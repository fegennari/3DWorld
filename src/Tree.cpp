// 3D World - OpenGL CS184 Computer Graphics Project
// by Hiral Patel and Frank Gennari
// 3/15/02

#include "3DWorld.h"
#include "mesh.h"
#include "tree_3dw.h"
#include "tree_leaf.h"
#include "explosion.h"
#include "physics_objects.h"
#include "gl_ext_arb.h"
#include "shaders.h"


float const BURN_RADIUS      = 0.2;
float const BURN_DAMAGE      = 400.0;
float const BEAM_DAMAGE      = 0.06;
float const LEAF_DAM_SCALE   = 0.00002;
float const LEAF_HEAL_RATE   = 0.00001;
float const TREE_SIZE        = 0.005; // 0.0015
float const REL_LEAF_SIZE    = 3.5;
float const TREE_MIN_H       = -0.15; // 0.15
float const TREE_MAX_H       = 0.6;
int   const LEAF_GEN_RAND1   = 10; //500
int   const LEAF_GEN_RAND2   = 200000; // larger is fewer leaves falling
float const DIST_C_SCALE     = 0.01;
float const LEAF_SHADOW_VAL  = 0.25;
float const BR_SHADOW_VAL    = 0.6;
float const TREE_DEPTH       = 0.1;
float const BASE_LEN_SCALE   = 0.8; // determines base cylinder overlap
float const BRANCH_LEN_SCALE = 0.9; // determines branch cylinder overlap
int const TREE_4TH_BRANCHES  = 0; // 1 = branches, 2 = branches + leaves
int const TREE_COLL          = 2;
int const DISABLE_LEAVES     = 0;
int const ENABLE_CLIP_LEAVES = 1;
int const TEST_RTREE_COBJS   = 0; // draw cobjs instead of tree (slow)

bool const USE_LEAF_GEOM_SHADER   = 0;
bool const USE_BRANCH_GEOM_SHADER = 0;


reusable_mem<tree_cylin >   tree::cylin_cache [CYLIN_CACHE_ENTRIES ];
reusable_mem<tree_branch>   tree::branch_cache[BRANCH_CACHE_ENTRIES];
reusable_mem<tree_branch *> tree::branch_ptr_cache;


// tree helper methods
void rotate_all(point const &rotate, float angle, float x, float y, float z);
int generate_next_cylin(int cylin_num, float branch_curveness, int ncib, bool branch_just_created, bool &branch_deflected);
void process_leaf(vector<tree_leaf> &leaves, point start, point rotate, float b_deg_rot, float tsize);


// tree_mode: 0 = no trees, 1 = large only, 2 = small only, 3 = both large and small
unsigned max_leaves(0), max_unique_trees(0);
int tree_mode(1), tree_coll_level(TREE_COLL);
float tree_temp_matrix[3], i_matrix[9], re_matrix[3];
float leaf_color_coherence(0.5), tree_color_coherence(0.2), tree_deadness(-1.0), nleaves_scale(1.0), leaf_size(0.0);
colorRGBA leaf_base_color(BLACK);
point leaf_points[4]; // z = 0.0 -> -0.05
tree_cont_t t_trees;


extern bool has_snow, no_sun_lpos_update;
extern int shadow_detail, island, num_trees, do_zoom, begin_motion, display_mode, animate2, iticks, draw_model;
extern int xoff2, yoff2, rand_gen_index, gm_blast, game_mode, leaf_color_changed, scrolling, dx_scroll, dy_scroll;
extern float zmin, zmax_est, zbottom, water_plane_z, mesh_scale, mesh_scale2, temperature, fticks, tree_size, vegetation;
extern point ocean;
extern lightning l_strike;
extern blastr latest_blastr;
extern texture_t textures[];
extern coll_obj_group coll_objects;
extern rand_gen_t global_rand_gen;


void calc_leaf_points() {

	leaf_size = REL_LEAF_SIZE*TREE_SIZE/(sqrt(nleaves_scale)*mesh_scale*mesh_scale2);
	leaf_points[0].assign(-2.0*leaf_size, 0.0, 0.0);
	leaf_points[1].assign(-2.0*leaf_size, 0.0, 4.0*leaf_size);
	leaf_points[2].assign( 2.0*leaf_size, 0.0, 4.0*leaf_size);
	leaf_points[3].assign( 2.0*leaf_size, 0.0, 0.0);
}


//tree related variables
extern GLUquadricObj* quadric;
int next_branch_num;


//designed to handle only positive numbers correctly
//return a random number between start and end
inline int rand_gen(int start, int end) {
	return (rand2()%(end - start + 1) + start);
}

float get_tree_z_bottom(float z, point const &pos) {
	return (is_over_mesh(pos) ? max(zbottom, (z - TREE_DEPTH)) : (z - TREE_DEPTH));
}


void mult_leaf_points_by(float val) {

	for (unsigned p = 0; p < 4; ++p) {
		leaf_points[p] *= val;
	}
}


inline colorRGBA get_leaf_base_color(int type) {

	colorRGBA color(tree_types[type].leafc);

	for (unsigned i = 0; i < 3; ++i) {
		color[i] = CLIP_TO_01(color[i] + leaf_base_color[i]);
	}
	return color;
}


bool tree::is_over_mesh() const {

	int const x1(get_xpos(sphere_center.x - sphere_radius)), y1(get_ypos(sphere_center.y - sphere_radius));
	int const x2(get_xpos(sphere_center.x + sphere_radius)), y2(get_ypos(sphere_center.y + sphere_radius));
	return (x1 < MESH_X_SIZE && y1 < MESH_Y_SIZE && x2 >= 0 && y2 >= 0); // completely off the mesh
}


void tree::gen_tree_shadows(unsigned light_sources) {

	if (shadow_detail < 2 || !(tree_mode & 1) || !created) return;
	// Note: not entirely correct since an off mesh tree can still cast a shadow on the mesh
	if (!is_over_mesh()) return; // optimization
	if (!enable_shadow_envelope(sphere_center, sphere_radius, light_sources, 1)) return;

	for (unsigned i = 0; i < all_cylins.size(); i++) {
		draw_cylin const &c(all_cylins[i]);
		if ((c.level + 2) > shadow_detail) return;
		cylinder_shadow(c.p1, c.p2, c.r1, c.r2, light_sources, 0, 0, (c.level < 2));
	}
	if (shadow_detail < 6) return;
	int const ltid(tree_types[type].leaf_tex);
	
	for (unsigned i = 0; i < leaves.size(); i++) { // loop through leaves
		polygon_shadow(leaves[i].pts, leaves[i].norm, 4, 0.0, light_sources, 0, 0, 0, ltid);
	}
	disable_shadow_envelope(light_sources);
}


void tree::add_tree_collision_objects() {

	//RESET_TIME;
	if (!(tree_mode & 1) || !tree_coll_level || !created) return;
	remove_collision_objects();
	if (!is_over_mesh()) return; // optimization
	int const btid(tree_types[type].bark_tex);
	colorRGBA const bcolor(tree_types[type].barkc);
	cobj_params cp(0.8, bcolor, TEST_RTREE_COBJS, 0, NULL, 0, btid, 4.0, 1, 0);
	cp.shadow = 0; // will be handled by gen_tree_shadows()

	for (unsigned i = 0; i < all_cylins.size(); i++) {
		draw_cylin &c(all_cylins[i]);

		if (c.level < min(tree_coll_level, 4)) {
			c.coll_index = add_coll_cylinder(c.p1, c.p2, c.r1, c.r2, cp, -1, ((c.level == 0) ? 0 : 1));
		}
		else {
			c.coll_index = -1;
		}
	}
	if (tree_coll_level >= 4) {
		int const ltid(tree_types[type].leaf_tex);
		colorRGBA lcolor(get_leaf_base_color(type).modulate_with(texture_color(ltid)));
		lcolor.alpha = 1.0;
		cobj_params cpl(0.3, lcolor, TEST_RTREE_COBJS, 0, NULL, 0, (TEST_RTREE_COBJS ? -1 : ltid), 1.0, 0, 0);
		cpl.shadow         = 0;
		cpl.is_destroyable = 1; // so that truly_static() returns false
		point const xlate(all_zeros); // for now

		for (unsigned i = 0; i < leaves.size(); i++) { // loop through leaves
			// Note: line collisions with leaves will use the texture alpha component for a more exact test
			leaves[i].coll_index = add_coll_polygon(leaves[i].pts, 4, cpl, 0.0, xlate, -1, 2);
			coll_objects[leaves[i].coll_index].is_billboard = 1;
		}
	}
	//PRINT_TIME("Tree Cobjs");
}


void tree::remove_collision_objects() {

	if (!created) return;

	for (unsigned i = 0; i < all_cylins.size(); i++) {
		remove_reset_coll_obj(all_cylins[i].coll_index);
	}
	if (tree_coll_level >= 5 || tree_coll_level >= 4) {
		for (unsigned i = 0; i < leaves.size(); ++i) {
			remove_reset_coll_obj(leaves[i].coll_index);
		}
	}
}


void remove_tree_cobjs() {

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].remove_collision_objects();
	}
}


void draw_trees_bl(shader_t const &s, bool draw_branches, bool draw_leaves, bool shadow_only) {

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].draw_tree(s, draw_branches, draw_leaves, shadow_only);
	}
}


void set_leaf_shader(shader_t &s, float min_alpha, bool gen_tex_coords, bool use_geom_shader) {

	s.set_prefix("#define USE_LIGHT_COLORS", 0); // VS
	if (gen_tex_coords) s.set_prefix("#define GEN_QUAD_TEX_COORDS", 0); // VS
	s.setup_enabled_lights(2);
	s.set_frag_shader("linear_fog.part+simple_texture");

	if (use_geom_shader) {
		s.set_vert_shader("ads_lighting.part*+leaf_lighting.part+tree_leaves_as_pts");
		s.set_geom_shader("output_textured_quad.part+point_to_quad", GL_POINTS, GL_TRIANGLE_STRIP, 4);
		s.begin_shader();
		//setup_wind_for_shader(s); // FIXME: add wind?
	}
	else {
		s.set_vert_shader("ads_lighting.part*+leaf_lighting.part+tree_leaves");
		s.begin_shader();
	}
	s.setup_fog_scale();
	s.add_uniform_float("min_alpha", min_alpha);
	set_multitex(0);
	s.add_uniform_int("tex0", 0);
	check_gl_error(301);
}


void check_leaf_shadow_change() {

	static point last_lpos(all_zeros);
	point const lpos(get_light_pos());
		
	if (!no_sun_lpos_update && lpos != last_lpos) {
		for (unsigned i = 0; i < t_trees.size(); ++i) {
			t_trees[i].calc_leaf_shadows();
		}
	}
	last_lpos = lpos;
}


void draw_trees(bool shadow_only) {

	//glFinish(); // testing
	//RESET_TIME;

	if (tree_mode & 2) { // small trees
		draw_small_trees(shadow_only);
		//PRINT_TIME("Small Trees");
	}
	if (tree_mode & 1) { // trees
		if (!shadow_only) check_leaf_shadow_change();
		set_fill_mode();

		// draw branches, then leaves: much faster for distant trees, slightly slower for near trees
		shader_t bs, ls;

		// draw branches
		bool const branch_smap(1 && !shadow_only); // looks better, but slower
		set_tree_branch_shader(bs, !shadow_only, !shadow_only, branch_smap, USE_BRANCH_GEOM_SHADER);
		draw_trees_bl(bs, 1, 0, shadow_only);
		bs.end_shader();
		disable_multitex_a();

		// draw leaves
		set_leaf_shader(ls, 0.75, !USE_LEAF_GEOM_SHADER, USE_LEAF_GEOM_SHADER);
		
		if (draw_model == 0) { // solid fill
			enable_blend();
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.75);
		}
		set_lighted_sides(2);
		set_specular(0.1, 10.0);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_NORMALIZE);
		draw_trees_bl(ls, 0, 1, shadow_only);
		set_lighted_sides(1);
		ls.end_shader();
		glDisable(GL_COLOR_MATERIAL);
		glEnable(GL_NORMALIZE);
		disable_blend();
		set_specular(0.0, 1.0);
		glDisable(GL_ALPHA_TEST);
		glDisable(GL_TEXTURE_2D);
		//glFinish(); // testing
		//PRINT_TIME(((tree_mode & 2) ? "Large + Small Trees" : "Large Trees"));
	}
}


void tree::mark_leaf_changed(unsigned i) {

	// Note: could keep track of the min/max changed leaf index and update a sub-range,
	//       but removing leaves will swap/resize the leaf data and invalidate the ranges,
	//       and many range updates may actually be slower than a full update
	leaves_changed = 1;
}


void tree::copy_color(colorRGB const &color, unsigned i) {

	if (USE_LEAF_GEOM_SHADER) {
		assert(i < leaf_data2.size());
		leaf_data2[i].set_c3(color);
	}
	else {
		assert(i < leaf_data.size());
	
		for (unsigned j = 0; j < 4; ++j) {
			leaf_data[j+(i<<2)].set_c3(color);
		}
	}
	if (i < leaves.size() && leaves[i].coll_index >= 0) { // update cobj color so that leaf water reflection is correct
		assert((unsigned)leaves[i].coll_index < coll_objects.size());
		coll_objects[leaves[i].coll_index].cp.color = colorRGBA(color).modulate_with(texture_color(tree_types[type].leaf_tex));
	}
	mark_leaf_changed(i);
}


void tree::change_leaf_color(colorRGBA &base_color, unsigned i) {

	if (!has_leaf_data()) return;
	colorRGB color3;
	colorRGBA const leafc(get_leaf_base_color(type));

	for (unsigned j = 0; j < 3; ++j) {
		color3[j] = leaves[i].color*(leaf_color_coherence*leafc[j] +
			(1.0 - leaf_color_coherence)*((j==0) ? leaves[i].lred : leaves[i].lgreen)) + base_color[j];
	}
	copy_color(color3, i);
}


void tree::remove_leaf(unsigned i, bool update_data) {

	assert(i < leaves.size());
	int const cix(leaves[i].coll_index);
	
	if (cix >= 0) {
		assert((unsigned)cix < coll_objects.size());
		remove_coll_object(cix);
	}
	leaves[i] = leaves.back();
	leaves.pop_back();
	if (!update_data || !has_leaf_data()) return;

	if (USE_LEAF_GEOM_SHADER) {
		unsigned const i4(i << 2), tnl4((unsigned)leaves.size() << 2);
		assert(leaves.size() <= leaf_data2.size());
		leaf_data2[i] = leaf_data2[leaves.size()]; // shift vertex array (last one is now invalid)
	}
	else {
		unsigned const i4(i << 2), tnl4((unsigned)leaves.size() << 2);
		assert(4*leaves.size() <= leaf_data.size());

		for (unsigned j = 0; j < 4; ++j) { // shift vertex array (last one is now invalid)
			leaf_data[j+i4] = leaf_data[j+tnl4];
		}
	}
	mark_leaf_changed(i);
}


void tree::burn_leaves() {

	float const max_t(get_max_t(LEAF));
	if (temperature < max_t || leaves.empty()) return;
	unsigned const num_burn(max(1U, min(5U, unsigned(5*(temperature - max_t)/max_t))));
	damage += ((1.0 - damage)*num_burn)/leaves.size();

	for (unsigned i = 0; i < num_burn && !leaves.empty(); ++i) {
		unsigned const index(rand()%leaves.size());
		leaves[index].color = max(0.0f, (leaves[index].color - 0.25f));
		change_leaf_color(base_color, index);
		if (rand()&1) gen_smoke(leaves[index].pts[0]);
		if (leaves[index].color <= 0.0) remove_leaf(index, 1);
	}
}


bool tree::damage_leaf(unsigned i, float damage_done) {

	if (damage_done == 0.0) return 0;
	assert(i < leaves.size());
	float &lcolor(leaves[i].color);
	colorRGB const black(0.0, 0.0, 0.0);

	if (damage_done > 1.0 || (damage_done > 0.2 && lcolor == 0.0)) {
		if (lcolor > 0.0) damage += damage_scale;
		lcolor = -1.0;
		if ((rand()&3) == 0) gen_leaf_at(leaves[i].pts, leaves[i].norm, type, black); // 50/50 chance of the burned leaf falling
		
		if (has_leaf_data()) { // remove leaf i
			remove_leaf(i, 1);
			return 1;
		}
	}
	else {
		float const llc(lcolor);
		lcolor = max(0.0, (lcolor - 0.3*damage_done));
		if (lcolor == 0.0 && llc > 0.0) damage += damage_scale;
		if (has_leaf_data()) change_leaf_color(base_color, i);
	}
	return 0;
}


void tree::blast_damage(blastr const *const blast_radius) {

	assert(blast_radius);
	float const bradius(blast_radius->cur_size), bdamage(LEAF_DAM_SCALE*blast_radius->damage);
	if (bdamage == 0.0) return;
	point const &bpos(blast_radius->pos);
	if (p2p_dist_sq(bpos, sphere_center) > (bradius + sphere_radius)*(bradius + sphere_radius)) return;
	float const bradius_sq(bradius*bradius);

	for (unsigned i = 0; i < leaves.size(); i++) {
		if (leaves[i].color < 0.0) continue;
		float const dist(p2p_dist_sq(bpos, leaves[i].pts[0]));
		if (dist < TOLERANCE || dist > bradius_sq) continue;
		float const blast_damage(bdamage*InvSqrt(dist));
		if (damage_leaf(i, blast_damage)) --i; // force reprocess of this leaf, wraparound to -1 is OK
	} // for i
	damage = min(1.0f, damage);
}


void tree::lightning_damage(point const &ltpos) {

	blastr const br(0, ETYPE_FIRE, -2, BURN_RADIUS, BURN_DAMAGE, ltpos, plus_z, LITN_C, LITN_C);
	blast_damage(&br);
}


void tree::drop_leaves() {

	if (damage >= 1.0 || leaves.empty()) return; // too damaged
	bool const llc(leaves_changed);
	unsigned const nleaves((unsigned)leaves.size());
	float const temp0(max(1.0f, min(0.3f, (20.0f-temperature)/30.0f)));
	int const rgen(min(LEAF_GEN_RAND2/10, int(rand_uniform(0.5, 1.5)*temp0*LEAF_GEN_RAND2/fticks)));

	for (unsigned i = (rand()%LEAF_GEN_RAND1); i < nleaves; i += LEAF_GEN_RAND1) {
		if ((rand()%rgen) == 0) {
			gen_leaf_at(leaves[i].pts, leaves[i].norm, type, get_leaf_color(i));
			// create a new leaf with a different color (and orient?)
			leaves[i].create_init_color(0);
			if (has_leaf_data()) copy_color(leaves[i].calc_leaf_color(leaf_color, base_color), i);
		}
	}
	leaves_changed = llc; // reset to original value so we don't do an update just because of this
}


void tree::gen_leaf_color() {

	float dc_scale(1.0 - 0.95*damage);
	colorRGBA const leafc(get_leaf_base_color(type));
	bcolor = tree_types[type].barkc;
	base_color.alpha = 1.0;

	for (unsigned i = 0; i < 3; ++i) {
		float cscale((i<2) ? 0.5 : 0.125);
		base_color[i] = rand_uniform2(-cscale*tree_color_coherence, cscale*tree_color_coherence);
		leaf_color[i] = leaf_color_coherence*leafc[i];
		bcolor[i]    *= min(1.0f, dc_scale*color[i]);
	}
}


colorRGB tree::get_leaf_color(unsigned i) const {

	assert(i < leaves.size());
	if (!has_leaf_data()) return leaves[i].calc_leaf_color(leaf_color, base_color);

	if (USE_LEAF_GEOM_SHADER) {
		assert(i < leaf_data2.size());
		return leaf_data2[i].get_c3();
	}
	else {
		assert(4*i < leaf_data.size());
		return leaf_data[i<<2].get_c3(); // return color of first vertex since they all should be the same
	}
}


inline colorRGB tree_leaf::calc_leaf_color(colorRGBA const &leaf_color, colorRGBA const &base_color) const {

	float const ilch(1.0 - leaf_color_coherence);
	return colorRGB(color*(leaf_color.R + ilch*lred) + base_color.R, color*(leaf_color.G + ilch*lgreen) + base_color.G, 0.0);
}


inline float tree_leaf::get_norm_scale(unsigned pt_ix) const {

	return ((shadow_bits & (1 << pt_ix)) ? LEAF_SHADOW_VAL : 1.0);
}


void tree::clear_vbo() {

	delete_vbo(branch_vbo);
	delete_vbo(branch_ivbo);
	delete_vbo(leaf_vbo);
	branch_vbo = branch_ivbo = leaf_vbo = 0;
	num_branch_quads = num_unique_pts   = 0;
}


bool tree::is_visible_to_camera() const {

	int const level((island || get_camera_pos().z > ztop) ? 1 : 2); // do we want to test the mesh in non-island mode?
	return sphere_in_camera_view(sphere_center, 1.1*sphere_radius, level);
}


void tree::draw_tree(shader_t const &s, bool draw_branches, bool draw_leaves, bool shadow_only) {

	if (!created) return;

	if (shadow_only) {
		if (!is_over_mesh()) return;
	
		if (draw_branches && branch_vbo > 0) { // draw branches (untextured)
			bind_vbo(branch_vbo, 0);

			if (USE_BRANCH_GEOM_SHADER) {
				draw_branches_as_lines(s, all_cylins.size());
			}
			else {
				size_t const branch_stride(sizeof(branch_vert_type_t));
				bind_vbo(branch_ivbo, 1);
				set_array_client_state(1, 0, 0, 0); // vertices only
				glVertexPointer(3, GL_FLOAT, branch_stride, 0);
				glDrawRangeElements(GL_QUADS, 0, num_unique_pts, num_branch_quads, GL_UNSIGNED_SHORT, 0);
				bind_vbo(0, 1);
			}
		}
		if (draw_leaves && leaf_vbo > 0 && !leaves.empty()) { // draw leaves
			select_texture(tree_types[type].leaf_tex);
			bind_vbo(leaf_vbo);

			if (USE_LEAF_GEOM_SHADER) {
				draw_leaves_as_points(s, leaves.size());
			}
			else {
				assert(leaf_data.size() >= 4*leaves.size());
				leaf_vert_type_t::set_vbo_arrays(); // could also disable normals and colors, but that doesn't seem to help much
				glDrawArrays(GL_QUADS, 0, 4*(unsigned)leaves.size());
			}
		}
		bind_vbo(0, 0);
		return;
	}
	set_rand2_state(trseed1, trseed2);
	gen_leaf_color();
	
	if (draw_leaves && deadness < 1.0 && init_deadness < 1.0) {
		burn_leaves();
		if (l_strike.enabled == 1)    lightning_damage(l_strike.end);
		if (game_mode && gm_blast)    blast_damage(&latest_blastr);
		if (begin_motion && animate2) drop_leaves();
	}
	if (TEST_RTREE_COBJS) return;
	if (!draw_leaves || draw_branches) {not_visible = !is_visible_to_camera();} // second pass only
	if (not_visible) return;
	bool const use_vbos(setup_gen_buffers());
	assert(use_vbos);
	float const size_scale((do_zoom ? ZOOM_FACTOR : 1.0)*base_radius/(distance_to_camera(sphere_center)*DIST_C_SCALE));
	if (draw_branches) draw_tree_branches(s, size_scale);
	if (draw_leaves && !has_no_leaves()) draw_tree_leaves(s, size_scale);
}


void tree::draw_branches_as_lines(shader_t const &s, unsigned num) {

	size_t const branch_stride(sizeof(branch_node_t));
	vert_norm::set_vbo_arrays(branch_stride);
	int const loc(s.get_attrib_loc("radius"));
	assert(loc > 0);
	glEnableVertexAttribArray(loc);
	glVertexAttribPointer(loc, 1, GL_FLOAT, GL_FALSE, branch_stride, (void *)sizeof(vert_norm));
	glDrawArrays(GL_LINES, 0, 2*num);
	glDisableVertexAttribArray(loc);
}


void tree::draw_tree_branches(shader_t const &s, float size_scale) {

	point const camera(get_camera_pos());
	unsigned const numcylin((unsigned)all_cylins.size());
	select_texture(tree_types[type].bark_tex);
	set_color(bcolor);
	BLACK.do_glColor();

	if (USE_BRANCH_GEOM_SHADER) {
		if (branch_vbo == 0) { // create vbo
			branch_vbo = create_vbo(); // do we need indexed vertices? no?
			assert(branch_vbo > 0);
			bind_vbo(branch_vbo, 0);
			vector<branch_node_t> data(2*numcylin);
			
			for (unsigned i = 0; i < numcylin; ++i) {
				bool const prev_connect(i   > 0        && all_cylins[i].can_merge(all_cylins[i-1]));
				bool const next_connect(i+1 < numcylin && all_cylins[i].can_merge(all_cylins[i+1]));
				data[(i<<1)+0].r = all_cylins[i].r1;
				data[(i<<1)+1].r = all_cylins[i].r2;
				data[(i<<1)+0].v = all_cylins[i].p1;
				data[(i<<1)+1].v = all_cylins[i].p2;
				data[(i<<1)+0].n = -all_cylins[i-prev_connect].get_norm_dir_vect();
				data[(i<<1)+1].n =  all_cylins[i+next_connect].get_norm_dir_vect();
			}
			upload_vbo_data(&data.front(), data.size()*sizeof(branch_node_t), 0);
		}
		else {
			assert(branch_vbo > 0);
			bind_vbo(branch_vbo, 0); // use vbo for rendering
		}
		unsigned const num(min(numcylin, max((numcylin/8), unsigned(18.0*numcylin*size_scale)))); // branch LOD
		draw_branches_as_lines(s, num);
		bind_vbo(0, 0);
	}
	else { // draw with branch vbos
		if (branch_vbo == 0) { // create vbos
			assert(branch_ivbo == 0);
			assert(num_branch_quads == 0 && num_unique_pts == 0);

			for (unsigned i = 0; i < numcylin; ++i) { // determine required data size
				bool const prev_connect(i > 0 && all_cylins[i].can_merge(all_cylins[i-1]));
				unsigned const ndiv(all_cylins[i].get_num_div());
				num_branch_quads += ndiv;
				num_unique_pts   += (prev_connect ? 1 : 2)*ndiv;
			}
			typedef unsigned short index_t;
			assert(num_unique_pts < (1 << 8*sizeof(index_t))); // cutting it close with 4th order branches
			vector<branch_vert_type_t> data;
			vector<index_t> idata;
			idata.reserve(4*num_branch_quads);
			data.reserve(num_unique_pts);
			unsigned cylin_id(0), data_pos(0), quad_id(0);

			for (unsigned i = 0; i < numcylin; i++) {
				draw_cylin const &cylin(all_cylins[i]);
				unsigned const ndiv(cylin.get_num_div());
				point const ce[2] = {cylin.p1, cylin.p2};
				float const ndiv_inv(1.0/ndiv);
				vector3d v12; // (ce[1] - ce[0]).get_norm()
				vector_point_norm const &vpn(gen_cylinder_data(ce, cylin.r1, cylin.r2, ndiv, v12, NULL, 0.0, 1.0, 0));
				bool const prev_connect(i > 0 && cylin.can_merge(all_cylins[i-1]));

				if (!prev_connect) { // new cylinder section
					data_pos = (unsigned)data.size();
					quad_id  = cylin_id = 0;
				}
				for (unsigned j = prev_connect; j < 2; ++j) { // create vertex data
					for (unsigned S = 0; S < ndiv; ++S) { // first cylin: 0,1 ; other cylins: 1
						float const tx(2.0*fabs(S*ndiv_inv - 0.5));
						// FIXME: something is still wrong - twisted branch segments due to misaligned or reversed starting points
						data.push_back(branch_vert_type_t(vpn.p[(S<<1)+j], vpn.n[S], tx, float(cylin_id + j))); // average normals?
					}
				}
				for (unsigned S = 0; S < ndiv; ++S) { // create index data
					bool const last_edge((quad_id % ndiv) == ndiv-1);
					unsigned const ix(data_pos + quad_id++);
					idata.push_back(ix);
					idata.push_back(ix+ndiv);
					idata.push_back(last_edge ? ix+1 : ix+ndiv+1);
					idata.push_back(last_edge ? ix+1-ndiv : ix+1);
				}
				++cylin_id;
			} // for i
			assert(data.size()  == data.capacity());
			assert(idata.size() == idata.capacity());
			branch_vbo  = create_vbo();
			branch_ivbo = create_vbo();
			assert(branch_vbo > 0 && branch_ivbo > 0);
			bind_vbo(branch_vbo,  0);
			bind_vbo(branch_ivbo, 1);
			upload_vbo_data(&data.front(),  data.size() *sizeof(branch_vert_type_t), 0); // ~350KB
			upload_vbo_data(&idata.front(), idata.size()*sizeof(index_t),            1); // ~75KB (with 16-bit index)
		} // end create vbo
		else {
			assert(branch_vbo > 0 && branch_ivbo > 0);
			bind_vbo(branch_vbo,  0); // use vbo for rendering
			bind_vbo(branch_ivbo, 1);
		}
		vert_norm_comp_tc::set_vbo_arrays();
		unsigned const num(4*min(num_branch_quads, max((num_branch_quads/8), unsigned(1.5*num_branch_quads*size_scale)))); // branch LOD
		glDrawRangeElements(GL_QUADS, 0, num_unique_pts, num, GL_UNSIGNED_SHORT, 0);
		bind_vbo(0, 0);
		bind_vbo(0, 1);
	}
}


void tree::draw_leaves_as_points(shader_t const &s, unsigned nl) {

	unsigned const leaf_stride(sizeof(leaf_node_t));
	vert_norm_color::set_vbo_arrays(leaf_stride);
	int const loc(s.get_attrib_loc("tangent"));
	assert(loc > 0);
	glEnableVertexAttribArray(loc);
	glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE, leaf_stride, (void *)sizeof(vert_norm_color));
	glDrawArrays(GL_POINTS, 0, nl);
	glDisableVertexAttribArray(loc);
}


void tree::draw_tree_leaves(shader_t const &s, float size_scale) {

	assert(leaves.size() <= max_leaves);
	bool const gen_arrays(!has_leaf_data()), create_leaf_vbo(leaf_vbo == 0);
	int const leaf_dynamic_en(!has_snow && (display_mode & 0x0100) != 0);

	if (gen_arrays && deadness > 0) {
		for (unsigned i = 0; i < leaves.size(); i++) { // process deadness stuff
			if (deadness > rand_float2()) {remove_leaf(i, 0);}
		}
	}
	if (!gen_arrays && leaf_dynamic_en && size_scale > 0.5) {update_leaf_orients();}
	unsigned nleaves((unsigned)leaves.size());
	unsigned nl(nleaves);
	if (ENABLE_CLIP_LEAVES) nl = min(nl, max((nl/8), unsigned(4.0*nl*size_scale))); // leaf LOD
	if (create_leaf_vbo) leaf_vbo = create_vbo();
	assert(leaf_vbo > 0);
	bind_vbo(leaf_vbo);
	unsigned const num_dlights(enable_dynamic_lights(sphere_center, sphere_radius));
	s.add_uniform_int("num_dlights", num_dlights);
	select_texture((draw_model == 0) ? tree_types[type].leaf_tex : WHITE_TEX); // what about texture color mod?

	if (USE_LEAF_GEOM_SHADER) {
		if (gen_arrays) leaf_data2.resize(nleaves);
		
		if (gen_arrays || (reset_leaves && !leaf_dynamic_en)) {
			for (unsigned i = 0; i < nleaves; i++) { // process leaf points - reset to default positions and normals
				leaf_data2[i].set_from_leaf(leaves[i]);
			}
			reset_leaves   = 0;
			leaves_changed = 1;
		}
		if (gen_arrays || leaf_color_changed) {copy_all_leaf_colors();}
		if (gen_arrays) {calc_leaf_shadows();}
		assert(leaf_data2.size() >= leaves.size());
		unsigned const leaf_stride(sizeof(leaf_node_t));

		if (create_leaf_vbo) {
			upload_vbo_data(&leaf_data2.front(), leaf_data2.size()*leaf_stride);
		}
		else if (leaves_changed) {
			upload_vbo_sub_data(&leaf_data2.front(), 0, leaf_data2.size()*leaf_stride);
		}
		draw_leaves_as_points(s, nl);
	}
	else {
		if (gen_arrays) {leaf_data.resize(4*nleaves);}

		if (gen_arrays || (reset_leaves && !leaf_dynamic_en)) {
			for (unsigned i = 0; i < nleaves; i++) { // process leaf points - reset to default positions and normals
				for (unsigned j = 0; j < 4; ++j) {
					leaf_data[j+(i<<2)].v = leaves[i].pts[j];
					leaf_data[j+(i<<2)].set_norm(leaves[i].norm*leaves[i].get_norm_scale(j));
				}
			}
			reset_leaves   = 0;
			leaves_changed = 1;
		}
		if (gen_arrays || leaf_color_changed) {copy_all_leaf_colors();}
		if (gen_arrays) {calc_leaf_shadows();}
		assert(leaf_data.size() >= 4*leaves.size());
		unsigned const leaf_stride(sizeof(leaf_vert_type_t));

		if (create_leaf_vbo) {
			upload_vbo_data(&leaf_data.front(), leaf_data.size()*leaf_stride); // ~150KB
		}
		else if (leaves_changed) {
			upload_vbo_sub_data(&leaf_data.front(), 0, leaf_data.size()*leaf_stride);
		}
		leaf_vert_type_t::set_vbo_arrays();
		glDrawArrays(GL_QUADS, 0, 4*nl);
	}
	disable_dynamic_lights(num_dlights);
	bind_vbo(0);
	leaves_changed = 0;
}


void leaf_node_t::set_from_leaf(tree_leaf const &l) {

	v = l.pts[0]; // origin is at LLC (TC 0,0)
	n = l.norm*l.get_norm_scale(0);
	t = (l.pts[3] - l.pts[0]); // not normalized
}


void tree::copy_all_leaf_colors() {

	for (unsigned i = 0; i < leaves.size(); i++) {
		copy_color(leaves[i].calc_leaf_color(leaf_color, base_color), i);
	}
}


void tree::update_leaf_orients() { // leaves move in wind or when struck by an object (somewhat slow)

	int last_xpos(0), last_ypos(0);
	vector3d local_wind;

	for (unsigned i = 0; i < leaves.size(); i++) { // process leaf wind and collisions
		tree_leaf const &l(leaves[i]);
		int const xpos(get_xpos(l.pts[0].x)), ypos(get_ypos(l.pts[0].y));
			
		// Note: should check for similar z-value, but z is usually similar within the leaves of a single tree
		if (i == 0 || xpos != last_xpos || ypos != last_ypos) {
			local_wind = get_local_wind(l.pts[0]); // slow
			last_xpos  = xpos;
			last_ypos  = ypos;
		}
		float hit_angle(0.0);

		if (l.coll_index >= 0) { // rotate leaves when hit by an object
			assert(l.coll_index < (int)coll_objects.size());
			unsigned char &last_coll(coll_objects[l.coll_index].last_coll);
			unsigned char &coll_type(coll_objects[l.coll_index].coll_type);
				
			if (last_coll > 0) {
				if (coll_type == BEAM) { // do burn damage
					if (damage_leaf(i, BEAM_DAMAGE)) {--i;}
				}
				else {
					if (coll_type == PROJECTILE) {
						if ((rand()&3) == 0) { // shoot off leaf
							gen_leaf_at(l.pts, l.norm, type, get_leaf_color(i));
							remove_leaf(i, 1);
							--i;
							continue;
						}
						coll_type = IMPACT; // reset to impact after first hit
					}
					hit_angle += PI_TWO*last_coll/TICKS_PER_SECOND; // 90 degree max rotate
				}
				last_coll = ((last_coll < iticks) ? 0 : (last_coll - iticks));
			}
		}
		float const wscale(dot_product(local_wind, l.norm));
		float const angle(0.5*PI*max(-1.0f, min(1.0f, wscale)) + hit_angle); // not physically correct, but it looks good
		point const p1((l.pts[1] + l.pts[2])*0.5); // tip
		point const p2((l.pts[0] + l.pts[3])*0.5); // base
		vector3d const orig_dir(p1 - p2); // vector from base to tip
		vector3d const new_dir(orig_dir*cos(angle) + l.norm*(orig_dir.mag()*sin(angle))); // s=orig_dir.get_norm(), t=l.norm
		vector3d const delta(new_dir - orig_dir);
		vector3d normal(cross_product(new_dir, (l.pts[3] - l.pts[0])).get_norm());
		
		if (USE_LEAF_GEOM_SHADER) {
			leaf_data2[i].n = normal*l.get_norm_scale(0);
		}
		else {
			unsigned const ix(i<<2);
			leaf_data[ix+1].v = l.pts[1] + delta;
			leaf_data[ix+2].v = l.pts[2] + delta;

			if (l.shadow_bits == 0) {
				norm_comp nc(normal);
				UNROLL_4X((norm_comp)leaf_data[i_+ix] = nc;)
			}
			else { // update the normals, even though this slows the algorithm down
				UNROLL_4X(leaf_data[i_+ix].set_norm(normal*l.get_norm_scale(i_));)
			}
		}
		if (LEAF_HEAL_RATE > 0.0 && l.color > 0.0 && l.color < 1.0) { // leaf heal
			leaves[i].color = min(1.0f, (l.color + LEAF_HEAL_RATE*fticks));
			copy_color(l.calc_leaf_color(leaf_color, base_color), i);
		}
		reset_leaves = 1; // Do we want to update the normals and collision objects as well?
		mark_leaf_changed(i);
	} // for i
}


void tree::calc_leaf_shadows() { // process leaf shadows/normals

	if (!has_leaf_data()) return; // leaf data not yet created (can happen if called when light source changes)
	int const light(get_light());

	for (unsigned i = 0; i < leaves.size(); i++) {
		tree_leaf &l(leaves[i]);

		if (USE_LEAF_GEOM_SHADER) {
			l.shadow_bits = (l.coll_index >= 0 && !is_visible_to_light_cobj(l.get_center(), light, 0.0, l.coll_index, 1));
			leaf_data2[i].n = l.norm*l.get_norm_scale(0);
		}
		else {
			l.shadow_bits = 0;

			for (unsigned j = 0; j < 4; ++j) {
				bool const shadowed(l.coll_index >= 0 && !is_visible_to_light_cobj(l.pts[j], light, 0.0, l.coll_index, 1));
				l.shadow_bits |= (int(shadowed) << j);
				leaf_data[j+(i<<2)].set_norm(l.norm*l.get_norm_scale(j));
			}
		}
	}
	leaves_changed = 1;
}


void delete_trees() {

	unsigned deleted(0);

	for (unsigned i = (unsigned)t_trees.size(); i > 0; --i) { // delete backwards (pop collision stack)
		if (t_trees[i-1].delete_tree()) ++deleted;
	}
	if (tree_coll_level && !t_trees.empty()) purge_coll_freed(1); // MUST do this for collision detection to work
}


int tree::delete_tree() {

	clear_vbo();
	if (!created)  return 0;
	if (tree_coll_level) remove_collision_objects();
	if (no_delete) return 0;
	all_cylins.clear();
	leaf_data.clear();
	leaf_data2.clear();
	created = 0;
	return 1;
}


inline void update_radius(point const &end, point const &sc, float &radius) {

	float const distance = (end.x-sc.x)*(end.x-sc.x) + (end.y-sc.y)*(end.y-sc.y) + (end.z-sc.z)*(end.z-sc.z);
	if (distance > radius*radius) radius = sqrt(distance);
}


void tree::process_cylins(tree_cylin const *const cylins, unsigned num) {

	for (unsigned i = 0; i < num; ++i) {
		assert(cylins[i].r1 > 0.0 || cylins[i].r2 > 0.0);
		assert(cylins[i].p1 != cylins[i].p2);
		update_radius(cylins[i].p2, sphere_center, sphere_radius);
		all_cylins.push_back(cylins[i]);

		if (leaves.capacity() > 0) { // leaves was reserved
			if (cylins[i].level > 1 && (cylins[i].level < 4 || TREE_4TH_BRANCHES > 1)) { // leaves will still be allocated
				float tsize((mesh_scale*base_radius/(TREE_SIZE*tree_size) + 10.0/mesh_scale2)/18.0);
				if (type == 3) tsize *= 1.5; // this tree type has larger leaves
				add_leaves_to_cylin(cylins[i], tsize);
			}
		}
	}
}


void tree::create_leaves_and_one_branch_array() {

	assert(all_cylins.empty());

	//tree leaf variables
	if (deadness < 1.0 && !DISABLE_LEAVES) {
		unsigned nl(unsigned((1.0 - init_deadness)*num_leaves_per_occ*all_cylins.size()) + 1); // determine the number of leaves
		leaves.reserve(nl);
	}

	//set the bounding sphere center
	assert(base_num_cylins > 0);
	sphere_center.z = base.cylin[max(0, (base_num_cylins - 2))].p2.z;
	sphere_radius   = 0.0;

	//process cylinders
	unsigned num_total_cylins(base_num_cylins); // start with trunk cylinders

	for (int i = 0; i < num_1_branches; i++) {
		for (int k = 0; k < (branches[i][0].num_branches + 1); k++) {
			num_total_cylins += branches[i][k].num_cylins;
		}
	}
	for (unsigned w = 0; w < 2; ++w) {
		for (int i = 0; i < num_34_branches[w]; i++) {
			num_total_cylins += branches_34[w][i].num_cylins;
		}
	}
	all_cylins.reserve(num_total_cylins);
	process_cylins(base.cylin, base_num_cylins);

	for (int i = 0; i < num_1_branches; i++) { //add the first and second order branches
		for (int k = 0; k < (branches[i][0].num_branches + 1); k++) {
			process_cylins(branches[i][k].cylin, branches[i][k].num_cylins);
		}
	}
	for (unsigned w = 0; w < 2; ++w) {
		for (int i = 0; i < num_34_branches[w]; i++) { //add the third and fourth order branches
			process_cylins(branches_34[w][i].cylin,  branches_34[w][i].num_cylins);
		}
	}
	assert(all_cylins.size() == all_cylins.capacity());

	//now delete the branches b/c they are unneccessary
	cylin_cache [2].reusable_free(branches_34[0][0].cylin);
	branch_cache[1].reusable_free(branches_34[0]);
	cylin_cache [1].reusable_free(branches[0][0].cylin);
	branch_cache[0].reusable_free(branches[0]);
	branch_ptr_cache.reusable_free(branches);
	cylin_cache [0].reusable_free(base.cylin);

	if (TREE_4TH_BRANCHES) {
		cylin_cache [3].reusable_free(branches_34[1][0].cylin);
		branch_cache[2].reusable_free(branches_34[1]);
	}
	max_leaves = max(max_leaves, unsigned(leaves.size()));

	/*for (unsigned i = 0; i < leaves.size(); ++i) { // scramble leaves so that LOD is unbiased/randomly sampled
		swap(leaves[i], leaves[(i + 1572869)%leaves.size()]);
	}*/
	reverse(leaves.begin(), leaves.end()); // order leaves so that LOD removes from the center first, which is less noticeable
	if (!leaves.empty()) damage_scale = 1.0/leaves.size();
}


inline void add_rotation(point &dest, point const &src, float mult) {

	for (unsigned i = 0; i < 3; ++i) {
		dest[i] = src[i] + mult*re_matrix[i];
	}
}


inline void setup_rotate(vector3d &rotate, float rotate_start, float temp_deg) {

	float const angle(rotate_start/TO_DEG + temp_deg);
	rotate.assign(cosf(angle), sinf(angle), 0.0);
}


inline void rotate_around_axis(tree_cylin const &c) {

	rotate_all(c.rotate, c.deg_rotate/TO_DEG, 0.0, 0.0, c.length);
}


inline void rotate_pts_around_axis(point const &p, point &rotation_v, float deg_rotate) {

	rotate_all(rotation_v, deg_rotate/TO_DEG, p.x, p.y, p.z);
}


inline void rotate_cylin(tree_cylin &c) {

	//rotate_around_axis(c);
	rotate_all(c.rotate, c.deg_rotate/TO_DEG, 0.0, 0.0, c.length);
	add_rotation(c.p2, c.p1, 1.0);
}


inline void rotate_it() {

	re_matrix[0] = i_matrix[0]*tree_temp_matrix[0] + i_matrix[1]*tree_temp_matrix[1];
	re_matrix[1] = i_matrix[3]*tree_temp_matrix[0] + i_matrix[4]*tree_temp_matrix[1];
	re_matrix[2] = tree_temp_matrix[2];
	tree_temp_matrix[0] = re_matrix[0];
	tree_temp_matrix[1] = re_matrix[1];
}

inline void rotate_itv2() {

	re_matrix[0] = tree_temp_matrix[0];
	re_matrix[1] = i_matrix[4]*tree_temp_matrix[1] + i_matrix[5]*tree_temp_matrix[2];
	re_matrix[2] = i_matrix[7]*tree_temp_matrix[1] + i_matrix[8]*tree_temp_matrix[2];
	tree_temp_matrix[1] = re_matrix[1];
	tree_temp_matrix[2] = re_matrix[2];
}


inline void turn_into_i(float *m) {

	m[0] = 1.0; m[1] = 0.0; m[2] = 0.0;
	m[3] = 0.0; m[4] = 1.0; m[5] = 0.0;
	m[6] = 0.0; m[7] = 0.0; m[8] = 1.0;
}


void rotate_all(point const &rotate, float angle, float x, float y, float z) {

	static float langle(0.0), sin_term(0.0), cos_term(1.0);

	if (angle != langle) {
		cos_term = cos(angle);
		sin_term = sin(angle);
		langle   = angle;
	}

	//create the starting point
	tree_temp_matrix[0] = x;
	tree_temp_matrix[1] = y;
	tree_temp_matrix[2] = z;

	//rotate around z
	turn_into_i(i_matrix);
	i_matrix[0] = -rotate.x;
	i_matrix[1] =  rotate.y;
	i_matrix[3] = -rotate.y;
	i_matrix[4] = -rotate.x;
	rotate_it(); //do matrix multiplication

	//rotate around x;
	turn_into_i(i_matrix);
	i_matrix[4] =  cos_term;
	i_matrix[5] = -sin_term;
	i_matrix[7] =  sin_term;
	i_matrix[8] =  cos_term;
	rotate_itv2(); //do matrix multiplication

	//do_reverse, inverse rotate around z
	turn_into_i(i_matrix);
	i_matrix[0] =  rotate.x;
	i_matrix[1] = -rotate.y;
	i_matrix[3] =  rotate.y;
	i_matrix[4] =  rotate.x;
	rotate_it(); //do matrix multiplication
}


void gen_cylin_rotate(vector3d &rotate, vector3d &lrotate, float rotate_start) {

	float temp_deg(safe_acosf(lrotate.x));
	if (lrotate.y < 0.0) temp_deg *= -1.0;
	setup_rotate(rotate, rotate_start, temp_deg);
}



void tree::gen_tree(point &pos, int size, int ttype, int calc_z, bool add_cobjs) {

	sphere_center = pos; // z value will be reset later
	if (calc_z) pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
	leaf_data.clear();
	leaf_data2.clear();

	//fixed tree variables
	type     = ((ttype < 0) ? rand2() : ttype)%NUM_TREE_TYPES; // maybe should be an error if > NUM_TREE_TYPES
	created  = 1;
	deadness = 0.0;
	damage   = 0.0;
	damage_scale = 0.0;
	leaf_vbo = branch_vbo = branch_ivbo = 0;
	color.alpha = 1.0;

	for (unsigned i = 0; i < 3; ++i) {
		color[i] = 0.01*rand_gen(80, 120);
	}
	if (tree_deadness >= 0.0) {
		init_deadness = tree_deadness;
	}
	else {
		int const num(rand_gen(1, 100));
		init_deadness = ((num > 94) ? min(1.0f, float(num - 94)/8.0f) : 0.0);
	}
	ncib                 = 10;
	base_num_cylins      = 5;
	num_cylin_factor     = 10.0;
	base_cylin_factor    = 5.0; //(type == TREE_A) ? 2.0 : 5.0;
	base_break_off       = 3;
	num_1_branches       = 8;
	num_big_branches_min = 3;
	num_big_branches_max = 4;
	num_2_branches_min   = 4;
	num_2_branches_max   = 6;
	num_34_branches[0]   = 0; // this should start out 0
	num_3_branches_min   = 6;
	num_3_branches_max   = 10;
	num_34_branches[1]   = 0;
	if (size <= 0) size  = rand_gen(40, 80); // tree size
	base_radius            = size * (0.1*tree_size*TREE_SIZE*tree_types[type].size/(mesh_scale*mesh_scale2));
	num_leaves_per_occ     = 0.01*nleaves_scale*(rand_gen(30, 60) + size);
	base_length_min        = rand_gen(4, 6) * base_radius;
	base_length_max        = base_length_min * 1.5;
	angle_rotate           = 60.0;
	base_curveness         = 10.0;
	tree_slimness          = 0;
	tree_wideness          = 70;
	branch_curveness       = 90.0;
	base_var               = 0.8;
	branch_1_var           = 100.0*0.85;
	branch_1_rad_var       = 100.0*0.64;
	branch_1_start         = 0.45;
	branch_2_start         = 0.9;
	branch_2_var           = 100.0*0.85;
	branch_2_rad_var       = 100.0*0.64;
	branch_4_max_radius    = 0.0002;
	rotate_factor          = 1.0; // how much to changes the horizontal rotation 
	branch_distribution    = 1.0;
	branch_1_distribution  = 1.0;
	branch_upwardness      = 0.9;
	branch_min_angle       = 20.0; // minimim angle of deflection from the current branch creating it
	branch_max_angle       = 40.0;
	max_2_angle_rotate     = 50.0;
	max_3_angle_rotate     = 50.0;
	branch_1_random_rotate = 40.0;
	trseed1                = rand2();
	trseed2                = rand2();

	//temporary variables
	int num_b_so_far(0);
	int const nbr(num_1_branches*(num_2_branches_max+1)), nbranches(nbr*num_3_branches_max);

	//allocate all the memory required for the tree------------------------------------------------
	cylin_cache [0].reusable_malloc(base.cylin,              base_num_cylins+1);
	branch_ptr_cache.reusable_malloc(branches,                num_1_branches);
	branch_cache[0].reusable_malloc(branches[0],             nbr);
	branch_cache[1].reusable_malloc(branches_34[0],          nbranches);
	cylin_cache [1].reusable_malloc(branches[0][0].cylin,    nbr*ncib);
	cylin_cache [2].reusable_malloc(branches_34[0][0].cylin, nbranches*ncib);

	for (int i = 1; i < num_1_branches; i++) { // start at 1
		branches[i] = branches[0] + i*(num_2_branches_max+1);
	}
	for (int i = 0; i < num_1_branches; i++) {
		for (int j = 0; j < (num_2_branches_max+1); j++) {
			branches[i][j].cylin = branches[0][0].cylin + (i*(num_2_branches_max+1) + j)*ncib;
		}
	}
	for (int i = 1; i < nbranches; i++) { // start at 1
		branches_34[0][i].cylin = branches_34[0][0].cylin + i*ncib;
	}

	//create tree base -------------------------------------------------------------
	int const init_num_cyl(base_num_cylins - base_break_off + 1);
	float const length((base_length_max - base_length_min)/2.0);
	float const cos_ar(cosf(angle_rotate/TO_DEG)), sin_ar(sinf(angle_rotate/TO_DEG));

	for (int i = 0; i < base_num_cylins; i++) {
		tree_cylin &cylin(base.cylin[i]);
		int num_b_to_create(0);

		if (i == 0) {
			cylin.assign_params(0, 0, base_radius, base_radius*base_var, length*(num_cylin_factor/base_num_cylins), 0.0);
			cylin.p1 = pos;
			cylin.rotate.assign(cos_ar, sin_ar, 0.0);
		}
		else {
			tree_cylin &lcylin(base.cylin[i-1]);
			cylin.assign_params(0, 0, lcylin.r2, lcylin.r2*base_var, lcylin.length*base_var*(base_cylin_factor/base_num_cylins),
				int(sinf(-PI_TWO + i*PI/base_num_cylins)*base_curveness));
			cylin.rotate.assign(lcylin.rotate.x, lcylin.rotate.y, 0.0);
			rotate_cylin(lcylin);
			add_rotation(cylin.p1, lcylin.p1, BASE_LEN_SCALE);
		}
		if (i == (base_num_cylins-1)) rotate_cylin(cylin); // last cylin

		if (base_break_off <= (i+1)) {
			float const temp_num(branch_distribution*((float)num_1_branches)/init_num_cyl*(float(i+2-base_break_off))/(float(num_b_so_far+1)));
			float rotate_start(rand_gen(0, 259));

			if (temp_num >= 1.0) {
				float const temp2(((float)num_1_branches - num_b_so_far)/(float(base_num_cylins - i)));
				num_b_to_create = ((temp2 <= 3.0) ? max(1, int(ceil(temp2))) : int(temp_num + 0.5));
			}
			if (num_b_to_create > num_1_branches) num_b_to_create = num_1_branches;
			
			for (int j = 0; j < num_b_to_create; j++) {
				if (num_b_so_far < num_1_branches) create_1_order_branch(i, rotate_start, num_b_so_far++);
				rotate_start += 360.0/((float)num_b_to_create) + (3 - 2*rand_gen(1,2))*rand_gen(0,int(branch_1_random_rotate));
			}
			if (num_b_to_create < 1) {
				if ((num_1_branches > 0) && ((i+1) == base_break_off)) {
					if (num_b_so_far < num_1_branches) create_1_order_branch(i, rotate_start, num_b_so_far++);
				}
			}
		}
	}
	if (num_1_branches > num_b_so_far) num_1_branches = num_b_so_far; // *** FIXED: Must update #branches to what has been allocated. ***
	
	//done with base ------------------------------------------------------------------
	if (TREE_4TH_BRANCHES) create_4th_order_branches();

	// cylinder from base into ground
	bool const mesh_disabled(is_mesh_disabled(get_xpos(pos.x), get_ypos(pos.y)));
	float const lmp(mesh_disabled ? (pos.z - TREE_DEPTH) : get_tree_z_bottom(lowest_mesh_point(pos, base_radius), pos));

	if (pos.z > lmp) { // add the bottom cylinder section (possibly to the bottom of the mesh)
		tree_cylin &cylin(base.cylin[base_num_cylins]);
		cylin.assign_params(0, 1, base_radius, base_radius, (pos.z - lmp), 180.0);
		cylin.p1 = cylin.p2 = pos;
		cylin.rotate.assign(1.0, 0.0, 0.0);
		rotate_cylin(cylin);
		++base_num_cylins;
	}
	create_leaves_and_one_branch_array();
	if (add_cobjs) add_tree_collision_objects();
}


inline float tree::gen_bc_size(float branch_var) {

	return (rand_gen((int)branch_var - 5, (int)branch_var + 5)/100.0)*(num_cylin_factor/ncib);
}


inline float tree::gen_bc_size2(float branch_var) {

	return rand_gen((int)branch_var - 5, (int)branch_var + 5)/100.0;
}


void tree::gen_next_cylin(tree_cylin &cylin, tree_cylin &lcylin, float var, float rad_var, int level, int branch_id, bool rad_var_test) {

	cylin.assign_params(level, branch_id, lcylin.r2, (rad_var_test ? lcylin.r1*gen_bc_size(rad_var) : 0.0),
		((level == 4) ? branch_4_length*var : lcylin.length*gen_bc_size(var)), lcylin.deg_rotate);
	rotate_around_axis(lcylin);
	add_rotation(cylin.p1, lcylin.p1, BRANCH_LEN_SCALE);
}


void tree::gen_first_cylin(tree_cylin &cylin, tree_cylin &src_cylin, float bstart, float rad_var, float rotate_start, int level, int branch_id) {

	rotate_around_axis(src_cylin);
	add_rotation(cylin.p1, src_cylin.p1, BRANCH_LEN_SCALE);
	float const radius1(bstart*src_cylin.r2);
	float deg_rotate(rand_gen(0, int(branch_max_angle)));
	
	if (rotate_start < branch_min_angle && deg_rotate < branch_min_angle) {
		deg_rotate = rand_gen(int(branch_min_angle - rotate_start), int(branch_max_angle));
	}
	if (rand2() & 1) deg_rotate *= -1.0;
	if (level == 2 && deg_rotate < 0.0) deg_rotate *= 0.5;
	deg_rotate += src_cylin.deg_rotate;
	cylin.assign_params(level, branch_id, radius1, radius1*gen_bc_size2(rad_var), bstart*src_cylin.length*(num_cylin_factor/ncib), deg_rotate);
}


void tree::create_1_order_branch(int base_cylin_num, float rotate_start, int branch_num) {

	bool branch_just_created(false), branch_deflected(false);
	int rotation((rand2() & 1) ? -1 : 1);

	//first cylinders of the first order branches
	tree_branch &branch(branches[branch_num][0]);
	tree_cylin &cylin(branch.cylin[0]);
	branch.num_cylins = ncib;
	rotate_around_axis(base.cylin[base_cylin_num]);
	add_rotation(cylin.p1, base.cylin[base_cylin_num].p1, BRANCH_LEN_SCALE);
	setup_rotate(cylin.rotate, rotate_start, 0.0);
	float const radius1(base_radius*branch_1_start);
	cylin.assign_params(1, branch_num, radius1, radius1*gen_bc_size(branch_1_rad_var),
		radius1*6*((float)rand_gen(7,10)/10)*(num_cylin_factor/ncib),
		(tree_slimness + (tree_wideness-tree_slimness)*rand_gen(6,10)/10.0*(num_1_branches - branch_num)/num_1_branches));

	//generate stats for the second order branches
	branch.num_branches = rand_gen(num_2_branches_min, num_2_branches_max);

	//create rest of the cylinders
	int num_2_branches_created(0);
	int const temp_num_big_branches(rand_gen(num_big_branches_min, num_big_branches_max));
	float const temp_num((float)temp_num_big_branches/((num_2_branches_created+1)*(float(0.5*branch.num_cylins+1))) +
		                  float(temp_num_big_branches/((float(branch.num_cylins/2+1))*(num_2_branches_created+1))));

	if (temp_num*branch_1_distribution >= 1.0) {
		if (num_2_branches_created < branch.num_branches) {
			create_2nd_order_branch(branch_num, ++num_2_branches_created, 0, branch_deflected, rotation);
			branch_just_created = true;
		}
	}
	if (branch.num_cylins < 2) rotate_cylin(cylin);
	
	for (int j = 1; j < branch.num_cylins; j++) {
		tree_cylin &cylin(branch.cylin[j]), &lcylin(branch.cylin[j-1]);
		branch_just_created = false;
		gen_next_cylin(cylin, lcylin, branch_1_var, branch_1_rad_var, 1, branch_num, (j < branch.num_cylins-1));
		float temp_num2;

		//calculate p1, based on the end point of lasy cylinder--------------------
		if (j < 0.5*branch.num_cylins) {
			temp_num2 = float((j+1)*temp_num_big_branches)/((num_2_branches_created+1)*(float(0.5*branch.num_cylins+1))) +
				        float(temp_num_big_branches/((float(branch.num_cylins/2+1))*(num_2_branches_created+1)));
		}
		else {
			temp_num2 = (float(((j+1)*branch.num_branches)/((num_2_branches_created+1)*branch.num_cylins))) +
				        (float(branch.num_branches/(branch.num_cylins*(num_2_branches_created+1))));
		}
		if (temp_num2*branch_1_distribution >= 1.0 && num_2_branches_created < branch.num_branches) branch_just_created = true;
		int const deg_added(generate_next_cylin(j, branch_curveness, branch.num_cylins, branch_just_created, branch_deflected));
		cylin.deg_rotate += deg_added;
		cylin.deg_rotate *= branch_upwardness;

		//set rotate ------------
		gen_cylin_rotate(cylin.rotate, lcylin.rotate, deg_added*rotate_factor);
		add_rotation(lcylin.p2, lcylin.p1, 1.0);
		if (j == (branch.num_cylins - 1)) rotate_cylin(cylin);
		
		if (branch_just_created) {
			if (j == branch.num_cylins-1) {
				branches[branch_num][num_2_branches_created+1].clear_num();
			}
			else {
				rotation = -rotation;
				create_2nd_order_branch(branch_num, ++num_2_branches_created, j, branch_deflected, rotation);
			}
		}
	}
}


void tree::create_2nd_order_branch(int i, int j, int cylin_num, bool branch_deflected, int rotation) {
	
	int index(0), num_3_branches_created(0);
	float const rotate_start((float)rand_gen(0, int(max_2_angle_rotate)));
	bool branch_just_created(false);
	
	//first cylinders of the second order branches
	tree_branch &branch(branches[i][j]);
	tree_cylin &cylin(branch.cylin[0]), &src_cylin(branches[i][0].cylin[cylin_num]);
	branch.num_cylins = branches[i][0].num_cylins - cylin_num;
	gen_first_cylin(cylin, src_cylin, branch_2_start, branch_2_rad_var, rotate_start, 2, j);
	gen_cylin_rotate(cylin.rotate, src_cylin.rotate, rotate_start*rotation);

	//generate stats for the third order branches
	branch.num_branches = rand_gen(num_3_branches_min, num_3_branches_max);
	int const temp_num_big_branches(branch.num_branches);
	float const temp_num((float)temp_num_big_branches/((num_3_branches_created+1)*((float)branch.num_cylins)) +
		                  float(temp_num_big_branches/(((float)branch.num_cylins)*(num_3_branches_created+1))));

	if (temp_num*branch_1_distribution >= 1.0) {
		if (num_3_branches_created < branch.num_branches && cylin.r2 > branch_4_max_radius) {
			create_3rd_order_branch(i, j, 0, ((num_3_branches_created++) + num_34_branches[0]), branch_deflected, rotation);
			branch_just_created = true;
		}
	}
	if (branch.num_cylins < 2) rotate_cylin(cylin);

	//create rest of the cylinders
	for (index = 1; index < branch.num_cylins; index++) {
		tree_cylin &cylin(branch.cylin[index]), &lcylin(branch.cylin[index-1]);
		gen_next_cylin(cylin, lcylin, branch_2_var, branch_2_rad_var, 2, j, (index < branch.num_cylins-1));
		float const temp_num2(float((index+1)*temp_num_big_branches)/((num_3_branches_created+1)*((float)branch.num_cylins)) +
			float(temp_num_big_branches/(((float)branch.num_cylins)*(num_3_branches_created+1))));
		if (temp_num*branch_1_distribution >= 1.0 && num_3_branches_created < branch.num_branches) branch_just_created = true;
		int const deg_added(generate_next_cylin(index, branch_curveness, branch.num_cylins, branch_just_created, branch_deflected));
		gen_cylin_rotate(cylin.rotate, lcylin.rotate, deg_added*rotate_factor);
		cylin.rotate.z    = lcylin.rotate.z;
		cylin.deg_rotate += ((branch.cylin[0].deg_rotate < 0.0) ? deg_added : -deg_added);
		cylin.deg_rotate *= branch_upwardness;
		add_rotation(lcylin.p2, lcylin.p1, 1.0);
		if (index == (branch.num_cylins - 1)) rotate_cylin(cylin);
		
		if (branch_just_created && cylin.deg_rotate < branch_min_angle) {
			int const bnum(num_3_branches_created + num_34_branches[0]);

			if (index == branch.num_cylins-1) {
				branches_34[0][bnum].clear_num();
			}
			else {
				rotation = -rotation;
				create_3rd_order_branch(i, j, index, bnum, branch_deflected, rotation);
				num_3_branches_created++;
			}
		}
	}
	num_34_branches[0] += num_3_branches_created;
}


void tree::create_3rd_order_branch(int i, int j, int cylin_num, int branch_num, bool branch_deflected, int rotation) {
	
	int index(0);
	float const rotate_start((float)rand_gen(0, int(max_3_angle_rotate)));
	bool branch_just_created(false);

	//first cylinders of the third order branches
	tree_branch &branch(branches_34[0][branch_num]);
	tree_cylin &cylin(branch.cylin[0]), &src_cylin(branches[i][j].cylin[cylin_num]);
	branch.num_cylins = branches[i][j].num_cylins - cylin_num;
	gen_first_cylin(cylin, src_cylin, branch_2_start, branch_2_rad_var, rotate_start, 3, branch_num);
	gen_cylin_rotate(cylin.rotate, src_cylin.rotate, rotate_start*rotation);

	//generate stats for the third order branches
	branch.num_branches = rand_gen(num_3_branches_min, num_3_branches_max);
	if (branch.num_cylins < 2) rotate_cylin(cylin);

	//create rest of the cylinders
	for (index = 1; index < branch.num_cylins; index++) {
		tree_cylin &cylin(branch.cylin[index]), &lcylin(branch.cylin[index-1]);
		gen_next_cylin(cylin, lcylin, branch_2_var, branch_2_rad_var, 3, branch_num, (index < branch.num_cylins-1));
		cylin.rotate = lcylin.rotate;
		int const deg_added(generate_next_cylin(index, branch_curveness, branch.num_cylins, branch_just_created, branch_deflected));
		cylin.deg_rotate += ((branch.cylin[0].deg_rotate < 0.0) ? deg_added : -deg_added);
		add_rotation(lcylin.p2, lcylin.p1, 1.0);
		if (index == (branch.num_cylins - 1)) rotate_cylin(cylin);
	}
}


void tree::gen_b4(tree_branch &branch, int &branch_num, int i, int k) {

	int ncib(branch.num_cylins);

	for (int j = 0; j < ncib-1; j++) {
		if ((float(((float)j)/((float)ncib))) >= branch_4_distribution) {
			tree_cylin &cylin(branch.cylin[j]);
			float temp_deg(safe_acosf(cylin.rotate.x)), rotate_start(0.0);
			if (cylin.rotate.y < 0.0) temp_deg *= -1.0;

			for (int l = 0; l < num_4_branches_per_occurance; l++) {
				rotate_start += 360.0/num_4_branches_per_occurance*l + rand_gen(1,30);
				if (rotate_start > 360.0) rotate_start -= 360.0;
				generate_4th_order_branch(branch, j, rotate_start, temp_deg, branch_num++);
			}
		}
	}
}


void tree::create_4th_order_branches() {

	num_34_branches[1]           = 2000;
	branch_4_distribution        = 0.2;
	num_4_branches_per_occurance = 2;
	num_4_cylins                 = 10;
	branch_4_rad_var             = 100.0*0.85;
	branch_4_var                 = 0.70;
	branch_4_length              = 0.006; //0.03;
	branch_4_max_radius          = 0.008;
	int branch_num(0);
	assert(num_34_branches[1] > 0);
	branch_cache[2].reusable_malloc(branches_34[1],          num_34_branches[1]);
	cylin_cache [3].reusable_malloc(branches_34[1][0].cylin, num_34_branches[1]*ncib);

	for (int i = 1; i < num_34_branches[1]; i++) { // skipping first branch
		branches_34[1][i].cylin = branches_34[1][0].cylin + i*ncib;
	}
	for (int i = 0; i < num_1_branches; i++) { // b1->b4
		for (int k = 0; k < (branches[i][0].num_branches + 1); k++) {
			gen_b4(branches[i][k], branch_num, i, k); // b2->b4
		}
	}
	for (int i = 0; i < num_34_branches[0]; i++) { // b2->b4
		gen_b4(branches_34[0][i], branch_num, i, 0);
	}
	num_34_branches[1] = branch_num;
}


void tree::generate_4th_order_branch(tree_branch &src_branch, int j, float rotate_start, float temp_deg, int branch_num) {
	
	int index(0);
	bool branch_deflected(false);
	tree_branch &branch(branches_34[1][branch_num]);
	tree_cylin &cylin(branch.cylin[0]);
	setup_rotate(cylin.rotate, rotate_start, temp_deg);
	tree_cylin &src_cylin(src_branch.cylin[j]);
	branch.num_cylins = num_4_cylins;
	float const radius1(src_branch.cylin[int(0.65*src_branch.num_cylins)].r2*0.9);
	cylin.assign_params(4, branch_num, radius1, radius1*gen_bc_size2(branch_4_rad_var), branch_4_length,
		src_cylin.deg_rotate + ((src_cylin.deg_rotate > 0.0) ? 1 : -1)*rand_gen(0,60));
	rotate_around_axis(src_cylin);
	add_rotation(cylin.p1, src_cylin.p1, BRANCH_LEN_SCALE);
	if (branch.num_cylins < 2) rotate_cylin(cylin);
	
	for (index = 1; index < branch.num_cylins; index++) {
		tree_cylin &cylin(branch.cylin[index]), &lcylin(branch.cylin[index-1]);
		gen_next_cylin(cylin, lcylin, branch_4_var, branch_4_rad_var, 4, branch_num, (index < branch.num_cylins-1));
		int const deg_added(generate_next_cylin(index, branch_curveness, branch.num_cylins, false, branch_deflected));
		gen_cylin_rotate(cylin.rotate, lcylin.rotate, deg_added*rotate_factor);
		cylin.rotate.z = lcylin.rotate.z;
		cylin.deg_rotate += ((branch.cylin[0].deg_rotate < 0.0) ? deg_added : -deg_added);
		add_rotation(lcylin.p2, lcylin.p1, 1.0);
		if (index == (branch.num_cylins - 1)) rotate_cylin(cylin);
	}
}


void tree_leaf::create_init_color(bool deterministic) {

	color  = 1.0;
	lred   = (deterministic ? rand2d() : rand_float());
	lgreen = (deterministic ? rand2d() : rand_float());
}


void process_leaf(vector<tree_leaf> &leaves, point start, point rotate, float b_deg_rot, float tsize) {

	tree_leaf leaf;
	int const val(rand_gen(0,60));
	float const deg_rotate(b_deg_rot + ((b_deg_rot > 0.0) ? val : -val));
	float const lsize(tsize*(0.7*rand2d() + 0.3));
	leaf.create_init_color(1);

	for (int p = 0; p < 4; ++p) {
		point lpts(leaf_points[p]*lsize);
		rotate_pts_around_axis(lpts, rotate, deg_rotate);
		add_rotation(leaf.pts[p], start, 1.0);
	}
	point lpts(0.0, 1.0, 0.0); // normal starts off in y-direction
	rotate_pts_around_axis(lpts, rotate, deg_rotate);
	leaf.norm.assign(re_matrix[0], re_matrix[1], re_matrix[2]);
	leaf.norm.normalize(); // should already be normalized
	leaves.push_back(leaf);
}


void tree::add_leaves_to_cylin(tree_cylin const &cylin, float tsize) {

	static float acc(0.0);
	point start;
	vector3d rotate;
	acc += num_leaves_per_occ;
	int const temp((int)acc);
	acc -= (float)temp;
	float const temp_deg(((cylin.rotate.y < 0.0) ? -1.0 : 1.0)*safe_acosf(cylin.rotate.x));
	float rotate_start(0.0);

	for (int l = 0; l < temp; l++) {
		if (init_deadness > 0 && init_deadness > rand_float2()) continue;
		rotate_start += 360.0/temp*l + rand_gen(1,30);
		if (rotate_start > 360.0) rotate_start -= 360.0;
		setup_rotate(rotate, rotate_start, temp_deg);
		rotate_around_axis(cylin);
		add_rotation(start, cylin.p1, 0.9);
		process_leaf(leaves, start, rotate, cylin.deg_rotate, tsize);
	}
}


int generate_next_cylin(int cylin_num, float branch_curveness, int ncib, bool branch_just_created, bool &branch_deflected) {

	//vars used in generating deg_rotate for cylinders
	float const PI_16(PI/16.0);
	float const t_start(TWO_PI/rand_gen(3,8)); //start in the first pi
	float const t_end(((rand_gen(1,3) == 1) ? 1.0 : 5.0)*PI_TWO + rand_gen(2,8)*PI_16); //either PI/2 to PI or 5*PI/2 to 3*PI - controls branch droopiness
	int add_deg_rotate(int(0.01*cylin_num*sinf(t_start+(t_end-t_start)*cylin_num/ncib)*branch_curveness)); //how much wavy the branch will be --in degrees
	int rg[2];

	if (cylin_num < int(ncib/3)) { //how much to start the starting deg_scale
		rg[0] = 5; rg[1] = 10;
	}
	else if (cylin_num < int(ncib*2/3)) { //scale for trig function -- to scale the middle of a branch
		rg[0] = 5; rg[1] = 10;
	}
	else { //how much to end the ending deg_scale
		rg[0] = 1; rg[1] = 5;
	}
	add_deg_rotate  *= rand_gen(rg[0], rg[1]);
	branch_deflected = false;

	if (branch_just_created) {
		if (rand_gen(0,10) < 6) return add_deg_rotate;  //don't deflect
		branch_deflected = true; //deflect
		return -add_deg_rotate;
	}
	return add_deg_rotate;
}


void regen_trees(bool recalc_shadows, bool keep_old) {

	cout << "vegetation: " << vegetation << endl;
	RESET_TIME;
	float const min_tree_h(island ? TREE_MIN_H : (water_plane_z + 0.01*zmax_est));
	float const max_tree_h(island ? TREE_MAX_H : 1.8*zmax_est);
	int const ext_x1(get_ext_x1()), ext_x2(get_ext_x2()), ext_y1(get_ext_y1()), ext_y2(get_ext_y2());
	static int init(0), last_rgi(0), last_xoff2(0), last_yoff2(0);
	static float last_ms(0.0), last_ms2(0.0);
	if (tree_mode && recalc_shadows) reset_shadows(OBJECT_SHADOW);
	
	if (tree_mode & 2) {
		gen_small_trees();
	}
	else {
		//remove_small_tree_cobjs();
	}
	if ((tree_mode & 1) && num_trees > 0) {
		if (keep_old && init && last_rgi == rand_gen_index && last_xoff2 == xoff2 && last_yoff2 == yoff2 &&
			last_ms == mesh_scale && last_ms2 == mesh_scale2)
		{ // keep old trees
			if (recalc_shadows) calc_visibility(SUN_SHADOW | MOON_SHADOW | TREE_ONLY);
			add_tree_cobjs();
			PRINT_TIME(" gen tree fast");
			return;
		}
		unsigned nkeep(0);

		if (scrolling) {
			vector3d const vd(-dx_scroll*DX_VAL, -dy_scroll*DY_VAL, 0.0);

			for (unsigned i = 0; i < t_trees.size(); ++i) { // keep any tree that's still in the scene
				point const gen_pos(t_trees[i].get_center());
				int const xp(get_xpos(gen_pos.x) - dx_scroll), yp(get_ypos(gen_pos.y) - dy_scroll);
				bool const keep(xp >= ext_x1 && xp <= ext_x2 && yp >= ext_y1 && yp <= ext_y2);
				t_trees[i].set_no_delete(keep);

				if (keep) { // shift it - don't have to recreate it
					t_trees[i].shift_tree(vd);
					++nkeep;
				}
			}
		}
		delete_trees();

		if (nkeep > 0) {
			for (unsigned i = 0; i < t_trees.size(); ++i) {
				if (t_trees[i].get_no_delete()) {
					t_trees[i].add_tree_collision_objects();
					t_trees[i].set_no_delete(0);
				}
				else { // remove this tree from the vector
					swap(t_trees[i], t_trees.back()); // deep copy = bad?
					t_trees.pop_back();
					--i;
				}
			}
		}
		else {
			max_leaves = 0;
			t_trees.resize(0);
		}
		PRINT_TIME(" Delete Trees");
		unsigned const smod(3.321*XY_MULT_SIZE+1), tree_prob(max(1, XY_MULT_SIZE/num_trees));
		unsigned const skip_val(max(1, int(1.0/sqrt((mesh_scale*mesh_scale2))))); // similar to deterministic gen in scenery.cpp

		for (int i = ext_y1; i < ext_y2; i += skip_val) {
			for (int j = ext_x1; j < ext_x2; j += skip_val) {
				if (scrolling) {
					int const ox(j + dx_scroll), oy(i + dy_scroll); // positions in original coordinate system
					if (ox >= ext_x1 && ox <= ext_x2 && oy >= ext_y1 && oy <= ext_y2) continue; // use orignal tree from last position
				}
				global_rand_gen.rseed1 = 805306457*(i + yoff2) + 12582917*(j + xoff2) + 100663319*rand_gen_index;
				global_rand_gen.rseed2 = 6291469  *(j + xoff2) + 3145739 *(i + yoff2) + 1572869  *rand_gen_index;
				rand2_mix();
				unsigned const val(((unsigned)rand2_seed_mix())%smod);
				if (val <= 100)         continue; // scenery
				if (val%tree_prob != 0) continue; // not selected
				if ((global_rand_gen.rseed1&127)/128.0 >= vegetation) continue;
				point pos((get_xval(j) + 0.5*DX_VAL*rand2d()), (get_yval(i) + 0.5*DY_VAL*rand2d()), 0.0);
				// Note: pos.z will be slightly different when calculated within vs. outside the mesh bounds
				pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
				if (pos.z > max_tree_h || pos.z < min_tree_h) continue;
				if (tree_mode == 3 && get_tree_type_from_height(pos.z) != 2) continue; // use a small (simple) tree here
				
				if (max_unique_trees > 0) {
					global_rand_gen.rseed1 = global_rand_gen.rseed1 % max_unique_trees;
					global_rand_gen.rseed2 = 12345;

					if (tree_coll_level == 0) {
						// FIXME: unique the trees so they can be reused
					}
				}
				t_trees.push_back(tree());
				t_trees.back().regen_tree(pos, 0); // use random function #2 for trees
			}
		}
		if (!scrolling) cout << "Num trees = " << t_trees.size() << endl;
		last_rgi   = rand_gen_index;
		last_xoff2 = xoff2;
		last_yoff2 = yoff2;
		last_ms    = mesh_scale;
		last_ms2   = mesh_scale2;
		init       = 1;
	}
	if (recalc_shadows) calc_visibility(SUN_SHADOW | MOON_SHADOW | TREE_ONLY);
	PRINT_TIME(" Gen Trees");
}


void tree::regen_tree(point &pos, int recalc_shadows) {

	gen_tree(pos, 0, -1, 1, 1);
	if (recalc_shadows) gen_tree_shadows((SUN_SHADOW | MOON_SHADOW));
}


void tree::shift_tree(vector3d const &vd) { // tree has no single pos, have to translate every coordinate!!!

	sphere_center += vd;

	for (unsigned j = 0; j < all_cylins.size(); ++j) {
		all_cylins[j].p1 += vd;
		all_cylins[j].p2 += vd;
	}
	for (unsigned j = 0; j < leaves.size(); ++j) {
		for (unsigned k = 0; k < 4; ++k) {
			leaves[j].pts[k] += vd;
		}
	}
	for (unsigned i = 0; i < leaf_data.size(); ++i) {
		leaf_data[i].v += vd;
	}
	for (unsigned i = 0; i < leaf_data2.size(); ++i) {
		leaf_data2[i].v += vd;
	}
}


void shift_trees(vector3d const &vd) {

	if (num_trees > 0) return; // dynamically created, not placed

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].shift_tree(vd);
	}
}


void add_tree_cobjs() {

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].add_tree_collision_objects();
	}
}


void clear_tree_vbos() {

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].clear_vbo();
	}
}




