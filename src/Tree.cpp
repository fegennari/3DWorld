// 3D World - OpenGL CS184 Computer Graphics Project
// by Hiral Patel and Frank Gennari
// 3/15/02

#include "GL/glew.h" // must be included first
#include "3DWorld.h"
#include "mesh.h"
#include "tree_3dw.h"
#include "tree_leaf.h"
#include "explosion.h"
#include "physics_objects.h"
#include "gl_ext_arb.h"


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
float const DIST_C_SCALE     = 0.25;
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
bool const USE_VBOS          = 1;


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


extern bool has_snow;
extern int shadow_detail, island, num_trees, do_zoom, begin_motion, display_mode, animate2, iticks, draw_model;
extern int xoff2, yoff2, rand_gen_index, gm_blast, game_mode, leaf_color_changed, scrolling, dx_scroll, dy_scroll;
extern long rseed1, rseed2;
extern float zmin, zmax_est, zbottom, water_plane_z, mesh_scale, mesh_scale2, temperature, fticks, tree_size, vegetation;
extern point ocean;
extern lightning l_strike;
extern blastr latest_blastr;
extern texture textures[];
extern vector<coll_obj> coll_objects;


void calc_leaf_points() {

	leaf_size = REL_LEAF_SIZE*TREE_SIZE/sqrt(nleaves_scale);
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


void tree::gen_tree_shadows(char light_sources, int index) {

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


void tree::add_tree_collision_objects(int ix) {

	//RESET_TIME;
	if (!(tree_mode & 1) || !tree_coll_level || !created) return;
	remove_collision_objects();
	if (!is_over_mesh()) return; // optimization
	int const btid(tree_types[type].bark_tex);
	colorRGBA const bcolor(tree_types[type].barkc);
	cobj_params cp(0.8, bcolor, TEST_RTREE_COBJS, 0, NULL, ix, btid, 4.0, 1, 0);
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
		cobj_params cpl(0.3, lcolor, TEST_RTREE_COBJS, 0, NULL, ix, (TEST_RTREE_COBJS ? -1 : ltid), 1.0, 0, 0);
		cpl.shadow = 0;
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


void remove_tree_cobjs(tree_cont_t &t_trees) {

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].remove_collision_objects();
	}
}


void draw_trees_bl(tree_cont_t &ts, bool lpos_change, bool draw_branches, bool near_draw_leaves, bool draw_far_leaves) {

	for (unsigned i = 0; i < ts.size(); ++i) {
		ts[i].draw_tree(lpos_change, draw_branches, near_draw_leaves, draw_far_leaves);
	}
}


void set_leaf_shader(float min_alpha, bool use_wind) {

	set_shader_prefix("#define USE_LIGHT_COLORS", 0); // VS
	setup_enabled_lights(2);
	unsigned p(0);
	
	if (use_wind) { // wind on leaves
		// FIXME: leaves are drawn as quads, not triangles
		p = set_shader_prog("ads_lighting.part*+tree_leaves", "linear_fog.part+simple_texture", "wind.part*+tri_wind", GL_TRIANGLES, GL_TRIANGLE_STRIP, 3);
		setup_wind_for_shader(p);
	}
	else {
		p = set_shader_prog("ads_lighting.part*+tree_leaves", "linear_fog.part+simple_texture");
	}
	setup_fog_scale(p);
	add_uniform_float(p, "min_alpha", min_alpha);
	set_multitex(0);
	add_uniform_int(p, "tex0", 0);
}


void draw_trees(tree_cont_t &ts) {

	//glFinish(); // testing
	//RESET_TIME;

	if (tree_mode & 2) {
		draw_small_trees();
		//PRINT_TIME("Small Trees");
	}
	if (tree_mode & 1) {
		static point last_lpos(all_zeros);
		point const lpos(get_light_pos());
		bool const lpos_change(lpos != last_lpos);
		if (lpos_change) update_cobj_tree();

		// draw branches, then leaves: much faster for distant trees, slightly slower for near trees
		colorRGBA const orig_fog_color(setup_smoke_shaders(0.0, 0, 0, 0, 1, 1, 0)); // dynamic lights, but no smoke (yet)
		draw_trees_bl(ts, lpos_change, 1, 0, 0); // branches
		end_smoke_shaders(orig_fog_color);
		set_leaf_shader(0.75, 0);
		draw_trees_bl(ts, lpos_change, 0, 0, 1); // far  leaves
		//unset_shader_prog();
		//set_leaf_shader(0.75, 1);
		draw_trees_bl(ts, lpos_change, 0, 1, 0); // near leaves
		unset_shader_prog();
		last_lpos = lpos;
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

	assert(i < leaf_data.size());
	
	for (unsigned j = 0; j < 4; ++j) {
		leaf_data[j+(i<<2)].set_c3(color);
	}
	if (i < leaves.size() && leaves[i].coll_index >= 0) { // update cobj color so that leaf water reflection is correct
		assert((unsigned)leaves[i].coll_index < coll_objects.size());
		coll_objects[leaves[i].coll_index].cp.color = colorRGBA(color).modulate_with(texture_color(tree_types[type].leaf_tex));
	}
	mark_leaf_changed(i);
}


void tree::change_leaf_color(colorRGBA &base_color, unsigned i) {

	if (leaf_data.empty()) return;
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
		vector<int> indices; // unused
		coll_objects[cix].update_shadowed_cobjs(coll_objects, indices, cix);
		remove_coll_object(cix);
	}
	leaves[i] = leaves.back();
	leaves.pop_back();

	if (update_data && !leaf_data.empty()) {
		unsigned const i4(i << 2), tnl4(leaves.size() << 2);
		assert(4*leaves.size() <= leaf_data.size());

		for (unsigned j = 0; j < 4; ++j) { // shift vertex array (last one is now invalid)
			leaf_data[j+i4] = leaf_data[j+tnl4];
		}
		mark_leaf_changed(i);
	}
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
		
		if (!leaf_data.empty()) { // remove leaf i
			remove_leaf(i, 1);
			return 1;
		}
	}
	else {
		float const llc(lcolor);
		lcolor = max(0.0, (lcolor - 0.3*damage_done));
		if (lcolor == 0.0 && llc > 0.0) damage += damage_scale;
		if (!leaf_data.empty()) change_leaf_color(base_color, i);
	}
	return 0;
}


void tree::blast_damage(blastr const *const blast_radius) {

	assert(blast_radius);
	float const bradius(blast_radius->cur_size), bdamage(LEAF_DAM_SCALE*blast_radius->damage);
	if (bdamage == 0.0) return;
	point const &blast_pos(blast_radius->pos);
	if (p2p_dist_sq(blast_pos, sphere_center) > (bradius + sphere_radius)*(bradius + sphere_radius)) return;
	float const bradius_sq(bradius*bradius);

	for (unsigned i = 0; i < leaves.size(); i++) {
		if (leaves[i].color < 0.0) continue;
		float const dist(p2p_dist_sq(blast_pos, leaves[i].pts[0]));
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
	unsigned const nleaves(leaves.size());
	float const temp0(max(1.0f, min(0.3f, (20.0f-temperature)/30.0f)));
	int const rgen(min(LEAF_GEN_RAND2/10, int(rand_uniform(0.5, 1.5)*temp0*LEAF_GEN_RAND2/fticks)));

	for (unsigned i = (rand()%LEAF_GEN_RAND1); i < nleaves; i += LEAF_GEN_RAND1) {
		if ((rand()%rgen) == 0) {
			gen_leaf_at(leaves[i].pts, leaves[i].norm, type, get_leaf_color(i));
			// create a new leaf with a different color (and orient?)
			leaves[i].create_init_color(0);
			if (!leaf_data.empty()) copy_color(leaves[i].calc_leaf_color(leaf_color, base_color), i);
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
	if (leaf_data.empty()) return leaves[i].calc_leaf_color(leaf_color, base_color);
	assert(4*i < leaf_data.size());
	return leaf_data[i<<2].get_c3(); // return color of first vertex since they all should be the same
}


inline colorRGB tree_leaf::calc_leaf_color(colorRGBA const &leaf_color, colorRGBA const &base_color) const {

	float const ilch(1.0 - leaf_color_coherence);
	return colorRGB(color*(leaf_color.red   + ilch*lred)   + base_color.red,
		            color*(leaf_color.green + ilch*lgreen) + base_color.green,
					0.0);
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


void tree::draw_tree(bool invalidate_norms, bool draw_branches, bool draw_near_leaves, bool draw_far_leaves) {

	if (!created) return;
	rseed1 = trseed1;
	rseed2 = trseed2;
	gen_leaf_color();

	float const mscale(mesh_scale*mesh_scale2);
	float const dist_c(mscale*DIST_C_SCALE/(do_zoom ? ZOOM_FACTOR : 1.0));
	float const dist_cs(distance_to_camera(sphere_center)*dist_c*(0.04/base_radius));
	int const leaf_dynamic_en(!has_snow && (display_mode & 0x0100) != 0);
	int const leaf_dynamic(leaf_dynamic_en && mscale/dist_cs > 0.5);
	bool const draw_leaves(leaf_dynamic ? draw_near_leaves : draw_far_leaves);
	
	if (draw_leaves && deadness < 1.0 && init_deadness < 1.0) {
		burn_leaves();
		if (l_strike.enabled == 1)    lightning_damage(l_strike.end);
		if (game_mode && gm_blast)    blast_damage(&latest_blastr);
		if (begin_motion && animate2) drop_leaves();
	}
	if (TEST_RTREE_COBJS) return;

	if (!draw_leaves || draw_branches) { // second pass only
		int const level((island || get_camera_pos().z > ztop) ? 1 : 2); // do we want to test the mesh in non-island mode?
		not_visible = !sphere_in_camera_view(sphere_center, 1.1*sphere_radius, level);
	}
	if (not_visible) return;
	bool const use_vbos(USE_VBOS && setup_gen_buffers());
	if (draw_branches) draw_tree_branches(mscale, dist_c, dist_cs, use_vbos);
	if (leaves.empty() || deadness >= 1.0 || init_deadness >= 1.0) return;
	if (draw_leaves  ) draw_tree_leaves(invalidate_norms, mscale, dist_cs, use_vbos, (leaf_dynamic + leaf_dynamic_en));
}


void tree::draw_tree_branches(float mscale, float dist_c, float dist_cs, bool use_vbos) {

	point const camera(get_camera_pos());
	assert(quadric);
	set_fill_mode();
	select_texture(tree_types[type].bark_tex);
	set_color(bcolor);
	BLACK.do_glColor();

	if (use_vbos) { // draw with branch vbos
		size_t const branch_stride(sizeof(vert_norm_tc));

		if (branch_vbo == 0) { // create vbo
			assert(branch_ivbo == 0);
			unsigned const numcylin(all_cylins.size());
			assert(num_branch_quads == 0 && num_unique_pts == 0);

			for (unsigned i = 0; i < numcylin; i++) { // determine required data size
				bool const prev_connect(i > 0 && all_cylins[i].can_merge(all_cylins[i-1]));
				unsigned const ndiv(all_cylins[i].get_num_div());
				num_branch_quads += ndiv;
				num_unique_pts   += (prev_connect ? 1 : 2)*ndiv;
			}
			typedef unsigned short index_t;
			assert(num_unique_pts < (1 << 8*sizeof(index_t))); // cutting it close with 4th order branches
			vector<vert_norm_tc> data;
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
					data_pos = data.size();
					quad_id  = cylin_id = 0;
				}
				for (unsigned j = prev_connect; j < 2; ++j) { // create vertex data
					for (unsigned S = 0; S < ndiv; ++S) { // first cylin: 0,1 ; other cylins: 1
						float const tx(2.0*fabs(S*ndiv_inv - 0.5));
						// FIXME: something is still wrong - twisted branch segments due to misaligned or reversed starting points
						data.push_back(vert_norm_tc(vpn.p[(S<<1)+j], vpn.n[S], tx, float(cylin_id + j))); // average normals?
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
			upload_vbo_data(&data.front(),  data.size()*branch_stride,    0); // ~350KB
			upload_vbo_data(&idata.front(), idata.size()*sizeof(index_t), 1); // ~75KB (with 16-bit index)
		} // end create vbo
		assert(branch_vbo > 0 && branch_ivbo > 0);
		bind_vbo(branch_vbo,  0); // use vbo for rendering
		bind_vbo(branch_ivbo, 1);
		set_array_client_state(1, 1, 1, 0);
		glVertexPointer(  3, GL_FLOAT, branch_stride, 0);
		glNormalPointer(     GL_FLOAT, branch_stride, (void *)(sizeof(point)));
		glTexCoordPointer(2, GL_FLOAT, branch_stride, (void *)(sizeof(point) + sizeof(vector3d)));
		unsigned const num(4*min(num_branch_quads, max((num_branch_quads/8), unsigned(1.5*num_branch_quads*mscale/dist_cs)))); // branch LOD
		glDrawRangeElements(GL_QUADS, 0, num_unique_pts, num, GL_UNSIGNED_SHORT, 0);
		bind_vbo(0, 0);
		bind_vbo(0, 1);
	}
	else { // draw branches
		float const bs_scale(350.0*mscale/dist_c);
		gluQuadricTexture(quadric, GL_TRUE);

		for (unsigned i = 0; i < all_cylins.size(); i++) {
			draw_cylin const &cylin(all_cylins[i]);
			float const bsize((cylin.r1 + cylin.r2)*InvSqrt(p2p_dist_sq(camera, cylin.p1)));
			unsigned ndiv(min(cylin.get_num_div(), unsigned(bs_scale*bsize)));
			if (cylin.level == 0) ndiv = max(ndiv, 3U);
			if (ndiv == 0) continue;

			if (ndiv <= 3) { // draw as a quad
				int npts(0);
				point pts[4];
				vector3d const t_to_c((camera - cylin.p1).get_norm());
				vector3d view_dir(t_to_c), v2;
				orthogonalize_dir(view_dir, (cylin.p2 - cylin.p1), view_dir, 0);
				view_dir.do_glNormal();
				cylinder_quad_projection(pts, cylin, t_to_c, v2, npts);
				glBegin((npts == 4) ? GL_QUADS : GL_TRIANGLES);

				for (int i = 0; i < npts; ++i) {
					glTexCoord2f((i&1), (i>>1)); pts[i].do_glVertex(); // textures?
				}
				glEnd();
			}
			else { // draw as a cylinder
				draw_fast_cylinder(cylin.p1, cylin.p2, cylin.r1, cylin.r2, ndiv, 1); // slow, but shouldn't be used anymore
			}
		} // for i
		gluQuadricTexture(quadric, GL_FALSE);
	}
	glDisable(GL_TEXTURE_2D);
}


void tree::draw_tree_leaves(bool invalidate_norms, float mscale, float dist_cs, bool use_vbos, int leaf_dynamic) {

	unsigned nleaves(leaves.size());
	assert(nleaves <= max_leaves);
	bool const gen_arrays(leaf_data.empty());
	unsigned const leaf_stride(sizeof(vert_norm_tc_color));

	if (gen_arrays) { // extra memory = 176*nleaves
		if (deadness > 0) {
			for (unsigned i = 0; i < nleaves; i++) { // process deadness stuff
				if (deadness > rand_float2()) {
					remove_leaf(i, 0);
					--nleaves;
				}
			}
		}
		leaf_data.resize(4*nleaves);
		gen_quad_tex_coords(&(leaf_data.front().t[0]), nleaves, leaf_stride/sizeof(float));
	}
	unsigned nl(nleaves);
	if (ENABLE_CLIP_LEAVES) nl = min(nl, max((nl/8), unsigned((use_vbos ? 4.0 : 2.0)*nl*mscale/dist_cs))); // leaf LOD
	
	if (gen_arrays || (reset_leaves && !leaf_dynamic)) {
		for (unsigned i = 0; i < nleaves; i++) { // process leaf points - reset to default positions and normals
			for (unsigned j = 0; j < 4; ++j) {
				leaf_data[j+(i<<2)].v = leaves[i].pts[j];
				leaf_data[j+(i<<2)].n = leaves[i].norm*leaves[i].get_norm_scale(j);
			}
		}
		reset_leaves   = 0;
		leaves_changed = 1;
	}
	if (leaf_dynamic > 1) { // leaves move in wind or when struck by an object (somewhat slow)
		int last_xpos(0), last_ypos(0);
		vector3d local_wind;

		for (unsigned i = 0; i < nl; i++) { // process leaf wind and collisions
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
						if (damage_leaf(i, BEAM_DAMAGE)) {--i; --nl; --nleaves;}
					}
					else {
						if (coll_type == PROJECTILE) {
							if ((rand()&3) == 0) { // shoot off leaf
								gen_leaf_at(l.pts, l.norm, type, get_leaf_color(i));
								remove_leaf(i, 1);
								--i; --nl; --nleaves;
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
			leaf_data[(i<<2)+1].v = l.pts[1] + delta;
			leaf_data[(i<<2)+2].v = l.pts[2] + delta;
			vector3d normal(cross_product(new_dir, (l.pts[3] - l.pts[0])).get_norm());

			for (unsigned j = 0; j < 4; ++j) { // update the normals, even though this slows the algorithm down
				leaf_data[j+(i<<2)].n = normal*l.get_norm_scale(j);
			}
			if (LEAF_HEAL_RATE > 0.0 && l.color > 0.0 && l.color < 1.0) { // leaf heal
				leaves[i].color = min(1.0f, (l.color + LEAF_HEAL_RATE*fticks));
				copy_color(l.calc_leaf_color(leaf_color, base_color), i);
			}
			reset_leaves = 1; // Do we want to update the normals and collision objects as well?
			mark_leaf_changed(i);
		} // for i
	}
	if (gen_arrays || leaf_color_changed) { // process leaf colors
		for (unsigned i = 0; i < nleaves; i++) {
			copy_color(leaves[i].calc_leaf_color(leaf_color, base_color), i);
		}
	}
	if (gen_arrays || (invalidate_norms && tree_coll_level >= 4)) { // process leaf shadows/normals
		int const light(get_light());

		for (unsigned i = 0; i < nleaves; i++) {
			tree_leaf &l(leaves[i]);
			l.shadow_bits = 0;

			for (unsigned j = 0; j < 4; ++j) {
				bool const shadowed(l.coll_index >= 0 && !is_visible_to_light_cobj(l.pts[j], light, 0.0, l.coll_index, 1));
				l.shadow_bits |= (int(shadowed) << j);
				leaf_data[j+(i<<2)].n = l.norm*l.get_norm_scale(j);
			}
		}
		leaves_changed = 1;
	}
	assert(leaf_data.size() >= 4*leaves.size());
	bool const draw_as_points(0); // testing

	if (use_vbos) {
		if (leaf_vbo == 0) {
			leaf_vbo = create_vbo();
			assert(leaf_vbo > 0);
			bind_vbo(leaf_vbo);
			upload_vbo_data(&leaf_data.front(), leaf_data.size()*leaf_stride); // ~150KB
		}
		else {
			bind_vbo(leaf_vbo);
			if (leaves_changed) upload_vbo_sub_data(&leaf_data.front(), 0, leaf_data.size()*leaf_stride);
		}
		vert_norm_tc_color::set_vbo_arrays(draw_as_points ? 4 : 1);
	}
	else {
		leaf_data.front().set_state(draw_as_points ? 4 : 1);
	}
	if (draw_model == 0) { // solid fill
		enable_blend();
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.75);
	}
	unsigned const num_dlights(enable_dynamic_lights(sphere_center, sphere_radius));
	add_uniform_int(0, "num_dlights", num_dlights);
	set_lighted_sides(2);
	set_specular(0.1, 10.0);
	set_fill_mode();
	if (!draw_as_points) select_texture((draw_model == 0) ? tree_types[type].leaf_tex : WHITE_TEX); // what about texture color mod?
	glEnable(GL_COLOR_MATERIAL);
	glDisable(GL_NORMALIZE);
	glDrawArrays((draw_as_points ? GL_POINTS : GL_QUADS), 0, (draw_as_points ? 1 : 4)*nl);
	glDisable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	glDisable(GL_TEXTURE_2D);
	disable_blend();
	set_specular(0.0, 1.0);
	glDisable(GL_ALPHA_TEST);
	set_lighted_sides(1);
	disable_dynamic_lights(num_dlights);
	if (leaf_vbo > 0) bind_vbo(0);
	leaves_changed = 0;
}


void delete_trees(tree_cont_t &ts) {

	unsigned deleted(0);

	for (unsigned i = ts.size(); i > 0; --i) { // delete backwards (pop collision stack)
		if (ts[i-1].delete_tree()) ++deleted;
	}
	if (tree_coll_level && !ts.empty()) purge_coll_freed(1); // MUST do this for collision detection to work
}


int tree::delete_tree() {

	clear_vbo();
	if (!created)  return 0;
	if (tree_coll_level) remove_collision_objects();
	if (no_delete) return 0;
	all_cylins.clear();
	leaf_data.clear();
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



void tree::gen_tree(point &pos, int size, int ttype, int calc_z, bool add_cobjs, int ix) {

	sphere_center = pos; // z value will be reset later
	if (calc_z) pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
	leaf_data.clear();

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
	if (add_cobjs) add_tree_collision_objects(ix);
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


void regen_trees(tree_cont_t &t_trees, bool recalc_shadows, bool keep_old) {

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
			add_tree_cobjs(t_trees);
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
		delete_trees(t_trees);

		if (nkeep > 0) {
			for (unsigned i = 0; i < t_trees.size(); ++i) {
				if (t_trees[i].get_no_delete()) {
					t_trees[i].add_tree_collision_objects(i);
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
				rseed1 = 805306457*(i + yoff2) + 12582917*(j + xoff2) + 100663319*rand_gen_index;
				rseed2 = 6291469  *(j + xoff2) + 3145739 *(i + yoff2) + 1572869  *rand_gen_index;
				rand2_mix();
				unsigned const val(((unsigned)rand2_seed_mix())%smod);
				if (val <= 100)         continue; // scenery
				if (val%tree_prob != 0) continue; // not selected
				if ((rseed1&127)/128.0 >= vegetation) continue;
				point pos((get_xval(j) + 0.5*DX_VAL*rand2d()), (get_yval(i) + 0.5*DY_VAL*rand2d()), 0.0);
				// Note: pos.z will be slightly different when calculated within vs. outside the mesh bounds
				pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
				if (pos.z > max_tree_h || pos.z < min_tree_h) continue;
				if (tree_mode == 3 && get_tree_type_from_height(pos.z) != 2) continue; // use a small (simple) tree here
				
				if (max_unique_trees > 0) {
					rseed1 = rseed1 % max_unique_trees;
					rseed2 = 12345;
					// *** unique the trees so they can be reused? ***
				}
				t_trees.push_back(tree());
				t_trees.back().regen_tree(pos, 0, t_trees.size()-1); // use random function #2 for trees
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


void tree::regen_tree(point &pos, int recalc_shadows, int index) {

	gen_tree(pos, 0, -1, 1, 1, index);
	if (recalc_shadows) gen_tree_shadows((SUN_SHADOW | MOON_SHADOW), index);
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
}


void shift_trees(tree_cont_t &t_trees, vector3d const &vd) {

	if (num_trees > 0) return; // dynamically created, not placed

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].shift_tree(vd);
	}
}


void add_tree_cobjs(tree_cont_t &t_trees) {

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].add_tree_collision_objects(i);
	}
}


void clear_tree_vbos(tree_cont_t &t_trees) {

	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].clear_vbo();
	}
}




