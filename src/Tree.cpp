// 3D World - OpenGL CS184 Computer Graphics Project
// by Hiral Patel and Frank Gennari
// 3/15/02

#include "3DWorld.h"
#include "mesh.h"
#include "tree_3dw.h"
#include "tree_leaf.h"
#include "explosion.h"
#include "physics_objects.h"
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
bool const FORCE_TREE_TYPE   = 1;
unsigned const CYLINS_PER_ROOT = 3;


vector<tree_cylin >   tree_builder_t::cylin_cache;
vector<tree_branch>   tree_builder_t::branch_cache;
vector<tree_branch *> tree_builder_t::branch_ptr_cache;


// tree helper methods
void rotate_all(point const &rotate, float angle, float x, float y, float z);


// tree_mode: 0 = no trees, 1 = large only, 2 = small only, 3 = both large and small
unsigned max_unique_trees(0);
int tree_mode(1), tree_coll_level(TREE_COLL);
float tree_temp_matrix[3], i_matrix[9], re_matrix[3];
float leaf_color_coherence(0.5), tree_color_coherence(0.2), tree_deadness(-1.0), nleaves_scale(1.0), leaf_size(0.0);
colorRGBA leaf_base_color(BLACK);
point leaf_points[4]; // z = 0.0 -> -0.05
tree_data_manager_t tree_data_manager;
tree_cont_t t_trees(tree_data_manager);


extern bool has_snow, no_sun_lpos_update, has_dl_sources, gen_tree_roots;
extern int shadow_detail, island, num_trees, do_zoom, begin_motion, display_mode, animate2, iticks, draw_model, frame_counter;
extern int xoff2, yoff2, rand_gen_index, gm_blast, game_mode, leaf_color_changed, scrolling, dx_scroll, dy_scroll, window_width, window_height;
extern float zmin, zmax_est, zbottom, water_plane_z, tree_scale, temperature, fticks, vegetation;
extern point ocean;
extern lightning l_strike;
extern blastr latest_blastr;
extern texture_t textures[];
extern coll_obj_group coll_objects;
extern rand_gen_t global_rand_gen;


struct texture_pair_t {

	unsigned tids[2]; // color, normal

	texture_pair_t() {tids[0] = tids[1] = 0;}

	void free_context() {
		for (unsigned d = 0; d < 2; ++d) {free_texture(tids[d]);}
	}
	void bind_textures() const {
		for (unsigned d = 0; d < 2; ++d) {
			assert(tids[d]);
			bind_2d_texture(tids[d]);
			set_multitex(d);
		}
		set_multitex(0);
	}
	static void ensure_tid(unsigned &tid, unsigned tsize) {
		if (tid) return; // already created
		setup_texture(tid, GL_MODULATE, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tsize, tsize, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	}
	void ensure_tids(unsigned tsize) {
		for (unsigned d = 0; d < 2; ++d) {ensure_tid(tids[d], tsize);}
	}
};


struct render_to_texture_t {

	unsigned tsize;

	render_to_texture_t(unsigned tsize_) : tsize(tsize_) {}
	virtual ~render_to_texture_t() {}

	void render(texture_pair_t &tpair, point const &eye, point const &center, float radius, vector3d const &up_dir) {
		assert(eye != center);
		assert(radius > 0.0);
		assert(tsize > 0);

		// setup matrices
		glViewport(0, 0, tsize, tsize);
		glClearColor(0.0, 0.0, 0.0, 0.0); // transparent
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
		glPushMatrix();
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		float const dist(p2p_dist(eye, center));
		assert(dist > radius);
		float const angle(asinf(radius/dist)), near_clip(dist - radius), far_clip(dist + radius);
		assert(near_clip > 0);
		gluPerspective(2.0*angle/TO_RADIANS, 1.0, near_clip, far_clip); // gluOrtho2D()?
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up_dir.x, up_dir.y, up_dir.z);

		// render
		tpair.ensure_tids(tsize);
		glDisable(GL_LIGHTING);

		for (unsigned d = 0; d < 2; ++d) {
			unsigned fbo_id(0);
			enable_fbo(fbo_id, tpair.tids[d], 0); // FIXME: too slow to create and free fbos every time?
			draw_geom(d != 0);
			free_fbo(fbo_id);
		}

		// restore state
		//glMatrixMode(GL_TEXTURE);
		//glLoadIdentity();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glEnable(GL_LIGHTING);
		disable_fbo();
		glViewport(0, 0, window_width, window_height);
	}

	virtual void draw_geom(bool is_normal_pass) = 0;
};


struct render_tree_to_texture_t : public render_to_texture_t {

	shader_t branch_shader[2], leaf_shader[2]; // color, normal
	tree *cur_tree;

	render_tree_to_texture_t(unsigned tsize_) : render_to_texture_t(tsize_), cur_tree(NULL) {}

	void free_context() {
		for (unsigned d = 0; d < 2; ++d) {
			branch_shader[d].end_shader();
			leaf_shader[d].end_shader();
		}
	}
	void set_leaf_shader_consts(bool ix) {
		leaf_shader[ix].begin_shader(1);
		leaf_shader[ix].add_uniform_float("min_alpha", 0.75);
		leaf_shader[ix].add_uniform_int("tex0", 0);
		leaf_shader[ix].add_uniform_int("tc_start_ix", 3);
		leaf_shader[ix].disable();
	}
	void setup_shaders() {
		if (!branch_shader[0].is_setup()) { // colors
			branch_shader[0].set_vert_shader("no_lighting_tex_coord");
			branch_shader[0].set_frag_shader("simple_texture");
			branch_shader[0].begin_shader(1);
			branch_shader[0].add_uniform_int("tex0", 0);
			branch_shader[0].disable();
		}
		if (!branch_shader[1].is_setup()) { // normals
			branch_shader[1].set_vert_shader("write_normal");
			branch_shader[1].set_frag_shader("write_normal");
			branch_shader[1].begin_shader(0); // begin, but don't enable
		}
		if (!leaf_shader[0].is_setup()) { // colors
			leaf_shader[0].set_vert_shader("tc_by_vert_id.part+tree_leaves_no_lighting");
			leaf_shader[0].set_frag_shader("simple_texture");
			set_leaf_shader_consts(0);
		}
		if (!leaf_shader[1].is_setup()) { // normals
			leaf_shader[1].set_vert_shader("tc_by_vert_id.part+tree_leaves_no_lighting");
			leaf_shader[1].set_frag_shader("write_normal_textured");
			set_leaf_shader_consts(1);
		}
	}
	virtual void draw_geom(bool is_normal_pass) {
		assert(cur_tree);
		BLACK.do_glColor();
		branch_shader[is_normal_pass].enable();
		cur_tree->draw_tree(branch_shader[is_normal_pass], 1, 0, is_normal_pass, zero_vector); // draw branches
		branch_shader[is_normal_pass].disable();
		disable_multitex_a();
		tree_cont_t::pre_leaf_draw(leaf_shader[is_normal_pass]);
		cur_tree->draw_tree(leaf_shader  [is_normal_pass], 0, 1, is_normal_pass, zero_vector); // draw leaves
		tree_cont_t::post_leaf_draw(leaf_shader[is_normal_pass]);
	}
	void render_tree(tree &t, texture_pair_t &tpair, vector3d const &view_dir, vector3d const &up_dir) {
		setup_shaders();
		cur_tree = &t;
		float const view_dist_mult = 20.0; // ???
		float const radius(t.get_radius());
		point const center(t.sphere_center());
		point const eye(center - view_dist_mult*radius*view_dir.get_norm());
		render(tpair, eye, center, radius, up_dir);
	}
	void render_tree_side_views(tree &t, vector<texture_pair_t> &tpairs, unsigned num) {
		assert(num > 0);

		for (unsigned i = 0; i < num; ++i) {
			float const angle(TWO_PI*((float)i/(float)num));
			vector3d const view_dir(cosf(angle), sinf(angle), 0.0); // should be normalized
			tpairs.push_back(texture_pair_t());
			render_tree(t, tpairs.back(), view_dir, plus_z);
		}
	}
	void render_tree_top_view(tree &t, texture_pair_t &tpair) {
		render_tree(t, tpair, vector3d(0.0, 0.0, -1.0), plus_x);
	}
};


void calc_leaf_points() {

	leaf_size = REL_LEAF_SIZE*TREE_SIZE/(sqrt(nleaves_scale)*tree_scale);
	leaf_points[0].assign(-2.0*leaf_size, 0.0, 0.0);
	leaf_points[1].assign(-2.0*leaf_size, 0.0, 4.0*leaf_size);
	leaf_points[2].assign( 2.0*leaf_size, 0.0, 4.0*leaf_size);
	leaf_points[3].assign( 2.0*leaf_size, 0.0, 0.0);
}


//designed to handle only positive numbers correctly
//return a random number between start and end
inline int rand_gen(int start, int end) {
	return (rand2()%(end - start + 1) + start);
}

float get_tree_z_bottom(float z, point const &pos) {
	return ((world_mode == WMODE_GROUND && is_over_mesh(pos)) ? max(zbottom, (z - TREE_DEPTH)) : (z - TREE_DEPTH));
}


colorRGBA get_leaf_base_color(int type) {

	colorRGBA color(tree_types[type].leafc);
	UNROLL_3X(color[i_] = CLIP_TO_01(color[i_] + leaf_base_color[i_]);)
	return color;
}


bool tree::is_over_mesh() const {

	//if (world_mode != WMODE_GROUND) return 1; // ???
	float const r(tdata().sphere_radius);
	int const x1(get_xpos(tree_center.x - r)), y1(get_ypos(tree_center.y - r));
	int const x2(get_xpos(tree_center.x + r)), y2(get_ypos(tree_center.y + r));
	return (x1 < MESH_X_SIZE && y1 < MESH_Y_SIZE && x2 >= 0 && y2 >= 0); // completely off the mesh
}


void tree::gen_tree_shadows(unsigned light_sources) {

	if (shadow_detail < 2 || !physics_enabled()) return;
	// Note: not entirely correct since an off mesh tree can still cast a shadow on the mesh
	if (!is_over_mesh()) return; // optimization
	tree_data_t const &td(tdata());
	if (!enable_shadow_envelope(sphere_center(), td.sphere_radius, light_sources, 1)) return;
	vector<draw_cylin> const &cylins(td.get_all_cylins());

	for (unsigned i = 0; i < cylins.size(); i++) {
		draw_cylin const &c(cylins[i]);
		if ((c.level + 2) > shadow_detail) return;
		cylinder_shadow(c.p1+tree_center, c.p2+tree_center, c.r1, c.r2, light_sources, 0, 0, (c.level < 2));
	}
	if (shadow_detail < 6) return;
	int const ltid(tree_types[type].leaf_tex);
	vector<tree_leaf> const &leaves(td.get_leaves());
	
	for (unsigned i = 0; i < leaves.size(); i++) { // loop through leaves
		point pts[4];
		get_abs_leaf_pts(pts, i);
		polygon_shadow(pts, leaves[i].norm, 4, 0.0, light_sources, 0, 0, 0, ltid);
	}
	disable_shadow_envelope(light_sources);
}


void tree::add_tree_collision_objects() {

	//RESET_TIME;
	if (!tree_coll_level || !physics_enabled()) return;
	remove_collision_objects();
	if (!is_over_mesh()) return; // optimization
	assert(type < NUM_TREE_TYPES);
	int const btid(tree_types[type].bark_tex), branch_coll_level(min(tree_coll_level, 4));
	cobj_params cp(0.8, tree_types[type].barkc, TEST_RTREE_COBJS, 0, NULL, 0, btid, 4.0, 1, 0);
	cp.shadow = 0; // will be handled by gen_tree_shadows()
	assert(branch_cobjs.empty());
	vector<draw_cylin> const &cylins(tdata().get_all_cylins());

	for (unsigned i = 0; i < cylins.size(); i++) {
		draw_cylin const &c(cylins[i]);

		if (c.level < branch_coll_level) {
			branch_cobjs.push_back(add_coll_cylinder(c.p1+tree_center, c.p2+tree_center, c.r1, c.r2, cp, -1, ((c.level == 0) ? 0 : 1)));
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
		vector<tree_leaf> const &leaves(tdata().get_leaves());
		leaf_cobjs.resize(leaves.size());

		for (unsigned i = 0; i < leaves.size(); i++) { // loop through leaves
			// Note: line collisions with leaves will use the texture alpha component for a more exact test
			point pts[4];
			get_abs_leaf_pts(pts, i);
			leaf_cobjs[i] = add_coll_polygon(pts, 4, cpl, 0.0, xlate, -1, 2);
			coll_objects[leaf_cobjs[i]].is_billboard = 1;
		}
	}
	//PRINT_TIME("Tree Cobjs");
}


void tree::remove_collision_objects() {

	for (unsigned i = 0; i < branch_cobjs.size(); i++) {
		remove_reset_coll_obj(branch_cobjs[i]);
	}
	for (unsigned i = 0; i < leaf_cobjs.size(); i++) {
		remove_reset_coll_obj(leaf_cobjs[i]);
	}
	branch_cobjs.clear();
	leaf_cobjs.clear();
}


void tree_cont_t::remove_cobjs() {

	for (iterator i = begin(); i != end(); ++i) {
		i->remove_collision_objects();
	}
}


void tree_cont_t::draw_branches_and_leaves(shader_t const &s, bool draw_branches, bool draw_leaves, bool shadow_only, vector3d const &xlate) {

	BLACK.do_glColor();

	for (iterator i = begin(); i != end(); ++i) {
		i->draw_tree(s, draw_branches, draw_leaves, shadow_only, xlate);
	}
}


void set_leaf_shader(shader_t &s, float min_alpha, bool gen_tex_coords, bool use_geom_shader, unsigned tc_start_ix) {

	s.set_prefix("#define USE_LIGHT_COLORS", 0); // VS - actually ignored due to custom lighting, but useful to have for reference
	if (!has_dl_sources)                 {s.set_prefix("#define NO_LEAF_DLIGHTS",     0);} // VS optimization
	if (gen_tex_coords)                  {s.set_prefix("#define GEN_QUAD_TEX_COORDS", 0);} // VS
	if (world_mode == WMODE_INF_TERRAIN) {s.set_prefix("#define USE_QUADRATIC_FOG",   1);} // FS
	s.setup_enabled_lights(2);
	s.set_frag_shader("linear_fog.part+textured_with_fog");

	if (use_geom_shader) { // unused
		s.set_vert_shader("ads_lighting.part*+leaf_lighting.part+tree_leaves_as_pts");
		s.set_geom_shader("output_textured_quad.part+point_to_quad", GL_POINTS, GL_TRIANGLE_STRIP, 4);
		s.begin_shader();
		//setup_wind_for_shader(s, 1); // FIXME: add wind?
	}
	else {
		s.set_vert_shader("ads_lighting.part*+leaf_lighting.part+tc_by_vert_id.part+tree_leaves");
		s.begin_shader();
	}
	s.setup_fog_scale();
	s.add_uniform_float("min_alpha", min_alpha);
	set_multitex(0);
	s.add_uniform_int("tex0", 0);
	s.add_uniform_int("tc_start_ix", tc_start_ix);
	check_gl_error(301);
}


void tree_cont_t::check_leaf_shadow_change() {

	static point last_lpos(all_zeros);
	point const lpos(get_light_pos());
		
	if (!no_sun_lpos_update && lpos != last_lpos) {
		for (iterator i = begin(); i != end(); ++i) {
			i->calc_leaf_shadows();
		}
	}
	last_lpos = lpos;
}


void tree_cont_t::pre_leaf_draw(shader_t &shader) {

	if (draw_model == 0) { // solid fill
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.75);
	}
	set_lighted_sides(2);
	set_specular(0.1, 10.0);
	glEnable(GL_COLOR_MATERIAL);
	glDisable(GL_NORMALIZE);
	
	if (shader.is_setup()) {
		shader.enable();
	}
	else {
		set_leaf_shader(shader, 0.75, 1, 0, 3);
	}
}


void tree_cont_t::post_leaf_draw(shader_t &shader) {

	shader.disable();
	set_lighted_sides(1);
	glDisable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	set_specular(0.0, 1.0);
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_TEXTURE_2D);
}


void tree_cont_t::draw(bool shadow_only) {

	// draw branches, then leaves: much faster for distant trees, slightly slower for near trees
	// draw branches
	shader_t bs;
	bool const branch_smap(1 && !shadow_only); // looks better, but slower
	set_tree_branch_shader(bs, !shadow_only, !shadow_only, branch_smap);
	draw_branches_and_leaves(bs, 1, 0, shadow_only, zero_vector);
	bs.add_uniform_vector3d("world_space_offset", zero_vector); // reset
	bs.end_shader();
	disable_multitex_a();

	// draw leaves
	shader_t ls;
	pre_leaf_draw(ls);
	draw_branches_and_leaves(ls, 0, 1, shadow_only, zero_vector);
	post_leaf_draw(ls);
}


void draw_trees(bool shadow_only) {

	//glFinish(); // testing
	//RESET_TIME;

	if (tree_mode & 2) { // small trees
		draw_small_trees(shadow_only);
		//PRINT_TIME("Small Trees");
	}
	if (tree_mode & 1) { // trees
		if (!shadow_only) {t_trees.check_leaf_shadow_change();}
		set_fill_mode();
		t_trees.draw(shadow_only);
		//glFinish(); // testing
		//PRINT_TIME(((tree_mode & 2) ? "Large + Small Trees" : "Large Trees"));
	}
	leaf_color_changed = 0; // Note: only visible trees will be updated
}


void tree_data_t::make_private_copy(tree_data_t &dest) const {
	
	dest = *this;
	dest.clear_vbo_ixs();
	//dest.ref_count = 1;
}


void tree::make_private_tdata_copy() {

	if (!tree_data || td_is_private()) return; // tree pointer is NULL or already private
	tree_data->make_private_copy(priv_tree_data);
	//tree_data->dec_ref_count();
	tree_data = NULL;
}


void tree::bind_to_td(tree_data_t *td) {

	assert(td);
	assert(!created); // too strong?
	assert(leaf_cobjs.empty() && branch_cobjs.empty());
	if (td_is_private()) {tdata().clear_data();}
	tree_data = td;
}


coll_obj &tree::get_leaf_cobj(unsigned i) const {

	assert(i < leaf_cobjs.size());
	assert(leaf_cobjs[i] < (int)coll_objects.size());
	return coll_objects[leaf_cobjs[i]];
}


void tree_data_t::gen_leaf_color() {

	leaf_color = get_leaf_base_color(tree_type) * leaf_color_coherence;
}


inline colorRGB tree_leaf::calc_leaf_color(colorRGBA const &leaf_color, colorRGBA const &base_color) const {

	float const ilch(1.0 - leaf_color_coherence);
	return colorRGB(color*(leaf_color.R + ilch*lred  ) + base_color.R*tree_color_coherence,
		            color*(leaf_color.G + ilch*lgreen) + base_color.G*tree_color_coherence, 0.0);
}


colorRGB tree_data_t::get_leaf_color(unsigned i) const {

	if (leaf_data_allocated()) {
		assert((i<<2) < leaf_data.size());
		return leaf_data[i<<2].get_c3(); // return color of first vertex since they all should be the same
	}
	assert(i < leaves.size());
	return leaves[i].calc_leaf_color(leaf_color, base_color);
}


void tree_data_t::update_leaf_color(unsigned i, bool no_mark_changed) {

	assert(i < leaves.size());
	colorRGB const color(leaves[i].calc_leaf_color(leaf_color, base_color));
	UNROLL_4X(leaf_data[i_+(i<<2)].set_c3(color);)
	if (!no_mark_changed) {leaves_changed = 1;}
}


void tree_data_t::update_all_leaf_colors() {

	for (unsigned i = 0; i < leaves.size(); i++) {update_leaf_color(i);}
}


void tree::update_leaf_cobj_color(unsigned i) {

	// update cobj color so that leaf water reflection is correct (Note: not always correct with tree instancing)
	get_leaf_cobj(i).cp.color = colorRGBA(tdata().get_leaf_color(i)).modulate_with(texture_color(tree_types[type].leaf_tex));
}


void tree::copy_color(unsigned i, bool no_mark_changed) {

	if (!has_leaf_data()) return;
	tdata().update_leaf_color(i, no_mark_changed);
	if (!leaf_cobjs.empty()) {update_leaf_cobj_color(i);}
}


void tree_data_t::remove_leaf_ix(unsigned i, bool update_data) {

	assert(i < leaves.size());
	leaves[i] = leaves.back();
	leaves.pop_back();
	if (!update_data) return;
	unsigned const i4(i << 2), tnl4((unsigned)leaves.size() << 2);
	assert(4*leaves.size() <= leaf_data.size());
	// shift vertex array (last one is now invalid)
	UNROLL_4X(leaf_data[i_+i4] = leaf_data[i_+tnl4];)
	leaves_changed = 1;
}


void tree::remove_leaf(unsigned i, bool update_data) {
	
	if (!leaf_cobjs.empty()) {
		assert(i < leaf_cobjs.size());
		remove_coll_object(leaf_cobjs[i]);
		leaf_cobjs[i] = leaf_cobjs.back();
		leaf_cobjs.pop_back();
	}
	make_private_tdata_copy();
	update_data &= has_leaf_data();
	tdata().remove_leaf_ix(i, update_data);
}


void tree::burn_leaves() {

	if (!physics_enabled()) return;
	tree_data_t &td(tdata());
	vector<tree_leaf> &leaves(td.get_leaves());
	float const max_t(get_max_t(LEAF));
	if (temperature < max_t || leaves.empty()) return;
	float const burn_amt(0.25*((td_is_private() || t_trees.empty()) ? 1.0 : float(tree_data_manager.size())/float(t_trees.size())));
	unsigned const num_burn(max(1U, min(5U, unsigned(5*(temperature - max_t)/max_t))));
	damage += ((1.0 - damage)*num_burn)/leaves.size();

	for (unsigned i = 0; i < num_burn && !leaves.empty(); ++i) {
		unsigned const index(rand()%leaves.size());
		leaves[index].color = max(0.0f, (leaves[index].color - burn_amt));
		copy_color(index);
		if (rand()&1) {gen_smoke(leaves[index].pts[0] + tree_center);}
		if (td_is_private() && leaves[index].color <= 0.0) {remove_leaf(index, 1);} // Note: if we modify shared data, leaves.size() must be dynamic
	}
}


bool tree::physics_enabled() const {

	return (created && (tree_mode & 1) && world_mode == WMODE_GROUND);
}


void tree::get_abs_leaf_pts(point pts[4], unsigned ix) const {
	
	vector<tree_leaf> const &leaves(tdata().get_leaves());
	UNROLL_4X(pts[i_] = leaves[ix].pts[i_]+tree_center;)
}


void tree::create_leaf_obj(unsigned ix) const {

	point pts[4];
	get_abs_leaf_pts(pts, ix);
	gen_leaf_at(pts, tdata().get_leaves()[ix].norm, type, tdata().get_leaf_color(ix));
}


bool tree::damage_leaf(unsigned i, float damage_done) {

	if (damage_done == 0.0) return 0;
	make_private_tdata_copy();
	vector<tree_leaf> &leaves(tdata().get_leaves());
	assert(i < leaves.size());
	float &lcolor(leaves[i].color);
	colorRGB const black(0.0, 0.0, 0.0);

	if (damage_done > 1.0 || (damage_done > 0.2 && lcolor == 0.0)) {
		if (lcolor > 0.0) damage += damage_scale;
		lcolor = -1.0;
		if ((rand()&3) == 0) {create_leaf_obj(i);} // 50/50 chance of the burned leaf falling
		
		if (has_leaf_data()) { // remove leaf i
			remove_leaf(i, 1);
			return 1;
		}
	}
	else {
		float const llc(lcolor);
		lcolor = max(0.0, (lcolor - 0.3*damage_done));
		if (lcolor == 0.0 && llc > 0.0) damage += damage_scale;
		copy_color(i);
	}
	return 0;
}


void tree::blast_damage(blastr const *const blast_radius) {

	if (!physics_enabled()) return;
	assert(blast_radius);
	float const bradius(blast_radius->cur_size), bdamage(LEAF_DAM_SCALE*blast_radius->damage);
	if (bdamage == 0.0) return;
	point const &bpos(blast_radius->pos);
	float const radius(bradius + tdata().sphere_radius), bradius_sq(bradius*bradius);
	if (p2p_dist_sq(bpos, sphere_center()) > radius*radius) return;
	unsigned nleaves(tdata().get_leaves().size());

	for (unsigned i = 0; i < nleaves; i++) {
		tree_leaf const &l(tdata().get_leaves()[i]);
		if (l.color < 0.0) continue;
		float const dist(p2p_dist_sq(bpos, l.pts[0]+tree_center));
		if (dist < TOLERANCE || dist > bradius_sq) continue;
		float const blast_damage(bdamage*InvSqrt(dist));
		if (damage_leaf(i, blast_damage)) {--i; --nleaves;} // force reprocess of this leaf, wraparound to -1 is OK
	} // for i
	damage = min(1.0f, damage);
}


void tree::lightning_damage(point const &ltpos) {

	blastr const br(0, ETYPE_FIRE, -2, BURN_RADIUS, BURN_DAMAGE, ltpos, plus_z, LITN_C, LITN_C);
	blast_damage(&br);
}


void tree::drop_leaves() {

	if (!physics_enabled()) return;
	tree_data_t &td(tdata());
	vector<tree_leaf> &leaves(td.get_leaves());
	unsigned const nleaves(leaves.size());
	if (damage >= 1.0 || nleaves == 0) return; // too damaged
	float const temp0(max(1.0f, min(0.3f, (20.0f-temperature)/30.0f)));
	int const rgen(min(LEAF_GEN_RAND2/10, int(rand_uniform(0.5, 1.5)*temp0*LEAF_GEN_RAND2/fticks)));

	for (unsigned i = (rand()%LEAF_GEN_RAND1); i < nleaves; i += LEAF_GEN_RAND1) {
		if ((rand()%rgen) != 0) continue;
		create_leaf_obj(i);
		
		if (td_is_private()) { // create a new leaf with a different color (and orient?)
			leaves[i].create_init_color(0);
			copy_color(i, 1); // don't set leaves_changed (too slow)
		}
	}
}


inline float tree_leaf::get_norm_scale(unsigned pt_ix) const {

	return ((shadow_bits & (1 << pt_ix)) ? LEAF_SHADOW_VAL : 1.0);
}


void tree_data_t::clear_vbo_ixs() {

	leaf_vbo = num_branch_quads = num_unique_pts = 0;
	branch_vbo_manager.reset_vbos_to_zero();
}


void tree_data_t::clear_vbos() {

	branch_vbo_manager.clear_vbos();
	delete_vbo(leaf_vbo);
	clear_vbo_ixs();
}


unsigned tree_data_t::get_gpu_mem() const {

	return (branch_vbo_manager.gpu_mem + (leaf_vbo ? leaf_data.size()*sizeof(leaf_vert_type_t) : 0));
}


bool tree::is_visible_to_camera(vector3d const &xlate) const {

	int const level((island || get_camera_pos().z > ztop) ? 1 : 2); // do we want to test the mesh in non-island mode?
	return sphere_in_camera_view((sphere_center() + xlate), 1.1*tdata().sphere_radius, level);
}


void tree_data_t::draw_tree_shadow_only(bool draw_branches, bool draw_leaves, int tree_type) const {

	if (draw_branches && branch_vbo_manager.vbo) { // draw branches (untextured)
		branch_vbo_manager.pre_render();
		set_array_client_state(1, 0, 0, 0); // vertices only
		glVertexPointer(3, GL_FLOAT, sizeof(branch_vert_type_t), 0);
		glDrawRangeElements(GL_QUADS, 0, num_unique_pts, num_branch_quads, GL_UNSIGNED_SHORT, 0);
		branch_vbo_manager.post_render();
	}
	if (draw_leaves && leaf_vbo > 0 && !leaves.empty()) { // draw leaves
		select_texture(tree_types[tree_type].leaf_tex);
		bind_vbo(leaf_vbo);
		assert(leaf_data.size() >= 4*leaves.size());
		leaf_vert_type_t::set_vbo_arrays(); // could also disable normals and colors, but that doesn't seem to help much
		glDrawArrays(GL_QUADS, 0, 4*(unsigned)leaves.size());
	}
}


void tree::draw_tree(shader_t const &s, bool draw_branches, bool draw_leaves, bool shadow_only, vector3d const &xlate) {

	if (!created) return;
	tree_data_t &td(tdata());

	if (shadow_only) {
		if (!is_over_mesh()) return;
		glPushMatrix();
		translate_to(tree_center + xlate);
		td.draw_tree_shadow_only(draw_branches, draw_leaves, type);
		glPopMatrix();
		return;
	}
	float dc_scale(1.0 - 0.95*damage);
	bcolor = tree_types[type].barkc;
	UNROLL_3X(bcolor[i_] *= min(1.0f, dc_scale*tree_color[i_]);)
	td.gen_leaf_color();
	bool const has_leaves(!td.get_leaves().empty());
	
	if (draw_leaves && has_leaves) {
		burn_leaves();
		if (l_strike.enabled == 1)    lightning_damage(l_strike.end);
		if (game_mode && gm_blast)    blast_damage(&latest_blastr);
		if (begin_motion && animate2) drop_leaves();
	}
	if (TEST_RTREE_COBJS) return;
	if (!draw_leaves || draw_branches) {not_visible = !is_visible_to_camera(xlate);} // second pass only
	if (not_visible) return;
	bool const use_vbos(setup_gen_buffers());
	assert(use_vbos);
	float const size_scale((do_zoom ? ZOOM_FACTOR : 1.0)*td.base_radius/(distance_to_camera((sphere_center() + xlate))*DIST_C_SCALE));
	if (draw_branches) draw_tree_branches(s, size_scale, xlate);
	if (draw_leaves && has_leaves) draw_tree_leaves(s, size_scale, xlate);
}


void tree_data_t::create_indexed_quads_for_branches(vector<branch_vert_type_t> &data, vector<branch_index_t> &idata) const {

	unsigned const numcylin((unsigned)all_cylins.size());
	unsigned num_quads(0), num_pts(0), cylin_id(0), data_pos(0), quad_id(0);

	for (unsigned i = 0; i < numcylin; ++i) { // determine required data size
		bool const prev_connect(i > 0 && all_cylins[i].can_merge(all_cylins[i-1]));
		unsigned const ndiv(all_cylins[i].get_num_div());
		num_quads += ndiv;
		num_pts   += (prev_connect ? 1 : 2)*ndiv;
	}
	assert(data.size() + num_pts < (1 << 8*sizeof(branch_index_t))); // cutting it close with 4th order branches
	data.reserve(num_pts);
	idata.reserve(4*num_quads);

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
				vector3d const n(0.5*vpn.n[S] + 0.5*vpn.n[(S+ndiv-1)%ndiv]); // average face normals to get vert normals
				data.push_back(branch_vert_type_t(vpn.p[(S<<1)+j], n, tx, float(cylin_id + j)));
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
}


void tree_data_t::draw_branches(float size_scale) {

	if (branch_vbo_manager.vbo == 0) { // vbos not yet created
		assert(branch_vbo_manager.ivbo == 0);
		assert(num_branch_quads == 0 && num_unique_pts == 0);
		vector<branch_vert_type_t> data;
		vector<branch_index_t> idata;
		create_indexed_quads_for_branches(data, idata);
		assert(data.size()  == data.capacity());
		assert(idata.size() == idata.capacity());
		num_branch_quads = idata.size()/4;
		num_unique_pts   = data.size();
		branch_vbo_manager.create_and_upload(data, idata); // ~350KB data + ~75KB idata (with 16-bit index)
	}
	branch_vbo_manager.pre_render();
	vert_norm_comp_tc::set_vbo_arrays();
	unsigned const num(4*min(num_branch_quads, max((num_branch_quads/8), unsigned(1.5*num_branch_quads*size_scale)))); // branch LOD
	glDrawRangeElements(GL_QUADS, 0, num_unique_pts, num, GL_UNSIGNED_SHORT, 0); // draw with branch vbos
	branch_vbo_manager.post_render();
}


void tree::draw_tree_branches(shader_t const &s, float size_scale, vector3d const &xlate) {

	if (size_scale < 0.05) return; // too far away, don't draw any branches
	select_texture(tree_types[type].bark_tex);
	set_color(bcolor);
	if (s.is_setup()) {s.add_uniform_vector3d("world_space_offset", (tree_center + xlate));}
	glPushMatrix();
	translate_to(tree_center + xlate);
	tdata().draw_branches(size_scale);
	glPopMatrix();
}


void tree_data_t::update_normal_for_leaf(unsigned i) {

	assert(i < leaves.size());
	UNROLL_4X(leaf_data[i_+(i<<2)].set_norm(leaves[i].norm*leaves[i].get_norm_scale(i_));)
	leaves_changed = 1;
}


void tree_data_t::reset_leaf_pos_norm() {

	for (unsigned i = 0; i < (unsigned)leaves.size(); i++) { // process leaf points - reset to default positions and normals
		UNROLL_4X(leaf_data[i_+(i<<2)].v = leaves[i].pts[i_];)
		update_normal_for_leaf(i);
	}
	reset_leaves = 0;
}


void tree_data_t::draw_leaves(float size_scale) {

	bool const create_leaf_vbo(leaf_vbo == 0);
	unsigned nl(leaves.size());
	assert(leaf_data.size() >= 4*nl);
	create_bind_vbo_and_upload(leaf_vbo, leaf_data, 0);

	if (!create_leaf_vbo && leaves_changed) {
		upload_vbo_sub_data(&leaf_data.front(), 0, leaf_data.size()*sizeof(leaf_vert_type_t)); // FIXME: track and update a sub-range?
	}
	leaves_changed = 0;
	leaf_vert_type_t::set_vbo_arrays();
	if (ENABLE_CLIP_LEAVES) {nl = min(nl, max((nl/8), unsigned(4.0*nl*size_scale)));} // leaf LOD
	glDrawArrays(GL_QUADS, 0, 4*nl);
	bind_vbo(0);
}


bool tree_data_t::leaf_draw_setup(bool leaf_dynamic_en) {

	bool const gen_arrays(!leaf_data_allocated());
	if (gen_arrays) {alloc_leaf_data();}
	if (gen_arrays || (reset_leaves && !leaf_dynamic_en)) {reset_leaf_pos_norm();}
	// may do duplicate work, but should at least keep the cobj colors in sync
	if (gen_arrays || leaf_color_changed) {update_all_leaf_colors();}
	return gen_arrays;
}


void tree::draw_tree_leaves(shader_t const &s, float size_scale, vector3d const &xlate) {

	tree_data_t &td(tdata());
	bool const leaf_dynamic_en(!has_snow && (display_mode & 0x0100) != 0);
	bool const gen_arrays(td.leaf_draw_setup(leaf_dynamic_en));
	if (!gen_arrays && leaf_dynamic_en && size_scale > 0.5) {update_leaf_orients();}

	if (gen_arrays || leaf_color_changed) {
		for (unsigned i = 0; i < leaf_cobjs.size(); ++i) {update_leaf_cobj_color(i);}
	}
	if (gen_arrays) {calc_leaf_shadows();}
	unsigned num_dlights(0);
	bool const enable_dlights(has_dl_sources && world_mode == WMODE_GROUND && s.is_setup());

	if (enable_dlights) {
		num_dlights = enable_dynamic_lights((sphere_center() + xlate), td.sphere_radius);
		s.add_uniform_int("num_dlights", num_dlights);
	}
	select_texture((draw_model == 0) ? tree_types[type].leaf_tex : WHITE_TEX); // what about texture color mod?
	glPushMatrix();
	translate_to(tree_center + xlate);
	td.draw_leaves(size_scale);
	if (enable_dlights) {disable_dynamic_lights(num_dlights);}
	glPopMatrix();
}


void tree_data_t::bend_leaf(unsigned i, float angle) {

	assert(i < leaves.size());
	tree_leaf const &l(leaves[i]);
	point const p1((l.pts[1] + l.pts[2])*0.5); // tip
	point const p2((l.pts[0] + l.pts[3])*0.5); // base
	vector3d const orig_dir(p1 - p2); // vector from base to tip
	vector3d const new_dir(orig_dir*cos(angle) + l.norm*(orig_dir.mag()*sin(angle))); // s=orig_dir.get_norm(), t=l.norm
	vector3d const delta(new_dir - orig_dir);
	vector3d normal(cross_product(new_dir, (l.pts[3] - l.pts[0])).get_norm());
	unsigned const ix(i<<2);
	leaf_data[ix+1].v = l.pts[1] + delta;
	leaf_data[ix+2].v = l.pts[2] + delta;

	if (l.shadow_bits == 0) {
		norm_comp nc(normal);
		UNROLL_4X((norm_comp)leaf_data[i_+ix] = nc;)
	}
	else { // update the normals, even though this slows the algorithm down
		UNROLL_4X(leaf_data[i_+ix].set_norm(normal*l.get_norm_scale(i_));) // almost update_normal_for_leaf(i);
	}
	leaves_changed = 1;
	reset_leaves   = 1; // do we want to update the normals as well?
}


bool tree_data_t::check_if_needs_updated() {

	bool const do_update(last_update_frame < frame_counter);
	last_update_frame = frame_counter;
	return do_update;
}


void tree::update_leaf_orients() { // leaves move in wind or when struck by an object (somewhat slow)

	int last_xpos(0), last_ypos(0);
	vector3d local_wind;
	tree_data_t &td(tdata());
	bool const do_update(td.check_if_needs_updated()), priv_data(td_is_private());
	if (!do_update && !physics_enabled()) return;
	unsigned nleaves(td.get_leaves().size());

	for (unsigned i = 0; i < nleaves; i++) { // process leaf wind and collisions
		float wscale(0.0), hit_angle(0.0);

		if (do_update) {
			tree_leaf const &leaf(td.get_leaves()[i]);
			point p0(leaf.pts[0]);
			if (priv_data) {p0 += tree_center;}
			int const xpos(get_xpos(p0.x)), ypos(get_ypos(p0.y));
			
			// Note: should check for similar z-value, but z is usually similar within the leaves of a single tree
			if (i == 0 || xpos != last_xpos || ypos != last_ypos) {
				local_wind = get_local_wind(p0, !priv_data); // slow
				last_xpos  = xpos;
				last_ypos  = ypos;
			}
			wscale = dot_product(local_wind, leaf.norm);
		}
		if (!leaf_cobjs.empty()) { // rotate leaves when hit by an object
			coll_obj &cobj(get_leaf_cobj(i));
				
			if (cobj.last_coll > 0) {
				bool removed(0);

				// Note: the code below can cause the tree_data_t to be deep copied, invalidating any iterators
				if (cobj.coll_type == BEAM) { // do burn damage
					removed = damage_leaf(i, BEAM_DAMAGE);
				}
				else {
					if (cobj.coll_type == PROJECTILE) {
						if ((rand()&3) == 0) { // shoot off leaf
							create_leaf_obj(i);
							remove_leaf(i, 1);
							removed = 1;
						}
						cobj.coll_type = IMPACT; // reset to impact after first hit
					}
					hit_angle += PI_TWO*cobj.last_coll/TICKS_PER_SECOND; // 90 degree max rotate
				}
				if (removed) {
					--i; --nleaves;
					continue;
				}
				cobj.last_coll = ((cobj.last_coll < iticks) ? 0 : (cobj.last_coll - iticks));
			}
		}
		if (wscale != 0.0 || hit_angle != 0.0) {
			float const angle(0.5*PI*max(-1.0f, min(1.0f, wscale)) + hit_angle); // not physically correct, but it looks good
			td.bend_leaf(i, angle); // do we want to update collision objects as well?
		}
		float &lcolor(td.get_leaves()[i].color);

		if (LEAF_HEAL_RATE > 0.0 && lcolor > 0.0 && lcolor < 1.0) { // leaf heal
			if (do_update) {lcolor = min(1.0f, (lcolor + LEAF_HEAL_RATE*fticks));}
			copy_color(i);
		}
	} // for i
}


void tree::calc_leaf_shadows() { // process leaf shadows/normals

	if (!physics_enabled()) return;
	tree_data_t &td(tdata());
	if (!td.leaf_data_allocated() || !td_is_private()) return; // leaf data not yet created (can happen if called when light source changes), or shared data
	int const light(get_light());
	vector<tree_leaf> &leaves(td.get_leaves());

	for (unsigned i = 0; i < leaves.size(); i++) {
		tree_leaf &l(leaves[i]);
		l.shadow_bits = 0;

		for (unsigned j = 0; j < 4; ++j) {
			bool const shadowed(!leaf_cobjs.empty() && !is_visible_to_light_cobj(l.pts[j]+tree_center, light, 0.0, leaf_cobjs[i], 1));
			l.shadow_bits |= (int(shadowed) << j);
		}
		td.update_normal_for_leaf(i);
	}
}


unsigned tree_cont_t::delete_all() {

	unsigned deleted(0);

	for (reverse_iterator i = rbegin(); i != rend(); ++i) { // delete backwards (pop collision stack)
		if (i->delete_tree()) ++deleted;
	}
	generated = 0;
	return deleted;
}


void delete_trees() {

	t_trees.delete_all();
	if (tree_coll_level && !t_trees.empty()) {purge_coll_freed(1);} // MUST do this for collision detection to work
}


void tree::clear_vbo() {
	
	if (td_is_private()) {tdata().clear_vbos();}
}


int tree::delete_tree() {

	clear_vbo();
	if (!created)  return 0;
	if (tree_coll_level) {remove_collision_objects();}
	if (no_delete) return 0;
	if (td_is_private()) {tdata().clear_data();}
	created = 0;
	return 1;
}


void tree_data_t::clear_data() {
	
	clear_vbos();
	all_cylins.clear(); remove_excess_cap(all_cylins);
	leaf_data.clear();  remove_excess_cap(leaf_data);
	leaves.clear();     remove_excess_cap(leaves); // Note: not present in original delete_trees()
}


void tree_builder_t::process_cylins(tree_cylin const *const cylins, unsigned num, int tree_type, float deadness,
	vector<draw_cylin> &all_cylins, vector<tree_leaf> &leaves) const
{
	float const leaf_size(tree_types[tree_type].leaf_size*(tree_scale*base_radius/TREE_SIZE + 10.0)/18.0);

	for (unsigned i = 0; i < num; ++i) {
		assert(cylins[i].r1 > 0.0 || cylins[i].r2 > 0.0);
		assert(cylins[i].p1 != cylins[i].p2);
		all_cylins.push_back(cylins[i]);

		if (deadness < 1.0) { // leaves was reserved
			if (cylins[i].level > 1 && (cylins[i].level < 4 || TREE_4TH_BRANCHES > 1)) { // leaves will still be allocated
				add_leaves_to_cylin(cylins[i], leaf_size, deadness, leaves);
			}
		}
	}
}


void tree_builder_t::create_all_cylins_and_leaves(int tree_type, float deadness, vector<draw_cylin> &all_cylins, vector<tree_leaf> &leaves) {

	//process cylinders
	assert(all_cylins.empty());
	unsigned num_total_cylins(base_num_cylins + root_num_cylins); // start with trunk and root cylinders

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

	// tree leaf variables
	if (deadness < 1.0) {
		unsigned nl(unsigned((1.0 - deadness)*num_leaves_per_occ*num_total_cylins) + 1); // determine the number of leaves
		leaves.reserve(nl);
	}
	process_cylins(base.cylin,  base_num_cylins, tree_type, deadness, all_cylins, leaves);
	process_cylins(roots.cylin, root_num_cylins, tree_type, deadness, all_cylins, leaves);

	for (int i = 0; i < num_1_branches; i++) { // add the first order branches
		process_cylins(branches[i][0].cylin, branches[i][0].num_cylins, tree_type, deadness, all_cylins, leaves);
	}
	for (int i = 0; i < num_1_branches; i++) { // add second order branches
		for (int k = 1; k <= branches[i][0].num_branches; k++) {
			process_cylins(branches[i][k].cylin, branches[i][k].num_cylins, tree_type, deadness, all_cylins, leaves);
		}
	}
	for (unsigned w = 0; w < 2; ++w) {
		for (int i = 0; i < num_34_branches[w]; i++) { // add the third and fourth order branches
			process_cylins(branches_34[w][i].cylin, branches_34[w][i].num_cylins, tree_type, deadness, all_cylins, leaves);
		}
	}

	// now NULL (pseudo free) the individual branches because they are unneccessary
	branches_34[0] = branches_34[1] = NULL;
	branches = NULL;
	roots.cylin = base.cylin = NULL;
	assert(all_cylins.size() == all_cylins.capacity());
}


float tree_builder_t::get_bsphere_center_zval() const {

	assert(base_num_cylins > 0);
	assert(base.cylin);
	return base.cylin[max(0, (base_num_cylins - 2))].p2.z;
}


inline void add_rotation(point &dest, point const &src, float mult) {

	UNROLL_3X(dest[i_] = src[i_] + mult*re_matrix[i_];)
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


float get_default_tree_depth() {
	return TREE_DEPTH*(0.5 + 0.5/tree_scale);
}


//gen_tree(pos, size, ttype>=0, calc_z, 0, 1);
//gen_tree(pos, 0, -1, 1, 1, 0);
void tree::gen_tree(point const &pos, int size, int ttype, int calc_z, bool add_cobjs, bool user_placed) {

	assert(calc_z || user_placed);
	tree_center = pos;
	created     = 1;
	tree_color.alpha = 1.0;
	UNROLL_3X(tree_color[i_] = 0.01*rand_gen(80, 120);) // rand_gen() called outside gen_tree_data()
	tree_data_t &td(tdata());

	if (td.is_created()) { // pre-allocated, shared tree
		assert(!td_is_private());
		assert(!user_placed);
		assert(size <= 0); // too strong?
		int const td_type(tdata().get_tree_type());
		
		if (ttype < 0) {
			type = (FORCE_TREE_TYPE ? td_type : (rand2() % NUM_TREE_TYPES));
		}
		else if (FORCE_TREE_TYPE) {
			type = ttype;
			// Note: we could relax this restruction and let type be different, which gives us more tree variation,
			// but then the leaf and branch sizes/colors couldn't be different per tree
			assert(type == td_type); // up to the caller to ensure this
		}
	}
	else {
		type = ((ttype < 0) ? rand2() : ttype) % NUM_TREE_TYPES; // maybe should be an error if > NUM_TREE_TYPES
		float tree_depth(get_default_tree_depth());
	
		if (calc_z) {
			tree_center.z = interpolate_mesh_zval(tree_center.x, tree_center.y, 0.0, 1, 1);
			if (user_placed) {tree_depth = tree_center.z - get_tree_z_bottom(tree_center.z, tree_center);} // more accurate
		}
		td.gen_tree_data(type, size, tree_depth); // create the tree here
	}
	assert(type < NUM_TREE_TYPES);
	unsigned const nleaves(td.get_leaves().size());
	damage_scale = (nleaves ? 1.0/nleaves : 0.0);
	damage       = 0.0;
	if (add_cobjs) {add_tree_collision_objects();}
}


void tree_data_t::gen_tree_data(int tree_type_, int size, float tree_depth) {

	tree_type = tree_type_;
	assert(tree_type < NUM_TREE_TYPES);
	calc_leaf_points(); // required for placed trees
	leaf_data.clear();
	clear_vbo_ixs();
	float deadness(DISABLE_LEAVES ? 1.0 : tree_deadness);

	if (deadness < 0.0) {
		int const num(rand_gen(1, 100));
		deadness = ((num > 94) ? min(1.0f, float(num - 94)/8.0f) : 0.0);
	}
	tree_builder_t builder;
	base_radius = builder.create_tree_branches(tree_type, size, tree_depth, base_color);
	
	// create leaves and all_cylins
	sphere_center_zoff = builder.get_bsphere_center_zval();
	builder.create_all_cylins_and_leaves(tree_type, deadness, all_cylins, leaves);

	// set the bounding sphere center
	sphere_radius = 0.0;

	for (vector<draw_cylin>::const_iterator i = all_cylins.begin(); i != all_cylins.end(); ++i) {
		sphere_radius = max(sphere_radius, p2p_dist_sq(i->p2, vector3d(0.0, 0.0, sphere_center_zoff)));
	}
	sphere_radius = sqrt(sphere_radius);

	/*for (unsigned i = 0; i < leaves.size(); ++i) { // scramble leaves so that LOD is unbiased/randomly sampled
		swap(leaves[i], leaves[(i + 1572869)%leaves.size()]);
	}*/
	reverse(leaves.begin(), leaves.end()); // order leaves so that LOD removes from the center first, which is less noticeable
}


float tree_builder_t::create_tree_branches(int tree_type, int size, float tree_depth, colorRGBA &base_color) {

	//fixed tree variables
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
	num_34_branches[1]   = (TREE_4TH_BRANCHES ? 2000 : 0);
	if (size <= 0) size  = rand_gen(40, 80); // tree size
	base_radius            = size * (0.1*TREE_SIZE*tree_types[tree_type].branch_size/tree_scale);
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
	float const branch_1_random_rotate = 40.0;
	int const min_num_roots(10), max_num_roots(12);
	base_color = colorRGBA(0.5*signed_rand_float2(), 0.5*signed_rand_float2(), 0.0, 1.0); // no blue

	//temporary variables
	int num_b_so_far(0);
	int const nbr(num_1_branches*(num_2_branches_max+1)), nbranches(nbr*num_3_branches_max);
	
	//allocate all the memory required for the tree------------------------------------------------
	branch_ptr_cache.resize(num_1_branches);
	branches = &branch_ptr_cache.front();
	unsigned const tot_branches(nbr + nbranches + max(num_34_branches[1], 1)); // allocate at least one branches_34[1] so we can set the first cylin
	branch_cache.resize(tot_branches);
	branches[0]    = &branch_cache.front();
	branches_34[0] = branches   [0] + nbr;
	branches_34[1] = branches_34[0] + nbranches;
	unsigned const tot_cylins(base_num_cylins+1 + CYLINS_PER_ROOT*max_num_roots + ncib*(nbr + nbranches + num_34_branches[1]));
	cylin_cache.resize(tot_cylins);
	base.cylin  = &cylin_cache.front();
	roots.cylin = base.cylin + base_num_cylins+1;

	for (int i = 0; i < num_1_branches; i++) {
		branches[i] = branches[0] + i*(num_2_branches_max+1);

		for (int j = 0; j < (num_2_branches_max+1); j++) {
			branches[i][j].cylin = roots.cylin + CYLINS_PER_ROOT*max_num_roots + (i*(num_2_branches_max+1) + j)*ncib;
		}
	}
	for (int i = 0; i < nbranches; i++) {
		branches_34[0][i].cylin = branches[0][0].cylin + (nbr + i)*ncib;
	}

	//create tree base -------------------------------------------------------------
	for (int i = 0; i < base_num_cylins; i++) {
		tree_cylin &cylin(base.cylin[i]);

		if (i == 0) {
			float const length((base_length_max - base_length_min)/2.0);
			cylin.assign_params(0, 0, base_radius, base_radius*base_var, length*(num_cylin_factor/base_num_cylins), 0.0);
			cylin.p1 = all_zeros;
			cylin.rotate.assign(cosf(angle_rotate/TO_DEG), sinf(angle_rotate/TO_DEG), 0.0);
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
			int const init_num_cyl(base_num_cylins - base_break_off + 1);
			float const temp_num(branch_distribution*((float)num_1_branches)/init_num_cyl*(float(i+2-base_break_off))/(float(num_b_so_far+1)));
			float rotate_start(rand_gen(0, 259));
			int num_b_to_create(0);

			if (temp_num >= 1.0) {
				float const temp2(((float)num_1_branches - num_b_so_far)/(float(base_num_cylins - i)));
				num_b_to_create = min(num_1_branches, ((temp2 <= 3.0) ? max(1, int(ceil(temp2))) : int(temp_num + 0.5)));
			}
			for (int j = 0; j < num_b_to_create; j++) {
				if (num_b_so_far < num_1_branches) create_1_order_branch(i, rotate_start, num_b_so_far++);
				rotate_start += 360.0/((float)num_b_to_create) + (3 - 2*rand_gen(1,2))*rand_gen(0,int(branch_1_random_rotate));
			}
		}
	}
	
	//done with base ------------------------------------------------------------------
	if (TREE_4TH_BRANCHES) {create_4th_order_branches(nbranches);}

	// FIXME: hack to prevent roots from being generated in scenes where trees are user-placed and affect the lighting voxel sparsity
	root_num_cylins = (gen_tree_roots ? CYLINS_PER_ROOT*rand_gen(min_num_roots, max_num_roots) : 0);

	for (int i = 0; i < root_num_cylins; i += CYLINS_PER_ROOT) { // add roots
		tree_cylin &cylin1(roots.cylin[i]), &cylin2(roots.cylin[i+1]), &cylin3(roots.cylin[i+2]);
		float const root_radius(rand_uniform2(0.38, 0.45)*base_radius);
		float const theta((TWO_PI*(i + 0.3*signed_rand_float2()))/root_num_cylins);
		float const deg_rot(180.0+rand_uniform2(40.0, 50.0));
		vector3d const dir(sin(theta), cos(theta), 0.0);
		cylin1.assign_params(1, i, root_radius, 0.75*root_radius, 1.0*base_radius, deg_rot); // level 1, with unique branch_id's
		cylin1.p1     = cylin1.p2 = point(0.0, 0.0, 0.75*base_radius);
		cylin1.rotate = dir;
		rotate_cylin(cylin1);
		cylin1.p1    += (0.3*base_radius/cylin1.length)*(cylin1.p2 - cylin1.p1); // move away from the tree centerline
		cylin1.p1    *= 2.0;
		cylin1.p2    *= 1.3;

		cylin2.assign_params(1, i, 0.75*root_radius, 0.5*root_radius, 2.0*base_radius, deg_rot+10.0);
		cylin2.p1     = cylin1.p2 = cylin1.p2;
		cylin2.rotate = cylin1.rotate;
		rotate_cylin(cylin2);

		cylin3.assign_params(1, i, 0.5*root_radius, 0.0, 4.0*base_radius, deg_rot-24.0);
		cylin3.p1     = cylin2.p2 = cylin2.p2;
		cylin3.rotate = cylin2.rotate;
		rotate_cylin(cylin3);

		if (tree_depth > 0.0) { // keep roots above the tree depth
			for (unsigned j = i; j < i + CYLINS_PER_ROOT; ++j) {
				roots.cylin[j].p1.z = max(roots.cylin[j].p1.z, -tree_depth+roots.cylin[j].r1);
				roots.cylin[j].p2.z = max(roots.cylin[j].p2.z, -tree_depth+roots.cylin[j].r2);
			}
		}
	}
	if (tree_depth > 0.0 && root_num_cylins == 0) { // add the bottom cylinder section from the base into the ground
		tree_cylin &cylin(base.cylin[base_num_cylins]);
		cylin.assign_params(0, 1, base_radius, base_radius, tree_depth, 180.0);
		cylin.p1 = cylin.p2 = all_zeros;
		cylin.rotate = plus_x;
		rotate_cylin(cylin);
		++base_num_cylins;
	}
	return base_radius;
}


inline float tree_builder_t::gen_bc_size(float branch_var) {

	return (rand_gen((int)branch_var - 5, (int)branch_var + 5)/100.0)*(num_cylin_factor/ncib);
}


inline float tree_builder_t::gen_bc_size2(float branch_var) {

	return rand_gen((int)branch_var - 5, (int)branch_var + 5)/100.0;
}


void tree_builder_t::gen_next_cylin(tree_cylin &cylin, tree_cylin &lcylin, float var, float rad_var, int level, int branch_id, bool rad_var_test) {

	cylin.assign_params(level, branch_id, lcylin.r2, (rad_var_test ? lcylin.r1*gen_bc_size(rad_var) : 0.0),
		((level == 4) ? branch_4_length*var : lcylin.length*gen_bc_size(var)), lcylin.deg_rotate);
	rotate_around_axis(lcylin);
	add_rotation(cylin.p1, lcylin.p1, BRANCH_LEN_SCALE);
}


void tree_builder_t::gen_first_cylin(tree_cylin &cylin, tree_cylin &src_cylin, float bstart, float rad_var, float rotate_start, int level, int branch_id) {

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


void tree_builder_t::create_1_order_branch(int base_cylin_num, float rotate_start, int branch_num) {

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
	float const temp_num(2.0*temp_num_big_branches/(0.5*branch.num_cylins+1));

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
			temp_num2 = (j+2)*temp_num_big_branches/((num_2_branches_created+1)*(0.5*branch.num_cylins+1));
		}
		else {
			temp_num2 = (j+1)*branch.num_branches/((num_2_branches_created+1)*branch.num_cylins) +
				              branch.num_branches/((num_2_branches_created+1)*branch.num_cylins);
		}
		if (temp_num2*branch_1_distribution >= 1.0 && num_2_branches_created < branch.num_branches) branch_just_created = true;
		int const deg_added(generate_next_cylin(j, branch.num_cylins, branch_just_created, branch_deflected));
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


void tree_builder_t::create_2nd_order_branch(int i, int j, int cylin_num, bool branch_deflected, int rotation) {
	
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
	float const temp_num(2.0*branch.num_branches/((float)branch.num_cylins));

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
		float const temp_num2(float((index+2)*branch.num_branches)/((num_3_branches_created+1)*((float)branch.num_cylins)));
		if (temp_num*branch_1_distribution >= 1.0 && num_3_branches_created < branch.num_branches) branch_just_created = true;
		int const deg_added(generate_next_cylin(index, branch.num_cylins, branch_just_created, branch_deflected));
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


void tree_builder_t::create_3rd_order_branch(int i, int j, int cylin_num, int branch_num, bool branch_deflected, int rotation) {
	
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
		int const deg_added(generate_next_cylin(index, branch.num_cylins, branch_just_created, branch_deflected));
		cylin.deg_rotate += ((branch.cylin[0].deg_rotate < 0.0) ? deg_added : -deg_added);
		add_rotation(lcylin.p2, lcylin.p1, 1.0);
		if (index == (branch.num_cylins - 1)) rotate_cylin(cylin);
	}
}


void tree_builder_t::gen_b4(tree_branch &branch, int &branch_num, int num_4_branches, int i, int k) {

	int ncib(branch.num_cylins);
	float const branch_4_distribution = 0.2;

	for (int j = 0; j < ncib-1; j++) {
		if ((float(((float)j)/((float)ncib))) >= branch_4_distribution) {
			tree_cylin &cylin(branch.cylin[j]);
			float temp_deg(safe_acosf(cylin.rotate.x)), rotate_start(0.0);
			if (cylin.rotate.y < 0.0) temp_deg *= -1.0;

			for (int l = 0; l < num_4_branches; l++) {
				rotate_start += 360.0/num_4_branches*l + rand_gen(1,30);
				if (rotate_start > 360.0) rotate_start -= 360.0;
				generate_4th_order_branch(branch, j, rotate_start, temp_deg, branch_num++);
			}
		}
	}
}


void tree_builder_t::create_4th_order_branches(int nbranches) {

	int num_4_branches  = 2;
	branch_4_length     = 0.006; //0.03;
	branch_4_max_radius = 0.008;
	assert(num_34_branches[1] > 0);
	int branch_num(0);

	for (int i = 0; i < num_34_branches[1]; i++) {
		branches_34[1][i].cylin = branches_34[0][0].cylin + (nbranches + i)*ncib;
	}
	for (int i = 0; i < num_1_branches; i++) { // b1->b4
		for (int k = 0; k < (branches[i][0].num_branches + 1); k++) {
			gen_b4(branches[i][k], branch_num, num_4_branches, i, k); // b2->b4
		}
	}
	for (int i = 0; i < num_34_branches[0]; i++) { // b2->b4
		gen_b4(branches_34[0][i], branch_num, num_4_branches, i, 0);
	}
	num_34_branches[1] = branch_num;
}


void tree_builder_t::generate_4th_order_branch(tree_branch &src_branch, int j, float rotate_start, float temp_deg, int branch_num) {
	
	float const branch_4_rad_var = 85.0;
	float const branch_4_var     = 0.70;
	int index(0);
	bool branch_deflected(false);
	tree_branch &branch(branches_34[1][branch_num]);
	tree_cylin &cylin(branch.cylin[0]);
	setup_rotate(cylin.rotate, rotate_start, temp_deg);
	tree_cylin &src_cylin(src_branch.cylin[j]);
	branch.num_cylins = 10; // number of cylinders
	float const radius1(src_branch.cylin[int(0.65*src_branch.num_cylins)].r2*0.9);
	cylin.assign_params(4, branch_num, radius1, radius1*gen_bc_size2(branch_4_rad_var), branch_4_length,
		src_cylin.deg_rotate + ((src_cylin.deg_rotate > 0.0) ? 1 : -1)*rand_gen(0,60));
	rotate_around_axis(src_cylin);
	add_rotation(cylin.p1, src_cylin.p1, BRANCH_LEN_SCALE);
	if (branch.num_cylins < 2) rotate_cylin(cylin);
	
	for (index = 1; index < branch.num_cylins; index++) {
		tree_cylin &cylin(branch.cylin[index]), &lcylin(branch.cylin[index-1]);
		gen_next_cylin(cylin, lcylin, branch_4_var, branch_4_rad_var, 4, branch_num, (index < branch.num_cylins-1));
		int const deg_added(generate_next_cylin(index, branch.num_cylins, false, branch_deflected));
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


void tree_builder_t::add_leaves_to_cylin(tree_cylin const &cylin, float tsize, float deadness, vector<tree_leaf> &leaves) const {

	static float acc(0.0);
	point start;
	vector3d rotate;
	acc += num_leaves_per_occ;
	int const temp((int)acc);
	acc -= (float)temp;
	float const temp_deg(((cylin.rotate.y < 0.0) ? -1.0 : 1.0)*safe_acosf(cylin.rotate.x));
	float rotate_start(0.0);

	for (int l = 0; l < temp; l++) {
		if (deadness > 0 && deadness > rand_float2()) continue;
		rotate_start += 360.0/temp*l + rand_gen(1,30);
		if (rotate_start > 360.0) rotate_start -= 360.0;
		setup_rotate(rotate, rotate_start, temp_deg);
		rotate_around_axis(cylin);
		add_rotation(start, cylin.p1, 0.9);
		tree_leaf leaf;
		int const val(rand_gen(0,60));
		float const deg_rotate(cylin.deg_rotate + ((cylin.deg_rotate > 0.0) ? val : -val));
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
}


int tree_builder_t::generate_next_cylin(int cylin_num, int ncib, bool branch_just_created, bool &branch_deflected) {

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


unsigned tree_cont_t::scroll_trees(int ext_x1, int ext_x2, int ext_y1, int ext_y2) {

	unsigned nkeep(0);
	vector3d const vd(-dx_scroll*DX_VAL, -dy_scroll*DY_VAL, 0.0);

	for (iterator i = begin(); i != end(); ++i) { // keep any tree that's still in the scene
		point const &gen_pos(i->get_center());
		int const xp(get_xpos(gen_pos.x) - dx_scroll), yp(get_ypos(gen_pos.y) - dy_scroll);
		bool const keep(xp >= ext_x1 && xp <= ext_x2 && yp >= ext_y1 && yp <= ext_y2);
		i->set_no_delete(keep);

		if (keep) { // shift it - don't have to recreate it
			i->shift_tree(vd);
			++nkeep;
		}
	}
	return nkeep;
}


void tree_cont_t::post_scroll_remove() {

	for (unsigned i = 0; i < size(); ++i) {
		if (at(i).get_no_delete()) {
			at(i).add_tree_collision_objects();
			at(i).set_no_delete(0);
		}
		else { // remove this tree from the vector
			std::swap(at(i), back()); // deep copy = bad?
			pop_back();
			--i;
		}
	}
}


void tree_cont_t::gen_deterministic(int x1, int y1, int x2, int y2, float vegetation_) {

	bool const NONUNIFORM_TREE_DEN = 1; // based on world_mode?
	float const min_tree_h(island ? TREE_MIN_H : (water_plane_z + 0.01*zmax_est));
	float const max_tree_h(island ? TREE_MAX_H : 1.8*zmax_est);
	unsigned const mod_num_trees(num_trees/(NONUNIFORM_TREE_DEN ? TREE_DEN_THRESH : 1.0));
	unsigned const smod(3.321*XY_MULT_SIZE+1), tree_prob(max(1U, XY_MULT_SIZE/mod_num_trees));
	unsigned const skip_val(max(1, int(1.0/sqrt(tree_scale)))); // similar to deterministic gen in scenery.cpp
	shared_tree_data.ensure_init();
	mesh_xy_grid_cache_t density_gen[NUM_TREE_TYPES+1];

	if (NONUNIFORM_TREE_DEN) {
		for (unsigned i = 0; i <= NUM_TREE_TYPES; ++i) {
			float const tds(TREE_DIST_SCALE*(XY_MULT_SIZE/16384.0)*(i==0 ? 1.0 : 0.1)), xscale(tds*DX_VAL*DX_VAL), yscale(tds*DY_VAL*DY_VAL);
			density_gen[i].build_arrays(xscale*(x1 + xoff2 + 1000*i), yscale*(y1 + yoff2 - 1500*i), xscale, yscale, (x2-x1), (y2-y1));
		}
	}
	for (int i = y1; i < y2; i += skip_val) {
		for (int j = x1; j < x2; j += skip_val) {
			if (scrolling) {
				int const ox(j + dx_scroll), oy(i + dy_scroll); // positions in original coordinate system
				if (ox >= x1 && ox <= x2 && oy >= y1 && oy <= y2) continue; // use orignal tree from last position
			}
			global_rand_gen.rseed1 = 805306457*(i + yoff2) + 12582917*(j + xoff2) + 100663319*rand_gen_index;
			global_rand_gen.rseed2 = 6291469  *(j + xoff2) + 3145739 *(i + yoff2) + 1572869  *rand_gen_index;
			rand2_mix();
			unsigned const val(((unsigned)rand2_seed_mix())%smod);
			if (val <= 100)         continue; // scenery
			if (val%tree_prob != 0) continue; // not selected
			if ((global_rand_gen.rseed1&127)/128.0 >= vegetation_) continue;
			point pos((get_xval(j) + 0.5*DX_VAL*rand2d()), (get_yval(i) + 0.5*DY_VAL*rand2d()), 0.0);
			// Note: pos.z will be slightly different when calculated within vs. outside the mesh bounds
			pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
			if (pos.z > max_tree_h || pos.z < min_tree_h) continue;
			if (tree_mode == 3 && get_tree_class_from_height(pos.z, 0) != TREE_CLASS_DECID) continue; // use a pine tree here (or no tree)
			int ttype(-1), tree_id(-1);

			if (NONUNIFORM_TREE_DEN) {
				float const dist_test(get_rel_height(density_gen[0].eval_index(j-x1, i-y1, 1), -zmax_est, zmax_est));
				if (dist_test > TREE_DEN_THRESH) continue; // density function test
				float max_val(0.0);

				for (unsigned tt = 0; tt < NUM_TREE_TYPES; ++tt) {
					float const den_val(density_gen[tt+1].eval_index(j-x1, i-y1, 1));
					if (max_val == 0.0 || den_val > max_val) {max_val = den_val; ttype = tt;}
				}
				max_val = get_rel_height(max_val, -zmax_est, zmax_est);
			}
			if (!shared_tree_data.empty()) {
				if (ttype >= 0) {
					unsigned const num_per_type(max(1U, shared_tree_data.size()/NUM_TREE_TYPES));
					tree_id = min(unsigned((global_rand_gen.rseed2 % num_per_type) + ttype*num_per_type), shared_tree_data.size()-1);
					if (shared_tree_data[tree_id].is_created()) {ttype = shared_tree_data[tree_id].get_tree_type();} // in case there weren't enough generated to get the requested type
				}
				else {
					tree_id = (global_rand_gen.rseed2 % shared_tree_data.size());
					ttype   = tree_id % NUM_TREE_TYPES;
				}
				//cout << "selected tree " << tree_id << " of " << shared_tree_data.size() << " type " << ttype << endl;
			}
			push_back(tree());
			if (tree_id >= 0) {back().bind_to_td(&shared_tree_data[tree_id]);}
			back().gen_tree(pos, 0, ttype, 1, 1, 0);
		}
	}
	//sort(begin(), end()); // doesn't help
	generated = 1;
}


void regen_trees(bool recalc_shadows, bool keep_old) {

	cout << "vegetation: " << vegetation << endl;
	RESET_TIME;
	calc_leaf_points();
	static int init(0), last_rgi(0), last_xoff2(0), last_yoff2(0);
	static float last_ts(0.0);
	if (tree_mode && recalc_shadows) {reset_shadows(OBJECT_SHADOW);}
	
	if (tree_mode & 2) {
		gen_small_trees();
	}
	else {
		//remove_small_tree_cobjs();
	}
	if ((tree_mode & 1) && num_trees > 0) {
		if (keep_old && init && last_rgi == rand_gen_index && last_xoff2 == xoff2 && last_yoff2 == yoff2 && last_ts == tree_scale)
		{ // keep old trees
			if (recalc_shadows) calc_visibility(SUN_SHADOW | MOON_SHADOW | TREE_ONLY);
			add_tree_cobjs();
			PRINT_TIME(" gen tree fast");
			return;
		}
		int const ext_x1(1), ext_x2(MESH_X_SIZE-1), ext_y1(1), ext_y2(MESH_Y_SIZE-1);

		if (scrolling && t_trees.scroll_trees(ext_x1, ext_x2, ext_y1, ext_y2)) {
			t_trees.post_scroll_remove();
		}
		else {
			t_trees.resize(0);
		}
		PRINT_TIME(" Delete Trees");
		t_trees.gen_deterministic(ext_x1, ext_y1, ext_x2, ext_y2, vegetation);
		if (!scrolling) {cout << "Num trees = " << t_trees.size() << endl;}
		last_rgi   = rand_gen_index;
		last_xoff2 = xoff2;
		last_yoff2 = yoff2;
		last_ts    = tree_scale;
		init       = 1;
	}
	if (recalc_shadows) calc_visibility(SUN_SHADOW | MOON_SHADOW | TREE_ONLY);
	PRINT_TIME(" Gen Trees");
}


void tree_data_manager_t::ensure_init() {

	if (max_unique_trees > 0 && empty()) {
		resize(max_unique_trees);
	}
	else if (tree_scale != last_tree_scale || rand_gen_index != last_rgi) {
		for (iterator i = begin(); i != end(); ++i) {i->clear_data();}
		last_tree_scale = tree_scale;
		last_rgi        = rand_gen_index;
	}
}

void tree_data_manager_t::clear_vbos() {
	for (iterator i = begin(); i != end(); ++i) {i->clear_vbos();}
}

unsigned tree_data_manager_t::get_gpu_mem() const {
	unsigned mem(0);
	for (const_iterator i = begin(); i != end(); ++i) {mem += i->get_gpu_mem();}
	return mem;
}


unsigned tree_cont_t::get_gpu_mem() const {
	unsigned mem(0);
	for (const_iterator i = begin(); i != end(); ++i) {mem += i->get_gpu_mem();}
	return mem;
}

float tree_cont_t::get_rmax() const {
	float rmax(0.0);
	for (const_iterator i = begin(); i != end(); ++i) {rmax = max(rmax, i->get_radius());}
	return rmax;
}

void tree_cont_t::update_zmax(float &tzmax) const {
	for (const_iterator i = begin(); i != end(); ++i) {tzmax = max(tzmax, (i->get_center().z + i->get_radius()));}
}

void tree_cont_t::shift_by(vector3d const &vd) {
	for (iterator i = begin(); i != end(); ++i) {i->shift_tree(vd);}
}

void tree_cont_t::add_cobjs() {
	for (iterator i = begin(); i != end(); ++i) {i->add_tree_collision_objects();}
}

void tree_cont_t::clear_vbos() {
	for (iterator i = begin(); i != end(); ++i) {i->clear_vbo();}
}


void shift_trees(vector3d const &vd) {
	if (num_trees > 0) return; // dynamically created, not placed
	t_trees.shift_by(vd);
}

void add_tree_cobjs   () {t_trees.add_cobjs();}
void remove_tree_cobjs() {t_trees.remove_cobjs();}

void clear_tree_vbos() {
	t_trees.clear_vbos();
	tree_data_manager.clear_vbos();
}




