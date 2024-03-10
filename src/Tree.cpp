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
#include "sinf.h"
#include "cobj_bsp_tree.h"
#include "draw_utils.h"

float const BURN_RADIUS      = 0.2;
float const BURN_DAMAGE      = 80.0;
float const BEAM_DAMAGE      = 0.06;
float const LEAF_DAM_SCALE   = 0.0001;
int   const LEAF_HEAL_RATE   = 25; // Note: in units of 1000; only applied every 64 frames
float const TREE_SIZE        = 0.005;
int   const LEAF_GEN_RAND1   = 16;
int   const LEAF_GEN_RAND2   = 200000; // larger is fewer leaves falling
float const DIST_C_SCALE     = 0.01;
float const BASE_LEN_SCALE   = 0.8; // determines base cylinder overlap
float const BRANCH_LEN_SCALE = 0.9; // determines branch cylinder overlap
float const REL_LEAF_SIZE    = 3.5;
float const LEAF_4TH_SCALE   = 0.4;
int const TREE_4TH_LEAVES    = 1;
int const DISABLE_LEAVES     = 0;
int const ENABLE_CLIP_LEAVES = 1;
int const TLEAF_START_TUID   = 8; // trees use texture units 8-12
bool const FORCE_TREE_TYPE   = 1;
unsigned const CYLINS_PER_ROOT     = 3;
unsigned const TREE_BILLBOARD_SIZE = 256;


// bark_tex, leaf_tex, branch_size, branch_radius, leaf_size, leaf_x_ar, height_scale, branch_break_off, branch_tscale, branch_color_var, bush_prob, barkc, leafc
tree_type tree_types[NUM_TREE_TYPES] = { // not const - config file can override parameters (such as bush_prob)
	tree_type(BARK3_TEX, LEAF_TEX,     1.0, 0.7, 1.0, 1.00, 1.0, 1.0, 1.0, 0.1,  0.0, colorRGBA(0.7, 0.7,  0.5,  1.0), WHITE),
	tree_type(BARK4_TEX, LIVE_OAK_TEX, 1.0, 1.0, 1.0, 0.63, 1.0, 1.0, 1.5, 0.1,  0.0, colorRGBA(1.0, 0.9,  0.8,  1.0), WHITE),
	tree_type(BARK1_TEX, LEAF2_TEX,    1.0, 1.0, 1.0, 0.82, 2.0, 0.5, 1.0, 0.1,  0.0, colorRGBA(0.8, 0.5,  0.3,  1.0), WHITE),
	tree_type(BARK5_TEX, LEAF3_TEX,    1.0, 0.7, 1.5, 0.81, 1.0, 1.0, 0.3, 0.01, 0.0, colorRGBA(0.8, 0.75, 0.65, 1.0), WHITE), // birch bark
	tree_type(BARK6_TEX, PAPAYA_TEX,   1.0, 1.0, 1.0, 1.00, 2.0, 2.0, 0.5, 0.1,  0.0, colorRGBA(0.7, 0.6,  0.5,  1.0), WHITE)
};

vector<tree_cylin >   tree_builder_t::cylin_cache;
vector<tree_branch>   tree_builder_t::branch_cache;
vector<tree_branch *> tree_builder_t::branch_ptr_cache;


bool has_any_billboard_coll(0), next_has_any_billboard_coll(0), tree_4th_branches(0);
unsigned max_unique_trees(0);
int tree_mode(1), tree_coll_level(2); // tree_mode: 0 = no trees, 1 = large only, 2 = small only, 3 = both large and small
float leaf_color_coherence(0.5), tree_color_coherence(0.2), tree_deadness(-1.0), tree_dead_prob(0.0), nleaves_scale(1.0), branch_radius_scale(1.0), tree_height_scale(1.0);
float tree_lod_scales[4] = {0, 0, 0, 0}; // branch_start, branch_end, leaf_start, leaf_end
colorRGBA leaf_base_color(BLACK);
tree_data_manager_t tree_data_manager;
tree_cont_t t_trees(tree_data_manager);
tree_placer_t tree_placer;


extern bool has_snow, no_sun_lpos_update, has_dl_sources, gen_tree_roots, tt_lightning_enabled, tree_indir_lighting, begin_motion, enable_grass_fire, rotate_trees;
extern int num_trees, do_zoom, display_mode, animate2, iticks, draw_model, frame_counter;
extern int xoff2, yoff2, rand_gen_index, leaf_color_changed, scrolling, dx_scroll, dy_scroll, window_width, window_height;
extern unsigned smoke_tid;
extern float zmin, zmax, zmax_est, zbottom, water_plane_z, tree_scale, temperature, fticks, vegetation, tree_density_thresh, tree_slope_thresh;
extern double sim_ticks;
extern vector3d wind;
extern lightning_t l_strike;
extern coll_obj_group coll_objects;

void set_indir_color(shader_t &s);
bool is_csm_active();


bool enable_rotate_trees() {return (rotate_trees && world_mode == WMODE_INF_TERRAIN);}

inline colorRGBA get_leaf_base_color(int type) {
	colorRGBA color(tree_types[type].leafc);
	UNROLL_3X(color[i_] = CLIP_TO_01(color[i_] + leaf_base_color[i_]);)
	return color;
}
colorRGBA get_leaf_texture_color(unsigned type) {
	// alpha is always 1.0 - texture alpha is handled by check_poly_billboard_alpha()
	return colorRGBA(texture_color(tree_types[type].leaf_tex), 1.0);
}
colorRGBA get_avg_leaf_color(unsigned type) {
	return get_leaf_base_color(type).modulate_with(get_leaf_texture_color(type));
}


// better place for this?
struct render_to_texture_shader_t : public render_to_texture_t {

	shader_t shaders[2]; // color, normal

	render_to_texture_shader_t(unsigned tsize_) : render_to_texture_t(tsize_) {}
	~render_to_texture_shader_t() {free_context();}
	
	void free_context() {
		for (unsigned d = 0; d < 2; ++d) {shaders[d].end_shader();}
	}
	void setup_shader(string const &vs, string const &fs, bool ix) {
		shaders[ix].set_vert_shader(vs);
		shaders[ix].set_frag_shader(fs);
		shaders[ix].begin_shader(1);
		shaders[ix].add_uniform_float("min_alpha", 0.75);
		shaders[ix].add_uniform_int("tex0", 0);
		set_other_shader_consts(ix);
		shaders[ix].disable();
	}
	virtual void set_other_shader_consts(bool ix) {}
};

struct render_tree_to_texture_t : public render_to_texture_shader_t {

	tree_data_t *cur_tree;
	render_tree_to_texture_t(unsigned tsize_) : render_to_texture_shader_t(tsize_), cur_tree(NULL) {}

	tree_type const &get_tree_type() const {
		assert(cur_tree);
		return tree_types[cur_tree->get_tree_type()];
	}
};

struct render_tree_leaves_to_texture_t : public render_tree_to_texture_t {

	render_tree_leaves_to_texture_t(unsigned tsize_) : render_tree_to_texture_t(tsize_) {}

	void set_other_shader_consts(bool ix) {
		shaders[ix].add_uniform_int("tc_start_ix", 3);
	}
	virtual void draw_geom(bool is_normal_pass) {
		select_texture(get_tree_type().leaf_tex);
		shader_t &s(shaders[is_normal_pass]);
		s.enable();
		tree_data_t::pre_leaf_draw(s);
		cur_tree->gen_leaf_color();
		cur_tree->leaf_draw_setup(1);
		cur_tree->draw_leaves(0.0); // sets color
		tree_data_t::post_leaf_draw();
		s.disable();
	}
	void render_tree(tree_data_t &t, tree_bb_tex_t &ttex, vector3d const &view_dir=plus_y) {
		if (!shaders[0].is_setup()) {setup_shader("texture_gen.part+tree_leaves_no_lighting", "simple_texture",        0);} // colors
		if (!shaders[1].is_setup()) {setup_shader("texture_gen.part+tree_leaves_no_lighting", "write_normal_textured", 1);} // normals
		cur_tree = &t;
		colorRGBA const leaf_bkg_color(get_avg_leaf_color(t.get_tree_type()), 0.0); // transparent
		bool const use_depth_buffer(1), mipmap(0); // Note: for some reason mipmaps are slow and don't look any better
		render(ttex, t.lr_x, t.lr_z, vector3d(0.0, 0.0, t.lr_z_cent), view_dir, leaf_bkg_color, use_depth_buffer, mipmap);
	}
};

struct render_tree_branches_to_texture_t : public render_tree_to_texture_t {

	render_tree_branches_to_texture_t(unsigned tsize_) : render_tree_to_texture_t(tsize_) {}

	virtual void draw_geom(bool is_normal_pass) {
		shader_t &s(shaders[is_normal_pass]);
		s.enable();
		s.set_cur_color(WHITE); // branch color will be applied to the billboard later
		select_texture(get_tree_type().bark_tex);
		tree_data_t::pre_branch_draw(s, 0);
		cur_tree->draw_branches(0.0, 0);
		tree_data_t::post_branch_draw(0);
		s.disable();
	}
	void render_tree(tree_data_t &t, tree_bb_tex_t &ttex, vector3d const &view_dir=plus_y) {
		if (!shaders[0].is_setup()) {setup_shader("no_lighting_tex_coord",     "simple_texture",        0);} // colors
		if (!shaders[1].is_setup()) {setup_shader("tree_branches_no_lighting", "write_normal_textured", 1);} // normals
		cur_tree = &t;
		t.ensure_branch_vbo(); // Note: for some reason, this *must* be called before we get into draw_geom()
		colorRGBA const branch_bkg_color(texture_color(get_tree_type().bark_tex), 0.0); // transparent
		bool const use_depth_buffer(1), mipmap(0); // Note: for some reason mipmaps are slow and don't look any better
		render(ttex, t.br_x, t.br_z, t.get_center(), view_dir, branch_bkg_color, use_depth_buffer, mipmap);
	}
};


void tree_lod_render_t::finalize() {
	sort(leaf_vect.begin(),   leaf_vect.end());
	sort(branch_vect.begin(), branch_vect.end());
}

void tree_lod_render_t::render_billboards(shader_t &s, bool render_branches) const { // branches or leaves

	vector<entry_t> const &data(render_branches ? branch_vect : leaf_vect);
	if (data.empty()) return;
	tree_data_t const *last_td(nullptr);
	s.add_uniform_vector3d("camera_pos", get_camera_pos());
	s.add_uniform_vector3d("up_vector",  up_vector);
	static vector<vert_tc_color> pts; // reused across frames

	for (vector<entry_t>::const_iterator i = data.begin(); i != data.end(); ++i) {
		if (i->td != last_td) {
			assert(i->td);
			last_td = i->td;
			draw_and_clear_verts(pts, GL_POINTS);
			(render_branches ? i->td->get_render_branch_texture() : i->td->get_render_leaf_texture()).bind_texture();
		}
		pts.emplace_back(i->pos, (render_branches ? i->td->br_x : i->td->lr_x), (render_branches ? i->td->br_z : i->td->lr_z), i->cw.c);
	} // for i
	assert(!pts.empty());
	draw_and_clear_verts(pts, GL_POINTS);
	pts.clear();
	bind_vbo(0);
}


float get_tree_z_bottom(float z, point const &pos) {
	return ((world_mode == WMODE_GROUND && is_over_mesh(pos)) ? max(zbottom, (z - TREE_DEPTH)) : (z - TREE_DEPTH));
}


bool tree::is_over_mesh() const {

	//if (world_mode != WMODE_GROUND) return 1; // ???
	float const r(tdata().sphere_radius);
	int const x1(get_xpos(tree_center.x - r)), y1(get_ypos(tree_center.y - r));
	int const x2(get_xpos(tree_center.x + r)), y2(get_ypos(tree_center.y + r));
	return (x1 < MESH_X_SIZE && y1 < MESH_Y_SIZE && x2 >= 0 && y2 >= 0); // completely off the mesh
}


void tree::add_tree_collision_objects() {

	//RESET_TIME;
	if (!tree_coll_level || !physics_enabled()) return;
	remove_collision_objects();
	if (!is_over_mesh()) return; // optimization
	assert(type < NUM_TREE_TYPES);
	int const btid(tree_types[type].bark_tex), branch_coll_level(min(tree_coll_level, 4));
	cobj_params cp(0.8, tree_types[type].barkc, 0, 0, NULL, 0, btid, 4.0, 1, 0);
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
		colorRGBA const lcolor(get_avg_leaf_color(type)); // will be reset in update_leaf_cobj_color()
		cobj_params cpl(0.3, lcolor, 0, 0, NULL, 0, ltid, 1.0, 0, 0);
		cpl.flags |= COBJ_DESTROYABLE; // so that truly_static() returns false
		point const xlate(all_zeros); // for now
		vector<tree_leaf> const &leaves(tdata().get_leaves());
		leaf_cobjs.resize(leaves.size());

		for (unsigned i = 0; i < leaves.size(); i++) { // loop through leaves
			// Note: line collisions with leaves will use the texture alpha component for a more exact test
			point pts[4];
			get_abs_leaf_pts(pts, i);
			for (int p = 0; p < 4; ++p) {pts[p] += xlate;}
			leaf_cobjs[i] = add_simple_coll_polygon(pts, 4, cpl, leaves[i].norm, 2);
			coll_objects[leaf_cobjs[i]].is_billboard = 1;
		}
	}
	//PRINT_TIME("Tree Cobjs");
}


void tree::remove_collision_objects() {
	for (unsigned i = 0; i < branch_cobjs.size(); i++) {remove_reset_coll_obj(branch_cobjs[i]);}
	for (unsigned i = 0; i < leaf_cobjs  .size(); i++) {remove_reset_coll_obj(leaf_cobjs  [i]);}
	branch_cobjs.clear();
	leaf_cobjs  .clear();
}

void tree_cont_t::remove_cobjs() {
	for (iterator i = begin(); i != end(); ++i) {i->remove_collision_objects();}
}


bool tree::check_sphere_coll(point &center, float radius) const { // collision with the trunk
	// 10.7% of CPU time, 4.5x slower than scenery coll
	tree_data_t const &td(tdata());
	if (!td.branches_bcube.contains_pt_exp((center - tree_center), radius)) return 0; // optimization
	float const trunk_radius(0.9*td.br_scale*td.base_radius);
	float const trunk_height(max(td.sphere_radius, td.sphere_center_zoff)); // very approximate
	cylinder_3dw const cylin(tree_center, tree_center+vector3d(0.0, 0.0, trunk_height), trunk_radius, trunk_radius);
	return sphere_vert_cylin_intersect(center, radius, cylin);
}
bool tree_cont_t::check_sphere_coll(point &center, float radius) const {
	if (!all_bcube.is_zero_area() && !sphere_cube_intersect(center, radius, all_bcube)) return 0; // optimization
	bool coll(0);
	for (const_iterator i = begin(); i != end(); ++i) {coll |= i->check_sphere_coll(center, radius);}
	return coll;
}

bool tree::check_cube_int(cube_t const &c) const {
	cube_t const c_tree_rel(c - tree_center);
	tree_data_t const &td(tdata());
	if (!td.leaves_bcube.intersects(c_tree_rel) && !td.branches_bcube.intersects(c_tree_rel)) return 0;
	return sphere_cube_intersect(td.get_center(), td.sphere_radius, c_tree_rel);
}
bool tree_cont_t::check_cube_int(cube_t const &c) const {
	if (!all_bcube.is_zero_area() && !all_bcube.intersects(c)) return 0; // optimization
	
	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->check_cube_int(c)) return 1;
	}
	return 0;
}


int tree_cont_t::draw_branches_and_leaves(shader_t &s, tree_lod_render_t &lod_renderer,
	bool draw_branches, bool draw_leaves, bool shadow_only, bool reflection_pass, vector3d const &xlate)
{
	assert(draw_branches != draw_leaves); // must enable exactly one
	if (!all_bcube.is_zero_area() && !camera_pdu.cube_visible(all_bcube + xlate)) return -1; // VFC
	int const wsoff_loc(shadow_only ? (is_csm_active() ? s.get_uniform_loc("xlate_scale") : -1) : s.get_uniform_loc("world_space_offset")); // not needed in non-CSM shadow mode
	bool const tt_shadow_mode(world_mode == WMODE_INF_TERRAIN && shadow_only);

	if (draw_branches) {
		tree_data_t::pre_branch_draw(s, shadow_only);

		for (iterator i = begin(); i != end(); ++i) {
			if (tt_shadow_mode && !dist_less_than((i->sphere_center() + xlate), camera_pdu.pos, camera_pdu.far_)) continue; // Note: intentionally excludes tree radius
			i->draw_branches_top(s, lod_renderer, shadow_only, reflection_pass, xlate, wsoff_loc);
		}
		tree_data_t::post_branch_draw(shadow_only);
	}
	else { // draw_leaves
		to_update_leaves.clear();
		tree_data_t::pre_leaf_draw(s);

		if (tt_shadow_mode) {
			for (tree &t : *this) {
				point const center(t.sphere_center() + xlate);
				// still need VFC in tiled terrain mode since the shadow volume doesn't include the entire scene (especially for local city light shadows)
				if (!dist_less_than(center, camera_pdu.pos, camera_pdu.far_))    continue; // Note: intentionally excludes tree radius
				if (!camera_pdu.sphere_visible_test(center, 1.1*t.get_radius())) continue;
				t.draw_leaves_top(s, lod_renderer, shadow_only, reflection_pass, xlate, wsoff_loc, -1, to_update_leaves);
			}
			tree_data_t::post_leaf_draw();
		}
		else {
			int const tex0_loc(shadow_only ? -1 : s.get_uniform_loc("tex0")); // not needed in shadow mode
			point const local_camera(get_camera_pos() - xlate);
			
			if (shadow_only || (world_mode == WMODE_INF_TERRAIN && !all_bcube.is_all_zeros() &&
				!all_bcube.closest_dist_less_than(local_camera, 1.0*(X_SCENE_SIZE + Y_SCENE_SIZE))))
			{ // direct draw
				for (tree &t : *this) {t.draw_leaves_top(s, lod_renderer, shadow_only, reflection_pass, xlate, wsoff_loc, tex0_loc, to_update_leaves);}
			}
			else { // sorted draw - when close to the player
				sorted.clear();
				for (unsigned i = 0; i < size(); ++i) {sorted.emplace_back(p2p_dist_sq(operator[](i).sphere_center(), local_camera), i);}
				sort(sorted.begin(), sorted.end()); // sort front to back for better early z culling

				for (auto i = sorted.begin(); i != sorted.end(); ++i) {
					operator[](i->second).draw_leaves_top(s, lod_renderer, shadow_only, reflection_pass, xlate, wsoff_loc, tex0_loc, to_update_leaves);
				}
			}
			tree_data_t::post_leaf_draw();
			int const num_to_update(to_update_leaves.size());
#pragma omp parallel for num_threads(max(1, min(4, num_to_update))) schedule(static) if (num_to_update > 1)
			for (int i = 0; i < num_to_update; ++i) {to_update_leaves[i]->update_leaf_orients_wind();}
		}
	}
	return wsoff_loc;
}


float get_plant_leaf_wind_mag(bool shadow_only) {
	//if (shadow_only || animate2) return 0.0; // faster, but looks odd
	return (has_snow ? 0.0 : 0.001*min(2.0f, wind.mag())/tree_scale); // Note: animate2 is 0 for shadow pass, so can't check it
}

void setup_leaf_wind(shader_t &s, float wind_mag, bool underwater) {

	if (wind_mag == 0.0) return;
	s.add_uniform_float("wind_mag",   wind_mag);
	s.add_uniform_float("wind_scale", 1.0);
	s.add_uniform_float("wind_time",  (underwater ? 0.025 : 0.1)*sim_ticks); // lower frequency movement for underwater seaweed (but applies to all plants)
	s.add_uniform_float("wind_freq",  80.0*tree_scale);
}

void set_leaf_shader(shader_t &s, float min_alpha, unsigned tc_start_ix, bool enable_opacity, bool no_dlights, float wind_mag,
	bool underwater, bool use_fs_smap, bool enable_smap, bool enable_tex_coord_weight, bool shadow_only)
{
	cube_t lights_bcube(all_zeros);

	if (world_mode == WMODE_INF_TERRAIN && !shadow_only) {
		if (!no_dlights) {
			lights_bcube = get_city_lights_bcube();
			if (lights_bcube.is_all_zeros()) {no_dlights = 1;}
		}
		setup_tt_fog_pre(s); // FS
	}
	int const shader_type(use_fs_smap ? 1 : 0); // VS/FS (for lighting)
	float const water_depth(setup_underwater_fog(s, 0)); // VS
	bool const use_indir(tree_indir_lighting && smoke_tid && !shadow_only);
	bool const use_smap(enable_smap && shadow_map_enabled());
	bool const for_trees(tc_start_ix == 3); // hack
	s.set_prefix(make_shader_bool_prefix("indir_lighting", use_indir), shader_type);
	if (wind_mag > 0.0)          {s.set_prefix("#define ENABLE_WIND",             0);} // VS
	if (enable_tex_coord_weight) {s.set_prefix("#define ENABLE_TEX_COORD_WEIGHT", 0);} // VS
	s.check_for_fog_disabled();
	s.setup_enabled_lights(2, (1<<shader_type));
	s.set_prefix(make_shader_bool_prefix("enable_light2", (world_mode == WMODE_INF_TERRAIN && tt_lightning_enabled)), shader_type);
	set_dlights_booleans(s, !no_dlights, shader_type, 1); // no_dl_smap=1
	s.set_prefix(make_shader_bool_prefix("use_shadow_map", use_smap), shader_type);
	if (use_smap) {s.set_prefix("#define ENABLE_LEAF_SMAP", shader_type);} // need to set this as well to avoid even using the shadow texture in the shader on ATI cards
	if (enable_rotate_trees() && for_trees) {s.set_prefix("#define ENABLE_ROTATIONS", 0);} // VS
	string const ls_str("ads_lighting.part*+shadow_map.part*+dynamic_lighting.part*+leaf_lighting_comp.part*+tree_leaf_lighting.part*");

	if (shader_type == 0) { // VS
		s.set_frag_shader(enable_opacity ? "linear_fog.part+noise_dither.part+textured_with_fog_opacity" : "linear_fog.part+textured_with_fog");
		s.set_vert_shader(ls_str + "+world_space_offset_rot.part+leaf_lighting.part+texture_gen.part+leaf_wind.part+tree_leaves");
	}
	else { // FS
		if (enable_opacity) {s.set_prefix("#define ENABLE_OPACITY", 1);} // FS
		s.set_frag_shader(ls_str + (enable_opacity ? "+noise_dither.part" : "") + "+linear_fog.part+leaf_lighting_ppl");
		s.set_vert_shader("texture_gen.part+leaf_wind.part+world_space_offset_rot.part+tree_leaves_ppl");
	}
	s.begin_shader();
	s.setup_scene_bounds(); // also sets camera_pos

	if (!shadow_only) {
		if (world_mode == WMODE_GROUND && use_smap) {set_smap_shader_for_all_lights(s);}
		if (!no_dlights) {setup_dlight_textures(s, 0);} // no dlight smap
		if (world_mode == WMODE_INF_TERRAIN) {setup_tt_fog_post(s);} else {s.setup_fog_scale();}
		if (!lights_bcube.is_all_zeros()) {set_city_lighting_shader_opts(s, lights_bcube, !no_dlights, use_smap);} // will reset some values
	}
	if (enable_tex_coord_weight) {s.add_uniform_float("tex_coord_weight", 0.0);}
	s.add_uniform_float("water_depth", water_depth);
	s.add_uniform_float("min_alpha",   min_alpha);
	s.add_uniform_int("tex0", 0);
	s.add_uniform_int("tc_start_ix", tc_start_ix);
	s.add_uniform_vector4d("world_space_offset", vector4d()); // reset

	if (use_indir) {
		bind_texture_tu(smoke_tid, 1);
		s.add_uniform_int("smoke_and_indir_tex", 1);
		set_indir_color(s);
	}
	setup_leaf_wind(s, wind_mag, underwater);
	check_gl_error(301);
}


void tree_cont_t::pre_leaf_draw(shader_t &shader, bool enable_opacity, bool shadow_only, bool use_fs_smap, bool enable_smap, bool enable_dlights) {
	
	if (shader.is_setup()) {shader.enable();}
	else { // Note: disabling leaf wind when shadow_only is faster but looks odd
		float const wind_mag((has_snow || !animate2 /*|| shadow_only*/) ? 0.0f : 0.05f*REL_LEAF_SIZE*TREE_SIZE/(sqrt(nleaves_scale)*tree_scale)*min(2.0f, wind.mag()));
		if (enable_smap) {shader.set_prefix("#define NO_SHADOW_PCF", (use_fs_smap ? 1 : 0));} // faster shadows
		set_leaf_shader(shader, 0.75, 3, enable_opacity, (shadow_only || !enable_dlights), wind_mag, 0, use_fs_smap, enable_smap, 0, shadow_only); // no underwater trees
		for (int i = 0; i < NUM_TREE_TYPES; ++i) {select_texture(((draw_model == 0) ? tree_types[i].leaf_tex : WHITE_TEX), TLEAF_START_TUID+i);}
	}
	if (!shadow_only) {
		shader.set_specular(0.2, 20.0); // small amount of specular
		set_multisample(0); // disable AA to prevent bright pixel artifacts and to improve frame rate
	}
}

void tree_cont_t::post_leaf_draw(shader_t &shader, bool shadow_only) {
	if (!shadow_only) {
		set_multisample(1); // re-enable AA
		shader.clear_specular();
	}
	shader.disable();
}


void tree_cont_t::draw(bool shadow_only, bool reflection_pass) {

	if (empty()) return;
	if (!all_bcube.is_zero_area() && !camera_pdu.cube_visible(all_bcube)) return;
	tree_lod_render_t lod_renderer(0); // disabled

	// draw leaves, then branches: much faster for distant trees, slightly slower for near trees
	shader_t ls;
	pre_leaf_draw(ls, 0, shadow_only);
	draw_branches_and_leaves(ls, lod_renderer, 0, 1, shadow_only, reflection_pass, zero_vector);
	post_leaf_draw(ls, shadow_only);

	// draw branches
	shader_t bs;
	if (shadow_only) {bs.begin_shadow_map_shader(0, 1);} // use_alpha_mask=0, enable_xlate_scale=1
	else {set_tree_branch_shader(bs, 1, 1, 1);} // direct_lighting=1, dlights=1, smap=1
	int const wsoff_loc(draw_branches_and_leaves(bs, lod_renderer, 1, 0, shadow_only, reflection_pass, zero_vector));
	if (!shadow_only) {bs.set_uniform_vector4d(wsoff_loc, vector4d());} // reset
	else              {bs.set_uniform_vector4d(wsoff_loc, vector4d(0.0, 0.0, 0.0, 1.0));} // needed for CSMs
	bs.end_shader();
}


void draw_trees(bool shadow_only, bool reflection_pass) {

	if (tree_mode & 2) {draw_small_trees(shadow_only, reflection_pass);} // small trees

	if (tree_mode & 1) { // trees
		if (shadow_only) {t_trees.draw(1, reflection_pass);}
		else {
			next_has_any_billboard_coll = 0; // reset for this frame
			t_trees.draw(0, reflection_pass);
			has_any_billboard_coll = next_has_any_billboard_coll; // keep only if some leaf already has a coll
		}
	}
	if (!shadow_only) {leaf_color_changed = 0;}
}

void register_leaf_color_change() {
	leaf_color_changed = 1;
	tree_data_manager.on_leaf_color_change();
}


void tree_data_t::make_private_copy(tree_data_t &dest) const {
	dest = *this;
	dest.clear_vbo_ixs();
}


void tree::make_private_tdata_copy() {

	if (!tree_data || td_is_private()) return; // tree pointer is NULL or already private
	tree_data->make_private_copy(priv_tree_data);
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
	return coll_objects.get_cobj(leaf_cobjs[i]);
}


void tree_data_t::mark_leaf_changed(unsigned ix) {

	assert(ix < leaves.size());
	leaf_change_start = min(ix,   leaf_change_start);
	leaf_change_end   = max(ix+1, leaf_change_end);
}


void tree_data_t::gen_leaf_color() {leaf_color = get_leaf_base_color(tree_type)*leaf_color_coherence;}

inline colorRGB tree_leaf::calc_leaf_color(colorRGBA const &leaf_color, colorRGBA const &base_color) const {

	float const color(lcolor/1000.0), ilch((1.0 - leaf_color_coherence)/255.0); // scale by 255 to convert uchar to 0-1 range
	return colorRGB(max(0.0f, (color*(leaf_color.R + ilch*lred  ) + base_color.R*tree_color_coherence)),
		            max(0.0f, (color*(leaf_color.G + ilch*lgreen) + base_color.G*tree_color_coherence)), 0.0);
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
	color_wrapper cw(leaves[i].calc_leaf_color(leaf_color, base_color));
	UNROLL_4X(leaf_data[i_+(i<<2)].copy_color(cw);)
	if (!no_mark_changed) {mark_leaf_changed(i);}
}

void tree_data_t::update_all_leaf_colors() {
	for (unsigned i = 0; i < leaves.size(); i++) {update_leaf_color(i);}
}


void tree::update_leaf_cobj_color(unsigned i) {

	// update cobj color so that leaf water reflection is correct (Note: not always correct with tree instancing)
	get_leaf_cobj(i).cp.color = colorRGBA(tdata().get_leaf_color(i)).modulate_with(get_leaf_texture_color(type));
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
	if (i < leaves.size()) {mark_leaf_changed(i);} // not the last leaf
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


bool tree::spraypaint_leaves(point const &pos, float radius, colorRGBA const &color) {

	if (!td_is_private()) { // using shared tdata, can't just modify the leaves, need to make a private copy if any leaf colors will be modified
		if (!tdata().spraypaint_leaves((pos - tree_center), radius, color, 1)) return 0; // check returns false, we're done
		make_private_tdata_copy(); // check returned true, make a private copy then go back and modify the leaf colors
	}
	return tdata().spraypaint_leaves((pos - tree_center), radius, color, 0); // only this call actually changes leaf colors
}


bool tree_data_t::spraypaint_leaves(point const &pos, float radius, colorRGBA const &color, bool check_only) {

	if (!leaf_data_allocated()) return 0;
	bool changed(0);

	// Note: may apply to multiple trees if instancing is used
	for (unsigned i = 0; i < leaves.size(); ++i) {
		point const center(leaves[i].get_center());
		bool const center_close(dist_less_than(pos, center, radius));

		for (unsigned j = 0; j < 4; ++j) { // use closest point - leaf corner or center
			if (!center_close && !dist_less_than(pos, leaves[i].pts[j], radius)) continue;
			if (check_only) return 1; // leaf color will be changed
			float const blend_val(color.alpha*(1.0 - min(p2p_dist(pos, leaves[i].pts[j]), p2p_dist(pos, center))/radius));
			unsigned const ix((i<<2)+j);
			assert(ix < leaf_data.size());
			UNROLL_3X(leaf_data[ix].c[i_] = (unsigned char)((1.0 - blend_val)*leaf_data[ix].c[i_] + 255.0*blend_val*CLIP_TO_01(color[i_]));)
			mark_leaf_changed(i);
			changed = 1;
		}
	}
	return changed;
}


void tree::burn_leaves() {

	if (!physics_enabled()) return;
	tree_data_t &td(tdata());
	vector<tree_leaf> &leaves(td.get_leaves());
	float const max_t(get_max_t(LEAF));
	if (temperature < max_t || leaves.empty()) return;
	float const burn_amt(0.25*((td_is_private() || t_trees.empty()) ? 1.0 : float(tree_data_manager.size())/float(t_trees.size())));
	int const burn_amt_scaled(1000*burn_amt);
	unsigned const num_burn(max(1U, min(5U, unsigned(5*(temperature - max_t)/max_t))));
	damage += ((1.0 - damage)*num_burn)/leaves.size();

	for (unsigned i = 0; i < num_burn && !leaves.empty(); ++i) {
		unsigned const index(rand()%leaves.size());
		leaves[index].lcolor = max(0, (leaves[index].lcolor - burn_amt_scaled));
		copy_color(index);
		if (rand()&1) {gen_smoke(leaves[index].pts[0] + tree_center);}
		if (td_is_private() && leaves[index].lcolor <= 0) {remove_leaf(index, 1);} // Note: if we modify shared data, leaves.size() must be dynamic
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


bool tree::damage_leaf(unsigned i, float damage_done, rand_gen_t &rgen) {

	make_private_tdata_copy();
	vector<tree_leaf> &leaves(tdata().get_leaves());
	assert(i < leaves.size());
	short &lcolor(leaves[i].lcolor);

	if (damage_done > 4.0 || (lcolor == 0 && (damage_done > 0.4 || (rgen.rand()&31) == 0))) {
		if (lcolor > 0) {damage += damage_scale;}
		lcolor = -1000; // tag as destroyed
		if ((rgen.rand()&3) == 0) {tdata().update_leaf_color(i); create_leaf_obj(i);} // 25% chance of the burned leaf falling
		
		if (has_leaf_data()) { // remove leaf i
			remove_leaf(i, 1);
			return 1;
		}
	}
	else {
		short const llc(lcolor);
		lcolor = max(0, (lcolor - int(300.0*damage_done)));
		if (lcolor == 0 && llc > 0) {damage += damage_scale;} // turned pure black
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
	float const radius(bradius + tdata().sphere_radius);
	if (p2p_dist_sq(bpos, sphere_center()) > radius*radius) return;
	burn_leaves_within_radius(bpos, bradius, bdamage);
	add_fire(bpos, 0.25*bradius, blast_radius->damage);
}

void tree::burn_leaves_within_radius(point const &bpos, float bradius, float damage) {

	if (damage == 0.0) return;
	float const bradius_sq(bradius*bradius);
	make_private_tdata_copy(); // copy before taking leaves reference, assuming we're going to modify the tree
	vector<tree_leaf> const &leaves(tdata().get_leaves());
	unsigned nleaves(leaves.size());
	point const rel_pos(bpos - tree_center);
	rand_gen_t rgen;
	rgen.set_state(leaf_burn_ix, nleaves);

	for (unsigned i = 0; i < nleaves; ++i) {
		tree_leaf const &l(leaves[i]);
		if (fabs(rel_pos.z - l.pts[0].z) > bradius) continue; // test z dist first
		if (l.lcolor < 0) continue; // already black
		float const dist_sq(p2p_dist_sq(rel_pos, l.pts[0]));
		if (dist_sq > bradius_sq || dist_sq < TOLERANCE) continue;
		if (damage_leaf(i, damage*InvSqrt(dist_sq), rgen)) {--i; --nleaves;} // force reprocess of this leaf, wraparound to -1 is OK
	} // for i
	//damage = min(1.0f, damage);
}

void tree::apply_fire_damage(vector<fire_damage_t> const &fire_damage, unsigned skipval) {

	if (fire_damage.empty()) return;
	assert(skipval > 0);
	vector<tree_leaf> const &leaves(tdata().get_leaves());
	unsigned nleaves(leaves.size());
	cube_t bcube;
	bcube.set_from_sphere(fire_damage.front());
	for (auto f = fire_damage.begin()+1; f != fire_damage.end(); ++f) {bcube.union_with_sphere(*f);}
	unsigned const block_ix((leaf_burn_ix++)%skipval), block_sz(1+nleaves/skipval), start(block_ix*block_sz), end(min(nleaves, (block_ix+1)*block_sz));
	rand_gen_t rgen;
	rgen.set_state(leaf_burn_ix, nleaves);

	for (unsigned i = start; i < min(end, nleaves); ++i) {
		tree_leaf const &l(leaves[i]);
		if (l.lcolor < 0) continue; // already black
		point const lpos(l.pts[0] + tree_center);
		if (!bcube.contains_pt(lpos)) continue;
		float damage(0.0);

		for (auto f = fire_damage.begin(); f != fire_damage.end(); ++f) {
			if (fabs(lpos.z - f->pos.z) > f->radius) continue; // test z dist first
			float const dist_sq(p2p_dist_sq(lpos, f->pos));
			if (dist_sq < f->radius*f->radius && dist_sq > TOLERANCE) {damage += f->damage*InvSqrt(dist_sq);}
		}
		if (damage > 0.0 && damage_leaf(i, damage, rgen)) {--i; --nleaves;} // force reprocess of this leaf, wraparound to -1 is OK
	}
	damage = min(1.0f, damage);
}


void tree::lightning_damage(point const &ltpos) {

	blastr const br(0, ETYPE_FIRE, NO_SOURCE, BURN_RADIUS, BURN_DAMAGE, ltpos, plus_z, LITN_C, LITN_C);
	blast_damage(&br);
}


void tree::drop_leaves() {

	if (!physics_enabled() || !enable_leaf_wind) return;
	tree_data_t &td(tdata());
	vector<tree_leaf> &leaves(td.get_leaves());
	unsigned const nleaves(leaves.size());
	if (damage >= 1.0 || nleaves == 0) return; // too damaged
	static rand_gen_t rgen;
	float const temp0(max(1.0f, min(0.3f, (20.0f-temperature)/30.0f)));
	int const rmod(min(LEAF_GEN_RAND2/10, int(rgen.rand_uniform(0.5, 1.5)*temp0*LEAF_GEN_RAND2/fticks)));
	int const rstep((td.get_has_4th_branches() ? 10 : 1)*LEAF_GEN_RAND1);

	for (unsigned i = (rgen.rand()%rstep); i < nleaves; i += rstep) {
		if ((rgen.rand()%rmod) != 0) continue;
		create_leaf_obj(i);
		
		if (td_is_private()) { // create a new leaf with a different color (and orient?)
			leaves[i].create_init_color(rgen);
			copy_color(i, 1); // don't call td.mark_leaf_changed(i) (too slow)
		}
	}
}


void tree_data_t::clear_vbo_ixs() {
	leaf_vbo = num_branch_quads = num_unique_pts = 0;
	branch_manager.reset_vbos_to_zero();
}
void tree_data_t::clear_context() {

	render_leaf_texture.free_context();
	render_branch_texture.free_context();
	branch_manager.clear_vbos();
	delete_vbo(leaf_vbo);
	clear_vbo_ixs();
}
void tree_data_t::on_leaf_color_change() {
	render_leaf_texture.free_context();
}


unsigned tree_data_t::get_gpu_mem() const {

	unsigned mem(branch_manager.gpu_mem + (leaf_vbo ? get_leaf_data_mem() : 0));
	unsigned const bbsz(TREE_BILLBOARD_SIZE*TREE_BILLBOARD_SIZE*8); // 8 bytes per pixel
	if (render_leaf_texture.is_valid  ()) {mem += bbsz;}
	if (render_branch_texture.is_valid()) {mem += bbsz;}
	return mem;
}


bool tree::is_visible_to_camera(vector3d const &xlate) const {
	int const level((get_camera_pos().z > max(ztop, czmax)) ? 0 : 2); // test cobjs and mesh unless camera is in the air
	return sphere_in_camera_view((sphere_center() + xlate), 1.1*tdata().sphere_radius, level);
}
void tree::add_bounds_to_bcube(cube_t &bcube) const {
	bcube.assign_or_union_with_cube(tdata().branches_bcube + tree_center);
	bcube.union_with_cube(tdata().leaves_bcube + tree_center);
}

void tree::shift_tree(vector3d const &vd) {
	tree_center += vd;
	if (tree_fire) {tree_fire->shift(vd);}
}

void tree_data_t::check_render_textures() {

	if (!render_leaf_texture.is_valid() && !leaves.empty()) {
#ifdef USE_TREE_BB_TEX_ATLAS
		render_leaf_texture.nx = 2;
#endif
		render_tree_leaves_to_texture_t renderer(TREE_BILLBOARD_SIZE);
		renderer.render_tree(*this, render_leaf_texture);
	}
	if (!render_branch_texture.is_valid() && !all_cylins.empty()) {
#ifdef USE_TREE_BB_TEX_ATLAS
		render_branch_texture.nx = 2;
#endif
		render_tree_branches_to_texture_t renderer(TREE_BILLBOARD_SIZE);
		renderer.render_tree(*this, render_branch_texture);
	}
}


void tree_data_t::pre_branch_draw(shader_t &s, bool shadow_only) {
	glEnable(GL_CULL_FACE);
	s.enable_vnct_atribs(1, !shadow_only, !shadow_only, 0); // vertices only in shadow_only mode
}
void tree_data_t::pre_leaf_draw(shader_t &s) {
	s.enable_vnct_atribs(1, 0, 1, 1);
}
void tree_data_t::post_branch_draw(bool shadow_only) {
	glDisable(GL_CULL_FACE);
	if (!shadow_only) {bind_vbo(0, 1);} // branch index vbo
	bind_vbo(0, 0);
}
void tree_data_t::post_leaf_draw() {
	bind_vbo(0, 0);
}


void tree_data_t::draw_leaf_quads_from_vbo(unsigned max_leaves) const {

	check_bind_vbo(leaf_vbo);
	leaf_vert_type_t::set_vbo_arrays(0);
	assert(max_leaves <= leaves.size() && leaf_data.size() >= 4*leaves.size());
	draw_quads_as_tris(4*max_leaves);
}


void tree_data_t::draw_leaves_shadow_only(float size_scale) {

	if (leaves.empty()) return;
	select_texture(tree_types[tree_type].leaf_tex);
	draw_leaves(size_scale); // could disable normals and colors, but that doesn't seem to help much
}


float tree::calc_size_scale(point const &draw_pos) const {

	float const dist_sq(distance_to_camera_sq(draw_pos));

	if (world_mode == WMODE_INF_TERRAIN) {
		float const dmax(get_draw_tile_dist());
		if (dist_sq > dmax*dmax) return 0.0; // too far away to draw
	}
	return (do_zoom ? ZOOM_FACTOR : 1.0f)*tdata().base_radius*InvSqrt(dist_sq)/DIST_C_SCALE;
}

float tree_data_t::get_size_scale_mult() const {return (has_4th_branches ? LEAF_4TH_SCALE : 1.0);}

float tree::pre_transform(vector3d const &tree_xlate) const { // returns rotation angle
	fgPushMatrix();
	translate_to(tree_xlate);
	float rot_angle(0.0);

	// rotate trees in tiled terrain only; ground mode has fewer trees, and they're often all unique anyway;
	// also, ground mode trees often have collisions, dropping leaves, wind, fires, indir lighting, etc. that would be wrong when rotated
	if (enable_rotate_trees()) {
		float const xy_mult(1.0/tdata().sphere_radius); // need enough random variation between adjacent trees
		rot_angle = TWO_PI*fract(xy_mult*tree_xlate.x) + fract(xy_mult*tree_xlate.y); // random angle based on pos
		fgRotateRadians(rot_angle, 0.0, 0.0, 1.0); // rotate around Z axis
	}
	return rot_angle;
}
void tree::post_transform() const {fgPopMatrix();}

void tree::draw_branches_top(shader_t &s, tree_lod_render_t &lod_renderer, bool shadow_only, bool reflection_pass, vector3d const &xlate, int wsoff_loc) {

	if (!created || (!shadow_only && not_visible)) return;
	tree_data_t &td(tdata());
	vector3d const tree_xlate(tree_center + xlate);
	if (!camera_pdu.cube_visible_likely(td.branches_bcube + tree_xlate)) return;
	bool const ground_mode(world_mode == WMODE_GROUND), wind_enabled(ground_mode && (display_mode & 0x0100) != 0);

	if (shadow_only) {
		if (ground_mode && !is_over_mesh()) return;
		pre_transform(tree_xlate);
		if (wsoff_loc >= 0) {s.set_uniform_vector4d(wsoff_loc, vector4d((tree_xlate - get_camera_coord_space_xlate()), 1.0));} // needed for CSMs; no support for rotate; scale=1.0
		td.draw_branches((ground_mode ? (wind_enabled ? last_size_scale : 0.0) : 1.0), 1, 1); // draw branches (untextured), low_detail=1, shadow_pass=1
		post_transform();
		return;
	}
	point const draw_pos(sphere_center() + xlate);
	float const size_scale(calc_size_scale(draw_pos));
	last_size_scale = size_scale;
	if (size_scale < 0.05) return; // if too far away, don't draw any branches
	colorRGBA bcolor(tree_types[type].barkc);
	if (!ground_mode) {bcolor *= 0.8;} // darken slightly in TT mode to account for lack of shadowing on tree branches/trunk
	float const dval(1.0f - 0.95f*damage);
	UNROLL_3X(bcolor[i_] *= min(1.0f, dval*tree_color[i_]);)

	if (lod_renderer.is_enabled()) {
		float const lod_start(tree_lod_scales[0]), lod_end(tree_lod_scales[1]), lod_denom(lod_start - lod_end);
		float geom_opacity(1.0);

		if (td.get_render_branch_texture().is_valid() && size_scale < lod_start) {
			geom_opacity = ((lod_denom == 0.0) ? 0.0 : CLIP_TO_01((size_scale - lod_end)/lod_denom));
			lod_renderer.add_branches(&td, draw_pos, (1.0 - geom_opacity), bcolor);
		}
		if (geom_opacity == 0.0) return;
		s.set_uniform_float(lod_renderer.branch_opacity_loc, geom_opacity);
	}
	select_texture(tree_types[type].bark_tex);
	float const rot_angle(pre_transform(tree_xlate));
	s.set_cur_color(bcolor);
	s.set_uniform_vector4d(wsoff_loc, vector4d((tree_xlate - get_camera_coord_space_xlate()), rot_angle));
	td.draw_branches(size_scale, reflection_pass);
	post_transform();
}


void tree::draw_leaves_top(shader_t &s, tree_lod_render_t &lod_renderer, bool shadow_only, bool reflection_pass, vector3d const &xlate,
	int wsoff_loc, int tex0_off, vector<tree *> &to_update_leaves)
{
	if (!created) return;
	tree_data_t &td(tdata());
	bool const ground_mode(world_mode == WMODE_GROUND), wind_enabled(ground_mode && (display_mode & 0x0100) != 0);
	vector3d const tree_xlate(tree_center + xlate);

	if (shadow_only) {
		if (ground_mode && !is_over_mesh()) return;
		pre_transform(tree_xlate);
		td.leaf_draw_setup(1);
		// Note: since the shadow map is updated every frame when wind is enabled, we can use dynamic LOD without locking in a low-LOD static shadow map
		td.draw_leaves_shadow_only((ground_mode ? (wind_enabled ? last_size_scale : 0.0) : 0.5));
		post_transform();
		return;
	}
	bool const has_leaves(!td.get_leaves().empty());
	
	if (has_leaves && ground_mode && !reflection_pass) {
		burn_leaves();
		if (l_strike.is_enabled() && animate2) {lightning_damage(l_strike.get_hit_pos());}
		if (begin_motion && animate2) {drop_leaves();}
	}
	not_visible = !is_visible_to_camera(xlate); // first pass only
	if (not_visible && !leaf_color_changed) return; // if leaf_color_changed=1, we always draw the leaves as that forces the leaf color update
	if (!has_leaves) return; // only after not_visible is calculated
	if (!leaf_color_changed && !camera_pdu.cube_visible_likely(td.leaves_bcube + tree_xlate)) return;
	point const draw_pos(sphere_center() + xlate);
	float const size_scale(calc_size_scale(draw_pos));
	last_size_scale = size_scale;
	if (size_scale == 0.0) return;
	td.gen_leaf_color();

	if (lod_renderer.is_enabled()) {
		float const lod_start(tree_lod_scales[2]), lod_end(tree_lod_scales[3]), lod_denom(lod_start - lod_end);
		float geom_opacity(1.0);

		if (td.get_render_leaf_texture().is_valid() && size_scale < lod_start) {
			geom_opacity = ((lod_denom == 0.0) ? 0.0 : CLIP_TO_01((size_scale - lod_end)/lod_denom));
			lod_renderer.add_leaves(&td, (draw_pos + vector3d(0.0, 0.0, (td.lr_z_cent - td.sphere_center_zoff))), (1.0 - geom_opacity));
		}
		if (geom_opacity == 0.0) return;
		s.set_uniform_float(lod_renderer.leaf_opacity_loc, geom_opacity);
	}
	bool const leaf_dynamic_en(!has_snow && enable_leaf_wind && wind_enabled && !reflection_pass), gen_arrays(td.leaf_draw_setup(leaf_dynamic_en));
	if (!gen_arrays && leaf_dynamic_en && size_scale*td.get_size_scale_mult() > (leaf_orients_valid ? 0.75f : 0.2f)) {update_leaf_orients_all(to_update_leaves);}

	if (gen_arrays || leaf_color_changed) {
		for (unsigned i = 0; i < leaf_cobjs.size(); ++i) {update_leaf_cobj_color(i);}
	}
	float const rot_angle(pre_transform(tree_xlate));

	if ((has_dl_sources || (tree_indir_lighting && smoke_tid)) && (ground_mode || !get_city_lights_bcube().is_all_zeros())) {
		s.set_uniform_vector4d(wsoff_loc, vector4d((tree_xlate - get_camera_coord_space_xlate()), rot_angle));
	}
	s.set_uniform_int(tex0_off, TLEAF_START_TUID+type); // what about texture color mod?
	td.draw_leaves(size_scale);
	post_transform();
}


template <typename T> void add_cylin_indices_tris(vector<T> &idata, unsigned ndiv, unsigned ix_start, unsigned &idix, unsigned step) {

	for (unsigned S = 0; S < ndiv; S += step) {
		bool const last_edge(S == ndiv-step);
		unsigned const ix0(ix_start + S), ixs[4] = {0, ndiv, (last_edge ? step : step+ndiv), (last_edge ? step-ndiv : step)};
		for (unsigned i = 0; i < 6; ++i) {idata[idix++] = ix0 + ixs[quad_to_tris_ixs[i]];}
	}
}

template<typename branch_index_t> void tree_data_t::create_branch_vbo() {

	vector<branch_vert_type_t> data; // static/reused?
	vector<branch_index_t> idata;
	unsigned const numcylin((unsigned)all_cylins.size());
	unsigned cylin_id(0), data_pos(0), quad_id(0), dix(0), idix(0);
	assert(numcylin > 0);
	branch_index_bytes = sizeof(branch_index_t);
	assert(num_unique_pts < (1ULL << 8*branch_index_bytes));
	data.resize(num_unique_pts);
	idata.resize(9*num_branch_quads); // quads + quads/2 for LOD
	unsigned idix2(6*num_branch_quads);

	for (unsigned i = 0; i < numcylin; i++) {
		draw_cylin const &cylin(all_cylins[i]);
		unsigned const ndiv(cylin.get_num_div());
		point const ce[2] = {cylin.p1, cylin.p2};
		float const ndiv_inv(1.0/ndiv);
		vector3d v12; // (ce[1] - ce[0]).get_norm()
		vector_point_norm const &vpn(gen_cylinder_data(ce, cylin.r1, cylin.r2, ndiv, v12, NULL, 0.0, 1.0, 0)); // force_dim=0
		bool const prev_connect(i > 0 && cylin.can_merge(all_cylins[i-1]));

		if (!prev_connect) { // new cylinder section
			data_pos = dix;
			quad_id  = cylin_id = 0;
		}
		for (unsigned j = prev_connect; j < 2; ++j) { // create vertex data
			float const ty(b_tex_scale*(cylin_id + j)); // Note: not in [0.0, 1.0] range

			for (unsigned S = 0; S < ndiv; ++S) { // first cylin: 0,1 ; other cylins: 1
				float const tx(2.0f*fabs(S*ndiv_inv - 0.5f));
				vector3d const n(0.5f*(vpn.n[S] + vpn.n[(S+ndiv-1)%ndiv])); // average face normals to get vert normals, don't need to normalize
				data[dix++] = branch_vert_type_t(vpn.p[(S<<1)+j], n, tx, ty);
			}
		}
		add_cylin_indices_tris(idata, ndiv, (data_pos + quad_id), idix,  1); // create index data
		add_cylin_indices_tris(idata, ndiv, (data_pos + quad_id), idix2, 2); // create index data for next LOD level
		quad_id += ndiv;
		++cylin_id;
	} // for i
	assert(idix  == 6*num_branch_quads);
	assert(dix   == data.size());
	assert(idix2 == idata.size());
	branch_manager.create_and_upload(data, idata); // ~350KB data + ~75KB idata (with 16-bit index)
}

void tree_data_t::ensure_branch_vbo() {

	if (branch_manager.vbo != 0) return; // vbos already created
	assert(branch_manager.ivbo == 0);
	assert(num_branch_quads == 0 && num_unique_pts == 0);
	num_branch_quads = 0;

	for (unsigned i = 0; i < (unsigned)all_cylins.size(); ++i) { // determine required data size
		bool const prev_connect(i > 0 && all_cylins[i].can_merge(all_cylins[i-1]));
		unsigned const ndiv(all_cylins[i].get_num_div());
		num_branch_quads += ndiv;
		num_unique_pts   += (prev_connect ? 1 : 2)*ndiv;
	}
	// determine 16 vs. 32 bits for branch index buffer based on number of branch points
	if (num_unique_pts >= (1 << 16)) {create_branch_vbo<unsigned>();} else {create_branch_vbo<unsigned short>();}
}


void tree_data_t::draw_branches(float size_scale, bool force_low_detail, bool shadow_pass) {

	unsigned const num((size_scale == 0.0) ? num_branch_quads : min(num_branch_quads, max((num_branch_quads/40), unsigned(1.5*num_branch_quads*size_scale*get_size_scale_mult())))); // branch LOD
	bool low_detail(force_low_detail || (size_scale > 0.0 && size_scale < 2.0));
	// need high detail when using 4th order branches, since otherwise they would only have ndiv=2 (zero width); but this is okay for the shadow pass
	if (!shadow_pass && has_4th_branches && num > num_branch_quads/5) {low_detail = 0;}
	ensure_branch_vbo();
	branch_manager.pre_render();
	vert_norm_comp_tc::set_vbo_arrays(0);
	int const index_type((branch_index_bytes == 2) ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT);
	unsigned const idata_sz(6*num_branch_quads*branch_index_bytes);
	glDrawRangeElements(GL_TRIANGLES, 0, num_unique_pts, (low_detail ? 3 : 6)*num, index_type, (void *)size_t(low_detail ? idata_sz : 0));
	++num_frame_draw_calls;
	branch_manager.post_render();
}


void tree_data_t::update_normal_for_leaf(unsigned i) {

	assert(i < leaves.size());
	vector3d const &normal(leaves[i].norm); // standard leaf plane normal
	//vector3d const &normal(leaves[i].get_center().get_norm()); // use position of leaf relative to tree center (could also mix in leaf normal)
	norm_comp nc; nc.set_norm_no_clamp(normal); // already normalized, no need to clamp
	UNROLL_4X(leaf_data[i_+(i<<2)].set_norm(nc);)
	mark_leaf_changed(i);
}


void tree_data_t::reset_leaf_pos_norm() {

	for (unsigned i = 0; i < (unsigned)leaves.size(); i++) { // process leaf points - reset to default positions and normals
		UNROLL_4X(leaf_data[i_+(i<<2)].v = leaves[i].pts[i_];)
		update_normal_for_leaf(i);
	}
	reset_leaves = 0;
}


void tree_data_t::ensure_leaf_vbo() {

	leaf_change_end = min(leaf_change_end, (unsigned)leaves.size()); // in case a leaf was removed after a previous update this frame

	if (leaf_vbo == 0) {
		create_vbo_and_upload(leaf_vbo, leaf_data, 0, 0, 1); // dynamic draw, due to wind updates, collision, burn damage, etc.
	}
	else if (leaf_change_start < leaf_change_end) {
		assert(leaf_change_end <= leaves.size());
		bind_vbo(leaf_vbo);
		unsigned const per_leaf_stride(4*sizeof(leaf_vert_type_t));
		upload_vbo_sub_data((leaf_data.data() + 4*leaf_change_start), leaf_change_start*per_leaf_stride, (leaf_change_end - leaf_change_start)*per_leaf_stride);
	}
	leaf_change_start = leaves.size();
	leaf_change_end   = 0;
}


void tree_data_t::draw_leaves(float size_scale) {

	ensure_leaf_vbo();
	unsigned nl(leaves.size());
	if (ENABLE_CLIP_LEAVES && size_scale > 0.0) {nl = min(nl, max((nl/8), unsigned(4.0*nl*size_scale*get_size_scale_mult())));} // leaf LOD
	draw_leaf_quads_from_vbo(nl);
}


bool tree_data_t::leaf_draw_setup(bool no_leaf_reset) {

	bool const gen_arrays(!leaf_data_allocated());
	if (gen_arrays) {alloc_leaf_data();}
	if (gen_arrays || (reset_leaves && !no_leaf_reset)) {reset_leaf_pos_norm();}
	// may do duplicate work, but should at least keep the cobj colors in sync
	if (gen_arrays || leaf_color_changed) {update_all_leaf_colors();}
	return gen_arrays;
}


void tree_data_t::bend_leaf(unsigned i, float angle) { // Note: slow

	assert(i < leaves.size());
	tree_leaf const &l(leaves[i]);
	vector3d const orig_dir(l.pts[1] - l.pts[0]); // vector from base to tip
	vector3d const new_dir(orig_dir*COSF(angle) + l.norm*(orig_dir.mag()*SINF(angle))); // s=orig_dir.get_norm(), t=l.norm
	vector3d const delta(new_dir - orig_dir);
	vector3d const normal(cross_product(new_dir, (l.pts[3] - l.pts[0])).get_norm());
	unsigned const ix(i<<2);
	leaf_data[ix+1].v = l.pts[1] + delta;
	leaf_data[ix+2].v = l.pts[2] + delta;
	norm_comp nc; nc.set_norm_no_clamp(normal); // already normalized, no need to clamp
	UNROLL_4X(leaf_data[i_+ix].set_norm(nc);) // similar to update_normal_for_leaf()
	mark_leaf_changed(i);
	reset_leaves = 1; // do we want to update the normals as well?
}


bool tree_data_t::check_if_needs_updated() {

	bool const do_update(last_update_frame < frame_counter);
	last_update_frame = frame_counter;
	return do_update;
}


void tree::update_leaf_orients_wind() { // leaves move in wind (somewhat slow)

	tree_data_t &td(tdata());
	vector<tree_leaf> const &leaves(td.get_leaves());
	rand_gen_t rgen;
	rgen.set_state(frame_counter, leaves.size());
	bool const priv_data(td_is_private());
	bool const heal_pass(priv_data && LEAF_HEAL_RATE > 0 && world_mode == WMODE_GROUND && (rgen.rand()&7) == 0); // only update healed color every 8 frames
	int last_xpos(0), last_ypos(0);
	vector3d local_wind(zero_vector);

	for (unsigned i = 0; i < leaves.size(); ++i) { // process leaf wind and collisions
		point p0(leaves[i].pts[0]);
		if (priv_data) {p0 += tree_center;}
		int const xpos(get_xpos(p0.x)), ypos(get_ypos(p0.y));
			
		// Note: should check for similar z-value, but z is usually similar within the leaves of a single tree
		if (i == 0 || xpos != last_xpos || ypos != last_ypos) {
			local_wind = get_local_wind(xpos, ypos, p0.z, !priv_data); // slow
			last_xpos  = xpos;
			last_ypos  = ypos;
		}
		if (local_wind != zero_vector) {
			float const angle(PI_TWO*max(-1.0f, min(1.0f, dot_product(local_wind, leaves[i].norm)))); // not physically correct, but it looks good
			td.bend_leaf(i, angle);
		}
		if (heal_pass && (rgen.rand()&63) == 0) { // leaf heals every 64 frames
			short &lcolor(td.get_leaves()[i].lcolor); // non-const, can't use <leaves>

			if (lcolor > 0 && lcolor < 1000) { // partially damaged
				lcolor = min(1000, (lcolor + int(LEAF_HEAL_RATE*fticks)));
				copy_color(i);
			}
		}
	} // for i
	leaf_orients_valid = 1;
}

void tree::update_leaf_orients_all(vector<tree *> &to_update_leaves) {

	tree_data_t &td(tdata());
	if (td.check_if_needs_updated() || !leaf_orients_valid) {to_update_leaves.push_back(this);}
	if (!has_any_billboard_coll || leaf_cobjs.empty()) return; // no coll movement
	// leaves move when struck by an object
	vector<tree_leaf> const &leaves(td.get_leaves());
	rand_gen_t rgen;
	rgen.set_state(frame_counter, leaves.size());
	unsigned num(min(leaves.size(), leaf_cobjs.size())); // should be the same size

	for (unsigned i = 0; i < num; ++i) {
		coll_obj &cobj(get_leaf_cobj(i));
		if (cobj.last_coll == 0) continue; // no collision
		bool removed(0);
		float hit_angle(0.0);
		if (cobj.coll_type == BEAM) {removed = damage_leaf(i, BEAM_DAMAGE, rgen);} // do burn damage
		else {
			if (cobj.coll_type == PROJECTILE) {
				if ((rgen.rand()&3) == 0) { // shoot off leaf
					create_leaf_obj(i);
					remove_leaf(i, 1);
					removed = 1;
				}
				cobj.coll_type = IMPACT; // reset to impact after first hit
			}
			hit_angle += PI_TWO*cobj.last_coll/TICKS_PER_SECOND; // 90 degree max rotate
		}
		if (removed) {--i; --num; continue;}
		cobj.last_coll = ((hit_angle == 0.0 || cobj.last_coll < iticks) ? 0 : (cobj.last_coll - iticks));
		if (cobj.last_coll) {next_has_any_billboard_coll = 1;}
		if (hit_angle != 0.0) {td.bend_leaf(i, hit_angle);} // do we want to update collision objects as well?
	} // for i
}


unsigned tree_cont_t::delete_all() {

	unsigned deleted(0);

	for (reverse_iterator i = rbegin(); i != rend(); ++i) { // delete backwards (pop collision stack)
		if (i->delete_tree()) {++deleted;}
	}
	generated = 0;
	return deleted;
}

void delete_trees() {
	t_trees.delete_all();
	if (tree_coll_level && !t_trees.empty()) {purge_coll_freed(1);} // MUST do this for collision detection to work
}

void tree::clear_context() {
	if (td_is_private()) {tdata().clear_context();}
}

int tree::delete_tree() {

	clear_context();
	if (!created)  return 0;
	if (tree_coll_level) {remove_collision_objects();}
	if (no_delete) return 0;
	if (td_is_private()) {tdata().clear_data();}
	tree_fire.reset();
	created = 0;
	return 1;
}

void tree_data_t::clear_data() {
	
	clear_context();
	clear_cont(all_cylins);
	clear_cont(leaf_data);
	clear_cont(leaves); // Note: not present in original delete_trees()
}


void copy_cylins(tree_cylin *start_cylin, int num, tree_cylin *&cur_cylin) {
	for (int i = 0; i < num; ++i) {*cur_cylin = *(start_cylin+i); ++cur_cylin;}
}


void tree_builder_t::create_all_cylins_and_leaves(vector<draw_cylin> &all_cylins, vector<tree_leaf> &leaves,
	int tree_type, float deadness, float br_scale, float nl_scale, bool has_4th_branches, int tree_size)
{
	// compact the cylinders into a contiguous block
	assert(all_cylins.empty());
	tree_cylin *cylins(&cylin_cache.front()), *cur_cylin(cylins);
	assert(base.cylin == cur_cylin);
	copy_cylins(base.cylin,  base_num_cylins, cur_cylin);
	copy_cylins(roots.cylin, root_num_cylins, cur_cylin);

	for (int i = 0; i < num_1_branches; i++) { // add the first order branches
		copy_cylins(branches[i][0].cylin, branches[i][0].num_cylins, cur_cylin);
	}
	for (int i = 0; i < num_1_branches; i++) { // add second order branches
		for (int k = 1; k <= branches[i][0].num_branches; k++) {
			copy_cylins(branches[i][k].cylin, branches[i][k].num_cylins, cur_cylin);
		}
	}
	for (unsigned w = 0; w < 2; ++w) {
		for (int i = 0; i < num_34_branches[w]; i++) { // add the third and fourth order branches
			copy_cylins(branches_34[w][i].cylin, branches_34[w][i].num_cylins, cur_cylin);
		}
	}

	// create the all_cylins vector
	unsigned const num_total_cylins(cur_cylin - cylins);
	all_cylins.reserve(num_total_cylins);

	for (unsigned i = 0; i < num_total_cylins; ++i) {
		tree_cylin &c(cylins[i]);
		assert(c.r1 > 0.0 || c.r2 > 0.0);
		assert(c.p1 != c.p2);
		c.r1 *= br_scale;
		c.r2 *= br_scale;
		
		if (clip_cube && c.level > 0) { // clip non-trunk non-roots
			if (!clip_cube->contains_pt(c.p1)) {c.p2 = c.p1; continue;} // stars outside bounds - make invalid and skip it
			if (!clip_cube->contains_pt(c.p2)) {c.r2 = 0.0;} // ends outside bounds - close off the end with a point
			clip_cube->clamp_pt(c.p1);
			clip_cube->clamp_pt(c.p2);
		}
		all_cylins.push_back(c);
	}

	// add leaves
	if (deadness < 1.0) {
		float const leaf_size(REL_LEAF_SIZE*TREE_SIZE*(has_4th_branches ? LEAF_4TH_SCALE : 1.0f)/(sqrt(nl_scale*nleaves_scale)*tree_scale));
		float rel_leaf_size(2.0*leaf_size*tree_types[tree_type].leaf_size*(tree_scale*base_radius/TREE_SIZE + 10.0)/18.0);
		if (tree_size > 100) {rel_leaf_size *= sqrt(tree_size/100.0);} // scale leaf size to match branch size for large custom trees
		unsigned nl(unsigned((1.0 - deadness)*num_leaves_per_occ*num_total_cylins) + 1); // determine the number of leaves
		//cout << "branch cylinders: " << num_total_cylins << ", leaves: " << nl << endl;
		leaves.reserve(nl);

		for (unsigned i = 0; i < num_total_cylins; ++i) {
			if (cylins[i].p2 != cylins[i].p1 && cylins[i].level > 1 && (TREE_4TH_LEAVES || cylins[i].level < 4)) { // leaves will still be allocated
				add_leaves_to_cylin(i, tree_type, rel_leaf_size, deadness, leaves);
			}
		}
		// randomly reorder leaves so that LOD produces more even partial coverage
		rand_gen_t rgen2;
		unsigned const nleaves(leaves.size());
		for (unsigned i = 0; i < nleaves; ++i) {swap(leaves[i], leaves[rgen2.rand()%nleaves]);}
		//remove_excess_cap(leaves);
	}

	// now NULL (pseudo free) the individual branches because they are unneccessary
	branches_34[0] = branches_34[1] = NULL;
	branches = NULL;
	roots.cylin = base.cylin = NULL;
	assert(clip_cube || all_cylins.size() == all_cylins.capacity());
}


void tree_xform_t::rotate_cylin(tree_cylin &c) {
	rotate_around_axis(c);
	c.p2 = c.p1 + re_matrix;
}

void tree_xform_t::set_sin_cos_terms(float deg_rotate) {

	if (deg_rotate == last_deg_rotate) return; // already set
	float const angle(deg_rotate/TO_DEG);
	cos_term = cos(angle);
	sin_term = sin(angle);
	last_deg_rotate = deg_rotate;
}

void tree_xform_t::rotate_pts_around_axis(point const &p, point const &rotate, float deg_rotate) {

	set_sin_cos_terms(deg_rotate);
	// rotate around z
	re_matrix.x = -rotate.x*p.x + rotate.y*p.y;
	re_matrix.y = -rotate.y*p.x - rotate.x*p.y;
	// rotate around x
	float const xv(re_matrix.x), yv(re_matrix.y);
	re_matrix.y = cos_term*yv - sin_term*p.z;
	re_matrix.z = sin_term*yv + cos_term*p.z;
	// do_reverse, inverse rotate around z
	float const yv2(re_matrix.y);
	re_matrix.x = -rotate.x*xv - rotate.y*yv2;
	re_matrix.y =  rotate.y*xv - rotate.x*yv2;
}

void tree_xform_t::rotate_around_axis(tree_cylin const &c) {

	set_sin_cos_terms(c.deg_rotate);
	float const yv(-sin_term*c.length);
	re_matrix.assign(-c.rotate.y*yv, c.rotate.x*yv, cos_term*c.length);
}

void tree_xform_t::gen_cylin_rotate(vector3d &rotate, vector3d &lrotate, float rotate_start) {

	float temp_deg(safe_acosf(lrotate.x));
	if (lrotate.y < 0.0) {temp_deg *= -1.0;}
	setup_rotate(rotate, rotate_start, temp_deg);
}


float get_default_tree_depth() {return TREE_DEPTH*(0.5 + 0.5/tree_scale);}

float get_tree_size_scale(int tree_type, bool create_bush) {
	return TREE_SIZE*tree_types[tree_type].branch_size/tree_scale * (create_bush ? 0.7 : 1.0);
}

void get_tile_z_minmax_for_area(tile_t const &tile, point const &pos, float radius, float &zmin, float &zmax);

bool adjust_tree_zval(point &pos, int size, int type, bool create_bush, tile_t const *const cur_tile=nullptr) {

	int const size_est((size == 0) ? 60 : size); // use average of 40-80
	float const size_scale(get_tree_size_scale(type, create_bush)), base_radius(size_est*(0.1*size_scale)), radius(2.0*base_radius);
	float mzmax(pos.z);

	if (cur_tile) { // use tile to more quickly evaluate (estimate) min/max height
		get_tile_z_minmax_for_area(*cur_tile, (pos + vector3d(xoff2*DX_VAL, yoff2*DY_VAL, 0.0)), 0.5*radius, pos.z, mzmax); // offset by current mesh transform
		return ((mzmax - pos.z) < tree_slope_thresh*radius); // drop trees on steep slopes
	}
	for (int dy = -1; dy <= 1; dy += 2) { // take max in 4 directions to prevent intersections with the terrain on steep slopes
		for (int dx = -1; dx <= 1; dx += 2) {
			float const z(interpolate_mesh_zval((pos.x + dx*radius), (pos.y + dy*radius), 0.0, 1, 1));
			min_eq(pos.z, z); // use min zval to ensure roots are below the soil
			max_eq(mzmax, z);
			if ((mzmax - pos.z) > tree_slope_thresh*radius) return 0; // drop trees on steep slopes
		}
	}
	return 1;
}


//gen_tree(pos, size, ttype>=0, calc_z, 0, 1);
//gen_tree(pos, 0, -1, 1, 1, 0);
void tree::gen_tree(point const &pos, int size, int ttype, int calc_z, bool add_cobjs, bool user_placed, rand_gen_t &rgen,
	float height_scale, float br_scale_mult, float nl_scale, bool has_4th_branches, bool allow_bushes, bool force_bushes)
{
	//assert(calc_z || user_placed);
	tree_center      = pos;
	created          = 1;
	tree_color.alpha = 1.0;
	tree_nl_scale    = nl_scale;
	vector3d const color_var(rgen.signed_rand_vector()); // rand_gen() called outside gen_tree_data()
	tree_data_t &td(tdata());

	if (td.is_created()) { // pre-allocated, shared tree
		assert(!td_is_private());
		assert(!user_placed);
		//assert(size <= 0); // too strong, but note that size is only used when creating new trees
		int const td_type(tdata().get_tree_type());
		
		if (ttype < 0) {
			type = (FORCE_TREE_TYPE ? td_type : (rgen.rand() % NUM_TREE_TYPES));
		}
		else if (FORCE_TREE_TYPE) {
			type = ttype;
			// Note: we could relax this restriction and let type be different, which gives us more tree variation,
			// but then the leaf and branch sizes/colors couldn't be different per tree
			assert(type == td_type); // up to the caller to ensure this
		}
		UNROLL_3X(tree_color[i_] = 1.0f + tree_types[type].branch_color_var*color_var[i_];)
	}
	else {
		type = ((ttype < 0) ? rgen.rand() : ttype) % NUM_TREE_TYPES; // maybe should be an error if > NUM_TREE_TYPES
		bool const create_bush(force_bushes ? 1 : (allow_bushes && rgen.rand_probability(tree_types[type].bush_prob)));
		if (create_bush) {type = (type + 1) % NUM_TREE_TYPES;} // mix up the tree types so that bushes stand out from trees
		tree_type const &treetype(tree_types[type]);
		UNROLL_3X(tree_color[i_] = 1.0f + treetype.branch_color_var*color_var[i_];)
		float tree_depth(get_default_tree_depth());
	
		if (calc_z) {
			tree_center.z = interpolate_mesh_zval(tree_center.x, tree_center.y, 0.0, 1, 1);
			adjust_tree_zval(tree_center, size, type, create_bush); // ignore return value, as it's too late to drop the tree due to slope here
			if (user_placed) {tree_depth = tree_center.z - get_tree_z_bottom(tree_center.z, tree_center);} // more accurate
		}
		cube_t const cc(clip_cube - tree_center); // Note: cc is only used if use_clip_cube==1
		float const hscale((height_scale == 1.0) ? treetype.height_scale : height_scale);
		float const br_scale((br_scale_mult == 1.0) ? treetype.branch_radius : br_scale_mult);
		float const bbo_scale((height_scale == 1.0) ? treetype.branch_break_off : 1.0);
		td.gen_tree_data(type, size, tree_depth, hscale, br_scale, nl_scale, bbo_scale, has_4th_branches, (use_clip_cube ? &cc : NULL), create_bush, rgen); // create the tree here
	}
	assert(type < NUM_TREE_TYPES);
	unsigned const nleaves(td.get_leaves().size());
	damage_scale = (nleaves ? 1.0/nleaves : 0.0);
	damage       = 0.0;
	if (add_cobjs) {add_tree_collision_objects();}
}


void tree_data_t::gen_tree_data(int tree_type_, int size, float tree_depth, float height_scale, float br_scale_mult,
	float nl_scale, float bbo_scale, bool has_4th_branches_, cube_t const *clip_cube, bool create_bush, rand_gen_t &rgen)
{
	//RESET_TIME;
	tree_type = tree_type_;
	has_4th_branches = has_4th_branches_;
	assert(tree_type < NUM_TREE_TYPES);
	leaf_data.clear();
	clear_vbo_ixs();
	float deadness(DISABLE_LEAVES ? 1.0 : tree_deadness);

	if (deadness < 0.0) {
		int const num(rgen.rand_int(1, 100));
		deadness = ((num > 94) ? min(1.0f, float(num - 94)/8.0f) : 0.0);
	}
	if (deadness < 1.0 && rgen.rand_probability(tree_dead_prob)) {deadness = 1.0;}
	tree_builder_t builder(clip_cube, rgen);
	
	// create leaves and all_cylins
	br_scale    = br_scale_mult*branch_radius_scale;
	b_tex_scale = tree_types[tree_type].branch_tscale*height_scale/br_scale;
	base_radius = builder.create_tree_branches(tree_type, size, tree_depth, base_color, height_scale, br_scale, nl_scale, bbo_scale, has_4th_branches, create_bush);
	builder.create_all_cylins_and_leaves(all_cylins, leaves, tree_type, deadness, br_scale, nl_scale, has_4th_branches, size);

	// set the bounding sphere center
	assert(!all_cylins.empty());
	branches_bcube.set_from_point(all_cylins.front().p1);

	for (vector<draw_cylin>::const_iterator i = all_cylins.begin(); i != all_cylins.end(); ++i) {
		branches_bcube.union_with_pt(i->p2); // only need to check branch (and root) end points
	}
	br_x = max(branches_bcube.d[0][1], -branches_bcube.d[0][0]);
	br_y = max(branches_bcube.d[1][1], -branches_bcube.d[1][0]);
	float const bzmin(branches_bcube.d[2][0]), bzmax(branches_bcube.d[2][1]);
	sphere_center_zoff = 0.5f*(bzmin + bzmax);
	br_z               = 0.5f*(bzmax - bzmin);
	sphere_radius      = 0.0f;
	point const center(get_center());
	leaves_bcube.set_from_point(leaves.empty() ? center : leaves.front().pts[0]);

	for (vector<draw_cylin>::const_iterator i = all_cylins.begin(); i != all_cylins.end(); ++i) {
		sphere_radius = max(sphere_radius, p2p_dist_sq(i->p2, center)); // only need to check branch (and root) end points
	}
	for (vector<tree_leaf>::const_iterator i = leaves.begin(); i != leaves.end(); ++i) {
		UNROLL_4X(leaves_bcube.union_with_pt(i->pts[i_]);)
	}
	lr_x  = max(leaves_bcube.d[0][1], -leaves_bcube.d[0][0]);
	lr_y  = max(leaves_bcube.d[1][1], -leaves_bcube.d[1][0]);
	float const lr_z1(leaves_bcube.d[2][0]), lr_z2(leaves_bcube.d[2][1]);
	sphere_radius = sqrt(sphere_radius);
	lr_z_cent     = 0.5f*(lr_z1 + lr_z2);
	lr_z          = 0.5f*(lr_z2 - lr_z1);
	reverse(leaves.begin(), leaves.end()); // order leaves so that LOD removes from the center first, which is less noticeable
	//PRINT_TIME("Gen Tree");
}


float tree_builder_t::create_tree_branches(int tree_type, int size, float tree_depth, colorRGBA &base_color, float height_scale,
	float br_scale, float nl_scale, float bbo_scale, bool has_4th_branches, bool create_bush)
{
	// fixed tree variables
	// bush: minimal trunk, no roots, starts below the ground, smaller, fewer leaves
	ncib                 = 10;
	base_num_cylins      = (create_bush ? 1 : round_fp(5 * max(1.0f, tree_height_scale)));
	num_cylin_factor     = 10.0;
	base_cylin_factor    = 1.0*base_num_cylins; //(type == TREE_A) ? 0.4 : 1.0;
	base_break_off       = max(1, min(base_num_cylins-1, round_fp(bbo_scale*(2 + base_num_cylins/5.0)))); // controls the point at which branches start on the base
	num_1_branches       = 8; // max value
	num_big_branches_min = 3;
	num_big_branches_max = 4;
	num_2_branches_min   = 4;
	num_2_branches_max   = 6;
	num_34_branches[0]   = 0; // this should start out 0
	num_3_branches_min   = 6;
	num_3_branches_max   = 10;
	num_34_branches[1]   = (has_4th_branches ? 2000 : 0);
	if (size <= 0) {size = rgen.rand_int(40, 80);} // tree size
	float const size_scale(get_tree_size_scale(tree_type, create_bush));
	base_radius            = size * (0.1*size_scale);
	num_leaves_per_occ     = 0.01f*nl_scale*nleaves_scale*(rgen.rand_int(30, 60) + size) * (create_bush ? 0.3f : 1.0f);
	base_length_min        = rgen.rand_int(4, 6) * height_scale * base_radius * tree_height_scale * (create_bush ? 0.05f : 1.0f); // short trunk for bushes
	base_length_max        = base_length_min * 1.5;
	angle_rotate           = 60.0;
	base_curveness         = 10.0;
	tree_slimness          = 0;
	tree_wideness          = 70;
	branch_curveness       = 90.0;
	branch_1_var           = 100.0*0.85;
	branch_1_rad_var       = 100.0*0.64; // Note: to be physically correct, branching should preserve cross-sectional area
	branch_1_start         = 0.45;
	branch_2_start         = 0.9;
	branch_2_var           = 100.0*0.85;
	branch_2_rad_var       = 100.0*0.64;
	branch_4_max_radius    = 0.0002;
	rotate_factor          = 1.0; // how much to change the horizontal rotation 
	branch_distribution    = 1.0;
	branch_1_distribution  = 1.0;
	branch_upwardness      = 0.9;
	branch_min_angle       = 20.0; // minimim angle of deflection from the current branch creating it
	branch_max_angle       = 40.0;
	max_2_angle_rotate     = 50.0;
	max_3_angle_rotate     = 50.0;
	leaf_acc               = 0.0;
	float const base_rad_var = (1.0 - 0.75/base_num_cylins);
	float const base_len_var = (1.0 - 1.0 /base_num_cylins);
	float const height_offset(create_bush ? -2.0*height_scale*base_radius*tree_height_scale : 0);
	float const branch_1_random_rotate = 40.0;
	int const min_num_roots(10), max_num_roots(12);
	base_color = colorRGBA(0.5*rgen.signed_rand_float(), 0.5*rgen.signed_rand_float(), 0.0, 1.0); // no blue
	// temporary variables
	int const nbr(num_1_branches*(num_2_branches_max+1)), nbranches(nbr*num_3_branches_max);
	
	// allocate all the memory required for the tree ------------------------------------------------
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
	tree_cylin *cur_ptr(roots.cylin + CYLINS_PER_ROOT*max_num_roots);

	for (int i = 0; i < num_1_branches; i++) {
		branches[i] = branches[0] + i*(num_2_branches_max+1);
		branches[i][0].cylin = cur_ptr;
		cur_ptr += ncib;
	}
	for (int i = 0; i < num_1_branches; i++) {
		for (int j = 1; j < (num_2_branches_max+1); j++) {
			branches[i][j].cylin = cur_ptr;
			cur_ptr += ncib;
		}
	}
	for (int i = 0; i < nbranches; i++) {
		branches_34[0][i].cylin = cur_ptr;
		cur_ptr += ncib;
	}

	// create tree base -------------------------------------------------------------
	int num_b_so_far(0);

	for (int i = 0; i < base_num_cylins; i++) {
		tree_cylin &cylin(base.cylin[i]);
		float const deg_rot(sinf(-PI_TWO + i*PI/base_num_cylins + ((tree_height_scale > 1.0) ? 100.0*base_color.R : 0))*base_curveness); // "hash" the color to get a unique curve

		if (i == 0) {
			float const length((base_length_max - base_length_min)/2.0f);
			cylin.assign_params(0, 0, base_radius, base_radius*base_rad_var, length*(num_cylin_factor/base_num_cylins), deg_rot);
			cylin.p1 = point(0.0, 0.0, height_offset);
			cylin.rotate.assign(cosf(angle_rotate/TO_DEG), sinf(angle_rotate/TO_DEG), 0.0);
		}
		else {
			tree_cylin &lcylin(base.cylin[i-1]);
			cylin.assign_params(0, 0, lcylin.r2, lcylin.r2*base_rad_var, lcylin.length*base_len_var*(base_cylin_factor/base_num_cylins), deg_rot);
			cylin.rotate.assign(lcylin.rotate.x, lcylin.rotate.y, 0.0);
			rotate_cylin(lcylin);
			cylin.p1 = lcylin.p1 + BASE_LEN_SCALE*re_matrix;
		}
		if (i == (base_num_cylins-1)) {rotate_cylin(cylin);} // last cylin

		if (base_break_off <= (i+1) && num_b_so_far < num_1_branches) {
			int const init_num_cyl(base_num_cylins - base_break_off + 1);
			float const temp_num(branch_distribution*float(num_1_branches)/init_num_cyl*float(i+2-base_break_off)/float(num_b_so_far+1));
			float rotate_start(rgen.rand_int(0, 259));
			int num_b_to_create(0);

			if (temp_num >= 1.0) {
				float const temp2((float(num_1_branches) - num_b_so_far)/float(base_num_cylins - i));
				num_b_to_create = min(num_1_branches-num_b_so_far, ((temp2 <= 3.0) ? max(1, int(ceil(temp2))) : int(temp_num + 0.5)));
			}
			for (int j = 0; j < num_b_to_create; j++) {
				assert(num_b_so_far < num_1_branches);
				create_1_order_branch(i, rotate_start, num_b_so_far++);
				rotate_start += 360.0f/float(num_b_to_create) + (3 - 2*rgen.rand_int(1,2))*rgen.rand_int(0,int(branch_1_random_rotate));
			}
		}
	} // for i
	num_1_branches = num_b_so_far; // clamp to the number actually created
	
	// done with base ------------------------------------------------------------------
	if (has_4th_branches) {create_4th_order_branches(nbranches, size_scale*br_scale);}
	root_num_cylins = ((gen_tree_roots && !create_bush) ? CYLINS_PER_ROOT*rgen.rand_int(min_num_roots, max_num_roots) : 0);

	for (int i = 0; i < root_num_cylins; i += CYLINS_PER_ROOT) { // add roots
		tree_cylin &cylin1(roots.cylin[i]), &cylin2(roots.cylin[i+1]), &cylin3(roots.cylin[i+2]);
		float const root_radius(rgen.rand_uniform(0.38, 0.45)*base_radius);
		float const theta((TWO_PI*(i + 0.3*rgen.signed_rand_float()))/root_num_cylins);
		float const deg_rot(180.0+rgen.rand_uniform(40.0, 50.0));
		vector3d const dir(sin(theta), cos(theta), 0.0);
		cylin1.assign_params(1, i, root_radius, 0.75*root_radius, 1.0*base_radius, deg_rot); // level 1, with unique branch_id's
		cylin1.p1     = cylin1.p2 = point(0.0, 0.0, 0.75*br_scale*base_radius);
		cylin1.rotate = dir;
		rotate_cylin(cylin1);
		cylin1.p1 += (0.3*br_scale*base_radius/cylin1.length)*(cylin1.p2 - cylin1.p1); // move away from the tree centerline
		cylin1.p1 *= 1.6; // 2.0 is okay for non-curved trunks
		cylin1.p2 *= 1.3;

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
	if (tree_depth > 0.0 && root_num_cylins == 0 && !create_bush) { // add the bottom cylinder section from the base into the ground
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
	return (rgen.rand_int((int)branch_var - 5, (int)branch_var + 5)/100.0)*(num_cylin_factor/ncib);
}

inline float tree_builder_t::gen_bc_size2(float branch_var) {
	return rgen.rand_int((int)branch_var - 5, (int)branch_var + 5)/100.0;
}


void tree_builder_t::gen_next_cylin(tree_cylin &cylin, tree_cylin &lcylin, float var, float rad_var, int level, int branch_id, bool rad_var_test) {

	cylin.assign_params(level, branch_id, lcylin.r2, (rad_var_test ? lcylin.r1*gen_bc_size(rad_var) : 0.0f),
		((level == 4) ? branch_4_length*var : lcylin.length*gen_bc_size(var)), lcylin.deg_rotate);
	rotate_around_axis(lcylin);
	cylin.p1 = lcylin.p1 + BRANCH_LEN_SCALE*re_matrix;
}


void tree_builder_t::gen_first_cylin(tree_cylin &cylin, tree_cylin &src_cylin, float bstart, float rad_var, float rotate_start, int level, int branch_id) {

	rotate_around_axis(src_cylin);
	cylin.p1 = src_cylin.p1 + BRANCH_LEN_SCALE*re_matrix;
	float const radius1(bstart*src_cylin.r2);
	float deg_rotate(rgen.rand_int(0, int(branch_max_angle)));
	
	if (rotate_start < branch_min_angle && deg_rotate < branch_min_angle) {
		deg_rotate = rgen.rand_int(int(branch_min_angle - rotate_start), int(branch_max_angle));
	}
	if (rgen.rand() & 1) deg_rotate *= -1.0;
	if (level == 2 && deg_rotate < 0.0) deg_rotate *= 0.5;
	deg_rotate += src_cylin.deg_rotate;
	cylin.assign_params(level, branch_id, radius1, radius1*gen_bc_size2(rad_var), bstart*src_cylin.length*(num_cylin_factor/ncib), deg_rotate);
}


void tree_builder_t::create_1_order_branch(int base_cylin_num, float rotate_start, int branch_num) {

	bool branch_just_created(false), branch_deflected(false);
	int rotation((rgen.rand() & 1) ? -1 : 1);

	//first cylinders of the first order branches
	tree_branch &branch(branches[branch_num][0]);
	tree_cylin &cylin(branch.cylin[0]);
	branch.num_cylins = ncib;
	rotate_around_axis(base.cylin[base_cylin_num]);
	cylin.p1 = base.cylin[base_cylin_num].p1 + BRANCH_LEN_SCALE*re_matrix;
	setup_rotate(cylin.rotate, rotate_start, 0.0);
	float const radius1(base_radius*branch_1_start);
	cylin.assign_params(1, branch_num, radius1, radius1*gen_bc_size(branch_1_rad_var),
		radius1*6*((float)rgen.rand_int(7,10)/10)*(num_cylin_factor/ncib),
		(tree_slimness + (tree_wideness-tree_slimness)*rgen.rand_int(6,10)/10.0f*(num_1_branches - branch_num)/num_1_branches));

	//generate stats for the second order branches
	branch.num_branches = rgen.rand_int(num_2_branches_min, num_2_branches_max);

	//create rest of the cylinders
	int num_2_branches_created(0);
	int const temp_num_big_branches(rgen.rand_int(num_big_branches_min, num_big_branches_max));
	float const temp_num(2.0*temp_num_big_branches/(0.5*branch.num_cylins+1));

	if (temp_num*branch_1_distribution >= 1.0f) {
		if (num_2_branches_created < branch.num_branches) {
			create_2nd_order_branch(branch_num, ++num_2_branches_created, 0, branch_deflected, rotation);
		}
	}
	if (branch.num_cylins < 2) {rotate_cylin(cylin);}
	
	for (int j = 1; j < branch.num_cylins; j++) {
		tree_cylin &cylin(branch.cylin[j]), &lcylin(branch.cylin[j-1]);
		branch_just_created = false;
		gen_next_cylin(cylin, lcylin, branch_1_var, branch_1_rad_var, 1, branch_num, (j < branch.num_cylins-1));
		float temp_num2;

		//calculate p1, based on the end point of last cylinder--------------------
		if (j < 0.5*branch.num_cylins) {
			temp_num2 = (j+2)*temp_num_big_branches/((num_2_branches_created+1)*(0.5f*branch.num_cylins+1));
		}
		else {
			int const denom(((num_2_branches_created+1)*branch.num_cylins));
			temp_num2 = float(int((j+1)*int(branch.num_branches)/denom) + int(branch.num_branches)/denom);
		}
		if (temp_num2*branch_1_distribution >= 1.0f && num_2_branches_created < branch.num_branches) branch_just_created = true;
		int const deg_added(generate_next_cylin(j, branch.num_cylins, branch_just_created, branch_deflected));
		cylin.deg_rotate += deg_added;
		cylin.deg_rotate *= branch_upwardness;

		//set rotate ------------
		gen_cylin_rotate(cylin.rotate, lcylin.rotate, deg_added*rotate_factor);
		lcylin.p2 = lcylin.p1 + re_matrix;
		if (j == (branch.num_cylins - 1)) {rotate_cylin(cylin);}
		
		if (branch_just_created) {
			if (j == branch.num_cylins-1) {
				branches[branch_num][num_2_branches_created+1].clear_num();
			}
			else {
				rotation = -rotation;
				create_2nd_order_branch(branch_num, ++num_2_branches_created, j, branch_deflected, rotation);
			}
		}
	} // for j
}


void tree_builder_t::create_2nd_order_branch(int i, int j, int cylin_num, bool branch_deflected, int rotation) {
	
	int index(0), num_3_branches_created(0);
	float const rotate_start((float)rgen.rand_int(0, int(max_2_angle_rotate)));
	bool branch_just_created(false);
	
	//first cylinders of the second order branches
	tree_branch &branch(branches[i][j]);
	tree_cylin &cylin(branch.cylin[0]), &src_cylin(branches[i][0].cylin[cylin_num]);
	branch.num_cylins = branches[i][0].num_cylins - cylin_num;
	gen_first_cylin(cylin, src_cylin, branch_2_start, branch_2_rad_var, rotate_start, 2, j);
	gen_cylin_rotate(cylin.rotate, src_cylin.rotate, rotate_start*rotation);

	//generate stats for the third order branches
	branch.num_branches = rgen.rand_int(num_3_branches_min, num_3_branches_max);
	float const temp_num(2.0*branch.num_branches/((float)branch.num_cylins));

	if (temp_num*branch_1_distribution >= 1.0f) {
		if (num_3_branches_created < branch.num_branches && cylin.r2 > branch_4_max_radius) {
			create_3rd_order_branch(i, j, 0, ((num_3_branches_created++) + num_34_branches[0]), branch_deflected, rotation);
			branch_just_created = true;
		}
	}
	if (branch.num_cylins < 2) {rotate_cylin(cylin);}

	//create rest of the cylinders
	for (index = 1; index < branch.num_cylins; index++) {
		tree_cylin &cylin(branch.cylin[index]), &lcylin(branch.cylin[index-1]);
		gen_next_cylin(cylin, lcylin, branch_2_var, branch_2_rad_var, 2, j, (index < branch.num_cylins-1));
		//float const temp_num2(float((index+2)*branch.num_branches)/((num_3_branches_created+1)*((float)branch.num_cylins)));
		if (temp_num*branch_1_distribution >= 1.0f && num_3_branches_created < branch.num_branches) {branch_just_created = true;}
		int const deg_added(generate_next_cylin(index, branch.num_cylins, branch_just_created, branch_deflected));
		gen_cylin_rotate(cylin.rotate, lcylin.rotate, deg_added*rotate_factor);
		cylin.rotate.z    = lcylin.rotate.z;
		cylin.deg_rotate += ((branch.cylin[0].deg_rotate < 0.0) ? deg_added : -deg_added);
		cylin.deg_rotate *= branch_upwardness;
		lcylin.p2 = lcylin.p1 + re_matrix;
		if (index == (branch.num_cylins - 1)) {rotate_cylin(cylin);}
		
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
	float const rotate_start((float)rgen.rand_int(0, int(max_3_angle_rotate)));
	bool branch_just_created(false);

	//first cylinders of the third order branches
	tree_branch &branch(branches_34[0][branch_num]);
	tree_cylin &cylin(branch.cylin[0]), &src_cylin(branches[i][j].cylin[cylin_num]);
	branch.num_cylins = branches[i][j].num_cylins - cylin_num;
	gen_first_cylin(cylin, src_cylin, branch_2_start, branch_2_rad_var, rotate_start, 3, branch_num);
	gen_cylin_rotate(cylin.rotate, src_cylin.rotate, rotate_start*rotation);

	//generate stats for the third order branches
	branch.num_branches = rgen.rand_int(num_3_branches_min, num_3_branches_max);
	if (branch.num_cylins < 2) {rotate_cylin(cylin);}

	//create rest of the cylinders
	for (index = 1; index < branch.num_cylins; index++) {
		tree_cylin &cylin(branch.cylin[index]), &lcylin(branch.cylin[index-1]);
		gen_next_cylin(cylin, lcylin, branch_2_var, branch_2_rad_var, 3, branch_num, (index < branch.num_cylins-1));
		cylin.rotate = lcylin.rotate;
		int const deg_added(generate_next_cylin(index, branch.num_cylins, branch_just_created, branch_deflected));
		cylin.deg_rotate += ((branch.cylin[0].deg_rotate < 0.0) ? deg_added : -deg_added);
		lcylin.p2 = lcylin.p1 + re_matrix;
		if (index == (branch.num_cylins - 1)) {rotate_cylin(cylin);}
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
				rotate_start += 360.0/num_4_branches*l + rgen.rand_int(1,30);
				generate_4th_order_branch(branch, j, rotate_start, temp_deg, branch_num++);
			}
		}
	}
}


void tree_builder_t::create_4th_order_branches(int nbranches, float branch_scale) {

	int num_4_branches  = 2;
	branch_4_length     = 6.0*branch_scale; //0.03;
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
		src_cylin.deg_rotate + ((src_cylin.deg_rotate > 0.0) ? 1 : -1)*rgen.rand_int(0,60));
	rotate_around_axis(src_cylin);
	cylin.p1 = src_cylin.p1 + BRANCH_LEN_SCALE*re_matrix;
	if (branch.num_cylins < 2) {rotate_cylin(cylin);}
	
	for (index = 1; index < branch.num_cylins; index++) {
		tree_cylin &cylin(branch.cylin[index]), &lcylin(branch.cylin[index-1]);
		gen_next_cylin(cylin, lcylin, branch_4_var, branch_4_rad_var, 4, branch_num, (index < branch.num_cylins-1));
		int const deg_added(generate_next_cylin(index, branch.num_cylins, false, branch_deflected));
		gen_cylin_rotate(cylin.rotate, lcylin.rotate, deg_added*rotate_factor);
		cylin.rotate.z = lcylin.rotate.z;
		cylin.deg_rotate += ((branch.cylin[0].deg_rotate < 0.0) ? deg_added : -deg_added);
		lcylin.p2 = lcylin.p1 + re_matrix;
		if (index == (branch.num_cylins - 1)) {rotate_cylin(cylin);}
	}
}


void tree_leaf::create_init_color(rand_gen_t &rgen) {

	lcolor = 1000;
	lred   = rgen.rand() & 255;
	lgreen = rgen.rand() & 255;
}


void tree_builder_t::add_leaves_to_cylin(unsigned cylin_ix, int tree_type, float rel_leaf_size, float deadness, vector<tree_leaf> &leaves) {

	assert(cylin_ix < cylin_cache.size());
	tree_cylin const &cylin(cylin_cache[cylin_ix]);
	leaf_acc += num_leaves_per_occ;
	int const temp((int)leaf_acc);
	leaf_acc -= (float)temp;
	float const temp_deg(((cylin.rotate.y < 0.0) ? -1.0 : 1.0)*safe_acosf(cylin.rotate.x)), angle_step(360.0/temp);

	for (int l = 0; l < temp; l++) {
		if (deadness > 0 && deadness > rgen.rand_float()) continue;
		rotate_around_axis(cylin); // rotate leaf starting point
		point const start(cylin.p1 + 0.9*re_matrix);
		float const rotate_start(angle_step*l);
		vector3d rotate;
		setup_rotate(rotate, (rotate_start + rgen.rand_int(1,60)), temp_deg); // rotate is in XY plane (z=0)
		float const deg_rotate(cylin.deg_rotate + ((cylin.deg_rotate > 0.0) ? 1 : -1)*rgen.rand_int(0,60));
		float const lsize(rel_leaf_size*tree_types[tree_type].leaf_x_ar*(0.7*rgen.randd() + 0.3));
		tree_leaf leaf;

		for (int p = 0; p < 4; ++p) {
			// calculate re_matrix = offset from leaf starting point for this vertex; deg_rotate is X rotation applied after Z rotation
			rotate_pts_around_axis(leaf_points[p]*lsize, rotate, deg_rotate);
			leaf.pts[p] = start + re_matrix;
		}
		point const base_pt((leaf.pts[0] + leaf.pts[3])*0.5);
		vector3d const xlate(cylin.r1*(base_pt - cylin.p1).get_norm());
		for (unsigned i = 0; i < 4; ++i) {leaf.pts[i] += xlate;} // move away from the branch centerline by radius
		leaf.create_init_color(rgen);
		leaf.norm = cross_product((leaf.pts[1] - leaf.pts[0]), (leaf.pts[3] - leaf.pts[0])).get_norm();
		leaves.push_back(leaf);
	} // for l
}


int tree_builder_t::generate_next_cylin(int cylin_num, int ncib, bool branch_just_created, bool &branch_deflected) {

	//vars used in generating deg_rotate for cylinders
	float const PI_16(PI/16.0);
	float const t_start(TWO_PI/rgen.rand_int(3,8)); //start in the first pi
	float const t_end(((rgen.rand_int(1,3) == 1) ? 1.0f : 5.0f)*PI_TWO + rgen.rand_int(2,8)*PI_16); //either PI/2 to PI or 5*PI/2 to 3*PI - controls branch droopiness
	int add_deg_rotate(int(0.01*cylin_num*sinf(t_start+(t_end-t_start)*cylin_num/ncib)*branch_curveness)); //how much wavy the branch will be --in degrees
	int rg[2] = {};

	if (cylin_num < int(ncib/3)) { //how much to start the starting deg_scale
		rg[0] = 5; rg[1] = 10;
	}
	else if (cylin_num < int(ncib*2/3)) { //scale for trig function -- to scale the middle of a branch
		rg[0] = 5; rg[1] = 10;
	}
	else { //how much to end the ending deg_scale
		rg[0] = 1; rg[1] = 5;
	}
	add_deg_rotate  *= rgen.rand_int(rg[0], rg[1]);
	branch_deflected = false;

	if (branch_just_created) {
		if (rgen.rand_int(0,10) < 6) return add_deg_rotate;  //don't deflect
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


void tree_cont_t::gen_deterministic(int x1, int y1, int x2, int y2, float vegetation_, float mesh_dz, tile_t const *const cur_tile) { // default full tile generation function
	gen_trees_tt_within_radius(x1, y1, x2, y2, all_zeros, 0.0, 0, mesh_dz, cur_tile, vegetation_, 1); // not using bounding sphere
	//cout << TXT(mod_num_trees) << TXT(size()) << endl;
}

void tree_cont_t::add_new_tree(rand_gen_t &rgen, int &ttype) {

	push_back(tree());
	if (shared_tree_data.empty()) return; // no fixed ID
	int tree_id(-1);

	if (ttype >= 0) {
		unsigned const num_per_type(max(1U, (unsigned)shared_tree_data.size()/NUM_TREE_TYPES));
		tree_id = min(unsigned((((rgen.rseed1 >> 7) + rgen.rseed2) % num_per_type) + ttype*num_per_type), (unsigned)shared_tree_data.size()-1);
	}
	else {
		tree_id = (rgen.rseed2 % shared_tree_data.size());
		ttype   = tree_id % NUM_TREE_TYPES;
	}
	if (shared_tree_data[tree_id].is_created()) {ttype = shared_tree_data[tree_id].get_tree_type();} // in case there weren't enough generated to get the requested type
	//cout << "selected tree " << tree_id << " of " << shared_tree_data.size() << " type " << ttype << endl;
	if (tree_id >= 0) {back().bind_to_td(&shared_tree_data[tree_id]);}
}

bool tree_placer_t::have_small_trees() const {return (world_mode == WMODE_INF_TERRAIN && !tree_placer.blocks   .empty());}
bool tree_placer_t::have_decid_trees() const {return (world_mode == WMODE_INF_TERRAIN && !tree_placer.sm_blocks.empty());}

void tree_placer_t::add(point const &pos, float size, int type, bool allow_bush, bool add_bush, bool is_sm_tree) {
	if ((is_sm_tree ? sm_blocks : blocks).empty()) {begin_block(is_sm_tree);} // begin a new block in case user isn't creating the blocks themselves
	tree_block &block((is_sm_tree ? sm_blocks : blocks).back());
	if (block.trees.empty()) {block.bcube.set_from_point(pos);} else {block.bcube.union_with_pt(pos);}
	(is_sm_tree ? sm_bcube : bcube).assign_or_union_with_pt(pos);
	block.trees.emplace_back(pos, size, type, allow_bush, add_bush);
}

void tree_cont_t::gen_trees_tt_within_radius(int x1, int y1, int x2, int y2, point const &center, float radius, bool is_square,
	float mesh_dz, tile_t const *const cur_tile, float vegetation_, bool use_density)
{
	//timer_t timer("Gen Trees");
	bool const NONUNIFORM_TREE_DEN = 1; // based on world_mode?
	unsigned const mod_num_trees(num_trees/(NONUNIFORM_TREE_DEN ? sqrt(tree_density_thresh) : 1.0));
	float const min_tree_h(water_plane_z + 0.01*zmax_est), max_tree_h(1.8*zmax_est);
	rand_gen_t rgen;
	generated = 1;

	if (tree_placer.have_decid_trees()) { // now add pre-placed trees within the city (TT mode)
		shared_tree_data.ensure_init();
		vector3d const xlate(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0);
		cube_t const bounds(get_xval(x1), get_xval(x2), get_yval(y1), get_yval(y2), min_tree_h, max_tree_h); // Note: zvals are unused

		if (bounds.intersects_xy(tree_placer.bcube + xlate)) { // test full bcube
			for (auto b = tree_placer.blocks.begin(); b != tree_placer.blocks.end(); ++b) {
				if (!bounds.intersects_xy(b->bcube + xlate)) continue;

				for (auto t = b->trees.begin(); t != b->trees.end(); ++t) {
					point const pos(t->pos + xlate);
					if (!bounds.contains_pt_xy(pos)) continue; // tree not within this tile
					int ttype(t->type);
					if (ttype >= 0) {ttype %= NUM_TREE_TYPES;} // make sure it maps to a valid tree type if specified
					add_new_tree(rgen, ttype);
					// Note: can't be user placed + instanced
					back().gen_tree(pos, int(t->size), ttype, 1, 1, 0, rgen, 1.0, 1.0, 1.0, tree_4th_branches, t->allow_bush, t->force_bush);
				} // for t
			} // for b
		}
	}
	if (mod_num_trees == 0 || !(tree_mode & 1)) {calc_bcube(); return;} // no generated trees
	float const height_thresh(get_median_height(tree_density_thresh));
	unsigned const smod(3.321*XY_MULT_SIZE+1), tree_prob(max(1U, XY_MULT_SIZE/mod_num_trees));
	unsigned const skip_val(max(1, int(1.0/tree_scale))); // similar to deterministic gen in scenery.cpp
	shared_tree_data.ensure_init();
	mesh_xy_grid_cache_t density_gen[NUM_TREE_TYPES+1];

	if (NONUNIFORM_TREE_DEN) { // i==0 is the coverage density map, i>0 are the per-tree type coverage maps
#pragma omp parallel for schedule(dynamic) num_threads(2)
		for (int i = (use_density ? 0 : 1); i <= NUM_TREE_TYPES; ++i) { // Note: i should be signed
			float const tds(TREE_DIST_SCALE*(XY_MULT_SIZE/16384.0)*(i==0 ? 1.0 : 0.1)), xscale(tds*DX_VAL*DX_VAL), yscale(tds*DY_VAL*DY_VAL);
			density_gen[i].build_arrays(xscale*(x1 + xoff2 + 1000*i), yscale*(y1 + yoff2 - 1500*i), xscale, yscale, (x2-x1), (y2-y1), 0, 1); // force_sine_mode=1
		}
	}
	for (int i = y1; i < y2; i += skip_val) {
		float const yval(get_yval(i));
		if (radius > 0.0 && fabs(yval - center.y) > radius) continue;

		for (int j = x1; j < x2; j += skip_val) {
			if (radius > 0.0) {
				float const xval(get_xval(j));
				if (fabs(xval - center.x) > radius) continue;
				point const tpos(xval, yval, 0.0);
				if (!dist_xy_less_than(center, tpos, radius)) continue; // Note: uses mesh xy center, not actual tree pos (for simplicity and efficiency)
			}
			if (scrolling) {
				int const ox(j + dx_scroll), oy(i + dy_scroll); // positions in original coordinate system
				if (ox >= x1 && ox <= x2 && oy >= y1 && oy <= y2) continue; // use orignal tree from last position
			}
			rgen.rseed1 = 805306457*(i + yoff2) + 12582917*(j + xoff2) + 100663319*rand_gen_index;
			rgen.rseed2 = 6291469  *(j + xoff2) + 3145739 *(i + yoff2) + 1572869  *rand_gen_index;
			rgen.rand_mix();
			unsigned const val(((unsigned)rgen.rand_seed_mix())%smod);
			if (val <= 100)         continue; // scenery
			if (val%tree_prob != 0) continue; // not selected
			if ((rgen.rseed1&127)/128.0 >= vegetation_) continue;
			point pos((get_xval(j) + 0.5*DX_VAL*rgen.randd()), (yval + 0.5*DY_VAL*rgen.randd()), 0.0);
			// Note: pos.z will be slightly different when calculated within vs. outside the mesh bounds
			pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
			if (pos.z > max_tree_h || pos.z < min_tree_h) continue;
			if (tree_mode == 3 && get_tree_class_from_height(pos.z, 0) != TREE_CLASS_DECID) continue; // use a pine tree here (or no tree)
			if (point_inside_voxel_terrain(pos)) continue; // don't create trees that start inside voxels (but what about trees that grow into voxels?)
			int ttype(-1);

			if (NONUNIFORM_TREE_DEN) {
				if (use_density && density_gen[0].eval_index(j-x1, i-y1) > height_thresh) continue; // density function test
				float max_val(0.0);

				for (unsigned tt = 0; tt < NUM_TREE_TYPES; ++tt) {
					float den_val(density_gen[tt+1].eval_index(j-x1, i-y1)); // no glaciate
					float const jitter(extract_low_bits_01(den_val, 100.0)); // randomize to make the border more fuzzy
					den_val += 0.8*jitter*jitter; // quadratic distribution
					if (max_val == 0.0 || den_val > max_val) {max_val = den_val; ttype = tt;}
				}
			}
			if (!check_valid_scenery_pos((pos + vector3d(0.0, 0.0, 0.3*tree_scale)), 0.4*tree_scale, 1)) continue; // approximate bsphere; is_tall=1
			// if terrain is highly varying in height, calculate surrounding mesh zvals and use the min for tree height; if slope is too high, skip this tree
			if (mesh_dz < 0.0 || mesh_dz > 1.0) {
				if (!adjust_tree_zval(pos, 0, ttype, 0, cur_tile)) continue; // create_bush=0
			}
			add_new_tree(rgen, ttype);
			back().gen_tree(pos, 0, ttype, 0, 1, 0, rgen, 1.0, 1.0, 1.0, tree_4th_branches, 1); // allow bushes
		} // for j
	} // for i
	calc_bcube();
}


void regen_trees(bool keep_old) {

	if (!scrolling) {cout << "vegetation: " << vegetation << endl;}
	RESET_TIME;
	static int init(0), last_rgi(0), last_xoff2(0), last_yoff2(0);
	static float last_ts(0.0);
	if (tree_mode & 2) {gen_small_trees();}
	//else {remove_small_tree_cobjs();}

	if ((tree_mode & 1) && num_trees > 0) {
		if (keep_old && init && last_rgi == rand_gen_index && last_xoff2 == xoff2 && last_yoff2 == yoff2 && last_ts == tree_scale) {
			add_tree_cobjs(); // keep old trees
			return;
		}
		int const border(1), ext_x1(border), ext_x2(MESH_X_SIZE-border), ext_y1(border), ext_y2(MESH_Y_SIZE-border);
		if (scrolling && t_trees.scroll_trees(ext_x1, ext_x2, ext_y1, ext_y2)) {t_trees.post_scroll_remove();}
		else {t_trees.resize(0);}
		t_trees.gen_deterministic(ext_x1, ext_y1, ext_x2, ext_y2, vegetation, /*(zmax - zmin)*/-1.0); // don't use mesh_dz (may cause problems with scrolling)
		if (!scrolling) {cout << "Num trees = " << t_trees.size() << endl;}
		last_rgi   = rand_gen_index;
		last_xoff2 = xoff2;
		last_yoff2 = yoff2;
		last_ts    = tree_scale;
		init       = 1;
	}
	if (!scrolling) {PRINT_TIME(" Gen Trees");}
}


void tree::write_to_cobj_file(std::ostream &out) const {

	tree_data_t const &td(tdata());
	tree_type const &tt(tree_types[type]);
	float const size(td.base_radius/(0.1*TREE_SIZE*tt.branch_size/tree_scale)); // reverse engineer these parameters from known values 
	out << "g " << td.b_tex_scale*td.br_scale/tt.branch_tscale << " " << td.br_scale/branch_radius_scale << " " << tree_nl_scale << " " << enable_leaf_wind << endl;

	if (use_clip_cube) {
		assert(!tree_4th_branches);
		// 'H': place hedges: xstart ystart dx dy nsteps size, type [cx1 cx2 cy1 cy2 cz1 cz2]
		out << "H " << tree_center.x << " " << tree_center.y << " 0.0 0.0 1 " << size << " " << type << " " << clip_cube.raw_str() << endl;
	}
	else {
		// 'E': // place tree: xpos ypos size type [zpos [tree_4th_branches]], type: TREE_MAPLE = 0, TREE_LIVE_OAK = 1, TREE_A = 2, TREE_B = 3, 4 = TREE_PAPAYA
		//fscanf(fp, "%f%f%f%i%f%i", &pos.x, &pos.y, &fvals[0], &ivals[0], &pos.z, local_tree_4th_branches)
		out << "E " << tree_center.x << " " << tree_center.y << " " << size << " " << type << " " << tree_center.z << " " << tree_4th_branches << endl;
	}
}
void write_trees_to_cobj_file(std::ostream &out) {
	for (auto i = t_trees.begin(); i != t_trees.end(); ++i) {i->write_to_cobj_file(out);}
}


bool tree_cont_t::update_zvals(int x1, int y1, int x2, int y2) {

	bool updated(0);

	for (iterator i = begin(); i != end(); ++i) {
		point const &pos(i->get_center());
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
		if (xpos < x1 || xpos > x2 || ypos < y1 || ypos > y2) continue;
		float const new_z(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1)); // use-real_equation=0
		if (fabs(pos.z - new_z) < 0.01*i->get_radius()) continue;
		i->remove_collision_objects();
		i->shift_tree(point(0.0, 0.0, (new_z - pos.z)));
		if (tree_mode & 1) {i->add_tree_collision_objects();}
		updated = 1;
	}
	return updated;
}


void tree_cont_t::spraypaint_leaves(point const &pos, float radius, colorRGBA const &color) {

	for (iterator i = begin(); i != end(); ++i) {
		if (dist_less_than(pos, i->sphere_center(), (radius + i->get_radius()))) {i->spraypaint_leaves(pos, radius, color);}
	}
}


void tree_data_manager_t::ensure_init() {

	if (max_unique_trees > 0 && empty()) {
		resize(max_unique_trees);
	}
	else if (tree_scale != last_tree_scale || rand_gen_index != last_rgi) {
		for (iterator i = begin(); i != end(); ++i) {i->clear_data();}
	}
	last_tree_scale = tree_scale;
	last_rgi        = rand_gen_index;
}

void tree_data_manager_t::clear_context() {
	for (iterator i = begin(); i != end(); ++i) {i->clear_context();}
}

void tree_data_manager_t::on_leaf_color_change() {
	for (iterator i = begin(); i != end(); ++i) {i->on_leaf_color_change();}
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

unsigned tree_cont_t::get_closest_tree_type(point const &pos) const {
	int closest_type(TREE_MAPLE); // default value for the no-trees-found case
	float min_dist_sq(0.0);
	
	for (const_iterator i = begin(); i != end(); ++i) {
		int const type(i->get_type());
		if (type < 0) continue; // tree type not set
		float const dist_sq(p2p_dist_sq(pos, i->get_center()));
		
		if (min_dist_sq == 0.0 || dist_sq < min_dist_sq) {
			min_dist_sq = dist_sq; closest_type = type;
			float const r(i->get_radius());
			if (dist_sq < 2.0*r*r) break; // close enough, early terminate
		}
	} // for i
	assert(closest_type >= 0 && closest_type < NUM_TREE_TYPES);
	return closest_type;
}

void tree_cont_t::update_zmax(float &tzmax) const {
	for (const_iterator i = begin(); i != end(); ++i) {tzmax = max(tzmax, (i->get_center().z + i->get_radius()));}
}

void tree_cont_t::shift_by(vector3d const &vd) {
	for (iterator i = begin(); i != end(); ++i) {i->shift_tree(vd);}
}

void tree_cont_t::add_cobjs() {
	for (iterator i = begin(); i != end(); ++i) {i->add_tree_collision_objects();}
	calc_bcube();
}

void tree_cont_t::calc_bcube() {
	all_bcube.set_to_zeros();
	for (iterator i = begin(); i != end(); ++i) {i->add_bounds_to_bcube(all_bcube);}
}

void tree_cont_t::clear_context() {
	for (iterator i = begin(); i != end(); ++i) {i->clear_context();}
}

void tree_cont_t::check_render_textures() {
	//timer_t timer("Check Render Textures"); // 1463 total, 239 max
	for (iterator i = begin(); i != end(); ++i) {i->check_render_textures();}
}

void tree_cont_t::apply_exp_damage(point const &epos, float damage, float bradius, int type) {
	blastr const br(0, ETYPE_FIRE, NO_SOURCE, bradius, damage, epos, plus_z, YELLOW, RED);
	for (iterator i = begin(); i != end(); ++i) {i->blast_damage(&br);}
}


void shift_trees(vector3d const &vd) {
	if (num_trees > 0) return; // dynamically created, not placed
	t_trees.shift_by(vd);
}

bool update_decid_tree_zvals(int x1, int y1, int x2, int y2) {
	return t_trees.update_zvals(x1, y1, x2, y2);
}

void spraypaint_tree_leaves(point const &pos, float radius, colorRGBA const &color) {
	t_trees.spraypaint_leaves(pos, radius, color);
}

void exp_damage_trees(point const &epos, float damage, float bradius, int type) {
	t_trees.apply_exp_damage(epos, damage, bradius, type);
}

void add_tree_cobjs   () {t_trees.add_cobjs();}
void remove_tree_cobjs() {t_trees.remove_cobjs();}

void clear_tree_context() {
	t_trees.clear_context();
	tree_data_manager.clear_context();
}


// tree fires

unsigned const MAX_BRANCH_BURN_LEVEL = 3;

void tree::add_fire(point const &pos, float radius, float val, bool spread_mode) {

	if (spread_mode && tree_fire != nullptr && tree_fire->is_burning()) return; // this tree is already on fire
	if (val < 100.0 && !spread_mode) return; // too small optimization
	if (!enable_grass_fire || !physics_enabled()) return;
	if (!dist_less_than(pos, sphere_center(), (radius + get_radius()))) return; // not within tree bounding sphere
	//make_private_tdata_copy(); // required if branches and/or leaves are modified (burned)
	if (tree_fire == nullptr) {tree_fire.reset(new tree_fire_t(tdata().get_all_cylins(), tree_center, tdata().base_radius));} // create if needed
	tree_fire->add_fire(pos, radius, 0.1*val);
}
void tree::next_fire_frame() {
	if (tree_fire != nullptr) {tree_fire->next_frame(*this);}
}
void tree::draw_fire(shader_t &s) const {
	if (tree_fire != nullptr && is_visible_to_camera()) {tree_fire->draw(s);}
}


void tree_cont_t::apply_fire(point const &pos, float radius, float val, bool spread_mode) {
	for (iterator i = begin(); i != end(); ++i) {i->add_fire(pos, radius, val, spread_mode);}
}

void tree_cont_t::next_fire_frame() {
	for (iterator i = begin(); i != end(); ++i) {i->next_fire_frame();}
}
void tree_cont_t::draw_fire(shader_t &s) const {
	//timer_t timer("Draw Tree Fires");
	for (const_iterator i = begin(); i != end(); ++i) {i->draw_fire(s);}
}
bool tree_cont_t::has_any_fire() const {
	for (const_iterator i = begin(); i != end(); ++i) {if (i->is_on_fire()) return 1;}
	return 0;
}

void apply_tree_fire(point const &pos, float radius, float val, bool spread_mode) {t_trees.apply_fire(pos, radius, val, spread_mode);}
void next_frame_tree_fires() {t_trees.next_fire_frame();}
void draw_tree_fires(shader_t &s) {t_trees.draw_fire(s);}
bool any_trees_on_fire() {return t_trees.has_any_fire();}


tree_fire_t::tree_fire_t(vector<draw_cylin> const &branches_, point const &tree_center_, float tree_base_radius) :
	branches(branches_), tree_center(tree_center_), update_ix(0), has_fire(0)
{
	rand_gen_t rgen;
	unsigned num(branches.size());

	for (vector<draw_cylin>::const_iterator i = branches.begin(); i != branches.end(); ++i) {
		if (i->level > MAX_BRANCH_BURN_LEVEL) {num = i - branches.begin(); break;} // end when we reach branches of a high enough level
	}
	fires.resize(num);

	for (unsigned i = 0; i < num; ++i) { // gen random params
		fires[i].fuel = 1.0E6*rgen.rand_uniform(1.0, 1.5)*branches[i].get_volume() + 0.01; // ~500 for trunks, 0.01 for tiny branches
		fires[i].hp   = min(100.0, 1.0E4*rgen.rand_uniform(0.5, 1.0)*branches[i].get_avg_radius()); // ~300 for trunks, 2.0 for tiny branches
		fires[i].branch_bradius = branches[i].get_bounding_radius();
	}
	fire_radius = max(0.5*HALF_DXY, 3.0*tree_base_radius); // scale to tree size
}

void tree_fire_t::shift(vector3d const &vd) {
	for (auto i = fires.begin(); i != fires.end(); ++i) {i->pos += vd;}
}

void update_dist_to_fire(point const &pos, float dist_mult);

void tree_fire_t::next_frame(tree &t) {

	if (!has_fire) return;
	//timer_t timer("Tree Fire Next Frame");
	float const burn_rate(fire_elem_t::get_burn_rate()), spread_rate(1.2*fticks*burn_rate);
	unsigned const num_fires(fires.size());
	unsigned const counter_offset(frame_counter + uintptr_t(this)/sizeof(void *)); // make offset unique per tree
	assert(num_fires <= branches.size());
	has_fire = 0;
	rand_gen_t rgen;

	for (unsigned i = 0; i < num_fires; ++i) {
		tree_fire_elem_t &elem(fires[i]);
		if (elem.burn_amt == 0.0) continue; // optimization
		elem.next_frame(4.0*burn_rate, branches[i].get_surface_area(), 0.2);
		if (elem.burn_amt == 0.0) continue; // burned out
		draw_cylin const &cylin(branches[i]);
		if (is_underwater(cylin.p1 + tree_center)) {elem.burn_amt = 0.0; continue;} // skip roots that are underwater
		has_fire = 1;
		update_dist_to_fire(elem.pos, 1.0);
		if (elem.burn_amt < 0.5) continue; // not large enough to spread or do damage
		int const counter(i + counter_offset);
		if ((counter&3) != 0) continue; // update every 4 frames as an optimization
		bool const trunk(cylin.level == 0); // trunk fire spreads more quickly
		float const radius(elem.burn_amt*fire_radius*rgen.rand_uniform(0.8, 1.3)), burn_radius(radius + elem.branch_bradius);
		vector3d const dir(rgen.signed_rand_vector_spherical().get_norm() + 0.2*wind + vector3d(0, 0, 0.5)); // add minor wind influence; spread is biased upward
		point const pos(elem.pos + radius*dir);

		if (elem.sleep_time > 0) {--elem.sleep_time;} // sleeping
		else {
			int const ret(add_fire(pos, burn_radius*(trunk ? 1.5 : 1.0), spread_rate*elem.burn_amt*(trunk ? 2.0 : 1.0))); // expand fire to cover the entire branch
			if (ret == 0) {elem.sleep_time = 8;} // sleep for 8 frames
		}
		if ((counter&15  ) == 0) {fire_damage.emplace_back(pos, 1.5*burn_radius, 0.01*fticks*elem.burn_amt);} // update every 16 frames as an optimization
		if ((counter&1023) == 0) {gen_smoke(elem.pos, 1.0, 1.0, colorRGBA(0.3, 0.3, 0.3, 0.4), 1);} // no_lighting=1

		if (trunk || (cylin.level == 1 && elem.pos.z < interpolate_mesh_zval(elem.pos.x, elem.pos.y, 0.0, 0, 1))) { // trunk or below the mesh
			add_ground_fire(elem.pos, radius, 20.0);
		}
		if ((counter&127) == 0) {apply_tree_fire(elem.pos, radius, 100.0*spread_rate*elem.burn_amt, 1);} // occasionally spread to other trees
	} // for i
	t.apply_fire_damage(fire_damage, 4); // skipval=4
	fire_damage.clear();
}

int tree_fire_t::add_fire(point const &pos, float radius, float val) { // returns: 0=no change, 1=damage only, 2=create fire

	if (!animate2) return 0;
	if (val == 0.0 || radius == 0.0) return 0; // no fire
	assert(radius > 0.0);
	point const rel_pos(pos - tree_center); // translate to tree local coords
	float const heat_amt(0.5*val), radius_sq(radius*radius);
	unsigned const skipval((val > 50.0) ? 1 : 2); // optimization
	unsigned const start_ix(update_ix%skipval), num_fires(fires.size());
	unsigned best_ix(0);
	float dmin_sq(0.0);
	bool damaged(0);
	assert(num_fires <= branches.size());
	++update_ix;

	for (unsigned i = start_ix; i < num_fires; i += skipval) { // iterate over branches
		float &hp(fires[i].hp);
		if (hp <= 0.0 || fires[i].fuel == 0.0) continue; // already burning/burned
		draw_cylin const &cylin(branches[i]);
		float const dist_sq(p2p_dist_sq(cylin.get_center(), rel_pos)); // min distance
		if (dist_sq < radius_sq) {hp -= min(heat_amt, 0.5f*hp);} // apply heat damage to nearby branches
		damaged = 1;
		if (hp > val) continue; // not enough heat to burn, skip
		float const sep_dist(radius + fires[i].branch_bradius);
		if (dist_sq > sep_dist*sep_dist) continue; // bounding sphere test
		
		if ((dmin_sq == 0.0 || dist_sq < dmin_sq) && sphere_intersect_cylinder(rel_pos, radius, cylin.p1, cylin.p2, cylin.r1, cylin.r2, 1)) { // check_ends=1
			best_ix = i;
			dmin_sq = dist_sq; // choose to burn the nearest branch
		}
	} // for i
	if (dmin_sq == 0.0) return damaged; // no fire was started
	tree_fire_elem_t &elem(fires[best_ix]);
	//cout << TXT(radius) << TXT(val) << TXT(dmin_sq) << TXT(best_ix) << TXT(elem.hp) << TXT(elem.fuel);
	if (!elem.burn(val)) return 1;
	rand_gen_t rgen;
	rgen.set_state(best_ix, update_ix);
	draw_cylin const &branch(branches[best_ix]);
	float const t(rgen.rand_float()), branch_radius(branch.r1 + (branch.r2 - branch.r1)*t);
	vector3d const branch_dir(branch.p2 - branch.p1);
	elem.pos = branch.p1 + branch_dir*t + tree_center; // random point along branch centerline
	vector3d dir;
	orthogonalize_dir(rgen.signed_rand_vector_spherical(), branch_dir, dir, 1);
	if (dir.z < -0.2) {dir.z = -dir.z;} // don't point too much downward
	elem.pos += branch_radius*dir;
	has_fire = 1;
	return 2; // fire was started
}

void tree_fire_t::draw(shader_t &s) const {

	if (!has_fire) return;
	//timer_t timer("Draw Tree Fire");
	fire_drawer_t fire_drawer; // cache inside tree_fire_t for reuse?

	for (vector<tree_fire_elem_t>::const_iterator i = fires.begin(); i != fires.end(); ++i) {
		if (i->burn_amt == 0.0) continue; // not burning
		if (!camera_pdu.sphere_visible_test(i->pos, 2.0*fire_radius)) continue; // VFC
		int const frame_ix(123*(i - fires.begin()));
		float const radius(0.5*fire_radius*i->burn_amt);
		fire_drawer.add_fire((i->pos + vector3d(0.0, 0.0, 0.75*radius)), radius, frame_ix, 1.0); // one flame per fire
	}
	fire_drawer.draw(s);
}

