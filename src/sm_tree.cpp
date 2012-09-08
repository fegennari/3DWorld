// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 2/13/03

#include "3DWorld.h"
#include "mesh.h"
#include "small_tree.h"
#include "gl_ext_arb.h"
#include "shaders.h"


float const SM_TREE_SIZE    = 0.05;
float const TREE_DIST_SCALE = 100.0;
float const TREE_DIST_RAND  = 0.2;
float const TREE_DIST_MH_S  = 100.0;
float const SM_TREE_AMT     = 0.55;
float const SM_TREE_QUALITY = 1.0;
float const LINE_THRESH     = 700.0;
bool const SMALL_TREE_COLL  = 1;
bool const DRAW_COBJS       = 0; // for debugging
bool const USE_BUMP_MAP     = 0;
int  const NUM_SMALL_TREES  = 40000;
unsigned const N_PT_LEVELS  = 6;
unsigned const N_PT_RINGS   = 5;


small_tree_group small_trees;
pt_line_drawer tree_scenery_pld;

extern bool disable_shaders, has_dir_lights;
extern int window_width, shadow_detail, draw_model, island, num_trees, do_zoom, tree_mode, xoff2, yoff2;
extern int rand_gen_index, display_mode, force_tree_class;
extern float zmin, zmax_est, water_plane_z, tree_scale, vegetation, OCEAN_DEPTH;
extern GLUquadricObj* quadric;


struct sm_tree_type {

	int leaf_tid, bark_tid;
	float w2, ws, h, ss;
	colorRGBA c;

	sm_tree_type(float w2_, float ws_, float h_, float ss_, colorRGBA const &c_, int ltid, int btid)
		: leaf_tid(ltid), bark_tid(btid), w2(w2_), ws(ws_), h(h_), ss(ss_), c(c_) {}
};

sm_tree_type const stt[NUM_ST_TYPES] = { // w2, ws, h, ss, c, tid
	sm_tree_type(0.00, 0.10, 0.35, 0.4, PTREE_C, PINE_TEX,      BARK2_TEX), // T_PINE
	sm_tree_type(0.13, 0.15, 0.75, 0.8, TREE_C,  TREE_HEMI_TEX, BARK3_TEX), // T_DECID // GRASS_TEX?
	sm_tree_type(0.13, 0.15, 0.75, 0.7, TREE_C,  GROUND_TEX,    BARK1_TEX), // T_TDECID
	sm_tree_type(0.00, 0.15, 0.00, 0.8, WHITE,   GROUND_TEX,    BARK4_TEX), // T_BUSH NOTE: bark texture is not used in trees, but is used in logs
	sm_tree_type(0.03, 0.15, 1.00, 0.6, TREE_C,  PALM_TEX,      BARK1_TEX), // T_PALM FIXME: PALM_BARK_TEX?
	sm_tree_type(0.00, 0.07, 0.00, 0.4, PTREE_C, PINE_TEX,      BARK2_TEX), // T_SH_PINE
};

int get_bark_tex_for_tree_type(int type) {
	assert(type < NUM_ST_TYPES);
	return stt[type].bark_tid;
}

colorRGBA get_tree_trunk_color(int type, bool modulate_with_texture) {
	assert(type < NUM_ST_TYPES);
	colorRGBA c(stt[type].c);
	
	if (modulate_with_texture && stt[type].bark_tid >= 0) {
		return c.modulate_with(texture_color(stt[type].bark_tid));
	}
	return c;
}


void small_tree_group::add_tree(small_tree &st) {

	if (st.is_pine_tree()) {
		max_pt_radius = max(max_pt_radius, st.get_pine_tree_radius());
		++num_pine_trees;
	}
	push_back(st);
}


void small_tree_group::calc_trunk_pts() {

	if (!trunk_pts.empty()) return; // already calculated
	trunk_pts.reserve(2*size());

	for (iterator i = begin(); i != end(); ++i) {
		i->add_trunk_as_line(trunk_pts);
	}
}


void small_tree_group::finalize(bool low_detail, bool pri_dim) {

	if (empty()) return;
	assert(!is_uploaded(low_detail));
	vbo_manager[low_detail].clear();
	vbo_manager[low_detail].reserve_pts(4*(low_detail ? 1 : N_PT_LEVELS*N_PT_RINGS)*num_pine_trees);

	for (iterator i = begin(); i != end(); ++i) {
		i->calc_points(vbo_manager[low_detail], low_detail, pri_dim);
	}
	last_pri_dim = pri_dim;
}


void small_tree_group::finalize_upload_and_clear_pts(bool low_detail, bool pri_dim) {

	if (low_detail && pri_dim != last_pri_dim) {clear_vbo_and_ids_if_needed(low_detail);}
	if (empty() || is_uploaded(low_detail)) return;
	//RESET_TIME;
	finalize(low_detail, pri_dim);
	//if (!low_detail) {PRINT_TIME("Finalize");}
	vbo_manager[low_detail].upload_and_clear_points();
	//if (!low_detail) {PRINT_TIME("Finalize + Upload");}
}


void small_tree_group::add_trunk_pts(point const &xlate, vector<point> &pts) const {

	for (vector<point>::const_iterator i = trunk_pts.begin(); i != trunk_pts.end(); ++i) {
		pts.push_back(*i + xlate);
	}
}


void small_tree_group::clear_vbos() {
	
	for (unsigned i = 0; i < 2; ++i) {
		vbo_manager[i].clear_vbo();
	}
}


void small_tree_group::clear_vbo_manager(int which) {
	
	for (unsigned d = 0; d < 2; ++d) {
		if (which & (1<<d)) {vbo_manager[d].clear();}
	}
}


void small_tree_group::clear_vbo_manager_and_ids(int which) {
	
	if (which & 1) { // high detail
		for (iterator i = begin(); i != end(); ++i) {i->clear_vbo_mgr_ix();}
	}
	clear_vbo_manager(which);
}


void small_tree_group::clear_vbo_and_ids_if_needed(bool low_detail) {
	
	if (is_uploaded(low_detail)) {clear_vbo_manager_and_ids(low_detail ? 2 : 1);}
}


void small_tree_group::clear_all() {

	clear();
	clear_vbo_manager_and_ids();
	trunk_pts.clear();
	generated     = 0;
	max_pt_radius = 0.0;
}


void small_tree_group::add_cobjs() {

	if (empty() || !SMALL_TREE_COLL) return;
	cobj_params cp      (0.65, GREEN, DRAW_COBJS, 0, NULL, 0, -1);
	cobj_params cp_trunk(0.9, TREE_C, DRAW_COBJS, 0, NULL, 0, -1);
	cp.shadow       = (shadow_detail >= 5);
	cp_trunk.shadow = (shadow_detail >= 6);

	for (iterator i = begin(); i != end(); ++i) {
		i->add_cobjs(cp, cp_trunk);
	}
}


void small_tree_group::remove_cobjs() {

	for (iterator i = begin(); i != end(); ++i) {
		i->remove_cobjs();
	}
}


void small_tree_group::translate_by(vector3d const &vd) {

	for (iterator i = begin(); i != end(); ++i) {
		i->translate_by(vd);
	}
}


void small_tree_group::draw_branches(bool shadow_only, vector3d const xlate, vector<point> *points) const {

	BLACK.do_glColor();

	for (const_iterator i = begin(); i != end(); ++i) {
		i->draw(1, shadow_only, 0, vbo_manager[0], xlate, points);
	}
}


void small_tree_group::draw_pine_leaves(bool shadow_only, bool low_detail, bool draw_all_pine, vector3d const xlate) const {

	vbo_vnc_block_manager_t const &vbomgr(vbo_manager[low_detail]);
	vbomgr.begin_render(1);
	select_texture((draw_model != 0) ? WHITE_TEX : stt[T_PINE].leaf_tid);

	if (draw_all_pine) {
		vbomgr.render_all(GL_QUADS);
	}
	else {
		assert(!low_detail);

		for (const_iterator i = begin(); i != end(); ++i) {
			int const type(i->get_type());
			if (type != T_PINE && type != T_SH_PINE) continue;
			i->draw(2, shadow_only, (size() < 100), vbomgr, xlate); // only cull pine tree leaves if there aren't too many
		}
	}
	vbomgr.end_render();
}


void small_tree_group::draw_non_pine_leaves(bool shadow_only, vector3d const xlate) const {

	enable_blend();

	for (const_iterator i = begin(); i != end(); ++i) {
		int const type(i->get_type());
		if (i->is_pine_tree()) continue;
		if (i == begin() || (i-1)->get_type() != type) {select_texture(stt[type].leaf_tid);} // first of this type
		i->draw(2, shadow_only, (size() < 100), vbo_manager[0], xlate); // only cull pine tree leaves if there aren't too many
	}
	disable_blend();
}


float calc_tree_scale() {return (Z_SCENE_SIZE*tree_scale)/16.0;}
float calc_tree_size () {return SM_TREE_SIZE*Z_SCENE_SIZE/calc_tree_scale();}


void small_tree_group::gen_trees(int x1, int y1, int x2, int y2) {

	if (vegetation == 0.0) return;
	generated = 1;
	float const tscale(calc_tree_scale()), tsize(calc_tree_size()); // random tree generation based on transformed mesh height function
	int const ntrees(int(min(1.0f, vegetation*tscale*tscale/8.0f)*NUM_SMALL_TREES));
	if (ntrees == 0) return;
	assert(x1 < x2 && y1 < y2);
	float const x0(TREE_DIST_MH_S*X_SCENE_SIZE + xoff2*DX_VAL), y0(TREE_DIST_MH_S*Y_SCENE_SIZE + yoff2*DY_VAL);
	float const tds(TREE_DIST_SCALE*(XY_MULT_SIZE/16384.0));
	int const tree_prob(max(1, XY_MULT_SIZE/ntrees)), trees_per_block(max(1, ntrees/XY_MULT_SIZE)), skip_val(max(1, int(1.0/sqrt(tree_scale))));
	float const xv(DX_VAL*(tds*(get_xval(x1)+x0) - (MESH_X_SIZE >> 1))), yv(DY_VAL*(tds*(get_yval(y1)+y0) - (MESH_Y_SIZE >> 1)));
	mesh_xy_grid_cache_t height_gen;
	height_gen.build_arrays(xv, yv, tds*DX_VAL*DX_VAL, tds*DY_VAL*DY_VAL, (x2-x1), (y2-y1));
	//PRINT_TIME("Delete");
	
	for (int i = y1; i < y2; i += skip_val) {
		for (int j = x1; j < x2; j += skip_val) {
			set_rand2_state((657435*(i + yoff2) + 243543*(j + xoff2) + 734533*rand_gen_index),
							(845631*(j + xoff2) + 667239*(i + yoff2) + 846357*rand_gen_index));
			if ((rand2_seed_mix()%tree_prob) != 0) continue; // not selected
			//float const val(eval_one_surface_point((tds*(xpos+x0)-xoff2), (tds*(ypos+y0)-yoff2)));
			float const val(height_gen.eval_index(j-x1, i-y1, 1));
			float const dist_test(get_rel_height(val, -zmax_est, zmax_est));

			for (int n = 0; n < trees_per_block; ++n) {
				if (dist_test > (SM_TREE_AMT*(1.0 - TREE_DIST_RAND) + TREE_DIST_RAND*rand_float2())) continue; // tree density function test
				rand2_mix();
				float const xpos(get_xval(j) + 0.5*skip_val*DX_VAL*signed_rand_float2());
				float const ypos(get_yval(i) + 0.5*skip_val*DY_VAL*signed_rand_float2());
				float const zpos(interpolate_mesh_zval(xpos, ypos, 0.0, 1, 1));
				int const ttype(get_tree_type_from_height(zpos));
				if (ttype == TREE_NONE) continue;
				float const height(rand_uniform2(0.4*tsize, tsize)), width(rand_uniform2(0.25*height, 0.35*height));
				small_tree st(point(xpos, ypos, zpos-0.1*height), height, width, ttype, 0);
				st.setup_rotation();
				add_tree(st);
			}
		}
	}
	sort_by_type();
}


int get_tree_class_from_height(float zpos) {

	if (island) {
		if (zpos > 0.14*Z_SCENE_SIZE || zpos < max(water_plane_z, (zmin + OCEAN_DEPTH))) return TREE_CLASS_NONE; // none
	}
	else {
		if (zpos < water_plane_z) return TREE_CLASS_NONE;
	}
	if (force_tree_class >= 0) {
		assert(force_tree_class < NUM_TREE_CLASSES);
		return force_tree_class;
	}
	if (world_mode == WMODE_INF_TERRAIN) return TREE_CLASS_PINE;
	float const relh(get_rel_height(zpos, -zmax_est, zmax_est));
	if (relh > 0.6) return TREE_CLASS_PINE;

	if (island) {
		if (zpos < 0.85*(zmin + OCEAN_DEPTH)) return TREE_CLASS_PALM;
	}
	else {
		if (zpos < 0.85*water_plane_z && relh < 0.435) return TREE_CLASS_PALM;
	}
	return TREE_CLASS_DECID;
}


int get_tree_type_from_height(float zpos) {

	//return T_DECID; // TESTING
	switch (get_tree_class_from_height(zpos)) {
	case TREE_CLASS_NONE: return TREE_NONE;
	case TREE_CLASS_PINE: return ((rand2()%10 == 0) ? T_SH_PINE : T_PINE);
	case TREE_CLASS_PALM: return T_PALM;
	case TREE_CLASS_DECID:
		//if (tree_mode == 3) return TREE_NONE; // use a large (complex) tree here
		return T_DECID + rand2()%3;
	default: assert(0);
	}
	return TREE_NONE; // never gets here
}


int add_small_tree(point const &pos, float height, float width, int tree_type, bool calc_z) {

	assert(height > 0.0 && width > 0.0);
	small_trees.add_tree(small_tree(pos, height, width, (abs(tree_type)%NUM_ST_TYPES), calc_z)); // could have a type error
	small_trees.back().calc_points(small_trees.vbo_manager[0], 0, 0);
	return 1; // might return zero in some case
}


colorRGBA colorgen(float r1, float r2, float g1, float g2, float b1, float b2) {
	return colorRGBA(rand_uniform2(r1, r2), rand_uniform2(g1, g2), rand_uniform2(b1, b2), 1.0);
}


void gen_small_trees() {

	if (num_trees == 0) return;
	//RESET_TIME;
	
	if (SMALL_TREE_COLL && !small_trees.empty()) {
		remove_small_tree_cobjs();
		purge_coll_freed(1);
	}
	small_trees.clear_all();
	small_trees = small_tree_group(); // really force a clear
	small_trees.gen_trees(1, 1, MESH_X_SIZE-1, MESH_Y_SIZE-1);
	small_trees.finalize(0, 0);
	//PRINT_TIME("Gen");
	small_trees.add_cobjs();
	//PRINT_TIME("Cobj");
	cout << "small trees: " << small_trees.size() << endl;
}


void clear_sm_tree_vbos() {
	small_trees.clear_vbos();
}

void add_small_tree_coll_objs() { // doesn't handle rotation angle
	small_trees.add_cobjs();
}

void remove_small_tree_cobjs() {
	small_trees.remove_cobjs();
}

void shift_small_trees(vector3d const &vd) {
	if (num_trees > 0) return; // dynamically created, not placed
	small_trees.translate_by(vd);
}


void draw_small_trees(bool shadow_only) {

	//RESET_TIME;
	if (small_trees.empty() || !(tree_mode & 2)) return;
	if (small_trees.size() < 100) {small_trees.sort_by_dist_to_camera();} // shadow_only?
	shader_t s;
	bool const can_use_shaders(!disable_shaders), v(!shadow_only);
	colorRGBA orig_fog_color;

	if (can_use_shaders) {
		orig_fog_color = setup_smoke_shaders(s, 0.0, 0, 0, 0, v, v, 0, 0, v, (USE_BUMP_MAP && v), 0, v); // dynamic lights, but no smoke
		s.add_uniform_float("tex_scale_t", 5.0);

		if (USE_BUMP_MAP) {
			vector4d const tangent(0.0, 0.0, 1.0, 1.0); // FIXME: set based on tree trunk direction?
			int const tangent_loc(s.get_attrib_loc("tangent"));
			if (tangent_loc >= 0) glVertexAttrib4fv(tangent_loc, &tangent.x);
			set_multitex(5);
			select_texture(BARK2_NORMAL_TEX, 0);
			set_multitex(0);
		}
	}
	set_fill_mode();
	small_trees.draw_branches(shadow_only);
	if (s.is_setup()) end_smoke_shaders(s, orig_fog_color);
	small_trees.vbo_manager[0].upload();

	if (can_use_shaders) {
		s.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
		orig_fog_color = setup_smoke_shaders(s, 0.75, 3, 0, 0, v, v, 0, 0, v, 0, 0, v, v); // dynamic lights, but no smoke
		s.add_uniform_float("base_color_scale", 0.0); // hack to force usage of material properties instead of color
	}
	set_lighted_sides(2);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.75);
	small_trees.draw_pine_leaves(shadow_only);
	
	if (s.is_setup()) {
		s.add_uniform_float("base_color_scale", 1.0);
		end_smoke_shaders(s, orig_fog_color);
	}
	small_trees.draw_non_pine_leaves(shadow_only); // not using shaders
	glDisable(GL_ALPHA_TEST);
	set_lighted_sides(1);
	glDisable(GL_TEXTURE_2D);
	tree_scenery_pld.draw_and_clear();
	//PRINT_TIME("small tree draw");
}


small_tree::small_tree(point const &p, float h, float w, int t, bool calc_z) :
	type(t), height(h), width(w), r_angle(0.0), rx(0.0), ry(0.0), pos(p)
{
	clear_vbo_mgr_ix();
	for (unsigned i = 0; i < 3; ++i) rv[i] = rand2d();
	if (calc_z) pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1) - 0.1*height;

	switch (type) {
	case T_PINE: // pine tree
		width  *= 1.1;
		height *= 1.2;
		color   = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35);
		break;
	case T_SH_PINE: // pine tree
		width  *= 1.2;
		height *= 0.8;
		color   = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35);
		break;
	case T_PALM: // palm tree
		color   = colorgen(0.6, 1.0, 0.6, 1.0, 0.6, 1.0);
		width  *= 1.4;
		height *= 2.0;
		break;
	case T_DECID: // decidious tree
		color = colorgen(0.6, 0.9, 0.7, 1.0, 0.4, 0.7);
		break;
	case T_TDECID: // tall decidious tree
		color = colorgen(0.3, 0.6, 0.7, 1.0, 0.1, 0.2);
		break;
	case T_BUSH: // bush
		color  = colorgen(0.6, 0.9, 0.7, 1.0, 0.1, 0.2);
		pos.z += 0.3*height;
		if (rand2()%100 < 50) pos.z -= height*rand_uniform2(0.0, 0.2);
		break;
	default: assert(0);
	}
	color.alpha = 1.0;
}


void small_tree::setup_rotation() {

	if (!(rand2()&3) && (type == T_DECID || type == T_TDECID || type == T_PALM)) {
		r_angle = ((type == T_PALM) ? rand_uniform2(-15.0, 15.0) : rand_uniform2(-7.5, 7.5));
		rx = rand_uniform2(-1.0, 1.0);
		ry = rand_uniform2(-1.0, 1.0);
	}
}


vector3d small_tree::get_rot_dir() const {

	vector3d dir(plus_z);
	if (r_angle != 0.0) rotate_vector3d(vector3d(rx, ry, 0.0), -r_angle/TO_DEG, dir); // oops, rotation is backwards
	return dir;
}


void small_tree::add_cobjs(cobj_params &cp, cobj_params &cp_trunk) {

	if (!is_over_mesh(pos)) return; // not sure why, but this makes drawing slower
	vector3d const dir(get_rot_dir());

	if (!DRAW_COBJS) {
		cp.tid       = stt[type].leaf_tid;
		cp_trunk.tid = stt[type].bark_tid;
	}
	if (type != T_BUSH && type != T_SH_PINE) {
		float const hval((type == T_PINE) ? 1.0 : 0.75), zb(pos.z - 0.2*width), zbot(get_tree_z_bottom(zb, pos)), len(hval*height + (zb - zbot));
		point const p1((pos + dir*(zbot - pos.z)));
		coll_id.push_back(add_coll_cylinder(p1, (p1 + dir*len), stt[type].ws*width, stt[type].w2*width, cp_trunk, -1, 1));
	}
	switch (type) {
	case T_PINE: // pine tree
	case T_SH_PINE: // short pine tree
		// Note: smaller than the actual pine tree render at the base so that it's easier for the player to walk under pine trees
		coll_id.push_back(add_coll_cylinder((pos + ((type == T_PINE) ? dir*(0.35*height) : all_zeros)), (pos + dir*height), get_pine_tree_radius(), 0.0, cp, -1, 1));
		break;
	case T_DECID: // decidious tree
		coll_id.push_back(add_coll_sphere((pos + dir*(0.75*height)), 1.3*width, cp, -1, 1));
		break;
	case T_TDECID: // tall decidious tree
		coll_id.push_back(add_coll_sphere((pos + dir*(0.8*height)), 0.8*width, cp, -1, 1));
		break;
	case T_BUSH: // bush
		coll_id.push_back(add_coll_sphere(pos, 1.2*width, cp, -1, 1));
		break;
	case T_PALM: // palm tree
		coll_id.push_back(add_coll_sphere((pos + dir*(0.65*height)), 0.5*width, cp, -1, 1)); // small sphere
		break;
	default: assert(0);
	}
}


void small_tree::remove_cobjs() {

	for (unsigned j = 0; j < coll_id.size(); ++j) {
		remove_reset_coll_obj(coll_id[j]);
	}
}


float small_tree::get_pine_tree_radius() const {
	float const height0(((type == T_PINE) ? 0.75 : 1.0)*height);
	return 0.35*(height0 + 0.03/tree_scale);
}


void small_tree::calc_points(vbo_vnc_block_manager_t &vbo_manager, bool low_detail, bool pri_dim) {

	if (type != T_PINE && type != T_SH_PINE) return; // only for pine trees
	float const height0(((type == T_PINE) ? 0.75 : 1.0)*height), sz_scale(SQRT2*get_pine_tree_radius());
	point const center(pos + point(0.0, 0.0, ((type == T_PINE) ? 0.35*height : 0.0)));
	vector<vert_norm> &points(vbo_manager.temp_points);

	if (!low_detail) { // high detail
		points.resize(4*N_PT_LEVELS*N_PT_RINGS);
		float const rd(0.5), theta0((int(1.0E6*height0)%360)*TO_RADIANS);

		for (unsigned j = 0, ix = 0; j < N_PT_LEVELS; ++j) {
			float const sz(sz_scale*(N_PT_LEVELS - j - 0.4)/(float)N_PT_LEVELS);
			float const z((j + 1.8)*height0/(N_PT_LEVELS + 2.8) - rd*sz);
			vector3d const scale(sz, sz, sz);

			for (unsigned k = 0; k < N_PT_RINGS; ++k) {
				float const theta(TWO_PI*(3.3*j + k/(float)N_PT_RINGS) + theta0);
				add_rotated_quad_pts(&points.front(), ix, theta, rd, z, center, scale); // bounds are (sz, sz, rd*sz+z)
			}
		}
		vbo_mgr_ix = vbo_manager.add_points_with_offset(points, color);
	}
	else { // low detail billboard
		float const nz_avg(0.816), r1(1.7*sz_scale), z1(center.z - 0.55*sz_scale - 0.3*height), z2(center.z + 1.45*sz_scale - 0.1*height);
		unsigned const d(pri_dim);
		unsigned ix(0);
		points.resize(4);
		vert_norm vn;
		vn.n[!d] = 0.0;
		vn.n[d]  = -sqrt(1.0 - nz_avg*nz_avg); // partially facing up and partially facing towards the camera
		vn.n[2]  = nz_avg;
		vn.v[!d] = pos[!d];
		vn.v[2]  = z2;
		vn.v[d]  = pos[d] + r1;
		points[ix++] = vn;
		vn.v[d]  = pos[d] - r1;
		points[ix++] = vn;
		vn.v[2]  = z1;
		points[ix++] = vn;
		vn.v[d]  = pos[d] + r1;
		points[ix++] = vn;
		vbo_manager.add_points(points, color);
	}
}


float small_tree::get_zmax() const {return (pos.z + height);} // approximate


void small_tree::add_trunk_as_line(vector<point> &points) const {

	float const hval(is_pine_tree() ? 1.0 : 0.75);
	points.push_back(pos);
	points.push_back(pos + get_rot_dir()*(hval*height*((stt[type].w2 == 0.0) ? 0.7 : 1.0))); // slightly shorter for distant pine trees
}


colorRGBA small_tree::get_bark_color() const {

	colorRGBA tcolor(stt[type].c);
	UNROLL_3X(tcolor[i_] = min(1.0f, tcolor[i_]*(0.85f + 0.3f*rv[i_]));)
	return tcolor;
}


void small_tree::draw(int mode, bool shadow_only, bool do_cull, vbo_vnc_block_manager_t const &vbo_manager,
	vector3d const xlate, vector<point> *points) const
{
	if (!(tree_mode & 2)) return; // disabled
	if (type == T_BUSH && !(mode & 2)) return; // no bark
	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	bool const pine_tree(is_pine_tree()), cobj_cull(pine_tree && do_cull);
	if (shadow_only ? !is_over_mesh(pos2) : !sphere_in_camera_view(pos2, max(1.5*width, 0.5*height), (cobj_cull ? 2 : 0))) return;
	float const dist(distance_to_camera(pos + xlate));
	float const zoom_f(do_zoom ? ZOOM_FACTOR : 1.0), size(zoom_f*SM_TREE_QUALITY*stt[type].ss*width*window_width/dist);
	int const max_sides(N_CYL_SIDES);

	if ((mode & 1) && (shadow_only || size > 1.0) && type != T_BUSH) { // trunk
		float const cz(camera_origin.z - cview_radius*cview_dir.z), vxy(1.0 - (cz - pos.z)/dist);

		if (pine_tree || dist < 0.2 || vxy >= 0.2*width/stt[type].h) { // if trunk not obscured by leaves
			float const hval(pine_tree ? 1.0 : 0.75), w1(stt[type].ws*width), w2(stt[type].w2*width);
			float const zb(pos.z - 0.2*width), zbot(get_tree_z_bottom(zb, pos)), len(hval*height + (zb - zbot));
			vector3d const dir(get_rot_dir());
			point const p1((pos + dir*(zbot - pos.z)));

			if (shadow_only) {
				cylinder_3dw const cylin(p1, (p1 + dir*len), w1, w2);
				draw_cylin_quad_proj(cylin, ((cylin.p1 + cylin.p2)*0.5 - get_camera_pos()), -1, 1);
			}
			else if (LINE_THRESH*zoom_f*(w1 + w2) < dist) { // draw as line
				point const p2(p1 + dir*(len*((w2 == 0.0) ? 0.8 : 1.0)));
				
				if (points) {
					points->push_back(p1 + xlate);
					points->push_back(p2 + xlate);
				}
				else {
					tree_scenery_pld.add_textured_line(p1+xlate, p2+xlate, get_bark_color(), stt[type].bark_tid);
				}
			}
			else { // draw as cylinder
				if (world_mode == WMODE_GROUND) {set_color(get_bark_color());}
				select_texture(stt[type].bark_tid);
				glPushMatrix();
				translate_to(pos);
				if (r_angle != 0.0) glRotatef(r_angle, rx, ry, 0.0);
				glTranslatef(0.0, 0.0, (zbot - pos.z));
				int const nsides2(max(3, min(2*max_sides/3, int(0.25*size))));
				draw_cylin_fast(w1, w2, len, nsides2, 1, 0, 1); // trunk (draw quad if small?)
				glPopMatrix();
			}
		}
	}
	if (mode & 2) { // leaves
		if (pine_tree) { // 30 quads per tree
			assert(vbo_mgr_ix >= 0);
			vbo_manager.render_range(GL_QUADS, vbo_mgr_ix, vbo_mgr_ix+1); // draw textured quad if far away?
		}
		else { // palm or decidious
			set_color(color);
			glPushMatrix();
			translate_to(pos);
			if (r_angle != 0.0) glRotatef(r_angle, rx, ry, 0.0);

			switch (type) { // draw leaves
			case T_DECID: // decidious tree
				glTranslatef(0.0, 0.0, 0.75*height);
				glScalef(1.2, 1.2, 0.8);
				break;
			case T_TDECID: // tall decidious tree
				glTranslatef(0.0, 0.0, 1.0*height);
				glScalef(0.7, 0.7, 1.6);
				break;
			case T_BUSH: // bush
				glScalef((0.1*height+0.8*width)/width, (0.1*height+0.8*width)/width, 1.0);
				break;
			case T_PALM: // palm tree
				glTranslatef(0.0, 0.0, 0.71*height-0.5*width);
				glScalef(1.2, 1.2, 0.5);
				break;
			}
			int const nsides(max(6, min(max_sides, (shadow_only ? get_smap_ndiv(width) : (int)size))));

			/*if (type == T_BUSH && nsides >= 24) {
				draw_cube_map_sphere(all_zeros, width, nsides/2, 1, 1); // slower, but looks better
			}
			else*/ {
				draw_sphere_dlist(all_zeros, width, nsides, 1, (type == T_PALM));
			}
			glPopMatrix();
		} // end pine else
	} // end mode
}



