// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 2/13/03

#include "3DWorld.h"
#include "mesh.h"
#include "small_tree.h"
#include "gl_ext_arb.h"
#include "shaders.h"


float const SM_TREE_SIZE    = 0.05;
float const TREE_DIST_RAND  = 0.2;
float const LINE_THRESH     = 700.0;
bool const SMALL_TREE_COLL  = 1;
bool const DRAW_COBJS       = 0; // for debugging
bool const USE_BUMP_MAP     = 1;
int  const NUM_SMALL_TREES  = 40000;
unsigned const N_PT_LEVELS  = 6;
unsigned const N_PT_RINGS   = 5;


small_tree_group small_trees;
pt_line_drawer tree_scenery_pld;

extern bool has_dir_lights;
extern int window_width, shadow_detail, draw_model, island, num_trees, do_zoom, tree_mode, xoff2, yoff2;
extern int rand_gen_index, display_mode, force_tree_class;
extern float zmin, zmax_est, water_plane_z, tree_scale, sm_tree_density, vegetation, OCEAN_DEPTH;


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


void small_tree_group::finalize(bool low_detail) {

	if (empty()) return;
	assert(!is_uploaded(low_detail));
	vbo_manager[low_detail].clear();
	vbo_manager[low_detail].reserve_pts(4*(low_detail ? 1 : N_PT_LEVELS*N_PT_RINGS)*num_pine_trees);

	for (iterator i = begin(); i != end(); ++i) {
		i->calc_points(vbo_manager[low_detail], low_detail);
	}
}


void small_tree_group::finalize_upload_and_clear_pts(bool low_detail) {

	if (empty() || is_uploaded(low_detail)) return;
	//RESET_TIME;
	finalize(low_detail);
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


void small_tree_group::add_cobjs_range(iterator b, iterator e) {

	if (b == e || !SMALL_TREE_COLL) return;
	assert(b < e && e <= end());
	cobj_params cp      (0.65, GREEN, DRAW_COBJS, 0, NULL, 0, -1);
	cobj_params cp_trunk(0.9, TREE_C, DRAW_COBJS, 0, NULL, 0, -1);
	cp.shadow       = (shadow_detail >= 5);
	cp_trunk.shadow = (shadow_detail >= 6);

	for (iterator i = b; i < e; ++i) {
		i->add_cobjs(cp, cp_trunk);
	}
}


void small_tree_group::remove_cobjs() {

	for (iterator i = begin(); i != end(); ++i) {
		i->remove_cobjs();
	}
}


bool small_tree_group::check_sphere_coll(point &center, float radius) const {

	bool coll(0);

	for (const_iterator i = begin(); i != end(); ++i) {
		coll |= i->check_sphere_coll(center, radius);
	}
	return coll;
}


void small_tree_group::translate_by(vector3d const &vd) {

	for (iterator i = begin(); i != end(); ++i) {
		i->translate_by(vd);
	}
}


void small_tree_group::draw_branches(bool shadow_only, vector3d const &xlate, vector<point> *points) const {

	BLACK.do_glColor();

	for (const_iterator i = begin(); i != end(); ++i) {
		i->draw(1, shadow_only, vbo_manager[0], xlate, points);
	}
}


void small_tree_group::draw_pine_leaves(bool shadow_only, bool low_detail, bool draw_all_pine, vector3d const &xlate) const {

	vbo_vnc_block_manager_t const &vbomgr(vbo_manager[low_detail]);
	vbomgr.begin_render(1);
	select_texture((draw_model != 0) ? WHITE_TEX : (low_detail ? PINE_TREE_TEX : stt[T_PINE].leaf_tid));

	if (draw_all_pine) {
		vbomgr.render_all(GL_QUADS);
	}
	else {
		assert(!low_detail);

		for (const_iterator i = begin(); i != end(); ++i) {
			if (i->is_pine_tree()) {i->draw(2, shadow_only, vbomgr, xlate);}
		}
	}
	vbomgr.end_render();
}


void small_tree_group::draw_non_pine_leaves(bool shadow_only, vector3d const &xlate) const {

	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->is_pine_tree()) continue;
		int const type(i->get_type());
		if (i == begin() || (i-1)->get_type() != type) {select_texture(stt[type].leaf_tid);} // first of this type
		i->draw(2, shadow_only, vbo_manager[0], xlate);
	}
}


float calc_tree_scale() {return (Z_SCENE_SIZE*tree_scale)/16.0;}
float calc_tree_size () {return SM_TREE_SIZE*Z_SCENE_SIZE/calc_tree_scale();}


void small_tree_group::gen_trees(int x1, int y1, int x2, int y2, float vegetation_) {

	generated = 1; // mark as generated if we got here, even if there are no actual trees generated
	if (vegetation_ == 0.0) return;
	float const tscale(sm_tree_density*calc_tree_scale()), tsize(calc_tree_size()); // random tree generation based on transformed mesh height function
	int const ntrees(int(min(1.0f, vegetation_*tscale*tscale/8.0f)*NUM_SMALL_TREES));
	if (ntrees == 0) return;
	assert(x1 < x2 && y1 < y2);
	int const tree_prob(max(1, XY_MULT_SIZE/ntrees)), trees_per_block(max(1, ntrees/XY_MULT_SIZE)), skip_val(max(1, int(1.0/(sm_tree_density*sqrt(tree_scale)))));
	float const tds(TREE_DIST_SCALE*(XY_MULT_SIZE/16384.0)), xscale(tds*DX_VAL*DX_VAL), yscale(tds*DY_VAL*DY_VAL);
	mesh_xy_grid_cache_t density_gen;
	density_gen.build_arrays(xscale*(x1 + xoff2), yscale*(y1 + yoff2), xscale, yscale, (x2-x1), (y2-y1));
	
	for (int i = y1; i < y2; i += skip_val) {
		for (int j = x1; j < x2; j += skip_val) {
			rgen.set_state((657435*(i + yoff2) + 243543*(j + xoff2) + 734533*rand_gen_index),
						   (845631*(j + xoff2) + 667239*(i + yoff2) + 846357*rand_gen_index));
			if ((rgen.rand_seed_mix()%tree_prob) != 0) continue; // not selected
			float const dist_test(get_rel_height(density_gen.eval_index(j-x1, i-y1, 1), -zmax_est, zmax_est));

			for (int n = 0; n < trees_per_block; ++n) {
				if (dist_test > (TREE_DEN_THRESH*(1.0 - TREE_DIST_RAND) + TREE_DIST_RAND*rgen.rand_float())) continue; // tree density function test
				rgen.rand_mix();
				float const xpos(get_xval(j) + 0.5*skip_val*DX_VAL*rgen.signed_rand_float());
				float const ypos(get_yval(i) + 0.5*skip_val*DY_VAL*rgen.signed_rand_float());
				float const zpos(interpolate_mesh_zval(xpos, ypos, 0.0, 1, 1));
				int const ttype(get_tree_type_from_height(zpos, rgen));
				if (ttype == TREE_NONE) continue;
				float const height(tsize*rgen.rand_uniform(0.4, 1.0)), width(height*rgen.rand_uniform(0.25, 0.35));
				small_tree st(point(xpos, ypos, zpos-0.1*height), height, width, ttype, 0, rgen);
				st.setup_rotation(rgen);
				add_tree(st);
			}
		}
	}
	if (world_mode == WMODE_GROUND) {sort_by_type();}
}


bool small_tree_group::update_zvals(int x1, int y1, int x2, int y2) {

	assert(trunk_pts.empty());
	bool updated(0);

	for (iterator i = begin(); i != end(); ++i) {
		point const &pos(i->get_pos());
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
		if (xpos < x1 || xpos > x2 || ypos < y1 || ypos > y2) continue;
		float const new_z(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1) - 0.1*i->get_height()); // use-real_equation=0
		if (fabs(pos.z - new_z) < 0.01*i->get_height()) continue;
		i->remove_cobjs();
		i->translate_by(point(0.0, 0.0, (new_z - pos.z)));
		i->update_points_vbo(vbo_manager[0], 0); // high detail only
		if (tree_mode & 2) {add_cobjs_range(i, i+1);}
		updated = 1;
	}
	return updated;
}


void small_tree_group::update_zmax(float &tzmax) const {
	for (const_iterator i = begin(); i != end(); ++i) {tzmax = max(tzmax, i->get_zmax());}
}


int get_tree_class_from_height(float zpos, bool pine_trees_only) {

	if (zpos < water_plane_z) return TREE_CLASS_NONE;
	if (island && (zpos > 0.14*Z_SCENE_SIZE || zpos < (zmin + OCEAN_DEPTH))) return TREE_CLASS_NONE;
	
	if (force_tree_class >= 0) {
		assert(force_tree_class < NUM_TREE_CLASSES);
		return force_tree_class;
	}
	float const relh(get_rel_height(zpos, -zmax_est, zmax_est));
	if (relh > 0.9) return TREE_CLASS_NONE; // too high
	if (relh > 0.6) return TREE_CLASS_PINE;
	if (pine_trees_only) {return ((tree_mode == 3) ? TREE_CLASS_NONE : TREE_CLASS_PINE);}

	if (island) {
		if (zpos < 0.85*(zmin + OCEAN_DEPTH)) return TREE_CLASS_PALM;
	}
	else {
		if (zpos < 0.85*water_plane_z && relh < 0.435) return TREE_CLASS_PALM;
	}
	return TREE_CLASS_DECID;
}


int get_tree_type_from_height(float zpos, rand_gen_t &rgen) {

	//return T_DECID; // TESTING
	switch (get_tree_class_from_height(zpos, (world_mode == WMODE_INF_TERRAIN))) {
	case TREE_CLASS_NONE: return TREE_NONE;
	case TREE_CLASS_PINE: return ((rgen.rand()%10 == 0) ? T_SH_PINE : T_PINE);
	case TREE_CLASS_PALM: return T_PALM;
	case TREE_CLASS_DECID:
		//if (tree_mode == 3) return TREE_NONE; // use a large (complex) tree here
		return T_DECID + rgen.rand()%3;
	default: assert(0);
	}
	return TREE_NONE; // never gets here
}


int add_small_tree(point const &pos, float height, float width, int tree_type, bool calc_z) {

	assert(height > 0.0 && width > 0.0);
	small_trees.add_tree(small_tree(pos, height, width, (abs(tree_type)%NUM_ST_TYPES), calc_z, small_trees.rgen)); // could have a type error
	small_trees.back().calc_points(small_trees.vbo_manager[0], 0);
	return 1; // might return zero in some case
}


colorRGBA colorgen(float r1, float r2, float g1, float g2, float b1, float b2, rand_gen_t &rgen) {
	return colorRGBA(rgen.rand_uniform(r1, r2), rgen.rand_uniform(g1, g2), rgen.rand_uniform(b1, b2), 1.0);
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
	small_trees.gen_trees(1, 1, MESH_X_SIZE-1, MESH_Y_SIZE-1, vegetation);
	small_trees.finalize(0);
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

bool update_small_tree_zvals(int x1, int y1, int x2, int y2) {
	return small_trees.update_zvals(x1, y1, x2, y2);
}


void draw_small_trees(bool shadow_only) {

	//RESET_TIME;
	if (small_trees.empty() || !(tree_mode & 2)) return;
	if (small_trees.size() < 100) {small_trees.sort_by_dist_to_camera();} // shadow_only?
	shader_t s;
	bool const v(!shadow_only), use_bump_map(USE_BUMP_MAP && v);

	// draw trunks
	setup_smoke_shaders(s, 0.0, 0, 0, 0, v, v, 0, 0, v, use_bump_map, 0, v); // dynamic lights, but no smoke
	s.add_uniform_float("tex_scale_t", 5.0);

	if (use_bump_map) {
		set_active_texture(5);
		select_texture(BARK2_NORMAL_TEX, 0);
		set_active_texture(0);
	}
	set_fill_mode();
	small_trees.draw_branches(shadow_only);
	s.add_uniform_float("tex_scale_t", 1.0);
	s.end_shader();

	// draw leaves
	small_trees.vbo_manager[0].upload();
	setup_smoke_shaders(s, 0.75, 3, 0, 0, v, v, 0, 0, v, 0, 0, v, v, 1); // dynamic lights, but no smoke, use light colors
	set_lighted_sides(2);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.75);
	small_trees.draw_pine_leaves(shadow_only);
	s.end_shader();
	small_trees.draw_non_pine_leaves(shadow_only); // not using shaders
	glDisable(GL_ALPHA_TEST);
	set_lighted_sides(1);
	glDisable(GL_TEXTURE_2D);
	tree_scenery_pld.draw_and_clear();
	//PRINT_TIME("small tree draw");
}


small_tree::small_tree(point const &p, float h, float w, int t, bool calc_z, rand_gen_t &rgen) :
	type(t), height(h), width(w), r_angle(0.0), rx(0.0), ry(0.0), pos(p)
{
	clear_vbo_mgr_ix();
	for (unsigned i = 0; i < 3; ++i) rv[i] = rgen.randd();
	if (calc_z) pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1) - 0.1*height;

	switch (type) {
	case T_PINE: // pine tree
		width  *= 1.1;
		height *= 1.2;
		color   = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35, rgen);
		break;
	case T_SH_PINE: // pine tree
		width  *= 1.2;
		height *= 0.8;
		color   = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35, rgen);
		break;
	case T_PALM: // palm tree
		color   = colorgen(0.6, 1.0, 0.6, 1.0, 0.6, 1.0, rgen);
		width  *= 1.4;
		height *= 2.0;
		break;
	case T_DECID: // decidious tree
		color = colorgen(0.6, 0.9, 0.7, 1.0, 0.4, 0.7, rgen);
		break;
	case T_TDECID: // tall decidious tree
		color = colorgen(0.3, 0.6, 0.7, 1.0, 0.1, 0.2, rgen);
		break;
	case T_BUSH: // bush
		color  = colorgen(0.6, 0.9, 0.7, 1.0, 0.1, 0.2, rgen);
		pos.z += 0.3*height;
		if (rgen.rand()%100 < 50) pos.z -= height*rgen.rand_uniform(0.0, 0.2);
		break;
	default: assert(0);
	}
	color.alpha = 1.0;
}


void small_tree::setup_rotation(rand_gen_t &rgen) {

	if (!(rgen.rand()&3) && (type == T_DECID || type == T_TDECID || type == T_PALM)) {
		r_angle = ((type == T_PALM) ? rgen.rand_uniform(-15.0, 15.0) : rgen.rand_uniform(-7.5, 7.5));
		rx = rgen.rand_uniform(-1.0, 1.0);
		ry = rgen.rand_uniform(-1.0, 1.0);
	}
}


vector3d small_tree::get_rot_dir() const {

	vector3d dir(plus_z);
	if (r_angle != 0.0) rotate_vector3d(vector3d(rx, ry, 0.0), -r_angle/TO_DEG, dir); // oops, rotation is backwards
	return dir;
}


cylinder_3dw small_tree::get_trunk_cylin() const {

	if (type == T_BUSH) {return cylinder_3dw();}
	vector3d const dir(get_rot_dir());
	float const hval(is_pine_tree() ? 1.0 : 0.75), zb(pos.z - 0.2*width), zbot(get_tree_z_bottom(zb, pos)), len(hval*height + (zb - zbot));
	point const p1((pos + dir*(zbot - pos.z)));
	return cylinder_3dw(p1, (p1 + dir*len), stt[type].ws*width, stt[type].w2*width);
}


void small_tree::add_cobjs(cobj_params &cp, cobj_params &cp_trunk) {

	if (!is_over_mesh(pos)) return; // not sure why, but this makes drawing slower

	if (!DRAW_COBJS) {
		cp.tid       = stt[type].leaf_tid;
		cp_trunk.tid = stt[type].bark_tid;
	}
	if (type != T_BUSH && type != T_SH_PINE) {coll_id.push_back(add_coll_cylinder(get_trunk_cylin(), cp_trunk, -1, 1));}
	vector3d const dir(get_rot_dir());

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


bool small_tree::check_sphere_coll(point &center, float radius) const { // very simple check against trunk only

	if (type == T_BUSH) return 0; // no trunk, not yet handled
	return sphere_vert_cylin_intersect(center, radius, get_trunk_cylin());
}


float small_tree::get_pine_tree_radius() const {
	float const height0(((type == T_PINE) ? 0.75 : 1.0)*height);
	return 0.35*(height0 + 0.03/tree_scale);
}


void small_tree::calc_points(vbo_vnc_block_manager_t &vbo_manager, bool low_detail, bool update_mode) {

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
		if (update_mode) {
			assert(vbo_mgr_ix >= 0);
			vbo_manager.update_range(points, color, vbo_mgr_ix, vbo_mgr_ix+1);
		}
		else {
			vbo_mgr_ix = vbo_manager.add_points_with_offset(points, color);
		}
	}
	else { // low detail billboard
		assert(!update_mode);
		points.resize(4);
		vert_norm vn(pos, vector3d(1.5*sz_scale/calc_tree_size(), 0.0, 0.816));
		vn.v.z = center.z + 1.45*sz_scale + 0.1*height;
		points[0] = points[1] = vn; // top two vertices
		vn.v.z = center.z - 0.55*sz_scale - 0.2*height;
		points[2] = points[3] = vn; // bottom two vertices
		vbo_manager.add_points(points, color);
	}
}


void small_tree::update_points_vbo(vbo_vnc_block_manager_t &vbo_manager, bool low_detail) {
	if (vbo_mgr_ix >= 0) {calc_points(vbo_manager, low_detail, 1);}
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


void small_tree::draw_pine(vbo_vnc_block_manager_t const &vbo_manager) const { // 30 quads per tree

	assert(is_pine_tree());
	assert(vbo_mgr_ix >= 0);
	vbo_manager.render_range(GL_QUADS, vbo_mgr_ix, vbo_mgr_ix+1); // draw textured quad if far away?
}


void small_tree::draw(int mode, bool shadow_only, vbo_vnc_block_manager_t const &vbo_manager, vector3d const &xlate, vector<point> *points) const {

	if (!(tree_mode & 2)) return; // disabled
	if (type == T_BUSH && !(mode & 2)) return; // no bark
	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	if (shadow_only ? !is_over_mesh(pos2) : !sphere_in_camera_view(pos2, max(1.5*width, 0.5*height), 0)) return;
	float const zoom_f(do_zoom ? ZOOM_FACTOR : 1.0), size_scale(zoom_f*stt[type].ss*width*window_width);
	bool const pine_tree(is_pine_tree());

	if ((mode & 1) && type != T_BUSH) {
		float const dist(distance_to_camera(pos + xlate)), size(size_scale/dist);

		if (shadow_only || size > 1.0) { // trunk
			float const cz(camera_origin.z - cview_radius*cview_dir.z), vxy(1.0 - (cz - pos.z)/dist);

			if (pine_tree || dist < 0.2 || vxy >= 0.2*width/stt[type].h) { // if trunk not obscured by leaves
				cylinder_3dw const cylin(get_trunk_cylin());

				if (!shadow_only && LINE_THRESH*zoom_f*(cylin.r1 + cylin.r2) < dist) { // draw as line
					point const p2((cylin.r2 == 0.0) ? (0.2*cylin.p1 + 0.8*cylin.p2) : cylin.p2);
				
					if (points) {
						points->push_back(cylin.p1 + xlate);
						points->push_back(p2 + xlate);
					}
					else {
						tree_scenery_pld.add_textured_line(cylin.p1+xlate, p2+xlate, get_bark_color(), stt[type].bark_tid);
					}
				}
				else { // draw as cylinder
					if (world_mode == WMODE_GROUND) {set_color(get_bark_color());}
					if (!shadow_only) {select_texture(stt[type].bark_tid);}
					int const nsides2(max(3, min(N_CYL_SIDES, int(0.25*size))));
					draw_fast_cylinder(cylin.p1, cylin.p2, cylin.r1, cylin.r2, nsides2, 1);
				}
			}
		}
	}
	if (mode & 2) { // leaves
		if (pine_tree) {
			draw_pine(vbo_manager);
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
			int const nsides(max(6, min(N_SPHERE_DIV, (shadow_only ? get_smap_ndiv(width) : (int)(size_scale/distance_to_camera(pos + xlate))))));

			/*if (type == T_BUSH && nsides >= 24) {
				draw_cube_map_sphere(all_zeros, width, N_SPHERE_DIV/2, 1); // slower, but looks better
			}
			else*/ {
				draw_sphere_vbo(all_zeros, width, nsides, 1, (type == T_PALM));
			}
			glPopMatrix();
		} // end pine else
	} // end mode
}



