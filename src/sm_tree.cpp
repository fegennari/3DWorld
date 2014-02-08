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
small_tree_group tree_instances;
pt_line_drawer tree_scenery_pld;

extern int window_width, shadow_detail, draw_model, island, num_trees, do_zoom, tree_mode, xoff2, yoff2;
extern int rand_gen_index, display_mode, force_tree_class;
extern unsigned max_unique_trees;
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
	sm_tree_type(0.13, 0.15, 0.75, 0.8, TREE_C,  TREE_HEMI_TEX, BARK3_TEX), // T_DECID // HEDGE_TEX?
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


unsigned get_pine_tree_inst_gpu_mem() {return tree_instances.get_gpu_mem();}


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

	if (!low_detail) { // high detail, create using multiple threads
		#pragma omp parallel for schedule(static,1)
		for (int i = 0; i < (int)size(); ++i) {
			operator[](i).calc_points(vbo_manager[low_detail], low_detail);
		}
	}
	else {
		for (iterator i = begin(); i != end(); ++i) {
			i->calc_points(vbo_manager[low_detail], low_detail);
		}
	}
}


void small_tree_group::finalize_upload_and_clear_pts(bool low_detail) {

	if (empty() || is_uploaded(low_detail)) return;

	if (instanced && !low_detail) {
		tree_instances.finalize_upload_and_clear_pts(0); // high detail
		return;
	}
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
		coll |= i->check_sphere_coll(center, radius); // Note: calculates/updates center, so we can't early terminate in case there are multiple collisions
	}
	return coll;
}


// if t is non-NULL, calculate the point of intersection
bool small_tree_group::line_intersect(point const &p1, point const &p2, float *t) const { // if t != NULL, it should start at 1.0 (the end of the line)

	bool coll(0);

	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->line_intersect(p1, p2, t)) {
			coll = 1;
			if (t == NULL) break;
		}
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
		i->draw(1, shadow_only, xlate, points);
	}
}


void small_tree_group::get_back_to_front_ordering(vector<pair<float, unsigned> > &to_draw, vector3d const &xlate) const {

	point const ref_pos(get_camera_pos() - xlate);

	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->is_visible_pine(xlate)) {to_draw.push_back(make_pair(p2p_dist_sq(i->get_pos(), ref_pos), i-begin()));}
	}
	sort(to_draw.begin(), to_draw.end()); // sort front to back for early Z culling
}


void small_tree_group::draw_pine_leaves(bool shadow_only, bool low_detail, bool draw_all_pine, bool sort_front_to_back, vector3d const &xlate, int xlate_loc) {

	if (empty()) return;
	select_texture((draw_model != 0) ? WHITE_TEX : (low_detail ? PINE_TREE_TEX : stt[T_PINE].leaf_tid));

	if (instanced && !low_detail) { // sort_front_to_back not supported
		assert(xlate_loc >= 0);

		if (!draw_all_pine || insts.size() != size()) { // recompute insts
			insts.clear();

			for (const_iterator i = begin(); i != end(); ++i) {
				if (draw_all_pine || i->is_visible_pine(xlate)) {insts.push_back(pine_tree_inst_t(i->get_inst_id(), i->get_pos()));}
			}
			if (insts.empty()) return; // nothing to draw
			sort(insts.begin(), insts.end());
		}
		static vector<point> pts; // class member?
		pts.resize(insts.size()); // force max possible size so that pts won't be resized later, and the pointer should remain valid
		glVertexAttribPointer(xlate_loc, 3, GL_FLOAT, GL_FALSE, sizeof(point), &pts.front());
		vbo_vnc_block_manager_t const &vbomgr(tree_instances.vbo_manager[0]); // high detail
		vbomgr.begin_render(1);

		for (vector<pine_tree_inst_t>::const_iterator i = insts.begin(); i != insts.end();) { // Note: no increment
			unsigned const inst_id(i->id);
			unsigned ix(0);
			for (; i != insts.end() && i->id == inst_id; ++i) {pts[ix++] = i->pt;}
			assert(inst_id < tree_instances.size());
			tree_instances[inst_id].draw_pine(vbomgr, ix);
		}
		vbomgr.end_render();
	}
	else {
		vbo_vnc_block_manager_t const &vbomgr(vbo_manager[low_detail]);
		vbomgr.begin_render(1);

		if (draw_all_pine) {
			vbomgr.render_all(GL_QUADS);
		}
		else if (sort_front_to_back) {
			assert(!low_detail);
			vector<pair<float, unsigned> > to_draw;
			get_back_to_front_ordering(to_draw, xlate);

			for (unsigned i = 0; i < to_draw.size(); ++i) {
				operator[](to_draw[i].second).draw_pine_leaves(vbomgr, xlate);
			}
		}
		else {
			for (const_iterator i = begin(); i != end(); ++i) {
				i->draw_pine_leaves(vbomgr, xlate);
			}
		}
		vbomgr.end_render();
	}
}


void small_tree_group::draw_non_pine_leaves(bool shadow_only, vector3d const &xlate) const {

	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->is_pine_tree()) continue;
		assert(!instanced); // only pine trees can be instanced
		int const type(i->get_type());
		if (i == begin() || (i-1)->get_type() != type) {select_texture(stt[type].leaf_tid);} // first of this type
		i->draw(2, shadow_only, xlate);
	}
}


float calc_tree_scale() {return (Z_SCENE_SIZE*tree_scale)/16.0;}
float calc_tree_size () {return SM_TREE_SIZE*Z_SCENE_SIZE/calc_tree_scale();}


void create_pine_tree_instances() {

	if (tree_instances.size() == max_unique_trees) return; // already done
	assert(tree_instances.empty());
	rand_gen_t rgen;

	for (unsigned i = 0; i < max_unique_trees; ++i) { // Note: pine trees not rotated
		int const ttype((rgen.rand()%10 == 0) ? T_SH_PINE : T_PINE);
		float const height(rgen.rand_uniform(0.4, 1.0)), width(rgen.rand_uniform(0.25, 0.35));
		tree_instances.add_tree(small_tree(all_zeros, height, width, ttype, 0, rgen));
	}
}


// density = x1,y1 x2,y1 x1,y2 x2,y2
void small_tree_group::gen_trees(int x1, int y1, int x2, int y2, float const density[4]) {

	generated = 1; // mark as generated if we got here, even if there are no actual trees generated
	if (sm_tree_density == 0.0 || vegetation == 0.0) return;
	if (density[0] == 0.0 && density[1] == 0.0 && density[2] == 0.0 && density[3] == 0.0) return;
	assert(x1 < x2 && y1 < y2);
	float const tscale(calc_tree_scale()), tsize(calc_tree_size()), ntrees_mult(vegetation*sm_tree_density*tscale*tscale/8.0f);
	int const skip_val(max(1, int(1.0/(sqrt(sm_tree_density*tree_scale)))));
	bool const use_hmap_tex(using_tiled_terrain_hmap_tex());
	bool const approx_zval(world_mode == WMODE_INF_TERRAIN && !use_hmap_tex); // faster, but lower z-value accuracy, and only works for tiled terrain mode
	float const tds(TREE_DIST_SCALE*(XY_MULT_SIZE/16384.0)), xscale(tds*DX_VAL*DX_VAL), yscale(tds*DY_VAL*DY_VAL);
	float const zval_adj((world_mode == WMODE_INF_TERRAIN) ? 0.0 : -0.1);
	mesh_xy_grid_cache_t density_gen, height_gen; // random tree generation based on transformed mesh height function
	density_gen.build_arrays(xscale*(x1 + xoff2), yscale*(y1 + yoff2), xscale, yscale, (x2-x1), (y2-y1));
	if (approx_zval) {height_gen.build_arrays(DX_VAL*(x1 + xoff2 - (MESH_X_SIZE >> 1) + 0.5), DY_VAL*(y1 + yoff2 - (MESH_Y_SIZE >> 1) + 0.5), DX_VAL, DY_VAL, (x2-x1), (y2-y1));}
	float const dxv(skip_val/(x2 - x1 - 1.0)), dyv(skip_val/(y2 - y1 - 1.0));
	float xv(0.0), yv(0.0);

	for (int i = y1; i < y2; i += skip_val, yv += dyv) {
		for (int j = x1; j < x2; j += skip_val, xv += dxv) {
			float const cur_density(yv*(xv*density[3] + (1.0-xv)*density[2]) + (1.0-yv)*(xv*density[1] + (1.0-xv)*density[0]));
			int const ntrees(int(min(1.0f, cur_density*ntrees_mult)*NUM_SMALL_TREES));
			if (ntrees == 0) continue;
			int const tree_prob(max(1, XY_MULT_SIZE/ntrees));
			rgen.set_state((657435*(i + yoff2) + 243543*(j + xoff2) + 734533*rand_gen_index),
						   (845631*(j + xoff2) + 667239*(i + yoff2) + 846357*rand_gen_index));
			rgen.rand(); // increase randomness (slower, but less regular - which is important for pine trees in tiled terrain mode)
			if ((rgen.rand_seed_mix()%tree_prob) != 0) continue; // not selected
			float const dist_test(get_rel_height(density_gen.eval_index(j-x1, i-y1, 1), -zmax_est, zmax_est));
			int const trees_per_block(max(1, ntrees/XY_MULT_SIZE));

			for (int n = 0; n < trees_per_block; ++n) {
				if (dist_test > (TREE_DEN_THRESH*(1.0 - TREE_DIST_RAND) + TREE_DIST_RAND*rgen.rand_float())) continue; // tree density function test
				rgen.rand_mix();
				float const xval(get_xval(j) + 0.5*skip_val*DX_VAL*rgen.signed_rand_float());
				float const yval(get_yval(i) + 0.5*skip_val*DY_VAL*rgen.signed_rand_float());
				float const zpos(approx_zval ? height_gen.eval_index(j-x1, i-y1, 1) : interpolate_mesh_zval(xval, yval, 0.0, 1, 1));
				int const ttype(get_tree_type_from_height(zpos, rgen));
				if (ttype == TREE_NONE) continue;
				if (use_hmap_tex && get_tiled_terrain_height_tex_norm(j+xoff2, i+yoff2).z < 0.8) {continue;}

				if (instanced) {
					assert(ttype == T_SH_PINE || ttype == T_PINE);
					assert(zval_adj == 0.0); // must be inf terrain mode
					assert(max_unique_trees > 0);
					add_tree(small_tree(point(xval, yval, zpos), (rgen.rand()%max_unique_trees)));
				}
				else {
					float const height(tsize*rgen.rand_uniform(0.4, 1.0)), width(height*rgen.rand_uniform(0.25, 0.35));
					small_tree st(point(xval, yval, zpos+zval_adj*height), height, width, ttype, 0, rgen);
					st.setup_rotation(rgen);
					add_tree(st);
				}
			}
		}
		xv = 0.0; // reset for next y iter
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
	float const density[4] = {1.0, 1.0, 1.0, 1.0};
	small_trees.clear_all();
	small_trees = small_tree_group(); // really force a clear
	small_trees.gen_trees(1, 1, MESH_X_SIZE-1, MESH_Y_SIZE-1, density);
	small_trees.finalize(0);
	//PRINT_TIME("Gen");
	small_trees.add_cobjs();
	//PRINT_TIME("Cobj");
	cout << "small trees: " << small_trees.size() << endl;
}


void clear_sm_tree_vbos() {

	small_trees.clear_vbos();
	tree_instances.clear_vbos();
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


// instanced constructor
small_tree::small_tree(point const &p, unsigned instance_id) {

	*this = tree_instances.get_tree(instance_id);
	assert(coll_id.empty());
	clear_vbo_mgr_ix();
	inst_id = instance_id;
	pos     = p;
	float const tsize(calc_tree_size());
	width  *= tsize;
	height *= tsize;
}


small_tree::small_tree(point const &p, float h, float w, int t, bool calc_z, rand_gen_t &rgen) :
	type(t), inst_id(-1), height(h), width(w), r_angle(0.0), rx(0.0), ry(0.0), pos(p)
{
	clear_vbo_mgr_ix();
	for (unsigned i = 0; i < 3; ++i) {rv[i] = rgen.randd();} // for bark color
	if (calc_z) {pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1) - 0.1*height;}

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
		if (rgen.rand()%100 < 50) {pos.z -= height*rgen.rand_uniform(0.0, 0.2);}
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

	if (r_angle == 0.0) return plus_z;
	vector3d dir(plus_z);
	rotate_vector3d(vector3d(rx, ry, 0.0), -r_angle/TO_DEG, dir); // oops, rotation is backwards
	return dir;
}


cylinder_3dw small_tree::get_trunk_cylin() const {

	if (type == T_BUSH) {return cylinder_3dw();}
	vector3d const dir(get_rot_dir());
	bool const is_pine(is_pine_tree());
	float const hval(is_pine ? 1.0 : 0.75), zb(pos.z - 0.2*width), zbot(get_tree_z_bottom(zb, pos)), len(hval*height + (zb - zbot));
	float const mod_width(width*(is_pine ? 0.8*len/(hval*height) : 1.0));
	point const p1((pos + dir*(zbot - pos.z)));
	return cylinder_3dw(p1, (p1 + dir*len), stt[type].ws*mod_width, stt[type].w2*mod_width);
}


void small_tree::add_cobjs(cobj_params &cp, cobj_params &cp_trunk) {

	if (!is_over_mesh(pos)) return; // not sure why, but this makes drawing slower

	if (!DRAW_COBJS) {
		cp.tid       = stt[type].leaf_tid;
		cp_trunk.tid = stt[type].bark_tid;
	}
	if (type != T_BUSH && type != T_SH_PINE) {coll_id.push_back(add_coll_cylinder(get_trunk_cylin(), cp_trunk, -1, 1));}
	vector3d const dirh(get_rot_dir()*height);

	switch (type) {
	case T_PINE: // pine tree
	case T_SH_PINE: // short pine tree
		// Note: smaller than the actual pine tree render at the base so that it's easier for the player to walk under pine trees
		coll_id.push_back(add_coll_cylinder((pos + ((type == T_PINE) ? 0.35*dirh : all_zeros)), (pos + dirh), get_pine_tree_radius(), 0.0, cp, -1, 1));
		break;
	case T_DECID: // decidious tree
		coll_id.push_back(add_coll_sphere((pos + 0.75*dirh), 1.3*width, cp, -1, 1));
		break;
	case T_TDECID: // tall decidious tree
		coll_id.push_back(add_coll_sphere((pos + 0.8*dirh), 0.8*width, cp, -1, 1));
		break;
	case T_BUSH: // bush
		coll_id.push_back(add_coll_sphere(pos, 1.2*width, cp, -1, 1));
		break;
	case T_PALM: // palm tree
		coll_id.push_back(add_coll_sphere((pos + 0.65*dirh), 0.5*width, cp, -1, 1)); // small sphere
		break;
	default: assert(0);
	}
}


void small_tree::remove_cobjs() {

	for (unsigned j = 0; j < coll_id.size(); ++j) {
		remove_reset_coll_obj(coll_id[j]);
	}
}


// very simple check against trunk only, for collisions with a player walking on the ground
bool small_tree::check_sphere_coll(point &center, float radius) const {

	if (type == T_BUSH) return 0; // no trunk, not yet handled
	return sphere_vert_cylin_intersect(center, radius, get_trunk_cylin());
}


bool small_tree::line_intersect(point const &p1, point const &p2, float *t) const {

	assert(is_pine_tree()); // Note: can work on other tree types, but it's more complex, and we only need pine trees in tiled terrain mode
	vector3d const dirh(get_rot_dir()*height);
	cylinder_3dw const cylins[2] = {get_trunk_cylin(), cylinder_3dw((pos + ((type == T_PINE) ? 0.35*dirh : all_zeros)), (pos + dirh), get_pine_tree_radius(), 0.0)};
	bool coll(0);
	
	for (unsigned i = 0; i < 2; ++i) {
		if (t == NULL) {
			coll |= line_intersect_cylinder(p1, p2, cylins[i], 0); // no check ends
			if (coll) break;
		}
		else {
			float t_new(0.0);
			
			if (line_intersect_trunc_cone(p1, p2, cylins[i].p1, cylins[i].p2, cylins[i].r1, cylins[i].r2, 0, t_new) && t_new < *t) {
				*t   = t_new; // would be more efficient if we could pass *t in as t_max
				coll = 1;
			}
		}
	}
	return coll;
}


float small_tree::get_pine_tree_radius() const {

	float const height0(((type == T_PINE) ? 0.75 : 1.0)*height);
	return 0.35*(height0 + 0.03/tree_scale);
}


void small_tree::calc_points(vbo_vnc_block_manager_t &vbo_manager, bool low_detail, bool update_mode) {

	if (type != T_PINE && type != T_SH_PINE) return; // only for pine trees
	float const sz_scale(SQRT2*get_pine_tree_radius()), dz((type == T_PINE) ? 0.35*height : 0.0);

	if (!low_detail) { // high detail
		unsigned const npts(4*N_PT_LEVELS*N_PT_RINGS);
		float const rd(0.5), height0(((type == T_PINE) ? 0.75 : 1.0)*height), theta0((int(1.0E6*height0)%360)*TO_RADIANS);
		point const center(pos + point(0.0, 0.0, dz));
		vert_norm points[npts];

		for (unsigned j = 0, ix = 0; j < N_PT_LEVELS; ++j) {
			float const sz(sz_scale*(N_PT_LEVELS - j - 0.4)/(float)N_PT_LEVELS);
			float const z((j + 1.8)*height0/(N_PT_LEVELS + 2.8) - rd*sz);

			for (unsigned k = 0; k < N_PT_RINGS; ++k) {
				float const theta(TWO_PI*(3.3*j + k/(float)N_PT_RINGS) + theta0);
				add_rotated_quad_pts(points, ix, theta, z, center, sz, rd*sz); // bounds are (sz, sz, rd*sz+z)
			}
		}
		if (update_mode) {
			assert(vbo_mgr_ix >= 0);
			vbo_manager.update_range(points, npts, color, vbo_mgr_ix, vbo_mgr_ix+1);
		}
		else { // we only get into this case when running in parallel
			#pragma omp critical(pine_tree_vbo_update)
			vbo_mgr_ix = vbo_manager.add_points_with_offset(points, npts, color);
		}
	}
	else { // low detail billboard
		assert(!update_mode);
		vert_norm points[4];
		vert_norm vn(pos, vector3d(1.5*sz_scale/calc_tree_size(), 0.0, 0.0)); // ranges from around 0.25 to 0.75
		vn.v.z = pos.z + dz + 1.45*sz_scale + 0.1*height;
		points[0] = points[1] = vn; // top two vertices
		vn.v.z = pos.z + dz - 0.55*sz_scale - 0.2*height;
		points[2] = points[3] = vn; // bottom two vertices
		vbo_manager.add_points(points, 4, color);
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


void small_tree::draw_pine(vbo_vnc_block_manager_t const &vbo_manager, unsigned num_instances) const { // 30 quads per tree

	assert(is_pine_tree());
	assert(vbo_mgr_ix >= 0);
	vbo_manager.render_range(GL_QUADS, vbo_mgr_ix, vbo_mgr_ix+1, num_instances);
}


bool small_tree::is_visible_pine(vector3d const &xlate) const {

	return camera_pdu.sphere_visible_test((pos + xlate + point(0.0, 0.0, 0.5*height)), max(1.5*width, 0.5*height));
}


void small_tree::draw_pine_leaves(vbo_vnc_block_manager_t const &vbo_manager, vector3d const &xlate) const {

	if (is_pine_tree() && is_visible_pine(xlate)) {draw_pine(vbo_manager);}
}


void small_tree::draw(int mode, bool shadow_only, vector3d const &xlate, vector<point> *points) const {

	if (!(tree_mode & 2)) return; // disabled
	if (type == T_BUSH && !(mode & 2)) return; // no bark
	if (shadow_only ? !is_over_mesh(pos + xlate + point(0.0, 0.0, 0.5*height)) : !is_visible_pine(xlate)) return;
	float const zoom_f(do_zoom ? ZOOM_FACTOR : 1.0), size_scale(zoom_f*stt[type].ss*width*window_width);

	if ((mode & 1) && type != T_BUSH) { // trunk
		float const dist(distance_to_camera(pos + xlate));

		if (shadow_only || size_scale > dist) {
			if (is_pine_tree() || dist < 0.2 || (1.0 - ((camera_origin.z - cview_radius*cview_dir.z) - pos.z)/dist)*stt[type].h >= 0.2*width) { // if trunk not obscured by leaves
				cylinder_3dw const cylin(get_trunk_cylin()); // cache in the tree?

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
					int const nsides2(max(3, min(N_CYL_SIDES, int(0.25*size_scale/dist))));
					draw_fast_cylinder(cylin.p1, cylin.p2, cylin.r1, cylin.r2, nsides2, 1);
				}
			}
		}
	}
	if (mode & 2) { // leaves
		assert(!is_pine_tree()); // handled through draw_pine_leaves()
		// palm or decidious
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
	} // end mode
}

