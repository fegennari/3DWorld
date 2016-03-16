// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 2/13/03

#include "3DWorld.h"
#include "mesh.h"
#include "small_tree.h"
#include "gl_ext_arb.h"
#include "shaders.h"


float const SM_TREE_SIZE    = 0.05;
float const TREE_DIST_RAND  = 0.125;
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

extern bool tree_indir_lighting, only_pine_palm_trees;
extern int window_width, draw_model, num_trees, do_zoom, tree_mode, xoff2, yoff2;
extern int rand_gen_index, display_mode, force_tree_class;
extern unsigned max_unique_trees;
extern float zmin, zmax_est, water_plane_z, tree_scale, sm_tree_density, vegetation, tree_density_thresh, tree_height_scale, CAMERA_RADIUS;


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
	sm_tree_type(0.13, 0.15, 0.75, 0.7, TREE_C,  HEDGE_TEX,     BARK1_TEX), // T_TDECID
	sm_tree_type(0.00, 0.15, 0.00, 0.8, TREE_C,  HEDGE_TEX,     BARK1_TEX), // T_BUSH NOTE: bark texture is not used in trees, but is used in logs
	sm_tree_type(0.03, 0.15, 1.00, 0.6, TREE_C,  PALM_FROND_TEX,PALM_BARK_TEX), // T_PALM
	sm_tree_type(0.00, 0.07, 0.00, 0.4, PTREE_C, PINE_TEX,      BARK2_TEX), // T_SH_PINE
};

bool is_pine_tree_type(int type) {return (type == T_PINE || type == T_SH_PINE);}

int get_bark_tex_for_tree_type(int type) {
	assert(type < NUM_ST_TYPES);
	return stt[type].bark_tid;
}

colorRGBA get_tree_trunk_color(int type, bool modulate_with_texture) {
	int const bark_tid(get_bark_tex_for_tree_type(type));
	colorRGBA c(stt[type].c);
	if (modulate_with_texture && bark_tid >= 0) {return c.modulate_with(texture_color(bark_tid));}
	return c;
}


unsigned get_pine_tree_inst_gpu_mem() {return tree_instances.get_gpu_mem();}


void small_tree_group::add_tree(small_tree const &st) {

	if (st.is_pine_tree()) {
		max_pt_radius = max(max_pt_radius, st.get_pine_tree_radius());
		++num_pine_trees;
	}
	else if (st.get_type() == T_PALM) {++num_palm_trees;}
	push_back(st);
}


void small_tree_group::calc_trunk_pts() {

	if (!trunk_pts.empty()) return; // already calculated
	trunk_pts.reserve(2*size());
	for (iterator i = begin(); i != end(); ++i) {i->add_trunk_as_line(trunk_pts);}
}


void small_tree_group::finalize(bool low_detail) {

	if (empty()) return;
	assert(!is_uploaded(low_detail));
	vbo_manager[low_detail].clear();
	vbo_manager[low_detail].reserve_pts(4*(low_detail ? 1 : N_PT_LEVELS*N_PT_RINGS)*num_pine_trees);

	#pragma omp parallel for schedule(static,1) if (!low_detail)
	for (int i = 0; i < (int)size(); ++i) {
		operator[](i).calc_points(vbo_manager[low_detail], low_detail);
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


void small_tree_group::add_trunk_pts(point const &xlate, vector<vert_wrap_t> &pts) const {
	for (vector<point>::const_iterator i = trunk_pts.begin(); i != trunk_pts.end(); ++i) {pts.push_back(*i + xlate);}
}


void small_tree_group::clear_vbos() {
	for (unsigned i = 0; i < 2; ++i) {vbo_manager[i].clear_vbo();}
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
	for (iterator i = b; i < e; ++i) {i->add_cobjs(cp, cp_trunk);}
}

void small_tree_group::remove_cobjs() {
	for (iterator i = begin(); i != end(); ++i) {i->remove_cobjs();}
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
	for (iterator i = begin(); i != end(); ++i) {i->translate_by(vd);}
}


void small_tree_group::draw_trunks(bool shadow_only, vector3d const &xlate, vector<vert_wrap_t> *points) const {

	static vector<vert_norm_tc> cylin_verts; // class member?
	for (const_iterator i = begin(); i != end(); ++i) {i->draw_trunks(shadow_only, xlate, points, &cylin_verts);}

	if (!cylin_verts.empty()) {
		if (!shadow_only) {
			get_tree_trunk_color(T_PINE, 0).set_for_cur_shader(); // all a constant color
			select_texture(get_bark_tex_for_tree_type(T_PINE));
		}
		// make this work for non-pine trees (where we have to deal with colors, texture changes, and tris vs. strips)?
		draw_and_clear_verts(cylin_verts, GL_TRIANGLES); // inf terrain mode pine tree trunks
	}
}

void small_tree_group::sort_by_dist_to_camera() {

	point const camera(get_camera_pos());
	if (dist_less_than(camera, last_cpos, CAMERA_RADIUS)) return; // no sort needed
	last_cpos = camera;
	sort(begin(), end(), small_tree::comp_by_type_dist(camera));
}

void small_tree_group::get_back_to_front_ordering(vector<pair<float, unsigned> > &to_draw, vector3d const &xlate) const { // for leaves

	point const ref_pos(get_camera_pos() - xlate);

	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->are_leaves_visible(xlate)) {to_draw.push_back(make_pair(p2p_dist_sq(i->get_pos(), ref_pos), i-begin()));}
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
				if (draw_all_pine || i->are_leaves_visible(xlate)) {insts.push_back(pine_tree_inst_t(i->get_inst_id(), i->get_pos()));}
			}
			if (insts.empty()) return; // nothing to draw
			sort(insts.begin(), insts.end());
		}
		vbo_vnc_block_manager_t const &vbomgr(tree_instances.vbo_manager[0]); // high detail
		vbomgr.begin_render();
		inst_pts.resize(insts.size());
		for (vector<pine_tree_inst_t>::const_iterator i = insts.begin(); i != insts.end(); ++i) {inst_pts[i-insts.begin()] = i->pt;}
		void const *vbo_ptr(get_dynamic_vbo_ptr(&inst_pts.front(), inst_pts.size()*sizeof(point)));
		unsigned ptr_offset(0), ix(0);

		for (vector<pine_tree_inst_t>::const_iterator i = insts.begin(); i != insts.end();) { // Note: no increment
			unsigned const inst_id(i->id);
			assert(inst_id < tree_instances.size());
			for (ix = 0; i != insts.end() && i->id == inst_id; ++i, ++ix) {}
			glVertexAttribPointer(xlate_loc, 3, GL_FLOAT, GL_FALSE, sizeof(point), (void const *)((point const *)vbo_ptr + ptr_offset));
			tree_instances[inst_id].draw_pine(vbomgr, ix);
			ptr_offset += ix;
		}
		vbomgr.end_render();
	}
	else {
		vbo_vnc_block_manager_t const &vbomgr(vbo_manager[low_detail]);
		vbomgr.begin_render();

		if (draw_all_pine) {
			vbomgr.render_all();
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
			for (const_iterator i = begin(); i != end(); ++i) {i->draw_pine_leaves(vbomgr, xlate);}
		}
		vbomgr.end_render();
	}
}


void small_tree_group::draw_non_pine_leaves(bool shadow_only, bool draw_palm, bool draw_non_palm, int xlate_loc, int scale_loc, vector3d const &xlate) const {

	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->is_pine_tree()) continue;
		assert(!instanced); // only pine trees can be instanced
		int const type(i->get_type());
		if (!((type == T_PALM) ? draw_palm : draw_non_palm)) continue;
		if (i == begin() || (i-1)->get_type() != type) {select_texture(stt[type].leaf_tid);} // first of this type
		i->draw_leaves(shadow_only, xlate_loc, scale_loc, xlate);
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
	density_gen.build_arrays(xscale*(x1 + xoff2), yscale*(y1 + yoff2), xscale, yscale, (x2-x1), (y2-y1), 0, 1); // force_sine_mode=1
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
			int const trees_per_block(max(1, ntrees/XY_MULT_SIZE));
			float const hval(density_gen.eval_index(j-x1, i-y1, 0));

			for (int n = 0; n < trees_per_block; ++n) {
				if (hval > get_median_height(tree_density_thresh - TREE_DIST_RAND*rgen.rand_float())) continue; // tree density function test
				rgen.rand_mix();
				float const xval(get_xval(j) + 0.5*skip_val*DX_VAL*rgen.signed_rand_float());
				float const yval(get_yval(i) + 0.5*skip_val*DY_VAL*rgen.signed_rand_float());
				float const zpos(approx_zval ? height_gen.eval_index(j-x1, i-y1, 1) : interpolate_mesh_zval(xval, yval, 0.0, 1, 1));
				int const ttype(get_tree_type_from_height(zpos, rgen));
				if (ttype == TREE_NONE) continue;
				if (use_hmap_tex && get_tiled_terrain_height_tex_norm(j+xoff2, i+yoff2).z < 0.8) continue;

				if (instanced) { // check voxel terrain? or does this mode only get used for tiled terrain?
					assert(ttype == T_SH_PINE || ttype == T_PINE);
					assert(zval_adj == 0.0); // must be inf terrain mode
					assert(max_unique_trees > 0);
					add_tree(small_tree(point(xval, yval, zpos), (rgen.rand()%max_unique_trees)));
				}
				else {
					float const height(tsize*rgen.rand_uniform(0.4, 1.0)), width(height*rgen.rand_uniform(0.25, 0.35));
					point const pos(xval, yval, zpos+zval_adj*height);
					if (point_inside_voxel_terrain(pos)) continue; // don't create trees that start inside voxels (but what about trees that grow into voxels?)
					add_tree(small_tree(pos, height, width, ttype, 0, rgen, 1)); // allow_rotation=1
				}
			} // for n
		} // for j
		xv = 0.0; // reset for next y iter
	} // for i
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
	
	if (force_tree_class >= 0) {
		assert(force_tree_class < NUM_TREE_CLASSES);
		return force_tree_class;
	}
	float const relh(get_rel_height(zpos, -zmax_est, zmax_est));
	if (relh > 0.9) return TREE_CLASS_NONE; // too high
	if (relh > 0.6) return TREE_CLASS_PINE;
	//if (zpos < 0.85*water_plane_z && relh < 0.435) return TREE_CLASS_PALM;
	if (pine_trees_only) {return ((tree_mode == 3) ? TREE_CLASS_NONE : TREE_CLASS_PINE);}
	if (zpos < 0.85*water_plane_z && relh < 0.435) return TREE_CLASS_PALM;
	return (only_pine_palm_trees ? TREE_CLASS_PINE : TREE_CLASS_DECID);
}


int get_tree_type_from_height(float zpos, rand_gen_t &rgen) {

	//return T_DECID; // TESTING
	switch (get_tree_class_from_height(zpos, (world_mode == WMODE_INF_TERRAIN && (tree_mode & 2)))) { // force pine trees in small tree mode
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
	if (!shadow_only && small_trees.size() < 100) {small_trees.sort_by_dist_to_camera();} // not in shadow pass, since trees usually don't overlap in z
	shader_t s;
	bool const all_pine(small_trees.num_pine_trees == small_trees.size());
	bool const v(!shadow_only), use_bump_map(USE_BUMP_MAP && !shadow_only && all_pine); // bump maps only work with pine tree trunks

	// draw trunks
	if (shadow_only) {
		s.begin_color_only_shader();
	}
	else {
		setup_smoke_shaders(s, 0.0, 0, 0, tree_indir_lighting, 1, 1, 0, 0, 1, use_bump_map, 0, 1, 0, 0.0, 0.0, 0, 0, 1); // dynamic lights, but no smoke, is_outside=1
		s.add_uniform_float("tex_scale_t", 5.0);
	}
	if (use_bump_map) {select_multitex(BARK2_NORMAL_TEX, 5, 1);}
	small_trees.draw_trunks(shadow_only);
	if (!shadow_only) {s.add_uniform_float("tex_scale_t", 1.0);}
	s.end_shader();

	// draw leaves
	float const wind_mag(get_plant_leaf_wind_mag(shadow_only));

	if (small_trees.num_pine_trees > 0) { // pine trees
		small_trees.vbo_manager[0].upload();
		if (wind_mag > 0.0) {s.set_prefix("#define ENABLE_WIND", 0);} // VS
		s.set_prefix("#define NO_SPECULAR", 1); // FS - disable rain effect
		setup_smoke_shaders(s, 0.5, 3, 0, (v && tree_indir_lighting), v, v, 0, 0, v, 0, 0, v, v, 0.0, 0.0, 0, 0, 1); // dynamic lights, but no smoke, texgen, is_outside=1
		setup_leaf_wind(s, wind_mag, 0);
		small_trees.draw_pine_leaves(shadow_only);
		s.end_shader();
	}
	if (!all_pine) { // non-pine trees
		if (small_trees.num_pine_trees + small_trees.num_palm_trees < small_trees.size()) { // deciduous trees
			s.begin_simple_textured_shader(0.75, !shadow_only); // with lighting, unless shadow_only
			bind_draw_sphere_vbo(1, 1); // texture, even in shadow pass, to handle alpha masking of some tree types
			int const xlate_loc(s.get_uniform_loc("xlate")), scale_loc(s.get_uniform_loc("scale"));
			small_trees.draw_non_pine_leaves(shadow_only, 0, 1, xlate_loc, scale_loc);
			s.set_uniform_vector3d(xlate_loc, zero_vector);
			s.set_uniform_vector3d(scale_loc, vector3d(1,1,1));
			bind_vbo(0);
			s.end_shader();
		}
		if (small_trees.num_palm_trees > 0) { // palm trees
			if (wind_mag > 0.0) {s.set_prefix("#define ENABLE_WIND", 0);} // VS
			setup_smoke_shaders(s, 0.75, 0, 0, 0, v, v, 0, 0, v, 0, 0, 0, 1); // dynamic lights, but no smoke (slow, but looks better)
			setup_leaf_wind(s, wind_mag, 0);
			small_trees.draw_non_pine_leaves(shadow_only, 1, 0);
			s.end_shader();
		}
	}
	if (!tree_scenery_pld.empty()) {
		shader_t s;
		s.begin_untextured_lit_glcolor_shader();
		tree_scenery_pld.draw_and_clear();
		s.end_shader();
	}
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
	trunk_cylin = get_trunk_cylin();
}


small_tree::small_tree(point const &p, float h, float w, int t, bool calc_z, rand_gen_t &rgen, bool allow_rotation) :
	type(t), inst_id(-1), height(h), width(w), r_angle(0.0), rx(0.0), ry(0.0), pos(p)
{
	height *= tree_height_scale;
	clear_vbo_mgr_ix();
	bark_color = stt[type].c;
	if (!is_pine_tree()) {UNROLL_3X(bark_color[i_] = min(1.0, bark_color[i_]*(0.85 + 0.3f*rgen.randd()));)} // gen bark color for decid trees
	if (calc_z) {pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1) - 0.1*height;}

	switch (type) {
	case T_PINE: // pine tree
		width  *= 1.1;
		height *= 1.2;
		color   = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35, rgen);
		break;
	case T_SH_PINE: // short pine tree
		width  *= 1.2;
		height *= 0.8;
		color   = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35, rgen);
		break;
	case T_PALM: // palm tree
		color   = colorgen(0.9, 1.0, 0.9, 1.0, 0.9, 1.0, rgen);
		width  *= 1.4;
		height *= 2.0;
		break;
	case T_DECID: // decidious tree
		color = colorgen(0.5, 0.8, 0.7, 1.0, 0.4, 0.7, rgen);
		break;
	case T_TDECID: // tall decidious tree
		color = colorgen(0.4, 0.7, 0.8, 1.0, 0.3, 0.5, rgen);
		break;
	case T_BUSH: // bush
		color  = colorgen(0.5, 0.8, 0.7, 1.0, 0.3, 0.5, rgen);
		pos.z += 0.3*height;
		if (rgen.rand()%100 < 50) {pos.z -= height*rgen.rand_uniform(0.0, 0.2);}
		break;
	default: assert(0);
	}
	color.alpha = 1.0;
	if (allow_rotation) {setup_rotation(rgen);}
	trunk_cylin = get_trunk_cylin();
}


void small_tree::setup_rotation(rand_gen_t &rgen) {

	if (!(rgen.rand()&3) && (type == T_DECID || type == T_TDECID || type == T_PALM)) {
		r_angle = ((type == T_PALM) ? rgen.rand_uniform(-15.0, 15.0) : rgen.rand_uniform(-7.5, 7.5));
		rx = rgen.rand_uniform(-1.0, 1.0);
		ry = rgen.rand_uniform(-1.0, 1.0);
	}
}

vector3d small_tree::get_rot_dir() const {
	return ((r_angle == 0.0) ? plus_z : trunk_cylin.get_norm_dir_vect());
}

cylinder_3dw small_tree::get_trunk_cylin() const {

	if (type == T_BUSH) {return cylinder_3dw();}
	vector3d dir(plus_z);
	rotate_vector3d(vector3d(rx, ry, 0.0), -r_angle/TO_DEG, dir); // oops, rotation is backwards
	bool const is_pine(is_pine_tree());
	float const hval(is_pine ? 1.0 : 0.8), zb(pos.z - 0.2*width), zbot(get_tree_z_bottom(zb, pos)), len(hval*height + (zb - zbot));
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
	cp.color       = color;
	cp_trunk.color = bark_color;
	if (type != T_BUSH && type != T_SH_PINE) {coll_id.push_back(add_coll_cylinder(trunk_cylin, cp_trunk, -1, 1));}
	vector3d const dirh(get_rot_dir()*height);

	switch (type) {
	case T_PINE: // pine tree
	case T_SH_PINE: // short pine tree
		// Note2: smaller than the actual pine tree render at the base so that it's easier for the player to walk under pine trees
		// Note2: could use the actual pine tree branch polygons, but that may be too conservative and too slow
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
		assert(palm_verts != nullptr && !palm_verts->empty());
		for (unsigned i = 0; i < palm_verts->size(); i += 4) { // iterate over the quads
			point pts[4];
			for (unsigned p = 0; p < 4; ++p) {pts[p] = (*palm_verts)[i+p].v;}
			cp.color = color.modulate_with((*palm_verts)[i].get_c4());
			coll_id.push_back(add_coll_polygon(pts, 4, cp, 0.0, -1, 1));
		}
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
	return sphere_vert_cylin_intersect(center, radius, trunk_cylin);
}


bool small_tree::line_intersect(point const &p1, point const &p2, float *t) const {

	assert(is_pine_tree()); // Note: can work on other tree types, but it's more complex, and we only need pine trees in tiled terrain mode
	vector3d const dirh(get_rot_dir()*height);
	cylinder_3dw const cylins[2] = {trunk_cylin, cylinder_3dw((pos + ((type == T_PINE) ? 0.35*dirh : all_zeros)), (pos + dirh), get_pine_tree_radius(), 0.0)};
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


void small_tree::calc_palm_tree_points() {

	if (palm_verts != nullptr) return;
	unsigned const num_fronds = 20;
	float const frond_l(0.3*height), frond_hw(0.2*width), frond_dz(-0.2*width);
	vector3d const trunk_dir(get_rot_dir());
	palm_verts.reset(new vector<vert_norm_tc_color>(8*num_fronds));
	vector<vert_norm_tc_color> &verts(*palm_verts);
	rand_gen_t rgen;
	rgen.set_state(long(1000*color.R), long(1000*color.G)); // seed random number generator with the tree color, which is deterministic

	for (unsigned n = 0, vix = 0; n < num_fronds; ++n) {
		vector3d const dir(rgen.signed_rand_vector_norm(1.0));
		float const dz(-0.4*width*rgen.rand_float());
		vector3d const binorm(cross_product(dir, plus_z).get_norm());
		vector3d const pa(trunk_cylin.p2 + trunk_dir*dz), pb(pa + vector3d(0.0, 0.0, frond_dz)), dx(frond_hw*binorm), dy(frond_l*dir);
		point const p0(pb-dx), p14(pa), p27(pa+dy), p3(pb-dx+dy), p5(pb+dx), p6(pb+dx+dy);
		vector3d const n1(cross_product(p14-p0, p27-p14).get_norm()), n2(cross_product(p5-p14, p6-p5).get_norm());
		float const brownness(0.5*rgen.rand_float() + 0.5*max(0.0f, -dir.z)); // fronds pointing down are browner
		color_wrapper_ctor cw(color.modulate_with(BROWN*brownness + WHITE*(1.0-brownness))); // random per-frond color
		verts[vix++].assign(p27, n1, 0.5, 1.0, cw.c); // 2
		verts[vix++].assign(p3,  n1, 0.0, 1.0, cw.c); // 3
		verts[vix++].assign(p0,  n1, 0.0, 0.0, cw.c); // 0
		verts[vix++].assign(p14, n1, 0.5, 0.0, cw.c); // 1
		verts[vix++].assign(p6,  n2, 1.0, 1.0, cw.c); // 6
		verts[vix++].assign(p27, n2, 0.5, 1.0, cw.c); // 7
		verts[vix++].assign(p14, n2, 0.5, 0.0, cw.c); // 4
		verts[vix++].assign(p5,  n2, 1.0, 0.0, cw.c); // 5
	}
}


float small_tree::get_pine_tree_radius() const {

	float const height0(((type == T_PINE) ? 0.75 : 1.0)*height/tree_height_scale);
	return 0.35*(height0 + 0.03/tree_scale);
}

void small_tree::calc_points(vbo_vnc_block_manager_t &vbo_manager, bool low_detail, bool update_mode) {

	if (type == T_PALM) {calc_palm_tree_points(); return;} // palm tree
	if (!is_pine_tree()) return; // only for pine trees
	float const sz_scale(SQRT2*get_pine_tree_radius()), dz((type == T_PINE) ? 0.35*height : 0.0);

	if (!low_detail) { // high detail
		unsigned const npts(4*N_PT_LEVELS*N_PT_RINGS);
		float const rd(0.45), height0(((type == T_PINE) ? 0.75 : 1.0)*height), theta0((int(1.0E6*height0)%360)*TO_RADIANS);
		float level_dz(height0/(N_PT_LEVELS + 1.2)), rd_scale(1.0), height_off(1.8*level_dz);
		point const center(pos + point(0.0, 0.0, dz));
		vert_norm points[npts];

		for (unsigned j = 0, ix = 0; j < N_PT_LEVELS; ++j) {
			level_dz *= 0.9; rd_scale *= 1.2; // higher slope, closer spacing near the top levels
			float const sz(sz_scale*(N_PT_LEVELS - j - 0.4)/(float)N_PT_LEVELS), z(height_off - rd*sz);
			height_off += level_dz;

			for (unsigned k = 0; k < N_PT_RINGS; ++k) {
				float const theta(TWO_PI*(3.3*j + k/(float)N_PT_RINGS) + theta0);
				add_rotated_quad_pts(points, ix, theta, z, center, sz, sz, sz, rd_scale*rd*sz); // bounds are (sz, sz, rd*sz+z)
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
		vert_norm vn(pos, vector3d(2.25*sz_scale/calc_tree_size(), 0.0, 0.0)); // ranges from around 0.25 to 0.75
		vn.v.z = pos.z + dz + 1.8*sz_scale + 0.1*height;
		points[0] = points[1] = vn; // top two vertices
		vn.v.z = pos.z + dz - 0.1*sz_scale - 0.2*height;
		points[2] = points[3] = vn; // bottom two vertices
		vbo_manager.add_points(points, 4, color);
	}
}

void small_tree::update_points_vbo(vbo_vnc_block_manager_t &vbo_manager, bool low_detail) {
	if (vbo_mgr_ix >= 0) {calc_points(vbo_manager, low_detail, 1);}
}


// approximate, trunk_cylin part is generally smaller so may not need to be added
float small_tree::get_zmax() const {return max((pos.z + height), trunk_cylin.p2.z);}


void small_tree::add_trunk_as_line(vector<point> &points) const {

	float const hval(is_pine_tree() ? 1.0 : 0.75);
	points.push_back(pos);
	points.push_back(pos + get_rot_dir()*(hval*height*((stt[type].w2 == 0.0) ? 0.7 : 1.0))); // slightly shorter for distant pine trees
}


void small_tree::draw_pine(vbo_vnc_block_manager_t const &vbo_manager, unsigned num_instances) const { // 30 quads per tree

	assert(is_pine_tree());
	assert(vbo_mgr_ix >= 0);
	vbo_manager.render_range(vbo_mgr_ix, vbo_mgr_ix+1, num_instances); // use glDrawElementsIndirect()?
}


bool small_tree::are_leaves_visible(vector3d const &xlate) const {

	if (type == T_PALM) { // slower, use occlusion culling
		return sphere_in_camera_view((trunk_cylin.p2 - 0.2*width*get_rot_dir() + xlate), (0.3*height + 0.2*width), 2);
	}
	else {
		return camera_pdu.sphere_visible_test((pos + 0.5*height*get_rot_dir() + xlate), max(1.5*width, 0.5*height));
	}
}


void small_tree::draw_pine_leaves(vbo_vnc_block_manager_t const &vbo_manager, vector3d const &xlate) const {
	if (is_pine_tree() && are_leaves_visible(xlate)) {draw_pine(vbo_manager);}
}


void small_tree::draw_trunks(bool shadow_only, vector3d const &xlate, vector<vert_wrap_t> *points, vector<vert_norm_tc> *cylin_verts) const {

	if (!(tree_mode & 2) || type == T_BUSH) return; // disabled, or no trunk/bark
	
	if (shadow_only) {
		if (!is_over_mesh(pos + xlate + height*get_rot_dir())) return;
	}
	else {
		if (!camera_pdu.sphere_visible_test((trunk_cylin.get_center() + xlate), (trunk_cylin.r1 + 0.5*trunk_cylin.get_length()))) return;
	}
	float const zoom_f(do_zoom ? ZOOM_FACTOR : 1.0), size_scale(zoom_f*stt[type].ss*width*window_width);
	float const dist(distance_to_camera(pos + xlate));
	if (!shadow_only && size_scale < dist) return; // too small/far

	if (!shadow_only && LINE_THRESH*zoom_f*(trunk_cylin.r1 + trunk_cylin.r2) < dist) { // draw as line
		point const p2((trunk_cylin.r2 == 0.0) ? (0.2*trunk_cylin.p1 + 0.8*trunk_cylin.p2) : trunk_cylin.p2);
				
		if (points) {
			points->push_back(trunk_cylin.p1 + xlate);
			points->push_back(p2 + xlate);
		}
		else {
			tree_scenery_pld.add_textured_line(trunk_cylin.p1+xlate, p2+xlate, bark_color, stt[type].bark_tid);
		}
	}
	else { // draw as cylinder
		int const nsides(max(3, min(N_CYL_SIDES, int(0.25*size_scale/dist))));

		if (cylin_verts && is_pine_tree()) {
			assert(trunk_cylin.r2 == 0.0); // cone
			point const ce[2] = {trunk_cylin.p1, trunk_cylin.p2};
			vector3d v12;
			gen_cone_triangles(*cylin_verts, gen_cylinder_data(ce, trunk_cylin.r1, trunk_cylin.r2, nsides, v12));
		}
		else {
			if (!shadow_only) {
				bark_color.set_for_cur_shader();
				select_texture(stt[type].bark_tid);
			}
			if (type == T_PALM && !shadow_only) {
				point const pa(0.92*trunk_cylin.p2 + 0.08*trunk_cylin.p1);
				draw_fast_cylinder(trunk_cylin.p1, pa, trunk_cylin.r1, trunk_cylin.r2, nsides, 1, 0, 0, nullptr, 2.0);

				if (nsides >= 8) {
					color.set_for_cur_shader(); // palm frond color
					point const pb(0.98*trunk_cylin.p2 + 0.02*trunk_cylin.p1), pc(1.04*trunk_cylin.p2 - 0.04*trunk_cylin.p1);
					draw_fast_cylinder(pa, pb, trunk_cylin.r2, 0.8*trunk_cylin.r2, nsides, 1);
					draw_fast_cylinder(pb, pc, 0.8*trunk_cylin.r2, 0.0, nsides, 1);
				}
			}
			else {
				draw_fast_cylinder(trunk_cylin.p1, trunk_cylin.p2, trunk_cylin.r1, trunk_cylin.r2, nsides, !shadow_only);
			}
		}
	}
}


void small_tree::draw_leaves(bool shadow_only, int xlate_loc, int scale_loc, vector3d const &xlate) const {

	if (!(tree_mode & 2)) return; // disabled
	assert(!is_pine_tree()); // handled through draw_pine_leaves()
	if (shadow_only ? !is_over_mesh(pos + xlate + point(0.0, 0.0, 0.5*height)) : !are_leaves_visible(xlate)) return;

	if (type == T_PALM) {
		assert(palm_verts != nullptr);
		assert(xlate == zero_vector);
		draw_vect_quads(*palm_verts);
		return;
	}
	float const size_scale((do_zoom ? ZOOM_FACTOR : 1.0)*stt[type].ss*width*window_width);
	vector3d scale;
	point xl(all_zeros);

	switch (type) { // draw leaves (palm or deciduous)
	case T_DECID: // decidious tree
		xl.z  = 0.75*height;
		scale = width*vector3d(1.2, 1.2, 0.8);
		break;
	case T_TDECID: // tall decidious tree
		xl.z  = 1.0*height;
		scale = width*vector3d(0.7, 0.7, 1.6);
		break;
	case T_BUSH: // bush
		scale = vector3d((0.1*height+0.8*width), (0.1*height+0.8*width), width);
		break;
	}
	int const nsides(max(6, min(N_SPHERE_DIV, (shadow_only ? get_def_smap_ndiv(width) : (int)(size_scale/distance_to_camera(pos + xlate))))));

	if (r_angle != 0.0) {
		fgPushMatrix();
		translate_to(pos);
		fgRotate(r_angle, rx, ry, 0.0);
	}
	else {xl += pos;}
	assert(xlate_loc >= 0 && scale_loc >= 0);
	shader_t::set_uniform_vector3d(xlate_loc, xl);
	shader_t::set_uniform_vector3d(scale_loc, scale);
	if (!shadow_only) {color.set_for_cur_shader();}
	draw_sphere_vbo_pre_bound(nsides, !shadow_only, (type == T_PALM));
	if (r_angle != 0.0) {fgPopMatrix();}
}


void small_tree::write_to_cobj_file(std::ostream &out) const {
	// 'F': // place small tree: xpos ypos height width type [zpos], type: T_PINE = 0, T_DECID = 1, T_TDECID = 2, T_BUSH = 3, T_PALM = 4, T_SH_PINE = 5
	//add_small_tree(pos, xf.scale*fvals[0], xf.scale*fvals[1], ivals[0], !use_z)
	out << "F " << pos.x << " " << pos.y << " " << height << " " << width << " " << int(type) << " " << pos.z << endl;
}
void write_small_trees_to_cobj_file(std::ostream &out) {
	for (auto i = small_trees.begin(); i != small_trees.end(); ++i) {i->write_to_cobj_file(out);}
}
