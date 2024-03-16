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

unsigned const PINE_TREE_NPTS = 4*N_PT_LEVELS*N_PT_RINGS;


small_tree_group small_trees;
small_tree_group tree_instances;
pt_line_drawer tree_scenery_pld;

extern bool tree_indir_lighting, only_pine_palm_trees;
extern int window_width, draw_model, num_trees, do_zoom, tree_mode, xoff2, yoff2;
extern int rand_gen_index, display_mode, force_tree_class, mesh_gen_mode;
extern unsigned max_unique_trees;
extern float zmin, zmax_est, water_plane_z, tree_scale, sm_tree_density, vegetation, tree_density_thresh, tree_height_scale, sm_tree_scale, CAMERA_RADIUS, tree_type_rand_zone, pine_tree_radius_scale;
extern tree_placer_t tree_placer;


struct sm_tree_type {

	int leaf_tid, bark_tid;
	float w2, ws, h, ss, width_scale, height_scale;
	colorRGB bc; // trunk color

	sm_tree_type(float w2_, float ws_, float h_, float ss_, float wscale, float hscale, colorRGB const &bc_, int ltid, int btid)
		: leaf_tid(ltid), bark_tid(btid), w2(w2_), ws(ws_), h(h_), ss(ss_), width_scale(wscale), height_scale(hscale), bc(bc_) {}
};

sm_tree_type const stt[NUM_ST_TYPES] = { // w2, ws, h, ss, wscale, hscale, c, tid
	sm_tree_type(0.00, 0.14, 0.35, 0.4, 1.0, 1.2, PTREE_C, PINE_TEX,      BARK2_TEX), // T_PINE
	sm_tree_type(0.13, 0.15, 0.75, 0.8, 1.0, 1.0, TREE_C,  TREE_HEMI_TEX, BARK3_TEX), // T_DECID // HEDGE_TEX?
	sm_tree_type(0.13, 0.15, 0.75, 0.7, 1.0, 1.0, TREE_C,  HEDGE_TEX,     BARK1_TEX), // T_TDECID
	sm_tree_type(0.00, 0.15, 0.00, 0.8, 1.0, 1.0, TREE_C,  HEDGE_TEX,     BARK1_TEX), // T_BUSH NOTE: bark texture is not used in trees, but is used in logs
	sm_tree_type(0.03, 0.15, 1.00, 0.6, 1.4, 2.0, TREE_C,  PALM_FROND_TEX,PALM_BARK_TEX), // T_PALM
	sm_tree_type(0.00, 0.08, 0.00, 0.4, 1.2, 0.8, PTREE_C, PINE_TEX,      BARK2_TEX), // T_SH_PINE
};

bool is_pine_tree_type(int type) {return (type == T_PINE || type == T_SH_PINE);}
bool small_trees_enabled() {return ((tree_mode & 2) || tree_placer.have_small_trees());}
bool are_trees_enabled  () {return ((tree_mode & 1) || tree_placer.have_decid_trees());}

int get_bark_tex_for_tree_type(int type) {
	assert(type < NUM_ST_TYPES);
	return stt[type].bark_tid;
}

colorRGBA get_tree_trunk_color(int type, bool modulate_with_texture) {
	int const bark_tid(get_bark_tex_for_tree_type(type));
	colorRGBA c(stt[type].bc);
	if (modulate_with_texture && bark_tid >= 0) {return c.modulate_with(texture_color(bark_tid));}
	return c;
}

unsigned get_tree_inst_gpu_mem() {return tree_instances.get_gpu_mem();}


void small_tree_group::add_tree(small_tree const &st) {

	if      (st.is_pine_tree()      ) {++num_pine_trees;}
	else if (st.get_type() == T_PALM) {++num_palm_trees;}
	max_tree_radius = max(max_tree_radius, st.get_radius()); // includes all tree types
	push_back(st);
}


void small_tree_group::finalize(bool low_detail) {

	palm_vbo_mem = 0;
	if (empty()) return;
	assert(!is_uploaded(low_detail));
	vbo_vnc_block_manager_t &vbo_mgr(vbo_manager[low_detail]);
	vbo_mgr.clear(0); // clear_pts_mem = 0
	vbo_mgr.reserve_pts(num_pine_trees*(low_detail ? 1 : PINE_TREE_NPTS));
	if (!low_detail) {vbo_mgr.reserve_offsets(num_pine_trees);}
#pragma omp parallel for schedule(static,1) num_threads(4) if (!low_detail)
	for (int i = 0; i < (int)size(); ++i) {operator[](i).calc_points(vbo_mgr, low_detail);}

	if (num_pine_trees > 0) {
		for (const_iterator i = begin(); i != end(); ++i) {palm_vbo_mem += i->get_palm_mem();}
	}
}


void small_tree_group::finalize_upload_and_clear_pts(bool low_detail) { // called in tiled terrain mode

	if (empty() || is_uploaded(low_detail)) return;
	//timer_t timer("Tree Finalize");

	if (instanced && !low_detail) {
		tree_instances.finalize_upload_and_clear_pts(0); // high detail
		return;
	}
	//RESET_TIME;
	static vector<vert_norm_comp_color> reused_pts;
	vbo_vnc_block_manager_t &vbo_mgr(vbo_manager[low_detail]);
	vbo_mgr.swap_points(reused_pts); // use the possibly pre-allocated points rather than allocating a new vector
	finalize(low_detail);
	//if (!low_detail) {PRINT_TIME("Finalize");}
	vbo_mgr.upload();
	vbo_mgr.swap_points(reused_pts); // make reused_pts available to other tiles
	//if (!low_detail) {PRINT_TIME("Finalize + Upload");}
}


void small_tree_group::draw_trunk_pts() {

	if (!trunk_pts_vbo.vbo_valid()) { // create trunk points
		vector<point> trunk_pts;
		trunk_pts.reserve(2*size());
		for (iterator i = begin(); i != end(); ++i) {i->add_trunk_as_line(trunk_pts);}
		trunk_pts_vbo.create_and_upload(trunk_pts);
		num_trunk_pts = trunk_pts.size();
	}
	trunk_pts_vbo.pre_render();
	vert_wrap_t::set_vbo_arrays();
	glDrawArrays(GL_LINES, 0, num_trunk_pts);
	++num_frame_draw_calls;
	bind_vbo(0);
}


void small_tree_group::clear_vbos() {

	trunk_pts_vbo.clear_vbo();
	for (unsigned i = 0; i < 2; ++i) {vbo_manager[i].clear_vbo();}
	if (num_palm_trees == 0) return; // no palm tree vbos to clear
	for (iterator i = begin(); i != end(); ++i) {i->clear_vbo();}
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

	clear_vbos(); // required to clear palm tree VBOs
	clear();
	clear_vbo_manager_and_ids();
	generated       = 0;
	max_tree_radius = 0.0;
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


void small_tree_group::calc_bcube() {
	all_bcube.set_to_zeros();
	for (iterator i = begin(); i != end(); ++i) {i->add_bounds_to_bcube(all_bcube);}
}


bool small_tree_group::check_sphere_coll(point &center, float radius) const {

	if (!all_bcube.is_zero_area() && !sphere_cube_intersect(center, radius, all_bcube)) return 0;
	bool coll(0);

	for (const_iterator i = begin(); i != end(); ++i) {
		coll |= i->check_sphere_coll(center, radius); // Note: calculates/updates center, so we can't early terminate in case there are multiple collisions
	}
	return coll;
}


// if t is non-NULL, calculate the point of intersection
bool small_tree_group::line_intersect(point const &p1, point const &p2, float *t) const { // if t != NULL, it should start at 1.0 (the end of the line)

	if (!all_bcube.is_zero_area() && !check_line_clip(p1, p2, all_bcube.d)) return 0;
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


bool small_tree_group::draw_trunks(bool shadow_only, bool all_visible, bool skip_lines, vector3d const &xlate) const {

	if (!all_visible && !(shadow_only && world_mode == WMODE_GROUND) && !camera_pdu.cube_visible(all_bcube + xlate)) return 0; // VFC
	static vector<vert_norm_tc> cylin_verts; // class member?
	bool all_drawn(1);

	for (const_iterator i = begin(); i != end(); ++i) {
		if (!i->draw_trunk(shadow_only, all_visible, skip_lines, xlate, &cylin_verts)) {all_drawn = 0;}
	}
	if (!cylin_verts.empty()) {
		if (!shadow_only) {
			get_tree_trunk_color(T_PINE, 0).set_for_cur_shader(); // all a constant color
			select_texture(get_bark_tex_for_tree_type(T_PINE));
		}
		// make this work for non-pine trees (where we have to deal with colors, texture changes, and tris vs. strips)?
		draw_and_clear_verts(cylin_verts, GL_TRIANGLES); // inf terrain mode pine tree trunks
	}
	return all_drawn;
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
		if (i->are_leaves_visible(xlate)) {to_draw.emplace_back(p2p_dist_sq(i->get_pos(), ref_pos), i-begin());}
	}
	sort(to_draw.begin(), to_draw.end()); // sort front to back for early Z culling
}


void small_tree_group::draw_tree_insts(shader_t &s, bool draw_all, vector3d const &xlate, int xlate_loc, vector<tree_inst_t> &insts, bool is_pine) {

	assert(instanced);
	unsigned const num_of_this_type(is_pine ? num_pine_trees : num_palm_trees);
	if (num_of_this_type == 0) return;

	if (!draw_all || insts.size() != num_of_this_type) { // recompute insts
		insts.clear();
		if (!draw_all && !camera_pdu.cube_visible(all_bcube + xlate)) return; // VFC
		insts.reserve(num_of_this_type);

		for (const_iterator i = begin(); i != end(); ++i) {
			if ((is_pine && !i->is_pine_tree()) || (!is_pine && i->get_type() != T_PALM)) continue; // only pine/plam trees are instanced
			if (draw_all || i->are_leaves_visible(xlate)) {insts.emplace_back(i->get_inst_id(), i->get_pos());}
		}
		sort(insts.begin(), insts.end());
	}
	if (insts.empty()) return; // no insts for this tree type
	inst_pts.resize(insts.size());
	for (auto i = insts.begin(); i != insts.end(); ++i) {inst_pts[i-insts.begin()] = i->pt;}
	vbo_vnc_block_manager_t const &vbomgr(tree_instances.vbo_manager[0]); // high detail, only used for pine trees
	if (is_pine) {vbomgr.begin_render();}
	select_texture((draw_model != 0) ? WHITE_TEX : stt[is_pine ? (unsigned)T_PINE : (unsigned)T_PALM].leaf_tid);
	void const *vbo_ptr(get_dynamic_vbo_ptr(&inst_pts.front(), inst_pts.size()*sizeof(point)));
	unsigned ptr_offset(0), ix(0);
	assert(xlate_loc >= 0);

	for (auto i = insts.begin(); i != insts.end();) { // Note: no increment
		unsigned const inst_id(i->id);
		assert(inst_id < tree_instances.size());
		for (ix = 0; i != insts.end() && i->id == inst_id; ++i, ++ix) {}
		if (!is_pine) {bind_dynamic_vbo();} // need to rebind VBO since draw_palm_leaves() will bind a different VBO
		glVertexAttribPointer(xlate_loc, 3, GL_FLOAT, GL_FALSE, sizeof(point), (void const *)((point const *)vbo_ptr + ptr_offset));
		if (is_pine) {tree_instances[inst_id].draw_pine(vbomgr, ix);} else {tree_instances[inst_id].draw_palm_leaves(ix);}
		ptr_offset += ix;
	}
	if (is_pine) {vbomgr.end_render();}
}

void small_tree_group::draw_pine_leaves(shader_t &s, bool shadow_only, bool low_detail, bool draw_all, bool sort_front_to_back, vector3d const &xlate, int xlate_loc) {

	if (empty()) return;

	if (instanced && !low_detail) { // high detail instanced; sort_front_to_back not supported
		draw_pine_insts(s, draw_all, xlate, xlate_loc);
		return;
	}
	if (num_pine_trees == 0) return;
	if (!draw_all && !camera_pdu.cube_visible(all_bcube + xlate)) return; // VFC
	select_texture((draw_model != 0) ? WHITE_TEX : (low_detail ? PINE_TREE_TEX : stt[T_PINE].leaf_tid));
	vbo_vnc_block_manager_t const &vbomgr(vbo_manager[low_detail]);// non-instanced (pine trees only)
	vbomgr.begin_render();

	if (draw_all) {
		if (low_detail) {vbo_manager[1].set_prim_type(GL_POINTS);}
		vbomgr.render_all();
	}
	else if (sort_front_to_back) {
		assert(!low_detail);
		vector<pair<float, unsigned> > to_draw;
		get_back_to_front_ordering(to_draw, xlate);
		for (unsigned i = 0; i < to_draw.size(); ++i) {operator[](to_draw[i].second).draw_pine_leaves(vbomgr, xlate);}
	}
	else {
		for (const_iterator i = begin(); i != end(); ++i) {i->draw_pine_leaves(vbomgr, xlate);}
	}
	vbomgr.end_render();
}


void small_tree_group::draw_non_pine_leaves(bool shadow_only, bool draw_palm, bool draw_non_palm, int xlate_loc, int scale_loc, vector3d const &xlate) const {

	if (!draw_non_palm && num_palm_trees == 0) return; // no palm trees to draw
	if (!shadow_only && !camera_pdu.cube_visible(all_bcube + xlate)) return; // VFC
	
	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->is_pine_tree()) continue;
		int const type(i->get_type());
		if (!((type == T_PALM) ? draw_palm : draw_non_palm)) continue;
		if (i == begin() || (i-1)->get_type() != type) {select_texture(stt[type].leaf_tid);} // first of this type
		i->draw_leaves(shadow_only, xlate_loc, scale_loc, xlate);
	}
}


float calc_tree_scale() {return (Z_SCENE_SIZE*tree_scale)/16.0f;}
float calc_tree_size () {return 16.0f*SM_TREE_SIZE/tree_scale;}

float rand_tree_height(rand_gen_t &rgen) {return rgen.rand_uniform(0.4, 1.0);}
float rand_tree_width (rand_gen_t &rgen) {return rgen.rand_uniform(0.25, 0.35);}

struct tree_inst_range_t {
	unsigned start, end;
	tree_inst_range_t() : start(0), end(0) {}
	unsigned select_inst(rand_gen_t &rgen) const {
		assert(start < end);
		return (start + rgen.rand()%(end - start));
	}
};
tree_inst_range_t num_insts_per_type[NUM_ST_TYPES];


void create_pine_tree_instances() {

	unsigned const num_pine_unique(max(1U, max_unique_trees/2)); // don't need as many insts as decid trees
	unsigned const num_palm_unique(max(1U, max_unique_trees/2)); // don't need as many insts as decid trees
	if (tree_instances.size() == num_pine_unique+num_palm_unique) return; // already done
	assert(tree_instances.empty());
	tree_instances.reserve(num_pine_unique+num_palm_unique);
	rand_gen_t rgen;

	for (unsigned i = 0; i < num_pine_unique; ++i) { // Note: pine trees not rotated
		int const ttype((rgen.rand()%10 == 0) ? (int)T_SH_PINE : (int)T_PINE);
		float const height(rand_tree_height(rgen));
		tree_instances.add_tree(small_tree(all_zeros, height, height*rand_tree_width(rgen), ttype, 0, rgen));
	}
	num_insts_per_type[T_PINE].end = num_insts_per_type[T_PALM].start = tree_instances.size();

	for (unsigned i = 0; i < num_palm_unique; ++i) { // Note: palm trees not rotated
		float const height(rand_tree_height(rgen));
		tree_instances.add_tree(small_tree(all_zeros, height, height*rand_tree_width(rgen), T_PALM, 0, rgen));
	}
	num_insts_per_type[T_PALM].end = tree_instances.size();
	num_insts_per_type[T_SH_PINE]  = num_insts_per_type[T_PINE]; // share the same instance range
}


int small_tree_group::get_ntrees_for_mesh_xy(int i, int j, float ntrees_mult_density) {

	int const ntrees(int(min(1.0f, ntrees_mult_density)*NUM_SMALL_TREES));
	if (ntrees == 0) return 0;
	
	if (XY_MULT_SIZE >= 2*ntrees) { // otherwise we will be modding with a value <= 1
		rgen.set_state((657435*(i + yoff2) + 243543*(j + xoff2) + 734533*rand_gen_index), (845631*(j + xoff2) + 667239*(i + yoff2) + 846357*rand_gen_index));
		rgen.rand(); // increase randomness (slower, but less regular - which is important for pine trees in tiled terrain mode)
		if ((rgen.rand_seed_mix()%(XY_MULT_SIZE/ntrees)) != 0) return 0; // not selected based on tree probability
	}
	return max(1, ntrees/XY_MULT_SIZE);
}

void small_tree_group::maybe_add_tree(int i, int j, float zpos_in, float tsize, int skip_val, bool check_hmap_normal) {

	// Note: faster to do this check when terrain is rough
	if (check_hmap_normal && get_tiled_terrain_height_tex_norm(j+xoff2, i+yoff2).z < 0.8) return; // 0.7ms
	rgen.rand_mix();
	float const xval(get_xval(j) + 0.5*skip_val*DX_VAL*rgen.signed_rand_float());
	float const yval(get_yval(i) + 0.5*skip_val*DY_VAL*rgen.signed_rand_float());
	float const zpos((zpos_in != 0.0) ? zpos_in : interpolate_mesh_zval(xval, yval, 0.0, 1, 1)); // 0.9ms
	int const ttype(get_tree_type_from_height(zpos, rgen, 0)); // 0.3ms
	if (ttype == TREE_NONE) return;
	small_tree tree;
	point pos(xval, yval, zpos);

	if (instanced) { // check voxel terrain? or does this mode only get used for tiled terrain?
		assert(ttype == T_SH_PINE || ttype == T_PINE || ttype == T_PALM);
		assert(world_mode == WMODE_INF_TERRAIN); // must be inf terrain mode
		tree = small_tree(pos, num_insts_per_type[ttype].select_inst(rgen)); // 0.5ms
	}
	else {
		float const height(tsize*rand_tree_height(rgen));
		if (world_mode != WMODE_INF_TERRAIN) {pos.z -= 0.1*height;} // move down slightly to ensure the trunk starts under the terrain
		if (point_inside_voxel_terrain(pos)) return; // don't create trees that start inside voxels (but what about trees that grow into voxels?)
		tree = small_tree(pos, height, height*rand_tree_width(rgen), ttype, 0, rgen, 1); // allow_rotation=1
	}
	if (!check_valid_scenery_pos(pos, 2.5*tree.get_radius(), 1)) return; // conservative bsphere; is_tall=1 1.4ms
	add_tree(tree); // 0.4ms
}

// density = x1,y1 x2,y1 x1,y2 x2,y2
void small_tree_group::gen_trees(int x1, int y1, int x2, int y2, float const density[4]) {

	//timer_t timer("Gen Trees");
	generated = 1; // mark as generated if we got here, even if there are no actual trees generated

	if (tree_placer.have_decid_trees()) { // now add pre-placed trees within the city (TT mode)
		vector3d const xlate(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0);
		cube_t const bounds(get_xval(x1), get_xval(x2), get_yval(y1), get_yval(y2), 0.0, 0.0); // Note: zvals are unused
		float const tsize(calc_tree_size());
		int const allowed_types[3] = {T_PINE, T_SH_PINE, T_PALM};

		if (bounds.intersects_xy(tree_placer.sm_bcube + xlate)) { // test full bcube
			for (auto b = tree_placer.sm_blocks.begin(); b != tree_placer.sm_blocks.end(); ++b) {
				if (!bounds.intersects_xy(b->bcube + xlate)) continue;

				for (auto t = b->trees.begin(); t != b->trees.end(); ++t) {
					point const pos(t->pos + xlate);
					if (!bounds.contains_pt_xy(pos)) continue; // tree not within this tile
					int const ttype(allowed_types[t->type%3]);

					if (instanced) {
						add_tree(small_tree(pos, num_insts_per_type[ttype].select_inst(rgen)));
					}
					else {
						float height(tsize*rand_tree_height(rgen));
						if (t->size > 0.0) {height *= t->size;} // treat size as a size scale
						add_tree(small_tree(pos, height, height*rand_tree_width(rgen), ttype, 0, rgen, 1, t->pine_xy_sz)); // allow_rotation=1
					}
				} // for t
			} // for b
		}
	}
	if (sm_tree_density == 0.0 || vegetation == 0.0 || !(tree_mode & 2)) {calc_bcube(); return;}
	if (density[0] == 0.0 && density[1] == 0.0 && density[2] == 0.0 && density[3] == 0.0) {calc_bcube(); return;}
	assert(x1 < x2 && y1 < y2);
	float const tscale(calc_tree_scale()), tsize(calc_tree_size()), ntrees_mult(vegetation*sm_tree_density*tscale*tscale/8.0f);
	int const skip_val(max(1, int(1.0/(sqrt(sm_tree_density*tree_scale)))));
	bool const use_hmap_tex(using_tiled_terrain_hmap_tex());
	// faster, but lower z-value accuracy, and only works for tiled terrain mode
	bool const approx_zval(world_mode == WMODE_INF_TERRAIN && !use_hmap_tex && (mesh_gen_mode == MGEN_SINE || ntrees_mult > 0.025));
	float const tds(TREE_DIST_SCALE*(XY_MULT_SIZE/16384.0)), xscale(tds*DX_VAL*DX_VAL), yscale(tds*DY_VAL*DY_VAL);
	mesh_xy_grid_cache_t density_gen, height_gen; // random tree generation based on transformed mesh height function
	density_gen.build_arrays(xscale*(x1 + xoff2), yscale*(y1 + yoff2), xscale, yscale, (x2-x1), (y2-y1), 0, 1); // force_sine_mode=1
	
	if (approx_zval) {
		height_gen.build_arrays(DX_VAL*(x1 + xoff2 - (MESH_X_SIZE >> 1) + 0.5f), DY_VAL*(y1 + yoff2 - (MESH_Y_SIZE >> 1) + 0.5f), DX_VAL, DY_VAL, (x2-x1), (y2-y1));
		height_gen.enable_glaciate();
	}
	float const dxv(skip_val/(x2 - x1 - 1.0f)), dyv(skip_val/(y2 - y1 - 1.0f));
	float xv(0.0), yv(0.0);

	for (int i = y1; i < y2; i += skip_val, yv += dyv) {
		for (int j = x1; j < x2; j += skip_val, xv += dxv) {
			float const cur_density(yv*(xv*density[3] + (1.0f-xv)*density[2]) + (1.0f-yv)*(xv*density[1] + (1.0f-xv)*density[0]));
			int const trees_this_xy(get_ntrees_for_mesh_xy(i, j, cur_density*ntrees_mult));
			if (trees_this_xy == 0) continue;
			float const hval(density_gen.eval_index(j-x1, i-y1));

			for (int n = 0; n < trees_this_xy; ++n) {
				if (hval > get_median_height(tree_density_thresh - TREE_DIST_RAND*rgen.rand_float())) continue; // tree density function test
				maybe_add_tree(i, j, (approx_zval ? height_gen.eval_index(j-x1, i-y1) : 0.0), tsize, skip_val, use_hmap_tex); // 4.3ms in this call
			} // for n
		} // for j
		xv = 0.0; // reset for next y iter
	} // for i
	if (world_mode == WMODE_GROUND) {sort_by_type();}
	calc_bcube();
}

// Note: for user placed trees in tiled terrain mode with heightmap texture; ignores vegetation, tree density functions, slope, etc.
void small_tree_group::gen_trees_tt_within_radius(int x1, int y1, int x2, int y2, point const &pos, float radius, bool is_square, float mesh_dz, tile_t const *const cur_tile) {

	generated = 1; // may already be set
	assert(x1 < x2 && y1 < y2);
	float const tscale(calc_tree_scale()), tsize(calc_tree_size()), ntrees_mult(sm_tree_density*tscale*tscale/8.0f);
	int const skip_val(max(1, int(1.0/(sqrt(sm_tree_density*tree_scale)))));

	for (int i = y1; i < y2; i += skip_val) {
		float const yval(get_yval(i));
		if (fabs(yval - pos.y) > radius) continue;

		for (int j = x1; j < x2; j += skip_val) {
			float const xval(get_xval(j));
			if (fabs(xval - pos.x) > radius) continue;
			point const tpos(xval, yval, 0.0);
			if (!is_square && !dist_xy_less_than(pos, tpos, radius)) continue; // Note: uses mesh xy center, not actual tree pos (for simplicity and efficiency)
			int const trees_this_xy(get_ntrees_for_mesh_xy(i, j, ntrees_mult));
			
			for (int n = 0; n < trees_this_xy; ++n) {
				rgen.rand_float(); // to match the get_median_height() call in gen_trees()
				maybe_add_tree(i, j, 0.0, tsize, skip_val, 0); // Note: mesh_dz is unused, but could be used for terrain slope estimates
			}
		} // for j
	} // for i
	calc_bcube();
}


bool small_tree_group::update_zvals(int x1, int y1, int x2, int y2) {

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


// use low bits of v to apply a random bias to v within the window [ref_pt +/- zone_width]
float val_signed_rand_bias_zone(float v, float ref_pt, float zone_width) {
	if (zone_width == 0.0) return v;
	float const dist(abs(v - ref_pt)), range(zone_width - dist);
	if (range <= 0.0) return v; // outside the zone, return the original value
	return (v + range*extract_low_bits_pm1(dist, 100.0/zone_width));
}
bool rel_height_check(float v, float thresh, float zw_scale=1.0) {
	return (val_signed_rand_bias_zone(v, thresh, zw_scale*tree_type_rand_zone) > thresh);
}

int get_tree_class_from_height(float zpos, bool pine_trees_only) {

	if (zpos < water_plane_z) return TREE_CLASS_NONE;
	
	if (force_tree_class >= 0) {
		assert(force_tree_class < NUM_TREE_CLASSES);
		return force_tree_class;
	}
	float relh(get_rel_height(zpos, -zmax_est, zmax_est));
	if (rel_height_check(relh, 0.9)) return TREE_CLASS_NONE; // too high
	if (rel_height_check(relh, 0.6)) return TREE_CLASS_PINE;
	//bool const allow_palm_trees(!pine_trees_only);
	bool const allow_palm_trees(tree_mode == 3);
	if (allow_palm_trees && zpos < 0.85*water_plane_z && !rel_height_check(relh, 0.435, 0.2)) return TREE_CLASS_PALM;
	if (pine_trees_only) {return ((tree_mode == 3) ? (int)TREE_CLASS_NONE : (int)TREE_CLASS_PINE);}
	return (only_pine_palm_trees ? (int)TREE_CLASS_PINE : (int)TREE_CLASS_DECID);
}

int get_tree_type_from_height(float zpos, rand_gen_t &rgen, bool for_scenery) {

	//return T_DECID; // TESTING
	switch (get_tree_class_from_height(zpos, (world_mode == WMODE_INF_TERRAIN && (tree_mode == 2 || (!for_scenery && tree_mode == 3))))) { // force pine trees in small tree mode
	case TREE_CLASS_NONE: return TREE_NONE;
	case TREE_CLASS_PINE: return ((rgen.rand()%10 == 0) ? (int)T_SH_PINE : (int)T_PINE);
	case TREE_CLASS_PALM: return T_PALM;
	case TREE_CLASS_DECID:
		//if (tree_mode == 3) return TREE_NONE; // use a large (complex) tree here
		return T_DECID + rgen.rand()%3;
	default: assert(0);
	}
	return TREE_NONE; // never gets here
}

bool can_have_pine_palm_trees_in_zrange(float z_min, float z_max, bool skip_range_check_if_manually_placed) {

	if (!small_trees_enabled()) return 0;
	if (z_max < water_plane_z)  return 0; // underwater
	if (force_tree_class >= 0) {return (force_tree_class != TREE_CLASS_NONE);}
	if (skip_range_check_if_manually_placed && tree_placer.have_small_trees()) return 1;
	float relh1(get_rel_height(z_min, -zmax_est, zmax_est)), relh2(get_rel_height(z_max, -zmax_est, zmax_est));
	if (relh1 - tree_type_rand_zone > 0.9f) return 0; // too high
	if (relh2 + tree_type_rand_zone > 0.6f) return 1; // can have pine trees
	if (tree_mode != 3) return 1; // no decid trees
	return (z_min < 0.85*water_plane_z && relh1 - 0.2*tree_type_rand_zone < 0.435); // can have palm trees
}

bool can_have_decid_trees_in_zrange(float z_min, float z_max, bool skip_range_check_if_manually_placed) {

	if (!are_trees_enabled())  return 0;
	if (z_max < water_plane_z) return 0; // underwater
	if (skip_range_check_if_manually_placed && tree_placer.have_decid_trees()) return 1;
	float relh1(get_rel_height(z_min, -zmax_est, zmax_est));
	if (relh1 - tree_type_rand_zone > 0.6f) return 0; // must have pine trees, or too high
	return 1; // not worth checking plam tree case, since the rand zone overlaps with the water plane
}

int add_small_tree(point const &pos, float height, float width, int tree_type, bool calc_z) {

	assert(height > 0.0 && width > 0.0);
	small_trees.add_tree(small_tree(pos, height, width, (abs(tree_type)%NUM_ST_TYPES), calc_z, small_trees.rgen)); // could have a type error
	small_trees.back().calc_points(small_trees.vbo_manager[0], 0);
	small_trees.back().add_bounds_to_bcube(small_trees.all_bcube);
	return 1; // might return zero in some case
}

colorRGBA colorgen(float r1, float r2, float g1, float g2, float b1, float b2, rand_gen_t &rgen) {
	return colorRGBA(rgen.rand_uniform(r1, r2), rgen.rand_uniform(g1, g2), rgen.rand_uniform(b1, b2), 1.0);
}


void gen_small_trees() { // called in ground mode

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

void add_small_tree_coll_objs() {small_trees.add_cobjs();} // doesn't handle rotation angle
void remove_small_tree_cobjs () {small_trees.remove_cobjs();}

void shift_small_trees(vector3d const &vd) {
	if (num_trees > 0) return; // dynamically created, not placed
	small_trees.translate_by(vd);
}

bool update_small_tree_zvals(int x1, int y1, int x2, int y2) {
	return small_trees.update_zvals(x1, y1, x2, y2);
}

void draw_small_trees(bool shadow_only, int reflection_pass) {small_trees.draw(shadow_only, reflection_pass);}

void small_tree_group::draw(bool shadow_only, int reflection_pass) {

	//RESET_TIME;
	if (empty() || !small_trees_enabled()) return;
	if (!all_bcube.is_zero_area() && !camera_pdu.cube_visible(all_bcube)) return;
	if (!shadow_only && !reflection_pass && size() < 100) {sort_by_dist_to_camera();} // not in shadow pass, since trees usually don't overlap in z
	shader_t s;
	bool const all_pine(num_pine_trees == size());
	bool const v(!shadow_only), use_bump_map(USE_BUMP_MAP && !shadow_only && all_pine); // bump maps only work with pine tree trunks

	// draw trunks
	if (shadow_only) {
		s.begin_shadow_map_shader();
	}
	else {
		setup_smoke_shaders(s, 0.0, 0, 0, tree_indir_lighting, 1, 1, 0, 0, 2, use_bump_map, 0, 1, 0, 0.0, 0.0, 0, 0, 1); // dynamic lights, but no smoke, is_outside=1
		s.add_uniform_float("tex_scale_t", 5.0);
	}
	if (use_bump_map) {select_texture(BARK2_NORMAL_TEX, 5);}
	draw_trunks(shadow_only);
	if (!shadow_only) {s.add_uniform_float("tex_scale_t", 1.0);}
	s.end_shader();

	// draw leaves; not drawn with the shadow map shader in the shadow pass because it doesn't support leaf wind
	float const wind_mag(get_plant_leaf_wind_mag(shadow_only));

	if (num_pine_trees > 0) { // pine trees
		vbo_manager[0].upload();
		if (wind_mag > 0.0) {s.set_prefix("#define ENABLE_WIND", 0);} // VS
		s.set_prefix("#define NO_SPECULAR", 1); // FS - disable rain effect
		setup_smoke_shaders(s, 0.5, 3, 0, (v && tree_indir_lighting), v, v, 0, 0, (v ? 2 : 1), 0, 0, v, v, 0.0, 0.0, 0, 0, 1); // dynamic lights, but no smoke, texgen, is_outside=1
		setup_leaf_wind(s, wind_mag, 0);
		draw_pine_leaves(s, shadow_only);
		s.end_shader();
	}
	if (!all_pine) { // non-pine trees
		if (num_pine_trees + num_palm_trees < size()) { // deciduous trees
			s.begin_simple_textured_shader(0.75, !shadow_only); // with lighting, unless shadow_only
			bind_draw_sphere_vbo(1, 1); // texture, even in shadow pass, to handle alpha masking of some tree types
			int const xlate_loc(s.get_uniform_loc("xlate")), scale_loc(s.get_uniform_loc("scale"));
			draw_non_pine_leaves(shadow_only, 0, 1, xlate_loc, scale_loc);
			s.set_uniform_vector3d(xlate_loc, zero_vector);
			s.set_uniform_vector3d(scale_loc, vector3d(1,1,1));
			bind_vbo(0);
			s.end_shader();
		}
		if (num_palm_trees > 0) { // palm trees
			if (wind_mag > 0.0) {s.set_prefix("#define ENABLE_WIND", 0);} // VS
			setup_smoke_shaders(s, 0.75, 4, 0, 0, v, v, 0, 0, (v ? 2 : 1), 0, 0, 0, 1); // dynamic lights, but no smoke (slow, but looks better); use bent quad texgen mode 4
			setup_leaf_wind(s, wind_mag, 0);
			draw_non_pine_leaves(shadow_only, 1, 0);
			s.end_shader();
		}
	}
	if (!tree_scenery_pld.empty()) { // not drawn in the shadow pass
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
	assert(coll_id[0] < 0 && coll_id[1] < 0);
	init();
	inst_id = instance_id;
	pos     = p;
	float const tsize(calc_tree_size());
	width  *= tsize;
	height *= tsize;
	trunk_cylin = get_trunk_cylin();
}


small_tree::small_tree(point const &p, float h, float w, int t, bool calc_z, rand_gen_t &rgen, bool allow_rotation, float bxys) :
	type(t), inst_id(-1), height(h), width(w), branch_xy_scale(bxys), pos(p)
{
	init();
	height    *= tree_height_scale*sm_tree_scale;
	bark_color = stt[type].bc;
	if (!is_pine_tree()) {UNROLL_3X(bark_color[i_] = min(1.0, bark_color[i_]*(0.85 + 0.3f*rgen.randd()));)} // gen bark color for decid trees
	if (calc_z) {pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1) - 0.1*height;}
	width  *= stt[type].width_scale;
	height *= stt[type].height_scale;

	switch (type) {
	case T_PINE: // pine tree
		leaf_color = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35, rgen);
		break;
	case T_SH_PINE: // short pine tree
		leaf_color = colorgen(0.05, 0.1, 0.3, 0.6, 0.15, 0.35, rgen);
		break;
	case T_PALM: // palm tree
		leaf_color = colorgen(0.7, 0.8, 0.9, 1.0, 0.6, 0.7, rgen); // r1, r2, g1, g2, b1, b2
		break;
	case T_DECID: // decidious tree
		leaf_color = colorgen(0.5, 0.8, 0.7, 1.0, 0.4, 0.7, rgen);
		break;
	case T_TDECID: // tall decidious tree
		leaf_color = colorgen(0.4, 0.7, 0.8, 1.0, 0.3, 0.5, rgen);
		break;
	case T_BUSH: // bush
		leaf_color = colorgen(0.5, 0.8, 0.7, 1.0, 0.3, 0.5, rgen);
		pos.z += 0.3*height;
		if (rgen.rand()%100 < 50) {pos.z -= height*rgen.rand_uniform(0.0, 0.2);}
		break;
	default: assert(0);
	}
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
	if (r_angle != 0.0) {rotate_vector3d(vector3d(rx, ry, 0.0), -r_angle/TO_DEG, dir);} // oops, rotation is backwards
	bool const is_pine(is_pine_tree());
	float const hval(is_pine ? 1.0 : 0.8), zb(pos.z - 0.2*width), zbot(get_tree_z_bottom(zb, pos)), len(hval*height + (zb - zbot));
	float const mod_width(width*(is_pine ? 0.8f*len/(hval*height) : 1.0f));
	point const p1((pos + dir*(zbot - pos.z)));
	return cylinder_3dw(p1, (p1 + dir*len), stt[type].ws*mod_width, stt[type].w2*mod_width);
}


void small_tree::add_cobjs(cobj_params &cp, cobj_params &cp_trunk) {

	if (!is_over_mesh(pos)) return; // not sure why, but this makes drawing slower

	if (!DRAW_COBJS) {
		cp.tid       = stt[type].leaf_tid;
		cp_trunk.tid = stt[type].bark_tid;
	}
	cp.color       = leaf_color;
	cp_trunk.color = bark_color;
	if (type != T_BUSH && type != T_SH_PINE) {coll_id[0] = add_coll_cylinder(trunk_cylin, cp_trunk, -1, 1);}
	vector3d const dirh(get_rot_dir()*height);

	switch (type) {
	case T_PINE: // pine tree
	case T_SH_PINE: // short pine tree
		// Note2: smaller than the actual pine tree render at the base so that it's easier for the player to walk under pine trees
		// Note2: could use the actual pine tree branch polygons, but that may be too conservative and too slow
		coll_id[1] = add_coll_cylinder((pos + ((type == T_PINE) ? 0.35*dirh : all_zeros)), (pos + dirh), get_radius(), 0.0, cp, -1, 1);
		break;
	case T_DECID: // decidious tree
		coll_id[1] = add_coll_sphere((pos + 0.75*dirh), 1.3*width, cp, -1, 1);
		break;
	case T_TDECID: // tall decidious tree
		coll_id[1] = add_coll_sphere((pos + 0.8*dirh), 0.8*width, cp, -1, 1);
		break;
	case T_BUSH: // bush
		coll_id[1] = add_coll_sphere(pos, 1.2*width, cp, -1, 1);
		break;
	case T_PALM: // palm tree
		{
			colorRGBA const palm_tex_color(texture_color(PALM_FROND_TEX));
			assert(palm_verts != nullptr && !palm_verts->v.empty());
			for (unsigned i = 0; i < palm_verts->v.size(); i += 4) { // iterate over the quads
				point pts[4];
				for (unsigned p = 0; p < 4; ++p) {pts[p] = palm_verts->v[i+p].v;}
				cp.color = palm_tex_color.modulate_with(palm_verts->v[i].get_c4());
				palm_verts->coll_id.push_back(add_coll_polygon(pts, 4, cp, 0.0, -1, 1));
			}
			break;
		}
	default: assert(0);
	}
}


void small_tree::remove_cobjs() {

	if (palm_verts != nullptr) {
		for (unsigned j = 0; j < palm_verts->coll_id.size(); ++j) {remove_reset_coll_obj(palm_verts->coll_id[j]);}
		palm_verts->coll_id.clear();
	}
	for (unsigned d = 0; d < 2; ++d) {remove_reset_coll_obj(coll_id[d]);}
}


void small_tree::add_bounds_to_bcube(cube_t &bcube) const {
	bcube.assign_or_union_with_sphere(trunk_cylin.get_center(), get_xy_radius()); // approximate/conservative
}


// very simple check against trunk only, for collisions with a player walking on the ground
bool small_tree::check_sphere_coll(point &center, float radius) const {
	if (type == T_BUSH) return 0; // no trunk, not yet handled
	if (!dist_xy_less_than(center, pos, (radius + trunk_cylin.r1))) return 0; // optimization
	return sphere_vert_cylin_intersect(center, radius, trunk_cylin);
}


bool small_tree::line_intersect(point const &p1, point const &p2, float *t) const { // for tiled terrain mode

	if (!is_pine_tree()) return 0; // Note: can work on other tree types such as palms, but it's more complex
	vector3d const dirh(get_rot_dir()*height);
	cylinder_3dw const cylins[2] = {trunk_cylin, cylinder_3dw((pos + ((type == T_PINE) ? 0.35*dirh : all_zeros)), (pos + dirh), get_radius(), 0.0)};
	bool coll(0);
	
	for (unsigned i = 0; i < 2; ++i) {
		if (t == NULL) {
			coll |= line_intersect_cylinder(p1, p2, cylins[i], 0); // no check ends
			if (coll) break;
		}
		else {
			float t_new(0.0);
			
			if (line_intersect_trunc_cone(p1, p2, cylins[i].p2, cylins[i].p1, cylins[i].r2, cylins[i].r1, 0, t_new) && t_new < *t) { // r2 > r1
				*t   = t_new; // would be more efficient if we could pass *t in as t_max
				coll = 1;
			}
		}
	} // for i
	return coll;
}


void small_tree::calc_palm_tree_points() {

	if (palm_verts != nullptr) return;
	unsigned const num_fronds = 20;
	float const frond_l(0.3*height), frond_hw(0.2*width), frond_dz(-0.2*width);
	vector3d const trunk_dir(get_rot_dir());
	palm_verts = std::make_shared<palm_verts_t>();
	vector<vert_norm_comp_color> &verts(palm_verts->v);
	verts.resize(8*num_fronds);
	rand_gen_t rgen;
	rgen.set_state(long(1000*leaf_color.R), long(1000*leaf_color.G)); // seed random number generator with the tree color, which is deterministic

	for (unsigned n = 0, vix = 0; n < num_fronds; ++n) {
		vector3d const dir(rgen.signed_rand_vector_norm(1.0));
		float const dz(-0.4*width*rgen.rand_float());
		vector3d const binorm(cross_product(dir, plus_z).get_norm());
		vector3d const pa(trunk_cylin.p2 + trunk_dir*dz), pb(pa + vector3d(0.0, 0.0, frond_dz)), dx(frond_hw*binorm), dy(frond_l*dir);
		point const p0(pb-dx), p14(pa), p27(pa+dy), p3(pb-dx+dy), p5(pb+dx), p6(pb+dx+dy);
		vector3d const n1(cross_product(p14-p0, p27-p14).get_norm()), n2(cross_product(p5-p14, p6-p5).get_norm());
		float const brownness(0.5*rgen.rand_float() + 0.5*max(0.0f, -dir.z)); // fronds pointing down are browner
		color_wrapper_ctor cw(leaf_color.modulate_with(BROWN*brownness + WHITE*(1.0-brownness))); // random per-frond color
		verts[vix++].assign(p27, n1, cw.c); // 2
		verts[vix++].assign(p3,  n1, cw.c); // 3
		verts[vix++].assign(p0,  n1, cw.c); // 0
		verts[vix++].assign(p14, n1, cw.c); // 1
		verts[vix++].assign(p6,  n2, cw.c); // 6
		verts[vix++].assign(p27, n2, cw.c); // 7
		verts[vix++].assign(p14, n2, cw.c); // 4
		verts[vix++].assign(p5,  n2, cw.c); // 5
	} // for n
	// Note: vbo_manager is unused, but we could put the palm verts into it, though perf seems okay with individual VBOs
}


float small_tree::get_pine_tree_radius() const { // Note: doesn't include branch_xy_scale
	float const height0(((type == T_PINE) ? 0.75 : 1.0)*height/(tree_height_scale*sm_tree_scale));
	return 0.35*pine_tree_radius_scale*(height0 + 0.03/tree_scale);
}

void small_tree::alloc_pine_tree_pts(vbo_vnc_block_manager_t &vbo_manager) {
	if (is_pine_tree()) {vbo_mgr_ix = vbo_manager.alloc_points_with_offset(PINE_TREE_NPTS);}
}

void small_tree::calc_points(vbo_vnc_block_manager_t &vbo_manager, bool low_detail, bool update_mode) {

	if (type == T_PALM) {calc_palm_tree_points(); return;} // palm tree
	if (!is_pine_tree()) return; // only for pine trees
	float const sz_scale(SQRT2*get_pine_tree_radius()), dz((type == T_PINE) ? 0.35*height : 0.0);

	if (!low_detail) { // high detail
		rand_gen_t rgen;
		rgen.set_state(long(10000*height), long(10000*leaf_color.B));
		float const rd(0.45), height0(((type == T_PINE) ? 0.75f : 1.0f)*height), theta0((int(1.0E6f*(height0 + leaf_color.B))%360)*TO_RADIANS);
		float const branch_start_height((num_trees > 0) ? 0.4*height0*rgen.rand_float() : 0.0); // only used with generated trees, not placed trees
		float level_dz((height0 - branch_start_height)/(N_PT_LEVELS + 1.2f)), rd_scale(1.0), height_off(branch_start_height + 1.8f*level_dz);
		point const center(pos + point(0.0, 0.0, dz));
		vert_norm_comp points[PINE_TREE_NPTS];

		for (unsigned j = 0, ix = 0; j < N_PT_LEVELS; ++j) {
			level_dz *= 0.9; rd_scale *= 1.2; // higher slope, closer spacing near the top levels
			float const sz(sz_scale*(N_PT_LEVELS - j - 0.4)/(float)N_PT_LEVELS), z(height_off - rd*sz);
			height_off += level_dz;

			for (unsigned k = 0; k < N_PT_RINGS; ++k) {
				float const theta(TWO_PI*(3.3*j + k/(float)N_PT_RINGS) + theta0 + 0.5*rgen.signed_rand_float());
				float const zz(z + 0.4*sz*rgen.signed_rand_float()), xy_sz(branch_xy_scale*sz); // reducing xy_sz makes tree tall+thin with more vertical branches
				add_rotated_quad_pts(points, ix, theta, zz, center, xy_sz, xy_sz, xy_sz, rd_scale*rd*sz); // bounds are (xy_sz, xy_sz, rd*sz+z)
			}
		} // for j
		if (update_mode) {
			assert(vbo_mgr_ix >= 0);
			vbo_manager.update_range(points, PINE_TREE_NPTS, leaf_color, vbo_mgr_ix, vbo_mgr_ix+1);
		}
		/*else if (vbo_mgr_ix >= 0) { // already allocated, just copy the points (requires additional code to be enabled in small_tree_group::finalize())
			vbo_manager.fill_pts_from(points, PINE_TREE_NPTS, leaf_color, vbo_mgr_ix);
		}*/
		else { // we only get into this case when running in parallel
#pragma omp critical(pine_tree_vbo_update)
			vbo_mgr_ix = vbo_manager.add_points_with_offset(points, PINE_TREE_NPTS, leaf_color);
		}
	}
	else { // low detail billboard
		assert(!update_mode);
		float const zv1(0.75*dz); // shift slightly down to account for sparse tree texture image
		// 0.9x to prevent clipping above 1.0
		vbo_manager.add_point(vert_norm_comp_color((pos + point(0.0, 0.0, zv1)), vector3d(2.0f*sz_scale/calc_tree_size(), 0.9f*(height - zv1), 0.0f), leaf_color));
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
	vbo_manager.render_single(vbo_mgr_ix, num_instances); // use glDrawElementsIndirect()?
}


bool small_tree::are_leaves_visible(vector3d const &xlate) const {

	if (type == T_PALM) { // slower, use occlusion culling
		return sphere_in_camera_view((trunk_cylin.p2 - 0.2*width*get_rot_dir() + xlate), (0.3*height + 0.2*width), 2);
	}
	else if (r_angle == 0.0) { // vertical trunk - common case
		return camera_pdu.sphere_visible_test((pos + vector3d(0.0, 0.0, 0.5*height) + xlate), get_xy_radius());
	}
	else {
		return camera_pdu.sphere_visible_test((pos + 0.5*height*get_rot_dir() + xlate), get_xy_radius());
	}
}


void small_tree::draw_pine_leaves(vbo_vnc_block_manager_t const &vbo_manager, vector3d const &xlate) const {
	if (is_pine_tree() && are_leaves_visible(xlate)) {draw_pine(vbo_manager);}
}


bool small_tree::draw_trunk(bool shadow_only, bool all_visible, bool skip_lines, vector3d const &xlate, vector<vert_norm_tc> *cylin_verts) const {

	if (!small_trees_enabled() || type == T_BUSH) return 1; // disabled, or no trunk/bark
	
	if (all_visible) {}
	else if (shadow_only && world_mode == WMODE_GROUND) {
		if (!is_over_mesh(pos + xlate + height*get_rot_dir())) return 1;
	}
	else {
		if (!camera_pdu.sphere_visible_test((trunk_cylin.get_center() + xlate), get_trunk_bsphere_radius())) return 1;
	}
	float const zoom_f(do_zoom ? ZOOM_FACTOR : 1.0), size_scale(zoom_f*stt[type].ss*width*window_width);
	float const dist(distance_to_camera(pos + xlate));
	if (!shadow_only && size_scale < dist) return 1; // too small/far

	if (!shadow_only && LINE_THRESH*zoom_f*(trunk_cylin.r1 + trunk_cylin.r2) < dist) { // draw as line
		if (skip_lines) return 0; // skipped - this is the only case that returns 0
		point const p2((trunk_cylin.r2 == 0.0) ? (0.2*trunk_cylin.p1 + 0.8*trunk_cylin.p2) : trunk_cylin.p2);
		tree_scenery_pld.add_textured_line(trunk_cylin.p1+xlate, p2+xlate, bark_color, stt[type].bark_tid);
	}
	else { // draw as cylinder
		int const nsides(max(3, min(N_CYL_SIDES, int(0.25*size_scale/dist))));

		if (cylin_verts && is_pine_tree()) { // flatten the tip (truncated cone)?
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
					leaf_color.set_for_cur_shader(); // palm frond color
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
	return 1;
}


void small_tree::draw_palm_leaves(unsigned num_instances) const {

	assert(type == T_PALM);
	assert(palm_verts != nullptr);
	palm_verts->vbo.create_and_upload(palm_verts->v, 0, 0);
	palm_verts->vbo.pre_render();
	draw_quad_verts_as_tris((vert_norm_comp_color *)nullptr, palm_verts->v.size(), 0, num_instances);
	palm_verts->vbo.post_render();
}

void small_tree::draw_leaves(bool shadow_only, int xlate_loc, int scale_loc, vector3d const &xlate) const {

	if (!small_trees_enabled()) return; // disabled
	assert(!is_pine_tree()); // handled through draw_pine_leaves()
	if (shadow_only ? !is_over_mesh(pos + xlate + point(0.0, 0.0, 0.5*height)) : !are_leaves_visible(xlate)) return;
	if (type == T_PALM) {draw_palm_leaves(); return;}
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
		scale.assign((0.1*height+0.8*width), (0.1*height+0.8*width), width);
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
	if (!shadow_only) {leaf_color.set_for_cur_shader();}
	draw_sphere_vbo_pre_bound(nsides, !shadow_only, (type == T_PALM));
	if (r_angle != 0.0) {fgPopMatrix();}
}


void small_tree::write_to_cobj_file(std::ostream &out) const {
	// 'F': // place small tree: xpos ypos height width type [zpos], type: T_PINE = 0, T_DECID = 1, T_TDECID = 2, T_BUSH = 3, T_PALM = 4, T_SH_PINE = 5
	//add_small_tree(pos, xf.scale*fvals[0], xf.scale*fvals[1], ivals[0], !use_z)
	out << "F " << pos.x << " " << pos.y << " " << height/stt[type].height_scale << " " << width/stt[type].width_scale << " " << int(type) << " " << pos.z << endl;
}
void write_small_trees_to_cobj_file(std::ostream &out) {
	for (auto i = small_trees.begin(); i != small_trees.end(); ++i) {i->write_to_cobj_file(out);}
}
