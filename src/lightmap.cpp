// 3D World - lighting code, incuding static and dynamic lights, profile generation, and flow calculation
// by Frank Gennari
// 1/16/06
#include "3DWorld.h"
#include "mesh.h"
#include "csg.h"
#include "lightmap.h"
#include "gl_ext_arb.h"
#include "shaders.h"


bool const CAMERA_CANDLE_LT  = 0;
bool const CAMERA_FLASH_LT   = 0; // looks cool
bool const POLY_XY_OVER_CHK  = 0;
bool const DYNAMIC_LT_FLOW   = 1;
bool const SHOW_STAT_LIGHTS  = 0; // debugging
bool const SHOW_DYNA_LIGHTS  = 0; // debugging
unsigned const NUM_LT_SMOOTH = 2; // nominally 2
unsigned const NUM_XY_PASSES = 2; // nominally 2
unsigned const NUM_RAND_LTS  = 0;
unsigned const FLOW_CACHE_BS = 17;
unsigned const FLOW_CACHE_SZ = (1 << FLOW_CACHE_BS);

int      const START_LIGHT   = GL_LIGHT2;
int      const END_LIGHT     = GL_LIGHT7 + 1;
unsigned const MAX_LIGHTS    = unsigned(END_LIGHT - START_LIGHT);

float const CTHRESH          = 0.025;
float const MIN_LIGHT        = 0.0;
float const MAX_LIGHT        = 1.0;
float const Z_WT_SCALE       = 1.0;
float const XY_WT_SCALE      = 1.0;
float const LIGHT_SCALE      = 1.0;
float const LIGHT_OFFSET     = 0.02;
float const LIGHT_SPREAD     = 0.4;
float const PASS_WEIGHT_ATT  = 0.4;
float const UNDER_M_RECOVER  = 0.25;
float const Z_LT_ATTEN       = 0.99;
float const XY_LT_ATTEN      = 0.94;
float const DZ_VAL_SCALE     = 2.0;
float const SHIFT_VAL        = 0.5; // hack to fix some offset problem
float const SLT_LINE_TEST_WT = 0.5; // SLT_LINE_TEST_WT + SLT_FLOW_TEST_WT == 1.0
float const SLT_FLOW_TEST_WT = 0.5;
float const DLIGHT_AMBIENT   = 0.25; // in range [0.0, 1.0]
float const DLIGHT_DIFFUSE   = 0.75; // in range [0.0, 1.0], DLIGHT_AMBIENT + DLIGHT_DIFFUSE should be close to 1.0
float const LT_DIR_FALLOFF   = 0.005;
float const LT_DIR_FALLOFF_INV(1.0/LT_DIR_FALLOFF);
float const DARKNESS_THRESH  = 0.1;


bool using_lightmap(0), lm_alloc(0), has_dl_sources(0), has_dir_lights(0), use_dense_voxels(0);
unsigned dl_tid(0), elem_tid(0), gb_tid(0), flow_tid(0);
float DZ_VAL_INV2(DZ_VAL_SCALE/DZ_VAL), SHIFT_DX(SHIFT_VAL*DX_VAL), SHIFT_DY(SHIFT_VAL*DY_VAL);
float czmin0(0.0), lm_dz_adj(0.0);
float dlight_bb[3][2] = {0}, SHIFT_DXYZ[3] = {SHIFT_DX, SHIFT_DY, 0.0};
dls_cell **ldynamic = NULL;
vector<light_source> light_sources, dl_sources, dl_sources2; // static, dynamic {cur frame, next frame}
flow_cache_e flow_cache[FLOW_CACHE_SZ]; // 2MB
lmap_manager_t lmap_manager;


extern bool disable_shaders;
extern int animate2, display_mode, frame_counter, read_light_files[], write_light_files[];
extern unsigned num_vpls;
extern float czmin, czmax, fticks, zbottom, ztop, XY_SCENE_SIZE;
extern vector<coll_obj> coll_objects;
extern vector<light_source> enabled_lights;
extern vector<cube_light_source> global_cube_lights;


// *** USEFUL INLINES ***

inline bool add_cobj_ok(coll_obj const &cobj) { // skip small things like tree leaves and such
	
	return (cobj.fixed && !cobj.disabled() && cobj.volume > 0.0001); // cobj.type == COLL_CUBE
}


inline void reset_lighting(float &mesh_light, float &vals, float &vscale) {

	mesh_light = min(1.0f, (mesh_light + UNDER_M_RECOVER*(128.0f/MESH_X_SIZE))); // mesh accumulates ambient light
	vals       = mesh_light;
	vscale     = XY_WT_SCALE;
}


// *** LIGHT_SOURCE IMPLEMENTATION ***


// radius == 0.0 is really radius == infinity (no attenuation)
light_source::light_source(float sz, point const &p, colorRGBA const &c, bool id, vector3d const &d, float bw, float ri, int gl_id) :
	dynamic(id), gl_light_id(gl_id), radius(sz), radius_inv((radius == 0.0) ? 0.0 : 1.0/radius),
	r_inner(ri), bwidth(bw), center(p), dir(d.get_norm()), color(c)
{
	assert(bw > 0.0 && bw <= 1.0);
	calc_cent();
}


void light_source::calc_cent() {

	for (unsigned i = 0; i < 3; ++i) {
		cent[i] = max(0, min(MESH_SIZE[i]-1, get_dim_pos(center[i], i))); // clamp to mesh bounds
	}
}


void light_source::add_color(colorRGBA const &c) {

	color = color*color.alpha + c*c.alpha;
	color.alpha = 1.0;
}


float light_source::get_intensity_at(point const &pos) const {

	if (radius == 0.0)                   return color[3]; // no falloff
	if (fabs(pos.z - center.z) > radius) return 0.0; // fast test
	float const dist_sq(p2p_dist_sq(pos, center));
	if (dist_sq > radius*radius)         return 0.0;
	float const rscale((radius - sqrt(dist_sq))*radius_inv);
	return rscale*rscale*color[3]; // quadratic 1/r^2 attenuation
}


bool light_source::lights_polygon(point const &pc, float rsize, vector3d const* const norm) const {
	
	if (norm && dot_product_ptv(*norm, center, pc) <= 0.0) return 0;
	return (radius == 0.0 || dist_less_than(pc, center, (radius + rsize)));
}


float light_source::get_dir_intensity(vector3d const &obj_dir) const {

	if (bwidth == 1.0) return 1.0;
	float const dp(dot_product(obj_dir, dir));
	if (dp >= 0.0 && (bwidth + LT_DIR_FALLOFF) < 0.5) return 0.0;
	float const dp_norm(0.5*(-dp*InvSqrt(obj_dir.mag_sq()) + 1.0)); // dp = -1.0 to 1.0, bw = 0.0 to 1.0
	return CLIP_TO_01(2.0f*(dp_norm + bwidth + LT_DIR_FALLOFF - 1.0f)*LT_DIR_FALLOFF_INV);
}


void light_source::get_bounds(point bounds[2], int bnds[3][2], float thresh) const {

	if (radius == 0.0) {
		for (unsigned d = 0; d < 3; ++d) {
			bounds[0][d] = -SCENE_SIZE[d];
			bounds[1][d] =  SCENE_SIZE[d];
			bnds[d][0]   = 0;
			bnds[d][1]   = MESH_SIZE[d]-1;
		}
	}
	else {
		float const rb(radius*(1.0 - sqrt(thresh)));

		for (unsigned d = 0; d < 3; ++d) {
			for (unsigned j = 0; j < 2; ++j) {
				bounds[j][d] = center[d] + (j ? rb : -rb);
				bnds[d][j]   = max(0, min(MESH_SIZE[d]-1, get_dim_pos(bounds[j][d], d)));
			}
		}
	}
}


bool light_source::is_visible() const {

	return (radius == 0.0 || sphere_in_camera_view(center, radius, 0)); // max_level?
}


void light_source::shift_by(vector3d const &vd) {
	
	center += vd;

	for (unsigned i = 0; i < 3; ++i) { // z-shift may not be correct
		cent[i] = max(0, min(MESH_SIZE[i]-1, (cent[i] + int((2.0*SCENE_SIZE[i]/MESH_SIZE[i])*vd[i])))); // clamp to mesh bounds
	}
}


void light_source::combine_with(light_source const &l) {

	assert(radius > 0.0);
	float const w1(radius*radius*radius), w2(l.radius*l.radius*l.radius), wsum(w1 + w2), wa(w1/wsum), wb(w2/wsum);
	radius     = pow(wsum, (1.0f/3.0f));
	radius_inv = 1.0/radius;
	center    *= wa;
	center    += l.center*wb; // weighted average
	blend_color(color, color, l.color, wa, 1);
	calc_cent();
}


void light_source::draw(int ndiv) const {

	if (radius == 0.0) return;
	set_color(color);
	draw_sphere_at(center, 0.05*radius, ndiv);
}


void shift_light_sources(vector3d const &vd) {

	for (unsigned i = 0; i < light_sources.size(); ++i) {
		light_sources[i].shift_by(vd);
	}
}


void dls_cell::get_close_sources(point const &pos, float radius, vector<unsigned> &dlights) const {

	unsigned const lsz(lsrc.size());

	for (unsigned l = 0; l < lsz; ++l) { // slow, unfinished
		light_source &ls(dl_sources[lsrc[l]]);
		//if (!ls.check_counter()) continue;
		if (!dist_less_than(pos, ls.get_center(), (radius + ls.get_radius()))) continue;
		dlights.push_back(lsrc[l]);
	}
}


// *** R_PROFILE IMPLEMENTATION ***


void r_profile::reset_bbox(float const bb_[2][2]) {
	
	clear();
	bb       = rect(bb_);
	tot_area = bb.area();
	assert(tot_area > 0.0);
}


void r_profile::clear() {
	
	rects.resize(0);
	filled    = 0;
	avg_alpha = 1.0;
}


void r_profile::add_rect_int(rect const &r) {

	for (unsigned i = 0; i < rects.size(); ++i) { // try rect merge
		if (rects[i].merge_with(r)) return;
	}
	rects.push_back(r);
}


bool r_profile::add_rect(float const d[3][2], unsigned d0, unsigned d1, float alpha=1.0) {
	
	if (filled || alpha == 0.0) return 0;
	rect r(d, d0, d1);
	if (!r.nonzero() || !r.overlaps(bb.d)) return 0;
	r.clip_to(bb.d);
	//if (r.is_near_zero_area())  return 0;
	
	if (r.equal(bb.d)) { // full containment
		if (alpha < 1.0) {
			// *** WRITE ***
		}
		else {
			avg_alpha = 1.0;
		}
		rects.resize(0);
		add_rect_int(r);
		filled = 1;
		return 1;
	}
	if (rects.empty()) { // single rect performance optimization
		add_rect_int(r);
		avg_alpha = alpha;
		return 1;
	}
	unsigned const nrects(rects.size());

	for (unsigned i = 0; i < nrects; ++i) { // check if contained in any rect
		if (rects[i].contains(r.d)) return 1; // minor performance improvement
	}
	pend.push_back(r);

	while (!pend.empty()) { // merge new rect into working set while removing overlaps
		bool bad_rect(0);
		rect rr(pend.front());
		pend.pop_front();

		for (unsigned i = 0; i < nrects; ++i) { // could start i at the value of i where rr was inserted into pend
			if (rects[i].overlaps(rr.d)) { // split rr
				rects[i].subtract_from(rr, pend);
				bad_rect = 1;
				break;
			}
		}
		if (!bad_rect) add_rect_int(rr);
	}
	if (rects.size() > nrects) avg_alpha = 1.0; // at least one rect was added, *** FIX ***
	return 1;
}


float r_profile::clipped_den_inv(float const c[2]) const { // clip by first dimension
	
	if (filled)        return (1.0 - avg_alpha);
	if (rects.empty()) return 1.0;
	bool const no_clip(c[0] == bb.d[0][0] && c[1] == bb.d[0][1]);
	unsigned const nrects(rects.size());
	float a(0.0);

	if (no_clip) {
		for (unsigned i = 0; i < nrects; ++i) {a += rects[i].area();}
	}
	else {
		for (unsigned i = 0; i < nrects; ++i) {a += rects[i].clipped_area(c);}
	}
	if (a == 0.0) return 1.0;
	a *= avg_alpha;
	float const area(no_clip ? tot_area : (c[1] - c[0])*(bb.d[1][1] - bb.d[1][0]));
	if (a > area + TOLER) cout << "a = " << a << ", area = " << area << ", size = " << rects.size() << endl;
	assert(a <= area + TOLER);
	return (area - a)/area;
}


void r_profile::clear_within(float const c[2]) {

	float const rd[2][2] = {{c[0], c[1]}, {bb.d[1][0], bb.d[1][1]}};
	rect const r(rd);
	bool removed(0);

	for (unsigned i = 0; i < rects.size(); ++i) {
		if (!r.overlaps(rects[i].d)) continue;
		r.subtract_from(rects[i], pend);
		swap(rects[i], rects.back());
		rects.pop_back();
		removed = 1;
		--i; // wraparound OK
	}
	copy(pend.begin(), pend.end(), back_inserter(rects));
	pend.resize(0);
	if (removed) filled = 0;
	//avg_alpha = 1.0; // *** FIX ***
}


// *** MAIN LIGHTMAP CODE ***


void reset_cobj_counters() {

	unsigned const ncobjs(coll_objects.size());

	for (unsigned i = 0; i < ncobjs; ++i) {
		coll_objects[i].counter = -1;
	}
}


void reset_flow_cache() {

	for (unsigned i = 0; i < FLOW_CACHE_SZ; ++i) {
		flow_cache[i].reset();
	}
}


lmcell *lmap_manager_t::get_lmcell(point const &p) {

	int const x(get_xpos(p.x - SHIFT_DX)), y(get_ypos(p.y - SHIFT_DY)), z(get_zpos(p.z));
	return (is_valid_cell(x, y, z) ? &vlmap[y][x][z] : NULL);
}


void lmap_manager_t::alloc(unsigned nbins, unsigned zsize, unsigned char **need_lmcell) {

	assert(need_lmcell != NULL);
	if (vlmap == NULL) matrix_gen_2d(vlmap);
	lm_zsize = zsize;
	vldata_alloc.resize(nbins);
	unsigned cur_v(0);

	// initialize lightmap
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (!need_lmcell[i][j]) {
				vlmap[i][j] = NULL;
				continue;
			}
			assert(cur_v + lm_zsize <= vldata_alloc.size());
			vlmap[i][j] = &vldata_alloc[cur_v];
			cur_v      += lm_zsize;
		}
	}
	assert(cur_v == vldata_alloc.size());
}


void lmap_manager_t::normalize_light_val(float min_light, float max_light, float light_scale, float light_off) {

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			lmcell *vlm(vlmap[i][j]);
			if (vlm == NULL) continue;

			for (unsigned v = 0; v < lm_zsize; ++v) {
				vlm[v].v = max(min_light, min(max_light, (light_scale*vlm[v].v + light_off)));
			}
		}
	}
}


#define GET_DIST_FROM(A, B) {int const val(dv[A]*(from[B] - dvt[B]) - dv[B]*(from[A] - dvt[A])); dist += val*val;}


// Note: assumes input ranges are valid
// from = center, to = query point
float get_flow_val(CELL_LOC_T const from[3], CELL_LOC_T const to[3], bool use_flow_cache) { // from should be the light source

	CELL_LOC_T const dv[3] = {(to[0] - from[0]), (to[1] - from[1]), (to[2] - from[2])};
	if (dv[0] == 0 && dv[1] == 0 && dv[2] == 0) return 1.0;
	flow_cache_e ce(from, to), &cached(use_flow_cache ? (flow_cache[ce.hash() & (FLOW_CACHE_SZ-1)]) : flow_cache[0]);
	if (use_flow_cache && cached == ce) return cached.val; // check if this value is cached
	CELL_LOC_T cur[3] = {from[0], from[1], from[2]};
	float mult_flow(1.0), max_flow(1.0), val(1.0);

	for (unsigned i = 0; i < 3; ++i) {
		assert(from[i] >= 0 && from[i] < MESH_SIZE[i]);
		assert(to[i]   >= 0 && to[i]   < MESH_SIZE[i]);
	}
	while (val > CTHRESH) {
		unsigned di(0);
		int dmin(-1);
	
		for (unsigned d = 0; d < 3; ++d) {
			if (cur[d] == to[d]) continue;
			int dist(0);
			CELL_LOC_T dvt[3] = {cur[0], cur[1], cur[2]};
			dvt[d] += ((to[d] > cur[d]) ? 1 : -1);
			GET_DIST_FROM(0, 1);
			GET_DIST_FROM(1, 2);
			GET_DIST_FROM(2, 0);
			if (dmin == -1 || dist < dmin) {dmin = dist; di = d;}
		}
		if (dmin < 0) break;
		bool const positive(to[di] > cur[di]);
		if (!positive) --cur[di]; // step backwards first then calc flow
		
		if (lmap_manager.vlmap[cur[1]][cur[0]]) {
			float const fval(lmap_manager.vlmap[cur[1]][cur[0]][cur[2]].lflow[di]/255.0);
			mult_flow *= fval;
			max_flow   = min(max_flow, fval);
		}
		if ( positive) ++cur[di]; // calc flow first then step forward
		val = 0.5*(max_flow + mult_flow); // not sure which is best - max or mult, so choose their average
	}
	if (use_flow_cache) {
		ce.val = val;
		cached = ce;
	}
	return val;
}


bool has_fixed_cobjs(int x, int y) {

	bool has_fixed(0);
	assert(!point_outside_mesh(x, y));
	vector<int> const &cvals(v_collision_matrix[y][x].cvals);
	unsigned const ncv(cvals.size());

	for (unsigned k = 0; k < ncv && !has_fixed; ++k) {
		if (coll_objects[cvals[k]].fixed && coll_objects[cvals[k]].status == COLL_STATIC) has_fixed = 1;
	}
	return has_fixed;
}


void shift_lightmap(vector3d const &vd) {

	// *** FIX ***
	regen_lightmap();
}


void regen_lightmap() {

	assert(lmap_manager.vlmap != NULL);
	clear_lightmap();
	//assert(lmap_manager.vlmap == NULL);
	build_lightmap(0);
	assert(lmap_manager.vlmap != NULL);
}


void clear_lightmap() {

	if (lmap_manager.vlmap == NULL) return;
	if (using_lightmap) reset_flow_cache();
	lmap_manager.clear();
	using_lightmap = 0;
	lm_alloc       = 0;
	czmin0         = czmin;
}


void calc_flow_for_xy(r_profile flow_prof[2][3], int **z_light_depth, int i, int j,
	bool proc_cobjs, bool calc_lval, float zstep, float z_atten, float light_off)
{
	assert(zstep > 0.0);
	if (z_light_depth) z_light_depth[i][j] = MESH_SIZE[2];
	lmcell *vldata(lmap_manager.vlmap[i][j]);
	if (vldata == NULL) return;
	float const bbz[2][2] = {{get_xval(j), get_xval(j+1)}, {get_yval(i), get_yval(i+1)}}; // X x Y
	coll_cell const &cell(v_collision_matrix[i][j]);
	unsigned const ncv(cell.cvals.size());
	float val(1.0), vscale(Z_WT_SCALE);
	r_profile prof(bbz);
	vector<pair<float, unsigned> > cobj_z;
	bool alpha1(1);

	if (proc_cobjs) {
		for (unsigned q = 0; q < ncv; ++q) {
			unsigned const cid(cell.cvals[q]);
			assert(cid < coll_objects.size());
			coll_obj const &cobj(coll_objects[cid]);
			if (cobj.status != COLL_STATIC) continue;
			if (cobj.type == COLL_CYLINDER_ROT && !line_is_axis_aligned(cobj.points[0], cobj.points[1])) continue; // bounding cube is too conservative, skip
			if (cobj.d[2][1] < zbottom)     continue; // below the mesh
			rect const r_cobj(cobj.d, 0, 1);
			if (!r_cobj.nonzero())          continue; // zero Z cross section (vertical polygon)
			float cztop;
					
			if (r_cobj.overlaps(bbz) && add_cobj_ok(cobj) && cobj.clip_in_2d(bbz, cztop, 0, 1, 1)) {
				cobj_z.push_back(make_pair(cztop, cid)); // still incorrect for coll polygon since x and y aren't clipped
				if (cobj.cp.color.alpha < 1.0) alpha1 = 0;
			}
		}
		sort(cobj_z.begin(), cobj_z.end(), std::greater<pair<float, unsigned> >()); // max to min z
	}
	unsigned const ncv2(cobj_z.size());
	unsigned c(0);

	for (int v = MESH_SIZE[2]-1; v >= 0; --v) { // top to bottom
		float const old_val(val);
		float zb(czmin0 + v*zstep), zt(zb + zstep); // cell Z bounds
		float flow_val[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
				
		if (Z_WT_SCALE == 0.0 || zt < mesh_height[i][j]) { // under mesh
			val = 0.0;
		}
		else if (!proc_cobjs) {
			for (unsigned d = 0; d < 3; ++d) {
				flow_val[0][d] = flow_val[1][d] = 1.0;
			}
			val = 1.0;
		}
		else {
			if (val > 0.0) {
				for (; val > 0.0 && c < ncv2 && cobj_z[c].first >= zb; ++c) {
					coll_obj const &cobj(coll_objects[cobj_z[c].second]);

					if (cobj.cp.color.alpha > MIN_SHADOW_ALPHA && cobj_z[c].first < zt) { // test for ztop in current voxel //(v == MESH_SIZE[2]-1) ||
						prof.add_rect(cobj.d, 0, 1, cobj.cp.color.alpha);
						if (prof.is_filled()) val = 0.0;
					}
				}
				float const new_val(prof.den_inv());
				assert(new_val > -TOLER && new_val <= (val + TOLER));
						
				if (new_val == 1.0) {
					if (z_light_depth) z_light_depth[i][j] = v;
					vscale = Z_WT_SCALE;
				}
				val = new_val;
			} // if val > 0.0
			float const bb[3][2]  = {{bbz[0][0], bbz[0][1]}, {bbz[1][0], bbz[1][1]}, {zb, zt}};
			float const bbx[2][2] = {{bb[1][0], bb[1][1]}, {zb, zt}}; // YxZ
			float const bby[2][2] = {{zb, zt}, {bb[0][0], bb[0][1]}}; // ZxX

			for (unsigned d = 0; d < unsigned(1 + !alpha1); ++d) {
				flow_prof[d][0].reset_bbox(bbx);
				flow_prof[d][1].reset_bbox(bby);
				flow_prof[d][2].reset_bbox(bbz);
			}
			for (unsigned c2 = 0; c2 < ncv2; ++c2) { // could make this more efficient
				coll_obj &cobj(coll_objects[cobj_z[c2].second]);
				if (cobj.d[0][0] >= bb[0][1] || cobj.d[0][1]     <= bb[0][0]) continue; // no intersection
				if (cobj.d[1][0] >= bb[1][1] || cobj.d[1][1]     <= bb[1][0]) continue;
				if (cobj.d[2][0] >= bb[2][1] || cobj_z[c2].first <= bb[2][0]) continue;
				float const cztop(cobj.d[2][1]), alpha(cobj.cp.color.alpha);
				cobj.d[2][1] = cobj_z[c2].first;
						
				for (unsigned d = 0; d < 3; ++d) { // critical path
					if (alpha > MIN_SHADOW_ALPHA) flow_prof[0][d].add_rect(cobj.d, (d+1)%3, (d+2)%3, alpha);
					if (!alpha1)                  flow_prof[1][d].add_rect(cobj.d, (d+1)%3, (d+2)%3, 1.0);
				}
				cobj.d[2][1] = cztop;
			} // for c2
			for (unsigned d = 0; d < unsigned(1 + !alpha1); ++d) {
				for (unsigned e = 0; e < 3; ++e) {
					flow_val[d][e] = flow_prof[d][e].den_inv();
					assert(flow_val[d][e] > -TOLER);
				}
			}
		} // if above mesh
		for (unsigned d = 0; d < 3; ++d) {
			vldata[v].lflow[d] = (unsigned char)(255.5*CLIP_TO_01(flow_val[0][d]));
			vldata[v].pflow[d] = (unsigned char)(255.5*CLIP_TO_01(flow_val[!alpha1][d]));
		}
		if (calc_lval) {
			vldata[v].v = 0.5*vscale*(val + old_val) + light_off; // what about colors?
			vscale     *= z_atten;
			UNROLL_3X(vldata[v].ac[i_] = 1.0;)
		}
	} // for v
}


float calc_czspan() {return max(0.0f, ((czmax + lm_dz_adj) - czmin0 + TOLER));}


void build_lightmap(bool verbose) {

	if (lm_alloc) return; // what about recreating the lightmap if the scene has changed?
	if (verbose) cout << "Building lightmap" << endl;
	RESET_TIME;
	unsigned nonempty(0);
	unsigned char **need_lmcell = NULL;
	matrix_gen_2d(need_lmcell);
	bool has_fixed(0);
	
	// determine where we will need lmcells
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			bool const fixed(!coll_objects.empty() && has_fixed_cobjs(j, i));
			need_lmcell[i][j] = (use_dense_voxels || fixed);
			has_fixed        |= fixed;
			if (need_lmcell[i][j]) ++nonempty;
		}
	}

	// add cells surrounding static scene lights
	// Note: this isn't really necessary when using ray casting for lighting,
	//       but it helps ensure there are lmap cells around light sources to light the dynamic objects
	for (unsigned i = 0; i < light_sources.size(); ++i) {
		point bounds[2];
		int bnds[3][2];
		light_sources[i].get_bounds(bounds, bnds, CTHRESH);

		for (int y = bnds[1][0]; y <= bnds[1][1]; ++y) {
			for (int x = bnds[0][0]; x <= bnds[0][1]; ++x) {
				if (!need_lmcell[y][x]) ++nonempty;
				need_lmcell[y][x] |= 2;
			}
		}
	}

	// prevent the z range from being empty/denormalized when there are no cobjs
	if (use_dense_voxels) {
		czmin = min(czmin, zbottom);
		czmax = max(czmax, (czmin + Z_SCENE_SIZE - 0.5f*DZ_VAL));
	}
	else if (czmin >= czmax) {
		czmin = min(czmin, zbottom);
		czmax = max(czmax, ztop);
	}

	// determine allocation and voxel grid sizes
	reset_cobj_counters();
	assert(DZ_VAL > 0.0 && Z_LT_ATTEN > 0.0 && Z_LT_ATTEN <= 1.0 && XY_LT_ATTEN > 0.0 && XY_LT_ATTEN <= 1.0);
	DZ_VAL_INV2   = DZ_VAL_SCALE/DZ_VAL;
	SHIFT_DXYZ[0] = SHIFT_DX = SHIFT_VAL*DX_VAL;
	SHIFT_DXYZ[1] = SHIFT_DY = SHIFT_VAL*DY_VAL;
	czmin0        = czmin;//max(czmin, zbottom);
	assert(lm_dz_adj >= 0.0);
	float const czspan(calc_czspan()), dz(DZ_VAL_INV2*czspan);
	assert(dz >= 0.0);
	assert(coll_objects.empty() || !has_fixed || dz > 0.0); // too strict (all cobjs can be shifted off the mesh)
	unsigned zsize(unsigned(dz + 1));
	
	if ((int)zsize > MESH_Z_SIZE) {
		cout << "* Warning: Scene height extends beyond the specified z range. Clamping zsize of " << zsize << " to " << MESH_Z_SIZE << "." << endl;
		zsize = MESH_Z_SIZE;
	}
	unsigned const nbins(nonempty*zsize);
	MESH_SIZE[2] = zsize; // override MESH_SIZE[2]
	float const zstep(czspan/zsize), scene_scale(MESH_X_SIZE/128.0);
	if (verbose) cout << "zsize= " << zsize << ", nonempty= " << nonempty << ", bins= " << nbins << ", czmin= " << czmin0 << ", czmax= " << czmax << endl;
	assert(zstep > 0.0);
	float const z_atten(1.0 - (1.0 - Z_LT_ATTEN)/scene_scale), xy_atten(1.0 - (1.0 - XY_LT_ATTEN)/scene_scale);
	if (!ldynamic) matrix_gen_2d(ldynamic, MESH_X_SIZE, MESH_Y_SIZE);
	lmap_manager.alloc(nbins, zsize, need_lmcell);
	using_lightmap = (nonempty > 0);
	lm_alloc       = 1;
	assert(ldynamic && lmap_manager.vlmap);
	int **z_light_depth = NULL;
	matrix_gen_2d(z_light_depth);
	if (verbose) PRINT_TIME(" Lighting Setup");
	bool raytrace_lights[3];
	UNROLL_3X(raytrace_lights[i_] = (read_light_files[i_] || write_light_files[i_]););
	bool const no_comp_light(raytrace_lights[LIGHTING_SKY]);
	float const light_off(no_comp_light ? 0.0f : LIGHT_OFFSET);

	// process vertical (Z) light projections
	r_profile flow_prof[2][3]; // {light, particle} x {x, y, z}

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			bool const proc_cobjs(need_lmcell[i][j] & 1);
			calc_flow_for_xy(flow_prof, z_light_depth, i, j, proc_cobjs, !no_comp_light, zstep, z_atten, light_off);
		} // for j
	} // for i
	if (verbose) PRINT_TIME(" Lighting Z + Flow");
	int const bnds[2][2] = {{0, MESH_X_SIZE-1}, {0, MESH_Y_SIZE-1}};
	float const fbnds[2] = {(X_SCENE_SIZE - TOLER), (Y_SCENE_SIZE - TOLER)};
	int counter(0);
	float pass_weight(1.0);

	// process lateral (X and Y) light projections
	for (unsigned pass = 0; pass < NUM_XY_PASSES; ++pass) {
		if (XY_WT_SCALE == 0.0 || nbins == 0 || dz == 0.0 || no_comp_light) continue;

		for (unsigned dim = 0; dim < 2; ++dim) { // x/y
			for (unsigned dir = 0; dir < 2; ++dir) { // +/-
				for (int s = bnds[dim][0]; s < bnds[dim][1]; ++s) { // outer iteration (direction doesn't matter)
					float bb[2][2] = {{czmin0, (czmax + lm_dz_adj)}, {0, 0}}; // {z, x/y}

					for (unsigned d = 0; d < 2; ++d) {
						bb[1][d] = get_dim_val((s+d), dim);
					}
					r_profile prof(bb);
					vector<float> vals(zsize, 1.0); // start at full light
					vector<float> vscale(zsize, XY_WT_SCALE*pass_weight);
					vector<float> mesh_light(zsize, 1.0);
					float const tmax(fbnds[!dim]);
					int const dt(dir ? -1 : 1);
					++counter;

					for (int t = bnds[!dim][dir]; t != (bnds[!dim][!dir] + dt); t += dt) { // inner iteration
						int i[2] = {(dim ? t : s), (dim ? s : t)}; // dim = 0 => {s, t} = {x, y}
						assert(!point_outside_mesh(i[0], i[1]));
						coll_cell const &cell(v_collision_matrix[i[1]][i[0]]);
						unsigned const ncv(cell.cvals.size());
						lmcell *vlm(lmap_manager.vlmap[i[1]][i[0]]);

						if (vlm == NULL) { // no cobjs
							for (int v = zsize-1; v >= 0; --v) {
								reset_lighting(mesh_light[v], vals[v], vscale[v]);
							}
							prof.clear();
							continue;
						}
						if (ncv > 0) {
							float t_bnds[2];
							
							for (unsigned d = 0; d < 2; ++d) {
								t_bnds[d] = get_dim_val((t+d), !dim); // bounds of current t row
							}
							for (unsigned c = 0; c < ncv; ++c) { // skip entire column somehow?
								coll_obj &cobj(coll_objects[cell.cvals[c]]);
								if (cobj.status != COLL_STATIC) continue;
								if (cobj.type == COLL_CUBE && cobj.counter == counter) continue; // already processed
								if (cobj.d[2][1] < zbottom || cobj.cp.color.alpha <= MIN_SHADOW_ALPHA || !add_cobj_ok(cobj)) continue;
								float xyval(cobj.d[!dim][dir]);
								//if (!cobj.clip_in_2d(bb, xyval, 2, dim, dir)) continue;
								float const tv(min(tmax, max(-tmax, xyval))); // clip to scene

								if ((tv < t_bnds[1] && tv >= t_bnds[0]) || (POLY_XY_OVER_CHK && cobj.type == COLL_POLYGON &&
									cobj.d[!dim][0] < t_bnds[1] && cobj.d[!dim][1] >= t_bnds[0]))
								{ // test for intersection or polygon bbox overlap with current column
									cobj.counter = counter;
									prof.add_rect(cobj.d, 2, dim, cobj.cp.color.alpha); // somewhat critical path
									if (prof.is_filled()) break;
								}
							} // for c
						} // if ncv > 0
						//assert(vlm <= &vldata_alloc[lmap_manager.size() - zsize]);
						float const m_height(mesh_height[i[1]][i[0]]);
						int const zl_depth(z_light_depth[i[1]][i[0]]);

						if (zl_depth < int(zsize)) {
							float const clipval[2] = {(czmin0 + zl_depth*zstep), (czmin0 + zsize*zstep)};
							prof.clear_within(clipval);
						}
						for (int v = zsize-1; v >= 0; --v) {
							float const zb(czmin0 + v*zstep), zt(zb + zstep); // cell Z bounds
							float old_val(vals[v]);

							if (zl_depth <= v) { // lit from above
								reset_lighting(mesh_light[v], vals[v], vscale[v]);
								old_val = vals[v];
							}
							else if (vals[v] > 0.0) { // not already zero
								float const clipval[2] = {zb, zt};
								float const new_val(prof.clipped_den_inv(clipval));
								assert(new_val > -TOLER);
								vals[v] = min(old_val, new_val); // can't increase due to density
							}
							if (old_val > 0.0) {
								if (zt < m_height) { // under mesh
									mesh_light[v] = 0.0; // reset mesh light
									vals[v]       = 0.0; // can't get lit if under mesh
								}
								vlm[v].v = min(1.0f, (vlm[v].v + vscale[v]*0.5f*(vals[v] + old_val))); // merge XY and Z values
							}
							vscale[v] *= xy_atten;
						} // for v
					} // for t
				} // for s
			} // for dir
		} // for dim
		pass_weight *= PASS_WEIGHT_ATT;
	} // for pass
	if (verbose) PRINT_TIME(" Lighting XY");

	// add in static light sources
	if (!raytrace_lights[LIGHTING_LOCAL]) {
		for (unsigned i = 0; i < light_sources.size(); ++i) {
			light_source &ls(light_sources[i]);
			point const &lpos(ls.get_center());
			if (!is_over_mesh(lpos)) continue;
			ls.calc_cent();
			colorRGBA const &lcolor(ls.get_color());
			point bounds[2];
			int bnds[3][2], cobj(-1), last_cobj(-1);
			CELL_LOC_T const *const cent(ls.get_cent());
			ls.get_bounds(bounds, bnds, CTHRESH);
			if (SLT_LINE_TEST_WT > 0.0) check_coll_line(lpos, lpos, cobj, -1, 1, 2); // check cobj containment and ignore that shape

			for (int y = bnds[1][0]; y <= bnds[1][1]; ++y) {
				for (int x = bnds[0][0]; x <= bnds[0][1]; ++x) {
					assert(lmap_manager.vlmap[y][x]);
					float const xv(get_xval(x)), yv(get_yval(y));

					for (int z = bnds[2][0]; z <= bnds[2][1]; ++z) {
						assert(unsigned(z) < zsize);
						point const p(xv, yv, get_zval(z));
						float cscale(ls.get_intensity_at(p));
						if (cscale < CTHRESH) {if (z > cent[2]) break; else continue;}
						CELL_LOC_T const cur_loc[3] = {x, y, z};
				
						if (ls.is_directional()) {
							vector3d const vlp(lpos, p);
							cscale *= ls.get_dir_intensity(vlp);
							if (cscale < CTHRESH) continue;
						}
						float flow[2] = {1.0, 1.0};

						if (SLT_LINE_TEST_WT > 0.0) { // slow
							if ((last_cobj >= 0 && coll_objects[last_cobj].line_intersect(lpos, p)) ||
								check_coll_line(p, lpos, last_cobj, cobj, 1, 2)) flow[0] = 0.0;
						}
						if (SLT_FLOW_TEST_WT > 0.0) flow[1] = get_flow_val(cent, cur_loc, 0); // determine flow value to this lmcell
						cscale *= SLT_LINE_TEST_WT*flow[0] + SLT_FLOW_TEST_WT*flow[1];
						if (cscale < CTHRESH) continue;
						lmcell &lmc(lmap_manager.vlmap[y][x][z]);
						UNROLL_3X(lmc.c[i_] = min(1.0f, (lmc.c[i_] + cscale*lcolor[i_]));) // what about diffuse/normals?
					} // for z
				} // for x
			} // for y
		} // for i
		if (verbose) PRINT_TIME(" Light Source Addition");
	}
	float const lscales[4] = {1.0/SQRT3, 1.0/SQRT2, 1.0, 0.0};

	// smoothing passes
	if (LIGHT_SPREAD > 0.0 && (!raytrace_lights[LIGHTING_SKY] || !raytrace_lights[LIGHTING_LOCAL])) {
		for (unsigned n = 0; n < NUM_LT_SMOOTH; ++n) { // low pass filter light, bleed to adjacent cells
			for (int i = 0; i < MESH_Y_SIZE; ++i) {
				for (int j = 0; j < MESH_X_SIZE; ++j) {
					lmcell *vlm(lmap_manager.vlmap[i][j]);
					if (vlm == NULL) continue;

					for (int v = 0; v < (int)zsize; ++v) {
						float vlmv(vlm[v].v);
						float const *const color(vlm[v].c);
						if (vlmv == 0.0 && color[0] == 0.0 && color[1] == 0.0 && color[2] == 0.0) continue;
						vlmv *= LIGHT_SPREAD*scene_scale;

						for (int k = max(0, i-1); k < min(MESH_Y_SIZE-1, i+1); ++k) {
							for (int l = max(0, j-1); l < min(MESH_X_SIZE-1, j+1); ++l) {
								if (lmap_manager.vlmap[k][l] == NULL) continue;

								for (int m = max(0, v-1); m < min((int)zsize-1, v+1); ++m) {
									unsigned const nsame((k == i) + (l == j) + (m == v));
									if (nsame != 2) continue; // for now, to make flow calculations simpler
									unsigned const dim((k!=i) ? 1 : ((l!=j) ? 0 : 2)); // x, y, z
									unsigned const ix((l>j)?(j+1):j), iy((k>i)?(i+1):i), iz((m>v)?(v+1):v);
									bool const oob((int)ix >= MESH_X_SIZE || (int)iy >= MESH_Y_SIZE || iz >= zsize); // out of bounds
									float flow(oob ? 1.0 : lmap_manager.vlmap[iy][ix][iz].lflow[dim]/255.0);
									if (flow <= 0.0) continue;
									lmcell &lm(lmap_manager.vlmap[k][l][m]);
									flow *= lscales[nsame];
									lm.v += flow*vlmv;
									UNROLL_3X(lm.c[i_] = min(1.0f, (lm.c[i_] + flow*color[i_]));) // assumes color.alpha is always either 0.0 or 1.0
								} // for m
							} // for l
						} // for k
					} // for v
				} // for j
			} // for i
		} // for n
		if (verbose) PRINT_TIME(" Lighting Smooth");
	} // if LIGHT_SPREAD

	if (nbins > 0) {
		if (raytrace_lights[LIGHTING_SKY]) {
			compute_ray_trace_lighting_sky();
			if (verbose) PRINT_TIME(" Sky Lightmap Ray Trace");
		}
		if (raytrace_lights[LIGHTING_GLOBAL]) {
			compute_ray_trace_lighting_global();
			if (verbose) PRINT_TIME(" Global Lightmap Ray Trace");
		}
		if (raytrace_lights[LIGHTING_LOCAL]) {
			compute_ray_trace_lighting_local();
			if (verbose) PRINT_TIME(" Local Lightmap Ray Trace");
		}
	}
	// normalize final light value to [MIN_LIGHT, MAX_LIGHT]
	lmap_manager.normalize_light_val(MIN_LIGHT, MAX_LIGHT, LIGHT_SCALE, light_off);
	reset_cobj_counters();
	matrix_delete_2d(z_light_depth);
	matrix_delete_2d(need_lmcell);
	PRINT_TIME(" Lighting");
}


void update_flow_for_voxels(cube_t const &cube) {

	if (!lm_alloc || !lmap_manager.vlmap) return;
	int const cx1(max(0, get_xpos(cube.d[0][0]))), cx2(min(MESH_X_SIZE-1, get_xpos(cube.d[0][1])));
	int const cy1(max(0, get_ypos(cube.d[1][0]))), cy2(min(MESH_Y_SIZE-1, get_ypos(cube.d[1][1])));
	float const zstep(calc_czspan()/MESH_SIZE[2]);
	r_profile flow_prof[2][3];

	for (int y = cy1; y <= cy2; ++y) {
		for (int x = cx1; x <= cx2; ++x) {
			assert(!point_outside_mesh(x, y));
			bool const fixed(!coll_objects.empty() && has_fixed_cobjs(x, y));
			bool const proc_cobjs(use_dense_voxels || fixed);
			calc_flow_for_xy(flow_prof, NULL, y, x, proc_cobjs, 0, zstep, 1.0, 0.0);
		} // for x
	} //for y
}


// *** Dynamic Lights Code ***


void setup_2d_texture(unsigned &tid) {

	assert(tid == 0);
	glGenTextures(1, &tid);
	bind_2d_texture(tid);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}


unsigned upload_voxel_flow_texture() {

	unsigned const zsize(MESH_SIZE[2]), sz(MESH_X_SIZE*MESH_Y_SIZE*zsize), ncomp(3);
	vector<unsigned char> data(ncomp*sz, 255); // start at max flow

	for (int y = 0; y < MESH_Y_SIZE; ++y) {
		for (int x = 0; x < MESH_X_SIZE; ++x) {
			lmcell const *const vlm(lmap_manager.vlmap[y][x]);
			if (vlm == NULL) continue;
			unsigned const off(zsize*(y*MESH_X_SIZE + x));

			for (unsigned z = 0; z < zsize; ++z) {
				unsigned const off2(ncomp*(off + z));
				UNROLL_3X(data[off2+i_] = vlm[z].lflow[i_];)
			}
		}
	}
	return create_3d_texture(zsize, MESH_X_SIZE, MESH_Y_SIZE, ncomp, data, GL_LINEAR);
}


// texture units used:
// 0: reserved for object textures
// 1: reserved for indirect sky lighting and smoke (if enabled)
// 2: dynamic light data
// 3: dynamic light element array
// 4: dynamic light grid bag
// 5: voxel flow (not yet enabled) / reserved for bump maps
// 6: reserved for shadow map sun
// 7: reserved for shadow map moon
// 8: reserved for specular maps
void upload_dlights_textures() {

	RESET_TIME;
	assert(lm_alloc && lmap_manager.vlmap);
	if (disable_shaders) return;
	static int supports_tex_int(2); // starts at unknown
	static bool last_dlights_empty(0);
	bool const cur_dlights_empty(dl_sources.empty());
	if (cur_dlights_empty && last_dlights_empty && dl_tid != 0 && elem_tid != 0 && gb_tid != 0) return; // no updates
	last_dlights_empty = cur_dlights_empty;
	
	if (supports_tex_int == 2) {
		supports_tex_int = has_extension("GL_EXT_texture_integer");
		if (!supports_tex_int) cout << "Error: GL_EXT_texture_integer extension not supported. Dynamic lighting will not work correctly." << endl;
	}
	if (!supports_tex_int) {
		dl_tid = elem_tid = gb_tid = 0; // should already be 0
		return;
	}

	// step 1: the light sources themselves
	unsigned const max_dlights      = 1024; // must agree with value in shader
	unsigned const floats_per_light = 12;
	float dl_data[max_dlights*floats_per_light] = {0.0};
	unsigned const ndl(min(max_dlights, dl_sources.size()));
	unsigned const ysz(floats_per_light/4);
	float const radius_scale(1.0/X_SCENE_SIZE);
	vector3d const poff(-X_SCENE_SIZE, -Y_SCENE_SIZE, get_zval(0));
	vector3d const pscale(0.5/X_SCENE_SIZE, 0.5/Y_SCENE_SIZE, 1.0/(get_zval(MESH_SIZE[2]) - poff.z));
	has_dir_lights = 0;

	for (unsigned i = 0; i < ndl; ++i) {
		float *data(dl_data + i*floats_per_light);
		dl_sources[i].pack_to_floatv(data); // {center,radius, color, dir,beamwidth}
		UNROLL_3X(data[i_] = (data[i_] - poff[i_])*pscale[i_];) // scale to [0,1] range
		UNROLL_3X(data[i_+4] *= 0.1;) // scale color down
		data[3] *= radius_scale;
		has_dir_lights |= dl_sources[i].is_directional();
	}
	if (dl_tid == 0) {
		setup_2d_texture(dl_tid);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16, ysz, max_dlights, 0, GL_RGBA, GL_FLOAT, dl_data); // 2 x M
	}
	else {
		bind_2d_texture(dl_tid);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ysz, ndl, GL_RGBA, GL_FLOAT, dl_data);
	}

	// step 2: grid bag entries
	vector<unsigned> gb_data(XY_MULT_SIZE, 0);
	unsigned const elem_tex_sz = 256; // must agree with value in shader
	unsigned const max_gb_entries(elem_tex_sz*elem_tex_sz);
	unsigned short elem_data[max_gb_entries] = {0};
	unsigned elix(0);

	for (int y = 0; y < MESH_Y_SIZE && elix < max_gb_entries; ++y) {
		for (int x = 0; x < MESH_X_SIZE && elix < max_gb_entries; ++x) {
			unsigned const gb_ix(x + y*MESH_X_SIZE); // {start, end, unused}
			gb_data[gb_ix] = elix; // start_ix
			vector<unsigned short> const &ixs(ldynamic[y][x].get_src_ixs());
			unsigned const num_ixs(min(ixs.size(), 256U)); // max of 256 lights per bin
			
			for (unsigned i = 0; i < num_ixs && elix < max_gb_entries; ++i) { // end if exceed max entries
				if (ixs[i] >= ndl) continue; // dlight index is too high, skip
				elem_data[elix++] = (unsigned short)ixs[i];
			}
			gb_data[gb_ix] += (elix << 16); // end_ix
		}
	}
	if (elem_tid == 0) {
		setup_2d_texture(elem_tid);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE16UI_EXT, elem_tex_sz, elem_tex_sz, 0, GL_LUMINANCE_INTEGER_EXT, GL_UNSIGNED_SHORT, elem_data);
	}
	else {
		bind_2d_texture(elem_tid);
		unsigned const height(min(elem_tex_sz, (elix/elem_tex_sz+1))); // approximate ceiling
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, elem_tex_sz, height, GL_LUMINANCE_INTEGER_EXT, GL_UNSIGNED_SHORT, elem_data);
	}

	// step 3: grid bag(s)
	if (gb_tid == 0) {
		setup_2d_texture(gb_tid);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE32UI_EXT, MESH_X_SIZE, MESH_Y_SIZE, 0, GL_LUMINANCE_INTEGER_EXT, GL_UNSIGNED_INT, &gb_data.front()); // Nx x Ny
	}
	else {
		bind_2d_texture(gb_tid);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, MESH_X_SIZE, MESH_Y_SIZE, GL_LUMINANCE_INTEGER_EXT, GL_UNSIGNED_INT, &gb_data.front());
	}

	// step 4: voxel flow
#if 0
	set_multitex(5); // texture unit 5

	if (flow_tid == 0) {
		flow_tid = upload_voxel_flow_texture();
	}
	else { // no dynamic updates
		bind_3d_texture(flow_tid);
	}
#endif
	//PRINT_TIME("Dlight Texture Upload");
	//cout << "ndl: " << ndl << ", elix: " << elix << ", gb_sz: " << XY_MULT_SIZE << endl;
}


void set_one_texture(shader_t &s, unsigned tid, unsigned tu_id, const char *const name) {

	set_multitex(tu_id); // texture unit
	bind_2d_texture(tid);
	s.add_uniform_int(name, tu_id);
}


void setup_dlight_textures(shader_t &s) {

	set_one_texture(s, dl_tid,   2, "dlight_tex");
	set_one_texture(s, elem_tid, 3, "dlelm_tex");
	set_one_texture(s, gb_tid,   4, "dlgb_tex");
	//set_one_texture(p, flow_tid, 5, "flow_tex");
}


colorRGBA gen_fire_color(float &cval, float &inten, float rate) {

	inten = max(0.6f, min(1.0f, (inten + 0.04f*rate*fticks*signed_rand_float())));
	cval  = max(0.0f, min(1.0f, (cval  + 0.02f*rate*fticks*signed_rand_float())));
	colorRGBA color(1.0, 0.9, 0.7);
	blend_color(color, color, colorRGBA(1.0, 0.6, 0.2), cval, 0);
	return color;
}


void add_camera_candlelight() {

	static float cval(0.5), inten(0.75);
	add_dynamic_light(1.5*inten, get_camera_pos(), gen_fire_color(cval, inten));
}


void add_camera_flashlight() {

	add_dynamic_light(4.0, get_camera_pos(), SUN_C, cview_dir, 0.02);
}


void add_dynamic_light(float sz, point const &p, colorRGBA const &c, vector3d const &d, float bw) {

	if (!animate2) return;
	float const sz_scale(sqrt(0.1*XY_SCENE_SIZE));
	dl_sources2.push_back(light_source(sz_scale*sz, p, c, 1, d, bw));
}


void add_line_light(point const &p1, point const &p2, colorRGBA const &color, float size, float intensity) {

	if (!animate2) return;
	point p[2] = {p1, p2};
	if (!do_line_clip_scene(p[0], p[1], zbottom, max(ztop, czmax))) return;
	vector3d const dir(p[1] - p[0]);
	float const length(dir.mag());

	for (float d = 0.0; d <= length; d += 0.5*size) {
		add_dynamic_light(size*intensity, (p[0] + dir*(d/max(length, TOLERANCE))), color);
	}
}


inline void dls_cell::clear() {

	if (lsrc.empty()) return;
	if (lsrc.capacity() > INIT_CCELL_SIZE) lsrc.clear(); else lsrc.resize(0);
	z1 =  FAR_CLIP;
	z2 = -FAR_CLIP;
}


inline void dls_cell::add_light(unsigned ix, float zmin, float zmax) {
	
	if (lsrc.capacity() == 0) lsrc.reserve(INIT_CCELL_SIZE);
	lsrc.push_back(ix);
	z1 = min(z1, zmin);
	z2 = max(z2, zmax);
}


bool dls_cell::check_add_light(unsigned ix) const {

	if (empty()) return 1;
	assert(ix < dl_sources.size());
	light_source const &ls(dl_sources[ix]);

	for (unsigned i = 0; i < lsrc.size(); ++i) {
		unsigned const ix2(lsrc[i]);
		assert(ix2 < dl_sources.size());
		assert(ix2 != ix);
		if (ls.try_merge_into(dl_sources[ix2])) return 0;
	}
	return 1;
}


bool light_source::try_merge_into(light_source &ls) const {

	if (ls.radius < radius) return 0; // shouldn't get here because of radius sort
	if (!dist_less_than(center, ls.center, 0.2*max(HALF_DXY, radius)))         return 0;
	if (ls.bwidth != bwidth || ls.r_inner != r_inner || ls.dynamic != dynamic) return 0;
	if (bwidth < 1.0 && dot_product(dir, ls.dir) < 0.95)                       return 0;
	colorRGBA lcolor(color);
	float const rr(radius/ls.radius);
	lcolor.alpha *= rr*rr; // scale by radius ratio squared
	ls.add_color(lcolor);
	return 1;
}


void add_vpls() {

	if (num_vpls == 0 || (display_mode & 0x08)) return;
	float const scene_radius(get_scene_radius());

	for (vector<cube_light_source>::const_iterator i = global_cube_lights.begin(); i != global_cube_lights.end(); ++i) {
		// use the ztop plane of the cube for vpls - FIXME: make more dynamic
		cube_t const &c(i->bounds);
		float const zval(c.d[2][1]), dx(c.d[0][1] - c.d[0][0]), dy(c.d[1][1] - c.d[1][0]), ar(dy/dx);
		unsigned const nx(max(1U, (unsigned)(sqrt(num_vpls/ar) + 0.5))), ny(max(1U, (unsigned)(num_vpls/nx)));
		float const dx_step(dx/nx), dy_step(dy/ny);
		float const light_size(2.0*(dx_step + dy_step));
		//float const light_size(12.0*HALF_DXY);
		point lpos;
	
		for (int l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
			unsigned const gl_light(GL_LIGHT0 + l);
			if (!light_valid(0xFF, l, lpos) || !glIsEnabled(gl_light)) continue;
			colorRGBA ambient, diffuse;
			glGetLightfv(gl_light, GL_AMBIENT, &ambient.red);
			glGetLightfv(gl_light, GL_DIFFUSE, &diffuse.red);
			colorRGBA base_color(ambient + diffuse);
			base_color.set_valid_color();
			
			for (unsigned y = 0; y < ny; ++y) {
				for (unsigned x = 0; x < nx; ++x) {
					point const pt((c.d[0][0] + (x + 0.5)*dx_step), (c.d[1][0] + (y + 0.5)*dy_step), zval);
					vector3d const dir((pt - lpos).get_norm());
					point const target(pt + dir*scene_radius);
					point cpos;
					vector3d cnorm;
					int cindex;
					
					if (check_coll_line_exact(lpos, target, cpos, cnorm, cindex, 0.0, -1, 1, 1, 1)) { // fast=1
						float const dp(-dot_product(cnorm, dir));
						if (dp < 0.0) continue; // can this happen?
						// add small bias to prevent fp errors in the fragment shader when the light source is on top of a fragment
						cpos += cnorm*(light_size*0.01);
						assert(cindex >= 0);
						coll_obj const &cobj(coll_objects[cindex]);
						colorRGBA color(base_color);
						color = color.modulate_with(cobj.cp.color);
						if (cobj.cp.tid >= 0) color = color.modulate_with(texture_color(cobj.cp.tid));
						add_dynamic_light(dp*light_size, cpos, color, cnorm, 0.5); // 180 degree point spotlight
						//set_color(color); draw_sphere_at(cpos, dp*light_size, N_SPHERE_DIV);
					}
				}
			}
		}
	}
}


void clear_dynamic_lights() { // slow for large lights

	//if (!animate2) return;
	if (dl_sources.empty()) return; // only clear if light pos/size has changed?
	assert(ldynamic);
	
	for (int y = 0; y < MESH_Y_SIZE; ++y) {
		for (int x = 0; x < MESH_X_SIZE; ++x) ldynamic[y][x].clear();
	}
	dl_sources.resize(0);
}


void add_dynamic_lights() {

	//RESET_TIME;
	if (!animate2) return;
	assert(ldynamic);
	clear_dynamic_lights();
	if (CAMERA_CANDLE_LT) add_camera_candlelight();
	if (CAMERA_FLASH_LT)  add_camera_flashlight();
	add_vpls();
	dl_sources.swap(dl_sources2);

	for (unsigned i = 0; i < NUM_RAND_LTS; ++i) { // add some random lights (omnidirectional)
		dl_sources.push_back(light_source(0.94, gen_rand_scene_pos(), BLUE, 1));
	}
	// Note: do we want to sort by y/x position to minimize cache misses?
	sort(dl_sources.begin(), dl_sources.end(), std::greater<light_source>()); // sort by largest to smallest radius
	unsigned const ndl(dl_sources.size());
	has_dl_sources = (ndl > 0);
	bool first(1);

	for (unsigned i = 0; i < ndl; ++i) {
		light_source &ls(dl_sources[i]);
		if (!ls.is_visible()) continue; // view culling
		point const &center(ls.get_center());
		if ((center.z - ls.get_radius()) > max(ztop, czmax)) continue; // above everything, rarely occurs
		int xcent(get_xpos(center.x)), ycent(get_ypos(center.y));
		if (!point_outside_mesh(xcent, ycent) && !ldynamic[ycent][xcent].check_add_light(i)) continue;
		point bounds[2];
		int bnds[3][2];
		unsigned const ix(i);
		ls.get_bounds(bounds, bnds, 0.0);
		
		for (unsigned j = 0; j < 3; ++j) {
			dlight_bb[j][0] = (first ? bounds[0][j] : min(dlight_bb[j][0], bounds[0][j]));
			dlight_bb[j][1] = (first ? bounds[1][j] : max(dlight_bb[j][1], bounds[1][j]));
		}
		first = 0;
		int const xsize(bnds[0][1]-bnds[0][0]), ysize(bnds[1][1]-bnds[1][0]);
		int const radius((max(xsize, ysize)>>1)+2), rsq(radius*radius);

		for (int y = bnds[1][0]; y <= bnds[1][1]; ++y) {
			int const y_sq((y-ycent)*(y-ycent));

			for (int x = bnds[0][0]; x <= bnds[0][1]; ++x) {
				if (rsq == 1 || ((x-xcent)*(x-xcent) + y_sq) <= rsq) {
					ldynamic[y][x].add_light(ix, bounds[0][2], bounds[1][2]); // could do flow clipping here?
				}
			}
		}
	}
	if (SHOW_STAT_LIGHTS) {
		for (unsigned i = 0; i < light_sources.size(); ++i) {
			light_sources[i].draw(16);
		}
	}
	if (SHOW_DYNA_LIGHTS) {
		for (unsigned i = 0; i < dl_sources.size(); ++i) {
			dl_sources[i].draw(16);
		}
	}
	//PRINT_TIME("Dynamic Light Add");
}


bool is_shadowed_lightmap(point const &p) {

	if (p.z <= czmin0)  return is_under_mesh(p);
	lmcell const *const lmc(lmap_manager.get_lmcell(p));
	return (lmc ? (lmc->v < 1.0) : 0);
}


void light_source::pack_to_floatv(float *data) const {

	// store light_source as: {center.xyz, radius}, {color.rgba}, {dir, bwidth}
	// Note: we don't really need to store the z-component of dir because we can calculate it from sqrt(1 - x*x - y*y),
	//       but doing this won't save us any texture data so it's not worth the trouble
	assert(data);
	UNROLL_3X(*(data++) = center[i_];)
	*(data++) = radius;
	UNROLL_3X(*(data++) = 0.5*(1.0 + color[i_]);) // map [-1,1] => [0,1] for negative light support
	*(data++) = color[3];
	UNROLL_3X(*(data++) = 0.5*(1.0 + dir[i_]);) // map [-1,1] to [0,1]
	*(data++) = bwidth; // [0,1]
}


inline float add_specular(point const &p, vector3d ldir, vector3d const &norm, float const *const spec) {

	if (!spec || spec[0] == 0.0) return 0.0;
	ldir.normalize();
	vector3d c2p(get_camera_pos(), p);
	c2p.normalize();
	c2p += ldir;
	float dp(dot_product(norm, c2p));
	if (dp <= 0.0) return 0.0;
	// Note: we assume exponent (spec[1]) is high, so we can skip if dp is small enough
	dp *= InvSqrt(c2p.mag_sq());
	return ((dp > 0.5) ? spec[0]*pow(dp, spec[1]) : 0.0);
}


bool is_in_darkness(point const &pos, float radius, int cobj) {

	colorRGBA c(WHITE);
	float const val(get_indir_light(c, WHITE, pos, 0, 1, NULL, NULL)); // this is faster so do it first
	if ((c.red + c.green + c.blue) > DARKNESS_THRESH) return 0;

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (is_visible_to_light_cobj(pos, l, radius, cobj, 1)) return 0;
	}
	return 1;
}


bool get_dynamic_light(int x, int y, int z, point const &p, float lightscale, float *ls,
					   vector3d const *const norm, float const *const spec)
{
	if (dl_sources.empty()) return 0;
	assert(!point_outside_mesh(x, y));
	dls_cell const &ldv(ldynamic[y][x]);
	if (!ldv.check_z(p[2])) return 0;
	unsigned const lsz(ldv.size());
	CELL_LOC_T const cl[3] = {x, y, z}; // what about SHIFT_VAL?
	bool added(0);
	unsigned index(0);

	for (unsigned l = 0; l < lsz; ++l) {
		unsigned const ls_ix(ldv.get(l));
		assert(ls_ix < dl_sources.size());
		light_source const &lsrc(dl_sources[ls_ix]);
		float cscale(lightscale*lsrc.get_intensity_at(p));
		if (cscale < CTHRESH) continue;
		bool const directional(lsrc.is_directional());
		point const &lpos(lsrc.get_center());
		
		if (norm || directional) {
			vector3d const dir(lpos, p);

			if (directional) {
				cscale *= lsrc.get_dir_intensity(dir);
				if (cscale < CTHRESH) continue;
			}
			if (norm) { // ambient + diffuse + specular lighting
				float const dp(dot_product(*norm, dir));
				if (dp <= 0.0)        continue; // back facing
				cscale *= (DLIGHT_AMBIENT + DLIGHT_DIFFUSE*dp*InvSqrt(dir.mag_sq()) + (spec ? add_specular(p, dir, *norm, spec) : 0.0));
				if (cscale < CTHRESH) continue;
			}
		}
		if (DYNAMIC_LT_FLOW && using_lightmap && z >= 0) { // slow for large lights, and somewhat inaccurate
			CELL_LOC_T const *const c(lsrc.get_cent());
			CELL_LOC_T c1[3], c2[3];
			unsigned equal(0);
			
			for (unsigned d = 0; d < 3; ++d) { // take one step towards each side in each direction
				c1[d] = max(0, min(MESH_SIZE[d]-1, int(cl[d])));
				c2[d] = max(0, min(MESH_SIZE[d]-1, int(c [d])));
				if (c1[d] > c2[d]) --c1[d]; else if (c1[d] < c2[d]) ++c1[d];
				if (c2[d] > c1[d]) --c2[d]; else if (c2[d] < c1[d]) ++c2[d]; else ++equal;
			}
			if (equal < 3) {
				cscale *= get_flow_val(c2, c1, 1);
				if (cscale < CTHRESH) continue;
			}
		}
		colorRGBA const &lsc(lsrc.get_color());
		UNROLL_3X(ls[i_] += lsc[i_]*cscale;)
		added = 1;
	}
	return added;
}


// used on mesh and water
void get_sd_light(int x, int y, int z, point const &p, bool no_dynamic, float lightscale, float *ls, vector3d const *const norm, float const *const spec) {

	assert(lm_alloc && lmap_manager.vlmap);

	if (using_lightmap && !light_sources.empty() && lmap_manager.is_valid_cell(x, y, z)) {
		float const *const color(lmap_manager.vlmap[y][x][z].c);
		ADD_LIGHT_CONTRIB(color, ls);
	}
	if (!no_dynamic && !dl_sources.empty()) get_dynamic_light(x, y, z, p, lightscale, ls, norm, spec);
}


float get_indir_light(colorRGBA &a, colorRGBA cscale, point const &p, bool no_dynamic, bool shadowed, vector3d const *const norm, float const *const spec) {

	//UNROLL_3X(a[i_] *= cscale[i_];) return 1.0;
	assert(lm_alloc && lmap_manager.vlmap);
	bool const global_lighting(read_light_files[LIGHTING_SKY] || write_light_files[LIGHTING_SKY]);
	float val(MAX_LIGHT);
	bool outside_mesh(0);
	colorRGBA ls(BLACK);
	point const p_adj((norm && !global_lighting) ? (p + (*norm)*(0.25*HALF_DXY)) : p);
	int const x(get_xpos(p_adj.x - SHIFT_DX)), y(get_ypos(p_adj.y - SHIFT_DY)), z(get_zpos(p_adj.z));
	
	if (point_outside_mesh(x, y)) {
		outside_mesh = 1; // outside the range
	}
	else if (p.z <= czmin0) {
		outside_mesh = is_under_mesh(p);
		if (outside_mesh) val = 0.0; // under all collision objects (is this correct?)
	}
	else if (using_lightmap && p.z < czmax && lmap_manager.vlmap[y][x] != NULL) { // not above all collision objects and not empty cell
		lmcell const &lmc(lmap_manager.vlmap[y][x][z]);
		
		if (shadowed) { // Note: this test is optional
			val = lmc.v; // Note: could interpolate between voxels here, but it's slow and doesn't look any better
			if (val > 0.0 && global_lighting) {UNROLL_3X(cscale[i_] *= lmc.ac[i_];)} // add indirect color
		}
		ADD_LIGHT_CONTRIB(lmc.c, ls);
	}
	if (!no_dynamic && !outside_mesh && !dl_sources.empty() && p.z < dlight_bb[2][1] && p.z > dlight_bb[2][0]) {
		get_dynamic_light(x, y, z, p, 1.0, (float *)&ls, norm, spec);
	}
	UNROLL_3X(a[i_] *= (cscale[i_]*val + ls[i_]);) // unroll the loop
	return val;
}


// within a sphere, unless radius == 0.0
unsigned enable_dynamic_lights(point const center, float radius) {

	point const camera(get_camera_pos());
	vector<pair<float, unsigned> > vis_lights;

	for (unsigned i = 0; i < dl_sources.size(); ++i) {
		light_source const &ls(dl_sources[i]);
		float const ls_radius(ls.get_radius());
		point const &ls_center(ls.get_center());
		if (ls_radius == 0.0) continue; // not handling zero radius lights yet
		if (radius > 0.0 && !dist_less_than(center, ls_center, (radius + ls_radius))) continue;
		if (!sphere_in_camera_view(ls_center, ls_radius, 0)) continue;
		float weight(p2p_dist(ls_center, camera)); // distance from light to camera
		if (radius > 0.0) weight += p2p_dist(ls_center, center); // distance from light to object center
		vis_lights.push_back(make_pair(weight/ls_radius, i));
	}
	sort(vis_lights.begin(), vis_lights.end());
	unsigned const num_dlights(min(vis_lights.size(), MAX_LIGHTS));

	for (unsigned i = 0; i < num_dlights; ++i) {
		int const gl_light(START_LIGHT+i);
		light_source const &ls(dl_sources[vis_lights[i].second]);
		float udiffuse[4] = {0}; // diffuse = 0 because we don't have correct normals
		set_colors_and_enable_light(gl_light, (float *)(&ls.get_color()), udiffuse);
		glLightf(gl_light, GL_CONSTANT_ATTENUATION,  1.0);
		glLightf(gl_light, GL_LINEAR_ATTENUATION,    0.0);
		glLightf(gl_light, GL_QUADRATIC_ATTENUATION, 10.0/(ls.get_radius()*ls.get_radius()));
		set_gl_light_pos(gl_light, ls.get_center(), 1.0); // point light source position
	}
	return num_dlights;
}


void disable_dynamic_lights(unsigned num_dlights) {

	for (int i = START_LIGHT; i < int(START_LIGHT+num_dlights); ++i) {
		glDisable(i);
	}
}


