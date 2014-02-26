// 3D World - lighting code, incuding static and dynamic lights, profile generation, and flow calculation
// by Frank Gennari
// 1/16/06
#include "3DWorld.h"
#include "mesh.h"
#include "csg.h"
#include "lightmap.h"
#include "gl_ext_arb.h"
#include "shaders.h"


bool const SHOW_STAT_LIGHTS  = 0; // debugging
bool const SHOW_DYNA_LIGHTS  = 0; // debugging
unsigned const NUM_RAND_LTS  = 0;

int      const START_LIGHT   = 2;
int      const END_LIGHT     = 8; // one past the end
unsigned const MAX_LIGHTS    = unsigned(END_LIGHT - START_LIGHT);

float const CTHRESH          = 0.025;
float const SQRT_CTHRESH     = sqrt(CTHRESH);
float const DZ_VAL_SCALE     = 2.0;
float const SHIFT_VAL        = 0.5; // hack to fix some offset problem
float const LT_DIR_FALLOFF   = 0.005;
float const LT_DIR_FALLOFF_INV(1.0/LT_DIR_FALLOFF);
float const DARKNESS_THRESH  = 0.1;
float const DEF_SKY_GLOBAL_LT= 0.25; // when ray tracing is not used


bool using_lightmap(0), lm_alloc(0), has_dl_sources(0), has_spotlights(0), has_line_lights(0), use_dense_voxels(0), has_indir_lighting(0);
unsigned dl_tid(0), elem_tid(0), gb_tid(0);
float DZ_VAL2(DZ_VAL/DZ_VAL_SCALE), DZ_VAL_INV2(1.0/DZ_VAL2), SHIFT_DX(SHIFT_VAL*DX_VAL), SHIFT_DY(SHIFT_VAL*DY_VAL);
float czmin0(0.0), lm_dz_adj(0.0), dlight_add_thresh(0.0);
float dlight_bb[3][2] = {0}, SHIFT_DXYZ[3] = {SHIFT_DX, SHIFT_DY, 0.0};
dls_cell **ldynamic = NULL;
vector<light_source> light_sources, dl_sources, dl_sources2; // static, dynamic {cur frame, next frame}
lmap_manager_t lmap_manager;


extern int animate2, display_mode, frame_counter, camera_coll_id, read_light_files[], write_light_files[];
extern unsigned create_voxel_landscape;
extern float czmin, czmax, fticks, zbottom, ztop, XY_SCENE_SIZE, indir_light_exp, light_int_scale[];
extern colorRGBA cur_ambient, cur_diffuse;
extern coll_obj_group coll_objects;
extern vector<light_source> enabled_lights;


// *** USEFUL INLINES ***

inline bool add_cobj_ok(coll_obj const &cobj) { // skip small things like tree leaves and such
	
	return (cobj.fixed && !cobj.disabled() && cobj.volume > 0.0001); // cobj.type == COLL_CUBE
}


// *** LIGHT_SOURCE IMPLEMENTATION ***


// radius == 0.0 is really radius == infinity (no attenuation)
light_source::light_source(float sz, point const &p, point const &p2, colorRGBA const &c, bool id, vector3d const &d, float bw, float ri) :
	dynamic(id), radius(sz), radius_inv((radius == 0.0) ? 0.0 : 1.0/radius),
	r_inner(ri), bwidth(bw), pos(p), pos2(p2), dir(d.get_norm()), color(c)
{
	assert(bw > 0.0 && bw <= 1.0);
	assert(!(is_directional() && is_line_light())); // can't be both
}


void light_source::add_color(colorRGBA const &c) {

	color = color*color.alpha + c*c.alpha;
	color.alpha = 1.0;
}


float light_source::get_intensity_at(point const &p, point &updated_lpos) const {

	if (radius == 0.0) return color[3]; // no falloff

	if (is_line_light()) {
		vector3d const L(pos2 - pos);
		updated_lpos = pos + L*CLIP_TO_01(dot_product((p - pos), L)/L.mag_sq());
	}
	else {
		updated_lpos = pos;
	}
	if (fabs(p.z - updated_lpos.z) > radius) return 0.0; // fast test
	float const dist_sq(p2p_dist_sq(updated_lpos, p));
	if (dist_sq > radius*radius) return 0.0;
	float const rscale((radius - sqrt(dist_sq))*radius_inv);
	return rscale*rscale*color[3]; // quadratic 1/r^2 attenuation
}


float light_source::get_dir_intensity(vector3d const &obj_dir) const {

	if (bwidth == 1.0) return 1.0;
	float const dp(dot_product(obj_dir, dir));
	if (dp >= 0.0 && (bwidth + LT_DIR_FALLOFF) < 0.5) return 0.0;
	float const dp_norm(0.5*(-dp*InvSqrt(obj_dir.mag_sq()) + 1.0)); // dp = -1.0 to 1.0, bw = 0.0 to 1.0
	return CLIP_TO_01(2.0f*(dp_norm + bwidth + LT_DIR_FALLOFF - 1.0f)*LT_DIR_FALLOFF_INV);
}


void light_source::get_bounds(point bounds[2], int bnds[3][2], float sqrt_thresh, vector3d const &bounds_offset) const {

	if (radius == 0.0) { // global light source
		for (unsigned d = 0; d < 3; ++d) {
			bounds[0][d] = -SCENE_SIZE[d];
			bounds[1][d] =  SCENE_SIZE[d];
			bnds[d][0]   = 0;
			bnds[d][1]   = MESH_SIZE[d]-1;
		}
	}
	else {
		float const rb(radius*(1.0 - sqrt_thresh));

		for (unsigned d = 0; d < 3; ++d) {
			bounds[0][d] = min(pos[d], pos2[d]) - rb; // lower
			bounds[1][d] = max(pos[d], pos2[d]) + rb; // upper

			for (unsigned j = 0; j < 2; ++j) {
				bnds[d][j] = max(0, min(MESH_SIZE[d]-1, get_dim_pos((bounds[j][d] + bounds_offset[d]), d)));
			}
		}
	}
}


bool light_source::is_visible() const {

	if (radius == 0.0 || sphere_in_camera_view(pos, radius, 0)) return 1; // max_level?
	if (is_line_light() && sphere_in_camera_view(0.5*(pos + pos2), (radius + 0.5*p2p_dist(pos, pos2)), 0)) return 1; // use bounding sphere
	return 0;
}


void light_source::combine_with(light_source const &l) {

	assert(radius > 0.0);
	float const w1(radius*radius*radius), w2(l.radius*l.radius*l.radius), wsum(w1 + w2), wa(w1/wsum), wb(w2/wsum);
	radius     = pow(wsum, (1.0f/3.0f));
	radius_inv = 1.0/radius;
	pos       *= wa;
	pos       += l.pos*wb; // weighted average
	blend_color(color, color, l.color, wa, 1);
}


void light_source::draw(int ndiv) const {

	if (radius == 0.0) return;
	set_color(color);
	draw_sphere_vbo(pos, 0.05*radius, ndiv, 0);
	if (pos2 != pos) {draw_sphere_vbo(pos2, 0.05*radius, ndiv, 0);} // line light source, draw both points
}


void shift_light_sources(vector3d const &vd) {

	for (unsigned i = 0; i < light_sources.size(); ++i) {
		light_sources[i].shift_by(vd);
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
		else {avg_alpha = 1.0;}
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
	unsigned const nrects((unsigned)rects.size());

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
	unsigned const nrects((unsigned)rects.size());
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

	for (unsigned i = 0; i < (unsigned)coll_objects.size(); ++i) {
		coll_objects[i].counter = -1;
	}
}


void lmcell::get_final_color(colorRGB &color, float max_indir, float indir_scale, float extra_ambient) const {

	float const max_s(max(sc[0], max(sc[1], sc[2])));
	float const max_g(max(gc[0], max(gc[1], gc[2])));
	float const sv_scaled((max_s > 0.0 && sv > 0.0) ? min(1.0f, sv*light_int_scale[LIGHTING_SKY   ])/max_s : 0.0);
	float const gv_scaled((max_g > 0.0 && gv > 0.0) ? min(1.0f, gv*light_int_scale[LIGHTING_GLOBAL])/max_g : 0.0);

	UNROLL_3X(float indir_term((sv_scaled*sc[i_] + extra_ambient)*cur_ambient[i_] + gv_scaled*gc[i_]*cur_diffuse[i_]); \
			  if (indir_term > 0.0 && indir_light_exp != 1.0) indir_term = pow(indir_term, indir_light_exp); \
			  color[i_] = min(max_indir, indir_scale*indir_term) + min(1.0f, lc[i_]*light_int_scale[LIGHTING_LOCAL]);)
}


void lmcell::set_outside_colors() {

	sv = 1.0;
	gv = 0.0;
	UNROLL_3X(sc[i_] = gc[i_] = 1.0; lc[i_] = 0.0;)
}


bool lmap_manager_t::is_valid_cell(int x, int y, int z) const {
	return (z >= 0 && z < MESH_SIZE[2] && !point_outside_mesh(x, y) && vlmap[y][x] != NULL);
}


lmcell *lmap_manager_t::get_lmcell(point const &p) {

	int const x(get_xpos(p.x - SHIFT_DX)), y(get_ypos(p.y - SHIFT_DY)), z(get_zpos(p.z));
	return (is_valid_cell(x, y, z) ? &vlmap[y][x][z] : NULL);
}


template<typename T> void lmap_manager_t::alloc(unsigned nbins, unsigned zsize, T **nonempty_bins, lmcell const &init_lmcell) {

	assert(nonempty_bins != NULL);
	if (vlmap == NULL) {matrix_gen_2d(vlmap);} // create column headers once
	lm_zsize = zsize;
	vldata_alloc.resize(max(nbins, 1U), init_lmcell); // make size at least 1, even if there are no bins, so we can test on emptiness
	unsigned cur_v(0);

	// initialize lightmap
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (!nonempty_bins[i][j]) {
				vlmap[i][j] = NULL;
				continue;
			}
			assert(cur_v + lm_zsize <= vldata_alloc.size());
			vlmap[i][j] = &vldata_alloc[cur_v];
			cur_v      += lm_zsize;
		}
	}
	assert(cur_v == nbins);
}


void lmap_manager_t::init_from(lmap_manager_t const &src) {

	//assert(!is_allocated());
	//clear_cells(); // probably unnecessary
	alloc(src.vldata_alloc.size(), src.lm_zsize, src.vlmap, lmcell());
	copy_data(src);
}


// *this = blend_weight*dest + (1.0 - blend_weight)*(*this)
void lmap_manager_t::copy_data(lmap_manager_t const &src, float blend_weight) {

	assert(vlmap && src.vlmap);
	assert(src.lm_zsize == lm_zsize);
	assert(src.vldata_alloc.size() == vldata_alloc.size());
	assert(blend_weight >= 0.0);
	if (blend_weight == 0.0) return; // keep existing dest

	if (blend_weight == 1.0) {
		vldata_alloc = src.vldata_alloc; // deep copy all lmcell data
		return;
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) { // openmp?
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (!vlmap[i][j]) {assert(!src.vlmap[i][j]); continue;}
			assert(src.vlmap[i][j]);
			
			for (unsigned z = 0; z < lm_zsize; ++z) {
				vlmap[i][j][z].mix_lighting_with(src.vlmap[i][j][z], blend_weight);
			}
		}
	}
}


// *this = val*lmc + (1.0 - val)*(*this)
void lmcell::mix_lighting_with(lmcell const &lmc, float val) {

	float const omv(1.0 - val); // Note: we ignore the flow values and smoke for now
	sv = val*lmc.sv + omv*sv;
	gv = val*lmc.gv + omv*gv;
	UNROLL_3X(sc[i_] = val*lmc.sc[i_] + omv*sc[i_];)
	UNROLL_3X(gc[i_] = val*lmc.gc[i_] + omv*gc[i_];)
	UNROLL_3X(lc[i_] = val*lmc.lc[i_] + omv*lc[i_];)
}


bool has_fixed_cobjs(int x, int y) {

	assert(!point_outside_mesh(x, y));
	vector<int> const &cvals(v_collision_matrix[y][x].cvals);

	for (vector<int>::const_iterator i = cvals.begin(); i != cvals.end(); ++i) {
		if (coll_objects[*i].fixed && coll_objects[*i].status == COLL_STATIC) {return 1;}
	}
	return 0;
}

void regen_lightmap() {

	if (MESH_Z_SIZE == 0) return; // not using lmap
	assert(lmap_manager.is_allocated());
	clear_lightmap();
	assert(!lmap_manager.is_allocated());
	build_lightmap(0);
	assert(lmap_manager.is_allocated());
}


void clear_lightmap() {

	if (!lmap_manager.is_allocated()) return;
	kill_current_raytrace_threads(); // kill raytrace threads and wait for them to finish since they are using the current lightmap
	lmap_manager.clear_cells();
	using_lightmap = 0;
	lm_alloc       = 0;
	czmin0         = czmin;
}


void calc_flow_for_xy(r_profile flow_prof[3], int i, int j, bool proc_cobjs, float zstep) {

	assert(zstep > 0.0);
	lmcell *vldata(lmap_manager.get_column(j, i));
	if (vldata == NULL) return;
	float const bbz[2][2] = {{get_xval(j), get_xval(j+1)}, {get_yval(i), get_yval(i+1)}}; // X x Y
	coll_cell const &cell(v_collision_matrix[i][j]);
	unsigned const ncv((unsigned)cell.cvals.size());
	vector<pair<float, unsigned> > cobj_z;

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
			}
		}
		sort(cobj_z.begin(), cobj_z.end(), std::greater<pair<float, unsigned> >()); // max to min z
	}
	unsigned const ncv2((unsigned)cobj_z.size());

	for (int v = MESH_SIZE[2]-1; v >= 0; --v) { // top to bottom
		float zb(czmin0 + v*zstep), zt(zb + zstep); // cell Z bounds
				
		if (zt < mesh_height[i][j]) { // under mesh
			UNROLL_3X(vldata[v].pflow[i_] = 0;) // all zeros
		}
		else if (!proc_cobjs) {
			UNROLL_3X(vldata[v].pflow[i_] = 255;) // all ones
		}
		else { // above mesh case
			float const bb[3][2]  = {{bbz[0][0], bbz[0][1]}, {bbz[1][0], bbz[1][1]}, {zb, zt}};
			float const bbx[2][2] = {{bb[1][0], bb[1][1]}, {zb, zt}}; // YxZ
			float const bby[2][2] = {{zb, zt}, {bb[0][0], bb[0][1]}}; // ZxX
			flow_prof[0].reset_bbox(bbx);
			flow_prof[1].reset_bbox(bby);
			flow_prof[2].reset_bbox(bbz);
			
			for (unsigned c2 = 0; c2 < ncv2; ++c2) { // could make this more efficient
				coll_obj &cobj(coll_objects[cobj_z[c2].second]);
				if (cobj.d[0][0] >= bb[0][1] || cobj.d[0][1]     <= bb[0][0]) continue; // no intersection
				if (cobj.d[1][0] >= bb[1][1] || cobj.d[1][1]     <= bb[1][0]) continue;
				if (cobj.d[2][0] >= bb[2][1] || cobj_z[c2].first <= bb[2][0]) continue;
				float const cztop(cobj.d[2][1]);
				cobj.d[2][1] = cobj_z[c2].first;
						
				for (unsigned d = 0; d < 3; ++d) { // critical path
					flow_prof[d].add_rect(cobj.d, (d+1)%3, (d+2)%3, 1.0);
				}
				cobj.d[2][1] = cztop; // restore original value
			} // for c2
			for (unsigned e = 0; e < 3; ++e) {
				float const fv(flow_prof[e].den_inv());
				assert(fv > -TOLER);
				vldata[v].pflow[e] = (unsigned char)(255.5*CLIP_TO_01(fv));
			}
		} // if above mesh
	} // for v
}


float calc_czspan() {return max(0.0f, ((czmax + lm_dz_adj) - czmin0 + TOLER));}


void build_lightmap(bool verbose) {

	if (lm_alloc) return; // what about recreating the lightmap if the scene has changed?

	// prevent the z range from being empty/denormalized when there are no cobjs
	if (use_dense_voxels) {
		czmin = min(czmin, zbottom);
		czmax = max(czmax, (czmin + Z_SCENE_SIZE - 0.5f*DZ_VAL));
	}
	else if (czmin >= czmax) {
		czmin = min(czmin, zbottom);
		czmax = max(czmax, ztop);
	}

	// calculate and allocate some data we need even if the lmap is not used
	assert(DZ_VAL > 0.0);
	DZ_VAL2       = DZ_VAL/DZ_VAL_SCALE;
	DZ_VAL_INV2   = 1.0/DZ_VAL2;
	SHIFT_DXYZ[0] = SHIFT_DX = SHIFT_VAL*DX_VAL;
	SHIFT_DXYZ[1] = SHIFT_DY = SHIFT_VAL*DY_VAL;
	czmin0        = czmin;//max(czmin, zbottom);
	assert(lm_dz_adj >= 0.0);
	if (!ldynamic) matrix_gen_2d(ldynamic, MESH_X_SIZE, MESH_Y_SIZE);
	if (MESH_Z_SIZE == 0) return;

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
			has_fixed        |= fixed; // only used in an assertion below
			if (need_lmcell[i][j]) ++nonempty;
		}
	}

	// add cells surrounding static scene lights
	// Note: this isn't really necessary when using ray casting for lighting,
	//       but it helps ensure there are lmap cells around light sources to light the dynamic objects
	for (unsigned i = 0; i < light_sources.size(); ++i) {
		point bounds[2]; // unused
		int bnds[3][2];
		light_sources[i].get_bounds(bounds, bnds, SQRT_CTHRESH);

		for (int y = bnds[1][0]; y <= bnds[1][1]; ++y) {
			for (int x = bnds[0][0]; x <= bnds[0][1]; ++x) {
				if (!need_lmcell[y][x]) ++nonempty;
				need_lmcell[y][x] |= 2;
			}
		}
	}

	// determine allocation and voxel grid sizes
	reset_cobj_counters();
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
	if (verbose) {cout << "Lightmap zsize= " << zsize << ", nonempty= " << nonempty << ", bins= " << nbins << ", czmin= " << czmin0 << ", czmax= " << czmax << endl;}
	assert(zstep > 0.0);
	bool raytrace_lights[3];
	UNROLL_3X(raytrace_lights[i_] = (read_light_files[i_] || write_light_files[i_]););
	has_indir_lighting = (raytrace_lights[LIGHTING_SKY] || raytrace_lights[LIGHTING_GLOBAL] || create_voxel_landscape);
	lmcell init_lmcell;

	if (!has_indir_lighting) { // set a default value that isn't all black
		init_lmcell.sv = init_lmcell.gv = DEF_SKY_GLOBAL_LT;
		UNROLL_3X(init_lmcell.sc[i_] = init_lmcell.gc[i_] = 1.0;)
	}
	lmap_manager.alloc(nbins, zsize, need_lmcell, init_lmcell);
	assert(ldynamic && lmap_manager.is_allocated());
	using_lightmap = (nonempty > 0);
	lm_alloc       = 1;

	// process vertical (Z) light projections
	r_profile flow_prof[3]; // particle {x, y, z}

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			bool const proc_cobjs(need_lmcell[i][j] & 1);
			calc_flow_for_xy(flow_prof, i, j, proc_cobjs, zstep);
		}
	}
	int const bnds[2][2] = {{0, MESH_X_SIZE-1}, {0, MESH_Y_SIZE-1}};

	// add in static light sources
	if (!raytrace_lights[LIGHTING_LOCAL]) {
		for (unsigned i = 0; i < light_sources.size(); ++i) {
			light_source &ls(light_sources[i]);
			assert(!ls.is_line_light()); // not supported here
			point lpos(ls.get_pos()); // may be updated for line lights (if they're ever supported)
			if (!is_over_mesh(lpos)) continue;
			colorRGBA const &lcolor(ls.get_color());
			point bounds[2];
			int bnds[3][2], cobj(-1), last_cobj(-1);
			CELL_LOC_T cent[3];
			
			for (unsigned i = 0; i < 3; ++i) {
				cent[i] = max(0, min(MESH_SIZE[i]-1, get_dim_pos(lpos[i], i))); // clamp to mesh bounds
			}
			ls.get_bounds(bounds, bnds, SQRT_CTHRESH);
			check_coll_line(lpos, lpos, cobj, -1, 1, 2, 1); // check cobj containment and ignore that shape (ignore voxels)

			for (int y = bnds[1][0]; y <= bnds[1][1]; ++y) {
				for (int x = bnds[0][0]; x <= bnds[0][1]; ++x) {
					assert(lmap_manager.get_column(x, y));
					float const xv(get_xval(x)), yv(get_yval(y));

					for (int z = bnds[2][0]; z <= bnds[2][1]; ++z) {
						assert(unsigned(z) < zsize);
						point const p(xv, yv, get_zval(z));
						float cscale(ls.get_intensity_at(p, lpos));
						if (cscale < CTHRESH) {if (z > cent[2]) break; else continue;}
				
						if (ls.is_directional()) {
							cscale *= ls.get_dir_intensity(lpos - p);
							if (cscale < CTHRESH) continue;
						}
						point const lpos_ext(lpos + HALF_DXY*(p - lpos).get_norm()); // extend away from light to account for light fixtures
						if ((last_cobj >= 0 && coll_objects[last_cobj].line_intersect(lpos_ext, p)) ||
							check_coll_line(p, lpos_ext, last_cobj, cobj, 1, 3)) {continue;}
						lmcell &lmc(lmap_manager.get_lmcell(x, y, z));
						UNROLL_3X(lmc.lc[i_] = min(1.0f, (lmc.lc[i_] + cscale*lcolor[i_]));) // what about diffuse/normals?
					} // for z
				} // for x
			} // for y
		} // for i
	}
	if (verbose) PRINT_TIME(" Lighting Setup + XYZ Passes");

	if (nbins > 0) {
		// Note: sky and global lighting use the same data structure for reading/writing, so they should have the same filename if used together
		string const type_names[NUM_LIGHTING_TYPES] = {" Sky", " Global", " Local"};

		for (unsigned ltype = 0; ltype < NUM_LIGHTING_TYPES; ++ltype) {
			if (raytrace_lights[ltype]) {
				compute_ray_trace_lighting(ltype);
				if (verbose) {PRINT_TIME((type_names[ltype] + " Lighting Load/Ray Trace").c_str());}
			}
		}
	}
	reset_cobj_counters();
	matrix_delete_2d(need_lmcell);
	PRINT_TIME(" Lighting Total");
}


void update_flow_for_voxels(cube_t const &cube) {

	if (!lm_alloc || !lmap_manager.is_allocated()) return;
	int const cx1(max(0, get_xpos(cube.d[0][0]))), cx2(min(MESH_X_SIZE-1, get_xpos(cube.d[0][1])));
	int const cy1(max(0, get_ypos(cube.d[1][0]))), cy2(min(MESH_Y_SIZE-1, get_ypos(cube.d[1][1])));
	float const zstep(calc_czspan()/MESH_SIZE[2]);
	r_profile flow_prof[3];

	for (int y = cy1; y <= cy2; ++y) {
		for (int x = cx1; x <= cx2; ++x) {
			assert(!point_outside_mesh(x, y));
			bool const fixed(!coll_objects.empty() && has_fixed_cobjs(x, y));
			bool const proc_cobjs(use_dense_voxels || fixed);
			calc_flow_for_xy(flow_prof, y, x, proc_cobjs, zstep);
		} // for x
	} //for y
}


// *** Dynamic Lights Code ***


void setup_2d_texture(unsigned &tid) {

	setup_texture(tid, 0, 0, 0, 0, 0, 1);
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
void upload_dlights_textures(cube_t const &bounds) {

	//RESET_TIME;
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
	unsigned const ndl(min(max_dlights, (unsigned)dl_sources.size()));
	unsigned const ysz(floats_per_light/4);
	float const radius_scale(1.0/(0.5*(bounds.d[0][1] - bounds.d[0][0]))); // bounds x radius inverted
	vector3d const poff(bounds.get_llc()), psize(bounds.get_urc() - poff);
	vector3d const pscale(1.0/psize.x, 1.0/psize.y, 1.0/psize.z);
	has_spotlights = has_line_lights = 0;

	for (unsigned i = 0; i < ndl; ++i) {
		bool const line_light(dl_sources[i].is_line_light());
		float *data(dl_data + i*floats_per_light);
		dl_sources[i].pack_to_floatv(data); // {center,radius, color, dir,beamwidth}
		UNROLL_3X(data[i_] = (data[i_] - poff[i_])*pscale[i_];) // scale to [0,1] range
		UNROLL_3X(data[i_+4] *= 0.1;) // scale color down
		if (line_light) {UNROLL_3X(data[i_+8] = (data[i_+8] - poff[i_])*pscale[i_];)} // scale to [0,1] range
		data[3] *= radius_scale;
		has_spotlights  |= dl_sources[i].is_directional();
		has_line_lights |= line_light;
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
	static vector<unsigned> gb_data;
	gb_data.resize(XY_MULT_SIZE, 0);
	unsigned const elem_tex_sz = 256; // must agree with value in shader
	unsigned const max_gb_entries(elem_tex_sz*elem_tex_sz);
	unsigned short elem_data[max_gb_entries] = {0};
	unsigned elix(0);

	for (int y = 0; y < MESH_Y_SIZE && elix < max_gb_entries; ++y) {
		for (int x = 0; x < MESH_X_SIZE && elix < max_gb_entries; ++x) {
			unsigned const gb_ix(x + y*MESH_X_SIZE); // {start, end, unused}
			gb_data[gb_ix] = elix; // start_ix
			vector<unsigned short> const &ixs(ldynamic[y][x].get_src_ixs());
			unsigned const num_ixs(min((unsigned)ixs.size(), 256U)); // max of 256 lights per bin
			
			for (unsigned i = 0; i < num_ixs && elix < max_gb_entries; ++i) { // end if exceed max entries
				if (ixs[i] >= ndl) continue; // dlight index is too high, skip
				elem_data[elix++] = (unsigned short)ixs[i];
			}
			gb_data[gb_ix] += (elix << 16); // end_ix
		}
	}
	if (elix > 0.9*max_gb_entries) {
		if (elix >= max_gb_entries) {std::cerr << "Warning: Exceeded max number of indexes in dynamic light texture upload" << endl;}
		//cout << "elix: " << elix << ", dlight_add_thresh: " << dlight_add_thresh << endl;
		dlight_add_thresh = min(0.25, (dlight_add_thresh + 0.005)); // increase thresh to clip the dynamic lights to a smaller radius
	}
	if (elem_tid == 0) {
		setup_2d_texture(elem_tid);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R16UI, elem_tex_sz, elem_tex_sz, 0, GL_RED_INTEGER, GL_UNSIGNED_SHORT, elem_data);
	}
	else {
		bind_2d_texture(elem_tid);
		unsigned const height(min(elem_tex_sz, (elix/elem_tex_sz+1))); // approximate ceiling
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, elem_tex_sz, height, GL_RED_INTEGER, GL_UNSIGNED_SHORT, elem_data);
	}

	// step 3: grid bag(s)
	if (gb_tid == 0) {
		setup_2d_texture(gb_tid);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, MESH_X_SIZE, MESH_Y_SIZE, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &gb_data.front()); // Nx x Ny
	}
	else {
		bind_2d_texture(gb_tid);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, MESH_X_SIZE, MESH_Y_SIZE, GL_RED_INTEGER, GL_UNSIGNED_INT, &gb_data.front());
	}
	//PRINT_TIME("Dlight Texture Upload");
	//cout << "ndl: " << ndl << ", elix: " << elix << ", gb_sz: " << XY_MULT_SIZE << endl;
}


void set_one_texture(shader_t &s, unsigned tid, unsigned tu_id, const char *const name) {

	set_active_texture(tu_id); // texture unit
	bind_2d_texture(tid);
	s.add_uniform_int(name, tu_id);
}


void setup_dlight_textures(shader_t &s) {

	assert(dl_tid > 0 && elem_tid > 0 && gb_tid > 0 );
	set_one_texture(s, dl_tid,   2, "dlight_tex");
	set_one_texture(s, elem_tid, 3, "dlelm_tex");
	set_one_texture(s, gb_tid,   4, "dlgb_tex");
	set_active_texture(0);
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

	float const bwidth = 0.02;
	float const radius = 4.0;
	point const camera(get_camera_pos());
	add_dynamic_light(radius, camera, SUN_C, cview_dir, bwidth);

	if (display_mode & 0x10) { // add one bounce of indirect lighting
		unsigned const NUM_VPLS = 32;
		float const theta(acosf(1.0f - bwidth /*- 0.5*LT_DIR_FALLOFF*/)); // flashlight beam angle
		float const rad_per_len(tan(theta));
		vector3d vab[2];
		get_ortho_vectors(cview_dir, vab);

		for (unsigned i = 0; i < NUM_VPLS; ++i) {
			float const a(TWO_PI*i/NUM_VPLS);
			vector3d const delta((sin(a)*vab[0] + cos(a)*vab[1]).get_norm()); // already normalized?
			vector3d const dir(cview_dir + rad_per_len*delta);
			int cindex;
			point cpos;
			vector3d cnorm;
			
			if (check_coll_line_exact(camera, (camera + 0.5*radius*dir), cpos, cnorm, cindex, 0.0, camera_coll_id, 1, 0, 0)) {
				cpos -= 0.0001*radius*cnorm; // move behind the collision plane so as not to multiply light
				assert(cindex >= 0);
				colorRGBA const color(SUN_C.modulate_with(coll_objects[cindex].get_avg_color()));
				add_dynamic_light(0.1*radius, cpos, color*0.15, cnorm, 0.5); // wide angle (almost hemisphere)
			}
		}
	}
}


void add_dynamic_light(float sz, point const &p, colorRGBA const &c, vector3d const &d, float bw, point *line_end_pos) {

	if (!animate2) return;
	float const sz_scale((world_mode == WMODE_UNIVERSE) ? 1.0 : sqrt(0.1*XY_SCENE_SIZE));
	dl_sources2.push_back(light_source(sz_scale*sz, p, (line_end_pos ? *line_end_pos : p), c, 1, d, bw));
}


void add_line_light(point const &p1, point const &p2, colorRGBA const &color, float size, float intensity) {

	if (!animate2) return;
	point p[2] = {p1, p2};
	if (!do_line_clip_scene(p[0], p[1], zbottom, max(ztop, czmax))) return;
	float const radius(size*intensity), pt_offset((1.0 - SQRTOFTWOINV)*radius);

	if (dist_less_than(p1, p2, radius)) { // short segment, use a single point light
		add_dynamic_light(radius, p[0], color);
	}
	else { // add a real line light
		vector3d const dir((p[1] - p[0]).get_norm());
		p[0] += dir*pt_offset; // shrink line slightly for a better effect
		p[1] -= dir*pt_offset;
		add_dynamic_light(radius, p[0], color, dir, 1.0, &p[1]);
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
	if (!dist_less_than(pos, ls.pos, 0.2*max(HALF_DXY, radius))) return 0;
	if (ls.bwidth != bwidth || ls.r_inner != r_inner || ls.dynamic != dynamic) return 0;
	if (bwidth < 1.0 && dot_product(dir, ls.dir) < 0.95) return 0;
	if (pos != pos2 || ls.pos != ls.pos2) return 0; // don't merge line lights
	colorRGBA lcolor(color);
	float const rr(radius/ls.radius);
	lcolor.alpha *= rr*rr; // scale by radius ratio squared
	ls.add_color(lcolor);
	return 1;
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


void add_dynamic_lights_ground() {

	//RESET_TIME;
	if (!animate2) return;
	assert(ldynamic);
	clear_dynamic_lights();
	dl_sources.swap(dl_sources2);

	for (unsigned i = 0; i < NUM_RAND_LTS; ++i) { // add some random lights (omnidirectional)
		point const pos(gen_rand_scene_pos());
		dl_sources.push_back(light_source(0.94, pos, pos, BLUE, 1));
	}
	// Note: do we want to sort by y/x position to minimize cache misses?
	sort(dl_sources.begin(), dl_sources.end(), std::greater<light_source>()); // sort by largest to smallest radius
	unsigned const ndl((unsigned)dl_sources.size());
	has_dl_sources = (ndl > 0);
	dlight_add_thresh *= 0.99;
	bool first(1);
	float const sqrt_dlight_add_thresh(sqrt(dlight_add_thresh));
	point const dlight_shift(-SHIFT_DX, -SHIFT_DY, 0.0);

	for (unsigned i = 0; i < ndl; ++i) {
		light_source &ls(dl_sources[i]);
		if (!ls.is_visible()) continue; // view culling
		float const ls_radius(ls.get_radius());
		if ((min(ls.get_pos().z, ls.get_pos2().z) - ls_radius) > max(ztop, czmax)) continue; // above everything, rarely occurs
		point const &lpos(ls.get_pos()), &lpos2(ls.get_pos2());
		bool const line_light(ls.is_line_light());
		int const xcent(get_xpos(lpos.x)), ycent(get_ypos(lpos.y));
		if (!line_light && !point_outside_mesh(xcent, ycent) && !ldynamic[ycent][xcent].check_add_light(i)) continue;
		point bounds[2];
		int bnds[3][2];
		unsigned const ix(i);
		ls.get_bounds(bounds, bnds, sqrt_dlight_add_thresh, dlight_shift);
		
		for (unsigned j = 0; j < 3; ++j) {
			dlight_bb[j][0] = (first ? bounds[0][j] : min(dlight_bb[j][0], bounds[0][j]));
			dlight_bb[j][1] = (first ? bounds[1][j] : max(dlight_bb[j][1], bounds[1][j]));
		}
		first = 0;
		int const xsize(bnds[0][1]-bnds[0][0]), ysize(bnds[1][1]-bnds[1][0]);
		int const radius((max(xsize, ysize)>>1)+2), rsq(radius*radius);
		float const line_rsq((ls_radius + HALF_DXY)*(ls_radius + HALF_DXY));

		for (int y = bnds[1][0]; y <= bnds[1][1]; ++y) {
			int const y_sq((y-ycent)*(y-ycent));

			for (int x = bnds[0][0]; x <= bnds[0][1]; ++x) {
				if (rsq > 1) {
					if (line_light) {
						float const px(get_xval(x)), py(get_yval(y)), lx(lpos2.x - lpos.x), ly(lpos2.y - lpos.y);
						float const cp_mag(lx*(lpos.y - py) - ly*(lpos.x - px));
						if (cp_mag*cp_mag > line_rsq*(lx*lx + ly*ly)) continue;
					}
					else if (((x-xcent)*(x-xcent) + y_sq) > rsq) {continue;} // skip
				}
				ldynamic[y][x].add_light(ix, bounds[0][2], bounds[1][2]); // could do flow clipping here?
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


void light_source::pack_to_floatv(float *data) const {

	// store light_source as: {pos.xyz, radius}, {color.rgba}, {dir.xyz|pos2.xyz, bwidth}
	// Note: we don't really need to store the z-component of dir because we can calculate it from sqrt(1 - x*x - y*y),
	//       but doing this won't save us any texture data so it's not worth the trouble
	assert(data);
	UNROLL_3X(*(data++) = pos[i_];)
	*(data++) = radius;
	UNROLL_3X(*(data++) = 0.5*(1.0 + color[i_]);) // map [-1,1] => [0,1] for negative light support
	*(data++) = color[3];

	if (is_line_light()) { // FIXME: cache this value in light_source for efficiency?
		UNROLL_3X(*(data++) = pos2[i_];)
		*(data++) = 0.0; // pack bwidth as 0 to indicate a line light
	}
	else {
		UNROLL_3X(*(data++) = 0.5*(1.0 + dir  [i_]);) // map [-1,1] to [0,1]
		*(data++) = bwidth; // [0,1]
	}
}


bool is_in_darkness(point const &pos, float radius, int cobj) {

	colorRGBA c(WHITE);
	get_indir_light(c, pos); // this is faster so do it first
	if ((c.R + c.G + c.B) > DARKNESS_THRESH) return 0;

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (is_visible_to_light_cobj(pos, l, radius, cobj, 1)) return 0;
	}
	return 1;
}


void get_dynamic_light(int x, int y, int z, point const &p, float lightscale, float *ls) {

	if (dl_sources.empty()) return;
	assert(!point_outside_mesh(x, y));
	dls_cell const &ldv(ldynamic[y][x]);
	if (!ldv.check_z(p[2])) return;
	unsigned const lsz((unsigned)ldv.size());
	CELL_LOC_T const cl[3] = {x, y, z}; // what about SHIFT_VAL?

	for (unsigned l = 0; l < lsz; ++l) {
		unsigned const ls_ix(ldv.get(l));
		assert(ls_ix < dl_sources.size());
		light_source const &lsrc(dl_sources[ls_ix]);
		point lpos;
		float cscale(lightscale*lsrc.get_intensity_at(p, lpos));
		if (cscale < CTHRESH) continue;
		
		if (lsrc.is_directional()) {
			cscale *= lsrc.get_dir_intensity(lpos - p);
			if (cscale < CTHRESH) continue;
		}
		colorRGBA const &lsc(lsrc.get_color());
		UNROLL_3X(ls[i_] += lsc[i_]*cscale;)
	}
}


// used on mesh and water
void get_sd_light(int x, int y, int z, float *ls) {

	if (lm_alloc && using_lightmap && !light_sources.empty() && lmap_manager.is_valid_cell(x, y, z)) {
		assert(lmap_manager.is_allocated());
		float const *const lcolor(lmap_manager.get_lmcell(x, y, z).lc);
		ADD_LIGHT_CONTRIB(lcolor, ls);
	}
}


float get_indir_light(colorRGBA &a, point const &p) { // Note: return value is unused

	float val(get_voxel_terrain_ao_lighting_val(p)); // shift p?
	if (!lm_alloc) return val;
	assert(lmap_manager.is_allocated());
	bool outside_mesh(0);
	colorRGB cscale(cur_ambient);
	int const x(get_xpos(p.x - SHIFT_DX)), y(get_ypos(p.y - SHIFT_DY)), z(get_zpos(p.z));
	
	if (point_outside_mesh(x, y)) {
		outside_mesh = 1; // outside the range
	}
	else if (p.z <= czmin0) {
		outside_mesh = is_under_mesh(p);
		if (outside_mesh) val = 0.0; // under all collision objects (is this correct?)
	}
	else if (using_lightmap && p.z < czmax && lmap_manager.get_column(x, y) != NULL) { // not above all collision objects and not empty cell
		lmcell const &lmc(lmap_manager.get_lmcell(x, y, z));
		lmc.get_final_color(cscale, 0.5, val);
		val *= lmc.sv + lmc.gv;
	}
	else if (val < 1.0) {
		cscale *= val;
	}
	if (!outside_mesh && !dl_sources.empty() && p.z < dlight_bb[2][1] && p.z > dlight_bb[2][0]) {
		get_dynamic_light(x, y, z, p, 1.0, (float *)&cscale);
	}
	UNROLL_3X(a[i_] *= min(1.0f, cscale[i_]);)
	return val;
}


unsigned enable_dynamic_lights(point const &center, float radius) { // used for tree leaves

	assert(radius > 0.0);
	point const camera(get_camera_pos());
	vector<pair<float, unsigned> > vis_lights;

	for (unsigned i = 0; i < dl_sources.size(); ++i) { // Note: could use ldynamic for faster queries
		light_source const &ls(dl_sources[i]);
		//if (ls.is_directional()) continue; // not correctly handled by GL point lights
		float const ls_radius(ls.get_radius());
		if (ls_radius == 0.0) continue; // not handling zero radius lights yet
		point const &lpos(ls.get_pos());
		float odist_sq(p2p_dist_sq(lpos, center));
		// not correctly handled by GL point lights - we choose the closest end point only
		if (ls.is_line_light()) {odist_sq = min(odist_sq, p2p_dist_sq(ls.get_pos2(), center));}
		if (odist_sq > (radius + ls_radius)*(radius + ls_radius)) continue;
		if (ls.is_directional() && ls.get_dir_intensity(lpos - center) == 0.0) continue; // wrong direction
		if (!ls.is_visible()) continue;
		float cdist_sq(p2p_dist(lpos, camera));
		if (ls.is_line_light()) {cdist_sq = min(cdist_sq, p2p_dist_sq(ls.get_pos2(), camera));} // not entirely correct, but good enough
		float const weight(sqrt(cdist_sq) + sqrt(odist_sq)); // distance from light to camera + distance from light to object center
		vis_lights.push_back(make_pair(weight/ls_radius, i));
	}
	sort(vis_lights.begin(), vis_lights.end());
	unsigned const num_dlights(min((unsigned)vis_lights.size(), MAX_LIGHTS));

	for (unsigned i = 0; i < num_dlights; ++i) {
		int const gl_light(GL_LIGHT0+START_LIGHT+i);
		light_source const &ls(dl_sources[vis_lights[i].second]);
		colorRGBA const &dcolor(ls.get_color()), acolor(dcolor*0.2); // 20% ambient
		set_colors_and_enable_light(gl_light, (float *)(&acolor), (float *)(&dcolor));
		setup_gl_light_atten(gl_light, 1.0, 0.0, 10.0/(ls.get_radius()*ls.get_radius()));
		point lpos(ls.get_pos());
		if (ls.is_line_light() && p2p_dist_sq(ls.get_pos2(), center) < p2p_dist_sq(lpos, center)) {lpos = ls.get_pos2();}
		set_gl_light_pos(gl_light, lpos, 1.0); // point light source position
	}
	return num_dlights;
}


void disable_dynamic_lights(unsigned num_dlights) {

	for (int i = START_LIGHT; i < int(START_LIGHT+num_dlights); ++i) {disable_light(i);}
}


