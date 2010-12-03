// 3D World - Dynamic Surface Subdivision, Lighting, and Shadows
// by Frank Gennari
// 4/10/05
#include "3DWorld.h"
#include "mesh.h"
#include "subdiv.h"
#include "lightmap.h"
#include "dynamic_particle.h"
#include "physics_objects.h"
#include "gl_ext_arb.h"


bool const DO_ROTATE         = 1;
bool const VERTEX_LIGHTING   = 1;
bool const VERBOSE_DYNAMIC   = 0;
bool const TEST_DS_TIME      = 0;
bool const MERGE_STRIPS      = 1; // makes display significantly faster but looks a little worse
bool const BETTER_QUALITY    = 1; // slightly slower, but looks better (except for popping artifacts)
bool const FULL_MAP_LOOKUP   = 1;
bool const USE_DLIST         = 1;
bool const ENABLE_DL_LOD     = 1; // faster but lower visual quality
bool const ALLOW_SPECULAR    = 1;
bool const DO_CIRC_PROJ_CLIP = 1;
bool const RENDER_PART_INT   = 1; // if 1, render only the visible portions of the surface using the slow algorithm, if 0 render all of the surface using the fast algorithm
bool const LOD_QUAD_TRIS     = 0; // slightly faster but somewhat lower quality (better off when dlists are enabled)
bool const LOD_NO_SUBDIV     = 0; // faster but causes color popping effects when switching LOD levels on non subdivided shapes
bool const FAST_SHAPE_DRAW   = 0; // disable lighting and subdivision of shapes
int  const USE_MESH_INT      = 1; // 0 = none, 1 = high res, 2 = low res
int  const VERBOSE           = 0; // 0, 1, 2, 3

unsigned const MAX_CALC_PER_FRAME = 1000;

float const MAX_QUAD_LT_SIZE = 0.005;
float const DIST_CUTOFF      = 2.0;
float const DIV_VAL          = 5.0; // larger = larger polygons (less subdivision)
float const TOLER_           = 1.0E-6;
float const DYNAM_RES_SCALE  = 1.0; // smaller subdivision size for dynamic lights (better quality but slower and has popping artifacts)
float const SHAPE_SPLIT_FACT = 0.0125; // smaller = more top level splitting

float const subdiv_size_inv(1.0/MAX_QUAD_LT_SIZE);


int cube_time(0), invalid_shadows(0);
unsigned ncubes(0), npolys(0), dlists(0), nverts(0), nquads(0), nsurfaces(0), nvlight(0), ntested(0);
unsigned ALL_LT[5] = {0};
unsigned long long max_lighted(0);
float L1_SUBDIV_SIZE(1.0);


extern unsigned cobj_counter;
extern int coll_border, begin_motion, num_groups, camera_coll_id, spectate;
extern int display_mode, camera_mode, camera_view, do_zoom, xoff2, yoff2;
extern float max_proj_rad, subdiv_size_mult, ztop, zbottom, zmax, zmin;
extern float DX_VAL, DY_VAL, XY_SCENE_SIZE, czmin, czmax, SHIFT_DX, SHIFT_DY;
extern double camera_zh;
extern point up_vector;
extern vector<int> weap_cobjs;
extern vector<coll_obj> coll_objects;
extern platform_cont platforms; // only needed for empty test
extern vector<light_source> enabled_lights;
extern obj_type object_types[];
extern obj_group obj_groups[];


class shadow_sphere {

public:
	char lighted, ctype;
	int cid;
	float radius;
	point pos;

	shadow_sphere() : lighted(0), ctype(COLL_CUBE), cid(-1), radius(0.0) {}

	shadow_sphere(point const &pos0, float radius0, int cid0, bool lighted0) :
		lighted(lighted0), cid(cid0), radius(radius0), pos(pos0)
	{
		if (cid < 0) {
			ctype = COLL_SPHERE; // sphere is the default
		}
		else {
			assert(size_t(cid) < coll_objects.size());
			ctype = coll_objects[cid].type;
		}
	}
	inline bool line_intersect(point const &p1, point const &p2) const {
		if (!line_sphere_intersect(p1, p2, pos, radius)) return 0;
		if (ctype == COLL_SPHERE) return 1;
		assert(cid >= 0 && cid < (int)coll_objects.size());
		return coll_objects[cid].line_intersect(p1, p2);
	}
	inline bool test_volume(point const *const pts, unsigned npts, point const &lpos) const {
		if (cid < 0) return 1;
		coll_obj const &c(coll_objects[cid]);
		if (ctype == COLL_SPHERE && (pos != c.points[0] || radius != c.radius)) return 1; // camera sphere != pos
		return !c.cobj_plane_side_test(pts, npts, lpos);
	}
};

vector<shadow_sphere> shadow_objs;


void init_draw_stats() {

	L1_SUBDIV_SIZE = min(SHAPE_SPLIT_FACT*(X_SCENE_SIZE + Y_SCENE_SIZE), subdiv_size_mult*HALF_DXY);

	if (L1_SUBDIV_SIZE > max_proj_rad) {
		cout << "***** Changing max_proj_rad from " << max_proj_rad << " to " << L1_SUBDIV_SIZE << endl;
		set_coll_rmax(L1_SUBDIV_SIZE - min(DX_VAL, DY_VAL)); // may be too late
	}
	cube_time = glutGet(GLUT_ELAPSED_TIME);
	ncubes    = npolys = dlists = nverts = nquads = nsurfaces = nvlight = ntested = 0;
}


void show_draw_stats() {

	if (VERBOSE && (VERBOSE >= 2 || nvlight > 0 || ncubes > 0 || npolys > 0 || nquads > 0)) {
		cout << "Time = " << (glutGet(GLUT_ELAPSED_TIME) - cube_time) << endl;

		if (VERBOSE >= 2 || nvlight > 0) {
			cout << "cubes= " << ncubes << ", polys= " << npolys << ", ndl = " << dlists << ", nverts = " << nverts
				 << ", quads= " << nquads << ", surfs= " << nsurfaces << ", vcalls= " << nvlight << endl
				 << "cobjs = " << coll_objects.size() << ", dynamic spheres = " << shadow_objs.size() << ", tested= " << ntested << endl;
		}
		if (VERBOSE >= 3) {
			// status: COLL_UNUSED COLL_FREED COLL_PENDING  COLL_STATIC COLL_DYNAMIC      COLL_NEGATIVE
			// type  : COLL_NULL   COLL_CUBE  COLL_CYLINDER COLL_SPHERE COLL_CYLINDER_ROT COLL_POLYGON COLL_INVALID
			unsigned tcounts[7] = {0}, scounts[6] = {0}, nsdep(0);

			for (unsigned i = 0; i < coll_objects.size(); ++i) {
				++scounts[unsigned(coll_objects[i].status)];
				if (coll_objects[i].status == COLL_STATIC) ++tcounts[unsigned(coll_objects[i].type)];
				nsdep += coll_objects[i].shadow_depends.size();
			}
			cout << "scounts = ";
			for (unsigned i = 0; i < 6; ++i) {
				cout << scounts[i] << "  ";
			}
			cout << endl << "tcounts = ";
			for (unsigned i = 0; i < 7; ++i) {
				cout << tcounts[i] << "  ";
			}
			cout << endl << "nsdep = " << nsdep << endl;
		}
	}
}


void calc_params(point const *const pts, vector3d *dirs, float *len, float *ninv, unsigned *n, float ss_inv) {

	for (unsigned d = 0; d < 2; ++d) {
		dirs[d]  = pts[1+(d<<1)] - pts[0];
		len[d]   = dirs[d].mag(); // can have zero len when drawing sphere ends
		n[d]     = max(1U, unsigned(len[d]*ss_inv + 0.5));
		ninv[d]  = 1.0/n[d];
		dirs[d] *= ninv[d];
	}
}


struct dqt_params {

	bool no_shadow, double_sided, is_quadric;
	int cobj;
	float ts_x, ts_y, t_tx, t_ty;
	point origin;
	colorRGBA const &color;
	vector<vector<int> > const &stest;
	
	dqt_params(point const &o, bool ns, bool db, int cob, float tx, float ty, float tdx, float tdy, bool iq)
		: no_shadow(ns), double_sided(db), is_quadric(iq), cobj(cob), ts_x(tx), ts_y(ty), t_tx(tdx), t_ty(tdy),
		origin(o), color(coll_objects[cobj].cp.color), stest(coll_objects[cobj].sobjs) {}
};


struct dqd_params {

	bool double_sided;
	int cobj, tri;
	float ts_x, ts_y, t_tx, t_ty;
	point origin;
	vector3d normal;
	colorRGBA const &color;
	quad_div const &qd;
	vector<vector<int> > const &stest;
	unsigned n[2], lod_level, num_internal, num_subdiv, all_lighted;
	float spec[2], len[2], ninv[2];
	vector3d dirs[2];
	unsigned char all_shad, all_unshad;
	bool in_dlist, all_int_surf, no_shadow_calc, use_n;

	dqd_params(dqt_params const &q, int tr, int sub_tri, vector3d const &norm, quad_div const &qd_,
		vector<vector<int> > const &s, bool idl, float const sp[2], unsigned ll, unsigned ns, unsigned al, bool nsc, bool usen) :
	    double_sided(q.double_sided), cobj(q.cobj), tri(tr), ts_x(q.ts_x), ts_y(q.ts_y), t_tx(q.t_tx), t_ty(q.t_ty),
		origin(q.origin), normal(norm), color(q.color), qd(qd_), stest(s), lod_level(ll), num_internal(0), num_subdiv(ns),
		all_lighted(al), all_shad(0xFF), all_unshad(0xFF), in_dlist(idl), all_int_surf(1), no_shadow_calc(nsc), use_n(usen)
	{
		spec[0] = sp[0]; spec[1] = sp[1];
	}

	void do_calc_params(point const *const pts) {
		calc_params(pts, dirs, len, ninv, n, subdiv_size_inv);
	}

	void inherit_flags_from(dqd_params const &p) {
		all_shad     &= p.all_shad;
		all_unshad   &= p.all_unshad;
		all_int_surf &= p.all_int_surf;
	}
};


struct pt_pair { // size = 28

	bool used;
	point p[2];
	pt_pair() : used(0) {}
};


struct vertex_t : public color_wrapper { // size = 32

	point p;
	vector3d n;
};


inline bool light_source::lights_polygon(point const &pc, float rsize, vector3d const* const norm) const {
	
	if (norm && dot_product_ptv(*norm, center, pc) <= 0.0) return 0;
	return (radius == 0.0 || dist_less_than(pc, center, (radius + rsize)));
}


void draw_verts(vector<vertex_t> &verts, unsigned const *ix, int npts, unsigned char shadowed, dqd_params const &p) {

	for (int i = 0; i < npts; ++i) { // much of the frame time is spent here
		assert(ix[i] < verts.size());
		vertex_t &v(verts[ix[i]]);
		
		if (v.c[3] == 0) { // Note: shadowed should agree across all uses of this vertex
			colorRGBA color;
			get_vertex_color(color, p.color, v.p, shadowed, v.n, p.spec, p.in_dlist);
			v.set_c4(color);
		}
		glColor4ubv(v.c);
		v.p.do_glVertex();
	}
	nverts += npts;
	++nsurfaces;
}


#define INTERP_1D(v, s, t, n, i) t*(s*v[2]i + (1.0-s)*v[n-1]i) + (1.0-t)*(s*v[1]i + (1.0-s)*v[0]i)

// s => n1 - n0, t => n3 - n0
template<typename T> inline T interpolate_3d(T const *v, unsigned npts, float s, float t) {

	return T(INTERP_1D(v, s, t, npts, [0]), INTERP_1D(v, s, t, npts, [1]), INTERP_1D(v, s, t, npts, [2]));
}


bool is_above_mesh(point const &pos) {

	if (pos.z > ztop)    return 1;
	if (pos.z < zbottom) return 0;
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

	if (!point_outside_mesh(xpos, ypos)) {
		if (pos.z > mesh_height[ypos][xpos] + DX_VAL) return 1;
		if (pos.z < mesh_height[ypos][xpos] - DX_VAL) return 0;
	}
	return (pos.z > interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1));
}


#define DO_SCALE(tri, s, n) ((tri) ? 1.0-(s)/(float(n)) : 1.0)


bool check_face_containment(point const *const pts, unsigned npts, int dim, int dir, int cobj) { // what about under mesh?

	if (coll_objects[cobj].type != COLL_CUBE) return 0;
	assert((dim >= 0 && dim <= 2) && (dir == 0 || dir == 1));
	point const cent(get_center(pts, npts));
	int const x(get_xpos(cent.x)), y(get_ypos(cent.y));
	if (point_outside_mesh(x, y)) return 0;
	coll_cell const &cell(v_collision_matrix[y][x]);
	unsigned const ncv(cell.cvals.size()), d(2-dim);

	for (unsigned i = 0; i < ncv; ++i) { // test for internal faces to be removed
		unsigned const cid(cell.cvals[i]);
		coll_obj const &c(coll_objects[cid]);
		if (c.type != COLL_CUBE || !c.fixed || c.status != COLL_STATIC || c.platform_id >= 0) continue;
		if ((int)cid == cobj || c.is_semi_trans()) continue;
		if (fabs(c.d[d][!dir] - cent[d]) > TOLER_) continue;
		bool contained(1);

		for (unsigned j = 0; j < npts && contained; ++j) {
			for (unsigned k = 0; k < 2; ++k) {
				unsigned const dk((d+k+1)%3);
				if (pts[j][dk] < (c.d[dk][0]-TOLER_) || pts[j][dk] > (c.d[dk][1]+TOLER_)) contained = 0;
			}
		}
		if (contained) return 1;
	}
	return 0;
}


inline unsigned total_size(vector<vector<int> > const &s) {

	unsigned num(0);
	unsigned const sz(s.size());
	for (unsigned i = 0; i < sz; ++i) num += s[i].size();
	return num;
}


inline void clear_vectors(vector<vector<int> > &s) {

	unsigned const sz(s.size());
	for (unsigned i = 0; i < sz; ++i) s[i].resize(0);
}


inline float scaled_view_dist(point const &camera, point const &center) {

	return p2p_dist(camera, center)/get_zoom_scale();
}


inline void gen_stepsize(point const &camera, point const &center, int step[2], unsigned n[2], float len[2]) {

	float const dist(scaled_view_dist(camera, center));
	if (dist <= DIST_CUTOFF) return;

	for (unsigned d = 0; d < 2; ++d) {
		step[d] = max(int(len[d]*dist*subdiv_size_inv/(n[d]*DIST_CUTOFF) + 0.5), 1);
	}
}


unsigned char test_all_light_val(unsigned lighted, unsigned val) {

	unsigned const num_lights(enabled_lights.size());
	unsigned char res(0);
	assert((val & ~0xF) == 0);

	for (unsigned i = 0; i < num_lights; ++i) {
		if (get_light_val(lighted, i) == val) res |= (1 << i);
	}
	return res;
}


bool check_lv_req(unsigned lighted, unsigned val1, unsigned val2, bool inv) {

	unsigned const num_lights(enabled_lights.size());

	for (unsigned i = 0; i < num_lights; ++i) {
		unsigned const lval(get_light_val(lighted, i));
		if ((lval != val1 && lval != val2) ^ inv) return inv;
	}
	return !inv;
}


float get_mesh_zmax(point const *const pts, unsigned npts) {

	float mesh_ztop(zmin);
	int xmin(MESH_X_SIZE-1), ymin(MESH_Y_SIZE-1), xmax(0), ymax(0);

	for (unsigned i = 0; i < npts; ++i) { // get xy bbox
		int const xv(get_xpos_clamp(pts[i].x)), yv(get_ypos_clamp(pts[i].y));
		xmin = min(xmin, xv); ymin = min(ymin, yv); xmax = max(xmax, xv); ymax = max(ymax, yv);
	}
	for (int i = ymin; i <= ymax; ++i) { // calculate highest mesh point for this quad/triangle
		for (int j = xmin; j <= xmax; ++j) {
			mesh_ztop = max(mesh_ztop, mesh_height[i][j]);
		}
	}
	return mesh_ztop;
}


unsigned determine_shadow_matrix(point const *const pts, vector<unsigned char> &norms, dqd_params &p) {

	unsigned const npts(3 + !p.tri);
	bool under_mesh(1), above_mesh(1), outside_mesh(1), calced_mesh_z(0);
	float mesh_ztop(zmin), qzmin(mesh_ztop + Z_SCENE_SIZE), rsize;

	for (unsigned i = 0; i < npts; ++i) {
		if (is_over_mesh(pts[i]))  outside_mesh = 0;
		if (is_above_mesh(pts[i])) under_mesh   = 0; else above_mesh = 0; // check if vertex is above mesh
		qzmin = min(qzmin, pts[i].z);
	}
	if (outside_mesh) return ALL_LT[2]; // not overlapping mesh - default is unshadowed
	if (under_mesh)   return ALL_LT[1]; // always shadowed

	// find shadowing collision objects
	bool do_mesh_intersect(!above_mesh);
	++nvlight;
	unsigned const n0(p.n[0]), n1(p.n[1]), xstride(n1+1), ystride(n0+1), numverts(xstride*ystride);
	unsigned const num_lights(enabled_lights.size());
	unsigned lighted(0);
	point center;
	polygon_bounding_sphere(pts, npts, 0.0, center, rsize);
	for (unsigned i = 0; i < numverts; ++i) norms[i] = 0; // initialize to 0 (unshadowed)
	assert(num_lights > 0);
	lv_val old_val(0, NULL);
	quad_div qd_old(p.qd);
	qd_old.tag |= QD_TAG_OLD;
	lvmap::iterator const it(coll_objects[p.cobj].lightmap.find(qd_old));
	
	if (it != coll_objects[p.cobj].lightmap.end()) {
		old_val = it->second; // computation can be skipped if light source has not changed since last time
		assert(old_val.status > 0 && old_val.status < max_lighted);
		assert(old_val.n0 == n0 && old_val.n1 == n1);
		coll_objects[p.cobj].lightmap.erase(it);
	}
	for (unsigned L = 0; L < num_lights; ++L) {
		point const &lpos(enabled_lights[L].get_center());
		unsigned const L_mask(1 << L);
		unsigned cur_lighted(get_light_val(old_val.status, L));
		static vector<int> cobjs;
		cobjs.resize(0);
		bool const visible(enabled_lights[L].lights_polygon(center, rsize, (p.double_sided ? NULL : &p.normal)));

		if (cur_lighted == 3) { // partially shadowed cached old value
			assert(old_val.nvals);
			for (unsigned i = 0; i < numverts; ++i) norms[i] |= (old_val.nvals[i] & L_mask);
		}
		else if (cur_lighted == 0 && visible && enabled_lights[L].is_dynamic()) {
			cur_lighted = 2; // all lighted
		}
		else if (cur_lighted == 1 || !visible || coll_pt_vis_test_large(lpos, center, cobjs, p.cobj, rsize, 1, 1, pts, npts, lpos)) {
			cur_lighted = 1; // completely shadowed
			for (unsigned i = 0; i < numverts; ++i) norms[i] |= L_mask;
		}
		else { // calculate variables for mesh shadow tests
			for (unsigned i = 0; i < npts && USE_MESH_INT && !do_mesh_intersect; ++i) { // not always correct
				if (!is_visible_from_light(pts[i], lpos, USE_MESH_INT)) do_mesh_intersect = 1; // test mesh
			}
			if (!calced_mesh_z) {
				mesh_ztop     = get_mesh_zmax(pts, npts);
				calced_mesh_z = 1;
			}
			if (above_mesh && cobjs.empty() && (!USE_MESH_INT || qzmin >= mesh_ztop)) cur_lighted = 2; // block is completely unshadowed

			if (cur_lighted == 0 || do_mesh_intersect) { // block is at least partially shadowed
				int lci(-1);
				unsigned const ncobjs(cobjs.size());
				unsigned const step(min(n1, 6U));
				ntested += ncobjs;
				vector<pt_pair> bounds(ncobjs);

				for (unsigned j = 0; j < ncobjs; ++j) { // line-plane int check optimization
					coll_obj &c(coll_objects[cobjs[j]]);
					point cube_pts[8], *c_pts(NULL);
					unsigned c_npts(0);
					
					if (c.type == COLL_POLYGON && c.thickness <= MIN_POLY_THICK2) {
						c_pts  = c.points;
						c_npts = c.npoints;
					}
					else if (c.type == COLL_CUBE) {
						get_cube_points(c.d, cube_pts);
						c_pts  = cube_pts;
						c_npts = 8;
					}
					if (c_npts == 0) continue;

					for (unsigned k = 0; k < c_npts; ++k) {
						point p_int;
						float t; // unused

						if (line_int_plane(c_pts[k], lpos, pts[0], p.normal, p_int, t, 1)) {
							if (!bounds[j].used) {
								bounds[j].p[0] = bounds[j].p[1] = p_int;
								bounds[j].used = 1;
							}
							else {
								for (unsigned d = 0; d < 3; ++d) {
									bounds[j].p[0][d] = min(bounds[j].p[0][d], p_int[d]);
									bounds[j].p[1][d] = max(bounds[j].p[1][d], p_int[d]);
								}
							}
						}
					} // for k
					for (unsigned d = 0; d < 3; ++d) { // slight adjustment to account for fp errors in intersection math
						bounds[j].p[0][d] -= 0.01*rsize;
						bounds[j].p[1][d] += 0.01*rsize;
					}
				} // for j
				for (unsigned s0 = 0; s0 <= n0; ++s0) {
					int const offset(s0*xstride);
					float const scale(DO_SCALE(p.tri, s0, n0));
					// shift test points by this amount from edges to make polygon meshes not shadow each other's corners
					float const toler(0.001);
					vector3d const vstep(p.dirs[1]*scale*(1.0 - toler)), vstep2(vstep*step);
					point const pos0(pts[0] + p.dirs[0]*(s0*(1.0 - toler) + toler));
					point pos(pos0 + vstep*toler);

					for (unsigned s1 = 0; s1 <= n1; ++s1, pos += vstep) {
						bool coll(0);

						if ((!above_mesh || qzmin < mesh_ztop) && !is_above_mesh(pos)) {
							coll = 1; // under mesh
						}
						else {
							for (unsigned j = 0; j < ncobjs; ++j) { // like playing Battleships
								if (bounds[j].used) { // line-plane int check optimization
									if (pos.x < bounds[j].p[0].x || pos.x > bounds[j].p[1].x) continue;
									if (pos.y < bounds[j].p[0].y || pos.y > bounds[j].p[1].y) continue;
									if (pos.z < bounds[j].p[0].z || pos.z > bounds[j].p[1].z) continue;
								}
								int const cindex(cobjs[j]);
								assert(cindex != p.cobj);
								coll_obj &c(coll_objects[cindex]);

								// test for adjacent shapes from the same plane (different subdivided polygon triangles)?
								if (c.line_intersect(pos, lpos)) { // intersection
									if (j > 0) swap(cobjs[0], cobjs[j]); // move this cobj to the beginning to increase its priority
									if (j > 0) swap(bounds[0], bounds[j]);
									if (cindex != lci) c.shadow_depends.insert(p.cobj);
									lci  = cindex;
									coll = 1;

									for (s1 += step, pos += vstep2; s1 <= n1; s1 += step, pos += vstep2) { // 25% speedup
										if (!c.line_intersect(pos, lpos)) break;
										for (unsigned q = s1-step; q < s1; ++q) norms[offset+q] |= L_mask;
									}
									pos -= vstep2;
									s1  -= step;
									break;
								}
							}
						} // above mesh
						if (do_mesh_intersect && !coll) { // test mesh
							if (!is_visible_from_light(pos, lpos, USE_MESH_INT)) coll = 1; // doesn't work quite right if under mesh
						}
						if (coll) {
							cur_lighted      |= 1;
							norms[offset+s1] |= L_mask;
						}
						else {
							cur_lighted      |= 2;
						}
					} // for s1
				} // for s0
			} // lighted == 0
		} // not completely shadowed
		set_light_val(lighted, cur_lighted, L);
	} // for L
	assert(lighted != 0);
	return lighted;
}


struct vert_color_comp {

	colorRGBA c[4];
};


unsigned draw_quad_div(vector<vertex_t> &verts, unsigned const *ix, dqd_params &p, bool &in_strip) {

	if (!VERTEX_LIGHTING) return ALL_LT[2]; // all unshadowed
	++nquads;
	point pts[4];
	vector3d normals[4];

	for (unsigned i = 0; i < 4; ++i) {
		pts[i]     = verts[ix[i]].p;
		normals[i] = verts[ix[i]].n;
	}
	p.do_calc_params(pts);
	unsigned char *nvals(NULL);
	unsigned lighted(0);
	unsigned const npts(3 + !p.tri), n0(p.n[0]), n1(p.n[1]);
	unsigned const xstride(n1+1), ystride(n0+1), numverts(xstride*ystride), num_lights(enabled_lights.size());
	static vector<unsigned char> norms, nval_buffer;
	if (numverts > norms.size()) norms.resize(numverts);
	bool const lighting_known(require_either_lightval(p.all_lighted, 1, 2));

	if (lighting_known) { // minor optimization
		lighted = p.all_lighted;
	}
	else {
		assert(size_t(p.cobj) < coll_objects.size());
		coll_obj &c_obj(coll_objects[p.cobj]);
		lvmap::const_iterator const it(c_obj.lightmap.find(p.qd));
		bool const cached(it != c_obj.lightmap.end());

		// lighted: 0 = unknown, 1 = shadowed, 2 = unshadowed, 3 = partially shadowed, 4 = hidden
		if (cached) {
			lighted = it->second.status;
			nvals   = it->second.nvals; // NULL if no lights are 0 and no lights are 3
			if (lighting_known) assert(!nvals);
			assert(lighted > 0 && lighted < max_lighted);
			assert(it->second.n0 == n0 && it->second.n1 == n1);
		} // cached
		else { // new quad
			if (p.no_shadow_calc) {
				lighted = ALL_LT[2]; // mark as unshadowed for now
			}
			else {
				if (!p.tri && !c_obj.is_semi_trans() && check_face_containment(pts, npts, p.qd.dim, p.qd.dir, p.cobj)) {
					lighted = ALL_LT[4];
				}
				else {
					lighted = determine_shadow_matrix(pts, norms, p);
					//if (get_light_val(lighted, 0) == 2 || get_light_val(lighted, 0) == 3) {} // add light to this cell

					if (is_partial_shadow(lighted)) {
						if (c_obj.status == COLL_STATIC) {
							nvals = new unsigned char[numverts];
						}
						else { // can we ever get here?
							nval_buffer.resize(numverts);
							nvals = &nval_buffer[0];
						}
						for (unsigned i = 0; i < numverts; ++i) nvals[i] = norms[i];
					}
				} // face not contained
				if (c_obj.status == COLL_STATIC) {
					c_obj.lightmap[p.qd] = lv_val(lighted, n0, n1, nvals);
					if (require_any_lightval(lighted, 2, 3)) c_obj.lighted = COBJ_LIT_TRUE; // some part of the surface is lit
				}
			}
		} // !cached
		if (FULL_MAP_LOOKUP && !p.no_shadow_calc) {
			bool const all_hidden(lighted == ALL_LT[4]);

			if (!all_hidden || (RENDER_PART_INT && 2*p.num_internal++ > p.num_subdiv)) {
				for (unsigned i = 0; i < num_lights; ++i) {
					unsigned const light_val(get_light_val(lighted, i));
					if (light_val != 1) p.all_shad   &= ~(1 << i);
					if (light_val != 2) p.all_unshad &= ~(1 << i);
				}
			}
			if (!all_hidden) p.all_int_surf = 0;
		}
	} // !lighting_known
	assert(lighted < max_lighted);
	bool has_dynamic(0), can_return(1);

	for (unsigned L = 0; L < num_lights; ++L) {
		unsigned const val(get_light_val(lighted, L));
		assert(val != 0);

		if ((val == 2 || val == 3) && L < p.stest.size() && !p.stest[L].empty()) { // dynamic lights
			has_dynamic = 1;
			can_return  = 0;
			break;
		}
		else if (val == 3) {
			can_return  = 0; // partial shadow
		}
	}
	if (can_return) return lighted;
	int step[2] = {p.lod_level, p.lod_level};
	
	if (!USE_DLIST || (ENABLE_DL_LOD && !p.in_dlist)) { // could use multiple display lists for LOD
		gen_stepsize(get_camera_pos(), pts[0], step, p.n, p.len);
		if (USE_DLIST) {for (unsigned i = 0; i < 2; ++i) step[i] = min(3, step[i]);}
	}

	// check for dynamic shadows
	if (has_dynamic && num_lights > 0) {
		point center;
		float rsize, t; // t unused
		polygon_bounding_sphere(pts, npts, 0.0, center, rsize);
		bool const cp_clip(DO_CIRC_PROJ_CLIP && n0 > 1 && n1 > 0);
		unsigned char const* const old_nvals(nvals);
		nvals = &norms[0]; // use a static array to avoid memory allocation issues
		for (unsigned i = 0; i < numverts; ++i) norms[i] = 0;

		for (unsigned L = 0; L < num_lights; ++L) {
			unsigned const L_mask(1 << L);
			unsigned cur_lighted(get_light_val(lighted, L));

			if (cur_lighted == 1 || cur_lighted == 3) {
				for (unsigned i = 0; i < numverts; ++i) {
					norms[i] |= ((cur_lighted == 1) ? L_mask : (old_nvals[i] & L_mask));
				}
			}
			assert(L < p.stest.size());
			unsigned const num_spheres(p.stest[L].size());
			if (num_spheres == 0) continue;
			if (cur_lighted == 1 || cur_lighted == 4)                        continue; // already shadowed or hidden
			if (!enabled_lights[L].lights_polygon(center, rsize, &p.normal)) continue;
			point const &lpos(enabled_lights[L].get_center());
			vector3d const v1(center, lpos);

			for (unsigned i = 0; i < num_spheres; ++i) { // find all potential sphere intersections
				shadow_sphere const &ss(shadow_objs[p.stest[L][i]]);
				point const &spos(ss.pos);
				float const sz(rsize + ss.radius);
				if (!sphere_test_comp(center, spos, v1, sz*sz)) continue;
				if (!ss.test_volume(pts, npts, lpos))           continue;

				// project bounding sphere into polygon plane and use projected circle as a clipping volume
				bool ssval(0);
				point p_int;
				vector3d const v2(lpos, spos);
				float const s2l_dist(v2.mag()), dist_to_light(v1.mag() + ss.radius);
				float const dp(fabs(dot_product(p.normal, v1)/dist_to_light));
				float const expand(max(1.0f, dist_to_light/max(1.0E-6f, (s2l_dist - sz)))/dp);
				float const sr_adj(expand*ss.radius + max(p.len[0]/n0, p.len[1]/n1));
				float const rsq(sr_adj*sr_adj), rtot(expand*rsize + sr_adj);
				unsigned s1_s(0), s1_e(n1);
				
				if (cp_clip) {
					if (!line_int_plane(spos, lpos, pts[0], p.normal, p_int, t, 1)) continue; // this call is required to calculate p_int
					if (!dist_less_than(p_int, center, rtot)) continue;

					if (!p.tri) { // minor performance improvement
						swap(s1_s, s1_e);

						for (unsigned s1 = 0, end_loop1 = 0; !end_loop1; s1 += step[1]) {
							if (s1 >= n1) {s1 = n1; end_loop1 = 1;}
							point const pos(pts[0] + p.dirs[1]*s1), p_end(pos + p.dirs[0]*n0);
							vector3d const vll(p_end, pos), cp(cross_product(vll, vector3d(pos - p_int)));
							
							if (cp.mag_sq() <= rsq*vll.mag_sq()) { // point-line distance => good scanline
								s1_s = min(s1, s1_s);
								s1_e = max(s1, s1_e);
							}
						}
						if (s1_s > s1_e) continue; // skip this light source
						s1_e = min(n1, s1_e+1);
					}
				}
				for (unsigned s0 = 0, end_loop0 = 0; !end_loop0; s0 += step[0]) { // fill in mesh values
					if (s0 >= n0) {s0 = n0; end_loop0 = 1;}
					float const scale(DO_SCALE(p.tri, s0, n0));
					point const pos0(pts[0] + p.dirs[0]*s0);
					vector3d const vstep(p.dirs[1]*scale);

					if (cp_clip && scale > 0.0) {
						vector3d const vll((pos0 + vstep*n1), pos0), cp(cross_product(vll, vector3d(pos0 - p_int)));
						if (cp.mag_sq() > rsq*vll.mag_sq()) continue; // point-line distance => bad scanline
					}
					bool sskip(0), sval(0);
					int const offset(s0*xstride);

					for (unsigned s1 = s1_s, end_loop1 = 0; !end_loop1; s1 += step[1]) {
						if (s1 >= s1_e) {s1 = s1_e; end_loop1 = 1;}
						if (norms[offset+s1] & L_mask) {sskip = 1; continue;} // already shadowed
						point const pos(pos0 + vstep*s1);
						if (cp_clip && p2p_dist_sq(pos, p_int) > rsq) continue;

						if (ss.line_intersect(pos, lpos)) { // lpos is constant - optimize?
							norms[offset+s1] |= L_mask;
							cur_lighted      |= 1;
							sval              = 1;
						}
						else if (sval) break; // reached the end of the sphere in s0 (convex)
					} // for s1
					if (!sval &&  ssval && !sskip) break; // reached the end of the sphere in s1 (convex)
					if ( sval && !ssval) ssval = 1;
				} // for s0
			} // for i
			assert(cur_lighted <= 4);
			set_light_val(lighted, cur_lighted, L);
		} // for L
		assert(lighted < max_lighted);
	} // end dynamic shadows
	// Note: could check if object is completely shadowed, but that rarely happens so it's probably not worth the trouble
	if (!is_partial_shadow(lighted)) return lighted; // all shadowed or none shadowed
	assert(nvals);
	if (in_strip) {glEnd(); in_strip = 0;}
	//if (p.in_dlist) {} // create some textures instead of drawing - QD_TAG_TEXTURE
	bool const more_strips(BETTER_QUALITY && !p.in_dlist && has_dynamic &&
		scaled_view_dist(get_camera_pos(), pts[0]) <= DIST_CUTOFF);
	static vert_color_comp ccomps[256]; // doesn't have to be static
	bool created[256] = {0};

	// render the quads - is there a way to render these using a single textured quad?
	for (unsigned s0 = 0; s0 < n0; s0 += step[0]) {
		unsigned const off1(min(n0, s0+step[0])), o0(off1*xstride), o1(s0*xstride);
		unsigned s_end(off1);

		if (MERGE_STRIPS) { // find the next set of strips that are identical to this one and merge them
			for (bool end_loop = 0; s_end < n0 && !end_loop; s_end += step[0]) {
				int const o2(s_end*xstride);

				for (unsigned s1 = 0; s1 <= n1; s1 += step[1]) {
					if (nvals[o1+s1] != nvals[o2+s1]) {end_loop = 1; break;} // differing element, stop
					if ((s1+step[1]) > n1 && s1 < n1) s1 = n1 - step[1]; // force the last s1 == n1
				}
			}
			s_end = min(s_end, n0);
			if (more_strips && s_end < n0 && s_end > off1) --s_end;
		}
		glBegin(GL_QUAD_STRIP);
		float const scale[2] = {DO_SCALE(p.tri, s_end, n0), DO_SCALE(p.tri, s0, n0)};
		float const s_[2]    = {s_end*p.ninv[0], s0*p.ninv[0]};
		point pos[2] = {(pts[0] + p.dirs[0]*s_end), (pts[0] + p.dirs[0]*s0)}; // s0max, s0min, dirs[0] => p[1] - p[0]

		for (unsigned s1 = 0, end_loop = 0; !end_loop; s1 += step[1]) {
			if (s1 >= n1) {s1 = n1; end_loop = 1;}
			unsigned const s1n(s1-step[1]), s1p(min(n1, s1+step[1]));

			if (s1 == 0 || s1 == n1 || n1 < 3 ||
				nvals[o0+s1] != nvals[o0+s1n] || nvals[o0+s1] != nvals[o0+s1p] ||
				nvals[o1+s1] != nvals[o1+s1n] || nvals[o1+s1] != nvals[o1+s1p])
			{
				float const t_(s1*p.ninv[1]);
				nverts += 2;

				for (int i = 1; i >= 0; --i) {
					vector3d const n(p.use_n ? interpolate_3d(normals, npts, s_[i], t_) : p.normal);
					point const v(pos[i] + p.dirs[1]*(s1*scale[i])); // dirs[1] => p[3] - p[0]
					unsigned char const shadowed(nvals[(i ? o1 : o0)+s1]);

					if (!created[shadowed]) {
						vert_color_comp &cc(ccomps[shadowed]);
						created[shadowed] = 1;

						for (unsigned r = 0; r < npts; ++r) {
							get_vertex_color(cc.c[r], p.color, pts[r], shadowed, normals[r], p.spec, p.in_dlist);
						}
					}
					vert_color_comp const &cc(ccomps[shadowed]);
					colorRGBA a(interpolate_3d(cc.c, npts, s_[i], t_));
					a.alpha = INTERP_1D(cc.c, s_[i], t_, npts, .alpha);
					a.do_glColor();
					v.do_glVertex();
				}
				++nsurfaces;
			} // else continue strip
		}
		s0 = s_end - step[0]; // advance past merged strips
		glEnd();
	}
	return lighted;
}


// currently only handles parallelograms and triangles
// lighted: 0 = unknown, 1 = shad, 2 = unshad, 3 = partially shad, 4 = hidden int surf
//          values per light source (per nibble of unsigned), so max light sources is 32/4 = 8
void draw_quad_tri(point const *pts0, vector3d const *normals0, int npts, int dim, int dir, unsigned face, dqt_params const &q) {

	if (q.color.alpha == 0.0) return; // transparent
	assert(size_t(q.cobj) < coll_objects.size());
	assert(npts == 3 || npts == 4);
	coll_obj &c_obj(coll_objects[q.cobj]);
	point const camera(get_camera_pos());
	if ((display_mode & 0x08) && is_occluded(c_obj.occluders, pts0, npts, camera)) return; // cull the entire face
	unsigned const num_lights(enabled_lights.size());
	quad_div qd2(dim, dir, QD_TAG_GLOBAL, face, 0); // set face to some magic value
	bool no_shadow(q.no_shadow), no_shadow_edge(no_shadow), first_render(0);
	unsigned all_lighted(0), orig_all(0);
	bool const no_shadow_calc(c_obj.platform_id < 0 && nvlight > MAX_CALC_PER_FRAME);

	if (FULL_MAP_LOOKUP && !no_shadow_calc) {
		lvmap::const_iterator const it(c_obj.lightmap.find(qd2));

		if (it != c_obj.lightmap.end()) { // cached
			all_lighted = orig_all = it->second.status;
			assert(all_lighted < max_lighted);
			if (all_lighted == ALL_LT[4]) return; // entire surface is hidden
		}
		else if (check_face_containment(pts0, npts, dim, dir, q.cobj)) { // back facing and hidden
			c_obj.lightmap[qd2] = lv_val(ALL_LT[4]);
			return;
		}
		else { // placeholder - tells us nothing but avoids duplicate face containment check
			c_obj.lightmap[qd2] = lv_val(ALL_LT[0]);
			first_render = 1;
		}
		for (unsigned i = 0; i < num_lights; ++i) {
			unsigned const lval(get_light_val(all_lighted, i));

			if (lval == 0 || (lval == 2 && i < q.stest.size() && !q.stest[i].empty())) {
				set_light_val(all_lighted, 3, i); // unknown or no static shadow, dynamic shadows still possible
			}
		}
		no_shadow_edge = require_either_lightval(all_lighted, 1, 2); // all or none shadowed
	}
	++npolys;
	point pos;
	float rsize;
	polygon_bounding_sphere(pts0, npts, 0.0, pos, rsize);
	assert(rsize > 0.0);
	bool const tri(npts == 3);
	float const subdiv_size_inv2(1.0/L1_SUBDIV_SIZE); // do LOD here???
	bool const is_black(q.color.red == 0.0 && q.color.green == 0.0 && q.color.blue == 0.0);
	bool const use_norms(normals0 != NULL); // for quadrics: spheres and cylinders
	unsigned shifti(0);
	point pts[4];
	vector3d normals[4];
	bool const double_sided(q.double_sided); // || c_obj.is_semi_trans()

	if (DO_ROTATE && npts == 3) { // rotate points on triangles to improve subdivision direction
		float minlen(0.0);

		for (int i = 0; i < npts; ++i) {
			float const len(p2p_dist_sq(pts0[i], pts0[(i+1+(dir ? npts-2 : 0))%npts]));
			if (i == 0 || len < minlen) {minlen = len; shifti = i;}
		}
		shifti = ((shifti+1)%npts);
	}
	for (int i = 0; i < npts; ++i) {
		unsigned const ix((i+npts-shifti)%npts);
		pts[i] = pts0[ix];
		if (use_norms) normals[i] = normals0[ix];
	}
	if (use_norms) { // rotate and mirror the polygon points and normals appropriately
		swap(pts[0], pts[1]);
		swap(normals[0], normals[1]);
		if (!tri) swap(pts[2], pts[3]);
		if (!tri) swap(normals[2], normals[3]);
		swap(normals[1], normals[npts-1]);
		if ( tri) normals[3] = normals[2];
	}
	else if (!dir) {
		swap(pts[1], pts[npts-1]); // acount for reversed point traversal
	}
	if (tri) pts[3] = pts[2];
	unsigned lighted(ALL_LT[0]), n[2];
	float len[2], ninv[2];
	assert(dir == 0 || dir == 1);
	vector3d dirs[2]; // edge step lengths/directions
	calc_params(pts, dirs, len, ninv, n, subdiv_size_inv2);
	vector3d normal(cross_product(dirs[1], dirs[0]).get_norm());

	unsigned shift_bits(1);
	for (unsigned val = n[0]*n[1]+1; val > 0; val >>= 1, ++shift_bits) {}
	assert(face < (1U<<(32-shift_bits)));
	face <<= shift_bits;
	
	if (double_sided && dot_product_ptv(normal, camera, pos) < 0.0) { // viewing the back side
		if (use_norms) {
			for (int i = 0; i < npts; ++i) {
				normals[i].negate();
			}
		}
		normal.negate();
	}
	vector3d const norm(use_norms ? get_center(normals, npts) : normal);
	bool const bf_test(!no_shadow_edge && !double_sided);
	bool back_facing(bf_test), dg_lights(0), has_d_shad(0);
	unsigned char shad(0);

	for (unsigned L = 0; L < num_lights; ++L) {
		if (!enabled_lights[L].lights_polygon(pos, rsize)) continue;
		bool const bf(dot_product_ptv(norm, enabled_lights[L].get_center(), pos) <= 0.0);
		if (!bf && enabled_lights[L].is_dynamic()) dg_lights = 1;
		if (!bf_test) continue;

		if (bf) { // back facing
			shad |= (1 << L);
		}
		else {
			back_facing = 0;
		}
	}
	if (back_facing) no_shadow = no_shadow_edge = 1; // back facing to all light sources
	static vector<vector<int> > stest;
	stest.resize(num_lights);
	
	if (!no_shadow_edge) { // if not facing away from the light (already shadowed)
		if (q.stest.size() == num_lights) { // nonempty dynamic shadow objects
			for (unsigned L = 0; L < num_lights; ++L) {
				vector<int> const &qstest(q.stest[L]);
				if (qstest.empty() || !enabled_lights[L].lights_polygon(pos, rsize, &normal)) continue;
				point const &lpos(enabled_lights[L].get_center());

				for (unsigned s = 0; s < qstest.size(); ++s) { // find all potential sphere intersections
					unsigned const sts(qstest[s]);
					assert(sts < shadow_objs.size());
					point const &spos(shadow_objs[sts].pos);
					float const radius(rsize + shadow_objs[sts].radius);

					if ((double_sided || (dot_product_ptv(norm, spos, pos) > 0.0)) &&
						line_sphere_int_cont(pos, lpos, spos, radius))
					{ // in front of quad
						stest[L].push_back(sts); // point inside test might not be necessary
						has_d_shad = 1;
					}
				}
			}
		}
		else {
			unsigned const tot(total_size(q.stest));
			if (tot > 0) {
				cout << "num_lights = " << num_lights << ", q.stest.size() = " << q.stest.size() << ", tot_sz = " << tot << endl;
				assert(0);
			}
		}
	}
	// Requirements for using a display list:
	// 1. Must have no dynamic shadows
	// 2. Must have no dynamic lights
	// 3. Must not be both specular and lit
	// 4. Must not be part of a quadric (too many LODs)
	float const spec[2] = {c_obj.cp.specular, c_obj.cp.shine};
	bool const is_specular(ALLOW_SPECULAR && spec[0] > 0.0 && !back_facing && all_lighted != ALL_LT[1]);
	bool const use_dlist(USE_DLIST && !first_render && !no_shadow_calc && !is_specular && !q.is_quadric && !dg_lights &&
		!has_d_shad && ((display_mode & 0x10) || !has_dynamic_lights(pts, npts)));
	bool const no_subdiv(no_shadow_edge || is_black);
	unsigned lod_level(1);
	float lod_scale(1.0);

	if (use_dlist) {
		if ((!no_subdiv || LOD_NO_SUBDIV) && (n[0] > 1 || n[1] > 1)) { // might create dlists for the same LOD if there are no shadow edges
			float const dist(scaled_view_dist(camera, pos));
			unsigned const lod(unsigned(len[0]*dist*subdiv_size_inv2/(n[0]*DIST_CUTOFF) + 0.5));
			lod_level = ((lod >= 4) ? 2 : 1); // two LOD levels for display lists
			if (no_subdiv) lod_scale = 1.0/float(lod_level);
		}
		quad_div const qddl(dim, dir, QD_TAG_DLIST, (face + lod_level), shift_bits);
		lvmap::const_iterator const it(c_obj.lightmap.find(qddl));
		bool const found(it != c_obj.lightmap.end());

		if (!found || it->second.status == 0) {
			if (found) c_obj.lightmap.erase(qddl);
			unsigned const dlist(glGenLists(1)); // shouldn't return 0
			assert(glIsList(dlist));
			c_obj.lightmap[qddl] = lv_val(dlist);
			glNewList(dlist, GL_COMPILE_AND_EXECUTE); // mapx has 4067 dlists, 88.5K verts
		}
		else {
			++dlists;
			glCallList(it->second.status);
			return;
		}
	}
	else if (no_subdiv) {
		// use dot_product((camera - pos), norm)/p2p_dist(camera, pos)?
		if (lod_level > 1) { // nothing
		}
		else if (LOD_QUAD_TRIS && (n[0] > 1 || n[1] > 1) && scaled_view_dist(camera, pos) > 2.4*DIST_CUTOFF) {
			lod_scale = 0.5; // half the resolution
		}
		else if (DYNAM_RES_SCALE > 1.0 && (!(display_mode & 0x10) && has_dynamic_lights(pts, npts))) {
			lod_scale = DYNAM_RES_SCALE; // LOD
		}
	}
	++cobj_counter;
	if (lod_scale != 1.0) calc_params(pts, dirs, len, ninv, n, lod_scale*subdiv_size_inv2); // smaller n
	bool const occlusion_test(!use_dlist && (display_mode & 0x08) && !c_obj.occluders.empty());
	bool in_strip(0), occluded(0);
	unsigned const num_subdiv(n[0]*n[1]);
	quad_div qd(dim, dir, QD_TAG_QUAD, face, shift_bits);
	dqd_params params(q, tri, 0, normal, qd, stest, use_dlist, spec, lod_level, num_subdiv, orig_all, no_shadow_calc, use_norms);
	unsigned const nv0(n[0]+1), nv1(n[1]+1);
	static vector<vertex_t> verts;
	verts.resize(nv0*nv1);

	for (unsigned s0 = 0; s0 < nv0; ++s0) {
		point const pt(pts[0] + dirs[0]*s0); // pts[1]
		vector3d_d const vs(dirs[1]*DO_SCALE(tri, s0, n[0])); // pts[3]

		for (unsigned s1 = 0; s1 < nv1; ++s1) {
			unsigned const ix(s0*nv1 + s1);
			verts[ix].c[3] = 0; // initialize
			verts[ix].p    = pt + vs*s1;
			verts[ix].n    = (use_norms ? interpolate_3d(normals, npts, s1*ninv[1], s0*ninv[0]) : normal);
		}
	}
	if (no_subdiv) {
		for (unsigned i = 0; i < num_lights; ++i) {
			if (get_light_val(all_lighted, i) == 1) shad |= (1 << i);
		}
	}
	else if (tri) {
		glBegin(GL_TRIANGLES);
	}
	for (unsigned s0 = 0; s0 < n[0]; ++s0) { // pre-subdivide large quads
		unsigned const ix0(s0*nv1);

		if (occlusion_test && n[0] > 1) { // test this row for occlusion
			point const pts2[4] = {verts[ix0].p, verts[ix0+nv1-1].p, verts[ix0+nv1].p, verts[ix0+nv1+nv1-1].p};

			if (is_occluded(c_obj.occluders, pts2, 4, camera)) {
				occluded = 1;
				continue; // cull a full strip
			}
		}
		if (no_subdiv) {
			glBegin(GL_QUAD_STRIP); // is there any way to safely skip vertices for LOD draw?

			for (unsigned s1 = 0; s1 <= n[1]; ++s1) {
				unsigned const ix(s0*nv1 + s1), ixs[2] = {ix, ix+nv1};
				draw_verts(verts, ixs, 2, shad, params);
			}
			glEnd();
			continue;
		}
		unsigned last_lighted(0), last_ixs[2];

		for (unsigned s1 = 0; s1 < n[1]; ++s1) {
			unsigned const ix(ix0 + s1);
			unsigned ixs[4] = {ix, ix+1, ix+nv1+1, ix+nv1};
			qd.face = face + ix;

			if (tri) {
				bool tri_begin(1);

				if (s0 < n[0]-1) { // not a single triangle at the top
					qd.tag |=  QD_TAG_TRIANGLE;
					unsigned const ixs2[4] = {ixs[0], ixs[2], ixs[3], ixs[3]};
					dqd_params params_tri(q, tri, 1, normal, qd, stest, use_dlist, spec, lod_level, num_subdiv, orig_all, no_shadow_calc, use_norms);
					lighted = draw_quad_div(verts, ixs2, params_tri, tri_begin);
					params.inherit_flags_from(params_tri);
					qd.tag &= ~QD_TAG_TRIANGLE;

					if (!tri_begin) {
						glBegin(GL_TRIANGLES);
						tri_begin = 1;
					}
					else if (lighted != ALL_LT[4]) {
						draw_verts(verts, ixs2, npts, test_all_light_val(lighted, 1), params_tri);
					}
				}
				ixs[3] = ixs[2]; // have to recompute since dir is not constant due to scale
				dqd_params params2(q, tri, 0, normal, qd, stest, use_dlist, spec, lod_level, num_subdiv, orig_all, no_shadow_calc, use_norms);
				lighted = draw_quad_div(verts, ixs, params2, tri_begin);
				if (!tri_begin) {glBegin(GL_TRIANGLES); tri_begin = 1;}
				params.inherit_flags_from(params2);
			}
			else { // quad
				lighted = draw_quad_div(verts, ixs, params, in_strip);
			}
			assert(lighted != 0);

			if (lighted != ALL_LT[4] && !is_partial_shadow(lighted)) {
				unsigned char const cur_shadowed(test_all_light_val(lighted, 1));

				if (tri) {
					draw_verts(verts, ixs, npts, cur_shadowed, params);
				}
				else { // technically, should never be able to directly transition between 1 and 0 unless the shadow falls exactly at the boundary
					unsigned const ixs2[4] = {ixs[3], ixs[0], ixs[2], ixs[1]};

					if (in_strip) {
						if (lighted != last_lighted) { // stutter the last point at the new lighting value so that there's no transition
							draw_verts(verts, last_ixs, 2, cur_shadowed, params);
						}
						draw_verts(verts, ixs2+2, 2, cur_shadowed, params);
					}
					else {
						glBegin(GL_QUAD_STRIP);
						in_strip = 1;
						draw_verts(verts, ixs2, npts, cur_shadowed, params);
					}
					last_lighted = lighted;
					last_ixs[0]  = ixs[2];
					last_ixs[1]  = ixs[1];
				}
			} // not partial shadow && lighted != ALL_LT[4]
		} // for s1
		if (in_strip) {glEnd(); in_strip = 0;}
	} // for s0
	assert(!tri || !in_strip);
	if (!no_subdiv && tri) glEnd();
	if (use_dlist) glEndList();

	if (FULL_MAP_LOOKUP && VERTEX_LIGHTING && !no_shadow_calc && all_lighted == ALL_LT[0] && !occluded) {
		if (!no_subdiv) {
			if (params.all_int_surf) {
				all_lighted = ALL_LT[4]; // free display list?
			}
			else {
				for (unsigned i = 0; i < num_lights; ++i) {
					bool const all_shad((params.all_shad   & (1 << i)) != 0);
					bool const all_unsh((params.all_unshad & (1 << i)) != 0);
					assert(!all_shad || !all_unsh);
					if (all_shad) set_light_val(all_lighted, 1, i);
					if (all_unsh) set_light_val(all_lighted, 2, i);
				}
			}
		}
		c_obj.lightmap[qd2] = lv_val(all_lighted);
	}
	clear_vectors(stest);
}


void draw_quad(float d1a, float d1b, float d2a, float d2b, float d0, int dim, int dir, dqt_params const &q) {

	double const dv(DIV_VAL*L1_SUBDIV_SIZE);
	unsigned const n1(max((unsigned)1, unsigned(fabs(d1b - d1a)/dv)));
	unsigned const n2(max((unsigned)1, unsigned(fabs(d2b - d2a)/dv)));
	double const d1((d1b - d1a)/n1), d2((d2b - d2a)/n2);
	int const dim1((dim+1)%3), dim2((dim+2)%3);
	point pts[4];
	assert(d1 > 0.0 && d2 > 0.0);

	for (unsigned s1 = 0; s1 < n1; ++s1) {
		double const s1d1(d1a + s1*d1);

		for (unsigned s2 = 0; s2 < n2; ++s2) {
			double const s2d2(d2a + s2*d2);

			double v[4][3] = {
				{s1d1,    s2d2,    d0},
				{s1d1,    s2d2+d2, d0},
				{s1d1+d1, s2d2+d2, d0},
				{s1d1+d1, s2d2,    d0}};
			for (unsigned i = 0; i < 4; ++i) {
				pts[i].assign(v[i][dim], v[i][dim1], v[i][dim2]);
			}
			draw_quad_tri(pts, NULL, 4, dim, dir, (s1*n2 + s2), q);
		}
	}
}


void draw_coll_cube(float ar, int do_fill, int cobj, int tid) {

	assert(size_t(cobj) < coll_objects.size());
	coll_obj const &c(coll_objects[cobj]);
	int const sides((int)c.cp.surfs);
	if (sides == EF_ALL) return; // all sides hidden
	bool const back_face_cull(!c.is_semi_trans()); // no alpha
	point const pos(c.points[0]), camera(get_camera_pos());
	bool inside(!back_face_cull);
	bool const textured(tid >= 0);
	float const tscale[2] = {c.cp.tscale, get_tex_ar(tid)*c.cp.tscale};

	if (!inside) { // check if in the camera's view
		float const dist(NEAR_CLIP + CAMERA_RADIUS);
		inside = 1;

		for (unsigned i = 0; i < 3 && inside; ++i) {
			if (camera[i] <= c.d[i][0]-dist || camera[i] >= c.d[i][1]+dist) inside = 0;
		}
	}
	++ncubes;
	assert((size_t)cobj < coll_objects.size());
	dqt_params const q(c.get_ref_pt(), 0, 0, cobj, c.cp.tscale, ar*c.cp.tscale, c.cp.tdx, c.cp.tdy, 0);
	pair<float, unsigned> faces[6];
	for (unsigned i = 0; i < 6; ++i) faces[i].second = i;
	vector3d tex_delta(xoff2*DX_VAL, yoff2*DY_VAL, 0.0);

	if (c.platform_id >= 0) { // make texture scroll with platform
		assert(c.platform_id < (int)platforms.size());
		tex_delta -= platforms[c.platform_id].get_delta();
	}
	if (!back_face_cull) { // semi-transparent
		for (unsigned i = 0; i < 6; ++i) {
			unsigned const dim(i>>1), dir(i&1), d0((dim+1)%3), d1((dim+2)%3);
			point pos;
			pos[dim] = c.d[dim][dir];
			pos[d0]  = 0.5*(c.d[d0][0] + c.d[d0][1]);
			pos[d1]  = 0.5*(c.d[d1][0] + c.d[d1][1]);
			faces[i].first = -p2p_dist_sq(pos, camera); // draw ordered furthest to closest to camera
		}
		sort(faces, (faces+6));
	}
	for (unsigned i = 0; i < 6; ++i) {
		unsigned const fi(faces[i].second), dim(fi>>1), dir(fi&1);
		if ((sides & EFLAGS[dim][dir]) || (!inside && !((camera[dim] < c.d[dim][dir]) ^ dir))) continue;
		unsigned const d0((dim+1)%3), d1((dim+2)%3), t0((2-dim)>>1), t1(1+((2-dim)>0));

		if (textured) {
			float a[4] = {0.0}, b[4] = {0.0};
			a[t0] = tscale[0];
			b[t1] = tscale[1];
			a[3]  = tex_delta[t0]*tscale[0];
			b[3]  = tex_delta[t1]*tscale[1];
			glTexGenfv(GL_S, GL_EYE_PLANE, a);
			glTexGenfv(GL_T, GL_EYE_PLANE, b);
		}
		if (FAST_SHAPE_DRAW) {
			c.cp.color.do_glColor();
			glBegin(GL_QUADS);
			point p;
			p[d0 ] = c.d[d0][0];
			p[d1 ] = c.d[d1][0];
			p[dim] = c.d[dim][dir];
			p.do_glVertex();
			p[d0 ] = c.d[d0][1];
			p.do_glVertex();
			p[d1 ] = c.d[d1][1];
			p.do_glVertex();
			p[d0 ] = c.d[d0][0];
			p.do_glVertex();
			glEnd();
		}
		else {
			draw_quad(c.d[d0][0], c.d[d0][1], c.d[d1][0], c.d[d1][1], c.d[dim][dir], 2-dim, dir, q);
		}
	}
}


bool camera_back_facing(point const *const points, int npoints, vector3d const &normal) {

	return (dot_product_ptv(normal, get_camera_pos(), get_center(points, npoints)) >= 0.0);
}


bool camera_behind_polygon(point const *const points, int npoints, bool &cbf) {

	point const center(get_center(points, npoints)), camera(get_camera_pos());
	vector3d const dirs[2] = {vector3d(points[1], points[0]), vector3d(points[npoints-1], points[0])};
	vector3d const normal(cross_product(dirs[0], dirs[1]));
	cbf = (dot_product_ptv(normal, camera, center) >= 0.0);

	for (unsigned L = 0; L < enabled_lights.size(); ++L) {
		if (!enabled_lights[L].lights_polygon(center, max(dirs[0].mag(), dirs[1].mag()))) continue;
		bool const lbf(dot_product_ptv(normal, enabled_lights[L].get_center(), center) >= 0.0);
		if (!(lbf ^ cbf)) return 0; // camera is on the opposite side of the polygon as all lights
	}
	return 1;
}


void draw_polygon(point const *points, const vector3d *normals, int npoints, vector3d const &norm,
				  int id, int subpoly, dqt_params &q)
{
	if (FAST_SHAPE_DRAW || (npoints != 3 && npoints != 4)) {
		q.color.do_glColor();
		draw_simple_polygon(points, npoints, norm);
		return;
	}
	bool const temp_ns(q.no_shadow); // subdivide based on L1_SUBDIV_SIZE
	bool cbf(0);
	unsigned npts(3);
	q.no_shadow = (subpoly ? 0 : camera_behind_polygon(points, npoints, cbf));

	if (npoints == 4) { // minimize cut length	
		vector3d v01(points[0] - points[1]), v32(points[3] - points[2]);
		vector3d v12(points[1] - points[2]), v03(points[0] - points[3]);

		// test for parallelograms - can be inaccurate for small shapes
		if (!normals &&
			fabs(v01.mag_sq() - v32.mag_sq()) < 1.0E-5*(v01.mag_sq() + v32.mag_sq()) &&
			fabs(v12.mag_sq() - v03.mag_sq()) < 1.0E-5*(v12.mag_sq() + v03.mag_sq()))
		{
			npts = 4;
		}
		else {
			if (!normals) { // is all this complexity worth the trouble?
				v01.normalize(); v32.normalize(); v12.normalize(); v03.normalize();
				bool const p_01_32(p2p_dist_sq(v01, v32) < TOLERANCE), p_12_03(p2p_dist_sq(v12, v03) < TOLERANCE);

				if (p_01_32 || p_12_03) { // trapezoid
					assert(!(p_01_32 && p_12_03)); // should have gotten into the first if conditional in this case?
					point pts[5];

					if (p_01_32) { // rotate points by 1 so that p1p2 | p0p3
						for (unsigned i = 0; i < 4; ++i) pts[i] = points[(i+1)&3];
					}
					else {
						for (unsigned i = 0; i < 4; ++i) pts[i] = points[i];
					}
					if (p2p_dist_sq(pts[0], pts[3]) < p2p_dist_sq(pts[1], pts[2])) { // rotate pts by 2 so that p0p3 > p1p2
						swap(pts[0], pts[2]);
						swap(pts[1], pts[3]);
					}
					pts[4] = pts[1] + (pts[3] - pts[2]);
					draw_quad_tri((pts+1), normals, 4, 2+cbf, !subpoly, id, q); // parallelogram
					pts[2] = pts[4];
					draw_quad_tri( pts,    normals, 3, 0+cbf, !subpoly, id, q); // triangle
					return;
				}
			} // normals test
			float const len[2] = {p2p_dist_sq(points[0], points[2]), p2p_dist_sq(points[1], points[3])}; // diagonals
			point points2[3]   = {points[0], points[1], points[3]};
			point normals3[3], *normals2(NULL);

			if (normals != NULL) {
				normals2    = normals3; // set to actual data
				normals2[0] = normals[0];
				normals2[1] = normals[1];
				normals2[2] = normals[3];
			}
			if (len[0] < len[1]) { // 0-2 cut
				if (normals != NULL) normals2[1] = normals[2];
				points2[1] = points[2];
			}
			else { // 1-3 cut
				if (normals != NULL) ++normals;
				++points; // move pointer
			}
			draw_quad_tri(points2, normals2, 3, 0+cbf, !subpoly, id, q);
		} // parallelogram test
	}
	draw_quad_tri(points, normals, npts, 2+cbf, !subpoly, id, q);
	q.no_shadow = temp_ns;
}


void draw_extruded_polygon(float thick, point const *const points, const vector3d *normals,
						   int npoints, vector3d const &norm, int cobj, int tid)
{
	assert(size_t(cobj) < coll_objects.size());
	coll_obj const &c(coll_objects[cobj]);
	thick = fabs(thick);
	dqt_params q(c.get_ref_pt(), 0, 1, cobj, 0.0, 0.0, 0.0, 0.0, 0); // just set double_sided
	bool const textured(tid >= 0);
	float const tscale[2] = {c.cp.tscale, get_tex_ar(tid)*c.cp.tscale};
	
	if (thick <= MIN_POLY_THICK2) { // double_sided = 0, relies on points being specified in the correct CW/CCW order
		if (textured) setup_polygon_texgen(norm, tscale);
		draw_polygon(points, normals, npoints, norm, 0, 0, q);
		return;
	}
	assert(points != NULL && (npoints == 3 || npoints == 4));
	static vector<point> pts[2];
	gen_poly_planes(points, npoints, norm, thick, pts);
	bool const bfc(!c.is_semi_trans()), cbf(camera_back_facing(&(pts[1].front()), npoints, norm)), back_facing(bfc && cbf);
	unsigned const nsides(unsigned(npoints)+2);
	assert(nsides <= 6);
	pair<int, unsigned> faces[6];
	for (unsigned i = 0; i < nsides; ++i) faces[i] = make_pair(0, i);

	if (!bfc) { // sort by the number of centerlines crossing the surfaces
		point const camera(get_camera_pos());
		point centers[6];

		for (unsigned i = 0; i < 2; ++i) { // front and back
			centers[i] = get_center(&(pts[i].front()), npoints);
		}
		for (int i = 0; i < npoints; ++i) { // sides
			unsigned const ii((i+1)%npoints);
			point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
			centers[i+2] = get_center(side_pts, 4);
		}
		for (unsigned f = 0; f < nsides; ++f) {
			for (unsigned i = 0; i < 2; ++i) { // front and back
				if (i != f && line_poly_intersect((centers[f] - camera), camera, &(pts[i].front()), npoints)) {
					--faces[f].first;
					++faces[i].first;
				}
			}
			for (int i = 0; i < npoints; ++i) { // sides
				unsigned const ii((i+1)%npoints);
				point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
				
				if ((i+2) != f && line_poly_intersect((centers[f] - camera), camera, side_pts, 4)) {
					--faces[f].first;
					++faces[i+2].first;
				}
			}
		}
		sort(faces, (faces+nsides));
	}
	for (unsigned fi = 0; fi < nsides; ++fi) { // draw back to front
		unsigned const s(faces[fi].second);

		if (s < 2) { // draw front and back
			if (bfc && (back_facing ^ (s == 0))) continue;
			vector3d norm2(norm), n2[4];

			if (!s) {
				reverse(pts[s].begin(), pts[s].end());

				if (normals != NULL) {
					for (int i = 0; i < npoints; ++i) n2[i] = normals[npoints-i-1]; // reverse
				}
				norm2.negate();
			}
			if (textured) setup_polygon_texgen(norm2, tscale);
			draw_polygon(&(pts[s].front()), ((normals && !s) ? n2 : NULL), npoints, norm2, s, 0, q); // draw bottom surface
			if (!s) reverse(pts[s].begin(), pts[s].end());
		}
		else { // draw sides
			unsigned const i(s-2), ii((i+1)%npoints);
			point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
			bool cbf2;
			q.no_shadow = camera_behind_polygon(side_pts, 4, cbf2);

			if (!bfc || !cbf2) {
				if (textured) setup_polygon_texgen(get_poly_norm(side_pts), tscale);
				draw_quad_tri(side_pts, NULL, 4, q.no_shadow, 1, i+4, q); // back face cull?
			}
		}
	}
}


void draw_subdiv_cylinder(point const &p1, point const &p2, float radius1, float radius2, int nsides, int nstacks,
						  bool draw_ends, bool no_bfc, int cobj, bool no_lighting, int tid)
{
	assert(radius1 > 0.0 || radius2 > 0.0);
	assert(size_t(cobj) < coll_objects.size());
	coll_obj const &c(coll_objects[cobj]);

	if (FAST_SHAPE_DRAW || no_lighting) {
		c.cp.color.do_glColor();
		draw_fast_cylinder(p1, p2, radius1, radius2, nsides, (tid >= 0), draw_ends);
		return;
	}
	bool const no_clip(no_bfc || c.is_semi_trans());
	point camera(get_camera_pos());
	fix_nsides(nsides);
	point const ce[2] = {p1, p2};
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius1, radius2, nsides, v12));
	dqt_params q(c.get_ref_pt(), 0, no_bfc, cobj, 0.0, 0.0, 0.0, 0.0, 1); // no_bfc == double sided

	for (int S = 0; S < nsides; ++S) { // nsides can change
		int const index(S + (nsides<<8) + 1);
		int const prevS((S-1+nsides)%nsides), nextS((S+1)%nsides);
		point pts[4] = {vpn.p[S<<1], vpn.p[(S<<1)+1], vpn.p[(nextS<<1)+1], vpn.p[nextS<<1]};

		if (no_clip || dot_product_ptv(vpn.n[S], camera, get_center(pts, 4)) > 0.0) { // cull the sides facing away from the camera
			vector3d normals[4], pn_norms[2] = {vpn.n[prevS], vpn.n[nextS]};

			for (unsigned i = 0; i < 2; ++i) {
				normals[(i<<1)]   = vpn.n[S] + pn_norms[i]; // average of two normals
				normals[(i<<1)].normalize(); // is this necessary?
				normals[(i<<1)+1] = normals[(i<<1)];
			}
			if (radius1 == radius2) {
				draw_quad_tri(pts, normals, 4, 0, 0, index, q); // quad
			}
			else if (radius1 == 0.0) {
				draw_quad_tri(pts, normals, 3, 0, 0, index, q); // triangle
			}
			else if (radius2 == 0.0) {
				pts[2] = pts[3]; // duplicate the last point
				draw_quad_tri(pts, normals, 3, 0, 0, index, q); // triangle
			}
			else {
				draw_polygon(pts, normals, 4, up_vector, index, 1, q); // polygon
			}
		}
	} // for S
	if (draw_ends) { // not quite right due to long thin triangles
		if (tid >= 0) { // textured
			float const tscale[2] = {c.cp.tscale, get_tex_ar(tid)*c.cp.tscale};
			setup_polygon_texgen(v12, tscale);
		}
		bool ends_bf[2];
		unsigned const i_end(1 + (radius2 != 0.0));
		q.double_sided = 0;

		for (unsigned i = (radius1 == 0.0); i < i_end; ++i) {
			ends_bf[i] = ((i != 0) ^ (dot_product_ptv(v12, camera, ce[i]) > 0.0));
		}
		for (int S = 0; S < nsides; ++S) { // nsides can change
			int const index(S + (nsides<<8) + 1), nextS((S+1)%nsides); // prevS?

			for (unsigned i = (radius1 == 0.0); i < i_end; ++i) {
				if (!ends_bf[i]) { // facing towards the camera
					point const pts2[4] = {vpn.p[(S<<1)+i], vpn.p[(nextS<<1)+i], ce[i], ce[i]};
					draw_quad_tri(pts2, NULL, 3, 4, i, index, q);
				}
			}
		}
	} // if draw_ends
}


void draw_subdiv_sphere_at(point const &pos, float radius, int ndiv, int cobj, bool no_lighting, int tid) {

	assert(radius > 0.0);
	assert(size_t(cobj) < coll_objects.size());
	coll_obj const &c(coll_objects[cobj]);

	if (FAST_SHAPE_DRAW || no_lighting) {
		c.cp.color.do_glColor();
		draw_subdiv_sphere(pos, radius, ndiv, (tid >= 0), 1);
		return;
	}
	bool const bfc(!c.is_semi_trans());
	fix_nsides(ndiv);
	point camera(get_camera_pos());
	sd_sphere_d sd(pos, radius, ndiv, NULL);
	sd.gen_points_norms();
	point **points   = sd.get_points();
	vector3d **norms = sd.get_norms();
	dqt_params q(c.get_ref_pt(), 0, 0, cobj, 0.0, 0.0, 0.0, 0.0, 1);

	for (int s = 0; s < ndiv; ++s) {
		int const sn((s+1)%ndiv);

		for (int t = 0; t < ndiv; ++t) {
			int const tn(t+1);
			point    const pts[4]     = {points[s][t], points[sn][t], points[sn][tn], points[s][tn]};
			vector3d const normals[4] = {norms[s][t],  norms[sn][t],  norms[sn][tn],  norms[s][tn]};

			if (!bfc || dot_product_ptv(get_center(normals, 4), camera, get_center(pts, 4)) > 0.0) { // back face culling
				int const index(s + (t<<8) + (ndiv<<16) + 1);
				draw_polygon(pts, normals, 4, up_vector, index, 1, q); // polygon segment
			}
		}
	}
}


void clear_all_lightmaps(int mode, unsigned keep) {

	for (unsigned i = 0; i < coll_objects.size(); ++i) {
		coll_objects[i].clear_lightmap(mode, keep);
	}
}


void add_shadow_obj(point const &pos, float radius, int coll_id, bool light_in_front, bool lighted) {

	if (light_in_front && dot_product_ptv(cview_dir, pos, get_camera_pos()) < 0.0) return;
	shadow_objs.push_back(shadow_sphere(pos, radius, coll_id, lighted));
}


void add_coll_shadow_objs() {
	
	RESET_TIME;
	shadow_objs.resize(0);
	unsigned const num_lights(enabled_lights.size());
	point const camera(get_camera_pos());
	bool light_in_front(1);

	if ((camera_mode == 1 || camera_view == 0) && !has_invisibility(CAMERA_ID)) { // shadow the camera even when in the air (but not when dead)
		point camera_pos(camera);
		if (camera_mode == 1 && !spectate) camera_pos.z -= 0.5*camera_zh; // cancel out the z height that was previously added
		shadow_objs.push_back(shadow_sphere(camera_pos, CAMERA_RADIUS, camera_coll_id, 0));
	}
	for (unsigned L = 0; L < num_lights && light_in_front; ++L) {
		if (dot_product_ptv(cview_dir, enabled_lights[L].get_center(), camera) < 0.0) light_in_front = 0;
	}
	if (begin_motion) { // can ignore if behind camera and light in front of camera
		for (int i = 0; i < num_groups; ++i) { // can we simply use the collision objects for this?
			obj_group const &objg(obj_groups[i]);
			if (!objg.enabled || !objg.large_radius()) continue;
			float const radius(object_types[objg.type].radius);
				
			for (unsigned j = 0; j < objg.end_id; ++j) {
				dwobject const &obj(objg.get_obj(j));
				if (obj.disabled() || obj.status == !objg.obj_has_shadow(j)) continue;
				add_shadow_obj(obj.pos, radius, obj.coll_id, light_in_front);
			}
		}
	}
	for (unsigned i = 0; i < weap_cobjs.size(); ++i) {
		unsigned const cid(weap_cobjs[i]);
		if (cid < 0) continue;
		assert(cid < coll_objects.size());
		float brad;
		point center;
		coll_objects[cid].bounding_sphere(center, brad);
		add_shadow_obj(center, brad, cid, light_in_front);
	}
	if (!platforms.empty()) {
		for (unsigned i = 0; i < coll_objects.size(); ++i) {
			coll_obj const &c(coll_objects[i]);
			if (c.disabled() || !c.dynamic_shadows_only()) continue;
			point center;
			float radius;
			c.bounding_sphere(center, radius);
			add_shadow_obj(center, radius, i, light_in_front);
		}
	}
	if (display_mode & 0x0200) d_part_sys.add_cobj_shadows(light_in_front);

	if (TEST_DS_TIME) {
		for (unsigned q = 0; q < 10000; ++q) {
			add_shadow_obj(gen_rand_scene_pos(), object_types[BALL].radius, -1, light_in_front);
		}
	}
	if (VERBOSE_DYNAMIC || TEST_DS_TIME) {PRINT_TIME("Shadow Object Creation");}
	static bool test_all(0);
	static set<unsigned> shadowed;
	unsigned const ncobjs(coll_objects.size());

	if (test_all) {
		for (unsigned i = 0; i < ncobjs; ++i) {
			if (coll_objects[i].status == COLL_STATIC) clear_vectors(coll_objects[i].sobjs);
		}
	}
	else {
		for (set<unsigned>::const_iterator it = shadowed.begin(); it != shadowed.end(); ++it) {
			assert(*it < ncobjs);
			clear_vectors(coll_objects[*it].sobjs);
		}
	}
	test_all = 0;
	shadowed.clear();
	if (VERBOSE_DYNAMIC || TEST_DS_TIME) {PRINT_TIME("Sobjs Reset");}
	vector<int> cobjs;
	unsigned nadded(0);

	for (unsigned L = 0; L < num_lights; ++L) {
		point const &lpos(enabled_lights[L].get_center());
		float const lrad(enabled_lights[L].get_radius());
		float const light_dist((lrad > 0.0) ? lrad : FAR_CLIP);

		for (unsigned s = 0; s < shadow_objs.size(); ++s) {
			point start_pos(shadow_objs[s].pos), end_pos(start_pos);
			vector3d const vl(start_pos, lpos); // object to light vector
			float const v1_mag(vl.mag());
			if (v1_mag > TOLERANCE) end_pos += vl*(light_dist/v1_mag); // extend the end point
			coll_pt_vis_test_large2(start_pos, end_pos, cobjs, shadow_objs[s].cid, shadow_objs[s].radius, 2, 1); // skip dynamic objects and non-drawn objects
			unsigned const ncobjs(cobjs.size());

			for (unsigned i = 0; i < ncobjs; ++i) {
				coll_obj &cobj(coll_objects[cobjs[i]]);
				if (cobj.all_shadowed()) continue; // check for distance from light?
				cobj.sobjs.resize(num_lights);
				cobj.sobjs[L].push_back((int)s);
				
				if (!test_all) {
					shadowed.insert(cobjs[i]);
					if (++nadded > ncobjs) test_all = 1;
				}
			}
			//if (s > 0) draw_sphere_at(shadow_objs[s].pos, shadow_objs[s].radius, 16); // debugging
		}
	}
	if (VERBOSE_DYNAMIC || TEST_DS_TIME) {PRINT_TIME("Shadow Object Addition");}
}


struct light_sphere {

	point pos;
	float radius;
	
	light_sphere(point const &p=all_zeros, float r=0.0) : pos(p), radius(r) {}
	bool operator==(light_sphere const &ls) const {return (pos == ls.pos && radius == ls.radius);}
	bool operator!=(light_sphere const &ls) const {return !operator==(ls);}
};


void init_subdiv_lighting() {

	RESET_TIME;
	static vector<light_sphere> last_lposes;
	static vector<colorRGBA> last_colors;
	unsigned const nlights(enabled_lights.size()), nlights2(last_lposes.size());
	max_lighted = (16ULL << (max(0U, (nlights-1))<<2));
	assert(last_colors.size() == nlights2);
	vector<light_sphere> lposes(nlights);
	vector<colorRGBA> colors(nlights);
	
	for (unsigned L = 0; L < nlights; ++L) {
		bool const dynamic(enabled_lights[L].is_dynamic()); // dynamic lights can move and change color
		lposes[L] = (dynamic ? light_sphere() : light_sphere(enabled_lights[L].get_center(), enabled_lights[L].get_radius()));
		colors[L] = (dynamic ? WHITE          : enabled_lights[L].get_color());
	}
	if (lposes != last_lposes || invalid_shadows == 1) { // clear the lightmap
		unsigned keep_lights(0);
		
		if (!invalid_shadows) {
			for (unsigned L = 0; L < min(nlights, nlights2); ++L) { // can be different sizes (can also use max and test L)
				if (lposes[L] == last_lposes[L]) keep_lights |= (1 << L);
			}
		}
		invalid_shadows = 0;
		clear_all_lightmaps(0, keep_lights); // clear all
		if (VERBOSE_DYNAMIC) {PRINT_TIME("Lightmap Clear");}
	}
	else if (invalid_shadows == 2) {
		invalid_shadows = 0;
		clear_all_lightmaps(2); // update shadows
		PRINT_TIME("Lightmap Bitset and Dlist Clear");
	}
	else if (colors != last_colors) {
		clear_all_lightmaps(1); // clear only dlists/textures
		if (VERBOSE_DYNAMIC) {PRINT_TIME("Lightmap Dlist Clear");}
	}
	last_lposes = lposes;
	last_colors = colors;

	for (unsigned n = 0; n < sizeof(ALL_LT)/sizeof(unsigned); ++n) {
		ALL_LT[n] = 0;

		for (unsigned i = 0; i < nlights; ++i) {
			set_light_val(ALL_LT[n], n, i);
		}
	}
}


