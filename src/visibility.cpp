// 3D World - OpenGL CS184 Computer Graphics Project - Sphere visibility code and shadow driver
// by Frank Gennari
// 4/17/02

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"


// better when moving the sun/moon, with very large number of trees
// worse with sparse cobjs, and requires tree cobjs for tree shadows
bool const DISABLE_SHADOWS  = 0;
int const FAST_LIGHT_VIS    = 1;
float const NORM_VIS_EXTEND = 0.02;


point ocean;
pos_dir_up camera_pdu;

extern bool combined_gu;
extern int window_width, window_height, do_zoom, display_mode, camera_coll_id, DISABLE_WATER;
extern float zmin, zmax, czmin, czmax, zbottom, ztop, sun_rot, moon_rot, NEAR_CLIP, FAR_CLIP;
extern point sun_pos, moon_pos, litning_pos;
extern obj_type object_types[];



// if cobj_ix is not NULL, it will be set if and only if the light is present and the result is false
bool is_visible_to_light_cobj(point const &pos, int light, float radius, int cobj, int skip_dynamic, int *cobj_ix) { // slow...

	point lpos;
	if (!get_light_pos(lpos, light)) return 0;

	if (lpos.z < czmax || pos.z < czmax) {
		int index;

		if (!coll_pt_vis_test(pos, lpos, 1.5*radius, index, cobj, skip_dynamic, 3)) { // test alpha?
			if (cobj_ix) {*cobj_ix = index;}
			return 0; // test collision objects
		}
	}
	if (is_visible_from_light(pos, lpos, 1+FAST_LIGHT_VIS)) return 1; // test mesh
	if (cobj_ix) {*cobj_ix = -1;}
	return 0;
}


bool coll_pt_vis_test(point pos, point pos2, float dist, int &index, int cobj, int skip_dynamic, int test_alpha) {

	float const minz(skip_dynamic ? max(zbottom, czmin) : zbottom);
	float const maxz(skip_dynamic ? czmax : (ztop + Z_SCENE_SIZE));
	// *** note that this will not work with tree and other cobjs if regenerated ***
	if (!do_line_clip_scene(pos, pos2, minz, maxz)) return 1; // assumes pos is in the simulation region
	vector3d const vcf(pos2, pos);
	float const range(vcf.mag());
	if (range < TOLERANCE) return 1; // too close
	return (!check_coll_line((pos + vcf*(dist/range)), pos2, index, cobj, skip_dynamic, test_alpha));
}


void set_camera_pdu() {
	assert(NEAR_CLIP > 0.0 && NEAR_CLIP < FAR_CLIP);
	camera_pdu = pos_dir_up(get_camera_pos(), cview_dir, up_vector, 0.0, NEAR_CLIP, FAR_CLIP);
}

// Note: near and far clip aren't flat planes like a true frustum - they're curved like a spherical section
// this is close enough to a frustum in practice, and makes sphere intersection much easier and faster
// Note: angle is in radians
pos_dir_up::pos_dir_up(point const &p, vector3d const &d, vector3d const &u, float angle_, float n, float f, float a, bool no_zoom)
		: pos(p), dir(d), upv(u), angle(angle_), near_(n), far_(f), A(a), valid(1)
{
	assert(near_ >= 0.0 && far_ > 0.0 && far_ > near_);
	if (A == 0.0) {A = double(window_width)/double(window_height);} // not yet cacluated, use default window
	if (angle == 0.0) {angle = 0.5*TO_RADIANS*PERSP_ANGLE*((do_zoom && !no_zoom) ? 1.0/ZOOM_FACTOR : 1.0);} // not yet cacluated
	tterm = tanf(angle);
	sterm = sinf(angle);
	assert(tterm > 0.0); // angle < 90

	if (A == 1.0) {x_sterm = sterm;} // common/simple case
	else {
		float atan_val(atanf(A*tterm));
		if (atan_val < 0.0) {atan_val += PI;}
		x_sterm = (atan_val/angle)*sterm; // ratio of X-angle to Y-angle
	}
	tterm_sq2_inv = 2.0/(tterm*tterm*(1.0 + A*A));
	orthogonalize_up_dir();
}

void pos_dir_up::orthogonalize_up_dir() {

	assert(dir != zero_vector);
	orthogonalize_dir(upv, dir, upv_, 1);
	cross_product(dir, upv_, cp);
}

bool pos_dir_up::point_visible_test(point const &pos_) const { // simplified/optimized version of sphere_visible_test()

	if (!valid) return 1; // invalid - the only reasonable thing to do is return true for safety
	vector3d const pv(pos_, pos);
	if (dot_product(dir, pv) < 0.0) return 0; // point behind - optimization
	float const dist(pv.mag());
	if (fabs(dot_product(upv_, pv)) > dist*  sterm) return 0; // y-direction (up)
	if (fabs(dot_product(cp,   pv)) > dist*x_sterm) return 0; // x-direction
	return (dist > near_ && dist < far_); // Note: approximate/conservative but fast
}

// view frustum check: dir and upv must be normalized - checks view frustum
bool pos_dir_up::sphere_visible_test(point const &pos_, float radius) const {

	if (!valid) return (radius >= 0.0); // invalid - the only reasonable thing to do is return true for safety
	vector3d const pv(pos_, pos);
	if (dot_product(dir, pv) < 0.0) return (radius > 0.0 && pv.mag_sq() < max(1.0, A*A)*radius*radius*tterm_sq2_inv); // sphere behind - optimization
	float const dist(pv.mag());
	if (fabs(dot_product(upv_, pv)) > (dist*  sterm + radius)) return 0; // y-direction (up)
	if (fabs(dot_product(cp,   pv)) > (dist*x_sterm + radius)) return 0; // x-direction
	return ((dist + radius) > near_ && (dist - radius) < far_); // Note: approximate/conservative but fast
}

template<unsigned N> bool pos_dir_up::pt_set_visible(point const *const pts) const {

	bool npass(0), fpass(0); // near, far

	for (unsigned i = 0; i < N && (!npass || !fpass); ++i) {
		float const dp(dot_product(dir, vector3d(pts[i], pos)));
		npass = (dp > near_);
		fpass = (dp < far_);
	}
	if (!npass || !fpass) return 0;
	vector3d const v [2] = {upv_, cp};
	float    const aa[2] = {sterm*sterm, x_sterm*x_sterm};

	for (unsigned xy = 0; xy < 2; ++xy) { // y, x
		for (unsigned d = 0; d < 2; ++d) { // lo, hi
			float const w(d ? -1.0 : 1.0);
			bool pass(0);
		
			for (unsigned i = 0; i < N && !pass; ++i) {
				vector3d const pv(pts[i], pos);
				float const dp(w*dot_product(v[xy], pv));
				pass = (dp <= 0.0 || dp*dp <= aa[xy]*pv.mag_sq());
			}
			if (!pass) return 0;
		}
	}
	return 1;
}

bool pos_dir_up::cube_visible(cube_t const &cube) const {

	if (!valid) return 1; // invalid - the only reasonable thing to do is return true for safety
	point cube_pts[8];
	cube.get_points(cube_pts);
	return pt_set_visible<8>(cube_pts);
	// Note: if the above call returns true, we could perform a further check for the frustum (all points) to the outside of each plane of the cube
}

// approximate
bool pos_dir_up::projected_cube_visible(cube_t const &cube, point const &proj_pt) const {

	if (!valid) return 1; // invalid - the only reasonable thing to do is return true for safety
	point cube_pts[16];
	cube.get_points(cube_pts);
	float dmax(0.0);
	for (unsigned i = 0; i < 8; ++i) dmax = max(dmax, p2p_dist_sq(cube_pts[i], pos));
	dmax = sqrt(dmax) + far_; // upper bound: (max_dist_to_furstum_origin + far_clip_dist) usually > max_dist_to_frustum_corner

	for (unsigned i = 0; i < 8; ++i) {
		vector3d const dir(cube_pts[i] - proj_pt);
		cube_pts[i+8] = cube_pts[i] + dir*(dmax/dir.mag()); // projected point = dmax*normalized_dir
	}
	return pt_set_visible<16>(cube_pts);
}

bool pos_dir_up::sphere_and_cube_visible_test(point const &pos_, float radius, cube_t const &cube) const {

	if (!sphere_visible_test(pos_,  radius)) return 0; // none of the sphere is visible
	if ( sphere_visible_test(pos_, -radius)) return 1; // all  of the sphere is visible
	return cube_visible(cube);
}

void pos_dir_up::rotate(vector3d const &axis, float angle) { // unused, but could use in model3d_xform_t::apply_inv_xform_to_pdu()

	if (angle == 0.0) return;
	vector3d v[2] = {dir, upv};
	rotate_vector3d_multi(axis, TO_RADIANS*angle, v, 2); // angle is in degrees
	dir = v[0]; upv = v[1]; // renormalize?
	orthogonalize_up_dir();
	// FIXME: what about pos?
}

void pos_dir_up::apply_z_mirror(float zval) {

	pos.z = 2*zval - pos.z;
	dir.z = -dir.z; // mirror
	upv_  = -upv_;
	orthogonalize_up_dir();
}


bool sphere_cobj_occluded(point const &viewer, point const &sc, float radius) {

	if (!have_occluders()) return 0;
	if (radius*radius < 1.0E-6*p2p_dist_sq(viewer, sc)) {return cobj_contained(viewer, sc, &sc, 1, -1);} // small and far away
	point pts[8];
	
	for (unsigned i = 0; i < 8; ++i) { // really only need 4 points
		for (unsigned j = 0; j < 3; ++j) {
			pts[i][j] = sc[j] + ((i&(1<<j)) ? -1.0 : 1.0)*radius;
		}
	}
	return cobj_contained(viewer, sc, pts, 8, -1);
}


bool cube_cobj_occluded(point const &viewer, cube_t const &cube) {

	if (!have_occluders()) return 0;
	point pts[8];
	cube.get_points(pts);
	return cobj_contained(viewer, cube.get_cube_center(), pts, 8, -1);
}


/*
conditions under which an object is viewable:
	0. inside perspective view volume
	1. on correct side of surface mesh
	2. real visibility with mesh and cobjs
	3. cobj vis check with center
	4. cobj vis check with center and top
	5. cobj vis check with all 7 points
	6. cobj vis check with all 7 points, including dynamic objects
*/
// dir must be normalized
bool sphere_in_view(pos_dir_up const &pdu, point const &pos, float radius, int max_level, bool no_frustum_test) {

	if (!no_frustum_test && !pdu.sphere_visible_test(pos, radius)) return 0;
	if (max_level == 0 || world_mode != WMODE_GROUND)              return 1;
	if (!(display_mode & 0x08)) return 1; // occlusion culling disabled
	point const &viewer(pdu.pos);
	
	if (is_over_mesh(viewer)) { // mesh occlusion culling
		float const zmax(pos.z + radius);
		int const cxpos(get_xpos(viewer.x)), cypos(get_ypos(viewer.y));

		if (!point_outside_mesh(cxpos, cypos)) {
			float const above_mesh_dist(viewer.z - mesh_height[cypos][cxpos]);

			if (fabs(above_mesh_dist) > object_types[SMILEY].radius) { // don't use this test when the viewer is close to the mesh
				int const pxpos1(get_xpos(pos.x - radius)), pxpos2(get_xpos(pos.x + radius));
				int const pypos1(get_ypos(pos.y - radius)), pypos2(get_ypos(pos.y + radius));

				if (pxpos1 >= 0 && pypos1 >= 0 && pxpos2 < MESH_X_SIZE && pypos2 < MESH_Y_SIZE) {
					bool const c_above_m(above_mesh_dist > 0.0); // viewer above mesh
					int flag(1);

					for (int i = pypos1; flag && i <= pypos2; ++i) { // works better for small radius spheres
						for (int j = pxpos1; j <= pxpos2; ++j) {
							if ((zmax < mesh_height[i][j]) ^ c_above_m) {flag = 0; break;}
						}
					}
					if (flag) return 0; // case 6
				}
			}
		}
	}
	if (max_level == 1) return 1;
	if (sphere_cobj_occluded (viewer, pos, radius)) return 0; // intersect view with cobjs
	if (!sphere_visible_to_pt(viewer, pos, radius)) return 0; // intersect view with mesh
	if (max_level == 2) return 1;
	
	// do collision object visibility test (not guaranteed to be correct, typically used only with smileys)
	// *** might be unnecessary with real cobj tests ***
	float ext_dist(1.2*object_types[SMILEY].radius);
	if (dist_less_than(viewer, pos, ext_dist)) return 1; // too close to sphere
	int index;
	point qp[5];
	unsigned const nrays((radius == 0.0 || max_level == 3) ? 1 : ((max_level == 4) ? 2 : 5));
	int const skip_dynamic((max_level < 6) ? 1 : 0); // skip dynamic (what about non-drawn?)
	get_sphere_border_pts(qp, pos, viewer, radius, nrays);
	int const cid((pdu.pos == get_camera_pos()) ? camera_coll_id : -1); // what about smiley coll_ids?

	for (unsigned i = 0; i < nrays; ++i) { // can see through transparent objects
		if (coll_pt_vis_test(qp[i], viewer, ext_dist, index, cid, skip_dynamic, 1)) return 1;
	}
	return 0; // case 7
}


int get_light_pos(point &lpos, int light) {

	if      (light == LIGHT_SUN ) {lpos = sun_pos;  return 1;}
	else if (light == LIGHT_MOON) {lpos = moon_pos; return 1;}
	return 0;
}


inline bool back_face_test(int xpos, int ypos, point const &lpos) {

	point pos;
	get_matrix_point(xpos, ypos, pos);
	return (dot_product_ptv(vertex_normals[ypos][xpos], lpos, pos) < 0.0); // back-face culling
}


bool light_visible_from_vertex(int xpos, int ypos, point const &lpos, int fast) {

	point pos;
	get_matrix_point(xpos, ypos, pos);
	pos += vertex_normals[ypos][xpos]*NORM_VIS_EXTEND; // small offset so that line does not intersect with the starting surface
	return (!line_intersect_mesh(lpos, pos, fast));
}


class mesh_shadow_gen {

	float const *mh;
	unsigned char *smask;
	float const *sh_in_x, *sh_in_y;
	float *sh_out_x, *sh_out_y;
	float dist;
	int xsize, ysize;
	vector3d dir;

	void trace_shadow_path(point v1) {
		point v2(v1 + vector3d(dir.x*dist, dir.y*dist, 0.0));
		float const d[3][2] = {{-X_SCENE_SIZE, get_xval(xsize)}, {-Y_SCENE_SIZE, get_yval(ysize)}, {zmin, zmax}};
		if (!do_line_clip(v1, v2, d)) return; // edge case ([zmin, zmax] should contain 0.0)
		int const xa(get_xpos(v1.x)), ya(get_ypos(v1.y)), xb(get_xpos(v2.x)), yb(get_ypos(v2.y));
		int const dx(xb - xa), dy(yb - ya), steps(max(abs(dx), abs(dy)));
		bool const dim(fabs(dir.x) < fabs(dir.y));
		double const xinc(dx/(double)steps), yinc(dy/(double)steps), dir_ratio(dir.z/dir[dim]);
		double x(xa), y(ya);
		bool inited(0);
		point cur(all_zeros);
		int xstart(0), ystart(0), xend(xsize-1), yend(ysize-1);
		if (dir.x <= 0.0) {swap(xstart, xend);}
		if (dir.y <= 0.0) {swap(ystart, yend);}

		for (int k = 0; ; ++k) { // DDA algorithm
			int const xp((int)x), yp((int)y);
			
			if (xp >= 0 && yp >= 0 && xp < xsize && yp < ysize) {
				int const xp1(min(xsize-1, xp+1)), yp1(min(ysize-1, yp+1));
				float const xpi(x - (float)xp), ypi(y - (float)yp);
				float const mh00(mh[yp*xsize+xp]), mh01(mh[yp*xsize+xp1]), mh10(mh[yp1*xsize+xp]), mh11(mh[yp1*xsize+xp1]);
				float const mh((1.0 - xpi)*((1.0 - ypi)*mh00 + ypi*mh10) + xpi*((1.0 - ypi)*mh01 + ypi*mh11));
				point const pt((-X_SCENE_SIZE + DX_VAL*x), (-Y_SCENE_SIZE + DY_VAL*y), mh);

				// use starting shadow height value
				if (sh_in_y != NULL && xp == xstart && sh_in_y[yp] > MESH_MIN_Z) {
					cur.assign(pt.x, pt.y, sh_in_y[yp]);
					inited = 1;
				}
				else if (sh_in_x != NULL && yp == ystart && sh_in_x[xp] > MESH_MIN_Z) {
					cur.assign(pt.x, pt.y, sh_in_x[xp]);
					inited = 1;
				}
				float const shadow_z((pt[dim] - cur[dim])*dir_ratio + cur.z);

				if (inited && shadow_z > pt.z) { // shadowed
					smask[yp*xsize+xp] |= MESH_SHADOW;
					// set ending shadow height value
					if (sh_out_y != NULL && xp == xend) {sh_out_y[yp] = shadow_z;}
					if (sh_out_x != NULL && yp == yend) {sh_out_x[xp] = shadow_z;}
				}
				else {cur = pt;} // update point
				inited = 1;
			}
			else if (k > steps) {break;}
			x += xinc;
			y += yinc;
		} // for k
	}

public:
	mesh_shadow_gen(float const *const h, unsigned char *sm, int xsz, int ysz, float const *shix, float const *shiy, float *shox, float *shoy)
		: mh(h), smask(sm), sh_in_x(shix), sh_in_y(shiy), sh_out_x(shox), sh_out_y(shoy), xsize(xsz), ysize(ysz) {
		assert(mh != NULL && smask != NULL);
	}

	void run(point const &lpos) { // assumes light source directional/at infinity
		assert(smask != NULL);
		dir  = (all_zeros - lpos).get_norm();
		dist = 2.0*XY_SUM_SIZE/sqrt(dir.x*dir.x + dir.y*dir.y);
		#pragma omp parallel sections num_threads(2)
		{
			#pragma omp section
			{
				float const xval(get_xval((dir.x > 0) ? 0 : xsize));
				for (int y = 0; y < 2*ysize; ++y) { // half increments
					trace_shadow_path(point(xval, (-Y_SCENE_SIZE + 0.5*DY_VAL*y), 0.0));
				}
			}
			#pragma omp section
			{
				float const yval(get_yval((dir.y > 0) ? 0 : ysize));
				for (int x = 0; x < 2*xsize; ++x) { // half increments
					trace_shadow_path(point((-X_SCENE_SIZE + 0.5*DX_VAL*x), yval, 0.0));
				}
			}
		}
	}
};


void calc_mesh_shadows(unsigned l, point const &lpos, float const *const mh, unsigned char *smask, int xsize, int ysize,
					   float const *sh_in_x, float const *sh_in_y, float *sh_out_x, float *sh_out_y)
{
	bool const no_shadow(l == LIGHT_MOON && combined_gu);
	bool const all_shadowed(!no_shadow && lpos.z < zmin);
	unsigned char const val(all_shadowed ? MESH_SHADOW : 0);
	for (int i = 0; i < xsize*ysize; ++i) {smask[i] = val;}
	if (no_shadow || FAST_VISIBILITY_CALC == 3) return;
	if (lpos.x == 0.0 && lpos.y == 0.0) return; // straight down = no mesh shadows
	mesh_shadow_gen(mh, smask, xsize, ysize, sh_in_x, sh_in_y, sh_out_x, sh_out_y).run(lpos);
}


void calc_visibility(unsigned light_sources) {

	if (world_mode == WMODE_UNIVERSE) return;
	check_update_global_lighting(light_sources);
	update_sun_and_moon();
}



