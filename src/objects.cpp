// 3D World - OpenGL CS184 Computer Graphics Project
// obj_group and coll_obj class definitions
// by Frank Gennari
// 1/11/06

#include "3DWorld.h"
#include "mesh.h"
#include "transform_obj.h"
#include "physics_objects.h"
#include "collision_detect.h"
#include "shaders.h"
#include "subdiv.h"


float const NDIV_SCALE = 200.0;


extern bool group_back_face_cull, has_any_billboard_coll, begin_motion, fast_transparent_spheres;
extern int display_mode, destroy_thresh, xoff2, yoff2;
extern float temperature, rain_wetness, snow_cov_amt;
extern double tfticks;
extern unsigned ALL_LT[];
extern obj_type object_types[];
extern dwobject def_objects[];
extern vector<texture_t> textures;
extern coll_obj_group coll_objects;
extern platform_cont platforms;
extern vector<obj_draw_group> obj_draw_groups;


// ******************* COLL_OBJ MEMBERS ******************


void coll_obj::init() {

	cp.coll_func = NULL;
	cp.tscale    = 1.0;
	cp.surfs     = 0;
	cp.flags     = 0;
	cp.spec_color= BLACK;
	cp.shine     = 0.0;
	cp.is_emissive = 0;
	norm         = zero_vector;
	npoints      = 0;
	type         = COLL_NULL;
	status       = COLL_UNUSED;
	destroy      = NON_DEST;
	fixed        = 1;
	id           = -1;
	cgroup_id    = -1;
	dgroup_id    = -1;
}


void coll_obj::clear_internal_data() {

	fixed    = 0; // unfix it so that it's actually removed
	cp.surfs = 0;
	occluders.clear();
}


void coll_obj::print_bounds() const {

	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 2; ++j) {cout << d[i][j] << ",";}
		cout << " ";
	}
}


void coll_obj::bb_union(float bb[3][2], int init) {

	for (unsigned i = 0; i < 3; ++i) {
		bb[i][0] = (init ? d[i][0] : min(bb[i][0], d[i][0]));
		bb[i][1] = (init ? d[i][1] : max(bb[i][1], d[i][1]));
	}
}


void coll_obj::setup_internal_state() {

	if (type == COLL_CUBE) {points[0] = get_cube_center();} // set center point
	status = ((cp.flags & COBJ_DYNAMIC) ? (int)COLL_DYNAMIC : (int)COLL_STATIC);
	calc_volume();
	set_npoints();
	calc_bcube();
}


void coll_obj::calc_volume() {

	volume = 0.0;

	switch (type) {
	case COLL_CUBE: // use bbox
		volume = get_volume();
		break;
	case COLL_SPHERE:
		volume = (4.0f/3.0f)*PI*radius*radius*radius;
		break;
	case COLL_POLYGON:
		volume = polygon_area(points, npoints)*thickness; // Note: thickness is 0 for voxel polygons
		break;
	case COLL_CAPSULE:
		volume = (2.0f/3.0f)*PI*(radius*radius*radius + radius2*radius2*radius2); // two hemispheres
		// fallthrough
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		volume += PI*(radius*radius+radius*radius2+radius2*radius2)*p2p_dist(points[0], points[1])/3.0f;
		break;
	case COLL_TORUS:
		volume = (PI*radius2*radius2)*(2.0f*PI*radius);
		break;
	default: assert(0);
	}
}


void coll_obj::calc_bcube() {

	switch (type) {
	case COLL_CUBE:
		break; // already set
	case COLL_SPHERE:
		assert(radius > 0.0);
		set_from_sphere(points[0], radius);
		break;
	case COLL_POLYGON:
		if (norm == zero_vector) {norm = get_poly_norm(points);} // Note: doesn't check for non-coplanar or degenerate polygons
		assert(!is_nan(norm));
		set_from_points(points, npoints); // set cube_t

		for (unsigned p = 0; p < 3; ++p) {
			float const thick(0.5*thickness*fabs(norm[p]) + 1.0E-6);
			d[p][0] -= thick;
			d[p][1] += thick;
		}
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
	case COLL_TORUS:
		get_bounding_cylinder().calc_bcube(*this);
		break;
	case COLL_CAPSULE: {
		assert(radius > 0.0 && radius2 > 0.0);
		cube_t bcube2;
		set_from_sphere(points[0], radius);
		bcube2.set_from_sphere(points[1], radius2);
		union_with_cube(bcube2); // union of bboxes of two end spheres
		break;
	}
	default: assert(0);
	}
}


float coll_obj::calc_min_dim() const {

	float min_dim(min_len());
	if (type == COLL_POLYGON) min_dim = min(min_dim, thickness);
	return min_dim;
}


float coll_obj::get_min_dist_to_pt(point const &pt) const {

	switch (type) {
	case COLL_CUBE:
		return p2p_dist(pt, cube_t::closest_pt(pt));
	case COLL_SPHERE:
		return max(0.0f, (p2p_dist(pt, points[0]) - radius));
	case COLL_POLYGON: {
		point const center(get_center(points, npoints));
		vector3d delta(norm*(0.5*thickness));
		if (dot_product_ptv(delta, pt, center) < 0.0) {delta = -delta;} // choose the side closest to pt
		float dmin_sq(p2p_dist_sq(pt, (center + delta))); // check center point
		for (int i = 0; i < npoints; ++i) {dmin_sq = min(dmin_sq, p2p_dist_sq(pt, (points[i] + delta)));} // check each corner
		return sqrt(dmin_sq); // slightly approximate
	}
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
	case COLL_CAPSULE:
		return max(0.0f, (p2p_dist(pt, get_closest_pt_on_line(pt, points[0], points[1]) - 0.5f*(radius + radius2)))); // approximate
	case COLL_TORUS: {
		float const dp(fabs(dot_product_ptv(norm, pt, points[0])));
		return max(0.0f, (p2p_dist(pt, points[0]) - radius2 - (1.0f - dp)*radius)); // approximate (ri when parallel to axis and ri+ro when perpendicular to axis)
	}
	default: assert(0);
	}
	return 0.0; // never gets here
}


// Note: these functions are intended to be called on coll cubes that are part of platforms,
// but are okay to call on other shapes to get their bounding cube extents; they return the original bounding cube for non-platforms;
// if delta is not aligned with x/y/z axes then the boundary will be an over approximation, which is inefficient but ok
bool coll_obj::is_update_light_platform() const {
	return (platform_id >= 0 && platforms.get_cobj_platform(*this).get_update_light());
}

cube_t coll_obj::get_extended_platform_bcube() const {return (*this + platforms.get_cobj_platform(*this).get_range());}

cube_t coll_obj::get_platform_max_bcube() const {
	if (is_expanded_platform()) return *this; // already expanded to max bcube
	cube_t bcube(*this);
	if (platform_id >= 0) {bcube.union_with_cube(get_extended_platform_bcube());}
	return bcube;
}

cube_t coll_obj::get_platform_min_bcube() const {
	cube_t bcube(*this);

	if (platform_id >= 0) {
		if (!intersects(get_extended_platform_bcube())) {bcube.set_to_zeros();}
		else {bcube.intersect_with_cube(get_extended_platform_bcube());}
	}
	return bcube;
}

void coll_obj::expand_to_platform_max_bounds() {
	if (platform_id < 0 || type != COLL_CUBE || is_expanded_platform()) return; // only unexpanded platform cubes
	union_with_cube(get_extended_platform_bcube());
	cp.flags |= COBJ_EXPANDED_PLATFORM;
}

void coll_obj::unexpand_from_platform_max_bounds() {
	if (platform_id < 0 || type != COLL_CUBE || !is_expanded_platform()) return; // only expanded platform cubes
	platform const &pf(platforms.get_cobj_platform(*this));
	assert(!pf.is_rotation()); // rotations not supported
	vector3d const platform_range(pf.get_range());
	for (unsigned i = 0; i < 3; ++i) {d[i][platform_range[i] > 0.0] -= platform_range[i];} // undo the above transform
	cp.flags &= ~COBJ_EXPANDED_PLATFORM;
}

void expand_or_unexpand_update_light_platform_cobjs(bool unexpand) {
	for (cobj_id_set_t::const_iterator i = coll_objects.platform_ids.begin(); i != coll_objects.platform_ids.end(); ++i) {
		coll_obj &cobj(coll_objects.get_cobj(*i));
		if (!cobj.is_update_light_platform()) continue;
		if (unexpand) {cobj.unexpand_from_platform_max_bounds();} else {cobj.expand_to_platform_max_bounds();}
	}
}

void pre_rt_bvh_build_hook () {expand_or_unexpand_update_light_platform_cobjs(0);}
void post_rt_bvh_build_hook() {expand_or_unexpand_update_light_platform_cobjs(1);}


bool coll_obj::clip_in_2d(float const bb[2][2], float &val, int d1, int d2, int dir) const {

	assert(d1 >= 0 && d1 <= 3 && d2 >= 0 && d2 <= 3 && d1 != d2 && (dir == 0 || dir == 1));
	int d3(0);
	for (int i = 0; i < 3; ++i) {if (i != d1 && i != d2) d3 = i;}

	switch (type) {
	case COLL_CUBE:
		val = d[d3][dir];
		return 1;

	case COLL_SPHERE: {
		cube_t cube;

		for (unsigned i = 0; i < 2; ++i) {
			cube.d[d1][i] = bb[0][i];
			cube.d[d2][i] = bb[1][i];
			cube.d[d3][i] = (i ? SCENE_SIZE[d3] : -SCENE_SIZE[d3]);
		}
		if (!sphere_cube_intersect(points[0], radius, cube)) return 0; //val = d[d3][dir];
		val = points[0][d3] + pow((PI/6.0), (1.0/3.0))*(dir ? radius : -radius); // doesn't seem to help
		return 1;
	}
	case COLL_CYLINDER:
		if (d3 == 2) {
			val = d[d3][dir]; // can be inaccurate, especially if radius1 != radius2 (cone, etc.)
		}
		else { // unused, untested
			int const bb_z((d1 == 2) ? d1 : d2);
			float const z(0.5f*(bb[bb_z][0] + bb[bb_z][1]));
			float const t((z - points[0].z)/(points[1].x - points[0].z));
			float const r(radius + t*(radius2 - radius));
			val = points[0][d3] + sqrt(PI/4.0)*(dir ? r : -r); // approximate with a square of the same cross sectional area
		}
		return 1;

	case COLL_CYLINDER_ROT:
	case COLL_TORUS: // ???
	case COLL_CAPSULE: // ???
		val = d[d3][dir];
		// should really do something with this - can be very inaccurate
		return 1;

	case COLL_POLYGON: {
		bool in_poly(0);

		for (unsigned i = 0; i < (unsigned)npoints && !in_poly; ++i) {
			if (points[i][d1] > bb[0][0] && points[i][d1] < bb[0][1] && points[i][d2] > bb[1][0] && points[i][d2] < bb[1][1]) {
				in_poly = 1; // polygon has a point inside bb
			}
		}
		for (unsigned i = 0; i < 4 && !in_poly; ++i) {
			if (point_in_polygon_2d(bb[0][i>>1], bb[1][i&&(i<3)], points, npoints, d1, d2)) {
				in_poly = 1; // bb has a point inside polygon
			}
		}
		if (!in_poly) return 0;
		val = d[d3][dir];
			
		if (fabs(norm[d3]) > 0.01 && thickness <= MIN_POLY_THICK) { // doesn't work on thick or vertical polygons
			float const dval(-dot_product(norm, points[0]));
			float const cent[2] = {0.5f*float(bb[0][0] + bb[0][1]), 0.5f*float(bb[1][0] + bb[1][1])};
			val = -(cent[0]*norm[d1] + cent[1]*norm[d2] + dval)/norm[d3] + (dir ? 1.0 : -1.0)*0.5*fabs(norm[d3]*thickness);
		}
		return (val >= d[d3][0] && val <= d[d3][1]);
	}
	default: assert(0);
	}
	return 0;
}


void coll_obj::set_npoints() {

	switch (type) {
	case COLL_CUBE:         npoints = 1; break;
	case COLL_SPHERE:       npoints = 1; break;
	case COLL_TORUS:        npoints = 1; break;
	case COLL_CYLINDER:     npoints = 2; break;
	case COLL_CYLINDER_ROT: npoints = 2; break;
	case COLL_CAPSULE:      npoints = 2; break;
	case COLL_POLYGON:      break; // 3 or 4, should be set at another time
	default: assert(0);
	}
}


void coll_obj::set_from_pts(point const *const pts, unsigned npts) {

	//assert(type == COLL_POLYGON); // check the type?
	//type = COLL_POLYGON; // set the type?
	npoints = npts;
	for (unsigned i = 0; i < npts; ++i) {points[i] = pts[i];}
}


bool coll_obj::is_occluded_from_viewer(point const &viewer) const {

	if (!(display_mode & 0x08) || occluders.empty()) return 0;
	if (is_thin_poly()) {return is_occluded(occluders, points, npoints, viewer);}
	point pts[8];
	unsigned const ncorners(get_cube_corners(d, pts, viewer, 0)); // 8 corners allocated, but only 6 used
	return is_occluded(occluders, pts, ncorners, viewer);
}

bool coll_obj::is_cobj_visible() const {
	if (group_id >= 0) return 1; // grouped cobjs are always visible, since the etnire group is drawn at the same time
	return (check_pdu_visible(camera_pdu) && !is_occluded_from_camera());
}
bool coll_obj::check_pdu_visible(pos_dir_up const &pdu) const {

	if (!pdu.valid) return 1;
	point center;
	float brad;
	bounding_sphere(center, brad);
	return pdu.sphere_and_cube_visible_test(center, brad, *this);
}

void coll_obj::register_coll(unsigned char coll_time, unsigned char coll_type_) {
	last_coll = coll_time;
	coll_type = coll_type_;
	has_any_billboard_coll |= is_billboard; // set global state
}

void coll_obj::check_indoors_outdoors() {

	point test_pts[10];
	unsigned npts(0);

	switch (type) {
	case COLL_CUBE:
		test_pts[npts++].assign(0.5f*(d[0][0] + d[0][1]), 0.5f*(d[1][0] + d[1][1]), d[2][1]); // cube center

		if (radius > HALF_DXY) { // large cube
			for (unsigned i = 0; i < 4; ++i) {
				float const xv(0.99*d[0][i&1] + 0.01*d[0][(i&1)^1]), yv(0.99*d[1][i>>1] + 0.01*d[1][(i>>1)^1]);
				test_pts[npts++].assign(xv, yv, d[2][1]); // use points near the top 4 corners
				test_pts[npts  ] = test_pts[0]; // start at top center
				test_pts[npts++][i>>1] = 0.99*d[i>>1][i&1] + 0.01*d[i>>1][(i&1)^1];
			}
		}
		break;
	case COLL_SPHERE:
		test_pts[npts++] = points[0] + vector3d(0.0, 0.0, radius); // top point of the sphere
		break;
	case COLL_CYLINDER_ROT:
	case COLL_CAPSULE:
		if (!is_cylin_vertical()) { // assumed to be horizontal, but probably good enough for any non-vertical angle
			test_pts[npts++] = points[0] + vector3d(0.0, 0.0, radius ); // first end
			test_pts[npts++] = points[1] + vector3d(0.0, 0.0, radius2); // first end
			test_pts[npts++] = 0.5*(test_pts[0] + test_pts[1]); // center point
			break;
		}
		// else fallthrough to vertical cylinder case
	case COLL_CYLINDER:
		test_pts[npts++] = ((points[0].z < points[1].z) ? points[1] : points[0]); // center of top circle
		if (type == COLL_CAPSULE) {test_pts[0].z += ((points[0].z < points[1].z) ? radius2 : radius);}
		break;
	case COLL_TORUS:
		if (has_z_normal()) { // Note: +z torus only - use 4 points on top
			for (unsigned i = 0; i < 2; ++i) {
				test_pts[npts++] = points[0] + vector3d((i ? radius : -radius), 0.0, radius2);
				test_pts[npts++] = points[0] + vector3d(0.0, (i ? radius : -radius), radius2);
			}
		}
		else {
			vector3d top_dir;
			orthogonalize_dir(plus_z, norm, top_dir, 1); // Note: not valid for norm == +/- z
			test_pts[npts++] = points[0] + top_dir*radius + vector3d(0.0, 0.0, radius2); // should be the point with highest z-value
		}
		break;
	case COLL_POLYGON:
		assert(npoints <= 4);
		if (thickness > MIN_POLY_THICK) { // extruded polygon
			vector3d const delta(norm*(0.5*thickness));
			for (int i = 0; i < npoints; ++i) {
				test_pts[npts++] = points[i] - delta; // neg side
				test_pts[npts++] = points[i] + delta; // neg side
			}
		}
		else {
			test_pts[npts++] = get_center_pt();
			if (radius > HALF_DXY) { // large polygon, use all points
				for (int i = 0; i < npoints; ++i) {test_pts[npts++] = points[i];}
			}
		}
		break;
	default: assert(0);
	}
	bool all_blocked(1);
	int cindex(-1); // unused
	for (unsigned i = 0; i < npts && all_blocked; ++i) {all_blocked &= check_coll_line(test_pts[i], point(test_pts[i].x, test_pts[i].y, czmax), cindex, id, 1, 0, 1, 0);}
	if (all_blocked) {cp.flags |= COBJ_IS_INDOORS;} else {cp.flags &= ~COBJ_IS_INDOORS;} // only indoors if all rays are blocked
}

void coll_obj::setup_cobj_sc_texgen(vector3d const &dir, shader_t &shader) const { // dir does not need to be normalized

	float const tscale(((cp.tid >= 0) ? get_tex_ar(cp.tid) : 1.0)*cp.tscale);
	int const dim(::get_max_dim(dir));
	point p1, p2;
	
	for (int i = 0; i < 3; ++i) {
		p1[i] = (i == dim) ? cp.tscale : 0.0;
		p2[i] = (i == dim) ? 0.0       : tscale;
	}
	if (cp.swap_txy()) {swap(p1, p2);}
	setup_texgen_full(p1.x, p1.y, p1.z, dot_product(p1, texture_offset), p2.x, p2.y, p2.z, dot_product(p2, texture_offset), shader, 1);
}

void coll_obj::get_sphere_cylin_tparams(vector3d const &dir, texgen_params_t &tp, float tscale0_mult) const {

	float const tscale0(tscale0_mult*cp.tscale), tscale1(((cp.tid >= 0) ? get_tex_ar(cp.tid) : 1.0)*cp.tscale);
	int const dim(::get_max_dim(dir));
	point p[2];

	for (int i = 0; i < 3; ++i) {
		p[0][i] = (i == dim) ? tscale0 : 0.0;
		p[1][i] = (i == dim) ? 0.0     : tscale1;
	}
	for (unsigned i = 0; i < 2; ++i) {
		bool const d((i != 0) ^ cp.swap_txy());
		UNROLL_3X(tp.st[d][i_] = tscale1*p[i][i_];)
		tp.st[d][3] = tscale1*dot_product(texture_offset, p[i]);
	}
}

void draw_subdiv_sphere_back_to_front(point const &pos, float radius, int ndiv, int texture) {

	glEnable(GL_CULL_FACE);
	sd_sphere_d sd(pos, radius, ndiv);
	sd.gen_points_norms_static();
	sd.draw_subdiv_sphere(all_zeros, texture, 1, nullptr, nullptr, nullptr, 0.0, 0.0, 1.0, 0.0, 1.0, 1); // back_to_front=1
	glDisable(GL_CULL_FACE);
}

void coll_obj::draw_cobj(unsigned &cix, int &last_tid, int &last_group_id, shader_t &shader, cobj_draw_buffer &cdb, int reflection_pass) const {

	assert(!disabled());
	bool const in_group(group_id >= 0), same_group(group_id == last_group_id), start_group(in_group && !same_group);
	// we want everything to be textured for simplicity in code/shaders,
	// so if there is no texture specified just use a plain white texture
	int const tid(get_tid());

	// process groups
	if (same_group) {
		assert(!in_group || tid == last_tid); // Note: each texture must be in its own group (see where group_ids are set in build_world.cpp)
	}
	else { // group changed
		end_group(last_group_id);
	}
	float const prev_refract_ix(cdb.get_last_layer().refract_ix);
	cdb.on_new_obj_layer(cp); // may flush/draw

	if (!in_group || start_group) { // should be the same across groups
		shader.set_material(cp);
		if (cp.refract_ix != prev_refract_ix) {shader.set_refract_ix(cp.refract_ix);}
	}
	if (cp.is_emissive) {assert(!in_group); cdb.flush(); shader.add_uniform_float("emissive_scale", 1.0);} // Note: slow (causes flush)

	if (tid != last_tid) {
		bool const textured(select_texture(tid));
		assert(textured);
		last_tid = tid;
	}
	last_group_id = group_id;
	bool const draw_as_wet(is_wet() && cp.metalness < 0.5); // wetness doesn't affect metal material properties

	if (in_group) {
		if (start_group) { // start rendering a new group
			// Note: this can be done between the begin/end of a group and is more efficient in some cases, but less efficient in others
			if (group_back_face_cull) {glEnable(GL_CULL_FACE);}
			vector3d tex_dir(0,0,0);
			tex_dir[::get_max_dim(norm)] = 1.0;
			set_poly_texgen(tid, tex_dir, shader);
			assert((unsigned)group_id < obj_draw_groups.size());
			if (obj_draw_groups[group_id].begin_render(cix)) return; // cix may change here
		}
		if (!obj_draw_groups[group_id].vbo_enabled()) { // don't do view frustum culling when creating the vbo
			if (group_back_face_cull && dot_product_ptv(norm, get_camera_pos(), points[0]) < 0.0) return;
			if (!is_cobj_visible()) return;
		}
		assert(is_thin_poly()); // thin triangle/quad
		obj_draw_groups[group_id].add_draw_polygon(points, norm, npoints, cix);
		return;
	}
	else if (int(draw_as_wet) != cdb.is_wet) { // check raining state when not in a group
		cdb.flush();
		cdb.is_wet = draw_as_wet;
		shader.add_uniform_float("wet_effect",   (cdb.is_wet ? rain_wetness : 0.0)); // either it's wet or it's not
		shader.add_uniform_float("reflectivity", (cdb.is_wet ? rain_wetness : 1.0)); // either it's partially wet or dry and reflective
	}
	switch (type) {
	case COLL_CUBE:
		draw_coll_cube(tid, cdb);
		break;

	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
	case COLL_TORUS:
	case COLL_CAPSULE: {
		float const scale(NDIV_SCALE*get_zoom_scale()*(reflection_pass ? 0.5 : 1.0));
		float const size(scale*sqrt(((max(radius, radius2) + 0.002)/distance_to_camera(get_closest_pt_on_line(get_camera_pos(), points[0], points[1])))));
		int ndiv(min(N_CYL_SIDES, max(4, (int)size)));
		bool const use_tcs(cp.tscale == 0.0); // even if not textured
		vector3d const dir(points[1] - points[0]);
		if (!use_tcs) {setup_cobj_sc_texgen(dir, shader);} // use texgen

		if (type == COLL_CAPSULE) {
			ndiv &= (~1); // make sure ndiv is even
			draw_fast_cylinder(points[0], points[1], radius, radius2, ndiv, use_tcs, 0, 0, nullptr, 1.0, (use_tcs ? 0.5 : 0.0)); // Note: using texgen, not textured
			float const r[2] = {radius, radius2};

			if (is_cylin_vertical()) { // is the hemisphere optimization worth the trouble?
				bool const d(points[1].z > points[0].z);
				draw_subdiv_sphere_section(points[ d], r[ d], ndiv, use_tcs, 0.0, 1.0, 0.0, 0.5);
				draw_subdiv_sphere_section(points[!d], r[!d], ndiv, use_tcs, 0.0, 1.0, 0.5, 1.0);
			}
			else {
				for (unsigned d = 0; d < 2; ++d) {
					// if using tex coords, rotate the sphere so that textures line up in at least one dim
					// Note: when normal maps are enabled, rotation doesn't correclty rotate the normals
					if (use_tcs && cp.normal_map < 0) {
						fgPushMatrix();
						translate_to(points[d]);
						rotate_from_v2v(dir, plus_z);
						rotate_to_dir(dir, 90.0);
						draw_subdiv_sphere(all_zeros, r[d], ndiv, use_tcs, 1);
						fgPopMatrix();
					}
					else {draw_subdiv_sphere(points[d], r[d], ndiv, use_tcs, 1);}
				}
			}
		}
		else if (type == COLL_TORUS) {
			float const tscale((cp.tscale == 0.0) ? 1.0 : cp.tscale); // always textured, use tscale of 1.0 if not set
			draw_rot_torus(points[0], norm, radius2, radius, ndiv, 3*ndiv/2, tscale, tscale);
		}
		else { // cylinder (Note: use_tcs is usually 0 - using texgen, not textured)
			bool const draw_ends(!(cp.surfs & 1));
			draw_fast_cylinder(points[0], points[1], radius, radius2, ndiv, use_tcs, (draw_ends && use_tcs), !draw_ends);
			if (draw_ends && !use_tcs) {draw_cylin_ends(tid, ndiv, cdb);}
		}
		break;
	}
	case COLL_SPHERE: {
		float const scale(NDIV_SCALE*get_zoom_scale()*(reflection_pass ? 0.5 : 1.0)), size(scale*sqrt((radius + 0.002)/distance_to_camera(points[0])));
		int ndiv(min(N_SPHERE_DIV, max(5, (int)size)));
		if (is_reflective()) {ndiv *= 2;} // higher tessellation for reflective spheres; required for high quality refraction
		if (reflection_pass) {min_eq(ndiv, 3*N_SPHERE_DIV/2);} // limit to MAX_SPHERE_VBO_NDIV for reflection pass to ensure we can use VBOs
		bool const use_tcs(cp.tscale == 0.0);
		if (!use_tcs) {setup_cobj_sc_texgen(plus_z, shader);} // use texgen
		//draw_cube_mapped_sphere(points[0], radius, ndiv/2, use_tcs);
		if (cp.light_atten > 0.0) {glEnable (GL_CULL_FACE);}
		// Note: if using textures (and using MVM in the shader), draw_sphere_vbo() is faster and can be used
		// Note: semi-transparent spheres must not use tex coords because the translate in draw_sphere_vbo() produces the wrong indir lighting vpos;
		// the transparent pass is sorted back-to-front and not grouped by texgen vs. tex-coord shader;
		// however, we allow VBO drawing when fast_transparent_spheres=1 because this will force use_mvm=1 in the untextured/no texgen case
		bool const semi_trans(is_semi_trans());
		if ((use_tcs || fast_transparent_spheres) && !(semi_trans && cp.tid < 0)) {draw_sphere_vbo(points[0], radius, ndiv, use_tcs);}
		else if (semi_trans) { // draw back-to-front
			if (fast_transparent_spheres) {draw_sphere_vbo_back_to_front(points[0], radius, ndiv, use_tcs);}
			else {draw_subdiv_sphere_back_to_front(points[0], radius, ndiv, use_tcs);}
		}
		else {draw_subdiv_sphere(points[0], radius, ndiv, use_tcs, 1);}
		if (cp.light_atten > 0.0) {glDisable(GL_CULL_FACE);}
		break;
	}
	case COLL_POLYGON:
		draw_extruded_polygon(tid, cdb);
		break;
	default: assert(0);
	}
	if (cp.is_emissive) {cdb.flush(); shader.add_uniform_float("emissive_scale", 0.0);} // reset to default
}


void get_shadow_cube_triangle_verts(vector<vert_wrap_t> &verts, cube_t const &c, int eflags) {

	point pts[4];

	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d0(cube_dim_table[0][i]), d1(cube_dim_table[1][i]);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (eflags & EFLAGS[i][j]) continue; // back facing or disabled
			pts[0][i ] = pts[1][i ] = pts[2][i] = pts[3][i] = c.d[i][j];
			pts[0][d0] = pts[3][d0] = c.d[d0][0];
			pts[1][d0] = pts[2][d0] = c.d[d0][1];
			pts[0][d1] = pts[1][d1] = c.d[d1][0];
			pts[2][d1] = pts[3][d1] = c.d[d1][1];
			if (!j) {swap(pts[0], pts[3]); swap(pts[1], pts[2]);}
			for (unsigned k = 0; k < 6; ++k) {verts.push_back(pts[quad_to_tris_ixs[k]]);}
		}
	}
}


void get_cylinder_triangles(vector<vert_wrap_t> &verts, point const &p1, point const &p2, float radius1, float radius2, int ndiv, int draw_ends) {

	assert(ndiv > 0);
	assert(radius1 > 0.0 || radius2 > 0.0);
	point const ce[2] = {p1, p2};
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius1, radius2, ndiv));

	if (radius2 == 0.0) { // cone (triangles)
		for (unsigned s = 0; s < (unsigned)ndiv; ++s) {
			verts.push_back(vpn.p[(s<<1)+1]);
			verts.push_back(vpn.p[(((s+1)%ndiv)<<1)+0]);
			verts.push_back(vpn.p[(s<<1)+0]);
		}
	}
	else {
		for (unsigned S = 0; S <= (unsigned)ndiv; ++S) { // triangle strip
			for (unsigned d = 0; d < 2; ++d) {
				if ((2*S+d) > 2) {tri_strip_push(verts);}
				verts.push_back(vpn.p[((S%ndiv)<<1)+d]);
			}
		}
	}
	if (draw_ends) {
		float const r[2] = {radius1, radius2};

		for (unsigned i = 0; i < 2; ++i) { // triangle fan
			if (r[i] == 0.0) continue;
			unsigned const center_ix(verts.size());
			verts.push_back(ce[i]);

			for (unsigned S = 0; S <= (unsigned)ndiv; ++S) {
				if (S > 1) { // tri_fan_push
					verts.push_back(verts[center_ix]);
					verts.push_back(verts[verts.size()-2]);
				}
				verts.push_back(vpn.p[((S%ndiv)<<1)+i]);
			}
		}
	}
}


void get_torus_triangles(vector<vert_wrap_t> &verts, point const &center, vector3d const &dir, float ro, float ri, int ndiv) {

	assert(ro > 0.0 || ri > 0.0);
	float const ds(TWO_PI/ndiv), cds(cos(ds)), sds(sin(ds));
	vector<float> const &sin_cos(gen_torus_sin_cos_vals(ndiv));
	vector3d vab[2];
	get_ortho_vectors(dir, vab);
	
	for (unsigned s = 0; s < (unsigned)ndiv; ++s) { // outer
		float const theta(s*ds), ct(cos(theta)), st(sin(theta)), ct2(ct*cds - st*sds), st2(st*cds + ct*sds);
		point const pos [2] = {(vab[0]*ct + vab[1]*st), (vab[0]*ct2 + vab[1]*st2)};
		point const vpos[2] = {(center + pos[0]*ro), (center + pos[1]*ro)};

		for (unsigned t = 0; t <= (unsigned)ndiv; ++t) { // inner
			unsigned const t_((t == (unsigned)ndiv) ? 0 : t);
			float const cp(sin_cos[(t_<<1)+0]), sp(sin_cos[(t_<<1)+1]);
			
			for (unsigned i = 0; i < 2; ++i) {
				if ((2*t+i) > 2) {tri_strip_push(verts);}
				verts.push_back(vpos[1-i] + (pos[1-i]*sp + dir*cp)*ri);
			}
		}
	} // for s
}


void get_polygon_triangles(vector<vert_wrap_t> &verts, point const *const points, int npoints) {

	assert(npoints == 3 || npoints == 4);
	for (int i = 0; i < ((npoints == 3) ? 3 : 6); ++i) {verts.push_back(points[quad_to_tris_ixs[i]]);} // 1-2 triangles
}


void get_extruded_polygon_triangles(vector<vert_wrap_t> &verts, float thick, point const *const points, int npoints) {

	assert(points != NULL && (npoints == 3 || npoints == 4));
	thick = fabs(thick);
	
	if (thick <= MIN_POLY_THICK) {
		get_polygon_triangles(verts, points, npoints);
	}
	else {
		point pts[2][4];
		gen_poly_planes(points, npoints, get_poly_norm(points), thick, pts);
		for (unsigned d = 0; d < 2; ++d) {get_polygon_triangles(verts, pts[d], npoints);} // top/bottom surface
	
		for (int i = 0; i < npoints; ++i) { // sides
			int const ii((i+1)%npoints);
			point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
			get_polygon_triangles(verts, side_pts, 4);
		}
	}
}


void coll_obj::get_shadow_triangle_verts(vector<vert_wrap_t> &verts, int ndiv, bool skip_spheres, unsigned char eflags) const {

	switch (type) {
	case COLL_CUBE: {
		get_shadow_cube_triangle_verts(verts, *this, (cp.surfs | eflags));
		break;
	}
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		get_cylinder_triangles(verts, points[0], points[1], radius, radius2, ndiv, !(cp.surfs & 1));
		break;
	case COLL_TORUS:
		get_torus_triangles(verts, points[0], norm, radius, radius2, ndiv);
		break;
	case COLL_CAPSULE:
		if (!skip_spheres) {
			get_sphere_triangles(verts, points[0], radius,  ndiv); // some extra triangles (only need hemisphere), but okay
			get_sphere_triangles(verts, points[1], radius2, ndiv);
		}
		get_cylinder_triangles(verts, points[0], points[1], radius, radius2, ndiv, 0);
		break;
	case COLL_SPHERE:
		if (!skip_spheres) {get_sphere_triangles(verts, points[0], radius, ndiv);}
		break;
	case COLL_POLYGON:
		get_extruded_polygon_triangles(verts, thickness, points, npoints);
		break;
	default: assert(0);
	}
	assert((verts.size()%3) == 0);
}


bool obj_layer::has_alpha_texture() const {
	return (tid >= 0 && textures[tid].has_alpha());
}

bool coll_obj::is_player() const { // sort of a hack
	return ((cp.coll_func == camera_collision) || (cp.coll_func == smiley_collision));
}

bool coll_obj::is_invis_player() const { // sort of a hack
	
	if (cp.coll_func == smiley_collision && has_invisibility(cp.cf_index)) return 1;
	if (cp.coll_func == camera_collision && has_invisibility(-1))          return 1;
	return 0;
}

bool coll_obj::truly_static() const {

	if (may_be_dynamic())                      return 0;
	if (cp.flags & COBJ_DESTROYABLE)           return 0;
	if (cp.cobj_type == COBJ_TYPE_MODEL3D)     return 1;
	if (cp.cobj_type == COBJ_TYPE_VOX_TERRAIN) return 0; // ???
	return (destroy < max((int)SHATTERABLE, destroy_thresh+1));
}

void coll_obj::bounding_sphere(point &center, float &brad) const {

	switch (type) {
	case COLL_CUBE:
	case COLL_SPHERE:
		center = points[0];
		brad   = radius;
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		cylinder_bounding_sphere(points, radius, radius2, center, brad);
		break;
	case COLL_TORUS:
		center = points[0];
		brad   = radius + radius2;
		break;
	case COLL_CAPSULE:
		center = get_center_n2(points);
		brad   = p2p_dist(center, points[0]) + max(radius, radius2); // conservative
		break;
	case COLL_POLYGON:
		center = get_center(points, npoints);
		brad   = radius;
		break;
	default:
		assert(0);
	}
}

cylinder_3dw coll_obj::get_bounding_cylinder() const { // Note: only valid for torus, capsule, and cylinder

	assert(radius >= 0.0 && radius2 >= 0.0);
	
	switch (type) {
	case COLL_TORUS:
		return cylinder_3dw(points[0]-norm*radius2, points[0]+norm*radius2, radius+radius2, radius+radius2);
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
	case COLL_CAPSULE:
		return cylinder_3dw(points[0], points[1], radius, radius2);
	default: assert(0);
	}
	return cylinder_3dw(); // never gets here
}

bool coll_obj::use_tex_coords() const {
	if (was_a_cube() && cp.tid >= 0) return 1; // textured cube
	return (cp.tscale == 0.0 && (is_cylinder() || type == COLL_SPHERE || type == COLL_CAPSULE || type == COLL_TORUS)); // default curved shape
}

point coll_obj::get_center_pt() const {

	switch (type) {
	case COLL_CUBE:         return points[0];
	case COLL_SPHERE:       return points[0];
	case COLL_TORUS:        return points[0];
	case COLL_CYLINDER:     return point(points[0].x, points[0].y, 0.5f*(points[0].z + points[1].z));
	case COLL_CYLINDER_ROT: return get_center_n2(points);
	case COLL_CAPSULE:      return get_center_n2(points); // approximate
	case COLL_POLYGON:      return get_center(points, npoints);
	default: assert(0);
	}
	return all_zeros; // never gets here
}

float coll_obj::get_light_transmit(point v1, point v2) const {

	if (type != COLL_CUBE)        return 1.0; // only cubes are supported due to clipping issues
	if (cp.light_atten == 0.0)    return 1.0;
	if (!do_line_clip(v1, v2, d)) return 1.0;
	return exp(-cp.light_atten*p2p_dist(v1, v2));
}

bool coll_obj::check_poly_billboard_alpha(point const &p1, point const &p2, float t) const {
	return (!has_poly_billboard_alpha() || !is_billboard_texture_transparent(points, (p1 + (p2 - p1)*t), cp.tid));
}

bool sphere_int_capsule(point const &sc, float sr, point const &cp1, point const &cp2, float cr1, float cr2) {
	if (dist_less_than(sc, cp1, (cr1+sr)) || dist_less_than(sc, cp2, (cr2+sr))) return 1;
	return sphere_intersect_cylinder(sc, sr, cp1, cp2, cr1, cr2);
}

int coll_obj::contains_point(point const &pos) const {
	
	if (!contains_pt(pos)) return 0; // test bounding cube
	switch (type) {
	case COLL_CUBE:         return 1; // if bcube contains the point, we're done
	case COLL_SPHERE:       return dist_less_than   (pos, points[0], radius);
	case COLL_CYLINDER:     return dist_xy_less_than(pos, points[0], radius); // z has been checked above
	case COLL_CYLINDER_ROT: return sphere_intersect_cylinder(pos, 0.0, points[0], points[1], radius, radius2); // use a zero radius sphere
	case COLL_TORUS:        return sphere_torus_intersect(pos, 0.0, points[0], radius2, radius); // use a zero radius sphere
	case COLL_CAPSULE:      return sphere_int_capsule(pos, 0.0, points[0], points[1], radius, radius2); // use a zero radius sphere
	case COLL_POLYGON:
		if (is_thin_poly()) return 0; // thin polygons are effectively 2D and can't contain a point (not a volume)
		return sphere_ext_poly_intersect(points, npoints, norm, pos, 0.0, thickness, 0.0); // use a zero radius sphere
	}
	return 0; // never gets here
}

// Note: simplified version of coll_obj::intersects_cobj()
bool coll_obj::sphere_intersects(point const &sc, float sr) const {

	if (!sphere_cube_intersect(sc, sr, *this)) return 0; // optimization
	switch (type) {
	case COLL_CUBE:         return 1; // if bcube intersects the sphere, we're done
	case COLL_SPHERE:       return dist_less_than(sc, points[0], (radius+sr));
	case COLL_CYLINDER:     if (sphere_def_coll_vert_cylin(sc, sr, points[0], points[1], radius)) return 1; // fallthrough
	case COLL_CYLINDER_ROT: return sphere_intersect_cylinder(sc, sr, points[0], points[1], radius, radius2);
	case COLL_CAPSULE:      return sphere_int_capsule(sc, sr, points[0], points[1], radius, radius2);
	case COLL_TORUS:        return sphere_torus_intersect(sc, sr, points[0], radius2, radius);
	case COLL_POLYGON:      return sphere_ext_poly_intersect(points, npoints, norm, sc, sr, thickness, 0.0);
	}
	return 0; // never gets here
}


// ******************* OBJ_GROUP MEMBERS ******************


void obj_group::create(int obj_type_, unsigned max_objects_, unsigned init_objects_, unsigned app_rate_,
	bool init_enabled_, bool reorderable_, bool auto_max, bool predef_use_once_)
{
	flags        = JUST_INIT;
	type         = obj_type_;
	app_rate     = app_rate_;
	reorderable  = reorderable_;
	predef_use_once = predef_use_once_;
	init_objects = init_objects_;
	max_objs     = max_objects_;
	end_id       = 0;
	if (auto_max) flags |= APP_FROM_LT;

	if (type == PRECIP) {
		type   = RAIN;
		flags |= PRECIPITATION;
	}
	max_objs = get_updated_max_objs();
	assert(max_objs < 1000000); // sanity check
	set_enable(init_enabled_);
	init_group();
}


unsigned obj_group::get_updated_max_objs() const {

	if (!(flags & APP_FROM_LT)) return max(max_objs, init_objects);
	int const lifetime((flags & PRECIPITATION) ? max(object_types[RAIN].lifetime,
		max(object_types[HAIL].lifetime, object_types[SNOW].lifetime)) : object_types[type].lifetime);
	return unsigned(lifetime*app_rate + init_objects);
}


void obj_group::update_app_rate(float const val, unsigned min_app, unsigned max_app) {

	app_rate = max(min_app, unsigned(app_rate*val + 0.5));
	if (max_app > 0) {app_rate = min(app_rate, max_app);}
	unsigned const new_max_objs(get_updated_max_objs());

	if (enabled) { // change capacity (increase or decrease)
		assert(objects.size() == max_objs);
		objects.resize(new_max_objs);
		
		for (unsigned j = max_objs; j < new_max_objs; ++j) {
			objects[j] = def_objects[type]; // allocate new objects
			objects[j].status = 0;
		}
	}
	max_objs = new_max_objs;
}


void obj_group::init_group() {

	for (unsigned j = 0; j < max_objects(); ++j) {
		objects[j] = def_objects[type];
		if (j < init_objects) {gen_object_pos(objects[j].pos, object_types[type].flags);}
		else {objects[j].status = 0;}
	}
	flags &= (PRECIPITATION | APP_FROM_LT); // keep only this flag
}


// normally called before using objects, but can be called dynamically later
void obj_group::add_predef_obj(point const &pos, int type, int rtime) {
	
	predef_objs.emplace_back(pos, type, rtime);
	max_eq(max_objs, (unsigned)predef_objs.size());
	reorderable = 0; // need to unset reorderable so that predef_objs indexes remain correct
}


void obj_group::preproc_this_frame() {

	unsigned const nobjs((unsigned)max_objects());
	unsigned const max_used_id(min(nobjs, max(end_id, init_objects+1))); // one past the end
	end_id = nobjs; // likely will be reset to a smaller value below
	new_id = 0;
	if (!enabled) return;

	if (reorderable && begin_motion) { // some objects such as smileys are position dependent
		assert(predef_objs.empty());
		unsigned saw_id(0);
		for (unsigned j = 0; j < max_used_id; ++j) {if (objects[j].enabled()) {saw_id = j+1;}} // only need to check up to max_used_id
		sort(objects.begin(), (objects.begin() + saw_id));
		for (unsigned j = 0; j < nobjs; ++j) {if (!objects[j].enabled()) {end_id = j; break;}} // Note: will likely exit before j reaches max_used_id
	}
	if (!predef_use_once) {reap_predef_objs();}
}

void obj_group::reap_predef_objs(set<unsigned> const *const only_of_type) {

	for (vector<predef_obj>::iterator i = predef_objs.begin(); i != predef_objs.end(); ++i) {
		if (i->obj_used == -1) continue; // already unused/available
		assert((unsigned)i->obj_used < max_objects());
		if (only_of_type != nullptr && only_of_type->find(i->type) == only_of_type->end()) continue; // not a selected type
			
		if (objects[i->obj_used].disabled()) { // unused object
			i->obj_used = -1; // reset back to 'unused'
			i->cur_time = tfticks;
		}
	}
}


void obj_group::remove_reset_cobjs() {

	if (!large_radius()) return;

	for (unsigned j = 0; j < max_objects(); ++j) {
		if (objects[j].status == 0) continue; // is status check correct?
		if (objects[j].coll_id >= 0) {remove_reset_coll_obj(objects[j].coll_id);}
		if (type == MAT_SPHERE) {remove_mat_sphere(j);}
	}
}


void obj_group::create_object_at(unsigned i, point const &pos) {

	enable();
	assert(max_objects() == max_objs && i < max_objects());
	assert(type >= 0 && type < NUM_TOT_OBJS);
	assert(!is_nan(pos));
	if (objects[i].coll_id >= 0) {remove_reset_coll_obj(objects[i].coll_id);} // just in case
	objects[i]     = def_objects[type];
	objects[i].pos = pos;
}


int obj_group::choose_object(bool peek)  { // could return unsigned?
	
	enable();
	assert(max_objects() > 0); // enabled == 1 should be true after before the object is used
	if (!reorderable) return objects.choose_element(peek);
	assert(!peek);
	// Note: To guarantee correctness, the times of all objects should be updated by a constant amount per frame
	// Note: Assumes objects are sorted oldest to newest, with all unused objects at the end
	if (end_id <  max_objects()) return end_id++; // unused object
	if (new_id == max_objects()) new_id = 0; // wraparound (circular queue)
	return new_id++; // used, old object (increment so that the first object isn't reused in the same frame)
}


// return value: 0 = skip object, 1 = use predef object, 2 = gen new object
int obj_group::get_next_predef_obj(dwobject &obj, unsigned ix) {

	if (predef_objs.empty()) return 2;
	float min_time(0.0);
	int best_obj(-1);

	for (unsigned i = 0; i < predef_objs.size(); ++i) { // inefficient, but predef_objs should be small
		if (predef_objs[i].obj_used >= 0) continue; // already in use
		if (predef_objs[i].cur_time > 0.0 && (tfticks - predef_objs[i].cur_time) < predef_objs[i].regen_time) continue; // not regenerated
		if (best_obj >= 0 && predef_objs[i].cur_time >= min_time) continue; // not the oldest
		min_time = predef_objs[i].cur_time;
		best_obj = i;
	}
	if (best_obj == -1) return ((max_objects() > predef_objs.size()) ? 2 : 0); // no valid predef_obj found
	obj.pos       = predef_objs[best_obj].pos;
	obj.direction = (unsigned char)predef_objs[best_obj].type;
	obj.flags    |= USER_PLACED;
	predef_objs[best_obj].obj_used = ix;
	return 1;
}


bool obj_group::large_radius() const {

	assert(type >= 0 && type < NUM_TOT_OBJS);
	return (object_types[type].radius >= LARGE_OBJ_RAD);
}


void obj_group::enable() {

	if (enabled == 1) return;
	objects.resize(max_objs); // allocate
	unsigned const tflags(object_types[type].flags);

	if (tflags & (OBJ_ROLLS | VERTEX_DEFORM)) {
		if (!td) td.reset(new transform_data);
		if (tflags & OBJ_ROLLS)     {td->matrices.resize(max_objs);}
		if (tflags & VERTEX_DEFORM) {td->perturb_maps.resize(max_objs);}
	}
	end_id  = 0;
	enabled = 1;
}


void obj_group::disable() {

	if (enabled == 0) return;
	remove_reset_cobjs();
	objects.clear();
	if (!large_radius()) {remove_excess_cap(objects);} // free the memory
	if ((object_types[type].flags & (OBJ_ROLLS | VERTEX_DEFORM))) {td.reset();}
	for (vector<predef_obj>::iterator i = predef_objs.begin(); i != predef_objs.end(); ++i) {i->obj_used = -1;}
	end_id  = 0;
	enabled = 0;
}


void obj_group::free_objects() {

	remove_reset_cobjs();
	if (!objects.empty()) {reset_status(objects);}
	td.reset();
}


void obj_group::shift(vector3d const &vd) {

	if (!enabled) return;

	for (unsigned j = 0; j < max_objects(); ++j) {
		if (!objects[j].disabled()) {
			objects[j].pos += vd;
			if (type == SMILEY) shift_player_state(vd, j);
		}
	}
	for (unsigned i = 0; i < predef_objs.size(); ++i) {predef_objs[i].pos += vd;}
}


bool obj_group::obj_within_dist(unsigned i, point const &pos, float dist) const { // cache miss dominated

	assert(enabled && i < objects.size());
	if (objects[i].disabled()) return 0;
	if (fabs(objects[i].pos.x - pos.x) > dist || fabs(objects[i].pos.y - pos.y) > dist) return 0; // quick test
	return (p2p_dist_sq(objects[i].pos, pos) < dist*dist);
}

bool obj_group::temperature_ok() const {
	return ((flags & PRECIPITATION) || type == PRECIP || (temperature >= object_types[type].min_t && temperature < object_types[type].max_t));
}

bool obj_group::obj_has_shadow(unsigned obj_id) const {

	if (type != SMILEY) return 1;
	assert(obj_id < objects.size());
	return (!has_invisibility(obj_id));
}


// ******************* COBJ_GROUP_T MEMBERS ******************

void cobj_group_t::update_props() {

	volume = mass = 0.0;
	center_of_mass = all_zeros;

	for (auto i = begin(); i != end(); ++i) {
		coll_obj const &c(coll_objects.get_cobj(*i));
		volume += c.volume;
		mass   += c.get_mass();
		center_of_mass += c.get_mass()*c.get_center_of_mass(1); // weighted average (ignore_group=1)
	}
	center_of_mass /= mass; // normalize
	valid = 1;
}


// ******************* COBJ_DRAW_GROUPS MEMBERS ******************

unsigned cobj_draw_groups::new_group() {
	unsigned const group_id(groups.size());
	groups.push_back(cobj_draw_group());
	return group_id;
}

void cobj_draw_groups::add_to_group(coll_obj const &cobj) { // Note: cobj.dgroup_id should be set correctly
	get_group(cobj.dgroup_id).ids.push_back(dcobjs.size());
	dcobjs.push_back(cobj);
	dcobjs.back().setup_internal_state();
}

bool cobj_draw_groups::set_parent_or_add_cobj(coll_obj const &cobj) { // Note: cobj.dgroup_id should be set correctly

	assert(cobj.id >= 0);
	cobj_draw_group &group(get_group(cobj.dgroup_id));
	
	if (group.empty() && group.parent < 0) { // this is the parent cobj
		group.parent     = cobj.id;
		group.parent_pos = cobj.get_cube_center();
		return 1;
	}
	add_to_group(cobj); // not the parent
	return 0;
}

// Note: better to pass the parent in rather than using group.parent, which may be stale
// Note: non-const because it translates the objects (which is one reason why they can't be instanced)
vector<unsigned> const &cobj_draw_groups::get_draw_group(int group_id, coll_obj const &parent) {
	
	cobj_draw_group &group(get_group(group_id));
	point const ppos(parent.get_cube_center());
	
	if (ppos != group.parent_pos) { // translate the group
		vector3d const delta(ppos - group.parent_pos); // determine parent translation
		group.parent_pos = ppos; // update parent pos with new value

		for (auto i = group.ids.begin(); i != group.ids.end(); ++i) {
			assert(*i < dcobjs.size());
			dcobjs[*i].shift_by(delta);
		}
	}
	return group.ids;
}


// ******************* OBJ_DRAW_GROUP MEMBERS ******************

void obj_draw_group::free_vbo() {

	assert(!inside_beg_end);
	if (!use_vbo) return;
	if (vbo) {delete_and_zero_vbo(vbo);}
	verts.clear();
	start_cix = end_cix = num_verts = 0;
}


void obj_draw_group::draw_vbo() const {

	if (num_verts == 0) return;
	check_bind_vbo(vbo);
	draw_verts<vert_norm>(NULL, num_verts, GL_TRIANGLES);
	bind_vbo(0);
}


bool obj_draw_group::begin_render(unsigned &cix) { // returns true if rendering can be skipped

	assert(!inside_beg_end); // can't call begin_render() twice
	inside_beg_end = 1;
	assert(verts.empty());
	if (!use_vbo) return 0;

	if (!vbo) { // no vbo exists, so create it
		start_cix = cix;
		return 0;
	}
	// already created the vbo, so draw it
	assert(start_cix == cix); // must start at the same location as last time
	assert(start_cix < end_cix); // can't have an empty range
	cix = end_cix-1; // move the current index to the one before the end
	draw_vbo();
	return 1;
}


void obj_draw_group::end_render() {

	assert(inside_beg_end); // can't call end_render() without begin_render()
	inside_beg_end = 0;
	if (verts.empty()) return;

	if (use_vbo) {
		assert(!vbo);
		create_vbo_and_upload(vbo, verts);
		num_verts = verts.size();
		draw_vbo();
	}
	else {
		draw_verts(verts, GL_TRIANGLES);
	}
	verts.resize(0);
}


void obj_draw_group::add_draw_polygon(point const *const points, vector3d const &normal, unsigned npoints, unsigned cix) {

	assert(inside_beg_end);
	assert(npoints == 3 || npoints == 4);
	//normal = get_norm_camera_orient(normal, points[0]);
		
	for (unsigned i = 0; i < ((npoints == 3) ? 3U : 6U); ++i) { // triangle or quad (2 tris)
		verts.emplace_back(points[quad_to_tris_ixs[i]], normal);
	}
	if (use_vbo) {end_cix = cix+1;}
}


void free_cobj_draw_group_vbos() {
	for (vector<obj_draw_group>::iterator i = obj_draw_groups.begin(); i != obj_draw_groups.end(); ++i) {i->free_vbo();}
}

