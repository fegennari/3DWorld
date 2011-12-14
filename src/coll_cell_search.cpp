// 3D World - OpenGL CS184 Computer Graphics Project - collision detection coll cell iteration and search
// by Frank Gennari
// 2/6/06

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"


bool const CACHE_COBJ_LITES   = 0;
unsigned const QLP_CACHE_SIZE = 10000;


int cobj_counter(0);

extern int display_mode, begin_motion;
extern float zmin, zbottom, water_plane_z;
extern vector<coll_obj> coll_objects;


// returns 1 if there is no intersection
bool coll_obj::cobj_plane_side_test(point const *pts, unsigned npts, point const &lpos) const {

	vector<tquad_t> ppts;

	for (unsigned i = 0; i < npts; ++i) {
		point const spts[3] = {pts[i], pts[(i+1)%npts], lpos};
		point const center(get_center(spts, 3));
		vector3d const pts_norm(get_poly_norm(spts));
		point pt;

		if (type == COLL_POLYGON) {
			if (thickness > MIN_POLY_THICK) { // thick polygon
				if (ppts.empty()) ppts = thick_poly_to_sides(points, npoints, norm, thickness);

				for (unsigned i = 0; i < ppts.size(); ++i) {
					for (unsigned j = 0; j < ppts[i].npts; ++j) {
						if (dot_product_ptv(pts_norm, center, ppts[i][j]) > 0.0) return 0;
					}
				}
			}
			else { // thin polygon
				for (unsigned i = 0; i < unsigned(npoints); ++i) {
					if (dot_product_ptv(pts_norm, center, points[i]) > 0.0) return 0;
				}
			}
		}
		else { // use bounding cube (bad for pine tree cones)
			for (unsigned x = 0; x < 2; ++x) {
				pt[0] = d[0][x];
				for (unsigned y = 0; y < 2; ++y) {
					pt[1] = d[1][y];
					for (unsigned z = 0; z < 2; ++z) {
						pt[2] = d[2][z];
						if (dot_product_ptv(pts_norm, center, pt) > 0.0) return 0;
					}
				}
			}
		}
	}
	return 1;
}


// false negatives are OK except when called from check_coll_line()
bool coll_obj::line_intersect(point const &p1, point const &p2) const {

	if (!check_line_clip(p1, p2, d)) return 0;
	
	switch (type) {
		case COLL_CUBE:
			return 1;
		case COLL_SPHERE:
			return line_sphere_intersect(p1, p2, points[0], radius);
		case COLL_CYLINDER:
		case COLL_CYLINDER_ROT:
			return line_intersect_cylinder(p1, p2, cylinder_3dw(points[0], points[1], radius, radius2), !(cp.surfs & 1));
		case COLL_POLYGON: // must be coplanar
			assert(npoints >= 3);

			if (thickness > MIN_POLY_THICK) { // test extruded (3D) polygon
				point pts[2][4];
				gen_poly_planes(points, npoints, norm, thickness, pts);
				vector3d const v1(p2, p1);
				bool const test_side(dot_product(v1, norm) > 0.0);
				if (thick_poly_intersect(v1, p1, norm, pts, test_side, npoints)) return 1;
			}
			else { // test planar (2D) polygon
				float t;
				if (!line_poly_intersect(p1, p2, points, npoints, norm, t)) return 0;
				return check_poly_billboard_alpha(p1, p2, t);
			}
			break;
	}
	return 0;
}


bool coll_obj::line_int_exact(point const &p1, point const &p2, float &t, vector3d &cnorm, float tmin, float tmax) const {

	float clip_tmin, clip_tmax;
	if (type != COLL_POLYGON && (!get_line_clip(p1, p2, d, clip_tmin, clip_tmax)
		                         || clip_tmin > tmax || clip_tmax < tmin)) return 0;
	
	switch (type) {
		case COLL_CUBE:
			{
				t = clip_tmin;
				if (t > tmax || t < tmin) return 0;
				get_closest_cube_norm(d, (p1 + (p2 - p1)*t), cnorm);
				return 1;
			}
		case COLL_SPHERE:
			{
				point coll_pos;
				vector3d const v1((p2 - p1).get_norm());
				if (!line_sphere_int(v1, p1, points[0], radius, coll_pos, 0)) return 0;
				t = -1.0; // start at a bad value
				
				for (unsigned i = 0; i < 3; ++i) {
					if (fabs(p2[i] - p1[i]) > TOLERANCE) {
						t = (coll_pos[i] - p1[i])/(p2[i] - p1[i]);
						break;
					}
				}
				if (t > tmax || t < tmin) return 0;
				cnorm = (coll_pos - points[0]);
				if (!cnorm.normalize_test()) cnorm = plus_z; // arbitrary
				return 1;
			}
		case COLL_CYLINDER:
		case COLL_CYLINDER_ROT:
			{
				int const int_type(line_int_thick_cylinder(p1, p2, points[0], points[1], 0.0, 0.0, radius, radius2, 1, t));
				if (!int_type || t > tmax || t < tmin) return 0;
				
				if (int_type == 1) { // side intersection
					vector3d const cv(points[0] - points[1]);
					orthogonalize_dir((p1 - points[0]), cv, cnorm, 0);
					if (cnorm == zero_vector) orthogonalize_dir((p2 - points[0]), cv, cnorm, 0); // p1 is bad, so try p2

					if (radius != radius2) {
						if (!cnorm.normalize_test()) cnorm = plus_z; // arbitrary
						float const len(cv.mag());

						if (len > TOLERANCE) {
							float const dr(radius2 - radius), denom(sqrt(len*len + dr*dr));
							assert(denom > TOLERANCE);
							cnorm *= len/denom;
							cnorm += cv*(dr*len/denom);
						}
					}
				}
				else { // top/bottom intersection
					cnorm = (points[1] - points[0]);
					if ((int_type == 2) ^ (radius >= radius2)) cnorm.negate(); // r1 >= r2 => swap
				}
				if (!cnorm.normalize_test()) cnorm = plus_z; // arbitrary
				return 1;
			}
		case COLL_POLYGON: // must be coplanar
			{
				assert(npoints >= 3);

				if (thickness > MIN_POLY_THICK) { // test extruded (3D) polygon
					t = 2.0; // start at a bad value
					float tval;
					point pts[2][4];
					gen_poly_planes(points, npoints, norm, thickness, pts);
					bool const test_side(dot_product((p2 - p1), norm) > 0.0);
					point const *const points2(pts[test_side]);
					
					if (line_poly_intersect(p1, p2, points2, npoints, norm, tval) && (tval <= tmax && tval >= tmin)) {
						t     = tval;
						cnorm = get_poly_dir_norm(norm, p1, (p2 - p1), t);
					}
					for (int j = 0; j < npoints; ++j) { // now test the <npoints> sides
						unsigned const jnext((j+1)%npoints);
						point const side_pts[4] = {pts[0][j], pts[0][jnext], pts[1][jnext], pts[1][j]};
						vector3d const side_norm(get_poly_norm(side_pts));
					
						if (line_poly_intersect(p1, p2, side_pts, 4, side_norm, tval)) {
							if (tval < t && (tval <= tmax && tval >= tmin)) {
								t     = tval;
								cnorm = get_poly_dir_norm(side_norm, p1, (p2 - p1), t);
							}
						}
					}
					return (t <= tmax && t >= tmin);
				}
				if (!line_poly_intersect(p1, p2, points, npoints, norm, t) || t > tmax || t < tmin) return 0;
				if (!check_poly_billboard_alpha(p1, p2, t)) return 0;
				cnorm = get_poly_dir_norm(norm, p1, (p2 - p1), t);
				return 1;
			}
	}
	return 0;
}


class water_spalsh_search {

	point const &pos1, &pos2;
	float splash_val;

public:
	water_spalsh_search(point const &pos1_, point const &pos2_, float splash_val_) :
	   pos1(pos1_), pos2(pos2_), splash_val(splash_val_) {}

	bool do_iter() const {
		if (splash_val == 0.0) return 0;
		int const xa(get_xpos(pos1.x)), ya(get_ypos(pos1.y)), xb(get_xpos(pos2.x)), yb(get_ypos(pos2.y));
		int const dx(xb - xa), dy(yb - ya), steps(max(1, ((abs(dx) > abs(dy)) ? abs(dx): abs(dy))));
		double const xinc(dx/(double)steps), yinc(dy/(double)steps);
		double x(xa), y(ya);
		int bnds[2][2] = {{MESH_X_SIZE,0}, {MESH_Y_SIZE,0}}; // {x,y}{min,max}

		for (int k = 0; k <= steps; ++k) { // DDA algorithm
			int const xpos(int(x + 0.5)), ypos(int(y + 0.5));
			x += xinc;
			y += yinc;
			if (point_outside_mesh(xpos, ypos) || !mesh_is_underwater(xpos, ypos)) continue;
			float const wmz(water_matrix[ypos][xpos]);
			if ((pos1.z < wmz) ^ (pos2.z > wmz)) continue;
			float const tz((wmz - pos1.z)/(pos2.z - pos1.z));
			point const splash_pos(pos1 + (pos2 - pos1)*tz);
			float const mx(get_xval(xpos)), my(get_yval(ypos));
		
			if (splash_pos.x > (mx-DX_VAL) && splash_pos.x < (mx+DX_VAL) && splash_pos.y > (my-DY_VAL) && splash_pos.y < (my+DY_VAL)) {
				add_splash(xpos, ypos, 25.0, 0.01, 1); // dynamic water
				draw_splash(splash_pos.x, splash_pos.y, (wmz + 0.0001), splash_val);
				gen_line_of_bubbles(splash_pos, pos2);
				return 1;
			}
		} // for k
		return 0;
	}
};


bool check_coll_line(point pos1, point pos2, int &cindex, int cobj, int skip_dynamic, int test_alpha) {

	if (check_coll_line_tree(pos1, pos2, cindex, cobj, 0, test_alpha, (skip_dynamic >= 2))) return 1;

	if (!skip_dynamic && begin_motion) { // find dynamic cobj intersection
		update_cobj_tree(1, 0);
		if (check_coll_line_tree(pos1, pos2, cindex, cobj, 1, test_alpha)) return 1;
	}
	return 0;
}


bool check_coll_line_exact(point pos1, point pos2, point &cpos, vector3d &cnorm, int &cindex, float splash_val,
						   int ignore_cobj, bool fast, bool test_alpha, bool skip_dynamic)
{
	if (check_coll_line_exact_tree(pos1, pos2, cpos, cnorm, cindex, ignore_cobj, 0, test_alpha)) pos2 = cpos;

	if (!skip_dynamic && begin_motion) { // find dynamic cobj intersection
		update_cobj_tree(1, 0);
		int cindex2;
		if (check_coll_line_exact_tree(pos1, pos2, cpos, cnorm, cindex2, ignore_cobj, 1, test_alpha)) cindex = cindex2;
	}
	if (splash_val > 0.0) { // handle water splashes
		if (cindex >= 0) pos2 = cpos;
		
		if (do_line_clip_scene(pos1, pos2, zbottom, max(ztop, water_plane_z))) { // max of dynamic and static water
			water_spalsh_search wss(pos1, pos2, splash_val);
			wss.do_iter();
		}
	}
	return (cindex >= 0);
}


bool cobj_contained(point pos1, point center, const point *pts, unsigned npts, int cobj) {

	if (!have_occluders()) return 0;
	assert(npts > 0);
	static int last_cobj(-1);

	if (last_cobj >= 0 && last_cobj != cobj && !coll_objects[last_cobj].disabled()) {
		if (coll_objects[last_cobj].intersects_all_pts(pos1, pts, npts)) return 1;
	}
	return cobj_contained_tree(pos1, center, pos1, pts, npts, cobj, last_cobj);
}


bool get_coll_line_cobjs(point pos1, point pos2, int cobj, vector<int> &cobjs) {

	cobjs.resize(0);
	get_coll_line_cobjs_tree(pos1, pos2, cobj, cobjs, 0, 1);
	return (!cobjs.empty());
}


bool coll_obj::is_occluder() const {
	
	if (status != COLL_STATIC || !cp.draw || is_semi_trans()) return 0; // cp.might_be_drawn()?
	if (type == COLL_CUBE)    return 1;
	if (type != COLL_POLYGON) return 0;
	unsigned big_dims(0);
	UNROLL_3X(if ((d[i_][1] - d[i_][0]) > 0.2*SCENE_SIZE[i_]) ++big_dims;)
	return (big_dims >= 2);
}


bool coll_obj::intersects_all_pts(point const &pos, point const *const pts, unsigned npts) const {

	switch (type) {
	case COLL_CUBE:
		for (unsigned i = 0; i < npts; ++i) { // can almost skip two corners on a quad
			if (!check_line_clip(pos, pts[i], d)) return 0;
		}
		break;
	case COLL_POLYGON:
		for (unsigned i = 0; i < npts; ++i) {
			float t; // unused
			if (!line_poly_intersect(pos, pts[i], points, npoints, norm, t)) return 0;
		}
		break;
	default:
		return 0; // not supported
	}
	return 1;
}


bool is_occluded(vector<int> const &occluders, point const *const pts, int npts, point const &camera) {

	unsigned const nocc(occluders.size());

	for (unsigned i = 0; i < nocc; ++i) { // cache last occluder?, promote to the front if occluded?
		coll_obj const &cobj(coll_objects[occluders[i]]);
		if (cobj.status == COLL_STATIC && cobj.intersects_all_pts(camera, pts, npts)) return 1;
	}
	return 0;
}


void get_occluders() {

	RESET_TIME;
	if (!(display_mode & 0x08) || !have_occluders()) return;
	static unsigned startval(0), stopped_count(0);
	static bool first_run(1);
	unsigned const skipval(first_run ? 0 : 8); // spread update across many frames
	first_run = 0;
	if (++startval >= skipval) startval = 0;
	static point last_camera(FAR_CLIP, FAR_CLIP, FAR_CLIP);
	point const camera(get_camera_pos());
	
	if (p2p_dist(camera, last_camera) < 0.1*HALF_DXY) { // camera hasn't moved much
		if (skipval == 0 || stopped_count == skipval) return;
		++stopped_count;
	}
	else {
		stopped_count = 0;
		last_camera   = camera;
	}
	unsigned const ncobjs(coll_objects.size());

	//#pragma omp parallel for schedule(static,1) // helps in some cases and hurts in others
	for (unsigned i = startval; i < ncobjs; i += max(1U, skipval)) {
		coll_obj &cobj(coll_objects[i]);
		if (!cobj.fixed || cobj.group_id >= 0 || cobj.status != COLL_STATIC || !cobj.cp.draw || cobj.cp.surfs == EF_ALL) continue;
		get_coll_line_cobjs(camera, cobj.get_cube_center(), i, coll_objects[i].occluders);
	}
	if (skipval == 0) {PRINT_TIME("Occlusion Preprocessing");}
}



