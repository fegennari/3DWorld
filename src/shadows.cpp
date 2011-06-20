// 3D World - OpenGL CS184 Computer Graphics Project - Mesh Shadow Generation code
// by Frank Gennari
// 3/6/06

#include "3DWorld.h"
#include "mesh.h"
#include "tree_3dw.h"
#include "physics_objects.h"


int const SHADOW_BORDER     = 1;
int const CLIP_SHADOW_CUBES = 1;


struct shad_env { // shadow envelope

	int enabled, x1, y1, x2, y2;
	float z1, z2;

	shad_env() : enabled(0) {}
};


int dshadow_lights(0);
point light_pos;
shad_env s_env[NUM_LIGHT_SRC];

extern bool combined_gu, use_stencil_shadows;
extern int stencil_shadow, island, ground_effects_level;
extern float sun_rot, moon_rot, zmin, zmax, zbottom, ztop;
extern point sun_pos, moon_pos, mesh_origin;
extern vector3d up_norm;
extern lightning l_strike;
extern tree_cont_t t_trees;
extern vector<coll_obj> coll_objects;


int check_shadow_edge_clip(point const &pt, point const &lpos, int &xmin, int &xmax, int &ymin, int &ymax);



void create_shadows() {

	glClear(GL_STENCIL_BUFFER_BIT);
	glDisable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_STENCIL_TEST);

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!(stencil_shadow & (1 << l))) continue;
		if ((l == LIGHT_SUN && light_factor < 0.5) || (l == LIGHT_MOON && light_factor >= 0.6)) continue;
		float shadow_alpha;

		if (light_factor <= 0.4 || light_factor >= 0.6) {
			shadow_alpha = SHADOW_COLOR;
		}
		else {
			if (l == LIGHT_SUN) {
				shadow_alpha = 10.0*(light_factor - 0.5)*SHADOW_COLOR;
			}
			else if (l == LIGHT_MOON) {
				shadow_alpha = 5.0*(0.6 - light_factor)*SHADOW_COLOR;
			}
			else continue;
		}
		int inverts(0); // two pass stencil
		glColorMask(0, 0, 0, 0);
		glStencilFunc(GL_ALWAYS, 1, ~0);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK); // could use glFrontFace(GL_CW)
		glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);
		inverts = draw_shadowed_objects(l); // draw_world.cpp
		glCullFace(GL_FRONT);
		glStencilOp(GL_KEEP, GL_KEEP, GL_DECR);
		draw_shadowed_objects(l); // draw_world.cpp
		glCullFace(GL_BACK);
		glDisable(GL_CULL_FACE);
		glColorMask(1, 1, 1, 1);

		// draw a shadowing rectangle covering the entire screen
		glColor4f(0.0, 0.0, 0.0, shadow_alpha);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glStencilFunc(((inverts & 1) ? GL_EQUAL : GL_NOTEQUAL), 0, ~0); // test for invert stencil
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		glPushMatrix();
		glLoadIdentity();
		vector<camera_filter> cfs;
		cfs.push_back(camera_filter(colorRGBA(0.0, 0.0, 0.0, shadow_alpha), 1, -1));
		draw_camera_filters(cfs);
		glPopMatrix();
	}
	glDepthFunc(GL_LEQUAL);
	glDepthMask(GL_TRUE);
	glEnable(GL_LIGHTING);
	glDisable(GL_STENCIL_TEST);
	glShadeModel(GL_SMOOTH);
}


void update_sun_shadows() {

	if (!combined_gu && fabs(sun_rot/PI - 1.0) < 0.6 && fabs(sun_rot/PI - 1.0) > 0.4) {
		calc_visibility(MOON_SHADOW); // seems like this has to be first, but questionable
	}
	calc_visibility(SUN_SHADOW);
}


void update_sun_and_moon() {

	float const radius(0.6*(FAR_CLIP+X_SCENE_SIZE));
	sun_rot      = fix_angle(sun_rot);
	moon_rot     = fix_angle(moon_rot);
	light_factor = fabs(sun_rot/PI - 1.0);
	moon_pos     = mesh_origin + rtp_to_xyz(radius, MOON_THETA, moon_rot);
	sun_pos      = mesh_origin + rtp_to_xyz(radius,  SUN_THETA, sun_rot);
	light_pos    = ((light_factor >= 0.5 || combined_gu) ? sun_pos : moon_pos);
	up_norm      = light_pos.get_norm();
}


int light_valid(char light_sources, int l, point &lpos) {

	if (!(light_sources & (1 << l))) return 0;
	if ((l == LIGHT_SUN && light_factor < 0.4) || (l == LIGHT_MOON && light_factor > 0.6)) return 0;
	if (!get_light_pos(lpos, l) || lpos.z < zmin) return 0;
	return 1;
}


void coll_obj::add_shadow(char light_sources, bool dynamic) const {

	if (status != COLL_STATIC || cp.color.alpha == 0.0 || !cp.shadow) return;
	if (dynamic != dynamic_shadows_only()) return;

	switch (type) {
	case COLL_CUBE:
		cube_shadow(*this, light_sources, dynamic, 1);
		break;
	case COLL_SPHERE:
		sphere_shadow(points[0], radius, light_sources, dynamic, 1);
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		cylinder_shadow(points[0], points[1], radius, radius2, light_sources, !(cp.surfs & 1), dynamic, 1);
		break;
	case COLL_POLYGON:
		assert(npoints > 2);
		polygon_shadow(points, norm, npoints, thickness, light_sources, dynamic, 1, 0, (has_poly_billboard_alpha() ? cp.tid : -1));
		break;
	default: assert(0);
	}
}


void add_cobj_shadows(char light_sources) {

	if (ground_effects_level == 0) return; // disabled
	//RESET_TIME;
#if 0
	// three problems with this approach, though this is *much* simpler
	// 1. tree shadows are only added if tree cobjs are created
	// 2. requires cobj_tree to be efficient
	// 3. the test line wil not hit very small cobjs, so shadows will be missed
	point lpos;

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid(light_sources, l, lpos)) continue;
		
		for (int y = 0; y < MESH_Y_SIZE; ++y) {
			for (int x = 0; x < MESH_X_SIZE; ++x) {
				if (shadow_mask[l][y][x] & OBJECT_SHADOW) continue;
				point const pos(get_xval(x), get_yval(y), mesh_height[y][x]);
				int cindex; // unused
				if (!coll_pt_vis_test(pos, lpos, 0.0, cindex, -1, 1, 0)) shadow_mask[l][y][x] |= OBJECT_SHADOW;
			}
		}
	}
#else
	for (unsigned i = 0; i < t_trees.size(); ++i) {
		t_trees[i].gen_tree_shadows(light_sources, i);
	}
	for (unsigned i = 0; i < coll_objects.size(); ++i) {
		coll_objects[i].add_shadow(light_sources, 0);
	}
#endif
	//PRINT_TIME("Cobj Shadows");
}


// first point must be the light position
bool point_to_point_visibility(point const &pos1, point const &pos2, int &xpos, int &ypos, float &zval, int fast, int light) {

	// extend line to span z values of [zmin+?, zmax-?]
	// P = pos1 + t*(pos2 - pos1) = pos1 + t*v12
	// t = (P - pos1)/v12
	vector3d const v12(pos2, pos1);
	float z1(zmin), z2(zmax);

	if (light >= 0) {
		assert(light < NUM_LIGHT_SRC);

		if (s_env[light].enabled) {
			z1 = max(z1, s_env[light].z1);
			z2 = min(z2, s_env[light].z2);
		}
	}
	if (v12.z > 0.0) z1 = max(z1, pos2.z);
	if (v12.z < 0.0) z2 = min(z2, pos2.z);
	float const tmin((z1 - pos1.z)/v12.z), tmax((z2 - pos1.z)/v12.z);
	point pt1((pos1.x + tmin*v12.x), (pos1.y + tmin*v12.y), z1);
	point pt2((pos1.x + tmax*v12.x), (pos1.y + tmax*v12.y), z2);
	if (pos1.z > pos2.z) swap(pt1, pt2);
	return line_intersect_mesh(pt1, pt2, xpos, ypos, zval, fast, 1);
}


int camera_shadow(point const &camera) {

	return (has_invisibility(CAMERA_ID) ? 0 : sphere_shadow2(camera, CAMERA_RADIUS, CHECK_ALL_SHADOW, 1, 1));
}


int sphere_shadow2(point const &pos, float radius, char light_sources, int is_dynamic, int quality) {

	if (!is_dynamic || !use_stencil_shadows) {
		sphere_shadow(pos, radius, light_sources, is_dynamic, quality);
		return 0;
	}
	stencil_shadow = 0xFFFF;
	return stencil_shadow;
}


#define ADD_POINT_BB(xp, yp) \
	xmin = min(xmin, xp); \
	xmax = max(xmax, xp); \
	ymin = min(ymin, yp); \
	ymax = max(ymax, yp);


int get_shape_shadow_bb(point const *points, int npoints, int l, int quality, point const &lpos,
						int &xmin, int &ymin, int &xmax, int &ymax, int &ret_val, unsigned char stype)
{
	assert(points != NULL);
	int xp, yp, miss(0), ss2(0), saw(0), tot_pts(0);
	float zval, radius;
	point points2[1];
	xmin = MESH_X_SIZE - SHADOW_BORDER - 1;
	ymin = MESH_Y_SIZE - SHADOW_BORDER - 1;
	xmax = SHADOW_BORDER;
	ymax = SHADOW_BORDER;

	if (quality == 0) {
		polygon_bounding_sphere(points, npoints, 0.0, points2[0], radius);

		if (radius < HALF_DXY) { // small polygon - use a single point at its center
			npoints = 1;
			points  = points2;
		}
	}
	for (int i = 0; i < npoints; ++i) { // trace an outline with ray casting
		point const &p_cur(points[i]), &p_next(points[(i+1)%npoints]);
		float const len_thresh(16*HALF_DXY);
		vector3d delta(zero_vector);
		int nsteps(1);
		point pt(p_cur);

		if (!dist_less_than(p_cur, p_next, len_thresh)) {
			nsteps = int(p2p_dist(p_cur, p_next)/len_thresh) + 1;
			delta  = (p_next - p_cur)/float(nsteps);
		}
		for (int j = 0; j < nsteps; ++j) {
			if (!is_above_mesh(pt)) {
				xp = get_xpos(pt.x);
				yp = get_ypos(pt.y);

				if (!point_outside_mesh(xp, yp)) {
					ADD_POINT_BB(xp, yp);
					saw = 1;
				}
			}
			else if (point_to_point_visibility(lpos, pt, xp, yp, zval, 1, l)) {
				// pt.z = object z, lpos.z = light z, zval = mesh intersection z
				if ((pt.z > lpos.z && pt.z      < zval)    // mesh   => object => light    - can this happen?
					|| (pt.z < lpos.z && pt.z   > zval)    // light  => object => mesh
					|| (pt.z > zval   && lpos.z < zval)    // object => mesh   => light  ? - can this happen?
					|| (pt.z < zval   && lpos.z > zval)) { // light  => mesh   => object ?
					ADD_POINT_BB(xp, yp);
					saw = 1;
				}
			}
			else {
				ss2 |= check_shadow_edge_clip(pt, lpos, xmin, xmax, ymin, ymax);
				++miss;
			}
			++tot_pts;
			pt += delta;
		} // for j
	} // for i
	if (miss >= (tot_pts+1)/2) ret_val |= (1 << l);
	if (!saw && (!ss2 || xmin > xmax || ymin > ymax)) return 0;
	int const sborder((quality || npoints == 1) ? SHADOW_BORDER : 0);
	xmin = max(xmin - sborder, 0);
	ymin = max(ymin - sborder, 0);
	xmax = min(xmax + sborder, MESH_X_SIZE-1);
	ymax = min(ymax + sborder, MESH_Y_SIZE-1);
	if (xmin == xmax && ymin == ymax && !point_outside_mesh(xmin, ymin) && (shadow_mask[l][ymin][xmin] & stype)) return 0; // already shadowed
	return 1;
}


void get_sphere_points(point const &pos, float radius, point pts[8]) {

	float const radsq2(radius/SQRT2);
	for (int i = 0; i < 8; ++i) pts[i] = pos;
	pts[0].y += radius; pts[1].y -= radius;
	pts[2].x += radius; pts[3].x -= radius;
	pts[4].x += radsq2; pts[4].y += radsq2;
	pts[5].x += radsq2; pts[5].y -= radsq2;
	pts[6].x -= radsq2; pts[6].y += radsq2;
	pts[7].x -= radsq2; pts[7].y -= radsq2;
}


int enable_shadow_envelope(point const &pos, float radius, char light_sources, int is_dynamic) {

	unsigned char const SHADOW_TYPE(is_dynamic ? DYNAMIC_SHADOW : OBJECT_SHADOW);
	int shadowed(0), ret_val(0);
	point lpos, pts[8];
	get_sphere_points(pos, radius, pts);

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid(light_sources, l, lpos)) continue;
		//if (!is_visible_from_light(pos, lpos, 1)) continue;
		shad_env &se(s_env[l]);
		se.enabled = get_shape_shadow_bb(pts, 8, l, 1, lpos, se.x1, se.y1, se.x2, se.y2, ret_val, SHADOW_TYPE);

		if (se.enabled) {
			se.z1 = zmax;
			se.z2 = zmin;
			
			for (int i = se.y1; i <= se.y2; ++i) {
				for (int j = se.x1; j <= se.x2; ++j) {
					if (shadow_mask[l][i][j] & SHADOW_TYPE) continue;
					se.z1 = min(se.z1, mesh_height[i][j]);
					se.z2 = max(se.z2, mesh_height[i][j]);
				}
			}
			if (se.z1 > se.z2) { // no unshadowed vertices in range
				se.enabled = 0;
			}
			else {
				shadowed = 1;
			}
		}
	}
	return shadowed;
}


void disable_shadow_envelope(char light_sources) {

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (light_sources & (1 << l)) s_env[l].enabled = 0;
	}
}


int sphere_shadow(point const &pos, float radius, char light_sources, int is_dynamic, int quality) {

	int xmin, xmax, ymin, ymax, ret_val(0);
	float const rad2(radius*radius);
	unsigned char const SHADOW_TYPE(is_dynamic ? DYNAMIC_SHADOW : OBJECT_SHADOW);
	point pts[8], lpos;
	get_sphere_points(pos, radius, pts);

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid(light_sources, l, lpos)) continue;
		if (!is_visible_from_light(pos, lpos, 1)) continue;
		if (!get_shape_shadow_bb(pts, 8, l, quality, lpos, xmin, ymin, xmax, ymax, ret_val, SHADOW_TYPE)) continue;
		if (is_dynamic) dshadow_lights |= (1 << l);

		for (int i = ymin; i <= ymax; ++i) { // fast approximation - uses hybrid ray casting and intersection test area fill
			for (int j = xmin; j <= xmax; ++j) {
				if (shadow_mask[l][i][j] & SHADOW_TYPE) continue;
				point const pt(get_xval(j), get_yval(i), mesh_height[i][j]);
				if (sphere_test_comp(pt, pos, vector3d(pt, lpos), rad2)) shadow_mask[l][i][j] |= SHADOW_TYPE;
			}
		}
	}
	return ret_val;
}


// used for tree branches, etc.
int cylinder_shadow(point p1, point p2, float radius1, float radius2, char light_sources, int shadow_ends, int is_dynamic, int quality) {

	if (radius1 == 0.0 && radius2 == 0.0) return 0;
	int xmin, xmax, ymin, ymax, ret_val(0);
	float t; // unused
	unsigned char const SHADOW_TYPE(is_dynamic ? DYNAMIC_SHADOW : OBJECT_SHADOW);
	point pts[8], lpos;
	vector3d v2, norm;
	float const len(p2p_dist(p1, p2)), rfudge(0.2*HALF_DXY);
	if (radius1 > 0.0) radius1 += rfudge;
	if (radius2 > 0.0) radius2 += rfudge;

	if (p1.x == p2.x && p1.y == p2.y) { // vertical cylinder
		if (p1.z > p2.z) swap(p1.z, p2.z);
		int const xp(get_xpos(p1.x)), yp(get_ypos(p1.y));

		if (point_interior_to_mesh(xp, yp)) {
			float const mh(mesh_height[yp][xp]);
			if (p2.z < mh) return 0; // entire cylinder below mesh
			if (p1.z < mh) p1.z = mh + 0.05*(p2.z - p1.z); // move slightly above mesh
		}
	}
	cylinder_3dw const cylin(p1, p2, radius1, radius2);
	point const p[2]   = {p1, p2};
	float const r[2]   = {radius1, radius2};
	float const rsq[2] = {radius1*radius1, radius2*radius2};

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid(light_sources, l, lpos)) continue;

		if (max(len, max(radius1, radius2)) > HALF_DXY) {
			if (!is_visible_from_light(p1, lpos, 1) && !is_visible_from_light(p2, lpos, 1)) continue;
		}
		int npts(0);
		vector3d const v1(p1, lpos);
		cylinder_quad_projection(pts, cylin, v1, v2, npts);
		int const ncpts(npts);

		if (shadow_ends) { // usually false
			assert(radius1 > 0.0 || radius2 > 0.0);
			vector3d ortho2;
			orthogonalize_dir(v1, v2, ortho2, 1);
			ortho2.negate();
			for (unsigned d = 0; d < 2; ++d) gen_cylin_pts(pts, npts, p[d], r[d], ortho2);
		}
		if (!get_shape_shadow_bb(pts, npts, l, quality, lpos, xmin, ymin, xmax, ymax, ret_val, SHADOW_TYPE)) continue;
		if (is_dynamic) dshadow_lights |= (1 << l);
		get_normal(pts[0], pts[1], pts[2], norm, 0);

		for (int i = ymin; i <= ymax; ++i) { // fast approximation - uses hybrid ray casting and intersection test area fill
			for (int j = xmin; j <= xmax; ++j) {
				if (shadow_mask[l][i][j] & SHADOW_TYPE) continue;
				point const pt(get_xval(j), get_yval(i), mesh_height[i][j]);
				vector3d const v1(lpos, pt);
				bool shadowed(line_poly_intersect(v1, pt, pts, ncpts, norm));
				
				if (shadow_ends) { // use circle for shadow
					for (unsigned d = 0; d < 2 && !shadowed; ++d) {
						if (rsq[d] > 0.0 && circle_test_comp(pt, p[d], v1, v2, rsq[d], t)) shadowed = 1; // end pt
					}
				}
				if (shadowed) shadow_mask[l][i][j] |= SHADOW_TYPE;
			}
		}
	}
	return ret_val;
}


// used for cubes, leaves, etc.
int polygon_shadow(point const *points, vector3d const &norm, int npoints, float thick, char light_sources,
				   int is_dynamic, int quality, int is_cube, int tid)
{
	assert(points && npoints >= 3);
	if (npoints < 3) return 0;
	int xmin, xmax, ymin, ymax, ret_val(0);
	unsigned char const SHADOW_TYPE(is_dynamic ? DYNAMIC_SHADOW : OBJECT_SHADOW);
	bool test_side(0);
	point lpos;
	static vector<point> pts[2];
	bool const thick_test(thick > MIN_POLY_THICK2);
	point const center(get_center(points, npoints));
	if (thick_test) gen_poly_planes(points, npoints, norm, thick, pts); // generate top and bottom surfaces
	
	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid(light_sources, l, lpos)) continue;
		if (is_cube && dot_product_ptv(norm, lpos, center) <= 0.0) continue; // back-facing polygon (thick > 0)?

		if (thick_test) { // test top and bottom surfaces
			xmin = MESH_X_SIZE;
			ymin = MESH_Y_SIZE;
			xmax = 0;
			ymax = 0;
			int x1, y1, x2, y2;
			bool got_bb(0);

			for (unsigned i = 0; i < 2; ++i) {
				if (get_shape_shadow_bb(&(pts[i].front()), npoints, l, quality, lpos, x1, y1, x2, y2, ret_val, SHADOW_TYPE)) {
					got_bb = 1;
					xmin   = min(xmin, x1);
					ymin   = min(ymin, y1);
					xmax   = max(xmax, x2);
					ymax   = max(ymax, y2);
				}
			}
			if (!got_bb) continue;
		}
		else {
			if (!get_shape_shadow_bb(points, npoints, l, quality, lpos, xmin, ymin, xmax, ymax, ret_val, SHADOW_TYPE)) continue;
		}
		/*if (quality == 0 && xmin == xmax && ymin == ymax) {
			shadow_mask[l][ymin][xmin] |= SHADOW_TYPE;
			continue;
		}*/
		if (is_dynamic) dshadow_lights |= (1 << l);
		if (thick_test) test_side = (dot_product(vector3d(lpos, center), norm) > 0.0); // inexact but probably OK
		
		for (int i = ymin; i <= ymax; ++i) { // fast approximation - uses hybrid ray casting and intersection test area fill
			for (int j = xmin; j <= xmax; ++j) {
				if (shadow_mask[l][i][j] & SHADOW_TYPE) continue;
				point const pt(get_xval(j), get_yval(i), mesh_height[i][j]);
				vector3d const v1(lpos, pt);
				bool shadowed(0);

				if (thick_test) { // test extruded (3D) polygon
					if (thick_poly_intersect(v1, pt, norm, pts, test_side, npoints)) shadowed = 1;
				}
				else { // test planar (2D) polygon
					float t;
					if (line_poly_intersect(pt, lpos, points, npoints, norm, t)) {
						if (tid < 0 || npoints != 4 || !is_billboard_texture_transparent(points, (pt + v1*t), tid)) shadowed = 1;
					}
				}
				if (shadowed) shadow_mask[l][i][j] |= SHADOW_TYPE;
			}
		}
	}
	return ret_val;
}


int cube_shadow(cube_t const &cube, char light_sources, int is_dynamic, int quality) {

	cube_t c(cube);

	if (CLIP_SHADOW_CUBES) { // clip to xy scene - not perfect but better than nothing
		c.d[0][0] = max(-X_SCENE_SIZE, c.d[0][0]);
		c.d[1][0] = max(-Y_SCENE_SIZE, c.d[1][0]);
		c.d[0][1] = min( X_SCENE_SIZE, c.d[0][1]);
		c.d[1][1] = min( Y_SCENE_SIZE, c.d[1][1]);
	}
	int xmin, xmax, ymin, ymax, ret_val(0);
	unsigned char const SHADOW_TYPE(is_dynamic ? DYNAMIC_SHADOW : OBJECT_SHADOW);
	point lpos, pts[8];
	get_cube_points(c.d, pts);

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid(light_sources, l, lpos)) continue;
		if (!get_shape_shadow_bb(pts, 8, l, quality, lpos, xmin, ymin, xmax, ymax, ret_val, SHADOW_TYPE)) continue;
		
		if (quality == 0 && xmin == xmax && ymin == ymax) {
			shadow_mask[l][ymin][xmin] |= SHADOW_TYPE;
			continue;
		}
		if (is_dynamic) dshadow_lights |= (1 << l);

		for (int i = ymin; i <= ymax; ++i) { // fast approximation - uses hybrid ray casting and intersection test area fill
			for (int j = xmin; j <= xmax; ++j) {
				if (shadow_mask[l][i][j] & SHADOW_TYPE) continue;
				point const pt(get_xval(j), get_yval(i), mesh_height[i][j]);
				if (check_line_clip(pt, lpos, cube.d)) shadow_mask[l][i][j] |= SHADOW_TYPE;
			}
		}
	}
	return ret_val;
}


int check_shadow_edge_clip(point const &pt, point const &lpos, int &xmin, int &xmax, int &ymin, int &ymax) {

	vector3d const dir((pt - lpos).get_norm());
	point pts[2] = {lpos, (lpos + dir*(10.0*FAR_CLIP))}; // extend pt off to infinity (well, very far)

	if (do_line_clip_scene(pts[0], pts[1], zbottom, max(lpos.z, ztop))) {
		for (unsigned i = 1; i < 2; ++i) { // iteration bounds have a significant effect on shadow time/quality
			int const xp(get_xpos_clamp(pts[i].x)), yp(get_ypos_clamp(pts[i].y));
			ADD_POINT_BB(xp, yp);
		}
		return 1;
	}
	return 0;
}


void reset_shadows(unsigned char type) {

	unsigned char not_type(~type);

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (type == DYNAMIC_SHADOW && !(dshadow_lights & (1 << l))) continue;
		
		for (int i = 0; i < MESH_Y_SIZE; ++i) {
			for (int j = 0; j < MESH_X_SIZE; ++j) {
				shadow_mask[l][i][j] &= not_type;
			}
		}
	}
	if (type == DYNAMIC_SHADOW) dshadow_lights = 0;
}



