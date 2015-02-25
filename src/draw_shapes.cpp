// 3D World - Dynamic Surface Subdivision, Lighting, and Shadows
// by Frank Gennari
// 4/10/05
#include "3DWorld.h"
#include "mesh.h"
#include "dynamic_particle.h"
#include "physics_objects.h"
#include "shaders.h"


float const TOLER_ = 1.0E-6;


extern int begin_motion, num_groups, camera_coll_id, spectate, display_mode, camera_mode, camera_view;
extern float zmin, NEAR_CLIP;
extern double camera_zh;
extern vector<int> weap_cobjs;
extern vector<unsigned> falling_cobjs;
extern coll_obj_group coll_objects;
extern platform_cont platforms;
extern obj_type object_types[];
extern obj_group obj_groups[];

vector<shadow_sphere> shadow_objs;

bool no_sparse_smap_update();


shadow_sphere::shadow_sphere(point const &pos0, float radius0, int cid0) : sphere_t(pos0, radius0), cid(cid0) {

	if (cid < 0) {
		ctype = COLL_SPHERE; // sphere is the default
	}
	else {
		assert(size_t(cid) < coll_objects.size());
		ctype = coll_objects[cid].type;
	}
}


bool shadow_sphere::line_intersect_cobj(point const &p1, point const &p2) const {

	assert(cid >= 0 && cid < (int)coll_objects.size());
	return coll_objects[cid].line_intersect(p1, p2);
}


bool shadow_sphere::test_volume_cobj(point const *const pts, unsigned npts, point const &lpos) const {

	coll_obj const &c(coll_objects[cid]);
	if (ctype == COLL_SPHERE && (pos != c.points[0] || radius != c.radius)) return 1; // camera sphere != pos
	return !c.cobj_plane_side_test(pts, npts, lpos);
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


bool check_face_containment(cube_t const &cube, int dim, int dir, int cobj) {

	assert((dim >= 0 && dim <= 2) && (dir == 0 || dir == 1));

	if (dim == 2 && dir == 0) { // z bottom
		int const x1(get_xpos(cube.d[0][0])), x2(get_xpos(cube.d[0][1]));
		int const y1(get_ypos(cube.d[1][0])), y2(get_ypos(cube.d[1][1]));
		bool under_mesh(1);
			
		for (int y = max(0, y1); y <= min(MESH_Y_SIZE-1, y2) && under_mesh; ++y) {
			for (int x = max(0, x1); x <= min(MESH_X_SIZE-1, x2) && under_mesh; ++x) {
				under_mesh &= (cube.d[2][0] < mesh_height[y][x]);
			}
		}
		if (under_mesh) return 1;
	}
	point const cent(cube.get_cube_center());
	int const x(get_xpos(cent.x)), y(get_ypos(cent.y));
	if (point_outside_mesh(x, y)) return 0;
	coll_cell const &cell(v_collision_matrix[y][x]);
	unsigned const ncv((unsigned)cell.cvals.size());

	for (unsigned i = 0; i < ncv; ++i) { // test for internal faces to be removed
		coll_obj const &c(coll_objects[cell.cvals[i]]);
		if (c.type != COLL_CUBE || !c.fixed || c.status != COLL_STATIC || c.platform_id >= 0 || c.destroy >= SHATTERABLE) continue;
		if ((int)cell.cvals[i] == cobj || c.is_semi_trans() || fabs(c.d[dim][!dir] - cube.d[dim][dir]) > TOLER_) continue;
		bool contained(1);

		for (unsigned k = 0; k < 2 && contained; ++k) {
			unsigned const dk((dim+k+1)%3);
			if (cube.d[dk][0] < (c.d[dk][0]-TOLER_) || cube.d[dk][1] > (c.d[dk][1]+TOLER_)) contained = 0;
		}
		if (contained) return 1;
	}
	return 0;
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


void cobj_draw_buffer::on_new_obj_layer(obj_layer const &layer) {

	if (layer.color != last_layer.color || layer.tid != last_layer.tid ||
		layer.spec_color != last_layer.spec_color || layer.shine != last_layer.shine) {flush();}
	last_layer = layer;
}


void coll_obj::draw_coll_cube(int tid, cobj_draw_buffer &cdb) const {

	int const sides((int)cp.surfs);
	if (sides == EF_ALL) return; // all sides hidden
	bool const back_face_cull(!is_semi_trans()); // no alpha
	point const pos(points[0]), camera(get_camera_pos());
	bool inside(!back_face_cull);
	float const tscale[2] = {cp.tscale, get_tex_ar(tid)*cp.tscale};

	if (!inside) { // check if in the camera's view volume intersects the cube - if so we must render all faces
		float const dist(NEAR_CLIP + CAMERA_RADIUS);
		inside = 1;

		for (unsigned i = 0; i < 3 && inside; ++i) {
			if (camera[i] <= d[i][0]-dist || camera[i] >= d[i][1]+dist) inside = 0;
		}
	}
	pair<float, unsigned> faces[6];
	for (unsigned i = 0; i < 6; ++i) {faces[i].second = i;}

	if (!back_face_cull) { // semi-transparent
		for (unsigned i = 0; i < 6; ++i) {
			unsigned const dim(i>>1), dir(i&1), d0((dim+1)%3), d1((dim+2)%3);
			point pos;
			pos[dim] = d[dim][dir];
			pos[d0]  = 0.5*(d[d0][0] + d[d0][1]);
			pos[d1]  = 0.5*(d[d1][0] + d[d1][1]);
			faces[i].first = -p2p_dist_sq(pos, camera); // draw ordered furthest to closest to camera
		}
		sort(faces, (faces+6));
	}
	for (unsigned i = 0; i < 6; ++i) {
		unsigned const fi(faces[i].second), dim(fi>>1), dir(fi&1);
		if ((sides & EFLAGS[dim][dir]) || (!inside && !((camera[dim] < d[dim][dir]) ^ dir))) continue; // side disabled
		unsigned const d0((dim+1)%3), d1((dim+2)%3), t0((2-dim)>>1), t1(1+((2-dim)>0));
		point pts[4], p;
		p[dim] = d[dim][dir];
		p[d0 ] = d[d0][0];
		p[d1 ] = d[d1][0]; pts[0] = p;
		p[d0 ] = d[d0][1]; pts[1] = p;
		p[d1 ] = d[d1][1]; pts[2] = p;
		p[d0 ] = d[d0][0]; pts[3] = p;
		if (!dir) {swap(pts[0], pts[3]); swap(pts[1], pts[2]);}
		vector3d normal(zero_vector);
		normal[dim] = (dir ? 1.0 : -1.0);
		texgen_params_t tp;

		for (unsigned e = 0; e < 2; ++e) {
			unsigned const tdim(e ? t1 : t0);
			bool const s_or_t(cp.swap_txy() ^ (e != 0));
			float *tg(tp.st[s_or_t]);

			if (tscale[0] == 0) { // special value of tscale=0 will result in the texture being fit exactly to the cube (mapped from 0 to 1)
				tg[tdim] = 1.0/(d[tdim][1] - d[tdim][0]);
				tg[3]    = (-d[tdim][0] + texture_offset[tdim])*tg[tdim];
			}
			else {
				tg[tdim] = tscale[e];
				tg[3]    = texture_offset[tdim]*tscale[e];
			}
		}
		for (unsigned i = 0; i < 6; ++i) { // quads (2 triangles)
			cdb.add_vert(vert_norm(pts[quad_to_tris_ixs[i]], normal), tp);
		}
	} // for i
}


bool camera_back_facing(point const *const points, int npoints, vector3d const &normal) {

	return (dot_product_ptv(normal, get_camera_pos(), get_center(points, npoints)) >= 0.0);
}


bool camera_behind_polygon(point const *const points, int npoints) {

	point const center(get_center(points, npoints)), camera(get_camera_pos());
	vector3d const dirs[2] = {vector3d(points[1], points[0]), vector3d(points[npoints-1], points[0])};
	vector3d const normal(cross_product(dirs[0], dirs[1]));
	return (dot_product_ptv(normal, camera, center) >= 0.0);
}


void coll_obj::set_poly_texgen(int tid, vector3d const &normal, shader_t &shader) const {

	if (tid < 0) return; // texturing disabled
	float const tscale[2] = {cp.tscale, get_tex_ar(tid)*cp.tscale}, xlate[2] = {cp.tdx, cp.tdy};
	setup_polygon_texgen(normal, tscale, xlate, texture_offset, cp.swap_txy(), shader, 1);
}

void coll_obj::get_polygon_tparams(int tid, vector3d const &normal, texgen_params_t &tp) const {

	float const tscale[2] = {cp.tscale, get_tex_ar(tid)*cp.tscale}, xlate[2] = {cp.tdx, cp.tdy};
	vector3d v[2];
	get_poly_texgen_dirs(normal, v);
	
	for (unsigned i = 0; i < 2; ++i) {
		bool const d((i != 0) ^ cp.swap_txy());
		UNROLL_3X(tp.st[d][i_] = tscale[i]*v[i][i_];)
		tp.st[d][3] = xlate[i] + tscale[i]*dot_product(texture_offset, v[i]);
	}
}


void coll_obj::draw_polygon(int tid, point const *points, int npoints, vector3d normal, cobj_draw_buffer &cdb) const {

	normal = get_norm_camera_orient(normal, get_center(points, npoints));
	assert(npoints == 3 || npoints == 4);
	texgen_params_t tp;
	get_polygon_tparams(tid, normal, tp);

	for (int i = 0; i < ((npoints == 3) ? 3 : 6); ++i) { // 1-2 triangles
		cdb.add_vert(vert_norm(points[quad_to_tris_ixs[i]], normal), tp);
	}
}


void get_sorted_thick_poly_faces(point pts[2][4], pair<int, unsigned> faces[6], point const *points,
	unsigned npoints, vector3d const &norm, float thick, bool bfc)
{
	gen_poly_planes(points, npoints, norm, thick, pts);
	unsigned const nsides(unsigned(npoints)+2);
	assert(nsides <= 6);
	for (unsigned i = 0; i < nsides; ++i) {faces[i] = make_pair(0, i);}

	if (!bfc) { // sort by the number of centerlines crossing the surfaces
		point centers[6];
		point const camera(get_camera_pos());

		for (unsigned i = 0; i < 2; ++i) { // front and back
			centers[i] = get_center(pts[i], npoints);
		}
		for (unsigned i = 0; i < npoints; ++i) { // sides
			unsigned const ii((i+1)%npoints);
			point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
			centers[i+2] = get_center(side_pts, 4);
		}
		for (unsigned f = 0; f < nsides; ++f) {
			for (unsigned i = 0; i < 2; ++i) { // front and back
				if (i != f && line_poly_intersect((centers[f] - camera), camera, pts[i], npoints)) {
					--faces[f].first;
					++faces[i].first;
				}
			}
			for (unsigned i = 0; i < npoints; ++i) { // sides
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
}


void coll_obj::draw_extruded_polygon(int tid, cobj_draw_buffer &cdb) const {

	assert(points != NULL && (npoints == 3 || npoints == 4));
	float const thick(fabs(thickness));
	
	if (thick <= MIN_POLY_THICK) { // double_sided = 0, relies on points being specified in the correct CW/CCW order
		draw_polygon(tid, points, npoints, norm, cdb);
	}
	else {
		point pts[2][4];
		pair<int, unsigned> faces[6];
		bool const bfc(!is_semi_trans());
		get_sorted_thick_poly_faces(pts, faces, points, npoints, norm, thick, bfc);
		bool const back_facing(bfc && camera_back_facing(pts[1], npoints, norm));
		unsigned const nsides(unsigned(npoints)+2);

		for (unsigned fi = 0; fi < nsides; ++fi) { // draw back to front
			unsigned const s(faces[fi].second);

			if (s < 2) { // draw front and back
				if (bfc && (back_facing ^ (s == 0))) continue;
				if (!s) {std::reverse(pts[s], pts[s]+npoints);}
				draw_polygon(tid, pts[s], npoints, (s ? norm : -norm), cdb); // draw bottom surface
				if (!s) {std::reverse(pts[s], pts[s]+npoints);}
			}
			else { // draw sides
				unsigned const i(s-2), ii((i+1)%npoints);
				point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};

				if (!bfc || !camera_behind_polygon(side_pts, 4)) {
					draw_polygon(tid, side_pts, 4, get_poly_norm(side_pts), cdb);
				}
			}
		}
	} // end thick case
}


void coll_obj::draw_cylin_ends(int tid, int ndiv, cobj_draw_buffer &cdb) const {

	vector3d v12; // (ce[1] - ce[0]).get_norm()
	vector_point_norm const &vpn(gen_cylinder_data(points, radius, radius2, ndiv, v12));
	texgen_params_t tp;
	get_polygon_tparams(tid, v12, tp);
	float const r[2] = {radius, radius2};

	for (unsigned i = 0; i < 2; ++i) {
		if (r[i] == 0.0) continue;
		vector3d const normal(i ? v12 : -v12);

		for (unsigned s = 0; s < (unsigned)ndiv; ++s) {
			unsigned const sn((s+1)%ndiv);
			cdb.add_vert(vert_norm(points[i],        normal), tp); // center point
			cdb.add_vert(vert_norm(vpn.p[(s <<1)+i], normal), tp);
			cdb.add_vert(vert_norm(vpn.p[(sn<<1)+i], normal), tp);
		}
	}
}


void add_shadow_obj(point const &pos, float radius, int coll_id) {
	shadow_objs.push_back(shadow_sphere(pos, radius, coll_id));
}


void add_shadow_cobj(int cid) {

	if (cid < 0) return;
	assert((unsigned)cid < coll_objects.size());
	if (coll_objects[cid].disabled() || coll_objects[cid].cp.color.alpha < MIN_SHADOW_ALPHA) return;
	point center;
	float radius;
	coll_objects[cid].bounding_sphere(center, radius);
	add_shadow_obj(center, radius, cid);
}


void add_coll_shadow_objs() {
	
	//RESET_TIME;
	shadow_objs.resize(0);
	if (!shadow_map_enabled()) return;
	point const camera(get_camera_pos());

	if ((camera_mode == 1 || camera_view == 0) && !has_invisibility(CAMERA_ID) && no_sparse_smap_update()) { // shadow the camera even when in the air (but not when dead)
		point camera_pos(camera);
		if (camera_mode == 1 && !spectate) {camera_pos.z -= 0.5*camera_zh;} // cancel out the z height that was previously added
		add_shadow_obj(camera_pos, CAMERA_RADIUS, camera_coll_id);
	}
	if (begin_motion) { // can ignore if behind camera and light in front of camera
		for (int i = 0; i < num_groups; ++i) { // can we simply use the collision objects for this?
			obj_group const &objg(obj_groups[i]);
			if (!objg.enabled || !objg.large_radius())                      continue;
			if (object_types[objg.type].color.alpha < 0.5*MIN_SHADOW_ALPHA) continue; // too low? nothing fails this yet
			float const radius(object_types[objg.type].radius);
				
			for (unsigned j = 0; j < objg.end_id; ++j) {
				dwobject const &obj(objg.get_obj(j));
				if (obj.disabled() || !objg.obj_has_shadow(j)) continue;
				add_shadow_obj(obj.pos, radius, obj.coll_id);
			}
		}
	}
	for (unsigned i = 0; i < weap_cobjs.size(); ++i) {
		add_shadow_cobj(weap_cobjs[i]);
	}
	for (platform_cont::const_iterator i = platforms.begin(); i != platforms.end(); ++i) {
		for (vector<unsigned>::const_iterator j = i->cobjs.begin(); j != i->cobjs.end(); ++j) {
			add_shadow_cobj(*j);
		}
	}
	for (unsigned i = 0; i < falling_cobjs.size(); ++i) {
		add_shadow_cobj(falling_cobjs[i]);
	}
	if (display_mode & 0x0200) d_part_sys.add_cobj_shadows();
	//PRINT_TIME(" Shadow Object Creation");
}




