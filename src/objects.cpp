// 3D World - OpenGL CS184 Computer Graphics Project
// obj_group and coll_obj class definitions
// by Frank Gennari
// 1/11/06

#include "3DWorld.h"
#include "mesh.h"
#include "transform_obj.h"
#include "physics_objects.h"
#include "collision_detect.h"


float const NDIV_SCALE = 120.0;


extern int draw_model, display_mode, destroy_thresh, xoff2, yoff2;
extern float temperature, tfticks;
extern unsigned ALL_LT[];
extern obj_type object_types[];
extern dwobject def_objects[];
extern texture textures[];
extern platform_cont platforms;


// ******************* COLL_OBJ MEMBERS ******************


void coll_obj::init() {

	cp.coll_func = NULL;
	cp.tscale    = 1.0;
	cp.surfs     = 0;
	cp.specular  = 0.0;
	cp.shine     = 1.0;
	norm         = zero_vector;
	npoints      = 0;
	type         = COLL_NULL;
	status       = COLL_UNUSED;
	destroy      = 1;
	fixed        = 1;
	id           = -1;
}


void coll_obj::clear_internal_data() {

	fixed    = 0; // unfix it so that it's actually removed
	lighted  = 0;
	cp.surfs = 0;
	occluders.clear();
}


void coll_obj::print_bounds() const {

	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			cout << d[i][j] << ",";
		}
		cout << " ";
	}
}


void coll_obj::bb_union(float bb[3][2], int init) {

	for (unsigned i = 0; i < 3; ++i) {
		bb[i][0] = (init ? d[i][0] : min(bb[i][0], d[i][0]));
		bb[i][1] = (init ? d[i][1] : max(bb[i][1], d[i][1]));
	}
}


void coll_obj::calc_size() {

	switch (type) {
	case COLL_CUBE: // use bbox
		volume = get_volume();
		break;
	case COLL_SPHERE:
		volume = (4.0/3.0)*PI*radius*radius*radius;
		break;
	case COLL_POLYGON:
		assert(npoints == 3 || npoints == 4);
		volume = triangle_area(points);
		if (npoints == 4) {
			point const points2[3] = {points[2], points[3], points[0]};
			volume += triangle_area(points2); // other triangle
		}
		volume *= thickness;
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		volume = PI*(radius*radius+radius*radius2+radius2*radius2)*p2p_dist(points[0], points[1])/3.0;
		break;
	default: assert(0);
	}
}


float coll_obj::calc_min_dim() const {

	float min_dim(min_len());
	if (type == COLL_POLYGON) min_dim = min(min_dim, thickness);
	return min_dim;
}


bool coll_obj::clip_in_2d(float const bb[2][2], float &val, int d1, int d2, int dir) const {

	assert(d1 >= 0 && d1 <= 3 && d2 >= 0 && d2 <= 3 && d1 != d2 && (dir == 0 || dir == 1));
	int d3(0);
	for (int i = 0; i < 3; ++i) {if (i != d1 && i != d2) d3 = i;}

	switch (type) {
	case COLL_CUBE:
		val = d[d3][dir];
		return 1;

	case COLL_SPHERE:
		{
			cube_t cube;

			for (unsigned i = 0; i < 2; ++i) {
				cube.d[d1][i] = bb[0][i];
				cube.d[d2][i] = bb[1][i];
				cube.d[d3][i] = (i ? SCENE_SIZE[d3] : -SCENE_SIZE[d3]);
			}
			if (!sphere_cube_intersect(points[0], radius, cube)) return 0; //val = d[d3][dir];
			val = points[0][d3] + pow((PI/6.0), (1.0/3.0))*(dir ? radius : -radius); // doesn't seem to help
		}
		return 1;

	case COLL_CYLINDER:
		if (d3 == 2) {
			val = d[d3][dir]; // can be inaccurate, especially if radius1 != radius2 (cone, etc.)
		}
		else { // unused, untested
			int const bb_z((d1 == 2) ? d1 : d2);
			float const z(0.5*(bb[bb_z][0] + bb[bb_z][1]));
			float const t((z - points[0].z)/(points[1].x - points[0].z));
			float const r(radius + t*(radius2 - radius));
			val = points[0][d3] + sqrt(PI/4.0)*(dir ? r : -r); // approximate with a square of the same cross sectional area
		}
		return 1;

	case COLL_CYLINDER_ROT:
		val = d[d3][dir];
		// should really do something with this - can be very inaccurate
		return 1;

	case COLL_POLYGON:
		{
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
			
			if (fabs(norm[d3]) > 0.01 && thickness <= MIN_POLY_THICK2) { // doesn't work on thick or vertical polygons
				float const dval(-dot_product(norm, points[0]));
				float const cent[2] = {0.5*(bb[0][0] + bb[0][1]), 0.5*(bb[1][0] + bb[1][1])};
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
	case COLL_CUBE:
		npoints = 1; break;
	case COLL_SPHERE:
		npoints = 1; break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		npoints = 2; break;
	case COLL_POLYGON:
		break; // 3 or 4, should be set at another time
	default: assert(0);
	}
}


void setup_sphere_cylin_texgen(float s_scale, float t_scale, vector3d const &dir, bool as_attr) { // dir does not need to be normalized

	int const dim(get_max_dim(dir));
	point p1, p2;
	
	for (unsigned i = 0; i < 3; ++i) {
		p1[i] = (i == dim) ? s_scale : 0.0;
		p2[i] = (i == dim) ? 0.0     : t_scale;
	}
	setup_texgen_full(p1.x, p1.y, p1.z, 0.0, p2.x, p2.y, p2.z, 0.0, GL_EYE_LINEAR, as_attr);
}


void coll_obj::draw_cobj(unsigned i, int &last_tid, int &last_group_id, int &last_pri_dim) { // non-const: modifies shadow state

	if (no_draw()) return;
	assert(status != COLL_FREED && !disabled());
	point center;
	float brad;
	bounding_sphere(center, brad);
	if ((type == COLL_CUBE) ? !camera_pdu.cube_visible(*this) : !sphere_in_camera_view(center, brad, 0)) return;
	
	if ((display_mode & 0x08) && !occluders.empty()) {
		point pts[8];
		point const camera(get_camera_pos());
		unsigned const ncorners(get_cube_corners(d, pts, camera, 0)); // 8 corners allocated, but only 6 used
		if (is_occluded(occluders, pts, ncorners, camera)) return;
	}
	// we want everything to be textured for simplicity in code/shaders,
	// so if there is no texture specified just use a plain white texture
	int const tid((cp.tid >= 0) ? cp.tid : WHITE_TEX);
	float const ar(get_tex_ar(tid));
	if (lighted == COBJ_LIT_UNKNOWN) lighted = COBJ_LIT_FALSE;

	// process groups
	int const pri_dim(::get_max_dim(norm));
	bool const in_group(group_id >= 0), same_group(group_id == last_group_id && tid == last_tid && pri_dim == last_pri_dim);
	bool const start_group(in_group && !same_group), end_group(last_group_id >= 0 && !same_group);
	last_group_id = group_id;
	last_pri_dim  = pri_dim;
	
	if (end_group) {
		glEnd();
	}
	if (!in_group || start_group) { // should be the same across groups
		set_specular(cp.specular, cp.shine);
		set_color_d(cp.color); // set material diffuse
		colorRGBA(0.0, 0.0, 0.0, cp.color.alpha).do_glColor();
	}
	if (tid != last_tid) {
		bool const textured(select_texture(tid));
		assert(textured);
		last_tid = tid;
	}
	if (start_group) {
		float const tscale[2] = {cp.tscale, ar*cp.tscale}, xlate[2] = {cp.tdx, cp.tdy};
		vector3d tex_dir(0,0,0);
		tex_dir[pri_dim] = 1.0;
		setup_polygon_texgen(tex_dir, tscale, xlate, cp.swap_txy, USE_ATTR_TEXGEN);
		glBegin(GL_TRIANGLES);
	}
	if (in_group) { // FIXME: color bug when using dynamic lighting
		assert(type == COLL_POLYGON && thickness <= MIN_POLY_THICK2); // thin triangle/quad
		vector3d const normal(get_norm_camera_orient(norm, center));
		unsigned const ixs[6] = {0,1,2,0,2,3};
		unsigned const nix((npoints == 3) ? 3 : 6); // triangle or quad (2 tris)
		
		for (unsigned i = 0; i < nix; ++i) {
			normal.do_glNormal(); // FIXME: smooth?
			points[ixs[i]].do_glVertex();
		}
		return;
	}
	switch (type) {
	case COLL_CUBE:
		draw_coll_cube((draw_model == 0), tid);
		break;

	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		{
			float const scale(NDIV_SCALE*get_zoom_scale());
			float const size(scale*sqrt(((max(radius, radius2) + 0.002)/min(distance_to_camera(center),
				min(distance_to_camera(points[0]), distance_to_camera(points[1]))))));
			int const ndiv(min(N_CYL_SIDES, max(3, (int)size)));
			bool const draw_ends(!(cp.surfs & 1));
			setup_sphere_cylin_texgen(cp.tscale, ar*cp.tscale, (points[1] - points[0]), USE_ATTR_TEXGEN);
			draw_fast_cylinder(points[0], points[1], radius, radius2, ndiv, 0, (draw_ends && tid < 0)); // Note: using texgen, not textured

			if (draw_ends && tid >= 0) { // draw ends with different texture matrix
				float const tscale[2] = {cp.tscale, get_tex_ar(tid)*cp.tscale}, xlate[2] = {cp.tdx, cp.tdy};
				setup_polygon_texgen((points[1] - points[0]).get_norm(), tscale, xlate, 0, USE_ATTR_TEXGEN);
				// FIXME: Not exactly correct, we're redrawing the sides here as well but there are texgen issues if we don't
				draw_fast_cylinder(points[0], points[1], radius, radius2, ndiv, 0, 1); // Note: using texgen, not textured
			}
		}
		break;

	case COLL_SPHERE:
		{
			float const scale(0.7*NDIV_SCALE*get_zoom_scale()), size(scale*sqrt((radius + 0.002)/distance_to_camera(points[0])));
			int const ndiv(min(N_SPHERE_DIV, max(5, (int)size)));
			setup_sphere_cylin_texgen(cp.tscale, ar*cp.tscale, plus_z, USE_ATTR_TEXGEN);
			draw_subdiv_sphere(points[0], radius, ndiv, 0, 1); // Note: using texgen, not textured
		}
		break;

	case COLL_POLYGON:
		draw_extruded_polygon(tid);
		break;
	}
}


int coll_obj::simple_draw(int ndiv, int in_cur_prim, bool no_normals, bool in_dlist) const {

	switch (type) {
	case COLL_CUBE:
		in_cur_prim = draw_simple_cube(*this, 0, in_cur_prim, no_normals, cp.surfs);
		break;

	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		{
			bool const draw_ends(!(cp.surfs & 1));

			if (no_normals && !in_dlist && ndiv == 3 && !draw_ends) { // special case draw as quad
				in_cur_prim = draw_cylin_quad_proj(cylinder_3dw(points[0], points[1], radius, radius2),
					((points[0] + points[1])*0.5 - get_camera_pos()), in_cur_prim, 1);
			}
			else {
				if (in_cur_prim != PRIM_DISABLED) {
					if (in_cur_prim >= 0) glEnd();
					in_cur_prim = PRIM_UNSET;
				}
				draw_fast_cylinder(points[0], points[1], radius, radius2, ndiv, 0, draw_ends);
			}
		}
		break;

	case COLL_SPHERE:
		if (in_cur_prim != PRIM_DISABLED) {
			if (in_cur_prim >= 0) glEnd();
			in_cur_prim = PRIM_UNSET;
		}
		//if (no_normals && ndiv == 3) {} // draw as circle/texture?
		if (in_dlist) {
			draw_subdiv_sphere(points[0], radius, ndiv, 0, 1);
		}
		else {
			draw_sphere_dlist(points[0], radius, ndiv, 0); // faster, but doesn't work when in dlist
		}
		break;

	case COLL_POLYGON:
		if (thickness <= MIN_POLY_THICK2) {
			in_cur_prim = draw_simple_polygon(points, npoints, norm, in_cur_prim, no_normals);
		}
		else {
			in_cur_prim = draw_simple_extruded_polygon(thickness, points, npoints, in_cur_prim, no_normals); // pass in norm?
		}
		break;
	}
	return in_cur_prim;
}


bool coll_obj::is_semi_trans() const {

	return (cp.color.alpha < 1.0 || (cp.tid >= 0 && textures[cp.tid].ncolors == 4));
}


bool coll_obj::is_player() const { // sort of a hack
	
	return ((cp.coll_func == camera_collision) || (cp.coll_func == smiley_collision));
}


bool coll_obj::is_invis_player() const { // sort of a hack
	
	if (cp.coll_func == smiley_collision && has_invisibility(cp.cf_index)) return 1;
	if (cp.coll_func == camera_collision && has_invisibility(-1))          return 1;
	return 0;
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
	case COLL_POLYGON:
		center = get_center(points, npoints);
		brad   = radius;
		break;
	default:
		assert(0);
	}
}


bool coll_obj::truly_static() const {

	if (status != COLL_STATIC || destroy >= max((int)SHATTERABLE, destroy_thresh+1)) return 0;
	if (platform_id < 0) return 1;
	return 0;
	//assert((unsigned)platform_id < platforms.size());
	//return !platforms[platform_id].is_moving();
}


bool coll_obj::can_be_scorched() const {

	return (status == COLL_STATIC && !(cp.tid >= 0 && textures[cp.tid].ncolors == 4) && !no_draw()); // allow destroyable and transparent objects
}


point coll_obj::get_center_pt() const {

	switch (type) {
	case COLL_CUBE:
	case COLL_SPHERE:
		return points[0];
	case COLL_CYLINDER:
		return point(points[0].x, points[0].y, 0.5*(points[0].z + points[1].z));
	case COLL_CYLINDER_ROT:
		return get_center_n2(points);
	case COLL_POLYGON:
		return get_center(points, npoints);
	default:
		assert(0);
	}
	return all_zeros; // never gets here
}


float coll_obj::get_max_dim() const {

	float md(0.0);

	for (unsigned i = 0; i < 3; ++i) {
		md = max(md, fabs(d[i][1] - d[i][0]));
	}
	return md;
}


float coll_obj::get_light_transmit(point v1, point v2) const {

	if (type != COLL_CUBE)        return 1.0; // only cubes are supported due to clipping issues
	if (cp.light_atten == 0.0)    return 1.0;
	if (!do_line_clip(v1, v2, d)) return 1.0;
	return exp(-cp.light_atten*p2p_dist(v1, v2));
}


bool coll_obj::has_poly_billboard_alpha() const {

	if (!is_billboard || type != COLL_POLYGON || thickness > MIN_POLY_THICK2 || npoints != 4) return 0;
	if (cp.tid < 0 || textures[cp.tid].ncolors != 4) return 0; // no alpha channel texture
	return 1;
}


bool coll_obj::check_poly_billboard_alpha(point const &p1, point const &p2, float t) const {

	return (!has_poly_billboard_alpha() || !is_billboard_texture_transparent(points, (p1 + (p2 - p1)*t), cp.tid));
}


// ******************* OBJ_GROUP MEMBERS ******************


void obj_group::create(int obj_type_, unsigned max_objects_, unsigned init_objects_,
					   unsigned app_rate_, bool init_enabled_, bool reorderable_, bool auto_max)
{
	flags        = JUST_INIT;
	type         = obj_type_;
	app_rate     = app_rate_;
	reorderable  = reorderable_;
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
	if (max_app > 0) app_rate = min(app_rate, max_app);
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

		if (j < init_objects) {
			gen_object_pos(objects[j].pos, object_types[type].flags);
		}
		else {
			objects[j].status = 0;
		}
	}
	flags &= (PRECIPITATION | APP_FROM_LT); // keep only this flag
}


// normally called before using objects, but can be called dynamically later
void obj_group::add_predef_obj(point const &pos, int type, int rtime) {
	
	predef_objs.push_back(predef_obj(pos, type, rtime));
	reorderable = 0; // need to unset reorderable so that predef_objs indexes remain correct
}


void obj_group::preproc_this_frame() {

	unsigned const nobjs(max_objects());
	end_id = nobjs;
	new_id = 0;
	if (!enabled) return;

	if (reorderable) { // some objects such as smileys are position dependent
		assert(predef_objs.empty());
		unsigned saw_id(0);
		
		for (unsigned j = 0; j < nobjs; ++j) {
			if (objects[j].status != 0) saw_id = j+1;
		}
		sort(objects.begin(), (objects.begin() + saw_id));
		
		for (unsigned j = 0; j < nobjs; ++j) {
			if (objects[j].status == 0) {
				end_id = j;
				break;
			}
		}
	}
	for (vector<predef_obj>::iterator i = predef_objs.begin(); i != predef_objs.end(); ++i) {
		if (i->obj_used == -1) continue;
		assert((unsigned)i->obj_used < max_objects());
			
		if (objects[i->obj_used].disabled()) { // unused object
			i->obj_used = -1; // reset back to 'unused'
			i->cur_time = tfticks;
		}
	}
}


void obj_group::remove_reset_cobjs() {

	if (!large_radius()) return;

	for (unsigned j = 0; j < max_objects(); ++j) {
		if (objects[j].status != 0 && objects[j].coll_id >= 0) { // is status check correct?
			remove_reset_coll_obj(objects[j].coll_id);
		}
	}
}


void obj_group::create_object_at(unsigned i, point const &pos) {

	enable();
	assert(max_objects() == max_objs && i < max_objects());
	assert(type >= 0 && type < NUM_TOT_OBJS);
	assert(!is_nan(pos));
	if (objects[i].coll_id >= 0) remove_reset_coll_obj(objects[i].coll_id); // just in case
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
		if (!td) td = new transform_data;
		if (tflags & OBJ_ROLLS)     td->matrices.resize(max_objs);
		if (tflags & VERTEX_DEFORM) td->perturb_maps.resize(max_objs);
	}
	enabled = 1;
}


void obj_group::disable() {

	if (enabled == 0) return;
	remove_reset_cobjs();
	objects.clear();
	if (!large_radius()) remove_excess_cap(objects); // free

	if ((object_types[type].flags & (OBJ_ROLLS | VERTEX_DEFORM))) {
		delete td;
		td = NULL;
	}
	for (vector<predef_obj>::iterator i = predef_objs.begin(); i != predef_objs.end(); ++i) {
		i->obj_used = -1;
	}
	enabled = 0;
}


void obj_group::free() {

	remove_reset_cobjs();
	if (!objects.empty()) reset_status(objects);

	if (td) {
		delete td;
		td = NULL;
	}
}


void obj_group::shift(vector3d const &vd) {

	if (!enabled) return;

	for (unsigned j = 0; j < max_objects(); ++j) {
		if (!objects[j].disabled()) {
			objects[j].pos += vd;
			if (type == SMILEY) shift_player_state(vd, j);
		}
	}
	for (unsigned i = 0; i < predef_objs.size(); ++i) {
		predef_objs[i].pos += vd;
	}
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



