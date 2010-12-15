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


extern int draw_model, display_mode, destroy_thresh, do_zoom, xoff2, yoff2;
extern float temperature;
extern obj_type object_types[];
extern dwobject def_objects[];
extern texture textures[];


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


bool coll_obj::clear_lightmap_entry(lvmap::iterator it, int mode, unsigned keep, vector<pair<quad_div, lv_val> > *to_add) {

	unsigned char const tag(it->first.tag);
	unsigned &status(it->second.status);

	if (mode == 2) { // update shadows - can only add shadows
		if (require_either_lightval(status, 1, 4)) return 0; // all hidden or all shadowed
	}
	if (tag == QD_TAG_DLIST && status != 0) { // display list
		if (glIsList(status)) glDeleteLists(status, 1);
		status = 0; // erase?
	}
	if (tag == QD_TAG_TEXTURE && status != 0) { // texture
		free_texture(status); // erase?
	}
	if (mode != 1) {
		if (keep > 0 && (tag == QD_TAG_QUAD || tag == QD_TAG_TRIANGLE)) { // backup precomputed lighting values
			assert(to_add);
			assert(keep < (1 << 8)); // only support up to 8 light sources for now
			quad_div qd(it->first);
			lv_val lv(it->second);

			for (unsigned L = 0; L < 8; ++L) {
				if (!(keep & (1 << L))) set_light_val(lv.status, 0, L); // mark as reset
			}
			qd.tag |= QD_TAG_OLD;
			to_add->push_back(make_pair(qd, lv));
			return 1;
		}
		delete [] it->second.nvals; // clear all
		it->second.nvals = NULL; // just in case
		if (keep > 0) return 1;
	}
	return (mode == 2);
}


// mode: 0 = clear all, 1 = clear dlists only, 2 = update shadows (keep is only used in mode 0)
void coll_obj::clear_lightmap(int mode, unsigned keep, bool keep_depends) {

	vector<pair<quad_div, lv_val> > to_add;

	for (lvmap::iterator it = lightmap.begin(); it != lightmap.end();) {
		bool const deleted(clear_lightmap_entry(it, mode, keep, &to_add));
		lvmap::iterator to_del(it++);
		if (deleted) lightmap.erase(to_del);
	}
	if (mode == 0 && keep == 0) {
		assert(to_add.empty());
		if (!keep_depends) shadow_depends.clear();
		lightmap.clear();
		lighted = COBJ_LIT_UNKNOWN;
	}
	for (unsigned i = 0; i < to_add.size(); ++i) { // add back in old values
		lightmap.insert(to_add[i]);
	}
}


void coll_obj::clear_lightmap_if_lighted_eq(int shadowed, int partial) {

	for (lvmap::iterator it = lightmap.begin(); it != lightmap.end();) {
		int const lv(it->second.status); // lighted: 1 = shadowed, 2 = unshadowed, 3 = partially shadowed
		unsigned char const tag(it->first.tag);
		lvmap::iterator temp(it++);

		if (tag == QD_TAG_GLOBAL || tag == QD_TAG_DLIST || tag == QD_TAG_TEXTURE ||
			(lv == 1 && shadowed) || (lv == 2 && !shadowed) || (lv == 3 && partial))
		{
			clear_lightmap_entry(temp, 0, 0);
			lightmap.erase(temp);
		}
	}
}


void coll_obj::clear_dependent_cobjs_lightmaps(vector<coll_obj> &cobjs, unsigned ix) const {

	for (set<int>::const_iterator it = shadow_depends.begin(); it != shadow_depends.end(); ++it) {
		assert((size_t)*it < cobjs.size() && *it != (int)ix); // check for circular reference
		cobjs[*it].clear_lightmap(0, 0, 0);
	}
}


void coll_obj::update_shadowed_cobjs(vector<coll_obj> &cobjs, vector<int> const &indices, unsigned ix) const {

	for (set<int>::const_iterator it = shadow_depends.begin(); it != shadow_depends.end(); ++it) {
		assert((size_t)*it < cobjs.size() && *it != (int)ix); // check for circular reference
		coll_obj &cobj(cobjs[*it]);

		// optimization: check if target object is still shadowed by a fragment of the original source object
		if (cobj.type == COLL_CUBE && !indices.empty()) {
			bool shadowed(0);
			point pts[8];
			get_cube_points(cobj.d, pts);
			point const lpos(get_light_pos());

			for (unsigned m = 0; m < indices.size() && !shadowed; ++m) {
				assert((size_t)indices[m] < cobjs.size());
				shadowed = 1;

				for (unsigned j = 0; j < 8; ++j) {
					if (!cobjs[indices[m]].line_intersect(pts[j], lpos)) shadowed = 0;
				}
				if (shadowed) cobjs[indices[m]].shadow_depends.insert(*it);
			}
			if (shadowed) continue; // got into the above if (shadowed) {} statement
		}
		cobj.clear_lightmap_if_lighted_eq(1, 1); // clear if shadowed or partially shadowed
	}
}


void coll_obj::clear_internal_data(vector<coll_obj> &cobjs, vector<int> const &indices, unsigned ix) {

	update_shadowed_cobjs(cobjs, indices, ix);
	fixed    = 0; // unfix it so that it's actually removed
	lighted  = 0;
	cp.surfs = 0;
	shadow_depends.clear();
	clear_lightmap(0);
	occluders.clear();
	sobjs.clear();
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
	case COLL_POLYGON:
		volume = (4.0/3.0)*PI*radius*radius*radius;
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		volume = PI*(radius*radius+radius*radius2+radius2*radius2)*p2p_dist(points[0], points[1])/3.0;
		break;
	default: assert(0);
	}
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


void setup_sphere_cylin_texgen(float s_scale, float t_scale, vector3d const &dir) { // dir does not need to be normalized

	int const dim(get_max_dim(dir));
	point p1, p2;
	
	for (unsigned i = 0; i < 3; ++i) {
		p1[i] = (i == dim) ? s_scale : 0.0;
		p2[i] = (i == dim) ? 0.0     : t_scale;
	}
	setup_texgen_full(p1.x, p1.y, p1.z, 0.0, p2.x, p2.y, p2.z, 0.0, GL_EYE_LINEAR);
}


void coll_obj::draw_cobj(unsigned i, int &last_tid) { // non-const: modifies shadow state

	if (no_draw()) return;
	assert(status != COLL_FREED && !disabled());
	point center;
	float brad;
	bounding_sphere(center, brad);
	if (!sphere_in_camera_view(center, brad, 0)) return;
	
	if ((display_mode & 0x08) && !occluders.empty()) {
		point pts[8];
		point const camera(get_camera_pos());
		unsigned const ncorners(get_cube_corners(d, pts, camera, 0)); // 8 corners allocated, but only 6 used
		if (is_occluded(occluders, pts, ncorners, camera)) return;
	}
	//if (brad/distance_to_camera(center) < 0.01) return; // too far/small
	// we want everything to be textured for simplicity in code/shaders,
	// so if there is no texture specified just use a plain white texture
	int const tid((cp.tid >= 0) ? cp.tid : WHITE_TEX);
	float const ar(get_tex_ar(tid));
	bool const no_lighting(cp.color == BLACK && cp.specular == 0.0);
	if (lighted == COBJ_LIT_UNKNOWN) lighted = COBJ_LIT_FALSE;
	set_specular(cp.specular, cp.shine);
	set_color_d(cp.color); // set material ambient and diffuse
	colorRGBA(0.0, 0.0, 0.0, cp.color.alpha).do_glColor();

	if (tid != last_tid) {
		bool const textured(select_texture(tid));
		assert(textured);
		last_tid = tid;
	}
	switch (type) {
	case COLL_CUBE:
		draw_coll_cube(ar, (draw_model == 0), i, tid);
		break;

	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		{
			point const p1(points[0]), p2(points[1]);
			float const scale(NDIV_SCALE*get_zoom_scale());
			float const size(scale*sqrt(((max(radius, radius2) + 0.002)/min(distance_to_camera(center),
				min(distance_to_camera(p1), distance_to_camera(p2))))));
			int const nsides(min(N_CYL_SIDES, max(3, (int)size)));
			setup_sphere_cylin_texgen(cp.tscale, ar*cp.tscale, (p2 - p1));
			draw_subdiv_cylinder(p1, p2, radius, radius2, nsides, 1, !(cp.surfs & 1), (cp.surfs == 1), i, no_lighting, tid);
		}
		break;

	case COLL_SPHERE:
		{
			float const scale(0.7*NDIV_SCALE*get_zoom_scale()), size(scale*sqrt((radius + 0.002)/distance_to_camera(points[0])));
			int const nsides(min(N_SPHERE_DIV, max(5, (int)size)));
			setup_sphere_cylin_texgen(cp.tscale, ar*cp.tscale, plus_z);
			draw_subdiv_sphere_at(points[0], radius, nsides, i, no_lighting, tid);
		}
		break;

	case COLL_POLYGON:
		draw_extruded_polygon(thickness, points, NULL, npoints, norm, i, tid);
		break;
	}
}


bool coll_obj::is_semi_trans() const {

	return (cp.color.alpha < 1.0 || (cp.tid >= 0 && textures[cp.tid].ncolors == 4));
}


bool coll_obj::is_player() const { // sort of a hack
	
	return ((cp.coll_func == camera_collision) || (cp.coll_func == smiley_collision));
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
		center = ::get_center(points, npoints);
		brad   = radius;
		break;
	default:
		assert(0);
	}
}


bool coll_obj::truly_static() const {

	return (status == COLL_STATIC && (type != COLL_CUBE || destroy <= destroy_thresh) && platform_id < 0);
}


bool coll_obj::can_be_scorched() const {

	return (truly_static() && !is_semi_trans() && !no_draw());
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
		return ::get_center(points, npoints);
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


void obj_group::sort_and_calc_end() {

	unsigned const nobjs(max_objects());
	end_id = nobjs;
	new_id = 0;

	if (enabled && reorderable) { // some objects such as smileys are position dependent
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



