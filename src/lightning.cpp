// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/13/02

#include "3DWorld.h"
#include "mesh.h"
#include "openal_wrap.h"
#include "shaders.h"
#include "physics_objects.h" // for lightning_t
#include "cube_map_shadow_manager.h"


bool   const ENABLE_LT_SHADOWS= 1;
int    const PATH_FORK_MOD    = 15;
int    const PATH_END_MOD     = 15;
float  const FORK_ATTEN       = 0.8;
float  const L_STRENGTH_MULT  = 150.0;
float  const L_DAMAGE_MULT    = 80.0;
float  const LIGHTNING_PERIOD = 100.0; // average time between strikes in ticks
float  const STEP_VARIANCE    = 0.5;
int    const DISCHARGE_RAD    = 20;
double const CHARGE_HALF_D    = 5.0;
float  const EVAP_AMOUNT      = 10.0;
float  const LITNING_ICE_EFF  = 0.5;
unsigned const MAX_LITN_FORKS = 25;
unsigned const smap_obj_id    = 0; //  0


point litning_pos;
lightning_t l_strike;
cube_map_shadow_manager lightning_smgr; // should this be inside lightning_t?

extern int is_cloudy, world_mode, iticks, DISABLE_WATER, animate2;
extern float zmin, temperature, lt_green_int, water_plane_z, ztop, fticks, XY_SCENE_SIZE;
extern vector<valley> valleys;


void do_lightning_damage(point &pos, float damage, int hit_water);


inline float get_lit_h(int xpos, int ypos) {
	float h(h_collision_matrix[ypos][xpos]);
	if (!v_collision_matrix[ypos][xpos].cvals.empty()) {h = max(h, v_collision_matrix[ypos][xpos].zmax);}
	return h;
}


void lightning_t::disable() {
	enabled = 0;
	paths.clear();
	if (ENABLE_LT_SHADOWS) {lightning_smgr.remove_obj_light(smap_obj_id);}
}

void lightning_t::gen() {

	is_cloudy = 1;

	if (!animate2) { // physics stopped - draw but don't update
		draw();
		return;
	}
	int const rnum(int(LIGHTNING_PERIOD/fticks));

	if (rnum <= 1 || (rgen.rand()%rnum) != 0) { // update timestep randomly
		if (enabled) {
			if (time < int(0.1f*LITNING_TIME/fticks)) {
				draw();
				if (animate2) {time += iticks;}
			}
			else {disable();} // time is up
		}
		return;
	}
	// generate a new lightning strike
	disable();
	time = 0;
	float const cloud_zval(CLOUD_CEILING + ztop);
	int xpos(0), ypos(0);
	float max_e(0.0), charge(0);

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			float const dist_to_ground(max(0.5f*CLOUD_CEILING, (cloud_zval - get_lit_h(j, i)))); // make sure it's positive
			float const this_e(0.03*charge_dist[i][j]/dist_to_ground); // higher potential for discharge when ground is closer
			if (this_e > max_e) {max_e = this_e; xpos = j; ypos = i;}
		}
	}
	int const x0(max(0, (xpos - DISCHARGE_RAD))), x1(min(MESH_X_SIZE-1, (xpos + DISCHARGE_RAD)));
	int const y0(max(0, (ypos - DISCHARGE_RAD))), y1(min(MESH_Y_SIZE-1, (ypos + DISCHARGE_RAD)));

	for (int i = y0; i <= y1; ++i) {
		for (int j = x0; j <= x1; ++j) {
			float const distance(sqrt(float((xpos - j)*(xpos - j) + (ypos - i)*(ypos - i))));
			float const d_charge(charge_dist[i][j]/pow(2.0, double(distance/CHARGE_HALF_D)));
			charge            += d_charge;
			charge_dist[i][j] -= d_charge;
		}
	}
	float const d_charge(charge/XY_MULT_SIZE);

	for (int i = 0; i < MESH_Y_SIZE; ++i) { // redistribute charge
		for (int j = 0; j < MESH_X_SIZE; ++j) {charge_dist[i][j] += d_charge;}
	}
	enabled = 1;
	gen_recur(point(get_xval(xpos), get_yval(ypos), cloud_zval), -plus_z, max_e, 0);
	assert(!paths.empty());
	unsigned pri_path(0), min_len(0);

	for (auto p = paths.begin(); p != paths.end(); ++p) { // find primary path
		if (!p->full_path) continue;
		unsigned const len(p->get_len());
		assert(len >= 2);
		if (min_len == 0 || len < min_len) {min_len = len; pri_path = (p - paths.begin());}
	}
	for (auto p = paths.begin(); p != paths.end(); ++p) {
		if (!p->full_path) continue;
		unsigned const len(p->get_len());
		assert(len >= min_len);
		unsigned const to_trim(p->has_child ? 0 : (len - min_len)); // only trim leaf path segments
		if (to_trim > 0 && to_trim+2 <= p->points.size()-2) {p->points.resize(p->points.size() - to_trim);} // shorten non-primary paths to length of primary path
		else {do_lightning_damage(p->points.back(), p->damage, p->hit_water);} // do damage
	}
	paths[pri_path].width *= 2.0; // main branch
	litning_pos    = hit_pos = paths[pri_path].points.back();
	litning_pos.z += 0.01;
	play_thunder(litning_pos, 4.0, 0.0);

	if (ENABLE_LT_SHADOWS) {
		float const radius(7.8*LITNING_LINEAR_I*sqrt(0.1*XY_SCENE_SIZE)); // see add_dynamic_light()
		cube_map_lix_t lix(lightning_smgr.add_obj(smap_obj_id, 1));
		lix.add_cube_face_lights(litning_pos, radius, LITN_C, (DX_VAL + DY_VAL), 1); // outdoor_shadows=1
	}
	draw();
}


void lightning_t::gen_recur(point const &start, vector3d const &start_dir, float strength, unsigned parent_len) {

	bool hit_water(0), full_path(1), has_child(0);
	unsigned const path_id((unsigned)paths.size());
	if (path_id >= MAX_LITN_FORKS) return;
	if (path_id > 0 && !is_over_mesh(start)) return; // bad fork - skip
	paths.push_back(lseg_t());
	vector<point> points;
	unsigned const max_steps(MESH_X_SIZE + MESH_Y_SIZE);
	float const step_sz(2.0*HALF_DXY);
	vector3d delta(step_sz*start_dir);
	point pos(start);

	for (unsigned step = 0; step < max_steps; ++step) {
		if (step > 1 && (rgen.rand()%PATH_FORK_MOD) == 0) { // create a fork
			gen_recur(pos, delta.get_norm(), FORK_ATTEN*strength, (parent_len + step));
			has_child = 1;
		}
		float depth(0.0);
		
		if (is_underwater(pos, 0, &depth)) { // terminate path - hit water
			pos.z += depth;
			points.push_back(pos);
			hit_water = 1;
			break;
		}
		int xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

		if (point_outside_mesh(xpos, ypos)) { // outside mesh
			if (path_id == 0) {xpos = max(0, min(MESH_X_SIZE-1, xpos)); ypos = max(0, min(MESH_Y_SIZE-1, ypos));} // primary path - apply clamp
			else {full_path = 0; break;} // terminate path
		}
		float const zval(get_lit_h(xpos, ypos));

		if (pos.z <= zval) { // terminate path - hit something
			if (!points.empty() && (points.back().z - zval) > 0.5f*step_sz) { // hit top of an object, add a line segment down
				pos.z = zval;
				points.push_back(pos);
			} // else hit side of an object, end path here
			break;
		}
		points.push_back(pos);
		if (path_id > 0 && (rgen.rand()%PATH_END_MOD) == 0) {full_path = 0; break;} // end path early

		for (unsigned n = 0; n < 20; ++n) { // make several attempts to find a delta that keeps the line inside the mesh
			vector3d new_delta(delta);
			new_delta += rgen.signed_rand_vector_spherical(STEP_VARIANCE*step_sz); // add in random dir change
			float const mag(new_delta.mag());
			if (mag < TOLERANCE) continue;
			new_delta *= step_sz*rgen.rand_uniform(0.5, 1.0)/mag; // recompute length
			if (dot_product(new_delta, delta) < 0.0f) {new_delta = -new_delta;} // too sharp a turn, switch direction
			if (n+1 == 20) {delta = new_delta; break;} // last attempt, keep it

			if (is_over_mesh(pos)) { // starts inside the scene, so keep it inside the scene
				if (is_over_mesh(pos + new_delta)) {delta = new_delta; break;} // inside the mesh, keep it
			}
			else {
				point pos_clamp(pos);
				pos_clamp.x = max(-X_SCENE_SIZE, min(X_SCENE_SIZE, pos_clamp.x)); // clamp to scene x bounds
				pos_clamp.y = max(-Y_SCENE_SIZE, min(Y_SCENE_SIZE, pos_clamp.y)); // clamp to scene y bounds
				vector3d const target_dir(pos_clamp - pos);
				if (dot_product(new_delta, target_dir) > 0.0f) {delta = new_delta; break;} // moving toward interior of scene, keep it
			}
		} // for n
		delta.z = -fabs(delta.z); // must be negative
		pos += delta;
	} // for step
	assert(path_id < paths.size());
	if (full_path) {points.push_back(pos);} // add lightning path end point
	if (points.size() < 2) {paths.pop_back(); return;} // small segment, discard
	float const dist_ratio(1.0/distance_to_camera(points.back()));
	paths[path_id].points = points;
	paths[path_id].width  = max(rgen.rand_uniform(1.0, 2.0), min(4.0f, L_STRENGTH_MULT*strength*dist_ratio));
	paths[path_id].set_params(parent_len, L_DAMAGE_MULT*strength, full_path, hit_water, has_child);
}


void do_lightning_damage(point &pos, float damage, int hit_water) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos))          return; // object off edge
	if (pos.z - mesh_height[ypos][xpos] > 0.05f) return; // didn't strike the mesh
	
	if (!DISABLE_WATER && wminside[ypos][xpos] == 1) { // interior water
		int const wsi(watershed_matrix[ypos][xpos].wsi);
		float ice_eff((temperature <= W_FREEZE_POINT) ? LITNING_ICE_EFF : 1.0);
		valleys[wsi].fvol += (MEL_EVAP_RATIO/SNOW_ACC)*ice_eff*accumulation_matrix[ypos][xpos]; // melt snow
		accumulation_matrix[ypos][xpos] = 0.0;

		if (hit_water) {
			draw_splash(pos.x, pos.y, pos.z, 0.05);
			valleys[wsi].w_volume = max(0.0f, (valleys[wsi].w_volume - ice_eff*EVAP_AMOUNT));
		}
		else {
			surface_damage[ypos][xpos] += damage;
		}
	}
	else if (!DISABLE_WATER && wminside[ypos][xpos] == 2) { // exterior water
		if (hit_water) {draw_splash(pos.x, pos.y, pos.z, 0.05);}
	}
	else {
		surface_damage[ypos][xpos] += damage;
		modify_grass_at(pos, HALF_DXY, 0, 1); // burn
	}
	if (ypos > 0)             {surface_damage[ypos-1][xpos  ] += damage;}
	if (xpos > 0)             {surface_damage[ypos  ][xpos-1] += damage;}
	if (xpos > 0 && ypos > 0) {surface_damage[ypos-1][xpos-1] += damage;}
	add_crater_to_landscape_texture(pos.x, pos.y, 0.5*damage);
}


void lightning_t::draw() const {

	if (!is_enabled()) return;
	enable_blend();
	set_additive_blend_mode();
	glDepthMask(GL_FALSE); // no depth writing
	shader_t s;
	s.begin_simple_textured_shader();

	for (auto i = paths.begin(); i != paths.end(); ++i) {
		assert(i->points.size() >= 2);
		if (animate2) {add_dynamic_light(0.6*i->width*LITNING_LINEAR_I, i->points.back(), LITN_C);} // unshadowed
		i->draw_lines(1, (i+1 != paths.end())); // fade_ends=1
	}
	s.end_shader();
	glDepthMask(GL_TRUE);
	set_std_blend_mode();
	disable_blend();
	if (animate2 && !ENABLE_LT_SHADOWS) {add_dynamic_light(7.8*LITNING_LINEAR_I, litning_pos, LITN_C);}
}

