// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 11/2/03

#include "gameplay.h"
#include "transform_obj.h"
#include "player_state.h"
#include "physics_objects.h"
#include "textures_3dw.h"


bool const NO_SMILEY_ACTION       = 0;
bool const SMILEYS_LOOK_AT_TARGET = 1;
bool const LEAD_SHOTS             = 0; // doesn't seem to make too much difference (or doesn't work correctly)
bool const SMILEY_BUTT            = 1;
bool const WAYPTS_ALWAYS_VIS      = 0;
unsigned const SMILEY_COLL_STEPS  = 10;
unsigned const PMAP_SIZE          = (2*N_SPHERE_DIV)/3;
int const WP_RESET_FRAMES         = 100; // Note: in frames, not ticks, fix?
int const WP_RECENT_FRAMES        = 200;


float SSTEPS_PER_FRAME(0.0), smiley_speed(1.0), smiley_acc(0);
vector<point> app_spots;
vector<point> waypoints;


extern int island, iticks, num_smileys, free_for_all, teams, frame_counter;
extern int DISABLE_WATER, xoff, yoff, world_mode, spectate, camera_reset, camera_mode, following, game_mode;
extern int recreated, mesh_scale_change, UNLIMITED_WEAPONS;
extern float fticks, temperature, zmax, ztop, XY_SCENE_SIZE, ball_velocity, TIMESTEP, self_damage;
extern point ocean, orig_camera, orig_cdir;
extern int coll_id[];
extern obj_group obj_groups[];
extern obj_type object_types[];
extern player_state *sstates;
extern team_info *teaminfo;
extern vector<string> avail_smiley_names;


// ********** waypt_used_set **********


void waypt_used_set::clear() {

	used.clear();
	last_wp    = 0;
	last_frame = 0;
}

void waypt_used_set::insert(unsigned wp) { // called when a waypoint has been reached

	last_wp    = wp;
	last_frame = frame_counter;
	used[wp]   = frame_counter;
}

bool waypt_used_set::is_valid(unsigned wp) { // called to determine whether or not a waypoint is valid based on when it was last used

	if (last_frame > 0) {
		if ((frame_counter - last_frame) >= WP_RECENT_FRAMES) {
			clear(); // last_wp has expired, so all of the others must have expired as well
			return 1;
		}
		if (wp == last_wp) return 0; // too recent (lasts used)
	}
	map<unsigned, int>::iterator it(used.find(wp));
	if (it == used.end()) return 1; // new waypoint
	if ((frame_counter - it->second) < WP_RESET_FRAMES) return 0; // too recent
	used.erase(it); // lazy update - remove the waypoint when found to be expired
	return 1;
}


// ********** unreachable_pts **********


bool unreachable_pts::cant_reach(point const &pos) const {

	float const toler(SMALL_NUMBER*SMALL_NUMBER);

	for (unsigned j = 0; j < cant_get.size(); ++j) {
		if (p2p_dist_sq(pos, cant_get[j]) < toler) return 1;
	}
	return 0;
}


bool unreachable_pts::proc_target(point const &pos, point const &target, point const &last_target, bool can_reach) {

	float const toler(SMALL_NUMBER*SMALL_NUMBER);

	if (p2p_dist_sq(target, last_target) < toler) { // same as last target
		float const dist_sq(p2p_dist_sq(pos, target)); // current distance

		if (try_dist_sq > 0.0 && dist_sq >= try_dist_sq) { // no closer than before
			try_counts += max(1, iticks);

			if (try_counts > int(1.0*TICKS_PER_SECOND)) { // timeout after 1s - can't get object
				if (!can_reach) cant_get.push_back(target); // give up, can't reach it
				return 0;
			}
		}
		else { // closer
			reset_try(dist_sq);
		}
	}
	else { // reset - new objective
		reset_try();
	}
	return 1;
}


void unreachable_pts::shift_by(vector3d const &vd) {

	for (unsigned i = 0; i < cant_get.size(); ++i) {
		cant_get[i] += vd;
	}
}


// ********** destination_marker **********


 // (x1,y1) is cur pos, (x2,y2) is candidate pos
bool destination_marker::add_candidate(int x1, int y1, int x2, int y2, float depth, float radius) {

	if ((x1 == x2 && y1 == y2) || point_outside_mesh(x2, y2)) return 0;
	
	if (depth > radius) { // still underwater
		if ((valid && depth >= min_depth) || dmin_sq > 0) return 0; // deeper or already found land
		min_depth = (valid ? min(min_depth, depth) : depth);
		valid = 1;
		return 1;
	}
	int const dist_sq((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)); // must be nonzero
	if (valid && dmin_sq != 0 && dist_sq >= dmin_sq) return 0;
	dmin_sq = dist_sq;
	xpos    = x2;
	ypos    = y2;
	valid   = 1;
	return 1;
}


point destination_marker::get_pos() const {

	assert(!point_outside_mesh(xpos, ypos));
	return point(get_xval(xpos), get_yval(ypos), mesh_height[ypos][xpos]);
}


// ********** SMILEY AI CODE **********


bool proj_coll_test(point const &pos, point const &target_pos, vector3d const &orient,
	float target_dist, float radius, int weapon, int coll_id)
{
	int xpos(0), ypos(0), index(0);
	point coll_pos;
	int const test_alpha((weapon == W_LASER) ? 1 : 3);
	point const pos2(target_pos - orient*(1.2*radius));
	if (!coll_pt_vis_test(pos, pos2, 1.2*radius, index, coll_id, 0, test_alpha)) return 0; // cobj collision
		
	if (get_range_to_mesh(pos, orient, coll_pos, xpos, ypos)) {
		if (p2p_dist(pos, coll_pos) + 2.0*radius < target_dist) return 0; // mesh collision
	}
	return 1;
}


void smiley_fire_weapon(int smiley_id) {

	if (!game_mode) return;
	int cid(coll_id[SMILEY]), status;
	float target_dist;
	player_state &sstate(sstates[smiley_id]);
	dwobject const &smiley(obj_groups[cid].get_obj(smiley_id));
	int const &weapon(sstate.weapon);
	if (smiley.disabled()) return;
	if (sstate.target_visible != 1 && (weapon != W_LANDMINE || (rand()&3) != 0)) return;
	int const last_weapon(weapon);
	
	if (weapon == W_UNARMED || (!UNLIMITED_WEAPONS && sstate.no_weap_or_ammo())) {
		init_smiley_weapon(smiley_id); // out of ammo, switch weapons
		if (weapon != last_weapon) sstate.fire_frame = 0;
	}
	if (sstate.target_visible && self_damage > 0.0 && sstate.powerup != PU_SHIELD && weapons[weapon].self_damage) { // can damage self
		if (dist_less_than(sstate.target_pos, smiley.pos, (weapons[weapon].blast_radius + object_types[SMILEY].radius))) { // will damage self
			init_smiley_weapon(smiley_id); // switch weapons to avoid suicide
			if (weapon != last_weapon) sstate.fire_frame = 0;
		}
	}
	assert(!sstate.no_weap_or_ammo());
	if (weapon == W_BALL && (rand()&15) != 0) return; // wait to throw
	point pos(smiley.pos);
	if (temperature <= W_FREEZE_POINT && is_underwater(pos)) return; // under ice
	float const radius(object_types[SMILEY].radius);
	pos.z      += 0.1*radius; // shoot slightly higher than center
	vector3d orient(sstate.target_pos, pos);
	target_dist = orient.mag();

	if (LEAD_SHOTS && sstate.target_visible == 1 && smiley_acc > 0.0) { // should use smiley_acc, target_dist is not quite right
		float const vweap(weapons[weapon].v_add + ball_velocity*weapons[weapon].v_mult);

		if (vweap > TOLERANCE) { // not an instant hit weapon
			int const targ(sstate.target);
			assert(targ >= CAMERA_ID && targ < num_smileys);
			float const oz(orient.z);
			orient   = lead_target(pos, sstate.target_pos, sstate.velocity, sstates[targ].velocity, vweap);
			orient.z = oz; // ???
			orient.normalize();
		}
	}
	if (smiley_acc <= 0.0) {
		orient = smiley.orientation;
	}
	else if (smiley_acc < 1.0) { // add firing error (should this be before or after the range test?)
		vadd_rand(orient, 0.1*((game_mode == 2) ? 0.5 : 1.0)*(1.0 - smiley_acc)/max(0.2f, min(1.0f, target_dist)));
	}
	orient.normalize();
	
	if (weapon != W_LANDMINE && weapon != W_BBBAT && target_dist > 2.0*radius) {
		// make sure it has a clear shot (excluding invisible smileys)
		if (!proj_coll_test(pos, sstate.target_pos, orient, target_dist, radius, weapon, smiley.coll_id)) return;

		if (weapon == W_ROCKET || weapon == W_SEEK_D || weapon == W_PLASMA) { // large projectile
			assert(weapons[weapon].obj_id != UNDEF);
			float const proj_radius(object_types[weapons[weapon].obj_id].radius);
			point const pos2(pos + point(0.0, 0.0, -proj_radius));

			if (!proj_coll_test(pos2, sstate.target_pos, orient, target_dist, radius, weapon, smiley.coll_id)) {
				orient   *= target_dist;
				orient.z += min(proj_radius, 0.7f*radius); // shoot slightly upward
				orient.normalize();
			}
		}
	}
	int &ammo(sstate.p_ammo[weapon]);
	assert(ammo >= 0);
	if (weapon == W_GRENADE && (sstate.wmode&1) && ammo < 3) sstate.wmode = 0;
	int chosen;
	status = fire_projectile(pos, orient, smiley_id, chosen);

	if (status != 0 && !UNLIMITED_WEAPONS && !sstate.no_weap_or_ammo() && weapons[weapon].need_ammo) {
		ammo -= (weapon == W_GRENADE && (sstate.wmode&1)) ? 3 : 1; // check for cluster grenade
		assert(ammo >= 0);
		if (ammo == 0) sstate.fire_frame = 0; // could switch weapons
	}
}


void add_target(vector<od_data> &oddatav, pos_dir_up const &pdu, point const &pos2, float radius, int id, int hitter, int killer) {

	if (!sphere_in_view(pdu, pos2, radius, 0)) return; // view culling
	float dist_sq(p2p_dist_sq(pdu.pos, pos2));
	if (hitter == id) dist_sq *= 0.25; // perfer to attack your attacker
	if (killer == id) dist_sq *= 0.5;  // perfer to attack your last killer in revenge
	oddatav.push_back(od_data(((id == CAMERA_ID) ? CAMERA : SMILEY), id, dist_sq));
}


int find_nearest_enemy(point const &pos, point const &avoid_dir, int smiley_id, point &target, int &target_visible, float &min_dist) {

	assert(smiley_id < num_smileys);
	int min_i(NO_SOURCE), hitter(NO_SOURCE);
	float const radius(object_types[SMILEY].radius);
	int const cid(coll_id[SMILEY]), killer(sstates[smiley_id].killer);
	float tterm, sterm;
	point const camera(get_camera_pos());
	static vector<od_data> oddatav;
	if (sstates[smiley_id].was_hit) hitter = sstates[smiley_id].hitter;
	vector3d const sorient(obj_groups[cid].get_obj(smiley_id).orientation);
	calc_view_test_terms(tterm, sterm, 0);
	pos_dir_up const pdu(pos, sorient, plus_z, tterm, sterm, NEAR_CLIP, FAR_CLIP);
	min_dist = 0.0;

	if (free_for_all) { // smileys attack each other, not only the player
		assert((int)obj_groups[cid].max_objects() == num_smileys);

		for (int i = obj_groups[cid].max_objects()-1; i >= 0; --i) {
			dwobject const &obj(obj_groups[cid].get_obj(i));
			if (obj.disabled() || i == smiley_id || (obj.flags & IN_DARKNESS)) continue;
			if (sstates[i].powerup == PU_INVISIBILITY && hitter != i)          continue; // invisible
			if (same_team(smiley_id, i)) continue; // don't shoot a teammate
			add_target(oddatav, pdu, obj.pos, radius, i, hitter, killer);
		}
	}
	if (camera_mode != 0 && !spectate && !same_team(smiley_id, CAMERA_ID) && (sstates[CAMERA_ID].powerup != PU_INVISIBILITY || hitter == CAMERA_ID)) {
		add_target(oddatav, pdu, camera, radius, CAMERA_ID, hitter, killer); // camera IN_DARKNESS?
	}
	sort(oddatav.begin(), oddatav.end());

	for (unsigned i = 0; i < oddatav.size(); ++i) { // find closest visible target
		point const pos2(get_sstate_pos(oddatav[i].id));
		if (island && (pos2.z < ocean.z)) continue;
		if (avoid_dir != zero_vector && dot_product_ptv(pos2, pos, avoid_dir) > 0.0) continue; // need to avoid this direction
		float const dist(oddatav[i].dist);

		if (sphere_in_view(pdu, pos2, radius, 5)) {
			min_dist = sqrt(dist);
			min_i    = oddatav[i].id;
			assert(min_i >= CAMERA_ID);
			target         = pos2;
			target_visible = 1;
			break;
		}
	}
	oddatav.clear();
	return min_i;
}


struct type_wt_t {
	unsigned type;
	float weight;
	type_wt_t(unsigned t=0, float w=1.0) : type(t), weight(w) {}
};


int find_nearest_obj(point const &pos, point const &avoid_dir, int smiley_id, point &target, float &min_dist, vector<type_wt_t> types) { // health, powerup, etc.

	assert(smiley_id < num_smileys);
	int min_ic(-1);
	float sradius(object_types[SMILEY].radius), sr_sq(sradius*sradius), ra_smiley(C_STEP_HEIGHT*sradius);
	float tterm, sterm;
	vector3d const sorient(obj_groups[coll_id[SMILEY]].get_obj(smiley_id).orientation);
	player_state &sstate(sstates[smiley_id]);
	static vector<od_data> oddatav;
	calc_view_test_terms(tterm, sterm, 0);
	pos_dir_up const pdu(pos, sorient, plus_z, tterm, sterm, NEAR_CLIP, FAR_CLIP);
	min_dist = 0.0;

	for (unsigned i = 0; i < waypoints.size(); ++i) { // inefficient - use subdivision?
		point const wp(waypoints[i]);
		if (!is_over_mesh(wp) || !sstate.waypts_used.is_valid(i)) continue;
		
		if (p2p_dist_sq(wp, pos) < sr_sq) { // smiley has reached waypoint
			sstate.waypts_used.insert(i); // insert as the last used waypoint and remove from consideration
			continue;
		}
		if (!WAYPTS_ALWAYS_VIS && !sphere_in_view(pdu, wp, 0.0, 0)) continue; // view culling - more detailed query later
		oddatav.push_back(od_data(WAYPOINT, i, 4.0*p2p_dist_sq(pos, wp))); // add weight of 4.0 to prefer other objects
	}
	for (unsigned t = 0; t < types.size(); ++t) {
		unsigned const type(types[t].type);
		assert(t < NUM_TOT_OBJS);
		if (UNLIMITED_WEAPONS && (type == WEAPON || type == AMMO || type == WA_PACK)) continue;
		obj_group const &objg(obj_groups[coll_id[type]]);
		if (!objg.is_enabled()) continue;
		float const radius(object_types[type].radius), dmult(1.0/types[t].weight);

		for (unsigned i = 0; i < objg.end_id; ++i) {
			dwobject const &obj(objg.get_obj(i));
			if (obj.disabled() || (obj.flags & IN_DARKNESS))                 continue;
			if (!is_over_mesh(obj.pos) || (island && (obj.pos.z < ocean.z))) continue;
			if (!sphere_in_view(pdu, obj.pos, radius, 0))                    continue; // view culling (disable for predef object locations?)
			if (avoid_dir != zero_vector && dot_product_ptv(obj.pos, pos, avoid_dir) > 0.0) continue; // need to avoid this direction
			float cost(1.0);

			for (int j = 0; j < num_smileys; ++j) { // too slow?
				if (j != smiley_id && sstates[j].objective_pos == obj.pos) {
					float const dist_ratio(p2p_dist(pos, obj.pos)/p2p_dist(obj_groups[coll_id[SMILEY]].get_obj(j).pos, obj.pos));
					if (dist_ratio > 1.0) cost += dist_ratio; // icrease the cost since the other smiley will likely get there first
				}
			}
			oddatav.push_back(od_data(type, i, cost*dmult*p2p_dist_sq(pos, obj.pos)));
		}
	}
	if (oddatav.empty()) { // no objects
		sstate.unreachable.reset_try();
		return min_ic;
	}
	sort(oddatav.begin(), oddatav.end());

	for (unsigned i = 0; i < oddatav.size(); ++i) {
		int const type(oddatav[i].type);
		assert(type >= 0 && type < NUM_TOT_OBJS);
		point pos2;
		dwobject const *obj(NULL);

		if (type == WAYPOINT) {
			assert(size_t(oddatav[i].id) < waypoints.size());
			pos2 = waypoints[oddatav[i].id];
		}
		else {
			int const cid(coll_id[type]);
			assert(cid >= 0 && cid < NUM_TOT_OBJS);
			obj  = &obj_groups[cid].get_obj(oddatav[i].id);
			pos2 = obj->pos;
		}
		int const xpos(get_xpos(pos2.x)), ypos(get_ypos(pos2.y));
		if (point_outside_mesh(xpos, ypos))      continue;
		if (sstate.unreachable.cant_reach(pos2)) continue;
		float const oradius(object_types[type].radius);
		bool const ice(temperature <= W_FREEZE_POINT);
		bool const not_too_high(pos2.z < max((mesh_height[ypos][xpos] + sradius + oradius), (pos.z + ra_smiley + oradius - sradius)));
		bool const not_too_high2(pos2.z < interpolate_mesh_zval(pos2.x, pos2.y, oradius, 0, 0)     + ra_smiley + oradius);
		bool const ice_height_ok(ice && wminside[ypos][xpos] && (pos2.z < water_matrix[ypos][xpos] + ra_smiley + oradius));
		bool const in_motion(obj && obj->status == 1 && obj->velocity.mag_sq() > 9.0);
		bool const can_reach(not_too_high || ice_height_ok || in_motion);

		if (not_too_high2) { // almost in reach, ultimate reachability questionable
			if (!sstate.unreachable.proc_target(pos, pos2, sstate.objective_pos, can_reach)) continue;
		}
		if (can_reach || not_too_high2) { // not_too_high2 - may be incorrect
			if ((WAYPTS_ALWAYS_VIS && type == WAYPOINT) || sphere_in_view(pdu, pos2, oradius, ((type == BALL) ? 5 : 4))) {
				min_dist = sqrt(oddatav[i].dist); // find closest reachable/visible object
				min_ic   = type;
				target   = pos2;
				break;
			}
		}
	}
	oddatav.clear();
	return min_ic;
}


int is_good_smiley_pos(int xpos, int ypos) {

	if (xpos < 1 || ypos < 1 || xpos >= MESH_X_SIZE-1 || ypos >= MESH_Y_SIZE-1) return 0;
	float const radius(object_types[SMILEY].radius), zval(mesh_height[ypos][xpos] + radius);
	if (temperature > W_FREEZE_POINT && zval <= water_matrix[ypos][xpos])       return 0;
	if (island && zval <= ocean.z) return 0;
	int cindex;
	if (!check_legal_move(xpos, ypos, zval, radius, cindex)) return 0;
	return 1;
}


int check_smiley_status(dwobject &obj, int smiley_id) {

	assert(sstates != NULL);
	if (obj.disabled()) return 0;

	if (obj.time > object_types[SMILEY].lifetime) {
		if (sstates) drop_pack(sstates[smiley_id], obj.pos);
		obj.status = 0; // smiley dies of old age (or boredom) after a REALLY long time
		string const msg(sstates[smiley_id].name + " died of " + ((rand()&1) ? "boredom" : "old age"));
		print_text_onscreen(msg, YELLOW, 1.0, MESSAGE_TIME, 0);
		return 0;
	}
	if (obj.time < 1) { // smiley is born
		if (app_spots.empty() && !is_good_smiley_pos(get_xpos(obj.pos.x), get_ypos(obj.pos.y))) {
			obj.status = 0; // shouldn't occur any more
			return 0;
		}
		init_smiley(smiley_id);
	}
	return 1;
}


void smiley_select_target(dwobject &obj, int smiley_id) {

	if (obj.disabled()) return;
	int min_ie(NO_SOURCE), min_ih(-1);
	float diste, disth, health(obj.health);
	vector<type_wt_t> types;
	point targete(all_zeros), targeth(all_zeros);
	player_state &sstate(sstates[smiley_id]);
	int const last_target_visible(sstate.target_visible), last_target_type(sstate.target_type);
	sstate.target_visible = 0;
	sstate.target_type    = 0; // 0 = none, 1 = enemy, 2 = health/powerup
	float const health_eq(min(4.0f*health, (health + sstate.shields)));

	// look for landmines and avoid them
	point avoid_dir(zero_vector);
	obj_group const &objg(obj_groups[coll_id[LANDMINE]]);
	
	if (objg.is_enabled()) {
		float min_dist(weapons[LANDMINE].blast_radius);

		for (unsigned i = 0; i < objg.end_id; ++i) {
			dwobject const &obj2(objg.get_obj(i));
			
			if (!obj2.disabled() && !lm_coll_invalid(obj2) && obj2.source == smiley_id && dist_less_than(obj.pos, obj2.pos, min_dist)) {
				avoid_dir  = (obj2.pos - obj.pos);
				min_dist   = avoid_dir.mag();
				avoid_dir /= min_dist;
			}
		}
	}
	if (game_mode == 2 && !UNLIMITED_WEAPONS) { // want the ball
		if (sstate.p_ammo[W_BALL] > 0) { // already have a ball
			min_ie = find_nearest_enemy(obj.pos, avoid_dir, smiley_id, sstate.target_pos, sstate.target_visible, diste);
		}
		if (!sstate.target_visible) { // don't have a ball or no enemy in sight
			types.push_back(type_wt_t(BALL, 1.0));
			min_ih = find_nearest_obj(obj.pos, avoid_dir, smiley_id, sstate.target_pos, disth, types);
		}
	}
	else if (health_eq < 20.0) { // want health (or shields)
		types.push_back(type_wt_t(HEALTH, 1.5));
		types.push_back(type_wt_t(SHIELD, (1.0 - sstate.shields/MAX_SHIELDS)));
		min_ih = find_nearest_obj(obj.pos, avoid_dir, smiley_id, sstate.target_pos, disth, types);
		if (min_ih < 0) min_ie = find_nearest_enemy(obj.pos, avoid_dir, smiley_id, sstate.target_pos, sstate.target_visible, diste);
	}
	else { // choose health/attack based on distance
		// want to get powerup badly if you don't have one
		float const pu_wt((sstates[smiley_id].powerup == 0) ? 1.5 : (1.0 - float(sstate.powerup_time)/POWERUP_TIME));
		types.push_back(type_wt_t(POWERUP, pu_wt));
		types.push_back(type_wt_t(WEAPON,  0.8));
		types.push_back(type_wt_t(AMMO,    0.7));
		types.push_back(type_wt_t(WA_PACK, 1.0));
		types.push_back(type_wt_t(SHIELD,  1.2*(1.0 - sstate.shields/MAX_SHIELDS))); // always below max since it ticks down over time
		if (health < MAX_HEALTH) types.push_back(type_wt_t(HEALTH, 1.5*(1.0 - health/MAX_HEALTH)));
		min_ie = find_nearest_enemy(obj.pos, avoid_dir, smiley_id, targete, sstate.target_visible, diste);
		min_ih = find_nearest_obj(  obj.pos, avoid_dir, smiley_id, targeth, disth, types);

		if (!sstate.target_visible) { // can't find an enemy, choose health/pickup
			if (min_ih >= 0) sstate.target_pos = targeth;
		}
		else if (min_ih < 0) { // can't find health/pickup, choose attack
			if (sstate.target_visible) sstate.target_pos = targete;
		}
		else { // found both item and enemy
			float const dp(dot_product((obj.pos - targeth).get_norm(), (obj.pos - targete).get_norm()));

			if (diste <= disth || dp > 0.95) { // enemy is closer or in the direction of the desired pickup
				sstate.target_pos = targete;
				min_ih            = -1;
			}
			else { // health/pickup is closer
				sstate.target_pos     = targeth;
				sstate.target_visible = 0;
			}
		}
	}
	if (sstate.target_visible) { // targeting and enemy
		sstate.target_type   = 1;
		sstate.objective_pos = sstate.target_pos;
	}
	else if (sstates[smiley_id].was_hit && sstates[smiley_id].hit_dir != zero_vector) { // turn around if hit in the back
		sstate.target_pos    = obj.pos + sstates[smiley_id].hit_dir.get_norm();
		sstate.target_type   = 1; // enemy
		sstate.objective_pos = sstate.target_pos;
	}
	else {
		sstate.objective_pos = sstate.target_pos;
		
		if (min_ih >= 0) {
			sstate.target_type    = 2;
			sstate.target_visible = 2;
		}
		if (sstate.target_pos != obj.pos) { // look beyond the target
			vector3d const o(sstate.target_pos - obj.pos); // motion direction
			sstate.target_pos += o*(2.0*object_types[SMILEY].radius/o.mag());
		}
	}
#if 0 // need proper path finding for this to work correctly
	if (!sstate.target_visible && last_target_visible == 1) { // no targets - use last known target enemy location
		sstate.target_visible = last_target_visible;
		sstate.target_type    = last_target_type;
		sstate.objective_pos  = sstate.target_pos; // should still be valid
		min_ie = sstate.target;
	}
#endif
	if (sstate.target_visible == 1) assert(min_ie >= CAMERA_ID);
	sstate.target = min_ie;
}


vector3d get_orient_from_mesh(int x1, int y1, int x2, int y2) {

	assert(x1 != x2 || y1 != y2);
	float zval(0.0);
	if (!is_mesh_disabled(x1, y1) && !is_mesh_disabled(x2, y2)) zval = (mesh_height[y2][x2] - mesh_height[y1][x1]);
	return vector3d((x2 - x1)*DX_VAL, (y2 - y1)*DY_VAL, zval).get_norm();
}


int smiley_motion(dwobject &obj, int smiley_id) {

	if (NO_SMILEY_ACTION || obj.disabled()) return 0;
	player_state &sstate(sstates[smiley_id]);
	int step, move, moved(0), dir((int)obj.direction), ddir, cdir(0), chdir, ds_count;
	int chase_target_var(0), target_type(sstate.target_type);
	float const speed(sstate.get_rspeed_scale()), radius(object_types[SMILEY].radius);
	assert(radius >= 0.0);
	float const step_dist(smiley_speed*fticks*GROUND_SPEED), step_height(S_SH_SCALE*C_STEP_HEIGHT*radius);
	point const opos(obj.pos);
	point target(sstate.target_pos);
	vector3d v(0.0, 0.0, 0.0);
	int movements[9][2], dir_status[9];
	if (SSTEPS_PER_FRAME == 0.0) SSTEPS_PER_FRAME = (SSTEPS_PER_FRAME0/XY_SCENE_SIZE)*((float)XY_SUM_SIZE/256.0);
	float const health_eq(min(4.0f*obj.health, (obj.health + sstate.shields)));
	int xpos(get_xpos_clamp(opos.x)), ypos(get_ypos_clamp(opos.y)), xnew(xpos), ynew(ypos);
	obj.angle += speed*min(8.0f, fticks)*SSTEPS_PER_FRAME; // sigma-delta A/D conversion
	int const nsteps(max(1, min(MAX_XY_SIZE, (int)obj.angle)));
	obj.angle -= (float)nsteps;
	bool const has_flight(sstate.powerup == PU_FLIGHT), is_water_temp(temperature > W_FREEZE_POINT);
	int const xcpos(get_xpos(target.x)), ycpos(get_ypos(target.y));
	
	// STOP N S E W NE NW SE SW
	// 0    1 5 3 7 2  8  4  6
	for (step = 0; step < nsteps; ++step) {
		if (sstate.target_visible) {
			if (xcpos == xpos && ycpos == ypos) {
				cdir = 0;
			}
			else {
				point2d<float> tdir(xcpos - xpos, ycpos - ypos);
				tdir.normalize();
				float angle(atan2(tdir.y, tdir.x) + PI);
				int const cmap[9] = {7,6,5,4,3,2,1,8,7}; // can wraparound by 1
				int const val(int(angle/(PI/4) + 0.5));
				assert(val >= 0 && val < 9);
				cdir = cmap[val];
			}
			
			// if there has been bloodshed and the smiley can see and attack the target
			if      (target_type == 0) chase_target_var = 0;
			else if (target_type == 2) chase_target_var = 4 + (150 - min(150, int(health_eq)))/15;
			else if (health_eq < 10.0 && sstate.powerup != PU_DAMAGE) { // run away unless you have quad damage
				if ((rand()&3) == 0) {
					chase_target_var = 1;
					cdir = (cdir+4)%9; // reverse direction
				}
				else chase_target_var = rand()&1;
			}
			else if (health_eq < 15.0)  chase_target_var =     rand()%3; // 0, 1, 2 (1.0)
			else if (health_eq < 40.0)  chase_target_var = 1 + rand()%2; // 1, 2 (1.5)
			else if (health_eq < 100.0) chase_target_var = 1 + rand()%3; // don't care much, 1, 2, 3 (2.0)
			else                        chase_target_var = 2 + rand()%3; // attack, 2, 3, 4 (3.0)
		}
		for (unsigned i = 0; i < 9; ++i) {
			movements[i][0] = xpos;
			movements[i][1] = ypos;
		}
		int const edist(rand() & 1);

		if (ypos < MESH_Y_SIZE-edist-1) { // N
			++movements[8][1];
			++movements[1][1];
			++movements[2][1];
		}
		if (ypos > edist) { // S
			--movements[4][1];
			--movements[5][1];
			--movements[6][1];
		}
		if (xpos < MESH_X_SIZE-edist-1) { // E
			++movements[2][0];
			++movements[3][0];
			++movements[4][0];
		}
		if (xpos > edist) { // W
			--movements[6][0];
			--movements[7][0];
			--movements[8][0];
		}
		ds_count = 9;
		for (unsigned i = 0; i < 9; ++i) dir_status[i] = 1;
		unsigned tries(0);
		
		for (; tries < SMILEY_MAX_TRIES && ds_count > 0; ++tries) {
			dir = (int)obj.direction;
			if (tries == SMILEY_MAX_TRIES/2) dir = (dir + 1)%9; // force slow right turn
			if (dir == 0 && rand()%10 < 5)   dir = rand()%9; // stopped, choose new dir
			assert(dir >= 0 && dir < 9);

			if (ds_count == 1) {
				for (unsigned j = 0; j < 9; ++j) {
					if (dir_status[j] != 0) {dir = j; break;}
				}
			}
			else {
				do {
					bool bias(0);
					int const rval(rand());
					move = rval%100;

					if (move%20 < chase_target_var) { // move towards target
						bias = 1;

						if (dir_status[cdir] == 0) { // navagate around obstacle towards target
							move = (rval>>8)&1;
							dir  = (cdir+2*move-1+9)%9;

							if (dir == 0) dir = (cdir+4*move-2+9)%9;
							if (dir_status[dir] != 0) cdir = dir;
							else {
								dir = (cdir-2*move+1+9)%9;
								if (dir == 0) dir = (cdir-4*move+2+9)%9;
								if (dir_status[dir] != 0) cdir = dir;
								else bias = 0;
							}
						}
						else dir = cdir;
					}
					if (!bias && move < int(100.0*SMILEY_DIR_FACTOR)) {
						ddir  = 0;
						chdir = (rval>>10)%15;

						for (unsigned c = 0; c < 4; ++c) {
							if (chdir >= ((1<<c)-1)) ++ddir;
						}
						dir = (dir + (((rval>>16)&1) ? ddir : (9-ddir)))%9;
					}
				} while (dir_status[dir] == 0);
			}
			assert(dir >= 0 && dir < 9);

			if (dir != 0) {
				ddir = dir - (int)obj.direction;
				dir  = (int)obj.direction;
				assert(dir >= 0 && dir < 9);
				if (abs(ddir) == 4) dir += ((rand()&1) ? 1 : -1);
				else if (ddir != 0) dir += (((ddir > 0 && ddir < 4) || ddir < -4) ? 1 : -1);
				if        (dir < 1) dir += 8;
				else if   (dir > 8) dir -= 8;
			}
			assert(dir >= 0 && dir < 9);
			if (dir_status[dir] == 0) continue;
			xnew = movements[dir][0];
			ynew = movements[dir][1];
			bool bad_move(point_outside_mesh(xnew, ynew));

			if (!bad_move) {
				float const mhyx(mesh_height[ynew][xnew]);

				if ((mhyx - obj.pos.z) > step_height && !has_flight) {
					bad_move = 1; // too high to step
				}
				else if (!has_flight && (is_water_temp &&
					max(obj.pos.z, (mhyx + radius)) < water_matrix[ynew][xnew]) || (island && ((mhyx + radius) <= ocean.z)))
				{
					bad_move = 1; // don't go under water/blood
				}
				else if (xnew != xpos || ynew != ypos) {
					// *** FIXME: check for walking into your own landmine, falling off a cliff, and other stupid decisions ***
					obj.orientation = get_orient_from_mesh(xpos, ypos, xnew, ynew);
					xpos = xnew;
					ypos = ynew;
				}
				else if (rand()%10 < 5) {
					cdir = (cdir+4)%9; // reverse direction
					continue;
				}
			}
			if (bad_move) {
				dir_status[dir] = 0;
				--ds_count;
				continue;
			}
			assert(dir >= 0 && dir < 9);
			obj.direction = (char)dir;
			moved = 1;
			break;
		} // for tries
		int cindex(-1);
		assert(!point_outside_mesh(xpos, ypos));

		if (!moved && (dir_status[0] == 0 || tries == SMILEY_MAX_TRIES-1)) { // stuck - not sure if this helps
			int xpos2, ypos2; // usually takes one iteration, I think
			do {xpos2 = (xpos + ((rand()&3) - 1));} while (xpos2 < 0 || xpos2 >= MESH_X_SIZE);
			do {ypos2 = (ypos + ((rand()&3) - 1));} while (ypos2 < 0 || ypos2 >= MESH_Y_SIZE);
			xpos = xpos2;
			ypos = ypos2;
			break;
		}
	} // for step
	ypos = max((int)0, min(MESH_Y_SIZE-1, ypos));
	xpos = max((int)0, min(MESH_X_SIZE-1, xpos));
	assert(!point_outside_mesh(xnew, ynew));

	if (!moved && (xnew != xpos || ynew != ypos)) {
		obj.orientation = get_orient_from_mesh(xpos, ypos, xnew, ynew);
		assert(dir >= 0 && dir < 9);
		obj.direction = (char)dir;
	}
	if (!DISABLE_WATER && is_water_temp && wminside[ypos][xpos] &&
		water_matrix[ypos][xpos] >= max(obj.pos.z, (mesh_height[ypos][xpos] + radius)))
	{ // stuck - underwater
		sstate.dest_mark.update_dmin(xpos, ypos);

		for (unsigned i = 0; i < SMILEY_MAX_TRIES; ++i) { // find some randomly chosen dry spots to head for
			int const xt(rand() % MESH_X_SIZE), yt(rand() % MESH_Y_SIZE); // rand() is only 16 bits here?
			float const depth(wminside[yt][xt] ? (water_matrix[yt][xt] - mesh_height[yt][xt]) : 0.0);
			sstate.dest_mark.add_candidate(xpos, ypos, xt, yt, depth, radius);
		}
		if (sstate.dest_mark.valid) {
			xpos = sstate.dest_mark.xpos;
			ypos = sstate.dest_mark.ypos;
			sstate.target_visible = 0;
			sstate.target_type    = target_type = 1;
			sstate.objective_pos  = sstate.target_pos = target = sstate.dest_mark.get_pos();
		}
	}
	else {
		sstate.dest_mark.clear();
	}

	// do actual movement
	vector3d stepv((get_xval(xpos) - obj.pos.x), (get_yval(ypos) - obj.pos.y), 0.0);

	if (stepv.normalize_test()) {
		vadd_rand(stepv, 0.1);
		obj.pos += stepv*step_dist;
	}
	obj.velocity.assign(0.0, 0.0, -1.0);
	point const wanted_pos(obj.pos);
	bool stuck(0), no_up(0), no_down(0);
	int const coll(obj.multistep_coll(opos, smiley_id, SMILEY_COLL_STEPS));

	if (has_flight) {
		no_down = (obj.pos.z > wanted_pos.z);
		no_up   = (obj.pos.z < wanted_pos.z);
	}
	else {
		clip_to_scene(obj.pos);
	}
	if (obj.status == 0) { // killed
		remove_reset_coll_obj(obj.coll_id);
		return 0;
	}
	if (coll && p2p_dist(obj.pos, opos) < p2p_dist(obj.pos, wanted_pos)) {
		if (++sstate.stopped_time > 4) stuck = 1;
	}
	else {
		sstate.stopped_time = 0;
	}
	int ohval(set_true_obj_height(obj.pos, opos, C_STEP_HEIGHT, sstate.zvel, SMILEY, smiley_id, has_flight, 0));
	float const zval(obj.pos.z);

	if (recreated || mesh_scale_change || (zval - opos.z) < step_height || (has_flight && !no_down)) { // also gets them stuck in ice
		obj.pos.z = min(zval, (opos.z + (has_flight ? 0.0f : 0.4f*radius)));
	}
	else if (!has_flight) {
		obj.pos = opos; // reset to old position
		ohval   = 3;
	}
	if (stuck || ohval == 3 || (!point_interior_to_mesh(xnew, ynew) && (rand()&3) == 0)) {
		obj.direction = (((int)obj.direction) + 4)%9; // stuck, turn around
	}
	if (has_flight) {
		if ((target_type == 2 && target.z < opos.z && obj.pos.z < opos.z) || // want health/powerup below - go down
			(target_type == 1 && sstate.weapon == W_BBBAT)) // have to go down for baseball bat hit
		{
			obj.pos.z = max(max(target.z, obj.pos.z), (opos.z - (no_down ? 0.0f : float(radius*rand_uniform(0.2, 0.6)))));
		}
		else {
			obj.pos.z = min(min((zmax + 0.5f), (ztop + 1.0f)), max(obj.pos.z, opos.z));
			if (!no_up && (rand()&7) == 0) obj.pos.z += radius*rand_uniform(0.0, 0.5); // sometimes fly up
		}
		if (!point_outside_mesh(xpos, ypos)) {
			obj.pos.z = max(obj.pos.z, (mesh_height[ypos][xpos] + radius)); // just in case (i.e. mesh was changed when under ice)
		}
		obj.pos.z = min(obj.pos.z, (opos.z + step_height)); // limit max rise distance
	}
	else { // ensure smiley stays on the mesh
		for (unsigned i = 0; i < 2; ++i) {
			obj.pos[i] = max(-SCENE_SIZE[i]+radius, min(SCENE_SIZE[i]-radius, obj.pos[i]));
		}
		obj.pos.z = max(obj.pos.z, interpolate_mesh_zval(obj.pos.x, obj.pos.y, 0.0, 1, 1) + radius); // ignore ice
	}
	if (SMILEYS_LOOK_AT_TARGET && target_type != 0) {
		obj.orientation = ((target == obj.pos) ? signed_rand_vector_norm() : (target - obj.pos).get_norm());
	}
	if (obj.pos.x == opos.x && obj.pos.y == opos.y) {
		if (obj.direction > 0 || rand()%4 == 0) obj.direction = (obj.direction + 1)%9; // last resort - force slow right turn
	}
	if (spectate && smiley_id == 0 /*&& !camera_view*/) { // can only follow smiley 0 for now
		obj_groups[coll_id[SMILEY]].get_obj(smiley_id).flags |= CAMERA_VIEW;
		orig_camera   = camera_origin;
		camera_reset  = 0;
		following     = 1;
		camera_mode   = 1;
		orig_cdir     = cview_dir;
		cview_dir     = obj.orientation; // doesn't seem to actually change the viewing direction
		camera_origin = obj.pos;
		update_cpos();
	}
	if (is_water_temp && is_underwater(obj.pos, 1) && (rand()&1)) gen_bubble(obj.pos);
	sstate.velocity = (obj.pos - opos)/(TIMESTEP*fticks);
	sstate.ctv      = chase_target_var;
	obj.status      = 3;
	return 1;
}


void advance_smiley(dwobject &obj, int smiley_id) { // seems to slightly favor smileys with later ids

	assert(smiley_id < num_smileys);
	assert(obj.type == SMILEY);
	assert(obj_groups[coll_id[SMILEY]].enabled);
	if (!check_smiley_status(obj, smiley_id)) {sstates[smiley_id].fall_counter = 0; return;}
	smiley_select_target(obj, smiley_id);
	obj.time += iticks;
	if (!smiley_motion(obj, smiley_id)) {sstates[smiley_id].fall_counter = 0; return;}
	smiley_action(smiley_id);
}


void shift_player_state(vector3d const &vd, int smiley_id) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);
	player_state &sstate(sstates[smiley_id]);
	sstate.target_pos    += vd;
	sstate.objective_pos += vd;
	sstate.unreachable.shift_by(vd);
}


void add_damage_to_smiley_texture(vector3d const &dir, float size, int smiley_id, int type) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);
	assert(size > 0.0);
	colorRGBA color;
	vector3d sdir(obj_groups[coll_id[SMILEY]].get_obj(smiley_id).orientation);
	unsigned char *tdata(sstates[smiley_id].tdata);
	assert(tdata != NULL);
	int const tsize(int(size*SMILEY_TEX_SIZE/100.0 + 0.5)), radsq(4*tsize*tsize);
	int const tex_size(SMILEY_TEX_SIZE*SMILEY_TEX_SIZE);
	int tx, ty;
	get_tex_coord(dir, sdir, SMILEY_TEX_SIZE, SMILEY_TEX_SIZE, tx, ty, 0);
	int x1(tx - tsize), y1(ty - 2*tsize), x2(tx + tsize), y2(ty + 2*tsize);

	switch (type) {
		case IMPACT:
			color = PURPLE; break;
		case PLASMA: case FIRE: case BURNED: case LASER: case BLAST_RADIUS:
			color = BLACK;  break;
		default:
			color = RED;
	}
	unsigned char c[3];
	unpack_color(c, color);

	for (int yy = y1; yy < y2; ++yy) { // allow texture wrap
		int const y((yy + SMILEY_TEX_SIZE) % SMILEY_TEX_SIZE), offset(y*SMILEY_TEX_SIZE), yterm((yy-ty)*(yy-ty));

		for (int xx = x1; xx < x2; ++xx) { // allow texture wrap
			int const x((xx + SMILEY_TEX_SIZE) % SMILEY_TEX_SIZE);

			if (((xx-tx)*(xx-tx) << 2) + yterm <= radsq) {
				double const dist(sqrt(double((xx-tx)*(xx-tx) + yterm)));
				double const blend(1.0/(dist + 2.0));
				assert(offset+x < tex_size);
				int const index(3*(offset+x));
				BLEND_COLOR((tdata + index), c, (tdata + index), blend);
			}
		}
	}
	if (y1 < 0 || y2 > (int)SMILEY_TEX_SIZE) { // texture wraparound - just reload the whole thing for now
		y1 = 0;
		y2 = SMILEY_TEX_SIZE;
	}
	select_smiley_texture(smiley_id); // update a strip from y1 to y2
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, y1, SMILEY_TEX_SIZE, (y2-y1), GL_RGB, GL_UNSIGNED_BYTE, (tdata + 3*SMILEY_TEX_SIZE*y1));
	glDisable(GL_TEXTURE_2D);
}


void add_damage_to_smiley_surface(vector3d const &dir, float size, int smiley_id) {

	float const radius(object_types[SMILEY].radius);
	//obj_groups[coll_id[SMILEY]].get_td()->add_rand_perturb(smiley_id, 0.004*size*radius, -0.6*radius, 0.4*radius);
	transform_data *td(obj_groups[coll_id[SMILEY]].get_td());
	td->set_perturb_size(smiley_id, PMAP_SIZE);
	vector3d sdir(obj_groups[coll_id[SMILEY]].get_obj(smiley_id).orientation);
	int const tsize(int(size*PMAP_SIZE/60.0 + 0.5)), radsq(4*tsize*tsize);
	int tx, ty;
	get_tex_coord(dir, sdir, PMAP_SIZE, PMAP_SIZE, tx, ty, 1);
	int x1(tx - tsize), y1(ty - 2*tsize), x2(tx + tsize), y2(ty + 2*tsize);

	for (int yy = y1; yy < y2; ++yy) { // allow texture wrap
		int const y((yy + PMAP_SIZE) % PMAP_SIZE), yterm((yy-ty)*(yy-ty));

		for (int xx = x1; xx < x2; ++xx) { // allow texture wrap
			int const x((xx + PMAP_SIZE) % PMAP_SIZE);

			if (((xx-tx)*(xx-tx) << 2) + yterm <= radsq) {
				double const dist(sqrt(double((xx-tx)*(xx-tx) + yterm)));
				td->add_perturb_at(x, y, smiley_id, -0.02/(dist + 2.0), -0.5*radius, 0.5*radius);
			}
		}
	}
}


void add_damage_to_smiley(vector3d const &dir, float size, int smiley_id, int type) {

	add_damage_to_smiley_texture(dir, size, smiley_id, type);
	if (size > 4.0) add_damage_to_smiley_surface(dir, size, smiley_id);
}


void init_smiley_texture(int smiley_id) {

	assert(sstates != NULL);
	int const texels(SMILEY_TEX_SIZE*SMILEY_TEX_SIZE), tbytes(3*texels);
	unsigned char *tdata(sstates[smiley_id].tdata);
	if (tdata == NULL) tdata = sstates[smiley_id].tdata = new unsigned char[tbytes];
	memset(tdata, 255, tbytes*sizeof(unsigned char));

	if (SMILEY_BUTT) {
		for (unsigned y = SMILEY_TEX_SIZE/4; y < SMILEY_TEX_SIZE/2; ++y) {
			for (unsigned i = 0; i < 3; ++i) { // butt crack ;)
				tdata[3*(y*SMILEY_TEX_SIZE + (SMILEY_TEX_SIZE/2)) + i] = 0;
			}
		}
	}
	if (sstates[smiley_id].tid == 0) {
		setup_texture(sstates[smiley_id].tid, GL_MODULATE, 0, 0, 0);
	}
	else {
		glBindTexture(GL_TEXTURE_2D, sstates[smiley_id].tid);
	}
	glTexImage2D(GL_TEXTURE_2D, 0, 3, SMILEY_TEX_SIZE, SMILEY_TEX_SIZE, 0, GL_RGB, GL_UNSIGNED_BYTE, tdata);
}


void init_smiley(int smiley_id) {

	assert(sstates != NULL);
	assert(smiley_id >= 0 && smiley_id < num_smileys);
	int const cid(coll_id[SMILEY]);
	assert(cid >= 0 && cid < NUM_TOT_OBJS);
	obj_group &objg(obj_groups[cid]);
	if (!objg.enabled) return;
	init_smiley_texture(smiley_id);
	init_sstate(smiley_id, (game_mode == 1));
	init_smiley_weapon(smiley_id);
	dwobject &obj(objg.get_obj(smiley_id));
	obj.direction     = 0;
	obj.orientation   = signed_rand_vector(); // this is likely reset later
	obj.orientation.z = 0.0;
	obj.orientation.normalize();
	objg.get_td()->reset_perturb_if_set(smiley_id);

	/*if (SMILEY_BUTT) {
		transform_data *td(obj_groups[coll_id[SMILEY]].get_td());
		td->set_perturb_size(smiley_id, PMAP_SIZE);
		float const radius(object_types[SMILEY].radius);
		td->add_perturb_at(PMAP_SIZE/2, 14, smiley_id, -0.3*radius, -0.5*radius, 0.5*radius);
	}*/
}


void init_smiley_weapon(int smiley_id) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);
	int bbat_iter(0), bb_range(-1);
	player_state &sstate(sstates[smiley_id]);
	sstate.wmode = ((rand()&3) == 0);

	if (game_mode == 2) { // dodgeball mode
		sstate.weapon     = ((UNLIMITED_WEAPONS || sstate.p_ammo[W_BALL] > 0) ? W_BALL : W_UNARMED);
		sstate.fire_frame = 0;
		return;
	}
	do {
		sstate.weapon = 1 + (rand()%(NUM_WEAPONS-1)); // no cgrenade

		if (sstate.weapon == W_BBBAT && bbat_iter++ < 10 && (rand()%8) != 0) { // only use baseball bat when in close range
			if (bb_range == -1) {
				float dist(FAR_CLIP);
				if (sstate.target_visible) dist = p2p_dist(sstate.target_pos, obj_groups[coll_id[SMILEY]].get_obj(smiley_id).pos);
				bb_range = (dist < 6.0*object_types[SMILEY].radius);
			}
			if (!bb_range) sstate.weapon = W_UNARMED;
		}
		else if (sstate.weapon == W_SBALL && (rand()%5) != 0) {
			sstate.weapon = W_UNARMED; // prefers not to use a bouncy ball
		}
	} while (sstate.weapon == W_UNARMED || sstate.no_weap_or_ammo());

	//sstate.weapon     = W_GRENADE; // set to force this as weapon choice (if available)
	sstate.fire_frame = 0;
	if (sstate.weapon == W_PLASMA && sstate.wmode == 1 && (rand()%4) != 0) sstate.plasma_loaded = 1; // fire it up!
}


void smiley_action(int smiley_id) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);
	float depth(0.0);
	player_state &sstate(sstates[smiley_id]);
	int const weapon(sstate.weapon), chase(sstate.target_visible && sstate.ctv);
	dwobject &smiley(obj_groups[coll_id[SMILEY]].get_obj(smiley_id));

	if (chase && rand()%10 != 0) {
		float const range(weapons[weapon].range);
		if (range == 0.0 || dist_less_than(sstate.target_pos, smiley.pos, range)) smiley_fire_weapon(smiley_id);
	}
	if (sstate.powerup == PU_REGEN) smiley.health = min(MAX_REGEN_HEALTH, smiley.health + 0.1f*fticks);
	if (rand()%500 == 0) init_smiley_weapon(smiley_id); // change weapons
	if (sstate.was_hit > 0) --sstate.was_hit;
	++sstate.kill_time;
	check_underwater(smiley_id, depth);
}


colorRGBA get_smiley_team_color(int smiley_id) {

	colorRGBA const tcolors[] = {RED, BLUE, GREEN, MAGENTA, CYAN, BROWN, PINK, GRAY, LT_BLUE, ORANGE, WHITE};
	if (teams <= 1) return tcolors[0];
	return tcolors[min(((smiley_id+1)%teams), int(sizeof(tcolors)/sizeof(colorRGBA)-1))];
}


int gen_smiley_or_player_pos(point &pos, int index) {

	if (!app_spots.empty()) {
		float const dmin(2.0*object_types[SMILEY].radius);
		point pos0(pos);

		for (unsigned i = 0; i < SMILEY_MAX_TRIES; ++i) {
			pos0 = app_spots[rand()%app_spots.size()];
			//if (!is_good_smiley_pos(get_xpos(pos0.z), get_ypos(pos0.y))) continue; // bad starting pos
			if (camera_mode == 1 && !spectate && dist_less_than(pos0, get_camera_pos(), dmin)) continue;
			obj_group const &objg(obj_groups[coll_id[SMILEY]]);
			bool too_close(0);

			if (objg.is_enabled()) {
				for (unsigned i = 0; i < objg.end_id && !too_close; ++i) {
					if (!objg.get_obj(i).disabled() && dist_less_than(pos0, objg.get_obj(i).pos, dmin)) too_close = 1;
				}
			}
			if (!too_close) break;
		}
		pos = pos0;
		//return 0;
		return 1; // bad app spot, but return it anyway
	}
	int const team(teams ? (index + teams + 1)%teams : 0);
	float const xc(teaminfo[team].bb.x1 + teaminfo[team].bb.x2);
	float const yc(teaminfo[team].bb.y1 + teaminfo[team].bb.y2);
	float const xd(teaminfo[team].bb.x2 - teaminfo[team].bb.x1);
	float const yd(teaminfo[team].bb.y2 - teaminfo[team].bb.y1);
	unsigned const MAX_ITER(100);
	pos.assign(0.0, 0.0, 0.0);
	
	for (unsigned i = 0; i < MAX_ITER; ++i) {
		pos.x = 0.5*(xc + xd*signed_rand_float());
		pos.y = 0.5*(yc + yd*signed_rand_float());
		
		if (is_good_smiley_pos(get_xpos(pos.x), get_ypos(pos.y))) {
			pos.z = interpolate_mesh_zval(pos.x, pos.y, object_types[SMILEY].radius, 0, 0);
			return 1;
		}
	}
	return 0; // bad pos
}


void select_smiley_texture(int smiley_id) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);

	if (!glIsTexture(sstates[smiley_id].tid)) {
		free_texture(sstates[smiley_id].tid);
		init_smiley_texture(smiley_id);
	}
	assert(glIsTexture(sstates[smiley_id].tid));
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, sstates[smiley_id].tid);
}


void free_smiley_textures() {

	if (sstates == NULL) return;

	for (int i = 0; i < num_smileys; ++i) {
		free_texture(sstates[i].tid);
	}
}


