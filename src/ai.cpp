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
bool const LEAD_SHOTS             = 1; // doesn't seem to make too much difference (or doesn't work correctly)
bool const ACCURATE_GRAV_PREDICT  = 0;
bool const SMILEY_BUTT            = 1;
int const WAYPT_VIS_LEVEL[2]      = {1, 2}; // {unconnected, connected}: 0: visible in frustum, 1 = visible from the position, 2 = always visible
unsigned const SMILEY_COLL_STEPS  = 10;
unsigned const PMAP_SIZE          = (2*N_SPHERE_DIV)/3;
float const UNREACHABLE_TIME      = 0.75; // in seconds


float smiley_speed(1.0), smiley_acc(0);
vector<point> app_spots;


extern bool has_wpt_goal;
extern int island, iticks, num_smileys, free_for_all, teams, frame_counter, display_mode;
extern int DISABLE_WATER, xoff, yoff, world_mode, spectate, camera_reset, camera_mode, following, game_mode;
extern int recreated, mesh_scale_change, UNLIMITED_WEAPONS, camera_coll_id;
extern float fticks, tfticks, temperature, zmax, ztop, XY_SCENE_SIZE, TIMESTEP, self_damage, base_gravity;
extern double camera_zh;
extern point ocean, orig_camera, orig_cdir;
extern int coll_id[];
extern obj_group obj_groups[];
extern obj_type object_types[];
extern player_state *sstates;
extern team_info *teaminfo;
extern vector<string> avail_smiley_names;
extern waypoint_vector waypoints;



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

			if (try_counts > int(UNREACHABLE_TIME*TICKS_PER_SECOND)) { // timeout - can't get object
				if (!can_reach) add(target); // give up, can't reach it
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
	depth = max(0.0f, depth); // can't have negative depth
	
	if (depth > radius) { // still underwater
		if ((valid && depth >= min_depth) || dmin_sq > 0) return 0; // deeper or already found land
		min_depth = (valid ? min(min_depth, depth) : depth);
	}
	else {
		int const dist_sq((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)); // must be nonzero
		if (valid && dmin_sq != 0 && dist_sq >= dmin_sq) return 0;
		dmin_sq = dist_sq;
	}
	xpos    = x2;
	ypos    = y2;
	valid   = 1;
	return 1;
}


point destination_marker::get_pos() const {

	assert(!point_outside_mesh(xpos, ypos));
	return point(get_xval(xpos), get_yval(ypos), (mesh_height[ypos][xpos] + object_types[SMILEY].radius));
}


// ********** SMILEY AI CODE (player_state) **********


bool proj_coll_test(point const &pos, point const &target_pos, vector3d const &orient, float radius, int weapon, int coll_id) {

	int xpos(0), ypos(0), index(0);
	point coll_pos;
	int const test_alpha((weapon == W_LASER) ? 1 : 3);
	point const pos2(target_pos - orient*(1.2*radius));
	if (!coll_pt_vis_test(pos, pos2, 1.2*radius, index, coll_id, 0, test_alpha)) return 0; // cobj collision
		
	if (get_range_to_mesh(pos, orient, coll_pos, xpos, ypos)) {
		if ((p2p_dist(pos, coll_pos) + 2.0*radius) < p2p_dist(pos, target_pos)) return 0; // mesh collision
	}
	return 1;
}


pos_dir_up get_smiley_pdu(point const &pos, vector3d const &orient) {
	
	float tterm, sterm;
	calc_view_test_terms(tterm, sterm, 0);
	return pos_dir_up(pos, orient, plus_z, tterm, sterm, NEAR_CLIP, FAR_CLIP);
}


bool check_left_and_right(point const &pos, point const &tpos, vector3d const &orient,
	float check_radius, float radius, int weapon, int coll_id)
{
	vector3d const check_dir(cross_product(orient, plus_z).get_norm()*check_radius);

	for (unsigned d = 0; d < 2; ++d) { // test left and right
		point const pos1(pos  + check_dir*(d ? 1.0 : -1.0));
		point const pos2(tpos + check_dir*(d ? 1.0 : -1.0));
		if (!proj_coll_test(pos1, pos2, orient, radius, weapon, coll_id)) return 0;
	}
	return 1;
}


void player_state::smiley_fire_weapon(int smiley_id) {

	if (!game_mode) return;
	int cid(coll_id[SMILEY]), status;
	dwobject const &smiley(obj_groups[cid].get_obj(smiley_id));
	if (smiley.disabled()) return;
	if (target_visible != 1 && (weapon != W_LANDMINE || (rand()&3) != 0)) return;
	assert(target >= CAMERA_ID && target < num_smileys);
	int const last_weapon(weapon);
	
	if (weapon == W_UNARMED || (!UNLIMITED_WEAPONS && no_weap_or_ammo())) {
		check_switch_weapon(smiley_id); // out of ammo, switch weapons
		if (weapon != last_weapon) fire_frame = 0;
	}
	if (target_visible && self_damage > 0.0 && powerup != PU_SHIELD && weapons[weapon].self_damage) { // can damage self
		if (dist_less_than(target_pos, smiley.pos, (weapons[weapon].blast_radius + object_types[SMILEY].radius))) { // will damage self
			check_switch_weapon(smiley_id); // switch weapons to avoid suicide
			if (weapon != last_weapon) fire_frame = 0;
		}
	}
	assert(!no_weap_or_ammo());
	if (weapon == W_BALL && !UNLIMITED_WEAPONS && (rand()&15) != 0) return; // wait to throw
	point pos(smiley.pos);
	bool const underwater(is_underwater(pos));
	if (temperature <= W_FREEZE_POINT && underwater) return; // under ice
	weapon_t const &w(weapons[weapon]);
	if (!w.use_underwater && (underwater || is_underwater(target_pos))) return; // self or target underwater and not hittable
	float const radius(object_types[SMILEY].radius), vweap(w.get_fire_vel());
	point tpos(target_pos);
	vector3d orient;
	bool target_lead(0);

	if (LEAD_SHOTS && smiley_acc >= 0.9 && vweap > TOLERANCE) { // lead target shots
		//orient = lead_target(pos, tpos, velocity, sstates[target].velocity, vweap);
		float const dist(p2p_dist(pos, tpos)), hit_time(dist/vweap); // approximate because it ignores enemy vel
		vector3d const enemy_vel(sstates[target].velocity);
		tpos += enemy_vel*hit_time; // the predicted enemy location
		target_lead = 1;
	}

	// maybe aim up to account for gravity
	if (smiley_acc <= 0.0) {
		orient = smiley.orientation;
	}
	else if (smiley_acc >= 0.5 && vweap > TOLERANCE && weapon != W_LANDMINE && w.obj_id != UNDEF) { // shrapnel?
		float const rel_enemy_vel(target_lead ? 0.0 : get_rel_enemy_vel(pos));
		if (rel_enemy_vel > vweap) return; // should already have been tested
		vector3d const tdir((tpos - pos).get_norm());
		float const wvel(vweap - rel_enemy_vel), radius2(radius + object_types[w.obj_id].radius), gscale(object_types[w.obj_id].gravity);

		if (gscale > 0.0 && wvel > 0.0) {
			point const fpos(pos + tdir*(0.75*radius2));

			if (ACCURATE_GRAV_PREDICT) {
				orient = get_firing_dir(fpos, tpos, wvel, gscale); // more accurate
				if (orient == all_zeros) return; // out of range
			}
			else {
				float const dist(p2p_dist(fpos, tpos));
				float const len(gscale*base_gravity*GRAVITY * dist*dist / (2*wvel*wvel)); // simpler and more efficient
				orient = (tpos + point(0.0, 0.0, len) - pos).get_norm();
			}
			// test line of sight here before using orient to help exclude invalid trajectories
			if (!proj_coll_test(pos, tpos, tdir, radius, weapon, smiley.coll_id)) return;
			float const proj_radius(object_types[w.obj_id].radius);
			if (!check_left_and_right(pos, tpos, tdir, proj_radius, radius, weapon, smiley.coll_id)) return;
		}
		else {
			orient = tpos - pos;
		}
	}
	else {
		bool const using_shrapnel((wmode&1) && (weapon == W_SHOTGUN || weapon == W_M16));
		float aim_up_val(0.0);
		if      (using_shrapnel   ) aim_up_val = 0.5;
		else if (w.obj_id != UNDEF) aim_up_val = object_types[w.obj_id].gravity;
		pos.z += CLIP_TO_01(0.1f +  aim_up_val)*radius; // shoot slightly higher than center
		orient = tpos - pos;
	}
	float const target_dist(p2p_dist(tpos, pos));

	if (smiley_acc > 0.0 && smiley_acc < 1.0) { // add firing error (should this be before or after the range test?)
		vadd_rand(orient, 0.1*((game_mode == 2) ? 0.5 : 1.0)*(1.0 - smiley_acc)/max(0.2f, min(1.0f, target_dist)));
	}
	if (target == CAMERA_ID && weapon == W_LASER) {
		// less accurate when shooting at the player, to make it more fair?
	}
	orient.normalize();
	
	if (weapon != W_LANDMINE && weapon != W_BBBAT && target_dist > 2.0*radius) {
		// make sure it has a clear shot (excluding invisible smileys)
		if (!proj_coll_test(pos, tpos, orient, radius, weapon, smiley.coll_id)) return; // Note: inexact, fails to account for gravity

		// check if we need to fire above or to the side to avoid a projectile collision with an obstacle
		if (weapon == W_ROCKET || weapon == W_SEEK_D || weapon == W_PLASMA) { // large projectile
			assert(w.obj_id != UNDEF);
			float const proj_radius(object_types[w.obj_id].radius);
			point const pos1(pos  + point(0.0, 0.0, -proj_radius)); // proj_radius up (+z)
			point const pos2(tpos + point(0.0, 0.0, -proj_radius));

			if (!proj_coll_test(pos1, pos2, orient, radius, weapon, smiley.coll_id)) {
				orient   *= target_dist;
				orient.z += min(proj_radius, 0.7f*radius); // shoot slightly upward
				orient.normalize();
			}
			if (!check_left_and_right(pos, tpos, orient, proj_radius, radius, weapon, smiley.coll_id)) return;
		}
	}
	int &ammo(p_ammo[weapon]);
	assert(ammo >= 0);
	if (weapon == W_GRENADE && (wmode&1) && ammo < 3) wmode = 0;
	int chosen;
	status = fire_projectile(pos, orient, smiley_id, chosen);

	if (status != 0 && !UNLIMITED_WEAPONS && !no_weap_or_ammo() && w.need_ammo) {
		ammo -= (weapon == W_GRENADE && (wmode&1)) ? 3 : 1; // check for cluster grenade
		assert(ammo >= 0);
		if (ammo == 0) fire_frame = 0; // could switch weapons
	}
}


void add_target(vector<od_data> &oddatav, pos_dir_up const &pdu, point const &pos2, float radius, int id, int hitter, int killer) {

	if (!sphere_in_view(pdu, pos2, radius, 0)) return; // view culling
	float dist_sq(p2p_dist_sq(pdu.pos, pos2));
	if (hitter == id) dist_sq *= 0.25; // perfer to attack your attacker
	if (killer == id) dist_sq *= 0.5;  // perfer to attack your last killer in revenge
	oddatav.push_back(od_data(((id == CAMERA_ID) ? CAMERA : SMILEY), id, dist_sq));
}


int player_state::find_nearest_enemy(point const &pos, pos_dir_up const &pdu, point const &avoid_dir,
	int smiley_id, point &target, int &target_visible, float &min_dist) const
{
	assert(smiley_id < num_smileys);
	int min_i(NO_SOURCE);
	int const last_hitter(was_hit ? hitter : NO_SOURCE);
	float const radius(object_types[SMILEY].radius);
	int const cid(coll_id[SMILEY]);
	point const camera(get_camera_pos());
	static vector<od_data> oddatav;
	min_dist = 0.0;

	if (free_for_all) { // smileys attack each other, not only the player
		assert((int)obj_groups[cid].max_objects() == num_smileys);

		for (int i = obj_groups[cid].max_objects()-1; i >= 0; --i) {
			dwobject const &obj(obj_groups[cid].get_obj(i));
			if (obj.disabled() || i == smiley_id || same_team(smiley_id, i))      continue;
			if (last_hitter != i && sstates[i].powerup == PU_INVISIBILITY)        continue; // invisible
			if (last_hitter != i && is_in_darkness(obj.pos, radius, obj.coll_id)) continue; // too dark to be visible
			add_target(oddatav, pdu, obj.pos, radius, i, last_hitter, killer);
		}
	}
	if (camera_mode != 0 && !spectate && !same_team(smiley_id, CAMERA_ID) && (sstates[CAMERA_ID].powerup != PU_INVISIBILITY || last_hitter == CAMERA_ID)) {
		if (!is_in_darkness(camera, radius, camera_coll_id)) {
			add_target(oddatav, pdu, camera, radius, CAMERA_ID, last_hitter, killer); // camera IN_DARKNESS?
		}
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


void player_state::check_cand_waypoint(point const &pos, point const &avoid_dir, int smiley_id,
	vector<od_data> &oddatav, unsigned i, int curw, float dmult, pos_dir_up const &pdu, bool next, float max_dist_sq)
{
	assert(i < waypoints.size());
	point const &wp(waypoints[i].pos);
	float const dist_sq(p2p_dist_sq(pos, wp));
	if (max_dist_sq > 0.0 && dist_sq > max_dist_sq && i != curw)               return; // too far away
	bool const can_see(next || (on_waypt_path && i == curw));
	if (!is_over_mesh(wp) || is_underwater(wp))                                return; // invalid smiley location
	if (WAYPT_VIS_LEVEL[can_see] == 0 && !sphere_in_view(pdu, wp, 0.0, 0))     return; // view culling - more detailed query later
	if (avoid_dir != zero_vector && dot_product_ptv(wp, pos, avoid_dir) > 0.0) return; // need to avoid this directio
	unsigned other_smiley_targets(0);

	for (int s = 0; s < num_smileys; ++s) {
		if (s != smiley_id && sstates[s].last_waypoint == i) ++other_smiley_targets;
	}
	dmult *= (1.0 + 1.0*other_smiley_targets); // increase distance cost if other smileys are going for the same waypoint
	map<unsigned, count_t>::const_iterator it(blocked_waypts.find(i));
	if (it != blocked_waypts.end())     dmult *= (1.0 + (1 << it->second.c)); // exponential increase in cost for blocked waypoints
	if (!waypts_used.is_valid(i))       dmult *= 100.0;
	if (waypoints[i].next_wpts.empty()) dmult *= 10.0; // increase the cost of waypoints disconnected from the rest of the waypoint graph
	if (i == curw)                      dmult *= 1.0E-6; // prefer the current waypoint to avoid indecision and force next connections
	float const time_weight(tfticks - waypoints[i].get_time_since_last_visited(smiley_id));
	float const tot_weight(dmult*(0.5*time_weight + dist_sq)*rand_uniform(0.8, 1.2));
	oddatav.push_back(od_data(WAYPOINT, i, tot_weight, can_see)); // add high weight to prefer other objects
}


// health, shields, powerup, weapon, ammo, pack, waypoint
int player_state::find_nearest_obj(point const &pos, pos_dir_up const &pdu, point const &avoid_dir, int smiley_id,
	point &target_pt, float &min_dist, vector<type_wt_t> types, int last_target_visible, int last_target_type)
{
	assert(smiley_id < num_smileys);
	int min_ic(-1);
	float sradius(object_types[SMILEY].radius), ra_smiley(C_STEP_HEIGHT*sradius);
	static vector<od_data> oddatav;
	min_dist = 0.0;

	// process dynamic pickup objects and waypoints
	for (unsigned t = 0; t < types.size(); ++t) {
		unsigned const type(types[t].type);
		assert(t < NUM_TOT_OBJS);
		float const dmult(1.0/types[t].weight);

		if (type == WAYPOINT) { // process waypoints
			int curw(last_waypoint);
			int ignore_w(-1);
			// mode: 0: none, 1: user wpt, 2: placed item wpt, 3: goal wpt, 4: wpt index, 5: closest wpt, 6: closest visible wpt, 7: goal pos (new wpt)
			wpt_goal goal((has_wpt_goal ? 3 : 2), 0, all_zeros); // mode, wpt, goal_pos
			//wpt_goal goal(6, 0, get_camera_pos()-point(0.0, 0.0, camera_zh)); // closest wpt visible to camera

			if (last_target_visible && last_target_type != 3 && goal.mode <= 2) { // have a previous enemy/item target and no real goal
				goal.mode = 6; // closest visible waypoint
				goal.pos  = target_pos; // should still be valid
			}
			if (goal.mode <= 2) { // add waypoint to a team member engaging an enemy
				for (int i = 0; i < num_smileys; ++i) { // what about camera/player (CAMERA_ID)?
					if (i == smiley_id || !same_team(i, smiley_id)) continue;
					player_state const &ss(sstates[i]);
					if (!ss.target_visible || ss.target_type != 1 || ss.target == NO_SOURCE) continue;
					goal.mode = 6;
					goal.pos  = ss.target_pos;
				}
			}
			if (curw >= 0) { // currently targeting a waypoint
				assert((unsigned)curw < waypoints.size());

				if (dist_less_than(waypoints[curw].pos, pos, sradius)) { // smiley has reached waypoint
					//cout << "reached target waypoint " << curw << " at time " << tfticks << endl; // testing
					waypts_used.insert(curw); // insert as the last used waypoint and remove from consideration
					waypoints[curw].mark_visited_by_smiley(smiley_id);
					unreachable[1].clear();
					waypt_adj_vect const &next(waypoints[curw].next_wpts);

					if (!next.empty()) { // choose next waypoint from graph
						//cout << "choose next waypoint, curw: " << curw << endl;
						// FIXME: skip path waypoints that are in unreachable[1]?
						curw = find_optimal_next_waypoint(curw, goal); // can return -1
						//cout << "next curw: " << curw << endl;

						for (unsigned i = 0; i < next.size(); ++i) {
							check_cand_waypoint(pos, avoid_dir, smiley_id, oddatav, next[i], curw, dmult, pdu, 1, 0.0);
						}
						//cout << "size: " << oddatav.size() << endl;
						continue;
					}
					// disconntected waypoint - should rarely get here
					// don't pick this goal again until another waypoint is reached, otherwise we will get stuck here until another objective appears
					unreachable[1].add(waypoints[curw].pos);
					ignore_w = curw;
					curw     = -1;
				}
			}
			float const max_dist(0.25*(X_SCENE_SIZE + Y_SCENE_SIZE)), max_dist_sq(max_dist*max_dist);

			for (unsigned i = 0; i < waypoints.size(); ++i) { // inefficient - use subdivision?
				if (waypoints[i].disabled || i == ignore_w) continue;
				check_cand_waypoint(pos, avoid_dir, smiley_id, oddatav, i, curw, dmult, pdu, 0, max_dist_sq);
			}
			if (curw < 0) find_optimal_waypoint(pos, oddatav, goal);
		}
		else { // not a waypoint (pickup item)
			if (UNLIMITED_WEAPONS && (type == WEAPON || type == AMMO || type == WA_PACK)) continue;
			obj_group const &objg(obj_groups[coll_id[type]]);
			if (!objg.is_enabled()) continue;
			float const radius(object_types[type].radius);

			for (unsigned i = 0; i < objg.end_id; ++i) {
				dwobject const &obj(objg.get_obj(i));
				bool const placed((obj.flags & USER_PLACED) != 0);
				if (obj.disabled())                                              continue;
				if (!is_over_mesh(obj.pos) || (island && (obj.pos.z < ocean.z))) continue;
				if (!placed && !sphere_in_view(pdu, obj.pos, radius, 0))         continue; // view culling (disabled for predef object locations)
				if (avoid_dir != zero_vector && dot_product_ptv(obj.pos, pos, avoid_dir) > 0.0) continue; // need to avoid this direction
				if (is_in_darkness(obj.pos, radius, obj.coll_id))                continue;
				float cost((target_pos == obj.pos) ? 0.75 : 1.0); // favor original targets

				for (int j = 0; j < num_smileys; ++j) { // too slow?
					if (j != smiley_id && sstates[j].objective_pos == obj.pos) {
						float const dist_ratio(p2p_dist(pos, obj.pos)/p2p_dist(obj_groups[coll_id[SMILEY]].get_obj(j).pos, obj.pos));
						if (dist_ratio > 1.0) cost += dist_ratio; // icrease the cost since the other smiley will likely get there first
					}
				}
				oddatav.push_back(od_data(type, i, cost*dmult*p2p_dist_sq(pos, obj.pos)));
			}
		}
	} // for t
	if (oddatav.empty()) return min_ic;
	sort(oddatav.begin(), oddatav.end());

	for (unsigned i = 0; i < oddatav.size(); ++i) {
		int const type(oddatav[i].type), id(oddatav[i].id);
		assert(type >= 0 && type < NUM_TOT_OBJS);
		point pos2;
		dwobject const *obj(NULL);
		bool const is_wpt(type == WAYPOINT);
		bool no_frustum_test(0), skip_vis_test(0);

		if (is_wpt) {
			assert(size_t(id) < waypoints.size());
			bool const can_see(oddatav[i].val != 0 || (on_waypt_path && id == last_waypoint));
			pos2 = waypoints[id].pos;
			no_frustum_test = (WAYPT_VIS_LEVEL[can_see] == 1);
			skip_vis_test   = (WAYPT_VIS_LEVEL[can_see] == 2);
		}
		else {
			int const cid(coll_id[type]);
			assert(cid >= 0 && cid < NUM_TOT_OBJS);
			obj  = &obj_groups[cid].get_obj(id);
			pos2 = obj->pos;
			no_frustum_test = ((obj->flags & USER_PLACED) != 0);
		}
		int const xpos(get_xpos(pos2.x)), ypos(get_ypos(pos2.y));
		if (point_outside_mesh(xpos, ypos))       continue;
		if (unreachable[is_wpt].cant_reach(pos2)) continue;
		float const oradius(object_types[type].radius);
		bool const ice(temperature <= W_FREEZE_POINT);
		bool const not_too_high( pos2.z < max((mesh_height[ypos][xpos] + sradius + oradius), (pos.z + ra_smiley + oradius - sradius)));
		bool const not_too_high2(pos2.z < interpolate_mesh_zval(pos2.x, pos2.y, oradius, 0, 0)      + ra_smiley + oradius);
		bool const ice_height_ok(ice && wminside[ypos][xpos] && (pos2.z < water_matrix[ypos][xpos]  + ra_smiley + oradius));
		bool const in_motion(obj && obj->status == 1 && obj->velocity.mag_sq() > 9.0);
		bool const can_reach(not_too_high || ice_height_ok || in_motion);

		if (is_wpt || not_too_high2) { // almost in reach, ultimate reachability questionable
			// Note: if a waypoint is found to be unreachable, it may be be avoided until the smiley dies
			if (!unreachable[is_wpt].proc_target(pos, pos2, objective_pos, (can_reach && !is_wpt))) continue;
		}
		if (can_reach || not_too_high2 || (is_wpt && on_waypt_path)) { // not_too_high2 - may be incorrect
			int const max_vis_level(is_wpt ? 3 : ((type == BALL) ? 5 : 4));
			bool skip_path_comp(target_pos == pos2 && (rand()&15) != 0); // infrequent updates if same target

			if ((skip_vis_test || sphere_in_view(pdu, pos2, oradius, max_vis_level, no_frustum_test)) &&
				(powerup == PU_FLIGHT || skip_path_comp || is_valid_path(pos, pos2, !is_wpt)))
			{
				// select this object as our target and return
				if (is_wpt) {
					if (id != last_waypoint) on_waypt_path = (oddatav[i].val != 0);
					last_waypoint = id;
					last_wpt_dist = p2p_dist_xy(pos, pos2);
				}
				min_dist = sqrt(oddatav[i].dist); // find closest reachable/visible object
				min_ic    = type;
				target_pt = pos2;
				break;
			}
			else if (is_wpt && id == last_waypoint) {
				++blocked_waypts[id].c;
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


int player_state::check_smiley_status(dwobject &obj, int smiley_id) {

	assert(sstates != NULL);
	if (obj.disabled()) return 0;

	if (obj.time > object_types[SMILEY].lifetime) {
		drop_pack(obj.pos);
		obj.status = 0; // smiley dies of old age (or boredom) after a REALLY long time
		string const msg(name + " died of " + ((rand()&1) ? "boredom" : "old age"));
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


void player_state::drop_pack(point const &pos) {

	if (UNLIMITED_WEAPONS) return; // no pack
	int const ammo(p_ammo[weapon]);
	if (!weapons[weapon].need_weapon && (!weapons[weapon].need_ammo || ammo == 0)) return; // no weapon/ammo
	bool const dodgeball(game_mode == 2 && weapon == W_BALL); // drop their balls
	if (dodgeball) assert(ammo == balls.size());
	unsigned const num(dodgeball ? ammo : 1);
	int const type(dodgeball ? BALL : WA_PACK), cid(coll_id[type]);
	obj_group &objg(obj_groups[cid]);

	for (unsigned i = 0; i < num; ++i) {
		int const max_t_i(dodgeball ? balls[i] : objg.choose_object());
		dwobject &obj(objg.get_obj(max_t_i));
		if (dodgeball) assert(obj.status == OBJ_STAT_RES);
		objg.create_object_at(max_t_i, pos);
		obj.direction = (unsigned char)(weapon + 0.5);
		obj.angle     = (dodgeball ? 1.0 : (ammo + 0.5));
		obj.velocity  = gen_rand_vector(1.0, 6.0, PI_TWO);
	}
	if (dodgeball) balls.clear();
}


int player_state::drop_weapon(vector3d const &coll_dir, vector3d const &nfront, point const &pos, int index, float energy, int type) {

	if (game_mode != 1) return 0;
	if (type == BEAM || type == BLAST_RADIUS || (type >= DROWNED && type <= CRUSHED)) return 0;

	if ((rand()%20 == 0) && energy > 25.0 && (weapons[weapon].need_weapon || (weapons[weapon].need_ammo && p_ammo[weapon] > 0))) {
		float const frontv(dot_product(coll_dir, nfront)/(coll_dir.mag()*nfront.mag()));

		if (frontv > 0.95) {
			vector3d rv(signed_rand_float(), signed_rand_float(), 1.2*rand_float());
			point const dpos(pos + rv*(2.0*object_types[SMILEY].radius));
			drop_pack(dpos);
			p_weapons[weapon] = 0;
			p_ammo[weapon]    = 0;
			if (index == CAMERA_ID) switch_weapon(1, 0); else check_switch_weapon(index);
			return 1;
		}
	}
	return 0;
}


vector3d get_avoid_dir(point const &pos, int smiley_id, pos_dir_up const &pdu) {

	// look for landmines and [c]grenades and avoid them
	point avoid_dir(zero_vector);
	int const check_ids[3] = {  GRENADE,   CGRENADE,   LANDMINE};
	int const weap_ids [3] = {W_GRENADE, W_CGRENADE, W_LANDMINE};

	for (unsigned t = 0; t < 3; ++t) {
		int const type(check_ids[t]);
		obj_group const &objg(obj_groups[coll_id[type]]);
	
		if (objg.is_enabled()) {
			float min_dist(weapons[weap_ids[t]].blast_radius);
			if (type == W_LANDMINE) min_dist *= 0.5; // trigger radius is only about half the blast radius

			for (unsigned i = 0; i < objg.end_id; ++i) {
				dwobject const &obj2(objg.get_obj(i));
			
				if (!obj2.disabled() && obj2.source == smiley_id && (type != LANDMINE || !obj2.lm_coll_invalid()) &&
					dist_less_than(pos, obj2.pos, min_dist) && sphere_in_view(pdu, pos, object_types[type].radius, 0))
				{
					avoid_dir  = (obj2.pos - pos);
					min_dist   = avoid_dir.mag();
					avoid_dir /= min_dist;
					break; // can only handle one right now
				}
			}
		}
	}
	return avoid_dir;
}


void player_state::smiley_select_target(dwobject &obj, int smiley_id) {

	if (obj.disabled()) return;
	int min_ie(NO_SOURCE), min_ih(-1);
	float diste, disth, health(obj.health);
	vector<type_wt_t> types;
	point targete(all_zeros), targeth(all_zeros);
	int const last_target_visible(target_visible), last_target_type(target_type);
	target_visible = 0;
	target_type    = 0; // 0 = none, 1 = enemy, 2 = health/powerup, 3 = waypoint
	float const health_eq(min(4.0f*health, (health + shields)));
	bool const almost_dead(health_eq < 20.0);
	pos_dir_up const pdu(get_smiley_pdu(obj.pos, obj.orientation));
	vector3d const avoid_dir(get_avoid_dir(obj.pos, smiley_id, pdu));

	if (game_mode == 2 && !UNLIMITED_WEAPONS) { // want the ball
		if (p_ammo[W_BALL] > 0) { // already have a ball
			min_ie = find_nearest_enemy(obj.pos, pdu, avoid_dir, smiley_id, target_pos, target_visible, diste);
		}
		if (!target_visible) { // don't have a ball or no enemy in sight
			types.push_back(type_wt_t(BALL, 1.0));
			min_ih = find_nearest_obj(obj.pos, pdu, avoid_dir, smiley_id, target_pos, disth, types, last_target_visible, last_target_type);
		}
	}
	else { // choose health/attack based on distance
		// want to get powerup badly if you don't have one
		float const pu_wt((powerup == 0) ? 1.5 : (1.0 - float(powerup_time)/POWERUP_TIME));
		types.push_back(type_wt_t(POWERUP, pu_wt));
		types.push_back(type_wt_t(WEAPON,  0.8));
		types.push_back(type_wt_t(AMMO,    0.7));
		types.push_back(type_wt_t(WA_PACK, 1.0));
		types.push_back(type_wt_t(SHIELD,  (almost_dead ? 10 : 1.2)*(1.0 - shields/MAX_SHIELDS))); // always below max since it ticks down over time
		if (health < MAX_HEALTH) types.push_back(type_wt_t(HEALTH, (almost_dead ? 15 : 1.5)*(1.0 - health/MAX_HEALTH)));

		if (game_mode) {
			min_ie = find_nearest_enemy(obj.pos, pdu, avoid_dir, smiley_id, targete, target_visible, diste);
		}
		min_ih = find_nearest_obj(  obj.pos, pdu, avoid_dir, smiley_id, targeth, disth, types, last_target_visible, last_target_type);

		if (!target_visible) { // can't find an enemy, choose health/pickup
			if (min_ih >= 0) target_pos = targeth;
		}
		else if (min_ih < 0) { // can't find health/pickup, choose attack
			if (target_visible) target_pos = targete;
		}
		else { // found both item and enemy
			float const dp(dot_product((obj.pos - targeth).get_norm(), (obj.pos - targete).get_norm()));

			if (diste <= disth || dp > 0.95) { // enemy is closer or in the direction of the desired pickup
				target_pos = targete;
				min_ih     = -1;
			}
			else { // health/pickup is closer
				target_pos     = targeth;
				target_visible = 0;
			}
		}
	}
	if (target_visible) { // targeting and enemy
		target_type   = 1;
		objective_pos = target_pos;
	}
	else if (was_hit && hit_dir != zero_vector) { // turn around if hit in the back
		target_pos    = obj.pos + hit_dir.get_norm();
		target_type   = 1; // enemy
		objective_pos = target_pos;
	}
	else { // targeting an item or nothing
		if (min_ih < 0) { // no item - use a waypoint to move to a different place in the scene
			types.clear();
			types.push_back(type_wt_t(WAYPOINT, 1.0));
			min_ih = find_nearest_obj(obj.pos, pdu, avoid_dir, smiley_id, target_pos, disth, types, last_target_visible, last_target_type);
		}
		objective_pos = target_pos;
		
		if (min_ih >= 0) {
			target_type    = (min_ih == WAYPOINT) ? 3 : 2;
			target_visible = 2;
		}
		/*if (target_pos != obj.pos) { // look beyond the target
			vector3d const o(target_pos - obj.pos); // motion direction
			target_pos += o*(2.0*object_types[SMILEY].radius/o.mag());
		}*/
	}
	if (min_ih != WAYPOINT) {
		unreachable[1].reset_try(); // no target waypoint
		reset_wpt_state();
	}
	else if (min_ih < 0) {
		unreachable[0].reset_try(); // no target item
	}
	if (on_waypt_path)       blocked_waypts.clear();
	if (target_visible == 1) assert(min_ie >= CAMERA_ID);
	target = min_ie;
}


bool is_targeting_smiley(int targeter, int targetee, point const &targetee_pos) {

	assert(targeter >= CAMERA_ID && targeter < num_smileys);
	assert(targetee >= CAMERA_ID && targetee < num_smileys);
	assert(targeter != targetee);
	if (targeter == CAMERA_ID) return sphere_in_camera_view(targetee_pos, 0.0, 0); // camera is targeting if target is in view
	return (sstates[targeter].target_type == 1 && sstates[targeter].target == targetee);
}


float player_state::get_pos_cost(int smiley_id, point pos, point const &opos, pos_dir_up const &pdu,
	float radius, float step_height, bool check_dists)
{
	// target_type: 0=none, 1=enemy, 2=item, 3=waypoint
	if (!is_over_mesh(pos)) return 10.0; // off the mesh - high cost
	int xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos)) return 10.0; // off the mesh - high cost

	if (powerup != PU_FLIGHT) {
		if (!check_step_dz(pos, opos, radius)) return 9.0;

		if (!on_waypt_path && temperature > W_FREEZE_POINT && !is_mesh_disabled(xpos, ypos)) {
			float depth(0.0);
			if (is_underwater(pos, 0, &depth)) return 6.0 + 0.01*depth; // don't go under water/blood
		}
	}
	vector3d const avoid_dir(get_avoid_dir(pos, smiley_id, pdu));
	if (avoid_dir != zero_vector) return 5.0 + 0.1*dot_product(avoid_dir, (pos - opos).get_norm());
	if (powerup != PU_FLIGHT && !can_make_progress(pos, opos, !on_waypt_path)) return 4.0; // can't make progress in this direction

	if (target_type == 1 && target != NO_SOURCE) { // don't get too close to enemy with ranged weapons
		assert(target >= CAMERA_ID && target < num_smileys && target != smiley_id);
		int cindex(-1);
		if (check_coll_line(pos, target_pos, cindex, -1, 1, 0)) return 3.5; // target not visible / no line of sight
		float const dist(p2p_dist(pos, target_pos));
		float const range(weapon_range(1));
		if (check_dists && dist > range) return (3.0 + 0.01*dist); // out of weapon range
		// FIXME: use target_in_range()?

		if (is_targeting_smiley(target, smiley_id, pos)) { // enemy (or camera) is targeting me
			float const enemy_range(sstates[target].weapon_range(1)); // use sstates[target].target_in_range?
			if (range > 2.0*enemy_range && dist < enemy_range) return (2.5 + 0.01*(enemy_range - dist)); // too close to enemy weapon range
		}
		float const min_dist(min(0.25*range, 8.0*radius));
		if (weapon != W_BBBAT && dist < min_dist) return (2.0 + 0.01*(min_dist    - dist)); // too close to enemy
	}
	// enforce turn speed?
	return 0.0; // good
}


void add_rand_dir_change(vector3d &dir, float mag) {

	dir += vector3d(mag*signed_rand_float(), mag*signed_rand_float(), 0.0);
	dir.normalize();
}


struct dir_cost_t {
	float cost, dp;
	vector3d dir;

	dir_cost_t() : cost(0.0), dp(0.0) {}
	dir_cost_t(float cost_, vector3d const &dir_, vector3d const &opt_dir)
		: cost(cost_), dp(dot_product(dir_, opt_dir)), dir(dir_) {}

	bool operator<(dir_cost_t const &dct) const {
		return ((cost == dct.cost) ? (dp > dct.dp) : (cost < dct.cost));
	}
};


vector3d step_dist_scale(dwobject const &obj, vector3d const &dir) {

	float const dp(dot_product(obj.orientation, dir));
	if (dp >  0.5) return dir;
	if (dp < -0.5) return dir*BACKWARD_SPEED;
	return                dir*SIDESTEP_SPEED;
}


int player_state::smiley_motion(dwobject &obj, int smiley_id) {

	if (NO_SMILEY_ACTION || obj.disabled()) return 0;
	float const speed(get_rspeed_scale()), radius(object_types[SMILEY].radius);
	assert(radius >= 0.0);
	float const step_dist(speed*smiley_speed*fticks*GROUND_SPEED), step_height(C_STEP_HEIGHT*radius);
	point const opos(obj.pos);
	int xpos(get_xpos_clamp(opos.x)), ypos(get_ypos_clamp(opos.y));
	bool const has_flight(powerup == PU_FLIGHT), is_water_temp(temperature > W_FREEZE_POINT), underwater(is_underwater(opos));
	assert(!point_outside_mesh(xpos, ypos));
	bool stuck(0), in_ice(0), no_up(0), no_down(0), using_dest_mark(0);
	
	// check for stuck underwater
	if (is_water_temp && underwater && !on_waypt_path) { // ok if on a waypoint path
		dest_mark.update_dmin(xpos, ypos);
		int xt(dest_mark.xpos), yt(dest_mark.ypos); // start with the last dest_mark used, if any (to remember a valid dry spot)

		// find a valid dest mark
		for (unsigned i = 0; i < SMILEY_MAX_TRIES; ++i) { // find some randomly chosen dry spots to head for
			if (dest_mark.valid && dest_mark.min_depth == 0.0) break;

			if (i > 0 || !point_interior_to_mesh(xt, yt)) { // need to calculate a new one
				xt = 1 + (rand() % (MESH_X_SIZE-2));
				yt = 1 + (rand() % (MESH_Y_SIZE-2));
			}
			float const depth(has_water(xt, yt) ? (water_matrix[yt][xt] - mesh_height[yt][xt]) : 0.0);
			dest_mark.add_candidate(xpos, ypos, xt, yt, depth, radius);
		}
		if (dest_mark.valid) {
			target_visible = 1; // ???
			target_type    = 2; // like an item
			objective_pos  = target_pos = dest_mark.get_pos();
			using_dest_mark= 1;
		}
		if (dest_mark.valid) assert(dest_mark.xpos > 0);
	}
	else {
		dest_mark.clear();
	}

	// do actual movement
	if (!is_water_temp && underwater) {
		obj.velocity = zero_vector; // stuck in ice, can't move
		stuck = in_ice = 1;
	}
	else {
		if (target_visible && dist_less_than(obj.pos, target_pos, step_dist)) {
			obj.pos = target_pos; // close enough to hit the target exactly
		}
		else {
			vector3d stepv;
			
			if (target_visible) {
				stepv = (target_pos - obj.pos).get_norm();
			}
			else if (!(xpos > 1 && ypos > 1 && xpos < MESH_X_SIZE-1 && ypos < MESH_Y_SIZE-1)) {
				stepv = vector3d(-obj.pos.x, -obj.pos.y, 0.0).get_norm(); // move toward the center of the scene (0,0)
			}
			else {
				stepv = obj.orientation; // continue along the same orient if no target visible
			}
			if (!on_waypt_path && ((rand() % 200) == 0)) {
				add_rand_dir_change(stepv, 1.0); // "dodging"
				stepv.normalize();
			}
			obj.pos += step_dist_scale(obj, stepv)*step_dist;

			// check if movement was valid
			pos_dir_up const pdu(get_smiley_pdu(obj.pos, obj.orientation));
			float const start_cost(using_dest_mark ? 0.0 : get_pos_cost(smiley_id, obj.pos, opos, pdu, radius, step_height, 0));
		
			if (start_cost > 0.0) {
				//cout << "cost: " << start_cost << ", stepv: "; stepv.print(); cout << endl; // testing
				unsigned const ndirs(16);
				dir_cost_t best;

				for (unsigned i = 0; i < ndirs; ++i) {
					vector3d dir(stepv);
					rotate_vector3d_norm(plus_z, TWO_PI*i/ndirs, dir);
					float const cost(get_pos_cost(smiley_id, (opos + step_dist_scale(obj, dir)*step_dist), opos, pdu, radius, step_height, 1)); // FIXME: zval has not been set
					dir_cost_t const cur(cost, dir, stepv);
					if (i == 0 || cur < best) best = cur;
				}
				//if (best.dp < 1.0) {cout << "best: cost: " << best.cost << ", dp: " << best.dp << ", dir: "; best.dir.print(); cout << endl;}
				if (best.cost > 0.0) {} // FIXME: still not good, what to do?
				obj.pos = opos + step_dist_scale(obj, best.dir)*step_dist;
			}
		}
		obj.velocity.assign(0.0, 0.0, -1.0);
	}
	point const wanted_pos(obj.pos);

	// perform collision detection and mesh constraints to get z height
	int const coll(obj.multistep_coll(opos, smiley_id, SMILEY_COLL_STEPS));

	if (has_flight) {
		no_down = (obj.pos.z > wanted_pos.z);
		no_up   = (obj.pos.z < wanted_pos.z);
	}
	else {
		player_clip_to_scene(obj.pos);
	}
	if (obj.status == 0) { // killed
		remove_reset_coll_obj(obj.coll_id);
		return 0;
	}
	if (coll && p2p_dist_sq(obj.pos, opos) < p2p_dist_sq(obj.pos, wanted_pos)) {
		if (++stopped_time > 4) stuck = 1;
	}
	else {
		stopped_time = 0;
	}
	int ohval(set_true_obj_height(obj.pos, opos, C_STEP_HEIGHT, zvel, SMILEY, smiley_id, has_flight, 0));

	if (!has_flight && !recreated && !mesh_scale_change && (obj.pos.z - opos.z) >= step_height) { // stuck, can't step that high
		obj.pos = opos; // reset to old position
		ohval   = 3;
	}
	if (has_flight && !in_ice) {
		if ((target_type >= 2 && target_pos.z < opos.z && obj.pos.z < opos.z) || // want health/powerup/waypoint below - go down
			(target_type == 1 && weapon == W_BBBAT)) // have to go down for baseball bat hit
		{
			obj.pos.z = max(max(target_pos.z, obj.pos.z), (opos.z - (no_down ? 0.0f : float(radius*rand_uniform(0.2, 0.6)))));
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

	// set orientation
	if (in_ice) {
		obj.orientation = all_zeros; // nothing
	}
	else if (SMILEYS_LOOK_AT_TARGET && target_type != 0) {
		vector3d new_orient((target_pos - obj.pos).get_norm()); // target dir
		if (new_orient != zero_vector) obj.orientation = new_orient;
	}
	else {
		vector3d const new_orient((obj.pos - opos).get_norm()); // actual movement dir
		obj.orientation.x = new_orient.x;
		obj.orientation.y = new_orient.y;
		obj.orientation.z = 0.8*obj.orientation.z + 0.2*new_orient.z; // smooth z step
		float const mh(int_mesh_zval_pt_off(obj.pos, 1, 1));
		
		if (obj.pos.z < (mh + 1.1*radius)) { // walking on the mesh (approx)
			int const xpos(get_xpos(obj.pos.x)), ypos(get_ypos(obj.pos.y));

			if (!point_outside_mesh(xpos, ypos)) { // determine orient.z from mesh normal
				obj.orientation.z = -dot_product(vector3d(obj.orientation.x, obj.orientation.y, 0.0), surface_normals[ypos][xpos]);
			}
		}
		obj.orientation.normalize();
	}
	if (obj.orientation == zero_vector) {
		obj.orientation = vector3d(1,0,0);
	}
	else if ((stuck || ohval == 3) && !in_ice && target_type != 3 && !dist_less_than(obj.pos, target_pos, 2.1*radius)) { // not on a waypoint
		obj.orientation.negate(); // turn around - will this fix it or just oscillate orient when really stuck?
		add_rand_dir_change(obj.orientation, 0.25);
	}

	// set other state
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
	velocity   = (obj.pos - opos)/(TIMESTEP*fticks);
	obj.status = 3;
	return 1;
}


void advance_smiley(dwobject &obj, int smiley_id) {

	assert(smiley_id < num_smileys);
	sstates[smiley_id].advance(obj, smiley_id);
}


void player_state::advance(dwobject &obj, int smiley_id) { // seems to slightly favor smileys with later ids

	assert(obj.type == SMILEY);
	assert(obj_groups[coll_id[SMILEY]].enabled);
	if (!check_smiley_status(obj, smiley_id)) {fall_counter = 0; return;}
	smiley_select_target(obj, smiley_id);
	obj.time += iticks;
	if (!smiley_motion(obj, smiley_id)) {fall_counter = 0; return;}
	smiley_action(smiley_id);
}


void player_state::shift(vector3d const &vd) {

	target_pos    += vd;
	objective_pos += vd;
	for (unsigned i = 0; i < 2; ++i) unreachable[i].shift_by(vd);
}


void shift_player_state(vector3d const &vd, int smiley_id) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);
	sstates[smiley_id].shift(vd);
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
	if (sstates[smiley_id].tdata == NULL) sstates[smiley_id].tdata = new unsigned char[tbytes];
	unsigned char *const tdata(sstates[smiley_id].tdata);
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
	sstates[smiley_id].check_switch_weapon(smiley_id);
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


void player_state::check_switch_weapon(int smiley_id) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);
	int bbat_iter(0), bb_range(-1);
	wmode = ((rand()&3) == 0);

	if (game_mode == 2) { // dodgeball mode
		weapon     = ((UNLIMITED_WEAPONS || p_ammo[W_BALL] > 0) ? W_BALL : W_UNARMED);
		fire_frame = 0;
		return;
	}
	// we get here even when game_mode==0, but that's ok
	vector<pair<float, unsigned> > choices;
	point const &pos(obj_groups[coll_id[SMILEY]].get_obj(smiley_id).pos);
	float min_weight(0.0);
	unsigned chosen_weap(W_UNARMED);

	for (unsigned i = 1; i < NUM_WEAPONS-1; ++i) {
		weapon = i;
		if (no_weap_or_ammo()) continue;
		float weight(rand_float());
		if (!weapons[weapon].use_underwater && is_underwater(pos)) weight += 0.5;
		float const range(weapon_range(0));
		if (range > 0.0)  weight += ((target_in_range(pos) != 0) ? -0.2 : 0.8); // ranged weapon
		if (i == W_BBBAT) weight *= 1.5;
		if (i == W_SBALL) weight *= 1.2;
		
		if (chosen_weap == W_UNARMED || weight < min_weight) {
			min_weight  = weight;
			chosen_weap = weapon;
		}
	}
	fire_frame = 0;
	weapon     = chosen_weap;
	//weapon     = W_LASER; // set to force this as weapon choice (if available)
	if (weapon == W_PLASMA && wmode == 1 && (rand()%4) != 0) plasma_loaded = 1; // fire it up!
}


float player_state::get_rel_enemy_vel(point const &pos) const {

	if (target_type != 1 || target == NO_SOURCE) return 0.0;
	vector3d const enemy_vel(sstates[target].velocity);
	vector3d const enemy_dir((target_pos - pos).get_norm());
	return dot_product(enemy_vel, enemy_dir);
}


int player_state::target_in_range(point const &pos) const {

	if (!target_visible || target_type != 1 || target == NO_SOURCE) return 2;
	float range(weapon_range(0));
	
	if (weapons[weapon].obj_id == LANDMINE || weapons[weapon].obj_id == UNDEF) {
		return (range == 0.0 || dist_less_than(target_pos, pos, range));
	}
	assert(target >= CAMERA_ID && target < num_smileys);
	float const wvel(weapons[weapon].get_fire_vel()); // what about wmode==1?
	float const rel_enemy_vel(get_rel_enemy_vel(pos));
	assert(wvel > 0.0);
	if (rel_enemy_vel > wvel) return 0; // enemy is moving away faster than our projectile
	if (range == 0.0)         return 1;
	range *= (wvel - rel_enemy_vel)/wvel; // adjust range based on enemy velocity
	float const gravity(object_types[weapons[weapon].obj_id].gravity);
	if (gravity == 0.0) return dist_less_than(target_pos, pos, range); // no gravity, use simple distance test
	float const xy_dist_sq(p2p_dist_xy_sq(target_pos, pos)); // xy distance
	if (target_pos.z <= pos.z) return (xy_dist_sq < range*range); // shooting down, ignore the z (height) distance
	float const eff_dz((1.0 + gravity)*(target_pos.z - pos.z)); // increase the cost of the z distance due to gravity (approximate)
	return ((xy_dist_sq + eff_dz*eff_dz) < range*range);
}


void player_state::smiley_action(int smiley_id) {

	assert(smiley_id >= 0 && smiley_id < num_smileys);
	float depth(0.0);
	dwobject &smiley(obj_groups[coll_id[SMILEY]].get_obj(smiley_id));
	int const in_range(target_in_range(smiley.pos));
	if (in_range == 1) smiley_fire_weapon(smiley_id);
	if (powerup == PU_REGEN) smiley.health = min(MAX_REGEN_HEALTH, smiley.health + 0.1f*fticks);
	if ((rand()%((in_range == 0) ? 50 : 500)) == 0) check_switch_weapon(smiley_id); // change weapons
	if (was_hit > 0) --was_hit;
	kill_time += max(1, iticks);
	check_underwater(smiley_id, depth);
}


colorRGBA get_smiley_team_color(int smiley_id) {

	colorRGBA const tcolors[] = {RED, BLUE, GREEN, MAGENTA, CYAN, BROWN, PINK, GRAY, LT_BLUE, ORANGE, WHITE};
	if (teams <= 1) return tcolors[0];
	return tcolors[min(((smiley_id+1)%teams), int(sizeof(tcolors)/sizeof(colorRGBA)-1))];
}


int gen_smiley_or_player_pos(point &pos, int index) {

	float const radius(object_types[SMILEY].radius);

	if (!app_spots.empty()) {
		float const dmin(2.0*radius);
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
			pos.z = interpolate_mesh_zval(pos.x, pos.y, radius, 0, 0) + radius;
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


void init_sstate(int id, bool w_start) {

	assert(sstates != NULL && id >= CAMERA_ID && id < num_smileys);
	sstates[id].init(w_start);

	for (int i = CAMERA_ID; i < num_smileys; ++i) {
		if (sstates[i].target_visible == 1 && sstates[i].target == id) {
			sstates[i].target_visible = 0;
		}
	}
}


void init_smileys() {

	for (int i = 0; i < num_smileys; ++i) {
		init_smiley(i);
	}
}


bool has_invisibility(int id) {

	assert(id >= CAMERA_ID && id < num_smileys);
	if (!game_mode)      return 0;
	if (sstates == NULL) return 0; // not initialized - should this be an error?
	return (sstates[id].powerup == PU_INVISIBILITY);
}


void player_state::init(bool w_start) {

	assert(balls.empty());
	
	for (int i = 0; i < NUM_WEAPONS; ++i) {
		p_weapons[i] = 0;
		p_ammo[i]    = 0;
	}
	if (!UNLIMITED_WEAPONS) {
		if (w_start) {
			p_weapons[W_UNARMED] = 2;
			p_weapons[W_BBBAT]   = 2;
			p_weapons[W_SBALL]   = 1;
			p_ammo[W_SBALL]      = weapons[W_SBALL].def_ammo;
			weapon               = W_SBALL;
		}
		else {
			weapon = W_UNARMED;
		}
		wmode     = 0;
	}
	timer         = 0;
	init_frame    = frame_counter;
	fire_frame    = 0;
	was_hit       = 0;
	rot_counter   = 0;
	plasma_loaded = 0;
	uw_time       = 0;
	cb_hurt       = 0;
	target_visible= 0;
	target_type   = 0;
	target        = 0;
	plasma_size   = 1.0;
	zvel          = 0.0;
	stopped_time  = 0;
	fall_counter  = 0;
	last_dz       = 0.0;
	last_zvel     = 0.0;
	velocity      = zero_vector;
	kill_time     = 100*TICKS_PER_SECOND;

	if (game_mode == 1) {
		shields       = INIT_SHIELDS;
		powerup       = ((INIT_PU_SH_TIME > 0) ? PU_SHIELD : -1);
		powerup_time  = INIT_PU_SH_TIME;
	}
	else {
		shields       = 0.0;
		powerup       = -1;
		powerup_time  = 0;
	}
	//powerup = PU_FLIGHT; powerup_time = 1000000; // testing
	for (unsigned i = 0; i < 2; ++i) unreachable[i].clear();
	reset_wpt_state();
	waypts_used.clear();
	dest_mark.clear();
}


void player_state::reset_wpt_state() {

	last_waypoint = -1;
	on_waypt_path = 0;
	last_wpt_dist = 0.0;
	blocked_waypts.clear();
}


bool player_state::no_weap() const {

	if (!game_mode) return 0;
	assert(weapon < NUM_WEAPONS);
	assert(p_weapons[weapon] >= 0);
	return (!UNLIMITED_WEAPONS && weapons[weapon].need_weapon && p_weapons[weapon] == 0);
}


bool player_state::no_ammo() const {

	if (!game_mode) return 0;
	assert(weapon < NUM_WEAPONS);
	assert(p_ammo[weapon] >= 0);
	return (!UNLIMITED_WEAPONS && weapons[weapon].need_ammo && p_ammo[weapon] == 0);
}


float player_state::weapon_range(bool use_far_clip) const {

	assert(weapon >= 0 && weapon <= NUM_WEAPONS);
	float const range(weapons[weapon].range[wmode&1]);
	return ((use_far_clip && range == 0.0) ? FAR_CLIP : range);
}


void player_state::verify_wmode() {

	if (weapon == W_GRENADE && (wmode&1) && p_ammo[weapon] < int(weapons[W_CGRENADE].def_ammo) && !UNLIMITED_WEAPONS) wmode = 0;
}


