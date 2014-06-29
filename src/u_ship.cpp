// 3D World - u_ship implementation for universe mode
// by Frank Gennari
// 11/10/05

#include "ship.h"
#include "ship_util.h"
#include "explosion.h"
#include "openal_wrap.h"
#include "shaders.h"
#include <sstream>


bool const ENABLE_PARTS        = 1; // engine particle trail
bool const GEN_FRAGMENTS       = 1;
bool const SHOW_SHIELDS        = 1;
bool const NO_AI_FIRE          = 0;
bool const AI_MULTI_FIRE       = 0;
bool const MULTI_SEC_FIRE      = 0;
bool const INTER_RELOAD_FIRE   = 0;
bool const PREDICT_TARGETS     = 1; // turret direction
bool const PREDICT_TARGETS2    = 1; // ship direction
bool const PTARG_ALL_WEAPS     = 0; // questionable
bool const OBSTACLE_AVOID      = 1;
bool const TEST_SHIP_BOUNDS    = 0;
int  const COMMON_TARGETS      = 1; // 0 = no target sharing, 1 = only if within sensor dist, 2 = within sensor network

unsigned const SHIELDS_TIME    = unsigned(0.50*TICKS_PER_SECOND);
unsigned const TARGET_CTIME    = unsigned(0.50*TICKS_PER_SECOND);
unsigned const RETARG_DELAY    = unsigned(1.00*TICKS_PER_SECOND);
unsigned const SHIP_AI_DELAY   = unsigned(1.00*TICKS_PER_SECOND);
unsigned const NEUT_CHASE_T    = unsigned(30.0*TICKS_PER_SECOND);
unsigned const FFIRE_WAIT_T    = unsigned(10.0*TICKS_PER_SECOND); // wait time after all enemies are gone before firing on an attacking ally
unsigned const DISABLE_TIME    = unsigned(4.00*TICKS_PER_SECOND);
unsigned const ENG_REPAIR_TIME = unsigned(5.00*TICKS_PER_SECOND);

unsigned const MAX_N_PARTICLES = 80;
unsigned const BASE_NUM_PARTS  = 40;

float const MAX_COLL_DAMAGE    = 1000.0; // limit the damage done by hyperspeed hitting an object
float const MAX_PARTICLE_SIZE  = 0.0004;
float const SHIP_EXP_DAM_SCALE = 0.5;
float const DISABLE_ARMOR      = 0.15;
float const ENGINE_DOWN_ARMOR  = 0.45;
float const WEAPON_DOWN_ARMOR  = 0.4;
float const KILL_CREW_ARMOR    = 0.5;
float const OBS_AVOID_TIME     = 2.5*TICKS_PER_SECOND; // look 2.5 seconds ahead
float const CLOAK_RATE         = 0.5/TICKS_PER_SECOND;
float const SHIP_ROT_ATTEN     = 0.96;
float const ENERGY_XFER_EFF    = 0.5;
float const ORBITING_RD_SCALE  = 2.0;
float const FAST_TARG_DIST     = 1.0; // fast speed if target dist > this value
float const SHIP_GMAX          = 10.0*MAX_SOBJ_GRAVITY;


extern bool player_autopilot, player_auto_stop, player_enemy, regen_uses_credits, respawn_req_hw, hold_fighters, dock_fighters, build_any;
extern int frame_counter, iticks, begin_motion, onscreen_display, display_mode, animate2;
extern float fticks, urm_proj, global_regen, hyperspeed_mult, player_turn_rate, rand_spawn_ship_dmax;
extern unsigned alloced_fobjs[], team_credits[], init_credits[], ind_ships_used[];
extern exp_type_params et_params[];
extern vector<free_obj const *> a_targets, attackers;
extern vector<ship_explosion> exploding;
extern vector<us_class> sclasses;
extern vector<us_weapon> us_weapons;
extern vector<usw_ray> b_wrays;
extern vector<temp_source> temp_sources;
extern vector<unsigned> build_types[];
extern shader_t emissive_shader;


// ************ U_SHIP ************


u_ship::u_ship(unsigned sclass_, point const &pos0, unsigned align, unsigned ai_type_, unsigned target_mode_, bool rand_orient)
			   : free_obj(pos0), u_ship_base(sclass_), ai_type(ai_type_), curr_weapon(0), target_mode(target_mode_)
{
	if (specs().kamikaze) ai_type |= AI_KAMIKAZE;
	flags       = (OBJ_FLAGS_SHIP | OBJ_FLAGS_TARG);
	init_align  = align;
	is_flagship = 0;
	child_stray_dist = 0.0;
	reset();
	if (rand_orient) do_rotate(TWO_PI*rand_float(), TWO_PI*rand_float());
}


u_ship::~u_ship() {
	
	status = 1;
	free_weapons();
	surface_mesh.clear();
}


void u_ship::reset() {

	lhyper       = 0;
	damaged      = 0;
	last_hit     = 0;
	target_set   = 0;
	fire_primary = 0;
	has_obstacle = 0;
	captured     = 0;
	dest_override= 0;
	eflags       = 0;
	alignment    = init_align;
	retarg_time  = 0;
	exp_time     = 0;
	tup_time     = 0;
	disable_t    = 0;
	elapsed_on_t = 0;
	last_targ_t  = FFIRE_WAIT_T; // start off ready to fire on an attacking ally
	kills        = 0;
	docked       = 0;
	o_docked     = 0;
	exp_val      = 0.0;
	cloaked      = 0.0;
	fuel         = 1.0;
	tow_mass     = 0.0;
	used_cargo   = 0.0; // ???
	roll_val     = 0.0;
	pitch_r      = 0.0;
	yaw_r        = 0.0;
	roll_r       = 0.0;
	cached_rsv   = 0.0;
	size_scale   = 1.0;
	check_size_scale();
	shields      = get_max_shields();
	armor        = get_max_armor();
	ncrew        = specs().ncrew;
	tcent        = all_zeros;
	hit_dir      = zero_vector;
	obs_orient   = zero_vector;
	target_dir   = zero_vector;
	time         = rand() % (SHIP_AI_DELAY/2 + 1); // to randomize startups
	if (rand()&1) dest_mgr.clear(); // sometimes clear the destination (and possibly choose a new one) and sometimes keep it
	
	if (!is_orbiting()) {
		if (has_homeworld()) { // respawn at homeworld
			homeworld.update_pos(); // what if this fails?
			reset_pos = homeworld.get_pos() + signed_rand_vector()*homeworld.get_radius();
		}
		else {
			homeworld.clear(); // possibly unnecessary, but just to be safe
		}
	}
	set_max_sf(SLOW_SPEED_FACTOR);
	fighters.clear();
	free_weapons();
	copy_weapons_from_sclass();
	if (weapons.empty()) weapons.push_back(ship_weapon(UWEAP_NONE)); // default, so that weapons isn't empty
	if (player_ship_ptr == NULL || !is_player_ship()) curr_weapon = 0;
	free_obj::reset();
	calc_wpt_center();
}


void u_ship::create_from(u_ship_base const &base) { // related to u_ship_base::create_from()

	ncrew      = base.ncrew;
	shields    = base.shields;
	armor      = base.armor;
	energy     = base.energy;
	fuel       = base.fuel;
	size_scale = base.size_scale;
	ncredits   = base.ncredits;
	kills      = base.kills;
	tot_kills  = base.tot_kills;
	deaths     = base.deaths;
	used_cargo = base.used_cargo;
	weapons    = base.weapons;
	check_size_scale();
}


bool u_ship::regen_enabled() const {
	
	if (is_fighter()) return 0;
	if (respawn_req_hw && !homeworld.is_valid()) return 0;
	unsigned const reserve_credits(sclasses[USC_ARMED_COL].cost + max(sclasses[USC_DEFSAT].cost, sclasses[USC_ANTI_MISS].cost));
	return (!regen_uses_credits || alloc_resources_for(sclass, alignment, reserve_credits, 0.5)); // half price
}


bool u_ship::player_controlled() const {
	
	return (is_player_ship() && !player_autopilot);
}


string u_ship::get_name() const {
	
	return name + " [" + align_names[get_align()] + " " + specs().name + "]";
}


string u_ship::get_info() const { // destination, homeworld?
	
	std::ostringstream oss;
	if (specs().max_shields > 0.0) oss << "  Shields: " << unsigned(100.0*shields/specs().max_shields) << "%";
	if (specs().max_armor   > 0.0) oss << "  Armor: "   << unsigned(100.0*armor/specs().max_armor)     << "%";
	if (specs().ncrew       > 0  ) oss << "  Crew: "    << ncrew;
	if ((ai_type & AI_BASE_TYPE) == AI_ATT_ALL) oss << " (Pirate)";
	if (ai_type & AI_GUARDIAN) oss << " (Guardian)";
	if (is_flagship)           oss << " (Flagship)";
	if (is_exploding())        oss << " (Exploding)";
	else if (disabled_priv())  oss << " (Disabled)";
	if (target_obj)            oss << endl << " (Target: "      << target_obj->get_name() << ")";
	if (parent)                oss << endl << " (Parent: "      << parent->get_name()     << ")";
	if (dest_mgr.is_valid())   oss << endl << " (Destination: " << dest_mgr.get_name()    << ")";
	return oss.str();
}


void u_ship::move_by(point const &pos_) {
	
	tcent += pos_;
	dest_mgr.move_by(pos_);
	homeworld.move_by(pos_);
	free_obj::move_by(pos_);
}


inline void u_ship::no_ai_wait() {

	time = SHIP_AI_DELAY;
}


void u_ship::check_size_scale() {
	
	assert(size_scale > 0.0);
	radius   = size_scale*specs().radius;
	c_radius = specs().cr_scale*radius; // collision radius
}


void u_ship::check_ref_objs() {

	if (status != 0) return;
	u_ship_base::check_ref_objs(this);
	free_obj::check_ref_objs();
}


 // engine status (power) is equal to the number of functioning engines divided by the number of total engines
float u_ship::get_engine_status() const {

	if (specs().nengines == 0 || eflags == 0) return 1.0;
	unsigned const num_engines(specs().nengines);
	float estatus(0.0);

	for (unsigned i = 0; i < num_engines; ++i) {
		if (!(eflags & (1 << i))) estatus += 1.0; 
	}
	return (estatus/num_engines);
}


float u_ship::get_real_speed_val() const {

	float const over_mass_factor(get_mass()/u_ship_base::get_mass());
	float const cscale(min(1.0, (get_crew_scale() + 0.5*SHIP_REQ_CREW))); // slightly slower with fewer crew
	float const trm_scale(min(1.0f, get_true_rel_mass_scale())); // max is 1.0 in case ship has special weapons/ammo (player ship)
	assert(trm_scale > TOLERANCE);
	return get_engine_status()*cscale*over_mass_factor/(0.5 + 0.5*trm_scale); // use trm_scale in other mass calculations?
}


void u_ship::thrust(int tdir, float speed, bool hyperspeed) {

	if (invalid_or_disabled()) return;
	us_class const &sc(specs());
	if (!specs().has_hyper) hyperspeed = 0;

	if (tdir == MOVE_LEFT || tdir == MOVE_RIGHT) {
		roll_r = ((tdir == MOVE_LEFT) ? -1.0 : 1.0)*fticks*speed*sc.roll_rate;

		if (player_controlled()) { // make acceleration-based for the player
			rot_axis = (rot_axis*rot_rate + dir*(0.05*roll_r));
			rot_rate = min(max(rot_rate, fabs(roll_r)), rot_axis.mag());
			rot_axis.normalize();
		}
		else {
			tilt(roll_r); // rotate up vector along current view direction
		}
		return;
	}
	if (tdir == MOVE_FRONT && dot_product(velocity, dir) >= 0.9999*(hyperspeed ? hyperspeed_mult : 1.0)*get_max_speed()) {
		return; // already moving at max speed
	}
	bool stop(0);
	if (cached_rsv == 0.0) cached_rsv = get_real_speed_val();
	speed *= max_sfactor; // *speed_factor?
	if (tdir == MOVE_FRONT) use_fuel((hyperspeed ? 10.0 : 1.0)*speed);
	speed *= cached_rsv;
	if (!hyperspeed && lhyper && velocity.mag() > sc.max_speed) hyperspeed = 1; // can't switch out of hyperspeed until speed is low
	if (speed_factor < 1.0) hyperspeed = 0; // can't do hyperspeed speed if speed is limited
	lhyper = hyperspeed;
	
	if (tdir == MOVE_STOP && !sc.stoppable) {
		if (sc.reversible) {
			float const dp(dot_product(velocity, dir));

			if (dp > TOLERANCE) { // reverse thrust
				stop  = 1;
				tdir  = MOVE_FRONT;
				speed = -speed; // *** really want to turn to <velocity> and also control the throttle ***
			}
			else return;
		}
		else return; // *** turn around to <-velocity> and move forward? ***
	}
	if (tdir == MOVE_BACK) {
		if (sc.reversible) {
			tdir  = MOVE_FRONT;
			speed = -speed;
		}
		else if (sc.stoppable) {
			tdir = MOVE_STOP; // closest we can get to moving backwards
		}
		else return; // *** turn around to <-velocity> and move forward? ***
	}
	if (tdir == MOVE_FRONT || tdir == MOVE_STOP) {
		speed *= fticks;
		if (!sc.has_hyper) hyperspeed = 0;
		if (hyperspeed) speed *= hyperspeed_mult;
	}
	if (tdir == MOVE_FRONT) {
		if (stop && velocity.mag() < fabs(speed)*sc.accel) {
			velocity = zero_vector; // stopped
		}
		else {
			accelerate(speed, sc.accel);
			//set_ship_max_speed();
		}
	}
	else if (tdir == MOVE_STOP) {
		decelerate(speed, sc.decel);
	}
}


void u_ship::turn(vector3d delta) {

	if (invalid_or_disabled()) return;
	float mult(1.0);

	if (player_controlled()) {
		delta *= player_turn_rate; // turn into angular velocity
		mult  *= fticks;
	}
	if (cached_rsv == 0.0) cached_rsv = get_real_speed_val();
	delta *= mult*cached_rsv;
	delta.set_max_mag(specs().max_turn); // limit turning speed
	pitch_r = delta.y;
	yaw_r   = delta.x; // delta.z is unused
	do_rotate(pitch_r, yaw_r);
}


int u_ship::get_move_dir() {

	switch (ai_type & AI_BASE_TYPE) {
		case AI_IGNORE:    // do nothing, don't even move
			break;
		case AI_RETREAT:   // move in direction opposite closest enemy
			return -1;
		case AI_ATT_WAIT:  // attack if target_obj but not closest enemy
		case AI_ATT_ENEMY: // attack closest/all enemies
		case AI_ATT_ALL:   // attack closest ship (rogue)
			return 1;
		default:
			assert(0);
	}
	return 0;
}


vector3d u_ship::get_tot_vel_at(point const &cpos) const { // velocity may be too high in some cases

	vector3d vtot(free_obj::get_tot_vel_at(cpos));
	vtot += calc_angular_vel(cpos, cross_product(dir, upv), pitch_r);
	vtot += calc_angular_vel(cpos, upv, yaw_r);
	vtot += calc_angular_vel(cpos, dir, roll_r);
	assert(!is_nan(vtot));
	return vtot;
}


bool u_ship::do_multi_target() const { // add time delay/rand-mod?

	if (target_mode != TARGET_CLOSEST || !TEAM_ALIGNED(alignment) || !weap_turret(get_weapon_id())) return 0;
	if ((ai_type & AI_BASE_TYPE) == AI_ATT_WAIT) return 0;
	return 1;
}


free_obj const *u_ship::find_closest_target(point const &pos0, float min_dist, float max_dist, bool req_shields) const {

	bool const dir_pref(specs().max_turn > 0.0);

	if ((ai_type & AI_BASE_TYPE) == AI_ATT_ALL || alignment == ALIGN_PIRATE) { // everyone is your enemy
		return get_closest_ship(pos0, min_dist, max_dist, 1, 1, req_shields, 0, dir_pref);
	}
	else { // RETREAT, ENEMY
		assert(alignment < NUM_ALIGNMENT);

		switch (alignment) {
			case ALIGN_NEUTRAL:
			case ALIGN_GOV:
				return NULL;
			case ALIGN_PLAYER:
				if (!player_enemy) return NULL;
			default: // ALIGN_PIRATE, ALIGN_RED, ALIGN_BLUE, etc.
				if (COMMON_TARGETS && (rand()&7) == 0) { // every 8th frame
					free_obj const *friendly(get_closest_ship(pos0, min_dist, max_dist, 0, 0, 0, 0, 0));
					
					if (friendly) { // see if a friendly has chosen a target, and if so, then accept the target as our own
						free_obj const *targ(friendly->get_target());
						
						if (target_valid(targ) && targ != friendly->get_parent() && (!req_shields || targ->has_shields()) &&
							(target_obj == NULL  || p2p_dist(targ->get_pos(), pos) <= min_dist) &&
							(COMMON_TARGETS == 2 || p2p_dist(targ->get_pos(), pos) <= max_dist))
						{
							return targ;
						}
					}
				}
				return get_closest_ship(pos0, min_dist, max_dist, 1, 0, req_shields, 0, dir_pref);
		}
	}
	return NULL;
}


void u_ship::acquire_target(float min_dist) {

	unsigned const ai_base_type(ai_type & AI_BASE_TYPE);
	float const tdist((target_obj == NULL) ? 0.0 : p2p_dist(pos, target_obj->get_pos()));
	float search_dist(specs().sensor_dist);

	if (!can_move() && fighters.empty()) { // if can't move, then there is no point to acquiring a target out of weapons range
		float const weap_range(specs().get_weap_range());
		if (weap_range > 0.0) search_dist = min(search_dist, (1.1f*weap_range + c_radius));
	}
	if (target_obj != NULL && (target_obj->is_resetting() || target_obj->is_invisible() ||
		(COMMON_TARGETS < 2 && tdist > search_dist)))
	{
		target_obj = NULL; // don't target a ship that's out of sensor range or already dead
	}

	// RETREAT, WAIT, ENEMY, ALL
	if (ai_base_type != AI_ATT_WAIT) { // RETREAT, ENEMY, ALL
		if (target_obj == NULL || target_obj == parent || target_obj->invalid() ||
			time > (tup_time + TARGET_CTIME) || tdist > search_dist || tdist < min_dist)
		{
			tup_time = time + ((rand()%TARGET_CTIME) >> 1);  // update target every so often, randomize
			free_obj const *new_target_obj(NULL);
			bool find_closest(0);

			switch (target_mode) {
			case TARGET_CLOSEST:
				find_closest = (target_obj == NULL || retarg_time == 0);
				break;
			case TARGET_ATTACKER:
			case TARGET_LAST:
				find_closest = (target_obj == NULL);
				break;
			case TARGET_PARENT:
				if (parent != NULL && target_valid(parent->get_target()) && !parent->get_target()->is_invisible()) {
					target_obj = parent->get_target();
				}
				else {
					find_closest = 1;
				}
				break;
			default:
				assert(0);
			}
			bool const has_dest(dest_mgr.is_valid());
			if (has_dest && (rand()&3)) find_closest = 0; // every 4th frame if already have a destination
			
			if (find_closest) {
				if (target_obj != NULL) {
					if (alignment == ALIGN_NEUTRAL && ai_base_type == AI_ATT_ENEMY && (rand() % NEUT_CHASE_T) == 0) {
						target_obj = NULL; // give up the chase after awhile
					}
					if (tdist > 2.0*search_dist) target_obj = NULL; // (tdist < min_dist) is ignored for now, out of range
				}
				float eff_search_dist(search_dist);
				if (has_dest) eff_search_dist = min(search_dist, p2p_dist(pos, dest_mgr.get_pos()));
				if (target_obj != NULL && tdist >= min_dist) eff_search_dist = min(search_dist, 0.8f*tdist);
				new_target_obj = find_closest_target(pos, min_dist, eff_search_dist, 0);
				if (new_target_obj == NULL) new_target_obj = target_obj; // keep the same target

				if (new_target_obj == NULL && alignment != ALIGN_NEUTRAL) { // no target, choose to attack same target as teammates
					assert(alignment < a_targets.size());
					
					if (target_valid(a_targets[alignment]) && dist_less_than(pos, a_targets[alignment]->get_pos(), search_dist)) {
						new_target_obj = a_targets[alignment];
					}
				}
				if ((ai_type & AI_GUARDIAN) && new_target_obj == NULL) { // seek out the last attacker
					unsigned const start_i(rand() % NUM_ALIGNMENT); // don't show favoritism
					
					for (unsigned i = 0; i < NUM_ALIGNMENT; ++i) {
						unsigned const ii((start_i + i) % NUM_ALIGNMENT);
						assert(alignment < attackers.size());
						free_obj const *att(attackers[ii]);
						
						if (target_valid(att) && att->get_align() != ALIGN_GOV && att->get_align() != alignment &&
							dist_less_than(pos, att->get_pos(), search_dist)) // what about targets related to this target?
						{
							new_target_obj = att; // don't attack gov or friend ships
							break;
						}
					}
				}
			}
			if (new_target_obj != NULL) {
				target_obj  = new_target_obj;
				retarg_time = RETARG_DELAY;
				
				if (target_obj->is_player_ship() && target_obj != parent) {
					send_warning_message(string("Enemy Ship Detected: ") + get_name());
				}
			}
			if (target_obj != NULL && target_mode == TARGET_LAST) target_set = 1;
		}
	}
	if ((ai_type & AI_GUARDIAN) && target_obj != NULL && target_obj->get_align() == alignment) {
		target_obj = NULL; // don't attack a friendly
	}
	if (target_obj == NULL && parent != NULL && target_valid(parent->get_target())) {
		target_obj = parent->get_target(); // as a last resort, even if not TARGET_PARENT
	}
	if (!fighters.empty()) get_fighter_target(this);
	
	if (target_obj != NULL && target_obj != parent) {
		if (specs().for_boarding && !target_obj->can_board()) target_obj = NULL; // can't board this ship
		else if (target_obj->is_invisible())                  target_obj = NULL; // invisible (cloaked ship)
	}
	assert(target_obj != this);
}


uobject const *u_ship::setup_int_query(vector3d const &qdir, float qdist, free_obj *&fobj,
									   float &tdist, bool sobjs_only, float line_radius) const
{
	assert(qdist > 0.0);
	int qtype(get_line_query_obj_types(qdist));
	if (sobjs_only) qtype &= OBJ_TYPE_SOBJ;
	if (qtype == 0) return NULL;
	line_int_data li_data(pos, qdir, qdist, this, target_obj, 1, 0, line_radius);
	uobject const *obj(line_intersect_objects(li_data, fobj, qtype));
	tdist = li_data.dist;
	return obj;
}


void u_ship::calc_wpt_center() {

	wpt_center = all_zeros;

	for (unsigned i = 0; i < weapons.size(); ++i) {
		for (unsigned j = 0; j < weapons[i].weap_pts.size(); ++j) {
			wpt_center += weapons[i].weap_pts[j];
		}
	}
}


uobject const *u_ship::get_obstacle(float oa_time, float max_dist, free_obj *&fobj, float &tdist, bool sobjs_only) const {

	if (oa_time == 0.0) return NULL;
	float const vmag(velocity.mag());
	if (vmag < max(TOLERANCE, 0.1f*get_max_speed())) return NULL;
	float ldist(oa_time*min(2.5f, fticks)*vmag); // see if there is an obstacle directly ahead
	if (max_dist > c_radius) ldist = min(ldist, max_dist);
	return setup_int_query(velocity/vmag, (ldist + c_radius), fobj, tdist, sobjs_only, c_radius);
}


// this is really slow
bool u_ship::obstacle_avoid(vector3d &orient, float target_dist, bool sobjs_only) { // find an obstacle other than target_obj

	if (!OBSTACLE_AVOID || (flags & OBJ_FLAGS_DIST)) return 0;

	if ((time&3) != (sclass&3)) { // only check every 4th frame (for efficiency)
		if (has_obstacle) {
			assert(obs_orient != zero_vector);
			orient = obs_orient;
		}
		return has_obstacle;
	}
	has_obstacle = 0;
	float tdist(0.0);
	free_obj *fobj = NULL;
	uobject const *obstacle(get_obstacle(OBS_AVOID_TIME, target_dist, fobj, tdist, sobjs_only)); // sort of incomplete
	if (obstacle == NULL) obstacle = setup_int_query(orient, target_dist, fobj, tdist, sobjs_only, c_radius); // now see if there is an obstacle between us and our target
	// really need to check if the trajectories of any close ships will cause a future intersection
	// trajectory (cylinder-cylinder) intersections are difficult and slow though, and can't account for direction change
	//if (obstacle == NULL) obstacle = get_closest_ship(pos, c_radius, target_dist, 0, 1, 0, 0); // dir_pref?
	
	if (obstacle != NULL) { // obstacle in the way
		if (fobj != NULL && fobj->get_parent() == this && fobj->get_target() == this) return 0; // docking fighter
		
		if (velocity != zero_vector) {
			float const avoid_radius(min(0.95*tdist, (obstacle->get_bounding_radius() + 1.5*c_radius))); // make sure it works if close
			float const angle(asinf(avoid_radius/tdist)); // -angle?
			vector3d const v2obj(pos, obstacle->get_pos()), vrot(cross_product(v2obj, velocity).get_norm());
			if (vrot != zero_vector) rotate_vector3d_norm(vrot, angle, orient); // more complex than this - what if orient != vnorm?
		}
		obs_orient   = orient;
		has_obstacle = 1;
		return 1;
	}
	return 0;
}


ship_explosion u_ship::get_explosion() const {

	float const size(get_def_explode_damage());
	return ship_explosion(size*radius, EXPLOSION_DAMAGE*size*radius, pos, sclasses[sclass].exp_disint);
}


bool u_ship::avoid_explosions(vector3d &orient) const {

	float const min_damage(1.2*specs().damage_abs + 1.0);
	float max_damage(min_damage);

	for (unsigned i = 0; i < exploding.size(); ++i) { // inefficient
		vector3d const dist(pos, exploding[i].pos);
		float const dsq(dist.mag_sq()), br(exploding[i].bradius);
		if (dsq > br*br) continue;
		float const damage(exploding[i].damage*calc_damage_scale(sqrt(dsq), radius, br));

		if (damage > max_damage) {
			max_damage = damage;
			orient     = dist; // point away from the exploding ship
		}
	}
	if (max_damage > min_damage) {
		if (orient.mag() < TOLERANCE) orient = signed_rand_vector();
		orient.normalize();
		return 1;
	}
	return 0;
}


void u_ship::do_turn(vector3d const &orient) {

	float const turn_angle(fabs(get_angle(orient, dir)));

	if (turn_angle > 0.001) {
		calc_rotation_vectors();
		vector3d orient2(orient);
		rotate_vector3d     (rv1, ra1, orient2);
		rotate_vector3d_norm(rv2, ra2, orient2);
		vector3d const delta(-orient2.y, -orient2.x, 0.0); // AI doesn't have that problem with screen Y-values being backwards
		float const dmag(delta.mag());
		if (dmag > TOLERANCE) turn(delta*(turn_angle/dmag)); // delta is the actual turn angle broken up into x and y components
	}
}


bool u_ship::can_return_to_parent() const {
	
	if (parent == NULL || parent->invalid() || parent->disabled()) return 0; // no parent to return to or parent is disabled
	if (is_orbiting()) return 0;
	if (size_scale > 4.0 || (size_scale > 1.0 && get_mass() > 0.1*parent->get_mass())) return 0; // too large to return to parent

	// don't dock until fighters that can return have returned, otherwise they will be assigned to our parent
	for (set<u_ship *>::const_iterator i = fighters.begin(); i != fighters.end(); ++i) { // works for boarding shuttles as well
		if (has_space_for_fighter((*i)->sclass)) return 0; // must wait for this fighter to return
	}
	if (!parent->get_ship_base()->has_space_for_fighter(sclass)) return 0; // check if parent has space in its fighter bay
	return 1;
}


bool u_ship::check_return_to_parent() const {

	if (last_hit > 0)            return 0; // under attack
	if (!can_return_to_parent()) return 0;
	if (dock_fighters && parent != NULL && parent->is_player_ship()) return 1; // dock with parent
	if (specs().for_boarding && ncrew <= specs().ncrew/2)            return 1; // can no longer board
	u_ship_base const *parent_ship(parent->get_ship_base());
	assert(parent_ship);
	bool const dock_allowed(player_autopilot || !parent->is_player_ship());

	if (target_obj == NULL) { // no targets
		if (dock_allowed || !dist_less_than(pos, parent->get_pos(), (1.2*parent->get_c_radius() + 4.0*c_radius))) {
			return 1; // avoid automatically docking with player
		}
	}
	else if (is_fighter() && specs().max_speed > TOLERANCE) { // need to dock to regen weapons/shields/armor?
		assert(parent->is_ship());
		if (!dock_allowed) return 0; // not sure what to do about this case
		if (eflags != 0)   return 1; // engine(s) need rapaired
		float dock_value(0.0);

		for (unsigned i = 0; i < weapons.size(); ++i) {
			ship_weapon const &w(weapons[i]);

			if (w.get_usw().need_ammo() && w.ammo < w.init_ammo && parent_ship->has_ammo_for(w.wclass)) {
				dock_value += 1.0*(1.0 - ((float)w.ammo)/((float)w.init_ammo)); // weapon ammo bonus
			}
			dock_value += 1.5*w.ndamaged; // damaged weapon bonus
		}
		dock_value += 2.0*get_damage(); // recharge bonus
		dock_value -= 0.25*p2p_dist(pos, parent->get_pos())/specs().max_speed;
		if (dock_value > 1.0) return 1;
	}
	return 0;
}


bool sobj_manager::claim_object(free_obj *parent, bool homeworld) {

	return check_dest_ownership(uobj_id, pos, parent, 0, homeworld);
}


void sobj_manager::set_object(uobject const *obj) {

	assert(obj != NULL);
	pos     = obj->get_pos();
	radius  = obj->get_bounding_radius();
	uobj_id = obj->get_id();
	otype   = obj->get_type();
	owner   = obj->get_owner();
	name    = obj->get_name();
}


void sobj_manager::choose_dest(point const &p, unsigned align, float tmax) {

	uobject const *dest(choose_dest_world(p, old_uobj_id, align, tmax));

	if (dest != NULL) {
		set_object(dest);
		assert(uobj_id != old_uobj_id);
	}
}


bool sobj_manager::update_pos_if_close(uobject const *obj) {

	if (obj == NULL || obj->get_id() != uobj_id) return 0;
	point const new_pos(obj->get_pos());
	if (p2p_dist(pos, new_pos) > 0.1*CELL_SIZE)  return 0; // don't update if in a different cell, let shift handle it
	pos   = new_pos;
	owner = obj->get_owner();
	return 1;
}


bool sobj_manager::update_pos() {

	return (is_valid() && update_pos_if_close(get_closest_world_ptr(pos, otype)));
}


bool u_ship::choose_destination() {

	if ((is_fighter() && parent != NULL) || is_orbiting()) return 0;

	if (dest_mgr.is_valid()) {
		float const rad(c_radius + dest_mgr.get_radius());

		if (powered_priv() && dist_less_than(pos, dest_mgr.get_pos(), 1.6f*rad)) {
			claim_world(NULL); // or at least try to
			dest_mgr.at_dest();
		}
		else {
			bool renew_dest((rand()&15) == 0 && p2p_dist(pos, dest_mgr.get_pos()) > 0.2*CELL_SIZE);

			if (!renew_dest && dist_less_than(pos, dest_mgr.get_pos(), specs().sensor_dist) && !(rand()&15)) {
				if (!dest_mgr.update_pos()) { // update every 16 frames in case planet has moved
					if (dist_less_than(pos, dest_mgr.get_pos(), 5.0f*rad)) renew_dest = 1;
				}
			}
			if (!renew_dest) return 1; // already chosen - change it?
		}
	}
	dest_mgr.choose_dest(pos, alignment, 1.5*get_max_t()/FOBJ_TEMP_SCALE); // allow higher than tmax, since ships can approach in the shade of the planet/moon
	
	if (!dest_mgr.is_valid()) { // delete the ship? move towards the player?
		if (specs().has_fast_speed) set_max_sf(FAST_SPEED_FACTOR);
		thrust(MOVE_FRONT, 1.0, 0); // full speed ahead
		return 0;
	}
	return 1;
}


bool u_ship::claim_world(uobject const *uobj) { // doesn't do a whole lot yet

	if (!is_player_ship() && !dest_mgr.claim_object(this, !has_homeworld())) {return 0;} // already claimed

	if (!has_homeworld()) { // make this our homeworld
		if (uobj == NULL) {
			homeworld = dest_mgr; // use the current destination as a homeworld
		}
		else {
			homeworld.set_object(uobj);
		}
	}
	return 1;
}


u_ship const *u_ship::try_fighter_pickup() const {

	u_ship const *cur_targ(NULL);
	float fdist_sq(0.0);

	for (set<u_ship *>::const_iterator i = fighters.begin(); i != fighters.end(); ++i) {
		u_ship const *const f(*i);
		assert(f);
		if (f->invalid_or_disabled()) continue; // disabled ships can't dock anyway, so wait until they become un-disabled
		if (f->get_parent() != this || f->get_target() != this) continue; // not ready to dock
		if (f->specs().max_turn > 0.0 && specs().max_speed <= f->specs().max_speed) continue;
		if (!has_space_for_fighter(f->sclass)) continue; // no space
		float const fd_sq(p2p_dist_sq(pos, f->get_pos()));

		if (fdist_sq == 0.0 || fd_sq < fdist_sq) { // pickup this fighter
			cur_targ = f;
			fdist_sq = fd_sq;
		}
	}
	return cur_targ;
}


free_obj const *u_ship::try_orbital_regen(free_obj const *cur_targ, bool last_od, bool &targ_friend, bool &o_dock_close) {

	bool regen(cur_targ == NULL);
	float max_dist(specs().sensor_dist);

	if (!regen) {
		float const target_dist(p2p_dist(cur_targ->get_pos(), pos)), wrange(specs().get_weap_range());
		float min_targ_dist(wrange);
		if (last_od || get_damage() > 0.75) min_targ_dist = min(min_targ_dist, 10.0f*radius);
		
		if (target_dist > min_targ_dist) {
			regen = 1;
			if (target_dist < wrange) max_dist = min(max_dist, 0.25f*target_dist);
		}
	}
	if (regen) {
		free_obj *dock(get_closest_dock(max_dist)); // damaged, look for a dock (what if have fighters?)
		
		if (dock) {
			float const time_to_target(min_time_to_target(dock->get_pos()));
			float const damage_at_target(get_damage_after_time(time_to_target));

			if (last_od || damage_at_target > 0.25) {
				float const dock_dist(c_radius + dock->get_c_radius());
				cur_targ    = dock;
				targ_friend = 1;

				if (dist_less_than(pos, dock->get_pos(), 2.3*dock_dist)) {
					dock->orbital_dock(this);
					o_dock_close = dist_less_than(pos, dock->get_pos(), 1.7*dock_dist);
					o_docked     = 1;
				}
			}
		}
	}
	return cur_targ;
}


bool u_ship::roll_to_face_target(float &roll_amt) const {

	assert(target_dir != zero_vector);
	if (fabs(wpt_center.x) < 0.01 && fabs(wpt_center.y) < 0.01) return 0;
	vector3d tdir(target_dir);
	rotate_point(tdir);
	tdir.z = 0.0;
	float const td_mag(tdir.mag());

	if (td_mag > 0.01) {
		tdir /= td_mag;
		point wcenter(wpt_center);
		wcenter.z = 0.0;
		wcenter.normalize();
		float const cpz(cross_product(wcenter, tdir).z), angle(safe_acosf(dot_product(wcenter, tdir)));

		if (fabs(cpz) > 0.01 && fabs(angle) > 0.01) { // roll to face weapons toward target
			roll_amt = CLIP_TO_pm1(((cpz < 0.0) ? -1 : 1)*angle);
		}
	}
	return 1;
}


float u_ship::get_fast_target_dist(free_obj const *const target) const {

	float ftd(FAST_TARG_DIST);
	if (specs().weap_range > 0.0) ftd = min(ftd, specs().weap_range);

	if (target != NULL && target->get_velocity().mag() > specs().max_speed*SLOW_SPEED_FACTOR) { // target is at high speed
		if (specs().fire_speed != FSPEED_SLOW) ftd *= 0.5; // return 0.0?
		if (!weapons.empty() && us_weapons[get_weapon_id()].hyper_fire) ftd *= 0.5; // return 0.0?
	}
	return ftd;
}


bool u_ship::has_slow_fighters() const {

	for (set<u_ship *>::const_iterator i = fighters.begin(); i != fighters.end(); ++i) {
		if (!(*i)->invalid() && (*i)->get_parent() == this && !(*i)->specs().has_fast_speed) return 1;
	}
	return 0;
}


void u_ship::fire_at_target(free_obj const *const targ_obj, float min_dist) {

	if (!target_valid(targ_obj)) return;
	free_obj const *const orig_tobj(target_obj);
	target_obj = targ_obj;
	upos_point_type const &tpos(targ_obj->get_pos());
	target_dir = (tpos - pos).get_norm();
	ai_fire(target_dir, p2p_dist(tpos, pos), min_dist, 0);
	target_obj = orig_tobj; // restore original value
}


void u_ship::ai_action() {

	float const old_cloaked(cloaked);
	cloaked = 0.0;
	if (time < SHIP_AI_DELAY || invalid_or_disabled()) return;
	if (!player_controlled()) set_max_sf(SLOW_SPEED_FACTOR);
	bool const player_ship(is_player_ship());
	
	if (begin_motion || player_ship) {
		if (rot_rate != 0.0) rot_rate *= pow(SHIP_ROT_ATTEN, fticks); // stabilize (should this vary per ship class?)
		fire_point_defenses(); // always, even player's ship
	}
	if (player_controlled() && specs().stoppable) {
		if (player_auto_stop) thrust(MOVE_STOP, 1.0, 0);
		return;
	}
	if (!begin_motion) return;
	us_class const &sc(specs());
	bool const can_move_(can_move()), last_od(o_docked);
	float const max_turn(sc.max_turn);
	if (sc.roll_rate < 0.0 && !player_ship) thrust(MOVE_LEFT, 1.0, 0); // negative roll implies always rolling
	o_docked = 0;

	// too close to the sun or a hot object, or too much gravity from a black hole, move away at full speed
	if (can_move_) {
		vector3d orient(zero_vector);

		if (is_burning()) {
			if (get_over_temp_factor() > 0.0) { // actually taking damage
				orient = pos - tcent; // fly directly away from star
			}
			else { // fly on a tangent
				orthogonalize_dir(dir, (pos - tcent), orient, 0);
			}
		}
		else if (near_b_hole && gvect.mag() > SHIP_GMAX) { // gravity too high
			orient = -gvect; // fly directly away from black hole
		}
		if (orient != zero_vector) {
			orient.normalize();
			if (max_turn > TOLERANCE) do_turn(orient); // a little unstable
			if (dot_product(dir, orient) > 0.0) thrust(MOVE_FRONT, 1.0, 0);
			return;
		}
	}
	int move_dir(get_move_dir());
	if (!move_dir) return;

	if (hold_fighters && parent != NULL && parent->is_player_ship()) {
		thrust(MOVE_STOP, 1.0, 0); // hold at this position
		return;
	}
	bool const no_ammo(out_of_ammo(0)), boarding(sc.for_boarding && ncrew > sc.ncrew/2), kamikaze((ai_type & AI_KAMIKAZE) != 0);
	if (no_ammo && !kamikaze && !boarding && target_obj != parent) move_dir = -1; // out of ammo, run away
	vector3d avoid_orient(dir);
	float const min_attack(get_min_att_dist()), vmag(velocity.mag());
	float const min_dist((no_ammo || kamikaze || boarding) ? 0.0 : min_attack); // ram the enemy
	bool const avoid_exp(can_move_ && avoid_explosions(avoid_orient)), local_dest(dest_override);
	dest_override = 0;
	
	if (!is_orbiting() || (time&3) == 0) { // every 4th frame if orbiting
		acquire_target(min_dist); // slow
	}
	free_obj const *const acquired_target(target_obj);
	if (local_dest) target_obj = NULL;
	bool const parent_is_player(parent && !player_autopilot && parent->is_player_ship());
	bool const use_stray(parent && !parent_is_player && !parent->disabled() && parent->can_move());
	float stray_dist(0.0);
	
	if (use_stray) {
		float const psd(parent->get_child_stray_dist()), fsd(is_fighter() ? specs().stray_dist : 0.0f);
		stray_dist = ((psd > 0.0 && fsd > 0.0) ? min(psd, fsd) : max(psd, fsd));
	}
	bool const has_strayed(stray_dist > 0.0 && !dist_less_than(pos, parent->get_pos(), 2.0*stray_dist));

	if (check_return_to_parent() || has_strayed) {
		assert(parent != NULL && parent != this);
		target_obj = parent;
		
		if (sc.stoppable && dist_less_than(pos, parent->get_pos(), 4.0*TICKS_PER_SECOND*fticks*velocity.mag())) {
			thrust(MOVE_STOP, (is_fighter() ? 0.12 : 0.6), 0); // return to parent slowly
		}
	}
	free_obj const *cur_targ(target_obj);
	bool targ_friend(cur_targ != NULL && cur_targ == parent), o_dock_close(0);

	if (!local_dest && can_move_ && !last_hit && (!is_fighter() || !parent) && get_damage() > (last_od ? 0.05 : 0.5)) {
		cur_targ = try_orbital_regen(cur_targ, last_od, targ_friend, o_dock_close);
	}
	bool const has_target(cur_targ != NULL && !local_dest);

	if (acquired_target && !no_ammo && (!has_target || targ_friend || local_dest)) {
		fire_at_target(acquired_target, min_dist); // fire at a target while moving towards a destination or returning to parent
	}
	if (!has_target) { // no targets
		if (avoid_exp) {
			do_turn(avoid_orient);
			thrust(MOVE_FRONT, 1.0, 0); // move forward very quickly
		}
		else if (can_move_ && (local_dest || fighters.empty()) && choose_destination()) {
			// have a destination and not waiting for fighters
			if (sc.roll_rate > 0.0 && !player_ship) thrust(MOVE_LEFT, 0.25, 0);
			point dest(dest_mgr.get_pos());
			vector3d orient(dest, pos);
			float const dist(orient.mag()), min_dest_dist(c_radius + dest_mgr.get_radius());
			bool const use_high_speed(dist > 5.0*min_dest_dist && specs().has_fast_speed);
			orient /= dist;
			
			if (can_move_ && vmag > 0.1*get_max_speed()) { // no hope if travelling at high speed
				obstacle_avoid(orient, min(dist, 4.0f*TICKS_PER_SECOND*fticks*vmag), use_high_speed); // 4.0s lookahead
			}
			if (max_turn > TOLERANCE) do_turn(orient);

			if (dist < 1.5*min_dest_dist) {
				thrust(MOVE_STOP, 1.0, 0); // full stop
			}
			else {
				if (use_high_speed && dot_product(orient, dir) > 0.0) set_max_sf(FAST_SPEED_FACTOR);
				thrust(MOVE_FRONT, 1.0, 0); // full speed ahead
			}
		}
		else {
			if (can_move_) {
				cur_targ = try_fighter_pickup();
				if (cur_targ != NULL) targ_friend = 1;
			}
			if (cur_targ == NULL && sc.decel > 0.0) {
				bool const use_high_speed(max_sfactor >= FAST_SPEED_FACTOR); // no hope of avoiding a collision
				float tdist; // value unused
				free_obj *fobj(NULL);

				if (vmag > TOLERANCE) {
					vector3d const vnorm(velocity/vmag);
					if (p2p_dist(vnorm, dir) > 1.0E-6) do_turn(vnorm); // align ourselves with our velocity
				}
				if (parent != NULL || get_obstacle(0.75*OBS_AVOID_TIME, 0.0, fobj, tdist, use_high_speed) != NULL) {
					thrust(MOVE_STOP, 1.0, 0); // hard stop to avoid a collision
				}
				else {
					thrust(MOVE_STOP, 0.1, 0); // slow down slowly - keep from being slowly pulled into a gravitational object
				}
			}
		}
		cloaked = 0.0;
		if (cur_targ == NULL) return;
	}
	if (sc.has_cloak && has_target && (has_strayed || !targ_friend)) {
		cloaked = min(1.0f, (old_cloaked + fticks*CLOAK_RATE));
	}

	// determine target_dir
	assert(cur_targ != NULL && cur_targ != this);
	upos_point_type target(cur_targ->get_pos());
	if (o_docked) target += cur_targ->get_dir()*c_radius; // move away from a planet when docking
	target_dir = target - pos;
	float const target_dist((float)p2p_dist(point_d(target), point_d(pos))); // need more precision
	float const target_radius(cur_targ->get_c_radius());
	if (target_dist < TOLERANCE) {target_dir = zero_vector; return;} // what to do?
	target_dir /= target_dist;
	vector3d fire_dir(target_dir);

	// has some stability problems as velocity is not constant
	if (PREDICT_TARGETS2 && cur_targ != NULL && move_dir == 1 && max_turn > TOLERANCE && !last_od) {
		if (boarding || targ_friend || kamikaze) { // want a collision
			vector3d const tdir(predict_target_dir(pos, cur_targ));
			if (tdir != zero_vector) target_dir = tdir;
		}
		else if (has_target) {
			bool lead_shot(can_lead_shot_with(curr_weapon));
			unsigned targ_weap(curr_weapon);

			if (PTARG_ALL_WEAPS) {
				for (unsigned i = 0; i < weapons.size() && !lead_shot; ++i) { // see if can lead shots for any weapon
					lead_shot = can_lead_shot_with(i);
					targ_weap = i;
				}
			}
			if (lead_shot) {
				unsigned const wid(weapons[targ_weap].wclass);
				bool const s_turret(weap_turret(wid)); // curr_weapon can change
				vector3d const tdir(predict_target_dir(pos, cur_targ, (s_turret ? UWEAP_NONE : wid)));
				
				if (tdir != zero_vector && get_angle(target_dir, tdir) < 0.4) { // only lead shot if in a similar direction
					target_dir = tdir; // can catch the target
					if (!s_turret) fire_dir = target_dir;
				}
			}
		}
	}
	float const mappd(radius*sc.min_app_dist);
	float const min_target_dist(max(mappd, (target_radius + (1.0f + 0.2f*max(min_attack, mappd)/radius)*c_radius)));
	vector3d orient(target_dir);

	if (avoid_exp) {
		orient = avoid_orient;
	}
	else if (move_dir == -1) {
		orient.negate(); // run away (RETREAT)
	}
	else if (has_target && target_dist < min_target_dist && !targ_friend && !kamikaze && !specs().for_boarding) {
		if (get_max_speed() > 1.1*cur_targ->get_velocity().mag()) { // if you can get away
			orient.negate(); // back away, should slow on approach so as not to ram it, unless docking, boarding, etc.
		}
	}

	// turn
	bool obstacle(0);

	if (max_turn > TOLERANCE) { // otherwise dir can't change
		if (can_move_ && !avoid_exp) obstacle = obstacle_avoid(orient, target_dist, 0);
		do_turn(orient);
	}

	// roll
	if (sc.roll_rate > 0.0 && !o_dock_close) {
		float roll_amt(0.0);
		bool roll_set(has_target && !targ_friend && roll_to_face_target(roll_amt));
		
		if (!player_ship && !roll_set) {
			roll_val += fticks*rand_uniform(-0.05, 0.05); // roll is random
			roll_val  = CLIP_TO_pm1(roll_val);
			roll_amt  = roll_val;
		}
		if (roll_amt != 0.0) thrust(((roll_amt < 0.0) ? MOVE_LEFT : MOVE_RIGHT), fabs(roll_amt), 0);
	}

	// move
	if (specs().has_fast_speed && target_dist > get_fast_target_dist(cur_targ) &&
		dot_product(orient, dir) > 0.0 && !has_slow_fighters()) {
		set_max_sf(FAST_SPEED_FACTOR); // fast speed when far from the target
	}
	vector3d const delta_v(velocity - cur_targ->get_velocity());
	float const d_turn_ang_mv((!targ_friend && weap_turret(get_weapon_id())) ? 0.0 : get_angle(orient, dir));
	float const dist_per_sec(TICKS_PER_SECOND*fticks*delta_v.mag());
	bool const pickup_fighter(!has_target && cur_targ != NULL);
	bool const dock(!parent_is_player && ((is_fighter() && cur_targ == parent) || pickup_fighter));
	bool const very_close(target_dist < 1.3*min_target_dist || target_dist < 1.0f*dist_per_sec);
	bool const target_close(very_close || target_dist < 4.5f*dist_per_sec);
	if (dock) cloaked = 0.0;

	// multiple calls to thrust() are illegal since the max of the sums of all of the calls is not checked
	if (avoid_exp) {
		thrust(MOVE_FRONT, 1.0, 0); // move forward very quickly
	}
	else if (obstacle) {
		thrust(MOVE_FRONT, 0.8, 0); // move forward relatively quickly
	}
	else if (can_move_ && target_close && (dock || boarding || kamikaze)) {
		if (d_turn_ang_mv < 1.0 && (cur_targ != parent || can_return_to_parent())) {
			thrust(MOVE_FRONT, 0.8*(1.0 - d_turn_ang_mv), 0);
		}
		else {
			thrust(MOVE_STOP, 1.0, 0); // hard stop to slow down for docking/return
		}
	}
	else if (d_turn_ang_mv < 1.0 && can_move_) { // move forward
		if (o_dock_close) {
			thrust(MOVE_STOP, 1.0, 0); // orbital docking, full stop
		}
		else if (stray_dist > 0.0 && !dist_less_than(pos, parent->get_pos(), stray_dist) &&
			vmag > parent->get_velocity().mag() && dot_product(vector3d(pos - parent->get_pos()), parent->get_velocity()) > 0.0)
		{
			thrust(MOVE_STOP, 0.8, 0); // fairly hard stop
		}
		else if (!no_ammo && !kamikaze && move_dir == 1 && target_dist >= min_target_dist && very_close &&
			dot_product_ptv(delta_v, point(target), point(pos)) > 0.0)
		{ // could be better
			float const decel((cur_targ->radius < 0.75*radius) ? 0.04 : 0.4); // smaller target, OK to ram it (use mass?)
			thrust(MOVE_BACK, decel, 0); // better slow down/back up to avoid a collision
		}
		else { // what about getting too far from parent?
			thrust(MOVE_FRONT, (1.0 - d_turn_ang_mv), 0);
		}
	}
	else {
		thrust(MOVE_STOP, 0.2, 0); // slow down more
	}

	// fire
	if (!no_ammo && !pickup_fighter && !targ_friend) ai_fire(fire_dir, target_dist, min_dist, move_dir);
}


void u_ship::fire_point_defenses() { // point defense weapon - fires at incoming projectiles only

	if (urm_proj == 0.0) return; // no projectiles
	unsigned const curr(curr_weapon);

	for (unsigned i = 0; i < weapons.size(); ++i) {
		assert(weapons[i].wclass < us_weapons.size());
		us_weapon const &weap(us_weapons[weapons[i].wclass]);
		if (!weap.point_def || !check_fire_delay(i) || !check_fire_speed()) continue;
		free_obj const *inc_proj(check_for_incoming_proj(pos, alignment, weap.range));
		if (!inc_proj) continue;
		
		if (target_obj != NULL && target_obj != parent && inc_proj->damage_done() < weap.damage &&
			dist_less_than(target_obj->get_pos(), pos, weap.range))
		{
			continue; // better to shoot at the target than the projectile
		}
		vector3d const fdir(inc_proj->get_pos(), pos);
		float const target_dist(fdir.mag());
		if (target_dist < TOLERANCE) continue;
		curr_weapon = i;
		fire_weapon(fdir/target_dist, target_dist);
	}
	curr_weapon = curr;
}


bool u_ship::find_coll_enemy_proj(float dmax, point &p_int) const {

	free_obj const *const inc_proj(check_for_incoming_proj(pos, alignment, dmax));
	if (!inc_proj) return 0;
	vector3d const rvel((inc_proj->get_velocity() - velocity).get_norm());
	point const ppos(inc_proj->get_pos()), ppos2(ppos + rvel*(2.0*dmax));
	return (line_sphere_int(rvel, ppos, pos, c_radius, p_int, 1) && line_int_obj(ppos, ppos2));
}


void u_ship::ai_fire(vector3d const &targ_dir, float target_dist, float min_dist, int move_dir) {

	assert(target_obj != NULL);
	if (NO_AI_FIRE || time < unsigned(1.5*SHIP_AI_DELAY)) return;
	bool s_turret(weap_turret(get_weapon_id()));
	if (move_dir != 1 && !s_turret) return; // no attack (should s_turret be checked?)
	unsigned const nweap((unsigned)weapons.size());
	if (nweap == 0 || (nweap == 1 && weapons[0].wclass == UWEAP_NONE)) return; // no weapons
	if (target_dist > specs().fire_dist || target_dist < 0.5*min_dist) return; // target not in range
	assert(nweap > 0);
	float const dv(p2p_dist(velocity, target_obj->get_velocity()));
	float const cur_damage(get_damage()), d_angle(get_angle(targ_dir, dir));
	unsigned const target_crew(target_obj->get_ncrew());
	unsigned default_next_weap(curr_weapon);
	static vector<unsigned> good_weapons;
	good_weapons.resize(0);
	
	// multiple weapons fire in the same frame - player's ship works differently (primary/secondary fire) ?
	for (unsigned i = 0; i < nweap; ++i) { // eventually move the code below into here
		unsigned const w((curr_weapon + i) % nweap);
		ship_weapon const &sw(weapons[w]);
		if (sw.wclass == UWEAP_NONE || sw.wcount == 0) continue; // no weapon
		if (sw.no_ammo()) continue; // have no ammo for this weapon
		us_weapon const &uw(sw.get_usw());
		if (uw.cost == 0) continue; // these are not real weapons
		if (!uw.is_decoy  && uw.speed > 0.0 && uw.speed < dv)           continue; // can't hit target - weapon is too slow
		if (uw.is_decoy   && target_obj->disabled())                    continue; // decoy is useless (what about targets with no seeking weapons?)
		if (uw.is_fighter && sw.min_damage() > 0.3 && cur_damage < 0.8) continue; // fighter(s) are too damaged to release
		if (uw.is_fighter && sclasses[uw.ammo_type].for_boarding && !target_obj->can_board()) continue;
		if (!uw.is_beam   && !uw.is_fighter && !uw.is_decoy && uw.damage > 0.0 && uw.damage <= target_obj->get_min_damage()) continue; // can't damage target
		float const range(uw.range);
		
		if (!uw.is_decoy && !uw.is_fighter && range > 0.0 && (range + c_radius) > 0.0) {
			float tdist(target_dist);

			if (sw.weap_pts.size() == 1) {
				point wpos(sw.weap_pts[0]);
				xform_point_inv(wpos); // calculate more accurate distance based on actual weapon pos
				tdist = p2p_dist(wpos, target_obj->get_pos());
			}
			if (tdist > (uw.is_beam ? 1.0 : 1.05)*range) continue; // out of range
		}
		if (!weap_turret(sw.wclass) && bad_angle(d_angle, target_dist, sw.wclass)) { // target not aligned
			default_next_weap = w; // if the only useable weapon is unaligned, then make this the next cur_weapon so that we can align ourselves next frame
			continue;
		}
		if (uw.is_beam) {
			beam_weap_params const &bwp(uw.get_beam_params());
			if (bwp.temp_src && target_obj->get_max_t() > get_max_t()) continue; // can't damage target
			if (bwp.mind_control && (target_crew == 0 || target_crew > get_ncrew())) continue; // no crew/too many crew to mind control
			if (bwp.mind_control && alignment != ALIGN_NEUTRAL && alignment == target_obj->get_align()) continue; // don't mind control a friendly
			if (uw.damage == 0.0 && uw.force != 0.0 && !bwp.mind_control && target_obj->is_orbiting())  continue; // tractor beam on immovable ship
		}
		good_weapons.push_back(w); // acceptable weapon
	}
	if (good_weapons.empty()) { // no weapon for the given range, etc.
		curr_weapon = default_next_weap;
		return;
	}
	float const rsum(c_radius + target_obj->get_c_radius()); // can target_obj be NULL?
	float min_value(-2000.0); // set to 0 or positive to conserve ammo, etc.
	static vector<pair<float, unsigned> > wchoices;
	wchoices.resize(0);

	for (unsigned i = 0; i < good_weapons.size(); ++i) { // choose a weapon
		ship_weapon const &sw(weapons[good_weapons[i]]);
		us_weapon const &uw(sw.get_usw());
		float const range(uw.range + c_radius), overrange((uw.range > 0.0)*(range - target_dist));
		float value(uw.preference);
		value += 4.0*(AI_MULTI_FIRE && check_fire_delay(good_weapons[i])); // ready to fire bonus
		value += 2.5*(sw.init_ammo == 0); // infinite ammo bonus
		value += 0.8*min(1.5f, logf((float)sw.ammo + 1.0))/((float)sw.wcount + 1.0); // ammo bonus
		value -= 2.0*overrange/(0.5*range + rsum) + 8.0*overrange; // overrange penalty
		value -= 5.0*(sw.init_ammo - sw.ammo)/((float)sw.init_ammo + 1.0); // used ammo penalty
		value -= 8.0*(uw.range > 0.0 && target_dist > range); // out of range penalty
		value -= 6.0*d_angle*(1 - weap_turret(sw.wclass))*(uw.seeking ? 0.5 : 1.0); // bad aim penalty

		if (uw.is_beam) {
			beam_weap_params const &bwp(uw.get_beam_params());
			if (bwp.temp_src) {value += 5.0*TEMP_FACTOR*uw.damage*get_max_t();}
			if (bwp.temp_src) {value *= 0.5*(1.0 + get_max_t()/target_obj->get_max_t());}
			
			if (bwp.mind_control) { // prefers a mostly undamaged target
				value += 40.0*(1.0 - target_obj->get_damage())*(target_obj->get_cost() +
					10.0*(target_obj->offense() + target_obj->defense()))/max(1U, target_crew);
			}
		}
		value *= 1.0 + 0.5*uw.no_coll; // can hit multiple targets and can't be destroyed easily
		value *= sw.wcount;
		//cout << "[" << i << "]: " << uw.name << ", value = " << value << endl;
		if (value >= min_value) wchoices.push_back(make_pair(-value, good_weapons[i])); // want largest to smallest so negate w
	}
	if (wchoices.empty()) return; // probably won't get here
	sort(wchoices.begin(), wchoices.end());
	bool fired_secondary(0);

	for (unsigned i = 0; i < wchoices.size(); ++i) {
		curr_weapon = wchoices[i].second;
		unsigned const weapon_id(get_weapon_id());
		us_weapon const &weap(us_weapons[weapon_id]);
		if (!MULTI_SEC_FIRE && fired_secondary && weap.secondary) continue; // already fired secondary weapon
		if ( INTER_RELOAD_FIRE && !check_fire_delay(curr_weapon)) continue; // not reloaded yet
		if (weap.secondary) fired_secondary = 1;
		if (!INTER_RELOAD_FIRE && !check_fire_delay(curr_weapon)) continue; // not reloaded yet
		unsigned const wcount(get_weapon().wcount);
		assert(wcount > 0);
		s_turret = weap_turret(weapon_id); // update in case curr_weapon changed
		
		if (!weap.is_fighter && need_ammo()) { // conserve ammo (could divide by fticks)
			unsigned const mmult((weap.seeking && target_dist < weap.seek_dist) ? 2 : 4);
			unsigned const modval(unsigned(mmult*weap.fire_delay/wcount));
			if (modval > 1 && (rand() % modval) != 0) continue;
		}
		vector3d fire_dir((s_turret ? targ_dir : dir));

		if (PREDICT_TARGETS && s_turret && !weap.is_beam && weap.speed > 0.0) {
			fire_dir = predict_target_dir(pos, target_obj, weapon_id);
			if (fire_dir == zero_vector) continue;
		}
		free_obj *fobj;
		float line_radius(0.0), tdist; // tdist is unused
		if (!weap.is_beamlike() && !weap.is_fighter && !weap.no_ffire) line_radius = weap.radius;
		uobject const *line_of_fire(setup_int_query(fire_dir, target_dist, fobj, tdist, 0, line_radius));
		assert(line_of_fire == NULL || line_of_fire != target_obj);

		if (!weap.is_fighter && line_of_fire != NULL && (fobj == NULL || fobj->is_stationary() || (!is_enemy(fobj) &&
			(!(weap.no_ffire && weap.const_dam) || !is_related(fobj)))))
		{
			continue; // don't shoot if shot is blocked by a non-enemy or fixed object
		}
		//cout << weap.name << " " << target_obj->get_name() << endl; // testing
		if (target_dist < TOLERANCE || fire_dir.mag() < TOLERANCE) continue; // not sure how we can get here
		fire_weapon(fire_dir, target_dist);
	}
	curr_weapon = wchoices[0].second; // reset back to first
}


bool u_ship::test_self_intersect(point const &fpos, point const &tpos, vector3d const &fdir, float target_dist) {

	if (!specs().dynamic_cobjs && !weap_turret(get_weapon_id())) return 0; // fixed weapons can't self intersect unless dynamic cobjs

	if (line_sphere_intersect(fpos, tpos, pos, c_radius)) {
		point const f_start(fpos + fdir*((custom_wpt() ? 0.05 : 0.2)*radius)); // not sure about this 0.2
		point const f_end  (fpos + fdir*target_dist);
		if (line_int_obj(f_start, f_end)) return 1;
	}
	return 0;
}


bool u_ship::fire_weapon(vector3d const &fire_dir, float target_dist) {

	unsigned const wcount(get_weapon().wcount);
	if (reset_timer > 0 || wcount == 0)          return 0; // can't fire when dying or have no weapon
	unsigned const weapon_id(get_weapon_id());
	assert(weapon_id < us_weapons.size());
	us_weapon const &weap(us_weapons[weapon_id]);
	if (!weap.hyper_fire && !check_fire_speed()) return 0;
	ship_weapon &sw(weapons[curr_weapon]);
	sw.last_fframe = frame_counter;
	float const wrange(weap.range + c_radius), cscale(sqrt(get_crew_scale())); // firing error is slightly higher with fewer crew
	float const firing_error(weap.firing_error/max(0.001f, cscale*specs().stability));
	int const base_intersect_type(get_line_query_obj_types(wrange));
	int intersect_type(base_intersect_type);
	if (weap.hit_proj) intersect_type |= OBJ_TYPE_PROJ;
	if (weap.hit_all)  intersect_type |= OBJ_TYPE_FREE;
	unsigned num_shots(1), shots_per_round(1), remainder(0), used_ammo;
	free_obj *fobj;
	static vector<point> offsets;
	offsets.resize(0);
	point const target_pos(pos + fire_dir*target_dist); // recompute target pos since it may account for motion prediction
	unsigned const nwpts((unsigned)sw.weap_pts.size()), wspread(specs().weap_spread);

	if (weap.is_beam) {
		num_shots = wcount;
		used_ammo = iticks*num_shots;
	}
	else {
		if (weap.parallel_fire || specs().parallel_fire) {
			if (nwpts > 1) {
				num_shots = nwpts;
			}
			else if (wspread > 1) {
				num_shots = wspread;
			}
		}
		used_ammo = num_shots;
	}
	if (weap.need_ammo()) {
		num_shots = min(num_shots, get_weapon().ammo);
		used_ammo = min(used_ammo, get_weapon().ammo);
	}
	assert(num_shots > 0);

	// currently only works for beam weapons fired by the AI
	if (nwpts > 0 && num_shots >= nwpts) {
		shots_per_round = num_shots/nwpts;
		remainder       = num_shots - shots_per_round*nwpts;
		
		for (unsigned i = 0; i < nwpts; ++i) {
			point pos(sw.weap_pts[i]);
			rotate_point_inv(pos);
			offsets.push_back(pos);
		}
	}
	else if (wspread > 1 && num_shots > wspread && target_dist > 0.0) { // weapons fire from around the ship, not from a single point
		// generate <num_shots> points around a circle centered at (0,0,0) with radius=1 and facing in direction <target_dir>
		shots_per_round = num_shots/wspread;
		remainder       = num_shots - shots_per_round*wspread;
		vector3d const vxx(cross_product(upv, fire_dir).get_norm()), vxy(cross_product(vxx, fire_dir).get_norm());
		float const fire_r(0.9);
		
		for (unsigned i = 0; i < wspread; ++i) { // shouldn't generate self-collisions
			float const theta(TWO_PI*((float)i)/((float)wspread));
			offsets.push_back(vxx*(fire_r*sinf(theta))+ vxy*(fire_r*cosf(theta)));
		}
	}
	else {
		shots_per_round = num_shots;

		if (nwpts > 0) {
			sw.cur_wpt = (sw.cur_wpt+1)%nwpts;
			point pos(sw.weap_pts[sw.cur_wpt]);
			rotate_point_inv(pos);
			offsets.push_back(pos);
		}
		else {
			offsets.push_back(all_zeros); // change this for special weapons, and use this for targeting/leading shots, etc.
		}
	}
	assert(shots_per_round > 0);
	assert(!offsets.empty());
	assert(offsets.size() <= num_shots);
	bool const s_turret(weap_turret(weapon_id)), multi_target(do_multi_target());
	bool fired(0);

	for (unsigned o = 0; o < offsets.size(); ++o) {
		unsigned nshots(shots_per_round);
		if (o < remainder) ++nshots;
		assert(num_shots >= nshots);
		num_shots -= nshots;
		assert(nshots > 0);
		vector3d fdir(fire_dir);
		point fpos(pos);
		float cur_target_dist(target_dist);
		free_obj const *tobj(target_obj);
		
		if (cur_target_dist > TOLERANCE && offsets[o] != all_zeros) {
			point tpos(target_pos);
			fpos += offsets[o]*radius; // fire from this point

			if (multi_target) {
				free_obj const *new_tobj(find_closest_target(fpos, get_min_att_dist(), weap.range, weap.shield_d_only));
				
				if (new_tobj != NULL) {
					tobj = new_tobj;
					point const new_tpos(tobj->get_pos());
					vector3d const temp_fd((new_tpos - fpos).get_norm());
					float const new_target_dist(p2p_dist(new_tpos, fpos));

					if (temp_fd != all_zeros && (weap.is_fighter || !test_self_intersect(fpos, new_tpos, temp_fd, new_target_dist))) {
						tpos            = new_tpos; // new target line does not intersect ourselves
						cur_target_dist = new_target_dist;
					}
				}
			}
			fdir = (tpos - fpos).get_norm(); // recompute direction to target from new firing location

			if (PREDICT_TARGETS && multi_target && !weap.is_beam && weap.speed > 0.0) {
				vector3d const new_fdir(predict_target_dir(fpos, tobj, weapon_id));
				if (new_fdir != zero_vector) fdir = new_fdir;
			}
			if (!weap.is_fighter && test_self_intersect(fpos, tpos, fdir, cur_target_dist)) {
				if (used_ammo >= nshots) used_ammo -= nshots; else used_ammo = 0; // not quite correct?
				continue; // collision with self
			}
			if (fdir == zero_vector) continue; // really shouldn't get here...
			target_dir = fdir;
		}
		if (tobj && !player_controlled()) {
			if ((weap.shield_d_only && !tobj->has_shields()) ||
				(!weap.is_beam && !weap.is_fighter && !weap.is_decoy && weap.damage < tobj->get_damage_abs()))
			{
				if (used_ammo >= nshots) used_ammo -= nshots; else used_ammo = 0; // not quite correct?
				continue; // can't damage the target with this weapon
			}
		}
		if (firing_error > 0.0) { // add in firing error
			vadd_rand(fdir, (s_turret ? (2.0*firing_error + 0.05) : firing_error)*(radius + 0.01));
			fdir.normalize();
		}
		switch (weapon_id) {
			case UWEAP_NONE: // shouldn't get here
				return 0; // nothing

			case UWEAP_RENAME: // player only
			case UWEAP_TARGET:
			case UWEAP_QUERY:
				{
					assert(is_player_ship());
					bool const rename(weapon_id == UWEAP_RENAME), query(weapon_id == UWEAP_TARGET), target(weapon_id == UWEAP_TARGET);
					line_int_data li_data(fpos, fdir, wrange, this, NULL, 0, (query ? 2 : 0));
					uobject *sobj(line_intersect_objects(li_data, fobj, (target ? OBJ_TYPE_SHIP : intersect_type)));
					
					if (sobj != NULL && rename) {
						rename_obj(sobj, alignment);
						break;
					}
					if (fobj != NULL && !fobj->invalid()) { // can query projectiles for now
						string msg;
						if (target) {
							msg = string("Target Acquired: ") + fobj->get_name();
						}
						else {
							msg = fobj->get_name() + "\n" + fobj->get_info() + "  Dist: " + make_string(li_data.dist);
						}
						print_text_onscreen(msg, PURPLE, 0.5, TICKS_PER_SECOND, 1);

						if (target && fobj->is_target() && target_valid(fobj)) {
							target_obj = fobj;
							target_set = 1;
						}
					}
					else if (sobj != NULL) {
						string const msg(sobj->get_name() + " at distance " + make_string(li_data.dist) + "\n" + sobj->get_info());
						print_text_onscreen(msg, PURPLE, 0.5, TICKS_PER_SECOND, 1);
					}
				}
				break;

			case UWEAP_DESTROY:
				fire_planet_killer(this, fpos, fdir, wrange, intersect_type); // see universe.cpp
				break;

			case UWEAP_CHAFF:
				{
					vector3d v(fdir*weap.speed);
					point p(fpos + fdir*(1.5*c_radius));

					for (unsigned j = 0; j < weap.nshots; ++j) {
						v += signed_rand_vector()*(0.04*weap.speed*weap.firing_error);
						p += signed_rand_vector()*(0.5*weap.radius*weap.firing_error); // random walk
						//gen_particle(PTYPE_TRIANGLE, LT_GRAY, LT_GRAY, int(weap.lifetime), p, v, weap.radius, 0.0, alignment, 1);
						create_projectile(us_weapons[get_weapon_id()].ammo_type, this, alignment, p, v, signed_rand_vector(), upv);
					}
				}
				break;

			default:
				if (weap.is_beamlike()) { // beams, LRCPA, PT_DEF
					fire_beam(fpos, fdir, weapon_id, nshots, intersect_type);
				}
				else if (weap.is_fighter) {
					assert(weap.is_fighter);
					if (nshots > 1) used_ammo -= nshots-1; // can only launch one at a time
					unsigned const stype(weap.ammo_type);
					assert(stype < sclasses.size());
					float const f_radius(sclasses[stype].calc_cradius());
					float const launch_dist(custom_wpt() ? 0.0 : (1.1*c_radius + 2.0*f_radius)*rand_uniform(0.8, 1.2));
					line_int_data li_data(fpos, fdir, (launch_dist + f_radius), this, NULL, 0, 0);

					if (line_intersect_objects(li_data, fobj, base_intersect_type) != NULL) {
						if (used_ammo >= nshots) used_ammo -= nshots; else used_ammo = 0;
						continue; // can't launch a fighter - something is in the way
					}
					spawn_fighter((fpos + fdir*launch_dist), fdir);
				}
				else if (weap.ammo_type != UWEAP_NONE) {
					if (nshots > 1) used_ammo -= nshots-1; // can only fire one at a time
					fire_projectile(fpos, fdir);
				}
				else {
					cout << "Invalid weapon id: " << weapon_id << endl;
					assert(0);
				}
		} // switch (weapon_id)
		fired = 1;
	} // for o
	if (fired) {partial_uncloak(0.25);} // must de-cloak partially when firing
	if (used_ammo > 0) {dec_ammo(used_ammo);}

	// generate weapon fire sounds
	if (fired && (is_player_ship() || dist_less_than(get_camera_pos(), pos, 0.1))) {
		float const gain(is_player_ship() ? 1.0 : 0.1);
		switch (weapon_id) {
		//case UWEAP_PBEAM: case UWEAP_EBEAM: case UWEAP_REPULSER: case UWEAP_TRACTORB: case UWEAP_FUSCUT: case UWEAP_LITNING: case UWEAP_PARALYZE: case UWEAP_MIND_C: case UWEAP_ESTEAL
		case UWEAP_DESTROY: break;
		case UWEAP_LRCPA:   break;
		case UWEAP_ENERGY:  break;
		case UWEAP_ATOMIC:  gen_sound(SOUND_GUNSHOT,  pos, gain, 1.0);  break;
		case UWEAP_SHIELDD: gen_sound(SOUND_GUNSHOT,  pos, gain, 1.5);  break;
		case UWEAP_ROCKET:  gen_sound(SOUND_ROCKET,   pos, gain, 1.0);  break;
		case UWEAP_NUKEDEV: gen_sound(SOUND_ROCKET,   pos, gain, 0.75); break;
		case UWEAP_TORPEDO: gen_sound(SOUND_ITEM,     pos, gain, 1.0);  break;
		case UWEAP_THUNDER: gen_sound(SOUND_ITEM,     pos, gain, 0.75); break;
		case UWEAP_EMP:     gen_sound(SOUND_POWERUP,  pos, gain, 1.0);  break;
		case UWEAP_SEIGEC:  gen_sound(SOUND_SHOTGUN,  pos, gain, 1.0);  break;
		case UWEAP_RFIRE:   gen_sound(SOUND_FIREBALL, pos, gain, 1.0);  break;
		case UWEAP_INFERNO: gen_sound(SOUND_FIREBALL, pos, gain, 1.5);  break;

		case UWEAP_FIGHTER: case UWEAP_B_BAY: case UWEAP_CRU_BAY: case UWEAP_SOD_BAY: case UWEAP_BOARDING: case UWEAP_NM_BAY: // fallthrough
		case UWEAP_WRAI_BAY: case UWEAP_HUNTER: case UWEAP_DEATHORB: case UWEAP_SAUC_BAY: gen_sound(SOUND_HISS, pos, gain, 1.0); break;
		}
	}
	return 1;
}


void u_ship::fire_beam(point const &fpos, vector3d const &fdir, unsigned weapon_id, unsigned num_shots, int intersect_type) {

	assert(num_shots > 0);
	free_obj *fobj;
	us_weapon const &weap(us_weapons[weapon_id]);
	assert(weap.radius > 0.0 && weap.c_radius > 0.0);
	unsigned const check_parent(weap.no_ffire ? 2 : 1);
	float const wscale(sqrt((float)num_shots)), beamwidth(weap.radius*wscale), blastwidth(weap.c_radius*wscale);
	beam_weap_params const &bwp(weap.get_beam_params());
	float range(weap.range + c_radius);
	assert(beamwidth > 0.0 && blastwidth > 0.0);
	point p1(fpos), p2(fpos);
	if (!custom_wpt()) p1 += fdir*(0.5*radius); // move the starting point away from the ship
	p1 += velocity*fticks; // display hack: adjust for velocity of the next frame
	line_int_data li_data(fpos, fdir, range, this, NULL, 0, check_parent);
	uobject const *const hit_obj(line_intersect_objects(li_data, fobj, intersect_type));
	
	if (hit_obj != NULL) {
		unsigned const etype(weap.exp_type);
		float const hradius(hit_obj->get_bounding_radius());
		point const hpos(hit_obj->get_pos());
		if (!line_sphere_int(fdir, fpos, hpos, hradius, p2, 0)) return; // no intersection (more robust test) (Note: inexact for asteroids)
		float dscale((weap.is_beam && !bwp.paralyze && !bwp.temp_src) ? min(4.0f, fticks) : 1.0);

		if (fobj != NULL) { // update p2 if we have more accurate data
			float cobj_dscale(1.0);
			point end_pos(fpos);
			float const fdir_mag(fdir.mag());
			if (fdir_mag > TOLERANCE) end_pos += fdir*(range/fdir_mag);
			fobj->line_int_obj(fpos, end_pos, &p2, &cobj_dscale);
			dscale *= cobj_dscale;
			p2     += fobj->get_velocity()*fticks; // display hack: adjust for velocity of the next frame
		}
		float const br(max(wscale*weap.bradius, blastwidth));

		if (etype != ETYPE_NONE && !is_distant(p2, 0.5*br)) {
			int const time(max(int(TICKS_PER_SECOND/10), int(0.2*TICKS_PER_SECOND*(0.75 + 0.25*et_params[etype].duration))));
			add_blastr(p2, zero_vector, br, 0.0, time, alignment, bwp.brc[0], bwp.brc[1], etype, this); // fdir?
		}
		vector3d vcoll(fdir*weap.force); // positive for repel, negative for attract (tractor beam)
		
		if (fobj != NULL && (!weap.hit_all || fobj->get_parent() != this) && dscale > 0.0) { // hit a free_obj and did damage
			float const damage(dscale*weap.damage*num_shots);
			
			if (weap.force != 0.0 && weap.mass != 0.0) {
				fobj->coll_physics(p2, vcoll, weap.mass, weap.radius, NULL, 1.0); // apply force to target object
			}
			if (bwp.temp_src) { // no dscale
				float const sval(get_state_val());

				if (sval > 0.25) {
					float const temperature(damage*TEMP_FACTOR*sval*specs().max_t);
					temp_sources.push_back(temp_source(p2, weap.bradius, temperature, this));
				}
			}
			else if (bwp.paralyze) {
				fobj->try_paralyze(this, damage, p2);
			}
			else if (bwp.mind_control) {
				if (player_controlled() || fobj == target_obj) fobj->try_mind_control(this, num_shots);
			}
			else {
				float const damage_done(fobj->damage(damage, DAMAGE_HEAT, p2, this, weapon_id)); // beam is heat damage

				if (bwp.energy_drain) {
					energy += ENERGY_XFER_EFF*damage_done*(4.0f/max(4.0f, fticks)); // scale to reasonable values for slow framerate
				}
			}
		}
		if (weap.force > 0.0 && weap.f_inv > 0.0) { // inverse force applied to the shooter (like recoil)
			coll_physics(p2, vcoll*-weap.f_inv, weap.f_inv*weap.mass, weap.radius, NULL, 1.0);
		}
	}
	else {
		assert(range > 0.0);
		p2 += fdir*range;
	}
	float const length(p2p_dist(p1, p2));
	if (length < 0.1*radius) return; // rarely happens - only if the ship and target are intersecting (or the ship hits itself?)

	// draw beam (even for player)
	assert(beamwidth > 0.0 && radius > 0.0 && fticks > 0.0 && bwp.bw_escale > 0.0);
	float const sscale(is_player_ship() ? 0.2 : 1.0);

	if (!bwp.multi_segment) {
		b_wrays.push_back(usw_ray(sscale*beamwidth, bwp.bw_escale*beamwidth, p1, p2, bwp.beamc[0], bwp.beamc[1]));
	}
	else {
		unsigned const segments((rand()&7)+4);
		float const step(length/segments);
		vector3d deltas[12];
		point prev(p1);
		deltas[0] = deltas[segments] = zero_vector; // first and last segments are zero

		for (unsigned i = 1; i < segments; ++i) {
			do {
				deltas[i] = signed_rand_vector(0.1*step);
			} while (deltas[i] == deltas[i-1]);
		}
		for (unsigned i = 0; i < segments; ++i) {
			float bw[2];
			point pt[2];
			colorRGBA c[2];

			for (unsigned d = 0; d < 2; ++d) {
				float const val(float(i+d)/float(segments));
				bw[d] = beamwidth*(val*bwp.bw_escale + (1.0 - val)*sscale);
				pt[d] = (p2*val + p1*(1.0 - val)) + deltas[i+d];
				blend_color(c[d], bwp.beamc[1], bwp.beamc[0], val, 1);
			}
			if (pt[0] != pt[1]) {
				if (i > 0) {b_wrays.back().next = pt[1];}
				b_wrays.push_back(usw_ray(bw[0], bw[1], pt[0], pt[1], c[0], c[1]));
				if (i > 0) {b_wrays.back().prev = prev;}
				prev = pt[0];
			}
		}
	}
}


void u_ship::fire_projectile(point const &fpos, vector3d const &fire_dir) {

	us_weapon const &weap(us_weapons[get_weapon_id()]);
	unsigned const weapon_id(weap.ammo_type);
	assert(get_weapon_id() == weapon_id); // sanity check
	point const ppos(fpos + fire_dir*(1.0*((custom_wpt() ? 0.0 : radius) + weap.radius)));
	vector3d const pvel((weap.no_ship_vel ? zero_vector : velocity) + fire_dir*weap.speed);
	create_projectile(weapon_id, this, alignment, ppos, pvel, fire_dir, upv); // return value unused

	if (weap.f_inv > 0.0) { // inverse force applied to the shooter (recoil)
		coll_physics(ppos, pvel*-weap.f_inv, weap.f_inv*weap.mass, weap.radius, NULL, 1.0);
	}
}


void u_ship::spawn_fighter(point const &fpos, vector3d const &fire_dir) {

	assert(fire_dir != zero_vector);
	unsigned const wid(get_weapon_id()), type(us_weapons[wid].ammo_type);
	assert(type < sclasses.size() && type < us_weapons.size());
	bool const boarding(sclasses[type].for_boarding); //specs().for_boarding
	unsigned const targeting((boarding && !player_controlled()) ? TARGET_CLOSEST : TARGET_PARENT);
	u_ship *ship(create_ship(type, fpos, alignment, AI_ATT_ENEMY, targeting, 0, is_rand_spawn()));
	ship->target_obj = NULL;
	ship->velocity   = zero_vector; // probably unnecessary
	ship->dvel       = zero_vector;
	ship->set_dir(us_weapons[wid].auto_orient ? fire_dir : dir);
	ship->set_upv(upv);
	ship->no_ai_wait();
	ship->check_distant();
	weapons[curr_weapon].release_ship(ship, this);
	add_fighter(ship, this, 1);
}


bool u_ship::try_fire_weapon() { // player fire

	if (invalid_or_disabled()) return 0;
	unsigned const temp_cw(curr_weapon), nweap((unsigned)weapons.size());
	bool fired_secondary(0), fired(0);
	target_dir = dir;
	if (target_obj) target_dir = (target_obj->get_pos() - pos).get_norm();

	for (unsigned i = 0; i < nweap; ++i) {
		unsigned const ix((i + temp_cw) % nweap); // try curr_weapon first
		if (!fire_primary && ix != temp_cw) continue;
		curr_weapon = ix;
		unsigned const wid(get_weapon_id());
		us_weapon const &weap(us_weapons[wid]);
		if (ix != temp_cw && (weap.is_fighter || wid == UWEAP_DESTROY || wid == UWEAP_TARGET ||
			wid == UWEAP_QUERY || wid == UWEAP_RENAME)) continue;
		if (!INTER_RELOAD_FIRE && !check_fire_delay(curr_weapon)) continue; // not reloaded yet
		if (!MULTI_SEC_FIRE && fired_secondary && weap.secondary && curr_weapon != temp_cw) continue; // already fired secondary weapon
		if (get_weapon().wcount == 0 || out_of_ammo(1)) continue;
		vector3d fire_dir;
		free_obj const *ctarg = NULL;
		
		if (weap_turret(wid)) { // weapon is turreted
			ctarg = target_obj;
			float search_dist(specs().sensor_dist);
			if (weap.range > 0.0) search_dist = min(search_dist, (c_radius + weap.range)); // ???
			if (ctarg == NULL)    ctarg       = get_closest_ship(pos, 0.0, search_dist, 1, 0, 0); // no target selected
		}
		if (ctarg != NULL) { // have a target
			fire_dir = (ctarg->get_pos() - pos).get_norm();
		}
		else { // fire straight ahead
			fire_dir = get_dir();
		}
		// can't always know what target to converge on, so can only fire weapons in a straight line from the center of the ship
		float const target_dist((target_set && target_obj != NULL) ? p2p_dist_sq(pos, target_obj->get_pos()) : 0.0);

		if (fire_weapon(fire_dir, target_dist)) {
			if (weap.secondary) fired_secondary = 1;
			fired = 1;
		}
	}
	curr_weapon = temp_cw;
	return fired;
}


// UWEAP_NONE => calculate ship dir, not weapon dir
vector3d u_ship::predict_target_dir(point const &fpos, free_obj const *targ, unsigned wclass) const {

	// Note: can more accurately test weapon range this way
	assert(targ != NULL);
	point const pt(targ->get_pos());

	// get weapon projectile speed
	if (wclass != UWEAP_NONE && us_weapons[wclass].speed == 0.0) { // not a projectile weapon, just face the target
		return vector3d(pt, fpos).get_norm();
	}
	double const vweap((wclass == UWEAP_NONE) ? /*velocity.mag()*/get_max_speed() : us_weapons[wclass].speed);
	return lead_target(fpos, pt, velocity, targ->get_velocity(), vweap);
}


inline float u_ship::get_min_att_dist() const {
	
	float const mad(radius*specs().min_att_dist);
	if (specs().for_boarding) return mad;
	us_weapon const &usw(us_weapons[get_weapon_id()]);
	return max(mad, (1.5f*c_radius + usw.bradius + usw.radius));
}


bool u_ship::check_fire_speed() const {

	if (specs().max_speed == 0.0) return 1; // can't move, but can fire

	switch(specs().fire_speed) {
	case FSPEED_HYPER:
		return 1;
	case FSPEED_SLOW:
		if (get_max_sf() > SLOW_SPEED_FACTOR) return 0;
	case FSPEED_FAST:
		return (velocity.mag() <= specs().max_speed);
	default:
		assert(0);
	}
	return 0;
}


void u_ship::switch_weapon(bool prev) {

	if (prev) curr_weapon += ((unsigned)weapons.size()-1); else ++curr_weapon;
	curr_weapon = curr_weapon % weapons.size();
	show_weapon_name();
}


void u_ship::show_weapon_name() const {

	unsigned const weapon_id(get_weapon_id());
	assert(us_weapons.size() == NUM_UWEAP && weapon_id < us_weapons.size());
	print_text_onscreen(us_weapons[weapon_id].name.c_str(), WHITE, 1.0, TICKS_PER_SECOND, 10); // higher priority to override enemy ship detected
}


inline void u_ship::dec_ammo(unsigned num) {

	if (need_ammo()) {
		assert(num > 0 && get_weapon().ammo >= num);
		weapons[curr_weapon].ammo -= num;
	}
}


void u_ship::reset_target() {
	
	if (invalid_priv()) return;

	for (set<u_ship *>::const_iterator i = fighters.begin(); i != fighters.end(); ++i) {
		if ((*i)->get_parent() == this && (*i)->get_target() == target_obj) (*i)->reset_target();
	}
	target_set = 0;
	free_obj::reset_target();
}


bool u_ship::dock_fighter(u_ship *ship) {

	assert(ship != NULL);
	if (parent == ship) return ship->dock_fighter(this); // was called backwards
	if (ship->get_parent() != this || ship->get_target() != this) return 0;
	if (invalid_or_disabled() || ship->invalid_or_disabled())     return 0;
	if (fighters.find(ship) == fighters.end()) add_fighter(ship, this, 0); // rarely happens
	
	if (!ship->can_return_to_parent()) {
		ship->target_obj = NULL; // can't return to parent any more
		return 0;
	}
	unsigned const fsclass(ship->sclass);
	unsigned ix(0);

	if (find_weapon_atype(fsclass, 1, ix)) {
		ship_weapon &w(weapons[ix]);
		assert(w.space_for_fighter()); // no room in fighter bay - should have checked this already
		if (!us_weapons[w.wclass].hyper_fire && !check_fire_speed()) return 0; // can't dock at this speed
		++w.ammo;
		//assert(w.ammo <= w.init_ammo); // not valid if ship docks a dead fighter's fighters
		w.init_ammo = max(w.ammo, w.init_ammo);
		w.dock_ship(ship, this);
		return 1;
	}
	return 0; // else can't dock this type of fighter, no bay that holds its type
}


bool u_ship::orbital_dock(u_ship *ship) { // ship docks with orbital station

	assert(ship != this);
	if (invalid_priv() || ship->invalid_priv()) return 0;
	bool const o_dock(specs().orbiting_dock), orbiting(is_orbiting()), ship_orbiting(ship->is_orbiting());
	if (o_dock) assert(orbiting);
	if (ship_orbiting && orbiting) return 0; // both are orbiting - what to do?
	if (ship_orbiting) return ship->orbital_dock(this); // what if both are orbiting?
	if (!orbiting || disabled_priv() || ship->alignment != alignment || !o_dock) return 0;
	orbital_ship_regen(ship);
	return 1;
}


float u_ship::get_crew_strength() const {

	float val(ncrew + specs().mass); // larger, heavier ships have more automated defenses
	if (specs().for_boarding) val *= 4.0; // better prepared/equipped for battle
	val *= 0.5 + min(0.5f, get_armor()/max(1.0f, get_max_armor())); // ship is damaged and more difficult to defend

	for (unsigned i = 0; i < weapons.size(); ++i) {
		if (us_weapons[weapons[i].wclass].is_fighter && sclasses[us_weapons[weapons[i].wclass].ammo_type].for_boarding) {
			val += 2.0*weapons[i].ammo; // ship contains a boarding shuttle itself
		}
	}
	return val;
}


float u_ship::get_child_stray_dist() const {
	
	if (child_stray_dist > 0.0)              return child_stray_dist;
	if (can_move() && !specs().for_boarding) return specs().sensor_dist;
	return 0.0;
}


float u_ship::min_time_to_target(point const &targ_pos) const { // in seconds

	float const max_speed(specs().max_speed);
	if (max_speed == 0.0) return 0.0; // what to do - can't move?
	float const dist(p2p_dist(targ_pos, pos)), fast_tdist(specs().has_fast_speed ? get_fast_target_dist() : dist);
	float const slow_dist(min(dist, fast_tdist)), fast_dist(max(0.0f, (dist - fast_tdist)));
	return ((slow_dist/SLOW_SPEED_FACTOR) + (fast_dist/FAST_SPEED_FACTOR))/(TICKS_PER_SECOND*max_speed);
}


// this boards ship - boarding on collision ensures the ship is close enough
bool u_ship::board_ship(u_ship *ship) { // returns whether or not boarding was attempted

	assert(ship != NULL);
	if (invalid_or_disabled() || !specs().for_boarding || ship != target_obj || ship->invalid() || !ship->specs().can_board) return 0;
	if (ship->captured)   return 0; // not sure what to do with two captures in the same frame
	if (is_related(ship)) return 0; // boarding a friendly, what else can we do here? transfer crew/ammo/shields?
	ship->register_attacker(this); // counts as attacking the ship
	if (ship->shields_up() || ship->get_damage() < 0.5) return 1; // not damaged enough, maybe don't even target it in this case?
	//cout << this << " is boarding " << ship << ", strength: " << get_crew_strength() << " vs. " << ship->get_crew_strength() << endl; // testing
	assert(ncrew > 0);
	
	// capture crew and ammo?
	if (get_crew_strength()*((rand()%100) + 10) > ship->get_crew_strength()*((rand()%100) + 10)) { // at most a factor of 11
		//cout << get_name() << " " << this << " captured " << ship->get_name() << " " << ship << endl; // testing
		capture_ship(ship, 0); // attackers defeat defenders
	}
	ncrew = max(1U, ncrew/2); // lose half the crew, both sides take casualties?
	return 1;
}


bool u_ship::capture_ship(u_ship *ship, bool add_as_fighter) { // this captures ship - not entirely tested and has problems

	// too bad we can't remove ourself from our parent's fighter list or add our parent as the captured ship's parent
	assert(ship != NULL && ship != this);
	free_obj const *const old_parent(ship->parent);
	bool const player(ship->is_player_ship());
	if (player) destroy_player_ship(1);
	register_damage(sclass, ship->sclass, WCLASS_CAPTURE, (ship->get_shields() + ship->get_armor()), alignment, ship->get_align(), 1); // counts as a kill
	if (is_player_ship()) acknowledge_kill(); // counts as a kill - player only for now
	ship->register_destruction(this); // counts as destroying the ship
	ship->alignment  = alignment;
	ship->parent     = this; // or could set to NULL?
	ship->target_obj = ((target_obj == ship) ? NULL : target_obj);
	ship->captured   = 1;

	if (player) { // create a ship of the same type as the player's ship, but can't have the special player weapons
		u_ship *new_ship(create_ship(ship->sclass, ship->pos, ship->alignment, ship->ai_type, ship->target_mode, 0));
		new_ship->target_obj = ship->target_obj;
		new_ship->set_dir(ship->dir);
		new_ship->set_upv(ship->upv);
		new_ship->create_from(*ship);
		new_ship->weapons = sclasses[new_ship->sclass].weapons; // can't have the unfair player weapons
		if (add_as_fighter) add_fighter(new_ship, this, 0);
	}
	else if (add_as_fighter) {
		add_fighter(ship, this, 0);
	}
	target_obj = NULL; // no longer targeting this ship
	
	for (set<u_ship *>::iterator i = ship->fighters.begin(); i != ship->fighters.end(); ) {
		set<u_ship *>::iterator cur(i++);
		u_ship *fighter(*cur);
		assert(fighter != NULL && fighter != ship && fighter != this && fighter != parent);

		if (!fighter->invalid() && fighter->get_parent() == ship) {
			// reset alignment, target, etc? - does fighter ignore, accept capturer as parent's parent, attack capturer, attack parent???
			if (rand()&1) { // keep <ship> as a parent
				fighter->target_obj = NULL; // reset target
				fighter->alignment  = ship->alignment;
			}
			else { // no longer a child of <ship>
				assert(old_parent != fighter); // old_parent can be NULL
				fighter->parent = old_parent; // what if <ship> is a itself fighter?
				ship->fighters.erase(cur);
			}
		}
	}
	// what else? transfer crew? restore shields/armor?
	return 1; // always successful, for now
}


bool u_ship::try_paralyze(u_ship const *source, float intensity, point const &ppos) {

	assert(source);
	register_attacker(source);
	if (intensity < 1.0E-6) return 0; // hit an undamageable part of the ship?
	float const resistance((has_shields() ? 2.0 : 1.0)*get_mass()), tval(10.0*intensity/resistance);
	unsigned ptime(0);
				
	if (tval >= 0.5) {
		ptime = unsigned(max(1.0f, fticks)*tval + 0.5);
	}
	else if ((rand() % unsigned(1.0/tval + 0.5)) == 0) {
		ptime = 1;
	}
	if (ptime == 0) return 0;
	disable(ptime);
	add_blastr(ppos, dir, 3.0*radius, 0.0, min(max(iticks, 1), int(ptime)), source->get_align(),
		YELLOW, ORANGE, ETYPE_NONE, source); // colored light
	return 1;
}


bool u_ship::try_mind_control(u_ship *source, unsigned num_tries) { // source mind controls this

	assert(source);
	if (num_tries == 0 || iticks == 0 || get_ncrew() == 0)     return 0;
	if (source->is_player_ship() && alignment == ALIGN_PLAYER) return 0; // already player team, don't attack
	register_attacker(source);
	unsigned const rval((shields_up() ? 2 : 1)*((rand()%40) + get_ncrew())/(3*iticks*num_tries));
	
	if (rval <= 1 || (rand()%rval) == 0) {
		colorRGBA const color_b(alignment_colors[get_align()]);
		colorRGBA const color_a(alignment_colors[source->get_align()]);
		add_blastr(pos, dir, 5.0*radius, 0.0, TICKS_PER_SECOND/4, source->get_align(),
			color_a, color_a, ETYPE_NONE, source); // colored light

		for (unsigned i = 0; i < 50; ++i) {
			float const psize(radius*rand_uniform(0.1, 0.25));
			vector3d const vel(velocity*0.1 + gen_rand_vector_uniform(1.0E-4*rand_uniform(1.0, 1.5)));
			gen_particle(PTYPE_GLOW, color_b, color_a, TICKS_PER_SECOND, (pos + vel*radius), vel, psize, 0.0, alignment, 0);
		}
		source->capture_ship(this, 0);
		return 1;
	}
	return 0;
}


bool u_ship::hostile_explode() const { // called on a ship or a parent
	
	if (specs().suicides) return 1;

	for (unsigned i = 0; i < weapons.size(); ++i) {
		if (weapons[i].init_ammo == 0) continue;
		us_weapon const &uw(weapons[i].get_usw());
		// too hard to determine which fighter just died, so just assume it was one of this type
		if (uw.is_fighter && sclasses[uw.ammo_type].suicides) return 1;
	}
	return 0;
}


inline bool u_ship::is_enemy(free_obj const *obj) const {
	
	assert(obj != NULL);
	if (!TEAM_ALIGNED(alignment)) return 0;
	unsigned const oal(obj->get_align());
	if (!TEAM_ALIGNED(oal))       return 0;
	return (alignment != oal);
}


void u_ship::get_fighter_target(u_ship const *ship) { // recursive

	if (target_obj != NULL || ship == NULL || invalid_priv()) return; // base case 1
	
	if (ship != this && ship->target_obj != NULL) { // base case 2
		if (target_valid(ship->target_obj)) target_obj = ship->target_obj;
		return;
	}
	for (set<u_ship *>::const_iterator i = ship->fighters.begin(); i != ship->fighters.end() && target_obj == NULL; ++i) {
		assert(*i != NULL); // see if a fighter of yours has a target
		if (!(*i)->invalid() && (*i)->get_parent() == this) get_fighter_target(*i);
	}
}


void u_ship::set_temp(float temp, point const &tcenter, free_obj const *source) {
	
	tcent = tcenter;
	free_obj::set_temp(temp, tcenter, source);
}


void u_ship::apply_physics() {

	bool const exploding(is_exploding_priv());
	us_class const &sc(specs());
	if (!exploding && armor <= 0.0 && !invalid_priv() && !is_resetting()) destroy_ship(0.0);
	
	if (exploding && sc.uses_mesh2d) {
		if (sc.mesh_remove) surface_mesh.unset_rand_rmap(unsigned(2.0*fticks));
		if (sc.mesh_deform) surface_mesh.mult_by(1.0 + 0.01*fticks);
		if (sc.mesh_expand) surface_mesh.set_rand_expand(0.1*rand_float(), unsigned(12.0*fticks));
		if (sc.mu_expand  ) surface_mesh.expand += 0.01*fticks;
		
		if (sc.mesh_trans) {
			point const tp(velocity * -10.0*fticks*rand_float()/radius);
			surface_mesh.set_rand_translate(tp, unsigned(10.0*fticks));
		}
	}
	if (!is_ok()) {
		last_hit = 0; // maybe remove this?
	}
	else {
		if (exploding) {
			if (explode_now()) {
				exp_time = 0;
				do_explode();
				return;
			}
			else if (sc.cont_frag && (rand() & 7) == 0) {
				fragment(dir, 0.05, 1);
			}
			set_ship_max_speed(2.0); // allow it to move faster than it normally could (due to explosions), but still cap it
		}
		else if (!invalid_priv()) {
			if (!fighters.empty()) get_fighter_target(this); // see if a fighter of yours has a target
			if (last_hit    > (unsigned)iticks) last_hit    -= iticks; else last_hit    = 0;
			if (retarg_time > (unsigned)iticks) retarg_time -= iticks; else retarg_time = 0;
			if (target_obj == NULL)             last_targ_t += iticks; else last_targ_t = 0;
			set_ship_max_speed(1.0);
			regen(1.0, NULL);
			check_size_scale();
		}
	}
	if (disabled() && armor >= DISABLE_ARMOR*get_max_armor()) disable_t -= iticks; // un-disables faster
	
	if (eflags != 0 && iticks > 0 && !disabled()) { // repair any engines that are down
		float const cscale(get_crew_scale()); // engines are repaired faster with more crew

		if ((rand() % max(1U, (ENG_REPAIR_TIME/max(1, int(cscale*iticks))))) == 0) {
			for (unsigned i = 0; i < sc.nengines; ++i) {
				if (eflags & (1 << i)) {eflags &= ~(1 << i); break;}
			}
		}
	}
	cached_rsv = 0.0;
	free_obj::apply_physics();
	bool const is_powered(powered_priv());
	if (is_powered) elapsed_on_t += iticks;

	// rand() mod should be f(fticks), but leave as is for improved performance when framerate is low
	if (ENABLE_PARTS && animate2 && is_powered && sc.nengines > 0 && cloaked < 0.5 && (rand()%6) == 0 &&
		fticks*velocity.mag() > 0.1*radius && univ_sphere_vis(pos, 2.0*c_radius))
	{ 
		float const pscale(0.1*min(5.0f*MAX_PARTICLE_SIZE, radius));

		if (!is_distant(pos, 2.0*pscale)) {
			point const pos2(pos + gen_rand_vector_uniform(0.5*radius) - dir*(1.5*c_radius));
			float const psize(pscale*rand_uniform(0.5, 1.0));
			vector3d const vel(velocity*-0.2);
			gen_particle(PTYPE_GLOW, LT_GRAY, colorRGBA(0.75, 0.75, 0.75, 0.0), ((5*TICKS_PER_SECOND)/2), pos2, vel, psize, 0.0, alignment, 0);
		}
	}
	if (exploding && specs().exp_subtype != ETYPE_NONE) { // many small explosions
		unsigned const niters(unsigned(rand()%10)/6);

		for (unsigned i = 0; i < niters; ++i) {
			point const bpos(pos + gen_rand_vector_uniform(radius));
			float const brad(0.5*(radius + c_radius)*(((specs().death_delay > 0) ? 0.4*get_t_exp() : 0.2) + rand_uniform(0.2, 0.5)));
			int const btime(int(TICKS_PER_SECOND*rand_uniform(0.5, 0.9)));
			add_blastr(bpos, (bpos - pos).get_norm(), brad, 0.0, btime, alignment, YELLOW, RED, specs().exp_subtype, this);
		}
	}
	if (eflags != 0)      set_max_sf(SLOW_SPEED_FACTOR); // can only move at slow speed when an engine is disabled
	if ((time & 15) == 0) homeworld.update_pos();
}


void u_ship::next_frame() {
	
	captured = 0;
	pitch_r  = yaw_r = roll_r = 0.0;
}


void u_ship::set_ship_max_speed(float ms_scale) {

	float max_speed(ms_scale*get_max_speed()); // *fticks
	if (lhyper) {max_speed *= hyperspeed_mult;}
	set_max_speed(max_speed);
}


bool u_ship::collision(point const &copos, vector3d const &vcoll, float obj_mass, float obj_radius, free_obj *source, float elastic) {

	bool no_damage(0), odock(0); // orbiting objects are special - no coll damage

	if (source != NULL) {
		assert(source != this); // no self-collisions
		if (source->get_parent() == this && source->get_time() < W_SELF_ARM_T) return 0; // not yet armed

		if (source->is_ship()) { // ship, not a projectile
			if (source->dock_fighter(this) || source->do_boarding(this)) no_damage = 1;
			if (source->orbital_dock(this)) no_damage = odock = 1;
		}
	}
	if (is_orbiting()) return 1; // no collisions
	if (odock) elastic = 0.0; // no bounce off (stick to the dock?)
	float energy(coll_physics(copos, vcoll, obj_mass, obj_radius, source, elastic));

	if (!no_damage && energy > TOLERANCE) {
		apply_torque_force(copos, vcoll, energy, get_mass());

		if (energy > 1.0) {
			float const e_thresh(0.2*MAX_COLL_DAMAGE);
			if (energy > e_thresh) energy = e_thresh + sqrt(energy - e_thresh);

			if (source != NULL && (is_fighter() || source->is_fighter()) &&
				(parent == source || source->get_parent() == this || (parent != NULL && source->get_parent() == parent)))
			{
				energy *= 0.25; // reduce damage to 25% if parent/fighter or fighter/fighter of the same parent collision
			}
			damage(min(energy, MAX_COLL_DAMAGE), DAMAGE_COLL, copos, source, WCLASS_COLLISION);
		}
	}
	set_ship_max_speed(); // might not be necessary?

	if (specs().suicides && target_obj != NULL && target_obj != parent && source == target_obj) { // explode on impact
		exp_time = time + specs().death_delay; // suicide - no one will get credit for the kill
		exp_val  = armor;
	}
	return 1;
}


float u_ship::damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {

	assert(wc != WCLASS_CAPTURE);
	if (invalid_priv()) return 0.0;
	
	if (val <= 0.0) {
		if (type == DAMAGE_HEAT) partial_uncloak(0.25); // heat (even from tractor beam) decloaks
		return 0.0;
	}
	if (reset_timer > 0) return 0.0; // already destroyed

	if (type != DAMAGE_HEAT) {
		if (val <= specs().damage_abs) return 0.0; // no damage done
		val -= specs().damage_abs;
	}
	if (type == DAMAGE_COLL) {
		if (val <= specs().hull_str) return 0.0; // no damage done
		val -= specs().hull_str;
	}
	bool const valid_wc(wc >= 0 && unsigned(wc) < us_weapons.size());
	bool const ignores_shields(lhyper || (valid_wc && us_weapons[wc].ignores_shields)); // no shields in hyperspeed
	bool const shield_d_only(valid_wc && us_weapons[wc].shield_d_only);
	
	float const hit_points((shield_d_only ? 0.0 : armor) + (ignores_shields ? 0.0 : shields));
	if (VERIFY_REFS && source) source->verify_status();
	bool const is_kill(!shield_d_only && val > hit_points);
	int const src_sclass(source ? source->get_src_sclass() : SWCLASS_UNDEF);
	unsigned const src_align(source ? source->get_align()  : NUM_ALIGNMENT);
	// make sure to report the real damage done (before or after various subtractions?)
	register_damage(src_sclass, sclass, wc, min(val, hit_points), src_align, get_align(), is_kill);
	
	assert(shields >= 0.0 && armor >= 0.0);
	damaged = 1;
	partial_uncloak(0.5); // must de-cloak partially when hit
	float damage_amt(0.0);
	bool const friendly(source == this || ((type == DAMAGE_COLL || alignment == ALIGN_PLAYER) && src_align == alignment)); // accidental collision, etc.
	bool const attack(!friendly && (wc != WCLASS_EXPLODE || (source && source->hostile_explode())) &&
		(wc != WCLASS_HEAT || source->is_ship()) && register_attacker(source));
	if (attack && wc != WCLASS_COLLISION && is_player_ship()) send_warning_message("Under Attack");
	
	if (!ignores_shields) { // determine damage absorbed by shields
		if (shields >= val) {
			shields -= val;
			last_hit = SHIELDS_TIME;
			if (hit_pos != pos) hit_dir = (hit_pos - pos).get_norm();
			return val;
		}
		damage_amt += shields;
		val        -= shields;
		shields     = 0.0;
	}
	if (shield_d_only) return damage_amt; // shield damage only
	
	if (armor >= val) { // shields are down or ignored, armor takes damage, when low ship takes damage
		armor -= val;
		do_structure_damage(val);
		return (damage_amt + val);
	}
	damage_amt += armor;
	val        -= armor;
	armor       = 0.0;
	assert(is_kill);
	if (attack) register_destruction(source);
	if (source && source->source_is_player()) player_ship().acknowledge_kill();  // player only for now
	acknowledge_death();
	destroy_ship(val);
	return damage_amt;
}


void u_ship::destroy_ship(float val) {

	if (specs().regen_delay > 0 && regen_enabled() && !is_player_ship()) {
		reset_after(unsigned((is_orbiting() ? ORBITING_RD_SCALE : 1)*specs().regen_delay));
	}
	exp_time = time + specs().death_delay;
	exp_val  = val;
	
	if (specs().uses_mesh2d) {
		surface_mesh.set_size((2*N_SPHERE_DIV)/3);
		if (specs().mesh_deform) surface_mesh.add_random(0.05, -0.5, 0.2); // surface_mesh.reset_data()?
	}
}


bool u_ship::register_attacker(free_obj const *source) {

	if (source == NULL) return 0;
	bool attack(0);
	free_obj const *damager(source->get_src());
	if (damager == NULL || damager->invalid()) return 0;
	unsigned const d_align(damager->get_align());
	assert(alignment < NUM_ALIGNMENT && d_align < NUM_ALIGNMENT);

	switch (alignment) {
		case ALIGN_NEUTRAL: // attack the shooter except sometimes when GOV
			attack = (d_align != ALIGN_GOV || (rand()&1) == 0);
			break;
		case ALIGN_GOV: // attack the shooter if not GOV
			attack = (d_align != ALIGN_GOV);
			break;
		case ALIGN_PIRATE:
			attack = 1; // always attack
			break;
		case ALIGN_PLAYER:
		default: // ALIGN_RED, ALIGN_BLUE, etc. - attack the shooter, but only attack an ally if there are no enemies
			attack = (d_align != alignment || (target_obj == NULL && last_targ_t >
				(source->is_proj() ? source->get_time() : FFIRE_WAIT_T)) && !is_orbiting() && !damager->is_orbiting());
	}
	if (attack && (!target_set || target_obj == NULL) && target_valid(damager)) { // target_set keeps current target
		if (target_obj == NULL || target_mode == TARGET_ATTACKER || !is_enemy(damager) ||
			((rand()&3) != 0 && p2p_dist_sq(pos, damager->get_pos()) < 1.05*p2p_dist_sq(pos, target_obj->get_pos()))) {
			target_obj = damager; // target our attacker
		}
		if (d_align != ALIGN_GOV && d_align != alignment) { // don't team up on a government ship
			register_attack_from(damager, alignment); // team takes offense
		}
	}
	assert(alignment < attackers.size());
	attackers[alignment] = damager;
	return attack;
}


void u_ship::register_destruction(free_obj const *source) {

	if (source == NULL) return;
	free_obj const *damager(source->get_src());
	if (damager == NULL || damager->invalid()) return;
	unsigned const d_align(damager->get_align());

	if (d_align != ALIGN_GOV) { // don't team up on a government ship
		if (d_align != alignment || !(ai_type & AI_GUARDIAN) || !(damager->get_ai_type() & AI_GUARDIAN)) { // guardians on the same team don't attack each other at all
			register_attack_from(damager, alignment); // team takes offense
		}
	}
}


void u_ship::do_structure_damage(float val) {

	float const ar(((get_max_shields() > 0.0) ? 2.0 : 1.0)*armor/max(1.0f, max(armor, get_max_armor())));
	float const vr(val/(armor + 1.0f));

	if (ar < KILL_CREW_ARMOR && ncrew > specs().req_crew() && vr > 0.1) { // kill some crew
		float const kc_val(SHIP_REQ_CREW*rand_uniform(0.2, 0.5)*min(1.0f, vr));
		ncrew = max((ncrew - unsigned(kc_val*ncrew)), specs().req_crew());
	}
	if (ar < ENGINE_DOWN_ARMOR && specs().nengines > 0 && vr > 0.15) { // disable an engine
		eflags |= (1 << (rand() % specs().nengines)); // engine might already be down
	}
	if (ar < WEAPON_DOWN_ARMOR && !weapons.empty() && vr > 0.12 /*&& (rand()&1) == 0*/) { // disable a weapon
		unsigned const wid(rand() % weapons.size());
		ship_weapon &w(weapons[wid]);
		if (w.ammo > 1) w.ammo -= min(us_weapons[w.wclass].def_ammo, (rand()%(w.ammo>>1))); // ammo is destroyed (up to half)
		
		if (w.wcount > 0 && (w.wcount > 1 || (rand()&1) == 0) && (rand()&1) == 0) {
			--w.wcount; // weapon is damaged until repaired
			++w.ndamaged;
		}
	}
	if (ar < DISABLE_ARMOR && !specs().no_disable) {
		disable(DISABLE_TIME); // disable the ship for awhile
		if (is_player_ship()) disable_player_ship();
	}
}


float u_ship::get_def_explode_damage() const {
	
	return ((specs().exp_type == ETYPE_NONE) ? 0.0 : specs().exp_scale*(0.9f + max(min(0.5f, 0.01f*exp_val), 0.1f)));
}


void u_ship::do_explode() {

	vector3d edir(dir + gen_rand_vector_uniform(0.25));
	edir.normalize();
	int const etype(specs().exp_type);
	def_explode(get_def_explode_damage(), etype, edir, WCLASS_EXPLODE, alignment, EXP_FLAGS_SHIP, parent); // can we blame a ship's explosion on its parent?
	
	if (specs().mass >= 10.0) { // medium to large ship
		float const gain(0.001/max(1E-6f, distance_to_camera(pos)));
		if (gain > 0.001) {gen_sound(SOUND_EXPLODE, pos, gain);}
	}
	if (specs().mass >= 25.0 && (etype == ETYPE_NUCLEAR || etype == ETYPE_STARB)) { // fairly large ship
		colorRGBA ci((etype == ETYPE_NUCLEAR) ? WHITE : YELLOW), co((etype == ETYPE_NUCLEAR) ? BLUE : RED);
		add_uparticle_cloud(pos, 1.0*radius, 4.0*radius, ci, co, colorRGBA(ci, 0.0), colorRGBA(co, 0.0), 10*TICKS_PER_SECOND, 0.0, 0.3, 0.1);
	}
}


void u_ship::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass, int align, unsigned eflags, free_obj const *parent_) {

	assert(parent == NULL || parent_ == NULL || parent == parent_);
	uobject::explode(SHIP_EXP_DAM_SCALE*damage, bradius, etype, edir, exp_time, wclass, align,
		(eflags | EXP_FLAGS_SHIP), parent_); // death may be delayed
	fragment(edir, 1.0, 0);
	surface_mesh.clear();
	detonate_weapons();
}


void u_ship::detonate_weapons() {

	for (unsigned i = 0; i < weapons.size(); ++i) { // ammo ignites
		unsigned const wid(weapons[i].wclass);
		us_weapon const &weap(us_weapons[wid]);
		if (!weap.det_on_exp || weapons[i].ammo == 0) continue;

		for (unsigned n = 0; n < max(1U, weapons[i].ammo/2); ++n) {
			point const ppos(pos + signed_rand_vector(radius));

			if (weap.is_fighter) {
				spawn_fighter(ppos, signed_rand_vector_norm());
			}
			else {
				float const mag(min(0.5*weap.speed, max(0.02*weap.speed, 0.2*velocity.mag())));
				vector3d const pvel(velocity*0.25 + ((mag == 0.0) ? zero_vector : signed_rand_vector(mag)));
				create_projectile(wid, NULL, alignment, ppos, pvel, pvel.get_norm(), upv);
			}
		}
		weapons[i].ammo = 0;
	}
}


void u_ship::fragment(vector3d const &edir, float num, bool outside_cr) const { // send off ship parts

	if (!GEN_FRAGMENTS || (flags & OBJ_FLAGS_DIST)) return;
	static free_obj_allocator<uparticle> allocator;
	unsigned const etype(specs().exp_type);
	float const nscale((etype == ETYPE_NONE) ? 1.5 : 1.0), pscale(0.25*radius), psize(min(pscale, MAX_PARTICLE_SIZE*nscale));
	float const mvscale(max(0.5f, (2.0f*pscale/MAX_PARTICLE_SIZE)));
	float const rscale(is_player_ship() ? 4.0 : 1.0); // move particles further out if player ship so they don't block the screen
	unsigned const modval(min(unsigned(MAX_N_PARTICLES*nscale), max(2U, unsigned(BASE_NUM_PARTS*mvscale)))/2);
	unsigned const num_particles(unsigned(num*(modval + (rand()%modval))));
	bool const nuclear_distrib(num >= 0.5 && specs().exp_type == ETYPE_NUCLEAR);
	assert(etype < NUM_ETYPES);

	for (unsigned i = 0; i < num_particles; ++i) {
		unsigned const plifetime(unsigned(TICKS_PER_SECOND*rand_uniform(20.0, 30.0)));
		float const sz(psize*rand_uniform(0.25, 1.0));
		vector3d dpos(gen_rand_vector_uniform(1.0 + num));
		//if (outside_cr) dpos *= rand_uniform(1.0, 1.2)/dpos.mag();
		point const ppos(pos + dpos*(outside_cr ? c_radius : radius)*rscale);
		vector3d const dv(ppos, pos);
		float const dv_mag(dv.mag());
		vector3d vel(velocity*0.1);
		if (dv_mag > TOLERANCE) {vel + dv*(0.000015/dv_mag);}

		if (nuclear_distrib) {
			vector3d expdir(edir + gen_rand_vector_uniform(0.15));
			if (rand() & 1) expdir.negate();
			vel *= 0.2;
			vel += expdir*(0.002*(0.05 + radius));
		}
		if (etype == ETYPE_STARB) { // parent = this?
			gen_particle(PTYPE_GLOW, et_params[etype].c1, et_params[etype].c2, plifetime/2, ppos, vel, sz, 5.0, alignment, 1);
		}
		else {
			colorRGBA const pcolor((rand()&1) ? specs().base_color : alignment_colors[alignment]);
			unsigned const ptype((rand()&1) ? PTYPE_SPHERE : PTYPE_TRIANGLE);
			gen_particle(ptype, pcolor, pcolor, plifetime, ppos, vel, sz, 0.0, alignment, 1);
		}
	}
}


//#define TIME_SHIP_DRAW


void u_ship::draw_obj(uobj_draw_data &ddata) const { // front is in -z

#ifdef TIME_SHIP_DRAW
	static int times[NUM_US_CLASS] = {0};
	static unsigned count(0), counts[NUM_US_CLASS];
	RESET_TIME;
#endif
	us_class const &sc(specs());
	ddata.nengines = sc.nengines;
	ddata.eflags   = eflags;
	ddata.t_exp    = get_t_exp();
	ddata.on_time  = elapsed_on_t;

	if (sc.uses_tdir && target_dir != zero_vector) {
		ddata.tdir = target_dir;
		if (calc_rvs()) rotate_point(ddata.tdir);
	}
	assert(alignment < NUM_ALIGNMENT);
	colorRGBA color_a(alignment_colors[alignment]);
	colorRGBA color_b(sc.base_color);
	
	if (cloaked > 0.0) {
		assert(sc.has_cloak);
		ddata.cloakval = cloaked;
		ddata.color_a  = ddata.apply_cloak(color_a);
		ddata.color_b  = ddata.apply_cloak(color_b);
	}
	else {
		ddata.color_a  = color_a;
		ddata.color_b  = color_b;
	}
	ddata.set_color(color_a);
	if (sc.engine_lights && ddata.can_have_engine_lights()) {ddata.dlights = 1;} // has engine lights (sc.emits_light?)
	float over_temp(CLIP_TO_01(0.005f*get_over_temp_factor()));
	if (ddata.t_exp > 0.0 && sc.exp_type == ETYPE_FUSION) {over_temp = max(over_temp, (1.0f - ddata.t_exp));} // heats up during explosion
	
	if (over_temp > 0.0) {
		colorRGBA const ecolor((over_temp > 0.5) ? colorRGBA(1.0, (over_temp - 0.5), 0.0, 1.0) : colorRGBA(2.0*over_temp, 0.0, 0.0, 1.0));
		ddata.shader->set_color_e(ecolor);
	}
	if (cloaked < 1.0) {
		if (sclasses[sclass].exp_disint) {ddata.setup_exp_texture();}

		//RESET_TIME;
		switch (sclass) {
		case USC_FIGHTER:    ddata.draw_us_fighter();   break;
		case USC_X1EXTREME:  ddata.draw_x1_extreme();   break;
		case USC_FRIGATE:    ddata.draw_us_frigate();   break;
		case USC_DESTROYER:  ddata.draw_us_destroyer(); break;
		case USC_LCRUISER:   ddata.draw_us_lcruiser();  break;
		case USC_HCRUISER:   ddata.draw_us_hcruiser();  break;
		case USC_BCRUISER:   ddata.draw_us_bcruiser();  break;
		case USC_CARRIER:    ddata.draw_us_carrier();   break;
		case USC_SHADOW:     ddata.draw_us_shadow();    break;
		case USC_ENFORCER:   ddata.draw_us_enforcer();  break;
		case USC_DEFSAT:     ddata.draw_defsat();       break;
		case USC_STARBASE:   ddata.draw_starbase();     break;
		case USC_BCUBE:      ddata.draw_borg(1, 0);     break;
		case USC_BSPHERE:    ddata.draw_borg(0, 0);     break;
		case USC_BTCUBE:     ddata.draw_borg(1, 1);     break;
		case USC_BSPH_SM:    ddata.draw_borg(0, 1);     break;
		case USC_BSHUTTLE:   ddata.draw_bshuttle();     break;
		case USC_TRACTOR:    ddata.draw_tractor();      break;
		case USC_GUNSHIP:    ddata.draw_gunship();      break;
		case USC_NIGHTMARE:  ddata.draw_nightmare();    break;
		case USC_DWCARRIER:  ddata.draw_dwcarrier();    break;
		case USC_DWEXTERM:   ddata.draw_dwexterm();     break;
		case USC_WRAITH:     ddata.draw_wraith();       break;
		case USC_ABOMIN:     ddata.draw_abomination();  break;
		case USC_REAPER:     ddata.draw_reaper();       break;
		case USC_DEATHORB:   ddata.draw_death_orb();    break;
		case USC_SUPPLY:     ddata.draw_supply();       break;
		case USC_ANTI_MISS:  ddata.draw_anti_miss();    break;
		case USC_JUGGERNAUT: ddata.draw_juggernaut();   break;
		case USC_SAUCER:     ddata.draw_saucer(1, 0);   break;
		case USC_SAUCER_V2:  ddata.draw_saucer(0, 0);   break;
		case USC_MOTHERSHIP: ddata.draw_saucer(1, 1);   break;
		case USC_HUNTER:     ddata.draw_headhunter();   break;
		case USC_SEIGE:      ddata.draw_seige();        break;
		case USC_COLONY:     ddata.draw_colony(0,0,0);  break;
		case USC_ARMED_COL:  ddata.draw_colony(1,0,0);  break;
		case USC_HW_COL:     ddata.draw_colony(1,1,0);  break;
		case USC_STARPORT:   ddata.draw_colony(1,0,1);  break;
		case USC_HW_SPORT:   ddata.draw_colony(1,1,1);  break;
		case USC_ARMAGEDDON: ddata.draw_armageddon(surface_mesh); break;
		default: assert(0);
		}
		//if (GET_DELTA_TIME > 10) PRINT_TIME(get_name());
		if (sclasses[sclass].exp_disint) {ddata.end_exp_texture();}
	}
	if (over_temp > 0.0) {ddata.shader->clear_color_e();}

	if (ddata.final_pass && ddata.phase2) {
		glDisable(GL_STENCIL_TEST); // disable in case it's enabled

		if (TEST_SHIP_BOUNDS) {
			ddata.set_color(colorRGBA(WHITE, 0.25));
			draw_bounding_volume(ddata.ndiv);
		}
		if (SHOW_SHIELDS && show_shields()) { // draw shields if recently hit
			fgPopMatrix();
			fgPushMatrix();
			set_additive_blend_mode();
			glDepthFunc(GL_LEQUAL);
			// FIXME: disable alpha testing to avoid artifacts at the shields boundary? but then we have potential alpha sort order problems
			//ddata.shader.add_uniform_float("min_alpha", -1.0);
			assert(last_hit <= SHIELDS_TIME);
			colorRGBA color_alpha(disabled() ? YELLOW : GREEN);
			color_alpha.alpha = (0.5*last_hit)*min(1.0, 2.5*shields/get_max_shields())/SHIELDS_TIME; // decrease at 40% shields
			bool const has_hit_dir(hit_dir != zero_vector);
			int const ndiv(4*min(ddata.ndiv, N_SPHERE_DIV)/3);
			unsigned const ssects(sc.shield_sects);
			float ssize(1.1*c_radius/radius);

			if (ssects > 1) {
				float const ssize0(ssize);
				ssize *= 1.25/ssects;
				vector3d hdir_rot(hit_dir);
				if (has_hit_dir) {rotate_point(hdir_rot);}
				float const z_end(ssize0 - ssize), z_span(2.0*z_end), z_step(z_span/(ssects-1));
				unsigned const sid(unsigned((hdir_rot.z + 1.0)*(ssects-1)*0.5 + 0.5)); // -1.0,1.0 => 0.0,2.0 => 0,ssects-1
				translate_to(dir*(z_end - sid*z_step));
			}
			if (has_hit_dir) { // rotate so that shields appear at hit direction
				select_texture(SBLUR_TEX);
				rotate_sphere_tex_to_dir(hit_dir);
			}
			//if (ssects == 1) {} // scale to create tightly bounding ellipsoid?
			assert(radius > 0.0);
			emissive_shader.enable();
			ddata.set_color(color_alpha);
			draw_sphere_vbo_back_to_front(all_zeros, ssize, 3*ndiv/2, has_hit_dir); // partial sphere?
			glDisable(GL_CULL_FACE);
			glDepthFunc(GL_LESS);
			set_std_blend_mode();
			ddata.shader->enable();
			if (has_hit_dir) {end_texture();}
		} // show shields
	} // final pass/phase
#ifdef TIME_SHIP_DRAW
	unsigned const off((display_mode & 0x08) == 0);
	//display_mode ^= 0x08;
	times[sclass+off] += GET_DELTA_TIME;
	++counts[sclass+off];
	//if (frame_counter == 1000) exit(0); // testing
	
	if (((++count) & 1023) == 0) {
		int tot_time(0);
		unsigned tot_count(0);
		cout << "time count time/count" << endl;

		for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
			if (times[i] > 0) {
				cout << sclasses[i].name << ": " << times[i] << " " << counts[i] << " " << float(times[i])/float(counts[i]) << endl;
				tot_time  += times[i];
				tot_count += counts[i];
			}
		}
		cout << "total: " << tot_time << " " << tot_count << " " << float(tot_time)/float(tot_count) << endl << endl;
	}
#endif
}


// ************ ATTACHED_SHIP ************


void attached_ship::attach(free_obj *obj) {
		
	assert(obj);
	if (attached_to) unattach(); // ???
	attached_to = obj;
	att_pos     = attached_to->get_pos();
	attached_to->add_mass(get_mass());
}


void attached_ship::unattach() {

	assert(attached_to);
	attached_to->add_mass(-get_mass()); // what if mass has changed?
	attached_to = NULL;
	att_pos     = get_pos();
}


void attached_ship::apply_physics() { // might also want to override advance_time

	u_ship::apply_physics();

	if (attached_to) {
		velocity = attached_to->get_velocity();
		pos     += attached_to->get_pos() - att_pos;
		att_pos  = attached_to->get_pos();
	}
}


// ************ MULTIPART_SHIP ************


multipart_ship::multipart_ship(unsigned sclass_, point const &pos0, unsigned align, unsigned ai_type_, unsigned target_mode_, bool rand_orient)
	: u_ship(sclass_, pos0, align, ai_type_, target_mode_, rand_orient), state_val(0.0)
{
	assert(sclass < sclasses.size());
	cobj_vector_t const &sc_cobjs(sclasses[sclass].cobjs);
	cobjs.reserve(sc_cobjs.size());

	for (unsigned i = 0; i < sc_cobjs.size(); ++i) {
		cobjs.add(sc_cobjs[i]->clone()); // will create the correct ship_coll_obj type
		assert(cobjs[i]);
	}
}


void multipart_ship::apply_physics() {

	switch (sclass) { // have to special case this
	case USC_ABOMIN:
		{
			if (!powered() && !cobjs.empty()) break;
			float const aspeed(0.004), aspeed2(0.002), aspeed3(2.8*aspeed2);
			unsigned on_time(get_on_time()); // make more dynamic, controlled by ai_action
			
			if (on_time == 0) {
				seed_on_time(int(1.0/aspeed)); // start at random state
				on_time = get_on_time();
			}
			state_val = CLIP_TO_01(6.0f*fabs(aspeed*(int(on_time)%int(1.0/aspeed) - int(0.5/aspeed))) - 1.0f);
			float const val2(cosf(aspeed2*TWO_PI*(on_time%unsigned(1.0/aspeed2))));
			float const val3(aspeed3*TWO_PI*(on_time%unsigned(1.0/aspeed3))), cv3(0.05*cosf(val3)), sv3(0.05*sinf(val3));
			float const z_start(-2.5), dz(0.2*state_val);
			c_radius = max((z_start + 1.0f - 2.0f*dz), specs().cr_scale);
			cobjs.clear(); // could update or use a custom allocator (37us per abomination)
			cobjs.add(new ship_sphere(point(0.0, 0.0, (z_start + dz)), (1.0 - dz), 0.5)); // main sphere, eye armor
			cobjs.add(new ship_sphere(point(0.0, 0.0, z_start), 0.8, 5.0)); // white part of the eye, more damage
			point tpos[2] = {point(0.0, 0.0, (z_start + 0.15)), all_zeros};

			for (unsigned i = 0; i < 16; ++i) { // add tail cylinders
				float const rv(1.0 - 0.062*i), r1(0.7*rv*rv*rv + 0.25*rv + 0.05);
				float const r2((i == 15) ? 0.0 : (0.92*r1 - 0.01));
				point delta(i*fabs(val2)*cv3, i*fabs(val2)*sv3, (1.0 - 0.01*i*i*fabs(val2)));
				delta *= 0.5/delta.mag();
				tpos[1] = tpos[0] + delta*1.25;
				cobjs.add(new ship_cylinder(tpos[0], tpos[1], r1, r2, 0, 0.25)); // half damage
				c_radius = max(c_radius, max((tpos[0].mag() + SQRT2*r1), (tpos[1].mag() + SQRT2*r2)));
				tpos[0] += delta;
			}
			c_radius *= radius;
			break;
		}

	case USC_REAPER:
		if (cobjs.empty() || (powered() && !invalid_or_disabled())) {
			point p_int(pos + dir*c_radius);

			if (find_coll_enemy_proj(4.0*c_radius, p_int)) {
				// nothing
			}
			else if (get_last_hit_dir() != zero_vector) {
				p_int = pos + get_last_hit_dir()*c_radius;
			}
			else if (target_obj != NULL) {
				p_int = pos + (target_obj->get_pos() - pos).get_norm()*c_radius;
			}
			else if (!cobjs.empty()) {
				break;
			}
			// setup blocking shield
			xform_point(p_int);

			if (cobjs.size() == 2) {
				point const old_pos(cobjs[1]->get_center()/0.8);
				float const angle(get_angle(p_int.get_norm(), old_pos.get_norm()));
				if (fabs(angle) < TOLERANCE) break; // angle hasn't changed
				float const max_angle(2.0*specs().max_turn);

				if (fabs(angle) > max_angle) {
					vector3d const axis(cross_product(p_int, old_pos).get_norm());
					p_int = old_pos;
					rotate_vector3d(axis, ((angle < 0.0) ? -max_angle : max_angle), p_int);
				}
			}
			cobjs.clear();
			cobjs.add(new ship_sphere(all_zeros, 1.0, 1.0)); // main sphere
			cobjs.add(new ship_cylinder(p_int*0.9, p_int*0.7, 0.0, 0.45, 1, 0.0)); // takes no damage
			//weapons[1].weap_pts.resize(0);
			//weapons[1].weap_pts.push_back(p_int*0.9); // lightning
		}
		break;

	default:
		assert(0);
	}
	// change weap_pts?
	u_ship::apply_physics();
}


void multipart_ship::ai_action() {

	switch (sclass) { // have to special case this
	case USC_ABOMIN:
		if (state_val < 0.1) return; // eye mostly closed, can't see, no action
		// write - other stuff, hit with tail?
		break;
	case USC_REAPER:
		// WRITE?
		break;
	default:
		assert(0);
	}
	u_ship::ai_action();
}


void multipart_ship::check_size_scale() { // could update c_radius here instead of in apply_physics()

	assert(size_scale > 0.0);
	radius = size_scale*specs().radius;
}


// ************ ORBITING_SHIP ************


void orbiting_ship::ai_action() {

	assert(alignment < NUM_ALIGNMENT);
	assert(sclasses.size() == NUM_US_CLASS);

	if (begin_motion && sobj_liveable && specs().orbiting_dock && !invalid_or_disabled() && !is_player_ship() &&
		init_credits[alignment] > 0 && (time - last_build_time) > unsigned(2.0*TICKS_PER_SECOND*global_regen))
	{
		unsigned const reserve_credits(sclasses[USC_HW_SPORT].cost + 2*(sclasses[USC_DEFSAT].cost + sclasses[USC_ANTI_MISS].cost));
		vector<unsigned> const &btypes(build_types[alignment]); // hmmm, crashes if just copy btypes...
		bool const excess_credits(team_credits[alignment] > init_credits[alignment]);

		if (team_credits[alignment] >= reserve_credits &&
			((excess_credits && build_any) || ind_ships_used[alignment] < (excess_credits ? 2 : 1)*btypes.size())) // inefficient, but rarely called
		{
			vector<pair<int, unsigned> > scs; // {-cost, type}

			if (excess_credits) { // second game phase build - have excess credits
				set<unsigned> unique_sc;

				if (build_any) {
					for (unsigned i = 0; i < USC_COLONY; ++i) {
						unique_sc.insert(i);
					}
				}
				else {
					copy(btypes.begin(), btypes.end(), inserter(unique_sc, unique_sc.begin()));
				}
				for (set<unsigned>::const_iterator it = unique_sc.begin(); it != unique_sc.end(); ++it) {
					scs.push_back(make_pair(-int(sclasses[*it].cost), *it));
				}
				sort(scs.begin(), scs.end());
			}
			else {
				unsigned const sc(btypes[rand()%btypes.size()]); // choose a random ship type with the same distribution as init
				scs.push_back(make_pair(-int(sclasses[sc].cost), sc));
			}
			for (unsigned i = 0; i < scs.size(); ++i) { // try each suggested build sclass until one can be built
				unsigned const sc(scs[i].second);
				//cout << get_name() << " " << i << ": trying " << sclasses[sc].name << ", cost = " << -scs[i].first << endl;

				if (sclasses[sc].can_move() && alloc_resources_for(sc, alignment, reserve_credits)) {
					point const pos2(pos + dir*(c_radius + sclasses[sc].radius*sclasses[sc].cr_scale));
					u_ship *ship(create_ship(sc, pos2, alignment, AI_ATT_ENEMY, TARGET_CLOSEST, 1));
					assert(ship != NULL); // set parent?
					//cout << get_name() << " built " << ship->get_name() << endl;
					break;
				}
			}
		}
		last_build_time = time; // don't queue up builds - if can't build now then too bad, have to wait another iteration
	}
	u_ship::ai_action();
}


// ************ RAND_SPAWN_SHIP ************


// created at a random location, away from the player but not too far, out of sight
// when out of sight and far away, destroy (but not explode)
// if destroyed, respawn
// any fighters are also type rand_spawn_ship, except they don't respawn
// can't colonize or build anything
rand_spawn_ship::rand_spawn_ship(unsigned sclass_, point const &pos0, unsigned align, unsigned ai_type_, unsigned target_mode_, bool rand_orient, bool will_respawn_)
	: u_ship(sclass_, all_zeros, align, ai_type_, target_mode_, rand_orient), rand_spawn_mixin(pos, radius, rand_spawn_ship_dmax), will_respawn(will_respawn_)
{
	if (pos == all_zeros) {gen_valid_pos();} // otherwise we assume pos is where we want to start
	//cout << "spawn " << get_name() << endl;
}


void rand_spawn_ship::gen_valid_pos() {

	do {
		gen_rand_pos();
	} while (!is_valid_starting_ship_pos(pos, sclass));
}


void rand_spawn_ship::destroy_or_respawn() {

	if (will_respawn) { // respawn
		vector<unsigned> sclasses;
		choose_n_random_sclasses(sclasses, alignment, 1, 0, 1);
		assert(sclasses.size() == 1);
		sclass = sclasses.front();
		reset();
		gen_valid_pos();
		//cout << "respawn " << get_name() << endl;
	}
	else {
		status = 1; // mark for removal
	}
}


void rand_spawn_ship::apply_physics() {

	u_ship::apply_physics();

	if (is_ok() && needs_respawned()) {
		if (okay_to_respawn()) {
			destroy_or_respawn();
		}
		else {
			// what to do here? We have a ship far from the player, but the player is in hyperspeed or has left the system/galaxy
			// for now, just let our ship continue on its course
		}
	}
}


void rand_spawn_ship::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time,
	int wclass, int align, unsigned eflags, free_obj const *parent)
{
	u_ship::explode(damage, bradius, etype, edir, exp_time, wclass, align, eflags, parent);
	if (status == 1) {destroy_or_respawn();} // actually destroyed
}



