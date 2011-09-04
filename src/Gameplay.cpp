// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 5/19/02

#include "gameplay.h"
#include "explosion.h"
#include "transform_obj.h"
#include "player_state.h"
#include "physics_objects.h"
#include "shape_line3d.h"
#include "openal_wrap.h"


bool const SELF_LASER_DAMAGE = 1;
bool const LASER_PATH_LIGHT  = 1; // slow
bool const SMILEY_GAS        = 1;
float const SHADOW_COL_VAL   = 0.5;
float const UWATER_FERR_ADD  = 0.05;
float const UWATER_FERR_MUL  = 2.0;
float const LASER_REFL_ATTEN = 0.95;

// should put these in a file
string const all_smiley_names[] =
{"0-ZERO-0", "Biotch", "Grim Reaper", "Left Nut", "Psycho", "Yomamma", "Mr. Awesome", "IownU", "Fragmaster", "Archangel", "Smiley:)"};


struct text_message_params {

	int mtime, priority;
	float size;
	colorRGBA color;
	bool bitmap;
	text_message_params() : mtime(0), priority(0), size(0.0), color(WHITE), bitmap(0) {}
	text_message_params(int t, float s, colorRGBA const &c, int p, bool b) : mtime(t), priority(p), size(s), color(c), bitmap(b) {}
};


int frags(0), tot_frags(0), best_frags(-1);
int following(0), invalid_collision(0);
int flight(0), blood_spilled(0);
int fired(0), camera_invincible(0), br_source(0), UNLIMITED_WEAPONS(0);
float camera_health(100.0), team_damage(1.0), self_damage(1.0), player_damage(1.0), smiley_damage(1.0);
point orig_camera(all_zeros), orig_cdir(plus_z);
vector<spark_t> sparks;
vector<laser_beam> lasers;
text_message_params msg_params;
string message;
blood_spot blood_spots[NUM_BS];
player_state *sstates = NULL;
team_info *teaminfo = NULL;
vector<bbox> team_starts;


extern int game_mode, window_width, window_height, world_mode, fire_key, spectate, begin_motion, animate2;
extern int camera_reset, frame_counter, camera_mode, camera_coll_id, camera_surf_collide, b2down;
extern int ocean_set, num_groups, island, num_smileys, left_handed, iticks, DISABLE_WATER, spectate;
extern int free_for_all, teams, show_scores, camera_view, xoff, yoff, display_mode;
extern float temperature, ball_velocity, water_plane_z, zmin, zmax, ztop, zbottom, fticks, crater_size;
extern float max_water_height, XY_SCENE_SIZE, czmax, TIMESTEP, atmosphere, camera_shake, base_gravity;
extern point ocean, surface_pos, camera_last_pos;
extern int coll_id[];
extern obj_type object_types[];
extern obj_group obj_groups[];
extern char player_name[];
extern texture_t textures[];
extern vector<rock_shape3d> rock_shapes;
extern vector<coll_obj> coll_objects;
extern GLUquadricObj* quadric;



point get_sstate_pos(int id) {
	return ((id == CAMERA_ID) ? get_camera_pos() : obj_groups[coll_id[SMILEY]].get_obj(id).pos);
}

vector3d get_sstate_dir(int id) {
	return ((id == CAMERA_ID) ? cview_dir : obj_groups[coll_id[SMILEY]].get_obj(id).orientation);
}

float get_sstate_radius(int id) {
	return ((id == CAMERA_ID) ? CAMERA_RADIUS : object_types[SMILEY].radius);
}


int gen_game_obj(int type) {

	switch (type) {
		case POWERUP:
			//return PU_SPEED;
			return rand()%NUM_POWERUPS;
		case WEAPON:
			while (1) {
				//return W_LASER;
				int const id(rand()%NUM_WEAPONS);
				if (weapons[id].need_weapon) return id;
			}
		case AMMO:
			while (1) {
				//return W_LASER;
				int const id(rand()%NUM_WEAPONS);
				if (weapons[id].need_ammo) return id;
			}
		default: assert(0);
	}
	return -1;
}


int get_ammo_or_obj(int wid) {

	switch (wid) {
	case W_M16:     return SHELLC;
	case W_SHOTGUN: return PROJECTILE;
	case W_PLASMA:  return PLASMA;
	case W_LASER:   return BEAM;
	}
	return weapons[wid].obj_id; // could be -1 for invalid
}


int wid_need_weapon(int wid) { // used in draw_world.cpp
	return weapons[wid].need_weapon;
}

bool is_area_damage(int type) {
	return (type == DROWNED || type == FROZEN || type == SUFFOCATED || type == GASSED);
}


int same_team(int source, int target) {

	if (teams <= 1) return (!((source == CAMERA_ID) ^ (target == CAMERA_ID)));
	int source_team((source + teams)%teams);
	int target_team((target + teams)%teams);
	return (source_team == target_team);
}


void blood_on_camera(unsigned num_spots) {

	num_spots = min(num_spots, NUM_BS/4);
	float const xrand(0.012*((float)window_width/(float)window_height));

	for (unsigned i = 0, n = 0; i < NUM_BS && n < num_spots; ++i) {
		if (blood_spots[i].time <= 0 || ((NUM_BS - i) <= (num_spots - n))) {
			blood_spots[i].pos.assign(rand_uniform(-xrand, xrand), rand_uniform(-0.012, 0.012), -0.02);
			blood_spots[i].size = rand_uniform(3.0, 50.0);
			blood_spots[i].time = int(10 + 3*blood_spots[i].size);
			++n;
		}
	}
}


int compute_damage(float &energy, int type, int obj_index, int source, int target) {

	energy = max(energy, 0.0f);

	if (type >= 0 && type < NUM_TOT_OBJS && (type != BALL || obj_groups[coll_id[type]].get_obj(obj_index).status == 1 ||
		obj_groups[coll_id[type]].get_obj(obj_index).status == 2))
	{
		energy += object_types[type].damage;
	}
	if (type == SHIELD) sstates[target].shields = min(MAX_SHIELDS, (100.0f + sstates[target].shields));
	if (type == HEALTH || type == SHIELD || type == POWERUP) return 1;

	if (source == target && energy > 0.0) { // hit yourself
		if (self_damage == 0.0) return 0;
		if (type == BALL) energy = 0.0;
		else energy *= self_damage;
	}
	if (energy == 0.0) return 1;
	blood_spilled = 1;

	if (type == STAR5 || type == SHRAPNEL) {
		energy *= 0.8*fticks;
	}
	else if (type == FRAGMENT) {
		energy *= obj_groups[coll_id[type]].get_obj(obj_index).vdeform.x;
	}
	if (sstates[target].powerup == PU_SHIELD) {
		if (source == target && (type != LANDMINE && type != FELL && type != DROWNED && type != CRUSHED)) return 0;
		energy *= sstates[target].get_shield_scale();
	}
	if (type != FELL && type != DROWNED) energy *= sstates[source].get_damage_scale();
	energy *= ((target == CAMERA_ID) ? player_damage : smiley_damage);

	if (source != target && same_team(source, target)) {
		if (team_damage == 0.0) return 0;
		energy *= team_damage; // hit by a teammate
	}
	float const shield_damage(min(0.75f*HEALTH_PER_DAMAGE*energy, sstates[target].shields));
	sstates[target].shields -= shield_damage;
	energy -= shield_damage/HEALTH_PER_DAMAGE;
	return 1;
}


int self_coll_invalid(int type, int obj_index) {

	if (type == ROCKET || type == SEEK_D || type == PROJECTILE || type == LASER || type == STAR5 || type == GASSED) {
		if (type == ROCKET || type == SEEK_D || type == STAR5) invalid_collision = 1;
		return 1;
	}
	if ((type == GRENADE || type == CGRENADE || type == S_BALL || type == BALL || type == PLASMA || type == SHRAPNEL) &&
		obj_groups[coll_id[type]].get_obj(obj_index).time < 10)
	{
		invalid_collision = 1;
		return 1;
	}
	if ((type == STAR5 || (type == SHRAPNEL && obj_groups[coll_id[type]].get_obj(obj_index).angle == 0.0)) &&
		obj_groups[coll_id[type]].get_obj(obj_index).velocity.mag_sq() > 1.0)
	{
		invalid_collision = 1;
		return 1;
	}
	if (type == LANDMINE && obj_groups[coll_id[type]].get_obj(obj_index).time < (int)SMILEY_LM_ACT_TIME) {
		invalid_collision = 1;
		return 1;
	}
	return 0;
}


// smileys and camera
void gen_dead_smiley(int source, int target, float energy, point const &pos, vector3d const &velocity, vector3d const &coll_dir,
					 int damage_type, float health, float radius, float blood_v, bool burned, int type)
{
	// eyes, nose, and tongue
	int const pcid(coll_id[SFPART]), offset(4*(target+1));
	int end_p_loop(3);
	float part_v(5.0 + 0.5*sqrt(energy)), chunk_v(7.5 + 0.2*sqrt(energy));
	vector3d orient(0.0, 0.0, 0.0);
	player_state &sstate(sstates[target]);

	if (target == CAMERA_ID) {
		orient = cview_dir;
	}
	else {
		assert(target < num_smileys);
		orient = obj_groups[coll_id[SMILEY]].get_obj(target).orientation; // should be normalized
	}
	if (sstate.kill_time < int(2*TICKS_PER_SECOND) || sstate.powerup == PU_DAMAGE) {
		++end_p_loop; // tongue out
	}
	assert(unsigned(offset + end_p_loop) <= obj_groups[pcid].max_objs);

	for (int i = 0; i < end_p_loop; ++i) { // {eye, eye, nose, [tongue]}
		obj_groups[pcid].create_object_at((i+offset), pos);
		dwobject &obj(obj_groups[pcid].get_obj(i+offset));
		obj.pos += gen_rand_vector(radius*rand_uniform(0.7, 1.3), 1.0, PI);
		obj.direction = (i <= 1 ? 0 : i-1);
		gen_blood_velocity(obj.velocity, velocity, coll_dir, part_v, 0.3, 0.2, damage_type, health);
		if (i == 3) obj.orientation = orient; // tongue
	}

	// chunks
	int const chunk_start(NUM_CHUNK_BLOCKS*(SMILEY_NCHUNKS*(target+1) + sstate.chunk_index));
	int const ccid(coll_id[CHUNK]), end_c_loop(chunk_start + SMILEY_NCHUNKS);
	assert(unsigned(end_c_loop) <= obj_groups[ccid].max_objs);
	if (burned) chunk_v *= 0.4;

	for (int i = chunk_start; i < end_c_loop; ++i) { // chunks
		obj_groups[ccid].create_object_at(i, pos);
		dwobject &obj(obj_groups[ccid].get_obj(i));
		obj.pos       += gen_rand_vector(radius*rand_uniform(0.7, 1.3), 1.0, PI);
		obj.init_dir   = signed_rand_vector_norm();
		obj.init_dir.z = fabs(obj.init_dir.z);
		gen_blood_velocity(obj.velocity, velocity, coll_dir, chunk_v, 0.25, 0.22, damage_type, health);
		float const vmag(obj.velocity.mag());
		if (vmag > TOLERANCE) obj.pos += obj.velocity*(radius/vmag);
		if (burned) obj.flags |= TYPE_FLAG;
	}

	// skull
	if (burned || energy < 250.0) {
		obj_groups[coll_id[SKULL]].create_object_at(target+1, pos);
		dwobject &obj(obj_groups[coll_id[SKULL]].get_obj(target+1));
		obj.orientation = orient;
		gen_blood_velocity(obj.velocity, velocity, coll_dir, part_v, 0.15, 0.1, damage_type, health);
	}
	if (burned) {
		gen_fire(pos, 1.0, source);
		gen_smoke(pos);
	}
	else { // add blood
		add_color_to_landscape_texture(BLOOD_C, pos.x, pos.y, min(4.0, double(sqrt(blood_v)))*radius, 0);
	}
	sstate.chunk_index = (sstate.chunk_index + 1) % NUM_CHUNK_BLOCKS;
	
	if (type == DROWNED) {
		gen_sound(SOUND_DROWN, pos);
	}
	else if (type == CRUSHED || type == FELL) {
		gen_sound(SOUND_SQUISH, pos);
	}
	else if (burned) {
		gen_sound(SOUND_SCREAM2, pos, 1.0, rand_uniform(0.7, 1.3));
	}
	else {
		gen_sound(SOUND_SCREAM1, pos, 1.0, rand_uniform(0.7, 1.3));
	}
}


unsigned create_blood(int index, int amt_denom, point const &pos, float obj_radius, vector3d const &velocity,
					  vector3d const &coll_dir, float blood_v, int damage_type, float health, bool burned)
{
	assert(amt_denom > 0);
	int const cid(coll_id[burned ? CHARRED : BLOOD]);
	obj_group &objg(obj_groups[cid]);
	unsigned const blood_amt(objg.max_objs/(obj_groups[coll_id[SMILEY]].max_objs+1));
	unsigned const start_ix(blood_amt*index);
	unsigned const end_ix(min(objg.max_objs, unsigned(blood_amt*(index+1))));

	for (unsigned i = (start_ix + rand()%amt_denom); i < end_ix; i += amt_denom) {
		objg.create_object_at(i, pos);
		objg.get_obj(i).pos += gen_rand_vector(obj_radius*rand_uniform(0.7, 1.3), 1.0, PI);
		gen_blood_velocity(objg.get_obj(i).velocity, velocity, coll_dir, blood_v, 0.3, 0.3, damage_type, health);
	}
	return blood_amt;
}


inline float get_shrapnel_damage(float energy, int index) {

	return 0.5*energy + object_types[SHRAPNEL].damage*(1.0 -
		((float)obj_groups[coll_id[SHRAPNEL]].get_obj(index).time)/((float)object_types[SHRAPNEL].lifetime));
}


bool player_state::pickup_ball(int index) {

	if (game_mode != 2)    return 0;
	dwobject &obj(obj_groups[coll_id[BALL]].get_obj(index));
	if (obj.disabled())    return 0; // already picked up this frame?
	if (UNLIMITED_WEAPONS) return 0; //obj.disable()?
	assert(p_ammo[W_BALL] == balls.size());
	weapon            = W_BALL;
	p_weapons[W_BALL] = 1;
	++p_ammo[W_BALL];
	balls.push_back(index);
	obj.status = OBJ_STAT_RES; // reserved status
	return 1;
}


bool is_burned(int type, int br_source) {

	return (type == PLASMA || type == FIRE || type == BURNED || type == LASER || type == BEAM ||
		(type == BLAST_RADIUS && br_source == PLASMA));
}


void player_coll(int type, int obj_index) {

	if (type == STAR5) {
		obj_groups[coll_id[type]].get_obj(obj_index).time += int(200.0*fticks);
	}
	else if (type == SHRAPNEL) {
		obj_groups[coll_id[type]].get_obj(obj_index).time += int(10.0*fticks);
	}
}


void update_kill_health(float &health) {

	health = max(health, min(100.0f, (health + KILL_HEALTH)));
}


// ***********************************
// COLLISION DAMAGE CODE
// ***********************************


int proc_coll_types(int type, int obj_index, float &energy) {

	int const ocid(coll_id[type]);

	if (type == WEAPON || type == AMMO || type == WA_PACK || type == POWERUP || type == HEALTH || type == SHIELD) {
		obj_groups[ocid].get_obj(obj_index).disable();
		energy = 0.0; // no damage done
	}
	else if ((type == GRENADE || type == CGRENADE) && (rand()&3) == 0) {
		obj_groups[ocid].get_obj(obj_index).time = object_types[type].lifetime; // maybe explode on collision
	}
	if (type == POWERUP || type == WEAPON || type == AMMO || type == WA_PACK) {
		return (int)obj_groups[ocid].get_obj(obj_index).direction;
	}
	return 0;
}


void camera_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	if (type == CAMERA || type == SMILEY || invalid_collision || !camera_mode || !game_mode || spectate ||
		(type != SMILEY && !damage_done(type, obj_index)))
	{
		return;
	}
	if (camera_health < 0.0) return; // already dead
	int const source(get_damage_source(type, obj_index, CAMERA_ID));
	assert(source >= CAMERA_ID);
	if (source == CAMERA_ID && self_coll_invalid(type, obj_index)) return; // hit yourself
	int damage_type(0);
	player_state &sstate(sstates[CAMERA_ID]);
	colorRGBA cam_filter_color(RED);
	float const last_h(camera_health);
	int const wa_id(proc_coll_types(type, obj_index, energy));

	switch (type) {
	case POWERUP:
		sstate.powerup      = wa_id;
		sstate.powerup_time = POWERUP_TIME;
		cam_filter_color    = WHITE;
		print_text_onscreen(powerup_names[wa_id], GREEN, 1.0, 2*MESSAGE_TIME/3, 1);
		gen_sound(SOUND_POWERUP, position);
		break;

	case HEALTH:
		cam_filter_color = GREEN;
		print_text_onscreen("+50 Health", GREEN, 1.0, 2*MESSAGE_TIME/3, 1);
		gen_sound(SOUND_ITEM, position, 0.5);
		break;

	case SHIELD:
		cam_filter_color = YELLOW;
		print_text_onscreen("+100 Shields", GREEN, 1.0, 2*MESSAGE_TIME/3, 1);
		gen_sound(SOUND_ITEM, position, 0.5);
		break;

	case WEAPON:
		sstate.p_weapons[wa_id] = 1;
		sstate.p_ammo[wa_id]    = min(weapons[wa_id].max_ammo, sstate.p_ammo[wa_id]+weapons[wa_id].def_ammo);
		print_text_onscreen(weapons[wa_id].name, GREEN, 0.8, 2*MESSAGE_TIME/3, 1);
		cam_filter_color = BLUE;
		gen_sound(SOUND_ITEM, position, 0.25);
		break;

	case AMMO:
		sstate.p_ammo[wa_id] = min(weapons[wa_id].max_ammo, sstate.p_ammo[wa_id]+weapons[wa_id].def_ammo);
		print_text_onscreen((make_string(weapons[wa_id].def_ammo) + " " + weapons[wa_id].name + " ammo"), GREEN, 0.8, 2*MESSAGE_TIME/3, 1);
		cam_filter_color = BLUE;
		gen_sound(SOUND_ITEM, position, 0.25);
		break;

	case WA_PACK:
		{
			int const pickup_ammo((int)obj_groups[coll_id[type]].get_obj(obj_index).angle);
			sstate.p_weapons[wa_id] = 1;
			sstate.p_ammo[wa_id]    = min((int)weapons[wa_id].max_ammo, (sstate.p_ammo[wa_id] + pickup_ammo));
			print_text_onscreen((weapons[wa_id].name + " pack with ammo " + make_string(pickup_ammo)), GREEN, 0.8, 2*MESSAGE_TIME/3, 1);
			cam_filter_color = BLUE;
			gen_sound(SOUND_ITEM, position, 0.25);
		}
		break;

	case BALL:
		if (energy < 10.0 && sstate.pickup_ball(obj_index)) {
			print_text_onscreen("You have the ball", GREEN, 1.2, 2*MESSAGE_TIME/3, 1);
			gen_sound(SOUND_POWERUP, position, 0.5);
		}
		else cam_filter_color = RED;
		break;

	case LANDMINE:
		damage_type      = 1;
		cam_filter_color = RED;
		break;
	case SHRAPNEL:
		energy = get_shrapnel_damage(energy, obj_index);
		cam_filter_color = RED;
		break;
	case BLAST_RADIUS:
		if (br_source == LANDMINE) damage_type = 1;
		cam_filter_color = RED;
		break;
	case GASSED:
		print_text_onscreen("Poison Gas Detected", OLIVE, 1.2, MESSAGE_TIME, 1);
		cam_filter_color = ((frame_counter & 15) ? DK_RED : OLIVE);
		break;
	default:
		cam_filter_color = RED;
	}
	if (!compute_damage(energy, type, obj_index, source, CAMERA_ID)) return;
	if (energy > 0.0 && camera_invincible) return;
	camera_health -= HEALTH_PER_DAMAGE*energy;
	camera_health = min(camera_health, ((sstate.powerup == PU_REGEN) ? MAX_REGEN_HEALTH : max(last_h, MAX_HEALTH)));
	bool const is_blood(energy > 0.0 && !is_area_damage(type));
	player_coll(type, obj_index);
	point const camera(get_camera_pos());
	vector3d const coll_dir(get_norm_rand(vector3d(position, camera)));
	bool const burned(is_burned(type, br_source)), alive(camera_health >= 0.0);
	float const blood_v((energy > 0.0) ? (6.0 + 0.6*sqrt(energy)) : 0.0);
	if (is_blood) create_blood(0, (alive ? 30 : 1), camera, CAMERA_RADIUS, velocity,
		coll_dir, blood_v, damage_type, camera_health, burned);

	if (alive) {
		if (is_blood && cam_filter_color == RED) {
			if (sstate.drop_weapon(coll_dir, cview_dir, camera, CAMERA_ID, energy, type)) {
				print_text_onscreen("Oops, you dropped your weapon!", RED, 1.2, 2*MESSAGE_TIME/3, 3);
			}
			if (camera_shake == 0.0) camera_shake = 1.0;
		}
		if (energy > 1.0E-4 || cam_filter_color != RED) {
			cam_filter_color.alpha = 0.001*fabs(energy) + 0.25;
		}
		else {
			cam_filter_color.alpha = 0.0;
		}
		if (type == FELL) gen_sound(SOUND_SQUISH, camera, 0.3);
	}
	else { // dead
		if (!spectate) {
			camera_mode = !camera_mode;
			reset_camera_pos();
		}
		sstate.powerup         = -1;
		sstate.powerup_time    = 0;
		camera_health          = 100.0;
		cam_filter_color.alpha = 1.0;
		remove_reset_coll_obj(camera_coll_id);
		gen_dead_smiley(source, CAMERA_ID, energy, camera, velocity, coll_dir, damage_type, camera_health, CAMERA_RADIUS, blood_v, burned, type);
		assert(obj_groups[coll_id[SFPART]].max_objects() > 0);
		obj_groups[coll_id[SFPART]].get_obj(0).flags |= CAMERA_VIEW; // camera follows the eye
		orig_cdir = cview_dir;
		
		if (type == SMILEY) {
			string const str(string("Fragged by ") + sstates[source].name);

			if (same_team(source, CAMERA_ID)) {
				--sstates[source].kills;
				print_text_onscreen(str, RED, 1.0, MESSAGE_TIME, 3); // killed by your teammate
			}
			else {
				++sstates[source].kills;
				sstates[source].kill_time = 0;
				print_text_onscreen(str, ORANGE, 1.0, MESSAGE_TIME, 3); // killed by an enemy
			}
		}
		else {
			string str;

			if (source == CAMERA_ID) { // camera/player
				cam_filter_color = BLACK;

				switch (type) {
					case FELL:       str = "FELL";             break;
					case DROWNED:    str = "DROWNED";          break;
					case FIRE:       str = "SUICIDE by FIRE";  break;
					case BURNED:     str = "BURNED to DEATH";  break;
					case FROZEN:     str = "FROZE to DEATH";   break;
					case SUFFOCATED: str = "SUFFOCATED";       break;
					case CRUSHED:    str = "CRUSHED to DEATH"; break;
					case GASSED:     str = "Gassed to Death";  break;
					default:         str = string("SUICIDE with ") + obj_type_names[type];
				}
				print_text_onscreen(str, RED, 1.0, MESSAGE_TIME, 3);
				--sstates[CAMERA_ID].kills;
				--frags; // suicide
			}
			else {
				if (type == FIRE || type == BURNED) {
					str = string("BURNED to DEATH by ") + sstates[source].name;
				}
				else {
					str = string("Fragged by ") + sstates[source].name + "'s " +
						get_weapon_qualifier(type, index, source) + " " + obj_type_names[type];
				}
				if (same_team(source, CAMERA_ID)) {
					--sstates[source].kills;
					print_text_onscreen(str, RED, 1.0, MESSAGE_TIME, 3); // killed by your teammate
				}
				else {
					++sstates[source].kills;
					sstates[source].kill_time = 0;
					print_text_onscreen(str, ORANGE, 1.0, MESSAGE_TIME, 3); // killed by an enemy
				}
			}
		}
		if (!same_team(CAMERA_ID, source) && obj_groups[coll_id[SMILEY]].is_enabled()) {
			update_kill_health(obj_groups[coll_id[SMILEY]].get_obj(source).health);
		}
		if (frags > best_frags) best_frags = frags;
		if (is_blood && !burned) blood_on_camera(rand()%10);
		sstate.drop_pack(camera);
		remove_reset_coll_obj(camera_coll_id);
		init_sstate(CAMERA_ID, 0);
		sstate.killer = source;
		frags         = 0;
		camera_reset  = 0;
		b2down        = 0;
		following     = 0; // ???
		orig_camera   = camera_origin;
		++sstates[CAMERA_ID].deaths;
	}
	if (cam_filter_color.alpha > 0.0) add_camera_filter(cam_filter_color, CAMERA_SPHERE_TIME, -1, CAM_FILT_DAMAGE);
}


void smiley_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	if (index == CAMERA_ID) {
		camera_collision(type, obj_index, velocity, position, energy, type);
		return;
	}
	if (type == CAMERA) {
		camera_collision(type, index, velocity, position, energy, SMILEY);
		return;
	}
	int damage_type(0), cid(coll_id[SMILEY]);
	assert(cid >= 0);
	assert(obj_groups[cid].enabled);
	assert(index >= 0 && index < num_smileys);
	if (invalid_collision || (!game_mode || type == SMILEY)) return;
	if (obj_groups[cid].get_obj(index).disabled() || obj_groups[cid].get_obj(index).health < 0.0) return;
	if (!damage_done(type, obj_index)) return;
	int const source(get_damage_source(type, obj_index, index));
	assert(source >= CAMERA_ID);
	if (source == index && self_coll_invalid(type, obj_index)) return; // hit itself
	player_coll(type, obj_index);
	int const wa_id(proc_coll_types(type, obj_index, energy));
	player_state &sstate(sstates[index]);

	switch (type) {
	case POWERUP:
		sstate.powerup      = wa_id;
		sstate.powerup_time = POWERUP_TIME;
		print_text_onscreen((sstates[index].name + " has " + powerup_names[wa_id]),
			get_smiley_team_color(index), 0.8, 2*MESSAGE_TIME, 0);

	case HEALTH:
		obj_groups[cid].get_td()->get_mesh(index).mult_by(0.4);

	case SHIELD:
		energy = 0.0;
		break;

	case WEAPON:
		sstate.p_weapons[wa_id] = 1;
		sstate.p_ammo[wa_id]    = min(weapons[wa_id].max_ammo, sstate.p_ammo[wa_id]+weapons[wa_id].def_ammo);
		if (sstate.weapon == W_BBBAT || sstate.weapon == W_SBALL || rand()%10 > 4) {
			sstate.weapon = wa_id; // switch weapons
		}
		break;

	case AMMO:
		sstate.p_ammo[wa_id] = min(weapons[wa_id].max_ammo, sstate.p_ammo[wa_id]+weapons[wa_id].def_ammo);
		if ((!weapons[wa_id].need_weapon || (sstate.p_weapons[wa_id] != 0 && sstate.p_ammo[wa_id] == 0)) &&
			(sstate.weapon == W_BBBAT || wa_id != W_SBALL) &&
		    (sstate.weapon == W_BBBAT || sstate.weapon == W_SBALL || rand()%10 > 5)) {
			sstate.weapon = wa_id; // switch weapons
		}
		break;

	case WA_PACK:
		sstate.p_weapons[wa_id] = 1;
		sstate.p_ammo[wa_id]    = min((int)weapons[wa_id].max_ammo,
			sstate.p_ammo[wa_id]+(int)obj_groups[coll_id[type]].get_obj(obj_index).angle);
		if ((sstate.weapon == W_BBBAT || wa_id != W_SBALL) &&
			(sstate.weapon == W_BBBAT || sstate.weapon == W_SBALL || rand()%10 > 6)) {
			sstate.weapon = wa_id; // switch weapons
		}
		break;

	case BALL:
		if (energy < 10.0) sstate.pickup_ball(obj_index);
		break;
	case LANDMINE:
		damage_type = 1;
		break;
	case SHRAPNEL:
		energy = get_shrapnel_damage(energy, obj_index);
		break;
	case BLAST_RADIUS:
		if (br_source == LANDMINE) damage_type = 1;
		break;
	}
	if (!compute_damage(energy, type, obj_index, source, index)) return;
	dwobject &obji(obj_groups[cid].get_obj(index));
	obji.health = min((obji.health - HEALTH_PER_DAMAGE*energy), (sstate.powerup == PU_REGEN) ? MAX_REGEN_HEALTH : MAX_HEALTH);
	int const alive(obji.health >= 0.0);
	if (energy <= 0.0 && alive) return;
	bool const burned(is_burned(type, br_source));
	float const radius(object_types[SMILEY].radius);
	point const obj_pos(obji.pos);
	float blood_v(6.0 + 0.6*sqrt(energy));
	vector3d coll_dir(get_norm_rand(vector3d(position, obj_pos)));

	if (alive) {
		if (type != FELL && type != CRUSHED && !is_area_damage(type)) {
			if (sstate.shields < 0.01) {
				float const damage_amt(max(2.0, min(8.0, 0.8*sqrt(energy))));
				add_damage_to_smiley(coll_dir, damage_amt, index, type);
			}
			if (energy > 0.1 && source != index && coll_dir != zero_vector && type != FIRE && type != BURNED &&
				type != ROCK && (type != SHRAPNEL || energy > 10.0*object_types[SHRAPNEL].damage))
			{
				sstate.was_hit = HIT_TIME;
				sstate.hitter  = source;
				sstate.hit_dir = coll_dir;
			}
		}
		sstate.drop_weapon(coll_dir, vector3d(sstate.target_pos, obj_pos), obj_pos, index, energy, type);
		blood_v  *= 0.5;
		coll_dir *= -4.0;
		if (type == FELL) gen_sound(SOUND_SQUISH, obj_pos, 0.2);
	}
	if (!burned && !is_area_damage(type)) {
		unsigned const blood_amt(create_blood(index+1, (alive ? 30 : 1), obj_pos, radius,
			velocity, coll_dir, blood_v, damage_type, obji.health, burned));
		float const cdist(distance_to_camera(obj_pos));

		if (cdist < 4.0*radius && sphere_in_camera_view(obj_pos, radius, 2)) {
			int const nblood(max(1, int(((cdist + 0.1)/radius)*blood_amt*max(0.5f, min(1.0f, (-obji.health/50.0f + 0.2f)))/2.0f)));
			blood_on_camera(rand()%(nblood+1));
			gen_sound(SOUND_SPLAT1, get_camera_pos());
		}
	}
	if (alive) return;

	// dead
	sstate.powerup      = -1;
	sstate.powerup_time = 0;
	gen_dead_smiley(source, index, energy, obj_pos, velocity, coll_dir, damage_type, obji.health, radius, blood_v, burned, type);
	player_state &ssource(sstates[source]);

	if (source == CAMERA_ID) { // camera/player
		if (camera_mode == 1 && camera_surf_collide) { // playing
			int type2(type);
			int const weapon((type == BLAST_RADIUS ? br_source : ssource.weapon));
			bool const head_shot(type != BLAST_RADIUS && type != FIRE && type != LANDMINE && !is_area_damage(type));
			if (game_mode == 2 && type == CAMERA) type2 = BALL; // dodgeball
			if (!same_team(index, source)) update_kill_health(camera_health);
			string const str(string("You fragged ") + sstates[index].name + " with " + get_weapon_qualifier(type, weapon, source) +
				" " + obj_type_names[type2] + (head_shot ? "\nHead Shot" : ""));

			if (same_team(source, index)) { // player killed a teammate
				print_text_onscreen(str, RED, 1.0, MESSAGE_TIME, 2);
				--frags;
				--sstates[CAMERA_ID].kills;
				--tot_frags;
			}
			else { // player killed an enemy
				print_text_onscreen(str, MAGENTA, 1.0, MESSAGE_TIME, 2);
				++frags;
				++sstates[CAMERA_ID].kills;
				++tot_frags;
				if (frags > best_frags) best_frags = frags; // not sure this is correct with suicides and team kills
			}
		}
	}
	else if (source == index) { // suicide
		assert(source < num_smileys);
		string str(sstates[source].name);

		switch (type) {
		case DROWNED:    str += " Drowned";              break;
		case BURNED:     str += " Burned to Death";      break;
		case FROZEN:     str += " Froze to Death";       break;
		case FELL:       str += " Fell to his Death";    break;
		case SUFFOCATED: str += " Suffocated";           break;
		case CRUSHED:    str += " was Crushed to Death"; break;
		case GASSED:     str += " was Gassed to Death";  break;
		default:         str += " suicided with " +
							 get_weapon_qualifier(type, (type == BLAST_RADIUS ? br_source : ssource.weapon), source) + " " + obj_type_names[type];
		}
		print_text_onscreen(str, CYAN, 1.0, MESSAGE_TIME/2, 0);
		//--ssource.kills;
	}
	else {
		assert(source < num_smileys && index < num_smileys);
		string const str(sstates[index].name + " was fragged by " + sstates[source].name + "'s " +
			get_weapon_qualifier(type, (type == BLAST_RADIUS ? br_source : ssource.weapon), source) + " " + obj_type_names[type]);
		
		if (same_team(source, index)) { // killed a teammate
			print_text_onscreen(str, PINK, 1.0, MESSAGE_TIME/2, 0);
			--ssource.kills;
		}
		else { // killed an enemy
			print_text_onscreen(str, YELLOW, 1.0, MESSAGE_TIME/2, 0);
			if (free_for_all) {
				++ssource.kills;
				ssource.kill_time = 0;
			}
		}
		if (!same_team(index, source)) update_kill_health(obj_groups[coll_id[SMILEY]].get_obj(source).health);
	}
	sstate.drop_pack(obj_pos);
	remove_reset_coll_obj(obji.coll_id);
	++sstate.deaths;
	sstate.killer = source;
	obji.status   = 0;
	if (game_mode != 2) gen_smoke(position);
}


int get_smiley_hit(vector3d &hdir, int index) {

	assert(index < num_smileys);
	if (sstates[index].was_hit == 0) return 0;
	hdir = sstates[index].hit_dir;
	return sstates[index].was_hit;
}


string get_weapon_qualifier(int type, int index, int source) {

	bool const qd(sstates[source].powerup == PU_DAMAGE); // quad damage
	if (index == CAMERA_ID) index = sstates[CAMERA_ID].weapon; // camera weapon
	string qstr;
	
	if ((type == IMPACT && index != IMPACT) || (type == PROJECTILE && index != PROJECTILE)) {
		assert(index >= 0 && index < NUM_WEAPONS);
		qstr = weapons[index].name;
	}
	else if (type == BLAST_RADIUS && index != BLAST_RADIUS) {
		assert(index >= 0 && index < NUM_TOT_OBJS);
		qstr = obj_type_names[index];
	}
	if (qd) return string("quad damage ") + qstr;
	return qstr;
}


bool dwobject::lm_coll_invalid() const {

	return (time < LM_ACT_TIME);
}


int damage_done(int type, int index) {

	int const cid(coll_id[type]);

	if (cid >= 0) {
		dwobject &obj(obj_groups[cid].get_obj(index));
		if ((obj.flags & (WAS_PUSHED | FLOATING)) && (type != BALL || game_mode != 2)) return 0; // floating on the water or pushed after stopping, can no longer do damage

		if (type == LANDMINE) {
			if (obj.status == 1 || obj.lm_coll_invalid()) return 0; // not activated
			obj.status = 0;
			return 1;
		}
	}
	return damage_done_obj[type];
}


void gen_blood_velocity(vector3d &vout, vector3d const &velocity, vector3d const &coll_dir, float blood_v, float md, float mv, int type, float health) {

	assert(!is_nan(coll_dir));
	float const hv(max(0.7, min(1.1, (-health/40.0 + 0.25)))), mag(rand_uniform(0.5*blood_v, blood_v));
	vout   = gen_rand_vector(mag, 2.0, 0.52*PI);
	vout.z = fabs(vout.z); // always up

	for (unsigned i = 0; i < 3; ++i) {
		vout[i] = hv*(-md*blood_v*coll_dir[i] + mv*velocity[i] + vout[i]);
		if (type == 1 && i < 2) vout[i] *= 0.2; // x and y
		assert(!is_nan(vout[i]));
	}
}


void default_obj_coll(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type, int cobj_type) {

	dwobject &obj(obj_groups[coll_id[cobj_type]].get_obj(index));

	if (type == CAMERA || type == SMILEY) {
		if (cobj_type == S_BALL) return;
		energy = get_coll_energy(zero_vector, obj.velocity, object_types[obj.type].mass);
		smiley_collision(((type == CAMERA) ? CAMERA_ID : obj_index), index, velocity, position, energy, cobj_type);
		if (!invalid_collision) obj.disable(); // return?
	}
	obj.elastic_collision(position, energy, type); // partially elastic collision
}


void landmine_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	if (type != CAMERA && type != SMILEY) return;
	if (type == CAMERA) obj_index = CAMERA_ID; else assert(obj_index < num_smileys);
	dwobject &obj(obj_groups[coll_id[LANDMINE]].get_obj(index));
	point const pos(get_sstate_pos(obj_index));

	if (sstates[obj_index].powerup == PU_FLIGHT && pos.z > object_types[SMILEY].radius + int_mesh_zval_pt_off(pos, 1, 0)) {
		invalid_collision = 1;
		return; // don't run into landmine when flying above it
	}
	if (obj_index == get_damage_source(LANDMINE, index, obj_index) && obj.time < (int)SMILEY_LM_ACT_TIME) {
		invalid_collision = 1;
		return; // camera/smiley ran into his own landmine
	}
	//smiley_collision(obj_index, index, velocity, obj.pos, energy, LANDMINE);
	blast_radius(obj.pos, LANDMINE, index, obj.source, 0);
	gen_smoke(obj.pos);
	obj.status = 0;
}


bool pushable_collision(int index, point const &position, float force, int type, int obj_type) {

	if (type == CAMERA || type == SMILEY) {
		dwobject &obj(obj_groups[coll_id[obj_type]].get_obj(index));

		if (obj.status != 1 && obj.status != 2) { // only if on ground or stopped
			// Note: an object resting on a destroyable static object will still have status 1 and will not be pushable
			if (obj.status == 4) obj.flags |= WAS_PUSHED;
			obj.elastic_collision(position, force, type); // add some extra energy so that we can push the skull
			return 1;
		}
	}
	return 0;
}

// something runs into a dodgeball
void dodgeball_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	if (type == CAMERA || type == SMILEY) {
		if (pushable_collision(index, position, 100.0, type, BALL)) {
			if (game_mode != 2) return; // doesn't seem to help
			energy = 0.0;
		}
		smiley_collision(((type == CAMERA) ? CAMERA_ID : obj_index), index, velocity, position, energy, BALL);
	}
	else {
		default_obj_coll(index, obj_index, velocity, position, energy, type, BALL);
	}
}


void skull_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	pushable_collision(index, position, 2000.0, type, SKULL);
}


void health_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	default_obj_coll(index, obj_index, velocity, position, energy, type, HEALTH);
}

void shield_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	default_obj_coll(index, obj_index, velocity, position, energy, type, SHIELD);
}


void powerup_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	default_obj_coll(index, obj_index, velocity, position, energy, type, POWERUP);
}


void weapon_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	default_obj_coll(index, obj_index, velocity, position, energy, type, WEAPON);
}


void ammo_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	default_obj_coll(index, obj_index, velocity, position, energy, type, AMMO);
}


void pack_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	default_obj_coll(index, obj_index, velocity, position, energy, type, WA_PACK);
}


void sball_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	default_obj_coll(index, obj_index, velocity, position, energy, type, S_BALL);
	pushable_collision(index, position, 20.0, type, S_BALL);
}


void rock_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {

	if (type != SEEK_D && type != ROCKET && type != IMPACT) return;
	float num(rand_uniform(0.0, 6.0));
	if (index == 0)          num *= 2.0; // large rock
	if      (type == SEEK_D) num *= 2.0;
	else if (type == ROCKET) num *= 1.5;
	int const shooter(get_damage_source(type, obj_index, index));
	float const p[7] = {2.5, 5.0, 4.0, 0.2, 1.0, 0.5, 0.5};
	gen_rubble(ROCK, int(num), position, shooter, p);
}


// ***********************************
// WEAPON/EFFECTS CODE
// ***********************************


void gen_rubble(int type, int num, point const &pos, int shooter, float const p[7]) {

	obj_group &objg(obj_groups[coll_id[type]]);

	for (int o = 0; o < num; ++o) {
		int const i(objg.choose_object());
		objg.create_object_at(i, pos);
		dwobject &obj(objg.get_obj(i));
		obj.init_dir      = signed_rand_vector_norm();
		obj.velocity      = gen_rand_vector(rand_uniform(p[0], p[1]), p[2], PI_TWO);
		obj.orientation.x = p[3] + p[4]*rand_float(); // size
		obj.orientation.y = p[5] + p[6]*rand_float(); // color
		obj.time          = int(0.5*rand_float()*(float)object_types[type].lifetime);
		obj.source        = shooter;
	}
}


void create_ground_rubble(point pos, int shooter, float hv, float close, int calc_hv) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	bool const outside(point_outside_mesh(xpos, ypos));
	if ((is_in_ice(xpos, ypos) && is_underwater(pos)) || is_mesh_disabled(xpos, ypos)) return;
	if (!outside && wminside[ypos][xpos] && (water_matrix[ypos][xpos] - mesh_height[ypos][xpos]) > 0.25*MAX_SPLASH_DEPTH) return;
	
	if (calc_hv) {
		hv = (outside ? 0.5 : get_rel_height(mesh_height[ypos][xpos], zmin, zmax));
	}
	if (island && hv < 0.25) { // sand
		int const num(int(close*(4.5*(0.2 - hv) + 0.1)*(50 + rand()%100)));
		float const params[7] = {4.0, 6.0, 4.0, 0.5, 0.5, 0.8, 0.2};
		gen_rubble(SAND, num, pos, shooter, params);
	}
	if (hv < 0.57 && (!island || hv > 0.2)) { // dirt
		int const num(int(close*(1.0 - 3.0*fabs(hv - 0.3))*(75 + rand()%150)));
		float const params[7] = {3.5, 5.5, 5.0, 0.3, 0.7, 0.6, 0.4};
		gen_rubble(DIRT, num, pos, shooter, params);
	}
	if (hv > 0.5 || (island && hv > 0.35)) { // rocks
		int const num(int(close*(3.0*min(0.3, (hv - 0.5 + 0.15*(island == 1))) + 0.1)*(30 + rand()%60)));
		float const params[7] = {3.0, 6.0, 4.5, 0.2, 1.0, 0.5, 0.5};
		gen_rubble(ROCK, num, pos, shooter, params);
	}
}


void dwobject::update_vel_from_damage(vector3d const &dv) {

	velocity += dv*min(1.0, 2.0*object_types[type].terminal_vel/dv.mag()); // 2x terminal velocity
}


void dwobject::damage_object(float damage, point const &dpos, point const &shoot_pos, int weapon) {

	point const pos0((weapon == W_BLADE) ? shoot_pos : dpos);
	float const damage2(0.1f*damage/max(1.0E-6f, (p2p_dist(pos, pos0)*sqrt(object_types[type].mass)))); // careful of divide by zero
	health -= HEALTH_PER_DAMAGE*damage;
	flags  &= ~STATIC_COBJ_COLL;

	if (health >= 0 && weapon != W_LASER) { // still alive, send it flying
		status = 1;
		if (type != LANDMINE) flags &= ~ALL_COLL_STOPPED;
		update_vel_from_damage((pos - pos0)*((weapon == W_BLADE) ? -1.6*damage2 : damage2)); // W_BLADE : pull object towards you
	}
}


void blast_radius(point const &pos, int type, int obj_index, int shooter, int chain_level) {

	point const temp1(camera_origin);
	vector3d const temp2(cview_dir);
	if (!game_mode) return;
	int const wtype(obj_weapons[type]);
	if (wtype < 0)  return;
	assert(type >= 0);
	if (BLAST_CHAIN_DELAY > 0) gen_smoke(pos);
	assert(wtype <= NUM_WEAPONS);
	float damage(weapons[wtype].blast_damage), size(weapons[wtype].blast_radius);
	dwobject const &obj(obj_groups[coll_id[type]].get_obj(obj_index));

	if (type == PLASMA) {
		float const psize(obj.init_dir.x);
		assert(psize >= 0.0);
		damage *= psize*psize;
		size   *= psize*psize;
	}
	if (following) {
		camera_origin = orig_camera;
		cview_dir     = orig_cdir;
	}
	if (type != IMPACT && type < NUM_TOT_OBJS && coll_id[type] > 0) shooter = obj.source;
	damage *= sstates[shooter].get_damage_scale();
	bool const cview((obj.flags & CAMERA_VIEW) != 0);
	create_explosion(pos, shooter, chain_level, damage, size, type, cview);

	if (following) {
		camera_origin = temp1;
		cview_dir     = temp2;
	}
}


// returns true if there are no objects blocking the explosion
bool check_explosion_damage(point const &p1, point const &p2, int cobj) {

	int cindex;
	if (line_intersect_mesh(p1, p2)) return 0;
	if (!check_coll_line(p1, p2, cindex, cobj, 1, 0)) return 1;
	assert((unsigned)cindex < coll_objects.size());
	return (coll_objects[cindex].destroy >= SHATTERABLE); // blocked by a non destroyable static object
}


void exp_damage_groups(point const &pos, int shooter, int chain_level, float damage, float size, int type, bool cview) { // slow

	float dist(distance_to_camera(pos)), ssq(size*size);

	if (dist <= size && (type != IMPACT || shooter != CAMERA_ID) && (type != SEEK_D || !cview)) {
		if (check_explosion_damage(pos, get_camera_pos(), camera_coll_id)) {
			br_source = type;
			camera_collision(type, shooter, zero_vector, pos, damage*(1.02 - dist/size), BLAST_RADIUS);
		}
	}
	for (int g = 0; g < num_groups; ++g) { // apply blast radius damage to objects
		obj_group &objg(obj_groups[g]);
		if (!objg.enabled) continue;
		int type2(objg.type);
		bool const large_obj(objg.large_radius());

		if (objg.flags & PRECIPITATION) { // rain isn't affected
			if (temperature <= RAIN_MIN_TEMP) type2 = PRECIP; else continue;
		}
		if (type2 == CAMERA || type2 == PLASMA || type2 == BLOOD || type2 == CHARRED) continue;
		assert(object_types[type2].mass > 0.0);
		bool const can_move(object_types[type2].friction_factor < 3.0*STICK_THRESHOLD);
		float const dscale(0.1/sqrt(object_types[type2].mass));
		unsigned const nobj(objg.end_id);
		assert(nobj <= objg.max_objects());

		for (unsigned i = 0; i < nobj; ++i) {
			if (!objg.obj_within_dist(i, pos, size)) continue; // size+radius?
			dwobject &obj(objg.get_obj(i));
			if (large_obj && !check_explosion_damage(pos, obj.pos, obj.coll_id)) continue; // blocked by an object
			float const damage2(damage*(1.02 - p2p_dist(obj.pos, pos)/size));
			
			if (type2 == SMILEY && (type != IMPACT || shooter != (int)i)) {
				br_source = type;
				smiley_collision(i, shooter, zero_vector, pos, damage2, BLAST_RADIUS);
			}
			else {
				obj.health -= HEALTH_PER_DAMAGE*damage2;
				
				if (can_move) {
					obj.flags &= ~STATIC_COBJ_COLL;

					if (obj.health >= 0.0 || (BLAST_CHAIN_DELAY > 0 && chain_level > 0)) {
						if (temperature > W_FREEZE_POINT || !(obj.flags & IN_WATER)) { // not stuck in ice
							obj.status = 1;
							obj.flags &= ~ALL_COLL_STOPPED;
							obj.update_vel_from_damage((obj.pos - pos)*(damage2*dscale)); // similar to damage_object()
						}
					}
				}
				if (obj.health < 0.0) {
					if ((object_types[type2].flags & OBJ_EXPLODES) && (type2 != LANDMINE || !obj.lm_coll_invalid())) {
						if (BLAST_CHAIN_DELAY == 0) {
							obj.status = 0;
							blast_radius(obj.pos, type2, i, obj.source, chain_level+1);
							gen_smoke(obj.pos);
						}
						else {
							obj.health = object_types[type2].health; // ???
							obj.time   = max(obj.time, object_types[type2].lifetime-((int)BLAST_CHAIN_DELAY)*(chain_level+1));
						}
					}
					else {
						if (type2 == PLASMA) gen_fire(obj.pos, obj.init_dir.x, obj.source);
						obj.status = 0;
					}
				} // health test
			} // SMILEY test
		} // for i
	} // for g
}


void create_explosion(point const &pos, int shooter, int chain_level, float damage, float size, int type, bool cview) {

	assert(damage >= 0.0 && size >= 0.0);
	assert(type != SMILEY);
	if (!game_mode || damage < TOLERANCE || size < TOLERANCE) return;
	//RESET_TIME;
	
	if (type == GRENADE || type == CGRENADE) {
		add_blastr(pos, (pos - get_camera_pos()), 0.9*size, damage, int(1.5*BLAST_TIME), shooter, YELLOW, RED, ETYPE_STARB);
	}
	else {
		int const time(((type == BLAST_RADIUS) ? 2 : 1)*BLAST_TIME);
		add_blastr(pos, signed_rand_vector_norm(), 0.7*size, damage, int(1.5*time), shooter, WHITE, WHITE, ETYPE_ANIM_FIRE);
		//add_blastr(pos, signed_rand_vector_norm(), 0.7*size, damage, time, shooter, YELLOW, RED, ETYPE_FIRE);
	}
	//exp_cobjs.push_back(add_coll_sphere(pos, size, cobj_params(0.0, WHITE, 0, 1, explosion_coll, exp_cobjs.size()))); // cobj for next frame
	exp_damage_groups(pos, shooter, chain_level, damage, size, type, cview);

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	float depth(0.0);
	point pos_zr(pos);
	pos_zr.z -= 0.5*size;
	bool const underwater(is_underwater(pos_zr, 0, &depth));

	if (underwater) {
		depth = min(0.25f, max(depth, 0.01f));
		assert(damage >= 0.0);
		add_splash(xpos, ypos, 0.002*damage/depth, (0.4 + 2.0*depth)*size);
	}
	if (damage > 500.0) { // everything except for plasma
		gen_sound((underwater? SOUND_SPLASH2 : SOUND_EXPLODE), pos);
	}
	if (type == GRENADE) { // shrapnel fragments
		unsigned const num(weapons[W_GRENADE].nfragments + rand()%(weapons[W_GRENADE].nfragments/4));
		obj_group &objg(obj_groups[coll_id[SHRAPNEL]]);

		for (unsigned o = 0; o < num; ++o) {
			int const i(objg.choose_object());
			objg.create_object_at(i, pos);
			dwobject &obj(objg.get_obj(i));
			obj.velocity    = gen_rand_vector(rand_uniform(6.0, 12.0), 2.5, 0.75*PI);
			obj.orientation = signed_rand_vector_norm();
			obj.angle       = 360.0*rand_float();
			obj.source      = shooter;
			obj.direction   = W_GRENADE;
		}
	}
	else if (type == CGRENADE) { // grenades
		unsigned const num(weapons[W_CGRENADE].nfragments);
		obj_group &objg(obj_groups[coll_id[GRENADE]]);

		for (unsigned o = 0; o < num; ++o) { // less efficient
			int const i(objg.choose_object());
			objg.create_object_at(i, pos);
			dwobject &obj(objg.get_obj(i));
			obj.time     = o;
			obj.velocity = gen_rand_vector(rand_uniform(8.0, 10.0), 2.0, PI_TWO);
			obj.init_dir = signed_rand_vector_norm();
			obj.source   = shooter;
		}
	}
	if (size > 0.2) gen_particles(pos, (rand() % int(25*size)));
	
	// large damage - throws up dirt and makes craters (later destroys trees)
	if ((type == IMPACT || damage > 1000.0) && is_over_mesh(pos) && !point_outside_mesh(xpos, ypos)) {
		float const zval(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1));
		float const damage2(5.0E-6*Z_SCENE_SIZE*crater_size*damage*(256.0/(float)XY_SUM_SIZE)), crater_dist(0.36*size);

		if (fabs(zval - pos.z) < crater_dist) { // on/close to ground
			int const crater(damage >= 1000.0 && point_interior_to_mesh(xpos, ypos) && type != PLASMA);
			float hv(0.5), close(1.0);

			if (crater) {
				close = 1.1 - fabs(zval - pos.z)/crater_dist;
				size *= close*sqrt(XY_SCENE_SIZE);
				hv    = add_crater_to_landscape_texture(pos.x, pos.y, size);
				update_mesh_height(xpos, ypos, int(crater_dist/HALF_DXY), damage2, 0.0, 0);
				modify_grass_at(pos, crater_dist, 0, 0, 0, 1); // update grass height
			}
			if ((h_collision_matrix[ypos][xpos] - mesh_height[ypos][xpos]) < SMALL_NUMBER) {
				create_ground_rubble(pos, shooter, hv, close, !crater);
			}
		}
	}
	if (damage > 100.0) destroy_coll_objs(pos, damage, shooter, (type == BLAST_RADIUS));
	//PRINT_TIME("Blast Radius");
}


void do_impact_damage(point const &fpos, vector3d const &dir, vector3d const &velocity, point const &shoot_pos, float radius, int shooter, int weapon, float mag) {

	float damage(mag*weapons[weapon].blast_damage);

	if (weapon == W_BLADE) {
		int const ammo(UNLIMITED_WEAPONS ? 1 : sstates[shooter].p_ammo[weapon]);
		damage *= sqrt((double)ammo + 1.0)/SQRT2;
	}
	point pos(fpos + dir*(1.25*radius));
	float const coll_radius(0.75*radius);

	for (int g = 0; g < num_groups; ++g) {
		obj_group &objg(obj_groups[g]);
		if (!objg.enabled || (!objg.large_radius() && (objg.type != FRAGMENT || weapon != W_BLADE))) continue;
		int const type(objg.type);
		float const robj(object_types[type].radius), rad(coll_radius + robj);
		
		for (unsigned i = 0; i < objg.end_id; ++i) {
			if (type == SMILEY && (int)i == shooter) continue; // this is the shooter
			if (!objg.obj_within_dist(i, pos, rad))  continue;

			if (type == SMILEY) {
				++sstates[shooter].cb_hurt;
				smiley_collision(i, shooter, velocity, pos, damage, IMPACT);
				gen_sound(SOUND_SQUISH2, pos, 0.5);
			}
			else if (type == FRAGMENT) {
				objg.get_obj(i).status = 0; // just destroy the fragment
			}
			else if (type != CAMERA) {
				objg.get_obj(i).damage_object(damage, pos, shoot_pos, weapon);
			}
		}
	}
	if (shooter >= 0) {
		point const camera(get_camera_pos());

		if (dist_less_than(pos, camera, (coll_radius + CAMERA_RADIUS))) {
			camera_collision(weapon, shooter, velocity, pos, damage, IMPACT);
			gen_sound(SOUND_SQUISH2, pos, 1.0);
		}
	}
	damage *= sstates[shooter].get_damage_scale();

	if (weapon != W_BLADE || rand()%20 == 0) {
		for (unsigned i = 0; i < rock_shapes.size(); ++i) {
			if (rock_shapes[i].do_impact_damage(pos, radius)) rock_collision(0, -1, zero_vector, pos, damage, IMPACT);
		}
	}
	if (weapon == W_BLADE && sstates[shooter].fire_frame > 0 && (rand()&7) == 0) {
		destroy_coll_objs(pos, 18.0*damage, shooter, 0);
	}
	if (weapon == W_BBBAT && sstates[shooter].fire_frame > 0) {
		destroy_coll_objs(pos, 0.5*damage, shooter, 0);
	}
	int const xpos(get_xpos(fpos.x)), ypos(get_ypos(fpos.y));

	if (has_water(xpos, ypos) && water_matrix[ypos][xpos] > (fpos.z - radius) &&
		water_matrix[ypos][xpos] < (fpos.z + radius + MAX_SPLASH_DEPTH))
	{
		add_splash(xpos, ypos, damage, 0.6*radius);
	}
}


// fire and gas damage
void do_area_effect_damage(point &pos, float effect_radius, float damage, int index, int source, int type) {

	vector3d velocity(zero_vector);
	float const radius(object_types[SMILEY].radius + effect_radius), radius_sq(radius*radius);
	
	if (p2p_dist_sq(pos, get_camera_pos()) < radius_sq) { // what if camera/player is a different size or height from smiley?
		camera_collision(CAMERA_ID, ((source == NO_SOURCE) ? CAMERA_ID : source), velocity, pos, damage, type);
	}
	obj_group const &objg(obj_groups[coll_id[SMILEY]]);
	
	if (objg.enabled) {		
		for (unsigned i = 0; i < objg.end_id; ++i) {
			if (!objg.get_obj(i).disabled() && p2p_dist_sq(pos, objg.get_obj(i).pos) < radius_sq) {
				// test for objects blocking the damage effects?
				smiley_collision(i, ((source == NO_SOURCE) ? i : source), velocity, pos, damage, type);
			}
		}
	}
}


void switch_player_weapon(int val) {

	if (sstates != NULL && game_mode) sstates[CAMERA_ID].switch_weapon(val, 1);
}


void player_state::switch_weapon(int val, int verbose) {

	if (game_mode == 2) {
		weapon = ((UNLIMITED_WEAPONS || p_ammo[W_BALL] > 0) ? W_BALL : W_UNARMED);
		return;
	}
	do {
		weapon = (weapon+NUM_WEAPONS+val)%NUM_WEAPONS;
	} while (!UNLIMITED_WEAPONS && no_weap_or_ammo());

	if (verbose) print_weapon(weapon);
	wmode      = 0; // maybe don't reset this?
	fire_frame = 0;
	cb_hurt    = 0;
}


void player_state::gamemode_fire_weapon() { // camera/player fire

	static int fire_frame(0);
	if (frame_counter == fire_frame) return; // to prevent two fires in the same frame
	fire_frame = frame_counter;

	if (!game_mode) { // flashlight/candlelight mode only
		if (wmode & 1) add_camera_candlelight(); else add_camera_flashlight();
		return;
	}
	if (!camera_reset) return;
	point const camera(get_camera_pos());
	if (temperature <= W_FREEZE_POINT && is_underwater(camera)) return; // under ice
	
	if (following) {
		following    = 0;
		camera_reset = 1;
		return;
	}
	if (!UNLIMITED_WEAPONS && weapon != W_UNARMED && no_weap_or_ammo()) {
		if (weapon != W_ROCKET && weapon != W_SEEK_D && weapon != W_PLASMA && weapon != W_GRENADE) { // this test is questionable
			switch_weapon(1, 1);
			if (weapon == W_BBBAT)   switch_weapon( 1, 1);
			if (weapon == W_UNARMED) switch_weapon(-1, 1);
			if (game_mode == 2 && weapon == W_BBBAT) switch_weapon(1, 1);
			return; // no weapon/out of ammo
		}
	}
	else {
		verify_wmode();
		bool const fmode2(wmode & 1);
		int const psize((int)ceil(plasma_size));
		int chosen;
		int const status(fire_projectile(camera, cview_dir, CAMERA_ID, chosen));
		fired = 1;

		if (status == 2 && weapon == W_SEEK_D && fmode2) { // follow the seek and destroy
			if (chosen >= 0) obj_groups[coll_id[SEEK_D]].get_obj(chosen).flags |= CAMERA_VIEW;
			orig_camera  = camera_origin; // camera is actually at the SEEK_D location, not quite right, but interesting
			orig_cdir    = cview_dir;
			camera_reset = 0;
			following    = 1;
		}
		int &pammo(p_ammo[weapon]);

		if (status != 0 && !UNLIMITED_WEAPONS && !no_weap() && pammo > 0) {
			if (weapon == W_PLASMA && psize > 1) {
				pammo = max(0, pammo-psize); // large plasma burst
			}
			else if (weapon == W_GRENADE && fmode2 && (pammo >= int(weapons[W_CGRENADE].def_ammo) || UNLIMITED_WEAPONS)) {
				pammo -= weapons[W_CGRENADE].def_ammo; // cluster grenade
			}
			else if (weapon != W_BLADE) {
				--pammo;
			}
		}
	}
	//if (frame_counter == fire_frame) return;
	int const dtime(int(get_fspeed_scale()*(frame_counter - timer)*fticks));

	if ((game_mode == 2 && weapon == W_BBBAT) || (!UNLIMITED_WEAPONS && no_ammo() &&
		(weapon == W_LASER || dtime >= int(weapons[weapon].fire_delay))))
	{
		switch_weapon(1, 1);
		if (weapon == W_UNARMED) switch_weapon(1, 1);
		if (weapon == W_BBBAT)   switch_weapon(1, 1);
		if (weapon == W_UNARMED) switch_weapon(1, 1);
	}
}


void add_laser_beam(laser_beam const &laser) {

	lasers.push_back(laser);
	if (LASER_PATH_LIGHT) add_line_light(laser.pts[0], laser.pts[1], RED, 0.4, min(1.0f, sqrt(laser.intensity)));
}


void create_shell_casing(point const &fpos, vector3d const &dir, int shooter, float radius, unsigned char type) {

	obj_group &objg(obj_groups[coll_id[SHELLC]]);
	int const i(objg.choose_object());
	objg.create_object_at(i, fpos);
	dwobject &obj(objg.get_obj(i));
	obj.pos.z      -= 0.8*radius;
	obj.pos        += dir*(1.2*radius);
	obj.velocity    = gen_rand_vector(rand_uniform(0.0, 0.25*(type + 1)), 3.0, PI) + vector3d(dir.y*0.7, -dir.x*0.7, 5.0);
	obj.orientation = signed_rand_vector_norm();
	obj.angle       = rand_uniform(0.0, TWO_PI);
	//obj.init_dir    = dir;
	obj.init_dir    = vector3d(PI*signed_rand_float(), 0.0, 0.0); // angle
	obj.time        = -1;
	obj.source      = shooter;
	obj.direction   = type; // 0 = M16, 1 = shotgun
}


void create_shrapnel(point const &pos, vector3d const &dir, float firing_error, unsigned nshots, int shooter, int weapon) {

	unsigned const num(2*nshots + 2);
	unsigned nobj(0);
	obj_group &objg(obj_groups[coll_id[SHRAPNEL]]);

	for (unsigned o = 0; o < num; ++o) {
		int const i(objg.choose_object());
		objg.create_object_at(i, pos);
		dwobject &obj(objg.get_obj(i));
		obj.angle       = 1.0;
		obj.velocity    = dir;
		vadd_rand(obj.velocity, firing_error);
		obj.velocity   *= rand_uniform(65.0, 80.0);
		obj.orientation = signed_rand_vector_norm();
		obj.angle       = 360.0*rand_float();
		obj.source      = shooter;
		obj.direction   = weapon;
	}
}


float weapon_t::get_fire_vel() const {

	return (v_add + ball_velocity*v_mult);
}


int player_state::fire_projectile(point fpos, vector3d dir, int shooter, int &chosen_obj) {

	chosen_obj = -1;
	float damage_scale(1.0), range;
	vector3d velocity(zero_vector);
	assert(UNLIMITED_WEAPONS || !no_weap_or_ammo());
	int const weapon_id((weapon == W_GRENADE && (wmode&1)) ? W_CGRENADE : weapon);
	if (weapon_id == W_M16 && (wmode&1) == 1) ++rot_counter;
	int const dtime(int(get_fspeed_scale()*(frame_counter - timer)*fticks));
	bool const rapid_fire(weapon_id == W_ROCKET && (wmode&1));
	weapon_t const &w(weapons[weapon_id]);
	int fire_delay((int)w.fire_delay);
	if (UNLIMITED_WEAPONS && shooter != CAMERA_ID && weapon_id == W_LANDMINE) fire_delay *= 2; // avoid too many landmines
	if (rapid_fire) fire_delay /= 3;
	float const radius(get_sstate_radius(shooter));

	if (weapon_id == W_LASER) { // always fires
		assert(fire_delay > 0);
		damage_scale = fticks/fire_delay;
	}
	else if (dtime < fire_delay) {
		return 0;
	}
	timer = frame_counter;
	bool const underwater(is_underwater(fpos));
	float firing_error(w.firing_error*(underwater ? UWATER_FERR_MUL : 1.0));
	if (rapid_fire) firing_error *= 20.0;
	dir.normalize();
	point pos(fpos + dir*(0.1*radius));
	fire_frame = max(1, fire_delay);
	float const damage(damage_scale*w.blast_damage), vel(w.get_fire_vel());

	switch (weapon_id) {
		case W_M16:     add_dynamic_light(1.0, fpos, YELLOW); break;
		case W_SHOTGUN: add_dynamic_light(1.3, fpos, YELLOW); break;
		case W_LASER:   add_dynamic_light(0.6, fpos, RED);    break;
		case W_PLASMA:  add_dynamic_light(0.9, fpos, ORANGE); break;
	}
	switch (weapon_id) {
	case W_M16: // line of sight damage
		gen_sound(SOUND_GUNSHOT, fpos, 0.5);

		if ((wmode&1) != 1) { // not firing shrapnel
			if (dtime > 10) firing_error *= 0.1;
			if (underwater) firing_error += UWATER_FERR_ADD;
			projectile_test(fpos, dir, firing_error, damage, shooter, range);
			create_shell_casing(fpos, dir, shooter, radius, 0);
			return 1;
		} // fallthrough to shotgun case
	case W_SHOTGUN:
		if ((wmode&1) == 1) { // shrapnel cannon/chaingun
			create_shrapnel(pos, dir, firing_error, w.nshots, shooter, weapon_id);
		}
		else { // normal 12-gauge/M16
			if (underwater) firing_error += UWATER_FERR_ADD;

			for (int i = 0; i < int(w.nshots); ++i) { // can be slow if trees are involved
				projectile_test(fpos, dir, firing_error, damage, shooter, range);
			}
		}
		if (weapon_id == W_SHOTGUN) {
			for (unsigned i = 0; i < 2; ++i) {
				create_shell_casing(fpos, dir, shooter, radius, 1);
			}
			gen_sound(SOUND_SHOTGUN, fpos);
		}
		return 1;

	case W_BBBAT: // baseball bat
		do_impact_damage(fpos, dir, velocity, fpos, radius, shooter, weapon_id, 1.0);
		//cobj.register_coll() ???
		return 1;

	case W_PLASMA:
		if ((wmode&1) == 1) {
			plasma_loaded = !plasma_loaded;

			if (plasma_loaded) {
				plasma_size = 1.0;
				return 0; // loaded but not fired
			}
		}
		else plasma_loaded = 0;
		gen_sound(SOUND_FIREBALL, fpos);
		break;

	case W_LASER: // line of sight damage
		{
			projectile_test(fpos, dir, firing_error, damage, shooter, range);
			laser_beam const laser((range >= 0.9*FAR_CLIP), shooter, (fpos + dir*radius), (fpos + dir*range), RED);
			add_laser_beam(laser); // might not need to actually add laser itself for camera/player
		}
		break;

	case W_SEEK_D:
		fpos += dir*radius; // fire from in front of shooter, but then there is no recoil
	case W_ROCKET: // fallthrough from above
		gen_sound(SOUND_ROCKET, fpos);
		break;

	case W_GASSER:
		{
			float const r(radius + w.blast_radius);
			point start_pos(fpos + dir*(0.5*r));
			
			if (is_underwater(start_pos)) {
				colorRGBA color(DK_GREEN);
				color.alpha = 0.75;

				for (unsigned n = 0; n < 4; ++n) {
					start_pos += signed_rand_vector(0.05*radius);
					gen_bubble(start_pos, 0.0, color);
				}
			}
			else {
				start_pos   += dir*(0.5*r);
				start_pos.z -= 0.25*r;
				vector3d dir2(dir);
				if (firing_error != 0.0) vadd_rand(dir2, firing_error);
				bool const is_fire(wmode & 1);
				vector3d const gas_vel(dir2*vel + vector3d(0.0, 0.0, 0.2));
				colorRGBA const color(is_fire ? colorRGBA(1.0, 0.75, 0.0) : DK_GREEN);
				int const smoke_type (is_fire ? FIRE : GASSED);
				float const density(0.5*rand_uniform(0.5, 1.0));
				float const darkness(0.6*rand_uniform(0.7, 1.0));
				float const radius(w.blast_radius*rand_uniform(0.8, 1.2));
				gen_arb_smoke(start_pos, color, gas_vel, radius, density, darkness, w.blast_damage, shooter, smoke_type, 0);
				if (is_fire) gen_sound(SOUND_FIREBALL, start_pos);
			}
		}
		return 1;

	case W_BLADE:
		gen_sound(SOUND_DRILL, fpos, 0.3);
		break;
	}
	int type(w.obj_id);
	if (type < 0) return 3;
	int const cid(coll_id[type]);
	assert(cid >= 0 && cid < NUM_TOT_OBJS);
	obj_group &objg(obj_groups[cid]);
	assert(objg.max_objs > 0);
	float const rdist(0.75 + ((weapon_id == W_PLASMA) ? 0.5*(plasma_size - 1.0) : 0.0)); // change?
	float const radius2(radius + object_types[type].radius);
	assert(w.nshots <= objg.max_objs);
	bool const dodgeball(game_mode == 2 && weapon_id == W_BALL && !UNLIMITED_WEAPONS);
	if (dodgeball) assert(w.nshots <= balls.size());

	for (unsigned shot = 0; shot < w.nshots; ++shot) {
		int const chosen(dodgeball ? balls.back() : objg.choose_object());
		if (dodgeball) balls.pop_back();
		chosen_obj = chosen;
		assert(chosen >= 0); // make sure there is an object available
		vector3d dir2(dir);

		if (firing_error != 0.0) {
			vadd_rand(dir2, firing_error);
			dir2.normalize();
		}
		velocity = dir2*vel;
		objg.create_object_at(chosen, fpos);
		dwobject &obj(objg.get_obj(chosen));
		obj.pos      += dir2*(rdist*radius2);
		obj.velocity  = velocity;
		obj.init_dir  = dir2;
		obj.init_dir.negate();
		obj.time      = -1;
		obj.source    = shooter;
		obj.direction = rapid_fire;
		if (rapid_fire) obj.time = int((1.0 - rand_uniform(0.1, 0.9)*rand_uniform(0.1, 0.9))*object_types[type].lifetime);

		switch (weapon_id) {
		case W_PLASMA:
			obj.init_dir.x  = float(pow(double(plasma_size), 0.75)); // psize
			obj.pos.z      += 0.2*radius2;
			plasma_size     = 1.0;
			break;
		case W_BALL:
			obj.pos.z += 0.2*radius2;
			break;
		case W_STAR5:
			obj.init_dir += gen_rand_vector(0.1, 1.0, PI);
			obj.init_dir.normalize();
			break;
		}
	}
	if ((wmode&1) == 0 && weapon_id == W_PLASMA) plasma_size = 1.0;
	return 2;
}


int get_range_to_mesh(point const &pos, vector3d const &vcf, point &coll_pos, int &xpos, int &ypos) {

	bool ice_int(0);
	int ixpos(0), iypos(0);
	float zval;
	vector3d const vca(pos + vcf*FAR_CLIP);
	point ice_coll_pos(vca);

	// compute range to ice surface (currently only at the default water level)
	if (temperature <= W_FREEZE_POINT && !island) { // firing into ice
		ixpos = get_xpos(pos.x);
		iypos = get_ypos(pos.y);

		if (has_water(ixpos, iypos) && pos.z < water_matrix[iypos][ixpos]) {
			coll_pos = pos;
			xpos     = ixpos;
			ypos     = iypos;
			return 2; // can't get any range if under ice
		}
		if (vcf.z < -1.0E-6) { // firing down
			float const t(-(pos.z - water_plane_z)/vcf.z), xval(pos.x + t*vcf.x), yval(pos.y + t*vcf.y);
			ixpos = get_xpos(xval);
			iypos = get_ypos(yval);

			if (mesh_is_underwater(ixpos, iypos)) {
				ice_coll_pos.assign(xval, yval, water_matrix[iypos][ixpos]);
				ice_int = 1;
			}
		}
	}
	if (fabs(vcf.z) > TOLERANCE && line_intersect_mesh(pos, vca, xpos, ypos, zval)) { // skip if dir has no z component
		assert(!point_outside_mesh(xpos, ypos));
		//coll_pos.assign(get_xval(xpos), get_yval(ypos), (mesh_height[ypos][xpos] + SMALL_NUMBER));
		coll_pos = (pos + (vca - pos)*(zval - pos.z)/(vca.z - pos.z));
		if (!ice_int || p2p_dist_sq(pos, ice_coll_pos) > p2p_dist_sq(pos, coll_pos)) return 1; // collides with mesh
	}
	if (ice_int) { // collides with ice before mesh
		coll_pos = ice_coll_pos;
		xpos     = ixpos;
		ypos     = iypos;
		return 2;
	}
	return 0;
}


void add_laser_beam_segment(point const &start_pos, point coll_pos, vector3d const &vref, int coll, bool distant, float intensity) {

	if (!coll || distant) {
		coll_pos = start_pos + vref*FAR_CLIP;
	}
	else if (start_pos == coll_pos) {
		return;
	}
	add_laser_beam(laser_beam(distant, NO_SOURCE, start_pos, coll_pos, RED, intensity));
}


point projectile_test(point const &pos, vector3d const &vcf_, float firing_error, float damage,
					  int shooter, float &range, float intensity, int ignore_cobj)
{
	assert(!is_nan(damage));
	assert(intensity <= 1.0 && LASER_REFL_ATTEN < 1.0);
	int closest(-1), closest_t(0), coll(0), xpos(0), ypos(0), cindex(-1);
	point coll_pos, res_pos(pos);
	vector3d vcf(vcf_), vca(pos), coll_norm(plus_z);
	float const vcf_mag(vcf.mag()), MAX_RANGE(min((float)FAR_CLIP, 2.0f*(X_SCENE_SIZE + Y_SCENE_SIZE + Z_SCENE_SIZE)));
	float specular(0.0), luminance(0.0), alpha(1.0), refract_ix(1.0);
	player_state const &sstate(sstates[shooter]);
	int const wtype(sstate.weapon), wmode(sstate.wmode);
	float const radius(get_sstate_radius(shooter));
	bool const is_laser(wtype == W_LASER);
	range = MAX_RANGE;
	vcf  /= vcf_mag;

	if (firing_error != 0.0) {
		vcf += signed_rand_vector_spherical(firing_error, 0);
		vcf.normalize();
	}
	vector3d const vcf0(vcf*vcf_mag);
	vca += vcf*FAR_CLIP;
	int const intersect(get_range_to_mesh(pos, vcf, coll_pos, xpos, ypos));

	if (intersect) {
		range = p2p_dist(pos, coll_pos);
		if (intersect == 1) coll_pos = pos + vcf*(range - 0.01); // not an ice intersection - simple and inexact, but seems OK
		coll_pos.z += SMALL_NUMBER;
	}

	// search for collisions with static objects (like trees)
	int const laser_m2(SELF_LASER_DAMAGE && is_laser && (wmode&1) && intensity < 1.0);
	int const proj_type(is_laser ? BEAM : PROJECTILE);
	range = get_projectile_range(pos, vcf, 0.01*radius, range, coll_pos, coll_norm, coll, cindex, shooter, !is_laser, ignore_cobj);
	if (cindex >= 0) assert(unsigned(cindex) < coll_objects.size());
	point cpos(pos), p_int(pos), lsip(pos);

	if (is_laser && cindex >= 0) {
		cobj_params const &cp(coll_objects[cindex].cp);
		get_lum_alpha(cp.color, cp.tid, luminance, alpha);
		specular    = cp.specular;
		refract_ix  = cp.refract_ix;
	}
	for (int g = 0; g < num_groups; ++g) { // collisions with dynamic group objects - this can be slow
		obj_group const &objg(obj_groups[g]);
		int const type(objg.type);
		obj_type const &otype(object_types[type]);
		if (!objg.enabled || !objg.large_radius() || type == PLASMA) continue;
		
		for (unsigned i = 0; i < objg.end_id; ++i) {
			if (type == SMILEY && (int)i == shooter && !laser_m2) continue; // this is the shooter
			if (!objg.obj_within_dist(i, pos, (range + otype.radius + SMALL_NUMBER))) continue;
			point const &apos(objg.get_obj(i).pos);

			if (line_sphere_int(vcf0, lsip, apos, otype.radius, p_int, 1)) {
				closest   = i;
				closest_t = type;
				cpos      = p_int;
				range     = p2p_dist(pos, cpos);
				
				if (is_laser && type != SMILEY && type != CAMERA) {
					get_lum_alpha(otype.color, otype.tid, luminance, alpha);
					specular = ((otype.flags & SPECULAR) ? 1.0 : ((otype.flags & LOW_SPECULAR) ? 0.5 : 0.0));
				}
			}
		}
	}
	if (shooter >= 0 || laser_m2) { // check camera/player
		point const camera(get_camera_pos());
		float const dist(distance_to_camera(pos));

		if (dist < (range + CAMERA_RADIUS + SMALL_NUMBER) && line_sphere_int(vcf0, lsip, camera, CAMERA_RADIUS, p_int, 1)) {
			closest   = 0;
			closest_t = CAMERA;
			cpos      = p_int;
			range     = p2p_dist(pos, cpos);
			camera_collision(wtype, shooter, zero_vector, cpos, damage, proj_type);
		}
	}
	if (closest_t == SMILEY) {
		smiley_collision(closest, shooter, zero_vector, cpos, damage, proj_type);
		if (is_laser && (rand()%6) == 0) gen_smoke(coll_pos);
	}
	if ((intersect || coll) && dist_less_than(coll_pos, pos, (X_SCENE_SIZE + Y_SCENE_SIZE))) { // spark
		if (!is_underwater(coll_pos)) {
			colorRGBA const scolor(is_laser ? RED : ORANGE);
			float const ssize((is_laser ? ((wmode&1) ? 0.015 : 0.020)*intensity : 0.025)*((closest_t == CAMERA) ? 0.5 : 1.0));
			sparks.push_back(spark_t(coll_pos, scolor, ssize));
			point const light_pos(coll_pos - vcf*(0.1*ssize));
			add_dynamic_light(0.6*CLIP_TO_01(sqrt(intensity)), light_pos, scolor);

			if (coll && intensity >= 1.0 && (!is_laser || (alpha == 1.0 && ((rand()&1) == 0)))) {
				gen_particles(light_pos, 1, 0.5, 1); // sparks
			}
		}
		res_pos = coll_pos;
	}
	else if (closest >= 0) {
		res_pos = cpos;
	}
	if (coll && cindex >= 0 && closest < 0) { // hit cobjs (like tree leaves)
		coll_objects[cindex].register_coll(TICKS_PER_SECOND/(is_laser ? 4 : 2), (is_laser ? BEAM : PROJECTILE));

		if (!is_laser && coll_objects[cindex].can_be_scorched()) {
			bool const is_glass(coll_objects[cindex].cp.color.alpha <= 0.5);
			gen_scorch_mark(coll_pos, 0.005, coll_norm, cindex, 1.0, (is_glass ? 1.0 : 0.0));
		}
	}
	
	// process laser reflections
	bool const reflects((intersect == 2 && !coll && closest == -1) || ((wmode&1) && coll));
	
	if (is_laser && ((coll && alpha < 1.0) || reflects) &&
		closest_t != CAMERA && closest_t != SMILEY && intensity > 0.01)
	{
		float range0(MAX_RANGE), range1(range0), atten(LASER_REFL_ATTEN);
		point end_pos(coll_pos);
		vector3d vref(vcf);

		if (coll) { // hit coll obj
			float reflect(alpha);

			if (reflects) {
				if (alpha < 1.0 && refract_ix != 1.0) { // semi-transparent - fresnel reflection
					reflect = get_reflected_weight(get_fresnel_reflection(vcf*-1, coll_norm, 1.0, refract_ix), alpha);
				}
				else { // specular + diffuse reflections
					reflect = CLIP_TO_01(alpha*(specular + (1.0f - specular)*luminance)); // could use red component
				}
				if (reflect > 0.01) { // reflected light
					calc_reflection_angle(vcf, vref, coll_norm);
					end_pos = projectile_test(coll_pos, vref, 0.0, reflect*damage, shooter, range0, reflect*intensity);
				}
			}
			if (alpha < 1.0) {
				float refract(1.0 - reflect);

				if (refract > 0.01) { // refracted light (index of refraction changes angle?)
					if (coll && cindex >= 0) refract *= coll_objects[cindex].get_light_transmit(coll_pos, (coll_pos + vcf*FAR_CLIP));
					refract = min(LASER_REFL_ATTEN, refract);

					if (refract > 0.01) {
						point const cp(projectile_test(coll_pos, vcf, 0.0, refract*damage, shooter, range1, refract*intensity, cindex));
						add_laser_beam_segment(coll_pos, cp, vcf, coll, (range1 > 0.9*MAX_RANGE), refract*intensity);
					}
				}
			}
		}
		else { // hit ice - ice may not actually be in +z (if it has frozen ripples)
			calc_reflection_angle(vcf, vref, wat_vert_normals[ypos][xpos]); // don't have water surface normals
			end_pos = projectile_test(coll_pos, vref, 0.0, atten*damage, shooter, range0, atten*intensity);
			coll    = (end_pos != coll_pos);
		}
		if (reflects && atten > 0.0) {
			add_laser_beam_segment(coll_pos, end_pos, vref, coll, (range0 > 0.9*MAX_RANGE), atten*intensity);
		}
	}
	float const dscale((shooter >= CAMERA_ID) ? sstate.get_damage_scale() : 1.0);

	if (closest < 0) {
		if (intersect == 1 && !coll) {
			surface_damage[ypos][xpos] += dscale;
			if (is_laser) modify_grass_at(res_pos, 0.25*HALF_DXY, 0, 1, 0, 0); // burn
		}
		return res_pos;
	}
	if (closest_t != SMILEY && closest_t != CAMERA) {
		obj_groups[coll_id[closest_t]].get_obj(closest).damage_object(dscale*damage, res_pos, pos, wtype);
		return res_pos;
	}
	return cpos;
}


// inefficient and inaccurate - this should be made better using real line/coll_obj intersection
float get_projectile_range(point const &pos, vector3d vcf, float dist, float range, point &coll_pos, vector3d &coll_norm,
						   int &coll, int &cindex, int source, int check_splash, int ignore_cobj)
{
	vcf.normalize();
	float const splash_val((!DISABLE_WATER && check_splash && (temperature > W_FREEZE_POINT || island)) ? SPLASH_BASE_SZ*100.0 : 0.0);
	point const pos1(pos + vcf*dist), pos2(pos + vcf*range);
	coll = 0;

	if (check_coll_line_exact(pos1, pos2, coll_pos, coll_norm, cindex, splash_val, ignore_cobj)) {
		assert(cindex >= 0 && (unsigned)cindex < coll_objects.size());
		coll_obj &cobj(coll_objects[cindex]);
		if (cobj.cp.coll_func) cobj.cp.coll_func(cobj.cp.cf_index, 0, zero_vector, pos, 0.0, PROJC); // apply collision function
		coll  = 1;
		range = p2p_dist(coll_pos, pos);
	}
	if (splash_val > 0.0 && is_underwater(pos1)) gen_line_of_bubbles(pos1, (pos + vcf*range));
	return range;
}


void do_cblade_damage_and_update_pos(point &pos, int shooter) {

	player_state &sstate(sstates[shooter]);
	int fframe(sstate.fire_frame);
	int const delay(max(1u, weapons[sstate.weapon].fire_delay));
	float const cradius(object_types[SMILEY].radius);
	vector3d const dir(get_sstate_dir(shooter));
	point const shoot_pos(pos);
	pos.z -= 0.05;

	if (fframe > 0) { // carnage blade extension
		point coll_pos;
		vector3d coll_norm; // unused
		int coll, xpos, ypos, cindex;
		int const fdir(fframe > (delay/2)), ff(fdir ? (delay - fframe) : fframe); // fdir = forward
		float range(get_projectile_range(pos, dir, 1.1*cradius, 1.5*cradius+CBLADE_EXT, coll_pos, coll_norm, coll, cindex, shooter, 0));
		bool cobj_coll(coll && cindex >= 0);

		if (get_range_to_mesh(pos, dir, coll_pos, xpos, ypos)) {
			float const mesh_range(p2p_dist(pos, coll_pos)-0.1f);

			if (mesh_range < range) {
				range     = mesh_range;
				cobj_coll = 0;
			}
			if (CBLADE_EXT_PT*ff > (range - 0.8f*cradius)) modify_grass_at(coll_pos, 0.75*cradius, 0, 0, 1, 0); // cut grass
		}
		if (cobj_coll) {
			assert(unsigned(cindex) < coll_objects.size());
			coll_objects[cindex].register_coll(TICKS_PER_SECOND/2, IMPACT);
		}
		sstate.dpos = max(0.0f, min(CBLADE_EXT_PT*ff, (range - 0.8f*cradius)));
		pos += dir*sstate.dpos;
	}
	// always doing damage
	//draw_sphere_at(pos, 0.01, 16);
	do_impact_damage(pos, dir, zero_vector, shoot_pos, cradius, shooter, W_BLADE, (1.0 + 0.25*(sstate.dpos > 0)));
	pos.z += 0.05;

	// throw up debris
	if ((rand()%6) == 0) {
		point pos2(pos + dir*(1.25*cradius));
		int const xpos(get_xpos(pos2.x)), ypos(get_ypos(pos2.y));
		
		if (!point_outside_mesh(xpos, ypos) && (pos.z - mesh_height[ypos][xpos]) < 0.5*cradius) {
			surface_damage[ypos][xpos] += 1.0;

			if (rand()%5 == 0) {
				pos2.z = mesh_height[ypos][xpos] + 0.25*cradius;
				create_ground_rubble(pos2, shooter, 0.0, 0.01, 1);
			}
		}
	}
	if (fframe > 0) pos -= dir*(1.0*cradius);
}


// ***********************************
// DRAWING CODE
// ***********************************


void show_user_stats() {

	bool const is_smiley0(spectate && num_smileys > 0 && obj_groups[coll_id[SMILEY]].enabled);
	colorRGBA color;
	player_state &sstate(sstates[is_smiley0 ? 0 : CAMERA_ID]);
	int nkills(sstate.kills), ndeaths(sstate.deaths);
	float chealth(is_smiley0 ? obj_groups[coll_id[SMILEY]].get_obj(0).health : camera_health);
	static char text[MAX_CHARS];
	static int *team_kills = NULL, *team_deaths = NULL;

	if (teams > 1) {
		if (team_kills  == NULL) team_kills  = new int[teams];
		if (team_deaths == NULL) team_deaths = new int[teams];
		
		for (unsigned i = 0; i < (unsigned)teams; ++i) {
			team_kills[i]  = 0;
			team_deaths[i] = 0;
		}
		team_kills[teams-1]  = nkills;
		team_deaths[teams-1] = ndeaths;
	}
	if (camera_mode == 1 && camera_surf_collide) {
		int const ammo((UNLIMITED_WEAPONS && weapons[sstate.weapon].need_ammo) ? -666 : sstate.p_ammo[sstate.weapon]);
		RED.do_glColor();
		sprintf(text, "%s %d  %s %d  %s %d  Frags %d  Best %d  Total %d  Deaths %d",
			((chealth < 25.0) ? "HEALTH" : "Health"), int(chealth + 0.5),
			((sstate.shields < 25.0) ? "SHIELDS" : "Shields"), int(sstate.shields + 0.5),
			((sstate.no_ammo()) ? "AMMO" : "Ammo"), ammo,
			frags, max(best_frags, -ndeaths), tot_frags, ndeaths);
		draw_text(-0.014, -0.012, -0.022, text);

		if (sstate.powerup_time > 0 && sstate.powerup >= 0) {
			colorRGBA const pc(get_powerup_color(sstate.powerup));
			pc.do_glColor();
			sprintf(text, "%is %s", int(sstate.powerup_time/TICKS_PER_SECOND + 0.5), powerup_names[sstate.powerup].c_str());
			draw_text(-0.015, -0.012, -0.025, text);
		}
	}
	if (show_scores) {
		for (unsigned i = 0; i < (unsigned)num_smileys; ++i) {
			color = get_smiley_team_color(i);
			color.do_glColor();
			sprintf(text, "%s: kills = %i, deaths = %i, score = %i\n",
				sstates[i].name.c_str(), sstates[i].kills, sstates[i].deaths, sstates[i].kills-sstates[i].deaths);
			draw_text(-0.008, 0.01-0.0014*i, -0.02, text);
			nkills  += sstates[i].kills;
			ndeaths += sstates[i].deaths;

			if (teams > 1) {
				team_kills[i%teams]  += sstates[i].kills;
				team_deaths[i%teams] += sstates[i].deaths;
			}
		}
		if (teams > 1) {
			for (unsigned i = 0; i < (unsigned)teams; ++i) {
				color = get_smiley_team_color(i);
				color.do_glColor();
				sprintf(text, "Team %u: kills = %i, deaths = %i, score = %i\n",
					i, team_kills[i], team_deaths[i], team_kills[i]-team_deaths[i]);
				draw_text(-0.008, 0.01-0.0014*(i+num_smileys)-0.0008, -0.02, text);
			}
		}
		WHITE.do_glColor();
		sprintf(text, "Total: kills = %i, deaths = %i, score = %i\n", nkills, ndeaths, nkills-ndeaths);
		draw_text(-0.008, 0.01-0.0014*(num_smileys+teams)-0.0016, -0.02, text);
	}
}


void show_other_messages() {

	if (msg_params.mtime <= 0) return;
	msg_params.color.do_glColor();
	draw_text(-0.008/msg_params.size, 0.005, -0.02/msg_params.size, message.c_str(), 1.0, msg_params.bitmap);
	msg_params.mtime -= iticks;
}


void print_text_onscreen(string const &text, colorRGBA const &color, float size, int time, int priority, bool bitmap) {

	if (msg_params.mtime > 0 && msg_params.priority > priority) return; // do this before the strcpy
	message    = text;
	msg_params = text_message_params(time, size, color, priority, bitmap);
}


void print_weapon(int weapon_id) {

	print_text_onscreen(weapons[weapon_id].name, WHITE, 1.0, MESSAGE_TIME/4, 1);
}


// ***********************************
// GAME CONTROL/QUERY CODE
// ***********************************


void init_game_state() {

	if (sstates != NULL) return; // make sure this isn't called more than once
	sstates  = new player_state[num_smileys+1]; // most of the parameters are initialized in the constructor
	teaminfo = new team_info[teams];
	++sstates; // advance pointer so that camera/player can be sstates[-1]
	sstates[CAMERA_ID].name = player_name;
	vector<string> avail_smiley_names;

	for (unsigned i = 0; i < sizeof(all_smiley_names)/sizeof(string); ++i) {
		avail_smiley_names.push_back(all_smiley_names[i]);
	}
	for (int i = 0; i < num_smileys; ++i) {
		init_smiley(i);

		if (!avail_smiley_names.empty()) {
			unsigned const nid(rand()%avail_smiley_names.size());
			sstates[i].name = avail_smiley_names[nid];
			avail_smiley_names.erase(avail_smiley_names.begin() + nid);
		}
		else {
			std::ostringstream oss;
			oss << "Smiley " << i;
			sstates[i].name = oss.str();
		}
	}
	for (int i = 0; i < teams; ++i) {
		teaminfo[i].bb.x1 = -X_SCENE_SIZE + DX_VAL;
		teaminfo[i].bb.y1 = -Y_SCENE_SIZE + DY_VAL;
		teaminfo[i].bb.x2 =  X_SCENE_SIZE - DX_VAL;
		teaminfo[i].bb.y2 =  Y_SCENE_SIZE - DY_VAL;
	}
	for (unsigned i = 0; i < team_starts.size(); ++i) {
		assert(team_starts[i].index >= 0 && team_starts[i].index < teams);
		if (team_starts[i].x1 > team_starts[i].x2) swap(team_starts[i].x1, team_starts[i].x2);
		if (team_starts[i].y1 > team_starts[i].y2) swap(team_starts[i].y1, team_starts[i].y2);
		teaminfo[team_starts[i].index].bb = team_starts[i];
	}
}


int get_damage_source(int type, int index, int questioner) {

	assert(questioner >= CAMERA_ID);
	if (index == NO_SOURCE) return questioner; // hurt/killed by nature, call it a suicide
	if (type == DROWNED || type == FELL) return questioner;
	assert(index >= CAMERA_ID);
	if (type == SMILEY || type == BURNED || type == FIRE) return index;
	if (type == CAMERA) return -1;
	assert(type >= 0);
	
	if (type < CAMERA) {
		int cid(coll_id[type]);

		if (cid < 0 || cid >= num_groups) {
			for (int i = 0; i < num_groups; ++i) {
				if (obj_groups[i].type == type) {cid = i; break;}
			}
		}
		if (cid >= 0 && cid < num_groups && (unsigned)index < obj_groups[cid].max_objects()) {
			int const source(obj_groups[cid].get_obj(index).source);
			if (source == NO_SOURCE) return questioner;
			assert(source >= CAMERA_ID);
			return source;
		}
		cout << "cid = " << cid << ", index = " << index << ", max = " << obj_groups[cid].max_objects() << endl;
		assert(0);
	}
	return index;
}


bool check_underwater(int who, float &depth) { // check if player is drowning

	assert(who >= CAMERA_ID && who < num_smileys);
	point const pos(get_sstate_pos(who));
	bool const underwater(is_underwater(pos, 1, &depth));
	if (sstates == NULL) return underwater; // assert(0)?

	if (underwater) {
		++sstates[who].uw_time;

		if (game_mode) {
			int dtime(int(sstates[who].uw_time - int((float)DROWN_TIME)/fticks));

			if (dtime > 0 && (sstates[who].uw_time%TICKS_PER_SECOND) == 0) {
				float const damage(2.0*fticks*dtime);
				smiley_collision(who, who, zero_vector, pos, damage, DROWNED);
			}
		}
	}
	else {
		sstates[who].uw_time = 0;
	}
	return underwater;
}


void player_fall(int id) { // smileys and the player (camera)

	if (!game_mode) return;
	assert(id >= CAMERA_ID && id < num_smileys);
	float const zvel(sstates[id].last_zvel), dz(sstates[id].last_dz);
	float const vel(-zvel - FALL_HURT_VEL), dz2(-dz - FALL_HURT_HEIGHT*CAMERA_RADIUS);
	if (dz2 < 0.0) return;
	smiley_collision(id, id, vector3d(0.0, 0.0, -zvel), get_sstate_pos(id), 5.0*vel*vel, FELL);
}


void update_camera_velocity(vector3d const &v) {

	if (sstates == NULL) return; // assert(0)?
	sstates[CAMERA_ID].velocity = v/(TIMESTEP*fticks);
}


void init_game_mode() {

	string const str(string("Playing ") + ((game_mode == 1) ? "Deathmatch" : "Dodgeball") + " as " + player_name);
	print_text_onscreen(str, WHITE, 2.0, MESSAGE_TIME, 4);
	if (!free_for_all) teams = 0;
	free_dodgeballs(1, 1);
	init_sstate(CAMERA_ID, (game_mode == 1));
	if (game_mode == 2) init_smileys();
	
	for (int i = CAMERA_ID; i < num_smileys; ++i) {
		sstates[i].killer = NO_SOURCE; // no one
		if (game_mode == 1) init_sstate(i, 1); // ???
	}
}


void update_game_frame() {

	assert(sstates != NULL);
	sstates[CAMERA_ID].update_camera_frame();
	
	for (int i = CAMERA_ID; i < num_smileys; ++i) {
		sstates[i].update_sstate_game_frame(i);
	}
}


void player_state::update_camera_frame() {

	if (powerup_time < 0.0)   print_text_onscreen("Powerup Expired", WHITE, 1.0, MESSAGE_TIME/2, 1);
	if (powerup == PU_REGEN ) camera_health = min(MAX_REGEN_HEALTH, camera_health + 0.1f*fticks);
	if (powerup == PU_FLIGHT) flight = 1;
	kill_time += max(1, iticks);
}


void player_state::update_sstate_game_frame(int i) {

	if (powerup_time == 0) {
		powerup = -1;
	}
	else if (animate2) {
		powerup_time -= iticks;
		if (powerup_time < 0) powerup_time = 0;
	}
	if (powerup == PU_REGEN && shields > 1.0) shields = min(MAX_SHIELDS, shields + 0.075f*fticks);
		
	if (plasma_loaded && weapon == W_PLASMA) {
		plasma_size += get_fspeed_scale()*fticks*PLASMA_SIZE_INCREASE;
		plasma_size  = min(plasma_size, MAX_PLASMA_SIZE);
	}
	fire_frame = max(0,    (fire_frame - iticks));
	shields    = max(0.0f, (shields    - 0.01f*fticks));

	// check temperature for too hot/too cold
	obj_group const &objg(obj_groups[coll_id[SMILEY]]);
	if (i != CAMERA_ID && (!begin_motion || !objg.enabled)) return;
	obj_type const &objt(object_types[SMILEY]);
	point const pos(get_sstate_pos(i));
	bool const obj_enabled((i == CAMERA_ID && camera_mode == 1) || (i != CAMERA_ID && !objg.get_obj(i).disabled()));

	if (temperature < 0.75*objt.min_t) {
		float const damage(1.0*fticks/max(0.001f, (objt.min_t - temperature)/objt.min_t));
		smiley_collision(i, -2, zero_vector, pos, damage, FROZEN);
	}
	if (temperature > 0.75*objt.max_t) {
		float const damage(2.0*fticks/max(0.001f, (objt.max_t - temperature)/objt.max_t));
		smiley_collision(i, -2, zero_vector, pos, damage, BURNED);
		if (obj_enabled && (rand()&3) == 0) gen_smoke(pos);
	}
	if (atmosphere < 0.2) {
		float const damage(1.0*fticks/max(atmosphere, 0.01f));
		smiley_collision(i, -2, zero_vector, pos, damage, SUFFOCATED);
	}
	if (powerup >= 0 && powerup_time > 0 && obj_enabled) {
		add_dynamic_light(1.3, pos, get_powerup_color(powerup));
	}
	if (SMILEY_GAS && game_mode == 1 && obj_enabled && powerup == PU_SHIELD && powerup_time > INIT_PU_SH_TIME && !(rand()&31)) {
		vector3d const dir(get_sstate_dir(i)), vel(velocity*0.5 - dir*1.2);
		point const spos(pos - dir*get_sstate_radius(i)); // generate gas
		gen_arb_smoke(spos, DK_GREEN, vel, rand_uniform(0.01, 0.05), rand_uniform(0.3, 0.7), rand_uniform(0.2, 0.6), 10.0, i, GASSED, 0);
	}
}


void player_state::free_balls() {

	if (balls.empty()) return;
	obj_group &objg(obj_groups[coll_id[BALL]]);

	if (objg.enabled) {
		for (unsigned j = 0; j < balls.size(); ++j) {
			unsigned const index(balls[j]);
			assert(objg.get_obj(index).disabled());
			objg.get_obj(index).status = 0;
		}
	}
	balls.clear();
	p_weapons[W_BALL] = 0;
	p_ammo[W_BALL]    = 0;
	if (weapon == W_BALL) weapon = W_UNARMED;
}


void free_dodgeballs(bool camera, bool smileys) {

	if (sstates == NULL) return;

	for (int i = (camera ? CAMERA_ID : 0); i < (smileys ? num_smileys : 0); ++i) {
		sstates[i].free_balls();
	}
}


void gamemode_rand_appear() {

	if (!game_mode || world_mode != WMODE_GROUND) return;
	gen_smiley_or_player_pos(surface_pos, CAMERA_ID);
	camera_last_pos = surface_pos;
	free_dodgeballs(1, 0);
	init_sstate(CAMERA_ID, (game_mode == 1));
}


void change_game_mode() {

	int types[] = {HEALTH, SHIELD, POWERUP, WEAPON, AMMO};
	unsigned const ntypes(UNLIMITED_WEAPONS ? 3 : 5);

	for (unsigned i = 0; i < ntypes; ++i) {
		obj_groups[coll_id[types[i]]].set_enable(game_mode == 1);
	}
	obj_groups[coll_id[BALL]].set_enable(game_mode == 2);

	if (game_mode == 2) {
		assert(sstates != NULL);
		sstates[CAMERA_ID].switch_weapon(1, 1); // player switch to dodgeball

		for (unsigned i = 0; i < ntypes; ++i) {
			int const group(coll_id[types[i]]);
			assert(group >= 0 && group < NUM_TOT_OBJS);
			obj_groups[group].free();
		}
	}
	else if (game_mode == 0) {
		assert(coll_id[BALL] >= 0 && coll_id[BALL] < NUM_TOT_OBJS);
		obj_groups[coll_id[BALL]].free();
	}
	if (game_mode) init_game_mode(); else free_dodgeballs(1, 1);
}



