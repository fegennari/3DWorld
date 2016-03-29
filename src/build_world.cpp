// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/13/02

#include "3DWorld.h"
#include "mesh.h"
#include "tree_3dw.h"
#include "lightmap.h"
#include "textures_3dw.h"
#include "dynamic_particle.h"
#include "physics_objects.h"
#include "model3d.h"
#include "subdiv.h"
#include "player_state.h"
#include <fstream>


bool const MORE_COLL_TSTEPS       = 1; // slow
bool const SHOW_PROC_TIME         = 0;
bool const FIXED_COBJS_SWAP       = 1; // attempt to swap fixed_cobjs with coll_objects to reduce peak memory
bool const EXPLODE_EVERYTHING     = 0; // for debugging/fun
unsigned const MAX_OBJ_ON_MESH    = 4;
unsigned const BLOOD_PER_SMILEY   = 300;
unsigned const LG_STEPS_PER_FRAME = 10;
unsigned const SM_STEPS_PER_FRAME = 1;
unsigned const SHRAP_DLT_IX_MOD   = 8;
float const STAR_INNER_RAD        = 0.4;
float const ROTATE_RATE           = 25.0;


// object variables
bool printed_ngsp_warning(0), using_model_bcube(0);
int num_groups(0), used_objs(0);
unsigned next_cobj_group_id(0);
float model_czmin(czmin), model_czmax(czmax);
obj_group obj_groups[NUM_TOT_OBJS];
dwobject def_objects[NUM_TOT_OBJS];
int coll_id[NUM_TOT_OBJS] = {0};
point star_pts[2*N_STAR_POINTS];
vector<user_waypt_t> user_waypoints;
coll_obj_group fixed_cobjs;
vector<portal> portals;
vector<teleporter> teleporters;
vector<obj_draw_group> obj_draw_groups;
cube_light_src_vect sky_cube_lights, global_cube_lights;

extern bool clear_landscape_vbo, use_voxel_cobjs, tree_4th_branches, lm_alloc, reflect_dodgeballs;
extern int camera_view, camera_mode, camera_reset, begin_motion, animate2, recreated, temp_change, preproc_cube_cobjs, precip_mode;
extern int is_cloudy, num_smileys, load_coll_objs, world_mode, start_ripple, has_snow_accum, has_accumulation, scrolling, num_items, camera_coll_id;
extern int num_dodgeballs, display_mode, game_mode, num_trees, tree_mode, has_scenery2, UNLIMITED_WEAPONS, ground_effects_level;
extern float temperature, zmin, TIMESTEP, base_gravity, orig_timestep, fticks, tstep, sun_rot, czmax, czmin, dodgeball_metalness;
extern point cpos2, orig_camera, orig_cdir;
extern unsigned create_voxel_landscape, scene_smap_vbo_invalid, num_dynam_parts, init_item_counts[];
extern obj_type object_types[];
extern string cobjs_out_fn;
extern coll_obj_group coll_objects;
extern cobj_groups_t cobj_groups;
extern cobj_draw_groups cdraw_groups;
extern platform_cont platforms;
extern lightning l_strike;
extern vector<int> hmv_coll_obj;
extern char *coll_obj_file;
extern vector<point> app_spots;
extern vector<light_source> light_sources_a;
extern vector<light_source_trig> light_sources_d;
extern indir_dlight_group_manager_t indir_dlight_group_manager;
extern tree_cont_t t_trees;
extern vector<texture_t> textures;
extern reflective_cobjs_t reflective_cobjs;


int create_group(int obj_type, unsigned max_objects, unsigned init_objects,
				 unsigned app_rate, bool init_enabled, bool reorderable, bool auot_max);
void add_all_coll_objects(const char *coll_obj_file, bool re_add);
int read_coll_objects(const char *coll_obj_file);
bool write_coll_objects_file(coll_obj_group const &cobjs, string const &fn);
void gen_star_points();
int gen_game_obj(int type);
point get_sstate_pos(int id);
void reset_smoke_tex_data();
void calc_uw_atten_colors();

void write_trees_to_cobj_file(ostream &out);
void write_small_trees_to_cobj_file(ostream &out);
void write_plants_to_cobj_file(ostream &out);


void create_object_groups() {

	static bool inited(0);
	if (inited) return; // prevent multiple inits
	inited = 1;
	unsigned const num_player_blocks(num_smileys + 1); // camera + smileys
	for (int i = 0; i < NUM_TOT_OBJS; ++i) {coll_id[i] = 0;} // will be offset to -1 at the end
	
	// type, max, init, rate, enabled, reorderable, auto_max
	coll_id[SMILEY]   = create_group(SMILEY,   num_smileys, 0, 1, 1, 0, 0);
	coll_id[PRECIP]   = create_group(PRECIP,   0, 0,  40, 0, 0, 1);
	coll_id[DROPLET]  = create_group(DROPLET,  20*MAX_SPLASH_DROP, 0, 0, 0, 1, 0);
	coll_id[WDROPLET] = create_group(WDROPLET, 4000,  0, 0, 0, 1, 0);
	coll_id[SAND]     = create_group(SAND,     1000,  0, 0, 0, 1, 0);
	coll_id[DIRT]     = create_group(DIRT,     1500,  0, 0, 0, 1, 0);
	coll_id[ROCK]     = create_group(ROCK,     500,   0, 0, 0, 1, 0);
	coll_id[BLOOD]    = create_group(BLOOD,    BLOOD_PER_SMILEY*num_player_blocks, 0, 0, 0, 0, 0);
	coll_id[CHARRED]  = create_group(CHARRED,  BLOOD_PER_SMILEY*num_player_blocks, 0, 0, 0, 0, 0);
	coll_id[SFPART]   = create_group(SFPART,   4*num_player_blocks, 0, 0, 0, 0, 0); // 2 eyes, 1 nose, 1 tongue
	coll_id[CHUNK]    = create_group(CHUNK,    SMILEY_NCHUNKS*NUM_CHUNK_BLOCKS*num_player_blocks, 0, 0, 0, 0, 0);
	coll_id[SKULL]    = create_group(SKULL,    num_player_blocks, 0, 0, 0, 0, 0);
	coll_id[BALL]     = create_group(BALL,     num_dodgeballs, 0, 1, 0, 0, 0);
	coll_id[S_BALL]   = create_group(S_BALL,   200,   0, 0, 0, 1, 0);
	coll_id[ROCKET]   = create_group(ROCKET,   100,   0, 0, 0, 1, 0);
	coll_id[LANDMINE] = create_group(LANDMINE, 100,   0, 0, 0, 1, 0);
	coll_id[SEEK_D]   = create_group(SEEK_D,   50,    0, 0, 0, 1, 0);
	coll_id[STAR5]    = create_group(STAR5,    200,   0, 0, 0, 1, 0);
	coll_id[GRENADE]  = create_group(GRENADE,  200,   0, 0, 0, 1, 0);
	coll_id[CGRENADE] = create_group(CGRENADE, 20,    0, 0, 0, 1, 0);
	coll_id[SHRAPNEL] = create_group(SHRAPNEL, 8000,  0, 0, 0, 1, 0);
	coll_id[SHELLC]   = create_group(SHELLC,   5*object_types[SHELLC].lifetime,  0, 0, 0, 1, 0);
	coll_id[LEAF]     = create_group(LEAF,     3000,  0, 0, 1, 1, 0);
	coll_id[HEALTH]   = create_group(HEALTH,   init_item_counts[0], 0, 1, 0, 1, 0);
	coll_id[SHIELD]   = create_group(SHIELD,   init_item_counts[1], 0, 1, 0, 1, 0);
	coll_id[POWERUP]  = create_group(POWERUP,  init_item_counts[2], 0, 1, 0, 1, 0);
	coll_id[WEAPON]   = create_group(WEAPON,   init_item_counts[3], 0, 1, 0, 1, 0);
	coll_id[AMMO]     = create_group(AMMO,     init_item_counts[4], 0, 1, 0, 1, 0);
	coll_id[WA_PACK]  = create_group(WA_PACK,  50,    0, 0, 0, 1, 0);
	coll_id[FRAGMENT] = create_group(FRAGMENT, 1200,  0, 0, 0, 1, 0);
	coll_id[PARTICLE] = create_group(PARTICLE, 800,   0, 0, 0, 1, 0);
	coll_id[SAWBLADE] = create_group(SAWBLADE, 50,    0, 0, 0, 1, 0);
	coll_id[PLASMA]   = create_group(PLASMA,   150,   0, 0, 0, 1, 0); // Note: create plasma group last since it uses a special shader during drawing
	for (int i = 0; i < NUM_TOT_OBJS; ++i) {coll_id[i] -= 1;} // offset by -1
}


int create_group(int obj_type, unsigned max_objects, unsigned init_objects,
				 unsigned app_rate, bool init_enabled, bool reorderable, bool auto_max)
{
	if (num_groups > NUM_TOT_OBJS) {
		cerr << "Error: Exceeded max of " << NUM_TOT_OBJS << " object groups." << endl;
		exit(1);
	}
	if (obj_type >= NUM_TOT_OBJS) {
		cerr << "Error: Illegal object type: " << obj_type << "." << endl;
		assert(0);
	}
	obj_groups[num_groups].create(obj_type, max_objects, init_objects, app_rate, init_enabled, reorderable, auto_max);
	return ++num_groups;
}


void dwobject::add_obj_dynamic_light(int index) const {

	switch(type) {
	case PLASMA:
		add_dynamic_light(min(3.5, 45.0*get_true_radius()), pos, get_plasma_color(init_dir.x));
		break;
	case ROCKET:
		add_dynamic_light(0.5, pos, ORANGE);
		break;
	case SEEK_D:
		add_dynamic_light(0.6, pos, colorRGBA(1.0, 0.25, 0.0, 1.0)); // red-orange
		break;
	case LANDMINE:
		if (time > 5) {
			float const radius(object_types[type].radius);
			float const sensor_height(get_landmine_sensor_height(radius, time) + 0.15*radius);
			add_dynamic_light(0.3, (pos + point(0.0, 0.0, sensor_height)), get_landmine_light_color(time));
		}
		break;
	case SHRAPNEL:
	case PARTICLE:
		{
			if (type == SHRAPNEL && direction == W_GRENADE && (index % SHRAP_DLT_IX_MOD) != 0) break; // optimization hack
			float stime;
			colorRGBA const scolor(get_glowing_obj_color(pos, time, object_types[type].lifetime, stime, (type == SHRAPNEL), 0));
			if (stime < 1.0) add_dynamic_light(0.2, pos, scolor);
		}
		break;
	case BALL:
		{
			colorRGBA colors[NUM_DB_TIDS] = {BLUE, colorRGBA(1.0, 0.5, 0.5, 1.0), colorRGBA(0.5, 1.0, 0.5, 1.0)};
			colorRGBA const color((game_mode == 2) ? colors[index%NUM_DB_TIDS] : colorRGBA(-1.0, -1.0, -1.0, 1.0));
			add_dynamic_light(0.8, pos, color);
		}
		break;
	case CHUNK:
		if (flags & TYPE_FLAG) { // charred
			float stime;
			colorRGBA const scolor(get_glowing_obj_color(pos, time, object_types[type].lifetime, stime, 1, 0));
			if (stime < 1.0) add_dynamic_light(0.5, pos, scolor);
		}
	}
}


bool is_rain_enabled() {return (temperature >  W_FREEZE_POINT && precip_mode != 0);}
bool is_snow_enabled() {return (temperature <= W_FREEZE_POINT && precip_mode != 0);}
int get_precip_type () {return ((temperature > RAIN_MIN_TEMP) ? RAIN : ((temperature > SNOW_MAX_TEMP) ? HAIL : SNOW));}

int obj_group::get_ptype() const {return ((flags & PRECIPITATION) ? get_precip_type() : type);}


void dwobject::update_precip_type() {

	int const ptype(get_precip_type());
	if (type == ptype) return;
	if (status == 4) {status = 1; flags  = 0;}
	type   = ptype;
	health = def_objects[type].health;
}


void process_platforms_falling_moving_and_light_triggers() {

	if (!animate2) return; // no updates
	int const start_i((camera_mode == 1) ? CAMERA_ID : 0);
	int const end_i(obj_groups[coll_id[SMILEY]].is_enabled() ? num_smileys : 0);

	if (!coll_objects.platform_ids.empty()) { // update platforms
		for (int i = start_i; i < end_i; ++i) {platforms.check_activate(get_sstate_pos(i), CAMERA_RADIUS, i);}
		platforms.advance_timestep();
	}
	proc_moving_cobjs(); // Note: depends on platforms and uses+modified the cobj BVHs
	check_falling_cobjs();
	build_static_moving_cobj_tree();

	for (auto l = light_sources_d.begin(); l != light_sources_d.end(); ++l) { // update scene lights
		for (int i = start_i; i < end_i; ++i) {l->check_activate(get_sstate_pos(i), CAMERA_RADIUS, i);}
		l->advance_timestep();
	}
}


void object_line_coll(dwobject &obj, point const &old_pos, float radius, unsigned obj_index, int &cindex) {

	vector3d cnorm(zero_vector);
	point cpos(obj.pos);

	//if (coll_objects[cindex].line_int_exact(old_pos, pos, t, cnorm)) {
	if (check_coll_line_exact(old_pos, obj.pos, cpos, cnorm, cindex)) { // slower, but more correct
		assert(cnorm != zero_vector);
		obj.flags |= OBJ_COLLIDED;
		obj.pos    = cpos; // move it to collision point
		bool coll(0);
		if (cindex >= 0) coll     = (obj.check_vert_collision(obj_index, 1, 0, NULL, all_zeros, 0, 0, cindex) != 0);
		if (!coll)       obj.pos += cnorm*(0.99*radius); // move so it only slightly collides
		assert(!is_nan(obj.pos));
	}
}


void set_global_state() {

	camera_view = 0;
	used_objs   = 0;
	is_cloudy   = (precip_mode > 0);
}


void process_groups() {

	if (animate2) {advance_physics_objects();}

	if (display_mode & 0x0200) {
		d_part_sys.create_particles(num_dynam_parts, 1);
		d_part_sys.apply_physics();
		d_part_sys.add_light();
	}
	set_global_state();
	if (num_groups == 0) return; // groups not enabled
	RESET_TIME;
	unsigned num_objs(0);
	static int camera_follow(0);
	static unsigned scounter(0);
	int const lcf(camera_follow);
	++scounter;
	camera_follow = 0;
	build_cobj_tree(1, 0); // could also do after group processing
	
	for (int i = 0; i < num_groups; ++i) {
		obj_group &objg(obj_groups[i]);
		objg.preproc_this_frame();
		if (!objg.enabled) continue;
		unsigned const flags(objg.flags);
		obj_type const &otype(object_types[objg.type]);
		bool const precip((flags & PRECIPITATION) != 0);
		int const type(objg.get_ptype());

		if (!objg.temperature_ok()) {
			if (temp_change) {
				if (type == SMILEY) free_dodgeballs(0, 1);
				objg.remove_reset_cobjs();
				objg.init_group();
			}
			continue;
		}
		if (!animate2) continue;

		if (!begin_motion) {
			objg.remove_reset_cobjs();
			if (flags & (JUST_INIT | WAS_ADVANCED)) objg.init_group();
			continue;
		}
		float const radius(otype.radius);
		bool const large_radius(objg.large_radius());
		collision_func coll_func(NULL);

		switch (type) {
		case CAMERA:   coll_func = camera_collision;    break;
		case SMILEY:   coll_func = smiley_collision;    break;
		case LANDMINE: coll_func = landmine_collision;  break;
		case HEALTH:   coll_func = health_collision;    break;
		case SHIELD:   coll_func = shield_collision;    break;
		case POWERUP:  coll_func = powerup_collision;   break;
		case WEAPON:   coll_func = weapon_collision;    break;
		case AMMO:     coll_func = ammo_collision;      break;
		case WA_PACK:  coll_func = pack_collision;      break;
		case S_BALL:   coll_func = sball_collision;     break;
		case BALL:     coll_func = dodgeball_collision; break;
		case SKULL:    coll_func = skull_collision;     break;
		case SAWBLADE: coll_func = sawblade_collision;  break;
		}
		//cout << "group %d %d\n", i, GET_DELTA_TIME);
		RESET_TIME;
		unsigned gen_count(0);
		num_objs = 0;
		float const fticks_max(min(4.0f, fticks)); // clamp effective fticks so that we don't slow the framerate down even more
		unsigned app_rate(unsigned(((float)objg.app_rate)*fticks_max));
		if (objg.app_rate > 0 && fticks > 0 && app_rate == 0) {app_rate = 1;}
		float const time(TIMESTEP*fticks_max), grav_dz(min(otype.terminal_vel*time, base_gravity*GRAVITY*time*time*otype.gravity));
		size_t const max_objs(objg.max_objects());
		bool const reflective(reflect_dodgeballs && type == BALL && enable_all_reflections()); // Note: cobjs only have a lifetime of one frame
		cobj_params cp(otype.elasticity, otype.color, reflective, 1, coll_func, -1, otype.tid, 1.0, 0, 0);
		if (reflective) {cp.metalness = dodgeball_metalness; cp.tscale = 0.0; cp.color = WHITE; cp.spec_color = WHITE; cp.shine = 100.0;} // reflective metal sphere

		for (size_t jj = 0; jj < max_objs; ++jj) {
			unsigned const j(unsigned((type == SMILEY) ? (jj + scounter)%objg.max_objects() : jj)); // handle smiley permutation
			dwobject &obj(objg.get_obj(j));
			if (large_radius && obj.coll_id >= 0) remove_reset_coll_obj(obj.coll_id);
			if (obj.status == OBJ_STAT_RES) continue; // ignore
			point &pos(obj.pos);

			if (obj.status == 0) {
				if (gen_count >= app_rate || !(flags & WAS_ADVANCED))      continue;
				if (type == BALL && (game_mode != 2 || UNLIMITED_WEAPONS)) continue; // not in dodgeball mode
				++gen_count;
				if (precip && temperature >= WATER_MAX_TEMP) continue; // skip it
				dwobject new_obj(def_objects[type]);
				int const ret(objg.get_next_predef_obj(new_obj, j));
				if (ret == 0) continue; // skip this object (no slots available)
				obj = new_obj;
				
				if (ret == 1) { // use a predefined object
					assert(type != SMILEY); // use an appearance spot for a smiley
				}
				else { // standard random generation
					if (type == SMILEY) {
						if (!gen_smiley_or_player_pos(pos, j)) {
							if (!printed_ngsp_warning) cout << "No good smiley pos." << endl;
							printed_ngsp_warning = 1;
						}
					}
					else {
						gen_object_pos(pos, otype.flags);
						assert(!is_nan(pos));
						vadd_rand(obj.velocity, 1.0);
					}
					if (type != SMILEY && (otype.flags & NO_FALL)) {
						pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 0) + radius;
					}
					if (type == POWERUP || type == WEAPON || type == AMMO) {
						obj.direction = (unsigned char)gen_game_obj(type);
					}
				}
				if (otype.flags & OBJ_IS_FLAT) {
					obj.init_dir = signed_rand_vector_norm();
					obj.angle    = signed_rand_float();
				}
				if (type == SNOW) {obj.angle = rand_uniform(0.7, 1.3);} // used as radius
			}
			if (precip) obj.update_precip_type();
			unsigned char const obj_flags(obj.flags);
			int const orig_status(obj.status);
			obj.flags &= ~PLATFORM_COLL;
			++used_objs;
			++num_objs;

			if (obj.health < 0.0) { // can get here for smileys?
				obj.status = 0;
			}
			else if (type == SMILEY) {
				advance_smiley(obj, j);
			}
			else {
				if (obj.time >= 0) {
					if (type == PLASMA && obj.velocity.mag_sq() < 1.0) {
						obj.disable(); // plasma dies when it stops
					}
					else {
						if (large_radius) {maybe_teleport_object(obj.pos, radius, NO_SOURCE);} // teleport!
						point const old_pos(pos); // after teleporting
						unsigned spf(1);
						int cindex(-1);

						// What about rolling objects (type_flags & OBJ_ROLLS) on the ground (status == 3)?
						if (obj.status == 1 && is_over_mesh(pos) && !((obj_flags & XY_STOPPED) && (obj_flags & Z_STOPPED))) {
							if (obj.flags & CAMERA_VIEW) { // smaller timesteps if camera view
								spf = 4*LG_STEPS_PER_FRAME;
							}
							else if (type == PLASMA || type == BALL || type == SAWBLADE) {
								spf = 3*LG_STEPS_PER_FRAME;
							}
							else if (type == ROCKET || type == SEEK_D) {
								spf = 2*LG_STEPS_PER_FRAME;
							}
							else if (large_radius /*|| type == STAR5 || type == SHELLC*/ || type == FRAGMENT) {
								spf = LG_STEPS_PER_FRAME;
							}
							else if (type == SHRAPNEL) {
								spf = max(1, min(((obj.direction == W_GRENADE) ? 4 : 20), int(0.2*obj.velocity.mag())));
							}
							else if (type == PRECIP || (flags & PRECIPITATION)) {
								spf = 1;
							}
							else {
								spf = SM_STEPS_PER_FRAME;
							}
							if (MORE_COLL_TSTEPS && obj.status == 1 && spf < LG_STEPS_PER_FRAME && pos.z < czmax && pos.z > czmin) {
								point pos2(pos + obj.velocity*time); // makes precipitation slower, but collision detection is more correct
								pos2.z -= grav_dz; // maybe want to try with and without this?
								// Note: we only do the line intersection test if the object moves by more than its radius this frame (static leaves don't)
								// Note: could also test pos.z > v_collision_matrix[y][x].zmax
								if (!dist_less_than(pos, pos2, radius)) {check_coll_line(pos, pos2, cindex, -1, 0, 0);} // return value is unused
							}
							assert(spf > 0);

							if (spf > 1) {
								assert(fticks > 0.0);
								orig_timestep = TIMESTEP; // incremental multistep object advance
								TIMESTEP     /= float(spf);
								tstep         = TIMESTEP*fticks;
								point const obj_pos(obj.pos);
								
								for (unsigned k = 0; k < spf; ++k) {
									obj.advance_object(!recreated, k, j);
									if (obj.status != 1)    break; // no longer airborne
									if (obj.pos == obj_pos) break; // stopped
								}
								TIMESTEP = orig_timestep;
								tstep    = time;
							}
						}
						if (spf == 1) obj.advance_object(!recreated, 0, j);
						obj.verify_data();
						
						if (!obj.disabled() && cindex >= 0 && !large_radius && spf < LG_STEPS_PER_FRAME) { // test collision with this cobj
							object_line_coll(obj, old_pos, radius, j, cindex);
						}
					} // not plasma
				} // obj.time < 0
				else obj.time = 0;
			} // not smiley
			if (!obj.disabled()) {
				update_deformation(obj);
				
				if ((otype.flags & OBJ_IS_FLAT) || type == SHELLC || type == SHRAPNEL || type == STAR5 || type == LEAF || type == SAWBLADE ||
					(type == FRAGMENT && (obj.flags & TYPE_FLAG))) // shatterable fragments
				{
					float const vmag(fabs(obj.velocity.z)); // rotate
					
					if (obj.status == 1 && vmag > 0.5 && !(obj.flags & (STATIC_COBJ_COLL | OBJ_COLLIDED))) {
						float const rr((type == LEAF) ? 0.25 : 1.0);
						obj.angle += fticks*(TIMESTEP/DEF_TIMESTEP)*rr*ROTATE_RATE*sqrt(vmag); // rotate while airborne
					}
				}
				if (large_radius) {
					if (type != LANDMINE) {
						float const r2((otype.flags & COLL_DESTROYS) ? 0.25*radius : radius);
						collision_detect_large_sphere(pos, r2, obj_flags);
					}
					if (type != CHUNK && (type != LANDMINE || !obj.lm_coll_invalid())) {
						if (type == BALL) {cp.tid = dodgeball_tids[(game_mode == 2) ? (j%NUM_DB_TIDS) : 0];}
						cp.cf_index = j;
						obj.coll_id = add_coll_sphere(pos, radius, cp, -1, 0, reflective);
					}
				}
				if (obj_flags & CAMERA_VIEW) {
					if (camera_reset) {
						obj.flags ^= CAMERA_VIEW;
						//camera_reset = 0;
					}
					else {
						cpos2         = pos;
						camera_view   = 1;
						camera_follow = 1;
						
						if (type == SEEK_D) { // follow the player's direction
							float const vmag(obj.velocity.mag());
							obj.orientation = cview_dir;
							obj.velocity    = cview_dir*vmag;
						}
					}
				}
				if (!reflective) {obj.add_obj_dynamic_light(j);}
			} // !obj.disabled()
			if (!recreated && orig_status != 0 && obj.status != 1 && obj.status != OBJ_STAT_RES) {
				if ((precip || type == BLOOD || type == WDROPLET) && obj.status == 0) {
					accumulate_object(pos, type, 1.0);
				}
				else if (type == SKULL && obj.status == 0) { // create skull fragments
					unsigned const num((rand()&7) + 12);
					bool const burned(obj.direction == 1);

					for (unsigned o = 0; o < num; ++o) {
						point const fpos(obj.pos + signed_rand_vector(radius));
						gen_fragment(fpos, zero_vector, 0.5, 0.0, (burned ? GRAY_BLACK : LT_GRAY), PLASTER_TEX, 1.0, obj.source, 0);
					}
				}
				else if ((otype.flags & EXPL_ON_COLL) || (obj.status == 0 && (otype.flags & OBJ_EXPLODES))) {
					obj.status = 0;
					if (otype.flags & EXPL_ON_COLL) collision_detect_large_sphere(pos, radius, flags);
					blast_radius(pos, type, j, obj.source, 0);
					gen_smoke(pos);
					gen_fire(pos, ((type == PLASMA) ? obj.init_dir.x : rand_uniform(0.4, 1.0)), obj.source);
					if (type == LANDMINE) {gen_landmine_scorch(obj.pos);}
				}
			}
			if (type == LANDMINE && obj.status == 1 && !(obj.flags & (STATIC_COBJ_COLL | PLATFORM_COLL))) obj.time = 0; // don't start time until it lands
		} // for jj
		objg.flags |= WAS_ADVANCED;
		if (num_objs > 0 && (SHOW_PROC_TIME /*|| type == SMILEY*/)) {cout << "type = " << type << ", num = " << num_objs << " "; PRINT_TIME("Process");}
	} // for i
	temp_change = 0;
	recreated   = 0;

	if (lcf == 1 && camera_follow == 0) {
		camera_reset  = 1;
		camera_origin = orig_camera;
		cview_dir     = orig_cdir;
		//reset_camera_pos();
	}
	if (SHOW_PROC_TIME) {PRINT_TIME("Final");}
}


void gen_scene(int generate_mesh, int gen_trees, int keep_sin_table, int update_zvals, int rgt_only) {
	
	cout << "Generating Scene..." << endl;
	RESET_TIME;
	static int st_valid(0);
	bool const inf_terrain(world_mode == WMODE_INF_TERRAIN);
	scene_smap_vbo_invalid = 2; // needed to force smap update - full rebuild of shadowers

	if (!st_valid) {
		keep_sin_table = 0;
		st_valid = 1;
	}
	l_strike.time = LITNING_TIME; // reset lightning
	kill_current_raytrace_threads();
	if (!keep_sin_table) {clear_tiled_terrain();}
	calc_uw_atten_colors();

	if (generate_mesh) {
		if (generate_mesh != 2) {
		  gen_mesh(0, keep_sin_table, update_zvals);
		  PRINT_TIME("Surface generation");
		}
		gen_tex_height_tables();
		clear_landscape_vbo = 1;
	}
	compute_matrices();
	PRINT_TIME("Matrix generation");
	
	if (generate_mesh) {
		create_landscape_texture();
		PRINT_TIME("Landscape Texture generation");
		has_snow_accum   = 0;
		has_accumulation = 0;
		start_ripple     = 0;
	}
	compute_volume_matrix();
	calc_motion_direction();
	PRINT_TIME("Volume+Motion matrix generation");
	
	if (generate_mesh && !inf_terrain && create_voxel_landscape == 1) {
		gen_voxel_landscape();
		PRINT_TIME("Voxel Landscape Generation");
	}
	if (num_trees > 0) {
		if (!inf_terrain && gen_trees) {
			regen_trees(1);
			PRINT_TIME("Tree generation");
		}
		else {
			delete_trees();
		}
	}
	if (!inf_terrain) {
		gen_scenery();
		PRINT_TIME("Scenery generation");
	}
	add_all_coll_objects(coll_obj_file, (num_trees == 0));
	PRINT_TIME("Collision object addition");

	if (!inf_terrain && !rgt_only) {
		calc_watershed();
		PRINT_TIME("Water generation");
	}
	if (!inf_terrain && !scrolling) {
		create_waypoints(user_waypoints);
		PRINT_TIME("Waypoint Creation");
	}
	reanimate_objects(); // allow stationary/stuck objects to move about the new terrain (fast so no timing)
	invalidate_snow_coverage();

	unsigned char sflags(0);
	float const lf(fabs(sun_rot/PI - 1.0)); // light_factor
	if (!scrolling || lf >= 0.4) {sflags |= SUN_SHADOW;}
	if (!scrolling || lf <= 0.6) {sflags |= MOON_SHADOW;}
	calc_visibility(sflags);
	PRINT_TIME("Visibility calculation");

	if (!inf_terrain) {
		if (generate_mesh) {gen_grass();}
		if (generate_mesh || gen_trees) {reset_smoke_tex_data();}
	}
	//triangle_counter(1).proc_cobjs(coll_objects);
}


void shift_point_vector(vector<point> &pts, vector3d const &vd) {
	for (unsigned i = 0; i < pts.size(); ++i) pts[i] += vd;
}

void shift_all_cobjs(vector3d const &vd) {
	for (unsigned i = 0; i < coll_objects.size(); ++i) {coll_objects[i].shift_by(vd);}
}


void shift_all_objs(vector3d const &vd) {

	shift_all_cobjs(vd);
	shift_hmv(vd);
	shift_trees(vd);
	shift_small_trees(vd);
	shift_scenery(vd);
	shift_water_springs(vd);
	shift_other_objs(vd);
	shift_light_sources(vd);
	platforms.shift_by(vd);
	shift_waypoints(vd); // is this correct?
	for (vector<user_waypt_t>::iterator i = user_waypoints.begin(); i != user_waypoints.end(); ++i) {i->pos += vd;}
	//shift_point_vector(app_spots, vd); // what if an appearance spot shifts off the map?

	if (begin_motion) {
		for (int i = 0; i < num_groups; ++i) {obj_groups[i].shift(vd);}
	}
}


void coll_obj::shift_by(vector3d const &vd, bool force, bool no_texture_offset) {

	if (!fixed && !force) return;
	for (unsigned j = 0; j < unsigned(npoints); ++j) {points[j] += vd;}
	cube_t::translate(vd);
	if (!no_texture_offset && cp.tscale != 0.0 && !(cp.flags & COBJ_WAS_CUBE)) {texture_offset -= vd;}
	if (cgroup_id >= 0) {cobj_groups.invalidate_group(cgroup_id);} // force recompute of center of mass, etc.
}

void coll_obj::move_cobj(vector3d const &vd, bool update_colls) {

	if (update_colls) {remove_coll_object(id, 0);}
	shift_by(vd); // move object
	if (update_colls) {re_add_coll_cobj(id, 0);}
	if (update_colls && (is_rain_enabled() || is_wet() || is_reflective())) {check_indoors_outdoors();} // update indoor/outdoor state if it's raining, reflective, or if already wet
}


void free_all_coll_objects() {

	// Note: all cobjs should have been removed from coll_objects/cobj_manager at the point,
	//       but the various scene objects could still reference them and need to be cleared
	free_scenery_cobjs();
	remove_small_tree_cobjs();
	remove_tree_cobjs();
	bool have_fixed_cobjs(0);
	
	for (unsigned i = 0; i < coll_objects.size(); ++i) {
		have_fixed_cobjs |= coll_objects[i].fixed;
		coll_objects[i].waypt_id = -1;
	}
	for (unsigned i = 0; i < hmv_coll_obj.size(); ++i) {
		remove_reset_coll_obj(hmv_coll_obj[i]);
	}
	if (begin_motion) {
		for (int i = 0; i < num_groups; ++i) {
			if (obj_groups[i].enabled) {obj_groups[i].remove_reset_cobjs();}
		}
	}
	remove_coll_object(camera_coll_id); // not necessary?
	purge_coll_freed(1);

	if (!have_fixed_cobjs) { // Note: if there are fixed cobjs, these maps will be nonempty
		assert(coll_objects.dynamic_ids.empty());
		assert(coll_objects.drawn_ids.empty());
		assert(coll_objects.platform_ids.empty());
	}
	if (!lm_alloc) { // if the lighting has already been computed, we can't change czmin/czmax/get_zval()/get_zpos()
		czmin = model_czmin; // reset zmin/zmax to original values before cobjs were added
		czmax = model_czmax;
	}
	//reflective_cobjs.clear(); // no longer needed or correct - reflective cobjs are tracked automatically and should not be manually cleared
	reflective_cobjs.mark_faces_invalid(); // force recreation of cube maps
}


void check_contained_cube_sides() {

	for (coll_obj_group::iterator i = coll_objects.begin(); i != coll_objects.end(); ++i) {
		if (!i->fixed || i->may_be_dynamic() || i->type != COLL_CUBE || i->is_semi_trans()) continue;

		for (unsigned dim = 0; dim < 3; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				if (i->cp.surfs & EFLAGS[dim][dir]) continue;
				if (check_face_containment(*i, dim, dir, (i - coll_objects.begin()))) {i->cp.surfs |= EFLAGS[dim][dir];} // set this flag bit
			}
		}
	}
}

void flag_cobjs_indoors_outdoors() {

	for (coll_obj_group::iterator i = coll_objects.begin(); i != coll_objects.end(); ++i) {
		if (i->fixed && !i->no_draw()) {i->check_indoors_outdoors();}
	}
}


void coll_obj_group::clear_ids() {
	
	dynamic_ids.clear();
	drawn_ids.clear();
	platform_ids.clear();
}


void coll_obj_group::clear() { // unused, but may be useful
	
	has_lt_atten = has_voxel_cobjs = 0;
	vector<coll_obj>::clear();
	clear_ids();
}


void coll_obj_group::finalize() {

	process_negative_shapes(); // must be first because requires an unmodified ordering of shapes
	bool has_cubes(0), any_drawn(0);

	for (coll_obj_group::const_iterator i = begin(); i != end(); ++i) {
		has_cubes |= (i->type == COLL_CUBE);
		any_drawn |= i->cp.draw;
	}
	if (has_cubes) { // Note: important to do this test on large polygon-only models
		remove_overlapping_cubes((preproc_cube_cobjs == 1) ? NON_DEST : SHATTERABLE); // always remove overlaps with >= shatterable cobjs

		if (preproc_cube_cobjs) {
			merge_cubes(); // and alpha sort
			subdiv_cubes();
		}
		check_cubes(); // sanity check, should be last
	}
	if (any_drawn) sort_cobjs_for_rendering();
}


void add_all_coll_objects(const char *coll_obj_file, bool re_add) {

	static int init(0);

	if (!init) {
		if (load_coll_objs) {
			if (!read_coll_objects(coll_obj_file)) {exit(1);}
			fixed_cobjs.finalize();
			bool const has_voxel_cobjs(gen_voxels_from_cobjs(fixed_cobjs));
			unsigned const ncobjs(fixed_cobjs.size());
			RESET_TIME;
			
			if (!FIXED_COBJS_SWAP || has_voxel_cobjs || !swap_and_set_as_coll_objects(fixed_cobjs)) {
				if (ncobjs > 2*coll_objects.size()) {reserve_coll_objects(coll_objects.size() + 1.1*ncobjs);} // reserve with 10% buffer
				
				for (unsigned i = 0; i < ncobjs; ++i) {
					if (fixed_cobjs[i].cp.cobj_type == COBJ_TYPE_VOX_TERRAIN) continue; // skip it
					fixed_cobjs[i].add_as_fixed_cobj(); // don't need to remove it
				}
			}
			PRINT_TIME(" Add Fixed Cobjs");
			fixed_cobjs.clear();
			remove_excess_cap(fixed_cobjs); // free the memory
			if (!cobjs_out_fn.empty()) {write_coll_objects_file(coll_objects, cobjs_out_fn);} // after fixed cobjs processing
		}
		init = 1;
	}
	else {
		for (unsigned i = 0; i < coll_objects.size(); ++i) {coll_objects[i].re_add_coll_cobj(i);}
	}
	purge_coll_freed(1);
	add_shape_coll_objs();
	
	if (re_add) {
		if (num_trees == 0) {
			if (tree_mode & 1) {add_tree_cobjs();} // multiple adds?
			if (tree_mode & 2) {add_small_tree_coll_objs();}
		}
		if (has_scenery2) {add_scenery_cobjs();}
	}
	bool const verbose(!scrolling);
	if (verbose) {cobj_stats();}
	build_cobj_tree(0, verbose);
	check_contained_cube_sides();
	flag_cobjs_indoors_outdoors();
}


int read_error(FILE *fp, const char *param, const char *filename) {
	cout << "*** Error reading " << param << " from file '" << filename << "'. ***" << endl;
	fclose(fp);
	return 0;
}
int read_error(FILE *fp, const char *param, const char *filename, unsigned line_num) {
	cout << "*** Error reading " << param << " from file '" << filename << "' on line " << line_num << ". ***" << endl;
	fclose(fp);
	return 0;
}

void check_layer(bool has_layer) {
	if (!has_layer) cout << "* Warning: Shape found before a layer specification in config file. Using default layer." << endl;
}


// Note: always called with cobjs == fixed_cobjs
void coll_obj::add_to_vector(coll_obj_group &cobjs, int type_) {

	type = type_;
	id   = (unsigned)cobjs.size();
	check_if_cube();
	set_npoints();
	if (type == COLL_POLYGON) {assert(npoints >= 3); norm = get_poly_norm(points);}

	if (dgroup_id >= 0) { // grouped cobj
		calc_bcube(); // may be needed to get center point if this is the parent
		bool const is_parent(cdraw_groups.set_parent_or_add_cobj(*this));
		if (!is_parent) return; // a child draw cobj, don't add it to cobjs
	}
	cobjs.push_back(*this);
}


void coll_obj::check_if_cube() {

	if (type != COLL_POLYGON || thickness == 0.0 || npoints != 4) return;
	cube_t bb(points, 4);
	unsigned zdim(0), nz(0);
	float const smax(bb.max_len());
	float const tolerance(1.0E-6*smax);

	for (unsigned j = 0; j < 3; ++j) {
		if ((bb.d[j][0] - bb.d[j][1]) < tolerance) {
			zdim = j;
			++nz;
		}
	}
	if (nz != 1) return;
	unsigned const dims[2] = {(zdim+1)%3, (zdim+2)%3};
	
	for (unsigned i = 0; i < 4 ; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			bool on_edge(0);
			
			for (unsigned k = 0; k < 2 && !on_edge; ++k) {
				if (fabs(points[i][dims[j]] - bb.d[dims[j]][k]) < tolerance) on_edge = 1;
			}
			if (!on_edge) return; // not a cube
		}
	}
	type = COLL_CUBE;
	copy_from(bb);
	d[zdim][0] -= 0.5*thickness;
	d[zdim][1] += 0.5*thickness;
}


void copy_polygon_to_cobj(polygon_t const &poly, coll_obj &cobj) {

	cobj.npoints = (short)poly.size();
	for (int j = 0; j < cobj.npoints; ++j) {cobj.points[j] = poly[j].v;}
}


void copy_tquad_to_cobj(coll_tquad const &tquad, coll_obj &cobj) {

	cobj.npoints = tquad.npts;
	for (int j = 0; j < cobj.npoints; ++j) {cobj.points[j] = tquad.pts[j];}
}


void maybe_reserve_fixed_cobjs(size_t size) {

	unsigned const ncobjs(fixed_cobjs.size());
	if (size > 2*ncobjs) {fixed_cobjs.reserve(ncobjs + (FIXED_COBJS_SWAP ? 1.1 : 1.0)*size);} // reserve to the correct size
}


void add_polygons_to_cobj_vector(vector<coll_tquad> const &ppts, coll_obj const &cobj, int const *group_ids, unsigned char cobj_type) {

	coll_obj poly(cobj);
	poly.cp.cobj_type = cobj_type;
	if (cobj_type != COBJ_TYPE_STD) poly.cp.draw = 0;
	maybe_reserve_fixed_cobjs(ppts.size());
	
	for (vector<coll_tquad>::const_iterator i = ppts.begin(); i != ppts.end(); ++i) {
		unsigned const npts(i->npts);
		assert(npts == 3 || npts == 4);
		if (!i->is_valid()) continue; // invalid zero area polygon - skip
		copy_tquad_to_cobj(*i, poly);
		vector3d const norm(get_poly_norm(poly.points));

		if (norm == zero_vector) {
			static bool had_zero_area_warning(0);

			if (!had_zero_area_warning) {
				cout << "* Warning: Ignoring zero area polygon." << endl;
				had_zero_area_warning = 1;
			}
			continue;
		}
		if (i->color.c[3] > 0) { // we have a valid override color, so set the color to this and reset the texture id
			poly.cp.color = i->color.get_c4();
			poly.cp.tid   = -1;
		}
		if (group_ids) poly.group_id = group_ids[get_max_dim(norm)];
		poly.add_to_vector(fixed_cobjs, COLL_POLYGON); // 3 or 4 point convex polygons only
	}
}


void create_xyz_groups(int *group_ids, bool use_vbo) {

	for (unsigned i = 0; i < 3; ++i) {
		group_ids[i] = (int)obj_draw_groups.size();
		obj_draw_groups.push_back(obj_draw_group(use_vbo));
	}
}


void read_or_calc_zval(FILE *fp, point &pos, float interp_rad, float radius, geom_xform_t const &xf) {

	pos.z = 0.0;
	bool const interpolate(fscanf(fp, "%f", &pos.z) != 1);
	xf.xform_pos(pos); // better not try to rotate z when interpolating
	if (interpolate) pos.z = interpolate_mesh_zval(pos.x, pos.y, interp_rad, 0, 0) + radius;
}


string read_filename(FILE *fp, unsigned &line_num) {

	assert(fp);
	string str;
	bool in_quote(0);

	while (1) {
		int const c(getc(fp));
		if (c == '\0' || c == EOF) {return str;} // end of file
		else if (c == '"') {in_quote ^= 1;} // quote
		else if (isspace(c)) { // whitespace character
			if (c == '\n') {++line_num;}
			if (in_quote) {str.push_back(c);} // part of quoted string
			else if (!str.empty()) {return str;} // not leading whitespace
		}
		else {str.push_back(c);}
	}
	assert(0); // never gets here
	return str;
}


bool read_texture(char const *const str, unsigned line_num, int &tid, bool is_normal_map, bool invert_y=0) {

	tid = get_texture_by_name(std::string(str), is_normal_map, invert_y);
	
	if (tid >= 0 && (unsigned)tid >= textures.size()) {
		cout << "Illegal texture on line " << line_num << ": " << tid << ", max is " << textures.size()-1 << endl;
		return 0;
	}
	return 1;
}


void add_model_polygons_to_cobjs(vector<coll_tquad> const &ppts, coll_obj &cobj, bool group_cobjs, bool use_vbo, unsigned cobj_type, bool has_layer, float scale) {

	int group_ids[3] = {-1, -1, -1}; // one for each primary dim
	if (group_cobjs) {create_xyz_groups(group_ids, use_vbo);}
	check_layer(has_layer);
	float const prev_thick(cobj.thickness);
	cobj.thickness *= scale;
	if (cobj.thickness == 0.0) {cobj.thickness = MIN_POLY_THICK;} // optional - will be set to this value later anyway
	else if (group_cobjs) {cobj.thickness = min(cobj.thickness, MIN_POLY_THICK);} // grouping code requires thin polygons
	add_polygons_to_cobj_vector(ppts, cobj, group_ids, cobj_type);
	cobj.group_id   = -1; // reset
	cobj.thickness  = prev_thick;
}

// returns error string
string add_loaded_model(vector<coll_tquad> const &ppts, coll_obj cobj, int group_cobjs_level, float voxel_spacing, float scale, bool has_layer, model3d_xform_t const &model_xf) {

	// group_cobjs_level: 0=no grouping, 1=simple grouping, 2=vbo grouping, 3=full 3d model, 4=no cobjs, 5=cubes from polygons (voxels), 6=cubes from edges
	bool const group_cobjs(group_cobjs_level >= 1);
	bool const use_vbo    (group_cobjs_level == 2);
	bool const use_model3d(group_cobjs_level >= 3);
	bool const no_cobjs   (group_cobjs_level >= 4);
	bool const use_cubes  (group_cobjs_level >= 5);
	bool const cube_edges (group_cobjs_level >= 6);
	
	if (!no_cobjs) { // add cobjs for collision detection
		add_model_polygons_to_cobjs(ppts, cobj, group_cobjs, use_vbo, (use_model3d ? COBJ_TYPE_MODEL3D : COBJ_TYPE_STD), has_layer, scale);
	}
	else if (use_cubes) {
		voxel_spacing *= scale;
		if (voxel_spacing <= 0.0) {return "model file voxel_spacing (scaled)";}
		vector<cube_t> cubes;
		if (cube_edges) {get_cur_model_edges_as_cubes(cubes, model_xf, voxel_spacing);}
		else {get_cur_model_as_cubes(cubes, model_xf, voxel_spacing);}
		check_layer(has_layer);
		coll_obj cur_cube(cobj); // color and tid left as-is for now
		cur_cube.cp.draw      = 0;
		cur_cube.cp.surfs     = 0; // clear
		cur_cube.type         = COLL_CUBE;
		cur_cube.cp.cobj_type = COBJ_TYPE_MODEL3D;
		//cur_cube.cp.color = BLUE; cur_cube.cp.draw = 1; // testing
		maybe_reserve_fixed_cobjs(cubes.size());

		for (vector<cube_t>::const_iterator i = cubes.begin(); i != cubes.end(); ++i) {
			cur_cube.copy_from(*i);
			cur_cube.id = (int)fixed_cobjs.size();
			fixed_cobjs.push_back(cur_cube);
		}
	}
	else if (use_model3d) {
		using_model_bcube = 1; // if any model is set (group_cobjs_level should generally agree across all models loaded)
	}
	return "";
}


// returns the number of values read
unsigned read_cube(FILE *fp, geom_xform_t const &xf, cube_t &c) {

	point pt[2];

	for (unsigned d = 0; d < 6; ++d) {
		if (fscanf(fp, "%f", &pt[d&1][d>>1]) != 1) return d;
	}
	for (unsigned i = 0; i < 2; ++i) {xf.xform_pos(pt[i]);}
	c = cube_t(pt[0], pt[1]);
	return 6;
}


bool is_end_of_string(int v) {return (v == '#' || v == EOF || v == 0 || isspace(v));}


int read_coll_obj_file(const char *coll_obj_file, geom_xform_t xf, coll_obj cobj, bool has_layer, colorRGBA lcolor) {

	assert(coll_obj_file != NULL);
	FILE *fp;
	if (!open_file(fp, coll_obj_file, "collision object")) return 0;
	char str[MAX_CHARS];
	unsigned line_num(1), npoints(0), indir_dlight_ix(0), prev_light_ix_start(0);
	int end(0), use_z(0), use_vel(0), ivals[3];
	float fvals[3];
	point pos(all_zeros);
	vector3d tv0, vel;
	polygon_t poly;
	vector<coll_tquad> ppts;
	// tree state
	float tree_br_scale_mult(1.0), tree_nl_scale(1.0), tree_height(1.0);
	bool enable_leaf_wind(1), remove_t_junctions(0);
	int reflective(0); // reflective: 0=none, 1=planar, 2=cube map (applies to cobjs and model3d)
	typedef map<string, cobj_params> material_map_t;
	material_map_t materials;
	multi_trigger_t triggers;
	
	while (!end) { // available: hou
		assert(fp != NULL);
		int letter(getc(fp));

		if (!is_end_of_string(letter)) {
			int const next_letter(getc(fp));

			if (!is_end_of_string(next_letter)) {
				string keyword;
				keyword.push_back(letter);
				letter = next_letter;
				while (!is_end_of_string(letter)) {keyword.push_back(letter); letter = getc(fp);}

				if (0) {}
				else if (keyword == "cube") {letter = 'B';}
				else if (keyword == "sphere") {letter = 'S';}
				else if (keyword == "cylinder") {letter = 'C';}
				else if (keyword == "capsule") {letter = 'k';}
				else if (keyword == "polygon") {letter = 'P';}
				else if (keyword == "torus") {letter = 'z';}
				else if (keyword == "trigger") {letter = 'K';}
				else if (keyword == "light"  ) {letter = 'L';}
				else if (keyword == "bind_light") {letter = 'V';}
				else if (keyword == "indir_dlight_group") {letter = 'U';}
				else if (keyword == "movable") {letter = 'd';}
				else if (keyword == "end") {letter = 'q';}
				else if (keyword == "density") {
					if (fscanf(fp, "%f", &cobj.cp.density) != 1) {return read_error(fp, "density", coll_obj_file);}
				}
				else if (keyword == "tj") {
					if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "remove t junctions", coll_obj_file);}
					remove_t_junctions = (ivals[0] != 0);
				}
				else if (keyword == "reflective") {
					if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "reflective", coll_obj_file);}
					reflective = ((ivals[0] != 0) ? 1 : 0);
					cobj.set_reflective_flag(reflective == 2); // only for cube maps
				}
				else if (keyword == "cube_map_ref") {
					if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "cube_map_ref", coll_obj_file);}
					reflective = ((ivals[0] != 0) ? 2 : 0);
					cobj.set_reflective_flag(reflective == 2); // only for cube maps
				}
				else if (keyword == "metalness") {
					if (fscanf(fp, "%f", &cobj.cp.metalness) != 1) {return read_error(fp, "metalness", coll_obj_file);}
				}
				else if (keyword == "start_cobj_group") {cobj.cgroup_id = cobj_groups.new_group();}
				else if (keyword == "end_cobj_group") {cobj.cgroup_id = -1;}
				else if (keyword == "start_draw_group") {cobj.dgroup_id = cdraw_groups.new_group();}
				else if (keyword == "end_draw_group") {cobj.dgroup_id = -1;}
				else if (keyword == "destroy_prob") {
					if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "destroy_prob", coll_obj_file);}
					cobj.cp.destroy_prob = (unsigned char)max(0, min(255, ivals[0])); // 0 = default
				}
				else {
					ostringstream oss;
					oss << "unrecognized keyword: '" << keyword << "' on line " << line_num;
					return read_error(fp, oss.str().c_str(), coll_obj_file);
				}
			}
		}
		switch (letter) {
		case 0:
		case EOF:
			end = 1;
			break;

		case '\n':
			++line_num;
		case '\t': case '\f': case '\r': case '\v': case ' ':
			break;

		case '#': // line comment
			do {letter = getc(fp);} while (letter != '\n' && letter != EOF && letter != 0);
			if (letter == '\n') {++line_num;}
			break;

		case 'i': // include file (only translation and scale are saved state)
			{
				string const fn(read_filename(fp, line_num));
				if (fn.empty()) {return read_error(fp, "include file", coll_obj_file);}
				if (!read_coll_obj_file(fn.c_str(), xf, cobj, has_layer, lcolor)) {return read_error(fp, "include file", coll_obj_file);}
			}
			break;

		case 'O': // load *.obj | *.3ds | *.model3d file: <filename> <group_cobjs_level> <recalc_normals/use_vertex_normals> <write_file> [<voxel_xy_spacing>]
			{
				string const fn(read_filename(fp, line_num));
				int recalc_normals(0), write_file(0);
				float voxel_xy_spacing(0.0);

				if (fn.empty() || fscanf(fp, "%i%i%i%f", &ivals[0], &recalc_normals, &write_file, &voxel_xy_spacing) < 3) {
					return read_error(fp, "load model file command", coll_obj_file);
				}
				if (ivals[0] < 0 || ivals[0] > 6) {return read_error(fp, "load model file command group_cobjs_level", coll_obj_file);}
				bool const use_model3d(ivals[0] >= 3), no_cobjs(ivals[0] >= 4);
				ppts.clear();
				RESET_TIME;
				
				if (!read_model_file(fn, (no_cobjs ? nullptr : &ppts), xf, cobj.cp.tid, cobj.cp.color, reflective, cobj.cp.metalness, use_model3d, (recalc_normals != 0), (write_file != 0), 1)) {
					return read_error(fp, "model file data", coll_obj_file);
				}
				string const error_str(add_loaded_model(ppts, cobj, ivals[0], voxel_xy_spacing, xf.scale, has_layer, model3d_xform_t()));
				if (!error_str.empty()) {return read_error(fp, error_str.c_str(), coll_obj_file);}
				PRINT_TIME("Model File Load/Process");
				break;
			}

		case 'Z': // add model3d transform: group_cobjs_level tx ty tz [scale [rx ry rz angle [<voxel_spacing>]]]
			{
				model3d_xform_t model_xf;
				float voxel_spacing(0.0);
				int const num_args(fscanf(fp, "%i%f%f%f%f%f%f%f%f%f", &ivals[0], &model_xf.tv.x, &model_xf.tv.y, &model_xf.tv.z, &model_xf.scale,
					&model_xf.axis.x, &model_xf.axis.y, &model_xf.axis.z, &model_xf.angle, &voxel_spacing));
				if (num_args != 4 && num_args != 5 && num_args != 9 && num_args != 10) {return read_error(fp, "model3d transform", coll_obj_file);}
				if (ivals[0] < 0 || ivals[0] > 6) {return read_error(fp, "add model transform command group_cobjs_level", coll_obj_file);}
				if (model_xf.scale == 0.0) {return read_error(fp, "model3d transform scale", coll_obj_file);} // what about negative scales?
				model_xf.material = cobj.cp; // copy base material from cobj
				add_transform_for_cur_model(model_xf);
				bool const no_cobjs(ivals[0] >= 4);
				if (!no_cobjs) {get_cur_model_polygons(ppts, model_xf);} // add cobjs for collision detection
				string const error_str(add_loaded_model(ppts, cobj, ivals[0], voxel_spacing, xf.scale, has_layer, model_xf));
				if (!error_str.empty()) {return read_error(fp, error_str.c_str(), coll_obj_file);}
			}
			break;

		case 'Q': // platform: enabled [fspeed rspeed sdelay rdelay ext_dist|rot_angle act_dist origin<x,y,z> dir|rot_axis<x,y,z> cont [is_rotation=0]]
			if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "platform", coll_obj_file);}
			assert(ivals[0] == 0 || ivals[0] == 1); // boolean
			
			if (ivals[0] == 0) { // disable platforms
				cobj.platform_id = -1;
			}
			else {
				cobj.platform_id = (short)platforms.size();
				if (!platforms.add_from_file(fp, xf, triggers)) {return read_error(fp, "platform", coll_obj_file);}
				assert(cobj.platform_id < (int)platforms.size());
			}
			break;

		case 'g': // set tree parameters state: height, branch scale, nleaves scale, enable wind
			if (fscanf(fp, "%f%f%f%i", &tree_height, &tree_br_scale_mult, &tree_nl_scale, &ivals[0]) != 4 || tree_height <= 0.0 || tree_br_scale_mult <= 0.0 || tree_nl_scale <= 0.0) {
				return read_error(fp, "tree_params", coll_obj_file);
			}
			enable_leaf_wind = (ivals[0] != 0);
			break;

		case 'E': // place tree: xpos ypos size type [zpos [tree_4th_branches]], type: TREE_MAPLE = 0, TREE_LIVE_OAK = 1, TREE_A = 2, TREE_B = 3, 4 = TREE_PAPAYA
			if (fscanf(fp, "%f%f%f%i", &pos.x, &pos.y, &fvals[0], &ivals[0]) != 4) {return read_error(fp, "tree", coll_obj_file);}
			assert(fvals[0] > 0.0);
			use_z = (fscanf(fp, "%f", &pos.z) == 1);
			
			if (num_trees > 0) {
				cout << "Must set ntrees to zero in order to add trees through collision objects file." << endl;
			}
			else {
				bool local_tree_4th_branches(tree_4th_branches);
				if (fscanf(fp, "%i", &ivals[1])) {local_tree_4th_branches = (ivals[1] != 0);}
				xf.xform_pos(pos);
				t_trees.push_back(tree(enable_leaf_wind));
				t_trees.back().gen_tree(pos, max(1, int(fvals[0]*xf.scale)), ivals[0], !use_z, 0, 1, tree_height, tree_br_scale_mult, tree_nl_scale, local_tree_4th_branches);
				tree_mode |= 1; // enable trees
			}
			break;

		case 'H': // place hedges: xstart ystart dx dy nsteps size, type [cx1 cx2 cy1 cy2 cz1 cz2]
			if (fscanf(fp, "%f%f%f%f%i%f%i", &pos.x, &pos.y, &fvals[0], &fvals[1], &ivals[0], &fvals[2], &ivals[1]) != 7 || ivals[0] <= 0) {
				return read_error(fp, "hedges", coll_obj_file);
			}
			if (num_trees > 0) {
				cout << "Must set ntrees to zero in order to add hedges through collision objects file." << endl;
			}
			else {
				cube_t clip_cube;
				unsigned const num_read(read_cube(fp, xf, clip_cube));
				bool const use_clip_cube(num_read == 6);
				if (num_read > 0 && !use_clip_cube) {return read_error(fp, "hedges bounding cube", coll_obj_file);}
				xf.xform_pos(pos);
				point cur(pos);
				vector3d delta(fvals[0], fvals[1], 0.0); // single hedge with zero delta is okay - this is just a tree with a clip cube
				xf.xform_pos_rms(delta);

				for (int i = 0; i < ivals[0]; ++i) {
					t_trees.push_back(tree(enable_leaf_wind));
					if (use_clip_cube) {t_trees.back().enable_clip_cube(clip_cube);}
					t_trees.back().gen_tree(cur, max(1, int(fvals[2]*xf.scale)), ivals[1], 1, 0, 1, tree_height, tree_br_scale_mult, tree_nl_scale, 0); // no 4th branches
					cur += delta;
				}
				tree_mode |= 1; // enable trees
			}
			break;

		case 'F': // place small tree: xpos ypos height width type [zpos], type: T_PINE = 0, T_DECID = 1, T_TDECID = 2, T_BUSH = 3, T_PALM = 4, T_SH_PINE = 5
			if (fscanf(fp, "%f%f%f%f%i", &pos.x, &pos.y, &fvals[0], &fvals[1], &ivals[0]) != 5) {return read_error(fp, "small tree", coll_obj_file);}
			assert(fvals[0] > 0.0 && fvals[1] > 0.0);
			use_z = (fscanf(fp, "%f", &pos.z) == 1);
			xf.xform_pos(pos);

			if (add_small_tree(pos, xf.scale*fvals[0], xf.scale*fvals[1], ivals[0], !use_z)) {
				tree_mode |= 2; // enable small trees
			}
			else {
				cout << "Error: Failed to add a small tree of type " << ivals[0] << endl;
			}
			break;

		case 'G': // place plant: xpos ypos height radius type [zpos], type: PLANT_MJ = 0, PLANT1, PLANT2, PLANT3, PLANT4
			if (fscanf(fp, "%f%f%f%f%i", &pos.x, &pos.y, &fvals[0], &fvals[1], &ivals[0]) != 5) {return read_error(fp, "plant", coll_obj_file);}
			assert(fvals[0] > 0.0 && fvals[1] > 0.0);
			use_z = (fscanf(fp, "%f", &pos.z) == 1);
			xf.xform_pos(pos);
			add_plant(pos, xf.scale*fvals[0], xf.scale*fvals[1], ivals[0], !use_z);
			break;

		case 'A': // appearance spot: xpos ypos [zpos]
			if (fscanf(fp, "%f%f", &pos.x, &pos.y) != 2) {return read_error(fp, "appearance spot", coll_obj_file);}
			{
				float const smiley_radius(object_types[SMILEY].radius);
				read_or_calc_zval(fp, pos, smiley_radius, smiley_radius, xf);
				app_spots.push_back(pos);
			}
			break;

		case 'L': // point/spot/line light: ambient_size diffuse_size xpos ypos zpos color [direction|pos2 [beamwidth=1.0 [inner_radius=0.0 [is_line_light=0 [use_shadow_map=0]]]]]
			// type: 0 = ambient/baked only, 1 = diffuse/dynamic only, 2 = both
			if (fscanf(fp, "%f%f%f%f%f%f%f%f%f", &fvals[0], &fvals[1], &pos.x, &pos.y, &pos.z, &lcolor.R, &lcolor.G, &lcolor.B, &lcolor.A) != 9) {
				return read_error(fp, "light source", coll_obj_file);
			}
			{
				prev_light_ix_start = light_sources_d.size(); // capture start of range so that we can allow binding of grouped lights (from cube map)
				xf.xform_pos(pos);
				float beamwidth(1.0), r_inner(0.0);
				vector3d dir(zero_vector); // defaults to (0,0,0)
				ivals[0] = 0; // default is point/spotlight
				point pos2(pos);
				int use_smap(0); // 0=none, 1=spotlight, 2=point light/cube map

				// direction|pos2 [beamwidth=1.0 [inner_radius=0.0 [is_line_light=0 [use_shadow_map=0]]]]
				if (fscanf(fp, "%f%f%f%f%f%i%i", &dir.x, &dir.y, &dir.z, &beamwidth, &r_inner, &ivals[0], &use_smap) >= 3) {
					if (use_smap < 0 || use_smap > 2) {return read_error(fp, "light source use_smap (must be 0, 1, or 2)", coll_obj_file);}
					if (ivals[0] != 0) {pos2 = dir; dir = zero_vector; beamwidth = 1.0; xf.xform_pos(pos2);} // line light
					else {xf.xform_pos_rm(dir);} // spotlight (or hemispherical light ray culling if beamwidth == 1.0)
				}
				for (unsigned d = 0; d < 2; ++d) { // {ambient, diffuse}
					if (fvals[d] == 0.0) continue;
					vector<light_source> lss;

					if (use_smap == 2 && ivals[0] == 0 && beamwidth == 1.0) { // special case for shadowed point light with cube map
						if (dir != zero_vector) {return read_error(fp, "light source point light smap (dir must be all zeros)", coll_obj_file);}
						beamwidth = 0.4; // 0.3 to 0.5 are okay

						for (unsigned ldim = 0; ldim < 3; ++ldim) {
							dir = zero_vector;
							for (unsigned ldir = 0; ldir < 2; ++ldir) {
								dir[ldim] = (ldir ? 1.0 : -1.0);
								lss.push_back(light_source(fvals[d], pos, pos2, lcolor, 0, dir, beamwidth, r_inner, (use_smap == 2))); // add lights for each cube face
							}
						}
					}
					else {
						lss.push_back(light_source(fvals[d], pos, pos2, lcolor, 0, dir, beamwidth, r_inner, (use_smap == 2)));
					}
					for (auto ls = lss.begin(); ls != lss.end(); ++ls) {
						if (d) { // diffuse
							if (cobj.platform_id >= 0) {platforms.get_cobj_platform(cobj).add_light(light_sources_d.size());}
							indir_dlight_group_manager.add_dlight_ix_for_tag_ix(indir_dlight_ix, light_sources_d.size());
							light_sources_d.push_back(light_source_trig(*ls, (use_smap != 0), cobj.platform_id, indir_dlight_ix));
							light_sources_d.back().add_triggers(triggers);
						}
						else {light_sources_a.push_back(*ls);} // ambient
					}
				}
			}
			break;

			// FIXME: should region option use first {x,y,z} as the region LLC and only specify the region URC?
		case 'K': // scene diffuse point light or platform trigger: x y z  activate_dist auto_on_time auto_off_time player_only requires_action [act_cube_region x1 x2 y1 y2 z1 z2]
			{
				trigger_t trigger;
				int const num_read(fscanf(fp, "%f%f%f%f%f%f%i%i", &trigger.act_pos.x, &trigger.act_pos.y, &trigger.act_pos.z,
					&trigger.act_dist, &trigger.auto_on_time, &trigger.auto_off_time, &ivals[0], &ivals[1]));
				if (num_read <= 0) {triggers.clear(); break;} // bare K, just reset params and disable the trigger, or EOF
				if (num_read != 8) {return read_error(fp, "light source trigger", coll_obj_file, line_num);}
				xf.xform_pos(trigger.act_pos);
				trigger.act_dist       *= xf.scale;
				trigger.player_only     = (ivals[0] != 0);
				trigger.requires_action = (ivals[1] != 0);
				cube_t act_region;
				int const num_read2(read_cube(fp, xf, act_region));
				if (num_read2 == 6) {trigger.set_act_region(act_region);}
				else if (num_read2 > 0) {return read_error(fp, "light source trigger activation cube", coll_obj_file);}
				triggers.push_back(trigger);
			}
			break;

		case 'V': // bind prev light source to cobj at location <x y z>
			if (fscanf(fp, "%f%f%f", &pos.x, &pos.y, &pos.z) != 3) {return read_error(fp, "light source", coll_obj_file);}
			if (prev_light_ix_start == 0 || prev_light_ix_start == light_sources_d.size()) {return read_error(fp, "light source <no previous dynamic light source to bind to>", coll_obj_file);}
			assert(prev_light_ix_start > FLASHLIGHT_LIGHT_ID);
			xf.xform_pos(pos);
			for (unsigned i = prev_light_ix_start; i < light_sources_d.size(); ++i) {light_sources_d[i].bind_to_pos(pos);} // bind the group of light sources
			break;

		case 'b': // cube volume light (for sky/global indirect): x1 x2 y1 y2 z1 z2  color.R color.G color.B  intensity num_rays ltype [disabled_edges_bits]
			{
				cube_light_src cls;
				if (read_cube(fp, xf, cls.bounds) != 6) {return read_error(fp, "cube volume global light", coll_obj_file);}
				
				if (fscanf(fp, "%f%f%f%f%u%i%u", &cls.color.R, &cls.color.G, &cls.color.B, &cls.intensity, &cls.num_rays, &ivals[0], &cls.disabled_edges) < 6) {
					return read_error(fp, "cube volume global light", coll_obj_file);
				}
				switch (ivals[0]) {
				case LIGHTING_SKY:    sky_cube_lights.push_back   (cls); break;
				case LIGHTING_GLOBAL: global_cube_lights.push_back(cls); break;
				default: return read_error(fp, "cube volume global light ltype param", coll_obj_file);
				}
			}
			break;

		case 'U': // indir dlight group: name [scale]
			fvals[0] = 1.0; // default scale
			if (fscanf(fp, "%255s%f", str, &fvals[0]) < 1) {return read_error(fp, "indir dlight group name", coll_obj_file);}
			indir_dlight_ix = indir_dlight_group_manager.get_ix_for_name(str, fvals[0]);
			break;

		case 'f': // place fire: size light_beamwidth intensity xpos ypos zpos
			if (fscanf(fp, "%f%f%f%f%f%f", &fvals[0], &fvals[1], &fvals[2], &pos.x, &pos.y, &pos.z) != 6) {
				return read_error(fp, "place fire", coll_obj_file);
			}
			xf.xform_pos(pos);
			gen_fire(pos, fvals[0], NO_SOURCE, 1, 1, fvals[1], fvals[2]);
			break;

		case 'p': // smiley path waypoint: type xpos ypos [zpos]
			if (fscanf(fp, "%i%f%f", &ivals[0], &pos.x, &pos.y) != 3) { // type: 0 = normal, 1 = goal
				return read_error(fp, "waypoint", coll_obj_file);
			}
			{
				read_or_calc_zval(fp, pos, SMALL_NUMBER, object_types[WAYPOINT].radius, xf);
				user_waypoints.push_back(user_waypt_t(ivals[0], pos));
			}
			break;

		case 'I': // items (health=28, shield=29, powerup=30, weapon=31, ammo=32)
			// obj_class obj_subtype regen_time(s) xpos ypos [zpos]
			if (fscanf(fp, "%i%i%f%f%f", &ivals[0], &ivals[1], &fvals[0], &pos.x, &pos.y) != 5) {
				return read_error(fp, "place item", coll_obj_file);
			}
			{
				assert(ivals[0] >= 0 && ivals[0] < NUM_TOT_OBJS);
				float const radius(object_types[ivals[0]].radius);
				read_or_calc_zval(fp, pos, radius, radius, xf);
				init_objects();
				create_object_groups();
				int const cid(coll_id[ivals[0]]);
				assert(cid < NUM_TOT_OBJS);
				assert(obj_groups[cid].type == ivals[0]); // make sure this group is properly initialized
				obj_groups[cid].add_predef_obj(pos, ivals[1], fvals[0]*TICKS_PER_SECOND);
			}
			break;

		case 'w': // water spring/source: xpos ypos rate [zpos] [vx vy vz diff]
			if (fscanf(fp, "%f%f%f", &pos.x, &pos.y, &fvals[0]) != 3) {return read_error(fp, "water source", coll_obj_file);}
			fvals[1] = 0.1;
			use_z    = (fscanf(fp, "%f", &pos.z) == 1);
			use_vel  = (fscanf(fp, "%f%f%f%f", &vel.x, &vel.y, &vel.z, &fvals[1]) == 4);
			if (use_vel) xf.xform_pos_rms(vel); // scale?
			xf.xform_pos(pos);
			add_water_spring(pos, vel, fvals[0], fvals[1], !use_z, !use_vel);
			break;

		case 'W': // water section
			{
				float x1, y1, x2, y2, zval, wvol;

				if (fscanf(fp, "%f%f%f%f%f%f", &x1, &x2, &y1, &y2, &zval, &wvol) != 6) {
					return read_error(fp, "water section", coll_obj_file);
				}
				add_water_section((xf.scale*x1+xf.tv[0]), (xf.scale*y1+xf.tv[1]), (xf.scale*x2+xf.tv[0]), (xf.scale*y2+xf.tv[1]), (xf.scale*zval+xf.tv[2]), wvol);
			}
			break;

		case 'e': // set shape edges to skip when drawing (cube and cylinder)
			if (fscanf(fp, "%i", &ivals[0]) != 1 || ivals[0] < 0 || ivals[0] > 255) {
				return read_error(fp, "shape edge skip draw", coll_obj_file);
			}
			cobj.cp.surfs = (unsigned char)ivals[0]; // all six sides (cube) or top/bottom (cylinder)
			break;

		case 'B': // cube: xmin xmax ymin ymax zmin zmax [corner_radius]
			{
				if (read_cube(fp, xf, cobj) != 6) {return read_error(fp, "collision cube", coll_obj_file);}
				if (cobj.is_zero_area()) {return read_error(fp, "collision cube: zero area cube", coll_obj_file);}
				float corner_radius(0.0);
				fscanf(fp, "%f", &corner_radius); // okay if fails
				check_layer(has_layer);
				cobj.radius2 = corner_radius*xf.scale;
				cobj.counter = (remove_t_junctions ? OBJ_CNT_REM_TJ : 0); // remove T-junctions
				cobj.add_to_vector(fixed_cobjs, COLL_CUBE);
				cobj.radius2 = 0.0; // reset
				cobj.counter = 0; // reset
				if (cobj.platform_id >= 0 && platforms[cobj.platform_id].is_rotation()) {return read_error(fp, "collision cube: cannot use a rotation platform", coll_obj_file);}
			}
			break;

		case 'S': // sphere: x y z radius
			if (fscanf(fp, "%f%f%f%f", &cobj.points[0].x, &cobj.points[0].y, &cobj.points[0].z, &cobj.radius) != 4) {
				return read_error(fp, "collision sphere", coll_obj_file);
			}
			check_layer(has_layer);
			cobj.radius *= xf.scale;
			xf.xform_pos(cobj.points[0]);
			cobj.add_to_vector(fixed_cobjs, COLL_SPHERE);
			break;

		case 'C': // cylinder: x1 y1 z1 x2 y2 z2 r1 r2
		case 'k': // capsule: x1 y1 z1 x2 y2 z2 r1 r2
			if (fscanf(fp, "%f%f%f%f%f%f%f%f", &cobj.points[0].x, &cobj.points[0].y, &cobj.points[0].z, &cobj.points[1].x, &cobj.points[1].y, &cobj.points[1].z, &cobj.radius, &cobj.radius2) != 8) {
				return read_error(fp, "collision cylinder/capsule", coll_obj_file);
			}
			assert(cobj.radius >  0.0 || cobj.radius2 >  0.0);
			assert(cobj.radius >= 0.0 && cobj.radius2 >= 0.0);
			check_layer(has_layer);
			cobj.radius  *= xf.scale;
			cobj.radius2 *= xf.scale;
			for (unsigned i = 0; i < 2; ++i) {xf.xform_pos(cobj.points[i]);}
			// cylinder surfs: 0 = draw ends + bfc (solid), 1 = no draw ends + no bfc (hollow), 3 = no draw ends + bfc (ends are hidden)
			cobj.add_to_vector(fixed_cobjs, ((letter == 'k') ? COLL_CAPSULE : COLL_CYLINDER));
			break;

		case 'z': // torus: x y z dir_x dir_y dir_z ro ri
			if (fscanf(fp, "%f%f%f%f%f%f%f%f", &cobj.points[0].x, &cobj.points[0].y, &cobj.points[0].z, &cobj.norm.x, &cobj.norm.y, &cobj.norm.z, &cobj.radius, &cobj.radius2) != 8) {
				return read_error(fp, "collision torus", coll_obj_file);
			}
			assert(cobj.radius > 0.0 && cobj.radius2 > 0.0);
			check_layer(has_layer);
			cobj.radius  *= xf.scale;
			cobj.radius2 *= xf.scale;
			xf.xform_pos(cobj.points[0]);
			xf.xform_pos_rm(cobj.norm);
			cobj.norm.normalize();
			cobj.add_to_vector(fixed_cobjs, COLL_TORUS);
			break;

		case 'P': // polygon: npts (x y z)* thickness
			if (fscanf(fp, "%u", &npoints) != 1) {return read_error(fp, "collision polygon npoints", coll_obj_file);}

			if (npoints < 3) {
				cout << "Error: Collision polygon must have at least 3 points: " << npoints << "." << endl;
				fclose(fp);
				return 0;
			}
			cobj.npoints = npoints;
			poly.resize(npoints);

			for (unsigned i = 0; i < npoints; ++i) {
				if (fscanf(fp, "%f%f%f", &poly[i].v.x, &poly[i].v.y, &poly[i].v.z) != 3) {
					cout << "Error reading collision polygon point " << i << " from file '" << coll_obj_file << "'." << endl;
					fclose(fp);
					return 0;
				}
				xf.xform_pos(poly[i].v);
			}
			if (fscanf(fp, "%f", &cobj.thickness) != 1) {return read_error(fp, "collision polygon", coll_obj_file);}
			check_layer(has_layer);
			cobj.thickness *= xf.scale;
			ppts.resize(0);
			split_polygon(poly, ppts, 0.99);
			add_polygons_to_cobj_vector(ppts, cobj, NULL, cobj.cp.cobj_type);
			break;

		case 'c': // hollow cylinder (multisided): x1 y1 z1  x2 y2 z2  ro ri  nsides [start_ix [end_ix]]
			{
				point pt[2];
				float ro, ri;

				if (fscanf(fp, "%f%f%f%f%f%f%f%f%u", &pt[0].x, &pt[0].y, &pt[0].z, &pt[1].x, &pt[1].y, &pt[1].z, &ro, &ri, &npoints) != 9) {
					return read_error(fp, "hollow cylinder", coll_obj_file);
				}
				if (npoints < 3 || ro <= 0.0 || ri < 0.0 || ro < ri || pt[0] == pt[1]) {
					return read_error(fp, "hollow cylinder values", coll_obj_file);
				}
				unsigned six(0), eix(npoints);
				fscanf(fp, "%u%u", &six, &eix); // optional
				if (six >= eix || eix > npoints) {return read_error(fp, "hollow cylinder start/end indices", coll_obj_file);}
				check_layer(has_layer);
				for (unsigned i = 0; i < 2; ++i) {xf.xform_pos(pt[i]);}
				cobj.thickness = xf.scale*(ro - ri);
				float const r(0.5*xf.scale*(ro + ri)), step(TWO_PI/float(npoints)), edist(0.5*cobj.thickness*tanf(0.5*step));
				vector3d const vc((pt[1] - pt[0]).get_norm());
				unsigned const dmin((vc.x < vc.y) ? ((vc.x < vc.z) ? 0 : 2) : ((vc.y < vc.z) ? 1 : 2));
				vector3d vn(zero_vector), dirs[2];
				vn[dmin] = 1.0;
				cross_product(vc, vn,      dirs[0]); // first ortho dir
				cross_product(vc, dirs[0], dirs[1]); // second ortho dir
				for (unsigned i = 0; i < 2; ++i) dirs[i].normalize();

				for (unsigned i = six; i < eix; ++i) {
					float const val[2] = {(i - 0.5), (i + 0.5)};
					vector3d deltas[2];

					for (unsigned j = 0; j < 2; ++j) {
						float const v(step*val[j]);
						deltas[j] = (dirs[0]*cosf(v) + dirs[1]*sinf(v))*r;
					}
					vector3d const extend((deltas[1] - deltas[0]).get_norm()*edist);
					deltas[0] -= extend;
					deltas[1] += extend;
					for (unsigned j = 0; j < 4; ++j) {cobj.points[j] = pt[j>>1] + deltas[(j>>1)^(j&1)];}
					cobj.npoints = 4; // have to reset every time in case it was a cube
					cobj.add_to_vector(fixed_cobjs, COLL_POLYGON);
				}
			}
			break;

		case 'N': // portal: xyz1 xyz2 xyz3 xyz4 [nx ny nz]
			{
				portal p;
				for (unsigned i = 0; i < 4; ++i) {
					if (fscanf(fp, "%f%f%f", &p.pts[i].x, &p.pts[i].y, &p.pts[i].z) != 3) {return read_error(fp, "portal", coll_obj_file);}
					xf.xform_pos(p.pts[i]);
				}
				fscanf(fp, "%f%f%f", &p.normal.x, &p.normal.y, &p.normal.z); // optional/can fail
				portals.push_back(p);
			}
			break;

		case 'x': // teleporter sx sy sz  dx dy dz  radius
			{
				teleporter tp;
				if (fscanf(fp, "%f%f%f%f%f%f%f", &tp.pos.x, &tp.pos.y, &tp.pos.z, &tp.dest.x, &tp.dest.y, &tp.dest.z, &tp.radius) != 7) {
					return read_error(fp, "teleporter", coll_obj_file);
				}
				xf.xform_pos(tp.pos);
				xf.xform_pos(tp.dest);
				tp.setup();
				teleporters.push_back(tp);
			}
			break;

		case 'D': // step delta (for stairs, etc.): dx dy dz num [dsx [dsy [dsz]]]
			if (cobj.type == COLL_NULL) {return read_error(fp, "step delta must appear after shape definition", coll_obj_file);}
			if (fscanf(fp, "%f%f%f%u", &pos.x, &pos.y, &pos.z, &npoints) != 4) {return read_error(fp, "step delta", coll_obj_file);}
			vel = zero_vector; // size delta
			fscanf(fp, "%f%f%f", &vel.x, &vel.y, &vel.z);
			if (pos == all_zeros && vel == zero_vector) {return read_error(fp, "step delta must have nonzero delta", coll_obj_file);}
			xf.xform_pos_rms(pos); // no translate
			xf.xform_pos_rms(vel); // no translate

			for (unsigned i = 0; i < npoints; ++i) {
				if (vel != zero_vector) {
					switch (cobj.type) {
					case COLL_CUBE:
						for (unsigned j = 0; j < 3; ++j) {cobj.d[j][1] += vel[j];}
						cobj.normalize();
						break;
					case COLL_CYLINDER:
					case COLL_CYLINDER_ROT:
					case COLL_CAPSULE:
						cobj.points[1] += vel; // just increase the length?
						break;
					case COLL_POLYGON: break; // nothing to translate by
					case COLL_SPHERE:  break; // nothing to translate by
					default: assert(0);
					}
				}
				cobj.shift_by(pos, 1, 1);
				cobj.add_to_vector(fixed_cobjs, cobj.type);
			}
			break;

		case 'l': // object layer/material: elasticity R G B A texture_id/texture_name [draw=1 [refract_ix=1.0 [light_atten=0.0 [emissive=0]]]]
			if (fscanf(fp, "%f%f%f%f%f%255s", &cobj.cp.elastic, &cobj.cp.color.R, &cobj.cp.color.G, &cobj.cp.color.B, &cobj.cp.color.A, str) != 6) {
				return read_error(fp, "layer/material properties", coll_obj_file);
			}
			if (!read_texture(str, line_num, cobj.cp.tid, 0)) {fclose(fp); return 0;}
			has_layer    = 1;
			cobj.cp.draw = (fscanf(fp, "%i", &ivals[0]) ? (ivals[0] != 0) : 1); // optional
			if (!fscanf(fp, "%f", &cobj.cp.refract_ix )) {cobj.cp.refract_ix  = 1.0;} // optional
			if (!fscanf(fp, "%f", &cobj.cp.light_atten)) {cobj.cp.light_atten = 0.0;} // optional
			if (!fscanf(fp, "%i", &ivals[0])) {ivals[0] = 0;} // optional
			cobj.cp.is_emissive = (ivals[0] != 0);
			break;

		case 'j': // restore material <name>
			if (fscanf(fp, "%255s", str) != 1) {return read_error(fp, "restore material name", coll_obj_file);}
			{
				material_map_t::const_iterator it(materials.find(str));
				if (it == materials.end()) {
					cout << "*** Error: material '" << str << "' not defined in file '" << coll_obj_file << "'. ***" << endl;
					fclose(fp);
					return 0;
				}
				cobj.cp = it->second;
			}
			break;
		case 'J': // save material <name>
			if (fscanf(fp, "%255s", str) != 1) {return read_error(fp, "save material name", coll_obj_file);}
			materials[str] = cobj.cp; // Note: okay to overwrite/redefine a material
			break;

		case 'X': // normal map texture id/name [invert_y=0 [swap_binorm_sign=0]]
			if (fscanf(fp, "%255s", str) != 1) {return read_error(fp, "normal map texture", coll_obj_file);}
			{
				int invert_y(0), swap_bns(0);
				fscanf(fp, "%i%i", &invert_y, &swap_bns); // optional
				cobj.cp.set_swap_tcs_flag(SWAP_TCS_NM_BS, (swap_bns != 0));
				if (!read_texture(str, line_num, cobj.cp.normal_map, 1, (invert_y != 0))) {fclose(fp); return 0;}
			}
			break;

		case 'r': // set specular: <specular intensity> <shininess> [R G B]
			{
				if (fscanf(fp, "%f%f", &fvals[0], &cobj.cp.shine) != 2) {return read_error(fp, "specular lighting", coll_obj_file);}
				int const num_read(fscanf(fp, "%f%f%f", &cobj.cp.spec_color.R, &cobj.cp.spec_color.G, &cobj.cp.spec_color.B));
				if (num_read > 0) {
					if (num_read != 3) {return read_error(fp, "specular lighting {R,G,B} values", coll_obj_file);}
					cobj.cp.spec_color *= fvals[0]; // multiply by intensity
				}
				else {
					cobj.cp.spec_color.set_to_val(fvals[0]);
				}
			}
			break;

		case 't': // relative translate
			if (fscanf(fp, "%f%f%f", &tv0.x, &tv0.y, &tv0.z) != 3) {return read_error(fp, "translate", coll_obj_file);}
			xf.tv += tv0;
			break;
		case 'T': // absolute translate
			if (fscanf(fp, "%f%f%f", &xf.tv.x, &xf.tv.y, &xf.tv.z) != 3) {return read_error(fp, "translate", coll_obj_file);}
			break;

		case 'm': // scale/magnitude
			if (fscanf(fp, "%f", &xf.scale) != 1) {return read_error(fp, "scale", coll_obj_file);}
			assert(xf.scale > 0.0);
			break;

		case 'M': // mirror <dim>, dim = [0,1,2] => [x,y,z]
			if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "mirror", coll_obj_file);}
			if (ivals[0] < 0 || ivals[0] > 2) {return read_error(fp, "mirror: dim must be in [0,2]", coll_obj_file);}
			xf.mirror[ivals[0]] ^= 1;
			break;
		case 's': // swap dimensions <dim1> <dim2>
			if (fscanf(fp, "%i%i", &ivals[0], &ivals[1]) != 2) {return read_error(fp, "swap dimensions", coll_obj_file);}
			if (ivals[0] == ivals[1] || ivals[0] < 0 || ivals[0] > 2 || ivals[1] < 0 || ivals[1] > 2) {
				return read_error(fp, "swap dimensions: dims must be different and in [0,2]", coll_obj_file);
			}
			xf.swap_dim[ivals[0]][ivals[1]] ^= 1;
			break;
		case 'R': // restore mirrors and swaps to default
			xf.restore_mirror_and_swap();
			break;

		case 'y': // texture scale
			if (fscanf(fp, "%f", &cobj.cp.tscale) != 1) {return read_error(fp, "texture scale", coll_obj_file);} // okay if tscale=0 for cubes
			break;

		case 'Y': // texture translate (cubes, polygons, cylinder ends only), swap xy (cubes/polygons only): <tdx> <tdy> [<swap_xy>]
			if (fscanf(fp, "%f%f", &cobj.cp.tdx, &cobj.cp.tdy) != 2) {return read_error(fp, "texture translate", coll_obj_file);}
			ivals[0] = 0;
			fscanf(fp, "%i", &ivals[0]); // optional
			cobj.cp.set_swap_tcs_flag(SWAP_TCS_XY, (ivals[0] != 0));
			break;

		case 'n': // toggle negative shape
			if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "negative shape", coll_obj_file);}
			set_bit_flag_to(cobj.status, COLL_NEGATIVE, (ivals[0] != 0));
			break;
		case 'a': // set destroyability
			if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "destroy shape", coll_obj_file);}
			cobj.destroy = (EXPLODE_EVERYTHING ? EXPLODEABLE : (char)ivals[0]);
			break;
		case 'd': // toggle movable
			if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "movable", coll_obj_file);}
			set_bit_flag_to(cobj.cp.flags, COBJ_MOVABLE, (ivals[0] != 0));
			break;

		case 'v': // set voxel mode
			if (fscanf(fp, "%i", &ivals[0]) != 1) {return read_error(fp, "set voxel mode", coll_obj_file);}
			cobj.cp.cobj_type = (ivals[0] ? COBJ_TYPE_VOX_TERRAIN : COBJ_TYPE_STD);
			break;

		case 'q': // quit reading object file
			fclose(fp);
			fp  = NULL;
			end = 1;
			break;

		default:
			cerr << "Error: unrecognized symbol in collision object file line " << line_num << ": " << char(letter) << " (ASCII " << letter << ")." << endl;
			fclose(fp);
			exit(1);
			return 0;
		}
	}
	if (fp != NULL) fclose(fp);
	return 1;
}

int read_coll_objects(const char *coll_obj_file) {

	geom_xform_t xf;
	coll_obj cobj;
	cobj.init();
	cobj.cp.elastic = 0.5; // default
	cobj.cp.draw    = 1;   // default
	if (EXPLODE_EVERYTHING) {cobj.destroy = EXPLODEABLE;}
	if (use_voxel_cobjs) {cobj.cp.cobj_type = COBJ_TYPE_VOX_TERRAIN;}
	if (!read_coll_obj_file(coll_obj_file, xf, cobj, 0, WHITE)) return 0;

	if (using_model_bcube) {
		cube_t const model_bcube(get_all_models_bcube());
		model_czmin = min(czmin, model_bcube.d[2][0]);
		model_czmax = max(czmax, model_bcube.d[2][1]);
		czmin = min(czmin, model_czmin);
		czmax = max(czmax, model_czmax);
	}
	if (tree_mode & 2) {add_small_tree_coll_objs();}
	if (has_scenery2)  {add_scenery_cobjs();}
	return 1;
}


string texture_str(int tid) {

	//ostringstream oss; oss << tid; return oss.str();
	if (tid < 0) {return "none";} // or -1
	assert((unsigned)tid < textures.size());
	return textures[tid].name;
}

void light_source::write_to_cobj_file(ostream &out, bool is_diffuse) const {

	// 'L'/"light": // point/spot/line light: ambient_size diffuse_size xpos ypos zpos color [direction|pos2 [beamwidth=1.0 [inner_radius=0.0 [is_line_light=0 [use_shadow_map=0]]]]]
	out << "light " << (is_diffuse ? 0.0 : radius) << " " << (is_diffuse ? radius : 0.0) << " " << pos.raw_str() << " " << color.raw_str() << " "
		<< (is_line_light() ? pos2 : dir).raw_str() << " " << bwidth << " " << r_inner << " " << is_line_light();
}
void light_source_trig::write_to_cobj_file(ostream &out, bool is_diffuse) const {

	triggers.write_to_cobj_file(out);
	indir_dlight_group_manager.write_entry_to_cobj_file(indir_dlight_ix, out);
	light_source::write_to_cobj_file(out, is_diffuse);
	out << " " << (use_smap != 0) << endl;
	if (bound) {out << "bind_light " << bind_pos.raw_str() << endl;} // 'V'/"bind_light": // bind prev light source to cobj at location <x y z>
	triggers.write_end_triggers_cobj_file(out);
	out << endl; // separate with a blank line
}

void indir_dlight_group_manager_t::write_entry_to_cobj_file(unsigned tag_ix, ostream &out) const {
	if (tag_ix == 0) {out << "indir_dlight_group none" << endl; return;}
	assert(tag_ix < groups.size());
	out << "indir_dlight_group " << groups[tag_ix].filename << " " << groups[tag_ix].scale << endl; // 'U'/"indir_dlight_group": name [scale]
}

void coll_obj::write_to_cobj_file(ostream &out, coll_obj &prev) const {

	if (type == COLL_NULL || !fixed) return; // unused/non-fixed cobj
	if (is_near_zero_area() || min_len() < 1.0E-5) return; // near zero area (use lower tolerance due to limited file write precision)
	if (!cp.draw) return; // don't write out non-drawn collision hulls as they're generally not useful by themselves
	bool const diff2(cp.draw != prev.cp.draw || cp.refract_ix != prev.cp.refract_ix || cp.light_atten != prev.cp.light_atten || cp.is_emissive != prev.cp.is_emissive);

	if (diff2 || cp.elastic != prev.cp.elastic || cp.color != prev.cp.color || cp.tid != prev.cp.tid) { // material parameters changed
		out << "l " << cp.elastic << " " << cp.color.raw_str() << " " << texture_str(cp.tid);
		if (diff2) {out << " " << cp.draw << " " << cp.refract_ix << " " << cp.light_atten << " " << cp.is_emissive;} // uncommon optional values
		out << endl;
	}
	if (cp.shine != prev.cp.shine || cp.spec_color != prev.cp.spec_color) {out << "r 1.0 " << cp.shine << " " << cp.spec_color.raw_str() << endl;}
	if (cp.density != prev.cp.density) {out << "density " << cp.density << endl;}
	if (cp.tscale  != prev.cp.tscale ) {out << "y " << cp.tscale << endl;}
	if (cp.tdx != prev.cp.tdx || cp.tdy != prev.cp.tdy || cp.swap_txy() != prev.cp.swap_txy()) {out << "Y " << cp.tdx << " " << cp.tdy << " " << cp.swap_txy() << endl;}
	if (cp.metalness != prev.cp.metalness) {out << "metalness " << cp.metalness << endl;}

	if (cp.normal_map != prev.cp.normal_map) {
		out << "X " << texture_str(cp.normal_map);
		if (cp.normal_map >= 0) {out << " " << (textures[textures.size()].invert_y != 0) << " " << cp.negate_nm_bns();}
		out << endl;
	}
	if (cp.surfs != prev.cp.surfs) {out << "e " << (unsigned)cp.surfs << endl;}
	if (is_movable() != prev.is_movable()) {out << "movable " << is_movable() << endl;} // or 'd'
	if (destroy != prev.destroy) {out << "a " << (unsigned)destroy << endl;}
	if (is_reflective()) {out << "cube_map_ref 1" << endl;} // cube map reflections only

	switch (type) {
	case COLL_CUBE: // 'B': cube: xmin xmax ymin ymax zmin zmax [corner_radius]
		out << "B " << d[0][0] << " " << d[0][1] << " " << d[1][0] << " " << d[1][1] << " " << d[2][0] << " " << d[2][1];
		if (radius2 > 0.0) {out << " " << radius2;} // rounded corner
		out << endl;
		break;
	case COLL_SPHERE: // 'S': sphere: x y z radius
		out << "S " << points[0].raw_str() << " " << radius << endl;
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT: // 'C': cylinder: x1 y1 z1 x2 y2 z2 r1 r2
		out << "C " << points[0].raw_str() << " " << points[1].raw_str() << " " << radius << " " << radius2 << endl;
		break;
	case COLL_CAPSULE: // 'k': capsule: x1 y1 z1 x2 y2 z2 r1 r2
		out << "k " << points[0].raw_str() << " " << points[1].raw_str() << " " << radius << " " << radius2 << endl;
		break;
	case COLL_TORUS: // 'z': x y z dir_x dir_y dir_z ro ri
		out << "z " << points[0].raw_str() << " " << norm.raw_str() << " " << radius << " " << radius2 << endl;
		break;
	case COLL_POLYGON: // 'P': polygon: npts (x y z)* thickness
		out << "P " << npoints << " ";
		for (int n = 0; n < npoints; ++n) {out << points[n].raw_str() << " ";}
		out << thickness << endl;
		break;
	default: assert(0);
	}
	if (is_reflective()) {out << "cube_map_ref 0" << endl;} // disable reflections
	prev = *this; // update previous cobj
}

bool write_coll_objects_file(coll_obj_group const &cobjs, string const &fn) { // call on fixed_cobjs

	ofstream out(fn);
	if (!out.good()) {cerr << "Error opening coll object file '" << fn << "' for output" << endl; return 0;}
	
	// add normal drawn cobjs
	coll_obj prev_cobj; // starts as default cobj

	for (auto c = cobjs.begin(); c != cobjs.end(); ++c) {
		if (c->platform_id >= 0) continue; // platforms are written out below
		c->write_to_cobj_file(out, prev_cobj);
	}
	out << endl;
	
	// add platforms
	platforms.add_current_cobjs(); // to ensure cobjs are valid for each platform

	for (auto p = platforms.begin(); p != platforms.end(); ++p) {
		p->write_to_cobj_file(out);

		for (auto c = p->cobjs.begin(); c != p->cobjs.end(); ++c) {
			assert(*c < cobjs.size());
			cobjs[*c].write_to_cobj_file(out, prev_cobj);
		}
		out << endl; // separate with a blank line for readability
		for (auto l = p->lights.begin(); l != p->lights.end(); ++l) {
			assert(*l < light_sources_d.size());
			light_sources_d[*l].write_to_cobj_file(out, 1);
		}
	}
	if (!platforms.empty()) {out << "Q 0" << endl << endl;} // end platforms section
	
	// add teleporters, portals, and appearance spots
	for (auto t = teleporters.begin(); t != teleporters.end(); ++t) {
		out << "x " << t->pos.raw_str() << " " << t->dest.raw_str() << " " << t->radius << endl;
	}
	out << endl;
	for (auto p = portals.begin(); p != portals.end(); ++p) {
		out << "N ";
		for (unsigned i = 0; i < 4; ++i) {out << p->pts[i].raw_str() << " ";}
		out << p->normal.raw_str() << endl;
	}
	out << endl;
	for (auto i = app_spots.begin(); i != app_spots.end(); ++i) {out << "A " << i->raw_str() << endl;}
	out << endl;
	
	// add light sources
	for (auto l = light_sources_a.begin(); l != light_sources_a.end(); ++l) {
		l->write_to_cobj_file(out, 0); out << endl;
	}
	out << endl;
	for (auto l = light_sources_d.begin()+FLASHLIGHT_LIGHT_ID+1; l != light_sources_d.end(); ++l) {
		if (l->has_bound_platform()) continue; // already output during the platform pass
		l->write_to_cobj_file(out, 1);
	}

	// 'b': // cube volume light (for sky/global indirect): x1 x2 y1 y2 z1 z2  color.R color.G color.B  intensity num_rays ltype [disabled_edges_bits]
	for (unsigned d = 0; d < 2; ++d) {
		auto lsv(d ? global_cube_lights : sky_cube_lights);
		for (auto i = lsv.begin(); i != lsv.end(); ++i ) {
			//read_cube(fp, xf, cls.bounds);
			//fscanf(fp, "%f%f%f%f%u%i%u", &cls.color.R, &cls.color.G, &cls.color.B, &cls.intensity, &cls.num_rays, &ivals[0], &cls.disabled_edges);
		}
	}
	write_trees_to_cobj_file(out);
	out << endl;
	write_small_trees_to_cobj_file(out);
	out << endl;
	write_plants_to_cobj_file(out);
	out << endl << "end" << endl;
	return 1;
}


void init_models() {

	build_hmv_shape();
	gen_star_points();
}


void free_models() {
	delete_hmv_shape();
}


void gen_star_points() {

	for (unsigned i = 0; i < 2*N_STAR_POINTS; ++i) {
		float const angle(TWO_PI*((float)i/(float)N_STAR_POINTS));
		
		if (!(i&1)) { // outer point
			star_pts[i].x = cosf(angle);
			star_pts[i].y = sinf(angle);
		}
		else { // inner point
			star_pts[i].x = STAR_INNER_RAD*cosf(angle);
			star_pts[i].y = STAR_INNER_RAD*sinf(angle);
		}
		star_pts[i].z = 0.0;
	}
}


