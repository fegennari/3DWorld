// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/10/02

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"


float    const KILL_DEPTH          = 12.0;
float    const RECOVER_DEPTH       = 1.0;
float    const MIN_BOUNCE_VEL      = 2.0;
float    const BOUNCE_CUTOFF       = MIN_BOUNCE_VEL*MIN_BOUNCE_VEL;
float    const WATER_SURF_FRICTION = 0.95;
float    const SURF_ADV_STEP       = 2.0;
float    const MAX_TEMP            = 37.0;
float    const MIN_TEMP            = -18.0;
float    const ICE_BOUNCE_ELAS     = 0.4;
float    const ICE_ELASTICITY      = 0.95;
float    const WATER_ELASTIC       = 0.6;
float    const LAND_ELASTICITY     = 0.8;
float    const SPILL_ELASTIC       = 0.8; // also multiplied by LAND_ELASTICITY
float    const WATER_DAMPING       = 0.1;
float    const CRITICAL_ANGLE      = 0.5; // in radians, for skipping objects on water
float    const BURN_DAMAGE         = 1200.0;
unsigned const MAX_FIRE_TIME       = 10000;
bool     const ball_camera_view    = 0;
bool     const PRINT_TIME_OF_DAY   = 1;


// Global Variables
bool no_sun_lpos_update(0);
int TIMESCALE2(TIMESCALE), I_TIMESCALE2(0);
float temperature(DEF_TEMPERATURE), max_obj_radius(0.0);
float TIMESTEP(DEF_TIMESTEP), orig_timestep(DEF_TIMESTEP); // temp in degrees C
float cloud_cover(0.0);
vector3d wind(0.4, 0.2, 0.0), total_wind(0.0, 0.0, 0.0);
point flow_source(0.0, 0.0, -2.0);
obj_type object_types[NUM_TOT_OBJS];

extern int num_groups, display_mode, frame_counter, game_mode, camera_coll_id, precip_mode;
extern int s_ball_id, world_mode, has_accumulation, has_snow_accum, iticks, auto_time_adv, DISABLE_WATER, enable_fsource, animate2;
extern float max_water_height, zmin, zmax, ztop, zbottom, zmax_est, base_gravity, tstep, fticks, water_plane_z;
extern float sun_rot, moon_rot, alt_temp, light_factor, XY_SCENE_SIZE, TWO_XSS, TWO_YSS, czmax, grass_length;
extern vector3d up_norm, orig_cdir;
extern vector<valley> valleys;
extern int coll_id[];
extern dwobject def_objects[];
extern obj_group obj_groups[];
extern obj_vector_t<bubble> bubbles;
extern obj_vector_t<particle_cloud> part_clouds;
extern obj_vector_t<fire> fires;
extern obj_vector_t<decal_obj> decals;
extern water_particle_manager water_part_man;
extern physics_particle_manager explosion_part_man;


int get_obj_zval(point &pt, float &dz, float z_offset);
int snow_height(point pos);


float get_max_t(int obj_type) {return object_types[obj_type].max_t;}


void init_objects() {

	static bool inited(0);
	assert(TIMESTEP > 0.0);
	if (!inited) inited = 1; else return; // prevent multiple inits
	float rmax(0.0);
	dwobject obj;
	obj.pos         = all_zeros;
	obj.velocity    = zero_vector;
	obj.status      = 1; // starts airborne
	obj.time        = 0;
	obj.flags       = 0;
	obj.coll_id     = -1;
	obj.direction   = 0;
	obj.source      = NO_SOURCE; // undefined/natural
	obj.angle       = 0.0;
	obj.orientation = up_norm; // up (+z) needed for spin, not yet used
	obj.init_dir    = up_norm;
	obj.vdeform     = all_ones;

	for (unsigned i = 0; i < NUM_TOT_OBJS; ++i) {
		object_types[i].tid             = -1; // texture enum, not gl texture
		object_types[i].radius          = 0.001;
		object_types[i].air_factor      = 0.001;
		object_types[i].friction_factor = 0.0;
		object_types[i].damage          = 0.0;
		object_types[i].min_t           = -100;
		object_types[i].max_t           = 200;
		object_types[i].density         = 1.0;
		object_types[i].gravity         = 1.0;
		object_types[i].elasticity      = 0.0;
		object_types[i].damage          = 0.0;
		object_types[i].deform          = 0.0;
		object_types[i].def_recover     = 0.0;
	}
	object_types[RAIN].air_factor          = 0.2;
	object_types[RAIN].friction_factor     = 0.005;
	object_types[RAIN].radius              = 0.0035;
	object_types[RAIN].lifetime            = int(50 + 40*Z_SCENE_SIZE);
	object_types[RAIN].min_t               = RAIN_MIN_TEMP;
	object_types[RAIN].max_t               = WATER_MAX_TEMP;
	object_types[RAIN].density             = 1.0;
	object_types[RAIN].health              = 5.0; // survivability from surface collision
	object_types[RAIN].color               = WATER;
	object_types[RAIN].color.alpha         = 0.8*WATER_ALPHA;
	object_types[RAIN].flags               = SEMI_TRANSPARENT | BLEND | FALL_EVERYWHERE | LOW_SPECULAR | TAIL_WHEN_FALL | IS_PRECIP | OBJ_IS_DROP;

	object_types[SNOW].air_factor          = 0.4;
	object_types[SNOW].friction_factor     = 3.0; //2.49; //2.2
	object_types[SNOW].radius              = 0.002;
	object_types[SNOW].lifetime            = int(150 + 70*Z_SCENE_SIZE);
	object_types[SNOW].min_t               = -1000;
	object_types[SNOW].max_t               = SNOW_MAX_TEMP;
	object_types[SNOW].density             = 0.2;
	object_types[SNOW].health              = 10.0;
	object_types[SNOW].color               = SNOW_COLOR;
	object_types[SNOW].color.alpha         = SNOW_ALPHA;
	object_types[SNOW].flags               = BLEND | FALL_EVERYWHERE | IS_PRECIP | OBJ_IS_DROP;
	object_types[SNOW].tid                 = SNOWFLAKE_TEX;

	object_types[HAIL].air_factor          = 0.175;
	object_types[HAIL].friction_factor     = 0.1;
	object_types[HAIL].radius              = 0.0045;
	object_types[HAIL].lifetime            = int(100 + 50*Z_SCENE_SIZE);
	object_types[HAIL].min_t               = SNOW_MAX_TEMP;
	object_types[HAIL].max_t               = RAIN_MIN_TEMP;
	object_types[HAIL].density             = 0.8;
	object_types[HAIL].elasticity          = 0.8;
	object_types[HAIL].health              = 30.0;
	object_types[HAIL].color               = TRANS;
	object_types[HAIL].color.alpha         = 0.9*ICE_ALPHA;
	object_types[HAIL].flags               = FALL_EVERYWHERE | IS_PRECIP;

	object_types[LEAF].air_factor          = 0.5;
	object_types[LEAF].friction_factor     = 0.9;
	object_types[LEAF].gravity             = 0.01;
	object_types[LEAF].radius              = 0.005;
	object_types[LEAF].lifetime            = 1200;
	object_types[LEAF].density             = 0.3;
	object_types[LEAF].health              = 4000.0;
	object_types[LEAF].color               = LEAF_C;
	object_types[LEAF].flags               = NO_WATER_DAMAGE | OBJ_IS_FLAT;
	object_types[LEAF].tid                 = LEAF_TEX;

	object_types[BALL].air_factor          = 0.05;
	object_types[BALL].friction_factor     = 0.02;
	object_types[BALL].gravity             = 0.3; // < half gravity for dodgeball
	object_types[BALL].radius              = 0.042;
	object_types[BALL].damage              = 2200.0;
	object_types[BALL].lifetime            = 1600;
	object_types[BALL].density             = 0.7;
	object_types[BALL].elasticity          = 0.92;
	object_types[BALL].health              = 20000.0;
	object_types[BALL].color               = WHITE;
	object_types[BALL].flags               = SPECULAR | SELECTABLE | OBJ_ROLLS;
	object_types[BALL].tid                 = SKULL_TEX;

	object_types[S_BALL].air_factor        = 0.035;
	object_types[S_BALL].friction_factor   = 0.02;
	object_types[S_BALL].radius            = 0.015;
	object_types[S_BALL].damage            = 8.0;
	object_types[S_BALL].lifetime          = 210;
	object_types[S_BALL].density           = 0.4;
	object_types[S_BALL].elasticity        = 0.9;
	object_types[S_BALL].deform            = 0.5;
	object_types[S_BALL].def_recover       = 0.5;
	object_types[S_BALL].health            = 10000.0;
	object_types[S_BALL].color             = ORANGE;
	object_types[S_BALL].flags             = SPECULAR | SELECTABLE | DEFORMABLE;

	object_types[SMILEY].air_factor        = 0.005;
	object_types[SMILEY].friction_factor   = 0.001;
	object_types[SMILEY].gravity           = 0.4;
	object_types[SMILEY].radius            = CAMERA_RADIUS;
	object_types[SMILEY].damage            = 50.0;
	object_types[SMILEY].lifetime          = 1000000;
	object_types[SMILEY].density           = 1.2;
	object_types[SMILEY].elasticity        = 0.4;
	object_types[SMILEY].health            = 100.0;
	object_types[SMILEY].color             = YELLOW;
	object_types[SMILEY].flags             = SELECTABLE | NO_FALL | VERTEX_DEFORM | NO_WATER_DAMAGE | SPECULAR;
	object_types[SMILEY].min_t             = -50.0;
	object_types[SMILEY].max_t             = 100.0;

	object_types[BLOOD].air_factor         = 0.15;
	object_types[BLOOD].friction_factor    = 0.01;
	object_types[BLOOD].radius             = 0.004;
	object_types[BLOOD].lifetime           = int(90 + 10*Z_SCENE_SIZE);
	object_types[BLOOD].max_t              = WATER_MAX_TEMP;
	object_types[BLOOD].density            = 1.0;
	object_types[BLOOD].health             = 15.0; // survivability from surface collision
	object_types[BLOOD].color              = BLOOD_C;
	object_types[BLOOD].flags              = LOW_SPECULAR | BLEND | OBJ_IS_DROP;

	object_types[CHARRED].air_factor       = 0.2;
	object_types[CHARRED].friction_factor  = 0.1;
	object_types[CHARRED].radius           = 0.003;
	object_types[CHARRED].lifetime         = 240;
	object_types[CHARRED].max_t            = 1000;
	object_types[CHARRED].density          = 0.9;
	object_types[CHARRED].elasticity       = 0.5;
	object_types[CHARRED].health           = 20.0;
	object_types[CHARRED].color            = colorRGBA(0.01, 0.01, 0.01, 1.0); // nearly BLACK, but that seems to hit a driver bug
	object_types[CHARRED].flags            = BLEND;

	object_types[CHUNK].air_factor         = 0.08;
	object_types[CHUNK].friction_factor    = 1.2;
	object_types[CHUNK].gravity            = 0.6;
	object_types[CHUNK].radius             = 0.011;
	object_types[CHUNK].lifetime           = 500;
	object_types[CHUNK].density            = 1.2;
	object_types[CHUNK].elasticity         = 0.9;
	object_types[CHUNK].deform             = 0.5;
	object_types[CHUNK].def_recover        = 0.0;
	object_types[CHUNK].health             = 250.0;
	object_types[CHUNK].color              = YELLOW;
	object_types[CHUNK].flags              = DEFORMABLE;

	object_types[SFPART].air_factor        = 0.1; // eye, nose, tongue
	object_types[SFPART].friction_factor   = 0.06;
	object_types[SFPART].radius            = CAMERA_RADIUS/6.0;
	object_types[SFPART].lifetime          = 1000;
	object_types[SFPART].density           = 0.5;
	object_types[SFPART].elasticity        = 0.9;
	object_types[SFPART].health            = 200.0;
	object_types[SFPART].color             = BLACK;
	object_types[SFPART].flags             = SPECULAR | NO_WATER_DAMAGE;

	object_types[DROPLET].air_factor       = 0.17; // blood/water
	object_types[DROPLET].friction_factor  = 0.005;
	object_types[DROPLET].radius           = 0.0035;
	object_types[DROPLET].lifetime         = int(40 + 10*Z_SCENE_SIZE);
	object_types[DROPLET].min_t            = W_FREEZE_POINT;
	object_types[DROPLET].max_t            = WATER_MAX_TEMP;
	object_types[DROPLET].density          = 1.0;
	object_types[DROPLET].health           = 5.0; // survivability from surface collision
	object_types[DROPLET].color            = WATER;
	object_types[DROPLET].color.alpha      = 1.0*WATER_ALPHA;
	object_types[DROPLET].flags            = SEMI_TRANSPARENT | BLEND | FALL_EVERYWHERE | LOW_SPECULAR | OBJ_IS_DROP;

	object_types[WDROPLET].air_factor      = 0.17; // water from springs
	object_types[WDROPLET].friction_factor = 0.005;
	object_types[WDROPLET].radius          = 0.0035;
	object_types[WDROPLET].lifetime        = int(50 + 20*Z_SCENE_SIZE);
	object_types[WDROPLET].min_t           = W_FREEZE_POINT;
	object_types[WDROPLET].max_t           = WATER_MAX_TEMP;
	object_types[WDROPLET].density         = 1.0;
	object_types[WDROPLET].health          = 10.0; //0.0; // survivability from surface collision (or not)
	object_types[WDROPLET].color           = WATER;
	object_types[WDROPLET].color.alpha     = 1.0*WATER_ALPHA;
	object_types[WDROPLET].flags           = SEMI_TRANSPARENT | BLEND | LOW_SPECULAR | OBJ_IS_DROP;

	object_types[SAND].air_factor          = 0.5;
	object_types[SAND].friction_factor     = 0.3;
	object_types[SAND].radius              = 0.005;
	object_types[SAND].lifetime            = 80;
	object_types[SAND].density             = 1.5;
	object_types[SAND].elasticity          = 0.3;
	object_types[SAND].health              = 300.0;
	object_types[SAND].color               = BROWN;
	object_types[SAND].flags               = BLEND;
	object_types[SAND].tid                 = SAND_TEX;

	object_types[DIRT].air_factor          = 0.25;
	object_types[DIRT].friction_factor     = 0.8;
	object_types[DIRT].radius              = 0.008;
	object_types[DIRT].lifetime            = 700;
	object_types[DIRT].density             = 1.1;
	object_types[DIRT].elasticity          = 0.6;
	object_types[DIRT].health              = 400.0;
	object_types[DIRT].color               = DK_BROWN;
	object_types[DIRT].flags               = BLEND;
	object_types[DIRT].tid                 = DIRT_TEX;

	object_types[ROCK].air_factor          = 0.1;
	object_types[ROCK].friction_factor     = 0.9;
	object_types[ROCK].radius              = 0.0099;
	object_types[ROCK].damage              = 1.0;
	object_types[ROCK].lifetime            = 800;
	object_types[ROCK].density             = 1.8;
	object_types[ROCK].elasticity          = 0.5;
	object_types[ROCK].health              = 1000.0;
	object_types[ROCK].color               = DK_GRAY;
	object_types[ROCK].flags               = BLEND;
	object_types[ROCK].tid                 = ROCK_TEX;

	object_types[FRAGMENT].air_factor      = 0.06;
	object_types[FRAGMENT].friction_factor = 0.6;
	object_types[FRAGMENT].radius          = 0.0098;
	object_types[FRAGMENT].damage          = 1.0;
	object_types[FRAGMENT].lifetime        = 400;
	object_types[FRAGMENT].density         = 2.2;
	object_types[FRAGMENT].elasticity      = 0.4;
	object_types[FRAGMENT].health          = 800.0;
	object_types[FRAGMENT].color           = WHITE;
	object_types[FRAGMENT].flags           = BLEND; //| OBJ_IS_FLAT;

	object_types[PARTICLE].air_factor      = 0.08;
	object_types[PARTICLE].friction_factor = 0.2;
	object_types[PARTICLE].radius          = 0.005;
	object_types[PARTICLE].damage          = 0.0;
	object_types[PARTICLE].lifetime        = 100;
	object_types[PARTICLE].max_t           = 1000.0;
	object_types[PARTICLE].density         = 1.6;
	object_types[PARTICLE].elasticity      = 0.6;
	object_types[PARTICLE].gravity         = 0.4;
	object_types[PARTICLE].health          = 1000.0;
	object_types[PARTICLE].color           = BLACK;
	object_types[PARTICLE].tid             = BLUR_TEX;
	object_types[PARTICLE].flags           = BLEND | SEMI_TRANSPARENT;

	object_types[CAMERA].air_factor        = 0.005;
	object_types[CAMERA].friction_factor   = 0.001;
	object_types[CAMERA].gravity           = 0.4;
	object_types[CAMERA].radius            = CAMERA_RADIUS;
	object_types[CAMERA].lifetime          = 0;
	object_types[CAMERA].density           = 1.2;
	object_types[CAMERA].elasticity        = 0.4;
	object_types[CAMERA].health            = 100.0;
	object_types[CAMERA].color             = WHITE;
	object_types[CAMERA].flags             = NO_FALL | NO_WATER_DAMAGE;
	object_types[CAMERA].min_t             = -50.0;
	object_types[CAMERA].max_t             = 100.0;

	object_types[ROCKET].air_factor        = 0.02;
	object_types[ROCKET].friction_factor   = 0.5;
	object_types[ROCKET].gravity           = 0.0;
	object_types[ROCKET].radius            = 0.022;
	object_types[ROCKET].damage            = 400.0;
	object_types[ROCKET].lifetime          = 300;
	object_types[ROCKET].density           = 1.3;
	object_types[ROCKET].health            = 10.0;
	object_types[ROCKET].color             = GRAY;
	object_types[ROCKET].flags             = SPECULAR | SELECTABLE | OBJ_EXPLODES | EXPL_ON_COLL | COLL_DESTROYS;

	object_types[LANDMINE].air_factor      = 0.01;
	object_types[LANDMINE].friction_factor = 3.5;
	object_types[LANDMINE].radius          = 0.03;
	object_types[LANDMINE].damage          = 1200.0;
	object_types[LANDMINE].lifetime        = 1400;
	object_types[LANDMINE].density         = 2.5;
	object_types[LANDMINE].health          = 80.0;
	object_types[LANDMINE].color           = WHITE;
	object_types[LANDMINE].flags           = SPECULAR | SELECTABLE | OBJ_EXPLODES;
	object_types[LANDMINE].tid             = CAMOFLAGE_TEX;

	object_types[SEEK_D].air_factor        = 0.015;
	object_types[SEEK_D].friction_factor   = 0.5;
	object_types[SEEK_D].gravity           = 0.0;
	object_types[SEEK_D].radius            = 0.026;
	object_types[SEEK_D].damage            = 1200.0;
	object_types[SEEK_D].lifetime          = 400;
	object_types[SEEK_D].density           = 1.4;
	object_types[SEEK_D].health            = 12.0;
	object_types[SEEK_D].color             = BLACK;
	object_types[SEEK_D].flags             = SPECULAR | SELECTABLE | OBJ_EXPLODES | EXPL_ON_COLL | COLL_DESTROYS;

	object_types[STAR5].air_factor         = 0.1;
	object_types[STAR5].friction_factor    = 2.49;
	object_types[STAR5].gravity            = 0.12;
	object_types[STAR5].radius             = 0.005;
	object_types[STAR5].damage             = 40.0;
	object_types[STAR5].lifetime           = 500;
	object_types[STAR5].density            = 1.8;
	object_types[STAR5].elasticity         = 0.1;
	object_types[STAR5].health             = 600.0;
	object_types[STAR5].color              = GRAY; // silver?
	object_types[STAR5].flags              = SPECULAR | SELECTABLE; // OBJ_IS_FLAT?

	object_types[PLASMA].air_factor        = 0.012;
	object_types[PLASMA].friction_factor   = 1.5;
	object_types[PLASMA].gravity           = 0.04;
	object_types[PLASMA].radius            = 0.02;
	object_types[PLASMA].damage            = 250.0;
	object_types[PLASMA].lifetime          = 700;
	object_types[PLASMA].min_t             = -1000;
	object_types[PLASMA].max_t             = 10000;
	object_types[PLASMA].density           = 0.12;
	object_types[PLASMA].elasticity        = 0.2;
	object_types[PLASMA].health            = 0.1;
	object_types[PLASMA].color             = YELLOW;
	object_types[PLASMA].color.alpha       = 0.6; // average color
	object_types[PLASMA].flags             = SELECTABLE | BLEND | SEMI_TRANSPARENT | OBJ_EXPLODES | EXPL_ON_COLL;
	object_types[PLASMA].tid               = PLASMA_TEX;

	object_types[GRENADE].air_factor       = 0.07;
	object_types[GRENADE].friction_factor  = 0.35;
	object_types[GRENADE].gravity          = 0.6;
	object_types[GRENADE].radius           = 0.02;
	object_types[GRENADE].damage           = 15.0;
	object_types[GRENADE].lifetime         = 35;
	object_types[GRENADE].density          = 1.2;
	object_types[GRENADE].elasticity       = 0.6;
	object_types[GRENADE].health           = 70.0;
	object_types[GRENADE].color            = COPPER_C;
	object_types[GRENADE].flags            = SPECULAR | SELECTABLE | OBJ_EXPLODES;

	object_types[CGRENADE].air_factor      = 0.06; // cluster of 10 grenades
	object_types[CGRENADE].friction_factor = 0.45;
	object_types[CGRENADE].gravity         = 0.5;
	object_types[CGRENADE].radius          = 0.04;
	object_types[CGRENADE].damage          = 40.0;
	object_types[CGRENADE].lifetime        = 25;
	object_types[CGRENADE].density         = 1.2;
	object_types[CGRENADE].elasticity      = 0.6;
	object_types[CGRENADE].health          = 50.0;
	object_types[CGRENADE].color           = GOLD;
	object_types[CGRENADE].flags           = SPECULAR | SELECTABLE | OBJ_EXPLODES;

	object_types[SHRAPNEL].air_factor      = 0.05;
	object_types[SHRAPNEL].friction_factor = 2.49;
	object_types[SHRAPNEL].radius          = 0.003;
	object_types[SHRAPNEL].damage          = 4.0;
	object_types[SHRAPNEL].lifetime        = 300;
	object_types[SHRAPNEL].density         = 4.0;
	object_types[SHRAPNEL].elasticity      = 0.5;
	object_types[SHRAPNEL].health          = 1000.0;
	object_types[SHRAPNEL].color           = BLACK;
	object_types[SHRAPNEL].flags           = 0; // OBJ_IS_FLAT?

	object_types[SHELLC].air_factor        = 0.12;
	object_types[SHELLC].friction_factor   = 0.9;
	object_types[SHELLC].gravity           = 0.5;
	object_types[SHELLC].radius            = 0.0015;
	object_types[SHELLC].lifetime          = 400;
	object_types[SHELLC].density           = 1.8;
	object_types[SHELLC].elasticity        = 0.8;
	object_types[SHELLC].health            = 1000.0;
	object_types[SHELLC].color             = BRASS_C;
	object_types[SHELLC].flags             = SPECULAR | OBJ_IS_CYLIN;

	object_types[PROJC].air_factor         = 0.0;
	object_types[PROJC].friction_factor    = 0.0;
	object_types[PROJC].gravity            = 0.0;
	object_types[PROJC].radius             = 0.005;
	object_types[PROJC].lifetime           = 1000;
	object_types[PROJC].density            = 1.0;
	object_types[PROJC].health             = 1.0;
	object_types[PROJC].color              = BLACK; // never actually drawn
	object_types[PROJC].flags              = 0;

	object_types[SKULL].friction_factor    = 0.2;
	object_types[SKULL].radius             = 0.7*CAMERA_RADIUS;
	object_types[SKULL].lifetime           = 600;
	object_types[SKULL].density            = 1.1;
	object_types[SKULL].elasticity         = 0.75;
	object_types[SKULL].health             = 25.0;
	object_types[SKULL].color              = LT_GRAY;
	object_types[SKULL].flags              = SPECULAR | SELECTABLE | NO_WATER_DAMAGE;
	object_types[SKULL].tid                = SMILEY_SKULL_TEX;

	object_types[HEALTH].color             = MAGENTA;
	object_types[SHIELD].color             = YELLOW;
	object_types[POWERUP].color            = WHITE; // should be variable
	object_types[WEAPON].color             = CYAN;
	object_types[AMMO].color               = BLUE;
	object_types[WA_PACK].color            = LT_BROWN;

	object_types[GASSED].gravity           = 0.0;
	object_types[GASSED].color             = OLIVE;
	object_types[GASSED].tid               = YUCK_TEX;
	object_types[GASSED].radius            = 0.035; // 0.5*weapons[W_GASSER].blast_radius
	object_types[GASSED].air_factor        = 0.5;
	object_types[GASSED].min_t             = -1000;
	object_types[GASSED].max_t             = 1000;
	object_types[GASSED].density           = 0.01;

	object_types[WAYPOINT].radius          = CAMERA_RADIUS;

	object_types[SAWBLADE].air_factor      = 0.01;
	object_types[SAWBLADE].friction_factor = 0.01;
	object_types[SAWBLADE].gravity         = 0.01;
	object_types[SAWBLADE].radius          = 0.04;
	object_types[SAWBLADE].damage          = 100.0;
	object_types[SAWBLADE].lifetime        = 400;
	object_types[SAWBLADE].density         = 0.4; // to reduce damage
	object_types[SAWBLADE].elasticity      = 2.0; // 100% elastic collisions with any surface that has elasticity >= 0.5
	object_types[SAWBLADE].health          = 500.0;
	object_types[SAWBLADE].color           = WHITE;
	object_types[SAWBLADE].flags           = SPECULAR | SELECTABLE | BLEND; // OBJ_IS_FLAT?

	for (unsigned i = HEALTH; i <= WA_PACK; ++i) { // all other physics are the same
		object_types[i].air_factor      = 0.05;
		object_types[i].friction_factor = 0.9;
		object_types[i].radius          = 0.025;
		object_types[i].lifetime        = 2400;
		object_types[i].density         = 0.05;
		object_types[i].elasticity      = 0.4;
		object_types[i].health          = 1600.0;
		object_types[i].color.alpha     = (i == WEAPON || i == AMMO || i == WA_PACK) ? 0.3 : 0.6;
		object_types[i].flags           = SEMI_TRANSPARENT | BLEND | SPECULAR | SELECTABLE;
	}
	object_types[HEALTH].damage   = -500.0;
	object_types[WA_PACK].density = 0.5; // so a blast doesn't send it flying
	object_types[PRECIP].flags   |= IS_PRECIP;

	for (unsigned i = 0; i < NUM_TOT_OBJS; ++i) { // i <= CAMERA?
		def_objects[i]        = obj;
		def_objects[i].type   = i;
		def_objects[i].health = object_types[i].health;
		float const radius(object_types[i].radius);
		rmax = max(rmax, radius);

		// objects are assumed to be spherical
		object_types[i].surface_area = 4.0*PI*radius*radius;
		object_types[i].volume       = (4.0/3.0)*PI*radius*radius*radius;
		object_types[i].terminal_vel = 1.0/max(1.0E-6f, object_types[i].air_factor);
		object_types[i].mass         = 150000.0*object_types[i].density*object_types[i].volume;
		object_types[i].lifetime     = int((0.01/TIMESTEP)*object_types[i].lifetime);
	}
	def_objects[SMILEY].orientation.assign(1.0, 0.0, 0.0); // start facing in +x
	object_types[BALL].mass             = 1.0; // fudge the mass so that the splashes look better
	object_types[ROCKET].terminal_vel   = 0.5;
	object_types[SEEK_D].terminal_vel   = 0.1;
	object_types[STAR5].terminal_vel    = 2.5;
	object_types[FIRE].friction_factor  = 2.0; // not a real object
	object_types[FIRE].terminal_vel     = 1.5;
	object_types[FIRE].gravity          = 0.2;
	object_types[DYNAM_PART].elasticity = 1.0; // not a real object

	for (unsigned i = 0; i < NUM_TOT_OBJS; ++i) {
		assert(object_types[i].radius  > TOLERANCE);
		assert(object_types[i].density > TOLERANCE);
		assert(object_types[i].mass    > TOLERANCE);
	}
	if (ball_camera_view) {
		def_objects[BALL].flags |= CAMERA_VIEW;
		orig_cdir = cview_dir;
	}
	set_coll_rmax(rmax);
}


void set_coll_rmax(float rmax) {

	max_obj_radius = rmax; // only used to cache the init value of rmax for use in later calls
	//cout << "rmax = " << max_obj_radius << ", DXY = " << max(DX_VAL, DY_VAL) << ", ratio: " << max_obj_radius/max(DX_VAL, DY_VAL) << endl;
}


void change_timestep(float mult_factor) {

	assert(mult_factor > 0.0);
	TIMESTEP *= mult_factor;
	float const its((TIMESTEP/DEF_TIMESTEP)*((int)I_TIMESCALE));

	if (its < 0.5) {
		I_TIMESCALE2 = 1;
		TIMESCALE2   = min(MAX_I_TIMESCALE, int(1.0/its));
	}
	else {
		TIMESCALE2   = 1;
		I_TIMESCALE2 = min(MAX_I_TIMESCALE, int(its + 0.5));
	}
	for (int i = 0; i < NUM_TOT_OBJS; ++i) {
		if (object_types[i].lifetime/mult_factor > 10) object_types[i].lifetime = int(object_types[i].lifetime/mult_factor);
	}
	orig_timestep = TIMESTEP;
}


vector3d get_flow_velocity(point pos) { // has same effect as wind

	int const do_tornado(0), do_swirl(0);
	if (!enable_fsource) return all_zeros;
	vector3d v(flow_source, pos);
	float const dist(v.mag()), dxy(do_swirl ? sqrt(v.x*v.x + v.y*v.y) : 0.0);
	float const vmag(10.0*(0.25 + 1.0/(dist + 0.7))); // direct attraction
	float vmagxy(vmag);
	v *= (vmag/dist);

	if (do_tornado) {
		if (v.z < 0)  v.z    = -0.8*v.z/(dist + 1.0);
		if (do_swirl) vmagxy = 10.0*(0.4 + 1.0/(dxy + 0.5));
	}
	if (do_swirl) {
		double swirl(min(1.0f, 0.1f/(dxy + 0.1f))); // v X up = (v.y, -v.x, 0.0)
		v.x = (1.0 - swirl)*v.x + swirl*v.y*(vmagxy/vmag);
		v.y = (1.0 - swirl)*v.y - swirl*v.x*(vmagxy/vmag);
	}
	return v;
}


vector3d get_local_wind(int xpos, int ypos, float zval, bool no_use_mesh) {

	float pressure(1.0), hval(0.5);
	vector3d local_wind(wind);

	if (!no_use_mesh && world_mode == WMODE_GROUND) {
		if (point_outside_mesh(xpos, ypos)) return wind;
		// calculate direction of wind based on mesh orientation
		float const mh(mesh_height[ypos][xpos]);
		if (zval < mh)    return all_zeros; // under the mesh - no wind
		float const szmax(max(ztop, czmax)); // scene zmax
		if (zval > szmax) return wind; // above the top of the mesh
		float const rel_height((zval - mh)/(szmax - mh)); // 0 at mesh level, 1 at scene_ztop
		vector3d v_ortho;
		orthogonalize_dir(wind, vertex_normals[ypos][xpos], v_ortho, 0);
		v_ortho.z *= 0.1; // z component of velocity is much smaller
		pressure   = min(2.0, 0.5*(zmax - zbottom)/(zmax - mh)); // pressure is higher at the top of hills
		hval       = (1.0 - rel_height)*(1.0 - rel_height); // at surface: 1.0, middle: 0.25, top: 0.0
		local_wind = v_ortho*hval + wind*(1.0 - hval); // wind follows the surface contour when close to the mesh
	}
	// calculate wind intensity
	float const tx((xpos + xoff2 - total_wind.x/TWO_XSS)/MESH_X_SIZE);
	float const ty((ypos + yoff2 - total_wind.y/TWO_YSS)/MESH_Y_SIZE);
	float const wind_intensity(CLIP_TO_01(1.0f - 2.0f*get_texture_component(WIND_TEX, tx, ty, 0))); // roughly (0.0, 0.5)
	return local_wind*(pressure*(hval*wind_intensity + (1.0 - hval)));
}

vector3d get_local_wind(point const &pt, bool no_use_mesh) {
	return get_local_wind(get_xpos(pt.x), get_ypos(pt.y), pt.z, no_use_mesh);
}


void dwobject::do_coll_damage() {

	if (type == LANDMINE) return;
	assert(type != SMILEY);
	
	if (health <= COLL_DAMAGE) {
		disable();
	}
	else {
		health -= COLL_DAMAGE;
	}
}


float dwobject::get_true_radius() const {

	float const radius(object_types[type].radius);
	
	switch (type) {
	case FRAGMENT: return radius*vdeform.x;
	case SAND:     return radius*orientation.x;
	case DIRT:     return radius*orientation.x;
	case ROCK:     return radius*orientation.x;
	case PLASMA:   return radius*init_dir.x;
	}
	return radius;
}


// 0 = out of range/expired, 1 = airborne, 2 = collision, 3 = moving on ground, 4 = motionless
void dwobject::advance_object(bool disable_motionless_objects, int iter, int obj_index) { // returns collision status

	assert(!disabled());
	if (temperature <= ABSOLUTE_ZERO) return;
	bool const coll_last_frame((flags & OBJ_COLLIDED) != 0);
	flags &= ~OBJ_COLLIDED;
	verify_data();
	obj_type const &otype(object_types[type]);

	if (status == 0 || pos.z < zmin || time > otype.lifetime || (type == PARTICLE && is_underwater(pos))) {
		assert(type != SMILEY);
		status = 0;
		return;
	}
	if (iter == 0) time += iticks;
	bool const frozen(temperature <= W_FREEZE_POINT);
	if (frozen && type == SHRAPNEL) flags &= ~IN_WATER;

	if (!frozen && (flags & IS_ON_ICE)) { // ice melted while object was on it
		flags  |= IN_WATER;
		flags  &= ~(XYZ_STOPPED | IS_ON_ICE);
		status  = 1;
	}
	if (disable_motionless_objects && status == 4) {
		if ((flags & IS_ON_ICE) || (!(flags & (FLOATING | STATIC_COBJ_COLL)) && object_still_stopped(obj_index))) {
			point const old_pos(pos);
			check_vert_collision(obj_index, 1, iter); // needed for gameplay (already tested in object_still_stopped()?)
			pos = old_pos;
			if (disabled() || check_water_collision(velocity.z)) return;
			if (pos.z < zmin || !is_over_mesh(pos)) status = 0;
			flags &= ~Z_STOPPED;
			return;
		}
		flags &= ~XY_STOPPED;
	}
	float const radius(get_true_radius()), friction(otype.friction_factor);

	if (status == 1 || type == LANDMINE) { // airborne
		if (type == ROCKET && direction == 1) { // rapid fire rocket
			rotate_vector3d(signed_rand_vector(), 0.02*fticks*signed_rand_float(), velocity);
		}
		float air_factor(0.0);

		if (!(flags & UNDERWATER)) {
			if (flags & FLOATING) {
				if (otype.flags & OBJ_IS_FLAT) {
					//init_dir.z = 0.0;
					int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
					vector3d const wnorm(has_water(xpos, ypos) ? wat_vert_normals[ypos][xpos] : plus_z);
					set_orient_for_coll(&wnorm);
				}
				if (WATER_SURF_FRICTION < 1.0) {air_factor = (1.0 - WATER_SURF_FRICTION)*otype.air_factor;}
			}
			else {
				air_factor = otype.air_factor;
			}
		}
		if (flags & Z_STOPPED) {
			int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

			if (!point_outside_mesh(xpos, ypos) && (pos.z - radius) > water_matrix[ypos][xpos] && ((friction < 2.0*STICK_THRESHOLD) || (friction < rand_uniform(2.0, 2.5)*STICK_THRESHOLD))) {
				flags &= ~Z_STOPPED;
			}
			else {
				velocity.z = 0.0;
			}
		}
		point const old_pos(pos);
		bool const collided(coll_last_frame || fabs(velocity.z) < 1.0E-6);
		vector3d v_flow(enable_fsource ? get_flow_velocity(pos) : velocity), vtot(v_flow);
		float const vz_old(velocity.z);
		vector3d const local_wind(get_local_wind(pos));
		
		if (iter == 0) {
			if (collided) vtot.z += local_wind.z; else vtot += local_wind;
		}
		if (!(flags & Z_STOPPED)) {
			double gscale((type == PLASMA && init_dir.x != 0.0) ? 1.0/sqrt(init_dir.x) : 1.0);
			if ((flags & IN_WATER) && otype.density > WATER_DENSITY) gscale *= (otype.density - WATER_DENSITY)/otype.density;

			if (enable_fsource) {
				double const grav_well(min(1.0f, 0.1f*v_flow.mag()));

				if (-velocity.z < otype.terminal_vel) {
					velocity.z -= (1.0 - grav_well)*base_gravity*gscale*GRAVITY*tstep*otype.gravity;
					velocity.z  = grav_well*velocity.z - (1.0 - grav_well)*min(-velocity.z, otype.terminal_vel);
				}
				if (fabs(air_factor*vtot.z) > fabs(velocity.z) || ((vtot.z < 0) != (velocity.z < 0))) {
					velocity.z = (1.0 - grav_well*air_factor)*velocity.z + air_factor*vtot.z; // wind?
				}
			}
			else {
				if (-velocity.z < otype.terminal_vel) {
					velocity.z -= base_gravity*gscale*GRAVITY*tstep*otype.gravity;
					velocity.z  = -min(-velocity.z, otype.terminal_vel);
				}
				if (fabs(air_factor*local_wind.z) > fabs(velocity.z) || ((local_wind.z < 0) != (velocity.z < 0))) {
					velocity.z += air_factor*local_wind.z;
				}
			}
		}
		if (!(flags & XY_STOPPED)) {
			for (unsigned d = 0; d < 2; ++d) {
				if (fabs(air_factor*vtot[d]) > fabs(velocity[d]) || ((vtot[d] < 0) != (velocity[d] < 0))) {
					velocity[d] = (1.0 - air_factor)*velocity[d] + air_factor*vtot[d];
				}
				if (collided && iter == 0 && !(flags | IN_WATER)) { // apply static friction
					bool const stopped(friction >= 2.0*STICK_THRESHOLD || fabs(velocity[d]) <= friction);
					velocity[d] = (stopped ? 0.0 : max(0.0f, (velocity[d] + ((velocity[d] > 0.0) ? -friction : friction))));
				}
				pos[d] += tstep*velocity[d]; // move object
			}
			if (flags & FLOATING) float_downstream(pos, radius);
		}
		assert(!is_nan(tstep));
		pos.z += tstep*velocity.z;
		verify_data();

		// check collisions
		float const r2((otype.flags & COLL_DESTROYS) ? 0.25*radius : radius);
		float dz;
		int val(get_obj_zval(pos, dz, r2)); // 0 = out of simulation region, 1 = airborne, 2 = on ground

		if (val == 2 && dz > radius && !is_over_mesh(old_pos) && old_pos.z < pos.z) { // hit side of simulation region
			status = 0;
			return;
		}
		if (val == 0) {
			if (pos.z < zmin || (flags & Z_STOPPED)) {status = 0;} // out of simulation region and underwater
			return;
		}
		int const wcoll(check_water_collision(vz_old));
		vector3d cnorm;
		bool const last_stat_coll((flags & STATIC_COBJ_COLL) != 0);
		int const coll(check_vert_collision(obj_index, 1, iter, &cnorm));
		if (disabled()) return;
		if (!coll) flags &= ~Z_STOPPED; // fix for landmine no longer stuck to cobj
		
		if (wcoll) {
			if (!frozen) status = 1;
			flags &= (frozen ? ~STATIC_COBJ_COLL : ~ALL_COLL_STOPPED);
			return;
		}
		if (val == 2 && !coll) { // collision with mesh surface but not vertical surface
			if (iter == 0) {surf_collide_obj();} // only supports blood and chunks for now
			
			if (object_bounce(0, cnorm, 0.0, radius)) {
				if (radius >= LARGE_OBJ_RAD) {modify_grass_at(pos, 2.0*radius, 1);} // crush grass a lot
				status = 1;
				return; // objects bounce on mesh but not on collision objects
			}
			do_coll_damage();
			if (status == 0) return;
			bool const stopped(otype.friction_factor >= STICK_THRESHOLD || (flags & XY_STOPPED) || velocity.mag_sq() < BOUNCE_CUTOFF);
			velocity *= (stopped ? 0.0 : 0.95); // apply some damping
		}
		if (coll) {
			if (type == SAWBLADE) {destroy_coll_objs(pos, 500.0, source, IMPACT);} // shatterable but not destroyable
			bool const stat_coll((flags & STATIC_COBJ_COLL) != 0);

			if (!stat_coll || !last_stat_coll) {
				do_coll_damage();
				if (status == 0) return;
			}
			if (stat_coll && (friction >= STICK_THRESHOLD || velocity.mag_sq() < BOUNCE_CUTOFF)) {
				velocity = zero_vector;
				val = 4;
			}
			else if (status == 4) { // was set to stopped status in check_vert_collision (possibly to to stick friction)
				val = 4;
			}
		}
		status = val;
	} // end in the air
	else { // on the ground
		if (!is_over_mesh(pos)) { // rolled off the mesh - destroy it
			status = 0;
			return;
		}
		if (otype.flags & COLL_DESTROYS) {assert(type != SMILEY); status = 0; return;}
		if (flags & STATIC_COBJ_COLL) return; // stuck on vertical collision surface
		if (check_water_collision(velocity.z) && (frozen || otype.density < WATER_DENSITY)) return;
		if (otype.flags & (OBJ_IS_FLAT | OBJ_IS_CYLIN)) set_orient_for_coll(NULL);
		point const old_pos(pos);
		int const val(surface_advance()); // move along ground

		if (val == 2) { // moved, recalculate velocity from position change
			status = 3;
			if (radius >= LARGE_OBJ_RAD) check_vert_collision(obj_index, 1, iter); // adds instability though
			assert(tstep > 0.0);
			if (radius >= LARGE_OBJ_RAD && velocity != zero_vector) {modify_grass_at(pos, radius, 1);} // crush grass
		}
		else if (val == 1) { // stopped
			if (otype.flags & IS_PRECIP) {
				int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

				if (!point_outside_mesh(xpos, ypos) && spillway_matrix[ypos][xpos] >= short(frame_counter-1)) {
					status = 0; // precipitation melts in pool
					return;
				}
			}
			if (status != 4) {
				check_vert_collision(obj_index, 0, iter); // one last time before the object is "stopped"???
				velocity = zero_vector;
				if (!disabled()) status = 4;
			}
		}
		else {
			assert(type != SMILEY);
			status = 0; // bad position
		}
	}
}


int get_obj_zval(point &pt, float &dz, float z_offset) { // 0 = out of bounds/error, 1 = airborne, 2 = on ground

	if (!is_over_mesh(pt))              return 0;
	int const xpos(get_xpos(pt.x)), ypos(get_ypos(pt.y));
	if (point_outside_mesh(xpos, ypos)) return 0;
	if ((pt.z - z_offset) > ztop)       return 1;
	if (is_in_ice(xpos, ypos) && pt.z < water_matrix[ypos][xpos]) return 1;
	float const zval(interpolate_mesh_zval(pt.x, pt.y, 0.0, 0, 0));
	if ((pt.z - z_offset) > zval) return 1;
	dz   = zval - pt.z;
	pt.z = zval + z_offset;
	return 2;
}


int dwobject::object_still_stopped(int obj_index) {

	float const zval(pos.z - object_types[type].radius);
	float const mh(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 0));

	if ((zval - SMALL_NUMBER) <= mh) {
		if (object_types[type].friction_factor >= STICK_THRESHOLD) return 1;
		if (status == 1 || status == 2) status = 3;
		return 0;
	}
	point const old_pos(pos);
	pos.z = zval;
	int const coll(check_vert_collision(obj_index, 0, 0)); // apply coll functions?
	pos   = old_pos;
	if (!disabled() && !coll) status = 1;
	return coll;
}


// 0 = error (bad position), 1 = stopped, 2 = moved
int dwobject::surface_advance() {

	obj_type const &otype(object_types[type]);
	
	if (otype.friction_factor >= STICK_THRESHOLD || (flags & XY_STOPPED)) { // stopped
		velocity = zero_vector;
		return 1;
	}
	int xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y)), val(0);
	if (point_outside_mesh(xpos, ypos)) return 0; // object off edge
	float const mesh_height(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 0)), radius(otype.radius), density(otype.density);

	if (pos.z < (mesh_height - RECOVER_DEPTH*radius)) { // below surface
		if (pos.z < (mesh_height - KILL_DEPTH*radius)) return 0; // far below surface, it's gone
		pos.z = mesh_height; // recover it
	}
	float const grass_friction(0.1*min(1.0f, (grass_length/radius))*get_grass_density(pos)*(is_rain_enabled() ? 0.5 : 1.0));
	float const friction(otype.friction_factor + grass_friction);
	
	if (friction >= STICK_THRESHOLD) { // stopped by grass
		velocity = zero_vector;
		return 1;
	}
	float const s((pos.x - get_xval(xpos))*DX_VAL_INV + 0.5), t((pos.y - get_yval(ypos))*DY_VAL_INV + 0.5);
	int const xpp1(min(xpos+1, MESH_X_SIZE-1)), ypp1(min(ypos+1, MESH_Y_SIZE-1));
	vector3d const &n00(vertex_normals[ypos][xpos]);
	vector3d const &n01(vertex_normals[ypp1][xpos]);
	vector3d const &n10(vertex_normals[ypos][xpp1]);
	vector3d const &n11(vertex_normals[ypp1][xpp1]);
	vector3d const snorm((n11*t + n10*(1.0-t))*s + (n01*t + n00*(1.0-t))*(1.0-s)); // interpolate across the quad
	float const dzn(sqrt(snorm.x*snorm.x + snorm.y*snorm.y));
	vector3d mesh_vel(zero_vector);

	if (dzn > TOLERANCE && dzn > friction) {
		float vel((SURF_ADV_STEP/XY_SCENE_SIZE)*dzn*(1.0 - 0.5*friction)/DEF_TIMESTEP);
		assert(density > 0.0);
		if ((flags & IN_WATER) && density >= WATER_DENSITY) vel *= (density - WATER_DENSITY)/density;

		if (vel > TOLERANCE) {
			mesh_vel.x = vel*DX_VAL*snorm.x/dzn;
			mesh_vel.y = vel*DY_VAL*snorm.y/dzn;
			val        = 1;
		}
	}
	float const vmult((otype.flags & OBJ_IS_DROP) ? 0.0 : pow((1.0f - friction), fticks)); // droplets stick - no momentum
	velocity = (mesh_vel*(1.0 - vmult) + velocity*vmult);
	pos.x   += velocity.x*tstep;
	pos.y   += velocity.y*tstep;
	pos.z    = mesh_height + radius;
	return val+1;
}


// NOTE: norm must be normalized
float get_terrain_rotation(vector3d &axis, point const &pos, vector3d const &norm, vector3d const *const forced_norm) {

	vector3d snorm(plus_z);

	if (forced_norm) {
		snorm = *forced_norm;
	}
	else {
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

		if (point_interior_to_mesh(xpos, ypos)
			&& (h_collision_matrix[ypos][xpos] < mesh_height[ypos][xpos] + SMALL_NUMBER)
			&& !is_mesh_disabled(xpos, ypos) && !mesh_is_underwater(xpos, ypos))
		{
			snorm = surface_normals[ypos][xpos];
		}
	}
	if (snorm == norm) { // normals are parallel
		axis.assign(1.0, 0.0, 0.0);
		return 0.0;
	}
	cross_product(norm, snorm, axis);
	axis.normalize();
	return TO_DEG*get_angle(norm, snorm);
}


void dwobject::set_orient_for_coll(vector3d const *const forced_norm) {

	float const r_angle(init_dir.x);
	vector3d const norm(-sinf(r_angle), cosf(r_angle), 0.0);
	angle = get_terrain_rotation(orientation, pos, norm, forced_norm);
}


int dwobject::check_water_collision(float vz_old) {

	obj_type const &otype(object_types[type]);
	float const radius(otype.radius);
	if ((pos.z - radius) > max_water_height) return 0; // quick check for efficiency
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos))      return 0; // off the mesh
	float water_height;
	vector3d old_v(velocity);
	if (!has_water(xpos, ypos)) return 0; // not over water
	water_height = water_matrix[ypos][xpos];
	if (!is_mesh_disabled(xpos, ypos) && water_height < mesh_height[ypos][xpos]) return 0;
	if (!(flags & IN_WATER) && (pos.z - radius) > water_height)                  return 0; // ???
	if (!is_mesh_disabled(xpos, ypos) && (pos.z + radius + SMALL_NUMBER) < mesh_height[ypos][xpos]) return 0;
	bool const exp_on_coll((otype.flags & EXPL_ON_COLL) != 0);
	
	if (temperature > W_FREEZE_POINT) { // water - object adds to water if precipitation
		assert(type != SMILEY);
		bool const no_water_damage((object_types[type].flags & NO_WATER_DAMAGE) != 0);
		bool splash(0);

		if (!no_water_damage && health <= WATER_DAMAGE) {
			status = 0; // melts
			time   = 2*otype.lifetime;
			old_v  = zero_vector; // actually, old and new velocity are switched
			splash = 1;
		}
		else {
			flags  |= IN_WATER;
			if (!no_water_damage) health -= fticks*WATER_FRAME_DAM;
			float const density(otype.density), v_tot_sq(velocity.mag_sq());
			float const ground_height(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1) + radius);

			if (v_tot_sq < BOUNCE_CUTOFF || (flags & Z_STOPPED)) {
				if (density < WATER_DENSITY || (density == WATER_DENSITY && velocity.z >= 0)) { // floats
					float const zpos(max(water_height + radius*(1.0f - 2.0f*density), ground_height));

					if ((zpos - pos.z) > 2.0*radius) { // under the surface
						velocity.z  = vz_old;
						velocity.z -= ((otype.density - WATER_DENSITY)/otype.density)*base_gravity*GRAVITY*tstep;
						flags      |= Z_STOPPED;
						if ((pos.z - radius) > water_height) splash = 1;
					}
					else {
						pos.z    = zpos + SMALL_NUMBER;
						velocity = zero_vector;
						flags   |= FLOATING;
					}
					if (pos.z > ground_height) {
						status = 1; // floating on water
					}
					else { // floating at edge of water
						if (status == 1 || status == 2) {
							status = 3; // still moving
							return 0;
						}
						else status = 4; // stationary
					}
				}
				else if (v_tot_sq < BOUNCE_CUTOFF) { // sinks
					if (pos.z > (ground_height + 0.00001)) { // sinking in water
						velocity *= density/(density + WATER_DENSITY);
					}
					else { // stationary/moving under water
						pos.z    = ground_height;
						velocity = zero_vector;
						status   = 3;
						return 0;
					}
				}
			}
			else { // collision with water
				bool const surf_coll(pos.z > (water_height - radius - MAX_SPLASH_DEPTH) && !(flags & (UNDERWATER | Z_STOPPED)));
				float const coll_angle(surf_coll ? get_water_coll_angle(velocity) : 0.0);
				float const den_ratio(WATER_DENSITY/(density + WATER_DENSITY));
				float const zpos(max(water_height + radius*(1.0f - 2.0f*density), ground_height));
				if (surf_coll) splash = 1; // first time
				vector3d norm(wat_vert_normals[ypos][xpos]);
				
				if (coll_angle < CRITICAL_ANGLE/den_ratio || (zpos - pos.z) > 6.0*radius || !object_bounce(2, norm, 0.0, 0.0)) {
					// object enters water
					velocity *= (1.0 - WATER_DAMPING*den_ratio);
			
					if (density >= WATER_DENSITY) {
						flags |= UNDERWATER; // will stay below the surface (what about water evaporating?)

						if (pos.z <= (ground_height + 0.00001)) {
							status = 3;
							return 0;
						}
					}
					else {
						velocity.z = 0.0;
						flags     |= (Z_STOPPED | FLOATING);
					}
				}
			}
		}
		if (splash) { // hit the pool of water
			float energy(get_coll_energy(old_v, (exp_on_coll ? zero_vector : velocity), otype.mass));

			if (energy > 0.0) {
				draw_splash(pos.x, pos.y, water_height, SPLASH_BASE_SZ*sqrt(energy));
				
				if (type != DROPLET) {
					if (type == SHRAPNEL) {
						if (rand()%10 < 6) {energy = 0.0;} else {energy *= 0.2;}
					}
					//else if (type == FRAGMENT) {energy *= 0.2;} // too many fragments adding energy gives too large of a splash
					if (energy > 0.0) {add_splash(pos, xpos, ypos, energy, radius, (radius >= LARGE_OBJ_RAD));}
				}
			}
		}
	}
	else { // stuck in ice
		do_coll_damage();
		if (status == 0) return 1;

		if (status != 4 && !(flags & IN_WATER)) {
			if (time <= 1) { // started out under the ice
				velocity = zero_vector;
				flags   |= (XYZ_STOPPED | IN_WATER);
				status   = 4;
			}
			else {
				pos.z = water_height + radius + SMALL_NUMBER; // sitting on ice surface
			}
		}
		vector3d norm(plus_z);

		if (otype.elasticity < ICE_BOUNCE_ELAS || !object_bounce(1, norm, 0.0, 0.0)) {
			velocity = zero_vector; // stuck to ice
			flags   |= (XYZ_STOPPED | IS_ON_ICE);
			status   = 4;
		}
	}
	if (exp_on_coll) disable();
	return 1;
}


void dwobject::surf_collide_obj() const {

	switch (type) {
	case CHUNK:
		if (flags & TYPE_FLAG) break; // charred, not blood
	case BLOOD:
		if (snow_height(pos)) { // in the snow
			add_color_to_landscape_texture(BLOOD_C, pos.x, pos.y, ((type == BLOOD) ? 4.0 : 2.2)*object_types[type].radius);
		}
		break;
	}
}


void dwobject::elastic_collision(point const &obj_pos, float energy, int obj_type) {

	if (disabled() || (object_types[type].flags & COLL_DESTROYS)) return; // self-propelled
	if (temperature <= W_FREEZE_POINT && (flags & IN_WATER))      return; // stuck in ice
	vector3d const vdir(pos, obj_pos);
	//float const elastic(object_types[otype].elasticity*object_types[obj_type].elasticity);
	float const elastic(object_types[type].elasticity), vdir_mag(vdir.mag());
	float const vmag(sqrt(2.0*elastic*energy/object_types[type].mass)); // E = 0.5*M*dV^2 => dV = sqrt(2*E/M)
	if (vdir_mag > TOLERANCE) velocity += vdir*(vmag/vdir_mag);
	status = 1; // re-animate
	flags &= ~ALL_COLL_STOPPED;
}


void reanimate_objects() {

	for (int i = 0; i < num_groups; ++i) {
		obj_group &objg(obj_groups[i]);
		if (!objg.enabled) continue;
		if (objg.type == SMILEY) continue; // doesn't apply to smileys, since they move independently from the physics system

		for (unsigned j = 0; j < objg.end_id; ++j) {
			dwobject &obj(objg.get_obj(j));
			if (obj.disabled()) continue;
			int const xpos(get_xpos(obj.pos.x)), ypos(get_ypos(obj.pos.y));
			int cindex; // unused
			float const radius(object_types[obj.type].radius);

			if (!(point_outside_mesh(xpos, ypos)) && !check_legal_move(xpos, ypos, obj.pos.z, radius, cindex)) {
				obj.status = 0; // within collision object
			}
			else {
				obj.flags = 0; // no longer stopped in x and y
				if (obj.status == 2 || obj.status == 4) obj.status = 1; // stationary => moving in air
			}
		}
	}
}


void seed_water_on_mesh(float amount) {

	for (int y = 0; y < MESH_Y_SIZE; ++y) {
		for (int x = 0; x < MESH_X_SIZE; ++x) {
			if (wminside[y][x] != 1) continue;
			int const wsi(watershed_matrix[y][x].wsi);
			assert(wsi < (int)valleys.size());
			valleys[wsi].w_volume += amount;
			has_accumulation = 1; // object inside volume
		}
	}
}


void accumulate_object(point const &pos, int type, float amount) {

	if (pos.x < -X_SCENE_SIZE+0.1*DX_VAL || pos.x > X_SCENE_SIZE-0.1*DX_VAL ||
		pos.y < -Y_SCENE_SIZE+0.1*DY_VAL || pos.y > Y_SCENE_SIZE-0.1*DY_VAL) return;
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos)) return;

	if (temperature > W_FREEZE_POINT) {
		if (!DISABLE_WATER && (display_mode & 0x04) && wminside[ypos][xpos] == 1) {
			int const wsi(watershed_matrix[ypos][xpos].wsi);
			assert(wsi < (int)valleys.size());
			float const wvol(valleys[wsi].w_volume);
			if (wvol >= 0) {valleys[wsi].blood_mix = (valleys[wsi].blood_mix*wvol + amount*(type == BLOOD))/(wvol + amount);}
			valleys[wsi].w_volume += amount;
			has_accumulation = 1; // object inside volume
		}
	}
	else if (type == SNOW) {
		float const acc(SNOW_ACC*amount*(1.0 + rand_float()));
		accumulation_matrix[ypos][xpos] += acc;
		add_snow_to_landscape_texture(pos, acc);
		has_snow_accum  |= (acc > 0.0);
		has_accumulation = 1; // object inside volume
	}
}


int dwobject::object_bounce(int coll_type, vector3d &norm, float elasticity2, float z_offset, vector3d const &obj_vel) {

	float elasticity(object_types[type].elasticity);
	if (elasticity == 0.0)      return 0;
	vector3d const delta_v(velocity - obj_vel);
	if (delta_v == zero_vector) return 0;
	vector3d const v_orig(velocity);
	vector3d bounce_v;
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

	if (point_outside_mesh(xpos, ypos)) { // object on/over edge
		status = 0;
		return 0;
	}
	switch (coll_type) {
	case 0: // mesh surface
		{
			float const zval(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 0));
			if ((pos.z - z_offset) < zval) pos.z = zval + z_offset;
			norm = surface_normals[ypos][xpos];
			elasticity *= LAND_ELASTICITY*(1.0 - 0.5*get_grass_density(pos)); // half elastic in dense grass
			if (spillway_matrix[ypos][xpos] >= short(frame_counter-1)) elasticity *= SPILL_ELASTIC;
		}
		break;
	case 1: // ice
		norm.assign(0.0, 0.0, -1.0);
		elasticity *= ICE_ELASTICITY;
		break;
	case 2: // water
		norm.assign(0.0, 0.0, -1.0);
		elasticity *= WATER_ELASTIC;
		break;
	case 3:
	default: // horizontal/vertical surface or other
		elasticity *= elasticity2;
		assert(norm != zero_vector); // norm must come in correct
	}
	elasticity = CLIP_TO_01(elasticity);
	assert(!is_nan(norm));
	calc_reflection_angle(delta_v, bounce_v, norm);
	assert(!is_nan(bounce_v));
	float const xy_elasticity(elasticity*(1.0 - object_types[type].air_factor));
	velocity.assign(xy_elasticity*bounce_v.x, xy_elasticity*bounce_v.y, elasticity*bounce_v.z);
	float const v_tot_sq(velocity.mag_sq());

	if (v_tot_sq >= BOUNCE_CUTOFF || type == DYNAM_PART) {
		if (type == PLASMA && (coll_type == 0 || coll_type == 3) && v_tot_sq >= 2.25*BOUNCE_CUTOFF && (rand()%10) < 8) {
			gen_fire(pos, rand_uniform(0.4, 1.2), source);
		}
		return 1;
	}
	velocity = v_orig;
	return 0;
}


void bubble::apply_physics(unsigned i) {

	if (!status) return;
	
	if (temperature <= W_FREEZE_POINT) {
		status = 0; // frozen
		return;
	}
	time  += iticks;
	pos.z += tstep*velocity;
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

	if (point_outside_mesh(xpos, ypos)) {
		status = 0; // out of simulation region
	}
	else if (pos.z >= water_matrix[ypos][xpos]) {
		draw_splash(pos.x, pos.y, water_matrix[ypos][xpos], 2.0*radius, color);
		status = 0; // pops at water surface
	}
}


void particle_cloud::apply_physics(unsigned i) {

	unsigned const num_smoke_advance(5); // must be > 0
	if (status == 0) return;

	if (pos.z >= (CLOUD_CEILING + zmax_est) || radius > MAX_PART_CLOUD_RAD || is_underwater(pos)) { // smoke dies
		status = 0;
		return;
	}
	vector3d v_flow(get_flow_velocity(pos));
	int coll(0);
	dwobject obj(SMOKE, pos, zero_vector, 1, 10000.0); // make a SMOKE object for collision detection
	object_types[SMOKE].radius = radius;
	unsigned const steps((czmax > -FAR_DISTANCE && pos.z > czmax) ? 1 : num_smoke_advance);

	for (unsigned j = 0; j < steps; ++j) {
		vector3d vel((get_local_wind(obj.pos) + v_flow)*0.5);
		vel.z   *= 0.5;
		obj.pos += (vel + init_vel)*(tstep/(double)num_smoke_advance);
		vector3d cnorm;
		
		if (obj.check_vert_collision(0, 0, j, &cnorm, all_zeros, 1, 1)) { // skip dynamic, only_drawn
			// destroy the smoke if it's not damaging and hits the bottom of a static drawn object (excludes trees and scenery)
			if (cnorm.z < 0.0 && damage == 0.0) { // <= 0.0?
				if (acc_smoke && time > 0) add_smoke(pos, 1.0);
				status = 0;
				return;
			}
			coll = 1;
			break;
		}
	}
	float const tstep_scale(TIMESTEP/DEF_TIMESTEP), rscale(get_rscale());
	pos       = obj.pos;
	time     += iticks;
	density  *= pow(0.97f, tstep_scale);
	darkness *= pow(0.98f, tstep_scale);
	radius   *= pow(1.03f, tstep_scale);
	if (density  < 0.0001) density  = 0.0;
	if (darkness < 0.0001) darkness = 0.0;
	
	if (damage > 0.0) {
		if (is_fire()) {modify_grass_at(pos, radius, 0, 1);} // burn grass
		do_area_effect_damage(pos, radius, damage*rscale, i, source, damage_type);
	}
	if (is_fire()) {
		colorRGBA color(base_color);
		color.G *= rscale;
		add_dynamic_light(3*radius, pos, color);
		if (coll && radius >= MAX_PART_CLOUD_RAD && (rand()&7) == 0) gen_fire(pos, 1.0, source); // will be destoyed next frame
	}
	if (damage_type == GASSED) { // check for gas ignition near fire
		for (unsigned i = 0; i < fires.size(); ++i) {
			if (fires[i].status != 0 && dist_less_than(fires[i].pos, pos, radius)) {
				create_explosion(pos, source, 0, 10*damage*rscale, 4*radius, BLAST_RADIUS, 0);
				status = 0;
				return;
			}
		}
	}
}


void fire::apply_physics(unsigned i) {

	if (status == 0) return;
	assert(radius > 0.0);
	float const damage(0.5*heat*radius);
	do_area_effect_damage(pos, 2.0*radius, BURN_DAMAGE*damage, i, source, FIRE);
	colorRGBA const fcolor(gen_fire_color(cval, inten));
	if (damage > 0.001) add_dynamic_light(64.0*inten*damage, pos, fcolor, plus_z, light_bwidth);

	if (!is_static) {
		int const rn(max(1, int(8.0 + 0.02/(0.1 + sqrt(radius*sqrt(heat))))));
		if (rand()%rn == 0) gen_smoke(pos);
	}
	if (animate2) time += iticks;
	point pos2(pos);
	pos2.z -= radius;
	bool underwater(is_underwater(pos2));

	if (underwater && temperature <= W_FREEZE_POINT) { // on/under ice
		radius    *= pow(0.95f, fticks); // slowly die out
		underwater = 0;
	}
	if ((!is_static && time > (int)MAX_FIRE_TIME) || radius < TOLERANCE || underwater) {
		extinguish();
		return;
	}
	if (!is_static) {
		object_types[FIRE].radius = 1.75*radius;
		//destroy_coll_objs(pos, 100000.0*damage, source, FIRE); // very, very slow

		if (status == 2) {
			velocity = zero_vector; // above ground - no movement
			float const zval(interpolate_mesh_zval(pos.x, pos.y, 0.0/*radius*/, 0, 0));

			if (pos.z - zval > 1.5*radius) {
				dwobject obj(FIRE, pos, zero_vector, 1, 10000.0); // make a FIRE object for collision detection
				obj.source = source;

				if (!obj.check_vert_collision(i, 1, 0)) {
					pos.z -= radius;
					status = 1; // re-animate
				}
			}
		}
		else {
			point const lpos(pos);
			vector3d const local_wind(get_local_wind(pos));
			vector3d const vel((local_wind.x + rand_uniform(-1.5, 1.5)), (local_wind.y + rand_uniform(-1.5, 1.5)), rand_uniform(-0.05, 0.0585));
			velocity *= pow(0.95f, fticks);
			velocity += vel*(0.005*tstep);
			pos.x    += fticks*velocity.x;
			pos.y    += fticks*velocity.y;
			set_true_obj_height(pos, lpos, FAR_DISTANCE, velocity.z, FIRE, 0, 0, 0, 1);
			pos.z    -= radius;
			pos.z     = 0.9*lpos.z + 0.1*pos.z; // slow movement
			//pos.z     = interpolate_mesh_zval(pos.x, pos.y, radius, 0, 0) + 0.6*radius;
		}
		radius += (0.02 + radius)*(rand_uniform(-0.02, 0.02) + 250.0*velocity.z);
		heat    = 0.8*heat + 0.2*rand_uniform(0.25, 1.2)/(0.9 + 2.0*radius);
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

		if (radius <= 0.0001 || !point_interior_to_mesh(xpos, ypos)) {
			extinguish();
			return;
		}
		if (status != 2) {surface_damage[ypos][xpos] += 20.0*radius*heat;} // near mesh

		if (radius > 0.04) { // split into smaller fires
			for (unsigned i = 0; i < 2; ++i) {pos2[i] = pos[i] + rand_uniform(-0.05, 0.05);}
			if (rand() & 1) {pos2.z = pos.z = interpolate_mesh_zval(pos.x, pos.y, radius, 0, 0) + 0.3*radius;}
			else {pos2.z = pos.z + rand_uniform(-0.05, 0.05);}
			gen_fire(pos2, 1.0, source, 1);
			radius -= 0.017;
		}
	} // !is_static
	if (animate2 && damage > 0.005 && (rand()%max(1, int(0.5/damage))) == 0) {gen_particles(pos, 1);}
}

void fire::extinguish() {

	status = 0;
	gen_smoke(pos + point(0.0, 0.0, radius));
}


void decal_obj::apply_physics(unsigned i) {

	if (!status) return;
	time += iticks;
	if (time > lifetime) {status = 0;}
}

float decal_obj::get_alpha() const {
	return alpha*CLIP_TO_01(2.0f - 2.0f*float(time)/float(lifetime)); // first half alpha=1, second half fade to 0
}


bool physics_particle_manager::is_pos_valid(point const &pos) const {

	if (!is_over_mesh(pos))    return 0; // outside simulation region
	if (pos.z < water_plane_z) return 0; // underwater
	if (!dist_less_than(pos, get_camera_pos(), 8.0)) return 0; // too far away
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos)) return 0; // can this fail?
	return (pos.z > mesh_height[ypos][xpos] && pos.z > water_matrix[ypos][xpos]); // above mesh and water
}

void physics_particle_manager::apply_physics(float gravity, float terminal_velocity) {

	if (parts.empty()) return;
	//RESET_TIME;
	unsigned o(0);
	float const g_acc(base_gravity*GRAVITY*tstep*gravity), xy_damp(pow(0.98f, fticks));

	for (unsigned i = 0; i < parts.size(); ++i) {
		part_t &part(parts[i]);
		part.v.z  = max(-terminal_velocity, (part.v.z - g_acc)); // apply gravity + terminal velocity
		part.v.x *= xy_damp;
		part.v.y *= xy_damp;
		//point const p0(part.p);
		part.p   += tstep*part.v; // add velocity to position
		int cindex;
		
		//if (check_coll_line(p0, part.p, cindex, -1, 1, 0)) { // skip dynamic
		if (check_point_contained_tree(part.p, cindex, 0)) { // skip dynamic
			continue; // destroy particle, don't bounce
			//part.p = cpos; vector3d bounce_v; calc_reflection_angle(part.v.get_norm(), bounce_v, cnorm); part.v = 0.9*part.v.mag()*bounce_v;
		}
		if (is_pos_valid(part.p)) {parts[o++] = part;} // above water and mesh - copy/compact
	}
	parts.resize(o);
	//PRINT_TIME("Particle Physics"); // 0.07ms average / 0.24ms with collisions
}

void water_particle_manager::apply_physics() {
	physics_particle_manager::apply_physics(object_types[DROPLET].gravity, object_types[DROPLET].terminal_vel);
}


template<typename T> void shift_objs(vector<T> &objs, vector3d const &vd) {

	for (unsigned i = 0; i < objs.size(); ++i) {
		if (objs[i].status) {objs[i].pos += vd;}
	}
}

template<typename T> void apply_obj_physics(vector<T> &objs) {
	for (unsigned i = 0; i < objs.size(); ++i) {objs[i].apply_physics(i);}
}


void shift_other_objs(vector3d const &vd) {

	shift_objs(bubbles,     vd);
	shift_objs(part_clouds, vd);
	shift_objs(fires,       vd);
	shift_objs(decals,      vd);
}


void advance_physics_objects() {

	apply_obj_physics(bubbles);
	apply_obj_physics(part_clouds);
	apply_obj_physics(fires);
	apply_obj_physics(decals);
	explosion_part_man.apply_physics(1.0, 4.0); // gravity=1.0, air_factor=0.25
	water_part_man.apply_physics();
	for (unsigned i = 0; i < decals.size(); ++i) {decals[i].check_cobj();}
}


void reset_other_objects_status() {

	reset_status(bubbles);
	reset_status(part_clouds);
	reset_status(fires);
	reset_status(decals);
}


// Note: In combined mode, times should be configured for planet/moon rotation and revolution periods
void auto_advance_time() { // T = 1 hour

	no_sun_lpos_update = 0;
	if (!auto_time_adv || !animate2) return;
	static int last_itime(0), precip_inited(0);
	static float ftime, precip(0.0), precip_app_rate(0.0), prate(0.0), ccover(0.0);
	static rand_gen_t rgen;
	if (last_itime == 0) {rgen.set_state(rand(), rand());} // random seed at the start

	// auto_time_adv = 0 => no change
	// auto_time_adv = 1 => 1 hour = 60 min = 1 hour
	// auto_time_adv = 2 => 1 hour = 6 min
	// auto_time_adv = 3 => 1 hour = 36 sec
	// auto_time_adv = 4 => 1 hour = 3.6 sec
	float const delta_time(fticks*pow(10.0, auto_time_adv-1)/56000.0);
	ftime += delta_time;
	int const itime((int)ftime); // half hours
	int const hrtime(itime/2 + 12), hrtime24(hrtime%24); // hours
	int const date(hrtime/24); // days
	int vis_recalc(0);

	// move sun (daily: dt=24, every t)
	sun_rot = (PI/24)*((itime%48) + (ftime - itime)); // 0 at 12:00 noon
	update_sun_and_moon();
	if (light_factor >= 0.4) {vis_recalc |= SUN_SHADOW;}

	// change wind
	for (unsigned d = 0; d < 3; ++d) {
		wind[d] += 1.0*delta_time*rgen.signed_rand_float();
		wind[d]  = CLIP_TO_pm1(wind[d]);
	}
	wind.z *= 0.99;
	
	if (itime == last_itime) {
		no_sun_lpos_update = 1;
		return;
	}
	last_itime = itime;

	// move moon (28 days: dt = 24*28 = 672, every 24t)
	if (hrtime24 == 12) {
		moon_rot += PI/14.0;
		if (moon_rot > TWO_PI) {moon_rot -= TWO_PI;}
		update_sun_and_moon();
		if (light_factor <= 0.6) {vis_recalc |= MOON_SHADOW;}
	}

	// change cloudiness (7 days: dt = 24*7 = 168, every t)
	ccover += 0.1*rgen.signed_rand_float();
	if (ccover < -0.2) {ccover = 0.2 - ccover;}
	cloud_cover = CLIP_TO_01(0.5f*max(-0.1f, precip) + 0.5f*min(1.0f, ccover));

	// change temperature (daily: dt = 24, seasonal, dt = 24*365.25 = 8766, every t)
	float const delta_t(MAX_TEMP - MIN_TEMP);
	float const random_temp(max(MIN_TEMP, min(MAX_TEMP, (temperature + 0.05f*delta_t*rgen.signed_rand_float()))));
	float const cloud_temp(MIN_TEMP + delta_t*(1.0 - cloud_cover));
	float const seasonal_temp(MIN_TEMP + delta_t*2.0*(0.5 - (float((hrtime24 + 4383)%8766))/8766.0));
	float const time_of_day_temp(MIN_TEMP + delta_t*(fabs(PI - sun_rot)/PI));
	temperature      = 0.22*random_temp + 0.1*cloud_temp + 0.28*seasonal_temp + 0.17*time_of_day_temp + 0.23*alt_temp;
	temperature      = min(MAX_TEMP, max(MIN_TEMP, temperature));
	//cout << "rtemp = " << random_temp << ", ctemp = " << cloud_temp << ", stemp = " << seasonal_temp << ", todtemp = " << time_of_day_temp << ", atemp = " << alt_temp << ", temp = " << temperature << endl;
	int const cid(coll_id[PRECIP]);

	// change precipitation (7 days: dt = 24*7 = 168, every t)
	if (!precip_inited) {
		precip_app_rate = (float)obj_groups[cid].app_rate;
		precip_inited   = 1;
	}
	prate += 0.2*rgen.signed_rand_float();
	if      (prate < -0.5) {prate = -1.0 - prate;}
	else if (prate >  0.5) {prate =  1.0 - prate;}
	precip      = max(-2.0f, min(1.0f, (precip + 0.4f*prate)));
	precip_mode = ((precip > 0.0) ? 1 : 0);
	obj_groups[cid].app_rate = ((precip_mode > 0) ? max(0, int(2.0*precip*precip_app_rate)) : 0.0);

	// change leaf color (seasonal, dt = 24*365.25 = 8766, every 168t)
	if (hrtime%168 == 0) {} // *** WRITE ***
	if (vis_recalc) {check_update_global_lighting(vis_recalc);}

	if (PRINT_TIME_OF_DAY) {
		if (hrtime24 < 12) {
			cout << "Time = " << (hrtime24==0  ? 12 : hrtime24   ) << ":" << (((itime&1) == 0) ? "00" : "30") << " AM";
		}
		else {
			cout << "Time = " << (hrtime24==12 ? 12 : hrtime24-12) << ":" << (((itime&1) == 0) ? "00" : "30") << " PM";
		}
		printf(", Day=%i, Temp=%.3f, Clouds=%.3f, Precip=%.3f/%u/%.3f.\n",
			date, temperature, cloud_cover, precip, obj_groups[cid].app_rate, prate);
		cout << "Wind = " << wind.x << ", " << wind.y << ", " << wind.z << "." << endl;
	}
}




