// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 11/2/03

#ifndef _GAMEPLAY_H_
#define _GAMEPLAY_H_


#include "3DWorld.h"
#include "mesh.h"
#include "player_state.h"


float const SSTEPS_PER_FRAME0     = 7.0; // 7.0-8.0: smileys about the same speed as player
unsigned const SMILEY_MAX_TRIES   = 100;
float const SMILEY_DIR_FACTOR     = 0.05;
float const S_SH_SCALE            = 2.0;
unsigned const HIT_TIME           = 6;

float const HEALTH_PER_DAMAGE     = 0.1;

unsigned const NUM_BS             = 100;
unsigned const BLAST_CHAIN_DELAY  = 3;

float const FALL_HURT_VEL         = 4.0;
float const FALL_HURT_HEIGHT      = 4.0;
float const MAX_HEALTH            = 150.0;
float const MAX_SHIELDS           = 150.0;
float const MAX_REGEN_HEALTH      = 200.0;
float const KILL_HEALTH           = 5.0;
float const INIT_SHIELDS          = 10.0;
unsigned const CAMERA_SPHERE_TIME = 4;

float const PLASMA_SIZE_INCREASE  = 0.004;
float const MAX_PLASMA_SIZE       = 3.5;
float const PLASMA_LT_DAMAGE      = 18.0;
float const BURN_DAMAGE           = 1200.0;
float const PGAS_DAMAGE           = 1.0;

unsigned const SMILEY_TEX_SIZE    = 256;

unsigned const MESSAGE_TIME       = unsigned(2 *TICKS_PER_SECOND);
unsigned const DROWN_TIME         = unsigned(10*TICKS_PER_SECOND);
unsigned const INIT_PU_SH_TIME    = unsigned(4 *TICKS_PER_SECOND);
unsigned const LM_ACT_TIME        = unsigned(2 *TICKS_PER_SECOND);
unsigned const SMILEY_LM_ACT_TIME = unsigned(8 *TICKS_PER_SECOND);


struct blood_spot { // size = 20

	int time;
	float size;
	point pos;

	blood_spot() : time(0), size(0.0), pos(all_zeros) {}
};


struct weapon { // size = 64

	bool self_damage, need_weapon, need_ammo;
	unsigned def_ammo, max_ammo, obj_id, fire_delay, nshots, nfragments;
	float v_mult, v_add, blast_damage, blast_radius, firing_error, range;
	string name;

	weapon(bool sd, bool nw, bool na, unsigned da, unsigned ma, unsigned oi, unsigned fd, unsigned ns, unsigned nf,
		float vm, float va, float bd, float br, float fe, float ra, string const &name_)
		: self_damage(sd), need_weapon(nw), need_ammo(na), def_ammo(da), max_ammo(ma), obj_id(oi), fire_delay(fd), nshots(ns),
		nfragments(nf), v_mult(vm), v_add(va), blast_damage(bd), blast_radius(br), firing_error(fe), range(ra), name(name_) {}
};


struct laser_beam { // size = 32

	int distant, shooter;
	float intensity;
	point pts[2];
	colorRGBA color;
	
	laser_beam(int dist, int shoot, point const &pt0, point const &pt1, colorRGBA const &c, float int_=1.0)
		: distant(dist), shooter(shoot), intensity(int_), color(c)
	{
		pts[0] = pt0;
		pts[1] = pt1;
	}
	void draw() const;
};


unsigned const UNDEF      = (unsigned)-1;
unsigned const CBFD       = 60;
float const CBLADE_EXT_PT = 0.04;
float const CBLADE_EXT(0.5*CBLADE_EXT_PT*CBFD);


// self_damage need_weapon need_ammo def_ammo max_ammo obj_id fire_delay nshots nfragments v_mult v_add blast_damage blast_radius firing_error range name
weapon const weapons[NUM_WEAPONS+1] = {
	weapon(0, 0, 0, 0,   0,   UNDEF,    0,   0,  0,   0.0,  0.0, 0.0,    0.0,  0.0,   0.0,  "Unarmed"         ),
	weapon(0, 0, 0, 0,   0,   UNDEF,    23,  1,  1,   0.0,  0.0, 500.0,  0.25, 0.0,   0.25, "Baseball Bat"    ),
	weapon(0, 0, 1, 1,   3,   BALL,     25,  1,  1,   1.5,  3.0, 0.0,    0.0,  0.05,  3.0,  "Dodgeball"       ),
	weapon(0, 0, 1, 30,  500, S_BALL,   18,  1,  1,   1.3,  3.3, 0.0,    0.0,  0.0,   1.5,  "Bouncy Ball"     ),
	weapon(1, 1, 1, 10,  100, ROCKET,   32,  1,  1,   0.7,  3.1, 1000.0, 0.42, 0.003, 0.0,  "Rocket Launcher" ),
	weapon(1, 0, 1, 5,   50,  LANDMINE, 30,  1,  1,   0.0,  2.0, 4000.0, 0.39, 0.0,   4.0,  "Proximity Mine"  ),
	weapon(1, 1, 1, 5,   50,  SEEK_D,   60,  1,  1,   0.5,  2.5, 2300.0, 0.50, 0.0,   0.0,  "Seek and Destroy"),
	weapon(0, 0, 1, 25,  500, STAR5,    10,  1,  1,   1.1,  3.0, 0.0,    0.0,  0.015, 2.0,  "Throwing Star"   ),
	weapon(0, 1, 1, 100, 600, UNDEF,    2,   1,  1,   0.0,  0.0, 70.0,   0.0,  0.02,  0.0,  "M16"             ),
	weapon(0, 1, 1, 12,  100, UNDEF,    27,  24, 1,   0.0,  0.0, 50.0,   0.0,  0.08,  5.0,  "Shotgun"         ),
	weapon(1, 0, 1, 12,  60,  GRENADE,  22,  1,  140, 1.0,  1.2, 700.0,  0.44, 0.01,  1.7,  "Grenade"         ),
	weapon(0, 1, 1, 200, 800, UNDEF,    1,   1,  1,   0.0,  0.0, 16.0,   0.0,  0.0,   0.0,  "Laser"           ),
	weapon(1, 1, 1, 20,  200, PLASMA,   13,  1,  1,   1.4,  3.5, 200.0,  0.43, 0.0,   3.8,  "Plasma Cannon"   ),
	weapon(0, 1, 0, 1,   10,  UNDEF,    CBFD,1,  1,   1.5,  4.0, 40.0,   0.2,  0.0, CBLADE_EXT, "Carnage Blade"),
	weapon(0, 1, 1, 50,  250, GASSED,   4,   1,  1,   1.2,  2.8, 100.0,  0.07, 0.1,   2.8,  "Gasser"),
	/* non-selectable */
	weapon(1, 0, 1, 3,   20,  CGRENADE, 80,  1,  8,   0.9,  1.1, 800.0,  0.45, 0.02,  1.8,  "Cluster Grenade" )
};

int const obj_weapons[NUM_TOT_OBJS] = {
	-1, -1, -1, -1, W_BALL, W_SBALL, -1, -1, -1, -1,
	-1, W_ROCKET, W_LANDMINE, W_SEEK_D, W_STAR5, W_PLASMA, W_GRENADE, W_CGRENADE, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1
};

bool const damage_done_obj[NUM_TOT_OBJS] = {
	0,0,0,0,1,1,0,0,0,0,
	0,1,0,1,1,1,1,1,1,0,
	0,0,0,0,0,1,1,0,1,1,
	1,1,1,1,0,0,1,1,1,1,
	1,1,1,1,1,1,1,1,1,1,
	0,0,0,0,0
};

string const obj_type_names[NUM_TOT_OBJS] = {
	"Rain", "Snow", "Hail", "Leaf", "Dodgeball", "Bouncy Ball", "Smiley", "Blood", "Charred Smiley", "Chunk",
	"Smiley Part", "Rocket", "Landmine", "Seek and Destroy", "Throwing Star", "Plasma", "Grenade", "Cluster Grenade", "Shrapnel", "Shell Casing",
	"Projectile Hit", "Droplet", "Water Droplet", "Sand", "Dirt", "Rock", "Fragment", "Particle", "Health", "Shields",
	"Powerup", "Weapon", "Ammo", "Pack", "Camera", "Precipitation", "Blast Radius", "Projectile", "laser beam", "Impact",
	"Plasma Lightning Damage", "Laser", "Drowned", "Burned", "Fire", "Fell", "Froze", "Suffocated", "Crushed", "Poison Gas",
	"Waypoint", "Smoke", "Dynamic Particle", "Skull", "Grass"
};

string const powerup_names[NUM_POWERUPS] =
	{"Quad Damage", "Regeneration", "Shielding", "Haste", "Flight", "Invisibility"};


point get_sstate_pos(int id);
vector3d get_sstate_dir(int id);

int same_team(int source, int target);
void drop_pack(player_state &sstate, point const &pos);
string get_weapon_qualifier(int type, int index, int source);
void gen_blood_velocity(vector3d &vout, vector3d const &velocity, vector3d const &coll_dir, float blood_v, float md, float mv, int type, float health);
int  damage_done(int type, int index);
void blood_on_camera(unsigned num_spots);
void init_sstate(int id, bool w_start);
int  get_range_to_mesh(point const &pos, vector3d const &vcf, point &coll_pos, int &xpos, int &ypos);
point projectile_test(point const &pos, vector3d const &vcf_, float firing_error, float damage,
					  int shooter, float &range, float intensity=1.0, int ignore_cobj=-1);
int fire_projectile(point fpos, vector3d dir, int shooter, int &chosen_obj);
float get_projectile_range(point const &pos, vector3d vcf, float dist, float range, point &coll_pos, vector3d &cnorm,
						   int &coll, int &cindex, int source, int check_splash, int ignore_cobj=-1);
void draw_plasmaball(point const &pos0, int shooter);
void init_smiley(int smiley_id);
void init_smiley_weapon(int smiley_id);
void smiley_fire_weapon(int smiley_id);
int  get_damage_source(int type, int index, int questioner);
void gen_rubble(int type, int num, point const &pos, int shooter, float const p[7]);

void add_damage_to_smiley(vector3d const &dir, float size, int smiley_id, int type);
void draw_plasma(point const &pos, point const &part_pos, float radius, float size, int ndiv, int shpere_tex, bool gen_parts, int time);
void do_cblade_damage_and_update_pos(point &pos, int shooter);


#endif


