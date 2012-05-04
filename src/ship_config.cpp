// 3D World - Ship models for universe mode
// by Frank Gennari
// 8/25/05

#include "ship.h"
#include "ship_util.h"
#include "explosion.h"


bool const SHOW_SHIP_RATINGS = 0;


bool allow_add_ship(0), regen_uses_credits(0), respawn_req_hw(0), player_enemy(0), build_any(0);
float spawn_dist(1.0), global_regen(0.0), hyperspeed_mult(200.0), ship_speed_scale(1.0), player_turn_rate(1.0);
unsigned gen_counts  [NUM_ALIGNMENT] = {0};
point ustart_pos(all_zeros);
u_ship *player_ship_ptr = NULL; // easiest to just make this global
unsigned ship_add_prob[2][NUM_ALIGNMENT][NUM_US_CLASS] = {0};
vector<unsigned> build_types[NUM_ALIGNMENT];
vector<string> ship_names;
vector<us_class> sclasses;
vector<us_weapon> us_weapons;
vector<ship_weapon> player_init_weapons;
map<unsigned, unsigned> sclasses_to_weapons;


extern int do_run;
extern unsigned team_credits[], init_credits[], alloced_fobjs[];
extern point player_death_pos, universe_origin, last_camera;
extern char *ship_def_file;


void init_ship_weapon_classes();


// ************ setup and helper functions ************


bool create_player_ship(unsigned sclass, unsigned align) {

	if (player_ship_ptr != NULL) return 0;
	player_ship_ptr = create_ship(sclass, ustart_pos, align, AI_ATT_ENEMY, TARGET_LAST, 0);
	player_ship_ptr->rename("Player");
	change_speed_mode(do_run);
	return 1;
}


void reset_player_ship() {

	player_ship().weapons = player_init_weapons;
}


void setup_ships() { // sets up all ship and weapons classes and default objects as well as creating the player's ship

	static bool inited(0);
	if (inited) return;
	inited = 1;
	cout << "Initializing ships..."; // should only happen once
	cout.flush();
	init_ship_weapon_classes(); // maybe should be somewhere more global
	cout << "Done." << endl;
}


void add_player_weap(unsigned weapon, unsigned num=1, unsigned ammo=0) {

	player_ship().add_weapon(ship_weapon(weapon, num, ammo));
}


void add_ship_weapon(unsigned sclass, unsigned weapon, unsigned num, unsigned ammo, vector<point> const &weap_pts) {

	assert(sclass < sclasses.size());
	sclasses[sclass].add_weapon(ship_weapon(weapon, num, ammo, weap_pts));
}


u_ship *add_ship(unsigned sclass, unsigned align, unsigned ai, unsigned targ, point const &pos, float spread) {

	// check for p2p_dist(pos, get_player_pos2()) > thresh ?
	assert(sclass < sclasses.size() && align < NUM_ALIGNMENT && (ai & AI_BASE_TYPE) < NUM_AI && targ < NUM_TARGET);
	assert(!sclasses[sclass].orbiting_dock); // can't add these - must be attached to a planetary body
	unsigned coll_ix(0);
	point spos;
	unsigned const MAXITER(10);

	// check for collisions (invalid starting locations)
	for (unsigned i = 0; i < MAXITER; ++i) { // Note: This only works for non-new objects so won't work with init objects
		spos    = pos;
		if (spread > 0.0) spos += signed_rand_vector_spherical(spread);
		coll_ix = check_for_obj_coll(spos, sclasses[sclass].calc_cradius()); // what about planet collisions?
		if (coll_ix == 0) break;
	}
	return create_ship(sclass, spos, align, ai, targ, 1);
}


// ************ ship_defs_file_reader ************


class ship_defs_file_reader {

	enum {CMD_GLOBAL_REGEN=0, CMD_RAND_SEED, CMD_SPAWN_DIST, CMD_START_POS, CMD_HYPERSPEED, CMD_SPEED_SCALE, CMD_PLAYER_TURN,
		CMD_SPAWN_HWORLD, CMD_PLAYER_ENEMY, CMD_BUILD_ANY, CMD_TEAM_CREDITS, CMD_SHIP, CMD_WEAP, CMD_WBEAM, CMD_SHIP_WEAP,
		CMD_ADD, CMD_WEAP_PT, CMD_PLAYER_WEAP, CMD_MESH_PARAMS, CMD_SHIP_CYLINDER, CMD_SHIP_CUBE, CMD_SHIP_SPHERE, CMD_SHIP_TORUS,
		CMD_SHIP_BCYLIN, CMD_SHIP_TRIANGLE, CMD_FLEET, CMD_SHIP_ADD_INIT, CMD_SHIP_ADD_GEN, CMD_SHIP_BUILD, CMD_ALIGN, CMD_SHIP_NAMES,
		CMD_ADD_SHIP, CMD_ADD_ASTEROID, CMD_BLACK_HOLE, CMD_PLAYER, CMD_LAST_PARENT, CMD_END};
	ifstream cfg;
	kw_map command_m, ship_m, weap_m, explosion_m, align_m, ai_m, target_m;
	u_ship *last_ship;
	vector<point> weap_pts;
	unsigned cur_ship_type, ship_add_mode;
	bool is_player, player_setup, saw_end, add_ship_enabled, ship_weap_set, last_parent;

	void print_kw_map(kw_map const &strmap) const;
	void strings_to_map(string const &str, kw_map &strmap) const;
	bool read_enum(kw_map &strmap, unsigned &val, string const &type);
	bool read_ship_type(unsigned &type);
	bool read_weap_type(unsigned &type);
	bool parse_command(unsigned cmd);
	bool read_pt(point &pt) {return (cfg >> pt.x >> pt.y >> pt.z) != 0;}
	void setup_keywords();

public:
	ship_defs_file_reader() : last_ship(NULL), cur_ship_type(NUM_US_CLASS), ship_add_mode(CMD_END),
		is_player(0), player_setup(0), saw_end(0), add_ship_enabled(0), ship_weap_set(0), last_parent(0) {}
	static bool read_string(ifstream &in, string &str);
	bool read_file(const char *fn);
};


void ship_defs_file_reader::print_kw_map(kw_map const &strmap) const {

	cout << "stringmap:" << endl;
	for (kw_map::const_iterator i = strmap.begin(); i != strmap.end(); ++i) {
		cout << i->first << ":" << i->second << " ";
	}
	cout << endl;
}


void ship_defs_file_reader::strings_to_map(string const &str, kw_map &strmap) const {

	unsigned const size((unsigned)str.size());

	for (unsigned pos = 0, pos_end = 0, index = 0; pos < size; pos = pos_end) {
		for( ; pos < size && str[pos] == ' '; ++pos) {}
		if (pos == size) break;
		pos_end = min((unsigned)str.find(' ', pos), size);
		assert(pos+1 < pos_end);
		string const substr((str.begin()+pos), (str.begin()+pos_end));
		strmap.insert(make_pair(substr, index++));
	}
	//print_kw_map(strmap);
}


bool ship_defs_file_reader::read_enum(kw_map &strmap, unsigned &val, string const &type) {

	string str;
	
	if (!cfg.good()) {
		cerr << "Error: Bad stream while reading " << type << " enumerated type." << endl;
		return 0;
	}
	if (!(cfg >> str)) {
		cerr << "Error reading " << type << " enumerated type from file." << endl;
		return 0;
	}
	if (strmap.find(str) == strmap.end()) {
		cerr << "Error: Unrecognized " << type << " enumerated type: " << str << "." << endl;
		return 0;
	}
	val = strmap[str];
	return 1;
}


bool ship_defs_file_reader::read_ship_type(unsigned &type) {

	if (!read_enum(ship_m, type, "ship")) return 0;
	assert(type < sclasses.size());
	return 1;
}


bool ship_defs_file_reader::read_weap_type(unsigned &type) {

	if (!read_enum(weap_m, type, "weapon")) return 0;
	assert(type < us_weapons.size());
	return 1;
}


bool ship_defs_file_reader::parse_command(unsigned cmd) {

	unsigned type;
	float dscale(1.0), read_dscale;

	switch (cmd) {
		case CMD_GLOBAL_REGEN: // <float regen_time(s)>
			if (!(cfg >> global_regen)) return 0;
			break;

		case CMD_RAND_SEED: // <unsigned seed>
			{
				unsigned rseed;
				if (!(cfg >> rseed)) return 0;
				srand(rseed);
			}
			break;

		case CMD_SPAWN_DIST: // <float dist>
			if (!(cfg >> spawn_dist)) return 0;
			break;

		case CMD_SPAWN_HWORLD: // <bool respawn_requires_homeworld>
			if (!(cfg >> respawn_req_hw)) return 0;
			break;

		case CMD_PLAYER_ENEMY: // <bool player_is_enemy_team>
			if (!(cfg >> player_enemy)) return 0;
			break;

		case CMD_BUILD_ANY: // <bool teams_can_build_any_ships>
			if (!(cfg >> build_any)) return 0;
			break;

		case CMD_START_POS: // <point start_pos>
			if (!read_pt(ustart_pos)) return 0;
			for (unsigned i = 0; i < 3; ++i) ustart_pos[i] *= CELL_SIZE;
			player_death_pos = universe_origin = last_camera = ustart_pos;
			break;

		case CMD_HYPERSPEED: // <float hyperspeed_mult>
			if (!(cfg >> hyperspeed_mult)) return 0;
			assert(hyperspeed_mult > 0.0);
			break;

		case CMD_SPEED_SCALE: // <float ship_speed_scale>
			if (!(cfg >> ship_speed_scale)) return 0;
			assert(ship_speed_scale > 0.0);
			break;

		case CMD_PLAYER_TURN: // <float player_turn_rate>
			if (!(cfg >> player_turn_rate)) return 0;
			assert(player_turn_rate > 0.0);
			player_turn_rate /= 1000.0;
			break;

		case CMD_TEAM_CREDITS: // <enum alignment> <unsigned num_team_credits(K)>
			{
				unsigned credits, align; // currently, all teams have the same number of credits
				if (!read_enum(align_m, align, "alignment")) return 0;
				if (!(cfg >> credits)) return 0;
				assert(align < NUM_ALIGNMENT);
				regen_uses_credits  = 1;
				credits            *= 1000; // convert from K
				team_credits[align] = init_credits[align] = credits;
			}
			break;

		case CMD_SHIP: // <enum ship_id> <string name> {params...} <enum explosion_t>
			if (!read_ship_type(type)) return 0;
			assert(!sclasses[type].inited);
			if (!sclasses[type].read_from_ifstream(cfg)) return 0; // make this function take the reader (*this)?
			if (!read_enum(explosion_m, sclasses[type].exp_type,    "explosion"    )) return 0;
			if (!read_enum(explosion_m, sclasses[type].exp_subtype, "sub_explosion")) return 0;
			break;

		case CMD_WEAP: // <enum weap_id> <string name> {params...} {<enum ship_id> | <enum weap_id>} <enum explosion_t>
			if (!read_weap_type(type)) return 0;
			assert(!us_weapons[type].inited);
			if (!us_weapons[type].read_from_ifstream(cfg)) return 0;
			if (us_weapons[type].is_fighter) {
				if (!read_ship_type(us_weapons[type].ammo_type)) return 0;
			}
			else {
				if (!read_weap_type(us_weapons[type].ammo_type)) return 0;
			}
			if (!read_enum(explosion_m, us_weapons[type].exp_type, "explosion")) return 0;
			break;

		case CMD_WBEAM: // <enum weap_id> <color brc1> <color brc2> <color beamc1> <color beamc2> <float bw_escale> <bool drains energy> <bool temperature source> <bool paralyzes> <bool mind control> <bool multi-segment>
			if (!read_weap_type(type)) return 0;
			if (!us_weapons[type].read_beam_params_from_ifstream(cfg)) return 0;
			break;

		case CMD_SHIP_WEAP: // <enum ship_id> {[$WEAP_PT] $ADD}*
			if (!read_ship_type(cur_ship_type)) return 0;
			is_player     = 0;
			ship_weap_set = 1;
			weap_pts.clear();
			break;

		case CMD_WEAP_PT: // <point pt>*
			{
				assert(ship_weap_set);
				weap_pts.clear();
				point pt;
				while (read_pt(pt)) {pt.z = -pt.z; weap_pts.push_back(pt);} // negate z
				cfg.clear();
			}
			break;

		case CMD_ADD: // <enum weap_id> [num=1] [ammo=default_ammo]
			{
				if (!read_weap_type(type)) return 0;
				unsigned num, ammo;
				
				if (cfg >> num) {
					if (!(cfg >> ammo)) {
						cfg.clear();
						ammo = 0;
					}
				}
				else {
					cfg.clear();
					num  = 1;
					ammo = 0;
				}
				if (is_player) {
					assert(weap_pts.empty());
					add_player_weap(type, num, ammo);
				}
				else {
					assert(ship_weap_set);
					add_ship_weapon(cur_ship_type, type, num, ammo, weap_pts);
				}
			}
			break;

		case CMD_PLAYER_WEAP: // {$ADD}*
			assert(player_setup);
			is_player = 1;
			weap_pts.clear();
			break;

		case CMD_MESH_PARAMS: // <enum ship_id> <bool deform> <bool remove> <bool expand> <bool uniform_expand> <bool trans>
			{
				if (!read_ship_type(type)) return 0;
				bool params[5];
				for (unsigned i = 0; i < 5; ++i) {
					if (!(cfg >> params[i])) return 0;
				}
				sclasses[type].set_mesh_params(params[0], params[1], params[2], params[3], params[4]);
			}
			break;

		case CMD_SHIP_CYLINDER: // <enum ship_id> <point p1> <point p2> <float r1> <float r2> <end_type> [<float dscale>]
			{
				if (!read_ship_type(type)) return 0;
				int etype;
				point p1, p2;
				float r1, r2;
				if (!read_pt(p1) || !read_pt(p2) || !(cfg >> r1 >> r2 >> etype)) return 0;
				assert(etype >= 0 && etype <= 2); // 0=none, 1=flat, 2=spherical
				if (cfg >> read_dscale) dscale = read_dscale; else cfg.clear();
				sclasses[type].add_bcylinder(ship_cylinder(p1, p2, r1, r2, (etype == 1), dscale));

				if (etype == 2) { // spherical ends - doesn't really work correctly for sphere intersections
					sclasses[type].add_bsphere(p1, r1, dscale);
					sclasses[type].add_bsphere(p2, r2, dscale);
				}
			}
			break;

		case CMD_SHIP_CUBE: // <enum ship_id> <float x1> <float x2> <float y1> <float y2> <float z1> <float z2> [<float dscale>]
			{
				if (!read_ship_type(type)) return 0;
				float x1, x2, y1, y2, z1, z2;
				if (!(cfg >> x1 >> x2 >> y1 >> y2 >> z1 >> z2)) return 0;
				if (cfg >> read_dscale) dscale = read_dscale; else cfg.clear();
				sclasses[type].add_bcube(x1, x2, y1, y2, z1, z2, dscale);
			}
			break;

		case CMD_SHIP_SPHERE: // <enum ship_id> <point center> <float radius> [<float dscale>]
			{
				if (!read_ship_type(type)) return 0;
				point center;
				float r;
				if (!read_pt(center) || !(cfg >> r)) return 0;
				if (cfg >> read_dscale) dscale = read_dscale; else cfg.clear();
				sclasses[type].add_bsphere(center, r, dscale);
			}
			break;

		case CMD_SHIP_TORUS: // <enum ship_id> <point center> <float r_inner> <float r_outer> [<float dscale>]
			{
				if (!read_ship_type(type)) return 0;
				point center;
				float ri, ro;
				if (!read_pt(center) || !(cfg >> ri >> ro)) return 0;
				if (cfg >> read_dscale) dscale = read_dscale; else cfg.clear();
				sclasses[type].add_btorus(center, ri, ro, dscale);
			}
			break;

		case CMD_SHIP_BCYLIN: // <enum ship_id> <point p1> <point p2> <float r1> <float r2>  <float x1> <float x2> <float y1> <float y2> <float z1> <float z2> [<float dscale>]
			{
				if (!read_ship_type(type)) return 0;
				point p1, p2;
				float r1, r2, x1, x2, y1, y2, z1, z2;
				if (!read_pt(p1) || !read_pt(p2) || !(cfg >> r1 >> r2  >> x1 >> x2 >> y1 >> y2 >> z1 >> z2)) return 0;
				if (cfg >> read_dscale) dscale = read_dscale; else cfg.clear();
				sclasses[type].add_bcylin_cube(ship_cylinder(p1, p2, r1, r2, 1, dscale), x1, x2, y1, y2, z1, z2); // check_ends=1
			}
			break;

		case CMD_SHIP_TRIANGLE: // <enum ship_id> <point p1> <point p2> <point p3>
			{
				if (!read_ship_type(type)) return 0;
				triangle tri;
				if (!read_pt(tri.pts[0]) || !read_pt(tri.pts[1]) || !read_pt(tri.pts[2])) return 0;
				sclasses[type].add_triangle(tri);
			}
			break;

		case CMD_FLEET: // <string name> <unsigned multiplier> <enum alignment> <enum ai_type> <enum target_type> <float rand_gen_dist> <point pos> <unsigned counts>+ [<float flagship_child_stray_dist> <enum ship_id>]
			{
				bool using_flagship(0);
				unsigned multiplier;
				unsigned align, ai_type, targ_type;
				float rgen_dist, fc_stray_dist(0.0);
				point pos;
				string name;
				unsigned counts[NUM_US_CLASS];
				if (!read_string(cfg, name)) return 0;
				if (!(cfg >> multiplier))    return 0;
				if (!read_enum(align_m,  align,     "alignment"))   return 0;
				if (!read_enum(ai_m,     ai_type,   "ai_type"))     return 0;
				if (!read_enum(target_m, targ_type, "target_type")) return 0;
				if (!(cfg >> rgen_dist) || !read_pt(pos))           return 0;
				
				for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
					if (!(cfg >> counts[i])) return 0;
					assert(multiplier*counts[i] < 100000);
					
					for (unsigned j = 0; j < multiplier*counts[i]; ++j) {
						build_types[align].push_back(i);
					}
				}
				if (cfg >> fc_stray_dist) { // optional flagship params
					if (!read_ship_type(type)) return 0;
					using_flagship = 1;
				} else cfg.clear();

				if (multiplier > 0) {
					us_fleet fleet(name, align, ai_type, targ_type, rgen_dist, (ustart_pos + pos), counts, multiplier);
					if (using_flagship) fleet.set_flagship(type, fc_stray_dist/1000.0);
					fleet.spawn();
				}
			}
			break;

		case CMD_SHIP_BUILD: // <unsigned counts>+
			for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
				unsigned count;
				if (!(cfg >> count)) return 0;
				
				for (unsigned align = 0; align < NUM_ALIGNMENT; ++align) { // common across all alignments
					for (unsigned j = 0; j < count; ++j) {
						build_types[align].push_back(i);
					}
				}
			}
			break;

		case CMD_SHIP_ADD_INIT: // <bool enabled> {$ALIGN}*
		case CMD_SHIP_ADD_GEN:  // <bool enabled> {$ALIGN}*
			{
				ship_add_mode = cmd;
				if (!(cfg >> add_ship_enabled)) return 0;
				if (cmd == CMD_SHIP_ADD_GEN) allow_add_ship = add_ship_enabled;
			}
			break;

		case CMD_ALIGN: // <enum alignment> <unsigned num> <unsigned counts>+
			{
				assert(ship_add_mode == CMD_SHIP_ADD_INIT || ship_add_mode == CMD_SHIP_ADD_GEN);
				unsigned align, num;
				if (!read_enum(align_m, align, "alignment")) return 0;
				if (!(cfg >> num)) return 0;
				bool const add_init(ship_add_mode == CMD_SHIP_ADD_INIT);
				unsigned *prob(ship_add_prob[add_init][align]);
				
				for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
					if (!(cfg >> prob[i])) return 0;
					assert(prob[i] == 0 || !sclasses[i].orbiting_dock);
				}
				if (add_ship_enabled && add_init && num > 0) add_other_ships(align, num, 1);
				if (!add_init) gen_counts[align] = num;
			}
			break;

		case CMD_SHIP_NAMES: // <unsigned num> <string names>*
			{
				unsigned num;
				if (!(cfg >> num)) return 0;
				string name;
				for (unsigned i = 0; i < num; ++i) {
					if (!read_string(cfg, name)) return 0;
					ship_names.push_back(name);
				}
			}
			break;

		case CMD_ADD_SHIP: // <enum ship_id> <unsigned num> <enum alignment> <enum ai_type> <enum target_type> <bool guardian> <bool kamikaze> [px py pz]
			{
				bool guardian, kamikaze;
				unsigned num, align, ai_type, targ_type;
				float spread(0.0);
				point pos;
				if (!read_ship_type(type)) return 0;
				if (!(cfg >> num))         return 0;
				if (!read_enum(align_m,  align,     "alignment"))   return 0;
				if (!read_enum(ai_m,     ai_type,   "ai_type"))     return 0;
				if (!read_enum(target_m, targ_type, "target_type")) return 0;
				if (!(cfg >> guardian >> kamikaze)) return 0;
				if (!read_pt(pos)) {
					pos    = all_zeros;
					spread = spawn_dist;
					cfg.clear();
				}
				if (guardian) ai_type |= AI_GUARDIAN;
				if (kamikaze) ai_type |= AI_KAMIKAZE;
				assert(num == 0 || !sclasses[type].orbiting_dock);
				u_ship *const parent(last_ship);

				for (unsigned i = 0; i < num; ++i) {
					last_ship = add_ship(type, align, ai_type, targ_type, (pos + ustart_pos), spread);
					assert(last_ship);
					
					if (last_parent) {
						assert(parent);
						last_ship->set_parent(parent);
					}
				}
				last_parent = 0;
			}
			break;

		case CMD_ADD_ASTEROID: // <unsigned num> <unsigned model> <float min_radius> <float max_radius> [px py pz]
			{
				unsigned model, num;
				float r1, r2, dist(0.0);
				point pos;
				if (!(cfg >> num >> model >> r1 >> r2)) return 0;
				assert(r1 <= r2 && model < 4);
				bool const pos_set(read_pt(pos));
				if (!pos_set) cfg.clear();
				assert(num == 1 || !pos_set);

				for (unsigned i = 0; i < num; ++i) { // could check for collisions
					point const pos(ustart_pos + (pos_set ? pos : signed_rand_vector_spherical(spawn_dist)));
					add_uobj(new uobj_asteroid(pos, rand_uniform(r1, r2), model));
				}
			}
			break;

		case CMD_BLACK_HOLE: // <point pos> <float radius>
			{
				point pos;
				float radius;
				if (!read_pt(pos) || !(cfg >> radius)) return 0;
				add_uobj(new stationary_obj(SO_BLACK_HOLE, (ustart_pos + pos), radius));
			}
			break;

		case CMD_PLAYER: // <enum ship_id> <enum alignment>
			{
				assert(!player_setup);
				unsigned align;
				if (!read_ship_type(type)) return 0;
				if (!read_enum(align_m, align, "alignment")) return 0;
				if (!create_player_ship(type, align)) {
					cerr << "Error: Player ship has already been created." << endl;
					return 0;
				}
				player_setup = 1;
			}
			break;

		case CMD_LAST_PARENT:
			last_parent = 1;
			break;

		case CMD_END:
			assert(!saw_end);
			saw_end = 1;
			break;

		default:
			assert(0);
	}
	return 1;
}


void ship_defs_file_reader::setup_keywords() {

	string const commands  ("$GLOBAL_REGEN $RAND_SEED $SPAWN_DIST $START_POS $HYPERSPEED $SPEED_SCALE $PLAYER_TURN $SPAWN_HWORLD $PLAYER_ENEMY $BUILD_ANY $TEAM_CREDITS $SHIP $WEAP $WBEAM $SHIP_WEAP $ADD $WEAP_PT $PLAYER_WEAP $MESH_PARAMS $SHIP_CYLINDER $SHIP_CUBE $SHIP_SPHERE $SHIP_TORUS $SHIP_BCYLIN $SHIP_TRIANGLE $FLEET $SHIP_ADD_INIT $SHIP_ADD_GEN $SHIP_BUILD $ALIGN $SHIP_NAMES $ADD_SHIP $ADD_ASTEROID $BLACK_HOLE $PLAYER $LAST_PARENT $END");
	string const ship_strs ("USC_FIGHTER USC_X1EXTREME USC_FRIGATE USC_DESTROYER USC_LCRUISER USC_HCRUISER USC_BCRUISER USC_ENFORCER USC_CARRIER USC_ARMAGEDDON USC_SHADOW USC_DEFSAT USC_STARBASE USC_BCUBE USC_BSPHERE USC_BTCUBE USC_BSPH_SM USC_BSHUTTLE USC_TRACTOR USC_GUNSHIP USC_NIGHTMARE USC_DWCARRIER USC_DWEXTERM USC_WRAITH USC_ABOMIN USC_REAPER USC_DEATH_ORB USC_SUPPLY USC_ANTI_MISS USC_JUGGERNAUT USC_SAUCER USC_SAUCER_V2 USC_MOTHERSHIP USC_HUNTER USC_SEIGE USC_COLONY USC_ARMED_COL USC_HW_COL USC_STARPORT USC_HW_SPORT");
	string const weap_strs ("UWEAP_NONE UWEAP_TARGET UWEAP_QUERY UWEAP_RENAME UWEAP_DESTROY UWEAP_PBEAM UWEAP_EBEAM UWEAP_REPULSER UWEAP_TRACTORB UWEAP_G_HOOK UWEAP_LRCPA UWEAP_ENERGY UWEAP_ATOMIC UWEAP_ROCKET UWEAP_NUKEDEV UWEAP_TORPEDO UWEAP_EMP UWEAP_PT_DEF UWEAP_DFLARE UWEAP_CHAFF UWEAP_FIGHTER UWEAP_B_BAY UWEAP_CRU_BAY UWEAP_SOD_BAY UWEAP_BOARDING UWEAP_NM_BAY UWEAP_RFIRE UWEAP_FUSCUT UWEAP_SHIELDD UWEAP_THUNDER UWEAP_ESTEAL UWEAP_WRAI_BAY UWEAP_STAR UWEAP_HUNTER UWEAP_DEATHORB UWEAP_LITNING UWEAP_INFERNO UWEAP_PARALYZE UWEAP_MIND_C UWEAP_SAUC_BAY UWEAP_SEIGEC");
	string const exp_strs  ("ETYPE_NONE ETYPE_FIRE ETYPE_NUCLEAR ETYPE_ENERGY ETYPE_ATOMIC ETYPE_PLASMA ETYPE_EMP ETYPE_STARB ETYPE_FUSION ETYPE_EBURST ETYPE_ESTEAL ETYPE_ANIM_FIRE ETYPE_SIEGE");
	string const align_strs("NEUTRAL PLAYER GOV PIRATE RED BLUE ORANGE PURPLE");
	string const ai_strs   ("AI_IGNORE AI_RETREAT AI_ATT_WAIT AI_ATT_ENEMY AI_ATT_ALL AI_SEEKING AI_NONE");
	string const targ_strs ("TARGET_CLOSEST TARGET_ATTACKER TARGET_LAST TARGET_PARENT");
	strings_to_map(commands,   command_m);
	strings_to_map(ship_strs,  ship_m);
	strings_to_map(weap_strs,  weap_m);
	strings_to_map(exp_strs,   explosion_m);
	strings_to_map(align_strs, align_m);
	strings_to_map(ai_strs,    ai_m);
	strings_to_map(targ_strs,  target_m);
}


bool ship_defs_file_reader::read_string(ifstream &in, string &str) { // can be called by outside functions

	str.clear();

	if (!in.good()) {
		cerr << "Error: Bad stream while reading string." << endl;
		return 0;
	}
	int c;
	do {c = in.get();} while (c == ' ' || c == '\t' || c == '\n');

	if (c != '"') {
		cerr << "Error: String does not begin with a double quote: " << char(c) << "." << endl;
		return 0;
	}
	while (1) {
		int const c(in.get());
		if (c == '"') break;
		
		if (c == 0 || c == EOF) {
			cerr << "Error: Unterminated string." << endl;
			return 0;
		}
		str.push_back(char(c));
	}
	return 1;
}


bool ship_defs_file_reader::read_file(const char *fn) {

	// try to open the file
	assert(fn);
	cfg.open(fn);

	if (!cfg.good()) {
		cerr << "Error: Could not open ship definitions file." << endl;
		return 0;
	}
	sclasses.resize(NUM_US_CLASS);
	us_weapons.resize(NUM_UWEAP);
	setup_keywords();

	// parse the file
	string str;

	while (!saw_end && (cfg >> str)) {
		//cout << "str = " << str << endl;
		if (!cfg.good()) {
			cerr << "Error: Bad file stream." << endl;
			return 0;
		}
		assert(!str.empty());
		
		if (str[0] == '#') { // comment
			while (1) {
				int const c(cfg.get());
				if (c == '\n' || c == 0 || c == EOF) break;
			}
		}
		else if (command_m.find(str) == command_m.end()) {
			cerr << "Error: unrecognized command keyword " << str << "." << endl;
			return 0;
		}
		else if (!parse_command(command_m[str])) {
			cerr << "Error reading command " << str << "." << endl;
			return 0;
		}
	}
	if (!cfg.good()) {
		cerr << "Error: Bad stream at end of read." << endl;
		return 0;
	}
	if (!saw_end) cout << "Warning: Missing end of file command." << endl;

	if (!player_setup) {
		cerr << "Error: $PLAYER command is required." << endl;
		return 0;
	}

	// postprocessing and error checking
	for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
		sclasses[i].setup(i);
	}
	for (unsigned i = 0; i < NUM_UWEAP; ++i) {
		us_weapons[i].setup(i);
	}
	return 1;
}


// ************ us_class ************


bool read_color(ifstream &in, colorRGBA &color) {

	return ((in >> color.R >> color.G >> color.B >> color.A) != 0);
}


bool us_class::read_from_ifstream(ifstream &in) {

	if (!ship_defs_file_reader::read_string(in, name)) return 0;
	float ddelay, rdelay, kcost;
	if (!(in >> kcost >> ncrew >> nengines >> radius >> cr_scale >> mass >> cargo >> exp_scale >> accel >> decel >> roll_rate
		>> max_speed >> max_turn >> stability >> max_shields >> max_armor >> shield_re >> armor_re >> max_t
		>> hull_str >> damage_abs >> min_att_dist >> min_app_dist >> sensor_dist >> fire_dist >> stray_dist >> reversible
		>> stoppable >> has_hyper >> has_fast_speed >> mpredict >> has_cloak >> regen_fighters >> regen_ammo >> regen_crew
		>> parallel_fire >> symmetric >> self_shadow >> cont_frag >> for_boarding >> can_board >> orbiting_dock
		>> dynamic_cobjs >> uses_tdir >> emits_light >> engine_lights >> suicides >> kamikaze >> no_disable >> uses_mesh2d
		>> turreted >> weap_spread >> shield_sects >> draw_passes >> ddelay >> rdelay)) return 0;
	if (!read_color(in, base_color)) return 0;
	death_delay = unsigned(TICKS_PER_SECOND*ddelay);
	regen_delay = ((rdelay > 0.0 || global_regen > 0.0) ? (death_delay + unsigned(TICKS_PER_SECOND*(rdelay + global_regen))) : 0);
	mesh_deform = mesh_remove = mesh_expand = mu_expand = mesh_trans = 0;
	accel      *= ship_speed_scale;
	decel      *= ship_speed_scale;
	max_speed  *= ship_speed_scale;
	radius     /= 1000.0;
	accel      /= 1000.0;
	decel      /= 1000.0;
	roll_rate  /= 1000.0;
	max_speed  /= 1000.0;
	max_turn   /= 1000.0;
	stray_dist /= 1000.0;
	shield_re  /= TICKS_PER_SECOND;
	armor_re   /= TICKS_PER_SECOND;
	cost        = unsigned(1000.0*kcost);
	if (for_boarding) can_board = 0; // can't board a boarding shuttle (for now)
	bnd_sphere  = ship_sphere(all_zeros, cr_scale);
	inited      = 1;
	return 1;
}


void us_class::setup(unsigned sclass_) {

	if (!inited) {cout << "Error: Ship ID " << sclass_ << " was not initialized" << endl; exit(1);}
	assert(exp_type < NUM_ETYPES);
	assert(regen_delay  >= 0.0 && death_delay >= 0.0);
	assert(cr_scale     >= 1.0); // might want to relax this later
	assert(stability    >  0.0);
	assert(max_t        > 0.0);
	assert(shield_sects >  0);
	assert(draw_passes  >  0);
	//assert(!reversible || !stoppable);
	assert(mass > 0.0);
	assert(!suicides      || kamikaze);
	assert(!for_boarding  || ncrew > 0);
	assert(!can_board     || ncrew > 0);
	assert(has_fast_speed || !has_hyper);
	sclass = sclass_;

	if (!cobj_triangles.empty()) {
		ship_triangle_list *tlist(new ship_triangle_list(bnd_sphere));
		
		for (vector<triangle>::const_iterator i = cobj_triangles.begin(); i != cobj_triangles.end(); ++i) {
			tlist->add_triangle(*i);
		}
		cobjs.push_back(tlist);
	}
	if (!dynamic_cobjs) { // check to make sure none of the turretted weapon points are within a bounding collision shape
		u_ship test_ship(sclass, all_zeros, ALIGN_NEUTRAL, AI_IGNORE, TARGET_CLOSEST, 0);

		for (unsigned w = 0; w < weapons.size(); ++w) {
			if (!test_ship.weap_turret(weapons[w].wclass)) continue;

			for (unsigned p = 0; p < weapons[w].weap_pts.size(); ++p) {
				assert(!test_ship.sphere_int_obj(weapons[w].weap_pts[p], 0.0)); // error checking
			}
		}
		test_ship.status = 2; // set bad flag so that the destructor can be called
	}
}


void us_class::clear_cobjs() {

	for (unsigned i = 0; i < cobjs.size(); ++i) {
		delete cobjs[i];
	}
	cobjs.clear();
	cobj_triangles.clear();
}


// ************ us_weapon ************


bool us_weapon::read_from_ifstream(ifstream &in) {

	if (!ship_defs_file_reader::read_string(in, name)) return 0;
	float bt;
	if (!(in >> cost >> ammo_cost >> radius >> c_radius >> bradius >> damage >> fire_delay >> firing_error >> regen_time
		>> range >> speed >> seek_dist >> def_ammo >> nshots >> lifetime >> bt >> max_t >> mass >> w_mass >> a_mass >> force
		>> f_inv >> armor >> preference >> hit_proj >> hit_all >> c2_flag >> no_coll >> no_exp_dam >> const_dam >> no_ffire
		>> is_beam >> secondary >> hyper_fire >> point_def >> is_decoy >> ignores_shields >> shield_d_only >> no_light
		>> parallel_fire >> turreted >> auto_orient >> no_ship_vel >> det_on_exp >> symmetric >> is_fighter >> do_regen)) return 0;
	btime       = int(TICKS_PER_SECOND*bt);
	radius     /= 1000.0;
	bradius    /= 1000.0;
	range      /= 1000.0;
	speed      /= 1000.0;
	force      /= 1000.0;
	fire_delay *= TICKS_PER_SECOND;
	lifetime   *= TICKS_PER_SECOND;
	regen_time *= TICKS_PER_SECOND;
	seek_dist  *= radius;
	c_radius   *= radius;
	seeking     = (seek_dist > 0.0);
	inited      = 1;
	return 1;
}


bool us_weapon::read_beam_params_from_ifstream(ifstream &in) {

	assert(inited);
	return bwp.read(in);
}


bool beam_weap_params::read(ifstream &in) {

	for (unsigned i = 0; i < 2; ++i) {
		if (!read_color(in,brc[i]))   return 0;
	}
	for (unsigned i = 0; i < 2; ++i) {
		if (!read_color(in, beamc[i])) return 0;
	}
	return ((in >> bw_escale >> energy_drain >> temp_src >> paralyze >> mind_control >> multi_segment) != 0);
}


void us_weapon::setup(unsigned wclass_) {

	if (!inited) {cout << "Error: Weapon ID " << wclass_ << " was not initialized" << endl; exit(1);}
	assert(damage == 0.0 || !ignores_shields || !shield_d_only); // otherwise there would be no damage
	assert(int(is_beam) + int(is_fighter) + int(is_decoy) <= 1); // can only have one of these flags set
	assert(!is_beam     || !det_on_exp);
	assert(!is_beam     || !no_ship_vel);
	assert(!do_regen    || !is_fighter); // too strict?
	assert(!det_on_exp  || ammo_type != UWEAP_NONE);
	assert(mass > 0.0   || ammo_type == UWEAP_NONE || is_fighter);
	assert(radius > 0.0 || ammo_type == UWEAP_NONE || is_fighter);
	if (range == 0.0) range = speed*lifetime;
	wclass = wclass_;
	
	if (is_fighter) {
		assert(a_mass == 0.0 && ammo_cost == 0 && max_t == 0.0);
		assert(ammo_type < sclasses.size());
		sclasses_to_weapons[ammo_type] = (unsigned)us_weapons.size();
		a_mass    = sclasses[ammo_type].mass;
		ammo_cost = sclasses[ammo_type].cost;
		max_t     = sclasses[ammo_type].max_t;
	}
	if (is_beam) { // can be only one of these
		assert(int(bwp.energy_drain) + int(bwp.temp_src) + int(bwp.paralyze) + int(bwp.mind_control) <= 1);
	}
	calc_preference();
}


void us_weapon::calc_preference() {

	float dscale(1.0);
	if (const_dam)                   dscale *= lifetime;
	if (is_beam && bwp.paralyze)     dscale *= 0.2;
	if (is_beam && bwp.energy_drain) dscale *= 1.5;
	preference += 1.0*is_beam; // beam bonus (fires every frame with instant hit)
	preference += 5.0*dscale*damage/(fire_delay + 1.0/TICKS_PER_SECOND); // weapon damage bonus (damage/delay)
	preference += 250.0*speed; // weapon speed bonus
	preference += 800.0*radius; // radius bonus
	preference += 1.5*(speed == 0.0); // infinite speed bonus
	preference += 400.0*bradius; // blast radius bonus
	preference += 1.5*seeking; // seeking bonus
	preference += 1.8*turreted; // turreted bonus
	preference += 2.0*(1 - secondary); // primary weapon bonus (multiple weapons can fire)
	preference += 1.4*do_regen; // regen bonus
	preference += 1000.0*is_fighter; // fighter bonus
	preference += 0.05*armor; // weapon armor bonus
	preference += 0.5*const_dam; // constant damage bonus (?)
	preference += 1.5*no_ffire; // no friendly fire bonus
	preference += 10.0*(exp_type == ETYPE_EBURST); // energy burst bonus
	preference += 1.7*is_decoy; // decoy bonus
	preference -= 2.0*firing_error; // firing error penalty
	preference -= 3.0*point_def; // point defense penalty (save for point/missile defense use)
	if (is_beam) preference += 5.0*bwp.temp_src + 4.0*bwp.paralyze + 5.0*bwp.mind_control;
}


// ************ other stuff ************


void print_ship_ratings() {

	assert(sclasses.size() == NUM_US_CLASS);
	assert(us_weapons.size() == NUM_UWEAP);
	int const cwidth(13);
	float norm_value(1.0), last_value(1.0);

	for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
		us_class const &sc(sclasses[i]);
		float const offense(int(TICKS_PER_SECOND*sc.offense_rating() + 0.5));
		float const defense(sc.defense_rating());
		float const value(int(0.001*offense*defense + 0.5));
		if (i == 0) norm_value = value;
		float const nv(0.1*(int(10.0*value/norm_value + 0.5)));
		float const lv(0.1*(int(10.0*value/last_value + 0.5)));
		cout << sc.name;
		print_n_spaces(cwidth - (int)sc.name.size());
		cout << " O: " << offense << " D: " << defense << " V: " << value << " NV: " << nv << " LV: " << lv
			 << " WG: " << sc.mass << " U: " << sc.used_mass() << endl;
		cout << "cost: " << sc.cost << ", weap_cost: " << sc.weap_cost() << ", ammo_cost: " << sc.ammo_cost() << endl;
		last_value = value;
	}
}


void init_ship_weapon_classes() {

	ship_defs_file_reader reader;

	if (!reader.read_file(ship_def_file)) {
		cerr << "Error reading ship definitions file '" << ship_def_file << "." << endl;
		exit(1);
	}
	assert(sclasses.size()   == NUM_US_CLASS);
	assert(us_weapons.size() == NUM_UWEAP);
	player_init_weapons = player_ship().weapons;
	if (SHOW_SHIP_RATINGS) print_ship_ratings(); // print out offense/defense ratings
}


void add_other_ships(int align, unsigned num, bool initial) {

	assert(align < NUM_ALIGNMENT);
	
	if (!initial) {
		if (!allow_add_ship) return;
		num = gen_counts[align];
	}
	if (num == 0) return;
	unsigned psum[NUM_US_CLASS];
	upos_point_type const &player(initial ? ustart_pos : get_player_pos2()); // can't call get_player_pos() while initializing
	unsigned const NS_NAMES(8);

	for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
		psum[i] = ship_add_prob[initial][align][i];
		if (i > 0) psum[i] += psum[i-1];
	}
	unsigned const ptot(psum[NUM_US_CLASS-1]);
	assert(ptot > 0);

	for (unsigned i = 0; i < num; ++i) { // *** generate based on random but deterministic function ***
		unsigned const randval(rand()%ptot);
		unsigned sclass, ai_type(AI_ATT_ENEMY);

		for (sclass = 0; sclass < NUM_US_CLASS; ++sclass) {
			if (randval < psum[sclass]) break;
		}
		assert(sclass < NUM_US_CLASS);
		u_ship *ship(add_ship(sclass, align, ai_type, TARGET_CLOSEST, player, spawn_dist)); // spawn near the player
		assert(ship != NULL);
		if (!ship_names.empty()) ship->rename(ship_names[rand()%ship_names.size()]);
	}
}


void merge_weapons(vector<ship_weapon> &weapons, ship_weapon const &w) {
	
	bool found(0);

	for (unsigned i = 0; i < weapons.size(); ++i) { // check to see if ship already has that weapon type
		if (weapons[i].wclass == w.wclass) { // same weapon, merge contents
			weapons[i].init_ammo += w.init_ammo;
			weapons[i].ammo      += w.ammo;
			weapons[i].wcount    += w.wcount;
			weapons[i].nregen    += w.nregen;
			
			if (w.docked != NULL) {
				weapons[i].check_docked();

				for (unsigned j = 0; j < w.docked->size(); ++j) {
					weapons[i].docked->push_back((*w.docked)[j]);
				}
			}
			found = 1;
			break;
		}
	}
	if (!found) weapons.push_back(w);
}


void clear_usc_cobjs() { // can call at the end to free up memory (but so far unused)

	for (unsigned i = 0; i< sclasses.size(); ++i) {
		sclasses[i].clear_cobjs();
	}
}



