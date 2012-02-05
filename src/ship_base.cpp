// 3D World - definitions for u_ship_base, us_class, us_weapon, ship_weapon
// by Frank Gennari
// 12/9/05

#include "ship.h"
#include "ship_util.h"


float const DAMAGED_WEAP_REPAIR = 8.0*TICKS_PER_SECOND;
float const DAMAGED_SHIELD_RE   = 0.1;
float const DOCKED_REGEN_RATE   = 10.0;
float const ORBIT_REGEN_RATE    = 25.0;
float const DOCK_REGEN_BONUS    = 0.25; // 0.0 to 1.0
float const ARMOR_ENERGY_EFF    = 0.5;
float const SHIELD_ENERGY_EFF   = 1.0;
float const GROW_ENERGY_EFF     = 4.0;
float const MAX_SHIP_SIZE_SCALE = 10.0;

bool const FIGHTER_WEAP_RANGE   = 1;


extern int iticks, frame_counter;
extern float fticks;
extern vector<us_class> sclasses;
extern vector<us_weapon> us_weapons;
extern map<unsigned, unsigned> sclasses_to_weapons;


// ************ U_SHIP_BASE ************


u_ship_base::u_ship_base(unsigned sclass_) : sclass(sclass_), ncredits(0), kills(0), tot_kills(0), deaths(0), regened(1),
	energy(0.0), fuel(1.0), used_cargo(0.0), size_scale(1.0)
{
	assert(sclass < sclasses.size());
	assert(specs().inited);
}


void u_ship_base::create_from(u_ship const *const ship) { // related to u_ship::create_from()

	ncrew      = ship->get_ncrew();
	shields    = ship->get_shields() + DOCK_REGEN_BONUS*(ship->get_max_shields() - ship->get_shields());
	armor      = ship->get_armor()   + DOCK_REGEN_BONUS*(ship->get_max_armor()   - ship->get_armor());
	energy     = ship->get_energy();
	size_scale = ship->size_scale;
	ncredits   = ship->ncredits;
	kills      = ship->kills;
	tot_kills  = ship->tot_kills;
	deaths     = ship->deaths;
	fuel       = ship->fuel;
	used_cargo = ship->used_cargo;
	weapons    = ship->weapons;
}


// Note: Can't be in u_ship_base since it can be copied
void u_ship_base::free_weapons() { // recursive

	for (unsigned i = 0; i < weapons.size(); ++i) {
		if (weapons[i].docked != NULL) {
			for (unsigned j = 0; j < weapons[i].docked->size(); ++j) { // free docked ships of the docked ships
				(*weapons[i].docked)[j].free_weapons();
			}
			delete weapons[i].docked;
		}
	}
	weapons.clear();
}


void u_ship_base::check_ref_objs(u_ship *cur_ship) {

	assert(cur_ship);

	for (set<u_ship *>::iterator i = fighters.begin(); i != fighters.end(); ) {
		set<u_ship *>::iterator cur(i++);
		u_ship *fighter(*cur);
		assert(fighter != NULL);

		if (fighter->invalid()) { // a fighter has been destroyed
			if (fighter->get_parent() == cur_ship || fighter->get_align() != cur_ship->get_align()) { // still our fighter (or captured)
				if (!fighter->docked) build_fighter(fighter->sclass);
				accept_fighters_from(fighter, cur_ship); // add the living fighters in the dead fighter tree below us to our own fighters
			}
			fighters.erase(cur);
		}
		else if (fighter->get_parent() != cur_ship) { // in case the fighter was captured
			fighters.erase(cur);
		}
	}
	if (VERIFY_REFS) {
		for (set<u_ship *>::const_iterator i = fighters.begin(); i != fighters.end(); ++i) {
			assert(*i && !(*i)->invalid() && !(*i)->to_be_removed());
			(*i)->verify_status();
		}
	}
}


bool u_ship_base::find_weapon_atype(unsigned atype, bool fighter, unsigned &ix) const {

	for (unsigned w = 0; w < weapons.size(); ++w) {
		if (fighter && !weapons[w].get_usw().is_fighter) continue;
		if (weapons[w].get_usw().ammo_type == atype) {ix = w; return 1;}
	}
	return 0;
}


bool u_ship_base::find_weapon_wclass(unsigned wclass, unsigned &ix) const {

	for (unsigned w = 0; w < weapons.size(); ++w) {
		if (weapons[w].wclass == wclass) {ix = w; return 1;}
	}
	return 0;
}


bool u_ship_base::build_fighter(unsigned fsclass) {

	assert(fsclass < sclasses.size());
	if (!specs().regen_fighters) return 0; // can't build a new one
	unsigned const req_crew(sclasses[fsclass].ncrew);
	unsigned j(0), w(0);
	if (req_crew > (ncrew - specs().req_crew())) return 0; // use crew if available, need half to run ship
	if (!find_weapon_atype(fsclass, 1, j))       return 0; // can't have multiple fighter weapons with the same type
	if (weapons[j].ammo != 0)                    return 0; // only when ammo == 0?
	unsigned tot_owned(weapons[j].ammo + weapons[j].nregen);

	for (set<u_ship *>::const_iterator i = fighters.begin(); i != fighters.end(); ++i) {
		if (!(*i)->invalid() && (*i)->get_sclass() == fsclass) ++tot_owned;
	}
	if (tot_owned >= weapons[w].init_ammo) return 0;
	vector<ship_weapon> const &sc_weap(sclasses[fsclass].weapons);

	for (unsigned k = 0; k < sc_weap.size(); ++k) {
		// if ammo is needed for this weapon, and it isn't automatically regenerated, and this weapon exists on the dock
		if (sc_weap[k].get_usw().need_ammo() && !sc_weap[k].get_usw().do_regen, find_weapon_wclass(sc_weap[k].wclass, w)) { // else no weapon, just create the ammo
			if (weapons[w].ammo < sc_weap[k].init_ammo) return 0; // no ammo for this weapon
		}
	}
	for (unsigned k = 0; k < sc_weap.size(); ++k) {
		if (sc_weap[k].get_usw().need_ammo() && !sc_weap[k].get_usw().do_regen && find_weapon_wclass(sc_weap[k].wclass, w)) {
			unsigned const req_ammo(sc_weap[k].init_ammo);
			assert(weapons[w].ammo >= req_ammo);
			weapons[w].ammo -= req_ammo;
		}
	}
	ncrew -= req_crew;
	++weapons[j].nregen;
	return 1;
}


void u_ship_base::orbital_ship_regen(u_ship *ship) {

	bool was_docked(ship->docked);
	ship->docked = 1;
	ship->regen(ORBIT_REGEN_RATE, this);
	ship->docked = was_docked;
}


void u_ship_base::use_fuel(float val) {

	if (val == 0.0) return;
	assert(val > 0.0);
	
	if (energy > 0.0) {
		float const escale(0.01*get_mass());
		float used_energy(min(energy, fticks*escale*val));
		energy -= used_energy;
		val    -= used_energy/escale;
	}
	if (val > 0.0) {
		float const fscale(0.0001);
		float used_fuel(min(fuel, fticks*fscale*val));
		fuel -= used_fuel;
		val  -= used_fuel/fscale;
	}
	if (fuel <= 0.0) { // out of fuel
		if (size_scale > 2.0) { // die
			armor   = 0.0;
			shields = 0.0;
		}
		else {
			//cout << "out of fuel" << endl;
			// *** WRITE ***
		}
	}
}


void u_ship_base::regen(float rate, u_ship_base *dock) {

	if (docked) assert(dock != NULL); else assert(dock == NULL);
	
	if (docked) {
		regened = 1;

		if (fuel < 1.0 && dock->fuel > 0.0) {
			float const fuel_ratio(dock->get_mass()/get_mass()); // > 1.0
			float const fuel_transfer(min((1.0f - fuel), fuel_ratio*dock->fuel));
			fuel       += fuel_transfer; // refill on fuel
			dock->fuel -= fuel_transfer/fuel_ratio;
		}
	}
	if (specs().regen_crew && ncrew < specs().ncrew && int(rand()%TICKS_PER_SECOND) < iticks) {
		ncrew = min(specs().ncrew, (ncrew + 1)); // regen one crew per second (randomly)
	}
	float max_armor(get_max_armor()), max_shields(get_max_shields()), max_energy(max_armor + max_shields);
	float const shield_rate((armor == max_armor) ? 1.0 : DAMAGED_SHIELD_RE);
	float const cscale(get_crew_scale()); // armor is repaired faster with more crew
	armor   = min(max_armor,   (armor   + fticks*rate*specs().armor_re*cscale));
	shields = min(max_shields, (shields + fticks*rate*specs().shield_re*shield_rate));

	if (max_energy > 0.0 && energy > 1.0) { // use stored energy to regenerate armor and shields
		energy = min(energy, max_energy);

		if (armor < max_armor) {
			float const armor_val(min(ARMOR_ENERGY_EFF*energy, (max_armor - armor)));
			armor   += armor_val;
			energy  -= armor_val/ARMOR_ENERGY_EFF;
		}
		if (shields < max_shields) {
			float const shield_val(min(SHIELD_ENERGY_EFF*energy, (max_shields - shields)));
			shields += shield_val;
			energy  -= shield_val/SHIELD_ENERGY_EFF;
		}
		if (energy > 0.0) {
			if (size_scale >= MAX_SHIP_SIZE_SCALE) { // too large - die
				armor   = 0.0;
				shields = 0.0;
				energy  = 0.0;
			}
			else { // grow
				float const gef(GROW_ENERGY_EFF/max_energy), grow_amt(min(gef*energy, (MAX_SHIP_SIZE_SCALE - size_scale)));
				unsigned const size_val((unsigned)size_scale), new_size_val(unsigned(size_scale + grow_amt));

				if (size_val > 0 && new_size_val > size_val) { // weapons upgrade
					for (unsigned i = 0; i < weapons.size(); ++i) {
						ship_weapon &w(weapons[i]);

						if (us_weapons[w.wclass].is_fighter) {
							if (new_size_val > 2) { // ammo (even fighter) upgrade
								++w.init_ammo;
								++w.ammo;
							}
						}
						else {
							w.init_ammo = (w.init_ammo*new_size_val)/size_val;
							w.ammo      = (w.ammo     *new_size_val)/size_val; // w.nregen?
							w.wcount    = (w.wcount   *new_size_val)/size_val;
						}
					}
				}
				ncrew       = unsigned((ncrew*new_size_val)/size_val);
				size_scale += grow_amt;
				energy     -= grow_amt/gef;
			}
			// Do other stuff...
		}
	}
	for (unsigned i = 0; i < weapons.size(); ++i) { // regenerate ammo for weapons
		weapons[i].regen_ammo(rate, this);

		if (docked && weapons[i].get_usw().need_ammo()) { // see if the dock can give us more ammo
			unsigned ammo_needed(weapons[i].init_ammo - weapons[i].ammo);
			unsigned j(0);

			if (ammo_needed > 0 && dock->find_weapon_wclass(weapons[i].wclass, j)) { // find compatible weapon (Note: could search for compatible ammo, but there could be multiple weapons)
				if (dock->weapons[j].ammo > 0) { // dock has ammo
					unsigned const transfer_ammo(min(ammo_needed, dock->weapons[j].ammo));
					assert(transfer_ammo > 0);
					dock->weapons[j].ammo -= transfer_ammo; // take from the dock
					weapons[i].ammo       += transfer_ammo; // give to the fighter
					ammo_needed           -= transfer_ammo; // now less ammo is needed
					//cout << "transfer of " << transfer_ammo << " ammo, left " << dock->weapons[j].ammo << endl;
				}
			}
		}
	}
}


void u_ship_base::accept_fighters_from(u_ship *ship, u_ship *cur_ship) { // recursive

	assert(ship != NULL && cur_ship != NULL && ship != cur_ship);
	if (cur_ship->invalid() || ship->get_parent() != cur_ship) return;
	
	for (set<u_ship *>::const_iterator i = ship->fighters.begin(); i != ship->fighters.end(); ++i) {
		u_ship *fighter(*i);
		assert(fighter != NULL);

		if (fighter->invalid()) {
			accept_fighters_from(fighter, cur_ship);
		}
		else {
			add_fighter(fighter, cur_ship, 0);
		}
	}
	ship->fighters.clear();
}


void u_ship_base::copy_weapons_from_sclass() {

	copy(specs().weapons.begin(), specs().weapons.end(), back_inserter(weapons));
}


void u_ship_base::add_fighter(u_ship *ship, u_ship *cur_ship, bool from_fighter_bay) {
	
	assert(ship != NULL && cur_ship != NULL);
	if (!cur_ship->invalid()) ship->set_parent(cur_ship); // can be invalid if fighters are launched on explosion
	if (from_fighter_bay)     ship->add_flag(OBJ_FLAGS_FITR);
	assert(fighters.find(ship) == fighters.end());
	fighters.insert(ship);
}


bool u_ship_base::has_ammo_for(unsigned wclass) const {

	assert(wclass < us_weapons.size());
	assert(us_weapons[wclass].need_ammo());
	unsigned i;
	return (find_weapon_wclass(wclass, i) && weapons[i].ammo > 0);
}


bool u_ship_base::has_space_for_fighter(unsigned sclass) const {

	unsigned ix;
	if (!find_weapon_atype(sclass, 1, ix)) return 0;
	assert(ix < weapons.size());
	return weapons[ix].space_for_fighter();
}


void u_ship_base::reset_ammo() { // used for player "cheat"

	for (unsigned i = 0; i < weapons.size(); ++i) {
		weapons[i].ammo = weapons[i].init_ammo;
	}
}


bool u_ship_base::can_attack() const {
	
	return (weapons.empty() || (weapons.size() == 1 && weapons[0].wclass == UWEAP_NONE));
}


bool u_ship_base::can_lead_shot_with(unsigned weap) const {

	assert(weap < weapons.size());
	if (!weapons[weap].can_lead_shot()) return 0;
	return (specs().mpredict || weap_turret(weapons[weap].wclass));
}


bool u_ship_base::need_ammo_for(unsigned wix) const {

	assert(wix < weapons.size());
	unsigned const weapon_id(weapons[wix].wclass);
	assert(weapon_id < us_weapons.size());
	return (us_weapons[weapon_id].need_ammo());
}


void u_ship_base::print_ammo() const {

	for (unsigned i = 0; i < weapons.size(); ++i) {
		if (us_weapons[weapons[i].wclass].need_ammo()) {
			cout << us_weapons[weapons[i].wclass].name << ":" << weapons[i].ammo << " ";
		}
	}
}


bool u_ship_base::check_fire_delay(unsigned wix) const {

	assert(wix < weapons.size());
	ship_weapon const &w(weapons[wix]);
	unsigned wcount(w.wcount);
	if (wcount == 0) return 0; // unlikely to get here
	us_weapon const &usw(us_weapons[w.wclass]);
	if (usw.is_beam) return 1;

	if (usw.parallel_fire || specs().parallel_fire) {
		if (w.weap_pts.size() > 1) {
			wcount = max(1U, wcount/(unsigned)w.weap_pts.size());
		}
		else if (specs().weap_spread > 1) {
			wcount = max(1U, wcount/specs().weap_spread);
		}
	}
	unsigned const fire_delay(unsigned(fticks*(frame_counter - w.last_fframe))); // not exactly right to use fticks
	unsigned const fdelay(unsigned(usw.fire_delay/(float)wcount)); // *** what about high count/fire rate? ***
	return (fire_delay >= fdelay);
}


bool u_ship_base::weap_turret(unsigned weapon_id) const {

	assert(weapon_id < us_weapons.size());
	return ((specs().turreted + us_weapons[weapon_id].turreted) >= 2);
}


bool u_ship_base::bad_angle(float const angle, float target_dist, unsigned weapon_id) const {

	assert(weapon_id < us_weapons.size());
	if (angle > PI_TWO) return 1; // > 90 degrees is always bad
	us_weapon const &weap(us_weapons[weapon_id]);
	if (weap.is_fighter) return 0;
	if (weap.seeking)    return (angle > 0.24 && target_dist > weap.seek_dist); // in seeking range?
	if (weap.is_decoy)   return (angle > 0.18);
	if (weap.const_dam)  return (angle > 0.14);
	return (angle > 0.08);
}


bool u_ship_base::out_of_ammo_for(unsigned wix, bool current_only) const {
	
	assert(wix < weapons.size());
	unsigned const num_weapons(weapons.size());

	for (unsigned i = 0; i < num_weapons; ++i) { // iterate over all weapons
		unsigned const w((wix + i) % num_weapons);
		if (weapons[w].wclass == UWEAP_NONE || weapons[w].wcount == 0) continue;
		//if (weapons[w].get_usw().cost == 0) continue; // not a real weapon
		if (!weapons[w].no_ammo()) return 0; // have ammo for this weapon
		if (current_only)          return 1;
	}
	return 1;
}


float u_ship_base::get_weap_ammo_mass(vector<ship_weapon> const &ws) const {

	float wmass(0.0);

	for (unsigned i = 0; i < ws.size(); ++i) {
		assert(ws[i].wclass < us_weapons.size());
		wmass += ws[i].wcount*us_weapons[ws[i].wclass].w_mass + ws[i].ammo*us_weapons[ws[i].wclass].a_mass;
	}
	return wmass;
}


float u_ship_base::get_damage_after_time(float time_seconds) const {
	
	float const max_shields(get_max_shields()), max_armor(get_max_armor());
	float fin_shields(shields), fin_armor(armor);
	assert(time_seconds >= 0.0);

	if (time_seconds > 0) {
		float const time_ticks(time_seconds*TICKS_PER_SECOND);
		fin_shields = min(max_shields, (shields + time_ticks*specs().shield_re));
		fin_armor   = min(max_armor,   (armor   + time_ticks*specs().armor_re));
	}
	return (1.0 - (fin_shields + fin_armor + 1.0)/(max_shields + max_armor + 1.0));
}


float u_ship_base::get_true_rel_mass_scale() const {

	assert(sclass < sclasses.size());
	assert(specs().mass > 0.0);
	float const sclass_mass(get_weap_ammo_mass(sclasses[sclass].weapons));
	float const cur_weap_mass(get_weap_ammo_mass(weapons));
	return (specs().mass + cur_weap_mass - sclass_mass)/specs().mass;
}


// ************ SHIP_WEAPON ************


ship_weapon::ship_weapon(unsigned weapon, unsigned num, unsigned ammo0, vector<point> const &weap_pts_)
	: wclass(weapon), init_ammo(ammo0), wcount(num), rtime(0), nregen(0), ndamaged(0), last_fframe(0), cur_wpt(0), docked(NULL)
{
	assert(wclass < us_weapons.size());
	if (init_ammo == 0) init_ammo = wcount*us_weapons[wclass].def_ammo; // get the default for that weapon
	ammo     = init_ammo;
	weap_pts = weap_pts_;
}


us_weapon const &ship_weapon::get_usw() const {

	assert(wclass < us_weapons.size());
	return us_weapons[wclass];
}


void ship_weapon::regen_ammo(float rate, u_ship_base *dock) {

	assert(dock != NULL);
	assert(rate > 0.0);

	if (docked != NULL) {
		for (unsigned i = 0; i < docked->size(); ++i) { // for docked fighters
			(*docked)[i].regen(max(rate, DOCKED_REGEN_RATE), dock);
		}
	}
	// can get fighters from a destroyed fighter while regenerating a fighter so that (ammo + nregen) > init_ammo
	assert(ammo <= init_ammo);
	//assert((ammo + nregen) <= init_ammo);

	if (ndamaged > 0) { // repair damaged weapons
		unsigned const dr(unsigned(DAMAGED_WEAP_REPAIR/rate));

		if (dr <= 1 || (rand() % dr) == 0) {
			++wcount;
			--ndamaged;
		}
	}
	if (ammo >= init_ammo || wcount == 0) return; // already at max ammo or have no weapons
	assert(wclass < us_weapons.size());
	us_weapon const &weap(us_weapons[wclass]);
	if (!weap.do_regen && (weap.is_fighter || !dock->specs().regen_ammo) && nregen == 0) return; // this weapon/ship combo doesn't regenerate
	float const regen_t(rate*weap.regen_time/(float)wcount);
	rtime += iticks;
	if (rtime < (unsigned)regen_t) return; // too early to regenerate
	rtime -= (unsigned)regen_t; //rtime = 0;
	unsigned const numregen(max(1U, unsigned(fticks/regen_t)));
	ammo = min(init_ammo, (ammo + numregen));
	if (!weap.do_regen) nregen = (unsigned)max(0, ((int)nregen - (int)numregen));
}


// Note: dock argument is unused
void ship_weapon::dock_ship(u_ship *ship, u_ship *dock) { // can't copy fighters though (unless docked)

	assert(ship != NULL);
	assert(us_weapons[wclass].is_fighter && us_weapons[wclass].ammo_type == ship->get_sclass());
	check_docked();
	docked->push_back(u_ship_base(ship->get_sclass())); // engines are repaired fully when docked
	docked->back().create_from(ship);
	docked->back().regened = 0;
	docked->back().docked  = 1;
	ship->docked = 1;
	ship->status = 1; // remove the fighter, must do this after increasing ammo so that the 'lost' fighter is not regenerated
	ship->weapons.clear(); // don't want to delete docked yet

	if (ship->get_parent() == dock && !ship->is_player_ship()) { // parent gets the fighter's credits
		dock->ncredits += ship->ncredits;
		ship->ncredits  = 0;
	}
}


void ship_weapon::release_ship(u_ship *ship, u_ship *dock) {

	assert(ship != NULL);
	if (docked == NULL || docked->empty()) return; // no docked ships - OK
	assert(us_weapons[wclass].is_fighter && us_weapons[wclass].ammo_type == ship->get_sclass());
	float min_damage(1.0);
	unsigned min_damaged(0);

	for (unsigned i = 0; i < docked->size(); ++i) {
		if (!(*docked)[i].regened) continue; // must regenerate each frame
		float const damage((*docked)[i].get_damage());
		
		if (damage < min_damage) {
			min_damage  = damage;
			min_damaged = i;
		}
	}
	ship->create_from((*docked)[min_damaged]);
	docked->erase(docked->begin() + min_damaged);

	if (dock != NULL && ship->get_ncrew() < ship->specs().ncrew && dock->get_ncrew() > dock->specs().req_crew()) {
		int const needed(ship->specs().ncrew - ship->get_ncrew());
		int const avail(dock->get_ncrew() - dock->specs().req_crew());
		int const transfer(min(needed, avail)); // transfer crew from dock to ship to replace the crew that died (> 0)
		ship->adj_crew( transfer);
		dock->adj_crew(-transfer);
	}
}


bool ship_weapon::space_for_fighter() const {

	return (ammo < max(init_ammo, wcount*us_weapons[wclass].def_ammo));
}


float ship_weapon::min_damage() const {

	if (docked == NULL || docked->empty()) return 0.0; // no docked ships - OK
	assert(us_weapons[wclass].is_fighter);
	float min_damage(1.0);

	for (unsigned i = 0; i < docked->size(); ++i) {
		if (!(*docked)[i].regened) continue;
		min_damage = min(min_damage, (*docked)[i].get_damage());
	}
	return min_damage;
}


// ************ US_CLASS ************


void us_class::set_mesh_params(bool deform, bool remove, bool expand, bool mu_exp, bool trans) {
	
	assert(inited && uses_mesh2d);
	mesh_deform = deform;
	mesh_remove = remove;
	mesh_expand = expand;
	mu_expand   = mu_exp;
	mesh_trans  = trans;
}


void us_class::add_bcube(float x1, float x2, float y1, float y2, float z1, float z2, float dscale) {

	assert(inited);
	assert(x1 < x2 && y1 < y2 && z1 < z2);
	cobjs.push_back(new ship_cube(x1, x2, y1, y2, z1, z2, dscale));
}


void us_class::add_bcylinder(ship_cylinder const &c) {
	
	assert(inited);
	assert(c.r1 >= 0.0 && c.r2 >= 0.0 && (c.r1 > 0.0 || c.r2 > 0.0));
	cobjs.push_back(new ship_cylinder(c));
}


void us_class::add_bsphere(point const &center, float r, float dscale) {

	assert(inited);
	assert(r > 0.0);
	cobjs.push_back(new ship_sphere(center, r, dscale));
}


void us_class::add_btorus(point const &center, float ri, float ro, float dscale) {

	assert(inited);
	assert(ri > 0.0 && ro > 0.0 && ri <= ro);
	cobjs.push_back(new ship_torus(center, ri, ro, dscale));
}


void us_class::add_bcylin_cube(ship_cylinder const &c, float x1, float x2, float y1, float y2, float z1, float z2) {

	assert(inited);
	assert(c.r1 >= 0.0 && c.r2 >= 0.0 && (c.r1 > 0.0 || c.r2 > 0.0));
	assert(x1 < x2 && y1 < y2 && z1 < z2);
	cobjs.push_back(new ship_bounded_cylinder(c, ship_cube(x1, x2, y1, y2, z1, z2, c.get_dscale())));
}


float us_class::offense_rating() const {

	if (offense >= 0.0) return offense; // cached
	offense = 0.0;
	unsigned nsec(0);

	for (unsigned i = 0; i < weapons.size(); ++i) {
		assert(weapons[i].wclass < us_weapons.size());
		if (us_weapons[weapons[i].wclass].secondary) ++nsec;
	}
	for (unsigned i = 0; i < weapons.size(); ++i) {
		us_weapon const &weap(us_weapons[weapons[i].wclass]);

		if (weap.is_fighter) {
			if (weapons[i].init_ammo > 0) {
				assert(weap.ammo_type != sclass); // recursion
				offense += 0.5*weapons[i].init_ammo*sclasses[weap.ammo_type].offense_rating();
			}
		}
		else {
			offense += ((nsec > 1 && weap.secondary) ? 0.75 : 1.0)*weapons[i].wcount*weap.offense_rating();
		}
	}
	if (kamikaze && exp_type != 0) offense += (suicides ? 1.0 : 0.25)*0.0005*EXPLOSION_DAMAGE*radius*exp_scale;
	return offense; // round to nearest int
}


float us_class::defense_rating() const {

	if (defense >= 0.0) return defense; // cached
	defense = (max_shields + max_armor);
	
	for (unsigned i = 0; i < weapons.size(); ++i) {
		us_weapon const &weap(us_weapons[weapons[i].wclass]);

		if (weap.is_fighter && weapons[i].init_ammo > 0) {
			float const mult((weap.do_regen || regen_fighters) ? 1.5 : 1.0);
			defense += 0.25*mult*weapons[i].init_ammo*sclasses[weap.ammo_type].defense_rating();
		}
		if (weap.point_def || weap.is_decoy) defense += 10.0*weapons[i].wcount*weap.damage;
	}
	if (has_cloak) defense *= 2.0;
	return defense;
}


unsigned us_class::weap_cost() const {

	unsigned c(0);

	for (unsigned i = 0; i < weapons.size(); ++i) {
		c += us_weapons[weapons[i].wclass].cost*weapons[i].wcount;
	}
	return c;
}


unsigned us_class::ammo_cost() const {

	unsigned c(0);

	for (unsigned i = 0; i < weapons.size(); ++i) {
		c += us_weapons[weapons[i].wclass].ammo_cost*weapons[i].init_ammo;
	}
	return c;
}


float us_class::used_mass() const {

	float used_w(0.0);

	for (unsigned i = 0; i < weapons.size(); ++i) {
		us_weapon const &weap(us_weapons[weapons[i].wclass]);
		float const tscale((!weap.is_fighter && (turreted + weap.turreted) >= 2) ? 2.0 : 1.0);
		used_w += weapons[i].wcount*weap.w_mass*tscale + weapons[i].init_ammo*weap.a_mass;
	}
	return used_w;
}


bool us_class::can_attack() const {
	
	return (weapons.empty() || (weapons.size() == 1 && weapons[0].wclass == UWEAP_NONE));
}


unsigned us_class::req_crew() const {

	return (unsigned)ceil(SHIP_REQ_CREW*ncrew);
}


float us_class::get_weap_range() const { // caches weap_range

	if (weap_range >= 0.0) return weap_range;
	weap_range = 0.0;
	
	for (unsigned i = 0; i < weapons.size(); ++i) {
		us_weapon const &uw(weapons[i].get_usw());
		float range(uw.range);
		
		if (FIGHTER_WEAP_RANGE && uw.is_fighter) {
			us_class const &sc(sclasses[uw.ammo_type]);
			range = max(sc.sensor_dist, sc.get_weap_range());
			if (sc.stray_dist > 0.0) range = min(range, (radius*cr_scale + 4.0f*sc.stray_dist + sc.get_weap_range()));
		}
		if (range == 0.0) {weap_range = 0.0; break;} // unranged weapon (fighter, etc.)
		weap_range = max(range, weap_range);
	}
	return weap_range;
}


float us_class::get_min_weap_range() const { // does not cache weap_range

	float min_weap_range(0.0);
	
	for (unsigned i = 0; i < weapons.size(); ++i) {
		us_weapon const &uw(weapons[i].get_usw());
		if (uw.range == 0.0) return 0.0; // infinite range (includes fighters)
		min_weap_range = ((min_weap_range == 0.0) ? uw.range : min(uw.range, min_weap_range));
	}
	return min_weap_range;
}


// ************ US_WEAPON ************


float us_weapon::offense_rating() const {

	if (is_fighter) return sclasses[ammo_type].offense_rating();
	float mult(1.0);
	if (def_ammo > 0) mult *= 0.5; // limited ammo penalty
	if (const_dam)    mult *= (1.0 + 0.4*lifetime); // constant damage bonus
	if (speed == 0.0) mult *= (2.0 + 10.0*range); // instant hit and range bonus
	return mult*damage/max(1.0f, fire_delay);
}


float us_weapon::defense_rating() const {

	return (is_fighter ? sclasses[ammo_type].defense_rating() : armor);
}


// ************ US_FLEET ************


us_fleet::us_fleet(string const &name_, unsigned align_, unsigned ai_, unsigned targ_, float spread_,
				   point const &pos_, unsigned counts[], unsigned multiplier) 
				   : name(name_), align(align_), ai(ai_), targ(targ_), spread(spread_), pos(pos_), flagship(NULL)
{
	assert(multiplier > 0);

	for (unsigned i = 0; i < NUM_US_CLASS; ++i) {
		ships.push_back(make_pair(i, multiplier*counts[i]));
	}
}


void us_fleet::set_flagship(unsigned sclass, float child_stray_dist) {

	flagship = add_ship(sclass, align, ai, targ, pos, spread);
	flagship->make_flagship(child_stray_dist);
}


void us_fleet::spawn() {

	for (unsigned i = 0; i < ships.size(); ++i) {
		for (unsigned j = 0; j < ships[i].second; ++j) {
			u_ship *ship(add_ship(ships[i].first, align, ai, targ, pos, spread));
			assert(ship != NULL);
			if (flagship) ship->set_parent(flagship);
		}
	}
}


void ship_cluster::update() { // unused

	center   = all_zeros;
	velocity = zero_vector;

	for (unsigned i = 0; i < ships.size(); ++i) {
		assert(ships[i] != NULL);

		if (ships[i]->invalid() || ships[i]->get_align() != align) { // remove this ship from the cluster
			swap(ships[i], ships.back());
			ships.pop_back();
			--i;
			continue;
		}
		center   += ships[i]->get_pos();
		velocity += ships[i]->get_velocity();
	}
	if (!ships.empty()) {
		center   /= ships.size();
		velocity /= ships.size();
	}
}



