// 3D World - Building Animals (rats, etc.)
// by Frank Gennari 1/16/22

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h" // for object_model_loader_t

extern float fticks;
extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;


float get_rat_height(float radius) {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT)); // 3878, 861, 801
	return 2.0*radius*sz.z/max(sz.x, sz.y); // use max of x/y size; the x/y size represents the bcube across rotations
}

cube_t rat_t::get_bcube() const {
	cube_t bcube(pos, pos);
	bcube.expand_by_xy(radius);
	bcube.z2() += get_rat_height(radius);
	return bcube;
}
cube_t rat_t::get_bcube_with_dir() const {
	bool const pri_dim(fabs(dir.x) < fabs(dir.y));
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT));
	cube_t bcube(pos, pos);
	bcube.expand_in_dim( pri_dim, radius); // larger dim
	bcube.expand_in_dim(!pri_dim, radius*min(sz.x, sz.y)/max(sz.x, sz.y)); // smaller dim
	bcube.z2() += get_rat_height(radius);
	return bcube;
}

void building_t::update_animals(unsigned building_ix) {
	if (is_rotated() || !has_room_geom() || interior->rooms.empty() || global_building_params.num_rats == 0) return;
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) return; // no rat model
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, mat_ix+1); // unique per building
	vector<rat_t> &rats(interior->room_geom->rats);

	if (rats.empty()) { // new building - place rats
		float const base_radius(0.1*get_window_vspace());
		rats.reserve(global_building_params.num_rats);

		for (unsigned n = 0; n < global_building_params.num_rats; ++n) {
			float const radius(base_radius*rgen.rand_uniform(0.8, 1.2));
			point const pos(gen_rat_pos(radius, rgen));
			if (pos == all_zeros) continue; // bad pos? skip this rat
			rats.emplace_back(pos, radius);
		}
	}
	for (auto &rat : rats) {update_rat(rat, rgen);}
}

point building_t::gen_rat_pos(float radius, rand_gen_t &rgen) const {
	for (unsigned n = 0; n < 100; ++n) { // make up to 100 tries
		unsigned const room_ix(rgen.rand() % interior->rooms.size());
		room_t const &room(interior->rooms[room_ix]);
		if (room.z1() > ground_floor_z1) continue; // not on the ground floor or basement
		cube_t place_area(room); // will represent the usable floor area
		place_area.expand_by_xy(-(radius + get_wall_thickness()));
		point pos(gen_xy_pos_in_area(place_area, radius, rgen));
		pos.z = place_area.z1() + get_fc_thickness(); // on top of the floor
		if (is_valid_ai_placement(pos, radius)) {return pos;} // success
	} // for n
	return all_zeros; // failed
}

void building_t::update_rat(rat_t &rat, rand_gen_t &rgen) const {
	if (rat.speed > 0.0) {
		// TODO: movement logic, collision detection
		rat.pos += rat.speed*rat.dir;
	}
	if (rat.fear > 0.0) {
		// TODO: hide (in opposite direction from fear_pos?)
		rat.fear = max(0.0, (rat.fear - 0.1*(fticks/TICKS_PER_SECOND))); // reduce fear over 10s
	}
	else {
		// explore
	}
	if (rat.dest != rat.pos) {
		// TODO: set speed
		rat.dir = (rat.dest - rat.pos).get_norm(); // TODO: slow turn (like people)
	}
}

void building_t::scare_animals(point const &scare_pos, float sight_amt, float sound_amt) {
	vector<rat_t> &rats(interior->room_geom->rats);
	if (rats.empty()) return; // nothing to do
	float const amount(min((sight_amt + sound_amt), 1.0f)); // for now we don't differentiate between sight and sound for rats
	assert(amount > 0.0);
	float const max_scare_dist(3.0*get_window_vspace()), scare_dist(max_scare_dist*amount);
	int const scare_room(get_room_containing_pt(scare_pos));
	if (scare_room < 0) return; // error?

	for (auto &rat : rats) {
		if (rat.fear > 0.99) continue; // already max fearful (optimization)
		float const dist(p2p_dist(rat.pos, scare_pos));
		if (dist < scare_dist) continue; // optimization
		int const rat_room(get_room_containing_pt(rat.pos));
		assert(rat_room >= 0);
		if (rat_room != scare_room) continue; // only scared if in the same room
		float fear(amount);

		if (sight_amt > 0.0) { // check line of sight
			bool is_visible(1);
			// TODO: visibility check
		}
		fear = fear*max_scare_dist - dist;
		if (fear <= 0.0) continue;
		rat.fear     = min(1.0f, (rat.fear + fear));
		rat.fear_pos = scare_pos;
	} // for rat
}
