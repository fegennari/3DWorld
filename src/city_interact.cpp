// 3D World - Player Interaction with Cities and Spectating Building AI/Pedestrians/Cars
// by Frank Gennari
// 11/23/21
#include "city.h"

extern bool camera_in_building, enable_mouse_look, no_tt_footsteps;
extern int spectate;
extern float CAMERA_RADIUS;
extern point surface_pos;
extern city_params_t city_params;


// spectate mode

class city_spectate_manager_t {
	enum {FOLLOW_NONE, FOLLOW_BAI, FOLLOW_PED, FOLLOW_CAR};
	unsigned spectate_mode;
	unsigned follow_ix; // building AI, pedestrian, or car index
	unsigned follow_id; // unique ID
	car_manager_t *car_manager;
	ped_manager_t *ped_manager;

	template<typename T> bool match_id(T const &actor) const {return (actor.get_unique_id() == follow_id);}

	int find_closest_person(vector<pedestrian_t> const &peds, point const &pos) const {
		float const dmax(4.0*CAMERA_RADIUS); // max distance for spectate
		float dmin_sq(dmax*dmax);
		int closest_ix(-1);
		
		for (auto i = peds.begin(); i != peds.end(); ++i) {
			float const dsq(p2p_dist_sq(pos, i->pos));
			if (dsq < dmin_sq) {closest_ix = (i - peds.begin()); dmin_sq = dsq;}
		}
		return closest_ix;
	}
	int find_closest_car(vector<car_t> const &cars, point const &pos) const {
		float const dmax(4.0*CAMERA_RADIUS); // max distance for spectate
		float dmin_sq(dmax*dmax);
		int closest_ix(-1);

		for (auto i = cars.begin(); i != cars.end(); ++i) {
			if (i->is_parked()) continue; // skip parked cars; should never get here?
			float const dsq(p2p_dist_sq(pos, i->bcube.closest_pt(pos)));
			if (dsq < dmin_sq) {closest_ix = (i - cars.begin()); dmin_sq = dsq;}
		}
		return closest_ix;
	}
	void set_camera_to_follow_person(vector<pedestrian_t> const &people) const {
		assert(unsigned(follow_ix) < people.size());
		pedestrian_t const &person(people[follow_ix]);
		surface_pos = person.get_eye_pos();
		spectate    = !enable_mouse_look; // 'V' key contols whether or not the player can move the camera
		if (spectate) {cview_dir = person.dir;}
	}
	template<typename T> bool update_ix_for_correct_agent(vector<T> const &agents) {
		assert(unsigned(follow_ix) < agents.size());
		unsigned ntest(1);
		if (match_id(agents[follow_ix])) {return 1;} // done/correct
		unsigned const search_dist = 100; // distance to search in backwards as first pass (optimization)
		auto search_start(agents.begin() + ((follow_ix < search_dist) ? 0 : (follow_ix - search_dist)));

		for (auto i = search_start; i != agents.end(); ++i, ++ntest) { // start searching around the old index
			if (match_id(*i)) {follow_ix = (i - agents.begin()); return 1;} // success
		}
		for (auto i = agents.begin(); i != search_start; ++i, ++ntest) { // search the remaining range starting at the first element
			if (match_id(*i)) {follow_ix = (i - agents.begin()); return 1;} // success
		}
		return 0; // agent with matching ID not found
	}
public:
	city_spectate_manager_t() : spectate_mode(FOLLOW_NONE), follow_ix(0), follow_id(0), car_manager(nullptr), ped_manager(nullptr) {}

	void init(car_manager_t &car_manager_, ped_manager_t &ped_manager_) {
		ped_manager = &ped_manager_;
		car_manager = &car_manager_;
	}
	void clear() {
		follow_ix       = follow_id = 0;
		spectate_mode   = FOLLOW_NONE;
		spectate        = 0;
		no_tt_footsteps = 0;
	}
	void toggle_enabled() {
		if (world_mode != WMODE_INF_TERRAIN) return;
		if (spectate_mode != FOLLOW_NONE) {clear(); return;} // stop spectating
		point const camera_bs(get_camera_building_space());

		if (camera_in_building) {
			if (ped_manager && city_params.num_building_peds > 0) {
				int const ix(find_closest_person(ped_manager->peds_b, camera_bs));
				if (ix >= 0) {
					assert(unsigned(ix) < ped_manager->peds_b.size());
					//if (ped_manager->peds_b[ix].dest_bldg != player_building) break; // wrong building; doesn't work, but the dmax test should be good enough
					follow_ix = ix;
					follow_id = ped_manager->peds_b[ix].get_unique_id();
					spectate_mode = FOLLOW_BAI;
				}
			}
			return;
		}
		if (!have_cities()) return; // no need to check peds or cars
		float dmin_sq(0.0);

		if (ped_manager && city_params.num_peds > 0) {
			int const ix(find_closest_person(ped_manager->peds, camera_bs));
			if (ix >= 0) {
				assert(unsigned(ix) < ped_manager->peds.size());
				dmin_sq   = p2p_dist_sq(camera_bs, ped_manager->peds[ix].pos);
				follow_ix = ix;
				follow_id = ped_manager->peds[ix].get_unique_id();
				spectate_mode = FOLLOW_PED;
			}
		}
		if (car_manager && city_params.num_cars > 0) {
			int const ix(find_closest_car(car_manager->cars, camera_bs));
			if (ix >= 0) {
				assert(unsigned(ix) < car_manager->cars.size());
				car_t const &car(car_manager->cars[ix]);

				if (dmin_sq == 0.0 || p2p_dist_sq(camera_bs, car.get_center()) < dmin_sq) { // no ped, or car is closer
					follow_ix = ix;
					follow_id = car.get_unique_id(); // should be constant and pseudo-unique across cars
					spectate_mode = FOLLOW_CAR;
				}
			}
		}
	}
	void next_frame() {
		if (world_mode != WMODE_INF_TERRAIN) {clear(); return;} // not spectating

		switch (spectate_mode) {
		case FOLLOW_NONE:
			return;
		case FOLLOW_BAI: {
			assert(ped_manager);
			set_camera_to_follow_person(ped_manager->peds_b); // index never changes (no peds_b sort), don't need to update
			break;
		}
		case FOLLOW_PED: {
			assert(ped_manager);
			if (!update_ix_for_correct_agent(ped_manager->peds)) {assert(0);} // should never fail
			assert(unsigned(follow_ix) < ped_manager->peds.size());
			if (ped_manager->peds[follow_ix].at_dest) {clear(); return;} // stop following when ped reaches the destination, before respawn
			set_camera_to_follow_person(ped_manager->peds);
			break;
		}
		case FOLLOW_CAR: {
			assert(car_manager);
			if (!update_ix_for_correct_agent(car_manager->cars)) {clear(); return;} // if this fails, maybe it's a car that parked in a driveway
			assert(unsigned(follow_ix) < car_manager->cars.size());
			car_t const &car(car_manager->cars[follow_ix]);
			point const center(car.get_center());
			surface_pos.assign(center.x, center.y, (center.z + 0.25*car.height)); // 75% of car height
			spectate = !enable_mouse_look; // 'V' key contols whether or not the player can move the camera
			
			if (spectate) { // calculate direction using code copied from car_draw_state_t::draw_car()
				point pb[8], pt[8]; // bottom and top sections
				car_draw_state_t::gen_car_pts(car, 0, pb, pt); // draw_top=0
				cview_dir = cross_product((pb[5] - pb[1]), (pb[0] - pb[1])).get_norm() * ((car.dim ^ car.dir) ? -1.0 : 1.0);
			}
			no_tt_footsteps = 1; // not walking
			break;
		}
		default: // undefined mode
			assert(0);
		}
		surface_pos += get_tiled_terrain_model_xlate(); // convert back to camera space
	}
	// functions used to disable drawing of the actor the player is following
	bool skip_bai_draw(pedestrian_t const &ped) const {return (spectate_mode == FOLLOW_BAI && match_id(ped));}
	bool skip_ped_draw(pedestrian_t const &ped) const {return (spectate_mode == FOLLOW_PED && match_id(ped));}
	bool skip_car_draw(car_t        const &car) const {return (spectate_mode == FOLLOW_CAR && match_id(car));} // almost exact ; id collision is rare
};

city_spectate_manager_t city_spectate_manager;

void init_city_spectate_manager(car_manager_t &car_manager, ped_manager_t &ped_manager) {city_spectate_manager.init(car_manager, ped_manager);}
void toggle_city_spectate_mode() {city_spectate_manager.toggle_enabled();} // key F8
void follow_city_actor() {city_spectate_manager.next_frame();}
bool skip_bai_draw(pedestrian_t const &bai) {return city_spectate_manager.skip_bai_draw(bai);}
bool skip_ped_draw(pedestrian_t const &ped) {return city_spectate_manager.skip_ped_draw(ped);}
bool skip_car_draw(car_t        const &car) {return city_spectate_manager.skip_car_draw(car);}

