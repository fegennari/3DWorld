// 3D World - Player Interaction with Cities and Spectating Building AI/Pedestrians/Cars
// by Frank Gennari
// 11/23/21
#include "city.h"

extern bool camera_in_building, enable_mouse_look;
extern int spectate;
extern float CAMERA_RADIUS;
extern point surface_pos;
extern city_params_t city_params;


// spectate mode

class city_spectate_manager_t {
	enum {FOLLOW_NONE, FOLLOW_BAI, FOLLOW_PED, FOLLOW_CAR};
	unsigned spectate_mode;
	unsigned follow_ix; // building AI, pedestrian, or car index
	unsigned person_ssn; // for building AI and pedestrians
	float car_speed; // pseudo-unique value per-car
	car_manager_t *car_manager;
	ped_manager_t *ped_manager;

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
	bool update_ix_for_correct_person(vector<pedestrian_t> const &people) {
		assert(unsigned(follow_ix) < people.size());
		if (people[follow_ix].ssn == person_ssn) return 1; // done/correct

		for (auto i = people.begin(); i != people.end(); ++i) { // slow, could maybe start searching around follow_ix
			if (i->ssn != person_ssn) continue; // wrong person
			follow_ix = (i - people.begin());
			return 1; // success
		}
		return 0; // person with matching SSN not found, error? or maybe reached dest?
	}
	void set_camera_to_follow_person(vector<pedestrian_t> const &people) const {
		assert(unsigned(follow_ix) < people.size());
		pedestrian_t const &person(people[follow_ix]);
		surface_pos = person.get_eye_pos();
		spectate    = !enable_mouse_look; // 'V' key contols whether or not the player can move the camera
		if (spectate) {cview_dir = person.dir;}
	}
	bool update_ix_for_correct_car() {
		assert(car_manager);
		assert(unsigned(follow_ix) < car_manager->cars.size());
		if (car_manager->cars[follow_ix].max_speed == car_speed) return 1; // done/correct

		for (auto i = car_manager->cars.begin(); i != car_manager->cars.end(); ++i) { // slow, could maybe start searching around follow_ix
			if (i->max_speed != car_speed) continue; // wrong car
			follow_ix = (i - car_manager->cars.begin());
			return 1; // success
		}
		return 0; // car with matching speed not found, maybe it parked
	}
public:
	city_spectate_manager_t() : spectate_mode(FOLLOW_NONE), follow_ix(0), person_ssn(0), car_speed(0.0), car_manager(nullptr), ped_manager(nullptr) {}

	void init(car_manager_t &car_manager_, ped_manager_t &ped_manager_) {
		ped_manager = &ped_manager_;
		car_manager = &car_manager_;
	}
	void clear() {
		follow_ix = person_ssn = 0;
		car_speed = 0.0;
		spectate_mode = FOLLOW_NONE;
		spectate = 0;
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
					follow_ix  = ix;
					person_ssn = ped_manager->peds_b[ix].ssn;
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
				dmin_sq    = p2p_dist_sq(camera_bs, ped_manager->peds[ix].pos);
				follow_ix  = ix;
				person_ssn = ped_manager->peds[ix].ssn;
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
					car_speed = car.max_speed; // should be constant and pseudo-unique across cars
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
			set_camera_to_follow_person(ped_manager->peds_b);
			break;
		}
		case FOLLOW_PED: {
			assert(ped_manager);
			if (!update_ix_for_correct_person(ped_manager->peds)) {clear(); return;}
			set_camera_to_follow_person(ped_manager->peds);
			break;
		}
		case FOLLOW_CAR: {
			assert(car_manager);
			if (!update_ix_for_correct_car()) {clear(); return;}
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
			break;
		}
		default: // undefined mode
			assert(0);
		}
	}
	// functions used to disable drawing of the actor the player is following
	//bool skip_bai_draw(unsigned ix) const {return (spectate_mode == FOLLOW_BAI && ix == follow_ix);}
	//bool skip_ped_draw(unsigned ix) const {return (spectate_mode == FOLLOW_PED && ix == follow_ix);}
	//bool skip_car_draw(unsigned ix) const {return (spectate_mode == FOLLOW_CAR && ix == follow_ix);}
	bool skip_bai_draw(pedestrian_t const &ped) const {return (spectate_mode == FOLLOW_BAI && ped.ssn == person_ssn);}
	bool skip_ped_draw(pedestrian_t const &ped) const {return (spectate_mode == FOLLOW_PED && ped.ssn == person_ssn);}
	bool skip_car_draw(car_t        const &car) const {return (spectate_mode == FOLLOW_CAR && car.max_speed == car_speed);} // not exact enough?
};

city_spectate_manager_t city_spectate_manager;

void init_city_spectate_manager(car_manager_t &car_manager, ped_manager_t &ped_manager) {city_spectate_manager.init(car_manager, ped_manager);}
void toggle_city_spectate_mode() {city_spectate_manager.toggle_enabled();} // key F8
void follow_city_actor() {city_spectate_manager.next_frame();}
bool skip_bai_draw(pedestrian_t const &bai) {return city_spectate_manager.skip_bai_draw(bai);}
bool skip_ped_draw(pedestrian_t const &ped) {return city_spectate_manager.skip_ped_draw(ped);}
bool skip_car_draw(car_t        const &car) {return city_spectate_manager.skip_car_draw(car);}

