// 3D World - Player Interaction with Cities and Spectating Building AI/Pedestrians/Cars
// by Frank Gennari
// 11/23/21
#include "city.h"



extern bool camera_in_building;
extern int spectate;
extern point surface_pos;
extern city_params_t city_params;


// spectate mode

class city_spectate_manager_t {
	enum {FOLLOW_NONE, FOLLOW_BAI, FOLLOW_PED, FOLLOW_CAR};
	unsigned spectate_mode;
	unsigned follow_ix; // building AI, pedestrian, or car index
	car_manager_t *car_manager;
	ped_manager_t *ped_manager;

	int find_nearest_bai(point const &pos) const {
		return -1; // TODO
	}
	int find_nearest_ped(point const &pos) const {
		return -1; // TODO
	}
	int find_nearest_car(point const &pos) const {
		return -1; // TODO
	}
public:
	city_spectate_manager_t() : spectate_mode(FOLLOW_NONE), follow_ix(0), car_manager(nullptr), ped_manager(nullptr) {}

	void init(car_manager_t &car_manager_, ped_manager_t &ped_manager_) {
		ped_manager = &ped_manager_;
		car_manager = &car_manager_;
	}
	void toggle_enabled() {
		if (world_mode != WMODE_INF_TERRAIN) return;
		if (spectate_mode != FOLLOW_NONE) {spectate_mode = FOLLOW_NONE; return;} // stop spectating
		point const camera_bs(get_camera_building_space());

		if (camera_in_building) {
			if (ped_manager && city_params.num_building_peds > 0) {
				int const ix(find_nearest_bai(camera_bs));
				if (ix >= 0) {
					assert(unsigned(ix) < ped_manager->peds_b.size());
					follow_ix = ix;
					spectate_mode = FOLLOW_BAI;
				}
			}
			return;
		}
		if (!have_cities()) return; // no need to check peds or cars
		float dmin_sq(0.0);

		if (ped_manager && city_params.num_peds > 0) {
			int const ix(find_nearest_ped(camera_bs));
			if (ix >= 0) {
				assert(unsigned(ix) < ped_manager->peds.size());
				dmin_sq   = p2p_dist_sq(camera_bs, ped_manager->peds[ix].pos);
				follow_ix = ix;
				spectate_mode = FOLLOW_PED;
			}
		}
		else if (car_manager && city_params.num_cars > 0) {
			int const ix(find_nearest_car(camera_bs));
			if (ix >= 0) {
				assert(unsigned(ix) < car_manager->cars.size());

				if (dmin_sq == 0.0 || p2p_dist_sq(camera_bs, car_manager->cars[ix].get_center()) < dmin_sq) { // no ped, or car is closer
					follow_ix = ix;
					spectate_mode = FOLLOW_CAR;
				}
			}
		}
	}
	void next_frame() {
		if (world_mode != WMODE_INF_TERRAIN) {spectate_mode = FOLLOW_NONE; return;} // not spectating

		switch (spectate_mode) {
		case FOLLOW_NONE:
			return;
		case FOLLOW_BAI: {
			assert(ped_manager);
			assert(unsigned(follow_ix) < ped_manager->peds_b.size());
			pedestrian_t const &bai(ped_manager->peds_b[follow_ix]);
			surface_pos.assign(bai.pos.x, bai.pos.y, bai.get_height()); // TODO: eye height?
			cview_dir = bai.dir;
			break;
		}
		case FOLLOW_PED: { // FIXME: handle ped index change after sorting somehow, maybe by searching for ped by position?
			assert(ped_manager);
			assert(unsigned(follow_ix) < ped_manager->peds.size());
			pedestrian_t const &ped(ped_manager->peds[follow_ix]);
			surface_pos.assign(ped.pos.x, ped.pos.y, ped.get_height()); // TODO: eye height?
			cview_dir = ped.dir;
			break;
		}
		case FOLLOW_CAR: { // FIXME: handle car index change after sorting somehow, maybe by searching for car by position?
			assert(car_manager);
			assert(unsigned(follow_ix) < car_manager->cars.size());
			car_t const &car(car_manager->cars[follow_ix]);
			point const center(car.get_center());
			surface_pos.assign(center.x, center.y, (center.z + 0.25*car.height)); // 75% of car height
			// don't set cview_dir; allow the player to turn their head
			break;
		}
		default:
			assert(0);
		}
		// set surface_pos and maybe cview_dir
	}
	// functions used to disable drawing of the actor the player is following
	bool skip_bai_draw(unsigned ix) const {return (spectate_mode == FOLLOW_BAI && ix == follow_ix);}
	bool skip_ped_draw(unsigned ix) const {return (spectate_mode == FOLLOW_PED && ix == follow_ix);}
	bool skip_car_draw(unsigned ix) const {return (spectate_mode == FOLLOW_CAR && ix == follow_ix);}
};

city_spectate_manager_t city_spectate_manager;

void init_city_spectate_manager(car_manager_t &car_manager, ped_manager_t &ped_manager) {city_spectate_manager.init(car_manager, ped_manager);}
void toggle_city_spectate_mode() {city_spectate_manager.toggle_enabled();} // key F8
void follow_city_actor() {city_spectate_manager.next_frame();}

