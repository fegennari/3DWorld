// 3D World - Building Points of Interest for AI
// by Frank Gennari 11/22/26

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for person_t
#include "openal_wrap.h"


extern building_t const *player_building;


unsigned point_of_interest_t::get_stand_areas(cube_t const &room_bounds, cube_t sas[4], float radius) const {
	cube_t avoid_area(look_area), clip_area(room_bounds);
	avoid_area.expand_by_xy( radius);
	clip_area .expand_by_xy(-radius);
	unsigned num(0);

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			float const act_edge(act_area.d[dim][dir]), avoid_edge(avoid_area.d[dim][dir]);
			if (act_edge == avoid_edge || ((avoid_edge < act_edge) ^ dir)) continue; // active edge doesn't extend beyond avoid edge
			cube_t &sa(sas[num]);
			sa = act_area;
			sa.d[dim][!dir] = avoid_edge;
			// trim off the far edge if there's space to put the person closer to the object
			if (radius > 0.0 && sa.get_sz_dim(dim) > radius) {sa.d[dim][dir] -= (dir ? 1.0 : -1.0)*radius;}
			if (!sa.intersects(clip_area)) {cout << "invalid stand area: " << TXTS(sa) << TXTS(clip_area) << endl; continue;}
			sa.intersect_with_cube(clip_area);
			++num;
		} // for dir
	} // for dim
	return num;
}

float building_t::get_poi_act_dist() const {return 0.4*get_window_vspace();} // about arm's reach distance

void building_t::add_poi(cube_t const &c, cube_t const &act_area, unsigned room_id, int obj_id) {
	assert(has_room_geom());
	interior->room_geom->pois.emplace_back(c, act_area, room_id, obj_id);
}
void building_t::add_poi_xy(cube_t const &c, unsigned room_id, int obj_id) {
	cube_t act_area(c);
	act_area.expand_by_xy(get_poi_act_dist());
	add_poi(c, act_area, room_id, obj_id);
}
void building_t::add_poi_dim(cube_t const &c, unsigned room_id, bool dim, int obj_id) {
	cube_t act_area(c);
	act_area.expand_in_dim(dim, get_poi_act_dist());
	add_poi(c, act_area, room_id, obj_id);
}
void building_t::add_poi_dim_dir(cube_t const &c, unsigned room_id, bool dim, bool dir, float dscale, int obj_id) {
	cube_t act_area(c);
	act_area.d[dim][dir] += (dir ? 1.0 : -1.0)*max(dscale*get_poi_act_dist(), 0.5f*get_min_front_clearance_inc_people());
	add_poi(c, act_area, room_id, obj_id);
}
void building_t::add_poi_object(unsigned obj_id, float dscale) {
	assert(has_room_geom());
	room_object_t const &obj(get_room_object_by_index(obj_id));
	add_poi_dim_dir(obj, obj.room_id, obj.dim, obj.dir, dscale, obj_id);
}
bool building_t::room_has_poi(unsigned room_id) const {
	return (has_room_geom() && !interior->room_geom->pois.empty() && get_room(room_id).get_has_pois());
}
void building_t::get_poi_stand_areas_for_room(unsigned room_id, float radius, float zval, vect_cube_t &stand_areas) const {
	if (!room_has_poi(room_id)) return;
	cube_t const room_area(get_walkable_room_bounds(get_room(room_id)));
	cube_t sa[4];

	for (point_of_interest_t const &p : interior->room_geom->pois) {
		if (p.room_id != room_id || p.act_area.z1() > zval || p.act_area.z2() < zval) continue;
		//if (p.obj_id < 0 || get_room_object_by_index(p.obj_id).type != TYPE_TOILET) continue; // TESTING
		unsigned const num(p.get_stand_areas(room_area, sa, radius));
		for (unsigned n = 0; n < num; ++n) {stand_areas.push_back(sa[n]);} // should we be weighting based on p.use_count here?
	}
}

struct poi_cache_t {
	vect_cube_t stand_areas;
	vector<unsigned> rpool;
	building_t const *last_bldg=nullptr;
	float last_zval=0.0, last_radius=0.0;
	int last_room_id=-1;

	bool select_stand_area(building_t const &b, unsigned room_id, float radius, point &pos, rand_gen_t &rgen) {
		float const split_len(20.0*radius);

		if (&b != last_bldg || (int)room_id != last_room_id || pos.z != last_zval || radius != last_radius) { // recompute POIs
			stand_areas.clear();
			rpool.clear();
			b.get_poi_stand_areas_for_room(room_id, radius, pos.z, stand_areas);
			if (stand_areas.empty()) return 0;

			for (unsigned i = 0; i < stand_areas.size(); ++i) {
				// longer stand areas are added multiple times for higher probability
				unsigned const num_add(max(1U, min(4U, unsigned(stand_areas[i].get_size_xy().get_max_val()/split_len))));
				for (unsigned n = 0; n < num_add; ++n) {rpool.push_back(i);}
			}
			last_bldg = &b; last_room_id = room_id; last_zval = pos.z; last_radius = radius;
		}
		if (rpool.empty()) return 0; // shouldn't fail
		unsigned const ix((rpool.size() == 1) ? 0 : (rgen.rand() % rpool.size()));
		gen_xy_pos_in_cube(pos, stand_areas[rpool[ix]], rgen);
		return 1;
	}
};
poi_cache_t poi_cache;

bool building_t::select_person_poi(unsigned room_id, float radius, point &pos, rand_gen_t &rgen) const {
	// some other ideas: select least used or least recently used POI? use sink after toilet in bathroom? return false sometimes if there's only one?
	return poi_cache.select_stand_area(*this, room_id, radius, pos, rgen);
}
bool building_t::select_poi_dest(point_of_interest_t const &poi, float radius, point &pos, rand_gen_t &rgen) const {
	cube_t sa[4];
	unsigned const num(poi.get_stand_areas(get_room(poi.room_id), sa, radius));
	if (num == 0) return 0; // no valid stand area; error?
	gen_xy_pos_in_cube(pos, sa[(num == 1) ? 0 : (rgen.rand() % num)], rgen); // select a random stand area
	return 1;
}

point_of_interest_t const *building_t::maybe_interact_with_poi(person_t &person, bool arrived) { // and set look dir
	if (!arrived) { // departed
		if (!person.at_poi) return nullptr; // not at a POI, done
		person.at_poi = 0; // no longer at a POI
	}
	if (!has_room_geom() || person.cur_room < 0) return nullptr;
	float dmin_sq(0.0);
	point_of_interest_t const *poi(nullptr);

	for (point_of_interest_t const &p : interior->room_geom->pois) {
		if ((int)p.room_id != person.cur_room)   continue; // wrong room
		if (!p.act_area.contains_pt(person.pos)) continue;
		float const dist_sq(p2p_dist_xy_sq(person.pos, p.look_area.closest_pt(person.pos)));
		if (!poi || dmin_sq < dist_sq) {poi = &p; dmin_sq = dist_sq;} // choose the closest
	}
	if (poi == nullptr) return nullptr; // no POI
	
	if (arrived) { // look at the closest point if outside the POI bounds, otherwise look at the center
		point const look_pt(poi->look_area.contains_pt(person.pos) ? poi->look_area.get_cube_center() : poi->look_area.closest_pt(person.pos));
		person.look_dir.assign((look_pt.x - person.pos.x), (look_pt.y - person.pos.y));
		person.look_dir.normalize();
		person.at_poi = 1;
		++poi->use_count;
	}
	if (this != player_building) return nullptr; // not the player's building; no object interaction
	if (poi->obj_id < 0        ) return nullptr; // POI not associated with a room object
	point_of_interest_t const *next_poi(nullptr);
	room_object_t &obj(get_room_object_by_index(poi->obj_id));
	float const floor_spacing(get_window_vspace());
	bool const on_player_floor(fabs(get_camera_pos().z - person.pos.z) < 0.5*floor_spacing);

	// Note: a removed object (TYPE_BLOCKER or TYPE_DRAIN) won't match any of the cases below and will be ignored; maybe the POI should be removed?
	if (obj.is_sink_type()) {
		bool use_sink_action(on_player_floor);

		if (use_sink_action && obj.type == TYPE_KSINK) { // kitchen sink is only part of the counter
			cube_t const sink(get_sink_cube(obj));
			float const val(person.pos[!obj.dim]);
			if (val < sink.d[!obj.dim][0] || val > sink.d[!obj.dim][1]) {use_sink_action = 0;} // not in front of the sink
		}
		if (use_sink_action) {
			if (arrived) {obj.flags |=  RO_FLAG_IS_ACTIVE;} // turn water on; water doesn't accumulate in the sink, this only happens for player activate
			else         {obj.flags &= ~RO_FLAG_IS_ACTIVE;} // turn water off
		}
	}
	else if (obj.type == TYPE_TOILET || obj.type == TYPE_URINAL) { // flush toilet/urinal when departing
		if (arrived) { // close bathroom door if close; should be done earlier
			room_t const &room(get_room(poi->room_id));
			
			if (room.get_room_type_for_zval(person.pos.z, floor_spacing) == RTYPE_BATH) { // not shared men's/women's room
				float const dz_max(0.5*floor_spacing);
				vect_door_stack_t &doorways(get_doorways_for_room(room, person.pos.z));

				for (auto &ds : doorways) {
					if (ds.get_auto_close()) continue; // no need to close this door
					assert(ds.first_door_ix < interior->doors.size());

					for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
						door_t &door(interior->doors[dix]);
						if (!ds.is_same_stack(door)) break; // moved to a different stack, done
						if (!door.open || fabs(door.zc() - person.pos.z) > dz_max) continue; // skip closed doors and doors on other floors
						if (!dist_xy_less_than(door.get_cube_center(), person.pos, 1.0*floor_spacing)) continue; // too far away to close from here
						toggle_door_state(dix, 1, 0, person.pos); // player_in_this_building=1; by_player=0
					}
				} // for ds
			}
		}
		else { // departing
			if (on_player_floor) {play_flush_sound(obj);}

			// find sinks in the same room; they could be placed before or after this POI, so search from the beginning
			for (point_of_interest_t const &p : interior->room_geom->pois) {
				if (p.room_id != poi->room_id || p.obj_id < 0 || fabs(p.act_area.z1() - poi->act_area.z1()) > 0.5*floor_spacing) continue;
				room_object_t &obj2(get_room_object_by_index(p.obj_id));
				if (!obj2.is_sink_type()) continue;
				if (!next_poi || p.use_count < next_poi->use_count) {next_poi = &p;}
			}
		}
	}
	else if (obj.type == TYPE_STOVE) {
		point const local_center(local_to_camera_space(obj.get_cube_center()));

		if (arrived) { // turn a burner on
			obj.item_flags |= (1 << (person.ssn % 3)); // select a random burner; should it be one with a pot on it?
			if (on_player_floor) {gen_sound_thread_safe(SOUND_HISS, local_center, 0.25, 0.6);}
		}
		else { // turn all burners off
			obj.item_flags = 0;
			if (on_player_floor) {gen_sound_thread_safe(SOUND_CLICK, local_center, 0.75);}
		}
	}
	else if (obj.type == TYPE_VENDING || obj.type == TYPE_CHECKOUT || obj.type == TYPE_CASHREG) { // beep
		if (arrived && on_player_floor) {gen_sound_thread_safe(SOUND_BEEP, local_to_camera_space(obj.get_cube_center()), 0.25, 0.85);}
	}
	return next_poi;
}

void building_room_geom_t::remove_obj_update_poi(unsigned obj_id) {
	for (auto i = pois.begin(); i != pois.end(); ++i) {
		if (i->obj_id == (int)obj_id) {pois.erase(i); break;}
	}
}
void building_room_geom_t::move_obj_update_poi(unsigned obj_id, vector3d const &delta) { // should be XY only
	for (point_of_interest_t &p : pois) {
		if (p.obj_id != (int)obj_id) continue;
		p.look_area += delta;
		p.act_area  += delta;
		break;
	}
}

// does this go here? It's similar to points of interest
void building_t::place_stationary_people(float radius, unsigned max_people, rand_gen_t &rgen) const {
	assert(has_room_geom());
	// restaurant entrances, retail checkouts, mall checkouts and food counters, and convenience stores
	vector<person_place_t> cands(interior->room_geom->people_place); // deep copy so that we can modify it
	if (cands.size() > max_people) {vector_random_shuffle(cands, rgen);} // too many cands, choose a random subset
	unsigned num_added(0);

	for (person_place_t const &pp : cands) {
		person_t person(radius);
		person.pos = pp.pos;
		person.dir = pp.dir;
		person.is_stationary = 1;
		add_person(person, rgen);
		if (++num_added == max_people) break;
	} // for pp
}

