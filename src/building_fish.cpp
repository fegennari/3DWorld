// 3D World - Building Fish
// by Frank Gennari 6/23/24

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h" // for animation_state_t
//#include "profiler.h"

extern int animate2, frame_counter, player_in_water;
extern float fticks;

bool check_ramp_collision(room_object_t const &c, point &pos, float radius, vector3d *cnorm);
void draw_animated_fish_model(shader_t &s, vector3d const &pos, float radius, vector3d const &dir, float anim_time, colorRGBA const &color);
bool play_attack_sound(point const &pos, float gain, float pitch, rand_gen_t &rgen);
void draw_bubbles(vector<sphere_t> const &bubbles, shader_t &s, bool in_fishtank);


int get_future_frame(float min_secs, float max_secs, rand_gen_t &rgen) {
	float const secs((min_secs == max_secs) ? min_secs : rgen.rand_uniform(min_secs, max_secs));
	return (frame_counter + round_fp(secs*TICKS_PER_SECOND));
}

class fish_manager_t {
	struct fish_t {
		float radius=0.0, mspeed=0.0, tspeed=0.0, speed_mult=1.0; // movement speed and tail speed
		bool attacks_player=0;
		unsigned id=0;
		int next_splash_frame=0, next_turn_frame=0, next_speed_frame=0, next_alert_frame=0;
		point pos; // p_last?
		vector3d dir, target_dir; // Note: fish only rotated about Z and remains in the XY plane
		colorRGBA color=WHITE;

		void draw(shader_t &s, animation_state_t &anim_state, float anim_time) const {
			anim_state.anim_time = 0.5*sqrt(speed_mult)*tspeed*anim_time; // tail nominally moves twice per second
			if (anim_state.enabled) {anim_state.set_animation_id_and_time(s, 0, 1.0);}
			draw_animated_fish_model(s, pos, radius, dir, anim_state.anim_time, color);
		}
		bool can_splash(rand_gen_t &rgen) {
			if (frame_counter < next_splash_frame) return 0;
			next_splash_frame = get_future_frame(1.0, 3.0, rgen); // once every 2s on average
			return 1;
		}
		bool has_target() const {return (target_dir != zero_vector);}
	};

	class fish_cont_t {
	protected:
		cube_t bcube, valid_area;
		rand_gen_t rgen;
		vector<fish_t> fish;
	public:
		bool present=0;

		~fish_cont_t() {}
		bool empty() const {return fish.empty();}
		cube_t const &get_bcube() const {return bcube;}
		
		void clear() {
			fish.clear();
			bcube.set_to_zeros();
			present = 0;
		}
		void init(cube_t const &bcube_, cube_t const &valid_area_, int rseed, unsigned min_num, unsigned max_num) {
			bcube      = bcube_;
			valid_area = valid_area_;
			rgen.set_state(rseed+1, rseed+1);
			for (unsigned n = 0; n < 2; ++n) {rgen.rand_mix();}
			fish.resize(rgen.rand_uniform_uint(min_num, max_num));
		}
		void populate(float max_fish_radius, float speed_mult) {
			assert(max_fish_radius > 0.0);
			unsigned id(0);

			for (fish_t &f : fish) {
				f.radius = max_fish_radius*rgen.rand_uniform(0.67, 1.0);
				f.dir    = assign_fish_dir();
				f.tspeed = rgen.rand_uniform(0.8, 1.2)/TICKS_PER_SECOND;
				f.mspeed = speed_mult*f.tspeed; // movement speed is also modulated by the tail animation speed
				f.id     = id++;
				f.attacks_player = rgen.rand_bool(); // 50% of the time; only applies to swimming pool and flooded basement fish
				for (unsigned d = 0; d < 3; ++d) {f.color[d] *= rgen.rand_uniform(0.8, 1.0);} // slight color variation
				cube_t const valid_area(get_valid_area(1.05*f.radius)); // slightly larger radius so that we don't start out intersecting

				for (unsigned n = 0; n < 100; ++n) { // 100 tries
					gen_xyz_pos_in_cube(f.pos, valid_area, rgen);
					point coll_pos; // unused
					if (!check_fish_coll(f.pos, f.radius, f.id, coll_pos)) break; // success
				}
			} // for f
		}
		void next_frame(float speed_mult_max) {
			if (!animate2) return;
			float const delta_dir(0.6*(1.0 - pow(0.7f, fticks)));

			for (fish_t &f : fish) {
				if (f.has_target()) { // fish is turning in place
					f.dir = delta_dir*f.target_dir + (1.0 - delta_dir)*f.dir;
					if (f.dir == zero_vector) {f.dir = assign_fish_dir();} // error?
					else {f.dir.normalize();}
					if (dot_product(f.dir, f.target_dir) < 0.95) continue; // still turning
					f.target_dir = zero_vector; // turn complete
				}
				cube_t const valid_area(get_valid_area(f.radius));
				assert(valid_area.contains_pt(f.pos));
				point const prev_pos(f.pos);
				f.pos += (f.mspeed*f.speed_mult*fticks)*f.dir;
				vector3d coll_dir; // fish => target
				point coll_pos;

				if (!valid_area.contains_pt(f.pos)) { // hit the container edge
					coll_dir = valid_area.closest_side_dir(prev_pos);
				}
				else if (check_fish_coll(f.pos, 0.7*f.radius, f.id, coll_pos)) { // hit another fish or a collider (use a smaller radius)
					coll_dir = (coll_pos - f.pos).get_norm();
				}
				if (coll_dir != zero_vector) {
					f.pos        = prev_pos;
					f.target_dir = assign_fish_dir();
					if (dot_product(f.target_dir, coll_dir) > 0.0) {f.target_dir.negate();}
					f.next_alert_frame = get_future_frame(1.0, 2.0, rgen); // ignore the player for a bit
				}
				else { // no collision
					if (frame_counter >= f.next_turn_frame) { // time to turn
						f.target_dir      = assign_fish_dir();
						f.next_turn_frame = get_future_frame(10.0, 30.0, rgen);
					}
					if (frame_counter >= f.next_speed_frame) { // time to change speed
						f.speed_mult      += rgen.rand_uniform(-0.5, 0.5);
						f.next_speed_frame = get_future_frame(5.0, 10.0, rgen);
					}
				}
				f.speed_mult = max(0.5f, min(speed_mult_max, f.speed_mult)); // clamp to a reasonable range
			} // for f
		}
		void draw(shader_t &s, animation_state_t &anim_state, float anim_time) const {
			for (fish_t const &f : fish) {f.draw(s, anim_state, anim_time);}
		}
	protected:
		vector3d assign_fish_dir() { // match the aspect ratio of the container, smaller in Z, but limit zval aspect ratio to 10%
			vector3d const sz(bcube.get_size()), dir_weight(sz.x, sz.y, 0.5*max(sz.z, 0.1f*max(sz.x, sz.y)));
			return (rgen.signed_rand_vector_spherical() * dir_weight).get_norm();
		}
		cube_t get_valid_area(float pad) const {
			cube_t va(valid_area);
			va.expand_by(-pad);
			assert(va.is_strictly_normalized());
			return va;
		}
		virtual bool check_fish_coll(point const &pos, float radius, unsigned id, point &coll_pos) const { // default fish-fish collision check
			for (fish_t const &f : fish) {
				if (f.radius == 0.0) continue; // not yet setup
				if (f.id >= id     ) continue; // skip ourself and higher ID fish
				if (dist_less_than(pos, f.pos, (radius + f.radius))) {coll_pos = f.pos; return 1;}
			}
			return 0;
		}
	}; // fish_cont_t

	class fishtank_t : public fish_cont_t {
		vector<sphere_t> bubbles;
	public:
		unsigned obj_id=0; // used as a unique identifier
		bool visible=0;

		fishtank_t(room_object_t const &obj, bool visible_=0) : obj_id(obj.obj_id), visible(visible_) {
			present = 1;
			cube_t tank_inner(obj);
			tank_inner.expand_by(-vector3d(0.04, 0.04, 0.1)*obj.get_height()); // subtract off the glass and some area from the top and bottom
			init(obj, tank_inner, obj_id, 1, 4); // 1-4 fish
			float const max_fish_radius(0.125*(1.0 + 1.0/fish.size())*bcube.min_len()); // more fish = smaller size
			populate(max_fish_radius, 0.0016);
		}
		void clear() {
			fish_cont_t::clear();
			bubbles.clear();
		}
		void next_frame() {
			present = visible = 0; // mark as not present or visible until it's seen
			fish_cont_t::next_frame(1.5); // speed_mult_max=1.5
			if (!animate2) return;
			// update bubbles
			float const elapsed_secs(fticks/TICKS_PER_SECOND), bubble_rate(elapsed_secs/4.0), rise_rate(0.25*elapsed_secs*bcube.dz()); // avg one bubble per 4s

			if (bubbles.size() < fish.size()) { // don't allow more bubbles than fish
				for (fish_t const &f : fish) {
					if (rgen.rand_float() > bubble_rate) continue;
					vector3d const bdir(vector3d(f.dir.x, f.dir.y, 0.0).get_norm()); // XY
					bubbles.emplace_back((f.pos + 1.4*f.radius*bdir), rgen.rand_uniform(0.06, 0.12)*f.radius); // spawn bubble in front
				}
			}
			auto i(bubbles.begin()), o(i);

			for (; i != bubbles.end(); ++i) {
				i->pos.z += rise_rate;
				if (i->pos.z + 0.5*i->radius < valid_area.z2()) {*o++ = *i;} // remove bubble when it reaches the water surface
			}
			bubbles.erase(o, bubbles.end());
		}
		void update_object(room_object_t const &obj) { // handle movement when the table is pushed
			present = 1;
			if (obj == bcube) return; // no update
			vector3d const delta(obj.get_llc() - bcube.get_llc());
			bcube = obj;
			valid_area += delta;
			for (fish_t &f : fish) {f.pos += delta;}
		}
		void draw(shader_t &s, animation_state_t &anim_state, float anim_time) const {
			if (!visible) return;
			fish_cont_t::draw(s, anim_state, anim_time);
			draw_bubbles(bubbles, s, 1); // in_fishtank=1
		}
	}; // fishtank_t

	class area_fish_cont_t : public fish_cont_t { // base class for swimming pools and flooded basements
	protected:
		vect_cube_t obstacles;
		bool test_line_of_sight=0;
	public:
		void next_frame() {
			vector3d const xlate(get_camera_coord_space_xlate());
			fish_cont_t::next_frame(3.0); // speed_mult_max=3.0

			for (fish_t &f : fish) { // handle water splashes if fish is visible
				if (f.pos.z < valid_area.z2() - 2.0*f.radius) continue; // too low to splash
				if (f.can_splash(rgen) && camera_pdu.sphere_visible_test((f.pos + xlate), f.radius)) {add_water_splash(f.pos, 2.0*f.radius, 0.25);}
			}
			if (player_in_water) { // handle player interaction
				float const player_radius(building_t::get_scaled_player_radius()), player_height(get_bldg_player_height());
				point const camera_bs(get_camera_pos() - xlate);

				for (fish_t &f : fish) {
					if (f.has_target()) continue; // skip if recently collided
					if (frame_counter < f.next_alert_frame) continue; // distracted by a collision; required to avoid getting stuck against a wall
					if (f.pos.z < camera_bs.z - player_height || f.pos.z > camera_bs.z + player_radius) continue; // no Z overlap

					if (f.attacks_player && in_building_gameplay_mode()) {
						if (!dist_xy_less_than(camera_bs, f.pos, 10.0*(CAMERA_RADIUS + f.radius))) continue; // not close enough to attack
						float const coll_dist(player_radius + f.radius);
						point const test_pt(camera_bs.x, camera_bs.y, f.pos.z); // at fish zval

						if (dist_xy_less_than(test_pt, f.pos, coll_dist)) { // collided with the player
							f.pos = test_pt + coll_dist*(f.pos - test_pt).get_norm(); // move to not collide
							get_valid_area(f.radius).clamp_pt(f.pos);
							play_attack_sound((f.pos + get_camera_coord_space_xlate()), 1.0, 1.0, rgen);
							bool const player_dead(player_take_damage(0.001)); // small damage
							if (player_dead) {register_achievement("Fish Food");} // damage over time; achievement if the player dies
							// Note: no blood decal or dropping of inventory items in the water
						}
						bool visible(1); // check for player visibilty

						for (cube_t const &c : obstacles) {
							if (c.line_intersects(camera_bs, f.pos)) {visible = 0; break;}
						}
						if (!visible) continue;
						f.target_dir = (test_pt - f.pos).get_norm(); // move toward the player
						max_eq(f.speed_mult, 2.0f); // faster
					}
					else { // avoid player
						float const alert_dist(4.0*(CAMERA_RADIUS + f.radius));
						if (!dist_xy_less_than(camera_bs, f.pos, alert_dist)) continue; // not close enough
						point const test_pt(camera_bs.x, camera_bs.y, f.pos.z); // at fish zval
						if (dot_product_ptv(f.dir, test_pt, f.pos) < 0.0) continue; // moving away from the player
						if (test_line_of_sight && line_int_cubes(test_pt, f.pos, obstacles, cube_t(camera_bs, f.pos))) continue; // not visible
						f.target_dir = (f.pos - test_pt).get_norm(); // move away from the player
						max_eq(f.speed_mult, 5.0f*(1.2f - p2p_dist_xy(test_pt, f.pos)/alert_dist)); // higher speed
					}
				} // for f
			}
		}
		virtual bool check_fish_coll(point const &pos, float radius, unsigned id, point &coll_pos) const {
			if (fish_cont_t::check_fish_coll(pos, radius, id, coll_pos)) return 1;

			for (cube_t const &c : obstacles) {
				if (sphere_cube_intersect(pos, 1.8*radius, c)) {coll_pos = c.get_cube_center(); return 1;} // extra radius for clearance
			}
			return 0;
		}
		void draw(shader_t &s, animation_state_t &anim_state, float anim_time) const { // override fish_cont_t::draw() to include VFC
			vector3d const xlate(get_camera_coord_space_xlate());

			for (fish_t const &f : fish) {
				if (camera_pdu.sphere_visible_test((f.pos + xlate), f.radius)) {f.draw(s, anim_state, anim_time);} // VFC, assuming non-rotated building
			}
		}
	}; // player_int_fish_cont_t

	class swimming_pool_t : public area_fish_cont_t {
		room_object_t pool_ramp;
	public:
		void init(indoor_pool_t const &pool, room_object_t const &pool_ramp_, cube_t const &water_bcube, vect_cube_t const &obstacles_, int rseed) {
			obstacles = obstacles_;
			pool_ramp = pool_ramp_;
			present   = 1;
			fish_cont_t::init(pool, water_bcube, rseed, 2, 5); // 2-5 fish
			float const max_fish_radius(0.05*min(min(pool.dx(), pool.dy()), pool.dz()));
			populate(max_fish_radius, 0.003);
		}
		virtual bool check_fish_coll(point const &pos, float radius, unsigned id, point &coll_pos) const {
			if (area_fish_cont_t::check_fish_coll(pos, radius, id, coll_pos)) return 1;
			point cpos(pos); // since we can't modify pos
			
			if (!pool_ramp.is_all_zeros() && check_ramp_collision(pool_ramp, cpos, radius, nullptr)) { // check pool ramp
				coll_pos = pos - radius*plus_z; // collision is always below
				return 1;
			}
			return 0;
		}
	}; // swimming_pool_t

	class flooded_basement_t : public area_fish_cont_t {
		unsigned fixed_obstacles_end_ix=0;
	public:
		void init(cube_t const &water_bc, float floor_spacing, vect_cube_t const &obstacles_, int rseed) {
			obstacles = obstacles_;
			present   = 1;
			fixed_obstacles_end_ix = obstacles.size();
			test_line_of_sight     = 1;
			fish_cont_t::init(water_bc, water_bc, rseed, 10, 20); // 10-20 fish
			float const max_fish_radius(0.25*min(bcube.dz(), 0.25f*floor_spacing));
			populate(max_fish_radius, 0.004);
		}
		void register_doors_state(vect_door_stack_t const &door_stacks, vect_door_t const &doors) {
			if (!present || fish.empty()) return;
			assert(fixed_obstacles_end_ix <= obstacles.size());
			obstacles.resize(fixed_obstacles_end_ix); // remove any previously added doors

			for (door_stack_t const &ds : door_stacks) { // add doors
				if (!valid_area.intersects(ds)) continue;
				assert(ds.first_door_ix < doors.size()); // should be only one door
				if (doors[ds.first_door_ix].open_amt < 1.0) {obstacles.push_back(doors[ds.first_door_ix]);} // add if not fully open
			}
		}
	}; // flooded_basement_t

	vector<fishtank_t> fishtanks;
	swimming_pool_t swimming_pool;
	flooded_basement_t flooded_basement;
	building_t const *prev_building=nullptr;
	unsigned fishtank_ix=0, rseed_ix=0;
	float anim_time=0.0;
	bool had_warning=0;
public:
	bool empty() const {return (fishtanks.empty() && swimming_pool.empty() && flooded_basement.empty());}
	bool has_swimming_pool   () const {return swimming_pool   .present;}
	bool has_flooded_basement() const {return flooded_basement.present;}

	void clear() {
		fishtanks       .clear();
		swimming_pool   .clear();
		flooded_basement.clear();
	}
	void next_frame(building_t const &building) {
		if (!building.interior) return; // error?

		if (&building != prev_building || building.interior->rgen_seed_ix != rseed_ix) { // new building or new seed
			clear();
			prev_building = &building;
			rseed_ix      = building.interior->rgen_seed_ix;
			anim_time     = 0.0; // reset
		}
		for (fishtank_t &ft : fishtanks) {ft.next_frame();}
		swimming_pool   .next_frame();
		flooded_basement.next_frame();
		if (animate2) {anim_time += fticks;}
		fishtank_ix = 0;
	}
	void register_fishtank(room_object_t const &obj, bool is_visible) {
		if (fishtanks.size() <= fishtank_ix) { // fishtank not yet added
			fishtanks.emplace_back(obj, is_visible);
		}
		else { // find existing fishtank
			bool found(0);

			for (unsigned pass = 0; pass < 2 && !found; ++pass) { // first try to match bcube, then try to match obj_id + size to handle movement
				for (fishtank_t &ft : fishtanks) {
					if ((pass == 0) ? (ft.get_bcube() != obj) : (ft.obj_id != obj.obj_id || ft.get_bcube().get_size() != obj.get_size())) continue;
					assert(!ft.present);
					ft.update_object(obj);
					ft.visible = is_visible;
					found = 1;
					break; // can only have one
				} // for ft
			} // for pass
			if (!found && !had_warning) {cout << "Failed to find previous fishtank at " << obj.str() << endl; had_warning = 1;} // print once
		}
		++fishtank_ix;
	}
	void register_swimming_pool(indoor_pool_t const &pool, room_object_t const &pool_ramp, cube_t const &water_bcube, vect_cube_t const &obstacles, int rseed) {
		if (!has_swimming_pool()) {swimming_pool.init(pool, pool_ramp, water_bcube, obstacles, rseed);}
	}
	void register_flooded_basement(cube_t const &water_bc, float floor_spacing, vect_cube_t const &obstacles, int rseed) {
		if (!has_flooded_basement()) {flooded_basement.init(water_bc, floor_spacing, obstacles, rseed);}
	}
	void register_doors_state(vect_door_stack_t const &door_stacks, vect_door_t const &doors) {
		if (has_flooded_basement()) {flooded_basement.register_doors_state(door_stacks, doors);}
	}
	void end_objs() { // remove any fishtanks that are no longer present
		unsigned const pre_num_fishtanks(fishtanks.size());
		fishtanks.erase(remove_if(fishtanks.begin(), fishtanks.end(), [](fishtank_t const &ft) {return !ft.present;}), fishtanks.end());
		
		if (fishtank_ix != fishtanks.size()) {
			cout << "Error: Fishtanks changed for building: " << TXT(pre_num_fishtanks) << TXT(fishtanks.size()) << TXT(fishtank_ix);
			if (prev_building) {cout << " name: " << prev_building->name;}
			cout << endl;
			fishtanks.clear(); // fishtanks changed due to nondeterministic building regen; rebuild on the next frame
		}
	}
	void draw_fish(shader_t &s, bool inc_pools_and_fb) const {
		if (empty() || (!inc_pools_and_fb && fishtanks.empty())) return;
		glDisable(GL_CULL_FACE); // fix for looking through the fish's mouth
		animation_state_t anim_state(1, ANIM_ID_FISH_TAIL); // enabled=1; animation_scale and model_delta_height are not set or used
		for (fishtank_t const &f : fishtanks) {f.draw(s, anim_state, anim_time);}
		if (inc_pools_and_fb) {swimming_pool    .draw(s, anim_state, anim_time);}
		if (inc_pools_and_fb) {flooded_basement .draw(s, anim_state, anim_time);}
		anim_state.clear_animation_id(s); // clear animations
		check_mvm_update();
		glEnable(GL_CULL_FACE);
	}
}; // fish_manager_t

fish_manager_t fish_manager;

void register_fishtank(room_object_t const &obj, bool is_visible) {fish_manager.register_fishtank(obj, is_visible);}

void end_fish_draw(shader_t &s, bool inc_pools_and_fb) {
	fish_manager.end_objs();
	fish_manager.draw_fish(s, inc_pools_and_fb);
}


bool building_t::begin_fish_draw() const { // returns true of pool or basement water is visible
	assert(interior);
	fish_manager.next_frame(*this);
	if (!water_visible_to_player()) return 0;
	// handle underwater fish in pools and basements
	int const rseed(interior->rooms.size());
	vect_cube_t obstacles;
	
	if (has_pool()) {
		if (!fish_manager.has_swimming_pool()) { // not yet setup
			float const tile_thickness(get_flooring_thick());
			cube_t water_bcube(interior->pool);
			water_bcube.z1() += tile_thickness; // subtract off tile width on bottom
			water_bcube.z2()  = interior->water_zval; // water level is below the top of the pool
			water_bcube.expand_by_xy(-tile_thickness); // subtract off tile width on sides
			room_object_t pool_ramp;

			if (has_room_geom()) { // should always be true
				vect_room_object_t const &objs(interior->room_geom->objs);
				unsigned const ramp_ix(interior->room_geom->pool_ramp_obj_ix), upper_ix(ramp_ix + 1), ladder_ix(ramp_ix + 3); // {ramp, upper, diving board, ladder}
				unsigned const stairs_ix(interior->room_geom->pool_stairs_start_ix);
			
				if (ramp_ix > 0) { // add ramp, upper part, and pool ladder
					assert(ramp_ix < objs.size());
					if (objs[ramp_ix].type == TYPE_POOL_TILE) {pool_ramp = objs[ramp_ix];}
					if (upper_ix  < objs.size() && objs[upper_ix ].type == TYPE_POOL_TILE) {obstacles.push_back(objs[upper_ix ]);}
					if (ladder_ix < objs.size() && objs[ladder_ix].type == TYPE_POOL_LAD ) {obstacles.push_back(objs[ladder_ix]);}
				}
				if (stairs_ix > 0) { // has stairs; should always be true
					assert(stairs_ix < objs.size());

					for (auto i = objs.begin()+stairs_ix; i != objs.end(); ++i) { // add stairs and railings as obstacles
						if (i->type != TYPE_STAIR && i->type != TYPE_RAILING) break; // done
						obstacles.push_back(*i);
					}
				}
			}
			fish_manager.register_swimming_pool(interior->pool, pool_ramp, water_bcube, obstacles, rseed);
		}
	}
	else if (has_ext_basement()) { // flooded backrooms basement
		cube_t water_bcube(get_water_cube());
		water_bcube.z1() += get_fc_thickness();
		min_eq(water_bcube.z2(), water_bcube.z1() + get_floor_ceil_gap()); // constraint to the lowest level
		if (water_bcube.z2() <= water_bcube.z1()) return 0; // shouldn't happen?
		if (get_camera_pos().z > water_bcube.z2() + get_window_vspace()) return 0; // player on the floor above the water, fish likely not visible
		water_bcube.expand_by_xy(-get_wall_thickness()); // subtract off exterior walls

		if (!fish_manager.has_flooded_basement()) { // not yet setup
			for (stairwell_t const &s : interior->stairwells) { // add stairs
				if (water_bcube.intersects(s)) {obstacles.push_back(s);}
			}
			if (has_room_geom()) { // add room objects
				vect_room_object_t const &objs(interior->room_geom->objs);
				unsigned const objs_start(interior->room_geom->backrooms_start);

				if (objs_start > 0) { // add walls, pillars, and other placed room objects
					assert(objs_start < objs.size());
				
					for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) { // add obstacles
						if (i->type == TYPE_BLOCKER || i->type == TYPE_FLOORING || i->type == TYPE_WALL_TRIM || i->type == TYPE_OUTLET || i->type == TYPE_VENT) continue; // ignore
						if (water_bcube.intersects(*i)) {obstacles.push_back(*i);}
					}
				}
			}
			fish_manager.register_flooded_basement(water_bcube, get_window_vspace(), obstacles, rseed);
		}
		fish_manager.register_doors_state(interior->door_stacks, interior->doors);
	}
	return 1;
}

