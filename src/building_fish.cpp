// 3D World - Building Fish
// by Frank Gennari 6/23/24

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h" // for animation_state_t

extern int animate2, frame_counter;
extern float fticks;

bool check_ramp_collision(room_object_t const &c, point &pos, float radius, vector3d *cnorm);
bool add_water_splash(point const &pos, float radius, float size);
void draw_animated_fish_model(shader_t &s, vector3d const &pos, float radius, vector3d const &dir, float anim_time, colorRGBA const &color);


class fish_manager_t {
	struct fish_t {
		float radius=0.0, mspeed=0.0, tspeed=0.0; // movement speed and tail speed
		unsigned id=0;
		int next_splash_frame=0;
		point pos; // p_last?
		vector3d dir, target_dir; // Note: fish only rotated about Z and remains in the XY plane
		colorRGBA color=WHITE;

		void draw(shader_t &s, animation_state_t &anim_state, float anim_time) const {
			anim_state.anim_time = 0.5*tspeed*anim_time; // tail moves twice per second
			if (anim_state.enabled) {anim_state.set_animation_id_and_time(s, 0, 1.0);}
			draw_animated_fish_model(s, pos, radius, dir, anim_state.anim_time, color);
		}
		bool can_splash(rand_gen_t &rgen) {
			if (frame_counter < next_splash_frame) return 0;
			next_splash_frame = frame_counter + (rgen.rand() % (2*TICKS_PER_SECOND)); // once every 2s on average
			return 1;
		}
	};

	class fish_cont_t {
	protected:
		cube_t bcube, valid_area;
		rand_gen_t rgen;
		vector<fish_t> fish;
	public:
		~fish_cont_t() {}
		bool empty() const {return fish.empty();}
		void clear() {fish.clear(); bcube.set_to_zeros();}

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
				for (unsigned d = 0; d < 3; ++d) {f.color[d] *= rgen.rand_uniform(0.8, 1.0);} // slight color variation
				cube_t const valid_area(get_valid_area(1.05*f.radius)); // slightly larger radius so that we don't start out intersecting

				for (unsigned n = 0; n < 100; ++n) { // 100 tries
					for (unsigned d = 0; d < 3; ++d) {f.pos[d] = rgen.rand_uniform(valid_area.d[d][0], valid_area.d[d][1]);}
					point coll_pos; // unused
					if (!check_fish_coll(f.pos, f.radius, f.id, coll_pos)) break; // success
				}
			} // for f
		}
		void next_frame() {
			if (!animate2) return;
			float const delta_dir(0.6*(1.0 - pow(0.7f, fticks)));

			for (fish_t &f : fish) {
				if (f.target_dir != zero_vector) { // fish is turning in place
					f.dir = delta_dir*f.target_dir + (1.0 - delta_dir)*f.dir;
					if (f.dir == zero_vector) {f.dir = assign_fish_dir();} // error?
					else {f.dir.normalize();}
					if (dot_product(f.dir, f.target_dir) > 0.95) {f.target_dir = zero_vector;} // turn complete
					continue;
				}
				cube_t const valid_area(get_valid_area(f.radius));
				point const prev_pos(f.pos);
				assert(valid_area.contains_pt(prev_pos));
				f.pos += f.mspeed*f.dir*fticks;
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
				}
			} // for f
		}
		void draw(shader_t &s, animation_state_t &anim_state, float anim_time) const {
			for (fish_t const &f : fish) {f.draw(s, anim_state, anim_time);}
		}
	protected:
		vector3d assign_fish_dir() {
			return (rgen.signed_rand_vector_spherical() * bcube.get_size() * vector3d(1.0, 1.0, 0.5)).get_norm(); // match the aspect ratio of the container, smaller in Z
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
			} // for f
			return 0;
		}
	}; // fish_cont_t

	class fishtank_t : public fish_cont_t {
	public:
		unsigned obj_id=0; // used as a unique identifier
		bool present=0, visible=0;

		fishtank_t(room_object_t const &obj) : obj_id(obj.obj_id) {
			cube_t tank_inner(obj);
			tank_inner.expand_by(-vector3d(0.04, 0.04, 0.1)*obj.get_height()); // subtract off the glass and some area from the top and bottom
			init(obj, tank_inner, obj_id, 1, 4); // 1-4 fish
			float const max_fish_radius(0.125*(1.0 + 1.0/fish.size())*bcube.min_len()); // more fish = smaller size
			populate(max_fish_radius, 0.0016);
		}
		void next_frame() {
			present = visible = 0; // mark as not present or visible until it's seen
			fish_cont_t::next_frame();
		}
		void update_object(room_object_t const &obj) { // handle movement when the table is pushed
			if (obj == bcube) return; // no update
			vector3d const delta(obj.get_llc() - bcube.get_llc());
			bcube = obj;
			for (fish_t &f : fish) {f.pos += delta;}
		}
		void draw(shader_t &s, animation_state_t &anim_state, float anim_time) const {
			if (visible) {fish_cont_t::draw(s, anim_state, anim_time);}
		}
	}; // fishtank_t

	class area_fish_cont_t : public fish_cont_t {
	protected:
		vect_cube_t obstacles;
	public:
		void next_frame() {
			// TODO: handle player interaction
			fish_cont_t::next_frame();

			for (fish_t &f : fish) { // handle water splashes
				if (f.pos.z < valid_area.z2() - 2.0*f.radius) continue; // too low to splash
				if (f.can_splash(rgen)) {add_water_splash(f.pos, 2.0*f.radius, 0.25);}
			}
		}
		virtual bool check_fish_coll(point const &pos, float radius, unsigned id, point &coll_pos) const {
			if (fish_cont_t::check_fish_coll(pos, radius, id, coll_pos)) return 1;

			for (cube_t const &c : obstacles) {
				if (sphere_cube_intersect(pos, 1.5*radius, c)) {coll_pos = c.get_cube_center(); return 1;} // extra radius for clearance
			}
			return 0;
		}
	}; // player_int_fish_cont_t

	class swimming_pool_t : public area_fish_cont_t {
		//indoor_pool_t const &pool;
		room_object_t pool_ramp;
	public:
		void init(indoor_pool_t const &pool, room_object_t const &pool_ramp_, cube_t const &water_bcube, vect_cube_t const &obstacles_, int rseed) {
			obstacles = obstacles_;
			pool_ramp = pool_ramp_;
			fish_cont_t::init(pool, water_bcube, rseed, 10, 40); // 1-4 fish FIXME
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
	public:
		void init(cube_t const &water_bc, float floor_spacing, vect_cube_t const &obstacles_, int rseed) {
			obstacles = obstacles_;
			fish_cont_t::init(water_bc, water_bc, rseed, 10, 20); // 10-20 fish
			float const max_fish_radius(0.25*min(bcube.dz(), 0.25f*floor_spacing));
			populate(max_fish_radius, 0.004);
		}
	}; // flooded_basement_t

	vector<fishtank_t> fishtanks;
	swimming_pool_t swimming_pool;
	flooded_basement_t flooded_basement;
	building_t const *prev_building=nullptr;
	unsigned fishtank_ix=0;
	float anim_time=0.0;
public:
	bool empty() const {return (fishtanks.empty() && swimming_pool.empty() && flooded_basement.empty());}

	void clear() {
		fishtanks       .clear();
		swimming_pool   .clear();
		flooded_basement.clear();
	}
	void next_frame(building_t const &building) {
		if (&building != prev_building) { // new building
			clear();
			prev_building = &building;
			anim_time     = 0.0; // reset
		}
		for (fishtank_t &ft : fishtanks) {ft.next_frame();}
		swimming_pool   .next_frame();
		flooded_basement.next_frame();
		if (animate2) {anim_time += fticks;}
		fishtank_ix = 0;
	}
	void register_fishtank(room_object_t const &obj, bool is_visible) {
		if (fishtanks.size() <= fishtank_ix) {fishtanks.push_back(fishtank_t(obj));} // fishtank not yet added

		for (fishtank_t &ft : fishtanks) {
			if (ft.obj_id == obj.obj_id) {
				assert(!ft.present); // check for duplicate obj_id
				ft.update_object(obj);
				ft.present = 1;
				ft.visible = is_visible;
			}
		}
		++fishtank_ix;
	}
	void register_swimming_pool(indoor_pool_t const &pool, room_object_t const &pool_ramp, cube_t const &water_bcube, vect_cube_t const &obstacles, int rseed) {
		if (swimming_pool.empty()) {swimming_pool.init(pool, pool_ramp, water_bcube, obstacles, rseed);}
	}
	void register_flooded_basement(cube_t const &water_bc, float floor_spacing, vect_cube_t const &obstacles, int rseed) {
		if (flooded_basement.empty()) {flooded_basement.init(water_bc, floor_spacing, obstacles, rseed);}
	}
	void end_objs() { // remove any fishtanks that are no longer present
		fishtanks.erase(remove_if(fishtanks.begin(), fishtanks.end(), [](fishtank_t const &ft) {return !ft.present;}), fishtanks.end());
		assert(fishtank_ix == fishtanks.size());
	}
	void draw_fish(shader_t &s, bool inc_pools_and_fb) const {
		if (empty() || (!inc_pools_and_fb && fishtanks.empty())) return;
		glDisable(GL_CULL_FACE); // fix for looking through the fish's mouth
		animation_state_t anim_state(1, ANIM_ID_FISH_TAIL); // enabled=1; animation_scale and model_delta_height are not set or used
		for (fishtank_t const &f : fishtanks) {f.draw(s, anim_state, anim_time);}
		if (inc_pools_and_fb) {swimming_pool    .draw(s, anim_state, anim_time);}
		if (inc_pools_and_fb) {flooded_basement .draw(s, anim_state, anim_time);}
		anim_state.clear_animation_id(s); // clear animations
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
	bool inc_pools_and_fb(0);
	fish_manager.next_frame(*this);
	if (!water_visible_to_player()) return 0;
	// handle underwater fish in pools and basements
	int const rseed(interior->rooms.size());
	vect_cube_t obstacles;
	
	if (has_pool()) {
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
	else if (has_ext_basement()) { // flooded backrooms basement
		// TODO: backrooms walls, doors, stairs, larger objects, etc.
		cube_t water_bcube(get_water_cube());
		water_bcube.z1() += get_fc_thickness();
		min_eq(water_bcube.z2(), water_bcube.z1() + get_floor_ceil_gap()); // constraint to the lowest level
		if (water_bcube.z2() <= water_bcube.z1()) return 0; // shouldn't happen?
		fish_manager.register_flooded_basement(water_bcube, get_window_vspace(), obstacles, rseed);
	}
	return 1;
}

