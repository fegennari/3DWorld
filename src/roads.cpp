// 3D World - Roads for Procedural Cities
// by Frank Gennari
// 11/20/18
#include "city.h"
#include "lightmap.h"

float const STREETLIGHT_BEAMWIDTH       = 0.225;
float const SLIGHT_DIST_TO_CORNER_SCALE = 3.0; // larger is closer to the road surface
float const BRIDGE_HEIGHT_TO_LEN        = 0.3;

extern bool tt_fire_button_down;
extern int frame_counter, game_mode, display_mode;
extern float fticks, FAR_CLIP;
extern vector<light_source> dl_sources;
extern city_params_t city_params;


string gen_random_first_name(rand_gen_t &rgen); // from pedestrians.cpp
void add_cylin_as_tris(vector<vert_norm_tc_color> &verts, point const ce[2], float r1, float r2, color_wrapper const &cw,
	unsigned ndiv, unsigned draw_top_bot, float tst=1.0, float tss=1.0, bool swap_ts_tt=0);
bool line_quad_intersect(point const &p1, point const &p2, point const *const pts, float &t);

class road_name_gen_t {
	static string get_numbered_street_name(unsigned num) {
		assert(num > 0 && num < 100);
		string const names_1_to_20[21] = {"Zeroth", "First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth", "Ninth", "Tenth",
			"Eleventh", "Twelfth", "Thirteenth", "Fourteenth", "Fifteenth", "Sixteenth", "Seventeenth", "Eighteenth", "Nineteenth", "Twentieth"};
		string const decade_names[10] = {"", "Ten", "Twent", "Thirt", "Fort", "Fift", "Sixt", "Sevent", "Eight", "Ninet"}; // missing the 'y'
		if (num <= 20) {return names_1_to_20[num];}
		if ((num%10) == 0) {return decade_names[num/10] + "ieth";}
		return decade_names[num/10] + "y " + names_1_to_20[num%10];
	}
	static string get_dir_name(bool dim, bool dir) {
		return (dim ? (dir ? "North" : "South") : (dir ? "East" : "West"));
	}
	static string select_road_name(rand_gen_t &rgen) { // use either a randomly generated name or a person's first name
		if (rgen.rand_bool()) {return gen_random_name(rgen, 4);}
		string const name(gen_random_first_name(rgen));
		return ((name.size() < 4) ? gen_random_name(rgen, 4) : name); // 2-3 letter names look bad on signs (too stretched)
	}
	string gen_name_int(road_t const &road, unsigned city_ix) const {
		if (road.is_all_zeros()) { // city connector road; only city_ix is used
			rand_gen_t rgen;
			rgen.set_state(city_ix+123, city_ix+321); // should be unique per road
			rgen.rand_mix();
			string const suffix[3] = {"St", "Rd", "Expy"};
			return select_road_name(rgen) + " " + suffix[rgen.rand()%3];
		}
		if (road.dim == 1 && road.road_ix < 99 && (city_ix & 1)) { // can use numbered roads, for odd cities
			string const suffix[2] = {"St", "Ave"};
			return get_numbered_street_name(road.road_ix+1) + " " + suffix[(city_ix & 2) >> 1]; // alternate suffixes between cities (may be unused)
		}
		rand_gen_t rgen;
		rgen.set_state(road.road_ix+1, city_ix+1); // should be unique per road
		rgen.rand_mix();
		string prefix;

		if ((rgen.rand() & 3) == 0) { // consider adding a direction 25% of the time; would make more sense if this name applies to half a road
			bool const dim(!road.dim);
			cube_t const city_bcube(get_city_bcube(city_ix));
			float const t_not_dim((road.get_center_dim(dim) - city_bcube.d[dim][0])/city_bcube.get_sz_dim(dim));
			if (t_not_dim > 0.7 || t_not_dim < 0.3) {prefix = get_dir_name(dim, (t_not_dim > 0.5)) + " ";} // use dir name if far from center
		}
		string const suffix[13] = {"St", "St", "St", "Ave", "Ave", "Rd", "Rd", "Rd", "Dr", "Blvd", "Ln", "Way", "Ct"}; // more common suffixes are duplicated
		return prefix + select_road_name(rgen) + " " + suffix[rgen.rand()%13];
	}
	map<road_t, string> road_name_cache; // only need to check road x/y
public:
	string gen_name(road_t const &road, unsigned city_ix) {
		string &name(road_name_cache[road]);
		if (name.empty()) {name = gen_name_int(road, city_ix);} // generate if needed
		return name;
	}
};
road_name_gen_t road_name_gen;
string road_t::get_name(unsigned city_ix) const {return road_name_gen.gen_name(*this, city_ix);}

void road_mat_mgr_t::ensure_road_textures() {
	if (inited) return;
	string const img_names[NUM_RD_TIDS] = {"sidewalk.jpg", "straight_road.jpg", "bend_90.jpg", "int_3_way.jpg", "int_4_way.jpg",
		                                   "parking_lot.png", "rail_tracks.jpg", "grass_park.jpg", "concrete.jpg", "asphalt.jpg"};
	float const aniso   [NUM_RD_TIDS] = {4.0, 16.0, 8.0, 8.0, 8.0, 4.0, 16.0, 16.0};
	int   const wrap_mir[NUM_RD_TIDS] = {1, 1, 0, 0, 0, 1, 1, 1, 1, 1}; // bend and intersections are clamped, and the others are wrapped
	for (unsigned i = 0; i < NUM_RD_TIDS; ++i) {tids[i] = get_texture_by_name(("roads/" + img_names[i]), 0, 0, wrap_mir[i], aniso[i]);}
	sl_tid = get_texture_by_name("roads/traffic_light.png");
	inited = 1;
}

void road_mat_mgr_t::set_texture(unsigned type) {
	assert(type < NUM_RD_TIDS);
	ensure_road_textures();
	select_texture(tids[type]);
}
void road_mat_mgr_t::set_stoplight_texture() {
	ensure_road_textures();
	select_texture(sl_tid);
}

tex_range_t get_uniform_tscale_ar(float length, float width, float base_tscale, bool dim) { // for driveways and road skirts
	float txy[2] = {base_tscale, base_tscale};
	txy[dim] *= length/width; // ensure 1:1 texture scale
	return tex_range_t(0.0, 0.0, txy[0], txy[1], 0, 0); // since the asphalt driveway texture isn't directional, we don't need to worry about swap_xy
}

road_t::road_t(point const &s, point const &e, float width, bool dim_, bool slope_, unsigned road_ix_) : road_ix(road_ix_), dim(dim_), slope(slope_) {
	assert(s != e);
	assert(width > 0.0);
	vector3d const dw(0.5*width*cross_product((e - s), plus_z).get_norm());
	point const pts[4] = {(s - dw), (s + dw), (e + dw), (e - dw)};
	set_from_points(pts, 4);
}

void road_t::register_bridge_or_tunnel(cube_t const &bc, bool is_bridge) {
	(is_bridge ? has_bridge : has_tunnel) = 1;
	bt_lo = bc.d[dim][0];
	bt_hi = bc.d[dim][1];
}

// specialized here for sloped roads (road segments and railroad tracks)
void road_t::add_road_quad(quad_batch_draw &qbd, colorRGBA const &color, float ar, bool add_skirt) const {
	if (z1() == z2()) {
		if (!add_skirt) {add_flat_city_quad(*this, qbd, color, ar);} // no skirt for flat roads
		return;
	}
	if (add_skirt && has_bridge) {
		for (unsigned d = 0; d < 2; ++d) { // split into two sections around the bridge
			road_t sub_road(*this);
			sub_road.has_bridge = 0;
			sub_road.d[dim][!d] = (d ? (bt_hi + 0.25*city_params.road_width) : (bt_lo - 0.25*city_params.road_width));
			float const len(get_length()), sub_len(sub_road.get_length()), delta_z(dz()*(len - sub_len)/len);
			if (sub_len <= 0.0) continue; // no segment
			if (slope ^ bool(d)) {sub_road.z1() += delta_z;} else {sub_road.z2() -= delta_z;} // adjust zvals
			sub_road.add_road_quad(qbd, color, ar, add_skirt);
		} // for d
		return;
	}
	bool const s(slope ^ dim);
	point pts[4] = {point(x1(), y1(), d[2][!s]), point(x2(), y1(), d[2][!s]), point(x2(), y2(), d[2][s]), point(x1(), y2(), d[2][s])};
	if (!dim) {swap(pts[0].z, pts[2].z);}

	if (add_skirt) {
		float const road_length(get_length()), skirt_width(0.2*get_width()), skirt_height(1.0*skirt_width);

		for (unsigned d = 0; d < 2; ++d) { // each side
			point skirt[4];
			if (dim) {skirt[0] = skirt[3] = pts[d ? 1 : 0]; skirt[1] = skirt[2] = pts[d ? 2 : 3];} // y
			else     {skirt[0] = skirt[3] = pts[d ? 2 : 0]; skirt[1] = skirt[2] = pts[d ? 3 : 1];} // x

			for (unsigned e = 0; e < 2; ++e) {
				skirt[e+2][!dim] += (d ? 1.0 : -1.0)*skirt_width; // expand outward
				skirt[e+2].z     -= skirt_height; // move downward
			}
			vector3d normal(cross_product((skirt[2] - skirt[1]), (skirt[0] - skirt[1])).get_norm());
			if (normal.z < 0.0) {reverse(skirt, skirt+4); normal.negate();} // make CCW
			qbd.add_quad_pts(skirt, color, normal, get_uniform_tscale_ar(road_length, skirt_width, 1.0, 0));
		} // for d
	}
	else { // add road quad
		vector3d const normal(cross_product((pts[2] - pts[1]), (pts[0] - pts[1])).get_norm());
		qbd.add_quad_pts(pts, color, normal, get_tex_range(ar));
	}
}


tex_range_t parking_lot_t::get_tex_range(float ar) const { // ar is unused
	bool const d(!dim); // Note: R90
	float const xscale(1.0/(2.0*PARK_SPACE_WIDTH *city_params.get_nom_car_size().y));
	float const yscale(1.0/(1.0*PARK_SPACE_LENGTH*city_params.get_nom_car_size().x)*((dim^dir) ? -1.0 : 1.0)); // swap texture in Y for some dim/dir
	float const tx(0.24), ty(0.0); // x=cols, y=rows
	return tex_range_t(tx, ty, (xscale*(d ? dy() : dx()) + tx), (yscale*(d ? dx() : dy()) + ty), 0, d);
}

void driveway_t::mark_ped_this_frame() const {last_ped_frame = frame_counter;} // can/must be const; last_ped_frame is mutable
bool driveway_t::has_recent_ped() const {return (frame_counter <= (int)last_ped_frame+1);} // allow one frame lag so that it doesn't matter which thread updates first

tex_range_t driveway_t::get_tex_range(float ar) const { // ar is unused
	return get_uniform_tscale_ar(get_length(), get_width(), 2.0, dim); // base_tscale=2.0
}
cube_t driveway_t::extend_across_road() const {
	cube_t dw_ext(*this); // includes driveway and the road adjacent to it
	dw_ext.d[dim][dir] += (dir ? 1.0 : -1.0)*city_params.road_width;
	return dw_ext;
}

namespace stoplight_ns {

	rand_gen_t stoplight_rgen;

	float stoplight_max_height() {return 9.2*(0.03*city_params.road_width);} // 2.0*h + 1.2*h*num_segs
	
	void stoplight_t::advance_state() {
		if (num_conn == 4) {next_state();} // all states are valid for 4-way intersections
		else { // 3-way intersection
			assert(num_conn == 3);
			while (1) {
				next_state();
				bool valid(0);
				switch (conn) {
				case 7 : {bool const allow[6] = {0,1,1,1,0,0}; valid = allow[cur_state]; break;} // no +y / S
				case 11: {bool const allow[6] = {1,1,0,0,0,1}; valid = allow[cur_state]; break;} // no -y / N
				case 13: {bool const allow[6] = {1,0,0,1,1,0}; valid = allow[cur_state]; break;} // no +x / W
				case 14: {bool const allow[6] = {0,0,1,0,1,1}; valid = allow[cur_state]; break;} // no -x / E
				default: assert(0);
				}
				if (valid) break;
			} // end while
		}
		cur_state_ticks = 0.0; // reset for this state
	}

	bool stoplight_t::is_any_car_waiting_at_this_state() const {
		if (num_conn == 2) return 0; // 2-way intersection, no cross traffic
		return ((left_orient_masks[cur_state] & car_waiting_left) || (st_r_orient_masks[cur_state] & car_waiting_sr));
	}

	void stoplight_t::find_state_with_waiting_car() {
		uint8_t const prev_state(cur_state);

		while (1) {
			// cycle through all directions and find one where a car is waiting to go; this allows the traffic light to skip some green lights that aren't needed;
			// TODO: however, it can lead to a light transitioning from green => yellow => green without going through red; this is difficult to fix with the current architecture
			advance_state();
			if (any_blocked()) break; // if some car is blocking the intersection in some dir, force all states (no skipped lights) to guarantee we can make progress
			if (is_any_car_waiting_at_this_state()) break; // car is waiting at this state
			if (cur_state == prev_state) {advance_state(); break;} // wrapped around, leave at the next valid state after the prev state
		}
		car_waiting_sr = car_waiting_left = 0; // waiting bits have been used, clear for next state change
	}

	void stoplight_t::run_update_logic() {
		assert(cur_state < NUM_STATE);
		//if (cur_state_ticks > get_cur_state_time_secs()) {advance_state();} // time to update to next state
		if (cur_state_ticks > get_cur_state_time_secs()) {find_state_with_waiting_car();} // time to update to next state
	}

	void stoplight_t::ffwd_to_future(float time_secs) {
		cur_state_ticks += TICKS_PER_SECOND*time_secs;
		run_update_logic();
	}

	void stoplight_t::init(uint8_t num_conn_, uint8_t conn_) {
		num_conn = num_conn_; conn = conn_;
		if (num_conn == 2) return; // nothing else to do
		cur_state = stoplight_rgen.rand() % NUM_STATE; // start at a random state
		advance_state(); // make sure cur_state is valid
		cur_state_ticks = get_cur_state_time_secs()*stoplight_rgen.rand_float(); // start at a random time within this state
	}

	void stoplight_t::next_frame() {
		reset_blocked();
		cw_in_use = 0;
		if (num_conn == 2) return; // nothing else to do
		cur_state_ticks += fticks;
		run_update_logic();
	}

	void stoplight_t::notify_waiting_car(bool dim, bool dir, unsigned turn) const {
		unsigned const orient(2*dim + dir); // {W, E, S, N}
		((turn == TURN_LEFT) ? car_waiting_left : car_waiting_sr) |= (1 << orient);
	}

	bool stoplight_t::red_light(bool dim, bool dir, unsigned turn) const {
		assert(cur_state < NUM_STATE);
		assert(turn == TURN_LEFT || turn == TURN_RIGHT || turn == TURN_NONE);
		if (num_conn == 2) return 0; // 2-way intersection, no cross traffic
		unsigned const mask(((turn == TURN_LEFT) ? left_orient_masks : st_r_orient_masks)[cur_state]);
		unsigned const orient(2*dim + dir); // {W, E, S, N}
		return (((1<<orient) & mask) == 0);
	}

	unsigned stoplight_t::get_light_state(bool dim, bool dir, unsigned turn) const { // 0=green, 1=yellow, 2=red
		if (red_light(dim, dir, turn)) return RED_LIGHT;
		if (num_conn == 2) return GREEN_LIGHT;
		stoplight_t future_self(*this);
		future_self.ffwd_to_future(2.0); // yellow light time = 2.0s
		return (future_self.red_light(dim, dir, turn) ? (unsigned)YELLOW_LIGHT : (unsigned)GREEN_LIGHT);
	}
	unsigned stoplight_t::get_future_light_state(bool dim, bool dir, unsigned turn, float future_seconds) const {
		if (future_seconds == 0.0) {return get_light_state(dim, dir, turn);}
		assert(future_seconds > 0.0); // not in the past
		stoplight_t future_self(*this);
		future_self.ffwd_to_future(future_seconds);
		return future_self.get_light_state(dim, dir, turn);
	}

	bool stoplight_t::can_walk(bool dim, bool dir) const { // Note: symmetric across the two sides of the street
		// Note: ignores blocked state and right-on-red
		unsigned const orient(2*dim + dir); // {W, E, S, N}
		assert(conn & (1<<orient));
		// check for cars entering the isec from this side
		if ((1<<orient) & (left_orient_masks[cur_state] | st_r_orient_masks[cur_state])) return 0; // any turn dir
		// check for cars entering the isec from the other side; Note that we don't check conn mask here because unconnected dirs don't get green lights
		if (st_r_orient_masks[cur_state] & (1<<other_lane[orient])) return 0; // opposing traffic going straight
		if (left_orient_masks[cur_state] & (1<<to_right  [orient])) return 0; // traffic to the right turning left
		//if (st_r_orient_masks[cur_state] & (1<<to_left   [orient])) return 0; // traffic to the left turning right (skip)
		return 1;
	}

	unsigned stoplight_t::get_crosswalk_state(bool dim, bool dir) const {
		if (!can_walk(dim, dir)) return CW_STOP;
		stoplight_t future_self(*this);
		future_self.ffwd_to_future(2.0); // crosswalk warn time = 2.0s
		if (!future_self.can_walk(dim, dir)) return CW_WARN;
		return CW_WALK;
	}

	bool stoplight_t::check_int_clear(unsigned orient, unsigned turn_dir) const { // check for cars on other lanes blocking the intersection
		switch (turn_dir) { // orient: {W, E, S, N}
		case TURN_NONE:  return (!blocked[to_right[orient]] && !blocked[to_left[orient]]); // straight
		case TURN_LEFT:  return (!blocked[to_right[orient]] && !blocked[to_left[orient]] && !blocked[other_lane[orient]]);
		case TURN_RIGHT: return (!blocked[to_right[orient]]);
		}
		return 1;
	}

	bool stoplight_t::can_turn_right_on_red(car_base_t const &car) const { // check for legal right on red (no other lanes turning into the road to our right)
		if (car.turn_dir != TURN_RIGHT) return 0;
		unsigned const orient(car.get_orient());  // {W, E, S, N}
		if (!red_light(!car.dim, (to_right  [orient] & 1), TURN_NONE)) return 0; // traffic to our left has a green or yellow light for going straight
		if (!red_light( car.dim, (other_lane[orient] & 1), TURN_LEFT)) return 0; // opposing traffic has a green or yellow light for turning left
		// Note: there are no U-turns, so we don't need to worry about cars coming from the other dir of our dest road
		if (cw_in_use & (1 << other_lane[orient])) return 0; // ped in the crosswalk, wait
		return 1; // can turn right on red to_left[orient]
	}

	string stoplight_t::str() const {
		std::ostringstream oss;
		oss << TXTi(num_conn) << TXTi(conn) << TXTi(cur_state) << TXT(at_conn_road) << TXT(cur_state_ticks) << TXTi(car_waiting_sr) << TXTi(car_waiting_left)
			<< "blocked: " << blocked[0] << blocked[1] << blocked[2] << blocked[3];
		return oss.str();
	}

	string stoplight_t::label_str() const {
		std::ostringstream oss;
		oss << TXTi(num_conn) << TXTin(conn) << TXTi(cur_state) << TXTn(cur_state_ticks) << TXTn(at_conn_road) << TXTi(car_waiting_sr) << TXTin(car_waiting_left)
			<< "blocked: " << blocked[0] << blocked[1] << blocked[2] << blocked[3];
		return oss.str();
	}

} // end stoplight_ns


namespace streetlight_ns {

	float get_streetlight_height     () {return light_height*city_params.road_width;}
	float get_streetlight_pole_radius() {return pole_radius *city_params.road_width;}

	point streetlight_t::get_lpos() const {
		float const height(get_streetlight_height());
		return (pos + vector3d(0.0, 0.0, 1.1*height) + 0.4*height*dir);
	}

	void streetlight_t::draw(road_draw_state_t &dstate, bool shadow_only, bool is_local_shadow, bool always_on) const {
		// Note: translate has already been applied as a transform
		float const height(get_streetlight_height());
		point const center(pos + dstate.xlate + vector3d(0.0, 0.0, 0.5*height));
		if (shadow_only && is_local_shadow && !dist_less_than(camera_pdu.pos, center, 0.8*camera_pdu.far_)) return;
		if (!camera_pdu.sphere_visible_test(center, height)) return; // VFC
		bool const is_on(!shadow_only && is_lit(always_on));
		point const lpos(get_lpos());
		float camera_dist(0.0), dist_val(0.0);

		if (shadow_only) {dist_val = (is_local_shadow ? 0.12 : 0.06);}
		else {
			// draw an emissive circular blur texture area on the ground below the streetlight if it's on but outside our city lights bcube
			if (is_on && !get_city_lights_bcube().contains_pt_xy(lpos)) {
				point below_pos(lpos);
				below_pos.z = pos.z + 0.01*height; // shift up slightly to prevent Z-fighting
				point pts[4];
				set_z_plane_square_pts(below_pos, 1.2*city_params.road_width, pts);
				dstate.qbd_emissive.add_quad_pts(pts, colorRGBA(light_color, 0.25), plus_z);
			}
			camera_dist = p2p_dist(camera_pdu.pos, center);
			dist_val    = camera_dist/dstate.draw_tile_dist;
			if (dist_val > 0.2) return; // too far
		}
		float const pradius(get_streetlight_pole_radius()), lradius(light_radius*city_params.road_width);
		int const ndiv(shadow_only ? (is_local_shadow ? 4 : 8) : max(4, min(N_SPHERE_DIV, int(0.5/dist_val))));
		point const top(pos + vector3d(0.0, 0.0, 0.96*height)), arm_end(lpos + vector3d(0.0, 0.0, 0.025*height) - 0.06*height*dir);
		point const vpost_ce[2] = {pos, (pos + height*plus_z)};
		add_cylin_as_tris(dstate.qbd_untextured.verts, vpost_ce, pradius, 0.7*pradius, pole_color, min(ndiv, 24), 0); // untextured, no ends

		if (dist_val <= 0.12) {
			point const hbar_ce[2] = {top, arm_end};
			add_cylin_as_tris(dstate.qbd_untextured.verts, hbar_ce, 0.5*pradius, 0.4*pradius, pole_color, min(ndiv, 16), 0); // untextured, no ends
		}
		if (shadow_only && is_local_shadow) return; // top part never projects a shadow on a visible object near the ground

		if (!shadow_only) {
			if (!is_on && dist_val > 0.15) return; // too far
			if (is_on) {dstate.s.set_color_e(light_color);} else {dstate.s.set_cur_color(light_color);} // emissive when lit
			
			if (is_on && dist_val > 0.01 && (camera_pdu.pos.z - center.z) < 0.5*camera_dist) { // streetlight bloom, drawn if camera is not too high above
				// should be elliptical, and needs tighter alpha mask
				dstate.add_light_flare((lpos - 0.6*lradius*plus_z), zero_vector, light_color, min(1.0f, 12.0f*dist_val), 4.0*lradius); // non-directional
			}
		}
		if (!shadow_only && dist_val < 0.05) {try_bind_tile_smap_at_point(center, dstate.s);} // bind the correct shadow map when near the player
		fgPushMatrix();
		translate_to(lpos);
		scale_by(lradius*vector3d(1.0+fabs(dir.x), 1.0+fabs(dir.y), 1.0)); // scale 2x in dir
		bind_draw_sphere_vbo(0, 1); // untextured
		draw_sphere_vbo_pre_bound(ndiv, 0); // untextured
		if (!shadow_only && is_on) {dstate.s.clear_color_e();}

		if (!shadow_only && dist_val < 0.12) {
			fgTranslate(0.0, 0.0, 0.1); // translate up slightly and draw top cap of light
			dstate.s.set_cur_color(pole_color);
			draw_sphere_vbo_pre_bound(ndiv, 0); // untextured
		}
		bind_vbo(0);
		fgPopMatrix();
	}

	void streetlight_t::add_dlight(vector3d const &xlate, cube_t &lights_bcube, bool always_on) const {
		if (!is_lit(always_on)) return;
		float const ldist(light_dist*city_params.road_width);
		if (!lights_bcube.contains_pt_xy(pos)) return; // not contained within the light volume
		point const lpos(get_lpos());
		if (!camera_pdu.sphere_visible_test((lpos + xlate), ldist)) return; // VFC
		min_eq(lights_bcube.z1(), (lpos.z - ldist));
		max_eq(lights_bcube.z2(), (lpos.z + 0.1f*ldist)); // pointed down - don't extend as far up
		//float const far_clip(1.2*get_streetlight_height()); // slightly beyond the height above the ground - faster, but misses shadows due to incorrect far plane test
		float const far_clip(0.0); // default radius far clip
		dl_sources.emplace_back(ldist, lpos, lpos, light_color, 0, -plus_z, STREETLIGHT_BEAMWIDTH, 0.0, 0, 0.0, far_clip); // points down
		
		// cache shadow maps if there are no dynamic cars or pedestrians (player doesn't cast a shadow);
		// if this streetligt is on a bridge or tunnel, there are no pedestrians since they don't go here
		if (city_params.num_cars == 0 && (on_bridge_or_tunnel || city_params.num_peds == 0)) {
			if (cached_smap) {dl_sources.back().assign_smap_id(uintptr_t(this)/sizeof(void *));} // cache on second frame
			cached_smap = 1;
		} else {cached_smap = 0;}
	}

	bool streetlight_t::proc_sphere_coll(point &center, float radius, vector3d const &xlate, vector3d *cnorm) const {
		point const p2(pos + xlate);
		float const pradius(get_streetlight_pole_radius());
		if (!dist_xy_less_than(p2, center, (pradius + radius))) return 0;
		return sphere_vert_cylin_intersect(center, radius, cylinder_3dw(p2, (p2 + vector3d(0.0, 0.0, get_streetlight_height())), pradius, 0.7*pradius), cnorm);
	}

	bool streetlight_t::line_intersect(point const &p1, point const &p2, float &t) const {
		float const pradius(get_streetlight_pole_radius());
		float t_new(0.0);
		if (line_int_cylinder(p1, p2, pos, (pos + vector3d(0.0, 0.0, get_streetlight_height())), pradius, 0.7*pradius, 0, t_new) && t_new < t) {t = t_new; return 1;}
		return 0;
	}
} // streetlight_ns


void streetlights_t::draw_streetlights(road_draw_state_t &dstate, bool shadow_only, bool always_on) const {
	if (streetlights.empty()) return;
	//timer_t t("Draw Streetlights");
	select_texture(WHITE_TEX);
	bool const is_local_shadow(camera_pdu.far_ < 0.1*FAR_CLIP); // Note: somewhat of a hack, but I don't have a better way to determine this
	for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {i->draw(dstate, shadow_only, is_local_shadow, always_on);}
}

void streetlights_t::add_streetlight_dlights(vector3d const &xlate, cube_t &lights_bcube, bool always_on) const {
	if (!always_on && !is_night(STREETLIGHT_ON_RAND)) return; // none of them can be on
	for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {i->add_dlight(xlate, lights_bcube, always_on);}
}

bool streetlights_t::proc_streetlight_sphere_coll(point &center, float radius, vector3d const &xlate, vector3d *cnorm) const {
	float const pradius(streetlight_ns::get_streetlight_pole_radius()), rtot(pradius + radius), y_end(center.y + rtot - xlate.y);

	for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {
		if (i->pos.y > y_end) break; // streetlights are sorted by y-val; if this holds, we can no longer intersect
		if (i->proc_sphere_coll(center, radius, xlate, cnorm)) return 1;
	}
	return 0;
}

bool streetlights_t::check_streetlight_sphere_coll_xy(point const &center, float radius, cube_t &coll_cube) const {
	float const pradius(streetlight_ns::get_streetlight_pole_radius()), rtot(pradius + radius), y_end(center.y + rtot);

	for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {
		if (i->pos.y > y_end) break; // streetlights are sorted by y-val; if this holds, we can no longer intersect
		if (!dist_xy_less_than(i->pos, center, rtot)) continue;
		coll_cube.set_from_point(i->pos);
		coll_cube.expand_by_xy(pradius);
		coll_cube.z2() += streetlight_ns::get_streetlight_height();
		return 1;
	}
	return 0;
}

bool streetlights_t::line_intersect_streetlights(point const &p1, point const &p2, float &t) const {
	bool ret(0);
	for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
	return ret;
}


road_isec_t::road_isec_t(cube_t const &c, int rx, int ry, unsigned char conn_, bool at_conn_road, bool has_stoplight_, short conn_to_city_) :
	cube_t(c), has_stoplight(has_stoplight_), num_conn(0), conn(conn_), conn_to_city(conn_to_city_), stoplight(at_conn_road)
{
	rix_xy[0] = rix_xy[1] = rx; rix_xy[2] = rix_xy[3] = ry; conn_ix[0] = conn_ix[1] = conn_ix[2] = conn_ix[3] = 0;
	if (conn == 15) {num_conn = 4;} // 4-way
	else if (conn == 7 || conn == 11 || conn == 13 || conn == 14) {num_conn = 3;} // 3-way
	else if (conn == 5 || conn == 6  || conn == 9  || conn == 10) {num_conn = 2;} // 2-way
	else {assert(0);}
	if (has_stoplight) {init_stoplights();}
}

tex_range_t road_isec_t::get_tex_range(float ar) const {
	switch (conn) {
	case 5 : return tex_range_t(1.0, 0.0, 0.0, 1.0, 0, 0); // 2-way: MX
	case 6 : return tex_range_t(0.0, 0.0, 1.0, 1.0, 0, 0); // 2-way: R0
	case 9 : return tex_range_t(1.0, 1.0, 0.0, 0.0, 0, 0); // 2-way: MXMY
	case 10: return tex_range_t(0.0, 1.0, 1.0, 0.0, 0, 0); // 2-way: MY
	case 7 : return tex_range_t(0.0, 0.0, 1.0, 1.0, 0, 0); // 3-way: R0
	case 11: return tex_range_t(1.0, 1.0, 0.0, 0.0, 0, 0); // 3-way: MY
	case 13: return tex_range_t(0.0, 1.0, 1.0, 0.0, 0, 1); // 3-way: R90MY
	case 14: return tex_range_t(1.0, 0.0, 0.0, 1.0, 0, 1); // 3-way: R90MX
	case 15: return tex_range_t(0.0, 0.0, 1.0, 1.0, 0, 0); // 4-way: R0
	default: assert(0);
	}
	return tex_range_t(0.0, 0.0, 1.0, 1.0); // never gets here
}

void road_isec_t::init_stoplights() {
	stoplight.init(num_conn, conn);
	max_eq(z2(), (z1() + max(get_street_sign_height(), stoplight_ns::stoplight_max_height()))); // include stoplights in Z-range
}
void road_isec_t::make_stop_signs() {
	has_stopsign = 1;
	max_eq(z2(), (z1() + max(get_street_sign_height(), get_stop_sign_height()))); // include stoplights in Z-range
}
void road_isec_t::make_4way(unsigned conn_to_city_) {
	num_conn = 4; conn = 15;
	assert(conn_to_city < 0); conn_to_city = conn_to_city_;
	has_stopsign = 0; has_stoplight = 1; // replace 3 stop signs with 4 stoplights
	init_stoplights(); // re-init with new conn
}
void road_isec_t::next_frame() {
	if (has_stoplight) {stoplight.next_frame();}
}
void road_isec_t::init_ssign_state(car_t const &car, ssign_state_t &ss, bool is_entering) const {
	ss.arrive_frame = frame_counter;
	ss.turn_dir     = car.turn_dir;
	ss.dest_orient  = get_dest_orient_for_car_in_isec(car, is_entering);
	ss.is_truck     = car.is_truck;
	ss.is_emergency = car.is_emergency;
	ss.in_use       = 1;
}
void road_isec_t::notify_waiting_car(car_t const &car) const { // can be called every frame, when waiting or in the intersection
	if (has_stoplight) {stoplight.notify_waiting_car(car.dim, car.dir, car.turn_dir);}
	
	if (has_stopsign) {
		ssign_state_pair_t &ss(get_ssign_state(car));

		if (car.stopped_at_light) { // stopped and waiting (at stop sign)
			if (!ss.waiting.in_use) {init_ssign_state(car, ss.waiting, 1);} // not yet marked as waiting; is_entering=1
		}
		else { // driving through intersection
			if (!ss.entering.in_use) { // not yet entered intersection
				//assert(ss.waiting.in_use); // car must have been waiting - unless it spawned here?
				init_ssign_state(car, ss.entering, 0); // is_entering=0
				ss.waiting.in_use = 0; // car is no longer waiting
			}
		}
	}
}
void road_isec_t::notify_leaving_car(car_t const &car) const { // called once per exit
	if (has_stopsign) {
		ssign_state_t &ss(get_ssign_state(car).entering);
		assert(ss.in_use); // must have been in the intersection
		ss = ssign_state_t(); // clear to zeros
	}
}
void road_isec_t::notify_turned_car(car_t const &car) const { // called once per turn completed
	if (has_stopsign) {init_ssign_state(car, get_ssign_state(car).entering, 0);} // is_entering=0; no in_use checking or setting
}
void road_isec_t::mark_crosswalk_in_use(bool dim, bool dir) const { // Note: const because in_use flag is mutable
	if (has_stoplight) {stoplight.mark_crosswalk_in_use(dim, dir);}
}
bool road_isec_t::can_go_based_on_light(car_base_t const &car) const {
	return (!has_stoplight || !red_or_yellow_light(car) || stoplight.can_turn_right_on_red(car)); // green light or right turn on red
}
bool road_isec_t::is_orient_currently_valid(unsigned orient, unsigned turn_dir) const {
	if (!(conn & (1<<orient))) return 0; // no road connection in this orient
	//if (!can_go_based_on_light(car)) return 1; // Note: can't call this without more car info
	//return stoplight.check_int_clear(orient, turn_dir); // intersection not clear (Note: too strong of a check, blocking car may be exiting the intersection)
	return 1;
}

unsigned road_isec_t::get_dest_orient_for_car_in_isec(car_base_t const &car, bool is_entering) const {
	unsigned const orient_in(car.get_orient_in_isec()); // invert dir (incoming, not outgoing)
	//cout << TXT(car.rot_z) << TXT(car.turn_val) << TXT(unsigned(car.turn_dir)) << TXT(car.dim) << TXT(car.dir) << TXT(orient_in) << hex << unsigned(conn) << dec << endl;
	if (is_entering) {assert(conn & (1<<orient_in));} // car must come from an enabled orient
	unsigned new_orient(0);
	switch (car.turn_dir) {
	case TURN_NONE:  new_orient = car.get_orient(); break;
	case TURN_LEFT:  new_orient = stoplight_ns::conn_left [orient_in]; break;
	case TURN_RIGHT: new_orient = stoplight_ns::conn_right[orient_in]; break;
	default: assert(0);
	}
	assert(conn & (1<<new_orient)); // car must go to an enabled orient
	return new_orient;
}

bool check_paths_cross(ssign_state_t const &ss, unsigned ss_cur_orient, unsigned cur_orient, unsigned dest_orient, unsigned turn_dir, bool is_truck) {
	if (ss.dest_orient == dest_orient) return 1; // same destination
	
	// check if at least one is a truck, and moving in the same general direction
	if ((is_truck || ss.is_truck) && (ss.dest_orient == cur_orient || ss_cur_orient == dest_orient || ss_cur_orient == cur_orient)) {
		// trucks are larger and may clip through each other when one is making a right turn opposite another making a left turn
		if (turn_dir == TURN_LEFT && ss.turn_dir == TURN_RIGHT) return 1;
		if (ss.turn_dir == TURN_LEFT && turn_dir == TURN_RIGHT) return 1;
	}
	if (ss.turn_dir == TURN_RIGHT || turn_dir == TURN_RIGHT) return 0; // can't cross if either is making a right turn
	if (ss.turn_dir == TURN_LEFT  || turn_dir == TURN_LEFT ) return 1; // otherwise, must cross if either is making a left turn
	return ((unsigned(ss.dest_orient)>>1) != (dest_orient>>1)); // both cars going straight; they cross if their dims are different
}
void update_arrive_frame(int &arrive_frame, bool is_emergency) {
	if (is_emergency) {arrive_frame -= 60*TICKS_PER_SECOND;} // emergency vehicles have priority over non-emergency vehicles; subtract 60s from their arrive frame
}
bool road_isec_t::can_go_now(car_t const &car) const {
	if (has_stopsign) { // run stop sign logic
		if (!car.is_emergency) { // emergency vehicles don't need to stop at stop signs
			if (!car.stopped_for_ssign) return 0; // must stop at stop sign first
			bool const allow_rolling_stop(bool(car.model_id & 1) ^ bool(car.color_id & 1) ^ car.dim ^ car.dir); // random-ish, but consistent for each car + stopsign pair
			if (!allow_rolling_stop && car.get_wait_time_secs() < 0.25) return 0; // must wait a quarter second before going
		}
		// run logic to check other cars at other stop signs or in the intersection
		unsigned const cur_orient(car.get_orient_in_isec()), dest_orient(get_dest_orient_for_car_in_isec(car, 1)); // is_entering=1
		assert(conn & (1<<cur_orient));
		int arrive_frame(ssign_state[cur_orient].waiting.arrive_frame);
		update_arrive_frame(arrive_frame, car.is_emergency);

		for (unsigned n = 0; n < 4; ++n) { // test each intersection connection
			if (!(conn & (1<<n))) continue; // no connection in this orient, skip
			if (n == cur_orient)  continue; // skip our own orient, assuming no two cars can be at the same place
			ssign_state_pair_t const &ss(ssign_state[n]);
			// check entering and waiting slots; if another car is waiting, check for earlier arrival time; if tied, use the orient as a tie breaker
			if (ss.entering.in_use && check_paths_cross(ss.entering, n, cur_orient, dest_orient, car.turn_dir, car.is_truck)) return 0;
			int waiting_arrive_trame(ss.waiting.arrive_frame);
			update_arrive_frame(waiting_arrive_trame, ss.waiting.is_emergency);
			if (ss.waiting .in_use && (waiting_arrive_trame < arrive_frame || (waiting_arrive_trame == arrive_frame && n < cur_orient)) &&
				check_paths_cross(ss.waiting, n, cur_orient, dest_orient, car.turn_dir, car.is_truck)) return 0;
		} // for n
		return 1;
	}
	if (!has_stoplight) return 1;
	// emergency vehicles can run red lights; can this cause collisions? usually the vehicle is stuck behind another one stopped at the light
	if (!car.is_emergency && !can_go_based_on_light(car)) return 0;
	if (stoplight.check_int_clear(car)) return 1;
	if ((frame_counter&15) == 0) {car.honk_horn_if_close_and_fast();} // honk every so often
	return 0; // intersection not clear
}

bool road_isec_t::has_left_turn_signal(unsigned orient) const {
	if (!has_stoplight) return 0; // never
	if (num_conn == 4) return 1; // always
	assert(num_conn == 3);
	switch (conn) {
	case 7 : return (orient == 1 || orient == 2);
	case 11: return (orient == 0 || orient == 3);
	case 13: return (orient == 0 || orient == 2);
	case 14: return (orient == 1 || orient == 3);
	default: assert(0);
	}
	return 0;
}

cube_t road_isec_t::get_stoplight_cube(unsigned n) const { // Note: mostly duplicated with draw_stoplights(), but difficult to factor the code out and share it
	assert(has_stoplight);
	assert(conn & (1<<n));
	float const sz(0.03*city_params.road_width), h(1.0*sz);
	bool const dim((n>>1) != 0), dir((n&1) == 0), side((dir^dim^1) != 0); // Note: dir is inverted here to represent car dir
	float const zbot(z1() + 2.0*h), dim_pos(d[dim][!dir] + SLIGHT_DIST_TO_CORNER_SCALE*(dir ? sz : -sz)); // pos in road dim
	float const side_len(side ? -sz : sz), v1(d[!dim][side] + (SLIGHT_DIST_TO_CORNER_SCALE - 1.0)*side_len), v2(v1 + side_len); // pos in other dim
	unsigned const num_segs(has_left_turn_signal(n) ? 6 : 3);
	float const sl_top(zbot + 1.2*h*num_segs), sl_lo(min(v1, v2) - 0.25*sz), sl_hi(max(v1, v2) + 0.25*sz);
	cube_t c;
	c.z1() = z1(); c.z2() = sl_top;
	c.d[ dim][0] = dim_pos - (dir ? -0.04 : 0.5)*sz; c.d[dim][1] = dim_pos + (dir ? 0.5 : -0.04)*sz;
	c.d[!dim][0] = sl_lo; c.d[!dim][1] = sl_hi;
	return c;
}

point road_isec_t::get_stop_sign_pos(unsigned n) const {
	float const dist_to_edge(0.78), dist_to_center(1.0 - dist_to_edge);
	bool const dim((n>>1) != 0), dir((n&1) == 0), side((dir^dim^1) != 0); // Note: dir is inverted here to represent car dir
	point pos(xc(), yc(), z1());
	pos[ dim] = dist_to_center*pos[ dim] + dist_to_edge*d[ dim][!dir];
	pos[!dim] = dist_to_center*pos[!dim] + dist_to_edge*d[!dim][side];
	return pos;
}

bool road_isec_t::check_sphere_coll(point const &pos, float radius, cube_t &coll_cube) const { // used for peds
	//if (has_stopsign) {} // not handled here
	if (!has_stoplight) return 0; // no stoplights
	if (!sphere_cube_intersect_xy(pos, radius, *this)) return 0;

	for (unsigned n = 0; n < 4; ++n) {
		if (!(conn & (1 << n))) continue; // no road in this dir
		cube_t const sl_cube(get_stoplight_cube(n));
		if (!sphere_cube_intersect(pos, radius, sl_cube)) continue;
		coll_cube = sl_cube;
		return 1;
	}
	return 0;
}

void road_isec_t::add_stoplight_bcubes_in_region(cube_t const &region, vect_cube_t &bcubes) const {
	//if (has_stopsign) {} // not handled here
	if (!has_stoplight) return; // no stoplights
	if (!intersects_xy(region)) return;

	for (unsigned n = 0; n < 4; ++n) {
		if (!(conn & (1 << n))) continue; // no road in this dir
		cube_t const sl_cube(get_stoplight_cube(n));
		if (sl_cube.intersects_xy(region)) {bcubes.push_back(sl_cube);}
	}
}

bool road_isec_t::proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, float dist, vector3d *cnorm) const {
	//if (has_stopsign) {} // not handled here
	if (!has_stoplight) return 0; // no stoplights
	if (!moving_sphere_cube_intersect_xy(pos, p_last, (*this + xlate), dist, radius)) return 0;

	for (unsigned n = 0; n < 4; ++n) {
		if (!(conn & (1<<n))) continue; // no road in this dir
		if (sphere_cube_int_update_pos(pos, radius, (get_stoplight_cube(n) + xlate), p_last, 0, cnorm)) return 1; // typically only one coll, just return the first one
	}
	return 0;
}

bool check_line_clip_update_t(point const &p1, point const &p2, float &t, cube_t const &c) {
	float tmin(0.0), tmax(1.0);
	if (get_line_clip(p1, p2, c.d, tmin, tmax) && tmin < t) {t = tmin; return 1;}
	return 0;
}

bool road_isec_t::line_intersect(point const &p1, point const &p2, float &t) const {
	//if (has_stopsign) {} // not handled here
	if (!has_stoplight) return 0; // no stoplights
	if (!line_intersects(p1, p2)) return 0;
	bool ret(0);

	for (unsigned n = 0; n < 4; ++n) {
		if (conn & (1<<n)) {ret |= check_line_clip_update_t(p1, p2, t, get_stoplight_cube(n));}
	}
	return ret;
}

void road_isec_t::draw_sl_block(quad_batch_draw &qbd, draw_state_t &dstate, point p[4], float h, unsigned state, bool draw_unlit, float flare_alpha, vector3d const &n, tex_range_t const &tr) const {
	for (unsigned j = 0; j < 3; ++j) {
		colorRGBA const &color(stoplight_ns::stoplight_colors[j]);

		if (j == state) {
			qbd.add_quad_pts(p, color, n, tr);
			if (flare_alpha > 0.0) {dstate.add_light_flare(0.25*(p[0] + p[1] + p[2] + p[3]), n, color, flare_alpha, 2.0*h);}
		}
		else if (draw_unlit) {
			qbd.add_quad_pts(p, (color + WHITE)*0.05, n, tr);
		}
		for (unsigned e = 0; e < 4; ++e) {p[e].z += 1.2*h;}
	}
}

// should this be a function of road_draw_state_t that takes the road_isec_t instead?
void road_isec_t::draw_stoplights_and_street_signs(road_draw_state_t &dstate, vector<road_t> const &roads, unsigned cur_city, bool shadow_only) const {
	if (!(has_stoplight || has_stopsign)) return; // no stoplights or street signs
	if (!dstate.check_cube_visible(*this, 0.16)) return; // dist_scale=0.16
	point const center(get_cube_center() + dstate.xlate);
	float const dist_val(shadow_only ? 0.0 : p2p_dist(camera_pdu.pos, center)/dstate.draw_tile_dist);
	bool const draw_street_sign(dist_val < 0.06);
	if (!(has_stoplight || draw_street_sign)) return; // nothing to draw
	vector3d const cview_dir(camera_pdu.pos - center);
	float const sz(0.03*city_params.road_width), h(1.0*sz);
	float const ss_height(get_stop_sign_height()), sign_z1(z1() + get_street_sign_height());
	color_wrapper cw(BLACK);

	for (unsigned n = 0; n < 4; ++n) { // {-x, +x, -y, +y} = {W, E, S, N} facing = car traveling {E, W, N, S}
		if (!(conn & (1<<n))) continue; // no road in this dir
		bool const dim((n>>1) != 0), dir((n&1) == 0), side((dir^dim^1) != 0); // Note: dir is inverted here to represent car dir
		float const side_len(side ? -sz : sz);
		cube_t sc; // stoplight/sign cube

		if (has_stoplight) {
			float const zbot(z1() + 2.0*h), dim_pos(d[dim][!dir] + SLIGHT_DIST_TO_CORNER_SCALE*(dir ? sz : -sz)); // pos in road dim
			float const v1(d[!dim][side] + (SLIGHT_DIST_TO_CORNER_SCALE - 1.0)*side_len), v2(v1 + side_len); // pos in other dim
			// draw base
			unsigned const num_segs(has_left_turn_signal(n) ? 6 : 3);
			float const sl_height(1.2*h*max(float(num_segs), 4.2f)); // use 4.2 segs rather than 3 so that the street sign isn't blocked by the crosswalk sign
			float const sl_top(zbot + sl_height), sl_lo(min(v1, v2) - 0.25*sz), sl_hi(max(v1, v2) + 0.25*sz);

			if (!draw_street_sign) { // draw front face only
				point pts[4];
				pts[0][dim]  = pts[1][ dim] = pts[2][dim] = pts[3][dim] = dim_pos; // plane of the face
				pts[0][!dim] = pts[3][!dim] = sl_lo;
				pts[1][!dim] = pts[2][!dim] = sl_hi;
				pts[0].z = pts[1].z = z1();
				pts[2].z = pts[3].z = sl_top;
				dstate.qbd_sl.add_quad_pts(pts, cw, (dim ? (dir ? plus_y : -plus_y) : (dir ? plus_x : -plus_x))); // Note: normal doesn't really matter since color is black
			}
			else {
				set_cube_zvals(sc, z1(), sl_top);
				sc.d[ dim][0] = dim_pos - (dir ? -0.04 : 0.5)*sz; sc.d[dim][1] = dim_pos + (dir ? 0.5 : -0.04)*sz;
				sc.d[!dim][0] = sl_lo; sc.d[!dim][1] = sl_hi;
				dstate.draw_cube(dstate.qbd_sl, sc, cw, 1); // skip_bottom=1; Note: uses traffic light texture, but color is black so it's all black anyway

				if (!shadow_only && tt_fire_button_down && game_mode != GAME_MODE_FPS) { // player debug visualization
					point const p1(camera_pdu.pos - dstate.xlate), p2(p1 + camera_pdu.dir*FAR_CLIP);
					if (sc.line_intersects(p1, p2)) {dstate.set_label_text(stoplight.label_str(), (sc.get_cube_center() + dstate.xlate));}
				}
				if (sc.z2() < sign_z1) { // sign is too low, add a metal pole to extend it upward
					sc.expand_in_dim(dim, -0.1*sc.get_sz_dim(dim)); // small shrink
					set_wall_width(sc, sc.get_center_dim(!dim), 0.5*sc.get_sz_dim(dim), !dim); // make it centered and square
					set_cube_zvals(sc, sc.z2(), sign_z1);
					dstate.draw_cube(dstate.qbd_untextured, sc, GRAY);
				}
			}
			bool const draw_detail(dist_val < 0.05); // add flares only when very close

			// Note skip crosswalk signal for road to the right of an unconnected 3-way int (T-junction)
			// because it's not needed and will never get to walk (due to opposing left turn light), and there's no stoplight to attach it to
			if (dist_val <= 0.06 && (conn & (1<<stoplight_ns::to_right[2*dim + (!dir)])) != 0) { // draw crosswalk signal
				int const cw_state(shadow_only ? 0 : stoplight.get_crosswalk_state(dim, !dir)); // Note: backwards dir
				colorRGBA cw_color(stoplight_ns::crosswalk_colors[cw_state] * 0.6); // less bright

				if (cw_state == stoplight_ns::CW_WARN) { // flashing
					float const flash_period = 0.5; // in seconds
					double const time(fract(tfticks/(double(flash_period)*TICKS_PER_SECOND)));
					if (time > 0.5) {cw_color = BLACK;}
				}
				float const cw_dim_pos(dim_pos + 0.7*(dir ? sz : -sz)); // location in road dim
				cube_t c;
				c.z1() = zbot + 1.2*h + 1.4*h*dim; c.z2() = c.z1() + 1.6*h;
				c.d[dim][0] = cw_dim_pos - 0.5*sz; c.d[dim][1] = cw_dim_pos + 0.5*sz;

				for (unsigned i = 0; i < 2; ++i) { // opposite sides of the road
					float const ndp(d[!dim][i] - (SLIGHT_DIST_TO_CORNER_SCALE + 0.2)*(i ? sz : -sz));
					c.d[!dim][0] = ndp - 0.1*sz; c.d[!dim][1] = ndp + 0.1*sz;
					dstate.draw_cube(dstate.qbd_sl, c, cw, shadow_only); // skip_bottom=shadow_only

					if (!shadow_only && cw_color != BLACK) { // draw light
						point p[4];
						p[0][!dim] = p[1][!dim] = p[2][!dim] = p[3][!dim] = c.d[!dim][!i] + 0.01*(i ? -sz : sz);
						p[0][dim] = p[3][dim] = c.d[dim][0]; p[1][dim] = p[2][dim] = c.d[dim][1];
						p[0].z = p[1].z = c.z1(); p[2].z = p[3].z = c.z2();
						vector3d const normal(vector_from_dim_dir(!dim, !i));
						point const cw_center(0.25*(p[0] + p[1] + p[2] + p[3]));
						vector3d const cw_cview_dir(camera_pdu.pos - (cw_center + dstate.xlate));
						float dp(dot_product(normal, cw_cview_dir));

						if (dp > 0.0) { // if back facing, don't draw the lights
							float const mag(CLIP_TO_01(2.5f*dp/cw_cview_dir.mag() - 1.0f)); // normalize and strengthen slope

							if (mag > 0.0) {
								bool const is_walk(cw_state == stoplight_ns::CW_WALK);
								dstate.qbd_sl.add_quad_pts(p, WHITE*mag, normal, tex_range_t((is_walk ? 0.25 : 0.0), 0.0, (is_walk ? 0.5 : 0.25), 0.5)); // color stored in texture
								if (draw_detail) {dstate.add_light_flare(cw_center, normal, cw_color*mag, 1.0, 0.8*h);}
							}
						}
					}
				} // for i
			} // end crosswalk signal
			if (!shadow_only && dist_val < 0.1) { // draw light decals when not in the shadow pass
				vector3d const normal(vector_from_dim_dir(dim, !dir));
				
				if (dot_product(normal, cview_dir) > 0.0) { // only draw lights if front facing
					// draw straight/line turn light
					point p[4];
					p[0][dim] = p[1][dim] = p[2][dim] = p[3][dim] = dim_pos;
					p[0][!dim] = p[3][!dim] = v1; p[1][!dim] = p[2][!dim] = v2;
					p[0].z = p[1].z = zbot; p[2].z = p[3].z = zbot + h;
					draw_sl_block(dstate.qbd_sl, dstate, p, h, stoplight.get_light_state(dim, dir, TURN_NONE), draw_detail, draw_detail, normal, tex_range_t(0.0, 0.5, 0.5, 1.0));

					if (has_left_turn_signal(n)) { // draw left turn light (upper light section)
						draw_sl_block(dstate.qbd_sl, dstate, p, h, stoplight.get_light_state(dim, dir, TURN_LEFT), draw_detail, 0.5*draw_detail, normal, tex_range_t(1.0, 0.5, 0.5, 1.0));
					}
				}
			}
		} // end stoplight drawing
		else { // stop sign case; add a vertical pole
			point spos(get_stop_sign_pos(n));
			spos[dim] += (dir ? 1.0 : -1.0)*0.07*ss_height; // move a bit behind the stop sign; must match stopsign_t::get_bird_bcube()
			sc.set_from_point(spos);
			sc.expand_by_xy(0.02*ss_height);
			set_cube_zvals(sc, z1(), sign_z1);
			if (draw_street_sign) {dstate.draw_cube(dstate.qbd_untextured, sc, GRAY);} // draw the sign pole if close
#if 0 // code for debugging stop sign logic
			colorRGBA const color(ssign_state[n].entering.in_use, ssign_state[n].waiting.in_use, 0, 1); // entering=RED, waiting=GREEN, both=YELLOW, neither=BLACK
			cube_t c(sc);
			c.expand_by_xy(0.1*ss_height);
			set_cube_zvals(c, sc.z2(), sc.z2()+0.2*ss_height);
			dstate.draw_cube(dstate.qbd_untextured, c, color);
#endif
		}
		if (draw_street_sign) { // draw street sign
			// draw green sign cube
			assert(sc.is_strictly_normalized());
			cube_t sign(sc); // start with extended (pole) stoplight or stop sign cube
			set_cube_zvals(sign, sc.z2()-0.7*sz, sc.z2()-0.1*sz); // high up
			sign.d[ dim][ dir ] -= (dir ? 1.0 : -1.0)*0.8*sc.get_sz_dim(dim); // shrink
			sign.d[!dim][ side]  = sc.d[!dim][!side]; // flush with the side of the stoplight body or stop sign pole
			sign.d[!dim][!side] += 5.0*side_len; // extend into the road
			dstate.draw_cube(dstate.qbd_untextured, sign, colorRGBA(0.0, 0.6, 0.0, 1.0));

			if (!shadow_only && dist_val < 0.05) { // draw sign text when very close
				unsigned const to_the_right[4] = {2,3,1,0};
				int const road_ix(rix_xy[to_the_right[n]]); // map from the current road to the one on the right
				string name;

				if (road_ix >= 0) {
					assert((unsigned)road_ix < roads.size());
					name = roads[road_ix].get_name(cur_city);
				}
				else { // city connector road
					assert(conn_to_city >= 0);
					unsigned const city1(min(cur_city, (unsigned)conn_to_city)), city2(max(cur_city, (unsigned)conn_to_city));
					name = road_t().get_name(city1 + (city2 << 16)); // use the canonical pair of city indices so that the road name matches at each city end
				}
				bool const text_dir(sign.d[dim][dir] + dstate.xlate[dim] < camera_pdu.pos[dim]); // draw text on the side facing the player
				// split the road name from the extension and draw them in different sizes, like in real street signs
				auto space_start(name.find_last_of(" "));
				assert(space_start != string::npos);
				string const road_name(name.substr(0, space_start)), ext(name.substr(space_start+1));
				cube_t sign_main(sign), sign_ext(sign);
				bool const ext_dir(text_dir ^ dim);
				float const signed_sz(ext_dir ? sz : -sz);
				sign_main.d[!dim][ ext_dir] -= 1.0*signed_sz; // leave space for extension
				sign_ext .d[!dim][!ext_dir] += 4.0*signed_sz; // leave space for street name (full size is 5.0x)
				sign_ext .d[!dim][ ext_dir] -= 0.4*signed_sz; // gap to pole
				sign_ext.z1() += 0.15*sz; // upper part / smaller text
				add_sign_text_verts(road_name, sign_main, dim, text_dir, WHITE, dstate.text_verts);
				add_sign_text_verts(ext,       sign_ext,  dim, text_dir, WHITE, dstate.text_verts); // smaller text
				//add_sign_text_verts(name, sign, dim, text_dir, WHITE, dstate.text_verts); // single font size text
			}
			if (hospital_dir & (3 << 2*n)) { // maybe draw hospital sign on the stoplight pole
				point sign_center;
				sign_center[ dim] = sc.d[dim][!dir] + (dir ? -1.0 : 1.0)*0.02*sz;
				sign_center[!dim] = sc.get_center_dim(!dim);
				sign_center.z     = sc.z2() - 0.85*sz;
				vector3d const normal(vector_from_dim_dir(dim, !dir)), dx(0.5*sz*vector_from_dim_dir(!dim, !side)), dy(0.8*sz*plus_z);

				if (!shadow_only && dot_product((camera_pdu.pos - (sign_center + dstate.xlate)), normal) > 0.0) { // front facing non-shadow
					bool const arrow_dir(hospital_dir & (2 << 2*n)); // check upper bit
					dstate.qbd_hospital.add_quad_dirs(sign_center, (arrow_dir ? -dx : dx), dy, WHITE, normal);
				}
				else { // back facing or shadow pass
					dstate.qbd_untextured.add_quad_dirs(sign_center, dx, dy, GRAY, -normal);
				}
			}
		} // end street sign
	} // for n
}


float road_connector_t::get_player_zval(point const &center, cube_t const &c) const {
	float const t((center[dim] - c.d[dim][0])/c.get_sz_dim(dim));
	float const za(slope ? c.z2() : c.z1()), zb(slope ? c.z1() : c.z2()), zval(za + (zb - za)*t); // z-value at x/y location
	return zval - src_road.get_z_adj() - ROAD_HEIGHT; // place exactly on mesh under the road/bridge/tunnel surface
}

void road_connector_t::add_streetlights(unsigned num_per_side, bool staggered, float dn_shift_mult, float za, float zb) {
	streetlights.reserve(2*num_per_side);
	float const dsz(get_sz_dim(dim)), dnsz(get_sz_dim(!dim)), dn_shift(dn_shift_mult*dnsz);
	float const dz(-(0.1 + 2.0*get_slope_val())*streetlight_ns::get_streetlight_pole_radius());
	za += dz; zb += dz; // shift down slightly since we may be on a slope and don't want to leave a gap
	if (staggered) {num_per_side *= 2;}

	for (unsigned n = 0; n < num_per_side; ++n) {
		float const v((n + 0.5)/num_per_side); // 1/8, 3/8, 5/8, 7/8
		point pos1, pos2;
		pos1[ dim] = pos2[dim] = d[dim][0] + v*dsz;
		pos1[!dim] = d[!dim][0] - dn_shift; pos2[!dim] = d[!dim][1] + dn_shift;
		pos1.z = pos2.z = za + v*(zb - za);
		vector3d dir(zero_vector); dir[!dim] = 1.0;
		if (!staggered || (n&1) == 0) {streetlights.emplace_back(pos1,  dir, 1);} // on_bridge_or_tunnel=1
		if (!staggered || (n&1) == 1) {streetlights.emplace_back(pos2, -dir, 1);} // on_bridge_or_tunnel=1
	} // for n
	sort_streetlights_by_yx();
}


bool bridge_t::proc_sphere_coll(point &center, point const &prev, float sradius, float prev_frame_zval, vector3d const &xlate, vector3d *cnorm) const {
	if (proc_streetlight_sphere_coll(center, sradius, xlate, cnorm)) return 1;
	float const exp(0.5*sradius);
	bool const prev_int(contains_pt_xy_exp((prev - xlate), exp));
	if (!prev_int && !src_road.contains_pt_xy_exp((center - xlate), exp)) return 0; // no x/y containment for road segment (cur or prev)
	cube_t const c(src_road + xlate);
	float const zval(get_player_zval(center, c));
	// Note: we need to use prev_frame_zval for camera collisions because center.z will always start at the mesh zval;
	// we can't just always place the camera on the bridge because the player may be walking on the ground under the bridge
	if (center.z - sradius > zval || prev_frame_zval + sradius < zval) return 0; // completely above or below the road/bridge
	center.z = zval + sradius; // place exactly on the road/bridge surface
	if (prev_int) {center[!dim] = min(c.d[!dim][1], max(c.d[!dim][0], center[!dim]));} // keep the sphere from falling off the bridge (approximate; assumes bridge has walls)
	if (cnorm) {*cnorm = plus_z;} // approximate - assume collision with bottom surface of road (intended for player)
	return 1;
}

bool bridge_t::line_intersect(point const &p1, point const &p2, float &t) const {
	float const size_scale(city_params.road_width), wall_width(0.05*size_scale), wall_height(0.1*size_scale);
	cube_t road_surface(*this); // start from bcube
	// set z1/z2 to include the road surface and guardrails; approximate since the bridge may be slightly sloped
	set_cube_zvals(road_surface, (min(z1(), z2()) - wall_width - 0.25f*ROAD_HEIGHT), (max(z1(), z2()) + wall_height));
	road_surface.expand_by_xy(wall_width);
	return check_line_clip_update_t(p1, p2, t, road_surface);
}


void tunnel_t::init(point const &start, point const &end, float radius_, float end_length, bool dim) {
	radius = radius_;
	height = (1.0 + TUNNEL_WALL_THICK)*radius;
	vector3d const dir((end - start).get_norm());
	point const extend[2] = {(start + dir*end_length), (end - dir*end_length)}; // extend toward the tunnel interior
	ends[0].set_from_point(start);
	ends[1].set_from_point(end);

	for (unsigned d = 0; d < 2; ++d) {
		ends[d].z2() += height; // height
		ends[d].d[!dim][0] -= radius; // width
		ends[d].d[!dim][1] += radius; // width
		ends[d].union_with_pt(extend[d]); // length
		facade_height[d] = 2.0*TUNNEL_WALL_THICK*radius; // min value - will likely be increased later
	}
	get_bcube() = ends[0]; get_bcube().union_with_cube(ends[1]); // bounding cube of both ends - overwrite road bcube
}

cube_t tunnel_t::get_tunnel_bcube() const {
	cube_t bcube(*this);
	bcube.expand_by_xy(radius);
	bcube.z2() += max(facade_height[0], facade_height[1]); // should be close enough
	return bcube;
}

void tunnel_t::calc_top_bot_side_cubes(cube_t cubes[4]) const { // {bottom, top, left, right}
	float const wall_thick(TUNNEL_WALL_THICK*city_params.road_width), end_ext(2.0*(dim ? DY_VAL : DX_VAL));
	cube_t tc(*this);
	tc.d[dim][0] -= end_ext; tc.d[dim][1] += end_ext; // extend to cover the gaps in the mesh
	cubes[0] = cubes[1] = cubes[2] = cubes[3] = tc;
	cubes[0].z1() = -wall_thick; // bottom cube extends below the road surface
	cubes[0].z2() = -ROAD_HEIGHT; // top of bottom cube is just below the road surface
	cubes[1].z1() = height - wall_thick; // top cube
	cubes[1].z2() = height;
	cubes[2].z1() = cubes[3].z1() = cubes[0].z2(); // bottom of sides meet bottom cube
	cubes[2].z2() = cubes[3].z2() = cubes[1].z1(); // top of sides meet top cube
	cubes[2].d[!dim][1] = cubes[2].d[!dim][0] + wall_thick; // left side
	cubes[3].d[!dim][0] = cubes[3].d[!dim][1] - wall_thick; // right side
}

bool tunnel_t::proc_sphere_coll(point &center, point const &prev, float sradius, float prev_frame_zval, vector3d const &xlate, vector3d *cnorm) const {
	if (proc_streetlight_sphere_coll(center, sradius, xlate, cnorm)) return 1;
	cube_t tunnel_ext(*this); // including ends
	tunnel_ext.expand_in_dim(dim, get_end_ext());
	bool const prev_int(tunnel_ext.contains_pt_xy(prev - xlate));
	if (!prev_int && !tunnel_ext.contains_pt_xy(center - xlate)) return 0; // no x/y containment for road segment (cur or prev)
	cube_t const c(src_road + xlate);
	float const zval(get_player_zval(center, c));
	// Note: we need to use prev_frame_zval for camera collisions because center.z will always start at the mesh zval;
	// we can't just always place the camera on the tunnel because the player may be walking on the ground above the tunnel
	if (prev_frame_zval - sradius > zval + height) return 0; // completely above the tunnel
	center.z = zval + sradius; // place exactly on mesh under the road/tunnel surface
	if (prev_int) {center[!dim] = min(c.d[!dim][1], max(c.d[!dim][0], center[!dim]));} // keep the sphere inside the tunnel (approximate)
	if (cnorm) {*cnorm = plus_z;} // approximate - assume collision with bottom surface of road (intended for player)
	return 1;
}

float tunnel_t::get_end_ext() const {return (2.0*(dim ? DY_VAL : DX_VAL));}

void tunnel_t::calc_zvals_and_eext(float &zf, float &zb, float &end_ext) const {
	zf = ends[0].z1(); zb = ends[1].z1();
	end_ext = get_end_ext();
	float const dz_ext(end_ext*(zb - zf)/get_length());
	zf -= dz_ext; zb += dz_ext; // adjust zvals for extension
}

bool tunnel_t::line_intersect(point const &p1, point const &p2, float &t) const {
	cube_t const bcube(get_tunnel_bcube());
	if (!check_line_clip(p1, p2, bcube.d)) return 0;
	cube_t cubes[4]; // {bottom, top, left, right}
	calc_top_bot_side_cubes(cubes);
	float const wall_thick(TUNNEL_WALL_THICK*city_params.road_width), width(max(0.5*get_width(), 2.0*(dim ? DX_VAL : DY_VAL)));
	float zf, zb, end_ext;
	calc_zvals_and_eext(zf, zb, end_ext);
	bool const d(dim);
	bool ret(0);

	// check tunnel {bottom, top, left, right} as tilted cubes; this really isn't needed since line intersection tests are broken when under the terrain anyway
	for (unsigned i = 0; i < 4; ++i) {
		cube_t const &c(cubes[i]);
		point p[8];
		draw_state_t::set_cube_pts(c, c.z1()+zf, c.z1()+zb, c.z2()+zf, c.z2()+zb, d, 0, p);
		if (i == 0) {ret |= line_quad_intersect(p1, p2, p+4, t);} // top
		if (i == 1) {ret |= line_quad_intersect(p1, p2, p+0, t);} // bottom
		if (i == 3) {point const pts1[4] = {p[1], p[2], p[6], p[5]}; ret |= line_quad_intersect(p1, p2, pts1, t);} // left
		if (i == 2) {point const pts2[4] = {p[3], p[0], p[4], p[7]}; ret |= line_quad_intersect(p1, p2, pts2, t);} // right
	} // for i
	for (unsigned n = 0; n < 2; ++n) { // check tunnel facades (Note: similar to code in road_draw_state_t::draw_tunnel())
		cube_t const &tend(ends[n]);
		cube_t c(tend);
		c.z1() += height - 0.5*wall_thick; // tunnel ceiling
		c.z2() += facade_height[n]; // high enough to cover the hole in the mesh
		c.d[d][0] -= 0.9*end_ext; c.d[d][1] += 0.9*end_ext; // extend to cover the gaps in the mesh (both dirs) - slightly less than interior so that it sticks out
		ret |= check_line_clip_update_t(p1, p2, t, c);
		c.z1() = ends[n].z1() - 2.0*wall_thick; // extend below the bottom to cover any gaps at the corners
		c.d[!d][0] = ends[n].d[!d][0] - width; c.d[!d][1] = tend.d[!d][0]; // left side
		ret |= check_line_clip_update_t(p1, p2, t, c);
		c.d[!d][1] = ends[n].d[!d][1] + width; c.d[!d][0] = tend.d[!d][1]; // right side
		ret |= check_line_clip_update_t(p1, p2, t, c);
	} // for n
	return ret;
}


road_mat_mgr_t road_mat_mgr;

void road_draw_state_t::draw_city_region_int(quad_batch_draw &cache, unsigned type_ix) {
	if (cache.verts.empty()) return; // nothing to draw

	if (emit_now) { // draw shadow blocks directly
		road_mat_mgr.set_texture(type_ix);
		cache.draw();
	} else {qbd_batched[type_ix].add_quads(cache);} // add non-shadow blocks for drawing later
}

void road_draw_state_t::draw_city_skirt(cube_t const &bcube) {
	cube_t skirt(bcube);
	skirt.z1() -= 0.2*city_params.road_width; // extend into the ground
	draw_cube(qbd_skirt, skirt, WHITE, 1, 16.0, 4); // skip zvals
	for (auto &v : qbd_skirt.verts) {v.set_norm(-plus_z);} // hack to avoid shadow artifacts: point the normal downward so that there is no light
	road_mat_mgr.set_texture(TYPE_ROAD_SKIRT);
	qbd_skirt.draw_and_clear();
}

void road_draw_state_t::pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only, bool always_setup_shader) {
	draw_state_t::pre_draw(xlate_, use_dlights_, shadow_only, always_setup_shader);
	ar = city_params.get_road_ar();
}

void road_draw_state_t::draw_unshadowed() {
	for (unsigned i = 0; i < NUM_RD_TIDS; ++i) { // only unshadowed blocks
		road_mat_mgr.set_texture(i);
		qbd_batched[i].draw_and_clear();
	}
}

void road_draw_state_t::post_draw() {
	draw_state_t::post_draw();
	if (qbd_sl.empty() && qbd_emissive.empty()) return;
	shader_t s;
	
	if (shadow_only) {s.begin_shadow_map_shader();}
	else {
		s.begin_simple_textured_shader(); // no lighting
		road_mat_mgr.set_stoplight_texture();
	}
	if (!qbd_sl.empty()) { // have stoplights to draw
		set_std_depth_func_with_eq(); // helps prevent Z-fighting
		qbd_sl.draw_and_clear();
		set_std_depth_func();
	}
	if (!qbd_emissive.empty()) { // have streetlight emissive light spots to draw
		assert(!shadow_only);
		s.set_color_e(streetlight_ns::light_color);
		draw_and_clear_blur_qbd(qbd_emissive);
		s.clear_color_e();
	}
	s.end_shader();
}

void road_draw_state_t::end_cur_tile() {
	if (!qbd_hospital.empty()) { // draw hospital signs
		select_texture(get_texture_by_name("roads/hospital_arrow.png"));
		s.add_uniform_float("min_alpha", 0.9);
		enable_blend();
		qbd_hospital.draw_and_clear();
		disable_blend();
		s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA);
	}
	if (!shadow_only && !qbd_untextured.empty()) {select_texture(WHITE_TEX);}
	qbd_untextured.draw_and_clear();
	
	if (!text_verts.empty()) { // draw street names on signs; must be after qbd_untextured due to alpha blending
		assert(!shadow_only);
		text_drawer::bind_font_texture();
		enable_blend();
		draw_quad_verts_as_tris(text_verts);
		text_verts.clear();
		disable_blend();
	}
}

struct plot_region_t : public cube_t {
	tex_range_t tr;
	plot_region_t(cube_t const &c, tex_range_t const &tr_) : cube_t(c), tr(tr_) {}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(ar*tr.x1, ar*tr.y1, ar*tr.x2, ar*tr.y2);}
};

void road_draw_state_t::add_city_quad(road_seg_t  const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool) { // road segment or skirt
	bool const add_skirt(type_ix == TYPE_ROAD_SKIRT);
	r.add_road_quad(qbd, color, ar, add_skirt);
}
void road_draw_state_t::add_city_quad(road_t      const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool) { // tracks
	r.add_road_quad(qbd, color, ar/TRACKS_WIDTH);
}
void road_draw_state_t::add_city_quad(road_plot_t const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool draw_all) { // plots and parks
	if (!draw_all && (type_ix == TYPE_PARK) != r.is_park) return;
	
	if (!plot_cuts.empty()) {
		float const dx_inv(1.0/r.dx()), dy_inv(1.0/r.dy());
		vect_cube_t parts, temp;
		subtract_cubes_from_cube((cube_t)r, plot_cuts, parts, temp, 1); // zval_mode=1 (ignore zvals)

		for (cube_t const &p : parts) {
			tex_range_t const tr((p.x1() - r.x1())*dx_inv, (p.y1() - r.y1())*dy_inv, (p.x2() - r.x1())*dx_inv, (p.y2() - r.y1())*dy_inv);
			add_flat_city_quad(plot_region_t(p, tr), qbd, color, ar);
		}
	}
	else {add_flat_city_quad(r, qbd, color, ar);}
}

cube_t bridge_t::get_drawn_bcube() const {
	cube_t bcube(*this);
	float const scale(1.0*city_params.road_width);
	float const l_expand(2.0*(d ? DY_VAL : DY_VAL)); // slight expand along road dim so that we're sure to cover the entire gap
	float const w_expand(0.25*scale); // expand width to add space for supports
	bcube.d[ dim][0] -= l_expand; bcube.d[ dim][1] += l_expand;
	bcube.d[!dim][0] -= w_expand; bcube.d[!dim][1] += w_expand;
	max_eq(bcube.d[dim][0], src_road.d[dim][0]); // clamp to orig road segment length
	min_eq(bcube.d[dim][1], src_road.d[dim][1]);
	bcube.z2() += BRIDGE_HEIGHT_TO_LEN*bcube.get_sz_dim(dim); // make it higher
	return bcube; // approximate
}

void road_draw_state_t::draw_bridge(bridge_t const &bridge, bool shadow_only) { // Note: called rarely, so doesn't need to be efficient
	//timer_t timer("Draw Bridge"); // 0.065ms - 0.11ms
	cube_t const bcube(bridge.get_drawn_bcube());
	if (!check_cube_visible(bcube, 1.0)) return; // VFC/too far
	float const scale(1.0*city_params.road_width), w_expand(0.25*scale);
	unsigned const d(bridge.dim);
	point const cpos(camera_pdu.pos - xlate);
	float const center(bcube.get_center_dim(!d)), len(bcube.get_sz_dim(d)), peak_height(BRIDGE_HEIGHT_TO_LEN*len);
	point p1, p2; // centerline end points
	p1.z = bridge.get_start_z();
	p2.z = bridge.get_end_z();
	p1[!d] = p2[!d] = center;
	p1[ d] = bcube.d[d][0];
	p2[ d] = bcube.d[d][1];
	vector3d const delta(p2 - p1);
	colorRGBA const main_color(WHITE), cables_color(LT_GRAY), concrete_color(LT_GRAY);
	color_wrapper const cw_main(main_color), cw_cables(cables_color), cw_concrete(concrete_color);
	float const thickness(0.2*scale), conn_thick(0.25*thickness), cable_thick(0.1*thickness), wall_width(0.25*thickness), wall_height(0.5*thickness);
	point const closest_pt((bridge + xlate).closest_pt(camera_pdu.pos));
	float const dist_val(shadow_only ? 1.0 : p2p_dist(camera_pdu.pos, closest_pt)/draw_tile_dist);
	int const cable_ndiv(min(24, max(4, int(0.4/dist_val))));
	unsigned const num_segs(max(16U, min(48U, unsigned(ceil(2.5*len/scale))))); // scale to fit the gap, with reasonable ranges
	float const step_sz(1.0/num_segs), delta_d(step_sz*delta[d]), delta_z(step_sz*delta.z);
	float zvals[48+1] = {}, cur_zpos(p1.z), cur_dval(p1[d]); // resize zvals based on max num_segs
	vector<float> sm_split_pos;
	uint64_t prev_tile_id(0);
	point query_pt(bridge.get_cube_center());
	ensure_shader_active(); // needed for use_smap=0 case
	if (!shadow_only) {select_texture(WHITE_TEX);}

	for (unsigned n = 0; n <= num_segs; ++n) { // populate zvals and dvals
		float const t(n*step_sz), v(2.0f*fabs(t - 0.5f)), zpos(p1.z + delta.z*t);
		zvals[n] = zpos + peak_height*(1.0f - v*v) - 0.5f*scale;
	}
	for (unsigned n = 0; n < num_segs; ++n) { // add arches
		float const zval(zvals[n+1]), next_dval(cur_dval + delta_d);
		point pts[4], conn_pts[2];
		pts[0][d] = pts[3][d] = cur_dval;
		pts[1][d] = pts[2][d] = next_dval;

		if (!shadow_only) {
			query_pt[d] = 0.5f*(cur_dval + next_dval); // center point of this segment
			uint64_t const tile_id(get_tile_id_containing_point(query_pt + xlate));

			if (n == 0 || tile_id != prev_tile_id) { // first segment, or new tile for this segment
				qbd_bridge.draw_and_clear(); // flush
				begin_tile(query_pt, 1); // will_emit_now=1
				prev_tile_id = tile_id;
				if (n > 0) {sm_split_pos.push_back(query_pt[d]);} // record split point for splitting road surface and guardrails
			}
		}
		for (unsigned e = 0; e < 2; ++e) { // two sides
			float const ndv(bcube.d[!d][e]);
			pts[0][!d] = pts[1][!d] = (e ? ndv-w_expand : ndv);
			pts[2][!d] = pts[3][!d] = (e ? ndv : ndv+w_expand);

			for (unsigned f = 0; f < 2; ++f) { // top/bottom surfaces
				float const dz(f ? 0.0 : thickness);
				pts[0].z = pts[3].z = zvals[n] + dz;
				pts[1].z = pts[2].z = zval + dz;
				// construct previous and next pts in order to calculate the face normals below
				point prev_pt(pts[0]); prev_pt[d] -= delta_d; if (n   > 0       ) {prev_pt.z = zvals[n-1] + dz;}
				point next_pt(pts[1]); next_pt[d] += delta_d; if (n+1 < num_segs) {next_pt.z = zvals[n+2] + dz;}
				vector3d const v_shared((pts[2] - pts[1])*((f^d) ? -1.0 : 1.0));
				vector3d const fn(cross_product((pts[1] - pts[0]), v_shared));

				if (dot_product((cpos - pts[0]), fn) > 0) { // BFC
					vector3d const fn_prev(cross_product((pts[0] - prev_pt), v_shared)), fn_next(cross_product((next_pt - pts[1]), v_shared));
					vector3d const vn_prev((fn + fn_prev).get_norm()), vn_next((fn + fn_next).get_norm()); // average of both face normals
					vector3d const normals[4] = {vn_prev, vn_next, vn_next, vn_prev};
					qbd_bridge.add_quad_pts_vert_norms(pts, normals, cw_main);
				}
			}
			for (unsigned f = 0; f < 2; ++f) { // side surfaces
				unsigned const i0(f ? 3 : 0), i1(f ? 2 : 1);
				point pts2[4] = {pts[i0], pts[i1], pts[i1], pts[i0]};
				pts2[2].z += thickness; pts2[3].z += thickness; // top surface
				add_bridge_quad(pts2, cw_main, ((f^d) ? -1.0 : 1.0));
			}
			if (zval > cur_zpos) { // vertical cables
				conn_pts[e] = 0.5*(pts[1] + pts[2]);
				conn_pts[e].z += 0.5*thickness;

				if (!shadow_only && dist_val < 0.1) { // use high detail vertical cylinders
					s.set_cur_color(cables_color);
					point const &p(conn_pts[e]);
					draw_fast_cylinder(point(p.x, p.y, cur_zpos+0.15*thickness), point(p.x, p.y, zval), cable_thick, cable_thick, cable_ndiv, 0); // no ends
				}
				else { // use lower detail cubes; okay for shadow pass
					cube_t c(conn_pts[e]);
					c.z1() = cur_zpos;
					c.z2() = zval;
					c.expand_by(cable_thick);
					draw_cube(qbd_bridge, c, cw_cables, 1); // skip_bottom=1
				}
			}
		} // for e
		if (zval > cur_zpos + 0.5*scale) { // add horizontal connectors if high enough
			cube_t c(conn_pts[0], conn_pts[1]);
			vector3d exp(zero_vector);
			exp.z  = 0.75*conn_thick;
			exp[d] = conn_thick;
			c.expand_by(exp);
			draw_cube(qbd_bridge, c, cw_main, 0); // skip_bottom=0
		}
		cur_dval  = next_dval;
		cur_zpos += delta_z;
	} // for s
	// bottom surface
	float tscale(0.0);

	if (!shadow_only) {
		qbd_bridge.draw_and_clear(); // flush before texture change
		select_texture(get_texture_by_name("roads/asphalt.jpg"));
		tscale = 1.0/scale; // scale texture to match road width
	}
	bool const dir(bridge.slope), invert_normals((d!=0) ^ dir);
	float const dz_scale(bridge.dz()/(bridge.get_sz_dim(d)));
	cur_dval = bcube.d[d][0]; // reset to start
	point pts[8];

	for (unsigned n = 0; n <= sm_split_pos.size(); ++n) {
		float const next_dval((n == sm_split_pos.size()) ? bcube.d[d][1] : sm_split_pos[n]);
		cube_t bot_bc(bcube);
		bot_bc.d[ d][0]  = cur_dval; // one tile slice
		bot_bc.d[ d][1]  = next_dval;
		bot_bc.d[!d][0] += 0.4*w_expand;
		bot_bc.d[!d][1] -= 0.4*w_expand;
		float extend_dz1(dz_scale*(bridge.d[d][0] - bot_bc.d[d][0])), extend_dz2(dz_scale*(bot_bc.d[d][1] - bridge.d[d][1]));
		if (dir) {swap(extend_dz1, extend_dz2);}
		float const z1(bridge.z1() - extend_dz1 - 0.25*ROAD_HEIGHT), z2(bridge.z2() + extend_dz2 - 0.25*ROAD_HEIGHT); // move slightly downward
		point bot_center(bot_bc.get_cube_center());
		bot_center.z = 0.5f*(z1 + z2) - 0.5f*wall_width;

		if (!sm_split_pos.empty()) { // multiple tiles, must select a new shadow map set
			query_pt[d] = cur_dval;
			begin_tile(query_pt, 1); // will_emit_now=1
		}
		// add bottom road/bridge surface
		set_cube_pts(bot_bc, z1-wall_width, z2-wall_width, z1, z2, (d != 0), dir, pts);
		draw_cube(qbd_bridge, cw_concrete, bot_center, pts, 0, invert_normals, tscale); // skip_bottom=0

		// add guardrails/walls
		for (unsigned e = 0; e < 2; ++e) { // two sides
			cube_t side_bc(bot_bc);
			side_bc.d[!d][!e] = side_bc.d[!d][e] + (e ? -wall_width : wall_width);
			point side_center(side_bc.get_cube_center());
			side_center.z = 0.5f*(z1 + z2) + 0.5f*wall_height;
			set_cube_pts(side_bc, z1, z2, z1+wall_height, z2+wall_height, (d != 0), dir, pts);
			draw_cube(qbd_bridge, cw_concrete, side_center, pts, 1, invert_normals, tscale); // skip_bottom=1
		}
		qbd_bridge.draw_and_clear(); // flush
		cur_dval = next_dval;
	} // for n
}

void road_draw_state_t::add_bridge_quad(point const pts[4], color_wrapper const &cw, float normal_scale) {
	vector3d const normal(cross_product((pts[1] - pts[0]), (pts[3] - pts[0]))*normal_scale);
	if (dot_product(((camera_pdu.pos - xlate) - pts[0]), normal) < 0) return; // BFC
	qbd_bridge.add_quad_pts(pts, cw, normal.get_norm());
}

void road_draw_state_t::draw_tunnel(tunnel_t const &tunnel, bool shadow_only) { // Note: called rarely, so doesn't need to be efficient
	cube_t const bcube(tunnel.get_tunnel_bcube());
	if (!check_cube_visible(bcube, 1.0)) return; // VFC/too far
	ensure_shader_active(); // needed for use_smap=0 case
	quad_batch_draw &qbd(qbd_bridge); // use same qbd as bridges
	color_wrapper cw_concrete(LT_GRAY);
	bool const d(tunnel.dim), invert_normals(d);
	float const scale(1.0*city_params.road_width), wall_thick(TUNNEL_WALL_THICK*city_params.road_width), width(max(0.5*tunnel.get_width(), 2.0*(d ? DX_VAL : DY_VAL)));
	float zf, zb, end_ext;
	tunnel.calc_zvals_and_eext(zf, zb, end_ext);
	cube_t cubes[4]; // {bottom, top, left, right}
	tunnel.calc_top_bot_side_cubes(cubes);
	float const tile_sz(d ? MESH_Y_SIZE*DY_VAL : MESH_X_SIZE*DX_VAL), xy1(cubes[0].d[d][0]), xy2(cubes[0].d[d][1]), length(xy2 - xy1);
	unsigned const num_segs(ceil(length/tile_sz));
	float const dt(1.0/num_segs);
	point pts[8];
	float tscale(0.0);

	if (!shadow_only) {
		s.add_uniform_float("hemi_lighting_scale", 0.1); // mostly disable hemispherical lighting, which doesn't work for tunnel interiors
		select_texture(get_texture_by_name("roads/asphalt.jpg"));
		tscale = 1.0/scale; // scale texture to match road width
	}
	for (unsigned s = 0; s < num_segs; ++s) { // split into segments, one per tile shadow map
		float const t1(s*dt), t2((s+1)*dt), tmid(0.5f*(t1 + t2)); // in range [0.0, 1.0]
		point center(tunnel.get_cube_center());
		center[d] = xy1 + tmid*length; // halfway between the end points
		begin_tile(center, 1);

		for (unsigned i = 0; i < 4; ++i) { // add tunnel top, bottom, and sides
			cube_t c(cubes[i]); // deep copy so we can modify it
			c.d[d][0] = xy1 + t1*length;
			c.d[d][1] = xy1 + t2*length;
			float const dz(zb - zf), zft(zf + t1*dz), zbt(zf + t2*dz);
			point center(c.get_cube_center());
			center.z += 0.5f*(zft + zbt);
			// Note: could use draw_cylindrical_section() or draw_circle_normal() for cylindrical tunnel
			set_cube_pts(c, c.z1()+zft, c.z1()+zbt, c.z2()+zft, c.z2()+zbt, d, 0, pts); // dir=0 here
			draw_cube(qbd, cw_concrete, center, pts, (!shadow_only && i != 1), invert_normals, tscale); // skip_bottom=1 for all but the top cube unless shadowed
		} // for i
		qbd.draw_and_clear();
	} // for s
	if (!shadow_only) {
		s.add_uniform_float("hemi_lighting_scale", 0.5); // set back to the default of 0.5
		select_texture(get_texture_by_name("cblock2.jpg"));
		select_texture(get_texture_by_name("normal_maps/cblock2_NRM.jpg", 1), 5); // set normal map
		tscale *= 4.0;
	}
	for (unsigned n = 0; n < 2; ++n) { // add tunnel facades
		cube_t const &tend(tunnel.ends[n]);
		cube_t c(tend);
		c.z1() += tunnel.height - 0.5*wall_thick; // tunnel ceiling
		c.z2() += tunnel.facade_height[n]; // high enough to cover the hole in the mesh
		c.d[d][0] -= 0.9*end_ext; c.d[d][1] += 0.9*end_ext; // extend to cover the gaps in the mesh (both dirs) - slightly less than interior so that it sticks out
		begin_tile(c.get_cube_center(), 1); // required for long tunnels where facade is in a different tile from the tunnel center
		draw_cube(qbd, c, cw_concrete, 0, tscale); // skip_bottom=0
		c.z1() = tunnel.ends[n].z1() - 2.0*wall_thick; // extend below the bottom to cover any gaps at the corners
		c.d[!d][0] = tunnel.ends[n].d[!d][0] - width; c.d[!d][1] = tend.d[!d][0]; // left side
		draw_cube(qbd, c, cw_concrete, 1, tscale); // skip_bottom=1
		c.d[!d][1] = tunnel.ends[n].d[!d][1] + width; c.d[!d][0] = tend.d[!d][1]; // right side
		draw_cube(qbd, c, cw_concrete, 0, tscale); // skip_bottom=0 in case there are overhangs due to steep cliffs
		qbd.draw_and_clear();
	} // for n
	if (!shadow_only) {bind_default_flat_normal_map();} // restore flat normal map
}

void road_draw_state_t::draw_stoplights_and_street_signs(vector<road_isec_t> const &isecs, vector<road_t> const &roads, range_pair_t const &rp, unsigned cur_city, bool shadow_only) {
	for (unsigned i = rp.s; i < rp.e; ++i) {isecs[i].draw_stoplights_and_street_signs(*this, roads, cur_city, shadow_only);}
}

// not really related to roads, but I guess this goes here
void road_draw_state_t::draw_transmission_line_wires(point const &p1, point const &p2, point const wire_pts1[3], point const wire_pts2[3], float radius) {
	if (shadow_only) return; // not shadow casters
	if (p1 == p2)    return; // zero length wire segment?
	point const pts[6] = {p1+wire_pts1[0], p2+wire_pts2[0], p1+wire_pts1[1], p2+wire_pts2[1], p1+wire_pts1[2], p2+wire_pts2[2]};
	cube_t wires_bcube(pts, 6);
	if (!check_cube_visible(wires_bcube, 0.2)) return; // VFC/too far
	s.set_cur_color(BLACK);
	unsigned const ndiv = 4;
	for (unsigned n = 0; n < 3; ++n) {draw_fast_cylinder(pts[2*n], pts[2*n+1], radius, radius, ndiv, 0);} // no ends
}
void road_draw_state_t::draw_transmission_line(transmission_line_t const &tline) {
	float const wire_radius(0.0024*city_params.road_width), tower_radius(0.035*city_params.road_width);
	float const tower_bar_len(0.25*city_params.road_width), tower_bar_radius(0.4*tower_radius), end_rscale(0.5), bar_extend(1.04);
	float const standoff_radius(0.12*tower_radius), standoff_height(1.2*tower_radius), standoff_dmax(0.05*draw_tile_dist);
	float const standoff_len(standoff_height - bar_extend*end_rscale*tower_bar_radius); // actual height
	point const camera_bs(camera_pdu.pos - xlate);
	vector3d const wire_sep_dir(cross_product((tline.p2 - tline.p1), plus_z).get_norm()); // should be a straight line through all towers
	point wire_pts[3]; // relative to the top of the tower; this is the bottom of the standoffs

	for (unsigned n = 0; n < 3; ++n) {
		wire_pts[n].z = -((n + 1)*0.07*tline.tower_height); // shifted down
		wire_pts[n]  += (((n&1) ? -1.0 : 1.0)*tower_bar_len)*wire_sep_dir; // shift away from the tower
	}
	ensure_shader_active(); // is this necessary?
	point cur_pt; // starts as zero offset
	uint64_t last_tile_id(0);

	for (auto const &p : tline.tower_pts) {
		point bot_pt(p - vector3d(0.0, 0.0, 1.25*tline.tower_height)); // extend 25% further into the ground in case it's on a steep slope
		cube_t tower_bcube(p, bot_pt);
		tower_bcube.expand_by_xy(bar_extend*tower_bar_len);

		if (check_cube_visible(tower_bcube, 0.5)) {
			if (!shadow_only && use_smap) {
				uint64_t const tile_id(get_tile_id_containing_point_no_xyoff(p));

				if (last_tile_id == 0 || tile_id != last_tile_id) { // moved to a new tile
					try_bind_tile_smap_at_point((p + xlate), s);
					last_tile_id = tile_id;
				}
			}
			if (!shadow_only) {s.set_cur_color(GRAY);}
			bool const draw_top(camera_bs.z > tower_bcube.z2());
			draw_fast_cylinder(bot_pt, p, tower_radius, end_rscale*tower_radius, 20, 0, (draw_top ? 4 : 0)); // draw sides and maybe top
			tower_bcube.z1() = p.z + wire_pts[2].z - tower_bar_radius; // bottom of lowest bar

			if (check_cube_visible(tower_bcube, 0.3)) {
				for (unsigned n = 0; n < 3; ++n) { // draw 3 horizontal bars to support each wire, top to bottom
					point const bar_start(p + vector3d(0.0, 0.0, wire_pts[n].z+standoff_height)); // top of standoffs
					point const bar_end  (bar_start + bar_extend*vector3d(wire_pts[n].x, wire_pts[n].y, 0.0)); // extend a bit further out
					draw_fast_cylinder(bar_start, bar_end, tower_bar_radius, end_rscale*tower_bar_radius, 12, 0, 4); // draw sides and end
				}
				if (!shadow_only && tower_bcube.closest_dist_less_than(camera_bs, standoff_dmax)) {
					s.set_cur_color(LT_GRAY);

					for (unsigned n = 0; n < 3; ++n) { // draw 3 standoffs
						point const bot(p + wire_pts[n]);

						if (tower_bcube.closest_dist_less_than(camera_bs, 0.4*standoff_dmax)) { // draw as a stack of cones
							unsigned const num_segs = 8;
							float const len_per_seg(standoff_len/num_segs), hlen_per_seg(0.5*len_per_seg);
							float const r1(1.25*standoff_radius), r2(0.2*standoff_radius);

							for (unsigned n = 0; n < num_segs; ++n) {
								point pb(bot), pm(bot), pt(bot); // bottom, middle, top
								pb.z += n*len_per_seg;
								pm.z  = pb.z + hlen_per_seg;
								pt.z  = pm.z + hlen_per_seg;
								draw_fast_cylinder(pm, pb, r1, r2, 16, 0, 0); // truncated cone with no end
								draw_fast_cylinder(pm, pt, r1, r2, 16, 0, 0);
							} // for n
						}
						else { // draw as a single cylinder
							point const top(bot + vector3d(0.0, 0.0, standoff_len));
							bool const draw_top(camera_bs.z > 0.5f*(bot.z + top.z));
							draw_fast_cylinder(bot, top, standoff_radius, standoff_radius, 16, 0, (draw_top ? 4 : 3)); // draw sides and one end
						}
					}
				}
			}
		}
		draw_transmission_line_wires(cur_pt, p, ((cur_pt == all_zeros) ? tline.p1_wire_pts : wire_pts), wire_pts, wire_radius);
		cur_pt = p;
	} // for p
	draw_transmission_line_wires(cur_pt, all_zeros, wire_pts, tline.p2_wire_pts, wire_radius); // final segment
}

