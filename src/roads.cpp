// 3D World - Roads for Procedural Cities
// by Frank Gennari
// 11/20/18
#include "city.h"
#include "lightmap.h"

float const STREETLIGHT_BEAMWIDTH       = 0.25;
float const SLIGHT_DIST_TO_CORNER_SCALE = 3.0; // larger is closer to the road surface

extern bool tt_fire_button_down;
extern int frame_counter, game_mode, display_mode;
extern float FAR_CLIP;
extern vector<light_source> dl_sources;
extern city_params_t city_params;


void road_mat_mgr_t::ensure_road_textures() {
	if (inited) return;
	string const img_names[NUM_RD_TIDS] = {"sidewalk.jpg", "straight_road.jpg", "bend_90.jpg", "int_3_way.jpg", "int_4_way.jpg",
		                                   "parking_lot.png", "rail_tracks.jpg", "grass_park.jpg", "asphalt.jpg"};
	float const aniso[NUM_RD_TIDS] = {4.0, 16.0, 8.0, 8.0, 8.0, 4.0, 16.0};
	for (unsigned i = 0; i < NUM_RD_TIDS; ++i) {tids[i] = get_texture_by_name(("roads/" + img_names[i]), 0, 0, 1, aniso[i]);}
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

road_t::road_t(point const &s, point const &e, float width, bool dim_, bool slope_, unsigned road_ix_) : road_ix(road_ix_), dim(dim_), slope(slope_) {
	assert(s != e);
	assert(width > 0.0);
	vector3d const dw(0.5*width*cross_product((e - s), plus_z).get_norm());
	point const pts[4] = {(s - dw), (s + dw), (e + dw), (e - dw)};
	set_from_points(pts, 4);
}

void road_t::add_road_quad(quad_batch_draw &qbd, colorRGBA const &color, float ar) const { // specialized here for sloped roads (road segments and railroad tracks)
	if (z1() == z2()) {add_flat_city_quad(*this, qbd, color, ar); return;}
	bool const s(slope ^ dim);
	point pts[4] = {point(x1(), y1(), d[2][!s]), point(x2(), y1(), d[2][!s]), point(x2(), y2(), d[2][ s]), point(x1(), y2(), d[2][ s])};
	if (!dim) {swap(pts[0].z, pts[2].z);}
	vector3d const normal(cross_product((pts[2] - pts[1]), (pts[0] - pts[1])).get_norm());
	qbd.add_quad_pts(pts, color, normal, get_tex_range(ar));
}


tex_range_t parking_lot_t::get_tex_range(float ar) const { // ar is unused
	bool const d(!dim); // Note: R90
	float const xscale(1.0/(2.0*PARK_SPACE_WIDTH *city_params.get_nom_car_size().y));
	float const yscale(1.0/(1.0*PARK_SPACE_LENGTH*city_params.get_nom_car_size().x));
	float const tx(0.24), ty(0.0); // x=cols, y=rows
	return tex_range_t(tx, ty, (xscale*(d ? dy() : dx()) + tx), (yscale*(d ? dx() : dy()) + ty), 0, d);
}

tex_range_t driveway_t::get_tex_range(float ar) const { // ar is unused
	float txy[2] = {2.0, 2.0};
	txy[dim] *= get_sz_dim(dim)/get_sz_dim(!dim); // ensure 1:1 texture scale
	return tex_range_t(0.0, 0.0, txy[0], txy[1], 0, 0); // since the asphalt driveway texture isn't directional, we don't need to worry about swap_xy
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
			// however, it can lead to a light transitioning from green => yellow => green without going through red; this is difficult to fix with the current architecture
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

	void streetlight_t::draw(draw_state_t &dstate, bool shadow_only, bool is_local_shadow, bool always_on) const { // Note: translate has already been applied as a transform
		float const height(get_streetlight_height());
		point const center(pos + dstate.xlate + vector3d(0.0, 0.0, 0.5*height));
		if (shadow_only && is_local_shadow && !dist_less_than(camera_pdu.pos, center, 0.8*camera_pdu.far_)) return;
		if (!camera_pdu.sphere_visible_test(center, height)) return; // VFC
		float dist_val(0.0);

		if (shadow_only) {dist_val = (is_local_shadow ? 0.12 : 0.06);}
		else {
			dist_val = p2p_dist(camera_pdu.pos, center)/get_draw_tile_dist();
			if (dist_val > 0.2) return; // too far
		}
		float const pradius(get_streetlight_pole_radius()), lradius(light_radius*city_params.road_width);
		int const ndiv(shadow_only ? (is_local_shadow ? 4 : 8) : max(4, min(N_SPHERE_DIV, int(0.5/dist_val))));
		point const top(pos + vector3d(0.0, 0.0, 0.96*height)), lpos(get_lpos()), arm_end(lpos + vector3d(0.0, 0.0, 0.025*height) - 0.06*height*dir);
		if (!shadow_only) {dstate.s.set_cur_color(pole_color);}
		draw_fast_cylinder(pos, pos+vector3d(0.0, 0.0, height), pradius, 0.7*pradius, min(ndiv, 24), 0, 0); // vertical post, untextured, no ends
		if (dist_val <= 0.12) {draw_fast_cylinder(top, arm_end, 0.5*pradius, 0.4*pradius, min(ndiv, 16), 0, 0);} // untextured, no ends
		if (shadow_only && is_local_shadow) return; // top part never projects a shadow on a visible object near the ground
		bool const is_on(is_lit(always_on));

		if (!shadow_only) {
			if (!is_on && dist_val > 0.15) return; // too far
			if (is_on) {dstate.s.set_color_e(light_color);} else {dstate.s.set_cur_color(light_color);} // emissive when lit
			// streetlight bloom: should be elliptical, disable when viewed from above, use tighter alpha mask
			//if (is_on && dist_val > 0.05) {dstate.add_light_flare((lpos - vector3d(0.0, 0.0, 0.6*lradius)), zero_vector, light_color, min(1.0f, 10.0f*dist_val), 2.5*lradius);} // non-directional
		}
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
		dl_sources.emplace_back(ldist, lpos, lpos, light_color, 0, -plus_z, STREETLIGHT_BEAMWIDTH); // points down
		
		// cache shadow maps if there are no dynamic cars or pedestrians (player doesn't cast a shadow)
		if (city_params.num_cars == 0 && city_params.num_peds == 0) {
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


void streetlights_t::draw_streetlights(draw_state_t &dstate, bool shadow_only, bool always_on) const {
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

bool streetlights_t::check_streetlight_sphere_coll_xy(point const &center, float radius) const {
	float const pradius(streetlight_ns::get_streetlight_pole_radius()), rtot(pradius + radius), y_end(center.y + rtot);

	for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {
		if (i->pos.y > y_end) break; // streetlights are sorted by y-val; if this holds, we can no longer intersect
		if (dist_xy_less_than(i->pos, center, rtot)) return 1;
	}
	return 0;
}

bool streetlights_t::line_intersect_streetlights(point const &p1, point const &p2, float &t) const {
	bool ret(0);
	for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
	return ret;
}


road_isec_t::road_isec_t(cube_t const &c, int rx, int ry, unsigned char conn_, bool at_conn_road, short conn_to_city_) :
	cube_t(c), num_conn(0), conn(conn_), conn_to_city(conn_to_city_), stoplight(at_conn_road)
{
	rix_xy[0] = rix_xy[1] = rx; rix_xy[2] = rix_xy[3] = ry; conn_ix[0] = conn_ix[1] = conn_ix[2] = conn_ix[3] = 0;
	if (conn == 15) {num_conn = 4;} // 4-way
	else if (conn == 7 || conn == 11 || conn == 13 || conn == 14) {num_conn = 3;} // 3-way
	else if (conn == 5 || conn == 6  || conn == 9  || conn == 10) {num_conn = 2;} // 2-way
	else {assert(0);}
	stoplight.init(num_conn, conn);
}

tex_range_t road_isec_t::get_tex_range(float ar) const {
	switch (conn) {
	case 5 : return tex_range_t(0.0, 0.0, -1.0,  1.0, 0, 0); // 2-way: MX
	case 6 : return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 2-way: R0
	case 9 : return tex_range_t(0.0, 0.0, -1.0, -1.0, 0, 0); // 2-way: MXMY
	case 10: return tex_range_t(0.0, 0.0,  1.0, -1.0, 0, 0); // 2-way: MY
	case 7 : return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 3-way: R0
	case 11: return tex_range_t(0.0, 0.0, -1.0, -1.0, 0, 0); // 3-way: MY
	case 13: return tex_range_t(0.0, 0.0,  1.0, -1.0, 0, 1); // 3-way: R90MY
	case 14: return tex_range_t(0.0, 0.0, -1.0,  1.0, 0, 1); // 3-way: R90MX
	case 15: return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 4-way: R0
	default: assert(0);
	}
	return tex_range_t(0.0, 0.0, 1.0, 1.0); // never gets here
}

void road_isec_t::make_4way(unsigned conn_to_city_) {
	num_conn = 4; conn = 15;
	assert(conn_to_city < 0); conn_to_city = conn_to_city_;
	stoplight.init(num_conn, conn); // re-init with new conn
}

bool road_isec_t::can_go_based_on_light(car_base_t const &car) const {
	return (!red_or_yellow_light(car) || stoplight.can_turn_right_on_red(car)); // green light or right turn on red
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
	assert(conn & (1<<new_orient)); // car mustto go an enabled orient
	return new_orient;
}

bool road_isec_t::can_go_now(car_t const &car) const {
	if (!can_go_based_on_light(car)) return 0;
	if (stoplight.check_int_clear(car)) return 1;
	if ((frame_counter&15) == 0) {car.honk_horn_if_close_and_fast();} // honk every so often
	return 0; // intersection not clear
}

bool road_isec_t::has_left_turn_signal(unsigned orient) const {
	if (num_conn == 2) return 0; // never
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

bool road_isec_t::check_sphere_coll(point const &pos, float radius) const { // used for peds
	if (num_conn == 2) return 0; // no stoplights
	if (!sphere_cube_intersect_xy(pos, radius, *this)) return 0;

	for (unsigned n = 0; n < 4; ++n) {
		if (!(conn & (1 << n))) continue; // no road in this dir
		if (sphere_cube_intersect(pos, radius, get_stoplight_cube(n))) return 1;
	}
	return 0;
}

bool road_isec_t::proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, float dist, vector3d *cnorm) const {
	if (num_conn == 2) return 0; // no stoplights
	if (!sphere_cube_intersect_xy(pos, (radius + dist), (*this + xlate))) return 0;

	for (unsigned n = 0; n < 4; ++n) {
		if (!(conn & (1<<n))) continue; // no road in this dir
		if (sphere_cube_int_update_pos(pos, radius, (get_stoplight_cube(n) + xlate), p_last, 1, 0, cnorm)) return 1; // typically only one coll, just return the first one
	}
	return 0;
}

bool check_line_clip_update_t(point const &p1, point const &p2, float &t, cube_t const &c) {
	float tmin(0.0), tmax(1.0);
	if (get_line_clip(p1, p2, c.d, tmin, tmax) && tmin < t) {t = tmin; return 1;}
	return 0;
}

bool road_isec_t::line_intersect(point const &p1, point const &p2, float &t) const {
	if (num_conn == 2) return 0; // no stoplights
	cube_t c(*this); // deep copy
	c.z2() += stoplight_ns::stoplight_max_height();
	if (!c.line_intersects(p1, p2)) return 0;
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

void road_isec_t::draw_stoplights(quad_batch_draw &qbd, draw_state_t &dstate, bool shadow_only) const {
	if (num_conn == 2) return; // no stoplights
	cube_t sl_bcube(*this);
	sl_bcube.z2() += 0.276*city_params.road_width; // add max stoplight height
	if (!dstate.check_cube_visible(sl_bcube, 0.16, shadow_only)) return; // dist_scale=0.16
	point const center(get_cube_center() + dstate.xlate);
	float const dist_val(shadow_only ? 0.0 : p2p_dist(camera_pdu.pos, center)/get_draw_tile_dist());
	vector3d const cview_dir(camera_pdu.pos - center);
	float const sz(0.03*city_params.road_width), h(1.0*sz);
	color_wrapper cw(BLACK);

	for (unsigned n = 0; n < 4; ++n) { // {-x, +x, -y, +y} = {W, E, S, N} facing = car traveling {E, W, N, S}
		if (!(conn & (1<<n))) continue; // no road in this dir
		bool const dim((n>>1) != 0), dir((n&1) == 0), side((dir^dim^1) != 0); // Note: dir is inverted here to represent car dir
		float const zbot(z1() + 2.0*h), dim_pos(d[dim][!dir] + SLIGHT_DIST_TO_CORNER_SCALE*(dir ? sz : -sz)); // pos in road dim
		float const side_len(side ? -sz : sz), v1(d[!dim][side] + (SLIGHT_DIST_TO_CORNER_SCALE - 1.0)*side_len), v2(v1 + side_len); // pos in other dim
		// draw base
		unsigned const num_segs(has_left_turn_signal(n) ? 6 : 3);
		float const sl_top(zbot + 1.2*h*num_segs), sl_lo(min(v1, v2) - 0.25*sz), sl_hi(max(v1, v2) + 0.25*sz);

		if (dist_val > 0.06) { // draw front face only
			point pts[4];
			pts[0][dim]  = pts[1][dim]  = pts[2][dim] = pts[3][dim] = dim_pos;
			pts[0][!dim] = pts[3][!dim] = sl_lo;
			pts[1][!dim] = pts[2][!dim] = sl_hi;
			pts[0].z = pts[1].z = z1();
			pts[2].z = pts[3].z = sl_top;
			qbd.add_quad_pts(pts, cw,  (dim ? (dir ? plus_y : -plus_y) : (dir ? plus_x : -plus_x))); // Note: normal doesn't really matter since color is black
		}
		else {
			cube_t c;
			c.z1() = z1(); c.z2() = sl_top;
			c.d[ dim][0] = dim_pos - (dir ? -0.04 : 0.5)*sz; c.d[dim][1] = dim_pos + (dir ? 0.5 : -0.04)*sz;
			c.d[!dim][0] = sl_lo; c.d[!dim][1] = sl_hi;
			dstate.draw_cube(qbd, c, cw, 1); // skip_bottom=1; Note: uses traffic light texture, but color is black so it's all black anyway

			if (!shadow_only && tt_fire_button_down && game_mode != 1) {
				point const p1(camera_pdu.pos - dstate.xlate), p2(p1 + camera_pdu.dir*FAR_CLIP);
				if (c.line_intersects(p1, p2)) {dstate.set_label_text(stoplight.label_str(), (c.get_cube_center() + dstate.xlate));}
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
				dstate.draw_cube(qbd, c, cw, shadow_only); // skip_bottom=shadow_only

				if (!shadow_only && cw_color != BLACK) { // draw light
					point p[4];
					p[0][!dim] = p[1][!dim] = p[2][!dim] = p[3][!dim] = c.d[!dim][!i] + 0.01*(i ? -sz : sz);
					p[0][dim] = p[3][dim] = c.d[dim][0]; p[1][dim] = p[2][dim] = c.d[dim][1];
					p[0].z = p[1].z = c.z1(); p[2].z = p[3].z = c.z2();
					vector3d normal(zero_vector);
					normal[!dim] = (i ? -1.0 : 1.0);
					point const cw_center(0.25*(p[0] + p[1] + p[2] + p[3]));
					vector3d const cw_cview_dir(camera_pdu.pos - (cw_center + dstate.xlate));
					float dp(dot_product(normal, cw_cview_dir));

					if (dp > 0.0) { // if back facing, don't draw the lights
						float const mag(CLIP_TO_01(2.5f*dp/cw_cview_dir.mag() - 1.0f)); // normalize and strengthen slope

						if (mag > 0.0) {
							bool const is_walk(cw_state == stoplight_ns::CW_WALK);
							qbd.add_quad_pts(p, WHITE*mag, normal, tex_range_t((is_walk ? 0.25 : 0.0), 0.0, (is_walk ? 0.5 : 0.25), 0.5)); // color is stored in the texture
							if (draw_detail) {dstate.add_light_flare(cw_center, normal, cw_color*mag, 1.0, 0.8*h);}
						}
					}
				}
			} // for i
		} // end crosswalk signal
		if (shadow_only)    continue; // no lights in shadow pass
		if (dist_val > 0.1) continue; // too far away
		vector3d normal(zero_vector);
		normal[dim] = (dir ? -1.0 : 1.0);
		if (dot_product(normal, cview_dir) < 0.0) continue; // if back facing, don't draw the lights
		// draw straight/line turn light
		point p[4];
		p[0][dim] = p[1][dim] = p[2][dim] = p[3][dim] = dim_pos;
		p[0][!dim] = p[3][!dim] = v1; p[1][!dim] = p[2][!dim] = v2;
		p[0].z = p[1].z = zbot; p[2].z = p[3].z = zbot + h;
		draw_sl_block(qbd, dstate, p, h, stoplight.get_light_state(dim, dir, TURN_NONE), draw_detail, draw_detail, normal, tex_range_t(0.0, 0.5, 0.5, 1.0));

		if (has_left_turn_signal(n)) { // draw left turn light (upper light section)
			draw_sl_block(qbd, dstate, p, h, stoplight.get_light_state(dim, dir, TURN_LEFT), draw_detail, 0.5*draw_detail, normal, tex_range_t(1.0, 0.5, 0.5, 1.0));
		}
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
	if (staggered) {num_per_side *= 2;}

	for (unsigned n = 0; n < num_per_side; ++n) {
		float const v((n + 0.5)/num_per_side); // 1/8, 3/8, 5/8, 7/8
		point pos1, pos2;
		pos1[ dim] = pos2[dim] = d[dim][0] + v*dsz;
		pos1[!dim] = d[!dim][0] - dn_shift; pos2[!dim] = d[!dim][1] + dn_shift;
		pos1.z = pos2.z = za + v*(zb - za);
		vector3d dir(zero_vector); dir[!dim] = 1.0;
		if (!staggered || (n&1) == 0) {streetlights.emplace_back(pos1,  dir);}
		if (!staggered || (n&1) == 1) {streetlights.emplace_back(pos2, -dir);}
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

void tunnel_t::calc_top_bot_side_cubes(cube_t cubes[4]) const {
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
	bool const prev_int(contains_pt_xy(prev - xlate));
	if (!prev_int && !contains_pt_xy(center - xlate)) return 0; // no x/y containment for road segment (cur or prev)
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

void tunnel_t::calc_zvals_and_eext(float &zf, float &zb, float &end_ext) const {
	zf = ends[0].z1(); zb = ends[1].z1();
	end_ext = (2.0*(dim ? DY_VAL : DX_VAL));
	float const dz_ext(end_ext*(zb - zf)/get_length());
	zf -= dz_ext; zb += dz_ext; // adjust zvals for extension
}

bool tunnel_t::line_intersect(point const &p1, point const &p2, float &t) const {
	cube_t const bcube(get_tunnel_bcube());
	if (!check_line_clip(p1, p2, bcube.d)) return 0;
	cube_t cubes[4];
	calc_top_bot_side_cubes(cubes);
	float const wall_thick(TUNNEL_WALL_THICK*city_params.road_width), width(max(0.5*get_width(), 2.0*(dim ? DX_VAL : DY_VAL)));
	float zf, zb, end_ext;
	calc_zvals_and_eext(zf, zb, end_ext);
	bool const d(dim);
	bool ret(0);

	for (unsigned i = 0; i < 4; ++i) { // check tunnel top, bottom, and sides
		//cube_t const &c(cubes[i]);
		//set_cube_pts(c, c.z1()+zf, c.z1()+zb, c.z2()+zf, c.z2()+zb, d, 0, pts); // TODO: tilted cube case
	}
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
	if (qbd_sl.empty()) return; // no stoplights to draw
	glDepthFunc(GL_LEQUAL); // helps prevent Z-fighting
	shader_t s;

	if (shadow_only) {s.begin_color_only_shader();}
	else {
		s.begin_simple_textured_shader(); // Note: no lighting
		road_mat_mgr.set_stoplight_texture();
	}
	qbd_sl.draw_and_clear();
	s.end_shader();
	glDepthFunc(GL_LESS);
}

void road_draw_state_t::draw_bridge(bridge_t const &bridge, bool shadow_only) { // Note: called rarely, so doesn't need to be efficient
	//timer_t timer("Draw Bridge"); // 0.065ms - 0.11ms
	cube_t bcube(bridge);
	float const scale(1.0*city_params.road_width);
	bcube.z2() += 2.0*scale; // make it higher
	unsigned const d(bridge.dim);
	float const l_expand(2.0*(d ? DY_VAL : DY_VAL)); // slight expand along road dim so that we're sure to cover the entire gap
	float const w_expand(0.25*scale); // expand width to add space for supports
	bcube.d[ d][0] -= l_expand; bcube.d[ d][1] += l_expand;
	bcube.d[!d][0] -= w_expand; bcube.d[!d][1] += w_expand;
	max_eq(bcube.d[d][0], bridge.src_road.d[d][0]); // clamp to orig road segment length
	min_eq(bcube.d[d][1], bridge.src_road.d[d][1]);
	if (!check_cube_visible(bcube, 1.0, shadow_only)) return; // VFC/too far
	point const cpos(camera_pdu.pos - xlate);
	float const center(bcube.get_center_dim(!d)), len(bcube.get_sz_dim(d));
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
	float const dist_val(shadow_only ? 1.0 : p2p_dist(camera_pdu.pos, closest_pt)/get_draw_tile_dist());
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
		zvals[n] = zpos + 0.3f*len*(1.0f - v*v) - 0.5f*scale;
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
					draw_fast_cylinder(point(p.x, p.y, cur_zpos), point(p.x, p.y, zval), cable_thick, cable_thick, cable_ndiv, 0); // no ends
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
	if (!check_cube_visible(bcube, 1.0, shadow_only)) return; // VFC/too far
	ensure_shader_active(); // needed for use_smap=0 case
	quad_batch_draw &qbd(qbd_bridge); // use same qbd as bridges
	color_wrapper cw_concrete(LT_GRAY);
	bool const d(tunnel.dim), invert_normals(d);
	float const scale(1.0*city_params.road_width), wall_thick(TUNNEL_WALL_THICK*city_params.road_width), width(max(0.5*tunnel.get_width(), 2.0*(d ? DX_VAL : DY_VAL)));
	float zf, zb, end_ext;
	tunnel.calc_zvals_and_eext(zf, zb, end_ext);
	cube_t cubes[4];
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
		select_multitex(get_texture_by_name("normal_maps/cblock2_NRM.jpg", 1), 5); // TODO: normal maps not enabled in the shader, this is future work
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
	if (!shadow_only) {select_multitex(FLAT_NMAP_TEX, 5);} // restore flat normal map
}

void road_draw_state_t::draw_stoplights(vector<road_isec_t> const &isecs, range_pair_t const &rp, bool shadow_only) {
	for (unsigned i = rp.s; i < rp.e; ++i) {isecs[i].draw_stoplights(qbd_sl, *this, shadow_only);}
}

