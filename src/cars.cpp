// 3D World - Cars for Procedural Cities
// by Frank Gennari
// 11/19/18
#include "city.h"
#include "file_utils.h"
#include "openal_wrap.h"
#include "explosion.h" // for add_blastr()
#include "lightmap.h" // for light_source
#include <cfloat> // for FLT_MAX

float const MIN_CAR_STOP_SEP = 0.25; // in units of car lengths

extern bool tt_fire_button_down;
extern int display_mode, game_mode, map_mode, animate2;
extern float FAR_CLIP;
extern point pre_smap_player_pos;
extern vector<light_source> dl_sources;
extern city_params_t city_params;


bool city_model_t::read(FILE *fp) { // filename body_material_id fixed_color_id xy_rot dz lod_mult shadow_mat_ids

	assert(fp);
	if (!read_string(fp, fn)) return 0;
	if (!read_int(fp, body_mat_id)) return 0;
	if (!read_int(fp, fixed_color_id)) return 0;
	if (!read_float(fp, xy_rot)) return 0;
	if (!read_float(fp, scale)) return 0;
	if (!read_float(fp, lod_mult) || lod_mult <= 0.0) return 0;
	unsigned shadow_mat_id;
	while (read_uint(fp, shadow_mat_id)) {shadow_mat_ids.push_back(shadow_mat_id);}
	return 1;
}

bool city_model_t::check_filename() {
	// if converting files, check if user specified the post-converted model3d filename rather than the input obj file
	if (city_params.convert_model_files && get_file_extension(fn, 0, 1) == "model3d" && !check_file_exists(fn)) {
		string const fn_as_obj(fn.substr(0, (fn.size() - 7)) + "obj"); // strip off 'model3d' and add 'obj'
		if (check_file_exists(fn_as_obj)) {fn = fn_as_obj; return 1;}
	}
	return check_file_exists(fn); // try to open model file for reading, but don't actually read anything; also, let the caller handle error printing
}


float car_t::get_max_lookahead_dist() const {return (get_length() + city_params.road_width);} // extend one car length + one road width in front
float car_t::get_turn_rot_z(float dist_to_turn) const {return (1.0 - CLIP_TO_01(4.0f*fabs(dist_to_turn)/city_params.road_width));}

bool car_t::headlights_on() const { // no headlights when parked
	return (!is_parked() && (in_tunnel || ((light_factor < (0.5 + HEADLIGHT_ON_RAND)) && is_night(HEADLIGHT_ON_RAND*signed_rand_hash(height + max_speed)))));
}

void car_t::apply_scale(float scale) {
	if (scale == 1.0) return; // no scale
	float const prev_height(height);
	height *= scale;
	point const pos(get_center());
	bcube.z2() += height - prev_height; // z1 is unchanged
	float const dx(bcube.x2() - pos.x), dy(bcube.y2() - pos.y);
	bcube.x1() = pos.x - scale*dx; bcube.x2() = pos.x + scale*dx;
	bcube.y1() = pos.y - scale*dy; bcube.y2() = pos.y + scale*dy;
}

void car_t::destroy() { // Note: not calling create_explosion(), so no chain reactions
	point const pos(get_center() + get_tiled_terrain_model_xlate());
	float const length(get_length());
	static rand_gen_t rgen;

	for (int n = 0; n < rgen.rand_int(3, 5); ++n) {
		vector3d off(rgen.signed_rand_vector_spherical()*(0.5*length));
		off.z = abs(off.z); // not into the ground
		point const exp_pos(pos + off);
		float const radius(rgen.rand_uniform(1.0, 1.5)*length), time(rgen.rand_uniform(0.3, 0.8));
		add_blastr(exp_pos, (exp_pos - get_camera_pos()), radius, 0.0, time*TICKS_PER_SECOND, CAMERA_ID, YELLOW, RED, ETYPE_ANIM_FIRE, nullptr, 1);
		gen_smoke(exp_pos, 1.0, rgen.rand_uniform(0.4, 0.6));
	} // for n
	gen_delayed_from_player_sound(SOUND_EXPLODE, pos, 1.0);
	park();
	destroyed = 1;
}

float car_t::get_min_sep_dist_to_car(car_t const &c, bool add_one_car_len) const {
	float const avg_len(0.5f*(get_length() + c.get_length())); // average length of the two cars
	float const min_speed(max(0.0f, (min(cur_speed, c.cur_speed) - 0.1f*max_speed))); // relative to max speed of 1.0, clamped to 10% at bottom end for stability
	return avg_len*(MIN_CAR_STOP_SEP + 1.11*min_speed + (add_one_car_len ? 1.0 : 0.0)); // 25% to 125% car length, depending on speed (2x on connector roads)
}

string car_t::str() const {
	std::ostringstream oss;
	oss << "Car " << TXT(dim) << TXT(dir) << TXT(cur_city) << TXT(cur_road) << TXT(cur_seg) << TXT(dz) << TXT(max_speed) << TXT(cur_speed)
		<< TXTi(cur_road_type) << TXTi(color_id) << " bcube=" << bcube.str();
	return oss.str();
}

string car_t::label_str() const {
	std::ostringstream oss;
	oss << TXT(dim) << TXTn(dir) << TXT(cur_city) << TXT(cur_road) << TXTn(cur_seg) << TXT(dz) << TXTn(turn_val) << TXT(max_speed) << TXTn(cur_speed)
		<< "wait_time=" << get_wait_time_secs() << "\n" << TXTin(cur_road_type)
		<< TXTn(stopped_at_light) << TXTn(in_isect()) << "cars_in_front=" << count_cars_in_front() << "\n" << TXT(dest_city) << TXTn(dest_isec);
	oss << "car=" << this << " car_in_front=" << car_in_front << endl; // debugging
	return oss.str();
}

void car_t::move(float speed_mult) {
	prev_bcube = bcube;
	if (destroyed || stopped_at_light || is_stopped()) return;
	assert(speed_mult >= 0.0 && cur_speed > 0.0 && cur_speed <= CONN_ROAD_SPEED_MULT*max_speed); // Note: must be valid for connector road => city transitions
	float dist(cur_speed*speed_mult);
	if (dz != 0.0) {dist *= min(1.25, max(0.75, (1.0 - 0.5*dz/get_length())));} // slightly faster down hills, slightly slower up hills
	min_eq(dist, 0.25f*city_params.road_width); // limit to half a car length to prevent cars from crossing an intersection in a single frame
	move_by(dir ? dist : -dist);
	// update waiting state
	float const cur_pos(bcube.d[dim][dir]);
	if (fabs(cur_pos - waiting_pos) > get_length()) {waiting_pos = cur_pos; reset_waiting();} // update when we move at least a car length
}

void car_t::maybe_accelerate(float mult) {
	if (car_in_front) {
		float const dist_sq(p2p_dist_xy_sq(get_center(), car_in_front->get_center())), length(get_length());

		if (dist_sq > length*length) { // if cars are colliding, let the collision detection system handle it
			float const dmin(get_min_sep_dist_to_car(*car_in_front, 1)); // add_one_car_len=1; space between the two car centers
			if (dist_sq < dmin*dmin) {decelerate(mult); return;} // too close to the car in front - decelerate instead
		}
	}
	accelerate(mult);
}

point car_base_t::get_front(float dval) const {
	point car_front(get_center());
	car_front[dim] += (dir ? dval : -dval)*get_length(); // half length
	return car_front;
}

bool car_t::front_intersects_car(car_t const &c) const {
	return (c.bcube.contains_pt(get_front(0.25)) || c.bcube.contains_pt(get_front(0.5))); // check front-middle and very front
}

void car_t::honk_horn_if_close() const {
	point const pos(get_center());
	if (dist_less_than((pos + get_tiled_terrain_model_xlate()), get_camera_pos(), 1.0)) {gen_sound(SOUND_HORN, pos);}
}

void car_t::honk_horn_if_close_and_fast() const {
	if (cur_speed > 0.25*max_speed) {honk_horn_if_close();}
}

void car_t::on_alternate_turn_dir(rand_gen_t &rgen) {
	honk_horn_if_close();
	if ((rgen.rand()&3) == 0) {dest_valid = 0;} // 25% chance of choosing a new destination rather than driving in circles; will be in current city
}

void car_t::register_adj_car(car_t &c) {
	if (car_in_front != nullptr && p2p_dist_xy_sq(get_center(), c.get_center()) > p2p_dist_xy_sq(get_center(), car_in_front->get_center())) return; // already found a closer car
	cube_t cube(bcube);
	cube.d[dim][!dir] = cube.d[dim][dir];
	cube.d[dim][dir] += (dir ? 1.0 : -1.0)*get_max_lookahead_dist();
	if (cube.intersects_xy(c.bcube)) {car_in_front = &c;} // projected cube intersects other car
}

unsigned car_t::count_cars_in_front(cube_t const &range) const { // Note: currently only used for debug printouts, so the iteration limit is acceptable
	unsigned num(0);
	car_t const *cur_car(this);

	for (unsigned i = 0; i < 50; ++i) { // limit iterations
		cur_car = cur_car->car_in_front;
		if (!cur_car || (!range.is_all_zeros() && !range.contains_pt_xy(cur_car->get_center()))) break;
		if (cur_car->dim != dim || cur_car->dir == dir) {++num;} // include if not going in opposite direction
	}
	return num;
}

float car_t::get_sum_len_space_for_cars_in_front(cube_t const &range) const {
	float len(0.0);
	car_t const *cur_car(this);

	// Note: should exit once we reach the end of the line of cars, or once we go off the current road segment;
	// this iteration may be very long for cars stopped on long, congested connector roads;
	// however, it should only be queried by one other car per frame, which means this is overall constant time per frame
	while (1) { // avg len = city_params.get_nom_car_size().x
		if (cur_car->dim != dim || cur_car->dir == dir) {len += cur_car->get_length();} // include if not going in opposite direction
		cur_car = cur_car->car_in_front;
		if (!cur_car || !range.contains_pt_xy(cur_car->get_center())) break;
	}
	return len * (1.0 + MIN_CAR_STOP_SEP); // car length + stopped space (including one extra space for the car behind us)
}

bool car_t::proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, vector3d *cnorm) const {
	return sphere_cube_int_update_pos(pos, radius, (bcube + xlate), p_last, 1, 0, cnorm); // Note: approximate when car is tilted or turning
	//return sphere_sphere_int((bcube.get_cube_center() + xlate), pos, bcube.get_bsphere_radius(), radius, cnorm, pos); // Note: handle cnorm in if using this
}

bool car_t::check_collision(car_t &c, road_gen_base_t const &road_gen) {

	if (c.dim != dim) { // turning in an intersection, etc. (Note: may not be needed, but at least need to return here)
		car_t *to_stop(nullptr);
		if (c.front_intersects_car(*this)) {to_stop = &c;}
		else if (front_intersects_car(c))  {to_stop = this;}
		if (!to_stop) return 0;
		to_stop->decelerate_fast(); // attempt to prevent one car from T-boning the other
		to_stop->bcube = to_stop->prev_bcube;
		to_stop->honk_horn_if_close_and_fast();
		return 1;
	}
	if (dir != c.dir) return 0; // traveling on opposite sides of the road
	float const sep_dist(get_min_sep_dist_to_car(c));
	float const test_dist(0.999*sep_dist); // slightly smaller than separation distance
	cube_t bcube_ext(bcube);
	bcube_ext.d[dim][0] -= test_dist; bcube_ext.d[dim][1] += test_dist; // expand by test_dist distance
	if (!bcube_ext.intersects_xy(c.bcube)) return 0;
	float const front(bcube.d[dim][dir]), c_front(c.bcube.d[dim][dir]);
	bool const move_c((front < c_front) ^ dir); // move the car that's behind
	// Note: we could slow the car in behind, but that won't work for initial placement collisions when speed == 0
	car_t &cmove(move_c ? c : *this); // the car that will be moved
	car_t const &cstay(move_c ? *this : c); // the car that won't be moved
	//cout << "Collision between " << cmove.str() << " and " << cstay.str() << endl;
	if (cstay.is_stopped()) {cmove.decelerate_fast();} else {cmove.decelerate();}
	float const dist(cstay.bcube.d[dim][!dir] - cmove.bcube.d[dim][dir]); // signed distance between the back of the car in front, and the front of the car in back
	point delta(all_zeros);
	delta[dim] += dist + (cmove.dir ? -sep_dist : sep_dist); // force separation between cars
	cube_t const &bcube(road_gen.get_bcube_for_car(cmove));
	if (cstay.max_speed < cmove.max_speed) {cmove.front_car_turn_dir = cstay.turn_dir;} // record the turn dir of this slow car in front of us so we can turn a different way

	if (!bcube.contains_cube_xy(cmove.bcube + delta)) { // moved outside its current road segment bcube
		//if (cmove.bcube == cmove.prev_bcube) {return 1;} // collided, but not safe to move the car (init pos or second collision)
		if (cmove.bcube != cmove.prev_bcube) { // try resetting to last frame's position
			cmove.bcube  = cmove.prev_bcube; // restore prev frame's pos
			//cmove.honk_horn_if_close_and_fast();
			return 1; // done
		}
		else { // keep the car from moving outside its current segment (init collision case)
			if (cmove.dir) {max_eq(delta[dim], min(0.0f, 0.999f*(bcube.d[cmove.dim][0] - cmove.bcube.d[cmove.dim][0])));}
			else           {min_eq(delta[dim], max(0.0f, 0.999f*(bcube.d[cmove.dim][1] - cmove.bcube.d[cmove.dim][1])));}
		}
	}
	cmove.bcube += delta;
	return 1;
}


bool comp_car_road_then_pos::operator()(car_t const &c1, car_t const &c2) const { // sort spatially for collision detection and drawing
	if (c1.cur_city != c2.cur_city) return (c1.cur_city < c2.cur_city);
	if (c1.is_parked() != c2.is_parked()) {return c2.is_parked();} // parked cars last
	if (c1.cur_road != c2.cur_road) return (c1.cur_road < c2.cur_road);

	if (c1.is_parked()) { // sort parked cars back to front relative to camera so that alpha blending works
		return (p2p_dist_xy_sq(c1.bcube.get_cube_center(), camera_pos) > p2p_dist_xy_sq(c2.bcube.get_cube_center(), camera_pos));
	}
	return (c1.bcube.d[c1.dim][c1.dir] < c2.bcube.d[c2.dim][c2.dir]); // compare front end of car (used for collisions)
}


unsigned car_model_loader_t::num_models() const {return city_params.car_model_files.size();}

city_model_t const &car_model_loader_t::get_model(unsigned id) const {
	assert(id < num_models());
	return city_params.car_model_files[id];
}

bool city_model_loader_t::is_model_valid(unsigned id) {
	assert(id < num_models());
	ensure_models_loaded(); // I guess we have to load the models here to determine if they're valid
	assert(id < models_valid.size());
	return (models_valid[id] != 0);
}

void city_model_loader_t::load_models() {
	models_valid.resize(num_models(), 1); // assume valid

	for (unsigned i = 0; i < num_models(); ++i) {
		string const &fn(get_model(i).fn);
		bool const recalc_normals = 1;

		if (!load_model_file(fn, *this, geom_xform_t(), -1, WHITE, 0, 0.0, recalc_normals, 0, city_params.convert_model_files, 1)) {
			cerr << "Error: Failed to read model file '" << fn << "'; Skipping this model (will use default low poly model)." << endl;
			push_back(model3d(fn, tmgr)); // add a placeholder dummy model
			models_valid[i] = 0;
		}
	} // for i
}

void city_model_loader_t::draw_model(shader_t &s, vector3d const &pos, cube_t const &obj_bcube, vector3d const &dir, colorRGBA const &color,
	point const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail, bool enable_animations)
{
	assert(is_model_valid(model_id));
	assert(size() == num_models()); // must be loaded
	city_model_t const &model_file(get_model(model_id));
	model3d &model(at(model_id));

	if (!is_shadow_pass && model_file.body_mat_id >= 0 && color.A != 0.0) { // use custom color for body material
		material_t &body_mat(model.get_material(model_file.body_mat_id));
		body_mat.ka = body_mat.kd = color;
	}
	model.bind_all_used_tids();
	cube_t const &bcube(model.get_bcube());
	point const orig_camera_pos(camera_pdu.pos);
	camera_pdu.pos += bcube.get_cube_center() - pos - xlate; // required for distance based LOD
	bool const camera_pdu_valid(camera_pdu.valid);
	camera_pdu.valid = 0; // disable VFC, since we're doing custom transforms here
	// Note: in model space, front-back=z, left-right=x, top-bot=y
	float const sz_scale(obj_bcube.get_size().sum() / bcube.get_size().sum());
	float const z_offset(0.5*bcube.dy() - (pos.z - obj_bcube.z1())/sz_scale); // translate required to map bottom of model to bottom of obj_bcube post transform
	
	if (enable_animations) {
		s.add_uniform_float("animation_scale",    model_file.scale/sz_scale); // Note: determined somewhat experimentally
		s.add_uniform_float("model_delta_height", (0.1*bcube.dy() + bcube.y1()));
	}
	fgPushMatrix();
	translate_to(pos + vector3d(0.0, 0.0, z_offset*sz_scale)); // z_offset is in model space, scale to world space
	if (fabs(dir.y) > 0.001) {rotate_to_plus_x(dir);} // orient facing front
	else if (dir.x < 0.0) {fgRotate(180.0, 0.0, 0.0, 1.0);}
	if (dir.z != 0.0) {fgRotate(TO_DEG*asinf(-dir.z), 0.0, 1.0, 0.0);} // handle cars on a slope
	if (model_file.xy_rot != 0.0) {fgRotate(model_file.xy_rot, 0.0, 0.0, 1.0);} // apply model rotation about z/up axis
	fgRotate(90.0, 1.0, 0.0, 0.0); // swap Y and Z dirs; models have up=Y, but we want up=Z
	uniform_scale(sz_scale); // scale from model space to the world space size of our target cube, using a uniform scale based on the averages of the x,y,z sizes
	translate_to(-bcube.get_cube_center()); // cancel out model local translate

	if ((low_detail || is_shadow_pass) && !model_file.shadow_mat_ids.empty()) { // low detail pass, normal maps disabled
		if (!is_shadow_pass && use_model3d_bump_maps()) {model3d::bind_default_flat_normal_map();} // still need to set the default here in case the shader is using it
		// TODO: combine shadow materials into a single VBO and draw with one call when is_shadow_pass==1; this is complex and may not yield a significant improvement
		for (auto i = model_file.shadow_mat_ids.begin(); i != model_file.shadow_mat_ids.end(); ++i) {model.render_material(s, *i, is_shadow_pass, 0, 2, 0);}
	}
	else {
		model.render_materials(s, is_shadow_pass, 0, 0, 2, 3, 3, model.get_unbound_material(), rotation_t(),
			nullptr, nullptr, is_shadow_pass, model_file.lod_mult, (is_shadow_pass ? 10.0 : 0.0)); // enable_alpha_mask=2 (both)
	}
	fgPopMatrix();
	camera_pdu.valid = camera_pdu_valid;
	camera_pdu.pos   = orig_camera_pos;
	select_texture(WHITE_TEX); // reset back to default/untextured
}


void ao_draw_state_t::draw_ao_qbd() {
	if (ao_qbd.empty()) return;
	enable_blend();
	select_texture(BLUR_CENT_TEX);
	ao_qbd.draw_and_clear();
	select_texture(WHITE_TEX); // reset back to default/untextured
	disable_blend();
}

void ao_draw_state_t::occlusion_checker_t::set_camera(pos_dir_up const &pdu) {
	if ((display_mode & 0x08) == 0) {state.building_ids.clear(); return;} // testing
	pos_dir_up near_pdu(pdu);
	near_pdu.far_ = 2.0*city_params.road_spacing; // set far clipping plane to one city block
	get_building_occluders(near_pdu, state);
	//cout << "occluders: " << state.building_ids.size() << endl;
}

bool ao_draw_state_t::occlusion_checker_t::is_occluded(cube_t const &c) {
	if (state.building_ids.empty()) return 0;
	float const z(c.z2()); // top edge
	point const corners[4] = {point(c.x1(), c.y1(), z), point(c.x2(), c.y1(), z), point(c.x2(), c.y2(), z), point(c.x1(), c.y2(), z)};
	return check_pts_occluded(corners, 4, state);
}

void ao_draw_state_t::pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_) {
	draw_state_t::pre_draw(xlate_, use_dlights_, shadow_only_, 1); // always_setup_shader=1 (required for model drawing)
	if (!shadow_only) {occlusion_checker.set_camera(camera_pdu);}
}

/*static*/ float car_draw_state_t::get_headlight_dist() {return 3.5*city_params.road_width;} // distance headlights will shine

colorRGBA car_draw_state_t::get_headlight_color(car_t const &car) const {
	return colorRGBA(1.0, 1.0, (1.0 + 0.8*(fract(1000.0*car.max_speed) - 0.5)), 1.0); // slight yellow-blue tinting using max_speed as a hash
}

void car_draw_state_t::pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_) {
	//set_enable_normal_map(use_model3d_bump_maps()); // used only for some car models, and currently doesn't work
	ao_draw_state_t::pre_draw(xlate_, use_dlights_, shadow_only_);
	select_texture(WHITE_TEX);
}

void car_draw_state_t::draw_unshadowed() {
	qbds[0].draw_and_clear();
	draw_ao_qbd();
}

void car_draw_state_t::add_car_headlights(vector<car_t> const &cars, vector3d const &xlate_, cube_t &lights_bcube) {
	xlate = xlate_; // needed earlier in the flow
	for (auto i = cars.begin(); i != cars.end(); ++i) {add_car_headlights(*i, lights_bcube);}
}

void car_draw_state_t::gen_car_pts(car_t const &car, bool include_top, point pb[8], point pt[8]) const {
	point const center(car.get_center());
	cube_t const &c(car.bcube);
	float const z1(center.z - 0.5*car.height), z2(center.z + 0.5*car.height), zmid(center.z + (include_top ? 0.1 : 0.5)*car.height), length(car.get_length());
	bool const dim(car.dim), dir(car.dir);
	set_cube_pts(c, z1, zmid, dim, dir, pb); // bottom

	if (include_top) {
		cube_t top_part(c);
		top_part.d[dim][0] += (dir ? 0.25 : 0.30)*length; // back
		top_part.d[dim][1] -= (dir ? 0.30 : 0.25)*length; // front
		set_cube_pts(top_part, zmid, z2, dim, dir, pt); // top
	}
	if (car.dz != 0.0) { // rotate all points about dim !d
		float const sine_val((dir ? 1.0f : -1.0f)*car.dz/length), cos_val(sqrt(1.0f - sine_val*sine_val));
		rotate_pts(center, sine_val, cos_val, dim, 2, pb);
		if (include_top) {rotate_pts(center, sine_val, cos_val, dim, 2, pt);}
	}
	if (car.rot_z != 0.0) { // turning about the z-axis: rot_z of [0.0, 1.0] maps to angles of [0.0, PI/2=90 degrees]
		float const sine_val(sinf(0.5f*PI*car.rot_z)), cos_val(sqrt(1.0f - sine_val*sine_val));
		rotate_pts(center, sine_val, cos_val, 0, 1, pb);
		if (include_top) {rotate_pts(center, sine_val, cos_val, 0, 1, pt);}
	}
}

bool sphere_in_light_cone_approx(pos_dir_up const &pdu, point const &center, float radius) {
	float const dist(p2p_dist(pdu.pos, center)), radius_at_dist(dist*pdu.sterm), rmod(radius_at_dist + radius);
	return pt_line_dist_less_than(center, pdu.pos, (pdu.pos + pdu.dir), rmod);
}

void car_draw_state_t::draw_car(car_t const &car, bool is_dlight_shadows) { // Note: all quads
	if (car.destroyed) return;

	if (is_dlight_shadows) { // dynamic spotlight shadow
		point const center(car.get_center());
		if (!dist_less_than(camera_pdu.pos, center, 0.6*camera_pdu.far_)) return; // optimization
		// since we know the dlight is a spotlight with a cone shape rather than a frustum, we can do a tighter visibility test
		if (!sphere_in_light_cone_approx(camera_pdu, center, car.bcube.get_xy_bsphere_radius())) return;
		cube_t bcube(car.bcube);
		bcube.expand_by(0.1*car.height);
		if (bcube.contains_pt(camera_pdu.pos)) return; // don't self-shadow
	}
	if (!check_cube_visible(car.bcube, (shadow_only ? 0.0 : 0.75))) return; // dist_scale=0.75
	point const center(car.get_center());
	begin_tile(center); // enable shadows
	colorRGBA const &color(car.get_color());
	float const tile_draw_dist(get_draw_tile_dist()), dist_val(p2p_dist(camera_pdu.pos, (center + xlate))/tile_draw_dist);
	bool const is_truck(car.height > 1.2*city_params.get_nom_car_size().z); // hack - truck has a larger than average size
	bool const draw_top(dist_val < 0.25 && !is_truck), dim(car.dim), dir(car.dir);
	bool const draw_model(car_model_loader.num_models() > 0 &&
		(is_dlight_shadows ? dist_less_than(pre_smap_player_pos, center, 0.05*tile_draw_dist) : (shadow_only || dist_val < 0.05)));
	float const sign((dim^dir) ? -1.0 : 1.0);
	point pb[8], pt[8]; // bottom and top sections
	gen_car_pts(car, draw_top, pb, pt);

	if (draw_model && car_model_loader.is_model_valid(car.model_id)) {
		if (is_occluded(car.bcube)) return; // only check occlusion for expensive car models
		vector3d const front_n(cross_product((pb[5] - pb[1]), (pb[0] - pb[1])).get_norm()*sign);
		car_model_loader.draw_model(s, center, car.bcube, front_n, color, xlate, car.model_id, shadow_only, (dist_val > 0.035));
	}
	else { // draw simple 1-2 cube model
		quad_batch_draw &qbd(qbds[emit_now]);
		color_wrapper cw(color);
		draw_cube(qbd, cw, center, pb, 1, (dim^dir)); // bottom (skip_bottom=1)
		if (draw_top) {draw_cube(qbd, cw, center, pt, 1, (dim^dir));} // top (skip_bottom=1)
		if (emit_now) {qbds[1].draw_and_clear();} // shadowed (only emit when tile changes?)
	}
	if (shadow_only) return; // shadow pass - done

	if (dist_val < 0.04 && fabs(car.dz) < 0.01) { // add AO planes when close to the camera and on a level road
		float const length(car.get_length());
		point pao[4];

		for (unsigned i = 0; i < 4; ++i) {
			point &v(pao[i]);
			v = pb[i] - center;
			v[ dim] += 0.1*length*SIGN(v[ dim]); // increase length slightly
			v[!dim] += 0.1*length*SIGN(v[!dim]); // increase width  slightly
			v   += center;
			v.z += 0.02*car.height; // shift up slightly to avoid z-fighting
		}
		/*if (!car.headlights_on()) { // daytime, adjust shadow to match sun pos
			vector3d const sun_dir(0.5*length*(center - get_sun_pos()).get_norm());
			vector3d const offset(sun_dir.x, sun_dir.y, 0.0);
			for (unsigned i = 0; i < 4; ++i) {pao[i] += offset;} // problems: double shadows, non-flat surfaces, buildings, texture coords/back in center, non-rectangular
		}*/
		ao_qbd.add_quad_pts(pao, colorRGBA(0, 0, 0, 0.9), plus_z);
	}
	if (dist_val > 0.3)  return; // to far - no lights to draw
	if (car.is_parked()) return; // no lights when parked
	vector3d const front_n(cross_product((pb[5] - pb[1]), (pb[0] - pb[1])).get_norm()*sign);
	unsigned const lr_xor(((camera_pdu.pos[!dim] - xlate[!dim]) - center[!dim]) < 0.0f);
	bool const brake_lights_on(car.is_almost_stopped() || car.stopped_at_light), headlights_on(car.headlights_on());

	if (headlights_on && dist_val < 0.3) { // night time headlights
		colorRGBA const hl_color(get_headlight_color(car));

		for (unsigned d = 0; d < 2; ++d) { // L, R
			unsigned const lr(d ^ lr_xor ^ 1);
			point const pos((lr ? 0.2 : 0.8)*(0.2*pb[0] + 0.8*pb[4]) + (lr ? 0.8 : 0.2)*(0.2*pb[1] + 0.8*pb[5]));
			add_light_flare(pos, front_n, hl_color, 2.0, 0.65*car.height); // pb 0,1,4,5
		}
	}
	if ((brake_lights_on || headlights_on) && dist_val < 0.2) { // brake lights
		for (unsigned d = 0; d < 2; ++d) { // L, R
			unsigned const lr(d ^ lr_xor);
			point const pos((lr ? 0.2 : 0.8)*(0.2*pb[2] + 0.8*pb[6]) + (lr ? 0.8 : 0.2)*(0.2*pb[3] + 0.8*pb[7]));
			add_light_flare(pos, -front_n, colorRGBA(1.0, 0.1, 0.05, 1.0), (brake_lights_on ? 1.0 : 0.5), 0.5*car.height); // near red; pb 2,3,6,7
		}
	}
	if (car.turn_dir != TURN_NONE && car.cur_city != CONN_CITY_IX && dist_val < 0.1) { // turn signals (not on connector road bends)
		float const ts_period = 1.5; // in seconds
		double const time(fract((tfticks + 1000.0*car.max_speed)/(double(ts_period)*TICKS_PER_SECOND))); // use car max_speed as seed to offset time base

		if (time > 0.5) { // flash on and off
			bool const tdir((car.turn_dir == TURN_LEFT) ^ dim ^ dir); // R=1,2,5,6 or L=0,3,4,7
			vector3d const side_n(cross_product((pb[6] - pb[2]), (pb[1] - pb[2])).get_norm()*sign*(tdir ? 1.0 : -1.0));

			for (unsigned d = 0; d < 2; ++d) { // B, F
				point const pos(0.3*pb[tdir ? (d ? 1 : 2) : (d ? 0 : 3)] + 0.7*pb[tdir ? (d ? 5 : 6) : (d ? 4 : 7)]);
				add_light_flare(pos, (side_n + (d ? 1.0 : -1.0)*front_n).get_norm(), colorRGBA(1.0, 0.75, 0.0, 1.0), 1.5, 0.3*car.height); // normal points out 45 degrees
			}
		}
	}
}

void car_draw_state_t::add_car_headlights(car_t const &car, cube_t &lights_bcube) {
	if (!car.headlights_on()) return;
	float const headlight_dist(get_headlight_dist());
	cube_t bcube(car.bcube);
	bcube.expand_by(headlight_dist);
	if (!lights_bcube.contains_cube_xy(bcube))   return; // not contained within the light volume
	if (!camera_pdu.cube_visible(bcube + xlate)) return; // VFC
	float const sign((car.dim^car.dir) ? -1.0 : 1.0);
	point pb[8], pt[8]; // bottom and top sections
	gen_car_pts(car, 0, pb, pt); // draw_top=0
	vector3d const front_n(cross_product((pb[5] - pb[1]), (pb[0] - pb[1])).get_norm()*sign);
	vector3d const dir((0.5*front_n - 0.5*plus_z).get_norm()); // point slightly down
	colorRGBA const color(get_headlight_color(car));
	float const beamwidth = 0.08;
	min_eq(lights_bcube.z1(), bcube.z1());
	max_eq(lights_bcube.z2(), bcube.z2());

	// TODO: merge headlights into a single spotlight for distant cars? share shadow maps between headlights?
	for (unsigned d = 0; d < 2; ++d) { // L, R
		point const pos((d ? 0.2 : 0.8)*(0.2*pb[0] + 0.8*pb[4]) + (d ? 0.8 : 0.2)*(0.2*pb[1] + 0.8*pb[5]));
		dl_sources.push_back(light_source(headlight_dist, pos, pos, color, 1, dir, beamwidth));
	}
}


void car_manager_t::remove_destroyed_cars() {
	remove_destroyed(cars);
	car_destroyed = 0;
}

void car_manager_t::init_cars(unsigned num) {
	if (num == 0) return;
	timer_t timer("Init Cars");
	cars.reserve(num);
	for (unsigned n = 0; n < num; ++n) {add_car();}
	cout << "Dynamic Cars: " << cars.size() << endl;
}

void car_manager_t::add_parked_cars(vector<car_t> const &new_cars, vect_cube_t const &garages) {
	first_parked_car = cars.size(); // Note: sort may invalidate this, but okay for use in finalize_cars()
	cars.reserve(cars.size() + new_cars.size() + garages.size());
	vector_add_to(new_cars, cars);
	first_garage_car = cars.size(); // Note: sort may invalidate this, but okay for use in finalize_cars()
	if (garages.empty()) return; // done
	vector3d const nom_car_size(city_params.get_nom_car_size());
	car_t car; // no cur_city/cur_road/cur_segq
	car.park();
	car.cur_city      = NO_CITY_IX; // special value
	car.cur_road_type = TYPE_BUILDING; // garage
	rand_gen_t rgen;
	
	for (auto i = garages.begin(); i != garages.end(); ++i) {
		if ((rgen.rand()&3) == 0) continue; // 25% of garages have no car
		vector3d car_sz(nom_car_size);
		car.dim    = (i->dx() < i->dy()); // long dim
		car.dir    = rgen.rand_bool(); // Note: ignores garage dir because some cars and backed in and some are pulled in
		car.height = car_sz.z;
		if (car.dim) {swap(car_sz.x, car_sz.y);}
		car.bcube.set_from_point(i->get_cube_center());
		car.bcube.expand_by(0.5*car_sz);
		assert(i->contains_cube(car.bcube));
		car.bcube.z1() = i->z1(); car.bcube.z2() = i->z1() + car.height;
		cars.push_back(car);
		garages_bcube.assign_or_union_with_cube(car.bcube);
	} // for i
}

void car_manager_t::finalize_cars() {
	if (empty()) return;
	unsigned const num_models(car_model_loader.num_models());

	for (auto i = cars.begin(); i != cars.end(); ++i) {
		int fixed_color(-1);

		if (num_models > 0) {
			for (unsigned n = 0; n < 20; ++n) {
				if (FORCE_MODEL_ID >= 0) {i->model_id = (unsigned char)FORCE_MODEL_ID;}
				else {i->model_id = ((num_models > 1) ? (rgen.rand() % num_models) : 0);}
				city_model_t const &model(car_model_loader.get_model(i->model_id));
				// if there are multiple models to choose from, and this car is in a garage, try for a model that's not scaled up (the truck)
				if (FORCE_MODEL_ID < 0 && num_models > 1 && (i-cars.begin()) >= first_garage_car && n+1 < 20 && model.scale > 1.0) continue;
				fixed_color = model.fixed_color_id;
				i->apply_scale(model.scale);
				break;
			} // for n
		}
		i->color_id = ((fixed_color >= 0) ? fixed_color : (rgen.rand() % NUM_CAR_COLORS));
		assert(i->is_valid());
	} // for i
	cout << "Total Cars: " << cars.size() << endl; // 4000 on the road + 4372 parked + 433 garage (out of 594) = 8805
}

void car_city_vect_t::clear_cars() {
	for (unsigned d = 0; d < 2; ++d) {cars[d][0].clear(); cars[d][1].clear();}
}

void car_manager_t::extract_car_data(vector<car_city_vect_t> &cars_by_city) const {
	if (cars.empty()) return;
	//timer_t timer("Extract Car Data");
	// create parked cars vectors on first call; this is used for pedestrian navigation within parking lots;
	// it won't be rebuilt on car destruction, but that should be okay
	bool const add_parked_cars(cars_by_city.empty());
	for (auto i = cars_by_city.begin(); i != cars_by_city.end(); ++i) {i->clear_cars();} // clear prev frame's state

	for (auto i = cars.begin(); i != cars.end(); ++i) {
		if (i->cur_city >= cars_by_city.size()) {cars_by_city.resize(i->cur_city+1);}
		auto &dest(cars_by_city[i->cur_city]);
		if (!i->is_parked()) {dest.cars[i->dim][i->dir].push_back(*i);} // moving on road
		else if (add_parked_cars) {dest.parked_car_bcubes.push_back(i->bcube);} // parked, not yet updated
	}
}

bool car_manager_t::proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm) const {
	vector3d const xlate(get_camera_coord_space_xlate());
	float const dist(p2p_dist(pos, p_last));

	for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
		cube_t const city_bcube(get_cb_bcube(*cb) + xlate);
		if (pos.z - radius > city_bcube.z2() + city_params.get_max_car_size().z) continue; // above the cars
		if (!sphere_cube_intersect_xy(pos, (radius + dist), city_bcube)) continue;
		cube_t sphere_bc; sphere_bc.set_from_sphere((pos - xlate), radius);
		unsigned start(0), end(0);
		get_car_ix_range_for_cube(cb, sphere_bc, start, end);

		for (unsigned c = start; c != end; ++c) {
			if (cars[c].proc_sphere_coll(pos, p_last, radius, xlate, cnorm)) return 1;
		}
	} // for cb
	return 0;
}

void car_manager_t::destroy_cars_in_radius(point const &pos_in, float radius) {
	vector3d const xlate(get_camera_coord_space_xlate());
	point const pos(pos_in - xlate);
	bool const is_pt(radius == 0.0);

	for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
		cube_t const city_bcube(get_cb_bcube(*cb));
		if (pos.z - radius > city_bcube.z2() + city_params.get_max_car_size().z) continue; // above the cars
		if (is_pt ? !city_bcube.contains_pt_xy(pos) : !sphere_cube_intersect_xy(pos, radius, city_bcube)) continue;
		unsigned const start(cb->start), end((cb+1)->start); // Note: shouldnt be called frequently enough to need road/parking lot acceleration
		assert(end <= cars.size() && start <= end);

		for (unsigned c = start; c != end; ++c) {
			car_t &car(cars[c]);

			if (is_pt ? car.bcube.contains_pt(pos) : dist_less_than(car.get_center(), pos, radius)) { // destroy if within the sphere
				car.destroy();
				car_destroyed = 1;
				// invalidate tile shadow map for destroyed parked cars
				if (city_params.car_shadows && car.is_parked()) {invalidate_tile_smap_at_pt((car.get_center() + xlate), 0.5*car.get_length());} // radius = length/2
			}
		} // for c
	} // for cb
}

bool car_manager_t::get_color_at_xy(point const &pos, colorRGBA &color, int int_ret) const { // Note: pos in local TT space
	if (cars.empty()) return 0;
	if (int_ret != INT_ROAD && int_ret != INT_PARKING) return 0; // not a road or a parking lot - no car intersections

	for (auto cb = car_blocks_by_road.begin(); cb+1 < car_blocks_by_road.end(); ++cb) { // use cars_by_road to accelerate query
		if (!get_cb_bcube(*cb).contains_pt_xy(pos)) continue; // skip
		unsigned start(cb->start), end((cb+1)->start);
		if      (int_ret == INT_ROAD)    {end   = cb->first_parked;} // moving cars only (beginning of range)
		else if (int_ret == INT_PARKING) {start = cb->first_parked;} // parked cars only (end of range)
		assert(start <= end);
		assert(end < cars_by_road.size()); // strictly less

		for (unsigned i = start; i != end; ++i) {
			cube_with_ix_t const &v(cars_by_road[i]);
			if (!v.contains_pt_xy(pos)) continue; // skip
			unsigned const ix_end(cars_by_road[i+1].ix);
			assert(ix_end <= cars.size());

			for (unsigned c = v.ix; c != ix_end; ++c) {
				if (cars[c].bcube.contains_pt_xy(pos)) {color = cars[c].get_color(); return 1;}
			}
		}
	} // for cb
	return 0;
}

car_t const *car_manager_t::get_car_at_pt(point const &pos, bool is_parked) const {
	for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
		if (!get_cb_bcube(*cb).contains_pt_xy(pos)) continue; // skip
		unsigned start(cb->start), end((cb+1)->start);
		if (!is_parked) {end   = cb->first_parked;} // moving cars only (beginning of range)
		else            {start = cb->first_parked;} // parked cars only (end of range)
		if (start > end || end > cars.size()) {cout << TXT(start) << TXT(end) << TXT(cars.size()) << TXT(is_parked) << endl;}
		assert(start <= end && end <= cars.size());

		for (unsigned c = start; c != end; ++c) {
			if (cars[c].bcube.contains_pt_xy(pos)) {return &cars[c];}
		}
	} // for cb
	return nullptr; // no car found
}

car_t const *car_manager_t::get_car_at(point const &p1, point const &p2) const { // Note: p1/p2 in local TT space
	for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
		if (!get_cb_bcube(*cb).line_intersects(p1, p2)) continue; // skip
		unsigned start(cb->start), end((cb+1)->start);
		assert(start <= end && end <= cars.size());

		for (unsigned c = start; c != end; ++c) { // Note: includes parked cars
			if (cars[c].bcube.line_intersects(p1, p2)) {return &cars[c];}
		}
	} // for cb
	return nullptr; // no car found
}

bool car_manager_t::line_intersect_cars(point const &p1, point const &p2, float &t) const { // Note: p1/p2 in local TT space
	bool ret(0);

	for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
		if (!get_cb_bcube(*cb).line_intersects(p1, p2)) continue; // skip
		unsigned start(cb->start), end((cb+1)->start);
		assert(start <= end && end <= cars.size());

		for (unsigned c = start; c != end; ++c) { // Note: includes parked cars
			ret |= check_line_clip_update_t(p1, p2, t, cars[c].bcube);
		}
	} // for cb
	return ret;
}

int car_manager_t::find_next_car_after_turn(car_t &car) {
	road_isec_t const &isec(get_car_isec(car));
	if (car.turn_dir == TURN_NONE && !isec.is_global_conn_int()) return -1; // car not turning, and not on connector road isec: should be handled by sorted car_in_front logic
	unsigned const dest_orient(isec.get_dest_orient_for_car_in_isec(car, 0)); // Note: may be before, during, or after turning
	int road_ix(isec.rix_xy[dest_orient]), seg_ix(isec.conn_ix[dest_orient]);
	unsigned city_ix(car.cur_city);
	//cout << TXT(car.get_orient()) << TXT(dest_orient) << TXT(city_ix) << TXT(road_ix) << TXT(seg_ix) << endl;
	assert((road_ix < 0) == (seg_ix < 0));

	if (road_ix < 0) { // goes to connector road
		city_ix = CONN_CITY_IX;
		road_ix = decode_neg_ix(road_ix);
		seg_ix  = decode_neg_ix(seg_ix );
	}
	point const car_center(car.get_center());
	float dmin(car.get_max_lookahead_dist()), dmin_sq(dmin*dmin);
	// include normal sorted order car; this is needed when going straight through connector road 4-way intersections where cur_road changes within the intersection
	if (car.car_in_front && car.car_in_front->get_orient() != dest_orient) {car.car_in_front = 0;} // not the correct car (turning a different way)
	if (car.turn_dir == TURN_NONE && car.car_in_front) {min_eq(dmin_sq, p2p_dist_sq(car_center, car.car_in_front->get_center()));}
	int ret_car_ix(-1);

	for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
		if (cb->cur_city != city_ix) continue; // incorrect city - skip
		unsigned const start(cb->start), end(cb->first_parked);
		assert(end <= cars.size() && start <= end);
		auto range_end(cars.begin()+end);
		car_t ref_car; ref_car.cur_road = road_ix;
		auto it(std::lower_bound(cars.begin()+start, range_end, ref_car, comp_car_road())); // binary search acceleration
		float prev_dist_sq(FLT_MAX);

		for (; it != range_end; ++it) {
			if (&(*it) == &car) continue; // skip self
			assert(it->cur_city == city_ix); // must be same city
			if (it->cur_road != road_ix) break; // different road, done

			if (it->cur_road_type == TYPE_RSEG) { // road segment
				if (it->cur_seg != seg_ix) continue; // on a different segment, skip
			}
			else if (it->cur_road_type != car.cur_road_type || it->cur_seg != car.cur_seg) continue; // in a different intersection
			if (it->get_orient() != dest_orient) continue; // wrong orient
			float const dist_sq(p2p_dist_sq(car_center, it->get_center()));
			if (p2p_dist_sq(car_center, it->get_front()) < dist_sq) continue; // front is closer than back - this car is not in front of us (waiting on other side of isect?)
			//cout << TXT(dmin_sq) << TXT(dist_sq) << (dist_sq < dmin_sq) << endl;

			if (dist_sq < dmin_sq) { // new closest car
				if (&(*it) != car.car_in_front) {ret_car_ix = (it - cars.begin());} // record index if set to a new value
				car.car_in_front = &(*it);
				dmin_sq = dist_sq;
			}
			else if (dist_sq > prev_dist_sq) break; // we're moving too far away from the car
			prev_dist_sq = dist_sq;
		} // for it
	} // for cb
	return ret_car_ix;
}

bool car_manager_t::check_car_for_ped_colls(car_t &car) const {
	if (car.cur_city >= peds_crossing_roads.peds.size())  return 0; // no peds in this city (includes connector road network)
	if (car.turn_val != 0.0 || car.turn_dir != TURN_NONE) return 0; // for now, don't check for cars when turning as this causes problems with blocked intersections
	auto const &peds_by_road(peds_crossing_roads.peds[car.cur_city]);
	if (car.cur_road >= peds_by_road.size()) return 0; // no peds in this road
	auto const &peds(peds_by_road[car.cur_road]);
	if (peds.empty()) return 0;
	cube_t coll_area(car.bcube);
	coll_area.d[car.dim][!car.dir] = coll_area.d[car.dim][car.dir]; // exclude the car itself
	coll_area.d[car.dim][car.dir] += (car.dir ? 1.25 : -1.25)*car.get_length(); // extend the front
	coll_area.d[!car.dim][0] -= 0.5*car.get_width();
	coll_area.d[!car.dim][1] += 0.5*car.get_width();
	static rand_gen_t rgen;

	for (auto i = peds.begin(); i != peds.end(); ++i) {
		if (coll_area.contains_pt_xy_exp(i->pos, i->radius)) {
			car.decelerate_fast();
			if ((rgen.rand()&3) == 0) {car.honk_horn_if_close_and_fast();}
			return 1;
		}
	} // for i
	return 0;
}

bool car_manager_t::choose_dest_parked_car(unsigned city_id, unsigned &plot_id, unsigned &car_ix, point &car_center, rand_gen_t &rgen) const { // for pedestrians
	if (cars.empty()) return 0; // no cars
	
	for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
		if (cb->cur_city != city_id) continue; // incorrect city - skip
		unsigned start(cb->first_parked), end((cb+1)->start);
		if (start == end) continue; // no parked cars
		assert(start < end && end <= cars.size());
		car_ix  = rgen.rand_int(start, end-1);
		assert(cars[car_ix].is_parked());
		plot_id    = cars[car_ix].cur_road; // plot_id is stored in cur_road for parked cars
		car_center = cars[car_ix].get_center();
		return 1;
	}
	return 0;
}

void car_manager_t::next_frame(ped_manager_t const &ped_manager, float car_speed) {
	if (cars.empty() || !animate2) return;
	// Warning: not really thread safe, but should be okay; the ped state should valid at all points (thought maybe inconsistent) and we don't need it to be exact every frame
	ped_manager.get_peds_crossing_roads(peds_crossing_roads);
	//timer_t timer("Update Cars"); // 4K cars = 0.7ms / 2.1ms with destinations + navigation
#pragma omp critical(modify_car_data)
	{
		if (car_destroyed) {remove_destroyed_cars();} // at least one car was destroyed in the previous frame - remove it/them
		sort(cars.begin(), cars.end(), comp_car_road_then_pos(camera_pdu.pos - dstate.xlate)); // sort by city/road/position for intersection tests and tile shadow map binds
	}
	entering_city.clear();
	car_blocks.clear();
	float const speed(CAR_SPEED_SCALE*car_speed*fticks);
	bool saw_parked(0);
	//unsigned num_on_conn_road(0);

	for (auto i = cars.begin(); i != cars.end(); ++i) { // move cars
		unsigned const cix(i - cars.begin());
		i->car_in_front = nullptr; // reset for this frame

		if (car_blocks.empty() || i->cur_city != car_blocks.back().cur_city) {
			if (!saw_parked && !car_blocks.empty()) {car_blocks.back().first_parked = cix;} // no parked cars in prev city
			saw_parked = 0; // reset for next city
			car_blocks.emplace_back(cix, i->cur_city);
		}
		if (i->is_parked()) {
			if (!saw_parked) {car_blocks.back().first_parked = cix; saw_parked = 1;}
			continue; // no update for parked cars
		}
		i->move(speed);
		if (i->entering_city) {entering_city.push_back(cix);} // record for use in collision detection
		if (!i->stopped_at_light && i->is_almost_stopped() && i->in_isect()) {get_car_isec(*i).stoplight.mark_blocked(i->dim, i->dir);} // blocking intersection
		register_car_at_city(*i);
	} // for i
	if (!saw_parked && !car_blocks.empty()) {car_blocks.back().first_parked = cars.size();} // no parked cars in final city
	car_blocks.emplace_back(cars.size(), 0); // add terminator

	for (auto i = cars.begin(); i != cars.end(); ++i) { // collision detection
		if (i->is_parked()) continue; // no collisions for parked cars
		bool const on_conn_road(i->cur_city == CONN_CITY_IX);
		float const length(i->get_length()), max_check_dist(max(3.0f*length, (length + i->get_max_lookahead_dist()))); // max of collision dist and car-in-front dist

		for (auto j = i+1; j != cars.end(); ++j) { // check for collisions with cars on the same road (can't test seg because they can be on diff segs but still collide)
			if (i->cur_city != j->cur_city || i->cur_road != j->cur_road) break; // different cities or roads
			if (!on_conn_road && i->cur_road_type == j->cur_road_type && abs((int)i->cur_seg - (int)j->cur_seg) > (on_conn_road ? 1 : 0)) break; // diff road segs or diff isects
			check_collision(*i, *j);
			i->register_adj_car(*j);
			j->register_adj_car(*i);
			if (!dist_xy_less_than(i->get_center(), j->get_center(), max_check_dist)) break;
		}
		if (on_conn_road) { // on connector road, check before entering intersection to a city
			for (auto ix = entering_city.begin(); ix != entering_city.end(); ++ix) {
				if (*ix != unsigned(i - cars.begin())) {check_collision(*i, cars[*ix]);}
			}
			//++num_on_conn_road;
		}
		if (i->in_isect()) {
			int const next_car(find_next_car_after_turn(*i)); // Note: calculates in i->car_in_front
			if (next_car >= 0) {check_collision(*i, cars[next_car]);} // make sure we collide with the correct car
		}
		if (!peds_crossing_roads.peds.empty()) {check_car_for_ped_colls(*i);}
	} // for i
	update_cars(); // run update logic

	if (map_mode) { // create cars_by_road
		// cars have moved since the last sort and may no longer be in city/road order, but this algorithm doesn't require that;
		// out-of-order cars will end up in their own blocks, which is less efficient but still correct
		car_blocks_by_road.clear();
		cars_by_road.clear();
		unsigned cur_city(1<<31), cur_road(1<<31); // start at invalid values
		bool saw_parked(0);

		for (auto i = cars.begin(); i != cars.end(); ++i) {
			bool const new_city(i->cur_city != cur_city), new_parked(!saw_parked && i->is_parked());
			unsigned const cbr_ix(cars_by_road.size());
			if (new_parked) {car_blocks_by_road.back().first_parked = cbr_ix; saw_parked = 1;}

			if (new_city || new_parked || i->cur_road != cur_road) { // new city/road
				if (new_city) { // new city
					if (!saw_parked && !car_blocks_by_road.empty()) {car_blocks_by_road.back().first_parked = cbr_ix;} // no parked cars in prev city
					saw_parked = 0; // reset for next city
					car_blocks_by_road.emplace_back(cbr_ix, i->cur_city);
				}
				cars_by_road.emplace_back(i->bcube, (i - cars.begin())); // start a new block
				cur_city = i->cur_city;
				cur_road = i->cur_road;
			}
			else {cars_by_road.back().union_with_cube(i->bcube);}
		} // for i
		if (!saw_parked && !car_blocks_by_road.empty()) {car_blocks_by_road.back().first_parked = cars_by_road.size();} // no parked cars in final city
		car_blocks_by_road.emplace_back(cars_by_road.size(), 0); // add terminator
		cars_by_road.emplace_back(cube_t(), cars.size()); // add terminator
	}
	//cout << TXT(cars.size()) << TXT(entering_city.size()) << TXT(in_isects.size()) << TXT(num_on_conn_road) << endl; // TESTING
}

void car_manager_t::draw(int trans_op_mask, vector3d const &xlate, bool use_dlights, bool shadow_only, bool is_dlight_shadows) {
	if (cars.empty()) return;

	if (trans_op_mask & 1) { // opaque pass, should be first
		if (is_dlight_shadows && !city_params.car_shadows) return;
		bool const only_parked(shadow_only && !is_dlight_shadows); // sun/moon shadows are precomputed and cached, so only include static objects such as parked cars
		//timer_t timer(string("Draw Cars") + (shadow_only ? " Shadow" : "")); // 10K cars = 1.5ms / 2K cars = 0.33ms
		dstate.xlate = xlate;
		fgPushMatrix();
		translate_to(xlate);
		dstate.pre_draw(xlate, use_dlights, shadow_only);

		for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
			if (!camera_pdu.cube_visible(get_cb_bcube(*cb) + xlate)) continue; // city not visible - skip
			unsigned const end((cb+1)->start);
			assert(end <= cars.size());

			for (unsigned c = cb->start; c != end; ++c) {
				if (only_parked && !cars[c].is_parked()) continue; // skip non-parked cars
				dstate.draw_car(cars[c], is_dlight_shadows);
			}
		} // for cb
		dstate.post_draw();
		fgPopMatrix();

		if (tt_fire_button_down && !game_mode) {
			point const p1(get_camera_pos() - xlate), p2(p1 + cview_dir*FAR_CLIP);
			car_t const *car(get_car_at(p1, p2));
			if (car != nullptr) {dstate.set_label_text(car->label_str(), (car->get_center() + xlate));} // car found
		}
	}
	if ((trans_op_mask & 2) && !shadow_only) {dstate.draw_and_clear_light_flares();} // transparent pass; must be done last for alpha blending, and no translate
	dstate.show_label_text();
}

