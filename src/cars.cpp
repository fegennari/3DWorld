// 3D World - Cars for Procedural Cities
// by Frank Gennari
// 11/19/18
#include "city.h"
#include "file_utils.h"
#include "openal_wrap.h"
#include "explosion.h" // for add_blastr()

float const MIN_CAR_STOP_SEP = 0.25; // in units of car lengths

extern city_params_t city_params;


bool car_model_t::read(FILE *fp) { // filename body_material_id fixed_color_id xy_rot dz lod_mult shadow_mat_ids

	assert(fp);
	if (!read_string(fp, fn)) return 0;
	if (!read_int(fp, body_mat_id)) return 0;
	if (!read_int(fp, fixed_color_id)) return 0;
	if (!read_float(fp, xy_rot)) return 0;
	if (!read_float(fp, dz)) return 0;
	if (!read_float(fp, scale)) return 0;
	if (!read_float(fp, lod_mult) || lod_mult <= 0.0) return 0;
	unsigned shadow_mat_id;
	while (read_uint(fp, shadow_mat_id)) {shadow_mat_ids.push_back(shadow_mat_id);}
	return 1;
}


float car_t::get_max_lookahead_dist() const {return (get_length() + city_params.road_width);} // extend one car length + one road width in front
float car_t::get_turn_rot_z(float dist_to_turn) const {return (1.0 - CLIP_TO_01(4.0f*fabs(dist_to_turn)/city_params.road_width));}

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
	float const avg_len(0.5*(get_length() + c.get_length())); // average length of the two cars
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
	if (fabs(cur_pos - waiting_pos) > get_length()) {waiting_pos = cur_pos; waiting_start = tfticks;} // update when we move at least a car length
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

point car_t::get_front(float dval) const {
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

unsigned car_t::count_cars_in_front(cube_t const &range) const {
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

	for (unsigned i = 0; i < 50; ++i) { // limit iterations; avg len = city_params.get_nom_car_size().x (FIXME: 50 not high enough for connector roads)
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
		return (p2p_dist_sq((c1.bcube.get_cube_center() + xlate), camera_pdu.pos) > p2p_dist_sq((c2.bcube.get_cube_center() + xlate), camera_pdu.pos));
	}
	return (c1.bcube.d[c1.dim][c1.dir] < c2.bcube.d[c2.dim][c2.dir]); // compare front end of car (used for collisions)
}


/*static*/ unsigned car_model_loader_t::num_models() {return city_params.car_model_files.size();}

bool car_model_loader_t::is_model_valid(unsigned id) {
	assert(id < num_models());
	ensure_models_loaded(); // I guess we have to load the models here to determine if they're valid
	assert(id < models_valid.size());
	return (models_valid[id] != 0);
}

car_model_t const &car_model_loader_t::get_model(unsigned id) const {
	assert(id < num_models());
	return city_params.car_model_files[id];
}

void car_model_loader_t::load_car_models() {
	models_valid.resize(num_models(), 1); // assume valid

	for (unsigned i = 0; i < num_models(); ++i) {
		string const &fn(get_model(i).fn);
		bool const recalc_normals = 1;

		if (!load_model_file(fn, *this, geom_xform_t(), -1, WHITE, 0, 0.0, recalc_normals, 0, 0, 1)) {
			cerr << "Error: Failed to read model file '" << fn << "'; Skipping this model (will use default box model)." << endl;
			push_back(model3d(fn, tmgr)); // add a placeholder dummy model
			models_valid[i] = 0;
		}
	} // for i
}

void car_model_loader_t::draw_car(shader_t &s, vector3d const &pos, cube_t const &car_bcube, vector3d const &dir, colorRGBA const &color,
	point const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail)
{
	assert(is_model_valid(model_id));
	assert(size() == num_models()); // must be loaded
	car_model_t const &model_file(get_model(model_id));
	model3d &model(at(model_id));

	if (!is_shadow_pass && model_file.body_mat_id >= 0) { // use custom color for body material
		material_t &body_mat(model.get_material(model_file.body_mat_id));
		body_mat.ka = body_mat.kd = color;
		//if (model_id == 5) {cout << body_mat.name << endl;}
	}
	model.bind_all_used_tids();
	cube_t const &bcube(model.get_bcube());
	point const orig_camera_pos(camera_pdu.pos);
	camera_pdu.pos += bcube.get_cube_center() - pos - xlate; // required for distance based LOD
	bool const camera_pdu_valid(camera_pdu.valid);
	camera_pdu.valid = 0; // disable VFC, since we're doing custom transforms here
							// Note: in model space, front-back=z, left-right=x, top-bot=y
	float const sz_scale(car_bcube.get_size().sum() / bcube.get_size().sum());
	fgPushMatrix();
	translate_to(pos + vector3d(0.0, 0.0, model_file.dz*sz_scale));
	if (fabs(dir.y) > 0.001) {rotate_to_plus_x(dir);}
	else if (dir.x < 0.0) {fgRotate(180.0, 0.0, 0.0, 1.0);}
	fgRotate(TO_DEG*asinf(-dir.z), 0.0, 1.0, 0.0);
	if (model_file.xy_rot != 0.0) {fgRotate(model_file.xy_rot, 0.0, 0.0, 1.0);}
	fgRotate(90.0, 1.0, 0.0, 0.0);
	uniform_scale(sz_scale);
	translate_to(-bcube.get_cube_center()); // cancel out model local translate

	if ((low_detail || is_shadow_pass) && !model_file.shadow_mat_ids.empty()) {
		for (auto i = model_file.shadow_mat_ids.begin(); i != model_file.shadow_mat_ids.end(); ++i) {model.render_material(s, *i, is_shadow_pass);}
	}
	else {
		auto const &unbound_mat(model.get_unbound_material());

		for (unsigned sam_pass = 0; sam_pass < (is_shadow_pass ? 2U : 1U); ++sam_pass) {
			model.render_materials(s, is_shadow_pass, 0, 0, (sam_pass == 1), 3, 3, unbound_mat, rotation_t(),
				nullptr, nullptr, is_shadow_pass, model_file.lod_mult, (is_shadow_pass ? 10.0 : 0.0));
		}
	}
	fgPopMatrix();
	camera_pdu.valid = camera_pdu_valid;
	camera_pdu.pos   = orig_camera_pos;
	select_texture(WHITE_TEX); // reset back to default/untextured
}

