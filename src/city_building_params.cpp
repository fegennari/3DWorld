// 3D World - City and Building Parameter Parsing
// by Frank Gennari
// 09/04/20

#include "city.h"
#include "buildings.h"


extern building_params_t global_building_params;
bool parse_buildings_option(FILE *fp) {return global_building_params.parse_buildings_option(fp);}


string const model_opt_names[NUM_OBJ_MODELS] =
{"toilet_model", "sink_model", "tub_model", "fridge_model", "stove_model", "tv_model", ""/*monitor*/, "couch_model", "office_chair_model", "urinal_model",
"lamp_model", "washer_model", "dryer_model", "key_model", "hanger_model", "clothing_model", "fire_escape_model", "wall_lamp_model", "cup_model", "toaster_model",
"hood_model", "rocking_chair_model", "silverware_model", "toy_model", "ceiling_fan_model", "fire_ext_model", "folded_shirt_model", "plant_model", "pool_table_model",
"pool_ladder_model", "bar_stool_model", "padlock_model", "cash_register_model", "water_fountain_model", "banana_model", "banana_peel_model", "phone_model", ""/*gbike*/, ""/*xfmr*/,
/*animal models*/ "rat_model", "roach_model",
/*building non-room objects*/ "door_handle_model",
/*city models*/ "fire_hydrant_model", "substation_model", "mailbox_model", "umbrella_model", "pigeon_model", "fountain_model", "bird_animated_model", "flag_model",
"bicycle_model", "swingset_model", "trampoline_model", "dumpster_model", "big_umbrella_model", "flower_model", "deck_chair_model"};

void city_params_t::init_kw_maps() {
	kwmu.add("num_cities",     num_cities);
	kwmu.add("num_rr_tracks",  num_rr_tracks);
	kwmu.add("num_samples",    num_samples);
	kwmu.add("num_conn_tries", num_conn_tries);
	kwmu.add("plots_to_parks_ratio", park_rate);
	kwmu.add("city_border", city_border); // in heightmap texels/mesh grids
	kwmu.add("road_border", road_border); // in heightmap texels/mesh grids
	kwmu.add("slope_width", slope_width); // in heightmap texels/mesh grids
	kwmb.add("assign_house_plots",     assign_house_plots);
	kwmb.add("new_city_conn_road_alg", new_city_conn_road_alg);
	kwmb.add("convert_model_files",    convert_model_files); // to model3d format; applies to cars, people, building items, etc.
	kwmr.add("road_width",   road_width,   FP_CHECK_NONNEG);
	kwmr.add("road_spacing", road_spacing, FP_CHECK_NONNEG);
	kwmr.add("road_spacing_rand",   road_spacing_rand,   FP_CHECK_NONNEG);
	kwmr.add("road_spacing_xy_add", road_spacing_xy_add, FP_CHECK_NONNEG);
	kwmr.add("conn_road_seg_len",   conn_road_seg_len,   FP_CHECK_POS);
	kwmr.add("max_road_slope",   max_road_slope,   FP_CHECK_POS);
	kwmr.add("max_track_slope",  max_track_slope,  FP_CHECK_POS);
	kwmr.add("model_anim_scale", model_anim_scale, FP_CHECK_POS);
	kwmr.add("residential_probability", residential_probability, FP_CHECK_01);
	kwmb.add("add_skyways", add_skyways);
	// cars
	kwmu.add("num_cars", num_cars);
	kwmb.add("enable_car_path_finding", enable_car_path_finding);
	kwmb.add("cars_use_driveways",  cars_use_driveways);
	kwmr.add("car_speed",           car_speed,           FP_CHECK_NONNEG);
	kwmr.add("traffic_balance_val", traffic_balance_val, FP_CHECK_01);
	kwmr.add("new_city_prob",       new_city_prob,       FP_CHECK_01);
	// pedestrians
	kwmu.add("num_peds", num_peds);
	kwmb.add("ped_respawn_at_dest", ped_respawn_at_dest);
	kwmb.add("use_animated_people", use_animated_people);
	kwmr.add("ped_speed",           ped_speed, FP_CHECK_NONNEG);
	// parking lots / trees / detail objects
	kwmu.add("min_park_spaces", min_park_spaces); // with default road parameters, can be up to 28
	kwmu.add("min_park_rows",   min_park_rows  ); // with default road parameters, can be up to 8
	kwmu.add("max_trees_per_plot",   max_trees_per_plot);
	kwmu.add("max_benches_per_plot", max_benches_per_plot);
	kwmr.add("min_park_density", min_park_density, FP_CHECK_01);
	kwmr.add("max_park_density", max_park_density, FP_CHECK_01);
	kwmr.add("tree_spacing",     tree_spacing,     FP_CHECK_POS);
	// lighting
	kwmu.add("max_lights",      max_lights);
	kwmu.add("max_shadow_maps", max_shadow_maps);
	kwmb.add("car_shadows",     car_shadows);
}
bool city_params_t::read_option(FILE *fp) {

	char strc[MAX_CHARS] = {0};
	if (!read_str(fp, strc)) return 0;
	string const str(strc);
	if (kwmb.maybe_set_from_fp(str, fp)) return !read_error_flag;
	if (kwmu.maybe_set_from_fp(str, fp)) return !read_error_flag;
	if (kwmr.maybe_set_from_fp(str, fp)) return !read_error_flag;

	if (str == "city_size_min") {
		if (!read_uint(fp, city_size_min)) {return read_error(str);}
		if (city_size_max == 0) {city_size_max = city_size_min;}
		if (city_size_max < city_size_min) {return read_error(str);}
	}
	else if (str == "city_size_max") {
		if (!read_uint(fp, city_size_max)) {return read_error(str);}
		if (city_size_min == 0) {city_size_min = city_size_max;}
		if (city_size_max < city_size_min) {return read_error(str);}
	}
	else if (str == "make_4_way_ints") {
		if (!read_uint(fp, make_4_way_ints) || make_4_way_ints > 3) {return read_error(str);}
	}
	else if (str == "add_transmission_lines") {
		if (!read_uint(fp, add_tlines) || add_tlines > 2) {return read_error(str);}
	}
	// cars/pedestrian/helicopter models
	else if (str == "car_model") { // multiple car models
		city_model_t car_model;
		if (!car_model.read(fp)) {return read_error(str);}
		if (!car_model.check_filename()) {cerr << "Error: car_model file '" << car_model.fn << "' does not exist; skipping" << endl; return 1;} // nonfatal
		car_model_files.push_back(car_model);
		max_eq(max_car_scale, car_model.scale);
	}
	else if (str == "helicopter_model") { // multiple helicopter models
		city_model_t hc_model;
		if (!hc_model.read(fp, 1)) {return read_error(str);} // is_helicopter=1
		if (!hc_model.check_filename()) {cerr << "Error: helicopter_model file '" << hc_model.fn << "' does not exist; skipping" << endl; return 1;} // nonfatal
		hc_model_files.push_back(hc_model);
	}
	else if (str == "ped_model") {
		city_model_t ped_model;
		if (!ped_model.read(fp, 0, 1)) {return read_error(str);} // is_helicopter=0, is_person=1
		if (!ped_model.check_filename()) {cerr << "Error: ped_model file '" << ped_model.fn << "' does not exist; skipping" << endl; return 1;} // nonfatal
		ped_model.default_anim_name = default_anim_name;
		ped_model_files.push_back(ped_model); // Note: no ped_model_scale
	}
	else if (str == "ped_model_add_anim") { // city ped_model_add_anim <filename> <animation name>
		string const fn(read_quoted_string(fp)), anim_name(read_quoted_string(fp));
		if (ped_model_files.empty()) {cerr << "Error: Can't use ped_model_add_anim without first declaring a ped_model" << endl; return 1;} // nonfatal
		if (!fn.empty()) {ped_model_files.back().anim_fns.emplace_back(fn, anim_name);} // ignore if empty filename, but I guess an empty anim_name is valid?
	}
	else if (str == "default_anim_name") {default_anim_name = read_quoted_string(fp);} // applies to all models loaded after this point
	// other
	else if (str == "smap_size") {
		if (!read_uint(fp, smap_size) || smap_size > 4096) {return read_error(str);}
	}
	else {
		for (unsigned i = 0; i < NUM_OBJ_MODELS; ++i) { // check for object models
			if (str != model_opt_names[i]) continue;
			if (!add_model(i, fp)) {return read_error(str);}
			return 1; // done
		}
		cout << "Unrecognized city keyword in input file: " << str << endl;
		return 0;
	}
	return 1;
}


void buildings_file_err(string const &str, int &error) {
	cout << "Error reading buildings config option " << str << "." << endl;
	error = 1;
}

int building_params_t::read_building_texture(FILE *fp, string const &str, bool is_normal_map, int &error, bool check_filename, bool *no_cracks) {
	char strc[MAX_CHARS] = {0};
	if (!read_str(fp, strc)) {buildings_file_err(str, error);}

	if (check_filename && !check_texture_file_exists(strc)) {
		std::cerr << "Warning: Skipping texture '" << strc << "' that can't be loaded" << endl;
		return -1; // texture filename doesn't exist
	}
	string const name(strc);
	int const ret(get_texture_by_name(name, is_normal_map, tex_inv_y, get_wrap_mir()));
	if (no_cracks != nullptr) {*no_cracks = (ret >= 0 && name.find("carpet") != string::npos);} // carpet floor textures have no cracks
	//cout << "texture filename: " << str << ", ID: " << ret << endl;
	return ret;
}
void building_params_t::read_texture_and_add_if_valid(FILE *fp, string const &str, int &error, vector<unsigned> &tids) {
	// Note: this version doesn't accept numbered texture IDs, but it also doesn't fail on missing files
	int const tid(read_building_texture(fp, str, 0, error, 1)); // is_normal_map=0, check_filename=1
	if (tid >= 0) {tids.push_back(tid);}
}
void read_building_tscale(FILE *fp, tid_nm_pair_t &tex, string const &str, int &error) {
	if (!read_float(fp, tex.tscale_x)) {buildings_file_err(str, error);}
	tex.tscale_y = tex.tscale_x; // uniform
}
void read_building_mat_specular(FILE *fp, string const &str, tid_nm_pair_t &tex, int &error) {
	float mag(0.0), shine(0.0);
	if (read_float(fp, mag) && read_float(fp, shine)) {tex.set_specular(mag, shine);} else {buildings_file_err(str, error);}
}

void building_params_t::init_kw_maps() {
	// global parameters
	kwmb.add("flatten_mesh", flatten_mesh);
	kwmu.add("num_place", num_place);
	kwmu.add("num_tries", num_tries);
	kwmu.add("rand_seed", buildings_rand_seed);
	kwmu.add("max_shadow_maps", max_shadow_maps);
	kwmu.add("max_ext_basement_hall_branches", max_ext_basement_hall_branches);
	kwmu.add("max_ext_basement_room_depth",    max_ext_basement_room_depth);
	kwmu.add("max_room_geom_gen_per_frame",    max_room_geom_gen_per_frame);
	kwmu.add("max_office_basement_floors",     max_office_basement_floors);
	kwmb.add("add_office_backroom_basements",  add_office_br_basements);
	kwmf.add("ao_factor", ao_factor);
	kwmf.add("sec_extra_spacing", sec_extra_spacing);
	kwmf.add("player_coll_radius_scale", player_coll_radius_scale);
	kwmf.add("max_floorplan_window_xscale", max_fp_wind_xscale);
	kwmf.add("max_floorplan_window_yscale", max_fp_wind_yscale);
	kwmf.add("interior_view_dist_scale", interior_view_dist_scale);
	kwmb.add("tt_only", tt_only);
	kwmb.add("infinite_buildings", infinite_buildings);
	kwmb.add("add_secondary_buildings", add_secondary_buildings);
	kwmb.add("cities_all_bldg_mats",    cities_all_bldg_mats);
	kwmb.add("small_city_buildings",    small_city_buildings);
	kwmb.add("add_office_basements",    add_office_basements);
	kwmb.add("put_doors_in_corners",    put_doors_in_corners);
	kwmb.add("add_door_handles",        add_door_handles    );
	kwmb.add("add_basement_tunnels",    add_basement_tunnels);
	kwmr.add("two_floor_retail_prob",   two_floor_retail_prob, FP_CHECK_01);
	kwmr.add("split_prob", cur_mat.split_prob, FP_CHECK_01);
	kwmr.add("cube_prob",  cur_mat.cube_prob,  FP_CHECK_01);
	kwmr.add("round_prob", cur_mat.round_prob, FP_CHECK_01);
	kwmr.add("alt_step_factor_prob", cur_mat.asf_prob,  FP_CHECK_01);
	kwmr.add("glass_floor_alpha",    glass_floor_alpha, FP_CHECK_01);
	// consistency probabilities of houses for cities and blocks
	kwmr.add("house_same_mat_prob",       house_same_mat_prob,       FP_CHECK_01);
	kwmr.add("house_same_size_prob",      house_same_size_prob,      FP_CHECK_01);
	kwmr.add("house_same_geom_prob",      house_same_geom_prob,      FP_CHECK_01);
	kwmr.add("house_same_per_city_prob",  house_same_per_city_prob,  FP_CHECK_01);
	kwmr.add("office_same_mat_prob",      office_same_mat_prob,      FP_CHECK_01);
	kwmr.add("office_same_size_prob",     office_same_size_prob,     FP_CHECK_01);
	kwmr.add("office_same_geom_prob",     office_same_geom_prob,     FP_CHECK_01);
	kwmr.add("office_same_per_city_prob", office_same_per_city_prob, FP_CHECK_01);
	// material parameters
	kwmf.add("place_radius", cur_mat.place_radius);
	kwmf.add("max_delta_z", cur_mat.max_delta_z);
	kwmf.add("min_level_height", cur_mat.min_level_height);
	kwmu.add("min_levels", cur_mat.min_levels);
	kwmu.add("max_levels", cur_mat.max_levels);
	kwmf.add("min_flat_side_amt", cur_mat.min_fsa);
	kwmf.add("max_flat_side_amt", cur_mat.max_fsa);
	kwmf.add("min_alt_step_factor", cur_mat.min_asf);
	kwmf.add("max_alt_step_factor", cur_mat.max_asf);
	kwmf.add("min_altitude",  cur_mat.min_alt);
	kwmf.add("max_altitude",  cur_mat.max_alt);
	kwmf.add("max_rot_angle", cur_mat.max_rot_angle);
	kwmb.add("dome_roof",   dome_roof);
	kwmb.add("onion_roof",  onion_roof);
	kwmb.add("no_city",     cur_mat.no_city);
	kwmb.add("no_walkways", cur_mat.no_walkways);
	kwmr.add("house_prob",  cur_mat.house_prob, FP_CHECK_01);
	kwmr.add("apartment_prob", cur_mat.apartment_prob, FP_CHECK_01);
	// material textures / colors
	kwmb.add("texture_mirror", tex_mirror);
	kwmb.add("texture_inv_y",  tex_inv_y);
	kwmf.add("side_tscale_x",  cur_mat.side_tex.tscale_x);
	kwmf.add("side_tscale_y",  cur_mat.side_tex.tscale_y);
	kwmf.add("side_color_grayscale_rand", cur_mat.side_color.grayscale_rand);
	kwmf.add("roof_color_grayscale_rand", cur_mat.roof_color.grayscale_rand);
	kwmc.add("side_color_min", cur_mat.side_color.cmin);
	kwmc.add("side_color_max", cur_mat.side_color.cmax);
	kwmc.add("roof_color_min", cur_mat.roof_color.cmin);
	kwmc.add("roof_color_max", cur_mat.roof_color.cmax);
	// windows (mostly per-material)
	kwmr.add("window_width",  window_width,  FP_CHECK_01);
	kwmr.add("window_height", window_height, FP_CHECK_01);
	kwmr.add("window_xspace", window_xspace, FP_CHECK_01);
	kwmr.add("window_yspace", window_yspace, FP_CHECK_01);
	kwmr.add("window_xscale", cur_mat.wind_xscale, FP_CHECK_POS);
	kwmr.add("window_yscale", cur_mat.wind_yscale, FP_CHECK_POS);
	kwmf.add("window_xoff",   cur_mat.wind_xoff);
	kwmf.add("window_yoff",   cur_mat.wind_yoff);
	kwmf.add("wall_split_thresh", wall_split_thresh);
	kwmb.add("add_windows",       cur_mat.add_windows);
	kwmb.add("add_window_lights", cur_mat.add_wind_lights);
	kwmc.add("window_color", cur_mat.window_color);
	kwmc.add("wall_color",   cur_mat.wall_color);
	kwmc.add("ceil_color",   cur_mat.ceil_color);
	kwmc.add("floor_color",  cur_mat.floor_color);
	kwmc.add("house_ceil_color",  cur_mat.house_ceil_color);
	kwmc.add("house_floor_color", cur_mat.house_floor_color);
	// people/AI logic
	kwmb.add("enable_people_ai", enable_people_ai);
	kwmu.add("ai_opens_doors",      ai_opens_doors); // 0=don't open doors, 1=only open if player closed door after path selection; 2=always open doors
	kwmb.add("ai_target_player",    ai_target_player);
	kwmb.add("ai_follow_player",    ai_follow_player);
	kwmu.add("ai_player_vis_test",  ai_player_vis_test); // 0=no test, 1=LOS, 2=LOS+FOV, 3=LOS+FOV+lit
	kwmu.add("ai_sees_player_hide", ai_sees_player_hide); // 0=doesn't see the player, 1=sees the player and waits outside the hiding spot, 2=opens the door and comes in
	kwmu.add("people_per_office_min", people_per_office_min);
	kwmu.add("people_per_office_max", people_per_office_max);
	kwmu.add("people_per_house_min",  people_per_house_min);
	kwmu.add("people_per_house_max",  people_per_house_max);
	kwmf.add("ai_retreat_time",       ai_retreat_time);
	kwmr.add("people_min_alpha",      people_min_alpha, FP_CHECK_01);
	kwmu.add("player_model_ix",       player_model_ix);
	kwmb.add("show_player_model",     show_player_model);
	// AI elevators
	kwmb.add("allow_elevator_line",         allow_elevator_line);
	kwmb.add("no_coll_enter_exit_elevator", no_coll_enter_exit_elevator);
	kwmu.add("elevator_capacity",           elevator_capacity);
	kwmf.add("elevator_wait_time",          elevator_wait_time);
	kwmr.add("use_elevator_prob",           use_elevator_prob,         FP_CHECK_01);
	kwmr.add("elevator_wait_recall_prob",   elevator_wait_recall_prob, FP_CHECK_01);
	// rats
	kwmu.add("num_rats_min",    num_rats_min);
	kwmu.add("num_rats_max",    num_rats_max);
	kwmu.add("min_attack_rats", min_attack_rats);
	kwmf.add("rat_speed",       rat_speed);
	kwmr.add("rat_size_min",    rat_size_min, FP_CHECK_POS);
	kwmr.add("rat_size_max",    rat_size_max, FP_CHECK_POS);
	// spiders
	kwmu.add("num_spiders_min", num_spiders_min);
	kwmu.add("num_spiders_max", num_spiders_max);
	kwmf.add("spider_speed",    spider_speed);
	kwmr.add("spider_size_min", spider_size_min, FP_CHECK_POS);
	kwmr.add("spider_size_max", spider_size_max, FP_CHECK_POS);
	kwmr.add("spider_drawer_prob", spider_drawer_prob, FP_CHECK_01);
	// snakes
	kwmu.add("num_snakes_min", num_snakes_min);
	kwmu.add("num_snakes_max", num_snakes_max);
	kwmf.add("snake_speed",    snake_speed);
	kwmr.add("snake_size_min", snake_size_min, FP_CHECK_POS);
	kwmr.add("snake_size_max", snake_size_max, FP_CHECK_POS);
	// insects
	kwmu.add("num_insects_min", num_insects_min);
	kwmu.add("num_insects_max", num_insects_max);
	kwmf.add("insect_speed",    insect_speed);
	kwmr.add("insect_size_min", insect_size_min, FP_CHECK_POS);
	kwmr.add("insect_size_max", insect_size_max, FP_CHECK_POS);
	// gameplay state
	kwmr.add("open_door_prob",       open_door_prob,       FP_CHECK_01);
	kwmr.add("locked_door_prob",     locked_door_prob,     FP_CHECK_01);
	kwmr.add("basement_prob_house",  basement_prob_house,  FP_CHECK_01);
	kwmr.add("basement_prob_office", basement_prob_office, FP_CHECK_01);
	kwmr.add("ball_prob",            ball_prob,            FP_CHECK_01);
	kwmr.add("split_stack_floorplan_prob", split_stack_floorplan_prob, FP_CHECK_01);
	kwmr.add("retail_floorplan_prob",      retail_floorplan_prob,      FP_CHECK_01);
	kwmf.add("player_weight_limit",  player_weight_limit);
	// building water
	kwmr.add("basement_water_level_min", basement_water_level_min); // negative is allowed for no water
	kwmr.add("basement_water_level_max", basement_water_level_max, FP_CHECK_NONNEG); // > 1.0 is allowed for more than one floor of water
	// special commands
	kwmu.add("probability",              cur_prob); // for building materials
	kwmb.add("add_city_interiors",       add_city_interiors);
	kwmb.add("gen_building_interiors",   gen_building_interiors);
	kwmb.add("enable_rotated_room_geom", enable_rotated_room_geom);
	kwmb.add("use_voronoise_cracks",     use_voronoise_cracks);
}

bool building_params_t::parse_buildings_option(FILE *fp) {

	char strc[MAX_CHARS] = {0};
	if (!read_str(fp, strc)) return 0;
	string const str(strc);
	if (kwmb.maybe_set_from_fp(str, fp)) return !read_error;
	if (kwmu.maybe_set_from_fp(str, fp)) return !read_error;
	if (kwmf.maybe_set_from_fp(str, fp)) return !read_error;
	if (kwmc.maybe_set_from_fp(str, fp)) return !read_error;
	if (kwmr.maybe_set_from_fp(str, fp)) return !read_error;

	// material parameters
	if (str == "range_translate") { // x,y only
		if (!(read_float(fp, range_translate.x) && read_float(fp, range_translate.y))) {buildings_file_err(str, read_error);}
	}
	else if (str == "pos_range") {
		if (!read_cube(fp, cur_mat.pos_range, 1)) {buildings_file_err(str, read_error);}
	}
	else if (str == "min_sides") {
		if (!read_uint(fp, cur_mat.min_sides)) {buildings_file_err(str, read_error);}
		if (cur_mat.min_sides < 3) {buildings_file_err(str+" (< 3)", read_error);}
	}
	else if (str == "max_sides") {
		if (!read_uint(fp, cur_mat.max_sides)) {buildings_file_err(str, read_error);}
		if (cur_mat.max_sides < 3) {buildings_file_err(str+" (< 3)", read_error);}
	}
	else if (str == "size_range") {
		if (!read_cube(fp, cur_mat.sz_range)) {buildings_file_err(str, read_error);}
	}
	else if (str == "house_scale_range") { // per-material
		if (!read_float(fp, cur_mat.house_scale_min) || !read_float(fp, cur_mat.house_scale_max)) {buildings_file_err(str, read_error);}
	}
	// material colors / textures
	else if (str == "side_color") {
		if (!read_color(fp, cur_mat.side_color.cmin)) {buildings_file_err(str, read_error);}
		cur_mat.side_color.cmax = cur_mat.side_color.cmin; // same
	}
	else if (str == "roof_color") {
		if (!read_color(fp, cur_mat.roof_color.cmin)) {buildings_file_err(str, read_error);}
		cur_mat.roof_color.cmax = cur_mat.roof_color.cmin; // same
	}
	else if (str == "side_tscale")  {read_building_tscale(fp, cur_mat.side_tex,  str, read_error);} // both X and Y
	else if (str == "roof_tscale")  {read_building_tscale(fp, cur_mat.roof_tex,  str, read_error);} // both X and Y
	else if (str == "wall_tscale")  {read_building_tscale(fp, cur_mat.wall_tex,  str, read_error);} // both X and Y
	else if (str == "ceil_tscale")  {read_building_tscale(fp, cur_mat.ceil_tex,  str, read_error);} // both X and Y
	else if (str == "floor_tscale") {read_building_tscale(fp, cur_mat.floor_tex, str, read_error);} // both X and Y
	else if (str == "house_ceil_tscale")  {read_building_tscale(fp, cur_mat.house_ceil_tex,  str, read_error);} // both X and Y
	else if (str == "house_floor_tscale") {read_building_tscale(fp, cur_mat.house_floor_tex, str, read_error);} // both X and Y
	else if (str == "basement_floor_tscale") {read_building_tscale(fp, cur_mat.basement_floor_tex, str, read_error);} // both X and Y
	// building textures
	// Warning: setting options such as tex_inv_y for textures that have already been loaded will have no effect!
	else if (str == "side_tid"    ) {cur_mat.side_tex.tid     = read_building_texture(fp, str, 0, read_error);}
	else if (str == "side_nm_tid" ) {cur_mat.side_tex.nm_tid  = read_building_texture(fp, str, 1, read_error);}
	else if (str == "roof_tid"    ) {cur_mat.roof_tex.tid     = read_building_texture(fp, str, 0, read_error);}
	else if (str == "roof_nm_tid" ) {cur_mat.roof_tex.nm_tid  = read_building_texture(fp, str, 1, read_error);}
	// interiors
	else if (str == "wall_tid"    ) {cur_mat.wall_tex.tid     = read_building_texture(fp, str, 0, read_error);}
	else if (str == "wall_nm_tid" ) {cur_mat.wall_tex.nm_tid  = read_building_texture(fp, str, 1, read_error);}
	else if (str == "floor_tid"   ) {cur_mat.floor_tex.tid    = read_building_texture(fp, str, 0, read_error, 0, &cur_mat.floor_tex.no_cracks);}
	else if (str == "floor_nm_tid") {cur_mat.floor_tex.nm_tid = read_building_texture(fp, str, 1, read_error);}
	else if (str == "ceil_tid"    ) {cur_mat.ceil_tex.tid     = read_building_texture(fp, str, 0, read_error);}
	else if (str == "ceil_nm_tid" ) {cur_mat.ceil_tex.nm_tid  = read_building_texture(fp, str, 1, read_error);}
	else if (str == "house_floor_tid"   ) {cur_mat.house_floor_tex.tid    = read_building_texture(fp, str, 0, read_error, 0, &cur_mat.house_floor_tex.no_cracks);}
	else if (str == "house_floor_nm_tid") {cur_mat.house_floor_tex.nm_tid = read_building_texture(fp, str, 1, read_error);}
	else if (str == "house_ceil_tid"    ) {cur_mat.house_ceil_tex.tid     = read_building_texture(fp, str, 0, read_error);}
	else if (str == "house_ceil_nm_tid" ) {cur_mat.house_ceil_tex.nm_tid  = read_building_texture(fp, str, 1, read_error);}
	else if (str == "basement_floor_tid") {cur_mat.basement_floor_tex.tid = read_building_texture(fp, str, 0, read_error, 0, &cur_mat.basement_floor_tex.no_cracks);}
	else if (str == "basement_floor_nm_tid") {cur_mat.basement_floor_tex.nm_tid = read_building_texture(fp, str, 1, read_error);}
	// specular
	else if (str == "side_specular" ) {read_building_mat_specular(fp, str, cur_mat.side_tex,  read_error);}
	else if (str == "roof_specular" ) {read_building_mat_specular(fp, str, cur_mat.roof_tex,  read_error);}
	else if (str == "wall_specular" ) {read_building_mat_specular(fp, str, cur_mat.wall_tex,  read_error);}
	else if (str == "ceil_specular" ) {read_building_mat_specular(fp, str, cur_mat.ceil_tex,  read_error);}
	else if (str == "floor_specular") {read_building_mat_specular(fp, str, cur_mat.floor_tex, read_error);}
	else if (str == "house_ceil_specular" ) {read_building_mat_specular(fp, str, cur_mat.house_ceil_tex,  read_error);}
	else if (str == "house_floor_specular") {read_building_mat_specular(fp, str, cur_mat.house_floor_tex, read_error);}
	// room objects/textures
	else if (str == "add_rug_texture"    ) {read_texture_and_add_if_valid(fp, str, read_error, rug_tids    );}
	else if (str == "add_picture_texture") {read_texture_and_add_if_valid(fp, str, read_error, picture_tids);}
	else if (str == "add_desktop_texture") {read_texture_and_add_if_valid(fp, str, read_error, desktop_tids);}
	else if (str == "add_sheet_texture"  ) {read_texture_and_add_if_valid(fp, str, read_error, sheet_tids  );}
	else if (str == "add_paper_texture"  ) {read_texture_and_add_if_valid(fp, str, read_error, paper_tids  );}
	else if (str == "add_metal_texture"  ) {read_texture_and_add_if_valid(fp, str, read_error, metal_tids  );}
	else if (str == "add_flag_texture"   ) {read_texture_and_add_if_valid(fp, str, read_error, flag_tids   );} // not a room object, but fits with this loading system
	else if (str == "add_food_box_texture") {
		read_texture_and_add_if_valid(fp, str, read_error, food_box_tids);
		string const food_name(read_quoted_string(fp));
		if (food_name.empty()) {buildings_file_err(str, read_error);}
		food_box_names.push_back(food_name);
	}
	// special commands
	else if (str == "add_material") {add_cur_mat();}
	else {
		cout << "Unrecognized buildings keyword in input file: " << str << endl;
		read_error = 1;
	}
	return !read_error;
}

void building_params_t::add_cur_mat() {
	unsigned const mat_ix(materials.size());

	for (unsigned n = 0; n < cur_prob; ++n) { // add more references to this mat for higher probability
		mat_gen_ix.push_back(mat_ix);
		if (cities_all_bldg_mats || ((!cur_mat.no_city) ^ small_city_buildings)) {mat_gen_ix_city.push_back(mat_ix);}
		if (cur_mat.no_city) {mat_gen_ix_nocity.push_back(mat_ix);}
		if (cur_mat.house_prob > 0.0) {mat_gen_ix_res.push_back(mat_ix);}
	}
	materials.push_back(cur_mat);
	materials.back().finalize();
	materials.back().update_range(range_translate);
	has_normal_map |= cur_mat.has_normal_map();
}
vector<unsigned> const &building_params_t::get_mat_list(bool city_only, bool non_city_only, bool residential) const {
	if (residential  ) return mat_gen_ix_res;
	if (city_only    ) return mat_gen_ix_city;
	if (non_city_only) return mat_gen_ix_nocity;
	return mat_gen_ix; // all materials
}
unsigned building_params_t::choose_rand_mat(rand_gen_t &rgen, bool city_only, bool non_city_only, bool residential) const {
	vector<unsigned> const &mat_ix_list(get_mat_list(city_only, non_city_only, residential));
	assert(!mat_ix_list.empty());
	return mat_ix_list[rgen.rand()%mat_ix_list.size()];
}
float building_params_t::get_max_house_size() const {
	float max_sz(0.0);

	for (auto m = materials.begin(); m != materials.end(); ++m) {
		if (m->house_prob > 0.0) {max_eq(max_sz, m->house_scale_max*max(m->sz_range.x2(), m->sz_range.y2()));} // take max of x/y size upper bounds
	}
	return max_sz;
}
void building_params_t::set_pos_range(cube_t const &pos_range) {
	//cout << "pos_range: " << pos_range.str() << endl;
	cur_mat.set_pos_range(pos_range);
	for (auto i = materials.begin(); i != materials.end(); ++i) {i->set_pos_range(pos_range);}
}
void building_params_t::restore_prev_pos_range() {
	cur_mat.restore_prev_pos_range();
	for (auto i = materials.begin(); i != materials.end(); ++i) {i->restore_prev_pos_range();}
}
void building_params_t::finalize() {
	//if (materials.empty()) {add_cur_mat();} // add current (maybe default) material - seems not to be needed
}

void building_mat_t::update_range(vector3d const &range_translate) {
	if (place_radius > 0.0) { // clip range to place_radius
		point const center(pos_range.get_cube_center());

		for (unsigned d = 0; d < 2; ++d) { // x,y
			max_eq(pos_range.d[d][0], (center[d] - place_radius));
			min_eq(pos_range.d[d][1], (center[d] + place_radius));
		}
	}
	pos_range += range_translate;
}

void color_range_t::gen_color(colorRGBA &color, rand_gen_t &rgen) const {
	if (cmin == cmax) {color = cmin;} // single exact color
	else {UNROLL_4X(color[i_] = rgen.rand_uniform(cmin[i_], cmax[i_]);)}
	if (grayscale_rand > 0.0) {
		float const v(grayscale_rand*rgen.rand_float());
		UNROLL_3X(color[i_] += v;)
	}
}

// windows are scaled to make the texture look correct; this is fine for exterior building wall materials that have no windows, since we can place the windows however we want;
// but some office buildings have windows spaced too close together, and we don't have control over it here;
// so instead we use a larger window space for floorplanning, since windows aren't cut out of the walls anyway
float building_mat_t::get_window_tx() const {return wind_xscale*global_building_params.get_window_tx();}
float building_mat_t::get_window_ty() const {return wind_yscale*global_building_params.get_window_ty();}

void building_mat_t::finalize() { // compute and cache spacing values
	if (!global_building_params.windows_enabled()) return; // don't set the variables below
	float tx(get_window_tx()), ty(get_window_ty());
	if (global_building_params.max_fp_wind_yscale > 0.0) {min_eq(ty, global_building_params.max_fp_wind_yscale*global_building_params.get_window_ty());}
	if (global_building_params.max_fp_wind_xscale > 0.0) {min_eq(tx, global_building_params.max_fp_wind_xscale*global_building_params.get_window_tx());}
	floor_spacing         = 1.0/(2.0*ty);
	floorplan_wind_xscale = 2.0*tx;
}

