// 3D World - City and Building Parameter Parsing
// by Frank Gennari
// 09/04/20

#include "city.h"
#include "buildings.h"
#include "file_utils.h"


extern building_params_t global_building_params;


string const model_opt_names[NUM_OBJ_MODELS] =
{"toilet_model", "sink_model", "tub_model", "fridge_model", "stove_model", "tv_model", ""/*monitor*/, "couch_model", "office_chair_model", "urinal_model",
"lamp_model", "washer_model", "dryer_model", "key_model", "hanger_model", "clothing_model", "fire_escape_model", "cup_model",
"fire_hydrant_model", "substation_model", "umbrella_model"};

bool city_params_t::read_option(FILE *fp) {

	char strc[MAX_CHARS] = {0};
	if (!read_str(fp, strc)) return 0;
	string const str(strc);

	if (str == "num_cities") {
		if (!read_uint(fp, num_cities)) {return read_error(str);}
	}
	else if (str == "num_rr_tracks") {
		if (!read_uint(fp, num_rr_tracks)) {return read_error(str);}
	}
	else if (str == "num_samples") {
		if (!read_uint(fp, num_samples) || num_samples == 0) {return read_error(str);}
	}
	else if (str == "num_conn_tries") {
		if (!read_uint(fp, num_conn_tries) || num_conn_tries == 0) {return read_error(str);}
	}
	else if (str == "plots_to_parks_ratio") {
		if (!read_uint(fp, park_rate)) {return read_error(str);}
	}
	else if (str == "city_size_min") {
		if (!read_uint(fp, city_size_min)) {return read_error(str);}
		if (city_size_max == 0) {city_size_max = city_size_min;}
		if (city_size_max < city_size_min) {return read_error(str);}
	}
	else if (str == "city_size_max") {
		if (!read_uint(fp, city_size_max)) {return read_error(str);}
		if (city_size_min == 0) {city_size_min = city_size_max;}
		if (city_size_max < city_size_min) {return read_error(str);}
	}
	else if (str == "city_border") {
		if (!read_uint(fp, city_border)) {return read_error(str);}
	}
	else if (str == "road_border") {
		if (!read_uint(fp, road_border)) {return read_error(str);}
	}
	else if (str == "slope_width") {
		if (!read_uint(fp, slope_width)) {return read_error(str);}
	}
	else if (str == "road_width") {
		if (!read_float(fp, road_width) || road_width < 0.0) {return read_error(str);}
	}
	else if (str == "road_spacing") {
		if (!read_float(fp, road_spacing) || road_spacing < 0.0) {return read_error(str);}
	}
	else if (str == "road_spacing_rand") {
		if (!read_float(fp, road_spacing_rand) || road_spacing_rand < 0.0) {return read_error(str);}
	}
	else if (str == "road_spacing_xy_add") {
		if (!read_float(fp, road_spacing_xy_add) || road_spacing_xy_add < 0.0) {return read_error(str);}
	}
	else if (str == "conn_road_seg_len") {
		if (!read_float(fp, conn_road_seg_len) || conn_road_seg_len <= 0.0) {return read_error(str);}
	}
	else if (str == "max_road_slope") {
		if (!read_float(fp, max_road_slope) || max_road_slope <= 0.0) {return read_error(str);}
	}
	else if (str == "max_track_slope") {
		if (!read_float(fp, max_track_slope) || max_track_slope <= 0.0) {return read_error(str);}
	}
	else if (str == "make_4_way_ints") {
		if (!read_uint(fp, make_4_way_ints) || make_4_way_ints > 3) {return read_error(str);}
	}
	else if (str == "add_transmission_lines") {
		if (!read_uint(fp, add_tlines) || add_tlines > 2) {return read_error(str);}
	}
	else if (str == "residential_probability") {
		if (!read_float(fp, residential_probability)) {return read_error(str);}
	}
	else if (str == "assign_house_plots") {
		if (!read_bool(fp, assign_house_plots)) {return read_error(str);}
	}
	else if (str == "new_city_conn_road_alg") {
		if (!read_bool(fp, new_city_conn_road_alg)) {return read_error(str);}
	}
	// cars
	else if (str == "num_cars") {
		if (!read_uint(fp, num_cars)) {return read_error(str);}
	}
	else if (str == "car_speed") {
		if (!read_float(fp, car_speed) || car_speed < 0.0) {return read_error(str);}
	}
	else if (str == "traffic_balance_val") {
		if (!read_float(fp, traffic_balance_val) || traffic_balance_val > 1.0 || traffic_balance_val < 0.0) {return read_error(str);}
	}
	else if (str == "new_city_prob") {
		if (!read_float(fp, new_city_prob) || new_city_prob < 0.0) {return read_error(str);}
	}
	else if (str == "enable_car_path_finding") {
		if (!read_bool(fp, enable_car_path_finding)) {return read_error(str);}
	}
	else if (str == "convert_model_files") {
		if (!read_bool(fp, convert_model_files)) {return read_error(str);}
	}
	else if (str == "cars_use_driveways") {
		if (!read_bool(fp, cars_use_driveways)) {return read_error(str);}
	}
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
	// pedestrians
	else if (str == "num_peds") {
		if (!read_uint(fp, num_peds)) {return read_error(str);}
	}
	else if (str == "num_building_peds") {
		if (!read_uint(fp, num_building_peds)) {return read_error(str);}
	}
	else if (str == "ped_speed") {
		if (!read_float(fp, ped_speed) || ped_speed < 0.0) {return read_error(str);}
	}
	else if (str == "ped_model") {
		city_model_t ped_model;
		if (!ped_model.read(fp)) {return read_error(str);}
		if (!ped_model.check_filename()) {cerr << "Error: ped_model file '" << ped_model.fn << "' does not exist; skipping" << endl; return 1;} // nonfatal
		ped_model_files.push_back(ped_model); // Note: no ped_model_scale
	}
	else if (str == "ped_respawn_at_dest") {
		if (!read_bool(fp, ped_respawn_at_dest)) {return read_error(str);}
	}
	// parking lots
	else if (str == "min_park_spaces") { // with default road parameters, can be up to 28
		if (!read_uint(fp, min_park_spaces)) {return read_error(str);}
	}
	else if (str == "min_park_rows") { // with default road parameters, can be up to 8
		if (!read_uint(fp, min_park_rows)) {return read_error(str);}
	}
	else if (str == "min_park_density") {
		if (!read_float(fp, min_park_density)) {return read_error(str);}
	}
	else if (str == "max_park_density") {
		if (!read_float(fp, max_park_density) || max_park_density < 0.0) {return read_error(str);}
	}
	// lighting
	else if (str == "max_lights") {
		if (!read_uint(fp, max_lights)) {return read_error(str);}
	}
	else if (str == "max_shadow_maps") {
		if (!read_uint(fp, max_shadow_maps)) {return read_error(str);}
	}
	else if (str == "smap_size") {
		if (!read_uint(fp, smap_size) || smap_size > 4096) {return read_error(str);}
	}
	else if (str == "car_shadows") {
		if (!read_bool(fp, car_shadows)) {return read_error(str);}
	}
	// trees
	else if (str == "max_trees_per_plot") {
		if (!read_uint(fp, max_trees_per_plot)) {return read_error(str);}
	}
	else if (str == "tree_spacing") {
		if (!read_float(fp, tree_spacing) || tree_spacing <= 0.0) {return read_error(str);}
	}
	// detail objects
	else if (str == "max_benches_per_plot") {
		if (!read_uint(fp, max_benches_per_plot)) {return read_error(str);}
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
bool check_01(float v) {return (v >= 0.0 && v <= 1.0);}

bool check_texture_file_exists(string const &filename);

int read_building_texture(FILE *fp, string const &str, int &error, bool check_filename=0) {
	char strc[MAX_CHARS] = {0};
	if (!read_str(fp, strc)) {buildings_file_err(str, error);}

	if (check_filename && !check_texture_file_exists(strc)) {
		std::cerr << "Warning: Skipping texture '" << strc << "' that can't be loaded" << endl;
		return -1; // texture filename doesn't exist
	}
	int const ret(get_texture_by_name(std::string(strc), 0, global_building_params.tex_inv_y, global_building_params.get_wrap_mir()));
	//cout << "texture filename: " << str << ", ID: " << ret << endl;
	return ret;
}
void read_texture_and_add_if_valid(FILE *fp, string const &str, int &error, vector<unsigned> &tids) {
	// Note: this version doesn't accept numbered texture IDs, but it also doesn't fail on missing files
	int const tid(read_building_texture(fp, str, error, 1)); // check_filename=1
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

bool parse_buildings_option(FILE *fp) {

	int error(0);
	char strc[MAX_CHARS] = {0};
	if (!read_str(fp, strc)) return 0;
	string const str(strc);
	bool unmatched(0);

	// global parameters
	if (str == "flatten_mesh") {
		if (!read_bool(fp, global_building_params.flatten_mesh)) {buildings_file_err(str, error);}
	}
	else if (str == "num_place") {
		if (!read_uint(fp, global_building_params.num_place)) {buildings_file_err(str, error);}
	}
	else if (str == "num_tries") {
		if (!read_uint(fp, global_building_params.num_tries)) {buildings_file_err(str, error);}
	}
	else if (str == "rand_seed") {
		if (!read_uint(fp, global_building_params.buildings_rand_seed)) {buildings_file_err(str, error);}
	}
	else if (str == "max_shadow_maps") {
		if (!read_uint(fp, global_building_params.max_shadow_maps)) {buildings_file_err(str, error);}
	}
	else if (str == "ao_factor") {
		if (!read_zero_one_float(fp, global_building_params.ao_factor)) {buildings_file_err(str, error);}
	}
	else if (str == "sec_extra_spacing") {
		if (!read_float(fp, global_building_params.sec_extra_spacing)) {buildings_file_err(str, error);}
	}
	else if (str == "player_coll_radius_scale") {
		if (!read_float(fp, global_building_params.player_coll_radius_scale)) {buildings_file_err(str, error);}
	}
	else if (str == "max_floorplan_window_xscale") {
		if (!read_float(fp, global_building_params.max_fp_wind_xscale)) {buildings_file_err(str, error);}
	}
	else if (str == "max_floorplan_window_yscale") {
		if (!read_float(fp, global_building_params.max_fp_wind_yscale)) {buildings_file_err(str, error);}
	}
	else if (str == "interior_view_dist_scale") {
		if (!read_float(fp, global_building_params.interior_view_dist_scale)) {buildings_file_err(str, error);}
	}
	else if (str == "tt_only") {
		if (!read_bool(fp, global_building_params.tt_only)) {buildings_file_err(str, error);}
	}
	else if (str == "infinite_buildings") {
		if (!read_bool(fp, global_building_params.infinite_buildings)) {buildings_file_err(str, error);}
	}
	else if (str == "add_secondary_buildings") {
		if (!read_bool(fp, global_building_params.add_secondary_buildings)) {buildings_file_err(str, error);}
	}
	else if (str == "add_office_basements") {
		if (!read_bool(fp, global_building_params.add_office_basements)) {buildings_file_err(str, error);}
	}
	else if (str == "enable_people_ai") {
		if (!read_bool(fp, global_building_params.enable_people_ai)) {buildings_file_err(str, error);}
	}
	// material parameters
	else if (str == "range_translate") { // x,y only
		if (!(read_float(fp, global_building_params.range_translate.x) &&
			read_float(fp, global_building_params.range_translate.y))) {buildings_file_err(str, error);}
	}
	else if (str == "pos_range") {
		if (!read_cube(fp, global_building_params.cur_mat.pos_range, 1)) {buildings_file_err(str, error);}
	}
	else if (str == "place_radius") {
		if (!read_float(fp, global_building_params.cur_mat.place_radius)) {buildings_file_err(str, error);}
	}
	else if (str == "max_delta_z") {
		if (!read_float(fp, global_building_params.cur_mat.max_delta_z)) {buildings_file_err(str, error);}
	}
	else if (str == "min_level_height") {
		if (!read_float(fp, global_building_params.cur_mat.min_level_height)) {buildings_file_err(str, error);}
	}
	else if (str == "max_rot_angle") {
		if (!read_float(fp, global_building_params.cur_mat.max_rot_angle)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.max_rot_angle *= TO_RADIANS; // specified in degrees, stored in radians
	}
	else if (str == "split_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.split_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "cube_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.cube_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "round_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.round_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "alt_step_factor_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.asf_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "min_levels") {
		if (!read_uint(fp, global_building_params.cur_mat.min_levels)) {buildings_file_err(str, error);}
	}
	else if (str == "max_levels") {
		if (!read_uint(fp, global_building_params.cur_mat.max_levels)) {buildings_file_err(str, error);}
	}
	else if (str == "min_sides") {
		if (!read_uint(fp, global_building_params.cur_mat.min_sides)) {buildings_file_err(str, error);}
		if (global_building_params.cur_mat.min_sides < 3) {buildings_file_err(str+" (< 3)", error);}
	}
	else if (str == "max_sides") {
		if (!read_uint(fp, global_building_params.cur_mat.max_sides)) {buildings_file_err(str, error);}
		if (global_building_params.cur_mat.max_sides < 3) {buildings_file_err(str+" (< 3)", error);}
	}
	else if (str == "min_flat_side_amt") {
		if (!read_float(fp, global_building_params.cur_mat.min_fsa)) {buildings_file_err(str, error);}
	}
	else if (str == "max_flat_side_amt") {
		if (!read_float(fp, global_building_params.cur_mat.max_fsa)) {buildings_file_err(str, error);}
	}
	else if (str == "min_alt_step_factor") {
		if (!read_float(fp, global_building_params.cur_mat.min_asf)) {buildings_file_err(str, error);}
	}
	else if (str == "max_alt_step_factor") {
		if (!read_float(fp, global_building_params.cur_mat.max_asf)) {buildings_file_err(str, error);}
	}
	else if (str == "size_range") {
		if (!read_cube(fp, global_building_params.cur_mat.sz_range)) {buildings_file_err(str, error);}
	}
	else if (str == "min_altitude") {
		if (!read_float(fp, global_building_params.cur_mat.min_alt)) {buildings_file_err(str, error);}
	}
	else if (str == "max_altitude") {
		if (!read_float(fp, global_building_params.cur_mat.max_alt)) {buildings_file_err(str, error);}
	}
	else if (str == "dome_roof") {
		if (!read_bool(fp, global_building_params.dome_roof)) {buildings_file_err(str, error);}
	}
	else if (str == "onion_roof") {
		if (!read_bool(fp, global_building_params.onion_roof)) {buildings_file_err(str, error);}
	}
	else if (str == "no_city") {
		if (!read_bool(fp, global_building_params.cur_mat.no_city)) {buildings_file_err(str, error);}
	}
	// material textures
	else if (str == "texture_mirror") {
		if (!read_bool(fp, global_building_params.tex_mirror)) {buildings_file_err(str, error);}
	}
	else if (str == "texture_inv_y") {
		if (!read_bool(fp, global_building_params.tex_inv_y)) {buildings_file_err(str, error);}
	}
	else if (str == "side_tscale") {read_building_tscale(fp, global_building_params.cur_mat.side_tex, str, error);} // both X and Y
	else if (str == "side_tscale_x") {
		if (!read_float(fp, global_building_params.cur_mat.side_tex.tscale_x)) {buildings_file_err(str, error);}
	}
	else if (str == "side_tscale_y") {
		if (!read_float(fp, global_building_params.cur_mat.side_tex.tscale_y)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_tscale")  {read_building_tscale(fp, global_building_params.cur_mat.roof_tex,  str, error);} // both X and Y
	else if (str == "wall_tscale")  {read_building_tscale(fp, global_building_params.cur_mat.wall_tex,  str, error);} // both X and Y
	else if (str == "ceil_tscale")  {read_building_tscale(fp, global_building_params.cur_mat.ceil_tex,  str, error);} // both X and Y
	else if (str == "floor_tscale") {read_building_tscale(fp, global_building_params.cur_mat.floor_tex, str, error);} // both X and Y
	else if (str == "house_ceil_tscale")  {read_building_tscale(fp, global_building_params.cur_mat.house_ceil_tex,  str, error);} // both X and Y
	else if (str == "house_floor_tscale") {read_building_tscale(fp, global_building_params.cur_mat.house_floor_tex, str, error);} // both X and Y
	else if (str == "basement_floor_tscale") {read_building_tscale(fp, global_building_params.cur_mat.basement_floor_tex, str, error);} // both X and Y
	// building textures
	// Warning: setting options such as tex_inv_y for textures that have already been loaded will have no effect!
	else if (str == "side_tid"    ) {global_building_params.cur_mat.side_tex.tid     = read_building_texture(fp, str, error);}
	else if (str == "side_nm_tid" ) {global_building_params.cur_mat.side_tex.nm_tid  = read_building_texture(fp, str, error);}
	else if (str == "roof_tid"    ) {global_building_params.cur_mat.roof_tex.tid     = read_building_texture(fp, str, error);}
	else if (str == "roof_nm_tid" ) {global_building_params.cur_mat.roof_tex.nm_tid  = read_building_texture(fp, str, error);}
	// interiors
	else if (str == "wall_tid"    ) {global_building_params.cur_mat.wall_tex.tid     = read_building_texture(fp, str, error);}
	else if (str == "wall_nm_tid" ) {global_building_params.cur_mat.wall_tex.nm_tid  = read_building_texture(fp, str, error);}
	else if (str == "floor_tid"   ) {global_building_params.cur_mat.floor_tex.tid    = read_building_texture(fp, str, error);}
	else if (str == "floor_nm_tid") {global_building_params.cur_mat.floor_tex.nm_tid = read_building_texture(fp, str, error);}
	else if (str == "ceil_tid"    ) {global_building_params.cur_mat.ceil_tex.tid     = read_building_texture(fp, str, error);}
	else if (str == "ceil_nm_tid" ) {global_building_params.cur_mat.ceil_tex.nm_tid  = read_building_texture(fp, str, error);}
	else if (str == "house_floor_tid"   ) {global_building_params.cur_mat.house_floor_tex.tid    = read_building_texture(fp, str, error);}
	else if (str == "house_floor_nm_tid") {global_building_params.cur_mat.house_floor_tex.nm_tid = read_building_texture(fp, str, error);}
	else if (str == "house_ceil_tid"    ) {global_building_params.cur_mat.house_ceil_tex.tid     = read_building_texture(fp, str, error);}
	else if (str == "house_ceil_nm_tid" ) {global_building_params.cur_mat.house_ceil_tex.nm_tid  = read_building_texture(fp, str, error);}
	else if (str == "basement_floor_tid"   ) {global_building_params.cur_mat.basement_floor_tex.tid    = read_building_texture(fp, str, error);}
	else if (str == "basement_floor_nm_tid") {global_building_params.cur_mat.basement_floor_tex.nm_tid = read_building_texture(fp, str, error);}
	else if (str == "open_door_prob") {
		if (!read_float(fp, global_building_params.open_door_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "locked_door_prob") {
	if (!read_float(fp, global_building_params.locked_door_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "basement_prob") {
		if (!read_float(fp, global_building_params.basement_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "ball_prob") {
		if (!read_float(fp, global_building_params.ball_prob)) {buildings_file_err(str, error);}
	}
	// material colors
	else if (str == "side_color") {
		if (!read_color(fp, global_building_params.cur_mat.side_color.cmin)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.side_color.cmax = global_building_params.cur_mat.side_color.cmin; // same
	}
	else if (str == "side_color_min") {
		if (!read_color(fp, global_building_params.cur_mat.side_color.cmin)) {buildings_file_err(str, error);}
	}
	else if (str == "side_color_max") {
		if (!read_color(fp, global_building_params.cur_mat.side_color.cmax)) {buildings_file_err(str, error);}
	}
	else if (str == "side_color_grayscale_rand") {
		if (!read_float(fp, global_building_params.cur_mat.side_color.grayscale_rand)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_color") {
		if (!read_color(fp, global_building_params.cur_mat.roof_color.cmin)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.roof_color.cmax = global_building_params.cur_mat.roof_color.cmin; // same
	}
	else if (str == "roof_color_min") {
		if (!read_color(fp, global_building_params.cur_mat.roof_color.cmin)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_color_max") {
		if (!read_color(fp, global_building_params.cur_mat.roof_color.cmax)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_color_grayscale_rand") {
		if (!read_float(fp, global_building_params.cur_mat.roof_color.grayscale_rand)) {buildings_file_err(str, error);}
	}
	// specular
	else if (str == "side_specular" ) {read_building_mat_specular(fp, str, global_building_params.cur_mat.side_tex,  error);}
	else if (str == "roof_specular" ) {read_building_mat_specular(fp, str, global_building_params.cur_mat.roof_tex,  error);}
	else if (str == "wall_specular" ) {read_building_mat_specular(fp, str, global_building_params.cur_mat.wall_tex,  error);}
	else if (str == "ceil_specular" ) {read_building_mat_specular(fp, str, global_building_params.cur_mat.ceil_tex,  error);}
	else if (str == "floor_specular") {read_building_mat_specular(fp, str, global_building_params.cur_mat.floor_tex, error);}
	else if (str == "house_ceil_specular" ) {read_building_mat_specular(fp, str, global_building_params.cur_mat.house_ceil_tex,  error);}
	else if (str == "house_floor_specular") {read_building_mat_specular(fp, str, global_building_params.cur_mat.house_floor_tex, error);}
	// windows
	// Note: this should be an else-if, but I have to split this if/else tree due to a MSVS compiler limit
	else {unmatched = 1;}
	/*else*/ if (str == "window_width") {
		if (!read_float(fp, global_building_params.window_width) || !check_01(global_building_params.window_width)) {buildings_file_err(str, error);}
	}
	else if (str == "window_height") {
		if (!read_float(fp, global_building_params.window_height) || !check_01(global_building_params.window_height)) {buildings_file_err(str, error);}
	}
	else if (str == "window_xspace") {
		if (!read_float(fp, global_building_params.window_xspace) || !check_01(global_building_params.window_xspace)) {buildings_file_err(str, error);}
	}
	else if (str == "window_yspace") {
		if (!read_float(fp, global_building_params.window_yspace) || !check_01(global_building_params.window_yspace)) {buildings_file_err(str, error);}
	}
	else if (str == "window_xscale") {
		if (!read_float(fp, global_building_params.cur_mat.wind_xscale) || global_building_params.cur_mat.wind_xscale < 0.0) {buildings_file_err(str, error);}
	}
	else if (str == "window_yscale") {
		if (!read_float(fp, global_building_params.cur_mat.wind_yscale) || global_building_params.cur_mat.wind_yscale < 0.0) {buildings_file_err(str, error);}
	}
	else if (str == "window_xoff") {
		if (!read_float(fp, global_building_params.cur_mat.wind_xoff)) {buildings_file_err(str, error);}
	}
	else if (str == "window_yoff") {
		if (!read_float(fp, global_building_params.cur_mat.wind_yoff)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.wind_yoff *= -1.0; // invert Y
	}
	else if (str == "wall_split_thresh") {
		if (!read_float(fp, global_building_params.wall_split_thresh)) {buildings_file_err(str, error);}
	}
	else if (str == "add_windows") { // per-material
		if (!read_bool(fp, global_building_params.cur_mat.add_windows)) {buildings_file_err(str, error);}
	}
	else if (str == "add_window_lights") { // per-material
		if (!read_bool(fp, global_building_params.cur_mat.add_wind_lights)) {buildings_file_err(str, error);}
	}
	else if (str == "house_prob") { // per-material
		if (!read_float(fp, global_building_params.cur_mat.house_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "house_scale_range") { // per-material
		if (!read_float(fp, global_building_params.cur_mat.house_scale_min) || !read_float(fp, global_building_params.cur_mat.house_scale_max)) {buildings_file_err(str, error);}
	}
	else if (str == "window_color") { // per-material
		if (!read_color(fp, global_building_params.cur_mat.window_color)) {buildings_file_err(str, error);}
	}
	else if (str == "wall_color") { // per-material
		if (!read_color(fp, global_building_params.cur_mat.wall_color)) {buildings_file_err(str, error);}
	}
	else if (str == "ceil_color") { // per-material
		if (!read_color(fp, global_building_params.cur_mat.ceil_color)) {buildings_file_err(str, error);}
	}
	else if (str == "floor_color") { // per-material
		if (!read_color(fp, global_building_params.cur_mat.floor_color)) {buildings_file_err(str, error);}
	}
	else if (str == "house_ceil_color") { // per-material
		if (!read_color(fp, global_building_params.cur_mat.house_ceil_color)) {buildings_file_err(str, error);}
	}
	else if (str == "house_floor_color") { // per-material
		if (!read_color(fp, global_building_params.cur_mat.house_floor_color)) {buildings_file_err(str, error);}
	}
	// room objects/textures
	else if (str == "add_rug_texture"    ) {read_texture_and_add_if_valid(fp, str, error, global_building_params.rug_tids    );}
	else if (str == "add_picture_texture") {read_texture_and_add_if_valid(fp, str, error, global_building_params.picture_tids);}
	else if (str == "add_desktop_texture") {read_texture_and_add_if_valid(fp, str, error, global_building_params.desktop_tids);}
	else if (str == "add_sheet_texture"  ) {read_texture_and_add_if_valid(fp, str, error, global_building_params.sheet_tids  );}
	else if (str == "add_paper_texture"  ) {read_texture_and_add_if_valid(fp, str, error, global_building_params.paper_tids  );}
	// AI logic
	else if (str == "ai_opens_doors") { // 0=don't open doors, 1=only open if player closed door after path selection; 2=always open doors
		if (!read_uint(fp, global_building_params.ai_opens_doors)) {buildings_file_err(str, error);}
	}
	else if (str == "ai_target_player") {
		if (!read_bool(fp, global_building_params.ai_target_player)) {buildings_file_err(str, error);}
	}
	else if (str == "ai_follow_player") {
		if (!read_bool(fp, global_building_params.ai_follow_player)) {buildings_file_err(str, error);}
	}
	else if (str == "ai_player_vis_test") { // 0=no test, 1=LOS, 2=LOS+FOV, 3=LOS+FOV+lit
		if (!read_uint(fp, global_building_params.ai_player_vis_test)) {buildings_file_err(str, error);}
	}
	// animals
	else if (str == "num_rats") {
		if (!read_uint(fp, global_building_params.num_rats)) {buildings_file_err(str, error);}
	}
	// gameplay state
	else if (str == "player_weight_limit") {
		if (!read_float(fp, global_building_params.player_weight_limit)) {buildings_file_err(str, error);}
	}
	// special commands
	else if (str == "probability") {
		if (!read_uint(fp, global_building_params.cur_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "add_material") {global_building_params.add_cur_mat();}
	else if (str == "add_city_interiors") {
		if (!read_bool(fp, global_building_params.add_city_interiors)) {buildings_file_err(str, error);}
	}
	else if (str == "gen_building_interiors") {
		if (!read_bool(fp, global_building_params.gen_building_interiors)) {buildings_file_err(str, error);}
	}
	else if (str == "enable_rotated_room_geom") {
		if (!read_bool(fp, global_building_params.enable_rotated_room_geom)) {buildings_file_err(str, error);}
	}
	else if (unmatched) {
		cout << "Unrecognized buildings keyword in input file: " << str << endl;
		error = 1;
	}
	return !error;
}


void building_params_t::add_cur_mat() {
	unsigned const mat_ix(materials.size());

	for (unsigned n = 0; n < cur_prob; ++n) { // add more references to this mat for higher probability
		mat_gen_ix.push_back(mat_ix);
		(cur_mat.no_city ? mat_gen_ix_nocity : mat_gen_ix_city).push_back(mat_ix);
		if (cur_mat.house_prob > 0.0) {mat_gen_ix_res.push_back(mat_ix);}
	}
	materials.push_back(cur_mat);
	materials.back().finalize();
	materials.back().update_range(range_translate);
	has_normal_map |= cur_mat.has_normal_map();
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
	floor_spacing = 1.0/(2.0*ty);
	floorplan_wind_xscale = 2.0f*tx;
}

