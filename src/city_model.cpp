// 3D World - Loading and Drawing of 3D Models for Cities
// by Frank Gennari
// 6/5/2020
#include "city.h"
#include "file_utils.h"

extern city_params_t city_params;


bool city_model_t::read(FILE *fp, bool is_helicopter) {

	// filename recalc_normals two_sided centered body_material_id fixed_color_id xy_rot swap_xy scale lod_mult <blade_mat_id for helicopter> [shadow_mat_ids]
	assert(fp);
	unsigned line_num(0), swap_xyz(0), shadow_mat_id(0); // Note: line_num is unused
	fn = read_quoted_string(fp, line_num);
	if (fn.empty()) return 0;
	if (!read_int  (fp, recalc_normals)) return 0; // 0,1,2
	if (!read_bool (fp, two_sided))      return 0;
	if (!read_int  (fp, centered))       return 0; // bit mask for {X, Y, Z}
	if (!read_int  (fp, body_mat_id))    return 0;
	if (!read_int  (fp, fixed_color_id)) return 0;
	if (!read_float(fp, xy_rot))         return 0;
	if (!read_uint(fp, swap_xyz))        return 0; // {swap none, swap Y with Z, swap X with Z}
	if (!read_float(fp, scale))          return 0;
	if (!read_float(fp, lod_mult) || lod_mult < 0.0)  return 0;
	if (is_helicopter && !read_int(fp, blade_mat_id)) return 0;
	shadow_mat_ids.clear();
	while (read_uint(fp, shadow_mat_id)) {shadow_mat_ids.push_back(shadow_mat_id);}
	swap_xz = bool(swap_xyz & 2);
	swap_yz = bool(swap_xyz & 1);
	valid   = 1; // success
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


model3d &city_model_loader_t::get_model3d(unsigned id) {
	city_model_t const &model(get_model(id));
	assert(model.is_loaded());
	assert((size_t)model.model3d_id < size());
	return operator[](model.model3d_id);
}
vector3d city_model_loader_t::get_model_world_space_size(unsigned id) { // Note: may need to load model
	if (!is_model_valid(id)) return zero_vector; // error?
	city_model_t const &model_file(get_model(id));
	vector3d sz(get_model3d(id).get_bcube().get_size());
	if (model_file.swap_xz) {std::swap(sz.x, sz.z);}
	if (model_file.swap_yz) {std::swap(sz.y, sz.z);}
	if (round_fp(model_file.xy_rot/90.0) & 1) {std::swap(sz.x, sz.y);} // swap x/y for 90 and 270 degree rotations
	return sz;
}
colorRGBA city_model_loader_t::get_avg_color(unsigned id) {
	if (!is_model_valid(id)) return BLACK; // error?
	return get_model3d(id).get_avg_color();
}
bool city_model_loader_t::is_model_valid(unsigned id) {
	ensure_models_loaded(); // I guess we have to load the models here to determine if they're valid
	return get_model(id).is_loaded();
}

void city_model_loader_t::load_models() {
	for (unsigned i = 0; i < num_models(); ++i) {load_model_id(i);}
}
void city_model_loader_t::load_model_id(unsigned id) {
	unsigned const num_sub_models(get_num_sub_models(id));

	for (unsigned sm = 0; sm < num_sub_models; ++sm) { // load all sub-models
		city_model_t &model(get_model(id + (sm << 8)));
		if (model.is_loaded()) continue; // already loaded (should never get here?)
		bool skip_model(!have_buildings() && id < OBJ_MODEL_FHYDRANT); // building model, but no buildings, don't need to load
		skip_model |= (id == OBJ_MODEL_UMBRELLA && city_params.num_peds == 0); // don't need to load the umbrella model if there are no pedestrians
		if (skip_model || model.fn.empty()) continue;
		int const def_tid(-1); // should this be a model parameter?
		colorRGBA const def_color(WHITE); // should this be a model parameter?
		model.model3d_id = size(); // set before adding the model

		if (!load_model_file(model.fn, *this, geom_xform_t(), def_tid, def_color, 0, 0.0, model.recalc_normals, 0, city_params.convert_model_files, 1)) {
			cerr << "Error: Failed to read model file '" << model.fn << "'; Skipping this model";
			if (has_low_poly_model()) {cerr << " (will use default low poly model)";}
			cerr << "." << endl;
			model.model3d_id = -1; // invalid
			continue;
		}
		if (model.shadow_mat_ids.empty()) { // empty shadow_mat_ids, create the list from all materials
			model3d const &m(back());
			unsigned const num_materials(max(m.num_materials(), size_t(1))); // max with 1 for unbound material
			for (unsigned j = 0; j < num_materials; ++j) {model.shadow_mat_ids.push_back(j);} // add them all
		}
	} // for sm
}

void city_model_loader_t::draw_model(shader_t &s, vector3d const &pos, cube_t const &obj_bcube, vector3d const &dir, colorRGBA const &color,
	vector3d const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail, bool enable_animations, unsigned skip_mat_mask)
{
	bool const is_valid(is_model_valid(model_id));
	assert(is_valid); // must be loaded
	city_model_t const &model_file(get_model(model_id));
	model3d &model(get_model3d(model_id));
	if (!is_shadow_pass && model_file.body_mat_id >= 0 && color.A != 0.0) {model.set_color_for_material(model_file.body_mat_id, color);} // use custom color for body material
	model.bind_all_used_tids();
	cube_t const &bcube(model.get_bcube());
	point const orig_camera_pos(camera_pdu.pos), bcube_center(bcube.get_cube_center());
	camera_pdu.pos += bcube_center - pos - xlate; // required for distance based LOD
	bool const camera_pdu_valid(camera_pdu.valid);
	camera_pdu.valid = 0; // disable VFC, since we're doing custom transforms here
	// Note: in model space, front-back=z, left-right=x, top-bot=y (for model_file.swap_yz=1)
	float const sz_scale(obj_bcube.get_size().sum() / bcube.get_size().sum());
	float const height(model_file.swap_xz ? bcube.dx() : (model_file.swap_yz ? bcube.dy() : bcube.dz()));
	float const z_offset(0.5*height - (pos.z - obj_bcube.z1())/sz_scale); // translate required to map bottom of model to bottom of obj_bcube post transform
	
	if (enable_animations) {
		s.add_uniform_float("animation_scale",    model_file.scale/sz_scale); // Note: determined somewhat experimentally
		s.add_uniform_float("model_delta_height", (0.1*height + (model_file.swap_xz ? bcube.x1() : (model_file.swap_yz ? bcube.y1() : bcube.z1()))));
	}
	fgPushMatrix();
	translate_to(pos + vector3d(0.0, 0.0, z_offset*sz_scale)); // z_offset is in model space, scale to world space
	if (fabs(dir.y) > 0.001) {rotate_to_plus_x(dir);} // orient facing front
	else if (dir.x < 0.0) {fgRotate(180.0, 0.0, 0.0, 1.0);}
	if (dir.z != 0.0) {fgRotate(TO_DEG*asinf(-dir.z), 0.0, 1.0, 0.0);} // handle cars on a slope
	if (model_file.xy_rot != 0.0) {fgRotate(model_file.xy_rot, 0.0, 0.0, 1.0);} // apply model rotation about z/up axis (in degrees)
	if (model_file.swap_xz) {fgRotate(90.0, 0.0, 1.0, 0.0);} // swap X and Z dirs; models have up=X, but we want up=Z
	if (model_file.swap_yz) {fgRotate(90.0, 1.0, 0.0, 0.0);} // swap Y and Z dirs; models have up=Y, but we want up=Z
	uniform_scale(sz_scale); // scale from model space to the world space size of our target cube, using a uniform scale based on the averages of the x,y,z sizes
	point center(all_zeros);
	UNROLL_3X(if (!(model_file.centered & (1<<i_))) {center[i_] = bcube_center[i_];}); // use centered bit mask to control which component is centered vs. translated
	translate_to(-center); // cancel out model local translate
	bool const disable_cull_face_ths_obj(/*!is_shadow_pass &&*/ model_file.two_sided && glIsEnabled(GL_CULL_FACE));
	if (disable_cull_face_ths_obj) {glDisable(GL_CULL_FACE);}

	if (skip_mat_mask > 0) { // draw select materials
		for (unsigned i = 0; i < model.num_materials(); ++i) {
			if (skip_mat_mask & (1<<i)) continue; // skip this material
			model.render_material(s, i, is_shadow_pass, 0, 2, 0);
		}
	}
	else if (low_detail || is_shadow_pass) { // low detail pass, normal maps disabled
		if (!is_shadow_pass && use_model3d_bump_maps()) {model3d::bind_default_flat_normal_map();} // still need to set the default here in case the shader is using it
		// TODO: combine shadow materials into a single VBO and draw with one call when is_shadow_pass==1; this is complex and may not yield a significant improvement
		for (auto i = model_file.shadow_mat_ids.begin(); i != model_file.shadow_mat_ids.end(); ++i) {model.render_material(s, *i, is_shadow_pass, 0, 2, 0);}
	}
	else { // draw all materials
		float lod_mult(model_file.lod_mult); // should model_file.lod_mult always be multiplied by sz_scale?
		if (model_file.lod_mult == 0.0) {lod_mult = 400.0*sz_scale;} // auto select lod_mult
		model.render_materials(s, is_shadow_pass, 0, 0, 2, 3, 3, model.get_unbound_material(), rotation_t(),
			nullptr, nullptr, is_shadow_pass, lod_mult, (is_shadow_pass ? 10.0 : 0.0)); // enable_alpha_mask=2 (both)
	}
	if (disable_cull_face_ths_obj) {glEnable(GL_CULL_FACE);} // restore previous value
	fgPopMatrix();
	camera_pdu.valid = camera_pdu_valid;
	camera_pdu.pos   = orig_camera_pos;
	select_texture(WHITE_TEX); // reset back to default/untextured
}

unsigned get_model_id(unsigned id) { // first 8 bits = model_id, second 8 bits = sub_model_id
	unsigned const model_id(id & 0xFF);
	assert(model_id < NUM_OBJ_MODELS);
	return model_id;
}
unsigned get_sub_model_id(unsigned id) {return (id >> 8);}

unsigned car_model_loader_t       ::num_models() const {return city_params.car_model_files.size();}
unsigned helicopter_model_loader_t::num_models() const {return city_params.hc_model_files .size();}
unsigned object_model_loader_t    ::num_models() const {return NUM_OBJ_MODELS;}

city_model_t const &car_model_loader_t::get_model(unsigned id) const {
	assert(id < num_models());
	return city_params.car_model_files[id];
}
city_model_t &car_model_loader_t::get_model(unsigned id) {
	assert(id < num_models());
	return city_params.car_model_files[id];
}

city_model_t const &helicopter_model_loader_t::get_model(unsigned id) const {
	assert(id < num_models());
	return city_params.hc_model_files[id];
}
city_model_t &helicopter_model_loader_t::get_model(unsigned id) {
	assert(id < num_models());
	return city_params.hc_model_files[id];
}

unsigned object_model_loader_t::get_num_sub_models(unsigned id) const {
	return city_params.building_models[get_model_id(id)].size();
}
city_model_t const &object_model_loader_t::get_model(unsigned id) const {
	auto const &models(city_params.building_models[get_model_id(id)]);
	if (models.empty()) return null_model;
	return models[get_sub_model_id(id) % models.size()]; // sel_ix will wrap around if too large, which allows rand() to be passed in
}
city_model_t &object_model_loader_t::get_model(unsigned id) {
	auto &models(city_params.building_models[get_model_id(id)]);
	if (models.empty()) return null_model; // can get here for OBJ_MODEL_MONITOR
	return models[get_sub_model_id(id) % models.size()]; // sel_ix will wrap around if too large, which allows rand() to be passed in
}

bool city_params_t::add_model(unsigned id, FILE *fp) {
	assert(id < NUM_OBJ_MODELS);
	city_model_t model;
	if (!model.read(fp)) return 0;
	bool const filename_valid(model.check_filename());
	if (!filename_valid) {cerr << "Error: model file '" << model.fn << "' does not exist; skipping" << endl;} // nonfatal
	if (filename_valid || building_models[id].empty()) {building_models[id].push_back(model);} // add if valid or the first model
	return 1;
}

