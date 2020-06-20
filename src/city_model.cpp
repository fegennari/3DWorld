// 3D World - Loading and Drawing of 3D Models for Cities
// by Frank Gennari
// 6/5/2020
#include "city.h"
#include "file_utils.h"

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
	valid = 1; // success
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
		int const recalc_normals = 1; // 0=no, 1=yes, 2=face_weight_avg

		if (!load_model_file(fn, *this, geom_xform_t(), -1, WHITE, 0, 0.0, recalc_normals, 0, city_params.convert_model_files, 1)) {
			cerr << "Error: Failed to read model file '" << fn << "'; Skipping this model (will use default low poly model)." << endl;
			push_back(model3d(fn, tmgr)); // add a placeholder dummy model
			models_valid[i] = 0;
		}
	} // for i
}

void city_model_loader_t::draw_model(shader_t &s, vector3d const &pos, cube_t const &obj_bcube, vector3d const &dir, colorRGBA const &color,
	vector3d const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail, bool enable_animations)
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

unsigned car_model_loader_t::num_models() const {return city_params.car_model_files.size();}

city_model_t const &car_model_loader_t::get_model(unsigned id) const {
	assert(id < num_models());
	return city_params.car_model_files[id];
}

city_model_t const &object_model_loader_t::get_model(unsigned id) const {
	switch (id) {
	case OBJ_MODEL_TOILET: return city_params.toilet_model;
	default: assert(0);
	}
	return null_model; // never gets here
}

