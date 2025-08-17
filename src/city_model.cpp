// 3D World - Loading and Drawing of 3D Models for Cities
// by Frank Gennari
// 6/5/2020
#include "city.h"
#include "file_utils.h"
#include "format_text.h"

extern city_params_t city_params;

bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, string const &anim_name, int recalc_normals, bool verbose);


bool read_keyword(FILE *fp, string &str) {
	str.clear();
	
	while (1) {
		int c(getc(fp));
		if (is_EOF(c) || c == '\n') break;

		if (isspace(c)) {
			if (str.empty()) continue; // leading whitespace
			break; // end of string
		}
		if (c == '=') {
			if (str.empty()) {
				cerr << "Error: Found stray '=' with no keyword in config file" << endl;
				return 0;
			}
			return 1; // done/success
		}
		if (c == '#' || (!isalpha(c) && c != '_')) {ungetc(c, fp); break;} // end of string
		str.push_back(c);
	} // end while
	if (str.empty()) return 1;
	cerr << "Error: Found keyword string in config file with missing '=': " << str << endl;
	return 0;
}

bool city_model_t::read(FILE *fp, bool is_helicopter, bool is_person) {

	// filename recalc_normals two_sided centered body_material_id fixed_color_id xy_rot swap_xy scale lod_mult <blade_mat_id for helicopter> [keywords] [shadow_mat_ids]
	assert(fp);
	unsigned swap_xyz(0), shadow_mat_id(0);
	fn = read_quoted_string(fp);
	if (fn.empty()) return 0;
	if (!read_int  (fp, recalc_normals)) return 0; // 0,1,2
	if (!read_bool (fp, two_sided))      return 0;
	if (!read_int  (fp, centered))       return 0; // bit mask for {X, Y, Z}
	if (!read_int  (fp, body_mat_id))    return 0; // -1=all
	if (!read_int  (fp, fixed_color_id)) return 0; // -1=none, -2=specified RGBA, -3=auto-from-model
	if (fixed_color_id == -2 && !read_color(fp, custom_color)) return 0; // read RGBA color
	if (!read_float(fp, xy_rot))         return 0;
	if (!read_uint(fp, swap_xyz))        return 0; // {swap none, swap Y with Z, swap X with Z}
	if (!read_float(fp, scale))          return 0;
	if (!read_float(fp, lod_mult) || lod_mult < 0.0)  return 0;
	
	if (is_helicopter) { // helicopter model is special because it has a blade_mat_id
		if (!read_int(fp, blade_mat_id)) return 0;
	}
	if (is_person) { // read animation data, etc.
		if (!read_float(fp, anim_speed)) return 0;
		if (!read_bool (fp, is_zombie))  return 0;
	}
	// read any keywords; must be before reading shadow_mat_ids, before the newline is encountered
	string keyword;

	while (1) {
		if (!read_keyword(fp, keyword)) return 0;
		if (keyword.empty()) break;

		if (keyword == "reverse_winding") {
			if (!read_uint64(fp, rev_winding_mask)) return 0;
		}
		else if (keyword == "mirrored") {
			if (!read_bool(fp, mirrored)) return 0;
		}
		else {
			cerr << "Error: Unrecognized keyword " << keyword << " in config file" << endl;
			return 0;
		}
	} // end while
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
bool city_model_t::has_animation(string const &anim_name) const {
	if (anim_name == default_anim_name) return 1;

	for (model_anim_t const &a : anim_fns) {
		if (a.anim_name == anim_name) return 1;
	}
	return 0;
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
	vector3d sz(get_model3d(id).get_bcube().get_size()*model_file.scale);
	if (model_file.swap_xz) {std::swap(sz.x, sz.z);}
	if (model_file.swap_yz) {std::swap(sz.y, sz.z);}
	if (round_fp(model_file.xy_rot/90.0) & 1) {std::swap(sz.x, sz.y);} // swap x/y for 90 and 270 degree rotations
	return sz;
}
colorRGBA city_model_loader_t::get_avg_color(unsigned id, bool area_weighted) {
	if (!is_model_valid(id)) return BLACK; // error?
	return get_model3d(id).get_and_cache_avg_color(area_weighted);
}
bool city_model_loader_t::model_filename_contains(unsigned id, string const &str, string const &str2) const {
	string const &fn(get_model(id).fn);
	if (fn.find(str) != string::npos) return 1;
	return (!str2.empty() && fn.find(str2) != string::npos);
}
bool city_model_loader_t::is_model_valid(unsigned id) {
	city_model_t &model(get_model(id));
	if (!model.tried_to_load) {load_model_id(id);} // load the model if needed
	return model.is_loaded();
}

void city_model_loader_t::load_model_id(unsigned id) { // currently up to 72 models loaded with people and cars
	unsigned const num_sub_models(get_num_sub_models(id));

	for (unsigned sm = 0; sm < num_sub_models; ++sm) { // load all sub-models
		city_model_t &model(get_model(combine_model_submodel_id(id, sm)));
		if (model.tried_to_load) continue; // already tried to load (should never get here?)
		if (can_skip_model(id) || model.fn.empty()) continue;
		int const def_tid(-1); // should this be a model parameter?
		colorRGBA const def_color(WHITE); // should this be a model parameter?
		model.tried_to_load   = 1; // flag, even if load fails
		model.model3d_id      = size(); // set before adding the model
		bool const verbose    = 1;
		float const lod_scale = 1.0; // or use model.lod_mult directly and not use it during drawing?

		if (!load_model_file(model.fn, *this, geom_xform_t(), model.default_anim_name, def_tid, def_color, 0, 0.0,
			lod_scale, model.recalc_normals, 0, city_params.convert_model_files, verbose, model.rev_winding_mask))
		{
			cerr << format_red("Error: Failed to read model file '" + model.fn + "'; Skipping this model");
			if (has_low_poly_model()) {cerr << " (will use default low poly model)";}
			cerr << "." << endl;
			model.model3d_id = -1; // invalid
			continue;
		}
		assert(!empty()); // must have loaded a model
		model3d &cur_model(back());

		if (model.shadow_mat_ids.empty()) { // empty shadow_mat_ids, create the list from all materials
			unsigned const num_materials(max(cur_model.num_materials(), size_t(1))); // max with 1 for unbound material
			for (unsigned j = 0; j < num_materials; ++j) {model.shadow_mat_ids.push_back(j);} // add them all
		}
		for (city_model_t::model_anim_t const &anim : model.anim_fns) {
			model3d anim_data(anim.fn, tmgr); // Note: texture manager is passed in, even though there should be no loaded textures; however, this isn't checked
			
			if (!read_assimp_model(anim.fn, anim_data, geom_xform_t(), anim.anim_name, model.recalc_normals, verbose)) {
				cerr << format_red("Error: Failed to read model animation file '" + anim.fn + "'; Skipping this animation") << endl;
			}
			else {cur_model.merge_animation_from(anim_data);}
		} // for anim
		city_params.any_model_has_animations |= cur_model.has_animations();
	} // for sm
}

bool object_model_loader_t::can_skip_model(unsigned id) const {
	if (!have_buildings() && id < OBJ_MODEL_FHYDRANT) return 1; // building model, but no buildings, don't need to load
	if (id == OBJ_MODEL_UMBRELLA && city_params.num_peds == 0) return 1; // don't need to load the umbrella model if there are no pedestrians
	return 0;
}

/*static*/ void city_model_loader_t::rotate_model_from_plus_x_to_dir(vector3d const &dir) {
	if (fabs(dir.y) > 0.001) {rotate_to_plus_x(dir);} // orient facing front
	else if (dir.x < 0.0) {fgRotate(180.0, 0.0, 0.0, 1.0);}
}

void enable_animations_for_shader(shader_t &s) {
	if (city_params.use_animated_people && city_params.any_model_has_animations) {s.set_prefix("#define USE_BONE_ANIMATIONS", 0);} // VS
	s.add_property("animation_shader", "pedestrian_animation.part+"); // this shader part now contains model bone animations as well
}
void set_anim_id(shader_t &s, bool enable_animations, int animation_id, unsigned model_anim_id, unsigned model_anim_id2, bool has_bone_animations) {
	if (!enable_animations) return;

	if (city_params.use_animated_people && has_bone_animations && animation_id == ANIM_ID_WALK) { // select bone animation rather than walking
		animation_id = ANIM_ID_SKELETAL; // used in the shader to select skeletal animation
		assert(model_anim_id < NUM_MODEL_ANIMS);
		s.add_property("animation_name", animation_names[model_anim_id]);

		if (model_anim_id2 != model_anim_id) { // blended animation
			assert(model_anim_id2 < NUM_MODEL_ANIMS);
			s.add_property("animation_name2", animation_names[model_anim_id2]);
		}
	}
	s.add_uniform_int("animation_id", animation_id);
}

// walking animations used by people use animation blending and shader animation_name properties; other animations use the stored model_anim_id
int animation_state_t::get_anim_id_for_setup_bone_transforms () const {return ((anim_id == ANIM_ID_WALK) ? -1 : model_anim_id );}
int animation_state_t::get_anim_id2_for_setup_bone_transforms() const {return ((anim_id == ANIM_ID_WALK) ? -1 : model_anim_id2);}

void animation_state_t::set_animation_id_and_time(shader_t &s, bool has_bone_animations, float anim_speed) const {
	set_anim_id(s, enabled, anim_id, model_anim_id, model_anim_id2, has_bone_animations);
	if (enabled && !(city_params.use_animated_people && has_bone_animations)) {s.add_uniform_float("animation_time", anim_speed*anim_time);} // only for custom animations
}
void animation_state_t::clear_animation_id(shader_t &s) const {set_anim_id(s, enabled, ANIM_ID_NONE, 0, 0, 0);} // has_bone_animations not needed here

bool city_model_loader_t::check_anim_wrapped(unsigned model_id, unsigned model_anim_id, float old_time, float new_time) {
	bool const is_valid(is_model_valid(model_id));
	assert(is_valid); // caller should have checked this previously
	return get_model3d(model_id).check_anim_wrapped(model_anim_id, old_time, new_time);
}
float city_model_loader_t::get_anim_duration(unsigned model_id, unsigned model_anim_id) { // in seconds
	bool const is_valid(is_model_valid(model_id));
	assert(is_valid); // caller should have checked this previously
	return get_model3d(model_id).get_anim_duration(model_anim_id);
}

void city_model_loader_t::draw_model(shader_t &s, vector3d const &pos, cube_t const &obj_bcube, vector3d const &dir, colorRGBA const &color,
	vector3d const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail, animation_state_t *anim_state, unsigned skip_mat_mask, bool untextured,
	bool force_high_detail, bool upside_down, bool emissive, bool do_local_rotate, int mirror_dim, bool using_custom_tid, int swap_xy_mode)
{
	assert(!(low_detail && force_high_detail));
	bool const is_valid(is_model_valid(model_id)); // first 8 bits is model ID, last 8 bits is sub-model ID
	if (!is_valid) {cerr << "Invalid model ID: " << model_id << endl;}
	assert(is_valid); // must be loaded
	city_model_t const &model_file(get_model(model_id));
	model3d &model(get_model3d(model_id));
	bool const is_animated(model.has_animations());
	bool const have_body_mat(model_file.body_mat_id >= 0);
	bool const use_custom_color   (!is_shadow_pass && have_body_mat && color.A != 0.0);
	bool const use_custom_emissive(!is_shadow_pass && have_body_mat && emissive);
	bool const use_custom_texture (!is_shadow_pass && have_body_mat && (untextured || using_custom_tid));
	colorRGBA orig_color;
	int orig_tid(-1);
	if (use_custom_color   ) {orig_color = model.set_color_for_material(model_file.body_mat_id, color);} // use custom color for body material
	if (use_custom_emissive) {model.set_material_emissive_color(model_file.body_mat_id, color);} // use custom color for body material

	if (use_custom_texture) {
		// Note: this indexes into the model's texture_manager textures vector; it's not a global texture ID
		int const untex_tid(using_custom_tid ? -2 : -1); // -2: texture is pre-bound, don't set; -1: use white texture
		if (using_custom_tid) {assert(model_file.body_mat_id == 0);} // caller must bind the texture to TU0; only works when body_mat_id is the first material
		orig_tid = model.set_texture_for_material(model_file.body_mat_id, untex_tid);
	}
	model.bind_all_used_tids();
	cube_t const bcube(model.get_bcube() * model_file.model_anim_scale);
	point const orig_camera_pos(camera_pdu.pos), bcube_center(bcube.get_cube_center());
	camera_pdu.pos += bcube_center - pos - xlate; // required for distance based LOD
	bool const camera_pdu_valid(camera_pdu.valid);
	camera_pdu.valid = 0; // disable VFC, since we're doing custom transforms here
	// Note: in model space, front-back=z, left-right=x, top-bot=y (for model_file.swap_yz=1)
	bool const swap_draw_xz(swap_xy_mode == 1), swap_draw_yz(swap_xy_mode == 2);
	float const model_height(model_file.swap_xz ? bcube    .dx() : (model_file.swap_yz ? bcube    .dy() : bcube    .dz()));
	float const obj_height  (swap_draw_xz       ? obj_bcube.dx() : (swap_draw_yz       ? obj_bcube.dy() : obj_bcube.dz()));
	// animated models don't have valid bcubes because they sometimes start in a bind pose, so use the height as the size scale since it's more likely to be accurate
	assert(obj_bcube.is_strictly_normalized());
	float sz_scale(0.0);
	if (is_animated) {sz_scale = (obj_height / model_height);} // use zsize only for scale
	else             {sz_scale = (obj_bcube.get_size().sum() / bcube.get_size().sum());} // use average XYZ size for scale
	// translate required to map bottom of model to bottom of obj_bcube post transform; up to the caller to handle this for swap_xy_mode > 0
	float const z_offset(swap_xy_mode ? 0.0 : (0.5*model_height - (pos.z - obj_bcube.z1())/sz_scale));
	
	if (anim_state && anim_state->enabled) {
		// skip expensive animations if low detail; may cause the model to T-pose, but it should be far enough that the user can't tell
		bool const use_bone_animations(city_params.use_animated_people && is_animated && !low_detail &&
			(anim_state->anim_id == ANIM_ID_WALK || anim_state->anim_id == ANIM_ID_SKELETAL)); // people use ANIM_ID_WALK; animals/others use ANIM_ID_SKELETAL
		float const anim_speed(anim_state->fixed_anim_speed ? 1.0 : model_file.anim_speed);
		anim_state->set_animation_id_and_time(s, use_bone_animations, anim_speed);

		if (use_bone_animations) {
			int const anim_id(anim_state->get_anim_id_for_setup_bone_transforms());
			float const speed_mult(SKELETAL_ANIM_TIME_CONST*anim_speed), anim_time(speed_mult*anim_state->anim_time);

			if (anim_state->blend_factor > 0.0) { // enable animation blending
				model.setup_bone_transforms_blended(s, anim_time, speed_mult*anim_state->anim_time2,
					anim_state->blend_factor, anim_id, anim_state->get_anim_id2_for_setup_bone_transforms());
			}
			else if (anim_state->cached) { // single cached animation
				model.setup_bone_transforms_cached(*anim_state->cached, s, anim_time, anim_id);
			}
			else { // single animation
				model.setup_bone_transforms(s, anim_time, anim_id);
			}
		}
		else {
			s.add_uniform_float("animation_scale",    model_file.scale/sz_scale); // Note: determined somewhat experimentally
			s.add_uniform_float("model_delta_height", (0.1*model_height + (model_file.swap_xz ? bcube.x1() : (model_file.swap_yz ? bcube.y1() : bcube.z1()))));
		}
	}
	bool const do_mirror(mirror_dim < 3);
	vector3d const local_rotate(do_local_rotate ? model_file.rotate_about : zero_vector);
	fgPushMatrix();
	translate_to(pos + z_offset*sz_scale*plus_z - local_rotate); // z_offset is in model space, scale to world space
	rotate_model_from_plus_x_to_dir(dir); // typically rotated about the Z axis
	if (local_rotate != all_zeros) {translate_to(local_rotate);}
	if (dir.z != 0.0      ) {fgRotate(TO_DEG*asinf(-dir.z), 0.0, 1.0, 0.0);} // handle cars on a slope
	if (do_mirror         ) {fgScale(((mirror_dim == 0) ? -1.0 : 1.0), ((mirror_dim == 1) ? -1.0 : 1.0), ((mirror_dim == 2) ? -1.0 : 1.0));}
	if (model_file.xy_rot != 0.0) {fgRotate(model_file.xy_rot, 0.0, 0.0, 1.0);} // apply model rotation about z/up axis (in degrees)
	if (swap_xy_mode      ) {fgRotate(-90.0, 1.0, 0.0, 0.0);} // swap Y and Z dirs - in both cases?
	if (model_file.swap_xz) {fgRotate( 90.0, 0.0, 1.0, 0.0);} // swap X and Z dirs; models have up=X, but we want up=Z
	if (model_file.swap_yz) {fgRotate( 90.0, 1.0, 0.0, 0.0);} // swap Y and Z dirs; models have up=Y, but we want up=Z
	if (upside_down       ) {fgRotate(180.0, 1.0, 0.0, 0.0);} // R180 about X to flip over
	uniform_scale(sz_scale); // scale from model space to the world space size of our target cube, using a uniform scale based on the averages of the x,y,z sizes
	point center(bcube_center);
	UNROLL_3X(if (model_file.centered & (1<<i_)) {center[i_] = 0.0;}); // use centered bit mask to control which component is centered vs. translated
	translate_to(-center); // cancel out model local translate
	bool const disable_cull_face_this_obj(/*!is_shadow_pass &&*/ model_file.two_sided && glIsEnabled(GL_CULL_FACE));
	if (disable_cull_face_this_obj) {glDisable(GL_CULL_FACE);}
	if (do_mirror) {glFrontFace(GL_CW);}

	if (!force_high_detail && (low_detail || is_shadow_pass)) { // low detail pass, normal maps disabled
		if (!is_shadow_pass && use_model3d_bump_maps()) {bind_default_flat_normal_map();} // still need to set the default here in case the shader is using it
		// combine shadow materials into a single VBO and draw with one call when is_shadow_pass==1? this is complex and may not yield a significant improvement
		for (auto i = model_file.shadow_mat_ids.begin(); i != model_file.shadow_mat_ids.end(); ++i) {
			if (skip_mat_mask & (1<<*i)) continue; // skip this material
			model.render_material(s, *i, is_shadow_pass, 0, 2, 0, nullptr, 1);
		}
	}
	else { // draw all materials
		// should model_file.lod_mult always be multiplied by sz_scale? this would make the config file values more consistent, but requires a lot of manual updating
		float lod_mult(model_file.lod_mult);
		if (!is_shadow_pass && model_file.lod_mult == 0.0) {lod_mult = 400.0*sz_scale;} // auto select lod_mult
		float const fixed_lod_dist((is_shadow_pass && !force_high_detail) ? 10.0 : 0.0);
		model.render_materials(s, is_shadow_pass, 0, 0, 2, 3, 3, model.get_unbound_material(), rotation_t(),
			nullptr, nullptr, is_shadow_pass, lod_mult, fixed_lod_dist, 0, 1, 1, skip_mat_mask); // enable_alpha_mask=2 (both), is_scaled=1, no_set_min_alpha=1
	}
	if (do_mirror) {glFrontFace(GL_CCW);}
	if (disable_cull_face_this_obj) {glEnable(GL_CULL_FACE);} // restore previous value
	fgPopMatrix();
	camera_pdu.valid = camera_pdu_valid;
	camera_pdu.pos   = orig_camera_pos;
	select_texture(WHITE_TEX); // reset back to default/untextured
	if (use_custom_color   ) {model.set_color_for_material     (model_file.body_mat_id, orig_color);} // restore original color
	if (use_custom_texture ) {model.set_texture_for_material   (model_file.body_mat_id, orig_tid  );} // restore original texture
	if (use_custom_emissive) {model.set_material_emissive_color(model_file.body_mat_id, BLACK     );} // reset
}

unsigned get_model_id(unsigned id) { // first 8 bits = model_id, second 8 bits = sub_model_id
	unsigned const model_id(id & 0xFF);
	assert(model_id < NUM_OBJ_MODELS);
	return model_id;
}

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
int object_model_loader_t::get_valid_sub_model_id(unsigned id, vector<city_model_t> const &models) const {
	unsigned const sub_model_id(get_sub_model_id(id)); // shift sub-model ID bits back to LSB; needed when models.size() is a power of 2

	for (unsigned i = 0; i < models.size(); ++i) { // check all models starting with the selected one and return the first valid
		unsigned const cand((sub_model_id + i) % models.size()); // index will wrap around if too large, which allows rand() to be passed in
		city_model_t const &model(models[cand]);
		if (model.valid && (!model.tried_to_load || model.is_loaded())) return cand;
	}
	return -1; // no valid models found
}
city_model_t const &object_model_loader_t::get_model(unsigned id) const {
	auto const &models(city_params.building_models[get_model_id(id)]);
	int const model_id(get_valid_sub_model_id(id, models));
	if (model_id < 0) return null_model;
	return models[model_id];
}
city_model_t &object_model_loader_t::get_model(unsigned id) {
	auto &models(city_params.building_models[get_model_id(id)]);
	int const model_id(get_valid_sub_model_id(id, models));
	if (model_id < 0) return null_model;
	return models[model_id];
}

bool city_params_t::add_model(unsigned id, FILE *fp) {
	if (id >= NUM_OBJ_MODELS) {cout << TXT(id) << TXT(NUM_OBJ_MODELS) << endl;}
	assert(id < NUM_OBJ_MODELS);
	city_model_t model;
	if (!model.read(fp)) return 0;
	model.default_anim_name = default_anim_name; // needed for birds
	model.model_anim_scale  = model_anim_scale ; // needed for birds
	bool const filename_valid(model.check_filename());
	if (!filename_valid) {cerr << format_red("Error: model file '" + model.fn + "' does not exist; skipping") << endl;} // nonfatal
	if (filename_valid || building_models[id].empty()) {building_models[id].push_back(model);} // add if valid or the first model
	return 1;
}

