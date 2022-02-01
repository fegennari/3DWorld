// 3D World - Screenshot/Picture System for Buildings
// by Frank Gennari
// 10/31/20
#include "3DWorld.h"
#include "mesh.h"
#include "textures.h"
#include "buildings.h"
#include "openal_wrap.h" // for teleport sound


extern bool def_tex_compress;
extern int window_width, window_height;
extern float def_tex_aniso;
extern vector<texture_t> textures;

void set_camera_pos_dir(point const &pos, vector3d const &dir);

vector3d get_global_camera_space_offset() {return vector3d((xoff2 - xoff)*DX_VAL, (yoff2 - yoff)*DY_VAL, 0.0);}

// Note: currently used for pictures hanging on building room walls, but could be made more generally useful
class screenshot_manager_t {
	struct screenshot_info_t {
		unsigned tid; // index into textures array
		point cpos; // camera position at the time the screenshot was taken, in global space
		vector3d cdir; // camera view direction at the time the screenshot was taken
		screenshot_info_t(unsigned tid_) : tid(tid_), cpos(get_camera_pos() + get_global_camera_space_offset()), cdir(cview_dir) {}
	};
	vector<screenshot_info_t> screenshots;
public:
	unsigned size() const {return screenshots.size();}
	bool empty   () const {return screenshots.empty();}
	void clear() {screenshots.clear();}
	int get_rand_tid(unsigned rand_ix) const {return (empty() ? -1 : screenshots[rand_ix%size()].tid);}

	void restore_camera_for_rand_ix(unsigned rand_ix) const {
		if (empty()) return;
		screenshot_info_t const &si(screenshots[rand_ix%size()]);
		set_camera_pos_dir((si.cpos - get_global_camera_space_offset()), si.cdir); // convert back to local camera space
	}
	void add_screenshot() {
		print_text_onscreen("Screenshot Saved as Texture", WHITE, 1.0, 2*TICKS_PER_SECOND, 0);
		unsigned tid(0);
		frame_buffer_to_texture(tid, 0);
		unsigned const tex_ix(textures.size());
		std::ostringstream oss;
		oss << "screenshot_" << screenshots.size();
		std::string const name(oss.str());
		// type format width height wrap_mir ncolors use_mipmaps name [invert_y=0 [do_compress=1 [anisotropy=1.0 [mipmap_alpha_weight=1.0 [normal_map=0]]]]]
		texture_t new_tex(0, 9, window_width, window_height, 0, 3, 0, name, 0, def_tex_compress, def_tex_aniso, 1.0, 0);
		new_tex.set_existing_tid(tid, WHITE); // not sure what to set the color to
		textures.push_back(new_tex);
		screenshots.emplace_back(tex_ix);
	}
};

screenshot_manager_t screenshot_manager;

void take_screenshot_texture() {screenshot_manager.add_screenshot();}
int get_rand_screenshot_texture(unsigned rand_ix) {return screenshot_manager.get_rand_tid(rand_ix);}
unsigned get_num_screenshot_tids() {return screenshot_manager.size();}

bool building_t::maybe_teleport_to_screenshot() const {
	if (!has_room_geom()) return 0;
	if (!is_house) return 0; // currently only houses have pictures hanging on walls
	if (get_num_screenshot_tids() == 0) return 0;
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip trim/buttons/stairs/elevators
	vector3d const xlate(get_tiled_terrain_model_xlate());
	point const camera_bs(get_camera_building_space());
	int closest_obj_id(-1);
	float dmin_sq(0.0);

	// find nearest screenshot picture
	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_PICTURE) continue;
		if (!camera_pdu.cube_visible(*i + xlate)) continue; // skip if invisible
		float const dist_sq(p2p_dist_sq(i->get_cube_center(), camera_bs));
		if (dmin_sq == 0.0 || dist_sq < dmin_sq) {closest_obj_id = i->obj_id; dmin_sq = dist_sq;}
	} // for i
	if (closest_obj_id < 0) return 0;
	screenshot_manager.restore_camera_for_rand_ix(closest_obj_id);
	gen_sound(SOUND_POWERUP, get_camera_pos(), 1.0, 0.6); // play teleport sound
	return 1;
}


