// 3D World - Screenshot/Picture System for Buildings
// by Frank Gennari
// 10/31/20
#include "3DWorld.h"
#include "mesh.h"
#include "textures.h"
#include "buildings.h"
#include "openal_wrap.h" // for teleport sound


extern bool def_tex_compress;
extern int window_width, window_height, animate2, frame_counter;
extern float def_tex_aniso, NEAR_CLIP;
extern double tfticks;
extern vector<texture_t> textures;

void set_camera_pos_dir(point const &pos, vector3d const &dir);
void get_security_camera_info(room_object_t const &c, point &lens_pt, point &rot_pt, vector3d &camera_dir, vector3d &rot_axis, float &rot_angle);

vector3d get_global_camera_space_offset() {return vector3d((xoff2 - xoff)*DX_VAL, (yoff2 - yoff)*DY_VAL, 0.0);}


// *** Pictures and Screenshots ***

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
		std::string const name("screenshot_" + std::to_string(screenshots.size()));
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
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
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


// *** Security Cameras and Monitors ***

unsigned const SEC_CAMERA_XSIZE(640), SEC_CAMERA_YSIZE(480); // or use half the screen resolution? or the monitor aspect ratio?

class video_camera_manager_t {
	struct camera_t {
		point pos;
		vector3d dir;
		unsigned obj_ix, tid=0;
		bool valid=1;
		camera_t(point const &p, vector3d const &d, unsigned ix) : pos(p), dir(d), obj_ix(ix) {}
	};
	struct monitor_t {
		unsigned obj_ix;
		int camera_ix; // -1 is unset
		monitor_t(unsigned ix, int cix=-1) : obj_ix(ix), camera_ix(cix) {}
	};
	bool update_this_frame=0;
	unsigned update_ix=0;
	building_t *cur_building=nullptr;
	rand_gen_t rgen; // for noise texture offset
	vector<camera_t > cameras;
	vector<monitor_t> monitors;
	vector<unsigned> tid_free_list, cur_monitors;

	unsigned alloc_tid() {
		if (!tid_free_list.empty()) {
			unsigned const tid(tid_free_list.back());
			tid_free_list.pop_back();
			return tid;
		}
		unsigned tid(0);
		setup_texture(tid, 0, 0, 0); // no mipmap or wrap
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, SEC_CAMERA_XSIZE, SEC_CAMERA_YSIZE, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
		assert(tid > 0);
		return tid;
	}
public:
	void register_building(building_t &b, room_t const &room) {
		update_this_frame = 1; // enable update_cameras() after this call
		if (&b == cur_building) return; // same building - nothing to do
		cur_building = &b;
		update_ix    = 0; // reset the camera cycle

		for (camera_t const &camera : cameras) { // store allocated tids on the free list before clearing cameras
			if (camera.tid > 0) {tid_free_list.push_back(camera.tid);}
		}
		cameras .clear();
		monitors.clear();
		auto objs_end(b.interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
		vect_room_object_t const &objs(b.interior->room_geom->objs);

		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->type == TYPE_CAMERA) { // register camera
				point lens_pt, rot_pt;
				vector3d camera_dir, rot_axis;
				float rot_angle(0.0);
				get_security_camera_info(*i, lens_pt, rot_pt, camera_dir, rot_axis, rot_angle);
				rotate_vector3d(rot_axis, rot_angle, camera_dir);
				lens_pt -= rot_pt;
				rotate_vector3d(rot_axis, rot_angle, lens_pt); // lens rotates about rot_pt
				lens_pt += rot_pt;
				//lens_pt += 0.1*NEAR_CLIP*camera_dir; // move slightly in front of the camera; is this necessary?
				cameras.emplace_back(lens_pt, camera_dir, (i - objs.begin()));
			}
			else if (i->type == TYPE_MONITOR && room.contains_cube(*i)) { // register monitor
				monitors.emplace_back(i - objs.begin());
			}
		} // for i
		assert(!cameras.empty()); // too strong?
		if (cameras.empty()) return;
		// round robin assign cameras to monitors; if there are more monitors than cameras, some cameras will be duplicated;
		// if there are more cameras than monitors, so cameras (on upper floors) won't be shown - or we could only show one per floor rather than each end of the hallway?
		for (unsigned i = 0; i < monitors.size(); ++i) {monitors[i].camera_ix = (i % cameras.size());}
	}
	void update_cameras() {
		if (!update_this_frame) return; // player not in security room; this also keeps us from using an invalid/deleted building
		update_this_frame = 0; // update once, until register_building() is called again
		if (!animate2) return; // updates disabled
		if (cur_building == nullptr) return; // error?
		if (monitors.empty() || cameras.empty()) return; // nothing to do
		vector3d const xlate(get_tiled_terrain_model_xlate()); // could pass this in

		for (unsigned n = 0; n < cameras.size(); ++n) { // try to find a valid camera and valid + visible monitor to update
			unsigned const camera_ix(update_ix);
			assert(camera_ix < cameras.size());
			update_ix = (update_ix + 1) % cameras.size();
			camera_t &camera(cameras[camera_ix]);
			vect_room_object_t &objs(cur_building->interior->room_geom->objs);
			assert(camera.obj_ix < objs.size());
			room_object_t const &camera_obj(objs[camera.obj_ix]);
			cur_monitors.clear();

			for (monitor_t const &m : monitors) {
				if (m.camera_ix != camera_ix) continue;
				assert(m.obj_ix < objs.size());
				room_object_t const &monitor(objs[m.obj_ix]);
				if (monitor.type != TYPE_MONITOR) continue; // taken by the player? should we remove this monitor from the list?
				if (monitor.is_broken() || (monitor.obj_id & 1)) continue; // broken or turned off
				if (!cur_building->is_rot_cube_visible(monitor, xlate)) continue; // monitor not visible to the player
				// Note: no occlusion culling since player is likely already in this room
				cur_monitors.push_back(m.obj_ix);
			} // for m
			if (cur_monitors.empty()) continue; // no monitors observing this camera
			if (camera_obj.type != TYPE_CAMERA) {camera.valid = 0;} // taken by the player? should we remove this camera from the list?
			else if (!camera_obj.is_powered())  {camera.valid = 0;} // no power
			else {camera.valid = 1;}
			//cout << TXT(n) << TXT(monitors.size()) << TXT(cameras.size()) << TXT(camera.valid) << TXT(cur_monitors.size()) << endl; // TESTING

			if (camera.valid) { // update camera texture
				if (camera.tid == 0) {camera.tid = alloc_tid();} // allocate the texture the first time this camera is rendered to
				// TODO: render cur_building to texture with camera.tid to update this camera
			}
			for (unsigned ix : cur_monitors) {
				room_object_t &monitor(objs[ix]);
				monitor.item_flags = camera_ix;
				monitor.flags     |= RO_FLAG_IS_ACTIVE; // flag as having an active texture
			}
			break; // done
		} // for n
	}
	void setup_monitor_screen_draw(room_object_t const &monitor, rgeom_mat_t &mat) {
		assert(monitor.type == TYPE_MONITOR);
		if (!monitor.is_active()) {select_texture(BLACK_TEX); return;} // error?
		assert(monitor.item_flags < cameras.size());
		camera_t const &camera(cameras[monitor.item_flags]);

		if (!camera.valid) { // draw noise
			mat.tex.tscale_x = 6.0/monitor.get_width ();
			mat.tex.tscale_y = 6.0/monitor.get_height();
			mat.tex.txoff    = rgen.rand_float(); // animate it
			mat.tex.tyoff    = rgen.rand_float();
			select_texture(PS_NOISE_TEX);
			return;
		}
		mat.tex.tscale_x = mat.tex.tscale_y = 0.0; // map texture to quad
		//select_texture((frame_counter + monitor.item_flags)%20); // TESTING
		assert(camera.tid > 0);
		bind_2d_texture(camera.tid);
	}
	unsigned get_gpu_mem() const {
		unsigned num_tids(tid_free_list.size());
		for (camera_t const &camera : cameras) {num_tids += (camera.tid > 0);}
		return 3*SEC_CAMERA_XSIZE*SEC_CAMERA_YSIZE*num_tids; // or 4x?
	}
};

video_camera_manager_t video_camera_manager; // reused across buildings

void building_t::update_security_cameras(point const &camera_bs) {
	return; // TODO: enable when this is working
	if (!has_room_geom()) return;
	if (interior->security_room_ix < 0) return; // no security room set, or not yet populated
	room_t const &sec_room(get_room(interior->security_room_ix)); // Note: security room is on the ground floor
	assert(sec_room.get_room_type(0) == RTYPE_SECURITY);
	if (!sec_room.contains_pt(camera_bs)) return;
	if (camera_bs.z > ground_floor_z1 + get_window_vspace()) return; // player not the ground floor (which contains the security room)
	video_camera_manager.register_building(*this, sec_room);
}
void update_security_camera_image() {video_camera_manager.update_cameras();}
void setup_monitor_screen_draw(room_object_t const &monitor, rgeom_mat_t &mat) {video_camera_manager.setup_monitor_screen_draw(monitor, mat);}
unsigned get_building_textures_gpu_mem() {return video_camera_manager.get_gpu_mem();}


