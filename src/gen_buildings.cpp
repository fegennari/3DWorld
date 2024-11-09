// 3D World - Building Generation
// by Frank Gennari
// 5/22/17

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "buildings.h"
#include "mesh.h"
#include "draw_utils.h" // for point_sprite_drawer_sized
#include "subdiv.h" // for sd_sphere_d
#include "tree_3dw.h" // for tree_placer_t
#include "profiler.h"
#include "lightmap.h" // for light_source

using std::string;

bool const ADD_ROOM_SHADOWS        = 1; // for room lights
bool const DRAW_EXT_REFLECTIONS    = 1; // draw building exteriors in mirror reflections; slower, but looks better; not shadowed
bool const DRAW_WALKWAY_INTERIORS  = 1;
float const WIND_LIGHT_ON_RAND     = 0.08;
unsigned const NO_SHADOW_WHITE_TEX = BLACK_TEX; // alias to differentiate shadowed    vs. unshadowed untextured objects
unsigned const SHADOW_ONLY_TEX     = RED_TEX;   // alias to differentiate shadow only vs. other      untextured objects

bool camera_in_building(0), interior_shadow_maps(0), player_is_hiding(0), player_in_unlit_room(0), player_in_walkway(0), player_in_int_elevator(0), player_on_house_stairs(0);
bool building_has_open_ext_door(0), sec_camera_shadow_mode(0), player_in_ww_elevator(0), player_in_skyway(0), player_on_moving_ww(0), player_on_escalator(0);
bool player_in_tunnel(0), player_in_mall(0);
int player_in_basement(0); // 0=no, 1=below ground level, 2=in basement and not on stairs, 3=in extended basement
int player_in_closet  (0); // uses flags RO_FLAG_IN_CLOSET (player in closet), RO_FLAG_LIT (closet light is on), RO_FLAG_OPEN (closet door is open)
int player_in_water   (0); // 0=no, 1=standing in water, 2=head underwater
int player_in_attic   (0); // 0=no, 1=attic with windows, 2=windowless attic
float building_bcube_expand(0.0), building_ambient_scale(0.0);
point player_candle_pos;
vector3d cur_camera_pos_xlate;
cube_t building_occluder;
building_params_t global_building_params;
building_t const *player_building(nullptr);
building_t const *vis_conn_bldg  (nullptr); // non-player building visible through extended basement connector room

extern bool start_in_inf_terrain, draw_building_interiors, flashlight_on, enable_use_temp_vbo, toggle_room_light;
extern bool teleport_to_screenshot, enable_dlight_bcubes, can_do_building_action, mirror_in_ext_basement;
extern unsigned room_mirror_ref_tid;
extern int rand_gen_index, display_mode, window_width, window_height, camera_surf_collide, animate2, building_action_key, player_in_elevator, frame_counter;
extern float CAMERA_RADIUS, fticks, NEAR_CLIP, FAR_CLIP;
extern colorRGB cur_ambient, cur_diffuse;
extern point pre_smap_player_pos, actual_player_pos;
extern cube_t smap_light_clip_cube;
extern vector<light_source> dl_sources;
extern vector<point> enabled_bldg_lights;
extern tree_placer_t tree_placer;
extern shader_t reflection_shader;


void bind_default_sun_moon_smap_textures();
void get_all_model_bcubes(vector<cube_t> &bcubes); // from model3d.h
cube_t get_building_indir_light_bounds(); // from building_lighting.cpp
float get_power_pole_height();
void register_player_not_in_building();
bool player_holding_lit_candle();
void parse_universe_name_str_tables();
void try_join_city_building_ext_basements(vect_building_t &buildings);
void add_sign_text_verts_both_sides(string const &text, cube_t const &sign, bool dim, bool dir, vect_vnctcc_t &verts);
void draw_candle_flames();
void update_security_camera_image();
void get_pedestrians_in_area(cube_t const &area, int building_ix, vector<point> &pts);
void setup_puddles_texture(shader_t &s);

float get_door_open_dist   () {return 3.5*CAMERA_RADIUS;}
bool player_in_ext_basement() {return (player_in_basement == 3 && player_building != nullptr);}

void tid_nm_pair_dstate_t::set_for_shader(float new_bump_map_mag) {
	if (new_bump_map_mag == bump_map_mag) return; // no change
	bump_map_mag = new_bump_map_mag;
	if (bmm_loc == -1) {bmm_loc = s.get_uniform_loc("bump_map_mag");} // set on the first call
	s.set_uniform_float(bmm_loc, bump_map_mag);
}
tid_nm_pair_dstate_t::~tid_nm_pair_dstate_t() { // restore to default if needed
	if (bmm_loc && bump_map_mag != 1.0) {s.set_uniform_float(bmm_loc, 1.0);} // bmm_loc should have been set
}

tid_nm_pair_t tid_nm_pair_t::get_scaled_version(float scale) const {
	tid_nm_pair_t tex(*this);
	tex.tscale_x *= scale;
	tex.tscale_y *= scale;
	return tex;
}
float tid_nm_pair_t::get_emissive_val() const {
	if (tid == RED_TEX) {return ((fract(tfticks/(1.5*TICKS_PER_SECOND)) < 0.5) ? 1.0 : 0.0);} // camera light flashes on and off with a period of 1.5s
	return emissive;
}
void tid_nm_pair_t::set_specular_color(colorRGB const &color, float mag, float shine) {
	if (shine == 0.0) {assert(color == WHITE);} // can't set zero shininess with a colored specular
	float max_comp(max(color.R, max(color.G, color.B)));
	if (max_comp == 0.0) {spec_color.set_c3(colorRGB(mag, mag, mag));} // black material has white specular; avoid divide-by-zero
	else {spec_color.set_c3(color*(mag/max_comp));} // extract color value normalized to largest component and multiply by mag; will cancel out any lighting
	shininess = (unsigned char)max(1, min(255, round_fp(shine)));
}
void tid_nm_pair_t::set_gl(tid_nm_pair_dstate_t &state) const {
	if (state.no_set_texture) {} // nothing to do
	else if (tid == FONT_TEXTURE_ID) {text_drawer::bind_font_texture();}
	else if (tid == REFLECTION_TEXTURE_ID) {
		if (bind_reflection_shader()) return;
	}
	else if (tid == NO_SHADOW_WHITE_TEX || tid == SHADOW_ONLY_TEX) {select_texture(WHITE_TEX);}
	else {select_texture(tid);}
	bool const has_normal_map(get_nm_tid() != FLAT_NMAP_TEX);
	if (has_normal_map) {select_texture(get_nm_tid(), 5);} // else we set bump_map_mag=0.0
	state.set_for_shader(has_normal_map ? 1.0 : 0.0); // enable or disable normal map (only ~25% of calls have a normal map)
	float const e_val(get_emissive_val());
	if (e_val     > 0.0) {state.s.add_uniform_float("emissive_scale", e_val);} // enable emissive
	if (shininess > 0  ) {state.s.set_specular_color(spec_color.get_c3(), shininess);} // colored specular
	if (no_cracks && state.crack_weight > 0.0) {state.s.add_uniform_float("crack_weight", 0.0);}
}
void tid_nm_pair_t::unset_gl(tid_nm_pair_dstate_t &state) const {
	if (tid == REFLECTION_TEXTURE_ID && room_mirror_ref_tid != 0) {state.s.make_current(); return;}
	bool const has_normal_map(get_nm_tid() != FLAT_NMAP_TEX);
	if (has_normal_map) {bind_default_flat_normal_map();} // reset back to flat normal map
	if (get_emissive_val() > 0.0) {state.s.add_uniform_float("emissive_scale", 0.0);} // disable emissive
	if (shininess          > 0  ) {state.s.clear_specular();} // clear specular
	if (no_cracks && state.crack_weight > 0.0) {state.s.add_uniform_float("crack_weight", state.crack_weight);} // restore original value
}
void tid_nm_pair_t::toggle_transparent_windows_mode() { // hack
	if      (tid == BLDG_WINDOW_TEX    ) {tid = BLDG_WIND_TRANS_TEX;}
	else if (tid == BLDG_WIND_TRANS_TEX) {tid = BLDG_WINDOW_TEX;}
}

void room_object_t::check_normalized() const {
	if (!is_strictly_normalized()) {std::cerr << "denormalized object of type " << unsigned(type) << " at " << str() << endl; assert(0);}
}

bool room_object_t::enable_rugs    () {return !global_building_params.rug_tids    .empty();}
bool room_object_t::enable_pictures() {return !global_building_params.picture_tids.empty();}

int select_tid_from_list(vector<unsigned> const &tids, unsigned ix) {return (tids.empty() ? -1 : tids[ix % tids.size()]);}
int room_object_t::get_rug_tid         () const {return select_tid_from_list(global_building_params.rug_tids,     obj_id);}
int room_object_t::get_picture_tid     () const {return select_tid_from_list(global_building_params.picture_tids, obj_id);}
int room_object_t::get_tv_tid          () const {return select_tid_from_list(global_building_params.picture_tids, obj_id/2);} // divide by 2 because even obj_id is turned off
int room_object_t::get_comp_monitor_tid() const {return select_tid_from_list(global_building_params.desktop_tids, obj_id/2);} // divide by 2 because even obj_id is turned off
int room_object_t::get_sheet_tid       () const {return select_tid_from_list(global_building_params.sheet_tids,   obj_id);}
int room_object_t::get_paper_tid       () const {return select_tid_from_list(global_building_params.paper_tids,   obj_id);}
int room_object_t::get_food_box_tid    () const {return select_tid_from_list(global_building_params.food_box_tids,obj_id);}
int get_metal_texture(unsigned id)              {return select_tid_from_list(global_building_params.metal_tids,       id);}
int get_flag_texture (unsigned id)              {return select_tid_from_list(global_building_params.flag_tids,        id);} // food_box_names

string const &select_str_from_list(vector<string> const &strs, unsigned ix) {
	static string empty_str;
	return (strs.empty() ? empty_str : strs[ix % strs.size()]);
}
string const &room_object_t::get_food_box_name() const {return select_str_from_list(global_building_params.food_box_names, obj_id);}

void do_xy_rotate(float rot_sin, float rot_cos, point const &center, point &pos) {
	float const x(pos.x - center.x), y(pos.y - center.y); // translate to center
	pos.x = x*rot_cos - y*rot_sin + center.x;
	pos.y = y*rot_cos + x*rot_sin + center.y;
}
void do_xy_rotate_normal(float rot_sin, float rot_cos, point &pos) { // point rotate without the translate
	float const x(pos.x), y(pos.y);
	pos.x = x*rot_cos - y*rot_sin;
	pos.y = y*rot_cos + x*rot_sin;
}
void building_geom_t::do_xy_rotate    (point const &center, point &pos) const {::do_xy_rotate( rot_sin, rot_cos, center, pos);}
void building_geom_t::do_xy_rotate_inv(point const &center, point &pos) const {::do_xy_rotate(-rot_sin, rot_cos, center, pos);}
void building_geom_t::do_xy_rotate_normal    (point &n) const {::do_xy_rotate_normal( rot_sin, rot_cos, n);}
void building_geom_t::do_xy_rotate_normal_inv(point &n) const {::do_xy_rotate_normal(-rot_sin, rot_cos, n);}


class building_texture_mgr_t {
	int window_tid=-1, hdoor_tid=-1, odoor_tid=-1, bdoor_tid=-1, bdoor2_tid=-1, gdoor_tid=-1, mdoor_tid=-1, ac_unit_tid1=-1, ac_unit_tid2=-1, bath_wind_tid=-1, helipad_tid=-1,
		solarp_tid=-1, concrete_tid=-1, met_plate_tid=-1, mplate_nm_tid=-1, met_roof_tid=-1, tile_floor_tid=-1, tile_floor_nm_tid=-1, duct_tid=-1, vent_tid=-1;

	int ensure_tid(int &tid, const char *name, bool is_normal_map=0, bool invert_y=0) {
		if (tid < 0) {tid = get_texture_by_name(name, is_normal_map, invert_y);}
		if (tid < 0) {tid = (is_normal_map ? FLAT_NMAP_TEX : WHITE_TEX);} // failed to load texture - use a simple white texture/flat normal map
		return tid;
	}
public:
	int get_window_tid   () const {return window_tid;}
	int get_hdoor_tid    () {return ensure_tid(hdoor_tid,     "white_door.jpg");} // house door
	int get_odoor_tid    () {return ensure_tid(odoor_tid,     "buildings/office_door.jpg");} // office door
	int get_bdoor_tid    () {return ensure_tid(bdoor_tid,     "buildings/building_door.jpg");} // metal + glass building door
	int get_bdoor2_tid   () {return ensure_tid(bdoor2_tid,    "buildings/metal_door.jpg");} // metal building door
	int get_gdoor_tid    () {return ensure_tid(gdoor_tid,     "buildings/garage_door.jpg");} // garage door
	int get_mdoor_tid    () {return ensure_tid(mdoor_tid,     "buildings/modern_door.jpg");} // unused; for future use, maybe with house exterior doors
	int get_ac_unit_tid1 () {return ensure_tid(ac_unit_tid1,  "buildings/AC_unit1.jpg");} // AC unit (should this be a <d> loop?)
	int get_ac_unit_tid2 () {return ensure_tid(ac_unit_tid2,  "buildings/AC_unit2.jpg");} // AC unit
	int get_duct_tid     () {return ensure_tid(duct_tid,      "interiors/duct.jpg");} // duct
	int get_vent_tid     () {return ensure_tid(vent_tid,      "interiors/vent.jpg");} // vent
	int get_bath_wind_tid() {return ensure_tid(bath_wind_tid, "buildings/window_blocks.jpg");} // bathroom window
	int get_helipad_tid  () {return ensure_tid(helipad_tid,   "buildings/helipad.jpg");}
	int get_solarp_tid   () {return ensure_tid(solarp_tid,    "buildings/solar_panel.jpg");}
	int get_concrete_tid () {return ensure_tid(concrete_tid,  "roads/concrete.jpg");}
	int get_met_plate_tid() {return ensure_tid(met_plate_tid, "metal_plate.jpg");}
	int get_mplate_nm_tid() {return ensure_tid(mplate_nm_tid, "normal_maps/metal_plate_NRM.jpg", 1);} // is_normal_map=1
	int get_met_roof_tid () {return ensure_tid(met_roof_tid,  "buildings/metal_roof.jpg");}
	int get_tile_floor_tid   () {return ensure_tid(tile_floor_tid,    "interiors/mosaic_tiles.jpg");}
	int get_tile_floor_nm_tid() {return ensure_tid(tile_floor_nm_tid, "interiors/mosaic_tiles_normal.jpg");}

	bool check_windows_texture() {
		if (!global_building_params.windows_enabled()) return 0;
		if (window_tid >= 0) return 1; // already generated
		gen_building_window_texture(global_building_params.get_window_width_fract(), global_building_params.get_window_height_fract());
		window_tid = BLDG_WINDOW_TEX;
		return 1;
	}
	bool is_door_tid(int tid) const {return (tid >= 0 && (tid == hdoor_tid || tid == odoor_tid || tid == bdoor_tid || tid == bdoor2_tid || tid == gdoor_tid || tid == mdoor_tid));}
};
building_texture_mgr_t building_texture_mgr;

int get_rect_panel_tid() {return building_texture_mgr.get_gdoor_tid    ();} // use garage doors
int get_bath_wind_tid () {return building_texture_mgr.get_bath_wind_tid();}
int get_int_door_tid  () {return building_texture_mgr.get_hdoor_tid    ();}
int get_bldg_door_tid () {return building_texture_mgr.get_bdoor_tid    ();}
int get_off_door_tid  () {return building_texture_mgr.get_odoor_tid    ();}
int get_concrete_tid  () {return building_texture_mgr.get_concrete_tid ();}
int get_solarp_tid    () {return building_texture_mgr.get_solarp_tid   ();}
int get_met_plate_tid () {return building_texture_mgr.get_met_plate_tid();}
int get_mplate_nm_tid () {return building_texture_mgr.get_mplate_nm_tid();}
int get_ac_unit_tid   (unsigned ix) {return ((ix & 1) ? building_texture_mgr.get_ac_unit_tid1() : building_texture_mgr.get_ac_unit_tid2());}

void set_tile_floor_texture() {
	select_texture(building_texture_mgr.get_tile_floor_tid   ());
	select_texture(building_texture_mgr.get_tile_floor_nm_tid(), 5);
}


class texture_id_mapper_t {
	vector<unsigned> tid_to_slot_ix;
	vector<int> tid_to_nm_tid;
	set<unsigned> ext_wall_tids, roof_tids;
	unsigned next_slot_ix;

	void register_tid(int tid) {
		if (tid < 0) return; // not allocated
		if (tid >= (int)tid_to_slot_ix.size()) {tid_to_slot_ix.resize(tid+1, 0);}
		if (tid_to_slot_ix[tid] == 0) {tid_to_slot_ix[tid] = next_slot_ix++;}
		//cout << "register " << tid << " slot " << tid_to_slot_ix[tid] << endl;
	}
	void register_tex(tid_nm_pair_t const &tex) {
		register_tid(tex.tid);

		if (tex.tid > 0 && tex.nm_tid > 0 && tex.nm_tid != FLAT_NMAP_TEX) {
			if (tex.tid >= (int)tid_to_nm_tid.size()) {tid_to_nm_tid.resize(tex.tid+1, -1);}
			tid_to_nm_tid[tex.tid] = tex.nm_tid;
		}
	}
public:
	texture_id_mapper_t() : next_slot_ix(1) {} // slots start at 1; slot 0 is for untextured

	void init() {
		if (!tid_to_slot_ix.empty()) return; // already inited
		// register all textures that will be used here, before we get into the OMP parallel block
		unsigned const num_special_tids = 7;
		int const special_tids[num_special_tids] = {WHITE_TEX, NO_SHADOW_WHITE_TEX, SHADOW_ONLY_TEX, FENCE_TEX, PANELING_TEX, TILE_TEX, WOOD_TEX}; // for elevators, etc.
		tid_to_slot_ix.push_back(0); // untextured case
		register_tid(building_texture_mgr.get_window_tid());
		register_tid(building_texture_mgr.get_hdoor_tid());
		register_tid(building_texture_mgr.get_odoor_tid());
		register_tid(building_texture_mgr.get_bdoor_tid());
		register_tid(building_texture_mgr.get_bdoor2_tid());
		register_tid(building_texture_mgr.get_gdoor_tid());
		//register_tid(building_texture_mgr.get_mdoor_tid()); // enable when this door type is used
		register_tid(building_texture_mgr.get_ac_unit_tid1());
		register_tid(building_texture_mgr.get_ac_unit_tid2());
		register_tid(building_texture_mgr.get_duct_tid());
		register_tid(building_texture_mgr.get_vent_tid());
		register_tid(building_texture_mgr.get_helipad_tid());
		register_tid(building_texture_mgr.get_solarp_tid());
		register_tid(building_texture_mgr.get_concrete_tid());
		register_tid(building_texture_mgr.get_met_plate_tid());
		register_tid(building_texture_mgr.get_mplate_nm_tid());
		register_tid(building_texture_mgr.get_met_roof_tid());
		register_tid(building_texture_mgr.get_tile_floor_tid());
		register_tid(building_texture_mgr.get_tile_floor_nm_tid());
		register_tid(get_plywood_tid()); // for attics
		register_tid(FONT_TEXTURE_ID); // for roof signs
		for (unsigned i = 0; i < num_special_tids; ++i) {register_tid(special_tids[i]);}

		for (auto i = global_building_params.materials.begin(); i != global_building_params.materials.end(); ++i) {
			register_tex(i->side_tex);
			register_tex(i->roof_tex);
			register_tex(i->wall_tex);
			register_tex(i->ceil_tex);
			register_tex(i->floor_tex);
			register_tex(i->house_ceil_tex);
			register_tex(i->house_floor_tex);
			ext_wall_tids.insert(i->side_tex.tid);
			roof_tids    .insert(i->roof_tex.tid);
		} // for i
		cout << "Used " << (next_slot_ix-1) << " slots for texture IDs up to " << (tid_to_slot_ix.size()-1) << endl;
	}
	unsigned get_slot_ix(int tid) const {
		if (tid < 0) return 0; // untextured - slot 0
		assert(tid < (int)get_num_slots());
		assert(tid_to_slot_ix[tid] > 0);
		return tid_to_slot_ix[tid];
	}
	int get_slot_ix_if_exists(int tid) const { // returns -1 if tid is not found
		if (tid < 0) return 0; // untextured - slot 0
		if (tid >= (int)get_num_slots()) return -1; // not found
		if (tid_to_slot_ix[tid] == 0)    return -1; // empty slot
		return tid_to_slot_ix[tid];
	}
	int get_normal_map_for_tid(int tid) const {
		if (tid < 0 || (unsigned)tid >= tid_to_nm_tid.size()) return -1; // no normal map
		return tid_to_nm_tid[tid];
	}
	unsigned get_num_slots() const {return tid_to_slot_ix.size();}
	bool is_ext_wall_tid(unsigned tid) const {return (ext_wall_tids.find(tid) != ext_wall_tids.end());}
	bool is_roof_tid    (unsigned tid) const {return (roof_tids    .find(tid) != roof_tids    .end());}
};
texture_id_mapper_t tid_mapper;

int get_normal_map_for_bldg_tid(int tid) {return tid_mapper.get_normal_map_for_tid(tid);}

class tid_vert_counter_t {
	vector<unsigned> counts;
public:
	tid_vert_counter_t() {counts.resize(tid_mapper.get_num_slots(), 0);} // resized to max tid
	void update_count(int tid, unsigned num) {
		if (tid < 0) return;
		assert((unsigned)tid < counts.size());
		counts[tid] += num;
	}
	unsigned get_count(int tid) const {
		if (tid < 0) return 0;
		assert((unsigned)tid < counts.size());
		return counts[tid];
	}
};


class indir_tex_mgr_t {
	unsigned tid; // Note: owned by building_indir_light_mgr, not us
public:
	indir_tex_mgr_t() : tid(0) {}
	bool enabled() const {return (tid > 0);}

	bool create_for_building(building_t const &b, unsigned bix, point const &target) {
		b.create_building_volume_light_texture(bix, target, tid);
		return 1;
	}
	bool setup_for_building(shader_t &s) const {
		if (!enabled()) return 0; // no texture set
		cube_t const lighting_bcube(get_building_indir_light_bounds());
		float const dx(lighting_bcube.dx()/MESH_X_SIZE), dy(lighting_bcube.dy()/MESH_Y_SIZE), dxy_offset(0.5f*(dx + dy));
		bind_texture_tu(tid, 1); // indir texture uses TU_ID=1
		s.add_uniform_vector3d("alt_scene_llc",   lighting_bcube.get_llc());
		s.add_uniform_vector3d("alt_scene_scale", lighting_bcube.get_size());
		s.add_uniform_float("half_dxy", dxy_offset);
		return 1;
	}
};
indir_tex_mgr_t indir_tex_mgr;

bool player_in_dark_room() {return (player_in_unlit_room || (player_in_closet && !(player_in_closet & (RO_FLAG_OPEN | RO_FLAG_LIT))));}

struct building_lights_manager_t : public city_lights_manager_t {
	void setup_building_lights(vector3d const &xlate, bool sec_camera_mode=0) {
		//highres_timer_t timer("Building Dlights Setup"); // 1.9/1.9
		float const light_radius(0.1*light_radius_scale*get_tile_smap_dist()); // distance from the camera where lights are drawn
		if (!begin_lights_setup(xlate, light_radius, dl_sources)) return;
		// include the building and it's extended basement and underground rooms in the lights_bcube; needed for malls
		if (player_building != nullptr) {lights_bcube.union_with_cube_xy(player_building->get_bcube_inc_extensions());}
		// no room lights if player is hiding in a closed closet/windowless room with light off (prevents light leakage)
		if (sec_camera_mode || !player_in_dark_room()) {add_building_interior_lights(xlate, lights_bcube, sec_camera_mode);}
		if (flashlight_on && !sec_camera_mode) {add_player_flashlight(0.12);} // add player flashlight, even when outside of building so that flashlight can shine through windows
		if (camera_in_building && !sec_camera_mode && player_holding_lit_candle()) {add_player_candle_light(xlate);}
		clamp_to_max_lights(xlate, dl_sources);
		tighten_light_bcube_bounds(dl_sources); // clip bcube to tight bounds around lights for better dlights texture utilization (possible optimization)
		
		if (ADD_ROOM_SHADOWS) {
			sec_camera_shadow_mode = sec_camera_mode; // optimization
			setup_shadow_maps(dl_sources, (camera_pdu.pos - xlate), global_building_params.max_shadow_maps, sec_camera_mode);
			sec_camera_shadow_mode = 0; // restore
		}
		finalize_lights(dl_sources);
	}
	void add_player_candle_light(vector3d const &xlate) {
		static float cval(0.5), inten(0.75);
		float const radius(10.0*CAMERA_RADIUS); // based on floor spacing?
		point const pos((player_candle_pos == all_zeros) ? (camera_pdu.pos - xlate) : player_candle_pos); // use player_candle_pos if valid, otherwise camera pos
		dl_sources.emplace_back(inten*radius, pos, pos, gen_fire_color(cval, inten, 1.0));
		dl_sources.back().disable_shadows(); // shadows not needed / not valid for point lights
		min_eq(lights_bcube.z1(), (pos.z - radius));
		max_eq(lights_bcube.z2(), (pos.z + radius));
	}
	virtual bool enable_lights() const {return (draw_building_interiors || flashlight_on);}
};

building_lights_manager_t building_lights_manager;

void setup_building_lights(vector3d const &xlate, bool sec_camera_mode=0) {
	interior_shadow_maps = 1; // set state so that above call will know that it was called recursively from here and should draw interior shadow maps
	enable_dlight_bcubes = 1; // needed around this call so that light bcubes are sent to the GPU
	building_lights_manager.setup_building_lights(xlate, sec_camera_mode);
	enable_dlight_bcubes = 0; // disable when creating the reflection image (will be set when we re-enter multi_draw())
	interior_shadow_maps = 0;
}

void interpolate_over_time(float &val, float target_val, float transition_secs, int &last_frame) {
	if (frame_counter == last_frame) return; // update once per frame
	last_frame = frame_counter;
	float const delta_val(fticks/(transition_secs*TICKS_PER_SECOND));
	if      (val > target_val) {val = max(target_val, (val - delta_val));} // decrease
	else if (val < target_val) {val = min(target_val, (val + delta_val));} // increase
}
void set_interior_lighting(shader_t &s, bool have_indir) {
	bool const player_in_mall(player_in_basement == 3 && player_building && player_building->has_mall());
	float const light_scale(0.5), target_blscale((player_in_basement || player_in_attic) ? (player_in_mall ? 0.8 : 0.0) : (player_in_walkway ? 2.0 : 1.0));
	float const target_ascale(player_in_tunnel ? 0.0 : 1.0); // brighter ambient unless in tunnel
	static float blscale(1.0), ascale(1.0);
	static int lu_frame1(0), lu_frame2(0);
	interpolate_over_time(blscale, target_blscale, 0.5, lu_frame1); // indir/ambient lighting slowly transitions when entering or leaving the basement or walkway
	interpolate_over_time(ascale,  target_ascale,  0.5, lu_frame2);
	float ambient_scale(0.5f*(ascale + blscale)*light_scale);
	float diffuse_scale(0.2f*blscale*light_scale); // reduce diffuse and specular lighting for sun/moon
	float hemi_scale(   0.2f*blscale*light_scale); // reduced hemispherical lighting

	if (have_indir || player_in_dark_room()) { // using indir lighting, or player in a closed closet/windowless room with the light off
		s.add_uniform_float("SHADOW_LEAKAGE", 0.0); // no light leakage
		
		if (have_indir) { // set ambient color to use with indir lookups outside the current building
			// since we can't add proper diffuse, make 50% of diffuse the ambient color assuming 50% of surfaces are diffusely lit
			s.add_uniform_color("out_range_indir_color", (cur_ambient*ambient_scale + cur_diffuse*(0.5*diffuse_scale)));
		}
		ambient_scale = ((!have_indir && player_in_dark_room()) ? 0.1 : 0.0); // no ambient for indir; slight ambient for closed closet/windowless room with light off
		diffuse_scale = hemi_scale = 0.0; // no diffuse or hemispherical from sun/moon
	}
	else if (player_in_basement && !player_in_mall) {
		s.add_uniform_float("SHADOW_LEAKAGE", 0.0); // make basements darker and avoid lights leaking through parking garage ceilings
	}
	s.add_uniform_float("diffuse_scale",       diffuse_scale);
	s.add_uniform_float("ambient_scale",       ambient_scale);
	s.add_uniform_float("hemi_lighting_scale", hemi_scale);
	building_ambient_scale = ambient_scale; // cache so that we can reset back to this value when drawing bubbles, etc.
}
void reset_interior_lighting(shader_t &s) {
	s.add_uniform_float("diffuse_scale",       1.0 ); // re-enable diffuse and specular lighting for sun/moon
	s.add_uniform_float("ambient_scale",       1.0 ); // reset to default
	s.add_uniform_float("hemi_lighting_scale", 0.5 ); // reset to default
	s.add_uniform_float("SHADOW_LEAKAGE",      0.05); // reset to default
}
void reset_interior_lighting_and_end_shader(shader_t &s) {
	reset_interior_lighting(s);
	s.end_shader();
}
bool have_building_indir_lighting() {
	return indir_tex_mgr.enabled() && enable_building_indir_lighting();
}
void setup_building_draw_shader_post(shader_t &s, bool have_indir) {
	set_interior_lighting(s, have_indir);
	if (have_indir) {indir_tex_mgr.setup_for_building(s);}
}
void setup_building_draw_shader(shader_t &s, float min_alpha, bool enable_indir, bool force_tsl, int use_texgen, float water_damage, float crack_damage) { // for building interiors
	float const pcf_scale = 0.2;
	if (player_building == nullptr) {water_damage = crack_damage = 0.0;} // water damage and cracks only apply to player building; this can fail on the exterior walls pass
	// disable indir if the player is in a closed closet
	bool const have_indir(enable_indir && have_building_indir_lighting() && !(player_in_closet && !(player_in_closet & RO_FLAG_OPEN)));
	bool const add_vorocracks(global_building_params.use_voronoise_cracks && crack_damage > 0.0);
	int const use_bmap(global_building_params.has_normal_map), interior_use_smaps(ADD_ROOM_SHADOWS ? 2 : 1); // dynamic light smaps only
	cube_t const lights_bcube(building_lights_manager.get_lights_bcube());
	if (have_indir) {s.set_prefix("#define ENABLE_OUTSIDE_INDIR_RANGE",  1);} // FS
	if (water_damage > 0.0) {s.set_prefix("#define ENABLE_WATER_DAMAGE", 1);} // FS
	if (crack_damage > 0.0) {s.set_prefix("#define ADD_CRACKS",          1);} // FS
	if (add_vorocracks    ) {s.set_prefix("#define USE_VOROCRACKS",      1);} // FS
	s.set_prefix("#define LINEAR_DLIGHT_ATTEN", 1); // FS; improves room lighting (better light distribution vs. framerate trade-off)
	city_shader_setup(s, lights_bcube, 1, interior_use_smaps, use_bmap, min_alpha, force_tsl, pcf_scale, use_texgen, have_indir, 0); // use_dlights=1, is_outside=0
	setup_building_draw_shader_post(s, have_indir);
	if (water_damage > 0.0 || crack_damage > 0.0) {setup_puddles_texture(s);} // 3D texture is used for both water damage and cracks

	if (water_damage > 0.0) {
		// Note: applies to basements only; needed for player building, but interiors are drawn by tile, and the other building basements aren't visible anyway
		s.add_uniform_float("wet_effect",   water_damage);
		s.add_uniform_float("puddle_scale", 0.5);
		s.add_uniform_float("water_damage_zmax", player_building->ground_floor_z1); // water damage is only in the basement
		s.add_uniform_float("water_damage_zscale", 0.25); // stretch out vertically on walls
	}
	if (crack_damage > 0.0) {
		s.add_uniform_float("crack_scale", 1.0); // Note: crack_weight will be set to crack_damage later
		s.add_uniform_float("crack_zmax",  player_building->ground_floor_z1); // cracks are only in the basement
		// disable cracks on on ceilings (-z) since they may be wood; carpet is special cased to not have cracks; are cracks on particle board ceilings okay?
		s.add_uniform_float("crack_normal_zmax", (player_building->is_house ? -0.5 : -2.0));
	}
}


/*static*/ void building_draw_utils::calc_normals(building_geom_t const &bg, vector<vector3d> &nv, unsigned ndiv) {

	assert(bg.flat_side_amt >= 0.0 && bg.flat_side_amt < 0.5); // generates a flat side
	assert(bg.alt_step_factor >= 0.0 && bg.alt_step_factor < 1.0);
	if (bg.flat_side_amt > 0.0) {assert(ndiv > 4);} // should be at least 5 sides, 6-8 is better
	float const ndiv_inv(1.0/ndiv), css(TWO_PI*ndiv_inv*(1.0f - bg.flat_side_amt));
	float sin_ds[2] = {}, cos_ds[2] = {};

	if (bg.alt_step_factor > 0.0) { // alternate between large and small steps (cube with angled corners, etc.)
		assert(!(ndiv&1));
		float const css_v[2] = {css*(1.0f + bg.alt_step_factor), css*(1.0f - bg.alt_step_factor)};
		UNROLL_2X(sin_ds[i_] = sin(css_v[i_]); cos_ds[i_] = cos(css_v[i_]);)
	}
	else { // uniform side length
		sin_ds[0] = sin_ds[1] = sin(css);
		cos_ds[0] = cos_ds[1] = cos(css);
	}
	float sin_s(0.0), cos_s(1.0), angle0(bg.start_angle); // start at 0.0
	if (bg.half_offset) {angle0 = 0.5*css;} // for cube
	if (angle0 != 0.0) {sin_s = sin(angle0); cos_s = cos(angle0);} // uncommon case
	nv.resize(ndiv);

	for (unsigned S = 0; S < ndiv; ++S) { // build normals table
		bool const d(S&1);
		float const s(sin_s), c(cos_s);
		nv[S].assign(s, c, 0.0);
		sin_s = s*cos_ds[d] + c*sin_ds[d];
		cos_s = c*cos_ds[d] - s*sin_ds[d];
	}
}

/*static*/ void building_draw_utils::calc_poly_pts(building_geom_t const &bg, cube_t const &bcube, cube_t const &part, vect_point &pts) {

	calc_normals(bg, pts, bg.num_sides);
	vector3d const sz(part.get_size());
	point const cc(part.get_cube_center());
	float const rx(0.5*sz.x), ry(0.5*sz.y);

	if (bg.is_rotated() && part != bcube) {
		// the building is rotated around the bcube center, but the part itself is rotated around its own center, so we have to adjust the points correctly
		point const rot_pos(part.get_cube_center()), inv_rot_pos(bcube.get_cube_center());

		for (point &pt : pts) {
			pt.assign((cc.x + rx*pt.x), (cc.y + ry*pt.y), 0.0); // convert normals to points
			bg.do_xy_rotate(rot_pos, pt);
			bg.do_xy_rotate_inv(inv_rot_pos, pt);
		}
	}
	else {
		for (point &pt : pts) {pt.assign((cc.x + rx*pt.x), (cc.y + ry*pt.y), 0.0);} // convert normals to points
	}
}

// Note: invert_tc only applies to doors
void add_tquad_to_verts(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex,
	colorRGBA const &color, vect_vnctcc_t &verts, bool invert_tc_x, bool exclude_frame, bool no_tc, bool no_rotate, bool swap_tc_xy)
{
	assert(tquad.npts == 3 || tquad.npts == 4); // triangles or quads
	bool const do_rotate(bg.is_rotated() && !no_rotate);
	point const center(do_rotate ? bcube.get_cube_center() : all_zeros); // rotate about bounding cube / building center
	vert_norm_comp_tc_color vert;
	float tsx(0.0), tsy(0.0), tex_off(0.0);
	bool dim(0);

	if (tquad.type == tquad_with_ix_t::TYPE_WALL) { // side/wall
		tsx = tex.get_drawn_tscale_x(); tsy = tex.get_drawn_tscale_y();
		dim = (tquad.pts[0].x == tquad.pts[1].x);
		if (world_mode != WMODE_INF_TERRAIN) {tex_off = (dim ? yoff2*DY_VAL : xoff2*DX_VAL);}
		tex_off -= (dim ? bcube.y1() : bcube.x1()); // normalize to building LLC to keep tex coords small
	}
	else if (tquad.is_roof() || tquad.type == tquad_with_ix_t::TYPE_ROOF_ACC) { // roof cap
		float const denom(0.5f*(bcube.dx() + bcube.dy()));
		tsx = tex.tscale_x/denom; tsy = tex.tscale_y/denom;
	}
	vert.set_c4(color);
	vector3d normal(tquad.get_norm());
	if (do_rotate) {bg.do_xy_rotate_normal(normal);}
	vert.set_norm(normal);
	invert_tc_x   ^= tquad.is_inside_face(); // invert interior/office door, inner/back face
	exclude_frame &= (tquad.is_interior_door() || tquad.is_exterior_door());

	for (unsigned i = 0; i < tquad.npts; ++i) {
		vert.v = tquad.pts[i];

		if (no_tc) { // untextured, for door edges
			vert.t[0] = vert.t[1] = 0.0;
		}
		else if (tquad.type == tquad_with_ix_t::TYPE_WALL) { // side/wall
			vert.t[0] = (vert.v[dim] + tex_off)*tsx; // use nonzero width dim
			vert.t[1] = (vert.v.z - bcube.z1())*tsy;
		}
		else if (tquad.is_roof()) { // roof cap
			vert.t[0] = (vert.v.x - bcube.x1())*tsx; // varies from 0.0 and bcube x1 to 1.0 and bcube x2
			vert.t[1] = (vert.v.y - bcube.y1())*tsy; // varies from 0.0 and bcube y1 to 1.0 and bcube y2
		}
		else if (tquad.type == tquad_with_ix_t::TYPE_ROOF_ACC) { // roof access cover
			if (fabs(normal.z) > 0.5) {vert.t[0] = vert.v.x*tsx; vert.t[1] = vert.v.y*tsy;} // facing up, use XY plane
			else {vert.t[0] = (vert.v.x + vert.v.y)*tsx; vert.t[1] = vert.v.z*tsy;} // facing to the side, use XZ or YZ plane
		}
		else if (tquad.is_exterior_door() || tquad.type == tquad_with_ix_t::TYPE_HELIPAD || tquad.type == tquad_with_ix_t::TYPE_SOLAR) { // textured from (0,0) to (1,1)
			vert.t[0] = float((i == 1 || i == 2) ^ invert_tc_x);
			vert.t[1] = float((i == 2 || i == 3));
			if      (tquad.type == tquad_with_ix_t::TYPE_SOLAR) {vert.t[0] *= 4.0; vert.t[1] *= 4.0;} // 4 reptitions in each dimension
			else if (tquad.is_rooftop_door()) { // only draw half of the door for rooftop doors; slightly more to pick up the frame if closed
				vert.t[0] *= ((!exclude_frame && tquad.type == tquad_with_ix_t::TYPE_RDOOR) ? 0.52 : 0.5);
				if (exclude_frame) {vert.t[1] *= 0.97;} // trim off the top door frame
			}
		}
		else if (tquad.is_interior_door()) { // interior door textured/stretched in Y
			vert.t[0]  = tex.tscale_x*((i == 1 || i == 2) ^ invert_tc_x);
			vert.t[1]  = tex.tscale_y*((i == 2 || i == 3));
			vert.t[1] *= 0.97; // trim off the top door frame
		}
		else if (tquad.type == tquad_with_ix_t::TYPE_TRIM) {} // untextured - no tex coords
		else {assert(0);}
		if (exclude_frame) {vert.t[0] = DOOR_FRAME_WIDTH + (1.0 - 2.0*DOOR_FRAME_WIDTH)*vert.t[0];}
		if (do_rotate) {bg.do_xy_rotate(center, vert.v);}
		if (swap_tc_xy) {swap(vert.t[0], vert.t[1]);}
		verts.push_back(vert);
	} // for i
}


#define EMIT_VERTEX() \
	vert.v = pt*sz + llc; \
	vert.t[ st] = tscale[ st]*(vert.v[d] + tex_vert_off[d]); \
	vert.t[!st] = tscale[!st]*(vert.v[i] + tex_vert_off[i]); \
	vert.t[0] += tex.txoff; \
	vert.t[1] += tex.tyoff; \
	if (apply_ao) {vert.copy_color(cw[pt.z > 0.5]);} \
	if (bg.is_rotated()) {bg.do_xy_rotate(center, vert.v);} \
	verts.push_back(vert);

#define EMIT_VERTEX_SIMPLE() \
	vert.v = pt*sz + llc; \
	vert.t[ st] = tscale[ st]*(ws_texture ? vert.v[d] : pt[d]); \
	vert.t[!st] = tscale[!st]*(ws_texture ? vert.v[i] : pt[i]); \
	verts.push_back(vert);

class building_draw_t {

	static vbo_cache_t vbo_cache; // shared across all bdraws/tiles/blocks/buildings

	class draw_block_t {
		struct vert_ix_pair {
			unsigned qix, tix; // {quads, tris}
			vert_ix_pair(unsigned qix_, unsigned tix_) : qix(qix_), tix(tix_) {}
			bool operator==(vert_ix_pair const &v) const {return (qix == v.qix && tix == v.tix);}
		};
		indexed_vao_manager_with_shadow_t vao_mgr; // Note: not using the indexed part
		vector<vert_ix_pair> pos_by_tile; // {quads, tris}
		unsigned tri_vbo_off, vert_vbo_sz;
		unsigned start_num_verts[2] = {0}; // for quads and triangles
	public:
		bool no_shadows;
		tid_nm_pair_t tex;
		vect_vnctcc_t quad_verts, tri_verts;

		draw_block_t() : tri_vbo_off(0), vert_vbo_sz(0), no_shadows(0) {}
		void record_num_verts() {start_num_verts[0] = num_quad_verts(); start_num_verts[1] = num_tri_verts();}

		void draw_geom_range(tid_nm_pair_dstate_t &state, bool shadow_only, vert_ix_pair const &vstart, vert_ix_pair const &vend) { // use VBO rendering
			if (vstart == vend)                  return; // empty range - no verts for this tile
			if (shadow_only && no_shadows)       return; // no shadows on this material
			if (!shadow_only && tex.shadow_only) return; // material is only drawn in the shadow pass
			int depth_write_disabled(0);
			
			if (tex.tid == FONT_TEXTURE_ID) {
				if (shadow_only) return; // no shadows for text
				enable_blend();
				glGetIntegerv(GL_DEPTH_WRITEMASK, &depth_write_disabled);
				if (depth_write_disabled) {glDepthMask(GL_FALSE);} // disable depth writing if it was enabled
			}
			if (!shadow_only) {tex.set_gl(state);}
			assert(vao_mgr.vbo_valid());
			vao_mgr.create_from_vbo<vert_norm_comp_tc_color>(shadow_only, 1, 1); // setup_pointers=1, always_bind=1

			if (vstart.qix != vend.qix) { // usually this is nonempty
				assert(vstart.qix < vend.qix);
				draw_quads_as_tris((vend.qix - vstart.qix), vstart.qix);
			}
			if (vstart.tix != vend.tix) { // this is empty over half the time; merging this with the quads draw call likely has little runtime effect
				assert(vstart.tix < vend.tix);
				glDrawArrays(GL_TRIANGLES, (vstart.tix + tri_vbo_off), (vend.tix - vstart.tix));
				++num_frame_draw_calls;
			}
			if (tex.tid == FONT_TEXTURE_ID) {
				if (depth_write_disabled) {glDepthMask(GL_TRUE);} // re-enable depth writing if needed
				disable_blend();
			}
			if (!shadow_only) {tex.unset_gl(state);}
			vao_manager_t::post_render();
		}
		void draw_all_geom(tid_nm_pair_dstate_t &state, bool shadow_only, bool direct_draw_no_vbo, vertex_range_t const *const exclude=nullptr) {
			if (shadow_only && no_shadows) return; // no shadows on this material

			if (direct_draw_no_vbo) {
				enable_use_temp_vbo = 1; // hack to fix missing wall artifacts when not using a core context
				assert(!exclude); // not supported in this mode
				bool const use_texture(!shadow_only && (!quad_verts.empty() || !tri_verts.empty()));
				if (use_texture) {tex.set_gl(state);} // Note: colors are not disabled here
				if (!quad_verts.empty()) {draw_quad_verts_as_tris(quad_verts, 0, 1, 1);}
				if (!tri_verts .empty()) {draw_verts(tri_verts, GL_TRIANGLES, 0, 1);}
				if (use_texture) {tex.unset_gl(state);}
				enable_use_temp_vbo = 0;
			}
			else {
				if (pos_by_tile.empty()) return; // nothing to draw for this block/texture
				vert_ix_pair const &start(pos_by_tile.front()), end(pos_by_tile.back());
				if (!exclude) {draw_geom_range(state, shadow_only, start, end); return;} // non-exclude case
				assert(exclude->start >= start.qix && exclude->start < exclude->end && exclude->end <= end.qix); // exclude (start, end) must be a subset of (start.qix, end.qix)
				draw_geom_range(state, shadow_only, start, vert_ix_pair(exclude->start, end.tix)); // first block of quads and all tris
				draw_geom_range(state, shadow_only, vert_ix_pair(exclude->end, end.tix), end); // second block of quads and no tris
			}
		}
		void draw_quad_geom_range(tid_nm_pair_dstate_t &state, vertex_range_t const &range, bool shadow_only=0) { // no tris; empty range is legal
			draw_geom_range(state, shadow_only, vert_ix_pair(range.start, 0), vert_ix_pair(range.end, 0));
		}
		void draw_tri_geom_range(tid_nm_pair_dstate_t &state, vertex_range_t const &range, bool shadow_only=0) { // no quads; empty range is legal
			draw_geom_range(state, shadow_only, vert_ix_pair(0, range.start), vert_ix_pair(0, range.end));
		}
		void draw_geom_tile(tid_nm_pair_dstate_t &state, unsigned tile_id, bool shadow_only) {
			if (pos_by_tile.empty()) return; // nothing to draw for this block/texture
			assert(tile_id+1 < pos_by_tile.size()); // tile and next tile must be valid indices
			draw_geom_range(state, shadow_only, pos_by_tile[tile_id], pos_by_tile[tile_id+1]); // shadow_only=0
		}
		void upload_to_vbos() {
			assert((num_quad_verts()%4) == 0);
			assert((num_tri_verts ()%3) == 0);
			tri_vbo_off = quad_verts.size(); // triangles start after quads
			vector_add_to(tri_verts, quad_verts);
			clear_cont(tri_verts); // no longer needed
			
			if (!quad_verts.empty()) {
				assert(!vao_mgr.vbo_valid());
				unsigned const verts_sz(quad_verts.size()*sizeof(vect_vnctcc_t::value_type));
				auto vret(vbo_cache.alloc(verts_sz, 0));
				vao_mgr.vbo = vret.vbo;
				check_bind_vbo(vao_mgr.vbo);

				if (vret.size == 0) { // newly created
					vert_vbo_sz = verts_sz;
					upload_vbo_data(quad_verts.data(), verts_sz);
				}
				else { // existing
					vert_vbo_sz = vret.size;
					assert(verts_sz <= vert_vbo_sz);
					upload_vector_to_vbo(quad_verts);
				}
				bind_vbo(0);
			}
			clear_cont(quad_verts); // no longer needed
		}
		void register_tile_id(unsigned tid) {
			if (tid+1 == pos_by_tile.size()) return; // already saw this tile
			assert(tid >= pos_by_tile.size()); // tid must be strictly increasing
			pos_by_tile.resize(tid+1, vert_ix_pair(num_quad_verts(), num_tri_verts())); // push start of new range back onto all previous tile slots
		}
		void finalize(unsigned num_tiles) {
			if (pos_by_tile.empty()) return; // nothing to do
			register_tile_id(num_tiles); // add terminator
			remove_excess_cap(pos_by_tile);
		}
		void clear_verts() {quad_verts.clear(); tri_verts.clear(); pos_by_tile.clear();}
		
		void clear_vbos() {
			vbo_cache.free(vao_mgr.vbo, vert_vbo_sz, 0);
			vao_mgr.clear_vaos(); // Note: VAOs not reused
			vert_vbo_sz = 0;
		}
		void clear() {clear_vbos(); clear_verts();}
		bool empty() const {return (quad_verts.empty() && tri_verts.empty());}
		bool has_drawn() const {return !pos_by_tile.empty();}
		unsigned num_quad_verts() const {return quad_verts.size();}
		unsigned num_tri_verts () const {return tri_verts .size();}
		unsigned num_verts() const {return (num_quad_verts() + num_tri_verts());}
		unsigned num_tris () const {return (num_quad_verts()/2 + num_tri_verts()/3);} // Note: 1 quad = 4 verts = 2 triangles
		unsigned start_quad_vert() const {return start_num_verts[0];}
		unsigned start_tri_vert () const {return start_num_verts[1];}
	}; // end draw_block_t
	vector<draw_block_t> to_draw; // one per texture, assumes tids are dense

public:
	vect_vnctcc_t &get_verts(tid_nm_pair_t const &tex, bool quads_or_tris=0) { // default is quads
		unsigned const ix(get_to_draw_ix(tex));
		if (ix >= to_draw.size()) {to_draw.resize(ix+1);}
		draw_block_t &block(to_draw[ix]);
		block.register_tile_id(cur_tile_id);
		if (block.empty()) {block.tex = tex;} // copy material first time
		else {
			assert(block.tex.tid == tex.tid);
			int const bnm(block.tex.get_nm_tid()), tnm(tex.get_nm_tid());

			if (bnm != tnm) { // else normal maps must agree
				if (bnm == FLAT_NMAP_TEX) {block.tex.nm_tid = tnm;} // assume this normal map is correct and assign it to the block
				else if (tnm != FLAT_NMAP_TEX) { // allow if if block has normal map but tex does not - block will override the texture
					std::cerr << "mismatched normal map for texture ID " << block.tex.tid << " in slot " << ix << ": " << bnm << " vs. " << tnm << endl;
					assert(0);
				}
			}
			// if new texture has specular and block does not, copy specular parameters from new texture; this is needed for house wood floors
			if (tex.shininess && !block.tex.shininess) {block.tex.spec_color = tex.spec_color; block.tex.shininess = tex.shininess;}
		}
		return (quads_or_tris ? block.tri_verts : block.quad_verts);
	}
private:
	static void setup_ao_color(colorRGBA const &color, float bcz1, float ao_bcz2, float z1, float z2, color_wrapper cw[2], vert_norm_comp_tc_color &vert, bool no_ao) {
		if (!no_ao && global_building_params.ao_factor > 0.0) {
			min_eq(z1, ao_bcz2); min_eq(z2, ao_bcz2); // clamp zvals to AO zmax
			float const dz_mult(global_building_params.ao_factor/(ao_bcz2 - bcz1));
			UNROLL_2X(cw[i_].set_c4(color*((1.0f - global_building_params.ao_factor) + dz_mult*((i_ ? z2 : z1) - bcz1)));)
		} else {vert.set_c4(color);} // color is shared across all verts
	}
	vector<vector3d> normals; // reused across add_cylinder() calls
	vector<vert_norm_tc> sphere_verts; // reused
	point cur_camera_pos;
	bool is_city;

	struct wall_seg_t {
		float dlo, dhi, ilo, ihi;
		wall_seg_t() : dlo(0.0), dhi(1.0), ilo(0.0), ihi(1.0) {}

		wall_seg_t(float dlo_, float dhi_, float ilo_, float ihi_) : dlo(dlo_), dhi(dhi_), ilo(ilo_), ihi(ihi_) {
			assert(dlo <= dhi && ilo <= ihi && dlo >= 0.0f && dhi <= 1.0f && ilo >= 0.0f && ihi <= 1.0f); // should be (dlo < dhi && ilo < ihi), but can fail due to FP error
		}
		bool is_normalized() const {return (dlo < dhi && ilo < ihi);}
	};
	vector<wall_seg_t> segs;
	vect_cube_t faces;

public:
	unsigned cur_tile_id;
	vect_cube_t temp_cubes, temp_cubes2;
	vector<float> temp_wall_edges;
	vect_tquad_with_ix_t temp_tquads;

	building_draw_t(bool is_city_=0) : cur_camera_pos(zero_vector), is_city(is_city_), cur_tile_id(0) {}
	void init_draw_frame() {cur_camera_pos = get_camera_pos();} // capture camera pos during non-shadow pass to use for shadow pass
	bool empty() const {return to_draw.empty();}
	void reserve_verts(tid_nm_pair_t const &tex, size_t num, bool quads_or_tris=0) {get_verts(tex, quads_or_tris).reserve(num);}
	unsigned get_to_draw_ix(tid_nm_pair_t const &tex) const {return tid_mapper.get_slot_ix(tex.tid);}
	int      get_to_draw_ix_if_exists(tid_nm_pair_t const &tex) const {return tid_mapper.get_slot_ix_if_exists(tex.tid);} // returns -1 if not found
	unsigned get_num_verts (tid_nm_pair_t const &tex, bool quads_or_tris=0) {return get_verts(tex, quads_or_tris).size    ();}
	unsigned get_cap_verts (tid_nm_pair_t const &tex, bool quads_or_tris=0) {return get_verts(tex, quads_or_tris).capacity();}
	vect_vnctcc_t &get_text_verts() {return get_verts(tid_nm_pair_t(FONT_TEXTURE_ID, 1.0, false, true));} // quads, unshadowed, transparent

	void print_stats() const {
		for (draw_block_t const &b : to_draw) {
			if (!b.empty()) {cout << "S=" << b.num_verts() << " " << get_texture_by_id(b.tex.tid).name << endl;}
		}
	}
	void get_all_mat_verts(vect_vnctcc_t &verts, bool triangles) const {
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {vector_add_to((triangles ? i->tri_verts : i->quad_verts), verts);}
	}
	void begin_draw_range_capture() {
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->record_num_verts();}
	}
	void end_draw_range_capture(draw_range_t &r) const { // capture quads added since begin_draw_range_capture() call across to_draw
		for (unsigned i = 0, rix = 0; i < to_draw.size(); ++i) { // quads
			unsigned const start(to_draw[i].start_quad_vert()), end(to_draw[i].num_quad_verts());
			if (start == end) continue; // empty, skip
			assert(start < end);
			assert(rix < MAX_DRAW_BLOCKS); // make sure we have enough slots for this entry
			r.vrq[rix++] = vertex_range_t(start, end, i);
		} // for i
		for (unsigned i = 0, rix = 0; i < to_draw.size(); ++i) { // triangles
			unsigned const start(to_draw[i].start_tri_vert()), end(to_draw[i].num_tri_verts());
			if (start == end) continue; // empty, skip
			assert(start < end);
			assert(rix < MAX_DRAW_BLOCKS); // make sure we have enough slots for this entry
			r.vrt[rix++] = vertex_range_t(start, end, i);
		} // for i
	}
	void toggle_transparent_windows_mode() {
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->tex.toggle_transparent_windows_mode();}
	}
	void set_no_shadows_for_tex(tid_nm_pair_t const &tex) { // must call get_verts() on this tex first
		if (to_draw.empty()) return; // no geometry; can get here with calls of tex.tid=-1 for empty building tiles
		int const ix(get_to_draw_ix_if_exists(tex));
		if (ix < 0) return; // tex doesn't exist - ignore it
		assert((unsigned)ix < to_draw.size());
		to_draw[ix].no_shadows = 1;
	}
	void add_cylinder(building_t const &bg, cube_t const &cube, tid_nm_pair_t const &tex,
		colorRGBA const &color, unsigned dim_mask, bool skip_bottom, bool skip_top, bool no_ao, bool clip_windows)
	{
		cube_t part(cube); // assume full part
		
		// skip this step for complex floorplan buildings because it doesn't work since multiple parts can contain the cube, and it's not needed anyway
		if (dim_mask == 4 && !bg.has_complex_floorplan) { // only for floors and ceilings
			// find the part containing this cube to determine if we need to clip the cylinder; needed for cutting holes in ceilings and floors for building interiors
			cube_t test_cube(cube);
			test_cube.expand_by(-0.1*cube.dz()); // shrink slightly to avoid failing due to FP error in clipping
			part = bg.get_part_containing_cube(test_cube);
			assert(part.is_strictly_normalized()); // must be found
		}
		//float const rscale(0.5*((bg.num_sides <= 8) ? SQRT2 : 1.0)); // larger for triangles/cubes/hexagons/octagons (to ensure overlap/connectivity), smaller for cylinders
		float const rscale(0.5); // use shape contained in bcube so that bcube tests are correct, since we're not creating L/T/U shapes for this case
		// get the bounds from the part, but clip to the cube
		point const pos(part.xc(), part.yc(), cube.z1());
		float const height(cube.dz()), rx(rscale*part.dx()), ry(rscale*part.dy());
		unsigned ndiv(bg.num_sides); // Note: no LOD
		assert(ndiv >= 3);
		bool const smooth_normals(ndiv >= 16); // cylinder vs. N-gon
		float const bcz1(bg.bcube.z1()), z_top(pos.z + height); // adjust for local vs. global space change
		float const ts_factor((dim_mask == 4) ? 1.0 : 2.0), tscale_x(ts_factor*tex.tscale_x), tscale_y(ts_factor*tex.tscale_y);
		bool const apply_ao(!no_ao && global_building_params.ao_factor > 0.0);
		vert_norm_comp_tc_color vert;
		color_wrapper cw[2];
		setup_ao_color(color, bcz1, bg.ao_bcz2, pos.z, z_top, cw, vert, no_ao);
		float tex_pos[2] = {0.0, 1.0};
		building_draw_utils::calc_normals(bg, normals, ndiv);
		UNROLL_2X(tex_pos[i_] = ((i_ ? z_top : pos.z) - bcz1););

		if (dim_mask & 3) { // draw sides
			auto &verts(get_verts(tex)); // Note: cubes are drawn with quads, so we want to emit quads here
			float tot_perim(0.0), cur_perim[2] = {0.0, 0.0};
			for (unsigned S = 0; S < ndiv; ++S) {tot_perim += p2p_dist(normals[S], normals[(S+1)%ndiv]);}
			float const tscale_mult(TWO_PI*sqrt((rx*rx + ry*ry)/2.0f)/tot_perim);
				
			for (unsigned S = 0; S < ndiv; ++S) { // generate vertex data quads
				vector3d const &n1(normals[S]), &n2(normals[(S+1)%ndiv]);
				cur_perim[0]  = cur_perim[1];
				cur_perim[1] += p2p_dist(n1, n2);
				vector3d normal(n1 + n2); normal.x *= ry; normal.y *= rx; // average the two vertex normals for the flat face normal
				if (bg.is_rotated()) {bg.do_xy_rotate_normal(normal);}
				bool const cur_smooth_normals(smooth_normals && (bg.flat_side_amt == 0.0 || S+1 != ndiv)); // flat side of cylindrical building is not smooth
				if (!cur_smooth_normals) {vert.set_norm(normal.get_norm());}

				for (unsigned d = 0; d < 2; ++d) {
					vector3d const &n(d ? n2 : n1);
					vert.t[0] = tscale_x*cur_perim[d]*tscale_mult + tex.txoff; // Note: could try harder to ensure an integer multiple to fix seams, but not a problem in practice
					
					if (cur_smooth_normals) {
						vector3d normal(n); normal.x *= ry; normal.y *= rx; // scale normal by radius (swapped)
						if (bg.is_rotated()) {bg.do_xy_rotate_normal(normal);}
						vert.set_norm(normal.get_norm());
					}
					vert.v.assign((pos.x + rx*n.x), (pos.y + ry*n.y), 0.0);
					if (bg.is_rotated()) {bg.do_xy_rotate(pos, vert.v);}

					for (unsigned e = 0; e < 2; ++e) { // top/bottom
						vert.v.z = ((d^e) ? z_top : pos.z);
						vert.t[1] = tscale_y*tex_pos[d^e] + tex.tyoff;
						if (apply_ao) {vert.copy_color(cw[d^e]);}
						verts.push_back(vert);
					}
					if (clip_windows) {clip_low_high_tc(verts[verts.size()-1].t[1], verts[verts.size()-2].t[1]);} // is this necessary?
				} // for d
			} // for S
		} // end draw sides
		if (dim_mask & 4) { // draw end(s) / ceiling/floor/roof
			auto &tri_verts(get_verts(tex, 1));
			// convert normals to vertices
			vector<point> &verts(normals);
			for (point &v : verts) {v.assign((pos.x + rx*v.x), (pos.y + ry*v.y), pos.z);}
			float tsx(tscale_x*(part.dx()/bg.bcube.dx())/rx), tsy(tscale_y*(part.dy()/bg.bcube.dy())/ry);

			if (part != cube) {
				// clip verts to cube and update ndiv; this isn't the cleanest or most efficient solution, but it's the simplest that uses existing math functions
				vector<point> verts2;
				verts2.reserve(verts.size());
				clip_polygon_xy(verts, cube, verts2);
				verts.swap(verts2);
				ndiv = verts.size();
			}
			for (unsigned d = 0; d < 2; ++d) { // bottom, top
				if (ndiv < 3) continue; // shouldn't happen, but maybe can due to FP error
				if (d ? skip_top : skip_bottom) continue;
				if (is_city && pos.z == bcz1 && d == 0) continue; // skip bottom
				vert.set_ortho_norm(2, d); // +/- z
				if (apply_ao) {vert.copy_color(cw[d]);}
				float const zval(pos.z + (d ? height : 0.0));
				// first vertex is shared across all triangles
				unsigned const start(tri_verts.size());

				for (unsigned S = 0; S < ndiv-2; ++S) { // generate vertex data triangles
					for (unsigned n = 0; n < 3; ++n) {
						if (S > 0 && n == 0) {tri_verts.push_back(tri_verts[start]); continue;} // reuse shared vertex
						if (S > 0 && n == 1) {tri_verts.push_back(tri_verts[tri_verts.size()-2]); continue;} // reuse prev vertex
						vert.v    = verts[S + n];
						vert.v.z  = zval;
						vert.t[0] = tsx*(vert.v.x - pos.x); vert.t[1] = tsy*(vert.v.y - pos.y);
						if (bg.is_rotated()) {bg.do_xy_rotate(pos, vert.v);}
						tri_verts.push_back(vert);
					} // for e
				} // for S
				if (d == 1) {std::reverse(tri_verts.begin()+start, tri_verts.end());} // winding order is wrong, but it's easier to reverse it than change all of the indexing logic
			} // for d
		} // end draw end(s)
	}

	void add_tquad(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex, colorRGBA const &color,
		bool invert_tc_x=0, bool exclude_frame=0, bool no_tc=0, bool swap_tc_xy=0)
	{
		add_tquad_to_verts(bg, tquad, bcube, tex, color, get_verts(tex, (tquad.npts == 3)), invert_tc_x, exclude_frame, no_tc, 0, swap_tc_xy); // 0=quads, 1=tris
	}

	static void set_rotated_normal(vert_norm_comp_tc_color &vert, building_t const &b, vector3d &norm, unsigned n, bool dir) {
		norm.z = 0.0; // likely doesn't need to be set, but okay to set to 0
		if (n == 0) {norm.x =  b.rot_cos; norm.y = b.rot_sin;} // X
		else        {norm.x = -b.rot_sin; norm.y = b.rot_cos;} // Y
		vert.set_norm(dir ? norm : -norm);
	}

	// clip_windows: 0=no clip, 1=clip for building, 2=clip for house
	// dim_mask bits: enable dims: 1=x, 2=y, 4=z | disable cube faces: 8=x1, 16=x2, 32=y1, 64=y2, 128=z1, 256=z2
	void add_section(building_t const &bg, bool clip_to_other_parts, cube_t const &cube, tid_nm_pair_t const &tex,
		colorRGBA const &color, unsigned dim_mask, bool skip_bottom, bool skip_top, bool no_ao, int clip_windows,
		float door_ztop=0.0, unsigned door_sides=0, float offset_scale=1.0, bool invert_normals=0, cube_t const *const clamp_cube=nullptr)
	{
		assert(bg.num_sides >= 3); // must be nonzero volume

		if ((clip_to_other_parts || dim_mask == 4) && !bg.is_cube()) {
			// not a cube, use cylinder; applies to exterior walls (clip_to_other_parts=1) and ceilings/floors (dim_mask == 4)
			//assert(door_ztop == 0.0); // not supported / ignored for testing purposes
			add_cylinder(bg, cube, tex, color, dim_mask, skip_bottom, skip_top, no_ao, clip_windows);
			return;
		}
		// else draw as a cube (optimized flow)
		bool const is_rotated(bg.is_rotated());
		point const center(!is_rotated ? all_zeros : bg.bcube.get_cube_center()); // rotate about bounding cube / building center
		vector3d const sz(cube.get_size()), llc(cube.get_llc()); // move origin from center to min corner
		auto &verts(get_verts(tex, bg.is_pointed)); // bg.is_pointed ? tris : quads
		vert_norm_comp_tc_color vert;
		if (bg.is_pointed) {dim_mask &= 3;} // mask off z-dim since pointed objects (antenna) have no horizontal surfaces
		float const tscale[2] = {tex.get_drawn_tscale_x(), tex.get_drawn_tscale_y()};
		bool const apply_ao(!no_ao && global_building_params.ao_factor > 0.0);
		color_wrapper cw[2];
		setup_ao_color(color, bg.bcube.z1(), bg.ao_bcz2, cube.z1(), cube.z2(), cw, vert, no_ao);
		vector3d tex_vert_off(((world_mode == WMODE_INF_TERRAIN) ? zero_vector : vector3d(xoff2*DX_VAL, yoff2*DY_VAL, 0.0)));
		point const bcube_llc(bg.bcube.get_llc());
		// don't adjust X/Y pos for windows, because other code needs to know where windows are placed; see tc_xlate code below
		if (clip_windows) {tex_vert_off.z -= bcube_llc.z;}
		else {tex_vert_off -= bcube_llc;} // normalize to building LLC to keep tex coords small
		if (is_city && cube.z1() == bg.bcube.z1()) {skip_bottom = 1;} // skip bottoms of first floor parts drawn in cities
		float const window_vspacing(bg.get_window_vspace()), window_h_border(0.75*bg.get_window_h_border()), offset_val(offset_scale*bg.get_door_shift_dist());
		vector3d norm; // used for rotated buildings

		for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
			unsigned const n((i+2)%3), d((i+1)%3), st(i&1); // n = dim of normal, i/d = other dims
			if (!(dim_mask & (1<<n))) continue; // check for enabled dims
			bool const do_xy_clip(clip_windows && n < 2 && !is_rotated), clip_d(d != 2); // clip the non-z dim;

			for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
				if (skip_bottom && n == 2 && j == 0) continue; // skip bottom side
				if (skip_top    && n == 2 && j == 1) continue; // skip top    side
				unsigned const dir(bool(j) ^ invert_normals);
				if (dim_mask & (1<<(2*n+dir+3)))     continue; // check for disabled faces
				if (clamp_cube != nullptr && (cube.d[n][dir] < clamp_cube->d[n][0] || cube.d[n][dir] > clamp_cube->d[n][1])) continue; // outside clamp cube, drop this face
				if (n < 2 && is_rotated) {set_rotated_normal(vert, bg, norm, n, j);} // XY only
				else {vert.n[i] = 0; vert.n[d] = 0; vert.n[n] = (j ? 127 : -128);} // -1.0 or 1.0
				point pt; // parameteric position within cube in [vec3(0), vec3(1)]
				pt[n] = dir; // our cube face, in direction of normal

				if (bg.is_pointed) { // antenna triangle; parts clipping doesn't apply to this case since there are no opposing cube faces
					unsigned const ix(verts.size()); // first vertex of this triangle
					assert(door_ztop == 0.0); // not supported
					pt[!n] = j^n^1; pt.z = 0;
					EMIT_VERTEX(); // bottom low
					pt[!n] = j^n;
					EMIT_VERTEX(); // bottom high
					pt[!n] = 0.5; pt[n] = 0.5; pt.z = 1;
					EMIT_VERTEX(); // top
					vector3d normal;
					get_normal(verts[ix+0].v, verts[ix+1].v, verts[ix+2].v, normal, 1); // update with correct normal
					vert.set_norm(normal);
					UNROLL_3X(verts[ix+i_].set_norm(vert);)
					continue; // no windows/clipping
				}
				segs.clear();

				if (bg.has_interior() && clip_to_other_parts && bg.real_num_parts > 1 && n != 2) {
					// clip walls XY to remove intersections; this applies to both walls and windows
					cube_t face(cube);
					face.d[n][!j] = face.d[n][j]; // shrink to zero thickness face
					faces.clear();
					faces.push_back(face);
					float const sz_d_inv(1.0/sz[d]), sz_i_inv(1.0/sz[i]);
					float num_windows(0.0);

					if (do_xy_clip) { // side of non-rotated building
						unsigned const cdim(clip_d ? d : i); // clip dim, horizontal (X or Y)
						float tlo(tscale[0]*(face.d[cdim][0] + tex_vert_off[cdim]) + tex.txoff), thi(tscale[0]*(face.d[cdim][1] + tex_vert_off[cdim]) + tex.txoff); // TCx
						clip_low_high_tc(tlo, thi); // TCx of unclipped wall
						num_windows = (thi - tlo); // for unclipped wall
						assert(num_windows >= 0.0);
						if (num_windows < 0.5) continue; // no space for a window before clip
					}
					for (auto p = bg.parts.begin(); p != bg.get_real_parts_end(); ++p) {
						if (p->contains_cube(cube)) continue; // skip ourself (including door part)
						if (cube.d[d][1] <= p->d[d][0] || cube.d[d][0] >= p->d[d][1] || cube.d[i][1] <= p->d[i][0] || cube.d[i][0] >= p->d[i][1]) continue; // check for overlap in the two quad dims
						subtract_cube_from_cubes(*p, faces, nullptr, 1, 1); // no holes, clip_in_z=1, include_adj=1
					}
					for (unsigned f = 0; f < faces.size(); ++f) { // convert from cubes to parametric coordinates in [0.0, 1.0] range
						cube_t const &F(faces[f]);
						if (F.get_sz_dim(d) == 0.0 || F.get_sz_dim(i) == 0.0) continue; // adjacent part zero area strip, skip
						// don't copy/enable the top segment for house windows because houses always have a sloped roof section on top that will block the windows
						if (clip_windows == 2 && F.z1() > cube.z1()) continue;
						wall_seg_t seg((F.d[d][0] - llc[d])*sz_d_inv, (F.d[d][1] - llc[d])*sz_d_inv, (F.d[i][0] - llc[i])*sz_i_inv, (F.d[i][1] - llc[i])*sz_i_inv);

						if (do_xy_clip) { // clip the horizontal vertices of the quad that was just added to make the windows line up with the other side
							float &lo(clip_d ? seg.dlo : seg.ilo), &hi(clip_d ? seg.dhi : seg.ihi);
							if (lo > 0.01) {lo = ceil (lo*num_windows - window_h_border)/num_windows;} // if clipped on lo edge, round up   to nearest whole window
							if (hi < 0.99) {hi = floor(hi*num_windows + window_h_border)/num_windows;} // if clipped on hi edge, round down to nearest whole window
							if (hi - lo < 0.01f) continue; // no space for a window after clip (optimization)
						}
						segs.emplace_back(seg);
					} // for f
				}
				else {
					segs.push_back(wall_seg_t()); // single seg
				}
				if (clip_windows && bg.has_chimney == 2) { // remove windows blocked by the chimney; there can be at most one segment
					cube_t block(bg.get_chimney());
					block.union_with_cube(bg.get_fireplace()); // merge with fireplace

					if (fabs(cube.d[n][j] - block.d[n][!j]) < 0.01*window_vspacing) { // chimney is adjacent to this side, or close enough
						// remove vertical columns of whole windows that overlap the chimney or fireplace
						unsigned const cdim(clip_d ? d : i); // clip dimension
						float const ts(tscale[bool(st) ^ clip_d ^ 1]), toff(cdim ? tex.tyoff : tex.txoff); // texture scale and offset for clip dim
						float t_lo(ts*(cube.d[cdim][0] + tex_vert_off[cdim]) + toff), t_hi(ts*(cube.d[cdim][1] + tex_vert_off[cdim]) + toff);
						clip_low_high_tc(t_lo, t_hi);
						float const num_windows(round_fp(t_hi - t_lo)); // window count for this wall; should be an exact integer

						if (num_windows > 0) { // we have at least one window
							float const edge_buffer(0.9*bg.get_window_h_border()/num_windows); // slightly smaller than the space to the sides of a window, in parametric space
							float c_lo((block.d[cdim][0] - llc[cdim])/sz[cdim]), c_hi((block.d[cdim][1] - llc[cdim])/sz[cdim]); // parametric bounds of chimney
							c_lo += edge_buffer; c_hi -= edge_buffer; // space to the left and right of the window is allowed to overlap the chimney/fireplace
							c_lo  = floor(num_windows*c_lo)/num_windows; // round down to an exact window boundary
							c_hi  = ceil (num_windows*c_hi)/num_windows; // round up   to an exact window boundary

							for (auto s = segs.begin(); s != segs.end(); ++s) {
								float &lo(clip_d ? s->dlo : s->ilo), &hi(clip_d ? s->dhi : s->ihi);
								if (lo > c_hi || hi < c_lo) continue; // doesn't overlap the chimney (side clipped to door may partially overlap the chimney)
								wall_seg_t s2(*s);
								hi = c_lo; // *s becomes lo seg
								(clip_d ? s2.dlo : s2.ilo) = c_hi; // s2 becomes hi seg
								if (s2.is_normalized()) {segs.push_back(s2);} // Note: s->is_normalized() check is done in the segs iteration below
								break; // s is invalidated, have to break; there shouldn't be any other segs spanning chimney anyway
							} // for s
						}
					}
				}
				for (auto s = segs.begin(); s != segs.end(); ++s) {
					wall_seg_t const &seg(*s);
					if (!seg.is_normalized()) continue; // zero area or denormalized - skip (can get here due to chimney clipping)
					unsigned const ix(verts.size()); // first vertex of this quad
					pt[d] = seg.dlo;
					pt[i] = (j ? seg.ilo : seg.ihi); // need to orient the vertices differently for each side
					//if (bg.roof_recess > 0.0 && n == 2 && j == 1) {pt.z -= bg.roof_recess*cube.dz();}
					EMIT_VERTEX(); // 0 !j
					pt[i] = (j ? seg.ihi : seg.ilo);
					EMIT_VERTEX(); // 0 j
					pt[d] = seg.dhi;
					EMIT_VERTEX(); // 1 j
					pt[i] = (j ? seg.ilo : seg.ihi);
					EMIT_VERTEX(); // 1 !j
					float const offset((j ? 1.0 : -1.0)*offset_val);

					if (((i == 2) ? seg.ilo : seg.dlo) == 0.0 && (door_sides & (1 << (2*n + j)))) {
						// clip zval to exclude door z-range (except for top segment); doesn't work with walkway doors
						for (unsigned k = ix; k < ix+4; ++k) {
							auto &v(verts[k]);
							float const delta(door_ztop - v.v.z);
							if (v.v.z < door_ztop) {v.v.z = door_ztop;} // make all windows start above the door
							// move slightly away from the building wall to avoid z-fighting (vertex is different from building and won't have same depth)
							if (is_rotated) {v.v += offset*norm;} else {v.v[n] += offset;}
							if (delta > 0.0) {v.t[1] += tscale[1]*delta;} // recalculate tex coord
						}
					}
					else if (clip_windows) { // move slightly away from the building wall to avoid z-fighting
						for (unsigned k = ix; k < ix+4; ++k) {
							if (is_rotated) {verts[k].v += offset*norm;} else {verts[k].v[n] += offset;}
						}
					}
					if (clip_windows && n < 2) { // clip the texture coordinates of the quad that was just added (side of building)
						clip_low_high_tc(verts[ix+0].t[!st], verts[ix+1].t[!st]);
						clip_low_high_tc(verts[ix+2].t[!st], verts[ix+3].t[!st]);
						clip_low_high_tc(verts[ix+0].t[ st], verts[ix+3].t[ st]);
						clip_low_high_tc(verts[ix+1].t[ st], verts[ix+2].t[ st]);
						// shift texture coords tp a local reference point near the building origin to prevent noise in shader texture sampling for large numbers
						float const tc_xlate[2] = {floor(verts[ix].t[0]), floor(verts[ix].t[1])}; // use first point for local origin

						for (unsigned n = 0; n < 4; ++n) {
							for (unsigned d = 0; d < 2; ++d) {verts[ix+n].t[d] -= tc_xlate[d];}
						}
						// Note: if we're drawing windows, and either of the texture coords have zero ranges, we can drop this quad; but this is uncommon and maybe not worth the trouble
					}
					if (clamp_cube != nullptr && *clamp_cube != cube && n < 2 && !is_rotated) { // x/y dims only; can't apply to rotated building
						unsigned const dim((i == 2) ? d : i); // x/y
						unsigned const sec_vix(ix + (st ? 1 : 3)); // opposite vertex in this dim
						float const delta_tc(verts[sec_vix].t[0]   - verts[ix].t[0]  );
						float const delta_v (verts[sec_vix].v[dim] - verts[ix].v[dim]);
						assert(delta_v != 0.0);
						float const fin_tscale(delta_tc/delta_v); // tscale[0] post-clip

						for (unsigned k = ix; k < ix+4; ++k) {
							auto &v(verts[k]);
							float &val(v.v[dim]);
							float const dlo(clamp_cube->d[dim][0] - val), dhi(val - clamp_cube->d[dim][1]); // dists outside clamp_cube
							if (dlo > 0.0) {val += dlo; v.t[0] += fin_tscale*dlo;} // shift pos
							if (dhi > 0.0) {val -= dhi; v.t[0] -= fin_tscale*dhi;} // shift neg
						} // for k
					}
				} // for seg s
			} // for j
		} // for i
	}

	void add_cube(building_t const &bg, cube_t const &cube, tid_nm_pair_t const &tex, colorRGBA const &color,
		bool swap_txy=0, unsigned dim_mask=7, bool skip_bottom=0, bool skip_top=0, bool ws_texture=0)
	{
		assert(dim_mask != 0); // must draw at least some face
		auto &verts(get_verts(tex));
		vector3d const sz(cube.get_size()), llc(cube.get_llc()); // move origin from center to min corner
		vert_norm_comp_tc_color vert;
		vert.set_c4(color);
		float const tscale[2] = {tex.tscale_x, tex.tscale_y};
		unsigned const verts_start(verts.size());
		bool const is_rotated(bg.is_rotated());
		vector3d norm; // used for rotated buildings

		for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
			unsigned const n((i+2)%3), d((i+1)%3), st(bool(i&1) ^ swap_txy); // n = dim of normal, i/d = other dims
			if (!(dim_mask & (1<<n))) continue; // check for enabled dims

			for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
				if (skip_bottom && n == 2 && j == 0) continue; // skip bottom side
				if (skip_top    && n == 2 && j == 1) continue; // skip top    side
				if (n < 2 && is_rotated) {set_rotated_normal(vert, bg, norm, n, j);} // XY only
				else {vert.n[i] = 0; vert.n[d] = 0; vert.n[n] = (j ? 127 : -128);} // -1.0 or 1.0
				point pt; // parameteric position within cube in [vec3(0), vec3(1)]
				pt[n] = j; // our cube face, in direction of normal
				pt[d] = 0.0;
				pt[i] = !j; // need to orient the vertices differently for each side
				EMIT_VERTEX_SIMPLE(); // 0 !j
				pt[i] = j;
				EMIT_VERTEX_SIMPLE(); // 0 j
				pt[d] = 1.0;
				EMIT_VERTEX_SIMPLE(); // 1 j
				pt[i] = !j;
				EMIT_VERTEX_SIMPLE(); // 1 !j
			} // for j
		} // for i
		if (is_rotated) {
			point const center(bg.bcube.get_cube_center());
			for (auto i = (verts.begin() + verts_start); i != verts.end(); ++i) {bg.do_xy_rotate(center, i->v);}
		}
	}

	void add_fence(building_t const &bg, cube_t const &fence, tid_nm_pair_t const &tex, colorRGBA const &color, bool mult_sections) {
		bool const dim(fence.dy() < fence.dx()); // smaller/separating dim
		float const length(fence.get_sz_dim(!dim)), height(fence.dz());
		float const post_width(fence.get_sz_dim(dim)), post_hwidth(0.5*post_width), beam_hwidth(0.5*post_hwidth), beam_hheight(1.0*post_hwidth);
		unsigned const num_sections(ceil(0.3*length/height)), num_posts(num_sections + 1), num_beams(2); // could also use 3 beams
		float const post_spacing((length - post_width)/num_sections), beam_spacing(height/(num_beams+0.5f));
		cube_t post(fence); // copy dim and Z values
		cube_t beam(fence); // copy dim and !dim values
		beam.expand_in_dim(!dim, -post_width); // remove overlap with end posts at both ends
		beam.expand_in_dim( dim, (beam_hwidth - post_hwidth));
		unsigned skip_ix(num_posts); // start at an invalid value

		if (mult_sections && dim == 0) { // skip end post on the corner in dim=0 because it's duplicated with the corner post in the other dim
			float const dmin(post_width + ((bg.has_chimney == 2) ? bg.get_chimney().get_sz_dim(!dim) : 0.0)); // include chimney width, can increase the house bcube beyond the fence
			if      (fabs(bg.bcube.d[!dim][1] - fence.d[!dim][1]) < dmin) {beam.d[!dim][1] += post_hwidth; skip_ix = num_posts-1;} // skip last post
			else if (fabs(bg.bcube.d[!dim][0] - fence.d[!dim][0]) < dmin) {beam.d[!dim][0] -= post_hwidth; skip_ix = 0;} // skip first post
		}
		for (unsigned i = 0; i < num_posts; ++i) { // add posts
			set_wall_width(post, (fence.d[!dim][0] + post_hwidth + i*post_spacing), post_hwidth, !dim);
			if (i != skip_ix) {add_cube(bg, post, tex, color, 0, 7, 1, 0, 1);} // skip bottom, ws_texture=1
		}
		for (unsigned i = 0; i < num_beams; ++i) { // add beams
			set_wall_width(beam, (fence.z1() + (i+1)*beam_spacing), beam_hheight, 2); // set beam zvals
			add_cube(bg, beam, tex, color, 0, (4U + (1U<<unsigned(dim))), 0, 0, 1); // skip !dim sides, ws_texture=1
		}
	}

	void add_roof_dome(point const &pos, float rx, float ry, tid_nm_pair_t const &tex, colorRGBA const &color, bool onion) { // Note: no rotation needed, no bg argument
		auto &verts(get_verts(tex));
		color_wrapper cw(color);
		unsigned const ndiv(N_SPHERE_DIV);
		float const ravg(0.5f*(rx + ry)), t_end(onion ? 1.0 : 0.5);
		point center(pos);
		if (onion) {center.z += 0.5*ravg; rx *= 1.2; ry *= 1.2;} // move up slightly and increase radius
		float const arx(rx/ravg), ary(ry/ravg);
		sphere_verts.clear();
		sd_sphere_d sd(all_zeros, ravg, ndiv);
		sd.gen_points_norms_static(0.0, 1.0, 0.0, t_end); // top half hemisphere dome
		sd.get_quad_points(sphere_verts, nullptr, 0, 0.0, 1.0, 0.0, t_end); // quads
			
		for (auto i = sphere_verts.begin(); i != sphere_verts.end(); ++i) {
			i->v.y *= ary;
			i->v.x *= arx;
			if (onion && i->v.z > 0.0) {i->v.z += 0.05f*ravg*(1.0f/(1.01f - i->v.z/ravg) - 1.0f);} // form a point at the top
			verts.emplace_back(vert_norm_comp_tc((i->v + center), i->n, i->t[0]*tex.tscale_x, i->t[1]*tex.tscale_y), cw);
		}
	}

	void add_vert_cylinder(point const &center, float z1, float z2, float radius, float tscale_x, float tscale_y, unsigned ndiv, colorRGBA const &color, vect_vnctcc_t &qverts) {
		float const ndiv_inv(1.0/ndiv);
		color_wrapper const cw(color);
		point const ce[2] = {point(center.x, center.y, z1), point(center.x, center.y, z2)};
		vector3d v12;
		vector_point_norm const &vpn(gen_cylinder_data(ce, radius, radius, ndiv, v12));

		for (unsigned i = 0; i < ndiv; ++i) { // similar to gen_cylinder_quads(), but with a color
			for (unsigned j = 0; j < 2; ++j) {
				unsigned const S(i + j), s(S%ndiv);
				norm_comp const normal((vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]).get_norm());
				float const ts(4.0*S*tscale_x*ndiv_inv);
				qverts.emplace_back(vpn.p[(s<<1)+ j], normal, ts, (1.0-j)*tscale_y, cw);
				qverts.emplace_back(vpn.p[(s<<1)+!j], normal, ts, j      *tscale_y, cw);
			}
		} // for i
	}
	void add_water_tower(building_t const &bg, cube_t const &wtc) {
		tid_nm_pair_t const side_tex(building_texture_mgr.get_met_plate_tid(), building_texture_mgr.get_mplate_nm_tid(), 1.0, 1.0);
		tid_nm_pair_t const base_tex(building_texture_mgr.get_met_roof_tid ()); // no normal map
		tid_nm_pair_t const roof_tex(WHITE_TEX); // untextured
		auto &tverts(get_verts(roof_tex, 1)), &qverts(get_verts(side_tex, 0)); // triangle and quad verts
		unsigned const ndiv(N_CYL_SIDES);
		float const radius(0.25*(wtc.dx() + wtc.dy())); // should be equal size in X vs. Y
		float const ndiv_inv(1.0/ndiv), height(wtc.dz());
		float const base_z1(wtc.z1() + 0.5*height - 0.5*radius), cylin_z1(base_z1 + 0.01*height), cylin_z2(wtc.z2() - 0.12*height), cone_z2(wtc.z2());
		color_wrapper const roof_cw(colorRGBA(0.15, 0.12, 0.10, 1.0));
		// draw base
		cube_t base(wtc);
		set_cube_zvals(base, base_z1, cylin_z1);
		add_cube(bg, base, base_tex, GRAY); // draw all sides
		// draw legs
		float const leg_width(0.08*radius);
		cube_t legs(wtc);
		legs.expand_by_xy(-0.5*leg_width); // shrink slightly
		legs.z2() = base_z1;

		for (unsigned y = 0; y < 2; ++y) {
			for (unsigned x = 0; x < 2; ++x) {
				cube_t leg(legs);
				leg.d[0][!x] = leg.d[0][x] + (x ? -1.0 : 1.0)*leg_width;
				leg.d[1][!y] = leg.d[1][y] + (y ? -1.0 : 1.0)*leg_width;
				add_cube(bg, leg, roof_tex, DK_GRAY, 0, 3); // skip top and bottom; untextured like roof
			} // for x
		} // for y
		// draw side quads
		vector3d center(wtc.get_cube_center()), v12;
		if (bg.is_rotated()) {bg.do_xy_rotate(bg.bcube.get_cube_center(), center);}
		add_vert_cylinder(center, cylin_z1, cylin_z2, radius, 2.0, 2.0, ndiv, WHITE, qverts); // tscale=2.0/2.0
		// draw top cone triangles
		point const ce[2] = {point(center.x, center.y, cylin_z2), point(center.x, center.y, cone_z2)};
		vector_point_norm const &vpn(gen_cylinder_data(ce, 1.1*radius, 0.0, ndiv, v12)); // slightly wider at the bottom

		for (unsigned i = 0; i < ndiv; ++i) { // similar to gen_cylinder_quads(), but with a color
			unsigned const ip((i+ndiv-1)%ndiv), in((i+1)%ndiv);
			float const ts(1.0 - i*ndiv_inv);
			tverts.emplace_back(vpn.p[(i <<1)+1],  vpn.n[i],                         0.5,             1.0, roof_cw);
			tverts.emplace_back(vpn.p[(in<<1)+0], (vpn.n[i] + vpn.n[in]).get_norm(), (ts - ndiv_inv), 0.0, roof_cw);
			tverts.emplace_back(vpn.p[(i <<1)+0], (vpn.n[i] + vpn.n[ip]).get_norm(), ts,              0.0, roof_cw);
		} // for i
		// draw pipe through the center going down into the roof
		add_vert_cylinder(center, wtc.z1(), base_z1, 0.1*radius, 1.0, 4.0, ndiv/2, WHITE, qverts); // tscale=1.0/4.0
	}

	unsigned num_verts() const {
		unsigned num(0);
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {num += i->num_verts();}
		return num;
	}
	unsigned num_tris() const {
		unsigned num(0);
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {num += i->num_tris();}
		return num;
	}
	void upload_to_vbos() {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->upload_to_vbos();}}
	void clear_vbos    () {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->clear_vbos();}}
	void clear         () {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->clear();}}
	unsigned get_num_draw_blocks() const {return to_draw.size();}
	void finalize(unsigned num_tiles) {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->finalize(num_tiles);}}
	
	// tex_filt_mode: 0=draw everything, 1=draw exterior walls only, 2=draw roofs and exterior doors, 3=draw roofs only
	void draw(shader_t &s, bool shadow_only, bool direct_draw_no_vbo=0, int tex_filt_mode=0, vertex_range_t const *const exclude=nullptr) {
		tid_nm_pair_dstate_t state(s);

		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {
			if (tex_filt_mode == 1 && !tid_mapper.is_ext_wall_tid(i->tex.tid)) continue; // exterior walls only
			if (tex_filt_mode == 2 && !tid_mapper.is_roof_tid    (i->tex.tid) && !building_texture_mgr.is_door_tid(i->tex.tid)) continue; // roofs and exterior doors only
			if (tex_filt_mode == 3 && !tid_mapper.is_roof_tid    (i->tex.tid)) continue; // roofs only
			bool const use_exclude(exclude && exclude->draw_ix == int(i - to_draw.begin()));
			i->draw_all_geom(state, shadow_only, direct_draw_no_vbo, (use_exclude ? exclude : nullptr));
		}
	}
	void draw_tile(shader_t &s, unsigned tile_id, bool shadow_only=0, float crack_weight=0.0) {
		tid_nm_pair_dstate_t state(s, 0, crack_weight); // no_set_texture=0
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->draw_geom_tile(state, tile_id, shadow_only);}
	}
	void draw_block(shader_t &s, unsigned ix, bool shadow_only, vertex_range_t const *const exclude=nullptr) {
		if (ix >= to_draw.size()) return;
		tid_nm_pair_dstate_t state(s);
		to_draw[ix].draw_all_geom(state, shadow_only, 0, exclude);
	}
	void draw_for_draw_range(shader_t &s, draw_range_t const &draw_range, bool shadow_only=0) {
		tid_nm_pair_dstate_t state(s);
		
		for (unsigned i = 0; i < MAX_DRAW_BLOCKS; ++i) {
			vertex_range_t const &qrange(draw_range.vrq[i]), &trange(draw_range.vrt[i]);
			if (qrange.draw_ix >= 0 && (unsigned)qrange.draw_ix < to_draw.size()) {to_draw[qrange.draw_ix].draw_quad_geom_range(state, qrange, shadow_only);}
			if (trange.draw_ix >= 0 && (unsigned)trange.draw_ix < to_draw.size()) {to_draw[trange.draw_ix].draw_tri_geom_range (state, trange, shadow_only);}
		}
	}
}; // end building_draw_t

/*static*/ vbo_cache_t building_draw_t::vbo_cache;


// *** Drawing ***

int get_building_door_tid(unsigned type) { // exterior doors, and interior doors of limited types
	switch(type) {
	case tquad_with_ix_t::TYPE_HDOOR : return building_texture_mgr.get_hdoor_tid ();
	case tquad_with_ix_t::TYPE_ODOOR : return building_texture_mgr.get_odoor_tid ();
	case tquad_with_ix_t::TYPE_RDOOR : return building_texture_mgr.get_bdoor2_tid();
	case tquad_with_ix_t::TYPE_BDOOR : return building_texture_mgr.get_bdoor_tid ();
	case tquad_with_ix_t::TYPE_BDOOR2: return building_texture_mgr.get_bdoor2_tid();
	case tquad_with_ix_t::TYPE_GDOOR : return building_texture_mgr.get_gdoor_tid ();
	//case tquad_with_ix_t::TYPE_MDOOR : return building_texture_mgr.get_mdoor_tid ();
	default: assert(0);
	}
	return -1; // never gets here
}

tid_nm_pair_t get_concrete_texture(float tscale=16.0) {return tid_nm_pair_t(get_concrete_tid(), tscale);}

void add_driveway_or_porch(building_draw_t &bdraw, building_t const &building, cube_t const &c, colorRGBA const &color, bool skip_bottom) {
	if (c.is_all_zeros()) return;
	tid_nm_pair_t const tex(get_concrete_texture());
	bdraw.add_section(building, 0, c, tex, color, 7, skip_bottom, 0, 1, 0); // all dims, no AO
}

tid_nm_pair_t building_t::get_basement_wall_texture() const { // okay to call if there's no basement
	if ((mat_ix + parts.size()) & 1) { // randomly select one of two textures
		return tid_nm_pair_t(CBLOCK_TEX, get_texture_by_name("normal_maps/cblock_NRM.jpg", 1), 1.0, 1.0);
	}
	else {
		return tid_nm_pair_t(get_texture_by_name("cblock2.jpg"), get_texture_by_name("normal_maps/cblock2_NRM.jpg", 1), 1.0, 1.0);
	}
}

tid_nm_pair_t building_t::get_attic_texture() const {
	if (!has_attic()) return tid_nm_pair_t();
	// plywood 50% of the time, boards 50% of the time
	int const tid((interior->rooms.size() & 1) ? get_plywood_tid() : FENCE_TEX);
	building_mat_t const &mat(get_material());
	return tid_nm_pair_t(tid, get_normal_map_for_bldg_tid(tid), 0.25*mat.house_floor_tex.tscale_x, 0.25*mat.house_floor_tex.tscale_y);
}
colorRGBA building_t::get_floor_tex_and_color(cube_t const &floor_cube, tid_nm_pair_t &tex) const {
	if (has_attic() && floor_cube.z2() > interior->attic_access.z1()) { // attic floor
		tex = get_attic_texture();
		return WHITE; // always white
	}
	bool const in_basement(floor_cube.z2() < ground_floor_z1);
	building_mat_t const &mat(get_material());

	if (is_house) {
		if (has_sec_bldg() && get_sec_bldg().contains_cube(floor_cube)) {tex = get_concrete_texture();} // garage or shed
		else if (in_basement) {tex = mat.basement_floor_tex;} // basement
		else {tex = mat.house_floor_tex;}
	}
	else { // office building
		bool const in_ext_basement(in_basement && !get_basement().contains_cube_xy(floor_cube));

		if ((in_ext_basement && has_mall()) || (has_retail() && floor_cube.z1() == ground_floor_z1)) { // retail or mall
			float const tscale(0.125*mat.floor_tex.tscale_x); // stretch the texture out for large tiles
			tex = tid_nm_pair_t(building_texture_mgr.get_tile_floor_tid(), building_texture_mgr.get_tile_floor_nm_tid(), tscale, tscale);
		}
		else if (in_basement && (has_parking_garage || in_ext_basement)) {tex = get_concrete_texture();} // parking garage or extended basement
		else {tex = mat.floor_tex;} // office block
	}
	return (is_house ? mat.house_floor_color : mat.floor_color);
}
colorRGBA building_t::get_ceil_tex_and_color(cube_t const &ceil_cube, tid_nm_pair_t &tex) const {
	bool const in_basement(ceil_cube.z1() < ground_floor_z1), in_ext_basement(in_basement && !get_basement().contains_cube_xy(ceil_cube));
	building_mat_t const &mat(get_material());
	colorRGBA color;

	if (is_house && in_basement && has_basement_pipes && !in_ext_basement) { // draw wood flooring for basement ceiling (not ext basement)
		tex = mat.house_floor_tex;
		return (is_house ? mat.house_floor_color : mat.floor_color);
	}
	else if (!is_house && in_ext_basement && !has_mall()) { // use concrete for office building extended basements except for malls
		tex = get_concrete_texture();
		return WHITE;
	}
	else if (in_basement && (is_house || (has_parking_garage && !in_ext_basement))) { // use wall texture for basement/parking garage ceilings, not ceiling texture
		tex = mat.wall_tex;
		return WHITE; // basement walls are always white
	}
	// normal ceiling texture
	bool const residential(is_residential()); // apartments and hotels use house ceiling textures and colors
	tex =  (residential ? mat.house_ceil_tex   : mat.ceil_tex  );
	return (residential ? mat.house_ceil_color : mat.ceil_color);
}

void draw_building_ext_door(building_draw_t &bdraw, tquad_with_ix_t const &door, building_t const &b) {
	colorRGBA const &dcolor((door.type == tquad_with_ix_t::TYPE_GDOOR) ? WHITE : b.door_color); // garage doors are always white
	bdraw.add_tquad(b, door, b.bcube, tid_nm_pair_t(get_building_door_tid(door.type), -1, 1.0, 1.0), dcolor);
}
void building_t::get_all_drawn_exterior_verts(building_draw_t &bdraw) { // exterior building parts
	if (!is_valid()) return; // invalid building
	building_mat_t const &mat(get_material());
	bool const need_top_roof(roof_type == ROOF_TYPE_FLAT || roof_type == ROOF_TYPE_DOME || roof_type == ROOF_TYPE_ONION);
	if (detail_color == BLACK) {detail_color = roof_color;} // use roof color if not set
		
	for (auto i = parts.begin(); i != parts.end(); ++i) { // multiple cubes/parts/levels - no AO for houses
		if (is_basement(i)) continue; // don't need to draw the basement exterior walls since they should be underground
		float const chimney_tscale = 1.0; // smaller bricks?
			
		if (has_chimney == 2 && (i+2 == parts.end())) { // fireplace; skip inside face (optimization); needs to be draw as part of the interior instead
			unsigned dim_mask(3); // disable faces: 8=x1, 16=x2, 32=y1, 64=y2
			if      (i->x1() <= bcube.x1()) {dim_mask |= 16;} // Note: may not work on rotated buildings
			else if (i->x2() >= bcube.x2()) {dim_mask |=  8;}
			else if (i->y1() <= bcube.y1()) {dim_mask |= 64;}
			else if (i->y2() >= bcube.y2()) {dim_mask |= 32;}
			tid_nm_pair_t const tp(mat.side_tex.get_scaled_version(chimney_tscale));
			bdraw.add_section(*this, 0, *i, tp, side_color, dim_mask, 0, 0, is_house, 0); // XY exterior walls
			bdraw.add_section(*this, 0, *i, tp, side_color, 4, 1, 0, 1, 0); // draw top of fireplace exterior, even if not wider than the chimney - should it be sloped?
			continue;
		}
		else if (has_chimney && (i+1 == parts.end())) { // chimney
			tid_nm_pair_t const tp(mat.side_tex.get_scaled_version(chimney_tscale));
			auto &verts(bdraw.get_verts(tp));
			unsigned const verts_start(verts.size());
			bdraw.add_section(*this, 0, *i, tp, side_color, 3, 0, 0, is_house, 0); // XY exterior walls
			float const wall_width(0.25*min(i->dx(), i->dy()));
			cube_t hole(*i);
			hole.expand_by_xy(-wall_width);
			cube_t sides[4]; // {-y, +y, -x, +x}
			subtract_cube_xy(*i, hole, sides);
			unsigned const face_masks[4] = {127-64, 127-32, 127-16, 127-8}; // enable XYZ but skip all but {+y, -y, +x, -x} in XY
			for (unsigned n = 0; n < 4; ++n) {bdraw.add_section(*this, 0, sides[n], tp, side_color, face_masks[n], 1, 0, 1, 0);} // skip bottom

			if (has_chimney == 1 && has_attic()) { // clip interior chimney verts to top of roof
				for (auto v = verts.begin()+verts_start; v != verts.end(); ++v) {
					if (v->v.z == i->z2()) continue; // top surface is unchanged

					for (auto const &tq : roof_tquads) {
						if (!is_attic_roof(tq, 1)) continue; // type_roof_only=1
						if (!point_in_polygon_2d(v->v.x, v->v.y, tq.pts, tq.npts)) continue; // check 2D XY point containment
						vector3d const normal(tq.get_norm());
						if (normal.z == 0.0) continue; // skip vertical sides
						float const rdist(dot_product_ptv(normal, v->v, tq.pts[0]));
						if (rdist >= 0) continue; // above the roof
						float const dz(rdist/normal.z);
						v->v.z  -= dz; // move exactly to the roof
						v->t[1] -= 2.0f*tp.tscale_y*dz;
						break;
					} // for tq
				} // for v
			}
			continue;
		}
		bdraw.add_section(*this, 1, *i, mat.side_tex, side_color, 3, 0, 0, is_house, 0); // XY exterior walls
		bool skip_top((!need_top_roof && (is_house || i+1 == parts.end())) || is_basement(i)); // don't add the flat roof for the top part in this case
		skip_top |= (has_retail() && i == parts.begin()); // skip drawing the roof between the retail area and the office above
		// skip the bottom of stacked cubes (not using ground_floor_z1); need to draw the porch roof, so test i->dz()
		bool const is_stacked(i->z1() > bcube.z1() && i->dz() > 0.5f*get_window_vspace()), is_stacked_cube(is_stacked && is_cube());
		if (is_stacked_cube && skip_top) continue; // no top/bottom to draw

		// add roof quads
		if (!is_house && !skip_top && interior && clip_part_ceiling_for_stairs(*i, bdraw.temp_cubes, bdraw.temp_cubes2)) {
			for (auto c = bdraw.temp_cubes.begin(); c != bdraw.temp_cubes.end(); ++c) { // add floors after removing stairwells
				bdraw.add_section(*this, 0, *c, mat.roof_tex, roof_color, 4, 1, 0, is_house, 0); // only top surface
			}
			// still need to draw in shadow pass
			tid_nm_pair_t shadow_mat(SHADOW_ONLY_TEX);
			shadow_mat.shadow_only = 1;
			bdraw.add_section(*this, 1, *i, shadow_mat, WHITE, 4, 1, 0, 1, 0); // only draw the top
			skip_top = 1;
			if (is_stacked_cube) continue; // no top/bottom to draw
		}
		if (is_stacked && !is_cube()) { // handle bottom of stacked parts from non-cube buildings that may overhang
			cube_t test_cube(*i);
			// clip_part_ceiling_for_stairs() expects a ceiling, shift zval to the top floor of the part below
			set_cube_zvals(test_cube, (i->z1() - get_window_vspace() - get_floor_thickness()), i->z1());

			if (clip_part_ceiling_for_stairs(test_cube, bdraw.temp_cubes, bdraw.temp_cubes2)) {
				for (auto c = bdraw.temp_cubes.begin(); c != bdraw.temp_cubes.end(); ++c) { // add floors after removing stairwells
					set_cube_zvals(*c, i->z1(), i->z2()); // set back to the correct zvals
					bdraw.add_section(*this, 0, *c, mat.roof_tex, roof_color, 4, 0, 1, is_house, 0); // only bottom surface
				}
				continue;
			}
		}
		bdraw.add_section(*this, 1, *i, mat.roof_tex, roof_color, 4, is_stacked_cube, skip_top, is_house, 0); // only Z dim
	} // for i
	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
		if (i->type == tquad_with_ix_t::TYPE_HELIPAD) {
			bdraw.add_tquad(*this, *i, bcube, building_texture_mgr.get_helipad_tid(), WHITE);
		}
		else if (i->type == tquad_with_ix_t::TYPE_SOLAR) {
			bdraw.add_tquad(*this, *i, bcube, building_texture_mgr.get_solarp_tid(), colorRGBA(0.6, 0.6, 0.6)); // panel is too bright compared to the roof, use a darker color
		}
		else if (i->type == tquad_with_ix_t::TYPE_TRIM) { // solar panel edges
			bdraw.add_tquad(*this, *i, bcube, tid_nm_pair_t(NO_SHADOW_WHITE_TEX), LT_GRAY); // untextured, no shadows
		}
		else if (is_house && (i->type == tquad_with_ix_t::TYPE_ROOF_PEAK || i->type == tquad_with_ix_t::TYPE_ROOF_SLOPE) && i->npts == 4) {
			// house peaked/sloped trapezoid roof: extend lower zvals out a bit
			tquad_with_ix_t tq(*i);
			float const extend(0.3*get_doorway_width());
			unsigned top_dim(2), bot_dim(2); // start at invalid values
			float top_lo(0.0), top_hi(0.0), bot_lo(0.0), bot_hi(0.0);
			
			// find the horizontal top and bottom edges
			for (unsigned n = 0; n < tq.npts; ++n) {
				point const &cur(i->pts[n]), &prev(i->pts[n ? n-1 : tq.npts-1]), &next(i->pts[(n == tq.npts-1) ? 0 : n+1]);
				if (cur.z < prev.z && cur.z == next.z) {bot_dim = (cur.x == next.x); bot_lo = min(cur[bot_dim], next[bot_dim]); bot_hi = max(cur[bot_dim], next[bot_dim]);}
				if (cur.z > prev.z && cur.z == next.z) {top_dim = (cur.x == next.x); top_lo = min(cur[top_dim], next[top_dim]); top_hi = max(cur[top_dim], next[top_dim]);}
			}
			if (top_dim == 2 || top_dim != bot_dim) {
				cout << "Bad house roof: " << TXT(top_dim) << TXT(bot_dim) << TXT(bcube.str()) << endl; // -464, -28 | -446, -13
				assert(0);
			}
			if (top_lo == bot_lo || top_hi == bot_hi) { // peaked, maybe clipped at one end, not hipped
				cube_t const tq_bcube(i->get_bcube());
				cube_t tq_bcube_lower(tq_bcube);
				tq_bcube_lower.z1() -= tq_bcube_lower.dz(); // extend downward so that it intersects the part below
				bool is_above_part(0);

				for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) { // includes garages and sheds, but not porches
					if (p->intersects_no_adj(tq_bcube_lower)) {is_above_part = 1; break;}
				}
				if (is_above_part) { // only extend if above a part (exclude porch roof)
					for (unsigned n = 0; n < tq.npts; ++n) { // extend edges downward and out first
						point &cur(tq.pts[n]);
						point const &orig_pt(i->pts[n]);
						point const other[2] = {i->pts[n ? n-1 : tq.npts-1], i->pts[(n == tq.npts-1) ? 0 : n+1]}; // prev and next; compare to unmodified points

						for (unsigned d = 0; d < 2; ++d) {
							if (cur.z     >= other[d].z) continue; // not along the bottom edge
							if (orig_pt.z == other[d].z) continue; // skip edge parallel to roofline
							vector3d ext_down_out; // extend down and out, but not along edge because it may be clipped diagonally to another roof tquad
							ext_down_out.z         = orig_pt.z         - other[d].z;
							ext_down_out[!top_dim] = orig_pt[!top_dim] - other[d][!top_dim];
							vector3d const delta((extend/ext_down_out.mag())*(orig_pt - other[d])); // use down + out dist for normalization, but move diagonally
							cur += delta;
							cube_t const new_bcube(tq.get_bcube());

							// if extended bcube intersects a part that the orig bcube didn't, undo the movement to avoid the roof clipping through another part
							for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
								if (!p->intersects_no_adj(tq_bcube_lower) && p->intersects_no_adj(new_bcube)) {cur -= delta; break;}
							}
						} // for d
					} // for n
					cube_t bot_edge_bcube;
					bool either_end_extended(0);

					for (unsigned n = 0; n < tq.npts; ++n) { // extend ends outward second
						point &cur(tq.pts[n]);
						float const top_lh[2] = {top_lo, top_hi}, bot_lh[2] = {bot_lo, bot_hi};

						for (unsigned d = 0; d < 2; ++d) {
							float const edge(top_lh[d]);
							if (edge != bot_lh[d] || cur[top_dim] != edge) continue; // vert not at end of roof
							float const extend_signed((d ? 1.0 : -1.0)*extend);

							if (edge != bcube.d[top_dim][d]) {
								cube_t test_cube(tq_bcube);
								test_cube.d[top_dim][ d] = edge + extend_signed;
								test_cube.d[top_dim][!d] = edge + 0.1*extend_signed; // minor shift to avoid intersecting the part this roof is placed on
								if (cube_int_parts_no_sec(test_cube)) continue;
								
								if (has_attic()) { // test the attic
									cube_t attic(get_attic_part());
									attic.z2() = bcube.z2();
									if (attic.intersects_no_adj(test_cube)) continue;
								}
							}
							cur[top_dim] += extend_signed; // extend out away from house
							either_end_extended = 1;
						} // for d
						// calculate bcube of points along bottom edge of roof section; required for intersecting/clipped roofs where bottom edge is shorter than top edge
						if (tq.pts[n].z < tq_bcube.zc()) {bot_edge_bcube.assign_or_union_with_pt(cur);}
					} // for n
					bot_edge_bcube.expand_in_dim(top_dim, -0.001*extend); // small inward bias to prevent Z-fighting with interior
					assert(!bot_edge_bcube.is_all_zeros()); // must have at least one point
					cube_t const new_bcube(tq.get_bcube());

					// add trim along the underside and edges of the roof to create rain gutters
					for (unsigned d = 0; d < 2; ++d) {
						float const old_edge(tq_bcube.d[!top_dim][d]), new_edge(new_bcube.d[!top_dim][d]);
						if (fabs(old_edge - new_edge) < 0.5*extend) continue; // not extended in this dir (can only extend in one dir per tquad)
						tquad_with_ix_t bot_surf(4, tquad_with_ix_t::TYPE_TRIM);
						UNROLL_4X(bot_surf.pts[i_].z = new_bcube.z1(););
						bot_surf.pts[0][!top_dim] = bot_surf.pts[1][!top_dim] = old_edge;
						bot_surf.pts[2][!top_dim] = bot_surf.pts[3][!top_dim] = new_edge;
						bot_surf.pts[0][ top_dim] = bot_surf.pts[3][ top_dim] = bot_edge_bcube.d[top_dim][0];
						bot_surf.pts[1][ top_dim] = bot_surf.pts[2][ top_dim] = bot_edge_bcube.d[top_dim][1];

						if (either_end_extended) { // currently always true, but could be false later if 3+ part houses are added
							// draw inside edge, which may be visible from above
							tquad_with_ix_t inside(bot_surf); // capture geometry before reverse()
							// move outward slightly to prevent Z-fighting with the interior wall
							for (unsigned n = 0; n < 2; ++n) {inside.pts[n][!top_dim] += 0.02*(new_edge - old_edge);}
							inside.pts[2] = inside.pts[1]; inside.pts[3] = inside.pts[0]; // move both points to the inside edge
							inside.pts[2].z = inside.pts[3].z = tq_bcube.z1() - 0.02*(tq_bcube.z1() - new_bcube.z1()); // top of edge, shifted down slightly
							if (d ^ top_dim ^ 1) {std::reverse(inside.pts, inside.pts+4);} // reverse to get the correct winding order
							tid_nm_pair_t const bot_tex(NO_SHADOW_WHITE_TEX); // untextured, no shadows
							bdraw.add_tquad(*this, inside, bcube, bot_tex, WHITE);
						}
						if (d ^ top_dim) {std::reverse(bot_surf.pts, bot_surf.pts+4);} // reverse to get the correct winding order
						tid_nm_pair_t const bot_tex(NO_SHADOW_WHITE_TEX); // untextured, no shadows
						bdraw.add_tquad(*this, bot_surf, bcube, bot_tex, WHITE);

						for (unsigned e = 0; e < 2; ++e) { // add triangle end caps
							tquad_with_ix_t end_cap(3, tquad_with_ix_t::TYPE_TRIM);
							UNROLL_3X(end_cap.pts[i_][top_dim] = bot_edge_bcube.d[top_dim][e];); // end
							end_cap.pts[0][!top_dim] = new_edge;
							end_cap.pts[1][!top_dim] = end_cap.pts[2][!top_dim] = old_edge;
							end_cap.pts[0].z = end_cap.pts[1].z = new_bcube.z1(); // bottom
							end_cap.pts[2].z = tq_bcube.z1(); // top
							if (d ^ e ^ top_dim ^ 1) {swap(end_cap.pts[0], end_cap.pts[1]);} // swap to get the correct winding order
							bdraw.add_tquad(*this, end_cap, bcube, bot_tex, WHITE);
						} // for e
						gutters.emplace_back(bot_surf.get_bcube(), (2*(!top_dim) + d)); // store 2*dim+dir in index - this is the roof edge the gutter connects to
					} // for d
				} // end is_above_part
			} // end peaked roof
			bdraw.add_tquad(*this, tq, bcube, mat.roof_tex.get_scaled_version(2.0), roof_color); // use roof texture
		} // end house roof quad
		else if (i->is_roof() || i->type == tquad_with_ix_t::TYPE_ROOF_ACC) {
			bdraw.add_tquad(*this, *i, bcube, mat.roof_tex.get_scaled_version(2.0), roof_color); // use roof texture
		}
		else if (i->type == tquad_with_ix_t::TYPE_BDOOR2 || i->type == tquad_with_ix_t::TYPE_RDOOR2) {
			bdraw.add_tquad(*this, *i, bcube, building_texture_mgr.get_bdoor2_tid(), WHITE);
		}
		else if (i->type == tquad_with_ix_t::TYPE_WALL) { // use wall texture
			bdraw.add_tquad(*this, *i, bcube, mat.side_tex, side_color);
		}
		else {assert(0);} // unsupported type
	} // for i
	for (auto i = details.begin(); i != details.end(); ++i) { // draw roof details
		if (i->type == DETAIL_OBJ_COLLIDER || i->type == DETAIL_OBJ_COLL_SHAD || i->type == DETAIL_OBJ_SHAD_ONLY) continue; // only drawn in the shadow pass

		if (i->type == ROOF_OBJ_AC) {
			bool const swap_st(i->dx() > i->dy());
			unsigned const tex_id(details.size() + parts.size() + mat_ix); // somewhat of a hash of various things; deterministic
			int const ac_tid(get_ac_unit_tid(tex_id));
			bdraw.add_cube(*this, *i, tid_nm_pair_t(ac_tid, -1, (swap_st ? 1.0 : -1.0), 1.0), WHITE, swap_st, 4, 1, 0, 0); // Z, skip bottom, ws_texture=0
			bdraw.add_cube(*this, *i, tid_nm_pair_t(ac_tid, -1, 0.3, 1.0), WHITE, 0, 3, 1, 0, 0); // XY with stretched texture, ws_texture=0
			continue;
		}
		if (i->type == ROOF_OBJ_DUCT) {
			int const duct_tid(building_texture_mgr.get_duct_tid()), vent_tid(building_texture_mgr.get_vent_tid());
			bool const duct_dim(i->dx() < i->dy()), swap_st_z(!duct_dim), has_vent(i->get_sz_dim(duct_dim) < 2.0*i->get_sz_dim(!duct_dim)); // add vent if short
			float ts[3] = {1.0, 1.0, 1.0};
			if (!has_vent) {ts[duct_dim] = 0.5;} // half texture/single section for the end

			for (unsigned dim = 0; dim < 3; ++dim) { // {x, y, z}
				bool const swap_st((dim == 2) ? swap_st_z : 0);
				int const tid((has_vent && bool(dim) == duct_dim) ? vent_tid : duct_tid);
				bdraw.add_cube(*this, *i, tid_nm_pair_t(tid, -1, ts[dim], 1.0), GRAY, swap_st, (1 << dim), 1, 0, 0); // skip_bottom, ws_texture=0
			}
			continue;
		}
		if (i->type == ROOF_OBJ_WTOWER) {
			bdraw.add_water_tower(*this, *i);
			continue;
		}
		bool const skip_bot(i->type != ROOF_OBJ_SCAP && i->type != ROOF_OBJ_SIGN && i->type != ROOF_OBJ_SIGN_CONN);
		bool const pointed(i->type == ROOF_OBJ_ANT); // draw antenna as a point
		building_t b(building_geom_t(4, rot_sin, rot_cos, pointed)); // cube
		b.bcube = bcube;
		tid_nm_pair_t tex;
		colorRGBA color;

		if (i->type == ROOF_OBJ_WALL && mat.add_windows) { // wall of brick/block building, use side color
			tex   = mat.side_tex;
			color = side_color;
		}
		else if (i->type == ROOF_OBJ_SIGN || i->type == ROOF_OBJ_SIGN_CONN) {
			color = WHITE; // also untextured, but casts shadows
		}
		else { // otherwise use roof color
			tex   = mat.roof_tex.get_scaled_version(1.5);
			color = detail_color*(pointed ? 0.5 : 1.0);
		}
		bdraw.add_section(b, 0, *i, tex, color, 7, skip_bot, 0, 1, 0); // all dims, no AO

		if (i->type == ROOF_OBJ_SIGN) { // add sign text
			bool const dim(i->dy() < i->dx());
			float const center_dim(i->get_center_dim(dim));
			bool dir(bcube.get_center_dim(dim) < center_dim); // choose dir based on center of building; inexact, will check parts below
			cube_t test_cube(*i);
			test_cube.expand_by(i->get_sz_dim(!dim)); // expand by sign width to include adjacency

			for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find part attached to the sign
				if (p->intersects(*i)) {dir = (p->get_center_dim(dim) < center_dim); break;}
			}
			vect_vnctcc_t &text_verts(bdraw.get_text_verts());
			unsigned text_verts_start(text_verts.size());
			cube_t text_bcube(*i);
			text_bcube.expand_in_dim(dim, 0.1*i->get_sz_dim(dim)); // expand outward a bit to reduce Z-fighting; doesn't seem to help much
			add_sign_text_verts_both_sides(name, text_bcube, dim, dir, text_verts);

			if (is_rotated()) {
				point const center(bcube.get_cube_center());
				// what about rotating the normal? seems difficult since normals are compressed, and maybe not very noticeable for sign text?
				for (auto i = text_verts.begin()+text_verts_start; i != text_verts.end(); ++i) {do_xy_rotate(center, i->v);}
			}
		}
	} // for i
	for (auto d = doors.begin(); d != doors.end(); ++d) {draw_building_ext_door(bdraw, *d, *this);} // draw exterior doors

	for (auto i = fences.begin(); i != fences.end(); ++i) {
		bdraw.add_fence(*this, *i, tid_nm_pair_t(WOOD_TEX, 0.4f/min(i->dx(), i->dy())), WHITE, (fences.size() > 1));
	}
	bool const skip_bottom(is_in_city); // skip_bottom=0, since it may be visible when extended over the terrain; okay to skip bottom for city driveways
	add_driveway_or_porch(bdraw, *this, driveway, LT_GRAY, skip_bottom);
	add_driveway_or_porch(bdraw, *this, porch,    LT_GRAY, 1); // skip_bottom=1

	if (roof_type == ROOF_TYPE_DOME || roof_type == ROOF_TYPE_ONION) {
		cube_t const &top(parts.back()); // top/last part
		point const center(top.get_cube_center());
		float const dx(top.dx()), dy(top.dy()), tscale(4.0f/(dx + dy));
		tid_nm_pair_t tex(mat.roof_tex); // use a different dome texture?
		tex.tscale_x *= tscale; tex.tscale_y *= tscale;
		bdraw.add_roof_dome(point(center.x, center.y, top.z2()), 0.5*dx, 0.5*dy, tex, roof_color*1.5, (roof_type == ROOF_TYPE_ONION));
	}
}

void building_t::get_detail_shadow_casters(building_draw_t &bdraw) {
	for (auto const &i : details) {
		if (i.type == DETAIL_OBJ_COLL_SHAD || i.type == DETAIL_OBJ_SHAD_ONLY) {bdraw.add_cube(*this, i, tid_nm_pair_t(), WHITE);} // always a cube
	}
}

void building_t::get_all_drawn_ext_wall_verts(building_draw_t &bdraw) {
	//if (interior == nullptr) return; // only needed if building has an interior?
	if (!is_valid()) return; // invalid building
	building_mat_t const &mat(get_material());
	ext_side_qv_range.draw_ix = bdraw.get_to_draw_ix(mat.wall_tex);
	ext_side_qv_range.start   = bdraw.get_num_verts (mat.wall_tex);

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // multiple cubes/parts/levels, room parts only, no AO
		if (!is_basement(i)) {bdraw.add_section(*this, 1, *i, mat.wall_tex, wall_color, 3, 0, 0, 1, 0);} // XY
	}
	ext_side_qv_range.end = bdraw.get_num_verts(mat.wall_tex);
	get_basement_ext_wall_verts(bdraw);
}
void building_t::get_basement_ext_wall_verts(building_draw_t &bdraw) const {
	if (!has_basement()) return;
	tid_nm_pair_t const tp(get_basement_wall_texture());
	// find basement door and exclude it
	// dim_mask bits: enable dims: 1=x, 2=y, 4=z | disable cube faces: 8=x1, 16=x2, 32=y1, 64=y2, 128=z1, 256=z2
	unsigned dim_mask(3); // XY
	cube_t const &basement(get_basement());

	if (has_basement_door) { // remove a section of wall around the basement door; can't use stencil test associated with ext_side_qv_range
		assert(interior);
		door_t const &door(interior->get_ext_basement_door());
		unsigned const this_face(1 << (2*interior->extb_wall_dim + interior->extb_wall_dir + 3));
		dim_mask |= this_face; // skip this face for the full basement call below

		for (unsigned d = 0; d < 2; ++d) {
			cube_t side(basement);
			side.d[!door.dim][!d] = door.d[!door.dim][d];
			bdraw.add_section(*this, 0, side, tp, WHITE, ~(this_face | 4), 0, 0, 1, 0); // single face, always white
		}
		if (door.z2() < basement.z2()) { // door shorter than basement; can happen with multi-level parking garages
			cube_t top(basement);
			for (unsigned d = 0; d < 2; ++d) {top.d[!door.dim][d] = door.d[!door.dim][d];} // same range as door in !door.dim
			top.z1() = door.z2();
			bdraw.add_section(*this, 0, top, tp, WHITE, ~(this_face | 4), 0, 0, 1, 0); // single face, always white
		}
	}
	bdraw.add_section(*this, 0, basement, tp, WHITE, dim_mask, 0, 0, 1, 0); // XY, always white
}

void set_skip_faces_for_nearby_cube_edge(cube_t const &c, cube_t const &C, float dist, bool dim, unsigned &dim_mask) {
	for (unsigned dir = 0; dir < 2; ++dir) { // skip faces along the edges of the building bcube or along an extended basement exterior facing wall
		if (fabs(c.d[dim][dir] - C.d[dim][dir]) < dist) {dim_mask |= (1<<(2*dim+dir+3));}
	}
}
void building_t::get_all_drawn_interior_verts(building_draw_t &bdraw) {
	if (!is_valid() || interior == nullptr) return; // invalid building or no interior
	building_mat_t const &mat(get_material());
	auto const parts_end(get_real_parts_end());
	float const floor_thickness(get_floor_thickness()), fc_thickness(get_fc_thickness());
	bdraw.begin_draw_range_capture();

	for (auto i = interior->floors.begin(); i != interior->floors.end(); ++i) { // 600K T
		tid_nm_pair_t tex;
		colorRGBA const color(get_floor_tex_and_color(*i, tex));
		// expand_by_xy(-get_trim_thickness()) to prevent z-fighting when AA is disabled? but that will leave small gaps where floors from adjacent parts meet
		bdraw.add_section(*this, 0, *i, tex, color, 4, 1, 0, 1, 0); // no AO; skip_bottom; Z dim only
	}
	for (auto i = interior->ceilings.begin(); i != interior->ceilings.end(); ++i) { // 600K T
		// skip top surface of all but top floor ceilings if the roof is sloped;
		// if this is an office building, the ceiling could be at a lower floor with a flat roof even if the highest floor has a sloped roof, so we must skip it
		bool skip_top(skip_top_of_ceilings());

		if (!skip_top) { // check if this is a top ceiling; needed for light occlusion
			float const toler(floor_thickness);
			skip_top = 1; // assume it's not

			for (auto p = parts.begin(); p != parts_end; ++p) { // Note: excludes garages and sheds
				if (!is_basement(p) && p->contains_cube_xy(*i) && fabs(i->z2() - p->z2()) < toler) {skip_top = 0; break;}
			}
		}
		tid_nm_pair_t tex;
		colorRGBA const color(get_ceil_tex_and_color(*i, tex));
		bdraw.add_section(*this, 0, *i, tex, color, 4, 0, skip_top, 1, 0); // no AO; Z dim only
	} // for i
	// minor optimization: don't need shadows for ceilings because lights only point down; assumes ceil_tex is only used for ceilings; not true for all houses/apts/hotels
	if (!is_residential()) {bdraw.set_no_shadows_for_tex(mat.ceil_tex);}
	float const wall_thickness(get_wall_thickness()), extb_wall_thresh(1.1*wall_thickness); // extb_wall_thresh uses wall thickness + tolerance

	for (unsigned dim = 0; dim < 2; ++dim) { // Note: can almost pass in (1U << dim) as dim_filt, if it wasn't for door cutouts (2.2M T)
		for (auto i = interior->walls[dim].begin(); i != interior->walls[dim].end(); ++i) {
			//unsigned const dim_mask(1 << dim); // doesn't work with office building hallway intersection corners and door frame shadows
			unsigned dim_mask(3); // XY
			set_skip_faces_for_nearby_cube_edge(*i, bcube, wall_thickness, !dim, dim_mask); // easy case: skip faces along the edges of the building bcube
			bool const in_basement(i->z1() < ground_floor_z1);
			bool const in_ext_basement(in_basement && i >= (interior->walls[dim].begin() + interior->extb_walls_start[dim]));
			
			// check rooms; skip this for above ground complex floorplans because they may have unexpected wall ends visible at non-rectangular rooms
			if (!has_complex_floorplan || in_basement) {
				unsigned const extb_room_start((interior->ext_basement_hallway_room_id >= 0) ? interior->ext_basement_hallway_room_id : interior->rooms.size());
				unsigned const rooms_start(in_ext_basement ? extb_room_start : 0);
				unsigned const rooms_end  (in_ext_basement ? interior->rooms.size() : extb_room_start);
			
				for (auto r = interior->rooms.begin()+rooms_start; r != interior->rooms.begin()+rooms_end; ++r) {
					if (!r->intersects(*i)) continue; // wall doesn't intersect this room
					// office hallways can have outside corners, and we need to draw the walls there
					if (!is_house && has_pri_hall() && !in_basement && r->is_hallway) {dim_mask = 3; break;} // force all 4 sides
					set_skip_faces_for_nearby_cube_edge(*i, *r, wall_thickness, !dim, dim_mask); // wall ends

					// ext basement rooms don't need to have their exterior wall surfaces drawn, but only valid for walls not shared between hallways and connected rooms
					if (in_ext_basement && !interior->has_mall) { // also, skip for malls, because this doesn't work with store separators
						if (fabs(i->d[!dim][0] - r->d[!dim][0]) < extb_wall_thresh && fabs(i->d[!dim][1] - r->d[!dim][1]) < extb_wall_thresh) {
							// use slightly more than half wall_thickness here so that we pick up the edges of the current room but not nearby adjacent rooms;
							// see building_t::is_basement_room_placement_valid() wall_expand_toler
							set_skip_faces_for_nearby_cube_edge(*i, *r, 0.6*wall_thickness, dim, dim_mask); // wall side edges
						}
					}
				} // for r
			}
			if (is_cube() && !in_basement && dim_mask != (1U << dim)) { // disable interior walls at building exteriors for cube buildings if we still have ends enabled
				for (auto p = parts.begin(); p != parts_end; ++p) {
					if (!p->contains_cube(*i)) continue;

					for (unsigned dir = 0; dir < 2; ++dir) { // skip faces along part exteriors
						if (fabs(i->d[!dim][dir] - p->d[!dim][dir]) > wall_thickness) continue;
						point test_pt;
						test_pt.z     = i->z1() + wall_thickness; // ground floor, since this will be the widest point
						test_pt[ dim] = i->get_center_dim(dim); // wall centerline
						test_pt[!dim] = p->d[!dim][dir] + (dir ? 1.0 : -1.0)*wall_thickness; // extend slightly away from the exterior wall
						bool contained(0);

						for (auto p2 = parts.begin(); p2 != parts_end; ++p2) { // only needed if has_small_part?
							if (p2 != p && p2->contains_pt(test_pt)) {contained = 1; break;}
						}
						if (!contained) {dim_mask |= (1<<(2*(!dim)+dir+3));}
					} // for dir
				} // for p
			}
			if (check_skylight_intersection(*i)) {dim_mask |= 4;} // draw top surface if under skylight
			colorRGBA color(in_basement ? WHITE : wall_color); // basement walls are always white
			tid_nm_pair_t tex;
			
			if (!is_house && in_ext_basement) {
				if (has_mall()) {tex = mat.wall_tex; color = colorRGBA(1.0, 0.9, 0.7);} // light brown/yellow-ish
				else            {tex = get_concrete_texture();} // office building extended basement
			}
			else {tex = mat.wall_tex;}
			bdraw.add_section(*this, 0, *i, tex, color, dim_mask, 1, 0, 1, 0); // no AO; X and/or Y dims only, skip bottom, only draw top if under skylight
		} // for i
	} // for dim
	// add the walkway interiors; these are outside the building and may get culled when outside the view frustum
	if (DRAW_WALKWAY_INTERIORS) {
		for (building_walkway_t &w : walkways) {get_walkway_interior_verts(bdraw, w);}
	}
	// Note: stair/elevator landings can probably be drawn in room_geom along with stairs, though I don't think there would be much benefit in doing so
	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) { // added per-floor (530K T)
		if (i->in_mall) continue; // mall stairs are open; skip drawing of landing interior faces
		unsigned dim_mask(3); // disable faces: 8=x1, 16=x2, 32=y1, 64=y2
		if (i->for_elevator) {dim_mask |= (120 - (1 << (i->get_face_id() + 3)));} // disable all but the open side of the elevator
		else if (i->for_ramp) {
			cube_t const &basement(get_basement());
			if (i->dim == 0) {dim_mask |= ((i->yc() < basement.yc()) ? 32 : 64);} // X, skip Y face closest to building wall
			else             {dim_mask |= ((i->xc() < basement.xc()) ?  8 : 16);} // Y, skip X face closest to building wall
		}
		if (!is_house && has_ext_basement() && point_in_extended_basement_not_basement(i->get_cube_center())) { // concrete edges for office ext basements
			bdraw.add_section(*this, 0, *i, get_concrete_texture(), WHITE, dim_mask, 0, 0, 1, 0, 0.0, 0, 1.0, 1); // no AO; X/Y dims only, inverted normals
		}
		else { // else use the wall texture
			bdraw.add_section(*this, 0, *i, mat.wall_tex, mat.wall_color, dim_mask, 0, 0, 1, 0, 0.0, 0, 1.0, 1); // no AO; X/Y dims only, inverted normals
		}
	} // for i
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		bool const dim(i->dim), dir(i->dir);
		float const spacing(i->get_wall_thickness()), frame_width(i->get_frame_width()); // space between inner/outer walls + frame around door
		unsigned dim_mask(3); // x and y dims enabled
		dim_mask |= (1 << (i->get_door_face_id() + 3)); // disable the face for the door opening
		cube_t shaft(*i);
		shaft.z2() -= ELEVATOR_Z2_SHIFT*fc_thickness*(i->under_skylight ? 1.0 : 0.25); // avoid clipping through skylights
		bdraw.add_section(*this, 0, shaft, mat.wall_tex, mat.wall_color, dim_mask, 0, 0, 1, 0); // outer elevator is textured like the walls
		cube_t entrance(shaft);
		entrance.d[dim][!dir] = entrance.d[dim][dir] + (dir ? -1.0f : 1.0f)*spacing; // set correct thickness

		for (unsigned d = 0; d < 2; ++d) { // add frame on both sides of the door opening
			cube_t frame(entrance); // one side
			frame.d[!dim][d] = frame.d[!dim][!d] + (d ? 1.0f : -1.0f)*frame_width; // set position
			unsigned dim_mask2(3); // x and y dims enabled
			dim_mask2 |= (1 << (2*(!dim) + (!d) + 3)); // 3 faces drawn
			bdraw.add_section(*this, 0, frame, mat.wall_tex, mat.wall_color, dim_mask2, 0, 0, 1, 0);
		}
		cube_t inner_cube(shaft);
		inner_cube.expand_by_xy(-spacing);
		// add interior of elevator by drawing the inside of the cube with a slightly smaller size, with invert_normals=1; normal mapped?
		tid_nm_pair_t wall_panel_tex(FENCE_TEX, -1, 16.0, 16.0);
		wall_panel_tex.set_specular(0.1, 50.0);
		bdraw.add_section(*this, 0, inner_cube, wall_panel_tex, WHITE, dim_mask, 0, 0, 1, 0, 0.0, 0, 1.0, 1);

		if (i->under_skylight ) { // under skylight, draw the top
			bdraw.add_section(*this, 0, shaft, mat.wall_tex, WHITE, 4, 1, 0, 1, 0); // top surface only
		}
		// Note elevator doors are dynamic and are drawn as part of room_geom
	} // for i
	if (has_attic()) {
		// add inside surface of attic access hole; could be draw as room geom if needed
		bdraw.add_section(*this, 0, interior->attic_access, mat.wall_tex, mat.wall_color, 3, 0, 0, 1, 0, 0.0, 0, 1.0, 1); // no AO; X/Y dims only, inverted normals
	}
	// Note: interior doors are drawn as part of room_geom
	bdraw.end_draw_range_capture(interior->draw_range); // 80MB, 394MB, 836ms
}

void building_t::get_walkway_interior_verts(building_draw_t &bdraw, building_walkway_t &w) {
	w.windows.clear(); // in case this gets called multiple times on the same walkway
	w.frames .clear();
	if (!w.is_owner) return;
	assert(w.bcube.z1() >= ground_floor_z1); // must be above ground
	building_mat_t const &mat(get_material());
	auto const parts_end(get_real_parts_end());
	float const floor_spacing(get_window_vspace()), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	unsigned const bot_floor(round_fp((w.bcube.z1() - ground_floor_z1)/floor_spacing)), top_floor(round_fp((w.bcube.z2() - ground_floor_z1)/floor_spacing));
	assert(bot_floor < top_floor); // must be at least one floor
	cube_t ww_floor(w.bcube), ww_ceil(w.bcube);
	vect_cube_t ww_walls;

	for (unsigned f = bot_floor; f < top_floor; ++f) {
		float const zval(ground_floor_z1 + f*floor_spacing), next_zval(zval + floor_spacing);
		set_cube_zvals(ww_floor, zval, zval+fc_thick);
		set_cube_zvals(ww_ceil,  next_zval-fc_thick, next_zval);
		bdraw.add_section(*this, 0, ww_floor, mat.floor_tex, mat.floor_color, 4, 1, 0, 1, 0); // no AO; top only
		bdraw.add_section(*this, 0, ww_ceil,  mat.ceil_tex,  mat.ceil_color,  4, 0, 1, 1, 0); // no AO; bottom only; applies to apartments and hotels as well
	}
	// walls on all 4 sides; walls extend through all floors; unlike normal exterior walls, these are windowless and have thickness
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned d = 0; d < 2; ++d) { // dir
			bool const is_end_dim(bool(dim) == w.dim), add_ext_door(is_end_dim && w.has_ext_door(d));
			if (is_end_dim && w.open_ends[d]) continue; // no wall (or door) at this end
			float const wall_thick(add_ext_door ? /*get_door_shift_dist()*/0.5*wall_thickness : wall_thickness); // flush with exterior door? then there's a light gap
			cube_t wall(w.bcube);
			wall.d[dim][!d] = w.bcube.d[dim][d] + (d ? -1.0 : 1.0)*wall_thick;
			unsigned const base_dim_mask(1 << unsigned(dim));
			unsigned dim_mask(base_dim_mask | (1<<(2*dim+d+3))); // only inside face in dim !w.dim should be visible

			if (is_end_dim) { // draw ends
				//if (has_int_windows()) {dim_mask = base_dim_mask;} // draw opposite of ends to block unwanted interior windows; but then we can get narrow strips

				if (add_ext_door) { // draw half wall to either side of door
					for (unsigned s = 0; s < 2; ++s) { // {left, right} side
						cube_t side(wall);
						side.d[!dim][!s] = w.door_bounds[d][s];
						bdraw.add_section(*this, 0, side, mat.wall_tex, wall_color, dim_mask, 1, 1, 1, 0); // no AO; skip bottom and top
					}
				}
				else { // draw full wall
					bdraw.add_section(*this, 0, wall, mat.wall_tex, wall_color, dim_mask, 1, 1, 1, 0); // no AO; skip bottom and top
				}
			}
			else { // draw sides with (real) cutouts for windows
				ww_walls.clear();
				ww_walls.push_back(wall);

				if (!w.elevator_cut.is_all_zeros() && w.elevator_cut.intersects(wall)) {
					subtract_cube_from_cubes(w.elevator_cut, ww_walls);
					// add frame around each floor of the elevator
					float const frame_thickness(1.5*get_trim_thickness()), frame_width(1.5*get_trim_height()), frame_dz(fc_thick + frame_thickness);
					cube_t cutout(wall);
					cutout.intersect_with_cube(w.elevator_cut);

					for (unsigned f = bot_floor; f <= top_floor; ++f) { // add horizontal frames for each floor and at the top
						float const zval(ground_floor_z1 + f*floor_spacing);
						cube_t frame(cutout);
						set_cube_zvals(frame, max(cutout.z1(), zval-frame_dz), min(cutout.z2(), zval+frame_dz));
						frame.expand_in_dim(w.dim, -frame_width); // clip off vertical frame width
						w.frames.push_back(frame);
					}
					for (unsigned d = 0; d < 2; ++d) {
						cube_t frame(cutout);
						frame.d[w.dim][!d] = cutout.d[w.dim][d] + (d ? -1.0 : 1.0)*frame_width;
						w.frames.push_back(frame);
					}
				}
				for (cube_t &wseg : ww_walls) {
					float const length(wseg.get_sz_dim(w.dim));
					unsigned num_segs(length/floor_spacing);

					if (num_segs == 0) { // short segment with no window
						float const zval(ground_floor_z1 + bot_floor*floor_spacing);
						set_cube_zvals(wseg, zval, (zval + (top_floor - bot_floor)*floor_spacing));
						bdraw.add_section(*this, 0, wseg, mat.wall_tex, wall_color, dim_mask, 1, 1, 1, 0); // no AO; skip bottom and top
						continue;
					}
					if (!(num_segs & 1)) {++num_segs;} // must be odd
					float const seg_len(length/num_segs), first_seg_end(wseg.d[!dim][0] + seg_len);
					cube_t seg(wseg);
					seg.d[!dim][1] = first_seg_end;
					dim_mask |= (1 << unsigned(!dim)); // draw edges

					for (unsigned n = 0; n < num_segs; n += 2) { // alternate seg - window - seg
						bdraw.add_section(*this, 0, seg, mat.wall_tex, wall_color, dim_mask, 1, 1, 1, 0); // no AO; skip bottom and top
						seg.translate_dim(!dim, 2.0*seg_len);
					}
					for (unsigned f = bot_floor; f < top_floor; ++f) { // add top and bottom strips for each floor
						float const zval(ground_floor_z1 + f*floor_spacing), next_zval(zval + floor_spacing);
						float const wz1(zval+0.30*floor_spacing), wz2(next_zval -0.15*floor_spacing);
						cube_t bot(wseg), top(wseg);
						set_cube_zvals(bot, zval, wz1);
						set_cube_zvals(top, wz2,  next_zval);
						bdraw.add_section(*this, 0, bot, mat.wall_tex, wall_color, (dim_mask | 4), 1, 0, 1, 0); // no AO; skip bottom
						bdraw.add_section(*this, 0, top, mat.wall_tex, wall_color, (dim_mask | 4), 0, 1, 1, 0); // no AO; skip top
						// add windows
						cube_with_ix_t window(wseg, (2*dim + d));
						set_cube_zvals(window, wz1, wz2);
						window.d[!dim][0] = first_seg_end;
						window.d[!dim][1] = first_seg_end + seg_len;

						for (unsigned n = 0; n < num_segs/2; ++n) {
							w.windows.push_back(window);
							window.translate_dim(!dim, 2.0*seg_len);
						}
					} // for f
				} // for wseg
			} // end drawn sides
		} // for d
	} // for dim
}

template<typename T> void building_t::add_door_verts(cube_t const &D, T &drawer, door_rotation_t &drot, uint8_t door_type, bool dim,
	bool dir, float open_amt, bool opens_out, bool exterior, bool on_stairs, bool hinge_side, bool open_min_amt, bool draw_top_edge) const
{
	bool const is_rooftop_door(door_type == tquad_with_ix_t::TYPE_RDOOR);
	int type(tquad_with_ix_t::TYPE_IDOOR); // use interior door type, even for exterior door, because we're drawing it in 3D inside the building
	if (is_rooftop_door) {type = door_type;} // roof door is special because we only draw half of the texture when open
	bool const opened(open_amt > 0.0), opens_up(door_type == tquad_with_ix_t::TYPE_GDOOR);
	// exclude the frame on open interior doors
	bool const exclude_frame((door_type == tquad_with_ix_t::TYPE_HDOOR || door_type == tquad_with_ix_t::TYPE_ODOOR || is_rooftop_door) && (!exterior || opened));
	unsigned const num_edges(opens_up ? 4 : 2);
	int const tid(get_building_door_tid(door_type));
	float const half_thickness(opens_up ? 0.01*D.dz() : 0.5*DOOR_THICK_TO_WIDTH*D.get_sz_dim(!dim));
	unsigned const num_sides((door_type == tquad_with_ix_t::TYPE_BDOOR || door_type == tquad_with_ix_t::TYPE_BDOOR2) ? 2 : 1); // double doors for office building exterior
	tid_nm_pair_t const tp(tid, -1, 1.0f/num_sides, 1.0); // map full texture in Y
	colorRGBA const &color((exterior && !opens_up) ? door_color : WHITE); // garage doors are always white

	for (unsigned side = 0; side < num_sides; ++side) { // {right, left}
		cube_t dc(D);
		if (num_sides == 2) {dc.d[!dim][bool(side) ^ dim ^ dir ^ 1] = 0.5f*(D.d[!dim][0] + D.d[!dim][1]);} // split door in half
		// we don't want to draw the open stairs door because it may get in the way, but we need to overwrite the previous verts, so make it zero area
		if (opened && on_stairs) {dc.z2() = dc.z1();}
		bool const int_other_side(exterior ? 0 : hinge_side), swap_sides(exterior ? (side == 0) : hinge_side); // swap sides for right half of exterior door
		// 0,1: bottom, 2,3: top; we pass in the same drot for both sides because the value is only filled in and used for interior doors, which have only one side
		tquad_with_ix_t const door(set_door_from_cube(dc, dim, dir, type, 0.0, exterior, open_amt, opens_out, opens_up, swap_sides, open_min_amt, drot));
		vector3d const normal(door.get_norm());
		tquad_with_ix_t door_edges[4] = {door, door, door, door}; // most doors will only use 2 of these

		for (unsigned d = 0; d < 2; ++d) { // draw front and back sides
			tquad_with_ix_t door_side(door);
			vector3d const offset((d ? -1.0 : 1.0)*half_thickness*normal);
			for (unsigned n = 0; n < 4; ++n) {door_side.pts[n] += offset;}

			for (unsigned e = 0; e < num_edges; ++e) {
				unsigned const ixs[4][2] = {{1, 2}, {3, 0}, {0, 1}, {2, 3}};
				door_edges[e].pts[2*d+1] = door_side.pts[ixs[e][ d]];
				door_edges[e].pts[2*d+0] = door_side.pts[ixs[e][!d]];
			}
			if (d == 1) { // back/inside face of house/office door
				swap(door_side.pts[0], door_side.pts[1]);
				swap(door_side.pts[2], door_side.pts[3]);
				door_side.type = (is_house ? (unsigned)tquad_with_ix_t::TYPE_IDOOR_IN :
					(is_rooftop_door ? (unsigned)tquad_with_ix_t::TYPE_RDOOR_IN : (unsigned)tquad_with_ix_t::TYPE_ODOOR_IN));
			}
			drawer.add_tquad(*this, door_side, bcube, tp, color, (int_other_side && !opened), exclude_frame, 0); // invert_tc_x=xxx, no_tc=0
		} // for d
		if (opened || on_stairs) { // add untextured door edges; only needed for open doors or doors at the bottom of basement stairs
			for (unsigned e = 0; e < num_edges; ++e) {
				drawer.add_tquad(*this, door_edges[e], bcube, tp, color, 0, 0, 1); // invert_tc_x=0, exclude_frame=0, no_tc=1, use single texel from corner of door texture
			}
		}
		if (opened && !exterior && !opens_up && num_sides == 1 && (draw_top_edge || check_skylight_intersection(door.get_bcube()))) {
			// open interior door at skylight or tall room; draw top edge of door
			tquad_with_ix_t top_edge(4, door.type);
			top_edge.pts[0] = door_edges[0].pts[0];
			top_edge.pts[1] = door_edges[0].pts[3];
			top_edge.pts[2] = door_edges[1].pts[2];
			top_edge.pts[3] = door_edges[1].pts[1];
			drawer.add_tquad(*this, top_edge, bcube, tp, color, 0, 0, 1); // invert_tc_x=0, exclude_frame=0, no_tc=1, use single texel from corner of door texture
		}
	} // for side
}

// explicit template specialization
template void building_t::add_door_verts(cube_t const &D, building_room_geom_t &drawer, door_rotation_t &drot, uint8_t door_type, bool dim,
	bool dir, float open_amt, bool opens_out, bool exterior, bool on_stairs, bool hinge_side, bool open_min_amt, bool draw_top_edge) const;

// Note: this is actually the geometry of walls that have windows, not the windows themselves
void building_t::get_all_drawn_window_verts(building_draw_t &bdraw, bool lights_pass, float offset_scale,
	point const *only_cont_pt_in, bool no_skylights, bool draw_int_windows) const
{
	if (!is_valid()) return; // invalid building

	if (!no_skylights && !lights_pass) { // draw skylight glass
		for (cube_t const &skylight : skylights) {
			tid_nm_pair_t tp;
			tp.transparent = 1; // doesn't do anything?
			cube_t glass(skylight);
			float const ceil_thickness(glass.dz());
			glass.z1() += 0.50*ceil_thickness; // glass pane is only 25% of ceiling thickness
			glass.z2() -= 0.25*ceil_thickness;
			bdraw.add_cube(*this, glass, tp, colorRGBA(WHITE, 0.1), 0, 4, 0, 0, 0); // top and bottom only, untextured
		} // for skylight
	}
	point const only_cont_pt(only_cont_pt_in ? get_inv_rot_pos(*only_cont_pt_in) : all_zeros);
	bool const cut_door_holes(only_cont_pt_in != nullptr); // needed even for non-cube buildings, so capture before clearing
	if (!is_cube()) {only_cont_pt_in = nullptr;} // not needed for current non-cube buildings (optimization)
	building_mat_t const &mat(get_material());
	bool const draw_windows(draw_int_windows ? has_int_windows() : has_windows());

	// Note: city office buildings don't have add_windows set because windows don't align with their interior rooms/floors,
	// which means the player currently can't see into or out of the building; but we can still cut out window holes on the interior side
	if (!global_building_params.windows_enabled() || (lights_pass ? !mat.add_wind_lights : !draw_windows)) {
		// no windows for this material (set in building_materials.txt)
		if (cut_door_holes) {cut_holes_for_ext_doors(bdraw, only_cont_pt, 0xFFFF);} // still need to draw holes for doors
		return;
	}
	int const window_tid(building_texture_mgr.get_window_tid());
	if (window_tid < 0) return; // not allocated - error?
	if (mat.wind_xscale == 0.0 || mat.wind_yscale == 0.0) return; // no windows for this material?
	tid_nm_pair_t tex;
	colorRGBA color;

	if (draw_int_windows && !has_windows()) { // calculate interior window spacing that aligns with actual floor spacing
		float const window_tx(0.5*mat.floorplan_wind_xscale), window_ty(0.5/get_window_vspace());
		tex = tid_nm_pair_t(window_tid, -1, window_tx, window_ty);
	}
	else {
		tex = tid_nm_pair_t(window_tid, -1, mat.get_window_tx(), mat.get_window_ty(), mat.wind_xoff, -mat.wind_yoff); // Note: wind_yoff is negated
	}
	if (lights_pass) { // slight yellow-blue tinting using bcube x1/y1 as a hash
		float const tint(0.2*fract(100.0f*(bcube.x1() + bcube.y1())));
		color = colorRGBA((1.0 - tint), (1.0 - tint), (0.8 + tint), 1.0);
	} else {color = mat.window_color;}
	int const clip_windows(draw_windows ? (is_house ? 2 : 1) : 0);
	float const floor_spacing(get_window_vspace()), first_floor_z1(ground_floor_z1 + floor_spacing);
	float const gf_door_ztop(doors.empty() ? 0.0f : (EXACT_MULT_FLOOR_HEIGHT ? first_floor_z1 : doors.front().pts[2].z));
	unsigned draw_parts_mask(0);
	bool room_with_stairs(0);
	cube_t cont_part; // part containing the point

	if (only_cont_pt_in) {
		cont_part        = get_part_containing_pt(only_cont_pt);
		room_with_stairs = room_or_adj_room_has_stairs(get_room_containing_pt(only_cont_pt), only_cont_pt.z, 1, 1); // inc_adj_rooms=1, check_door_open=1
	}
	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // multiple cubes/parts/levels, excluding chimney/porch/etc.
		if (is_basement(i)) continue; // skip the basement
		bool const split_per_floor(i == parts.begin() && floor_ext_door_mask > 1); // for multi-family houses
		unsigned const num_splits(split_per_floor ? calc_num_floors(*i, floor_spacing, get_floor_thickness()) : 1);

		for (unsigned f = 0; f < num_splits; ++f) {
			float const floor_offset(f*floor_spacing), slice_z1(i->z1() + floor_offset), door_ztop(gf_door_ztop + floor_offset);
			cube_t part(*i), draw_part;
			cube_t const *clamp_cube(nullptr);
			set_cube_zvals(part, slice_z1, (split_per_floor ? (slice_z1 + floor_spacing) : i->z2()));

			if (only_cont_pt_in && *i != cont_part && !i->contains_pt(only_cont_pt)) { // not the part containing the point
				float const z_exp(get_fc_thickness()); // allow a bit of extra Z overlap, which helps when the player is on the stairs

				if (i->contains_pt_xy(only_cont_pt) && only_cont_pt.z > i->z1()-z_exp && only_cont_pt.z < i->z2()+z_exp) {} // okay, can draw unsplit in this case
				else if (room_with_stairs && are_parts_stacked(*i, cont_part)) { // windows may be visible through stairs in rooms with stacked parts
					draw_part  = cont_part;
					draw_part.intersect_with_cube_xy(part);
					clamp_cube = &draw_part;
				}
				else {
					if (i->z2() < only_cont_pt.z || i->z1() > only_cont_pt.z) continue; // z-range not contained, skip
					bool skip(0);

					for (unsigned d = 0; d < 2; ++d) {
						if (i->d[ d][0] != cont_part.d[ d][1] && i->d[ d][1] != cont_part.d[ d][0]) continue; // not adj in dim d
						if (i->d[!d][0] >= cont_part.d[!d][1] || i->d[!d][1] <= cont_part.d[!d][0]) continue; // no overlap in dim !d
						if (i->d[!d][1] < only_cont_pt[!d] || i->d[!d][0] > only_cont_pt[!d]) {skip = 1; break;} // other dim range not contained, skip
						draw_part = part; // deep copy
						max_eq(draw_part.d[!d][0], cont_part.d[!d][0]); // clamp to contained part in dim !d
						min_eq(draw_part.d[!d][1], cont_part.d[!d][1]);
						clamp_cube = &draw_part;
						break;
					} // for d
					if (skip || clamp_cube == nullptr) continue; // skip if adj in neither dim, always skip (but could check chained adj case)
				}
			}
			// skip windows on sides with doors, but only for buildings with windows
			unsigned const part_ix(i - parts.begin()), dsides((part_ix < 4 && draw_windows && i->z1() == ground_floor_z1) ? door_sides[part_ix] : 0);
			bdraw.add_section(*this, 1, part, tex, color, 3, 0, 0, 1, clip_windows, door_ztop, dsides, offset_scale, 0, clamp_cube); // XY, no_ao=1
			draw_parts_mask |= (1 << part_ix);

			// add ground floor windows next to doors
			if (!is_cube())  continue; // below logic is only correct for cube-shaped buildings; other shapes are generally convex anyway
			if (dsides == 0) continue; // no doors
			float const space(0.25*floor_spacing), toler(0.1*floor_spacing);

			for (unsigned dim = 0; dim < 2; ++dim) {
				unsigned const num_windows(get_num_windows_on_side(i->d[!dim][0], i->d[!dim][1]));
				if (num_windows <= 1) continue; // no space to split the windows on this wall
				float const window_spacing(i->get_sz_dim(!dim)/num_windows), side_lo(i->d[!dim][0]), side_hi(i->d[!dim][1]);

				for (unsigned dir = 0; dir < 2; ++dir) {
					if (!(dsides & (1 << (2*dim + dir)))) continue; // no door on this side
					unsigned const dim_mask((1 << dim) + (1 << (3 + 2*dim + (1-dir)))); // enable only this dim but disable the other dir
					float const wall_pos(i->d[dim][dir]);
					vector<float> &wall_edges(bdraw.temp_wall_edges);
					wall_edges.clear();

					for (auto d = doors.begin(); d != doors.end(); ++d) {
						cube_t const c(d->get_bcube());
						if (!is_house && c.z1() > first_floor_z1) continue; // not ground floor door and not house upper door; walkway door
						if ((c.dy() < c.dx()) != dim)             continue; // wrong dim
						if (c.d[dim][0]-toler > wall_pos || c.d[dim][1]+toler < wall_pos) continue; // door not on this wall
						float const door_lo(c.d[!dim][0]), door_hi(c.d[!dim][1]);
						if (door_lo > side_hi || door_hi < side_lo)     continue; // door not on this part
						if (c.z1() >= part.z2() || c.z2() <= part.z1()) continue; // door not on this floor slice
						// align to an exact multiple of window period so that bottom floor windows line up with windows on the floors above and no walls are clipped
						if (wall_edges.empty()) {wall_edges.push_back(side_lo); wall_edges.push_back(side_hi);} // first wall, add end points
						wall_edges.push_back(door_lo - space); // low
						wall_edges.push_back(door_hi + space); // high
					} // for d
					if (wall_edges.empty()) { // no door, could be a non-main door (roof access, garage, shed) or slice with a door above or below on this wall
						// draw the full wall; does this mean there are no windows on this exterior wall?
						bdraw.add_section(*this, 1, part, tex, color, dim_mask, 0, 0, 1, clip_windows, door_ztop, 0, offset_scale, 0, clamp_cube); // no_ao=1
						continue;
					}
					assert(!(wall_edges.size() & 1)); // must be an even number
					sort(wall_edges.begin(), wall_edges.end());

					for (unsigned e = 0; e < wall_edges.size(); e += 2) { // each pair of points should be the {left, right} edge of a wall section
						bool const first(e == 0), last(e+2 == wall_edges.size());
						cube_t c(part);
						c.d[!dim][0] = (first ? side_lo : (window_spacing*ceil ((wall_edges[e  ] - side_lo)/window_spacing) + side_lo)); // lo, clamped to whole windows
						c.d[!dim][1] = (last  ? side_hi : (window_spacing*floor((wall_edges[e+1] - side_lo)/window_spacing) + side_lo)); // hi, clamped to whole windows
						float const wall_len(c.get_sz_dim(!dim));
						if (wall_len < 0.5*window_spacing) continue; // wall too small to add here
						c.z2() = door_ztop;
						assert(c.is_strictly_normalized());
						tid_nm_pair_t tex2(tex);
						tex2.tscale_x = 0.5f*round_fp(wall_len/window_spacing)/wall_len;
						tex2.txoff    = -2.0*tex2.tscale_x*c.d[!dim][0];
						bdraw.add_section(*this, 1, c, tex2, color, dim_mask, 0, 0, 1, clip_windows, door_ztop, 0, offset_scale, 0, clamp_cube); // no_ao=1
					} // for e
				} // for dir
			} // for dim
		} // for f
	} // for i (parts)
	// draw attic windows
	bdraw.temp_tquads.clear();
	get_attic_windows(bdraw.temp_tquads, offset_scale);
	if (bdraw.temp_tquads.empty()) {has_attic_window = 0;}
	for (tquad_with_ix_t const &window : bdraw.temp_tquads) {bdraw.add_tquad(*this, window, bcube, tex, color);}
	// if camera is inside this building, cut out holes so that the exterior doors show through
	if (cut_door_holes) {cut_holes_for_ext_doors(bdraw, only_cont_pt, draw_parts_mask);}
}

void building_t::get_all_drawn_window_verts_as_quads(vect_vnctcc_t &verts) const { // for interior drawing
	building_draw_t bdraw; // should this be a static variable?
	get_all_drawn_window_verts(bdraw, 0, 1.0, nullptr, 1, 1); // lights_pass=0, no_skylights=1, draw_int_windows=1
	bdraw.get_all_mat_verts(verts, 0); // triangles=0; combine quad verts across materials (should only be one)
	assert((verts.size() & 3) == 0); // must be a multiple of 4
}

void building_t::cut_holes_for_ext_doors(building_draw_t &bdraw, point const &contain_pt, unsigned draw_parts_mask) const {
	if (doors.empty()) return;
	float const floor_spacing(get_window_vspace());
	vector3d const xlate(get_camera_coord_space_xlate());
	auto const parts_end(get_real_parts_end_inc_sec());

	for (auto d = doors.begin(); d != doors.end(); ++d) { // cut a hole for each door
		if (!camera_pdu.cube_visible(d->get_bcube() + xlate)) continue; // VFC
		tquad_with_ix_t door(*d);
		move_door_to_other_side_of_wall(door, 0.3, 0); // move a bit in front of the normal interior door (0.3 vs. 0.2)
		cube_t const door_bcube(door.get_bcube());
		bool contained(0);

		for (auto i = parts.begin(); i != parts_end; ++i) {
			if (!i->intersects(door_bcube)) continue;
			contained = ((draw_parts_mask & (1<<(i-parts.begin()))) != 0);
			if (contain_pt.z > (door_bcube.z1() + floor_spacing) && !i->contains_pt(contain_pt)) {contained = 0;} // camera in a different part on a floor above the door
			break;
		}
		if (!contained) continue; // part skipped, skip door as well

		if (draw_parts_mask == 0xFFFF) { // windowless case - check for exterior walls blocking the door
			point const end_pt(door_bcube.get_cube_center());

			for (auto i = parts.begin(); i != parts_end; ++i) {
				float tmin(0.0), tmax(1.0);
				if (!get_line_clip(contain_pt, end_pt, i->d, tmin, tmax)) continue; // no intersection
				float const t(tmax + 0.001); // slightly past the intersection
				if (t > 1.0) continue; // past the door, skip
				point const p(contain_pt + t*(end_pt - contain_pt));
				contained = 0;

				for (auto j = parts.begin(); j != parts_end; ++j) {
					if (j != i && j->contains_pt(p)) {contained = 1; break;} // inside the building
				}
				if (!contained) break; // outside the building
			} // for i
			if (!contained) continue;
		}
		clip_door_to_interior(door);
		bdraw.add_tquad(*this, door, bcube, tid_nm_pair_t(WHITE_TEX), WHITE);
	} // for d
}

bool building_t::get_nearby_ext_door_verts(building_draw_t &bdraw, shader_t &s, point const &pos, vector3d const &view_dir, float dist, bool update_state, bool only_open) {
	tquad_with_ix_t door;
	int const door_ix(find_ext_door_close_to_point(door, pos, dist));
	if (update_state) {register_open_ext_door_state(door_ix);}
	if (door_ix < 0) return 0; // no nearby door
	move_door_to_other_side_of_wall(door, -1.01, 0); // move a bit further away from the outside of the building to make it in front of the orig door
	clip_door_to_interior(door);
	bdraw.add_tquad(*this, door, bcube, tid_nm_pair_t(WHITE_TEX), WHITE);
	// draw the opened door
	building_draw_t door_draw;
	door_rotation_t drot; // return value is unused
	vector3d const normal(door.get_norm());
	bool const opens_outward(!is_house), dim(fabs(normal.x) < fabs(normal.y)), dir(normal[dim] < 0.0);
	add_door_verts(door.get_bcube(), door_draw, drot, door.type, dim, dir, 1.0, opens_outward, 1, 0); // open_amt=1.0, exterior=1, on_stairs=0
	// draw other exterior doors as closed in case they're visible through the open door; is this needed for pedestrians?
	if (!only_open) {get_ext_door_verts(door_draw, pos, view_dir, door_ix);}
	door_draw.draw(s, 0, 1); // direct_draw_no_vbo=1
	return 1;
}
void building_t::get_ext_door_verts(building_draw_t &bdraw, point const &viewer, vector3d const &view_dir, int skip_door_ix) const {
	for (auto d = doors.begin(); d != doors.end(); ++d) {
		if (int(d - doors.begin()) == skip_door_ix) continue; // skip this door
		vector3d const normal2(d->get_norm());
		if (dot_product_ptv(normal2, viewer, d->pts[0]) > 0.0) continue; // facing exterior of door rather than interior, skip
		if (view_dir != zero_vector && dot_product(view_dir, normal2) < 0.0) continue; // not visible
		tquad_with_ix_t door_rev(*d);
		std::swap(door_rev.pts[0], door_rev.pts[1]); // reverse winding order
		std::swap(door_rev.pts[2], door_rev.pts[3]);
		draw_building_ext_door(bdraw, door_rev, *this);
	}
}
bool building_t::get_all_nearby_ext_door_verts(building_draw_t &bdraw, shader_t &s, vector<point> const &pts, float dist) { // for pedestrians
	for (auto const &p : pts) {
		// we currently only support drawing one open door, so stop when we find one; future work is to use a bit mask to keep track of which doors are open
		if (get_nearby_ext_door_verts(bdraw, s, p, zero_vector, dist, 0, 1)) return 1; // no view_dir, update_state=0, only_open=1
	}
	return 0;
}

void building_t::get_split_int_window_wall_verts(building_draw_t &bdraw_front, building_draw_t &bdraw_back, point const &only_cont_pt_in, bool make_all_front) const {
	if (!is_valid()) return; // invalid building
	point const only_cont_pt(get_inv_rot_pos(only_cont_pt_in));
	building_mat_t const &mat(get_material());
	cube_t const cont_part(get_part_containing_pt(only_cont_pt)); // part containing the point
	// complex floorplan buildings can have odd exterior wall geometry where this splitting approach doesn't work well,
	// but if the building is windowless, then we can at least make the walls all front so that exterior doors are drawn properly
	if (has_complex_floorplan && !has_int_windows()) {make_all_front = 1;}
	
	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // multiple cubes/parts/levels; include house garage/shed
		if (is_basement(i)) continue; // skip basement walls because they have no windows
		
		if (building_has_open_ext_door) { // skip drawing wall in front of door if the camera is within NEAR_CLIP of it
			vector3d const offset(NEAR_CLIP*vector3d(cview_dir.x, cview_dir.y, 0.0));
			if (i->contains_pt_xy(only_cont_pt) && !i->contains_pt_xy(only_cont_pt + offset)) continue;
		}
		if (make_all_front || *i == cont_part || i->contains_pt(only_cont_pt) || // part containing the point
			are_parts_stacked(*i, cont_part)) // stacked building parts, contained, draw as front in case player can see through stairs
		{
			bdraw_front.add_section(*this, 1, *i, mat.wall_tex, wall_color, 3, 0, 0, 1, 0); // XY
			continue;
		}
		unsigned back_dim_mask(3), front_dim_mask(0); // enable dims: 1=x, 2=y, 4=z | disable cube faces: 8=x1, 16=x2, 32=y1, 64=y2, 128=z1, 256=z2
		cube_t front_clip_cube(*i);

		for (unsigned d = 0; d < 2; ++d) {
			if (i->d[ d][0] != cont_part.d[ d][1] && i->d[ d][1] != cont_part.d[ d][0]) continue; // not adj in dim d
			if (i->d[!d][0] >= cont_part.d[!d][1] || i->d[!d][1] <= cont_part.d[!d][0]) continue; // no overlap in dim !d
			// if we get here, *i and cont_part are adjacent in dim d
			if (i->d[!d][1] < only_cont_pt[!d] || i->d[!d][0] > only_cont_pt[!d]) break; // point not contained in other dim range, draw part as back
			
			for (unsigned e = 0; e < 2; ++e) { // check for coplanar sides (wall extensions)
				unsigned const disable_bit(1 << (3 + 2*(1-d) + e));
				if (i->d[!d][e] != cont_part.d[!d][e] && ((i->d[!d][e] < cont_part.d[!d][e]) ^ e)) {front_dim_mask |= disable_bit; continue;} // not coplanar, disable edge from front
				front_dim_mask |= (1<<(1-d)); // coplanar, make other edge dim a front dim
				back_dim_mask  |= disable_bit; // disable this edge from back
			}
			for (unsigned e = 0; e < 2; ++e) { // check for extensions outside cont_part where back walls could be viewed through a window and split them out
				if (i->d[!d][e] != cont_part.d[!d][e] && ((i->d[!d][e] < cont_part.d[!d][e]) ^ e)) {
					cube_t back_clip_cube(*i);
					front_clip_cube.d[!d][e] = back_clip_cube.d[!d][!e] = cont_part.d[!d][e]; // split point
					bdraw_back.add_section(*this, 1, *i, mat.wall_tex, wall_color, back_dim_mask, 0, 0, 1, 0, 0.0, 0, 1.0, 0, &back_clip_cube);
				}
			}
			back_dim_mask &= ~(1<<d); front_dim_mask |= (1<<d); // draw only the other dim as back and this dim as front
			break;
		} // for d
		if (back_dim_mask  > 0) {bdraw_back .add_section(*this, 1, *i, mat.wall_tex, wall_color, back_dim_mask,  0, 0, 1, 0);}
		if (front_dim_mask > 0) {bdraw_front.add_section(*this, 1, *i, mat.wall_tex, wall_color, front_dim_mask, 0, 0, 1, 0, 0.0, 0, 1.0, 0, &front_clip_cube);}
	} // for i
}

void building_t::get_ext_wall_verts_no_sec(building_draw_t &bdraw) const { // used for blocking room shadows between parts
	if (real_num_parts == 1) return; // one part, light can't leak
	float const clip_cube_dist_thresh(2.0*get_wall_thickness());
	building_mat_t const &mat(get_material());

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		if (p->z1() < ground_floor_z1) continue; // not needed for basement and extended basement
		unsigned const part_ix(p - parts.begin());
		unsigned dim_mask(3); // start with XY only
		bool draw_any(0);

		for (unsigned d = 0; d < 4; ++d) { // 4 sides of this part
			bool const dim(d >> 1), dir(d & 1);
			float const side_pos(p->d[dim][dir]);
			bool skip_this_side(side_pos == bcube.d[dim][dir]); // exterior wall is on the edge of the bcube and can't shadow anything
			// skip drawing of walls on sides that are already clipped by the light; optimization, and helps with drawing of inner faces of walkway exterior doors
			if (!smap_light_clip_cube.is_all_zeros()) {skip_this_side |= (fabs(side_pos - smap_light_clip_cube.d[dim][dir]) < clip_cube_dist_thresh);}
			// houses and building courtyards can have exterior doors not along the bcube that aren't handled by the above case;
			// drawing the wall containing this door in the shadow map will cause lighting artifacts, so skip this wall;
			// it should be okay because we only need one of two walls intersecting an interior->exterior->interior light ray to suppress it,
			// and buildings generally won't have two doors on adjacent interior sides
			if (part_ix < 4) {skip_this_side |= bool(door_sides[part_ix] & (1<<d));} // only check base parts
			if (skip_this_side) {dim_mask |= (1<<(d+3));} // disable cube faces: 8=x1, 16=x2, 32=y1, 64=y2
			else {draw_any = 1;}
		} // for d
		if (!draw_any) continue; // nothing to draw (optimization)
		// Note: this can cause shadows over walkway doors for buildings with walkways connecting to recessed part edges, which is rare
		bdraw.add_section(*this, 1, *p, mat.side_tex, side_color, dim_mask, 0, 0, 1, 0); // Note: ignores windows and door cutouts
	} // for p
}

void building_t::get_walkway_end_verts(building_draw_t &bdraw, point const &pos) const {
	float const room_exp(2.0*get_window_vspace());

	for (building_walkway_t const &w : walkways) {
		if (w.bcube.contains_pt(pos)) return; // light inside walkway - end shadow not needed
		if (!w.bcube_inc_rooms.contains_pt_exp_xy_only(pos, room_exp)) continue; // expand to include nearby lights
		bool const dir(w.bcube.get_center_dim(w.dim) < pos[w.dim]);
		float const wall_thickness(get_wall_thickness());
		cube_t wall_cube(w.bcube);
		wall_cube.d[w.dim][0] = wall_cube.d[w.dim][1] = w.bcube.d[w.dim][dir] - (dir ? 1.0 : -1.0)*0.5*wall_thickness; // shrink to zero width near the wall
		cube_t wall_cube_exp(wall_cube);
		wall_cube_exp.expand_in_dim(w.dim, wall_thickness);
		static vect_cube_t cubes;
		cubes.clear();
		cubes.push_back(wall_cube);
		tid_nm_pair_t tp; // untextured
		
		for (unsigned b = 0; b < 2; ++b) { // check doors for both buildings
			if (b && w.conn_bldg == nullptr) continue; // no connected building
			if (!(b ? w.conn_bldg->bcube : bcube).intersects(wall_cube_exp)) continue; // wrong building

			for (tquad_with_ix_t const &door : (b ? w.conn_bldg->doors : doors)) { // check for open door
				cube_t const door_bc(door.get_bcube());
				if (door_bc.z2() < pos.z || door_bc.z1() > pos.z) continue; // wrong floor
				if (!door_bc.intersects(wall_cube_exp)) continue; // wrong door
				if (!door_bc.contains_pt_exp(pre_smap_player_pos, get_door_open_dist())) continue; // skip closed door
				cube_t door_bc_exp(door_bc);
				door_bc_exp.expand_in_dim(w.dim, wall_thickness);
				swap_cube_dims(wall_cube,   w.dim, 2); // swap so that subtract can be done in the XY plane
				swap_cube_dims(door_bc_exp, w.dim, 2);
				cubes.clear();
				subtract_cube_from_cube(wall_cube, door_bc_exp, cubes);
				for (cube_t &c : cubes) {swap_cube_dims(c, w.dim, 2);} // swap back
				// draw open side doors
				float const door_width(door_bc.get_sz_dim(!w.dim)), door_hwidth(0.5*door_width);
				cube_t door_side(door_bc);
				door_side.d[w.dim][!dir] -= (dir ? 1.0 : -1.0)*door_hwidth; // extend into walkway

				for (unsigned d = 0; d < 2; ++d) { // left, right doors
					set_wall_width(door_side, door_bc.d[!w.dim][d], 0.8*wall_thickness, !w.dim);
					bdraw.add_cube(*this, door_side, tp, BLACK, 0, 3); // draw all sides
				}
			} // for door
		} // for b
		for (cube_t const &c : cubes) {bdraw.add_cube(*this, c, tp, BLACK, 0, (1 << unsigned(w.dim)));} // draw ends
		bdraw.add_cube(*this, w.bcube, tp, BLACK, 0, (1 << unsigned(!w.dim))); // draw sides of walkway to prevent light leaks from diagonal/adjacent rooms
	} // for w
}

void building_t::add_split_roof_shadow_quads(building_draw_t &bdraw) const {
	if (!interior || is_house || real_num_parts == 1) return; // no a stacked case
	float const light_zval(get_camera_pos().z); // light pos is stored in camera_pos during the shadow pass

	for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
		if (is_basement(i))        continue; // skip the basement
		if (i->z2() == bcube.z2()) continue; // skip top roof
		if (i->z2() > light_zval ) continue; // if roof is above the light, and lights all point down, then it can't cast a shadow

		if (clip_part_ceiling_for_stairs(*i, bdraw.temp_cubes, bdraw.temp_cubes2)) {
			for (auto c = bdraw.temp_cubes.begin(); c != bdraw.temp_cubes.end(); ++c) { // add floors after removing stairwells
				bdraw.add_section(*this, 1, *c, tid_nm_pair_t(), BLACK, 4, 1, 0, 1, 0); // only Z dim
			}
		}
	} // for i
}

// writes to the depth buffer only to prevent the terrain from being drawn in the basement
// to be called when the player is inside this building; when outside the building, the exterior walls/windows will write to the depth buffer instead
void draw_basement_entrance_cap(cube_t const &c, float z) {
	// only draw top surface - bottom surface of terrain is not drawn when player is in the basement
	vert_wrap_t const verts[4] = {point(c.x1(), c.y1(), z), point(c.x2(), c.y1(), z), point(c.x2(), c.y2(), z), point(c.x1(), c.y2(), z)};
	draw_verts(verts, 4, GL_TRIANGLE_FAN); // single quad
}
void building_t::write_basement_entrance_depth_pass(shader_t &s) const {
	if (!interior || !has_basement()) return;
	float const zval(get_basement().z2()), camera_z(get_camera_pos().z);
	if (camera_z < zval) return; // below upper basement level
	if (!is_house && camera_z > ground_floor_z1 + 2.0*get_window_vspace()) return; // floor 3+ of office building (can be visible through house L-shaped stairs)
	float const z(zval + BASEMENT_ENTRANCE_SCALE*get_floor_thickness()); // offset is required to prevent Z-fighting
	bool const depth_clamp_enabled(glIsEnabled(GL_DEPTH_CLAMP));
	s.set_cur_color(ALPHA0); // fully transparent
	select_texture(WHITE_TEX);
	enable_blend();
	glEnable(GL_CULL_FACE);
	if (!depth_clamp_enabled) {glEnable(GL_DEPTH_CLAMP);}

	for (auto i = interior->stairwells.begin(); i != interior->stairwells.end(); ++i) {
		if (i->z1() < zval && !i->in_ext_basement) {draw_basement_entrance_cap(*i, z);} // draw if this is a basement stairwell (not extended basement stairs)
	}
	if (has_pg_ramp() && !interior->ignore_ramp_placement) { // add opening for the ramp onto the ground floor
		draw_basement_entrance_cap(interior->pg_ramp, z);
	}
	if (!depth_clamp_enabled) {glDisable(GL_DEPTH_CLAMP);}
	glDisable(GL_CULL_FACE);
	disable_blend();
}


class building_creator_t {

	bool use_smap_this_frame=0, has_interior_geom=0, is_city=0, vbos_created=0;
	unsigned grid_sz=1, gpu_mem_usage=0;
	vector3d range_sz, range_sz_inv, max_extent;
	cube_t range, buildings_bcube;
	rand_gen_t rgen, ai_rgen;
	vect_building_t buildings;
	vect_bldg_walkway_t all_walkways; // walkways connecting city buildings
	vector<vector<unsigned>> bix_by_plot; // cached for use with pedestrian collisions
	// dynamic verts, static exterior verts, windows, window lights, interior walls/ceilings/floors, interior exterior walls
	building_draw_t building_draw, building_draw_vbo, building_draw_windows, building_draw_wind_lights, building_draw_interior, building_draw_int_ext_walls;
	point_sprite_drawer_sized building_lights;
	vector<point> points; // reused temporary

	struct grid_elem_t {
		vector<cube_with_ix_t> bc_ixs;
		vect_cube_t road_segs; // or driveways, porches, etc.
		cube_t bcube, extb_bcube; // base building, extended basement
		bool has_room_geom=0;

		bool empty() const {return (bc_ixs.empty() && road_segs.empty());}

		void add_bcube(cube_t const &c, unsigned ix, bool is_road_seg=0) {
			bcube.assign_or_union_with_cube(c);
			if (is_road_seg) {road_segs.push_back(c);} else {bc_ixs.emplace_back(c, ix);}
		}
		void add_building(building_t const &b, unsigned ix) {
			add_bcube(b.bcube, ix);
			if (b.has_ext_basement()) {extb_bcube.assign_or_union_with_cube(b.interior->basement_ext_bcube);}

			if (DRAW_WALKWAY_INTERIORS) {
				for (building_walkway_t const &w : b.walkways) {
					if (w.is_owner) {bcube.assign_or_union_with_cube(w.bcube);}
				}
			}
		}
		cube_t const &get_vis_bcube() const {return (player_in_ext_basement() ? extb_bcube : bcube);}
	};
	vector<grid_elem_t> grid, grid_by_tile;

	grid_elem_t &get_grid_elem(unsigned gx, unsigned gy) {
		assert(gx < grid_sz && gy < grid_sz && !grid.empty());
		return grid[gy*grid_sz + gx];
	}
	grid_elem_t const &get_grid_elem(unsigned gx, unsigned gy) const {
		assert(gx < grid_sz && gy < grid_sz && !grid.empty());
		return grid[gy*grid_sz + gx];
	}
	struct bix_by_x1 {
		vector<building_t> const &buildings;
		bix_by_x1(vector<building_t> const &buildings_) : buildings(buildings_) {}
		bool operator()(unsigned const a, unsigned const b) const {return (buildings[a].bcube.x1() < buildings[b].bcube.x1());}
	};
	unsigned get_grid_ix(point pos) const {
		range.clamp_pt_xy(pos);
		unsigned gxy[2] = {};
		for (unsigned d = 0; d < 2; ++d) {
			float const v((pos[d] - range.d[d][0])*range_sz_inv[d]);
			gxy[d] = unsigned(v*(grid_sz-1));
			assert(gxy[d] < grid_sz);
		}
		return (gxy[1]*grid_sz + gxy[0]);
	}
	void get_grid_range(cube_t const &bcube, unsigned ixr[2][2], bool expand_by_one=0) const { // {lo,hi}x{x,y}
		point llc(bcube.get_llc()), urc(bcube.get_urc());
		range.clamp_pt_xy(llc);
		range.clamp_pt_xy(urc);
		for (unsigned d = 0; d < 2; ++d) {
			float const v1((llc[d] - range.d[d][0])*range_sz_inv[d]), v2((urc[d] - range.d[d][0])*range_sz_inv[d]);
			ixr[0][d] = unsigned(v1*(grid_sz-1));
			ixr[1][d] = unsigned(v2*(grid_sz-1));
			assert(ixr[0][d] < grid_sz && ixr[1][d] < grid_sz);
			if (expand_by_one && ixr[0][d]   > 0      ) {--ixr[0][d];}
			if (expand_by_one && ixr[1][d]+1 < grid_sz) {++ixr[1][d];}
		}
	}
	void add_to_grid(cube_t const &bcube, unsigned bix, bool is_road_seg) {
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {get_grid_elem(x, y).add_bcube(bcube, bix, is_road_seg);}
		}
	}
	bool check_for_overlaps(vector<unsigned> const &ixs, cube_t const &test_bc, building_t const &b, float expand_rel, float expand_abs) const {
		for (auto i = ixs.begin(); i != ixs.end(); ++i) {
			building_t const &ob(get_building(*i));
			if (test_bc.intersects_xy(ob.bcube) && ob.check_bcube_overlap_xy(b, expand_rel, expand_abs)) return 1;
		}
		return 0;
	}
	bool check_for_overlaps(vector<cube_with_ix_t> const &bc_ixs, cube_t const &test_bc, building_t const &b, float expand_rel, float expand_abs) const {
		for (auto i = bc_ixs.begin(); i != bc_ixs.end(); ++i) {
			if (test_bc.intersects_xy(*i) && get_building(i->ix).check_bcube_overlap_xy(b, expand_rel, expand_abs)) return 1;
		}
		return 0;
	}

	void add_building_to_grid(building_t const &b, unsigned gix, unsigned bix) {
		grid_by_tile[gix].add_building(b, bix);
		if (b.enable_driveway_coll() && !b.driveway.is_all_zeros()) {grid_by_tile[gix].add_bcube(b.driveway, bix, 1);} // add driveway as well, but not porch
	}
	void build_grid_by_tile(bool single_tile) {
		grid_by_tile.clear();

		if (single_tile || world_mode != WMODE_INF_TERRAIN) { // not used in this mode - add all buildings to the first tile
			grid_by_tile.resize(1);
			grid_by_tile[0].bc_ixs.reserve(buildings.size());

			for(unsigned bix = 0; bix < buildings.size(); ++bix) {
				building_t const &b(buildings[bix]);
				if (!b.bcube.is_all_zeros()) {add_building_to_grid(b, 0, bix);} // skip invalid buildings
			}
			return;
		}
		//timer_t timer("build_grid_by_tile");
		map<uint64_t, unsigned> tile_to_gbt;

		for(unsigned bix = 0; bix < buildings.size(); ++bix) {
			building_t const &b(buildings[bix]);
			if (b.bcube.is_all_zeros()) continue; // skip invalid buildings
			uint64_t const tile_id(get_tile_id_for_cube(b.bcube));
			auto it(tile_to_gbt.find(tile_id));
			unsigned gix;

			if (it == tile_to_gbt.end()) { // new element
				gix = grid_by_tile.size();
				grid_by_tile.push_back(grid_elem_t());
				tile_to_gbt[tile_id] = gix;
			}
			else { // existing element
				gix = it->second;
				assert(gix < grid_by_tile.size());
			}
			add_building_to_grid(b, gix, bix);
		} // for bix
	}

	bool check_valid_building_placement(building_params_t const &params, building_t const &b, vect_cube_t const &avoid_bcubes, cube_t const &avoid_bcubes_bcube,
		float min_building_spacing, unsigned plot_ix, bool non_city_only, bool use_city_plots, bool check_plot_coll)
	{
		float const expand_val(b.is_rotated() ? 0.05 : 0.1); // expand by 5-10% (relative - multiplied by building size)
		vector3d const b_sz(b.bcube.get_size());
		vector3d expand(expand_val*b_sz);
		for (unsigned d = 0; d < 2; ++d) {max_eq(expand[d], min_building_spacing);} // ensure the min building spacing (only applies to the current building)
		cube_t test_bc(b.bcube);
		test_bc.expand_by_xy(expand);

		if (use_city_plots) { // use city blocks
			assert(plot_ix < bix_by_plot.size());
			if (check_for_overlaps(bix_by_plot[plot_ix], test_bc, b, expand_val, min_building_spacing)) return 0;
			bix_by_plot[plot_ix].push_back(buildings.size());
		}
		else if (check_plot_coll && !avoid_bcubes.empty() && avoid_bcubes_bcube.intersects_xy(test_bc) &&
			has_bcube_int_xy(test_bc, avoid_bcubes, params.sec_extra_spacing)) // extra expand val
		{
			return 0;
		}
		else if (non_city_only && check_city_tline_cube_intersect_xy(test_bc)) {return 0;} // check transmission lines
		else {
			float const extra_spacing(non_city_only ? params.sec_extra_spacing : 0.0); // absolute value of expand
			test_bc.expand_by_xy(extra_spacing);
			unsigned ixr[2][2];
			get_grid_range(test_bc, ixr);

			for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
				for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
					grid_elem_t const &ge(get_grid_elem(x, y));
					if (!test_bc.intersects_xy(ge.bcube)) continue;
					if (check_for_overlaps(ge.bc_ixs, test_bc, b, expand_val, max(min_building_spacing, extra_spacing))) {return 0;}
				} // for x
			} // for y
		}
		return 1;
	}
	bool own_this_building(building_t const *building) const {return (building >= buildings.data() && building < buildings.data()+buildings.size());}

	struct building_cand_t : public building_t {
		vect_cube_t &temp_parts;
		building_cand_t(vect_cube_t &temp_parts_) : temp_parts(temp_parts_) {temp_parts.clear(); parts.swap(temp_parts);} // parts takes temp_parts memory
		~building_cand_t() {parts.swap(temp_parts);} // memory returned from parts to temp_parts
	};

public:
	building_creator_t(bool is_city_=0) : is_city(is_city_), max_extent(zero_vector), building_draw(is_city), building_draw_vbo(is_city) {}
	bool empty() const {return buildings.empty();}
	bool get_is_city() const {return is_city;}
	bool has_interior_to_draw() const {return (has_interior_geom && !building_draw_interior.empty());}

	void clear() {
		buildings.clear();
		grid.clear();
		grid_by_tile.clear();
		bix_by_plot.clear();
		clear_vbos();
		buildings_bcube = cube_t();
		gpu_mem_usage = 0;
	}
	unsigned get_num_buildings() const {return buildings.size();}
	unsigned get_gpu_mem_usage() const {return gpu_mem_usage;}
	vector3d const &get_max_extent() const {return max_extent;}
	building_t const &get_building(unsigned ix) const {assert(ix < buildings.size()); return buildings[ix];}
	building_t       &get_building(unsigned ix)       {assert(ix < buildings.size()); return buildings[ix];} // non-const version; not intended to be used to change geometry
	cube_t const &get_building_bcube(unsigned ix) const {return get_building(ix).bcube;}
	
	bool get_building_door_pos_closest_to(unsigned ix, point const &target_pos, point &door_pos, bool inc_garage_door) const {
		return get_building(ix).get_building_door_pos_closest_to(target_pos, door_pos, inc_garage_door);
	}
	cube_t register_deck_and_get_part_bounds(unsigned ix, cube_t const &deck) {return get_building(ix).register_deck_and_get_part_bounds(deck);}
	cube_t const &get_bcube() const {return buildings_bcube;}
	bool is_visible(vector3d const &xlate) const {return (!empty() && camera_pdu.cube_visible(buildings_bcube + xlate));}
	bool is_single_tile() const {return (grid_by_tile.size() == 1);}
	
	bool get_building_hit_color(point const &p1, point const &p2, colorRGBA &color) const { // exterior only; p1 and p2 are in building space
		if (p1.x == p2.x && p1.y == p2.y && player_in_basement >= 3 && player_building != nullptr) {
			if (player_building->get_interior_color_at_xy(p1, color)) return 1; // handle interior of extended basement the player is in
		}
		float t(1.0); // unused
		unsigned hit_bix(0);
		unsigned const ret(check_line_coll(p1, p2, t, hit_bix, 0, 1, 1)); // ret_any_pt=0, no_coll_pt=1, check_non_coll=1; returns type of surface that was hit
		if (ret == BLDG_COLL_NONE) return 0;
		building_t const &b(get_building(hit_bix));

		if (p1.x == p2.x && p1.y == p2.y && b.get_interior_color_at_xy(p1, color)) {
			return 1; // vertical line (from map mode) hit the roof of a building with an interior
		}
		switch (ret) {
		case BLDG_COLL_SIDE    : color = b.get_avg_side_color  (); break;
		case BLDG_COLL_ROOF    : color = b.get_avg_roof_color  (); break;
		case BLDG_COLL_DRIVEWAY: color = LT_GRAY ; break;
		case BLDG_COLL_FENCE   : color = LT_BROWN; break;
		case BLDG_COLL_SKYLIGHT: color = LT_BLUE ; break;
		case BLDG_COLL_DETAIL  :
			color = b.get_avg_detail_color();
			if (color == b.get_avg_roof_color()) {color *= 1.5;} // lighten if the same color as the roof so that details stand out
			break;
		default: assert(0);
		}
		return 1;
	}

	struct city_prob_t {
		bool inited=0, same_mat_per_block[2]={0,0}, same_size_per_block[2]={0,0}, same_geom_per_mat[2]={0,0}, same_houses_citywide[2]={0,0}; // {office, house}

		void init(building_params_t const &params, rand_gen_t &rgen) {
			if (inited) return;
			same_mat_per_block  [1] = rgen.rand_probability(params.house_same_mat_prob );
			same_size_per_block [1] = rgen.rand_probability(params.house_same_size_prob);
			same_geom_per_mat   [1] = rgen.rand_probability(params.house_same_geom_prob);
			same_houses_citywide[1] = rgen.rand_probability(params.house_same_per_city_prob);
			same_mat_per_block  [0] = rgen.rand_probability(params.office_same_mat_prob );
			same_size_per_block [0] = rgen.rand_probability(params.office_same_size_prob);
			same_geom_per_mat   [0] = rgen.rand_probability(params.office_same_geom_prob);
			same_houses_citywide[0] = rgen.rand_probability(params.office_same_per_city_prob);
			inited = 1;
		}
	};
	struct vect_city_prob_t {
		vector<city_prob_t> cps;
		vector<unsigned> city_for_building;
		city_prob_t def_prob;
		bool enabled;

		vect_city_prob_t(bool enabled_) : enabled(enabled_) {}

		city_prob_t const &get(building_params_t const &params, rand_gen_t &rgen, unsigned city_ix) {
			if (!enabled) return def_prob;
			if (city_ix >= cps.size()) {cps.resize(city_ix+1);} // allocate a new entry
			city_prob_t &ret(cps[city_ix]);
			ret.init(params, rgen);
			return ret;
		}
		city_prob_t const &get(unsigned bix) {
			if (!enabled) return def_prob;
			assert(bix < city_for_building.size());
			assert(city_for_building[bix] < cps.size());
			return cps[city_for_building[bix]];
		}
		void next_building(unsigned city_ix) {
			if (enabled) {city_for_building.push_back(city_ix);} // only if in use
		}
	};

	void gen(building_params_t const &params, bool city_only, bool non_city_only, bool is_tile, bool allow_flatten, int rseed=123) {
		assert(!(city_only && non_city_only));
		clear();
		if (params.tt_only && world_mode != WMODE_INF_TERRAIN)    return;
		if (params.gen_inf_buildings() && !is_tile && !city_only) return; // secondary buildings - not added here
		vect_city_zone_t city_plot_bcubes;
		vector<unsigned> valid_city_plot_ixs;
		bool maybe_residential(0); // Note: may not be correct for a mix of residential and non-residential, but this only matters if the material nonemptiness or ranges differ

		if (city_only) {
			unsigned num_residential(0), num_non_residential(0);
			get_city_plot_zones(city_plot_bcubes); // Note: assumes approx equal area for placement distribution

			for (auto i = city_plot_bcubes.begin(); i != city_plot_bcubes.end(); ++i) {
				if (i->is_park) continue; // skip parks
				valid_city_plot_ixs.push_back(i - city_plot_bcubes.begin()); // record non-park plots
				if (i->is_residential) {++num_residential;} else {++num_non_residential;}
			}
			//assert(!valid_city_plot_ixs.empty()); // too strong? what happens if this is true?
			maybe_residential = (num_residential > num_non_residential); // consider this city residential if the majority of non-park plots are residential
		}
		vector<unsigned> const &mat_ix_list(params.get_mat_list(city_only, non_city_only, maybe_residential));
		if (params.materials.empty() || mat_ix_list.empty()) return; // no materials
		timer_t timer("Gen Buildings", !is_tile);
		float const def_water_level(get_water_z_height()), min_building_spacing(get_min_obj_spacing());
		vector3d const offset(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0);
		vector3d const xlate((world_mode == WMODE_INF_TERRAIN) ? offset : zero_vector); // cancel out xoff2/yoff2 translate
		vector3d const delta_range((world_mode == WMODE_INF_TERRAIN) ? zero_vector : offset);
		range = params.materials[mat_ix_list.front()].pos_range; // range is union over all material ranges
		for (auto i = mat_ix_list.begin()+1; i != mat_ix_list.end(); ++i) {range.union_with_cube(params.materials[*i].pos_range);}
		range     += delta_range;
		range_sz   = range.get_size(); // Note: place_radius is relative to range cube center
		max_extent = zero_vector;
		assert(range_sz.x > 0.0 && range_sz.y > 0.0);
		UNROLL_2X(range_sz_inv[i_] = 1.0/range_sz[i_];) // xy only
		if (!is_tile) {buildings.reserve(params.num_place);}
		grid_sz = (is_tile ? 4 : 32); // tiles are small enough that they don't need grids
		grid.resize(grid_sz*grid_sz); // square
		unsigned num_tries(0), num_gen(0), num_skip(0);
		if (rseed == 0) {rseed = 123;} // 0 is a bad value
		rseed += params.buildings_rand_seed; // add in rand_seed from the config file
		rgen.set_state(rand_gen_index, rseed); // update when mesh changes, otherwise determinstic
		vect_cube_t avoid_bcubes;
		cube_t avoid_bcubes_bcube;

		if (non_city_only) {
			get_city_bcubes(avoid_bcubes);
			get_city_road_bcubes(avoid_bcubes, 1); // connector roads only
			get_all_model_bcubes(avoid_bcubes);
			expand_cubes_by_xy(avoid_bcubes, get_road_max_width());
			avoid_bcubes_bcube = get_bcubes_union(avoid_bcubes);
		}
		bool const use_city_plots(!valid_city_plot_ixs.empty()), check_plot_coll(!avoid_bcubes.empty());
		vect_city_prob_t city_prob(use_city_plots); // calculated and reused once per city
		bix_by_plot.resize(city_plot_bcubes.size());
		point center(all_zeros);
		unsigned num_consec_fail(0), max_consec_fail(0);
		vect_cube_t temp_parts;

		for (unsigned i = 0; i < params.num_place; ++i) {
			bool success(0);

			for (unsigned n = 0; n < params.num_tries; ++n) { // 10 tries to find a non-overlapping building placement
				unsigned plot_ix(0), city_block_ix(0), pref_dir(0);
				int city_ix(-1);
				bool residential(0), no_residential(0);

				if (use_city_plots) { // select a random plot, if available
					bool success(0);

					for (unsigned N = 0; N < 10; ++N) { // 10 tries to choose a plot that has the capacity for a new building
						plot_ix = valid_city_plot_ixs[rgen.rand() % valid_city_plot_ixs.size()];
						assert(plot_ix < city_plot_bcubes.size());
						if (!city_plot_bcubes[plot_ix].is_full()) {success = 1; break;}
					}
					if (!success) break; // all candidate plots were full
					if (city_plot_bcubes[plot_ix].is_residential) {residential = 1;} else {no_residential = 1;}
					if (residential && params.mat_gen_ix_res.empty()) break; // no residential buildings available, break from n loop (but retry i loop with new plot)
				}
				cube_t pos_range;
				float border_scale(1.0);
				unsigned max_floors(0); // starts at unlimited
				building_cand_t b(temp_parts);
				
				if (use_city_plots) { // select a random plot, if available
					city_zone_t const &plot(city_plot_bcubes[plot_ix]);
					pos_range       = plot;
					center.z        = plot.zval; // optimization: take zval from plot rather than calling get_exact_zval()
					b.assigned_plot = plot; // only really needed for residential sub-plots
					b.address       = plot.address;
					city_block_ix   = ((plot.parent_plot_ix >= 0) ? plot.parent_plot_ix : plot_ix);
					city_ix         = plot.city_ix;
					max_floors      = plot.max_floors;
					if (residential) {pref_dir = plot.street_dir;}
					// force min spacing between building and edge of plot, but make sure the plot remains normalized after shrinking
					pos_range.expand_by_xy(-min(min_building_spacing, 0.45f*min(plot.dx(), plot.dy())));
					if (plot.capacity == 1) {border_scale *= 2.0;} // use smaller border scale since the individual building plots should handle borders
				}
				city_prob_t const &CP(city_prob.get(params, rgen, city_ix));
				rand_gen_t group_rgen;
				group_rgen.set_state(rseed, (CP.same_houses_citywide[residential] ? city_ix : city_block_ix)+1); // varies per city block
				group_rgen.rand_mix();
				rand_gen_t &rgen_mat(CP.same_mat_per_block [residential] ? group_rgen : rgen); // for material and color
				rand_gen_t &rgen_sz (CP.same_size_per_block[residential] ? group_rgen : rgen); // for size, height, and orient
				b.mat_ix = params.choose_rand_mat(rgen_mat, city_only, non_city_only, residential); // set material
				building_mat_t const &mat(b.get_material());
				if (!use_city_plots) {pos_range = mat.pos_range + delta_range;} // select pos range by material
				vector3d const pos_range_sz(pos_range.get_size());
				assert(pos_range_sz.x > 0.0 && pos_range_sz.y > 0.0);
				point const place_center(pos_range.get_cube_center());
				bool const is_residential_block(residential && !pref_dir);
				float const min_center_dist(is_residential_block ? 0.3*min(pos_range_sz.x, pos_range_sz.y) : 0.0);
				bool keep(0);
				++num_tries;

				for (unsigned m = 0; m < params.num_tries; ++m) {
					for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(pos_range.d[d][0], pos_range.d[d][1]);} // x,y
					// place residential buildings around the edges of the plot (pos_range) unless the street dir has already been assigned (individual house plot) / keep out of the center
					if (min_center_dist > 0.0 && dist_xy_less_than(center, place_center, min_center_dist)) continue;
					if (is_tile || mat.place_radius == 0.0 || dist_xy_less_than(center, place_center, mat.place_radius)) {keep = 1; break;} // place_radius ignored for tiles
				}
				if (!keep) continue; // placement failed, skip

				if (!no_residential) {
					b.is_house = (mat.house_prob > 0.0 && (residential || rgen.rand_probability(mat.house_prob))); // force a house if residential and houses are enabled
				}
				float const size_scale(b.is_house ? mat.gen_house_size_scale(rgen_sz) : 1.0);
				
				for (unsigned d = 0; d < 2; ++d) { // x,y
					float const size_cap(border_scale*pos_range_sz[d]*(b.is_house ? 0.8 : 1.0)); // size cap relative to plot size
					float const sz(0.5*rgen_sz.rand_uniform(min(size_scale*mat.sz_range.d[d][0], 0.3f*size_cap),
						                                    min(size_scale*mat.sz_range.d[d][1], 0.5f*size_cap))); // use pos range size for max
					b.bcube.d[d][0] = center[d] - sz;
					b.bcube.d[d][1] = center[d] + sz;
				}
				if (is_residential_block && b.is_house && b.bcube.contains_pt_xy(place_center)) continue; // house should not contain the center point of the plot
				if ((use_city_plots || is_tile) && !pos_range.contains_cube_xy(b.bcube)) continue; // not completely contained in plot/tile (pre-rot)
				if (!use_city_plots) {b.gen_rotation(rgen_sz);} // city plots are Manhattan (non-rotated) - must rotate before bcube checks below
				if (is_tile && !pos_range.contains_cube_xy(b.bcube)) continue; // not completely contained in tile
				if (start_in_inf_terrain && b.bcube.contains_pt_xy(get_camera_pos())) continue; // don't place a building over the player appearance spot
				if (!check_valid_building_placement(params, b, avoid_bcubes, avoid_bcubes_bcube, min_building_spacing,
					city_block_ix, non_city_only, use_city_plots, check_plot_coll)) continue; // check overlap (use city plot_ix rather than sub-plot ix)
				++num_gen;
				if (!use_city_plots) {center.z = get_exact_zval(center.x+xlate.x, center.y+xlate.y);} // only calculate when needed
				float const z_sea_level(center.z - def_water_level);
				if (z_sea_level < 0.0) break; // skip underwater buildings, failed placement
				if (z_sea_level < mat.min_alt || z_sea_level > mat.max_alt) break; // skip bad altitude buildings, failed placement
				float const hmin(use_city_plots ? pos_range.z1() : 0.0), hmax(use_city_plots ? pos_range.z2() : 1.0);
				assert(hmin <= hmax);
				float const height_range(mat.sz_range.dz());
				assert(height_range >= 0.0);
				float const z_size_scale(size_scale*(b.is_house ? rgen_sz.rand_uniform(0.6, 0.8) : 1.0)); // make houses slightly shorter on average to offset extra height added by roof
				float height_val(0.5f*z_size_scale*(mat.sz_range.z1() + height_range*rgen_sz.rand_uniform(hmin, hmax)));
				if (max_floors > 0) {min_eq(height_val, max_floors*b.get_window_vspace());} // limit height based on max floors
				assert(height_val > 0.0);
				b.set_z_range(center.z, (center.z + height_val));
				assert(b.bcube.is_strictly_normalized());
				mat.side_color.gen_color(b.side_color, rgen_mat);
				mat.roof_color.gen_color(b.roof_color, rgen_mat);
				if (use_city_plots) {b.street_dir = (pref_dir ? pref_dir : get_street_dir(b.bcube, pos_range));}
				if (city_only     ) {b.is_in_city = 1; b.city_ix = city_ix;}
				add_to_grid(b.bcube, buildings.size(), 0);
				vector3d const sz(b.bcube.get_size());
				float const mult[3] = {0.5, 0.5, 1.0}; // half in X,Y and full in Z
				UNROLL_3X(max_extent[i_] = max(max_extent[i_], mult[i_]*sz[i_]);)
				buildings.push_back(b);
				city_prob.next_building(city_ix);
				if (use_city_plots) {++city_plot_bcubes[plot_ix].nbuildings;}
				success = 1;
				break; // done
			} // for n
			if (success) {num_consec_fail = 0;}
			else {
				++num_consec_fail;
				max_eq(max_consec_fail, num_consec_fail);

				if (num_consec_fail >= (is_tile ? 50U : 5000U)) { // too many failures - give up
					if (!is_tile) {cout << "Failed to place a building after " << num_consec_fail << " tries, giving up after " << i << " iterations" << endl;}
					break;
				}
			}
		} // for i
		if (buildings.capacity() > 2*buildings.size()) {buildings.shrink_to_fit();}
		// after this point buildings should no longer be resized and their pointers can be used without worrying about invalidation, at least within this buildings block
		bix_by_x1 cmp_x1(buildings);
		for (auto i = bix_by_plot.begin(); i != bix_by_plot.end(); ++i) {sort(i->begin(), i->end(), cmp_x1);}
		if (!is_tile) {timer.end();} // use a single timer for tile mode
		parse_universe_name_str_tables(); // must do this here because it's not legal to call in MT code below

		if (params.flatten_mesh && !use_city_plots) { // not needed for city plots, which are already flat
			timer_t timer("Gen Building Zvals", !is_tile);
			bool const do_flatten(allow_flatten && using_tiled_terrain_hmap_tex()); // can't always flatten terrain when using tiles

#pragma omp parallel for schedule(static,1) if (!is_tile)
			for (int i = 0; i < (int)buildings.size(); ++i) {
				building_t &b(buildings[i]);

				if (do_flatten) { // flatten the mesh under the bcube to a height of mesh_zval
					//assert(!b.is_rotated()); // too strong?
					flatten_hmap_region(b.bcube);
				}
				else { // extend building bottom downward to min mesh height
					bool const shift_top(1); // shift is required to preserve height for floor alignment of building interiors
					float &zmin(b.bcube.z1()); // Note: grid bcube z0 value won't be correct, but will be fixed conservatively below
					float const orig_zmin(zmin);
					unsigned num_below(0);
					
					for (int d = 0; d < 4; ++d) {
						float const zval(get_exact_zval(b.bcube.d[0][d&1]+xlate.x, b.bcube.d[1][d>>1]+xlate.y)); // approximate for rotated buildings
						min_eq(zmin, zval);
						num_below += (zval < def_water_level);
					}
					max_eq(zmin, def_water_level); // don't go below the water
					float const dz(orig_zmin - zmin), max_dz(b.get_material().max_delta_z);
					if (shift_top) {b.bcube.z2() -= dz;} // shift top down as well to keep the height constant
					if (num_below > 2 || // more than 2 corners underwater
						(max_dz > 0.0 && dz > max_dz)) // too steep of a slope
					{
						b.bcube.set_to_zeros();
						++num_skip;
					}
					else if (!b.parts.empty()) {
						b.parts.back().z1() = b.bcube.z1(); // update base z1
						if (shift_top) {b.parts.back().z2() -= dz;} // shift top down as well
						assert(b.parts.back().dz() > 0.0);
					}
				}
			} // for i
			if (do_flatten) { // use conservative zmin for grid
				for (auto i = grid.begin(); i != grid.end(); ++i) {i->bcube.z1() = def_water_level;}
			}
		} // if flatten_mesh
		{ // open a scope
			has_office_chair_model(); // must call this to load models here, since it's called inside building_t::gen_geometry() and is not thread safe
			timer_t timer2("Gen Building Geometry", !is_tile); // 120ms/700ms => 160ms/900ms
			bool const gen_interiors(global_building_params.gen_building_interiors);
			bool const use_mt(!is_tile || gen_interiors); // only single threaded for tiles with no interiors, which is a fast case anyway

			// split into houses and office buildings run on 2 threads
#pragma omp parallel for schedule(static) num_threads(2) if (use_mt)
			for (int is_house=0; is_house < 2; ++is_house) {
				for (unsigned i = 0; i < buildings.size(); ++i) {
					building_t &b(buildings[i]);
					if (b.is_house != bool(is_house)) continue; // wrong pass
					unsigned const rs_ix(city_prob.get(i).same_geom_per_mat[b.is_house] ? b.mat_ix : i); // same material, maybe from same block/city; could also use city_ix
					b.gen_geometry(rs_ix, 1337*rs_ix+rseed);
				}
			} // for is_house
			if (city_only && gen_interiors && global_building_params.max_ext_basement_room_depth > 0) {try_join_city_building_ext_basements(buildings);}
		} // close the scope
		if (0 && non_city_only) { // perform room graph analysis
			timer_t timer3("Building Room Graph Analysis");
			for (auto b = buildings.begin(); b != buildings.end(); ++b) {
				if (b->has_complex_floorplan) continue; // room graph isn't really valid for this building type
				//if (b->is_house) continue;
				unsigned num_comp(b->count_connected_room_components());
				if (b->has_sec_bldg()) {--num_comp;} // exclude garage/shed
				//cout << num_comp;
				if (num_comp > 1) {cout << num_comp << ": " << b->bcube.get_cube_center().str() << endl;}
			}
			cout << endl;
		}
		// re-generate grid based on new building bcubes that include things like roofs and chimneys;
		// since bcubes should only increase in size, we don't need to reset grid bcubes
		for (auto g = grid.begin(); g != grid.end(); ++g) {g->bc_ixs.clear();}

		for (auto b = buildings.begin(); b != buildings.end(); ++b) { // add driveways, porches, etc. and calculate has_interior_geom
			unsigned const bix(b - buildings.begin());

			if (b->enable_driveway_coll()) {
				if (!b->driveway.is_all_zeros()) {
					cube_t driveway_ext(b->driveway);
					driveway_ext.expand_by_xy(0.2*b->get_window_vspace()); // expand so that grass is excluded at the edges; determined experimentally
					add_to_grid(driveway_ext, bix, 1);
				}
				if (!b->porch.is_all_zeros()) {add_to_grid(b->porch, bix, 1);}

				for (auto const &d : b->doors) { // handle steps for exterior doors
					if (d.type == tquad_with_ix_t::TYPE_GDOOR) continue; // already handled by driveway
					cube_t step(b->get_step_for_ext_door(d));
					step.translate_dim(2, -b->get_fc_thickness()); // shift down to make player coll smoother
					if (step.z1() > b->ground_floor_z1) continue; // not on the ground floor
					add_to_grid(step, bix, 1);
				}
			}
			add_to_grid(b->bcube, bix, 0);
			buildings_bcube.assign_or_union_with_cube(b->bcube);
			has_interior_geom |= b->has_interior();
		} // for b
		if (!is_tile && (!city_only || maybe_residential)) {place_building_trees(rgen);}

		if (!is_tile) {
			cout << "WM: " << world_mode << " MCF: " << max_consec_fail << " Buildings: " << params.num_place << " / " << num_tries << " / " << num_gen
				 << " / " << buildings.size() << " / " << (buildings.size() - num_skip) << endl;
			building_stats_t s;
			for (auto b = buildings.begin(); b != buildings.end(); ++b) {b->update_stats(s);}
			cout << TXT(s.nbuildings) << TXT(s.nparts) << TXT(s.ndetails) << TXT(s.ntquads) << TXT(s.ndoors) << TXT(s.ninterior)
				 << TXT(s.nrooms) << TXT(s.nceils) << TXT(s.nfloors) << TXT(s.nwalls) << TXT(s.nrgeom) << TXT(s.nobjs) << TXT(s.nverts) << endl;
		}
		if (city_only) { // connect with walkways here
			vect_cube_t city_bcubes;
			get_city_bcubes(city_bcubes);
			for (cube_t const &c : city_bcubes) {connect_buildings_with_walkways(c);}
		}
		build_grid_by_tile(is_tile);
		if (!city_only) {create_vbos(is_tile);} // city VBOs are created later, after skyways are added
	} // end gen()

	void place_building_trees(rand_gen_t &rgen) {
		if (!has_city_trees()) return;
		vector<point> placements;

		for (auto b = buildings.begin(); b != buildings.end(); ++b) {
			if (b->tree_pos != all_zeros) {placements.push_back(b->tree_pos);}
		}
		if (placements.empty()) return;
		sort(placements.begin(), placements.end(), [](point const &a, point const &b) {return (a.x < b.x);}); // sort by xval
		float const block_xsize(X_SCENE_SIZE);
		float cur_xval(0.0);
		unsigned num_blocks(0);
		bool const allow_bush = 0; // no bushes for now
		bool const is_sm_tree = 0; // deciduous trees only
		bool const add_bush   = 0;

		for (auto p = placements.begin(); p != placements.end(); ++p) {
			if (p == placements.begin() || p->x > (cur_xval + block_xsize)) {
				tree_placer.begin_block(is_sm_tree);
				cur_xval = p->x;
				++num_blocks;
			}
			int const ttype(rgen.rand()%100); // Note: okay to leave at -1; also, don't have to set to a valid tree type
			tree_placer.add(*p, 0, ttype, allow_bush, add_bush, is_sm_tree);
		} // for p
		cout << "Num Placed Trees: " << placements.size() << ", Blocks: " << num_blocks << endl;
	}

	void get_all_helipads(vect_cube_t &helipads) const {
		for (auto b = buildings.begin(); b != buildings.end(); ++b) {
			if (b->has_helipad) {helipads.push_back(b->get_helipad_bcube());}
		}
	}
	void add_building_signs(cube_t const &region_bcube, vector<sign_t> &signs) const {
		// Note: region_bcube is currently only used to select the buildings within a city, and there are a small number of cities,
		// so it should be okay to iterate rather than using a grid query
		for (building_t const &b : buildings) {
			if (!region_bcube.intersects_xy(b.bcube)) continue; // wrong region/city
			b.add_signs(signs);
		}
	}
	void add_building_flags(cube_t const &region_bcube, vector<city_flag_t> &flags) { // non-const because flags are cached in buildings
		for (building_t &b : buildings) { // same note as in add_building_signs
			if (!region_bcube.intersects_xy(b.bcube)) continue; // wrong region/city
			b.add_flags(flags);
		}
	}
	void update_ai_state(float delta_dir) { // called once per frame
		if (!global_building_params.building_people_enabled()) return;
		point const camera_bs(get_camera_building_space());
		float const dmax(1.5f*(X_SCENE_SIZE + Y_SCENE_SIZE));
		if (!get_bcube().closest_dist_less_than(camera_bs, dmax)) return; // too far away
		buildings.ai_room_update(delta_dir, dmax, camera_bs, ai_rgen);
	}

	static void select_person_shadow_shader(shader_t &person_shader) {
		if (!person_shader.is_setup()) {
			enable_animations_for_shader(person_shader);
			setup_smoke_shaders(person_shader, 0.0, 0, 0, 0, 0, 0, 0);
		} else {person_shader.make_current();}
	}

	static void multi_draw_shadow(vector3d const &xlate, vector<building_creator_t *> const &bcs) {
		DebugScope scope("building_multi_draw_shadow");
		//timer_t timer("Draw Buildings Shadow");
		push_scene_xlate(xlate); // drawn in building space
		shader_t s, amask_shader, person_shader;
		s.begin_shadow_map_shader();
		glEnable(GL_CULL_FACE); // slightly faster for interior shadow maps
		vector<point> points; // reused temporary
		static building_draw_t ext_parts_draw; // roof and exterior walls; reused across calls
		bool const sec_camera_mode(pre_smap_player_pos != actual_player_pos); // hack to determine if this is the shadow for a security camera light
		bool is_house(0);

		for (auto i = bcs.begin(); i != bcs.end(); ++i) {
			if (interior_shadow_maps) { // draw interior shadow maps
				occlusion_checker_noncity_t oc(**i);
				point const lpos(get_camera_pos() - xlate); // Note: camera_pos is actually the light pos
				bool found_building(0);

				// draw interior for the building containing the light
				for (auto g = (*i)->grid_by_tile.begin(); g != (*i)->grid_by_tile.end(); ++g) {
					if (!g->get_vis_bcube().contains_pt_xy(lpos)) continue; // wrong tile (note that z test is skipped to handle skylights)

					for (auto bi = g->bc_ixs.begin(); bi != g->bc_ixs.end(); ++bi) {
						building_t &b((*i)->get_building(bi->ix));
						if (!b.interior) continue; // no interior
						point lpos_clamped(lpos);
						// include skylight light sources, which are above the building; buildings can't stack vertically, so the light can't belong to a different building
						if (!b.skylights.empty()) {min_eq(lpos_clamped.z, b.bcube.z2());}
						bool const camera_in_walkway(b.check_pt_in_or_near_walkway(pre_smap_player_pos, 1, 0, 0)); // owned_only=1, inc_open_door=0, inc_conn_room=0
						if (!b.point_in_building_or_basement_bcube(lpos_clamped) && !camera_in_walkway) continue; // wrong building
						(*i)->building_draw_interior.draw_for_draw_range(s, b.interior->draw_range, 1); // shadow_only=1
						b.add_split_roof_shadow_quads(ext_parts_draw);
						// no batch draw for shadow pass since textures aren't used; draw everything, since shadow may be cached
						bool camera_in_this_building(b.check_point_or_cylin_contained(pre_smap_player_pos, 0.0, points, 1, 1, 0)); // inc_attic=1, inc_ext_basement=1, inc_roof_acc=0
						camera_in_this_building |= b.interior_visible_from_other_building_ext_basement(xlate, 1); // check conn building as well; expand_for_light=1
						camera_in_this_building |= camera_in_walkway;
						// generate interior detail objects during the shadow pass when the player is in the building so that it can be done in parallel with small static geom gen
						// skip drawing small object shadows for secondary camera (security camera) as an optimization
						int const inc_small(sec_camera_shadow_mode ? 0 : (camera_in_this_building ? 3 : 1));
						b.draw_room_geom(nullptr, s, amask_shader, oc, xlate, bi->ix, 1, 0, inc_small, 1); // shadow_only=1, player_in_building=1
						bool const basement_light(lpos.z < b.ground_floor_z1);
						if (!basement_light) {b.get_ext_wall_verts_no_sec(ext_parts_draw);} // add exterior walls to prevent light leaking between adjacent parts, if not basement
						else if (b.has_ext_basement()) {b.get_basement_ext_wall_verts(ext_parts_draw);} // draw basement exterior walls to block light from entering ext basement
						if (!basement_light) {b.get_walkway_end_verts(ext_parts_draw, lpos);}
						b.draw_cars_in_building(s, xlate, 1, 1); // player_in_building=1, shadow_only=1
						is_house |= b.is_house;
						bool const in_retail_room(b.check_pt_in_retail_room(lpos));
						float const player_smap_dist((in_retail_room ? RETAIL_SMAP_DSCALE : 1.0)*camera_pdu.far_);
						bool const viewer_close(dist_less_than(lpos, pre_smap_player_pos, player_smap_dist)); // Note: pre_smap_player_pos already in building space
						bool const add_player_shadow(camera_surf_collide && camera_in_this_building && viewer_close && !sec_camera_mode &&
							(actual_player_pos.z - get_bldg_player_height()) < lpos.z);
						bool const add_people_shadow((camera_in_this_building || viewer_close) && b.has_people());
						bool const enable_animations(global_building_params.enable_people_ai);

						if (add_people_shadow || add_player_shadow) {
							shader_t &shader(enable_animations ? person_shader : s);
							if (enable_animations) {select_person_shadow_shader(person_shader);}
							if (add_people_shadow) {gen_and_draw_people_in_building(ped_draw_vars_t(b, oc, shader, xlate, bi->ix, 1, 0, in_retail_room));}
							if (add_player_shadow) {draw_player_model(shader, xlate, 1);} // shadow_only=1
							if (enable_animations) {s.make_current();} // switch back to normal building shader
						}
						// since two buildings can have overlapping extended basement bcubes, we can only exit these loops if the light is in the main building itself (inc basement)
						if (b.get_part_ix_containing_pt(lpos_clamped) >= 0) {found_building = 1;}
						if (found_building) break; // only one building can contain the shadow
					} // for bi
					if (found_building) break; // only one building can contain the shadow
				} // for g
			}
			else { // draw exterior shadow maps
				for (auto g = (*i)->grid_by_tile.begin(); g != (*i)->grid_by_tile.end(); ++g) { // draw only visible tiles
					if (!building_grid_visible(xlate, g->bcube)) continue; // VFC; use exterior bcube
					(*i)->building_draw_vbo.draw_tile(s, (g - (*i)->grid_by_tile.begin()), 1);

					// draw shadow casters such as balconies that are added later as buildings come into view; won't show up until shadow map is regenerated
					for (auto bi = g->bc_ixs.begin(); bi != g->bc_ixs.end(); ++bi) {
						building_t &b((*i)->get_building(bi->ix));
						if (!b.interior || !camera_pdu.cube_visible(b.bcube + xlate)) continue; // no interior or not visible
						b.get_detail_shadow_casters(ext_parts_draw);
					}
				}
				//(*i)->building_draw_vbo.draw(s, 1); // less CPU time but more GPU work, in general seems to be slower
			}
		} // for i
		// need to draw back faces of exterior parts to handle shadows on blinds; only needed for houses, and causes problems with walkway doors
		bool const enable_back_faces(interior_shadow_maps && is_house);
		if ( enable_back_faces) {glDisable(GL_CULL_FACE);}
		ext_parts_draw.draw(s, 1, 1); // shadow_only=1, direct_draw_no_vbo=1
		ext_parts_draw.clear();
		if (!enable_back_faces) {glDisable(GL_CULL_FACE);}
		s.end_shader();
		pop_scene_xlate();
	}
	static bool check_tile_smap(bool shadow_only) {
		return (!shadow_only && world_mode == WMODE_INF_TERRAIN && shadow_map_enabled());
	}
	static bool building_grid_visible(vector3d const &xlate, cube_t const &grid_bcube, pos_dir_up const &pdu=camera_pdu) {
		return pdu.sphere_and_cube_visible_test((grid_bcube.get_cube_center() + xlate), grid_bcube.get_bsphere_radius(), (grid_bcube + xlate));
	}

	void add_interior_lights(vector3d const &xlate, cube_t &lights_bcube, bool sec_camera_mode) { // Note: non const because this caches light bcubes
		if (!draw_building_interiors || !has_interior_geom) return; // no interior
		point const camera(get_camera_pos()), camera_bs(camera - xlate);
		vector<point> points; // reused temporary
		vect_cube_with_ix_t ped_bcubes; // reused temporary
		occlusion_checker_noncity_t oc(*this, 1); // for_light=1
		bool is_first_building(1);
		//highres_timer_t timer("Add Interior Lights");

		for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
			cube_t const &grid_bcube(g->get_vis_bcube());
			if (!lights_bcube.intersects_xy  (grid_bcube)) continue; // not within light volume (too far from camera)
			if (!building_grid_visible(xlate, grid_bcube)) continue; // VFC

			for (auto bi = g->bc_ixs.begin(); bi != g->bc_ixs.end(); ++bi) {
				building_t &b(get_building(bi->ix));
				if (!b.has_room_geom()) continue; // no interior room geom, skip
				if (!lights_bcube.intersects_xy(b.bcube) && &b != player_building) continue; // not within light volume (too far from camera); allow if player building (extb)
				bool const camera_in_this_building(b.check_point_or_cylin_contained(camera_bs, 0.0, points, 1, 1, 0)); // inc_attic=1, inc_ext_basement=1, inc_roof_acc=0
				if (sec_camera_mode && !camera_in_this_building) continue; // security cameras only show lights in their building
				// limit room lights to when the player is in a building because we can restrict them to a single floor, otherwise it's too slow
				if (!camera_in_this_building && !camera_pdu.cube_visible(b.bcube + xlate) &&
					!b.interior_visible_from_other_building_ext_basement(xlate, 1) && !b.check_pt_in_or_near_walkway(camera_bs, 1, 1, 0)) continue; // VFC
				if (is_first_building) {oc.set_camera(camera_pdu, sec_camera_mode);} // setup occlusion culling on the first visible building; cur_building_only=sec_camera_mode
				is_first_building = 0;
				oc.set_exclude_bix(bi->ix);
				b.add_room_lights(xlate, bi->ix, camera_in_this_building, sec_camera_mode, oc, ped_bcubes, lights_bcube);
			} // for bi
		} // for g
	}

	void add_exterior_lights(vector3d const &xlate, cube_t &lights_bcube) const {
		for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
			if (!lights_bcube.intersects_xy  (g->bcube)) continue; // not within light volume (too far from camera)
			if (!building_grid_visible(xlate, g->bcube)) continue; // VFC; use exterior bcube

			for (auto bi = g->bc_ixs.begin(); bi != g->bc_ixs.end(); ++bi) {
				building_t const &b(get_building(bi->ix));
				if (b.ext_lights.empty()) continue;
				if (!lights_bcube.intersects_xy(b.bcube)) continue; // not within light volume (too far from camera)
				
				for (colored_sphere_t const &light : b.ext_lights) {
					if (!lights_bcube.contains_pt_xy(light.pos)) continue; // not within light volume (too far from camera)
					if (!camera_pdu.sphere_visible_test((light.pos + xlate), light.radius)) continue; // VFC
					min_eq(lights_bcube.z1(), light.pos.z - light.radius);
					max_eq(lights_bcube.z2(), light.pos.z + light.radius);
					dl_sources.push_back(light_source(light.radius, light.pos, light.pos, light.color));
				}
			} // for bi
		} // for g
	}

	struct defer_ped_draw_vars_t {
		building_t *building=nullptr;
		building_creator_t const *bc=nullptr;
		unsigned bix=0;
		void assign(building_t *b, building_creator_t const *c, unsigned ix) {assert(b); assert(c); building = b; bc = c; bix = ix;}
		bool valid() const {return (building != nullptr);}
	};

	static void ensure_city_lighting_setup(bool reflection_pass, vector3d const &xlate, bool &is_setup) {
		if (is_setup) return;
		if (!reflection_pass) {setup_city_lights(xlate);}
		is_setup = 1;
	}
	static void enable_holes_shader(shader_t &s) {
		if (!s.is_setup()) {setup_smoke_shaders(s, 0.9, 0, 0, 0, 0, 0, 0);} // min_alpha=0.9 for depth test
		else {s.enable();}
	}
	static void enable_city_shader(shader_t &s, bool use_dlights, int use_bmap, float min_alpha) {
		if (!s.is_setup()) {city_shader_setup(s, get_city_lights_bcube(), use_dlights, 1, use_bmap, min_alpha);} // use_smap=1
		else {s.enable();}
	}

	static void create_building_reflections(vector3d const &xlate) {
		bind_default_sun_moon_smap_textures();
		update_security_camera_image();
		setup_building_lights(xlate); // setup lights on first (opaque) non-shadow pass
		create_mirror_reflection_if_needed(vis_conn_bldg, xlate);
	}
	static void push_scene_xlate(vector3d const &xlate) {
		fgPushMatrix();
		translate_to(xlate);
		cur_camera_pos_xlate = xlate; // needed for correct dlights specular
	}
	static void pop_scene_xlate() {
		fgPopMatrix();
		cur_camera_pos_xlate = zero_vector;
	}

	// reflection_pass: 0 = not reflection pass, 1 = reflection for room with exterior wall,
	// 2 = reflection for room no exterior wall (can't see outside windows), 3 = reflection from mirror in a house (windows and doors need to be drawn)
	static void multi_draw(int shadow_only, int reflection_pass, vector3d const &xlate, vector<building_creator_t *> const &bcs) {
		if (bcs.empty()) return;

		if (shadow_only) {
			assert(!reflection_pass);
			multi_draw_shadow(xlate, bcs);
			return;
		}
		DebugScope scope("building_multi_draw");
		bind_default_sun_moon_smap_textures(); // bind default sun/moon smap textures
		building_t const *new_player_building(nullptr);
		//timer_t timer("Draw Buildings"); // 0.57ms (2.6ms with glFinish(), 6.3ms with building interiors)
		point const camera(get_camera_pos()), camera_bs(camera - xlate);
		int const use_bmap(global_building_params.has_normal_map);
		bool const night(is_night(WIND_LIGHT_ON_RAND)), use_city_dlights(!reflection_pass);
		bool const ref_pass_int_only(reflection_pass & REF_PASS_INT_ONLY), ref_pass_interior(reflection_pass & REF_PASS_INTERIOR); // REF_PASS_HOUSE is no longer used here
		bool const ref_pass_water(reflection_pass & REF_PASS_WATER), ref_pass_extb(reflection_pass & REF_PASS_EXTB);
		bool const not_mirror(reflection_pass & REF_PASS_NO_MIRROR), swap_front_back(reflection_pass && !not_mirror); // for mirror reflection, but not security cameras
		bool const ref_glass_floor(reflection_pass & REF_PASS_GLASS_FLOOR);
		// check for sun or moon; also need the smap pass for drawing with dynamic lights at night, so basically it's always enabled
		bool const use_tt_smap(check_tile_smap(0)); // && (night || light_valid_and_enabled(0) || light_valid_and_enabled(1)));
		bool have_windows(0), have_wind_lights(0), have_interior(0), is_city_lighting_setup(0), this_frame_camera_in_building(0), this_frame_player_in_mall(0);
		int this_frame_player_in_basement(0), this_frame_player_in_water(0), this_frame_player_in_attic(0);
		unsigned max_draw_ix(0);
		shader_t s, amask_shader, holes_shader, city_shader;

		for (auto i = bcs.begin(); i != bcs.end(); ++i) {
			assert(*i);
			have_windows     |= !(*i)->building_draw_windows.empty();
			have_wind_lights |= !(*i)->building_draw_wind_lights.empty();
			have_interior    |= (draw_building_interiors && (*i)->has_interior_geom);
			max_eq(max_draw_ix, (*i)->building_draw_vbo.get_num_draw_blocks());
			if (night) {(*i)->ensure_window_lights_vbos();}
			
			if ((*i)->is_single_tile()) { // only for tiled buildings
				(*i)->use_smap_this_frame = (use_tt_smap && try_bind_tile_smap_at_point(((*i)->grid_by_tile[0].bcube.get_cube_center() + xlate), s, 1)); // check_only=1
			}
		}
		bool const draw_interior((have_windows || global_building_params.add_city_interiors) && draw_building_interiors);
		bool const v(world_mode == WMODE_GROUND), indir(v), dlights(v), use_smap(v);
		float const min_alpha = 0.0; // 0.0 to avoid alpha test
		enable_dlight_bcubes  = 1; // using light bcubes is both faster and more correct when shadow maps are not enabled
		push_scene_xlate(xlate);
		float water_damage(0.0), crack_damage(0.0);
		building_draw_t interior_wind_draw, ext_door_draw;
		vector<building_draw_t> int_wall_draw_front, int_wall_draw_back;
		vector<vertex_range_t> per_bcs_exclude;
		building_t const *building_cont_player(nullptr);
		defer_ped_draw_vars_t defer_ped_draw_vars;
		vector<building_t *> buildings_with_cars;
		vector<point> pts;
		static brg_batch_draw_t bbd; // allocated memory is reused across building interiors
		bool const defer_people_draw_for_player_building(global_building_params.people_min_alpha > 0.0);
		vis_conn_bldg = nullptr;

		// draw building interiors with standard shader and no shadow maps; must be drawn first before windows depth pass
		if (have_interior) {
			//timer_t timer2("Draw Building Interiors");
			float const interior_draw_dist(global_building_params.interior_view_dist_scale*2.0f*(X_SCENE_SIZE + Y_SCENE_SIZE));
			float const room_geom_draw_dist(0.4*interior_draw_dist), room_geom_clear_dist(1.05*room_geom_draw_dist), room_geom_sm_draw_dist(0.14*interior_draw_dist);
			float const room_geom_int_detail_draw_dist(0.045*interior_draw_dist), room_geom_ext_detail_draw_dist(0.08*interior_draw_dist), z_prepass_dist(0.25*interior_draw_dist);
			glEnable(GL_CULL_FACE); // back face culling optimization, helps with expensive lighting shaders
			glCullFace(swap_front_back ? GL_FRONT : GL_BACK);

			// draw lit interiors; use z-prepass to reduce time taken for shading
			setup_smoke_shaders(s, 0.0, 0, 0, 0, 0, 0, 0); // everything disabled, but same shader so that vertex transforms are identical
			glPolygonOffset(1.0, 1.0);
			if (swap_front_back) {glEnable(GL_POLYGON_OFFSET_FILL);} // not sure why, but a polygon offset is required for the mirror reflection pass
			glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color rendering, we only want to write to the Z-Buffer
				
			for (auto i = bcs.begin(); i != bcs.end(); ++i) { // draw interior for the tile containing the camera
				float const ddist_scale((*i)->building_draw_windows.empty() ? 0.1 : 1.0), zpp_dist_scale(ddist_scale*z_prepass_dist);

				for (auto g = (*i)->grid_by_tile.begin(); g != (*i)->grid_by_tile.end(); ++g) {
					cube_t const &grid_bcube(g->get_vis_bcube());

					if (reflection_pass ? grid_bcube.contains_pt_xy(camera_bs) : grid_bcube.closest_dist_xy_less_than(camera_bs, zpp_dist_scale)) {
						if (!building_grid_visible(xlate, grid_bcube)) continue; // VFC
						(*i)->building_draw_interior.draw_tile(s, (g - (*i)->grid_by_tile.begin()));
					}
				}
			}
			if (swap_front_back) {glDisable(GL_POLYGON_OFFSET_FILL);}
			glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
			s.end_shader();
			set_std_depth_func_with_eq();

			if (!enabled_bldg_lights.empty()) { // used for debugging
				s.begin_color_only_shader(RED);
				begin_sphere_draw(0); // textured=0
				for (point const &lpos : enabled_bldg_lights) {draw_sphere_vbo(lpos, CAMERA_RADIUS, 16, 0);}
				end_sphere_draw();
				s.end_shader();
			}
			// Note: the best I can come up with is applying animations to both buildings and people, making sure to set animation_time to 0.0 for buildings;
			// otherwise, we would need to switch between two different shaders every time we come across a building with people in it; not very clean, but seems to work
			bool const enable_animations(global_building_params.enable_people_ai && draw_interior);
			if (enable_animations) {enable_animations_for_shader(s);}

			if (camera_in_building && player_building != nullptr) { // handle damage effects
				if (camera_bs.z < player_building->ground_floor_z1 || player_building->point_on_basement_stairs(camera_bs)) { // entering or in basement
					water_damage = player_building->water_damage;
					crack_damage = player_building->crack_damage;
					float const player_feet_zval(camera_bs.z - get_bldg_player_height());
					// incremental transition when entering/exiting the basement
					float const damage_weight(CLIP_TO_01(1.25f*(player_building->ground_floor_z1 + player_building->get_fc_thickness() -
						player_feet_zval)/player_building->get_window_vspace()));
					water_damage *= damage_weight;
					crack_damage *= damage_weight;
				}
			}
			setup_building_draw_shader(s, min_alpha, 1, 0, 0, water_damage, crack_damage); // enable_indir=1, force_tsl=0, use_texgen=0
			vector<point> points; // reused temporary
			bbd.clear_obj_models();
			int indir_bcs_ix(-1), indir_bix(-1);

			if (draw_interior) {
				per_bcs_exclude    .resize(bcs.size());
				int_wall_draw_front.resize(bcs.size());
				int_wall_draw_back .resize(bcs.size());
			}
			for (auto i = bcs.begin(); i != bcs.end(); ++i) { // draw only nearby interiors
				unsigned const bcs_ix(i - bcs.begin());
				float const door_open_dist(get_door_open_dist());
				// if there are no windows, we can wait until the player is very close to draw the interior;
				// this is generally okay when the player is flying over, and necessary for performance;
				// however, details won't show up when the player is on the ground, so use a larger scale in that case
				bool const int_not_visible((*i)->building_draw_windows.empty());
				float const ddist_scale(int_not_visible ? (camera_surf_collide ? 0.1 : 0.05) : 1.0), ddist_scale_sq(ddist_scale*ddist_scale);
				float const int_draw_dist_sq(ddist_scale_sq*interior_draw_dist*interior_draw_dist);
				float const rgeom_clear_dist_sq(ddist_scale_sq*room_geom_clear_dist*room_geom_clear_dist);
				float const rgeom_draw_dist_sq(ddist_scale_sq*room_geom_draw_dist*room_geom_draw_dist);
				float const rgeom_sm_draw_dist_sq(ddist_scale_sq*room_geom_sm_draw_dist*room_geom_sm_draw_dist);
				float const rgeom_int_detail_dist_sq(ddist_scale_sq*room_geom_int_detail_draw_dist*room_geom_int_detail_draw_dist);
				float const rgeom_ext_detail_dist_sq(ddist_scale_sq*room_geom_ext_detail_draw_dist*room_geom_ext_detail_draw_dist);
				occlusion_checker_noncity_t oc(**i);
				bool is_first_tile(1), can_break_from_loop(0);

				for (auto g = (*i)->grid_by_tile.begin(); g != (*i)->grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
					cube_t const &grid_bcube(g->get_vis_bcube());
					// for the reflection pass, we only need to look at the grid containing the building with the mirror, which must be the player's building
					if (reflection_pass && !grid_bcube.contains_pt_xy(camera_bs)) continue; // not the correct tile
					float const gdist_sq(p2p_dist_sq(camera_bs, grid_bcube.closest_pt(camera_bs)));

					if (!reflection_pass && gdist_sq > rgeom_clear_dist_sq && g->has_room_geom) { // need to clear room geom
						for (auto bi = g->bc_ixs.begin(); bi != g->bc_ixs.end(); ++bi) {(*i)->get_building(bi->ix).clear_room_geom();}
						g->has_room_geom = 0;
					}
					if (gdist_sq > int_draw_dist_sq)               continue; // too far
					if (!building_grid_visible(xlate, grid_bcube)) continue; // VFC
					if (is_first_tile) {(*i)->ensure_interior_geom_vbos();} // we need the interior geom at this point, even if it's the reflection pass
					if (crack_damage > 0.0) {s.add_uniform_float("crack_weight", crack_damage);} // crack damage for interior
					// Note: in cases where another building's extended basement is visible, we can't set crack_weight and wet_effect per-building, only per-tile
					(*i)->building_draw_interior.draw_tile(s, (g - (*i)->grid_by_tile.begin()), 0, crack_damage); // shadow_only=0
					// iterate over nearby buildings in this tile and draw interior room geom, generating it if needed
					if (gdist_sq > rgeom_draw_dist_sq) continue; // too far
					if (crack_damage > 0.0) {s.add_uniform_float("crack_weight", 0.0);} // no crack damage for room objects
					if (is_first_tile && !reflection_pass) {oc.set_camera(camera_pdu);} // setup occlusion culling on the first visible tile
					if (!ref_pass_interior) {bbd.next_tile(g->bcube);} // only needed for exterior geom; always uses main/exterior bcube
					is_first_tile = 0;

					for (auto bi = g->bc_ixs.begin(); bi != g->bc_ixs.end(); ++bi) {
						building_t &b((*i)->get_building(bi->ix));
						if (!b.interior) continue; // no interior, skip
						float const bdist_sq(p2p_dist_sq(camera_bs, b.bcube.closest_pt(camera_bs)));
						//if (bdist_sq > rgeom_clear_dist_sq) {b.clear_room_geom(); continue;} // too far away - is this useful?
						if (bdist_sq > rgeom_draw_dist_sq) continue; // too far away
						bool player_in_building_bcube(b.bcube.contains_pt_xy(camera_bs) || b.point_in_extended_basement(camera_bs)); // player within building's bcube
						bool const ext_basement_conn_visible(b.interior_visible_from_other_building_ext_basement(xlate));
						if (reflection_pass && !player_in_building_bcube && !ext_basement_conn_visible) continue; // not the correct building
						bool const debug_draw(0 && b.interior->has_backrooms); // TESTING
						
						if (b.check_pt_in_or_near_walkway(camera_bs, 1, 1, 1)) { // owned_only=1, inc_open_door=1, inc_conn_room=1
							if (toggle_room_light) {b.toggle_walkway_lights(camera_bs);}
							player_in_building_bcube = 1; // walkways count as in building bcube
						}
						if (!debug_draw && !player_in_building_bcube && !ext_basement_conn_visible && !camera_pdu.cube_visible(b.bcube + xlate)) continue; // VFC
						b.maybe_gen_chimney_smoke();
						bool const camera_near_building(player_in_building_bcube || (!b.doors.empty() && b.bcube.contains_pt_xy_exp(camera_bs, door_open_dist)));
						bool cant_see_inside(0);

						if (!debug_draw && !ext_basement_conn_visible) {
							// check if player is outside a windowless building (city office building); need to account for open doors
							if (!player_in_building_bcube && !b.has_windows()) {
								if (!b.point_near_ext_door(camera_bs, 20.0*door_open_dist)) continue; // too far away (use larger dist for door steps and ext door signs)
								if (!camera_near_building) {cant_see_inside = 1;} // can see exterior objects, but not interiors
							}
							else if ((display_mode & 0x08) && !player_in_building_bcube && b.is_entire_building_occluded(camera_bs, oc)) continue; // check occlusion
						}
						// draw interior detail objects if player is in the building (inc ext basement), even if far from the building center
						unsigned inc_small(bdist_sq < rgeom_sm_draw_dist_sq);
						if      (cant_see_inside)                                  {inc_small = 4;} // only exterior detail objects
						else if (player_in_building_bcube)                         {inc_small = 3;} // include interior and exterior detail objects
						else if (ext_basement_conn_visible)                        {inc_small = 3;} // include interior and exterior detail objects
						else if (inc_small && bdist_sq < rgeom_int_detail_dist_sq) {inc_small = 3;} // include interior and exterior detail objects
						else if (inc_small && bdist_sq < rgeom_ext_detail_dist_sq) {inc_small = 2;} // include exterior detail objects
						if (debug_draw) {inc_small = 3;} // TESTING
						bool const player_in_bldg(debug_draw || player_in_building_bcube);
						if (ext_basement_conn_visible) {s.add_uniform_float("wet_effect", 0.0);} // disable for non-player building
						b.gen_and_draw_room_geom(&bbd, s, amask_shader, oc, xlate, bi->ix, 0, reflection_pass, inc_small, player_in_bldg, ext_basement_conn_visible); // shadow_only=0
						if (ext_basement_conn_visible) {s.add_uniform_float("wet_effect", water_damage);}
						g->has_room_geom = 1;
						if (!draw_interior) continue;
						
						// when player is in the building (not attic or ext basement), draw people later so that alpha blending of hair against ext walls and windows works properly
						if (defer_people_draw_for_player_building && player_in_building_bcube && b.has_people() && b.check_point_or_cylin_contained(camera_bs, 0.0, points, 0, 0, 0)) {
							defer_ped_draw_vars.assign(&b, *i, bi->ix);
						}
						else {gen_and_draw_people_in_building(ped_draw_vars_t(b, oc, s, xlate, bi->ix, 0, reflection_pass));} // draw people in this building
						// there currently shouldn't be any parked cars visible in mirrors or security cameras, so skip them in the reflection pass
						if (!reflection_pass && b.has_cars_to_draw(player_in_building_bcube)) {buildings_with_cars.push_back(&b);}					

						if ((*i)->get_is_city()) { // check for nearby pedestrians in city buildings and open doors for them
							float const ped_od(0.4*door_open_dist); // smaller than player dist
							pts.clear();
							cube_t door_test_cube(b.bcube);
							door_test_cube.expand_by_xy(ped_od);
							get_pedestrians_in_area(door_test_cube, bi->ix, pts); // is this thread safe?
							b.get_all_nearby_ext_door_verts(ext_door_draw, s, pts, ped_od);
						}
						// check the bcube rather than check_point_or_cylin_contained() so that it works with roof doors that are outside any part?
						if (!camera_near_building && !ext_basement_conn_visible) { // camera not near building or ext basement conn
							if (!reflection_pass) {b.player_not_near_building();}
							continue;
						}
						if (ref_pass_interior) continue; // interior room, don't need to draw windows and exterior doors
						// and draw opened door; update_state if not ref pass
						bool const had_open_door(b.get_nearby_ext_door_verts(ext_door_draw, s, camera_bs, cview_dir, door_open_dist, !reflection_pass, 0)); // only_open=0
						bool const camera_in_this_building(b.check_point_or_cylin_contained(camera_bs, 0.0, points, 1, 1, 1)); // inc_attic=1, inc_ext_basement=1, inc_roof_acc=1
						bool const player_in_bldg_bc_or_door(player_in_building_bcube || had_open_door);
						
						if (!reflection_pass && (camera_in_this_building || !this_frame_camera_in_building)) { // player in this building, or near but not inside another
							// disable grass in building part(s) containing the player
							b.update_grass_exclude_at_pos(camera_bs, xlate, camera_in_this_building);
						}
						if (!reflection_pass && player_in_bldg_bc_or_door) {b.update_animals(camera_bs, bi->ix);}
						
						// Note: if we skip this check and treat all walls/windows as front/containing part, this almost works, but will skip front faces of other buildings
						if (!camera_in_this_building) { // camera not in building
							if (ext_basement_conn_visible && animate2) {b.update_player_interact_objects(camera_bs);} // need to at least update door open/close state
							if (ext_basement_conn_visible && !reflection_pass) {vis_conn_bldg = &b;} // for now we only support one visible connected building

							if (had_open_door && !reflection_pass && b.glass_floor_visible(xlate, 1)) { // from_outside_building=1
								b.draw_glass_surfaces(xlate);
								s.make_current();
							}
							continue;
						}
						// we should get here for at most one building
						// pass in camera pos to only include the part that contains the camera to avoid drawing artifacts when looking into another part of the building
						// neg offset to move windows on the inside of the building's exterior wall;
						// since there are no basement windows, we should treat the player as being in the part above so that windows are drawn correctly through the basement stairs
						assert(bcs_ix < int_wall_draw_front.size() && bcs_ix < int_wall_draw_back.size());
						point pt_ag(camera_bs);
						max_eq(pt_ag.z, (b.ground_floor_z1 + b.get_floor_thickness()));
						b.get_all_drawn_window_verts(interior_wind_draw, 0, -0.1, &pt_ag, 0, 1); // lights_pass=0, no_skylights=0, draw_int_windows=1
						b.get_split_int_window_wall_verts(int_wall_draw_front[bcs_ix], int_wall_draw_back[bcs_ix], pt_ag, 0);
						building_cont_player    = &b; // there can be only one
						if (!interior_wind_draw.empty() && !ref_pass_interior) {per_bcs_exclude[bcs_ix] = b.ext_side_qv_range;} // only if there are drawn windows
						if (reflection_pass) continue; // don't execute the code below
						if (display_mode & 0x20) {b.debug_people_in_building(s, camera_bs);} // debug visualization
						float const basement_z_adj(2.0*BASEMENT_ENTRANCE_SCALE*b.get_floor_thickness()); // adjust to prevent problems when camera is close to the plane
						this_frame_camera_in_building = 1;
						this_frame_player_in_basement =   b.check_player_in_basement(camera_bs - basement_z_adj*plus_z); // set once
						this_frame_player_in_mall     =   b.point_in_mall(camera_bs);
						this_frame_player_in_attic    =  (b.point_in_attic(camera_bs) ? (b.has_attic_window ? 1 : 2) : 0);
						this_frame_player_in_water    =   b.point_in_water_area(camera_bs, 1); // full_room_height=1
						if (this_frame_player_in_water && b.point_in_water_area(camera_bs, 0)) {this_frame_player_in_water = 2;} // full_room_height=0; test for underwater
						
						if (!camera_surf_collide) { // handle player clipping/flying into or out of an elevator
							if (!b.point_in_elevator(camera_bs, 1)) {player_in_elevator = 0;} // check_elevator_car=1
							else {max_eq(player_in_elevator, 1);} // at least in an elevator
						}
						// player can only be in one basement or attic, except for extended basement connector rooms;
						// be conservative and don't break if the player is in the basement and this building has any connections to other basements
						can_break_from_loop |= ((this_frame_player_in_basement >= 2 && !b.has_conn_info()) || this_frame_player_in_attic == 2);
						new_player_building = &b;
						b.register_player_in_building(camera_bs, bi->ix); // required for AI following logic
						if (enable_building_indir_lighting()) {indir_bcs_ix = bcs_ix; indir_bix = bi->ix;} // compute indirect lighting for this building
						// run any player interaction logic here
						b.update_security_cameras(camera_bs);
						if (toggle_room_light  ) {b.toggle_room_light(camera_bs);}
						if (building_action_key) {b.apply_player_action_key(camera_bs, cview_dir, (building_action_key-1), 0);} // check_only=0
						else {can_do_building_action = b.apply_player_action_key(camera_bs, cview_dir, 0, 1);} // mode=0, check_only=1
						b.player_pickup_object(camera_bs, cview_dir);
						if (teleport_to_screenshot) {b.maybe_teleport_to_screenshot();}
						if (animate2) {b.update_player_interact_objects(camera_bs);} // update dynamic objects if the player is in the building
						building_occluder = b.get_best_occluder(camera_bs);
					} // for bi
					if (can_break_from_loop) break; // done
				} // for g
				if (can_break_from_loop) break; // done
			} // for i
			bbd.draw_and_clear(s);
			set_std_depth_func(); // restore
			glDisable(GL_CULL_FACE);

			// draw lower part of player model if not flying; doesn't work well when crouching
			if (!reflection_pass && this_frame_camera_in_building && camera_surf_collide && global_building_params.show_player_model) {
				glDisable(GL_DEPTH_CLAMP);
				draw_player_model(s, xlate, 0);
				setup_depth_clamp(); // restore
			}
			if (!reflection_pass) { // update once; non-interior buildings (such as city buildings) won't update this
				camera_in_building = this_frame_camera_in_building;
				player_in_basement = this_frame_player_in_basement;
				player_in_mall     = this_frame_player_in_mall;
				player_in_attic    = this_frame_player_in_attic;
				player_in_water    = this_frame_player_in_water;
				building_has_open_ext_door = !ext_door_draw.empty();
			}
			reset_interior_lighting_and_end_shader(s);
			reflection_shader.clear();

			// update indir lighting using ray casting
			if (indir_bcs_ix >= 0 && indir_bix >= 0) {indir_tex_mgr.create_for_building(bcs[indir_bcs_ix]->get_building(indir_bix), indir_bix, camera_bs);}
			else if (!reflection_pass) {end_building_rt_job();}
			
			if (draw_interior && !interior_wind_draw.empty() && !ref_pass_interior) {
				// draw interior windows to cut out holes; write to stencil buffer, use stencil test for back facing building walls
				enable_holes_shader(holes_shader);
				setup_stencil_buffer_write();
				glStencilOpSeparate((swap_front_back ? GL_BACK : GL_FRONT), GL_KEEP, GL_KEEP, GL_KEEP); // ignore front faces
				glStencilOpSeparate((swap_front_back ? GL_FRONT : GL_BACK), GL_KEEP, GL_KEEP, GL_INCR); // mark stencil on back faces
				glDepthMask(GL_FALSE);
				interior_wind_draw.draw(holes_shader, 0, 1); // draw back facing windows; direct_draw_no_vbo=1
				glDepthMask(GL_TRUE);
				end_stencil_write();
				holes_shader.disable();
			}
		} // end have_interior
		if (!reflection_pass) { // update player_building state
			if (new_player_building == nullptr) {register_player_not_in_building();}

			if (new_player_building != player_building) { // building transition
				if (new_player_building) {new_player_building->register_player_enter_building();}
				if (player_building    ) {player_building    ->register_player_exit_building (new_player_building != nullptr);}
				player_building = new_player_building;
			}
			toggle_room_light = teleport_to_screenshot = 0; building_action_key = 0; // reset these even if the player wasn't in a building
		}
		if (draw_interior) {
			if (!ref_pass_extb) { // skip for extended basement room reflections
				// draw back faces of buildings, which will be interior walls
				cube_t player_part;
				
				// since walls are mostly XY axis aligned, we can use both axes for the texture 's' component scale and Z for the 't' component scale;
				// this doesn't really work for non-cube buildings though, in particular on near 45 degree edges where the delta_x cancels with the delta_y;
				// so set special texgen mode so that X and Y are of opposite signs and don't cancel at near 45 degree edges; this fixes cylinders but not 5-6 sides
				if (player_building != nullptr && player_building->num_sides > 8) { // player in non-cube cylinder-like building
					player_part = player_building->get_part_containing_pt(camera_bs);
				}
				bool const diag_texgen_mode(!player_part.is_all_zeros());
				// water damage is needed here as well to apply to exterior basement walls, but cracks don't apply to concrete blocks
				setup_building_draw_shader(s, min_alpha, 1, 1, (diag_texgen_mode ? 2 : 1), water_damage, 0.0); // enable_indir=1, force_tsl=1, use_texgen=1|2
				glEnable(GL_CULL_FACE);
				glCullFace(swap_front_back ? GL_BACK : GL_FRONT); // draw back faces

				for (auto i = bcs.begin(); i != bcs.end(); ++i) {
					if ((*i)->empty() || !(*i)->has_interior_to_draw()) continue; // no buildings or no interiors
					unsigned const bcs_ix(i - bcs.begin());
					vertex_range_t const *exclude(nullptr);
					building_mat_t const &mat((*i)->buildings.front().get_material()); // Note: assumes all wall textures have a consistent tscale
					vector3d texgen_origin;

					if (diag_texgen_mode) { // use building player part center
						texgen_origin.assign(player_part.xc(), player_part.yc(), 0.0);
					}
					else { // translate texture near the camera to get better tex coord resolution; make a multiple of tscale to avoid visible shift
						texgen_origin.assign(xoff2*DX_VAL, yoff2*DY_VAL, 0.0);
						for (unsigned d = 0; d < 2; ++d) {texgen_origin[d] = mat.wall_tex.tscale_x*int(texgen_origin[d]/mat.wall_tex.tscale_x);}
					}
					float const side_tscale(2.0f*mat.wall_tex.tscale_x);
					s.add_uniform_vector3d("texgen_origin", texgen_origin);
					setup_texgen_full(side_tscale, side_tscale, 0.0, 0.0, 0.0, 0.0, 2.0f*mat.wall_tex.tscale_y, 0.0, s, 0);
				
					if (!per_bcs_exclude.empty()) { // draw this range using stencil test but the rest of the buildings without stencil test
						vertex_range_t const &vr(per_bcs_exclude[bcs_ix]);

						if (vr.draw_ix >= 0) { // nonempty
							exclude = &vr; // use this exclude
							glEnable(GL_STENCIL_TEST);
							glStencilFunc(GL_EQUAL, 0, ~0U); // keep if stencil bit has not been set by above pass
							glStencilOpSeparate(GL_FRONT_AND_BACK, GL_KEEP, GL_KEEP, GL_KEEP);
							int_wall_draw_front[bcs_ix].draw(s, 0, 1); // draw back facing walls for front part of building with    stencil test
							glDisable(GL_STENCIL_TEST);
							int_wall_draw_back [bcs_ix].draw(s, 0, 1); // draw back facing walls for back  part of building without stencil test
						}
					}
					(*i)->building_draw_int_ext_walls.draw(s, 0, 0, 0, exclude); // exterior walls only, no stencil test
					s.add_uniform_vector3d("texgen_origin", zero_vector);
				} // for i
				reset_interior_lighting_and_end_shader(s);

				// draw parked cars in building parking garages or house garages
				if (!buildings_with_cars.empty()) {
					glDisable(GL_CULL_FACE); // no back face culling for cars
					for (auto const &b : buildings_with_cars) {b->draw_cars_in_building(s, xlate, camera_in_building, 0);} // shadow_only=0
					if (s.is_setup()) {reset_interior_lighting_and_end_shader(s);}
					glEnable(GL_CULL_FACE);
				}
				if (DRAW_EXT_REFLECTIONS || !reflection_pass) {
					// if we're not by an exterior door, draw the back sides of exterior doors as closed; always draw roof geometry
					// this is required for drawing objects such as the underside of roof overhangs and roof doors and their covers
					int const tex_filt_mode(ext_door_draw.empty() ? 2 : 3);
					bool const enable_indir(camera_in_building); // need to enable indir lighting when drawing the back sides of exterior doors
					setup_building_draw_shader(s, min_alpha, enable_indir, 1, 0); // force_tsl=1, use_texgen=0, damage=0.0
					for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw_vbo.draw(s, 0, 0, tex_filt_mode);}
					reset_interior_lighting_and_end_shader(s);
				}
			} // end !ref_pass_extb
			glCullFace(swap_front_back ? GL_FRONT : GL_BACK); // draw front faces

			// draw people in the player's building here with alpha mask enabled
			if (defer_ped_draw_vars.valid() || (reflection_pass && !ref_pass_water)) {
				if (global_building_params.enable_people_ai) {enable_animations_for_shader(s);}
				setup_building_draw_shader(s, global_building_params.people_min_alpha, 1, 0, 0); // enable_indir=1, force_tsl=0, use_texgen=0, damage=0.0

				if (defer_ped_draw_vars.valid()) {
					occlusion_checker_noncity_t oc(*defer_ped_draw_vars.bc);
					if (!reflection_pass) {oc.set_camera(camera_pdu);} // setup occlusion culling
					gen_and_draw_people_in_building(ped_draw_vars_t(*defer_ped_draw_vars.building, oc, s, xlate, defer_ped_draw_vars.bix, 0, reflection_pass));
				}
				// draw player reflection last so that alpha blending of hair works properly; not visible in water reflections or glass floor (due to low Fresnel term)
				if (reflection_pass && !ref_pass_water && !not_mirror && !ref_glass_floor) {draw_player_model(s, xlate, 0);} // shadow_only=0
				reset_interior_lighting_and_end_shader(s);
			}
			if (!ref_pass_interior && have_buildings_ext_paint()) { // draw spraypaint/markers on building exterior walls/windows, if needed
				glDisable(GL_CULL_FACE);
				setup_building_draw_shader(s, DEF_CITY_MIN_ALPHA, 1, 1, 0); // alpha test, enable_indir=1, force_tsl=1, use_texgen=0, damage=0.0
				draw_buildings_ext_paint(s);
				reset_interior_lighting_and_end_shader(s);
				glEnable(GL_CULL_FACE);
			}
			if (!reflection_pass && player_in_ext_basement()) {
				player_building->draw_water(xlate);
				if (vis_conn_bldg) {vis_conn_bldg->draw_water(xlate);} // check any visible building as well
			}
			if (!ref_pass_interior && bbd.has_ext_geom()) { // skip for interior room reflections
				glDisable(GL_CULL_FACE);
				ensure_city_lighting_setup(reflection_pass, xlate, is_city_lighting_setup); // needed for dlights to work
				glEnable(GL_CULL_FACE); // above call may create shadow maps and disable face culling, so make sure it's re-enabled
				enable_city_shader(city_shader, use_city_dlights, use_bmap, min_alpha);
				bbd.draw_and_clear_ext_tiles(city_shader, xlate); // draw after ext walls but before windows so that alpha blending works properly
				city_shader.disable();
			}
			bbd.clear_ext_tiles(); // required, even if there's no ext_geom(), because tile bboxes may still be nonempty

			if (!reflection_pass) { // draw windows and doors in depth pass to create holes
				enable_holes_shader(holes_shader); // need same shader to avoid z-fighting
				glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color writing, we only want to write to the Z-Buffer
				if (!ext_door_draw.empty()) {glDisable(GL_DEPTH_CLAMP);}
				for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw_windows.draw(holes_shader, 0);} // draw windows on top of other buildings
				glEnable(GL_DEPTH_CLAMP); // make sure holes are not clipped by the near plane
				ext_door_draw.draw(holes_shader, 0, 1); // direct_draw_no_vbo=1
				setup_depth_clamp(); // restore
				glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
				holes_shader.end_shader();
			}
			glDisable(GL_CULL_FACE);
		} // end draw_interior
		draw_candle_flames();

		// everything after this point is part of the building exteriors and uses city lights rather than building room lights;
		// when the player is in the extended basement we still need to draw the exterior wall and door
		if ((reflection_pass && (!DRAW_EXT_REFLECTIONS || ref_pass_int_only)) || player_cant_see_outside_building()) {
			// early exit for player fully in basement or attic, or house reflections, if enabled
			pop_scene_xlate();
			enable_dlight_bcubes = 0;
			return;
		}
		// main/batched draw pass
		ensure_city_lighting_setup(reflection_pass, xlate, is_city_lighting_setup);
		// Note: indir and dlights are set for ground mode only
		bool const keep_alpha = 1; // required for fog on windows
		setup_smoke_shaders(s, min_alpha, 0, keep_alpha, indir, 1, dlights, 0, 0, (use_smap ? 2 : 1), use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1

		if (!reflection_pass) { // don't want to do this in the reflection pass
			for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw.init_draw_frame();}
		}
		glEnable(GL_CULL_FACE);
		glCullFace(swap_front_back ? GL_FRONT : GL_BACK);
		if (!ext_door_draw.empty()) {glDisable(GL_DEPTH_CLAMP);} // if an exterior door was drawn, make sure we don't clamp the walls over the holes

		// draw front faces of buildings, even in the reflection pass; nearby shadowed buildings will be drawn later
		for (unsigned ix = 0; ix < max_draw_ix; ++ix) {
			for (auto i = bcs.begin(); i != bcs.end(); ++i) {
				if (!(*i)->use_smap_this_frame) {(*i)->building_draw_vbo.draw_block(s, ix, 0);}
			}
		}
		set_std_depth_func_with_eq();
		glPolygonOffset(-1.0, -1.0); // useful for avoiding z-fighting on building windows

		if (have_windows) { // draw windows, front facing only (not viewed from interior)
			enable_blend();
			glDepthMask(GL_FALSE); // disable depth writing
			glEnable(GL_POLYGON_OFFSET_FILL);

			for (auto i = bcs.begin(); i != bcs.end(); ++i) { // draw windows on top of other buildings
				// need to swap opaque window texture with transparent texture for this draw pass
				bool const transparent_windows(draw_interior && (*i)->has_interior_to_draw() && !reflection_pass);
				if (transparent_windows) {(*i)->building_draw_windows.toggle_transparent_windows_mode();}
				(*i)->building_draw_windows.draw(s, 0);
				if (transparent_windows) {(*i)->building_draw_windows.toggle_transparent_windows_mode();}
			}
			//interior_wind_draw.draw(s, 0, 1); // draw opaque front facing windows of building the player is in; direct_draw_no_vbo=1
			glDisable(GL_POLYGON_OFFSET_FILL);
			glDepthMask(GL_TRUE); // re-enable depth writing
			disable_blend();
		}
		glDisable(GL_CULL_FACE);
		bool const use_smap_pass(use_tt_smap && !reflection_pass);
		if (!use_smap_pass) {bbd.draw_obj_models(s, xlate, 0);} // draw models here if the code below won't be run; shadow_only=0
		if (building_cont_player) {building_cont_player->write_basement_entrance_depth_pass(s);} // drawn last
		s.end_shader();

		// post-pass to render building exteriors in nearby tiles that have shadow maps; shadow maps don't work right when using reflections
		if (use_smap_pass) {
			//timer_t timer2("Draw Buildings Smap"); // 0.3
			enable_city_shader(city_shader, use_city_dlights, use_bmap, min_alpha);
			float const draw_dist(get_tile_smap_dist() + 0.5f*(X_SCENE_SIZE + Y_SCENE_SIZE));
			int const norm_bias_scale_loc(city_shader.get_uniform_loc("norm_bias_scale"));
			assert(norm_bias_scale_loc >= 0);
			// cull back faces to avoid lighting/shadows on inside walls of building interiors;
			// disable when there are no interiors so that the bottom surfaces of roof overhangs/gutters are drawn
			if (draw_interior) {glEnable(GL_CULL_FACE);}

			for (auto i = bcs.begin(); i != bcs.end(); ++i) {
				bool const single_tile((*i)->is_single_tile()), no_depth_write(!single_tile), transparent_windows(draw_interior && (*i)->has_interior_to_draw());
				if (single_tile && !(*i)->use_smap_this_frame) continue; // optimization
				if (no_depth_write) {glDepthMask(GL_FALSE);} // disable depth writing

				for (auto g = (*i)->grid_by_tile.begin(); g != (*i)->grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
					if (single_tile && (*i)->use_smap_this_frame) {} // not drawn in main/nonshadow pass, so must be drawn here
					else if (!g->bcube.closest_dist_less_than(camera_bs, draw_dist)) continue; // too far; uses exterior bcube
					if (!building_grid_visible(xlate, g->bcube)) continue; // VFC; use exterior bcube
					unsigned lod_level(0);
					if (!try_bind_tile_smap_at_point((g->bcube.get_cube_center() + xlate), city_shader, 0, &lod_level)) continue; // no shadow maps - not drawn in this pass
					// increase bias with smap texture LOD to make it constant per texel to avoid artifacts on distant tiles using low resolution smaps
					city_shader.set_uniform_float(norm_bias_scale_loc, DEF_NORM_BIAS_SCALE*(1.0 + max(0, ((int)lod_level-1))));
					unsigned const tile_id(g - (*i)->grid_by_tile.begin());
					// Note: we could skip detail materials like trim for tiles that are further, but it's unclear if that would make much difference
					(*i)->building_draw_vbo.draw_tile(city_shader, tile_id);

					if (!(*i)->building_draw_windows.empty()) {
						enable_blend();
						glEnable(GL_POLYGON_OFFSET_FILL);
						if (!no_depth_write) {glDepthMask(GL_FALSE);} // always disable depth writing
						if (transparent_windows) {(*i)->building_draw_windows.toggle_transparent_windows_mode();}
						(*i)->building_draw_windows.draw_tile(city_shader, tile_id); // draw windows on top of other buildings
						if (transparent_windows) {(*i)->building_draw_windows.toggle_transparent_windows_mode();}
						if (!no_depth_write) {glDepthMask(GL_TRUE);} // always re-enable depth writing
						glDisable(GL_POLYGON_OFFSET_FILL);
						disable_blend();
					}
				} // for g
				if (no_depth_write) {glDepthMask(GL_TRUE);} // re-enable depth writing
			} // for i
			city_shader.set_uniform_float(norm_bias_scale_loc, DEF_NORM_BIAS_SCALE); // restore the default
			if (draw_interior) {glDisable(GL_CULL_FACE);}
			bbd.draw_obj_models(city_shader, xlate, 0); // shadow_only=0
			city_shader.end_shader();
		}
		if (night && have_wind_lights) { // add night time random lights in windows
			enable_blend();
			glDepthMask(GL_FALSE); // disable depth writing
			float const low_v(0.5 - WIND_LIGHT_ON_RAND), high_v(0.5 + WIND_LIGHT_ON_RAND), lit_thresh_mult(1.0 + 2.0*CLIP_TO_01((light_factor - low_v)/(high_v - low_v)));
			s.set_vert_shader("window_lights");
			s.set_frag_shader("linear_fog.part+window_lights");
			s.set_prefix("#define FOG_FADE_TO_TRANSPARENT", 1);
			setup_tt_fog_pre(s);
			s.begin_shader();
			s.add_uniform_float("lit_thresh_mult", lit_thresh_mult); // gradual transition of lit window probability around sunset
			setup_tt_fog_post(s);
			for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw_wind_lights.draw(s, 0);} // add bloom?
			glDepthMask(GL_TRUE); // re-enable depth writing
			disable_blend();
		}
		if (!ext_door_draw.empty()) {setup_depth_clamp();} // restore
		glCullFace(GL_BACK);
		set_std_depth_func();
		pop_scene_xlate();
		enable_dlight_bcubes = 0;
	}

	static void draw_player_building_transparent(int reflection_pass, vector3d const &xlate) {
		// draw glass materials such as floors for the player's building
		if (reflection_pass || !draw_building_interiors) return;
		if (player_building == nullptr || !player_building->glass_floor_visible(xlate)) return;
		push_scene_xlate(xlate);
		player_building->draw_glass_surfaces(xlate);
		pop_scene_xlate();
	}

	void draw_building_lights(vector3d const &xlate) { // add night time lights to buildings; non-const because it modifies building_lights
		if (empty() || !is_night(WIND_LIGHT_ON_RAND)) return;
		//timer_t timer("Building Lights"); // 0.06ms
		set_additive_blend_mode();
		enable_blend();
		glDepthMask(GL_FALSE); // disable depth writing
		vector3d const max_extent(get_buildings_max_extent());
		float const draw_dist(20.0*max_extent.mag());
		point const camera(get_camera_pos() - xlate); // in building space
		colorRGBA const light_colors[16] = {RED,RED,RED,RED,RED,RED,RED,RED, BLUE,BLUE,BLUE,BLUE, WHITE,WHITE, YELLOW, GREEN};

		for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) {
			if (!g->bcube.closest_dist_less_than(camera, draw_dist)) continue; // too far away; use exterior bcube
			if (!camera_pdu.cube_visible(g->bcube + xlate)) continue;

			for (auto i = g->bc_ixs.begin(); i != g->bc_ixs.end(); ++i) {
				building_t const &b(get_building(i->ix));
				if (b.details.empty()) continue;
				if (!is_night((((321*i->ix) & 7)/7.0)*WIND_LIGHT_ON_RAND)) continue; // gradually turn on
				if (!b.bcube.closest_dist_less_than(camera, draw_dist)) continue; // too far away
				if (!camera_pdu.cube_visible(b.bcube + xlate)) continue;
				
				for (auto j = b.details.begin(); j != b.details.end(); ++j) {
					if (j->type != ROOF_OBJ_ANT) continue; // not an antenna
					unsigned const num_segs(max(1U, (((123*i->ix) & 3) + unsigned(6.0*j->dz()/max_extent.z)))); // some mix of height and randomness
					point const center(j->get_cube_center());
					point pos(point(center.x, center.y, j->z2()) + xlate);
					float const radius(1.2f*(j->dx() + j->dy())), z_step(0.6*j->dz()/num_segs);
					float const alpha(min(1.0f, 1.5f*(1.0f - p2p_dist(camera, center)/draw_dist))); // fade with distance
					colorRGBA const color(light_colors[i->ix & 15], alpha);

					for (unsigned n = 0; n < num_segs; ++n) { // distribute lights along top half of antenna
						building_lights.add_pt(sized_vert_t<vert_color>(vert_color(pos, color), radius));
						pos.z -= z_step;
					}
				} // for j
			} // for i
		} // for g
		building_lights.draw_and_clear(BLUR_TEX, 0.0, 0, 1, 0.005); // use geometry shader for unlimited point size
		glDepthMask(GL_TRUE); // re-enable depth writing
		disable_blend();
		set_std_blend_mode();
	}

	void get_all_window_verts(building_draw_t &bdraw, bool light_pass) { // for exterior drawing
		bdraw.clear();

		for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
			bdraw.cur_tile_id = (g - grid_by_tile.begin());
			for (auto i = g->bc_ixs.begin(); i != g->bc_ixs.end(); ++i) {get_building(i->ix).get_all_drawn_window_verts(bdraw, light_pass);}
		}
		bdraw.finalize(grid_by_tile.size());
	}
	void get_interior_drawn_verts() {
		// pre-allocate interior wall, celing, and floor verts, assuming all buildings have the same materials
		tid_vert_counter_t vert_counter;

		for (auto b = buildings.begin(); b != buildings.end(); ++b) {
			if (!b->interior) continue; // no interior
			unsigned const num_elevators(b->interior->elevators.size()), ceil_nverts(b->skip_top_of_ceilings() ? 4 : 8);
			// Note: here we use 14 verts per wall rather than the expected 16 due to estimated hidden surface culling
			unsigned const nv_wall(14*(b->interior->walls[0].size() + b->interior->walls[1].size()) + 16*b->interior->landings.size() + 16*b->has_attic() + 36*num_elevators);
			vert_counter.update_count(b->get_material().wall_tex.tid, nv_wall);
			vert_counter.update_count(FENCE_TEX, 12*num_elevators);

			for (cube_t const &f : b->interior->floors) {
				tid_nm_pair_t tex;
				b->get_floor_tex_and_color(f, tex);
				vert_counter.update_count(tex.tid, 4);
			}
			for (cube_t const &c : b->interior->ceilings) {
				tid_nm_pair_t tex;
				b->get_ceil_tex_and_color(c, tex);
				vert_counter.update_count(tex.tid, ceil_nverts);
			}
			if (b->has_attic()) {
				tid_nm_pair_t const attic_tex(b->get_attic_texture());

				for (tquad_with_ix_t const &i : b->roof_tquads) {
					if (b->is_attic_roof(i, 0)) {vert_counter.update_count(attic_tex.tid, i.npts);}
				}
			}
		}
		for (unsigned i = 0; i < tid_mapper.get_num_slots(); ++i) {
			unsigned const count(vert_counter.get_count(i));
			if (count > 0) {building_draw_interior.reserve_verts(tid_nm_pair_t(i), count);}
		}
		// generate vertex data
		building_draw_interior.clear();
		building_draw_int_ext_walls.clear();

		for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
			building_draw_interior.cur_tile_id = (g - grid_by_tile.begin());

			for (auto i = g->bc_ixs.begin(); i != g->bc_ixs.end(); ++i) {
				get_building(i->ix).get_all_drawn_interior_verts(building_draw_interior); // interior
				get_building(i->ix).get_all_drawn_ext_wall_verts(building_draw_int_ext_walls); // interior of exterior walls
			}
		}
#if 0
		for (unsigned i = 0; i < tid_mapper.get_num_slots(); ++i) { // walls: 14269482 / 16553020 => 15887564 => 12184164 => 11461212
			unsigned const count(vert_counter.get_count(i));
			if (count == 0) continue;
			cout << i << ": R=" << count << " S=" << building_draw_interior.get_num_verts(tid_nm_pair_t(i)) << " C="
				 << building_draw_interior.get_cap_verts(tid_nm_pair_t(i)) << " " << get_texture_by_id(i).name << endl;
		}
#endif
		building_draw_interior     .finalize(grid_by_tile.size());
		building_draw_int_ext_walls.finalize(grid_by_tile.size());
	}
	void get_all_drawn_verts(bool is_tile) { // Note: non-const; building_draw is modified
		if (buildings.empty()) return;
		//timer_t timer("Get Building Verts"); // 140/670
		int const num_passes(is_tile ? 2 : 3); // skip interior pass for tiles; these verts will be generated later when they're needed for drawing

#pragma omp parallel for schedule(static) num_threads(num_passes)
		for (int pass = 0; pass < num_passes; ++pass) { // parallel loop doesn't help much because pass 0 takes most of the time
			if (pass == 0) { // exterior pass
				building_draw_vbo.clear();

				for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
					building_draw_vbo.cur_tile_id = (g - grid_by_tile.begin());
					for (auto i = g->bc_ixs.begin(); i != g->bc_ixs.end(); ++i) {get_building(i->ix).get_all_drawn_exterior_verts(building_draw_vbo);} // exterior
				}
				// disable shadows for materials that don't need them
				building_draw_vbo.set_no_shadows_for_tex(tid_nm_pair_t(NO_SHADOW_WHITE_TEX)); // for roof and solar panel trim
				building_draw_vbo.set_no_shadows_for_tex(building_texture_mgr.get_helipad_tid());
				//building_draw_vbo.set_no_shadows_for_tex(building_texture_mgr.get_solarp_tid()); // should solar panels cast shadows? they're pretty thin for shadow casters
				building_draw_vbo.set_no_shadows_for_tex(building_texture_mgr.get_bdoor2_tid());
				building_draw_vbo.finalize(grid_by_tile.size());
			}
			else if (pass == 1) { // windows pass
				get_all_window_verts(building_draw_windows, 0);
				if (is_night(WIND_LIGHT_ON_RAND)) {get_all_window_verts(building_draw_wind_lights, 1);} // only generate window verts at night
			}
			else if (pass == 2) { // interior pass; skip for is_tile case
				get_interior_drawn_verts();
			}
		} // for pass
	}
	void update_mem_usage(bool is_tile) {
		unsigned const num_everts(building_draw_vbo.num_verts() + building_draw_windows.num_verts() + building_draw_wind_lights.num_verts());
		unsigned const num_etris( building_draw_vbo.num_tris () + building_draw_windows.num_tris () + building_draw_wind_lights.num_tris ());
		unsigned const num_iverts(building_draw_interior.num_verts() + building_draw_int_ext_walls.num_verts());
		unsigned const num_itris( building_draw_interior.num_tris () + building_draw_int_ext_walls.num_tris ());
		gpu_mem_usage += (num_everts + num_iverts)*sizeof(vert_norm_comp_tc_color);
		
		if (!is_tile && num_everts > 0) {
			cout << "Building V: " << num_everts << ", T: " << num_etris << ", interior V: " << num_iverts << ", T: " << num_itris << ", mem: " << gpu_mem_usage << endl;
		}
	}
	void create_vbos(bool is_tile=0) {
		if (vbos_created) return; // already created
		vbos_created = 1;
		building_texture_mgr.check_windows_texture();
		tid_mapper.init();
		timer_t timer("Create Building VBOs", !is_tile);
		get_all_drawn_verts(is_tile);
		update_mem_usage(is_tile);
		building_draw_vbo          .upload_to_vbos();
		building_draw_windows      .upload_to_vbos();
		building_draw_wind_lights  .upload_to_vbos(); // Note: may be empty if not night time
		building_draw_interior     .upload_to_vbos();
		building_draw_int_ext_walls.upload_to_vbos();
	}
	void ensure_interior_geom_vbos() { // only for is_tile case
		if (!has_interior_geom) return; // no interior geom, nothing to do
		if (!building_draw_interior.empty()) return; // already created
		//timer_t timer("Create Building Interiors VBOs");
		get_interior_drawn_verts();
		update_mem_usage(1); // is_tile=1
		building_draw_interior.upload_to_vbos();
		building_draw_int_ext_walls.upload_to_vbos();
	}
	void ensure_window_lights_vbos() {
		if (!building_draw_wind_lights.empty()) return; // already calculated
		building_texture_mgr.check_windows_texture();
		get_all_window_verts(building_draw_wind_lights, 1);
		building_draw_wind_lights.upload_to_vbos();
	}
	void clear_vbos() {
		building_draw.clear_vbos();
		building_draw_vbo.clear_vbos();
		building_draw_windows.clear_vbos();
		building_draw_wind_lights.clear_vbos();
		building_draw_interior.clear_vbos();
		building_draw_int_ext_walls.clear_vbos();
		for (auto i = buildings.begin(); i != buildings.end(); ++i) {i->clear_room_geom();} // likely required for tiled buildings
		gpu_mem_usage = 0;
	}
	bool check_point_coll_xy(point const &pos) const { // Note: pos is in camera space
		if (empty()) return 0;
		vector3d const xlate(get_camera_coord_space_xlate());
		point const p1x(pos - xlate); // convert back to building space
		if (!range.contains_pt_xy(p1x)) return 0; // outside buildings bcube
		unsigned const gix(get_grid_ix(p1x));
		grid_elem_t const &ge(grid[gix]);
		if (ge.empty()) return 0; // skip empty grid
		if (!ge.bcube.contains_pt_xy(p1x)) return 0; // no intersection - skip this grid
		bool const xy_only = 1;
		point pos2(pos); // make a non-const copy
		point const p_last(pos);

		for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) {
			if (!b->contains_pt_xy(p1x)) continue;
			if (get_building(b->ix).check_sphere_coll(pos2, p_last, xlate, 0.0, xy_only)) return 1;
		}
		return check_road_seg_sphere_coll(ge, pos2, p_last, xlate, 0.0, xy_only, nullptr);
	}
	// Note: pos is in camera space; assumes only player collision queries set check_interior=1
	bool check_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm=nullptr, bool check_interior=0) const {
		if (empty()) return 0;
		bool const xy_only = 0;
		vector3d const xlate(get_camera_coord_space_xlate());
		cube_t bcube;
		bcube.set_from_sphere((pos - xlate), (radius + building_bcube_expand)); // expand to handle AC units, balconies, fire escapes, etc.
		bool saw_player_building(0);
		
		if (range.intersects_xy(bcube)) { // inside buildings bcube
			unsigned ixr[2][2];
			get_grid_range(bcube, ixr);
			float const dist(p2p_dist(pos, p_last));

			for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
				for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
					grid_elem_t const &ge(get_grid_elem(x, y));
					if (ge.empty()) continue; // skip empty grid
					if (!sphere_cube_intersect(pos, (radius + dist), (ge.bcube + xlate))) continue; // Note: makes little difference

					// Note: assumes buildings are separated so that only one sphere collision can occur
					for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) {
						if (!b->intersects_xy(bcube)) continue;
						building_t const &building(get_building(b->ix));
						if (building.check_sphere_coll(pos, p_last, xlate, radius, xy_only, cnorm, check_interior)) return 1;
						saw_player_building |= (check_interior && &building == player_building);
					} // for b
					if (check_interior && player_in_basement == 3) continue; // hack to keep player from popping from extended basement to top of driveway
					if (check_road_seg_sphere_coll(ge, pos, p_last, xlate, radius, xy_only, cnorm)) return 1; // check driveways
				} // for x
			} // for y
		}
		// hack to handle player in extended basement, which may be outside the building or even grid bbox:
		// if player is in the basement, and we haven't checked the player's building, and the player's building is in our range of buildings, check it now
		if (check_interior && !saw_player_building && player_in_basement && own_this_building(player_building)) {
			if (player_building->check_sphere_coll(pos, p_last, xlate, radius, xy_only, cnorm, check_interior)) return 1;
		}
		return 0;
	}
	// used for extended basement intersection checks; Note: bcube is in local building space
	bool check_cube_coll(cube_t const &bcube, bool xy_only, bool inc_basement, building_t const *exclude1, building_t const *exclude2) const {
		if (empty() || !range.intersects_xy(bcube)) return 0; // no buildings, or outside buildings bcube
		unsigned ixr[2][2];
		// buildings can extend outside grid bcubes, so we need to look in adjacent grids;
		// note that buildings are actually added to each grid they overlap, but that's their bcube only, and doesn't include their extended basement (which is added later);
		// extended basements are limited to the grid containing the building's center, but the grid bcube can extend outside the grid itself
		// due to other buildings that extend off the grid, even if the current building is completely contained and isn't itself in the adjacent grid
		bool const expand_by_one(inc_basement); // example: (-1.12, -15.7)
		get_grid_range(bcube, ixr, expand_by_one);

		// Note: can't check driveways/road_segs because they may not have been created yet
		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.empty()) continue; // skip empty grid
				if (!bcube.intersects_xy(ge.bcube)) continue; // Note: no need to check z-range

				for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) {
					building_t const &building(get_building(b->ix));
					if (&building == exclude1 || &building == exclude2) continue;
					if (inc_basement && building.cube_intersects_extb_room(bcube)) return 1; // extended basement intersection
					if (!bcube.intersects_xy(*b)) continue; // no intersection
						
					if (!xy_only) {
						if (bcube.z1() >= b->z2()) continue; // above the building, doesn't intersect (I guess attics are skipped?)
						// if parts has been allocated, the basement should be known, and the building's z1 should be valid; otherwise, extend the bcube z1 down
						float basement_z1(building.bcube.z1());
						
						if (!building.parts_generated) { // 1 basement level for house, 1+ for office
							basement_z1 -= (building.is_house ? 1 : global_building_params.max_office_basement_floors)*building.get_window_vspace();
						}
						if (bcube.z2() <= basement_z1) continue;
					}
					if (!building.parts_generated) return 1; // parts not yet allocated, assume it intersects

					for (auto p = building.parts.begin(); p != building.get_real_parts_end_inc_sec(); ++p) {
						if (xy_only ? p->intersects_xy(bcube) : p->intersects(bcube)) return 1; // intersects this part
					}
				} // for b
			} // for x
		} // for y
		return 0;
	}
	bool check_road_seg_sphere_coll(grid_elem_t const &ge, point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only, vector3d *cnorm) const {
		for (auto r = ge.road_segs.begin(); r != ge.road_segs.end(); ++r) {
			cube_t const cube(*r + xlate); // convert to camera space to agree with pos
			if (!cube.contains_pt_xy(pos) || (pos.z - radius) > cube.z2()) continue; // no collision - test top surface only (even when xy_only=1)
			pos.z = cube.z2() + radius;
			if (cnorm) {*cnorm = plus_z;}
			return 1;
		}
		return 0;
	}
	// Note: region is in building space, out is in camera space; no building rotation applied
	void get_road_segs_in_region(cube_t const &region, vect_cube_t &out) const {
		if (empty()) return;
		if (!range.intersects_xy(region)) return; // outside buildings bcube
		vector3d const xlate(get_camera_coord_space_xlate());
		unsigned ixr[2][2];
		get_grid_range(region, ixr);

		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.road_segs.empty() || !ge.bcube.intersects(region)) continue; // skip empty or non-intersecting grids

				for (cube_t const &c : ge.road_segs) {
					if (c.intersects_xy(region)) {out.push_back(c + xlate);}
				}
			} // for x
		} // for y
	}

	// Note: p1 and p2 are in building space; returns type of surface that was hit
	unsigned check_line_coll(point const &p1, point const &p2, float &t, unsigned &hit_bix, bool ret_any_pt, bool no_coll_pt, bool check_non_coll=0) const {
		if (empty()) return 0;

		if (p1.x == p2.x && p1.y == p2.y) { // vertical line special case optimization (for example map mode)
			if (!get_bcube().contains_pt_xy(p1)) return 0;
			unsigned const gix(get_grid_ix(p1));
			grid_elem_t const &ge(grid[gix]);
			if (ge.bc_ixs.empty()) return 0; // skip empty grid
			if (!ge.bcube.contains_pt_xy(p1)) return 0; // no intersection - skip this grid

			for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) {
				if (!b->contains_pt_xy(p1)) continue;
				unsigned const ret(get_building(b->ix).check_line_coll(p1, p2, t, 0, ret_any_pt, no_coll_pt, check_non_coll));
				if (ret) {hit_bix = b->ix; return ret;} // can only intersect one building
			} // for b
			for (cube_t const &seg : ge.road_segs) { // check driveways; not guaranteed to be correct if driveway is in another grid than building?
				if (seg.contains_pt_xy(p1)) return BLDG_COLL_DRIVEWAY;
			}
			return 0; // no coll
		}
		cube_t bcube(p1, p2);
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		unsigned coll(BLDG_COLL_NONE);
		point end_pos(p2);

		// for now, just do a slow iteration over every grid element within the line's bbox in XY
		// Note: should probably iterate over the grid in XY order from the start to the end of the line, or better yet use a line drawing algorithm
		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.bc_ixs.empty()) continue; // skip empty grid
				if (!check_line_clip(p1, end_pos, ge.bcube.d)) continue; // no intersection - skip this grid

				for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) { // Note: okay to check the same building more than once
					if (!b->intersects(bcube)) continue;
					float t_new(t);
					unsigned const ret(get_building(b->ix).check_line_coll(p1, p2, t_new, 0, ret_any_pt, no_coll_pt, check_non_coll));

					if (ret && t_new <= t) { // closer hit pos, update state
						t = t_new; hit_bix = b->ix; coll = ret;
						end_pos = p1 + t*(p2 - p1);
						if (ret_any_pt) return coll;
					}
				} // for b
			} // for x
		} // for y
		return coll;
	}

	// Note: p1 and p2 are in building space
	void update_zmax_for_line(point const &p1_in, point const &p2_in, float radius, float house_extra_zval, float &cur_zmax) const {
		if (empty()) return;
		point p1(p1_in.x, p1_in.y, cur_zmax), p2(p2_in.x, p2_in.y, cur_zmax); // use current/starting zmax for line end points
		cube_t bcube(p1, p2);
		bcube.expand_by_xy(radius);
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		
		// for now, just do a slow iteration over every grid element within the line's bbox in XY
		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.bc_ixs.empty() || (ge.bcube.z2() + house_extra_zval) < cur_zmax) continue; // skip empty grid or grids below the line
				cube_t grid_bc(ge.bcube);
				grid_bc.z2() += house_extra_zval; // assume it could be a house
				grid_bc.expand_by_xy(radius);
				if (!check_line_clip(p1, p2, grid_bc.d)) continue; // no intersection - skip this grid

				for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) { // Note: okay to check the same building more than once
					cube_t bbc(*b);
					building_t const &building(get_building(b->ix));
					if (building.is_house) {bbc.z2() += house_extra_zval;}
					if (bbc.z2() < cur_zmax) continue; // building is not tall enough
					bbc.expand_by_xy(radius);
					if (!bcube.intersects(bbc) || !check_line_clip(p1, p2, bbc.d)) continue;
					float t(1.0); // result is unused
					// Note: unclear if we need to do a detailed line collision check, maybe testing the building bbox is good enough?
					// if radius is passed in as nonzero, then simply assume it intersects because check_line_coll() can't easily take a radius value
					if (radius > 0.0 || building.check_line_coll(p1, p2, t, 0, 1, 1)) {
						max_eq(cur_zmax, bbc.z2());
						p1.z = p2.z = cur_zmax; // update line end points to match this new elevation so that line clipping works properly
					}
				} // for b
			} // for x
		} // for y
	}

	// Note: we can get building_id by calling check_ped_coll() or get_building_bcube_at_pos(); p1 and p2 are in building space
	bool check_line_coll_building(point const &p1, point const &p2, unsigned building_id) const { // Note: not thread safe due to static points
		assert(building_id < buildings.size());
		float t_new(1.0);
		return buildings[building_id].check_line_coll(p1, p2, t_new, 0, 1, 1); // occlusion_only=0, ret_any_pt=1, no_coll_pt=1
	}
	bool check_sphere_coll_building(point const &pos, float radius, bool xy_only, unsigned building_id) const {
		assert(building_id < buildings.size());
		return buildings[building_id].check_sphere_coll(pos, radius, xy_only);
	}
	bool check_building_point_or_cylin_contained(point const &pos, float radius, bool inc_details, unsigned building_id) const {
		static vector<point> points; // reused across calls
		assert(building_id < buildings.size());
		return buildings[building_id].check_point_or_cylin_contained(pos, radius, points, 0, 0, 0, inc_details); // attic=0, extb=0, roof=0
	}

	int get_building_bcube_contains_pos(point const &pos) { // Note: not thread safe due to static points
		if (empty()) return -1;
		unsigned const gix(get_grid_ix(pos));
		grid_elem_t const &ge(grid[gix]);
		if (ge.bc_ixs.empty() || !ge.bcube.contains_pt(pos)) return -1; // skip empty or non-containing grid
		static vector<point> points; // reused across calls

		for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) {
			if (b->contains_pt(pos)) {return b->ix;} // found
		}
		return -1;
	}
	cube_t get_grid_bcube_for_building(building_t const &b) const {
		if (empty() || !own_this_building(&b)) return cube_t();
		return grid[get_grid_ix(b.bcube.get_cube_center())].bcube;
	}

	// Note: not thread safe due to static points
	// return value: 0=no cont, 1=part, 2=attic, 3=ext basement, 4=roof access, 5=detail
	int check_ped_coll(point const &pos, float bcube_radius, float detail_radius, unsigned plot_id, unsigned &building_id, cube_t *coll_cube) const {
		if (empty()) return 0;
		assert(plot_id < bix_by_plot.size());
		vector<unsigned> const &bixes(bix_by_plot[plot_id]); // should be populated in gen()
		if (bixes.empty()) return 0;
		cube_t bcube; bcube.set_from_sphere(pos, bcube_radius);
		static vector<point> points; // reused across calls

		// Note: assumes buildings are separated so that only one ped collision can occur
		for (auto b = bixes.begin(); b != bixes.end(); ++b) {
			building_t const &building(get_building(*b));
			if (building.bcube.x1() > bcube.x2())     break; // no further buildings can intersect (sorted by x1)
			if (!building.bcube.intersects_xy(bcube)) continue;
			// inc_attic=0, inc_ext_basement=0, inc_roof_acc=0, inc_details=1
			int const ret(building.check_point_or_cylin_contained(pos, detail_radius, points, 0, 0, 0, 1, coll_cube));
			if (ret) {building_id = *b; return ret;}
		}
		return 0;
	}
	bool select_building_in_plot(unsigned plot_id, unsigned rand_val, unsigned &building_id) const {
		if (bix_by_plot.empty()) return 0; // not setup / no buildings
		assert(plot_id < bix_by_plot.size());
		vector<unsigned> const &bixes(bix_by_plot[plot_id]);
		if (bixes.empty()) return 0;
		building_id = bixes[rand_val % bixes.size()];
		return 1;
	}

	// Note: called on init, don't need to use get_camera_coord_space_xlate()
	template<typename RET> void query_for_cube(cube_t const &query_cube, vector<RET> &cubes, int query_mode) const {
		if (empty()) return; // nothing to do
		unsigned ixr[2][2];
		get_grid_range(query_cube, ixr);

		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.bc_ixs.empty() || !query_cube.intersects_xy(ge.bcube)) continue;

				for (cube_with_ix_t const &b : ge.bc_ixs) {
					if (!query_cube.intersects_xy(b)) continue;
					if (get_grid_ix(query_cube.get_llc().max(b.get_llc())) != (y*grid_sz + x)) continue; // add only if in home grid (to avoid duplicates)
					if      (query_mode == 0) {cubes.push_back(b);} // return building bcube
					else if (query_mode == 1) { // return house driveways
						building_t const &B(get_building(b.ix));
						// B.city_driveway will be the added driveway; if a valid driveway already connects the garage to the road, it will remain all zeros
						if (B.maybe_add_house_driveway(query_cube, b.ix) && !B.city_driveway.is_all_zeros()) {cubes.push_back(B.city_driveway);} // ix not set
					}
					else {assert(0);} // invalid mode/not implemented
				}
			} // for x
		} // for y
	}
	void get_overlapping_bcubes      (cube_t const &xy_range, vect_cube_t         &bcubes) const {return query_for_cube(xy_range, bcubes,    0);}
	void get_overlapping_bcubes      (cube_t const &xy_range, vect_cube_with_ix_t &bcubes) const {return query_for_cube(xy_range, bcubes,    0);}
	void add_house_driveways_for_plot(cube_t const &plot,     vect_cube_t      &driveways) const {return query_for_cube(plot,     driveways, 1);}

	void get_building_ext_basement_bcubes(cube_t const &city_bcube, vect_cube_t &bcubes) const {
		vector<cube_with_ix_t> cand_bldgs;
		get_overlapping_bcubes(city_bcube, cand_bldgs);

		for (cube_with_ix_t const &b : cand_bldgs) {
			building_t const &building(get_building(b.ix));
			if (!building.has_ext_basement() || !building.interior) continue;
			for (auto r = building.interior->ext_basement_rooms_start(); r != building.interior->rooms.end(); ++r) {bcubes.push_back(*r);}
			for (tunnel_seg_t const &t : building.interior->tunnels) {bcubes.push_back(t.bcube);}
		}
	}

	// walkways
private:
	float get_walkway_buildings_and_max_sz(cube_t const &city_bcube, vector<cube_with_ix_t> &city_bldgs, vector<cube_with_ix_t> &ww_bldgs) {
		get_overlapping_bcubes(city_bcube, city_bldgs);
		float max_xy_sz(0.0);

		for (cube_with_ix_t const &b : city_bldgs) {
			building_t const &building(get_building(b.ix));
			if (building.is_house || !building.is_cube() || building.is_rotated()) continue; // walkways not supported for this building
			ww_bldgs.push_back(b);
			max_eq(max_xy_sz, max(building.bcube.dx(), building.bcube.dy()));
		}
		return max_xy_sz;
	}
	bool check_if_blocked_by_building(cube_t const &cand, vect_cube_with_ix_t const &city_bldgs, unsigned bix1, unsigned bix2) const {
		for (cube_with_ix_t const &b : city_bldgs) {
			if (b.ix == bix1 || b.ix == bix2 || !b.intersects_xy(cand)) continue;
			building_t const &building(get_building(b.ix));
			if (b.z2() > cand.z1() && building.cube_int_parts_no_sec(cand)) return 1;
			if (building.has_helipad && building.get_helipad_bcube().intersects_xy(cand)) return 1; // check for crossing above helipad as well
		}
		return 0;
	}
	int choose_rand_walkway_side_mat_ix(rand_gen_t &rgen) const {
		if (global_building_params.mat_gen_ix_city.empty()) return -1; // should never fail?

		for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to choose a valid walkway material
			unsigned const mat_ix(global_building_params.choose_rand_mat(rgen, 1, 0, 0)); // city_only=1, non_city_only=0, residential=0
			if (!global_building_params.get_material(mat_ix).no_walkways) return mat_ix;
		}
		return -1; // not found
	}
public:
	void connect_buildings_with_walkways(cube_t const &city_bcube) {
		vector<cube_with_ix_t> city_bldgs, ww_bldgs;
		float const max_xy_sz(get_walkway_buildings_and_max_sz(city_bcube, city_bldgs, ww_bldgs));
		if (ww_bldgs.size() < 2) return; // no buildings to connect
		float const max_walkway_len(1.5*max_xy_sz), road_width(get_road_max_width()), pp_height(get_power_pole_height());
		vect_cube_t blocked; // walkways currently placed for this city
		rand_gen_t rgen;

		for (auto i1 = ww_bldgs.begin(); i1 != ww_bldgs.end(); ++i1) {
			building_t &b1(get_building(i1->ix));
			float const min_ww_width(2.0*b1.get_office_ext_doorway_width()), floor_spacing(b1.get_window_vspace()); // should be the same for all buildings
			float const bot_z_add((DRAW_CITY_INT_WINDOWS ? 0.125 : 0.25)*floor_spacing); // reduce if there are interior windows so that we don't see inside the walkway bottom
			float const power_pole_clearance(1.25*bot_z_add);
			unsigned const min_floors_above_power_pole(unsigned((pp_height + power_pole_clearance)/floor_spacing) + 1U); // for crossing roads; take ceil
			float const walkway_zmin_short(b1.ground_floor_z1 + ((bot_z_add > 0.0) ? 2.0 : 1.0)*floor_spacing); // two floors up, to account for bot_z_add
			float const walkway_zmin_long (b1.ground_floor_z1 + min_floors_above_power_pole    *floor_spacing); // N floors up
			float const edge_pad(0.25*b1.get_wall_thickness());

			for (auto i2 = i1+1; i2 != ww_bldgs.end(); ++i2) {
				building_t &b2(get_building(i2->ix));
				assert(!b1.bcube.intersects_xy(b2.bcube)); // sanity check
				assert(b1.ground_floor_z1 == b2.ground_floor_z1); // must be at the same elevation
				if (b2.get_window_vspace() != floor_spacing) continue; // floor spacing differs, can't connect (optional, can instead not connect interiors)
				bool connected(0);

				for (unsigned dim = 0; dim < 2; ++dim) { // connection dim
					bool const dir(b1.bcube.get_center_dim(dim) < b2.bcube.get_center_dim(dim)); // dir of b2 relative to b1: 0=to the left, 1=to the right
					// first, check that the bcubes have a large enough projection and small enough gap
					float const lo(max(b1.bcube.d[!dim][0], b2.bcube.d[!dim][0]) + edge_pad), hi(min(b1.bcube.d[!dim][1], b2.bcube.d[!dim][1]) - edge_pad); // projection range
					if (hi - lo < min_ww_width)   continue; // projection too small
					float const length(fabs(b1.bcube.d[dim][dir] - b2.bcube.d[dim][!dir]));
					if (length > max_walkway_len) continue; // buildings are too far apart
					// check for other buildings in between that block this walkway
					cube_t cand;
					cand.d[ dim][ dir] = b2.bcube.d[dim][!dir];
					cand.d[ dim][!dir] = b1.bcube.d[dim][ dir];
					cand.d[!dim][0] = lo; cand.d[!dim][1] = hi;
					set_cube_zvals(cand, walkway_zmin_short, min(b1.bcube.z2(), b2.bcube.z2()));
					if (check_if_blocked_by_building(cand, city_bldgs, i1->ix, i2->ix)) continue; // Note: uses *all* buildings

					for (auto P1 = b1.parts.begin(); P1 != b1.parts.end(); ++P1) {
						for (auto P2 = b2.parts.begin(); P2 != b2.parts.end(); ++P2) {
							cube_t const &p1(*P1), &p2(*P2);
							float const lo(max(p1.d[!dim][0], p2.d[!dim][0]) + edge_pad), hi(min(p1.d[!dim][1], p2.d[!dim][1]) - edge_pad), width(hi - lo); // projection range
							if (width < min_ww_width)     continue; // projection too small
							float const length(fabs(p1.d[dim][dir] - p2.d[dim][!dir]));
							if (length > max_walkway_len) continue; // buildings are too far apart
							// we can't easily check if the walkway crosses a road since we have neither the roads nor the plots here, so be conservative and check for min length
							bool const is_long(length > road_width);
							float const walkway_zmin(is_long ? walkway_zmin_long : walkway_zmin_short);
							float const z2_min_test(walkway_zmin + 0.5*floor_spacing); // test z2 against center of next floor
							if (p1.z2() < z2_min_test || p2.z2() < z2_min_test) continue; // too short
							float const zlo(max(walkway_zmin, max(p1.z1(), p2.z1()))), zhi(min(p1.z2(), p2.z2()));
							if (zhi - zlo < 0.9*floor_spacing) continue; // no overlap of at least a floor (give or take)
							cube_t walkway;
							walkway.d[ dim][ dir] = p2.d[dim][!dir];
							walkway.d[ dim][!dir] = p1.d[dim][ dir];
							walkway.d[!dim][0] = lo; walkway.d[!dim][1] = hi;
							float const target_width(min_ww_width*rgen.rand_uniform(1.0, 1.8));
							if (width > target_width) {walkway.expand_in_dim(!dim, -0.5*(width - target_width));} // shrink the width if needed
							set_cube_zvals(walkway, zlo, zhi);
							unsigned num_floors_max(1 + (rgen.rand()%3)); // 1-3 floors
							float const z2_max(zlo + num_floors_max*floor_spacing); // limit height

							if (z2_max < walkway.z2()) { // reduce walkway height
								unsigned const num_floors_above(round_fp((walkway.z2() - z2_max)/floor_spacing));
								walkway.z2() = z2_max;

								if (num_floors_above > 1) { // reduced by at least one floor
									unsigned const num_floors_raise(rgen.rand() % num_floors_above);
									if (num_floors_raise > 0) {walkway.translate_dim(2, num_floors_raise*floor_spacing);} // raise it up
								}
							}
							assert(walkway.dz() > 0.0);
							cube_t const walkway_interior(walkway); // capture before applying bot_z_add
							walkway.z1() -= bot_z_add; // add extra space at the bottom for support; can't add to the top in case we're at the top building floor
							assert(walkway.is_strictly_normalized());
							// check for other parts or walkways blocking the walkway
							if (b1.cube_int_parts_no_sec(walkway) || b2.cube_int_parts_no_sec(walkway)) continue;
							bool ww_blocked(0);
							for (cube_t const &w : blocked) {ww_blocked |= w.intersects(walkway);} // check other walkways
							if (ww_blocked) continue;
							bool const mat1_valid(!b1.get_material().no_walkways), mat2_valid(!b2.get_material().no_walkways);
							bool owner_is_b1(0);
							int side_mat_ix(-1);
							if (mat1_valid && mat2_valid) {owner_is_b1 = rgen.rand_bool();} // both are value, choose a building randomly
							else if (mat1_valid) {owner_is_b1 = 1;}
							else if (mat2_valid) {owner_is_b1 = 0;}
							else { // neither is valid
								owner_is_b1 = rgen.rand_bool(); // choose a random building to get the roof and side color from
								side_mat_ix = choose_rand_walkway_side_mat_ix(rgen);
							}
							building_t const &ww_owner(owner_is_b1 ? b1 : b2);
							if (side_mat_ix < 0) {side_mat_ix = ww_owner.mat_ix;} // if side_mat_ix wasn't set above, use the parent building's material
							all_walkways.emplace_back(walkway, dim, side_mat_ix, ww_owner.mat_ix, ww_owner.side_color, ww_owner.roof_color, floor_spacing);
							blocked.push_back(walkway);
							building_walkway_geom_t bwg(walkway_interior, dim);

							if (ADD_WALKWAY_EXT_DOORS) { // add exterior doors connected to walkways
								b1.add_walkway_door(bwg,  dir, (P1 - b1.parts.begin()));
								b2.add_walkway_door(bwg, !dir, (P2 - b2.parts.begin()));
							}
							b1.walkways.emplace_back(bwg,  owner_is_b1, &b2);
							b2.walkways.emplace_back(bwg, !owner_is_b1, &b1);
							connected = 1;
							break; // only need one connection
						} // for p2
						if (connected) break;
					} // for p1
					if (connected) break;
				} // for dim
			} // for i2
		} // for i1
	}
	void get_walkways_for_city(cube_t const &city_bcube, vect_bldg_walkway_t &walkways) const {
		for (bldg_walkway_t const &w : all_walkways) {
			if (city_bcube.contains_cube_xy(w)) {walkways.push_back(w);}
		}
	}
private:
	struct walkway_cand_t {
		cube_t bcube;
		unsigned bix, pix;
		bool dir;
		walkway_cand_t(cube_t const &bc, unsigned b, unsigned p, bool d) : bcube(bc), bix(b), pix(p), dir(d) {}
	};
public:
	bool connect_buildings_to_skyway(cube_t &m_bcube, bool m_dim, cube_t const &city_bcube, vector<skyway_conn_t> &ww_conns) {
		vector<cube_with_ix_t> city_bldgs, ww_bldgs;
		float const max_xy_sz(get_walkway_buildings_and_max_sz(city_bcube, city_bldgs, ww_bldgs)), max_walkway_len(1.5*max_xy_sz);
		if (ww_bldgs.size() < 2) return 0; // not enough buildings
		bool const conn_dim(!m_dim);
		float const centerline(m_bcube.get_center_dim(conn_dim));
		cube_t conn_area(m_bcube);
		conn_area.expand_by_xy(max_walkway_len);
		cube_t all_conn_bc;
		vector<walkway_cand_t> cands;
		float ww_zmin(m_bcube.z2());

		for (auto i = ww_bldgs.begin(); i != ww_bldgs.end(); ++i) {
			building_t &b(get_building(i->ix));
			assert(!b.bcube.intersects_xy(m_bcube)); // sanity check
			if (!b.bcube.intersects_xy(conn_area)) continue; // too far from skyway
			float const min_ww_width(2.0*b.get_office_ext_doorway_width()), floor_spacing(b.get_window_vspace()); // should be the same for all buildings
			if (b.bcube.z2() < m_bcube.z1() + floor_spacing) continue; // too short to connect to skyway

			for (auto P = b.parts.begin(); P != b.parts.end(); ++P) {
				cube_t const &p(*P);
				if (!p.intersects_xy(conn_area))           continue; // too far from skyway
				if (p.z2() < m_bcube.z1() + floor_spacing) continue; // too short to connect to skyway
				if (p.z1() > m_bcube.z1())                 continue; // starts above skyway
				bool const dir(centerline < p.get_center_dim(conn_dim));
				cube_t conn_area(p);
				conn_area.expand_in_dim(!conn_dim, -0.25*b.get_wall_thickness()); // shrink slightly to prevent Z-fighting with inside edges of parts
				conn_area.d[conn_dim][ dir] = p      .d[conn_dim][!dir]; // flush with part
				conn_area.d[conn_dim][!dir] = m_bcube.d[conn_dim][ dir]; // connect to skyway
				// clamp to shared range in the skyway dim; really should not change the part width since the skyway should run the entire length of the city
				max_eq(conn_area.d[!conn_dim][0], m_bcube.d[!conn_dim][0]);
				min_eq(conn_area.d[!conn_dim][1], m_bcube.d[!conn_dim][1]);
				float const width(conn_area.get_sz_dim(!conn_dim));
				if (width < min_ww_width) continue; // too narrow
				unsigned const floor_ix(ceil((m_bcube.z1() - p.z1())/floor_spacing)); // round up
				conn_area.z1() = p.z1() + floor_ix*floor_spacing;
				conn_area.z2() = conn_area.z1()  + floor_spacing; // one floor in height
				if (conn_area.z2() > p.z2()) continue; // too high
				float const target_width(min_ww_width*rgen.rand_uniform(1.0, 1.8)), space_to_sides(0.5*(width - target_width));
				if (space_to_sides > 0.0) {conn_area.expand_in_dim(!conn_dim, -space_to_sides);} // shrink the width if needed
				assert(conn_area.is_strictly_normalized());
				// first try placing centered on the part; if there's extra space to the sides, and centered fails, then try to a random side, then the other side
				float const first_xlate((rgen.rand_bool() ? 1.0 : -1.0)*space_to_sides), xlate_vals[3] = {0.0f, first_xlate, -2.0f*first_xlate};
				bool success(0);

				for (unsigned n = 0; n < ((space_to_sides > 0.0) ? 3U : 1U); ++n) {
					conn_area.translate_dim(!conn_dim, xlate_vals[n]);
					if (check_if_blocked_by_building(conn_area, city_bldgs, i->ix, i->ix)) continue; // Note: uses *all* buildings
					if (b.cube_int_parts_no_sec(conn_area)) continue; // intersects another part
					cands.emplace_back(conn_area, i->ix, (P - b.parts.begin()), dir);
					all_conn_bc.assign_or_union_with_cube(conn_area);
					min_eq(ww_zmin, conn_area.z1());
					success = 1;
					break;
				} // for n
				if (success) break; // only connect the first valid part
			} // for p
		} // for i
		if (cands.size() < 2) return 0; // need at least two connected buildings
		if (all_conn_bc.get_sz_dim(m_dim) < 0.5*m_bcube.get_sz_dim(m_dim)) return 0; // less than half the length is connected: fail
		float const ww_conn_width(2.5*m_bcube.get_sz_dim(conn_dim));
		for (unsigned d = 0; d < 2; ++d) {m_bcube.d[m_dim][d] = all_conn_bc.d[m_dim][d];} // clip to shared connection sub-length
		m_bcube.expand_in_dim(m_dim, 0.25*m_bcube.get_sz_dim(conn_dim)); // extend slightly
		m_bcube.translate_dim(2, (ww_zmin - m_bcube.z1() - 0.05*m_bcube.dz())); // translate to the bottom of the lowest walkway; walkways are often all the same zval
		ww_conns.reserve(cands.size());
		
		for (walkway_cand_t const &cand : cands) { // actually add walkways
			building_t &b(get_building(cand.bix));
			int side_mat_ix(choose_rand_walkway_side_mat_ix(rgen));
			if (side_mat_ix < 0) {side_mat_ix = b.mat_ix;} // if side_mat_ix wasn't set above, use the building's material
			all_walkways.emplace_back(cand.bcube, conn_dim, side_mat_ix, b.mat_ix, b.side_color, b.roof_color, b.get_window_vspace());
			building_walkway_geom_t bwg(cand.bcube, conn_dim);
			if (ADD_WALKWAY_EXT_DOORS) {b.add_walkway_door(bwg, !cand.dir, cand.pix);}
			b.walkways.emplace_back(bwg, 1, nullptr); // owned, no conn_bldg
			b.walkways.back().open_ends[!cand.dir] = all_walkways.back().open_ends[!cand.dir] = 1; // flag end connected to skyway as open
			cube_t skyway_conn(m_bcube);
			max_eq(skyway_conn.d[m_dim][0], cand.bcube.d[m_dim][0]-ww_conn_width);
			min_eq(skyway_conn.d[m_dim][1], cand.bcube.d[m_dim][1]+ww_conn_width);
			b.walkways.back().skyway_conn = skyway_conn;
			cube_t conn(cand.bcube);
			conn.d[conn_dim][cand.dir] = conn.d[conn_dim][!cand.dir];
			ww_conns.emplace_back(conn, conn_dim, !cand.dir, &b);
		} // for cand
		return 1; // success
	}

	void get_city_building_walkways(cube_t const &city_bcube, vector<building_walkway_t *> &bwws) {
		vector<cube_with_ix_t> city_bldgs;
		get_overlapping_bcubes(city_bcube, city_bldgs);

		for (cube_with_ix_t const &b : city_bldgs) {
			building_t &building(get_building(b.ix));
			
			for (building_walkway_t &ww : building.walkways) {
				if (ww.is_owner) {bwws.push_back(&ww);} // add walkway only once, for the building owner
			}
		}
	}

	void get_power_points(cube_t const &xy_range, vector<point> &ppts) const { // similar to above function, but returns points rather than cubes
		if (empty()) return; // nothing to do
		unsigned ixr[2][2];
		get_grid_range(xy_range, ixr);

		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.bc_ixs.empty() || !xy_range.intersects_xy(ge.bcube)) continue;

				for (auto b = ge.bc_ixs.begin(); b != ge.bc_ixs.end(); ++b) {
					if (!xy_range.intersects_xy(*b)) continue;
					if (get_grid_ix(xy_range.get_llc().max(b->get_llc())) != (y*grid_sz + x)) continue; // add only if in home grid (to avoid duplicates)
					get_building(b->ix).get_power_point(ppts);
				}
			} // for x
		} // for y
	}

	void get_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state, bool cur_building_only=0) const {
		state.init(pdu.pos, get_camera_coord_space_xlate());
		if (cur_building_only) return; // no grid/buildings iteration
		
		for (auto g = grid.begin(); g != grid.end(); ++g) {
			if (g->bc_ixs.empty()) continue;
			if (!building_grid_visible(state.xlate, g->bcube, pdu)) continue; // VFC; use exterior bcube; pass in our custom pdu with lowered near clip plane
			
			for (auto b = g->bc_ixs.begin(); b != g->bc_ixs.end(); ++b) {
				if ((int)b->ix == state.exclude_bix) continue; // excluded
				cube_t const c(*b + state.xlate); // check far clipping plane first because that's more likely to reject buildings
				
				// if player is inside this building, skip occlusion so that objects are visible through windows
				if (state.skip_cont_camera && !(player_in_basement || player_in_attic) && c.contains_pt(pdu.pos)) {
					building_t const &bldg(get_building(b->ix));
					if (bldg.has_int_windows() || bldg.point_near_ext_door((state.pos - state.xlate), get_door_open_dist())) continue;
				}
				if (dist_less_than(pdu.pos, c.closest_pt(pdu.pos), pdu.far_) && pdu.cube_visible(c)) {state.building_ids.push_back(*b);}
			} // for b
		} // for g
	}
	bool check_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t const &state) const { // pts are in building space
		point const pos_bs(state.pos - state.xlate);

		for (auto b = state.building_ids.begin(); b != state.building_ids.end(); ++b) {
			if ((int)b->ix == state.exclude_bix) continue;
			if (get_region(pos_bs, b->d) & get_region(pts[0], b->d)) continue; // line outside - early reject optimization
			if (!b->line_intersects(pos_bs, pts[0])) continue; // early reject optimization
			building_t const &building(get_building(b->ix));
			bool occluded(1);

			for (unsigned i = 0; i < npts; ++i) {
				float t(1.0); // start at end of line
				if (!building.check_line_coll(pos_bs, pts[i], t, 1)) {occluded = 0; break;}
			}
			if (occluded) return 1;
		} // for b
		return 0;
	}
	bool single_cube_visible_check(point const &pos, cube_t const &c) const {
		if (empty()) return 1;
		float const z(c.z2()); // top edge
		point const pts[4] = {point(c.x1(), c.y1(), z), point(c.x2(), c.y1(), z), point(c.x2(), c.y2(), z), point(c.x1(), c.y2(), z)};
		cube_t query_region(c);
		query_region.union_with_pt(pos);

		for (auto g = grid.begin(); g != grid.end(); ++g) {
			if (g->bc_ixs.empty() || !g->bcube.intersects(query_region)) continue;

			for (auto b = g->bc_ixs.begin(); b != g->bc_ixs.end(); ++b) {
				if (!b->intersects(query_region)) continue;
				building_t const &building(get_building(b->ix));
				bool occluded(1);

				for (unsigned i = 0; i < 4; ++i) {
					float t(1.0); // start at end of line
					if (!building.check_line_coll(pos, pts[i], t, 1)) {occluded = 0; break;}
				}
				if (occluded) return 0;
			} // for b
		} // for g
		return 1;
	}
}; // building_creator_t


class building_tiles_t {
	typedef pair<int, int> xy_pair;
	typedef map<xy_pair, building_creator_t> tile_map_t;
	tile_map_t tiles; // key is {x, y} pair
	//set<xy_pair> generated; // only used in heightmap terrain mode, and generally limited to the size of the heightmap in tiles
	vector3d max_extent;

	tile_map_t::const_iterator get_tile_by_pos_cs(point const &pos) const { // Note: pos is in camera space
		vector3d const xlate(get_camera_coord_space_xlate());
		int const x(round_fp(0.5f*(pos.x - xlate.x)/X_SCENE_SIZE)), y(round_fp(0.5f*(pos.y - xlate.y)/Y_SCENE_SIZE));
		return tiles.find(make_pair(x, y));
	}
	tile_map_t::const_iterator get_tile_by_pos_bs(point const &pos) const { // Note: pos is in building space
		int const x(round_fp(0.5f*pos.x/X_SCENE_SIZE)), y(round_fp(0.5f*pos.y/Y_SCENE_SIZE));
		return tiles.find(make_pair(x, y));
	}
public:
	building_tiles_t() : max_extent(zero_vector) {}
	bool     empty() const {return tiles.empty();}
	unsigned size()  const {return tiles.size();}
	vector3d get_max_extent() const {return max_extent;}

	int create_tile(int x, int y, bool allow_flatten) { // return value: 0=already exists, 1=newly generaged, 2=re-generated
		xy_pair const loc(x, y);
		auto it(tiles.find(loc));
		if (it != tiles.end()) return 0; // already exists
		//cout << "Create building tile " << x << "," << y << ", tiles: " << tiles.size() << endl; // 299 tiles
		building_creator_t &bc(tiles[make_pair(x, y)]); // insert it
		assert(bc.empty());
		int const border(allow_flatten ? 1 : 0); // add a 1 pixel border around the tile to avoid creating a seam when an adjacent tile's edge height is modified
		cube_t bcube(all_zeros);
		bcube.x1() = get_xval(x*MESH_X_SIZE + border);
		bcube.y1() = get_yval(y*MESH_Y_SIZE + border);
		bcube.x2() = get_xval((x+1)*MESH_X_SIZE - border);
		bcube.y2() = get_yval((y+1)*MESH_Y_SIZE - border);
		global_building_params.set_pos_range(bcube);
		int const rseed(x + (y << 16) + 12345); // should not be zero
		bc.gen(global_building_params, 0, have_cities(), 1, allow_flatten, rseed); // if there are cities, then tiles are non-city/secondary buildings
		global_building_params.restore_prev_pos_range();
		max_extent = max_extent.max(bc.get_max_extent());
		//if (allow_flatten) {return (generated.insert(loc).second ? 1 : 2);} // Note: caller no longer uses this value, so don't need to maintain generated
		return 1;
	}
	bool remove_tile(int x, int y) {
		auto it(tiles.find(make_pair(x, y)));
		if (it == tiles.end()) return 0; // not found
		//cout << "Remove building tile " << x << "," << y << ", tiles: " << tiles.size() << endl;
		it->second.clear_vbos(); // free VBOs/VAOs
		tiles.erase(it);
		return 1;
	}
	void clear_vbos() {
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {i->second.clear_vbos();}
	}
	void clear() {
		clear_vbos();
		tiles.clear();
	}
	bool check_point_coll_xy(point const &pos) const { // Note: pos is in camera space
		if (empty()) return 0;
		auto it(get_tile_by_pos_cs(pos)); // single point, use map lookup optimization (for example for grass)
		if (it == tiles.end()) return 0;
		return it->second.check_point_coll_xy(pos);
	}
	bool check_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm=nullptr, bool check_interior=0) const { // Note: pos is in camera space
		if (empty()) return 0;

		for (auto const &t : tiles) {
			if (t.second.check_sphere_coll(pos, p_last, radius, cnorm, check_interior)) return 1;
		}
		return 0;
	}
	bool check_cube_coll(cube_t const &bcube, bool xy_only, bool inc_basement, building_t const *exclude1, building_t const *exclude2) const {
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {
			if (i->second.check_cube_coll(bcube, xy_only, inc_basement, exclude1, exclude2)) return 1;
		}
		return 0;
	}
	void get_road_segs_in_region(cube_t const &region, vect_cube_t &out) const {
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {i->second.get_road_segs_in_region(region, out);}
	}
	bool get_building_hit_color(point const &p1, point const &p2, colorRGBA &color) const { // Note: p1/p2 are in building space
		if (empty()) return 0;

		if (p1.x == p2.x && p1.y == p2.y) { // vertical line, use map lookup optimization (for overhead map mode)
			auto it(get_tile_by_pos_bs(p1));
			if (it == tiles.end()) return 0;
			return it->second.get_building_hit_color(p1, p2, color);
		}
		cube_t const line_bcube(p1, p2);

		for (auto i = tiles.begin(); i != tiles.end(); ++i) {
			if (!i->second.get_bcube().intersects(line_bcube)) continue; // optimization
			if (i->second.get_building_hit_color(p1, p2, color)) return 1; // line is generally pointed down and can only intersect one building; return the first hit
		}
		return 0;
	}
	unsigned check_line_coll(point const &p1, point const &p2, float &t, bool ret_any_pt, bool no_coll_pt) const {
		if (empty()) return 0;
		cube_t const line_bcube(p1, p2);
		unsigned ret(BLDG_COLL_NONE), hit_bix(0); // internal tile index, can't return

		for (auto i = tiles.begin(); i != tiles.end(); ++i) { // no vertical line test
			if (!i->second.get_bcube().intersects(line_bcube)) continue; // optimization
			unsigned const tile_ret(i->second.check_line_coll(p1, p2, t, hit_bix, ret_any_pt, no_coll_pt));
			if (tile_ret) {ret = tile_ret;}
		}
		return ret;
	}
	void update_zmax_for_line(point const &p1, point const &p2, float radius, float house_extra_zval, float &cur_zmax) const { // Note p1/p2 are in building space
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {
			if (!check_line_clip(p1, p2, i->second.get_bcube().d)) continue; // optimization
			i->second.update_zmax_for_line(p1, p2, radius, house_extra_zval, cur_zmax);
		}
	}
	cube_t get_grid_bcube_for_building(building_t const &b) const {
		auto it(get_tile_by_pos_bs(b.bcube.get_cube_center()));
		return ((it == tiles.end()) ? cube_t() : it->second.get_grid_bcube_for_building(b));
	}
	void add_drawn(vector3d const &xlate, vector<building_creator_t *> &bcs) {
		float const draw_dist(get_draw_tile_dist());
		point const camera(get_camera_pos() - xlate);

		for (auto i = tiles.begin(); i != tiles.end(); ++i) {
			//if (!i->second.get_bcube().closest_dist_xy_less_than(camera, draw_dist)) continue; // distance test (conservative)
			if (!dist_xy_less_than(camera, i->second.get_bcube().get_cube_center(), draw_dist)) continue; // distance test (aggressive)
			if (i->second.is_visible(xlate)) {bcs.push_back(&i->second);}
		}
	}
	void add_interior_lights(vector3d const &xlate, cube_t &lights_bcube, bool sec_camera_mode) {
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {
			cube_t const &bcube(i->second.get_bcube());
			if (!lights_bcube.intersects_xy(bcube))      continue; // not within light volume (too far from camera)
			if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
			i->second.add_interior_lights(xlate, lights_bcube, sec_camera_mode);
		}
	}
	void get_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state, bool cur_building_only) const {
		auto it(get_tile_by_pos_cs(pdu.pos));
		if (it != tiles.end()) {it->second.get_occluders(pdu, state, cur_building_only);}
	}
	bool check_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t const &state) const {
		auto it(get_tile_by_pos_cs(state.pos));
		return ((it == tiles.end()) ? 0 : it->second.check_pts_occluded(pts, npts, state));
	}
	unsigned get_tot_num_buildings() const {
		unsigned num(0);
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {num += i->second.get_num_buildings();}
		return num;
	}
	unsigned get_gpu_mem_usage() const {
		unsigned mem(0);
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {mem += i->second.get_gpu_mem_usage();}
		return mem;
	}
	void update_ai_state(float delta_dir) { // called once per frame
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {i->second.update_ai_state(delta_dir);}
	}
}; // end building_tiles_t


void occlusion_checker_noncity_t::set_camera(pos_dir_up const &pdu, bool cur_building_only) {
	if ((display_mode & 0x08) == 0) {state.building_ids.clear(); return;} // occlusion culling disabled
	pos_dir_up near_pdu(pdu);
	near_pdu.far_ = 0.5f*(X_SCENE_SIZE + Y_SCENE_SIZE); // set far clipping plane to half a tile (currently 4.0)
	bc.get_occluders(near_pdu, state, cur_building_only);
	//cout << "buildings: " << bc.get_num_buildings() << ", occluders: " << state.building_ids.size() << endl;
}
bool occlusion_checker_noncity_t::is_occluded(cube_t const &c) const {
	if (state.building_ids.empty()) return 0;
	float const z(c.z2()); // top edge
	point const corners[4] = {point(c.x1(), c.y1(), z), point(c.x2(), c.y1(), z), point(c.x2(), c.y2(), z), point(c.x1(), c.y2(), z)};
	return bc.check_pts_occluded(corners, 4, state);
}


building_creator_t building_creator(0), building_creator_city(1);
building_tiles_t building_tiles;

int create_buildings_tile(int x, int y, bool allow_flatten) { // return value: 0=already exists, 1=newly generaged, 2=re-generated
	if (!global_building_params.gen_inf_buildings()) return 0;
	return building_tiles.create_tile(x, y, allow_flatten);
}
bool remove_buildings_tile(int x, int y) {
	if (!global_building_params.gen_inf_buildings()) return 0;
	return building_tiles.remove_tile(x, y);
}

vector3d get_tt_xlate_val() {return ((world_mode == WMODE_INF_TERRAIN) ? vector3d(xoff*DX_VAL, yoff*DY_VAL, 0.0) : zero_vector);}

void gen_buildings() {
	global_building_params.finalize();
	update_sun_and_moon(); // need to update light_factor from sun to know if we need to generate window light geometry

	if (world_mode == WMODE_INF_TERRAIN && have_cities()) {
		building_creator_city.gen(global_building_params, 1, 0, 0, 1); // city buildings
		global_building_params.restore_prev_pos_range(); // hack to undo clip to city bounds to allow buildings to extend further out
		if (global_building_params.add_secondary_buildings) {building_creator.gen(global_building_params, 0, 1, 0, 1);} // non-city secondary buildings
	} else {building_creator .gen(global_building_params, 0, 0, 0, 1);} // mixed/non-city buildings
}
void draw_buildings(int shadow_only, int reflection_pass, vector3d const &xlate) {
	//if (!shadow_only && !reflection_pass && !building_tiles.empty()) {cout << "Building Tiles: " << building_tiles.size() << " Tiled Buildings: " << building_tiles.get_tot_num_buildings() << endl;} // debugging
	if (world_mode != WMODE_INF_TERRAIN) {building_tiles.clear();}
	building_creator_city.create_vbos(); // create VBOs for city buildings (after adding skyways, etc.), if needed
	vector<building_creator_t *> bcs;
	// don't draw city buildings for interior shadows
	bool const draw_city(world_mode == WMODE_INF_TERRAIN && (shadow_only != 2 || !interior_shadow_maps || global_building_params.add_city_interiors));
	bool const draw_sec ((shadow_only != 2 || interior_shadow_maps)); // don't draw secondary buildings for exterior dynamic shadows
	if (draw_city && building_creator_city.is_visible(xlate)) {bcs.push_back(&building_creator_city);}
	if (draw_sec  && building_creator     .is_visible(xlate)) {bcs.push_back(&building_creator     );}
	building_tiles.add_drawn(xlate, bcs);
	building_creator_t::multi_draw(shadow_only, reflection_pass, xlate, bcs);
}
void create_building_reflections() {
	building_creator_city.create_vbos(); // create VBOs for city buildings (after adding skyways, etc.), if needed
	building_creator_t::create_building_reflections(get_tiled_terrain_model_xlate());
}
void draw_player_building_transparent(int reflection_pass, vector3d const &xlate) {building_creator_t::draw_player_building_transparent(reflection_pass, xlate);}

void draw_building_lights(vector3d const &xlate) {
	building_creator_city.draw_building_lights(xlate);
	//building_creator.draw_building_lights(xlate); // only city buildings for now
}
bool proc_buildings_sphere_coll(point &pos, point const &p_int, float radius, vector3d *cnorm, bool check_interior, bool exclude_city) { // pos is in camera space
	if (check_interior) { // only called for the player
		player_in_closet       = 0; // reset for this call
		player_is_hiding       = 0;
		player_in_elevator     = 0;
		player_on_escalator    = 0;
		player_on_house_stairs = 0;
	}
	// we generally won't intersect more than one of these categories, so we can return true without checking all cases
	return ((!exclude_city && building_creator_city.check_sphere_coll(pos, p_int, radius, cnorm, check_interior)) ||
		                           building_creator.check_sphere_coll(pos, p_int, radius, cnorm, check_interior) ||
		                             building_tiles.check_sphere_coll(pos, p_int, radius, cnorm, check_interior));
}
bool check_buildings_no_grass(point const &pos) { // for tiled terrain mode; pos is in camera
	if (building_creator.check_point_coll_xy(pos)) return 1; // secondary buildings only
	if (building_tiles  .check_point_coll_xy(pos)) return 1;
	return 0;
}
bool check_buildings_cube_coll(cube_t const &c, bool xy_only, bool inc_basement, building_t const *exclude1, building_t const *exclude2) {
	return (building_creator_city.check_cube_coll(c, xy_only, inc_basement, exclude1, exclude2) ||
		building_creator.check_cube_coll(c, xy_only, inc_basement, exclude1, exclude2) ||
		building_tiles.check_cube_coll(c, xy_only, inc_basement, exclude1, exclude2));
}
void get_road_segs_in_region(cube_t const &region, vect_cube_t &out) { // for tiled terrain mode; pos is in local space
	building_creator.get_road_segs_in_region(region, out);
	building_tiles  .get_road_segs_in_region(region, out);
}
unsigned check_buildings_line_coll(point const &p1, point const &p2, float &t, unsigned &hit_bix, bool ret_any_pt) { // for line_intersect_city(); p1/p2 are in camera space
	vector3d const xlate(get_camera_coord_space_xlate());
	point const p1x(p1 - xlate), p2x(p2 - xlate); // convert from camera to building space
	unsigned const coll1(building_creator_city.check_line_coll(p1x, p2x, t, hit_bix, ret_any_pt, 0));
	if (coll1 && ret_any_pt) return coll1;
	unsigned const coll2(building_creator.check_line_coll(p1x, p2x, t, hit_bix, ret_any_pt, 1));
	if (coll2 && ret_any_pt) return coll2;
	unsigned const coll3(building_tiles.check_line_coll(p1x, p2x, t, ret_any_pt, 1)); // Note: does't take/set hit_bix
	return (coll3 ? coll3 : (coll2 ? coll2 : coll1));
}
bool check_city_building_line_coll_bs(point const &p1, point const &p2, point &p_int) { // Note: p1/p2 are in building space
	float t(1.0);
	unsigned hit_bix(0); // unused
	if (!building_creator_city.check_line_coll(p1, p2, t, hit_bix, 0, 0)) return 0; // ret_any_pt=0
	p_int = p1 + t*(p2 - p1);
	return 1;
}
bool check_city_building_line_coll_bs_any(point const &p1, point const &p2) { // Note: p1/p2 are in building space
	float t(1.0); // unused
	unsigned hit_bix(0); // unused
	return building_creator_city.check_line_coll(p1, p2, t, hit_bix, 1, 1); // ret_any_pt=1, no_coll_pt=1
}
void update_buildings_zmax_for_line(point const &p1, point const &p2, float radius, float house_extra_zval, float &cur_zmax) {
	building_creator_city.update_zmax_for_line(p1, p2, radius, house_extra_zval, cur_zmax);
	building_creator     .update_zmax_for_line(p1, p2, radius, house_extra_zval, cur_zmax);
	building_tiles       .update_zmax_for_line(p1, p2, radius, house_extra_zval, cur_zmax);
}
bool get_buildings_line_hit_color(point const &p1, point const &p2, colorRGBA &color) { // Note: p1 and p2 are in camera space
	vector3d const xlate(get_camera_coord_space_xlate());
	point const p1x(p1 - xlate), p2x(p2 - xlate); // convert from camera to building space
	if (world_mode == WMODE_INF_TERRAIN && building_creator_city.get_building_hit_color(p1x, p2x, color)) return 1;
	if (building_tiles.get_building_hit_color(p1x, p2x, color)) return 1;
	return building_creator.get_building_hit_color(p1x, p2x, color);
}
bool have_city_buildings() {return !building_creator_city.empty();}
bool have_secondary_buildings() {return (global_building_params.add_secondary_buildings && global_building_params.num_place > 0);}
bool have_buildings() {return (!building_creator.empty() || !building_creator_city.empty() || !building_tiles.empty());} // for postproc effects
bool no_grass_under_buildings() {return (world_mode == WMODE_INF_TERRAIN && !(building_creator.empty() && building_tiles.empty()) && global_building_params.flatten_mesh);}
unsigned get_buildings_gpu_mem_usage() {return (building_creator.get_gpu_mem_usage() + building_creator_city.get_gpu_mem_usage() + building_tiles.get_gpu_mem_usage());}
void add_city_building_signs(cube_t const &region_bcube, vector<sign_t     > &signs) {building_creator_city.add_building_signs(region_bcube, signs);}
void add_city_building_flags(cube_t const &region_bcube, vector<city_flag_t> &flags) {building_creator_city.add_building_flags(region_bcube, flags);}

vector3d get_buildings_max_extent() { // used for TT shadow bounds + map mode
	return building_creator.get_max_extent().max(building_creator_city.get_max_extent()).max(building_tiles.get_max_extent());
}
cube_t get_grid_bcube_for_building(building_t const &b) {
	cube_t ret(building_creator_city.get_grid_bcube_for_building(b));
	if (!ret.is_all_zeros()) return ret; // city building
	ret = building_creator.get_grid_bcube_for_building(b);
	if (!ret.is_all_zeros()) return ret; // secondary building
	ret = building_tiles.get_grid_bcube_for_building(b);
	return ret;
}
void clear_building_vbos() {
	building_creator     .clear_vbos();
	building_creator_city.clear_vbos();
	building_tiles       .clear_vbos();
}

// city interface
void set_buildings_pos_range(cube_t const &pos_range) {global_building_params.set_pos_range(pos_range);}
// Note: no xlate applied for any of these four queries below
void get_building_bcubes(cube_t const &xy_range, vect_cube_with_ix_t &bcubes  ) {building_creator_city.get_overlapping_bcubes(xy_range, bcubes);}
void get_building_bcubes(cube_t const &xy_range, vect_cube_t         &bcubes  ) {building_creator_city.get_overlapping_bcubes(xy_range, bcubes);}
void get_building_ext_basement_bcubes(cube_t const &city_bcube, vect_cube_t &bcubes) {building_creator_city.get_building_ext_basement_bcubes(city_bcube, bcubes);}
void get_walkways_for_city(cube_t const &city_bcube, vect_bldg_walkway_t &walkways ) {building_creator_city.get_walkways_for_city(city_bcube, walkways);}
void get_building_power_points(cube_t const &xy_range, vector<point> &ppts    ) {building_creator_city.get_power_points(xy_range, ppts);}
void add_house_driveways_for_plot(cube_t const &plot, vect_cube_t &driveways  ) {building_creator_city.add_house_driveways_for_plot(plot, driveways);}
void add_buildings_exterior_lights(vector3d const &xlate, cube_t &lights_bcube) {building_creator_city.add_exterior_lights(xlate, lights_bcube);}
float get_max_house_size() {return global_building_params.get_max_house_size();}

bool connect_buildings_to_skyway(cube_t &m_bcube, bool m_dim, cube_t const &city_bcube, vector<skyway_conn_t> &ww_conns) {
	return building_creator_city.connect_buildings_to_skyway(m_bcube, m_dim, city_bcube, ww_conns);
}
void get_city_building_walkways(cube_t const &city_bcube, vector<building_walkway_t *> &bwws) {
	building_creator_city.get_city_building_walkways(city_bcube, bwws);
}
void add_building_interior_lights(point const &xlate, cube_t &lights_bcube, bool sec_camera_mode) {
	//highres_timer_t timer("Add building interior lights"); // 0.97/0.37
	building_creator     .add_interior_lights(xlate, lights_bcube, sec_camera_mode);
	building_creator_city.add_interior_lights(xlate, lights_bcube, sec_camera_mode);
	building_tiles       .add_interior_lights(xlate, lights_bcube, sec_camera_mode);
}
// cars + peds
void get_city_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state) {building_creator_city.get_occluders(pdu, state);}
bool check_city_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t const &state) {return building_creator_city.check_pts_occluded(pts, npts, state);}
bool city_single_cube_visible_check(point const &pos, cube_t const &c) {return building_creator_city.single_cube_visible_check(pos, c);}
cube_t get_building_lights_bcube() {return building_lights_manager.get_lights_bcube();}
// used for pedestrians in cities
cube_t get_building_bcube(unsigned building_id) {return building_creator_city.get_building_bcube(building_id);}

bool get_building_door_pos_closest_to(unsigned building_id, point const &target_pos, point &door_pos, bool inc_garage_door) { // for city buildings only
	return building_creator_city.get_building_door_pos_closest_to(building_id, target_pos, door_pos, inc_garage_door);
}
cube_t register_deck_and_get_part_bounds(unsigned building_id, cube_t const &deck) {
	return building_creator_city.register_deck_and_get_part_bounds(building_id, deck);
}
bool check_sphere_coll_building(point const &pos, float radius, bool xy_only, unsigned building_id) {
	return building_creator_city.check_sphere_coll_building(pos, radius, xy_only, building_id);
}
bool check_building_point_or_cylin_contained(point const &pos, float radius, bool inc_details, unsigned building_id) {
	return building_creator_city.check_building_point_or_cylin_contained(pos, radius, inc_details, building_id);
}
int check_buildings_ped_coll(point const &pos, float bcube_radius, float detail_radius, unsigned plot_id, unsigned &building_id, cube_t *coll_cube) {
	return building_creator_city.check_ped_coll(pos, bcube_radius, detail_radius, plot_id, building_id, coll_cube);
}
bool check_line_coll_building(point const &p1, point const &p2, unsigned building_id) {return building_creator_city.check_line_coll_building(p1, p2, building_id);}
int get_building_bcube_contains_pos(point const &pos) {return building_creator_city.get_building_bcube_contains_pos(pos);}
bool select_building_in_plot(unsigned plot_id, unsigned rand_val, unsigned &building_id) {return building_creator_city.select_building_in_plot(plot_id, rand_val, building_id);}

// used for people in buildings
cube_t get_sec_building_bcube(unsigned building_id) {return building_creator.get_building_bcube(building_id);} // unused
bool enable_building_people_ai() {return global_building_params.enable_people_ai;}

void update_building_ai_state(float delta_dir) { // Note: each creator will manage its own range of people
	if (!global_building_params.enable_people_ai || !draw_building_interiors || !animate2) return;
	building_creator     .update_ai_state(delta_dir);
	building_creator_city.update_ai_state(delta_dir);
	building_tiles       .update_ai_state(delta_dir);
}

void get_all_city_helipads(vect_cube_t &helipads) {building_creator_city.get_all_helipads(helipads);} // city only for now

bool is_pos_in_player_building(point const &pos) { // pos is in global space
	if (!camera_in_building || player_building == nullptr) return 0;
	//static vector<point> points; // reused across calls
	//return player_building->check_point_or_cylin_contained(pos, 0.0, points);
	return player_building->check_point_xy_in_part(pos); // don't need to draw if above the building either, since there are no skylights
}
cube_t get_cur_basement() {return ((player_building != nullptr && player_building->has_basement()) ? player_building->get_basement() : cube_t());}

