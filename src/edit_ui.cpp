// 3D World - User Interface for Mesh Editing
// by Frank Gennari
// 10/26/13

#include "3DWorld.h"
#include "function_registry.h"
#include "heightmap.h" // for hmap_brush_t
#include "voxels.h" // for voxel_brush_params_t
#include "sphere_materials.h" // for sphere_mat_t
#include "shaders.h"

using namespace std;

float const MENU_TEXT_SIZE = 1.0;

extern bool camera_in_building;

bool is_shift_key_pressed();

class keyboard_menu_t {

protected:
	string title;
	unsigned num_controls, cur_control;
	mutable text_drawer text_draw;

	void draw_one_control_text(unsigned control_ix, string const &name, string const &cur_value, float slider_pos,
		bool reversed=0, float scale=1.0, colorRGBA const &color=ALPHA0, int tid=-1) const
	{
		if (reversed) {control_ix = num_controls - control_ix - 1;}
		//assert(slider_pos >= 0.0 && slider_pos <= 1.0);
		slider_pos = CLIP_TO_01(slider_pos);
		assert(control_ix < num_controls);
		bool const selected(control_ix == cur_control);
		unsigned const ndiv = 20;
		unsigned const pos(round_fp((ndiv-1)*slider_pos));
		ostringstream oss;
		for (unsigned n = 0; n < pos; ++n) {oss << '-';}
		oss << '+';
		for (unsigned n = pos+1; n < ndiv; ++n) {oss << '-';}
		oss << "  " << name << ": " << cur_value;
		float const mag(0.01f/scale), y_step(0.0014f*scale), yval(1.0f*mag - y_step*(num_controls - control_ix)), zval(-2.0f*mag);
		text_draw.add_text((selected ? ORANGE : YELLOW), -1.0*mag, yval, zval, oss.str().c_str(), MENU_TEXT_SIZE);
		if (color.A > 0.0) {text_draw.add_text(color, -1.06*mag, yval, zval, "O", MENU_TEXT_SIZE);} // draw color indicator

		if (tid >= 0) {
			float const sz(0.75*y_step);
			text_draw.flush();
			text_draw.set_color(WHITE);
			select_texture(tid);
			draw_one_tquad(-1.08*mag-sz, yval-0.5*sz, -1.08*mag+sz, yval+1.5*sz, zval);
			text_draw.bind_font_texture(); // restore font texture
		}
	}
	virtual void draw_one_control(unsigned control_ix) const = 0;

public:
	// Note: the max number of controls we can fit on the screen is around 13
	keyboard_menu_t(unsigned num_controls_, string const &title_) : title(title_), num_controls(num_controls_), cur_control(0) {assert(num_controls > 0);}
	virtual ~keyboard_menu_t() {}
	virtual bool is_enabled() const = 0;
	void next_control() {cur_control = (cur_control == num_controls-1) ? 0  : (cur_control+1);}
	void prev_control() {cur_control = (cur_control == 0) ? (num_controls-1): (cur_control-1);}
	virtual void change_value(int delta) = 0;

	void draw_controls() const {
		text_draw.begin_draw();
		if (!title.empty()) {text_draw.add_text(YELLOW, -0.01, 0.01, -0.02, title.c_str(), MENU_TEXT_SIZE);}
		for (unsigned i = 0; i < num_controls; ++i) {draw_one_control(i);}
		text_draw.end_draw();
	}
};


// ************ Tiled Terrain Heightmap Brush/Editing ************

enum {HMAP_DELAY=0, HMAP_CTR_SHAPE, HMAP_CTR_RADIUS, HMAP_CTR_DELTA, HMAP_NUM_CTR};
string const hmap_ctr_names[HMAP_NUM_CTR] = {"Placement Delay", "Brush Shape", "Brush Radius", "Brush Delta"};
string const brush_shapes[NUM_BSHAPES] = {"Constant Square", "Constant Circle", "Linear Circle", "Quadratic Circle", "Cosine Circle", "Sine Circle", "Flat Square", "Flat Circle"};

extern int show_scores, world_mode, game_mode;
extern unsigned inf_terrain_fire_mode;
extern hmap_brush_param_t cur_brush_param;

unsigned get_tile_size();


class hmap_kbd_menu_t : public keyboard_menu_t {

	hmap_brush_param_t &brush_param;
	unsigned max_radius_exp;

	virtual void draw_one_control(unsigned control_ix) const {
		assert(control_ix < HMAP_NUM_CTR);
		ostringstream value;
		float spos(0.0);

		switch (control_ix) {
		case HMAP_DELAY:
			spos = (brush_param.delay/10.0);
			value << brush_param.delay;
			break;
		case HMAP_CTR_SHAPE:
			spos = (brush_param.shape/float(BSHAPE_SINE));
			assert(brush_param.shape < NUM_BSHAPES);
			value << brush_shapes[brush_param.shape];
			break;
		case HMAP_CTR_RADIUS:
			spos = ((brush_param.radius_exp + 1)/float(max_radius_exp + 1));
			value << (brush_param.get_radius() ? brush_param.get_radius() : 0.5);
			break;
		case HMAP_CTR_DELTA:
			spos = (brush_param.delta_exp/9.0);
			value << brush_param.get_delta_mag();
			break;
		default: assert(0);
		}
		draw_one_control_text(control_ix, hmap_ctr_names[control_ix], value.str(), spos);
	}

public:
	hmap_kbd_menu_t(hmap_brush_param_t &bp) : keyboard_menu_t(HMAP_NUM_CTR, "Heightmap Edit"), brush_param(bp), max_radius_exp(0) {
		for (unsigned sz = 1; sz < get_tile_size(); sz <<= 1) {++max_radius_exp;}
	}
	virtual bool is_enabled() const {return (show_scores && !game_mode && inf_terrain_fire_mode && world_mode == WMODE_INF_TERRAIN);}

	virtual void change_value(int delta) {
		switch (cur_control) {
		case HMAP_DELAY:
			brush_param.delay = max(0, min(10, ((int)brush_param.delay + delta)));
			break;
		case HMAP_CTR_SHAPE:
			brush_param.shape = max(0, min((int)BSHAPE_SINE, (brush_param.shape + delta)));
			break;
		case HMAP_CTR_RADIUS:
			brush_param.radius_exp = max(-1, min((int)max_radius_exp, ((int)brush_param.radius_exp + delta)));
			break;
		case HMAP_CTR_DELTA:
			brush_param.delta_exp = max(0, min(9, ((int)brush_param.delta_exp + delta)));
			break;
		default: assert(0);
		}
	}
};


// ************ Voxel Editing ************

enum {VOXEL_DELAY=0, VOXEL_SHAPE, VOXEL_RADIUS, VOXEL_WEIGHT, NUM_VOXEL_CONT};
string const voxel_ctr_names[NUM_VOXEL_CONT] = {"Brush Delay", "Brush Shape", "Brush Radius", "Brush Weight"};
string const vb_shape_names [NUM_VB_SHAPES ] = {"Constant Cube", "Constant Sphere", "Linear Sphere", "Quadratic Sphere"};

unsigned const MAX_VB_RADIUS = 20; // in voxel dx size units
int const MAX_VB_WEIGHT_EXP  = 4;

extern int voxel_editing;
extern voxel_brush_params_t voxel_brush_params;

float get_voxel_brush_step();

class voxel_edit_kbd_menu_t : public keyboard_menu_t {

	voxel_brush_params_t &brush;

	virtual void draw_one_control(unsigned control_ix) const {
		assert(control_ix < NUM_VOXEL_CONT);
		ostringstream value;
		float spos(0.0);

		switch (control_ix) {
		case VOXEL_DELAY:
			spos = (brush.delay/10.0);
			value << brush.delay;
			break;
		case VOXEL_SHAPE:
			spos = (brush.shape/float(NUM_VB_SHAPES-1));
			assert(brush.shape < NUM_VB_SHAPES);
			value << vb_shape_names[brush.shape];
			break;
		case VOXEL_RADIUS:
			spos = float(brush.radius-1)/float(MAX_VB_RADIUS-1);
			value.precision(1);
			value << fixed; // fixed precision in units of 0.1
			value << brush.get_world_space_radius() << " (" << brush.radius << " units)";
			break;
		case VOXEL_WEIGHT:
			spos = 0.5f*(brush.weight_exp + MAX_VB_WEIGHT_EXP)/MAX_VB_WEIGHT_EXP;
			value << pow(2.0f, brush.weight_exp) * ((voxel_editing == 2) ? -1.0 : 1.0);
			break;
		default: assert(0);
		}
		draw_one_control_text(control_ix, voxel_ctr_names[control_ix], value.str(), spos);
	}

public:
	voxel_edit_kbd_menu_t(voxel_brush_params_t &brush_) : keyboard_menu_t(NUM_VOXEL_CONT, "Voxel Edit"), brush(brush_) {}
	virtual bool is_enabled() const {return (show_scores && !game_mode && voxel_editing && world_mode == WMODE_GROUND);}

	virtual void change_value(int delta) {
		switch (cur_control) {
		case VOXEL_DELAY:
			brush.delay = max(0, min(10, ((int)brush.delay + delta)));
			break;
		case VOXEL_SHAPE:
			brush.shape = max(0, min(NUM_VB_SHAPES-1, (brush.shape + delta)));
			break;
		case VOXEL_RADIUS:
			brush.radius = max(1, min((int)MAX_VB_RADIUS, ((int)brush.radius + delta)));
			break;
		case VOXEL_WEIGHT:
			brush.weight_exp = max(-MAX_VB_WEIGHT_EXP, min(MAX_VB_WEIGHT_EXP, (brush.weight_exp + delta)));
			break;
		default: assert(0);
		}
	}
};


// ************ Sun/Tree Leaf/Grass Colors/Properties ************

enum {GRASS_WIDTH=0, GRASS_LENGTH, TREE_COLOR_VAR, LEAF_COLOR_VAR, LEAF_RED_COMP, LEAF_GREEN_COMP, LEAF_BLUE_COMP, SUN_RED, SUN_GREEN, SUN_BLUE, NUM_LEAF_CONT};
string const leaf_ctr_names[NUM_LEAF_CONT] = {"Grass Width", "Grass Length", "Tree Color Variance", "Leaf Color Variance",
	"Leaf Red Component", "Leaf Green Component", "Leaf Blue Component", "Sun Red Color", "Sun Green Color", "Sun Blue Color"};

extern int leaf_color_changed;
extern float leaf_color_coherence, tree_color_coherence; // [0.0, 1.0] in steps of 0.1
extern float grass_length, grass_width;
extern colorRGBA leaf_base_color; // can set R and G in [-1.0, 1.0] in steps of 0.1
extern colorRGBA sunlight_color;

void update_grass_length_width(float new_grass_length, float new_grass_width);
void update_tiled_grass_colors();
void register_leaf_color_change();

class leaf_color_kbd_menu_t : public keyboard_menu_t {

	virtual void draw_one_control(unsigned control_ix) const {
		assert(control_ix < NUM_LEAF_CONT);
		ostringstream value;
		value.precision(1);
		value << fixed; // fixed precision in units of 0.1
		float spos(0.0);
		colorRGBA color(ALPHA0);

		switch (control_ix) {
		case TREE_COLOR_VAR:
			spos = tree_color_coherence;
			value << spos;
			break;
		case LEAF_COLOR_VAR:
			spos = (1.0 - leaf_color_coherence);
			value << spos;
			break;
		case GRASS_LENGTH:
			value.precision(3);
			spos = 10.0*grass_length - 0.05;
			value << grass_length;
			break;
		case GRASS_WIDTH:
			value.precision(4);
			spos = 100.0*grass_width - 0.05;
			value << grass_width;
			break;
		case LEAF_RED_COMP:
		case LEAF_GREEN_COMP:
		case LEAF_BLUE_COMP:
			spos  = leaf_base_color[control_ix-LEAF_RED_COMP];
			value << spos;
			spos  = 0.5*(spos + 1.0); // center around 0.0
			color = (leaf_base_color + WHITE)*0.5;
			break;
		case SUN_RED:
		case SUN_GREEN:
		case SUN_BLUE:
			spos  = sunlight_color[control_ix-SUN_RED];
			value << spos;
			color = sunlight_color;
			break;
		default: assert(0);
		}
		draw_one_control_text(control_ix, leaf_ctr_names[control_ix], value.str(), spos, 0, 1.0, color);
	}

public:
	leaf_color_kbd_menu_t() : keyboard_menu_t(NUM_LEAF_CONT, "Tree Leaves and Grass") {}
	virtual bool is_enabled() const {return (show_scores && !game_mode && (world_mode == WMODE_GROUND || world_mode == WMODE_INF_TERRAIN));}

	virtual void change_value(int delta) {
		switch (cur_control) {
		case TREE_COLOR_VAR:
			tree_color_coherence = max(0.0f, (tree_color_coherence + 0.1f*delta)); // no maximum
			break;
		case LEAF_COLOR_VAR:
			leaf_color_coherence = CLIP_TO_01(leaf_color_coherence - 0.1f*delta); // delta is backwards
			update_tiled_grass_colors();
			break;
		case GRASS_LENGTH:
			update_grass_length_width(max(0.005f, (grass_length + 0.005f*delta)), grass_width); // steps of 0.005
			break;
		case GRASS_WIDTH:
			update_grass_length_width(grass_length, max(0.0005f, (grass_width + 0.0005f*delta))); // steps of 0.0005
			break;
		case LEAF_RED_COMP:
		case LEAF_GREEN_COMP:
		case LEAF_BLUE_COMP:
			leaf_base_color[cur_control-LEAF_RED_COMP] = CLIP_TO_pm1(leaf_base_color[cur_control-LEAF_RED_COMP] + 0.1f*delta); // steps of 0.1
			update_tiled_grass_colors();
			break;
		case SUN_RED:
		case SUN_GREEN:
		case SUN_BLUE:
			sunlight_color[cur_control-SUN_RED] = CLIP_TO_01(sunlight_color[cur_control-SUN_RED] + 0.05f*delta); // steps of 0.05
			break;
		default: assert(0);
		}
		register_leaf_color_change();
	}
};


// ************ Water Colors/Properties ************

colorRGB uw_atten_max(WHITE), uw_atten_scale(BLACK);
water_params_t water_params;
extern bool water_is_lava;

void water_params_t::set_def_water() {alpha=0.5; mud=0.0; algae=0.0; bright=1.0; reflect=1.0; green=0.16; wave_amp=1.0;}
void water_params_t::set_def_lava () {alpha=1.5; mud=0.0; algae=0.0; bright=1.0; reflect=0.5; green=0.02; wave_amp=1.5;}

colorRGB get_underwater_atten_color(float mud_amt, float algae_amt) {
	colorRGBA const c(blend_color(colorRGB(1.5, 0.5, 0.9), colorRGB(1.5, 0.9, 0.5), algae_amt, 0));
	return blend_color(colorRGB(0.9, 1.0, 1.5), c, mud_amt, 0);
}
void calc_uw_atten_colors() {
	if (water_is_lava) {uw_atten_scale = colorRGB(1.0, 1.5, 2.2);} // red-orange
	else               {uw_atten_scale = get_underwater_atten_color(water_params.mud, water_params.algae);} // blend in mud color
	//UNROLL_3X(uw_atten_max[i_] = CLIP_TO_01(1.0f - 0.03f/uw_atten_scale[i_]);)
	uw_atten_scale *= 0.05 + 0.95*water_params.alpha;
}

enum {WATERP_ALPHA=0, WATERP_MUD, WATERP_BRIGHT, WATERP_REFLECT, WATERP_GREEN, WATERP_WAVE_AMP, NUM_WATER_CONT};
string const water_ctr_names[NUM_WATER_CONT] = {"Alpha Scale", "Mud Content", "Brightness", "Reflectivity", "Green Hue", "Wave Amplitude"};

class water_color_kbd_menu_t : public keyboard_menu_t {

	virtual void draw_one_control(unsigned control_ix) const {
		assert(control_ix < NUM_WATER_CONT);
		ostringstream value;
		value.precision(2);
		value << fixed; // fixed precision in units of 0.01
		float spos(0.0);

		switch (control_ix) {
		case WATERP_ALPHA:
			value << water_params.alpha;
			spos = water_params.alpha/1.5; // 0.0 to 1.5
			break;
		case WATERP_MUD:
			value << water_params.mud;
			spos = water_params.mud; // 0.0 to 1.0
			break;
		case WATERP_BRIGHT:
			value << water_params.bright;
			spos = 0.5*water_params.bright; // 0.0 to 2.0
			break;
		case WATERP_REFLECT:
			value << water_params.reflect;
			spos = water_params.reflect; // 0.0 to 1.0
			break;
		case WATERP_GREEN:
			value << water_params.green;
			spos = 2.5*water_params.green; // 0.0 to 0.5
			break;
		case WATERP_WAVE_AMP:
			value << water_params.wave_amp;
			spos = 0.25*water_params.wave_amp; // 0.0 to 4.0
			break;
		default:
			assert(0);
		}
		draw_one_control_text(control_ix, water_ctr_names[control_ix], value.str(), spos);
	}

public:
	water_color_kbd_menu_t() : keyboard_menu_t(NUM_WATER_CONT, "Water Colors") {}
	virtual bool is_enabled() const {return (show_scores && !game_mode && world_mode == WMODE_INF_TERRAIN);}

	virtual void change_value(int delta) {
		switch (cur_control) {
		case WATERP_ALPHA:
			water_params.alpha = max(0.0f, min(1.5f, (water_params.alpha + 0.05f*delta))); // 0.0 to 1.5 in steps of 0.05
			break;
		case WATERP_MUD:
			water_params.mud = CLIP_TO_01(water_params.mud + 0.05f*delta); // 0.0 to 1.0 in steps of 0.05
			break;
		case WATERP_BRIGHT:
			water_params.bright = max(0.0f, min(2.0f, (water_params.bright + 0.1f*delta))); // 0.0 to 2.0 in steps of 0.1
			break;
		case WATERP_REFLECT:
			water_params.reflect = CLIP_TO_01(water_params.reflect + 0.05f*delta); // 0.0 to 1.0 in steps of 0.05
			break;
		case WATERP_GREEN:
			water_params.green = max(0.0f, min(0.5f, (water_params.green + 0.02f*delta))); // 0.0 to 0.4 in steps of 0.02
			break;
		case WATERP_WAVE_AMP:
			water_params.wave_amp = max(0.0f, min(4.0f, (water_params.wave_amp + 0.2f*delta))); // 0.0 to 4.0 in steps of 0.2
			break;
		default:
			assert(0);
		}
		calc_uw_atten_colors();
	}
};


// ************ Physics and Weather ************

extern int temp_change, def_cube_map_reflect_mipmap_level;
extern float base_gravity, temperature, cloud_cover, sun_rot, moon_rot, ball_velocity, TIMESTEP, player_speed, def_atmosphere, def_vegetation, ambient_scale;
extern vector3d wind;

enum {PW_SPEED=0, PW_GRAVITY, PW_TEMP, PW_WATER, PW_VEG, PW_ATMOS, PW_CLOUD, PW_PRECIP, PW_WIND_X, PW_WIND_Y, PW_SUN_POS, PW_MOON_POS, PM_AMB_SCALE, PW_TIMESTEP, PW_WVEL, CM_MIP_BIAS, NUM_PW_CONT};
string const phys_weather_names[NUM_PW_CONT] = {"Player Speed", "Gravity", "Temperature", "Water Level", "Vegetation", "Atmosphere", "Cloudiness", "Precipitation",
												"Wind X", "Wind Y", "Sun Angle", "Moon Angle", "Ambient Light Scale", "Physics Timestep", "Weapon Velocity", "Cube Map Mipmap Bias"};

class phys_weather_kbd_menu_t : public keyboard_menu_t {

	virtual void draw_one_control(unsigned control_ix) const {
		assert(control_ix < NUM_PW_CONT);
		ostringstream value;
		value.precision(2);
		value << fixed; // fixed precision in units of 0.01
		float spos(0.0);

		switch (control_ix) {
		case PW_SPEED:
			value << player_speed;
			spos = player_speed/4.0; // 0.0 to 4.0
			break;
		case PW_GRAVITY:
			value << base_gravity;
			spos = base_gravity/2.0; // 0.0 to 2.0
			break;
		case PW_TEMP: // Note: doesn't update in tiled terrain mode
			value << temperature;
			spos = (temperature + 40.0)/160.0; // -40.0 to 120.0
			break;
		case PW_WATER:
			value << get_rel_wpz();
			spos = get_rel_wpz(); // 0.0 to 1.0
			break;
		case PW_VEG:
			value << def_vegetation;
			spos = def_vegetation; // 0.0 to 1.0
			break;
		case PW_ATMOS:
			value << def_atmosphere;
			spos = def_atmosphere; // 0.0 to 1.0
			break;
		case PW_CLOUD:
			value << cloud_cover;
			spos = (cloud_cover + 0.5)/1.5; // -0.5 to 1.0
			break;
		case PW_PRECIP:
			value << get_precip_rate();
			spos = get_precip_rate()/1000.0; // 0.0 to 1.0
			break;
		case PW_WIND_X:
			value << wind.x;
			spos = 0.25*(wind.x + 2.0); // -2.0 to 2.0
			break;
		case PW_WIND_Y:
			value << wind.y;
			spos = 0.25*(wind.y + 2.0); // -2.0 to 2.0
			break;
		case PW_SUN_POS:
			value << sun_rot;
			spos = sun_rot/TWO_PI; // 0 to 2*PI
			break;
		case PW_MOON_POS:
			value << moon_rot;
			spos = moon_rot/TWO_PI; // 0 to 2*PI
			break;
		case PM_AMB_SCALE:
			value << ambient_scale;
			spos = 0.5*ambient_scale; // 0.0 to 2.0
			break;
		case PW_TIMESTEP:
			value << 100.0*TIMESTEP;
			spos = 100.0*TIMESTEP; // 0.0 to 0.01
			break;
		case PW_WVEL:
			value << ball_velocity;
			spos = ball_velocity/40.0; // 0.0 to 40.0
			break;
		case CM_MIP_BIAS:
			value << def_cube_map_reflect_mipmap_level;
			spos = (def_cube_map_reflect_mipmap_level + 10)/20.0f; // -10 to 10
			break;
		default:
			assert(0);
		}
		draw_one_control_text(control_ix, phys_weather_names[control_ix], value.str(), spos, 0, 0.9);
	}

public:
	phys_weather_kbd_menu_t() : keyboard_menu_t(NUM_PW_CONT, "Physics and Weather") {}
	virtual bool is_enabled() const {return (show_scores && !game_mode && world_mode != WMODE_UNIVERSE);}

	virtual void change_value(int delta) {
		int regen_mesh(0);

		switch (cur_control) {
		case PW_SPEED:
			player_speed = max(0.0f, (player_speed + 0.2f*delta)); // >= 0.0 in steps of 0.2
			break;
		case PW_GRAVITY:
			base_gravity = max(0.0f, (base_gravity + 0.05f*delta)); // >= 0.0 in steps of 0.1
			break;
		case PW_TEMP:
			temperature += TEMP_INCREMENT*delta; // in steps of 10.0
			temp_change  = 1;
			//regen_mesh   = 1; // regen texture
			break;
		case PW_WATER: {
				float const water_level(max(-1.0f, min(1.0f, (get_rel_wpz() + 0.05f*delta)))); // -1.0 to 1.0 in steps of 0.05
				change_water_level(water_level); // 0.0 to 1.0
			}
			regen_mesh = 2; // regen texture
			break;
		case PW_VEG:
			def_vegetation = CLIP_TO_01(def_vegetation + 0.1f*delta); // 0.0 to 1.0 in steps of 0.1
			regen_mesh     = 1; // regen texture
			break;
		case PW_ATMOS:
			def_atmosphere = CLIP_TO_01(def_atmosphere + 0.1f*delta); // 0.0 to 1.0 in steps of 0.1
			break;
		case PW_CLOUD:
			cloud_cover = max(-0.5f, min(1.0f, (cloud_cover + 0.1f*delta))); // -0.5 to 1.0 in steps of 0.1
			break;
		case PW_PRECIP:
			update_precip_rate((delta > 0.0) ? 1.5 : 1.0/1.5);
			break;
		case PW_WIND_X:
			wind.x = (wind.x + WIND_ADJUST*delta); // -2.0 to 2.0 in steps of 0.2
			break;
		case PW_WIND_Y:
			wind.y = (wind.y + WIND_ADJUST*delta); // -2.0 to 2.0 in steps of 0.2
			break;
		case PW_SUN_POS:
			sun_rot += 0.5*LIGHT_ROT_AMT*delta; // -PI to PI in steps of 0.025
			if (sun_rot < 0.0) {sun_rot += TWO_PI;} else if (sun_rot > TWO_PI) {sun_rot -= TWO_PI;}
			update_sun_shadows();
			break;
		case PW_MOON_POS:
			moon_rot += 0.5*LIGHT_ROT_AMT*delta; // -PI to PI in steps of 0.025
			if (moon_rot < 0.0) {moon_rot += TWO_PI;} else if (moon_rot > TWO_PI) {moon_rot -= TWO_PI;}
			calc_visibility(MOON_SHADOW);
			break;
		case PM_AMB_SCALE:
			ambient_scale = max(0.0f, (ambient_scale + 0.1f*delta)); // 0.0 to 2.0+ in steps of 0.1
			break;
		case PW_TIMESTEP:
			change_timestep((delta > 0.0) ? 1.5 : 1.0/1.5);
			break;
		case PW_WVEL:
			ball_velocity = max(0.0f, (ball_velocity + 2.0f*delta)); // >= 0.0 in steps of 2.0
			break;
		case CM_MIP_BIAS:
			def_cube_map_reflect_mipmap_level = max(-10, min(10, def_cube_map_reflect_mipmap_level+delta)); // -10 to 10 in steps of 1
			break;
		default:
			assert(0);
		} // end switch
		if (regen_mesh) {
			init_terrain_mesh(); // Note: texture not actually regenerated
			clear_tiled_terrain(regen_mesh == 1);
		}
	}
};


// ************ Sphere Materials ************

extern unsigned spheres_mode;
extern float sphere_mat_fire_delay;
extern int coll_id[];
extern vector<texture_t> textures;

enum {SM_MAT_NAME=0, SM_TEXTURE, SM_FDELAY, SM_EMISS, SM_REFLECT, SM_DESTROY, SM_RSCALE, SM_HARDNESS, SM_DENSITY, SM_METAL, SM_ALPHA, SM_SPEC_MAG, SM_SHINE,
	SM_REFRACT_IX, SM_LIGHT_ATTEN, SM_LIGHT_RADIUS, SM_LIGHT_SHADOW, SM_DIFF_R, SM_DIFF_G, SM_DIFF_B, SM_SPEC_R, SM_SPEC_G, SM_SPEC_B, NUM_SM_CONT};
string const sphere_mode_names[NUM_SM_CONT] =
{"Material Name", "Texture", "Fire Delay", "Emissive ", "Reflective", "Destroyable", "Radius Scale", "Hardness ", "Density   ", "Metalness ", "Alpha    ", "Specular Mag", "Shininess  ",
"Refract Ix  ", "Light Atten ", "Light Radius", "Light Shadow", "Diffuse Red ", "Diffuse Green", "Diffuse Blue ", "Specular Red ", "Specular Green", "Specular Blue "};
string const destroy_tags[3] = {"Indestructible", "Shatters", "Explodes"};

class sphere_mat_kbd_menu_t : public keyboard_menu_t {

	virtual void draw_one_control(unsigned control_ix) const {
		assert(control_ix < NUM_SM_CONT);
		ostringstream value;
		value.precision(2);
		value << fixed; // fixed precision in units of 0.01
		float spos(0.0);
		sphere_mat_t const &mat(get_cur_sphere_mat());
		colorRGBA color(ALPHA0);
		int tid(-1);

		switch (control_ix) {
		case SM_MAT_NAME:     value << mat.get_name(); break; // spos stays at 0
		case SM_TEXTURE:      value << mat.tid << " " << ((mat.tid < 0) ? "None" : textures[mat.tid].name); tid = mat.tid; break; // spos stays at 0
		case SM_FDELAY:       value << sphere_mat_fire_delay; spos = sphere_mat_fire_delay; break; // 0.0 to 1.0 seconds
		case SM_EMISS:        value << mat.emissive;     spos = mat.emissive;         break; // 0/1
		case SM_REFLECT:      value << mat.reflective;   spos = mat.reflective;       break; // 0/1
		case SM_DESTROY:      value << destroy_tags[mat.destroyable]; spos = float(mat.destroyable/2); break; // 0/1
		case SM_RSCALE:       value << mat.radius_scale; spos = (mat.radius_scale-0.2)/1.8; break; // 0.2 to 2.0
		case SM_HARDNESS:     value << mat.hardness;     spos = mat.hardness;         break; // 0.05 to 1.0
		case SM_DENSITY:      value << mat.density;      spos = mat.density/4.0;      break; // 0.1 to 4.0
		case SM_METAL:        value << mat.metal;        spos = mat.metal;            break; // 0.0 to 1.0
		case SM_ALPHA:        value << mat.alpha;        spos = mat.alpha;            break; // 0.05 to 1.0
		case SM_SPEC_MAG:     value << mat.spec_mag;     spos = mat.spec_mag;         break; // 0.0 to 1.0
		case SM_SHINE:        value << int(mat.shine);   spos = log2(mat.shine)/10.0; break; // 1.0 to 1024.0 (log2 0 to 10)
		case SM_REFRACT_IX:   value << mat.refract_ix;   spos = (mat.refract_ix-1.0)/3.0; break; // 1.0 to 4.0
		case SM_LIGHT_ATTEN:  value << mat.light_atten;  spos = mat.light_atten/25.0; break; // 0.0 to 25.0
		case SM_LIGHT_RADIUS: value << mat.light_radius; spos = mat.light_radius/2.0; break; // 0.0 to 2.5
		case SM_LIGHT_SHADOW: value << mat.shadows;      spos = mat.shadows;          break; // 0/1
		case SM_DIFF_R:       value << mat.diff_c.R;     spos = mat.diff_c.R;         color = mat.diff_c; break; // 0.0 to 1.0
		case SM_DIFF_G:       value << mat.diff_c.G;     spos = mat.diff_c.G;         color = mat.diff_c; break; // 0.0 to 1.0
		case SM_DIFF_B:       value << mat.diff_c.B;     spos = mat.diff_c.B;         color = mat.diff_c; break; // 0.0 to 1.0
		case SM_SPEC_R:       value << mat.spec_c.R;     spos = mat.spec_c.R;         color = mat.spec_c; break; // 0.0 to 1.0
		case SM_SPEC_G:       value << mat.spec_c.G;     spos = mat.spec_c.G;         color = mat.spec_c; break; // 0.0 to 1.0
		case SM_SPEC_B:       value << mat.spec_c.B;     spos = mat.spec_c.B;         color = mat.spec_c; break; // 0.0 to 1.0
		default: assert(0);
		}
		draw_one_control_text(control_ix, sphere_mode_names[control_ix], value.str(), spos, 1, 0.75, color, tid); // draw in reverse order at 80% scale
	}
	static void change_texture(int &tid, int &nm_tid, int delta) {
		do {
			tid += delta;
			if (tid == -1) break; // first texture => no texture
			if (tid >= (int)textures.size()) {tid = -1; break;} // last texture => no texture
			if (tid < 0) {tid = textures.size()-1;} // no texture => last texture
		} while (textures[tid].normal_map || textures[tid].type > 0); // skip normal maps and generated textures
		nm_tid = ((tid < 0) ? -1 : textures[tid].bump_tid); // usually unset/-1
	}

public:
	sphere_mat_kbd_menu_t() : keyboard_menu_t(NUM_SM_CONT, "Sphere Materials") {cur_control = num_controls-1;} // start at material name
	virtual bool is_enabled() const {return (show_scores && spheres_mode != 0 && !game_mode);}

	virtual void change_value(int delta) {
		sphere_mat_t &mat(get_cur_sphere_mat());
		unsigned const ix(num_controls - cur_control - 1);

		switch (ix) { // reverse order
		case SM_MAT_NAME:     change_sphere_material(delta, 1);           break;
		case SM_TEXTURE:      change_texture(mat.tid, mat.nm_tid, delta); break; // reset normal map because it won't go with the texture when changed
		case SM_FDELAY:       sphere_mat_fire_delay = CLIP_TO_01(sphere_mat_fire_delay + 0.05f*delta);  break; // 0.0 to 1.0 in steps of 0.05
		case SM_EMISS:        mat.emissive     = ((delta < 0) ? 0 : 1);   break; // 0/1
		case SM_REFLECT:      mat.reflective   = ((delta < 0) ? 0 : 1);   break; // 0/1
		case SM_DESTROY:      mat.destroyable  = min(2, max(0, mat.destroyable+delta)); break; // 0-2
		case SM_RSCALE:       mat.radius_scale = max(0.2f, min(2.0f, (mat.radius_scale + 0.1f*delta))); break; // 0.2 to 2.0 in steps of 0.1
		case SM_HARDNESS:     mat.hardness     = max(0.05f, min(1.0f, (mat.hardness + 0.05f*delta)));   break; // 0.05 to 1.0 in steps of 0.05
		case SM_DENSITY:      mat.density      = max(0.1f,  min(4.0f, (mat.density  + 0.1f* delta)));   break; // 0.1 to 4.0 in steps of 0.1
		case SM_METAL:        mat.metal        = CLIP_TO_01(mat.metal    + 0.05f*delta);                break; // 0.0 to 1.0 in steps of 0.05
		case SM_ALPHA:        mat.alpha        = CLIP_TO_01(mat.alpha    + 0.05f*delta);                break; // 0.0 to 1.0 in steps of 0.05 (alpha of 0 is okay for refractive or scattering objects)
		case SM_SPEC_MAG:     mat.spec_mag     = CLIP_TO_01(mat.spec_mag + 0.05f*delta);                break; // 0.0 to 1.0 in steps of 0.05
		case SM_SHINE:        mat.shine        = pow(2.0, max(0, min(10, int(log2(mat.shine)+delta)))); break; // 1.0 to 1024.0 (log2 0 to 10 in steps of 1)
		case SM_REFRACT_IX:   mat.refract_ix   = max(1.0f, min(4.0f, (mat.refract_ix   + 0.1f*delta))); break; // 1.0 to 4.0 in steps of 0.1
		case SM_LIGHT_ATTEN:  mat.light_atten  = max(0.0f, min(25.0f,(mat.light_atten  + 1.0f*delta))); break; // 0.0 to 25.0 in steps of 1.0
		case SM_LIGHT_RADIUS: mat.light_radius = max(0.0f, min(2.5f, (mat.light_radius + 0.1f*delta))); if (mat.light_radius < 0.01) {mat.light_radius = 0.0;} break; // 0.0 to 2.5 in steps of 0.1
		case SM_LIGHT_SHADOW: mat.shadows      = ((delta < 0) ? 0 : 1); break; // 0/1
		case SM_DIFF_R:       mat.diff_c.R     = CLIP_TO_01(mat.diff_c.R + 0.05f*delta); break; // 0.0 to 1.0 in steps of 0.05
		case SM_DIFF_G:       mat.diff_c.G     = CLIP_TO_01(mat.diff_c.G + 0.05f*delta); break; // 0.0 to 1.0 in steps of 0.05
		case SM_DIFF_B:       mat.diff_c.B     = CLIP_TO_01(mat.diff_c.B + 0.05f*delta); break; // 0.0 to 1.0 in steps of 0.05
		case SM_SPEC_R:       mat.spec_c.R     = CLIP_TO_01(mat.spec_c.R + 0.05f*delta); break; // 0.0 to 1.0 in steps of 0.05
		case SM_SPEC_G:       mat.spec_c.G     = CLIP_TO_01(mat.spec_c.G + 0.05f*delta); break; // 0.0 to 1.0 in steps of 0.05
		case SM_SPEC_B:       mat.spec_c.B     = CLIP_TO_01(mat.spec_c.B + 0.05f*delta); break; // 0.0 to 1.0 in steps of 0.05
		default: assert(0);
		} // end switch
		if (ix == SM_RSCALE || ix == SM_HARDNESS || ix == SM_DENSITY) {
			int const gix(coll_id[MAT_SPHERE]);
			if (gix >= 0) {reanimate_group(gix, 0);}
		}
	}
};


// ************ Top-Level UI Hooks ************

unsigned selected_menu_ix(0);
hmap_kbd_menu_t hmap_menu(cur_brush_param);
voxel_edit_kbd_menu_t voxel_edit_menu(voxel_brush_params);
leaf_color_kbd_menu_t leaf_color_menu;
water_color_kbd_menu_t water_color_menu;
phys_weather_kbd_menu_t phys_weather_kbd_menu;
sphere_mat_kbd_menu_t sphere_mat_kbd_menu;


keyboard_menu_t *kbd_menus[] = {&hmap_menu, &voxel_edit_menu, &sphere_mat_kbd_menu, &leaf_color_menu, &water_color_menu, &phys_weather_kbd_menu};
unsigned const NUM_KBD_MENUS = sizeof(kbd_menus)/sizeof(kbd_menus[0]);


void next_selected_menu_ix() {

	unsigned num_enabled(0);
	for (unsigned i = 0; i < NUM_KBD_MENUS; ++i) {if (kbd_menus[i]->is_enabled()) {++num_enabled;}}
	if (num_enabled > 1) {++selected_menu_ix;} else {selected_menu_ix = 0;} // increment if there are multiple choices
	if (selected_menu_ix >= num_enabled) {selected_menu_ix = 0;} // wraparound
}


keyboard_menu_t *get_enabled_menu() {

	keyboard_menu_t *kbd_menu(NULL);
	unsigned num_enabled(0);
	
	for (unsigned i = 0; i < NUM_KBD_MENUS; ++i) {
		assert(kbd_menus[i]);

		if (kbd_menus[i]->is_enabled()) {
			kbd_menu = kbd_menus[i];
			if ((num_enabled++) == selected_menu_ix) {return kbd_menu;} // use this one (Nth enabled menu)
		}
	}
	if (kbd_menu == NULL) {return NULL;}
	assert(num_enabled > 0);
	selected_menu_ix = num_enabled - 1; // set to the last one (which was used)
	return kbd_menu;
}


bool ui_intercept_keyboard(unsigned char key, bool is_special) {

	if (!is_special) return 0; // only using special keys right now
	keyboard_menu_t *const kbd_menu(get_enabled_menu());
	if (kbd_menu == NULL) return 0; // no keyboard menu enabled
	bool const shift_pressed(is_shift_key_pressed());
	int const lr_mag(shift_pressed ? 10 : 1), ud_mag(shift_pressed ? 5 : 1); // move 10x/5x if shift is set
	
	switch (key) {
	case GLUT_KEY_UP:    for (int i = 0; i < ud_mag; ++i) {kbd_menu->next_control();} return 1;
	case GLUT_KEY_DOWN:  for (int i = 0; i < ud_mag; ++i) {kbd_menu->prev_control();} return 1;
	case GLUT_KEY_LEFT:  kbd_menu->change_value(-1*lr_mag); return 1;
	case GLUT_KEY_RIGHT: kbd_menu->change_value( 1*lr_mag); return 1;
	}
	return 0;
}

// if is_up_down=0, button and state and invalid
bool ui_intercept_mouse(int button, int state, int x, int y, bool is_up_down) {
	return 0; // do nothing (for now)
}


void show_onscreen_building_help() { // non-interactive menu
	ostringstream oss;
	oss << "wasd: Move" << endl << "e: Take Item" << endl << "q: Interact/Push Object" << endl << "r: Pull Object" << endl << "E or <Space>: Use Inventory Item" << endl
		<< "S: Toggle Room Light" << endl << "Left/Right Mouse Button: Flashlight" << endl << "Mouse Wheel: Switch Inventory Item" << endl << "Tab: Toggle Help";
	text_drawer text_draw;
	text_draw.begin_draw();
	text_draw.add_text(YELLOW, -0.01, 0.01, -0.02, oss.str().c_str(), MENU_TEXT_SIZE);
	text_draw.end_draw();
}

void draw_enabled_ui_menus() {
	if (show_scores && camera_in_building) {show_onscreen_building_help(); return;}
	keyboard_menu_t const *const kbd_menu(get_enabled_menu());
	if (kbd_menu) {kbd_menu->draw_controls();}
}


