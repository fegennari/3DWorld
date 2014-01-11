// 3D World - User Interface for Mesh Editing
// by Frank Gennari
// 10/26/13

#include "3DWorld.h"
#include "function_registry.h"
#include "heightmap.h" // for hmap_brush_t
#include "voxels.h" // for voxel_brush_params_t

using namespace std;

bool  const MENU_BITMAP_TEXT = 0;
float const MENU_TEXT_SIZE   = 1.0;


class keyboard_menu_t {

protected:
	string title;
	unsigned num_controls, cur_control;

	void draw_one_control_text(unsigned control_ix, string const &name, string const &cur_value, float slider_pos) const {
		//assert(slider_pos >= 0.0 && slider_pos <= 1.0);
		slider_pos = CLIP_TO_01(slider_pos);
		assert(control_ix < num_controls);
		bool const selected(control_ix == cur_control);
		(selected ? ORANGE : YELLOW).do_glColor();
		unsigned const ndiv = 20;
		unsigned const pos(round_fp((ndiv-1)*slider_pos));
		ostringstream oss;
		for (unsigned n = 0; n < pos; ++n) {oss << "-";}
		oss << "+";
		for (unsigned n = pos+1; n < ndiv; ++n) {oss << "-";}
		oss << "  " << name << ": " << cur_value;
		draw_text(-0.01, 0.01-0.0014*(num_controls - control_ix), -0.02, oss.str().c_str(), MENU_TEXT_SIZE, MENU_BITMAP_TEXT);
	}
	virtual void draw_one_control(unsigned control_ix) const = 0;

public:
	keyboard_menu_t(unsigned num_controls_, string const &title_) : num_controls(num_controls_), cur_control(0), title(title_) {assert(num_controls > 0);}
	virtual ~keyboard_menu_t() {}
	virtual bool is_enabled() const = 0;
	void next_control() {cur_control = (cur_control == num_controls-1) ? 0  : (cur_control+1);}
	void prev_control() {cur_control = (cur_control == 0) ? (num_controls-1): (cur_control-1);}
	virtual void change_value(int delta) = 0;

	void draw_controls() const {
		if (!title.empty()) {
			YELLOW.do_glColor();
			draw_text(-0.01, 0.01, -0.02, title.c_str(), MENU_TEXT_SIZE, MENU_BITMAP_TEXT);
		}
		for (unsigned i = 0; i < num_controls; ++i) {draw_one_control(i);}
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
			value << brush_param.get_radius();
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
			value << get_voxel_brush_step()*brush.radius << " (" << brush.radius << " units)";
			break;
		case VOXEL_WEIGHT:
			spos = 0.5*(brush.weight_exp + MAX_VB_WEIGHT_EXP)/MAX_VB_WEIGHT_EXP;
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


// ************ Leaf Colors ************

enum {TREE_COLOR_VAR=0, LEAF_COLOR_VAR, LEAF_RED_COMP, LEAF_GREEN_COMP, LEAF_BLUE_COMP, NUM_LEAF_CONT};
string const leaf_ctr_names[NUM_LEAF_CONT] = {"Tree Color Variance", "Leaf Color Variance", "Leaf Red Component", "Leaf Green Component", "Leaf Blue Component"};

extern int leaf_color_changed;
extern float leaf_color_coherence, tree_color_coherence; // [0.0, 1.0] in steps of 0.1
extern colorRGBA leaf_base_color; // can set R and G in [-1.0, 1.0] in steps of 0.1

class leaf_color_kbd_menu_t : public keyboard_menu_t {

	virtual void draw_one_control(unsigned control_ix) const {
		assert(control_ix < NUM_LEAF_CONT);
		ostringstream value;
		value.precision(1);
		value << fixed; // fixed precision in units of 0.1
		float spos(0.0);

		switch (control_ix) {
		case TREE_COLOR_VAR:
			spos = tree_color_coherence;
			value << spos;
			break;
		case LEAF_COLOR_VAR:
			spos = (1.0 - leaf_color_coherence);
			value << spos;
			break;
		default:
			spos = leaf_base_color[control_ix-LEAF_RED_COMP];
			value << spos;
			spos = 0.5*(spos + 1.0); // center around 0.0
		}
		draw_one_control_text(control_ix, leaf_ctr_names[control_ix], value.str(), spos);
	}

public:
	leaf_color_kbd_menu_t() : keyboard_menu_t(NUM_LEAF_CONT, "Leaf Colors") {}
	virtual bool is_enabled() const {return (show_scores && !game_mode && !voxel_editing && world_mode == WMODE_GROUND);}

	virtual void change_value(int delta) {
		switch (cur_control) {
		case TREE_COLOR_VAR:
			tree_color_coherence = max(0.0f, (tree_color_coherence + 0.1f*delta)); // no maximum
			break;
		case LEAF_COLOR_VAR:
			leaf_color_coherence = CLIP_TO_01(leaf_color_coherence - 0.1f*delta); // delta is backwards
			break;
		default:
			leaf_base_color[cur_control-LEAF_RED_COMP] = CLIP_TO_pm1(leaf_base_color[cur_control-LEAF_RED_COMP] + 0.1f*delta);
			break;
		}
		leaf_color_changed = 1;
	}
};


// ************ Water Colors/Properties ************

colorRGB uw_atten_max(WHITE), uw_atten_scale(BLACK);
water_params_t water_params;

void calc_uw_atten_colors() {

	blend_color(uw_atten_scale, colorRGB(0.9, 1.0, 1.5), colorRGB(1.5, 0.9, 0.5), water_params.mud); // blend in mud color
	//UNROLL_3X(uw_atten_max[i_] = CLIP_TO_01(1.0f - 0.03f/uw_atten_scale[i_]);)
	uw_atten_scale *= 0.05 + 0.95*water_params.alpha;
}

enum {WATERP_ALPHA=0, WATERP_MUD, WATERP_BRIGHT, WATERP_REFLECT, WATERP_GREEN, NUM_WATER_CONT};
string const water_ctr_names[NUM_WATER_CONT] = {"Alpha Scale", "Mud Content", "Brightness", "Reflectivity", "Green Hue"};

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
		default:
			assert(0);
		}
		draw_one_control_text(control_ix, water_ctr_names[control_ix], value.str(), spos);
	}

public:
	water_color_kbd_menu_t() : keyboard_menu_t(NUM_LEAF_CONT, "Water Colors") {}
	virtual bool is_enabled() const {return (show_scores && !game_mode && !inf_terrain_fire_mode && world_mode == WMODE_INF_TERRAIN);}

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
		default:
			assert(0);
		}
		calc_uw_atten_colors();
	}
};


// ************ Top-Level UI Hooks ************

hmap_kbd_menu_t hmap_menu(cur_brush_param);
voxel_edit_kbd_menu_t voxel_edit_menu(voxel_brush_params);
leaf_color_kbd_menu_t leaf_color_menu;
water_color_kbd_menu_t water_color_menu;


keyboard_menu_t *kbd_menus[] = {&hmap_menu, &voxel_edit_menu, &leaf_color_menu, &water_color_menu};
unsigned const NUM_KBD_MENUS = sizeof(kbd_menus)/sizeof(kbd_menus[0]);


bool ui_intercept_keyboard(unsigned char key, bool is_special) {

	if (!is_special) return 0; // only using special keys right now
	keyboard_menu_t *kbd_menu(NULL);
	
	for (unsigned i = 0; i < NUM_KBD_MENUS; ++i) { // Note: doesn't handle more than one enabled menu
		assert(kbd_menus[i]);
		if (kbd_menus[i]->is_enabled()) {assert(!kbd_menu); kbd_menu = kbd_menus[i];}
	}
	if (kbd_menu == NULL) return 0; // no keyboard menu enabled
	int const change_mag((glutGetModifiers() & GLUT_ACTIVE_SHIFT) ? 10 : 1); // move 10x if shift is set
	
	switch (key) {
	case GLUT_KEY_UP:    kbd_menu->next_control(); return 1;
	case GLUT_KEY_DOWN:  kbd_menu->prev_control(); return 1;
	case GLUT_KEY_LEFT:  kbd_menu->change_value(-1*change_mag); return 1;
	case GLUT_KEY_RIGHT: kbd_menu->change_value( 1*change_mag); return 1;
	}
	return 0;
}


// if is_up_down=0, button and state and invalid
bool ui_intercept_mouse(int button, int state, int x, int y, bool is_up_down) {

	return 0; // do nothing (for now)
}


void draw_enabled_ui_menus() {

	bool drawn(0);

	for (unsigned i = 0; i < NUM_KBD_MENUS; ++i) { // Note: doesn't handle more than one enabled menu
		assert(kbd_menus[i]);
		if (!kbd_menus[i]->is_enabled()) continue;
		assert(!drawn);
		kbd_menus[i]->draw_controls();
		drawn = 1;
	}
}


