// 3D World - User Interface for Mesh Editing
// by Frank Gennari
// 10/26/13

#include "3DWorld.h"
#include "function_registry.h"
#include "heightmap.h" // for hmap_brush_t

using namespace std;


class keyboard_menu_t {

protected:
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
		bool const bitmap = 0;
		float const size = 1.0;
		draw_text(-0.01, 0.01-0.0014*(num_controls - control_ix), -0.02, oss.str().c_str(), size, bitmap);
	}
	virtual void draw_one_control(unsigned control_ix) const = 0;

public:
	keyboard_menu_t(unsigned num_controls_) : num_controls(num_controls_), cur_control(0) {assert(num_controls > 0);}
	virtual ~keyboard_menu_t() {}
	virtual bool is_enabled() const = 0;
	void next_control() {cur_control = (cur_control == num_controls-1) ? 0  : (cur_control+1);}
	void prev_control() {cur_control = (cur_control == 0) ? (num_controls-1): (cur_control-1);}
	virtual void change_value(int delta) = 0;

	void draw_controls() const {
		for (unsigned i = 0; i < num_controls; ++i) {draw_one_control(i);}
	}
};


// ************ Tiled Terrain Heightmap Brush/Editing ************

enum {HMAP_DELAY=0, HMAP_CTR_SHAPE, HMAP_CTR_RADIUS, HMAP_CTR_DELTA, HMAP_NUM_CTR};
string const hmap_ctr_names[HMAP_NUM_CTR] = {"Placement Delay", "Brush Shape", "Brush Radius", "Brush Delta"};
string const brush_shapes[NUM_BSHAPES] = {"Constant Square", "Constant Circle", "Linear Circle", "Quadratic Circle", "Cosine Circle", "Flat Square", "Flat Circle"};

extern int show_scores, world_mode;
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
			spos = (brush_param.shape/float(BSHAPE_COSINE));
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
	hmap_kbd_menu_t(hmap_brush_param_t &bp) : keyboard_menu_t(HMAP_NUM_CTR), brush_param(bp), max_radius_exp(0) {
		for (unsigned sz = 1; sz < get_tile_size(); sz <<= 1) {++max_radius_exp;}
	}
	virtual bool is_enabled() const {return (show_scores && inf_terrain_fire_mode && world_mode == WMODE_INF_TERRAIN);}

	virtual void change_value(int delta) {
		switch (cur_control) {
		case HMAP_DELAY:
			brush_param.delay = max(0, min(10, ((int)brush_param.delay + delta)));
			break;
		case HMAP_CTR_SHAPE:
			brush_param.shape = max(0, min((int)BSHAPE_COSINE, (brush_param.shape + delta)));
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


// ************ Leaf Colors ************

enum {TREE_COLOR_VAR=0, LEAF_COLOR_VAR, LEAF_RED_COMP, LEAF_GREEN_COMP, LEAF_BLUE_COMP, NUM_LEAF_CONT};
string const leaf_ctr_names[NUM_LEAF_CONT] = {"Tree Color Variance", "Leaf Color Variance", "Leaf Red Component", "Leaf Green Component", "Leaf Blue Component"};

extern int leaf_color_changed, game_mode;
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
	leaf_color_kbd_menu_t() : keyboard_menu_t(NUM_LEAF_CONT) {}
	virtual bool is_enabled() const {return (show_scores && !game_mode && world_mode == WMODE_GROUND);}

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


// ************ Top-Level UI Hooks ************

hmap_kbd_menu_t hmap_menu(cur_brush_param);
leaf_color_kbd_menu_t leaf_color_menu;


unsigned const NUM_KBD_MENUS = 2;
keyboard_menu_t *kbd_menus[NUM_KBD_MENUS] = {&hmap_menu, &leaf_color_menu};


bool ui_intercept_keyboard(unsigned char key, bool is_special) {

	if (!is_special) return 0; // only using special keys right now
	keyboard_menu_t *kbd_menu(NULL);
	
	for (unsigned i = 0; i < NUM_KBD_MENUS; ++i) { // Note: doesn't handle more than one enabled menu
		assert(kbd_menus[i]);
		if (kbd_menus[i]->is_enabled()) {assert(!kbd_menu); kbd_menu = kbd_menus[i];}
	}
	if (kbd_menu == NULL) return 0; // no keyboard menu enabled
	
	switch (key) {
	case GLUT_KEY_UP:    kbd_menu->next_control(  ); return 1;
	case GLUT_KEY_DOWN:  kbd_menu->prev_control(  ); return 1;
	case GLUT_KEY_LEFT:  kbd_menu->change_value(-1); return 1; // more than one if modifier (shift?) is set?
	case GLUT_KEY_RIGHT: kbd_menu->change_value( 1); return 1; // more than one if modifier (shift?) is set?
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


