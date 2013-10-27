// 3D World - User Interface for Mesh Editing
// by Frank Gennari
// 10/26/13

#include "3DWorld.h"
#include "heightmap.h" // for hmap_brush_t


class keyboard_menu_t {

protected:
	bool enabled;
	unsigned num_controls, cur_control;

public:
	keyboard_menu_t(unsigned num_controls_) : enabled(0), num_controls(num_controls_), cur_control(0) {assert(num_controls > 0);}
	virtual ~keyboard_menu_t() {}
	bool is_enabled() const {return enabled;}
	void next_control() {cur_control = (cur_control == num_controls-1) ? 0  : (cur_control+1);}
	void prev_control() {cur_control = (cur_control == 0) ? (num_controls-1): (cur_control-1);}
	virtual void change_value(int delta) = 0;
};


enum {HMAP_DELAY=0, HMAP_CTR_SHAPE, HMAP_CTR_RADIUS, HMAP_CTR_DELTA, HMAP_NUM_CTR};

class hmap_kbd_menu_t : public keyboard_menu_t {

	tex_mod_map_manager_t::hmap_brush_t hmap_brush;
	unsigned hmap_delay_ms;

public:
	hmap_kbd_menu_t() : keyboard_menu_t(HMAP_NUM_CTR), hmap_delay_ms(0) {}

	virtual void change_value(int delta) {
		switch (cur_control) {
		case HMAP_DELAY:
			// modify hmap_delay
			break;
		case HMAP_CTR_SHAPE:
			// modify hmap_brush.shape
			break;
		case HMAP_CTR_RADIUS:
			// modify hmap_brush.radius
			break;
		case HMAP_CTR_DELTA:
			// modify hmap_brush.delta
			break;
		default: assert(0);
		}
		// FIXME: WRITE
	}
};

hmap_kbd_menu_t hmap_menu;


bool ui_intercept_keyboard(unsigned char key, bool is_special) {

	if (!is_special) return 0; // only using special keys right now
	keyboard_menu_t *kbd_menu(NULL);
	if (hmap_menu.is_enabled()) {kbd_menu = &hmap_menu;}
	// add more handlers here
	if (kbd_menu == NULL) return 0; // no keyboard menu enabled
	
	switch (key) {
	case GLUT_KEY_UP:    kbd_menu->prev_control(  ); return 1;
	case GLUT_KEY_DOWN:  kbd_menu->next_control(  ); return 1;
	case GLUT_KEY_LEFT:  kbd_menu->change_value(-1); return 1; // more than one if modifier (shift?) is set?
	case GLUT_KEY_RIGHT: kbd_menu->change_value( 1); return 1; // more than one if modifier (shift?) is set?
	}
	return 0;
}


// if is_up_down=0, button and state and invalid
bool ui_intercept_mouse(int button, int state, int x, int y, bool is_up_down) {

	return 0; // do nothing (for now)
}


