// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/17/02


#ifndef _U_EVENT_H_
#define _U_EVENT_H_


#include <stdio.h>
#include <stdlib.h>

// STL include
#include <vector>
#include <assert.h>


unsigned const NUM_UE_TYPES     = 9;
unsigned const UE_MAX_PARAMS    = 4;
unsigned const MAX_U_EVENTS     = 100000;
unsigned const MAX_EVENT_FRAMES = 1000000;
char *const UEL_SAVE_NAME       = "ueventlist";

enum {UE_SRAND = 0, UE_RESIZE, UE_MBUTTON, UE_MMOTION, UE_KEYBOARD, UE_BREAK, UE_GOTO, UE_NULL, UE_KEYBOARD_SPECIAL, UE_KEYBOARD_UP};


int const ue_nparams[NUM_UE_TYPES] = {1, 2, 4, 2, 3, 0, 2, 0, 3};


struct uevent {

	int type, frame, params[UE_MAX_PARAMS];

	uevent() {}
	uevent(int type_, int frame_) : type(type_), frame(frame_) {}
};


// function prototypes - main loop
void display(void);
void resize(int x, int y);
void mouseButton(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void mousePassiveMotion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void keyboard2(int key, int x, int y);
void keyboard_up(unsigned char key, int x, int y);
void keyboard2_up(int key, int x, int y);
void proc_kbd_events();


// function prototypes - user events
int  read_ueventlist(char *arg);
int  save_ueventlist();
void uevent_advance_frame();
void add_uevent_srand(int rseed);
void add_uevent_resize(int x, int y);
void add_uevent_mbutton(int button, int state, int x, int y);
void add_uevent_mmotion(int x, int y);
void add_uevent_keyboard(unsigned char key, int x, int y);
void add_uevent_keyboard_up(unsigned char key, int x, int y);
void add_uevent_keyboard_special(int key, int x, int y);
int  check_event_ok();


#endif
