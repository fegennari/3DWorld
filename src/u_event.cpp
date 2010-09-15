// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/17/02


#include "u_event.h"
#include "iostream"

using std::vector;
using std::cout;
using std::endl;

// Global Variables
int read_eventlist(0), make_eventlist(0), curr_event(0), n_events(0), n_frames(0), frame_counter(0);
vector<uevent> eventlist;

bool open_file(FILE *&fp, char const *const fn, std::string const &file_type, char const *const mode="r");



int read_ueventlist(char *arg) {

	if (strcmp(arg, "-uel") == 0) { // record user event list
		make_eventlist = 1;
		return 2;
	}
	uevent event;
	FILE *fp;
	if (!open_file(fp, arg, "input user eventlist")) return 0;

	if (fscanf(fp, "%u%u", &n_events, &n_frames) != 2) {
		cout << "Error reading user event list header." << endl;
		fclose(fp);
		return 0;
	}
	while (fscanf(fp, "%i %i", &(event.type), &(event.frame)) == 2) {
		if (event.type < (int)NUM_UE_TYPES) {
			unsigned const nparams(std::min((unsigned)ue_nparams[event.type], UE_MAX_PARAMS));

			for (unsigned i = 0; i < nparams; ++i) {
				if (fscanf(fp, "%u", &(event.params[i])) != 1) {
					cout << "Error reading user event file: event # " << eventlist.size() << ", type " << event.type << ", param " << i << endl;
					eventlist.clear();
					fclose(fp);
					return 0;
				}
			}
			eventlist.push_back(event);
		}
	}
	read_eventlist = 1;
	fclose(fp);
	return 1;
}


int save_ueventlist() {

	if (make_eventlist == 0) return 0;
	FILE *fp;
	if (!open_file(fp, UEL_SAVE_NAME, "output user eventlist"), "w") return 0;

	if (!fprintf(fp, "%zi %u\n", eventlist.size(), frame_counter)) {
		cout << "Error writing user event list header." << endl;
		fclose(fp);
		return 0;
	}
	for (unsigned i = 0; i < eventlist.size(); ++i) {
		unsigned const nparams(std::min((unsigned)ue_nparams[eventlist[i].type], UE_MAX_PARAMS));
		fprintf(fp, "%i %i", eventlist[i].type, eventlist[i].frame);

		for (unsigned j = 0; j < nparams; ++j) {
			fprintf(fp, " %i", eventlist[i].params[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	cout << "Ueventlist saved as '" << UEL_SAVE_NAME << "'." << endl;
	return 1;
}


void uevent_advance_frame() {

	if (read_eventlist == 0) {
		++frame_counter;
		return;
	}
	while (curr_event < (int)eventlist.size()) {
		if (eventlist[curr_event].frame > frame_counter) {
			++frame_counter;
			return;
		}
		if(eventlist[curr_event].type < (int)NUM_UE_TYPES) {
			int *params(eventlist[curr_event].params);

			switch (eventlist[curr_event].type) {
			case UE_SRAND:
				srand(params[0]);
				break;
			case UE_RESIZE:
				resize(params[0], params[1]);
				break;
			case UE_MBUTTON:
				mouseButton(params[0], params[1], params[2], params[3]);
				break;
			case UE_MMOTION:
				mouseMotion(params[0], params[1]);
				break;
			case UE_KEYBOARD:
				keyboard((unsigned char)(params[0]), params[1], params[2]);
				break;
			case UE_KEYBOARD_UP:
				keyboard_up((unsigned char)(params[0]), params[1], params[2]);
				break;
			case UE_KEYBOARD_SPECIAL:
				keyboard2((unsigned char)(params[0]), params[1], params[2]);
				break;
			case UE_BREAK:
				cout << "User break at frame " << frame_counter << endl;
				assert(0);
				break;
			case UE_GOTO: // must be manually added by user
				if (params[0] == 0) { // goto frame
					if (frame_counter == params[1]) break;
					frame_counter = params[1];

					if (frame_counter > n_frames + 1) {
						curr_event = (int)eventlist.size();
						return;
					}
					unsigned i;
					for (i = 0; i < eventlist.size(); ++i) {
						if (eventlist[i].frame == frame_counter) break;
						if (eventlist[i].frame > frame_counter) {--i; break;}
					}
					curr_event = (int)i;
					return;
				}
				else { // goto event
					if (params[1] >= (int)eventlist.size()) {
						curr_event = (int)eventlist.size();
						break;
					}
					curr_event    = params[1];
					frame_counter = eventlist[curr_event].frame;
					return;
				}
				break;
			case UE_NULL:
				break;
			}
		}
		++curr_event;
	}
	++frame_counter;
}


void add_uevent_srand(int rseed) {

	if (read_eventlist == 1) return;
	srand(rseed);
	if (!check_event_ok())   return;
	uevent event(UE_SRAND, frame_counter);
	event.params[0] = rseed;
	eventlist.push_back(event);
}


void add_uevent_resize(int x, int y) {

	if (!check_event_ok()) return;
	uevent event(UE_RESIZE, frame_counter);
	event.params[0] = x;
	event.params[1] = y;
	eventlist.push_back(event);
}


void add_uevent_mbutton(int button, int state, int x, int y) {

	if (!check_event_ok()) return;
	uevent event(UE_MBUTTON, frame_counter);
	event.params[0] = button;
	event.params[1] = state;
	event.params[2] = x;
	event.params[3] = y;
	eventlist.push_back(event);
}


void add_uevent_mmotion(int x, int y) {

	if (!check_event_ok()) return;
	uevent event(UE_MMOTION, frame_counter);
	event.params[0] = x;
	event.params[1] = y;
	eventlist.push_back(event);
}


void add_uevent_keyboard_def(unsigned char key, int x, int y, int type) {

	if (!check_event_ok()) return;
	uevent event(type, frame_counter);
	event.params[0] = (int)key;
	event.params[1] = x;
	event.params[2] = y;
	eventlist.push_back(event);
}


void add_uevent_keyboard(unsigned char key, int x, int y) {

	add_uevent_keyboard_def(key, x, y, UE_KEYBOARD);
}


void add_uevent_keyboard_up(unsigned char key, int x, int y) {

	add_uevent_keyboard_def(key, x, y, UE_KEYBOARD_UP);
}


void add_uevent_keyboard_special(int key, int x, int y) {

	add_uevent_keyboard_def(key, x, y, UE_KEYBOARD_SPECIAL);
}


int check_event_ok() {

	return (make_eventlist && frame_counter <= (int)MAX_EVENT_FRAMES && eventlist.size() <= (size_t)MAX_U_EVENTS);
}





