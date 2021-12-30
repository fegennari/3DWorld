// 3D World - Building Gameplay Logic
// by Frank Gennari 12/29/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "openal_wrap.h"

float const THROW_VELOCITY = 0.0050;
float const ALERT_THRESH   = 0.08; // min sound alert level for AIs

bool do_room_obj_pickup(0), use_last_pickup_object(0), show_bldg_pickup_crosshair(0), player_near_toilet(0), player_in_elevator(0);
bool city_action_key(0), can_do_building_action(0);
int can_pickup_bldg_obj(0);
float office_chair_rot_rate(0.0), cur_building_sound_level(0.0);
carried_item_t player_held_object;
bldg_obj_type_t bldg_obj_types[NUM_ROBJ_TYPES];
vector<sphere_t> cur_sounds; // radius = sound volume

extern bool camera_in_building, player_is_hiding;
extern int window_width, window_height, display_framerate, display_mode, game_mode, building_action_key;
extern float fticks, CAMERA_RADIUS;
extern double tfticks, camera_zh;
extern building_params_t global_building_params;


void place_player_at_xy(float xval, float yval);
room_object_t get_dresser_middle(room_object_t const &c);
room_object_t get_desk_drawers_part(room_object_t const &c);
void show_key_icon();
bool is_shirt_model(room_object_t const &obj);
bool is_pants_model(room_object_t const &obj);

bool in_building_gameplay_mode() {return (game_mode == 2);} // replaces dodgeball mode

// object types/pickup

void setup_bldg_obj_types() {
	static bool was_setup(0);
	if (was_setup) return; // nothing to do
	was_setup = 1;
	// player_coll, ai_coll, pickup, attached, is_model, lg_sm, value, weight, name [capacity]
	//                                                pc ac pu at im ls value  weight  name
	bldg_obj_types[TYPE_TABLE     ] = bldg_obj_type_t(1, 1, 1, 0, 0, 1, 70.0,  40.0,  "table");
	bldg_obj_types[TYPE_CHAIR     ] = bldg_obj_type_t(0, 1, 1, 0, 0, 1, 50.0,  25.0,  "chair"); // skip player collisions because they can be in the way and block the path in some rooms
	bldg_obj_types[TYPE_STAIR     ] = bldg_obj_type_t(1, 0, 0, 1, 0, 1, 0.0,   0.0,   "stair");
	bldg_obj_types[TYPE_STAIR_WALL] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "stairs wall");
	bldg_obj_types[TYPE_ELEVATOR  ] = bldg_obj_type_t(1, 1, 0, 1, 0, 0, 0.0,   0.0,   "elevator");
	bldg_obj_types[TYPE_LIGHT     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 0, 40.0,  5.0,   "light");
	bldg_obj_types[TYPE_RUG       ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 50.0,  20.0,  "rug");
	bldg_obj_types[TYPE_PICTURE   ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 100.0, 1.0,   "picture"); // should be random value
	bldg_obj_types[TYPE_WBOARD    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 50.0,  25.0,  "whiteboard");
	bldg_obj_types[TYPE_BOOK      ] = bldg_obj_type_t(0, 0, 1, 0, 0, 3, 10.0,  1.0,   "book");
	bldg_obj_types[TYPE_BCASE     ] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 150.0, 100.0, "bookcase"); // Note: can't pick up until bookcase can be expanded and books taken off
	bldg_obj_types[TYPE_TCAN      ] = bldg_obj_type_t(0, 1, 1, 0, 0, 2, 12.0,  2.0,   "trashcan"); // skip player collisions because they can be in the way and block the path in some rooms
	bldg_obj_types[TYPE_DESK      ] = bldg_obj_type_t(1, 1, 0, 0, 0, 3, 100.0, 80.0,  "desk"); // drawers are small items
	bldg_obj_types[TYPE_BED       ] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 300.0, 200.0, "bed"); // pillows are small, and the rest is large
	bldg_obj_types[TYPE_WINDOW    ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 0.0,   0.0,   "window");
	bldg_obj_types[TYPE_BLOCKER   ] = bldg_obj_type_t(0, 0, 0, 0, 0, 0, 0.0,   0.0,   "<blocker>");  // not a drawn object; block other objects, but not the player or AI
	bldg_obj_types[TYPE_COLLIDER  ] = bldg_obj_type_t(1, 1, 0, 0, 0, 0, 0.0,   0.0,   "<collider>"); // not a drawn object; block the player and AI
	bldg_obj_types[TYPE_CUBICLE   ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 500.0, 250.0, "cubicle"); // skip collisions because they have their own colliders
	bldg_obj_types[TYPE_STALL     ] = bldg_obj_type_t(1, 1, 1, 1, 0, 1, 40.0,  20.0,  "bathroom divider"); // can pick up short sections of bathroom stalls (urinal dividers)
	bldg_obj_types[TYPE_SIGN      ] = bldg_obj_type_t(0, 0, 1, 0, 0, 3, 10.0,  1.0,   "sign");
	bldg_obj_types[TYPE_COUNTER   ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "kitchen counter");
	bldg_obj_types[TYPE_CABINET   ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 0.0,   0.0,   "kitchen cabinet");
	bldg_obj_types[TYPE_KSINK     ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "kitchen sink");
	bldg_obj_types[TYPE_BRSINK    ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "bathroom sink"); // for office building bathrooms
	bldg_obj_types[TYPE_PLANT     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 3, 18.0,  8.0,   "potted plant");
	bldg_obj_types[TYPE_DRESSER   ] = bldg_obj_type_t(1, 1, 0, 0, 0, 3, 120.0, 110.0, "dresser"); // Note: can't pick up until drawers can be opened and items removed from them
	bldg_obj_types[TYPE_NIGHTSTAND] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 60.0,  45.0,  "nightstand");
	bldg_obj_types[TYPE_FLOORING  ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 0.0,   0.0,   "flooring");
	// closets can't be picked up, but they can block a pickup; marked as large because small objects are not modified; marked as is_model because closets can contain lamps
	bldg_obj_types[TYPE_CLOSET    ] = bldg_obj_type_t(1, 1, 1, 1, 1, 1, 0.0,   0.0,   "closet");
	bldg_obj_types[TYPE_WALL_TRIM ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "wall trim");
	bldg_obj_types[TYPE_RAILING   ] = bldg_obj_type_t(1, 0, 0, 1, 0, 2, 0.0,   0.0,   "railing");
	bldg_obj_types[TYPE_CRATE     ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 10.0,  12.0,  "crate"); // should be random value
	bldg_obj_types[TYPE_BOX       ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 5.0,   8.0,   "box");   // should be random value
	bldg_obj_types[TYPE_MIRROR    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 40.0,  15.0,  "mirror");
	bldg_obj_types[TYPE_SHELVES   ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 0.0,   0.0,   "shelves");
	bldg_obj_types[TYPE_KEYBOARD  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 15.0,  2.0,   "keyboard");
	bldg_obj_types[TYPE_SHOWER    ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "shower");
	bldg_obj_types[TYPE_RDESK     ] = bldg_obj_type_t(1, 1, 0, 0, 0, 1, 800.0, 300.0, "reception desk");
	bldg_obj_types[TYPE_BOTTLE    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 1.0,   1.0,   "bottle");
	bldg_obj_types[TYPE_WINE_RACK ] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 75.0,  40.0,  "wine rack");
	bldg_obj_types[TYPE_COMPUTER  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 500.0, 20.0,  "computer");
	bldg_obj_types[TYPE_MWAVE     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 100.0, 50.0,  "microwave oven");
	bldg_obj_types[TYPE_PAPER     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.0,   0.0,   "sheet of paper"); // will have a random value that's often 0
	bldg_obj_types[TYPE_BLINDS    ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 50.0,  7.0,   "window blinds");
	bldg_obj_types[TYPE_PEN       ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.10,  0.02,  "pen");
	bldg_obj_types[TYPE_PENCIL    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.10,  0.02,  "pencil");
	bldg_obj_types[TYPE_PAINTCAN  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 12.0,  8.0,   "paint can");
	bldg_obj_types[TYPE_LG_BALL   ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 15.0,  1.2,   "ball");
	bldg_obj_types[TYPE_HANGER_ROD] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 10.0,  5.0,   "hanger rod");
	bldg_obj_types[TYPE_DRAIN     ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "drain pipe");
	bldg_obj_types[TYPE_MONEY     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 20.0,  0.0,   "pile of money"); // $20 bills
	bldg_obj_types[TYPE_PHONE     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 200.0, 0.1,   "cell phone");
	bldg_obj_types[TYPE_TPROLL    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.25,  0.1,   "TP roll", 200);
	bldg_obj_types[TYPE_SPRAYCAN  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 2.0,   1.0,   "spray paint", 5000);
	bldg_obj_types[TYPE_MARKER    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.20,  0.05,  "marker",      10000);
	bldg_obj_types[TYPE_BUTTON    ] = bldg_obj_type_t(0, 0, 1, 1, 0, 2, 1.0,   0.05,  "button");
	bldg_obj_types[TYPE_CRACK     ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "crack");
	bldg_obj_types[TYPE_SWITCH    ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "switch");
	bldg_obj_types[TYPE_PLATE     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 6.0,   0.25,  "plate");
	bldg_obj_types[TYPE_LAPTOP    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 600.0, 8.0,   "laptop");
	bldg_obj_types[TYPE_FPLACE    ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   2000.0,"fireplace");
	bldg_obj_types[TYPE_LBASKET   ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 12.0,  2.0,   "laundry basket");
	bldg_obj_types[TYPE_WHEATER   ] = bldg_obj_type_t(1, 1, 0, 1, 0, 2, 300.0, 500.0, "water heater");
	bldg_obj_types[TYPE_TAPE      ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 2.0,   0.4,   "duct tape", 1000);
	// player_coll, ai_coll, pickup, attached, is_model, lg_sm, value, weight, name [capacity]
	// 3D models
	bldg_obj_types[TYPE_TOILET    ] = bldg_obj_type_t(1, 1, 1, 1, 1, 0, 120.0, 88.0,  "toilet");
	bldg_obj_types[TYPE_SINK      ] = bldg_obj_type_t(1, 1, 1, 1, 1, 0, 80.0,  55.0,  "sink");
	bldg_obj_types[TYPE_TUB       ] = bldg_obj_type_t(1, 1, 0, 1, 1, 1, 250.0, 200.0, "bathtub");
	bldg_obj_types[TYPE_FRIDGE    ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 700.0, 300.0, "refrigerator"); // no pickup, too large and may want to keep it for future hunger bar
	bldg_obj_types[TYPE_STOVE     ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 400.0, 150.0, "stove");
	bldg_obj_types[TYPE_TV        ] = bldg_obj_type_t(1, 1, 1, 0, 1, 1, 400.0, 70.0,  "TV");
	bldg_obj_types[TYPE_MONITOR   ] = bldg_obj_type_t(1, 1, 1, 0, 1, 1, 250.0, 15.0,  "computer monitor");
	bldg_obj_types[TYPE_COUCH     ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 600.0, 300.0, "couch");
	bldg_obj_types[TYPE_OFF_CHAIR ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 150.0, 60.0,  "office chair");
	bldg_obj_types[TYPE_URINAL    ] = bldg_obj_type_t(1, 1, 1, 1, 1, 0, 100.0, 80.0,  "urinal");
	bldg_obj_types[TYPE_LAMP      ] = bldg_obj_type_t(0, 0, 1, 0, 1, 0, 25.0,  12.0,  "lamp");
	bldg_obj_types[TYPE_WASHER    ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 300.0, 150.0, "washer");
	bldg_obj_types[TYPE_DRYER     ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 300.0, 160.0, "dryer");
	// keys are special because they're potentially either a small object or an object model (in a drawer)
	bldg_obj_types[TYPE_KEY       ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.0,   0.05,  "room key"); // drawn as an object, not a model
	bldg_obj_types[TYPE_HANGER    ] = bldg_obj_type_t(0, 0, 1, 0, 1, 2, 0.25,  0.05,  "clothes hanger");
	bldg_obj_types[TYPE_CLOTHES   ] = bldg_obj_type_t(0, 0, 1, 0, 1, 2, 10.0,  0.25,  "clothes"); // teeshirt, shirt, pants, etc.
	//                                                pc ac pu at im ls value  weight  name [capacity]
}

bldg_obj_type_t const &get_room_obj_type(room_object_t const &obj) {
	assert(obj.type < NUM_ROBJ_TYPES);
	return bldg_obj_types[obj.type];
}
float carried_item_t::get_remaining_capacity_ratio() const {
	unsigned const capacity(get_room_obj_type(*this).capacity);
	return ((capacity == 0) ? 1.0 : (1.0 - float(use_count)/float(capacity))); // Note: zero capacity is unlimited and ratio returned is always 1.0
}

bldg_obj_type_t get_taken_obj_type(room_object_t const &obj) {
	if (obj.type == TYPE_PICTURE && (obj.flags & RO_FLAG_TAKEN1)) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 20.0, 6.0, "picture frame");} // second item to take from picture
	if (obj.type == TYPE_TPROLL  && (obj.flags & RO_FLAG_TAKEN1)) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 6.0,  0.5, "toilet paper holder");} // second item to take from tproll

	if (obj.type == TYPE_BED) { // player_coll, ai_coll, pickup, attached, is_model, lg_sm, value, weight, name
		if (obj.flags & RO_FLAG_TAKEN2) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 250.0, 80.0, "mattress"  );} // third item to take from bed
		if (obj.flags & RO_FLAG_TAKEN1) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 80.0,  4.0,  "bed sheets");} // second item to take from bed
		return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 20.0, 1.0, "pillow"); // first item to take from bed
	}
	if (obj.type == TYPE_PLANT && !(obj.flags & RO_FLAG_ADJ_BOT)) { // plant not on a table/desk
		if (obj.flags & RO_FLAG_TAKEN2) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 10.0, 10.0, "plant pot");} // third item to take
		if (obj.flags & RO_FLAG_TAKEN1) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 1.0,  10.0, "dirt"     );} // second item to take
		return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 25.0, 5.0, "plant"); // first item to take
	}
	if (obj.type == TYPE_COMPUTER && (obj.flags & RO_FLAG_WAS_EXP)) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 100.0, 20.0, "old computer");}
	if (obj.type == TYPE_BOX      && (obj.flags & RO_FLAG_OPEN   )) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2,   0.0, 0.05, "opened box"  );}
	if (obj.type == TYPE_CRATE    && (obj.flags & RO_FLAG_OPEN   )) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2,   2.0, 0.5,  "opened crate");}
	if (obj.type == TYPE_TV       && (obj.flags & RO_FLAG_BROKEN )) {return bldg_obj_type_t(1, 1, 1, 0, 1, 1,  20.0, 70.0, "broken TV"   );}
	if (obj.type == TYPE_MONITOR  && (obj.flags & RO_FLAG_BROKEN )) {return bldg_obj_type_t(1, 1, 1, 0, 1, 1,  10.0, 15.0, "broken computer monitor");}

	if (obj.type == TYPE_BOTTLE) {
		bottle_params_t const &bparams(bottle_params[obj.get_bottle_type()]);
		bldg_obj_type_t type(0, 0, 1, 0, 0, 2,  bparams.value, 1.0, bparams.name);

		if (obj.is_bottle_empty()) {
			type.name    = "empty " + type.name;
			type.weight *= 0.25;
			type.value   = 0.0;
		}
		return type;
	}
	// default value, name may be modified below
	bldg_obj_type_t type(get_room_obj_type(obj));

	if (obj.type == TYPE_LG_BALL) {
		type.name = ((obj.item_flags & 1) ? "basketball" : "soccer ball"); // use a more specific type name; all other fields are shared across balls
	}
	else if (obj.type == TYPE_CLOTHES) {
		if      (is_shirt_model(obj)) {type.name = "shirt";}
		else if (is_pants_model(obj)) {type.name = "pants";}
	}
	return type;
}
rand_gen_t rgen_from_obj(room_object_t const &obj) {
	rand_gen_t rgen;
	rgen.set_state(12345*abs(obj.x1()), 67890*abs(obj.y1()));
	return rgen;
}
float get_obj_value(room_object_t const &obj) {
	float value(get_taken_obj_type(obj).value);
	if (obj.type == TYPE_CRATE || obj.type == TYPE_BOX) {value *= (1 + (rgen_from_obj(obj).rand() % 20));}
	else if (obj.type == TYPE_PAPER) {
		rand_gen_t rgen(rgen_from_obj(obj));
		if (rgen.rand_float() < 0.25) { // 25% of papers have some value
			float const val_mult((rgen.rand_float() < 0.25) ? 10 : 1); // 25% of papers have higher value
			value = val_mult*(2 + (rgen.rand()%10))*(1 + (rgen.rand()%10));
		}
	}
	else if (obj.type == TYPE_MONEY) {
		unsigned const num_bills(round_fp(obj.dz()/(0.01*obj.get_sz_dim(obj.dim))));
		value *= num_bills;
	}
	if (obj.flags & RO_FLAG_USED) {value = 0.01*floor(50.0*value);} // used objects have half value, rounded down to the nearest cent
	return value;
}
float get_obj_weight(room_object_t const &obj) {
	return get_taken_obj_type(obj).weight; // constant per object type, for now, but really should depend on object size/volume
}
bool is_consumable(room_object_t const &obj) {
	return (in_building_gameplay_mode() && obj.type == TYPE_BOTTLE && !obj.is_bottle_empty() && !(obj.flags & RO_FLAG_NO_CONS));
}

void show_weight_limit_message() {
	std::ostringstream oss;
	oss << "Over weight limit of " << global_building_params.player_weight_limit << " lbs";
	print_text_onscreen(oss.str(), RED, 1.0, 1.5*TICKS_PER_SECOND, 0);
}

class phone_manager_t {
	bool is_enabled=0, is_ringing = 0, is_on=0;
	double stop_ring_time=0.0, next_ring_time=0.0, next_cycle_time=0.0, auto_off_time=0.0, next_button_time=0.0;
	rand_gen_t rgen;

	void schedule_next_ring() {next_ring_time = tfticks + (double)rgen.rand_uniform(30.0, 120.0)*TICKS_PER_SECOND;} // 30s to 2min
public:
	bool is_phone_ringing() const {return is_ringing;}
	bool is_phone_on     () const {return is_on     ;}

	void next_frame() {
		if (!is_enabled || !camera_in_building) {} // do nothing
		else if (is_ringing) {
			if (tfticks > stop_ring_time) {is_ringing = 0; schedule_next_ring();} // stop automatically
			else if (tfticks > next_cycle_time) { // start a new ring cycle
				gen_sound_thread_safe_at_player(get_sound_id_for_file("phone_ring.wav"), 1.0);
				register_building_sound_at_player(1.0);
				next_cycle_time += 4.2*TICKS_PER_SECOND; // 4.2s between rings
			}
		}
		else {
			if (tfticks > next_ring_time) { // start a new ring cycle
				is_ringing      = 1;
				stop_ring_time  = tfticks + (double)rgen.rand_uniform(12.0, 24.0)*TICKS_PER_SECOND; // 10-20s into the future
				next_cycle_time = tfticks; // cycle begins now
			}
			if (is_on && tfticks > auto_off_time) {is_on = 0;} // auto off
		}
	}
	void enable() {
		if (is_enabled) return; // already enabled
		is_enabled = 1;
		is_ringing = is_on = 0;
		schedule_next_ring();
	}
	void disable() {
		is_enabled = is_ringing = is_on = 0;
	}
	void player_action() {
		if (tfticks < next_button_time) return; // skip if pressed immediately after the last press (switch debouncer)
		next_button_time = tfticks + 0.25*TICKS_PER_SECOND;

		if (is_ringing) { // switch off
			is_ringing     = is_on = 0;
			stop_ring_time = tfticks; // now
			schedule_next_ring();
		}
		else {
			is_on ^= 1;
			if (is_on) {auto_off_time = tfticks + 4.0*TICKS_PER_SECOND;} // 4s auto off delay
		}
		gen_sound_thread_safe_at_player(SOUND_CLICK, 0.5);
	}
};

phone_manager_t phone_manager;

struct tape_manager_t {
	bool in_use;
	float last_toggle_time;
	vector<point> points;
	point last_pos;
	room_object_t tape;
	building_t *cur_building;

	tape_manager_t() : in_use(0), last_toggle_time(0.0), cur_building(nullptr) {}

	void toggle_use(room_object_t const &tape_, building_t *building) {
		if ((tfticks - last_toggle_time) < 0.5*TICKS_PER_SECOND) return; // don't toggle too many times per frame
		last_toggle_time = tfticks;
		if (in_use) {clear();} // end use
		else {tape = tape_;	cur_building = building; in_use = 1;} // begin use
	}
	void clear() {
		if (cur_building != nullptr) {cur_building->maybe_update_tape(last_pos, 1);} // end_of_tape=1
		cur_building = nullptr;
		points.clear();
		tape = room_object_t();
		in_use = 0;
	}
};

tape_manager_t tape_manager;

class player_inventory_t { // manages player inventory, health, and other stats
	vector<carried_item_t> carried; // interactive items the player is currently carrying
	float cur_value, cur_weight, tot_value, tot_weight, damage_done, best_value, player_health, drunkenness, bladder, bladder_time, prev_player_zval;
	bool prev_in_building, has_key;

	void register_player_death(unsigned sound_id, std::string const &why) {
		point const xlate(get_camera_coord_space_xlate());
		place_player_at_xy(xlate.x, xlate.y); // move back to the origin/spawn location
		gen_sound_thread_safe_at_player(sound_id);
		print_text_onscreen(("You Have Died" + why), RED, 2.0, 2*TICKS_PER_SECOND, 10);
		clear(); // respawn
	}
public:
	player_inventory_t() {clear_all();}

	void clear() { // called on player death
		max_eq(best_value, tot_value);
		cur_value     = cur_weight = tot_value = tot_weight = damage_done = 0.0;
		drunkenness   = bladder = bladder_time = prev_player_zval = 0.0;
		player_health = 1.0; // full health
		prev_in_building = has_key = 0;
		phone_manager.disable();
		carried.clear();
		tape_manager.clear();
	}
	void clear_all() { // called on game mode init
		tot_value = best_value = 0.0;
		clear();
	}
	void take_damage(float amt) {player_health -= amt*(1.0f - 0.75f*min(drunkenness, 1.0f));} // up to 75% damage reduction when drunk
	void record_damage_done(float amt) {damage_done += amt;}
	void return_object_to_building(room_object_t const &obj) {damage_done -= get_obj_value(obj);}
	bool check_weight_limit(float weight) const {return ((cur_weight + weight) <= global_building_params.player_weight_limit);}
	bool can_pick_up_item(room_object_t const &obj) const {return check_weight_limit(get_obj_weight(obj));}
	float get_carry_weight_ratio() const {return min(1.0f, cur_weight/global_building_params.player_weight_limit);}
	float get_speed_mult () const {return (1.0f - 0.4f*get_carry_weight_ratio())*((bladder > 0.9) ? 0.6 : 1.0);} // 40% reduction for heavy load, 40% reduction for full bladder
	float get_drunkenness() const {return drunkenness;}
	bool  player_is_dead () const {return (player_health <= 0.0);}
	bool  player_has_key () const {return has_key;}

	bool can_open_door(door_t const &door) const {
		if (door.is_closed_and_locked() && !has_key) {
			print_text_onscreen("Door is locked", RED, 1.0, 2.0*TICKS_PER_SECOND, 0);
			gen_sound_thread_safe_at_player(SOUND_CLICK, 1.0, 0.6);
			return 0;
		}
		if (door.blocked) {
			print_text_onscreen("Door is blocked", RED, 1.0, 2.0*TICKS_PER_SECOND, 0);
			gen_sound_thread_safe_at_player(SOUND_DOOR_CLOSE, 1.0, 0.6);
			return 0;
		}
		return 1;
	}
	void switch_item(bool dir) { // Note: current item is always carried.back()
		if (carried.size() <= 1) return; // no other item to switch to
		if (dir) {std::rotate(carried.begin(), carried.begin()+1, carried.end());}
		else     {std::rotate(carried.begin(), carried.end  ()-1, carried.end());}
		gen_sound_thread_safe_at_player(SOUND_CLICK, 0.5);
		tape_manager.clear();
	}
	void add_item(room_object_t const &obj) {
		float health(0.0), drunk(0.0); // add these fields to bldg_obj_type_t?
		bool const bladder_was_full(bladder >= 0.9);
		float const value(get_obj_value(obj));
		damage_done += value;
		colorRGBA text_color(GREEN);
		std::ostringstream oss;
		oss << get_taken_obj_type(obj).name;

		if (is_consumable(obj)) { // nonempty bottle, consumable
			switch (obj.get_bottle_type()) {
			case 0: health =  0.25; break; // water
			case 1: health =  0.50; break; // Coke
			case 2: drunk  =  0.25; break; // beer
			case 3: drunk  =  0.50; break; // wine (entire bottle)
			case 4: health = -0.50; break; // poison - take damage
			default: assert(0);
			}
		}
		if (health > 0.0) { // heal
			player_health = min(1.0f, (player_health + health));
			oss << ": +" << round_fp(100.0*health) << "% Health";
		}
		if (health < 0.0) { // take damage
			player_health += health;
			oss << ": " << round_fp(100.0*health) << "% Health";
			text_color = RED;
			add_camera_filter(colorRGBA(RED, 0.25), 4, -1, CAM_FILT_DAMAGE); // 4 ticks of red damage
		}
		if (obj.type == TYPE_KEY) {
			has_key = 1; // mark as having the key, but it doesn't go into the inventory or contribute to weight or value
		}
		else if (health == 0.0 && drunk == 0.0) { // print value and weight if item is not consumed
			float const weight(get_obj_weight(obj));
			cur_value  += value;
			cur_weight += weight;
			
			if (obj.is_interactive()) {
				carried.push_back(obj);
				room_object_t &co(carried.back());

				if (obj.type == TYPE_BOOK) { // clear dim and dir for books
					float const dx(co.dx()), dy(co.dy()), dz(co.dz());

					if (dz > min(dx, dy)) { // upright book from a bookcase, put it on its side facing the player
						co.x2() = co.x1() + dz;
						co.y2() = co.y1() + max(dx, dy);
						co.z2() = co.z1() + min(dx, dy);
					}
					else if (co.dim) { // swap aspect ratio to make dim=0
						co.x2() = co.x1() + dy;
						co.y2() = co.y1() + dx;
					}
					co.dim = co.dir = 0;
					co.flags &= ~RO_FLAG_RAND_ROT; // remove the rotate bit
				}
				else if (obj.type == TYPE_PHONE) {
					if (co.dim) { // swap aspect ratio to make dim=0
						float const dx(co.dx()), dy(co.dy());
						co.x2() = co.x1() + dy;
						co.y2() = co.y1() + dx;
					}
					co.dim = co.dir = 0; // clear dim and dir
					phone_manager.enable();
				}
				tape_manager.clear();
			}
			oss << ": value $";
			if (value < 1.0 && value > 0.0) {oss << ((value < 0.1) ? "0.0" : "0.") << round_fp(100.0*value);} // make sure to print the leading/trailing zero for cents
			else {oss << value;}
			oss << " weight " << get_obj_weight(obj) << " lbs";
		}
		else { // add one drink to the bladder, 25% of capacity
			bladder = min(1.0f, (bladder + 0.25f));
		}
		if (drunk > 0.0) {
			drunkenness += drunk;
			oss << ": +" << round_fp(100.0*drunk) << "% Drunkenness";
			if (drunkenness > 0.99f && (drunkenness - drunk) <= 0.99f) {oss << "\nYou are drunk"; text_color = DK_GREEN;}
		}
		if (!bladder_was_full && bladder >= 0.9f) {oss << "\nYou need to use the bathroom"; text_color = YELLOW;}
		print_text_onscreen(oss.str(), text_color, 1.0, 3*TICKS_PER_SECOND, 0);
	}
	bool take_person(bool &person_has_key, float person_height) {
		if (drunkenness < 1.5) { // not drunk enough
			print_text_onscreen("Not drunk enough", RED, 1.0, 2.0*TICKS_PER_SECOND, 0);
			return 0;
		}
		float const value(1000), weight((person_height > 0.025) ? 180.0 : 80.0); // always worth $1000; use height to select man vs. girl
		if (!check_weight_limit(weight)) {show_weight_limit_message(); return 0;}
		has_key    |= person_has_key; person_has_key = 0; // steal their key
		cur_value  += value;
		cur_weight += weight;
		std::ostringstream oss;
		oss << "zombie: value $" << value << " weight " << weight << " lbs";
		print_text_onscreen(oss.str(), GREEN, 1.0, 4*TICKS_PER_SECOND, 0);
		return 1; // success
	}
	bool try_use_last_item(room_object_t &obj) {
		if (carried.empty()) return 0; // no interactive carried item
		obj = carried.back(); // deep copy
		if (!obj.has_dstate()) {return obj.can_use();} // not a droppable/throwable item(ball); should always return 1
		remove_last_item(); // drop the item - remove it from our inventory
		return 1;
	}
	bool update_last_item_use_count(int val) { // val can be positive or negative
		if (val == 0) return 1; // no change
		assert(!carried.empty());
		carried_item_t &obj(carried.back());
		obj.flags |= RO_FLAG_USED;
		unsigned const capacity(get_room_obj_type(obj).capacity);

		if (capacity > 0) {
			max_eq(val, -int(obj.use_count)); // can't go negative
			obj.use_count += val;
			if (obj.use_count >= capacity) {remove_last_item(); return 0;} // remove after too many uses
		}
		return 1;
	}
	bool mark_last_item_used() {
		return update_last_item_use_count(1);
	}
	void remove_last_item() {
		assert(!carried.empty());
		room_object_t const &obj(carried.back());
		cur_value  -= get_obj_value (obj);
		cur_weight -= get_obj_weight(obj);
		cur_value   = 0.01*round_fp(100.0*cur_value ); // round to nearest cent
		cur_weight  = 0.01*round_fp(100.0*cur_weight); // round to nearest 0.01 lb
		assert(cur_value > -0.01 && cur_weight > -0.01); // sanity check for math
		max_eq(cur_value,  0.0f); // handle FP rounding error
		max_eq(cur_weight, 0.0f);
		carried.pop_back(); // Note: invalidates obj
		tape_manager.clear();
	}
	void collect_items(bool keep_interactive) { // called when player exits a building
		if (!keep_interactive) {has_key = 0;} // key only good for current building
		phone_manager.disable(); // phones won't ring when taken out of their building, since the player can't switch to them anyway
		tape_manager.clear();
		if (carried.empty() && cur_weight == 0.0 && cur_value == 0.0) return; // nothing to add
		float keep_value(0.0), keep_weight(0.0);

		if (keep_interactive) { // carried items don't contribute to collected value and weight; their value and weight remain in the inventory after collection
			for (auto i = carried.begin(); i != carried.end(); ++i) {
				keep_value  += get_obj_value (*i);
				keep_weight += get_obj_weight(*i);
			}
			cur_value  -= keep_value;
			cur_weight -= keep_weight;
			cur_value   = 0.01*round_fp(100.0*cur_value); // round to nearest cent
		}
		else {
			carried.clear();
		}
		if (cur_weight > 0.0 || cur_value > 0.0) { // we have some change in value or weight to print
			std::ostringstream oss;
			oss << "Added value $" << cur_value << " Added weight " << cur_weight << " lbs\n";
			tot_value  += cur_value;
			tot_weight += cur_weight;
			tot_value   = 0.01*round_fp(100.0*tot_value); // round to nearest cent
			oss << "Total value $" << tot_value << " Total weight " << tot_weight << " lbs";
			print_text_onscreen(oss.str(), GREEN, 1.0, 4*TICKS_PER_SECOND, 0);
		}
		cur_value  = keep_value;
		cur_weight = keep_weight;
	}
	void show_stats() const {
		if (!carried.empty()) {
			player_held_object = carried.back(); // deep copy last pickup object if usable
			
			if (player_held_object.type == TYPE_PHONE) {
				if (phone_manager.is_phone_ringing()) {player_held_object.flags |= RO_FLAG_EMISSIVE;} // show ring screen
				else if (phone_manager.is_phone_on()) {player_held_object.flags |= RO_FLAG_OPEN    ;} // show lock screen
			}
		}
		if (display_framerate) { // controlled by framerate toggle
			float const aspect_ratio((float)window_width/(float)window_height);

			if (cur_weight > 0.0 || tot_weight > 0.0 || best_value > 0.0) { // don't show stats until the player has picked something up
				std::ostringstream oss;
				oss << "Cur $" << cur_value << " / " << cur_weight << " lbs  Total $" << tot_value << " / " << tot_weight
					<< " lbs  Best $" << best_value << "  Damage $" << damage_done;
				
				if (!carried.empty()) {
					unsigned const capacity(get_room_obj_type(player_held_object).capacity);
					oss << "  [" << get_taken_obj_type(player_held_object).name; // print the name of the throwable object
					if (capacity > 0) {oss << " " << (capacity - carried.back().use_count) << "/" << capacity;} // print use/capacity
					oss << "]";
				}
				draw_text(GREEN, -0.010*aspect_ratio, -0.011, -0.02, oss.str(), 0.8); // size=0.8
			}
			if (phone_manager.is_phone_ringing() && player_held_object.type == TYPE_PHONE) { // player is holding a ringing phone, give them a hint
				draw_text(LT_BLUE, -0.001*aspect_ratio, -0.009, -0.02, "[Press space to silence phone]");
			}
			if (in_building_gameplay_mode()) {
				float const lvl(min(cur_building_sound_level, 1.0f));
				unsigned const num_bars(round_fp(20.0*lvl));

				if (num_bars > 0) { // display sound meter
					colorRGBA const color(lvl, (1.0 - lvl), 0.0, 1.0); // green => yellow => orange => red
					draw_text(color, -0.005*aspect_ratio, -0.010, -0.02, std::string(num_bars, '#'));
				}
				if (player_is_hiding && !phone_manager.is_phone_ringing()) {draw_text(LT_BLUE, -0.001*aspect_ratio, -0.009, -0.02, "[Hiding]");}
			}
		}
		if (in_building_gameplay_mode()) {
			// Note: shields is used for drunkenness; values are scaled from 0-1 to 0-100; powerup values are for bladder fullness
			draw_health_bar(100.0*player_health, 100.0*drunkenness, bladder, YELLOW, get_carry_weight_ratio(), WHITE);
			if (has_key) {show_key_icon();}
		}
	}
	void next_frame() {
		show_stats();
		phone_manager.next_frame(); // even if not in gameplay mode?
		if (!in_building_gameplay_mode()) return;
		// handle player fall damage logic
		point const camera_pos(get_camera_pos());
		float const fall_damage_start(3.0*CAMERA_RADIUS); // should be a function of building floor spacing?
		float const player_zval(camera_pos.z), delta_z(prev_player_zval - player_zval);
		
		if (camera_in_building != prev_in_building) {prev_in_building = camera_in_building;}
		else if (prev_player_zval != 0.0 && delta_z > fall_damage_start && camera_in_building) {
			// only take fall damage when inside the building (no falling off the roof for now)
			player_health -= 1.0f*(delta_z - fall_damage_start)/fall_damage_start;
			if (player_is_dead()) {register_player_death(SOUND_SQUISH, " of a fall"); return;} // dead
			gen_sound_thread_safe_at_player(SOUND_SQUISH, 0.5);
			add_camera_filter(colorRGBA(RED, 0.25), 4, -1, CAM_FILT_DAMAGE); // 4 ticks of red damage
			register_building_sound_at_player(1.0);
		}
		prev_player_zval = player_zval;
		// handle death events
		if (player_is_dead() ) {register_player_death(SOUND_SCREAM3, ""); return;} // dead
		if (drunkenness > 2.0) {register_player_death(SOUND_DROWN,   " of alcohol poisoning"); return;}
		// update state for next frame
		drunkenness = max(0.0f, (drunkenness - 0.0001f*fticks)); // slowly decrease over time
		
		if (player_near_toilet) { // empty bladder
			if (bladder > 0.9) {gen_sound_thread_safe_at_player(SOUND_GASP);} // urinate
			if (bladder > 0.0) { // toilet flush
#pragma omp critical(gen_sound)
				gen_delayed_sound(1.0, SOUND_FLUSH, camera_pos); // delay by 1s
				register_building_sound_at_player(0.5);
			}
			bladder = 0.0;
		}
		else if (bladder > 0.9) {
			bladder_time += fticks;

			if (bladder_time > 5.0*TICKS_PER_SECOND) { // play the "I have to go" sound
				gen_sound_thread_safe_at_player(SOUND_HURT);
				bladder_time = 0.0;
			}
		}
		player_near_toilet = 0;
	}
}; // end player_inventory_t

player_inventory_t player_inventory;

float get_player_drunkenness() {return player_inventory.get_drunkenness();}
float get_player_building_speed_mult() {return player_inventory.get_speed_mult();}
bool player_can_open_door(door_t const &door) {return player_inventory.can_open_door(door);}

void register_building_sound_for_obj(room_object_t const &obj, point const &pos) {
	float const weight(get_obj_weight(obj)), volume((weight <= 1.0) ? 0.0 : min(1.0f, 0.01f*weight)); // heavier objects make more sound
	register_building_sound(pos, volume);
}

void register_broken_object(room_object_t const &obj) {player_inventory.record_damage_done(get_obj_value(obj));}

bool register_player_object_pickup(room_object_t const &obj, point const &at_pos) {
	bool const can_pick_up(player_inventory.can_pick_up_item(obj));

	if (!do_room_obj_pickup) { // player has not used the pickup key, but we can still use this to notify the player that an object can be picked up
		can_pickup_bldg_obj = (can_pick_up ? 1 : 2);
		return 0;
	}
	if (!can_pick_up) {
		show_weight_limit_message();
		return 0;
	}
	if (is_consumable(obj)) {gen_sound_thread_safe_at_player(SOUND_GULP, 1.0 );}
	else                    {gen_sound_thread_safe_at_player(SOUND_ITEM, 0.25);}
	register_building_sound_for_obj(obj, at_pos);
	do_room_obj_pickup = 0; // no more object pickups
	return 1;
}

bool building_t::player_pickup_object(point const &at_pos, vector3d const &in_dir) {
	if (!has_room_geom()) return 0;
	return interior->room_geom->player_pickup_object(*this, at_pos, in_dir);
}
bool building_room_geom_t::player_pickup_object(building_t &building, point const &at_pos, vector3d const &in_dir) {
	point at_pos_rot(at_pos);
	vector3d in_dir_rot(in_dir);
	building.maybe_inv_rotate_pos_dir(at_pos_rot, in_dir_rot);
	float const range(3.0*CAMERA_RADIUS), drawer_range_max(2.5*CAMERA_RADIUS);
	float drawer_range(drawer_range_max), obj_dist(0.0);
	int const obj_id(find_nearest_pickup_object(building, at_pos_rot, in_dir_rot, range, obj_dist));
	if (obj_id >= 0) {min_eq(drawer_range, obj_dist);} // only include drawers that are closer than the pickup object
	if (open_nearest_drawer(building, at_pos_rot, in_dir_rot, drawer_range_max, 1, 0)) return 1; // try objects in drawers; pickup_item=1
	if (obj_id < 0) return 0; // no object to pick up
	room_object_t &obj(get_room_object_by_index(obj_id));

	if (obj.type == TYPE_SHELVES || (obj.type == TYPE_WINE_RACK && !(obj.flags & RO_FLAG_EXPANDED))) { // shelves or unexpanded wine rack
		assert(!(obj.flags & RO_FLAG_EXPANDED)); // should not have been expanded
		expand_object(obj);
		bool const picked_up(player_pickup_object(building, at_pos, in_dir)); // call recursively on contents
		// if we picked up an object, assume the VBOs have already been updated; otherwise we need to update them to expand this object
		if (!picked_up) {create_small_static_vbos(building);} // assumes expanded objects are all "small"
		return picked_up;
	}
	if (obj.type == TYPE_BCASE) {
		static vector<room_object_t> books;
		books.clear();
		get_bookcase_books(obj, books);
		int closest_obj_id(-1);
		float dmin_sq(0.0);
		point const p2(at_pos + in_dir*range);

		for (auto i = books.begin(); i != books.end(); ++i) {
			point p1c(at_pos), p2c(p2);
			if (!do_line_clip(p1c, p2c, i->d))  continue; // test ray intersection vs. bcube
			float const dsq(p2p_dist(at_pos, p1c)); // use closest intersection point
			if (dmin_sq > 0.0 && dsq > dmin_sq) continue; // not the closest
			closest_obj_id = (i - books.begin()); // valid pickup object
			dmin_sq = dsq; // this object is the closest
		} // for i
		if (dmin_sq == 0.0) return 0; // no book to pick up
		room_object_t &book(books[closest_obj_id]);
		if (!register_player_object_pickup(book, at_pos)) return 0;
		obj.set_combined_flags(obj.get_combined_flags() | (1<<(book.item_flags&31))); // set flag bit to remove this book from the bookcase
		player_inventory.add_item(book);
		update_draw_state_for_room_object(book, building, 1);
		return 1;
	}
	if (!register_player_object_pickup(obj, at_pos)) return 0;
	remove_object(obj_id, building);
	return 1;
}

void building_t::register_player_enter_building() const {
	// nothing to do yet
}
void building_t::register_player_exit_building() const {
	// only collect items in gameplay mode where there's a risk the player can lose them; otherwise, let the player carry items between buildings
	player_inventory.collect_items(!in_building_gameplay_mode());
}

bool is_obj_in_or_on_obj(room_object_t const &parent, room_object_t const &child) {
	if (parent.type == TYPE_WINE_RACK && parent.contains_pt(child.get_cube_center()))     return 1; // check for wine bottles left in wine rack
	if (fabs(child.z1() - parent.z2()) < 0.05*parent.dz() && child.intersects_xy(parent)) return 1; // zval test
	if (parent.type == TYPE_BOX && parent.is_open() && parent.contains_cube(child))       return 1; // open box with an object inside
	return 0;
}
bool object_has_something_on_it(room_object_t const &obj, vector<room_object_t> const &objs, vector<room_object_t>::const_iterator objs_end) {
	// only these types can have objects placed on them (what about TYPE_SHELF?)
	if (obj.type != TYPE_TABLE && obj.type != TYPE_DESK && obj.type != TYPE_COUNTER && obj.type != TYPE_DRESSER && obj.type != TYPE_NIGHTSTAND &&
		obj.type != TYPE_BOX && obj.type != TYPE_CRATE && obj.type != TYPE_WINE_RACK && obj.type != TYPE_BOOK) return 0;

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type == TYPE_BLOCKER)      continue; // ignore blockers (from removed objects)
		if (*i == obj)                    continue; // skip self (bcube check)
		if (is_obj_in_or_on_obj(obj, *i)) return 1;
	}
	return 0;
}

cube_t get_true_obj_bcube(room_object_t const &obj) {
	if (obj.type == TYPE_PEN || obj.type == TYPE_PENCIL) {
		cube_t obj_bcube(obj);
		obj_bcube.expand_in_dim(!obj.dim, obj.get_sz_dim(!obj.dim));
		return obj_bcube;
	}
	if (obj.type == TYPE_BED) { // do more accurate check with various parts of the bed
		cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
		get_bed_cubes(obj, cubes);
		return cubes[3]; // check mattress only, since we can only take the mattress, sheets, and pillows
	}
	if (obj.is_obj_model_type()) {
		cube_t obj_bcube(obj);
		obj_bcube.expand_by(-0.1*obj.get_size()); // since models don't fill their bcubes, shrink them a bit when doing a ray query
		return obj_bcube;
	}
	return obj; // unmodified
}

bool obj_has_open_drawers(room_object_t const &obj) {return ((obj.type == TYPE_NIGHTSTAND || obj.type == TYPE_DRESSER || obj.type == TYPE_DESK) && obj.drawer_flags);}

int building_room_geom_t::find_nearest_pickup_object(building_t const &building, point const &at_pos, vector3d const &in_dir, float range, float &obj_dist) const {
	int closest_obj_id(-1);
	float dmin_sq(0.0);
	point const p2(at_pos + in_dir*range);

	for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
		auto const &obj_vect((vect_id == 1) ? expanded_objs : objs), &other_obj_vect((vect_id == 1) ? objs : expanded_objs);
		unsigned const obj_id_offset((vect_id == 1) ? objs.size() : 0);
		auto objs_end((vect_id == 1) ? expanded_objs.end() : get_stairs_start()); // skip stairs and elevators
		auto other_objs_end((vect_id == 1) ? get_stairs_start() : expanded_objs.end());

		for (auto i = obj_vect.begin(); i != objs_end; ++i) {
			if (!get_room_obj_type(*i).pickup) continue; // this object type can't be picked up
			cube_t const obj_bcube(get_true_obj_bcube(*i));
			point p1c(at_pos), p2c(p2);
			if (!do_line_clip(p1c, p2c, obj_bcube.d)) continue; // test ray intersection vs. bcube
			float const dsq(p2p_dist(at_pos, p1c)); // use closest intersection point
			if (dmin_sq > 0.0 && dsq > dmin_sq)       continue; // not the closest
		
			if (i->type == TYPE_CLOSET || (i->type == TYPE_STALL && i->shape != SHAPE_SHORT)) { // can only take short stalls (separating urinals)
				if (!i->is_open() && !i->contains_pt(at_pos)) { // stalls/closets block the player from taking toilets/boxes unless open, or the player is inside
					closest_obj_id = -1;
					dmin_sq = dsq;
				}
				continue;
			}
			if (i->type == TYPE_CHAIR) { // separate back vs. seat vs. legs check for improved accuracy
				cube_t cubes[3]; // seat, back, legs_bcube
				get_chair_cubes(*i, cubes);
				bool intersects(0);
				for (unsigned n = 0; n < 3 && !intersects; ++n) {intersects |= cubes[n].line_intersects(p1c, p2c);}
				if (!intersects) continue;
			}
			if (i->type == TYPE_HANGER_ROD && i->item_flags > 0) { // nonempty hanger rod
				// search for hangers and don't allow hanger rod to be taken until the hangers are all taken
				bool has_hanger(0);

				for (auto j = i+1; j != (i + i->item_flags); ++j) { // iterate over all objects hanging on the hanger rod and look for untaken hangers
					if (j->type == TYPE_HANGER) {has_hanger = 1; break;}
				}
				if (has_hanger) continue;
			}
			if (i->type == TYPE_HANGER  && (i->flags & RO_FLAG_HANGING) && (i+1) != objs_end && (i+1)->type == TYPE_CLOTHES) continue; // hanger with clothes - must take clothes first
			if (i->type == TYPE_MIRROR  && !i->is_house())                continue; // can only pick up mirrors from houses, not office buildings
			if (i->type == TYPE_TABLE   && i->shape == SHAPE_CUBE)        continue; // can only pick up short (TV) tables and cylindrical tables
			if (i->type == TYPE_BED     && (i->flags & RO_FLAG_TAKEN3))   continue; // can only take pillow, sheets, and mattress - not the frame
			if (i->type == TYPE_SHELVES && (i->flags & RO_FLAG_EXPANDED)) continue; // shelves are already expanded, can no longer select this object
			if (obj_has_open_drawers(*i))                                 continue; // can't take if any drawers are open
			if (object_has_something_on_it(*i,       obj_vect, objs_end)) continue; // can't remove a table, etc. that has something on it
			if (object_has_something_on_it(*i, other_obj_vect, other_objs_end)) continue; // check the other one as well
			if (building.check_for_wall_ceil_floor_int(at_pos, p1c))      continue; // skip if it's on the other side of a wall, ceiling, or floor
			closest_obj_id = (i - obj_vect.begin()) + obj_id_offset; // valid pickup object
			dmin_sq = dsq; // this object is the closest, even if it can't be picked up
		} // for i
	} // for vect_id
	obj_dist = sqrt(dmin_sq);
	return closest_obj_id;
}

bool is_counter(room_object_t const &obj) {return (obj.type == TYPE_COUNTER || obj.type == TYPE_KSINK);}

bool building_room_geom_t::open_nearest_drawer(building_t &building, point const &at_pos, vector3d const &in_dir, float range, bool pickup_item, bool check_only) {
	int closest_obj_id(-1);
	float dmin_sq(0.0);
	point const p2(at_pos + in_dir*range);
	auto objs_end(get_std_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		bool const is_counter_type(is_counter(*i) || i->type == TYPE_CABINET);
		if (!(i->type == TYPE_DRESSER || i->type == TYPE_NIGHTSTAND || i->type == TYPE_DESK || // drawers that can be opened or items picked up from
			(!pickup_item && is_counter_type))) continue; // doors that can be opened (no item pickup)
		cube_t bcube(*i);
		float &front_face(bcube.d[i->dim][i->dir]);
		if ((front_face < at_pos[i->dim]) ^ i->dir) continue; // can't open from behind
		if (!is_counter_type) {front_face += 0.75*(i->dir ? 1.0 : -1.0)*i->get_sz_dim(i->dim);} // expand outward to include open drawers
		point p1c(at_pos), p2c(p2);
		if (!do_line_clip(p1c, p2c, bcube.d)) continue; // test ray intersection vs. bcube
		float const dsq(p2p_dist_sq(at_pos, p1c)); // use closest intersection point
		if (dmin_sq > 0.0 && dsq > dmin_sq) continue; // not the closest
		if (building.check_for_wall_ceil_floor_int(at_pos, p1c)) continue; // skip if it's on the other side of a wall, ceiling, or floor
		closest_obj_id = (i - objs.begin());
		dmin_sq = dsq; // this object is the closest, even if it can't be picked up
	} // for i
	if (closest_obj_id < 0) return 0; // no object
	room_object_t &obj(objs[closest_obj_id]);
	room_object_t drawers_part;
	vect_cube_t drawers; // or doors
	bool const has_doors(is_counter(obj) || obj.type == TYPE_CABINET);
	float drawer_extend(0.0); // signed, for drawers only

	// Note: this is a messy solution and must match the drawing code, but it's unclear how else we can get the location of the drawers
	if (has_doors) {get_cabinet_or_counter_doors(obj, drawers);}
	else {
		if (obj.type == TYPE_DESK) {
			if (!(obj.room_id & 3)) return 0; // no drawers for this desk
			drawers_part = get_desk_drawers_part(obj);
			bool const side(obj.obj_id & 1);
			drawers_part.d[!obj.dim][side] -= (side ? 1.0 : -1.0)*0.85*get_tc_leg_width(obj, 0.06);
		}
		else {
			drawers_part = get_dresser_middle(obj);
			drawers_part.expand_in_dim(!obj.dim, -0.5*get_tc_leg_width(obj, 0.10));
		}
		drawer_extend = get_drawer_cubes(drawers_part, drawers, 0, 1); // front_only=0, inside_only=1
	}
	dmin_sq        = 0.0;
	closest_obj_id = -1;

	for (auto i = drawers.begin(); i != drawers.end(); ++i) {
		point p1c(at_pos), p2c(p2);
		if (!do_line_clip(p1c, p2c, i->d)) continue; // test ray intersection vs. drawer
		float const dsq(p2p_dist_sq(at_pos, p1c)); // use closest intersection point
		if (dmin_sq == 0.0 || dsq < dmin_sq) {closest_obj_id = (i - drawers.begin()); dmin_sq = dsq;} // update if closest
	}
	if (closest_obj_id < 0) return 0; // no drawer
	cube_t const &drawer(drawers[closest_obj_id]); // Note: drawer cube is the interior part
	
	if (pickup_item && !has_doors) { // pick up item in drawer rather than opening drawer; no pickup items behind doors yet
		if (!(obj.drawer_flags & (1U << closest_obj_id))) return 0; // drawer is not open
		room_object_t const item(get_item_in_drawer(obj, drawer, closest_obj_id));
		if (item.type == TYPE_NONE) return 0; // no item
		if (check_only) return 1;
		if (!register_player_object_pickup(item, at_pos)) return 0;
		obj.item_flags |= (1U << closest_obj_id); // flag item as taken
		player_inventory.add_item(item);
		update_draw_state_for_room_object(item, building, 1);
	}
	else { // open or close the drawer/door
		cube_t c_test(drawer);
		if (has_doors) {c_test.d[obj.dim][obj.dir] += (obj.dir ? 1.0 : -1.0)*drawer.get_sz_dim(!obj.dim);} // expand outward by the width of the door
		else {c_test.d[obj.dim][obj.dir] += drawer_extend;} // drawer
		if (cube_intersects_moved_obj(c_test, closest_obj_id)) return 0; // blocked, can't open; ignore this object
		if (check_only) return 1;
		obj.drawer_flags ^= (1U << (unsigned)closest_obj_id); // toggle flag bit for selected drawer
		point const drawer_center(drawer.get_cube_center());
		if (has_doors) {building.play_door_open_close_sound(drawer_center, obj.is_open(), 0.5, 1.5);}
		else {gen_sound_thread_safe(SOUND_SLIDING, building.local_to_camera_space(drawer_center), 0.5);}
		register_building_sound(drawer_center, 0.4);
		
		if (has_doors) { // expand any items in the cabinet so that the player can pick them up
			// find any cabinets adjacent to this one in the other dim (inside corner) and ensure the opposing door is closed so that they don't intersect
			for (auto i = objs.begin(); i != objs_end; ++i) {
				if ((is_counter(obj) && !is_counter(*i)) || (obj.type == TYPE_CABINET && i->type != TYPE_CABINET)) continue; // wrong object type
				if (i->dim == obj.dim) continue; // not opposing dim (also skips obj itself)
				float const dir_sign(i->dir ? 1.0 : -1.0);
				cube_t i_exp(*i);
				i_exp.d[i->dim][i->dir] += dir_sign*i->get_sz_dim(i->dim); // expand other counter/cabinet to account for open doors
				if (!i_exp.intersects(c_test)) continue;
				get_cabinet_or_counter_doors(*i, drawers);

				for (auto j = drawers.begin(); j != drawers.end(); ++j) {
					cube_t drawer_exp(*j);
					drawer_exp.d[i->dim][i->dir] += dir_sign*j->get_sz_dim(!i->dim); // expand outward by the width of the door
					if (drawer_exp.intersects(c_test)) {i->drawer_flags &= ~(1U << (j - drawers.begin()));} // make sure any intersecting doors are closed
				}
			} // for i
			// Note: expanding cabinets by opening a single door will allow the player to take items from anywhere in the cabinet, even if behind a closed door
			expand_object(obj);
			update_draw_state_for_room_object(obj, building, 0);
		}
		else {
			create_small_static_vbos(building); // only need to update small objects for drawers
		}
	}
	return 1;
}

room_object_t &building_room_geom_t::get_room_object_by_index(unsigned obj_id) {
	if (obj_id < objs.size()) {return objs[obj_id];}
	unsigned const exp_obj_id(obj_id - objs.size());
	assert(exp_obj_id < expanded_objs.size());
	return expanded_objs[exp_obj_id];
}

void building_room_geom_t::remove_object(unsigned obj_id, building_t &building) {
	room_object_t &obj(get_room_object_by_index(obj_id));
	room_object_t const old_obj(obj); // deep copy
	assert(obj.type != TYPE_ELEVATOR); // elevators require special updates for drawing logic and cannot be removed at this time
	player_inventory.add_item(obj);
	bldg_obj_type_t const type(get_taken_obj_type(obj)); // capture type before updating obj
	bool const is_light(obj.type == TYPE_LIGHT);

	if (obj.type == TYPE_PICTURE && !(obj.flags & RO_FLAG_TAKEN1)) {obj.flags |= RO_FLAG_TAKEN1;} // take picture, leave frame
	else if (obj.type == TYPE_TPROLL && !(obj.flags & (RO_FLAG_TAKEN1 | RO_FLAG_WAS_EXP))) {obj.flags |= RO_FLAG_TAKEN1;} // take toilet paper roll, leave holder; not for expanded TP rolls
	else if (obj.type == TYPE_BED) {
		if      (obj.flags & RO_FLAG_TAKEN2) {obj.flags |= RO_FLAG_TAKEN3;} // take mattress
		else if (obj.flags & RO_FLAG_TAKEN1) {obj.flags |= RO_FLAG_TAKEN2;} // take sheets
		else {obj.flags |= RO_FLAG_TAKEN1;} // take pillow(s)
	}
	else if (obj.type == TYPE_PLANT && !(obj.flags & RO_FLAG_ADJ_BOT)) { // plant not on a table/desk
		if      (obj.flags & RO_FLAG_TAKEN2) {obj.type = TYPE_BLOCKER;} // take pot - gone
		else if (obj.flags & RO_FLAG_TAKEN1) {obj.flags |= RO_FLAG_TAKEN2;} // take dirt
		else {obj.flags |= RO_FLAG_TAKEN1;} // take plant
	}
	else if (obj.type == TYPE_TOILET || obj.type == TYPE_SINK) { // leave a drain in the floor
		cube_t drain;
		drain.set_from_point(point(obj.xc(), obj.yc(), obj.z1()));
		drain.expand_by_xy(0.065*obj.dz());
		drain.z2() += 0.02*obj.dz();
		obj = room_object_t(drain, TYPE_DRAIN, obj.room_id, 0, 0, RO_FLAG_NOCOLL, obj.light_amt, SHAPE_CYLIN, DK_GRAY);
		create_small_static_vbos(building);
	}
	else { // replace it with an invisible blocker that won't collide with anything
		obj.type  = TYPE_BLOCKER;
		obj.flags = (RO_FLAG_NOCOLL | RO_FLAG_INVIS);
	}
	if (is_light) {clear_and_recreate_lights();}
	update_draw_state_for_room_object(old_obj, building, 1);
}

void building_room_geom_t::update_draw_state_for_obj_type_flags(bldg_obj_type_t const &type, building_t &building) {
	if (type.lg_sm & 2) {create_small_static_vbos(building);} // small object
	if (type.lg_sm & 1) {create_static_vbos      (building);} // large object
	if (type.is_model ) {create_obj_model_insts  (building);} // 3D model
	//if (type.ai_coll  ) {building.invalidate_nav_graph();} // removing this object should not affect the AI navigation graph
}
// Note: called when adding, removing, or moving objects
void building_room_geom_t::update_draw_state_for_room_object(room_object_t const &obj, building_t &building, bool was_taken) {
	// reuild necessary VBOs and other data structures
	if (obj.is_dynamic()) {mats_dynamic.clear();} // dynamic object
	else if (obj.type == TYPE_BUTTON && (obj.flags & RO_FLAG_IN_ELEV)) {update_dynamic_draw_data();} // interior elevator buttons are drawn as dynamic objects
	else {update_draw_state_for_obj_type_flags((was_taken ? get_taken_obj_type(obj) : get_room_obj_type(obj)), building);} // static object
	modified_by_player = 1; // flag so that we avoid re-generating room geom if the player leaves and comes back
}

int building_room_geom_t::find_avail_obj_slot() const {
	auto objs_end(get_stairs_start()); // skip stairs and elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type == TYPE_BLOCKER) {return int(i - objs.begin());} // blockers are used as temporaries for room object placement and to replace removed objects
	}
	return -1; // no slot found
}
void building_room_geom_t::add_expanded_object(room_object_t const &obj) {
	for (auto i = expanded_objs.begin(); i != expanded_objs.end(); ++i) {
		if (i->type == TYPE_BLOCKER) {*i = obj; return;} // found a slot - done
	}
	expanded_objs.push_back(obj); // not found - in this case we can add a new object
}
bool building_room_geom_t::add_room_object(room_object_t const &obj, building_t &building, bool set_obj_id, vector3d const &velocity) {
	assert(obj.type != TYPE_LIGHT && get_room_obj_type(obj).pickup); // currently must be a pickup object, and not a light

	if (!set_obj_id && (obj.flags & RO_FLAG_WAS_EXP)) { // if object was expanded, and it's not a dynamic object, use an expanded slot (books, etc.)
		assert(velocity == zero_vector);
		add_expanded_object(obj);
	}
	else {
		int const obj_id(find_avail_obj_slot());
		if (obj_id < 0) return 0; // no slot found
		room_object_t &added_obj(get_room_object_by_index(obj_id));
		added_obj = obj; // overwrite with new object
		if (set_obj_id) {added_obj.obj_id = (uint16_t)(obj.has_dstate() ? allocate_dynamic_state() : obj_id);}
		if (velocity != zero_vector) {get_dstate(added_obj).velocity = velocity;}
	}
	update_draw_state_for_room_object(obj, building, 0);
	return 1;
}

bool is_movable(room_object_t const &obj) {
	if (obj.no_coll() || obj.type == TYPE_BLOCKER) return 0; // no blockers
	bldg_obj_type_t const &bot(get_room_obj_type(obj));
	return (bot.weight >= 40.0 && !bot.attached); // heavy non-attached objects, including tables
}
bool building_t::move_nearest_object(point const &at_pos, vector3d const &in_dir, float range, int mode) { // mode: 0=normal, 1=pull
	assert(has_room_geom());
	int closest_obj_id(-1);
	float dmin_sq(0.0);
	point const p2(at_pos + in_dir*range);
	vector<room_object_t> &objs(interior->room_geom->objs), &expanded_objs(interior->room_geom->expanded_objs);
	auto objs_end(interior->room_geom->get_stairs_start()); // skip stairs and elevators

	// determine which object the player may be choosing to move
	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->no_coll() || i->type == TYPE_BLOCKER) continue; // not interactive
		cube_t const obj_bcube(get_true_obj_bcube(*i));
		point p1c(at_pos), p2c(p2);
		if (!do_line_clip(p1c, p2c, obj_bcube.d)) continue; // test ray intersection vs. bcube
		float const dsq(p2p_dist(at_pos, p1c)); // use closest intersection point
		if (dmin_sq > 0.0 && dsq > dmin_sq)       continue; // not the closest
		if (obj_has_open_drawers(*i))             continue; // can't move if any drawers are open
		if (check_for_wall_ceil_floor_int(at_pos, p1c)) continue; // skip if it's on the other side of a wall, ceiling, or floor
		closest_obj_id = (i - objs.begin()); // valid pickup object
		dmin_sq = dsq; // this object is the closest
	} // for i
	if (closest_obj_id < 0) return 0;

	// determine move direction and distance
	room_object_t &obj(objs[closest_obj_id]);
	if (!is_movable(obj))        return 0; // closest object isn't movable
	if (obj.contains_pt(at_pos)) return 0; // player is inside this object?
	float const move_dist(rand_uniform(0.5, 1.0)*CAMERA_RADIUS*(100.0f/max(75.0f, get_room_obj_type(obj).weight))); // heavier objects move less; add some global randomness
	vector3d delta(obj.closest_pt(at_pos) - at_pos);
	delta.z = 0.0; // XY only
	delta.normalize();
	if (mode == 1) {delta.negate();} // changes push to pull ('r' key vs 'e' key)
	cube_t player_bcube;
	player_bcube.set_from_sphere(at_pos, get_scaled_player_radius());
	player_bcube.z1() -= camera_zh;

	// attempt to move the object
	for (unsigned mdir = 0; mdir < 3; ++mdir) { // X+Y, closer dim, further dim
		vector3d move_vector(zero_vector);
		if (mdir == 0) {move_vector = delta*move_dist;} // move diag in XY
		else { // move in one dim
			if (delta.x == 0.0 || delta.y == 0.0) break; // no more dims to try (only one mdir iteration)
			bool const dim(fabs(delta.x) < fabs(delta.y));
			move_vector[dim] = delta[dim]*move_dist;
		}
		for (unsigned n = 0; n < 5; ++n, move_vector *= 0.5) { // move in several incrementally smaller steps
			room_object_t moved_obj(obj);
			moved_obj += move_vector; // only the position changes
			if (player_bcube.intersects(moved_obj)) continue; // don't intersect the player - applies to pull mode
			if (!is_obj_pos_valid(moved_obj, 1))    continue; // try a smaller movement; keep_in_room=1
			bool bad_placement(0);

			for (auto i = objs.begin(); i != objs_end && !bad_placement; ++i) {
				if (i == objs.begin() + closest_obj_id)      continue; // skip self
				if (i->no_coll() || i->type == TYPE_BLOCKER) continue; // skip non-colliding objects and blockers that add clearance between objects as these won't block this object
				
				if (i->type == TYPE_CLOSET && i->is_open() && i->is_small_closet()) { // check open closet door collision
					cube_t cubes[5];
					get_closet_cubes(*i, cubes, 1); // get cubes for walls and door; for_collision=1
					for (unsigned n = 0; n < 5; ++n) {bad_placement |= cubes[n].intersects(moved_obj);}
				}
				else {bad_placement = i->intersects(moved_obj);}
			} // for i
			// Note: okay to skip expanded_objs because these should already be on/inside some other object; this allows us to move wine racks containing wine
			if (bad_placement) continue; // intersects another object, try a smaller movement
			bldg_obj_type_t type_flags(get_room_obj_type(obj));

			// move objects inside or on top of this one
			for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
				auto &obj_vect((vect_id == 1) ? expanded_objs : objs);
				auto obj_vect_end((vect_id == 1) ? expanded_objs.end() : objs_end); // skip stairs and elevators

				for (auto i = obj_vect.begin(); i != obj_vect_end; ++i) {
					if (i->type == TYPE_BLOCKER || *i == obj) continue; // ignore blockers and self
					if (!is_obj_in_or_on_obj(obj, *i))        continue;
					*i += move_vector; // move this object as well
					if (i->is_dynamic()) {interior->room_geom->mats_dynamic.clear();} // dynamic object
					else {
						bldg_obj_type_t const &type(get_room_obj_type(*i));
						type_flags.lg_sm    |= type.lg_sm;
						type_flags.is_model |= type.is_model;
						type_flags.ai_coll  |= type.ai_coll;
					}
				} // for i
			} // for vect_id
			// mark doors as blocked
			room_t const &room(get_room(obj.room_id));

			for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) { // check for door intersection
				if (i->open || !door_opens_inward(*i, room)) continue; // if the door is already open, or opens in the other direction, it can't be blocked
				bool const inc_open(0), check_dirs(i->get_check_dirs());
				if (is_cube_close_to_door(moved_obj, 0.0, inc_open, *i, check_dirs))              {i->blocked = 1; interior->door_state_updated = 1;} // newly blocked
				else if (i->blocked && is_cube_close_to_door(obj, 0.0, inc_open, *i, check_dirs)) {i->blocked = 0; interior->door_state_updated = 1;} // newly unblocked
			}
			// update this object
			obj = moved_obj; // keep this placement
			if (!obj.was_moved()) {interior->room_geom->moved_obj_ids.push_back(closest_obj_id);} // add to moved_obj_ids on first movement
			obj.flags |= RO_FLAG_MOVED;
			interior->room_geom->update_draw_state_for_obj_type_flags(type_flags, *this);
			interior->room_geom->modified_by_player = 1; // flag so that we avoid re-generating room geom if the player leaves and comes back
			gen_sound_thread_safe_at_player(SOUND_SLIDING);
			register_building_sound_at_player(0.7);
			return 1; // success
		} // for n
	} // for mdir
	return 0; // failed
}

void play_obj_fall_sound(room_object_t const &obj, point const &player_pos) {
	gen_sound_thread_safe(SOUND_OBJ_FALL, (get_camera_pos() + (obj.get_cube_center() - player_pos)));
	register_building_sound_for_obj(obj, player_pos);
}

bool building_t::maybe_use_last_pickup_room_object(point const &player_pos) {
	if (player_in_elevator) return 0; // can't use items in elevators
	assert(has_room_geom());
	static bool delay_use(0);
	static double last_use_time(0.0);
	if (delay_use && (tfticks - last_use_time) < 0.5*TICKS_PER_SECOND) return 0; // half second delay on prev item use or switch
	delay_use = 0;
	room_object_t obj;
	if (!player_inventory.try_use_last_item(obj)) return 0;

	if (obj.has_dstate()) { // it's a dynamic object (ball), throw it; only activated with use_object/'E' key
		float const cradius(get_scaled_player_radius());
		point dest(player_pos + (1.2f*(cradius + obj.get_radius()))*cview_dir);
		dest.z -= 0.5*cradius; // slightly below the player's face
		obj.translate(dest - point(obj.xc(), obj.yc(), obj.z1()));
		obj.flags |= RO_FLAG_DYNAMIC; // make it dynamic, assuming it will be dropped/thrown
		if (!interior->room_geom->add_room_object(obj, *this, 1, THROW_VELOCITY*cview_dir)) return 0;
		player_inventory.return_object_to_building(obj); // re-add this object's value
		play_obj_fall_sound(obj, player_pos);
		delay_use = 1;
	}
	else if (obj.can_use()) { // active with either use_object or fire key
		if (obj.type == TYPE_TPROLL) {
			point const dest(player_pos + (1.5f*get_scaled_player_radius())*cview_dir);
			if (!apply_toilet_paper(dest, cview_dir, 0.5*obj.dz())) return 0;
			player_inventory.mark_last_item_used();
		}
		else if (obj.type == TYPE_SPRAYCAN || obj.type == TYPE_MARKER) { // spraypaint or marker
			if (!apply_paint(player_pos, cview_dir, obj.color, obj.type)) return 0;
			player_inventory.mark_last_item_used();
		}
		else if (obj.type == TYPE_TAPE) {
			tape_manager.toggle_use(obj, this);
		}
		else if (obj.type == TYPE_BOOK) {
			float const half_width(0.5*max(max(obj.dx(), obj.dy()), obj.dz()));
			point dest(player_pos + (1.2f*(get_scaled_player_radius() + half_width))*cview_dir);
			if (!get_zval_for_obj_placement(dest, half_width, dest.z, 0)) return 0; // no suitable placement found; add_z_bias=0
			// orient based on the player's primary direction
			bool const place_dim(fabs(cview_dir.y) < fabs(cview_dir.x));

			if (obj.dim != place_dim) {
				float const dx(obj.dx()), dy(obj.dy());
				obj.x2() = obj.x1() + dy;
				obj.y2() = obj.y1() + dx;
			}
			obj.dim    = place_dim;
			obj.dir    = ((cview_dir[!place_dim] > 0) ^ place_dim);
			obj.flags |= (RO_FLAG_TAKEN1 | RO_FLAG_WAS_EXP);
			obj.translate(dest - point(obj.xc(), obj.yc(), obj.z1()));
			if (!interior->room_geom->add_room_object(obj, *this)) return 0;
			player_inventory.return_object_to_building(obj); // re-add this object's value
			player_inventory.remove_last_item(); // used
			play_obj_fall_sound(obj, player_pos);
			delay_use = 1;
		}
		else if (obj.type == TYPE_PHONE) {phone_manager.player_action();}
		else {assert(0);}
	}
	else {assert(0);}
	last_use_time = tfticks;
	return 1;
}

// adds two back-to-back quads for two sided lighting
void add_tape_quad(point const &p1, point const &p2, float width, color_wrapper const &color, quad_batch_draw &qbd, vector3d const &wdir=plus_z) {
	vector3d const dir(p2 - p1), wvect(0.5*width*wdir);
	vector3d normal(cross_product(dir, wdir).get_norm());
	point pts[4] = {(p1 - wvect), (p1 + wvect), (p2 + wvect), (p2 - wvect)};
	qbd.add_quad_pts(pts, color,  normal);
	swap(pts[1], pts[3]); // swap winding order and draw with reversed normal for two sided lighting
	qbd.add_quad_pts(pts, color, -normal);
}

void building_t::play_tape_sound(point const &sound_pos, float sound_gain) const {
	gen_sound_thread_safe(get_sound_id_for_file("tape.wav"), local_to_camera_space(sound_pos), sound_gain);
	register_building_sound(sound_pos, 0.35*sound_gain);
}

bool building_t::maybe_update_tape(point const &player_pos, bool end_of_tape) {
	if (!tape_manager.in_use) return 0;
	assert(has_room_geom());
	auto &decal_mgr(interior->room_geom->decal_manager);
	room_object_t const &obj(tape_manager.tape);
	float const thickness(obj.dz()), pad_dist(0.1*thickness);
	point const pos(player_pos + (1.5f*get_scaled_player_radius())*cview_dir);
	point sound_pos;
	float sound_gain(0.0);

	if (end_of_tape) { // add final tape
		if (tape_manager.points.empty()) return 0; // no tape
		decal_mgr.pend_tape_qbd.clear();
		point const end_pos(interior->find_closest_pt_on_obj_to_pos(*this, pos, pad_dist, 1)); // no_ceil_floor=1
		add_tape_quad(tape_manager.points.back(), end_pos, thickness, obj.color, decal_mgr.tape_qbd); // add final segment
		sound_pos = end_pos; sound_gain = 1.0;
	}
	else if (tape_manager.points.empty()) { // first point
		point const start_pos(interior->find_closest_pt_on_obj_to_pos(*this, pos, pad_dist, 1)); // starting position for tape; no_ceil_floor=1
		tape_manager.points.push_back(start_pos);
		decal_mgr.commit_pend_tape_qbd(); // commit any previous tape
		interior->room_geom->modified_by_player = 1; // make sure tape stays in this building
		sound_pos = start_pos; sound_gain = 1.0;
	}
	else {
		point last_pt(tape_manager.points.back()), p_int;

		if (!dist_less_than(last_pt, pos, thickness) && interior->line_coll(*this, last_pt, pos, p_int)) { // no short segments
			p_int += 0.5*thickness*(pos - last_pt).get_norm(); // move past the object to avoid an intersection at the starting point on the next call
			tape_manager.points.push_back(p_int);
			add_tape_quad(last_pt, p_int, thickness, obj.color, decal_mgr.tape_qbd);
			last_pt = p_int;
			//sound_pos = p_int; sound_gain = 0.1; // too noisy?
		}
		decal_mgr.pend_tape_qbd.clear();
		add_tape_quad(last_pt, pos, thickness, obj.color, decal_mgr.pend_tape_qbd);
		// update use count based on length change
		float const prev_dist(p2p_dist(last_pt, tape_manager.last_pos)), cur_dist(p2p_dist(last_pt, pos)), delta(cur_dist - prev_dist);
		int const delta_use_count(round_fp(0.5f*delta/thickness));
		if (!player_inventory.update_last_item_use_count(delta_use_count)) {tape_manager.clear();} // check if we ran out of tape
	}
	if (sound_gain > 0.0) {play_tape_sound(sound_pos, sound_gain);}
	tape_manager.last_pos = pos;
	return 1;
}

// returns the index of the first quad vertex, or -1 if no intersection found
int tape_quad_batch_draw::moving_vert_cyilin_int_tape(point &cur_pos, point const &prev_pos, float z1, float z2, float radius, float slow_amt) const {
	if (verts.empty()) return -1;
	if (cur_pos == prev_pos) return -1; // stopped, no effect
	assert(!(verts.size() % 12)); // must be a multiple of 12 (pairs of quads formed from two triangles)
	assert(slow_amt >= 0.0 && slow_amt <= 1.0); // 0.0 => no change to cur_pos, 1.0 => move back to prev_pos
	cylinder_3dw const cylin(point(cur_pos.x, cur_pos.y, z1), point(cur_pos.x, cur_pos.y, z2), radius, radius);

	for (unsigned i = 0; i < verts.size(); i += 12) { // iterate over pairs of back-to-back quads
		point const &p1(verts[i].v), &p2(verts[i+1].v); // get two opposite corners of the first quad, which approximates the line of the tape
		if (dist_xy_less_than(p1, p2, 2.0*radius)) continue; // skip short lines; ignore zval to skip already split vertical tape segments
		if (!line_intersect_cylinder(p1, p2, cylin, 1)) continue; // check_ends=1
		if (pt_line_dist(prev_pos, p1, p2) < pt_line_dist(cur_pos, p1, p2)) continue; // new point is further from the line - moving away, skip
		if (slow_amt > 0.0) {cur_pos = slow_amt*prev_pos + (1.0 - slow_amt)*cur_pos;}
		// hack to avoid intersection with just-placed tape: skip if prev pos also intersects when not slowing
		else if (line_intersect_cylinder(p1, p2, cylinder_3dw(point(prev_pos.x, prev_pos.y, z1), point(prev_pos.x, prev_pos.y, z2), radius, radius), 1)) continue;
		return i;
	}
	return -1; // not found
}
void tape_quad_batch_draw::split_tape_at(unsigned first_vert, point const &pos, float min_zval) {
	assert(first_vert+12 <= verts.size());
	assert(!(first_vert % 12)); // must be the start of a quad pair
	point p1(verts[first_vert].v), p2(verts[first_vert+1].v); // get two opposite corners of the first quad; copy to avoid invalid reference
	// find the lengths of both hanging segments, clipped to min_zval (the floor)
	float const d1(p2p_dist(p1, pos)), d2(p2p_dist(p2, pos)), t(d1/(d1 + d2)); // find split point
	float const width(p2p_dist(p1, verts[first_vert+2].v));
	float const len(p2p_dist(p1, p2)), len1(min(t*len, (p1.z - min_zval))), len2(min((1.0f - t)*len, (p2.z - min_zval))); // should be positive?
	// find the bottom points of the two hanging segments
	vector3d const dir((p2 - p1).get_norm());
	//p1 += 0.1*width*dir; p2 -= 0.1*width*dir; // shift all points slightly inward to avoid z-fighting
	point const p1b(p1.x, p1.y, p1.z-len1), p2b(p2.x, p2.y, p2.z-len2);
	vector3d wdir(cross_product(dir, plus_z).get_norm()); // segments hang with normal oriented toward the split point/along the original dir
	color_wrapper const cw(verts[first_vert]);

	if (len1 > 0.0) {
		unsigned const verts_start(verts.size());
		add_tape_quad(p1, p1b, width, cw, *this, wdir);
		assert(verts.size() == verts_start+12);
		for (unsigned i = 0; i < 12; ++i) {verts[first_vert+i] = verts[verts_start+i];} // overwrite old verts with new verts
		verts.resize(verts_start); // remove new verts
	}
	else { // zero length segment, overwrite with zeros to remove this quad
		for (unsigned i = 0; i < 12; ++i) {verts[first_vert+i].v = zero_vector;}
	}
	if (len2 > 0.0) {add_tape_quad(p2, p2b, width, cw, *this, -wdir);} // add the other segment
}

// Note: cur_pos.z should be between z1 and z2
void building_t::handle_vert_cylin_tape_collision(point &cur_pos, point const &prev_pos, float z1, float z2, float radius) const {
	if (!has_room_geom()) return;
	tape_quad_batch_draw &tape_qbd(interior->room_geom->decal_manager.tape_qbd); // Note: technically, this violates const-ness of this function
	// first, test if tape is very close, and if so, break it; otherwise, slow down the player/AI when colliding with tape
	int const vert_ix(tape_qbd.moving_vert_cyilin_int_tape(cur_pos, prev_pos, z1, z2, 0.2*radius, 0.0)); // 20% radius, slow_amt=0.0
	
	if (vert_ix >= 0) {
		tape_qbd.split_tape_at(vert_ix, cur_pos, z1); // min_zval=z1
		play_tape_sound(cur_pos, 1.0);
	}
	else {tape_qbd.moving_vert_cyilin_int_tape(cur_pos, prev_pos, z1, z2, radius, 0.85);} // slow_amt=0.85
}

// spraypaint, markers, and decals

bool line_int_cube_get_t(point const &p1, point const &p2, cube_t const &cube, float &tmin) {
	float tmin0(0.0), tmax0(1.0);
	if (get_line_clip(p1, p2, cube.d, tmin0, tmax0) && tmin0 < tmin) {tmin = tmin0; return 1;}
	return 0;
}
bool line_int_cubes_get_t(point const &p1, point const &p2, vect_cube_t const &cubes, float &tmin, cube_t &target) {
	bool had_int(0);

	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (line_int_cube_get_t(p1, p2, *c, tmin)) {target = *c; had_int = 1;}
	}
	return had_int;
}
vector3d get_normal_for_ray_cube_int_xy(point const &p, cube_t const &c, float tolerance) {
	vector3d n(zero_vector);

	for (unsigned d = 0; d < 2; ++d) { // find the closest intersecting cube XY edge, which will determine the normal vector
		if (fabs(p[d] - c.d[d][0]) < tolerance) {n[d] = -1.0; break;} // test low  edge
		if (fabs(p[d] - c.d[d][1]) < tolerance) {n[d] =  1.0; break;} // test high edge
	}
	return n;
}

// decals

class paint_manager_t : public paint_draw_t { // for paint on exterior walls/windows, viewed from inside the building
	building_t const *paint_bldg = nullptr;
public:
	bool have_paint_for_building() const { // only true if the building contains the player
		return (paint_bldg && !(qbd[0].empty() && qbd[1].empty()) && paint_bldg->bcube.contains_pt(get_camera_building_space()));
	}
	quad_batch_draw &get_paint_qbd(building_t const *const building, bool is_marker) {
		if (building != paint_bldg) { // paint switches to this building
			for (unsigned d = 0; d < 2; ++d) {qbd[d].clear();}
			paint_bldg = building;
		}
		return qbd[is_marker];
	}
};

paint_manager_t ext_paint_manager;
bool have_buildings_ext_paint() {return ext_paint_manager.have_paint_for_building();}
void draw_buildings_ext_paint() {ext_paint_manager.draw_paint();}

bool building_t::apply_paint(point const &pos, vector3d const &dir, colorRGBA const &color, room_object const obj_type) const { // spraypaint or marker
	bool const is_spraypaint(obj_type == TYPE_SPRAYCAN), is_marker(obj_type == TYPE_MARKER);
	assert(is_spraypaint || is_marker); // only these two are supported
	// find intersection point and normal; assumes pos is inside the building
	assert(has_room_geom());
	float const max_dist((is_spraypaint ? 16.0 : 3.0)*CAMERA_RADIUS), tolerance(0.01*get_wall_thickness());
	point const pos2(pos + max_dist*dir);
	float tmin(1.0);
	vector3d normal;
	cube_t target;
	
	for (unsigned d = 0; d < 2; ++d) {
		if (line_int_cubes_get_t(pos, pos2, interior->walls[d], tmin, target)) {
			normal    = zero_vector;
			normal[d] = -SIGN(dir[d]); // normal is opposite of ray dir in this dim
		}
	}
	if (line_int_cubes_get_t(pos, pos2, interior->floors  , tmin, target)) {normal =  plus_z;}
	if (line_int_cubes_get_t(pos, pos2, interior->ceilings, tmin, target)) {normal = -plus_z;}
	
	// include exterior walls; okay to add spraypaint and markers over windows
	cube_t const part(get_part_containing_pt(pos));
	float tmin0(0.0), tmax0(1.0);
	bool exterior_wall(0);

	if (get_line_clip(pos, pos2, part.d, tmin0, tmax0) && tmax0 < tmin) { // part edge is the closest intersection point
		// check other parts to see if ray continues into them; if not, it exited the building; this implementation isn't perfect but should be close enough
		point const cand_p_int(pos + tmax0*(pos2 - pos));
		bool found(0);

		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (*p == part || !p->contains_pt_exp(cand_p_int, tolerance)) continue; // ray does not continue into this new part
			if (check_line_clip(cand_p_int, pos2, p->d)) {found = 1; break;} // ray continues into this part
		}
		if (!found) { // ray has exited the building
			vector3d const n(-get_normal_for_ray_cube_int_xy(cand_p_int, part, tolerance)); // negate the normal because we're looking for the exit point from the cube
			if (n != zero_vector) {tmin = tmax0; normal = n; target = part; exterior_wall = 1;}
		}
	}
	// check for rugs, pictures, and whiteboards, which can all be painted over; also check for walls from closets
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
	bool const is_wall(normal.x != 0.0 || normal.y != 0.0), is_floor(normal == plus_z);
	bool walls_blocked(0);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
 		if ((is_wall && (i->type == TYPE_PICTURE || i->type == TYPE_WBOARD || i->type == TYPE_MIRROR)) || (is_floor && (i->type == TYPE_RUG || i->type == TYPE_FLOORING))) {
			if (line_int_cube_get_t(pos, pos2, *i, tmin)) {target = *i;} // Note: return value is ignored, we only need to update tmin and target; normal should be unchanged
		}
		else if (i->type == TYPE_CLOSET && line_int_cube_get_t(pos, pos2, *i, tmin)) {
			point const cand_p_int(pos + tmin*(pos2 - pos));
			normal = get_normal_for_ray_cube_int_xy(cand_p_int, *i, tolerance); // should always return a valid normal
			target = *i;
		}
		else if (i->type == TYPE_STALL || i->type == TYPE_CUBICLE) {
			cube_t c(*i);

			if (i->type == TYPE_STALL && i->shape != SHAPE_SHORT) { // toilet stall, clip cube to wall height
				float const dz(c.dz());
				c.z2() -= 0.35*dz; c.z1() += 0.15*dz;
			}
			float tmin0(tmin);
			if (!line_int_cube_get_t(pos, pos2, c, tmin0)) continue;
			if (i->contains_pt(pos)) continue; // inside stall/cubicle, can't paint the exterior
			vector3d const n(get_normal_for_ray_cube_int_xy((pos + tmin0*(pos2 - pos)), c, tolerance)); // should always return a valid normal
			if (n[i->dim] != 0) {walls_blocked = 1; continue;} // only the side walls count
			tmin = tmin0; normal = n; target = c;
		}
	} // for i
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		float tmin0(tmin);
		if (!line_int_cube_get_t(pos, pos2, *i, tmin0)) continue;
		if (i->contains_pt(pos)) {walls_blocked = 1; continue;} // can't spraypaint the outside of the elevator when standing inside it
		vector3d const n(get_normal_for_ray_cube_int_xy((pos + tmin0*(pos2 - pos)), *i, tolerance)); // should always return a valid normal
		if (n[i->dim] == (i->dir ? 1.0 : -1.0)) {walls_blocked = 1; continue;} // skip elevator opening, even if not currently open
		tmin = tmin0; normal = n; target = *i;
	}
	for (auto i = interior->stairwells.begin(); i != interior->stairwells.end(); ++i) {
		if (i->shape == SHAPE_STRAIGHT) continue; // no walls, skip
		// expand by wall half-width; see building_t::add_stairs_and_elevators()
		float const step_len_pos(i->get_sz_dim(i->dim)/i->get_num_stairs());
		cube_t c(*i);
		c.expand_by_xy(0.15*step_len_pos); // wall half width
		float tmin0(tmin);
		if (!line_int_cube_get_t(pos, pos2, c, tmin0)) continue;
		if (c.contains_pt(pos)) {walls_blocked = 1; continue;} // can't spraypaint the outside of the stairs when standing inside them
		vector3d const n(get_normal_for_ray_cube_int_xy((pos + tmin0*(pos2 - pos)), c, tolerance)); // should always return a valid normal

		if (i->shape == SHAPE_U) {
			if (n[i->dim] == (i->dir ? -1.0 : 1.0)) {walls_blocked = 1; continue;} // skip stairs opening
		}
		else if (i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES) {
			// Note: we skip the end for SHAPE_WALLED and only check the sides because it depends on the floor we're on
			if (n[i->dim] != 0) {walls_blocked = 1; continue;} // skip stairs opening, either side
		}
		else {assert(0);} // unsupported stairs type
		tmin = tmin0; normal = n; target = c;
	}
	if (normal == zero_vector) return 0; // no walls, ceilings, floors, etc. hit
	if (walls_blocked && normal.z == 0.0) return 0; // can't spraypaint walls through elevator, stairs, etc.
	point p_int(pos + tmin*(pos2 - pos));
	if (check_line_intersect_doors(pos, p_int)) return 0; // blocked by door, no spraypaint; can't add spraypaint over door in case door is opened
	float const max_radius((is_spraypaint ? 2.0 : 0.035)*CAMERA_RADIUS);
	float const dist(p2p_dist(pos, p_int)), radius(is_spraypaint ? min(max_radius, max(0.05f*max_radius, 0.1f*dist)) : max_radius); // modified version of get_spray_radius()
	float const alpha((is_spraypaint && radius > 0.5*max_radius) ? (1.0 - (radius - 0.5*max_radius)/max_radius) : 1.0); // 0.5 - 1.0
	p_int += 0.01*radius*normal; // move slightly away from the surface
	assert(bcube.contains_pt(p_int));
	unsigned const dim(get_max_dim(normal)), d1((dim+1)%3), d2((dim+2)%3);

	// check that entire circle is contained in the target
	for (unsigned e = 0; e < 2; ++e) {
		unsigned const d(e ? d2 : d1);
		if (p_int[d] - 0.9*radius < target.d[d][0] || p_int[d] + 0.9*radius > target.d[d][1]) return 0; // extends outside the target surface in this dim
	}
	static point last_p_int(all_zeros);
	if (dist_less_than(p_int, last_p_int, 0.25*radius)) return 1; // too close to previous point, skip (to avoid overlapping sprays at the same location); still return 1
	last_p_int = p_int;
	vector3d dir1, dir2; // unit vectors
	dir1[d1] = 1.0; dir2[d2] = 1.0;
	float const winding_order_sign(-SIGN(normal[dim])); // make sure to invert the winding order to match the normal sign
	// Note: interior spraypaint draw uses back face culling while exterior draw does not; invert the winding order for exterior quads so that they show through windows correctly
	vector3d const dx(radius*dir1*winding_order_sign*(exterior_wall ? -1.0 : 1.0));
	interior->room_geom->decal_manager.paint_draw[exterior_wall].qbd[is_marker].add_quad_dirs(p_int, dx, radius*dir2, colorRGBA(color, alpha), normal); // add interior/exterior paint
	if (exterior_wall) {ext_paint_manager.get_paint_qbd(this, is_marker).add_quad_dirs(p_int, dx, radius*dir2, colorRGBA(color, alpha), normal);} // add exterior paint only
	static double next_sound_time(0.0);

	if (tfticks > next_sound_time) { // play sound if sprayed/marked, but not too frequently; marker has no sound
		gen_sound_thread_safe_at_player((is_spraypaint ? (int)SOUND_SPRAY : (int)SOUND_SQUEAK), 0.25);
		if (is_spraypaint) {register_building_sound(pos, 0.1);}
		next_sound_time = tfticks + double(is_spraypaint ? 0.5 : 0.25)*TICKS_PER_SECOND;
	}
	player_inventory.record_damage_done(is_spraypaint ? 1.0 : 0.1); // spraypaint does more damage than markers
	return 1;
}

bool room_object_t::can_use() const { // excludes dynamic objects
	return (type == TYPE_SPRAYCAN || type == TYPE_MARKER || type == TYPE_TPROLL || type == TYPE_BOOK || type == TYPE_PHONE || type == TYPE_TAPE);
}
bool room_object_t::can_place_onto() const {
	return (type == TYPE_TABLE || type == TYPE_DESK || type == TYPE_DRESSER || type == TYPE_NIGHTSTAND || type == TYPE_COUNTER || type == TYPE_KSINK ||
		type == TYPE_BRSINK || type == TYPE_BED || type == TYPE_BOX || type == TYPE_CRATE || type == TYPE_KEYBOARD || type == TYPE_BOOK);
}

bool building_t::apply_toilet_paper(point const &pos, vector3d const &dir, float half_width) {
	// for now, just drop a square of TP on the floor; could do better; should the TP roll shrink in size as this is done?
	assert(has_room_geom());
	static point last_tp_pos;
	if (dist_xy_less_than(pos, last_tp_pos, 1.5*half_width)) return 0; // too close to prev pos
	last_tp_pos = pos;
	float zval(pos.z);
	if (!get_zval_for_obj_placement(pos, half_width, zval, 1)) return 0; // no suitable placement found; add_z_bias=1
	vector3d d1(dir.x, dir.y, 0.0);
	if (d1 == zero_vector) {d1 = plus_x;} else {d1.normalize();}
	vector3d d2(cross_product(d1, plus_z));
	if (d2 == zero_vector) {d2 = plus_y;} else {d2.normalize();}
	interior->room_geom->decal_manager.tp_qbd.add_quad_dirs(point(pos.x, pos.y, zval), d1*half_width, d2*half_width, WHITE, plus_z);
	interior->room_geom->modified_by_player = 1; // make sure TP stays in this building
	// Note: no damage done for TP
	return 1;
}

void building_t::add_blood_decal(point const &pos) {
	assert(has_room_geom());
	float const radius(get_scaled_player_radius());
	float zval(pos.z);
	if (!get_zval_of_floor(pos, radius, zval)) return; // no suitable floor found
	tex_range_t const tex_range(tex_range_t::from_atlas((rand()&1), (rand()&1), 2, 2)); // 2x2 texture atlas
	interior->room_geom->decal_manager.blood_qbd.add_quad_dirs(point(pos.x, pos.y, zval), -plus_x*radius, plus_y*radius, WHITE, plus_z, tex_range); // Note: never cleared
	interior->room_geom->modified_by_player = 1; // make sure blood stays in this building
	player_inventory.record_damage_done(100.0); // blood is a mess to clean up (though damage will be reset on player death anyway)
}

// sound/audio tracking

void register_building_sound(point const &pos, float volume) {
	if (volume == 0.0 || !(show_bldg_pickup_crosshair || in_building_gameplay_mode())) return; // only when in gameplay/item pickup mode
	assert(volume > 0.0); // can't be negative
#pragma omp critical(building_sounds_update)
	{ // since this can be called by both the draw thread and the AI update thread, it should be in a critical section
		if (volume > ALERT_THRESH && cur_sounds.size() < 100) { // cap at 100 sounds in case they're not being cleared
			float const max_merge_dist(0.5*CAMERA_RADIUS);
			bool merged(0);

			for (auto i = cur_sounds.begin(); i != cur_sounds.end(); ++i) { // attempt to merge with an existing nearby sound
				if (dist_less_than(pos, i->pos, max_merge_dist)) {i->radius += volume; merged = 1;}
			}
			if (!merged) {cur_sounds.emplace_back(pos, volume);} // Note: volume is stored in radius field of sphere_t
		}
		cur_building_sound_level += volume;
	}
}
void register_building_sound_at_player(float volume) {
	register_building_sound(get_camera_building_space(), 1.0);
}

bool get_closest_building_sound(point const &at_pos, point &sound_pos, float floor_spacing) {
	if (cur_sounds.empty()) return 0;
	float max_vol(0.0); // 1.0 at a sound=1.0 volume at a distance of floor_spacing

	for (auto i = cur_sounds.begin(); i != cur_sounds.end(); ++i) {
		float vol(i->radius/max(0.01f*floor_spacing, p2p_dist(i->pos, at_pos)));
		if (fabs(i->pos.z - at_pos.z) > 0.75f*floor_spacing) {vol *= 0.5;} // half the volume when the sound comes from another floor
		if (vol > max_vol) {max_vol = vol; sound_pos = i->pos;}
	} // for i
	//cout << TXT(cur_sounds.size()) << TXT(max_vol) << endl;
	return (max_vol*floor_spacing > 0.06f);
}

void maybe_play_zombie_sound(point const &sound_pos_bs, unsigned zombie_ix, bool alert_other_zombies, bool high_priority) {
	unsigned const NUM_ZSOUNDS = 5;
	static rand_gen_t rgen;
	static double next_time_all(0.0), next_times[NUM_ZSOUNDS] = {};
	if (!high_priority && tfticks < next_time_all) return; // don't play any sound too frequently
	if (!high_priority && (rgen.rand()&3) != 0)    return; // only generate a sound 25% of the time (each frame), to allow more than one zombie to get a chance
	unsigned const sound_id(zombie_ix%NUM_ZSOUNDS); // choose one of the zombie sounds, determined by the current zombie
	double &next_time(next_times[sound_id]);
	if (!high_priority && tfticks < next_time) return; // don't play this particular sound too frequently
	next_time_all = tfticks + double(rgen.rand_uniform(1.0, 2.0))*TICKS_PER_SECOND; // next sound of any  type can play between 0.8 and 2.0s in the future
	next_time     = tfticks + double(rgen.rand_uniform(2.5, 5.0))*TICKS_PER_SECOND; // next sound of this type can play between 2.5 and 5.0s in the future
	gen_sound_thread_safe((SOUND_ZOMBIE1 + sound_id), (sound_pos_bs + get_camera_coord_space_xlate()));
	if (alert_other_zombies) {register_building_sound(sound_pos_bs, 0.4);}
}

void water_sound_manager_t::register_running_water(room_object_t const &obj, building_t const &building) {
	if (!(obj.flags & RO_FLAG_IS_ACTIVE)) return; // not turned on
	if (fabs(obj.z2() - camera_bs.z) > building.get_window_vspace()) return; // on the wrong floor
	point const pos(obj.get_cube_center());
	float const dsq(p2p_dist_sq(pos, camera_bs));
	if (dmin_sq == 0.0 || dsq < dmin_sq) {closest = building.local_to_camera_space(pos); dmin_sq = dsq;}
}
void water_sound_manager_t::finalize() {
	if (dmin_sq == 0.0) return; // no water found
	static point prev_closest;
	bool const skip_if_already_playing(closest == prev_closest); // don't reset sound loop unless it moves to a different sink
	prev_closest = closest;
	gen_sound_thread_safe(SOUND_SINK, closest, 1.0, 1.0, 0.06, skip_if_already_playing); // fast distance falloff; will loop at the end if needed
}

// gameplay logic

bool player_has_room_key() {return player_inventory.player_has_key();}

// return value: 0=no effect, 1=player is killed, 2=this person is killed
int register_ai_player_coll(bool &has_key, float height) {
	if (do_room_obj_pickup && player_inventory.take_person(has_key, height)) {
		gen_sound_thread_safe_at_player(SOUND_ITEM, 0.5);
		do_room_obj_pickup = 0; // no more object pickups
		return 2;
	}
	static double last_coll_time(0.0);
	
	if (tfticks - last_coll_time > 2.0*TICKS_PER_SECOND) {
		gen_sound_thread_safe_at_player(SOUND_SCREAM1);
		last_coll_time = tfticks;
	}
	add_camera_filter(colorRGBA(RED, 0.25), 1, -1, CAM_FILT_DAMAGE); // 1 tick of red damage
	player_inventory.take_damage(0.04*fticks); // take damage over time
	
	if (player_inventory.player_is_dead()) {
		if (player_has_room_key()) {has_key = 1;}
		return 1;
	}
	return 0;
}

void building_gameplay_action_key(int mode, bool mouse_wheel) {
	if (camera_in_building) { // building interior action
		if (mouse_wheel) {player_inventory.switch_item(mode != 0);}
		else if (mode == 0) {building_action_key    = 1;} // 'q'
		else if (mode == 1) {do_room_obj_pickup     = 1;} // 'e'
		else if (mode == 2) {use_last_pickup_object = 1;} // 'E'
		else if (mode == 3) {building_action_key    = 2;} // 'r'
		else {assert(0);} // unsupported key/mode
		show_bldg_pickup_crosshair = 1; // show crosshair on first pickup, too difficult to pick up objects without it
	}
	else { // building exterior/city/road/car action
		if (mode == 1 || mode == 2) {city_action_key = 1;} // 'e'/'E'
		else {} // 'q'/'r'
	}
}

void attenuate_rate(float &v, float rate) {
	if (v != 0.0) {
		v *= exp(-rate*fticks); // exponential slowdown
		if (fabs(v) < 0.001) {v = 0.0;} // stop moving
	}
}

void building_gameplay_next_frame() {
	attenuate_rate(office_chair_rot_rate, 0.05); // update office chair rotation

	if (in_building_gameplay_mode()) { // run gameplay update logic
		show_bldg_pickup_crosshair = 1;
		// update sounds used by AI
		auto i(cur_sounds.begin()), o(i);

		for (; i != cur_sounds.end(); ++i) {
			i->radius *= exp(-0.04*fticks);
			if (i->radius > ALERT_THRESH) {*(o++) = *i;} // keep if above thresh
		}
		cur_sounds.erase(o, cur_sounds.end());
	}
	player_held_object = carried_item_t();
	player_inventory.next_frame();
	// reset state for next frame
	cur_building_sound_level = min(1.2f, max(0.0f, (cur_building_sound_level - 0.01f*fticks))); // gradual decrease
	can_pickup_bldg_obj = 0;
	do_room_obj_pickup  = city_action_key = can_do_building_action = 0;
}

void enter_building_gameplay_mode() {player_inventory.clear_all();}

