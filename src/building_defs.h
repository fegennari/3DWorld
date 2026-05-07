// 3D World - Buildings Shared Definitions
// by Frank Gennari
// 5-6-26
#pragma once

#include "3DWorld.h"
#include "gl_ext_arb.h" // for vbo_wrap_t, etc.
#include "lightmap.h"
#include "file_utils.h" // for kw_to_val_map_t

bool const EXACT_MULT_FLOOR_HEIGHT = 1;
bool const ENABLE_MIRROR_REFLECTIONS = 1;
bool const ENABLE_GLASS_FLOOR_REF  = 1;
bool const DRAW_CITY_INT_WINDOWS   = 1; // requires having different window x/y size/space/offset values for interior vs. exterior windows
bool const ADD_WALKWAY_EXT_DOORS   = 1; // requires DRAW_WALKWAY_INTERIORS=1 in gen_buildings.cpp
unsigned const MAX_CYLIN_SIDES     = 36;
unsigned const MAX_DRAW_BLOCKS     = 10; // for building interiors only; currently have floor, ceiling, walls, and doors with a few different materials
unsigned const NUM_STAIRS_PER_FLOOR= 12;
unsigned const NUM_STAIRS_PER_FLOOR_U = 16;
unsigned const NUM_STAIRS_PER_FLOOR_L = 12;
unsigned const NUM_STAIRS_PER_FLOOR_ESC = 12;
unsigned const MAX_SHELVES         = 4;
float const FLOOR_THICK_VAL_HOUSE  = 0.10; // 10% of floor spacing
float const FLOOR_THICK_VAL_OFFICE = 0.11; // thicker for office buildings
float const FLOOR_THICK_VAL_WINDOWLESS = 0.12; // even thicker for windowless office buildings
float const RAMP_THICKNESS_SCALE   = 0.11; // thickness to height ratio
float const WINDOW_BORDER_MULT     = 0.94; // account for the frame part of the window texture, which is included in the interior cutout of the window
float const WALL_THICK_VAL         = 0.05; // 5% of floor spacing
float const DOOR_THICK_TO_WIDTH    = 0.04; // ratio of door thickness to width for doors opening to the side
float const DEF_CITY_MIN_ALPHA     = 0.01;
float const DOOR_WIDTH_SCALE       = 0.5;
float const DOOR_WIDTH_SCALE_OFFICE= 0.7; // wider than house doors
float const STAIRS_WALL_WIDTH_MULT = 0.15; // relative to the depth of a stair
float const ELEVATOR_Z2_SHIFT      = 0.6; // shift downward, relative to ceiling thickness
float const ELEVATOR_STAND_DIST    = 0.4; // distance that people stand from the elevator door relative to floor spacing
float const ELEVATOR_FRAME_WIDTH   = 0.2; // frame width to each side of door relative to elevator width
float const DOOR_FRAME_WIDTH       = 0.07; // for door texture, relative to door width
float const EXT_BASEMENT_JOIN_DIST = 4.0; // relative to floor spacing
float const BALCONY_PILLAR_SCALE   = 0.15; // relative to depth
float const BASEMENT_ENTRANCE_SCALE= 0.33;
float const CHAIR_LEG_WIDTH        = 0.15; // relative to chair width
float const CHAIR_LEG_WIDTH_MALL   = 0.07; // relative to chair width
float const BED_HEAD_WIDTH         = 0.04; // headboard width; relative to bed width
float const SHELF_RACK_HEIGHT_FS   = 0.85*(1.0 - FLOOR_THICK_VAL_OFFICE);
float const WHEATER_PIPE_SPACING   = 0.65; // relative to radius
float const WHEATER_PIPE_H_DIST    = 0.92; // relative to radius
float const PLANT_POT_RADIUS       = 0.75; // relative to plant bcube radius
float const DEF_NORM_BIAS_SCALE    = 10.0; // see shadow_map.part
float const CITY_BIAS_SCALE        = 0.1;
float const RETAIL_SMAP_DSCALE     = 0.7; // player shadow caster distance for retail rooms and malls relative to light radius
float const PERSON_INT_SMAP_DSCALE = 0.8; // building person shadow caster distance relative to light radius
float const ESCALATOR_SPEED        = 1.0/TICKS_PER_SECOND; // in steps per tick
float const MALL_FLOOR_HEIGHT      = 2.0; // as a multiple of normal building floor height
float const CEILING_BEAM_THICK     = 2.5; // as a multiple of wall thickness
float const BACKSPLASH_HEIGHT      = 0.33; // relative to cabinet height
float const LOCKER_BOT_SHELF_HEIGHT= 0.67; // relative to full height
float const GLASS_IOR              = 1.6;
float const RESTROOM_WIN_TSCALE_X  = 12.0; // use a consistent scale to avoid single long windows

unsigned const BOTTLE_EMPTY_MASK= 192; // empty if both bits 6 and 7 are set
unsigned const NUM_CHAIR_COLORS = 13;
unsigned const MAX_BCASE_BOOKS  = 48; // limited by available bit flags
unsigned const NUM_BOOK_COLORS  = 16;
unsigned const NUM_PAPER_COLORS = 6;
unsigned const NUM_SPCAN_COLORS = 11;
unsigned const NUM_LAMP_COLORS  = 6;
unsigned const NUM_TCAN_COLORS  = 6;
unsigned const NUM_TAPE_COLORS  = 7;
unsigned const NUM_SHIRT_COLORS = 14;
unsigned const NUM_STAPLER_COLORS = 5;
unsigned const NUM_TSHIRT_COLORS  = 9;
unsigned const NUM_TOASTER_COLORS = 7;
unsigned const NUM_TRASH_COLORS   = 8;
unsigned const NUM_PFLOAT_COLORS  = 6;
unsigned const NUM_HARDHAT_COLORS = 9;
unsigned const NUM_HANDLE_COLORS  = 4;
unsigned const NUM_SOAP_COLORS    = 5;
unsigned const NUM_SPICE_COLORS   = 7;
unsigned const NUM_FOOD_TUB_COLORS= 3;
unsigned const NUM_RAND_SIGN_BKG_COLORS = 8; // Note: 0 is auto
enum {SIGN_BKG_COLOR_BLACK=NUM_RAND_SIGN_BKG_COLORS, SIGN_BKG_COLOR_RED, NUM_SIGN_BKG_COLORS};
unsigned const NUM_MALL_CHAIR_COLORS  = 5;
unsigned const NUM_SP_EMISSIVE_COLORS = 2;
colorRGBA const GD_SP_COLOR(0.5, 1.0, 1.0); // used for glow-in-the-dark spraypaint
colorRGBA const chair_colors[NUM_CHAIR_COLORS] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN, DK_BROWN};
colorRGBA const book_colors [NUM_BOOK_COLORS ] = {GRAY_BLACK, WHITE, LT_GRAY, GRAY, DK_GRAY, DK_BLUE, BLUE, LT_BLUE, DK_RED, RED, ORANGE, YELLOW, DK_GREEN, LT_BROWN, BROWN, DK_BROWN};
colorRGBA const spcan_colors[NUM_SPCAN_COLORS] = {GD_SP_COLOR, WHITE, RED, GREEN, BLUE, YELLOW, PINK, ORANGE, PURPLE, BROWN, BLACK};
colorRGBA const sp_emissive_colors[NUM_SP_EMISSIVE_COLORS] = {colorRGBA(0.2, 1.0, 0.2), colorRGBA(0.2, 0.5, 1.0)}; // light green, greenish blue
colorRGBA const lamp_colors[NUM_LAMP_COLORS]   = {WHITE, GRAY_BLACK, BROWN, LT_BROWN, DK_BROWN, OLIVE};
colorRGBA const cream(0.9, 0.9, 0.8), vlt_yellow(1.0, 1.0, 0.5);
colorRGBA const paper_colors[NUM_PAPER_COLORS] = {WHITE, WHITE, WHITE, cream, cream, vlt_yellow};
colorRGBA const pen_colors    [4] = {WHITE, BLACK, colorRGBA(0.2, 0.4, 1.0), RED};
colorRGBA const pencil_colors [2] = {colorRGBA(1.0, 0.75, 0.25), colorRGBA(1.0, 0.5, 0.1)};
colorRGBA const marker_colors [8] = {BLACK, RED, BLACK, BLUE, BLACK, GREEN, RED, PURPLE};
colorRGBA const railing_colors[3] = {GOLD, DK_BROWN, BLACK};
colorRGBA const tcan_colors   [NUM_TCAN_COLORS   ] = {BLUE, DK_GRAY, LT_GRAY, GRAY, BLUE, WHITE};
colorRGBA const tape_colors   [NUM_TAPE_COLORS   ] = {GRAY, GRAY, GRAY, GRAY, BKGRAY, BKGRAY, colorRGBA(0.2, 0.2, 1.0)}; // gray duct tape is the most common
colorRGBA const shirt_colors  [NUM_SHIRT_COLORS  ] = {WHITE, WHITE, WHITE, BKGRAY, GRAY, RED, DK_RED, LT_BLUE, BLUE, DK_BLUE, DK_GREEN, DK_BROWN, BROWN, ORANGE};
colorRGBA const stapler_colors[NUM_STAPLER_COLORS] = {BLACK, RED, BLACK, BLUE, BLACK};
colorRGBA const TSHIRT_COLORS [NUM_TSHIRT_COLORS ] = {WHITE, BKGRAY, GRAY, RED, GREEN, BLUE, YELLOW, ORANGE, WHITE};
colorRGBA const toaster_colors[NUM_TOASTER_COLORS] = {WHITE, LT_GRAY, GRAY, DK_GRAY, GRAY_BLACK, colorRGBA(0.0, 0.0, 0.5), colorRGBA(0.5, 0.0, 0.0)};
colorRGBA const trash_colors  [NUM_TRASH_COLORS  ] = {WHITE, WHITE, WHITE, WHITE, cream, vlt_yellow, LT_GRAY, colorRGBA(0.8, 0.6, 0.4)};
colorRGBA const pfloat_colors [NUM_PFLOAT_COLORS ] = {WHITE, YELLOW, PINK, GREEN, ORANGE, LT_BLUE};
colorRGBA const hardhat_colors[NUM_HARDHAT_COLORS] = {YELLOW, YELLOW, YELLOW, YELLOW, ORANGE, ORANGE, RED, WHITE, colorRGBA(0.25, 0.25, 1.0, 1.0)};
colorRGBA const handle_colors [NUM_HANDLE_COLORS ] = {DK_RED, colorRGBA(0.1, 0.2, 0.4), colorRGBA(0.05, 0.3, 0.05), BKGRAY};
colorRGBA const soap_colors   [NUM_SOAP_COLORS   ] = {WHITE, cream, vlt_yellow, colorRGBA(1.0, 0.8, 0.6), colorRGBA(0.7, 1.0, 0.7)};
colorRGBA const spice_colors  [NUM_SPICE_COLORS  ] = {WHITE, BLACK, LT_BROWN, BROWN, DK_BROWN, DK_GREEN, OLIVE};
colorRGBA const food_tub_colors[NUM_FOOD_TUB_COLORS] = {WHITE, cream, colorRGBA(0.8, 0.7, 0.5)}; // tan
colorRGBA const sign_bkg_colors[NUM_SIGN_BKG_COLORS] = {WHITE, CYAN, MAGENTA, YELLOW, LT_RED, LT_BLUE, LT_GREEN, LT_BROWN, BLACK, RED};
colorRGBA const mall_chair_colors[NUM_MALL_CHAIR_COLORS] = {WHITE, LT_GRAY, GRAY, ORANGE, LT_BROWN};
colorRGBA const LAMP_COLOR(1.0, 0.8, 0.6); // soft white
colorRGBA const WALL_LAMP_COLOR(1.0, 0.9, 0.8);
colorRGBA const WOOD_COLOR(0.9, 0.7, 0.5); // light brown, multiplies wood texture color; typical value to use
colorRGBA const DUCT_COLOR(WHITE);
colorRGBA const rat_color(GRAY); // make the rat's fur darker
colorRGBA const candle_color(0.95, 0.9, 0.75, 1.0); // cream
colorRGBA const GLASS_COLOR(0.8, 1.0, 0.9, 0.25);
colorRGBA const BREAKER_PANEL_COLOR(0.50, 0.60, 0.70);
colorRGBA const WATER_HEATER_COLOR (0.50, 0.55, 0.60); // slightly blue-green tinted gray
colorRGBA const DARK_BRASS_C(0.4, 0.35, 0.15, 1.0);
colorRGBA const mall_tc_legs_color(BKGRAY);

unsigned const NUM_LOCK_COLORS = 8;
unsigned const MAX_LOCK_INDEX  = NUM_LOCK_COLORS + 2;
colorRGBA   const lock_colors     [NUM_LOCK_COLORS] = {WHITE, BLACK, RED, GREEN, BLUE, YELLOW, ORANGE, BROWN};
std::string const lock_color_names[NUM_LOCK_COLORS] = {"silver", "black", "red", "green", "blue", "yellow", "orange", "brown"};


// building types/functions; these are for primary buildings, not basements/rooms (such as malls or parking garages)
enum {BTYPE_UNSET=0, BTYPE_HOUSE, BTYPE_MULT_FAM, BTYPE_OFFICE, BTYPE_APARTMENT, BTYPE_HOTEL, BTYPE_HOSPITAL, BTYPE_PARKING, BTYPE_MALL, BTYPE_FACTORY,
	BTYPE_WAREHOUSE, BTYPE_POWERPLANT, BTYPE_SCHOOL, BTYPE_POLICE, BTYPE_FIRE_STAT, BTYPE_PRISON, BTYPE_RESTAURANT, BTYPE_CONV_STORE, BTYPE_RESTROOM, NUM_BUILDING_TYPES};
std::string const btype_names[NUM_BUILDING_TYPES] =
{"", "House", "Multi-Family House", "Office", "Apartments", "Hotel", "Hospital", "Parking", "Mall", "Factory", "Warehouse", "Power Plant",
"School", "Police Station", "Fire Station", "Prison", "Restaurant", "Store", "Restroom"};
colorRGBA const  btype_colors[NUM_BUILDING_TYPES] =
{WHITE, WHITE, YELLOW,               WHITE,    GREEN,        GREEN,   BLUE,       BROWN,     ORANGE, RED,       RED,         RED,
PURPLE,   MAGENTA,          MAGENTA,        BLACK,    CYAN,         WHITE,   PINK};
typedef uint8_t building_type_t;

enum { // room object types
	TYPE_NONE=0, TYPE_TABLE, TYPE_CHAIR, TYPE_STAIR, TYPE_STAIR_WALL, TYPE_ELEVATOR, TYPE_LIGHT, TYPE_RUG, TYPE_PICTURE, TYPE_WBOARD,
	TYPE_BOOK, TYPE_BCASE, TYPE_TCAN, TYPE_DESK, TYPE_BED, TYPE_WINDOW, TYPE_BLOCKER, TYPE_COLLIDER, TYPE_CUBICLE, TYPE_STALL,
	TYPE_SIGN, TYPE_COUNTER, TYPE_CABINET, TYPE_KSINK, TYPE_BRSINK, TYPE_PLANT, TYPE_DRESSER, TYPE_NIGHTSTAND, TYPE_FLOORING, TYPE_CLOSET,
	TYPE_WALL_TRIM, TYPE_RAILING, TYPE_CRATE, TYPE_BOX, TYPE_MIRROR, TYPE_SHELVES, TYPE_KEYBOARD, TYPE_SHOWER, TYPE_RDESK, TYPE_BOTTLE,
	TYPE_WINE_RACK, TYPE_COMPUTER, TYPE_MWAVE, TYPE_PAPER, TYPE_BLINDS, TYPE_PEN, TYPE_PENCIL, TYPE_PAINTCAN, TYPE_LG_BALL, TYPE_HANGER_ROD,
	TYPE_DRAIN, TYPE_MONEY, TYPE_PHONE, TYPE_TPROLL, TYPE_SPRAYCAN, TYPE_MARKER, TYPE_BUTTON, TYPE_VENT_HOOD, TYPE_SWITCH, TYPE_PLATE,
	TYPE_LAPTOP, TYPE_FPLACE, TYPE_LBASKET, TYPE_WHEATER, TYPE_TAPE, TYPE_OUTLET, TYPE_PG_WALL, TYPE_PG_PILLAR, TYPE_PG_BEAM, TYPE_PARK_SPACE,
	TYPE_RAMP, TYPE_PIPE, TYPE_CURB, TYPE_BRK_PANEL, TYPE_VENT, TYPE_BREAKER, TYPE_FURNACE, TYPE_ATTIC_DOOR, TYPE_CHIMNEY, TYPE_DUCT,
	TYPE_TOY, TYPE_DRESS_MIR, TYPE_PAN, TYPE_VASE, TYPE_URN, TYPE_FCABINET, TYPE_STAPLER, TYPE_WIND_SILL, TYPE_BALCONY, TYPE_SPRINKLER,
	TYPE_FEXT_MOUNT, TYPE_FEXT_SIGN, TYPE_PIZZA_BOX, TYPE_PIZZA_TOP, TYPE_TEESHIRT, TYPE_PANTS, TYPE_BLANKET, TYPE_SERVER, TYPE_EXT_STEP, TYPE_DBG_SHAPE,
	TYPE_POOL_BALL, TYPE_POOL_CUE, TYPE_WALL_MOUNT, TYPE_POOL_TILE, TYPE_POOL_FLOAT, TYPE_BENCH, TYPE_DIV_BOARD, TYPE_FALSE_DOOR, TYPE_FLASHLIGHT, TYPE_CANDLE,
	TYPE_CAMERA, TYPE_CLOCK, TYPE_DOWNSPOUT, TYPE_SHELFRACK, TYPE_CHIM_CAP, TYPE_FOOD_BOX, TYPE_SAFE, TYPE_LADDER, TYPE_CO_COUNTER, TYPE_FISHTANK,
	TYPE_LAVALAMP, TYPE_SHOWERTUB, TYPE_TRASH, TYPE_VALVE, TYPE_METAL_BAR, TYPE_OFF_PILLAR, TYPE_DRINK_CAN, TYPE_CONF_TABLE, TYPE_INT_WINDOW, TYPE_INT_LADDER,
	TYPE_MACHINE, TYPE_BUCKET, TYPE_SPIWEB, TYPE_TREE, TYPE_THEFT_SENS, TYPE_ELEC_WIRE, TYPE_ERASER, TYPE_DWASHER, TYPE_PET_CAGE, TYPE_IBEAM,
	TYPE_CATWALK, TYPE_VANITY, TYPE_CHEM_TANK, TYPE_HVAC_UNIT, TYPE_WARN_LIGHT, TYPE_GAUGE, TYPE_PALLET, TYPE_SHELF_WALL, TYPE_VENDING, TYPE_MED_CAB,
	TYPE_LOCKER, TYPE_TESTTUBE, TYPE_HARDHAT, TYPE_TOPHAT, TYPE_COMP_MOUSE, TYPE_PARK_GATE, TYPE_CONV_BELT, TYPE_JAIL_BARS, TYPE_STICK_NOTE, TYPE_GYM_WEIGHT,
	TYPE_FOOD_TRAY, TYPE_BAR_SOAP, TYPE_COAT_RACK, TYPE_O_SHOWER, TYPE_CARD_DECK, TYPE_CIGARETTE, TYPE_BULLETS, TYPE_CEIL_TILE, TYPE_WALL_GAP, TYPE_MUSHROOM,
	TYPE_SHELL_CASE, TYPE_PAN_SHELF, TYPE_JAR, TYPE_FOOD_TUB, TYPE_COM_FRIDGE,
	/* these next ones are all 3D models - see logic in room_object_t::is_obj_model_type() */
	TYPE_TOILET, TYPE_SINK, TYPE_TUB, TYPE_FRIDGE, TYPE_STOVE, TYPE_TV, TYPE_MONITOR, TYPE_COUCH, TYPE_OFF_CHAIR, TYPE_URINAL,
	TYPE_LAMP, TYPE_WASHER, TYPE_DRYER, TYPE_KEY, TYPE_HANGER, TYPE_CLOTHES, TYPE_FESCAPE, TYPE_WALL_LAMP, TYPE_CUP, TYPE_TOASTER,
	TYPE_HOOD, TYPE_RCHAIR, TYPE_SILVER, TYPE_TOY_MODEL, TYPE_CEIL_FAN, TYPE_FIRE_EXT, TYPE_FOLD_SHIRT, TYPE_PLANT_MODEL, TYPE_POOL_TABLE, TYPE_POOL_LAD,
	TYPE_BAR_STOOL, TYPE_PADLOCK, TYPE_CHECKOUT, TYPE_WFOUNTAIN, TYPE_BANANA, TYPE_BAN_PEEL, TYPE_CONF_PHONE, TYPE_SHOE, TYPE_SHOEBOX, TYPE_VENT_FAN,
	TYPE_HOSP_BED, TYPE_HOSP_CURT, TYPE_FORKLIFT, TYPE_WHEELCHAIR, TYPE_OP_TABLE, TYPE_TROLLEY, TYPE_STRETCHER, TYPE_APPLE, TYPE_EX_MACHINE, TYPE_VIS_PHONE,
	TYPE_JUMPSUIT, TYPE_HANDGUN, TYPE_SHOP_CART, TYPE_CASHREG, TYPE_FOOD_FISH, TYPE_KITCH_APP, TYPE_MILK, TYPE_RADIATOR, TYPE_RAD_FAN, TYPE_SURG_TOOLS,
	TYPE_TOWEL_DISP, TYPE_SOAP_DISP, TYPE_HAND_DRYER,
	/* shared with city objects */
	TYPE_GBIKE, TYPE_XFORMER, TYPE_US_FLAG, TYPE_BLDG_FOUNT,
	/* animals; bird is only used for pet stores */
	TYPE_RAT, TYPE_ROACH, TYPE_SPIDER, TYPE_SNAKE, TYPE_INSECT, TYPE_FISH, TYPE_BIRD,
	NUM_ROBJ_TYPES};
typedef uint8_t room_object;

// room object and stairs shapes
enum {SHAPE_CUBE=0, SHAPE_CYLIN, SHAPE_SPHERE, SHAPE_STAIRS_U, SHAPE_STAIRS_L, SHAPE_STAIRS_FAN, SHAPE_TALL, SHAPE_SHORT, SHAPE_ANGLED, SHAPE_VERT_TORUS, SHAPE_ROUNDED_CUBE};
typedef uint8_t room_obj_shape;

// room types; any new types must also be added to the three tables below
enum {RTYPE_NOTSET=0, RTYPE_HALL, RTYPE_STAIRS, RTYPE_OFFICE, RTYPE_BATH, RTYPE_MENS, RTYPE_WOMENS, RTYPE_BED, RTYPE_KITCHEN, RTYPE_LIVING,
	RTYPE_DINING, RTYPE_STUDY, RTYPE_ENTRY, RTYPE_LIBRARY, RTYPE_STORAGE, RTYPE_GARAGE, RTYPE_SHED, RTYPE_LOBBY, RTYPE_LAUNDRY, RTYPE_CARD,
	RTYPE_PLAY, RTYPE_ART, RTYPE_UTILITY, RTYPE_PARKING, RTYPE_RAMP_EXIT, RTYPE_ATTIC, RTYPE_MASTER_BED, RTYPE_UNFINISHED, RTYPE_SERVER, RTYPE_POOL,
	RTYPE_SWIM, RTYPE_SECURITY, RTYPE_LOUNGE, RTYPE_COMMON, RTYPE_BACKROOMS, RTYPE_RETAIL, RTYPE_ELEVATOR, RTYPE_CONF, RTYPE_MACHINE, RTYPE_INTERR,
	RTYPE_ELEV_EQUIP, RTYPE_STORE, RTYPE_MALL, RTYPE_RESTAURANT, RTYPE_FACTORY, RTYPE_WAREHOUSE, RTYPE_HOS_BED, RTYPE_HOS_OR, RTYPE_HOS_EXAM, RTYPE_CLASS,
	RTYPE_WAITING, RTYPE_LAB, RTYPE_CAFETERIA, RTYPE_LOCKER, RTYPE_JAIL, RTYPE_JAIL_CELL, RTYPE_GYM, RTYPE_VISIT, RTYPE_SHOWER, RTYPE_INFIRMARY,
	RTYPE_CAVE,
	NUM_RTYPES};
typedef uint8_t room_type;

inline bool is_bathroom (room_type   const rtype) {return (rtype == RTYPE_BATH   || rtype == RTYPE_MENS || rtype == RTYPE_WOMENS);}
inline bool is_bedroom  (room_type   const rtype) {return (rtype == RTYPE_BED    || rtype == RTYPE_MASTER_BED);}
inline bool is_jail_room(room_type   const rtype) {return (rtype == RTYPE_JAIL   || rtype == RTYPE_JAIL_CELL);}
inline bool is_ball_type(room_object const type ) {return (type  == TYPE_LG_BALL || type  == TYPE_POOL_BALL );}

// T-shirts are colored, jeans are always white
inline colorRGBA const &gen_teeshirt_color(rand_gen_t &rgen) {return TSHIRT_COLORS[rgen.rand()%NUM_TSHIRT_COLORS];}
inline colorRGBA const &gen_shirt_pants_color(unsigned type, rand_gen_t &rgen) {return ((type == TYPE_TEESHIRT) ? gen_teeshirt_color(rgen) : WHITE);}

// full room names for UI display
std::string const room_names[NUM_RTYPES] =
{"Unset", "Hallway", "Stairs", "Office", "Bathroom", "Men's Restroom", "Women's Restroom", "Bedroom", "Kitchen", "Living Room",
"Dining Room", "Study", "Entryway", "Library", "Storage Room", "Garage", "Shed", "Lobby", "Laundry Room", "Card Room",
"Play Room", "Art Room", "Utility Room", "Parking Garage", "Ramp Exit", "Attic", "Master Bedroom", "Unfinished Room", "Server Room", "Pool Room",
"Swimming Pool Room", "Security Room", "Lounge", "Common Room", "Backrooms", "Retail", "Elevator", "Conference Room", "Machine Room", "Interrogation Room",
"Elev Equip Room", "Store", "Mall Concourse", "Restaurant", "Factory Floor", "Warehouse", "Hospital Bedroom", "Operating Room", "Exam Room", "Classroom",
"Waiting Room", "Laboratory", "Cafeteria", "Locker Room", "Jailroom", "Jail Cell", "Gym", "Visitation", "Shower", "Infirmary"
};
// short room names for breaker panel labels (should be <= 8 characters)
std::string const room_names_short[NUM_RTYPES] =
{"", "Hall", "Stairs", "Office", "Bath", "Men", "Women", "Bed", "Kitchen", "Living",
"Dining", "Study", "Entry", "Library", "Storage", "Garage", "Shed", "Lobby", "Laundry", "Card",
"Play", "Art", "Utility", "Garage", "Ramp", "Attic", "Bed", "", "Server", "Pool",
"Swim", "Security", "Lounge", "Common", "Basement", "Retail", "Elevator", "Conference", "Machine", "Dungeon",
"Equipment", "Store", "Mall", "Restaurant", "Factory", "Warehouse", "Bedroom", "OR", "Exam", "Class",
"Waiting", "Lab", "Cafeteria", "Locker", "Jail", "Cells", "Gym", "Visit", "Shower", "Infirm"
};

unsigned const room_priorities[NUM_RTYPES] = { // for breaker labels; higher values have higher priority
	0, 2, 1, 1, 2, 2, 2, 3, 3, 3,
	3, 1, 3, 3, 2, 3, 3, 2, 2, 2,
	2, 2, 3, 3, 0, 3, 3, 0, 4, 3,
	4, 4, 4, 0, 1, 2, 1, 3, 2, 2,
	1, 3, 4, 3, 1, 1, 2, 3, 3, 2,
	2, 3, 4, 3, 3, 1, 4, 3, 3, 4,
	4
};

// store types, for use with object placement and naming
enum {STORE_OTHER=0, STORE_CLOTHING, STORE_FOOD, STORE_BOOK, STORE_RETAIL, STORE_FURNITURE, STORE_PETS, STORE_APPLIANCE, STORE_SHOE, NUM_STORE_TYPES};
enum {RETAIL_BOXED=0, RETAIL_FOOD, RETAIL_HOUSE_GOODS, RETAIL_KITCHEN, RETAIL_ELECTRONICS, NUM_RETAIL_CAT};
std::string const store_type_strs [NUM_STORE_TYPES] = {"", "clothing", "food", "book", "retail", "furniture", "pet", "appliance", "shoe"};
std::string const srack_categories[NUM_RETAIL_CAT ] = {"boxed items", "food", "household goods", "kitchen", "electronics"};

enum {SHAPE_STRAIGHT=0, SHAPE_U, SHAPE_WALLED, SHAPE_WALLED_SIDES, SHAPE_RAMP, SHAPE_L, SHAPE_FAN, SHAPE_SPIRAL}; // stairs shapes; SHAPE_FAN is used for mall entrances
typedef uint8_t stairs_shape;

enum {ROOM_WALL_INT=0, ROOM_WALL_SEP, ROOM_WALL_EXT, ROOM_WALL_BASEMENT};
enum {FLOORING_MARBLE=0, FLOORING_TILE, FLOORING_CONCRETE, FLOORING_CARPET, FLOORING_WOOD, FLOORING_LGTILE, FLOORING_RUBBER, FLOORING_ROCK, NUM_FLOORING_TYPES}; // not all are used
enum {MAT_TYPE_STATIC=0, MAT_TYPE_SMALL, MAT_TYPE_DYNAMIC, MAT_TYPE_DETAIL, MAT_TYPE_DOORS, MAT_TYPE_LIGHTS, MAT_TYPE_TEXT}; // building_room_geom_t material types; max is 8
enum {FTYPE_NONE=0, FTYPE_BASEMENT, FTYPE_ATTIC}; // for furnace
enum {ATTIC_TYPE_RAFTERS=0, ATTIC_TYPE_FIBERGLASS, ATTIC_TYPE_WOOD, ATTIC_TYPE_PLASTER, NUM_ATTIC_TYPES};
enum {BLDG_COLL_NONE=0, BLDG_COLL_SIDE, BLDG_COLL_ROOF, BLDG_COLL_DETAIL, BLDG_COLL_DRIVEWAY, BLDG_COLL_FENCE, BLDG_COLL_SKYLIGHT};
enum {PLACED_TOILET=1, PLACED_SINK=2, PLACED_TUB=4, PLACED_SHOWER=8}; // for bathroom objects

// Note: when adding an entry here, must also add a string to model_opt_names in city_building_params.cpp
enum {/*building models*/ OBJ_MODEL_TOILET=0, OBJ_MODEL_SINK, OBJ_MODEL_TUB, OBJ_MODEL_FRIDGE, OBJ_MODEL_STOVE, OBJ_MODEL_TV, OBJ_MODEL_MONITOR/*unused*/, OBJ_MODEL_COUCH,
	OBJ_MODEL_OFFICE_CHAIR, OBJ_MODEL_URINAL, OBJ_MODEL_LAMP, OBJ_MODEL_WASHER, OBJ_MODEL_DRYER, OBJ_MODEL_KEY, OBJ_MODEL_HANGER, OBJ_MODEL_CLOTHES, OBJ_MODEL_FESCAPE,
	OBJ_MODEL_WALL_LAMP, OBJ_MODEL_CUP, OBJ_MODEL_TOASTER, OBJ_MODEL_HOOD, OBJ_MODEL_RCHAIR, OBJ_MODEL_SILVER, OBJ_MODEL_TOY, OBJ_MODEL_CEIL_FAN, OBJ_MODEL_FIRE_EXT,
	OBJ_MODEL_FOLD_SHIRT, OBJ_MODEL_PLANT, OBJ_MODEL_POOL_TABLE, OBJ_MODEL_POOL_LAD, OBJ_MODEL_BAR_STOOL, OBJ_MODEL_PADLOCK, OBJ_MODEL_CHECKOUT, OBJ_MODEL_WFOUNTAIN,
	OBJ_MODEL_BANANA, OBJ_MODEL_BAN_PEEL, OBJ_MODEL_PHONE, OBJ_MODEL_SHOE, OBJ_MODEL_SHOEBOX, OBJ_MODEL_VENT_FAN, OBJ_MODEL_HOSP_BED, OBJ_MODEL_HOSP_CURT, OBJ_MODEL_FORKLIFT,
	OBJ_MODEL_WHEELCHAIR, OBJ_MODEL_OP_TABLE, OBJ_MODEL_TROLLEY, OBJ_MODEL_STRETCHER, OBJ_MODEL_APPLE, OBJ_MODEL_EX_MACHINE, OBJ_MODEL_VIS_PHONE, OBJ_MODEL_JUMPSUIT,
	OBJ_MODEL_HANDGUN, OBJ_MODEL_SHOP_CART, OBJ_MODEL_CASHREG, OBJ_MODEL_FISH, OBJ_MODEL_CK_APP, OBJ_MODEL_MILK, OBJ_MODEL_RADIATOR, OBJ_MODEL_RAD_FAN, OBJ_MODEL_SURG_TOOLS,
	OBJ_MODEL_TOWEL_DISP, OBJ_MODEL_SOAP_DISP, OBJ_MODEL_HAND_DRYER,
	OBJ_MODEL_GBIKE/*unused*/, OBJ_MODEL_XFMR/*unused*/, OBJ_MODEL_US_FLAG/*unused*/, OBJ_MODEL_BLDG_FOUNT/*unused*/,
	/*animal models*/ OBJ_MODEL_RAT, OBJ_MODEL_ROACH,
	/*building non-room objects*/ OBJ_MODEL_DOOR_HANDLE,
	/*city models*/ OBJ_MODEL_FHYDRANT, OBJ_MODEL_SUBSTATION, OBJ_MODEL_MAILBOX, OBJ_MODEL_UMBRELLA, OBJ_MODEL_PIGEON, OBJ_MODEL_FOUNTAIN, OBJ_MODEL_BIRD_ANIM, OBJ_MODEL_FLAG,
	OBJ_MODEL_BICYCLE, OBJ_MODEL_SWINGSET, OBJ_MODEL_TRAMPOLINE, OBJ_MODEL_DUMPSTER, OBJ_MODEL_BIG_UMBRELLA, OBJ_MODEL_FLOWER, OBJ_MODEL_DECK_CHAIR, OBJ_MODEL_PICNIC,
	OBJ_MODEL_WIND_TUR, OBJ_MODEL_BB_HOOP, OBJ_MODEL_GAS_PUMP, OBJ_MODEL_STATUE, NUM_OBJ_MODELS};

enum {PART_EFFECT_NONE=0, PART_EFFECT_SPARK, PART_EFFECT_CLOUD, PART_EFFECT_SMOKE, PART_EFFECT_SPLASH, PART_EFFECT_BUBBLE, PART_EFFECT_DROPLET, PART_EFFECT_STEAM,
	PART_EFFECT_FLASH, NUM_PART_EFFECTS};
enum {PIPE_TYPE_SEWER=0, PIPE_TYPE_CW, PIPE_TYPE_HW, PIPE_TYPE_GAS, NUM_PIPE_TYPES};

enum {BIRD_STATE_FLYING=0, BIRD_STATE_GLIDING, BIRD_STATE_LANDING, BIRD_STATE_STANDING, BIRD_STATE_TAKEOFF, NUM_BIRD_STATES};

enum {WOOD_TYPE_DARK=0, WOOD_TYPE_OAK, WOOD_TYPE_PLYWOOD};

struct ball_type_t {
	std::string name, tex_fname, nm_fname;
	float radius, density, value, weight, spec, shine, elastic, friction; // radius in inches, value in dollars, weight in pounds
	bool can_kick, hurts_zombie, breaks_glass;
	ball_type_t(std::string const &name_, std::string const &fn, std::string const &nm, float r, float d, float v, float w, bool ck, bool hz, bool bg,
		float spec_, float shine_, float e, float f) :
		name(name_), tex_fname(fn), nm_fname(nm), radius(r), density(d), value(v), weight(w), spec(spec_), shine(shine_),
		elastic(e), friction(f), can_kick(ck), hurts_zombie(hz), breaks_glass(bg) {}
};
enum {BALL_TYPE_SOCCER=0, BALL_TYPE_BASKET, BALL_TYPE_SOFT, BALL_TYPE_TENNIS, BALL_TYPE_BEACH, NUM_BALL_TYPES};

// name tex_fname nm_fname radius density value weight can_kick hurts_zombie breaks_glass spec shine elastic friction
ball_type_t const ball_types[NUM_BALL_TYPES] = {
	ball_type_t("soccer ball", "balls/soccer_ball_diffuse.png", "balls/soccer_ball_normal.png", 4.4, 0.50, 12.0, 0.90, 1, 1, 1, 0.4, 60.0, 1.0, 1.0),
	ball_type_t("basketball",  "balls/basketball.png",          "",                             4.7, 0.50, 15.0, 1.38, 1, 1, 1, 0.2, 40.0, 1.0, 1.0),
	ball_type_t("softball",    "balls/softball.jpg",            "",                             1.9, 1.20,  5.0, 0.40, 0, 1, 1, 0.1, 20.0, 0.8, 1.0), // balls/softball_bump.jpg    bad format
	ball_type_t("tennis ball", "balls/tennis_ball.jpg",         "",                             1.3, 0.75,  2.0, 0.13, 0, 1, 0, 0.0,  0.0, 1.0, 2.0), // balls/tennis_ball_bump.jpg bad format
	ball_type_t("beach ball",  "balls/beachball.jpg",           "",                            10.0, 0.01, 10.0, 0.10, 1, 0, 0, 0.5, 80.0, 0.8, 1.0)
};
ball_type_t const pool_ball_type("pool ball", "balls/pool_balls.png", "",                     1.125, 1.70,  2.0, 0.37, 0, 1, 1, 0.9, 100.0,0.5, 2.5);

// object flags
unsigned const RO_FLAG_LIT     = 0x01; // light is on
unsigned const RO_FLAG_TOS     = 0x02; // at top of stairs; used for railings and lights; used for objects on open top shelves of shelf racks
unsigned const RO_FLAG_IN_WH   = 0x02; // in warehouse; aliased with RO_FLAG_TOS
unsigned const RO_FLAG_RSTAIRS = 0x04; // in a room with stairs
unsigned const RO_FLAG_INVIS   = 0x08; // invisible
unsigned const RO_FLAG_NOCOLL  = 0x10; // no collision detection
unsigned const RO_FLAG_OPEN    = 0x20; // open, for elevators, closet doors, bathroom stalls, phones, etc.
unsigned const RO_FLAG_NODYNAM = 0x40; // for light shadow maps
unsigned const RO_FLAG_INTERIOR= 0x80; // applies to containing room
// object flags, second byte
unsigned const RO_FLAG_EMISSIVE= 0x0100; // for signs, lights, and phones
unsigned const RO_FLAG_HANGING = 0x0200; // for signs, blinds, hangers, shirts, teeshirts, beams, walls, ladders, balconies, and boxes on mesh; treated as "folding" for closet doors
unsigned const RO_FLAG_ADJ_LO  = 0x0400; // for kitchen counters/closets/door trim/blinds/railings
unsigned const RO_FLAG_ADJ_HI  = 0x0800; // for kitchen counters/closets/door trim/blinds/railings
unsigned const RO_FLAG_ADJ_BOT = 0x1000; // for door trim/railings/ext steps/etc.
unsigned const RO_FLAG_ADJ_TOP = 0x2000; // for door trim/railings
unsigned const RO_FLAG_IS_HOUSE= 0x4000; // used for mirror reflections, shelves, tables, desks, beds, closets, and false doors
unsigned const RO_FLAG_RAND_ROT= 0x8000; // random rotation; used for office chairs, papers, pictures, cups, balls, and ceiling tiles
unsigned const RO_FLAG_UNTEXTURED= 0x1000; // for shirts and desks, aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_FROM_SET  = 0x1000; // for books, aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_HAS_VOL_IX= 0x2000; // for books, aliased with RO_FLAG_ADJ_TOP
unsigned const RO_FLAG_FOR_CAR   = 0x1000; // for car blockers, aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_WALKWAY   = 0x1000; // for walkway objects (outside of buildings), aliased with RO_FLAG_ADJ_BOT
// object flags, third byte, for pickup/interact state
unsigned const RO_FLAG_IN_HALLWAY= 0x010000; // for attic doors, trashcans, and hanging wires
unsigned const RO_FLAG_IN_MALL   = 0x010000; // for mall chairs, tables, benches, trashcans, etc.; aliased with RO_FLAG_IN_HALLWAY and RO_FLAG_IN_FACTORY
unsigned const RO_FLAG_IN_FACTORY= 0x010000; // for machines, ladders, etc.; aliased with RO_FLAG_IN_HALLWAY and RO_FLAG_IN_MALL
unsigned const RO_FLAG_IN_ATTIC  = 0x020000; // in attic
unsigned const RO_FLAG_HAS_EXTRA = 0x040000; // used for counter backsplash, ext wall trim, desks with comp monitors and kbds, books on glass tables, hotel closets, paper towels, teeshirts, pants
unsigned const RO_FLAG_EXTERIOR  = 0x080000; // for signs, window trim, etc.
unsigned const RO_FLAG_EXPANDED  = 0x100000; // for shelves, closets, boxes, and mirrors
unsigned const RO_FLAG_WAS_EXP   = 0x200000; // for objects in/on shelves, closets, drawers, cabinets, shelfracks, and books
unsigned const RO_FLAG_ROTATING  = 0x400000; // for office chairs, clothes on hangers, ceiling fans, and broken lights hanging from and edge
unsigned const RO_FLAG_IN_CLOSET = 0x800000; // for closet lights and light switches, aliased with RO_FLAG_ON_FLOOR/RO_FLAG_ON_SRACK
unsigned const RO_FLAG_ON_SRACK  = 0x800000; // on shelf rack; aliased with RO_FLAG_IN_CLOSET/RO_FLAG_ON_FLOOR
unsigned const RO_FLAG_NONEMPTY  = 0x040000; // for microwaves, shelves, shelfracks, and cups, aliased with RO_FLAG_HAS_EXTRA
unsigned const RO_FLAG_ON_FLOOR  = 0x800000; // for books, fallen objects, upper floor shelfracks, etc., aliased with RO_FLAG_IN_CLOSET/RO_FLAG_ON_SRACK
unsigned const RO_FLAG_BROKEN2   = 0x040000; // for lights that are completely broken, aliased with RO_FLAG_HAS_EXTRA and RO_FLAG_NONEMPTY
unsigned const RO_FLAG_PLCOLL    = 0x040000; // player collidable, for chairs, aliased with RO_FLAG_HAS_EXTRA
unsigned const RO_FLAG_IN_JAIL   = 0x040000; // for jail cell/prison objects such as bed and stall/shower; aliased with RO_FLAG_HAS_EXTRA
// object flags, fourth byte
unsigned const RO_FLAG_DYNAMIC  = 0x01000000; // dynamic object (balls, elevators, etc.)
unsigned const RO_FLAG_DSTATE   = 0x02000000; // this object has dynamic state
unsigned const RO_FLAG_NO_CONS  = 0x04000000; // this object is not consumable (bottles)
unsigned const RO_FLAG_NO_POWER = 0x04000000; // unpowered; related to circuit breakers, aliased with RO_FLAG_NO_CONS
unsigned const RO_FLAG_IS_ACTIVE= 0x08000000; // active, for sinks, tubs, buttons, pool balls, shower curtains, etc.
unsigned const RO_FLAG_USED     = 0x10000000; // used by the player (spraypaint/marker), parking spaces to indicate cars, school books, occupied beds, and unopenable boxes
unsigned const RO_FLAG_IN_ELEV  = 0x20000000; // for elevator lights, buttons, and flooring; aliased with RO_FLAG_BACKROOM and RO_FLAG_IN_POOL
unsigned const RO_FLAG_BACKROOM = 0x20000000; // in backrooms, for walls and pillars; aliased with RO_FLAG_IN_ELEV and RO_FLAG_IN_POOL
unsigned const RO_FLAG_IN_POOL  = 0x20000000; // for stairs, railings, and drains; aliased with RO_FLAG_IN_ELEV and RO_FLAG_BACKROOM
unsigned const RO_FLAG_BROKEN   = 0x40000000; // for TVs, monitors, flickering lights, and computers; used for stained trays/plates; maybe can use for windows
unsigned const RO_FLAG_MOVED    = 0x80000000; // for player push/pull

// reflection pass flags
unsigned const REF_PASS_ENABLED     = 0x01;
unsigned const REF_PASS_INTERIOR    = 0x02;
unsigned const REF_PASS_CUBE_MAP    = 0x04;
unsigned const REF_PASS_WATER       = 0x08;
unsigned const REF_PASS_EXTB        = 0x10;
unsigned const REF_PASS_NO_MIRROR   = 0x20;
unsigned const REF_PASS_INT_ONLY    = 0x40;
unsigned const REF_PASS_GLASS_FLOOR = 0x80;
unsigned const REF_PASS_CITY_ONLY   = 0x100;
unsigned const REF_PASS_EXT_ONLY    = 0x200;

typedef vector<point> vect_point;
typedef vector<sphere_t> vect_sphere_t;

inline point get_camera_building_space() {return (get_camera_pos() - get_tiled_terrain_model_xlate());}
inline void set_cube_zvals(cube_t &c, float z1, float z2) {c.z1() = z1; c.z2() = z2;}
inline void copy_zvals(cube_t &to, cube_t const &from) {to.z1() = from.z1(); to.z2() = from.z2();}
inline float get_tc_leg_width(cube_t const &c, float width) {return min(0.5f*width*(c.dx() + c.dy()), 0.2f*c.dz());} // make legs square
inline unsigned get_rgeom_sphere_ndiv(bool low_detail) {return (low_detail ? N_SPHERE_DIV/2 : N_SPHERE_DIV);}
inline point cube_bot_center(cube_t const &c) {return point(c.xc(), c.yc(), c.z1());}
inline point cube_top_center(cube_t const &c) {return point(c.xc(), c.yc(), c.z2());}
inline point get_cube_center_zval(cube_t const &c, float zval) {return point(c.xc(), c.yc(), zval);}
inline bool is_known_metal_color(colorRGBA const &c) {return (c == COPPER_C || c == BRASS_C || c == DARK_BRASS_C || c == BRONZE_C || c == GOLD);}
inline colorRGBA get_specular_color(colorRGBA const &c) {return (is_known_metal_color(c) ? c : WHITE);}
inline colorRGBA gen_box_color(rand_gen_t &rgen) {return colorRGBA(rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0));} // add minor color variation


struct oriented_cube_t : public cube_t {
	bool dim=0, dir=0;
	oriented_cube_t() {}
	oriented_cube_t(cube_t const &c, bool dim_, bool dir_) : cube_t(c), dim(dim_), dir(dir_) {}
	float get_length() const {return get_sz_dim( dim);}
	float get_width () const {return get_sz_dim(!dim);}
	float get_height() const {return dz();}
};

struct building_occlusion_state_t {
	int exclude_bix=-1;
	bool skip_cont_camera=0;
	point pos;
	vector3d xlate;
	vector<cube_with_ix_t> building_ids;

	void init(point const &pos_, vector3d const &xlate_) {
		pos   = pos_;
		xlate = xlate_;
		building_ids.clear();
	}
};

class occlusion_checker_t {
	building_occlusion_state_t state;
public:
	vect_cube_t occluders;
	void set_exclude_bix(int exclude_bix) {state.exclude_bix = exclude_bix;}
	void set_exclude_camera_building() {state.skip_cont_camera = 1;}
	void set_camera(pos_dir_up const &pdu);
	bool is_occluded(cube_t const &c) const;
};

struct building_t;
class building_creator_t;

class occlusion_checker_noncity_t {
	building_occlusion_state_t state;
	building_creator_t const &bc;
public:
	bool query_is_for_light, for_shadows, extra_occ_dim=0;
	vect_cube_t const *extra_occluders=nullptr; // for mall stores

	occlusion_checker_noncity_t(building_creator_t const &bc_, bool for_light=0, bool for_shadows_=0) : bc(bc_), query_is_for_light(for_light), for_shadows(for_shadows_) {}
	void set_exclude_bix(int exclude_bix) {state.exclude_bix = exclude_bix;}
	void set_camera(pos_dir_up const &pdu, bool cur_building_only=0);
	bool is_occluded(cube_t const &c) const;
	bool check_custom_occluder_cull(cube_t const &c, point const &viewer) const;
	vector3d const &get_xlate() const {return state.xlate;}
};

struct ped_draw_vars_t {
	building_t &building; // Note: building ref is non-const because we may update rendering data for the people inside it
	occlusion_checker_noncity_t &oc;
	shader_t &s;
	vector3d const &xlate;
	unsigned bix;
	bool shadow_only, reflection_pass, in_retail_room;

	ped_draw_vars_t(building_t &b, occlusion_checker_noncity_t &oc_, shader_t &s_, vector3d const &x, unsigned bix_, bool so, bool rp, bool irr=0)
		: building(b), oc(oc_), s(s_), xlate(x), bix(bix_), shadow_only(so), reflection_pass(rp), in_retail_room(irr) {}
};

struct city_zone_t : public cube_t {
	float zval=0.0;
	bool is_park=0, is_residential=0;
	uint8_t street_dir=0; // encoded as 2*dim + dir + 1; 0 is unassigned
	unsigned nbuildings=0, capacity=0; // in number of buildings; 0 is unlimited
	unsigned max_floors=0; // 0=unlimited
	unsigned city_ix=0, street_num=0;
	int parent_plot_ix=-1; // if this is a sub-plot; -1 otherwise
	std::string address; // or just store road_name?

	city_zone_t() {}
	city_zone_t(cube_t const &c, float zval_, bool p, bool r, unsigned sdir, unsigned cap, int ppix, int cix, unsigned mf) :
		cube_t(c), zval(zval_), is_park(p), is_residential(r), street_dir(sdir), capacity(cap), max_floors(mf), city_ix(cix), parent_plot_ix(ppix) {}
	bool is_full() const {return (capacity > 0 && nbuildings >= capacity);}
};

typedef vector<city_zone_t> vect_city_zone_t;

struct tid_nm_pair_dstate_t {
	shader_t &s;
	depth_write_tracker_t dwt;
	int bmm_loc=-1;
	float bump_map_mag=1.0, crack_weight=0.0;
	bool no_set_texture=0;

	tid_nm_pair_dstate_t(shader_t &s_, bool no_set_texture_=0, float crack_weight_=0.0) : s(s_), crack_weight(crack_weight_), no_set_texture(no_set_texture_) {}
	void set_for_shader(float new_bump_map_mag);
	~tid_nm_pair_dstate_t();
};

struct tid_nm_pair_t { // size=48

	int tid=-1, nm_tid=-1; // Note: assumes each tid has only one nm_tid
	float tscale_x=1.0, tscale_y=1.0, txoff=0.0, tyoff=0.0, emissive=0.0, metalness=0.0, refract_ix=1.0;
	color_wrapper spec_color;
	unsigned char shininess=0; // Note: spec_mag is divided by 255.0
	bool shadowed   =0; // Note: doesn't directly affect rendering, only used for uniquing/operator==()
	bool shadow_only=0; // only drawn in the shadow pass; used for simplified shadow casters
	bool transparent=0; // used to draw batched alpha blended materials last
	bool no_cracks  =0; // for basement crack effects
	bool no_reflect =0; // no cube map reflections (optimization)

	tid_nm_pair_t() {}
	tid_nm_pair_t(int tid_, float txy=1.0, bool shadowed_=0, bool transparent_=0, bool no_reflect_=0) : tid(tid_), nm_tid(FLAT_NMAP_TEX),
		tscale_x(txy), tscale_y(txy), shadowed(shadowed_), transparent(transparent_), no_reflect(no_reflect_) {} // non-normal mapped 1:1 texture AR
	tid_nm_pair_t(int tid_, int nm_tid_, float tx, float ty, float xo=0.0, float yo=0.0, bool shadowed_=0, bool transparent_=0, bool no_reflect_=0) :
		tid(tid_), nm_tid(nm_tid_), tscale_x(tx), tscale_y(ty), txoff(xo), tyoff(yo), shadowed(shadowed_), transparent(transparent_), no_reflect(no_reflect_) {}
	tid_nm_pair_t(tid_nm_pair_t const &t, bool shadowed_) : tid_nm_pair_t(t) {shadowed = shadowed_;} // constructor to change the shadowed flag
	void set_shininess(float shine) {shininess = (unsigned char)max(1, min(255, round_fp(shine)));}
	void set_specular(float mag, float shine) {set_specular_color(WHITE, mag, shine);}
	void set_specular(float mag, float shine, float metalness_) {set_specular(mag, shine); metalness = metalness_;}
	void set_specular_color(colorRGB const &color, float mag, float shine);
	void set_metal_specular(colorRGB const &color=WHITE, float mag=0.8, float shine=60.0) {set_specular_color(color, mag, shine);}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}

	bool is_compat_ignore_shadowed(tid_nm_pair_t const &t) const {
		return (tid == t.tid && nm_tid == t.nm_tid && emissive == t.emissive && metalness == t.metalness && refract_ix == t.refract_ix &&
			shininess == t.shininess && transparent == t.transparent && spec_color == t.spec_color && no_reflect == t.no_reflect);
	}
	bool is_compatible(tid_nm_pair_t const &t) const {return (is_compat_ignore_shadowed(t) && shadowed == t.shadowed && shadow_only == t.shadow_only);}
	bool operator==(tid_nm_pair_t const &t) const {return (is_compatible(t) && tscale_x == t.tscale_x && tscale_y == t.tscale_y && txoff == t.txoff && tyoff == t.tyoff);}
	bool operator!=(tid_nm_pair_t const &t) const {return !operator==(t);}
	int get_nm_tid() const {return ((nm_tid < 0) ? FLAT_NMAP_TEX : nm_tid);}
	float get_drawn_tscale_x() const {return 2.0f*tscale_x;} // adjust for local vs. global space change
	float get_drawn_tscale_y() const {return 2.0f*tscale_y;} // adjust for local vs. global space change
	float get_emissive_val  () const;
	colorRGBA get_avg_color () const {return texture_color(tid);}
	tid_nm_pair_t get_scaled_version(float scale) const;
	static bool bind_reflection_shader();
	void   set_gl(tid_nm_pair_dstate_t &state) const;
	void unset_gl(tid_nm_pair_dstate_t &state) const;
	void toggle_transparent_windows_mode();
};

struct building_tex_params_t {
	tid_nm_pair_t side_tex, roof_tex; // exterior
	tid_nm_pair_t wall_tex, ceil_tex, floor_tex, house_ceil_tex, house_floor_tex, basement_floor_tex; // interior

	bool has_normal_map() const {return (side_tex.nm_tid >= 0 || roof_tex.nm_tid >= 0 || wall_tex.nm_tid >= 0 || ceil_tex.nm_tid >= 0 ||
		floor_tex.nm_tid >= 0 || house_ceil_tex.nm_tid >= 0 || house_floor_tex.nm_tid >= 0 || basement_floor_tex.nm_tid >= 0);}
};

struct ug_elev_info_t {
	cube_t entrance;
	float top_floor_z2;
	bool dim, dir;
	ug_elev_info_t(cube_t const &e, float tfz2, bool dim_, bool dir_) : entrance(e), top_floor_z2(tfz2), dim(dim_), dir(dir_) {}
};
typedef vector<ug_elev_info_t> vect_ug_elev_info_t;

struct walkway_base_t {
	bool open_ends[2]={};
	float floor_spacing;
	unsigned side_mat_ix, roof_mat_ix; // matches building material
	colorRGBA side_color, roof_color;
	walkway_base_t(unsigned smix, unsigned rmix, colorRGBA const &sc, colorRGBA const &rc, float floor_spacing_) :
		floor_spacing(floor_spacing_), side_mat_ix(smix), roof_mat_ix(rmix), side_color(sc), roof_color(rc) {}
};
struct bldg_walkway_t : public cube_t, public walkway_base_t {
	bool dim;
	bldg_walkway_t(cube_t const &c, bool d, unsigned smix, unsigned rmix, colorRGBA const &sc, colorRGBA const &rc, float fs) :
		cube_t(c), walkway_base_t(smix, rmix, sc, rc, fs), dim(d) {}
};
typedef vector<bldg_walkway_t> vect_bldg_walkway_t;

struct building_walkway_geom_t {
	cube_t bcube;
	bool dim;
	float door_bounds[2][2]={}; // {dir=0, dir=1} x {lo, hi} for exterior doors, one pair per end
	uint8_t has_door=0; // for false doors; one bit per floor

	building_walkway_geom_t(cube_t const &c, bool dim_) : bcube(c), dim(dim_) {}
	bool has_ext_door(bool dir) const {return (door_bounds[dir][0] < door_bounds[dir][1]);}
	float get_length() const {return bcube.get_sz_dim(dim);}
};
struct building_walkway_t : public building_walkway_geom_t { // "owned" walkway, one per connected building
	bool is_owner, open_ends[2]={};
	building_t *conn_bldg;
	vect_cube_with_ix_t windows;
	vect_cube_t frames;
	cube_t bcube_inc_rooms, skyway_conn, elevator_bcube, elevator_cut;

	building_walkway_t(building_walkway_geom_t const &g, bool owner, building_t *b) : building_walkway_geom_t(g), is_owner(owner), conn_bldg(b) {bcube_inc_rooms = bcube;}
	cube_t get_bcube_inc_open_door() const;
	bool has_skyway_conn() const {return !skyway_conn.is_all_zeros();}
	void attach_elevator(cube_t const &e);
};

struct skyway_conn_t : public cube_t {
	bool dim, dir;
	building_t const *const building;
	skyway_conn_t(cube_t const &c, bool dim_, bool dir_, building_t const *b) : cube_t(c), dim(dim_), dir(dir_), building(b) {}
	string get_building_name() const;
	float  get_doorway_width() const;
};

class city_lights_manager_t {
protected:
	cube_t lights_bcube;
	float light_radius_scale=1.0, dlight_add_thresh=0.0;
	bool prev_had_lights=0;
public:
	virtual ~city_lights_manager_t() {}
	cube_t get_lights_bcube() const {return lights_bcube;}
	void add_player_flashlight(float radius_scale);
	void tighten_light_bcube_bounds(vector<light_source> const &lights);
	void clamp_to_max_lights(vector3d const &xlate, vector<light_source> &lights);
	bool begin_lights_setup(vector3d const &xlate, float light_radius, vector<light_source> &lights);
	void finalize_lights(vector<light_source> &lights);
	void setup_shadow_maps(vector<light_source> &light_sources, point const &cpos, unsigned max_smaps, bool sec_camera_mode=0);
	virtual bool enable_lights() const = 0;
};

struct color_range_t {
	float grayscale_rand=0.0;
	colorRGBA cmin=WHITE, cmax=WHITE; // alpha is unused?
	void gen_color(colorRGBA &color, rand_gen_t &rgen) const;
};

struct building_mat_t : public building_tex_params_t {

	bool no_city=0, add_windows=0, add_wind_lights=0, no_walkways=0;
	unsigned min_levels=1, max_levels=1, min_sides=4, max_sides=4;
	float place_radius=0.0, max_delta_z=0.0, max_rot_angle=0.0, min_level_height=0.0, min_alt=-1000, max_alt=1000, house_prob=0.0, house_scale_min=1.0, house_scale_max=1.0;
	float split_prob=0.0, cube_prob=1.0, round_prob=0.0, asf_prob=0.0, min_fsa=0.0, max_fsa=0.0, min_asf=0.0, max_asf=0.0;
	float wind_xscale=1.0, wind_yscale=1.0, wind_xoff=0.0, wind_yoff=0.0;
	float floor_spacing=0.0, floorplan_wind_xscale=0.0; // these are derived values
	float apartment_prob=0.0;
	cube_t pos_range, prev_pos_range, sz_range; // pos_range z is unused?
	color_range_t side_color, roof_color; // exterior
	colorRGBA window_color=GRAY, wall_color=WHITE, ceil_color=WHITE, floor_color=LT_GRAY, house_ceil_color=WHITE, house_floor_color=WHITE;

	building_mat_t() : pos_range(-100,100,-100,100,0,0), sz_range(1,1,1,1,1,1) {}
	float gen_house_size_scale(rand_gen_t &rgen) const {return ((house_scale_min == house_scale_max) ? house_scale_min : rgen.rand_uniform(house_scale_min, house_scale_max));}
	void update_range(vector3d const &range_translate);
	void set_pos_range(cube_t const &new_pos_range) {prev_pos_range = pos_range; pos_range = new_pos_range;}
	void restore_prev_pos_range() {
		if (!prev_pos_range.is_all_zeros()) {pos_range = prev_pos_range;}
	}
	void finalize();
	float get_window_tx() const;
	float get_window_ty() const;
};

struct building_params_t {

	bool flatten_mesh=0, has_normal_map=0, tex_mirror=0, tex_inv_y=0, tt_only=0, infinite_buildings=0, dome_roof=0, onion_roof=0;
	bool gen_building_interiors=1, add_city_interiors=0, enable_rotated_room_geom=0, add_secondary_buildings=0, add_office_basements=0, add_office_br_basements=0;
	bool put_doors_in_corners=0, cities_all_bldg_mats=0, small_city_buildings=0, add_door_handles=0, use_voronoise_cracks=0, add_basement_tunnels=0, no_retail_and_mall=0;
	unsigned num_place=0, num_tries=10, cur_prob=1, max_shadow_maps=32, buildings_rand_seed=0, max_ext_basement_hall_branches=4, max_ext_basement_room_depth=4;
	unsigned max_room_geom_gen_per_frame=1, max_office_basement_floors=2, max_mall_levels=2;
	float ao_factor=0.0, sec_extra_spacing=0.0, player_coll_radius_scale=1.0, interior_view_dist_scale=1.0;
	float window_width=0.0, window_height=0.0, window_xspace=0.0, window_yspace=0.0; // windows
	float wall_split_thresh=4.0, max_fp_wind_xscale=0.0, max_fp_wind_yscale=0.0, basement_water_level_min=0.0, basement_water_level_max=0.0; // interiors
	float open_door_prob=1.0, locked_door_prob=0.0, basement_prob_house=0.5, basement_prob_office=0.5, ball_prob=0.3, two_floor_retail_prob=0.0; // interior probabilities
	float split_stack_floorplan_prob=0.0, retail_floorplan_prob=0.0, mall_prob=0.5; // floorplan probabilities
	float glass_floor_alpha=GLASS_COLOR.A, max_altitude_all=0.0;
	// consistency probabilities of houses for cities and blocks
	float house_same_mat_prob =0.0, house_same_size_prob =0.0, house_same_geom_prob =0.0, house_same_per_city_prob =0.0;
	float office_same_mat_prob=0.0, office_same_size_prob=0.0, office_same_geom_prob=0.0, office_same_per_city_prob=0.0;
	// building people/AI params
	bool enable_people_ai=0, ai_target_player=1, ai_follow_player=0, allow_elevator_line=1, no_coll_enter_exit_elevator=1, show_player_model=0;
	unsigned ai_opens_doors=1; // 0=don't open doors, 1=only open if player closed door after path selection; 2=always open doors
	unsigned ai_player_vis_test=0; // 0=no test, 1=LOS, 2=LOS+FOV, 3=LOS+FOV+lit
	unsigned ai_sees_player_hide=2; // 0=doesn't see the player, 1=sees the player and waits outside the hiding spot, 2=opens the door and comes in
	unsigned people_per_office_min=0, people_per_office_max=0, people_per_house_min=0, people_per_house_max=0, elevator_capacity=1;
	unsigned player_model_ix=0;
	float ai_retreat_time=4.0, elevator_wait_time=60.0, use_elevator_prob=0.25, elevator_wait_recall_prob=0.5;
	float people_min_alpha=0.0, zombie_fov=0.0;
	// building animal params
	unsigned num_rats_min=0, num_rats_max=0, min_attack_rats=0, num_spiders_min=0, num_spiders_max=0, num_snakes_min=0, num_snakes_max=0, num_insects_min=0, num_insects_max=0;
	float rat_speed   =0.0, rat_size_min   =0.5, rat_size_max   =1.0; // rats
	float spider_speed=0.0, spider_size_min=0.5, spider_size_max=1.0, spider_drawer_prob=0.0; // spiders
	float snake_speed =0.0, snake_size_min =0.5, snake_size_max =1.0; // snakes
	float insect_speed=0.0, insect_size_min=0.5, insect_size_max=1.0; // snakes
	// gameplay state
	float player_weight_limit=100.0;
	// materials
	vector3d range_translate; // used as a temporary to add to material pos_range
	building_mat_t cur_mat;
	vector<building_mat_t> materials;
	vector<unsigned> mat_gen_ix, mat_gen_ix_city, mat_gen_ix_nocity, mat_gen_ix_res; // {any, city_only, non_city, residential}
	vector<unsigned> rug_tids, picture_tids, desktop_tids, sheet_tids, paper_tids, food_box_tids, flag_tids, metal_tids;
	vector<std::string> food_box_names; // same size as food_box_tids
	map<unsigned, unsigned> tid_to_nmap_tid;
	int last_read_tid=-1;
	// use for option reading
	int read_error=0;
	kw_to_val_map_t<bool     >  kwmb;
	kw_to_val_map_t<unsigned >  kwmu;
	kw_to_val_map_t<float    >  kwmf;
	kw_to_val_map_t<colorRGBA>  kwmc;
	kw_to_val_map_float_check_t kwmr;

	building_params_t(unsigned num=0) : num_place(num),
		kwmb(read_error, "buildings"), kwmu(read_error, "buildings"), kwmf(read_error, "buildings"), kwmc(read_error, "buildings"), kwmr(read_error, "buildings") {init_kw_maps();}
	bool parse_buildings_option(FILE *fp);
	int get_wrap_mir() const {return (tex_mirror ? 2 : 1);}
	bool building_people_enabled() const {return (people_per_office_max > 0 || people_per_house_max > 0);}
	bool windows_enabled  () const {return (window_width > 0.0 && window_height > 0.0 && window_xspace > 0.0 && window_yspace);} // all must be specified as nonzero
	bool gen_inf_buildings() const {return (infinite_buildings && world_mode == WMODE_INF_TERRAIN);}
	float get_window_width_fract () const {assert(windows_enabled()); return window_width /(window_width  + window_xspace);}
	float get_window_height_fract() const {assert(windows_enabled()); return window_height/(window_height + window_yspace);}
	float get_window_tx() const {assert(windows_enabled()); return 1.0f/(window_width  + window_xspace);}
	float get_window_ty() const {assert(windows_enabled()); return 1.0f/(window_height + window_yspace);}
	void add_cur_mat();
	void finalize();

	building_mat_t const &get_material(unsigned mat_ix) const {
		assert(mat_ix < materials.size());
		return materials[mat_ix];
	}
	vector<unsigned> const &get_mat_list(bool city_only, bool non_city_only, bool residential) const;
	unsigned choose_rand_mat(rand_gen_t &rgen, bool city_only, bool non_city_only, bool residential) const;
	float get_max_house_size() const;
	void set_pos_range(cube_t const &pos_range);
	void restore_prev_pos_range();
	int get_nm_tid_for(unsigned tid) const;
private:
	void init_kw_maps();
	int read_building_texture(FILE *fp, std::string const &str, bool is_normal_map, int &error, bool check_filename=0, bool *no_cracks=nullptr);
	void read_texture_and_add_if_valid(FILE *fp, std::string const &str, int &error, vector<unsigned> &tids);
};

inline uint64_t get_tile_id_for_cube(cube_t const &c) {return get_tile_id_containing_point_no_xyoff(c.get_cube_center());}

struct cmp_by_tile { // not the most efficient solution, but no memory overhead
	bool operator()(cube_t const &a, cube_t const &b) const {return (get_tile_id_for_cube(a) < get_tile_id_for_cube(b));}
};


template<typename T> static bool check_vect_cube_contains_pt(vector<T> const &cubes, point const &pos);

inline void copy_dim(cube_t &c1, cube_t const &c2, unsigned dim) {
	for (unsigned d = 0; d < 2; ++d) {c1.d[dim][d] = c2.d[dim][d];}
}
inline void intersect_dim(cube_t &c1, cube_t const &c2, unsigned dim) {
	max_eq(c1.d[dim][0], c2.d[dim][0]);
	min_eq(c1.d[dim][1], c2.d[dim][1]);
}
inline void union_dim(cube_t &c1, cube_t const &c2, unsigned dim) {
	min_eq(c1.d[dim][0], c2.d[dim][0]);
	max_eq(c1.d[dim][1], c2.d[dim][1]);
}
inline void clip_low_high_tc(float &t0, float &t1) {
	if (fabs(t0 - t1) < 0.5) {t0 = t1 = 0.0;} // too small to have a window
	else {t0 = round_fp(t0); t1 = round_fp(t1);} // Note: round() is much faster than nearbyint(), and round_fp() is faster than round()
}
template<typename T> cube_t get_cube_height_radius(point const &center, T radius, float height) { // T can be float or vector3d
	cube_t c(center);
	c.expand_by_xy(radius);
	c.z2() += height;
	return c;
}
inline bool check_bcube_sphere_coll(cube_t const &bcube, point const &sc, float radius, bool xy_only) {
	return (xy_only ? sphere_cube_intersect_xy(sc, radius, bcube) : sphere_cube_intersect(sc, radius, bcube));
}

template<typename T> bool has_bcube_int(cube_t const &bcube, vector<T> const &cubes) { // T must derive from cube_t
	for (cube_t const &c : cubes) {if (c.intersects(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int(cube_t const &bcube, vector<T> const &cubes, unsigned start_ix) { // T must derive from cube_t
	assert(start_ix <= cubes.size());
	for (auto i = cubes.begin()+start_ix; i < cubes.end(); ++i) {if (i->intersects(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_no_adj(cube_t const &bcube, vector<T> const &cubes) { // T must derive from cube_t
	for (cube_t const &c : cubes) {if (c.intersects_no_adj(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_xy(cube_t const &bcube, vector<T> const &cubes, float pad_dist=0.0) { // T must derive from cube_t
	cube_t tc(bcube);
	tc.expand_by_xy(pad_dist);
	for (cube_t const &c : cubes) {if (c.intersects_xy(tc)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_xy_no_adj(cube_t const &bcube, vector<T> const &cubes) { // T must derive from cube_t
	for (cube_t const &c : cubes) {if (c.intersects_xy_no_adj(bcube)) return 1;}
	return 0;
}
template<typename T> static bool check_vect_cube_contains_pt_xy(vector<T> const &cubes, point const &pos) {
	for (cube_t const &c : cubes) {if (c.contains_pt_xy(pos)) return 1;}
	return 0;
}
inline bool point_in_cubes_xy_exp(vect_cube_t const &cubes, point const &p, float dist) {
	for (cube_t const &c : cubes) {if (c.contains_pt_xy_exp(p, dist)) return 1;}
	return 0;
}
template<typename T> static bool check_vect_cube_contains_pt(vector<T> const &cubes, point const &pos) {
	for (cube_t const &c : cubes) {if (c.contains_pt(pos)) return 1;}
	return 0;
}
template<typename T> cube_t get_bcubes_union(vector<T> const &cubes) {
	if (cubes.size() == 1) {return cubes.front();}
	cube_t bcube;
	for (cube_t const &c : cubes) {bcube.assign_or_union_with_cube(c);}
	return bcube;
}
inline void swap_cube_dims(cube_t &c, unsigned d1, unsigned d2) {
	for (unsigned d = 0; d < 2; ++d) {swap(c.d[d1][d], c.d[d2][d]);}
}
inline void expand_cube_zvals(cube_t &c, float z1, float z2) {
	min_eq(c.z1(), z1);
	max_eq(c.z2(), z2);
}
inline void set_wall_width(cube_t &wall, float pos, float half_thick, unsigned dim) {
	wall.d[dim][0] = pos - half_thick;
	wall.d[dim][1] = pos + half_thick;
}
template<typename T> void vector_random_shuffle(vector<T> &v, rand_gen_t &rgen) {std::shuffle(v.begin(), v.end(), rand_gen_wrap_t(rgen));}

struct cube_by_sz { // sort cube by size in dim descending
	bool dim;
	cube_by_sz(bool dim_) : dim(dim_) {}
	bool operator()(cube_t const &a, cube_t const &b) const {return (b.get_sz_dim(dim) < a.get_sz_dim(dim));}
};

template<typename T> void add_to_and_clear(T &src, T &dest) {
	vector_add_to(src, dest);
	src.clear();
}
template<typename T> void add_inverted_triangles(T &verts, vector<unsigned> &indices, unsigned verts_start, unsigned ixs_start, bool replace_mode=0);
template<typename T> void reserve_extra(vector<T> &v, unsigned num) {v.reserve(v.size() + num);}

template<typename T> void subtract_cube_from_cube(T const &c, cube_t const &s, vector<T> &out, bool clear_out=0);
template<typename T> void subtract_cube_from_cube_inplace(cube_t const &s, vector<T> &cubes, unsigned &ix, unsigned &iter_end);
template<typename T> void subtract_cubes_from_cube(cube_t const &c, vector<T> const &sub, vect_cube_t &out, vect_cube_t &out2, int zval_mode=0);
template<typename T> bool subtract_cube_from_cubes(cube_t const &s, vector<T> &cubes, vect_cube_t *holes=nullptr, bool clip_in_z=0, bool include_adj=0, bool no_z_test=0);
template<typename T> bool line_int_cubes(point const &p1, point const &p2, vector<T> const &cubes, cube_t const &line_bcube);
void subtract_cube_from_cube_transpose(cube_t c, cube_t s, vect_cube_t &out);

template<typename T> void subtract_cubes_from_cubes(T const &sub, vect_cube_t &cubes, bool clip_in_z=0) {
	for (auto const &s : sub) {subtract_cube_from_cubes(s, cubes, nullptr, clip_in_z);} // no holes
}
void get_city_plot_zones(vect_city_zone_t &zones);
void get_city_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state);
bool check_city_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t const &state);
bool city_single_cube_visible_check(point const &pos, cube_t const &c);
cube_t get_building_lights_bcube();
cube_t get_grid_bcube_for_building(building_t const &b);
bool get_building_door_pos_closest_to(unsigned building_id, point const &target_pos, point &door_pos, bool inc_garage_door=0, int mf_pref=2);
colorRGBA get_light_color_temp(float t);
vector3d get_nom_car_size();
float get_scaled_player_radius();
float get_bldg_player_height();
float get_player_eye_height();
bool check_city_tline_cube_intersect_xy(cube_t const &c);
void get_closest_dim_dir_xy(cube_t const &inner, cube_t const &outer, bool &dim, bool &dir);
template<typename T> void add_sign_text_verts(std::string const &text, cube_t const &sign, bool dim, bool dir, colorRGBA const &color,
	vector<T> &verts_out, float first_char_clip_val=0.0, float last_char_clip_val=0.0, bool include_space_chars=0, bool invert_z=0);
void enable_animations_for_shader(shader_t &s);
int get_normal_map_for_bldg_tid(int tid);
void apply_building_fall_damage(float delta_z);
unsigned get_face_mask(unsigned dim, bool dir);
unsigned get_skip_mask_for_xy(bool dim);
// functions in building_interact.cc and building_gameplay.cc
void gen_sound_thread_safe(unsigned id, point const &pos, float gain=1.0, float pitch=1.0, float gain_scale=1.0, bool skip_if_already_playing=0);
inline void gen_sound_thread_safe_at_player(unsigned id, float gain=1.0, float pitch=1.0, bool skip_if_already_playing=0) {
	gen_sound_thread_safe(id, get_camera_pos(), gain, pitch, 1.0, skip_if_already_playing);
}
bool in_building_gameplay_mode();
int get_rect_panel_tid();
int get_bath_wind_tid ();
int get_int_door_tid  ();
int get_off_door_tid  ();
int get_bldg_door_tid ();
int get_concrete_tid  ();
int get_plywood_tid   ();
int get_insulation_tid();
int get_met_plate_tid ();
int get_mplate_nm_tid ();
// functions in building_room_obj_expand.cc
point gen_xy_pos_in_area(cube_t const &S, vector3d const &sz, rand_gen_t &rgen, float zval=0.0);
point gen_xy_pos_in_area(cube_t const &S, float radius, rand_gen_t &rgen, float zval=0.0);
void gen_xy_pos_in_cube (point &pos, cube_t const &c, rand_gen_t &rgen);
void gen_xyz_pos_in_cube(point &pos, cube_t const &c, rand_gen_t &rgen);
void gen_xy_pos_for_cube_obj(cube_t &C, cube_t const &S, vector3d const &sz, float height, rand_gen_t &rgen, bool place_at_z1=0);
void gen_xy_pos_for_round_obj(cube_t &C, cube_t const &S, float radius, float height, float spacing, rand_gen_t &rgen, bool place_at_z1=0);
cube_t place_cylin_object(rand_gen_t &rgen, cube_t const &place_on, float radius, float height, float dist_from_edge,  bool place_at_z1=0);
// from Universe_name.cpp
std::string gen_random_name(rand_gen_t &rgen, unsigned min_len=0, bool for_universe=0);

