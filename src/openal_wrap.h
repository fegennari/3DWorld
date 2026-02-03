// 3D World - OpenAL Interface Code Header
// by Frank Gennari
// 9/1/11
#include "3DWorld.h"
#include "trigger.h"

#pragma once

// static sounds
enum {SOUND_BURNING=0, SOUND_RAIN1, SOUND_WIND1, SOUND_UNDERWATER, SOUND_EXPLODE, SOUND_GUNSHOT, SOUND_SHOTGUN, SOUND_FIREBALL, SOUND_DROWN, SOUND_SCREAM1,
	SOUND_SCREAM2, SOUND_GLASS, SOUND_DRILL, SOUND_ROCKET, SOUND_ITEM, SOUND_POWERUP, SOUND_ALERT, SOUND_SQUISH, SOUND_SQUISH2, SOUND_SPLAT1,
	SOUND_SPLASH1, SOUND_SPLASH2, SOUND_WATER, SOUND_THUNDER, SOUND_THUNDER2, SOUND_BOING, SOUND_SWING, SOUND_HISS, SOUND_DOH, SOUND_HURT,
	SOUND_DEATH, SOUND_AGONY, SOUND_SCARED, SOUND_GASP, SOUND_SCREAM3, SOUND_SQUEAL, SOUND_RICOCHET, SOUND_ROCK_FALL, SOUND_SPRAY, SOUND_CLICK,
	SOUND_SHELLC, SOUND_SH_DROP, SOUND_WATER_DROP, SOUND_SLIDING, SOUND_OBJ_FALL, SOUND_WOOD_CRACK, SOUND_FOOTSTEP, SOUND_SNOW_STEP, SOUND_ICE_CRACK, SOUND_RELOAD,
	SOUND_FALLING, SOUND_HORN, SOUND_DOOR_OPEN, SOUND_DOOR_CLOSE, SOUND_KICK_BALL, SOUND_FLUSH, SOUND_GULP, SOUND_ZOMBIE1, SOUND_ZOMBIE2, SOUND_ZOMBIE3,
	SOUND_ZOMBIE4, SOUND_ZOMBIE5, SOUND_SQUEAK, SOUND_BEEP, SOUND_SINK, SOUND_METAL_DOOR, SOUND_DOORBELL, SOUND_HELICOPTER, SOUND_RAT_SQUEAK, SOUND_HURT2,
	SOUND_FLY_BUZZ, SOUND_EATING, SOUND_BUBBLE, SOUND_NEON_SIGN, SOUND_SM_SPLAT, SOUND_POLICE, SOUND_ALARM, SOUND_SCRATCH, SOUND_HANDGUN, NUM_SOUNDS};

// looping sounds
enum {SOUND_LOOP_FIRE, SOUND_LOOP_RAIN, SOUND_LOOP_WIND, SOUND_LOOP_UNDERWATER, NUM_LOOP_SOUNDS};


struct openal_orient {

	vector3d at, up;
	openal_orient() {}
	openal_orient(vector3d const &at_, vector3d const &up_) : at(at_), up(up_) {}
};


struct sound_params_t {

	point pos;
	float gain, pitch;
	int sound_id;
	bool rel_to_listener;

	sound_params_t() : gain(1.0), pitch(1.0), sound_id(-1), rel_to_listener(0) {}
	sound_params_t(point const &P, unsigned sid, float g, float p, bool r) : pos(P), gain(g), pitch(p), sound_id(sid), rel_to_listener(r) {}
	float get_loudness() const;
	bool read_from_file(FILE *fp);
	void write_to_cobj_file(std::ostream &out) const;
};


struct delayed_sound_t : public sound_params_t {

	int id, time;

	delayed_sound_t() : id(-1), time(0) {}
	delayed_sound_t(sound_params_t const &p, int i, int t) : sound_params_t(p), id(i), time(t) {}
};


class openal_buffer {

	unsigned buffer;
	float time;

public:
	openal_buffer(unsigned buffer_=0) : buffer(buffer_), time(0.0) {}
	//~openal_buffer() {free_buffer();}
	bool load_check();
	bool load_from_file(std::string const &fn);
	bool load_from_file_std_path(std::string const &fn);
	bool load_from_memory(void const *data, size_t length);
	bool load_from_waveform(int waveshape, float frequency, float phase, float duration);
	bool is_valid() const {return (buffer > 0);}
	unsigned get_buffer_ix() const {return buffer;}
	void alloc();
	void free_buffer();
};


class buffer_manager_t {

	vector<openal_buffer> buffers;

public:
	openal_buffer &get_buffer(unsigned id) {assert(id < buffers.size()); return buffers[id];}

	size_t add_buffer(bool alloc) {
		size_t const ix(buffers.size());
		buffers.push_back(openal_buffer());
		if (alloc) {buffers.back().alloc();}
		return ix;
	}
	unsigned add_file_buffer(std::string const &fn);

	void clear() {
		for (unsigned i = 0; i < buffers.size(); ++i) {buffers[i].free_buffer();}
		buffers.clear();
	}
};


class openal_source {

	unsigned source;
	sound_params_t params;

public:
	openal_source(unsigned source_=0) : source(source_) {}
	//~openal_source() {free_source();}
	bool is_valid  () const {return (source > 0);}
	bool is_playing() const;
	bool is_active () const;
	bool is_playing_sound(unsigned sid) const {return ((int)sid == params.sound_id && is_playing());}
	float get_loudness() const {return (is_active() ? params.get_loudness() : 0.0);}
	
	void alloc();
	void free_source();
	void setup(openal_buffer const &buffer, point const &pos, unsigned sound_id, float gain=1.0, float pitch=1.0,
		bool looping=0, bool rel_to_listener=0, vector3d const &vel=zero_vector);
	void set_pos(point const &pos);
	void set_gain(float gain);
	void set_buffer(openal_buffer const &buffer) {set_buffer_ix(buffer.get_buffer_ix());}
	void set_buffer_ix(unsigned buffer_ix);
	void play_if_not_playing() const;
	void play()   const;
	void stop()   const;
	void pause()  const;
	void rewind() const;
	bool check_for_active_sound(point const &pos, float radius, float min_gain=0.0) const;
};


class source_manager_t {

	vector<openal_source> sources;
	unsigned next_source;

public:
	std::set<unsigned> used_this_frame;

	source_manager_t() : next_source(0) {}
	void create_channels(unsigned num_channels);
	unsigned new_source();
	openal_source &get_least_loud_source();
	openal_source &get_oldest_source();
	openal_source &get_inactive_source();
	openal_source &get_source(unsigned id) {assert(id < sources.size()); return sources[id];}
	openal_source const &get_source(unsigned id) const {assert(id < sources.size()); return sources[id];}
	void clear();
	bool is_playing   (unsigned id) const {return get_source(id).is_playing();}
	void play_source  (unsigned id) const {get_source(id).play  ();}
	void stop_source  (unsigned id) const {get_source(id).stop  ();}
	void pause_source (unsigned id) const {get_source(id).pause ();}
	void rewind_source(unsigned id) const {get_source(id).rewind();}
	bool is_playing_sound(unsigned sid) const;
	bool check_for_active_sound(point const &pos, float radius, float min_gain=0.0) const;
};


struct placed_sound_t {

	int sound_id;
	sound_params_t params;
	sensor_t sensor;
	//multi_trigger_t triggers;

	placed_sound_t() : sound_id(-1) {}
	placed_sound_t(unsigned id, sound_params_t const &params_, sensor_t const &sensor_=sensor_t()) : sound_id(id), params(params_), sensor(sensor_) {}
	bool enabled() const {return (sound_id >= 0);}
	void next_frame();
	void write_to_cobj_file(std::ostream &out, std::string const &name) const;
};


unsigned get_sound_id_for_file(std::string const &fn);
std::string const &get_sound_name(unsigned id);
void set_sound_loop_state(unsigned id, bool play, float volume=0.0);
bool check_for_active_sound(point const &pos, float radius, float min_gain=0.0);
void add_placed_sound(std::string const &fn, sound_params_t const &params, sensor_t const &sensor=sensor_t());
void write_placed_sounds_to_cobj_file(std::ostream &out);
void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient);
void set_openal_listener_as_player();
void gen_sound(unsigned id, point const &pos, float gain=1.0, float pitch=1.0, bool rel_to_listener=0, vector3d const &vel=zero_vector, bool skip_if_already_playing=0);
void gen_sound_random_var(unsigned id, point const &pos, float gain=1.0, float pitch=1.0);
void gen_delayed_sound(float delay, unsigned id, point const &pos, float gain=1.0, float pitch=1.0, bool rel_to_listener=0); // no vel
void gen_delayed_from_player_sound(unsigned id, point const &pos, float gain=1.0, float pitch=1.0);
void proc_delayed_and_placed_sounds();
void play_thunder(point const &pos, float gain, float delay);
void play_switch_weapon_sound();
void play_switch_wmode_sound();
void init_openal(int &argc, char** argv);
void exit_openal();
int read_sound_file(std::string const &name);

