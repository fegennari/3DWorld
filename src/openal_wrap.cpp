// 3D World - OpenAL Interface Code
// by Frank Gennari
// 8/28/11
// Sounds from http://www.findsounds.com
#include "openal_wrap.h"
#include "function_registry.h"
#include <iostream>
#include <assert.h>
#include <thread>
#include <chrono>

using namespace std;

extern float CAMERA_RADIUS;

#define ENABLE_OPENAL // comment this out to disable OpenAL sound support

#ifdef ENABLE_OPENAL
#if defined(_WIN32) && !defined(__MINGW32__) // Note: the Windows OpenAL 1.1 SDK doesn't have the AL directory, but the openal-soft GitHub repo does
#include <al.h>
#include <alc.h>
#else
#include <AL/al.h>
#include <AL/alc.h>
#endif
#include <AL/alut.h>

unsigned const NUM_CHANNELS = 8;
string const sounds_path("sounds/");

extern bool disable_sound;
extern int frame_counter, iticks;

float const loop_sound_gains  [NUM_LOOP_SOUNDS] = {0.5, 0.1, 0.1, 0.1};
float const loop_sound_pitches[NUM_LOOP_SOUNDS] = {1.0, 1.0, 1.0, 1.0};

void sleep_for_ms(unsigned milliseconds) {
	std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
}


void placed_sound_t::write_to_cobj_file(ostream &out, string const &name) const {

	sensor.write_to_cobj_file(out);
	out << "place_sound " << name << " ";
	params.write_to_cobj_file(out); // includes endl
	sensor.write_end_sensor_to_cobj_file(out);
}


class sound_manager_t {

	buffer_manager_t sounds;
	source_manager_t sources, looping_sources;
	vector<delayed_sound_t> delayed_sounds;
	vector<placed_sound_t> placed_sounds;
	map<string, unsigned> name_to_id_map;
	vector<string> sound_names;
	int last_frame;

public:
	sound_manager_t() : last_frame(0) {}
	openal_source &get_least_loud_source() {return sources.get_least_loud_source();}
	openal_buffer &get_buffer(unsigned id) {return sounds.get_buffer(id);}
	void add_delayed_sound(sound_params_t const &params, unsigned id, int delay_time) {delayed_sounds.push_back(delayed_sound_t(params, id, delay_time));}
	string const &get_name(unsigned id) const {assert(id < sound_names.size()); return sound_names[id];}
	bool is_playing_sound(unsigned sid) const {return (sources.is_playing_sound(sid) || looping_sources.is_playing_sound(sid));}

	unsigned add_new_sound(string const &fn) {
		unsigned const ix(sounds.add_file_buffer(fn));
		bool const did_ins(name_to_id_map.insert(make_pair(fn, ix)).second);
		assert(did_ins); // all sound filenames must be unique
		assert(ix == sound_names.size());
		sound_names.push_back(fn);
		return ix;
	}
	unsigned find_or_add_sound(string const &fn) {
		auto it(name_to_id_map.find(fn));
		if (it != name_to_id_map.end()) return it->second; // found
		return add_new_sound(fn);
	}
	void add_placed_sound(string const &fn, sound_params_t const &params, sensor_t const &sensor) {
		placed_sounds.push_back(placed_sound_t(find_or_add_sound(fn), params, sensor));
	}
	void write_placed_sounds_to_cobj_file(ostream &out) const {
		for (auto i = placed_sounds.begin(); i != placed_sounds.end(); ++i) {i->write_to_cobj_file(out, get_name(i->sound_id));}
	}

	// supported: au, wav
	// not supported: mp3, aif
	void setup_sounds() {

		cout << endl << "Loading Sounds"; cout.flush();
		// loop sounds
		add_new_sound("burning.wav"    ); // SOUND_BURNING
		add_new_sound("rain1.wav"      ); // SOUND_RAIN1
		add_new_sound("wind1.wav"      ); // SOUND_WIND1
		add_new_sound("underwater.wav" ); // SOUND_UNDERWATER
		// regular sounds
		add_new_sound("explosion1.au"  ); // SOUND_EXPLODE
		add_new_sound("gunshot.wav"    ); // SOUND_GUNSHOT
		add_new_sound("shotgun.wav"    ); // SOUND_SHOTGUN
		add_new_sound("fireball.wav"   ); // SOUND_FIREBALL
		add_new_sound("drown.wav"      ); // SOUND_DROWN
		add_new_sound("scream1.wav"    ); // SOUND_SCREAM1
		add_new_sound("scream2.wav"    ); // SOUND_SCREAM2
		add_new_sound("glass_break.wav"); // SOUND_GLASS
		add_new_sound("drill.wav"      ); // SOUND_DRILL
		add_new_sound("rocket.au"      ); // SOUND_ROCKET
		add_new_sound("item.wav"       ); // SOUND_ITEM
		add_new_sound("powerup.wav"    ); // SOUND_POWERUP
		add_new_sound("alert.wav"      ); // SOUND_ALERT
		add_new_sound("squish.wav"     ); // SOUND_SQUISH
		add_new_sound("squish2.wav"    ); // SOUND_SQUISH2
		add_new_sound("splat1.wav"     ); // SOUND_SPLAT1
		add_new_sound("splash1.wav"    ); // SOUND_SPLASH1
		add_new_sound("splash2.wav"    ); // SOUND_SPLASH2
		add_new_sound("water.wav"      ); // SOUND_WATER
		add_new_sound("thunder.wav"    ); // SOUND_THUNDER
		add_new_sound("thunder2.wav"   ); // SOUND_THUNDER2
		add_new_sound("boing.wav"      ); // SOUND_BOING
		add_new_sound("swing.wav"      ); // SOUND_SWING
		add_new_sound("hiss.wav"       ); // SOUND_HISS
		add_new_sound("doh.wav"        ); // SOUND_DOH
		add_new_sound("hurt.wav"       ); // SOUND_HURT
		add_new_sound("death.wav"      ); // SOUND_DEATH
		add_new_sound("agony.au"       ); // SOUND_AGONY
		add_new_sound("scared.wav"     ); // SOUND_SCARED
		add_new_sound("gasp.wav"       ); // SOUND_GASP
		add_new_sound("scream3.wav"    ); // SOUND_SCREAM3
		add_new_sound("squeal.wav"     ); // SOUND_SQUEAL
		add_new_sound("ricochet.wav"   ); // SOUND_RICOCHET
		add_new_sound("rockfall.wav"   ); // SOUND_ROCK_FALL // bricksfall.wav
		add_new_sound("spray.wav"      ); // SOUND_SPRAY
		add_new_sound("click.wav"      ); // SOUND_CLICK
		add_new_sound("shell_casing.wav");// SOUND_SHELLC
		add_new_sound("short_drop.wav" ); // SOUND_SH_DROP (unused)
		add_new_sound("water_drop.wav" ); // SOUND_WATER_DROP (unused, cartoony)
		add_new_sound("sliding.wav"    ); // SOUND_SLIDING
		add_new_sound("obj_fall.wav"   ); // SOUND_OBJ_FALL
		add_new_sound("wood_crack.wav" ); // SOUND_WOOD_CRACK
		add_new_sound("footstep.wav"   ); // SOUND_FOOTSTEP
		add_new_sound("snow_step.wav"  ); // SOUND_SNOW_STEP
		add_new_sound("ice_crack.wav"  ); // SOUND_ICE_CRACK
		add_new_sound("reload.wav"     ); // SOUND_RELOAD
		add_new_sound("falling.wav"    ); // SOUND_FALLING
		add_new_sound("horn.wav"       ); // SOUND_HORN
		add_new_sound("door_open.wav"  ); // SOUND_DOOR_OPEN
		add_new_sound("door_close.wav" ); // SOUND_DOOR_CLOSE
		add_new_sound("kickball.wav"   ); // SOUND_KICK_BALL
		add_new_sound("toilet_flush.wav");// SOUND_FLUSH
		add_new_sound("gulp.wav"       ); // SOUND_GULP
		add_new_sound("zombie1.wav"    ); // SOUND_ZOMBIE1
		add_new_sound("zombie2.wav"    ); // SOUND_ZOMBIE2
		add_new_sound("zombie3.wav"    ); // SOUND_ZOMBIE3
		add_new_sound("zombie4.wav"    ); // SOUND_ZOMBIE4
		add_new_sound("zombie5.wav"    ); // SOUND_ZOMBIE5
		add_new_sound("squeak.wav"     ); // SOUND_SQUEAK
		add_new_sound("beep.wav"       ); // SOUND_BEEP
		add_new_sound("flowingwater.wav");// SOUND_SINK
		add_new_sound("metal_door.wav" ); // SOUND_METAL_DOOR
		add_new_sound("doorbell.wav"   ); // SOUND_DOORBELL
		add_new_sound("helicopter.wav" ); // SOUND_HELICOPTER
		add_new_sound("rat_squeak.wav" ); // SOUND_RAT_SQUEAK
		add_new_sound("hurt2.wav"      ); // SOUND_HURT2
		add_new_sound("fly_buzz.wav"   ); // SOUND_FLY_BUZZ
		add_new_sound("eating.wav"     ); // SOUND_EATING
		add_new_sound("bubble.au"      ); // SOUND_BUBBLE
		add_new_sound("neon_sign_sm.wav"); // SOUND_NEON_SIGN
		add_new_sound("small_splat.wav"); // SOUND_SM_SPLAT
		add_new_sound("police.wav"     ); // SOUND_POLICE
		add_new_sound("alarm.wav"      ); // SOUND_ALARM
		cout << endl;

		// create sources
		sources.create_channels(NUM_CHANNELS);
		looping_sources.create_channels(NUM_LOOP_SOUNDS);

		for (unsigned i = 0; i < NUM_LOOP_SOUNDS; ++i) { // looping_source i is bound to sounds buffer i
			openal_source &source(looping_sources.get_source(i));
			source.setup(sounds.get_buffer(i), all_zeros, i, loop_sound_gains[i], loop_sound_pitches[i], 1, 1);
		}
	}

	void set_loop_state(unsigned id, bool play, float volume) { // volume=0.0 => use previous value

		if (disable_sound) return;
		assert(id < NUM_LOOP_SOUNDS);
		bool const playing(looping_sources.is_playing(id));
		if (play && volume > 0.0) {looping_sources.get_source(id).set_gain(CLIP_TO_01(volume)*loop_sound_gains[id]);}
		if      ( play && !playing) {looping_sources.play_source(id);} // start
		else if (!play &&  playing) {looping_sources.stop_source(id);} // stop
	}

	bool check_for_duplicate(int id) {

		if (frame_counter != last_frame) { // start new frame
			sources.used_this_frame.clear();
			last_frame = frame_counter;
		}
		if (sources.used_this_frame.find(id) != sources.used_this_frame.end()) return 1; // duplicate sound this frame
		sources.used_this_frame.insert(id);
		return 0;
	}

	void proc_delayed() {

		for (unsigned i = 0; i < delayed_sounds.size(); ++i) {
			delayed_sound_t &ds(delayed_sounds[i]);
			ds.time -= iticks;
			if (ds.time > 0) continue; // continue to delay
			gen_sound(ds.id, ds.pos, ds.gain, ds.pitch, ds.rel_to_listener);
			swap(delayed_sounds[i], delayed_sounds.back());
			delayed_sounds.pop_back();
			--i; // wraparound ok
		}
	}
	void proc_placed() {
		for (auto i = placed_sounds.begin(); i != placed_sounds.end(); ++i) {i->next_frame();}
	}

	bool check_for_active_sound(point const &pos, float radius, float min_gain=0.0) const {
		return sources.check_for_active_sound(pos, radius, min_gain);
	}
};

sound_manager_t sound_manager;

unsigned get_sound_id_for_file(string const &fn) {return sound_manager.find_or_add_sound(fn);}
string const &get_sound_name(unsigned id) {return sound_manager.get_name(id);}
void set_sound_loop_state(unsigned id, bool play, float volume) {sound_manager.set_loop_state(id, play, volume);}
bool check_for_active_sound(point const &pos, float radius, float min_gain) {return sound_manager.check_for_active_sound(pos, radius, min_gain);}
void add_placed_sound(string const &fn, sound_params_t const &params, sensor_t const &sensor) {sound_manager.add_placed_sound(fn, params, sensor);}
void write_placed_sounds_to_cobj_file(ostream &out) {sound_manager.write_placed_sounds_to_cobj_file(out);}

bool had_al_error  () {return (alGetError  () != AL_NO_ERROR);}
bool had_alut_error() {return (alutGetError() != AL_NO_ERROR);}

bool check_and_print_alut_error(const char *msg) { // returns 1 on error

	ALenum const error_id(alutGetError());

	if (error_id != AL_NO_ERROR) {
		cerr << "alut error in " << msg << ": " << alutGetErrorString(error_id) << endl;
		return 1;
	}
	return 0;
}
bool check_and_print_al_error(const char *msg) { // returns 1 on error

	ALenum error_id(alGetError());

	if (error_id != AL_NO_ERROR) {
		cerr << "OpenAL error in " << msg << ": ID " << error_id << endl; // Note: alutGetErrorString doesn't work on all AL error codes
		return 1;
	}
	return 0;
}


// openal_buffer
void openal_buffer::alloc() {
	assert(!is_valid());
	alGenBuffers(1, &buffer);
    
	if (had_al_error()) {
		cerr << "Error creating OpenAL buffers" << endl;
		exit(1);
	}
	assert(is_valid());
}

void openal_buffer::free_buffer() {
	if (is_valid()) {alDeleteBuffers(1, &buffer);}
	buffer = 0;
	time   = 0.0;
}

bool openal_buffer::load_check() {
	if (check_and_print_alut_error("load_check")) {
		free_buffer();
		return 0;
	}
	return 1;
}

bool openal_buffer::load_from_file(string const &fn) {
	assert(!fn.empty());
	buffer = alutCreateBufferFromFile(fn.c_str());
	return load_check();
}

bool openal_buffer::load_from_file_std_path(std::string const &fn) {
	assert(!fn.empty());
	return load_from_file(sounds_path + fn);
}

bool openal_buffer::load_from_memory(void const *data, size_t length) {
	assert(data);
	assert(length > 0);
	buffer = alutCreateBufferFromFileImage(data, (ALsizei)length);
	return load_check();
}

bool openal_buffer::load_from_waveform(int waveshape, float frequency, float phase, float duration) {
	assert(frequency > 0.0 && duration > 0.0);
	buffer = alutCreateBufferWaveform(waveshape, frequency, phase, duration);
	time   = duration;
	return load_check();
}


unsigned buffer_manager_t::add_file_buffer(std::string const &fn) {

	unsigned const ix((unsigned)buffers.size());
	buffers.push_back(openal_buffer());
	openal_buffer &buf(buffers.back());
	// check here because an existing OpenAL error will cause loading to fail
	check_and_print_al_error("pre_load");

	// check sounds directory first, and current directory second
	if (!buf.load_from_file_std_path(fn) && !buf.load_from_file(fn)) {
		cerr << "Failed to load sound file: " << fn << endl;
		exit(1);
	}
	if (frame_counter <= 1) {cout << "."; cout.flush();} // only show progress in pre-load
	return ix;
}


// openal_source
void openal_source::alloc() {
	assert(!is_valid());
	alGenSources(1, &source);

	if (had_al_error()) {
		cerr << "Error creating OpenAL sources" << endl;
		exit(1);
	}
	assert(is_valid());
}

void openal_source::free_source() {
	if (is_valid()) {
		alDeleteSources(1, &source);
		source = 0;
	}
}

void openal_source::setup(openal_buffer const &buffer, point const &pos, unsigned sound_id, float gain, float pitch,
	bool looping, bool rel_to_listener, vector3d const &vel)
{
	assert(is_valid() && buffer.is_valid());
	set_pos(pos);
	set_gain(gain);
	alSourcef (source, AL_PITCH,    pitch);
	alSourcefv(source, AL_VELOCITY, &vel.x);
	alSourcei (source, AL_LOOPING,  looping);
	alSourcei (source, AL_SOURCE_RELATIVE, rel_to_listener);
	set_buffer_ix(buffer.get_buffer_ix());
	params = sound_params_t(pos, sound_id, gain, pitch, rel_to_listener);
}
void openal_source::set_pos(point const &pos) {
	alSourcefv(source, AL_POSITION, &pos.x);
}
void openal_source::set_gain(float gain) {
	alSourcef(source, AL_GAIN, gain);
}

void openal_source::set_buffer_ix(unsigned buffer_ix) {alSourcei(source, AL_BUFFER, buffer_ix);}

void openal_source::play_if_not_playing() const {
	if (!is_playing()) {play();}
}
void openal_source::play()   const {alSourcePlay  (source);}
void openal_source::stop()   const {alSourceStop  (source);}
void openal_source::pause()  const {alSourcePause (source);}
void openal_source::rewind() const {alSourceRewind(source);}

void openal_source::blocking_play() const {
	play();
	do {sleep_for_ms(10);} while (is_active()); // sleep 10ms
}

int get_source_state(unsigned source) {
	int state;
	alGetSourcei(source, AL_SOURCE_STATE, &state);
	return state;
}

bool openal_source::is_active() const {
	if (!is_valid()) return 0;
	int const state(get_source_state(source));
	return (state == AL_PLAYING || state == AL_PAUSED);
}

bool openal_source::is_playing() const {
	if (!is_valid()) return 0;
	int const state(get_source_state(source));
	return (state == AL_PLAYING);
}

bool openal_source::check_for_active_sound(point const &pos, float radius, float min_gain) const {
	if (!is_valid()) return 0;
	if (params.gain < min_gain) return 0;
	if (!dist_less_than(params.pos, pos, radius)) return 0;
	return is_playing(); // do this test last to avoid the library call
}


// source_manager_t

void source_manager_t::create_channels(unsigned num_channels) {
	clear();
	for (unsigned i = 0; i < num_channels; ++i) {new_source();}
}

unsigned source_manager_t::new_source() {
	unsigned const ix((unsigned)sources.size());
	sources.push_back(openal_source());
	sources.back().alloc();
	return ix;
}

openal_source &source_manager_t::get_least_loud_source() {
	assert(!sources.empty());
	float least_loud(0.0);
	unsigned least_loud_source(0);

	for (unsigned i = 0; i < sources.size(); ++i) {
		float const loudness(sources[i].get_loudness());
		if (loudness == 0.0) return sources[i]; // not active

		if (least_loud == 0.0 || loudness < least_loud) {
			least_loud = loudness;
			least_loud_source = i;
		}
	}
	return sources[least_loud_source];
}

openal_source &source_manager_t::get_oldest_source() { // round robin
	assert(!sources.empty());
	if (next_source >= sources.size()) next_source = 0; // wraparound
	return sources[next_source++];
}

openal_source &source_manager_t::get_inactive_source() {
	for (unsigned i = 0; i < sources.size(); ++i) {
		if (!sources[i].is_active()) {
			if (next_source == i) {next_source++;} // move past this source
			return sources[i];
		}
	}
	return get_oldest_source(); // all were active
}

void source_manager_t::clear() {
	for (auto i = sources.begin(); i != sources.end(); ++i) {i->free_source();}
	sources.clear();
	next_source = 0;
}

bool source_manager_t::is_playing_sound(unsigned sid) const {
	for (auto i = sources.begin(); i != sources.end(); ++i) { // loop over every source and check if it's currently playing this sound
		if (i->is_playing_sound(sid)) return 1;
	}
	return 0;
}

bool source_manager_t::check_for_active_sound(point const &pos, float radius, float min_gain) const {
	for (auto i = sources.begin(); i != sources.end(); ++i) {
		if (i->check_for_active_sound(pos, radius, min_gain)) return 1;
	}
	return 0;
}


void placed_sound_t::next_frame() {
	if (disable_sound || !enabled()) return;
	if (check_for_active_sound(params.pos, 0.1*CAMERA_RADIUS)) return; // assume this sound is already playing (not a perfect check, but it's easy)
	if (sensor.enabled() && !sensor.check_active()) return;
	gen_sound(sound_id, params.pos, params.gain, params.pitch, params.rel_to_listener); // no velocity
}


// listner code
void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient) {
	if (disable_sound) return;
	alListenerfv(AL_POSITION,    &pos.x);
    alListenerfv(AL_VELOCITY,    &vel.x);
    alListenerfv(AL_ORIENTATION, &orient.at.x);
}


// non-blocking
void gen_sound(unsigned id, point const &pos, float gain, float pitch, bool rel_to_listener, vector3d const &vel, bool skip_if_already_playing) {

	//RESET_TIME;
	if (disable_sound) return;
	point const listener(get_camera_pos());
	float const dist(distance_to_camera(pos));
	bool const close(dist < CAMERA_RADIUS), is_underwater_sound(id == SOUND_DROWN || id == SOUND_SPLASH1 || id == SOUND_SPLASH2 || id == SOUND_WATER);
	if (!close && !is_underwater_sound && world_mode == WMODE_GROUND && (is_underwater(pos) || is_underwater(listener))) return; // can't hear sounds that are under water/when under water
	if (skip_if_already_playing && sound_manager.is_playing_sound(id)) return; // this sound is already playing; is it possible to update the position and gain?
#if 0
	openal_source &source(sources.get_inactive_source());
	if (!close && source.is_playing()) return; // already playing - don't stop it
#else
	openal_source &source(sound_manager.get_least_loud_source());
	float const loudness(gain/max(SMALL_NUMBER, dist));
	if (loudness < max(0.01f, source.get_loudness())) return; // too soft
#endif
	if (sound_manager.check_for_duplicate(id)) return; // duplicate sound this frame

	if (!close && world_mode == WMODE_GROUND) {
		int cindex;
		bool const line_of_sight(!check_coll_line(pos, listener, cindex, -1, 1, 0));
		if (!line_of_sight) {gain *= 0.25;} // attenuate by 4x if there is no line of sight between source and listener
	}
	if (source.is_active()) {source.stop();} // stop if already playing
	set_openal_listener_as_player();
	source.setup(sound_manager.get_buffer(id), pos, id, gain, pitch, 0, rel_to_listener, vel); // not looping
	source.play();
	//PRINT_TIME("Play Sound");
}
void gen_sound_random_var(unsigned id, point const &pos, float gain, float pitch) { // with minor random variationin gain and pitch
	static rand_gen_t rgen;
	gen_sound(id, pos, gain*rgen.rand_uniform(0.75, 1.3), pitch*rgen.rand_uniform(0.9, 1.11));
}

void gen_delayed_sound(float delay, unsigned id, point const &pos, float gain, float pitch, bool rel_to_listener) { // delay in seconds
	if (disable_sound) return;

	if (delay < 0.01) {gen_sound(id, pos, gain, pitch, rel_to_listener);} // less than 10ms - play now
	else {
		assert(delay > 0.0);
		sound_manager.add_delayed_sound(sound_params_t(pos, id, gain, pitch, rel_to_listener), id, round_fp(delay*TICKS_PER_SECOND)); // round to the nearest tick
	}
}

void proc_delayed_and_placed_sounds() {
	sound_manager.proc_delayed();
	sound_manager.proc_placed();
}


void init_openal(int &argc, char** argv) {
	if (disable_sound) return;

	if (!alutInit(&argc, argv)) {
		check_and_print_alut_error("init");
		cerr << "alutInit failed" << endl;
		exit(1);
	}
	alGetError(); // ignore any previous errors
	sound_manager.setup_sounds();
	//cout << "Supported OpenAL types: " << alutGetMIMETypes(ALUT_LOADER_BUFFER) << endl;
	//openal_hello_world();
}

void exit_openal() {
	if (!alutExit()) {check_and_print_alut_error("exit");}
}

#else // !ENABLE_OPENAL

unsigned get_sound_id_for_file(std::string const &fn) {return 0;}
string const empty_str;
string const &get_sound_name(unsigned id) {return empty_str;} // for writing to cobj file
void set_sound_loop_state(unsigned id, bool play, float volume) {}
bool check_for_active_sound(point const &pos, float radius, float min_gain) {return 0;}
void add_placed_sound(string const &fn, sound_params_t const &params, sensor_t const &sensor) {}
void write_placed_sounds_to_cobj_file(ostream &out) {}
void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient) {}
void gen_sound(unsigned id, point const &pos, float gain, float pitch, bool rel_to_listener, vector3d const &vel, bool skip_if_already_playing) {}
void gen_delayed_sound(float delay, unsigned id, point const &pos, float gain, float pitch, bool rel_to_listener) {}
void proc_delayed_and_placed_sounds() {}
void init_openal(int &argc, char** argv) {}
void exit_openal() {}

#endif // ENABLE_OPENAL

float sound_params_t::get_loudness() const {
	return gain/max(SMALL_NUMBER, distance_to_camera(pos));
}
bool sound_params_t::read_from_file(FILE *fp) {
	gain = 1.0; pitch = 1.0;
	return (fscanf(fp, "%f%f%f%f%f", &pos.x, &pos.y, &pos.z, &gain, &pitch) >= 3); // pos.x pos.y pos.z [gain [pitch]]
}
void sound_params_t::write_to_cobj_file(std::ostream &out) const {
	out << pos.raw_str() << " " << gain << " " << pitch << endl;
}

void play_switch_weapon_sound() {gen_sound(SOUND_CLICK,  get_camera_pos(), 1.0, 0.7);}
void play_switch_wmode_sound () {gen_sound(SOUND_RELOAD, get_camera_pos(), 1.0, 1.0);}

void gen_delayed_from_player_sound(unsigned id, point const &pos, float gain, float pitch) {
	float const SPEED_OF_SOUND = 16.0; // in world units per second
	float const dist(distance_to_camera(pos));
	gen_delayed_sound(((dist < CAMERA_RADIUS) ? 0.0 : dist/SPEED_OF_SOUND), id, pos, gain, pitch);
}

void play_thunder(point const &pos, float gain, float delay) {
	float const pitch(rand_uniform(0.8, 1.2)); // variable pitch
	gen_delayed_sound(delay, ((rand()&1) ? (unsigned)SOUND_THUNDER : (unsigned)SOUND_THUNDER2), pos, gain, pitch);
}

void set_openal_listener_as_player() {
	// Note: velocity is zero for now because the actual player velocity is not constant,
	// and passing it into the listener setup as a constant is incorrect
	vector3d const cvel(zero_vector);
	setup_openal_listener(get_camera_all(), cvel, openal_orient(get_vdir_all(), get_upv_all()));
}

int read_sound_file(string const &name) {
	if (name == "none" || name == "null") return -1; // no sound
	return get_sound_id_for_file(name);
}
