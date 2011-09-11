// 3D World - OpenAL Interface Code
// by Frank Gennari
// 8/28/11
// Sounds from http://www.findsounds.com
#include "openal_wrap.h"
#include <iostream>
#include <assert.h>
#include <al.h>
#include <alc.h>
#include <AL/alut.h>

using namespace std;


unsigned const NUM_CHANNELS = 8;
string const sounds_path("sounds/");

buffer_manager_t sounds;
source_manager_t sources, looping_sources;
vector<delayed_sound_t> delayed_sounds;

extern int frame_counter, iticks;


// supported: au, wav
// not supported: mp3, aif
void setup_sounds() {

	cout << endl << "Loading Sounds"; cout.flush();
	sounds.add_file_buffer("burning.wav"    ); // SOUND_BURNING
	sounds.add_file_buffer("rain1.wav"      ); // SOUND_RAIN1
	sounds.add_file_buffer("wind1.wav"      ); // SOUND_WIND1
	sounds.add_file_buffer("explosion1.au"  ); // SOUND_EXPLODE
	sounds.add_file_buffer("gunshot.wav"    ); // SOUND_GUNSHOT
	sounds.add_file_buffer("shotgun.wav"    ); // SOUND_SHOTGUN
	sounds.add_file_buffer("fireball.wav"   ); // SOUND_FIREBALL
	sounds.add_file_buffer("drown.wav"      ); // SOUND_DROWN
	sounds.add_file_buffer("scream1.wav"    ); // SOUND_SCREAM1
	sounds.add_file_buffer("scream2.wav"    ); // SOUND_SCREAM2
	sounds.add_file_buffer("glass_break.wav"); // SOUND_GLASS
	sounds.add_file_buffer("drill.wav"      ); // SOUND_DRILL
	sounds.add_file_buffer("rocket.au"      ); // SOUND_ROCKET
	sounds.add_file_buffer("item.wav"       ); // SOUND_ITEM
	sounds.add_file_buffer("powerup.wav"    ); // SOUND_POWERUP
	sounds.add_file_buffer("alert.wav"      ); // SOUND_ALERT
	sounds.add_file_buffer("squish.wav"     ); // SOUND_SQUISH
	sounds.add_file_buffer("squish2.wav"    ); // SOUND_SQUISH2
	sounds.add_file_buffer("splat1.wav"     ); // SOUND_SPLAT1
	sounds.add_file_buffer("splash1.wav"    ); // SOUND_SPLASH1
	sounds.add_file_buffer("splash2.wav"    ); // SOUND_SPLASH2
	sounds.add_file_buffer("water.wav"      ); // SOUND_WATER
	sounds.add_file_buffer("thunder.wav"    ); // SOUND_THUNDER
	sounds.add_file_buffer("boing.wav"      ); // SOUND_BOING
	sounds.add_file_buffer("swing.wav"      ); // SOUND_SWING
	sounds.add_file_buffer("hiss.wav"       ); // SOUND_HISS
	sounds.add_file_buffer("doh.wav"        ); // SOUND_DOH
	cout << endl;

	// create sources
	sources.create_channels(NUM_CHANNELS);
	looping_sources.create_channels(NUM_LOOP_SOUNDS);

	for (unsigned i = 0; i < NUM_LOOP_SOUNDS; ++i) { // looping_source i is bound to sounds buffer i
		openal_source &source(looping_sources.get_source(i));
		source.setup(sounds.get_buffer(i), all_zeros, 1.0, 1.0, 1, 1);
	}
}

void set_sound_loop_state(unsigned id, bool play) {
	bool const playing(looping_sources.is_playing(id));

	if (play && !playing) { // start
		looping_sources.play_source(id);
	}
	else if (!play && playing) { // stop
		looping_sources.stop_source(id);
	}
}

void alut_sleep(float seconds) {
	alutSleep(seconds);
}


bool had_al_error  () {return (alGetError  () != AL_NO_ERROR);}
bool had_alut_error() {return (alutGetError() != AL_NO_ERROR);}

bool check_and_print_alut_error() { // returns 1 on error

	ALenum const error_id(alutGetError());

	if (error_id != AL_NO_ERROR) {
		cerr << "alut error: " << alutGetErrorString(error_id) << endl;
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

void openal_buffer::free() {
	if (is_valid()) alDeleteBuffers(1, &buffer);
	buffer = 0;
	time   = 0.0;
}

bool openal_buffer::load_check() {
	if (check_and_print_alut_error()) {
		free();
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
	buffer = alutCreateBufferFromFileImage(data, length);
	return load_check();
}

bool openal_buffer::load_from_waveform(int waveshape, float frequency, float phase, float duration) {
	assert(frequency > 0.0 && duration > 0.0);
	buffer = alutCreateBufferWaveform(waveshape, frequency, phase, duration);
	time   = duration;
	return load_check();
}


unsigned buffer_manager_t::add_file_buffer(std::string const &fn) {

	unsigned const ix(buffers.size());
	buffers.push_back(openal_buffer());
	
	if (!buffers.back().load_from_file_std_path(fn)) {
		cerr << "Failed to load sound file: " << fn << endl;
		exit(1);
	}
	cout << "."; cout.flush();
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

void openal_source::free() {
	if (is_valid()) {
		alDeleteSources(1, &source);
		source = 0;
	}
}

void openal_source::setup(openal_buffer const &buffer, point const &pos, float gain, float pitch,
	bool looping, bool rel_to_listener, vector3d const &vel)
{
	assert(is_valid() && buffer.is_valid());
	alSourcef (source, AL_PITCH,    pitch);
	alSourcef (source, AL_GAIN,     gain);
	alSourcefv(source, AL_POSITION, &pos.x);
	alSourcefv(source, AL_VELOCITY, &vel.x);
	alSourcei (source, AL_LOOPING,  looping);
	alSourcei (source, AL_SOURCE_RELATIVE, rel_to_listener);
	set_buffer_ix(buffer.get_buffer_ix());
	params = sound_params_t(pos, gain, pitch, rel_to_listener);
}

void openal_source::set_buffer_ix(unsigned buffer_ix) {alSourcei(source, AL_BUFFER, buffer_ix);}

void openal_source::play_if_not_playing() const {if (!is_playing()) alSourcePlay(source);}
void openal_source::play()   const {alSourcePlay(source);}
void openal_source::stop()   const {alSourceStop  (source);}
void openal_source::pause()  const {alSourcePause (source);}
void openal_source::rewind() const {alSourceRewind(source);}

void openal_source::blocking_play() const {
	play();
	//if (buffer.time > 0.0) {alut_sleep(buffer.time); return;}

	do {
		alut_sleep(0.01); // sleep 10ms
	} while (is_active());
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


// source_manager_t

void source_manager_t::create_channels(unsigned num_channels) {
	clear();
	
	for (unsigned i = 0; i < num_channels; ++i) {
		new_source();
	}
}

unsigned source_manager_t::new_source() {
	unsigned const ix(sources.size());
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
			if (next_source == i) next_source++; // move past this source
			return sources[i];
		}
	}
	return get_oldest_source(); // all were active
}

void source_manager_t::clear() {
	for (unsigned i = 0; i < sources.size(); ++i) {
		sources[i].free();
	}
	sources.clear();
	next_source = 0;
}


// listner code
void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient) {

	alListenerfv(AL_POSITION,    &pos.x);
    alListenerfv(AL_VELOCITY,    &vel.x);
    alListenerfv(AL_ORIENTATION, &orient.at.x);
}


void set_openal_listener_as_player() {

	// Note: velocity is zero for now because the actual player velocity is not constant,
	// and passing it into the listener setup as a constant is incorrect
	vector3d const cvel(zero_vector);
	setup_openal_listener(get_camera_all(), cvel, openal_orient(get_vdir_all(), get_upv_all()));
}


// non-blocking
void gen_sound(unsigned id, point const &pos, float gain, float pitch, bool rel_to_listener, vector3d const &vel) {

	//RESET_TIME;
	point const listener(get_camera_pos());
	float const dist(distance_to_camera(pos));
	bool const close(dist < CAMERA_RADIUS);
	if (!close && id != SOUND_DROWN && id != SOUND_SPLASH1 && id != SOUND_SPLASH2 && id != SOUND_WATER && (is_underwater(pos) || is_underwater(listener))) return;
#if 0
	openal_source &source(sources.get_inactive_source());
	if (!close && source.is_playing()) return; // already playing - don't stop it
#else
	openal_source &source(sources.get_least_loud_source());
	float const loudness(gain/max(SMALL_NUMBER, dist));
	//cout << "f: " << frame_counter << " id: " << id << " loud: " << loudness << " s loud: " << source.get_loudness() << " keep: " << !(loudness < max(0.01f, source.get_loudness())) << endl; // testing
	if (loudness < max(0.01f, source.get_loudness())) return; // too soft
#endif
	static int last_frame(0);

	if (frame_counter != last_frame) { // start new frame
		sources.used_this_frame.clear();
		last_frame = frame_counter;
	}
	if (sources.used_this_frame.find(id) != sources.used_this_frame.end()) return; // duplicate sound this frame
	sources.used_this_frame.insert(id);

	if (!close) {
		int cindex;
		bool const line_of_sight(!check_coll_line(pos, listener, cindex, -1, 1, 0));
		if (!line_of_sight) gain *= 0.25; // attenuate by 4x if there is no line of sight between source and listener
	}
	if (source.is_active()) source.stop(); // stop if already playing
	set_openal_listener_as_player();
	source.setup(sounds.get_buffer(id), pos, gain, pitch, 0, rel_to_listener, vel); // not looping
	source.play();
	//PRINT_TIME("Play Sound");
}


void gen_delayed_sound(float delay, unsigned id, point const &pos, float gain, float pitch, bool rel_to_listener) {

	if (delay == 0.0) {
		gen_sound(id, pos, gain, pitch, rel_to_listener);
	}
	else {
		assert(delay > 0.0);
		int const delay_time(int(delay*TICKS_PER_SECOND + 0.5)); // round to the nearest tick
		sound_params_t params(pos, gain, pitch, rel_to_listener);
		delayed_sounds.push_back(delayed_sound_t(params, id, delay_time));
	}
}


void proc_delayed_sounds() {

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


void openal_hello_world() {

	openal_buffer buffer(alutCreateBufferHelloWorld());
	openal_source source;
	source.alloc();
	source.set_buffer(buffer);
	source.blocking_play();
	source.free();
	buffer.free();
}


void init_openal(int &argc, char** argv) {

	//return; // not yet ready

	if (!alutInit(&argc, argv)) {
		check_and_print_alut_error();
		cerr << "alutInit failed" << endl;
		exit(1);
	}
	alGetError(); // ignore any previous errors
	setup_sounds();
	//cout << "Supported OpenAL types: " << alutGetMIMETypes(ALUT_LOADER_BUFFER) << endl;
	//openal_hello_world();
}


void exit_openal() {

	if (!alutExit()) check_and_print_alut_error();
}


