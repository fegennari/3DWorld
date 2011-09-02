// 3D World - OpenAL Interface Code
// by Frank Gennari
// 8/28/11
#include "openal_wrap.h"
#include <iostream>
#include <assert.h>
#include <al.h>
#include <alc.h>
#include <AL/alut.h>

using namespace std;


buffer_manager_t buffer_manager;


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
	buffer = alutCreateBufferFromFile(fn.c_str());
	return load_check();
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

void openal_source::setup(openal_buffer const &buffer, point const &pos, vector3d const &vel, float pitch, float gain, bool looping) {
	assert(is_valid() && buffer.is_valid());
	alSourcef (source, AL_PITCH,    pitch);
	alSourcef (source, AL_GAIN,     gain);
	alSourcefv(source, AL_POSITION, &pos.x);
	alSourcefv(source, AL_VELOCITY, &vel.x);
	alSourcei (source, AL_LOOPING,  looping);
	set_buffer_ix(buffer.get_buffer_ix());
}

void openal_source::set_buffer_ix(unsigned buffer_ix) {alSourcei(source, AL_BUFFER, buffer_ix);}

void openal_source::play()   const {alSourcePlay  (source);}
void openal_source::stop()   const {alSourceStop  (source);}
void openal_source::pause()  const {alSourcePause (source);}
void openal_source::rewind() const {alSourceRewind(source);}


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


// test and init
void run_openal_test() {

	point    const listen_pos(0.0, 0.0, 4.0);
	vector3d const listen_vel(0.0, 0.0, 0.0);
	openal_orient const listen_orient(vector3d(0.0, 0.0, 1.0), vector3d(0.0, 1.0, 0.0));
	point    const source_pos(-2.0, 0.0, 0.0);
	vector3d const source_vel(0.0, 0.0, 0.0);
	setup_openal_listener(listen_pos, listen_vel, listen_orient);
	openal_buffer buffer;
	openal_source source;
	source.alloc();
	
	if (buffer.load_from_file("test.wav")) {
		source.setup(buffer, source_pos, source_vel, 1.0, 1.0, 0);
		source.play();
		alut_sleep(2.0);
	}
	source.free();
	buffer.free();
}


void openal_hello_world() {

	openal_buffer buffer(alutCreateBufferHelloWorld());
	openal_source source;
	source.alloc();
	source.set_buffer(buffer);
	source.play();
	alut_sleep(1.0);
	source.free();
	buffer.free();
}


void init_openal(int &argc, char** argv) {

	return; // not yet ready

	if (!alutInit(&argc, argv)) {
		check_and_print_alut_error();
		cerr << "alutInit failed" << endl;
		exit(1);
	}
	alGetError(); // ignore any previous errors
	//cout << "Supported OpenAL types: " << alutGetMIMETypes(ALUT_LOADER_BUFFER) << endl;
	openal_hello_world();
	run_openal_test();
}


void exit_openal() {

	if (!alutExit()) check_and_print_alut_error();
}


