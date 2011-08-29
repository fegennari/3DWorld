// 3D World - OpenAL Interface Code
// by Frank Gennari
// 8/28/11
//#include "3DWorld.h"
#include "3DWorld.h" // must be first
#include <iostream>
#include <assert.h>
#include <al.h>
#include <alc.h>
#include <AL/alut.h>

using namespace std;


void load_wave_file(string const &fn, unsigned buffer) {

#if 0 // old (deprecated) version
	ALsizei size, freq;
	ALenum format;
	ALboolean loop;
	ALvoid *data;
	alutLoadWAVFile((char *)fn.c_str(), &format, &data, &size, &freq, &loop);
    alBufferData(buffer, format, data, size, freq);
    alutUnloadWAV(format, data, size, freq);
#else

#endif
}


unsigned gen_openal_buffer() {

	unsigned buffer(0);
    alGenBuffers(1, &buffer);
    
    if (alGetError() != AL_NO_ERROR) {
        cerr << "Error creating OpenAL buffers" << endl;
        exit(1);
    }
	assert(buffer > 0);
	return buffer;
}


unsigned gen_openal_source() {

	unsigned source(0);
    alGenSources(1, &source);

    if (alGetError() != AL_NO_ERROR) {
       cerr << "Error creating OpenAL sources" << endl;
       exit(1);
    }
	assert(source > 0);
	return source;
}


void free_openal_buffer(unsigned &buffer) {

	assert(buffer > 0); // too strict?
	alDeleteBuffers(1, &buffer);
	buffer = 0;
}


void free_openal_source(unsigned &source) {

	assert(source > 0); // too strict?
	alDeleteSources(1, &source);
	source = 0;
}


void setup_openal_source(unsigned source, unsigned buffer, point const &pos, vector3d const &vel, float pitch=1.0, float gain=1.0, bool looping=0) {

	assert(source > 0 && buffer > 0);
	alSourcef (source, AL_PITCH,    pitch);
    alSourcef (source, AL_GAIN,     gain);
    alSourcefv(source, AL_POSITION, &pos.x);
    alSourcefv(source, AL_VELOCITY, &vel.x);
    alSourcei (source, AL_BUFFER,   buffer);
    alSourcei (source, AL_LOOPING,  looping);
}


struct openal_orient {
	vector3d v1, v2;
	openal_orient() {}
	openal_orient(vector3d const &v1_, vector3d const &v2_) : v1(v1_), v2(v2_) {}
};


void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient) {

	alListenerfv(AL_POSITION,    &pos.x);
    alListenerfv(AL_VELOCITY,    &vel.x);
    alListenerfv(AL_ORIENTATION, &orient.v1.x);
}


void run_openal_test() {

	point    const listen_pos(0.0, 0.0, 4.0);
	vector3d const listen_vel(0.0, 0.0, 0.0);
	openal_orient const listen_orient(vector3d(0.0, 0.0, 1.0), vector3d(0.0, 1.0, 0.0));
	point    const source_pos(-2.0, 0.0, 0.0);
	vector3d const source_vel(0.0, 0.0, 0.0);
	setup_openal_listener(listen_pos, listen_vel, listen_orient);
	
	unsigned buffer(gen_openal_buffer());
	unsigned source(gen_openal_source());
	load_wave_file("test.wav", buffer);
	setup_openal_source(source, buffer, source_pos, source_vel, 1.0, 1.0, 0);
	// write
	free_openal_buffer(buffer);
	free_openal_source(source);
}


void openal_hello_world() {

	unsigned buffer(alutCreateBufferHelloWorld());
	unsigned source(gen_openal_source());
	alSourcei(source, AL_BUFFER, buffer);
	alSourcePlay(source);
	//free_openal_buffer(buffer);
	//free_openal_source(source);
}


void init_openal(int &argc, char** argv) {

	return; // not yet ready
	alutInit(&argc, argv);
	//openal_hello_world();
}


