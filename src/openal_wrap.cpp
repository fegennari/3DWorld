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


struct openal_orient {

	vector3d at, up;
	openal_orient() {}
	openal_orient(vector3d const &at_, vector3d const &up_) : at(at_), up(up_) {}
};


class openal_buffer {

	unsigned buffer;

public:
	openal_buffer() : buffer(0) {}
	~openal_buffer() {free();}
	bool is_valid() const {return (buffer > 0);}
	unsigned get_buffer_ix() const {return buffer;}
	
	void alloc() {
		 alGenBuffers(1, &buffer);
    
		if (alGetError() != AL_NO_ERROR) {
			cerr << "Error creating OpenAL buffers" << endl;
			exit(1);
		}
		assert(is_valid());
	}
	void free() {
		if (is_valid()) {
			alDeleteBuffers(1, &buffer);
			buffer = 0;
		}
	}
};


class openal_source {

	unsigned source;

public:
	openal_source() : source(0) {}
	~openal_source() {free();}
	bool is_valid() const {return (source > 0);}
	
	void alloc() {
		alGenSources(1, &source);

		if (alGetError() != AL_NO_ERROR) {
		   cerr << "Error creating OpenAL sources" << endl;
		   exit(1);
		}
		assert(is_valid());
	}
	void free() {
		if (is_valid()) {
			alDeleteSources(1, &source);
			source = 0;
		}
	}
	void setup(openal_buffer const &buffer, point const &pos, vector3d const &vel, float pitch=1.0, float gain=1.0, bool looping=0) {
		assert(is_valid() && buffer.is_valid());
		alSourcef (source, AL_PITCH,    pitch);
		alSourcef (source, AL_GAIN,     gain);
		alSourcefv(source, AL_POSITION, &pos.x);
		alSourcefv(source, AL_VELOCITY, &vel.x);
		alSourcei (source, AL_LOOPING,  looping);
		set_buffer_ix(buffer.get_buffer_ix());
	}
	void set_buffer_ix(unsigned buffer_ix) {alSourcei(source, AL_BUFFER, buffer_ix);}
	void play()   const {alSourcePlay  (source);}
	void stop()   const {alSourceStop  (source);}
	void pause()  const {alSourcePause (source);}
	void rewind() const {alSourceRewind(source);}
};


void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient) {

	alListenerfv(AL_POSITION,    &pos.x);
    alListenerfv(AL_VELOCITY,    &vel.x);
    alListenerfv(AL_ORIENTATION, &orient.at.x);
}


void load_wave_file(string const &fn, openal_buffer &buffer) {

#if 0 // old (deprecated) version
	ALsizei size, freq;
	ALenum format;
	ALboolean loop;
	ALvoid *data;
	alutLoadWAVFile((char *)fn.c_str(), &format, &data, &size, &freq, &loop);
    alBufferData(buffer.get_buffer_ix(), format, data, size, freq);
    alutUnloadWAV(format, data, size, freq);
#else

#endif
}


void run_openal_test() {

	point    const listen_pos(0.0, 0.0, 4.0);
	vector3d const listen_vel(0.0, 0.0, 0.0);
	openal_orient const listen_orient(vector3d(0.0, 0.0, 1.0), vector3d(0.0, 1.0, 0.0));
	point    const source_pos(-2.0, 0.0, 0.0);
	vector3d const source_vel(0.0, 0.0, 0.0);
	setup_openal_listener(listen_pos, listen_vel, listen_orient);
	openal_buffer buffer;
	openal_source source;
	buffer.alloc();
	source.alloc();
	load_wave_file("test.wav", buffer);
	source.setup(buffer, source_pos, source_vel, 1.0, 1.0, 0);
	// write
}


void openal_hello_world() {

	unsigned buffer(alutCreateBufferHelloWorld());
	openal_source source;
	source.alloc();
	source.set_buffer_ix(buffer);
	source.play();
	alutSleep(1.0);
	//free_openal_buffer(buffer);
}


void init_openal(int &argc, char** argv) {

	return; // not yet ready
	alutInit(&argc, argv);
	alGetError(); // ignore any previous errors
	openal_hello_world();
}


