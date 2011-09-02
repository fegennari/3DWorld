// 3D World - OpenAL Interface Code Header
// by Frank Gennari
// 9/1/11
#include "3DWorld.h"


struct openal_orient {

	vector3d at, up;
	openal_orient() {}
	openal_orient(vector3d const &at_, vector3d const &up_) : at(at_), up(up_) {}
};


class openal_buffer {

	unsigned buffer;
	float time;

public:
	openal_buffer(unsigned buffer_=0) : buffer(buffer_), time(0.0) {}
	//~openal_buffer() {free();}
	bool load_check();
	bool load_from_file(std::string const &fn);
	bool load_from_memory(void const *data, size_t length);
	bool load_from_waveform(int waveshape, float frequency, float phase, float duration);
	bool is_valid() const {return (buffer > 0);}
	unsigned get_buffer_ix() const {return buffer;}
	void alloc();
	void free();
};


class buffer_manager_t {

	vector<openal_buffer> buffers;

public:
	openal_buffer &new_buffer(bool alloc) {
		buffers.push_back(openal_buffer());
		if (alloc) buffers.back().alloc();
		return buffers.back();
	}
	void clear() {
		for (unsigned i = 0; i < buffers.size(); ++i) buffers[i].free();
		buffers.clear();
	}
};


class openal_source {

	unsigned source;

public:
	openal_source(unsigned source_=0) : source(source_) {}
	//~openal_source() {free();}
	bool is_valid() const {return (source > 0);}
	
	void alloc();
	void free();
	void setup(openal_buffer const &buffer, point const &pos, vector3d const &vel, float pitch=1.0, float gain=1.0, bool looping=0);
	void set_buffer(openal_buffer const &buffer) {set_buffer_ix(buffer.get_buffer_ix());}
	void set_buffer_ix(unsigned buffer_ix);
	void play()   const;
	void stop()   const;
	void pause()  const;
	void rewind() const;
};


void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient);
void set_openal_listener_as_player();
void init_openal(int &argc, char** argv);


