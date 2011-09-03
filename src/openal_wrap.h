// 3D World - OpenAL Interface Code Header
// by Frank Gennari
// 9/1/11
#include "3DWorld.h"


enum {SOUND_EXPLODE=0, SOUND_GUNSHOT, SOUND_SHOTGUN, SOUND_FIREBALL, SOUND_BURNING, SOUND_DROWN, SOUND_SCREAM1, SOUND_SCREAM2, SOUND_GLASS,
	SOUND_DRILL, SOUND_ROCKET, SOUND_ITEM, SOUND_POWERUP, NUM_SOUNDS};


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
	bool load_from_file_std_path(std::string const &fn);
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
	openal_buffer &get_buffer(unsigned id) {assert(id < buffers.size()); return buffers[id];}

	unsigned add_buffer(bool alloc) {
		unsigned const ix(buffers.size());
		buffers.push_back(openal_buffer());
		if (alloc) buffers.back().alloc();
		return ix;
	}
	unsigned add_file_buffer(std::string const &fn);

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


class source_manager_t {

	vector<openal_source> sources;
	unsigned next_source;

public:
	source_manager_t() : next_source(0) {}

	void create_channels(unsigned num_channels) {
		clear();
		for (unsigned i = 0; i < num_channels; ++i) new_source();
	}
	unsigned new_source() {
		unsigned const ix(sources.size());
		sources.push_back(openal_source());
		sources.back().alloc();
		return ix;
	}
	openal_source &get_source() { // round robin
		assert(!sources.empty());
		if (next_source >= sources.size()) next_source = 0; // wraparound
		return sources[next_source++];
	}
	void clear() {
		for (unsigned i = 0; i < sources.size(); ++i) sources[i].free();
		sources.clear();
		next_source = 0;
	}
};


void setup_openal_listener(point const &pos, vector3d const &vel, openal_orient const &orient);
void set_openal_listener_as_player();
void gen_sound(unsigned id, point const &pos, vector3d const &vel=zero_vector, float pitch=1.0, float gain=1.0, bool looping=0);
void init_openal(int &argc, char** argv);


