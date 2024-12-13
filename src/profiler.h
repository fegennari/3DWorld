// 3D World - Performance Timing Profiler
// by Frank Gennari
// 5/23/14
#pragma once

#include <string>
#include <chrono>

using namespace std::chrono;

class highres_timer_t { // should this share a base class with timer_t?
	std::string name;
	bool enabled, no_loading_screen;
	high_resolution_clock::time_point timer1;
	high_resolution_clock clock;
public:
	highres_timer_t(char const *const name_,  bool enabled_=1, bool nls=0) : name(name_), enabled(enabled_), no_loading_screen(nls), timer1(clock.now()) {}
	highres_timer_t(std::string const &name_, bool enabled_=1, bool nls=0) : name(name_), enabled(enabled_), no_loading_screen(nls), timer1(clock.now()) {}
	~highres_timer_t() {end();}
	void end();
};

inline float get_delta_secs(high_resolution_clock::time_point b, high_resolution_clock::time_point a) {return duration_cast<duration<float>>(b - a).count();}

