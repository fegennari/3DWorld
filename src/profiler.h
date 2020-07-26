// 3D World - Shared shadow map related classes
// by Frank Gennari
// 5/23/14
#pragma once

#include <string>
#include <chrono>

using namespace std::chrono;

class highres_timer_t {
	std::string name;
	bool enabled;
	high_resolution_clock::time_point timer1;
	high_resolution_clock clock;
public:
	highres_timer_t(char const *const name_,  bool enabled_=1) : name(name_), timer1(clock.now()), enabled(enabled_) {}
	highres_timer_t(std::string const &name_, bool enabled_=1) : name(name_), timer1(clock.now()), enabled(enabled_) {}
	~highres_timer_t() {end();}
	void end();
};

