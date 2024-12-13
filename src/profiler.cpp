// 3D World - Performance Timing Profiler
// by Frank Gennari
// 4/20/13

#include "3DWorld.h"
#include "profiler.h"

using std::string;

void maybe_update_loading_screen(const char *str);
int omp_get_thread_num_3dw();


template <typename T> class timing_profiler {

	struct entry_t {
		unsigned count;
		T time, tmax;
		entry_t() : count(0), time(0), tmax(0) {}
		void add(T t) {++count; time += t; tmax = max(tmax, t);}
	};
	map<string, entry_t> entries;

public:
	bool enabled;

	timing_profiler() : enabled(0) {}
	void clear() {entries.clear();}

	void register_time(const char *str, T delta_time, bool no_loading_screen) {
#pragma omp critical(timer_update)
		{
			if (enabled) {entries[str].add(delta_time);}
			else {cout << str << " time = " << delta_time << endl;}
		}
		if (!no_loading_screen && delta_time > 0 && omp_get_thread_num_3dw() == 0) {maybe_update_loading_screen(str);} // only call on main thread when time has elapsed
	}
	void stats() const {
		if (entries.empty()) return;
		cout << "name count total max average" << endl;
		unsigned max_name(0);
		for (auto i = entries.begin(); i != entries.end(); ++i) {max_name = max(max_name, (unsigned)i->first.size());}

		for (auto i = entries.begin(); i != entries.end(); ++i) {
			string const spaces((max_name - i->first.size()), ' ');
			cout << i->first << spaces << ": " << i->second.count << "\t" << i->second.time << "\t"
					<< i->second.tmax << "\t" << float(i->second.time)/float(i->second.count) << endl;
		}
	}
};

timing_profiler<int> global_profiler;
timing_profiler<float> global_highres_profiler;

void toggle_timing_profiler() {global_profiler.enabled ^= 1; global_highres_profiler.enabled ^= 1;}
void register_timing_value(const char *str, int delta_time, bool no_loading_screen) {global_profiler.register_time(str, delta_time, no_loading_screen);}

void timing_profiler_stats() {
	global_profiler.stats();
	global_profiler.clear();
	global_highres_profiler.stats();
	global_highres_profiler.clear();
}

void highres_timer_t::end() {
	if (!enabled || name.empty()) return;
	float const elapsed(get_delta_secs(clock.now(), timer1));
	global_highres_profiler.register_time(name.c_str(), 1000.0*elapsed, no_loading_screen); // print in ms
	name.clear(); // make sure we don't double count this
}

