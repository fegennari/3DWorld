// 3D World - Performance Timing Profiler
// by Frank Gennari
// 4/20/13

#include "3DWorld.h"

using std::string;


class timing_profiler {

	struct entry_t {
		unsigned count;
		int time, tmax;
		entry_t() : count(0), time(0), tmax(0) {}
		void add(int t) {++count; time += t; tmax = max(tmax, t);}
	};

	map<string, entry_t> entries;

public:
	bool enabled;

	timing_profiler() : enabled(0) {}
	void clear() {entries.clear();}

	void register_time(const char *str, int delta_time) {
		if (enabled) {
			entries[str].add(delta_time);
		}
		else {
			cout << str << " time = " << delta_time << endl;
		}
	}
	void stats() const {
		cout << "name count total max average" << endl;
		unsigned max_name(0);
		for (map<string, entry_t>::const_iterator i = entries.begin(); i != entries.end(); ++i) {max_name = max(max_name, i->first.size());}

		for (map<string, entry_t>::const_iterator i = entries.begin(); i != entries.end(); ++i) {
			string const spaces((max_name - i->first.size()), ' ');
			cout << i->first << spaces << ": " << i->second.count << "\t" << i->second.time << "\t"
					<< i->second.tmax << "\t" << float(i->second.time)/float(i->second.count) << endl;
		}
	}
};

timing_profiler global_profiler;


void toggle_timing_profiler() {
	global_profiler.enabled ^= 1;
}

void register_timing_value(const char *str, int delta_time) {
	global_profiler.register_time(str, delta_time);
}

void timing_profiler_stats() {
	global_profiler.stats();
	global_profiler.clear();
}



