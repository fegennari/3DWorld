// 3D World - Logging
// by Frank Gennari
// 1/24/26

#include "3DWorld.h"
#include <fstream>

class logger_t {
	std::ofstream log;

	void open_log_file() { // open file on first write
		if (log.good()) return; // already open
		log.open("3DWorld.log");
	}
public:
	void log_str(std::string const &str, bool add_newline=1) {
		open_log_file();
		log << str;
		if (add_newline) {log << std::endl;}
	}
	template <typename T> logger_t& operator<<(const T& v) {
		open_log_file();
		log << v;
		return *this;
	}
	logger_t& operator<<(std::ostream& (*manip)(std::ostream&)) { // overload for manipulators like std::endl
		open_log_file();
		manip(log);
		return *this;
	}
	~logger_t() {log.close();}
};
logger_t global_logger;

void global_logger_log(string const &str) {global_logger.log_str(str);}

void log_location(point const &pos) {

	static std::ofstream out;
	static bool inited(0);

	if (!inited) {
		out.open("positions.log.txt");
		inited = 1;
	}
	assert(out.good());
	out << pos.x << " " << pos.y << " " << pos.z << endl;
}

