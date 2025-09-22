// 3D World - Simple File Logger
// by Frank Gennari
// 9/21/25
#pragma once

#include <string>
#include <iostream>
#include <fstream>

class logger_t {
	std::ofstream log;

	void open_log_file() { // open file on first write
		if (log.good()) return; // alredy open
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

