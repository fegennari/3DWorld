// 3D World - File Reader Base Class
// by Frank Gennari
// 11/1/14

#ifndef _FILE_READER_H_
#define _FILE_READER_H_

#include <fstream>
#include <iostream>
#include "3DWorld.h"

class base_file_reader {

protected:
	std::string filename;
	FILE *fp; // Note: we use a FILE* here instead of an ifstream because it's ~2.2x faster in MSVS
	static unsigned const MAX_CHARS = 1024;
	char buffer[MAX_CHARS];
	bool verbose;

	bool open_file();
	void close_file();
	int get_next_char() {assert(fp); return get_char(fp);}
	void unget_last_char(int c);
	bool fast_isspace(char c)  const {return (c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '/r');}
	bool fast_isdigit(char c)  const {return (c >= '0' && c <= '9');}
	int get_char(FILE *fp_)    const;
	int get_char(std::ifstream &in) const {return in.get();}

	template<typename T> void read_to_newline(T &stream, string *str=NULL) const {
		bool prev_was_escape(0);

		while (1) {
			int const c(get_char(stream));
			if ((!prev_was_escape && c == '\n') || c == '\0' || c == EOF) return;
			prev_was_escape = (c == '\\'); // handle escape character at end of line
			if (str && !prev_was_escape) {str->push_back(char(c));}
		}
		assert(0); // never gets here
	}
	template<typename T> void read_str_to_newline(T &stream, string &str) const {
		str.resize(0);

		while (1) {
			int const c(get_char(stream));
			if (c == '\n' || c == '\0' || c == EOF) break; // end of file or line
			if (!fast_isspace(c) || !str.empty()) str.push_back(c);
		}
		return;
	}
	int fast_atoi(char *str) const;
	bool read_int(int &v);
	bool read_string(char *s, unsigned max_len);

public:
	base_file_reader(std::string const &fn) : filename(fn), fp(NULL), verbose(0) {assert(!fn.empty());}
	~base_file_reader() {close_file();}
};

#endif // _FILE_READER_H_

