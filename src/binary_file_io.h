// 3D World - FILE Utility Functions
// by Frank Gennari
// 3/6/12

#pragma once

#include "3DWorld.h"
#include <zlib.h>
#include <fstream>

using std::string;
using std::istream;
using std::ostream;

class binary_file_io {
protected:
	FILE *fp;
	gzFile gzf;
public:
	binary_file_io() : fp(nullptr), gzf(nullptr) {}
	~binary_file_io() {close();}

	bool open(string const &filename, char const *const mode, string const &purpose) {
		if (filename.empty()) return 0;
		if (is_gz_file(filename)) {gzf = gzopen(filename.c_str(), mode);} else {fp = fopen(filename.c_str(), mode);}
		if (is_valid()) return 1;
		std::cerr << "Failed to open file " << filename << " for " << purpose << ".";
		return 0;
	}
	bool is_valid() const {return (fp || gzf);}
	void close() {
		if (fp ) {
			checked_fclose(fp);
			fp = nullptr;
		}
		if (gzf) {
			if (gzclose(gzf) != 0) {
				std::cerr << "Error: gzclose() call failed" << std::endl;
				exit(0); // how fatal should this be?
			}
			gzf = nullptr;
		}
	}
	static string get_extension(string const &filename) {return filename.substr(filename.find_last_of(".") + 1);}
	static bool   is_gz_file   (string const &filename) {return (get_extension(filename) == "gz");}
};

struct binary_file_reader : public binary_file_io {
	bool open(string const &filename) {return binary_file_io::open(filename, "rb", "reading");}

	bool read(void *ptr, size_t sz, size_t count) {
		if      (fp ) {return (fread (ptr, sz, count, fp) == count);}
		else if (gzf) {return (gzread(gzf, ptr, sz*count) == int(sz*count));}
		else {assert(0);} // no file opened
		return 0;
	}
};
struct binary_file_writer : public binary_file_io {
	bool open(string const &filename) {return binary_file_io::open(filename, "wb", "writing");}

	bool write(void const *ptr, size_t sz, size_t count) {
		if      (fp ) {return (fwrite (ptr, sz, count, fp) == count);}
		else if (gzf) {return (gzwrite(gzf, ptr, sz*count) == int(sz*count));}
		else {assert(0);} // no file opened
		return 0;
	}
};

// ************ read/write helpers ************

template<typename T> void write_val(ostream &out, T val) {
	out.write((const char *)&val, sizeof(T));
}
template<typename T> void read_val(istream &in, T &val) {
	in.read((char *)&val, sizeof(T));
}

inline void write_uint(ostream &out, unsigned val) {
	write_val(out, val);
}
inline unsigned read_uint(istream &in ) {
	unsigned val(0);
	read_val(in, val);
	return val;
}

template<typename V> void write_vector(ostream &out, V const &v) {
	write_uint(out, (unsigned)v.size());
	out.write((const char *)&v.front(), (std::streamsize)v.size()*sizeof(typename V::value_type));
}
template<typename V> void read_vector(istream &in, V &v) {
	v.clear();
	v.resize(read_uint(in));
	in.read((char *)&v.front(), (std::streamsize)v.size()*sizeof(typename V::value_type));
}

inline void write_string(ostream &out, string const &s) {
	write_uint(out, (unsigned)s.size());
	out.write(s.c_str(), s.size());
}
inline void read_string(istream &in, string &s) {
	vector<char> str;
	read_vector(in, str);
	s = string(str.begin(), str.end());
}

