// 3D World - FILE Utility Functions
// by Frank Gennari
// 3/6/12

#pragma once

#include "3DWorld.h"
#include <zlib.h>

using std::string;

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
		if (fp ) {fclose(fp); fp = nullptr;}
		if (gzf) {gzclose(gzf); gzf = nullptr;}
	}
	static string get_extension(string const &filename) {return filename.substr(filename.find_last_of(".") + 1);}
	static bool   is_gz_file   (string const &filename) {return (get_extension(filename) == "gz");}
};

struct binary_file_reader : public binary_file_io {
	bool open(string const &filename) {return binary_file_io::open(filename, "rb", "reading");}

	bool read(void *ptr, size_t sz, size_t count) {
		if      (fp ) {return (fread (ptr, sz, count, fp) == count);}
		else if (gzf) {return (gzread(gzf, ptr, sz*count) == sz*count);}
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

