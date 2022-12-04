// OpenGL Pattern Matcher Graphical Display - Image Export
// Must be compiled as C, not C++, to enable JPEG support
// by Frank Gennari 6/10/02

#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "gl_includes.h"

using namespace std;

int write_jpeg_data(unsigned width, unsigned height, FILE *fp, unsigned char const *const data, bool invert_y);
bool write_rgb_bmp_image(FILE *fp, string const &fn, unsigned char *data, unsigned width, unsigned height, unsigned ncolors);
void checked_fclose(FILE *fp);

void print_text_onscreen_default(string const &text);


void read_depth_buffer(unsigned window_width, unsigned window_height, vector<float> &depth, bool normalize) {

	depth.resize(window_width*window_height, 0.0);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, window_width, window_height, GL_DEPTH_COMPONENT, GL_FLOAT, &depth.front());
	if (!normalize) return; // done
	float dmin(1.0), dmax(0.0);
	for (auto i = depth.begin(); i != depth.end(); ++i) {dmin = min(*i, dmin); dmax = max(*i, dmax);}
	if (dmax == dmin) return; // all the same depth, can't normalize
	float const dmult(1.0f/(dmax - dmin));
	for (auto i = depth.begin(); i != depth.end(); ++i) {*i = (*i - dmin)*dmult;}
}

void read_pixels(unsigned window_width, unsigned window_height, vector<unsigned char> &buf) {

	unsigned const num_bytes(window_width*window_height*3);
	if (buf.size() < num_bytes) {buf.resize(num_bytes);}
#if 0 // for debugging depth buffer
	vector<float> depth;
	read_depth_buffer(window_width, window_height, depth, 1); // normalized
	for (unsigned i = 0; i < depth.size(); ++i) {buf[3*i+0] = buf[3*i+1] = buf[3*i+2] = (unsigned char)(255.0*depth[i]);}
#else
	glReadBuffer(GL_FRONT);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, &buf.front());
#endif
}


string get_screenshot_filename(char const *const file_path, string const &extension, unsigned &id) {
	ostringstream oss;
	if (file_path != nullptr) {oss << file_path;}
	oss << "screenshot" << id++ << "." << extension;
	string const fn(oss.str());
	cout << "Writing screenshot image " << fn << endl;
	print_text_onscreen_default("Writing screenshot " + fn);
	return fn;
}

int screenshot(unsigned window_width, unsigned window_height, char const *const file_path, bool write_bmp) {

	static unsigned ss_id(0);
	string const fn(get_screenshot_filename(file_path, (write_bmp ? "bmp" : "jpg"), ss_id));
	FILE* fp(fopen(fn.c_str(), "wb"));

	if (fp == nullptr) {
		printf("Error writing screenshot '%s'.\n", fn.c_str());
		return 0;
	}	
	vector<unsigned char> buf((window_width+1)*(window_height+1)*3); // allocate extra size for alignment purposes - is this really needed?
	read_pixels(window_width, window_height, buf);
	int ret(0);
	if (write_bmp) {ret = write_rgb_bmp_image(fp, "<screenshot>", &buf.front(), window_width, window_height, 3);} // RGB BMP (40ms)
	else           {ret = write_jpeg_data(window_width, window_height, fp, &buf.front(), 1);} // JPEG (110ms)
	checked_fclose(fp);
	return ret;
}

