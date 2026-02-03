// OpenGL Pattern Matcher Graphical Display - Image Export
// Must be compiled as C, not C++, to enable JPEG support
// by Frank Gennari 6/10/02

#include "globals.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "gl_includes.h"

using namespace std;

int write_jpeg_data(string const &fn, unsigned char const *const data, unsigned width, unsigned height, bool invert_y);
bool write_rgb_bmp_image(string const &fn, unsigned char *data, unsigned width, unsigned height, unsigned ncolors);

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

int screenshot(unsigned window_width, unsigned window_height, char const *const file_path, bool write_bmp) {

	static unsigned ss_id(0);
	ostringstream oss;
	if (file_path != nullptr) {oss << file_path;}
	oss << "screenshot" << ss_id++ << "." << (write_bmp ? "bmp" : "jpg");
	string const fn(oss.str());
	cout << "Writing screenshot image " << fn << endl;
	print_text_onscreen_default("Writing screenshot " + fn);
	vector<unsigned char> buf((window_width+1)*(window_height+1)*3); // allocate extra size for alignment purposes - is this really needed?
	read_pixels(window_width, window_height, buf);
	if (write_bmp) {return write_rgb_bmp_image(fn, &buf.front(), window_width, window_height, 3);} // BMP, RGB (ncolors=3) (40ms)
	else           {return write_jpeg_data    (fn, &buf.front(), window_width, window_height, 1);} // JPEG, invert_y=1 (110ms libjpeg, 226ms stb)
}

