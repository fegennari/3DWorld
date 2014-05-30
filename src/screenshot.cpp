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
#include "gl_includes.h"

using namespace std;

int write_jpeg_data(unsigned width, unsigned height, FILE *fp, unsigned char const *const data, bool invert_y);


FILE *open_screenshot_file(char *file_path, char *basename) {

	FILE *fp;
	if (basename == NULL) return NULL;

	if (file_path != NULL) {
		string const filename(string(file_path) + basename);
		fp = fopen(filename.c_str(), "wb");
	}
	else {
		fp = fopen(basename, "wb");
	}
	if (fp == NULL) {
		printf("Error writing screenshot '%s'.\n", basename);
		return NULL;
	}
	return fp;
}


void read_pixels(unsigned window_width, unsigned window_height, vector<unsigned char> &buf) {

	glReadBuffer(GL_FRONT);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, &buf.front());
}


int screenshot(unsigned window_width, unsigned window_height, char *file_path) {

	FILE *fp = open_screenshot_file(file_path, "screenshot.raw");
	if (fp == NULL) return 0;
	vector<unsigned char> buf(window_width*window_height*3);
	read_pixels(window_width, window_height, buf);

	for (unsigned i = 0; i < window_height; ++i) {
		unsigned const offset((window_height-i-1)*window_width);
		int const num_write(fwrite(&buf[3*offset], 3, window_width, fp));
		assert(num_write == window_width);
	}
	fclose(fp);
	return 1;
}


int write_jpeg(unsigned window_width, unsigned window_height, char *file_path) {

	FILE *fp = open_screenshot_file(file_path, "screenshot.jpg");
	if (fp == NULL) return 0;
	vector<unsigned char> buf((window_width+1)*(window_height+1)*3);
	read_pixels(window_width, window_height, buf);
	int const ret(write_jpeg_data(window_width, window_height, fp, &buf.front(), 1));
	return ret;
}
