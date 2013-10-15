// OpenGL Pattern Matcher Graphical Display - Image Export
// Must be compiled as C, not C++, to enable JPEG support
// by Frank Gennari 6/10/02

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gl_includes.h"


int write_jpeg_data(unsigned window_width, unsigned window_height, FILE *fp, unsigned char *data);


FILE *open_screenshot_file(char *file_path, char *basename) {

	FILE *fp;
	if (basename == NULL) return NULL;

	if (file_path != NULL) {
		char *filename = new char[strlen(file_path)+20];
		sprintf(filename, "%s%s", file_path, basename);
		fp = fopen(filename, "wb");
		delete [] filename;
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


int screenshot(unsigned window_width, unsigned window_height, char *file_path) {

	FILE *fp = open_screenshot_file(file_path, "screenshot.raw");
	if (fp == NULL) return 0;
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	unsigned const bufSize(((window_width)*(window_height)*3));
	unsigned char *buf(new unsigned char[bufSize]);
	glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	for (unsigned i = 0; i < window_height; ++i) {
		unsigned const offset((window_height-i-1)*window_width);

		for (unsigned j = 0; j < window_width; ++j) {
			fprintf(fp, "%c%c%c", buf[3*(j+offset)], buf[3*(j+offset)+1], buf[3*(j+offset)+2]);
		}
	}
	delete [] buf;
	fclose(fp);
	return 1;
}


int write_jpeg(unsigned window_width, unsigned window_height, char *file_path) {

	FILE *fp = open_screenshot_file(file_path, "screenshot.jpg");
	if (fp == NULL) return 0;
	unsigned char *buf(new unsigned char[(window_width+1)*(window_height+1)*3]);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, buf);
	int const ret(write_jpeg_data(window_width, window_height, fp, buf));
	delete [] buf;
	return ret;
}
