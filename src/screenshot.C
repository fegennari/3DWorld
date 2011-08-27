// OpenGL Pattern Matcher Graphical Display - Image Export
// Must be compiled as C, not C++, to enable JPEG support
// by Frank Gennari 6/10/02

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gl_includes.h"


FILE *open_screenshot_file(char *file_path, char *basename) {

	FILE *fp;
	if (basename == NULL) return NULL;

	if (file_path != NULL) {
		char *filename = (char *)malloc((strlen(file_path)+20)*sizeof(char));
		sprintf(filename, "%s%s", file_path, basename);
		fp = fopen(filename, "wb");
		free(filename);
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

	unsigned i, j, offset, bufSize;
	unsigned char *buf;
	FILE *fp = open_screenshot_file(file_path, "screenshot.raw");
	if (fp == NULL) return 0;
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	bufSize = ((window_width)*(window_height)*3*sizeof(unsigned char));
	buf     = (unsigned char *)malloc(bufSize);
	glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	for (i = 0; i < window_height; ++i) {
		offset = (window_height-i-1)*window_width;
		for (j = 0; j < window_width; ++j) {
			fprintf(fp, "%c%c%c", buf[3*(j+offset)], buf[3*(j+offset)+1], buf[3*(j+offset)+2]);
		}
	}
	free(buf);
	fclose(fp);
	return 1;
}


#ifdef ENABLE_JPEG

#include "../jpeg-6b/jpeglib.h"

int write_jpeg(unsigned window_width, unsigned window_height, char *file_path) {

	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	unsigned step_size;
	unsigned char *buf;
	JSAMPROW row_pointer[1];
	JSAMPLE *rgb_row;
	FILE *fp = open_screenshot_file(file_path, "screenshot.jpg");
	if (fp == NULL) return 0;
	buf = (unsigned char *)malloc((window_width+1)*(window_height+1)*3*sizeof(unsigned char));
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	step_size = 3*window_width;
	rgb_row   = (JSAMPLE *)malloc(sizeof(JSAMPLE)*step_size);
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width      = window_width;
	cinfo.image_height     = window_height;
	cinfo.input_components = 3;
	cinfo.in_color_space   = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_start_compress(&cinfo, TRUE);

	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &buf[(window_height-cinfo.next_scanline-1)*step_size];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	free(rgb_row);
	free(buf);
	fclose(fp);
	return 1;
}

#else

int write_jpeg(unsigned window_width, unsigned window_height, char *file_path) {

  printf("Error: JPEG support is not enabled.\n");
  return 0;
}

#endif
