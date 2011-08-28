// http://www.opengl.org/resources/features/KilgardTechniques/LensFlare/
/* texture.c - by David Blythe, SGI */
/* load_luminace is a simplistic routine for reading an SGI .bw image file. */

#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <cassert>

unsigned const MAX_IR_CHARS = 256;

FILE *open_texture_file(std::string filename);


struct ImageRec {

	unsigned short imagic, type, dim, xsize, ysize, zsize;
	unsigned min, max, wasteBytes;
	char name[MAX_IR_CHARS];
	unsigned long colorMap;
	FILE *file;
	unsigned char *tmp;
	unsigned long rleEnd;
	unsigned *rowStart;
	int *rowSize;
};


static void ConvertShort(unsigned short *array, unsigned length) {

	unsigned short b1, b2;
	unsigned char *ptr;
	ptr = (unsigned char *) array;

	while (length--) {
		b1 = *ptr++;
		b2 = *ptr++;
		*array++ = (b1 << 8) | (b2);
	}
}


static void ConvertUint(unsigned *array, unsigned length) {

	unsigned char *ptr = (unsigned char *) array;

	while (length--) {
		unsigned const b1(*ptr++);
		unsigned const b2(*ptr++);
		unsigned const b3(*ptr++);
		unsigned const b4(*ptr++);
		*array++ = (b1 << 24) | (b2 << 16) | (b3 << 8) | (b4);
	}
}


static ImageRec *ImageOpen(std::string const &filename) {

	union {
		int testWord;
		char testByte[4];
	} endianTest;

	endianTest.testWord = 1;
	bool const swapFlag(endianTest.testByte[0] == 1);
	ImageRec *image(new ImageRec);
	image->file = open_texture_file(filename);
	assert(image->file != NULL);
	fread(image, 1, 12, image->file);
	if (swapFlag) ConvertShort(&image->imagic, 6);
	image->tmp = new unsigned char[image->xsize * 256];

	if ((image->type & 0xFF00) == 0x0100) {
		int const x(image->ysize * image->zsize * (int) sizeof(unsigned));
		image->rowStart = new unsigned[x];
		image->rowSize  = new int[x];
		image->rleEnd = 512 + (2 * x);
		fseek(image->file, 512, SEEK_SET);
		fread(image->rowStart, 1, x, image->file);
		fread(image->rowSize, 1, x, image->file);

		if (swapFlag) {
			ConvertUint(image->rowStart, x / (int) sizeof(unsigned));
			ConvertUint((unsigned *) image->rowSize, x / (int) sizeof(int));
		}
	}
	else {
	  image->rowStart = NULL;
	  image->rowSize  = NULL;
	}
	return image;
}


static void ImageClose(ImageRec * image) {

	fclose(image->file);
	delete [] image->rowStart;
	delete [] image->rowSize;
	delete [] image->tmp;
	delete image;
}


static void ImageGetRow(ImageRec * image, unsigned char *buf, int y, int z) {

	if ((image->type & 0xFF00) == 0x0100) {
		fseek(image->file, (long) image->rowStart[y + z * image->ysize], SEEK_SET);
		fread(image->tmp, 1, (unsigned) image->rowSize[y + z * image->ysize], image->file);
		unsigned char *iPtr = image->tmp;
		unsigned char *oPtr = buf;

		for (;;) {
			unsigned char pixel(*iPtr++);
			int count((int) (pixel & 0x7F));

			if (!count) return;
			if (pixel & 0x80) {
				while (count--) *oPtr++ = *iPtr++;
			}
			else {
				pixel = *iPtr++;
				while (count--) *oPtr++ = pixel;
			}
		}
	}
	else {
		fseek(image->file, 512 + (y * image->xsize) + (z * image->xsize * image->ysize), SEEK_SET);
		fread(buf, 1, image->xsize, image->file);
	}
}


unsigned char *load_luminance(std::string const &filename, int *width, int *height, int *components) {

	ImageRec *image = ImageOpen(filename.c_str());
	if (!image) return NULL;
	
	if (image->zsize != 1) {
		ImageClose(image);
		return NULL;
	}
	*width      = image->xsize;
	*height     = image->ysize;
	*components = image->zsize;
	unsigned char *base = new unsigned char[image->xsize * image->ysize];
	if (base == NULL) return NULL;
	unsigned char *lptr = base;

	for (int y = 0; y < image->ysize; y++) {
		ImageGetRow(image, lptr, y, 0);
		lptr += image->xsize;
	}
	ImageClose(image);
	return base;
}

