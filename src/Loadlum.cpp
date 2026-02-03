// http://www.opengl.org/resources/features/KilgardTechniques/LensFlare/
/* texture.c - by David Blythe, SGI */
/* load_luminace is a simplistic routine for reading an SGI .bw image file. */

#include "globals.h"
#include <stdio.h>
#include <string>
#include <cassert>

unsigned const MAX_IR_CHARS = 256;

FILE *open_texture_file(std::string const &filename);
void checked_fseek_to(FILE *fp, long fpos);
void checked_fclose(FILE *fp);


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

	unsigned char *ptr = (unsigned char *) array;

	while (length--) {
		unsigned short b1 = *ptr++;
		unsigned short b2 = *ptr++;
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
	int num_read(fread(image, 1, 12, image->file));
	assert(num_read == 12);
	if (swapFlag) ConvertShort(&image->imagic, 6);
	image->tmp = new unsigned char[image->xsize * 256];

	if ((image->type & 0xFF00) == 0x0100) {
		int const x(image->ysize * image->zsize * (int) sizeof(unsigned));
		image->rowStart = new unsigned[x];
		image->rowSize  = new int[x];
		image->rleEnd = 512 + (2 * x);
		checked_fseek_to(image->file, 512);
		num_read = fread(image->rowStart, 1, x, image->file);
		assert(num_read == x);
		num_read = fread(image->rowSize, 1, x, image->file);
		assert(num_read == x);

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

	checked_fclose(image->file);
	delete [] image->rowStart;
	delete [] image->rowSize;
	delete [] image->tmp;
	delete image;
}


static void ImageGetRow(ImageRec * image, unsigned char *buf, int y, int z) {

	if ((image->type & 0xFF00) == 0x0100) {
		checked_fseek_to(image->file, (long)image->rowStart[y + z * image->ysize]);
		unsigned const row_sz(image->rowSize[y + z * image->ysize]);
		size_t const num_read(fread(image->tmp, 1, row_sz, image->file));
		assert(num_read == row_sz);
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
		checked_fseek_to(image->file, (512 + (y * image->xsize) + (z * image->xsize * image->ysize)));
		size_t const num_read(fread(buf, 1, image->xsize, image->file));
		assert(num_read == image->xsize);
	}
}


unsigned char *load_luminance(std::string const &filename, int *width, int *height, int *components) {

	ImageRec *image = ImageOpen(filename);
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

