// 3D World - Image I/O from texture_t
// by Frank Gennari
// 10/14/13
#include "targa.h"
#include "textures_3dw.h"

using std::string;
using std::cerr;

#ifdef ENABLE_JPEG
#include "jpeglib.h"
#endif

#ifdef ENABLE_PNG
#include "png.h"

void wrap_png_error(png_structp, png_const_charp) {	
	cerr << "Error reading PNG image file." << endl;
}
#endif

#ifdef ENABLE_TIFF
#include "tiffio.h"
#endif


std::string const texture_dir("textures");


FILE *open_texture_file(string const &filename) {

	FILE *fp = fopen((texture_dir + "/" + filename).c_str(), "rb");
	if (fp != NULL) return fp;

	// if not in the current directory, then look in the current directory
	fp = fopen(filename.c_str(), "rb");

	if (fp == NULL) {
		cout << "Error loading image " << filename << endl;
		exit(1);
	}
	return fp;
}


string get_file_extension(string const &filename, unsigned level, bool make_lower) {

	size_t const epos(filename.find_last_of('.'));
	size_t const spos[2] = {filename.find_last_of('\\'), filename.find_last_of('/')};
	size_t smax(0);
	string ext;

	for (unsigned i = 0; i < 2; ++i) {
		if (spos[i] != string::npos) smax = max(smax, spos[i]);
	}
	if (epos != string::npos && epos > smax) { // make sure the dot is after the last slash (part of the filename, not part of the path)
		ext = string(filename, epos+1, filename.length()-1);

		if (level > 0 && !ext.empty()) {
			string const fn2(string(filename, 0, epos));
			ext = get_file_extension(fn2, level-1, make_lower); // recursively strip off extensions
		}
	}
	unsigned const len((unsigned)ext.length());

	for (unsigned i = 0; i < len; ++i) { // convert upper case ext letters to lower case
		ext[i] = tolower(ext[i]);
	}
	return ext;
}


void texture_t::load(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale, bool ignore_word_alignment) {

	if (type > 0) { // generated texture
		alloc();
	}
	else {
		if (format == 7) { // auto
			// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel), 4: targa (*tga), 5: jpeg, 6: png, 7: auto, 8: tiff
			string const ext(get_file_extension(name, 0, 1));
		
			if (0) {}
			else if (ext == "raw") {
				format = ((ncolors == 4) ? 3 : 0);
			}
			else if (ext == "bmp") {
				format = 1;
			}
			else if (ext == "tga" || ext == "targa") {
				format = 4;
			}
			else if (ext == "jpg" || ext == "jpeg") {
				format = 5;
			}
			else if (ext == "png") {
				format = 6;
			}
			else if (ext == "tif" || ext == "tiff") {
				format = 8;
			}
			else {
				cerr << "Error: Unidentified image file format for autodetect: " << ext << " in filename " << name << endl;
				exit(1);
			}
		}
		unsigned want_alpha_channel(ncolors == 4);

		switch (format) {
		case 0: case 1: case 2: case 3: load_raw_bmp(index, allow_diff_width_height); break; // raw
		case 4: load_targa(index, allow_diff_width_height); break;
		case 5: load_jpeg (index, allow_diff_width_height); break;
		case 6: load_png  (index, allow_diff_width_height, allow_two_byte_grayscale); break;
		case 8: load_tiff (index, allow_diff_width_height, allow_two_byte_grayscale); break; // FIXME - tiff
		default:
			cerr << "Unsupported image format: " << format << endl;
			exit(1);
		}
		assert(is_allocated());
		if (invert_y) {do_invert_y();} // upside down

		if (want_alpha_channel && ncolors < 4) {
			add_alpha_channel();
		}
		else {
			try_compact_to_lum();
		}
		if (!ignore_word_alignment) {fix_word_alignment();}
	}
}


// http://paulbourke.net/dataformats/bmp/
struct bmp_header { // 14 bytes (may be padded to 16, but we only read 14)
   unsigned short int type;                 /* Magic identifier            */
   unsigned int size;                       /* File size in bytes          */
   unsigned short int reserved1, reserved2;
   unsigned int offset;                     /* Offset to image data, bytes */
};

struct bmp_infoheader { // 40 bytes
   unsigned int size;               /* Header size in bytes      */
   int width,height;                /* Width and height of image */
   unsigned short int planes;       /* Number of colour planes   */
   unsigned short int bits;         /* Bits per pixel            */
   unsigned int compression;        /* Compression type          */
   unsigned int imagesize;          /* Image size in bytes       */
   int xresolution,yresolution;     /* Pixels per meter          */
   unsigned int ncolours;           /* Number of colours         */
   unsigned int importantcolours;   /* Important colours         */
};


bool read_bmp_header(FILE *&fp, string const &fn, int &width, int &height, int &ncolors, bool allow_diff_width_height) {

	bmp_header header;
	bmp_infoheader infoheader;

	if (fread(&header, 14, 1, fp) != 1 || fread(&infoheader, 40, 1, fp) != 1) {
		cerr << "Error reading bitmap header/infoheader for file " << fn << endl;
		return 0;
	}
	int const img_ncolors(infoheader.bits >> 3);
	if (width   == 0 || allow_diff_width_height) {width  = infoheader.width;}
	if (height  == 0 || allow_diff_width_height) {height = infoheader.height;}
	if (ncolors == 0) {ncolors = img_ncolors;} // not reliable?

	if (width != infoheader.width || height != infoheader.height) { // check ncolors?
		cerr << "Error: bitmap file " << fn << " has incorrect width/height/ncolors: expected " << width << " " << height << " " << ncolors
			 << " but got " << infoheader.width << " " << infoheader.height << " " << img_ncolors << endl;
		return 0;
	}
	if (infoheader.compression != 0) {
		cerr << "Error: BMP compression mode " << infoheader.compression << " is not supported" << endl;
		return 0;
	}
	assert(width > 0 && height > 0 && ncolors > 0);

	if (ncolors == 1) { // read and discard color index table, and just use index values as grayscale values
		char data[1024];

		if (fread(data, 1024, 1, fp) != 1) {
			cerr << "Error reading bitmap data for file " << fn << endl;
			return 0;
		}
	}
	return 1;
}


// load an .RAW or .BMP file as a texture
// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel)
void texture_t::load_raw_bmp(int index, bool allow_diff_width_height) {

	assert(ncolors == 1 || ncolors == 3 || ncolors == 4);
	if (format == 3) assert(ncolors == 4);
	FILE *file(open_texture_file(name)); // open texture data
	assert(file != NULL);
	
	if (format == 1) { // BMP
		if (!read_bmp_header(file, name, width, height, ncolors, allow_diff_width_height)) {exit(1);}
	}
	unsigned const size(num_pixels()); // allocate buffer

	if (size == 0) {
		cerr << "Error loading texture image " << name << ": size not specified and not readable from this image format" << endl;
		exit(1);
	}
	assert(!is_allocated());
	alloc();

	// read texture data
	if (ncolors == 4 && format != 3) { // add alpha
		for (unsigned i = 0; i < size; ++i) {
			unsigned char buf[4];

			if (fread(buf, 3, 1, file) != 1) {
				cerr << "Error loading data from texture image " << name << endl;
				exit(1);
			}
			RGB_BLOCK_COPY((data+(i<<2)), buf);
		}
		auto_insert_alpha_channel(index);
	}
#if 0 // untested, enable if/when can be tested
	else if (format == 1 && (ncolors*width & 3)) { // not a multiple of 4 bytes - need to handle BMP padding
		unsigned const row_bytes(ncolors*width), stride(row_bytes + 4 - (row_bytes & 3)), nbytes(num_bytes());

		for (unsigned row = 0, pos = 0; row < (unsigned)height; ++row, pos += row_bytes) {
			assert(pos < nbytes);

			if (fread(data+pos, min(stride, nbytes-pos), 1, file) != 1) { // skip the padding bytes on the final scanline
				cerr << "Error loading data from texture image " << name << endl;
				exit(1);
			}
		}
	}
#endif
	else { // RGB or grayscale luminance
		if (fread(data, ncolors*size, 1, file) != 1) {
			cerr << "Error loading data from texture image " << name << endl;
			exit(1);
		}
	}
	if (format == 1 && (ncolors == 3 || ncolors == 4)) { // RGB/RGBA bitmap file
		for(unsigned i = 0; i < size; ++i) {
			swap(data[ncolors*i+0], data[ncolors*i+2]); // BGR[A] => RGB[A]
		}
	}
	fclose(file);
#if 0
	if (format == 0 && ncolors == 3) { // TESTING
		string fn(name);
		fn.erase(fn.begin()+fn.size()-4, fn.end());
		fn = texture_dir + "/gen/" + fn + ".jpg";
		cout << "Writing " << fn << endl;
		write_to_jpg(fn);
	}
#endif
}


void texture_t::load_targa(int index, bool allow_diff_width_height) {

	assert(!is_allocated());
	tga_image img;
	tga_result ret(tga_read(&img, (texture_dir + "/" + name).c_str())); // try textures directory
	//cout << "load texture" << name << endl;

	if (ret != TGA_NOERR) {
		ret = tga_read(&img, name.c_str()); // try current directory

		if (ret != TGA_NOERR) {
			cerr << "Error reading targa file " << name << ": " << tga_error(ret) << endl;
			exit(1);
		}
	}
	if (allow_diff_width_height || (width == 0 && height == 0)) {
		width  = img.width;
		height = img.height;
		assert(width > 0 && height > 0);
	}
	assert(img.width == width && img.height == height);
	alloc();
	//if (!tga_is_top_to_bottom(&img)) tga_flip_vert(&img);
	//if (tga_is_right_to_left(&img)) tga_flip_horiz(&img);

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			unsigned char const *const pixel(tga_find_pixel(&img, x, height-y-1)); // flip vert
			assert(pixel);
			unsigned char *d(data + ncolors*(x + y*width));
			tga_result const ret2(tga_unpack_pixel(pixel, img.pixel_depth, (ncolors>2 ? d+2 : 0), (ncolors>1 ? d+1 : 0), d, (ncolors>3 ? d+3 : 0)));
			assert(ret2 == TGA_NOERR);
		}
	}
	tga_free_buffers(&img);
}


void texture_t::load_jpeg(int index, bool allow_diff_width_height) {

#ifdef ENABLE_JPEG
	jpeg_decompress_struct cinfo;
	jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	FILE *fp(open_texture_file(name));

	if (fp == NULL) {
		cerr << "Error opening jpeg file " << name << " for read." << endl;
		exit(1);
	}
	jpeg_stdio_src(&cinfo, fp);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	if (allow_diff_width_height || (width == 0 && height == 0)) {
		width  = cinfo.output_width;
		height = cinfo.output_height;
		assert(width > 0 && height > 0);
	}
	assert(cinfo.output_width == width && cinfo.output_height == height);
	bool const want_alpha_channel(ncolors == 4 && cinfo.output_components == 3);
	ncolors = cinfo.output_components; // Note: can never be 4
	unsigned const scanline_size(ncolors*width);
	alloc();

	while (cinfo.output_scanline < cinfo.output_height) {
		JSAMPROW row_pointer[1] = {data + scanline_size*(cinfo.output_height - cinfo.output_scanline - 1)};
		jpeg_read_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(fp);

	if (want_alpha_channel) {
		add_alpha_channel();
		auto_insert_alpha_channel(index);
	}
#else
	cerr << "Error loading texture image file " << name << ": jpeg support has not been enabled." << endl;
	exit(1);
#endif
}


int write_jpeg_data(unsigned width, unsigned height, FILE *fp, unsigned char const *const data, bool invert_y) {

#ifdef ENABLE_JPEG
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	unsigned const step_size(3*width);
	JSAMPLE *rgb_row(new JSAMPLE[step_size]);
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width      = width;
	cinfo.image_height     = height;
	cinfo.input_components = 3;
	cinfo.in_color_space   = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_start_compress(&cinfo, TRUE);

	while (cinfo.next_scanline < cinfo.image_height) {
		JSAMPROW row_pointer[1];
		unsigned const yval(invert_y ? (height-cinfo.next_scanline-1) : cinfo.next_scanline);
		row_pointer[0] = (unsigned char *)&data[yval*step_size]; // cast away the const (we know the data won't be modified)
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	delete [] rgb_row;
	fclose(fp);
	return 1;
#else
  cerr << "Error: JPEG support is not enabled." << endl;
  return 0;
#endif
}


int texture_t::write_to_jpg(string const &fn) const {

	assert(ncolors == 3); // only supports RGB for now
	FILE *fp(fopen(fn.c_str(), "wb"));

	if (fp == NULL) {
		cerr << "Error opening jpg file " << fn << " for write." << endl;
		return 0;
	}
	return write_jpeg_data(width, height, fp, data, 0); // no invert
}


void texture_t::load_png(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale) {

#ifdef ENABLE_PNG
	FILE *fp(open_texture_file(name));

	if (fp == NULL) {
		cerr << "Error opening png file " << name << " for read." << endl;
		exit(1);
	}
	png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, (png_voidp)wrap_png_error, 0, 0);
	assert(png_ptr);
	png_infop info_ptr = png_create_info_struct(png_ptr);
	assert(info_ptr);
	png_infop end_info = png_create_info_struct(png_ptr);
	assert(end_info);
	png_init_io(png_ptr, fp);
	png_read_info(png_ptr, info_ptr);
	unsigned const w(png_get_image_width(png_ptr, info_ptr));
	unsigned const h(png_get_image_height(png_ptr, info_ptr));
	int const bit_depth(png_get_bit_depth(png_ptr, info_ptr));
	unsigned const png_ncolors(png_get_channels(png_ptr, info_ptr));

	if (allow_diff_width_height || (width == 0 && height == 0)) {
		width  = w;
		height = h;
		assert(width > 0 && height > 0);
	}
	assert(w == width && h == height);
	bool const want_alpha_channel(ncolors == 4 && png_ncolors == 3);
	ncolors = png_ncolors;

	if (allow_two_byte_grayscale && ncolors == 1 && bit_depth == 16) {
		ncolors        = 2; // change from 1 to 2 colors so that we can encode the high and low bytes into different channes to have 16-bit values
		is_16_bit_gray = 1;
		png_set_swap(png_ptr); // change big endian to little endian
	}
	else {
		if (bit_depth == 16) {png_set_strip_16(png_ptr);}
		if (bit_depth < 8)   {png_set_packing(png_ptr);}
	}
	vector<unsigned char *> rows(height);
	unsigned const scanline_size(ncolors*width);
	alloc();
	
	for (int i = 0; i < height; ++i) {
		rows[i] = data + (height - i - 1)*scanline_size;
	}
	png_read_image(png_ptr, &rows.front());
	png_read_end(png_ptr, end_info);
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	fclose(fp);

	if (want_alpha_channel) {
		add_alpha_channel();
		auto_insert_alpha_channel(index);
	}
#else
	cerr << "Error loading texture image file " << name << ": png support has not been enabled." << endl;
	exit(1);
#endif
}


void texture_t::load_tiff(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale) {

#ifdef ENABLE_TIFF
	TIFF* tif = TIFFOpen(name.c_str(), "r"); // FIXME: ignores texture directory
    
	if (tif == NULL) {
		cerr << "Error opening tiff file " << name << " for read." << endl;
		exit(1);
	}
	uint32 w(0), h(0);
	uint16 bit_depth(0), config(0);
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,    &w);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH,   &h);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bit_depth);
	TIFFGetField(tif, TIFFTAG_PLANARCONFIG,  &config);

	if (allow_diff_width_height || (width == 0 && height == 0)) {
		width  = w;
		height = h;
		assert(width > 0 && height > 0);
	}
	assert(w == width && h == height);
	
	if (allow_two_byte_grayscale && (ncolors == 0 || ncolors == 1) && bit_depth == 16) { // 16-bit grayscale
		ncolors        = 2; // change from 1 to 2 colors so that we can encode the high and low bytes into different channes to have 16-bit values
		is_16_bit_gray = 1;
		tmsize_t const sl_size(TIFFScanlineSize(tif));
		assert(sl_size == 2*width);
        tdata_t buf = _TIFFmalloc(sl_size);
		assert(buf);
		alloc();
		assert(config == PLANARCONFIG_CONTIG); // no support for PLANARCONFIG_SEPARATE, but could be added later

		for (int row = 0; row < height; row++) {
			TIFFReadScanline(tif, buf, row);

			for (int i = 0; i < sl_size; ++i) { // x-values
				data[sl_size*(height - row - 1) + i] = ((unsigned char const *)buf)[i]; // assumes little endian byte ordering, no swap required, may need to check this?
			}
		}
        _TIFFfree(buf);
	}
	else {
		if (ncolors == 0) {ncolors = 4;} // 3?
		uint32 *raster = (uint32 *)_TIFFmalloc(num_pixels() * sizeof(uint32));
		assert(raster != NULL);

		if (!TIFFReadRGBAImage(tif, width, height, raster, 0)) {
			cerr << "Error reading data from tiff file " << name << "." << endl;
			exit(1);
		}
		alloc();

		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				unsigned const ix(y*width + x);
				unsigned char const *d((unsigned char const *)(raster + ix));
				for (int i = 0; i < ncolors; ++i) {data[ncolors*ix+i] = d[i];} // not correct for lum+alpha textures?
			}
		}
		_TIFFfree(raster);
	}
	TIFFClose(tif);
#else
	cerr << "Error loading texture image file " << name << ": tiff support has not been enabled." << endl;
	exit(1);
#endif
}

