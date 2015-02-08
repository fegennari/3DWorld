// 3D World - Image I/O from texture_t
// by Frank Gennari
// 10/14/13
#include "targa.h"
#include "textures_3dw.h"

using std::string;
using std::cerr;

#ifdef ENABLE_JPEG
#define INT32 prev_INT32 // fix conflicting typedef used in freeglut
#include "jpeglib.h"
#undef prev_INT32
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

string append_texture_dir(string const &filename) {return (texture_dir + "/" + filename);}


FILE *open_texture_file(string const &filename) {

	FILE *fp = fopen(append_texture_dir(filename).c_str(), "rb");
	if (fp != NULL) return fp;

	// if not in the current directory, then look in the current directory
	fp = fopen(filename.c_str(), "rb");

	if (fp == NULL) {
		cerr << "Error loading image " << filename << endl;
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
	for (unsigned i = 0; i < len; ++i) {ext[i] = tolower(ext[i]);} // convert upper case ext letters to lower case
	return ext;
}


void texture_t::load(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale, bool ignore_word_alignment) {

	if (type > 0) { // generated texture
		alloc();
	}
	else {
		if (format == 7) { // auto
			// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel), 4: targa (*tga), 5: jpeg, 6: png, 7: auto, 8: tiff, 10: DDS
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
			else if (ext == "dds") {
				format = 10;
			}
			else {
				cerr << "Error: Unidentified image file format for autodetect: " << ext << " in filename " << name << endl;
				exit(1);
			}
		}
		unsigned const want_alpha_channel(ncolors == 4), want_luminance(ncolors == 1);

		switch (format) {
		case 0: case 1: case 2: case 3: load_raw_bmp(index, allow_diff_width_height, allow_two_byte_grayscale); break; // raw
		case 4: load_targa(index, allow_diff_width_height); break;
		case 5: load_jpeg (index, allow_diff_width_height); break;
		case 6: load_png  (index, allow_diff_width_height, allow_two_byte_grayscale); break;
		case 8: load_tiff (index, allow_diff_width_height, allow_two_byte_grayscale); break;
		case 10: load_dds (index); break;
		default:
			cerr << "Unsupported image format: " << format << endl;
			exit(1);
		}
		// defer this check until we actually need to access the data, in case we want to actually do the load on the fly later
		//assert(is_allocated());
		assert(is_loaded());
		if (invert_y) {do_invert_y();} // upside down

		if (want_alpha_channel && ncolors < 4) {
			add_alpha_channel();
		}
		else if (want_luminance && ncolors == 3) {
			try_compact_to_lum();
		}
		if (!ignore_word_alignment) {fix_word_alignment();}

		if (invert_alpha) {
			if (ncolors == 1 || ncolors == 3) { // if 3 colors, assume all are duplicate alpha channels
				assert(!is_16_bit_gray);
				unsigned const nbytes(num_bytes());
				for (unsigned i = 0; i < nbytes; ++i) {data[i] = (255 - data[i]);}
			}
			else {
				assert(ncolors == 4);
				unsigned const npixels(num_pixels());
				for (unsigned i = 0; i < npixels; ++i) {data[4*i+3] = (255 - data[4*i+3]);}
			}
		}
	}
#if 0
	if (name.size() > 4 && name.front() != '@') {
		string fn(name);
		fn.erase(fn.begin()+fn.size()-4, fn.end());
		fn += ".bmp";
		fn  = texture_dir + "/gen/";
		cout << "Writing " << fn << endl;
		write_to_bmp(fn);
	}
#endif
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


bool read_bmp_header(FILE *&fp, string const &fn, int &width, int &height, int &ncolors, bool allow_diff_width_height, bool allow_two_byte_grayscale, bool &is_16_bit_gray) {

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

	if (allow_two_byte_grayscale && ncolors == 1 && img_ncolors == 2) { // Note: not officially part of the BMP spec, but we allow it since we write in this format
		ncolors        = 2; // change from 1 to 2 colors so that we can encode the high and low bytes into different channes to have 16-bit values
		is_16_bit_gray = 1;
	}
	if (ncolors == 1) { // read and discard color index table, and just use index values as grayscale values
		char color_table[1024];

		if (fread(color_table, 1024, 1, fp) != 1) {
			cerr << "Error reading bitmap color table for file " << fn << endl;
			return 0;
		}
	}
	return 1;
}


// load an .RAW or .BMP file as a texture
// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel)
void texture_t::load_raw_bmp(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale) {

	assert(ncolors == 1 || ncolors == 3 || ncolors == 4);
	if (format == 3) assert(ncolors == 4);
	FILE *file(open_texture_file(name)); // open texture data
	assert(file != NULL);
	
	if (format == 1) { // BMP
		if (!read_bmp_header(file, name, width, height, ncolors, allow_diff_width_height, allow_two_byte_grayscale, is_16_bit_gray)) {exit(1);}
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
	if (format == 1) {maybe_swap_rb(data);}
	fclose(file);
}


void maybe_swap_rb(unsigned char *ptr, unsigned num_pixels, unsigned ncolors) {

	assert(ptr != NULL);
	if (ncolors != 3 && ncolors != 4) return;

	for(unsigned i = 0; i < num_pixels; ++i) {
		swap(ptr[ncolors*i+0], ptr[ncolors*i+2]); // BGR[A] => RGB[A]
	}
}


void texture_t::maybe_swap_rb(unsigned char *ptr) const { // ptr is assumed to be of size num_bytes()
	::maybe_swap_rb(ptr, num_pixels(), ncolors);
}


bool write_rgb_bmp_image(FILE *fp, string const &fn, unsigned char *data, unsigned width, unsigned height, unsigned ncolors) {

	maybe_swap_rb(data, width*height, ncolors); // Note: data not const because of this line
	unsigned char pad[4] = {0};
	unsigned const row_sz(width*ncolors), row_sz_mod(row_sz&3), row_pad(row_sz_mod ? 4-row_sz_mod : 0);
	bmp_header header = {0};
	header.type = 19778; // bitmap
	//header.size = 54 + (row_sz + row_pad)*height + ((ncolors == 1) ? 1024 : 0); // optional
	//header.offset = 54; // optional
	bmp_infoheader infoheader = {0};
	infoheader.width  = width;
	infoheader.height = height;
	infoheader.bits   = ncolors << 3;
	infoheader.planes = 1;
	infoheader.size   = 40;
	infoheader.xresolution = infoheader.yresolution = 1200; // arbitrary nonzero

	if (fwrite(&header, 14, 1, fp) != 1 || fwrite(&infoheader, 40, 1, fp) != 1) {
		cerr << "Error writing bitmap header/infoheader for file " << fn << endl;
		return 0;
	}
	if (ncolors == 1) { // add color index table
		char color_table[1024] = {0};
		for (unsigned i = 0; i < 256; ++i) {UNROLL_3X(color_table[(i<<2)+i_] = i;)}

		if (fwrite(color_table, 1024, 1, fp) != 1) {
			cerr << "Error writing bitmap color table for file " << fn << endl;
			return 0;
		}
	}
	for (unsigned i = 0; i < height; ++i) { // write one scanline at a time (could invert y if needed)
		if (fwrite(data, 1, row_sz, fp) != row_sz) { // row image data
			cerr << "Error writing bitmap data for file " << fn << endl;
			return 0;
		}
		if (row_pad > 0 && fwrite(pad, 1, row_pad, fp) != row_pad) { // maybe add padding
			cerr << "Error writing bitmap row padding for file " << fn << endl;
			return 0;
		}
		data += row_sz;
	}
	return 1;
}


int texture_t::write_to_bmp(string const &fn) const {

	FILE *fp(fopen(fn.c_str(), "wb"));

	if (fp == NULL) {
		cerr << "Error opening bmp file " << fn << " for write." << endl;
		return 0;
	}
	vector<unsigned char> data_swap_rb(data, data+num_bytes());
	bool const ret(write_rgb_bmp_image(fp, fn, &data_swap_rb.front(), width, height, ncolors));
	fclose(fp);
	return ret;
}


void texture_t::load_targa(int index, bool allow_diff_width_height) {

	assert(!is_allocated());
	tga_image img;
	tga_result ret(tga_read(&img, append_texture_dir(name).c_str())); // try textures directory
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
	if (img.width != width || img.height != height) {
		cerr << "Incorrect image size for " << name << ": expected " << width << "x" << height << ", got " << img.width << "x" << img.height << endl;
		exit(1);
	}
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
	if (cinfo.output_width != width || cinfo.output_height != height) {
		cerr << "Incorrect image size for " << name << ": expected " << width << "x" << height << ", got " << cinfo.output_width << "x" << cinfo.output_height << endl;
		exit(1);
	}
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
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width      = width;
	cinfo.image_height     = height;
	cinfo.input_components = 3;
	cinfo.in_color_space   = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality (&cinfo, 100, 0); // set highest quality
	jpeg_start_compress(&cinfo, TRUE);

	while (cinfo.next_scanline < cinfo.image_height) {
		JSAMPROW row_pointer[1];
		unsigned const yval(invert_y ? (height-cinfo.next_scanline-1) : cinfo.next_scanline);
		row_pointer[0] = (unsigned char *)&data[yval*step_size]; // cast away the const (we know the data won't be modified)
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	fclose(fp);
	return 1;
#else
  cerr << "Error: JPEG writing support is not enabled." << endl;
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
	return write_jpeg_data(width, height, fp, data, 0); // no invert, fclose(fp) is called within this function
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
	if (w != width || h != height) {
		cerr << "Incorrect image size for " << name << ": expected " << width << "x" << height << ", got " << w << "x" << h << endl;
		exit(1);
	}
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


int texture_t::write_to_png(string const &fn) const {

#ifdef ENABLE_PNG
	FILE *fp(fopen(fn.c_str(), "wb"));

	if (fp == NULL) {
		cerr << "Error opening png file " << fn << " for write." << endl;
		return 0;
	}
	// Initialize write structure
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	assert(png_ptr != NULL);

	// Initialize info structure
	png_infop info_ptr = png_create_info_struct(png_ptr);
	assert(info_ptr != NULL);
	png_init_io(png_ptr, fp);

	// Write header
	int color_type, bit_depth;

	if (is_16_bit_gray) {
		color_type = PNG_COLOR_TYPE_GRAY;
		bit_depth  = 16;
	}
	else {
		switch (ncolors) {
		case 1: color_type = PNG_COLOR_TYPE_GRAY; break;
		case 2: color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
		case 3: color_type = PNG_COLOR_TYPE_RGB; break;
		case 4: color_type = PNG_COLOR_TYPE_RGB_ALPHA; break;
		default: assert(0);
		}
		bit_depth = 8;
	}
	png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	png_write_info(png_ptr, info_ptr);

	// Write image data
	for (int y = 0; y < height; y++) {png_write_row(png_ptr, (data + y*width*ncolors));}

	// End write
	png_write_end(png_ptr, NULL);
	png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	fclose(fp);
	return 1;
#else
	cerr << "Error: PNG writing support is not enabled." << endl;
	return 0;
#endif
}


void texture_t::load_tiff(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale) {

#ifdef ENABLE_TIFF
	TIFF* tif = TIFFOpen(append_texture_dir(name).c_str(), "r"); // first try texture directory
	if (tif == NULL) {tif = TIFFOpen(name.c_str(), "r");} // not found, try current directory

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


#ifdef ENABLE_DDS
#include <gli.hpp>
#endif

void texture_t::load_dds(int index) {
	
#ifdef ENABLE_DDS
	defer_load_type = DEFER_TYPE_DDS;
#else
	cerr << "Error loading texture image file " << name << ": DDS support has not been enabled." << endl;
	exit(1);
#endif
}


void texture_t::deferred_load_and_bind() {

	defer_load();

	switch (defer_load_type) {
#ifdef ENABLE_DDS
	case DEFER_TYPE_DDS:
		{
			//cout << "Loading DDS image " << name << endl;
			gli::texture2D Texture(gli::load_dds(name.c_str()));
			bool const compressed(gli::is_compressed(Texture.format()));
			// here we assume the texture is upside down and flip it, if it's uncompressed and flippable
			if (!compressed) {Texture = flip(Texture);}
			assert(!Texture.empty());
			width   = Texture.dimensions().x;
			height  = Texture.dimensions().y;
			ncolors = component_count(Texture.format());
			assert(width > 0 && height > 0);
			glBindTexture(GL_TEXTURE_2D, tid);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, GLint(Texture.levels() - 1));
			glTexStorage2D(GL_TEXTURE_2D, Texture.levels(), gli::internal_format(Texture.format()), width, height);

			for (gli::texture2D::size_type Level = 0; Level < Texture.levels(); ++Level) {
				if (compressed) {
					glCompressedTexSubImage2D(GL_TEXTURE_2D, Level, 0, 0, Texture[Level].dimensions().x, Texture[Level].dimensions().y,
						gli::internal_format(Texture.format()), Texture[Level].size(), Texture[Level].data());
				}
				else {
					glTexSubImage2D(GL_TEXTURE_2D, Level, 0, 0, Texture[Level].dimensions().x, Texture[Level].dimensions().y,
						gli::external_format(Texture.format()), gli::type_format(Texture.format()), Texture[Level].data());
				}
			} 
		}
		break;
	default:
		cerr << "Unhandled texture defer type " << defer_load_type << endl;
		exit(1);
#endif
	}
}
