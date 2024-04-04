// 3D World - Image I/O from texture_t
// by Frank Gennari
// 10/14/13
#include "targa.h"
#include "textures.h"
//#include "profiler.h"
#include <fstream> // for filebuf

using namespace std;

// requires ENABLE_STB_IMAGE; no support for 16-bit PNG heightmaps, TIFF, or DDS (will default to other readers in those cases)
bool const USE_STB_IMAGE_LOAD = 0;

// Note: stb_image generally works as a replacement for ligpng/ligjpeg/libtiff, but it doesn't support 16-bit formats, so it can't fully replace libpng; similar runtime
#ifdef ENABLE_STB_IMAGE
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
bool stb_image_enabled() {return 1;}
#else
bool stb_image_enabled() {return 0;}
#endif

#ifdef ENABLE_PNG
#include "png.h"

void wrap_png_error(png_structp, png_const_charp) {cerr << "Error reading PNG image file." << endl;}
#endif

#ifdef ENABLE_DDS
#include <gli/gli.hpp> // Note: must be after tiffio.h include due to conflicting uint32_t typedef
#endif

#ifdef ENABLE_TIFF
#include "tiffio.h"
#endif


string const texture_dir("textures");

string append_texture_dir(string const &filename) {return (texture_dir + "/" + filename);}


size_t get_last_slash_pos(string const &filename) {
	size_t const spos[2] = {filename.find_last_of('\\'), filename.find_last_of('/')};
	size_t smax(0);

	for (unsigned i = 0; i < 2; ++i) {
		if (spos[i] != string::npos) {smax = max(smax, spos[i]);}
	}
	return smax;
}
string get_base_filename(string const &filename) {
	size_t const smax(get_last_slash_pos(filename));
	if (smax == 0) return filename; // no slash
	return string(filename, smax+1, filename.length()-1);
}
string get_file_extension(string const &filename, unsigned level, bool make_lower) {
	size_t const epos(filename.find_last_of('.'));
	size_t const smax(get_last_slash_pos(filename));
	string ext;

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

void checked_fclose(FILE *fp) {
	if (fclose(fp) != 0) {
		perror("Error: fclose() call failed");
		exit(1); // how fatal should this be?
	}
}

FILE *open_texture_file_no_check(string const &filename) {
	FILE *fp = fopen(append_texture_dir(filename).c_str(), "rb");
	if (fp != nullptr) return fp; // found
	// if not in the texture directory, look in the current directory
	return fopen(filename.c_str(), "rb");
}
FILE *open_texture_file(string const &filename) {
	FILE *fp(open_texture_file_no_check(filename));

	if (fp == nullptr) {
		cerr << endl << "Error loading image " << filename << endl;
		exit(1);
	}
	return fp;
}
bool check_texture_file_exists(string const &filename) {
	FILE *fp(open_texture_file_no_check(filename));
	if (fp == nullptr) return 0;
	checked_fclose(fp);
	return 1;
}

void texture_t::load(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale, bool ignore_word_alignment) {

	//timer_t timer("Load Texture " + name);

	if (type > 0) { // generated texture
		alloc();
		memset(data, 0, num_bytes()); // zero the values to make sure we don't accidentally use it uninitialized before the texture is generated
	}
	else {
		if (format == IMG_FMT_AUTO) { // auto
			string const ext(get_file_extension(name, 0, 1));
		
			if (0) {}
			else if (ext == "raw") {format = ((ncolors == 4) ? IMG_FMT_RAW_RGBA : IMG_FMT_RAW_RGB);}
			else if (ext == "bmp") {format = IMG_FMT_BMP;}
			else if (ext == "tga" || ext == "targa") {format = IMG_FMT_TGA;}
			else if (ext == "jpg" || ext == "jpeg" ) {format = IMG_FMT_JPG;}
			else if (ext == "png") {format = IMG_FMT_PNG;}
			else if (ext == "tif" || ext == "tiff") {format = IMG_FMT_TIFF;}
			else if (ext == "dds") {format = IMG_FMT_DDS;}
			else if (ext == "ppm") {format = IMG_FMT_PPM;}
			else if (ext == "hdr") {
				cerr << "Error: HDR texture format is not yet supported: " << name << endl;
				exit(1);
			}
			else { // unsupported native image format
#ifdef ENABLE_STB_IMAGE // try loading with stb_image; should work for GIF, etc.
				if (load_stb_image(index, allow_diff_width_height)) {format = IMG_FMT_OTHER;}
				else
#endif
				{
					cerr << "Error: Unidentified image file format for autodetect: " << ext << " in filename " << name << endl;
					exit(1);
				}
			}
		}
		unsigned const want_alpha_channel(ncolors == 4), want_luminance(ncolors == 1);
		//highres_timer_t timer("Load " + get_file_extension(name, 0, 1)); // 0.1s bmp, 6.3s jpeg, 4.4s png, 0.3s tga, 0.2s tiff

		// attempt to load image with STB if it's enabled; applies to BMP, TGA, JPEG, and PNG; 2-byte grayscale images are not supported
		if (USE_STB_IMAGE_LOAD && (format == IMG_FMT_BMP || format == IMG_FMT_TGA || format == IMG_FMT_JPG || format == IMG_FMT_PNG || format == IMG_FMT_OTHER) &&
			!allow_two_byte_grayscale && load_stb_image(index, allow_diff_width_height))
		{
			// loaded with stb_image
		}
		else {
			switch (format) {
			case IMG_FMT_RAW_RGB: case IMG_FMT_BMP: case IMG_FMT_RAW_INVY: case IMG_FMT_RAW_RGBA: load_raw_bmp(index, allow_diff_width_height, allow_two_byte_grayscale); break; // raw or BMP
			case IMG_FMT_TGA: load_targa(index, allow_diff_width_height); break;
			case IMG_FMT_JPG: load_jpeg (index, allow_diff_width_height); break;
			case IMG_FMT_PNG: load_png  (index, allow_diff_width_height, allow_two_byte_grayscale); break;
			case IMG_FMT_TIFF: load_tiff (index, allow_diff_width_height, allow_two_byte_grayscale); break;
			case IMG_FMT_DDS: load_dds (index); break;
			case IMG_FMT_PPM: load_ppm (index, allow_diff_width_height); break;
			case IMG_FMT_OTHER: break; // already loaded by stb_image
			default:
				cerr << "Unsupported image format: " << int(format) << endl;
				exit(1);
			}
		}
		//timer.end();
		// defer this check until we actually need to access the data, in case we want to actually do the load on the fly later
		//assert(is_allocated());
		assert(is_loaded());
		if (invert_y && format != IMG_FMT_DDS) {do_invert_y();} // upside down (not DDS)
		if (want_alpha_channel && ncolors < 4) {add_alpha_channel();}
		else if (want_luminance && ncolors == 3) {try_compact_to_lum();}
		//if (want_alpha_channel) {fill_transparent_with_avg_color();}
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
	} // end non-generated texture case
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

void texture_t::set_image_size(int w, int h, bool allow_diff_width_height) {
	if (allow_diff_width_height || (width == 0 && height == 0)) {
		width  = w;
		height = h;
		assert(width > 0 && height > 0);
	}
	else if (w != width || h != height) {
		cerr << "Incorrect image size for " << name << ": expected " << width << "x" << height << ", got " << w << "x" << h << endl;
		exit(1);
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
		is_16_bit_gray = 1; // AKA set_16_bit_grayscale()
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
		for (size_t i = 0; i < size; ++i) {
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
	checked_fclose(file);
}


void maybe_swap_rb(unsigned char *ptr, unsigned num_pixels, unsigned ncolors) {
	assert(ptr != NULL);
	if (ncolors != 3 && ncolors != 4) return;
	for(unsigned i = 0; i < num_pixels; ++i) {swap(ptr[ncolors*i+0], ptr[ncolors*i+2]);} // BGR[A] => RGB[A]
}


void texture_t::maybe_swap_rb(unsigned char *ptr) const { // ptr is assumed to be of size num_bytes()
	::maybe_swap_rb(ptr, num_pixels(), ncolors);
}


bool write_rgb_bmp_image(string const &fn, unsigned char *data, unsigned width, unsigned height, unsigned ncolors) {

	FILE *fp(fopen(fn.c_str(), "wb"));

	if (fp == NULL) {
		cerr << "Error opening BMP file " << fn << " for write." << endl;
		return 0;
	}
	maybe_swap_rb(data, width*height, ncolors); // Note: data not const because of this line
	unsigned char pad[4] = {0};
	unsigned const row_sz(width*ncolors), row_sz_mod(row_sz&3), row_pad(row_sz_mod ? 4-row_sz_mod : 0);
	bmp_header header = {};
	header.type = 19778; // bitmap
	//header.size = 54 + (row_sz + row_pad)*height + ((ncolors == 1) ? 1024 : 0); // optional
	//header.offset = 54; // optional
	bmp_infoheader infoheader = {};
	infoheader.width  = width;
	infoheader.height = height;
	infoheader.bits   = ncolors << 3;
	infoheader.planes = 1;
	infoheader.size   = 40;
	infoheader.xresolution = infoheader.yresolution = 1200; // arbitrary nonzero

	if (fwrite(&header, 14, 1, fp) != 1 || fwrite(&infoheader, 40, 1, fp) != 1) {
		cerr << "Error writing BMP header/infoheader for file " << fn << endl;
		checked_fclose(fp);
		return 0;
	}
	if (ncolors == 1) { // add color index table
		char color_table[1024] = {0};
		for (unsigned i = 0; i < 256; ++i) {UNROLL_3X(color_table[(i<<2)+i_] = i;)}

		if (fwrite(color_table, 1024, 1, fp) != 1) {
			cerr << "Error writing BMP color table for file " << fn << endl;
			checked_fclose(fp);
			return 0;
		}
	}
	for (unsigned i = 0; i < height; ++i) { // write one scanline at a time (could invert y if needed)
		if (fwrite(data, 1, row_sz, fp) != row_sz) { // row image data
			cerr << "Error writing BMP data for file " << fn << endl;
			checked_fclose(fp);
			return 0;
		}
		if (row_pad > 0 && fwrite(pad, 1, row_pad, fp) != row_pad) { // maybe add padding
			cerr << "Error writing BMP row padding for file " << fn << endl;
			checked_fclose(fp);
			return 0;
		}
		data += row_sz;
	}
	checked_fclose(fp);
	return 1;
}


int texture_t::write_to_bmp(string const &fn) const {
	vector<unsigned char> data_swap_rb(data, data+num_bytes());
	return write_rgb_bmp_image(fn, &data_swap_rb.front(), width, height, ncolors);
}


void texture_t::load_targa(int index, bool allow_diff_width_height) {

	assert(!is_allocated());
	tga_image img;
	tga_result ret(tga_read(&img, append_texture_dir(name).c_str())); // try textures directory

	if (ret != TGA_NOERR) {
		ret = tga_read(&img, name.c_str()); // try current directory

		if (ret != TGA_NOERR) {
			cerr << "Error reading targa file " << name << ": " << tga_error(ret) << endl;
			exit(1);
		}
	}
	set_image_size(img.width, img.height, allow_diff_width_height);
	// overwrite ncolors; if caller set ncolors=4 then assume it has an alpha channel
	if      (img.image_type == TGA_IMAGE_TYPE_BGR  || img.image_type == TGA_IMAGE_TYPE_BGR_RLE ) {ncolors = ((ncolors == 4) ? 4 : 3);} // RGB/RGBA
	else if (img.image_type == TGA_IMAGE_TYPE_MONO || img.image_type == TGA_IMAGE_TYPE_MONO_RLE) {ncolors = 1;} // Red
	// else leave ncolors unchanged from what was specified by the caller before load_targa()
	alloc();

#pragma omp parallel for schedule(static)
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
	//timer_t timer("Load Jpeg"); // 187 calls, 7063ms total
	if (!load_stb_image(index, allow_diff_width_height)) { // allow_two_byte_grayscale=0
		cerr << "Error loading JPEG texture image file " << name << endl;
		exit(1);
	}
}

class stb_buffered_writer {
	FILE *fp=nullptr;
	vector<unsigned char> buffer;
public:
	~stb_buffered_writer() {fclose(fp);}

	bool init(string const &fn) {
		fp = fopen(fn.c_str(), "wb");
		if (fp == NULL) {cerr << "Error opening image file " << fn << " for write." << endl; return 0;}
		return 1;
	}
	void add_data(void *data, int size) {buffer.insert(buffer.end(), (char *)data, (char *)data+size);}
	
	void flush() {
		assert(fp);
		if (buffer.empty()) return;
		fwrite(buffer.data(), 1, buffer.size(), fp);
		buffer.clear();
	}
};
void stb_buffered_writer_func(void *context, void *data, int size) {
	assert(context);
	((stb_buffered_writer *)context)->add_data(data, size);
}

int write_jpeg_data(string const &fn, unsigned char const *const data, unsigned width, unsigned height, bool invert_y) {
#ifdef ENABLE_STB_IMAGE
	stbi_flip_vertically_on_write(1);
	//timer_t timer("JPEG Write");
	// use a buffered writer because it's > 2x faster (100ms vs. 225ms for 1920x1080 image)
	//return stbi_write_jpg(fn.c_str(), width, height, 3, data, 95); // ncomp=3, quality=95
	stb_buffered_writer writer;
	if (!writer.init(fn)) return 0;
	if (!stbi_write_jpg_to_func(stb_buffered_writer_func, (void *)&writer, width, height, 3, data, 95)) return 0; // ncomp=3, quality=95
	writer.flush(); // flush once at the end, which writes the entire image as one fwrite() call
	return 1;
#else
	cerr << "Error: JPEG writing support is not enabled." << endl;
	return 0;
#endif
}
int texture_t::write_to_jpg(string const &fn) const {
	assert(ncolors == 3); // only supports RGB for now
	return write_jpeg_data(fn, data, width, height, 0); // no invert
}


void texture_t::set_16_bit_grayscale() {
	ncolors        = 2; // change from 1 to 2 colors so that we can encode the high and low bytes into different channes to have 16-bit values
	is_16_bit_gray = 1;
}
void texture_t::load_png(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale) {

	//timer_t timer("Load PNG", allow_two_byte_grayscale);
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
	int const color_type(png_get_color_type(png_ptr, info_ptr));
	
	if (color_type == PNG_COLOR_TYPE_PALETTE) {
		png_set_expand(png_ptr); // expand to RGB
		png_read_update_info(png_ptr, info_ptr);
	}
	unsigned const w(png_get_image_width(png_ptr, info_ptr));
	unsigned const h(png_get_image_height(png_ptr, info_ptr));
	int const bit_depth(png_get_bit_depth(png_ptr, info_ptr));
	unsigned const png_ncolors(png_get_channels(png_ptr, info_ptr));
	set_image_size(w, h, allow_diff_width_height);
	bool const want_alpha_channel(ncolors == 4 && png_ncolors == 3);
	ncolors = png_ncolors;

	if (allow_two_byte_grayscale && ncolors == 1 && bit_depth == 16) {
		set_16_bit_grayscale();
		png_set_swap(png_ptr); // change big endian to little endian
	}
	else {
		if (bit_depth == 16) {png_set_strip_16(png_ptr);}
		if (bit_depth < 8)   {png_set_packing(png_ptr);}
	}
	vector<unsigned char *> rows(height);
	unsigned const scanline_size(ncolors*width);
	alloc();
	for (int i = 0; i < height; ++i) {rows[i] = data + (height - i - 1)*scanline_size;}
	png_read_image(png_ptr, &rows.front());
	png_read_end(png_ptr, end_info);
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	checked_fclose(fp);

	if (want_alpha_channel) {
		add_alpha_channel();
		auto_insert_alpha_channel(index);
	}
#else
	if (!load_stb_image(index, allow_diff_width_height, allow_two_byte_grayscale)) {
		cerr << "Error loading texture image file " << name << ": png support has not been enabled." << endl;
		exit(1);
	}
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
	int color_type(0), bit_depth;

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
	if (is_16_bit_gray) {png_set_swap(png_ptr);} // change big endian to little endian
	// Write image data
	for (int y = 0; y < height; y++) {png_write_row(png_ptr, (data + y*width*ncolors));}
	// End write
	png_write_end(png_ptr, NULL);
	png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	checked_fclose(fp);
	return 1;
#else
#ifdef ENABLE_STB_IMAGE
	if (!is_16_bit_gray) { // 16-bit grayscale is not yet supported
		stb_buffered_writer writer;
		if (!writer.init(fn)) return 0;
		if (!stbi_write_png_to_func(stb_buffered_writer_func, (void *)&writer, width, height, 3, data, 0)) return 0; // ncomp=3, stride_bytes=0 (packed)
		writer.flush(); // flush once at the end, which writes the entire image as one fwrite() call
		return 1;
	}
#endif
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
	uint32_t w(0), h(0);
	uint16_t bit_depth(0), config(0);
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,    &w);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH,   &h);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bit_depth);
	TIFFGetField(tif, TIFFTAG_PLANARCONFIG,  &config);
	set_image_size(w, h, allow_diff_width_height);
	
	if (allow_two_byte_grayscale && (ncolors == 0 || ncolors == 1) && bit_depth == 16) { // 16-bit grayscale
		set_16_bit_grayscale();
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
		uint32_t *raster = (uint32_t *)_TIFFmalloc(num_pixels() * sizeof(uint32_t));
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


// compressed/deferred load texture formats

void texture_t::deferred_load_and_bind() {
	assert(defer_load());

	switch (defer_load_type) {
	case DEFER_TYPE_DDS:  deferred_load_dds (); break;
	default:
		cerr << "Unhandled texture defer type " << defer_load_type << endl;
		exit(1);
	}
}

void texture_t::load_dds(int index) {
#ifdef ENABLE_DDS
	defer_load_type = DEFER_TYPE_DDS;
#else
	cerr << "Error loading texture image file " << name << ": DDS support has not been enabled." << endl;
	exit(1);
#endif
}

void texture_t::deferred_load_dds() {
#ifdef ENABLE_DDS
	//cout << "Loading DDS image " << name << endl;
	gli::texture2d Texture(gli::load_dds(name.c_str()));
	bool const compressed(gli::is_compressed(Texture.format()));
	// here we assume the texture is upside down and flip it, if it's uncompressed and flippable
	if (!compressed && !invert_y) {Texture = flip(Texture);}
	assert(!Texture.empty());
	width   = Texture.extent().x;
	height  = Texture.extent().y;
	ncolors = component_count(Texture.format());
	assert(width > 0 && height > 0);
	gli::gl GL(gli::gl::PROFILE_GL33);
	gli::gl::format const Format(GL.translate(Texture.format(), Texture.swizzles()));
	glBindTexture(GL_TEXTURE_2D, tid);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, GLint(Texture.levels() - 1));
	glTexStorage2D(GL_TEXTURE_2D, Texture.levels(), Format.Internal, width, height);

	// handle mipmaps; the value of use_mipmaps is ignored here
	for (gli::texture2d::size_type Level = 0; Level < Texture.levels(); ++Level) {
		if (compressed) {
			glCompressedTexSubImage2D(GL_TEXTURE_2D, Level, 0, 0, Texture[Level].extent().x, Texture[Level].extent().y,
				Format.Internal, Texture[Level].size(), Texture[Level].data());
		}
		else {
			glTexSubImage2D(GL_TEXTURE_2D, Level, 0, 0, Texture[Level].extent().x, Texture[Level].extent().y,
				Format.External, Format.Type, Texture[Level].data());
		}
	}
#else
	assert(0);
#endif
}


string read_string_ignore_comment_line(istream &in) {
	string s;
	in >> s;
	while (s.front() == '#') {
		in.ignore(numeric_limits<streamsize>::max(), '\n');
		in >> s;
	}
	return s;
}

bool open_texture_filebuf(filebuf &fb, string const &name) {
	return (fb.open(name, (ios::in | ios::binary)) || fb.open(append_texture_dir(name), (ios::in | ios::binary)));
}

// from Deliot2019
// https://eheitzresearch.wordpress.com/738-2/
void texture_t::load_ppm(int index, bool allow_diff_width_height) {

	filebuf fb;

	if (!open_texture_filebuf(fb, name)) {
		cerr << "Error: Couldn't open ppm file " << name << endl;
		exit(1);
	}
	istream in(&fb);
	ncolors = 3; // Only RGB format supported

	// Header parsing while ignoring comment lines
	if (read_string_ignore_comment_line(in) != string("P6")) { // Read magic number
		cerr << "Error: PPM file's magic number needs to be P6 in file " << name << endl;
		exit(1);
	}
	// Read width and height
	int const w(stoi(read_string_ignore_comment_line(in)));
	int const h(stoi(read_string_ignore_comment_line(in)));
	set_image_size(w, h, allow_diff_width_height);

	if (read_string_ignore_comment_line(in) != string("255")) { // Read max value
		cerr << "Error: PPM file's max value needs to be 255 in file " << name << endl;
		exit(1);
	}
	alloc();
	in.read(reinterpret_cast<char*>(data), num_bytes()); // Pixel data reading

	if (!in) {
		cerr << "Error reading PPM file" << endl;
		exit(1);
	}
}

#ifdef ENABLE_STB_IMAGE
// custom implementation for returning 8 or 16 bit data
stbi_uc *stbi_load_8_or_16(char const *filename, int *x, int *y, int *comp, bool *is_16bit, int req_comp=0) {
	FILE *f = stbi__fopen(filename, "rb");
	if (!f) {return stbi__errpuc("can't fopen", "Unable to open file");}
	stbi__context s;
	stbi__start_file(&s,f);
	stbi__result_info ri;
	void *result = stbi__load_main(&s, x, y, comp, req_comp, &ri, 16);
	if (result == NULL) return NULL;
	STBI_ASSERT(ri.bits_per_channel == 8 || ri.bits_per_channel == 16);
	*is_16bit = (ri.bits_per_channel == 16);

	if (stbi__vertically_flip_on_load) {
		int channels = req_comp ? req_comp : *comp;
		stbi__vertical_flip(result, *x, *y, channels * (ri.bits_per_channel >> 3));
	}
	fclose(f);
	return (stbi_uc *)result;
}
#endif

bool texture_t::load_stb_image(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale, unsigned char const *const load_from_data, unsigned load_from_size) {

#ifdef ENABLE_STB_IMAGE
	int w(0), h(0), nc(0);
	unsigned char *file_data(nullptr);
	stbi_set_flip_vertically_on_load(1);

	if (load_from_data != nullptr) { // load from memory
		file_data = stbi_load_from_memory(load_from_data, load_from_size, &w, &h, &nc, 0);
	}
	else { // load from file on disk
		string filename(append_texture_dir(name)); // first, try looking in the texture directory
		FILE *const fp(fopen(filename.c_str(), "rb")); // see if we can open the file
		if (fp == nullptr) {filename = name;} // if not in the texture directory, look in the current directory
		else {checked_fclose(fp);}
		
		if (allow_two_byte_grayscale) {
			bool is_16bit(0);
			file_data = stbi_load_8_or_16(filename.c_str(), &w, &h, &nc, &is_16bit); // Note: ~10% faster than libpng
			
			if (nc != 1) { // shouldn't be reachable, unless a RGB/RGBA image is specified as a heightmap
				cerr << "Error: Expecting grayscale PNG image, but got " << nc << " color component image: " << filename << endl;
				return 0;
			}
			if (is_16bit) {set_16_bit_grayscale();}
		}
		else {
			file_data = stbi_load(filename.c_str(), &w, &h, &nc, 0); // or stbi_load_from_file(fp, &w, &h, &nc, 0);
		}
	}
	if (file_data == nullptr) { // still not found, or there was an error
		cerr << "Error: stbi_load() returned error: \"" << stbi_failure_reason() << "\" for file " << name << endl;
		return 0; // default to some other image reader
	}
	//cout << TXT(name) << TXT(width) << TXT(height) << TXT(w) << TXT(h) << TXT(ncolors) << TXT(nc) << endl;
	set_image_size(w, h, allow_diff_width_height);
	bool const want_alpha_channel(ncolors == 4 && nc == 3);
	if (!is_16_bit_gray) {ncolors = nc;}
	alloc();
	std::memcpy(data, file_data, num_bytes()); // Note: stb uses malloc()/free(), but we use new[]/delete[], so need to copy the data
	stbi_image_free(file_data);

	if (want_alpha_channel) {
		add_alpha_channel();
		auto_insert_alpha_channel(index);
	}
	return 1;
#else
	cerr << "Error: stb_image support has not been enabled" << endl;
	return 0; // disabled
#endif
}

