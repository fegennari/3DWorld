// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/25/02
#include "targa.h"
#include "3DWorld.h"
#include "mesh.h"
#include "sinf.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"
#include "shaders.h"

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


float const TEXTURE_SMOOTH        = 0.01;
float const TEXTURE_SMOOTH_I      = 0.1;
float const SPHERE_SECTION        = 0.75;
float const TEXTURE_SNOW_ACC_RATE = 0.05;
int   const SNOW_ACC_RADIUS       = 3;
float const LANDSCAPE_REGEN_AMT   = 0.05; // how much regen after LANDSCAPE_REGEN_MOD ticks
int   const LANDSCAPE_REGEN_MOD   = 32;   // ticks for entire mesh regen pass
int   const LANDSCAPE_REGEN_RATE  = 400;  // ticks per landscape texture update (slow)
int   const MAX_UPDATE_SIZE       = 16;
float const ADJ_VALUE             = 1.0;

bool const RELOAD_TEX_ON_HOLE  = 0;
bool const LANDSCAPE_MIPMAP    = 0; // looks better, but texture update doesn't recompute the mipmaps
bool const SHOW_TEXTURE_MEMORY = 0;
bool const INSTANT_LTEX_COL    = 1;
bool const COMPRESS_TEXTURES   = 1;
bool const CHECK_FOR_LUM       = 1;
float const SMOOTH_SKY_POLES   = 0.2;

std::string const texture_dir("textures");


struct lspot {

	int x, y;
	float mag;
};


//GROUND ROCK WATER WATER2 SKY SUN MOON EARTH ICE SNOW LEAF WOOD SAND ROCK2 CAMOFLAGE GRASS PALM SMOKE PLASMA GEN LANDSCAPE TREE_END TREE_SNOW TREE_HEMI ...
//0      1    2     3      4   5   6    7     8   9    10   11   12   13    14        15    16   17    18     19  20        21       22        23        ...

texture_t textures[NUM_TEXTURES] = { // 4 colors without wrap sometimes has a bad transparent strip on spheres
// type: 0 = read from file, 1 = generated, 2 generated and dynamically updated
// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel), 4: targa (*tga), 5: jpeg, 6: png, 7: auto
// use_mipmaps: 0 = none, 1 = standard OpenGL, 2 = openGL + CPU data, 3 = custom alpha OpenGL
// type format width height wrap ncolors use_mipmaps name [do_compress] [anisotropy], [mipmap_alpha_weight]
//texture_t(0, 0, 512,  512,  1, 3, 0, "ground.raw"),
texture_t(0, 0, 128,  128,  1, 3, 2, "grass29.raw"), // mipmap for small trees?
texture_t(0, 0, 256,  256,  1, 3, 1, "rock.raw"),
texture_t(0, 0, 512,  512,  1, 3, 1, "water.raw"),
texture_t(0, 0, 64,   64,   1, 3, 1, "water_sm.raw"), // WATER2_TEX is unused
texture_t(0, 0, 1024, 1024, 1, 4, 0, "sky.raw"),
texture_t(0, 0, 64,   64,   1, 3, 1, "sun.raw"),
texture_t(0, 0, 128,  128,  1, 3, 1, "moon.raw"),
texture_t(0, 0, 256,  256,  0, 3, 1, "earth.raw"),
texture_t(0, 0, 64,   64,   1, 3, 1, "ice.raw"), // marble?
//texture_t(0, 0, 256,  256,  1, 3, 2, "snow.raw"),
texture_t(0, 7, 0,  0,  1, 3, 2, "snow2.jpg"),
texture_t(0, 0, 128,  128,  0, 4, 3, "leaf.raw"),
texture_t(0, 0, 128,  128,  1, 3, 0, "bark.raw"), // mipmap?
texture_t(0, 0, 512,  512,  1, 3, 2, "desert_sand.raw"),
texture_t(0, 0, 256,  256,  1, 3, 2, "rock2.raw"),
texture_t(0, 0, 512,  512,  1, 3, 1, "camoflage.raw"),
texture_t(0, 0, 128,  128,  1, 3, 0, "grass4.raw"),
texture_t(0, 1, 512,  512,  1, 3, 1, "brick1.bmp", 1, 8.0),
texture_t(0, 2, 512,  512,  1, 3, 1, "manhole.bmp"),
texture_t(0, 0, 128,  128,  1, 4, 3, "palmtree.raw"),
texture_t(1, 0, 256,  256,  1, 4, 1, "@smoke.raw"),  // not real file
texture_t(1, 0, 64,   64,   1, 4, 1, "@plasma.raw"), // not real file
texture_t(1, 0, 128,  128,  0, 3, 0, "@gen.raw"),    // not real file - unused
texture_t(2, 0, 1024, 1024, 0, 3, LANDSCAPE_MIPMAP, "final1024.raw"), // for loading real landscape texture
texture_t(1, 0, 128,  128,  0, 3, 0, "@tree_end.raw"),  // not real file
texture_t(1, 0, 128,  128,  1, 4, 1, "@tree_hemi.raw"), // not real file, mipmap for trees?
texture_t(1, 1, 512,  512,  1, 3, 1, "@shingle.bmp", 1, 8.0), // not real file
texture_t(0, 0, 256,  256,  1, 3, 1, "paneling.raw", 1, 16.0),
texture_t(0, 0, 256,  256,  1, 3, 1, "cblock.raw", 1, 8.0),
texture_t(0, 0, 128,  128,  0, 4, 3, "mj_leaf.raw"),
texture_t(0, 0, 128,  128,  0, 4, 3, "live_oak.raw"),
texture_t(0, 0, 256,  256,  0, 4, 3, "leaf2.raw"),
texture_t(0, 0, 256,  256,  0, 4, 3, "leaf3c.raw"),
texture_t(0, 0, 256,  256,  0, 4, 3, "plant1.raw"),
texture_t(0, 0, 256,  256,  0, 4, 3, "plant2.raw"),
texture_t(0, 0, 256,  256,  0, 4, 3, "plant3.raw"),
texture_t(0, 0, 64,   64,   0, 4, 3, "hibiscus.raw"),
texture_t(1, 0, 256,  256,  1, 3, 1, "@fence.raw", 1, 8.0), // not real file, light paneling
texture_t(0, 2, 128,  128,  1, 3, 1, "skull.raw"),
texture_t(0, 0, 64,   64,   1, 3, 1, "radiation.raw"),
texture_t(0, 2, 128,  128,  1, 3, 1, "yuck.raw"),
texture_t(0, 0, 256,  256,  0, 4, 0, "sawblade.raw"),
texture_t(0, 0, 256,  256,  0, 4, 0, "sawblade_b.raw"),
texture_t(0, 0, 256,  256,  0, 4, 1, "blur.raw"),
texture_t(0, 0, 256,  256,  1, 4, 1, "blur_s.raw"),
texture_t(0, 0, 256,  256,  0, 4, 3, "pine.raw", 1, 1.0, 0.36),
texture_t(0, 0, 128,  128,  1, 3, 1, "noise.raw"),
texture_t(0, 0, 128,  128,  1, 3, 1, "wood.raw"),
texture_t(0, 0, 128,  128,  1, 3, 1, "hb_brick.raw", 1, 8.0),
texture_t(0, 0, 128,  128,  1, 3, 1, "particleb.raw", 1, 8.0),
texture_t(0, 0, 128,  128,  1, 3, 1, "plaster.raw"),
texture_t(0, 0, 256,  256,  1, 3, 1, "tile.raw", 1, 8.0),
texture_t(0, 2, 256,  32,   1, 3, 1, "CommandCAD.raw"),
texture_t(1, 0, 32,   32,   1, 4, 1, "@disint.raw"),   // not real file
texture_t(1, 0, 256,  256,  1, 4, 1, "@blur_inv.raw"), // not real file
texture_t(1, 0, 32,   32,   1, 3, 0, "@hstripe.raw", 1, 8.0), // not real file
texture_t(1, 0, 32,   32,   1, 3, 0, "@vstripe.raw", 1, 8.0), // not real file
texture_t(0, 0, 512,  512,  1, 3, 1, "bcube.raw"),
texture_t(0, 0, 512,  512,  0, 4, 1, "explosion.raw"),
texture_t(0, 0, 512,  512,  1, 3, 1, "shiphull.raw"),
texture_t(0, 0, 512,  512,  1, 3, 1, "bcube2.raw"),
texture_t(0, 0, 512,  512,  1, 3, 1, "bcube_tactical.raw"),
texture_t(0, 0, 512,  256,  1, 3, 1, "rock_sphere.raw"),
texture_t(0, 3, 256,  256,  0, 4, 3, "papaya_leaf.raw"),
texture_t(0, 3, 256,  256,  0, 4, 3, "coffee_leaf.raw"), // half the texture is wasted, but leaves must be square (for now)
texture_t(0, 0, 256,  256,  1, 4, 0, "smiley_skull.raw"),
texture_t(0, 0, 512,  512,  1, 3, 1, "ice.2.raw"),
texture_t(0, 0, 256,  256,  1, 3, 2, "rock.03.raw"),
texture_t(0, 0, 16,   16,   1, 3, 0, "black.raw"),
texture_t(0, 0, 16,   16,   1, 3, 0, "white.raw"),
texture_t(0, 2, 512,  512,  0, 4, 0, "fire.raw"),
texture_t(0, 0, 1024, 1024, 1, 4, 1, "sky.raw"),
texture_t(0, 0, 256,  256,  0, 4, 0, "snowflake.raw"),
texture_t(1, 0, 128,  128,  0, 4, 1, "@blur_center.raw"), // not real file
texture_t(1, 0, 1,    128,  1, 4, 0, "@gradient.raw"), // not real file
texture_t(0, 0, 1024, 128,  0, 3, 1, "grass_blade.raw"),
texture_t(1, 0, 1024, 1024, 1, 1, 1, "@wind_texture.raw"),  // not real file
texture_t(0, 5, 0,    0,    1, 3, 1, "mossy_rock.jpg"), // 500x500
// bark
texture_t(0, 5, 0,    0,    1, 4, 1, "bark/bark1.jpg"), // 600x600
texture_t(0, 5, 0,    0,    1, 4, 1, "bark/bark2.jpg"), // 512x512
texture_t(0, 5, 0,    0,    1, 4, 1, "bark/bark2-normal.jpg"), // 512x512
texture_t(0, 5, 0,    0,    1, 4, 1, "bark/bark_lendrick.jpg"), // 892x892
texture_t(0, 6, 0,    0,    1, 4, 1, "bark/bark_lylejk.png"), // 1024x768

texture_t(1, 0, 128,  128,  1, 1, 0, "@noise_gen.raw") // not real file
//texture_t(0, 4, 0,    0,    1, 3, 1, "../Sponza2/textures/spnza_bricks_a_diff.tga")
// type format width height wrap ncolors use_mipmaps name [do_compress]
};


// zval should depend on def_water_level and temperature
float h_sand[NTEX_SAND], h_dirt[NTEX_DIRT], clip_hs1, clip_hs2, clip_hd1;
std::set<int> ls_color_texels;
vector<colorRGBA> cached_ls_colors;
typedef map<string, unsigned> name_map_t;
name_map_t texture_name_map;

int landscape_changed(0), lchanged0(0), skip_regrow(0), ltx1(0), lty1(0), ltx2(0), lty2(0), ls0_invalid(1);
unsigned char *landscape0 = NULL;


extern unsigned smoke_tid, dl_tid, elem_tid, gb_tid, flow_tid, reflection_tid;
extern int world_mode, island, read_landscape, default_ground_tex, xoff2, yoff2, DISABLE_WATER;
extern int scrolling, dx_scroll, dy_scroll, display_mode, iticks, universe_only;
extern float zmax, zmin, zmax_est, glaciate_exp, relh_adj_tex, vegetation, fticks;


void gen_smoke_texture();
void gen_plasma_texture();
void gen_disintegrate_texture();
void gen_tree_snow_texture();
void gen_tree_hemi_texture();
void gen_shingle_texture();
void gen_plant_texture();
void gen_fence_texture();
void gen_blur_inv_texture();
void gen_stripe_texture(int tid, bool horiz);
void gen_tree_end_texture();
void gen_blur_cent_texture();
void gen_gradient_texture();
void gen_wind_texture();
void gen_noise_texture();
void regrow_landscape_texture_amt0();
void update_lt_section(int x1, int y1, int x2, int y2);
int get_bare_ls_tid(float zval);

void free_universe_textures();


bool is_tex_disabled(int i) {

	return (universe_only && (i == CLOUD_RAW_TEX || i == WIND_TEX || i == LANDSCAPE_TEX || i == TREE_END_TEX || i == TREE_HEMI_TEX));
}


void load_texture_names() {

	if (!texture_name_map.empty()) return; // already loaded

	for (int i = 0; i < NUM_TEXTURES; ++i) {
		if (is_tex_disabled(i)) continue; // skip
		//assert(texture_name_map.find(textures[i].name) == texture_name_map.end());
		texture_name_map[textures[i].name] = i; // multiply used textures such as sky.raw will be overwritten
	}
}


void load_textures() {

	cout << "loading textures";
	load_texture_names();
	if (read_landscape) textures[LANDSCAPE_TEX].type = 0; // loaded from file

	for (int i = 0; i < NUM_TEXTURES; ++i) {
		cout.flush();
		cout << ".";
		if (!is_tex_disabled(i)) textures[i].load(i);
	}
	cout << endl;
	gen_smoke_texture();
	gen_plasma_texture();
	gen_disintegrate_texture();
	gen_shingle_texture();
	gen_fence_texture();
	gen_blur_inv_texture(); // must be after BLUR_TEX
	gen_stripe_texture(HSTRIPE_TEX, 1);
	gen_stripe_texture(VSTRIPE_TEX, 0);
	gen_blur_cent_texture();
	gen_gradient_texture();
	gen_noise_texture();

	if (!universe_only) {
		gen_wind_texture();
		gen_tree_hemi_texture();
		gen_tree_end_texture();
	}
	for (int i = 0; i < NUM_TEXTURES; ++i) {
		if (is_tex_disabled(i)) continue; // skip
		textures[i].init();
	}
	setup_multitexture();

	int max_tc(0), max_tu(0), max_tiu(0), max_ctiu(0);
	glGetIntegerv(GL_MAX_TEXTURE_COORDS,      &max_tc);
	glGetIntegerv(GL_MAX_TEXTURE_UNITS,       &max_tu);
	glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &max_tiu);
	glGetIntegerv(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS, &max_ctiu);
	cout << "max TCs: " << max_tc << ", max TUs: " << max_tu << ", max FS TIUs: " << max_tiu << ", max combined TIUs: " << max_ctiu << endl;
}


int get_texture_by_name(string const &name) {

	int const ix(atoi(name.c_str()));
	if (ix > 0 || ix == -1 || name == "0") return ix; // a number was specified
	name_map_t::const_iterator it(texture_name_map.find(name));

	if (it == texture_name_map.end()) {
		std::cerr << "Invalid texture name: " << name << endl;
		exit(0);
	}
	return it->second;
}


void check_init_texture(int id) {

	textures[id].check_init();
	//assert(glIsTexture(textures[id].tid)); // glIsTexture is slow on some machines???
}


bool select_texture(int id, bool enable, bool white_tex_default) {

	if (id < 0) {
		if (white_tex_default) {
			check_init_texture(WHITE_TEX);
			textures[WHITE_TEX].bind_gl();
		}
		else {
			glBindTexture(GL_TEXTURE_2D, 0); // bind to none
		}
		return 0;
	}
	assert(id < NUM_TEXTURES);
	check_init_texture(id);
	if (enable) glEnable(GL_TEXTURE_2D);
	textures[id].bind_gl();
	return 1;
}


float get_tex_ar(int id) {

	if (id < 0) return 0;
	assert(id < NUM_TEXTURES);
	return (((double)textures[id].width)/((double)textures[id].height));
}


void free_textures() {

	for (int i = 0; i < NUM_TEXTURES; ++i) {
		textures[i].gl_delete();
	}
}


void reset_textures() {

	cout << "Recreating textures..." << endl;
	free_textures();
	free_smiley_textures(); // should this be guarded by a conditional?
	free_universe_textures();
	free_flare_textures();
	free_shadow_map_textures();
	free_cloud_textures();
	free_texture(smoke_tid);
	free_texture(dl_tid);
	free_texture(elem_tid);
	free_texture(gb_tid);
	free_texture(flow_tid);
	free_texture(reflection_tid);
}


void setup_landscape_tex_colors(colorRGBA const &c1, colorRGBA const &c2) { // c1 = high, c2 = low

	textures[ROCK_TEX].set_to_color(c1);
	textures[SAND_TEX].set_to_color(c2);
	textures[DIRT_TEX].set_to_color(c2);
}


void free_texture(unsigned &tid) {

	if (tid == 0) return; // avoid GL calls
	if (glIsTexture(tid)) glDeleteTextures(1, &tid);
	tid = 0;
}


void texture_t::alloc() {

	free();
	data = new unsigned char[num_bytes()];

	if (SHOW_TEXTURE_MEMORY) {
		static unsigned tmem(0);
		tmem += num_bytes();
		cout << "tex mem = " << tmem << endl;
	}
}


void texture_t::bind_gl() const {
	
	assert(tid > 0);
	//assert(glIsTexture(tid));
	glBindTexture(GL_TEXTURE_2D, tid);
}


void texture_t::free_mm_data() {

	delete [] mm_data;
	mm_data = NULL;
	mm_offsets.clear();
}

void texture_t::free() {

	gl_delete(); // ???
	if (orig_data    != data) delete [] orig_data;
	if (colored_data != data) delete [] colored_data;
	delete [] data;
	data = orig_data = colored_data = NULL;
	free_mm_data();
}

void texture_t::gl_delete() {

	free_texture(tid);
}


void texture_t::init() {

	calc_color();
	build_mipmaps();
}


GLenum texture_t::calc_internal_format() const {

	static int has_comp(2); // starts unknown
	if (has_comp == 2) has_comp = has_extension("GL_ARB_texture_compression"); // unknown, calculate it
	assert(ncolors >= 1 && ncolors <= 4);
	return get_internal_texture_format(ncolors, (COMPRESS_TEXTURES && has_comp && do_compress && type != 2));
}

GLenum texture_t::calc_format() const {
	return get_texture_format(ncolors);
}


void texture_t::do_gl_init() {

	if (SHOW_TEXTURE_MEMORY) {
		static unsigned tmem(0);
		unsigned tsize(num_bytes());
		if (use_mipmaps) tsize = 4*tsize/3;
		tmem += tsize;
		cout << "tex vmem = " << tmem << endl;
	}
	//cout << "bind texture " << name << " size " << width << "x" << height << endl;
	//RESET_TIME;
	assert(is_allocated() && width > 0 && height > 0);
	setup_texture(tid, GL_MODULATE/*GL_DECAL*/, (use_mipmaps != 0), wrap, wrap, 0, 0, 0, anisotropy);
	//if (use_mipmaps == 1 || use_mipmaps == 2) glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
	glTexImage2D(GL_TEXTURE_2D, 0, calc_internal_format(), width, height, 0, calc_format(), GL_UNSIGNED_BYTE, data);
	if (use_mipmaps == 1 || use_mipmaps == 2) gen_mipmaps();
	if (use_mipmaps == 3) create_custom_mipmaps();
	assert(glIsTexture(tid));
	//PRINT_TIME("Texture Init");
}


void texture_t::calc_color() {

	assert(is_allocated());
	float colors[4] = {0.0}, weight(0.0);
	unsigned const size(num_pixels());
	bool const has_alpha_comp(has_alpha());
	has_binary_alpha = 1;

	for(unsigned i = 0; i < size; ++i) {
		int const offset(i*ncolors);
		float const cscale(has_alpha_comp ? data[offset+3]/255.0 : 1.0); // alpha scale
		weight += cscale;

		if (ncolors == 1) { // grayscale luminance
			UNROLL_3X(colors[i_] += cscale*data[offset];);
		}
		else {
			UNROLL_3X(colors[i_] += cscale*data[offset+i_];);
			
			if (has_alpha_comp) {
				colors[3] += data[offset+3];
				if (data[offset+3] > 0 && data[offset+3] < 255) has_binary_alpha = 0;
			}
		}
	}
	UNROLL_3X(color[i_] = colors[i_]/(255.0*weight);)
	color.alpha = (has_alpha_comp ? colors[3]/(255.0*weight) : 1.0);
}


void texture_t::copy_alpha_from_texture(texture_t const &at, bool alpha_in_red_comp) {

	assert(at.is_allocated() && is_allocated()); // check that data is allocated in both textures
	assert(ncolors == 4); // check for alpha channel
	assert(!is_bound());  // check that texture isn't already bound
	assert(at.ncolors == 1 || at.ncolors == 4);
	
	if (at.width != width || at.height != height) {
		resize(at.width, at.height);
		assert(at.width == width && at.height == height);
	}
	unsigned const npixels(num_pixels()), alpha_offset(alpha_in_red_comp ? 0 : 3);
	bool const is_lum(at.ncolors == 1);

	for (unsigned i = 0; i < npixels; ++i) {
		data[4*i+3] = (is_lum ? at.data[i] : at.data[4*i+alpha_offset]); // copy alpha values
	}
}


void texture_t::build_mipmaps() {

	if (use_mipmaps != 2) return; // not enabled
	assert(width == height);
	if (!mm_offsets.empty()) {assert(mm_data); return;} // already built
	assert(mm_data == NULL);
	unsigned data_size(0);

	for (unsigned tsz = width/2; tsz >= 1; tsz /= 2) {
		mm_offsets.push_back(data_size);
		data_size += ncolors*tsz*tsz;
	}
	mm_data = new unsigned char[data_size];
	GLenum const format(calc_format());

	for (unsigned level = 0; level < mm_offsets.size(); ++level) {
		unsigned const tsz(width >> level);
		assert(tsz > 1);
		int const ret(gluScaleImage(format, tsz,   tsz,   GL_UNSIGNED_BYTE, get_mipmap_data(level),
			                                tsz/2, tsz/2, GL_UNSIGNED_BYTE, (mm_data + mm_offsets[level])));
		if (ret) cout << "GLU error during mipmap image scale: " << gluErrorString(ret) << "." << endl;
	}
}


unsigned char const *texture_t::get_mipmap_data(unsigned level) const {

	if (level == 0) return get_data(); // base texture
	assert(level-1 < mm_offsets.size());
	assert(mm_data != NULL);
	return (mm_data + mm_offsets[level-1]);
}


void texture_t::set_to_color(colorRGBA const &c) {

	assert(is_allocated());
	if (c == color) return; // already set
	if (c == ALPHA0 && (orig_data == NULL || data == orig_data)) return; // color disabled (but never enabled)
	color = c;
	gl_delete();

	if (c == ALPHA0) { // color disabled
		data = orig_data;
		free_mm_data();
		build_mipmaps();
		return;
	}
	color_wrapper c4;
	c4.set_c4(c);
	unsigned const size(num_pixels());
	float const cw_scale(1.0/(float(c4.c[0]) + float(c4.c[1]) + float(c4.c[2])));
	if (colored_data == NULL) colored_data = new unsigned char[num_bytes()];
	if (orig_data    == NULL) orig_data    = data; // make a copy
	data = colored_data;
	assert(data != NULL);

	for (unsigned i = 0; i < size; ++i) {
		unsigned const pos(i*ncolors);
		unsigned char *d(data + pos);
		float const cscale(min(1.0f, (unsigned(d[0]) + unsigned(d[1]) + unsigned(d[2]))*cw_scale));
		
		for (int n = 0; n < ncolors; ++n) {
			d[n] = (unsigned char)min(255.0, (0.5*cscale*c4.c[n] + 0.5*orig_data[pos+n]));
		}
	}
	free_mm_data();
	build_mipmaps();
}


colorRGBA texture_color(int tid) {

	assert(tid >= 0 && tid < NUM_TEXTURES);
	return textures[tid].get_avg_color();
}


unsigned get_texture_size(int tid, bool dim) {

	assert(tid >= 0 && tid < NUM_TEXTURES);
	return (dim ? textures[tid].height : textures[tid].width);
}


void get_lum_alpha(colorRGBA const &color, int tid, float &luminance, float &alpha) {

	luminance = color.get_luminance();
	alpha     = color.alpha;

	if (tid >= 0) {
		colorRGBA const tc(texture_color(tid));
		luminance *= tc.get_luminance();
		alpha     *= tc.alpha;
	}
}


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


void texture_t::load(int index) {

	if (type > 0) { // generated texture
		alloc();
	}
	else {
		if (format == 7) { // auto
			// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel), 4: targa (*tga), 5: jpeg, 6: png, 7: auto
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
			else {
				cerr << "Error: Unidentified image file format for autodetect: " << ext << " in filename " << name << endl;
				exit(1);
			}
		}
		unsigned want_alpha_channel(ncolors == 4);

		switch (format) {
		case 0: load_raw_bmp(index); break; // raw
		case 1: load_raw_bmp(index); break; // bmp
		case 2: load_raw_bmp(index); break; // raw
		case 3: load_raw_bmp(index); break; // raw
		case 4: load_targa(); break;
		case 5: load_jpeg(); break;
		case 6: load_png(); break;
		}
		if (want_alpha_channel && ncolors < 4) {
			add_alpha_channel();
		}
		else {
			try_compact_to_lum();
		}
		fix_word_alignment();
	}
}


// load an .RAW or .BMP file as a texture
// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel)
void texture_t::load_raw_bmp(int index) {

	assert(ncolors == 1 || ncolors == 3 || ncolors == 4);
	if (format == 3) assert(ncolors == 4);
	FILE *file = open_texture_file(name); // open texture data
	assert(file != NULL);
	if (format == 1 && !verify_bmp_header(file, 0)) exit(1);

	// allocate buffer
	unsigned const size(num_pixels());
	assert(!is_allocated());
	alloc();
	float const ssp_inv_sq((SMOOTH_SKY_POLES > 0.0) ? 1.0/(SMOOTH_SKY_POLES*SMOOTH_SKY_POLES) : 0.0);
	int alpha_white(0);
	unsigned char buf[4], alpha;

	// read texture data
	if (ncolors == 4 && format != 3) { // add alpha
		for(unsigned i = 0; i < size; ++i) {
			int const i4(i << 2);
			size_t const nread(fread(buf, 3, 1, file)); assert(nread == 1);

			if (index == BLUR_TEX || index == SBLUR_TEX || index == BLUR_CENT_TEX) { // could use grayscale texture
				RGBA_BLOCK_ASSIGN((data+i4), 255, 255, 255, buf[0]); // alpha - assumes buf[0] = buf[1] = buf[2]
				continue;
			}
			if (i == 0) { // key off of first (llc) pixel
				alpha_white = (index == SMILEY_SKULL_TEX) ? 0 : ((int)buf[0] + (int)buf[1] + (int)buf[2] > 400);
			}
			if (format == 1) { // BGR => RGB
				RGB_BLOCK_ASSIGN((data+i4), buf[2], buf[1], buf[0]);
			}
			else {
				RGB_BLOCK_COPY((data+i4), buf);
			}
			if (index == CLOUD_TEX || index == CLOUD_RAW_TEX) {
				// white -> alpha = 255
				// blue  -> alpha = 0
				float const val(float(buf[0]) + float(buf[1]));
				//tex_data[i*4+3] = (unsigned char)(0.5*val);
				alpha = ((val <= 340.0) ? 0 : ((unsigned char)1.0*(val - 340.0)));

				if (SMOOTH_SKY_POLES > 0.0 && index == CLOUD_TEX) {
					unsigned const y(i/width);
					float const dist(float(y)/float(height)), d2(min(dist, (1.0f - dist)));
					if (d2 < SMOOTH_SKY_POLES) alpha = (unsigned char)(d2*d2*ssp_inv_sq*float(alpha));
				}
			}
			else { // make white/black part transparent, for example leaves
				float const val(float(buf[0]) + float(buf[1]) + float(buf[2]));
				
				if (index == EXPLOSION_TEX || index == FIRE_TEX) {
					alpha = (unsigned char)(0.333*val);
				}
				else if (alpha_white) {
					float const thresh((index == LEAF3_TEX) ? 700.0 : 600.0);
					alpha = ((val > thresh) ? 0 : ((val < thresh-100.0) ? 255 : (unsigned char)(2.55*(thresh - val))));
				}
				else {
					alpha = ((val < ((index == PINE_TEX) ? 65 : 32)) ? 0 : 255);
				}
			}
			data[i4+3] = alpha;
		}
	}
	else if (ncolors == 1) { // grayscale luminance
		vector<unsigned char> td2(size);
		size_t const nread(fread(&td2.front(), size, 1, file)); assert(nread == 1);

		for(unsigned i = 0; i < size; ++i) {
			RGB_BLOCK_COPY((data+3*i), td2);
		}
	}
	else {
		size_t const nread(fread(data, ncolors*size, 1, file)); assert(nread == 1);

		if (format == 1) {
			for(unsigned i = 0; i < size; ++i) {
				swap(data[3*i+0], data[3*i+2]); // BGR => RGB
			}
		}
	}
	if (format == 2) { // upside down
		unsigned const h2(height >> 1), wc(ncolors*width);
		
		for(unsigned i = 0; i < h2; ++i) {
			unsigned const off1(i*wc), off2((height-i-1)*wc);
			
			for(unsigned j = 0; j < wc; ++j) {
				swap(data[off1+j], data[off2+j]); // invert y
			}
		}
	}
	fclose(file);
}


void texture_t::load_targa() {

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
	if (width == 0 && height == 0) {
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


void texture_t::load_jpeg() {

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

	if (width == 0 && height == 0) {
		width  = cinfo.output_width;
		height = cinfo.output_height;
		assert(width > 0 && height > 0);
	}
	assert(cinfo.output_width == width && cinfo.output_height == height);
	ncolors = cinfo.output_components;
	unsigned const scanline_size(ncolors*width);
	alloc();

	while (cinfo.output_scanline < cinfo.output_height) {
		JSAMPROW row_pointer[1] = {data + scanline_size*(cinfo.output_height - cinfo.output_scanline - 1)};
		jpeg_read_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(fp);
#else
	cerr << "Error loading texture image file " << name << ": jpeg support has not been enabled." << endl;
	exit(1);
#endif
}


void texture_t::load_png() {

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

	if (width == 0 && height == 0) {
		width  = w;
		height = h;
		assert(width > 0 && height > 0);
	}
	assert(w == width && h == height);
	ncolors = png_get_channels(png_ptr, info_ptr);
	if (bit_depth == 16) png_set_strip_16(png_ptr);
	if (bit_depth < 8)   png_set_packing(png_ptr);
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
#else
	cerr << "Error loading texture image file " << name << ": png support has not been enabled." << endl;
	exit(1);
#endif
}


void texture_t::fix_word_alignment() {

	assert(is_allocated());
	unsigned const byte_align = 4;
	if ((ncolors*width & (byte_align-1)) == 0) return; // nothing to do
	float const ar(float(width)/float(height));
	int const new_w(width - (width&(byte_align-1)) + byte_align); // round up to next highest multiple of byte_align
	int const new_h(int(new_w/ar + 0.5)); // preserve aspect ratio
	resize(new_w, new_h);
}


void texture_t::add_alpha_channel() {

	assert(is_allocated());
	if (ncolors == 4) return; // alread has an alpha channel
	assert(ncolors == 1 || ncolors == 3); // only supported for RGB => RGBA for now
	bool const luminance(ncolors == 1);
	ncolors = 4; // add alpha channel
	unsigned char *new_data(new unsigned char[num_bytes()]);
	unsigned const npixels(num_pixels());

	for (unsigned i = 0; i < npixels; ++i) {
		UNROLL_3X(new_data[4*i+i_] = (luminance ? data[i] : data[3*i+i_]););
		new_data[4*i+3] = 255; // alpha of 255
	}
	free();
	data = new_data;
}


void texture_t::resize(int new_w, int new_h) {

	assert(is_allocated());
	if (new_w == width && new_h == height) return; // already correct size
	assert(data != NULL && width > 0 && height > 0 && new_w > 0 && new_h > 0);
	unsigned char *new_data(new unsigned char[new_w*new_h*ncolors]);
	int const ret(gluScaleImage(calc_format(), width, height, GL_UNSIGNED_BYTE, data, new_w, new_h, GL_UNSIGNED_BYTE, new_data));
	if (ret) cout << "GLU error during image scale: " << gluErrorString(ret) << "." << endl;
	free();
	data   = new_data;
	width  = new_w;
	height = new_h;
}


void texture_t::try_compact_to_lum() {

	assert(is_allocated());
	if (!CHECK_FOR_LUM || ncolors != 3) return;
	// determine if it's really a luminance texture
	unsigned const npixels(num_pixels());
	bool is_lum(1);

	for (unsigned i = 0; i < npixels && is_lum; ++i) {
		is_lum &= (data[3*i+1] == data[3*i] && data[3*i+2] == data[3*i]);
	}
	if (!is_lum) return;
	// RGB equal, make it a single color (luminance) channel
	ncolors = 1; // add alpha channel
	unsigned char *new_data(new unsigned char[num_bytes()]);

	for (unsigned i = 0; i < npixels; ++i) {
		new_data[i] = data[3*i];
	}
	free();
	data = new_data;
}


void texture_t::make_normal_map() {

	assert(is_allocated());
	if (ncolors == 3) return; // already a normal map
	
	if (ncolors == 4) { // better not have an alpha component
		cout << "Error: Skipping RGBA Bump/Normal map " << name << "." << endl;
		return;
	}
	assert(ncolors == 1); // grayscale heightmap
	ncolors = 3; // convert to RGB
	unsigned char *new_data(new unsigned char[num_bytes()]);
	int max_delta(1);

	for (int y = 0; y < height; ++y) { // assume texture wraps
		int const ym1((y-1+height) % height), yp1((y+1) % height);

		for (int x = 0; x < width; ++x) {
			int const xm1((x-1+width) % width), xp1((x+1) % width), off(3*(x + y*width));
			max_delta = max(max_delta, abs((int)data[xp1+y*width] - (int)data[xm1+y*width]));
			max_delta = max(max_delta, abs((int)data[x+yp1*width] - (int)data[x+ym1*width]));
		}
	}
	float max_delta_f(max_delta);

	for (int y = 0; y < height; ++y) { // assume texture wraps
		int const ym1((y-1+height) % height), yp1((y+1) % height);

		for (int x = 0; x < width; ++x) {
			int const xm1((x-1+width) % width), xp1((x+1) % width), off(3*(x + y*width));
			vector3d n(-((int)data[xp1+y*width] - (int)data[xm1+y*width])/max_delta_f,
				        ((int)data[x+yp1*width] - (int)data[x+ym1*width])/max_delta_f, 1.0);
			n.normalize();
			UNROLL_3X(new_data[off+i_] = unsigned char(127.5*(n[i_]+1.0)););
		}
	}
	free();
	data = new_data;
}


void texture_t::create_custom_mipmaps() {

	assert(is_allocated());
	GLenum const format(calc_format());
	unsigned const tsize(num_bytes());
	vector<unsigned char> idata, odata;
	idata.resize(tsize);
	memcpy(&idata.front(), data, tsize);

	for (unsigned w = width, h = height, level = 1; w > 1 || h > 1; w >>= 1, h >>= 1, ++level) {
		unsigned const w1(max(w,    1U)), h1(max(h,    1U));
		unsigned const w2(max(w>>1, 1U)), h2(max(h>>1, 1U));
		unsigned const xinc((w2 < w1) ? ncolors : 0), yinc((h2 < h1) ? ncolors*w1 : 0);
		odata.resize(ncolors*w2*h2);

		for (unsigned y = 0; y < h2; ++y) {
			for (unsigned x = 0; x < w2; ++x) {
				unsigned const ix1(ncolors*(y*w2+x)), ix2(ncolors*((y<<1)*w1+(x<<1)));

				if (ncolors == 1) {
					odata[ix1] = (unsigned char)(((unsigned)idata[ix2] + idata[ix2+xinc] + idata[ix2+yinc] + idata[ix2+yinc+xinc]) >> 2);
				}
				else if (ncolors == 3) {
					UNROLL_3X(odata[ix1+i_] = (unsigned char)(((unsigned)idata[ix2+i_] + idata[ix2+xinc+i_] + idata[ix2+yinc+i_] + idata[ix2+yinc+xinc+i_]) >> 2);)
				}
				else { // custom alpha mipmaps
					assert(ncolors == 4);
					unsigned const a1(idata[ix2+3]), a2(idata[ix2+xinc+3]), a3(idata[ix2+yinc+3]), a4(idata[ix2+yinc+xinc+3]);
					unsigned const a_sum(max(1U, (a1 + a2 + a3 + a4))); // no div by 0
					UNROLL_3X(odata[ix1+i_] = (unsigned char)((a1*idata[ix2+i_] + a2*idata[ix2+xinc+i_] + a3*idata[ix2+yinc+i_] + a4*idata[ix2+yinc+xinc+i_]) / a_sum);)
					odata[ix1+3] = min(255U, min(max(max(a1, a2), max(a3, a4)), unsigned(mipmap_alpha_weight*a_sum)));
				}
			}
		}
		glTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, format, GL_UNSIGNED_BYTE, &odata.front());
		idata.swap(odata);
	}
}


void bind_1d_texture(unsigned tid) {

	glBindTexture(GL_TEXTURE_1D, tid);
	assert(glIsTexture(tid));
}


void bind_2d_texture(unsigned tid) {

	glBindTexture(GL_TEXTURE_2D, tid);
	assert(glIsTexture(tid));
}


void setup_texture(unsigned &tid, int type, bool mipmap, bool wrap_s, bool wrap_t, bool mirror_s, bool mirror_t, bool nearest, float anisotropy) {

	assert(tid == 0);
	assert(!nearest || !mipmap);
	glGenTextures(1, &tid);

	// select our current texture
	bind_2d_texture(tid);

	// select modulate to mix texture with color for shading (decal keeps texture color)
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, type);

	// when texture area is small, use linear filter (bilinear filter the closest mipmap)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (nearest ? GL_NEAREST : (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR))); // GL_LINEAR_MIPMAP_NEAREST?

	// when texture area is large, bilinear filter the first mipmap
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (nearest ? GL_NEAREST : GL_LINEAR));

	// enable anisotropic filtering (slower but higher quality)
	if (anisotropy > 1.0) glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, anisotropy);

	// if wrap is true,  the texture wraps over at the edges (repeat) or is mirrored
	// if wrap is false, the texture ends at the edges (clamp)
	int const mode_s(wrap_s ? (mirror_s ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE); // Note: clamp is more efficient than wrap
	int const mode_t(wrap_t ? (mirror_t ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode_s);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode_t);
}


void texture_t::gen_rand_texture(unsigned char val, unsigned char a_add, unsigned a_rand) {

	assert(ncolors == 4);
	unsigned const size(num_pixels());

	for (unsigned i = 0; i < size; ++i) {
		RGBA_BLOCK_ASSIGN((data+(i<<2)), val, val, val, (a_add + (unsigned char)(rand() % a_rand)));
	}
}


void gen_smoke_texture() {

	unsigned char const smoke_color(255);
	textures[SMOKE_TEX].gen_rand_texture(smoke_color, 0, 256); // same as PLASMA_TEX but larger
}


void gen_plasma_texture() {

	textures[PLASMA_TEX].gen_rand_texture(255, 0, 256);
}


void gen_disintegrate_texture() {

	textures[DISINT_TEX].gen_rand_texture(255, 230, 26);
}


void gen_tree_hemi_texture() {

	assert(SPHERE_SECTION >= 0.0 && SPHERE_SECTION <= 1.0);
	unsigned char const *grass_tex_data(textures[GRASS_TEX].get_data());
	texture_t &tex(textures[TREE_HEMI_TEX]);
	unsigned char *tex_data(tex.get_data());
	int const sphere_h(int((1.0 - SPHERE_SECTION)*tex.height));

	for (int i = 0; i < sphere_h; ++i) { // green
		int const iw(i*tex.width);
		for (int j = 0; j < tex.width; ++j) {
			int const offset((iw + j) << 2);
			RGBA_BLOCK_ASSIGN((tex_data+offset), 0, 0, 0, 0);
		}
	}
	for (int i = sphere_h; i < tex.height; ++i) { // alpha = 0.0 - transparent
		int const iw(i*tex.width);
		for (int j = 0; j < tex.width; ++j) {
			int const offset((iw + j) << 2), offset2(3*(iw + j));
			RGB_BLOCK_COPY((tex_data+offset), (grass_tex_data+offset2));
			tex_data[offset+3] = 255; // A
		}
	}
}


void gen_shingle_texture() {

	texture_t &tex(textures[SHINGLE_TEX]);
	unsigned char const *tex_data2(textures[BRICK_TEX].get_data());
	assert(tex_data2 != NULL && tex.width == textures[BRICK_TEX].width && tex.height == textures[BRICK_TEX].height);
	assert(tex.ncolors == 3 && textures[BRICK_TEX].ncolors == 3);
	unsigned char *tex_data(tex.get_data());
	unsigned const size(tex.num_pixels());

	for (unsigned i = 0; i < size; ++i) { // convert brick texture to grayscale
		unsigned const offset(3*i);
		unsigned char const val((unsigned char)(((unsigned)tex_data2[offset+0] + (unsigned)tex_data2[offset+1] + (unsigned)tex_data2[offset+2])/3));
		RGB_BLOCK_ASSIGN((tex_data+offset), val, val, val);
	}
}


void gen_fence_texture() {

	texture_t &tex(textures[FENCE_TEX]);
	unsigned char const *tex_data2(textures[PANELING_TEX].get_data());
	assert(tex_data2 != NULL && tex.width == textures[PANELING_TEX].width && tex.height == textures[PANELING_TEX].height);
	assert(tex.ncolors == 3 && textures[PANELING_TEX].ncolors == 3);
	unsigned char *tex_data(tex.get_data());
	unsigned const size(tex.num_bytes());

	for (unsigned i = 0; i < size; ++i) { // convert to lighter color
		tex_data[i] = (unsigned char)min((unsigned)255, ((unsigned)tex_data2[i]) << 2);
	}
}


void gen_blur_inv_texture() {

	texture_t &tex(textures[BLUR_TEX_INV]);
	assert(tex.ncolors == 4);
	unsigned char *tex_data(tex.get_data());
	memset(tex_data, 255, 4*tex.num_pixels()*sizeof(char));

	for (int i = 0; i < tex.height; ++i) {
		unsigned char val(max(64, min(255, (2*255*(tex.height-i-1))/tex.height)));
		
		for (int j = 0; j < tex.width; ++j) {
			tex_data[((i*tex.width+j)<<2)+3] = val;
		}
	}
}


void gen_stripe_texture(int tid, bool horiz) {

	texture_t &tex(textures[tid]);
	assert(tex.ncolors == 3);
	unsigned char *tex_data(tex.get_data());

	for (int i = 0; i < tex.height; ++i) {
		for (int j = 0; j < tex.width; ++j) {
			unsigned char const val(255*(((horiz ? i : j)&3) != 0));
			UNROLL_3X(tex_data[3*(i*tex.width+j)+i_] = val;)
		}
	}
}


void gen_tree_end_texture() {

	texture_t &tex(textures[TREE_END_TEX]);
	unsigned char *tex_data(tex.get_data());
	int const w2(tex.width >> 1), h2(tex.height >> 1);
	float const scale_vals[3] = {160.0, 130.0, 80.0};

	for (int i = 0; i < tex.height; ++i) {
		for (int j = 0; j < tex.width; ++j) {
			float const rsin(sin(sqrt(float((j - w2)*(j - w2) + (i - h2)*(i - h2)))));
			float const darkness(0.9 + 0.1*rsin);
			int const offset(3*(i*tex.width + j));
			UNROLL_3X(tex_data[offset+i_] = (unsigned char)(scale_vals[i_]*darkness + 20.0*signed_rand_float());)
		}
	}
}


void gen_blur_cent_texture() {

	texture_t &tex(textures[BLUR_CENT_TEX]);
	assert(textures[BLUR_CENT_TEX].ncolors == 4);
	unsigned char *tex_data(tex.get_data());
	int const w2(tex.width >> 1), h2(tex.height >> 1);
	float const scale(2.0/min(tex.width, tex.height));

	for (int i = 0; i < tex.height; ++i) {
		for (int j = 0; j < tex.width; ++j) {
			float const radius(sqrt(float((j - w2)*(j - w2) + (i - h2)*(i - h2)))*scale);
			int const offset(4*(i*tex.width + j));
			UNROLL_3X(tex_data[offset+i_] = 255;)
			tex_data[offset+3] = (unsigned char)(255.0*(1.0 - CLIP_TO_01(radius))); // linear scaling for alpha
		}
	}
}


void gen_gradient_texture() { // for horizon

	texture_t &tex(textures[GRADIENT_TEX]); // 1D
	assert(tex.width == 1 && tex.ncolors == 4);
	unsigned char *tex_data(tex.get_data());
	int const size(tex.num_pixels()); // Note: int, not unsigned

	for (int i = 0; i < size; ++i) {
		UNROLL_3X(tex_data[4*i+i_] = 255;)
		tex_data[4*i+3] = (unsigned char)max(0, (255 * 2 * (size/2 - abs(i - (int)size/2)) / size)); // linear gradient
	}
}


void gen_wind_texture() {

	texture_t &tex(textures[WIND_TEX]);
	unsigned char const *tex_data2(textures[CLOUD_RAW_TEX].get_data());
	assert(tex.ncolors == 1 && textures[CLOUD_RAW_TEX].ncolors == 4); // RGBA => grayscale luminance
	assert(tex_data2 != NULL && tex.width == textures[CLOUD_RAW_TEX].width && tex.height == textures[CLOUD_RAW_TEX].height);
	unsigned char *tex_data(tex.get_data());
	unsigned const size(tex.num_pixels());

	for (unsigned i = 0; i < size; ++i) {
		tex_data[i] = tex_data2[(i<<2)+3]; // put alpha in luminance
	}
}


void noise_fill(unsigned char *data, unsigned size) {
	for (unsigned i = 0; i < size; ++i) {data[i] = (rand() % 256);}
}


void gen_noise_texture() {

	texture_t &tex(textures[NOISE_GEN_TEX]);
	assert(tex.ncolors == 1);
	noise_fill(tex.get_data(), tex.num_pixels());
}


unsigned create_3d_noise_texture(unsigned size) {

	vector<unsigned char> data(size*size*size);
	noise_fill(&data.front(), data.size());
	return create_3d_texture(size, size, size, 1, data, GL_LINEAR, GL_REPEAT);
}


colorRGBA get_landscape_texture_color(int xpos, int ypos) {

	if (cached_ls_colors.empty()) cached_ls_colors.resize(XY_MULT_SIZE, ALPHA0);
	unsigned const ix(xpos + MESH_X_SIZE*ypos);
	assert(ix < cached_ls_colors.size());
	if (cached_ls_colors[ix].alpha != 0.0) return cached_ls_colors[ix];
	assert(!point_outside_mesh(xpos, ypos));
	texture_t &tex(textures[LANDSCAPE_TEX]);
	assert(tex.ncolors == 3);
	unsigned char const *tex_data(tex.get_data());
	unsigned const xstep(tex.width/MESH_X_SIZE), ystep(tex.height/MESH_Y_SIZE);
	unsigned const x0(xpos*xstep), y0(ypos*ystep), x1(x0 + xstep), y1(y0 + ystep);
	colorRGBA color(BLACK);

	for (unsigned y = y0; y < y1; ++y) { // like creating a mipmap
		for (unsigned x = x0; x < x1; ++x) {
			UNROLL_3X(color[i_] += tex_data[3*(tex.width*y + x) + i_];)
		}
	}
	color *= 1.0/(255.0*xstep*ystep);
	cached_ls_colors[ix] = color;
	return color;
}


inline int get_bare_ls_tid(float zval) {

	float const relh(relh_adj_tex + (zval - zmin)/(zmax - zmin));

	if (island) {
		if (relh > clip_hs1) { // deep into grass
			return ((relh > clip_hs2) ? ROCK_TEX : DIRT_TEX); // rock or dirt
		}
		else {
			return SAND_TEX; // sand
		}
	}
	else {
		return ((relh > clip_hd1) ? ROCK_TEX : DIRT_TEX); // rock or dirt
	}
	return 0;
}


void update_lttex_ix(int &ix) { // note: assumes lttex_dirt (no islands)

	if (DISABLE_WATER == 2 && lttex_dirt[ix].id == SNOW_TEX  ) --ix;
	if (vegetation == 0.0  && lttex_dirt[ix].id == GROUND_TEX) ++ix;
}


void get_tids(float relh, int NTEXm1, float const *const h_tex, int &k1, int &k2, float &t) {

	for (k1 = 0; k1 < NTEXm1 && relh >= h_tex[k1]; ++k1) {} // find first texture with height greater than relh
	float const blend_border(island ? TEXTURE_SMOOTH_I : TEXTURE_SMOOTH);

	if (k1 < NTEXm1 && (h_tex[k1] - relh) < blend_border) {
		t  = 1.0 - (h_tex[k1] - relh)/blend_border;
		k2 = k1+1;
		update_lttex_ix(k1);
		update_lttex_ix(k2);
	}
	else {
		update_lttex_ix(k1);
		k2 = k1;
	}
}


// multiple resolution texture passes?
void create_landscape_texture() {

	RESET_TIME;
	if (read_landscape || world_mode == WMODE_INF_TERRAIN) return;
	cached_ls_colors.clear();
	int tox0(0), toy0(0), scroll(0);
	texture_t &tex(textures[LANDSCAPE_TEX]);
	int const width(tex.width), height(tex.height);
	unsigned char *tex_data(tex.get_data());
	assert(tex.ncolors == 3);
	int x1(0), y1(0), x2(width), y2(height);
	static int tox(0), toy(0);
	if (!scrolling) cout << "Generating landscape terrain texture." << endl;

	if (scrolling) { // ensure texture alignment when scrolling
		tox0   = (dx_scroll*width) /MESH_X_SIZE;
		toy0   = (dy_scroll*height)/MESH_Y_SIZE;
		tox   += tox0;
		toy   += toy0;
		x1     = max(0, -tox0);
		y1     = max(0, -toy0);
		x2     = min(width,  width -tox0);
		y2     = min(height, height-toy0);
		scroll = (x1 < x2 && y1 < y2 && abs(tox0) < width && abs(toy0) < height);
	}
	float const dz(zmax - zmin), dz_inv(1.0/dz);
	float const *h_tex(island ? h_sand     : h_dirt);
	ttex  const *lttex(island ? lttex_sand : lttex_dirt);
	int   const   NTEX(island ? NTEX_SAND  : NTEX_DIRT);
	int const mxszm1(MESH_X_SIZE-1), myszm1(MESH_Y_SIZE-1), dxv(width/MESH_X_SIZE), dyv(height/MESH_Y_SIZE);
	int const NTEXm1(NTEX-1), def_id((default_ground_tex >= 0) ? default_ground_tex : GROUND_TEX);
	int const id0(lttex[NTEXm1].id);
	float const xscale(((float)MESH_X_SIZE)/((float)width)), yscale(((float)MESH_Y_SIZE)/((float)height));
	static char **tids = NULL;
	if (tids == NULL) matrix_gen_2d(tids);
	assert(NTEX < 128);
	
	for (int i = 0; i < MESH_Y_SIZE; ++i) { // makes a big performance improvement
		int const keepy(scroll && (i+dy_scroll) > 0 && (i+dy_scroll) < MESH_Y_SIZE-1);

		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (keepy && (j+dx_scroll) > 0 && (j+dx_scroll) < MESH_X_SIZE-1) continue;
			int const i1(min(myszm1, i+1)), j1(min(mxszm1, j+1));
			float const mh00(mesh_height[i][j]), mh01(mesh_height[i][j1]), mh10(mesh_height[i1][j]), mh11(mesh_height[i1][j1]);
			float const dist(fabs(mh01 - mh00) + fabs(mh00 - mh10) + fabs(mh01 - mh11) + fabs(mh10 - mh11));
			float const relh1(relh_adj_tex + (min(min(mh00, mh01), min(mh10, mh11)) - zmin)*dz_inv);
			float const relh2(relh_adj_tex + (max(max(mh00, mh01), max(mh10, mh11)) - zmin)*dz_inv);
			int k1a, k1b, k2a, k2b;
			float t;
			get_tids(relh1, NTEXm1, h_tex, k1a, k2a, t);
			get_tids(relh2, NTEXm1, h_tex, k1b, k2b, t);
			tids[i][j] = char((k1a == k2b) ? k1a : -1);
		}
	}
	int const wx(width - dxv), hy(height - dyv);
	int const i0((toy0 < 0) ? height-1 : 0), i1((toy0 < 0) ? -1 : height), di((toy0 < 0) ? -1 : 1);
	int const j0((tox0 < 0) ? width -1 : 0), j1((tox0 < 0) ? -1 : width ), dj((tox0 < 0) ? -1 : 1);
	int const wxtx(wx-tox0), wxtx3(3*wxtx), j00(max(0, -tox0)), j01(min(width, wx-tox0));
	
	for (int i = i0; i != i1; i += di) {
		int j10(j0);
		float const yp(yscale*(float)i);
		int const lly(i + toy0), offset(3*i*width), ypos(max(0, min(myszm1, (int)yp))), ypos1(min(myszm1, ypos+1));
		float const ypi(yp - (float)ypos);

		if (scroll && lly >= 0 && lly < hy) {
			memmove((tex_data+offset+3*j00), (tex_data+3*((j00+tox0)+lly*width)), 3*(j01-j00)); // range could be overlapping
			j10 = ((tox0 < 0) ? min(j0, -tox0-1) : max(j0, wx-tox0));
		}
		for (int j = j10; j != j1; j += dj) {
			float const xp(xscale*(float)j);
			int const o2(offset + 3*j), xpos(max(0, min(mxszm1, (int)xp))), xpos1(min(mxszm1, xpos+1));
			float const xpi(xp - (float)xpos);
			float t(0.0);
			int k1, k2, id, id2;

			if (default_ground_tex >= 0 || zmax == zmin) {
				id = id2 = def_id;
				k1 = k2  = 0;
			}
			else {
				if (tids[ypos][xpos] < 0) {
					float const mh00(mesh_height[ypos][xpos]), mh01(mesh_height[ypos][xpos1]), mh10(mesh_height[ypos1][xpos]), mh11(mesh_height[ypos1][xpos1]);
					float const mh((1.0 - xpi)*((1.0 - ypi)*mh00 + ypi*mh10) + xpi*((1.0 - ypi)*mh01 + ypi*mh11));
					float const relh(relh_adj_tex + (mh - zmin)*dz_inv);
					get_tids(relh, NTEXm1, h_tex, k1, k2, t);
					if (k1 != k2) assert(k2 == k1+1 || vegetation == 0.0);
				}
				else {
					k1 = k2 = tids[ypos][xpos];
				}
				id  = lttex[k1].id;
				id2 = lttex[k2].id;
			}
			texture_t const &t1(textures[id]);
			unsigned char const *t1_data(t1.get_data());
			int const tof(t1.ncolors*(((i+toy)&(t1.height-1))*t1.width + ((j+tox)&(t1.width-1))));

			if (k1 == k2) { // single texture
				RGB_BLOCK_COPY((tex_data + o2), (t1_data + tof));
			}
			else { // blend two textures - performance critical
				texture_t const &t2(textures[id2]);
				unsigned char const *t2_data(t2.get_data());
				int const tof2(t2.ncolors*(((i+toy)&(t2.height-1))*t2.width + ((j+tox)&(t2.width-1))));
				BLEND_COLOR((tex_data + o2), (t2_data + tof2), (t1_data + tof), t);
			}

			// handle steep slopes (dirt/rock texture replaces grass texture)
			bool const grass(id == GROUND_TEX || id2 == GROUND_TEX), snow(id2 == SNOW_TEX);
			if (!grass && !snow) continue;
			float const *const sti(sthresh[island][snow]);
			float const vnz00(vertex_normals[ypos][xpos].z);
			if (vnz00 > sti[1]+0.1) continue; // not steep enough
			float const vnz01(vertex_normals[ypos][xpos1].z), vnz10(vertex_normals[ypos1][xpos].z), vnz11(vertex_normals[ypos1][xpos1].z);
			float const vnz((1.0 - xpi)*((1.0 - ypi)*vnz00 + ypi*vnz10) + xpi*((1.0 - ypi)*vnz01 + ypi*vnz11));

			if (grass && vnz < sti[1]) { // ground/grass
				texture_t const &ta(textures[DIRT_TEX]);
				unsigned char const *ta_data(ta.get_data());
				int const tofa(ta.ncolors*(((i+toy)&(ta.height-1))*ta.width + ((j+tox)&(ta.width-1))));
				unsigned char temp[3];

				if (id == GROUND_TEX || id2 == ROCK_TEX) {
					texture_t const &tb(textures[ROCK_TEX]);
					unsigned char const *tb_data(tb.get_data());
					int const tofb(tb.ncolors*(((i+toy)&(tb.height-1))*tb.width + ((j+tox)&(tb.width-1))));
					BLEND_COLOR(temp, (tb_data+tofb), (ta_data+tofa), t);
				}
				else {
					RGB_BLOCK_COPY(temp, (ta_data+tofa));
				}
				float const val(CLIP_TO_01((vnz - sti[0])/(sti[1] - sti[0])));
				BLEND_COLOR((tex_data+o2), (tex_data+o2), temp, val);
			}
			else if (vnz < sti[1]) { // snow
				texture_t const &ta(textures[ROCK_TEX]);
				unsigned char const *ta_data(ta.get_data());
				int const tofa(ta.ncolors*(((i+toy)&(ta.height-1))*ta.width + ((j+tox)&(ta.width-1))));
				float const val(CLIP_TO_01(2.0f*(vnz - sti[0])/(sti[1] - sti[0])));
				BLEND_COLOR((tex_data+o2), (tex_data+o2), (ta_data+tofa), val);
			}
		} // for j
	} // for i
	PRINT_TIME(" Data Gen");

	if (!read_landscape) {
		if (landscape0 == NULL || !scroll) { // initialize/copy entire texture
			int const totsize(tex.num_bytes());
			if (landscape0 == NULL) landscape0 = new unsigned char[totsize];
			memcpy(landscape0, tex_data, totsize*sizeof(unsigned char));
			ls0_invalid = 0;
			PRINT_TIME(" Landscape0 Gen");
		}
		else if (lchanged0) { // create landscape0
			for (int i = i0; i != i1; i += di) {
				int const lly(i + toy0), off(3*i*width);

				if (lly >= 0 && lly < hy) { // copy part of this row from old landscape0
					memmove(landscape0+off+3*j00, landscape0+3*((j00+tox0)+lly*width), 3*(j01-j00)); // range could be overlapping
					if (-tox0 > 0)      memcpy(landscape0+off, tex_data+off, -3*tox0); // copy from beginning
					if (width-wxtx > 0) memcpy(landscape0+off+wxtx3, tex_data+off+wxtx3, 3*(width-wxtx)); // copy to end
				}
				else { // new section - copy entire row from tex_data
					memcpy(landscape0+off, tex_data+off, 3*width);
				}
			}
			ls0_invalid = 0;
			PRINT_TIME(" Landscape0 Copy");
		}
		else { // delay creation of landscape0 until it's needed
			ls0_invalid = 1;
		}
	}
	if (scrolling) { // supposedly more efficient
		tex.bind_gl();
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, tex.calc_format(), GL_UNSIGNED_BYTE, tex_data);
	}
	else {
		tex.gl_delete(); // should we try to update rather than recreating from scratch?
		tex.do_gl_init();
	}
	PRINT_TIME(" Final");
}


void regrow_landscape_texture_amt0() {

	//RESET_TIME;
	static int counter(0);
	if (read_landscape || iticks == 0) return; // is it too strong to use both iticks and fticks?
	if (ls0_invalid) create_landscape_texture();
	texture_t &tex(textures[LANDSCAPE_TEX]);
	assert(tex.is_allocated() && landscape0 != NULL);
	int const regen_bshift(max(1, min(7, int(log(1.0/max(0.0001f, fticks*LANDSCAPE_REGEN_AMT))/log(2.0)))));
	unsigned char *tex_data(tex.get_data());
	int const y1((counter%LANDSCAPE_REGEN_MOD)*(tex.height/LANDSCAPE_REGEN_MOD));
	int const y2(y1 + (tex.height/LANDSCAPE_REGEN_MOD));
	assert(y2 <= tex.height);

	if (LANDSCAPE_REGEN_AMT > 0.0 || skip_regrow) {
		for (int i = y1; i < y2; ++i) {
			int const i_step(i*tex.width);

			for (int j = 0; j < tex.width; ++j) { // performance critical
				int const offset(3*(i_step + j));
				UNROLL_3X(tex_data[offset+i_] += (unsigned char)(((int)landscape0[offset+i_] - (int)tex_data[offset+i_]) >> regen_bshift);)
			}
		}
	}
	++counter;
	skip_regrow = 0;
	check_init_texture(LANDSCAPE_TEX);
	tex.bind_gl();
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, y1, tex.width, (y2-y1), GL_RGB, GL_UNSIGNED_BYTE, (tex_data + 3*tex.width*y1));
	//PRINT_TIME("Regrow");
}


float add_crater_to_landscape_texture(float xval, float yval, float radius) {

	int const xpos(get_xpos(xval)), ypos(get_ypos(yval));
	if (point_outside_mesh(xpos, ypos)) return 0.0; // off the terrain area
	texture_t const &t1(textures[get_bare_ls_tid(mesh_height[ypos][xpos])]);
	texture_t &tex(textures[LANDSCAPE_TEX]);
	unsigned char const *data(t1.get_data());
	unsigned char *tex_data(tex.get_data());
	float const xscale(((float)MESH_X_SIZE)/((float)tex.width)), yscale(((float)MESH_Y_SIZE)/((float)tex.height));
	int const xpos2(int((xval + X_SCENE_SIZE)/(xscale*DX_VAL) + 0.5));
	int const ypos2(int((yval + Y_SCENE_SIZE)/(yscale*DY_VAL) + 0.5));
	if (xpos2 < 0 || ypos2 < 0 || xpos2 >= tex.width || ypos2 >= tex.height) return 0.0;
	int const rad(max(0, int(radius/((xscale + yscale)*(DX_VAL + DY_VAL))) - 1)), radsq(rad*rad); // size in texture space
	int const x1(max(0, xpos2-rad)), y1(max(0, ypos2-rad));
	int const x2(min(tex.width-1,  xpos2+rad)), y2(min(tex.height-1, ypos2+rad));

	for (int i = y1; i <= y2; ++i) {
		int const offset(i*tex.width), yterm((i - ypos2)*(i - ypos2));

		for (int j = x1; j <= x2; ++j) {
			if ((yterm + (j - xpos2)*(j - xpos2)) < radsq) {
				int const o2(3*(offset + j)), tof(t1.ncolors*((i&(t1.height-1))*t1.width + (j&(t1.width-1))));
				RGB_BLOCK_COPY((tex_data+o2), (data+tof));
			}
		}
	}
	lchanged0 = 1;
	cached_ls_colors.clear();
	update_lt_section(x1, y1, x2+1, y2+1);
	return get_rel_height(mesh_height[ypos][xpos], zmin, zmax);
}


void add_hole_in_landscape_texture(int xpos, int ypos, float blend) { // for water damage

	if (blend <= 0.0 || point_outside_mesh(xpos, ypos)) return; // off the terrain area
	blend = max(0.01f, min(1.0f, blend));
	texture_t const &t1(textures[get_bare_ls_tid(mesh_height[ypos][xpos])]);
	texture_t &tex(textures[LANDSCAPE_TEX]);
	unsigned char const *data(t1.get_data());
	unsigned char *tex_data(tex.get_data());
	float const xdiv(((((float)MESH_X_SIZE)/((float)tex.width))*DX_VAL)), ydiv(((((float)MESH_Y_SIZE)/((float)tex.height))*DY_VAL));
	int const xpos1(max(0,          int((get_xval(xpos)   + X_SCENE_SIZE)/xdiv + 0.5)));
	int const ypos1(max(0,          int((get_yval(ypos)   + Y_SCENE_SIZE)/ydiv + 0.5)));
	int const xpos2(min(tex.width,  int((get_xval(xpos+1) + X_SCENE_SIZE)/xdiv + 0.5)));
	int const ypos2(min(tex.height, int((get_yval(ypos+1) + Y_SCENE_SIZE)/ydiv + 0.5)));

	for (int i = ypos1; i < ypos2; ++i) {
		unsigned const iw(i*tex.width), ival((i&(t1.height-1))*t1.width);

		for (int j = xpos1; j < xpos2; ++j) {
			int const o2(3*(iw + j)), tof(t1.ncolors*(ival + (j&(t1.width-1))));
			BLEND_COLOR((tex_data + o2), (data + tof), (tex_data + o2), blend);
		}
	}
	lchanged0 = 1;

	if (RELOAD_TEX_ON_HOLE) {
		ltx1 = min(ltx1, xpos1);
		lty1 = min(lty1, ypos1);
		ltx2 = max(ltx2, xpos2);
		lty2 = max(lty2, ypos2);
		landscape_changed = 1;
	}
}


void add_color_to_landscape_texture(colorRGBA const &color, float xval, float yval, float radius, int check_unique) {

	int const xpos0(get_xpos(xval)), ypos0(get_ypos(yval));
	if (point_outside_mesh(xpos0, ypos0)) return; // off the terrain area
	if (wminside[ypos0][xpos0] && water_matrix[ypos0][xpos0] > mesh_height[ypos0][xpos0]) return; // underwater
	int const index(xpos0 + MESH_X_SIZE*ypos0);
	if (ls_color_texels.find(index) != ls_color_texels.end()) return;
	ls_color_texels.insert(index);
	texture_t &tex(textures[LANDSCAPE_TEX]);
	unsigned char *tex_data(tex.get_data());
	float const xscale(((float)MESH_X_SIZE)/((float)tex.width)), yscale(((float)MESH_Y_SIZE)/((float)tex.height));
	int const maxsize(tex.num_bytes()), tsize(tex.num_pixels()); // used only for assertions
	int const xpos(int((xval + X_SCENE_SIZE)/(((float)xscale)*DX_VAL) + 0.5));
	int const ypos(int((yval + Y_SCENE_SIZE)/(((float)yscale)*DY_VAL) + 0.5));
	if (xpos < 0 || ypos < 0 || xpos >= tex.width || ypos >= tex.height) return;
	int const rad(max(0, int(5.0*radius/((xscale + yscale)*(DX_VAL + DY_VAL))) - 1)), radsq(max(1, rad*rad)); // size in texture space
	int const x1(max(0, xpos-rad)), y1(max(0, ypos-rad)), x2(min(tex.width-1, xpos+rad)), y2(min(tex.height-1, ypos+rad));
	unsigned char color_i[3];
	unpack_color(color_i, color);

	for (int i = y1; i <= y2; ++i) {
		int const offset(i*tex.width), yterm((i - ypos)*(i - ypos));

		for (int j = x1; j <= x2; ++j) {
			int const dist_sq((yterm + (j - xpos)*(j - xpos)));

			if (dist_sq < radsq) {
				assert(offset + j <= tsize);
				float const blend((dist_sq == 0.0) ? 0.8 : min(0.8, 10.0/dist_sq));
				int const o2(3*(offset + j));
				assert(o2 + 3 <= maxsize);
				BLEND_COLOR((tex_data + o2), color_i, (tex_data + o2), blend); // test if this texel is enabled for draw?
			}
		}
	}
	if (INSTANT_LTEX_COL) {
		update_lt_section(x1, y1, x2+1, y2+1); // slow on certain hardware
	}
	else {
		lchanged0 = 1;
	}
}


void add_snow_to_landscape_texture(point const &pos, float acc) {

	if (acc <= 0.0) return;
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos)) return; // off the terrain area
	if (wminside[ypos][xpos] && water_matrix[ypos][xpos] > mesh_height[ypos][xpos]) return; // underwater
	texture_t &tex(textures[LANDSCAPE_TEX]);
	int const tx(int(((float)tex.width /(float)MESH_X_SIZE)*(pos.x + X_SCENE_SIZE)*DX_VAL_INV + 0.5));
	int const ty(int(((float)tex.height/(float)MESH_Y_SIZE)*(pos.y + Y_SCENE_SIZE)*DY_VAL_INV + 0.5));
	if (tx < 0 || tx >= tex.width || ty < 0 || ty >= tex.height) return;
	skip_regrow = 1;
	unsigned char *tex_data(tex.get_data());
	acc *= TEXTURE_SNOW_ACC_RATE;
	float const acc255(255.0*acc);
	int const rsq(SNOW_ACC_RADIUS*SNOW_ACC_RADIUS);
	float const frad((float)rsq), inv_frad(1.0/frad);
	int const x1(max(0, tx-SNOW_ACC_RADIUS)), x2(min(tex.width -1, tx+SNOW_ACC_RADIUS));
	int const y1(max(0, ty-SNOW_ACC_RADIUS)), y2(min(tex.height-1, ty+SNOW_ACC_RADIUS));

	for (int i = y1; i <= y2; ++i) {
		int const dy(i*tex.width), ysq((ty - i)*(ty - i));
		if (ysq > rsq) continue; // can never satisfy condition below

		for (int j = x1; j <= x2; ++j) {
			int const dist(ysq + (tx - j)*(tx - j));
			if (dist > rsq) continue;
			int const offset(3*(dy + j));

			if (acc >= 1.0) { // white
				RGB_BLOCK_ASSIGN((tex_data+offset), 255, 255, 255);
				continue;
			}
			// somewhat white
			float const mult((frad - (float)dist)*inv_frad), acc255s(acc255*mult), omacc(1.0 - acc*mult);
			UNROLL_3X(tex_data[offset+i_] = (unsigned char)min(255, int(tex_data[offset+i_]*omacc + acc255s));)
		}
	}
	lchanged0 = 1;
}


void update_landscape_texture() {

	//RESET_TIME;
	ls_color_texels.clear();
	if (landscape_changed)  lchanged0 = 1; // this is set if landscape ever changed
	if (lchanged0)          regrow_landscape_texture_amt0();
	if (!landscape_changed) return;
	
	if (ltx2 <= ltx1 || lty2 <= lty1) {
		landscape_changed = 0;
		return;
	}
	int y2(lty2);
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);

	if ((ltx2 - ltx1)*(lty2 - lty1) > width*height/MAX_UPDATE_SIZE) { // only update a small section per frame
		int const usize((width*height/MAX_UPDATE_SIZE)/(ltx2 - ltx1));
		y2   = min(height, (lty1 + usize)); // select a smaller section
		lty1 = y2; // advance lty1 to beginning of unprocessed section
	}
	else { // small section
		ltx1 = width;
		lty1 = height;
		ltx2 = lty2 = 0;
		landscape_changed = 0;
	}
	update_lt_section(ltx1, lty1, ltx2, y2);
	//PRINT_TIME("Update");
}


inline int get_align(int x1, int x2, unsigned am) {

	return ((1 + am - ((x2 - x1) & am)) & am);
}


void update_lt_section(int x1, int y1, int x2, int y2) {

	//RESET_TIME;
	if (x1 == x2 || y1 == y2) return;
	texture_t const &t1(textures[LANDSCAPE_TEX]);
	unsigned const nc(t1.ncolors);
	int const width(t1.width), height(t1.height);
	assert(nc == 3 && width > 0 && height > 0 && t1.is_allocated());

	if (!(x1 < x2 && y1 < y2 && x1 >= 0 && y1 >= 0 && x2 <= width && y2 <= height)) {
		cout << "x1 = " << x1 << ", y1 = " << y1 << ", x2 = " << x2 << ", y2 = " << y2 << endl;
		assert(0);
	}
	unsigned const align(4), am(align - 1);

	if (am > 0) { // force 4-byte alignment
		int xadd(get_align(x1, x2, am));
		x2   = min((x2 + xadd), width);
		xadd = get_align(x1, x2, am);
		x1   = max((x1 - xadd), 0);
		assert(((x2 - x1) & am) == 0);
	}
	int const offset(nc*(x1 + y1*width));
	unsigned char const *data(t1.get_data() + offset);
	unsigned char *copy(NULL);

	if ((x1 > 0 || x2 < width) && (y2 - y1) > 1) { // copy data if more than one row
		if ((x2 - x1) > width/2) { // just copy the entire strip
			x1 = 0;
			x2 = width;
		}
		else {
			unsigned const cw(nc*(x2 - x1));
			copy = new unsigned char[cw*(y2 - y1)];

			for (int i = 0; i < (y2 - y1); ++i) {
				memcpy((copy + cw*i), (data + nc*i*width), cw*sizeof(unsigned char));
			}
			data = copy;
		}
	}
	check_init_texture(LANDSCAPE_TEX);
	t1.bind_gl();
	glTexSubImage2D(GL_TEXTURE_2D, 0, x1, y1, (x2-x1), (y2-y1), t1.calc_format(), GL_UNSIGNED_BYTE, data);
	delete [] copy;
	//PRINT_TIME("LT Update");
}


int snow_height(point pos) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos)) return 0;
	double const relh(relh_adj_tex + (mesh_height[ypos][xpos] - zmin)/(zmax - zmin));
	return (island ? (relh > h_sand[2]) : (relh > h_dirt[2]));
}


void gen_tex_height_tables() {

	for (unsigned i = 0; i < NTEX_SAND; ++i) {
		h_sand[i] = pow(lttex_sand[i].zval, glaciate_exp);
	}
	for (unsigned i = 0; i < NTEX_DIRT; ++i) {
		h_dirt[i] = pow(lttex_dirt[i].zval, glaciate_exp);
	}
	clip_hs1 = (0.25*h_sand[1] + 0.75*h_sand[0]);
	clip_hs2 = (0.90*h_sand[1] + 0.10*h_sand[0]);
	clip_hd1 = (0.90*h_dirt[1] + 0.10*h_dirt[0]);
}


void set_texgen_vec4(float const v[4], bool s_or_t, bool enable_and_set_mode, shader_t *shader) {

	if (shader) {
		shader->add_attrib_float_array((s_or_t ? TEX0_T_ATTR : TEX0_S_ATTR), v, 4);
	}
	else {
		if (enable_and_set_mode) {
			glEnable(s_or_t ? GL_TEXTURE_GEN_T : GL_TEXTURE_GEN_S);
			glTexGeni((s_or_t ? GL_T : GL_S), GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
		}
		glTexGenfv((s_or_t ? GL_T : GL_S), GL_EYE_PLANE, v);
	}
}


void setup_texgen_full(float sx, float sy, float sz, float sw, float tx, float ty, float tz, float tw, int mode, shader_t *shader) {

	float const tex_param_s[4] = {sx, sy, sz, sw};
	float const tex_param_t[4] = {tx, ty, tz, tw};
	set_texgen_vec4(tex_param_s, 0, 1, shader);
	set_texgen_vec4(tex_param_t, 1, 1, shader);
}


void setup_texgen(float xscale, float yscale, float tx, float ty, float z_off, int mode) {

	assert(xscale != 0.0 && yscale != 0.0);
	setup_texgen_full(xscale, 0.0, z_off, tx, 0.0, yscale, z_off, ty, mode);
}


void disable_texgen() {

	glDisable(GL_TEXTURE_GEN_S);
	glDisable(GL_TEXTURE_GEN_T);
}


void disable_textures_texgen() {

	disable_multitex_a();
	disable_texgen();
	glDisable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0); // bind to none
}


void setup_polygon_texgen(vector3d const &norm, float const scale[2], float const xlate[2],
	vector3d const &offset, bool swap_txy, shader_t *shader)
{
	int const d0(get_min_dim(norm));
	vector3d v[2] = {all_zeros, all_zeros};
	v[0][d0] = 1.0;
	cross_product(norm, v[0], v[1]);
	cross_product(norm, v[1], v[0]);

	for (unsigned i = 0; i < 2; ++i) {
		float const tex_param[4] = {scale[i]*v[i].x, scale[i]*v[i].y, scale[i]*v[i].z, (xlate[i] + scale[i]*dot_product(offset, v[i]))};
		set_texgen_vec4(tex_param, ((i != 0) ^ swap_txy), 1, shader);
	}
}


void get_tex_coord(vector3d const &dir, vector3d const &sdir, unsigned txsize, unsigned tysize, int &tx, int &ty, bool invert) {

	double s, t;
	dir_to_sphere_s_t(dir, sdir, s, t);
	tx = int(txsize*s)%txsize; // size is not always a power of 2
	ty = int(tysize*t)%tysize;

	if (invert) {
		tx = int(txsize) - tx - 1; // backwards
		ty = int(tysize) - ty - 1;
	}
	assert(tx >= 0 && ty >= 0 && tx < int(txsize) && ty < int(tysize));
}


float texture_t::get_component(float xval, float yval, int comp) const {

	assert(comp < ncolors);
	int tx(int(width *xval) & (width -1)); // width and height are a power of 2
	int ty(int(height*yval) & (height-1));
	if (tx < 0) tx += width;
	if (ty < 0) ty += height;
	assert(tx >= 0 && ty >= 0 && tx < width && ty < height);
	return data[ncolors*(width*ty + tx) + comp]/255.0;
}


float get_texture_component(unsigned tid, float xval, float yval, int comp) {

	assert(tid < NUM_TEXTURES);
	return textures[tid].get_component(xval, yval, comp);
}


bool is_billboard_texture_transparent(point const *const points, point const &pos, int tid) {

	// ordered: (s,t) => (0,1), (0,0), (1,0), (1,1)
	float d[4]; // distance from coll point to quad edge

	for (unsigned i = 0; i < 4; ++i) {
		unsigned const in((i+1)&3);
		d[i] = cross_product((pos - points[i]), (pos - points[in])).mag()/p2p_dist(points[i], points[in]);
	}
	assert(d[0] + d[2] > 0.0);
	assert(d[1] + d[3] > 0.0);
	float const tx(d[0]/((d[0] + d[2]))), ty(d[1]/(d[1] + d[3])); // y is upside down
	assert(tx >= 0.0 && tx <= 1.0 && ty >= 0.0 && ty <= 1.0);
	return (get_texture_component(tid, tx, ty, 3) == 0.0);
}


void set_texture_specular(bool val) {

	static bool inited(0), has_ext(0);

	if (!inited) {
		has_ext = has_extension("GL_EXT_separate_specular_color");
		inited  = 1;
	}
	if (!has_ext) return; // check for extension enabled
	glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL_EXT, (val ? GL_SEPARATE_SPECULAR_COLOR_EXT : GL_SINGLE_COLOR_EXT));
}




