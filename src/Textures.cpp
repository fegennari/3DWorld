// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/25/02
#include "GL/glew.h" // must be included first
#include "3DWorld.h"
#include "mesh.h"
#include "sinf.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"

using std::string;


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
float const SMOOTH_SKY_POLES   = 0.2;

std::string const texture_dir("textures");


struct lspot {

	int x, y;
	float mag;
};


//GROUND ROCK WATER WATER2 SKY SUN MOON EARTH ICE SNOW LEAF WOOD SAND ROCK2 CAMOFLAGE GRASS PALM SMOKE PLASMA GEN LANDSCAPE TREE_END TREE_SNOW TREE_HEMI ...
//0      1    2     3      4   5   6    7     8   9    10   11   12   13    14        15    16   17    18     19  20        21       22        23        ...

texture textures[NUM_TEXTURES] = { // 4 colors without wrap sometimes has a bad transparent strip on spheres
// type: 0 = read from file, 1 = generated, 2 generated and dynamically updated
// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel)
// use_mipmaps: 0 = none, 1 = standard OpenGL, 2 = openGL + CPU data, 3 = custom alpha OpenGL
// type format width height wrap ncolors use_mipmaps ([data] name [id] [color])
//texture(0, 0, 512,  512,  1, 3, 0, "ground.raw"),
texture(0, 0, 128,  128,  1, 3, 2, "grass29.raw"), // mipmap for small trees?
texture(0, 0, 256,  256,  1, 3, 1, "rock.raw"),
texture(0, 0, 512,  512,  1, 3, 1, "water.raw"),
texture(0, 0, 64,   64,   1, 3, 1, "water_sm.raw"), // WATER2_TEX is unused
texture(0, 0, 1024, 1024, 1, 4, 1, "sky.raw"),
texture(0, 0, 64,   64,   1, 3, 1, "sun.raw"),
texture(0, 0, 128,  128,  1, 3, 1, "moon.raw"),
texture(0, 0, 256,  256,  0, 3, 1, "earth.raw"),
texture(0, 0, 64,   64,   1, 3, 1, "ice.raw"), // marble?
texture(0, 0, 256,  256,  1, 3, 2, "snow.raw"),
texture(0, 0, 128,  128,  0, 4, 3, "leaf.raw"),
texture(0, 0, 128,  128,  1, 3, 0, "bark.raw"),
texture(0, 0, 512,  512,  1, 3, 2, "desert_sand.raw"),
texture(0, 0, 256,  256,  1, 3, 2, "rock2.raw"),
texture(0, 0, 512,  512,  1, 3, 1, "camoflage.raw"),
texture(0, 0, 128,  128,  1, 3, 0, "grass4.raw"),
texture(0, 1, 512,  512,  1, 3, 1, "brick1.bmp"),
texture(0, 2, 512,  512,  1, 3, 1, "manhole.bmp"),
texture(0, 0, 128,  128,  1, 4, 3, "palmtree.raw"),
texture(1, 0, 256,  256,  1, 4, 1, "@smoke.raw"),  // not real file
texture(1, 0, 64,   64,   1, 4, 1, "@plasma.raw"), // not real file
texture(1, 0, 128,  128,  0, 3, 0, "@gen.raw"),    // not real file - unused
texture(2, 0, 1024, 1024, 0, 3, LANDSCAPE_MIPMAP, "final1024.raw"), // for loading real landscape texture
texture(1, 0, 128,  128,  0, 3, 0, "@tree_end.raw"),  // not real file
texture(1, 0, 128,  128,  1, 4, 1, "@tree_hemi.raw"), // not real file, mipmap for trees?
texture(1, 1, 512,  512,  1, 3, 1, "@shingle.bmp"),   // not real file
texture(0, 0, 256,  256,  1, 3, 1, "paneling.raw"),
texture(0, 0, 256,  256,  1, 3, 1, "cblock.raw"),
texture(0, 0, 128,  128,  0, 4, 3, "mj_leaf.raw"),
texture(0, 0, 128,  128,  0, 4, 3, "live_oak.raw"),
texture(0, 0, 256,  256,  0, 4, 3, "leaf2.raw"),
texture(0, 0, 256,  256,  0, 4, 3, "leaf3c.raw"),
texture(0, 0, 256,  256,  0, 4, 3, "plant1.raw"),
texture(0, 0, 256,  256,  0, 4, 3, "plant2.raw"),
texture(0, 0, 256,  256,  0, 4, 3, "plant3.raw"),
texture(0, 0, 64,   64,   0, 4, 3, "hibiscus.raw"),
texture(0, 0, 256,  256,  1, 3, 1, "@fence.raw"), // not real file, light paneling
texture(0, 2, 128,  128,  1, 3, 1, "skull.raw"),
texture(0, 0, 64,   64,   1, 3, 1, "radiation.raw"),
texture(0, 2, 128,  128,  1, 3, 1, "yuck.raw"),
texture(0, 0, 256,  256,  0, 4, 0, "sawblade.raw"),
texture(0, 0, 256,  256,  0, 4, 0, "sawblade_b.raw"),
texture(0, 0, 256,  256,  0, 4, 1, "blur.raw"),
texture(0, 0, 256,  256,  1, 4, 1, "blur_s.raw"),
texture(0, 0, 256,  256,  0, 4, 0, "pine.raw"),
texture(0, 0, 128,  128,  1, 3, 1, "noise.raw"),
texture(0, 0, 128,  128,  1, 3, 1, "wood.raw"),
texture(0, 0, 128,  128,  1, 3, 1, "hb_brick.raw"),
texture(0, 0, 128,  128,  1, 3, 1, "particleb.raw"),
texture(0, 0, 128,  128,  1, 3, 1, "plaster.raw"),
texture(0, 0, 256,  256,  1, 3, 1, "tile.raw"),
texture(0, 2, 256,  32,   1, 3, 1, "CommandCAD.raw"),
texture(1, 0, 32,   32,   1, 4, 1, "@disint.raw"),   // not real file
texture(1, 0, 256,  256,  1, 4, 1, "@blur_inv.raw"), // not real file
texture(1, 0, 32,   32,   1, 3, 0, "@hstripe.raw"),  // not real file
texture(1, 0, 32,   32,   1, 3, 0, "@vstripe.raw"),  // not real file
texture(0, 0, 512,  512,  1, 3, 1, "bcube.raw"),
texture(0, 0, 512,  512,  0, 4, 1, "explosion.raw"),
texture(0, 0, 512,  512,  1, 3, 1, "shiphull.raw"),
texture(0, 0, 512,  512,  1, 3, 1, "bcube2.raw"),
texture(0, 0, 512,  512,  1, 3, 1, "bcube_tactical.raw"),
texture(0, 0, 512,  256,  1, 3, 1, "rock_sphere.raw"),
texture(0, 3, 256,  256,  0, 4, 3, "papaya_leaf.raw"),
texture(0, 3, 256,  256,  0, 4, 3, "coffee_leaf.raw"), // half the texture is wasted, but leaves must be square (for now)
texture(0, 0, 256,  256,  1, 4, 0, "smiley_skull.raw"),
texture(0, 0, 512,  512,  1, 3, 1, "ice.2.raw"),
texture(0, 0, 256,  256,  1, 3, 2, "rock.03.raw"),
texture(0, 0, 16,   16,   1, 3, 0, "black.raw"),
texture(0, 0, 16,   16,   1, 3, 0, "white.raw"),
texture(0, 2, 512,  512,  0, 4, 0, "fire.raw"),
texture(0, 0, 1024, 1024, 1, 4, 1, "sky.raw"),
texture(0, 0, 256,  256,  0, 4, 0, "snowflake.raw"),
texture(1, 0, 128,  128,  0, 4, 1, "@blur_center.raw"), // not real file
texture(1, 0, 1,    128,  1, 4, 0, "@gradient.raw"), // not real file
texture(0, 0, 1024, 128,  0, 3, 1, "grass_blade.raw"),
texture(1, 0, 1024, 1024, 1, 1, 1, "@wind_texture.raw"),  // not real file
// type format width height wrap ncolors use_mipmaps ([data] name [id] [color])
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
extern int scrolling, dx_scroll, dy_scroll, display_mode;
extern float zmax, zmin, zmax_est, glaciate_exp, relh_adj_tex, vegetation;


unsigned char *LoadTextureRAW(texture const &t, int index);
void gen_smoke_texture();
void gen_plasma_texture();
void gen_disintegrate_texture();
void gen_tree_snow_texture();
void gen_tree_hemi_texture();
void gen_shingle_texture();
void gen_plant_texture();
void gen_fence_texture();
void gen_blur_inv_texture();
void gen_hstripe_texture();
void gen_vstripe_texture();
void gen_tree_end_texture();
void gen_blur_cent_texture();
void gen_gradient_texture();
void gen_wind_texture();
void regrow_landscape_texture_amt0();
void update_lt_section(int x1, int y1, int x2, int y2);
int get_bare_ls_tid(float zval);
void alloc_texture(int id);

void free_universe_textures();


void load_textures() {

	static int init(0);
	cout << "loading textures";
	glEnable(GL_TEXTURE_2D);

	for (int i = 0; i < NUM_TEXTURES; ++i) {
		cout.flush();
		cout << ".";
		
		switch (i) {
		case SMOKE_TEX:     gen_smoke_texture();        break;
		case PLASMA_TEX:    gen_plasma_texture();       break;
		case DISINT_TEX:    gen_disintegrate_texture(); break;
		case TREE_HEMI_TEX: gen_tree_hemi_texture();    break;
		case SHINGLE_TEX:   gen_shingle_texture();      break;
		case FENCE_TEX:     gen_fence_texture();        break;
		case BLUR_TEX_INV:  gen_blur_inv_texture();     break; // must be after BLUR_TEX
		case HSTRIPE_TEX:   gen_hstripe_texture();      break;
		case VSTRIPE_TEX:   gen_vstripe_texture();      break;
		case BLUR_CENT_TEX: gen_blur_cent_texture();    break;
		case GRADIENT_TEX:  gen_gradient_texture();     break;
		case WIND_TEX:      gen_wind_texture();         break;

		default:
			if (textures[i].type > 0) { // generated texture
				alloc_texture(i);
			}
			else {
				textures[i].data = LoadTextureRAW(textures[i], i);
			}
		} // switch
		textures[i].init();
		//assert(texture_name_map.find(textures[i].name) == texture_name_map.end());
		texture_name_map[textures[i].name] = i; // multiply used textures such as sky.raw will be overwritten
	}
	if (read_landscape) {
		int const i(LANDSCAPE_TEX);
		textures[i].data = LoadTextureRAW(textures[i], i);
	}
	cout << endl;
	gen_tree_end_texture();
	glDisable(GL_TEXTURE_2D);
	setup_multitexture();
	if (!init) init = 1;
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

	// glIsTexture is slow on some machines???
	if (textures[id].tid == 0) textures[id].do_gl_init();
	//assert(glIsTexture(textures[id].tid));
}


bool select_texture(int id, bool enable) {

	if (id < 0) return 0;
	assert(id < NUM_TEXTURES);
	check_init_texture(id);
	if (enable) glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, textures[id].tid);
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

	if (glIsTexture(tid)) glDeleteTextures(1, &tid);
	tid = 0;
}


void texture::alloc() {

	free();
	data = new unsigned char[width*height*ncolors];
}

void texture::free_mm_data() {

	delete [] mm_data;
	mm_data = NULL;
	mm_offsets.clear();
}

void texture::free() {

	if (orig_data    != data) delete [] orig_data;
	if (colored_data != data) delete [] colored_data;
	delete [] data;
	data = orig_data = colored_data = NULL;
	free_mm_data();
}

void texture::gl_delete() {

	free_texture(tid);
}


void texture::init() {

	calc_color();
	build_mipmaps();
}


GLenum texture::calc_internal_format() const {

	static int has_comp(2); // starts unknown
	if (has_comp == 2) has_comp = has_extension("GL_ARB_texture_compression"); // unknown, calculate it

	if (COMPRESS_TEXTURES && has_comp && type != 2) {
		switch (ncolors) {
		case 1: return GL_COMPRESSED_LUMINANCE;
		case 3: return GL_COMPRESSED_RGB;
		case 4: return GL_COMPRESSED_RGBA;
		default: assert(0);
		}
	}
	return ncolors;
}


GLenum texture::calc_format() const {
	
	switch (ncolors) {
	case 1: return GL_LUMINANCE;
	case 3: return GL_RGB;
	case 4: return GL_RGBA; // GL_BGRA is supposedly faster, but do we want to swap things here?
	default: assert(0);
	}
	return 0;
}


void texture::do_gl_init() {

	if (SHOW_TEXTURE_MEMORY) {
		static unsigned tmem(0);
		unsigned tsize(width*height*ncolors);
		if (use_mipmaps) tsize = 4*tsize/3;
		tmem += tsize;
		cout << "tmem = " << tmem << endl;
	}
	assert(width > 0 && height > 0 && data != NULL);
	setup_texture(tid, GL_MODULATE/*GL_DECAL*/, (use_mipmaps != 0), wrap, wrap);
	//if (use_mipmaps) glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
	glTexImage2D(GL_TEXTURE_2D, 0, calc_internal_format(), width, height, 0, calc_format(), GL_UNSIGNED_BYTE, data);
	if (use_mipmaps == 1 || use_mipmaps == 2) gen_mipmaps();
	if (use_mipmaps == 3) create_custom_mipmaps();
	assert(glIsTexture(tid));
}


void texture::calc_color() {

	float colors[4] = {0.0,0.0,0.0,0.0}, weight(0.0);
	int const size(width*height);
	if (data == NULL) {cout << "NULL texture: " << name << endl; assert(0);}
	bool const has_alpha(ncolors == 4);

	for(int i = 0; i < size; ++i) {
		int const offset(i*ncolors);
		float const cscale(has_alpha ? data[offset+3]/255.0 : 1.0); // alpha scale
		weight += cscale;

		for (int j = 0; j < ncolors; ++j) {
			colors[j] += ((j < 3) ? cscale : 1.0)*data[offset+j];
		}
	}
	for (int j = 0; j < ncolors; ++j) {
		color[j] = colors[j]/(255.0*weight);
	}
	if (!has_alpha) color.alpha = 1.0;
}


void texture::build_mipmaps() {

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
		gluScaleImage(format, tsz,   tsz,   GL_UNSIGNED_BYTE, get_mipmap_data(level),
			                  tsz/2, tsz/2, GL_UNSIGNED_BYTE, (mm_data + mm_offsets[level]));
	}
}


unsigned char const *texture::get_mipmap_data(unsigned level) const {

	if (level == 0) { // base texture
		assert(data != NULL);
		return data;
	}
	assert(level-1 < mm_offsets.size());
	assert(mm_data != NULL);
	return (mm_data + mm_offsets[level-1]);
}


void texture::set_to_color(colorRGBA const &c) {

	assert(data != NULL);
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
	unsigned const sz(unsigned(width*height));
	float const cw_scale(1.0/(float(c4.c[0]) + float(c4.c[1]) + float(c4.c[2])));
	if (colored_data == NULL) colored_data = new unsigned char[sz*ncolors];
	if (orig_data == NULL) orig_data = data; // make a copy
	data = colored_data;
	assert(data != NULL);

	for (unsigned i = 0; i < sz; ++i) {
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
	return textures[tid].color;
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


FILE *open_texture_file(string filename) {

	filename = texture_dir + "/" + filename;
	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp != NULL) return fp;

	// if not in the current directory, then look in the parent directory
	fp = fopen((string("../") + filename).c_str(), "rb");

	if (fp == NULL) {
		cout << "Error loading image " << filename << endl;
		exit(1);
	}
	return fp;
}



// load an .RAW or .BMP file as a texture
// format: 0 = RAW, 1 = BMP, 2 = RAW (upside down), 3 = RAW (alpha channel)
unsigned char *LoadTextureRAW(texture const &t, int index) {

	assert(t.ncolors == 1 || t.ncolors == 3 || t.ncolors == 4);
	if (t.format == 3) assert(t.ncolors == 4);
	int alpha_white(0);
	unsigned char buf[4], alpha;
	unsigned const size(t.width*t.height);
	FILE *file = open_texture_file(t.name); // open texture data
	assert(file != NULL);
	if (t.format == 1 && !verify_bmp_header(file, 0)) exit(1);

	// allocate buffer
	unsigned char *tex_data(new unsigned char[size*t.ncolors]);
	float const ssp_inv_sq((SMOOTH_SKY_POLES > 0.0) ? 1.0/(SMOOTH_SKY_POLES*SMOOTH_SKY_POLES) : 0.0);

	// read texture data
	if (t.ncolors == 4 && t.format != 3) { // add alpha
		for(unsigned i = 0; i < size; ++i) {
			int const i4(i << 2);
			size_t const nread(fread(buf, 3, 1, file)); assert(nread == 1);

			if (index == BLUR_TEX || index == SBLUR_TEX || index == BLUR_CENT_TEX) { // could use grayscale texture
				RGBA_BLOCK_ASSIGN((tex_data+i4), 255, 255, 255, buf[0]); // alpha - assumes buf[0] = buf[1] = buf[2]
				continue;
			}
			if (i == 0) { // key off of first (llc) pixel
				alpha_white = (index == SMILEY_SKULL_TEX) ? 0 : ((int)buf[0] + (int)buf[1] + (int)buf[2] > 400);
			}
			if (t.format == 1) { // BGR => RGB
				RGB_BLOCK_ASSIGN((tex_data+i4), buf[2], buf[1], buf[0]);
			}
			else {
				RGB_BLOCK_COPY((tex_data+i4), buf);
			}
			if (index == CLOUD_TEX || index == CLOUD_RAW_TEX) {
				// white -> alpha = 255
				// blue  -> alpha = 0
				float const val(float(buf[0]) + float(buf[1]));
				//tex_data[i*4+3] = (unsigned char)(0.5*val);
				alpha = ((val <= 340.0) ? 0 : ((unsigned char)1.0*(val - 340.0)));

				if (SMOOTH_SKY_POLES > 0.0 && index == CLOUD_TEX) {
					unsigned const y(i/t.width);
					float const dist(float(y)/float(t.height)), d2(min(dist, (1.0f - dist)));
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
			tex_data[i4+3] = alpha;
		}
	}
	else if (t.ncolors == 1) { // grayscale luminance
		vector<unsigned char> td2(size);
		size_t const nread(fread(&td2.front(), size, 1, file)); assert(nread == 1);

		for(unsigned i = 0; i < size; ++i) {
			RGB_BLOCK_COPY((tex_data+3*i), td2);
		}
	}
	else {
		size_t const nread(fread(tex_data, t.ncolors*size, 1, file)); assert(nread == 1);

		if (t.format == 1) {
			for(unsigned i = 0; i < size; ++i) {
				swap(tex_data[3*i+0], tex_data[3*i+2]); // BGR => RGB
			}
		}
	}
	if (t.format == 2) { // upside down
		unsigned const h2(t.height >> 1), wc(t.ncolors*t.width);
		
		for(unsigned i = 0; i < h2; ++i) {
			unsigned const off1(i*wc), off2((t.height-i-1)*wc);
			
			for(unsigned j = 0; j < wc; ++j) {
				swap(tex_data[off1+j], tex_data[off2+j]); // invert y
			}
		}
	}
	fclose(file);
	return tex_data;
}


void texture::create_custom_mipmaps() {

	GLenum const format(calc_format());
	unsigned const tsize(ncolors*width*height);
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
					odata[ix1+3] = max(max(a1, a2), max(a3, a4));
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


void setup_texture(unsigned &tid, int type, bool mipmap, bool wrap_s, bool wrap_t, bool mirror_s, bool mirror_t) {

	assert(tid == 0);
	glGenTextures(1, &tid);

	// select our current texture
	bind_2d_texture(tid);

	// select modulate to mix texture with color for shading (decal keeps texture color)
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, type);

	// when texture area is small, use linear filter (bilinear filter the closest mipmap)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR)); // GL_LINEAR_MIPMAP_NEAREST?

	// when texture area is large, bilinear filter the first mipmap
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	// if wrap is true,  the texture wraps over at the edges (repeat) or is mirrored
	// if wrap is false, the texture ends at the edges (clamp)
	int const mode_s(wrap_s ? (mirror_s ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE);
	int const mode_t(wrap_t ? (mirror_t ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode_s);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode_t);
}


void init_texture(int id) {

	assert(id < NUM_TEXTURES);
	textures[id].do_gl_init();
}


void gen_rand_texture(int tid, unsigned char val, unsigned char a_add=0, unsigned a_rand=256) {

	assert(textures[tid].ncolors == 4);
	int const width(textures[tid].width);
	int const height(textures[tid].height);
	unsigned char *tex_data(new unsigned char[width*height*4]);

	for (int i = 0; i < height; ++i) {
		int const iw(i*width);
		for (int j = 0; j < width; ++j) {
			int const offset((iw + j) << 2);
			RGBA_BLOCK_ASSIGN((tex_data+offset), val, val, val, (a_add + (unsigned char)(rand() % a_rand)));
		}
	}
	textures[tid].data = tex_data;
}


void gen_smoke_texture() {

	unsigned char const smoke_color(255);
	gen_rand_texture(SMOKE_TEX, smoke_color, 0, 256); // same as PLASMA_TEX but larger
}


void gen_plasma_texture() {

	gen_rand_texture(PLASMA_TEX, 255, 0, 256);
}


void gen_disintegrate_texture() {

	gen_rand_texture(DISINT_TEX, 255, 230, 26);
}


void gen_tree_hemi_texture() {

	unsigned char *grass_tex_data(textures[GRASS_TEX].data);
	assert(SPHERE_SECTION >= 0.0 && SPHERE_SECTION <= 1.0);
	assert(grass_tex_data != NULL);
	int const width(textures[TREE_HEMI_TEX].width), height(textures[TREE_HEMI_TEX].height);
	int const sphere_h(int((1.0 - SPHERE_SECTION)*height));
	unsigned char *tex_data(new unsigned char[width*height*4]);

	for (int i = 0; i < sphere_h; ++i) { // green
		int const iw(i*width);
		for (int j = 0; j < width; ++j) {
			int const offset((iw + j) << 2);
			RGBA_BLOCK_ASSIGN((tex_data+offset), 0, 0, 0, 0);
		}
	}
	for (int i = sphere_h; i < height; ++i) { // alpha = 0.0 - transparent
		int const iw(i*width);
		for (int j = 0; j < width; ++j) {
			int const offset((iw + j) << 2), offset2(3*(iw + j));
			RGB_BLOCK_COPY((tex_data+offset), (grass_tex_data+offset2));
			tex_data[offset+3] = 255; // A
		}
	}
	textures[TREE_HEMI_TEX].data = tex_data;
}


void gen_shingle_texture() {

	int const width(textures[SHINGLE_TEX].width), height(textures[SHINGLE_TEX].height);
	unsigned char *tex_data2(textures[BRICK_TEX].data);
	assert(tex_data2 != NULL && width == textures[BRICK_TEX].width && height == textures[BRICK_TEX].height);
	assert(textures[SHINGLE_TEX].ncolors == 3 && textures[BRICK_TEX].ncolors == 3);
	unsigned char *tex_data(new unsigned char[3*width*height]);

	for (int i = 0; i < height; ++i) { // convert brick texture to grayscale
		int const iw(i*width);
		for (int j = 0; j < width; ++j) {
			int const offset(3*(iw + j));
			unsigned char const val((unsigned char)(((unsigned)tex_data2[offset+0] + (unsigned)tex_data2[offset+1] + (unsigned)tex_data2[offset+2])/3));
			RGB_BLOCK_ASSIGN((tex_data+offset), val, val, val);
		}
	}
	textures[SHINGLE_TEX].data = tex_data;
}


void gen_fence_texture() {

	int const width(textures[FENCE_TEX].width), height(textures[FENCE_TEX].height), size(3*width*height);
	unsigned char *tex_data2(textures[PANELING_TEX].data);
	assert(tex_data2 != NULL && width == textures[PANELING_TEX].width && height == textures[PANELING_TEX].height);
	assert(textures[FENCE_TEX].ncolors == 3 && textures[PANELING_TEX].ncolors == 3);
	unsigned char *tex_data(new unsigned char[size]);

	for (int i = 0; i < size; ++i) { // convert to lighter color
		tex_data[i] = (unsigned char)min((unsigned)255, ((unsigned)tex_data2[i]) << 2);
	}
	textures[FENCE_TEX].data = tex_data;
}


void gen_blur_inv_texture() {

	int const width(textures[BLUR_TEX_INV].width), height(textures[BLUR_TEX_INV].height), size(4*width*height);
	assert(textures[BLUR_TEX_INV].ncolors == 4);
	unsigned char *tex_data(new unsigned char[size]);
	memset(tex_data, 255, size*sizeof(char));

	for (int i = 0; i < height; ++i) {
		unsigned char val(max(64, min(255, (2*255*(height-i-1))/height)));
		
		for (int j = 0; j < width; ++j) {
			tex_data[((i*width+j)<<2)+3] = val;
		}
	}
	textures[BLUR_TEX_INV].data = tex_data;
}


void gen_stripe_texture(int tid, bool horiz) {

	int const width(textures[tid].width), height(textures[tid].height), size(3*width*height);
	assert(textures[tid].ncolors == 3);
	unsigned char *tex_data(new unsigned char[size]);

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			for (unsigned k = 0; k < 3; ++k) {
				tex_data[3*(i*width+j)+k] = 255*(((horiz ? i : j)&3) != 0);
			}
		}
	}
	textures[tid].data = tex_data;
}


void gen_hstripe_texture() {

	gen_stripe_texture(HSTRIPE_TEX, 1);
}


void gen_vstripe_texture() {

	gen_stripe_texture(VSTRIPE_TEX, 0);
}


void gen_tree_end_texture() {

	int const width(textures[TREE_END_TEX].width), height(textures[TREE_END_TEX].height);
	unsigned char *tex_data(textures[TREE_END_TEX].data);
	int const w2(width >> 1), h2(height >> 1);
	float const scale_vals[3] = {160.0, 130.0, 80.0};

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			float const rsin(sin(sqrt(float((j - w2)*(j - w2) + (i - h2)*(i - h2)))));
			float const darkness(0.9 + 0.1*rsin);
			int const offset(3*(i*width + j));

			for (unsigned k = 0; k < 3; ++k) {
				tex_data[offset+k] = (unsigned char)(scale_vals[k]*darkness + 20.0*signed_rand_float());
			}
		}
	}
}


void gen_blur_cent_texture() {

	int const width(textures[BLUR_CENT_TEX].width), height(textures[BLUR_CENT_TEX].height);
	assert(textures[BLUR_CENT_TEX].ncolors == 4);
	unsigned char *tex_data(new unsigned char[4*width*height]);
	int const w2(width >> 1), h2(height >> 1);
	float const scale(2.0/min(width, height));

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			float const radius(sqrt(float((j - w2)*(j - w2) + (i - h2)*(i - h2)))*scale);
			int const offset(4*(i*width + j));
			for (unsigned k = 0; k < 3; ++k) tex_data[offset+k] = 255;
			tex_data[offset+3] = (unsigned char)(255.0*(1.0 - CLIP_TO_01(radius))); // linear scaling for alpha
		}
	}
	textures[BLUR_CENT_TEX].data = tex_data;
}


void gen_gradient_texture() { // for horizon

	assert(textures[GRADIENT_TEX].width == 1); // 1D
	int const size(textures[GRADIENT_TEX].height);
	assert(textures[GRADIENT_TEX].ncolors == 4);
	unsigned char *tex_data(new unsigned char[4*size]);

	for (int i = 0; i < size; ++i) {
		for (unsigned k = 0; k < 3; ++k) tex_data[4*i+k] = 255;
		tex_data[4*i+3] = (unsigned char)max(0, (255 * 2 * (size/2 - abs(i - (int)size/2)) / size)); // linear gradient
	}
	textures[GRADIENT_TEX].data = tex_data;
}


void gen_wind_texture() {

	int const width(textures[WIND_TEX].width), height(textures[WIND_TEX].height), size(width*height);
	unsigned char *tex_data2(textures[CLOUD_RAW_TEX].data);
	assert(textures[WIND_TEX].ncolors == 1 && textures[CLOUD_RAW_TEX].ncolors == 4); // RGBA => grayscale luminance
	assert(tex_data2 != NULL && width == textures[CLOUD_RAW_TEX].width && height == textures[CLOUD_RAW_TEX].height);
	unsigned char *tex_data(new unsigned char[size]);

	for (int i = 0; i < size; ++i) {
		tex_data[i] = tex_data2[(i<<2)+3]; // put alpha in luminance
	}
	textures[WIND_TEX].data = tex_data;
}


void alloc_texture(int id) {

	assert(id < NUM_TEXTURES);
	textures[id].alloc();
}


colorRGBA get_landscape_texture_color(int xpos, int ypos) {

	if (cached_ls_colors.empty()) cached_ls_colors.resize(XY_MULT_SIZE, ALPHA0);
	unsigned const ix(xpos + MESH_X_SIZE*ypos);
	assert(ix < cached_ls_colors.size());
	if (cached_ls_colors[ix].alpha != 0.0) return cached_ls_colors[ix];
	assert(!point_outside_mesh(xpos, ypos));
	assert(textures[LANDSCAPE_TEX].ncolors == 3);
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	unsigned char *tex_data(textures[LANDSCAPE_TEX].data);
	unsigned const xstep(width/MESH_X_SIZE), ystep(height/MESH_Y_SIZE);
	unsigned const x0(xpos*xstep), y0(ypos*ystep), x1(x0 + xstep), y1(y0 + ystep);
	colorRGBA color(BLACK);

	for (unsigned y = y0; y < y1; ++y) { // like creating a mipmap
		for (unsigned x = x0; x < x1; ++x) {
			UNROLL_3X(color[i_] += tex_data[3*(width*y + x) + i_];)
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
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	unsigned char *tex_data(textures[LANDSCAPE_TEX].data);
	assert(textures[LANDSCAPE_TEX].ncolors == 3);
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
	int const suc3(3*sizeof(unsigned char)), id0(lttex[NTEXm1].id);
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
			memmove((tex_data+offset+3*j00), (tex_data+3*((j00+tox0)+lly*width)), (j01-j00)*suc3); // range could be overlapping
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
			texture const &t1(textures[id]);
			int const tof(t1.ncolors*(((i+toy)&(t1.height-1))*t1.width + ((j+tox)&(t1.width-1))));

			if (k1 == k2) { // single texture
				RGB_BLOCK_COPY((tex_data + o2), (t1.data + tof));
			}
			else { // blend two textures - performance critical
				texture const &t2(textures[id2]);
				int const tof2(t2.ncolors*(((i+toy)&(t2.height-1))*t2.width + ((j+tox)&(t2.width-1))));
				BLEND_COLOR((tex_data + o2), (t2.data + tof2), (t1.data + tof), t);
			}

			// handle steep slopes (dirt/rock texture replaces grass texture)
			float const sthresh[2][2] = {{0.45, 0.7}, {0.3, 0.55}};
			float const *const sti(sthresh[island]);
			float const vnz00(vertex_normals[ypos][xpos].z);

			if (vnz00 < sti[1]+0.1) {
				float const vnz01(vertex_normals[ypos][xpos1].z), vnz10(vertex_normals[ypos1][xpos].z), vnz11(vertex_normals[ypos1][xpos1].z);
				float const vnz((1.0 - xpi)*((1.0 - ypi)*vnz00 + ypi*vnz10) + xpi*((1.0 - ypi)*vnz01 + ypi*vnz11));

				if ((id == GROUND_TEX || id2 == GROUND_TEX) && vnz < sti[1]) { // ground/grass
					texture const &ta(textures[DIRT_TEX]);
					int const tofa(ta.ncolors*(((i+toy)&(ta.height-1))*ta.width + ((j+tox)&(ta.width-1))));
					unsigned char temp[3];

					if (id == GROUND_TEX || id2 == ROCK_TEX) {
						texture const &tb(textures[ROCK_TEX]);
						int const tofb(tb.ncolors*(((i+toy)&(tb.height-1))*tb.width + ((j+tox)&(tb.width-1))));
						BLEND_COLOR(temp, (tb.data+tofb), (ta.data+tofa), t);
					}
					else {
						RGB_BLOCK_COPY(temp, (ta.data+tofa));
					}
					float const val(CLIP_TO_01((vnz - sti[0])/(sti[1] - sti[0])));
					BLEND_COLOR((tex_data+o2), (tex_data+o2), temp, val);
				}
				else if (id2 == SNOW_TEX && vnz < sti[1]) { // snow
					texture const &ta(textures[ROCK_TEX]);
					int const tofa(ta.ncolors*(((i+toy)&(ta.height-1))*ta.width + ((j+tox)&(ta.width-1))));
					float const val(CLIP_TO_01(2.0f*(vnz - sti[0])/(sti[1] - sti[0])));
					BLEND_COLOR((tex_data+o2), (tex_data+o2), (ta.data+tofa), val);
				}
			}
		} // for j
	} // for i
	PRINT_TIME(" Data Gen");

	if (!read_landscape) {
		if (landscape0 == NULL || !scroll) { // initialize/copy entire texture
			int const totsize(3*width*height);
			if (landscape0 == NULL) landscape0 = new unsigned char[totsize];
			memcpy(landscape0, tex_data, totsize*sizeof(unsigned char));
			ls0_invalid = 0;
			PRINT_TIME(" Landscape0 Gen");
		}
		else if (lchanged0) { // create landscape0
			for (int i = i0; i != i1; i += di) {
				int const lly(i + toy0), off(3*i*width);

				if (lly >= 0 && lly < hy) { // copy part of this row from old landscape0
					memmove(landscape0+off+3*j00, landscape0+3*((j00+tox0)+lly*width), (j01-j00)*suc3); // range could be overlapping
					if (-tox0 > 0)      memcpy(landscape0+off, tex_data+off, -tox0*suc3); // copy from beginning
					if (width-wxtx > 0) memcpy(landscape0+off+wxtx3, tex_data+off+wxtx3, (width-wxtx)*suc3); // copy to end
				}
				else { // new section - copy entire row from tex_data
					memcpy(landscape0+off, tex_data+off, width*suc3);
				}
			}
			ls0_invalid = 0;
			PRINT_TIME(" Landscape0 Copy");
		}
		else { // delay creation of landscape0 until it's needed
			ls0_invalid = 1;
		}
	}
	textures[LANDSCAPE_TEX].gl_delete(); // should we try to update rather than recreating from scratch?
	init_texture(LANDSCAPE_TEX); // performance bottleneck
	PRINT_TIME(" Final");
}


void regrow_landscape_texture_amt0() {

	//RESET_TIME;
	static int counter(0);
	if (read_landscape) return;
	if (ls0_invalid) create_landscape_texture();
	assert(textures[LANDSCAPE_TEX].data != NULL && landscape0 != NULL);
	int const regen_bshift(max(1, min(7, int(log(1.0/max(0.0001f, LANDSCAPE_REGEN_AMT))/log(2.0)))));
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	unsigned char *tex_data(textures[LANDSCAPE_TEX].data);
	int const y1((counter%LANDSCAPE_REGEN_MOD)*(height/LANDSCAPE_REGEN_MOD));
	int const y2(y1 + (height/LANDSCAPE_REGEN_MOD));
	assert(y2 <= height);

	if (LANDSCAPE_REGEN_AMT > 0.0 || skip_regrow) {
		for (int i = y1; i < y2; ++i) {
			int const i_step(i*width);

			for (int j = 0; j < width; ++j) { // performance critical
				int const offset(3*(i_step + j));
				UNROLL_3X(tex_data[offset+i_] += (unsigned char)(((int)landscape0[offset+i_] - (int)tex_data[offset+i_]) >> regen_bshift);)
			}
		}
	}
	++counter;
	skip_regrow = 0;
	check_init_texture(LANDSCAPE_TEX);
	glBindTexture(GL_TEXTURE_2D, textures[LANDSCAPE_TEX].tid);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, y1, width, (y2-y1), GL_RGB, GL_UNSIGNED_BYTE, (tex_data + 3*width*y1));
	//PRINT_TIME("Regrow");
}


float add_crater_to_landscape_texture(float xval, float yval, float radius) { // no grass underwater so don't check

	int const xpos(get_xpos(xval)), ypos(get_ypos(yval));
	if (point_outside_mesh(xpos, ypos)) return 0.0; // off the terrain area
	texture const &t1(textures[get_bare_ls_tid(mesh_height[ypos][xpos])]);
	int const tidw(t1.width), tidh(t1.height), ncol(t1.ncolors);
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	unsigned char *data(t1.data), *tex_data(textures[LANDSCAPE_TEX].data);
	float const xscale(((float)MESH_X_SIZE)/((float)width)), yscale(((float)MESH_Y_SIZE)/((float)height));
	int const xpos2(int((xval + X_SCENE_SIZE)/(xscale*DX_VAL) + 0.5));
	int const ypos2(int((yval + Y_SCENE_SIZE)/(yscale*DY_VAL) + 0.5));
	if (xpos2 < 0 || ypos2 < 0 || xpos2 >= width || ypos2 >= height) return 0.0;
	int const rad(max(0, int(radius/((xscale + yscale)*(DX_VAL + DY_VAL))) - 1)), radsq(rad*rad); // size in texture space
	int const x1(max(0, xpos2-rad)), y1(max(0, ypos2-rad));
	int const x2(min(width-1,  xpos2+rad)), y2(min(height-1, ypos2+rad));

	for (int i = y1; i <= y2; ++i) {
		int const offset(i*width), yterm((i - ypos2)*(i - ypos2));

		for (int j = x1; j <= x2; ++j) {
			if ((yterm + (j - xpos2)*(j - xpos2)) < radsq) {
				int const o2(3*(offset + j)), tof(ncol*((i&(tidh-1))*tidw + (j&(tidw-1))));
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
	texture const &t1(textures[get_bare_ls_tid(mesh_height[ypos][xpos])]);
	int const tidw(t1.width), tidh(t1.height), ncol(t1.ncolors);
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	unsigned char *data(t1.data), *tex_data(textures[LANDSCAPE_TEX].data);
	float const xdiv(((((float)MESH_X_SIZE)/((float)width))*DX_VAL)), ydiv(((((float)MESH_Y_SIZE)/((float)height))*DY_VAL));
	int const xpos1(max(0,      int((get_xval(xpos)   + X_SCENE_SIZE)/xdiv + 0.5)));
	int const ypos1(max(0,      int((get_yval(ypos)   + Y_SCENE_SIZE)/ydiv + 0.5)));
	int const xpos2(min(width,  int((get_xval(xpos+1) + X_SCENE_SIZE)/xdiv + 0.5)));
	int const ypos2(min(height, int((get_yval(ypos+1) + Y_SCENE_SIZE)/ydiv + 0.5)));

	for (int i = ypos1; i < ypos2; ++i) {
		unsigned const iw(i*width), ival((i&(tidh-1))*tidw);

		for (int j = xpos1; j < xpos2; ++j) {
			int const o2(3*(iw + j)), tof(ncol*(ival + (j&(tidw-1))));
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


void add_cutout_to_landscape_texture(int x1, int y1, int x2, int y2) { // unused

	if (read_landscape) return;
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	unsigned char *tex_data(textures[LANDSCAPE_TEX].data);
	x1 = (width*max(x1, 0))/MESH_X_SIZE;
	x2 = (width*min(x2, MESH_X_SIZE))/MESH_X_SIZE;
	y1 = (height*max(y1, 0))/MESH_Y_SIZE;
	y2 = (height*min(y2, MESH_Y_SIZE))/MESH_Y_SIZE;

	for (int i = y1; i < y2; ++i) {
		int const offset(i*width), yv((i*MESH_Y_SIZE)/height);

		for (int j = x1; j < x2; ++j) {
			int const xv((j*MESH_X_SIZE)/width);
			assert(xv < MESH_X_SIZE && yv < MESH_Y_SIZE);
			texture const &t1(textures[get_bare_ls_tid(mesh_height[yv][xv])]);
			int const o2(3*(offset + j)), tof(t1.ncolors*((i&(t1.height-1))*t1.width + (j&(t1.width-1))));
			RGB_BLOCK_COPY((tex_data+o2), (t1.data+tof));
		}
	}
}


void add_color_to_landscape_texture(colorRGBA const &color, float xval, float yval, float radius, int check_unique) {

	int const xpos0(get_xpos(xval)), ypos0(get_ypos(yval));
	if (point_outside_mesh(xpos0, ypos0)) return; // off the terrain area
	if (wminside[ypos0][xpos0] && water_matrix[ypos0][xpos0] > mesh_height[ypos0][xpos0]) return; // underwater
	int const index(xpos0 + MESH_X_SIZE*ypos0);
	if (ls_color_texels.find(index) != ls_color_texels.end()) return;
	ls_color_texels.insert(index);
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	unsigned char *tex_data(textures[LANDSCAPE_TEX].data);
	float const xscale(((float)MESH_X_SIZE)/((float)width)), yscale(((float)MESH_Y_SIZE)/((float)height));
	int const maxsize(3*width*height);
	int const xpos(int((xval + X_SCENE_SIZE)/(((float)xscale)*DX_VAL) + 0.5));
	int const ypos(int((yval + Y_SCENE_SIZE)/(((float)yscale)*DY_VAL) + 0.5));
	if (xpos < 0 || ypos < 0 || xpos >= width || ypos >= height) return;
	int const rad(max(0, int(5.0*radius/((xscale + yscale)*(DX_VAL + DY_VAL))) - 1)), radsq(max(1, rad*rad)); // size in texture space
	int const x1(max(0, xpos-rad)), y1(max(0, ypos-rad)), x2(min(width-1, xpos+rad)), y2(min(height-1, ypos+rad));
	unsigned char color_i[3];
	unpack_color(color_i, color);

	for (int i = y1; i <= y2; ++i) {
		int const offset(i*width), yterm((i - ypos)*(i - ypos));

		for (int j = x1; j <= x2; ++j) {
			int const dist_sq((yterm + (j - xpos)*(j - xpos)));

			if (dist_sq < radsq) {
				assert(offset + j <= width*height);
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
	int const width(textures[LANDSCAPE_TEX].width), height(textures[LANDSCAPE_TEX].height);
	int const tx(int(((float)width /(float)MESH_X_SIZE)*(pos.x + X_SCENE_SIZE)*DX_VAL_INV + 0.5));
	int const ty(int(((float)height/(float)MESH_Y_SIZE)*(pos.y + Y_SCENE_SIZE)*DY_VAL_INV + 0.5));
	if (tx < 0 || tx >= width || ty < 0 || ty >= height) return;
	skip_regrow = 1;
	unsigned char *tex_data(textures[LANDSCAPE_TEX].data);
	acc *= TEXTURE_SNOW_ACC_RATE;
	float const acc255(255.0*acc);
	int const rsq(SNOW_ACC_RADIUS*SNOW_ACC_RADIUS);
	float const frad((float)rsq), inv_frad(1.0/frad);
	int const x1(max(0, tx-SNOW_ACC_RADIUS)), x2(min(width -1, tx+SNOW_ACC_RADIUS));
	int const y1(max(0, ty-SNOW_ACC_RADIUS)), y2(min(height-1, ty+SNOW_ACC_RADIUS));

	for (int i = y1; i <= y2; ++i) {
		int const dy(i*width), ysq((ty - i)*(ty - i));
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
	texture const &t1(textures[LANDSCAPE_TEX]);
	unsigned const nc(t1.ncolors);
	assert(nc == 3);
	assert(t1.width > 0 && t1.height > 0 && t1.data != NULL);
	int const width(t1.width), height(t1.height);

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
	unsigned char *copy(NULL), *data(t1.data + offset);

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
	bind_2d_texture(t1.tid);
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


void set_texgen_vec4(float const v[4], bool s_or_t, bool as_attr, bool enable_and_set_mode) {

	if (as_attr) {
		add_attrib_float_array(1+s_or_t, v, 4);
	}
	else {
		if (enable_and_set_mode) {
			glEnable(s_or_t ? GL_TEXTURE_GEN_T : GL_TEXTURE_GEN_S);
			glTexGeni((s_or_t ? GL_T : GL_S), GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
		}
		glTexGenfv((s_or_t ? GL_T : GL_S), GL_EYE_PLANE, v);
	}
}


void setup_texgen_full(float sx, float sy, float sz, float sw, float tx, float ty, float tz, float tw, int mode, bool as_attr) {

	float const tex_param_s[4] = {sx, sy, sz, sw};
	float const tex_param_t[4] = {tx, ty, tz, tw};
	set_texgen_vec4(tex_param_s, 0, as_attr, 1);
	set_texgen_vec4(tex_param_t, 1, as_attr, 1);
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
}


void setup_polygon_texgen(vector3d const &norm, float const scale[2], float const xlate[2], bool swap_txy, bool as_attr) {

	int const d0(get_min_dim(norm));
	vector3d v[2] = {all_zeros, all_zeros};
	v[0][d0] = 1.0;
	cross_product(norm, v[0], v[1]);
	cross_product(norm, v[1], v[0]);

	for (unsigned i = 0; i < 2; ++i) { // ignoring xoff2/yoff2
		float const tex_param[4] = {scale[i]*v[i].x, scale[i]*v[i].y, scale[i]*v[i].z, xlate[i]};
		set_texgen_vec4(tex_param, ((i != 0) ^ swap_txy), as_attr, 1);
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


float get_texture_component(unsigned tid, float xval, float yval, int comp) {

	assert(tid < NUM_TEXTURES);
	texture const &tex(textures[tid]);
	assert(comp < tex.ncolors);
	int tx(int(tex.width *xval) & (tex.width -1)); // width and height are a power of 2
	int ty(int(tex.height*yval) & (tex.height-1));
	if (tx < 0) tx += tex.width;
	if (ty < 0) ty += tex.height;
	assert(tx >= 0 && ty >= 0 && tx < tex.width && ty < tex.height);
	return tex.data[tex.ncolors*(tex.width*ty + tx) + comp]/255.0;
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




