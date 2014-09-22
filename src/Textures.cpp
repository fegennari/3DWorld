// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/25/02
#include "3DWorld.h"
#include "mesh.h"
#include "sinf.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"
#include "shaders.h"


float const TEXTURE_SMOOTH        = 0.01;
float const TEXTURE_SMOOTH_I      = 0.1;
float const SPHERE_SECTION        = 0.75;
float const TEXTURE_SNOW_ACC_RATE = 0.05;
int   const SNOW_ACC_RADIUS       = 3;
float const LANDSCAPE_REGEN_AMT   = 0.02; // how much regen after LANDSCAPE_REGEN_MOD ticks
int   const LANDSCAPE_REGEN_MOD   = 32;   // ticks for entire mesh regen pass
int   const LANDSCAPE_REGEN_RATE  = 400;  // ticks per landscape texture update (slow)
int   const MAX_UPDATE_SIZE       = 16;
float const ADJ_VALUE             = 1.0;
float const LS_TEX_ANISO          = 2.0; // value of anisotropic texture filtering for landscape source textures

bool const RELOAD_TEX_ON_HOLE  = 0;
bool const LANDSCAPE_MIPMAP    = 0; // looks better, but texture update doesn't recompute the mipmaps
bool const SHOW_TEXTURE_MEMORY = 0;
bool const COMPRESS_TEXTURES   = 1;
bool const CHECK_FOR_LUM       = 1;
float const SMOOTH_SKY_POLES   = 0.2;



//GROUND ROCK WATER WATER2 SKY SUN MOON EARTH ICE SNOW LEAF WOOD SAND ROCK2 CAMOFLAGE GRASS PALM SMOKE PLASMA GEN LANDSCAPE TREE_END TREE_SNOW TREE_HEMI ...
//0      1    2     3      4   5   6    7     8   9    10   11   12   13    14        15    16   17    18     19  20        21       22        23        ...

texture_t textures[NUM_TEXTURES] = { // 4 colors without wrap sometimes has a bad transparent strip on spheres
// type: 0 = read from file, 1 = generated, 2 generated and dynamically updated
// format: 0 = RGB RAW, 1 = BMP, 2 = RGB RAW, 3 = RGBA RAW, 4: targa (*tga), 5: jpeg, 6: png, 7: auto, 8: tiff, 9: generate (not loaded from file), 10: DDS
// use_mipmaps: 0 = none, 1 = standard OpenGL, 2 = openGL + CPU data, 3 = custom alpha OpenGL
// type format width height wrap ncolors use_mipmaps name [invert_y=0 [do_compress=1 [anisotropy=1.0 [mipmap_alpha_weight=1.0]]]]
//texture_t(0, 6, 512,  512,  1, 3, 0, "ground.png"),
texture_t(0, 6, 128,  128,  1, 3, 2, "grass.png", 0, 1, LS_TEX_ANISO), // mipmap for small trees?
texture_t(0, 6, 256,  256,  1, 3, 1, "rock.png"),
texture_t(0, 5, 512,  512,  1, 3, 1, "water.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "stucco.jpg", 0, 0), // compression is slow
texture_t(0, 5, 0,    0,    1, 4, 0, "sky.jpg", 1), // 1024x1024
texture_t(0, 5, 0,    0,    1, 3, 1, "brick1.jpg"), // brick2?
texture_t(0, 5, 0,    0,    1, 3, 1, "moon.jpg"),
texture_t(0, 6, 256,  256,  0, 3, 1, "earth.png", 1),
texture_t(0, 5, 0,    0,    1, 3, 1, "marble.jpg", 0, 0), // or marble2.jpg, compression is slow
texture_t(0, 7, 0,    0,    1, 3, 2, "snow2.jpg", 0, 1, LS_TEX_ANISO),
texture_t(0, 5, 0,    0,    0, 4, 3, "leaf.jpg", 1, 1, 4.0), // 128x128
texture_t(0, 5, 0,    0,    1, 3, 1, "bark2.jpg"), // bark.jpg: 224x224, bark2.jpg: 400x400 (Note: must match baseball bat texture size)
texture_t(0, 5, 512,  512,  1, 3, 2, "desert_sand.jpg", 0, 1, LS_TEX_ANISO),
texture_t(0, 6, 256,  256,  1, 3, 2, "rock2.png", 0, 1, LS_TEX_ANISO),
texture_t(0, 5, 512,  512,  1, 3, 1, "camoflage.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "hedges.jpg", 0, 0), // 1024x1024, compression is slow
texture_t(0, 1, 512,  512,  1, 3, 1, "brick1.bmp", 0, 1, 8.0),
texture_t(0, 0, 512,  512,  1, 3, 1, "manhole.bmp", 1),
texture_t(0, 6, 128,  128,  1, 4, 3, "palmtree.png", 1),
texture_t(1, 9, 256,  256,  1, 4, 1, "@smoke"),  // not real file
texture_t(1, 9, 64,   64,   1, 4, 1, "@plasma"), // not real file
texture_t(1, 9, 128,  128,  0, 3, 0, "@gen"),    // not real file - unused
texture_t(2, 7, 1024, 1024, 0, 3, LANDSCAPE_MIPMAP, "@landscape_tex"), // for loading real landscape texture
texture_t(1, 9, 128,  128,  0, 3, 0, "@tree_end"),  // not real file
texture_t(1, 9, 1024, 1024, 1, 4, 1, "@tree_hemi", 0, 1), // not real file, compression is too slow, mipmap for trees?
texture_t(0, 5, 0  ,  0,    1, 3, 1, "shingles.jpg", 0, 0, 8.0), // compression is slow
texture_t(0, 6, 256,  256,  1, 3, 1, "paneling.png", 0, 1, 16.0),
texture_t(0, 6, 256,  256,  1, 3, 1, "cblock.png", 0, 1, 8.0),
texture_t(0, 5, 0,    0,    0, 4, 3, "mj_leaf.jpg", 1), // 128x128
texture_t(0, 5, 0,    0,    0, 4, 3, "live_oak.jpg", 1, 1, 4.0), // 80x128
texture_t(0, 5, 0,    0,    0, 4, 3, "leaf2.jpg", 1, 1, 4.0), // 212x256
texture_t(0, 5, 0,    0,    0, 4, 3, "leaf3c.jpg", 1, 1, 4.0), // 208x256
texture_t(0, 5, 0,    0,    0, 4, 3, "plant1.jpg", 1), // 256x256
texture_t(0, 6, 256,  256,  0, 4, 3, "plant2.png", 1),
//texture_t(0, 5, 0,    0,    0, 4, 3, "plant2.jpg", 1), // 160x256
texture_t(0, 6, 256,  256,  0, 4, 3, "plant3.png", 1),
//texture_t(0, 5, 0,    0,    0, 4, 3, "plant3.jpg", 1), // 176x256
texture_t(0, 5, 0,    0,    0, 4, 3, "hibiscus.jpg", 1), // 64x64
texture_t(0, 5, 0,    0,    1, 3, 1, "fence.jpg", 0, 0, 8.0), // 896x896, compression is slow
texture_t(0, 6, 128,  128,  1, 3, 1, "skull.png"),
texture_t(0, 6, 64,   64,   1, 3, 1, "radiation.png", 1),
texture_t(0, 6, 128,  128,  1, 3, 1, "yuck.png"),
texture_t(0, 6, 256,  256,  0, 4, 0, "sawblade.png", 1),
texture_t(0, 6, 256,  256,  0, 4, 0, "sawblade_b.png", 1),
texture_t(0, 6, 256,  256,  0, 4, 1, "blur.png"),
texture_t(0, 6, 256,  256,  1, 4, 1, "blur_s.png"),
texture_t(0, 6, 256,  256,  0, 4, 3, "pine.png", 1, 1, 1.0, 0.36),
//texture_t(0, 5, 0,    0,    0, 4, 3, "pine.jpg", 1, 1, 1.0, 0.36), // 184x256
texture_t(0, 6, 128,  128,  1, 3, 1, "noise.png"),
texture_t(0, 5, 0,    0,    1, 3, 1, "wood.jpg", 0, 0), // 768x768, compression is slow
texture_t(0, 6, 128,  128,  1, 3, 1, "hb_brick.png", 0, 1, 8.0),
texture_t(0, 6, 128,  128,  1, 3, 1, "particleb.png", 0, 1, 8.0),
texture_t(0, 6, 128,  128,  1, 3, 1, "plaster.png"),
texture_t(0, 6, 256,  256,  1, 3, 1, "tile.png", 0, 1, 8.0),
texture_t(0, 6, 256,  32,   1, 3, 1, "CommandCAD.png"),
texture_t(1, 9, 32,   32,   1, 1, 1, "@disint"),   // not real file
texture_t(1, 9, 256,  256,  1, 4, 1, "@blur_inv"), // not real file
texture_t(1, 9, 32,   32,   1, 3, 0, "@hstripe", 0, 1, 8.0), // not real file
texture_t(1, 9, 32,   32,   1, 3, 0, "@vstripe", 0, 1, 8.0), // not real file
texture_t(0, 5, 512,  512,  1, 3, 1, "bcube.jpg"),
texture_t(0, 6, 0  ,  0  ,  0, 4, 1, "atlas/explosion.png", 1),
texture_t(0, 5, 512,  512,  1, 3, 1, "shiphull.jpg"),
texture_t(0, 5, 512,  512,  1, 3, 1, "bcube2.jpg"),
texture_t(0, 5, 512,  512,  1, 3, 1, "bcube_tactical.jpg"),
texture_t(0, 6, 512,  256,  1, 3, 1, "rock_sphere.png"),
texture_t(0, 6, 0,    0,    0, 4, 3, "papaya_leaf.png", 1, 1, 4.0), // 256x256
//texture_t(0, 6, 256,  256,  0, 4, 3, "coffee_leaf.png"), // half the texture is wasted, but leaves must be square (for now)
texture_t(0, 6, 0,    0,    0, 4, 3, "coffee_leaf.png", 1), // 256x256
texture_t(0, 6, 256,  256,  1, 4, 0, "smiley_skull.png", 1),
texture_t(0, 5, 512,  512,  1, 3, 1, "ice.2.jpg"),
texture_t(0, 6, 256,  256,  1, 3, 2, "rock.03.png", 0, 1, LS_TEX_ANISO),
texture_t(0, 6, 16,   16,   1, 3, 0, "black.png"),
texture_t(0, 6, 16,   16,   1, 3, 0, "white.png"),
texture_t(0, 6, 0  ,  0  ,  0, 4, 0, "atlas/fire.png"),
texture_t(0, 5, 0,    0,    1, 4, 0, "sky.jpg", 1), // 1024x1024
texture_t(0, 6, 256,  256,  0, 4, 0, "snowflake.png", 1),
texture_t(1, 9, 128,  128,  0, 4, 1, "@blur_center"), // not real file
texture_t(1, 9, 1,    128,  1, 4, 0, "@gradient"), // not real file
texture_t(0, 6, 1024, 128,  0, 3, 1, "grass_blade.png", 1),
texture_t(1, 9, 1024, 1024, 1, 1, 1, "@wind_texture"),  // not real file
texture_t(0, 5, 0,    0,    1, 3, 1, "mossy_rock.jpg"), // 500x500
// bark
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark1.jpg"), // 600x600
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark2.jpg"), // 512x512
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark2-normal.jpg", 0, 0, 4.0), // 512x512
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark_lendrick.jpg", 0, 0), // 892x892, compression is slow
texture_t(0, 6, 0,    0,    1, 3, 1, "bark/bark_lylejk.png", 0, 0), // 1024x768, compression is slow
// normal/caustic maps
texture_t(0, 4, 0,    0,    1, 3, 1, "normal_maps/water_normal.tga", 0, 1, 8.0), // 512x512
texture_t(0, 6, 0,    0,    1, 3, 1, "normal_maps/ocean_water_normal.png", 0, 0, 4.0), // 1024x1024 (Note: compression disabled as it causes artifacts)
texture_t(0, 5, 0,    0,    1, 3, 1, "caustics.jpg"), // 512x512
// noise
texture_t(0, 6, 0,    0,    1, 1, 0, "perlin_simplex.png"), // 256x256
texture_t(1, 9, 128,  128,  1, 1, 0, "@noise_gen"), // not real file
texture_t(1, 9, 128,  128,  1, 1, 1, "@noise_gen_mipmap"), // not real file
texture_t(1, 9, 256,  256,  1, 1, 1, "@noise_gen_sparse"), // not real file
texture_t(1, 9, 400,  400,  1, 3, 1, "@player_bbb_tex"), // not real file (Note: must match bark texture size)
texture_t(0, 5, 0,    0,    0, 4, 3, "pine_tree_leaves.jpg", 1, 0, 1.0, 0.28), // 256x256
texture_t(0, 5, 0,    0,    0, 4, 1, "flare1.jpg"), // 384x384
texture_t(0, 5, 0,    0,    0, 4, 1, "flare2.jpg"), // 128x128 (Note: low resolution)
texture_t(0, 5, 0,    0,    0, 4, 1, "Flare3.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 4, 1, "flare4.jpg"), // 256x256
//texture_t(0, 5, 0,    0,    0, 4, 1, "flare5.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 4, 1, "flare5b.jpg"), // 416x416
texture_t(0, 5, 0,    0,    1, 3, 1, "foam1.jpg"), // 512x512

texture_t(0, 5, 0,    0,    0, 3, 0, "bullet_hole/bullet_diffuse.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 1, 0, "bullet_hole/bullet_alpha.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 3, 0, "bullet_hole/bullet_normal.jpg"), // 256x256
texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/sand_normal.jpg", 0, 0, 2.0), // (Note: compression disabled as it causes artifacts)

texture_t(0, 5, 0,    0,    1, 3, 1, "raindrop_dots.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "spaceship1.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "spaceship2.jpg"),
texture_t(0, 6, 0,    0,    0, 4, 1, "atlas/blood.png"),
texture_t(0, 5, 0,    0,    1, 3, 1, "lichen.jpg", 0, 0), // 1500x1500, compression is probably slow
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/palm_bark.jpg"), // 512x512
texture_t(0, 5, 0,    0,    0, 4, 3, "daisy.jpg", 0, 1, 4.0), // 1024x1024
//texture_t(0, 4, 0,    0,    1, 3, 1, "../Sponza2/textures/spnza_bricks_a_diff.tga")
// type format width height wrap ncolors use_mipmaps name [invert_y=0 [do_compress=1 [anisotropy=1.0 [mipmap_alpha_weight=1.0]]]]
};


// zval should depend on def_water_level and temperature
float h_dirt[NTEX_DIRT], clip_hd1;
std::set<int> ls_color_texels;
vector<colorRGBA> cached_ls_colors;

typedef map<pair<unsigned, unsigned>, unsigned> texture_map_t;
texture_map_t noise_tex_3ds;

typedef map<string, unsigned> name_map_t;
name_map_t texture_name_map;

int landscape_changed(0), lchanged0(0), skip_regrow(0), ltx1(0), lty1(0), ltx2(0), lty2(0), ls0_invalid(1);
unsigned char *landscape0 = NULL;


extern bool mesh_difuse_tex_comp, water_is_lava;
extern unsigned smoke_tid, dl_tid, elem_tid, gb_tid, reflection_tid;
extern int world_mode, read_landscape, default_ground_tex, xoff2, yoff2, DISABLE_WATER;
extern int scrolling, dx_scroll, dy_scroll, display_mode, iticks, universe_only;
extern float zmax, zmin, glaciate_exp, relh_adj_tex, vegetation, fticks;
extern char *mesh_diffuse_tex_fn;


void gen_smoke_texture();
void gen_plasma_texture();
void gen_disintegrate_texture();
void gen_tree_hemi_texture();
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


bool using_custom_landscape_texture() {return (read_landscape && mesh_diffuse_tex_fn != NULL);}


void set_landscape_texture_from_file() {

	assert(mesh_diffuse_tex_fn != NULL);
	texture_t &tex(textures[LANDSCAPE_TEX]);
	tex.name  = mesh_diffuse_tex_fn;
	tex.type  = 0; // loaded from file
	tex.width = tex.height = 0; // reset to 0 for autodetect (only works with some formats)
	tex.use_mipmaps = 1;
	tex.do_compress = mesh_difuse_tex_comp; // compression is slow (one time cost), but saves GPU memory
}


void load_textures() {

	cout << "loading textures";
	if (using_custom_landscape_texture()) {set_landscape_texture_from_file();} // must be done first
	load_texture_names();

	for (int i = 0; i < NUM_TEXTURES; ++i) {
		cout.flush();
		cout << ".";
		if (!is_tex_disabled(i)) textures[i].load(i);
	}
	cout << endl;
	textures[BULLET_D_TEX].merge_in_alpha_channel(textures[BULLET_A_TEX]);
	gen_smoke_texture();
	gen_plasma_texture();
	gen_disintegrate_texture();
	update_player_bbb_texture(0.0, 0);
	gen_blur_inv_texture(); // must be after BLUR_TEX
	gen_stripe_texture(HSTRIPE_TEX, 1);
	gen_stripe_texture(VSTRIPE_TEX, 0);
	gen_blur_cent_texture();
	gen_gradient_texture();
	gen_noise_texture();
	load_font_texture_atlas();

	if (!universe_only) {
		gen_wind_texture();
		gen_tree_hemi_texture();
		gen_tree_end_texture();
	}
	for (int i = 0; i < NUM_TEXTURES; ++i) {
		if (is_tex_disabled(i)) continue; // skip
		textures[i].init();
	}
	textures[TREE_HEMI_TEX].set_color_alpha_to_one();

	int max_tiu(0), max_ctiu(0);
	glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &max_tiu);
	glGetIntegerv(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS, &max_ctiu);
	cout << "max TIUs: " << max_tiu << ", max combined TIUs: " << max_ctiu << endl;
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
}


bool select_texture(int id) {

	bool const no_tex(id < 0);
	if (no_tex) {id = WHITE_TEX;} //glBindTexture(GL_TEXTURE_2D, 0); // bind to none
	assert(id < NUM_TEXTURES);
	check_init_texture(id);
	textures[id].bind_gl();
	return !no_tex;
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

	cout << "Freeing textures..." << endl; // only print if something was loaded?
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
	free_texture(reflection_tid);
	free_font_texture_atlas();

	for (texture_map_t::iterator i = noise_tex_3ds.begin(); i != noise_tex_3ds.end(); ++i) {
		free_texture(i->second);
	}
	noise_tex_3ds.clear();
}


void setup_landscape_tex_colors(colorRGBA const &c1, colorRGBA const &c2) { // c1 = high, c2 = low

	textures[ROCK_TEX].set_to_color(c1);
	textures[SAND_TEX].set_to_color(c2);
	textures[DIRT_TEX].set_to_color(c2);
}


void free_texture(unsigned &tid) {

	if (tid == 0) return; // avoid GL calls
	if (glIsTexture(tid)) {glDeleteTextures(1, &tid);}
	tid = 0;
}


void texture_t::alloc() {

	free_data();
	data = new unsigned char[num_bytes()];

	if (SHOW_TEXTURE_MEMORY) {
		static unsigned tmem(0);
		tmem += num_bytes();
		cout << "tex mem = " << tmem << endl;
	}
}


void texture_t::bind_gl() const {
	
	assert(tid > 0);
	bind_2d_texture(tid);
}


void texture_t::free_mm_data() {

	delete [] mm_data;
	mm_data = NULL;
	mm_offsets.clear();
}

void texture_t::free_data() {

	gl_delete();
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

	assert(ncolors >= 1 && ncolors <= 4);
	if (is_16_bit_gray) {return GL_R16;} // compressed?
	return get_internal_texture_format(ncolors, (COMPRESS_TEXTURES && do_compress && type != 2));
}

GLenum texture_t::calc_format() const {
	return (is_16_bit_gray ? GL_RED : get_texture_format(ncolors));
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
	setup_texture(tid, (use_mipmaps != 0 && !defer_load()), wrap, wrap, 0, 0, 0, anisotropy);

	if (defer_load()) {
		deferred_load_and_bind(); // FIXME: mipmaps?
	}
	else {
		assert(is_allocated() && width > 0 && height > 0);
		glTexImage2D(GL_TEXTURE_2D, 0, calc_internal_format(), width, height, 0, calc_format(), get_data_format(), data);
		if (use_mipmaps == 1 || use_mipmaps == 2) {gen_mipmaps();}
		if (use_mipmaps == 3) {create_custom_mipmaps();}
	}
	assert(glIsTexture(tid));
	//PRINT_TIME("Texture Init");
}


void texture_t::calc_color() { // incorrect in is_16_bit_gray mode

	if (defer_load() && !is_allocated()) {
		color = WHITE; // FIXME: can we do any better than this?
		return;
	}
	assert(is_allocated());
	float colors[4] = {0.0}, weight(0.0);
	unsigned const size(num_pixels());
	bool const has_alpha_comp(ncolors == 4);
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
				unsigned char const alpha(data[offset+3]);
				colors[3] += alpha;
				if (alpha > 0 && alpha < 255) {has_binary_alpha = 0;}
			}
		}
	}
	UNROLL_3X(color[i_] = colors[i_]/(255.0*weight);)
	color.alpha = (has_alpha_comp ? CLIP_TO_01(colors[3]/(255.0f*size)) : 1.0);
}


void texture_t::copy_alpha_from_texture(texture_t const &at, bool alpha_in_red_comp) {

	assert(at.is_allocated() && is_allocated()); // check that data is allocated in both textures
	assert(ncolors == 4); // check for alpha channel
	assert(!is_bound());  // check that texture isn't already bound
	assert(at.ncolors == 1 || at.ncolors == 3 || at.ncolors == 4);
	
	if (at.width != width || at.height != height) {
		resize(at.width, at.height);
		assert(at.width == width && at.height == height);
	}
	// alpha channel comes from either R or A in an RGBA texture, R in RGB texture, or R in grayscale texture
	unsigned const npixels(num_pixels()), alpha_offset((at.ncolors < 4 || alpha_in_red_comp) ? 0 : 3);

	for (unsigned i = 0; i < npixels; ++i) {
		data[4*i+3] = at.data[at.ncolors*i+alpha_offset]; // copy alpha values
	}
}


void texture_t::merge_in_alpha_channel(texture_t const &at) {

	assert(ncolors == 3 && at.ncolors == 1);
	assert(width == at.width && height == at.height);
	add_alpha_channel();
	copy_alpha_from_texture(at, 0);
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
		int const ret(gluScaleImage(format, tsz,   tsz,   get_data_format(), get_mipmap_data(level),
			                                tsz/2, tsz/2, get_data_format(), (mm_data + mm_offsets[level])));
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

	if (tid < 0) {return WHITE;} // ???
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


void texture_t::auto_insert_alpha_channel(int index) {

	int alpha_white(0);
	unsigned char alpha(255);
	unsigned const size(num_pixels());
	bool const is_alpha_mask(index == BLUR_TEX || index == SBLUR_TEX || index == BLUR_CENT_TEX || (index >= FLARE1_TEX && index <= FLARE5_TEX));
	bool const is_alpha_tex(index == EXPLOSION_TEX || index == FIRE_TEX || is_alpha_mask);
	bool has_zero_alpha(0);
	assert(is_allocated());

	for (unsigned i = 0; i < size; ++i) {
		int const i4(i << 2);
		unsigned char *buf(data+i4);

		if (index == CLOUD_TEX || index == CLOUD_RAW_TEX) {
			// white -> alpha = 255
			// blue  -> alpha = 0
			float const val(float(buf[0]) + float(buf[1]));
			alpha = ((val <= 340.0) ? 0 : ((unsigned char)1.0*(val - 340.0)));

			if (SMOOTH_SKY_POLES > 0.0 && index == CLOUD_TEX) {
				unsigned const y(i/width);
				float const dist(float(y)/float(height)), d2(min(dist, (1.0f - dist)));
				if (d2 < SMOOTH_SKY_POLES) {alpha = (unsigned char)(d2*d2*float(alpha)/(SMOOTH_SKY_POLES*SMOOTH_SKY_POLES));}
			}
		}
		else { // make white/black part transparent, for example leaves
			float const val(float(buf[0]) + float(buf[1]) + float(buf[2]));
				
			if (is_alpha_tex) { // animated/multipart textures
				alpha = (unsigned char)(0.333*val);
				if (is_alpha_mask) {buf[0] = buf[1] = buf[2] = 255;}
			}
			else {
				if (i == 0) { // key off of first (llc) pixel
					alpha_white = (index == SMILEY_SKULL_TEX) ? 0 : ((int)buf[0] + (int)buf[1] + (int)buf[2] > 400);
				}
				if (index == DAISY_TEX) {
					alpha = ((buf[0] == 255 && buf[1] == 255 && buf[2] == 255) ? 0 : 255); // all white = transparent
				}
				else if (alpha_white) {
					float const thresh((index == LEAF3_TEX) ? 700.0 : 600.0);
					alpha = ((val > thresh) ? 0 : ((val < thresh-100.0) ? 255 : (unsigned char)(2.55*(thresh - val))));
				}
				else {
					alpha = ((val < ((index == PINE_TEX || index == PINE_TREE_TEX) ? 65 : 32)) ? 0 : 255);
				}
				has_zero_alpha |= (alpha == 0);
			}
		}
		buf[3] = alpha;
	} // for i
	if (has_zero_alpha && !no_avg_color_alpha_fill) {
		calc_color(); // needed to ensure color is correct (may be recalculated later)
		unsigned char avg_rgb[3];
		UNROLL_3X(avg_rgb[i_] = (unsigned char)(255*color[i_]);)

		// set all alpha=0 texels to the average non-transparent color to improve mipmap quality
		for (unsigned i = 0; i < size; ++i) {
			int const i4(i << 2);
			if (data[i4+3] == 0) {RGB_BLOCK_COPY((data+i4), avg_rgb);} // alpha == 0
		}
	}
}


void texture_t::do_invert_y() {

	assert(is_allocated());
	unsigned const h2(height >> 1), wc(ncolors*width);
		
	for(unsigned i = 0; i < h2; ++i) {
		unsigned const off1(i*wc), off2((height-i-1)*wc);
			
		for(unsigned j = 0; j < wc; ++j) {
			swap(data[off1+j], data[off2+j]); // invert y
		}
	}
}


void texture_t::fix_word_alignment() {

	unsigned const byte_align = 4;
	if ((ncolors*width & (byte_align-1)) == 0) return; // nothing to do
	assert(is_allocated());
	float const ar(float(width)/float(height));
	int const new_w(width - (width&(byte_align-1)) + byte_align); // round up to next highest multiple of byte_align
	int const new_h(int(new_w/ar + 0.5)); // preserve aspect ratio
	resize(new_w, new_h);
}


void texture_t::add_alpha_channel() {

	if (ncolors == 4) return; // alread has an alpha channel
	assert(is_allocated());
	assert(ncolors == 1 || ncolors == 3); // only supported for RGB => RGBA for now
	bool const luminance(ncolors == 1);
	ncolors = 4; // add alpha channel
	unsigned char *new_data(new unsigned char[num_bytes()]);
	unsigned const npixels(num_pixels());

	for (unsigned i = 0; i < npixels; ++i) {
		UNROLL_3X(new_data[4*i+i_] = (luminance ? data[i] : data[3*i+i_]););
		new_data[4*i+3] = 255; // alpha of 255
	}
	free_data();
	data = new_data;
}


void texture_t::resize(int new_w, int new_h) {

	if (new_w == width && new_h == height) return; // already correct size
	assert(is_allocated());
	assert(width > 0 && height > 0 && new_w > 0 && new_h > 0);
	unsigned char *new_data(new unsigned char[new_w*new_h*ncolors]);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // required to handle calls from from fix_word_alignment()
	int const ret(gluScaleImage(calc_format(), width, height, get_data_format(), data, new_w, new_h, get_data_format(), new_data));
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	if (ret) cout << "GLU error during image scale: " << gluErrorString(ret) << "." << endl;
	free_data(); // only if size increases?
	data   = new_data;
	width  = new_w;
	height = new_h;
}


bool texture_t::try_compact_to_lum() {

	if (!CHECK_FOR_LUM || ncolors != 3) return 0;
	assert(is_allocated());
	// determine if it's really a luminance texture
	unsigned const npixels(num_pixels());
	bool is_lum(1);

	for (unsigned i = 0; i < npixels && is_lum; ++i) {
		is_lum &= (data[3*i+1] == data[3*i] && data[3*i+2] == data[3*i]);
	}
	if (!is_lum) return 0;
	//cout << "make luminance " << name << endl;
	// RGB equal, make it a single color (luminance) channel
	ncolors = 1; // add alpha channel
	unsigned char *new_data(new unsigned char[num_bytes()]);
	for (unsigned i = 0; i < npixels; ++i) {new_data[i] = data[3*i];}
	free_data();
	data = new_data;
	return 1;
}


void texture_t::make_normal_map() {

	if (ncolors == 3) return; // already a normal map
	
	if (ncolors == 4) { // better not have an alpha component
		cout << "Error: Skipping RGBA Bump/Normal map " << name << "." << endl;
		return;
	}
	assert(is_allocated());
	assert(ncolors == 1 && !is_16_bit_gray); // 8-bit grayscale heightmap
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
	free_data();
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
					unsigned const a_sum(a1 + a2 + a3 + a4);

					if (a_sum == 0) { // fully transparent - color is average of all 4 values
						UNROLL_3X(odata[ix1+i_] = (unsigned char)((idata[ix2+i_] + idata[ix2+xinc+i_] + idata[ix2+yinc+i_] + idata[ix2+yinc+xinc+i_]) / 4);)
						odata[ix1+3] = 0;
					}
					else {
						UNROLL_3X(odata[ix1+i_] = (unsigned char)((a1*idata[ix2+i_] + a2*idata[ix2+xinc+i_] + a3*idata[ix2+yinc+i_] + a4*idata[ix2+yinc+xinc+i_]) / a_sum);)
						odata[ix1+3] = min(255U, min(max(max(a1, a2), max(a3, a4)), unsigned(mipmap_alpha_weight*a_sum)));
					}
				}
			}
		}
		glTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, format, get_data_format(), &odata.front());
		idata.swap(odata);
	}
}


void texture_t::load_from_gl() { // also set tid?

	alloc();
	glGetTexImage(GL_TEXTURE_2D, 0, calc_format(), get_data_format(), data);
}


void bind_1d_texture(unsigned tid) {

	glBindTexture(GL_TEXTURE_1D, tid);
	assert(glIsTexture(tid)); // too slow?
}


void bind_2d_texture(unsigned tid) {

	glBindTexture(GL_TEXTURE_2D, tid);
	assert(glIsTexture(tid)); // too slow?
}


// 2D texture
void setup_texture(unsigned &tid, bool mipmap, bool wrap_s, bool wrap_t, bool mirror_s, bool mirror_t, bool nearest, float anisotropy) {

	assert(tid == 0);
	assert(!nearest || !mipmap);
	glGenTextures(1, &tid);

	// select our current texture
	bind_2d_texture(tid);

	// when texture area is small, use linear filter (bilinear filter the closest mipmap)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (nearest ? GL_NEAREST : (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR))); // GL_LINEAR_MIPMAP_NEAREST?

	// when texture area is large, bilinear filter the first mipmap
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (nearest ? GL_NEAREST : GL_LINEAR));

	// enable anisotropic filtering (slower but higher quality)
	if (anisotropy > 1.0) {glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, anisotropy);}

	// if wrap is true,  the texture wraps over at the edges (repeat) or is mirrored
	// if wrap is false, the texture ends at the edges (clamp)
	int const mode_s(wrap_s ? (mirror_s ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE); // Note: clamp is more efficient than wrap
	int const mode_t(wrap_t ? (mirror_t ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode_s);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode_t);
}


void setup_1d_texture(unsigned &tid, bool mipmap, bool wrap, bool mirror, bool nearest) {

	assert(tid == 0);
	assert(!nearest || !mipmap);
	glGenTextures(1, &tid);
	bind_1d_texture(tid);
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, (nearest ? GL_NEAREST : (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR))); // GL_LINEAR_MIPMAP_NEAREST?
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, (nearest ? GL_NEAREST : GL_LINEAR));
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, (wrap ? (mirror ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE));
}


void texture_t::gen_rand_texture(unsigned char val, unsigned char a_add, unsigned a_rand) {

	unsigned const size(num_pixels());
	assert(ncolors == 1 || ncolors == 4);

	for (unsigned i = 0; i < size; ++i) {
		if (ncolors == 1) { // alpha/luminance
			data[i] = (a_add + (unsigned char)(rand() % a_rand));
		}
		else {
			RGBA_BLOCK_ASSIGN((data+(i<<2)), val, val, val, (a_add + (unsigned char)(rand() % a_rand)));
		}
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
	textures[DISINT_TEX].gen_rand_texture(255, 1, 255);
}


void gen_tree_hemi_texture() {

	assert(SPHERE_SECTION >= 0.0 && SPHERE_SECTION <= 1.0);
	texture_t const &gtex(textures[HEDGE_TEX]);
	texture_t &tex(textures[TREE_HEMI_TEX]);
	assert(tex.width == gtex.width && tex.height == gtex.height);
	unsigned char const *grass_tex_data(gtex.get_data());
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


void update_player_bbb_texture(float extra_blood, bool recreate) {

	texture_t const &wood_tex(textures[WOOD_TEX]);
	texture_t &bbb_tex(textures[PLAYER_BBB_TEX]);
	assert(wood_tex.width == bbb_tex.width && wood_tex.height == bbb_tex.height && wood_tex.ncolors == bbb_tex.ncolors);
	unsigned char const *const wood_data(wood_tex.get_data());
	unsigned char *const bbb_data(bbb_tex.get_data());
	
	if (extra_blood == 0.0) { // copy/reset mode, just copy from WOOD_TEX
		unsigned const size(wood_tex.num_bytes());
		for (unsigned i = 0; i < size; ++i) {bbb_data[i] = wood_data[i];}
	}
	else { // add blood mode
		assert(bbb_tex.ncolors == 3);
		unsigned const npixels(wood_tex.num_pixels());
		unsigned const rval(max(1, round_fp(1.0/extra_blood)));

		for (unsigned i = 0; i < npixels; ++i) {
			if ((rand() % rval) == 0) {bbb_data[3*i+0] = 255; bbb_data[3*i+1] = 0; bbb_data[3*i+2] = 0;} // make all red
		}
	}
	if (recreate) {
		bbb_tex.gl_delete();
		bbb_tex.init(); // unnecessary?
	}
}


void gen_blur_inv_texture() {

	texture_t &tex(textures[BLUR_TEX_INV]);
	assert(tex.ncolors == 4);
	unsigned char *tex_data(tex.get_data());
	memset(tex_data, 255, 4*tex.num_pixels()*sizeof(char));

	for (int i = 0; i < tex.height; ++i) {
		unsigned char val(max(64, min(255, (2*255*(tex.height-i-1))/tex.height)));
		for (int j = 0; j < tex.width; ++j) {tex_data[((i*tex.width+j)<<2)+3] = val;}
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
	for (unsigned i = 0; i < size; ++i) {tex_data[i] = tex_data2[(i<<2)+3];} // put alpha in luminance
}


float wrap_tc(float v) { // unused
	float const fract(v - int(v));
	return CLIP_TO_01((v < 0.0) ? (1.0f + fract) : fract); // Note: clamped to allow for fp error
}


void noise_fill(unsigned char *data, unsigned size) {

	rand_gen_t rgen;
	for (unsigned i = 0; i < size; ++i) {data[i] = (rgen.rand() & 255);}
}


void gen_noise_texture() {

	for (unsigned i = NOISE_GEN_TEX; i <= NOISE_GEN_MIPMAP_TEX; ++i) {
		assert(textures[i].ncolors == 1);
		noise_fill(textures[i].get_data(), textures[i].num_pixels());
	}
	unsigned char *data(textures[SPARSE_NOISE_TEX].get_data());
	unsigned const size(textures[SPARSE_NOISE_TEX].num_pixels());
	rand_gen_t rgen;

	for (unsigned i = 0; i < size; ++i) {
		data[i] = (((rgen.rand()&7) == 0) ? 255 : 0);
		rgen.rand_mix();
	}
}


unsigned create_3d_noise_texture(unsigned size, unsigned ncomp) {

	vector<unsigned char> data(ncomp*size*size*size);
	noise_fill(&data.front(), data.size());
	return create_3d_texture(size, size, size, ncomp, data, GL_LINEAR, GL_REPEAT); // compressed?
}


unsigned get_noise_tex_3d(unsigned tsize, unsigned ncomp) {

	pair<texture_map_t::iterator, bool> ret(noise_tex_3ds.insert(make_pair(make_pair(tsize, ncomp), 0)));
	if (ret.second) {ret.first->second = create_3d_noise_texture(tsize, ncomp);}
	return ret.first->second;
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
	return ((relh > clip_hd1) ? ROCK_TEX : DIRT_TEX); // rock or dirt
}


void update_lttex_ix(int &ix) { // note: assumes lttex_dirt

	if ((water_is_lava || DISABLE_WATER == 2) && lttex_dirt[ix].id == SNOW_TEX) {--ix;}
	if (vegetation == 0.0  && lttex_dirt[ix].id == GROUND_TEX) {++ix;}
}


void get_tids(float relh, int NTEXm1, float const *const h_tex, int &k1, int &k2, float *t) {

	for (k1 = 0; k1 < NTEXm1 && relh >= h_tex[k1]; ++k1) {} // find first texture with height greater than relh

	if (k1 < NTEXm1 && (h_tex[k1] - relh) < TEXTURE_SMOOTH) {
		if (t) {*t = 1.0 - (h_tex[k1] - relh)/TEXTURE_SMOOTH;}
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
	if (using_custom_landscape_texture() || world_mode == WMODE_INF_TERRAIN) return;
	cached_ls_colors.clear();
	int tox0(0), toy0(0), scroll(0);
	texture_t &tex(textures[LANDSCAPE_TEX]);
	int const width(tex.width), height(tex.height);
	unsigned char *tex_data(tex.get_data());
	assert(tex.ncolors == 3);
	static int tox(0), toy(0);

	if (scrolling) { // ensure texture alignment when scrolling
		tox0 = (dx_scroll*width) /MESH_X_SIZE;
		toy0 = (dy_scroll*height)/MESH_Y_SIZE;
		tox += tox0;
		toy += toy0;
		int const x1(max(0, -tox0)), y1(max(0, -toy0)), x2(min(width, width-tox0)), y2(min(height, height-toy0));
		scroll = (x1 < x2 && y1 < y2 && abs(tox0) < width && abs(toy0) < height);
	}
	float const dz(zmax - zmin), dz_inv(1.0/dz);
	int const mxszm1(MESH_X_SIZE-1), myszm1(MESH_Y_SIZE-1), dxv(width/MESH_X_SIZE), dyv(height/MESH_Y_SIZE);
	int const NTEXm1(NTEX_DIRT-1), def_id((default_ground_tex >= 0) ? default_ground_tex : GROUND_TEX);
	float const xscale(((float)MESH_X_SIZE)/((float)width)), yscale(((float)MESH_Y_SIZE)/((float)height));
	static char **tids = NULL;
	if (tids == NULL) matrix_gen_2d(tids);
	assert(NTEX_DIRT < 128);
	
	for (int i = 0; i < MESH_Y_SIZE; ++i) { // makes a big performance improvement
		int const keepy(scroll && (i+dy_scroll) > 0 && (i+dy_scroll) < MESH_Y_SIZE-1);

		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (keepy && (j+dx_scroll) > 0 && (j+dx_scroll) < MESH_X_SIZE-1) continue;
			int const i1(min(myszm1, i+1)), j1(min(mxszm1, j+1));
			float const mh00(mesh_height[i][j]), mh01(mesh_height[i][j1]), mh10(mesh_height[i1][j]), mh11(mesh_height[i1][j1]);
			float const relh1(relh_adj_tex + (min(min(mh00, mh01), min(mh10, mh11)) - zmin)*dz_inv);
			float const relh2(relh_adj_tex + (max(max(mh00, mh01), max(mh10, mh11)) - zmin)*dz_inv);
			int k1a, k1b, k2a, k2b;
			get_tids(relh1, NTEXm1, h_dirt, k1a, k2a);
			get_tids(relh2, NTEXm1, h_dirt, k1b, k2b);
			tids[i][j] = char((k1a == k2b) ? k1a : -1);
		}
	}
	int const wx(width - dxv), hy(height - dyv);
	int const i0((toy0 < 0) ? height-1 : 0), i1((toy0 < 0) ? -1 : height), di((toy0 < 0) ? -1 : 1);
	int const j0((tox0 < 0) ? width -1 : 0), j1((tox0 < 0) ? -1 : width ), dj((tox0 < 0) ? -1 : 1);
	int const wxtx(wx-tox0), wxtx3(3*wxtx), j00(max(0, -tox0)), j01(min(width, wx-tox0));
	
	#pragma omp parallel for schedule(static,1)
	for (int ii = 0; ii < height; ++ii) {
		int const i((toy0 < 0) ? height-ii-1 : ii);
		float const yp(yscale*(float)i);
		int const lly(i + toy0), offset(3*i*width), ypos(max(0, min(myszm1, (int)yp))), ypos1(min(myszm1, ypos+1));
		float const ypi(yp - (float)ypos);
		int j10(j0);

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
					get_tids(relh, NTEXm1, h_dirt, k1, k2, &t);
					if (k1 != k2) assert(k2 == k1+1 || vegetation == 0.0);
				}
				else {
					k1 = k2 = tids[ypos][xpos];
				}
				id  = lttex_dirt[k1].id;
				id2 = lttex_dirt[k2].id;
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
			float const *const sti(sthresh[snow]);
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
	if (!using_custom_landscape_texture()) {
		if (landscape0 == NULL || !scroll) { // initialize/copy entire texture
			int const totsize(tex.num_bytes());
			if (landscape0 == NULL) landscape0 = new unsigned char[totsize];
			memcpy(landscape0, tex_data, totsize*sizeof(unsigned char));
			ls0_invalid = 0;
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
	PRINT_TIME(" Gen Landscape Texture");
}


void regrow_landscape_texture_amt0() {

	//RESET_TIME;
	static int counter(0);
	if (using_custom_landscape_texture() || iticks == 0) return; // is it too strong to use both iticks and fticks?
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

	if (using_custom_landscape_texture()) return 0.0;
	int const xpos(get_xpos(xval)), ypos(get_ypos(yval));
	if (point_outside_mesh(xpos, ypos))   return 0.0; // off the terrain area
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

	if (using_custom_landscape_texture()) return;
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


void add_color_to_landscape_texture(colorRGBA const &color, float xval, float yval, float radius) {

	if (using_custom_landscape_texture()) return;
	int const xpos0(get_xpos(xval)), ypos0(get_ypos(yval));
	if (point_outside_mesh(xpos0, ypos0)) return; // off the terrain area
	if ((display_mode & 0x04) && mesh_is_underwater(xpos0, ypos0)) return; // underwater + water enabled

	int const index(xpos0 + MESH_X_SIZE*ypos0);
	if (ls_color_texels.find(index) != ls_color_texels.end()) return; // don't update the same mesh location multiple times in the same frame
	ls_color_texels.insert(index);

	texture_t &tex(textures[LANDSCAPE_TEX]);
	unsigned char *tex_data(tex.get_data());
	float const xscale(((float)MESH_X_SIZE)/((float)tex.width)), yscale(((float)MESH_Y_SIZE)/((float)tex.height));
	int const maxsize(tex.num_bytes()), tsize(tex.num_pixels()); // used only for assertions
	int const xpos(int((xval + X_SCENE_SIZE)/(((float)xscale)*DX_VAL) + 0.5));
	int const ypos(int((yval + Y_SCENE_SIZE)/(((float)yscale)*DY_VAL) + 0.5));
	if (xpos < 0 || ypos < 0 || xpos >= tex.width || ypos >= tex.height) return;
	int const rad(max(0, int(5.0*radius/((xscale + yscale)*(DX_VAL + DY_VAL))) - 1)), rad_sq(max(1, rad*rad)); // size in texture space
	int const x1(max(0, xpos-rad)), y1(max(0, ypos-rad)), x2(min(tex.width-1, xpos+rad)), y2(min(tex.height-1, ypos+rad));
	float const blend_scale(1.0/float(rad_sq));
	unsigned char color_i[3];
	unpack_color(color_i, color);

	for (int i = y1; i <= y2; ++i) {
		int const offset(i*tex.width), yterm((i - ypos)*(i - ypos));

		for (int j = x1; j <= x2; ++j) {
			int const dist_sq((yterm + (j - xpos)*(j - xpos)));

			if (dist_sq < rad_sq) {
				assert(offset + j <= tsize);
				float const blend(0.8*color.alpha*(1.0 - float(dist_sq)*blend_scale));
				int const o2(3*(offset + j));
				assert(o2 + 3 <= maxsize);
				BLEND_COLOR((tex_data + o2), color_i, (tex_data + o2), blend); // test if this texel is enabled for draw?
			}
		}
	}
	update_lt_section(x1, y1, x2+1, y2+1); // can be slow
}


void add_snow_to_landscape_texture(point const &pos, float acc) {

	if (acc <= 0.0 || using_custom_landscape_texture()) return;
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
	assert(!using_custom_landscape_texture());
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
	return (relh > h_dirt[3]);
}


void gen_tex_height_tables() {

	for (unsigned i = 0; i < NTEX_DIRT; ++i) {h_dirt[i] = pow(lttex_dirt[i].zval, glaciate_exp);}
	clip_hd1 = (0.90*h_dirt[1] + 0.10*h_dirt[0]);
}


void set_texgen_vec4(vector4d const &v, bool s_or_t, shader_t &shader, int mode) {

	switch (mode) { // mode: 0 = shader uniform, 1 = shader attrib, 2 = shader uniform #2
	case 0: shader.add_uniform_vector4d((s_or_t ? "texgen_t"  : "texgen_s" ), v); break;
	case 2: shader.add_uniform_vector4d((s_or_t ? "texgen2_t" : "texgen2_s"), v); break;
	case 1: shader.add_attrib_float_array((s_or_t ? TEX0_T_ATTR : TEX0_S_ATTR), &v.x, 4); break;
	default: assert(0);
	}
}


void setup_texgen_full(float sx, float sy, float sz, float sw, float tx, float ty, float tz, float tw, shader_t &shader, int mode) {

	set_texgen_vec4(vector4d(sx, sy, sz, sw), 0, shader, mode);
	set_texgen_vec4(vector4d(tx, ty, tz, tw), 1, shader, mode);
}


void setup_texgen(float xscale, float yscale, float tx, float ty, float z_off, shader_t &shader, int mode) {

	assert(xscale != 0.0 && yscale != 0.0);
	setup_texgen_full(xscale, 0.0, z_off, tx, 0.0, yscale, z_off, ty, shader, mode);
}


void get_poly_texgen_dirs(vector3d const &norm, vector3d v[2]) { // similar to get_ortho_vectors()

	v[0] = all_zeros;
	v[0][get_min_dim(norm)] = 1.0;
	cross_product(norm, v[0], v[1]);
	cross_product(norm, v[1], v[0]);
}


void setup_polygon_texgen(vector3d const &norm, float const scale[2], float const xlate[2],
	vector3d const &offset, bool swap_txy, shader_t &shader, int mode)
{
	vector3d v[2];
	get_poly_texgen_dirs(norm, v);
	
	for (unsigned i = 0; i < 2; ++i) {
		vector4d const v(scale[i]*v[i].x, scale[i]*v[i].y, scale[i]*v[i].z, (xlate[i] + scale[i]*dot_product(offset, v[i])));
		set_texgen_vec4(v, ((i != 0) ^ swap_txy), shader, mode);
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


unsigned texture_t::get_texel_ix(float u, float v) const {

	u -= 1.0E-6f; v -= 1.0E-6f; // make 1.0 equal width-1 | height-1
	//int tx(int(width *u) & (width -1));
	//int ty(int(height*v) & (height-1)); // assumes width and height are a power of 2
	int tx(int(width *u) % width ); // width and height can be any nonzero value
	int ty(int(height*v) % height);
	if (tx < 0) tx += width;
	if (ty < 0) ty += height;
	assert(tx >= 0 && ty >= 0 && tx < width && ty < height);
	return (width*ty + tx);
}


colorRGBA texture_t::get_texel(unsigned ix) const {

	assert(ix < num_pixels());
	unsigned const off(ncolors*ix);

	switch (ncolors) {
	case 1: // grayscale
		{float const val(data[off]/255.0); return colorRGBA(val, val, val, 1.0);} // alpha=1.0
	case 2: // luminance + alpha
		{float const val(data[off]/255.0); return colorRGBA(val, val, val, data[off+1]/255.0);}
	case 3: // RGB
		return colorRGBA(data[off+0]/255.0, data[off+1]/255.0, data[off+2]/255.0, 1.0); // alpha=1.0
	case 4: // RGBA
		return colorRGBA(data[off+0]/255.0, data[off+1]/255.0, data[off+2]/255.0, data[off+3]/255.0);
	default: assert(0);
	}
	return BLACK; // never gets here
}


float texture_t::get_component(float u, float v, int comp) const {

	assert(comp < ncolors);
	return data[ncolors*get_texel_ix(u, v) + comp]/255.0;
}


float get_texture_component(unsigned tid, float u, float v, int comp) {

	assert(tid < NUM_TEXTURES);
	return textures[tid].get_component(u, v, comp);
}

colorRGBA get_texture_color(unsigned tid, float u, float v) {

	assert(tid < NUM_TEXTURES);
	return textures[tid].get_texel(u, v);
}


vector2d get_billboard_texture_uv(point const *const points, point const &pos) {

	assert(points != NULL);
	// ordered: (s,t) => (0,1), (0,0), (1,0), (1,1)
	float d[4]; // distance from coll point to quad edge

	for (unsigned i = 0; i < 4; ++i) {
		unsigned const in((i+1)&3);
		d[i] = cross_product((pos - points[i]), (pos - points[in])).mag()/p2p_dist(points[i], points[in]);
	}
	assert(d[0] + d[2] > 0.0);
	assert(d[1] + d[3] > 0.0);
	vector2d uv(d[0]/(d[0] + d[2]), d[1]/(d[1] + d[3])); // y is upside down
	uv.x = CLIP_TO_01(uv.x); uv.y = CLIP_TO_01(uv.y); // clamp uv to account for fp rounding errors (and incorrect pos?)
	//assert(uv.x >= 0.0 && uv.x <= 1.0 && uv.y >= 0.0 && uv.y <= 1.0); // more restrictive case
	return uv;
}


bool is_billboard_texture_transparent(point const *const points, point const &pos, int tid) {

	vector2d const uv(get_billboard_texture_uv(points, pos));
	return (get_texture_component(tid, uv.x, uv.y, 3) == 0.0);
}


void ensure_texture_loaded(unsigned &tid, unsigned txsize, unsigned tysize, bool mipmap, bool nearest) { // used with texture_pair_t/render-to-texture RGBA

	assert(txsize > 0 && tysize > 0);
	if (tid) return; // already created
	setup_texture(tid, mipmap, 0, 0, 0, 0, nearest);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, txsize, tysize, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
}


void build_texture_mipmaps(unsigned tid, unsigned dim) {

	bind_2d_texture(tid);
	gen_mipmaps(dim);
}


void texture_pair_t::free_context() {
	free_texture(tids[0]);
	free_texture(tids[1]);
}

void texture_pair_t::bind_texture() const {

	assert(is_valid());
	bind_2d_texture(tids[0]);
	set_active_texture(1);
	bind_2d_texture(tids[1]);
	set_active_texture(0);
}

void texture_pair_t::ensure_tid(unsigned tsize, bool mipmap) {
	ensure_texture_loaded(tids[0], tsize, tsize, mipmap, 0); // color
	ensure_texture_loaded(tids[1], tsize, tsize, mipmap, 0); // normal
}


void texture_atlas_t::free_context() {free_texture(tid);}

void texture_atlas_t::bind_texture() const {
	assert(tid);
	bind_2d_texture(tid);
}

void texture_atlas_t::ensure_tid(unsigned base_tsize, bool mipmap) {
	assert(nx > 0 && ny > 0);
	ensure_texture_loaded(tid, nx*base_tsize, ny*base_tsize, mipmap, 0);
}


