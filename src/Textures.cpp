// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/25/02
#include "3DWorld.h"
#include "mesh.h"
#include "sinf.h"
#include "textures.h"
#include "gl_ext_arb.h"
#include "shaders.h"


float const TEXTURE_SMOOTH        = 0.01;
float const SPHERE_SECTION        = 0.75;
float const TEXTURE_SNOW_ACC_RATE = 0.05;
int   const SNOW_ACC_RADIUS       = 3;
float const LANDSCAPE_REGEN_AMT   = 0.02; // how much regen after LANDSCAPE_REGEN_MOD ticks
int   const LANDSCAPE_REGEN_MOD   = 32;   // ticks for entire mesh regen pass
int   const MAX_UPDATE_SIZE       = 16;
float const LS_TEX_ANISO          = 2.0; // value of anisotropic texture filtering for landscape source textures

bool const RELOAD_TEX_ON_HOLE  = 0;
bool const LANDSCAPE_MIPMAP    = 1; // looks better, but texture update requires mipmap updates
bool const SHOW_TEXTURE_MEMORY = 0;
bool const COMPRESS_TEXTURES   = 1;
bool const ALLOW_SLOW_COMPRESS = 1;
bool const USE_STB_DXT         = 1;
bool const CHECK_FOR_LUM       = 1;
float const SMOOTH_SKY_POLES   = 0.2;


texture_t def_textures[NUM_PREDEF_TEXTURES] = { // 4 colors without wrap sometimes has a bad transparent strip on spheres
// type: 0 = read from file, 1 = generated, 2 generated and dynamically updated
// format: 0: RGB RAW, 1: BMP, 2: RGB RAW, 3: RGBA RAW, 4: targa (*tga), 5: jpeg, 6: png, 7: auto, 8: tiff, 9: generated (not loaded from file), 10: DDS, 11:ppm
// use_mipmaps: 0 = none, 1 = standard OpenGL, 2 = openGL + CPU data, 3 = custom alpha OpenGL, 4 = custom alpha OpenGL using average texture color for transparent pixels
// wrap_mir: 0 = clamp, 1 = wrap, 2 = mirror
// type format width height wrap_mir ncolors use_mipmaps name [invert_y=0 [do_compress=1 [anisotropy=1.0 [mipmap_alpha_weight=1.0 [normal_map=0]]]]]
//texture_t(0, 6, 512,  512,  1, 3, 0, "ground.png"),
texture_t(0, 6, 128,  128,  1, 3, 2, "grass.png", 0, 1, LS_TEX_ANISO), // mipmap for small trees?
//texture_t(0, 5, 0,    0,    1, 3, 2, "grass_new.jpg", 0, 1, LS_TEX_ANISO), // 1024x1024; has texture seams, not as bright as other grass
texture_t(0, 6, 256,  256,  1, 3, 1, "rock.png"),
texture_t(0, 5, 512,  512,  1, 3, 1, "water.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "stucco.jpg", 0, ALLOW_SLOW_COMPRESS),
texture_t(0, 5, 0,    0,    1, 4, 0, "sky.jpg", 1), // 1024x1024
texture_t(0, 5, 0,    0,    1, 3, 1, "brick1.jpg"), // brick2?
texture_t(0, 5, 0,    0,    1, 3, 1, "moon.jpg"),
texture_t(0, 6, 256,  256,  0, 3, 1, "earth.png", 1),
texture_t(0, 5, 0,    0,    1, 3, 1, "marble.jpg", 0, ALLOW_SLOW_COMPRESS), // or marble2.jpg
texture_t(0, 7, 0,    0,    1, 3, 2, "snow2.jpg", 0, 1, LS_TEX_ANISO),
texture_t(0, 5, 0,    0,    0, 4, 4, "leaves/green_maple_leaf.jpg", 1, 1, 4.0), // 960x744
//texture_t(0, 6, 0,    0,    0, 4, 4, "leaves/maple_leaf.png", 1, 1, 4.0), // 344x410
texture_t(0, 5, 0,    0,    1, 3, 1, "bark2.jpg"), // 512x512 (Note: must match baseball bat texture size)
texture_t(0, 5, 512,  512,  1, 3, 2, "desert_sand.jpg", 0, 1, LS_TEX_ANISO),
texture_t(0, 6, 256,  256,  1, 3, 2, "rock2.png", 0, 1, LS_TEX_ANISO),
texture_t(0, 5, 512,  512,  1, 3, 1, "camoflage.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "hedges.jpg", 0, ALLOW_SLOW_COMPRESS), // 1024x1024
texture_t(0, 1, 512,  512,  1, 3, 1, "brick1.bmp", 0, 1, 8.0),
texture_t(0, 5, 512,  512,  1, 3, 1, "manhole.jpg", 1),
texture_t(0, 5, 0,    0,    0, 4, 4, "leaves/palm_frond_diff.jpg", 0, 1, 4.0), // 512x1024
texture_t(1, 9, 256,  256,  1, 4, 1, "@smoke"),  // not real file
texture_t(1, 9, 64,   64,   1, 4, 1, "@plasma"), // not real file
texture_t(1, 9, 8,    8,    0, 3, 0, "@gen"),    // not real file - unused placeholder
texture_t(2, 7, 1024, 1024, 0, 3, LANDSCAPE_MIPMAP, "@landscape_tex"), // for loading real landscape texture
texture_t(1, 9, 128,  128,  0, 3, 0, "@tree_end"),  // not real file
texture_t(1, 9, 1024, 1024, 1, 4, 1, "@tree_hemi", 0, 1), // not real file, mipmap for trees?
texture_t(0, 5, 0  ,  0,    1, 3, 1, "shingles.jpg", 0, ALLOW_SLOW_COMPRESS, 8.0),
texture_t(0, 6, 256,  256,  1, 3, 1, "paneling.png", 0, 1, 16.0),
texture_t(0, 6, 256,  256,  1, 3, 1, "cblock.png", 0, 1, 8.0),
texture_t(0, 5, 0,    0,    0, 4, 3, "mj_leaf.jpg", 1), // 128x128
texture_t(0, 6, 0,    0,    0, 4, 4, "leaves/oak_leaf.png", 1, 1, 4.0), // 208x348
texture_t(0, 6, 0,    0,    0, 4, 4, "leaves/cherry_leaf.png", 1, 1, 4.0), // 576x1220
texture_t(0, 6, 0,    0,    0, 4, 4, "leaves/birch_leaf.png", 1, 1, 4.0), // 838x1372
texture_t(0, 5, 0,    0,    0, 4, 3, "plant1.jpg", 1), // 256x256
texture_t(0, 6, 256,  256,  0, 4, 3, "plant2.png", 1),
//texture_t(0, 5, 0,    0,    0, 4, 3, "plant2.jpg", 1), // 160x256
texture_t(0, 6, 256,  256,  0, 4, 3, "plant3.png", 1),
//texture_t(0, 5, 0,    0,    0, 4, 3, "plant3.jpg", 1), // 176x256
texture_t(0, 5, 0,    0,    0, 4, 4, "leaves/leaf_d.jpg", 1), // 200x500
texture_t(0, 5, 0,    0,    1, 3, 1, "fence.jpg", 0, ALLOW_SLOW_COMPRESS, 8.0), // 896x896
texture_t(0, 6, 128,  128,  1, 3, 1, "skull.png"),
texture_t(0, 6, 64,   64,   1, 3, 1, "radiation.png", 1),
texture_t(0, 6, 128,  128,  1, 3, 1, "yuck.png"),
texture_t(0, 6, 256,  256,  0, 4, 0, "sawblade.png", 1),
texture_t(0, 6, 256,  256,  0, 4, 0, "sawblade_b.png", 1),
texture_t(0, 6, 256,  256,  0, 4, 1, "blur.png", 0, 0), // disable compression - causes artifacts when drawing star flares
texture_t(0, 6, 256,  256,  1, 4, 1, "blur_s.png"),
texture_t(0, 5, 0,    0,    0, 4, 3, "pine2.jpg", 1, 1, 1.0, 0.5),
texture_t(0, 6, 128,  128,  1, 3, 1, "noise.png"),
texture_t(0, 5, 0,    0,    1, 3, 1, "wood.jpg", 0, ALLOW_SLOW_COMPRESS, 4.0), // 768x768
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
texture_t(0, 6, 0  ,  0  ,  0, 4, 1, "atlas/explosion.png", 1, 0), // no compression
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
texture_t(0, 6, 0  ,  0  ,  0, 4, 0, "atlas/fire.png", 0, 0), // no compression
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
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark2-normal.jpg", 0, 0, 4.0, 1.0, 1), // 512x512, no compress
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark_lendrick.jpg", 0, ALLOW_SLOW_COMPRESS), // 892x892
texture_t(0, 6, 0,    0,    1, 3, 1, "bark/bark_lylejk.png", 0, ALLOW_SLOW_COMPRESS), // 1024x768
// normal/caustic maps
texture_t(0, 4, 0,    0,    1, 3, 1, "normal_maps/water_normal.tga", 0, 0, 8.0, 1.0, 1), // 512x512, no compress
texture_t(0, 6, 0,    0,    1, 3, 1, "normal_maps/ocean_water_normal.png", 0, 0, 4.0, 1.0, 1), // 1024x1024 (Note: compression disabled as it causes artifacts)
texture_t(0, 5, 0,    0,    1, 3, 1, "caustics.jpg"), // 512x512
// noise
texture_t(0, 6, 0,    0,    1, 1, 0, "perlin_simplex.png"), // 256x256
texture_t(1, 9, 128,  128,  1, 1, 0, "@noise_gen"), // not real file
texture_t(1, 9, 128,  128,  1, 1, 1, "@noise_gen_mipmap"), // not real file
texture_t(1, 9, 256,  256,  1, 1, 1, "@noise_gen_sparse"), // not real file
texture_t(1, 9, 400,  400,  1, 3, 1, "@player_bbb_tex"), // not real file (Note: must match bark texture size)
texture_t(0, 5, 0,    0,    0, 4, 3, "pine_tree_leaves2.jpg", 1, 0, 1.0, 0.5), // 512x512
texture_t(0, 5, 0,    0,    0, 4, 1, "flare1.jpg"), // 384x384
texture_t(0, 5, 0,    0,    0, 4, 1, "flare2.jpg"), // 128x128 (Note: low resolution)
texture_t(0, 5, 0,    0,    0, 4, 1, "Flare3.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 4, 1, "flare4.jpg"), // 256x256
//texture_t(0, 5, 0,    0,    0, 4, 1, "flare5.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 4, 1, "flare5b.jpg"), // 416x416
texture_t(0, 5, 0,    0,    1, 3, 1, "foam1.jpg"), // 512x512

texture_t(0, 5, 0,    0,    0, 3, 0, "bullet_hole/bullet_diffuse.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 1, 0, "bullet_hole/bullet_alpha.jpg"), // 256x256
texture_t(0, 5, 0,    0,    0, 3, 0, "bullet_hole/bullet_normal.jpg", 0, 0, 1.0, 1.0, 1), // 256x256, no compress
texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/sand_normal.jpg", 0, 0, 2.0, 1.0, 1), // (Note: compression disabled as it causes artifacts)
//texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/sand_dunes.jpg", 0, 0, 2.0, 1.0, 1), // okay for sand, but not for other mesh textures/materials
//texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/test_normal.jpg", 0, 0, 2.0, 1.0, 1),

texture_t(0, 5, 0,    0,    1, 3, 1, "raindrop_dots.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "spaceship1.jpg"),
texture_t(0, 5, 0,    0,    1, 3, 1, "spaceship2.jpg"),
texture_t(0, 6, 0,    0,    0, 4, 1, "atlas/blood.png"),
texture_t(0, 5, 0,    0,    1, 3, 1, "lichen.jpg", 0, ALLOW_SLOW_COMPRESS), // 1500x1500
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/palm_bark.jpg"), // 512x512
texture_t(0, 5, 0,    0,    0, 4, 0, "daisy.jpg", 0, 1, 4.0), // 1024x1024 - no mipmap to avoid filtering artifacts making distant flowers look square (but not too bad with mode 3)
texture_t(0, 5, 0,    0,    1, 3, 1, "lava.jpg"), // 512x512
texture_t(0, 5, 0,    0,    0, 4, 1, "smoke_puff.jpg"), // 150x150
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark_birch.jpg", 0, 0), // 512x512, no compress
texture_t(0, 5, 0,    0,    1, 3, 1, "bark/bark6.jpg", 0, 0), // 894x894, no compress
texture_t(0, 6, 0,    0,    1, 4, 1, "ripple_map.png", 0, 0), // 256x256, mipmaps?, no compress?
texture_t(0, 6, 0,    0,    1, 4, 1, "starburst.png",  0, 0), // disable compression - causes artifacts
texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/rocks1.jpg", 0, 0, 2.0, 1.0, 1), // sand_normal = snow, rocks1 = rock, rocks2 = sand, rocks3 = dirt/grass
texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/rocks2.jpg", 0, 0, 2.0, 1.0, 1),
texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/rocks3.jpg", 0, 0, 2.0, 1.0, 1),
texture_t(0, 5, 0,    0,    1, 3, 1, "normal_maps/dirt_normal.jpg", 0, 0, 2.0, 1.0, 1),
texture_t(0, 6, 16,   16,   1, 3, 0, "cyan.png"), // for normal maps
texture_t(0, 6, 16,   16,   1, 3, 0, "red.png"), // for TT default sand weights texture
texture_t(0, 5, 0,    0,    1, 3, 1, "hazard_stripes.jpg", 0, 0, 4.0), // 500x500
texture_t(1, 9, 512,  512,  1, 4, 1, "@windows" , 0, 0, 4.0),  // not real file
texture_t(1, 9, 512,  512,  1, 4, 1, "@twindows", 0, 0, 4.0),  // not real file
texture_t(0, 6, 0,    0,    1, 3, 1, "keycard.png", 0, 1, 4.0), // 512x512
// type format width height wrap_mir ncolors use_mipmaps name [invert_y=0 [do_compress=1 [anisotropy=1.0 [mipmap_alpha_weight=1.0 [normal_map=0]]]]]
};

vector<texture_t> textures;


// zval should depend on def_water_level and temperature
int max_tius(0), max_ctius(0); // cached in case they are needed somewhere (for shadow map logic, etc.)
float h_dirt[NTEX_DIRT], clip_hd1;
std::set<int> ls_color_texels;
vector<colorRGBA> cached_ls_colors;

typedef map<pair<unsigned, unsigned>, unsigned> texture_map_t;
texture_map_t noise_tex_3ds;

typedef map<string, unsigned> name_map_t;
name_map_t texture_name_map;

bool textures_inited(0), def_tex_compress(1);
int landscape_changed(0), lchanged0(0), skip_regrow(0), ltx1(0), lty1(0), ltx2(0), lty2(0), ls0_invalid(1);
unsigned sky_zval_tid;
float def_tex_aniso(2.0);
unsigned char *landscape0 = NULL;


extern bool mesh_difuse_tex_comp, water_is_lava, invert_bump_maps, no_store_model_textures_in_memory;
extern unsigned smoke_tid, dl_tid, elem_tid, gb_tid, dl_bc_tid, reflection_tid, room_mirror_ref_tid, depth_tid, empty_smap_tid;
extern unsigned frame_buffer_RGB_tid, skybox_tid, skybox_cube_tid, univ_reflection_tid;
extern int world_mode, read_landscape, default_ground_tex, xoff2, yoff2, DISABLE_WATER;
extern int scrolling, dx_scroll, dy_scroll, display_mode, iticks, universe_only, window_width, window_height;
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
unsigned get_building_textures_gpu_mem();


bool is_tex_disabled(int i) {
	return (i == GEN_TEX || (universe_only && (i == CLOUD_RAW_TEX || i == WIND_TEX || i == LANDSCAPE_TEX || i == TREE_END_TEX || i == TREE_HEMI_TEX)));
}


void load_texture_names() {

	if (!texture_name_map.empty()) return; // already loaded
	textures.resize(sizeof(def_textures)/sizeof(texture_t));

	for (unsigned i = 0; i < textures.size(); ++i) {
		textures[i] = def_textures[i];
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

	timer_t timer("Texture Load");
	cout << "Loading " << textures.size() << " textures" << endl;
	if (using_custom_landscape_texture()) {set_landscape_texture_from_file();} // must be done first
	load_texture_names();

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < (int)textures.size(); ++i) {
		//cout << "."; cout.flush();
		if (!is_tex_disabled(i)) {textures[i].load(i, 0, 0, 1);} // ignore word alignment here, since resizing isn't thread safe
	}
	for (int i = 0; i < (int)textures.size(); ++i) {
		if (!is_tex_disabled(i)) {textures[i].fix_word_alignment();}
	}
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
	for (unsigned i = 0; i < textures.size(); ++i) {
		if (is_tex_disabled(i)) continue; // skip
		if (i == BLDG_WINDOW_TEX || i == BLDG_WIND_TRANS_TEX || i == LANDSCAPE_TEX) continue; // not yet generated
		textures[i].init();
	}
	textures[TREE_HEMI_TEX].set_color_alpha_to_one();
	textures_inited = 1;

	glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &max_tius);
	glGetIntegerv(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS, &max_ctius);
	cout << "max TIUs: " << max_tius << ", max combined TIUs: " << max_ctius << endl;
}


unsigned get_loaded_textures_cpu_mem() {
	unsigned mem(0);
	for (texture_t const &t : textures) {mem += t.get_cpu_mem();}
	return mem;
}
unsigned get_loaded_textures_gpu_mem() {
	unsigned mem(get_building_textures_gpu_mem());
	for (texture_t const &t : textures) {mem += t.get_gpu_mem();}
	return mem;
}
void print_texture_memory_usage() { // full cities scene: 142MB / 272MB (660MB uncompressed) (244MB with ALLOW_SLOW_COMPRESS=1)
	cout << "Texture Memory: " << get_loaded_textures_cpu_mem() << " CPU / " << get_loaded_textures_gpu_mem() << " GPU" << endl;
}
void print_texture_stats() {
	cout << "Textures: " << textures.size() << endl;
	print_texture_memory_usage();
}


int texture_lookup(string const &name) {
	name_map_t::const_iterator it(texture_name_map.find(name));
	return ((it != texture_name_map.end()) ? it->second : -1);
}

int get_texture_by_name(string const &name, bool is_normal_map, bool invert_y, int wrap_mir, float aniso,
	bool allow_compress, int use_mipmaps, unsigned ncolors, bool is_alpha_mask)
{
	if (name.empty()) return -1; // no texture
	int const ix(atoi(name.c_str()));
	if (ix > 0 || ix == -1 || name == "0") return ix; // a number was specified
	if (name == "none" || name == "null")  return -1; // no texture
	int tid(texture_lookup(name));
	if (tid >= 0) {assert((unsigned)tid < textures.size()); return tid;}
#ifndef ENABLE_PNG
	// no test; handled by stb_image as a fallback
#endif
#ifndef ENABLE_TIFF
	if (endswith(name, ".tif") || endswith(name, ".tiff")) return -1; // tiff format not enabled
#endif
#ifndef ENABLE_DDS
	if (endswith(name, ".dds")) return -1; // dds format not enabled
#endif
	//timer_t timer("Load Texture " + name);
	// try to load/add the texture directly from a file: assume it's RGB with wrap and mipmaps
	assert(omp_get_thread_num_3dw() == 0); // must be serial
	tid = textures.size();
	bool const do_compress(allow_compress && def_tex_compress && !is_normal_map);
	// type format width height wrap_mir ncolors use_mipmaps name [invert_y=0 [do_compress=1 [anisotropy=1.0 [mipmap_alpha_weight=1.0 [normal_map=0]]]]]
	texture_t new_tex(0, IMG_FMT_AUTO, 0, 0, wrap_mir, ncolors, use_mipmaps, name, invert_y, do_compress,
		((aniso > 0.0) ? aniso : def_tex_aniso), 1.0, is_normal_map);

	if (textures_inited) {
		new_tex.load(tid);
		if ((is_alpha_mask || ncolors == 1) && new_tex.ncolors == 4) {new_tex.fill_to_grayscale_color(255);} // alpha mask - fill color to white
		new_tex.init();
	}
	textures.push_back(new_tex);
	texture_name_map[name] = tid;
	return tid; // Note: city/buildings scene currently loads 348 textures
}

unsigned load_cube_map_texture(string const &name) { // used for sky boxes

	// Rather than specify each of the 6 textures (one per side), we're going to specify the top texture and try to guess the names of the others
	// There doesn't seem to be any well defined specs for which face corresponds to which direction, other than that Y is up (vs. Z being up in 3DWorld),
	// so we'll just hard-code the direction to face mapping that seems to look correct for each of our test image sets
	size_t const dot_pos(name.find_last_of("."));
	if (dot_pos == string::npos) {std::cerr << "Failed to find file extension in cube map texture name" << endl; return 0;}
	string const ext(name.substr(dot_pos)), prefix(name.substr(0, (name.size() - ext.size())));
	string names[6];

	if (endswith(prefix, "top")) {
		string const base(prefix.substr(0, (prefix.size() - 3)));
		string const dirs[6] = {"right", "left", "bottom", "top", "front", "back"};
		for (unsigned n = 0; n < 6; ++n) {names[n] = base + dirs[n] + ext;}
	}
	else if (endswith(prefix, "up")) {
		string const base(prefix.substr(0, (prefix.size() - 2)));
		string const dirs[6] = {"ft", "bk", "dn", "up", "rt", "lf"};
		for (unsigned n = 0; n < 6; ++n) {names[n] = base + dirs[n] + ext;}
	}
	else {
		std::cerr << "Expecting cube map texture filename to end with 'top' or 'up': " << prefix << endl;
		return 0;
	}
	for (unsigned n = 0; n < 6; ++n) {
		if (check_texture_file_exists(names[n])) continue;
		std::cerr << "Error: Failed to load cube map texture '" << names[n] << "'; Skipping file '" << name << "'" << endl;
		return 0; // fails if any of the cube side textures don't exist
	}
	unsigned tid(0), tex_size(0);
	bool const allocate(0), use_mipmaps(0), do_compress(1);
	setup_cube_map_texture(tid, tex_size, allocate, use_mipmaps, 1.0);
	texture_t texture(0, IMG_FMT_AUTO, 0, 0, 0, 3, use_mipmaps, "skybox", 0, do_compress); // type format width height wrap_mir ncolors use_mipmaps name [invert_y=0 [do_compress=1]]

	for (unsigned n = 0; n < 6; ++n) {
		texture.name = names[n];
		texture.load(0);
		texture.upload_cube_map_face(n);
		texture.free_data();
	}
	if (use_mipmaps) {gen_mipmaps(6);}
	return tid;
}


void check_init_texture(int id, bool free_after_upload) {textures[id].check_init(free_after_upload);}

bool select_texture(int id, unsigned tu_id) {
	bool const no_tex(id < 0);
	if (no_tex) {id = WHITE_TEX;} //glBindTexture(GL_TEXTURE_2D, 0); // bind to none
	assert((unsigned)id < textures.size());
	bool const free_after_upload(no_store_model_textures_in_memory && id >= NUM_PREDEF_TEXTURES); // free textures loaded by name only
	check_init_texture(id, free_after_upload);
	textures[id].bind_gl(tu_id);
	return !no_tex;
}

void bind_texture_tu_def_white_tex(unsigned tid, unsigned tu_id) {
	if (tid == 0) {select_texture(WHITE_TEX, tu_id);}
	else {bind_texture_tu(tid, tu_id);}
}

float get_tex_ar(int id) {
	if (id < 0) return 1.0;
	texture_t const &texture(get_texture_by_id(id));
	return (((double)texture.width)/((double)texture.height));
}

void free_textures() {
	for (unsigned i = 0; i < textures.size(); ++i) {textures[i].gl_delete();}
}


void reset_textures() {

	cout << "Freeing textures..." << endl; // only print if something was loaded?
	free_textures();
	free_smiley_textures(); // should this be guarded by a conditional?
	free_flare_textures();
	free_shadow_map_textures();
	free_cloud_textures();
	free_teleporter_textures();
	free_texture(smoke_tid);
	free_texture(dl_tid);
	free_texture(elem_tid);
	free_texture(gb_tid);
	free_texture(dl_bc_tid);
	free_texture(reflection_tid);
	free_texture(depth_tid);
	free_texture(sky_zval_tid);
	free_texture(empty_smap_tid);
	free_texture(frame_buffer_RGB_tid);
	free_texture(skybox_tid);
	free_texture(skybox_cube_tid);
	free_texture(univ_reflection_tid);
	free_texture(room_mirror_ref_tid);
	free_building_indir_texture();
	free_font_texture_atlas();
	for (texture_map_t::iterator i = noise_tex_3ds.begin(); i != noise_tex_3ds.end(); ++i) {free_texture(i->second);}
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

void texture_t::bind_gl(unsigned tu_id) const {
	assert(tid > 0);
	bind_texture_tu(tid, tu_id);
}
GLuint64 texture_t::get_bindless_handle(bool make_resident) const {
	assert(tid > 0);
	GLuint64 const handle(glGetTextureHandleARB(tid));
	if (make_resident) {glMakeTextureHandleResidentARB(handle);}
	return handle;
}
void texture_t::gl_delete() {
	free_texture(tid);
}

void texture_t::free_client_mem() {
	if (orig_data    != data) {delete [] orig_data;}
	if (colored_data != data) {delete [] colored_data;}
	delete [] data;
	data = orig_data = colored_data = NULL;
}

bool texture_t::is_texture_compressed() const {
	return (COMPRESS_TEXTURES && do_compress && type != 2); // no compress if dynamically updated
}
GLenum texture_t::calc_internal_format() const {
	assert(ncolors >= 1 && ncolors <= 4);
	if (is_16_bit_gray) {return GL_R16;} // compressed?
	return get_internal_texture_format(ncolors, is_texture_compressed(), 0); // linear_space=0
}
GLenum texture_t::calc_format() const {
	return (is_16_bit_gray ? GL_RED : get_texture_format(ncolors));
}


void texture_t::do_gl_init(bool free_after_upload) {

	//cout << "bind texture " << name << " size " << width << "x" << height << endl;
	//timer_t timer(("Load and Upload Texture " + name), 1, 1);
	setup_texture(tid, (use_mipmaps != 0 && !defer_load()), wrap, wrap, mirror, mirror, 0, anisotropy);

	if (SHOW_TEXTURE_MEMORY) {
		static unsigned tmem(0);
		tmem += get_gpu_mem();
		cout << "tex vmem = " << tmem << endl;
	}
	if (defer_load()) {deferred_load_and_bind();} // Note: mipmaps are stored in the DDS file and aren't controlled by the use_mipmaps option
	else {
		assert(is_allocated());
		assert(width > 0 && height > 0);
		bool const compressed(is_texture_compressed()), use_custom_compress(USE_STB_DXT && compressed && (ncolors == 3 || ncolors == 4));

		if (use_custom_compress) {compress_and_send_texture();} // compressed RGB or RGBA
		else { // font atlas and noise gen texture
			glTexImage2D(GL_TEXTURE_2D, 0, calc_internal_format(), width, height, 0, calc_format(), get_data_format(), data);
		}
		if (use_mipmaps == 1 || use_mipmaps == 2) {
			if (use_custom_compress) {create_compressed_mipmaps();}
			else {gen_mipmaps();}
		}
		else if (use_mipmaps == 3 || use_mipmaps == 4) {create_custom_mipmaps();}
	}
	if (free_after_upload) {free_client_mem();}
}

void texture_t::upload_cube_map_face(unsigned ix) {
	assert(ix < 6);
	assert(ncolors == 3);
	glTexImage2D((GL_TEXTURE_CUBE_MAP_POSITIVE_X + ix), 0, calc_internal_format(), width, height, 0, calc_format(), get_data_format(), data);
}

void texture_t::calc_color() { // incorrect in is_16_bit_gray mode

	if (normal_map) {color = WHITE; return;} // color not used for normal maps, set to white
	if (defer_load() && !is_allocated()) {color = WHITE; return;} // texture not loaded - this is the best we can do
	if (color != DEF_TEX_COLOR) return; // color already calculated; happens for leaf textures
	//highres_timer_t timer("Texture Color"); // 519ms
	assert(is_allocated());
	float colors[4] = {0.0}, weight(0.0);
	unsigned const size(num_pixels());
	has_binary_alpha = 1;

	if (ncolors == 4) { // RGBA - with alpha component
		for (unsigned i = 0; i < size; ++i) {
			int const offset(i*ncolors);
			unsigned char const alpha(data[offset+3]);
			float const cscale(alpha); // alpha scale
			UNROLL_3X(colors[i_] += cscale*data[offset+i_];);
			colors[3] += alpha;
			weight += cscale;
			if (alpha > 0 && alpha < 255) {has_binary_alpha = 0;}
		} // for i
		UNROLL_3X(color[i_] = colors[i_]/(255.0*weight);)
		color.alpha = CLIP_TO_01(colors[3]/(255.0f*size));
	}
	else {
		unsigned icolors[3] = {};

		if (ncolors == 1) { // grayscale luminance
			for (unsigned i = 0; i < size; ++i) {icolors[0] += data[i];}
			icolors[1] = icolors[2] = icolors[0]; // set G and B to the calculated value of R
		}
		else { // RGB
			for (unsigned i = 0; i < size; ++i) {UNROLL_3X(icolors[i_] += data[i*ncolors+i_];);}
		}
		UNROLL_3X(color[i_] = icolors[i_]/(255.0*size);) // all weights are 1.0
		color.alpha = 1.0;
	}
	if (color.A == 1.0 && (use_mipmaps == 3 || use_mipmaps == 4)) {use_mipmaps = 1;} // no need for custom alpha mipmaps (alpha channel is always 1)
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
	for (unsigned i = 0; i < npixels; ++i) {data[4*i+3] = at.data[at.ncolors*i+alpha_offset];} // copy alpha values
}


void texture_t::merge_in_alpha_channel(texture_t const &at) {

	assert(ncolors == 3 && at.ncolors == 1);
	assert(width == at.width && height == at.height);
	add_alpha_channel();
	copy_alpha_from_texture(at, 0);
}


void texture_t::set_to_color(colorRGBA const &c) {

	assert(is_allocated());
	if (c == color) return; // already set
	if (c == ALPHA0 && (orig_data == NULL || data == orig_data)) return; // color disabled (but never enabled)
	color = c;
	gl_delete();
	if (c == ALPHA0) {data = orig_data; return;} // color disabled
	assert(ncolors == 3 || ncolors == 4);
	color_wrapper c4;
	c4.set_c4(c);
	unsigned const size(num_pixels());
	float const cw_scale(1.0f/(float(c4.c[0]) + float(c4.c[1]) + float(c4.c[2])));
	if (colored_data == NULL) {colored_data = new unsigned char[num_bytes()];}
	if (orig_data    == NULL) {orig_data    = data;} // make a copy
	data = colored_data;
	assert(data != NULL);

	for (unsigned i = 0; i < size; ++i) {
		unsigned const pos(i*ncolors);
		unsigned char *d(data + pos);
		float const cscale(min(1.0f, (unsigned(d[0]) + unsigned(d[1]) + unsigned(d[2]))*cw_scale));
		
		for (int n = 0; n < ncolors; ++n) {
			d[n] = (unsigned char)min(255.0f, (0.5f*cscale*c4.c[n] + 0.5f*orig_data[pos+n]));
		}
	} // for i
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
	bool const is_alpha_mask(index == BLUR_TEX || index == SBLUR_TEX || index == BLUR_CENT_TEX || index == SMOKE_PUFF_TEX || (index >= FLARE1_TEX && index <= FLARE5_TEX));
	bool const is_alpha_tex(index == EXPLOSION_TEX || index == FIRE_TEX || is_alpha_mask);
	bool has_zero_alpha(0);
	assert(is_allocated());

	for (unsigned i = 0; i < size; ++i) {
		int const i4(i << 2);
		unsigned char *buf(data+i4);

		if (index == CLOUD_TEX || index == CLOUD_RAW_TEX) {
			// white -> alpha = 255
			// blue  -> alpha = 0
			float const val((float)buf[0] + (float)buf[1]);
			alpha = ((val <= 340.0) ? 0 : ((unsigned char)1.0*(val - 340.0)));

			if (SMOOTH_SKY_POLES > 0.0 && index == CLOUD_TEX) {
				unsigned const y(i/width);
				float const dist(float(y)/float(height)), d2(min(dist, (1.0f - dist)));
				if (d2 < SMOOTH_SKY_POLES) {alpha = (unsigned char)(d2*d2*float(alpha)/(SMOOTH_SKY_POLES*SMOOTH_SKY_POLES));}
			}
		}
		else { // make white/black part transparent, for example leaves
			float const val((float)buf[0] + (float)buf[1] + (float)buf[2]);
				
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
					float const thresh((index == PLANT4_TEX) ? 660 : 600.0);
					alpha = ((val > thresh) ? 0 : ((val < thresh-100.0) ? 255 : (unsigned char)(2.55f*(thresh - val))));
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
		unsigned char avg_rgb[3] = {};
		UNROLL_3X(avg_rgb[i_] = (unsigned char)(255*color[i_]);)

		// set all alpha=0 texels to the average non-transparent color to improve mipmap quality
		for (unsigned i = 0; i < size; ++i) {
			int const i4(i << 2);
			if (data[i4+3] == 0) {RGB_BLOCK_COPY((data+i4), avg_rgb);} // alpha == 0
		}
	}
}

void texture_t::fill_to_grayscale_color(unsigned char color_val) {
	assert(ncolors == 4);
	assert(is_allocated());
	unsigned const size(num_pixels());
	for(unsigned i = 0; i < size; ++i) {UNROLL_3X(data[(i<<2)+i_] = color_val;);}
}


void texture_t::fill_transparent_with_avg_color() { // unused

	assert(ncolors == 4);
	assert(is_allocated());
	unsigned const size(num_pixels());
	colorRGBA avg_color(0.0, 0.0, 0.0, 1.0);

	for(unsigned i = 0; i < size; ++i) {
		float const cscale(data[(i<<2)+3]/255.0); // alpha scale
		avg_color.A += cscale;
		UNROLL_3X(avg_color[i_] += cscale*data[(i<<2)+i_];);
	}
	UNROLL_3X(avg_color[i_] /= avg_color.A;);
	color_wrapper cw;
	cw.set_c3(avg_color);

	for(unsigned i = 0; i < size; ++i) {
		if (data[(i<<2)+3] == 0) {UNROLL_3X(data[(i<<2)+i_] = cw.c[i_];);} // reassign transparent pixels
	}
}


void texture_t::do_invert_y() {

	assert(is_allocated());
	unsigned const h2(height >> 1), wc(ncolors*width);
		
	for(unsigned i = 0; i < h2; ++i) {
		unsigned const off1(i*wc), off2((height-i-1)*wc);
		for(unsigned j = 0; j < wc; ++j) {swap(data[off1+j], data[off2+j]);} // invert y
	}
}


void texture_t::fix_word_alignment() {

	unsigned const byte_align = 4;
	if ((ncolors*width & (byte_align-1)) == 0) return; // nothing to do
	//timer_t timer("Resize " + name);
	assert(is_allocated());
	float const ar(float(width)/float(height));
	int const new_w(width - (width&(byte_align-1)) + byte_align); // round up to next highest multiple of byte_align
	int const new_h(int(new_w/ar + 0.5)); // preserve aspect ratio
	//cout << "resize " << name << " from " << width << " to " << new_w << ", invert_y: " << invert_y << endl;
	//timer_t timer("Texture Resize");
	resize(new_w, new_h);
	//if (ncolors == 3) {write_to_jpg(name);} // use this one when writing external textures; Warning: may need to set invert_y to get the orient to be correct
	//if (ncolors == 3) {write_to_jpg("textures\\"+name);} // auto-update of resized textures
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

void texture_t::expand_grayscale_to_rgb() {

	if (is_16_bit_gray) return; // 16-bit RGBA not supported; error?
	if (ncolors != 1  ) return; // not grayscale
	assert(is_allocated());
	unsigned const npixels(num_pixels());
	ncolors = 3; // make RGBA
	unsigned char *new_data(new unsigned char[num_bytes()]);
	for (unsigned i = 0; i < npixels; ++i) {UNROLL_3X(new_data[3*i+i_] = data[i];);}
	free_data();
	data = new_data;
}

void texture_t::resize(int new_w, int new_h) { // Note: not thread safe

	if (new_w == width && new_h == height) return; // already correct size
	
	if (omp_get_thread_num_3dw() != 0) {
		std::cerr << "Error: Can't resize texture '" << name << "' on an OpenMP worker thread" << endl;
		exit(1); // for now, this is fatal
		return; // don't scale it?
	}
	assert(is_allocated());
	assert(width > 0 && height > 0 && new_w > 0 && new_h > 0);
	unsigned char *new_data(new unsigned char[new_w*new_h*ncolors]);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // required to handle calls from from fix_word_alignment()
	int const ret(gluScaleImage(calc_format(), width, height, get_data_format(), data, new_w, new_h, get_data_format(), new_data));
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	if (ret) {cout << "GLU error during image scale: " << gluErrorString(ret) << "." << endl;}
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
	for (unsigned i = 0; i < npixels && is_lum; ++i) {is_lum &= (data[3*i+1] == data[3*i] && data[3*i+2] == data[3*i]);}
	if (!is_lum) return 0;
	// RGB equal, make it a single color (luminance) channel
	ncolors = 1; // make luminance
	unsigned char *new_data(new unsigned char[num_bytes()]);
	for (unsigned i = 0; i < npixels; ++i) {new_data[i] = data[3*i];}
	free_data();
	data = new_data;
	return 1;
}


void texture_t::make_normal_map() {

	normal_map = 1;
	if (ncolors == 3) return; // already a normal map
	
	if (ncolors == 4 || ncolors == 2) { // better not have an alpha component
		cout << "Error: Skipping RGBA Bump/Normal map " << name << ", which has " << ncolors << " color components." << endl;
		return;
	}
	assert(is_allocated());
	assert(ncolors == 1 && !is_16_bit_gray); // 8-bit grayscale heightmap
	ncolors = 3; // convert to RGB
	unsigned char *new_data(new unsigned char[num_bytes()]);
	int max_delta(1);

	for (int y = 0; y < height; ++y) { // assume texture wraps
		int const ym1((y == 0) ? height-1 : y-1), yp1((y+1 == height) ? 0 : y+1);

		for (int x = 0; x < width; ++x) {
			int const xm1((x == 0) ? width-1 : x-1), xp1((x+1 == width) ? 0 : x+1);
			max_delta = max(max_delta, abs((int)data[xp1+y*width] - (int)data[xm1+y*width]));
			max_delta = max(max_delta, abs((int)data[x+yp1*width] - (int)data[x+ym1*width]));
		}
	}
	float max_delta_inv(1.0/float(max_delta));

	for (int y = 0; y < height; ++y) { // assume texture wraps
		int const ym1((y == 0) ? height-1 : y-1), yp1((y+1 == height) ? 0 : y+1);

		for (int x = 0; x < width; ++x) {
			int const xm1((x == 0) ? width-1 : x-1), xp1((x+1 == width) ? 0 : x+1), off(3*(x + y*width));
			vector3d n(-((int)data[xp1+y*width] - (int)data[xm1+y*width])*max_delta_inv,
				        ((int)data[x+yp1*width] - (int)data[x+ym1*width])*max_delta_inv, 1.0);
			n.normalize();
			if (invert_bump_maps) {n.x = -n.x; n.y = -n.y;}
			UNROLL_3X(new_data[off+i_] = (unsigned char)(127.5f*(n[i_]+1.0)););
		}
	}
	free_data();
	data = new_data;
}


void texture_t::load_from_gl() { // also set tid?
	alloc();
	glGetTexImage(GL_TEXTURE_2D, 0, calc_format(), get_data_format(), data);
}


void bind_1d_texture(unsigned tid, bool is_array) {
	glBindTexture((is_array ? GL_TEXTURE_1D_ARRAY : GL_TEXTURE_1D), tid);
}

void bind_2d_texture(unsigned tid, bool is_array, bool multisample) {
	glBindTexture(get_2d_texture_target(is_array, multisample), tid);
	//assert(glIsTexture(tid)); // too slow?
}

void bind_cube_map_texture(unsigned tid, bool is_array) {
	glBindTexture((is_array ? GL_TEXTURE_CUBE_MAP_ARRAY : GL_TEXTURE_CUBE_MAP), tid);
	assert(glIsTexture(tid));
}


// 2D texture
void setup_texture(unsigned &tid, bool mipmap, bool wrap_s, bool wrap_t, bool mirror_s, bool mirror_t, bool nearest, float anisotropy, bool is_array, bool multisample) {

	assert(tid == 0);
	assert(!nearest || !mipmap);
	glGenTextures(1, &tid);
	check_gl_error(600);

	// select our current texture
	bind_2d_texture(tid, is_array, multisample);
	if (multisample) return; // don't set any other state
	int const target(get_2d_texture_target(is_array, multisample));

	// when texture area is small, use linear filter (bilinear filter the closest mipmap)
	glTexParameteri(target, GL_TEXTURE_MIN_FILTER, (nearest ? GL_NEAREST : (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR))); // GL_LINEAR_MIPMAP_NEAREST?
	// when texture area is large, bilinear filter the first mipmap
	glTexParameteri(target, GL_TEXTURE_MAG_FILTER, (nearest ? GL_NEAREST : GL_LINEAR));
	// enable anisotropic filtering (slower but higher quality)
	if (anisotropy > 1.0) {glTexParameterf(target, GL_TEXTURE_MAX_ANISOTROPY_EXT, anisotropy);}
	// this may be useful for model textures such as the Mineways texture atlas;
	// however, we can't use a global config here and would need to store an option per texture and pass it down through all of the function calls
	//if (mipmap) {glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LOD, max_mipmap_levels);}

	// if wrap is true,  the texture wraps over at the edges (repeat) or is mirrored
	// if wrap is false, the texture ends at the edges (clamp)
	int const mode_s(wrap_s ? (mirror_s ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE); // Note: clamp is more efficient than wrap
	int const mode_t(wrap_t ? (mirror_t ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE);
	glTexParameteri(target, GL_TEXTURE_WRAP_S, mode_s);
	glTexParameteri(target, GL_TEXTURE_WRAP_T, mode_t);
}


void setup_1d_texture(unsigned &tid, bool mipmap, bool wrap, bool mirror, bool nearest) {

	assert(tid == 0);
	assert(!nearest || !mipmap);
	glGenTextures(1, &tid);
	bind_1d_texture(tid);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, (nearest ? GL_NEAREST : (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR))); // GL_LINEAR_MIPMAP_NEAREST?
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, (nearest ? GL_NEAREST : GL_LINEAR));
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, (wrap ? (mirror ? GL_MIRRORED_REPEAT : GL_REPEAT) : GL_CLAMP_TO_EDGE));
}


void setup_cube_map_texture(unsigned &tid, unsigned tex_size, bool allocate, bool use_mipmaps, float aniso) { // Note: no mipmaps

	assert(tid == 0);
	glGenTextures(1, &tid);
	bind_cube_map_texture(tid);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, (use_mipmaps ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR));
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAX_ANISOTROPY_EXT, aniso);
	
	if (allocate) {
		assert(tex_size > 0);
		unsigned num_levels(1);

		if (use_mipmaps) { // determine the number of mipmap levels needed
			for (unsigned sz = tex_size; sz > 1; sz >>= 1) {++num_levels;}
		}
		glTexStorage2D(GL_TEXTURE_CUBE_MAP, num_levels, GL_RGB8, tex_size, tex_size);
	}
	//gluBuild2DMipmaps(GL_TEXTURE_CUBE_MAP_POSITIVE_X, GL_RGB8, tex_size, tex_size, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
	//glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
	check_gl_error(520);
}

void frame_buffer_to_texture(unsigned &tid, bool is_depth) {

	if (tid) {bind_2d_texture(tid);} else {setup_texture(tid, 0, 0, 0);}
	glTexImage2D(GL_TEXTURE_2D, 0, (is_depth ? GL_DEPTH_COMPONENT32F : GL_RGB8), window_width, window_height, 0,
		(is_depth ? GL_DEPTH_COMPONENT : GL_RGB), (is_depth ? GL_FLOAT : GL_UNSIGNED_BYTE), nullptr);
	glReadBuffer(GL_BACK);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, window_width, window_height);
}

void depth_buffer_to_texture    (unsigned &tid) {frame_buffer_to_texture(tid, 1);}
void frame_buffer_RGB_to_texture(unsigned &tid) {frame_buffer_to_texture(tid, 0);}


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
	if (!bbb_tex.is_allocated()) return; // not yet loaded/allocated, or texture not found; skip this update
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


void gen_building_window_texture(float width_frac, float height_frac) { // Note: generated when needed, not during load

	for (unsigned n = 0; n < 2; ++n) { // generate both opaque and transparent building window textures
		unsigned const tid(n ? (unsigned)BLDG_WIND_TRANS_TEX : (unsigned)BLDG_WINDOW_TEX);
		texture_t &tex(textures[tid]);
		assert(tex.ncolors == 4); // or 2 color grayscale + alpha?
		unsigned char *tex_data(tex.get_data());
		assert(width_frac  > 0.0 && width_frac  < 1.0);
		assert(height_frac > 0.0 && height_frac < 1.0);
		float const xspace(0.5*(1.0 - width_frac)), yspace(0.5*(1.0 - height_frac));
		int const w1(round_fp(xspace*tex.width)), w2(round_fp((1.0 - xspace)*tex.width)), h1(round_fp(yspace*tex.height)), h2(round_fp((1.0 - yspace)*tex.height));
		int const border(tex.width/32 + n), w1b(w1 - border), w2b(w2 + border), h1b(h1 - border), h2b(h2 + border); // window borders

		for (int i = 0; i < tex.height; ++i) {
			for (int j = 0; j < tex.width; ++j) {
				int const offset(4*(i*tex.width + j));
				bool const pane(i > h1 && i <= h2 && j > w1 && j <= w2), frame(i > h1b && i <= h2b && j > w1b && j <= w2b);
				unsigned char luminance(pane ? 128 : 255); // gray for window, white for border
				UNROLL_3X(tex_data[offset+i_] = luminance;)
				tex_data[offset+3] = (frame ? ((!pane || n==0) ? 255 : 32) : 0); // alpha: partially transparent inside window, opaque border, and transparent outside
			}
		}
		tex.init(); // called later in the flow, init it now
	} // for n
}


void noise_fill(unsigned char *data, unsigned size) {
	rand_gen_t rgen;
	for (unsigned i = 0; i < size; ++i) {data[i] = (rgen.rand() & 255);}
}

void noise_fill_01(vector<float> &data) {
	rand_gen_t rgen;
	for (auto i = data.begin(); i != data.end(); ++i) {*i = rgen.rand_float();}
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


unsigned create_3d_noise_texture(unsigned size, unsigned ncomp, unsigned bytes_per_pixel) {

	if (bytes_per_pixel == 4) { // special case - use floating-point texture
		assert(ncomp == 1); // grayscale only for now
		vector<float> data(ncomp*size*size*size);
		noise_fill_01(data);
		unsigned tid(0);
		setup_3d_texture(tid, GL_LINEAR, GL_REPEAT);
		glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, size, size, size, 0, get_texture_format(ncomp), GL_FLOAT, &data.front());
		return tid;
	}
	vector<unsigned char> data(ncomp*bytes_per_pixel*size*size*size);
	noise_fill(&data.front(), data.size());
	return create_3d_texture(size, size, size, ncomp, data, GL_LINEAR, GL_REPEAT, 0, bytes_per_pixel); // compressed?
}


unsigned get_noise_tex_3d(unsigned tsize, unsigned ncomp, unsigned bytes_per_pixel) {

	pair<texture_map_t::iterator, bool> ret(noise_tex_3ds.insert(make_pair(make_pair(tsize, (ncomp + 4*bytes_per_pixel)), 0)));
	if (ret.second) {ret.first->second = create_3d_noise_texture(tsize, ncomp, bytes_per_pixel);}
	return ret.first->second;
}


colorRGBA get_landscape_texture_color(int xpos, int ypos) {

	if (cached_ls_colors.empty()) {cached_ls_colors.resize(XY_MULT_SIZE, ALPHA0);}
	unsigned const ix(xpos + MESH_X_SIZE*ypos);
	assert(ix < cached_ls_colors.size());
	if (cached_ls_colors[ix].alpha != 0.0) {return cached_ls_colors[ix];}
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


int get_bare_ls_tid(float zval) {
	float const relh(relh_adj_tex + (zval - zmin)/(zmax - zmin));
	return ((relh > clip_hd1) ? (int)ROCK_TEX : (int)DIRT_TEX); // rock or dirt
}

void update_lttex_ix(int &ix) { // note: assumes lttex_dirt
	if ((water_is_lava || DISABLE_WATER == 2) && lttex_dirt[ix].id == SNOW_TEX) {--ix;} // convert snow to rock
	if (vegetation == 0.0 && lttex_dirt[ix].id == GROUND_TEX) {++ix;} // convert ground/grass to rock
}

void get_tids(float relh, int &k1, int &k2, float *t) {

	// find first texture with height greater than relh; unrolled for efficiency
	if      (relh < h_dirt[0]) {k1 = 0;}
	else if (relh < h_dirt[1]) {k1 = 1;}
	else if (relh < h_dirt[2]) {k1 = 2;}
	else if (relh < h_dirt[3]) {k1 = 3;}
	else                       {k1 = 4;}

	if (k1 < NTEX_DIRT-1 && (h_dirt[k1] - relh) < TEXTURE_SMOOTH) {
		if (t) {*t = 1.0 - (h_dirt[k1] - relh)/TEXTURE_SMOOTH;}
		k2 = k1+1;
		update_lttex_ix(k1);
		update_lttex_ix(k2);
	}
	else {
		update_lttex_ix(k1);
		k2 = k1;
	}
}


void clear_cached_ls_colors() {
	for (auto i = cached_ls_colors.begin(); i != cached_ls_colors.end(); ++i) {*i = ALPHA0;}
}


void create_landscape_texture() {

	RESET_TIME;
	if (using_custom_landscape_texture() || world_mode == WMODE_INF_TERRAIN) return;
	clear_cached_ls_colors();
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
	int const def_id((default_ground_tex >= 0) ? default_ground_tex : GROUND_TEX);
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
			get_tids(relh1, k1a, k2a);
			get_tids(relh2, k1b, k2b);
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
					float const mh((1.0f - xpi)*((1.0f - ypi)*mh00 + ypi*mh10) + xpi*((1.0f - ypi)*mh01 + ypi*mh11));
					float const relh(relh_adj_tex + (mh - zmin)*dz_inv);
					get_tids(relh, k1, k2, &t);
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
			float const vnz((1.0f - xpi)*((1.0f - ypi)*vnz00 + ypi*vnz10) + xpi*((1.0f - ypi)*vnz01 + ypi*vnz11));

			if (grass && vnz < sti[1]) { // ground/grass
				texture_t const &ta(textures[DIRT_TEX]);
				unsigned char const *ta_data(ta.get_data());
				int const tofa(ta.ncolors*(((i+toy)&(ta.height-1))*ta.width + ((j+tox)&(ta.width-1))));
				unsigned char temp[3] = {};

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
		tex.update_texture_data(0, 0, width, height);
	}
	else {
		tex.gl_delete(); // should we try to update rather than recreating from scratch?
		tex.do_gl_init();
		tex.calc_color();
	}
	if (!scrolling) {PRINT_TIME(" Gen Landscape Texture");}
}


inline int get_align(int x1, int x2, unsigned am) {
	return ((1 + am - ((x2 - x1) & am)) & am);
}

void texture_t::update_texture_data(int x1, int y1, int x2, int y2) {

	if (x1 == x2 || y1 == y2) return; // nothing to update
	assert(ncolors == 3 && width > 0 && height > 0 && is_allocated());

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
	check_init();
	bind_gl();
	glPixelStorei(GL_UNPACK_ROW_LENGTH, width);
	glTexSubImage2D(GL_TEXTURE_2D, 0, x1, y1, (x2-x1), (y2-y1), calc_format(), GL_UNSIGNED_BYTE, (get_data() + ncolors*(x1 + y1*width)));
	if (LANDSCAPE_MIPMAP) {gen_mipmaps(2);} // update mipmaps if needed (non-sparse)
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0); // reset to 0
}


void regrow_landscape_texture_amt0() {

	//RESET_TIME;
	static int counter(0);
	if (using_custom_landscape_texture() || iticks == 0) return; // is it too strong to use both iticks and fticks?
	if (ls0_invalid) create_landscape_texture();
	texture_t &tex(textures[LANDSCAPE_TEX]);
	assert(tex.is_allocated() && landscape0 != NULL);
	int const regen_bshift(max(1, min(7, int(log2(1.0/max(0.0001f, fticks*LANDSCAPE_REGEN_AMT))))));
	unsigned char *tex_data(tex.get_data());
	int const y1((counter%LANDSCAPE_REGEN_MOD)*(tex.height/LANDSCAPE_REGEN_MOD));
	int const y2(y1 + (tex.height/LANDSCAPE_REGEN_MOD));
	assert(y2 <= tex.height);

	//if (LANDSCAPE_REGEN_AMT > 0.0 || skip_regrow) {
		for (int i = y1; i < y2; ++i) {
			int const i_step(i*tex.width);

			for (int j = 0; j < tex.width; ++j) { // performance critical
				int const offset(3*(i_step + j));
				UNROLL_3X(tex_data[offset+i_] += (unsigned char)(((int)landscape0[offset+i_] - (int)tex_data[offset+i_]) >> regen_bshift);)
			}
		}
	//}
	++counter;
	skip_regrow = 0;
	tex.update_texture_data(0, y1, tex.width, y2);
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
	clear_cached_ls_colors();
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
	int const rad(max(0, int(5.0f*radius/((xscale + yscale)*(DX_VAL + DY_VAL))) - 1)), rad_sq(max(1, rad*rad)); // size in texture space
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
				float const blend(0.8f*color.alpha*(1.0f - float(dist_sq)*blend_scale));
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
			float const mult((frad - (float)dist)*inv_frad), acc255s(acc255*mult), omacc(1.0f - acc*mult);
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


void update_lt_section(int x1, int y1, int x2, int y2) {

	//RESET_TIME;
	if (x1 == x2 || y1 == y2) return;
	assert(!using_custom_landscape_texture());
	textures[LANDSCAPE_TEX].update_texture_data(x1, y1, x2, y2);
	//PRINT_TIME("LT Update");
}


int snow_height(point pos) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos)) return 0;
	float const relh(relh_adj_tex + (mesh_height[ypos][xpos] - zmin)/(zmax - zmin));
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

void get_poly_texgen_dirs(vector3d const &norm, vector3d v[2]) { // similar to get_ortho_vectors(), but not normalized

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
		vector4d const vv(scale[i]*v[i].x, scale[i]*v[i].y, scale[i]*v[i].z, (xlate[i] + scale[i]*dot_product(offset, v[i])));
		set_texgen_vec4(vv, ((i != 0) ^ swap_txy), shader, mode);
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
	int tx(int(width*u) % width), ty(int(height*v) % height); // width and height can be any nonzero value
	if (tx < 0) {tx += width;}
	if (ty < 0) {ty += height;}
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

unsigned texture_t::get_gpu_mem() const {
	if (!is_bound()) return 0;
	unsigned mem(0);
	if (is_texture_compressed()) {mem = (num_pixels()/16) * ((ncolors == 4) ? 16 : 8);} // assumes DXT1 8:1 RGB and DXT5 4:1 RGBA compression
	else {mem = num_bytes();}
	if (use_mipmaps) {mem += mem/3;} // 33% overhead
	return mem;
}

texture_t const &get_texture_by_id(unsigned tid) {
	assert(tid < textures.size());
	return textures[tid];
}
colorRGBA texture_color(int tid) {
	return ((tid < 0) ? WHITE : get_texture_by_id(tid).get_avg_color());
}
unsigned get_texture_size(int tid, bool dim) {
	assert(tid >= 0);
	return (dim ? get_texture_by_id(tid).height : get_texture_by_id(tid).width);
}
float get_texture_component(unsigned tid, float u, float v, int comp) {
	return get_texture_by_id(tid).get_component(u, v, comp);
}
float get_texture_component_grayscale_pow2(unsigned tid, float u, float v) {
	return get_texture_by_id(tid).get_component_grayscale_pow2(u, v);
}
colorRGBA get_texture_color(unsigned tid, float u, float v) {
	return get_texture_by_id(tid).get_texel(u, v);
}
int get_texture_normal_map_tid(unsigned tid) {
	return get_texture_by_id(tid).bump_tid;
}

void texture_t::write_pixel_16_bits(unsigned ix, float val) { // Note: no error checking
	unsigned char const high_bits(val); // high bits - truncate
	data[(ix<<1)+1] = high_bits;
	data[ix<<1]     = (unsigned char)(256.0f*(val - float(high_bits))); // low bits - remainder
}


vector2d get_billboard_texture_uv(point const *const points, point const &pos) {

	assert(points != NULL);
	// ordered: (s,t) => (0,1), (0,0), (1,0), (1,1)
	float d[4] = {}; // distance from coll point to quad edge

	for (unsigned i = 0; i < 4; ++i) {
		unsigned const in((i+1)&3);
		d[i] = cross_product((pos - points[i]), (pos - points[in])).mag()/p2p_dist(points[i], points[in]);
	}
	assert(d[0] + d[2] > 0.0f);
	assert(d[1] + d[3] > 0.0f);
	vector2d uv(d[0]/(d[0] + d[2]), d[1]/(d[1] + d[3])); // y is upside down
	uv.x = CLIP_TO_01(uv.x); uv.y = CLIP_TO_01(uv.y); // clamp uv to account for fp rounding errors (and incorrect pos?)
	//assert(uv.x >= 0.0 && uv.x <= 1.0 && uv.y >= 0.0 && uv.y <= 1.0); // more restrictive case
	return uv;
}


bool is_billboard_texture_transparent(point const *const points, point const &pos, int tid) {
	vector2d const uv(get_billboard_texture_uv(points, pos));
	return (get_texture_component(tid, uv.x, uv.y, 3) == 0.0);
}

void ensure_texture_loaded(unsigned &tid, unsigned txsize, unsigned tysize, bool mipmap, bool nearest, bool multisample) { // used with texture_pair_t/render-to-texture RGBA
	assert(txsize > 0 && tysize > 0);
	if (tid) return; // already created
	setup_texture(tid, mipmap, 0, 0, 0, 0, nearest, 1.0, 0, multisample);
	if (multisample) {glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, NUM_TEX_MS_SAMPLES, GL_RGBA8, txsize, tysize, false);}
	else {glTexImage2D(get_2d_texture_target(0, multisample), 0, GL_RGBA8, txsize, tysize, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);}
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
	bind_texture_tu(tids[0], 0);
	bind_texture_tu(tids[1], 1);
}
void texture_pair_t::ensure_tid(unsigned tsize, bool mipmap) {
	ensure_texture_loaded(tids[0], tsize, tsize, mipmap, 0, multisample); // color
	ensure_texture_loaded(tids[1], tsize, tsize, mipmap, 0, multisample); // normal
}

void texture_atlas_t::free_context() {free_texture(tid);}

void texture_atlas_t::bind_texture() const {
	assert(tid);
	bind_2d_texture(tid, 0, multisample);
}
void texture_atlas_t::ensure_tid(unsigned base_tsize, bool mipmap) {
	assert(nx > 0 && ny > 0);
	ensure_texture_loaded(tid, nx*base_tsize, ny*base_tsize, mipmap, 0, multisample);
}


