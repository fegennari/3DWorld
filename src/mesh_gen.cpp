// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/1/02

#include "3DWorld.h"
#include "mesh.h"
#include "textures.h"
#include "sinf.h"
#include "heightmap.h"
#include "shaders.h"
#include <glm/gtc/noise.hpp>


int      const NUM_FREQ_COMP      = 9;
float    const MESH_SCALE_Z_EXP   = 0.7;
int      const N_RAND_SIN2        = 10;
int      const FREQ_FILTER        = 2; // higher = smoother landscape
int      const MIN_FREQS          = 3;
float    const W_PLANE_Z          = 0.42;
float    const HEIGHT_SCALE       = 0.01;
unsigned const EST_RAND_PARAM     = 128;
float    const READ_MESH_H_SCALE  = 0.0008;
float    const MESH_SCALE_FACTOR  = 0.0007;
bool     const AUTOSCALE_HEIGHT   = 1;
bool     const DEF_GLACIATE       = 1;
float    const DEF_GLACIATE_EXP   = 3.0;
bool     const GEN_SCROLLING_MESH = 1;
float    const S_GEN_ATTEN_DIST   = 128.0;

int   const F_TABLE_SIZE = NUM_FREQ_COMP*N_RAND_SIN2;


// Global Variables
float MESH_START_MAG(0.02), MESH_START_FREQ(240.0), MESH_MAG_MULT(2.0), MESH_FREQ_MULT(0.5);
int cache_counter(1), start_eval_sin(0), GLACIATE(DEF_GLACIATE), mesh_gen_mode(MGEN_SINE), mesh_gen_shape(0), mesh_freq_filter(FREQ_FILTER);
float zmax, zmin, zmax_est, zcenter(0.0), zbottom(0.0), ztop(0.0), h_sum(0.0), alt_temp(DEF_TEMPERATURE);
float mesh_scale(1.0), tree_scale(1.0), mesh_scale_z(1.0), mesh_scale_z_inv(1.0), glaciate_exp(1.0), glaciate_exp_inv(1.0);
float mesh_height_scale(1.0), zmax_est2(1.0), zmax_est2_inv(1.0);
vector<float> sin_table;
float sinTable[F_TABLE_SIZE][5];

int const mesh_tids_dirt[NTEX_DIRT] = {SAND_TEX, DIRT_TEX, GROUND_TEX, ROCK_TEX, SNOW_TEX};
float const mesh_rh_dirt[NTEX_DIRT] = {0.40, 0.44, 0.60, 0.75, 1.0};
float sthresh[2][2] = {{0.68, 0.86}, {0.48, 0.72}}; // {grass, snow}, {lo, hi}
ttex lttex_dirt[NTEX_DIRT];
vector<float> height_histogram;
hmap_params_t hmap_params;


extern bool combined_gu;
extern int xoff, yoff, xoff2, yoff2, world_mode, rand_gen_index, mesh_rgen_index, mesh_scale_change, display_mode;
extern int read_heightmap, read_landscape, do_read_mesh, mesh_seed, scrolling, camera_mode, invert_mh_image;
extern unsigned erosion_iters;
extern double c_radius, c_phi, c_theta;
extern float water_plane_z, temperature, mesh_file_scale, mesh_file_tz, custom_glaciate_exp, MESH_HEIGHT, XY_SCENE_SIZE;
extern float water_h_off, water_h_off_rel, disabled_mesh_z, read_mesh_zmm, init_temperature, univ_temp;
extern point mesh_origin, surface_pos;
extern char *mh_filename, *mh_filename_tt, *mesh_file;


void glaciate();
void gen_terrain_map();
void estimate_zminmax(bool using_eq);
void set_zvals();
void update_temperature(bool verbose);
void compute_scale();

bool using_hmap_with_detail();



void create_sin_table() {

	if (!sin_table.empty()) return; // already setup
	sin_table.resize(2*TSIZE);

	for (unsigned i = 0; i < TSIZE; ++i) {
		sin_table[i]       = sinf(i/sscale);
		sin_table[i+TSIZE] = cosf(i/sscale);
	}
}


void matrix_min_max(float **matrix, float &minval, float &maxval) {

	assert(matrix);
	minval = maxval = matrix[0][0];

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			minval = min(minval, matrix[i][j]);
			maxval = max(maxval, matrix[i][j]);
		}
	}
}

void calc_zminmax() {
	matrix_min_max(mesh_height, zmin, zmax);
}


bool bmp_to_chars(char const *const fname, unsigned char **&data) { // Note: supports all image formats

	assert(fname != nullptr);
	texture_t texture(0, 7, MESH_X_SIZE, MESH_Y_SIZE, 0, 1, 0, fname, 0); // invert_y=0
	texture.load(-1, 0, 0, 1); // generates fatal internal errors if load() fails, so return is always 1
	
	if (texture.width != MESH_X_SIZE || texture.height != MESH_Y_SIZE || texture.ncolors != 1) {
		std::cerr << "Error loading BMP file '" << fname << "': Expected image to be " << MESH_X_SIZE << "x" << MESH_Y_SIZE << "x" << 1
			      << " but it was " << texture.width << "x" << texture.height << "x" << texture.ncolors << endl;
		exit(1);
	}
	matrix_gen_2d(data);
	memcpy(data[0], texture.get_data(), XY_MULT_SIZE*sizeof(char));
	texture.free_data();
	return 1;
}


float scale_mh_texture_val  (float val) {return (READ_MESH_H_SCALE*mesh_height_scale*mesh_file_scale*val + mesh_file_tz)*mesh_scale_z_inv;}
float unscale_mh_texture_val(float val) {return (mesh_scale_z*val - mesh_file_tz)/(READ_MESH_H_SCALE*mesh_height_scale*mesh_file_scale);}
float get_mh_texture_mult() {return READ_MESH_H_SCALE*mesh_height_scale*mesh_file_scale*mesh_scale_z_inv;}
float get_mh_texture_add () {return mesh_file_tz*mesh_scale_z_inv;}
float get_heightmap_scale() {return 256.0*READ_MESH_H_SCALE;}

void set_mesh_height_scales_for_zval_range(float min_z, float dz) {
	assert(dz > 0.0);
	// undo the transform: v = (h - min_z)/dz by applying h = get_mh_texture_mult()*v + get_mh_texture_add(); h = M*v + A; M = dz, A = min_z
	mesh_file_scale = dz/(READ_MESH_H_SCALE*mesh_height_scale*mesh_scale_z_inv);
	mesh_file_tz    = min_z/mesh_scale_z_inv;
}


// Note: only works for 8-bit heightmaps (higher precision textures are truncated)
bool read_mesh_height_image(char const *fn, bool allow_resize=1) {

	if (fn == NULL) {
		std::cerr << "Error: No mh_filename spcified in the config file." << endl;
		return 0;
	}
	//timer_t timer("Read Mesh Heightmap");
	cout << "Reading mesh heightmap " << mh_filename << endl;
	heightmap_t hmap(0, 7, MESH_X_SIZE, MESH_Y_SIZE, fn, (invert_mh_image != 0));
	hmap.load(-1, allow_resize, 1, 1); // allow 2-byte grayscale (currently only works for PNGs)
	if (allow_resize) {hmap.resize(MESH_X_SIZE, MESH_Y_SIZE);}
		
	if (hmap.width != MESH_X_SIZE || hmap.height != MESH_Y_SIZE) { // may be due to texture padding to make a multipe of 4
		std::cerr << "Error reading mesh height image: Expected size " << MESH_X_SIZE << "x" << MESH_Y_SIZE
				    << ", got size " << hmap.width << "x" << hmap.height << endl;
		return 0;
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			mesh_height[i][j] = scale_mh_texture_val(hmap.get_heightmap_value(j, i));
		}
	}
	hmap.free_data();
	return 1;
}


void set_zmax_est(float zval) {

	zmax_est      = zval;
	zmax_est2     = 2.0*zmax_est;
	zmax_est2_inv = 1.0/zmax_est2;
}


void update_disabled_mesh_height() {

	if (disabled_mesh_z == FAR_DISTANCE || mesh_draw == NULL) return;

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (is_mesh_disabled(j, i)) {mesh_height[i][j] = disabled_mesh_z;}
		}
	}
}


void clamp_to_mesh(int xy[2]) {
	for (unsigned i = 0; i < 2; ++i) {xy[i] = max(0, min(MESH_SIZE[i]-1, xy[i]));}
}


void gen_mesh_random(float height) {

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			float val(signed_rand_float2()*height);
			if (i + j == 0) val *= 10.0;
			if (i != 0)     val += 0.45*mesh_height[i-1][j];
			if (j != 0)     val += 0.45*mesh_height[i][j-1];
			mesh_height[i][j] = val;
		}
	}
}


void gen_mesh_sine_table(float **matrix, int x_offset, int y_offset, int xsize, int ysize) {

	assert(matrix);
	mesh_xy_grid_cache_t height_gen;
	height_gen.build_arrays(float(x_offset - xsize/2)*DX_VAL, float(y_offset - ysize/2)*DY_VAL, DX_VAL, DY_VAL, xsize, ysize);

	for (int i = 0; i < ysize; ++i) {
		for (int j = 0; j < xsize; ++j) {matrix[i][j] = height_gen.eval_index(j, i);}
	}
}


void apply_mesh_rand_seed(rand_gen_t &rgen) {
	if (mesh_seed != 0) {rgen.set_state(mesh_seed, 12345);}
	else if (mesh_gen_mode != MGEN_SINE) {rgen.set_state(mesh_rgen_index+1, 12345);}
}


void gen_rand_sine_table_entries(float scaled_height) {

	float xf_scale((float)MESH_Y_SIZE/(float)MESH_X_SIZE), yf_scale(1.0/xf_scale);
	if (X_SCENE_SIZE > Y_SCENE_SIZE) yf_scale *= (float)Y_SCENE_SIZE/(float)X_SCENE_SIZE;
	if (Y_SCENE_SIZE > X_SCENE_SIZE) xf_scale *= (float)X_SCENE_SIZE/(float)Y_SCENE_SIZE;
	float mags[NUM_FREQ_COMP] = {}, freqs[NUM_FREQ_COMP] = {};
	// Note: none of these config values are error checked, maybe should at least check >= 0.0
	freqs[0] = MESH_START_FREQ;
	mags [0] = MESH_START_MAG;

	for (int i = 1; i < NUM_FREQ_COMP; ++i) {
		freqs[i] = freqs[i-1]*MESH_FREQ_MULT;
		mags [i] = mags[i-1]*MESH_MAG_MULT;
	}
	h_sum = 0.0;
	for (int l = 0; l < NUM_FREQ_COMP; ++l) {h_sum += mags[l];}
	float const mesh_h(scaled_height/sqrt(0.1*N_RAND_SIN2));
	h_sum *= N_RAND_SIN2*scaled_height*HEIGHT_SCALE;
	static rand_gen_t rgen; // static so that later calls to this function will generate different values
	apply_mesh_rand_seed(rgen);

	for (int l = 0; l < NUM_FREQ_COMP; ++l) {
		int const offset(l*N_RAND_SIN2);
		float const x_freq(freqs[l]/((float)MESH_X_SIZE)), y_freq(freqs[l]/((float)MESH_Y_SIZE));
		float const mheight(mags[l]*mesh_h);

		for (int i = 0; i < N_RAND_SIN2; ++i) {
			int const index(offset + i);
			sinTable[index][0] = rgen.rand_uniform(0.2, 1.0)*mheight; // magnitude
			sinTable[index][1] = rgen.rand_float()*TWO_PI; // y phase
			sinTable[index][2] = rgen.rand_float()*TWO_PI; // x phase
			sinTable[index][3] = rgen.rand_uniform(0.1, 1.0)*x_freq*yf_scale; // y frequency
			sinTable[index][4] = rgen.rand_uniform(0.1, 1.0)*y_freq*xf_scale; // x frequency
		}
	}
}


void gen_mesh(int surface_type, int keep_sin_table, int update_zvals) {

	static bool init(0);
	vector<float> atten_table(MESH_X_SIZE);
	assert((read_landscape != 0) + (read_heightmap != 0) + (do_read_mesh != 0) <= 1);
	assert(surface_type != 2); // random+sin surface is no longer supported
	if (read_landscape) {surface_type = 3;} else if (read_heightmap) {surface_type = 4;} else if (do_read_mesh) {surface_type = 5;}
	bool surface_generated(0);
	bool const loaded_surface(surface_type >= 3);
	bool const gen_scroll_surface(GEN_SCROLLING_MESH && scrolling && loaded_surface);
	float const scaled_height(MESH_HEIGHT*mesh_height_scale);
	++cache_counter; // invalidate mesh cache

	if (surface_type != 5) {
		zmax = -LARGE_ZVAL;
		zmin =  LARGE_ZVAL;
	}
	compute_scale();
	matrix_clear_2d(mesh_height);

	// Note: we always create the sine table, even when using a heightmap, because it may be used for random tree distributions, etc.
	if (!keep_sin_table || !init) {
		surface_generated = 1;
		gen_rand_sine_table_entries(scaled_height);
	}
	if (surface_type == 0 || gen_scroll_surface) { // sine waves
		gen_mesh_sine_table(mesh_height, xoff2, yoff2, MESH_X_SIZE, MESH_Y_SIZE);
	}
	if (surface_type == 1) { // random
		zmax = -LARGE_ZVAL;
		zmin = -zmax;
		gen_mesh_random(0.1*scaled_height);
	}
	else if (surface_type >= 3) { // read mesh image file
		static float **read_mh = NULL;

		if (scrolling) {
			assert(read_mh != NULL);

			if (abs(xoff2) < MESH_X_SIZE && abs(yoff2) < MESH_Y_SIZE) { // otherwise the entire read area is scrolled off the mesh
				float mh_scale(1.0);
				float const atten_dist(S_GEN_ATTEN_DIST*XY_SUM_SIZE/256);

				if (GEN_SCROLLING_MESH) { // still not entirely correct
					float mzmin, mzmax;
					matrix_min_max(mesh_height, mzmin, mzmax);
					mh_scale = CLIP_TO_01(zmax_est/max(fabs(mzmax), fabs(mzmin)));
				}
				for (int i = 0; i < MESH_Y_SIZE; ++i) {
					for (int j = 0; j < MESH_X_SIZE; ++j) { // scroll by absolute offset
						int xy[2] = {(j + xoff2), (i + yoff2)};
						int const xy0[2] = {xy[0], xy[1]};
						clamp_to_mesh(xy);
						float const mh(read_mh[xy[1]][xy[0]]);

						if (GEN_SCROLLING_MESH && (xy[0] != xy0[0] || xy[1] != xy0[1])) {
							float const atten(CLIP_TO_01((abs(xy[0] - xy0[0]) + abs(xy[1] - xy0[1]))/atten_dist));
							mesh_height[i][j] = atten*mh_scale*mesh_height[i][j] + (1.0 - atten)*mh; // not yet correct
						}
						else {
							mesh_height[i][j] = mh;
						}
					}
				}
			} // xoff/yoff test
		}
		else if (read_mh != NULL) {
			memcpy(mesh_height[0], read_mh[0], XY_MULT_SIZE*sizeof(float));
		}
		else {
			if (!((surface_type == 5) ? read_mesh(mesh_file, read_mesh_zmm) : read_mesh_height_image(mh_filename))) {
				std::cerr << "Error reading landscape height data." << endl;
				exit(1);
			}
			if (read_mh == NULL) matrix_gen_2d(read_mh); // should this be deleted later? reread?
			memcpy(read_mh[0], mesh_height[0], XY_MULT_SIZE*sizeof(float));
		}
		//matrix_delete_2d(read_mh);
	} // end read
	update_disabled_mesh_height();
	if (surface_type != 5) {calc_zminmax();}

	if (surface_type != 5) {
		if (!keep_sin_table || !init || update_zvals) {
			if (AUTOSCALE_HEIGHT && world_mode == WMODE_GROUND) {mesh_origin.z = camera_origin.z = surface_pos.z = 0.5f*(zmin + zmax);}
			estimate_zminmax(surface_type < 3);
		}
		else {
			set_zvals();
		}
	}
	update_temperature(1);
	gen_terrain_map();

	if (GLACIATE && world_mode == WMODE_GROUND && (!keep_sin_table || !init || update_zvals) && AUTOSCALE_HEIGHT) {
		mesh_origin.z = camera_origin.z = surface_pos.z = 0.5f*(zbottom + ztop); // readjust camera height
	}
	if (surface_generated) init = 1;
}


float do_glaciate_exp(float value) {
	return ((custom_glaciate_exp == 0.0) ? value*value*value : pow(value, custom_glaciate_exp)); // DEF_GLACIATE_EXP = 3.0
}

float get_rel_wpz() {return CLIP_TO_01(W_PLANE_Z + water_h_off_rel);}

float get_volcano_height(float xi, float yi) {
	float const freq(mesh_scale/hmap_params.volcano_width), x(freq*xi), y(freq*yi), dist(sqrt(x*x + y*y));
	if (dist > 2.0) return 0.0; // too far, no effect (optimization)
	float const val(COSF(x)*COSF(y)), hole(max(0.0, 400.0*(val - 0.999)));
	//float const peak(150.0*(val - 0.99)); // cos function - discontinuous at the edges, repating for each island
	float const peak(0.08*val/max(0.04f, dist)); // sinc function - steep with flat top and wide base - only on center island
	return hmap_params.volcano_height*max(0.0f, (peak - hole))*mesh_scale_z_inv;
}

void apply_mesh_sine(float &zval, float x, float y) {
	if (hmap_params.sine_mag > 0.0) { // Note: snow thresh is still off when highly zoomed in
		float const freq(mesh_scale*hmap_params.sine_freq);
		zval += (hmap_params.sine_mag*COSF(x*freq)*COSF(y*freq) + hmap_params.sine_bias)*mesh_scale_z_inv;
		if (hmap_params.volcano_width > 0.0 && hmap_params.volcano_height > 0.0) {zval += get_volcano_height(x, y);}
	}
}
void apply_glaciate(float &zval) {
	if (GLACIATE) {
		float const relh((zval + zmax_est)*zmax_est2_inv);
		zval = do_glaciate_exp(relh)*zmax_est2 - zmax_est;
	}
}


void glaciate() {

	ztop             = -LARGE_ZVAL;
	zbottom          =  LARGE_ZVAL;
	glaciate_exp     = ((custom_glaciate_exp == 0.0) ? DEF_GLACIATE_EXP : custom_glaciate_exp);
	glaciate_exp_inv = 1.0/glaciate_exp;

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			float &zval(mesh_height[i][j]);
			apply_glaciate(zval);
			apply_mesh_sine(zval, float(j + xoff2 - MESH_X_SIZE/2), float(i + yoff2 - MESH_Y_SIZE/2));
			zbottom = min(zbottom, zval);
			ztop    = max(ztop, zval);
		}
	}
}

// should be named setup_lltex_sand_dirt() or something like that?
void init_terrain_mesh() {

	float const rel_wpz(get_rel_wpz());

	for (unsigned i = 0; i < NTEX_DIRT; ++i) { // move into gen_tex_height_tables()?
		lttex_dirt[i].id = mesh_tids_dirt[i];
		float const def_h(mesh_rh_dirt[i]);
		float h;

		if (mesh_rh_dirt[i] < W_PLANE_Z) { // below water
			h = def_h*rel_wpz/W_PLANE_Z;
		}
		else { // above water
			float const rel_h((def_h - W_PLANE_Z)/(1.0f - W_PLANE_Z));
			h = rel_wpz + rel_h*(1.0 - rel_wpz);
				
			if (mesh_tids_dirt[i] == SNOW_TEX) {
				h = min(h, def_h); // snow can't get lower when water lowers
				if (temperature > 40.0) h += 0.01*(temperature - 40.0); // less snow with increasing temperature
			}
		}
		lttex_dirt[i].zval = h;
	}
	gen_tex_height_tables();
}


void gen_terrain_map() {

	if (GLACIATE) {
		glaciate();
	}
	else {
		glaciate_exp     = 1.0;
		glaciate_exp_inv = 1.0;
	}
	apply_erosion(mesh_height[0], MESH_X_SIZE, MESH_Y_SIZE, zbottom, erosion_iters);
}


void estimate_zminmax(bool using_eq) {

	if (mesh_scale_change) {
		// don't include max() with -zmin or zmax when sine_bias has been set, as this creates sharp transitions in zmax_est
		// when the initial estimation doesn't include points that fall within the full range of sine values
		// TODO: better way to handle sine_bias? never max() with -zmin/zmax and just let height values fall outside the estimation?
		float zp(zmax_est);
		if (hmap_params.sine_bias >= 0.0) {zp = max(zp, -zmin);}
		if (hmap_params.sine_bias <= 0.0) {zp = max(zp,  zmax);}
		set_zmax_est(zp);
		mesh_scale_change = 0;
		set_zvals();
		return;
	}
	set_zmax_est(max(zmax, -zmin));

	if (using_eq && zmax == zmin) {
		set_zmax_est(zmax_est + 1.0E-6);
		return;
	}
	if (using_eq) {
		float const rm_scale(1000.0*XY_SCENE_SIZE/mesh_scale);
		mesh_xy_grid_cache_t height_gen;
		height_gen.build_arrays(0.0, 0.0, rm_scale, rm_scale, EST_RAND_PARAM, EST_RAND_PARAM);
		height_histogram.reserve(EST_RAND_PARAM*EST_RAND_PARAM/16); // 1024 values

		for (unsigned i = 0; i < EST_RAND_PARAM; ++i) {
			for (unsigned j = 0; j < EST_RAND_PARAM; ++j) {
				float const height(height_gen.eval_index(j, i)); // no glaciate
				zmax_est = max(zmax_est, float(fabs(height)));
				if (!(i&3) && !(j&3)) {height_histogram.push_back(height);} // only 1/16 of the values
			}
		}
		sort(height_histogram.begin(), height_histogram.end());
		if (mesh_gen_mode != MGEN_SINE) {zmax_est *= 1.2;} // perlin/simplex
	}
	set_zmax_est(1.1*zmax_est);
	set_zvals();
}

float get_median_height(float distribution_pos) {

	if (height_histogram.empty()) {return distribution_pos;} // ???
	return height_histogram[max(0, min((int)height_histogram.size()-1, int(height_histogram.size()*distribution_pos)))];
}


void set_zvals() {

	zcenter       = 0.5f*(zmax + zmin);
	zbottom       = zmin;
	ztop          = zmax;
	zmin          = -zmax_est;
	zmax          =  zmax_est;
	water_plane_z = get_water_z_height();
	CLOUD_CEILING = CLOUD_CEILING0*Z_SCENE_SIZE;
	LARGE_ZVAL    = max(LARGE_ZVAL, 100.0f*(CLOUD_CEILING + ztop));
}


float get_water_z_height() {

	float wpz(get_rel_wpz());
	if (GLACIATE) {wpz = do_glaciate_exp(wpz);}
	return wpz*zmax_est2 - zmax_est + water_h_off;
}


float get_cur_temperature() {return (combined_gu ? univ_temp : init_temperature);}


void update_temperature(bool verbose) {

	if (camera_mode != 1) return; // camera in air, don't use altitude temp

	// keep planet temperatures in combined landscape + universe
	alt_temp = get_cur_temperature();

	if (read_landscape || read_heightmap || do_read_mesh) {
		temperature = alt_temp;
		return;
	}
	float const cur_z((camera_mode == 1) ? get_camera_pos().z : 0.5f*(ztop + zbottom));
	float const rel_h((cur_z - zmin)/(zmax - zmin));
	//cout << "z: " << cur_z << ", rh: " << rel_h << ", wpz: " << water_plane_z << ", zmin: " << zmin << ", zmax: " << zmax << ", ztop: " << ztop << ", zbot: " << zbottom << endl;

	if (cur_z < water_plane_z) {
		float const znorm(min(1.0f, (water_plane_z - cur_z)/(water_plane_z - zmin)));
		alt_temp *= (1.0 - 0.9*znorm); // underwater
	}
	else if (rel_h > 0.2) {
		alt_temp -= 60.0*(rel_h - 0.2); // in snow covered mountains
	}
	if (verbose && !scrolling /*&& temperature != alt_temp*/) {cout << "Temperature = " << alt_temp << endl;}
	temperature = alt_temp;
}


void compute_scale() {

	int const iscale(int(log2(mesh_scale)));
	start_eval_sin = N_RAND_SIN2*max(0, min(NUM_FREQ_COMP-MIN_FREQS, (iscale+mesh_freq_filter)));
}

float get_hmap_scale(int mode) {
	float const scale((mode == MGEN_SIMPLEX || mode == MGEN_SIMPLEX_GPU || mode == MGEN_DWARP_GPU) ? 16.0 : 32.0); // simplex vs. perlin
	return scale*MESH_HEIGHT*mesh_height_scale*mesh_scale_z_inv;
}

void postproc_noise_zval(float &zval) {

	// Note: this applies to tree distributions as well, which is odd...
	hmap_params_t const &h(hmap_params);
	if (zval > h.plat_bot) {zval = h.plat_bot + h.plat_h*(zval - h.plat_bot) + min(h.plat_max, h.plat_s*(zval - h.plat_bot));} // plateau
	if (zval > h.crat_h  ) {zval = h.crat_h - h.crat_s*(zval - h.crat_h);} // craters
	if (zval > h.crack_lo && zval < h.crack_hi) {zval -= h.crack_d*min(zval-h.crack_lo, h.crack_hi-zval);} // cracks
}

void apply_noise_shape_final(float &noise, int shape) {
	switch (shape) {
	case 0: break; // linear - do nothing
	case 1: noise = fabs(noise) - 2.0; break; // billowy
	case 2: noise = 3.5 - fabs(noise); break; // ridged
	}
	postproc_noise_zval(noise);
}

void apply_noise_shape_per_term(float &noise, int shape) { // inputs and outputs [-1, 1]
	switch (shape) {
	case 0: break; // linear - do nothing
	case 1: noise = 2.0*fabs(noise) - 1.0; break; // billowy
	case 2: noise = 1.0 - 2.0*fabs(noise); break; // ridged
	}
}

void gen_rx_ry(float &rx, float &ry) {
	rand_gen_t rgen;
	apply_mesh_rand_seed(rgen);
	rx = rgen.rand_float() + 1.0;
	ry = rgen.rand_float() + 1.0;
}

bool mesh_xy_grid_cache_t::build_arrays(float x0, float y0, float dx, float dy,
	unsigned nx, unsigned ny, bool cache_values, bool force_sine_mode, bool no_wait)
{
	assert(nx > 0 && ny > 0);
	assert(start_eval_sin <= F_TABLE_SIZE);
	cur_nx = nx; cur_ny = ny; mx0 = x0; my0 = y0; mdx = dx; mdy = dy;
	gen_mode  = (force_sine_mode ? MGEN_SINE : mesh_gen_mode);
	gen_shape = (force_sine_mode ? 0 : mesh_gen_shape);
	do_glaciate = 0; // must set enable_glaciate() after this call if needed
	cached_vals.clear();

	if (gen_mode >= MGEN_SIMPLEX_GPU) { // GPU simplex noise - always cache values
		bool const is_running(cshader && cshader->get_is_running());
		if (!is_running) {run_gpu_simplex();} // launch the job
		if (no_wait && !is_running) return 0; // just started, results not yet available
		cache_gpu_simplex_vals();
		return 1; // results are available
	}
	yterms_start = nx*F_TABLE_SIZE;
	xyterms.resize((nx + ny)*F_TABLE_SIZE, 0.0);
	float const msx(mesh_scale*DX_VAL_INV), msy(mesh_scale*DY_VAL_INV), ms2(0.5*mesh_scale);

	for (int k = start_eval_sin; k < F_TABLE_SIZE; ++k) {
		float const x_mult(msx*sinTable[k][4]), y_mult(msy*sinTable[k][3]), y_scale(mesh_scale_z_inv*sinTable[k][0]);
		float const x_const(ms2*sinTable[k][4] + sinTable[k][2] + x_mult*x0), y_const(ms2*sinTable[k][3] + sinTable[k][1] + y_mult*y0);
		float const xmdx(x_mult*dx), ymdy(y_mult*dy);
		float *x_ptr(xyterms.data() + k), *y_ptr(x_ptr + yterms_start);

		for (unsigned i = 0; i < nx; ++i) {
			float sin_val(SINF(xmdx*i + x_const));
			//apply_noise_shape_per_term(sin_val, gen_shape);
			x_ptr[i*F_TABLE_SIZE] = sin_val;
		}
		for (unsigned i = 0; i < ny; ++i) {
			float sin_val(SINF(ymdy*i + y_const));
			//apply_noise_shape_per_term(sin_val, gen_shape);
			y_ptr[i*F_TABLE_SIZE] = y_scale*sin_val;
		}
	}
	if (cache_values) {
		cached_vals.resize(cur_nx*cur_ny);
		
#pragma omp parallel for schedule(static,1)
		for (int y = 0; y < (int)cur_ny; ++y) {
			for (unsigned x = 0; x < cur_nx; ++x) {
				cached_vals[y*cur_nx + x] = eval_index(x, y, 0, 0); // Note: no glaciate, min_start_sin=0, use_cache=0
			}
		}
	}
	return 1; // results are available
}

void mesh_xy_grid_cache_t::enable_glaciate() {

	do_glaciate = 1;
	if (hmap_params.sine_mag == 0.0) return;
	assert(cur_nx > 0 && cur_ny > 0); // build_arrays() must have been called first
	sine_mag_terms.resize(cur_nx*cur_ny);
	sine_offset = hmap_params.sine_bias*mesh_scale_z_inv;
	float const sm_scale(hmap_params.sine_mag*mesh_scale_z_inv), freq(mesh_scale*hmap_params.sine_freq);
	for (unsigned x = 0; x < cur_nx; ++x) {sine_mag_terms[x] = sm_scale*COSF((x*mdx + mx0)*DX_VAL_INV*freq);}
	for (unsigned y = 0; y < cur_ny; ++y) {sine_mag_terms[cur_nx + y] = COSF((y*mdy + my0)*DY_VAL_INV*freq);}
}

void mesh_xy_grid_cache_t::run_gpu_simplex() {

	//timer_t("GPU Mesh Gen");
	float const xy_scale(MESH_SCALE_FACTOR*mesh_scale), xscale(xy_scale*DX_VAL_INV), yscale(xy_scale*DY_VAL_INV);
	float rx, ry;
	gen_rx_ry(rx, ry);

	if (cshader == nullptr) {
		cshader = new grid_gen_shader_t("noise_2d_3d.part*+procedural_height_gen", cur_nx, cur_ny);
		if (gen_shape == 1) {cshader->set_comp_prefix("#define BILLOWY");}
		if (gen_shape == 2) {cshader->set_comp_prefix("#define RIDGED" );}
		if (gen_mode  == MGEN_DWARP_GPU) {cshader->set_comp_prefix("#define DOMAIN_WARP");}
		cshader->begin();
	}
	else {
		assert(cshader->get_xsize_req() == cur_nx && cshader->get_ysize_req() == cur_ny); // can't reuse with different sizes
		cshader->enable();
	}
	// returns heights in approximately [-1,1] range
	cshader->add_uniform_float("x0", (mx0 - 0.5*mdx)*xscale);
	cshader->add_uniform_float("y0", (my0 - 0.5*mdy)*yscale);
	cshader->add_uniform_float("dx", mdx*cur_nx*xscale);
	cshader->add_uniform_float("dy", mdy*cur_ny*yscale);
	cshader->add_uniform_float("rx", rx);
	cshader->add_uniform_float("ry", ry);
	cshader->add_uniform_float("zscale", /*zscale*/1.0);
	//cshader->gen_matrix_R32F(cached_vals, tid, 1, 1, 1); // reuse FBO
	cshader->setup_and_run(tid, 1, 1, 1); // R32F
	cshader->disable();
	cshader->prep_for_read_pixels(1);
}

void mesh_xy_grid_cache_t::cache_gpu_simplex_vals() {

	assert(cshader && cshader->get_is_running());
	cshader->read_float_vals(cached_vals, 1, 1); // reuse FBO
	float const zscale(get_hmap_scale(gen_mode));

	if (hmap_params.need_postproc()) {
		for (auto i = cached_vals.begin(); i != cached_vals.end(); ++i) {postproc_noise_zval(*i); *i *= zscale;}
	}
	else {
		for (auto i = cached_vals.begin(); i != cached_vals.end(); ++i) {*i *= zscale;}
	}
}

void mesh_xy_grid_cache_t::clear_context() { // for GPU-mode cached state
	free_texture(tid);
	if (cshader != nullptr) {cshader->end_shader(); free_cshader();}
}
void mesh_xy_grid_cache_t::free_cshader() { // don't make any GPU calls - hard reset, safe for destructors
	delete cshader; cshader = nullptr;
}


float gen_noise(float xv, float yv, int mode, int shape) {

	float zval(0.0), mag(1.0), freq(1.0), rx, ry;
	unsigned const end_octave(NUM_FREQ_COMP - start_eval_sin/N_RAND_SIN2);
	float const lacunarity(1.92), gain(0.5);
	gen_rx_ry(rx, ry);

	//#pragma omp parallel for schedule(static,1)
	for (unsigned i = 0; i < end_octave; ++i) {
		glm::vec2 const pos((freq*xv + rx), (freq*yv + ry));
		float noise((mode == MGEN_SIMPLEX || mode == MGEN_SIMPLEX_GPU || mode == MGEN_DWARP_GPU) ? glm::simplex(pos) : glm::perlin(pos));
		switch (shape) {
		case 0: break; // linear - do nothing
		case 1: noise = fabs(noise) - 0.40; break; // billowy
		case 2: noise = 0.45 - fabs(noise); break; // ridged
		//abs(0.5-abs(noise)*2.0)*2.0-0.5
		}
		zval += mag*noise;
		mag  *= gain;
		freq *= lacunarity;
		rx   *= 1.5;
		ry   *= 1.5;
	}
	return zval;
}

// mode: 0=sine tables, 1=simplex, 2=perlin, 3=GPU simplex, 4=GPU domain warp
// shape: 0=linear, 1=billowy, 2=ridged
float get_noise_zval(float xval, float yval, int mode, int shape) {

	assert(mode != MGEN_SINE); // mode 0 not supported by this function
	float const xy_scale(MESH_SCALE_FACTOR*mesh_scale);
	float xv(xy_scale*xval), yv(xy_scale*yval);

	if (mode == MGEN_DWARP_GPU) { // domain warping
		float const scale(0.2);
		float const dx1(gen_noise(xv+0.0, yv+0.0, mode, shape));
		float const dy1(gen_noise(xv+5.2, yv+1.3, mode, shape));
		float const dx2(gen_noise((xv + scale*dx1 + 1.7), (yv + scale*dy1 + 9.2), mode, shape));
		float const dy2(gen_noise((xv + scale*dx1 + 8.3), (yv + scale*dy1 + 2.8), mode, shape));
		xv += scale*dx2; yv += scale*dy2;
	}
	float zval(gen_noise(xv, yv, mode, shape));
	postproc_noise_zval(zval);
	return zval*get_hmap_scale(mode);
}


float mesh_xy_grid_cache_t::eval_index(unsigned x, unsigned y, int min_start_sin, bool use_cache) const {

	assert(x < cur_nx && y < cur_ny);
	float zval(0.0);

	if ((use_cache || gen_mode >= MGEN_SIMPLEX_GPU) && !cached_vals.empty()) {
		zval += cached_vals[y*cur_nx + x];
	}
	else if (gen_mode != MGEN_SINE) { // perlin/simplex
		float const xval((x*mdx + mx0)*DX_VAL_INV), yval((y*mdy + my0)*DY_VAL_INV);
		zval += get_noise_zval(xval, yval, gen_mode, gen_shape);
	}
	else { // sine tables
		float const *const xptr(xyterms.data() + x*F_TABLE_SIZE);
		float const *const yptr(xyterms.data() + yterms_start + y*F_TABLE_SIZE);
		int const start_ix(max(start_eval_sin, min_start_sin));
		// performance critical
		if (start_ix == 0) { // common case
			for (int i = 0; i < F_TABLE_SIZE; ++i) {zval += xptr[i]*yptr[i];}
		}
		else if (start_ix == 50) { // another common case
			for (int i = 50; i < F_TABLE_SIZE; ++i) {zval += xptr[i]*yptr[i];}
		}
		else {
			for (int i = start_ix; i < F_TABLE_SIZE; ++i) {zval += xptr[i]*yptr[i];}
		}
		apply_noise_shape_final(zval, gen_shape);
	}
	if (do_glaciate) {
		apply_glaciate(zval);
		
		if (hmap_params.sine_mag > 0.0) {
			assert(cur_nx + y < sine_mag_terms.size());
			zval += sine_mag_terms[x]*sine_mag_terms[cur_nx + y] + sine_offset;
			if (hmap_params.volcano_width > 0.0 && hmap_params.volcano_height > 0.0) {zval += get_volcano_height((x*mdx + mx0)*DX_VAL_INV, (y*mdy + my0)*DY_VAL_INV);}
		}
	}
	return zval;
}


// Note: called directly in tiled mesh and voxel code as a random number generator (not for mesh height);
// we always use sine tables here because get_noise_zval() is too slow
float eval_mesh_sin_terms(float xv, float yv) {

	float zval(0.0);

	for (int k = start_eval_sin; k < F_TABLE_SIZE; ++k) {
		float const *stk(sinTable[k]);
		zval += stk[0]*SINF(stk[3]*yv + stk[1])*SINF(stk[4]*xv + stk[2]); // performance critical
	}
	return zval;
}

float eval_mesh_sin_terms_scaled(float xval, float yval, float xy_scale) {

	float const xv(xy_scale*(xval - (MESH_X_SIZE >> 1))), yv(xy_scale*(yval - (MESH_Y_SIZE >> 1)));
	if (mesh_gen_mode != MGEN_SINE) {return get_noise_zval(xv, yv, mesh_gen_mode, mesh_gen_shape);}
	float val(eval_mesh_sin_terms(mesh_scale*xv, mesh_scale*yv)*mesh_scale_z_inv);
	apply_noise_shape_final(val, mesh_gen_shape);
	return val;
}


float get_exact_zval(float xval_in, float yval_in) {

	float xval((xval_in + X_SCENE_SIZE)*DX_VAL_INV + 0.5); // convert from real to index space, as in get_xpos()/get_ypos() but as FP
	float yval((yval_in + Y_SCENE_SIZE)*DY_VAL_INV + 0.5);

	if ((read_landscape || read_heightmap) && world_mode == WMODE_GROUND) { // interpolate from provided coords
		int xy[2] = {int(xval + 0.5), int(yval + 0.5)};
		clamp_to_mesh(xy);
		return mesh_height[xy[1]][xy[0]]; // could interpolate?
	}
	xval += xoff2; // offset by current mesh transform
	yval += yoff2;

	if (using_tiled_terrain_hmap_tex()) {
		float zval(get_tiled_terrain_height_tex(xval, yval));
		if (using_hmap_with_detail()) {zval += HMAP_DETAIL_MAG*eval_mesh_sin_terms_scaled(xval, yval, HMAP_DETAIL_SCALE);} // Note: agrees with tile_t::create_zvals()
		return zval;
	}
	if (world_mode == WMODE_INF_TERRAIN && mh_filename_tt != nullptr) {
		// using a heightmap texture, but it hasn't been loaded yet; this is where we get on the first frame of the loading screen
		// heightmap texture is 0-255, so use an intermediate value as the default/guess; this value works well for the heigtmap scene
		return (140*get_mh_texture_mult() + get_mh_texture_add());
	}
	float zval(eval_mesh_sin_terms_scaled(xval, yval, 1.0));
	apply_glaciate(zval);
	apply_mesh_sine(zval, (xval - float(MESH_X_SIZE >> 1)), (yval - float(MESH_Y_SIZE >> 1)));
	return zval;
}


void reset_offsets() {
	
	if (world_mode == WMODE_INF_TERRAIN) { // terrain
		point camera(get_camera_pos());
		xoff2 += int((get_xpos(camera.x) - MESH_X_SIZE/2));
		yoff2 += int((get_ypos(camera.y) - MESH_Y_SIZE/2));
	}
	else { // normal
		surface_pos = all_zeros;
	}
	xoff = yoff = 0;
}


void update_mesh(float dms, bool do_regen_trees) { // called when mesh_scale changes

	assert(dms > 0.0);
	++cache_counter;
	mesh_scale /= dms;
	tree_scale  = mesh_scale;
	xoff2       = int(xoff2*dms);
	yoff2       = int(yoff2*dms);
	set_zmax_est(zmax_est*pow(dms, (float)MESH_SCALE_Z_EXP));
	mesh_scale_z     = pow(mesh_scale, (float)MESH_SCALE_Z_EXP);
	mesh_scale_z_inv = 1.0/mesh_scale_z;

	if (world_mode == WMODE_INF_TERRAIN) {
		zmax = -LARGE_ZVAL;
		zmin =  LARGE_ZVAL;
		compute_scale();
		estimate_zminmax(1);
		gen_tex_height_tables();
	}
	else {
		gen_scene(1, do_regen_trees, 1, 1, 0);
		regen_lightmap();
	}
}


bool is_under_mesh(point const &p) {
	return ((world_mode == WMODE_GROUND && p.z < zbottom) || p.z < int_mesh_zval_pt_off(p, 0, 1, 1));
}


bool read_mesh(const char *filename, float zmm) {

	if (filename == NULL) {
		cout << "No input mesh file, using generated mesh." << endl;
		return 0;
	}
	FILE *fp;
	if (!open_file(fp, filename, "input mesh")) return 0;
	int xsize, ysize;
	float height;

	if (fscanf(fp, "%i%i", &xsize, &ysize) != 2) {
		cout << "Error reading size header in input file '" << filename << "'." << endl;
		checked_fclose(fp);
		return 0;
	}
	if (xsize != MESH_X_SIZE || ysize != MESH_Y_SIZE) {
		cout << "Error: Mesh size in file is " << xsize << "x" << ysize << " but should be " << MESH_X_SIZE << "x" << MESH_Y_SIZE << "." << endl;
		checked_fclose(fp);
		return 0;
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (fscanf(fp, "%f", &height) != 1) {
				cout << "Error reading mesh heights from file at position (" << j << ", " << i << ")." << endl;
				checked_fclose(fp);
				return 0;
			}
			mesh_height[i][j] = mesh_file_scale*height + mesh_file_tz;
		}
	}
	checked_fclose(fp);
	update_disabled_mesh_height();
	calc_zminmax();
	set_zmax_est((zmm != 0.0) ? zmm : max(-zmin, zmax));
	set_zvals();
	cout << "Read input mesh file '" << filename << "'." << endl;
	return 1;
}


bool write_mesh(const char *filename) {

	if (mesh_height == NULL) return 0;

	if (filename == NULL) {
		cout << "Error in filename: Cannot write mesh." << endl;
		return 0;
	}
	FILE *fp;
	if (!open_file(fp, filename, "mesh output", "w")) return 0;

	if (!fprintf(fp, "%i %i\n", MESH_X_SIZE, MESH_Y_SIZE)) {
		cout << "Error writing size header for mesh file '" << filename << "'." << endl;
		checked_fclose(fp);
		return 0;
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (!fprintf(fp, "%f ", mesh_height[i][j])) {
				cout << "Error writing mesh heights to file at position " << j << ", " << i << endl;
				checked_fclose(fp);
				return 0;
			}
		}
		fprintf(fp, "\n");
	}
	checked_fclose(fp);
	cout << "Wrote mesh file '" << filename << "'." << endl;
	return 1;
}


bool load_state(const char *filename) {

	if (filename == NULL) {
		cout << "Error in filename: Cannot load state." << endl;
		return 0;
	}
	FILE *fp;
	if (!open_file(fp, filename, "input state")) return 0;
	int v1, v2, v3, v4;
	++cache_counter;

	if (fscanf(fp, "%lf%lf%lf%f%f%f%f%f%f%i%i%i%i%i%li%li%i%i%i%i", &c_radius, &c_phi, &c_theta, &camera_origin.x,
		&camera_origin.y, &camera_origin.z, &surface_pos.x, &surface_pos.y, &surface_pos.z, &xoff, &yoff, &xoff2, &yoff2,
		&rand_gen_index, &global_rand_gen.rseed1, &global_rand_gen.rseed2, &v1, &v2, &v3, &v4) != 20)
	{
		cout << "Error reading state header." << endl;
		checked_fclose(fp);
		return 0;
	}
	if (v1 != MESH_X_SIZE || v2 != MESH_Y_SIZE || v3 != NUM_FREQ_COMP || v4 != N_RAND_SIN2) {
		cout << "Error: Saved state is incompatible with current configuration." << endl;
		checked_fclose(fp);
		return 0;
	}
	for (int i = 0; i < F_TABLE_SIZE; ++i) {
		for (int k = 0; k < 5; ++k) {
			if (fscanf(fp, "%f ", &(sinTable[i][k])) != 1) {
				cout << "Error reading state table entry (" << i << ", " << k << ")." << endl;
				checked_fclose(fp);
				return 0;
			}
		}
	}
	checked_fclose(fp);
	mesh_rgen_index = rand_gen_index;
	update_cpos();
	if (world_mode == WMODE_GROUND) {gen_scene(1, (world_mode == WMODE_GROUND), 1, 1, 0);}
	cout << "State file '" << filename << "' has been loaded." << endl;
	return 1;
}


bool save_state(const char *filename) {

	if (filename == NULL) {
		cout << "Error in filename: Cannot save state." << endl;
		return 0;
	}
	FILE *fp;
	if (!open_file(fp, filename, "output state", "w")) return 0;

	if (!fprintf(fp, "%lf %lf %lf %f %f %f %f %f %f %i %i %i %i %i %li %li\n%i %i %i %i\n",
		c_radius, c_phi, c_theta, camera_origin.x, camera_origin.y, camera_origin.z, surface_pos.x, surface_pos.y, surface_pos.z,
		xoff, yoff, xoff2, yoff2, rand_gen_index, global_rand_gen.rseed1, global_rand_gen.rseed2, MESH_X_SIZE, MESH_Y_SIZE, NUM_FREQ_COMP, N_RAND_SIN2))
	{
		cout << "Error writing state header." << endl;
		checked_fclose(fp);
		return 0;
	}
	for (unsigned i = 0; i < (unsigned)F_TABLE_SIZE; ++i) {
		for (unsigned k = 0; k < 5; ++k) {
			if (!fprintf(fp, "%f ", sinTable[i][k])) {
				cout << "Error writing state table entry (" << i << ", " << k << ")." << endl;
				checked_fclose(fp);
				return 0;
			}
		}
		fprintf(fp, "\n");
	}
	checked_fclose(fp);
	cout << "State file '" << filename << "' has been saved." << endl;
	return 1;
}


