// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari and Hiral Patel
// 3/15/02

#include "mesh.h"
#include "spillover.h"
#include "physics_objects.h"


float    const RIPPLE_DAMP1        = 0.95;
float    const RIPPLE_DAMP2        = 0.02;
float    const RIPPLE_MAT_ATTEN    = 0.965;
float    const MAX_RIPPLE_HEIGHT   = 1.0;
float    const MAX_SPLASH_SIZE     = 80.0;
float    const SMALL_DZ            = 0.001;
float    const WATER_WIND_EFF2     = 0.0006;
float    const MIN_SCALED_ALPHA    = 0.6;
float    const MAX_SCALED_ALPHA    = 1.1;
float    const SCALED_ALPHA_SLOPE  = 12.0;
float    const WLIGHT_SCALE        = 3.0; // 0.0 disables
float    const SNOW_THRESHOLD      = 8.0;
float    const MELT_RATE           = 10.0;
float    const SHADE_MELT          = 0.5;
float    const NIGHT_MELT          = 0.3;
float    const RAIN_VOLUME         = 0.05;
float    const BLOOD_VOLUME        = 0.01;
float    const FLOW_HEIGHT0        = 0.001;
float    const G_W_STOP_DEPTH      = 0.0;
float    const G_W_START_DEPTH     = 0.1;
float    const SQUARE_S_AREA       = 1.0;
float    const INIT_AREA           = 10.0;
float    const EVAPORATION         = 1.0E-5; // small for now
float    const EVAP_TEMP_POINT     = 0.0;
float    const W_TEX_STRETCH       = 2.0;
float    const MAX_WATER_ACC       = 25.0;
float    const DARK_WATER_ATTEN    = 1.0; // higher is lighter, lower is darker
float    const INT_WATER_ATTEN     = 0.5; // for interior water
unsigned const NUM_WATER_SPRINGS   = 2;
unsigned const MAX_RIPPLE_STEPS    = 2;
unsigned const UPDATE_STEP         = 8; // update water only every nth ripple computation
int      const EROSION_DIST        = 4;
float    const EROSION_RATE        = 0.0; //0.01
bool     const DEBUG_WATER_TIME    = 0; // DEBUGGING
bool     const DEBUG_RIPPLE_TIME   = 0;
bool     const FAST_WATER_RIPPLE   = 0;
bool     const NO_ICE_RIPPLES      = 0;
bool     const WATER_SHADOW        = 0; // uses mesh shadow, so is only accurate for shallow water
bool     const WATER_TEXTURE       = 1;
int      const UPDATE_UW_LANDSCAPE = 2;
int      const WATER_PRIM_TYPE     = GL_TRIANGLE_STRIP; // GL_QUAD_STRIP

float const w_spec[2][2] = {{0.3, 80.0}, {0.4, 70.0}};

enum {SPILL_NONE, SPILL_OUTSIDE, SPILL_INSIDE};



struct water_spring { // size = 40;

	int enabled;
	float dpf, diff, acc;
	point pos;
	vector3d vel;

	water_spring() {}

	water_spring(int en, float dpf_, float diff_, float acc_, point const &pos_, vector3d const &vel_)
		: enabled(en), dpf(dpf_), diff(diff_), acc(acc_), pos(pos_), vel(vel_) {}
};


struct water_section {

	float zval, wvol;
	int x1, y1, x2, y2;

	water_section() {}

	water_section(float x1_, float y1_, float x2_, float y2_, float zval_, float wvol_) : zval(zval_), wvol(wvol_),
		x1(get_xpos_clamp(x1_)), x2(get_xpos_clamp(x2_)), y1(get_ypos_clamp(y1_)), y2(get_ypos_clamp(y2_))
	{
		if (x1 > x2) swap(x1, x2);
		if (y1 > y2) swap(y1, y2);
	}
	bool is_nonzero() const {return (x1 < x2 && y1 < y2);}
};



// Global Variables
int total_watershed(0), w_acc(0), start_ripple(0), DISABLE_WATER(0), first_water_run(0), is_snow(0), added_wsprings(0);
float max_water_height, min_water_height, def_water_level;
vector<valley> valleys;
vector<water_spring> water_springs;
vector<water_section> wsections;
spillover spill;

extern bool mesh_invalidated, using_lightmap, has_dl_sources, has_snow;
extern int display_mode, frame_counter, game_mode, TIMESCALE2, I_TIMESCALE2, ocean_set;
extern int world_mode, island, rand_gen_index, begin_motion, animate, animate2, blood_spilled;
extern int landscape_changed, xoff2, yoff2, scrolling, dx_scroll, dy_scroll;
extern float temperature, zmax, zmin, zbottom, ztop, light_factor, water_plane_z, fticks, mesh_scale;
extern float TIMESTEP, TWO_XSS, TWO_YSS, XY_SCENE_SIZE;
extern point ocean;
extern vector3d up_norm, wind, total_wind;
extern int coll_id[];
extern obj_group obj_groups[];
extern vector<coll_obj> coll_objects;


void calc_water_normals();
void compute_ripples();
void update_valleys();
void update_water_volumes();
void draw_spillover(int i, int j, int si, int sj, int index, int vol_over, float blood_mix, float mud_mix);
int  calc_rest_pos(vector<int> &path_x, vector<int> &path_y, vector<char> &rp_set, int &x, int &y);
void calc_water_flow();
void init_water_springs(int nws);
void process_water_springs();
void add_waves();
void update_accumulation(int xpos, int ypos);
void shift_water_springs(vector3d const &vd);

void add_hole_in_landscape_texture(int xpos, int ypos, float blend);

void set_ocean_z();



inline bool cont_surf(int wsi, int wsi2, float zval) {

	assert(size_t(wsi) < valleys.size() && size_t(wsi2) < valleys.size());
	return (wsi == wsi2 || fabs(zval - valleys[wsi2].zval) < SMALL_DZ || spill.member2way(wsi, wsi2));
}


bool get_water_enabled(int x, int y) {

	if (water_enabled == NULL) return 1;
	int const xx(x + xoff2), yy(y + yoff2);
	return (point_outside_mesh(xx, yy) || water_enabled[yy][xx]);
}


bool has_water(int x, int y) {

	return (!point_outside_mesh(x, y) && wminside[y][x] && get_water_enabled(x, y));
}


bool mesh_is_underwater(int x, int y) {

	return (has_water(x, y) && mesh_height[y][x] < water_matrix[y][x]);
}


void select_water_ice_texture(colorRGBA &color, float *use_this_temp) {

	if (((use_this_temp != NULL) ? *use_this_temp : temperature) > W_FREEZE_POINT) { // water
		if (WATER_TEXTURE) select_texture(WATER_TEX); else glDisable(GL_TEXTURE_2D);
		set_specular(w_spec[0][0], w_spec[0][1]);
		color = WATER_C;
	}
	else { // ice
		if (WATER_TEXTURE) select_texture(ICE_TEX);   else glDisable(GL_TEXTURE_2D);
		set_specular(w_spec[1][0], w_spec[1][1]);
		color = ICE_C;
	}
	color *= DARK_WATER_ATTEN;
}


colorRGBA get_landscape_color(int xpos, int ypos) {

	assert(!point_outside_mesh(xpos, ypos));
	float const diffuse((shadow_mask[get_light()][ypos][xpos] & SHADOWED_ALL) ? 0.0 : get_cloud_shadow_atten(xpos, ypos));
	return get_landscape_texture_color(xpos, ypos)*(0.5*(1.0 + diffuse)); // half ambient and half diffuse
}


void get_object_color(int cindex, colorRGBA &color) {

	assert(unsigned(cindex) < coll_objects.size());
	coll_obj const &cobj(coll_objects[cindex]);

	// invisible simleys don't have a reflection
	if (cobj.cp.coll_func == smiley_collision && has_invisibility(cobj.cp.cf_index)) return;
	color = cobj.cp.color;
	if (cobj.cp.tid >= 0) color = color.modulate_with(texture_color(cobj.cp.tid));
}


class water_surface_draw {

	struct color_ix {
		colorRGBA c;
		int ix;
		color_ix() : ix(-1) {}
		color_ix(colorRGBA const &c_, int ix_) : c(c_), ix(ix_) {}
	};

	vector<color_ix> last_row_colors;
	bool big_water, inited;
	unsigned const bs, nx, ny;
	vector<unsigned char> underwater;

public:
	water_surface_draw() : big_water(0), inited(0), bs(5), nx(BITSHIFT_CEIL(MESH_X_SIZE, bs)), ny(BITSHIFT_CEIL(MESH_Y_SIZE, bs)) {
		last_row_colors.resize(MESH_X_SIZE);
	}

	void set_big_water() {big_water = 1;}

	void init() {

		if (inited) return;
		inited = 1;
		underwater.resize(nx*ny, 1);

		for (int y = 0; y < MESH_Y_SIZE; ++y) {
			for (int x = 0; x < MESH_X_SIZE; ++x) {
				if (max(v_collision_matrix[y][x].zmax, mesh_height[y][x]) > (water_matrix[y][x] - SMALL_NUMBER)) {
					unsigned const index((x >> bs) + nx*(y >> bs));
					assert(index < underwater.size());
					underwater[index] = 0; // not underwater
				}
			}
		}
	}

	void blend_reflection_color(point const &v, colorRGBA &color, vector3d const &n, point const &camera) {

		// Note: what about direct shadows on the water surface?
		assert(inited);
		vector3d const view_dir((v - camera).get_norm());
		
		if (dot_product(view_dir, n) >= 0.0) { // back facing water
			color = ALPHA0;
			return;
		}

		// add some green at shallow view angles
		blend_color(color, colorRGBA(0.0, 1.0, 0.5, 1.0), color, 0.25*(1.0 - fabs(dot_product(view_dir, n))), 0);

		vector3d dir;
		calc_reflection_angle(view_dir, dir, n);
		point vs(v), ve(v + dir*(4.0*XY_SCENE_SIZE/dir.mag()));
		bool skip_mesh(!do_line_clip_scene(vs, ve, zbottom, max(ztop, czmax)));
		point vs0(vs), ve0(ve);
		float const dval[2] = {DX_VAL, DY_VAL};
		int        uwv  [2] = {(get_xpos(vs.x) >> bs), (get_ypos(vs.y) >> bs)};
		int const  ve_xy[2] = {(get_xpos(ve.x) >> bs), (get_ypos(ve.y) >> bs)};

		// clip out regions of the ray that are entirely over the water
		// not entirely correct for dynamic objects since it uses czmin/czmax and skips regions over water
		while (!skip_mesh) {
			if (uwv[0] < 0 || uwv[0] >= int(nx) || uwv[1] < 0 || uwv[1] >= int(ny)) {skip_mesh = 1; break;} // all clipped
			if (!underwater[uwv[0] + nx*uwv[1]]) break;
			float tmin(1.0);
			unsigned dim0(0), dir0(0);

			for (unsigned dim = 0; dim < 2; ++dim) {
				float const delta(ve[dim] - vs[dim]);
				if (delta == 0.0) continue;
				bool const d(delta > 0.0);
				float const t((-SCENE_SIZE[dim] + dval[dim]*((uwv[dim] + d) << bs) - vs[dim])/delta);
				if (t < tmin) {tmin = t; dim0 = dim; dir0 = d;}
			}
			if (tmin == 1.0) break; // can this happen?
			if (uwv[0] == ve_xy[0] && uwv[1] == ve_xy[1]) {skip_mesh = 1; break;} // reached the end
			vs        += (ve - vs)*tmin;
			uwv[dim0] += (dir0 ? 1 : -1);
		} // while (1)
		int xpos, ypos;
		float mesh_zval;
		int const fast(2); // less accurate but faster
		colorRGBA rcolor(ALPHA0);
		bool mesh_int(0);

		if (!skip_mesh && line_intersect_mesh(vs, ve, xpos, ypos, mesh_zval, fast, 1) /*&& mesh_zval > vs.z*/) { // not sure about this last test
			float const t((mesh_zval - vs0.z)/(ve0.z - vs0.z));
			ve0      = (vs0 + (ve0 - vs0)*t);
			rcolor   = get_landscape_color(xpos, ypos);
			mesh_int = 1;
			
			if (dir.z < 0.0 && is_underwater(ve0)) {
				// if mesh int point is underwater, attenuate along light path and blend like we do when drawing the mesh
				// could calculate approx water intersection and call recursively, but there are too many problems with that approach
				water_color_atten_pt(&rcolor.red, xpos, ypos, ve0, v, get_light_pos());
				float const t2((water_matrix[ypos][xpos] - ve0.z)/(v.z - ve0.z));
				ve0 = (ve0 + (v - ve0)*t2); // updated approx water collision point
				blend_color(rcolor, color, rcolor, 0.5, 1); // add in a watery color
			}
		}
		int cindex;
		point cpos; // unused
		vector3d cnorm; // unused
		
		if (check_coll_line_exact(vs0, ve0, cpos, cnorm, cindex, 0.0, -1, 1, 0, 1)) { // skip_dynamic if !begin_motion?, !skip_mesh?
			get_object_color(cindex, rcolor);
			ve0 = cpos;
		}
		else if (!mesh_int) { // no mesh intersect and no cobj intersect
			vector3d const vdir((ve0 - vs0).get_norm());
			float const cloud_density(get_cloud_density(vs0, vdir)); // cloud_color vs. bkg_color
			blend_color(rcolor, get_cloud_color(), get_bkg_color(vs0, vdir), CLIP_TO_01(2.0f*cloud_density), 1);
		}
		if (begin_motion) { // find dynamic cobj intersection
			update_cobj_tree(1);
			if (check_coll_line_exact_tree(vs0, ve0, cpos, cnorm, cindex, -1, 1)) get_object_color(cindex, rcolor);
		}
		if (rcolor.alpha > 0.0) {
			float const r(get_fresnel_reflection(view_dir*-1, n, 1.0, WATER_INDEX_REFRACT));
			rcolor.alpha = 1.0; // transparent objects
			blend_color(color, rcolor, color, r*rcolor.alpha, 1);
		}
	}

	void draw_water_vertex(int i, int j, point const &v, colorRGBA color) { // i=y, j=x

		if (!sphere_in_camera_view(v, (DX_VAL + DY_VAL), 0)) return; // helps with reflections (is this correct?)
		vector3d const &n(wat_vert_normals[i][j]);
		assert(unsigned(j) < last_row_colors.size());
		
		if (last_row_colors[j].ix == i) { // gets here nearly half the time
			color = last_row_colors[j].c;
		}
		else {
			if (!(display_mode & 0x20) && !has_snow && v.z > mesh_height[i][j]) { // calculate water reflection and blend into color
				point const camera(get_camera_pos());
				if (camera.z > v.z) blend_reflection_color(v, color, n, camera); // below the camera
			}
			else if (WATER_SHADOW && shadow_mask[get_light()][i][j]) {
				color *= 0.5;
			}
			if (WLIGHT_SCALE > 0.0 && (using_lightmap || has_dl_sources)) {
				float const *const spec(w_spec[(temperature > W_FREEZE_POINT)]);
				get_sd_light(j, i, get_zpos(v.z), v, 0, WLIGHT_SCALE, (float *)&color, &n, spec);
			}
			last_row_colors[j] = color_ix(color, i);
		}
		color.do_glColor(); // note that the texture is blue
		n.do_glNormal();
		v.do_glVertex();
	}

	void draw_water_surface(int i, int j, colorRGBA const &color, int wsi) {

		float const x(get_xval(j)), y(get_yval(i));
		static float zval2(def_water_level);
		float const zval1(water_matrix[i][j] - SMALL_NUMBER); // what if wminside[i][j]==0 ?

		if (i < MESH_Y_SIZE-1) {
			if (big_water || watershed_matrix[i+1][j].wsi == wsi ||
				(wminside[i+1][j] == 1 && cont_surf(wsi, watershed_matrix[i+1][j].wsi, valleys[wsi].zval)))
			{
				zval2 = water_matrix[i+1][j] - SMALL_NUMBER;
			}
		}
		for (unsigned d = 0; d < 2; ++d) {
			draw_water_vertex(i+d, j, point(x, y+d*DY_VAL, (d ? zval2 : zval1)), color);
		}
	}
};


void draw_water() {

	RESET_TIME;
	int wsi(0), last_water(2), last_draw(0), lc0(landscape_changed);
	colorRGBA color(WHITE);
	static int wcounter(0);
	static float tdx(0.0), tdy(0.0);
	float const tx_scale(W_TEX_STRETCH/TWO_XSS), ty_scale(W_TEX_STRETCH/TWO_YSS);
	process_water_springs();
	add_waves();
	if (DISABLE_WATER || (island && !w_acc)) return;
	water_surface_draw wsd;
	set_fill_mode();
	point const camera(get_camera_pos());
	bool const is_ice(temperature <= W_FREEZE_POINT);

	if (!is_ice) {
		tdx -= WATER_WIND_EFF2*wind.x*fticks;
		tdy -= WATER_WIND_EFF2*wind.y*fticks;
	}
	float const tx_val(tx_scale*xoff2*DX_VAL + tdx), ty_val(ty_scale*yoff2*DX_VAL + tdy);

	if (!island) { // draw exterior water (oceans)
		if (camera.z >= water_plane_z) draw_water_sides(1);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_NORMALIZE);
		enable_blend();
		enable_point_specular();
		glNormal3f(0.0, 0.0, 1.0); // probably unnecessary
		select_water_ice_texture(color);
		setup_texgen(tx_scale, ty_scale, tx_val, ty_val);
		color.alpha *= 0.5;
		wsd.init();
		wsd.set_big_water();
		//glCullFace(GL_FRONT); // backwards?
		//glEnable(GL_CULL_FACE);

		for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
			for (int j = 0; j < MESH_X_SIZE; ++j) {
				if (j < MESH_X_SIZE-1 && (
					(wminside[i][j]     == 2 && water_matrix[i][j]     > mesh_height[i][j])   ||
					(wminside[i+1][j]   == 2 && water_matrix[i+1][j]   > mesh_height[i+1][j]) ||
					(wminside[i][j+1]   == 2 && water_matrix[i][j+1]   > mesh_height[i][j+1]) ||
					(wminside[i+1][j+1] == 2 && water_matrix[i+1][j+1] > mesh_height[i+1][j+1])))
				{
					if (!last_draw) glBegin(WATER_PRIM_TYPE);
					wsd.draw_water_surface(i, j, color, 0);
					last_draw = 1;
				}
				else if (last_draw) {
					wsd.draw_water_surface(i, j, color, 0);
					glEnd();
					last_draw = 0;
				}
			}
		}
		assert(!last_draw);
		//glDisable(GL_CULL_FACE);
		//glCullFace(GL_BACK); // backwards?
		disable_blend();
		disable_point_specular();
		set_specular(0.0, 1.0);
		disable_textures_texgen();
		glEnable(GL_NORMALIZE);
		glDisable(GL_COLOR_MATERIAL);
		last_draw = 0;
		if (DEBUG_WATER_TIME) {PRINT_TIME("Water Draw Fixed");}
		if (camera.z <  water_plane_z) draw_water_sides(1);
	}
	wsd.init();
	if (DEBUG_WATER_TIME) {PRINT_TIME("Water Init");}

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (wminside[i][j] == 1 && (rand()&255) == 0 && get_water_enabled(j, i) && !is_mesh_disabled(j, i)) {
				modify_grass_at(point(get_xval(j), get_yval(i), mesh_height[i][j]), HALF_DXY, 0, 0, 0, 0, 1);
			}
		}
	}
	if (DEBUG_WATER_TIME) {PRINT_TIME("Grass Update");}

	glEnable(GL_COLOR_MATERIAL);
	update_valleys();
	if (DEBUG_WATER_TIME) {PRINT_TIME("Water Valleys Update");}

	// call the function that computes the ripple effect
	assert(fticks != 0.0);
	//int const time_mod(int(((float)TIMESCALE2)/max(1.0f, fticks)));
	int const time_mod(1);

	if (animate2 && (time_mod <= 1 || (wcounter%time_mod) == 0)) {
		unsigned num_steps;

		if (FAST_WATER_RIPPLE || !(start_ripple || first_water_run) || is_ice || (!start_ripple && first_water_run)) {
			num_steps = 1;
		}
		else {
			assert(I_TIMESCALE2 > 0);
			num_steps = max(1U, min(MAX_RIPPLE_STEPS, unsigned(min(3.0f, fticks)*min(MAX_I_TIMESCALE, I_TIMESCALE2))));
		}
		for (unsigned i = 0; i < num_steps; ++i) {
			compute_ripples();
		}
		calc_water_normals();
	}
	if (DEBUG_WATER_TIME) {PRINT_TIME("Water Ripple Update");}
	glEnable(GL_TEXTURE_2D);
	setup_texgen(tx_scale, ty_scale, tx_val, ty_val);
	glDisable(GL_NORMALIZE);
	enable_blend();
	enable_point_specular();
	unsigned nin(0);
	int xin[4], yin[4], last_wsi(-1);
	bool const disp_snow((display_mode & 0x40) && temperature <= SNOW_MAX_TEMP);
	select_water_ice_texture(color);
	color *= INT_WATER_ATTEN; // attenuate for interior water
	colorRGBA wcolor(color);
	glBegin(WATER_PRIM_TYPE);
	
	// draw interior water (ponds)
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (!is_ice && is_snow) update_accumulation(i, j);
			nin = 0;

			if (wminside[i][j] == 1) {
				wsi = watershed_matrix[i][j].wsi;
				assert(size_t(wsi) < valleys.size());
				float const z_min(z_min_matrix[i][j]), zval(valleys[wsi].zval);

				if (zval >= z_min - G_W_STOP_DEPTH && (!island || zval > ocean.z)) {
					if (!is_ice && UPDATE_UW_LANDSCAPE && (rand()&63) == 0) {
						add_hole_in_landscape_texture(j, i, 1.2*fticks*(0.02 + zval - z_min));
					}
					float delta_area(1.0);

					if (point_interior_to_mesh(j, i)) { // more accurate but less efficient
						delta_area = 0.0;

						for (int y = -1; y <= 1; ++y) {
							for (int x = -1; x <= 1; ++x) {
								if (zval >= mesh_height[i+y][j+x]) delta_area += 1.0;
							}
						}
						delta_area /= 9.0;
					}
					valleys[wsi].area += SQUARE_S_AREA*delta_area;
					
					if (zval >= z_min && i < MESH_Y_SIZE-1 && j < MESH_X_SIZE-1) {
						xin[nin] = j; yin[nin] = i; ++nin;

						for (int ii = i; ii < i+2; ++ii) {
							for (int jj = j + (ii == i); jj < j+2; ++jj) {
								if (wminside[ii][jj] == 1 && cont_surf(wsi, watershed_matrix[ii][jj].wsi, zval)) {
									xin[nin] = jj; yin[nin] = ii; ++nin;
								}
							}
						}
						if (nin == 4) {
							if (disp_snow && accumulation_matrix[i][j] >= SNOW_THRESHOLD &&
								(mesh_height[i][j] <= zval || (i > 0 && j > 0 && mesh_height[i-1][j-1] <= zval)))
							{
								if (last_water != 0) { // snow
									glEnd();
									select_texture(SNOW_TEX);
									set_specular(0.6, 20.0);
									glBegin(WATER_PRIM_TYPE);
									last_water = 0;
									color      = WHITE;
								}
							}
							else if (last_water != 1) {
								glEnd();
								select_water_ice_texture(color);
								color *= INT_WATER_ATTEN; // attenuate for interior water
								glBegin(WATER_PRIM_TYPE);
								last_water = 1;
							}
							if (wsi != last_wsi) {
								float const blood_mix(valleys[wsi].blood_mix), mud_mix(valleys[wsi].mud_mix);
								wcolor = color;

								if (blood_mix > 0.99) { // all blood
									wcolor = RED;
								}
								else {
									wcolor.alpha *= 0.5;
									if (mud_mix   > 0.0) blend_color(wcolor, (!is_ice ? MUD_S_C : LT_BROWN), wcolor, mud_mix, 1);
									if (blood_mix > 0.0) blend_color(wcolor, RED, wcolor, blood_mix, 1);
								}
								last_wsi = wsi;
							}
							wsd.draw_water_surface(i, j, wcolor, wsi);
							last_draw = 1;
						} // nin == 4
					}
				}
			}
			if (nin < 4) {
				if (last_draw) wsd.draw_water_surface(i, j, wcolor, wsi);
				if (last_draw || nin == 3) glEnd();

				if (nin == 3) { // single triangle, color is already set from above draw
					glBegin(GL_TRIANGLES);

					for (unsigned p = 0; p < 3; ++p) {
						wat_vert_normals[yin[p]][xin[p]].do_glNormal();
						glVertex3f(get_xval(xin[p]), get_yval(yin[p]), (water_matrix[yin[p]][xin[p]] - SMALL_NUMBER));
					}
					glEnd();
				}
				if (last_draw || nin == 3) glBegin(WATER_PRIM_TYPE);
				last_draw = 0;
			}
		} // for j
	} // for i
	glEnd();
	disable_blend();
	disable_point_specular();
	set_specular(0.0, 1.0);
	disable_textures_texgen();
	glEnable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	update_water_volumes();
	if (!lc0 && rand()%5 != 0) landscape_changed = 0; // reset, only update landscape 20% of the time
	++wcounter;
	first_water_run = 0;
	if (DEBUG_WATER_TIME) {PRINT_TIME("Water Draw");}
}


void calc_water_normals() {

	if (DISABLE_WATER) return;
	vector3d *wsn0(wat_surf_normals[0]), *wsn1(wat_surf_normals[1]);

	if (temperature <= W_FREEZE_POINT) {
		for (int i = 0; i < MESH_Y_SIZE; ++i) {
			for (int j = 0; j < MESH_X_SIZE; ++j) {
				if (!point_interior_to_mesh(j, i) || wminside[i][j]) wsn1[j] = wat_vert_normals[i][j] = plus_z;
			} // for j
			swap(wsn0, wsn1);
		} // for i
		return;
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (!point_interior_to_mesh(j, i)) {
				wsn1[j] = wat_vert_normals[i][j] = plus_z;
				continue;
			}
			if (wminside[i][j] && water_matrix[i][j] >= z_min_matrix[i][j]) {
				vector3d nv(get_matrix_surf_norm(water_matrix, NULL, MESH_X_SIZE, MESH_Y_SIZE, i, j));
				wsn1[j] = nv;
				
				if (i > 0) {
					nv += wsn0[j];
					if (j > 0) nv += wsn0[j-1];
				}
				if (j > 0)  nv    += wsn1[j-1];
				if (i == 0) nv.z  += 2.0;
				if (j == 0) nv.z  += 2.0;
				wat_vert_normals[i][j] = nv*0.25;
			} // inside
		} // for j
		swap(wsn0, wsn1);
	} // for i
}


inline void update_water_edges(int i, int j) {

	if (i > 0 && wminside[i-1][j] == 1) {
		water_matrix[i][j] = valleys[watershed_matrix[i-1][j].wsi].zval;
	}
	else if (j > 0 && wminside[i][j-1] == 1) {
		water_matrix[i][j] = valleys[watershed_matrix[i][j-1].wsi].zval;
	}
	else {
		water_matrix[i][j] = def_water_level;
	}
}


void compute_ripples() {

	if (DISABLE_WATER) return;
	static unsigned dtime1(0), dtime2(0), counter(0);
	bool const update_iter((counter%UPDATE_STEP) == 0);
	RESET_TIME;

	if (temperature > W_FREEZE_POINT && (start_ripple || first_water_run)) {
		float const rm_atten(pow(RIPPLE_MAT_ATTEN, fticks)), rdamp1(pow(RIPPLE_DAMP1, fticks)), rdamp2(RIPPLE_DAMP2*fticks);
		start_ripple = 0;

		for (int i = 0; i < MESH_Y_SIZE; ++i) {
			for (int j = 0; j < MESH_X_SIZE; ++j) {
				if (!wminside[i][j] || water_matrix[i][j] < z_min_matrix[i][j]) continue;
				short const i8(watershed_matrix[i][j].inside8);
				fix_fp_mag(ripples[i][j].rval);
				float const rmij(ripples[i][j].rval);
				float &acc(ripples[i][j].acc);
				fix_fp_mag(acc);
				acc *= rm_atten;
				if (!start_ripple && fabs(acc) > 1.0E-6) start_ripple = 1;

				// 00 0- -0 0+ +0 -- +- ++ -+  22  11
				// 01 02 04 08 10 20 40 80 100 200 400
				if (point_interior_to_mesh(j, i)) { // fast mode
					float const d0( rmij - ripples[i  ][j-1].rval);
					if (i8 & 0x02)  ripples[i  ][j-1].acc += d0;
					float const d2( rmij - ripples[i  ][j+1].rval);
					if (i8 & 0x08)  ripples[i  ][j+1].acc += d2;
					float const d4((rmij - ripples[i-1][j-1].rval)*SQRTOFTWOINV);
					if (i8 & 0x20)  ripples[i-1][j-1].acc += d4;
					float const d1( rmij - ripples[i-1][j  ].rval);
					if (i8 & 0x04)  ripples[i-1][j  ].acc += d1;
					float const d7((rmij - ripples[i-1][j+1].rval)*SQRTOFTWOINV);
					if (i8 & 0x100) ripples[i-1][j+1].acc += d7;
					float const d5((rmij - ripples[i+1][j-1].rval)*SQRTOFTWOINV);
					if (i8 & 0x40)  ripples[i+1][j-1].acc += d5;
					float const d3( rmij - ripples[i+1][j  ].rval);
					if (i8 & 0x10)  ripples[i+1][j  ].acc += d3;
					float const d6((rmij - ripples[i+1][j+1].rval)*SQRTOFTWOINV);
					if (i8 & 0x80)  ripples[i+1][j+1].acc += d6;
					acc -= d0 + d1 + d2 + d3 + d4 + d5 + d6 + d7;
					fix_fp_mag(acc);
					continue;
				}
				if (j > 0) {
					float const dz(rmij - ripples[i][j-1].rval);
					if (i8 & 0x02) ripples[i][j-1].acc += dz;
					acc -= dz;

					if (i > 0) {
						float const dz((rmij - ripples[i-1][j-1].rval)*SQRTOFTWOINV);
						if (i8 & 0x20) ripples[i-1][j-1].acc += dz;
						acc -= dz;
					}
					if (i < MESH_Y_SIZE-1) {
						float const dz((rmij - ripples[i+1][j-1].rval)*SQRTOFTWOINV);
						if (i8 & 0x40) ripples[i+1][j-1].acc += dz;
						acc -= dz;
					}
				}
				if (i > 0) {
					float const dz(rmij - ripples[i-1][j].rval);
					if (i8 & 0x04) ripples[i-1][j].acc += dz;
					acc -= dz;
				}
				if (j < MESH_X_SIZE-1) {
					float const dz(rmij - ripples[i][j+1].rval);
					if (i8 & 0x08) ripples[i][j+1].acc += dz;
					acc -= dz;

					if (i < MESH_Y_SIZE-1) {
						float const dz((rmij - ripples[i+1][j+1].rval)*SQRTOFTWOINV);
						if (i8 & 0x80) ripples[i+1][j+1].acc += dz;
						acc -= dz;
					}
					if (i > 0) {
						float const dz((rmij - ripples[i-1][j+1].rval)*SQRTOFTWOINV);
						if (i8 & 0x100) ripples[i-1][j+1].acc += dz;
						acc -= dz;
					}
				}
				if (i < MESH_Y_SIZE-1) {
					float const dz(rmij - ripples[i+1][j].rval);
					if (i8 & 0x10) ripples[i+1][j].acc += dz;
					acc -= dz;
				}
				fix_fp_mag(acc);
			} // for j
		} // for i
		if (DEBUG_RIPPLE_TIME) dtime1 += GET_DELTA_TIME;
		
		for (int i = 0; i < MESH_Y_SIZE; ++i) {
			for (int j = 0; j < MESH_X_SIZE; ++j) {
				float ripple_zval(0.0);

				if (wminside[i][j]) {
					float const zval(rdamp1*(ripples[i][j].rval + rdamp2*ripples[i][j].acc)); // ripple wave height
					ripple_zval = ((fabs(zval) < TOLERANCE) ? 0.0 : zval); // prevent small floating point numbers
				}
				if (wminside[i][j] == 1) { // dynamic water
					int const wsi(watershed_matrix[i][j].wsi);
					assert(size_t(wsi) < valleys.size());

					if (water_matrix[i][j] < z_min_matrix[i][j] && fabs(ripples[i][j].rval) < 1.0E-4 && fabs(ripples[i][j].acc) < 1.0E-4) { // under ground - no ripple
						if (update_iter) water_matrix[i][j] = valleys[wsi].zval;
						continue;
					}
					float const depth(valleys[wsi].depth);

					if (depth < 0) {
						ripples[i][j].rval *= rm_atten;
						if (update_iter) water_matrix[i][j] = valleys[wsi].zval;
						continue;
					}
					float const zval(max(min(ripple_zval, depth), -depth)); // max ripple height equals water depth
					ripples[i][j].rval = rm_atten*zval;
					water_matrix[i][j] = valleys[wsi].zval + zval;
				}
				else if (wminside[i][j] == 2) { // fixed water
					ripples[i][j].rval = rm_atten*ripple_zval;
					water_matrix[i][j] = water_plane_z + min(MAX_RIPPLE_HEIGHT, ripple_zval);
					water_matrix[i][j] = max(water_matrix[i][j], zbottom);
				}
				else if (update_iter && get_water_enabled(j, i)) {
					update_water_edges(i, j);
				}
			} // for j
		} // for i
		if (DEBUG_RIPPLE_TIME) dtime2 += GET_DELTA_TIME;
	}
	else { // no ripple
		matrix_clear_2d(ripples);

		// must clear ripples at least once at the beginning
		if (NO_ICE_RIPPLES || counter == 0 || temperature > W_FREEZE_POINT) {
			for (int i = 0; i < MESH_Y_SIZE; ++i) {
				for (int j = 0; j < MESH_X_SIZE; ++j) {
					if (wminside[i][j] == 1) {
						water_matrix[i][j] = valleys[watershed_matrix[i][j].wsi].zval;
					}
					else {
						update_water_edges(i, j);
					}
				} // for j
			} // for i
		}
	} // ripple
	++counter;

	if (DEBUG_RIPPLE_TIME && (counter%20) == 0) {
		cout << "times = " << dtime1 << ", " << dtime2 << endl; // cumulative
		dtime1 = dtime2 = 0;
	}
}


void add_splash(int xpos, int ypos, float energy, float radius) {

	//energy *= 10.0;
	if (DISABLE_WATER || !(display_mode & 0x04))  return;
	if (temperature <= W_FREEZE_POINT && !island) return;
	int in_wmatrix(0), wsi(0);
	float water_mix(1.0), blood_mix(0.0), mud_mix(0.0);

	if (!point_outside_mesh(xpos, ypos) && wminside[ypos][xpos] == 1 && !is_mesh_disabled(xpos, ypos)) {
		if (water_matrix[ypos][xpos] < (mesh_height[ypos][xpos] - 0.5*radius)) return; // water is too low
		float const mh_r(mesh_height[ypos][xpos] + radius);

		if (h_collision_matrix[ypos][xpos] < mh_r && water_matrix[ypos][xpos] < (mh_r + MAX_SPLASH_DEPTH)) {
			wsi        = watershed_matrix[ypos][xpos].wsi;
			blood_mix  = valleys[wsi].blood_mix;
			mud_mix    = valleys[wsi].mud_mix;
			water_mix  = (1.0 - blood_mix);
			in_wmatrix = 1;
		}
	}
	assert(energy >= 0.0);
	energy = min(energy, 10000.0f);
	float const splash_size(min(0.006f*(0.5f + water_mix)*sqrt((4.0f/XY_SCENE_SIZE)*energy), MAX_SPLASH_SIZE));
	if (splash_size < TOLERANCE) assert(0); //return;
	float const rad(min(4, (int)ceil(0.5*radius/(DX_VAL + DY_VAL)))), radsq(rad*rad);
	int const droplet_id(coll_id[DROPLET]);

	if (energy > 2.0 && droplet_id >= 0 && temperature > W_FREEZE_POINT) {
		// splash size = 0.3 - rain, 2.0 - ball, 60.0 - rocket
		float const sqrt_energy(sqrt(energy));
		int ndrops(min(MAX_SPLASH_DROP, int((0.5 + 0.2*water_mix)*sqrt_energy)));

		if (ndrops > 0 && point_interior_to_mesh(xpos, ypos)) {
			if (water_matrix[ypos][xpos] >= z_min_matrix[ypos][xpos] && wminside[ypos][xpos]) {
				if (in_wmatrix) {
					valleys[wsi].mud_mix += 0.12*sqrt_energy/(valleys[wsi].w_volume + 1.0);
					valleys[wsi].mud_mix  = mud_mix = min(1.0f, valleys[wsi].mud_mix);
				}
				if (wminside[ypos][xpos] == 1) { // dynamic water
					ndrops = min(ndrops, int(0.25*valleys[watershed_matrix[ypos][xpos].wsi].w_volume));
				}
				if (ndrops > 0) {
					obj_group &objg(obj_groups[droplet_id]);
					//valleys[watershed_matrix[ypos][xpos].wsi].w_volume -= ndrops;
					float const vz(5.0 + (0.2 + 0.1*water_mix)*sqrt_energy);
					point const pos(get_xval(xpos), get_yval(ypos), (water_matrix[ypos][xpos] + 2.5*radius));
					
					for (int o = 0; o < ndrops; ++o) { // less efficient
						int const i(objg.choose_object());
						objg.create_object_at(i, pos);
						objg.get_obj(i).velocity = gen_rand_vector(vz*rand_uniform(0.05, 0.1), 20.0, PI_TWO);
						vadd_rand(objg.get_obj(i).pos, 1.0*radius);
					}
				}
			}
		}
	}
	int const irad((int)rad), x1(max((xpos - irad), 0)), y1(max((ypos - irad), 0));
	int const x2(min((xpos + irad), MESH_X_SIZE-1)), y2(min((ypos + irad), MESH_Y_SIZE-1));

	for (int i = y1; i <= y2; i++) {
		int const di_sq((i - ypos)*(i - ypos));

		for (int j = x1; j <= x2; j++) {
			if ((di_sq + (j - xpos)*(j - ypos)) <= radsq && wminside[i][j]) ripples[i][j].rval += splash_size;
		}
	}
	start_ripple = 1;
}


void add_waves() { // add waves due to wind

	//RESET_TIME;
	if (DISABLE_WATER || !(display_mode & 0x04) || !(display_mode & 0x0100)) return;
	if (temperature <= W_FREEZE_POINT || !animate2) return;
	float const wave_amplitude(0.006/mesh_scale), wave_freq(0.06), depth_scale(20.0);
	float const wind_freq(0.02), wind_amplitude(0.15/mesh_scale), wind_resist(0.75);
	float const wxoff(wind_resist*total_wind.x/TWO_XSS/MESH_X_SIZE);
	float const wyoff(wind_resist*total_wind.y/TWO_YSS/MESH_Y_SIZE);
	float const fticks_clamped(min(fticks, 10.0f));
	static float wave_time(0.0);
	wave_time += fticks_clamped;
	
	for (int y = 0; y < MESH_Y_SIZE; ++y) {
		for (int x = 0; x < MESH_X_SIZE; ++x) {
			if (!wminside[y][x]) continue; // only in water
			float const depth(water_matrix[y][x] - mesh_height[y][x]);
			if (depth < SMALL_NUMBER) continue; // not deep enough for waves
			point const p(get_xval(x), get_yval(y), water_matrix[y][x]);
			vector3d const local_wind(get_local_wind(p));
			float const lwmag(local_wind.mag());
			float const tx(min(0.2f, fabs(local_wind.y))*wind_freq*(x + xoff2)/lwmag - wxoff);
			float const ty(min(0.2f, fabs(local_wind.x))*wind_freq*(y + yoff2)/lwmag - wyoff);
			float const val(get_texture_alpha(CLOUD_RAW_TEX, tx, ty));
			float wval(wind_amplitude*min(2.5f, sqrt(lwmag))*val*min(depth, 0.1f));
			if (wminside[y][x] == 2) wval += wave_amplitude*fticks_clamped*sin(wave_freq*wave_time + depth_scale*depth); // outside water (oceans)
			ripples[y][x].rval += wval;
			start_ripple = 1;
		}
	}
	//PRINT_TIME("Add Waves");
}


// *** BEGIN VALLEYS/SPILLOVER ***


void check_spillover(int i, int j, int ii, int jj, int si, int sj, float zval, int wsi) { // does wsi overflow?

	float const z_over(zval - mesh_height[ii][jj]);
	if (z_over <= max(0.0f, valleys[wsi].sf.z_over)) return;

	if (wminside[i][j] != 1) { // spillover outside
		valleys[wsi].sf = valley::spill_func(-1, i, j, si, sj, SPILL_OUTSIDE, z_over);
		return;
	}
	int const index(watershed_matrix[i][j].wsi); // spills into this pool
	assert(size_t(index) < valleys.size());

	if (index != wsi && zval > valleys[index].zval) { // spillover inside
		if (!spill.member(index, wsi)) valleys[wsi].sf = valley::spill_func(index, i, j, si, sj, SPILL_INSIDE, z_over); // pools not combined
		valleys[wsi].spill_index = index;
	}
}


void sync_water_height(int wsi, int skip_ix, float zval, float z_over) {

	vector<unsigned> cc;
	spill.get_connected_components(wsi, cc);

	for (unsigned k = 0; k < cc.size(); ++k) {
		if (cc[k] == wsi || cc[k] == skip_ix) continue;
		valley &v(valleys[cc[k]]);
		v.zval      = zval;//max(zval, v.zval);
		v.dz       -= z_over;
		v.w_volume  = v.lwv; // ???
		v.mud_mix   = valleys[wsi].mud_mix; // is this necessary?
		v.blood_mix = valleys[wsi].blood_mix;
	}
}


void update_valleys() {

	vector<unsigned> cc;

	for (unsigned i = 0; i < valleys.size(); ++i) {
		valley &v(valleys[i]);
		assert(v.area >= 0.0);
		assert(v.fvol >= 0.0);
		assert(v.w_volume >= 0.0);
		v.area     += INIT_AREA*SQUARE_S_AREA;
		v.w_volume += v.fvol;
		v.fvol      = 0.0;
		v.blood_mix = CLIP_TO_01(v.blood_mix);
		v.mud_mix   = ((v.mud_mix < 0.0001) ? 0.0 : CLIP_TO_01(v.mud_mix)*pow(0.998f, fticks)); // slowly settles
		float const dv(min((v.w_volume - v.lwv), MAX_WATER_ACC));
		float delta_z(v.get_volume()*dv/v.area); // this is inaccurate because v.area is not constant
		//if (v.old_area > 0.0 && v.area > 0.0 && v.area != v.old_area) delta_z -= (1.0 - v.old_area/(0.5*(v.area + v.old_area)))*v.dz;
		v.dz    = delta_z;
		v.zval += v.dz;
		v.spill_integral = 0.9*v.spill_integral + 0.1*v.spill_vol;
		if (v.w_volume == 0.0) spill.remove_connected(i);
	}
	vector<bool> used(valleys.size(), 0);

	for (unsigned i = 0; i < valleys.size(); ++i) {
		if (used[i]) continue; // already processed
		spill.get_connected_components(i, cc);
		
		if (cc.empty()) { // unconnected (optimization)
			used[i] = 1;
			continue;
		}
		cc.push_back(i); // add current v to the combined set
		valley combined;

		for (vector<unsigned>::const_iterator j = cc.begin(); j != cc.end(); ++j) { // calculate sums
			assert(*j < valleys.size());
			assert(!used[*j]);
			used[*j] = 1;
			valley &v(valleys[*j]);
			combined.area      += v.area;
			combined.w_volume  += v.w_volume;
			combined.zval      += v.zval*v.area; // zval determined by surface area (weighed average)
			combined.blood_mix += v.w_volume*v.blood_mix; // composition determined by volume
			combined.mud_mix   += v.w_volume*v.mud_mix;
		}
		assert(combined.area     > 0.0);
		assert(combined.w_volume > 0.0);
		combined.zval     /= combined.area;
		combined.blood_mix = CLIP_TO_01(combined.blood_mix/combined.w_volume); // do we have to clip these?
		combined.mud_mix   = CLIP_TO_01(combined.mud_mix  /combined.w_volume);

		for (vector<unsigned>::const_iterator j = cc.begin(); j != cc.end(); ++j) { // update all connected valleys
			valley &v(valleys[*j]);
			float const delta_z(combined.zval - v.zval); // new - old
			v.w_volume  = max(TOLERANCE, (v.w_volume + delta_z*v.area/v.get_volume())); // not sure about this???
			v.dz       += delta_z;
			v.zval      = combined.zval; // do we want to update v.w_volume based on v.dz?
			v.blood_mix = combined.blood_mix;
			v.mud_mix   = combined.mud_mix;
		}
	}
	assert(used.size() == valleys.size());

	for (unsigned i = 0; i < valleys.size(); ++i) { // reset for next iteration
		valley &v(valleys[i]);
		v.area        = 0.0;
		v.spill_vol   = 0.0;
		v.has_spilled = 0;
		v.sf.spill    = 0;
		v.sf.z_over   = 0.0;
		v.spill_index = -1;
		v.depth       = v.zval - mesh_height[v.y][v.x];
	}

	// check for spillover offscreen or into another pool
	int const ijd[4][4] = {{0,1,0,1}, {0,-1,0,0}, {1,0,1,0}, {-1,0,0,0}};

	for (int i = 1; i < MESH_Y_SIZE-1; ++i) {
		for (int j = 1; j < MESH_X_SIZE-1; ++j) {
			if (wminside[i][j] != 1) continue;
			int const wsi(watershed_matrix[i][j].wsi);
			float const zval(valleys[wsi].zval);
			if (zval < z_min_matrix[i][j]) continue;

			for (unsigned k = 0; k < 4; ++k) {
				check_spillover(i+ijd[k][0], j+ijd[k][1], i+ijd[k][2], j+ijd[k][3], i, j, zval, wsi);
			}
		}
	}
	for (unsigned i = 0; i < valleys.size(); ++i) {
		valley &v(valleys[i]); // pool that may be spilling (source)
		valley::spill_func const &sf(v.sf);

		if (sf.spill != SPILL_NONE) { // this pool has spilled over its border
			if (sf.z_over <= 0.0) continue;
			// what happens when w_volume < lwv? what about dz < 0?
			float vol_over(max(0.0f, (v.w_volume - v.lwv)));
			if (v.dz > TOLERANCE) vol_over *= min(1.0f, (sf.z_over/v.dz)); // fraction of last amount that spilled
			
			if (sf.spill == SPILL_INSIDE && vol_over > 0.0) { // spilled into another pool (1 direction), send water to the target
				assert(watershed_matrix[sf.si][sf.sj].wsi == i);
				assert(sf.index != i);
				valley &vs(valleys[sf.index]); // pool that is spilled into (dest)
				float const nv_ratio(CLIP_TO_01(vol_over/max(TOLERANCE, vs.w_volume)));
				v.dz         -= sf.z_over;
				vs.spill_vol += vol_over;
				vs.blood_mix  = CLIP_TO_01(nv_ratio*v.blood_mix + (1.0f - nv_ratio)*vs.blood_mix);
				vs.mud_mix    = CLIP_TO_01(nv_ratio*v.mud_mix   + (1.0f - nv_ratio)*vs.mud_mix);
				spill.insert(i, sf.index); // source, dest
			}
			v.w_volume    = max(0.0f, (v.w_volume - vol_over)); // ???
			v.has_spilled = 1;
			float const zval(max((v.zval - sf.z_over), v.min_zval)); // zval at spill point (local minima)
			sync_water_height(i, sf.index, zval, sf.z_over);
			draw_spillover(sf.i, sf.j, sf.si, sf.sj, sf.index, int(vol_over), v.blood_mix, v.mud_mix);
			v.zval = zval;
		}
		else if (v.spill_index >= 0) {
			spill.insert(i, v.spill_index); // source, dest
		}
		/*else if (wminside[sf.si][sf.sj] == 1 && watershed_matrix[sf.si][sf.sj].wsi == i) { // current sf is valid
			if (v.zval + TOLERANCE < mesh_height[sf.si][sf.sj]) spill.remove_all_i(i); // water is below min point => no spill
		}*/
		max_water_height = max(max_water_height, valleys[i].zval);
	}
}


void update_water_volumes() {

	for (unsigned i = 0; i < valleys.size(); ++i) {
		valley &v(valleys[i]);
		bool const has_lwv(v.lwv > 0.0);

		if (v.has_spilled) {
			v.w_volume = v.lwv; // spilled into another pool, reset to last value
		}
		else {
			v.lwv += min((v.w_volume - v.lwv), MAX_WATER_ACC); // update (with thresh)
		}
		v.w_volume += v.spill_vol; // new water added from other pools this frame

		if (temperature > EVAP_TEMP_POINT && v.w_volume > 0.0) { // apply evaporation
			v.w_volume = max(0.0f, (v.w_volume - EVAPORATION*(temperature - EVAP_TEMP_POINT)*v.area));
		}
		if (v.w_volume < 0.0) cout << "wv: " << v.w_volume << ", sv: " << v.spill_vol << ", lwv: " << v.lwv << endl; // testing
		assert(v.w_volume >= 0.0);
		
		// don't let the graph get unlinked until the depth is 0 (should use z_spill or something?)
		//if (v.w_volume == 0.0 && v.depth > 0.0 && v.area > 0.0) v.w_volume = v.lwv = TOLERANCE;
		
		if (has_lwv && v.w_volume == 0.0) { // no more water (all evaporated)
			if (v.zval < v.min_zval) {
				v.zval = v.min_zval;
			}
			else {
				v.lwv = 1.0;
			}
			v.has_spilled = 0;
			spill.remove_all_i(i); //spill.remove_all_i_2way(i); // reset connectivity graph entry
		}
	}
}


// *** END VALLEYS/SPILLOVER ***


int draw_spill_section(int x1, int y1, int x2, int y2, float z1, float z2, float width, int volume, int index, float blood_mix, float mud_mix) {

	assert(abs(x2 - x1) <= 1 && abs(y2 - y1) <= 1);
	float const flow_height(FLOW_HEIGHT0*Z_SCENE_SIZE);
	assert(!point_outside_mesh(x1, y1));
	if (!get_water_enabled(x1, y1) && !get_water_enabled(x2, y2)) return 2;
	spillway_matrix[y1][x1] = (short)frame_counter;
	float xa(get_xval(x1)), ya(get_yval(y1));

	if (world_mode != WMODE_INF_TERRAIN && (z1+flow_height) < water_plane_z && (z2+flow_height) < water_plane_z) {
		return 0;
	}
	if (island && (z1+0.02) < ocean.z && (z2+0.02) < ocean.z) return 0;
	vector3d const &norm(vertex_normals[y1][x1]);
	norm.do_glNormal();

	if (x1 == x2 && y1 == y2) {
		for (unsigned i = 0; i < 2; ++i) glVertex3f(xa, ya, z1); // end at a point
		return 0;
	}
	float const xb(get_xval(x2)), yb(get_yval(y2));
	float const slope(fabs((z2 - z1)/sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya))));

	if (UPDATE_UW_LANDSCAPE > 1 && (rand()%5) == 0) {
		int xv[4] = {x1, x2, x2, x1}, yv[4] = {y1, y1, y2, y2};
		float const rate(0.01*fticks*slope*sqrt((float)volume));

		if (rate > 0) {
			for (unsigned i = 0; i < 4; ++i) {
				add_hole_in_landscape_texture(xv[i], yv[i], rate);

				if (EROSION_RATE > 0.0) { // erode terrain
					update_mesh_height(xv[i], yv[i], EROSION_DIST, EROSION_RATE*rate, 4.0, 1);
					// update grass?
				}
			}
		}
	}
	colorRGBA color(BLOOD_C);

	if (blood_mix < 0.99) {
		if (index < 0 && (min(x1, x2) <= 1 || min(y1, y2) <= 1 || max(x1, x2) >= MESH_X_SIZE-2 || max(y1, y2) >= MESH_Y_SIZE-2)) {
			color.assign(0.3, 0.3, 1.0, 0.22);
		}
		else if (z1 == z2) {
			color.assign(0.3, 0.3, 1.0, 0.2);
		}
		else {
			assert(slope != 0.0);
			color.blue  = min(1.0f, max(0.6f, (0.3f + 0.6f/slope)));
			color.alpha = min(0.6f, max(0.2f, (0.2f + 0.2f*slope)));
			color.red   = color.green = min(0.7*color.blue, max(0.3, (0.2 + 0.4*slope)));
		}
		if (mud_mix   > 0.01) blend_color(color, MUD_C,   color, mud_mix,   1);
		if (blood_mix > 0.01) blend_color(color, BLOOD_C, color, blood_mix, 1);
	}
	color.do_glColor();
	xa += flow_height*norm.x;
	ya += flow_height*norm.y;
	z1 += flow_height*norm.z;

	if (x1 != x2 && y1 != y2) { // diagonal
		assert(abs(x2 - x1) <= 1 && abs(y2 - y1) <= 1);
		width = width/SQRT2;
	}
	for (unsigned i = 0; i < 2; ++i) {
		float const xv(xa + ((y1 == y2) ? 0.0 : ((i ^ (y1 > y2)) ?   width : -width)));
		float const yv(ya + ((x1 == x2) ? 0.0 : ((i ^ (x1 > x2)) ?  -width :  width)));
		glVertex3f(xv, yv, z1);
	}
	return 1;
}


void draw_spillover(int i, int j, int si, int sj, int index, int vol_over, float blood_mix, float mud_mix) {

	if (vol_over <= 0) return;
	assert(!point_outside_mesh(j, i));
	bool const nov(index == -1);
	assert(nov || size_t(index) < valleys.size());
	int x1(j), x2(j), y1(i), y2(i);
	bool last_iteration(0);
	float z1(mesh_height[i][j]);
	float const zval(nov ? zmin : valleys[index].zval);
	float const width(min(0.012, 0.1*(DX_VAL + DY_VAL)*(sqrt((float)vol_over) + 1.0)));
	float const v_splash(min(10.0f, (float)vol_over));
	int const xs(nov ? -1 : valleys[index].x), ys(nov ? -1 : valleys[index].y);
	int count(0);
	enable_blend();
	glBegin(GL_QUAD_STRIP);
	up_norm.do_glNormal();

	if (!point_outside_mesh(sj, si)) {
		draw_spill_section(sj, si, x1, y1, mesh_height[si][sj], z1, width, vol_over, index, blood_mix, mud_mix);
	}
	for (; (count < XY_SUM_SIZE) && (x1 != xs || y1 != ys) && point_interior_to_mesh(x1, y1); ++count) {
		x2 = w_motion_matrix[y1][x1].x;
		y2 = w_motion_matrix[y1][x1].y;
		assert(!point_outside_mesh(x2, y2));
		float const z2(mesh_height[y2][x2]);
		int const draw_res(draw_spill_section(x1, y1, x2, y2, z1, z2, width, vol_over, index, blood_mix, mud_mix));
		if (last_iteration || draw_res == 0) add_splash(x1, y1, 1.5*v_splash, 0.002*v_splash); // hit fixed ocean/lake
		if (last_iteration || draw_res != 1 || (x2 == x1 && y2 == y1)) break; // edge, disabled, or valley
		last_iteration = (zval >= z2);
		z1 = (last_iteration ? zval : z2);
		x1 = x2;
		y1 = y2;
	}
	if (count == XY_SUM_SIZE) cout << "Error: Exceeded iteration limit in draw_spillover()." << endl;
	glEnd();
	disable_blend();
}


void float_downstream(point &pos, float radius) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos) || wminside[ypos][xpos] != 1) return; // not in an inside pool
	if (!mesh_is_underwater(xpos, ypos)) return; // no longer floating?
	if ((pos.z - radius) < (mesh_height[ypos][xpos] + 0.5*radius))   return; // too close to the bottom of the pool
	int const wsi(watershed_matrix[ypos][xpos].wsi);
	assert(size_t(wsi) < valleys.size());
	valley::spill_func const &sf(valleys[wsi].sf);
	if (sf.spill == SPILL_NONE || sf.z_over == 0.0 || valleys[wsi].spill_integral == 0.0) return; // not spilled
	point const spill_pt(get_xval(sf.sj), get_yval(sf.si), pos.z);
	float const dist(p2p_dist(pos, spill_pt));
	if (dist < SMALL_NUMBER) return; // avoid divide-by-zeros
	float const vel(5.0E-5*valleys[wsi].spill_integral/(DX_VAL + DY_VAL + dist));
	pos += (spill_pt - pos)*(min(0.005f, vel)/dist);
}


void calc_watershed() {

	int mode(0);

	if (DISABLE_WATER == 1) {
		for (int i = 0; i < MESH_Y_SIZE; ++i) {
			for (int j = 0; j < MESH_X_SIZE; ++j) {
				water_matrix[i][j] = zmin - 1.0;
			}
		}
		if (world_mode == WMODE_GROUND) def_water_level = water_plane_z = zmin; // wm == 3???
		init_water_springs(NUM_WATER_SPRINGS);
		return;
	}
	if (island) {
		set_ocean_z();
		def_water_level = water_plane_z = ocean.z;
		mode = 1;
	}
	else {
		ocean.z = water_plane_z;

		if (ztop < water_plane_z) { // all water
			def_water_level = water_plane_z;

			for (int i = 0; i < MESH_Y_SIZE; ++i) {
				for (int j = 0; j < MESH_X_SIZE; ++j) {
					water_matrix[i][j]     = water_plane_z;
					wminside[i][j]         = 2;
					wat_vert_normals[i][j] = plus_z;
				}
			}
			max_water_height = def_water_level;
			min_water_height = def_water_level;
			return;
		}
		else if (zbottom < water_plane_z) { // some water
			def_water_level = water_plane_z;
			mode            = 1;
		}
		else { // no water
			def_water_level = zmin;
		}
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			watershed_matrix[i][j].x = watershed_matrix[i][j].y = 0;
		}
	}
	max_water_height = def_water_level;
	min_water_height = def_water_level;
	vector<char> rp_set(XY_MULT_SIZE, 0);
	vector<int> path_x(XY_SUM_SIZE), path_y(XY_SUM_SIZE);

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (!get_water_enabled(j, i)) { // disabled
				wminside[i][j] = 0;
				continue;
			}
			int x(j), y(i);
			int const crp(point_interior_to_mesh(j, i) ? calc_rest_pos(path_x, path_y, rp_set, x, y) : 0);
			wminside[i][j] = ((mode == 1 && mesh_height[y][x] < water_plane_z) ? 2 : crp);
		}
	}
	calc_water_flow();
	init_water_springs(NUM_WATER_SPRINGS);
	matrix_clear_2d(ripples);
	first_water_run = 1;

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (wminside[i][j] == 1) { // dynamic water
				int const wsi(watershed_matrix[i][j].wsi);
				assert(wsi >= 0 && wsi < (int)valleys.size());
				water_matrix[i][j] = valleys[wsi].zval;
			}
			else if (wminside[i][j] == 2) { // fixed water
				water_matrix[i][j] = water_plane_z;
			}
			else { // no water
				water_matrix[i][j] = def_water_level; // this seems safe
			}
			short &i8(watershed_matrix[i][j].inside8);
			i8 = 0;
			// 00 0- -0 0+ +0 -- +- ++ -+  22  11
			// 01 02 04 08 10 20 40 80 100 200 400
			if (wminside[i][j])      i8 |= 0x01;
			if (wminside[i][j] == 2) i8 |= 0x200;
			if (wminside[i][j] == 1) i8 |= 0x400;
			
			if (j > 0 && wminside[i][j-1]) {
				i8 |= 0x02;
				if (i > 0 && wminside[i-1][j-1]) i8 |= 0x20;
			}
			if (i > 0 && wminside[i-1][j]) {
				i8 |= 0x04;
				if (j < MESH_X_SIZE-1 && wminside[i-1][j+1]) i8 |= 0x100;
			}
			if (j < MESH_X_SIZE-1 && wminside[i][j+1]) {
				i8 |= 0x08;
				if (i < MESH_Y_SIZE-1 && wminside[i+1][j+1]) i8 |= 0x80;
			}
			if (i < MESH_Y_SIZE-1 && wminside[i+1][j]) {
				i8 |= 0x10;
				if (j > 0 && wminside[i+1][j-1]) i8 |= 0x40;
			}
		}
	}
}


int calc_rest_pos(vector<int> &path_x, vector<int> &path_y, vector<char> &rp_set, int &x, int &y) { // return 0 if off the map

	int path_counter(0), x2(0), y2(0);
	bool found(0);

	do {
		if (path_counter < XY_SUM_SIZE) {
			path_x[path_counter] = x;
			path_y[path_counter] = y;
			++path_counter;
		}
		char const rps(rp_set[y*MESH_X_SIZE+x]);

		if (rps) {
			x2    = watershed_matrix[y][x].x;
			y2    = watershed_matrix[y][x].y;
			x     = x2;
			y     = y2;
			found = (rps == 1);
			break;
		}
		x2 = w_motion_matrix[y][x].x;
		y2 = w_motion_matrix[y][x].y;

		if (x2 == x && y2 == y) {
			found = 1;
			break;
		}
		x = x2;
		y = y2;
	} while (point_interior_to_mesh(x, y));

	for (--path_counter; path_counter >= 0; --path_counter) {
		int const x2(path_x[path_counter]), y2(path_y[path_counter]);
		watershed_matrix[y2][x2].x = x;
		watershed_matrix[y2][x2].y = y;
		rp_set[y2*MESH_X_SIZE+x2]  = (1 + (found == 0));
	}
	return found;
}


bool add_water_section(float x1, float y1, float x2, float y2, float zval, float wvol) {

	water_section ws(x1, y1, x2, y2, zval, wvol);
	if (!ws.is_nonzero()) return 0;
	wsections.push_back(ws);
	return 1;
}


void calc_water_flow() {

	if (DISABLE_WATER == 1) return;
	vector<int> minima;
	map<pair<int, int>, valley> pool_zvals;

	if (scrolling) { // save old pool zval state
		for (unsigned i = 0; i < valleys.size(); ++i) {
			pool_zvals.insert(make_pair(make_pair((valleys[i].x - dx_scroll), (valleys[i].y - dy_scroll)), valleys[i]));
		}
	}
	valleys.clear();
	spill.clear();

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (wminside[i][j] != 1) continue;
			int const wmij(watershed_matrix[i][j].x + MESH_X_SIZE*watershed_matrix[i][j].y);
			if (minima.empty() || wmij != minima.back()) minima.push_back(wmij);
		}
	}
	unsigned const msize(minima.size());
	if (msize == 0) return;
	std::sort(minima.begin(), minima.end());
	total_watershed = 0;

	for (unsigned i = 0; i < wsections.size(); ++i) {
		valleys.push_back(valley(((wsections[i].x1 + wsections[i].x2) >> 1), ((wsections[i].y1 + wsections[i].y2) >> 1)));
	}
	for (unsigned i = 0; i < msize; ++i) {
		if (i == msize-1 || minima[i] != minima[i+1]) {
			valley const v((minima[i]%MESH_X_SIZE), (minima[i]/MESH_X_SIZE));
			if (!(mesh_height[v.y][v.x] <= water_plane_z || !get_water_enabled(v.x, v.y))) valleys.push_back(v);
		}
	}
	if (valleys.size() > 1000) cout << "Warning: This landscape contains " << valleys.size() << " water pooling locations." << endl;
	
	if (valleys.size() > 32767) {
		cout << "Error: Too many water pools. Max is 32767." << endl;
		exit(1);
	}
	spill.init(valleys.size());
	vector<int> vmap(XY_MULT_SIZE, 0);

	for (unsigned k = 0; k < valleys.size(); ++k) {
		int const index(valleys[k].x + MESH_X_SIZE*valleys[k].y);
		assert(index < XY_MULT_SIZE);
		vmap[index] = k+1;
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			if (wminside[i][j] == 1) {
				int const index(watershed_matrix[i][j].x + MESH_X_SIZE*watershed_matrix[i][j].y);
				assert(index < XY_MULT_SIZE);

				if (vmap[index] == 0) {
					wminside[i][j] = 0;
				}
				else {
					watershed_matrix[i][j].wsi = vmap[index]-1;
				}
				++total_watershed;
			}
			else {
				watershed_matrix[i][j].wsi = -1; // invalid
			}
		}
	}
	for (unsigned i = 0; i < wsections.size(); ++i) {
		int const x1(max(0, wsections[i].x1)), y1(max(0, wsections[i].y1));
		int const x2(min(MESH_X_SIZE-1, wsections[i].x2)), y2(min(MESH_Y_SIZE-1, wsections[i].y2));
		assert(x1 <= x2 && y1 <= y2);

		for (int y = y1; y < y2; ++y) {
			for (int x = x1; x < x2; ++x) {
				wminside[y][x]             = 1;
				watershed_matrix[y][x].wsi = i;

				if (water_enabled != NULL) { // not sure about this
					if (y > y1 && y < y2-1 && x > x1 && x < x2-1) water_enabled[y][x] = 1;
				}
			}
		}
	}
	for (unsigned i = 0; i < valleys.size(); ++i) {
		valleys[i].create(i);
			
		if (scrolling) { // attempt to restore zval state
			map<pair<int, int>, valley>::const_iterator it(pool_zvals.find(make_pair(valleys[i].x, valleys[i].y)));
			
			if (it != pool_zvals.end()) { // not quite right - lacks spill/merge dependencies and other per-frame state, pool might be cut off by the mesh
				valleys[i].copy_state_from(it->second);
			}
		}
	}
	assert(valleys.size() >= wsections.size());

	for (unsigned i = 0; i < wsections.size(); ++i) {
		valleys[i].zval = valleys[i].min_zval = wsections[i].zval;
		valleys[i].fvol = wsections[i].wvol;
		valleys[i].lwv  = valleys[i].fvol;
	}
}


void init_water_springs(int nws) {

	if (added_wsprings || nws == 0) return;
	water_springs.clear();
	unsigned const smod(max(1, (XY_MULT_SIZE/nws)));

	for (int i = 1; i < MESH_Y_SIZE-1; ++i) {
		for (int j = 1; j < MESH_X_SIZE-1; ++j) {
			rseed1 = 54563  *(i + yoff2) + 23423  *rand_gen_index; // not prime
			rseed2 = 4365435*(j + xoff2) + 6456541*rand_gen_index; // not prime
			if ((rand2()%smod) != 0)   continue;
			point const pos(get_xval(j), get_yval(i), (interpolate_mesh_zval(get_xval(j), get_yval(i), 0.0, 0, 1) + 0.02));
			if (pos.z < water_plane_z) continue;
			water_springs.push_back(water_spring(1, rand_uniform2(1.5, 3.0), 0.1, 0.0, pos, gen_rand_vector2(5.0, 3.0, PI_TWO)));
		}
	}
}


void process_water_springs() {

	//if (DISABLE_WATER) return;
	int const cid(coll_id[WDROPLET]);
	if (world_mode != WMODE_GROUND || !begin_motion || !animate || !animate2 || temperature <= W_FREEZE_POINT) return;
	if (ztop < water_plane_z)  return; // all water
	if (water_springs.empty()) return;

	for (unsigned i = 0; i < water_springs.size(); ++i) {
		water_spring &ws(water_springs[i]);
		if (!ws.enabled) continue;
		
		if (temperature > get_max_t(WDROPLET)) { // water boils
			if (rand()&1) gen_smoke(ws.pos);
			continue;
		}
		ws.acc += fticks*(TIMESTEP/DEF_TIMESTEP)*ws.dpf;

		for (unsigned j = 0; j < ws.acc; ++j) {
			int const k(obj_groups[cid].choose_object());
			obj_groups[cid].create_object_at(k, ws.pos);
			vadd_rand(obj_groups[cid].get_obj(k).pos, 0.02*ws.diff);
			obj_groups[cid].get_obj(k).velocity = ws.vel + gen_rand_vector(ws.vel.mag()*ws.diff, 1.0, PI);
		}
		ws.acc -= (int)ws.acc;
	}
}


void add_water_spring(point const &pos, vector3d const &vel, float rate, float diff, int calc_z, int gen_vel) {

	vector3d const vel0(gen_vel ? gen_rand_vector(5.0, 3.0, PI_TWO) : vel);
	water_spring ws(1, rate, diff, 0.0, pos, vel0);
	if (calc_z) ws.pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1) + 0.02;
	water_springs.push_back(ws);
	++added_wsprings;
}


void shift_water_springs(vector3d const &vd) {

	if (added_wsprings == 0) return; // dynamically created, not placed

	for (unsigned i = 0; i < water_springs.size(); ++i) {
		water_springs[i].pos += vd;
	}
}


// update region is inclusive: [x1,x2]x[y1,y2]
void update_water_zvals(int x1, int y1, int x2, int y2) {

	for (int i = y1; i <= y2; ++i) {
		for (int j = x1; j <= x2; ++j) {
			// *** WRITE *** - add new pool if there is a new local mimima?
		}
	}
}


bool is_underwater(point const &pos, int check_bottom, float *depth) { // or under ice
	
	if (depth) *depth = 0.0;
	if (DISABLE_WATER) return 0; // water disabled

	if (world_mode == WMODE_INF_TERRAIN) {
		if (depth) *depth = max(0.0f, (water_plane_z - pos.z));
		return (pos.z < water_plane_z);
	}
	if (ocean_set && pos.z < ocean.z) {
		if (depth) *depth = ocean.z - pos.z;
		return 1;
	}
	assert(water_matrix && mesh_height);
	if (!(display_mode & 0x04) || !is_over_mesh(pos) || pos.z < zmin || pos.z > max_water_height) return 0;
	//if (pos.z < min_water_height) return 1;
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (!has_water(xpos, ypos)) return 0; // off the mesh or disabled
	
	if (pos.z <= water_matrix[ypos][xpos] && (!check_bottom || pos.z >= mesh_height[ypos][xpos])) {
		if (depth) *depth = water_matrix[ypos][xpos] - pos.z;
		return 1;
	}
	return 0;
}


inline void update_accumulation(int xpos, int ypos) {

	float acc(accumulation_matrix[ypos][xpos]);
	if (acc <= 0.0) return;
	w_acc  = 1; // melt snow
	float melted(((temperature - W_FREEZE_POINT)/MELT_RATE)*(NIGHT_MELT + (1.0 - NIGHT_MELT)*light_factor));
	if (light_factor >= 0.5 && (shadow_mask[LIGHT_SUN][ypos][xpos] & SHADOWED_ALL)) melted *= SHADE_MELT;
	accumulation_matrix[ypos][xpos] = max(0.0f, (acc - melted));

	if (!DISABLE_WATER && wminside[ypos][xpos] == 1) {
		int const wsi(watershed_matrix[ypos][xpos].wsi);
		valleys[wsi].fvol = min(4.0f*MAX_WATER_ACC, valleys[wsi].fvol + SNOW_ACC*melted);
	}
}


void select_liquid_color(colorRGBA &color, int xpos, int ypos) {

	if (point_outside_mesh(xpos, ypos) || wminside[ypos][xpos] != 1 || (island && mesh_height[ypos][xpos] <= ocean.z)) return;
	int const wsi(watershed_matrix[ypos][xpos].wsi);
	assert(size_t(wsi) < valleys.size());
	float const blood_mix(valleys[wsi].blood_mix);
	float const mud_mix(valleys[wsi].mud_mix);
	blend_color(color, MUD_C, color, mud_mix, 1);
	blend_color(color, BLOOD_C, color, blood_mix, 1);
}


void select_liquid_color(colorRGBA &color, point const &pos) {

	select_liquid_color(color, get_xpos(pos.x), get_ypos(pos.y));
}


void valley::copy_state_from(valley const &v) {

	zval      = v.zval;
	w_volume  = v.w_volume;
	blood_mix = v.blood_mix;
	mud_mix   = v.mud_mix;
	fvol      = v.fvol;
	dz        = v.dz;
	lwv       = v.lwv;
}


void valley::create(int wsi) {

	assert(!point_outside_mesh(x, y));
	zval = min_zval = max(def_water_level, (mesh_height[y][x] - (float)G_W_START_DEPTH));
	int pindex(watershed_matrix[y][x].wsi), px(x); // initialize start to center

	if (size_t(pindex) >= wsections.size() && pindex != wsi) {
		cout << "Error in valley::create(): Invalid pool index: " << wsi << " vs. " << pindex << "." << endl;
		assert(0);
	}
	if (!wminside[y][x]) {
		cout << "Error in valley::create(): Zero flag in watershed matrix." << endl;
		assert(0);
	}
	for (++px; pindex == wsi && px < MESH_X_SIZE && wminside[y][px]; ++px) { // find first point on edge in +x direction (not used yet)
		pindex = watershed_matrix[y][px].wsi;
	}
}


float valley::get_volume() const {

	return (blood_mix*BLOOD_VOLUME + (1.0 - blood_mix)*RAIN_VOLUME);
}



