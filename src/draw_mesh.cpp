// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/27/02

#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"
#include <GL/glext.h>


float const W_TEX_SCALE0     = 1.0;
float const WATER_WIND_EFF   = 0.0006;
float const VIEW_DIST0       = 4.0;
float const START_OFFSET0    = 0.05;
float const PRESP_ANGLE_ADJ  = 1.5;
float const VD_SCALE         = 1.0;
float const SURF_HEAL_RATE   = 0.005;
float const MAX_SURFD        = 20.0;
float const DLIGHT_SCALE     = 4.0;
float const WATER_COL_ATTEN  = 0.6;
unsigned const W_STEPS       = 16;
int   const DRAW_BORDER      = 3;

int   const SHOW_MESH_TIME   = 0;
int   const SHOW_NORMALS     = 0;
int   const DEBUG_COLLS      = 0;
int   const DISABLE_TEXTURES = 0;

int  const TILE_RADIUS       = 4; // WM0, in mesh sizes
int  const TILE_RADIUS_IT    = 5; // WM3, in mesh sizes
bool const DEBUG_TILES       = 0;


struct fp_ratio {
	float n, d;
	inline float get_val() const {return ((d == 0.0) ? 1.0 : n/d);}
};


// Global Variables
bool clear_landscape_vbo;
int island(0);
float lt_green_int(1.0), sm_green_int(1.0);
vector<fp_ratio> uw_mesh_lighting; // for water caustics

extern bool using_lightmap, has_dl_sources, combined_gu, has_snow, tiled_mesh_display;
extern unsigned num_jterms;
extern int draw_model, num_local_minima, world_mode, xoff, yoff, xoff2, yoff2, ocean_set, ground_effects_level;
extern int display_mode, frame_counter, resolution, verbose_mode, DISABLE_WATER, read_landscape, disable_inf_terrain;
extern float zmax, zmin, zmax_est, ztop, zbottom, light_factor, max_water_height, init_temperature;
extern float water_plane_z, temperature, fticks, mesh_scale, mesh_z_cutoff, TWO_XSS, TWO_YSS, XY_SCENE_SIZE;
extern point light_pos, litning_pos, sun_pos, moon_pos;
extern vector3d up_norm, wind;
extern float max_light[], h_dirt[];
extern texture textures[];


void draw_sides_and_bottom();



inline int test_draw_vertex(int i, int j, unsigned c) {

	if (c == 1) return 1;
	if (mesh_draw != NULL && (is_mesh_disabled(j, i) || is_mesh_disabled(j, i+1))) return 0;
	if (mesh_z_cutoff <= -FAR_CLIP) return 1;
	return (mesh_z_cutoff < max(mesh_height[i][j], mesh_height[i+1][j]));
}


float camera_min_dist_to_surface() { // min dist of four corners and center

	point pos;
	get_matrix_point(0, 0, pos);
	float dist(distance_to_camera(pos));
	get_matrix_point(MESH_X_SIZE-1, 0, pos);
	dist = min(dist, distance_to_camera(pos));
	get_matrix_point(0, MESH_Y_SIZE-1, pos);
	dist = min(dist, distance_to_camera(pos));
	get_matrix_point(MESH_X_SIZE-1, MESH_Y_SIZE-1, pos);
	dist = min(dist, distance_to_camera(pos));
	get_matrix_point(MESH_X_SIZE/2, MESH_Y_SIZE/2, pos);
	dist = min(dist, distance_to_camera(pos));
	return dist;
}


colorRGBA setup_mesh_lighting() {

	colorRGBA ambient_color(DEF_AMBIENT, DEF_AMBIENT, DEF_AMBIENT, 1.0);
	colorRGBA diffuse_color(DEF_DIFFUSE, DEF_DIFFUSE, DEF_DIFFUSE, 1.0);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &diffuse_color.red);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, &ambient_color.red);
	glEnable(GL_COLOR_MATERIAL);
	diffuse_color.do_glColor();
	set_fill_mode();
	enable_blend();
	return diffuse_color;
}


void setup_arrays(void *varr, void *narr, void *carr) {

	glDisable(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	if (!carr) glDisable(GL_COLOR_ARRAY);
	if (carr) glEnableClientState(GL_COLOR_ARRAY); else glDisableClientState(GL_COLOR_ARRAY);
	if (varr) glVertexPointer(3, GL_FLOAT, 0, varr);
	if (narr) glNormalPointer(   GL_FLOAT, 0, narr);
	if (carr) glColorPointer( 3, GL_FLOAT, 0, carr);
}


void run_post_mesh_draw() {
	
	glEnable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	disable_blend();
	glEnable(GL_TEXTURE_COORD_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glEnable(GL_COLOR_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
}


float integrate_water_dist(point const &targ_pos, point const &src_pos, float const water_z) {

	if (src_pos.z == targ_pos.z) return 0.0;
	float const t(min(1.0f, (water_z - targ_pos.z)/fabs(src_pos.z - targ_pos.z))); // min(1.0,...) for underwater case
	point p_int(targ_pos + (src_pos - targ_pos)*t);
	int const xp(get_xpos(targ_pos.x)), yp(get_ypos(targ_pos.y));
	if (!point_outside_mesh(xp, yp)) p_int.z = min(src_pos.z, water_matrix[yp][xp]); // account for ripples
	return p2p_dist(p_int, targ_pos)*mesh_scale;

}


inline void calc_norm(vector3d &norm, float *xv, float *yv, int xi, int yi, float z1, float z2, float z3, float z4) {

	float const dx(xv[xi-1] - xv[xi+1]), dy(yv[yi-1] - yv[yi+1]);
	norm.assign(dy*(z3 - z1), dx*(z4 - z2), dx*dy);
	//norm.normalize();
}


void water_color_atten_pt(float *c, int x, int y, point const &pos, point const &p1, point const &p2) {

	float const scale(WATER_COL_ATTEN*((wminside[y][x] == 2) ? 1.0 : 2.0)), wh(water_matrix[y][x]);
	float const dist(scale*(integrate_water_dist(pos, p1, wh) + integrate_water_dist(pos, p2, wh)));
	atten_by_water_depth(c, dist);
}


float get_cloud_shadow_atten(int x, int y) {

	point const pos(get_xval(x), get_yval(y), mesh_height[y][x]);
	float sval(0.0);
	
	// use the original sun/moon pos - it's wrong, but not too wrong and at least it's independent of the camera pos
	if (light_factor > 0.4) { // sun
		float const cloud_density(get_cloud_density(pos, (sun_pos - pos).get_norm()));
		sval += min(1.0,  5.0*(light_factor - 0.4))*(1.0 - CLIP_TO_01(1.7f*cloud_density));
	}
	if (light_factor < 0.6 && !combined_gu) { // moon
		float const cloud_density(get_cloud_density(pos, (moon_pos - pos).get_norm()));
		sval += min(1.0, -5.0*(light_factor - 0.6))*(1.0 - CLIP_TO_01(1.7f*cloud_density));
	}
	return sval;
};


inline float blend_light(bool has_sun, bool has_moon) {

	float const lfs(5.0*(light_factor - 0.4)); // 0: all moon, 1: all sun
	return lfs*has_sun + (1.0 - lfs)*has_moon;
}


class mesh_vertex_draw {

	static int fcount;
	float const healr;
	unsigned char **sml;
	vector<point>    varr;
	vector<vector3d> narr;
	vector<colorRGB> carr;

	struct norm_color_ix {
		vector3d n;
		colorRGB c;
		int ix;
		norm_color_ix() : ix(-1) {}
		norm_color_ix(vector3d const &n_, colorRGB const &c_, int ix_) : n(n_), c(c_), ix(ix_) {}
	};

	vector<norm_color_ix> last_rows;

	void update_color(int i, int j) {

		float color_scale(DEF_DIFFUSE);
		float &sd(surface_damage[i][j]);

		if (sd > 0.0) {
			sd = min(MAX_SURFD, max(0.0f, (sd - healr)));
			color_scale *= max(0.0f, (1.0f - sd));
		}
		carr[c].set_to_val(color_scale);

		if (DLIGHT_SCALE > 0.0 && (using_lightmap || has_dl_sources)) { // somewhat slow
			get_sd_light(j, i, get_zpos(varr[c].z), &varr[c][0], DLIGHT_SCALE, &carr[c].R, &surface_normals[i][j], NULL);
		}
	}

	void set_normal_array(int i, int j) {

		float light_scale(1.0);

		if (sml) { // sun or moon shadows
			light_scale = ((sml[i][j] & SHADOWED_ALL) ? 0.0 : 1.0);
		}
		else { // combined sun and moon shadows
			bool const no_sun ((shadow_mask[LIGHT_SUN ][i][j] & SHADOWED_ALL) != 0);
			bool const no_moon((shadow_mask[LIGHT_MOON][i][j] & SHADOWED_ALL) != 0);
			light_scale = blend_light(!no_sun, !no_moon);
		}

		// water light attenuation: total distance from sun/moon, reflected off bottom, to viewer
		if (!DISABLE_WATER && varr[c][2] < max_water_height && varr[c][2] < water_matrix[i][j]) {
			point const pos(get_xval(j), get_yval(i), mesh_height[i][j]);
			water_color_atten(((float *)&(carr[c].R)), j, i, pos);

			if (wminside[i][j] == 1) { // too slow?
				colorRGBA wc(WHITE);
				select_liquid_color(wc, j, i);
				UNROLL_3X(carr[c][i_] *= wc[i_];)
			}
			
			// water caustics: slow and low resolution, but conceptually interesting
			if (light_scale > 0.0 && !uw_mesh_lighting.empty()) {
				float const val(uw_mesh_lighting[i*MESH_X_SIZE + j].get_val());
				//light_scale *= val*val; // square to enhance the caustics effect
				light_scale = pow(val, 8);
			}
		}
		if (ground_effects_level >= 2 && !has_snow && light_scale > 0.0 && !(display_mode & 0x10)) {
			light_scale *= get_cloud_shadow_atten(j, i);
		}
		narr[c] = ((light_scale == 0.0) ? vector3d(0.0, 0.0, 0.0) : vertex_normals[i][j]*light_scale);
	}

public:
	unsigned c;

	mesh_vertex_draw() : healr(fticks*SURF_HEAL_RATE), sml(NULL),
		varr(2*(MAX_XY_SIZE+1)), narr(varr.size()), carr(varr.size()), c(0)
	{
		assert(shadow_mask != NULL);
		
		if (light_factor >= 0.6) { // sun shadows
			sml = shadow_mask[LIGHT_SUN];
		}
		else if (light_factor <= 0.4) { // moon shadows
			sml = shadow_mask[LIGHT_MOON];
		}
		last_rows.resize(MESH_X_SIZE+1);
		setup_arrays(&varr.front(), &narr.front(), &carr.front());
	}

	bool draw_mesh_vertex_pair(int i, int j, float x, float y) {

		if (!test_draw_vertex(i, j, c)) return 0;
		
		for (unsigned p = 0; p < 2; ++p, ++c) {
			int const iinc(min((MESH_Y_SIZE-1), int(i+p)));
			assert(c < varr.size());
			varr[c].assign(x, (y + p*DY_VAL), mesh_height[iinc][j]);
			assert(unsigned(j) < last_rows.size());
		
			if (last_rows[j].ix == iinc) { // gets here nearly half the time
				narr[c] = last_rows[j].n;
				carr[c] = last_rows[j].c;
			}
			else {
				update_color(iinc, j);
				set_normal_array(iinc, j);
				last_rows[j] = norm_color_ix(narr[c], carr[c], iinc);
			}
		}
		return 1;
	}
};

int mesh_vertex_draw::fcount(0);


void gen_uw_lighting() {

	uw_mesh_lighting.resize(MESH_X_SIZE*MESH_Y_SIZE);
	point const lpos(get_light_pos());
	float const ssize(X_SCENE_SIZE + Y_SCENE_SIZE + Z_SCENE_SIZE);
	vector<point> rows[2]; // {last, current} y rows
	for (unsigned i = 0; i < 2; ++i) rows[i].resize(MESH_X_SIZE, all_zeros);
	float const dxy_val_inv[2] = {DX_VAL_INV, DY_VAL_INV};

	for (vector<fp_ratio>::iterator i = uw_mesh_lighting.begin(); i != uw_mesh_lighting.end(); ++i) {
		i->n = i->d = 0.0; // initialize
	}
	for (int y = 0; y < MESH_Y_SIZE; ++y) {
		for (int x = 0; x < MESH_X_SIZE; ++x) {
			if (!mesh_is_underwater(x, y)) continue;
			point const p1(get_xval(x), get_yval(y), water_matrix[y][x]); // point on water surface
			vector3d const dir(p1 - lpos);
			vector3d v_refract(dir);
			calc_refraction_angle(dir, v_refract, wat_vert_normals[y][x], 1.0, WATER_INDEX_REFRACT);
			point const p2(p1 + v_refract.get_norm()*ssize); // distant point along refraction vector
			int xpos(0), ypos(0);
			float zval;

			if (p1.z == p2.z || !line_intersect_mesh(p1, p2, xpos, ypos, zval, 1, 1)) continue; // no intersection
			assert(!point_outside_mesh(xpos, ypos));
			float const t((zval - p1.z)/(p2.z - p1.z));
			point const cpos(p1 + (p2 - p1)*t); // collision point with underwater mesh
			rows[1][x] = cpos;
			if (x == 0 || y == 0) continue; // not an interior point
			if (rows[0][x].z == 0.0 || rows[0][x-1].z == 0.0 || rows[1][x-1].z == 0.0) continue; // incomplete block
			float rng[2][2] = {{cpos.x, cpos.y}, {cpos.x, cpos.y}}; // {min,max} x {x,y} - bounds of mesh surface light patch through this water patch
			int bnds[2][2]; // {min,max} x {x,y} - integer mesh index bounds

			for (unsigned d = 0; d < 2; ++d) { // x,y
				for (unsigned e = 0; e < 2; ++e) { // last,cur y (row)
					for (unsigned f = 0; f < 2; ++f) { // last,cur x
						rng[0][d] = min(rng[0][d], rows[e][x-f][d]);
						rng[1][d] = max(rng[1][d], rows[e][x-f][d]);
					}
				}
				assert(rng[0][d] < rng[1][d]);
				bnds[0][d] = max(0,              int(floor((rng[0][d] + SCENE_SIZE[d])*dxy_val_inv[d])));
				bnds[1][d] = min(MESH_SIZE[d]-1, int(ceil ((rng[1][d] + SCENE_SIZE[d])*dxy_val_inv[d])));
				assert(bnds[0][d] <= bnds[1][d]); // can this fail?
			}
			float const weight_n(1.0/((rng[1][0] - rng[0][0])*(rng[1][1] - rng[0][1]))); // weight of this patch of light
			float const weight_d(1.0/(DX_VAL*DY_VAL));
			float const init_cr[2] = {(-X_SCENE_SIZE + DX_VAL*bnds[0][0]), (-Y_SCENE_SIZE + DY_VAL*bnds[0][1])};
			float crng[2][2]; // {min,max} x {x,y} - range of this mesh quad
			crng[0][1] = init_cr[1];
			crng[1][1] = init_cr[1] + DY_VAL;

			for (int yy = bnds[0][1]; yy < bnds[1][1]; ++yy) {
				assert(yy >= 0 && yy < MESH_Y_SIZE);
				float const ysz(min(crng[1][1], rng[1][1]) - max(crng[0][1], rng[0][1])); // intersection: min(UB) - max(LB)
				assert(ysz > 0.0);
				crng[0][0] = init_cr[0];
				crng[1][0] = init_cr[0] + DX_VAL;

				for (int xx = bnds[0][0]; xx < bnds[1][0]; ++xx) {
					assert(xx >= 0 && xx < MESH_X_SIZE);
					float const xsz(min(crng[1][0], rng[1][0]) - max(crng[0][0], rng[0][0])); // intersection: min(UB) - max(LB)
					assert(xsz > 0.0);
					unsigned const ix(yy*MESH_X_SIZE + xx);
					uw_mesh_lighting[ix].n += weight_n*xsz*ysz; // amount of light through patch of water hitting this mesh quad
					uw_mesh_lighting[ix].d += weight_d*xsz*ysz;
					crng[0][0] += DX_VAL;
					crng[1][0] += DX_VAL;
				}
				crng[0][1] += DY_VAL;
				crng[1][1] += DY_VAL;
			}
		} // for x
		rows[0].swap(rows[1]);
		for (int x = 0; x < MESH_X_SIZE; ++x) rows[1][x] = all_zeros; // reset last row
	} // for y
}


void set_landscape_texgen(float tex_scale, int xoff, int yoff, int xsize, int ysize) {

	float const tx(tex_scale*(((float)xoff)/((float)xsize) + 0.5));
	float const ty(tex_scale*(((float)yoff)/((float)ysize) + 0.5));
	setup_texgen(tex_scale/TWO_XSS, tex_scale/TWO_YSS, tx, ty);
}


void draw_coll_vert(int i, int j) {

	(v_collision_matrix[i][j].cvz.empty() ? BLUE : RED).do_glColor();
	glVertex3f(get_xval(j), get_yval(i), max(czmin, v_collision_matrix[i][j].zmax));
}


void display_mesh() { // fast array version

	if (mesh_height == NULL) return; // no mesh to display
	// can't put the hole in the right place, so only draw the tiled terrain
	if (tiled_mesh_display && !island && (display_mode & 0x10) && (xoff2 != 0 || yoff2 != 0)) return;
	RESET_TIME;

	if ((display_mode & 0x80) && !DISABLE_WATER && !ocean_set && zmin < max_water_height) {
		gen_uw_lighting();
		if (SHOW_MESH_TIME) PRINT_TIME("Underwater Lighting");
	}
	else {
		uw_mesh_lighting.clear();
	}
	if (DEBUG_COLLS) {
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
			glBegin(GL_TRIANGLE_STRIP);
			
			for (int j = 0; j < MESH_X_SIZE; ++j) {
				draw_coll_vert(i+0, j);
				draw_coll_vert(i+1, j);
			}
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}
	setup_mesh_lighting();
	update_landscape_texture();
	if (SHOW_MESH_TIME) PRINT_TIME("Landscape Texture");
	set_landscape_texgen(1.0, xoff, yoff, MESH_X_SIZE, MESH_Y_SIZE);
	glDisable(GL_NORMALIZE);
	if (!DISABLE_TEXTURES) select_texture(LANDSCAPE_TEX);
	if (SHOW_MESH_TIME) PRINT_TIME("Preprocess");

	if (ground_effects_level == 0 && setup_gen_buffers_arb()) {
		static unsigned mesh_vbo(0);
		
		if (clear_landscape_vbo) {
			delete_vbo(mesh_vbo);
			mesh_vbo = 0;
			clear_landscape_vbo = 0;
		}
		if (mesh_vbo == 0) {
			vector<point> data; // vertex and normals

			for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
				for (int j = 0; j < MESH_X_SIZE; ++j) {
					for (unsigned k = 0; k < 2; ++k) {
						data.push_back(point(get_xval(j), get_yval(i+k), mesh_height[i+k][j]));
						data.push_back(vertex_normals[i+k][j]);
					}
				}
			}
			mesh_vbo = create_vbo();
			assert(mesh_vbo > 0);
			bind_vbo(mesh_vbo);
			upload_vbo_data(&data.front(), data.size()*sizeof(point));
		}
		else {
			bind_vbo(mesh_vbo);
		}
		glDisable(GL_TEXTURE_COORD_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisable(GL_COLOR_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glVertexPointer(3, GL_FLOAT, 2*sizeof(point), 0);
		glNormalPointer(   GL_FLOAT, 2*sizeof(point), (void *)sizeof(point));

		for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
			glDrawArrays(GL_QUAD_STRIP, 2*i*MESH_X_SIZE, 2*MESH_X_SIZE);
		}
		bind_vbo(0);
	}
	else {
		float y(-Y_SCENE_SIZE);
		mesh_vertex_draw mvd;

		for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
			float x(-X_SCENE_SIZE);
			mvd.c = 0;

			for (int j = 0; j < MESH_X_SIZE-1; ++j) {
				if (!mvd.draw_mesh_vertex_pair(i, j, x, y) && mvd.c > 0) {
					glDrawArrays(GL_TRIANGLE_STRIP, 0, mvd.c);
					mvd.c = 0;
				}
				x += DX_VAL;
			} // for j
			mvd.draw_mesh_vertex_pair(i, (MESH_X_SIZE - 1), x, y);
			if (mvd.c > 1) glDrawArrays(GL_TRIANGLE_STRIP, 0, mvd.c);
			y += DY_VAL;
		} // for i
	}
	if (SHOW_MESH_TIME) PRINT_TIME("Draw");
	disable_textures_texgen();
	glDisable(GL_COLOR_MATERIAL);
	if (!island && world_mode != WMODE_INF_TERRAIN) draw_sides_and_bottom();
	run_post_mesh_draw();

	if (SHOW_NORMALS) {
		set_color(RED);
		glNormal3f(0.0, 0.0, 1.0);
		glBegin(GL_LINES);

		for (int i = 1; i < MESH_Y_SIZE-2; ++i) {
			for (int j = 1; j < MESH_X_SIZE-1; ++j) {
				point const pos(get_xval(j), get_yval(i), mesh_height[i][j]);
				vector3d const norm(vertex_normals[i][j]*0.1);
				//vector3d const norm(get_local_wind(pos + vector3d(0.0, 0.0, DZ_VAL))*0.1);
				glVertex3f(pos.x, pos.y, pos.z);
				glVertex3f(pos.x+norm.x, pos.y+norm.y, pos.z+norm.z);
			}
		}
		glEnd();
	}
	if (SHOW_MESH_TIME) PRINT_TIME("Final");
}


int set_texture(float zval, int &tex_id) {

	int id;
	for (id = 0; id < NTEX_DIRT-1 && zval >= (h_dirt[id]*2.0*zmax_est - zmax_est); ++id) {}
	update_lttex_ix(id);
	bool const changed(tex_id != lttex_dirt[id].id);
	tex_id = lttex_dirt[id].id;
	return changed;
}


unsigned get_norm_texels() {
	return (world_mode == WMODE_INF_TERRAIN) ? 512 : get_texture_size(LANDSCAPE_TEX, 0);
}


float display_mesh3(int const *const hole_bounds) { // WM3 - infinite terrain

	if (tiled_mesh_display) {
		bool const add_hole((hole_bounds != NULL) && xoff2 == 0 && yoff2 == 0);
		return draw_tiled_terrain(add_hole);
	}

	//RESET_TIME;
	float const view_dist(Z_SCENE_SIZE*VIEW_DIST0);
	static float xv[DYNAMIC_MESH_SZ], yv[DYNAMIC_MESH_SZ], xv2[DYNAMIC_MESH_SZ], yv2[DYNAMIC_MESH_SZ], last_h[2][DYNAMIC_MESH_SZ];
	static vector3d last_n[DYNAMIC_MESH_SZ];
	int const ssize((int)pow(RES_STEP, resolution-1));
	float const step_size(STEP_SIZE*ssize);
	float const x0(-STEP_SIZE*X_SCENE_SIZE), y0(-STEP_SIZE*Y_SCENE_SIZE);
	float const xstep(step_size*DX_VAL), ystep(step_size*DY_VAL), ocxl(-cview_dir.x), ocyl(-cview_dir.y);
	float const sz_off(((float)ssize)*Z_SCENE_SIZE*START_OFFSET0), xc(sz_off*ocxl), yc(sz_off*ocyl);
	float const xo(xc - VD_SCALE*view_dist*ocxl), yo(yc - VD_SCALE*view_dist*ocyl);
	float const d_far(VD_SCALE*view_dist*tan(0.5*PRESP_ANGLE_ADJ*PERSP_ANGLE*TO_RADIANS));
	float const xfar1(xo - d_far*ocyl), xfar2(xo + d_far*ocyl), yfar1(yo + d_far*ocxl), yfar2(yo - d_far*ocxl);
	float const bbx1(min(xc, min(xfar1, xfar2))), bbx2(max(xc, max(xfar1, xfar2)));
	float const bby1(min(yc, min(yfar1, yfar2))), bby2(max(yc, max(yfar1, yfar2)));
	int const x1(int((bbx1 - x0)/xstep) - DRAW_BORDER), x2(int((bbx2 - x0)/xstep) + DRAW_BORDER);
	int const y1(int((bby1 - y0)/ystep) - DRAW_BORDER), y2(int((bby2 - y0)/ystep) + DRAW_BORDER);
	int const nx(x2 - x1), ny(y2 - y1), xoff3(ssize*(xoff2/ssize)), yoff3(ssize*(yoff2/ssize));
	float zmin2(FAR_CLIP), zmax2(-FAR_CLIP);
	if (nx == 0 || ny == 0) return zmin2;
	assert(nx > 2);
	if (xoff2 != xoff3 || yoff2 != yoff3) glTranslatef(DX_VAL*(xoff3 - xoff2), DY_VAL*(yoff3 - yoff2), 0.0);
	float x(x0 + (x1 + 1)*xstep), y(y0 + (y1 + 1)*ystep);

	for (int i = 0; i < ny; ++i) {
		yv[i]  = y;
		yv2[i] = y + (yoff3 + 0.5)*DY_VAL; // FIXME: why the +0.5, and why is it -0.5 in x?
		y     += ystep;
	}
	for (int j = 0; j < nx; ++j) {
		xv[j]  = x;
		xv2[j] = x + (xoff3 - 0.5)*DX_VAL;
		x     += xstep;
	}
	//PRINT_TIME("Init");
	build_xy_mesh_arrays(xv2, yv2, nx, ny);
	//PRINT_TIME("Array Build");

	for (int j = 1; j < nx-1; ++j) {
		last_h[0][j] = fast_eval_from_index(j, 1, 1);
		last_h[1][j] = fast_eval_from_index(j, 0, 1);
		calc_norm(last_n[j], xv, yv, j, 0, last_h[0][j-1], last_h[1][j], last_h[0][j+1], fast_eval_from_index(j, 2, 1));
	}
	vector3d norm;
	unsigned const nverts(2*(nx-2));
	vector<point>    varr(nverts);
	vector<vector3d> narr(nverts);

	int const tex_xoff(xoff + xoff3 - xoff2), tex_yoff(yoff + yoff3 - yoff2);
	unsigned const norm_texels(get_norm_texels());
	unsigned last_tsize(0);
	int tex_id(-1);
	set_texture(last_h[0][1], tex_id);
	last_tsize = get_texture_size(tex_id, 0);
	set_landscape_texgen(((float)norm_texels)/last_tsize, tex_xoff, tex_yoff, MESH_X_SIZE, MESH_Y_SIZE);
	select_texture(tex_id);
	setup_mesh_lighting();
	setup_arrays(&varr.front(), &narr.front(), NULL);
	//PRINT_TIME("Htable");
	glPushMatrix();
	glTranslatef(xoff*DX_VAL, yoff*DY_VAL, 0.0);
	
	for (int i = 1; i < ny-2; ++i) {
		unsigned c(0);

		for (int j = 1; j < nx-1; ++j) {
			float const zval(fast_eval_from_index(j, i+1, 1)); // FIXME: i+1?
			calc_norm(norm, xv, yv, j, i, last_h[0][j-1], last_h[1][j], last_h[0][j+1], zval);

			for (unsigned p = 0; p < 2; ++p, ++c) {
				varr[c].assign(xv[j], yv[i+p-1], last_h[!p][j]);
				narr[c] = (p ? norm : last_n[j]);
			}
			int const xpos(j + x1), ypos(i + y1);
			
			if (hole_bounds != NULL && ypos >= hole_bounds[2] && ypos < hole_bounds[3] && xpos+1 >= hole_bounds[0] && xpos+2 <= hole_bounds[1]) {
				if (c > 2) glDrawArrays(GL_TRIANGLE_STRIP, 0, c);
				c = 0;
			}
			else if (set_texture(last_h[0][j+1], tex_id)) {
				if (c > 2) glDrawArrays(GL_TRIANGLE_STRIP, 0, c);
				unsigned const new_tsize(get_texture_size(tex_id, 0));

				if (new_tsize != last_tsize) {
					set_landscape_texgen(((float)norm_texels)/new_tsize, tex_xoff, tex_yoff, MESH_X_SIZE, MESH_Y_SIZE);
					last_tsize = new_tsize;
				}
				select_texture(tex_id);
				
				for (unsigned q = 0; q < 2; ++q) {
					varr[q] = varr[c+q-2];
					narr[q] = narr[c+q-2];
				}
				c = 2;
			}
			last_h[1][j] = last_h[0][j];
			last_h[0][j] = zval;
			last_n[j]    = norm;
			zmin2        = min(zmin2, zval);
			zmax2        = max(zmax2, zval);
		}
		assert(c <= nverts);
		if (c > 2) glDrawArrays(GL_TRIANGLE_STRIP, 0, c);
	}
	//PRINT_TIME("Render");
	disable_textures_texgen();
	run_post_mesh_draw();

	if (SHOW_NORMALS) {
		set_color(RED);
		glNormal3f(0.0, 0.0, 1.0);
		glBegin(GL_LINES);

		for (int i = 1; i < ny-2; ++i) {
			for (int j = 1; j < nx-1; ++j) {
				float const zval(fast_eval_from_index(j, i, 1));
				vector3d norm;
				calc_norm(norm, xv, yv, j, i, fast_eval_from_index(j-1, i, 1), fast_eval_from_index(j, i-1, 1),
					fast_eval_from_index(j+1, i, 1), fast_eval_from_index(j, i+1, 1));
				glVertex3f(xv[j], yv[i], zval);
				glVertex3f(xv[j]+norm.x, yv[i]+norm.y, zval+norm.z);
			}
		}
		glEnd();
	}
	glPopMatrix();
	if (verbose_mode && frame_counter%100 == 0) cout << 2*(nx-1)*(ny-2) << " triangles." << endl;
	//PRINT_TIME("Done");
	return zmin2;
}


inline void draw_vertex(float x, float y, float z, bool in_y, float tscale=1.0) { // xz or zy

	glTexCoord2f(tscale*(in_y ? z : x), tscale*(in_y ? y : z));
	glVertex3f(x, y, z);
}


// NOTE: There is a buffer of one unit around the drawn area
void draw_sides_and_bottom() {

	if (!disable_inf_terrain && (display_mode & 0x10)) return;
	int const lx(MESH_X_SIZE-1), ly(MESH_Y_SIZE-1);
	float const botz(zbottom - 0.05), z_avg(0.5*(zbottom + ztop)), ts(4.0/(X_SCENE_SIZE + Y_SCENE_SIZE));
	float const x1(-X_SCENE_SIZE), y1(-Y_SCENE_SIZE), x2(X_SCENE_SIZE-DX_VAL), y2(Y_SCENE_SIZE-DY_VAL);
	int const texture((!read_landscape && get_rel_height(z_avg, zmin, zmax) > lttex_dirt[2].zval) ? ROCK_TEX : DIRT_TEX);
	set_color(WHITE);
	set_lighted_sides(2);
	set_fill_mode();
	
	if (!DISABLE_TEXTURES) select_texture(texture);
	glBegin(GL_QUADS);
	glNormal3f(0.0, 0.0, -1.0); // bottom surface
	draw_one_tquad(x1, y1, x2, y2, botz, 1, ts*x1, ts*y1, ts*x2, ts*y2);
	glEnd();
	{
		float xv(x1);
		glNormal3f(0.0, -1.0, 0.0);
		glBegin(GL_QUADS);
		unsigned steps(0);

		for (int i = 1; i < MESH_X_SIZE; ++i) { // y sides
			for (unsigned d = 0; d < 2; ++d) {
				int const xy_ix(d ? ly : 0);
				float const limit(d ? y2 : y1);
				draw_vertex(xv,        limit, botz, 0, ts);
				draw_vertex(xv+DX_VAL, limit, botz, 0, ts);
				draw_vertex(xv+DX_VAL, limit, mesh_height[xy_ix][i], 0, ts);
				draw_vertex(xv,        limit, mesh_height[xy_ix][i-1], 0, ts);
			}
			xv += DX_VAL;
		}
		glEnd();
	}
	{
		float yv(y1);
		glNormal3f(1.0, 0.0, 0.0);
		glBegin(GL_QUADS);

		for (int i = 1; i < MESH_Y_SIZE; ++i) { // x sides
			for (unsigned d = 0; d < 2; ++d) {
				int const xy_ix(d ? lx : 0);
				float const limit(d ? x2 : x1);
				draw_vertex(limit, yv,        botz, 1, ts);
				draw_vertex(limit, yv+DY_VAL, botz, 1, ts);
				draw_vertex(limit, yv+DY_VAL, mesh_height[i][xy_ix], 1, ts);
				draw_vertex(limit, yv,        mesh_height[i-1][xy_ix], 1, ts);
			}
			yv += DY_VAL;
		}
		glEnd();
	}
	glNormal3f(0.0, 0.0, 1.0);
	set_lighted_sides(1);
	glDisable(GL_TEXTURE_2D);
}


class water_renderer {

	int check_zvals;
	float tex_scale;

	void draw_x_sides(bool neg_edge) const;
	void draw_y_sides(bool neg_edge) const;
	void draw_sides(unsigned ix) const;

public:
	water_renderer(int ix, int iy, int cz) : check_zvals(cz), tex_scale(W_TEX_SCALE0/Z_SCENE_SIZE) {}
	void draw() const;
};


void water_renderer::draw_x_sides(bool neg_edge) const {

	int const end_val(neg_edge ? 0 : MESH_X_SIZE-1);
	float const limit(neg_edge ? -X_SCENE_SIZE : X_SCENE_SIZE-DX_VAL);
	float yv(-Y_SCENE_SIZE);
	glNormal3f(-1.0, 0.0, 0.0);
	glBegin(GL_QUADS);

	for (int i = 1; i < MESH_Y_SIZE; ++i) { // x sides
		float const mh1(mesh_height[i][end_val]), mh2(mesh_height[i-1][end_val]);
		float const wm1(water_matrix[i][end_val] - SMALL_NUMBER), wm2(water_matrix[i-1][end_val] - SMALL_NUMBER);

		if (!check_zvals || mh1 < wm1 || mh2 < wm2) {
			draw_vertex(limit, yv,        wm2,           1, tex_scale);
			draw_vertex(limit, yv+DY_VAL, wm1,           1, tex_scale);
			draw_vertex(limit, yv+DY_VAL, min(wm1, mh1), 1, tex_scale);
			draw_vertex(limit, yv,        min(wm2, mh2), 1, tex_scale);
		}
		yv += DY_VAL;
	}
	glEnd();
}

void water_renderer::draw_y_sides(bool neg_edge) const {

	int const end_val(neg_edge ? 0 : MESH_Y_SIZE-1);
	float const limit(neg_edge ? -Y_SCENE_SIZE : Y_SCENE_SIZE-DY_VAL);
	float xv(-X_SCENE_SIZE);
	glNormal3f(0.0, 1.0, 0.0);
	glBegin(GL_QUADS);
	
	for (int i = 1; i < MESH_X_SIZE; ++i) { // y sides
		float const mh1(mesh_height[end_val][i]), mh2(mesh_height[end_val][i-1]);
		float const wm1(water_matrix[end_val][i] - SMALL_NUMBER), wm2(water_matrix[end_val][i-1] - SMALL_NUMBER);

		if (!check_zvals || mh1 < wm1 || mh2 < wm2) {
			draw_vertex(xv,        limit, wm2,           0, tex_scale);
			draw_vertex(xv+DX_VAL, limit, wm1,           0, tex_scale);
			draw_vertex(xv+DX_VAL, limit, min(wm1, mh1), 0, tex_scale);
			draw_vertex(xv,        limit, min(wm2, mh2), 0, tex_scale);
		}
		xv += DX_VAL;
	}
	glEnd();
}


void water_renderer::draw_sides(unsigned ix) const {

	switch (ix) { // xn xp yn yp
		case 0: draw_x_sides(1); break;
		case 1: draw_x_sides(0); break;
		case 2: draw_y_sides(1); break;
		case 3: draw_y_sides(0); break;
		default: assert(0);
	}
}


void water_renderer::draw() const {

	colorRGBA color;
	select_water_ice_texture(color);
	set_color(color);
	set_fill_mode();
	set_lighted_sides(2);
	enable_blend();
	point const camera(get_camera_pos());
	float const pts[4][2] = {{-X_SCENE_SIZE, 0.0}, {X_SCENE_SIZE, 0.0}, {0.0, -Y_SCENE_SIZE}, {0.0, Y_SCENE_SIZE}};
	vector<pair<float, unsigned> > sides(4);

	for (unsigned i = 0; i < 4; ++i) {
		sides[i] = make_pair(-distance_to_camera_sq(point(pts[i][0], pts[i][1], water_plane_z)), i);
	}
	sort(sides.begin(), sides.end()); // largest to smallest distance

	for (unsigned i = 0; i < 4; ++i) {
		draw_sides(sides[i].second); // draw back to front
	}
	glNormal3f(0.0, 0.0, 1.0);
	disable_blend();
	set_specular(0.0, 1.0);
	set_lighted_sides(1);
	glDisable(GL_TEXTURE_2D);
}


void draw_water_sides(int check_zvals) {

	if (!disable_inf_terrain && display_mode & 0x10) return;
	water_renderer wr(resolution, resolution, check_zvals);
	wr.draw();
}


int get_tile_radius() {

	return ((world_mode == WMODE_INF_TERRAIN) ? TILE_RADIUS_IT : TILE_RADIUS);
}


void draw_water_plane(float zval, int const *const hole_bounds, bool disable_lighting, bool large_size) {

	if (DISABLE_WATER) return;
	float const tscale(W_TEX_SCALE0/Z_SCENE_SIZE);
	float const vd_scale(large_size ? (tiled_mesh_display ? 2.0*get_tile_radius() : VIEW_DIST0)*SQRT2 : X_SCENE_SIZE/(X_SCENE_SIZE + DX_VAL));
	float const dx(xoff*DX_VAL), dy(yoff*DY_VAL);
	float const vdx(vd_scale*X_SCENE_SIZE), vdy(vd_scale*Y_SCENE_SIZE);
	static float wxoff(0.0), wyoff(0.0);
	colorRGBA color;
	select_water_ice_texture(color, ((world_mode == WMODE_INF_TERRAIN) ? &init_temperature : &temperature));
	color.alpha *= ((display_mode & 0x20) ? 0.5 : 0.8);

	if (temperature > W_FREEZE_POINT) {
		wxoff -= WATER_WIND_EFF*wind.x*fticks;
		wyoff -= WATER_WIND_EFF*wind.y*fticks;
	}
	if (disable_lighting) {
		glDisable(GL_LIGHTING);
		color.do_glColor();
	}
	else {
		plus_z.do_glNormal();
		set_color(color);
	}
	set_fill_mode();
	set_lighted_sides(2);
	enable_blend();
	setup_texgen(tscale, tscale, (tscale*(xoff2 - xoff)*DX_VAL + wxoff), (tscale*(yoff2 - yoff)*DY_VAL + wyoff));
	glPushMatrix();
	glTranslatef(0.0, 0.0, zval);
	glBegin(GL_QUADS);

	if (hole_bounds) { // x1 x2 y1 y2
		assert(large_size);
		float const obnd[2][2] = {{dx-vdx, dx+vdx}, {dy-vdy, dy+vdy}}; // {x,y} x {1,2}
		float ibnd[2][2]; // {x,y} x {1,2}

		for (unsigned i = 0; i < 2; ++i) {
			ibnd[0][i] = get_xval(hole_bounds[i+0]);
			ibnd[1][i] = get_yval(hole_bounds[i+2]);
		}
		for (unsigned i = 0; i < 2; ++i) {
			glVertex2f(obnd[0][0], obnd[1][i]);
			glVertex2f(obnd[0][1], obnd[1][i]);
			glVertex2f(obnd[0][1], ibnd[1][i]);
			glVertex2f(obnd[0][0], ibnd[1][i]);

			glVertex2f(obnd[0][i], ibnd[1][0]);
			glVertex2f(ibnd[0][i], ibnd[1][0]);
			glVertex2f(ibnd[0][i], ibnd[1][1]);
			glVertex2f(obnd[0][i], ibnd[1][1]);
			if (!disable_lighting) glNormal3f(0.0, 0.0, -1.0); // invert the normal
		}
	}
	else {
		float const xinc(2.0*vdx/(float)W_STEPS); // about 9
		float const yinc(2.0*vdy/(float)W_STEPS);
		float yval(dy - vdy);

		for (unsigned i = 0; i < W_STEPS; ++i) {
			float xval(dx - vdx);

			for (unsigned j = 0; j < W_STEPS; ++j) {
				//if (hole_bounds) {}
				glVertex2f( xval,        yval);
				glVertex2f((xval+xinc),  yval);
				glVertex2f((xval+xinc), (yval+yinc));
				glVertex2f( xval,       (yval+yinc));
				xval += xinc;
			}
			yval += yinc;
		}
	}
	glEnd();
	glPopMatrix();
	disable_blend();
	set_specular(0.0, 1.0);
	set_lighted_sides(1);
	disable_textures_texgen();
	glEnable(GL_LIGHTING);
}


// ************* Tile Drawing Code *************


inline float get_water_atten_factor(float mh) {
	return 2.5*WATER_COL_ATTEN*(water_plane_z - mh)*mesh_scale;
}


struct tile_xy_pair {

	int x, y;
	tile_xy_pair(int x_=0, int y_=0) : x(x_), y(y_) {}
	bool operator<(tile_xy_pair const &t) const {return ((y == t.y) ? (x < t.x) : (y < t.y));}
};


class tile_t {

	int x1, y1, x2, y2, init_dxoff, init_dyoff;
	unsigned tid, vbo, ivbo, size, stride, zvsize, base_tsize, gen_tsize;
	float radius, mzmin, mzmax, xstart, ystart, xstep, ystep;
	vector<float> zvals;

public:
	tile_t() : tid(0), vbo(0), ivbo(0), size(0), stride(0), zvsize(0), gen_tsize(0) {}
	~tile_t() {clear_vbo_tid(1,1);}
	float get_zmin() const {return mzmin;}
	float get_zmax() const {return mzmax;}
	
	point get_center() const {
		return point(get_xval(((x1+x2)>>1) + (xoff - xoff2)), get_yval(((y1+y2)>>1) + (yoff - yoff2)), 0.5*(mzmin + mzmax));
	}

	void calc_start_step(int dx, int dy) {
		xstart = get_xval(x1 + dx);
		ystart = get_yval(y1 + dy);
		xstep  = (get_xval(x2 + dx) - xstart)/size;
		ystep  = (get_yval(y2 + dy) - ystart)/size;
	}
	
	tile_t(unsigned size_, int x, int y) : init_dxoff(xoff - xoff2), init_dyoff(yoff - yoff2),
		tid(0), vbo(0), ivbo(0), size(size_), stride(size+1), zvsize(stride+1), gen_tsize(0)
	{
		assert(size > 0);
		x1 = x*size;
		y1 = y*size;
		x2 = x1 + size;
		y2 = y1 + size;
		calc_start_step(0, 0);
		radius = 0.5*sqrt(xstep*xstep + ystep*ystep)*size; // approximate (lower bound)
		mzmin  = mzmax = get_camera_pos().z;
		base_tsize = get_norm_texels();
		
		if (DEBUG_TILES) {
			cout << "create " << size << ": " << x << "," << y << ", coords: " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
		}
	}

	unsigned get_gpu_memory() const {
		unsigned mem(0);
		if (vbo  > 0) mem += 2*stride*size*sizeof(vert_norm);
		if (ivbo > 0) mem += size*size*sizeof(unsigned short);
		if (tid  > 0) mem += 3*gen_tsize*gen_tsize;
		return mem;
	}

	void clear_tid() {
		if (tid > 0) glDeleteTextures(1, &tid);
		tid       = 0;
		gen_tsize = 0;
	}

	void clear_vbo_tid(bool vclear, bool tclear) {
		if (vclear) {
			delete_vbo(vbo);
			delete_vbo(ivbo);
			vbo = ivbo = 0;
		}
		if (tclear) {clear_tid();}
	}

	void create_zvals() {
		RESET_TIME;
		zvals.resize(zvsize*zvsize);
		static vector<float> xv, yv; // move somewhere else?
		xv.resize(zvsize);
		yv.resize(zvsize);
		calc_start_step(0, 0);
		mzmin =  FAR_CLIP;
		mzmax = -FAR_CLIP;

		for (unsigned i = 0; i < zvsize; ++i) { // not sure about this -x, +y thing
			xv[i] = (xstart + (i - 0.5)*xstep);
			yv[i] = (ystart + (i + 0.5)*ystep);
		}
		build_xy_mesh_arrays(&xv.front(), &yv.front(), zvsize, zvsize);

		for (unsigned y = 0; y < zvsize; ++y) {
			for (unsigned x = 0; x < zvsize; ++x) {
				float const zval(fast_eval_from_index(x, y, 0, 1));
				zvals[y*zvsize + x] = zval;
				mzmin = min(mzmin, zval);
				mzmax = max(mzmax, zval);
			}
		}
		assert(mzmin <= mzmax);
		radius = 0.5*sqrt((xstep*xstep + ystep*ystep)*size*size + (mzmax - mzmin)*(mzmax - mzmin));
		//PRINT_TIME("Create Zvals");
	}

	inline vector3d get_norm(unsigned ix) const {
		return vector3d(DY_VAL*(zvals[ix] - zvals[ix + 1]), DX_VAL*(zvals[ix] - zvals[ix + zvsize]), dxdy).get_norm();
	}

	void create_data(vector<vert_norm> &data, vector<unsigned short> &indices, vector<unsigned char> smask[]) {
		RESET_TIME;
		assert(zvals.size() == zvsize*zvsize);
		data.resize(stride*stride);
		indices.resize(4*size*size);
		calc_start_step(init_dxoff, init_dyoff);

		// FIXME: Need to shadow across adjacent tiles and update when new tiles are created
		bool const has_sun(light_factor >= 0.4), has_moon(light_factor <= 0.6);
		assert(has_sun || has_moon);
		
		if (has_sun)  {
			smask[LIGHT_SUN ].resize(zvals.size());
			calc_mesh_shadows(LIGHT_SUN,  sun_pos,  &zvals.front(), &smask[LIGHT_SUN ].front(), zvsize, zvsize);
		}
		if (has_moon) {
			smask[LIGHT_MOON].resize(zvals.size());
			calc_mesh_shadows(LIGHT_MOON, moon_pos, &zvals.front(), &smask[LIGHT_MOON].front(), zvsize, zvsize);
		}
		for (unsigned y = 0; y <= size; ++y) {
			for (unsigned x = 0; x <= size; ++x) {
				unsigned const ix(y*zvsize + x);
				float light_scale(1.0);

				if (has_sun && has_moon) {
					bool const no_sun( (smask[LIGHT_SUN ][ix] & SHADOWED_ALL) != 0);
					bool const no_moon((smask[LIGHT_MOON][ix] & SHADOWED_ALL) != 0);
					light_scale = blend_light(!no_sun, !no_moon);
				}
				else if (smask[has_sun ? LIGHT_SUN : LIGHT_MOON][ix] & SHADOWED_ALL) {
					light_scale = 0.0;
				}
				point const v((xstart + x*xstep), (ystart + y*ystep), zvals[ix]);
				data[y*stride + x].assign(v, get_norm(ix)*light_scale);
			}
		}
		for (unsigned y = 0; y < size; ++y) {
			for (unsigned x = 0; x < size; ++x) {
				unsigned const vix(y*stride + x), iix(4*(y*size + x));
				indices[iix+0] = vix;
				indices[iix+1] = vix + stride;
				indices[iix+2] = vix + stride + 1;
				indices[iix+3] = vix + 1;
			}
		}
		//PRINT_TIME("Create Data");
	}

	void create_texture(unsigned tex_bs) {
		assert(tid == 0);
		assert(!island);
		assert(zvals.size() == zvsize*zvsize);
		RESET_TIME;
		unsigned const tsize(base_tsize >> tex_bs), scale(tsize/size);
		float const dz_inv(1.0/(zmax - zmin)), fscale_inv(1.0/scale);
		unsigned char *data(new unsigned char[3*tsize*tsize]); // RGB
		int k1, k2, k3, k4;
		float t;

		for (unsigned y = 0; y < size; ++y) { // makes a big performance improvement
			for (unsigned x = 0; x < size; ++x) {
				unsigned const ix(y*zvsize + x);
				float const vnz00(get_norm(ix).z), vnz01(get_norm(ix+1).z);
				float const vnz10(get_norm(ix+zvsize).z), vnz11(get_norm(ix+zvsize+1).z);
				float const mh00(zvals[ix]), mh01(zvals[ix+1]), mh10(zvals[ix+zvsize]), mh11(zvals[ix+zvsize+1]);
				float const dist(fabs(mh01 - mh00) + fabs(mh00 - mh10) + fabs(mh01 - mh11) + fabs(mh10 - mh11));
				float const relh1((min(min(mh00, mh01), min(mh10, mh11)) - zmin)*dz_inv);
				float const relh2((max(max(mh00, mh01), max(mh10, mh11)) - zmin)*dz_inv);
				get_tids(relh1, NTEX_DIRT-1, h_dirt, k1, k2, t);
				get_tids(relh2, NTEX_DIRT-1, h_dirt, k3, k4, t);
				bool const same_tid(k1 == k4);
				k2 = k4;
				
				for (unsigned sy = 0; sy < scale; ++sy) {
					unsigned const ty(y*scale + sy);
					float const ypi(sy*fscale_inv);

					for (unsigned sx = 0; sx < scale; ++sx) {
						unsigned const tx(x*scale + sx), ix(y*zvsize + x), off(3*(ty*tsize + tx));
						float const xpi(sx*fscale_inv);
						float const mh((1.0 - xpi)*((1.0 - ypi)*mh00 + ypi*mh10) + xpi*((1.0 - ypi)*mh01 + ypi*mh11));

						if (!same_tid) {
							float const relh((mh - zmin)*dz_inv);
							get_tids(relh, NTEX_DIRT-1, h_dirt, k1, k2, t);
						}
						int const id(lttex_dirt[k1].id), id2(lttex_dirt[k2].id);
						texture const &t1(textures[id]);
						int const tof(t1.ncolors*(((ty<<tex_bs)&(t1.height-1))*t1.width + ((tx<<tex_bs)&(t1.width-1))));
						unsigned char *td(data + off);
						
						if (k1 == k2) { // single texture
							RGB_BLOCK_COPY(td, (t1.data + tof));
						}
						else { // blend two textures - performance critical
							texture const &t2(textures[id2]);
							int const tof2(t2.ncolors*(((ty<<tex_bs)&(t2.height-1))*t2.width + ((tx<<tex_bs)&(t2.width-1))));
							BLEND_COLOR(td, (t2.data + tof2), (t1.data + tof), t);
						}

						// handle steep slopes (dirt/rock texture replaces grass texture)
						float const sthresh[2] = {0.45, 0.7};
						float const vnz((1.0 - xpi)*((1.0 - ypi)*vnz00 + ypi*vnz10) + xpi*((1.0 - ypi)*vnz01 + ypi*vnz11));

						if (vnz < sthresh[1]) {
							if (id == GROUND_TEX || id2 == GROUND_TEX) { // ground/grass
								texture const &ta(textures[DIRT_TEX]);
								int const tofa(ta.ncolors*(((ty<<tex_bs)&(ta.height-1))*ta.width + ((tx<<tex_bs)&(ta.width-1))));
								unsigned char temp[3];

								if (id == GROUND_TEX || id2 == ROCK_TEX) {
									texture const &tb(textures[ROCK_TEX]);
									int const tofb(tb.ncolors*(((ty<<tex_bs)&(tb.height-1))*tb.width + ((tx<<tex_bs)&(tb.width-1))));
									BLEND_COLOR(temp, (tb.data+tofb), (ta.data+tofa), t);
								}
								else {
									RGB_BLOCK_COPY(temp, (ta.data+tofa));
								}
								float const val((vnz - sthresh[0])/(sthresh[1] - sthresh[0]));
								BLEND_COLOR(td, td, temp, CLIP_TO_01(val));
							}
							else if (id2 == SNOW_TEX) { // snow
								texture const &ta(textures[ROCK_TEX]);
								int const tofa(ta.ncolors*(((ty<<tex_bs)&(ta.height-1))*ta.width + ((tx<<tex_bs)&(ta.width-1))));
								float const val(2.0*(vnz - sthresh[0])/(sthresh[1] - sthresh[0]));
								BLEND_COLOR(td, td, (ta.data+tofa), CLIP_TO_01(val));
							}
						}

						// darken underwater regions
						if (!DISABLE_WATER && mh < water_plane_z) {
							float c[3] = {td[0], td[1], td[2]};
							atten_by_water_depth(c, get_water_atten_factor(mh));
							UNROLL_3X(td[i_] = (unsigned char)c[i_];)
						}
					} // for sx
				} // for sy
			} // for x
		} // for y
		//if (tex_bs == 0) {PRINT_TIME("Texture Gen");}
		bool const mipmaps(0); // mipmaps are too slow to build, not sure what wrap/mirror mode is best
		setup_texture(tid, GL_MODULATE, GL_LINEAR, mipmaps, 1, 1, 1, 1);
		assert(tid > 0);
		assert(glIsTexture(tid));

		if (mipmaps) {
			gluBuild2DMipmaps(GL_TEXTURE_2D, 3, tsize, tsize, GL_RGB, GL_UNSIGNED_BYTE, data);
		}
		else {
			//bool const has_comp(has_extension("GL_ARB_texture_compression")); GL_COMPRESSED_RGB - too slow
			glTexImage2D(GL_TEXTURE_2D, 0, 3, tsize, tsize, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
		}
		glDisable(GL_TEXTURE_2D);
		delete [] data;
		//if (tex_bs == 0) {PRINT_TIME("Texture Upload");}
	}

	void check_texture() {
		float dist(get_rel_dist_to_camera()*get_tile_radius()); // in tiles
		unsigned tex_bs(0);

		while (dist > 1.0 && (base_tsize >> tex_bs) > size) { // texture LOD
			dist /= 2;
			++tex_bs;
		}
		unsigned const new_tsize(base_tsize >> tex_bs);

		//if (new_tsize > gen_tsize) {
		//if (new_tsize != gen_tsize) { // power of 2 >= size
		if (new_tsize > gen_tsize || 2*new_tsize < gen_tsize) {
			clear_tid();
			gen_tsize = new_tsize; // power of 2 >= size
			create_texture(tex_bs);
		}
	}

	float get_rel_dist_to_camera() const {
		return (p2p_dist_xy(get_camera_pos(), get_center()) - radius)/(get_tile_radius()*(X_SCENE_SIZE + Y_SCENE_SIZE));
	}

	bool update_range() { // if returns 0, tile will be deleted
		float const dist(get_rel_dist_to_camera());
		if (dist > 1.5) clear_vbo_tid(1,1);
		return (dist < 2.0);
	}

	bool is_visible() const {
		return camera_pdu.sphere_visible_test(get_center(), radius);
	}

	void bind_vbos() {
		assert(vbo > 0 && ivbo > 0);
		bind_vbo(vbo,  0);
		bind_vbo(ivbo, 1);
	}

	bool draw(vector<vert_norm> &data, vector<unsigned short> &indices, vector<unsigned char> smask[]) { // make const or make vbo mutable?
		assert(size > 0);
		if (!is_visible()) return 0; // not visible to camera

		if (!DISABLE_TEXTURES) {
			check_texture();
			
			if (tid > 0) {
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, tid);
			}
		}
		glPushMatrix();
		glTranslatef(((xoff - xoff2) - init_dxoff)*DX_VAL, ((yoff - yoff2) - init_dyoff)*DY_VAL, 0.0);
		if (tid > 0) set_landscape_texgen(1.0, (-x1 - init_dxoff), (-y1 - init_dyoff), MESH_X_SIZE, MESH_Y_SIZE);
		unsigned ptr_stride(sizeof(vert_norm));

		if (vbo == 0) {
			create_data(data, indices, smask);
			vbo  = create_vbo();
			ivbo = create_vbo();
			bind_vbos();
			upload_vbo_data(&data.front(),    data.size()*ptr_stride,                0);
			upload_vbo_data(&indices.front(), indices.size()*sizeof(unsigned short), 1);
		}
		else {
			bind_vbos();
		}
		glVertexPointer(3, GL_FLOAT, ptr_stride, 0);
		glNormalPointer(   GL_FLOAT, ptr_stride, (void *)sizeof(point));
		glDrawElements(GL_QUADS, 4*size*size, GL_UNSIGNED_SHORT, 0);
		bind_vbo(0, 0);
		bind_vbo(0, 1);
		glPopMatrix();
		if (tid > 0) glDisable(GL_TEXTURE_2D);
		return 1;
	}
};


class tile_draw_t {

	typedef map<tile_xy_pair, tile_t*> tile_map;
	tile_map tiles;

public:
	tile_draw_t() {
		assert(MESH_X_SIZE == MESH_Y_SIZE && X_SCENE_SIZE == Y_SCENE_SIZE);
	}
	~tile_draw_t() {clear();}

	void clear() {
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			delete i->second;
		}
		tiles.clear();
	}

	void update() {
		//RESET_TIME;
		assert(MESH_X_SIZE == MESH_Y_SIZE); // limitation, for now
		point const camera(get_camera_pos() - point((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0));
		int const tile_radius(int(1.5*get_tile_radius()) + 1);
		int const toffx(int(0.5*camera.x/X_SCENE_SIZE)), toffy(int(0.5*camera.y/Y_SCENE_SIZE));
		int const x1(-tile_radius + toffx), y1(-tile_radius + toffy);
		int const x2( tile_radius + toffx), y2( tile_radius + toffy);
		unsigned const init_tiles(tiles.size());
		vector<tile_xy_pair> to_erase;

		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) { // free old tiles
			if (!i->second->update_range()) {
				to_erase.push_back(i->first);
				delete i->second;
			}
		}
		for (vector<tile_xy_pair>::const_iterator i = to_erase.begin(); i != to_erase.end(); ++i) {
			tiles.erase(*i);
		}
		for (int y = y1; y <= y2; ++y ) { // create new tiles
			for (int x = x1; x <= x2; ++x ) {
				tile_xy_pair const txy(x, y);
				if (tiles.find(txy) != tiles.end()) continue; // already exists
				tile_t tile(MESH_X_SIZE, x, y);
				if (tile.get_rel_dist_to_camera() >= 1.5) continue; // too far away to create
				tile_t *new_tile(new tile_t(tile));
				new_tile->create_zvals();
				tiles[txy] = new_tile;
			}
		}
		if (DEBUG_TILES && (tiles.size() != init_tiles || !to_erase.empty())) {
			cout << "update: tiles: " << init_tiles << " to " << tiles.size() << ", erased: " << to_erase.size() << endl;
		}
		//PRINT_TIME("Tiled Terrain Update");
	}

	float draw(bool add_hole) {
		bool const use_vbos(setup_gen_buffers_arb());
		assert(use_vbos);
		float zmin(FAR_CLIP);
		glDisable(GL_NORMALIZE);
		glDisable(GL_TEXTURE_COORD_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisable(GL_COLOR_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		setup_mesh_lighting();
		unsigned num_drawn(0);
		unsigned long long mem(0);
		vector<vert_norm> data;
		vector<unsigned short> indices;
		vector<unsigned char> smask[NUM_LIGHT_SRC];
		static point last_sun(all_zeros), last_moon(all_zeros);
		static bool last_water_en(1);
		bool const water_en((display_mode & 0x04) != 0);

		if (sun_pos != last_sun || moon_pos != last_moon) {
			clear_vbos_tids(1,0); // light source changed, clear vbos and build new shadow map
			last_sun  = sun_pos;
			last_moon = moon_pos;
		}
		if (water_en != last_water_en) { // should this be here?
			clear_vbos_tids(0,1); // water changed, recreate textures
			last_water_en = water_en;
		}
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			assert(i->second);
			if (DEBUG_TILES) mem += i->second->get_gpu_memory();
			if (add_hole && i->first.x == 0 && i->first.y == 0) continue;
			if (i->second->get_rel_dist_to_camera() > 1.4)      continue; // too far to draw
			zmin = min(zmin, i->second->get_zmin());
			num_drawn += i->second->draw(data, indices, smask);
		}
		if (DEBUG_TILES) cout << "tiles drawn: " << num_drawn << " of " << tiles.size() << ", gpu mem: " << mem/1024/1024 << endl;
		disable_textures_texgen();
		run_post_mesh_draw();
		return zmin;
	}

	void clear_vbos_tids(bool vclear, bool tclear) {
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			i->second->clear_vbo_tid(vclear, tclear);
		}
	}
};


tile_draw_t terrain_tile_draw;


void draw_vert_color(colorRGBA c, float x, float y, float z) {

	if (z < water_plane_z) atten_by_water_depth(&c.red, get_water_atten_factor(z));
	c.do_glColor();
	glVertex3f(x, y, z);
}


void fill_gap() {

	//RESET_TIME;
	colorRGBA const color(setup_mesh_lighting());
	set_landscape_texgen(1.0, xoff, yoff, MESH_X_SIZE, MESH_Y_SIZE);
	if (!DISABLE_TEXTURES) select_texture(LANDSCAPE_TEX);
	vector<float> xv(MESH_X_SIZE+1), yv(MESH_Y_SIZE+1);
	float const xstart(get_xval(xoff2)), ystart(get_yval(yoff2));

	for (int i = 0; i <= MESH_X_SIZE; ++i) { // not sure about this -x, +y thing
		xv[i] = (xstart + (i - 0.5)*DX_VAL);
	}
	for (int i = 0; i <= MESH_Y_SIZE; ++i) {
		yv[i] = (ystart + (i + 0.5)*DY_VAL);
	}

	// draw +x
	build_xy_mesh_arrays(&xv.front(), &yv[MESH_Y_SIZE], MESH_X_SIZE, 1);
	glBegin(GL_QUAD_STRIP);

	for (int x = 0; x < MESH_X_SIZE; ++x) {
		vertex_normals[MESH_Y_SIZE-1][x].do_glNormal();
		draw_vert_color(color, get_xval(x), Y_SCENE_SIZE-DY_VAL, mesh_height[MESH_Y_SIZE-1][x]);
		draw_vert_color(color, get_xval(x), Y_SCENE_SIZE       , fast_eval_from_index(x, 0, 0, 1));
	}
	glEnd();
	float const z_end_val(fast_eval_from_index(MESH_X_SIZE-1, 0, 0, 1));

	// draw +y
	build_xy_mesh_arrays(&xv[MESH_X_SIZE], &yv.front(), 1, MESH_Y_SIZE+1);
	glBegin(GL_QUAD_STRIP);

	for (int y = 0; y < MESH_X_SIZE; ++y) {
		vertex_normals[y][MESH_X_SIZE-1].do_glNormal();
		draw_vert_color(color, X_SCENE_SIZE-DX_VAL, get_yval(y), mesh_height[y][MESH_X_SIZE-1]);
		draw_vert_color(color, X_SCENE_SIZE       , get_yval(y), fast_eval_from_index(0, y, 0, 1));
	}

	// draw corner quad
	draw_vert_color(color, X_SCENE_SIZE-DX_VAL, Y_SCENE_SIZE, z_end_val);
	draw_vert_color(color, X_SCENE_SIZE       , Y_SCENE_SIZE, fast_eval_from_index(0, MESH_Y_SIZE, 0, 1));
	glEnd();

	disable_textures_texgen();
	glDisable(GL_COLOR_MATERIAL);
	run_post_mesh_draw();
	//PRINT_TIME("Fill Gap");
}


float draw_tiled_terrain(bool add_hole) {

	//RESET_TIME;
	terrain_tile_draw.update();
	float const zmin(terrain_tile_draw.draw(add_hole));
	if (add_hole) fill_gap(); // need to fill the gap on +x/+y
	//glFinish();
	//PRINT_TIME("Tiled Terrain Draw");
	return zmin;
}


void clear_tiled_terrain() {
	terrain_tile_draw.clear();
}

void reset_tiled_terrain_state() {
	terrain_tile_draw.clear_vbos_tids(1,1);
}

float get_inf_terrain_fog_dist() {
	return (tiled_mesh_display ? 2.5*get_tile_radius()*XY_SCENE_SIZE : VIEW_DIST0*XY_SCENE_SIZE);
}


