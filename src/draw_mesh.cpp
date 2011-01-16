// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/27/02

#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"


float const W_TEX_SCALE0     = 1.0;
float const WATER_WIND_EFF   = 0.0006;
float const START_OFFSET0    = 0.05;
float const PRESP_ANGLE_ADJ  = 1.5;
float const VD_SCALE         = 1.0;
float const SURF_HEAL_RATE   = 0.005;
float const MAX_SURFD        = 20.0;
float const DLIGHT_SCALE     = 4.0;
unsigned const W_STEPS       = 16;
int   const DRAW_BORDER      = 3;

int   const SHOW_MESH_TIME   = 0;
int   const SHOW_NORMALS     = 0;
int   const DEBUG_COLLS      = 0; // 0 = disabled, 1 = lines, 2 = cubes
int   const DISABLE_TEXTURES = 0;


struct fp_ratio {
	float n, d;
	inline float get_val() const {return ((d == 0.0) ? 1.0 : n/d);}
};


// Global Variables
bool clear_landscape_vbo;
int island(0);
float lt_green_int(1.0), sm_green_int(1.0);
vector<fp_ratio> uw_mesh_lighting; // for water caustics

extern bool using_lightmap, has_dl_sources, combined_gu, has_snow;
extern unsigned num_jterms;
extern int draw_model, num_local_minima, world_mode, xoff, yoff, xoff2, yoff2, ocean_set, ground_effects_level, animate2;
extern int display_mode, frame_counter, resolution, verbose_mode, DISABLE_WATER, read_landscape, disable_inf_terrain;
extern float zmax, zmin, zmax_est, ztop, zbottom, light_factor, max_water_height, init_temperature;
extern float water_plane_z, temperature, fticks, mesh_scale, mesh_z_cutoff, TWO_XSS, TWO_YSS, XY_SCENE_SIZE;
extern point light_pos, litning_pos, sun_pos, moon_pos;
extern vector3d up_norm, wind;
extern float h_dirt[];


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

	set_array_client_state(1, 0, 1, (carr != 0));
	if (varr) glVertexPointer(3, GL_FLOAT, 0, varr);
	if (narr) glNormalPointer(   GL_FLOAT, 0, narr);
	if (carr) glColorPointer( 3, GL_FLOAT, 0, carr);
}


void run_post_mesh_draw() {
	
	glEnable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	disable_blend();
}


float integrate_water_dist(point const &targ_pos, point const &src_pos, float const water_z) {

	if (src_pos.z == targ_pos.z) return 0.0;
	float const t(min(1.0f, (water_z - targ_pos.z)/fabs(src_pos.z - targ_pos.z))); // min(1.0,...) for underwater case
	point p_int(targ_pos + (src_pos - targ_pos)*t);
	int const xp(get_xpos(targ_pos.x)), yp(get_ypos(targ_pos.y));
	if (!point_outside_mesh(xp, yp)) p_int.z = min(src_pos.z, water_matrix[yp][xp]); // account for ripples
	return p2p_dist(p_int, targ_pos)*mesh_scale;

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

		if (!(display_mode & 0x08) && DLIGHT_SCALE > 0.0 && (using_lightmap || has_dl_sources)) { // somewhat slow
			get_sd_light(j, i, get_zpos(varr[c].z), varr[c], DLIGHT_SCALE, &carr[c].R, &surface_normals[i][j], NULL);
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
			light_scale = blend_light(light_factor, !no_sun, !no_moon);
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
		// Note: normal is never set to zero because we need it for dynamic light sources
		narr[c] = vertex_normals[i][j]*max(light_scale, 0.01f);
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
			bool const refracted(calc_refraction_angle(dir, v_refract, wat_vert_normals[y][x], 1.0, WATER_INDEX_REFRACT));
			assert(refracted); // can't have total internal reflection going into the water if the physics are sane
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
				assert(ysz >= 0.0);
				if (ysz <= 0.0) continue;
				crng[0][0] = init_cr[0];
				crng[1][0] = init_cr[0] + DX_VAL;

				for (int xx = bnds[0][0]; xx < bnds[1][0]; ++xx) {
					assert(xx >= 0 && xx < MESH_X_SIZE);
					float const xsz(min(crng[1][0], rng[1][0]) - max(crng[0][0], rng[0][0])); // intersection: min(UB) - max(LB)
					assert(xsz >= 0.0);
					if (xsz <= 0.0) continue;
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


void set_landscape_texgen(float tex_scale, int xoffset, int yoffset, int xsize, int ysize, bool use_detail_tex) {

	float const tx(tex_scale*(((float)xoffset)/((float)xsize) + 0.5));
	float const ty(tex_scale*(((float)yoffset)/((float)ysize) + 0.5));
	setup_texgen(tex_scale/TWO_XSS, tex_scale/TWO_YSS, tx, ty);

	if (use_detail_tex) { // blend in detail nose texture at 30x scale
		select_multitex(NOISE_TEX, 1);
		setup_texgen(30.0*tex_scale/TWO_XSS, 30.0*tex_scale/TWO_YSS, 0.0, 0.0);
	}
}


void draw_coll_vert(int i, int j) {

	(v_collision_matrix[i][j].cvz.empty() ? BLUE : RED).do_glColor();
	glVertex3f(get_xval(j), get_yval(i), max(czmin, v_collision_matrix[i][j].zmax));
}


void display_mesh() { // fast array version

	if (mesh_height == NULL) return; // no mesh to display
	// can't put the hole in the right place, so only draw the tiled terrain
	if (!island && (display_mode & 0x10) && (xoff2 != 0 || yoff2 != 0)) return;
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

		if (DEBUG_COLLS == 2) {
			enable_blend();
			colorRGBA(1.0, 0.0, 0.0, 0.1).do_glColor();

			for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
				for (int j = 0; j < MESH_X_SIZE; ++j) {
					if (v_collision_matrix[i][j].zmin < v_collision_matrix[i][j].zmax) {
						point const p1(get_xval(j+0), get_yval(i+0),v_collision_matrix[i][j].zmin);
						point const p2(get_xval(j+1), get_yval(i+1),v_collision_matrix[i][j].zmax);
						draw_cube((p1 + p2)*0.5, (p2.x - p1.x), (p2.y - p1.y), (p2.z - p1.z), 0, 1);
					}
				}
			}
			disable_blend();
		}
		else {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

			for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
				glBegin(GL_TRIANGLE_STRIP);
				
				for (int j = 0; j < MESH_X_SIZE; ++j) {
					draw_coll_vert(i+0, j);
					draw_coll_vert(i+1, j);
				}
				glEnd();
			}
		}
		glEnable(GL_LIGHTING);
	}
	setup_mesh_lighting();
	update_landscape_texture();
	if (SHOW_MESH_TIME) PRINT_TIME("Landscape Texture");
	glDisable(GL_NORMALIZE);

	if (!DISABLE_TEXTURES) {
		select_texture(LANDSCAPE_TEX);
		set_landscape_texgen(1.0, xoff, yoff, MESH_X_SIZE, MESH_Y_SIZE);
	}
	if (SHOW_MESH_TIME) PRINT_TIME("Preprocess");

	if (ground_effects_level == 0 && setup_gen_buffers()) { // simpler, more efficient mesh draw
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
		set_array_client_state(1, 0, 1, 0);
		glVertexPointer(3, GL_FLOAT, 2*sizeof(point), 0);
		glNormalPointer(   GL_FLOAT, 2*sizeof(point), (void *)sizeof(point));

		for (int i = 0; i < MESH_Y_SIZE-1; ++i) {
			glDrawArrays(GL_QUAD_STRIP, 2*i*MESH_X_SIZE, 2*MESH_X_SIZE);
		}
		bind_vbo(0);
	}
	else { // slower mesh draw with more features
		bool const use_shaders((display_mode & 0x08) != 0);

		if (use_shaders) {
			set_shader_prefix("#define USE_LIGHT_COLORS", 0); // VS
			setup_enabled_lights();
			set_dlights_booleans(1, 1); // FS
			unsigned const p(set_shader_prog("ads_lighting.part*+fog.part+texture_gen.part+draw_mesh", "dynamic_lighting.part*+linear_fog.part+draw_mesh"));
			setup_fog_scale(p);
			add_uniform_int(p, "tex0", 0);
			add_uniform_int(p, "tex1", 1);
			setup_scene_bounds(p);
			setup_dlight_textures(p);
		}
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
		if (use_shaders) unset_shader_prog();
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


float display_mesh3(int const *const hole_bounds) { // WM3 - infinite terrain

	bool const add_hole((hole_bounds != NULL) && xoff2 == 0 && yoff2 == 0);
	return draw_tiled_terrain(add_hole);
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
	bool in_quads(0);
	glNormal3f(-1.0, 0.0, 0.0);

	for (int i = 1; i < MESH_Y_SIZE; ++i) { // x sides
		float const mh1(mesh_height[i][end_val]), mh2(mesh_height[i-1][end_val]);
		float const wm1(water_matrix[i][end_val] - SMALL_NUMBER), wm2(water_matrix[i-1][end_val] - SMALL_NUMBER);

		if (!check_zvals || mh1 < wm1 || mh2 < wm2) {
			if (!in_quads) {glBegin(GL_QUADS); in_quads = 1;}
			draw_vertex(limit, yv,        wm2,           1, tex_scale);
			draw_vertex(limit, yv+DY_VAL, wm1,           1, tex_scale);
			draw_vertex(limit, yv+DY_VAL, min(wm1, mh1), 1, tex_scale);
			draw_vertex(limit, yv,        min(wm2, mh2), 1, tex_scale);
		}
		yv += DY_VAL;
	}
	if (in_quads) glEnd();
}


void water_renderer::draw_y_sides(bool neg_edge) const {

	int const end_val(neg_edge ? 0 : MESH_Y_SIZE-1);
	float const limit(neg_edge ? -Y_SCENE_SIZE : Y_SCENE_SIZE-DY_VAL);
	float xv(-X_SCENE_SIZE);
	bool in_quads(0);
	glNormal3f(0.0, 1.0, 0.0);
	
	for (int i = 1; i < MESH_X_SIZE; ++i) { // y sides
		float const mh1(mesh_height[end_val][i]), mh2(mesh_height[end_val][i-1]);
		float const wm1(water_matrix[end_val][i] - SMALL_NUMBER), wm2(water_matrix[end_val][i-1] - SMALL_NUMBER);

		if (!check_zvals || mh1 < wm1 || mh2 < wm2) {
			if (!in_quads) {glBegin(GL_QUADS); in_quads = 1;}
			draw_vertex(xv,        limit, wm2,           0, tex_scale);
			draw_vertex(xv+DX_VAL, limit, wm1,           0, tex_scale);
			draw_vertex(xv+DX_VAL, limit, min(wm1, mh1), 0, tex_scale);
			draw_vertex(xv,        limit, min(wm2, mh2), 0, tex_scale);
		}
		xv += DX_VAL;
	}
	if (in_quads) glEnd();
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

	if (!disable_inf_terrain && (display_mode & 0x10)) return;
	water_renderer wr(resolution, resolution, check_zvals);
	wr.draw();
}


void draw_water_plane(float zval, int const *const hole_bounds) {

	if (DISABLE_WATER) return;

	if (display_mode & 0x0100) { // add small waves
		static float time(0.0);
		if (animate2) time += fticks;
		zval += 0.01*sin(1.0*time/TICKS_PER_SECOND);
	}
	float const tscale(W_TEX_SCALE0/Z_SCENE_SIZE), vd_scale(2.0*get_tile_radius()*SQRT2);
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
	point const camera(get_camera_pos());
	vector3d(0.0, 0.0, ((camera.z < zval) ? -1.0 : 1.0)).do_glNormal();
	set_color(color);
	set_fill_mode();
	enable_blend();
	setup_texgen(tscale, tscale, (tscale*(xoff2 - xoff)*DX_VAL + wxoff), (tscale*(yoff2 - yoff)*DY_VAL + wyoff));

	//setup_enabled_lights();
	//unsigned const p(set_shader_prog("texture_gen.part+per_pixel_lighting_texgen", "ads_lighting.part*+per_pixel_lighting_textured")); // needs fog
	//add_uniform_int(p, "tex0", 0);
	
	glPushMatrix();
	glTranslatef(0.0, 0.0, zval);
	glBegin(GL_QUADS);

	if (hole_bounds) { // x1 x2 y1 y2
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
		}
	}
	else {
		float const xinc(2.0*vdx/(float)W_STEPS); // about 9
		float const yinc(2.0*vdy/(float)W_STEPS);
		float yval(dy - vdy);

		for (unsigned i = 0; i < W_STEPS; ++i) {
			float xval(dx - vdx);

			for (unsigned j = 0; j < W_STEPS; ++j) {
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
	//unset_shader_prog();
	disable_blend();
	set_specular(0.0, 1.0);
	disable_textures_texgen();
	glEnable(GL_LIGHTING);
}


float get_inf_terrain_fog_dist() {
	return 3.0*get_tile_radius()*XY_SCENE_SIZE;
}


