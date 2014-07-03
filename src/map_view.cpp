// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/2/02

#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "physics_objects.h"
#include "shaders.h"
#include "heightmap.h"


bool const MAP_VIEW_LIGHTING = 1;

int map_x(0), map_y(0);
float map_zoom(0.25);

extern int window_width, window_height, xoff2, yoff2, map_mode, map_color, begin_motion;
extern int world_mode, game_mode, display_mode, num_smileys, DISABLE_WATER, cache_counter;
extern float zmax_est, water_plane_z, glaciate_exp, glaciate_exp_inv, vegetation, relh_adj_tex;
extern int coll_id[];
extern obj_group obj_groups[];


void setup_height_gen(mesh_xy_grid_cache_t &height_gen, float x0, float y0, float dx, float dy, unsigned nx, unsigned ny, bool cache_values);
bool using_hmap_with_detail();


float get_mesh_height(mesh_xy_grid_cache_t const &height_gen, float xstart, float ystart, float scale, int i, int j) {

	if (using_tiled_terrain_hmap_tex()) {
		float zval(get_tiled_terrain_height_tex((xstart + X_SCENE_SIZE + j*scale)*DX_VAL_INV, (ystart + Y_SCENE_SIZE + i*scale)*DY_VAL_INV));
		if (using_hmap_with_detail()) {zval += HMAP_DETAIL_MAG*height_gen.eval_index(j, i, 0);}
		return zval;
	}
	return height_gen.eval_index(j, i, 1);
}


void draw_overhead_map() {

	//RESET_TIME
	unsigned tid(0);
	if (map_mode == 0) return;
	
	if (map_mode == 2) {
		map_mode = 1;
		return;
	}
	bool const no_water((DISABLE_WATER == 2) || !(display_mode & 0x04));
	int bx1(0), by1(0), bx2(0), by2(0), nx(1);
	while (min(window_width, window_height) > 2*nx) nx *= 2;
	if (nx < 4) return;
	int const ny(nx), nx2(nx/2), ny2(ny/2);
	float const zmax2(zmax_est*((map_color || no_water) ? 1.0 : 0.855)), hscale(0.5/zmax2);
	float const scale(2.0*map_zoom*X_SCENE_SIZE*DX_VAL), scale_val(scale/64);
	float x0(map_x + xoff2*DX_VAL), y0(map_y + yoff2*DY_VAL);
	float const relh_water(get_rel_height(water_plane_z, -zmax_est, zmax_est));
	point const camera(get_camera_pos());
	float map_heights[6];
	map_heights[0] = 0.9*lttex_dirt[3].zval  + 0.1*lttex_dirt[4].zval;
	map_heights[1] = 0.5*(lttex_dirt[2].zval + lttex_dirt[3].zval);
	map_heights[2] = 0.5*(lttex_dirt[1].zval + lttex_dirt[2].zval);
	map_heights[3] = 0.5*(lttex_dirt[0].zval + lttex_dirt[1].zval);
	map_heights[4] = relh_water;
	map_heights[5] = 0.5*map_heights[4];

	colorRGBA const map_colors[6] = {
		((DISABLE_WATER == 2) ? DK_GRAY : WHITE),
		GRAY,
		((vegetation == 0.0) ? colorRGBA(0.55,0.45,0.35,1.0) : GREEN),
		LT_BROWN,
		(no_water ? BROWN    : colorRGBA(0.3,0.2,0.6)),
		(no_water ? DK_BROWN : BLUE)};

	for (unsigned i = 0; i < 6; ++i) {
		map_heights[i] = pow(map_heights[i], glaciate_exp);
	}
	if (world_mode == WMODE_GROUND) {
		float const xv(-(camera.x + map_x)/X_SCENE_SIZE), yv(-(camera.y + map_y)/Y_SCENE_SIZE);
		float const xs(DX_VAL/scale_val), ys(DY_VAL/scale_val);
		x0 += camera.x;
		y0 += camera.y;
		bx1 = int(nx2 + xs*(xv - 1.0));
		by1 = int(ny2 + ys*(yv - 1.0));
		bx2 = int(nx2 + xs*(xv + 1.0));
		by2 = int(ny2 + ys*(yv + 1.0));
	}
	vector3d const dir(vector3d(cview_dir.x, cview_dir.y, 0.0).get_norm());
	int const cx(int(nx2 - map_x/scale)), cy(int(ny2 - map_y/scale));
	int const xx(cx + int(4*dir.x)), yy(cy + int(4*dir.y));
	float const xstart(x0 - nx2*scale), ystart(y0 - ny2*scale);
	mesh_xy_grid_cache_t height_gen;
	setup_height_gen(height_gen, xstart, ystart, scale, scale, nx, ny, 1); // cache_values=1
	vector<unsigned char> buf(nx*ny*3*sizeof(unsigned char));
	vector3d const light_dir(get_light_pos().get_norm()); // assume directional lighting to origin

	#pragma omp parallel for schedule(static,1)
	for (int i = 0; i < ny; ++i) {
		int const inx(i*nx), iyy((i - yy)*(i - yy)), icy((i - cy)*(i - cy));
		float last_height(0.0);

		for (int j = 0; j < nx; ++j) {
			int const offset(3*(inx + j));
			unsigned char *rgb(&buf[offset]);

			if (iyy + (j - xx)*(j - xx) <= 4) {
				rgb[0] = rgb[1] = rgb[2] = 0; // camera direction
			}
			else if (icy + (j - cx)*(j - cx) <= 9) {
				rgb[0] = 255;
				rgb[1] = rgb[2] = 0; // camera position
			}
			else if (world_mode == WMODE_GROUND && (((i == by1 || i == by2) && j >= bx1 && j < bx2) ||
				((j == bx1 || j == bx2) && i >= by1 && i < by2))) {
				rgb[0] = rgb[1] = rgb[2] = 0; // world boundary
			}
			else {
				float height(CLIP_TO_01(hscale*(get_mesh_height(height_gen, xstart, ystart, scale, i, j) + zmax2)));

				if (!map_color) { // grayscale
					rgb[0] = rgb[1] = rgb[2] = (unsigned char)(255.0*pow(height, glaciate_exp_inv)); // un-glaciate: slow
				}
				else {
					height += relh_adj_tex;

					if (height <= map_heights[5]) {
						unpack_color(rgb, map_colors[5]); // deep water
					}
					else if (height <= map_heights[3]) {
						unpack_color(rgb, map_colors[3]); // sand
					}
					else if (height >= map_heights[0]) {
						unpack_color(rgb, map_colors[0]); // snow
					}
					else {
						for (unsigned k = 0; k < 4; ++k) { // mixed
							if (height > map_heights[k+1]) {
								float const hval((height - map_heights[k+1])/(map_heights[k] - map_heights[k+1]));
								UNROLL_3X(rgb[i_] = (unsigned char)(255.0*(hval*map_colors[k][i_] + (1.0 - hval)*map_colors[k+1][i_]));)
								break;
							}
						}
					}
					if (height <= map_heights[4] && height > map_heights[5]) { // shallow water
						float const hval(0.5*(height - map_heights[5])/(map_heights[4] - map_heights[5]));
						UNROLL_3X(rgb[i_] = (unsigned char)(255.0*(1.0 - hval)*map_colors[5][i_] + hval*rgb[i_]);)
					}
					if (MAP_VIEW_LIGHTING) {
						vector3d normal(plus_z);

						if (height > map_heights[4]) {
							float const hx((j == 0) ? height : last_height);
							float const hy(CLIP_TO_01(hscale*(get_mesh_height(height_gen, xstart, ystart, scale, max(i-1, 0), j) + zmax2)));
							normal = vector3d(DY_VAL*(hx - height), DX_VAL*(hy - height), dxdy).get_norm();
						}
						last_height = height;
						float const light_val(0.2 + 0.8*max(0.0f, dot_product(light_dir, normal)));
						UNROLL_3X(rgb[i_] = (unsigned char)(light_val*rgb[i_]);)
					}
				}
			}
		}
	}
	if (begin_motion && obj_groups[coll_id[SMILEY]].enabled) { // game_mode?
		float const camx((world_mode == WMODE_GROUND) ? camera.x : 0.0), camy((world_mode == WMODE_GROUND) ? camera.y : 0.0);

		for (int s = 0; s < num_smileys; ++s) { // add in smiley markers
			point const spos(obj_groups[coll_id[SMILEY]].get_obj(s).pos);
			int const xpos(int(nx2 + ((-camx - map_x + spos.x)/X_SCENE_SIZE)*DX_VAL/scale_val));
			int const ypos(int(ny2 + ((-camy - map_y + spos.y)/Y_SCENE_SIZE)*DY_VAL/scale_val));
			colorRGBA const color(get_smiley_team_color(s));

			for (int i = max(0, ypos-1); i < min(ny, ypos+1); ++i) {
				for (int j = max(0, xpos-1); j < min(nx, xpos+1); ++j) {
					int const offset(3*(i*nx + j));
					unpack_color(&buf[offset], color);
				}
			}
		}
	}
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	shader_t s;
	s.begin_simple_textured_shader(0.0, 0, 0, &WHITE);
	setup_texture(tid, 0, 0, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, &buf.front());
	draw_tquad(0.58*((float)window_width)/((float)window_height), 0.58, -1.0);
	free_texture(tid);
	s.end_shader();
	//PRINT_TIME("draw map")
}
