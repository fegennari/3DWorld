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
bool const MAP_VIEW_SHADOWS  = 1;

int map_drag_x(0), map_drag_y(0);
float map_x(0.0), map_y(0.0), map_zoom(0.0);

extern bool water_is_lava;
extern int window_width, window_height, xoff2, yoff2, map_mode, map_color, begin_motion, read_landscape, read_heightmap, do_read_mesh;
extern int world_mode, game_mode, display_mode, num_smileys, DISABLE_WATER, cache_counter, default_ground_tex;
extern float zmax_est, water_plane_z, glaciate_exp, glaciate_exp_inv, vegetation, relh_adj_tex, temperature;
extern int coll_id[];
extern obj_group obj_groups[];
extern coll_obj_group coll_objects;


void setup_height_gen(mesh_xy_grid_cache_t &height_gen, float x0, float y0, float dx, float dy, unsigned nx, unsigned ny, bool cache_values);
bool using_hmap_with_detail();
void set_temp_clear_color(colorRGBA const &clear_color);


float get_mesh_height(mesh_xy_grid_cache_t const &height_gen, float xstart, float ystart, float xscale, float yscale, int i, int j) {

	if (using_tiled_terrain_hmap_tex()) {
		float zval(get_tiled_terrain_height_tex((xstart + X_SCENE_SIZE + j*xscale)*DX_VAL_INV, (ystart + Y_SCENE_SIZE + i*yscale)*DY_VAL_INV));
		if (using_hmap_with_detail()) {zval += HMAP_DETAIL_MAG*height_gen.eval_index(j, i, 0);}
		return zval;
	}
	return height_gen.eval_index(j, i, 1);
}

bool is_shadowed(point const &cpos, vector3d const &cnorm, point const &lpos, int &cindex) {

	if (!MAP_VIEW_SHADOWS)   return 0;
	if (display_mode & 0x20) return 0;
	point const cpos2(cpos + 0.001*cnorm);
	if (cindex >= 0 && coll_objects.get_cobj(cindex).line_intersect(cpos2, lpos)) return 1;
	return check_coll_line(cpos2, lpos, cindex, -1, 1, 3); // static cobj shadows only for performance
}


void draw_overhead_map() {

	//RESET_TIME
	unsigned tid(0);
	if (map_mode == 0) return;
	
	if (map_mode == 2) {
		map_mode = 1;
		return;
	}
	if (map_zoom == 0.0) {map_zoom = ((world_mode == WMODE_GROUND) ? 0.08 : 0.8);} // set reasonable defaults based on mode
	int bx1(0), by1(0), bx2(0), by2(0), nx(1);
	while (min(window_width, window_height) > 2*nx) {nx *= 2;}
	if (nx < 4) return;

	int const ny(nx), nx2(nx/2), ny2(ny/2);
	bool const no_water((DISABLE_WATER == 2) || !(display_mode & 0x04));
	bool const is_ice(((world_mode == WMODE_GROUND) ? temperature : get_cur_temperature()) <= W_FREEZE_POINT);
	float const zmax2(zmax_est*((map_color || no_water) ? 1.0 : 0.855)), hscale(0.5/zmax2);
	float const window_ar(float(window_width)/float(window_height)), scene_ar(X_SCENE_SIZE/Y_SCENE_SIZE);
	float const xscale(2.0*map_zoom*window_ar*HALF_DXY), yscale(2.0*map_zoom*scene_ar*HALF_DXY);
	float const xscale_val(xscale/64), yscale_val(yscale/64);

	// translate map_drag_x/y (screen pixel space) into map_x/y (world unit space)
	float const x_scale(nx*xscale/window_width), y_scale(ny*yscale/window_height);
	map_x += x_scale*map_drag_x; map_drag_x = 0;
	map_y += y_scale*map_drag_y; map_drag_y = 0;

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
	for (unsigned i = 0; i < 6; ++i) {map_heights[i] = pow(map_heights[i], glaciate_exp);}
	colorRGBA ground_color(BLACK);
	if (default_ground_tex >= 0) {ground_color = texture_color(default_ground_tex);}

	colorRGBA const map_colors[6] = {
		((water_is_lava || DISABLE_WATER == 2) ? DK_GRAY : WHITE),
		GRAY,
		((vegetation == 0.0) ? colorRGBA(0.55,0.45,0.35,1.0) : GREEN),
		LT_BROWN,
		(no_water ? BROWN    : (water_is_lava ? RED        : colorRGBA(0.3,0.2,0.6))),
		(no_water ? DK_BROWN : (water_is_lava ? LAVA_COLOR : (is_ice ? LT_BLUE : BLUE)))};

	if (world_mode == WMODE_GROUND) {
		float const xv(-(camera.x + map_x)/X_SCENE_SIZE), yv(-(camera.y + map_y)/Y_SCENE_SIZE);
		float const xs(DX_VAL/xscale_val), ys(DY_VAL/yscale_val);
		x0 += camera.x;
		y0 += camera.y;
		bx1 = int(nx2 + xs*(xv - 1.0));
		by1 = int(ny2 + ys*(yv - 1.0));
		bx2 = int(nx2 + xs*(xv + 1.0));
		by2 = int(ny2 + ys*(yv + 1.0));
	}
	vector3d const dir(vector3d(cview_dir.x, cview_dir.y, 0.0).get_norm());
	int const cx(int(nx2 - map_x/xscale)), cy(int(ny2 - map_y/yscale));
	int const xx(cx + int(4*dir.x)), yy(cy + int(4*dir.y));
	float const xstart(x0 - nx2*xscale), ystart(y0 - ny2*yscale);

	bool const uses_hmap(world_mode == WMODE_GROUND && (read_landscape || read_heightmap || do_read_mesh));
	mesh_xy_grid_cache_t height_gen;
	if (!uses_hmap) {setup_height_gen(height_gen, xstart, ystart, xscale, yscale, nx, ny, 1);} // cache_values=1
	vector<unsigned char> buf(nx*ny*3*sizeof(unsigned char));
	point const lpos(get_light_pos());
	vector3d const light_dir(lpos.get_norm()); // assume directional lighting to origin

	#pragma omp parallel for schedule(static,1)
	for (int i = 0; i < ny; ++i) {
		int const inx(i*nx), iyy((i - yy)*(i - yy)), icy((i - cy)*(i - cy));
		float last_height(0.0);
		point cpos;
		vector3d cnorm;
		int cindex(-1), cindex2(-1);

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
			else if (world_mode == WMODE_GROUND &&
				(((i == by1 || i == by2) && j >= bx1 && j < bx2) || ((j == bx1 || j == bx2) && i >= by1 && i < by2)))
			{
				rgb[0] = rgb[1] = rgb[2] = 0; // world boundary
			}
			else {
				float mh(0.0);
				bool mh_set(0), shadowed(0);

				if (world_mode == WMODE_GROUND) {
					float const xval((j - nx2)*xscale_val*(X_SCENE_SIZE/DX_VAL) + camera.x + map_x);
					float const yval((i - ny2)*yscale_val*(Y_SCENE_SIZE/DY_VAL) + camera.y + map_y);
					point p1(xval, yval, czmax);
					bool const over_mesh(is_over_mesh(p1));
					
					if (over_mesh || uses_hmap) { // if using a heightmap, clamp values to scene bounds
						mh = interpolate_mesh_zval(max(-X_SCENE_SIZE, min(X_SCENE_SIZE-DX_VAL, xval)), max(-Y_SCENE_SIZE, min(Y_SCENE_SIZE-DY_VAL, yval)), 0.0, 0, 1);
						mh_set = 1;
					}
					if (over_mesh && czmin < czmax) { // check cobjs
						point p2(xval, yval, max(mh, czmin));
						float t;
						int cindex0(-1);
						if (cindex >= 0 && coll_objects.get_cobj(cindex).line_int_exact(p1, p2, t, cnorm)) {cpos = p1 + t*(p2 - p1); p2 = cpos;} // previous cobj int
						else {cindex = -1;} // else reset
						if (check_coll_line_exact(p1, p2, cpos, cnorm, cindex0, 0.0, cindex, 1, 0, 0, 0, 0)) {cindex = cindex0;} // cobj intersection

						if (cindex >= 0) {
							colorRGBA const color(get_cobj_color_at_point(cindex, cpos, cnorm, 0));
							unpack_color(rgb, color*(is_shadowed(cpos, cnorm, lpos, cindex2) ? 0.5 : 1.0));
							continue;
						}
						if (mh_set) {shadowed = is_shadowed(point(xval, yval, mh), plus_z, lpos, cindex2);}
					}
				}
				if (default_ground_tex >= 0 && map_color) {
					unpack_color(rgb, ground_color*(shadowed ? 0.5 : 1.0));
					continue;
				}
				if (!mh_set) {mh = get_mesh_height(height_gen, xstart, ystart, xscale, yscale, i, j);} // calculate mesh height here if not yet set
				float height(CLIP_TO_01(hscale*(mh + zmax2)));

				if (!map_color) { // grayscale
					rgb[0] = rgb[1] = rgb[2] = (unsigned char)(255.0*pow(height, glaciate_exp_inv)); // un-glaciate: slow
				}
				else {
					height += relh_adj_tex;
					colorRGBA color;
					if      (height <= map_heights[5]) {color = map_colors[5];} // deep water
					else if (height <= map_heights[3]) {color = map_colors[3];} // sand
					else if (height >= map_heights[0]) {color = map_colors[0];} // snow
					else {
						for (unsigned k = 0; k < 4; ++k) { // mixed
							if (height > map_heights[k+1]) {
								float const h((height - map_heights[k+1])/(map_heights[k] - map_heights[k+1])), v(cubic_interpolate(h));
								blend_color(color, map_colors[k], map_colors[k+1], v);
								break;
							}
						}
					}
					if (height <= map_heights[4] && height > map_heights[5]) { // shallow water
						float const h(0.5*(height - map_heights[5])/(map_heights[4] - map_heights[5])), v(cubic_interpolate(h));
						blend_color(color, color, map_colors[5], v);
					}
					if (MAP_VIEW_LIGHTING && !uses_hmap && !(display_mode & 0x20)) {
						vector3d normal(plus_z);

						if (height > map_heights[4]) {
							float const hx((j == 0) ? height : last_height);
							float const hy(CLIP_TO_01(hscale*(get_mesh_height(height_gen, xstart, ystart, xscale, yscale, max(i-1, 0), j) + zmax2)));
							normal = vector3d(DY_VAL*(hx - height), DX_VAL*(hy - height), dxdy).get_norm();
						}
						last_height = height;
						color *= (0.2 + (shadowed ? 0.0 : 0.8)*max(0.0f, dot_product(light_dir, normal)));
						shadowed = 0; // handled correctly above
					}
					unpack_color(rgb, color*(shadowed ? 0.5 : 1.0));
				}
			}
		} // for j
	} // for i
	if (begin_motion && obj_groups[coll_id[SMILEY]].enabled) { // game_mode?
		float const camx((world_mode == WMODE_GROUND) ? camera.x : 0.0), camy((world_mode == WMODE_GROUND) ? camera.y : 0.0);

		for (int s = 0; s < num_smileys; ++s) { // add in smiley markers
			point const spos(obj_groups[coll_id[SMILEY]].get_obj(s).pos);
			int const xpos(int(nx2 + ((-camx - map_x + spos.x)/X_SCENE_SIZE)*DX_VAL/xscale_val));
			int const ypos(int(ny2 + ((-camy - map_y + spos.y)/Y_SCENE_SIZE)*DY_VAL/yscale_val));
			colorRGBA const color(get_smiley_team_color(s));

			for (int i = max(0, ypos-1); i < min(ny, ypos+1); ++i) {
				for (int j = max(0, xpos-1); j < min(nx, xpos+1); ++j) {
					int const offset(3*(i*nx + j));
					unpack_color(&buf[offset], color);
				}
			}
		}
	}
	set_temp_clear_color(BLACK);
	shader_t s;
	s.begin_simple_textured_shader(0.0, 0, 0, &WHITE);
	setup_texture(tid, 0, 0, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, &buf.front());
	draw_tquad(0.58*((float)window_width)/((float)window_height), 0.58, -1.0);
	free_texture(tid);
	s.end_shader();
	//PRINT_TIME("draw map")
}
