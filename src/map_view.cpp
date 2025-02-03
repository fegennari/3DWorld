// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/2/02

#include "3DWorld.h"
#include "mesh.h"
#include "textures.h"
#include "physics_objects.h"
#include "shaders.h"
#include "heightmap.h"
#include <cfloat> // for FLT_MAX


bool const MAP_VIEW_LIGHTING = 1;
bool const MAP_VIEW_SHADOWS  = 1;

int map_drag_x(0), map_drag_y(0);
float map_zoom(0.0);
double map_x(0.0), map_y(0.0);

extern bool water_is_lava, begin_motion;
extern int window_width, window_height, xoff2, yoff2, map_mode, map_color, read_landscape, read_heightmap, do_read_mesh;
extern int world_mode, display_mode, num_smileys, DISABLE_WATER, cache_counter, default_ground_tex, frame_counter;
extern unsigned show_map_view_fractal;
extern float zmax_est, zmin, zmax, water_plane_z, water_h_off, glaciate_exp, glaciate_exp_inv, vegetation, relh_adj_tex, temperature, mesh_height_scale, mesh_scale;
extern int coll_id[];
extern point surface_pos, camera_last_pos;
extern obj_group obj_groups[];
extern coll_obj_group coll_objects;


bool setup_height_gen(mesh_xy_grid_cache_t &height_gen, float x0, float y0, float dx, float dy, unsigned nx, unsigned ny, bool cache_values, bool no_wait=0);
bool using_hmap_with_detail();
void set_temp_clear_color(colorRGBA const &clear_color, bool clear_depth=0, bool clear_stencil=0);
float get_heightmap_scale();


struct complex_num {
	double r, i;
	complex_num() : r(0), i(0) {}
	complex_num(double r_, double i_) : r(r_), i(i_) {}
	complex_num operator+ (complex_num const &n) const {return complex_num((r+n.r), (i+n.i));}
	complex_num operator* (complex_num const &n) const {return complex_num((r*n.r - i*n.i), (r*n.i + i*n.r));}
	complex_num conjugaqte(complex_num const &n) const {return complex_num((r*n.r - i*n.i), -2.0*r*n.i);}
	complex_num abs() const {return complex_num(fabs(r), fabs(i));}
	double mag_sq()   const {return (r*r + i*i);}
	double mag   ()   const {return sqrt(mag_sq());}
};

double eval_fractal_set(complex_num const &c) {

	complex_num z(0.0, 0.0);
	unsigned val(0);
	
	switch (show_map_view_fractal) {
	case 0: // disabled, error?
		break;
	case 1: // madelbrot
		for (; val < 200; ++val) {
			if (z.mag_sq() > 4.0) break;
			z = z*z + c;
		}
		break;
	case 2: // tricorn
		for (; val < 200; ++val) {
			if (z.mag_sq() > 4.0) break;
			z = z.conjugaqte(z) + c;
		}
		break;
	case 3: // burning ship
		for (; val < 200; ++val) {
			if (z.mag_sq() > 4.0) break;
			z = z.abs();
			z = z*z + c;
		}
		break;
	} // end switch
	return (double(val) - log2(log2(z.mag_sq())) + 1.0)/200.0; // from http://www.iquilezles.org/www/articles/mset_smooth/mset_smooth.htm
	//return val/200.0;
}


float get_mesh_height(mesh_xy_grid_cache_t const &height_gen, float xstart, float ystart, float xscale, float yscale, int i, int j, bool nearest_texel=0) {

	if (using_tiled_terrain_hmap_tex()) {
		float zval(get_tiled_terrain_height_tex((xstart + X_SCENE_SIZE + j*xscale)*DX_VAL_INV, (ystart + Y_SCENE_SIZE + i*yscale)*DY_VAL_INV, nearest_texel));
		if (using_hmap_with_detail()) {zval += HMAP_DETAIL_MAG*height_gen.eval_index(j, i);}
		return zval;
	}
	return height_gen.eval_index(j, i);
}

bool is_shadowed(point const &cpos, vector3d const &cnorm, point const &lpos, int &cindex) {

	if (!MAP_VIEW_SHADOWS)   return 0;
	if (display_mode & 0x20) return 0;
	point const cpos2(cpos + 0.001*cnorm);
	if (cindex >= 0 && coll_objects.get_cobj(cindex).line_intersect(cpos2, lpos)) return 1;
	return check_coll_line(cpos2, lpos, cindex, -1, 1, 3); // static cobj shadows only for performance
}


void colorize(float val, unsigned char *rgb) {

	//rgb[0] = rgb[1] = rgb[2] = (unsigned char)(255.0*val); return;
	float const a(5*val), b(7*val), c(11*val);
	rgb[0] = (unsigned char)(255.0f*(a - int(a)));
	rgb[1] = (unsigned char)(255.0f*(b - int(b)));
	rgb[2] = (unsigned char)(255.0f*(c - int(c)));
}


void draw_overhead_map() {

	unsigned tid(0);
	if (map_mode == 0) return;
	if (map_mode == 2) {map_mode = 1; return;}
	if (map_zoom == 0.0) {map_zoom = ((world_mode == WMODE_GROUND) ? 0.08 : 0.8);} // set reasonable defaults based on mode
	int bx1(0), by1(0), bx2(0), by2(0), nx(1), ny(1);
	while (window_width  > 2*nx) {nx *= 2;}
	while (window_height > 2*ny) {ny *= 2;}
	//nx = (window_width & 0xFFFC); ny = (window_height & 0xFFFC); // looks nicer, but slower
	if (nx < 4 || ny < 4) return;

	//timer_t timer("Map Draw");
	int const nx2(nx/2), ny2(ny/2);
	bool const no_water((DISABLE_WATER == 2) || !(display_mode & 0x04));
	bool const is_ice(((world_mode == WMODE_GROUND) ? temperature : get_cur_temperature()) <= W_FREEZE_POINT);
	float const zmax2(zmax_est*((map_color || no_water) ? 1.0 : 0.855)), hscale(0.5/zmax2);
	float const window_ar((float(window_width)*ny)/(float(window_height)*nx)), scene_ar(X_SCENE_SIZE/Y_SCENE_SIZE);
	float const xscale(2.0*map_zoom*window_ar*HALF_DXY), yscale(2.0*map_zoom*scene_ar*HALF_DXY);
	float const xscale_val(xscale/64), yscale_val(yscale/64);

	// translate map_drag_x/y (screen pixel space) into map_x/y (world unit space)
	double const x_scale(nx*xscale/window_width), y_scale(ny*yscale/window_height);
	map_x += x_scale*map_drag_x; map_drag_x = 0;
	map_y += y_scale*map_drag_y; map_drag_y = 0;
	unsigned const tot_sz(nx*ny);
	vector<unsigned char> buf(tot_sz*3*sizeof(unsigned char));

	if (show_map_view_fractal) {
		double const y_scale(10.0*map_zoom), x_scale(window_ar*y_scale);
		double const i_scale(2.0*y_scale/ny), j_scale(2.0*x_scale/nx);
		double const x_off(-x_scale + 0.05*map_x), y_off(-y_scale + 0.05*map_y);
		//timer_t timer("Mandelbrot");

#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < ny; ++i) {
			double const my(i_scale*i + y_off);

			for (int j = 0; j < nx; ++j) {
				double const mx(j_scale*j + x_off);
				double const val(eval_fractal_set(complex_num(mx, -my)));
				colorize(val, &buf[3*(i*nx + j)]);
			}
		}
	}
	else {
		float x0((float)map_x + xoff2*DX_VAL), y0((float)map_y + yoff2*DY_VAL);
		float const relh_water(get_rel_height_no_clamp(water_plane_z, -zmax_est, zmax_est));
		point const camera(get_camera_pos());
		float map_heights[6] = {};
		map_heights[0] = 0.9f*lttex_dirt[3].zval  + 0.1f*lttex_dirt[4].zval;
		map_heights[1] = 0.5f*(lttex_dirt[2].zval + lttex_dirt[3].zval);
		map_heights[2] = 0.5f*(lttex_dirt[1].zval + lttex_dirt[2].zval);
		map_heights[3] = 0.5f*(lttex_dirt[0].zval + lttex_dirt[1].zval);
		map_heights[4] = relh_water; // Note: can be negative
		map_heights[5] = min(0.5f*relh_water, relh_water-0.01f); // handle negative case
		
		for (unsigned i = 0; i < 6; ++i) {
			if (map_heights[i] > 0.0) {map_heights[i] = pow(map_heights[i], glaciate_exp);} // handle negative case
		}
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
		float const xsv(xscale_val*(X_SCENE_SIZE/DX_VAL)), ysv(yscale_val*(Y_SCENE_SIZE/DY_VAL));
		float const max_building_dz(2.0*get_buildings_max_extent().z); // pad by 2x

		bool const uses_hmap(world_mode == WMODE_GROUND && (read_landscape || read_heightmap || do_read_mesh));
		mesh_xy_grid_cache_t height_gen;
		if (!uses_hmap && !show_map_view_fractal) {setup_height_gen(height_gen, xstart, ystart, xscale, yscale, nx, ny, 1);} // cache_values=1
		point const lpos(get_light_pos());
		vector3d const light_dir(lpos.get_norm()); // assume directional lighting to origin
		float const texels_per_pixel(mesh_scale*0.5f*(xscale*DX_VAL_INV + yscale*DY_VAL_INV));
		bool const nearest_texel(texels_per_pixel >= 1.0);

#pragma omp parallel for schedule(static,1)
		for (int i = 0; i < ny; ++i) {
			int const inx(i*nx);
			int64_t const iyy(((int64_t)i - (int64_t)yy)*((int64_t)i - (int64_t)yy)), icy(((int64_t)i - (int64_t)cy)*((int64_t)i - (int64_t)cy));
			float last_height(0.0);
			point cpos;
			vector3d cnorm;
			int cindex(-1), cindex2(-1);

			for (int j = 0; j < nx; ++j) {
				int const offset(3*(inx + j));
				unsigned char *rgb(&buf[offset]);
				int64_t const jxx((int64_t)j - (int64_t)xx), jcx((int64_t)j - (int64_t)cx);

				if (iyy + jxx*jxx <= 4) {
					rgb[0] = rgb[1] = rgb[2] = 0; // camera direction
				}
				else if (icy + jcx*jcx <= 9) {
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
					float const xval((j - nx2)*xsv + camera.x + map_x), yval((i - ny2)*ysv + camera.y + map_y);

					if (world_mode == WMODE_GROUND) {
						point p1(xval, yval, czmax);
						bool const over_mesh(is_over_mesh(p1));
						colorRGBA building_color;
					
						if (over_mesh || uses_hmap) { // if using a heightmap, clamp values to scene bounds
							mh = interpolate_mesh_zval(max(-X_SCENE_SIZE, min(X_SCENE_SIZE-DX_VAL, xval)), max(-Y_SCENE_SIZE, min(Y_SCENE_SIZE-DY_VAL, yval)), 0.0, 0, 1);
							mh_set = 1;
						}
						if (over_mesh && get_buildings_line_hit_color(point(xval, yval, mh+max_building_dz), point(xval, yval, mh), building_color)) {
							// if the player is inside a building, should we draw room labels with text? or is that too difficult and not worth the effort?
							unpack_color(rgb, building_color); // no shadows
							continue;
						}
						if (over_mesh && czmin < czmax) { // check cobjs
							// Note: as an optimization, can skip the cobj test if no cobjs at this pos, but it makes little difference and will miss dynamic objects
							//int const xpos(get_xpos(xval)), ypos(get_ypos(yval));
							//if (point_outside_mesh(xpos, ypos) || v_collision_matrix[ypos][xpos].zmin == v_collision_matrix[ypos][xpos].zmax) {}
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
					} // end ground mode
					else if (world_mode == WMODE_INF_TERRAIN && (have_cities() || have_buildings())) { // show cities and road networks
						colorRGBA city_color(BLACK);

						if (get_buildings_line_hit_color(point(xval, yval, zmax+max_building_dz), point(xval, yval, zmin), city_color)) {
							unpack_color(rgb, city_color); // no shadows
							continue;
						}
						if (get_city_color_at_xy(xval, yval, city_color)) {
							unpack_color(rgb, city_color); // no shadows
							continue;
						}
					}
					if (default_ground_tex >= 0 && map_color) {
						unpack_color(rgb, ground_color*(shadowed ? 0.5 : 1.0));
						continue;
					}
					if (!mh_set) {mh = get_mesh_height(height_gen, xstart, ystart, xscale, yscale, i, j, nearest_texel);} // calculate mesh height here if not yet set
					float height(min(1.0f, hscale*(mh + zmax2))); // can be negative

					if (!map_color) { // grayscale
						float const val(pow(height, glaciate_exp_inv)); // un-glaciate: slow
						//rgb[0] = rgb[1] = rgb[2] = (unsigned char)(255.0*val);
						// http://c0de517e.blogspot.com/2017/11/coder-color-palettes-for-data.html
						rgb[0] = (unsigned char)(255.0*(-0.121 + 0.893 * val + 0.276 * sin (1.94 - 5.69 * val)));
						rgb[1] = (unsigned char)(255.0*(0.07 + 0.947 * val));
						rgb[2] = (unsigned char)(255.0*(0.107 + (1.5 - 1.22 * val) * val));
					}
					else {
						height += relh_adj_tex;
						colorRGBA color;
						if      (height <= map_heights[5]) {color = map_colors[5];} // deep water
						else if (height <= map_heights[3]) {color = map_colors[3];} // sand
						else if (height >= map_heights[0]) {color = map_colors[0];} // snow
						else {
							color = BLACK;
							for (unsigned k = 0; k < 4; ++k) { // mixed
								if (height > map_heights[k+1]) {
									float const h((height - map_heights[k+1])/(map_heights[k] - map_heights[k+1])), v(cubic_interpolate(h));
									blend_color(color, map_colors[k], map_colors[k+1], v);
									break;
								}
							}
						}
						if (height <= map_heights[4] && height > map_heights[5]) { // shallow water
							float const h(0.5f*(height - map_heights[5])/(map_heights[4] - map_heights[5])), v(cubic_interpolate(h));
							blend_color(color, color, map_colors[5], v);
						}
						if (MAP_VIEW_LIGHTING && !uses_hmap && !(display_mode & 0x20)) {
							vector3d normal(plus_z);

							if (height > map_heights[4]) {
								float const mh2(get_mesh_height(height_gen, xstart, ystart, xscale, yscale, max(i-1, 0), j, nearest_texel));
								float const hx((j == 0) ? height : last_height), hy(CLIP_TO_01(hscale*(mh2 + zmax2)) + relh_adj_tex);
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
		if (begin_motion && obj_groups[coll_id[SMILEY]].enabled) {
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
	}
	set_temp_clear_color(BLACK, 1); // clear_depth=1
	shader_t s;
	s.begin_simple_textured_shader(0.0, 0, 0, &WHITE);
	setup_texture(tid, 0, 0, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, &buf.front());
	ensure_filled_polygons();
	draw_tquad(0.58*((float)window_width)/((float)window_height), 0.58, -1.0);
	reset_fill_mode();
	free_texture(tid);
	s.end_shader();
}


void place_player_at_xy(float xval, float yval) {
	surface_pos.assign(xval, yval, interpolate_mesh_zval(xval, yval, CAMERA_RADIUS, 0, 0));
}
void teleport_to_map_location() {
	static int last_update_frame(0);
	if ((frame_counter - last_update_frame) < 1.0f*TICKS_PER_SECOND) return; // teleport at most once per second if player holds down the key
	last_update_frame = frame_counter;
	float const xval(surface_pos.x + map_x), yval(surface_pos.y + map_y);
	place_player_at_xy(xval, yval);
	camera_last_pos = surface_pos; // avoid slow falling and rising on map teleport
	map_x = map_y = 0.0; // recenter on the new location
}


void get_heightmap_z_range(vector<float> const &heights, float &min_z, float &max_z) {
	min_z =  FLT_MAX;
	max_z = -FLT_MAX;

	for (unsigned i = 0; i < heights.size(); ++i) {
		min_eq(min_z, heights[i]);
		max_eq(max_z, heights[i]);
	}
}

void write_map_mode_heightmap_image() {

	float const window_ar((float(window_width)*window_height)/(float(window_height)*window_width));
	float const xscale(2.0*map_zoom*window_ar*HALF_DXY), yscale(2.0*map_zoom*(X_SCENE_SIZE/Y_SCENE_SIZE)*HALF_DXY);
	float const xstart((float)map_x + xoff2*DX_VAL - float(window_width/2)*xscale), ystart((float)map_y + yoff2*DY_VAL - float(window_height/2)*yscale);
	int const x1(get_xpos(xstart)), y1(get_ypos(ystart)), x2(get_xpos(xstart + window_width*xscale)), y2(get_ypos(ystart + window_height*yscale)), width(x2 - x1), height(y2 - y1);
	cout << "Heightmap image size: " << width << "x" << height << " = " << width*height/1024 << "K" << endl;
	if (width > 16384 || height > 16384) {std::cerr << "Error: heightmap image is too large, max size is 16384 pixels" << endl; return;} // fail

	string const fn("heightmap.png");
	texture_t texture(0, 6, width, height, 0, 2, 0, fn); // two bytes per pixel grayscale
	texture.set_16_bit_grayscale();
	texture.alloc();
	{ // open a scope
		timer_t timer("Heightmap Gen");
		vector<float> heights(texture.num_pixels());
		mesh_xy_grid_cache_t height_gen;
		setup_height_gen(height_gen, xstart, ystart, DX_VAL, DY_VAL, width, height, 1); // cache_values=1

	#pragma omp parallel for schedule(static,1)
		for (int i = 0; i < height; ++i) {
			int const off(width*(height - i - 1)); // invert yval
			for (int j = 0; j < width; ++j) {heights[off + j] = get_mesh_height(height_gen, xstart, ystart, DX_VAL, DY_VAL, i, j);}
		}
		float min_z(0), max_z(0);
		get_heightmap_z_range(heights, min_z, max_z);
		float const dz(max(TOLERANCE, (max_z - min_z))), height_scale(255.0/dz); // prevent divide-by-zero
		cout << "zval range: " << min_z << " to " << max_z << " total: " << dz << " scale: " << dz/(get_heightmap_scale()*mesh_height_scale) << endl;
		for (unsigned i = 0; i < heights.size(); ++i) {texture.write_pixel_16_bits(i, (heights[i] - min_z)*height_scale);}
	}
	cout << "Writing heightmap to image file " << fn << endl;
	timer_t timer("Heightmap Image Write");
	texture.write_to_png(fn);
}
