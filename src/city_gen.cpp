// 3D World - City Generation
// by Frank Gennari
// 2/10/18

#include "city.h"
#include "city_objects.h"
#include "city_terrain.h"
#include "heightmap.h"
#include "lightmap.h"
#include "buildings.h"
#include "profiler.h"
#include <cfloat> // for FLT_MAX

bool const CHECK_HEIGHT_BORDER_ONLY = 1; // choose building site to minimize edge discontinuity rather than amount of land that needs to be modified
float const CAR_LANE_OFFSET         = 0.15; // in units of road width
float const CITY_LIGHT_FALLOFF      = 0.2;


bool had_building_interior_coll(0), city_lights_custom_bcube(0);
vector2d actual_max_road_seg_len;
city_params_t city_params;
point pre_smap_player_pos(all_zeros), actual_player_pos(all_zeros); // Note: pre_smap_player_pos can be security cameras, but actual_player_pos is always the player

extern bool enable_dlight_shadows, dl_smap_enabled, enable_dlight_bcubes, flashlight_on, camera_in_building, have_indir_smoke_tex, disable_city_shadow_maps;
extern bool player_in_walkway, player_in_skyway, player_on_moving_ww, player_in_ww_elevator, player_in_tunnel;
extern int rand_gen_index, display_mode, animate2, draw_model, player_in_basement;
extern unsigned shadow_map_sz, cur_display_iter;
extern float cobj_z_bias, rain_wetness, NEAR_CLIP;
extern building_params_t global_building_params;
extern vector<light_source> dl_sources;


void add_dynamic_lights_city(cube_t const &scene_bcube, float &dlight_add_thresh, float falloff);
void add_buildings_exterior_lights(vector3d const &xlate, cube_t &lights_bcube);
void disable_shadow_maps(shader_t &s);
vector3d get_tt_xlate_val();
float get_max_house_size();
bool proc_buildings_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm=nullptr, bool check_interior=0, bool exclude_city=0);
void draw_player_building_transparent(int reflection_pass, vector3d const &xlate);


template<typename S, typename T> void get_all_bcubes(vector<T> const &v, S &bcubes) {
	for (auto i = v.begin(); i != v.end(); ++i) {bcubes.push_back(*i);}
}

bool is_night(float adj) {return (light_factor - adj < 0.5f);} // for car headlights and streetlights

void set_city_lighting_shader_opts(shader_t &s, cube_t const &lights_bcube, bool use_dlights, bool use_smap, float pcf_scale) {

	if (use_dlights) {
		s.setup_scene_bounds_from_bcube(lights_bcube); // reset with correct values
		s.add_uniform_float("LT_DIR_FALLOFF", CITY_LIGHT_FALLOFF); // smooth falloff for car headlights and streetlights

		if (flashlight_on) { // use a quicker falloff for flashlight beams
			s.add_uniform_float("LT_DIR_FALLOFF_SM", 0.02);
			s.add_uniform_float("LDIR_FALL_THRESH",  1.5*FLASHLIGHT_BW);
		}
	}
	if (use_smap) {
		s.add_uniform_float("z_bias", pcf_scale*cobj_z_bias); // I guess pcf_scale is really a light size scale and should apply to the z-bias as well
		s.add_uniform_float("shad_bias_scale", CITY_BIAS_SCALE); // fix for sun/moon shadows that are too far shifted on buildings
	}
}

// use_smap: 0=no, 1=sun/moon + dynamic lights; enable in shader and set shadow map uniforms, 2=dynamic lights only; disable in shader but set shadow map uniforms
void city_shader_setup(shader_t &s, cube_t const &lights_bcube, bool use_dlights, int use_smap, int use_bmap,
	float min_alpha, bool force_tsl, float pcf_scale, int use_texgen, bool indir_lighting, bool is_outside)
{
	use_dlights &= (lights_bcube.is_strictly_normalized() && !dl_sources.empty());
	have_indir_smoke_tex = indir_lighting; // assume someone is going to set the indir texture in this case; ***note that this breaks normal indir scene drawing***
	if (indir_lighting) {s.set_prefix("#define USE_ALT_SCENE_BOUNDS", 1);} // FS; need to use different scene_llc_scale for dynamic lighting vs. building indir lighting
	// Note: here use_texgen mode 5 is used as a hack so that the shader still has binding points for tex coords (can't optimize it out)
	// and we can share the same VAO between texgen and texcoords modes without having to worry about which mode we were in when the VAO was created;
	// use texgen mode 6 instead for cylinder buildings
	int const use_texgen_val(use_texgen ? ((use_texgen == 2) ? 6 : 5) : 0);
	bool const keep_alpha = 1; // required for fog on windows
	setup_smoke_shaders(s, min_alpha, use_texgen_val, keep_alpha, indir_lighting, 1, use_dlights, 0, 0,
		((use_smap == 1) ? 2 : 0), use_bmap, 0, (use_dlights || indir_lighting), force_tsl, 0.0, 0.0, 0, 0, is_outside); // use_spec_map=0
	set_city_lighting_shader_opts(s, lights_bcube, use_dlights, (use_smap != 0), pcf_scale);
	if (use_texgen) {s.add_uniform_float("tc_texgen_mix", 0.0);} // always uses texgen in this mode
	//if (use_smap) {bind_default_sun_moon_smap_textures();} // bind default sun/moon smap textures
}

void draw_state_t::begin_tile(point const &pos, bool will_emit_now, bool ensure_active) {
	if (ensure_active) {ensure_shader_active();} // needed for use_smap=0 case
	emit_now = (use_smap && try_bind_tile_smap_at_point((pos + xlate), s));
	if (will_emit_now && !emit_now) {disable_shadow_maps(s);} // not using shadow maps or second (non-shadow map) pass - disable shadow maps
}
void draw_state_t::pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_, bool always_setup_shader, bool enable_animations, bool enable_occlusion) {
	xlate       = xlate_;
	camera_bs   = camera_pdu.pos - xlate;
	shadow_only = shadow_only_;
	use_dlights = (use_dlights_ && !shadow_only);
	use_smap    = (shadow_map_enabled() && !shadow_only && !disable_city_shadow_maps);
	draw_tile_dist = get_draw_tile_dist();

	if (enable_occlusion && !shadow_only) {
		occlusion_checker.set_exclude_camera_building(); // if the player is inside a building, skip occlusion culling
		occlusion_checker.set_camera(camera_pdu);
	}
	if (!use_smap && !always_setup_shader) return;

	if (shadow_only) {
		if (enable_animations) {setup_smoke_shaders(s, 0.0, 0, 0, 0, 0, 0, 0);} // use main shader with all lighting and effects disabled to get animations
		else {s.begin_simple_textured_shader();}
	}
	else {
		bool const force_tsl = 0; // helps with hedges and flags, but causes problems with other models
		cube_t const &lights_bcube(use_building_lights ? get_building_lights_bcube() : get_city_lights_bcube());
		city_shader_setup(s, lights_bcube, use_dlights, use_smap, (use_bmap && !shadow_only), DEF_CITY_MIN_ALPHA, force_tsl, 0.5);
	}
}
void draw_state_t::end_draw() {
	emit_now = 0;
	if (s.is_setup()) {s.end_shader();} // use_smap case
}
/*virtual*/ void draw_state_t::post_draw() {
	end_draw();
	ensure_shader_active();
	draw_unshadowed();
	s.end_shader();
}
void draw_state_t::set_untextured_material() {
	if (shadow_only) return;
	select_texture(WHITE_TEX);
	if (normal_maps_enabled()) {s.add_uniform_float("bump_map_mag", 0.0);} // disable bump map
}
void draw_state_t::unset_untextured_material() {
	if (!shadow_only && normal_maps_enabled()) {s.add_uniform_float("bump_map_mag", 1.0);} // re-enable bump map
}
void draw_state_t::ensure_shader_active() {
	if (s.is_setup()) return; // already active
	if (shadow_only) {s.begin_shadow_map_shader();}
	else {city_shader_setup(s, get_city_lights_bcube(), use_dlights, 0, use_bmap, DEF_CITY_MIN_ALPHA, 0, 0.5);} // no smap
}
void draw_state_t::draw_and_clear_light_flares() {
	if (light_psd.empty()) return; // no lights to draw
	enable_blend();
	set_additive_blend_mode();
	glDepthMask(GL_FALSE); // disable depth write
	light_psd.draw_and_clear(BLUR_TEX, 0.0, 0, 1, 0.005); // use geometry shader for unlimited point size
	glDepthMask(GL_TRUE);
	set_std_blend_mode();
	disable_blend();
}
bool draw_state_t::check_cube_visible(cube_t const &bc, float dist_scale) const {
	if (!camera_pdu.valid) return 1;
	cube_t const bcx(bc + xlate);

	if (dist_scale > 0.0) {
		float const dmax(shadow_only ? camera_pdu.far_ : dist_scale*draw_tile_dist);
		if (!bcx.closest_dist_less_than(camera_pdu.pos, dmax)) return 0;
	}
	return camera_pdu.cube_visible(bcx);
}
/*static*/ void draw_state_t::set_cube_pts(cube_t const &c, float z1f, float z1b, float z2f, float z2b, bool d, bool D, point p[8]) {
	p[0][!d] = p[4][!d] = c.d[!d][1]; p[0][d] = p[4][d] = c.d[d][ D]; p[0].z = z1f; p[4].z = z2f; // front right
	p[1][!d] = p[5][!d] = c.d[!d][0]; p[1][d] = p[5][d] = c.d[d][ D]; p[1].z = z1f; p[5].z = z2f; // front left
	p[2][!d] = p[6][!d] = c.d[!d][0]; p[2][d] = p[6][d] = c.d[d][!D]; p[2].z = z1b; p[6].z = z2b; // back  left
	p[3][!d] = p[7][!d] = c.d[!d][1]; p[3][d] = p[7][d] = c.d[d][!D]; p[3].z = z1b; p[7].z = z2b; // back  right
}
/*static*/ void draw_state_t::rotate_pts(point const &center, float sine_val, float cos_val, int d, int e, point p[8]) {
	for (unsigned i = 0; i < 8; ++i) {
		point &v(p[i]); // rotate p[i]
		v -= center; // translate to origin
		float const a(v[d]*cos_val - v[e]*sine_val), b(v[e]*cos_val + v[d]*sine_val);
		v[d] = a; v[e] = b; // rotate
		v += center; // translate back
	}
}
void draw_state_t::draw_cube(quad_batch_draw &qbd, color_wrapper const &cw, point const &center, point const p[8],
	bool skip_bottom, bool invert_normals, float tscale, unsigned skip_dims) const
{
	vector3d const cview_dir(camera_bs - center);
	float const sign(invert_normals ? -1.0 : 1.0);
	vector3d const top_n  (cross_product((p[2] - p[1]), (p[0] - p[1]))*sign); // Note: normalization not needed
	vector3d const front_n(cross_product((p[5] - p[1]), (p[0] - p[1]))*sign);
	vector3d const right_n(cross_product((p[6] - p[2]), (p[1] - p[2]))*sign);
	tex_range_t tr_top, tr_front, tr_right;

	if (tscale > 0.0) { // compute texture s/t parameters from cube side lengths to get a 1:1 AR
		float const ts01(tscale*p2p_dist(p[0], p[1])), ts12(tscale*p2p_dist(p[1], p[2])), ts15(tscale*p2p_dist(p[1], p[5]));
		tr_top  .x2 = ts01; tr_top  .y2 = ts12;
		tr_front.x2 = ts01; tr_front.y2 = ts15;
		tr_right.x2 = ts12; tr_right.y2 = ts15;
	}
	if (!(skip_dims & 4)) { // Z
		if (dot_product(cview_dir, top_n) > 0) {qbd.add_quad_pts(p+4, cw,  top_n, tr_top);} // top
		else if (!skip_bottom)                 {qbd.add_quad_pts(p+0, cw, -top_n, tr_top);} // bottom - not always drawn
	}
	if (!(skip_dims & 1)) { // X
		if (dot_product(cview_dir, front_n) > 0) {point const pts[4] = {p[0], p[1], p[5], p[4]}; qbd.add_quad_pts(pts, cw,  front_n, tr_front);} // back
		else                                     {point const pts[4] = {p[2], p[3], p[7], p[6]}; qbd.add_quad_pts(pts, cw, -front_n, tr_front);} // front
	}
	if (!(skip_dims & 2)) { // Y
		if (dot_product(cview_dir, right_n) > 0) {point const pts[4] = {p[1], p[2], p[6], p[5]}; qbd.add_quad_pts(pts, cw,  right_n, tr_right);} // left
		else                                     {point const pts[4] = {p[3], p[0], p[4], p[7]}; qbd.add_quad_pts(pts, cw, -right_n, tr_right);} // right
	}
}
void draw_state_t::draw_cube(quad_batch_draw &qbd, cube_t const &c, color_wrapper const &cw, bool skip_bottom, float tscale, unsigned skip_dims,
	bool mirror_x, bool mirror_y, bool swap_tc_xy, float tscale_x, float tscale_y, float tscale_z, bool skip_top, bool no_cull) const
{
	point p[8];
	set_cube_pts(c, 0, 0, p);
	//draw_cube(qbd, cw, c.get_cube_center(), p, skip_bottom, 0, tscale, skip_dims); // customized for axis aligned cube below
	vector3d const cview_dir(no_cull ? zero_vector : (camera_bs - c.get_cube_center()));
	tex_range_t tr_top, tr_front, tr_right;
	if (swap_tc_xy) {tr_top.swap_xy = tr_front.swap_xy = tr_right.swap_xy = 1;}

	if (tscale != 0.0) { // tscale: compute texture s/t parameters from cube side lengths to get a 1:1 AR
		float const ts01(tscale*tscale_y*c.dy()), ts12(tscale*tscale_x*c.dx()), ts15(tscale*tscale_z*c.dz());
		tr_top  .x2 = ts01; tr_top  .y2 = ts12;
		tr_front.x2 = ts01; tr_front.y2 = ts15;
		tr_right.x2 = ts12; tr_right.y2 = ts15;
		if (swap_tc_xy) {swap(tr_top.x2, tr_top.y2); swap(tr_front.x2, tr_front.y2); swap(tr_right.x2, tr_right.y2);}
	}
	if (!(skip_dims & 4)) { // Z
		if (mirror_x) {tr_top.mirror_x();}
		if (mirror_y) {tr_top.mirror_y();}

		if (cview_dir.z >= 0.0) {
			if (!skip_top   ) {qbd.add_quad_pts(p+4, cw,  plus_z, tr_top);} // top - not always drawn
		}
		if (cview_dir.z <= 0.0) {
			if (!skip_bottom) {qbd.add_quad_pts(p+0, cw, -plus_z, tr_top);} // bottom - not always drawn
		}
	}
	if (!(skip_dims & 1)) { // X
		if (mirror_x) {tr_front.mirror_x();}
		if (mirror_y) {tr_front.mirror_y();}
		if (cview_dir.x <= 0.0) {point const pts[4] = {p[0], p[1], p[5], p[4]}; qbd.add_quad_pts(pts, cw, -plus_x, tr_front);} // back
		if (cview_dir.x >= 0.0) {point const pts[4] = {p[2], p[3], p[7], p[6]}; qbd.add_quad_pts(pts, cw,  plus_x, tr_front);} // front
	}
	if (!(skip_dims & 2)) { // Y
		if (mirror_x) {tr_right.mirror_x();}
		if (mirror_y) {tr_right.mirror_y();}
		if (cview_dir.y <= 0.0) {point const pts[4] = {p[1], p[2], p[6], p[5]}; qbd.add_quad_pts(pts, cw, -plus_y, tr_right);} // left
		if (cview_dir.y >= 0.0) {point const pts[4] = {p[3], p[0], p[4], p[7]}; qbd.add_quad_pts(pts, cw,  plus_y, tr_right);} // right
	}
}
bool draw_state_t::add_light_flare(point const &flare_pos, vector3d const &n, colorRGBA const &color, float alpha, float radius) {
	point pos(xlate + flare_pos);
	vector3d const view_dir((camera_pdu.pos - pos).get_norm());
	float dp((n == zero_vector) ? 1.0 : dot_product(n, view_dir)); // n == 0 => non-directional
	if (dp < 0.05) return 0; // back facing, skip
	pos += 0.75*radius*view_dir; // move toward the camera, away from the stoplight, to prevent clipping
	light_psd.add_pt(sized_vert_t<vert_color>(vert_color(pos, colorRGBA(color, dp*alpha)), radius));
	return 1;
}
void draw_state_t::show_label_text() {
	if (label_str.empty()) return;
	text_drawer_t text_drawer;
	text_drawer.strs.push_back(text_string_t(label_str, label_pos, 20.0, CYAN));
	text_drawer.draw();
	label_str.clear(); // drawn, clear for next frame
}


cube_t const &get_bcube(cube_t const &bcube) {return bcube;}

cube_t get_bcube(bridge_t const &bridge) {
	cube_t bcube(bridge);
	float const shrink(2.0*(bridge.dim ? DY_VAL : DX_VAL));
	bcube.d[bridge.dim][0] += shrink; bcube.d[bridge.dim][1] -= shrink;
	return bcube;
}
template<typename T> bool check_bcubes_sphere_coll(vector<T> const &bcubes, point const &sc, float radius, bool xy_only) {
	for (auto i = bcubes.begin(); i != bcubes.end(); ++i) {
		if (check_bcube_sphere_coll(get_bcube(*i), sc, radius, xy_only)) return 1;
	}
	return 0;
}
template bool check_bcubes_sphere_coll(vector<cube_t> const &bcubes, point const &sc, float radius, bool xy_only); // explicit instantiation

template<typename T> void get_bcubes_region_coll_xy(vector<T> const &bcubes, vect_cube_t &out, cube_t const &region, vector3d const &xlate) {
	for (T const &c : bcubes) {
		if (c.intersects_xy(region)) {out.push_back(get_bcube(c) + xlate);} // Note: out is in camera space
	}
}
template<typename T> cube_t calc_cubes_bcube(vector<T> const &cubes) {
	if (cubes.empty()) return cube_t(all_zeros);
	cube_t bcube(cubes.front()); // first cube
	for (auto r = cubes.begin()+1; r != cubes.end(); ++r) {bcube.union_with_cube(*r);} // skip first cube
	return bcube;
}

point rand_xy_pt_in_cube(cube_t const &c, float radius, rand_gen_t &rgen) {
	return point(rgen.rand_uniform(c.x1()+radius, c.x2()-radius), rgen.rand_uniform(c.y1()+radius, c.y2()-radius), c.z1());
}


class city_plot_gen_t : public heightmap_query_t {
protected:
	int last_rgi;
	rand_gen_t rgen;
	vector<rect_t> used;
	vect_cube_t plots; // same size as used
	cube_t bcube;

	bool overlaps_used(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		rect_t const cur(x1, y1, x2, y2);
		for (vector<rect_t>::const_iterator i = used.begin(); i != used.end(); ++i) {if (i->has_overlap(cur)) return 1;} // simple linear iteration
		return 0;
	}
	cube_t add_plot(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float elevation) {
		cube_t const c(get_cube_for_bounds(x1, y1, x2, y2, elevation));
		if (plots.empty()) {bcube = c;} else {bcube.union_with_cube(c);}
		plots.push_back(c);
		used.emplace_back(x1, y1, x2, y2);
		return c;
	}
	float get_avg_height(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		assert(is_normalized_region(x1, y1, x2, y2));
		float sum(0.0), denom(0.0);

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				sum   += get_height(x, y);
				denom += 1.0;
			}
		}
		return sum/denom;
	}
	float get_rms_height_diff(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		float const avg(get_avg_height(x1, y1, x2, y2));
		float diff(0.0);

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				float const delta(get_height(x, y) - avg);
				diff += delta*delta; // square the difference
			}
		}
		return diff;
	}
public:
	city_plot_gen_t() : last_rgi(0), bcube(all_zeros) {}
	void invalidate_heightmap() {heightmap = nullptr;}

	void init(float *heightmap_, unsigned xsize_, unsigned ysize_) {
		heightmap = heightmap_; xsize = xsize_; ysize = ysize_;
		assert(heightmap != nullptr);
		assert(xsize > 0 && ysize > 0); // any size is okay
		if (rand_gen_index != last_rgi) {rgen.set_state(rand_gen_index, 12345); last_rgi = rand_gen_index;} // only when rand_gen_index changes
	}
	bool find_best_city_location(unsigned wmin, unsigned hmin, unsigned wmax, unsigned hmax, unsigned border, unsigned slope_width, unsigned num_samples,
		unsigned &cx1, unsigned &cy1, unsigned &cx2, unsigned &cy2)
	{
		assert(num_samples > 0);
		if ((wmax + 2*border) >= xsize || (hmax + 2*border) >= ysize) return 0; // city can't fit in the map
		unsigned const num_iters(100*num_samples); // upper bound
		unsigned const xend(xsize - wmax - 2*border + 1), yend(ysize - hmax - 2*border + 1); // max rect LLC, inclusive
		assert(xend > 0 && yend > 0);
		unsigned num_cands(0);
		float best_diff(0.0);

		for (unsigned n = 0; n < num_iters; ++n) { // find min RMS height change across N samples
			unsigned const x1(border + (rgen.rand()%xend)), y1(border + (rgen.rand()%yend));
			unsigned const x2(x1 + ((wmin == wmax) ? wmin : rgen.rand_int(wmin, wmax)));
			unsigned const y2(y1 + ((hmin == hmax) ? hmin : rgen.rand_int(hmin, hmax)));
			if (overlaps_used (x1-slope_width, y1-slope_width, x2+slope_width, y2+slope_width)) continue; // skip if plot expanded by slope_width overlaps an existing city
			if (any_underwater(x1, y1, x2, y2, CHECK_HEIGHT_BORDER_ONLY)) continue; // skip
			float const diff(get_rms_height_diff(x1, y1, x2, y2));
			if (num_cands == 0 || diff < best_diff) {cx1 = x1; cy1 = y1; cx2 = x2; cy2 = y2; best_diff = diff;}
			if (++num_cands == num_samples) break; // done
		} // for n
		if (num_cands == 0) return 0;
		//cout << "City cands: " << num_cands << ", diff: " << best_diff << ", loc: " << (cx1+cx2)/2 << "," << (cy1+cy2)/2 << endl;
		return 1; // success
	}
	float flatten_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float const *const height=nullptr) {
		float const elevation(height ? *height : get_avg_height(x1, y1, x2, y2));
		flatten_region_to(x1, y1, x2, y2, slope_width, elevation);
		return elevation;
	}
	bool check_plot_sphere_coll(point const &pos, float radius, bool xy_only=1) const { // Note: cities, not blocks; pos is in camera space
		if (plots.empty()) return 0;
		point const query_pos(pos - get_camera_coord_space_xlate()); // convert from camera space to global city space
		radius += 0.5*city_params.road_width; // add extra padding around city plots
		if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
		return check_bcubes_sphere_coll(plots, query_pos, radius, xy_only);
	}
	void get_plots_in_region(cube_t const &region, vect_cube_t &out) const { // Note: region is in camera space, and out is returned in camera space
		if (plots.empty()) return;
		vector3d const xlate(get_camera_coord_space_xlate());
		cube_t const query_region(region - xlate); // in global city space
		if (!bcube.intersects_xy(query_region)) return;
		get_bcubes_region_coll_xy(plots, out, query_region, xlate);
	}
}; // city_plot_gen_t


void range_pair_t::update(unsigned v) {
	if (s == 0 && e == 0) {s = v;} // first insert
	else {assert(s < e && v >= e);} // v must strictly increase
	e = v+1; // one past the end
}

void plot_xy_t::gen_adj_plots(vector<road_plot_t> const &plots) {
	assert(plots.size() == num());
	adj_plots.clear();
	adj_plots.resize(plots.size()); // start at -1

	for (unsigned i = 0; i < plots.size(); ++i) {
		unsigned x(plots[i].xpos), y(plots[i].ypos), ix(y*nx + x);
		if (x+1 < nx) {adj_plots[ix+1 ].set_adj(0, i);} // i is to the west  of this plot
		if (x   > 0 ) {adj_plots[ix-1 ].set_adj(1, i);} // i is to the east  of this plot
		if (y+1 < ny) {adj_plots[ix+nx].set_adj(2, i);} // i is to the north of this plot
		if (y   > 0 ) {adj_plots[ix-nx].set_adj(3, i);} // i is to the south of this plot
	}
}


class road_network_t : public streetlights_t { // AKA city center

	vector<road_t> roads; // full overlapping roads with constant slope, for collisions, etc.
	vector<road_seg_t> segs; // non-overlapping road segments, for drawing with textures
	vector<cube_t> conn_roads; // connector road bounding cubes (contain multiple adjacent connected roads with different slopes/zvals)
	vector<road_isec_t> isecs[3]; // for drawing with textures: {2-way, 3-way, 4-way}
	vector<road_plot_t> plots; // plots of land that can hold buildings (city blocks)
	vector<bridge_t> bridges; // bridges, part of global road network
	vector<tunnel_t> tunnels; // tunnels, part of global road network
	vector<road_t> tracks, track_segs; // railroad tracks (for global road network)
	vect_cube_t parks;
	vect_cube_t plot_cuts; // for underground skylights, elevators, etc.
	vect_ug_elev_info_t uges;
	//vector<road_isec_t> track_turns; // for railroad tracks
	city_obj_placer_t city_obj_placer;
	cube_t bcube;
	set<unsigned> connected_to; // vector?
	map<uint64_t, unsigned> tile_to_block_map;
	map<unsigned, road_isec_t const *> cix_to_isec; // maps city_ix to intersection
	vector<vect_cube_t> plot_colliders;
	plot_xy_t plot_xy;
	unsigned city_id, cluster_id, plot_id_offset;
	//string city_name; // future work
	float tot_road_len;
	mutable unsigned num_cars; // Note: not counting parked cars; mutable so that car_manager can update this
	bool is_residential;

	// use only for the global road network
	struct city_id_pair_t {
		unsigned id[2]; // lo, hi
		city_id_pair_t(unsigned c1, unsigned c2) {id[0] = c1; id[1] = c2;}
	};
	vector<city_id_pair_t> road_to_city; // indexed by road ID
	vector<vector<unsigned>> city_to_seg; // maps city_id to set of road segments connecting to that city

	struct tile_block_t { // collection of road parts for a given tile
		range_pair_t ranges[NUM_RD_TYPES]; // {plot, seg, isec2, isec3, isec4, park_lot, tracks, park, driveway, road_skirt, building}
		quad_batch_draw quads[NUM_RD_TYPES];
		cube_t bcube;
		tile_block_t(cube_t const &bcube_) : bcube(bcube_) {}
	};
	vector<tile_block_t> tile_blocks;

	template<typename T> void add_tile_blocks(vector<T> &v, map<uint64_t, unsigned> &tile_to_block_map, unsigned type_ix) {
		assert(type_ix < NUM_RD_TYPES);
		sort(v.begin(), v.end(), cmp_by_tile());

		for (unsigned i = 0; i < v.size(); ++i) {
			uint64_t const tile_id(get_tile_id_for_cube(v[i]));
			auto it(tile_to_block_map.find(tile_id));
			unsigned block_id(0);
			
			if (it == tile_to_block_map.end()) { // not found, add new block
				tile_to_block_map[tile_id] = block_id = tile_blocks.size();
				tile_blocks.push_back(tile_block_t(v[i]));
			}
			else {block_id = it->second;}
			assert(block_id < tile_blocks.size());
			tile_blocks[block_id].ranges[type_ix].update(i);
			tile_blocks[block_id].bcube.union_with_cube(v[i]);
		} // for i
	}
	road_seg_t const &get_seg(unsigned seg_ix) const {
		assert(seg_ix < segs.size());
		return segs[seg_ix];
	}
	road_isec_t const &get_isec(unsigned type_ix, unsigned isec_ix) const {
		assert(type_ix < 3);
		assert(isec_ix < isecs[type_ix].size());
		return isecs[type_ix][isec_ix];
	}
public:
	road_network_t() : bcube(all_zeros), city_id(CONN_CITY_IX), cluster_id(0), plot_id_offset(0), tot_road_len(0.0), num_cars(0), is_residential(0) {} // global road network ctor
		
	road_network_t(cube_t const &bcube_, unsigned city_id_, bool is_residential_) :
		bcube(bcube_), city_id(city_id_), cluster_id(0), plot_id_offset(0), tot_road_len(0.0), num_cars(0), is_residential(is_residential_)
	{
		bcube.z2() += ROAD_HEIGHT; // make it nonzero size
	}
	bool get_is_residential() const {return is_residential;}
	cube_t const &get_bcube() const {return bcube;}
	cube_t const &get_plot_bcube(unsigned plot_ix) const {assert(plot_ix < plots.size()); return plots[plot_ix];}
	vector<power_pole_t> const &get_power_poles() const {return city_obj_placer.get_power_poles();} // used for city connectivity
	void set_bcube(cube_t const &bcube_) {bcube = bcube_;}
	unsigned num_roads() const {return roads.size();}
	vector<road_t> const &get_roads() const {return roads;} // used for connecting roads between cities with 4-way intersections
	bool empty() const {return roads.empty();}
	plot_xy_t const &get_plot_xy() const {return plot_xy;}
	bool has_tunnels() const {return !tunnels.empty();} // global connector road only
	void set_cluster(unsigned id) {cluster_id = id;}
	void register_connected_city(unsigned id) {connected_to.insert(id);}
	set<unsigned> const &get_connected() const {return connected_to;}
	bool is_connected_to(unsigned id) const {return (connected_to.find(id) != connected_to.end());}
	float get_traffic_density() const {return ((tot_road_len == 0.0) ? 0.0 : num_cars/tot_road_len);} // cars per unit road
	void register_car() const {++num_cars;} // Note: must be const; num_cars is mutable

	cube_t get_bcube_inc_stoplights_and_streetlights() const {
		cube_t c(bcube); // deep copy
		c.z2() += max(stoplight_ns::stoplight_max_height(), streetlight_ns::get_streetlight_height());
		return c;
	}
	bool gen_road_grid(float road_width, vector2d const &road_spacing) {
		if (road_width > 0.5*min(road_spacing.x, road_spacing.y)) {
			cerr << "Error: City road_width should not be set larger than half the road spacing" << endl;
			exit(1);
		}
		vector3d const size(bcube.get_size()); // use our bcube as the region to process
		assert(size.x > 0.0 && size.y > 0.0);
		float const half_width(0.5*road_width), zval(bcube.z1() + ROAD_HEIGHT);
		float const rx1(bcube.x1() + half_width), rx2(bcube.x2() - half_width), ry1(bcube.y1() + half_width), ry2(bcube.y2() - half_width); // shrink to include centerlines
		float road_pitch_x(road_width + road_spacing.x), road_pitch_y(road_width + road_spacing.y);
		int const num_x_roads((rx2 - rx1)/road_pitch_x), num_y_roads((ry2 - ry1)/road_pitch_y);
		road_pitch_x = 0.9999f*(rx2 - rx1)/num_x_roads; // auto-calculate, round down slightly to avoid FP error
		road_pitch_y = 0.9999f*(ry2 - ry1)/num_y_roads;
		max_eq(actual_max_road_seg_len.x, (road_pitch_x - road_width));
		max_eq(actual_max_road_seg_len.y, (road_pitch_y - road_width));

		// create a grid, for now; crossing roads will overlap
		for (float x = rx1; x < rx2; x += road_pitch_x) {
			roads.emplace_back(point(x, bcube.y1(), zval), point(x, bcube.y2(), zval), road_width, true, false, roads.size());
		}
		unsigned const num_x(roads.size());

		for (float y = ry1; y < ry2; y += road_pitch_y) {
			roads.emplace_back(point(bcube.x1(), y, zval), point(bcube.x2(), y, zval), road_width, false, false, roads.size());
		}
		unsigned const num_r(roads.size()), num_y(num_r - num_x);
		if (num_x <= 1 || num_y <= 1) {roads.clear(); return 0;} // not enough space for roads; this city will be removed
		bcube.x1() = roads[0      ].x1(); // actual bcube x1 from first x road
		bcube.x2() = roads[num_x-1].x2(); // actual bcube x2 from last  x road
		bcube.y1() = roads[num_x  ].y1(); // actual bcube y1 from first y road
		bcube.y2() = roads[num_r-1].y2(); // actual bcube y2 from last  y road

		// create road segments and intersections
		bool const add_stoplights(!is_residential); // add stoplights in cities, and stop signs in residential neighborhoods
		plot_xy.nx = num_x - 1; plot_xy.ny = num_y - 1;
		segs .reserve(num_x*(num_y-1) + (num_x-1)*num_y + 4); // X + Y segments, allocate one extra per side for connectors
		plots.reserve(plot_xy.num());

		if (num_x > 2 && num_y > 2) {
			isecs[0].reserve(4); // 2-way, always exactly 4 at each corner
			isecs[1].reserve(2*((num_x-2) + (num_y-2)) + 4); // 3-way, allocate one extra per side for connectors
			isecs[2].reserve((num_x-2)*(num_y-2) + 4); // 4-way, allocate one extra per side for connectors
		}
		for (unsigned x = 0; x < num_x; ++x) {
			for (unsigned y = num_x; y < num_r; ++y) {
				bool const FX(x == 0), FY(y == num_x), LX(x+1 == num_x), LY(y+1 == num_r);
				cube_t const &rx(roads[x]), &ry(roads[y]);
				unsigned const num_conn((!FX) + (!LX) + (!FY) + (!LY));
				if (num_conn < 2) continue; // error?
				uint8_t const conn(((!FX) << 0) | ((!LX) << 1) | ((!FY) << 2) | ((!LY) << 3)); // 1-15
				isecs[num_conn - 2].emplace_back(cube_t(rx.x1(), rx.x2(), ry.y1(), ry.y2(), zval, zval), y, x, conn, false, (add_stoplights && num_conn > 2)); // intersections
					
				if (!LX) { // skip last y segment
					cube_t const &rxn(roads[x+1]);
					segs.emplace_back(cube_t(rx.x2(), rxn.x1(), ry.y1(), ry.y2(), zval, zval), y, false); // y-segments
				}
				if (!LY) { // skip last x segment
					cube_t const &ryn(roads[y+1]);
					segs.emplace_back(cube_t(rx.x1(), rx.x2(), ry.y2(), ryn.y1(), zval, zval), x, true); // x-segments

					if (!LX) { // skip last y segment
						cube_t const &rxn(roads[x+1]);
						plots.emplace_back(cube_t(rx.x2(), rxn.x1(), ry.y2(), ryn.y1(), zval, zval), uint8_t(x), uint8_t(y - num_x), is_residential); // plots between roads
					}
				}
			} // for y
		} // for x
		plot_colliders.resize(plots.size());

		if (city_params.park_rate > 0) { // make some plots into parks
			rand_gen_t rgen;
			rgen.set_state(plots.size(), city_id+1);
			
			for (auto p = plots.begin(); p != plots.end(); ++p) {
				if ((rgen.rand() % city_params.park_rate) == 0) {p->is_park = 1; parks.push_back(*p);}
			}
		}
		return 1;
	}

	// called on road (road_seg_t) or track (road_t) segments
	template<typename T> bool handle_crossings(T begin, T end, point &p1, point &p2, float width, bool dim, vector<road_t> &to_re_flatten) {
		cube_t seg_bcube(p1, p2);
		seg_bcube.expand_in_dim(!dim, 0.5*width);

		for (T i = begin; i != end; ++i) {
			if (!i->intersects_xy(seg_bcube)) continue;
			// shouldn't return here because any intersection with a road parallel to the track should also intersect a city;
			// however, this is safer for the cases where a jog is allowed in a road
			if (i->dim == dim) return 0; // failed, invalid track pos
			float const prev_center(i->get_center_dim(dim)), cur_center(p1[!dim]);
			float const tt((prev_center - p1[dim])/(p2[dim] - p1[dim]));
			float const tr((cur_center  - i->d[!dim][0])/i->get_length());
			if (tt <= 0.0 || tt >= 1.0 || tr <= 0.0 || tr >= 1.0) continue; // int, but centerlines not crossing - some other seg intersects, or already handled this seg
			point p_int(p1);
			p_int[dim] = prev_center;
			float const z_cur(p1.z + tt*(p2.z - p1.z)); // interpolated from endpoints, assuming track is a constant slope
			// Note: we can't just call hq.get_road_zval_at_pt(p_int) here becase the terrain may not have been flattened yet (for prev placed tracks)
			float const z_prev(i->get_start_z() + tr*(i->get_end_z() - i->get_start_z())); // interpolated from endpoints, assuming road seg is a constant slope
			float const delta_z(z_cur - z_prev); // positive: track is above the road; negative: track is below the road
			if      (0 && delta_z >  4.0*width) {} // bridge?
			else if (0 && delta_z < -4.0*width) {} // tunnel?
			else { // adjust the tracks to meet the road
				p2   = p_int; // this segment will end at the road's center
				p2.z = z_prev;
				to_re_flatten.push_back(*i);
			}
			return 1; // assume we can hit at most one road
		} // for i
		return 1;
	}
	void gen_railroad_tracks(vect_cube_t const &blockers, heightmap_query_t &hq) { // global connector road only, for now
		float const width(TRACKS_WIDTH*city_params.road_width);
		cube_t const region(calc_cubes_bcube(blockers));
		assert(region.dx() > 0.0 && region.dy() > 0.0);
		if (region.dx() <= 2.0*width || region.dy() <= 2.0*width) return; // region too small (shouldn't happen)
		vect_cube_t dim_tracks[2]; // one in each dim, for collision detection with tracks going in the other dim
		float const step_sz(city_params.conn_road_seg_len), max_seg_len(city_params.road_spacing);
		unsigned const num_tries(2*max(city_params.num_conn_tries, 1U)); // twice as many tries as connector roads: try full length, then shorter lengths
		tracks.reserve(city_params.num_rr_tracks); // to avoid iterator invalidation
		rand_gen_t rgen;
		vector<road_t> to_re_flatten;

		for (unsigned n = 0; n < city_params.num_rr_tracks; ++n) {
			for (unsigned tries = 0; tries < num_tries; ++tries) {
				bool const dim(rgen.rand_bool());
				float const rv(rgen.rand_uniform(0.2, 0.8)); // use center area
				float const pos(region.d[!dim][0]*(1.0f - rv) + region.d[!dim][1]*rv);
				float const seg_start(region.d[dim][0]), seg_end(region.d[dim][1]);
				point p1, p2;
				p1[!dim] = p2[!dim] = pos;
				p1[dim]  = seg_start; p2[dim] = seg_end; // full segment for blockers check
				p1.z     = p2.z = region.z1();
				road_t track(p1, p2, width, dim, (p2.z < p1.z), n); // Note: zvals are at 0, but should be unused
				if (has_bcube_int_xy(track, blockers,            width)) continue; // check cities
				if (has_bcube_int_xy(track, dim_tracks[dim], 8.0*width)) continue; // check prev placed tracks in same dim
				p2[dim] = (p1[dim] + step_sz); // back to starting segment
				p2.z    = hq.get_road_zval_at_pt(p1); // will be used as the initial p1 zval
				unsigned const segs_start(track_segs.size());
				bool valid(1);

				while (p1[dim] < seg_end) { // split into per-tile segments
					p1.z = p2.z;
					p2.z = hq.get_road_zval_at_pt(p2);
					// check for collisions with previously placed tracks and connector roads, and handle them with intersections, bridges, or tunnels
					auto track_segs_end(track_segs.begin() + segs_start);
					if (!handle_crossings(track_segs.begin(), track_segs_end, p1, p2, width, dim, to_re_flatten)) {valid = 0; break;} // check prev placed tracks
					if (!handle_crossings(segs.begin(),       segs.end(),     p1, p2, width, dim, to_re_flatten)) {valid = 0; break;} // check connector roads
					float const seg_len(p2[dim] - p1[dim]), dz(fabs(p2.z - p1.z));
	
					if (dz/seg_len > city_params.max_track_slope) { // check the max slope
						if (tries < num_tries/2) {valid = 0; break;} // first half of tries: no slope violations are tolerated, tracks must span the entire region
						bool const is_halfway((p1[dim] - seg_start) > 0.5*(seg_end - seg_start));

						if (is_halfway) { // we've reached the halfway point
							track.d[dim][1] = p1[dim]; // end the tracks before this segment starts
							if (track.get_sz_dim(dim) < 0.5*region.get_sz_dim(dim)) {valid = 0;} // must be at least half the full length to be accepted
							break;
						}
						else { // not yet halfway; can get here multiple times
							track.d[dim][0] = p2[dim]; // begin the tracks at the end of this segment
							track_segs.resize(segs_start); // clear any previous tracks
						}
					}
					bool const slope(p2.z < p1.z);
					unsigned const num_segs(unsigned(ceil(seg_len/max_seg_len)));
					vector3d const step_delta((p2 - p1)/num_segs);

					for (unsigned s = 0; s < num_segs; ++s) { // split smaller so that segments are within the shadow map bounds
						point const p1s(p1 + s*step_delta);
						track_segs.emplace_back(p1s, (p1s + step_delta), width, dim, slope, n);
					}
					p1[dim] = p2[dim];
					p2[dim] = min((p1[dim] + step_sz), seg_end);
				} // end while
				if (!valid) {track_segs.resize(segs_start); continue;} // if not valid, clear any partial segs and try another location
				tracks.push_back(track);
				dim_tracks[dim].push_back(track); // add to dim_tracks, but not to blockers, since we want roads to cross tracks
				break; // success
			} // for tries
		} // for n
		for (unsigned pass = 0; pass < 2; ++pass) { // flatten mesh after placing all tracks: regular, decrease_only
			for (auto const &t : track_segs) {hq.flatten_for_road(t, TRACKS_WIDTH*city_params.road_border, 0, (pass == 1));}
		}
		for (auto const &r : to_re_flatten) { // flatten out the previous segments again to avoid piles of dirt
			hq.flatten_for_road(r, city_params.road_border, 0, 1); // decrease_only=1
		}
		cout << "tracks: " << tracks.size() << ", track segments: " << track_segs.size() << endl;
	}

	void calc_bcube_from_roads() { // Note: ignores isecs, plots, and bridges, which should be bounded by roads
		if (roads.empty()) return; // no roads (assumes also no tracks)
		bcube = calc_cubes_bcube(roads);
		for (auto t = tracks.begin(); t != tracks.end(); ++t) {bcube.union_with_cube(*t);}
	}
private:
	int find_conn_int_seg(cube_t const &c, bool dim, bool dir) const {
		float const min_seg_len(1.0*city_params.road_width);

		for (unsigned i = 0; i < segs.size(); ++i) {
			road_seg_t const &s(segs[i]);
			if (s.dim == dim) continue; // not perp dim
			if (s.d[dim][dir] != bcube.d[dim][dir]) continue; // not on edge of road grid
			if (s.d[!dim][1] < c.d[!dim][0] || s.d[!dim][0] > c.d[!dim][1]) continue; // no overlap/projection in other dim
			// c contained in segment in other dim with enough padding (min road width) on each side
			if (c.d[!dim][0] > s.d[!dim][0]+min_seg_len && c.d[!dim][1] < s.d[!dim][1]-min_seg_len) return i; // this is the one we want
			return -1; // partial overlap in other dim, can't split, fail
		} // for i
		return -1; // not found
	}
	int find_3way_int_at(cube_t const &c, bool dim, bool dir) const {
		float const cube_cent(c.get_center_dim(!dim));
		assert(bcube.d[!dim][1] > cube_cent && bcube.d[!dim][0] < cube_cent); // c must overlap bcube in !dim
		//float dmin(0.0);
		int ret(-1);

		for (unsigned i = 0; i < isecs[1].size(); ++i) {
			road_isec_t const &isec(isecs[1][i]);
			if (isec.d[dim][dir] != bcube.d[dim][dir]) continue; // not on edge of road grid
			if (isec.d[!dim][1] > cube_cent && isec.d[!dim][0] < cube_cent) return i; // this is the one we want (early terminate case)
			//float const int_cent(isec.get_center_dim(!dim)), dist(fabs(int_cent - cube_cent));
			//if (ret < 0 || dist < dmin) {ret = i; dmin = dist;} // update if closer
		} // for i
		return ret; // not found if ret is still at -1
	}
	template<typename T> static void do_road_align(vector<T> &v, float from, float to, bool dim) {
		for (auto i = v.begin(); i != v.end(); ++i) {
			for (unsigned d = 0; d < 2; ++d) { // low, high edge
				if (i->d[dim][d] == from) {i->d[dim][d] = to;}
			}
		} // for i
	}
	bool align_isec3_to(unsigned int3_ix, cube_t const &c, bool dim) {
		assert(int3_ix < isecs[1].size());
		cube_t const src(isecs[1][int3_ix]); // deep copy so that it's not changed below
		//cout << c.d[!dim][0] << " " << c.d[!dim][1] << " " << src.d[!dim][0] << " " << src.d[!dim][1] << " " << (c.d[!dim][0] - src.d[!dim][0]) << endl; // TESTING
		if (c.d[!dim][0] == src.d[!dim][0]) return 0; // already aligned - done

		for (unsigned d = 0; d < 2; ++d) { // low, high
			do_road_align(roads, src.d[!dim][d], c.d[!dim][d], !dim);
			do_road_align(segs,  src.d[!dim][d], c.d[!dim][d], !dim);
			do_road_align(plots, src.d[!dim][d], c.d[!dim][d], !dim);
			for (unsigned i = 0; i < 3; ++i) {do_road_align(isecs[i], src.d[!dim][d], c.d[!dim][d], !dim);}
		}
		return 1;
	}
	void make_4way_int(unsigned int3_ix, bool dim, bool dir, unsigned conn_to_city, int road_ix) { // turn a 3-way intersection into a 4-way intersection for a connector road
		assert(int3_ix < isecs[1].size());
		road_isec_t &isec(isecs[1][int3_ix]); // bbox doesn't change, only conn changes
		isec.make_4way(conn_to_city); // all connected
		isec.rix_xy[2*dim + dir] = road_ix; // Note: should be negative (connector road)
		isecs[2].push_back(isec); // add as 4-way intersection
		isecs[1][int3_ix] = isecs[1].back(); // remove original 3-way intersection
		isecs[1].pop_back();
	}

	struct road_ixs_t {
		vector<unsigned> seg_ixs, isec_ixs[3][2]; // {2-way, 3-way, 4-way} x {X, Y}
	};
	template<typename T> int search_for_adj(vector<T> const &v, vector<unsigned> const &ixs, cube_t const &bcube, bool dim, bool dir) const {
		for (auto i = ixs.begin(); i != ixs.end(); ++i) {
			assert(*i < v.size());
			cube_t const &c(v[*i]);
			if (c.d[dim][!dir] != bcube.d[dim][dir]) continue; // no shared edge
			if (c.d[!dim][0] != bcube.d[!dim][0] || c.d[!dim][1] != bcube.d[!dim][1]) continue; // no shared edge in other dim
			return *i; // there can be only one
		} // for i
		return -1; // not found
	}
	vector<unsigned> const &get_segs_connecting_to_city(unsigned city) const {
		assert(city < city_to_seg.size());
		return city_to_seg[city];
	}
public:
	void calc_ix_values(vector<road_network_t> const &road_networks, road_network_t const &global_rn, unsigned &global_plot_id) {
		plot_id_offset  = global_plot_id; // cache offset into global plots vector
		global_plot_id += plots.size();
		// now that the segments and intersections are in order, we can fill in the IDs; first, create a mapping from road to segments and intersections
		bool const is_global_rn(&global_rn == this);
		assert(road_to_city.size() == (is_global_rn ? roads.size() : 0));
		vector<road_ixs_t> by_ix(roads.size()); // maps road_ix to list of seg_ix values
		vector<unsigned> all_ixs;

		if (is_global_rn) {
			unsigned num_cities(0);
			for (auto r = road_to_city.begin(); r != road_to_city.end(); ++r) {
				for (unsigned d = 0; d < 2; ++d) {if (r->id[d] != CONN_CITY_IX) {max_eq(num_cities, r->id[d]+1);}}
			}
			city_to_seg.resize(num_cities);
		}
		for (unsigned i = 0; i < segs.size(); ++i) {
			unsigned const ix(segs[i].road_ix);
			assert(ix < by_ix.size());
			assert(roads[ix].dim == segs[i].dim);
			by_ix[ix].seg_ixs.push_back(i);
		}
		for (unsigned n = 0; n < 3; ++n) { // {2-way, 3-way, 4-way}
			for (unsigned i = 0; i < isecs[n].size(); ++i) {
				road_isec_t const &isec(isecs[n][i]);

				for (unsigned d = 0; d < 2; ++d) { // {x, y}
					for (unsigned e = 0; e < 2; ++e) {
						int ix(isec.rix_xy[2*d + e]);

						if (ix < 0) { // global connector road
							ix = decode_neg_ix(ix);
							assert((unsigned)ix < global_rn.roads.size()); // connector road, nothing else to do here?
						}
						else if (e == 1 && ix == isec.rix_xy[2*d]) {} // same local road in both dirs, skip
						else {
							assert((unsigned)ix < roads.size());
							assert(roads[ix].dim == (d != 0));
							by_ix[ix].isec_ixs[n][d].push_back(i);
						}
					} // for e
				} // for d
				if (isec.conn_to_city >= 0) {cix_to_isec[isec.conn_to_city] = &isec;}
			} // for i
		} // for n

		// next, connect segments and intersections together by index, using roads as an acceleration structure
		for (unsigned i = 0; i < segs.size(); ++i) {
			road_seg_t &seg(segs[i]);
			road_ixs_t const &rix(by_ix[seg.road_ix]);

			for (unsigned dir = 0; dir < 2; ++dir) { // dir
				bool found(0);
				int const seg_ix(search_for_adj(segs, rix.seg_ixs, seg, seg.dim, (dir != 0)));
				if (seg_ix >= 0) {assert(seg_ix != (int)i); seg.conn_ix[dir] = seg_ix; seg.conn_type[dir] = TYPE_RSEG; found = 1;} // found segment

				for (unsigned n = 0; n < 3; ++n) { // 2-way, 3-way, 4-way
					int const isec_ix(search_for_adj(isecs[n], rix.isec_ixs[n][seg.dim], seg, seg.dim, (dir != 0)));
					if (isec_ix >= 0) {assert(!found); seg.conn_ix[dir] = isec_ix; seg.conn_type[dir] = (TYPE_ISEC2 + n); found = 1;} // found intersection
				}
				if (is_global_rn && !found) { // connection to a city
					assert(seg.road_ix < road_to_city.size());
					unsigned const city(road_to_city[seg.road_ix].id[dir]);
					assert(city != CONN_CITY_IX); // internal segments should be connected and not get here
					assert(city < road_networks.size());
					assert(city < city_to_seg.size());
					city_to_seg[city].push_back(i); // add segment ID
					road_network_t const &rn(road_networks[city]);

					for (unsigned n = 1; n < 3; ++n) { // search 3-way and 4-way intersections
						all_ixs.resize(rn.isecs[n].size());
						for (unsigned m = 0; m < all_ixs.size(); ++m) {all_ixs[m] = m;} // all sequential index values
						int const isec_ix(rn.search_for_adj(rn.isecs[n], all_ixs, seg, seg.dim, (dir != 0)));
						if (isec_ix < 0) continue; // not be found
						seg.conn_ix  [dir] = isec_ix;
						seg.conn_type[dir] = TYPE_ISEC2 + n; // always connects to a 3-way or 4-way intersection within the city
						found = 1;
						break;
					} // for n
					assert(found); // must be found
					continue;
				}
				assert(found);
			} // for dir
		} // for i
		for (unsigned n = 0; n < 3; ++n) { // 2-way, 3-way, 4-way
			for (unsigned i = 0; i < isecs[n].size(); ++i) {
				road_isec_t &isec(isecs[n][i]);

				for (unsigned d = 0; d < 4; ++d) { // {-x, +x, -y, +y}
					if (!(isec.conn & (1<<d))) continue; // no connection in this position
					unsigned const dim(d>>1), dir(d&1);
					int const ix(isec.rix_xy[d]);

					if (ix < 0) { // global connector road
						vector<unsigned> const &seg_ids(global_rn.get_segs_connecting_to_city(city_id));
						assert(!seg_ids.empty());
						int const seg_ix(global_rn.search_for_adj(global_rn.segs, seg_ids, isec, (dim != 0), (dir != 0))); // global conn segment
						assert(seg_ix >= 0); // must be found
						if (seg_ix >= 0) {isec.conn_ix[d] = encode_neg_ix(seg_ix);}
					}
					else { // local segment
						int const seg_ix(search_for_adj(segs, by_ix[ix].seg_ixs, isec, (dim != 0), (dir != 0))); // always connects to a road segment
						assert(seg_ix >= 0); // must be found
						if (seg_ix >= 0) {isec.conn_ix[d] = seg_ix;}
					}
				} // for d
			} // for i
		} // for n
		for (auto r = roads.begin(); r != roads.end(); ++r) {tot_road_len += r->get_length();} // calculate tot_road_len
	}
	bool check_valid_conn_intersection(cube_t const &c, bool dim, bool dir, bool is_4_way) const {
		return (is_4_way ? (find_3way_int_at(c, dim, dir) >= 0) : (find_conn_int_seg(c, dim, dir) >= 0));
	}
	void insert_conn_intersection(cube_t const &c, bool dim, bool dir, unsigned grn_rix, unsigned dest_city_id, bool is_4_way) { // Note: dim is the dimension of the connector road
		assert(dest_city_id != city_id); // not connected to self

		if (is_4_way) {
			int const int3_ix(find_3way_int_at(c, dim, dir));
			if (int3_ix < 0) {cout << TXT(dim) << TXT(dir) << TXT(bcube.str()) << TXT(c.str()) << endl;}
			assert(int3_ix >= 0); // must be found
			align_isec3_to(int3_ix, c, dim);
			make_4way_int(int3_ix, dim, dir, dest_city_id, encode_neg_ix(grn_rix));
		}
		else {
			int const seg_id(find_conn_int_seg(c, dim, dir));
			assert(seg_id >= 0 && (unsigned)seg_id < segs.size());
			segs.push_back(segs[seg_id]); // clone the segment first
			road_seg_t &seg(segs[seg_id]);
			assert(seg.road_ix < roads.size() && roads[seg.road_ix].dim != dim); // sanity check
			seg        .d[!dim][1] = c.d[!dim][0]; // low part
			segs.back().d[!dim][0] = c.d[!dim][1]; // high part
			cube_t ibc(seg); // intersection bcube
			ibc.d[!dim][0] = c.d[!dim][0]; // copy width from c
			ibc.d[!dim][1] = c.d[!dim][1];
			uint8_t const conns[4] = {7, 11, 13, 14};
			int const other_rix(encode_neg_ix(grn_rix)); // make negative
			bool const add_stoplights = 1; // always true for connector roads or !is_residential?
			isecs[1].emplace_back(ibc, (dim ? seg.road_ix : (int)other_rix), (dim ? other_rix : (int)seg.road_ix), conns[2*(!dim) + dir], true, add_stoplights, dest_city_id); // 3-way
		}
	}

	// global connector road functions
	float create_connector_road(cube_t const &bcube1, cube_t const &bcube2, vect_cube_t &blockers, road_network_t *rn1, road_network_t *rn2, unsigned city1, unsigned city2,
		unsigned dest_city_id1, unsigned dest_city_id2, city_road_connector_t &crc, float road_width, float conn_pos, bool dim, bool check_only, bool is_4_way1, bool is_4_way2)
	{
		assert(city1 != city2 || city1 == CONN_CITY_IX); // only holds for 2-segment connector roads
		bool const dir(bcube1.d[dim][0] < bcube2.d[dim][0]);
		if (dir == 0) {swap(city1, city2);} // make {lo, hi}
		point p1, p2;
		p1.z = bcube1.z2();
		p2.z = bcube2.z2();
		p1[!dim] = p2[!dim] = conn_pos;
		p1[ dim] = bcube1.d[dim][ dir];
		p2[ dim] = bcube2.d[dim][!dir];
		bool const slope((p1.z < p2.z) ^ dir);
		road_t const road(p1, p2, road_width, dim, slope, roads.size());
		float const road_len(road.get_length()), delta_z(road.dz());
		assert(road_len > 0.0 && delta_z >= 0.0);
		heightmap_query_t &hq(crc.hq);

		if (delta_z/road_len > city_params.max_road_slope) { // slope is too high (split segments will have even higher slopes)
			if (!check_only) {cout << TXT(dim) << TXT(road_len) << TXT(delta_z) << TXT(bcube1.str()) << TXT(bcube2.str()) << TXT(p1.str()) << TXT(p2.str()) << endl;}
			assert(check_only);
			return -1.0;
		}
		unsigned const x1(hq.get_x_pos(road.x1())), y1(hq.get_y_pos(road.y1())), x2(hq.get_x_pos(road.x2())), y2(hq.get_y_pos(road.y2()));

		if (check_only) { // only need to do these checks in this case
			if (rn1 && !rn1->check_valid_conn_intersection(road, dim,  dir, is_4_way1)) return -1.0; // invalid, don't make any changes
			if (rn2 && !rn2->check_valid_conn_intersection(road, dim, !dir, is_4_way2)) return -1.0;

			for (auto b = blockers.begin(); b != blockers.end(); ++b) {
				if ((rn1 && b->contains_cube(bcube1)) || (rn2 && b->contains_cube(bcube2))) continue; // skip current cities
				// create an intersection if blocker is a road, and happens to be the same elevation?
				if (b->intersects_xy(road)) return -1.0; // bad intersection, fail
			}
			if (hq.any_underwater(x1, y1, x2+1, y2+1)) return -1.0; // underwater (Note: bounds check is done here)
		}
		if (!check_only) { // create intersections and add blocker
			unsigned const grn_rix(roads.size()); // may be wrong end of connector, but doesn't matter?
			if (rn1) {rn1->insert_conn_intersection(road, dim,  dir, grn_rix, dest_city_id2, is_4_way1);}
			if (rn2) {rn2->insert_conn_intersection(road, dim, !dir, grn_rix, dest_city_id1, is_4_way2);}
			float const blocker_padding(max(city_params.road_spacing, 2.0f*city_params.road_border*max(DX_VAL, DY_VAL)));
			blockers.push_back(road);
			blockers.back().expand_by(blocker_padding); // add extra padding
			conn_roads.push_back(road);
		}
		if (road_len <= city_params.conn_road_seg_len) { // simple single road segment case
			if (!check_only) {
				roads.push_back(road);
				road_to_city.emplace_back(city1, city2);
			}
			// Note: no bridges here, but could add them
			return hq.flatten_sloped_region(x1, y1, x2, y2, road.d[2][slope]-ROAD_HEIGHT, road.d[2][!slope]-ROAD_HEIGHT, dim, city_params.road_border, 0, 0, check_only);
		}
		if (!crc.segment_road(road, check_only)) return -1.0;
		float tot_dz(0.0);
		bool last_was_bridge(0), last_was_tunnel(0);
		vector<flatten_op_t> replay_fops;

		for (auto s = crc.segments.begin(); s != crc.segments.end(); ++s) {
			if (s->z2() < s->z1()) {swap(s->z2(), s->z1());} // swap zvals if needed
			assert(s->is_normalized());
			bridge_t bridge(*s);
			tunnel_t tunnel(*s);
			tot_dz += hq.flatten_for_road(*s, city_params.road_border, check_only, 0, (last_was_bridge ? nullptr : &bridge), (last_was_tunnel ? nullptr : &tunnel));
			replay_fops.push_back(hq.last_flatten_op);
				
			if (!check_only) {
				if (bridge.make_bridge) {bridges.push_back(bridge); s->register_bridge_or_tunnel(bridge, 1);}
				if (tunnel.enabled())   {tunnels.push_back(tunnel); s->register_bridge_or_tunnel(bridge, 0);}
				roads.push_back(*s);
				road_to_city.emplace_back(city1, city2); // Note: city index is specified even for internal (non-terminal) roads
			}
			last_was_bridge = bridge.make_bridge; // Note: conservative; used to prevent two consecutive bridges with no (or not enough) mesh in between
			last_was_tunnel = tunnel.enabled(); // same thing for tunnels
		} // for s
		if (!check_only) { // post-flatten pass to fix up dirt at road joints - doesn't help much
			for (auto f = replay_fops.begin(); f != replay_fops.end(); ++f) { // replay the same series of operations; Note that bridge and tunnel segments have been cached
				hq.flatten_sloped_region(f->x1, f->y1, f->x2, f->y2, f->z1, f->z2, f->dim, f->border, f->skip_six, f->skip_eix, 0, 1);
			}
		}
		return tot_dz; // success
	}
	void create_connector_bend(cube_t const &int_bcube, bool dx, bool dy, unsigned road_ix_x, unsigned road_ix_y) {
		uint8_t const conns[4] = {6, 5, 10, 9};
		isecs[0].emplace_back(int_bcube, road_ix_x, road_ix_y, conns[2*dy + dx], true, false); // add_stoplights=0
		//blockers.push_back(int_bcube); // ???
	}
	void split_connector_roads(float road_spacing) { // required for correct shadow maps, since default segments may be too long
		// Note: here we use segs, maybe 2-way isecs for bends, but not plots
		for (auto r = roads.begin(); r != roads.end(); ++r) {
			bool const d(r->dim), slope(r->slope);
			float const len(r->get_length());
			unsigned const rix(r->road_ix); // not (r - roads.begin())
			if (len <= road_spacing) {segs.emplace_back(*r, rix); continue;} // single segment road
			assert(len > 0.0);
			unsigned const num_segs(ceil(len/road_spacing));
			float const seg_len(len/num_segs), z1(r->d[2][slope]), z2(r->d[2][!slope]); // use fixed-length segments
			assert(seg_len <= road_spacing);
			road_t c(*r); // start by copying the road's bcube
				
			for (unsigned n = 0; n < num_segs; ++n) {
				c.d[d][1] = ((n+1 == num_segs) ? r->d[d][1] : (c.d[d][0] + seg_len)); // make sure it ends exactly at the correct location
				for (unsigned e = 0; e < 2; ++e) {c.d[2][e] = z1 + (z2 - z1)*((c.d[d][e] - r->d[d][0])/len);} // interpolate road height across segments
				if (c.z2() < c.z1()) {swap(c.z2(), c.z1());} // swap zvals if needed
				assert(c.is_normalized());
				segs.emplace_back(c, rix);
				if (c.has_bridge) {segs.back().has_bridge = has_bcube_int_xy(c, bridges);} // check if bridge overlaps this segment
				if (c.has_tunnel) {segs.back().has_tunnel = has_bcube_int_xy(c, tunnels);} // check if tunnel overlaps this segment
				c.d[d][0] = c.d[d][1]; // shift segment end point
			} // for n
		} // for r
	}
	void finalize_bridges_and_tunnels() {
		for (bridge_t &b : bridges) {b.add_streetlights();}
		for (tunnel_t &t : tunnels) {t.add_streetlights();}
	}

	void gen_tile_blocks() {
		tile_blocks.clear(); // should already be empty?
		tile_to_block_map.clear();
		add_tile_blocks(segs,       tile_to_block_map, TYPE_RSEG);
		add_tile_blocks(plots,      tile_to_block_map, TYPE_PLOT);
		add_tile_blocks(track_segs, tile_to_block_map, TYPE_TRACKS);
		for (unsigned i = 0; i < 3; ++i) {add_tile_blocks(isecs[i], tile_to_block_map, (TYPE_ISEC2 + i));}
		plot_xy.gen_adj_plots(plots);
		//cout << "tile_to_block_map: " << tile_to_block_map.size() << ", tile_blocks: " << tile_blocks.size() << endl;
	}
	void gen_parking_lots_and_place_objects(vector<car_t> &cars, bool have_cars, bool &have_plot_dividers) {
		city_obj_placer = city_obj_placer_t(); // clear; should be empty anyway, since city_obj_placer is not reused
		city_obj_placer.set_plot_subdiv_sz(get_plot_subdiv_sz());
		city_obj_placer.add_city_ug_elevator_entrances(uges); // Note: can clear ug_elev_entrances after this, but it's not needed
		city_obj_placer.gen_parking_and_place_objects(plots, plot_colliders, cars, roads, isecs, bcube, city_id, have_cars, is_residential, !streetlights.empty());
		add_tile_blocks(city_obj_placer.parking_lots, tile_to_block_map, TYPE_PARK_LOT); // need to do this later, after gen_tile_blocks()
		add_tile_blocks(city_obj_placer.driveways,    tile_to_block_map, TYPE_DRIVEWAY);
		city_obj_placer.remap_parking_lot_ixs(); // required after sorting parking_lots
		tile_to_block_map.clear(); // no longer needed
		city_obj_placer.finalize_streetlights_and_power(*this, plot_colliders);
		for (auto i = plot_colliders.begin(); i != plot_colliders.end(); ++i) {sort(i->begin(), i->end(), [](cube_t const &a, cube_t const &b) {return (a.x1() < b.x1());});}
		have_plot_dividers |= !city_obj_placer.has_plot_dividers();
	}
	void add_streetlights() {
		streetlights.clear();
		streetlights.reserve(4*plots.size()); // one on each side of each plot
		// spacing from light pos to plot edge, relative to plot size (placed just outside the plot, so spacing is negative)
		float const b(-(SIDEWALK_WIDTH*city_params.road_width - 0.5f*streetlight_ns::get_streetlight_pole_radius())/city_params.road_spacing), a(1.0 - b);
		assert(plot_colliders.size() == plots.size());

		for (auto i = plots.begin(); i != plots.end(); ++i) {
			unsigned const plot_ix(i - plots.begin());
			streetlights.emplace_back(point((a*i->x1() + b*i->x2()), (0.75*i->y1() + 0.25*i->y2()), i->z2()), -plus_x, 0, plot_ix); // left   edge one   quarter  up
			streetlights.emplace_back(point((a*i->x2() + b*i->x1()), (0.25*i->y1() + 0.75*i->y2()), i->z2()),  plus_x, 0, plot_ix); // right  edge three quarters up
			streetlights.emplace_back(point((0.25*i->x1() + 0.75*i->x2()), (a*i->y1() + b*i->y2()), i->z2()), -plus_y, 0, plot_ix); // bottom edge three quarters right
			streetlights.emplace_back(point((0.75*i->x1() + 0.25*i->x2()), (a*i->y2() + b*i->y1()), i->z2()),  plus_y, 0, plot_ix); // top    edge one   quarter  right
		}
		sort_streetlights_by_yx();
	}
	bool update_depth_if_underwater(point const &pos, float &depth) const { // Note: pos is in building space
		return city_obj_placer.update_depth_if_underwater(pos, depth);
	}
	void get_road_bcubes(vect_cube_t &bcubes) const {
		get_all_bcubes(roads,  bcubes);
		get_all_bcubes(tracks, bcubes);
	}
	float get_plot_subdiv_sz() const {return ((is_residential && city_params.assign_house_plots) ? 0.8*get_max_house_size() : 0.0);}

	void get_plot_zones(vect_city_zone_t &zones) const { // Note: z-values of cubes indicate building height ranges
		if (plots.empty()) return; // connector road city
		unsigned const start(zones.size());
		float const plot_subdiv_sz(get_plot_subdiv_sz());
		unsigned cur_global_plot_ix(plot_id_offset), max_floors(0); // max_floors of 0 is unlimited

		for (auto i = plots.begin(); i != plots.end(); ++i, ++cur_global_plot_ix) { // capture all plot zones, even parks (needed for pedestrians)
			if (plot_subdiv_sz > 0.0 && !i->is_park) { // split into smaller plots for each house
				if (city_obj_placer.subdivide_plot_for_residential(*i, roads, plot_subdiv_sz, cur_global_plot_ix, city_id, zones)) continue;
			}
			zones.emplace_back(*i, 0.0, i->is_park, is_residential, 0, 0, cur_global_plot_ix, city_id, max_floors); // cube, zval, park, res, sdir, capacity, ppix, cix, max_floors
		}
		vector3d const city_radius(0.5*bcube.get_size());
		point const city_center(bcube.get_cube_center());

		for (auto i = zones.begin()+start; i != zones.end(); ++i) { // set zvals to control building height range, higher in city center
			point const center(i->get_cube_center());
			float const dx(fabs(center.x - city_center.x)/city_radius.x); // 0 at city center, 1 at city perimeter
			float const dy(fabs(center.y - city_center.y)/city_radius.y); // 0 at city center, 1 at city perimeter
			float const hval(1.0 - max(dx*dx, dy*dy)); // square to give higher weight to larger height ranges
			i->zval = i->z2(); // capture z2 value for use in setting building height
			i->z1() = max(0.0, (hval - 0.25)); // bottom of height range
			i->z2() = min(1.0, (hval + 0.25)); // bottom of height range
		} // for i
	}
	bool check_road_sphere_coll(point const &pos, float radius, bool include_intersections, bool xy_only, bool exclude_bridges_and_tunnels) const {
		if (roads.empty()) return 0;
		point const query_pos(pos - get_camera_coord_space_xlate()); // convert from camera space to building space
		if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
			
		if (check_bcubes_sphere_coll(roads, query_pos, radius, xy_only)) { // collision with a road
			if (!exclude_bridges_and_tunnels ||
				!(check_bcubes_sphere_coll(bridges, query_pos, radius, xy_only) ||
				  check_bcubes_sphere_coll(tunnels, query_pos, radius, xy_only))) return 1; // ignore collisions with bridges and tunnels
		}
		if (include_intersections) { // used for global road network
			for (unsigned i = 0; i < 3; ++i) { // {2-way, 3-way, 4-way}
				if (check_bcubes_sphere_coll(isecs[i], query_pos, radius, xy_only)) return 1;
			}
		}
		if (check_bcubes_sphere_coll(tracks, query_pos, radius, xy_only)) return 1; // collision with a track
		return 0;
	}
	// Note: returns cubes in local pos space
	void get_roads_in_region(cube_t const &region, vect_cube_t &out, vect_cube_t &out_bt) const {
		if (roads.empty()) return;
		vector3d const xlate(get_camera_coord_space_xlate());
		cube_t const query_region(region - xlate);
		if (!bcube.intersects_xy(query_region)) return;
		get_bcubes_region_coll_xy(roads, out, query_region, xlate);	
		// include global road network intersections
		for (unsigned i = 0; i < 3; ++i) {get_bcubes_region_coll_xy(isecs[i], out, query_region, xlate);} // {2-way, 3-way, 4-way}
		get_bcubes_region_coll_xy(bridges, out_bt, query_region, xlate);
		get_bcubes_region_coll_xy(tunnels, out_bt, query_region, xlate);
		get_bcubes_region_coll_xy(tracks,  out,    query_region, xlate);
	}
	bool proc_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float dist, float radius, float prev_frame_zval, vector3d *cnorm) const { // pos in camera space
		if (!moving_sphere_cube_intersect_xy(pos, p_last, (bcube + xlate), dist, radius)) return 0;
		bool plot_coll(0);
			
		if (!plots.empty()) {
			float const max_obj_z(bcube.z1() + radius);
			if (pos.z < max_obj_z) {pos.z = max_obj_z; plot_coll = 1;} // make sure the sphere is above the city road/plot surface
		}
		for (unsigned n = 1; n < 3; ++n) { // intersections, possibly with stoplights (3-way, 4-way)
			for (auto i = isecs[n].begin(); i != isecs[n].end(); ++i) {
				if (i->proc_sphere_coll(pos, p_last, radius, xlate, dist, cnorm)) return 1;
			}
		}
		for (auto i = bridges.begin(); i != bridges.end(); ++i) {
			if (i->proc_sphere_coll(pos, p_last, radius, prev_frame_zval, xlate, cnorm)) return 1;
		}
		for (auto i = tunnels.begin(); i != tunnels.end(); ++i) {
			if (i->proc_sphere_coll(pos, p_last, radius, prev_frame_zval, xlate, cnorm)) return 1;
		}
		// test city objects before streetlights so that player doesn't fall through a walkway when over a streetlight
		if (city_obj_placer.proc_sphere_coll(pos, p_last, xlate, radius, cnorm)) return 1;

		if ((pos.z - xlate.z - radius) < (bcube.z2() + streetlight_ns::get_streetlight_height())) { // below the level of the streetlights
			if (proc_streetlight_sphere_coll(pos, radius, xlate, cnorm)) return 1;
		}
		if (0 && plot_coll) { // no other collisions - return collision with plot or road - doesn't work correctly for bouncing balls
			if (cnorm) {*cnorm = plus_z;}
			return 1;
		}
		return 0;
	}
	bool line_intersect(point const &p1, point const &p2, float &t) const { // Note: xlate has already been applied
		bool ret(0);
			
		if (get_bcube_inc_stoplights_and_streetlights().line_intersects(p1, p2)) { // z2 too small for streetlights?
			for (unsigned n = 1; n < 3; ++n) { // intersections, possibly with stoplights (3-way, 4-way)
				for (auto i = isecs[n].begin(); i != isecs[n].end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
			}
			for (auto i = bridges.begin(); i != bridges.end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
			for (auto i = tunnels.begin(); i != tunnels.end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
			ret |= line_intersect_streetlights(p1, p2, t);
		}
		ret |= city_obj_placer.line_intersect(p1, p2, t);
		return ret;
	}
	bool is_invalid_placement_for_cube(cube_t const &c) const { // Note: c is in global space
		if (!bcube.intersects_xy(c))    return 0; // not in this city
		if (has_bcube_int_xy(c, roads)) return 1; // intersects or above a road
		if (plots.empty()) return 0; // global connector road

		for (unsigned plot_ix = 0; plot_ix < plots.size(); ++plot_ix) {
			if (!plots[plot_ix].contains_cube_xy(c)) continue;
			if (plot_ix < plot_colliders.size() && has_bcube_int(c, plot_colliders[plot_ix])) return 1;
			if (city_obj_placer.intersects_parking_lot(c)) return 1; // parking lots and driveways may not have been added at this point
			return 0; // fully contained in this plot and not intersecting any colliders
		}
		return 1; // not contained in a plot
	}
	void add_plot_cut(cube_t const &cut) {
		if (bcube.intersects_xy(cut)) {plot_cuts.push_back(cut);}
	}
	void add_city_ug_elevator_entrance(ug_elev_info_t const &uge) {
		if (!bcube.intersects_xy(uge.entrance)) return; // wrong city
		assert(bcube.contains_cube_xy(uge.entrance)); // can't partially overlap
		uges.push_back(uge);
		add_cube_to_plot_colliders(uge.entrance);
	}
	void add_cube_to_plot_colliders(cube_t const &c) { // Note: must be before plot_colliders sorting in gen_parking_lots_and_place_objects()
		for (unsigned plot_ix = 0; plot_ix < plots.size(); ++plot_ix) {
			if (!plots[plot_ix].intersects_xy(c)) continue;
			plot_colliders[plot_ix].push_back(c); // or insert in the correct place?
			if (plots[plot_ix].contains_cube_xy(c)) break; // should always be true
		}
	}
	bool check_mesh_disable(point const &pos, float radius) const { // Note: pos is in camera space
		if (tunnels.empty()) return 0;
		point const query_pos(pos - get_camera_coord_space_xlate());
		cube_t query_region; query_region.set_from_sphere(query_pos, radius); // actually a cube, not a sphere
		if (!bcube.intersects_xy(query_region)) return 0;

		for (tunnel_t const &t : tunnels) {
			if (t.check_mesh_disable(query_region)) return 1;
		}
		return 0;
	}
	bool tile_contains_tunnel(cube_t const &tile_bcube) const { // Note: cube is in global space
		if (tunnels.empty() || !bcube.intersects_xy(tile_bcube)) return 0;

		for (tunnel_t const &t : tunnels) {
			if (t.intersects_xy(tile_bcube)) return 1;
		}
		return 0;
	}
	bool point_in_tunnel(point const &pos) const { // Note: pos is in global space
		if (tunnels.empty() || !bcube.contains_pt_xy(pos)) return 0;

		for (tunnel_t const &t : tunnels) {
			if (t.contains_pt(pos)) return 1; // Note: checks z-val
		}
		return 0;
	}
	bool cube_intersect_tunnel(cube_t const &c) const { // Note: cube is in global space
		if (tunnels.empty() || !bcube.intersects_xy(c)) return 0;

		for (tunnel_t const &t : tunnels) {
			if (t.get_tunnel_bcube().intersects(c)) return 1; // use outer bounding cube
		}
		return 0;
	}
	bool cube_int_underground_obj(cube_t const &c) const {
		return (c.intersects_xy(bcube)) && city_obj_placer.cube_int_underground_obj(c);
	}
	void get_ponds_in_xy_range(cube_t const &range, vect_cube_t &ponds) const {
		if (range.intersects_xy(bcube)) {city_obj_placer.get_ponds_in_xy_range(range, ponds);}
	}
	bool choose_pt_in_park(point &park_pos, rand_gen_t &rgen) const {
		if (parks.empty()) return 0;
		cube_t const &park(parks[rgen.rand() % parks.size()]); // select a random park
		park_pos = rand_xy_pt_in_cube(park, get_sidewalk_width(), rgen);
		return 1;
	}
	void add_manhole(point const &pos, float radius) {
		bool is_over_road(0);

		for (road_t const &road : roads) {
			if (road.contains_pt_xy_exp(pos, -radius)) {is_over_road = 1; break;} // must fully contain manhole
		}
		city_obj_placer.add_manhole(pos, radius, is_over_road);
	}
	template<typename T> bool check_tile_group_contains_pt_xy(vector<T> const &objs, point const &pos, unsigned type) const {
		assert(type < NUM_RD_TYPES);
		if (objs.empty()) return 0;

		for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
			if (!b->bcube.contains_pt_xy(pos)) continue;
			range_pair_t const &rp(b->ranges[type]);
			for (unsigned i = rp.s; i < rp.e; ++i) {if (objs[i].contains_pt_xy(pos)) return 1;}
		}
		return 0;
	}
	template<typename T> bool cube_overlaps_tile_group_xy(vector<T> const &objs, cube_t const &c, unsigned type) const {
		assert(type < NUM_RD_TYPES);
		if (objs.empty()) return 0;

		for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
			if (!b->bcube.intersects_xy(c)) continue;
			range_pair_t const &rp(b->ranges[type]);
			for (unsigned i = rp.s; i < rp.e; ++i) {if (objs[i].intersects_xy(c)) return 1;}
		}
		return 0;
	}
	bool cube_overlaps_pl_or_dw_xy(cube_t const &c) const {
		return (cube_overlaps_tile_group_xy(city_obj_placer.parking_lots, c, TYPE_PARK_LOT) || cube_overlaps_tile_group_xy(city_obj_placer.driveways, c, TYPE_DRIVEWAY));
	}
	int get_color_at_xy(point const &pos, colorRGBA &color) const { // Note: return value is currently unused, but it could be used for something in the future
		// Note: query results are mutually exclusive since there's no overlap, so can early terminate on true
		if (!bcube.contains_pt_xy(pos)) return 0;
		if (check_vect_cube_contains_pt_xy(tunnels, pos)) {color = BROWN; return INT_ROAD;}
		if (city_obj_placer.get_color_at_xy_pre_road(pos, color)) return INT_BUILDING; // walkway
			
		for (auto i = bridges.begin(); i != bridges.end(); ++i) {
			if (i->contains_pt_xy_exp(pos, 1.0*city_params.road_width)) {color = WHITE; return INT_ROAD;}
		}
		if (!conn_roads.empty()) { // global_rn connector roads - use this vector because we only care about XY projection (not Z), and conn_roads is smaller than roads
			if (check_vect_cube_contains_pt_xy(conn_roads, pos)) {color = GRAY; return INT_ROAD;}
		}
		else {
			if (check_vect_cube_contains_pt_xy(roads,      pos)) {color = GRAY; return INT_ROAD;}
		}
		for (auto i = tracks.begin(); i != tracks.end(); ++i) {
			if (!i->contains_pt_xy(pos)) continue;
			float const width(i->get_sz_dim(!i->dim)), edge_dist(fabs(pos[!i->dim] - i->get_center_dim(!i->dim))), val(edge_dist/width);
			color = ((val > 0.2 && val < 0.3) ? GRAY : LT_BROWN);
			return INT_TRACK;
		}
		if (plots.empty()) { // check connector road bends
			if (check_vect_cube_contains_pt_xy(isecs[0], pos)) {color = GRAY; return INT_ROAD;} // 2-way intersections
		}
		if (check_tile_group_contains_pt_xy(city_obj_placer.parking_lots, pos, TYPE_PARK_LOT)) {color = DK_GRAY; return INT_PARKING;}
		if (check_tile_group_contains_pt_xy(city_obj_placer.driveways,    pos, TYPE_DRIVEWAY)) {color = (is_residential ? LT_GRAY : colorRGBA(0.4, 0.4, 0.4)); return INT_PARKING;}
		if (city_obj_placer.get_color_at_xy(pos, color, 1)) {return INT_PLOT;} // hit a detail object, but still in a plot; skip objects in roads such as fire hydrants
			
		if (!plots.empty()) { // inside a city and not over a road - must be over a plot or park
			for (auto i = parks.begin(); i != parks.end(); ++i) {
				if (i->contains_pt_xy(pos)) {color = GREEN; return INT_PARK;}
			}
			color = (is_residential ? DK_GREEN : colorRGBA(0.65, 0.65, 0.65, 1.0)); // grass or concrete
			return INT_PLOT;
		}
		return INT_NONE;
	}
	bool cube_overlaps_road_xy(cube_t const &c) const {
		// can we use conn_roads here for global_rn?
		for (auto i = roads.begin(); i != roads.end(); ++i) {if (i->intersects(c)) return 1;}
		return 0;
	}
	void get_occluders(vect_cube_t &occluders, vector3d const &xlate) const {
		if (bcube.contains_pt_xy(camera_pdu.pos - xlate)) {city_obj_placer.get_occluders(camera_pdu, xlate, occluders);} // only add if this city contains the camera
	}
	vector<bridge_t> const &get_bridges() const {return bridges;}
	bool have_animations() const {return city_obj_placer.have_animations();}
	static void set_road_normal_map  () {select_texture(get_texture_by_name("normal_maps/dirt_normal.jpg", 1), 5);}
	static void reset_road_normal_map() {bind_default_flat_normal_map();} // no normal map

	void draw(road_draw_state_t &dstate, bool shadow_only, bool is_connector_road) {
		city_obj_placer.draw_detail_objects(dstate, shadow_only); // always drawn; does its own VFC and distance test
		if (!empty()) {draw_roads_and_plots(dstate, shadow_only, is_connector_road);}
		dstate.end_cur_tile(); // once for all tiles, to draw shadow casters and untextured streetlights
	}
	void draw_roads_and_plots(road_draw_state_t &dstate, bool shadow_only, bool is_connector_road) {
		if (!dstate.check_cube_visible(get_bcube_inc_stoplights_and_streetlights(), 1.0)) return; // VFC/too far

		if (shadow_only) {
			if (!is_connector_road) { // connector road has no stoplights to cast shadows
				// Note: we can store the contents of qbd_sl in a VBO to avoid recreating it every frame for the shadow pass
				for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
					if (!dstate.check_cube_visible(b->bcube, 0.16)) continue; // VFC/too far; dist_scale=0.16
					
					for (unsigned i = 1; i < 3; ++i) { // isecs with stoplights (3-way, 4-way)
						dstate.draw_stoplights_and_street_signs(isecs[i], roads, b->ranges[TYPE_ISEC2 + i], city_id, 1);
					}
				}
			}
		}
		else { // regular draw pass
			// it would be nice to have raised sidewalks here, or possibly curbs;
			// however, they would need to apply to all tiles, including city plots and curved road segments, which makes it difficult;
			// also, we would have to adjust the height of pedestrians, fire hydrants, etc.
			bool const use_road_normal_maps(rain_wetness > 0.0); // use dirt normal map texture for rain effects

			for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
				if (!dstate.check_cube_visible(b->bcube)) continue; // VFC/too far
				dstate.begin_tile(b->bcube.get_cube_center());

				// if the player is in the basement, don't draw the plot over the basement stairs; the player can't see any of this anyway
				if (!player_in_basement || is_connector_road) {
					if (!plots.empty()) { // draw plots if not global connector road network
						cube_t const plot_exclude(get_cur_basement()); // clip out basement
						dstate.plot_cuts = plot_cuts; // reset to static plot cuts
						
						if (!plot_exclude.is_all_zeros() && plot_exclude.intersects_xy(b->bcube)) {
							b->quads[TYPE_PLOT].clear(); // clear and rebuild plot cache
							dstate.plot_cuts.push_back(plot_exclude);
						}
						if (is_residential) { // draw all plots with grass, using the park materials
							city_obj_placer.get_plot_cuts(b->bcube, dstate.plot_cuts); // inground swimming pools, etc.
							dstate.draw_city_region(plots, b->ranges[TYPE_PLOT], b->quads[TYPE_PLOT], TYPE_PARK, 1); // draw_all=1
						}
						else {
							dstate.draw_city_region(plots, b->ranges[TYPE_PLOT], b->quads[TYPE_PLOT], TYPE_PLOT); // concrete
							dstate.draw_city_region(plots, b->ranges[TYPE_PLOT], b->quads[TYPE_PARK], TYPE_PARK); // grass parks (stored as plots)
						}
						dstate.plot_cuts.clear();
					}
					if (use_road_normal_maps) {set_road_normal_map();} // set normal maps for roads, parking lots, and driveways
					dstate.draw_city_region(segs, b->ranges[TYPE_RSEG], b->quads[TYPE_RSEG], TYPE_RSEG); // road segments
					dstate.draw_city_region(city_obj_placer.parking_lots, b->ranges[TYPE_PARK_LOT], b->quads[TYPE_PARK_LOT], TYPE_PARK_LOT); // parking lots

					if (is_connector_road) { // draw road skirts
						dstate.draw_city_region(segs, b->ranges[TYPE_RSEG], b->quads[TYPE_ROAD_SKIRT], TYPE_ROAD_SKIRT); // same ranges as the road, but a different type
					}
					if (!city_obj_placer.driveways.empty()) {
						glPolygonOffset(-1.0, -1.0); // useful for avoiding z-fighting with grassy ground under driveways
						glEnable(GL_POLYGON_OFFSET_FILL);
						dstate.draw_city_region(city_obj_placer.driveways, b->ranges[TYPE_DRIVEWAY], b->quads[TYPE_DRIVEWAY], TYPE_DRIVEWAY); // driveways
						glDisable(GL_POLYGON_OFFSET_FILL);
					}
					if (use_road_normal_maps) {reset_road_normal_map();}
					dstate.draw_city_region(track_segs, b->ranges[TYPE_TRACKS], b->quads[TYPE_TRACKS], TYPE_TRACKS); // railroad tracks
				}
				if (use_road_normal_maps) {set_road_normal_map();}

				for (unsigned i = 0; i < 3; ++i) { // intersections (2-way, 3-way, 4-way)
					dstate.draw_city_region(isecs[i], b->ranges[TYPE_ISEC2 + i], b->quads[TYPE_ISEC2 + i], (TYPE_ISEC2 + i));
				}
				if (use_road_normal_maps) {reset_road_normal_map();}

				for (unsigned i = 1; i < 3; ++i) { // intersections (3-way, 4-way)
					dstate.draw_stoplights_and_street_signs(isecs[i], roads, b->ranges[TYPE_ISEC2 + i], city_id, 0);
				}
				dstate.end_cur_tile();
			} // for b
		} // end !shadow_only
		if (!is_connector_road) { // draw exterior edges (skirts) of city to cover the gap above the terrain
			// it's okay to draw under connector roads, and these are vertical and generally have no light or shadows
			dstate.draw_city_skirt(bcube, shadow_only);
		}
		draw_streetlights(dstate, shadow_only, 0);
			
		// draw bridges and tunnels; only in connector road network; bridges and tunnels are sparse/uncommon, so don't need to be batched by blocks
		for (auto b = bridges.begin(); b != bridges.end(); ++b) {
			dstate.draw_bridge(*b, shadow_only);
			b->draw_streetlights(dstate, shadow_only, 0);
		}
		for (auto t = tunnels.begin(); t != tunnels.end(); ++t) {
			dstate.draw_tunnel(*t, shadow_only);
			t->draw_streetlights(dstate, shadow_only, 1); // always_on=1
		}
	}
	void draw_transparent(road_draw_state_t &dstate) {
		city_obj_placer.draw_transparent_objects(dstate);
	}
	void add_city_lights(vector3d const &xlate, cube_t &lights_bcube) const { // for now, the only light sources added by the road network are city block streetlights
		add_streetlight_dlights(xlate, lights_bcube, 0);
		for (auto b = bridges.begin(); b != bridges.end(); ++b) {b->add_streetlight_dlights(xlate, lights_bcube, 0);}
		for (auto t = tunnels.begin(); t != tunnels.end(); ++t) {t->add_streetlight_dlights(xlate, lights_bcube, 1);} // always_on=1
		city_obj_placer.add_lights(xlate, lights_bcube);
	}

	// cars/peds
	static float get_car_lane_offset() {return CAR_LANE_OFFSET*city_params.road_width;}

	static unsigned gen_rand_city(rand_gen_t &rgen, vector<road_network_t> const &road_networks) {
		assert(!road_networks.empty());
		return rgen.rand()%road_networks.size();
	}
	static bool add_car_to_rns(car_t &car, rand_gen_t &rgen, vector<road_network_t> const &road_networks) {
		car.cur_city = gen_rand_city(rgen, road_networks);
		return road_networks[car.cur_city].add_car(car, rgen);
	}
	static bool gen_ped_pos(pedestrian_t &ped, rand_gen_t &rgen, vector<road_network_t> const &road_networks) {
		ped.city = gen_rand_city(rgen, road_networks);
		return road_networks[ped.city].gen_ped_pos(ped, rgen);
	}

	bool add_car(car_t &car, rand_gen_t &rgen) const {
		if (segs.empty()) return 0; // no segments to place car on
		vector3d const nom_car_size(city_params.get_nom_car_size());

		for (unsigned n = 0; n < 10; ++n) { // make 10 tries
			unsigned const seg_ix(rgen.rand()%segs.size());
			road_seg_t const &seg(segs[seg_ix]); // chose a random segment
			car.dim = seg.dim;
			car.dir = rgen.rand_bool();
			car.choose_max_speed(rgen);
			car.cur_road = seg.road_ix;
			car.cur_seg  = seg_ix;
			car.cur_road_type = TYPE_RSEG;
			point pos;
			float val1(seg.d[seg.dim][0] + 0.6f*nom_car_size.x), val2(seg.d[seg.dim][1] - 0.6f*nom_car_size.x);
			if (val1 >= val2) continue; // failed, try again (connector road junction?)
			pos[!seg.dim]  = seg.get_center_dim(!seg.dim); // center of road
			pos[!seg.dim] += ((car.dir ^ car.dim) ? -1.0 : 1.0)*get_car_lane_offset(); // place in right lane
			pos[ seg.dim]  = rgen.rand_uniform(val1, val2); // place at random pos in segment
			pos.z = seg.z2(); // place above road surface
			car.set_bcube(pos, nom_car_size);
			assert(get_road_bcube_for_car(car).contains_cube_xy(car.bcube)); // sanity check
			return 1; // success
		} // for n
		return 0; // failed
	}
	bool gen_ped_pos(pedestrian_t &ped, rand_gen_t &rgen) const {
		if (plots.empty()) return 0; // no plots to place car on

		for (unsigned n = 0; n < 100; ++n) { // make 100 tries
			unsigned const plot_id(rgen.rand()%plots.size()); // chose a random plot
			if (ped.try_place_in_plot(plots[plot_id], plot_colliders[plot_id], (plot_id + plot_id_offset), rgen)) return 1; // success
		}
		return 0; // failed
	}
	unsigned decode_plot_id(unsigned global_plot_id) const { // global => local
		assert(global_plot_id >= plot_id_offset);
		unsigned const plot_id(global_plot_id - plot_id_offset);
		assert(plot_id < plots.size());
		return plot_id;
	}
	unsigned encode_plot_id(unsigned local_plot_id) const {return (local_plot_id + plot_id_offset);} // local => global
	road_plot_t const &get_plot_from_global_id(unsigned global_plot_id) const {return plots         [decode_plot_id(global_plot_id)];}
	vect_cube_t const &get_colliders_for_plot (unsigned global_plot_id) const {return plot_colliders[decode_plot_id(global_plot_id)];}

	int get_global_plot_id_for_pos(point const &pos) const { // pos is in building space
		for (unsigned i = 0; i < plots.size(); ++i) {
			if (plots[i].contains_pt_xy(pos)) return encode_plot_id(i);
		}
		return -1; // not found
	}
	// plot = current plot, dest_plot = final destination plot; returns next plot adj to cur plot on path to dest_plot
	unsigned get_next_plot(unsigned global_plot, unsigned global_dest_plot, int exclude_plot) const {
		if (global_plot == global_dest_plot) {return global_plot;} // identity, at destination, no change
		unsigned const plot(decode_plot_id(global_plot)), dest_plot(decode_plot_id(global_dest_plot)); // convert to local space
		assert(plot < plots.size() && dest_plot < plots.size());
		int const cur_x(plots[plot].xpos), cur_y(plots[plot].ypos), dest_x(plots[dest_plot].xpos), dest_y(plots[dest_plot].ypos); // use int to avoid signed problems
		int const dx(dest_x - cur_x), dy(dest_y - cur_y);
		bool move_dir(abs(dx) > abs(dy)); // technically !move_dir
		unsigned const dir(move_dir ? ((dx < 0) ? 0 : 1) : ((dy < 0) ? 2 : 3));
		int const next_plot(plot_xy.get_adj(cur_x, cur_y, dir));
		assert(next_plot >= 0 && next_plot <= int(plots.size())); // must be found
		road_plot_t const &np(plots[next_plot]); // this next part is only for error checking
		if (move_dir) {assert(np.xpos == cur_x + ((dx < 0) ? -1 : 1)); assert(np.ypos == cur_y);}
		else          {assert(np.ypos == cur_y + ((dy < 0) ? -1 : 1)); assert(np.xpos == cur_x);}
		unsigned global_next_plot(next_plot + plot_id_offset); // convert back to global space
			
		if (exclude_plot == int(global_next_plot)) { // the selected plot has been excluded, choose a different plot
			move_dir ^= 1; // move in the other dimension - don't move in the wrong direction
			unsigned dir(0);

			if (dx != 0 && dy != 0) { // moving in the other dimension makes progress
				dir = (move_dir ? ((dx < 0) ? 0 : 1) : ((dy < 0) ? 2 : 3));	
			}
			else { // take a detour in a random direction
				static rand_gen_t rgen;
				bool rand_dir(rgen.rand_bool());
				dir = (move_dir ? (rand_dir ? 0 : 1) : (rand_dir ? 2 : 3));
					
				if (plot_xy.get_adj(cur_x, cur_y, dir) < 0) { // that direction's not valid, choose the other one
					rand_dir ^= 1;
					dir = (move_dir ? (rand_dir ? 0 : 1) : (rand_dir ? 2 : 3));
				}
			}
			int const next_plot(plot_xy.get_adj(cur_x, cur_y, dir));
			assert(next_plot >= 0 && next_plot <= int(plots.size())); // must be found
			global_next_plot = next_plot + plot_id_offset; // convert back to global space
		}
		return global_next_plot;
	}
	bool choose_dest_building(unsigned &global_plot, unsigned &building, rand_gen_t &rgen) const { // for pedestrians
		if (plots.empty()) return 0; // no plots
		global_plot = (rgen.rand() % plots.size()) + plot_id_offset;
		if (!select_building_in_plot(global_plot, rgen.rand(), building)) return 0; // no buildings in plot (maybe it's a park)
		return 1;
	}
	bool check_future_isec_stop(car_t &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
		if (car.is_stopped() || car.stopped_at_light) return 0;
		if (car.cur_road_type != TYPE_RSEG)           return 0; // not on a road segment; maybe already stopped
		road_seg_t const &seg(get_car_seg(car));
		unsigned const next_road_type(seg.conn_type[car.dir]);
		if (!is_isect(next_road_type)) return 0; // skip if this is an adjacent segment (for global connector roads)
		float const stop_dist(0.6f*car.get_length()); // distance we should begin decelerating
		float const dist_to_end(max((car.dir ? 1.0f : -1.0f)*(seg.d[car.dim][car.dir] - car.bcube.d[car.dim][car.dir]), 0.0f));
		if (dist_to_end > stop_dist) return 0; // not close enough
		unsigned const cur_seg(seg.conn_ix[car.dir]), isec_type(next_road_type - TYPE_ISEC2);
		unsigned isec_city(city_id); // default to ourself

		if (!road_to_city.empty()) { // on connector road
			assert(seg.road_ix < road_to_city.size());
			unsigned const city_ix(road_to_city[seg.road_ix].id[car.dir]);
			if (city_ix != CONN_CITY_IX) {isec_city = city_ix;} // moving into a city
		}
		road_isec_t const &isec(get_city(isec_city, road_networks, global_rn).get_isec(isec_type, cur_seg));
		if (!isec.has_stopsign && isec.can_go_now(car)) return 0; // green light, no stop sign
		// it's difficult to calculate the correct decleration to stop at the right spot, especially when framerate is variable, so directly adjust the speed instead
		// use linear deceleration to 10% of max speed (for static vs. kinetic friction); should this vary per-car?
		float const max_speed(0.1f*car.max_speed + (dist_to_end/stop_dist)*car.get_max_speed());
		assert(max_speed > 0.0);
		if (car.cur_speed < max_speed) return 0;
		car.cur_speed  = max_speed; // clamp to the max
		car.is_braking = 1;
		return 1;
	}
	void find_car_next_seg(car_t &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
		if (car.cur_road_type == TYPE_RSEG) {
			road_seg_t const &seg(get_car_seg(car));
			car.cur_road_type = seg.conn_type[car.dir];
			car.cur_road      = seg.road_ix;
			car.cur_seg       = seg.conn_ix[car.dir];

			if (!road_to_city.empty()) { // on connector road
				assert(car.cur_road < road_to_city.size());
				unsigned const city_ix(road_to_city[car.cur_road].id[car.dir]);
					
				if (car.in_isect() && city_ix != CONN_CITY_IX) { // moving into a city
					assert(city_ix < road_networks.size());
					road_network_t const &rn(road_networks[city_ix]);
					road_isec_t const &isec(rn.get_isec(car.get_isec_type(), car.cur_seg)); // must be a 3-way or 4-way intersection
					unsigned orient(0);
						
					if (car.cur_road_type == TYPE_ISEC4 && car.turn_dir == TURN_NONE) { // straight through a 4-way isec (but may not have entered isec yet, so turn_dir isn't valid)
						orient = car.get_orient();
					}
					else {
						orient = 2*(!car.dim) + 0; // use the road in the other dim, since it must be within the new city (dir doesn't matter)
					}
					car.cur_city = city_ix;
					car.cur_road = isec.rix_xy[orient];
					assert(car.cur_road < rn.roads.size());
					car.entering_city = 1; // flag so that collision detection works
				}
			}
			assert(get_car_rn(car, road_networks, global_rn).get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
			return; // always within same city, no city update
		}
		if (car.cur_road_type == TYPE_DRIVEWAY) { // is this reachable? yes, on rare occasions
			find_and_set_car_road_and_seg(car);
			return;
		}
		// else it must be in an intersection
		road_isec_t const &isec(get_car_isec(car)); // conn_ix: {-x, +x, -y, +y}
		unsigned const orient(car.get_orient());
		int conn_ix(isec.conn_ix[orient]), rix(isec.rix_xy[car.get_orient()]);
		assert(isec.conn & (1<<orient));
		isec.notify_leaving_car(car);

		if (conn_ix < 0) { // city connector road case, use global_rn
			assert(rix < 0);
			conn_ix = decode_neg_ix(conn_ix);
			rix = global_rn.get_seg(conn_ix).road_ix;
			car.cur_city = CONN_CITY_IX; // move to global road network
		} else {assert(rix >= 0);} // local road
		car.cur_road = (unsigned)rix;
		car.cur_seg  = (unsigned)conn_ix;
		car.cur_road_type = TYPE_RSEG; // always connects to a road segment
		car.entering_city = 0;
		cube_t const road_bcube(get_road_bcube_for_car(car, global_rn));

		if (!road_bcube.intersects_xy(car.bcube)) { // sanity check
			cout << "bad intersection:" << endl << car.str() << endl << "bcube: " << road_bcube.str() << endl;
			assert(0);
		}
	}
	int get_next_seg(car_t const &car) const {
		if (car.in_isect()) { // use segment at exit of current intersection
			road_isec_t const &isec(get_car_isec(car));
			return isec.conn_ix[isec.get_dest_orient_for_car_in_isec(car, 1)];
		}
		else { // use current segment
			assert(car.cur_road_type == TYPE_RSEG);
			return car.cur_seg;
		}
	}
	bool car_can_fit_in_seg(car_t const &car, road_network_t const &global_rn) const {
		if (!car.car_in_front) return 1; // no car in front, assume we can fit (optimization)
		int seg_ix(get_next_seg(car));
		road_t const &seg((seg_ix < 0) ? global_rn.get_seg(decode_neg_ix(seg_ix)) : get_seg(seg_ix)); // handle global connector road
		cube_t region(seg);
		if (car.in_isect()) {region.union_with_cube(get_car_isec(car));} // include cars in the current intersection as well
		float const req_space(car.get_sum_len_space_for_cars_in_front(region)), avail_space(seg.get_length());
		//cout << "num_in_front: " << car.count_cars_in_front(region) << ", avail_space: " << avail_space << ", req_space: " << req_space << ", fits: " << (avail_space > req_space) << endl;
		return (avail_space > req_space); // check if there's enough space in straight segment
	}
	bool car_can_go_now(car_t const &car, road_network_t const &global_rn) const {
		if (!car.in_isect()) return 1; // not at an intersection
		if (car.cur_road_type != TYPE_ISEC2 && !get_car_isec(car).can_go_now(car)) return 0; // check stoplights, stopsigns, and blocked intersections
		return car_can_fit_in_seg(car, global_rn); // check if there's space, to avoid blocking the intersection
	}
private:
	void choose_another_dir(car_t &car, rand_gen_t &rgen, road_isec_t const &isec) const {
		unsigned char const orig_turn_dir(car.turn_dir);

		if (car.turn_dir == TURN_LEFT || car.turn_dir == TURN_RIGHT) { // turning left/right at intersection
			assert(isec.num_conn > 2); // must not be a bend (can't go straight, but can't be blocked)
			if (isec.is_orient_currently_valid(car.get_orient(), (unsigned)TURN_NONE)) {car.turn_dir = (unsigned char)TURN_NONE;} // give up on the left turn and go straight instead
			else {car.turn_dir = ((car.turn_dir == TURN_LEFT) ? (unsigned char)TURN_RIGHT : (unsigned char)TURN_LEFT);} // can't go straight - then go right/left instead
		}
		else if (car.turn_dir == TURN_NONE && car.turn_val != 0.0) { // was going straight, just entered the isec, turn not yet completed
			unsigned const isec_orient(car.get_orient_in_isec()); // invert dir (incoming, not outgoing)
			if      (isec.is_orient_currently_valid(stoplight_ns::conn_right[isec_orient], TURN_RIGHT)) {car.turn_dir = TURN_RIGHT;} // go right instead
			else if (isec.is_orient_currently_valid(stoplight_ns::conn_left [isec_orient], TURN_LEFT )) {car.turn_dir = TURN_LEFT ;} // go left instead
		}
		if (car.turn_dir != orig_turn_dir && isec.is_global_conn_int()) { // change of turn dir at global connector road intersection
			if (isec.rix_xy[isec.get_dest_orient_for_car_in_isec(car, 1)] < 0) {car.turn_dir = orig_turn_dir;} // abort turn change to avoid going to the wrong city
		}
		if (car.turn_dir != orig_turn_dir) {car.on_alternate_turn_dir(rgen);}
	}
	void stop_and_wait_car(car_t &car, rand_gen_t &rgen, vector<road_network_t> const &road_networks, road_isec_t const &isec) const {
		bool const yellow_light(isec.yellow_light(car)); // this logic is only triggered on yellow lights

		if (car.rot_z == 0.0 && yellow_light) { // not in the middle of a turn; helps with gridlock at connector roads
			// Note that right turns can't be blocked by cross traffic
			if (car.turn_dir != TURN_RIGHT && !isec.stoplight.check_int_clear(car))     {choose_another_dir(car, rgen, isec);} // light turned yellow and isec still blocked
			else if (car.car_in_front && car.car_in_front->get_wait_time_secs() > 60.0) {choose_another_dir(car, rgen, isec);} // car in front has been stopped for > 60s
		}
		car.stopped_at_light = 1;
		car.stop();
		float const wait_secs(car.get_wait_time_secs());

		if (wait_secs > 60.0 && yellow_light && isec.contains_pt_xy(car.get_center())) { // car in the intersection
			cout << "car waiting for " << wait_secs << " seconds" << endl;
			car.honk_horn_if_close();
			car_t const orig_car(car);
			car = car_t(); // reset default fields
			if (add_car_to_rns(car, rgen, road_networks)) {car.model_id = orig_car.model_id; car.color_id = orig_car.color_id; return;} // relocate car somewhere else
			car = orig_car; // failed (unlikely) - restore original car state
		}
	}
	void find_and_set_car_road_and_seg(car_t &car) const {
		// find the road and segment the car is currently on; this can be slow, but should be rare
		point const center(car.get_center());

		for (auto i = segs.begin(); i != segs.end(); ++i) {
			if (!i->contains_pt_xy(center)) continue;
			car.cur_seg       = (i - segs.begin());
			car.cur_road      = i->road_ix;
			car.cur_road_type = TYPE_RSEG;
			return;
		} // for i
		cout << "Warning: Car exiting driveway not centered on a road segment: " << TXT(car.str()) << endl;

		// cars exiting a driveway by backing up and turning right may end in the intersection rather than the road segment;
		// find the road segment the car overlaps - there should be only one; technically we should only need to check the front of the car, but this is safer
		for (auto i = segs.begin(); i != segs.end(); ++i) {
			if (!i->intersects_xy(car.bcube)) continue;
			car.cur_seg       = (i - segs.begin());
			car.cur_road      = i->road_ix;
			car.cur_road_type = TYPE_RSEG;
			return;
		} // for i
		cout << "Error: Car not overlapping a road segment" << endl;
		assert(0); // should never get here
	}
	int find_road_for_car(car_t const &car, bool dim) const {
		for (auto r = roads.begin(); r != roads.end(); ++r) { // similar to get_nearby_road_ix(), except checks cube overlap rather than point containment
			if (r->dim != dim) continue; // wrong direction
			if (r->intersects_xy(car.bcube)) {return (r - roads.begin());}
		}
		return -1; // road not found - error?
	}

	bool dest_driveway_in_this_city(car_t const &car) const {
		return (car.dest_driveway >= 0 && car.dest_city == city_id);
	}
	driveway_t const &get_driveway(unsigned dix) const {
		assert(dix < city_obj_placer.driveways.size());
		return city_obj_placer.driveways[dix];
	}
	cube_t get_driveway_cube_inc_parking_lot(driveway_t const &driveway) const { // Note: currently unused
		if (!driveway.is_parking_lot()) return driveway;
		cube_t bc(driveway);
		bc.union_with_cube(city_obj_placer.parking_lots[driveway.park_lot_ix]);
		return bc;
	}

	bool run_car_in_driveway_logic(car_t &car, vector<car_t> const &cars, rand_gen_t &rgen) const {
		assert(city_params.cars_use_driveways);
		driveway_t const &driveway(get_driveway(car.cur_seg));
		// if pedestrian(s) in driveway, stop and wait unless we've been waiting for > 60s
		if (driveway.has_recent_ped() && car.get_wait_time_secs() < 60.0) {car.sleep(rgen, 0.5); return 1;}

		if (car.dest_driveway == (int)car.cur_seg) { // entering driveway
			car.pull_into_driveway(driveway, rgen);
			return 1; // no other logic to run here
		}
		// else leaving driveway
		bool const ddim(driveway.dim);
		bool in_driveway(driveway.contains_cube_xy(car.bcube));
		// for cars in parking lots, we only need to check the bounds of the car inside the driveway dim, since parking spaces are to the side
		in_driveway |= (driveway.is_parking_lot() && car.bcube.d[ddim][0] >= driveway.d[ddim][0] && car.bcube.d[ddim][1] <= driveway.d[ddim][1]);
		
		if (in_driveway) { // car still in driveway or parking lot, continue to pull/back out
			car.back_or_pull_out_of_driveway(driveway);
			return 1;
		}
		int const road_ix(find_road_for_car(car, !ddim));

		if (road_ix < 0) {
			cerr << car.str() << TXT(driveway.str()) << endl;
			assert(0);
		}
		float const road_center(roads[road_ix].get_center_dim(car.dim));
		float const centerline(road_center + (driveway.dir ? -1.0 : 1.0)*get_car_lane_offset());
		if (!car.exit_driveway_to_road(cars, driveway, centerline, road_ix, rgen)) return 1; // still exiting driveway
		find_and_set_car_road_and_seg(car);
		return 0; // done, continue with car logic in update_car() below
	}
public:
	void update_car(car_t &car, vector<car_t> const &cars, rand_gen_t &rgen, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
		assert(car.cur_city == city_id);
		if (car.is_parked()) return; // stopped, no update (for now)
		car.is_braking = 0; // reset for this frame
			
		if (car.cur_road_type == TYPE_DRIVEWAY) { // moving in a driveway; could also be TYPE_PARK_LOT
			if (run_car_in_driveway_logic(car, cars, rgen)) return;
		}
		if (dest_driveway_in_this_city(car) && !car.in_isect() && !car.stopped_at_light) { // turning into an intersection is not the same as a driveway
			if (car.run_enter_driveway_logic(cars, get_driveway(car.dest_driveway))) return;
		}
		if (car.in_isect()) {
			road_isec_t const &isec(get_car_isec(car));
			isec.notify_waiting_car(car); // even if not stopped

			// if we started to pull into the intersection but aren't at least halfway inside, and we can't fit ourselves in the dest road,
			// then wait at the light rather than blocking the intersection
			if (isec.contains_pt_xy(car.get_front()) && !isec.contains_pt_xy(car.get_center()) && !car_can_fit_in_seg(car, global_rn)) {
				bool const was_stopped(car.stopped_at_light);
				stop_and_wait_car(car, rgen, road_networks, isec);
				if (!was_stopped) {isec.notify_leaving_car(car);} // unset the entering flag so that we don't keep ourselves from entering again (and deadlocking)
				return;
			}
		}
		if (car.stopped_at_light) {
			bool const was_stopped(car.is_stopped());
			if (car_can_go_now(car, global_rn)) {car.stopped_at_light = 0;} // can go now
			else if (car.in_isect()) {stop_and_wait_car(car, rgen, road_networks, get_car_isec(car));} // Note: is_isect test allows cars to coast through lights when decel is very low
			if (was_stopped) {car.stopped_for_ssign = 1;} // maybe at stop sign, or maybe at stoplight; after checking to ensure a full stop
			if (was_stopped) return; // no update needed
		}
		else if (check_future_isec_stop(car, road_networks, global_rn)) {} // no update
		else {car.maybe_accelerate();}

		cube_t const road_bcube(get_road_bcube_for_car(car));
		if (!bcube.intersects_xy(car.prev_bcube)) {cout << car.str() << endl << road_bcube.str() << endl; assert(0);} // sanity check
		bool const dim(car.dim);
		car.in_tunnel = point_in_tunnel(car.get_center());

		if (car.cur_road_type == TYPE_RSEG && road_bcube.z1() != road_bcube.z2()) { // car on connector road
			assert(car.cur_road_type == TYPE_RSEG);
			assert(car.cur_city == CONN_CITY_IX);
			bool const slope(get_car_seg(car).slope);
			float const car_pos(car.bcube.get_center_dim(dim)); // center of car in dim
			float const road_len(road_bcube.get_sz_dim(dim));
			assert(road_len > TOLERANCE);
			float const t((car_pos - road_bcube.d[dim][0])/road_len); // car pos along road in (0.0, 1.0)
			float const road_z(t*road_bcube.d[2][!slope] + (1.0 - t)*road_bcube.d[2][slope]);
			float const car_len(car.get_length());
			car.dz = ((slope ^ car.dir) ? 1.0 : -1.0)*road_bcube.dz()*(car_len/road_len);
			car.bcube.z1() = road_z - 0.5*fabs(car.dz);
			car.bcube.z2() = road_z + 0.5*fabs(car.dz) + car.height;
		}
		else if (car.dz != 0.0) { // car moving from connector road to level city
			float const road_z(road_bcube.z1());
			car.dz = 0.0;
			set_cube_zvals(car.bcube, road_z, (road_z + car.height));
		}
		if (city_params.cars_use_driveways && car.turn_dir != TURN_NONE && !car.in_isect()) {
			assert(dest_driveway_in_this_city(car));
			cout << car.str() << TXT(car.prev_bcube.str()) << TXT(get_driveway(car.dest_driveway).str()) << endl;
			car.turn_dir = TURN_NONE; // hack to handle misbehaving cars turning into driveways (maybe missed the turn because it was blocked?)
		}
		if (car.turn_dir != TURN_NONE) {
			assert(car.in_isect());
			float const isec_center(road_bcube.get_cube_center()[dim]);
			float const centerline(isec_center + (((car.turn_dir == TURN_RIGHT) ^ car.dir) ? 1.0 : -1.0)*get_car_lane_offset());
			car_t const car_pre_turn(car); // make a copy before turning that we can pass into isec.notify_leaving_car()
			car.maybe_apply_turn(centerline, 0); // for_driveway=0

			if (car.turn_dir == TURN_NONE) { // turn has been completed
				road_isec_t const &isec(get_car_isec(car));
				// special case for stop signs to handle orient change
				isec.notify_leaving_car(car_pre_turn); // leaving this orient, entering the exit orient
				isec.notify_turned_car (car); // register in new orient so that waiting cars don't see this slot as open for one frame
				
				if (isec.conn_ix[car.get_orient()] >= 0) {
					short const rix(isec.rix_xy[car.get_orient()]);
					assert(rix >= 0); // not connector road
					car.cur_road = rix; // switch to using road_ix in new dim
				}
			}
		}
		if (road_bcube.contains_cube_xy(car.bcube)) { // in same road seg/int
			assert(get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
			return; // done
		}
		point const car_front(car.get_front(0.375)); // near the front, so that we can stop at the intersection

		// car crossing the border of this bcube, update state
		if (!road_bcube.contains_pt_xy_inc_low_edge(car_front)) { // move to another road seg/int
			find_car_next_seg(car, road_networks, global_rn);
				
			if (car.in_isect()) { // moved into an intersection, choose direction
				road_network_t const &car_rn(get_car_rn(car, road_networks, global_rn));
				road_isec_t const &isec(car_rn.get_car_isec(car)); // Note: either == city_id, or just moved from global to a new city
				unsigned const orient_in(car.get_orient_in_isec()); // invert dir (incoming, not outgoing)
				assert(isec.conn & (1<<orient_in)); // car must come from an enabled orient
				unsigned orients[3] = {}; // {straight, left, right}
				orients[TURN_NONE ] = car.get_orient(); // straight
				orients[TURN_LEFT ] = stoplight_ns::conn_left [orient_in];
				orients[TURN_RIGHT] = stoplight_ns::conn_right[orient_in];

				// use dest_seg.car_count to estimate traffic and route around?
				if (car.dest_valid && car.cur_city != CONN_CITY_IX) { // Note: don't need to update dest logic on connector roads since there are no choices to make
					vector3d dest_dir;
						
					if (is_car_at_dest_isec(car)) { // this intersection is our destination
						if (dest_driveway_in_this_city(car)) { // drive toward the dest driveway
							driveway_t const &driveway(get_driveway(car.dest_driveway));
							point dest_pos(driveway.get_cube_center());
							dest_pos[driveway.dim] = driveway.get_edge_at_road();
							dest_dir = dest_pos - isec.get_cube_center();
						}
						else { // prefer to go straight, otherwise random
							dest_dir[ dim] = (car.dir ? 1.0 : -1.0); // straight
							dest_dir[!dim] = 0.1*(rgen.rand_bool() ? 1.0 : -1.0);
						}
					}
					else { // drive toward the destination intersection
						point const dest_pos(car_rn.get_car_dest_isec_center(car, road_networks, global_rn));
						dest_dir = dest_pos - car.get_center();
					}
					dest_dir.z = 0.0; // always level
					bool const pri_dim(fabs(dest_dir.x) < fabs(dest_dir.y)), pri_dir(dest_dir[pri_dim] > 0), sec_dir(dest_dir[!pri_dim] > 0);
					unsigned best_score(0);

					for (unsigned tdir = 0; tdir < 3; ++tdir) { // choose best scoring of all valid turn dirs from {none/straight, left, right}
						unsigned const orient(orients[tdir]);
						if (!isec.is_orient_currently_valid(orient, tdir)) continue; // can't turn in this dir

						if (isec.conn_to_city >= 0 && isec.conn_ix[orient] < 0) { // city connector isec
							if (isec.conn_to_city != car.dest_city) continue; // leads to incorrect city, skip
							car.turn_dir = tdir; // this is our destination - done
							best_score = 1; // set to avoid assertion failure below
							break;
						}
						bool const dim2((orient >> 1) != 0), dir2(orient & 1);
						unsigned score(1); // start at lowest valid score
						if      (dim2 == pri_dim && dir2 == pri_dir) {score = 3;} // best score
						else if (dim2 != pri_dim && dir2 == sec_dir) {score = 2;} // second best score
						if (score > best_score) {best_score = score; car.turn_dir = tdir;}
					} // for tdir
					assert(best_score > 0); // no dead end roads
				}
				else { // use random turn direction
					while (1) {
						unsigned new_turn_dir(0); // force turn on global conn road 75% of the time to get more cars traveling between cities
						bool const force_turn(isec.is_global_conn_int() && (rgen.rand()&3) != 0);
						int const rval(rgen.rand()%(force_turn ? 2 : 4));
						if      (rval == 0) {new_turn_dir = TURN_LEFT ;} // 25%
						else if (rval == 1) {new_turn_dir = TURN_RIGHT;} // 25%
						else                {new_turn_dir = TURN_NONE ;} // 50%
						if (new_turn_dir == car.front_car_turn_dir && (rgen.rand()%4) != 0) continue; // car in front is too slow, don't turn the same way as it
						if (isec.is_orient_currently_valid(orients[new_turn_dir], new_turn_dir)) {car.turn_dir = new_turn_dir; break;} // success
					} // end while
				}
				assert(isec.conn & (1<<orients[car.turn_dir]));
				car.front_car_turn_dir = TURN_UNSPEC; // reset state now that it's been used
				car.stopped_for_ssign  = 0; // reset for this intersection
				car.stopped_at_light   = (isec.red_or_yellow_light(car) || !car_rn.car_can_go_now(car, global_rn)); // is the red_or_yellow_light() call redundant?
				if (car.stopped_at_light) {car.stop();}
				if (car.turn_dir != TURN_NONE) {car.begin_turn();} // capture car centerline before the turn
			}
		} // end move to another road segment
		assert(get_car_rn(car, road_networks, global_rn).get_road_bcube_for_car(car, global_rn).intersects_xy(car.bcube)); // sanity check
	}
	bool is_car_at_dest_isec(car_t const &car) const {
		unsigned const isec_type(car.get_isec_type());
		unsigned flat_isec_ix(car.cur_seg);
		for (unsigned n = 0; n < isec_type; ++n) {flat_isec_ix += isecs[n].size();}
		return (car.dest_isec == flat_isec_ix); // dest_isec is in flat space, while the current isec is defined by {cur_road_type, cur_seg)
	}
	road_isec_t const &get_car_dest_isec(car_t const &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
		if (car.dest_city == city_id) {return get_isec_by_ix(car.dest_isec);} // local destination within the current city
		assert(car.dest_city < road_networks.size());
		road_isec_t const *const isec(find_isec_to_dest_city(car, road_networks[car.dest_city], global_rn)); // destination in another city
		assert(isec != nullptr); // path must exist, otherwise this city wouldn't have been chosen
		return *isec;
	}
	point get_car_dest_isec_center(car_t const &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
		return get_car_dest_isec(car, road_networks, global_rn).get_cube_center();
	}
	cube_t const &get_car_dest(car_t const &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn, bool isec_only) const { // isec or driveway
		if (!isec_only && dest_driveway_in_this_city(car)) {return get_driveway(car.dest_driveway);} // driveway in this city
		return get_car_dest_isec(car, road_networks, global_rn); // else assume we want andest intersection
	}
	int get_parking_lot_ix_for_car(car_t const &car) const {
		if (!car.in_parking_lot)    return -1; // caller should really check this
		if (car.dest_driveway >= 0) return get_driveway(car.dest_driveway).park_lot_ix; // have a cached driveway (optimization)
		// dest_driveway is reset after the car is parked, so we must re-query; this could be optimized by caching the parking lot index in the car
		return get_parking_lot_ix(car.get_center(), 1); // inc_driveways=1
	}
private:
	road_isec_t const *find_isec_to_dest_city(car_t const &car, road_network_t const &dest_rn, road_network_t const &global_rn) const {
		assert(car. cur_city == city_id);
		assert(car.dest_city == dest_rn.city_id);
		assert(dest_rn.city_id != city_id); // not ourself
		// Note: here we don't attempt to find shortcuts through other cities as this would be quite complex
		auto it(cix_to_isec.find(car.dest_city));
		if (it != cix_to_isec.end()) {return it->second;} // found
		return nullptr; // not found, caller can error check
	}
	bool select_avail_driveway_or_parking_space(car_t &car, rand_gen_t &rgen) const {
		// consider destination driveways, possibly including parking lot entrances
		if (city_obj_placer.driveways.empty()) return 0; // no driveways
		float const car_len(car.get_length()), car_width(car.get_width());
		
		for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to find a valid driveway
			unsigned const dix(rgen.rand()%city_obj_placer.driveways.size());
			driveway_t const &driveway(get_driveway(dix));
			if (driveway.in_use) continue;
			if (driveway.get_length() < 1.5*car_len  ) continue; // driveway is too short
			if (driveway.get_width () < 1.1*car_width) continue; // driveway is too narrow (mostly applies to trucks)

			if (!is_residential) { // not a residential city
				bool const allow_hcap(rgen.rand_float() < 0.25);
				int const psix(city_obj_placer.select_dest_parking_space(dix, allow_hcap, 1, car_len, rgen)); // reserve_spot=1
				if (psix < 0) continue;
				point const ps_center(city_obj_placer.get_parking_space_center(psix));
				car.park_space_cent.assign(ps_center.x, ps_center.y);

				if (driveway.park_lot_ix >= 0) { // should always be true
					// clamp point inside the parking lot so that large vehicles such as ambulances don't stick out; parking space should be large enough
					float const ends_pad(1.01*0.5*car_len);
					float &val(car.park_space_cent[!driveway.dim]);
					parking_lot_t const &parking_lot(city_obj_placer.parking_lots[driveway.park_lot_ix]);
					max_eq(val, parking_lot.d[!driveway.dim][0]+ends_pad);
					min_eq(val, parking_lot.d[!driveway.dim][1]-ends_pad);
				}
			}
			driveway.in_use   = 1; // temporarily in use
			car.dest_driveway = (unsigned short)dix;
			// find intersection before the driveway such that driving on the road exiting this intersection will encounter the driveway on the right
			bool const dim(driveway.dim), dir(driveway.dir), extend_dir(dim ^ dir); // extend_dir is the direction of the last intersection before our right turn
			float const dw_road_meet(driveway.d[dim][dir]); // point at which the road and driveway meet
			assert(!plots.empty());
			float const road_spacing(plots.front().get_sz_dim(!dim)); // use actual segment length (should be the same across segments)
			cube_t query_cube(driveway.extend_across_road()); // segment of road connected to driveway
			query_cube.d[ dim][!dir] = dw_road_meet;
			query_cube.d[!dim][extend_dir] += (extend_dir ? 1.0 : -1.0)*road_spacing; // extend out to include the adjacent intersection in this dim/dir, assuming driveway on a road seg
			car.dest_isec = 0; // used as a loop index

			// this is slow, but we need the correct index; should be rarely called
			for (unsigned n = 0; n < 3; ++n) { // {2-way, 3-way, 4-way}
				for (auto i = isecs[n].begin(); i != isecs[n].end(); ++i, ++car.dest_isec) {
					if (i->intersects_xy(query_cube)) return 1; // car.dest_isec is the correct value
				}
			}
			cout << TXT(dim) << TXT(dir) << TXT(extend_dir) << TXT(road_spacing) << TXT(dw_road_meet) << TXT(driveway.str()) << TXT(query_cube.str()) << endl;
			assert(0); // should never get here
		} // for n
		return 0; // failed
	}
public:
	bool choose_new_car_dest(car_t &car, rand_gen_t &rgen) const {
		// select a driveway if one is available and we're in the dest city; otherwise, select an intersection
		//assert(car.dest_driveway < 0); // generally okay, but could maybe fail due to floating-point error? better to reset below?
		car.dest_driveway = -1; // reset; if nonzero, that may mean this driveway is never used after this point
		if (city_params.cars_use_driveways && car.cur_city == car.dest_city && select_avail_driveway_or_parking_space(car, rgen)) return 1;
		unsigned const num_tot(isecs[0].size() + isecs[1].size() + isecs[2].size()); // include 2-way, 3-way, and 4-way intersections
		if (num_tot == 0) return 0; // no isecs to select
		car.dest_isec = (unsigned short)(rgen.rand() % num_tot);
		return 1;
	}
	bool car_at_dest(car_t const &car) const {
		if (car.cur_road_type == TYPE_DRIVEWAY) return 0; // driveway desination check is handled in update_car()
		if (!dest_driveway_in_this_city(car)) {return get_isec_by_ix(car.dest_isec).contains_pt_xy(car.get_center());} // dest isec
		cube_t const &driveway(get_driveway(car.dest_driveway));
		cube_t car_bcube(car.bcube);
		car_bcube.d[car.dim][!car.dir] -= (car.dir ? 1.0 : -1.0)*0.2*car.get_length(); // add some space behind the car
		return driveway.contains_cube_xy(car_bcube);
	}
	road_isec_t const &get_isec_by_ix(unsigned ix) const {
		for (unsigned n = 0; n < 3; ++n) {
			unsigned const sz(isecs[n].size());
			if (ix < sz) {return isecs[n][ix];}
			ix -= sz;
		}
		assert(0); // should never get here (invalid ix)
		return isecs[0][0]; // never gets here
	}
	void next_frame() {
		for (unsigned n = 1; n < 3; ++n) { // {2-way, 3-way, 4-way} - Note: 2-way can be skipped
			for (auto i = isecs[n].begin(); i != isecs[n].end(); ++i) {i->next_frame();} // update stoplight state
		}
		for (auto i = segs.begin(); i != segs.end(); ++i) {i->next_frame();}
		//cout << TXT(city_id) << TXT(tot_road_len) << TXT(num_cars) << TXT(get_traffic_density()) << endl;
		num_cars = 0;
		city_obj_placer.next_frame();
		city_obj_placer.play_sounds();
	}
	static road_network_t const &get_city(unsigned city_ix, vector<road_network_t> const &road_networks, road_network_t const &global_rn) {
		if (city_ix == CONN_CITY_IX) return global_rn;
		assert(city_ix < road_networks.size());
		return road_networks[city_ix];
	}
	static road_network_t const &get_car_rn(car_base_t const &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) {
		return get_city(car.cur_city, road_networks, global_rn);
	}
	void update_car_seg_stats(car_base_t const &car) const {
		if (car.cur_road_type == TYPE_RSEG) {++get_car_seg(car).car_count;}
	}
	road_seg_t const &get_car_seg(car_base_t const &car) const {
		assert(car.cur_road_type == TYPE_RSEG);
		return get_seg(car.cur_seg);
	}
	road_isec_t const &get_car_isec(car_base_t const &car) const {
		auto const &iv(isecs[car.get_isec_type()]);
		assert(car.cur_seg < iv.size());
		return iv[car.cur_seg];
	}
	cube_t get_road_bcube_for_car(car_base_t const &car) const {
		if (car.cur_road_type == TYPE_DRIVEWAY) {return get_driveway(car.cur_seg);} // is this case used?
		assert(car.cur_road < roads.size()); // Note: generally holds, but not required
		if (car.cur_road_type == TYPE_RSEG) {return get_car_seg(car);}
		return get_car_isec(car);
	}
	cube_t get_road_bcube_for_car(car_base_t const &car, road_network_t const &global_rn) const { // function variant that works with both global and local roads
		if (car.cur_city == city_id) {return get_road_bcube_for_car(car);}
		else if (car.cur_city == CONN_CITY_IX) {return global_rn.get_road_bcube_for_car(car);}
		else {assert(0); return cube_t();}
	}
	road_isec_t const *find_isec_containing_pt(point const &pos, unsigned ns, unsigned ne) const {
		assert(ns < ne && ne <= 3);
		for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
			if (!b->bcube.contains_pt_xy(pos)) continue;
			for (unsigned n = ns; n < ne; ++n) { // {2-way, 3-way, 4-way}
				range_pair_t const &range(b->ranges[TYPE_ISEC2 + n]);
				for (auto i = isecs[n].begin()+range.s; i != isecs[n].begin()+range.e; ++i) {
					if (i->contains_pt_xy(pos)) {return &(*i);}
				}
			}
		} // for b
		return nullptr;
	}
	dw_query_t get_nearby_driveway(unsigned global_plot_id, point const &pos, float dist) const {
		if (city_obj_placer.driveways.empty()) return dw_query_t(); // no driveways
		//if (!is_residential) return dw_query_t(); // logic doesn't really work well for city parking lot driveways as cars can get stuck, so skip
		unsigned const plot_id(decode_plot_id(global_plot_id));

		// plot_ix isn't required, but it's likely faster to check the plot first in the iteration
		for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) { // iterate by tile, which is likely faster than iterating over driveways
			if (!b->bcube.contains_pt_xy_exp(pos, dist)) continue;
			range_pair_t const &rp(b->ranges[TYPE_DRIVEWAY]);
				
			for (unsigned i = rp.s; i < rp.e; ++i) {
				driveway_t const &driveway(city_obj_placer.driveways[i]);
				if (driveway.plot_ix != plot_id) continue; // wrong plot (optimization)
				if (driveway.contains_pt_xy_exp(pos, dist)) return dw_query_t(&driveway, i); // found
			}
		} // for b
		return dw_query_t(); // driveway not found
	}
	bool mark_crosswalk_in_use(point const &pos, bool dim, bool dir) const {
		road_isec_t const *isec(find_isec_containing_pt(pos, 1, 3)); // 2-way can be skipped because there's no light/crosswalk
		if (isec == nullptr) return 0; // ped not at a crosswalk, maybe crossing in the middle of the street; this is okay for now, nothing else to do in this case
		isec->mark_crosswalk_in_use(dim, dir);
		return 1;
	}
	bool check_isec_sphere_coll(point const &pos, float radius, cube_t &coll_cube) const {
		road_isec_t const *isec(find_isec_containing_pt(pos, 1, 3)); // 2-way can be skipped because there's no light/crosswalk
		if (isec == nullptr) return 0;
		return isec->check_sphere_coll(pos, radius, coll_cube);
	}
	int get_nearby_road_ix(point const &pos, bool road_dim) const {
		for (auto r = roads.begin(); r != roads.end(); ++r) {
			if (r->dim != road_dim) continue;
			cube_t road_bcube(*r);
			road_bcube.expand_by_xy(0.01*city_params.road_width); // expand slightly to pick up roads that are adjacent to this point
			if (road_bcube.contains_pt_xy(pos)) {return (r - roads.begin());}
		}
		return -1; // should never get here, but occasionally can due to bad collisions between peds, floating-point error, etc.
	}
	int get_parking_lot_ix(point const &pos, bool inc_driveways) const {
		if (city_obj_placer.parking_lots.empty()) return -1;

		// iterate by tile, which is faster than iterating over parking lots and driveways; should we store maintain and use a plot_ix => parking_lot_ix mapping?
		for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
			if (!b->bcube.contains_pt_xy(pos)) continue;
			range_pair_t const &rpp(b->ranges[TYPE_PARK_LOT]);
			
			for (unsigned i = rpp.s; i < rpp.e; ++i) { // check parking lots
				if (city_obj_placer.parking_lots[i].contains_pt_xy(pos)) return i;
			}
			if (!inc_driveways) continue;
			range_pair_t const &rpd(b->ranges[TYPE_DRIVEWAY]);

			for (unsigned i = rpd.s; i < rpd.e; ++i) { // check driveways
				driveway_t const &d(city_obj_placer.driveways[i]);
				if (d.park_lot_ix >= 0 && d.contains_pt_xy(pos)) return d.park_lot_ix;
			}
		} // for b
		return -1; // not found
	}
}; // road_network_t

class city_road_gen_t : public road_gen_base_t {
	vector<road_network_t> road_networks; // one per city
	road_network_t global_rn; // connects cities together; no plots
	vector<transmission_line_t> transmission_lines;
	cube_t cities_bcube;
	road_draw_state_t dstate;
	rand_gen_t rgen;
	bool have_plot_dividers;

	static float rgen_uniform(float val1, float val2, rand_gen_t &rgen) {return (val1 + (val2 - val1)*rgen.rand_float());}

	void assign_city_clusters() {
		vector<unsigned char> used(road_networks.size(), 0);
		vector<unsigned> pend;
		unsigned cluster_id(0);

		for (unsigned i = 0; i < used.size(); ++i) {
			if (used[i]) continue; // already used
			pend.push_back(i);

			while (!pend.empty()) {
				unsigned const id(pend.back());
				pend.pop_back();
				road_networks[id].set_cluster(cluster_id);
				used[id] = 1;
				set<unsigned> const &conn(road_networks[id].get_connected());

				for (auto c = conn.begin(); c != conn.end(); ++c) {
					assert(*c < used.size());
					if (!used[*c]) {pend.push_back(*c);}
				}
			} // while
			++cluster_id;
		} // for i
		cout << "City clusters: " << cluster_id << endl;
	}
	road_network_t       &get_city_by_ix(unsigned ix)       {assert(ix < road_networks.size()); return road_networks[ix];}
	road_network_t const &get_city_by_ix(unsigned ix) const {assert(ix < road_networks.size()); return road_networks[ix];}

public:
	city_road_gen_t() : have_plot_dividers(0) {}
	bool empty() const {return road_networks.empty();}
	bool has_tunnels() const {return global_rn.has_tunnels();}
	bool point_in_tunnel(point const &pos) const {return global_rn.point_in_tunnel(pos);}
	road_network_t const &get_city(unsigned city_ix) const {return road_network_t::get_city(city_ix, road_networks, global_rn);} // call the static function
	cube_t const &get_city_bcube(unsigned city_ix) const {return get_city(city_ix).get_bcube();}
	cube_t const &get_city_plot_bcube(unsigned city_ix, unsigned plot_ix) const {return get_city(city_ix).get_plot_bcube(plot_ix);}
	vect_cube_t const &get_colliders_for_plot(unsigned city_ix, unsigned global_plot_id) const {return get_city(city_ix).get_colliders_for_plot(global_plot_id);}
	cube_t const &get_car_dest_bcube(car_t const &car, bool isec_only) const {return get_city(car.dest_city).get_car_dest(car, road_networks, global_rn, isec_only);}
	int get_parking_lot_ix_for_car(car_t const &car) const {return get_city(car.dest_city).get_parking_lot_ix_for_car(car);}

	bool cube_int_underground_obj(cube_t const &c) const {
		if (global_rn.cube_intersect_tunnel(c)) return 1;

		for (road_network_t const &rn : road_networks) {
			if (rn.cube_int_underground_obj(c)) return 1;
		}
		return 0;
	}
	bool is_invalid_placement_for_cube(cube_t const &c) const {
		if (global_rn.is_invalid_placement_for_cube(c)) return 1;

		for (road_network_t const &rn : road_networks) {
			if (rn.is_invalid_placement_for_cube(c)) return 1;
		}
		return 0;
	}
	void add_plot_cut(cube_t const &cut) {
		for (road_network_t &rn : road_networks) {rn.add_plot_cut(cut);}
	}
	void add_city_ug_elevator_entrance(ug_elev_info_t const &uge) {
		for (road_network_t &rn : road_networks) {rn.add_city_ug_elevator_entrance(uge);}
	}
	void get_ponds_in_xy_range(cube_t const &range, vect_cube_t &ponds) const {
		for (road_network_t const &rn : road_networks) {rn.get_ponds_in_xy_range(range, ponds);}
	}
	cube_t get_city_bcube_for_cars(unsigned city_ix) const {
		cube_t bcube(get_city_bcube(city_ix));
		bcube.expand_by_xy(city_params.get_max_car_size().x); // expand by car length to fully include cars that are partially inside connector road intersections
		return bcube;
	}
	cube_t get_city_bcube_at_pt(point const &pos) const { // skips global_rn
		for (road_network_t const &rn : road_networks) {
			if (rn.get_bcube().contains_pt_xy(pos)) return rn.get_bcube();
		}
		return cube_t(); // not found
	}
	cube_t get_city_bcube_overlapping(cube_t const &c) const { // skips global_rn
		for (road_network_t const &rn : road_networks) {
			if (rn.get_bcube().intersects_xy(c)) return rn.get_bcube(); // return the first overlapping cube; assumes c is small and doesn't overlap multiple cities
		}
		return cube_t(); // not found
	}
	void add_manhole(point const &pos, float radius) {
		for (road_network_t &rn : road_networks) {
			if (rn.get_bcube().contains_pt_xy(pos)) {rn.add_manhole(pos, radius); break;}
		}
	}
	bool cube_overlaps_road_xy    (cube_t const &c, unsigned city_ix) const {return get_city(city_ix).cube_overlaps_road_xy    (c);}
	bool cube_overlaps_pl_or_dw_xy(cube_t const &c, unsigned city_ix) const {return get_city(city_ix).cube_overlaps_pl_or_dw_xy(c);}

	void gen_roads_and_plots(cube_t const &region, float road_width, vector2d const &road_spacing, bool is_residential) {
		//timer_t timer("Gen Roads"); // ~0.5ms
		road_networks.push_back(road_network_t(region, road_networks.size(), is_residential));
		if (!road_networks.back().gen_road_grid(road_width, road_spacing)) {road_networks.pop_back(); return;}
		if (cities_bcube.is_all_zeros()) {cities_bcube = region;} else {cities_bcube.union_with_cube(region);} // Note: region is zero height
	}

	void get_all_conn_road_endpoints(road_network_t const &rn, cube_t const &dest_bcube, vector<road_endpoint_t> &rpts) {
		for (auto r = rn.get_roads().begin(); r != rn.get_roads().end(); ++r) {
			bool const dir(city_road_connector_t::get_closer_dir(*r, dest_bcube, r->dim));
			if (!rn.check_valid_conn_intersection(*r, r->dim, dir, 1)) continue; // not a valid 3-way intersection; is_4_way1=1
			point pt;
			pt[ r->dim] = r->d[r->dim][dir]; // extend to the outside of the intersection; should this be rn.get_bcube().d[r->dim][dir]?
			pt[!r->dim] = r->get_center_dim(!r->dim);
			pt.z        = r->z2();
			rpts.emplace_back(pt, r->dim, dir);
		} // for r
	}

	float connect_two_cities_new(unsigned city1, unsigned city2, vect_cube_t &blockers, city_road_connector_t &crc, float road_width) {
		// Note: ignores the value of city_params.make_4_way_ints because this will always start and end at a 4-way intersection
		road_network_t &rn1(get_city_by_ix(city1)), &rn2(get_city_by_ix(city2));
		cube_t const &bcube1(rn1.get_bcube()), &bcube2(rn2.get_bcube());
		heightmap_query_t &hq(crc.hq);
		float const road_hwidth(0.5*road_width);
		road_cand_t cand, best_cand;
		vect_cube_t active_blockers;
		
		// determine the set of active blockers, which exclude the two cities to be connected
		for (auto b = blockers.begin(); b != blockers.end(); ++b) {
			if (!b->contains_cube(bcube1) && !b->contains_cube(bcube2)) {active_blockers.push_back(*b);}
		}
		// try to extend all permutations of roads on the shared sides between cities
		vector<road_endpoint_t> rpts1, rpts2;
		get_all_conn_road_endpoints(rn1, bcube2, rpts1);
		get_all_conn_road_endpoints(rn2, bcube1, rpts2);

		for (auto r1 = rpts1.begin(); r1 != rpts1.end(); ++r1) {
			for (auto r2 = rpts2.begin(); r2 != rpts2.end(); ++r2) {
				cand.clear();
				cand.start_dim = r1->dim;
				cand.cost = crc.find_route_between_points(r1->pt, r2->pt, active_blockers, cand.pts, bcube1, bcube2, road_hwidth, r1->dim, r1->dir, r2->dim, r2->dir);
				if (cand.cost > 0.0 && (best_cand.cost == 0.0 || cand.cost < best_cand.cost)) {best_cand = cand;} // update best_can if valid and a lower cost
			} // for r2
		} // for r1
		if (!best_cand.valid()) return 0.0; // failed

		// create the road segments
		assert(best_cand.pts.size() > 1); // must be at least 1 segment/2 points (in practical cases there must be at least 2 segments/3 points)
		bool fdim(best_cand.start_dim), last_dir(0);
		cube_t cur_bcube(bcube1); // start at the first city
		vector<conn_isec_t> conn_isecs;
		unsigned last_road_ix(0);
		float tot_cost(0.0);
		
		for (auto p = best_cand.pts.begin(); p+1 != best_cand.pts.end(); ++p, fdim ^= 1) {
			bool const is_first(p == best_cand.pts.begin()), is_last(p+2 == best_cand.pts.end());
			point const &p1(*p), &p2(*(p+1));
			//cout << (p-best_cand.pts.begin()) << " " << best_cand.pts.size() << " " << TXT(fdim) << TXT(p1.str()) << TXT(p2.str()) << TXT(is_first) << TXT(is_last) << endl;
			assert(!bcube1.contains_pt(p2));
			assert(!bcube2.contains_pt(p1));
			assert(p1[!fdim] == p2[!fdim]); // must be a straight road in this dim
			assert(p1[ fdim] != p2[ fdim]); // must be a nonzero length road in this dim
			bool const dir(p1[fdim] < p2[fdim]);
			cube_t next_bcube;
			if (is_last) {next_bcube = bcube2;} // end at the second city
			else {next_bcube.set_from_point(p2); next_bcube.expand_by_xy(road_hwidth);} // intersection of this road and the next one
			hq.flatten_region_to(next_bcube, city_params.road_border); // do this first to improve flattening
			unsigned const road_ix(global_rn.num_roads());
			float const cost(global_rn.create_connector_road(cur_bcube, next_bcube, blockers,
				(is_first ? &rn1 : nullptr), (is_last ? &rn2 : nullptr), (is_first ? city1 : CONN_CITY_IX), (is_last ? city2 : CONN_CITY_IX),
				city1, city2, crc, road_width, p1[!fdim], fdim, 0, is_first, is_last)); // check_only=0, is_4way=at_city
			assert(cost >= 0.0);
			tot_cost += cost;
			// if not the last segment, cache state so that the loop below can create the intersection; can't create it here without the road_ix of the next segment
			if (!is_last) {conn_isecs.emplace_back(hq, next_bcube, road_ix, fdim, dir);}
			last_road_ix = road_ix;
			last_dir     = dir;
			cur_bcube    = next_bcube;
		} // for p
		for (auto i = conn_isecs.begin(); i != conn_isecs.end(); ++i) {
			bool const is_last(i+1 == conn_isecs.end()), next_dir(is_last ? last_dir : (i+1)->dir);
			unsigned const next_road_ix(is_last ? last_road_ix : (i+1)->road_ix);
			if (i->dim) {global_rn.create_connector_bend(i->int_bcube, (next_dir ^ i->dim), (i->dir ^ i->dim), next_road_ix, i->road_ix);} // Y
			else        {global_rn.create_connector_bend(i->int_bcube, (i->dir ^ i->dim), (next_dir ^ i->dim), i->road_ix, next_road_ix);} // X
			// decrease_only=1; remove any dirt that the prev road added
			hq.flatten_sloped_region(i->fop.x1, i->fop.y1, i->fop.x2, i->fop.y2, i->fop.z1, i->fop.z2, i->fop.dim, i->fop.border, i->fop.skip_six, i->fop.skip_eix, 0, 1);
			hq.flatten_region_to(i->int_bcube, city_params.road_border, 1); // one more pass to fix mesh that was raised above the intersection by a sloped road segment
		} // for i
		return tot_cost;
	}

	float connect_two_cities(unsigned city1, unsigned city2, vect_cube_t &blockers, city_road_connector_t &crc, float road_width) {
		assert(city1 != city2); // check for self reference

		if (city_params.new_city_conn_road_alg) { // run the new algorithm first; if that fails, run the old algorithm
			float const cost(connect_two_cities_new(city1, city2, blockers, crc, road_width));
			if (cost > 0.0) return cost;
		}
		road_network_t &rn1(get_city_by_ix(city1)), &rn2(get_city_by_ix(city2));
		cube_t const &bcube1(rn1.get_bcube()), &bcube2(rn2.get_bcube());
		assert(!bcube1.intersects_xy(bcube2));
		float const min_edge_dist(4.0*road_width), min_jog(2.0*road_width);
		// Note: cost function should include road length, number of jogs, total elevation change, and max slope

		for (unsigned d = 0; d < 2; ++d) { // try for single segment
			if (city_params.make_4_way_ints > 2) continue; // only allow connector roads that have 4-way intersections at both ends (single jog case below)
			float const shared_min(max(bcube1.d[d][0], bcube2.d[d][0])), shared_max(min(bcube1.d[d][1], bcube2.d[d][1]));
			
			if (shared_max - shared_min > min_edge_dist) { // can connect with single road segment in dim !d, if the terrain in between is passable
				float const val1(shared_min + 0.5*min_edge_dist), val2(shared_max - 0.5*min_edge_dist); // shrink by half of min_edge_dist
				float best_conn_pos(0.0), best_cost(-1.0);
				bool is_4way1(0), is_4way2(0);

				if (city_params.make_4_way_ints) { // currently only inserts connector roads that have 3-way intersections on one end and 4-way intersections on the other end
					for (unsigned r12 = 0; r12 < 2; ++r12) {
						vector<road_t> const &roads((r12 ? rn2 : rn1).get_roads());

						for (auto r = roads.begin(); r != roads.end(); ++r) {
							if (r->dim == (d != 0)) continue; // wrong dim
							if (r->d[d][0] < val1 || r->d[d][1] > val2) continue; // road not contained in placement range
							float const conn_pos(r->get_center_dim(d)); // road centerline position along the edge of the plot
							float const cost(0.5*global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2,
								city1, city2, city1, city2, crc, road_width, conn_pos, !d, 1, (r12==0), (r12!=0))); // check_only=1; half cost (prefer over 3-way intersection)
							
							if (cost >= 0.0 && (best_cost < 0.0 || cost < best_cost)) {
								best_conn_pos = conn_pos; best_cost = cost;
								is_4way1 = (r12==0); is_4way2 = (r12!=0); // make only one end a 4-way intersection
								//is_4way1 = is_4way2 = 1; // make both ends 4-way intersections and move the roads to make them connect (WIP)
							}
						} // for r
					} // for r12
				}
				if (city_params.make_4_way_ints < 2) { // include connector roads that have 3-way intersections on both ends
					for (unsigned n = 0; n < city_params.num_conn_tries; ++n) { // make up to num_tries attempts at connecting the cities with a straight line
						float const conn_pos(rgen_uniform(val1, val2, rgen)); // chose a random new connection point and try it
						float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, city1, city2, city1, city2, crc, road_width, conn_pos, !d, 1, 0, 0)); // check_only=1
						if (cost >= 0.0 && (best_cost < 0.0 || cost < best_cost)) {best_conn_pos = conn_pos; best_cost = cost; is_4way1 = is_4way2 = 0;}
					}
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2,
						city1, city2, city1, city2, crc, road_width, best_conn_pos, !d, 0, is_4way1, is_4way2)); // check_only=0; make change
					assert(cost >= 0.0);
					return cost;
				}
			}
		} // for d
		point const center1(bcube1.get_cube_center()), center2(bcube2.get_cube_center());
		bool const dx(center1.x < center2.x), dy(center1.y < center2.y);
		cube_t const bc[2] = {bcube1, bcube2};
		
		if ((bc[dx].x1() - bc[!dx].x1() > min_jog) && (bc[dy].y1() - bc[!dy].y1() > min_jog)) {
			// connect with two road segments using a jog: Note: assumes cities are all about the same size
			bool const inv_dim(rgen.rand_bool());
			cube_t bc1_conn(bcube1), bc2_conn(bcube2);
			bc1_conn.d[0][ dx] = bcube2.d[0][!dx];
			bc2_conn.d[0][!dx] = bcube1.d[0][ dx];
			bc1_conn.d[1][ dy] = bcube2.d[1][!dy];
			bc2_conn.d[1][!dy] = bcube1.d[1][ dy];

			for (unsigned d = 0; d < 2; ++d) { // x-then-y vs. y-then-x
				bool const fdim((d != 0) ^ inv_dim); // first segment dim: 0=x first, 1=y first
				bool const range_dir1(fdim ? dy : dx), range_dir2(fdim ? dx : dy);
				cube_t region1(bc1_conn), region2(bc2_conn); // regions of valid jog connectors
				region1.d[ fdim][!range_dir1] = region1.d[ fdim][ range_dir1]; // fixed edge
				region2.d[!fdim][ range_dir2] = region2.d[!fdim][!range_dir2]; // fixed edge
				region1.d[!fdim][0] += min_edge_dist; region1.d[!fdim][1] -= min_edge_dist; // variable edges
				region2.d[ fdim][0] += min_edge_dist; region1.d[ fdim][1] -= min_edge_dist; // variable edges
				float const xmin(region1.d[!fdim][0]), xmax(region1.d[!fdim][1]), ymin(region2.d[fdim][0]), ymax(region2.d[fdim][1]); // Note: x and y may be swapped!
				float best_xval(0.0), best_yval(0.0), best_cost(-1.0);
				cube_t best_int_cube;
				bool const is_4way(city_params.make_4_way_ints > 0);

				if (is_4way) {
					vector<road_t> const &roads1(rn1.get_roads()), &roads2(rn2.get_roads());
					bool const d1(!fdim), d2(fdim); // dims for {city1, city2}

					for (auto r1 = roads1.begin(); r1 != roads1.end(); ++r1) {
						if (r1->dim == (d1 != 0)) continue; // wrong dim
						if (r1->d[d1][0] < bcube1.d[d1][0]+min_edge_dist || r1->d[d1][1] > bcube1.d[d1][1]-min_edge_dist) continue; // not an interior road (edge road)
						if (r1->d[d1][0] < (fdim ? xmin : ymin) || r1->d[d1][1] > (fdim ? xmax : ymax)) continue; // road not contained in placement range
						float const cpos1(r1->get_center_dim(d1)); // fdim=0 => yval
						
						for (auto r2 = roads2.begin(); r2 != roads2.end(); ++r2) {
							if (r2->dim == (d2 != 0)) continue; // wrong dim
							if (r2->d[d2][0] < bcube2.d[d2][0]+min_edge_dist || r2->d[d2][1] > bcube2.d[d2][1]-min_edge_dist) continue; // not an interior road (edge road)
							if (r2->d[d2][0] < (fdim ? ymin : xmin) || r2->d[d2][1] > (fdim ? ymax : xmax)) continue; // road not contained in placement range
							float const cpos2(r2->get_center_dim(d2)); // fdim=0 => xval
							float const xval(fdim ? cpos1 : cpos2), yval(fdim ? cpos2 : cpos1);
							try_single_jog_conn_road(city1, city2, blockers, crc, road_width, fdim, xval, yval, 1, best_xval, best_yval, best_cost, best_int_cube);
						} // for r2
					} // for r1
				}
				else {
					for (unsigned n = 0; n < city_params.num_conn_tries; ++n) { // make up to num_tries attempts at connecting the cities with a single jog
						float xval(rgen_uniform(xmin, xmax, rgen)), yval(rgen_uniform(ymin, ymax, rgen));
						if (!fdim) {swap(xval, yval);}
						try_single_jog_conn_road(city1, city2, blockers, crc, road_width, fdim, xval, yval, 0, best_xval, best_yval, best_cost, best_int_cube);
					} // for n
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					crc.hq.flatten_region_to(best_int_cube, city_params.road_border); // do this first to improve flattening
					unsigned road_ix[2] = {};
					road_ix[ fdim] = global_rn.num_roads();
					float const cost1(global_rn.create_connector_road(bcube1, best_int_cube, blockers, &rn1, nullptr, city1,
						CONN_CITY_IX, city1, city2, crc, road_width, (fdim ? best_xval : best_yval),  fdim, 0, is_4way, 0)); // check_only=0
					assert(cost1 >= 0.0);
					flatten_op_t const fop(crc.hq.last_flatten_op); // cache for reuse later during decrease_only pass
					road_ix[!fdim] = global_rn.num_roads();
					float const cost2(global_rn.create_connector_road(best_int_cube, bcube2, blockers, nullptr, &rn2,
						CONN_CITY_IX, city2, city1, city2, crc, road_width, (fdim ? best_yval : best_xval), !fdim, 0, 0, is_4way)); // check_only=0
					assert(cost2 >= 0.0);
					global_rn.create_connector_bend(best_int_cube, (dx ^ fdim), (dy ^ fdim), road_ix[0], road_ix[1]);
					// decrease_only=1; remove any dirt that the prev road added
					crc.hq.flatten_sloped_region(fop.x1, fop.y1, fop.x2, fop.y2, fop.z1, fop.z2, fop.dim, fop.border, fop.skip_six, fop.skip_eix, 0, 1);
					crc.hq.flatten_region_to(best_int_cube, city_params.road_border, 1); // one more pass to fix mesh that was raised above the intersection by a sloped road segment
					return (cost1 + cost2);
				}
			} // for d
		} // end two segment/one jog case
		if (city_params.make_4_way_ints > 2) {
			// determine if the two cities project onto each other in one dim; if so, maybe we can use the new city connector algorithm to join them with a jog (3 segments + 2 bends)
			cube_t bc1(bcube1), bc2(bcube2);
			bc1.expand_by(-min_jog);
			bc2.expand_by(-min_jog);

			if ((bc1.x1() < bc2.x2() && bc1.x2() > bc2.x1()) || (bc1.y1() < bc2.y2() && bc1.y2() > bc2.y1())) {
				return connect_two_cities_new(city1, city2, blockers, crc, road_width);
			}
		}
		return 0.0; // failed to connect cities
	}

	bool try_add_transmission_line(unsigned city1, unsigned city2, point const &p1, point const &p2, city_road_connector_t &crc,
		vect_cube_t &blockers, float road_width, float road_spacing)
	{
		float const tower_height(2.0*road_width);
		transmission_line_t tline(city1, city2, tower_height, p1, p2);
		if (!crc.route_transmission_line(tline, blockers, road_width, road_spacing)) return 0;
		transmission_lines.push_back(tline);
		return 1;
	}
	bool connect_two_cities_with_power(unsigned city1, unsigned city2, city_road_connector_t &crc, vect_cube_t &blockers, float road_width, float road_spacing) {
		assert(city1 != city2); // check for self reference
		road_network_t &rn1(get_city_by_ix(city1)), &rn2(get_city_by_ix(city2));
		cube_t const &bcube1(rn1.get_bcube()), &bcube2(rn2.get_bcube());
		point const center1(bcube1.get_cube_center()), center2(bcube2.get_cube_center());
		vector3d const dir(center2 - center1);
		float tower_to_city_spacing(2.0*road_width);

		for (unsigned d = 0; d < 2; ++d) { // x/y
			float const seg_min(max(bcube1.d[d][0], bcube2.d[d][0])), seg_max(min(bcube1.d[d][1], bcube2.d[d][1]));
			if (seg_max <= seg_min) continue; // no shared edge
			bool const other_dir(dir[!d] > 0);
			float const spacing_from_city(tower_to_city_spacing*(other_dir ? 1.0 : -1.0));
			unsigned const num_samples(unsigned(ceil((seg_max - seg_min)/road_spacing)));

			for (unsigned n = 0; n < num_samples; ++n) { // take some random points in the shared range
				point p1, p2;
				p1[d] = p2[d] = rgen.rand_uniform(seg_min, seg_max);
				// connect to city bcube edges facing each other, shifted away from the city
				float const edge1(bcube1.d[!d][other_dir]), edge2(bcube2.d[!d][!other_dir]);
				if (fabs(edge2 - edge1) < 3.0*tower_to_city_spacing) continue; // too close; shouldn't happen, but if it does the results won't be good
				p1[!d] = edge1 + spacing_from_city;
				p2[!d] = edge2 - spacing_from_city;
				p1.z   = bcube1.z2(); p2.z = bcube2.z2(); // set to city zvals; will be reset later anyway
				if (try_add_transmission_line(city1, city2, p1, p2, crc, blockers, road_width, road_spacing)) return 1;
			} // for n
		} // for d
		// no projecting points could be connected, find closest pair of corners
		tower_to_city_spacing /= SQRT2; // decrease because corner has spacing in both dims
		point p1, p2;

		for (unsigned d = 0; d < 2; ++d) { // x/y
			bool const D(dir[d] > 0);
			float const edge1( bcube1.d[d][D]), edge2(bcube2.d[d][!D]);
			if (fabs(edge2 - edge1) < 3.0*tower_to_city_spacing) return 0; // too close; shouldn't happen, but if it does the results won't be good
			float const spacing_from_city(tower_to_city_spacing*(D ? 1.0 : -1.0));
			p1[d] = edge1 + spacing_from_city; // opposite corners
			p2[d] = edge2 - spacing_from_city;
		}
		p1.z = bcube1.z2(); p2.z = bcube2.z2(); // set to city zvals; will be reset later anyway
		return try_add_transmission_line(city1, city2, p1, p2, crc, blockers, road_width, road_spacing);
	}
private:
	void try_single_jog_conn_road(unsigned city1, unsigned city2, vect_cube_t &blockers, city_road_connector_t &crc, float road_width,
		bool fdim, float xval, float yval, bool is_4way, float &best_xval, float &best_yval, float &best_cost, cube_t &best_int_cube)
	{
		float const height(crc.hq.get_road_zval_at_pt(point(xval, yval, 0.0))), half_width(0.5*road_width);
		cube_t const int_cube(xval-half_width, xval+half_width, yval-half_width, yval+half_width, height, height); // the candidate intersection point
		if (has_bcube_int_xy(int_cube, blockers)) return; // bad intersection, fail
		road_network_t &rn1(get_city_by_ix(city1)), &rn2(get_city_by_ix(city2));
		float const cost1(global_rn.create_connector_road(rn1.get_bcube(), int_cube, blockers, &rn1, nullptr, city1,
			CONN_CITY_IX, city1, city2, crc, road_width, (fdim ? xval : yval), fdim, 1, is_4way, 0)); // check_only=1
		if (cost1 < 0.0) return; // bad segment
		if (best_cost > 0.0 && cost1 > best_cost) return; // bound - early terminate
		float const cost2(global_rn.create_connector_road(int_cube, rn2.get_bcube(), blockers, nullptr, &rn2,
			CONN_CITY_IX, city2, city1, city2, crc, road_width, (fdim ? yval : xval), !fdim, 1, 0, is_4way)); // check_only=1
		if (cost2 < 0.0) return; // bad segment
		float const cost(cost1 + cost2);
		if (best_cost < 0.0 || cost < best_cost) {best_xval = xval; best_yval = yval; best_int_cube = int_cube; best_cost = cost;}
	}
	static float city_dist_sq(road_network_t const &c1, road_network_t const &c2) {
		return p2p_dist_sq(c1.get_bcube().get_cube_center(), c2.get_bcube().get_cube_center()); // just compare centers
	}
public:
	void connect_all_cities(float *heightmap, unsigned xsize, unsigned ysize, float road_width, float road_spacing) {
		if (road_width == 0.0 || road_spacing == 0.0) return; // no roads
		unsigned const num_cities(road_networks.size());
		if (num_cities < 2) return; // no cities to connect
		timer_t timer("Connect Cities");
		assert(heightmap != nullptr); // must be called when heightmap is valid
		heightmap_query_t hq(heightmap, xsize, ysize);
		city_road_connector_t crc(hq);
		vector<unsigned> is_conn(num_cities, 0); // start with all cities unconnected (0=unconnected, 1=connected, 2=done/connect failed
		vector<pair<float, unsigned>> cands;
		vect_cube_t blockers; // existing cities and connector roads that we want to avoid intersecting
		// gather city blockers
		get_city_bcubes(blockers);
		expand_cubes_by_xy(blockers, road_spacing); // separate roads by at least this value
		vect_cube_t const city_blockers(blockers); // save before adding roads, for use with railroad tracks
		unsigned const city_bcubes_end(blockers.size());
		float tot_cost(0.0);
		unsigned num_conn(0);

		// full cross product road connectivity
		for (unsigned i = 0; i < num_cities; ++i) {
			for (unsigned j = i+1; j < num_cities; ++j) {
				float const cost(connect_two_cities(i, j, blockers, crc, road_width));
				if (cost == 0.0) continue;
				tot_cost += cost;
				road_networks[i].register_connected_city(j);
				road_networks[j].register_connected_city(i);
				++num_conn;
			} // for j
		} // for i
		assign_city_clusters();
		global_rn.split_connector_roads(road_spacing);
		// place railroad tracks after connector roads are placed and split into segs so that intersections will be properly handled
		global_rn.gen_railroad_tracks(city_blockers, hq);
		global_rn.calc_bcube_from_roads(); // must be after placing tracks
		global_rn.finalize_bridges_and_tunnels();
		
		// only connect cities with transmission lines if there are no secondary buildings in the way;
		// while buildings should now avoid tlines, it still looks a bit odd having them so close to houses, and the endpoints aren't checked;
		// but let the user override this by setting add_transmission_lines to 1
		if (city_params.add_tlines && (city_params.add_tlines == 1 || !have_secondary_buildings())) {
			for (auto i = blockers.begin(); i != blockers.begin() + city_bcubes_end; ++i) {i->expand_by_xy(-road_spacing);} // undo city expand
			connect_city_power_grids(crc, blockers, road_width, road_spacing);
		}
		timer.end();
		// old: 8, 12, 7, 19057 ; new: 8, 15, 7, 26085
		cout << "Cities: " << num_cities << ", connector roads: " << num_conn << ", transmission lines: " << transmission_lines.size() << ", total cost: " << tot_cost << endl;
	}

	struct conn_cand_t {
		unsigned c1, c2; // city index
		float dsq; // squared distance between cities
		conn_cand_t(unsigned c1_, unsigned c2_, float dsq_) : c1(c1_), c2(c2_), dsq(dsq_) {}
		bool operator<(conn_cand_t const &c) const {return (dsq < c.dsq);} // compare squared distance only as it will likely always be unique
	};
	void connect_city_power_grids(city_road_connector_t &crc, vect_cube_t &blockers, float road_width, float road_spacing) {
		unsigned const num_cities(road_networks.size());
		if (num_cities < 2) return; // no cities to connect
		vector<uint8_t> power_connected(num_cities, 0);
		vector<conn_cand_t> cands;

		// determine first two cities to connect
		for (unsigned i = 0; i < num_cities; ++i) {
			for (unsigned j = i+1; j < num_cities; ++j) {
				cands.emplace_back(i, j, city_dist_sq(road_networks[i], road_networks[j]));
			}
		}
		sort(cands.begin(), cands.end()); // sort min to max distance
		bool success(0);

		for (auto const &cand : cands) {
			if (connect_two_cities_with_power(cand.c1, cand.c2, crc, blockers, road_width, road_spacing)) {
				power_connected[cand.c1] = power_connected[cand.c2] = 1;
				success = 1;
				break;
			}
		}
		if (!success) return; // can't even connect two cities together, failed
		unsigned num_power_conn(2);

		// connect remaining cities
		while (num_power_conn < num_cities) {
			cands.clear();
			success = 0;

			for (unsigned i = 0; i < num_cities; ++i) {
				if (!power_connected[i]) continue; // not connected, skip

				for (unsigned j = 0; j < num_cities; ++j) {
					if (power_connected[j]) continue; // connected, skip (includes self)
					cands.emplace_back(i, j, city_dist_sq(road_networks[i], road_networks[j]));
				}
			} // for i
			sort(cands.begin(), cands.end()); // sort min to max distance

			for (auto const &cand : cands) {
				if (connect_two_cities_with_power(cand.c1, cand.c2, crc, blockers, road_width, road_spacing)) {
					power_connected[cand.c2] = 1;
					success = 1;
					break;
				}
			}
			if (!success) break; // failed to connect any more pairs, done (partially connected)
			++num_power_conn;
		} // end while()
	}
	void closest_edge_power_pole_conn_pts(point const &pt, unsigned city_ix, point conn_pts[3]) const {
		vector<power_pole_t> const &ppoles(get_city_by_ix(city_ix).get_power_poles());
		if (ppoles.empty()) {UNROLL_3X(conn_pts[i_] = pt;) return;} // set all conn pts to pt if there are no power poles
		float dmin_sq(0.0);

		for (auto const &p : ppoles) {
			if (!p.is_at_grid_edge()) continue; // not on edge
			point cand_pts[3];
			// Note: it's possible that the wires connecting from the power pole to the tline tower cross over each other or intersect the power pole,
			// depending on the orientations of the two sets of wires relative to each other; I have no idea how to avoid this
			p.get_top_wires_conn_pts(cand_pts);
			float const dsq(p2p_dist_xy_sq(pt, cand_pts[1])); // use center wire
			if (dmin_sq == 0.0 || dsq < dmin_sq) {dmin_sq = dsq; UNROLL_3X(conn_pts[i_] = cand_pts[i_];)}
		} // for p
	}
	void connect_power_poles_to_transmission_lines() {
		for (auto &t : transmission_lines) {
			closest_edge_power_pole_conn_pts(t.p1, t.city1, t.p1_wire_pts);
			closest_edge_power_pole_conn_pts(t.p2, t.city2, t.p2_wire_pts);
			t.calc_bcube();
		}
	}
	void add_streetlights() {
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->add_streetlights();}
	}
	void gen_tile_blocks() {
		timer_t timer("Gen Tile Blocks");
		global_rn.gen_tile_blocks(); // must be done first to fill in road_to_city and city_to_seg
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->gen_tile_blocks();}
		unsigned global_plot_id(0);
		global_rn.calc_ix_values(road_networks, global_rn, global_plot_id);
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->calc_ix_values(road_networks, global_rn, global_plot_id);}
	}
	void gen_parking_lots_and_place_objects(vector<car_t> &cars, bool have_cars) {
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->gen_parking_lots_and_place_objects(cars, have_cars, have_plot_dividers);}
	}
	void get_city_bcubes(vect_cube_t &bcubes) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {bcubes.push_back(r->get_bcube());}
	}
	int check_city_contains_overlaps(cube_t const &query) const { // XY only; 0=no intersection, 1=overlaps, 2=contains; ignores global connector RN
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {
			if (r->get_bcube().contains_cube_xy(query)) return 2;
			if (r->get_bcube().intersects_xy   (query)) return 1;
		}
		return 0;
	}
	bool check_inside_city(point const &pos, float radius) const { // Note: pos is in camera space
		cube_t query; query.set_from_sphere((pos - get_camera_coord_space_xlate()), radius);

		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {
			if (r->get_bcube().contains_cube_xy(query)) return 1;
		}
		return 0;
	}
	bool update_depth_if_underwater(point const &pos, float &depth) const { // Note: pos is in camera space
		point const pos_bs(pos - get_camera_coord_space_xlate());

		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {
			if (r->update_depth_if_underwater(pos_bs, depth)) return 1;
		}
		return 0;
	}
	void get_all_road_bcubes(vect_cube_t &bcubes, bool connector_only) const {
		global_rn.get_road_bcubes(bcubes); // not sure if this should be included
		if (connector_only) return;
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_road_bcubes(bcubes);}
	}
	void get_all_plot_zones(vect_city_zone_t &zones) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_plot_zones(zones);}
	}
	bool check_road_sphere_coll(point const &pos, float radius, bool xy_only, bool exclude_bridges_and_tunnels) const { // Note: pos is in camera space
		if (global_rn.check_road_sphere_coll(pos, radius, 1, xy_only, exclude_bridges_and_tunnels)) return 1;
		point const query_pos(pos - get_camera_coord_space_xlate()); // convert from camera space to global city space

		// since transmission lines run between cities, we can check the cities_bcube for early rejection
		if (!exclude_bridges_and_tunnels && xy_only && !transmission_lines.empty() && sphere_cube_intersect_xy(query_pos, radius, cities_bcube)) {
			// assume this is a valid scenery placement query and check transmission lines
			for (auto const &t : transmission_lines) {
				if (t.sphere_intersect_xy(query_pos, radius)) return 1;
			}
		}
		return 0;
	}
	bool check_tline_cube_intersect_xy(cube_t const &c) const {
		if (transmission_lines.empty())     return 0;
		if (!cities_bcube.intersects_xy(c)) return 0;

		for (auto const &t : transmission_lines) {
			if (t.cube_intersect_xy(c)) return 1;
		}
		return 0;
	}
	void get_roads_in_region(cube_t const &region, vect_cube_t &out, vect_cube_t &out_bt) const {
		global_rn.get_roads_in_region(region, out, out_bt);
	}
	bool choose_pt_in_park(point const &pos, point &park_pos, rand_gen_t &rgen) const {
		for (road_network_t const &r : road_networks) {
			if (!r.get_bcube().contains_pt_xy(pos))  continue; // point not in this city
			if (r.choose_pt_in_park(park_pos, rgen)) return 1;
		}
		return 0;
	}
	bool check_mesh_disable(point const &pos, float radius) const {return global_rn.check_mesh_disable(pos, radius);}
	bool tile_contains_tunnel(cube_t const &bcube) const {return global_rn.tile_contains_tunnel(bcube);}

	int get_color_at_xy(point const &pos, colorRGBA &color) const {
		for (road_network_t const &r : road_networks) {
			int const ret(r.get_color_at_xy(pos, color));
			if (ret) return ret;
		}
		return global_rn.get_color_at_xy(pos, color);
	}
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, float prev_frame_zval, vector3d *cnorm) const { // pos is in camera space
		vector3d const xlate(get_camera_coord_space_xlate());
		float const dist(p2p_dist(pos, p_last));

		for (road_network_t const &r : road_networks) {
			if (r.proc_sphere_coll(pos, p_last, xlate, dist, radius, prev_frame_zval, cnorm)) return 1;
		}
		return global_rn.proc_sphere_coll(pos, p_last, xlate, dist, radius, prev_frame_zval, cnorm); // needed for bridges and tunnels
	}
	bool line_intersect(point const &p1, point const &p2, float &t) const {
		bool ret(0);
		for (road_network_t const &r : road_networks) {ret |= r.line_intersect(p1, p2, t);}
		ret |= global_rn.line_intersect(p1, p2, t); // bridges and tunnels
		return ret;
	}
	void add_city_lights(vector3d const &xlate, cube_t &lights_bcube) const {
		global_rn.add_city_lights(xlate, lights_bcube); // no streetlights, but may need to add lights for bridges and tunnels
		for (road_network_t const &r : road_networks) {r.add_city_lights(xlate, lights_bcube);}
	}
	void get_occluders(vect_cube_t &occluders) const {
		for (road_network_t const &r : road_networks) {r.get_occluders(occluders, dstate.xlate);}
	}
	vector<bridge_t> const &get_bridges() const {return global_rn.get_bridges();}

	bool have_animations() const {
		for (road_network_t const &r : road_networks) {
			if (r.have_animations()) return 1;
		}
		return global_rn.have_animations();
	}
	void draw(int trans_op_mask, vector3d const &xlate, bool use_dlights, bool shadow_only) { // non-const because dstate/qbd is modified
		if (road_networks.empty() && global_rn.empty()) return;

		if (trans_op_mask & 1) { // opaque pass, should be first
			//highres_timer_t timer(shadow_only ? "Draw City Shadows" : "Draw City"); // 1.1ms / 0.42ms shadows
			fgPushMatrix();
			translate_to(xlate);
			set_std_depth_func_with_eq(); // helps prevent Z-fighting

			if (1 || have_plot_dividers) { // enable normal maps for fences and walls; also applies to tunnels and power poles
				dstate.set_enable_normal_map(1);
				bind_default_flat_normal_map(); // set flat normal map texture as the default
			}
			if (have_animations()) {enable_animations_for_shader(dstate.s);}
			if (!shadow_only) {enable_dlight_bcubes |= city_lights_custom_bcube;}
			dstate.pre_draw(xlate, use_dlights, shadow_only, 1); // always_setup_shader=1
			assert(dstate.s.is_setup());
			for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->draw(dstate, shadow_only, 0);}
			global_rn.draw(dstate, shadow_only, 1); // connector road may have bridges, and therefore needs shadows
			draw_transmission_lines();
			dstate.post_draw();
			if (!shadow_only) {enable_dlight_bcubes = 0;}
			set_std_depth_func();
			fgPopMatrix();
		}
		if (trans_op_mask & 2) {dstate.draw_and_clear_light_flares();} // transparent pass; must be done last for alpha blending, and no translate
	}
	void draw_transparent(vector3d const &xlate, bool use_dlights) { // non-const because dstate/qbd is modified
		if (road_networks.empty()) return;
		fgPushMatrix();
		translate_to(xlate);
		bind_default_flat_normal_map();
		enable_dlight_bcubes |= city_lights_custom_bcube;
		dstate.pre_draw(xlate, use_dlights, 0, 1); // shadow_only=0, always_setup_shader=1
		// bind a valid but empty shadow map texture for both the sun and moon in case it was left in an invalid state and the first tile doesn't have a valid smap
		for (unsigned i = 0; i < 2; ++i) {bind_texture_tu(get_empty_smap_tid(), 6+i);}
		for (road_network_t &r : road_networks) {r.draw_transparent(dstate);}
		dstate.post_draw();
		enable_dlight_bcubes = 0;
		fgPopMatrix();
	}
	void draw_transmission_lines() { // non-const because dstate is modified
		if (transmission_lines.empty()) return;
		//highres_timer_t timer("Draw Transmission Lines"); // 0.12ms
		dstate.set_untextured_material();
		for (auto const &i : transmission_lines) {dstate.draw_transmission_line(i);}
		dstate.unset_untextured_material();
	}
	void draw_label() {dstate.show_label_text();}

	// cars/peds
	void next_frame() {
		if (!animate2) return;
		//timer_t timer("Update Stoplights");
		//for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {cout << r->get_traffic_density() << " ";} cout << endl;
		for (road_network_t &r : road_networks) {r.next_frame();}
		global_rn.next_frame(); // not needed since there are no 3/4-way intersections/stoplights?
	}
	void register_car_at_city(unsigned city_id) const {get_city(city_id).register_car();} // Note: must be const
	
	bool add_car(car_t &car, rand_gen_t &rgen) const {
		if (road_networks.empty()) return 0; // no cities to add cars to
		return road_network_t::add_car_to_rns(car, rgen, road_networks);
	}
	bool gen_ped_pos(pedestrian_t &ped, rand_gen_t &rgen) const {
		if (road_networks.empty()) return 0; // no cities to add peds to
		return road_network_t::gen_ped_pos(ped, rgen, road_networks);
	}
	bool is_city_residential(unsigned city_id) const {return get_city(city_id).get_is_residential();}
	road_plot_t const &get_plot_from_global_id(unsigned city_id, unsigned global_plot_id) const {return get_city(city_id).get_plot_from_global_id(global_plot_id);}
	int get_global_plot_id_for_pos(unsigned city_id, point const &pos) const {return get_city(city_id).get_global_plot_id_for_pos(pos);}
	unsigned get_next_plot(unsigned city_id, unsigned plot, unsigned dest_plot, int exclude_plot) const {return get_city(city_id).get_next_plot(plot, dest_plot, exclude_plot);}
	bool choose_dest_building(unsigned city_id, unsigned &plot, unsigned &building, rand_gen_t &rgen) const {return get_city(city_id).choose_dest_building(plot, building, rgen);}
	
	bool update_car_dest(car_t &car) const {
		if (car.is_parked()) return 0; // no dest for parked cars
		if (car.dest_valid && !car_at_dest(car)) return 0; // not yet at destination, keep existing dest
		assert(!car.dest_valid || car.dest_city == car.cur_city); // sanity check
		static rand_gen_t rgen; // reused across calls
		choose_new_car_dest(car, rgen);
		return 1;
	}
	void choose_new_car_dest(car_t &car, rand_gen_t &rgen) const {
		if (road_networks.empty()) {assert(!car.dest_valid); return;} // no roads, no updates
		if (!car.dest_valid) {car.dest_city = car.cur_city;} // start in current city
		else {
			auto const &conn(get_city_by_ix(car.cur_city).get_connected());
			float const new_city_prob(city_params.new_city_prob*min(0.4f, 0.1f*conn.size())); // 10% to 40% chance, depending on the number of connecting cities (to reduce traffic congestion)

			if (rgen.rand_probability(new_city_prob)) { // select a different city when there are multiple cities
				if (rgen.rand_float() < city_params.traffic_balance_val) { // choose the connected city with the lowest traffic density
					float min_td(0.0);

					for (auto c = conn.begin(); c != conn.end(); ++c) {
						float const td(get_city_by_ix(*c).get_traffic_density()); // excludes global_rn
						if (min_td == 0.0 || td < min_td) {min_td = td; car.dest_city = *c;}
					}
				}
				else { // choose a randomly selected connected city
					vector<unsigned> const cands(conn.begin(), conn.end()); // copy set to vector; should not include car.dest_city
					car.dest_city = cands[rgen.rand() % cands.size()]; // choose a random connected (adjacent) city
				}
			}
		}
		car.dest_valid = get_city(car.dest_city).choose_new_car_dest(car, rgen);
	}
	bool car_at_dest(car_t const &car) const {
		if (!car.dest_valid) return 0;
		return get_city(car.dest_city).car_at_dest(car);
	}
	road_network_t const &get_car_rn(car_base_t const &car) const {return road_network_t::get_car_rn(car, road_networks, global_rn);}
	
	void update_car(car_t &car, vector<car_t> const &cars, rand_gen_t &rgen) const {
		if (car.cur_city == NO_CITY_IX) return; // not in a city (in a garage), nothing to update
		//update_car_seg_stats(car); // not needed - stats not yet used
		get_car_rn(car).update_car(car, cars, rgen, road_networks, global_rn);
		if (city_params.enable_car_path_finding) {update_car_dest(car);}
	}
	void update_car_seg_stats(car_base_t const &car) const {get_car_rn(car).update_car_seg_stats(car);}
	road_isec_t const &get_car_isec(car_base_t const &car) const {return get_car_rn(car).get_car_isec(car);}
	cube_t get_road_bcube_for_car(car_base_t const &car) const {return get_car_rn(car).get_road_bcube_for_car(car);}
	virtual cube_t get_bcube_for_car(car_base_t const &car) const {return get_road_bcube_for_car(car);}
	
	dw_query_t get_nearby_driveway(unsigned city_ix, unsigned global_plot_id, point const &pos, float dist) const {
		return get_city(city_ix).get_nearby_driveway(global_plot_id, pos, dist);
	}
}; // city_road_gen_t


// Note: the car_manager_t member functions that use road_gen are here rather than in cars.cpp
cube_t car_manager_t::get_cb_bcube(car_block_t const &cb )       const {return road_gen.get_city_bcube_for_cars(cb.cur_city);}
road_isec_t const &car_manager_t::get_car_isec(car_t const &car) const {return road_gen.get_car_isec(car);}
bool car_manager_t::check_collision(car_t &c1, car_t &c2)        const {return c1.check_collision(c2, road_gen);}
void car_manager_t::register_car_at_city(car_t const &car) {road_gen.register_car_at_city(car.cur_city);}
cube_t const &car_manager_t::get_car_dest_bcube(car_t const &car, bool isec_only) const {return road_gen.get_car_dest_bcube(car, isec_only);}
int car_manager_t::get_parking_lot_ix_for_car(car_t const &car) const {return road_gen.get_parking_lot_ix_for_car(car);}

void car_manager_t::add_car() {
	car_t car;
	if (road_gen.add_car(car, rgen)) {cars.push_back(car);}
}
void car_manager_t::update_cars() {
	for (auto i = cars.begin(); i != cars.end(); ++i) {road_gen.update_car(*i, cars, rgen);} // run update logic
}
void car_manager_t::get_car_ix_range_for_cube(vector<car_block_t>::const_iterator cb, cube_t const &bc, unsigned &start, unsigned &end) const {
	start = cb->start; end = (cb+1)->start;
	assert(end <= cars.size());
	if (!road_gen.cube_overlaps_pl_or_dw_xy(bc, cb->cur_city)) {end   = cb->first_parked;} // moving cars only (beginning of range)
	if (!road_gen.cube_overlaps_road_xy    (bc, cb->cur_city)) {start = cb->first_parked;} // parked cars only (end of range)
	assert(start <= end);
}
void car_manager_t::setup_occluders() {
	dstate.get_occluders().clear();
	if ((display_mode & 0x08) && !cars.empty()) {road_gen.get_occluders(dstate.get_occluders());}
}
vector<bridge_t> const &car_manager_t::get_bridges() const {return road_gen.get_bridges();}

struct cmp_light_source_sz_dist {
	point const &cpos;
	cmp_light_source_sz_dist(point const &p) : cpos(p) {}
	float get_value(light_source const &s) const {return s.get_beamwidth()*s.get_radius()*s.get_radius()/p2p_dist_sq(s.get_pos(), cpos);}
	bool operator()(light_source const &a, light_source const &b) const {return (get_value(a) > get_value(b));} // sort largest/closest to smallest/furthest
};

void sort_lights_by_dist_size(vector<light_source> &lights, point const &cpos) {
	stable_sort(lights.begin(), lights.end(), cmp_light_source_sz_dist(cpos));
}

void filter_dlights_to(vector<light_source> &lights, unsigned max_num, point const &cpos) {
	if (lights.size() <= max_num) return;
	sort_lights_by_dist_size(lights, cpos);
	lights.resize(max_num); // remove lowest scoring lights
}


// Note: these ped_manager_t functions are defined here because they use road_gen
road_plot_t const &ped_manager_t::get_city_plot_for_peds(unsigned city_ix, unsigned plot_ix) const {return road_gen.get_plot_from_global_id(city_ix, plot_ix);}
road_isec_t const &ped_manager_t::get_car_isec(car_base_t const &car) const {return road_gen.get_car_isec(car);}
int ped_manager_t::get_global_plot_id_for_pos(unsigned city_ix, point const &pos) const {return road_gen.get_global_plot_id_for_pos(city_ix, pos);}

dw_query_t ped_manager_t::get_nearby_driveway(unsigned city_ix, unsigned global_plot_ix, point const &pos, float dist) const {
	return road_gen.get_nearby_driveway(city_ix, global_plot_ix, pos, dist);
}
bool ped_manager_t::is_city_residential(unsigned city_ix) const {return road_gen.is_city_residential(city_ix);}
cube_t ped_manager_t::get_city_bcube_for_peds(unsigned city_ix) const {return road_gen.get_city_bcube(city_ix);}

cube_t ped_manager_t::get_expanded_city_bcube_for_peds(unsigned city_ix) const {
	cube_t bcube(get_city_bcube_for_peds(city_ix));
	expand_cube_for_ped(bcube);
	return bcube;
}
cube_t ped_manager_t::get_expanded_city_plot_bcube_for_peds(unsigned city_ix, unsigned plot_ix) const {
	cube_t bcube(road_gen.get_plot_from_global_id(city_ix, plot_ix));
	expand_cube_for_ped(bcube);
	bcube.expand_by_xy(city_params.road_width); // required to include pedestrians that are crossing the road and not contained in the plot bcube
	return bcube;
}
vect_cube_t const &ped_manager_t::get_colliders_for_plot(unsigned city_ix, unsigned plot_ix) const {return road_gen.get_colliders_for_plot(city_ix, plot_ix);}
bool ped_manager_t::gen_ped_pos(pedestrian_t &ped) {return road_gen.gen_ped_pos(ped, rgen);} // Note: non-const because rgen is modified

bool ped_manager_t::mark_crosswalk_in_use(pedestrian_t const &ped) {
	bool const dim(fabs(ped.dir.y) > fabs(ped.dir.x)), dir(ped.dir[dim] > 0); // something like this?
	return road_gen.get_city(ped.city).mark_crosswalk_in_use(ped.pos, dim, dir);
}
bool ped_manager_t::check_isec_sphere_coll(pedestrian_t const &ped, cube_t &coll_cube) const {
	return road_gen.get_city(ped.city).check_isec_sphere_coll(ped.pos, 0.6*ped.radius, coll_cube); // Note: no xlate is required since peds and city are in the same coord space
}
bool ped_manager_t::check_streetlight_sphere_coll(pedestrian_t const &ped, cube_t &coll_cube) const {
	return road_gen.get_city(ped.city).check_streetlight_sphere_coll_xy(ped.pos, ped.radius, coll_cube);
}
int ped_manager_t::get_road_ix_for_ped_crossing(pedestrian_t const &ped, bool road_dim) const { // returns -1 on failure (ped not in the road)
	return road_gen.get_city(ped.city).get_nearby_road_ix(ped.pos, road_dim);
}
int ped_manager_t::get_parking_lot_ix_for_ped(pedestrian_t const &ped, bool inc_driveways) const {
	return road_gen.get_city(ped.city).get_parking_lot_ix(ped.pos, inc_driveways);
}
void ped_manager_t::setup_occluders() {
	dstate.get_occluders().clear();
	if ((display_mode & 0x08) && !peds.empty()) {road_gen.get_occluders(dstate.get_occluders());}
}

// path finding
bool ped_manager_t::choose_dest_building_or_parked_car(pedestrian_t &ped) { // modifies rgen, non-const
	unsigned const prev_dest_plot(ped.dest_plot);
	ped.clear_current_dest(); // will choose a new dest

	if (city_params.num_cars == 0 || (rgen.rand() & 3) != 0) { // choose a dest building 75% of the time, 100% of the time if there are no cars
		ped.has_dest_bldg = road_gen.choose_dest_building(ped.city, ped.dest_plot, ped.dest_bldg, rgen);
	}
	if (city_params.num_cars > 0 && !ped.has_dest_bldg) { // chose a dest parked car 25% of the time, or if choosing a dest building failed
		ped.has_dest_car = choose_dest_parked_car(ped.city, ped.dest_plot, ped.dest_bldg, ped.dest_car_center);
		if (ped.has_dest_car) {ped.dest_plot = road_gen.get_city(ped.city).encode_plot_id(ped.dest_plot);}
	}
	bool const has_valid_dest(ped.has_dest_bldg || ped.has_dest_car);
	if (!has_valid_dest && prev_dest_plot > 0) {ped.dest_plot = prev_dest_plot;} // if a dest plot was selected, restore its value to prevent randomly jumping between plots each frame
	ped.next_plot = get_next_plot(ped);
	return has_valid_dest;
}
void ped_manager_t::choose_new_ped_plot_pos(pedestrian_t &ped) {
	if (city_params.ped_respawn_at_dest) { // respawn
		for (unsigned n = 0; n < 100; ++n) { // keep respawning until it's not visible by the camera
			float const prev_zval(ped.pos.z);
			bool const ret(road_gen.get_city(ped.city).gen_ped_pos(ped, rgen));
			ped.pos.z = prev_zval; // restore orig zval - don't want to change this (zval was set from ped radius post-model scale but should be pre-model scale)
			if (!ret) break; // failed to respawn, leave at current pos (should be very rare)
			float const draw_dist(500.0*get_ped_radius());
			if (!dist_less_than(get_camera_pos(), ped.pos, draw_dist) || !camera_pdu.sphere_visible_test(ped.pos, ped.radius)) break; // good pos
		}
		register_ped_new_plot(ped);
	}
	choose_dest_building_or_parked_car(ped);
}
unsigned ped_manager_t::get_next_plot(pedestrian_t &ped, int exclude_plot) const {return road_gen.get_next_plot(ped.city, ped.plot, ped.dest_plot, exclude_plot);}


void city_lights_manager_t::add_player_flashlight(float radius_scale) {
	if (!flashlight_on) return;
	add_player_flashlight_light_source(radius_scale);
	assert(!dl_sources.empty()); // must have been added
	float const zval(dl_sources.back().get_pos().z), light_radius(dl_sources.back().get_radius());
	min_eq(lights_bcube.z1(), (zval - light_radius));
	max_eq(lights_bcube.z2(), (zval + light_radius));
}

void city_lights_manager_t::tighten_light_bcube_bounds(vector<light_source> const &lights) {
	if (lights.empty()) return; // nothing to do
	cube_t tight_bcube;
	for (auto l = lights.begin(); l != lights.end(); ++l) {tight_bcube.assign_or_union_with_sphere(l->get_pos(), l->get_radius());}
	//cout << TXT(bcube.dx()) << TXT(bcube.dy()) << TXT(tight_bcube.dx()) << TXT(tight_bcube.dy()) << endl;
	if (!lights_bcube.intersects(tight_bcube)) {cerr << "Invalid bcubes: " << TXT(lights_bcube.str()) << TXT(tight_bcube.str()) << endl;}
	lights_bcube.intersect_with_cube(tight_bcube); // clip the original cube to the tight cube (better to just set to tight cube?)
}

void city_lights_manager_t::clamp_to_max_lights(vector3d const &xlate, vector<light_source> &lights) {
	unsigned const max_dlights(min(1024U, city_params.max_lights)); // Note: should be <= the value used in upload_dlights_textures()
	//cout << "dlights: " << lights.size() << ", bcube: " << lights_bcube.str() << endl; // 536/621/889

	if (lights.size() > max_dlights) {
		if (lights.size() > 4*max_dlights) { // too many lights, reduce light radius for next frame
			light_radius_scale *= 0.95;
			cout << "Too many city lights: " << lights.size() << ". Reducing light_radius_scale to " << light_radius_scale << endl;
		}
		filter_dlights_to(lights, max_dlights, (camera_pdu.pos - xlate));
	}
}

bool city_lights_manager_t::begin_lights_setup(vector3d const &xlate, float light_radius, vector<light_source> &lights) {
	assert(light_radius > 0.0);
	for (auto i = lights.begin(); i != lights.end(); ++i) {i->release_smap();} // must be done before clearing dlights
	clear_dynamic_lights();
	lights_bcube.set_to_zeros();
	dl_smap_enabled = 0; // here for safety, needed for buildings flow
	if (draw_model != 0) return 0; // skip shadow calculation in wireframe mode
	if (!enable_lights() && !prev_had_lights) return 0; // only have lights at night
	lights_bcube = cube_t(camera_pdu.pos - xlate);
	lights_bcube.expand_by(light_radius);
	lights_bcube.z1() =  FLT_MAX;
	lights_bcube.z2() = -FLT_MAX;
	return 1;
}

void city_lights_manager_t::finalize_lights(vector<light_source> &lights) { // Note: lights is always dl_sources and not passed into calls below
	add_dynamic_lights_city(lights_bcube, dlight_add_thresh, CITY_LIGHT_FALLOFF);
	upload_dlights_textures(lights_bcube, dlight_add_thresh);
	prev_had_lights = !lights.empty();
}

struct sel_smap_light_t {
	unsigned lix;
	bool matched_smap_id;
	sel_smap_light_t(unsigned lix_, bool matched) : lix(lix_), matched_smap_id(matched) {}
};
void city_lights_manager_t::setup_shadow_maps(vector<light_source> &light_sources, point const &cpos, unsigned max_smaps, bool sec_camera_mode) {
	unsigned const num_smaps(min((unsigned)light_sources.size(), min(max_smaps, MAX_DLIGHT_SMAPS)));
	dl_smap_enabled = 0;
	if (!enable_dlight_shadows || shadow_map_sz == 0 || num_smaps == 0) return;
	sort_lights_by_dist_size(light_sources, cpos); // Note: may already be sorted for enabled lights selection, but okay to sort again
	cmp_light_source_sz_dist sz_cmp(cpos);
	unsigned const smap_size(city_params.smap_size); // 0 = use default shadow map resolution
	// capture player pos in global coordinate space before replacing with light pos so it can be used for LOD during model drawing
	pre_smap_player_pos = cpos; // player or security camera (or maybe reflected pos in the future)
	if (!sec_camera_mode) {actual_player_pos = pre_smap_player_pos;} // actual_player_pos only applies to the player
	// Note: if using a dynamic (distance-based) sm_size, need to maintain a pool of different sm resolutions somehow
	check_gl_error(430);
	// Note: slow to recreate shadow maps every frame, but most city lights are either dynamic (headlights) or include dynamic shadow casters (cars) and need to be updated every frame anyway
	// Do we want to gradually fade in new shadow maps and fade out old ones? But how do we track which lights are associated with old shadow maps?
	// Tracking positions won't work for car headlights because they move. We don't have object pointers to track either. And what about lights that are no longer in our list?
	vector<sel_smap_light_t> selected;
	selected.reserve(num_smaps);

	for (auto i = light_sources.begin(); i != light_sources.end() && selected.size() < num_smaps; ++i) {
		if (i->has_no_shadows())       continue; // shadows not enabled for this light
		if (!i->is_very_directional()) continue; // not a spotlight
		if (sz_cmp.get_value(*i) < 0.002) break; // light influence is too low, skip even though we have enough shadow maps; can break because sort means all later lights also fail
		bool matched_smap_id(0);
		if (!i->alloc_shadow_map(matched_smap_id, smap_size)) break; // out of shadow maps, done
		selected.emplace_back((i - light_sources.begin()), matched_smap_id);
	} // for i
	// now that all smaps have been allocated, we can create them without worrying about the backing texture array getting resized and overwriting earlier shadow maps
	for (sel_smap_light_t const &s : selected) {light_sources[s.lix].update_shadow_map(s.matched_smap_id, CITY_LIGHT_FALLOFF);}
	dl_smap_enabled |= !selected.empty();
	check_gl_error(431);
}


class city_gen_t : public city_plot_gen_t, public city_lights_manager_t {

	city_road_gen_t road_gen;
	car_manager_t car_manager;
	ped_manager_t ped_manager;
	unsigned prev_city_lights_setup_frame;

public:
	city_gen_t() : car_manager(road_gen), ped_manager(road_gen, car_manager), prev_city_lights_setup_frame(-1) {}

	bool gen_city(city_params_t const &params, cube_t &cities_bcube) {
		unsigned x1(0), y1(0), x2(0), y2(0);
		if (!find_best_city_location(params.city_size_min, params.city_size_min, params.city_size_max, params.city_size_max,
			params.city_border, params.slope_width, max(params.num_samples, 1U), x1, y1, x2, y2)) return 0;
		float const elevation(flatten_region(x1, y1, x2, y2, params.slope_width));
		cube_t const pos_range(add_plot(x1, y1, x2, y2, elevation));
		if (cities_bcube.is_all_zeros()) {cities_bcube = pos_range;} else {cities_bcube.union_with_cube(pos_range);}
		
		if (params.roads_enabled()) {
			rand_gen_t rgen2; // don't use the built-in rgen to avoid affecting other state
			rgen2.set_state(x1, y1+y2);
			float const rp(params.residential_probability);
			bool const is_residential(rgen2.rand_probability(rp)), extra_spacing_dim(rgen2.rand_bool());
			vector2d road_spacing;

			for (unsigned d = 0; d < 2; ++d) {
				road_spacing[d] = params.road_spacing;
				if (params.road_spacing_rand   > 0.0) {road_spacing[d] *= 1.0f + params.road_spacing_rand*rgen2.rand_float();}
				if (params.road_spacing_xy_add > 0.0 && bool(d) == extra_spacing_dim) {road_spacing[d] *= 1.0 + params.road_spacing_xy_add;}
			}
			road_gen.gen_roads_and_plots(pos_range, params.road_width, road_spacing, is_residential);
		}
		return 1;
	}
	void gen_cities() { // Note: params better be city_params
		if (city_params.num_cities == 0) return;
		cube_t cities_bcube(all_zeros);
		{ // open a scope
			timer_t t("Choose City Location");
			for (unsigned n = 0; n < city_params.num_cities; ++n) {gen_city(city_params, cities_bcube);}
		}
		if (!cities_bcube.is_all_zeros()) {set_buildings_pos_range(cities_bcube);}
		road_gen.connect_all_cities(heightmap, xsize, ysize, city_params.road_width, city_params.road_spacing);
		road_gen.gen_tile_blocks();
		road_gen.add_streetlights();
		car_manager.init_cars(city_params.num_cars);
		init_city_spectate_manager(car_manager, ped_manager);
	}
	void gen_details() {
		if (road_gen.empty()) return; // nothing to do - no roads or cars
		vector<car_t> parked_cars;
		vect_cube_t hp_locs;
		bool const have_cars(!car_manager.empty());
		highres_timer_t timer("Gen City Details");
		road_gen.gen_parking_lots_and_place_objects(parked_cars, have_cars);
		road_gen.connect_power_poles_to_transmission_lines(); // must be after placing power poles
		if (city_params.has_helicopter_model()) {get_all_city_helipads(hp_locs);}
		timer.end(); // exclude the steps below, which are dominated by model load time
		car_manager.add_parked_cars(parked_cars);
		car_manager.finalize_cars();
		car_manager.add_helicopters(hp_locs);
		ped_manager.init(city_params.num_peds); // must be after buildings are placed
	}
	cube_t get_city_bcube(unsigned city_id)               const {return road_gen.get_city_bcube(city_id);}
	cube_t get_city_bcube_at_pt(point const &pos)         const {return road_gen.get_city_bcube_at_pt(pos);}
	cube_t get_city_bcube_overlapping(cube_t const &c)    const {return road_gen.get_city_bcube_overlapping(c);}
	void get_city_bcubes(vect_cube_t &bcubes)             const {road_gen.get_city_bcubes(bcubes);}
	int check_city_contains_overlaps(cube_t const &query) const {return road_gen.check_city_contains_overlaps(query);}
	bool update_depth_if_underwater(point const &pos, float &depth) const {return road_gen.update_depth_if_underwater(pos, depth);}
	void get_all_road_bcubes(vect_cube_t &bcubes, bool connector_only) const {road_gen.get_all_road_bcubes(bcubes, connector_only);}
	void get_all_plot_zones(vect_city_zone_t &zones) {road_gen.get_all_plot_zones(zones);} // caches plot_id_offset, so non-const

	// return: 0=no coll, 1=plot coll, 2=road coll, 3=both plot and road coll; pos is in camera space
	unsigned check_city_sphere_coll(point const &pos, float radius, bool xy_only, bool exclude_bridges_and_tunnels, bool ret_first_coll, unsigned check_mask) const {
		int ret(0);
		if ((check_mask & 1) && check_plot_sphere_coll(pos, radius, xy_only)) {ret |= 1;}
		if (ret_first_coll && ret) return ret;
		if ((check_mask & 2) && road_gen.check_road_sphere_coll(pos, radius, xy_only, exclude_bridges_and_tunnels)) {ret |= 2;}
		return ret;
	}
	bool check_tline_cube_intersect_xy(cube_t const &c) const {return road_gen.check_tline_cube_intersect_xy(c);}

	void get_grass_coll_cubes(cube_t const &region, vect_cube_t &out, vect_cube_t &out_bt) const { // region is in camera space
		get_plots_in_region(region, out);
		road_gen.get_roads_in_region(region, out, out_bt);
		get_road_segs_in_region((region - get_camera_coord_space_xlate()), out); // convert region from camera space to city/building space
	}
	bool proc_city_sphere_coll(point &pos, point const &p_last, float radius, float prev_frame_zval, bool inc_cars, vector3d *cnorm) const { // pos is in camera space
		if (road_gen.proc_sphere_coll(pos, p_last, radius, prev_frame_zval, cnorm)) return 1;
		if (!inc_cars) return 0;
		return car_manager.proc_sphere_coll(pos, p_last, radius, cnorm); // Note: doesn't really work well, at least for player collisions
	}
	bool line_intersect(point const &p1, point const &p2, float &t) const { // Note: p1 and p2 are in camera space
		vector3d const xlate(get_camera_coord_space_xlate()), p1x(p1 - xlate), p2x(p2 - xlate);
		bool ret(road_gen.line_intersect(p1x, p2x, t));
		ret |= car_manager.line_intersect_cars(p1x, p2x, t);
		ret |= ped_manager.line_intersect_peds(p1x, p2x, t);
		return ret;
	}
	bool choose_pt_in_park (point const &pos, point &park_pos, rand_gen_t &rgen) const {return road_gen.choose_pt_in_park(pos, park_pos, rgen);}
	bool check_mesh_disable(point const &pos, float radius ) const {return road_gen.check_mesh_disable(pos, radius);}
	bool check_inside_city (point const &pos, float radius ) const {return road_gen.check_inside_city (pos, radius);}
	bool tile_contains_tunnel(cube_t const &bcube) const {return road_gen.tile_contains_tunnel(bcube);}
	bool cube_int_underground_obj(cube_t const &c) const {return road_gen.cube_int_underground_obj(c);}
	bool is_invalid_placement_for_cube(cube_t const &c) const {return road_gen.is_invalid_placement_for_cube(c);}
	void add_plot_cut(cube_t const &cut) {road_gen.add_plot_cut(cut);}
	void add_city_ug_elevator_entrance(ug_elev_info_t const &uge) {road_gen.add_city_ug_elevator_entrance(uge);}
	void get_ponds_in_xy_range(cube_t const &range, vect_cube_t &ponds) const {road_gen.get_ponds_in_xy_range(range, ponds);}
	void add_manhole(point const &pos, float radius) {road_gen.add_manhole(pos, radius);}

	void destroy_in_radius(point const &pos, float radius) {
		car_manager.destroy_cars_in_radius(pos, radius);
		ped_manager.destroy_peds_in_radius(pos, radius);
	}
	bool get_color_at_xy(float x, float y, colorRGBA &color) const {
		if (!city_params.enabled()) return 0;
		point const pos(point(x, y, 0.0) - get_camera_coord_space_xlate());
		int const int_ret(road_gen.get_color_at_xy(pos, color)); // check roads/plots first to determine if we need to check cars
		if (int_ret == INT_NONE) return 0; // no road/plot intersection - done
		car_manager.get_color_at_xy(pos, color, int_ret); // check cars next, but override the color
		return 1;
	}
	void next_frame(bool use_threads_2_3) { // Note: threads: 0=draw, 1=roads and cars, 2=pedestrians
		if (!city_params.enabled()) return;

		if (!use_threads_2_3 || omp_get_thread_num_3dw() == 1) { // thread 1
			road_gen.next_frame(); // update stoplights; must be before car_manager next_frame() call
			car_manager.next_frame(ped_manager, city_params.car_speed);
		}
		if (!use_threads_2_3 || omp_get_thread_num_3dw() == 2) {ped_manager.next_frame();} // thread=2
	}
	void draw(int shadow_only, int reflection_pass, int trans_op_mask, vector3d const &xlate) { // shadow_only: 0=non-shadow pass, 1=sun/moon shadow, 2=dynamic shadow
		if (player_in_basement >= 2)            return; // player is fully in the basement, not on stairs - don't draw anything
		if (player_cant_see_outside_building()) return; // player can't see outside the building (in ext basement, parking garage, attic, or windowless building)
		if (!shadow_only && !reflection_pass && (trans_op_mask & 1)) {setup_city_lights(xlate);} // setup lights on first (opaque) non-shadow pass
		bool const use_dlights(enable_lights()), is_dlight_shadows(shadow_only == 2);
		if (reflection_pass == 0) {road_gen.draw(trans_op_mask, xlate, use_dlights, (shadow_only != 0));} // roads don't cast shadows/aren't reflected in water, but stoplights cast shadows
		car_manager.draw(trans_op_mask, xlate, use_dlights, (shadow_only != 0), is_dlight_shadows);
		if  (trans_op_mask & 1) {ped_manager.draw(xlate, use_dlights, (shadow_only != 0), is_dlight_shadows);} // opaque
		if ((trans_op_mask & 1) && !shadow_only) {road_gen.draw_label();} // after drawing cars so that it's in front
		// Note: buildings are drawn through draw_buildings()
	}
	void draw_transparent(vector3d const &xlate) {road_gen.draw_transparent(xlate, enable_lights());} // drawn at the very end
	void draw_roads(int trans_op_mask, vector3d const &xlate) {road_gen.draw(trans_op_mask, xlate, enable_lights(), 0);} // shadow_only=0
	void draw_car_in_pspace(car_t &car, shader_t &s, vector3d const &xlate, bool shadow_only) {car_manager.draw_car_in_pspace(car, s, xlate, shadow_only);}
	void set_car_model_color(car_t &car) {car_manager.set_car_model_color(car);}
	void gen_and_draw_people_in_building(ped_draw_vars_t const &pdv) {ped_manager.gen_and_draw_people_in_building(pdv);}
	void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only) {ped_manager.draw_player_model(s, xlate, shadow_only);}
	bool is_player_model_female() {return ped_manager.is_player_model_female();}

	void setup_city_lights(vector3d const &xlate) {
		if (world_mode != WMODE_INF_TERRAIN) return; // TT only
		if (prev_city_lights_setup_frame == cur_display_iter) return; // already called this frame
		prev_city_lights_setup_frame = cur_display_iter;
		city_lights_custom_bcube     = 0;
		//timer_t timer("City Dlights Setup");
		float const light_radius(1.0*light_radius_scale*get_tile_smap_dist()); // distance from the camera where headlights and streetlights are drawn
		if (!begin_lights_setup(xlate, light_radius, dl_sources)) return;
		car_manager.add_car_headlights(xlate, lights_bcube);
		road_gen.add_city_lights(xlate, lights_bcube);
		if (is_night()) {add_buildings_exterior_lights(xlate, lights_bcube);} // currently building lights are only on at night
		if (flashlight_on && !camera_in_building) {add_player_flashlight(0.25);} // add player flashlight
		clamp_to_max_lights(xlate, dl_sources);
		setup_shadow_maps(dl_sources, (camera_pdu.pos - xlate), city_params.max_shadow_maps);
		finalize_lights(dl_sources);
	}
	virtual bool enable_lights() const {return (is_night(max(STREETLIGHT_ON_RAND, HEADLIGHT_ON_RAND)) || road_gen.has_tunnels() || flashlight_on);} // only have lights at night
	void next_ped_animation() {ped_manager.next_animation();}
	void get_pedestrians_in_area(cube_t const &area, int building_ix, vector<point> &pts) const {ped_manager.get_pedestrians_in_area(area, building_ix, pts);}
	void free_context() {car_manager.free_context(); ped_manager.free_context();}
	unsigned get_model_gpu_mem() const {return (ped_manager.get_model_gpu_mem() + car_manager.get_model_gpu_mem());}
}; // city_gen_t

city_gen_t city_gen;


bool parse_city_option(FILE *fp) {return city_params.read_option(fp);}
bool have_cities() {return city_params.enabled();}
// Note: this is used for parallel car/pedestrian updates
bool have_city_models() {
	return ((have_cities() && (city_params.num_cars > 0 || city_params.num_peds > 0)) ||
		    (have_buildings() && enable_building_people_ai() && global_building_params.building_people_enabled()));
}
vector2d get_road_max_len() {return vector2d(max(city_params.road_spacing, actual_max_road_seg_len.x), max(city_params.road_spacing, actual_max_road_seg_len.y));}
float get_road_max_width () {return city_params.road_width;}
float get_min_obj_spacing() {return 4.0*ped_manager_t::get_ped_radius();} // allow a ped to walk between objects (two side-by-side)

void gen_cities(float *heightmap, unsigned xsize, unsigned ysize) {
	if (!have_cities()) return; // nothing to do
	city_gen.init(heightmap, xsize, ysize); // only need to call once for any given heightmap
	city_gen.gen_cities();
	city_gen.invalidate_heightmap();
}
void gen_city_details() {city_gen.gen_details();} // called after gen_buildings()
cube_t get_city_bcube(unsigned city_id) {return city_gen.get_city_bcube(city_id);}
cube_t get_city_bcube_at_pt(point const &pos) {return city_gen.get_city_bcube_at_pt(pos);}
cube_t get_city_bcube_overlapping(cube_t const &c) {return city_gen.get_city_bcube_overlapping(c);}
void get_city_bcubes(vect_cube_t &bcubes) {city_gen.get_city_bcubes(bcubes);}
int check_city_contains_overlaps(cube_t const &query) {return city_gen.check_city_contains_overlaps(query);}
bool update_depth_if_underwater(point const &pos, float &depth) {return city_gen.update_depth_if_underwater(pos, depth);}
void get_city_road_bcubes(vect_cube_t &bcubes, bool connector_only) {city_gen.get_all_road_bcubes(bcubes, connector_only);}
void get_city_plot_zones(vect_city_zone_t &zones) {city_gen.get_all_plot_zones(zones);}
void next_city_frame(bool use_threads_2_3) {city_gen.next_frame(use_threads_2_3);}
void draw_cities(int shadow_only, int reflection_pass, int trans_op_mask, vector3d const &xlate) {city_gen.draw(shadow_only, reflection_pass, trans_op_mask, xlate);}
void draw_city_roads(int trans_op_mask, vector3d const &xlate) {city_gen.draw_roads(trans_op_mask, xlate);}
void setup_city_lights(vector3d const &xlate) {city_gen.setup_city_lights(xlate);}
void add_city_manhole(point const &pos, float radius) {city_gen.add_manhole(pos, radius);}

void draw_transparent_city_bldg_geom(int reflection_pass, vector3d const &xlate) {
	if (!reflection_pass) {city_gen.draw_transparent(xlate);}
	draw_player_building_transparent(reflection_pass, xlate);
}
void gen_and_draw_people_in_building(ped_draw_vars_t const &pdv) {city_gen.gen_and_draw_people_in_building(pdv);}
void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only) {city_gen.draw_player_model(s, xlate, shadow_only);}
bool is_player_model_female() {return city_gen.is_player_model_female();}
float get_road_height() {return ROAD_HEIGHT;} // for rain splashes

// Note: pos is in global space for these next two calls
unsigned check_city_sphere_coll(point const &pos, float radius, bool exclude_bridges_and_tunnels, bool ret_first_coll, unsigned check_mask) {
	if (!have_cities()) return 0;
	return city_gen.check_city_sphere_coll((pos + get_tt_xlate_val()), radius, 1, exclude_bridges_and_tunnels, ret_first_coll, check_mask); // apply xlate for all static objects
}
void get_city_grass_coll_cubes(cube_t const &region, vect_cube_t &out, vect_cube_t &out_bt) { // Note: region and out are in camera space
	city_gen.get_grass_coll_cubes(region, out, out_bt);
}
// primarily used for player collision (typically check_interior=1), but also used for gameplay dynamic cobjs
bool proc_city_sphere_coll(point &pos, point const &p_last, float radius, float prev_frame_zval, bool inc_cars, vector3d *cnorm, bool check_interior) {
	if (check_interior) {player_in_walkway = player_in_skyway = player_on_moving_ww = player_in_ww_elevator = player_in_tunnel = 0;} // reset for next iter if this is the player
	had_building_interior_coll = 0;
	bool ret(proc_buildings_sphere_coll(pos, p_last, radius, cnorm, check_interior));
	if (ret && had_building_interior_coll) return ret; // skip city coll if player is in a building
	ret |= city_gen.proc_city_sphere_coll(pos, p_last, radius, prev_frame_zval, inc_cars, cnorm); // check city as well
	return ret;
}
bool line_intersect_city(point const &p1, point const &p2, float &t, bool ret_any_pt) {
	unsigned hit_bix(0); // unused
	bool ret(check_buildings_line_coll(p1, p2, t, hit_bix, ret_any_pt));
	ret |= city_gen.line_intersect(p1, p2, t);
	return ret;
}
bool line_intersect_city(point const &p1, point const &p2, point &p_int) {
	float t(1.0);
	if (!line_intersect_city(p1, p2, t)) return 0;
	p_int = p1 + t*(p2 - p1);
	return 1;
}
bool check_city_tline_cube_intersect_xy(cube_t const &c) {return city_gen.check_tline_cube_intersect_xy(c);}

class model_bcube_checker_t {
	vect_cube_t model_bcubes;
	cube_t all_bcube;
	vector3d max_sz;
	bool is_valid = 0;

	struct bcube_by_y2 {
		bool operator()(cube_t const &a, cube_t const &b) const {return (a.y2() < b.y2());}
	};

	void gather_bcubes() {
		if (is_valid) return; // nothing to do
		get_all_model_bcubes(model_bcubes);
		
		for (cube_t const &c : model_bcubes) {
			all_bcube.assign_or_union_with_cube(c);
			max_sz = max_sz.max(c.get_size());
		}
		sort(model_bcubes.begin(), model_bcubes.end(), bcube_by_y2());
		is_valid = 1;
	}
public:
	bool check_sphere_coll(point const &pos, float radius, bool xy_only) {
		if (!is_valid) {
#pragma omp critical(create_model_bcubes)
			gather_bcubes(); // will be run by the first thread to get here
		}
		if (model_bcubes.empty()) return 0;
		if (!check_bcube_sphere_coll(all_bcube, pos, radius, xy_only)) return 0;
		if (model_bcubes.size() == 1) return 1; // one bcube, must be equal to all_bcube

		if (model_bcubes.size() >= 100) { // large dataset, use a binary search over x2
			auto it(lower_bound(model_bcubes.begin(), model_bcubes.end(), cube_t(0,0,0,(pos.x - radius),0,0))); // only y2 is important
			float const y2_end(pos.y + radius + max_sz.y); // we can stop when the model bcube x1 value passes this point
			
			for (auto i = it; i != model_bcubes.end(); ++i) {
				if (i->y2() > y2_end) break; // done
				if (check_bcube_sphere_coll(*i, pos, radius, xy_only)) return 1;
			}
			return 0;
		}
		return check_bcubes_sphere_coll(model_bcubes, pos, radius, xy_only); // Note: can be somewhat slow for our Puget Sound 10K museums scene
	}
};

model_bcube_checker_t model_bcube_checker;

bool check_valid_scenery_pos(point const &pos, float radius, bool is_tall) {
	point const pos_cs(pos + get_tt_xlate_val()); // convert from city/building space to camera space
	point center(pos_cs); // make a copy; not actually modified
	if (proc_buildings_sphere_coll(center, pos_cs, radius, nullptr, 0, 1)) return 0; // check_interior=0, exclude_city=1 (since we're checking plots below)
	if (world_mode != WMODE_INF_TERRAIN) return 1; // the checks below are for tiled terrain mode only

	if (have_cities()) {
		if (city_gen.check_city_sphere_coll(pos_cs, radius, 1, !is_tall, 1, 3)) return 0; // check_mask=3 to include both plots and roads
		if (city_gen.check_mesh_disable(pos_cs, radius)) return 0;
	}
	if (model_bcube_checker.check_sphere_coll((pos_cs - get_tiled_terrain_model_xlate()), radius, 1)) return 0; // xy_only=1
	return 1;
}
bool check_mesh_disable(point const &pos, float radius) { // Note: pos is in global space
	if (!have_cities()) return 0;
	return city_gen.check_mesh_disable((pos + get_tt_xlate_val()), radius); // apply xlate for all static objects
}
bool check_inside_city(point const &pos, float radius) { // Note: pos is in global space
	if (!have_cities()) return 0;
	return city_gen.check_inside_city((pos + get_tt_xlate_val()), radius); // apply xlate for all static objects
}
bool cube_int_underground_obj(cube_t const &c) {return city_gen.cube_int_underground_obj(c);} // Note: cube is in global space
bool is_invalid_city_placement_for_cube(cube_t const &c) {return city_gen.is_invalid_placement_for_cube(c);}
void add_city_plot_cut(cube_t const &cut) {city_gen.add_plot_cut(cut);}
void add_city_ug_elevator_entrance(ug_elev_info_t const &uge) {city_gen.add_city_ug_elevator_entrance(uge);}
void get_ponds_in_xy_range(cube_t const &range, vect_cube_t &ponds) {city_gen.get_ponds_in_xy_range(range, ponds);} // Note: range is in global space
bool choose_pt_in_city_park(point const &pos, point &park_pos, rand_gen_t &rgen) {return city_gen.choose_pt_in_park(pos, park_pos, rgen);}
bool tile_contains_tunnel(cube_t const &bcube) {return city_gen.tile_contains_tunnel(bcube);}
void destroy_city_in_radius(point const &pos, float radius) {city_gen.destroy_in_radius(pos, radius);}
bool get_city_color_at_xy(float x, float y, colorRGBA &color) {return city_gen.get_color_at_xy(x, y, color);}
cube_t get_city_lights_bcube() {return city_gen.get_lights_bcube();}
unsigned get_city_model_gpu_mem() {return city_gen.get_model_gpu_mem();}
void next_pedestrian_animation() {city_gen.next_ped_animation();}
void get_pedestrians_in_area(cube_t const &area, int building_ix, vector<point> &pts) {city_gen.get_pedestrians_in_area(area, building_ix, pts);}
void free_city_context() {city_gen.free_context();}
bool has_city_trees() {return (city_params.max_trees_per_plot > 0);}
vector3d get_nom_car_size() {return city_params.get_nom_car_size();}
void draw_car_in_pspace(car_t &car, shader_t &s, vector3d const &xlate, bool shadow_only) {city_gen.draw_car_in_pspace(car, s, xlate, shadow_only);}
void set_car_model_color(car_t &car) {city_gen.set_car_model_color(car);}

