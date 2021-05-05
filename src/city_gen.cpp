// 3D World - City Generation
// by Frank Gennari
// 2/10/18

#include "city.h"
#include "mesh.h"
#include "heightmap.h"
#include "lightmap.h"
#include "buildings.h"
#include "tree_3dw.h"
#include "profiler.h"
#include <cfloat> // for FLT_MAX

using std::string;

bool const CHECK_HEIGHT_BORDER_ONLY = 1; // choose building site to minimize edge discontinuity rather than amount of land that needs to be modified
float const OUTSIDE_TERRAIN_HEIGHT  = 0.0;
float const CAR_LANE_OFFSET         = 0.15; // in units of road width
float const CITY_LIGHT_FALLOFF      = 0.2;


float city_dlight_pcf_offset_scale(1.0);
city_params_t city_params;
point pre_smap_player_pos(all_zeros);

extern bool enable_dlight_shadows, dl_smap_enabled, flashlight_on, camera_in_building, have_indir_smoke_tex, disable_city_shadow_maps;
extern int rand_gen_index, display_mode, animate2, draw_model;
extern unsigned shadow_map_sz, cur_display_iter, max_unique_trees;
extern float water_plane_z, shadow_map_pcf_offset, cobj_z_bias, fticks;
extern vector<light_source> dl_sources;
extern tree_placer_t tree_placer;
extern object_model_loader_t building_obj_model_loader;


void add_dynamic_lights_city(cube_t const &scene_bcube, float &dlight_add_thresh);
void disable_shadow_maps(shader_t &s);
vector3d get_tt_xlate_val();


template<typename S, typename T> void get_all_bcubes(vector<T> const &v, S &bcubes) {
	for (auto i = v.begin(); i != v.end(); ++i) {bcubes.push_back(*i);}
}

float smooth_interp(float a, float b, float mix) {
	mix = mix * mix * (3.0 - 2.0 * mix); // cubic Hermite interoplation (smoothstep)
	return mix*a + (1.0 - mix)*b;
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
		s.add_uniform_float("z_bias", pcf_scale*cobj_z_bias); // I guess pcf_scale is really some sort of light size scale and should apply to the z-bias as well
		s.add_uniform_float("pcf_offset", 8.0*pcf_scale*shadow_map_pcf_offset);
		s.add_uniform_float("dlight_pcf_offset", 0.0005*pcf_scale*city_dlight_pcf_offset_scale);
	}
}

// use_smap: 0=no, 1=sun/moon + dynamic lights; enable in shader and set shadow map uniforms, 2=dynamic lights only; disable in shader but set shadow map uniforms
void city_shader_setup(shader_t &s, cube_t const &lights_bcube, bool use_dlights, int use_smap, int use_bmap,
	float min_alpha, bool force_tsl, float pcf_scale, bool use_texgen, bool indir_lighting)
{
	use_dlights &= (!lights_bcube.is_zero_area() && !dl_sources.empty());
	have_indir_smoke_tex = indir_lighting; // assume someone is going to set the indir texture in this case
	if (indir_lighting) {s.set_prefix("#define USE_ALT_SCENE_BOUNDS", 1);} // FS; need to use different scene_llc_scale for dynamic lighting vs. building indir lighting
	//if ((display_mode & 0x10) && use_dlights) {s.set_prefix("#define LINEAR_DLIGHT_ATTEN", 1);} // FS
	// Note: here use_texgen mode 5 is used as a hack so that the shader still has binding points for tex coords (can't optimize it out)
	// and we can share the same VAO between texgen and texcoords modes without having to worry about which mode we were in when the VAO was created
	setup_smoke_shaders(s, min_alpha, (use_texgen ? 5 : 0), 0, indir_lighting, 1, use_dlights,
		0, 0, ((use_smap == 1) ? 2 : 0), use_bmap, 0, use_dlights, force_tsl, 0.0, 0.0, 0, 0, 1); // use_spec_map=0, is_outside=1
	set_city_lighting_shader_opts(s, lights_bcube, use_dlights, (use_smap != 0), pcf_scale);
	if (use_texgen) {s.add_uniform_float("tc_texgen_mix", 0.0);} // always uses texgen in this mode
}

void enable_animations_for_shader(shader_t &s) {s.add_property("animation_shader", "pedestrian_animation.part+");}

void draw_state_t::begin_tile(point const &pos, bool will_emit_now, bool ensure_active) {
	if (ensure_active) {ensure_shader_active();} // needed for use_smap=0 case
	emit_now = (use_smap && try_bind_tile_smap_at_point((pos + xlate), s));
	if (will_emit_now && !emit_now) {disable_shadow_maps(s);} // not using shadow maps or second (non-shadow map) pass - disable shadow maps
}
void draw_state_t::pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_, bool always_setup_shader) {
	xlate       = xlate_;
	shadow_only = shadow_only_;
	use_dlights = (use_dlights_ && !shadow_only);
	use_smap    = (shadow_map_enabled() && !shadow_only && !disable_city_shadow_maps);
	if (!use_smap && !always_setup_shader) return;
	if (shadow_only) {s.begin_simple_textured_shader();}
	else {
		cube_t const &lights_bcube(use_building_lights ? get_building_lights_bcube() : get_city_lights_bcube());
		city_shader_setup(s, lights_bcube, use_dlights, use_smap, (use_bmap && !shadow_only), 0.01, 0, 0.5);
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
void draw_state_t::ensure_shader_active() {
	if (s.is_setup()) return; // already active
	if (shadow_only) {s.begin_color_only_shader();}
	else {city_shader_setup(s, get_city_lights_bcube(), use_dlights, 0, use_bmap, 0.01, 0, 0.5);} // no smap
}
void draw_state_t::draw_and_clear_light_flares() {
	if (light_psd.empty()) return; // no lights to draw
	enable_blend();
	set_additive_blend_mode();
	light_psd.draw_and_clear(BLUR_TEX, 0.0, 0, 1, 0.005); // use geometry shader for unlimited point size
	set_std_blend_mode();
	disable_blend();
}
bool draw_state_t::check_cube_visible(cube_t const &bc, float dist_scale, bool shadow_only) const {
	if (!camera_pdu.valid) return 1;
	cube_t const bcx(bc + xlate);

	if (dist_scale > 0.0) {
		float const dmax(shadow_only ? camera_pdu.far_ : dist_scale*get_draw_tile_dist());
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
	vector3d const cview_dir((camera_pdu.pos - xlate) - center);
	float const sign(invert_normals ? -1.0 : 1.0);
	vector3d const top_n  (cross_product((p[2] - p[1]), (p[0] - p[1]))*sign); // Note: normalization not needed
	vector3d const front_n(cross_product((p[5] - p[1]), (p[0] - p[1]))*sign);
	vector3d const right_n(cross_product((p[6] - p[2]), (p[1] - p[2]))*sign);
	tex_range_t tr_top, tr_front, tr_right;

	if (tscale > 0.0) { // compute texture s/t parameters from cube side lengths to get a 1:1 AR
		tr_top   = tex_range_t(0.0, 0.0, tscale*(p[0] - p[1]).mag(), tscale*(p[2] - p[1]).mag());
		tr_front = tex_range_t(0.0, 0.0, tscale*(p[0] - p[1]).mag(), tscale*(p[5] - p[1]).mag());
		tr_right = tex_range_t(0.0, 0.0, tscale*(p[1] - p[2]).mag(), tscale*(p[6] - p[2]).mag());
	}
	if (!(skip_dims & 4)) { // Z
		if (dot_product(cview_dir, top_n) > 0) {qbd.add_quad_pts(p+4, cw,  top_n, tr_top);} // top
		else if (!skip_bottom)                 {qbd.add_quad_pts(p+0, cw, -top_n, tr_top);} // bottom - not always drawn
	}
	if (!(skip_dims & 1)) { // X
		if (dot_product(cview_dir, front_n) > 0) {point const pts[4] = {p[0], p[1], p[5], p[4]}; qbd.add_quad_pts(pts, cw,  front_n, tr_front);} // front
		else                                     {point const pts[4] = {p[2], p[3], p[7], p[6]}; qbd.add_quad_pts(pts, cw, -front_n, tr_front);} // back
	}
	if (!(skip_dims & 2)) { // Y
		if (dot_product(cview_dir, right_n) > 0) {point const pts[4] = {p[1], p[2], p[6], p[5]}; qbd.add_quad_pts(pts, cw,  right_n, tr_right);} // right
		else                                     {point const pts[4] = {p[3], p[0], p[4], p[7]}; qbd.add_quad_pts(pts, cw, -right_n, tr_right);} // left
	}
}
void draw_state_t::draw_cube(quad_batch_draw &qbd, cube_t const &c, color_wrapper const &cw, bool skip_bottom, float tscale, unsigned skip_dims) const {
	point pts[8];
	set_cube_pts(c, 0, 0, pts);
	draw_cube(qbd, cw, c.get_cube_center(), pts, skip_bottom, 0, tscale, skip_dims);
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

template<typename T> void get_bcubes_sphere_coll(vector<T> const &bcubes, vect_cube_t &out, point const &sc, float radius, bool xy_only, vector3d const &xlate) {
	for (auto i = bcubes.begin(); i != bcubes.end(); ++i) {
		if (check_bcube_sphere_coll(get_bcube(*i), sc, radius, xy_only)) {out.push_back(get_bcube(*i) + xlate);}
	}
}
template<typename T> cube_t calc_cubes_bcube(vector<T> const &cubes) {
	if (cubes.empty()) return cube_t(all_zeros);
	cube_t bcube(cubes.front()); // first cube
	for (auto r = cubes.begin()+1; r != cubes.end(); ++r) {bcube.union_with_cube(*r);} // skip first cube
	return bcube;
}

template<typename T> bool has_bcube_int_xy(cube_t const &bcube, vector<T> const &bcubes, float pad_dist) { // T must derive from cube_t
	cube_t tc(bcube);
	tc.expand_by_xy(pad_dist);

	for (auto c = bcubes.begin(); c != bcubes.end(); ++c) {
		if (c->intersects_xy(tc)) return 1; // intersection
	}
	return 0;
}
// explicit instantiations
template bool has_bcube_int_xy(cube_t const &bcube, vector<cube_t    > const &bcubes, float pad_dist);
template bool has_bcube_int_xy(cube_t const &bcube, vector<elevator_t> const &bcubes, float pad_dist);

point rand_xy_pt_in_cube(cube_t const &c, float radius, rand_gen_t &rgen) {
	return point(rgen.rand_uniform(c.x1()+radius, c.x2()-radius), rgen.rand_uniform(c.y1()+radius, c.y2()-radius), c.z1());
}


class heightmap_query_t {
protected:
	float *heightmap;
	unsigned xsize, ysize;

public:
	flatten_op_t last_flatten_op;

	heightmap_query_t() : heightmap(nullptr), xsize(0), ysize(0) {}
	heightmap_query_t(float *hmap, unsigned xsize_, unsigned ysize_) : heightmap(hmap), xsize(xsize_), ysize(ysize_) {}
	float get_x_value(int x) const {return get_xval(x - int(xsize)/2);} // convert from center to LLC
	float get_y_value(int y) const {return get_yval(y - int(ysize)/2);}
	int get_x_pos(float x) const {return (get_xpos(x) + int(xsize)/2);}
	int get_y_pos(float y) const {return (get_ypos(y) + int(ysize)/2);}
	float  get_height(unsigned x, unsigned y) const {return heightmap[y*xsize + x];} // Note: not bounds checked
	float &get_height(unsigned x, unsigned y)       {return heightmap[y*xsize + x];} // Note: not bounds checked
	bool is_normalized_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <  x2 && y1 <  y2 && x2 <= xsize && y2 <= ysize);}
	bool is_valid_region     (unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <= x2 && y1 <= y2 && x2 <= xsize && y2 <= ysize);}
	bool is_inside_terrain(int x, int y) const {return (x >= 0 && y >= 0 && x < (int)xsize && y < (int)ysize);}
	cube_t get_full_hmap_bcube() const {return get_cube_for_bounds(0, 0, xsize, ysize, 0.0);}
	
	cube_t get_cube_for_bounds(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float elevation) const {
		cube_t c;
		c.x1() = get_x_value(x1);
		c.x2() = get_x_value(x2);
		c.y1() = get_y_value(y1);
		c.y2() = get_y_value(y2);
		c.z1() = c.z2() = elevation;
		return c;
	}
	float get_height_at(float xval, float yval) const {
		int const x(get_x_pos(xval)), y(get_y_pos(yval));
		return (is_inside_terrain(x, y) ? get_height(x, y) : OUTSIDE_TERRAIN_HEIGHT);
	}
	bool any_underwater(unsigned x1, unsigned y1, unsigned x2, unsigned y2, bool check_border=0) const {
		min_eq(x2, xsize); min_eq(y2, ysize); // clamp upper bound
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (check_border && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				if (get_height(x, y) < water_plane_z) return 1;
			}
		}
		return 0;
	}
	void get_segment_end_pts(road_t const &r, unsigned six, unsigned eix, point &ps, point &pe) const {
		float const sv(r.dim ? get_y_value(six) : get_x_value(six)), ev(r.dim ? get_y_value(eix) : get_x_value(eix)); // start/end pos of bridge in road dim
		float const z1(r.get_start_z()), z2(r.get_end_z()), dz(z2 - z1), len(r.get_length()), v0(r.dim ? r.y1() : r.x1());
		ps[ r.dim] = sv;
		pe[ r.dim] = ev;
		ps[!r.dim] = pe[!r.dim] = r.get_center_dim(!r.dim);
		ps.z = z1 + dz*CLIP_TO_01((sv - v0)/len);
		pe.z = z1 + dz*CLIP_TO_01((ev - v0)/len);
	}
	void flatten_region_to(cube_t const c, unsigned slope_width, bool decrease_only=0) {
		flatten_region_to(get_x_pos(c.x1()), get_y_pos(c.y1()), get_x_pos(c.x2()), get_y_pos(c.y2()), slope_width, (c.z1() - ROAD_HEIGHT), decrease_only);
	}
	void flatten_region_to(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float elevation, bool decrease_only=0) {
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = max((int)y1-(int)slope_width, 0); y < min(y2+slope_width, ysize); ++y) {
			for (unsigned x = max((int)x1-(int)slope_width, 0); x < min(x2+slope_width, xsize); ++x) {
				float &h(get_height(x, y));
				if (decrease_only && h < elevation) continue; // don't increase

				if (slope_width > 0) {
					float const dx(max(0, max(((int)x1 - (int)x), ((int)x - (int)x2 + 1))));
					float const dy(max(0, max(((int)y1 - (int)y), ((int)y - (int)y2 + 1))));
					h = smooth_interp(h, elevation, min(1.0f, sqrt(dx*dx + dy*dy)/slope_width));
				} else {h = elevation;}
			} // for x
		} // for y
	}
	float flatten_sloped_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float z1, float z2, bool dim, unsigned border,
		unsigned skip_six=0, unsigned skip_eix=0, bool stats_only=0, bool decrease_only=0, bridge_t *bridge=nullptr, tunnel_t *tunnel=nullptr)
	{
		if (!stats_only) {last_flatten_op = flatten_op_t(x1, y1, x2, y2, z1, z2, dim, border);} // cache for later replay
		assert(is_valid_region(x1, y1, x2, y2));
		if (x1 == x2 || y1 == y2) return 0.0; // zero area
		float const run_len(dim ? (y2 - y1) : (x2 - x1)), denom(1.0f/max(run_len, 1.0f)), dz(z2 - z1), border_inv(1.0/border);
		int const pad(border + 1U); // pad an extra 1 texel to handle roads misaligned with the texture
		unsigned px1(x1), py1(y1), px2(x2), py2(y2), six(dim ? ysize : xsize), eix(0);
		float tot_dz(0.0), seg_min_dh(0.0);
		float const bridge_cost(0.0), bridge_dist_cost(0.0), tunnel_cost(0.0), tunnel_dist_cost(0.0); // Note: currently set to zero, but could be used
		unsigned const min_bridge_len(12), min_tunnel_len(12); // in mesh texels
		
		if (dim) {
			px1 = max((int)x1-pad, 0);
			px2 = min(x2+pad, xsize);
			py1 = max((int)y1-1, 0); // pad by 1 in road dim as well to blend with edge of city
			py2 = min(y2+1, ysize);
		}
		else {
			py1 = max((int)y1-pad, 0);
			py2 = min(y2+pad, ysize);
			px1 = max((int)x1-1, 0);
			px2 = min(x2+1, xsize);
		}
		if (!stats_only && !decrease_only && bridge != nullptr && fabs(bridge->get_slope_val()) < 0.1) { // determine if we should add a bridge here
			float added(0.0), removed(0.0), total(0.0);
			bool end_bridge(0);

			for (unsigned y = y1; y < y2; ++y) { // Note: not padded
				for (unsigned x = x1; x < x2; ++x) {
					float const t(((dim ? int(y - y1) : int(x - x1)) + ((dz < 0.0) ? 1 : -1))*denom); // bias toward the lower zval
					float const road_z(z1 + dz*t - ROAD_HEIGHT), h(get_height(x, y));
					
					if (road_z > h) {
						added += (road_z - h);

						if (!end_bridge && road_z > h + 1.0*city_params.road_width) { // higher than terrain by a significant amount
							min_eq(six, (dim ? y : x));
							max_eq(eix, (dim ? y : x));
						}
					}
					else {
						removed += (h - road_z);
						if (eix > 0) {end_bridge = 1;} // done with bridge - don't create bridge past high point
					}
					total += 1.0;
				} // for x
			} // for y
			max_eq(six, (dim ? y1+border : x1+border)); // keep away from segment end points (especially connector road jogs)
			min_eq(eix, (dim ? y2-border : x2-border));

			if (eix > six+min_bridge_len && added > 1.5*city_params.road_width*total && added > 2.0*removed) {
				point ps, pe;
				get_segment_end_pts(bridge->src_road, six, eix, ps, pe);
				bridge->d[dim][0] = ps[dim];
				bridge->d[dim][1] = pe[dim];
				bridge->z1() = min(ps.z, pe.z);
				bridge->z2() = max(ps.z, pe.z);
				bridge->make_bridge = 1;
				skip_six = six; skip_eix = eix; // mark so that mesh height isn't updated in this region
				tot_dz += bridge_cost + bridge_dist_cost*bridge->get_length();
			}
		} // end bridge logic
		if (!stats_only && tunnel != nullptr && skip_eix == 0 && fabs(tunnel->get_slope_val()) < 0.2) { // determine if we should add a tunnel here
			float const radius(1.0*city_params.road_width), min_height((1.0 + TUNNEL_WALL_THICK)*radius);
			float added(0.0), removed(0.0), total(0.0);
			bool end_tunnel(0);

			for (unsigned y = y1; y < y2; ++y) { // Note: not padded
				for (unsigned x = x1; x < x2; ++x) {
					float const t(((dim ? int(y - y1) : int(x - x1)) + ((dz < 0.0) ? 1 : -1))*denom); // bias toward the lower zval
					float const road_z(z1 + dz*t), h(get_height(x, y));

					if (road_z < h) {
						removed += (h - road_z);

						if (!end_tunnel && road_z + min_height < h) { // below terrain by a significant amount
							min_eq(six, (dim ? y : x));
							max_eq(eix, (dim ? y : x));
						}
					}
					else {
						added += (road_z - h);
						if (eix > 0) {end_tunnel = 1;} // done with tunnel - don't create tunnel past low point
					}
					total += 1.0;
				} // for x
			} // for y
			max_eq(six, (dim ? y1+border : x1+border)); // keep away from segment end points (especially connector road jogs)
			min_eq(eix, (dim ? y2-border : x2-border));

			if (eix > six+min_tunnel_len && removed > 1.0*city_params.road_width*total && removed > 2.0*added) {
				point ps, pe;
				get_segment_end_pts(*tunnel, six, eix, ps, pe);
				float const len(fabs(ps[dim] - pe[dim]));
				
				if (len > 4.0*radius) { // don't make the tunnel too short
					tunnel->init(ps, pe, radius, 0.5*radius, dim);
					unsigned const tunnel_border((unsigned)ceil(radius/(dim ? DY_VAL : DX_VAL)));
					seg_min_dh = tunnel->height + TUNNEL_WALL_THICK*radius; // add another wall thickness to account for sloped terrain minima
					skip_six = six; skip_eix = eix; // mark so that mesh height isn't updated in this region
					tot_dz += tunnel_cost + tunnel_dist_cost*tunnel->get_length();
					int const rwidth(ceil(city_params.road_width/(dim ? DX_VAL : DY_VAL)));

					for (int dxy = -rwidth; dxy <= rwidth; ++dxy) { // shifts in !dim
						unsigned qpt[2] = {(x1 + x2)/2, (y1 + y2)/2}; // start at center
						qpt[!dim] += dxy;

						for (unsigned n = 0; n < tunnel_border; ++n) { // take several samples and find the peak mesh height for the tunnel facades
							qpt[dim] = six + n;
							max_eq(tunnel->facade_height[0], (get_height(qpt[0], qpt[1]) - ps.z - radius)); // effectively adds an additional wall height (= tunnel->height - radius)
							qpt[dim] = eix - n;
							max_eq(tunnel->facade_height[1], (get_height(qpt[0], qpt[1]) - pe.z - radius));
						} // for n
					} // for dxy
				}
			}
		} // end tunnel logic
		if (!stats_only && skip_six < skip_eix) {last_flatten_op.skip_six = skip_six; last_flatten_op.skip_eix = skip_eix;} // clip to a partial range

		for (unsigned y = py1; y < py2; ++y) {
			for (unsigned x = px1; x < px2; ++x) {
				float const t(((dim ? int(y - y1) : int(x - x1)) + ((dz < 0.0) ? 1 : -1))*denom); // bias toward the lower zval
				float const road_z(z1 + dz*t - ROAD_HEIGHT);
				float &h(get_height(x, y));
				if (decrease_only && h < road_z) continue; // don't increase
				float new_h;
				unsigned dist(0);

				if (border > 0) {
					dist = (dim ? max(0, max(((int)x1 - (int)x - 1), ((int)x - (int)x2))) : max(0, max(((int)y1 - (int)y - 1), ((int)y - (int)y2))));
					new_h = smooth_interp(h, road_z, dist*border_inv);
				} else {new_h = road_z;}
				tot_dz += fabs(h - new_h);
				if (stats_only) continue; // no height update
				unsigned const dv(dim ? y : x);

				if (dv > skip_six && dv < skip_eix) { // don't modify mesh height at bridges or tunnels, but still count it toward the cost
					if (seg_min_dh > 0.0) { // clamp to roof of tunnel (Note: doesn't count toward tot_dz)
						float const zmin(road_z + seg_min_dh);
						if (h < zmin) {h = smooth_interp(h, zmin, dist*border_inv);}
					}
				}
				else {h = new_h;} // apply the height change
			} // for x
		} // for y
		return tot_dz;
	}
	float flatten_for_road(road_t const &road, unsigned border, bool stats_only=0, bool decrease_only=0, bridge_t *bridge=nullptr, tunnel_t *tunnel=nullptr) {
		float const z_adj(road.get_z_adj());
		unsigned const rx1(get_x_pos(road.x1())), ry1(get_y_pos(road.y1())), rx2(get_x_pos(road.x2())), ry2(get_y_pos(road.y2()));
		return flatten_sloped_region(rx1, ry1, rx2, ry2, road.d[2][road.slope]-z_adj, road.d[2][!road.slope]-z_adj, road.dim, border, 0, 0, stats_only, decrease_only, bridge, tunnel);
	}
}; // heightmap_query_t


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
	bool check_plot_sphere_coll(point const &pos, float radius, bool xy_only=1) const {
		if (plots.empty()) return 0;
		point const query_pos(pos - get_camera_coord_space_xlate());
		if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
		return check_bcubes_sphere_coll(plots, query_pos, radius, xy_only);
	}
	void get_plots_sphere_coll(point const &pos, float radius, bool xy_only, vect_cube_t &out) const {
		if (plots.empty()) return;
		vector3d const xlate(get_camera_coord_space_xlate());
		point const query_pos(pos - xlate);
		if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return;
		get_bcubes_sphere_coll(plots, out, query_pos, radius, xy_only, xlate);
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


class city_road_gen_t : public road_gen_base_t {

	struct city_obj_t : public sphere_t {
		cube_t bcube;
		city_obj_t() {}
		city_obj_t(point const &pos_, float radius_) : sphere_t(pos_, radius_) {}
		bool operator<(city_obj_t const &b) const {return (bcube.x1() < b.bcube.x1());} // sort by bcube x1
		static void post_draw(draw_state_t &dstate, bool shadow_only) {}
	};

	struct bench_t : public city_obj_t {
		bool dim, dir;

		bench_t() : dim(0), dir(0) {}

		void calc_bcube() {
			bcube.set_from_point(pos);
			bcube.expand_by(vector3d((dim ? 0.32 : 1.0), (dim ? 1.0 : 0.32), 0.0)*radius);
			bcube.z2() += 0.85*radius; // set bench height
		}
		static void pre_draw(draw_state_t &dstate, bool shadow_only) {
			if (!shadow_only) {select_texture(FENCE_TEX);} // normal map?
		}
		void draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const {
			if (!dstate.check_cube_visible(bcube, dist_scale, shadow_only)) return;

			cube_t cubes[] = { // Note: taken from mapx/bench.txt
				cube_t(-0.4, 0.0,  -5.0,   5.0,   1.6, 5.0), // back (straight)
				cube_t( 0.0, 4.0,  -5.35,  5.35,  1.6, 2.0), // seat
				cube_t( 0.3, 1.3,  -5.3,  -4.7,   0.0, 1.6), // legs
				cube_t( 2.7, 3.7,  -5.3,  -4.7,   0.0, 1.6),
				cube_t( 0.3, 1.3,   4.7,   5.3,   0.0, 1.6),
				cube_t( 2.7, 3.7,   4.7,   5.3,   0.0, 1.6),
				cube_t(-0.5, 3.8,  -5.4,  -4.5,   3.0, 3.2), // arms
				cube_t(-0.5, 3.8,   4.5,   5.4,   3.0, 3.2),
				cube_t( 0.8, 1.2,  -5.1,  -4.9,   2.0, 3.0), // arm supports
				cube_t( 2.8, 3.2,  -5.1,  -4.9,   2.0, 3.0),
				cube_t( 0.8, 1.2,   4.9,   5.1,   2.0, 3.0),
				cube_t( 2.8, 3.2,   4.9,   5.1,   2.0, 3.0),
			};
			point const center(pos + dstate.xlate);
			float const dist_val(shadow_only ? 0.0 : p2p_dist(camera_pdu.pos, center)/get_draw_tile_dist());
			cube_t bc; // bench bbox

			for (unsigned i = 0; i < 12; ++i) { // back still contributes to bbox
				if (dir)  {swap(cubes[i].d[0][0], cubes[i].d[0][1]); cubes[i].d[0][0] *= -1.0; cubes[i].d[0][1] *= -1.0;}
				if (!dim) {swap(cubes[i].d[0][0], cubes[i].d[1][0]); swap(cubes[i].d[0][1], cubes[i].d[1][1]);}
				if (i == 0) {bc = cubes[i];} else {bc.union_with_cube(cubes[i]);}
			}
			point const c1(bcube.get_cube_center()), c2(bc.get_cube_center());
			vector3d const scale(bcube.dx()/bc.dx(), bcube.dy()/bc.dy(), bcube.dz()/bc.dz()); // scale to fit to target cube
			color_wrapper const cw(WHITE);
			unsigned const num(shadow_only ? 6U : max(1U, min(6U, unsigned(0.2/dist_val)))); // simple distance-based LOD, in pairs
			for (unsigned i = 1; i < 2*num; ++i) {dstate.draw_cube(qbd, ((cubes[i] - c2)*scale + c1), cw, 1);} // skip back
			point pts[4] = {point(-1.0, -5.0, 5.0), point(-1.0, 5.0, 5.0), point(0.2, 5.0, 1.6), point(0.2, -5.0, 1.6)}; // Note: back not drawn
			point f[4], b[4];

			for (unsigned i = 0; i < 4; ++i) {
				if (dir)  {pts[i].x *= -1.0;}
				if (!dim) {swap(pts[i].x, pts[i].y);}
				pts[i] = ((pts[i] - c2)*scale + c1);
			}
			vector3d const normal(get_poly_norm(pts, 1)), delta((0.2*scale.x)*normal); // thickness = 0.4
			UNROLL_4X(f[i_] = pts[i_] + delta;);
			qbd.add_quad_pts(f, WHITE,  normal);
			UNROLL_4X(b[i_] = pts[i_] - delta;);
			qbd.add_quad_pts(b, WHITE, -normal);

			for (unsigned i = 0; i < 4; ++i) { // draw sides
				unsigned const j((i+1)&3); // next i
				point const s[4] = {f[i], b[i], b[j], f[j]};
				qbd.add_quad_pts(s, WHITE, get_poly_norm(s, 1));
			}
		}
		bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
			return sphere_cube_int_update_pos(pos_, radius_, (bcube + xlate), p_last, 1, 0, cnorm);
		}
	};

	struct tree_planter_t : public city_obj_t {
		tree_planter_t(point const &pos_, float radius_, float height) : city_obj_t(pos_, radius_) {
			bcube.set_from_point(pos);
			bcube.expand_by_xy(radius);
			bcube.z2() += height;
		}
		static void pre_draw(draw_state_t &dstate, bool shadow_only) {
			if (!shadow_only) {select_texture((dstate.pass_ix == 0) ? (int)DIRT_TEX : get_texture_by_name("roads/sidewalk.jpg"));}
		}
		void draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const {
			if (!dstate.check_cube_visible(bcube, dist_scale, shadow_only)) return;
			color_wrapper const cw(LT_GRAY);
			cube_t dirt(bcube);
			dirt.expand_by_xy(-0.1*dirt.get_size()); // shrink 10% on all XY sides

			if (dstate.pass_ix == 0) { // draw dirt
				dirt.z2() -= 0.25*bcube.dz(); // move down 25%
				dstate.draw_cube(qbd, dirt, cw, 1, 0.0, 3); // top only (skip X, Y, and bottom)
			}
			else { // draw stone
				cube_t walls[4] = {bcube, bcube, bcube, bcube}; // -X, +X, -Y, +Y
				walls[0].x2() = walls[2].x1() = walls[3].x1() = dirt.x1();
				walls[1].x1() = walls[2].x2() = walls[3].x2() = dirt.x2();
				walls[2].y2() = dirt.y1();
				walls[3].y1() = dirt.y2();
				float const tscale(40.0);
			
				for (unsigned d = 0; d < 2; ++d) {
					dstate.draw_cube(qbd, walls[d  ], cw, 1, tscale, 0); // X
					dstate.draw_cube(qbd, walls[d+2], cw, 1, tscale, 1); // Y, skip X dims
				}
			}
		}
	};

	struct fire_hydrant_t : public city_obj_t {
		float cylin_radius;
		vector3d orient;

		fire_hydrant_t(point const &pos_, float radius_, float height, vector3d const &orient_) : city_obj_t(pos_, radius_), cylin_radius(radius), orient(orient_) {
			bcube.set_from_sphere(*this);
			set_cube_zvals(bcube, pos.z, pos.z+height);
			pos.z += 0.5*height; // pos is bottom center point, make it the center
			max_eq(radius, 0.5f*height); // use a more accurate bounding sphere; Note: no cube root of (r*r + r*r + h*h)
		}
		static void pre_draw(draw_state_t &dstate, bool shadow_only) {
			if (!shadow_only) {select_texture(WHITE_TEX);}
			if (!shadow_only) {dstate.s.set_cur_color(colorRGBA(1.0, 0.75, 0.0));}
		}
		static void post_draw(draw_state_t &dstate, bool shadow_only) {
			if (!shadow_only) {dstate.s.set_cur_color(WHITE);} // restore to default color
		}
		void draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const { // Note: qbd is unused
			if (!dstate.check_cube_visible(bcube, dist_scale, shadow_only)) return;

			if (!shadow_only && building_obj_model_loader.is_model_valid(OBJ_MODEL_FHYDRANT)) {
				building_obj_model_loader.draw_model(dstate.s, pos, bcube, orient, WHITE, dstate.xlate, OBJ_MODEL_FHYDRANT, shadow_only);
			}
			else { // draw as a simple cylinder, untextured, top end only
				draw_fast_cylinder(point(pos.x, pos.y, bcube.z1()), point(pos.x, pos.y, bcube.z2()), 0.8*cylin_radius, 0.8*cylin_radius, (shadow_only ? 12 : N_CYL_SIDES), 0, 4);
			}
		}
		bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
			point const pos2(pos + xlate);
			float const r_sum(cylin_radius + radius_);
			if (!dist_less_than(pos_, pos2, r_sum)) return 0; // use sphere/vert cylinder instead?
			// since this is a cylinder, and we're not supposed to stand on top of it, assume collision normal is in the XY plane
			vector3d const coll_norm(vector3d((pos_.x - pos2.x), (pos_.y - pos2.y), 0.0).get_norm());
			pos_ += coll_norm*(r_sum - p2p_dist(pos_, pos2)); // move away from pos2
			if (cnorm) {*cnorm = coll_norm;}
			return 1;
		}
	};

	class city_obj_placer_t {
	public: // road network needs access to parking lots for drawing
		vector<parking_lot_t> parking_lots;
	private:
		vector<bench_t> benches;
		vector<tree_planter_t> planters;
		vector<fire_hydrant_t> fire_hydrants;
		vector<cube_with_ix_t> parking_lot_groups, bench_groups, planter_groups, fire_hydrant_groups; // index is last object in group
		quad_batch_draw qbd;
		unsigned num_spaces, filled_spaces;

		bool gen_parking_lots_for_plot(cube_t plot, vector<car_t> &cars, unsigned city_id, unsigned plot_ix, vect_cube_t &bcubes, vect_cube_t &colliders, rand_gen_t &rgen, bool is_new_tile) {
			vector3d const nom_car_size(city_params.get_nom_car_size()); // {length, width, height}
			float const space_width(PARK_SPACE_WIDTH *nom_car_size.y); // add 50% extra space between cars
			float const space_len  (PARK_SPACE_LENGTH*nom_car_size.x); // space for car + gap for cars to drive through
			float const pad_dist   (max(1.0f*nom_car_size.x, get_min_obj_spacing())); // one car length or min building spacing
			plot.expand_by_xy(-pad_dist);
			if (bcubes.empty()) return 0; // shouldn't happen, unless buildings are disabled; skip to avoid perf problems with an entire plot of parking lot
			unsigned const first_corner(rgen.rand()&3); // 0-3
			bool const car_dim(rgen.rand() & 1); // 0=cars face in X; 1=cars face in Y
			bool const car_dir(rgen.rand() & 1);
			float const xsz(car_dim ? space_width : space_len), ysz(car_dim ? space_len : space_width);
			bool has_parking(0);
			//cout << "max_row_sz: " << floor(plot.get_size()[!car_dim]/space_width) << ", max_num_rows: " << floor(plot.get_size()[car_dim]/space_len) << endl;
			car_t car;
			car.park();
			car.cur_city = city_id;
			car.cur_road = plot_ix; // store plot_ix in road field
			car.cur_road_type = TYPE_PLOT;

			for (unsigned c = 0; c < 4; ++c) { // generate 0-4 parking lots per plot, starting at the corners, in random order
				unsigned const cix((first_corner + c) & 3), xdir(cix & 1), ydir(cix >> 1), wdir(car_dim ? xdir : ydir), rdir(car_dim ? ydir : xdir);
				float const dx(xdir ? -xsz : xsz), dy(ydir ? -ysz : ysz), dw(car_dim ? dx : dy), dr(car_dim ? dy : dx); // delta-wdith and delta-row
				point const corner_pos(plot.d[0][xdir], plot.d[1][ydir], (plot.z1() + 0.1*ROAD_HEIGHT)); // shift up slightly to avoid z-fighting
				assert(dw != 0.0 && dr != 0.0);
				parking_lot_t cand(cube_t(corner_pos, corner_pos), car_dim, car_dir, city_params.min_park_spaces, city_params.min_park_rows); // start as min size at the corner
				cand.d[!car_dim][!wdir] += cand.row_sz*dw;
				cand.d[ car_dim][!rdir] += cand.num_rows*dr;
				if (!plot.contains_cube_xy(cand)) {continue;} // can't fit a min size parking lot in this plot, so skip it (shouldn't happen)
				if (has_bcube_int_xy(cand, bcubes, pad_dist)) continue; // intersects a building - skip (can't fit min size parking lot)
				cand.z2() += plot.dz(); // probably unnecessary
				parking_lot_t park(cand);

				// try to add more parking spaces in a row
				for (; plot.contains_cube_xy(cand); ++cand.row_sz, cand.d[!car_dim][!wdir] += dw) {
					if (has_bcube_int_xy(cand, bcubes, pad_dist)) break; // intersects a building - done
					park = cand; // success: increase parking lot to this size
				}
				cand = park;
				// try to add more rows of parking spaces
				for (; plot.contains_cube_xy(cand); ++cand.num_rows, cand.d[car_dim][!rdir] += dr) {
					if (has_bcube_int_xy(cand, bcubes, pad_dist)) break; // intersects a building - done
					park = cand; // success: increase parking lot to this size
				}
				assert(park.row_sz >= city_params.min_park_spaces && park.num_rows >= city_params.min_park_rows);
				assert(park.dx() > 0.0 && park.dy() > 0.0);
				car.cur_seg = (unsigned short)parking_lots.size(); // store parking lot index in cur_seg
				add_obj_to_group(park, park, parking_lots, parking_lot_groups, is_new_tile);
				bcubes.push_back(park); // add to list of blocker bcubes so that no later parking lots overlap this one
				//colliders.push_back(park); // added per-filled space below
				//parking_lots.back().expand_by_xy(0.5*pad_dist); // re-add half the padding for drawing (breaks texture coord alignment)
				unsigned const nspaces(park.row_sz*park.num_rows);
				num_spaces += nspaces;

				// fill the parking lot with cars
				vector<unsigned char> &used_spaces(parking_lots.back().used_spaces);
				used_spaces.resize(nspaces, 0); // start empty
				vector3d car_sz(nom_car_size);
				car.dim    = car_dim;
				car.dir    = car_dir;
				car.height = car_sz.z;
				if (car.dim) {swap(car_sz.x, car_sz.y);}
				point pos(corner_pos.x, corner_pos.y, (plot.z2() + 0.5*car_sz.z));
				pos[ car_dim] += 0.5*dr + (car_dim ? 0.15 : -0.15)*fabs(dr); // offset for centerline, biased toward the front of the parking space
				float const car_density(rgen.rand_uniform(city_params.min_park_density, city_params.max_park_density));

				for (unsigned row = 0; row < park.num_rows; ++row) {
					pos[!car_dim] = corner_pos[!car_dim] + 0.5*dw; // half offset for centerline
					bool prev_was_bad(0);

					for (unsigned col = 0; col < park.row_sz; ++col) { // iterate one past the end
						if (prev_was_bad) {prev_was_bad = 0;} // previous car did a bad parking job, leave this space empty
						else if (rgen.rand_float() < car_density) { // only half the spaces are filled on average
							point cpos(pos);
							cpos[ car_dim] += 0.05*dr*rgen.rand_uniform(-1.0, 1.0); // randomness of front amount
							cpos[!car_dim] += 0.12*dw*rgen.rand_uniform(-1.0, 1.0); // randomness of side  amount

							if (col+1 != park.row_sz && (rgen.rand()&15) == 0) {// occasional bad parking job
								cpos[!car_dim] += dw*rgen.rand_uniform(0.3, 0.35);
								prev_was_bad = 1;
							}
							car.bcube.set_from_point(cpos);
							car.bcube.expand_by(0.5*car_sz);
							cars.push_back(car);
							if ((rgen.rand()&7) == 0) {cars.back().dir ^= 1;} // pack backwards 1/8 of the time
							used_spaces[row*park.num_rows + col] = 1;
							++filled_spaces;
							has_parking = 1;
						}
						pos[!car_dim] += dw;
					} // for col
					pos[car_dim] += dr;
				} // for row
				// generate colliders for each group of used parking space columns
				cube_t cur_cube(park); // set zvals, etc.
				bool inside(0);

				for (unsigned col = 0; col <= park.row_sz; ++col) {
					// mark this space as blocked if any spaces in the row are blocked; this avoids creating diagonally adjacent colliders that cause dead ends and confuse path finding
					bool blocked(0);
					for (unsigned row = 0; col < park.row_sz && row < park.num_rows; ++row) {blocked |= (used_spaces[row*park.num_rows + col] != 0);}

					if (!inside && blocked) { // start a new segment
						cur_cube.d[!car_dim][0] = corner_pos[!car_dim] + col*dw;
						inside = 1;
					}
					else if (inside && !blocked) { // end the current segment
						cur_cube.d[!car_dim][1] = corner_pos[!car_dim] + col*dw;
						cur_cube.normalize();
						//assert(park.contains_cube(cur_cube)); // can fail due to floating-point precision
						colliders.push_back(cur_cube);
						inside = 0;
					}
				} // for col
			} // for c
			return has_parking;
		}
		static bool check_pt_and_place_blocker(point const &pos, vect_cube_t &blockers, float radius, float blocker_spacing) {
			cube_t bc(pos);
			if (has_bcube_int_xy(bc, blockers, radius)) return 0; // intersects a building or parking lot - skip
			bc.expand_by_xy(blocker_spacing);
			blockers.push_back(bc); // prevent trees and benches from being too close to each other
			return 1;
		}
		static bool try_place_obj(cube_t const &plot, vect_cube_t &blockers, rand_gen_t &rgen, float radius, float extra_spacing, float num_tries, point &pos) {
			for (unsigned t = 0; t < num_tries; ++t) {
				pos = rand_xy_pt_in_cube(plot, radius, rgen);
				if (check_pt_and_place_blocker(pos, blockers, radius, extra_spacing)) return 1; // success
			} // for t
			return 0;
		}
		static void place_tree(point const &pos, float radius, int ttype, vect_cube_t &colliders, vector<point> &tree_pos, bool allow_bush, bool is_sm_tree) {
			tree_placer.add(pos, 0, ttype, allow_bush, is_sm_tree); // use same tree type
			cube_t bcube; bcube.set_from_sphere(pos, 0.1*radius); // use 10% of the placement radius for collision
			bcube.z2() += radius; // increase cube height
			colliders.push_back(bcube);
			tree_pos.push_back(pos);
		}
		static void place_trees_in_plot(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> &tree_pos, rand_gen_t &rgen) {
			if (city_params.max_trees_per_plot == 0) return;
			float const radius(city_params.tree_spacing*city_params.get_nom_car_size().x); // in multiples of car length
			float const spacing(max(radius, get_min_obj_spacing())), radius_exp(2.0*spacing);
			vector3d const plot_sz(plot.get_size());
			if (min(plot_sz.x, plot_sz.y) < 2.0*radius_exp) return; // plot is too small for trees of this size
			unsigned num_trees(city_params.max_trees_per_plot);
			if (plot.is_park) {num_trees += (rgen.rand() % city_params.max_trees_per_plot);} // allow up to twice as many trees in parks

			for (unsigned n = 0; n < num_trees; ++n) {
				bool const is_sm_tree((rgen.rand()%3) == 0); // 33% of the time is a pine/palm tree
				int ttype(-1); // Note: okay to leave at -1; also, don't have to set to a valid tree type
				if (is_sm_tree) {ttype = (plot.is_park ? (rgen.rand()&1) : 2);} // pine/short pine in parks, palm in city blocks
				else {ttype = rgen.rand()%100;} // random type
				bool const is_palm(is_sm_tree && ttype == 2);
				bool const allow_bush(plot.is_park && max_unique_trees == 0); // can't place bushes if tree instances are enabled (generally true) because bushes may be instanced in non-parks
				float const bldg_extra_radius(is_palm ? 0.5f*radius : 0.0f); // palm trees are larger and must be kept away from buildings, but can overlap with other trees
				point pos;
				if (!try_place_obj(plot, blockers, rgen, (spacing + bldg_extra_radius), (radius - bldg_extra_radius), 10, pos)) continue; // 10 tries per tree, extra spacing for palm trees
				place_tree(pos, radius, ttype, colliders, tree_pos, allow_bush, is_sm_tree); // size is randomly selected by the tree generator using default values; allow bushes in parks
				if (plot.is_park) continue; // skip row logic and just place trees randomly throughout the park
				// now that we're here, try to place more trees at this same distance from the road in a row
				bool const dim(min((pos.x - plot.x1()), (plot.x2() - pos.x)) < min((pos.y - plot.y1()), (plot.y2() - pos.y)));
				bool const dir((pos[dim] - plot.d[dim][0]) < (plot.d[dim][1] - pos[dim]));
				float const step(1.25*radius_exp*(dir ? 1.0 : -1.0)); // positive or negative (must be > 2x radius spacing)
					
				for (; n < city_params.max_trees_per_plot; ++n) {
					pos[dim] += step;
					if (pos[dim] < plot.d[dim][0]+radius || pos[dim] > plot.d[dim][1]-radius) break; // outside place area
					if (!check_pt_and_place_blocker(pos, blockers, (spacing + bldg_extra_radius), (spacing - bldg_extra_radius))) break; // placement failed
					place_tree(pos, radius, ttype, colliders, tree_pos, plot.is_park, is_sm_tree); // use same tree type
				} // for n
			} // for n
		}
		template<typename T> void add_obj_to_group(T const &obj, cube_t const &bcube, vector<T> &objs, vector<cube_with_ix_t> &groups, bool &is_new_tile) {
			objs.push_back(obj);
			if (groups.empty() || is_new_tile) {groups.push_back(cube_with_ix_t(bcube));}
			else {groups.back().union_with_cube(bcube);}
			groups.back().ix = objs.size();
			is_new_tile = 0;
		}
		template<typename T> void sort_grouped_objects(vector<T> &objs, vector<cube_with_ix_t> const &groups) {
			unsigned start_ix(0);
			for (auto i = groups.begin(); i != groups.end(); start_ix = i->ix, ++i) {sort(objs.begin()+start_ix, objs.begin()+i->ix);}
		}
		// Note: blockers are used for placement of objects within this plot; colliders are used for pedestrian AI
		void place_detail_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> const &tree_pos, rand_gen_t &rgen, bool is_new_tile) {
			float const car_length(city_params.get_nom_car_size().x); // used as a size reference for other objects
			// FIXME: due to some problem I can't figure out, bench shadow maps don't work right unless is_new_bench_tile is reset for each plot;
			//        however, tree planters seem to work just fine without this, even though they use nearly the same code
			bool is_new_fh_tile(is_new_tile), is_new_bench_tile(1 || is_new_tile), is_new_planter_tile(is_new_tile);

			// place fire_hydrants; don't add fire hydrants in parks
			if (!plot.is_park) {
				float const radius(0.04*car_length), height(0.18*car_length), dist_from_road(radius);
				point pos(0.0, 0.0, plot.z2()); // XY will be assigned below

				for (unsigned dim = 0; dim < 2; ++dim) {
					pos[!dim] = plot.get_center_dim(!dim);

					for (unsigned dir = 0; dir < 2; ++dir) {
						pos[dim] = plot.d[dim][dir] - (dir ? 1.0 : -1.0)*dist_from_road; // move into the sidewalk along the road
						// Note: will skip placement if too close to a previously placed tree, but that should be okay as it is relatively uncommon
						if (!check_pt_and_place_blocker(pos, blockers, radius, 2.0*radius)) continue; // bad placement, skip
						vector3d orient(zero_vector);
						orient[!dim] = (dir ? 1.0 : -1.0); // oriented perpendicular to the road
						fire_hydrant_t const fire_hydrant(pos, radius, height, orient);
						add_obj_to_group(fire_hydrant, fire_hydrant.bcube, fire_hydrants, fire_hydrant_groups, is_new_fh_tile);
						colliders.push_back(fire_hydrant.bcube);
					} // for dir
				} // for dim
			}
			// place benches
			bench_t bench;
			bench.radius = 0.3*car_length;
			float const bench_spacing(max(bench.radius, get_min_obj_spacing()));

			for (unsigned n = 0; n < city_params.max_benches_per_plot; ++n) {
				if (!try_place_obj(plot, blockers, rgen, bench_spacing, 0.0, 1, bench.pos)) continue; // 1 try
				float dmin(0.0);

				for (unsigned dim = 0; dim < 2; ++dim) {
					for (unsigned dir = 0; dir < 2; ++dir) {
						float const dist(fabs(bench.pos[dim] - plot.d[dim][dir])); // find closest distance to road (plot edge) and orient bench that way
						if (dmin == 0.0 || dist < dmin) {bench.dim = !dim; bench.dir = !dir; dmin = dist;}
					}
				}
				bench.calc_bcube();
				add_obj_to_group(bench, bench.bcube, benches, bench_groups, is_new_bench_tile);
				colliders.push_back(bench.bcube);
			} // for n

			// place planters; don't add planters in parks
			if (!plot.is_park) {
				float const planter_height(0.05*car_length), planter_radius(0.25*car_length);

				for (auto i = tree_pos.begin(); i != tree_pos.end(); ++i) {
					tree_planter_t const planter(*i, planter_radius, planter_height);
					add_obj_to_group(planter, planter.bcube, planters, planter_groups, is_new_planter_tile); // no colliders for planters; pedestrians avoid the trees instead
				}
			}
		}
		template<typename T> void draw_objects(vector<T> const &objs, vector<cube_with_ix_t> const &groups,
			draw_state_t &dstate, float dist_scale, bool shadow_only, bool not_using_qbd=0)
		{
			if (objs.empty()) return;
			T::pre_draw(dstate, shadow_only);
			unsigned start_ix(0);
			assert(qbd.empty());

			for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
				if (!dstate.check_cube_visible(*g, dist_scale, shadow_only)) continue; // VFC/distance culling for group
				if (not_using_qbd) {dstate.begin_tile(g->get_cube_center(), 1, 1);} // must setup shader and tile shadow map before drawing

				for (unsigned i = start_ix; i < g->ix; ++i) {
					assert(i < objs.size());
					T const &obj(objs[i]);
					if (dstate.check_sphere_visible(obj.pos, obj.radius)) {obj.draw(dstate, qbd, dist_scale, shadow_only);}
				}
				if (!qbd.empty()) { // we have something to draw
					dstate.begin_tile(g->get_cube_center(), 1, 1); // will_emit_now=1, ensure_active=1
					qbd.draw_and_clear(); // draw this group with current smap
				}
			} // for g
			T::post_draw(dstate, shadow_only);
		}
	public:
		city_obj_placer_t() : num_spaces(0), filled_spaces(0) {}
		
		void clear() {
			parking_lots.clear(); parking_lot_groups.clear(); benches.clear(); planters.clear(); fire_hydrants.clear();
			bench_groups.clear(); planter_groups.clear(); fire_hydrant_groups.clear();
			num_spaces = filled_spaces = 0;
		}
		struct cube_by_x1 {
			bool operator()(cube_t const &a, cube_t const &b) const {return (a.x1() < b.x1());}
		};
		void gen_parking_and_place_objects(vector<road_plot_t> &plots, vector<vect_cube_t> &plot_colliders, vector<car_t> &cars, unsigned city_id, bool have_cars) {
			// Note: fills in plots.has_parking
			//timer_t timer("Gen Parking Lots and Place Objects");
			vect_cube_t bcubes; // local blockers for this plot; reused across calls
			vector<point> tree_pos;
			rand_gen_t rgen, detail_rgen;
			parking_lots.clear();
			rgen.set_state(city_id, 123);
			detail_rgen.set_state(3145739*(city_id+1), 1572869*(city_id+1));
			clear();
			if (city_params.max_trees_per_plot > 0) {tree_placer.begin_block(0); tree_placer.begin_block(1);} // both small and large trees
			bool const add_parking_lots(have_cars && city_params.min_park_spaces > 0 && city_params.min_park_rows > 0);
			uint64_t prev_tile_id(0);

			for (auto i = plots.begin(); i != plots.end(); ++i) {
				uint64_t const tile_id(road_network_t::get_tile_id_for_cube(*i));
				bool const is_new_tile(tile_id != prev_tile_id);
				bcubes.clear();
				tree_pos.clear();
				get_building_bcubes(*i, bcubes);
				size_t const plot_id(i - plots.begin());
				assert(plot_id < plot_colliders.size());
				vect_cube_t &colliders(plot_colliders[plot_id]);
				if (add_parking_lots && !i->is_park) {i->has_parking = gen_parking_lots_for_plot(*i, cars, city_id, plot_id, bcubes, colliders, rgen, is_new_tile);}
				place_trees_in_plot (*i, bcubes, colliders, tree_pos, detail_rgen);
				place_detail_objects(*i, bcubes, colliders, tree_pos, detail_rgen, is_new_tile);
				sort(colliders.begin(), colliders.end(), cube_by_x1());
				prev_tile_id = tile_id;
			} // for i
			sort_grouped_objects(benches,       bench_groups  );
			sort_grouped_objects(planters,      planter_groups);
			sort_grouped_objects(fire_hydrants, fire_hydrant_groups);

			if (add_parking_lots) {
				cout << "parking lots: " << parking_lots.size() << ", spaces: " << num_spaces << ", filled: " << filled_spaces << ", benches: " << benches.size() << endl;
			}
		}
		void draw_detail_objects(draw_state_t &dstate, bool shadow_only) {
			draw_objects(benches,       bench_groups,        dstate, 0.16, shadow_only, 0); // dist_scale=0.16
			draw_objects(fire_hydrants, fire_hydrant_groups, dstate, 0.07, shadow_only, 1); // dist_scale=0.12, not_using_qbd=1
			
			if (!shadow_only) { // low profile, not drawn in shadow pass
				draw_objects(planters, planter_groups, dstate, 0.1, shadow_only, 0); // dirt pass, dist_scale=0.1
				dstate.pass_ix = 1;
				draw_objects(planters, planter_groups, dstate, 0.1, shadow_only, 0); // stone pass, dist_scale=0.1
				dstate.pass_ix = 0;
			}
		}
		bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm) const {
			vector3d const xlate(get_camera_coord_space_xlate());

			for (auto i = benches.begin(); i != benches.end(); ++i) { // Note: could use bench_groups
				if (i->proc_sphere_coll(pos, p_last, radius, xlate, cnorm)) return 1;
			}
			for (auto i = fire_hydrants.begin(); i != fire_hydrants.end(); ++i) { // Note: could use fire_hydrant_groups
				if (i->proc_sphere_coll(pos, p_last, radius, xlate, cnorm)) return 1;
			}
			// Note: no coll with tree_planters because the tree coll should take care of it
			return 0;
		}
		bool line_intersect(point const &p1, point const &p2, float &t) const { // Note: nothing to do for parking lots or tree_planters
			bool ret(0);

			for (auto i = benches.begin(); i != benches.end(); ++i) {
				ret |= check_line_clip_update_t(p1, p2, t, i->bcube); // check bounding cube
			}
			for (auto i = fire_hydrants.begin(); i != fire_hydrants.end(); ++i) {
				ret |= check_line_clip_update_t(p1, p2, t, i->bcube); // check bounding cube; cylinder intersection may be more accurate, but likely doesn't matter much
			}
			return ret;
		}
		bool pt_in_parking_lot_xy(point const &pos) const {
			unsigned start_ix(0);

			for (auto i = parking_lot_groups.begin(); i != parking_lot_groups.end(); start_ix = i->ix, ++i) {
				if (!i->contains_pt_xy(pos)) continue;
				for (unsigned b = start_ix; b < i->ix; ++b) {if (parking_lots[b].contains_pt_xy(pos)) return 1;}
			}
			return 0;
		}
		bool cube_overlaps_parking_lot_xy(cube_t const &c) const {
			unsigned start_ix(0);

			for (auto i = parking_lot_groups.begin(); i != parking_lot_groups.end(); start_ix = i->ix, ++i) {
				if (!i->intersects_xy(c)) continue;
				for (unsigned b = start_ix; b < i->ix; ++b) {if (parking_lots[b].intersects_xy(c)) return 1;}
			}
			return 0;
		}
		bool get_color_at_xy(point const &pos, colorRGBA &color) const {
			unsigned start_ix(0);

			for (auto i = bench_groups.begin(); i != bench_groups.end(); start_ix = i->ix, ++i) {
				if (!i->contains_pt_xy(pos)) continue;
					
				for (auto b = benches.begin()+start_ix; b != benches.begin()+i->ix; ++b) {
					if (pos.x < b->bcube.x1()) break; // benches are sorted by x1, no bench after this can match
					if (b->bcube.contains_pt_xy(pos)) {color = texture_color(FENCE_TEX); return 1;}
				}
			} // for i
			float const expand(0.15*city_params.road_width), x_test(pos.x + expand); // expand to approx tree diameter
			start_ix = 0;

			for (auto i = planter_groups.begin(); i != planter_groups.end(); start_ix = i->ix, ++i) {
				if (!i->contains_pt_xy_exp(pos, expand)) continue;

				for (auto p = planters.begin()+start_ix; p != planters.begin()+i->ix; ++p) {
					if (x_test < p->bcube.x1()) break; // planters are sorted by x1, no planter after this can match
					if (!p->bcube.contains_pt_xy_exp(pos, expand)) continue;
					// treat this as a tree rather than a planter by testing against a circle, since trees aren't otherwise included
					if (dist_xy_less_than(pos, p->pos, (p->radius + expand))) {color = DK_GREEN; return 1;}
				}
			} // for i
			start_ix = 0;

			for (auto i = fire_hydrant_groups.begin(); i != fire_hydrant_groups.end(); start_ix = i->ix, ++i) {
				if (!i->contains_pt_xy(pos)) continue;

				for (auto b = fire_hydrants.begin()+start_ix; b != fire_hydrants.begin()+i->ix; ++b) {
					if (pos.x < b->bcube.x1()) break; // fire_hydrants are sorted by x1, no fire_hydrant after this can match
					if (dist_xy_less_than(pos, b->pos, b->radius)) {color = colorRGBA(1.0, 0.75, 0.0); return 1;} // orange/yellow color
				}
			} // for i
			return 0;
		}
	}; // city_obj_placer_t


	class road_network_t : public streetlights_t {

		vector<road_t> roads; // full overlapping roads with constant slope, for collisions, etc.
		vector<road_seg_t> segs; // non-overlapping road segments, for drawing with textures
		vector<cube_t> conn_roads; // connector road bounding cubes (contain multiple adjacent connected roads with different slopes/zvals)
		vector<road_isec_t> isecs[3]; // for drawing with textures: {2-way, 3-way, 4-way}
		vector<road_plot_t> plots; // plots of land that can hold buildings (city blocks)
		vector<bridge_t> bridges; // bridges, part of global road network
		vector<tunnel_t> tunnels; // tunnels, part of global road network
		vector<road_t> tracks, track_segs; // railroad tracks (for global road network)
		vect_cube_t parks;
		//vector<road_isec_t> track_turns; // for railroad tracks
		city_obj_placer_t city_obj_placer;
		cube_t bcube;
		vector<road_t> segments; // reused temporary
		set<unsigned> connected_to; // vector?
		map<uint64_t, unsigned> tile_to_block_map;
		map<unsigned, road_isec_t const *> cix_to_isec; // maps city_ix to intersection
		vector<vect_cube_t> plot_colliders;
		plot_xy_t plot_xy;
		unsigned city_id, cluster_id, plot_id_offset;
		//string city_name; // future work
		float tot_road_len;
		mutable unsigned num_cars; // Note: not counting parked cars; mutable so that car_manager can update this

		// use only for the global road network
		struct city_id_pair_t {
			unsigned id[2]; // lo, hi
			city_id_pair_t(unsigned c1, unsigned c2) {id[0] = c1; id[1] = c2;}
		};
		vector<city_id_pair_t> road_to_city; // indexed by road ID
		vector<vector<unsigned>> city_to_seg; // maps city_id to set of road segments connecting to that city

		struct cmp_by_tile { // not the most efficient solution, but no memory overhead
			bool operator()(cube_t const &a, cube_t const &b) const {return (get_tile_id_for_cube(a) < get_tile_id_for_cube(b));}
		};
		struct tile_block_t { // collection of road parts for a given tile
			range_pair_t ranges[NUM_RD_TYPES]; // {plot, seg, isec2, isec3, isec4, park_lot, tracks}
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
	public:
		road_network_t() : bcube(all_zeros), city_id(CONN_CITY_IX), cluster_id(0), plot_id_offset(0), tot_road_len(0.0), num_cars(0) {} // global road network ctor
		road_network_t(cube_t const &bcube_, unsigned city_id_) : bcube(bcube_), city_id(city_id_), cluster_id(0), plot_id_offset(0), tot_road_len(0.0), num_cars(0) {
			bcube.d[2][1] += ROAD_HEIGHT; // make it nonzero size
		}
		static uint64_t get_tile_id_for_cube(cube_t const &c) {return get_tile_id_containing_point_no_xyoff(c.get_cube_center());}
		cube_t const &get_bcube() const {return bcube;}
		cube_t const &get_plot_bcube(unsigned plot_ix) const {assert(plot_ix < plots.size()); return plots[plot_ix];}
		void set_bcube(cube_t const &bcube_) {bcube = bcube_;}
		unsigned num_roads() const {return roads.size();}
		vector<road_t> const &get_roads() const {return roads;} // used for connecting roads between cities with 4-way intersections
		bool empty() const {return roads.empty();}
		plot_xy_t const &get_plot_xy() const {return plot_xy;}
		bool has_tunnels() const {return !tunnels.empty();}
		void set_cluster(unsigned id) {cluster_id = id;}
		void register_connected_city(unsigned id) {connected_to.insert(id);}
		set<unsigned> const &get_connected() const {return connected_to;}
		bool is_connected_to(unsigned id) const {return (connected_to.find(id) != connected_to.end());}
		float get_traffic_density() const {return ((tot_road_len == 0.0) ? 0.0 : num_cars/tot_road_len);} // cars per unit road
		void register_car() const {++num_cars;} // Note: must be const; num_cars is mutable

		void clear() {
			roads.clear();
			segs.clear();
			conn_roads.clear();
			plots.clear();
			bridges.clear();
			tunnels.clear();
			tracks.clear();
			track_segs.clear();
			for (unsigned i = 0; i < 3; ++i) {isecs[i].clear();}
			streetlights.clear();
			city_obj_placer.clear();
			tile_blocks.clear();
			plot_colliders.clear();
		}
		bool gen_road_grid(float road_width, float road_spacing) {
			if (city_params.road_width > 0.5*city_params.road_spacing) {
				cerr << "Error: City road_width should not be set larger than half the road spacing" << endl;
				exit(1);
			}
			cube_t const &region(bcube); // use our bcube as the region to process
			vector3d const size(region.get_size());
			assert(size.x > 0.0 && size.y > 0.0);
			//rand_gen_t rgen; rgen.set_state(int(123.0*region.x1()), int(456.0*region.y1())); road_spacing *= rgen.rand_uniform(0.8, 1.2); // add some random variation?
			float const half_width(0.5*road_width), zval(region.z1() + ROAD_HEIGHT);
			float const rx1(region.x1() + half_width), rx2(region.x2() - half_width), ry1(region.y1() + half_width), ry2(region.y2() - half_width); // shrink to include centerlines
			float road_pitch_x(road_width + road_spacing), road_pitch_y(road_pitch_x);
			int const num_x_roads((rx2 - rx1)/road_pitch_x), num_y_roads((ry2 - ry1)/road_pitch_y);
			road_pitch_x = 0.9999f*(rx2 - rx1)/num_x_roads; // auto-calculate, round down slightly to avoid FP error
			road_pitch_y = 0.9999f*(ry2 - ry1)/num_y_roads;
			//cout << "road pitch: " << road_pitch_x/DX_VAL << " " << road_pitch_y/DY_VAL << " road width: " << road_width/DX_VAL << " " << road_width/DY_VAL << endl;
			//spacing = int((road_width + road_spacing)/DX_VAL)*DX_VAL - road_width;

			// create a grid, for now; crossing roads will overlap
			for (float x = rx1; x < rx2; x += road_pitch_x) {
				roads.emplace_back(point(x, region.y1(), zval), point(x, region.y2(), zval), road_width, true);
			}
			unsigned const num_x(roads.size());

			for (float y = ry1; y < ry2; y += road_pitch_y) {
				roads.emplace_back(point(region.x1(), y, zval), point(region.x2(), y, zval), road_width, false);
			}
			unsigned const num_r(roads.size()), num_y(num_r - num_x);
			if (num_x <= 1 || num_y <= 1) {clear(); return 0;} // not enough space for roads
			bcube.x1() = roads[0      ].x1(); // actual bcube x1 from first x road
			bcube.x2() = roads[num_x-1].x2(); // actual bcube x2 from last  x road
			bcube.y1() = roads[num_x  ].y1(); // actual bcube y1 from first y road
			bcube.y2() = roads[num_r-1].y2(); // actual bcube y2 from last  y road

			// create road segments and intersections
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
					isecs[num_conn - 2].emplace_back(cube_t(rx.x1(), rx.x2(), ry.y1(), ry.y2(), zval, zval), y, x, conn, false); // intersections
					
					if (!LX) { // skip last y segment
						cube_t const &rxn(roads[x+1]);
						segs.emplace_back(cube_t(rx.x2(), rxn.x1(), ry.y1(), ry.y2(), zval, zval), y, false); // y-segments
					}
					if (!LY) { // skip last x segment
						cube_t const &ryn(roads[y+1]);
						segs.emplace_back(cube_t(rx.x1(), rx.x2(), ry.y2(), ryn.y1(), zval, zval), x, true); // x-segments

						if (!LX) { // skip last y segment
							cube_t const &rxn(roads[x+1]);
							plots.emplace_back(cube_t(rx.x2(), rxn.x1(), ry.y2(), ryn.y1(), zval, zval), x, (y - num_x)); // plots between roads
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
		void gen_railroad_tracks(float width, unsigned num, cube_t const &region, vect_cube_t const &blockers, heightmap_query_t &hq) {
			rand_gen_t rgen;
			assert(region.dx() > 0.0 && region.dy() > 0.0);
			if (region.dx() <= 2.0*width || region.dy() <= 2.0*width) return; // region too small (shouldn't happen)
			vect_cube_t dim_tracks[2]; // one in each dim, for collision detection with tracks going in the other dim

			for (unsigned n = 0; n < num; ++n) {
				for (unsigned tries = 0; tries < city_params.num_conn_tries; ++tries) {
					bool const dim(rgen.rand_bool());
					float const rv(rgen.rand_uniform(0.2, 0.8)); // use center area
					float const pos(region.d[!dim][0]*(1.0f - rv) + region.d[!dim][1]*rv);
					float const step_sz(city_params.conn_road_seg_len);
					float const seg_end(region.d[dim][1]);
					point p1, p2;
					p1[!dim] = p2[!dim] = pos;
					p1[dim]  = region.d[dim][0]; p2[dim] = seg_end; // full segment for blockers check
					p1.z = p2.z = region.z1();
					cube_t tracks_bcube(p1, p2);
					if (has_bcube_int_xy(tracks_bcube, blockers,            width)) continue; // check cities
					if (has_bcube_int_xy(tracks_bcube, dim_tracks[dim], 8.0*width)) continue; // check prev placed tracks in same dim
					dim_tracks[dim].push_back(tracks_bcube); // add to dim_tracks, but not to blockers, since we want roads to cross tracks
					tracks.emplace_back(p1, p2, width, dim, (p2.z < p1.z), n); // Note: zvals are at 0, but should be unused
					p2[dim] = (p1[dim] + step_sz); // back to starting segment

					while (p1[dim] < seg_end) { // split into per-tile segments
						p1.z = hq.get_height_at(p1.x, p1.y) + ROAD_HEIGHT;
						p2.z = hq.get_height_at(p2.x, p2.y) + ROAD_HEIGHT;
						track_segs.emplace_back(p1, p2, width, dim, (p2.z < p1.z), n);
						p1[dim] += step_sz;
						p2[dim]  = min((p1[dim] + step_sz), seg_end);
					} // end while
					// TODO: check for collisions with roads and handle them with intersections, bridges, or tunnels
					// TODO: handle slopes that are too steep and have shadow artifacts
					break; // success
				} // for tries
			} // for n
			for (unsigned pass = 0; pass < 2; ++pass) { // flatten mesh after placing all tracks: regular, decrease_only
				for (auto i = track_segs.begin(); i != track_segs.end(); ++i) {hq.flatten_for_road(*i, city_params.road_border, 0, (pass == 1));}
			}
			cout << "track segments: " << track_segs.size() << endl;
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
						bool found2(0);

						for (unsigned n = 1; n < 3; ++n) { // search 3-way and 4-way intersections
							all_ixs.resize(rn.isecs[n].size());
							for (unsigned m = 0; m < all_ixs.size(); ++m) {all_ixs[m] = m;} // all sequential index values
							int const isec_ix(rn.search_for_adj(rn.isecs[n], all_ixs, seg, seg.dim, (dir != 0)));
							if (isec_ix < 0) continue; // not be found
							seg.conn_ix  [dir] = isec_ix;
							seg.conn_type[dir] = TYPE_ISEC2 + n; // always connects to a 3-way or 4-way intersection within the city
							found2 = 1;
							break;
						} // for n
						assert(found2); // must be found
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
				isecs[1].emplace_back(ibc, (dim ? seg.road_ix : (int)other_rix), (dim ? other_rix : (int)seg.road_ix), conns[2*(!dim) + dir], true, dest_city_id); // 3-way
			}
		}
		float create_connector_road(cube_t const &bcube1, cube_t const &bcube2, vect_cube_t &blockers, road_network_t *rn1, road_network_t *rn2, unsigned city1, unsigned city2,
			unsigned dest_city_id1, unsigned dest_city_id2, heightmap_query_t &hq, float road_width, float conn_pos, bool dim, bool check_only, bool is_4_way1, bool is_4_way2)
		{
			bool const dir(bcube1.d[dim][0] < bcube2.d[dim][0]);
			if (dir == 0) {swap(city1, city2);} // make {lo, hi}
			point p1, p2;
			p1.z = bcube1.d[2][1];
			p2.z = bcube2.d[2][1];
			p1[!dim] = p2[!dim] = conn_pos;
			p1[ dim] = bcube1.d[dim][ dir];
			p2[ dim] = bcube2.d[dim][!dir];
			bool const slope((p1.z < p2.z) ^ dir);
			road_t const road(p1, p2, road_width, dim, slope, roads.size());
			float const road_len(road.get_length()), delta_z(road.dz()), max_slope(city_params.max_road_slope);
			assert(road_len > 0.0 && delta_z >= 0.0);
			if (delta_z/road_len > max_slope) {assert(check_only); return -1.0;} // slope is too high (split segments will have even higher slopes)
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
				} // Note: no bridges here, but could add them
				return hq.flatten_sloped_region(x1, y1, x2, y2, road.d[2][slope]-ROAD_HEIGHT, road.d[2][!slope]-ROAD_HEIGHT, dim, city_params.road_border, 0, 0, check_only);
			}
			unsigned const num_segs(ceil(road_len/city_params.conn_road_seg_len));
			assert(num_segs > 0 && num_segs < 1000); // sanity check
			float const seg_len(road_len/num_segs);
			assert(seg_len <= city_params.conn_road_seg_len);
			road_t rs(road); // keep d[!dim][0], d[!dim][1], dim, and road_ix
			rs.z1() = road.d[2][slope];
			segments.clear();
			float tot_dz(0.0);
			bool last_was_bridge(0), last_was_tunnel(0);
			vector<flatten_op_t> replay_fops;

			for (unsigned n = 0; n < num_segs; ++n) {
				rs.d[dim][1] = ((n+1 == num_segs) ? road.d[dim][1] : (rs.d[dim][0] + seg_len)); // make sure it ends exactly at the correct location
				point pos;
				pos[ dim] = rs.d[dim][1];
				pos[!dim] = conn_pos;
				rs.z2()   = hq.get_height_at(pos.x, pos.y) + ROAD_HEIGHT; // terrain height at end of segment
				rs.slope  = (rs.z2() < rs.z1());
				
				if (fabs(rs.get_slope_val()) > max_slope) { // slope is too high, clamp z2 to max allowed value
					if (n+1 == num_segs) {
						// Note: the height of the first/last segment may change slightly after placing the bend,
						// which can make this slope check fail when check_only=0 while it passed when check_only=1;
						// returning here will create a disconnected road segment and fail an assert later, so instead we allow the high slope
						if (check_only) return -1.0;
					}
					else {rs.z2() = rs.z1() + seg_len*max_slope*SIGN(rs.dz());}
				}
				segments.push_back(rs);
				rs.d[dim][0] = rs.d[dim][1]; rs.z1() = rs.z2(); // shift segment end point
			} // for n
			for (auto s = segments.begin(); s != segments.end(); ++s) {
				if (s->z2() < s->z1()) {swap(s->z2(), s->z1());} // swap zvals if needed
				assert(s->is_normalized());
				bridge_t bridge(*s);
				tunnel_t tunnel(*s);
				tot_dz += hq.flatten_for_road(*s, city_params.road_border, check_only, 0, (last_was_bridge ? nullptr : &bridge), (last_was_tunnel ? nullptr : &tunnel));
				replay_fops.push_back(hq.last_flatten_op);
				
				if (!check_only) {
					roads.push_back(*s);
					road_to_city.emplace_back(city1, city2); // Note: city index is specified even for internal (non-terminal) roads
					if (bridge.make_bridge) {bridges.push_back(bridge);}
					if (tunnel.enabled()) {tunnels.push_back(tunnel);}
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
			isecs[0].emplace_back(int_bcube, road_ix_x, road_ix_y, conns[2*dy + dx], true);
			//blockers.push_back(int_bcube); // ???
		}
		void split_connector_roads(float road_spacing) {
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
				cube_t c(*r); // start by copying the road's bcube
				
				for (unsigned n = 0; n < num_segs; ++n) {
					c.d[d][1] = ((n+1 == num_segs) ? r->d[d][1] : (c.d[d][0] + seg_len)); // make sure it ends exactly at the correct location
					for (unsigned e = 0; e < 2; ++e) {c.d[2][e] = z1 + (z2 - z1)*((c.d[d][e] - r->d[d][0])/len);} // interpolate road height across segments
					if (c.z2() < c.z1()) {swap(c.z2(), c.z1());} // swap zvals if needed
					assert(c.is_normalized());
					segs.emplace_back(c, rix, d, r->slope);
					c.d[d][0] = c.d[d][1]; // shift segment end point
				} // for n
			} // for r
		}
		void finalize_bridges_and_tunnels() {
			for (auto b = bridges.begin(); b != bridges.end(); ++b) {b->add_streetlights();}
			for (auto b = tunnels.begin(); b != tunnels.end(); ++b) {b->add_streetlights();}
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
		void gen_parking_lots_and_place_objects(vector<car_t> &cars, bool have_cars) {
			city_obj_placer.gen_parking_and_place_objects(plots, plot_colliders, cars, city_id, have_cars);
			add_tile_blocks(city_obj_placer.parking_lots, tile_to_block_map, TYPE_PARK_LOT); // need to do this later, after gen_tile_blocks()
			tile_to_block_map.clear(); // no longer needed
		}
		void add_streetlights() {
			streetlights.clear();
			streetlights.reserve(4*plots.size()); // one on each side of each plot
			float const b(STREETLIGHT_DIST_FROM_PLOT_EDGE), a(1.0 - b); // spacing from light pos to plot edge (placed just outside the plot, so spacing is negative)

			for (auto i = plots.begin(); i != plots.end(); ++i) {
				streetlights.emplace_back(point((a*i->x1() + b*i->x2()), (0.75*i->y1() + 0.25*i->y2()), i->z2()), -plus_x); // left   edge one   quarter  up
				streetlights.emplace_back(point((a*i->x2() + b*i->x1()), (0.25*i->y1() + 0.75*i->y2()), i->z2()),  plus_x); // right  edge three quarters up
				streetlights.emplace_back(point((0.25*i->x1() + 0.75*i->x2()), (a*i->y1() + b*i->y2()), i->z2()), -plus_y); // bottom edge three quarters right
				streetlights.emplace_back(point((0.75*i->x1() + 0.25*i->x2()), (a*i->y2() + b*i->y1()), i->z2()),  plus_y); // top    edge one   quarter  right
			}
			sort_streetlights_by_yx();
		}
		void get_road_bcubes(vect_cube_t &bcubes) const {
			get_all_bcubes(roads,  bcubes);
			get_all_bcubes(tracks, bcubes);
		}
		void get_plot_bcubes(vect_cube_with_zval_t &bcubes) const { // Note: z-values of cubes indicate building height ranges
			if (plots.empty()) return; // connector road city
			unsigned const start(bcubes.size());

			for (auto i = plots.begin(); i != plots.end(); ++i) {
				bcubes.push_back(*i); // capture all plot bcubes, even parks (needed for pedestrians)
				bcubes.back().is_park = i->is_park;
			}
			vector3d const city_radius(0.5*bcube.get_size());
			point const city_center(bcube.get_cube_center());

			for (auto i = bcubes.begin()+start; i != bcubes.end(); ++i) { // set zvals to control building height range, higher in city center
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
			point const query_pos(pos - get_camera_coord_space_xlate());
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
		void get_roads_sphere_coll(point const &pos, float radius, bool include_intersections, bool xy_only, vect_cube_t &out, vect_cube_t *out_bt) const {
			if (roads.empty()) return;
			vector3d const xlate(get_camera_coord_space_xlate());
			point const query_pos(pos - xlate);
			if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return;
			get_bcubes_sphere_coll(roads, out, query_pos, radius, xy_only, xlate);
				
			if (include_intersections) { // used for global road network
				for (unsigned i = 0; i < 3; ++i) { // {2-way, 3-way, 4-way}
					get_bcubes_sphere_coll(isecs[i], out, query_pos, radius, xy_only, xlate);
				}
			}
			if (out_bt) {
				get_bcubes_sphere_coll(bridges, *out_bt, query_pos, radius, xy_only, xlate);
				get_bcubes_sphere_coll(tunnels, *out_bt, query_pos, radius, xy_only, xlate);
			}
			get_bcubes_sphere_coll(tracks, out, query_pos, radius, xy_only, xlate);
		}
		bool proc_sphere_coll(point &pos, point const &p_last, float radius, float prev_frame_zval, vector3d *cnorm) const {
			vector3d const xlate(get_camera_coord_space_xlate());
			float const dist(p2p_dist(pos, p_last));
			if (!sphere_cube_intersect_xy(pos, (radius + dist), (bcube + xlate))) return 0;
			bool plot_coll(0);
			
			if (!plots.empty()) {
				float const max_obj_z(bcube.z1() + radius);
				if (pos.z < max_obj_z) {pos.z = max_obj_z; plot_coll = 1;} // make sure the sphere is above the city road/plot surface
			}
			for (unsigned n = 1; n < 3; ++n) { // intersections with stoplights (3-way, 4-way)
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
			if ((pos.z - xlate.z - radius < bcube.z2()) + (streetlight_ns::get_streetlight_height())) { // below the level of the streetlights
				if (proc_streetlight_sphere_coll(pos, radius, xlate, cnorm)) return 1;
			}
			if (city_obj_placer.proc_sphere_coll(pos, p_last, radius, cnorm)) return 1;
			
			if (0 && plot_coll) { // no other collisions - return collision with plot or road - doesn't work correctly for bouncing balls
				if (cnorm) {*cnorm = plus_z;}
				return 1;
			}
			return 0;
		}
		bool line_intersect(point const &p1, point const &p2, float &t) const { // Note: xlate has already been applied
			cube_t c(bcube); // deep copy
			c.z2() += stoplight_ns::stoplight_max_height();
			if (!c.line_intersects(p1, p2)) return 0;
			bool ret(0);

			for (unsigned n = 1; n < 3; ++n) { // intersections with stoplights (3-way, 4-way)
				for (auto i = isecs[n].begin(); i != isecs[n].end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
			}
			for (auto i = bridges.begin(); i != bridges.end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
			for (auto i = tunnels.begin(); i != tunnels.end(); ++i) {ret |= i->line_intersect(p1, p2, t);}
			ret |= line_intersect_streetlights(p1, p2, t);
			ret |= city_obj_placer.line_intersect(p1, p2, t);
			return ret;
		}
		bool check_mesh_disable(point const &pos, float radius) const {
			if (tunnels.empty()) return 0;
			point const query_pos(pos - get_camera_coord_space_xlate());
			cube_t query_region; query_region.set_from_sphere(query_pos, radius); // actually a cube, not a sphere

			for (auto i = tunnels.begin(); i != tunnels.end(); ++i) {
				if (i->check_mesh_disable(query_region)) return 1;
			}
			return 0;
		}
		bool tile_contains_tunnel(cube_t const &bcube) const {
			for (auto i = tunnels.begin(); i != tunnels.end(); ++i) {
				if (i->intersects_xy(bcube)) return 1;
			}
			return 0;
		}
		bool point_in_tunnel(point const &pos) const {
			for (auto i = tunnels.begin(); i != tunnels.end(); ++i) {
				if (i->contains_pt(pos)) return 1; // Note: checks z-val
			}
			return 0;
		}
		int get_color_at_xy(point const &pos, colorRGBA &color) const {
			// Note: query results are mutually exclusive since there's no overlap, so can early terminate on true
			if (!bcube.contains_pt_xy(pos)) return 0;
			
			for (auto i = bridges.begin(); i != bridges.end(); ++i) {
				if (i->contains_pt_xy_exp(pos, 1.0*city_params.road_width)) {color = WHITE; return INT_ROAD;}
			}
			for (auto i = tunnels.begin(); i != tunnels.end(); ++i) {
				if (i->contains_pt_xy(pos)) {color = BROWN; return INT_ROAD;}
			}
			if (!conn_roads.empty()) { // global_rn connector roads - use this vector because we only care about XY projection (not Z), and conn_roads is smaller than roads
				for (auto i = conn_roads.begin(); i != conn_roads.end(); ++i) {
					if (i->contains_pt_xy(pos)) {color = GRAY; return INT_ROAD;}
				}
			}
			else {
				for (auto i = roads.begin(); i != roads.end(); ++i) {
					if (i->contains_pt_xy(pos)) {color = GRAY; return INT_ROAD;}
				}
			}
			for (auto i = tracks.begin(); i != tracks.end(); ++i) {
				if (i->contains_pt_xy(pos)) {color = LT_BROWN; return INT_ROAD;} // counts as road intersection (for now)
			}
			if (plots.empty()) { // connector road
				for (auto i = isecs[0].begin(); i != isecs[0].end(); ++i) { // 2-way intersections
					if (i->contains_pt_xy(pos)) {color = GRAY; return INT_ROAD;}
				}
			}
			if (city_obj_placer.pt_in_parking_lot_xy(pos)) {color = DK_GRAY; return INT_PARKING;}
			if (city_obj_placer.get_color_at_xy(pos, color)) {return INT_PLOT;} // hit a detail object, but still in a plot
			
			if (!plots.empty()) { // inside a city and not over a road - must be over a plot or park

				for (auto i = parks.begin(); i != parks.end(); ++i) {
					if (i->contains_pt_xy(pos)) {color = GREEN; return INT_PARK;}
				}
				color = colorRGBA(0.65, 0.65, 0.65, 1.0);
				return INT_PLOT;
			}
			return INT_NONE;
		}
		bool cube_overlaps_road_xy(cube_t const &c) const {
			// can we use conn_roads here for global_rn?
			for (auto i = roads.begin(); i != roads.end(); ++i) {if (i->intersects(c)) return 1;}
			return 0;
		}
		bool cube_overlaps_parking_lot_xy(cube_t const &c) const {return city_obj_placer.cube_overlaps_parking_lot_xy(c);}

		void draw(road_draw_state_t &dstate, bool shadow_only, bool is_connector_road) {
			if (empty()) return;
			if (!dstate.check_cube_visible(bcube, 1.0, shadow_only)) return; // VFC/too far

			if (shadow_only) {
				if (!is_connector_road) { // connector road has no stoplights to cast shadows
					// Note: we can store the contents of qbd_sl in a VBO to avoid recreating it every frame for the shadow pass
					for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
						if (!dstate.check_cube_visible(b->bcube, 0.16, 1)) continue; // VFC/too far; dist_scale=0.16
						for (unsigned i = 1; i < 3; ++i) {dstate.draw_stoplights(isecs[i], b->ranges[TYPE_ISEC2 + i], 1);} // intersections with stoplights (3-way, 4-way)
					}
				}
			}
			else {
				for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
					if (!dstate.check_cube_visible(b->bcube)) continue; // VFC/too far
					dstate.begin_tile(b->bcube.get_cube_center());
					dstate.draw_road_region(segs,       b->ranges[TYPE_RSEG  ], b->quads[TYPE_RSEG  ], TYPE_RSEG  ); // road segments
					dstate.draw_road_region(plots,      b->ranges[TYPE_PLOT  ], b->quads[TYPE_PLOT  ], TYPE_PLOT  ); // plots
					dstate.draw_road_region(plots,      b->ranges[TYPE_PLOT  ], b->quads[TYPE_PARK  ], TYPE_PARK  ); // parks (stored as plots)
					dstate.draw_road_region(track_segs, b->ranges[TYPE_TRACKS], b->quads[TYPE_TRACKS], TYPE_TRACKS); // railroad tracks
					dstate.draw_road_region(city_obj_placer.parking_lots, b->ranges[TYPE_PARK_LOT], b->quads[TYPE_PARK_LOT], TYPE_PARK_LOT); // parking lots
				
					for (unsigned i = 0; i < 3; ++i) { // intersections (2-way, 3-way, 4-way)
						dstate.draw_road_region(isecs[i], b->ranges[TYPE_ISEC2 + i], b->quads[TYPE_ISEC2 + i], (TYPE_ISEC2 + i));
						if (i > 0) {dstate.draw_stoplights(isecs[i], b->ranges[TYPE_ISEC2 + i], 0);}
					}
				} // for b
			}
			draw_streetlights(dstate, shadow_only, 0);
			
			// draw bridges and tunnels; only in connector road network; bridgesand tunnels are sparse/uncommon, so don't need to be batched by blocks
			for (auto b = bridges.begin(); b != bridges.end(); ++b) {
				dstate.draw_bridge(*b, shadow_only);
				b->draw_streetlights(dstate, shadow_only, 0);
			}
			for (auto t = tunnels.begin(); t != tunnels.end(); ++t) {
				dstate.draw_tunnel(*t, shadow_only);
				t->draw_streetlights(dstate, shadow_only, 1); // always_on=1
			}
			city_obj_placer.draw_detail_objects(dstate, shadow_only);
		}
		void add_city_lights(vector3d const &xlate, cube_t &lights_bcube) const { // for now, the only light sources added by the road network are city block streetlights
			add_streetlight_dlights(xlate, lights_bcube, 0);
			for (auto b = bridges.begin(); b != bridges.end(); ++b) {b->add_streetlight_dlights(xlate, lights_bcube, 0);}
			for (auto t = tunnels.begin(); t != tunnels.end(); ++t) {t->add_streetlight_dlights(xlate, lights_bcube, 1);} // always_on=1
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
				car.dim   = seg.dim;
				car.dir   = rgen.rand_bool();
				car.max_speed = rgen.rand_uniform(0.66, 1.0); // add some speed variation
				car.cur_road  = seg.road_ix;
				car.cur_seg   = seg_ix;
				car.cur_road_type = TYPE_RSEG;
				vector3d car_sz(nom_car_size); // {length, width, height} // Note: car models should all be the same size
				car.height = car_sz.z;
				point pos;
				float val1(seg.d[seg.dim][0] + 0.6f*car_sz.x), val2(seg.d[seg.dim][1] - 0.6f*car_sz.x);
				if (val1 >= val2) continue; // failed, try again (connector road junction?)
				pos[!seg.dim]  = seg.get_center_dim(!seg.dim); // center of road
				pos[!seg.dim] += ((car.dir ^ car.dim) ? -1.0 : 1.0)*get_car_lane_offset(); // place in right lane
				pos[ seg.dim]  = rgen.rand_uniform(val1, val2); // place at random pos in segment
				pos.z = seg.z2() + 0.5*car_sz.z; // place above road surface
				if (seg.dim) {swap(car_sz.x, car_sz.y);}
				car.bcube.set_from_point(pos);
				car.bcube.expand_by(0.5*car_sz);
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
		unsigned decode_plot_id(unsigned global_plot_id) const {
			assert(global_plot_id >= plot_id_offset);
			unsigned const plot_id(global_plot_id - plot_id_offset);
			assert(plot_id < plots.size());
			return plot_id;
		}
		unsigned encode_plot_id(unsigned local_plot_id) const {return (local_plot_id + plot_id_offset);}
		cube_t      const &get_plot_from_global_id(unsigned global_plot_id) const {return plots         [decode_plot_id(global_plot_id)];}
		vect_cube_t const &get_colliders_for_plot (unsigned global_plot_id) const {return plot_colliders[decode_plot_id(global_plot_id)];}

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
			if (!select_building_in_plot(global_plot, rgen.rand(), building)) return 0; // no buildings in plot
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
					
					if (car.in_isect() && city_ix != CONN_CITY_IX) {
						road_network_t const &rn(road_networks[city_ix]);
						vector<road_isec_t> const &isecs(rn.isecs[car.get_isec_type()]); // must be a 3-way or 4-way intersection
						car.cur_city = city_ix;
						assert(car.cur_seg  < isecs.size());
						car.cur_road = isecs[car.cur_seg].rix_xy[2*(!car.dim) + 0]; // use the road in the other dim, since it must be within the new city (dir doesn't matter)
						assert(car.cur_road < rn.roads.size());
						car.entering_city = 1; // flag so that collision detection works
					}
				}
				assert(get_car_rn(car, road_networks, global_rn).get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
				return; // always within same city, no city update
			}
			road_isec_t const &isec(get_car_isec(car)); // conn_ix: {-x, +x, -y, +y}
			unsigned const orient(car.get_orient());
			int conn_ix(isec.conn_ix[orient]), rix(isec.rix_xy[car.get_orient()]);
			assert(isec.conn & (1<<orient));

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
			cube_t const bcube(get_road_bcube_for_car(car, global_rn));

			if (!bcube.intersects_xy(car.bcube)) { // sanity check
				cout << "bad intersection:" << endl << car.str() << endl << "bcube: " << bcube.str() << endl;
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
			if (car.cur_road_type != TYPE_ISEC2 && !get_car_isec(car).can_go_now(car)) return 0; // check stoplights and blocked intersections
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
			car.decelerate_fast();
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
	public:
		void update_car(car_t &car, rand_gen_t &rgen, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
			assert(car.cur_city == city_id);
			if (car.is_parked()) return; // stopped, no update (for now)
			
			if (car.in_isect()) {
				road_isec_t const &isec(get_car_isec(car));
				isec.notify_waiting_car(car); // even if not stopped

				// unclear why this was needed (how was stopped_at_light not set earlier?)
				if (isec.contains_pt_xy(car.get_front()) && !isec.contains_pt_xy(car.get_center()) && !car_can_fit_in_seg(car, global_rn)) { // not yet in the isec - stop and wait
					stop_and_wait_car(car, rgen, road_networks, isec);
					return;
				}
			}
			if (car.stopped_at_light) {
				bool const was_stopped(car.is_stopped());
				if (car_can_go_now(car, global_rn)) {car.stopped_at_light = 0;} // can go now
				else if (car.in_isect()) {stop_and_wait_car(car, rgen, road_networks, get_car_isec(car));} // Note: is_isect test allows cars to coast through lights when decel is very low
				if (was_stopped) return; // no update needed
			} else {car.maybe_accelerate();}

			cube_t const bcube(get_road_bcube_for_car(car));
			if (!bcube.intersects_xy(car.prev_bcube)) {cout << car.str() << endl << bcube.str() << endl; assert(0);} // sanity check
			bool const dim(car.dim);
			float const road_dz(bcube.dz());
			car.in_tunnel = point_in_tunnel(car.get_center());

			if (road_dz != 0.0) { // car on connector road
				assert(car.cur_road_type == TYPE_RSEG);
				assert(car.cur_city == CONN_CITY_IX);
				bool const slope(get_car_seg(car).slope);
				float const car_pos(car.bcube.get_center_dim(dim)); // center of car in dim
				float const road_len(bcube.get_sz_dim(dim));
				assert(road_len > TOLERANCE);
				float const t((car_pos - bcube.d[dim][0])/road_len); // car pos along road in (0.0, 1.0)
				float const road_z(bcube.d[2][slope] + t*(bcube.d[2][!slope] - bcube.d[2][slope]));
				float const car_len(car.get_length());
				car.dz = ((slope ^ car.dir) ? 1.0 : -1.0)*road_dz*(car_len/road_len);
				car.bcube.z1() = road_z - 0.5*fabs(car.dz);
				car.bcube.z2() = road_z + 0.5*fabs(car.dz) + car.height;
			}
			else if (car.dz != 0.0) { // car moving from connector road to level city
				float const road_z(bcube.d[2][1]);
				car.dz = 0.0;
				car.bcube.z1() = road_z;
				car.bcube.z2() = road_z + car.height;
			}
			if (car.turn_dir != TURN_NONE) {
				assert(car.in_isect());
				bool const turn_dir(car.turn_dir == TURN_RIGHT); // 0=left, 1=right
				point const car_center(car.get_center()), prev_center(car.prev_bcube.get_cube_center());
				float const car_lane_offset(get_car_lane_offset());
				float const trad_mult((car.cur_road_type == TYPE_ISEC2) ? 2.0 : 1.0); // larger turn radius for 2-way intersections (bends)
				float const turn_radius((turn_dir ? 0.15 : 0.25)*trad_mult*city_params.road_width); // right turn has smaller radius
				float const isec_center(bcube.get_cube_center()[dim]);
				float const centerline(isec_center + (((car.turn_dir == TURN_LEFT) ^ car.dir) ? -1.0 : 1.0)*car_lane_offset);
				float const prev_val(prev_center[dim]), cur_val(car_center[dim]);
				float const dist_to_turn(fabs(cur_val - centerline));

				if (dist_to_turn < turn_radius) { // turn radius; Note: cars turn around their center points, not their front wheels, which looks odd
					float const dist_from_turn_start(turn_radius - dist_to_turn);
					float const dev(turn_radius - sqrt(turn_radius*turn_radius - dist_from_turn_start*dist_from_turn_start));
					float const new_center(car.turn_val + dev*((turn_dir^car.dir^dim) ? 1.0 : -1.0));
					float const adj(new_center - car_center[!dim]);
					float const frame_dist(p2p_dist_xy(car_center, prev_center)); // total XY distance the car is allowed to move
					car.rot_z = (turn_dir ? -1.0 : 1.0)*(1.0 - CLIP_TO_01(dist_to_turn/turn_radius));
					car.bcube.d[!dim][0] += adj; car.bcube.d[!dim][1] += adj;
					vector3d const move_dir(car.get_center() - prev_center); // total movement from car + turn
					float const move_dist(move_dir.mag());
					
					if (move_dist > TOLERANCE) { // avoid division by zero
						vector3d const delta(move_dir*(frame_dist/move_dist - 1.0)); // overshoot value due to turn
						car.bcube += delta;
					}
				}
				if (min(prev_val, cur_val) <= centerline && max(prev_val, cur_val) > centerline) { // crossed the lane centerline boundary
					car.move_by(centerline - cur_val); // align to lane centerline
					vector3d const car_sz(car.bcube.get_size());
					float const size_adj(0.5f*(car_sz[dim] - car_sz[!dim]));
					vector3d expand(zero_vector);
					expand[dim] -= size_adj; expand[!dim] += size_adj;
					car.bcube.expand_by(expand); // fix aspect ratio
					if ((dim == 0) ^ (car.turn_dir == TURN_LEFT)) {car.dir ^= 1;}
					car.dim     ^= 1;
					car.rot_z    = 0.0;
					car.turn_val = 0.0; // reset
					car.turn_dir = TURN_NONE; // turn completed
					car.entering_city = 0;
					road_isec_t const &isec(get_car_isec(car));
					
					if (isec.conn_ix[car.get_orient()] >= 0) {
						short const rix(isec.rix_xy[car.get_orient()]);
						assert(rix >= 0); // not connector road
						car.cur_road = rix; // switch to using road_ix in new dim
					}
				}
			}
			if (bcube.contains_cube_xy(car.bcube)) { // in same road seg/int
				assert(get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
				return; // done
			}
			point const car_front(car.get_front(0.375)); // near the front, so that we can stop at the intersection

			// car crossing the border of this bcube, update state
			if (!bcube.contains_pt_xy_inc_low_edge(car_front)) { // move to another road seg/int
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
						point const dest_pos(car_rn.get_car_dest_isec_center(car, road_networks, global_rn));
						vector3d const dest_dir(dest_pos - car.get_center());
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
						} // for d
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
					car.stopped_at_light   = (isec.red_or_yellow_light(car) || !car_rn.car_can_go_now(car, global_rn));
					if (car.stopped_at_light) {car.decelerate_fast();}
					if (car.turn_dir != TURN_NONE) {car.turn_val = car.get_center()[!dim];} // capture car centerline before the turn
				}
			}
			assert(get_car_rn(car, road_networks, global_rn).get_road_bcube_for_car(car, global_rn).intersects_xy(car.bcube)); // sanity check
		}
	private:
		point get_car_dest_isec_center(car_t &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
			if (car.dest_city == city_id) {return get_isec_by_ix(car.dest_isec).get_cube_center();} // local destination within the current city
			assert(car.dest_city < road_networks.size());
			road_isec_t const *const isec(find_isec_to_dest_city(car, road_networks[car.dest_city], global_rn)); // destination in another city
			assert(isec != nullptr); // path must exist, otherwise this city wouldn't have been chosen
			return isec->get_cube_center();
		}
		road_isec_t const *find_isec_to_dest_city(car_t &car, road_network_t const &dest_rn, road_network_t const &global_rn) const {
			assert(car. cur_city == city_id);
			assert(car.dest_city == dest_rn.city_id);
			assert(dest_rn.city_id != city_id); // not ourself
			// Note: here we don't attempt to find shortcuts through other cities as this would be quite complex
			auto it(cix_to_isec.find(car.dest_city));
			if (it != cix_to_isec.end()) {return it->second;} // found
			return nullptr; // not found, caller can error check
		}
	public:
		bool choose_new_car_dest(car_t &car, rand_gen_t &rgen) const {
			unsigned const num_tot(isecs[0].size() + isecs[1].size() + isecs[2].size());
			if (num_tot == 0) return 0; // no isecs to select
			car.dest_isec = (unsigned short)(rgen.rand() % num_tot);
			return 1;
		}
		bool car_at_dest(car_t const &car) const {
			return get_isec_by_ix(car.dest_isec).contains_pt_xy(car.get_center());
		}
		road_isec_t const &get_isec_by_ix(unsigned ix) const {
			for (unsigned n = 0; n < 3; ++n) {
				unsigned const sz(isecs[n].size());
				if (ix < sz) return isecs[n][ix];
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
		}
		static road_network_t const &get_car_rn(car_base_t const &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) {
			if (car.cur_city == CONN_CITY_IX) return global_rn;
			assert(car.cur_city < road_networks.size());
			return road_networks[car.cur_city];
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
		bool mark_crosswalk_in_use(point const &pos, bool dim, bool dir) const {
			road_isec_t const *isec(find_isec_containing_pt(pos, 1, 3)); // 2-way can be skipped because there's no light/crosswalk
			if (isec == nullptr) return 0; // ped not at a crosswalk, maybe crossing in the middle of the street; this is okay for now, nothing else to do in this case
			isec->stoplight.mark_crosswalk_in_use(dim, dir);
			return 1;
		}
		bool check_isec_sphere_coll(point const &pos, float radius) const {
			road_isec_t const *isec(find_isec_containing_pt(pos, 1, 3)); // 2-way can be skipped because there's no light/crosswalk
			if (isec == nullptr) return 0;
			return isec->check_sphere_coll(pos, radius);
		}
		int get_nearby_road_ix(point const &pos, bool road_dim) const {
			for (auto r = roads.begin(); r != roads.end(); ++r) {
				if (r->dim != road_dim) continue;
				cube_t bcube(*r);
				bcube.expand_by_xy(0.01*city_params.road_width); // expand slightly to pick up roads that are adjacent to this point
				if (bcube.contains_pt_xy(pos)) {return (r - roads.begin());}
			}
			return -1; // should never get here, but occasionally can due to bad collisions between peds, floating-point error, etc.
		}
	}; // road_network_t

	vector<road_network_t> road_networks; // one per city
	road_network_t global_rn; // connects cities together; no plots
	road_draw_state_t dstate;
	rand_gen_t rgen;

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

public:
	bool empty() const {return road_networks.empty();}
	bool has_tunnels() const {return global_rn.has_tunnels();}
	bool point_in_tunnel(point const &pos) const {return global_rn.point_in_tunnel(pos);}

	road_network_t const &get_city(unsigned city_ix) const {
		if (city_ix == CONN_CITY_IX) {return global_rn;}
		assert(city_ix < road_networks.size());
		return road_networks[city_ix];
	}
	cube_t const &get_city_bcube(unsigned city_ix) const {return get_city(city_ix).get_bcube();}
	cube_t const &get_city_plot_bcube(unsigned city_ix, unsigned plot_ix) const {return get_city(city_ix).get_plot_bcube(plot_ix);}
	vect_cube_t const &get_colliders_for_plot(unsigned city_ix, unsigned global_plot_id) const {return get_city(city_ix).get_colliders_for_plot(global_plot_id);}

	cube_t get_city_bcube_for_cars(unsigned city_ix) const {
		cube_t bcube(get_city_bcube(city_ix));
		bcube.expand_by_xy(city_params.get_max_car_size().x); // expand by car length to fully include cars that are partially inside connector road intersections
		return bcube;
	}
	bool cube_overlaps_road_xy(cube_t const &c, unsigned city_ix) const {return get_city(city_ix).cube_overlaps_road_xy(c);}
	bool cube_overlaps_parking_lot_xy(cube_t const &c, unsigned city_ix) const {return get_city(city_ix).cube_overlaps_parking_lot_xy(c);}

	void gen_roads(cube_t const &region, float road_width, float road_spacing) {
		//timer_t timer("Gen Roads"); // ~0.5ms
		road_networks.push_back(road_network_t(region, road_networks.size()));
		if (!road_networks.back().gen_road_grid(road_width, road_spacing)) {road_networks.pop_back(); return;}
		//cout << "Roads: " << road_networks.back().num_roads() << endl;
	}
	bool connect_two_cities(unsigned city1, unsigned city2, vect_cube_t &blockers, heightmap_query_t &hq, float road_width) {
		assert(city1 < road_networks.size() && city2 < road_networks.size());
		assert(city1 != city2); // check for self reference
		//cout << "Connect city " << city1 << " and " << city2 << endl;
		road_network_t &rn1(road_networks[city1]), &rn2(road_networks[city2]);
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
							float const conn_pos(r->get_center_dim(d));
							float const cost(0.5*global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2,
								city1, city2, city1, city2, hq, road_width, conn_pos, !d, 1, (r12==0), (r12!=0))); // check_only=1; half cost (prefer over 3-way intersection)
							
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
						float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, city1, city2, city1, city2, hq, road_width, conn_pos, !d, 1, 0, 0)); // check_only=1
						if (cost >= 0.0 && (best_cost < 0.0 || cost < best_cost)) {best_conn_pos = conn_pos; best_cost = cost; is_4way1 = is_4way2 = 0;}
					}
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Single segment dim: << "d " << cost: " << best_cost << endl;
					float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2,
						city1, city2, city1, city2, hq, road_width, best_conn_pos, !d, 0, is_4way1, is_4way2)); // check_only=0; make change
					assert(cost >= 0.0);
					return 1;
				}
			}
		} // for d
		point const center1(bcube1.get_cube_center()), center2(bcube2.get_cube_center());
		bool const dx(center1.x < center2.x), dy(center1.y < center2.y);
		cube_t const bc[2] = {bcube1, bcube2};
		
		if ((bc[dx].x1() - bc[!dx].x1() > min_jog) && (bc[dy].y1() - bc[!dy].y1() > min_jog)) {
			// connect with two road segments using a jog: Note: assumes cities are all the same size
			bool const inv_dim(rgen.rand_bool());
			//cout << "Try connect using jog in dim " << inv_dim << endl;
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
							try_single_jog_conn_road(city1, city2, blockers, hq, road_width, fdim, xval, yval, 1, best_xval, best_yval, best_cost, best_int_cube);
						} // for r2
					} // for r1
				}
				else {
					for (unsigned n = 0; n < city_params.num_conn_tries; ++n) { // make up to num_tries attempts at connecting the cities with a single jog
						float xval(rgen_uniform(xmin, xmax, rgen)), yval(rgen_uniform(ymin, ymax, rgen));
						if (!fdim) {swap(xval, yval);}
						try_single_jog_conn_road(city1, city2, blockers, hq, road_width, fdim, xval, yval, 0, best_xval, best_yval, best_cost, best_int_cube);
					} // for n
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Double segment cost: " << best_cost << " " << TXT(best_xval) << TXT(best_yval) << TXT(fdim) << ", int_cube: " << best_int_cube.str() << endl;
					hq.flatten_region_to(best_int_cube, city_params.road_border); // do this first to improve flattening
					unsigned road_ix[2] = {};
					road_ix[ fdim] = global_rn.num_roads();
					float const cost1(global_rn.create_connector_road(bcube1, best_int_cube, blockers, &rn1, nullptr, city1,
						CONN_CITY_IX, city1, city2, hq, road_width, (fdim ? best_xval : best_yval),  fdim, 0, is_4way, 0)); // check_only=0
					assert(cost1 >= 0.0);
					flatten_op_t const fop(hq.last_flatten_op); // cache for reuse later during decrease_only pass
					road_ix[!fdim] = global_rn.num_roads();
					float const cost2(global_rn.create_connector_road(best_int_cube, bcube2, blockers, nullptr, &rn2,
						CONN_CITY_IX, city2, city1, city2, hq, road_width, (fdim ? best_yval : best_xval), !fdim, 0, 0, is_4way)); // check_only=0
					assert(cost2 >= 0.0);
					global_rn.create_connector_bend(best_int_cube, (dx ^ fdim), (dy ^ fdim), road_ix[0], road_ix[1]);
					// decrease_only=1; remove any dirt that the prev road added
					hq.flatten_sloped_region(fop.x1, fop.y1, fop.x2, fop.y2, fop.z1, fop.z2, fop.dim, fop.border, fop.skip_six, fop.skip_eix, 0, 1);
					hq.flatten_region_to(best_int_cube, city_params.road_border, 1); // one more pass to fix mesh that was raised above the intersection by a sloped road segment
					return 1;
				}
			} // for d
		}
		return 0;
	}
	private:
	void try_single_jog_conn_road(unsigned city1, unsigned city2, vect_cube_t &blockers, heightmap_query_t &hq, float road_width,
		bool fdim, float xval, float yval, bool is_4way, float &best_xval, float &best_yval, float &best_cost, cube_t &best_int_cube)
	{
		float const height(hq.get_height_at(xval, yval) + ROAD_HEIGHT), half_width(0.5*road_width);
		cube_t const int_cube(xval-half_width, xval+half_width, yval-half_width, yval+half_width, height, height); // the candidate intersection point
		if (has_bcube_int_xy(int_cube, blockers)) return; // bad intersection, fail
		road_network_t &rn1(road_networks[city1]), &rn2(road_networks[city2]);
		float const cost1(global_rn.create_connector_road(rn1.get_bcube(), int_cube, blockers, &rn1, nullptr, city1,
			CONN_CITY_IX, city1, city2, hq, road_width, (fdim ? xval : yval), fdim, 1, is_4way, 0)); // check_only=1
		if (cost1 < 0.0) return; // bad segment
		if (best_cost > 0.0 && cost1 > best_cost) return; // bound - early terminate
		float const cost2(global_rn.create_connector_road(int_cube, rn2.get_bcube(), blockers, nullptr, &rn2,
			CONN_CITY_IX, city2, city1, city2, hq, road_width, (fdim ? yval : xval), !fdim, 1, 0, is_4way)); // check_only=1
		if (cost2 < 0.0) return; // bad segment
		float const cost(cost1 + cost2);
		if (best_cost < 0.0 || cost < best_cost) {best_xval = xval; best_yval = yval; best_int_cube = int_cube; best_cost = cost;}
	}
	public:
	void connect_all_cities(float *heightmap, unsigned xsize, unsigned ysize, float road_width, float road_spacing) {
		if (road_width == 0.0 || road_spacing == 0.0) return; // no roads
		unsigned const num_cities(road_networks.size());
		if (num_cities < 2) return; // not cities to connect
		timer_t timer("Connect Cities");
		heightmap_query_t hq(heightmap, xsize, ysize);
		vector<unsigned> is_conn(num_cities, 0); // start with all cities unconnected (0=unconnected, 1=connected, 2=done/connect failed
		vector<pair<float, unsigned>> cands;
		vect_cube_t blockers; // existing cities and connector roads that we want to avoid intersecting
		// gather city blockers
		get_city_bcubes(blockers);
		expand_cubes_by_xy(blockers, road_spacing); // separate roads by at least this value
		// place railroad tracks before roads so that roads will reset the mesh height; need to fix this later
		cube_t const tracks_region(calc_cubes_bcube(blockers));
		global_rn.gen_railroad_tracks(TRACKS_WIDTH*city_params.road_width, city_params.num_rr_tracks, tracks_region, blockers, hq);

		// full cross-product connectivity
		for (unsigned i = 0; i < num_cities; ++i) {
			for (unsigned j = i+1; j < num_cities; ++j) {
				bool const success(connect_two_cities(i, j, blockers, hq, road_width));
				//cout << "Trying to connect city " << i << " to city " << j << ": " << success << endl;
				if (!success) continue;
				//cout << i << " connected to " << j << endl;
				road_networks[i].register_connected_city(j);
				road_networks[j].register_connected_city(i);
			} // for j
		} // for i
		assign_city_clusters();
		global_rn.calc_bcube_from_roads();
		global_rn.split_connector_roads(road_spacing);
		global_rn.finalize_bridges_and_tunnels();
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
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->gen_parking_lots_and_place_objects(cars, have_cars);}
	}
	void get_city_bcubes(vect_cube_t &bcubes) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {bcubes.push_back(r->get_bcube());}
	}
	void get_all_road_bcubes(vect_cube_t &bcubes, bool connector_only) const {
		global_rn.get_road_bcubes(bcubes); // not sure if this should be included
		if (connector_only) return;
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_road_bcubes(bcubes);}
	}
	void get_all_plot_bcubes(vect_cube_with_zval_t &bcubes) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_plot_bcubes(bcubes);}
	}
	bool check_road_sphere_coll(point const &pos, float radius, bool xy_only, bool exclude_bridges_and_tunnels) const {
		return global_rn.check_road_sphere_coll(pos, radius, 1, xy_only, exclude_bridges_and_tunnels);
	}
	void get_roads_sphere_coll(point const &pos, float radius, bool include_intersections, bool xy_only, vect_cube_t &out, vect_cube_t *out_bt) const {
		global_rn.get_roads_sphere_coll(pos, radius, 1, xy_only, out, out_bt);
	}
	bool check_mesh_disable(point const &pos, float radius) const {return global_rn.check_mesh_disable(pos, radius);}
	bool tile_contains_tunnel(cube_t const &bcube) const {return global_rn.tile_contains_tunnel(bcube);}

	int get_color_at_xy(point const &pos, colorRGBA &color) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {
			int const ret(r->get_color_at_xy(pos, color));
			if (ret) return ret;
		}
		return global_rn.get_color_at_xy(pos, color);
	}
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, float prev_frame_zval, vector3d *cnorm) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {
			if (r->proc_sphere_coll(pos, p_last, radius, prev_frame_zval, cnorm)) return 1;
		}
		return global_rn.proc_sphere_coll(pos, p_last, radius, prev_frame_zval, cnorm); // needed for bridges and tunnels
	}
	bool line_intersect(point const &p1, point const &p2, float &t) const {
		bool ret(0);
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {ret |= r->line_intersect(p1, p2, t);}
		ret |= global_rn.line_intersect(p1, p2, t); // bridges and tunnels
		return ret;
	}
	void add_city_lights(vector3d const &xlate, cube_t &lights_bcube) const {
		global_rn.add_city_lights(xlate, lights_bcube); // no streetlights, but may need to add lights for bridges and tunnels
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->add_city_lights(xlate, lights_bcube);}
	}
	void draw(int trans_op_mask, vector3d const &xlate, bool use_dlights, bool shadow_only) { // non-const because qbd is modified
		if (road_networks.empty() && global_rn.empty()) return;

		if (trans_op_mask & 1) { // opaque pass, should be first
			//highres_timer_t timer(shadow_only ? "Draw City Shadows" : "Draw City"); // 1.1ms / 0.42ms shadows
			fgPushMatrix();
			translate_to(xlate);
			glDepthFunc(GL_LEQUAL); // helps prevent Z-fighting
			dstate.pre_draw(xlate, use_dlights, shadow_only, 1); // always_setup_shader=1
			assert(dstate.s.is_setup());
			for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->draw(dstate, shadow_only, 0);}
			global_rn.draw(dstate, shadow_only, 1); // connector road may have bridges, and therefore needs shadows
			dstate.post_draw();
			glDepthFunc(GL_LESS);
			fgPopMatrix();
		}
		if (trans_op_mask & 2) {dstate.draw_and_clear_light_flares();} // transparent pass; must be done last for alpha blending, and no translate
	}
	void draw_label() {dstate.show_label_text();}

	// cars/peds
	void next_frame() {
		if (!animate2) return;
		//timer_t timer("Update Stoplights");
		//for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {cout << r->get_traffic_density() << " ";} cout << endl;
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->next_frame();}
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
	cube_t const &get_plot_from_global_id(unsigned city_id, unsigned global_plot_id) const {return get_city(city_id).get_plot_from_global_id(global_plot_id);}
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
			auto const &conn(road_networks[car.cur_city].get_connected());
			float const new_city_prob(city_params.new_city_prob*min(0.4f, 0.1f*conn.size())); // 10% to 40% chance, depending on the number of connecting cities (to reduce traffic congestion)

			if (rgen.rand_float() < new_city_prob) { // select a different city when there are multiple cities
				if (rgen.rand_float() < city_params.traffic_balance_val) { // choose the connected city with the lowest traffic density
					float min_td(0.0);

					for (auto c = conn.begin(); c != conn.end(); ++c) {
						assert(*c < road_networks.size()); // excludes global_rn
						float const td(road_networks[*c].get_traffic_density());
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
		//cout << TXT(car.dest_valid) << TXT(car.dest_city) << TXT(car.dest_isec) << endl;
	}
	bool car_at_dest(car_t const &car) const {
		if (!car.dest_valid) return 0;
		return get_city(car.dest_city).car_at_dest(car);
	}
	road_network_t const &get_car_rn(car_base_t const &car) const {return road_network_t::get_car_rn(car, road_networks, global_rn);}
	
	void update_car(car_t &car, rand_gen_t &rgen) const {
		if (car.cur_city == NO_CITY_IX) return; // not in a city (in a garage), nothing to update
		//update_car_seg_stats(car); // not needed - stats not yet used
		get_car_rn(car).update_car(car, rgen, road_networks, global_rn);
		if (city_params.enable_car_path_finding) {update_car_dest(car);}
	}
	void update_car_seg_stats(car_base_t const &car) const {get_car_rn(car).update_car_seg_stats(car);}
	road_isec_t const &get_car_isec(car_base_t const &car) const {return get_car_rn(car).get_car_isec(car);}
	cube_t get_road_bcube_for_car(car_base_t const &car) const {return get_car_rn(car).get_road_bcube_for_car(car);}
	virtual cube_t get_bcube_for_car(car_base_t const &car) const {return get_road_bcube_for_car(car);}
}; // city_road_gen_t


// Note: the car_manager_t member functions that use road_gen are here rather than in cars.cpp
cube_t car_manager_t::get_cb_bcube(car_block_t const &cb ) const {
	if (cb.is_in_building()) {return garages_bcube;}
	return road_gen.get_city_bcube_for_cars(cb.cur_city);
}
road_isec_t const &car_manager_t::get_car_isec(car_t const &car) const {return road_gen.get_car_isec(car);}
bool car_manager_t::check_collision(car_t &c1, car_t &c2)        const {return c1.check_collision(c2, road_gen);}
void car_manager_t::register_car_at_city(car_t const &car) {road_gen.register_car_at_city(car.cur_city);}

void car_manager_t::add_car() {
	car_t car;
	if (road_gen.add_car(car, rgen)) {cars.push_back(car);}
}

void car_manager_t::update_cars() {
	for (auto i = cars.begin(); i != cars.end(); ++i) {road_gen.update_car(*i, rgen);} // run update logic
}

void car_manager_t::get_car_ix_range_for_cube(vector<car_block_t>::const_iterator cb, cube_t const &bc, unsigned &start, unsigned &end) const {
	start = cb->start; end = (cb+1)->start;
	assert(end <= cars.size());
	if (cb->is_in_building()) return; // cars parked in garages - keep full start/end range
	if (!road_gen.cube_overlaps_parking_lot_xy(bc, cb->cur_city)) {end   = cb->first_parked;} // moving cars only (beginning of range)
	if (!road_gen.cube_overlaps_road_xy       (bc, cb->cur_city)) {start = cb->first_parked;} // parked cars only (end of range)
	assert(start <= end);
}


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
cube_t const &ped_manager_t::get_city_plot_bcube_for_peds(unsigned city_ix, unsigned plot_ix) const {return road_gen.get_plot_from_global_id(city_ix, plot_ix);}
road_isec_t const &ped_manager_t::get_car_isec(car_base_t const &car) const {return road_gen.get_car_isec(car);}

cube_t ped_manager_t::get_expanded_city_bcube_for_peds(unsigned city_ix) const {
	cube_t bcube(road_gen.get_city_bcube(city_ix));
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
bool ped_manager_t::check_isec_sphere_coll(pedestrian_t const &ped) const {
	return road_gen.get_city(ped.city).check_isec_sphere_coll(ped.pos, 0.6*ped.radius); // Note: no xlate is required since peds and city are in the same coord space
}
bool ped_manager_t::check_streetlight_sphere_coll(pedestrian_t const &ped) const {
	return road_gen.get_city(ped.city).check_streetlight_sphere_coll_xy(ped.pos, ped.radius);
}
int ped_manager_t::get_road_ix_for_ped_crossing(pedestrian_t const &ped, bool road_dim) const { // returns -1 on failure (ped not in the road)
	return road_gen.get_city(ped.city).get_nearby_road_ix(ped.pos, road_dim);
}

// path finding
bool ped_manager_t::choose_dest_building_or_parked_car(pedestrian_t &ped) { // modifies rgen, non-const
	ped.has_dest_bldg = ped.has_dest_car = ped.at_dest = 0; // will choose a new dest

	if (city_params.num_cars == 0 || (rgen.rand() & 3) != 0) { // choose a dest building 75% of the time, 100% of the time if there are no cars
		ped.has_dest_bldg = road_gen.choose_dest_building(ped.city, ped.dest_plot, ped.dest_bldg, rgen);
	}
	if (!ped.has_dest_bldg) { // chose a dest parked car 25% of the time, or if choosing a dest building failed
		ped.has_dest_car = choose_dest_parked_car(ped.city, ped.dest_plot, ped.dest_bldg, ped.dest_car_center);
		if (!ped.has_dest_car) return 0;
		ped.dest_plot = road_gen.get_city(ped.city).encode_plot_id(ped.dest_plot);
	}
	ped.next_plot = get_next_plot(ped);
	return 1;
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
	//cout << "dlights: " << lights.size() << ", bcube: " << lights_bcube.str() << endl;

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
	add_dynamic_lights_city(lights_bcube, dlight_add_thresh);
	upload_dlights_textures(lights_bcube, dlight_add_thresh);
	prev_had_lights = !lights.empty();
}

void city_lights_manager_t::setup_shadow_maps(vector<light_source> &light_sources, point const &cpos, unsigned max_smaps) {
	unsigned const num_smaps(min((unsigned)light_sources.size(), min(max_smaps, MAX_DLIGHT_SMAPS)));
	dl_smap_enabled = 0;
	if (!enable_dlight_shadows || shadow_map_sz == 0 || num_smaps == 0) return;
	sort_lights_by_dist_size(light_sources, cpos); // Note: may already be sorted for enabled lights selection, but okay to sort again
	cmp_light_source_sz_dist sz_cmp(cpos);
	unsigned num_used(0);
	unsigned const smap_size(city_params.smap_size); // 0 = use default shadow map resolution
	// capture player pos in global coordinate space before replacing with light pos so it can be used for LOD during model drawing
	pre_smap_player_pos = get_camera_pos() - get_camera_coord_space_xlate();
	// Note: if using a dynamic (distance-based) sm_size, need to maintain a pool of different sm resolutions somehow
	check_gl_error(430);

	// Note: slow to recreate shadow maps every frame, but most lights are either dynamic (headlights) or include dynamic shadow casters (cars) and need to be updated every frame anyway
	// Do we want to gradually fade in new shadow maps and fade out old ones? But how do we track which lights are associated with old shadow maps?
	// Tracking positions won't work for car headlights because they move. We don't have object pointers to track either. And what about lights that are no longer in our list?
	for (auto i = light_sources.begin(); i != light_sources.end() && num_used < num_smaps; ++i) {
		if (i->has_no_shadows()) continue; // shadows not enabled for this light
		if (!i->is_very_directional()) continue; // not a spotlight
		if (sz_cmp.get_value(*i) < 0.002) break; // light influence is too low, skip even though we have enough shadow maps; can break because sort means all later lights also fail
		dl_smap_enabled |= i->setup_shadow_map(CITY_LIGHT_FALLOFF, 0, 0, 0, smap_size);
		++num_used;
	} // for i
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
			params.city_border, params.slope_width, params.num_samples, x1, y1, x2, y2)) return 0;
		float const elevation(flatten_region(x1, y1, x2, y2, params.slope_width));
		cube_t const pos_range(add_plot(x1, y1, x2, y2, elevation));
		if (cities_bcube.is_all_zeros()) {cities_bcube = pos_range;} else {cities_bcube.union_with_cube(pos_range);}
		if (params.roads_enabled()) {road_gen.gen_roads(pos_range, params.road_width, params.road_spacing);}
		return 1;
	}
	void gen_cities(city_params_t const &params) {
		if (params.num_cities == 0) return;
		cube_t cities_bcube(all_zeros);
		{ // open a scope
			timer_t t("Choose City Location");
			for (unsigned n = 0; n < params.num_cities; ++n) {gen_city(params, cities_bcube);}
		}
		if (!cities_bcube.is_all_zeros()) {set_buildings_pos_range(cities_bcube);}
		road_gen.connect_all_cities(heightmap, xsize, ysize, params.road_width, params.road_spacing);
		road_gen.add_streetlights();
		road_gen.gen_tile_blocks();
		car_manager.init_cars(city_params.num_cars);
	}
	void gen_details() {
		if (road_gen.empty()) return; // nothing to do - no roads or cars
		// generate parking lots
		vector<car_t> parked_cars;
		vect_cube_t garages, hp_locs;
		bool const have_cars(!car_manager.empty());
		road_gen.gen_parking_lots_and_place_objects(parked_cars, have_cars);
		if (have_cars) {get_all_garages(garages);}
		if (city_params.has_helicopter_model()) {get_all_city_helipads(hp_locs);}
		car_manager.add_parked_cars(parked_cars, garages);
		car_manager.finalize_cars();
		car_manager.add_helicopters(hp_locs);
		ped_manager.init(city_params.num_peds, city_params.num_building_peds); // must be after buildings are placed
	}
	void get_city_bcubes(vect_cube_t &bcubes) const {return road_gen.get_city_bcubes(bcubes);}
	void get_all_road_bcubes(vect_cube_t &bcubes, bool connector_only) const {road_gen.get_all_road_bcubes(bcubes, connector_only);}
	void get_all_plot_bcubes(vect_cube_with_zval_t &bcubes) {road_gen.get_all_plot_bcubes(bcubes);} // caches plot_id_offset, so non-const

	// return: 0=no coll, 1=plot coll, 2=road coll, 3=both plot and road coll
	unsigned check_city_sphere_coll(point const &pos, float radius, bool xy_only, bool exclude_bridges_and_tunnels, bool ret_first_coll, unsigned check_mask) const {
		int ret(0);
		if ((check_mask & 1) && check_plot_sphere_coll(pos, radius, xy_only)) {ret |= 1;}
		if (ret_first_coll && ret) return ret;
		if ((check_mask & 2) && road_gen.check_road_sphere_coll(pos, radius, xy_only, exclude_bridges_and_tunnels)) {ret |= 2;}
		return ret;
	}
	void get_sphere_coll_cubes(point const &pos, float radius, bool include_intersections, bool xy_only, vect_cube_t &out, vect_cube_t *out_bt) const {
		get_plots_sphere_coll(pos, radius, xy_only, out);
		road_gen.get_roads_sphere_coll(pos, radius, include_intersections, xy_only, out, out_bt);
		get_driveway_sphere_coll_cubes(pos, radius, xy_only, out);
	}
	bool proc_city_sphere_coll(point &pos, point const &p_last, float radius, float prev_frame_zval, bool inc_cars, vector3d *cnorm) const {
		if (road_gen.proc_sphere_coll(pos, p_last, radius, prev_frame_zval, cnorm)) return 1;
		if (!inc_cars) return 0;
		return car_manager.proc_sphere_coll(pos, p_last, radius, cnorm); // Note: doesn't really work well, at least for player collisions
	}
	bool line_intersect(point const &p1, point const &p2, float &t) const {
		vector3d const xlate(get_camera_coord_space_xlate()), p1x(p1 - xlate), p2x(p2 - xlate);
		bool ret(road_gen.line_intersect(p1x, p2x, t));
		ret |= car_manager.line_intersect_cars(p1x, p2x, t);
		ret |= ped_manager.line_intersect_peds(p1x, p2x, t);
		return ret;
	}
	bool check_mesh_disable(point const &pos, float radius ) const {return road_gen.check_mesh_disable(pos, radius);}
	bool tile_contains_tunnel(cube_t const &bcube) const {return road_gen.tile_contains_tunnel(bcube);}

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
		if (!shadow_only && !reflection_pass && (trans_op_mask & 1)) {setup_city_lights(xlate);} // setup lights on first (opaque) non-shadow pass
		bool const use_dlights(enable_lights()), is_dlight_shadows(shadow_only == 2);
		if (reflection_pass == 0) {road_gen.draw(trans_op_mask, xlate, use_dlights, (shadow_only != 0));} // roads don't cast shadows and aren't reflected in water, but stoplights cast shadows
		car_manager.draw(trans_op_mask, xlate, use_dlights, (shadow_only != 0), is_dlight_shadows, 0);
		if (trans_op_mask & 1) {ped_manager.draw(xlate, use_dlights, (shadow_only != 0), is_dlight_shadows);} // opaque
		if ((trans_op_mask & 1) && !shadow_only) {road_gen.draw_label();} // after drawing cars so that it's in front
		// Note: buildings are drawn through draw_buildings()
	}
	void draw_roads(int trans_op_mask, vector3d const &xlate) {road_gen.draw(trans_op_mask, xlate, enable_lights(), 0);} // shadow_only=0
	void draw_cars_in_garages(vector3d const &xlate, bool shadow_only) {car_manager.draw(1, xlate, 1, shadow_only, 0, 1);} // opaque + garages pass
	void draw_peds_in_building(int first_ped_ix, ped_draw_vars_t const &pdv) {ped_manager.draw_peds_in_building(first_ped_ix, pdv);}
	void get_ped_bcubes_for_building(int first_ped_ix, unsigned bix, vect_cube_t &bcubes, bool moving_only) const {ped_manager.get_ped_bcubes_for_building(first_ped_ix, bix, bcubes, moving_only);}
	void register_person_hit(unsigned person_ix, room_object_t const &obj, vector3d const &velocity) {ped_manager.register_person_hit(person_ix, obj, velocity);}
	void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only) {ped_manager.draw_player_model(s, xlate, shadow_only);}

	void setup_city_lights(vector3d const &xlate) {
		if (world_mode != WMODE_INF_TERRAIN) return; // TT only
		if (prev_city_lights_setup_frame == cur_display_iter) return; // already called this frame
		prev_city_lights_setup_frame = cur_display_iter;
		//timer_t timer("City Dlights Setup");
		float const light_radius(1.0*light_radius_scale*get_tile_smap_dist()); // distance from the camera where headlights and streetlights are drawn
		if (!begin_lights_setup(xlate, light_radius, dl_sources)) return;
		car_manager.add_car_headlights(xlate, lights_bcube);
		road_gen.add_city_lights(xlate, lights_bcube);
		if (flashlight_on && !camera_in_building) {add_player_flashlight(0.25);} // add player flashlight
		clamp_to_max_lights(xlate, dl_sources);
		setup_shadow_maps(dl_sources, (camera_pdu.pos - xlate), city_params.max_shadow_maps);
		finalize_lights(dl_sources);
	}
	virtual bool enable_lights() const {return (is_night(max(STREETLIGHT_ON_RAND, HEADLIGHT_ON_RAND)) || road_gen.has_tunnels() || flashlight_on);} // only have lights at night
	void next_ped_animation() {ped_manager.next_animation();}
	void free_context() {car_manager.free_context(); ped_manager.free_context();}
	unsigned get_model_gpu_mem() const {return (ped_manager.get_model_gpu_mem() + car_manager.get_model_gpu_mem());}
}; // city_gen_t

city_gen_t city_gen;


bool parse_city_option(FILE *fp) {return city_params.read_option(fp);}
bool have_cities() {return city_params.enabled();}
// Note: this is used for parallel car/pedestrian updates and does not include city_params.num_building_peds
bool have_city_models() {
	return ((have_cities() && (city_params.num_cars > 0 || city_params.num_peds > 0)) || (enable_building_people_ai() && city_params.num_building_peds > 0));
}
float get_road_max_len   () {return city_params.road_spacing;}
float get_road_max_width () {return city_params.road_width;}
float get_min_obj_spacing() {return 4.0*ped_manager_t::get_ped_radius();} // allow a ped to walk between objects (two side-by-side)

void gen_cities(float *heightmap, unsigned xsize, unsigned ysize) {
	if (!have_cities()) return; // nothing to do
	city_gen.init(heightmap, xsize, ysize); // only need to call once for any given heightmap
	city_gen.gen_cities(city_params);
}
void gen_city_details() {city_gen.gen_details();} // called after gen_buildings()
void get_city_bcubes(vect_cube_t &bcubes) {city_gen.get_city_bcubes(bcubes);}
void get_city_road_bcubes(vect_cube_t &bcubes, bool connector_only) {city_gen.get_all_road_bcubes(bcubes, connector_only);}
void get_city_plot_bcubes(vector<cube_with_zval_t> &bcubes) {city_gen.get_all_plot_bcubes(bcubes);}
void next_city_frame(bool use_threads_2_3) {city_gen.next_frame(use_threads_2_3);}
void draw_cities(int shadow_only, int reflection_pass, int trans_op_mask, vector3d const &xlate) {city_gen.draw(shadow_only, reflection_pass, trans_op_mask, xlate);}
void draw_city_roads(int trans_op_mask, vector3d const &xlate) {city_gen.draw_roads(trans_op_mask, xlate);}
void setup_city_lights(vector3d const &xlate) {city_gen.setup_city_lights(xlate);}

void draw_peds_in_building(int first_ped_ix, ped_draw_vars_t const &pdv) {city_gen.draw_peds_in_building(first_ped_ix, pdv);}
void get_ped_bcubes_for_building(int first_ped_ix, unsigned bix, vect_cube_t &bcubes, bool moving_only) {city_gen.get_ped_bcubes_for_building(first_ped_ix, bix, bcubes, moving_only);}
void register_person_hit(unsigned person_ix, room_object_t const &obj, vector3d const &velocity) {city_gen.register_person_hit(person_ix, obj, velocity);}
void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only) {city_gen.draw_player_model(s, xlate, shadow_only);}

unsigned check_city_sphere_coll(point const &pos, float radius, bool exclude_bridges_and_tunnels, bool ret_first_coll, unsigned check_mask) {
	if (!have_cities()) return 0;
	return city_gen.check_city_sphere_coll((pos + get_tt_xlate_val()), radius, 1, exclude_bridges_and_tunnels, ret_first_coll, check_mask); // apply xlate for all static objects
}
void get_city_sphere_coll_cubes(point const &pos, float radius, bool include_intersections, bool xy_only, vect_cube_t &out, vect_cube_t *out_bt) {
	city_gen.get_sphere_coll_cubes((pos + get_tt_xlate_val()), radius, include_intersections, xy_only, out, out_bt);
}
bool proc_city_sphere_coll(point &pos, point const &p_last, float radius, float prev_frame_zval, bool xy_only, bool inc_cars, vector3d *cnorm, bool check_interior) {
	if (proc_buildings_sphere_coll(pos, p_last, radius, xy_only, cnorm, check_interior)) return 1;
	return city_gen.proc_city_sphere_coll(pos, p_last, radius, prev_frame_zval, inc_cars, cnorm); // Note: no xy_only for cities
}
bool line_intersect_city(point const &p1, point const &p2, float &t, bool ret_any_pt) {
	unsigned hit_bix(0); // unused
	bool ret(check_buildings_line_coll(p1, p2, t, hit_bix, 0, ret_any_pt)); // apply_tt_xlate=0
	ret |= city_gen.line_intersect(p1, p2, t);
	return ret;
}
bool line_intersect_city(point const &p1, point const &p2, point &p_int) {
	float t(1.0);
	if (!line_intersect_city(p1, p2, t)) return 0;
	p_int = p1 + t*(p2 - p1);
	return 1;
}
bool check_valid_scenery_pos(point const &pos, float radius, bool is_tall) {
	if (check_buildings_sphere_coll(pos, radius, 1, 1, 0, 1)) return 0; // apply_tt_xlate=1, xy_only=1, check_interior=0, exclude_city=1 (since we're checking plots below)
	if (check_city_sphere_coll(pos, radius, !is_tall))        return 0; // exclude bridges if not tall
	if (check_mesh_disable(pos, (radius + 2.0*HALF_DXY)))     return 0; // check tunnels
	return 1;
}
bool check_mesh_disable(point const &pos, float radius) {
	if (!have_cities()) return 0;
	return city_gen.check_mesh_disable((pos + get_tt_xlate_val()), radius); // apply xlate for all static objects
}
bool tile_contains_tunnel(cube_t const &bcube) {return city_gen.tile_contains_tunnel(bcube + get_tt_xlate_val());}
void destroy_city_in_radius(point const &pos, float radius) {city_gen.destroy_in_radius(pos, radius);}
bool get_city_color_at_xy(float x, float y, colorRGBA &color) {return city_gen.get_color_at_xy(x, y, color);}
cube_t get_city_lights_bcube() {return city_gen.get_lights_bcube();}
unsigned get_city_model_gpu_mem() {return city_gen.get_model_gpu_mem();}
void next_pedestrian_animation() {city_gen.next_ped_animation();}
void free_city_context() {city_gen.free_context();}
bool has_city_trees() {return (city_params.max_trees_per_plot > 0);}
vector3d get_nom_car_size() {return city_params.get_nom_car_size();}
void draw_cars_in_garages(vector3d const &xlate, bool shadow_only) {city_gen.draw_cars_in_garages(xlate, shadow_only);}

