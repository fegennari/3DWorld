// 3D World - City Object Class Definitions
// by Frank Gennari
// 08/14/21

#include "city_objects.h"
#include "lightmap.h" // for light_source

extern bool player_in_walkway;
extern int animate2, display_mode;
extern float fticks;
extern double camera_zh;
extern city_params_t city_params;
extern object_model_loader_t building_obj_model_loader;
extern building_params_t global_building_params;

unsigned const q2t_ixs[6] = {0,2,1,0,3,2}; // quad => 2 tris

bool do_line_clip_xy(point &v1, point &v2, float const d[3][2]);

float get_power_pole_offset() {return 0.045*city_params.road_width;}
unsigned get_building_models_gpu_mem() {return building_obj_model_loader.get_gpu_mem();}
int get_solarp_tid();


void textured_mat_t::pre_draw(bool shadow_only) {
	if (shadow_only && !has_alpha_mask) return; // not textured

	if (tid < 0 && !tex_name.empty()) {
		unsigned const ncolors(has_alpha_mask ? 4 : 3);
		int const use_mipmaps (has_alpha_mask ? 0 : 1); // disable mipmaps if this texture has an alpha mask because we need to keep binary alpha
		tid = get_texture_by_name(tex_name, 0, 0, 1, 4.0, 1, use_mipmaps, ncolors); // load/lookup texture if needed, 4.0 aniso
	}
	if (nm_tid < 0 && !nm_tex_name.empty()) {nm_tid = get_texture_by_name(nm_tex_name, 1);} // load/lookup texture if needed
	select_texture(tid);
	if (nm_tid >= 0) {select_texture(nm_tid, 5);} // bind normal map if it was specified
}
void textured_mat_t::post_draw(bool shadow_only) {
	if (!shadow_only && nm_tid >= 0) {bind_default_flat_normal_map();} // restore default flat normal map
}


void city_obj_t::set_bcube_from_vcylin(point const &base, float height, float xy_radius) {
	bcube.set_from_point(base);
	bcube.expand_by_xy(xy_radius);
	bcube.z2() += height;
}
bool city_obj_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	return sphere_cube_int_update_pos(pos_, radius_, (bcube + xlate), p_last, 0, cnorm);
}

vector3d oriented_city_obj_t::get_orient_dir() const {
	vector3d orient(zero_vector);
	orient[dim] = (dir ? 1.0 : -1.0); // here dim is assumed to be 0=x or 1=y; 2=z has so far not been needed
	return orient;
}

bool sphere_city_obj_cylin_coll(point const &cpos, float cradius, point &spos, point const &p_last, float sradius, point const &xlate, vector3d *cnorm) {
	point const pos2(cpos + xlate);
	float const r_sum(cradius + sradius);
	if (!dist_less_than(spos, pos2, r_sum)) return 0; // use sphere/vert cylinder instead?
	// since this is a cylinder, and we're not supposed to stand on top of it, assume collision normal is in the XY plane
	vector3d const coll_norm(vector3d((spos.x - pos2.x), (spos.y - pos2.y), 0.0).get_norm());
	spos += coll_norm*(r_sum - p2p_dist(spos, pos2)); // move away from pos2
	if (cnorm) {*cnorm = coll_norm;}
	return 1;
}

void disable_hemi_lighting_pre_post(draw_state_t &dstate, bool shadow_only, bool is_post) {
	if (!shadow_only) {dstate.s.add_uniform_float("hemi_lighting_scale", (is_post ? 0.5 : 0.0));} // disable or restore
}

// model_city_obj_t

// can't call get_model_id() virtual, must pass model_id in
model_city_obj_t::model_city_obj_t(point const &pos_, float height, bool dim_, bool dir_, unsigned model_id, bool is_cylinder_) :
	oriented_city_obj_t(pos_, 0.5*height, dim_, dir_), is_cylinder(is_cylinder_) // radius = 0.5*height
{
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
	vector3d expand;
	expand[ dim] = height*sz.x/sz.z; // depth
	expand[!dim] = height*sz.y/sz.z; // width
	expand.z = height;
	pos.z += 0.5*height; // pos is on the ground, while we want the bsphere to be at the center
	bcube.set_from_point(pos);
	bcube.expand_by(0.5*expand);
	set_bsphere_from_bcube(); // recompute bsphere from bcube
}
void model_city_obj_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only, animation_state_t *anim_state) const {
	if (!dstate.is_visible_and_unoccluded(bcube, dist_scale)) return;
	if (min_alpha > 0.0) {dstate.s.add_uniform_float("min_alpha", 0.9);}
	building_obj_model_loader.draw_model(dstate.s, pos, bcube, get_orient_dir(), color, dstate.xlate, get_model_id(), shadow_only, 0, anim_state);
	if (min_alpha > 0.0) {dstate.s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA);} // restore to the default
}
bool model_city_obj_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (!is_cylinder) {return oriented_city_obj_t::proc_sphere_coll(pos_, p_last, radius_, xlate, cnorm);} // use default cube collision
	return sphere_city_obj_cylin_coll(pos, get_xy_coll_radius(), pos_, p_last, radius_, xlate, cnorm);
}

multi_model_city_obj_t::multi_model_city_obj_t(point const &pos_, float height, bool dim_, bool dir_, unsigned model_id, unsigned model_select, bool is_cylinder_) :
	model_city_obj_t(pos_, height, dim_, dir_, (model_id + (model_select << 8)), is_cylinder_), full_model_id(model_id + (model_select << 8)) {}

// swingsets

void swingset_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (anim_scale > 0.0) { // animated
		animation_state_t anim_state(1, ANIM_ID_SWINGS, anim_time); // enable_animations=1
		dstate.s.add_uniform_float("rotate_amount", anim_scale);
		model_city_obj_t::draw(dstate, qbds, dist_scale, shadow_only, &anim_state);
		dstate.s.add_uniform_float("rotate_amount", 1.0); // reset
		anim_state.clear_animation_id(dstate.s);
	}
	else { // static
		model_city_obj_t::draw(dstate, qbds, dist_scale, shadow_only);
	}
}
void swingset_t::next_frame(point const &camera_bs) {
	if (bcube.contains_pt_exp(camera_bs, 1.2*CAMERA_RADIUS)) {anim_scale = 1.0;}

	if (anim_scale > 0.01) {
		anim_time  += 0.0006*fticks;
		anim_scale *= (1.0f - min(1.0f, 0.0025f*fticks)); // dampen
	}
	else {anim_time = anim_scale = 0.0;}
}

// benches

bench_t::bench_t(point const &pos_, float radius_, bool dim_, bool dir_) : oriented_city_obj_t(pos_, radius_, dim_, dir_) {
	bcube.set_from_point(pos);
	bcube.expand_by(vector3d((dim ? 0.32 : 1.0), (dim ? 1.0 : 0.32), 0.0)*radius);
	bcube.z2() += 0.85*radius; // set bench height
	set_bsphere_from_bcube(); // calculate a more correct bsphere
}
cube_t bench_t::get_bird_bcube() const {
	cube_t top_place(bcube);
	top_place.expand_in_dim(!dim,  0.1*get_width ()); // expand the back outward a bit
	top_place.expand_in_dim( dim, -0.1*get_length()); // shrink a bit to account for the arms extending further to the sides than the back
	return top_place;
}
/*static*/ void bench_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {
		select_texture(FENCE_TEX);
		select_texture(get_normal_map_for_bldg_tid(FENCE_TEX), 5);
	}
}
/*static*/ void bench_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {bind_default_flat_normal_map();}
}
void bench_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	cube_t bcube_with_back(bcube);
	bcube_with_back.d[!dim][dir] += (dir ? 1.0 : -1.0)*0.2*radius; // adjust slightly since the back tilts outside the bcube
	if (!dstate.check_cube_visible(bcube_with_back, dist_scale)) return;

	cube_t cubes[] = { // Note: taken from mapx/bench.txt; front-back is in X
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
	float const dist_val(shadow_only ? 0.0 : p2p_dist(camera_pdu.pos, center)/dstate.draw_tile_dist);
	cube_t bc; // bench bounds

	for (unsigned i = 0; i < 12; ++i) { // back still contributes to bounds
		if (dir)  {swap(cubes[i].d[0][0], cubes[i].d[0][1]); cubes[i].d[0][0] *= -1.0; cubes[i].d[0][1] *= -1.0;}
		if (!dim) {swap(cubes[i].d[0][0], cubes[i].d[1][0]); swap(cubes[i].d[0][1], cubes[i].d[1][1]);}
		if (i == 0) {bc = cubes[i];} else {bc.union_with_cube(cubes[i]);}
	}
	point const c1(bcube.get_cube_center()), c2(bc.get_cube_center());
	vector3d const scale(bcube.dx()/bc.dx(), bcube.dy()/bc.dy(), bcube.dz()/bc.dz()); // scale to fit to target cube
	color_wrapper const cw(WHITE);
	unsigned const num(shadow_only ? 6U : max(1U, min(6U, unsigned(0.2/dist_val)))); // simple distance-based LOD, in pairs
	float const tscale(1.0/radius);
	bool const swap_tc_xy(dim);
	for (unsigned i = 1; i < 2*num; ++i) {dstate.draw_cube(qbds.qbd, ((cubes[i] - c2)*scale + c1), cw, 1, tscale, 0, 0, 0, swap_tc_xy);} // skip back
	// now draw the sloped back
	point pts[4] = {point(-1.0, -5.0, 5.0), point(-1.0, 5.0, 5.0), point(0.2, 5.0, 1.6), point(0.2, -5.0, 1.6)};
	point f[4], b[4];

	for (unsigned i = 0; i < 4; ++i) {
		if (dir)  {pts[i].x *= -1.0;}
		if (!dim) {swap(pts[i].x, pts[i].y);}
		pts[i] = ((pts[i] - c2)*scale + c1);
	}
	vector3d const normal(get_poly_norm(pts, 1)), delta((0.2f*scale.x)*normal); // thickness = 0.4
	UNROLL_4X(f[i_] = pts[i_] + delta;);
	qbds.qbd.add_quad_pts(f, WHITE,  normal);
	UNROLL_4X(b[i_] = pts[i_] - delta;);
	qbds.qbd.add_quad_pts(b, WHITE, -normal);

	for (unsigned i = 0; i < 4; ++i) { // draw sides
		unsigned const j((i+1)&3); // next i
		point const s[4] = {f[i], b[i], b[j], f[j]};
		qbds.qbd.add_quad_pts(s, WHITE, get_poly_norm(s, 1));
	}
}

// tree planters

tree_planter_t::tree_planter_t(point const &pos_, float radius_, float height) : city_obj_t(pos_, radius_) {
	set_bcube_from_vcylin(pos, height, radius);
}
/*static*/ void tree_planter_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture((dstate.pass_ix == 0) ? (int)DIRT_TEX : get_texture_by_name("roads/sidewalk.jpg"));}
}
void draw_xy_walls(cube_t const &bcube, cube_t const &hole, color_wrapper const &cw, float tscale, draw_state_t &dstate, quad_batch_draw &qbd) {
	cube_t walls[4] = {bcube, bcube, bcube, bcube}; // -X, +X, -Y, +Y
	walls[0].x2() = walls[2].x1() = walls[3].x1() = hole.x1();
	walls[1].x1() = walls[2].x2() = walls[3].x2() = hole.x2();
	walls[2].y2() = hole.y1();
	walls[3].y1() = hole.y2();

	for (unsigned d = 0; d < 2; ++d) {
		dstate.draw_cube(qbd, walls[d  ], cw, 1, tscale, 0); // X
		dstate.draw_cube(qbd, walls[d+2], cw, 1, tscale, 1); // Y, skip X dims
	}
}
void tree_planter_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	color_wrapper const cw(LT_GRAY);
	cube_t dirt(bcube);
	dirt.expand_by_xy(-0.1*dirt.get_size()); // shrink 10% on all XY sides

	if (dstate.pass_ix == 0) { // draw dirt
		dirt.z2() -= 0.25*bcube.dz(); // move down 25%
		dstate.draw_cube(qbds.qbd, dirt, cw, 1, 0.0, 3); // top only (skip X, Y, and bottom)
	}
	else { // draw stone
		draw_xy_walls(bcube, dirt, cw, 40.0, dstate, qbds.qbd);
	}
}

// trashcans

trashcan_t::trashcan_t(point const &pos_, float radius_, float height, bool is_cylin_) : city_obj_t(pos_, radius_), is_cylin(is_cylin_) {
	set_bcube_from_vcylin(pos, height, radius);
	set_bsphere_from_bcube(); // recompute bsphere from bcube
}
/*static*/ void trashcan_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (shadow_only) {} // nothing to do
	else if (dstate.pass_ix == 0) {select_texture(get_texture_by_name("roads/asphalt.jpg"));} // cube city/park
	else { // cylinder residential
		select_texture(get_texture_by_name("buildings/corrugated_metal.tif"));
		select_texture(get_texture_by_name("buildings/corrugated_metal_normal.tif", 1), 5);
		dstate.s.set_cur_color(GRAY);
	}
}
/*static*/ void trashcan_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only && dstate.pass_ix > 0) {bind_default_flat_normal_map();} // restore to default for cylindrical trashcan
	city_obj_t::post_draw(dstate, shadow_only);
}
void trashcan_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (is_cylin != (dstate.pass_ix == 1)) return; // wrong pass
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;

	if (is_cylin) { // cylindrical residential trashcan
		unsigned const ndiv(shadow_only ? 16 : max(4U, min(32U, unsigned(2.0f*dist_scale*dstate.draw_tile_dist/p2p_dist(dstate.camera_bs, pos)))));
		float const cylin_radius(get_cylin_radius()), lid_radius(1.08*cylin_radius), height(bcube.dz());
		point const rim_center(pos.x, pos.y, (bcube.z1() + 0.88*height)), lid_center(pos.x, pos.y, (bcube.z1() + 0.96*bcube.dz()));
		draw_fast_cylinder(point(pos.x, pos.y, bcube.z1()), rim_center, cylin_radius, cylin_radius, ndiv, 1, 0); // draw sides only

		if (!shadow_only && bcube.closest_dist_less_than(dstate.camera_bs, 0.4*dist_scale*dstate.draw_tile_dist)) { // draw the lid
			draw_fast_cylinder(rim_center, lid_center, lid_radius, lid_radius, ndiv, 1, 0, 0, nullptr, 0.1); // draw sides only with partial texture
			draw_fast_cylinder(lid_center, point(pos.x, pos.y, bcube.z2()), lid_radius, 0.001*lid_radius, ndiv, 1, 0); // draw sides only; not quite a cone

			if (bcube.closest_dist_less_than(dstate.camera_bs, 0.1*dist_scale*dstate.draw_tile_dist)) { // draw the handle
				color_wrapper const gray(GRAY);
				float const hlen(0.1*height), hwidth(0.02*height), thickness(0.005*height);
				cube_t top;
				set_wall_width(top, bcube.xc(), hwidth, 0);
				set_wall_width(top, bcube.yc(), hlen,   1);
				cube_t side(top); // copy xvals
				side.z1() = bcube.z2() - 0.03*height;
				side.z2() = top.z1() = bcube.z2() + 0.02*height;
				top .z2() = top.z1() + thickness;
				dstate.draw_cube(qbds.qbd, top, gray, 0, 0.001); // set tscale to map to a single texel

				for (unsigned s = 0; s < 2; ++s) {
					side.y1() = side.y2() = top.d[1][s];
					side.d[1][!s] += (s ? -1.0 : 1.0)*thickness;
					dstate.draw_cube(qbds.qbd, side, gray, 1, 0.001, 4); // set tscale to map to a single texel, skip top and bottom
				}
			}
		}
	}
	else { // cube city/park trashcan
		if (shadow_only) {dstate.draw_cube(qbds.qbd, bcube, WHITE, 1); return;} // draw a cube for the shadow
		color_wrapper const tan(colorRGBA(0.8, 0.6, 0.3, 1.0));
		cube_t hole(bcube);
		hole.expand_by_xy(-0.08*bcube.get_size()); // shrink on all XY sides
		draw_xy_walls(bcube, hole, tan, 25.0, dstate, qbds.qbd); // sides

		if (bcube.closest_dist_less_than(dstate.camera_bs, 0.4*dist_scale*dstate.draw_tile_dist)) {
			float const height(bcube.dz());
			cube_t bottom(hole);
			bottom.z2() -= 0.95*height;
			dstate.draw_cube(qbds.qbd, bottom, tan, 1, 30.0, 3); // inside bottom, top surface only
			cube_t top(hole);
			top.z1() += 0.92*height;
			top.z2() -= 0.02*height;
			cube_t top_hole(top);
			top_hole.expand_by_xy(-0.22*bcube.get_size()); // shrink on all XY sides
			draw_xy_walls(top, top_hole, color_wrapper(BROWN), 200.0, dstate, qbds.qbd); // brown top; technically don't need to draw some of the interior surfaces
		}
	}
}
bool trashcan_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (!is_cylin) {return city_obj_t::proc_sphere_coll(pos_, p_last, radius_, xlate, cnorm);} // use base class cube coll for cube trashcans
	return sphere_city_obj_cylin_coll(pos, get_cylin_radius(), pos_, p_last, radius_, xlate, cnorm);
}

// fire hydrants

fire_hydrant_t::fire_hydrant_t(point const &pos_, float radius_, float height, vector3d const &orient_) : city_obj_t(pos_, radius_), cylin_radius(radius), orient(orient_) {
	bcube.set_from_sphere(*this);
	set_cube_zvals(bcube, pos.z, pos.z+height);
	pos.z += 0.5*height; // pos is bottom center point, make it the center
	max_eq(radius, 0.5f*height); // use a more accurate bounding sphere; Note: no cube root of (r*r + r*r + h*h)
}
/*static*/ void fire_hydrant_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.set_cur_color(colorRGBA(1.0, 0.75, 0.0));} // override with custom color since the model color is black
	disable_hemi_lighting_pre_post(dstate, shadow_only, 0);
}
/*static*/ void fire_hydrant_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.set_cur_color(WHITE);} // restore to default color
	disable_hemi_lighting_pre_post(dstate, shadow_only, 1);
	city_obj_t::post_draw(dstate, shadow_only);
}
void fire_hydrant_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const { // Note: qbds are unused
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;

	if (!shadow_only) {
		if (dstate.is_occluded(bcube)) return;
		building_obj_model_loader.draw_model(dstate.s, pos, bcube, orient, WHITE, dstate.xlate, OBJ_MODEL_FHYDRANT, shadow_only);
	}
	else { // shadow pass: draw as a simple cylinder, untextured, top end only
		draw_fast_cylinder(point(pos.x, pos.y, bcube.z1()), point(pos.x, pos.y, bcube.z2()), 0.8*cylin_radius, 0.8*cylin_radius, 12, 0, 4); // ndiv=12
	}
}
bool fire_hydrant_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	return sphere_city_obj_cylin_coll(pos, cylin_radius, pos_, p_last, radius_, xlate, cnorm);
}

// plot dividers

plot_divider_type_t plot_divider_types[DIV_NUM_TYPES] = {
	plot_divider_type_t("cblock2.jpg", "normal_maps/cblock2_NRM.jpg", 0.50, 2.5, 1.0, 1, 0, WHITE, GRAY    ), // wall
	plot_divider_type_t("fence.jpg",   "normal_maps/fence_NRM.jpg",   0.15, 2.0, 1.0, 1, 0, WHITE, LT_BROWN), // fence
	plot_divider_type_t("hedges.jpg",  "", 1.00, 1.6, 1.0, 0, 0, GRAY, GREEN), // hedge - too short to be an occluder
	plot_divider_type_t("roads/chainlink_fence.png", "", 0.02, 1.55, 8.0, 0, 1, WHITE, GRAY) // chainlink fence with alpha mask; can't be taller than other fence types in case it intersects
};

/*static*/ void divider_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (dstate.pass_ix == DIV_NUM_TYPES) { // fence post pass - untextured
		if (!shadow_only) {
			select_texture(WHITE_TEX);
			dstate.s.set_specular(0.8, 60.0); // specular metal surface
		}
		return;
	}
	assert(dstate.pass_ix < DIV_NUM_TYPES);
	plot_divider_types[dstate.pass_ix].pre_draw(shadow_only);
}
/*static*/ void divider_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (dstate.pass_ix == DIV_NUM_TYPES) { // fence post pass
		if (!shadow_only) {dstate.s.clear_specular();}
	}
	else {plot_divider_types[dstate.pass_ix].post_draw(shadow_only);}
}
void divider_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (dstate.pass_ix == DIV_NUM_TYPES && type == DIV_CHAINLINK) { // add chainlink fence posts
		if (!dstate.check_cube_visible(bcube, 1.5*dist_scale)) return;
		float const length(get_width()), height(bcube.dz()), thickness(get_depth());
		float const post_hwidth(1.5*thickness), post_width(2.0*post_hwidth), top_width(1.5*thickness);
		unsigned const num_sections(ceil(0.3*length/height)), num_posts(num_sections + 1);
		float const post_spacing((length - post_width)/num_sections);
		color_wrapper cw(GRAY);
		cube_t post(bcube), top(bcube); // copy dim and Z values
		post.expand_in_dim(dim, 0.5f*(post_width - thickness)); // increase width to post_width
		top .expand_in_dim(dim, 0.5f*(top_width  - thickness));
		post.z2() += 0.02*height; // extend slightly above the top of the fence; bird feet may clip through the post when standing on it though
		set_wall_width(top, bcube.z2(), 0.5f*top_width, 2); // set height

		// for now we draw posts as cubes rather than cylinders since it's faster and easier, because we can use the existing qbd
		for (unsigned i = 0; i < num_posts; ++i) { // add posts
			set_wall_width(post, (bcube.d[!dim][0] + post_hwidth + i*post_spacing), post_hwidth, !dim);
			dstate.draw_cube(qbds.qbd, post, cw, 1);
		}
		dstate.draw_cube(qbds.qbd, top, cw, 1);
		return;
	}
	if (type != dstate.pass_ix) return; // this type not enabled in this pass
	if (type == DIV_CHAINLINK) {dist_scale *= 0.5;} // less visible
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	assert(dstate.pass_ix < DIV_NUM_TYPES);
	plot_divider_type_t const &pdt(plot_divider_types[dstate.pass_ix]);
	dstate.draw_cube(qbds.qbd, bcube, color_wrapper(pdt.color), 1, pdt.tscale/bcube.dz(), skip_dims); // skip bottom, scale texture to match the height

	if (!shadow_only && type == DIV_HEDGE && bcube.closest_dist_less_than(dstate.camera_bs, 0.35f*(X_SCENE_SIZE + Y_SCENE_SIZE))) {
		dstate.hedge_draw.add(bcube); // draw detailed leaves for nearby hedges
	}
}
bool divider_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	cube_t bcube_wide(bcube + xlate);
	bcube_wide.expand_in_dim(dim, max(0.0f, 0.5f*(0.5f*building_t::get_scaled_player_radius() - get_depth()))); // make sure it's at least half player radius in thickness
	return sphere_cube_int_update_pos(pos_, radius_, bcube_wide, p_last, 0, cnorm);
}

void hedge_draw_t::create(cube_t const &bc) {
	bcube = bc - bc.get_cube_center(); // centered on the origin
	unsigned const target_num_leaves(40000);
	vector3d const sz(bcube.get_size());
	float const leaf_sz(0.05*sz.z), surf_area(sz.x*sz.y + 2.0f*sz.z*(sz.x + sz.y));
	float const side_areas[5] = {sz.y*sz.z, sz.y*sz.z, sz.x*sz.z, sz.x*sz.z, sz.x*sz.y};
	rand_gen_t rgen;
	vector<vert_norm_comp_tc> verts;
	verts.reserve(4*target_num_leaves);

	for (unsigned n = 0; n < 5; ++n) { // {+X, -X, +Y, -Y, +Z} sides
		unsigned const dim(n>>1), dir(n&1), d1((dim+1)%3), d2((dim+2)%3);
		unsigned const num_this_face(target_num_leaves*side_areas[n]/surf_area);
		point pos;
		pos[dim] = bcube.d[dim][!dir];

		for (unsigned n = 0; n < num_this_face; ++n) {
			pos[d1] = rgen.rand_uniform(bcube.d[d1][0], bcube.d[d1][1]);
			pos[d2] = rgen.rand_uniform(bcube.d[d2][0], bcube.d[d2][1]);
			vector3d const normal(rgen.signed_rand_vector_spherical().get_norm());
			float const angle(TWO_PI*rgen.rand_float());
			vector3d tangent;
			rotate_vector3d(cross_product(normal, plus_x), normal, angle, tangent);
			vector3d const binormal(cross_product(normal, tangent)), dx(leaf_sz*tangent), dy(leaf_sz*binormal);
			point const pts[4] = {(pos - dx - dy), (pos + dx - dy), (pos + dx + dy), (pos - dx + dy)};
			float const ts [4] = {0.0, 1.0, 1.0, 0.0}, tt[4] = {0.0, 0.0, 1.0, 1.0};
			for (unsigned i = 0; i < 4; ++i) {verts.emplace_back(pts[i], normal, ts[i], tt[i]);}
		} // for n
	} // for s
	num_verts = verts.size();
	create_and_upload(verts, 0, 1);
}
void hedge_draw_t::draw_and_clear(shader_t &s) {
	if (empty()) return;
	if (!vbo_valid()) {create(to_draw.front());}
	select_texture(get_texture_by_name("pine2.jpg"));
	enable_blend(); // slightly smoother, but a bit of background shows through
	s.add_uniform_float("min_alpha", 0.9);
	s.add_uniform_int("two_sided_lighting", 1);
	pre_render();
	vector3d const sz_mult(bcube.get_size().inverse());
	// we can almost use an instance_render here, but that requires changes to the shader or a custom shader

	for (auto c = to_draw.begin(); c != to_draw.end(); ++c) {
		bool const swap_dims((c->dx() < c->dy()) ^ (bcube.dx() < bcube.dy())); // wrong dim, rotate 90 degrees
		vector3d sz(c->get_size());
		if (swap_dims) {swap(sz.x, sz.y);}
		fgPushMatrix();
		translate_to(c->get_cube_center()); // align the center
		if (swap_dims) {fgRotate(90.0, 0.0, 0.0, 1.0);}
		scale_by(sz_mult*sz); // scale to match the size
		s.upload_mvm();
		draw_quads_as_tris(num_verts);
		++num_frame_draw_calls;
		fgPopMatrix();
	} // for c
	post_render();
	s.add_uniform_int("two_sided_lighting", 0); // reset
	s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA); // restore to the default
	disable_blend();
	to_draw.clear();
}

// swimming pools

void begin_water_surface_draw() {
	if (0) {
		select_texture(WHITE_TEX);
		select_texture(get_texture_by_name("normal_maps/ocean_water_normal.png", 1), 5);
	}
	else {select_texture(get_texture_by_name("snow2.jpg"));}
	enable_blend(); // transparent water
	glDepthMask(GL_FALSE);
}
void end_water_surface_draw() {
	if (0) {bind_default_flat_normal_map();}
	disable_blend();
	glDepthMask(GL_TRUE);
}

// passes: 0=in-ground walls, 1=in-ground water, 2=above ground sides, 3=above ground water
/*static*/ void swimming_pool_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (shadow_only) {} // nothing
	else if (dstate.pass_ix == 2) {select_texture(WHITE_TEX);} // sides/untextured
	else if (dstate.pass_ix == 0) {select_texture(get_texture_by_name("bathroom_tile.jpg"));} // walls and maybe ladder
	else if (dstate.pass_ix == 1 || dstate.pass_ix == 3) {begin_water_surface_draw();} // water surface
	else if (dstate.pass_ix == 4) {} // caustics pass - handled by the caller
	else {assert(0);}
}
/*static*/ void swimming_pool_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only && dstate.pass_ix != 4) {dstate.s.set_cur_color(WHITE);} // restore to default color
	if (dstate.pass_ix == 1 || dstate.pass_ix == 3) {end_water_surface_draw();} // water surface
	city_obj_t::post_draw(dstate, shadow_only);
}
void swimming_pool_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	bool const caustics_pass(dstate.pass_ix == 4);
	if (((dstate.pass_ix > 1) ^ above_ground) && !caustics_pass) return; // not drawn in this pass
	if (caustics_pass && wcolor.A > 0.75) return; // water too murky, disable caustics
	float cull_dist_scale(dist_scale);
	if (caustics_pass) {min_eq(cull_dist_scale, 0.75f*get_inf_terrain_fog_dist()/dstate.draw_tile_dist);} // cull caustics at 75% fog scale since they don't have fog enabled
	if (!dstate.check_cube_visible(bcube, cull_dist_scale))      return;
	float const water_zval(get_water_zval());

	if (above_ground) { // cylindrical; bcube should be square in XY
		point const camera_bs(dstate.camera_bs);
		float const radius(get_radius()), xc(bcube.xc()), yc(bcube.yc()), dscale(dist_scale*dstate.draw_tile_dist), height(bcube.dz()), inner_bottom(bcube.z1() + 0.04f*height);
		unsigned const ndiv(shadow_only ? 24 : max(4U, min(64U, unsigned(6.0f*dscale/p2p_dist(camera_bs, pos)))));
		point const orig_cpos(camera_pos);
		camera_pos = dstate.camera_bs; // required for proper two sided cylinder normals

		if (dstate.pass_ix == 2) { // draw sides, bottom, and maybe ladder
			dstate.s.set_cur_color(color);
			draw_fast_cylinder(point(xc, yc, bcube.z1()), point(xc, yc, bcube.z2()), radius, radius, ndiv, 0, 0, 1); // untextured, no ends; two sided lighting

			if (!shadow_only) { // draw bottom, shifted slightly up
				dstate.s.set_cur_color(color*0.4); // darker due to light atten
				draw_circle_normal(0.0, radius, ndiv, 0, point(xc, yc, inner_bottom));
			}
			if (bcube.closest_dist_less_than(camera_bs, 0.5*dscale)) { // draw ladder
				unsigned const num_steps = 5;
				color_wrapper const step_color(LT_GRAY);
				cube_t ladder;
				float const side_pos(bcube.d[dim][dir]), swidth((dir ? 1.0 : -1.0)*0.065*radius); // ladder is on this side of the pool
				float const ladder_height(1.2*height), step_delta(ladder_height/(num_steps + 0.25)), step_offset(0.25*step_delta), step_height(0.14*step_delta);
				ladder.d[dim][!dir] = side_pos;
				ladder.d[dim][ dir] = side_pos + swidth;
				set_wall_width(ladder, (dim ? xc : yc), 0.16*radius, !dim);
				bool const is_close(bcube.closest_dist_less_than(camera_bs, 0.2*dscale));
				bool const is_very_close(is_close && bcube.closest_dist_less_than(camera_bs, 0.1*dscale));

				for (unsigned n = 0; n < num_steps; ++n) { // draw steps
					ladder.z1() = bcube .z1() + n*step_delta + step_offset;
					ladder.z2() = ladder.z1() + step_height;
					dstate.draw_cube(qbds.qbd, ladder, step_color, !is_very_close); // skip bottom if not close
				}
				if (is_close) { // draw bars
					float const bars_top(bcube.z1() + ladder_height), bar_radius(0.012*radius);
					bool const draw_top(is_very_close && camera_bs.z > bars_top);
					unsigned const bars_ndiv(max(3U, min(12U, ndiv/4)));
					point pt;
					pt.z    = bcube.z1();
					pt[dim] = side_pos + 0.4*swidth; // slightly off center
					dstate.s.set_cur_color(colorRGBA(0.33, 0.33, 0.33));

					for (unsigned n = 0; n < 2; ++n) {
						pt[!dim] = ladder.d[!dim][n] + (n ? -1.0 : 1.0)*2.0*bar_radius;
						draw_fast_cylinder(pt, point(pt.x, pt.y, bars_top), bar_radius, bar_radius, bars_ndiv, 0, (draw_top ? 4 : 0)); // untextured with top
					}
				}
			}
		}
		else if (!shadow_only && dstate.pass_ix == 3) { // draw water surface; not for the shadow pass
			dstate.s.set_cur_color(wcolor);
			draw_circle_normal(0.0, radius, ndiv, 0, point(xc, yc, water_zval)); // shift slightly below the top
		}
		else if (caustics_pass) { // draw sides and bottom caustics
			assert(!shadow_only);
			float const tscale = 2.0;
			int const bias_loc(dstate.s.get_uniform_loc("shad_bias_scale"));
			assert(bias_loc >= 0);
			dstate.s.set_uniform_float(bias_loc, 0.0); // disable - not needed, and looks slighlty better without this
			draw_circle_normal(0.0, 0.999*radius, ndiv, 0, point(xc, yc, inner_bottom), tscale, tscale); // draw bottom, shifted slightly up, and slightly smaller to avoid clipping
			dstate.s.set_uniform_float(bias_loc, CITY_BIAS_SCALE); // restore the default
			glEnable(GL_CULL_FACE); // inner surface only
			glCullFace(GL_FRONT);
			draw_fast_cylinder(point(xc, yc, inner_bottom), point(xc, yc, water_zval), radius, radius, ndiv, 1, 0, 1, nullptr, 0.3, 0.0, nullptr, 0, 6.0); // textured, sides
			glCullFace(GL_BACK);
			glDisable(GL_CULL_FACE);
		}
		camera_pos = orig_cpos;
	}
	else { // in-ground
		float const height(bcube.dz()), wall_thick(0.14*height), tscale((caustics_pass ? 0.1 : 0.5)/wall_thick);
		cube_t inner(bcube);
		inner.expand_by_xy(-wall_thick);

		if (dstate.pass_ix == 0 || caustics_pass) { // draw walls or caustics
			// draw sides
			color_wrapper const cw(caustics_pass ? WHITE : color);
			cube_t sides[4] = {bcube, bcube, bcube, bcube}; // {S, N, W center, E center}
			sides[0].y2() = sides[2].y1() = sides[3].y1() = inner.y1();
			sides[1].y1() = sides[2].y2() = sides[3].y2() = inner.y2();
			sides[2].x2() = inner.x1();
			sides[3].x1() = inner.x2();
			
			for (unsigned d = 0; d < 4; ++d) {
				if (caustics_pass) {sides[d].z2() = water_zval;} // clip to water level
				dstate.draw_cube(qbds.qbd, sides[d], cw, 1, tscale, ((d > 2) ? 2 : 0)); // skip_bottom=1
			}
			// draw bottom
			bool const dmax(bcube.dx() < bcube.dy());
			float const bz1(bcube.z1() + 0.5*wall_thick);
			point pts[4] = {point(bcube.x1(), bcube.y1(), bz1), point(bcube.x2(), bcube.y1(), bz1), point(bcube.x2(), bcube.y2(), bz1), point(bcube.x1(), bcube.y2(), bz1)};
			
			if (sloped) { // draw sloped bottom surface
				for (unsigned n = 0; n < 4; ++n) {
					if (pts[n][dmax] == bcube.d[dmax][dir]) {pts[n].z = (bz1 + 0.5*height);}
				}
			}
			qbds.qbd.add_quad_pts(pts, cw, get_poly_norm(pts, 1), tex_range_t(0.0, 0.0, tscale*bcube.dx(), tscale*bcube.dy()));
		}
		else if (dstate.pass_ix == 1) { // draw water surface
			inner.z2() = water_zval; // shift water surface a bit below the top
			dstate.draw_cube(qbds.qbd, inner, color_wrapper(wcolor), 1, 0.5*tscale, 3); // draw top water as semi-transparent
		}
	}
}
bool swimming_pool_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (above_ground) {
		if (!sphere_cube_intersect((pos_ - xlate), radius_, bcube)) return 0; // optimization
		float const radius(get_radius()), xc(bcube.xc() + xlate.x), yc(bcube.yc() + xlate.y), z1(bcube.z1() + xlate.z), z2(bcube.z2() + xlate.z);
		return sphere_vert_cylin_intersect(pos_, radius_, cylinder_3dw(point(xc, yc, z1), point(xc, yc, z2), radius, radius), cnorm); // checks sides
	}
	cube_t bcube_tall(bcube + xlate);
	bcube_tall.z2() += CAMERA_RADIUS + camera_zh; // extend upward so that player collision detection works better
	return sphere_cube_int_update_pos(pos_, radius_, bcube_tall, p_last, 0, cnorm);
}
bool swimming_pool_t::update_depth_if_underwater(point const &p, float &depth) const {
	if (!bcube.contains_pt(p)) return 0;
	float const water_zval(get_water_zval());
	if (p.z >= water_zval)     return 0;
	if (above_ground && !dist_xy_less_than(p, pos, get_radius())) return 0; // check circle containment
	depth = (water_zval - pos.z);
	return 1;
}

// pool decks

textured_mat_t pool_deck_mats[NUM_POOL_DECK_PASSES] = {
	textured_mat_t("fence.jpg",          "normal_maps/fence_NRM.jpg",    0, WHITE, LT_BROWN),
	textured_mat_t("roads/concrete.jpg", "",                             0, GRAY,  LT_GRAY ),
	textured_mat_t("shingles.jpg",       "normal_maps/shingles_NRM.jpg", 0, WHITE, GRAY    ), // roof
	textured_mat_t("wood.jpg",           "normal_maps/wood_NRM.jpg",     0, WHITE, BROWN   ), // pillars
};

pool_deck_t::pool_deck_t(cube_t const &base_, cube_t const &roof_, unsigned mat_id_, bool dim_, bool dir_) : // dim/dir is facing yard
	oriented_city_obj_t(base_, dim_, dir_), mat_id(mat_id_ % NUM_POOL_DECK_TYPES), base(base_), roof(roof_)
{
	if (roof.is_all_zeros()) return; // no roof
	bcube.union_with_cube(roof);
	set_bsphere_from_bcube(); // update
}
void pool_deck_t::calc_pillars(cube_t const &ladder) {
	if (!has_roof())      return; // no roof
	if (!pillars.empty()) return; // already calculated
	float const roof_width(roof.get_sz_dim(!dim)), roof_height(roof.z1() - base.z2());
	unsigned const num_pillars(2U + unsigned(0.5*roof_width/roof_height));
	float const pillar_width(1.8*roof.dz()), pillar_hwidth(0.5*pillar_width), pillar_edge_gap(0.1*pillar_width);
	float const span(roof_width - pillar_width - 2.0*pillar_edge_gap), spacing(span/(num_pillars - 1));
	cube_t pillar(roof), avoid(ladder);
	set_cube_zvals(pillar, base.z2(), roof.z1());
	pillar.d[dim][!dir] = roof.d[dim][dir] - (dir ? 1.0 : -1.0)*pillar_width;
	pillar.translate_dim(dim, -(dir ? 1.0 : -1.0)*pillar_edge_gap); // shift slightly toward the house
	if (!ladder.is_all_zeros()) {avoid.expand_by_xy(1.0*pillar_width);}
	pillars.reserve(num_pillars);

	for (unsigned n = 0; n < num_pillars; ++n) {
		set_wall_width(pillar, (roof.d[!dim][0] + pillar_hwidth + pillar_edge_gap + n*spacing), pillar_hwidth, !dim);
		if (!avoid.is_all_zeros() && pillar.intersects(avoid)) continue; // skip if close to the ladder
		pillars.push_back(pillar);
	}
}
/*static*/ void pool_deck_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	assert(dstate.pass_ix < NUM_POOL_DECK_PASSES);
	pool_deck_mats[dstate.pass_ix].pre_draw(shadow_only);
}
/*static*/ void pool_deck_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	assert(dstate.pass_ix < NUM_POOL_DECK_PASSES);
	pool_deck_mats[dstate.pass_ix].post_draw(shadow_only);
}
void pool_deck_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (dstate.pass_ix < NUM_POOL_DECK_TYPES) { // draw the deck
		if (mat_id != dstate.pass_ix) return; // this type not enabled in this pass
		dstate.draw_cube(qbds.qbd, base, pool_deck_mats[mat_id].color, 1, 1.2/get_depth(), 0, 0, 0, !dim); // skip bottom; what about the side against the house?
	}
	else if (dstate.pass_ix == NUM_POOL_DECK_TYPES) { // draw the roof; should the roof slope downward in dim/dir?
		if (!has_roof()) return; // no roof
		dstate.draw_cube(qbds.qbd, roof, pool_deck_mats[mat_id].color, 1, 1.5/get_depth(), 3, 0, dir, dim); // draw top surface only
	}
	else if (dstate.pass_ix == NUM_POOL_DECK_TYPES+1) { // draw the support pillars and inside/edges of the roof
		if (!has_roof()) return; // no roof
		if (!shadow_only && !roof.closest_dist_less_than(dstate.camera_bs, 0.3*dist_scale*dstate.draw_tile_dist)) return; // too far
		// Note: roof side against the house is drawn in case it's visible, but it's between floors and not visible from inside the house
		dstate.draw_cube(qbds.qbd, roof, pool_deck_mats[mat_id].color, 0, 4.0/get_depth(), 0, 0, dir, dim, 1.0, 1.0, 1.0, 1); // roof sides and bottom; skip_top=1
		float const tscale(4.0/(roof.z1() - base.z2()));
		for (cube_t const &pillar : pillars) {dstate.draw_cube(qbds.qbd, pillar, pool_deck_mats[mat_id].color, 1, tscale, 4, 0, 0, 1);} // skip top/bot; swap tc x/y
	}
	else {assert(0);}
}
bool pool_deck_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (pillars.empty() || !sphere_cube_intersect(pos_, radius_, (bcube + xlate))) return 0;
	bool coll(0);
	for (cube_t const &pillar : pillars) {coll |= sphere_cube_int_update_pos(pos_, radius_, (pillar + xlate), p_last, 0, cnorm);}
	return coll;
}

// newsracks

newsrack_t::newsrack_t(point const &pos_, float height, float width, float depth, bool dim_, bool dir_, unsigned style_, colorRGBA const &color_) :
	oriented_city_obj_t(pos_, 0.5*sqrt(height*height + width*width + depth*depth), dim_, dir_), color(color_), style(style_)
{
	bcube.set_from_point(pos_);
	bcube.expand_in_dim( dim, 0.5*depth);
	bcube.expand_in_dim(!dim, 0.5*width);
	bcube.z2() += height;
}
cube_t newsrack_t::get_bird_bcube() const {
	if ((style & 3) != 1) return bcube; // not cube with coin mechanism on top - use full bcube (only z2 face is needed)
	// return coin mechanism bcube
	vector3d const sz(bcube.get_size());
	cube_t cm(bcube);
	cm.z1() = bcube.z1() + 0.7*sz.z;
	cm.expand_in_dim(!dim, -0.3*sz[!dim]); // shrink width
	cm.d[dim][!dir] += (dir ? 1.0 : -1.0)*0.6*sz[dim]; // shift back inward
	return cm;
}
/*static*/ void newsrack_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture(get_texture_by_name("roads/fake_news.jpg"));}
	if (!shadow_only) {dstate.s.set_specular(0.33, 40.0);} // low specular
}
/*static*/ void newsrack_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.clear_specular();}
}
void newsrack_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!bcube.closest_dist_less_than(dstate.camera_bs, 0.45*dist_scale*dstate.draw_tile_dist)) { // far away, draw low detail single cube
		dstate.draw_cube(qbds.qbd, bcube, color, 1); // skip_bottom=1
		return;
	}
	// use a tiny texture scale so that we can use the newspaper texture for drawing the sides with the white texel in the LLC
	float const dir_sign(dir ? 1.0 : -1.0), llc_tscale(0.0001);
	bool const front_facing((dstate.camera_bs[dim] < bcube.get_center_dim(dim)) ^ dir);
	vector3d const sz(bcube.get_size());
	cube_t body(bcube);
	bool skip_bottom(1);

	switch (style & 3) {
	case 0: break; // simple cube
	case 1: { // cube with coin mechanism on top
		cube_t const cm(get_bird_bcube());
		body.z2() = cm.z1();
		dstate.draw_cube(qbds.qbd, cm, color, 1, llc_tscale); // coin mech; skip bottom

		if (front_facing) { // draw the lock bar
			cube_t bar(cm);
			bar.expand_in_dim(!dim, -0.1*sz[!dim] ); // shrink width
			bar.d[dim][!dir]  = cm.d[dim][dir]; // flush with front of coin mech
			bar.d[dim][ dir] += dir_sign*0.05*sz[dim]; // extend outward a bit
			bar.z2() -= 0.5*cm.dz();
			bar.z1() -= 0.3*cm.dz();
			dstate.draw_cube(qbds.qbd, bar, color, 0, llc_tscale); // skip face dim, !dir, or draw if player is in front?
		}
		break;
	}
	case 2: { // cube with extended front
		cube_t stand(bcube);
		stand.z2() = body.z1() = bcube.z1() + 0.35*sz.z;
		stand.d[dim][dir] -= dir_sign*0.1*sz[dim]; // shift front inward
		dstate.draw_cube(qbds.qbd, stand, color, 1, llc_tscale, 4); // stand, skip top and bottom
		skip_bottom = 0; // front of bottom may be visible
		break;
	}
	case 3: { // cube on narrow stand
		cube_t stand(bcube), base(bcube);
		stand.z2() = body .z1() = bcube.z1() + 0.50*sz.z;
		base .z2() = stand.z1() = bcube.z1() + 0.06*sz.z;
		stand.expand_by_xy(-0.3*min(sz.x, sz.y)); // shrink in XY
		dstate.draw_cube(qbds.qbd, stand, color, 1, llc_tscale, 4); // stand, skip top and bottom
		dstate.draw_cube(qbds.qbd, base,  color, 1, llc_tscale); // base, skip bottom
		skip_bottom = 0; // edges of bottom may be visible
		break;
	}
	} // end switch
	dstate.draw_cube(qbds.qbd, body, color, skip_bottom, llc_tscale); // main body
	
	if (front_facing) { // draw the door
		cube_t door(body);
		door.expand_in_dim( 2,   -0.1*body.dz()); // shrink height
		door.expand_in_dim(!dim, -0.1*sz[!dim] ); // shrink width
		max_eq(door.z1(), (door.z2() - 1.1f*door.get_sz_dim(!dim))); // shift up if it's too tall compared to width
		door.d[dim][!dir]  = body.d[dim][dir]; // flush with front of body
		door.d[dim][ dir] += dir_sign*0.02*sz[dim]; // extend outward a bit
		dstate.draw_cube(qbds.qbd, door, WHITE);
	}
}

// power poles

power_pole_t::power_pole_t(point const &base_, point const &center_, float pole_radius_, float height, float wires_offset_,
	float const pole_spacing_[2], uint8_t dims_, bool at_grid_edge_, bool const at_line_end_[2], bool residential_) :
	at_grid_edge(at_grid_edge_), residential(residential_), dims(dims_), pole_radius(pole_radius_), wires_offset(wires_offset_), base(base_), center(center_)
{
	UNROLL_2X(pole_spacing[i_] = pole_spacing_[i_]; at_line_end[i_] = at_line_end_[i_];)
	bcube.set_from_point(center);
	bcube.z2() += height;
	pos    = bcube.get_cube_center();
	radius = bcube.get_bsphere_radius();
	for (unsigned d = 0; d < 2; ++d) {bcube.expand_in_dim(d, (has_dim_set(!d) ? get_bar_extend() : pole_radius));} // add bar if dim bit not set
	bcube_with_wires = bcube; // cache for visibility query; could also recompute on each call

	for (unsigned d = 0; d < 2; ++d) {
		if (at_line_end[d] || !has_dim_set(d)) continue;
		bcube_with_wires.d[d][0] -= pole_spacing[d]; // extend to include wires
		bcube_with_wires.translate_dim(d, wires_offset); // Note: wires_offset should be 0 for poles that have both wire dims enabled
	}
	bsphere_radius = bcube_with_wires.furthest_dist_to_pt(pos);
}
cube_t power_pole_t::get_ped_occluder() const {
	cube_t occluder(base);
	occluder.z2() = bcube.z2();
	occluder.expand_by_xy(pole_radius/SQRT2); // take inner radius to reduce the occluder size for the pedestrian
	return occluder;
}
cube_t power_pole_t::calc_cbar(bool d) const {
	cube_t cbar;
	cbar.z1() = bcube.z1() + (d ? 0.90 : 0.96)*bcube.dz(); // stagger the heights in X vs. Y so that wires don't intersect on poles that have them in both dims
	cbar.z2() = cbar .z1() + 0.7*pole_radius;
	set_wall_width(cbar, (center[ d] + 1.3*pole_radius), 0.3*pole_radius, d); // offset
	set_wall_width(cbar,  center[!d], get_bar_extend(), !d);
	return cbar;
}
void power_pole_t::get_wires_conn_pts(point pts[3], bool d) const {
	float const wire_spacing(get_hwire_spacing()), wire_radius(get_wire_radius()), standoff_height(get_standoff_height());
	float const offsets[3] = {-wire_spacing, -0.3f*wire_spacing, wire_spacing}; // offset from the center to avoid intersecting the pole
	cube_t const cbar(calc_cbar(d));

	for (unsigned n = 0; n < 3; ++n) {
		pts[n][d]  = cbar.get_center_dim(d);
		pts[n][!d] = center[!d] + offsets[n]; // set wire offset
		pts[n].z   = cbar.z2()  + standoff_height + wire_radius; // resting on top of the standoff
	}
}
point power_pole_t::get_nearest_connection_point(point const &to_pos, bool near_power_pole) const {
	float dmin_sq(0.0);
	point ret(to_pos); // start at to_pos; will return this if no wires can be connected to

	for (unsigned d = 0; d < 2; ++d) {
		if (!has_dim_set(d) || at_line_end[d]) continue; // no wires in this dim
		point pw(center);
		pw[!d] = base[!d] + ((center[!d] != base[!d]) ? -1.0 : 1.0)*0.5*get_power_pole_offset(); // on the side of the pole
		pw.z  += 0.75*bcube.dz() - get_vwire_spacing(); // connect to middle wire; must match value used in drawing
		cube_t wires_bcube(pw, pw);
		wires_bcube.d[d][0] -= pole_spacing[d]; // extend to the adjacent pole
		wires_bcube.translate_dim(d, wires_offset); // offset wires in case the pole was moved to avoid a driveway
		point conn_pos(wires_bcube.closest_pt(to_pos)); // will be directly across with conn_pos[d] = to_pos[d]
		
		if (near_power_pole) { // off to the side of the pole
			// offset by an additional value based on distance so that wires coming from different points intersect at different places
			float const run_delta(to_pos[d] - base[d]);
			conn_pos[d] = base[d] + SIGN(run_delta)*pole_radius*(2.0 + 4.0*fabs(run_delta)/pole_spacing[d]);
		}
		float const dsq(p2p_dist_sq(to_pos, conn_pos));
		if (dmin_sq == 0.0 || dsq < dmin_sq) {ret = conn_pos; dmin_sq = dsq;}
	} // for d
	return ret;
}
void power_pole_t::get_top_wire_end_pts(point top_wires[2][3][2]) const { // {dim, wire_ix, end_ix}
	float const wire_radius(get_wire_radius());

	for (unsigned d = 0; d < 2; ++d) {
		if (!has_dim_set(d) || at_line_end[d]) continue; // no wires in this dim
		point conn_pts[3];
		get_wires_conn_pts(conn_pts, d);

		for (unsigned n = 0; n < 3; ++n) {
			conn_pts[n].z += wire_radius; // top of wire
			top_wires[d][n][0] = top_wires[d][n][1] = conn_pts[n];
			top_wires[d][n][1][d] -= pole_spacing[d]; // extend to the adjacent pole
		}
	} // for d
}
point power_pole_t::get_transformer_center() const { // not checking has_transformer()
	float const pole_height(bcube.dz()), tf_radius(2.0*pole_radius), tf_height(0.1*pole_height), y_sign(at_line_end[1] ? 1.0 : -1.0);
	point tf_pos(base);
	tf_pos.z += 0.77*pole_height + 0.5*tf_height;
	tf_pos.y += y_sign*(tf_radius + pole_radius);
	return tf_pos;
}
bool power_pole_t::add_wire(point const &p1, point const &p2, bool add_pole) { // Note: p1 connects to building or streetlight; p2 connects to wires on pole
	wire_t wire(p1, p2);
	
	if (add_pole) { // used for houses
		wire.pts[0]   .z += 0.040*bcube.dz(); // set the wire pole height
		wire.pole_base.z -= 0.006*bcube.dz(); // extend below the roof
		if (check_city_building_line_coll_bs_any(wire.pts[0], wire.pts[1])) return 0; // placement failed
	}
	wires.push_back(wire);
	for (unsigned d = 0; d < 2; ++d) {bcube_with_wires.union_with_sphere(wire.pts[d], get_wire_radius());} // okay to omit pole_base
	bsphere_radius = bcube_with_wires.furthest_dist_to_pt(pos); // recompute
	return 1;
}
/*static*/ void power_pole_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (shadow_only) return;
	select_texture(WOOD2_TEX);
	select_texture(get_texture_by_name("normal_maps/wood_NRM.jpg", 1), 5);
}
/*static*/ void power_pole_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {bind_default_flat_normal_map();} // restore to default
	city_obj_t::post_draw(dstate, shadow_only);
}

void add_cylin_as_tris(vector<vert_norm_tc_color> &verts, point const ce[2], float r1, float r2, color_wrapper const &cw,
	unsigned ndiv, unsigned draw_top_bot, float tst=1.0, float tss=1.0, bool swap_ts_tt=0)
{
	// added as individual triangles; would be more efficient to use indexed triangles
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, r1, r2, ndiv, v12));
	vert_norm_tc_color quad_pts[4];

	for (unsigned i = 0; i < ndiv; ++i) { // similar to gen_cylinder_quads(), but with a color
		for (unsigned j = 0; j < 2; ++j) {
			unsigned const S(i + j), s(S%ndiv);
			vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]); // normalize?
			float const ts(S*tss);
			quad_pts[2*j+0].assign(vpn.p[(s<<1)+!j], normal, (swap_ts_tt ?      j *tst : ts), (swap_ts_tt ? ts :      j *tst), cw.c);
			quad_pts[2*j+1].assign(vpn.p[(s<<1)+ j], normal, (swap_ts_tt ? (1.0-j)*tst : ts), (swap_ts_tt ? ts : (1.0-j)*tst), cw.c);
		}
		for (unsigned n = 0; n < 6; ++n) {verts.push_back(quad_pts[q2t_ixs[n]]);}

		for (unsigned d = 0; d < 2; ++d) { // draw bottom and top triangle(s)
			if (!(draw_top_bot & (1<<d))) continue;
			unsigned const I((i+1)%ndiv);
			vector3d const normal(d ? v12 : -v12);
			verts.emplace_back(ce[d], normal, 0.5, 0.5, cw);
			verts.emplace_back(vpn.p[(i<<1)+d], normal, 0.5*(1.0 + vpn.n[i].x), 0.5*(1.0 + vpn.n[i].y), cw);
			verts.emplace_back(vpn.p[(I<<1)+d], normal, 0.5*(1.0 + vpn.n[I].x), 0.5*(1.0 + vpn.n[I].y), cw);
		}
	} // for i
}
void draw_wire(point const *const pts, float radius, color_wrapper const &cw, quad_batch_draw &untex_qbd) { // pts is size 2
	unsigned const ndiv(4);
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(pts, radius, radius, ndiv, v12));

	for (unsigned i = 0; i < ndiv; ++i) { // similar to gen_cylinder_quads()
		unsigned const in((i+1)%ndiv);
		unsigned const pt_ixs[4] = {(i<<1)+1, (i<<1), (in<<1), (in<<1)+1};
		for (unsigned n = 0; n < 6; ++n) {untex_qbd.verts.emplace_back(vpn.p[pt_ixs[q2t_ixs[n]]], plus_z, 0, 0, cw.c);}
	}
}
void draw_ortho_wire(point const &p, float radius, float pole_spacing, bool d, color_wrapper const &cw, draw_state_t &dstate, quad_batch_draw &untex_qbd) {
	cube_t wire(p, p);
	wire.d[d][0] -= pole_spacing; // extend to the adjacent pole
	wire.expand_in_dim(!d, radius);
	wire.expand_in_dim(2,  radius);
	// black, don't need normals/tcs/colors; could use indexed triangles, but the time taken to draw these wires is insignificant (< 1% of total frame time)
	dstate.draw_cube(untex_qbd, wire, cw, 0); // since wires are black, and we can't see the ends, we can't even tell they're cubes rather than cylinders
}
void draw_standoff_geom(point const ce[2], float radius, float dmax, point const &camera_bs, color_wrapper const &cw, quad_batch_draw &untex_qbd) {
	unsigned const ndiv(max(4U, min(32U, unsigned(0.33f*dmax/p2p_dist(camera_bs, ce[1])))));

	if (ndiv <= 16) { // single truncated cone
		bool const draw_top(dot_product((ce[1] - ce[0]), (camera_bs - ce[1])) > 0.0);
		add_cylin_as_tris(untex_qbd.verts, ce, radius, 0.75*radius, cw, ndiv, (draw_top ? 2 : 0));
	}
	else { // multiple truncated cones
		unsigned const num_segs = 4;
		vector3d const step_delta((ce[1] - ce[0])/num_segs);
		point const ce_part[2] = {ce[0], ce[0]+step_delta};
		unsigned const verts_start(untex_qbd.verts.size());
		add_cylin_as_tris(untex_qbd.verts, ce_part, radius, 0.75*radius, cw, 16, 3); // truncated cone with top and bottom
		unsigned const verts_end(untex_qbd.verts.size());
	
		for (unsigned n = 1; n < num_segs; ++n) {
			vector3d const delta(n*step_delta);
		
			for (unsigned v = verts_start; v < verts_end; ++v) {
				untex_qbd.verts.push_back(untex_qbd.verts[v]);
				untex_qbd.verts.back().v += delta;
			}
		}
	}
}
void draw_vert_standoff(point const &p1, point const &camera_bs, float height, float radius, float delta_offset, float dmax, bool dim,
	bool is_first, unsigned verts_start, unsigned &verts_end, color_wrapper const &cw, quad_batch_draw &untex_qbd)
{
	if (is_first) { // first standoff, draw a truncated cone
		point ce[2] = {p1, p1};
		ce[1].z += height;
		draw_standoff_geom(ce, radius, dmax, camera_bs, cw, untex_qbd);
		verts_end = untex_qbd.verts.size();
	}
	else { // next standoff, copy and translate the previous truncated cone
		for (unsigned v = verts_start; v < verts_end; ++v) {
			untex_qbd.verts.push_back(untex_qbd.verts[v]);
			untex_qbd.verts.back().v[!dim] += delta_offset;
		}
	}
}

// Note: power line connectivity is all handled here in draw();
// the current grid is fully connected, forming loops on both the upper three high voltage lines and lower three low voltage lines;
// maybe this isn't realistic, but it does have a nice symmetry and higher apparent wiring complexity; the user will likely not notice
void power_pole_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	point const camera_bs(dstate.camera_bs);
	float const dmax(shadow_only ? camera_pdu.far_ : dist_scale*dstate.draw_tile_dist);
	if (!bcube.closest_dist_less_than(camera_bs, dmax)) return;
	if (!camera_pdu.cube_visible((shadow_only ? bcube : bcube_with_wires) + dstate.xlate)) return;
	color_wrapper const black(BLACK), white(colorRGBA(0.7, 0.7, 0.7)), gray(colorRGBA(0.4, 0.4, 0.4)), cw(LT_BROWN); // darken the wood color
	bool const pole_visible(camera_pdu.cube_visible(bcube + dstate.xlate));
	float const wire_radius(get_wire_radius()), pole_height(bcube.dz());
	cube_t tf_bcube;
	point conduit_top(all_zeros);
	quad_batch_draw &m_qbd(qbds.untex_qbd), &s_qbd(qbds.untex_spec_qbd); // {matte, specular}

	if (pole_visible) {
		unsigned const ndiv(shadow_only ? 16 : max(4U, min(32U, unsigned(1.5f*dmax/p2p_dist(camera_bs, pos)))));
		unsigned const pole_ndiv(min(ndiv, 24U));
		float const vert_tscale = 10.0;
		point ce[2] = {base, get_top()};
		bool const draw_top(ce[1].z < camera_bs.z);
		add_cylin_as_tris(qbds.qbd.verts, ce, pole_radius, pole_radius, cw, pole_ndiv, (draw_top ? 2 : 0), vert_tscale, 1.0/pole_ndiv, 1); // swap_ts_tt=1

		if (has_transformer() && (shadow_only || bcube.closest_dist_less_than(camera_bs, 0.7*dmax))) { // draw transformer, untextured
			float const tf_radius(2.0*pole_radius), tf_height(0.1*pole_height), y_sign(at_line_end[1] ? 1.0 : -1.0);
			ce[0].z = base.z  + 0.77*pole_height;
			ce[1].z = ce[0].z + tf_height;
			ce[0].x = ce[1].x = base.x;
			ce[0].y = ce[1].y = base.y + y_sign*(tf_radius + pole_radius); // offset in -y, +y at end (so that wires cross above it)
			bool const draw_top_bot(camera_bs.z > 0.5f*(ce[0].z + ce[1].z));
			add_cylin_as_tris(s_qbd.verts, ce, tf_radius, tf_radius, gray, ndiv, (draw_top_bot ? 2 : 1)); // specular
			tf_bcube.set_from_points(ce, 2);
			tf_bcube.expand_by_xy(tf_radius);

			if (!residential && ndiv > 4) { // draw conduit for wires that go into the ground
				float const cradius(3.5*wire_radius);
				conduit_top.assign(base.x, (base.y + y_sign*(0.5f*cradius + pole_radius)), (base.z + 0.6*pole_height)); // below the transformer
				point const cce[2] = {point(conduit_top.x, conduit_top.y, base.z), conduit_top};
				bool const draw_top(ndiv > 8 && camera_bs.z > conduit_top.z);
				add_cylin_as_tris(s_qbd.verts, cce, cradius, cradius, gray, min(ndiv, 16U), (draw_top ? 2 : 0)); // specular
			}
		}
	}
	float const wire_spacing(get_hwire_spacing()), vwire_spacing(get_vwire_spacing());
	float const standoff_height(get_standoff_height()), standoff_radius(0.25*pole_radius);
	point wire_pts[3][2];
	unsigned wire_mask(0);
	bool drew_wires(0);

	for (unsigned d = 0; d < 2; ++d) { // {x, y}
		if (!has_dim_set(d)) continue; // no wires in this dim
		float const offsets[3] = {-wire_spacing, -0.3f*wire_spacing, wire_spacing}; // offset from the center to avoid intersecting the pole
		cube_t const cbar(calc_cbar(d));
		point p1;
		p1[d] = cbar.get_center_dim(d);
		p1.z  = cbar.z2(); // resting on top of the bar
		
		if (pole_visible && (shadow_only || cbar.closest_dist_less_than(camera_bs, 0.5*dmax))) { // could test cbar cube visible, but unclear if that's faster
			dstate.draw_cube(qbds.qbd, cbar, cw, 0, 0.8/cbar.dz()); // draw all sides

			if (!shadow_only && cbar.closest_dist_less_than(camera_bs, 0.15*dmax)) { // draw insulator standoffs
				unsigned verts_start(s_qbd.verts.size()), verts_end(0);

				for (unsigned n = 0; n < 3; ++n) {
					p1[!d] = center[!d] + offsets[n]; // set wire offset
					wire_pts[n][d] = p1 + vector3d(0.0, 0.0, (standoff_height + wire_radius));

					if (n == 1 && d == 1) { // offset middle wire from standoff to avoid clipping through pole
						if      (!at_line_end[1]) {wire_pts[n][1].y = wire_pts[n][0].y;}
						else if (!at_line_end[0]) {wire_pts[n][0].x = wire_pts[n][1].x;}
						else {} // I guess it clips through the pole in this case
					}
					float const delta_offset(offsets[n] - offsets[0]);
					draw_vert_standoff(p1, camera_bs, standoff_height, standoff_radius, delta_offset, dmax, d, (n == 0), verts_start, verts_end, white, s_qbd); // specular
				} // for n
				wire_mask |= (1 << d); // mark wires as drawn in this dim
			}
			if (d == 1 && !tf_bcube.is_all_zeros()) { // connect wire to transformer if running in y dim
				float const spacing(0.3*tf_bcube.get_sz_dim(!d));
				point const tf_top_center(cube_top_center(tf_bcube));

				for (unsigned n = 0; n < 3; ++n) {
					point const pts[2] = {point(center.x+offsets[n], tf_top_center.y, p1.z+standoff_height+wire_radius),
						(tf_top_center + vector3d((n - 1.0)*spacing, 0.0, standoff_height-0.5f*wire_radius))}; // top wire, bottom transformer
					draw_wire(pts, wire_radius, black, m_qbd);
				}
				if (!shadow_only && tf_bcube.closest_dist_less_than(camera_bs, 0.1*dmax)) { // draw insulator standoffs
					unsigned verts_start(s_qbd.verts.size()), verts_end(0);

					for (unsigned n = 0; n < 3; ++n) {
						point const p2((tf_top_center.x + (n - 1.0)*spacing), tf_top_center.y, tf_top_center.z);
						draw_vert_standoff(p2, camera_bs, standoff_height, standoff_radius, n*spacing, dmax, d, (n == 0), verts_start, verts_end, white, s_qbd); // specular
					}
					wire_mask |= (1 << d); // mark wires as drawn in this dim
				}
			}
		}
		if (shadow_only) continue; // skip wires for shadow pass since they don't show up reliably
		bool const is_offset(center[!d] != base[!d]);
		float const sep_dist(0.5*get_power_pole_offset()), offset_sign(is_offset ? -1.0 : 1.0);
		float const bot_wire_zval(base.z + 0.75*pole_height), bot_wire_pos(base[!d] + offset_sign*sep_dist), thick_wire_delta_z(0.042*pole_height);

		if (!at_line_end[d]) { // no wires at end pole
			cube_t wires_bcube(cbar);
			min_eq(wires_bcube.z1(), (bot_wire_zval - 3.0f*vwire_spacing - thick_wire_delta_z)); // include lower wires
			UNROLL_2X(wires_bcube.d[d][i_] = bcube_with_wires.d[d][i_];)
			if (!wires_bcube.closest_dist_less_than(camera_bs, 0.45*dmax) || !camera_pdu.cube_visible(wires_bcube + dstate.xlate)) continue; // wires distance/VFC
			// draw the three top wires and one bottom wire in this dim
			p1[d] += wires_offset; // offset wires in case the pole was moved to avoid a driveway
			p1.z  += standoff_height + wire_radius; // resting on top of the standoff

			for (unsigned n = 0; n < 3; ++n) { // top wires
				p1[!d] = center[!d] + offsets[n]; // set wire spacing
				draw_ortho_wire(p1, wire_radius, pole_spacing[d], d, black, dstate, m_qbd);
			}
			// bottom 3 wires: split the difference between the offset corner and street power poles and attach to different sides of the poles
			float const wire_extend(sep_dist - pole_radius - 0.5*cbar.get_sz_dim(d)); // extend to fill gap between outer wire at standoffs and wire ending at cbar
			point pw;
			pw[!d] = bot_wire_pos; // on the side of the pole
			pw[d]  = p1[d] + wire_extend; // extend slightly to meet the crossing wire
			pw.z   = bot_wire_zval;
			float const bot_wire_extend(pole_spacing[d] + sep_dist + wire_radius); // overlap with next wire if there is one, extend past standoff at last pole
			float const cable_wire_radius(2.0*wire_radius);

			for (unsigned n = 0; n < 4; ++n) { // 3 bottom wires + thicker cable TV wire bundle as well
				if (n == 3) {pw.z -= thick_wire_delta_z;}
				draw_ortho_wire(pw, ((n == 3) ? cable_wire_radius : wire_radius), bot_wire_extend, d, black, dstate, m_qbd);
				if (n < 3) {pw.z -= vwire_spacing;}
			}
			// draw cable TV repeater/junction box
			float const box_hlen(1.8*vwire_spacing), box_radius(0.55*vwire_spacing);
			pw.z  -= (box_radius + cable_wire_radius); // just touching the bottom of the wire
			pw[d] -= (0.2 + 0.6*fract(12345*bcube.x1() + 54321*bcube.y1()))*pole_spacing[d]; // random spacing
			point epts[2] = {pw, pw};
			epts[0][d] -= box_hlen; epts[1][d] += box_hlen;
			cube_t rbcube(epts[0], epts[1]);
			rbcube.expand_by(box_radius);
			
			if (camera_pdu.cube_visible(rbcube + dstate.xlate)) {
				unsigned const ndiv(max(4U, min(24U, unsigned(0.75f*dmax/p2p_dist(camera_bs, pw)))));
				add_cylin_as_tris(s_qbd.verts, epts, box_radius, box_radius, color_wrapper(BKGRAY), ndiv, 3); // draw both ends; specular
			}
			drew_wires = 1;
		}
		if (pole_visible && bcube.closest_dist_less_than(camera_bs, 0.15*dmax)) {
			// bottom standoffs, even if there's no wire (because it will extend on end pole)
			point pb;
			pb[ d] = base[d]; // reset to the center of the pole
			pb[!d] = bot_wire_pos;
			pb.z   = bot_wire_zval;
			point ce[2] = {pb, pb};
			ce[1][!d] -= offset_sign*wire_radius; // end attached to the wire
			ce[0][!d]  = base[!d] + 0.96*offset_sign*pole_radius; // end attached to the pole, slightly offset into the pole
			unsigned verts_start(s_qbd.verts.size()), verts_end(0);

			for (unsigned n = 0; n < 4; ++n) { // 3 bottom wires + thicker cable TV wire bundle as well
				if (n == 0) { // first standoff, draw a truncated cone
					draw_standoff_geom(ce, standoff_radius, dmax, camera_bs, white, s_qbd);
					verts_end = s_qbd.verts.size();
				}
				else { // next standoff, copy and translate the previous truncated cone
					for (unsigned v = verts_start; v < verts_end; ++v) {
						s_qbd.verts.push_back(s_qbd.verts[v]);
						s_qbd.verts.back().v.z -= n*vwire_spacing;
						if (n == 3) {s_qbd.verts.back().v.z -= thick_wire_delta_z;}
					}
				}
			} // for n
			if (d == 1 && !tf_bcube.is_all_zeros()) { // vertical wire up to transformer
				point const tf_conn_pt(pb.x, (tf_bcube.yc() - 0.5f*tf_bcube.dy()), (tf_bcube.z1() + 0.9*tf_bcube.dz())); // along the side near the top
				cube_t wire(tf_conn_pt, tf_conn_pt);
				wire.z1()  = pb.z + wire_radius - 2.0*vwire_spacing; // meet the top of the lowest wire
				wire.z2() += wire_radius; // span entire standoff
				wire.expand_by_xy(wire_radius);
				dstate.draw_cube(m_qbd, wire, black, 1); // skip bottom
				// standoff for wire
				point ce[2] = {tf_conn_pt, tf_conn_pt};
				vector3d const so_dir(-1.0, 1.0, 0.0);
				ce[0] += 0.4*pole_radius*so_dir;
				ce[1] += 0.5*wire_radius*so_dir;
				draw_standoff_geom(ce, standoff_radius, dmax, camera_bs, white, s_qbd);

				if (conduit_top != all_zeros) { // draw wires to the top of the ground conduit
					point const conn_pt(tf_conn_pt.x, tf_conn_pt.y, (pb.z - 2.0*vwire_spacing)); // connect to bottom wire
					point const pts[2] = {(conduit_top - vector3d(0.0, 0.0, wire_radius)), conn_pt};
					draw_wire(pts, wire_radius, black, m_qbd);
				}
			}
		}
	} // for d
	if (drew_wires && wire_mask == 3 && bcube.closest_dist_less_than(camera_bs, 0.25*dmax)) { // both dims set, connect X and Y wires
		for (unsigned n = 0; n < 3; ++n) {draw_wire(wire_pts[n], wire_radius, black, m_qbd);}
	}
	if (!wires.empty() && bcube_with_wires.closest_dist_less_than(camera_bs, 0.3*dmax)) {
		for (auto &w : wires) { // represents all three wires tied together
			draw_wire(w.pts, wire_radius, black, m_qbd);
			// draw vertical wire segment connecting the three, which also represents all three power wires tied together
			cube_t wire(w.pts[1], w.pts[1]); // connection point to bottom horizontal wires
			wire.expand_in_dim(2, vwire_spacing); // connect to wires above and below
			wire.expand_by_xy(wire_radius);
			dstate.draw_cube(m_qbd, wire, black, 1, 0.0, 4); // skip top and bottom
			if (w.pole_base.z == w.pts[0].z || !dist_less_than(w.pts[0], camera_bs, 0.15*dmax)) continue; // no pole, or too far away
			point const ce[2] = {w.pole_base, (w.pts[0] + vector3d(0.0, 0.0, wire_radius))};
			float const radius(1.5f*wire_radius);
			unsigned const ndiv(max(4U, min(16U, unsigned(0.1f*dmax/p2p_dist(camera_bs, ce[1])))));
			bool const draw_top(camera_bs.z > 0.5f*(ce[0].z + ce[1].z));
			add_cylin_as_tris(m_qbd.verts, ce, radius, radius, gray, ndiv, (draw_top ? 2 : 1)); // this one is a cylinder
		} // for w
	}
}
bool power_pole_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (!sphere_cube_intersect((pos_ - xlate), radius_, bcube)) return 0; // optimization
	// for now we only check the pole itself rather than the cross bar
	return sphere_vert_cylin_intersect(pos_, radius_, cylinder_3dw((base + xlate), (get_top() + xlate), pole_radius, pole_radius), cnorm); // check pole using the base radius
}

// transmission lines

float get_tline_right_of_way() {return 0.3*city_params.road_width;} // include tower bar length (.25) plus some extra spacing

void transmission_line_t::calc_bcube() { // bcube of the right of way
	bcube = cube_t(p1, p2); // endpoints

	for (unsigned n = 0; n < 3; ++n) { // probably only need to include center point [1], but okay to include them all
		bcube.union_with_pt(p1_wire_pts[n]);
		bcube.union_with_pt(p2_wire_pts[n]);
	}
	bcube.z1() -= tower_height; // approximate
	bcube.expand_by_xy(get_tline_right_of_way());
}
bool transmission_line_t::sphere_intersect_xy(point const &pos, float radius) const {
	if (!sphere_cube_intersect_xy(pos, radius, bcube)) return 0;
	float const right_of_way(radius + get_tline_right_of_way());
	if (point_line_seg_dist_2d(pos, p1, p2            ) < right_of_way) return 1; // check the line of towers
	if (point_line_seg_dist_2d(pos, p1, p1_wire_pts[1]) < right_of_way) return 1; // check end points using center wire
	if (point_line_seg_dist_2d(pos, p2, p2_wire_pts[1]) < right_of_way) return 1; // check end points using center wire
	return 0;
}
bool transmission_line_t::cube_intersect_xy(cube_t const &c) const {
	if (bcube.is_all_zeros()) { // endpoints not yet connected; check p1-p2 only (for building generation)
		cube_t c_exp(c);
		c_exp.expand_by_xy(get_tline_right_of_way());
		return check_line_clip_xy(p1, p2, c_exp.d);
	}
	if (!bcube.intersects_xy(c)) return 0;
	cube_t c_exp(c);
	c_exp.expand_by_xy(get_tline_right_of_way());
	if (check_line_clip_xy(p1, p2, c_exp.d)) return 1; // check the line of towers
	if (check_line_clip_xy(p1, p1_wire_pts[1], c_exp.d)) return 1; // check end points using center wire
	if (check_line_clip_xy(p2, p2_wire_pts[1], c_exp.d)) return 1; // check end points using center wire
	return 0;
}

// handicap spaces

hcap_space_t::hcap_space_t(point const &pos_, float radius_, bool dim_, bool dir_) : oriented_city_obj_t(pos_, radius_, dim_, dir_) {
	assert(radius > 0.0);
	set_bcube_from_vcylin(pos, 0.01*radius, radius); // make nonzero height
}
/*static*/ void hcap_space_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	assert(!shadow_only); // not drawn in the shadow pass
	select_texture(get_texture_by_name("roads/handicap_parking.jpg"));
}
void hcap_space_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	float const x1(bcube.x1()), y1(bcube.y1()), x2(bcube.x2()), y2(bcube.y2()), z(bcube.z2());
	point const pts[4] = {point(x1, y1, z), point(x2, y1, z), point(x2, y2, z), point(x1, y2, z)};
	qbds.qbd.add_quad_pts(pts, WHITE, plus_z, tex_range_t(0.0, float(dir^dim), 1.0, float(dir^dim^1), 0, !dim));
}

hcap_with_dist_t::hcap_with_dist_t(hcap_space_t const &hs, cube_t const &plot, vect_cube_t &bcubes, unsigned bcubes_end) : hcap_space_t(hs), dmin_sq(plot.dx() + plot.dy()) {
	assert(bcubes_end <= bcubes.size());
	point const center(hs.bcube.get_cube_center());
		
	for (auto c = bcubes.begin(); c != bcubes.begin()+bcubes_end; ++c) {
		min_eq(dmin_sq, p2p_dist_xy_sq(center, c->get_cube_center())); // use distance to building center; c->closest_pt() has too many ties
	}
}

// manholes

manhole_t::manhole_t(point const &pos_, float radius_) : city_obj_t(pos_, radius_) {
	set_bcube_from_vcylin(pos, get_height(), radius);
}
/*static*/ void manhole_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	assert(!shadow_only); // not drawn in the shadow pass
	select_texture(MANHOLE_TEX);
	dstate.s.set_cur_color(colorRGBA(0.5, 0.35, 0.25, 1.0)); // gray-brown
}
void manhole_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	unsigned const ndiv(max(4U, min(32U, unsigned(1.0f*dist_scale*dstate.draw_tile_dist/p2p_dist(dstate.camera_bs, pos)))));
	draw_circle_normal(0.0, radius, ndiv, 0, point(pos.x, pos.y, pos.z+get_height()), -1.0); // draw top surface, invert texture coords
}

// trampolines

trampoline_t::trampoline_t(point const &pos_, float height, rand_gen_t &rgen) :
	model_city_obj_t(pos_, height, rgen.rand_bool(), rgen.rand_bool(), get_model_id(), 1) // random dim and dir; is_cylinder=1
{
	unsigned const NUM_COLORS = 5;
	colorRGBA const colors[NUM_COLORS] = {BLUE, RED, YELLOW, BLUE, GREEN};
	color = colors[rgen.rand() % NUM_COLORS];
}

// flowers: can use get_texture_by_name("sunflower.jpg")

// traffic cones

traffic_cone_t::traffic_cone_t(point const &pos_, float radius_) : city_obj_t(pos_, radius_) {
	set_bcube_from_vcylin(pos, 2.0*radius, 0.67*radius);
	pos.z += radius; // move to the center
}
/*static*/ void traffic_cone_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	dstate.s.add_uniform_int("two_sided_lighting", 1);
}
/*static*/ void traffic_cone_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	dstate.s.add_uniform_int("two_sided_lighting", 0);
}
void traffic_cone_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	cube_t base(bcube);
	base.z2() = bcube.z1() + 0.05*bcube.dz();
	color_wrapper const cw(colorRGBA(1.0, 0.25, 0.0, 1.0)); // dark orange
	unsigned const ndiv(max(4U, min(32U, unsigned(2.0f*dist_scale*dstate.draw_tile_dist/p2p_dist(dstate.camera_bs, pos)))));
	dstate.draw_cube(qbds.untex_qbd, base, color_wrapper(BKGRAY), 1); // skip bottom
	float const xc(bcube.xc()), yc(bcube.yc());
	float const cone_height(bcube.z2() - base.z2()), cone_r1(0.54*radius), cone_r2(0.13*radius), cone_dr(cone_r1 - cone_r2);
	point ce[2] = {point(xc, yc, base.z2()), point(xc, yc, bcube.z2())};
	add_cylin_as_tris(qbds.untex_qbd.verts, ce, cone_r1, cone_r2, cw, ndiv, 0); // draw sides, no end
	ce[0].z += 0.5*cone_height;
	ce[1].z -= 0.2*cone_height;
	add_cylin_as_tris(qbds.untex_qbd.verts, ce, 1.02*(cone_r1 - 0.5*cone_dr), 1.02*(cone_r2 + 0.2*cone_dr), color_wrapper(WHITE), ndiv, 0); // stripe
}

// ponds

pond_t::pond_t(point const &pos_, float x_radius, float y_radius, float depth) : city_obj_t(pos_, max(x_radius, y_radius)) {
	bcube.set_from_point(pos);
	bcube.expand_in_dim(0, x_radius);
	bcube.expand_in_dim(1, y_radius);
	bcube.z1() -= depth;
}
/*static*/ void pond_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	assert(!shadow_only);
	if      (dstate.pass_ix == 0) {select_texture(DIRT_TEX);} // dirt below
	else if (dstate.pass_ix == 1) {select_texture(BLUR_CENT_TEX); enable_blend();} // dark blur
	else {begin_water_surface_draw();} // water above
}
/*static*/ void pond_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if      (dstate.pass_ix == 0) {} // dirt below
	else if (dstate.pass_ix == 1) {disable_blend();} // dark blur
	else {end_water_surface_draw();} // water above
}
void pond_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	assert(!shadow_only);
	float const dist(p2p_dist(dstate.camera_bs, pos)), dz_off(max(0.0001f*bcube.dz(), 0.00025f*dist));
	unsigned const ndiv(max(4U, min(64U, unsigned(6.0f*dist_scale*dstate.draw_tile_dist/dist))));
	fgPushMatrix();
	translate_to(point(pos.x, pos.y, bcube.z2()));
	fgScale(bcube.dx(), bcube.dy(), 1.0); // set correct aspect ratio

	if (dstate.pass_ix == 0) { // dirt below
		float const tscale(4.0), bot_radius(0.5*0.5); // half the total radius
		point const bot_center(point(0.0, 0.0, -bcube.dz()));
		// bottom
		dstate.s.set_cur_color(GRAY_BLACK); // darker
		draw_circle_normal(0.0, bot_radius, ndiv, 0, bot_center, tscale, tscale); // invert_normals=0
		// sloped sides
		point const ce[2] = {bot_center, all_zeros};
		vector3d v12;
		vector_point_norm const &vpn(gen_cylinder_data(ce, bot_radius, 0.5, ndiv, v12));
		static vector<vert_norm_tc_color> verts;
		verts.resize(2U*(ndiv+1U));
		color_wrapper const cw_outer(GRAY), cw_inner(GRAY_BLACK);

		for (unsigned S = 0; S <= ndiv; ++S) {
			unsigned const s(S%ndiv), vix(2*S);
			vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]); // points down, must negate
			point const &p1(vpn.p[(s<<1)+0]), &p2(vpn.p[(s<<1)+1]);
			verts[vix+0].assign(p1, -normal, tscale*p1.x, tscale*p1.y, cw_inner.c);
			verts[vix+1].assign(p2, -normal, tscale*p2.x, tscale*p2.y, cw_outer.c);
		}
		draw_and_clear_verts(verts, GL_TRIANGLE_STRIP);
	}
	else if (dstate.pass_ix == 1) {
		float const tscale(1.0);
		dstate.s.set_cur_color(BLACK);
		draw_circle_normal(0.0, 0.5, ndiv, 0, point(0.0, 0.0, dz_off), tscale, tscale);
	}
	else { // water above
		float const tscale(2.0);
		dstate.s.set_cur_color(colorRGBA(0.2, 0.3, 0.5, 0.5)); // semi-transparent
		draw_circle_normal(0.0, 0.5, ndiv, 0, point(0.0, 0.0, 2.0*dz_off), tscale, tscale);
	}
	fgPopMatrix();
}
bool pond_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const { // pos_ is in camera space
	point const pos_bs(pos_ - xlate);
	if (!bcube.contains_pt_xy_exp(pos_bs, radius_)) return 0;
	vector3d const delta(pos_bs - pos);
	float const xv(delta.x/(radius_ + 0.5*bcube.dx())), yv(delta.y/(radius_ + 0.5*bcube.dy()));
	if (xv*xv + yv*yv > 1.0) return 0; // dist^2 > 1.0, outside the ellipse
	pos_ = p_last; // finding the actual intersection point requires solving a quartic equation, so simply revert to the last pos
	if (cnorm) {*cnorm = vector3d(delta.x, delta.y, 0.0).get_norm();} // assume collision normal is in the XY plane
	return 1;
}
bool point_in_ellipse(point const &p, cube_t const &c) {
	if (!c.contains_pt_xy(p)) return 0;
	float const xv((p.x - c.xc())/(0.5*c.dx())), yv((p.y - c.yc())/(0.5*c.dy()));
	return (xv*xv + yv*yv < 1.0);
}
bool pond_t::point_contains_xy(point const &p) const { // p is in global space
	return point_in_ellipse(p, bcube);
}
bool pond_t::update_depth_if_underwater(point const &p, float &depth) const {
	if (p.z <= bcube.z1() || p.z > bcube.z2() || !point_contains_xy(p)) return 0;
	depth = (bcube.z2() - p.z);
	return 1;
}

// walkways

walkway_t::walkway_t(bldg_walkway_t const &w) : oriented_city_obj_t(w, w.dim, 0), walkway_base_t(w) { // dir=0 (unused)
	assert(floor_spacing > 0.0);
	// use the roof because this is what's visible in overhead map mode
	map_mode_color = texture_color(global_building_params.get_material(roof_mat_ix).roof_tex.tid).modulate_with(roof_color);
	for (unsigned d = 0; d < 2; ++d) {open_ends[d] = w.open_ends[d];}
}
/*static*/ void walkway_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {bind_default_flat_normal_map();}
}
void walkway_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	tid_nm_pair_dstate_t state(dstate.s); // pass this in?

	if (dstate.camera_bs[!dim] < bcube.d[!dim][0] || dstate.camera_bs[!dim] > bcube.d[!dim][1]) { // camera not inside walkway projection - draw sides
		auto const &side_mat(global_building_params.get_material(side_mat_ix));
		side_mat.side_tex.set_gl(state);
		float const tsx(side_mat.side_tex.get_drawn_tscale_x()), tsy(side_mat.side_tex.get_drawn_tscale_y());
		dstate.draw_cube(qbds.qbd, bcube, side_color, 0, 1.0, ((1 << unsigned(dim)) | 4), 0, 0, 0, tsx, tsx, tsy); // sides; skip ends, top, and bottom
	}
	if (dstate.camera_bs.z < bcube.z1() || dstate.camera_bs.z > bcube.z2()) { // camera not inside walkway zvals - draw top and bottom
		auto const &roof_mat(global_building_params.get_material(roof_mat_ix));
		qbds.qbd.draw_and_clear(); // must draw here since texture was set dynamically
		roof_mat.roof_tex.set_gl(state);
		dstate.draw_cube(qbds.qbd, bcube, roof_color, 0, roof_mat.roof_tex.get_drawn_tscale_x(), 3); // top and bottom; skip ends and sides
	}
	qbds.qbd.draw_and_clear(); // must draw here since texture was set dynamically
}
bool walkway_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	// Note: if the player is coming from a building, this code won't be called until pos_ is outside the building bcube,
	// which can cause a one frame glitch during the transition; when coming from the walkway, the player isn't in the building until they're outside the walkway
	cube_t const bc(bcube + xlate);
	if (!bc.contains_pt_xy_exp(pos_, radius_)) return 0; // include radius so that we can walk from a building roof onto a walkway
	float const pzmax(max(pos_.z, p_last.z));
	if (pzmax < bc.z1()) return 0; // below the walkway

	if (pzmax > bc.z2()) { // above the walkway
		float const wwz(bc.z2() + radius_);
		if (pzmax > wwz) return 0; // in the air above the walkway
		pos_.z = wwz;
		if (cnorm) {*cnorm = plus_z;}
		return 1; // collision with roof
	}
	// else inside the walkway
	float const fc_thickness(0.5*FLOOR_THICK_VAL_WINDOWLESS*floor_spacing);
	float zval(bc.z2() - floor_spacing + fc_thickness); // bottom of upper walkway floor, assuming a windowless (city) office building
	assert(zval >= bc.z1());

	for (; zval >= bc.z1(); zval -= floor_spacing) {
		if (zval > (pzmax + radius_)) continue; // wrong floor - above top of sphere
		player_in_walkway = 1; // assumes only the player can be here
		float const wwz(zval + radius_);
		bool const floor_coll(pos_.z <= wwz);
		
		if (floor_coll) { // collision with floor
			pos_.z = wwz;
			if (cnorm) {*cnorm = plus_z;}
		}
		point const prev(pos_);
		max_eq(pos_[!dim], bc.d[!dim][0]+radius_); // clamp in !dim
		min_eq(pos_[!dim], bc.d[!dim][1]-radius_);
		if (pos_ == prev) return floor_coll; // not touching a side
		if (cnorm) {*cnorm = (pos_ - prev).get_norm();}
		return 1; // collision with side
	} // for zval
	return 0; // shouldn't get here
}
cube_t walkway_t::get_floor_occluder() const {
	// conservative: assumes player is on the bottom floor rather than an upper floor, since the actual floor is difficult to determine
	cube_t occluder(bcube);
	occluder.z2() = bcube.z1() + 0.3*floor_spacing; // matches bottom window zval calculation in building_t::get_all_drawn_interior_verts()
	return occluder;
}

// pillars

/*static*/ void pillar_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (shadow_only) {} // nothing to do
	else if (dstate.pass_ix == 0) {select_texture(get_texture_by_name("roads/concrete.jpg"));} // concrete cube
	else { // untextured steel cylinder
		select_texture(WHITE_TEX);
		dstate.s.set_specular(0.8, 60.0); // specular metal surface
		dstate.s.set_cur_color(WHITE);
	}
}
/*static*/ void pillar_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (shadow_only) {} // nothing to do
	else if (dstate.pass_ix == 0) {} // concrete cube
	else {dstate.s.clear_specular();} // untextured steel cylinder
}
void pillar_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (is_concrete != (dstate.pass_ix == 0)) return; // wrong pass
	if (is_concrete) {dstate.draw_cube(qbds.qbd, bcube, WHITE, 1, 16.0, 4);} // skip top and bottom
	else {
		unsigned const ndiv(shadow_only ? 8 : max(4U, min(32U, unsigned(2.0f*dist_scale*dstate.draw_tile_dist/p2p_dist(dstate.camera_bs, pos)))));
		float const cylin_radius(get_cylin_radius());
		draw_fast_cylinder(point(pos.x, pos.y, bcube.z1()), point(pos.x, pos.y, bcube.z2()), cylin_radius, cylin_radius, ndiv, 1, 0); // draw sides only
	}
}
bool pillar_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (is_concrete) {return city_obj_t::proc_sphere_coll(pos_, p_last, radius_, xlate, cnorm);} // use base class cube coll for concrete pillars
	return sphere_city_obj_cylin_coll(point(pos.x, pos.y, bcube.z1()), get_cylin_radius(), pos_, p_last, radius_, xlate, cnorm);
}

 // ww_elevator_t

/*static*/ void ww_elevator_t::pre_draw (draw_state_t &dstate, bool shadow_only) {enable_blend ();}
/*static*/ void ww_elevator_t::post_draw(draw_state_t &dstate, bool shadow_only) {disable_blend();}

void ww_elevator_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	dstate.draw_cube(qbds.untex_qbd, bcube, colorRGBA(0.8, 1.0, 0.9, 0.25), 1); // skip bottom
}
bool ww_elevator_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	// TODO: special logic to allow the player to enter and ride the elevator; this may go elsewhere
	return sphere_cube_int_update_pos(pos_, radius_, (bcube + xlate), p_last, 0, cnorm);
}

// parking lot solar roofs

parking_solar_t::parking_solar_t(cube_t const &c, bool dim_, bool dir_, unsigned ns, unsigned nr) : oriented_city_obj_t(c, dim_, dir_), num_spaces(ns), num_rows(nr) {
	set_bsphere_from_bcube();
}
/*static*/ void parking_solar_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture(get_solarp_tid());}
}
void parking_solar_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	// draw roof and solar panel
	float const height(bcube.dz());
	cube_t roof(bcube);
	roof.z1() = bcube.z2() - 0.05*height;
	cube_t panel(roof);
	panel.z2() += 0.004*height; // slightly above the roof
	panel.expand_by_xy(-0.04*height); // small shrink to create a border
	unsigned const vs1(qbds.qbd.verts.size()), vs2(qbds.untex_qbd.verts.size());
	dstate.draw_cube(qbds.qbd,       panel, WHITE,   1, 1.0/height, 3, 0, 0, 0, 1.0, 1.0, 1.0, 0, 1); // top surface only; no_cull=1 since it's rotated
	dstate.draw_cube(qbds.untex_qbd, roof,  LT_GRAY, 0, 0.0,        0, 0, 0, 0, 1.0, 1.0, 1.0, 0, 1); // draw all sides; no_cull=1 since it's rotated
	// rotate into place
	float const angle(0.25*(dim ? -1.0 : 1.0)*(height/get_length()));
	vector3d const axis(vector_from_dim_dir(!dim, dir));
	point about;
	about[ dim] = bcube.d[dim][dir];
	about[!dim] = bcube.get_center_dim(!dim);
	about.z     = roof.zc();
	rotate_verts(qbds.qbd.verts,       axis, angle, about, vs1);
	rotate_verts(qbds.untex_qbd.verts, axis, angle, about, vs2);
	// draw 4 untextured metal legs
	for (cube_t const &leg : get_legs()) {dstate.draw_cube(qbds.untex_qbd, leg, WHITE, 0, 0.0, 4);} // skip top and bottom
}
bool parking_solar_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	for (cube_t const &leg : get_legs()) {
		if (sphere_cube_int_update_pos(pos_, radius_, (leg + xlate), p_last, 0, cnorm)) return 1;
	}
	return 0;
}
vect_cube_t const &parking_solar_t::get_legs() const {
	float const height(bcube.dz()), border(0.01*height), leg_width(0.05*height), leg_hwidth(0.5*leg_width);
	cube_t inner(bcube), outer(bcube);
	inner.expand_in_dim( dim, -(leg_width  + border));
	inner.expand_in_dim(!dim, -(leg_hwidth + border));
	outer.expand_by_xy(-border);
	outer.z2() -= 0.025*height; // end mid-roof
	unsigned spans_per_side(1), spaces_per_leg(num_spaces);

	while (spaces_per_leg > 8 && (spaces_per_leg & 1) == 0) {
		spans_per_side *= 2; // double the number of spans
		spaces_per_leg /= 2; // halve the number of spaces per leg
	}
	while (spaces_per_leg > 12 && (spaces_per_leg % 3) == 0) {
		spans_per_side *= 3;
		spaces_per_leg /= 3;
	}
	float const step(inner.get_sz_dim(!dim)/spans_per_side);
	static vect_cube_t legs; // reused across calls
	legs.clear();

	for (unsigned side = 0; side < 2; ++side) {
		for (unsigned i = 0; i <= spans_per_side; ++i) {
			cube_t leg(outer);
			leg.d[dim][!side] = inner.d[dim][side];
			set_wall_width(leg, (inner.d[!dim][0] + i*step), leg_hwidth, !dim);
			if (bool(side) ^ dir) {leg.z2() -= 0.25*bcube.dz();} // shorter sloped side
			legs.push_back(leg);
		}
	} // for side
	return legs;
}

// birds/pigeons

// pos is at the feet; ignore dir.z
city_bird_base_t::city_bird_base_t(point const &pos_, float height, vector3d const &dir_, unsigned model_id) :
	city_obj_t(pos_, 0.0), dir(vector3d(dir_.x, dir_.y, 0.0).get_norm())
{
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
	float const hheight(0.5*height), xy_radius(hheight*sz.xy_mag()/sz.z); // assume model can be rotated in XY and take the max bounds
	set_bcube_from_vcylin(pos, height, xy_radius);
	radius = hheight*sz.mag()/sz.z; // max diagonal
	pos.z += hheight; // pos is on the ground, while we want the bsphere to be at the center
}

void pigeon_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!dstate.is_visible_and_unoccluded(bcube, dist_scale)) return;
	building_obj_model_loader.draw_model(dstate.s, pos, bcube, dir, WHITE, dstate.xlate, OBJ_MODEL_PIGEON, shadow_only);
}

// signs (for buildings)

sign_t::sign_t(cube_t const &bcube_, bool dim_, bool dir_, string const &text_, colorRGBA const &bc, colorRGBA const &tc,
	bool two_sided_, bool emissive_, bool small_, bool scrolling_, bool fs, bool in_skyway_, cube_t const &conn) :
	oriented_city_obj_t(dim_, dir_), two_sided(two_sided_), emissive(emissive_), small(small_), scrolling(scrolling_),
	free_standing(fs), in_skyway(in_skyway_), bkg_color(bc), text_color(tc), connector(conn)
{
	assert(!text_.empty());
	bcube = frame_bcube = text_bcube = bcube_; // excludes connector
	if (free_standing && !connector.is_all_zeros()) {bcube.union_with_cube(connector);}
	set_bsphere_from_bcube(); // recompute bsphere from bcube
	text  = (scrolling ? " "+text_+" " : text_); // pad with space on both sides if scrolling
	if (is_hospital_sign()) {text_bcube.z2() -= 0.75*text_bcube.dz();} // top 75% is hospital image and bottom 25% is hospital text

	if (scrolling) { // precompute text character offsets for scrolling effect
		unsigned const text_len(text.size());
		float const width(text_bcube.get_sz_dim(!dim)); // make text area a bit wider to account for the space padding
		text_bcube.expand_in_dim(!dim, 0.25*(float(text_len)/float(text_len-2) - 1.0)*width);
		text_bcube.expand_in_dim( dim, 0.1*get_depth()); // expand outward a bit to reduce Z-fighting; doesn't seem to help much
		vector<vert_norm_tc_color> verts;
		add_sign_text_verts(text, text_bcube, dim, dir, text_color, verts, 0.0, 0.0, 1); // include_space_chars=1, since we need offsets for all characters
		assert(verts.size() == 4*text_len);
		char_pos.resize(text_len);
		float const start_val(verts.front().v[!dim]);
		for (unsigned n = 0; n < text_len; ++n) {char_pos[n] = (verts[4*n+2].v[!dim] - start_val);}
		assert(char_pos.back() != 0.0);
		float const pos_scale(1.0/char_pos.back()); // can be positive or negative; result of divide should always be positive
		for (unsigned n = 0; n < text_len; ++n) {char_pos[n] *= pos_scale;} // scale so that full width is 1.0
	}
}
/*static*/ void sign_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {text_drawer::bind_font_texture();}
	if (!shadow_only) {enable_blend();}
}
/*static*/ void sign_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {disable_blend();}
}
void sign_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (small && shadow_only) return; // small signs have no shadows
	if (in_skyway && dstate.camera_bs.z < draw_zmin) return; // camera below the skyway
	float const dmax(dist_scale*dstate.draw_tile_dist);
	if (small && !bcube.closest_dist_less_than(dstate.camera_bs, 0.4*dmax)) return; // 40% view dist
	static quad_batch_draw temp_qbd;

	if (emissive && bkg_color == RED) { // special case for emissive emergency room sign
		dstate.draw_cube(temp_qbd, frame_bcube, bkg_color); // untextured, matte back
		select_texture(WHITE_TEX);
		dstate.s.add_uniform_float("emissive_scale", 1.0); // 100% emissive
		temp_qbd.draw_and_clear();
		dstate.s.add_uniform_float("emissive_scale", 0.0); // reset
		text_drawer::bind_font_texture(); // restore
	}
	else {
		dstate.draw_cube(qbds.untex_qbd, frame_bcube, bkg_color); // untextured, matte back
	}
	if (!connector.is_all_zeros()) { // draw connector; is this needed for the shadow pass?
		dstate.draw_cube(qbds.untex_qbd, connector, LT_GRAY, 0, 0.0, (free_standing ? 4 : 0)); // untextured, matte; skip top and bottom if free standing
	}
	if (free_standing) {} // connector is the base and sign bcube is the top
	if (shadow_only) return; // no text or images in shadow pass

	if (is_hospital_sign() && bcube.closest_dist_less_than(dstate.camera_bs, 0.4*dmax)) { // special case for hospital
		select_texture(get_texture_by_name("roads/hospital_sign.png"));
		cube_t top_part(frame_bcube);
		top_part.z1() = text_bcube.z2(); // above the text
		top_part.expand_in_dim(dim, 0.1*frame_bcube.get_sz_dim(dim)); // expand slightly to prevent Z-fighting
		dstate.draw_cube(temp_qbd, top_part, WHITE, 1, 0.0, (dim ? 5 : 6)); // only draw faces in <dim>
		temp_qbd.draw_and_clear();
		text_drawer::bind_font_texture(); // restore
	}
	if (!(emissive && is_night()) && !bcube.closest_dist_less_than(dstate.camera_bs, 0.9*(small ? 0.4 : 1.0)*dmax)) return; // too far to see the text in daytime

	if (scrolling && animate2) { // at the moment we can only scroll in integer characters, 4 per second
		// add pos x/y so that signs scroll at different points per building; make sure it's positive
		double const scroll_val(0.25*tfticks/TICKS_PER_SECOND + fabs(pos.x) + fabs(pos.y));
		float const scroll_val_mod(scroll_val - floor(scroll_val)); // take the fractional part
		assert(!char_pos.empty());
		auto it(std::lower_bound(char_pos.begin(), char_pos.end(), scroll_val_mod));
		unsigned const offset(it - char_pos.begin());
		assert(it != char_pos.end()); // can't be >= 1.0
		float const lo((offset == 0) ? 0.0f : char_pos[offset-1]), hi(char_pos[offset]), width(hi - lo), remainder(scroll_val_mod - lo);
		assert(width > 0.0);
		assert(remainder >= 0.0);
		assert(remainder <= width);
		float const first_char_clip_val(remainder/width), last_char_clip_val(1.0 - first_char_clip_val);
		string scroll_text(text);
		std::rotate(scroll_text.begin(), scroll_text.begin()+offset, scroll_text.end());
		scroll_text.push_back(scroll_text.front()); // duplicate first character for partial character scroll effect
		draw_text(dstate, qbds, scroll_text, first_char_clip_val, last_char_clip_val);
	}
	else {draw_text(dstate, qbds, text);}
}
void sign_t::draw_text(draw_state_t &dstate, city_draw_qbds_t &qbds, string const &text_to_draw, float first_char_clip_val, float last_char_clip_val) const {
	quad_batch_draw &qbd((emissive /*&& is_night()*/) ? qbds.emissive_qbd : qbds.qbd);
	bool const front_facing(((camera_pdu.pos[dim] - dstate.xlate[dim]) < text_bcube.d[dim][dir]) ^ dir);
	if (front_facing  ) {add_sign_text_verts(text_to_draw, text_bcube, dim,  dir, text_color, qbd.verts, first_char_clip_val, last_char_clip_val);} // draw the front side text
	else if (two_sided) {add_sign_text_verts(text_to_draw, text_bcube, dim, !dir, text_color, qbd.verts, first_char_clip_val, last_char_clip_val);} // draw the back  side text
}

stopsign_t::stopsign_t(point const &pos_, float height, float width, bool dim_, bool dir_, unsigned num_way_) :
	oriented_city_obj_t(pos_, max(width, height), dim_, dir_), num_way(num_way_)
{
	assert(num_way == 3 || num_way == 4);
	bcube.set_from_point(pos);
	bcube.expand_in_dim( dim, 0.05*width); // thickness
	bcube.expand_in_dim(!dim, 0.50*width); // width
	bcube.z2() += height;
}
cube_t stopsign_t::get_bird_bcube() const {
	cube_t top_place(bcube);
	if (SIGN_STOPSIGN_HEIGHT > 1.0) {top_place.z2() = top_place.z1() + SIGN_STOPSIGN_HEIGHT*bcube.dz();} // extend to top of street sign
	top_place.translate_dim(dim, (dir ? -1.0 : 1.0)*0.07*bcube.dz()); // translate to match street sign pole in draw_stoplights_and_street_signs()
	return top_place;
}
/*static*/ void stopsign_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	// Note: the texture is still needed in the shadow pass as an alpha mask
	int tid(-1);
	if      (dstate.pass_ix == 0) {tid = get_texture_by_name("roads/stop_sign.png",     0, 0, 0);} // wrap_mir=0
	else if (dstate.pass_ix == 1) {tid = get_texture_by_name("roads/white_octagon.png", 0, 0, 0, 0.0, 1, 1, 4, 1);} // wrap_mir=0, is_alpha_mask=1
	else if (!shadow_only)        {tid = get_texture_by_name("roads/stop_4_way.jpg",    0, 0, 0);} // wrap_mir=0
	if (tid >= 0) {select_texture(tid);}
	if (!shadow_only) {dstate.s.add_uniform_float("min_alpha", 0.25);} // fix mipmap drawing
}
/*static*/ void stopsign_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA);}
}
void stopsign_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	float const width(get_width()), thickness(get_depth()), sign_back(bcube.d[dim][dir] + (dir ? -1.0 : 1.0)*0.1*thickness);
	// draw the octagon with one side textured and the other not, in two passes
	bool const front_facing(((camera_pdu.pos[dim] - dstate.xlate[dim]) < bcube.d[dim][dir]) ^ dir);
	unsigned const skip_dims((1 << (1-dim)) | 4); // skip edges and top/bottom
	
	if (unsigned(!front_facing) == dstate.pass_ix) { // pass 0: front; pass 1: back
		cube_t sign(bcube);
		sign.z1() = bcube.z2() - width;
		sign.d[dim][!dir] = sign_back; // make it very thin
		dstate.draw_cube(qbds.qbd, sign, (front_facing ? WHITE : LT_GRAY), 0, 0.0, skip_dims);
	}
	if (num_way == 4 && dstate.pass_ix == 2) { // draw the 4-way sign part (only when close?)
		cube_t sign(bcube);
		set_cube_zvals(sign, (bcube.z2() - 1.3*width), (bcube.z2() - width)); // below the main octagon sign part
		sign.expand_in_dim(!dim, -0.2*width); // shrink
		sign.d[dim][!dir] = sign_back; // make it very thin
		dstate.draw_cube((front_facing ? qbds.qbd : qbds.untex_qbd), sign, (front_facing ? WHITE : LT_GRAY), 0, 0.0, skip_dims); // back side is untextured
	}
	if (dstate.pass_ix != 1) return; // no pole in this pass
	if (!shadow_only && !bcube.closest_dist_less_than(dstate.camera_bs, 0.4*dist_scale*dstate.draw_tile_dist)) return; // pole too far away to draw
	// draw the pole
	cube_t pole(bcube);
	pole.d[dim][dir] = sign_back; // fix for Z-fighting
	set_wall_width(pole, pos[!dim], 0.5*thickness, !dim);
	dstate.draw_cube(qbds.untex_qbd, pole, GRAY, 1); // skip_bottom=1
}

// city flags

int def_flag_tid(-1); // not yet loaded
int get_flag_texture(unsigned id);

city_flag_t::city_flag_t(cube_t const &flag_bcube_, bool dim_, bool dir_, point const &pole_base_, float pradius, int flag_id_) :
	oriented_city_obj_t(dim_, dir_), flag_bcube(flag_bcube_), pole_base(pole_base_), pole_radius(pradius), flag_id(flag_id_)
{
	draw_as_model = (is_horizontal() && building_obj_model_loader.is_model_valid(OBJ_MODEL_FLAG));

	if (draw_as_model) { // Note: pole_base and pole_radius are unused in this case
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FLAG)); // D, W, H
		flag_bcube.union_with_pt(pole_base); // make sure to include the pole
		// use bcube in !dim as the reference size; the flag sides face in dim and the flag sticks out in !dim
		float const sz_scale(flag_bcube.get_sz_dim(!dim)/sz.y), flag_z1(pole_base.z - 0.5*sz_scale*sz.z);
		set_wall_width(flag_bcube, flag_bcube.get_center_dim(dim), 0.5*sz_scale*sz.x, dim);
		set_cube_zvals(flag_bcube, flag_z1, (flag_z1 + sz_scale*sz.z)); // points upward starting at the pole base
		bcube = flag_bcube;
	}
	else { // calculate the bcube from the flag and pole
		bcube       = flag_bcube;
		bcube.z1()  = pole_base.z; // include pole Z-range
		bcube.z2() += 2.0*pole_radius; // account for the ball at the top
		bcube.union_with_pt(pole_base); // make sure to include the pole
		bcube.expand_by_xy(pradius); // include pole thickness
	}
	pos    = bcube.get_cube_center();
	radius = bcube.get_bsphere_radius();
}
/*static*/ void city_flag_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (shadow_only) return;
	def_flag_tid = get_texture_by_name("flags/american_flag_indexed.png"); // only if !draw_as_model (which we don't have here)?
	select_texture(def_flag_tid);
	disable_hemi_lighting_pre_post(dstate, shadow_only, 0);
}
void city_flag_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;

	if (draw_as_model) { // need to invert orient for dim=0 because dir applies to the side of the flag, not the direction of the pole
		if (dstate.is_occluded(bcube)) return;
		dstate.s.add_uniform_int("two_sided_lighting", 1); // must enable TSL because flags are not a closed volume
		building_obj_model_loader.draw_model(dstate.s, pos, bcube, get_orient_dir()*(dim ? 1.0 : -1.0), WHITE, dstate.xlate, OBJ_MODEL_FLAG, shadow_only);
		dstate.s.add_uniform_int("two_sided_lighting", 0);
		// def_flag_tid should still be set
		return;
	}
	bool const horizontal(is_horizontal());
	// draw the flag
	unsigned const skip_dims(4 + (1<<unsigned(!dim))); // top, bottom, and edges
	float const cview_dir((camera_pdu.pos[dim] - dstate.xlate[dim]) - pole_base[dim]);
	bool const visible_side((cview_dir < 0) ^ dir ^ dim), mirror_x(horizontal ? 1 : !visible_side), mirror_y(horizontal ? visible_side : 0);
	int tid(def_flag_tid);

	if (!shadow_only && flag_id >= 0) { // select custom flag texture
		int const flag_tid(get_flag_texture(flag_id));
		if (flag_tid >= 0) {tid = flag_tid;}
		
		if (tid != def_flag_tid && !qbds.qbd.empty()) { // flush verts drawn with default flag texture
			select_texture(def_flag_tid);
			qbds.qbd.draw_and_clear();
		}
	}
	dstate.draw_cube(qbds.qbd, flag_bcube, WHITE, 1, 0.0, skip_dims, mirror_x, mirror_y, horizontal); // swap_tc_xy=horizontal

	if (tid != def_flag_tid) {
		select_texture(tid);
		qbds.qbd.draw_and_clear(); // flush verts with this texture
		select_texture(def_flag_tid); // restore default flag texture
	}
	if (pole_radius == 0.0) return; // no pole to draw
	// draw the pole
	float const dmax(dist_scale*dstate.draw_tile_dist*(horizontal ? 0.7 : 1.0)); // horizontals flag poles are less visible
	if (!shadow_only && !bcube.closest_dist_less_than(dstate.camera_bs, 0.75*dmax)) return;
	unsigned const ndiv = 16;
	float const sphere_radius((horizontal ? 1.5 : 1.0)*pole_radius);
	point ce[2] = {pole_base, pole_base};
	if (horizontal) {ce[1][!dim] = bcube.d[!dim][dir] + (dir ? 1.0 : -1.0)*sphere_radius;}
	else {ce[1].z = bcube.z2() - sphere_radius;} // vertical pole
	add_cylin_as_tris(qbds.untex_qbd.verts, ce, pole_radius, (horizontal ? 1.0 : 0.5)*pole_radius, WHITE, ndiv, 0); // (truncated, if vertical) cone, sides only
	// draw the gold sphere at the top
	//if (shadow_only) return; // too small to cast a shadow?
	if (!shadow_only && !bcube.closest_dist_less_than(dstate.camera_bs, 0.4*dmax)) return;
	color_wrapper const cw(GOLD);
	dstate.temp_verts.clear();
	get_sphere_triangles(dstate.temp_verts, ce[1], sphere_radius, ndiv);
	for (vert_wrap_t const &v : dstate.temp_verts) {qbds.untex_qbd.verts.emplace_back(v.v, (v.v - ce[1]).get_norm(), 0.0, 0.0, cw);}
}
bool city_flag_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (sphere_cube_int_update_pos(pos_, radius_, (flag_bcube + xlate), p_last, 0, cnorm)) return 1; // flag coll
	if (pole_radius == 0.0 || draw_as_model) return 0; // no pole, skip
	return sphere_city_obj_cylin_coll(pole_base, pole_radius, pos_, p_last, radius_, xlate, cnorm); // pole coll
}

// park paths

void park_path_t::calc_bcube_bsphere() {
	assert(pts.size() >= 2);
	bcube.set_from_points(pts);
	bcube.expand_by_xy(hwidth); // conservative
	bcube.z2() += hwidth; // make sure it's nonzero height
	set_bsphere_from_bcube();
}
/*static*/ void park_path_t::pre_draw(draw_state_t &dstate, bool shadow_only) { // Note: not drawn in shadow pass
	select_texture(get_texture_by_name("roads/concrete.jpg"));
}
void park_path_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	assert(pts.size() >= 2);
	float const tscale(0.5/hwidth);
	vector3d prev_ortho;
	point prev_lo, prev_hi;
	vert_norm_tc_color vert;
	vert.set_norm(plus_z); // common across all verts
	vert.set_c4(color);

	for (unsigned i = 0; i < pts.size(); ++i) {
		point const &cur(pts[i]);
		vector3d const ortho((i+1 == pts.size()) ? prev_ortho : cross_product((pts[i+1] - cur), plus_z).get_norm());
		vector3d const v_side(hwidth*((i == 0) ? ortho : (ortho + prev_ortho).get_norm()));
		point const lo(cur - v_side), hi(cur + v_side);

		if (i > 0) { // emit a quad
			point qpts[4] = {prev_hi, prev_lo, lo, hi};

			if (i == 1 || i+1 == pts.size()) { // clip edges to plot
				do_line_clip_xy(qpts[0], qpts[3], plot.d); // hi
				do_line_clip_xy(qpts[1], qpts[2], plot.d); // lo
			}
			for (unsigned n = 0; n < 4; ++n) {
				vert.v = qpts[n];
				vert.t[0] = tscale*(vert.v.x - bcube.x1());
				vert.t[1] = tscale*(vert.v.y - bcube.y1());
				qbds.qbd.verts.push_back(vert);
			}
		}
		prev_ortho = ortho;
		prev_lo    = lo;
		prev_hi    = hi;
	} // for i
}
bool park_path_t::check_cube_coll_xy(cube_t const &c) const { // conservative
	if (!bcube.intersects_xy(c)) return 0;
	assert(pts.size() >= 2);
	cube_t test_cube(c);
	test_cube.expand_by_xy(hwidth);

	for (unsigned i = 0; i+1 < pts.size(); ++i) {
		if (check_line_clip_xy(pts[i], pts[i+1], test_cube.d)) return 1;
	}
	return 0;
}
bool park_path_t::check_point_contains_xy(point const &p) const {
	if (!bcube.contains_pt_xy(p)) return 0;
	cube_t c; c.set_from_point(p); // zero area
	return check_cube_coll_xy(c);
}

