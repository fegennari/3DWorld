// 3D World - City Object Placement Header
// by Frank Gennari
// 08/14/21

#include "city_objects.h"
#include "tree_3dw.h" // for tree_placer_t

extern int display_mode;
extern unsigned max_unique_trees;
extern double camera_zh;
extern tree_placer_t tree_placer;
extern city_params_t city_params;
extern object_model_loader_t building_obj_model_loader;

unsigned const q2t_ixs[6] = {0,2,1,0,3,2}; // quad => 2 tris

bool check_city_building_line_coll_bs_any(point const &p1, point const &p2);

float get_power_pole_offset() {return 0.045*city_params.road_width;}


struct plot_divider_type_t {
	bool is_occluder, has_alpha_mask;
	int tid, nm_tid;
	float wscale, hscale, tscale; // width, height, and texture scales
	colorRGBA color, map_color;
	string tex_name, nm_tex_name;

	plot_divider_type_t(string const &tn, string const &nm_tn, float ws, float hs, float ts, bool ic, bool ham, colorRGBA const &c, colorRGBA const &mc) :
		is_occluder(ic), has_alpha_mask(ham), tid(-1), nm_tid(-1), wscale(ws), hscale(hs), tscale(ts), color(c), map_color(mc), tex_name(tn), nm_tex_name(nm_tn) {}
	colorRGBA get_avg_color() const {return ((tid >= 0) ? texture_color(tid) : map_color).modulate_with(color);}

	void pre_draw(bool shadow_only) {
		if (shadow_only && !has_alpha_mask) return; // not textured

		if (tid < 0 && !tex_name.empty()) {
			unsigned const ncolors(has_alpha_mask ? 4 : 3);
			int const use_mipmaps (has_alpha_mask ? 0 : 1); // disable mipmaps if this texture has an alpha mask because we need to keep binary alpha
			tid = get_texture_by_name(tex_name, 0, 0, 1, 4.0, 1, use_mipmaps, ncolors); // load/lookup texture if needed, 4.0 aniso
		}
		if (nm_tid < 0 && !nm_tex_name.empty()) {nm_tid = get_texture_by_name(nm_tex_name, 1);} // load/lookup texture if needed
		select_texture(tid);
		if (nm_tid >= 0) {select_multitex(nm_tid, 5);} // bind normal map if it was specified
	}
	void post_draw(bool shadow_only) {
		if (!shadow_only && nm_tid >= 0) {select_multitex(FLAT_NMAP_TEX, 5);} // restore default flat normal map
	}
};
enum {DIV_WALL=0, DIV_FENCE, DIV_HEDGE, DIV_CHAINLINK, DIV_NUM_TYPES}; // types of plot dividers, with end terminator

plot_divider_type_t plot_divider_types[DIV_NUM_TYPES] = {
	plot_divider_type_t("cblock2.jpg", "normal_maps/cblock2_NRM.jpg", 0.50, 2.5, 1.0, 1, 0, WHITE, GRAY    ), // wall
	plot_divider_type_t("fence.jpg",   "normal_maps/fence_NRM.jpg",   0.15, 2.0, 1.0, 1, 0, WHITE, LT_BROWN), // fence
	plot_divider_type_t("hedges.jpg",  "", 1.00, 1.6, 1.0, 0, 0, GRAY, GREEN), // hedge - too short to be an occluder
	plot_divider_type_t("roads/chainlink_fence.png", "", 0.02, 1.55, 8.0, 0, 1, WHITE, GRAY) // chainlink fence with alpha mask; can't be taller than other fence types in case it intersects
};

void add_house_driveways_for_plot(cube_t const &plot, vect_cube_t &driveways);


bool city_obj_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	return sphere_cube_int_update_pos(pos_, radius_, (bcube + xlate), p_last, 1, 0, cnorm);
}

// benches

void bench_t::calc_bcube() {
	bcube.set_from_point(pos);
	bcube.expand_by(vector3d((dim ? 0.32 : 1.0), (dim ? 1.0 : 0.32), 0.0)*radius);
	bcube.z2() += 0.85*radius; // set bench height
}
/*static*/ void bench_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture(FENCE_TEX);} // normal map?
}
void bench_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;

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
	float const dist_val(shadow_only ? 0.0 : p2p_dist(camera_pdu.pos, center)/dstate.draw_tile_dist);
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
	vector3d const normal(get_poly_norm(pts, 1)), delta((0.2f*scale.x)*normal); // thickness = 0.4
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

// tree planters

tree_planter_t::tree_planter_t(point const &pos_, float radius_, float height) : city_obj_t(pos_, radius_) {
	bcube.set_from_point(pos);
	bcube.expand_by_xy(radius);
	bcube.z2() += height;
}
/*static*/ void tree_planter_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture((dstate.pass_ix == 0) ? (int)DIRT_TEX : get_texture_by_name("roads/sidewalk.jpg"));}
}
void tree_planter_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
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

// trashcans

trashcan_t::trashcan_t(point const &pos_, float radius_, float height) : city_obj_t(pos_, radius_) {
	bcube.set_from_point(pos);
	bcube.expand_by_xy(radius);
	bcube.z2() += height;
}
/*static*/ void trashcan_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture(get_texture_by_name("roads/asphalt.jpg"));}
}
void trashcan_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	color_wrapper const cw(LT_BROWN);
	cube_t hole(bcube);
	hole.expand_by_xy(-0.2*bcube.get_size()); // shrink 20% on all XY sides
	cube_t sides[4] = {bcube, bcube, bcube, bcube}; // -X, +X, -Y, +Y
	sides[0].x2() = sides[2].x1() = sides[3].x1() = hole.x1();
	sides[1].x1() = sides[2].x2() = sides[3].x2() = hole.x2();
	sides[2].y2() = hole.y1();
	sides[3].y1() = hole.y2();
	float const tscale(10.0);

	for (unsigned d = 0; d < 2; ++d) {
		dstate.draw_cube(qbd, sides[d  ], cw, 1, tscale, 0); // X
		dstate.draw_cube(qbd, sides[d+2], cw, 1, tscale, 1); // Y, skip X dims
	}
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
	if (!shadow_only) {dstate.s.add_uniform_float("hemi_lighting_scale", 0.0);} // disable hemispherical lighting
}
/*static*/ void fire_hydrant_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.set_cur_color(WHITE);} // restore to default color
	if (!shadow_only) {dstate.s.add_uniform_float("hemi_lighting_scale", 0.5);} // set hemispherical lighting back to the default
	city_obj_t::post_draw(dstate, shadow_only);
}
void fire_hydrant_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const { // Note: qbds are unused
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;

	if (!shadow_only) {
		building_obj_model_loader.draw_model(dstate.s, pos, bcube, orient, WHITE, dstate.xlate, OBJ_MODEL_FHYDRANT, shadow_only);
	}
	else { // shadow pass: draw as a simple cylinder, untextured, top end only
		draw_fast_cylinder(point(pos.x, pos.y, bcube.z1()), point(pos.x, pos.y, bcube.z2()), 0.8*cylin_radius, 0.8*cylin_radius, 12, 0, 4); // ndiv=12
	}
}
bool fire_hydrant_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	point const pos2(pos + xlate);
	float const r_sum(cylin_radius + radius_);
	if (!dist_less_than(pos_, pos2, r_sum)) return 0; // use sphere/vert cylinder instead?
	// since this is a cylinder, and we're not supposed to stand on top of it, assume collision normal is in the XY plane
	vector3d const coll_norm(vector3d((pos_.x - pos2.x), (pos_.y - pos2.y), 0.0).get_norm());
	pos_ += coll_norm*(r_sum - p2p_dist(pos_, pos2)); // move away from pos2
	if (cnorm) {*cnorm = coll_norm;}
	return 1;
}

// substations

substation_t::substation_t(cube_t const &bcube_, bool dim_, bool dir_) : dim(dim_), dir(dir_) {
	bcube = bcube_;
	*((sphere_t *)this) = bcube.get_bsphere();
}
/*static*/ void substation_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.add_uniform_float("hemi_lighting_scale", 0.0);} // disable hemispherical lighting
}
/*static*/ void substation_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.add_uniform_float("hemi_lighting_scale", 0.5);} // set hemispherical lighting back to the default
}
void substation_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	vector3d orient(zero_vector);
	orient[dim] = (dir ? 1.0 : -1.0);
	building_obj_model_loader.draw_model(dstate.s, pos, bcube, orient, WHITE, dstate.xlate, OBJ_MODEL_SUBSTATION, shadow_only);
}

// plot dividers

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
void divider_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const {
	if (dstate.pass_ix == DIV_NUM_TYPES && type == DIV_CHAINLINK) { // add chainlink fence posts
		if (!dstate.check_cube_visible(bcube, 1.5*dist_scale)) return;
		float const length(bcube.get_sz_dim(!dim)), height(bcube.dz()), thickness(bcube.get_sz_dim(dim));
		float const post_hwidth(1.5*thickness), post_width(2.0*post_hwidth), top_width(1.5*thickness);
		unsigned const num_sections(ceil(0.3*length/height)), num_posts(num_sections + 1);
		float const post_spacing((length - post_width)/num_sections);
		color_wrapper cw(GRAY);
		cube_t post(bcube), top(bcube); // copy dim and Z values
		post.expand_in_dim(dim, 0.5f*(post_width - thickness)); // increase width to post_width
		top .expand_in_dim(dim, 0.5f*(top_width  - thickness));
		post.z2() += 0.025*height; // extend slightly above the top of the fence
		set_wall_width(top, bcube.z2(), 0.5f*top_width, 2); // set height

		// for now we draw posts as cubes rather than cylinders since it's faster and easier, because we can use the existing qbd
		for (unsigned i = 0; i < num_posts; ++i) { // add posts
			set_wall_width(post, (bcube.d[!dim][0] + post_hwidth + i*post_spacing), post_hwidth, !dim);
			dstate.draw_cube(qbd, post, cw, 1);
		}
		dstate.draw_cube(qbd, top, cw, 1);
		return;
	}
	if (type != dstate.pass_ix) return; // this type not enabled in this pass
	if (type == DIV_CHAINLINK) {dist_scale *= 0.5;} // less visible
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	assert(dstate.pass_ix < DIV_NUM_TYPES);
	plot_divider_type_t const &pdt(plot_divider_types[dstate.pass_ix]);
	dstate.draw_cube(qbd, bcube, color_wrapper(pdt.color), 1, pdt.tscale/bcube.dz(), skip_dims); // skip bottom, scale texture to match the height

	if (!shadow_only && type == DIV_HEDGE && (bcube + dstate.xlate).closest_dist_less_than(camera_pdu.pos, 0.25f*(X_SCENE_SIZE + Y_SCENE_SIZE))) {
		dstate.hedge_draw.add(bcube); // draw detailed leaves for nearby hedges
	}
}
bool divider_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	cube_t bcube_wide(bcube + xlate);
	bcube_wide.expand_in_dim(dim, max(0.0f, 0.5f*(0.5f*building_t::get_scaled_player_radius() - bcube.get_sz_dim(dim)))); // make sure it's at least half player radius in thickness
	return sphere_cube_int_update_pos(pos_, radius_, bcube_wide, p_last, 1, 0, cnorm);
}

void hedge_draw_t::create(cube_t const &bc) {
	bcube = bc - bc.get_cube_center(); // centered on the origin
	unsigned const target_num_leaves(40000);
	vector3d const sz(bcube.get_size());
	float const leaf_sz(0.05*sz.z), surf_area(sz.x*sz.y + 2.0f*sz.z*(sz.x + sz.y));
	float const side_areas[5] = {sz.y*sz.z, sz.y*sz.z, sz.x*sz.z, sz.x*sz.z, sz.x*sz.y};
	rand_gen_t rgen;
	quad_batch_draw qbd;
	qbd.verts.reserve(6*target_num_leaves);

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
			vector3d const binormal(cross_product(normal, tangent));
			qbd.add_quad_dirs(pos, leaf_sz*tangent, leaf_sz*binormal, WHITE, normal);
		} // for n
	} // for s
	num_verts = qbd.verts.size();
	create_and_upload(qbd.verts, 0, 1);
}
void hedge_draw_t::draw_and_clear(shader_t &s) {
	if (empty()) return;
	if (!vbo_valid()) {create(to_draw.front());}
	select_texture(get_texture_by_name("pine2.jpg"));
	enable_blend(); // slightly smoother, but a bit of background shows through
	s.add_uniform_float("min_alpha", 0.5);
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
		glDrawArrays(GL_TRIANGLES, 0, num_verts);
		fgPopMatrix();
	} // for c
	post_render();
	s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA); // restore to the default
	disable_blend();
	to_draw.clear();
}

// swimming pools

// passes: 0=in-ground walls, 1=in-ground water, 2=above ground sides, 3=above ground water
/*static*/ void swimming_pool_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {
		if      (dstate.pass_ix == 2) {select_texture(WHITE_TEX);} // sides/untextured
		else if (dstate.pass_ix == 0) {select_texture(get_texture_by_name("bathroom_tile.jpg"));} // walls
		else if (dstate.pass_ix == 1 || dstate.pass_ix == 3) {select_texture(get_texture_by_name("snow2.jpg"));} // water surface
		else {assert(0);}
	}
}
void swimming_pool_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const {
	if ((dstate.pass_ix > 1) ^ above_ground) return; // not drawn in this pass
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;

	if (above_ground) { // cylindrical; bcube should be square in XY
		point const camera_bs(camera_pdu.pos - dstate.xlate);
		float const radius(get_radius()), xc(bcube.xc()), yc(bcube.yc()), dscale(dist_scale*dstate.draw_tile_dist);
		unsigned const ndiv(shadow_only ? 24 : max(4U, min(64U, unsigned(6.0f*dscale/p2p_dist(camera_bs, pos)))));

		if (dstate.pass_ix == 2) { // draw sides, and maybe ladder
			dstate.s.set_cur_color(color);
			draw_fast_cylinder(point(xc, yc, bcube.z1()), point(xc, yc, bcube.z2()), radius, radius, ndiv, 0, 0); // untextured, no ends
			point const camera_bs(camera_pdu.pos - dstate.xlate);

			if (bcube.closest_dist_less_than(camera_bs, 0.5*dscale)) { // draw ladder
				unsigned const num_steps = 5;
				color_wrapper const step_color(LT_GRAY);
				cube_t ladder;
				float const side_pos(bcube.d[dim][dir]), swidth((dir ? 1.0 : -1.0)*0.065*radius); // ladder is on this side of the pool
				float const height(1.2*bcube.dz()), step_delta(height/(num_steps + 0.25)), step_offset(0.25*step_delta), step_height(0.14*step_delta);
				ladder.d[dim][!dir] = side_pos;
				ladder.d[dim][ dir] = side_pos + swidth;
				set_wall_width(ladder, (dim ? xc : yc), 0.16*radius, !dim);
				bool const is_close(bcube.closest_dist_less_than(camera_bs, 0.2*dscale));
				bool const is_very_close(is_close && bcube.closest_dist_less_than(camera_bs, 0.1*dscale));

				for (unsigned n = 0; n < num_steps; ++n) { // draw steps
					ladder.z1() = bcube .z1() + n*step_delta + step_offset;
					ladder.z2() = ladder.z1() + step_height;
					dstate.draw_cube(qbd, ladder, step_color, !is_very_close); // skip bottom if not close
				}
				if (is_close) { // draw bars
					float const bars_top(bcube.z1() + height), bar_radius(0.012*radius);
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
		else if (dstate.pass_ix == 3) { // draw water surface
			dstate.s.set_cur_color(wcolor);
			draw_circle_normal(0.0, radius, ndiv, 0, point(xc, yc, (bcube.z2() - 0.1f*bcube.dz()))); // shift slightly below the top
		}
	}
	else { // in-ground
		float const dz(bcube.dz()), wall_thick(1.2*dz), tscale(0.5/wall_thick);
		cube_t inner(bcube);
		inner.expand_by_xy(-wall_thick);

		if (dstate.pass_ix == 0) { // draw walls
			color_wrapper const cw(color);
			cube_t sides[4] = {bcube, bcube, bcube, bcube}; // {S, N, W center, E center}
			sides[0].y2() = sides[2].y1() = sides[3].y1() = inner.y1();
			sides[1].y1() = sides[2].y2() = sides[3].y2() = inner.y2();
			sides[2].x2() = inner.x1();
			sides[3].x1() = inner.x2();
			for (unsigned d = 0; d < 4; ++d) {dstate.draw_cube(qbd, sides[d], cw, 1, tscale, ((d > 2) ? 2 : 0));}
		}
		else if (dstate.pass_ix == 1) { // draw water surface
			inner.z2() -= 0.5*dz; // reduce water height by 50%; can't make water below the mesh though
			dstate.draw_cube(qbd, inner, color_wrapper(wcolor), 1, 0.5*tscale, 3); // draw top water
		}
	}
}
/*static*/ void swimming_pool_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.set_cur_color(WHITE);} // restore to default color
	city_obj_t::post_draw(dstate, shadow_only);
}
bool swimming_pool_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (above_ground) {
		if (!sphere_cube_intersect((pos_ - xlate), radius_, bcube)) return 0; // optimization
		float const radius(get_radius()), xc(bcube.xc() + xlate.x), yc(bcube.yc() + xlate.y), z1(bcube.z1() + xlate.z), z2(bcube.z2() + xlate.z);
		return sphere_vert_cylin_intersect(pos_, radius_, cylinder_3dw(point(xc, yc, z1), point(xc, yc, z2), radius, radius), cnorm); // checks sides
	}
	cube_t bcube_tall(bcube + xlate);
	bcube_tall.z2() += CAMERA_RADIUS + camera_zh; // extend upward so that player collision detection works better
	return sphere_cube_int_update_pos(pos_, radius_, bcube_tall, p_last, 1, 0, cnorm);
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
	occluder.expand_by_xy(pole_radius/SQRT2); // take inner radius to reduce the occluder side for the pedestrian
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
	select_multitex(get_texture_by_name("normal_maps/wood_NRM.jpg", 1), 5);
}
/*static*/ void power_pole_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_multitex(FLAT_NMAP_TEX, 5);} // restore to default
	city_obj_t::post_draw(dstate, shadow_only);
}

void add_cylin_as_tris(vector<vert_norm_tc_color> &verts, point const ce[2], float r1, float r2, color_wrapper const &cw,
	unsigned ndiv, unsigned draw_top_bot, float tst=1.0, float tss=1.0, bool swap_ts_tt=0)
{
	// added as individual triangles; would be more efficient to use indexed triangles
	vector3d v12; // will be plus_z
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
			vector3d const normal(d ? plus_z : -plus_z);
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
void power_pole_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, quad_batch_draw &untex_qbd, float dist_scale, bool shadow_only) const {
	point const camera_bs(camera_pdu.pos - dstate.xlate);
	float const dmax(shadow_only ? camera_pdu.far_ : dist_scale*dstate.draw_tile_dist);
	if (!bcube.closest_dist_less_than(camera_bs, dmax)) return;
	if (!camera_pdu.cube_visible((shadow_only ? bcube : bcube_with_wires) + dstate.xlate)) return;
	color_wrapper const black(BLACK), white(colorRGBA(0.7, 0.7, 0.7)), gray(colorRGBA(0.4, 0.4, 0.4)), cw(LT_BROWN); // darken the wood color
	bool const pole_visible(camera_pdu.cube_visible(bcube + dstate.xlate));
	float const wire_radius(get_wire_radius());
	cube_t tf_bcube;
	point conduit_top(all_zeros);

	if (pole_visible) {
		unsigned const ndiv(shadow_only ? 16 : max(4U, min(32U, unsigned(1.5f*dmax/p2p_dist(camera_bs, pos)))));
		unsigned const pole_ndiv(min(ndiv, 24U));
		float const vert_tscale = 10.0;
		point ce[2] = {base, get_top()};
		bool const draw_top(ce[1].z < camera_bs.z);
		add_cylin_as_tris(qbd.verts, ce, pole_radius, pole_radius, cw, pole_ndiv, (draw_top ? 2 : 0), vert_tscale, 1.0/pole_ndiv, 1); // swap_ts_tt=1

		if (dims == 3 && (shadow_only || bcube.closest_dist_less_than(camera_bs, 0.7*dmax))) { // draw transformer, untextured
			float const tf_radius(2.0*pole_radius), pole_height(bcube.dz()), tf_height(0.1*pole_height), y_sign(at_line_end[1] ? 1.0 : -1.0);
			ce[0].z = base.z  + 0.77*pole_height;
			ce[1].z = ce[0].z + tf_height;
			ce[0].x = ce[1].x = base.x;
			ce[0].y = ce[1].y = base.y + y_sign*(tf_radius + pole_radius); // offset in -y, +y at end (so that wires across above it)
			bool const draw_top_bot(camera_bs.z > 0.5f*(ce[0].z + ce[1].z));
			add_cylin_as_tris(untex_qbd.verts, ce, tf_radius, tf_radius, gray, ndiv, (draw_top_bot ? 2 : 1));
			tf_bcube.set_from_points(ce, 2);
			tf_bcube.expand_by_xy(tf_radius);

			if (!residential && ndiv > 4) { // draw conduit for wires that go into the ground
				float const cradius(3.5*wire_radius);
				conduit_top.assign(base.x, (base.y + y_sign*(0.5f*cradius + pole_radius)), (base.z + 0.6*pole_height)); // below the transformer
				point const cce[2] = {point(conduit_top.x, conduit_top.y, base.z), conduit_top};
				bool const draw_top(ndiv > 8 && camera_bs.z > conduit_top.z);
				add_cylin_as_tris(untex_qbd.verts, cce, cradius, cradius, gray, min(ndiv, 16U), (draw_top ? 2 : 0));
			}
		}
	}
	float const wire_spacing(get_hwire_spacing()), vwire_spacing(get_vwire_spacing());
	float const standoff_height(get_standoff_height()), standoff_radius(0.25*pole_radius);
	point wire_pts[3][2];
	unsigned wire_mask(0);
	bool drew_wires(0);

	for (unsigned d = 0; d < 2; ++d) {
		if (!has_dim_set(d)) continue; // no wires in this dim
		float const offsets[3] = {-wire_spacing, -0.3f*wire_spacing, wire_spacing}; // offset from the center to avoid intersecting the pole
		cube_t const cbar(calc_cbar(d));
		point p1;
		p1[d] = cbar.get_center_dim(d);
		p1.z  = cbar.z2(); // resting on top of the bar
		
		if (pole_visible && (shadow_only || cbar.closest_dist_less_than(camera_bs, 0.5*dmax))) { // could test cbar cube visible, but unclear if that's faster
			dstate.draw_cube(qbd, cbar, cw, 0, 0.8/cbar.dz()); // draw all sides

			if (!shadow_only && cbar.closest_dist_less_than(camera_bs, 0.15*dmax)) { // draw insulator standoffs
				unsigned verts_start(untex_qbd.verts.size()), verts_end(0);

				for (unsigned n = 0; n < 3; ++n) {
					p1[!d] = center[!d] + offsets[n]; // set wire offset
					wire_pts[n][d] = p1 + vector3d(0.0, 0.0, (standoff_height + wire_radius));

					if (n == 1 && d == 1) { // offset middle wire from standoff to avoid clipping through pole
						if      (!at_line_end[1]) {wire_pts[n][1].y = wire_pts[n][0].y;}
						else if (!at_line_end[0]) {wire_pts[n][0].x = wire_pts[n][1].x;}
						else {} // I guess it clips through the pole in this case
					}
					float const delta_offset(offsets[n] - offsets[0]);
					draw_vert_standoff(p1, camera_bs, standoff_height, standoff_radius, delta_offset, dmax, d, (n == 0), verts_start, verts_end, white, untex_qbd);
				} // for n
				wire_mask |= (1 << d); // mark wires as drawn in this dim
			}
			if (d == 1 && !tf_bcube.is_all_zeros()) { // connect wire to transformer if running in y dim
				float const spacing(0.3*tf_bcube.get_sz_dim(!d));
				point const tf_top_center(tf_bcube.xc(), tf_bcube.yc(), tf_bcube.z2());

				for (unsigned n = 0; n < 3; ++n) {
					point const pts[2] = {point(center.x+offsets[n], tf_top_center.y, p1.z+standoff_height+wire_radius),
						(tf_top_center + vector3d((n - 1.0)*spacing, 0.0, standoff_height-0.5f*wire_radius))}; // top wire, bottom transformer
					draw_wire(pts, wire_radius, black, untex_qbd);
				}
				if (!shadow_only && tf_bcube.closest_dist_less_than(camera_bs, 0.1*dmax)) { // draw insulator standoffs
					unsigned verts_start(untex_qbd.verts.size()), verts_end(0);

					for (unsigned n = 0; n < 3; ++n) {
						point const p2((tf_top_center.x + (n - 1.0)*spacing), tf_top_center.y, tf_top_center.z);
						draw_vert_standoff(p2, camera_bs, standoff_height, standoff_radius, n*spacing, dmax, d, (n == 0), verts_start, verts_end, white, untex_qbd);
					}
					wire_mask |= (1 << d); // mark wires as drawn in this dim
				}
			}
		}
		if (shadow_only) continue; // skip wires for shadow pass since they don't show up reliably
		bool const is_offset(center[!d] != base[!d]);
		float const sep_dist(0.5*get_power_pole_offset()), offset_sign(is_offset ? -1.0 : 1.0);
		float const bot_wire_zval(base.z + 0.75*bcube.dz()), bot_wire_pos(base[!d] + offset_sign*sep_dist);

		if (!at_line_end[d]) { // no wires at end pole
			cube_t wires_bcube(cbar);
			min_eq(wires_bcube.z1(), bot_wire_zval); // include lower wires
			UNROLL_2X(wires_bcube.d[d][i_] = bcube_with_wires.d[d][i_];)
			if (!wires_bcube.closest_dist_less_than(camera_bs, 0.45*dmax) || !camera_pdu.cube_visible(wires_bcube + dstate.xlate)) continue; // wires distance/VFC
			// draw the three top wires and one bottom wire in this dim
			p1[d] += wires_offset; // offset wires in case the pole was moved to avoid a driveway
			p1.z  += standoff_height + wire_radius; // resting on top of the standoff

			for (unsigned n = 0; n < 3; ++n) { // top wires
				p1[!d] = center[!d] + offsets[n]; // set wire spacing
				draw_ortho_wire(p1, wire_radius, pole_spacing[d], d, black, dstate, untex_qbd);
			}
			// bottom 3 wires: split the difference between the offset corner and street power poles and attach to different sides of the poles
			float const wire_extend(sep_dist - pole_radius - 0.5*cbar.get_sz_dim(d)); // extend to fill gap between outer wire at standoffs and wire ending at cbar
			point pw;
			pw[!d] = bot_wire_pos; // on the side of the pole
			pw.z   = bot_wire_zval;
			pw[d]  = p1[d] + wire_extend; // extend slightly to meet the crossing wire
			float const bot_wire_extend(pole_spacing[d] + sep_dist + wire_radius); // overlap with next wire if there is one, extend past standoff at last pole

			for (unsigned n = 0; n < 3; ++n) { // 3 bottom wires as well
				draw_ortho_wire(pw, wire_radius, bot_wire_extend, d, black, dstate, untex_qbd);
				pw.z -= vwire_spacing;
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
			unsigned verts_start(untex_qbd.verts.size()), verts_end(0);

			for (unsigned n = 0; n < 3; ++n) {
				if (n == 0) { // first standoff, draw a truncated cone
					draw_standoff_geom(ce, standoff_radius, dmax, camera_bs, white, untex_qbd);
					verts_end = untex_qbd.verts.size();
				}
				else { // next standoff, copy and translate the previous truncated cone
					for (unsigned v = verts_start; v < verts_end; ++v) {
						untex_qbd.verts.push_back(untex_qbd.verts[v]);
						untex_qbd.verts.back().v.z -= n*vwire_spacing;
					}
				}
			}
			if (d == 1 && !tf_bcube.is_all_zeros()) { // vertical wire up to transformer
				point const tf_conn_pt(pb.x, (tf_bcube.yc() - 0.5f*tf_bcube.dy()), (tf_bcube.z1() + 0.9*tf_bcube.dz())); // along the side near the top
				cube_t wire(tf_conn_pt, tf_conn_pt);
				wire.z1()  = pb.z + wire_radius - 2.0*vwire_spacing; // meet the top of the lowest wire
				wire.z2() += wire_radius; // span entire standoff
				wire.expand_by_xy(wire_radius);
				dstate.draw_cube(untex_qbd, wire, black, 1); // skip bottom
				// standoff for wire
				point ce[2] = {tf_conn_pt, tf_conn_pt};
				vector3d const so_dir(-1.0, 1.0, 0.0);
				ce[0] += 0.4*pole_radius*so_dir;
				ce[1] += 0.5*wire_radius*so_dir;
				draw_standoff_geom(ce, standoff_radius, dmax, camera_bs, white, untex_qbd);

				if (conduit_top != all_zeros) { // draw wires to the top of the ground conduit
					point const conn_pt(tf_conn_pt.x, tf_conn_pt.y, (pb.z - 2.0*vwire_spacing)); // connect to bottom wire
					point const pts[2] = {(conduit_top - vector3d(0.0, 0.0, wire_radius)), conn_pt};
					draw_wire(pts, wire_radius, black, untex_qbd);
				}
			}
		}
	} // for d
	if (drew_wires && wire_mask == 3 && bcube.closest_dist_less_than(camera_bs, 0.25*dmax)) { // both dims set, connect X and Y wires
		for (unsigned n = 0; n < 3; ++n) {draw_wire(wire_pts[n], wire_radius, black, untex_qbd);}
	}
	if (!shadow_only && !wires.empty() && bcube_with_wires.closest_dist_less_than(camera_bs, 0.3*dmax)) {
		for (auto &w : wires) { // represents all three wires tied together
			draw_wire(w.pts, wire_radius, black, untex_qbd);
			// draw vertical wire segment connecting the three, which also represents all three power wires tied together
			cube_t wire(w.pts[1], w.pts[1]); // connection point to bottom horizontal wires
			wire.expand_in_dim(2, vwire_spacing); // connect to wires above and below
			wire.expand_by_xy(wire_radius);
			dstate.draw_cube(untex_qbd, wire, black, 1, 0.0, 4); // skip top and bottom
			if (w.pole_base.z == w.pts[0].z || !dist_less_than(w.pts[0], camera_bs, 0.15*dmax)) continue; // no pole, or too far away
			point const ce[2] = {w.pole_base, (w.pts[0] + vector3d(0.0, 0.0, wire_radius))};
			float const radius(1.5f*wire_radius);
			unsigned const ndiv(max(4U, min(16U, unsigned(0.1f*dmax/p2p_dist(camera_bs, ce[1])))));
			bool const draw_top(camera_bs.z > 0.5f*(ce[0].z + ce[1].z));
			add_cylin_as_tris(untex_qbd.verts, ce, radius, radius, gray, ndiv, (draw_top ? 2 : 1));
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


bool city_obj_placer_t::gen_parking_lots_for_plot(cube_t plot, vector<car_t> &cars, unsigned city_id, unsigned plot_ix, vect_cube_t &bcubes, vect_cube_t &colliders, rand_gen_t &rgen) {
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
		parking_lots.push_back(park);
		bcubes.push_back(park); // add to list of blocker bcubes so that no later parking lots overlap this one
		//parking_lots.back().expand_by_xy(0.5*pad_dist); // re-add half the padding for drawing (breaks texture coord alignment)
		unsigned const nspaces(park.row_sz*park.num_rows);
		num_spaces += nspaces;

		// fill the parking lot with cars
		vector<unsigned char> &used_spaces(parking_lots.back().used_spaces);
		used_spaces.resize(nspaces, 0); // start empty
		car.dim = car_dim; car.dir = car_dir;
		point pos(corner_pos.x, corner_pos.y, plot.z2());
		pos[car_dim] += 0.5*dr + (car_dim ? 0.15 : -0.15)*fabs(dr); // offset for centerline, biased toward the front of the parking space
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
					car.set_bcube(pos, nom_car_size);
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

// non-const because this sets driveway_t::car_ix through add_car()
void city_obj_placer_t::add_cars_to_driveways(vector<car_t> &cars, vector<road_plot_t> const &plots, vector<vect_cube_t> &plot_colliders, unsigned city_id, rand_gen_t &rgen) {
	car_t car;
	car.park();
	car.cur_city = city_id;
	car.cur_road_type = TYPE_DRIVEWAY;
	vector3d const nom_car_size(city_params.get_nom_car_size()); // {length, width, height}

	for (auto i = driveways.begin(); i != driveways.end(); ++i) {
		if (rgen.rand_float() < 0.5) continue; // no car in this driveway 50% of the time
		car.cur_road = (unsigned short)i->plot_ix; // store plot_ix in road field
		car.cur_seg  = (unsigned short)(i - driveways.begin()); // store driveway index in cur_seg
		cube_t const &plot(plots[i->plot_ix]);
		car.dim = (i->y1() == plot.y1() || i->y2() == plot.y2()); // check which edge of the plot the driveway is connected to, which is more accurate than the aspect ratio
		if (i->get_sz_dim(car.dim) < 1.6*nom_car_size.x || i->get_sz_dim(!car.dim) < 1.25*nom_car_size.y) continue; // driveway is too small to fit this car
		car.dir = rgen.rand_bool(); // randomly pulled in vs. backed in, since we don't know the direction to the house anyway
		float const pad_l(0.75*nom_car_size.x), pad_w(0.6*nom_car_size.y); // needs to be a bit larger to fit trucks
		point cpos(0.0, 0.0, i->z2());
		cpos[ car.dim] = rgen.rand_uniform(i->d[ car.dim][0]+pad_l, i->d[ car.dim][1]-pad_l);
		cpos[!car.dim] = rgen.rand_uniform(i->d[!car.dim][0]+pad_w, i->d[!car.dim][1]-pad_w); // not quite centered
		car.set_bcube(cpos, nom_car_size);
		// check if this car intersects another parked car; this can only happen if two driveways intersect, which should be rare
		bool intersects(0);

		for (auto c = cars.rbegin(); c != cars.rend(); ++c) {
			if (c->cur_road != i->plot_ix) break; // prev plot, done
			if (car.bcube.intersects(c->bcube)) {intersects = 1; break;}
		}
		if (intersects) continue; // skip
		i->in_use = 2; // permanently in use
		cars.push_back(car);
		plot_colliders[i->plot_ix].push_back(car.bcube); // prevent pedestrians from walking through this parked car
	} // for i
}

bool check_pt_and_place_blocker(point const &pos, vect_cube_t &blockers, float radius, float blocker_spacing) {
	cube_t bc(pos);
	if (has_bcube_int_xy(bc, blockers, radius)) return 0; // intersects a building or parking lot - skip
	bc.expand_by_xy(blocker_spacing);
	blockers.push_back(bc); // prevent trees and benches from being too close to each other
	return 1;
}
bool try_place_obj(cube_t const &plot, vect_cube_t &blockers, rand_gen_t &rgen, float radius, float blocker_spacing, unsigned num_tries, point &pos) {
	for (unsigned t = 0; t < num_tries; ++t) {
		pos = rand_xy_pt_in_cube(plot, radius, rgen);
		if (check_pt_and_place_blocker(pos, blockers, radius, blocker_spacing)) return 1; // success
	}
	return 0;
}
void place_tree(point const &pos, float radius, int ttype, vect_cube_t &colliders, vector<point> &tree_pos, bool allow_bush, bool is_sm_tree) {
	tree_placer.add(pos, 0, ttype, allow_bush, is_sm_tree); // use same tree type
	cube_t bcube; bcube.set_from_sphere(pos, 0.15*radius); // use 15% of the placement radius for collision (trunk + planter)
	bcube.z2() += radius; // increase cube height
	colliders.push_back(bcube);
	tree_pos.push_back(pos);
}

void city_obj_placer_t::place_trees_in_plot(road_plot_t const &plot, vect_cube_t &blockers,
	vect_cube_t &colliders, vector<point> &tree_pos, rand_gen_t &rgen, unsigned buildings_end)
{
	if (city_params.max_trees_per_plot == 0) return;
	float const radius(city_params.tree_spacing*city_params.get_nom_car_size().x); // in multiples of car length
	float const spacing(max(radius, get_min_obj_spacing())), radius_exp(2.0*spacing);
	vector3d const plot_sz(plot.get_size());
	if (min(plot_sz.x, plot_sz.y) < 2.0*radius_exp) return; // plot is too small for trees of this size
	unsigned num_trees(city_params.max_trees_per_plot);
	if (plot.is_park) {num_trees += (rgen.rand() % city_params.max_trees_per_plot);} // allow up to twice as many trees in parks
	assert(buildings_end <= blockers.size());
	// shrink non-building blockers (parking lots, driveways, fences, walls, hedges) to allow trees to hang over them; okay if they become denormalized
	unsigned const input_blockers_end(blockers.size());
	float const non_buildings_overlap(0.7*radius);
	for (auto i = blockers.begin()+buildings_end; i != blockers.end(); ++i) {i->expand_by_xy(-non_buildings_overlap);}

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
	for (auto i = blockers.begin()+buildings_end; i != blockers.begin()+input_blockers_end; ++i) {i->expand_by_xy(non_buildings_overlap);} // undo initial expand
}

template<typename T> void city_obj_groups_t::add_obj(T const &obj, vector<T> &objs) {
	by_tile[get_tile_id_for_cube(obj.bcube)].push_back(objs.size());
	objs.push_back(obj);
}
template<typename T> void city_obj_groups_t::create_groups(vector<T> &objs, cube_t &all_objs_bcube) {
	vector<T> new_objs;
	new_objs.reserve(objs.size());
	reserve(by_tile.size()); // the number of actual groups

	for (auto g = by_tile.begin(); g != by_tile.end(); ++g) {
		unsigned const group_start(new_objs.size());
		cube_with_ix_t group;

		for (auto i = g->second.begin(); i != g->second.end(); ++i) {
			assert(*i < objs.size());
			group.assign_or_union_with_cube(objs[*i].get_outer_bcube());
			new_objs.push_back(objs[*i]);
		}
		sort(new_objs.begin()+group_start, new_objs.end());
		group.ix = new_objs.size();
		push_back(group);
		all_objs_bcube.assign_or_union_with_cube(group);
	} // for g
	objs.swap(new_objs);
	by_tile.clear(); // no longer needed
}

// Note: blockers are used for placement of objects within this plot; colliders are used for pedestrian AI
void city_obj_placer_t::place_detail_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders,
	vector<point> const &tree_pos, rand_gen_t &rgen, bool is_residential, bool have_streetlights)
{
	float const car_length(city_params.get_nom_car_size().x); // used as a size reference for other objects

	// place fire_hydrants if the model has been loaded; don't add fire hydrants in parks
	if (!plot.is_park && building_obj_model_loader.is_model_valid(OBJ_MODEL_FHYDRANT)) {
		// we want the fire hydrant on the edge of the sidewalk next to the road, not next to the plot; this makes it outside the plot itself
		float const radius(0.04*car_length), height(0.18*car_length), dist_from_road(-0.5*radius - get_sidewalk_width());
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
				fhydrant_groups.add_obj(fire_hydrant, fhydrants);
				colliders.push_back(fire_hydrant.bcube);
			} // for dir
		} // for dim
	}
	// place benches in parks and non-residential areas
	if (!is_residential || plot.is_park) {
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
			bench_groups.add_obj(bench, benches);
			colliders.push_back(bench.bcube);
		} // for n
	}
	// place planters; don't add planters in parks or residential areas
	if (!is_residential && !plot.is_park) {
		float const planter_height(0.05*car_length), planter_radius(0.25*car_length);

		for (auto i = tree_pos.begin(); i != tree_pos.end(); ++i) {
			planter_groups.add_obj(tree_planter_t(*i, planter_radius, planter_height), planters); // no colliders for planters; pedestrians avoid the trees instead
		}
	}
	// place trashcans next to sidewalks in commercial cities
	if (!is_residential) {
	
	}
	// place power poles if there are houses or streetlights
	point corner_pole_pos(all_zeros);

	if ((is_residential && have_city_buildings()) || have_streetlights) {
		float const road_width(city_params.road_width), pole_radius(0.015*road_width), height(0.9*road_width);
		float const xspace(plot.dx() + road_width), yspace(plot.dy() + road_width); // == city_params.road_spacing?
		float const xyspace[2] = {0.5f*xspace, 0.5f*yspace};
		// we can move in toward the plot center so that they don't block pedestrians, but then they can block driveways;
		// if we move them into the road, then they block traffic light crosswalks;
		// so we move them toward the road in an assymetic way and allow the pole to be not centered with the wires
		float const offset(0.075*road_width), extra_offset(get_power_pole_offset()); // assymmetric offset to match the aspect ratio of stoplights
		float const pp_x(plot.x2() + offset + extra_offset), pp_y(plot.y2() + offset);
		point pts[3]; // one on the corner and two on each side: {corner, x, y}
		for (unsigned i = 0; i < 3; ++i) {pts[i].assign(pp_x, pp_y, plot.z2());} // start at plot upper corner
		pts[1].x -= 0.5*xspace;
		pts[2].y -= 0.5*yspace;
		unsigned const dims[3] = {3, 1, 2};
		unsigned const pp_start(ppoles.size());

		for (unsigned i = 0; i < 3; ++i) {
			point pos(pts[i]);
			float wires_offset(0.0);

			if (i > 0 && !driveways.empty()) {
				bool const dim(i == 2);
				float const prev_val(pos[dim]);
				move_to_not_intersect_driveway(pos, (pole_radius + get_sidewalk_width()), dim);
				wires_offset = prev_val - pos[dim];
			}
			point base(pos);
			if (i == 1) {base.y += extra_offset;} // shift the pole off the sidewalk and off toward the road to keep it out of the way of pedestrians
			bool const at_line_end[2] = {0, 0};
			bool const at_grid_edge(plot.xpos+1U == num_x_plots || plot.ypos+1U == num_y_plots);
			ppole_groups.add_obj(power_pole_t(base, pos, pole_radius, height, wires_offset, xyspace, dims[i], at_grid_edge, at_line_end, is_residential), ppoles);
			if (i == 0) {corner_pole_pos = base;}
		}
		if (plot.xpos == 0) { // no -x neighbor plot, but need to add the power poles there
			unsigned const pole_ixs[2] = {0, 2};

			for (unsigned i = 0; i < 2; ++i) {
				point pt(pts[pole_ixs[i]]);
				pt.x -= xspace;
				bool const at_line_end[2] = {1, 0};
				ppole_groups.add_obj(power_pole_t(pt, pt, pole_radius, height, 0.0, xyspace, dims[pole_ixs[i]], 1, at_line_end, is_residential), ppoles);
			}
		}
		if (plot.ypos == 0) { // no -y neighbor plot, but need to add the power poles there
			unsigned const pole_ixs[2] = {0, 1};

			for (unsigned i = 0; i < 2; ++i) {
				point pt(pts[pole_ixs[i]]);
				pt.y -= yspace;
				point base(pt);
				if (i == 1) {base.y += extra_offset;}
				bool const at_line_end[2] = {0, 1};
				ppole_groups.add_obj(power_pole_t(base, pt, pole_radius, height, 0.0, xyspace, dims[pole_ixs[i]], 1, at_line_end, is_residential), ppoles);
			}
		}
		if (plot.xpos == 0 && plot.ypos == 0) { // pole at the corner of the grid
			point pt(pts[0]);
			pt.x -= xspace;
			pt.y -= yspace;
			bool const at_line_end[2] = {1, 1};
			ppole_groups.add_obj(power_pole_t(pt, pt, pole_radius, height, 0.0, xyspace, dims[0], 1, at_line_end, is_residential), ppoles);
		}
		for (auto i = (ppoles.begin() + pp_start); i != ppoles.end(); ++i) {colliders.push_back(i->get_ped_occluder());}
	}
	// place substations in commercial cities, near the corner pole that routes power into the ground, if the model has been loaded
	if (!is_residential && corner_pole_pos != all_zeros && building_obj_model_loader.is_model_valid(OBJ_MODEL_SUBSTATION)) {
		bool const dim(0), dir(0); // hard-coded for now
		float const ss_height(0.08*city_params.road_width), dist_from_corner(0.12); // distance from corner relative to plot size
		vector3d const ss_center((1.0 - dist_from_corner)*corner_pole_pos + dist_from_corner*plot.get_cube_center());
		vector3d const model_sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SUBSTATION));
		vector3d bcube_exp;
		bcube_exp[ dim] = 0.5*ss_height*model_sz.x/model_sz.z;
		bcube_exp[!dim] = 0.5*ss_height*model_sz.y/model_sz.z;
		cube_t ss_bcube(ss_center, ss_center);
		ss_bcube.expand_by_xy(bcube_exp);
		ss_bcube.z2() += ss_height;
		
		if (!has_bcube_int_xy(ss_bcube, blockers, 0.2*ss_height)) { // skip if intersects a building or parking lot
			sstation_groups.add_obj(substation_t(ss_bcube, dim, dir), sstations);
			colliders.push_back(ss_bcube);
			blockers.push_back(ss_bcube);
		}
	}
}

bool is_placement_blocked(cube_t const &cube, vect_cube_t const &blockers, cube_t const &exclude, unsigned prev_blockers_end, float expand, bool exp_dim) {
	cube_t query_cube(cube);
	query_cube.expand_in_dim(exp_dim, expand);

	for (auto b = blockers.begin(); b != blockers.end()+prev_blockers_end; ++b) {
		if (*b != exclude && b->intersects_xy_no_adj(query_cube)) return 1;
	}
	return 0;
}

// dim=narrow dimension of fence; dir=house front/back dir in dimension dim; side=house left/right side in dimension !dim
float extend_fence_to_house(cube_t &fence, cube_t const &house, float fence_hwidth, float fence_height, bool dim, bool dir, bool side) {
	float &fence_end(fence.d[!dim][!side]);
	fence_end = house.d[!dim][side]; // adjacent to the house
	set_wall_width(fence, house.d[dim][dir], fence_hwidth, dim);
	// try to expand to the wall edge of two part houses by doing a line intersection query
	point p1, p2;
	p1.z     = p2.z    = fence.z1() + 0.25*fence_height; // slightly up from the bottom edge of the fence
	p1[ dim] = p2[dim] = fence.d[dim][!dir]; // use the side that overlaps the house bcube
	p1[!dim] = fence_end - (side ? -1.0 : 1.0)*fence_hwidth; // pull back slightly so that the start point isn't exactly at the house edge
	p2[!dim] = house.d[!dim][!side]; // end point is the opposite side of the house
	point p_int;
	if (!check_city_building_line_coll_bs(p1, p2, p_int)) return 0.0; // if this fails, house bcube must be wrong; should this be asserted?
	float const dist(fabs(fence_end - p_int[!dim]));
	fence_end = p_int[!dim];
	assert(fence.is_strictly_normalized());
	return dist;
}

void city_obj_placer_t::place_residential_plot_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, rand_gen_t &rgen) {
	assert(plot_subdiv_sz > 0.0);
	sub_plots.clear();
	if (plot.is_park) return; // no dividers in parks
	subdivide_plot_for_residential(plot, plot_subdiv_sz, 0, sub_plots); // parent_plot_ix=0, not needed
	if (sub_plots.size() <= 1) return; // nothing to divide
	if (rgen.rand_bool()) {std::reverse(sub_plots.begin(), sub_plots.end());} // reverse half the time so that we don't prefer a divider in one side or the other
	unsigned const shrink_dim(rgen.rand_bool()); // mostly arbitrary, could maybe even make this a constant 0
	float const sz_scale(0.06*city_params.road_width);
	unsigned const dividers_start(dividers.size()), prev_blockers_end(blockers.size());
	float const min_pool_spacing_to_plot_edge(0.5*city_params.road_width);
	colorRGBA const pool_side_colors[5] = {WHITE, WHITE, GRAY, LT_BROWN, LT_BLUE};
	vect_cube_t bcubes;

	for (auto i = sub_plots.begin(); i != sub_plots.end(); ++i) {
		// place plot dividers
		unsigned const type(rgen.rand()%DIV_NUM_TYPES); // use a consistent divider type for all sides of this plot
		if (type == DIV_CHAINLINK) continue; // chain link fence is not a primary divider - no divider for this plot; also, can't place a swimming pool here because it's not enclosed
		// should we remove or move house fences for divided sub-plots? I'm not sure how that would actually be possible at this point; or maybe skip dividers if the house has a fence?
		plot_divider_type_t const &pdt(plot_divider_types[type]);
		float const hwidth(0.5*sz_scale*pdt.wscale), z2(i->z1() + sz_scale*pdt.hscale);
		float const shrink_border(1.5*get_inner_sidewalk_width()); // needed for pedestrians to move along the edge of the plot; slightly larger to prevent collisions
		float translate_dist[2] = {0.0, 0.0};
		unsigned const prev_dividers_end(dividers.size());
		cube_t place_area(plot);
		place_area.expand_by_xy(-shrink_border);

		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				float const div_pos(i->d[dim][dir]);
				if (div_pos == plot.d[dim][dir]) continue; // sub-plot is against the plot border, don't need to add a divider
				bool const back_of_plot(i->d[!dim][0] != plot.d[!dim][0] && i->d[!dim][1] != plot.d[!dim][1]); // back of the plot, opposite the street
				unsigned const skip_dims(0); // can't make this (back_of_plot ? (1<<(1-dim)) : 0) because the edge may be showing at borders of different divider types
				cube_t c(*i);
				c.intersect_with_cube_xy(place_area);
				c.z2() = z2;
				set_wall_width(c, div_pos, hwidth, dim); // centered on the edge of the plot

				if (dim == shrink_dim) {
					translate_dist[dir] = (dir ? -1.0 : 1.0)*hwidth;
					c.translate_dim(dim, translate_dist[dir]); // move inside the plot so that edges line up
					// clip to the sides to remove overlap; may not line up with a neighboring divider of a different type/width, but hopefully okay
					for (unsigned d = 0; d < 2; ++d) {
						if (c.d[!dim][d] != plot.d[!dim][d]) {c.d[!dim][d] -= (d ? 1.0 : -1.0)*hwidth;}
					}
				}
				else {
					c.expand_in_dim(!dim, -0.001*hwidth); // fix for z-fighting
				}
				if (!back_of_plot) { // check for overlap of other plot dividers to the left and right
					cube_t test_cube(c);
					test_cube.expand_by_xy(4.0*hwidth); // expand so that adjacency counts as intersection
					bool overlaps(0);

					for (auto d = (dividers.begin()+dividers_start); d != (dividers.begin()+prev_dividers_end) && !overlaps; ++d) {
						overlaps |= (d->dim == bool(dim) && test_cube.contains_pt_xy(d->bcube.get_cube_center()));
					}
					if (overlaps) continue; // overlaps a previous divider, skip this one
				}
				divider_groups.add_obj(divider_t(c, type, dim, dir, skip_dims), dividers);
				colliders.push_back(c);
				blockers .push_back(c);
			} // for dir
		} // for dim

		// place swimming pools
		if (!i->is_residential || i->is_park || i->street_dir == 0) continue; // not a residential plot along a road
		if (rgen.rand_float() < 0.1) continue; // add pools 90% of the time
		bcubes.clear();
		if (have_buildings()) {get_building_bcubes(*i, bcubes);}
		if (bcubes.empty()) continue; // no house, skip adding swimming pool
		assert(bcubes.size() == 1); // there should be exactly one building/house in this sub-plot
		bool const dim((i->street_dir-1)>>1), dir((i->street_dir-1)&1); // direction to the road
		bool const above_ground(rgen.rand_bool());
		cube_t const &house(bcubes.front());
		cube_t pool_area(*i);
		pool_area.d[dim][dir] = house.d[dim][!dir]; // limit the pool to the back yard

		for (unsigned d = 0; d < 2; ++d) {
			if (i->d[!dim][d] == plot.d[!dim][d]) {pool_area.d[!dim][d] = house.d[!dim][d];} // adjacent to road - constrain to house projection so that side fence can be placed
		}
		float const dmin(min(pool_area.dx(), pool_area.dy())); // or should this be based on city_params.road_width?
		if (dmin < 0.75f*city_params.road_width) continue; // back yard is too small to add a pool

		for (unsigned d = 0; d < 2; ++d) { // keep pools away from the edges of plots; applies to sub-plots on the corners
			max_eq(pool_area.d[d][0], plot.d[d][0]+min_pool_spacing_to_plot_edge);
			min_eq(pool_area.d[d][1], plot.d[d][1]-min_pool_spacing_to_plot_edge);
		}
		pool_area.expand_by_xy(-0.05*dmin); // small shrink to keep away from walls, fences, and hedges
		vector3d pool_sz;
		pool_sz.z = (above_ground ? rgen.rand_uniform(0.08, 0.12)*city_params.road_width : 0.01f*dmin);

		for (unsigned d = 0; d < 2; ++d) {
			pool_sz[d] = ((above_ground && d == 1) ? pool_sz[0] : rgen.rand_uniform(0.5, 0.7)*dmin); // above_ground_cylin pools have square bcubes
			pool_area.d[d][1] -= pool_sz[d]; // shrink so that pool_area is where (x1, x2) can be placed
		}
		if (!pool_area.is_normalized()) continue; // pool area is too small; this can only happen due to shrink at plot edges
		point pool_llc;
		pool_llc.z = i->z2();

		for (unsigned n = 0; n < 20; ++n) { // make some attempts to generate a valid pool location
			for (unsigned d = 0; d < 2; ++d) {pool_llc[d] = rgen.rand_uniform(pool_area.d[d][0], pool_area.d[d][1]);}
			cube_t pool(pool_llc, (pool_llc + pool_sz));
			if (has_bcube_int_xy(pool, blockers, 0.08*dmin)) continue; // intersects some other object
			float const grayscale(rgen.rand_uniform(0.7, 1.0));
			float const water_white_comp(rgen.rand_uniform(0.1, 0.3)), extra_green(rgen.rand_uniform(0.2, 0.5)), lightness(rgen.rand_uniform(0.5, 0.8));
			colorRGBA const color(above_ground ? pool_side_colors[rgen.rand()%5]: colorRGBA(grayscale, grayscale, grayscale));
			colorRGBA const wcolor(lightness*water_white_comp, lightness*(water_white_comp + extra_green), lightness);
			
			// add fences along the sides of the house to separate the back yard from the front yard; if fences can't be added, then don't add the pool either
			plot_divider_type_t const &fence_pdt(plot_divider_types[DIV_CHAINLINK]);
			float const fence_hwidth(0.5*sz_scale*fence_pdt.wscale), fence_height(sz_scale*fence_pdt.hscale), fence_z2(i->z1() + fence_height);
			float const expand(-1.5*sz_scale*plot_divider_types[DIV_HEDGE].wscale); // shrink by widest divider to avoid false intersection with orthogonal dividers
			cube_t subplot_shrunk(*i);
			// translate so that fences line up with dividers; inexact if different width dividers are on each side
			for (unsigned d = 0; d < 2; ++d) {subplot_shrunk.d[shrink_dim][d] += translate_dist[d];}
			subplot_shrunk.expand_by_xy(-hwidth); // shrink by half width of surrounding dividers
			divider_t fences[2];
			bool bad_fence_place(0);

			for (unsigned side = 0; side < 2; ++side) { // left/right
				bool fence_dim(dim);
				cube_t fence(subplot_shrunk);
				fence.z2() = fence_z2;

				if (i->d[!dim][side] == plot.d[!dim][side]) { // at the edge of the plot, wrap the fence around in the back yard instead
					fence_dim ^= 1;
					extend_fence_to_house(fence, house, fence_hwidth, fence_height, !dim, side, !dir);

					if (is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {
						// Note: can't safely move the fence to the middle of the house if it intersects the pooly because it may intersect or block a door
						if ((house.d[!dim][!side] > pool.d[!dim][side]) ^ side) {bad_fence_place = 1; break;} // fence at back of house does not contain the pool
						extend_fence_to_house(fence, house, fence_hwidth, fence_height, !dim, !side, !dir); // try the back edge of the house
						if (is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {bad_fence_place = 1; break;} // blocked by a driveway, etc.
					}
				}
				else { // at the front of the house
					float const ext_dist(extend_fence_to_house(fence, house, fence_hwidth, fence_height, dim, dir, side));

					// check if front fence position is bad, or fence extension is too long (may block off the porch)
					if (ext_dist > 0.33*house.get_sz_dim(!dim) || is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {
						extend_fence_to_house(fence, house, fence_hwidth, fence_height, dim, !dir, side); // try the back edge of the house
						if (is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {bad_fence_place = 1; break;} // blocked by a driveway, etc.
					}
				}
				fences[side] = divider_t(fence, DIV_CHAINLINK, fence_dim, dir, 0); // Note: dir is unused in divider_t so doesn't have to be set correctly
			} // for side
			if (bad_fence_place) continue; // failed to fence off the pool, don't place it here
			pool_groups.add_obj(swimming_pool_t(pool, color, wcolor, above_ground, dim, dir), pools);
			pool.z2() += 0.1*city_params.road_width; // extend upward to make a better collider
			colliders.push_back(pool);
			blockers .push_back(pool);

			for (unsigned side = 0; side < 2; ++side) {
				divider_t const &fence(fences[side]);
				assert(fence.bcube.is_strictly_normalized());
				divider_groups.add_obj(fence, dividers);
				colliders.push_back(fence.bcube);
				blockers .push_back(fence.bcube);
			}
			break; // success
		} // for n
	} // for i
}

void city_obj_placer_t::add_house_driveways(road_plot_t const &plot, vect_cube_t &temp_cubes, rand_gen_t &rgen, unsigned plot_ix) {
	cube_t plot_z(plot);
	plot_z.z1() = plot_z.z2() = plot.z2() + 0.0002*city_params.road_width; // shift slightly up to avoid Z-fighting
	temp_cubes.clear();
	add_house_driveways_for_plot(plot_z, temp_cubes);

	for (auto i = temp_cubes.begin(); i != temp_cubes.end(); ++i) {
		bool dim(0), dir(0);
		get_closest_dim_dir_xy(*i, plot, dim, dir);
		driveways.emplace_back(*i, dim, dir, plot_ix);
	}
}

template<typename T> void city_obj_placer_t::draw_objects(vector<T> const &objs, city_obj_groups_t const &groups,
	draw_state_t &dstate, float dist_scale, bool shadow_only, bool has_immediate_draw)
{
	if (objs.empty()) return;
	T::pre_draw(dstate, shadow_only);
	unsigned start_ix(0);
	assert(qbd.empty() && untex_qbd.empty());

	for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
		if (!dstate.check_cube_visible(*g, dist_scale)) continue; // VFC/distance culling for group
		if (has_immediate_draw) {dstate.begin_tile(g->get_cube_center(), 1, 1);} // must setup shader and tile shadow map before drawing
		assert(start_ix <= g->ix && g->ix <= objs.size());

		for (unsigned i = start_ix; i < g->ix; ++i) {
			T const &obj(objs[i]);
			if (dstate.check_sphere_visible(obj.pos, obj.get_bsphere_radius(shadow_only))) {obj.draw(dstate, qbd, untex_qbd, dist_scale, shadow_only);}
		}
		if (!qbd.empty() || !untex_qbd.empty() || !dstate.hedge_draw.empty()) { // we have something to draw
			if (!has_immediate_draw) {dstate.begin_tile(g->get_cube_center(), 1, 1);} // will_emit_now=1, ensure_active=1
			qbd.draw_and_clear(); // draw this group with current smap
			bool must_restore_state(!dstate.hedge_draw.empty());
			dstate.hedge_draw.draw_and_clear(dstate.s);

			if (!untex_qbd.empty()) {
				dstate.set_untextured_material();
				untex_qbd.draw_and_clear();
				dstate.unset_untextured_material();
				must_restore_state = 1;
			}
			if (!shadow_only && must_restore_state) {T::pre_draw(dstate, shadow_only);} // re-setup for next tile
		}
	} // for g
	T::post_draw(dstate, shadow_only);
}

void city_obj_placer_t::clear() {
	parking_lots.clear(); benches.clear(); planters.clear(); trashcans.clear(); fhydrants.clear(); sstations.clear();
	driveways.clear(); dividers.clear(); pools.clear(); ppoles.clear();
	bench_groups.clear(); planter_groups.clear(); trashcan_groups.clear(); fhydrant_groups.clear(); sstation_groups.clear();
	divider_groups.clear(); pool_groups.clear(); ppole_groups.clear();
	all_objs_bcube.set_to_zeros();
	num_spaces = filled_spaces = 0;
}

void city_obj_placer_t::gen_parking_and_place_objects(vector<road_plot_t> &plots, vector<vect_cube_t> &plot_colliders, vector<car_t> &cars,
	unsigned city_id, bool have_cars, bool is_residential, bool have_streetlights)
{
	// Note: fills in plots.has_parking
	//timer_t timer("Gen Parking Lots and Place Objects");
	vect_cube_t bcubes, temp_cubes; // blockers, driveways
	vector<point> tree_pos;
	rand_gen_t rgen, detail_rgen;
	rgen.set_state(city_id, 123);
	detail_rgen.set_state(3145739*(city_id+1), 1572869*(city_id+1));
	if (city_params.max_trees_per_plot > 0) {tree_placer.begin_block(0); tree_placer.begin_block(1);} // both small and large trees
	bool const add_parking_lots(have_cars && !is_residential && city_params.min_park_spaces > 0 && city_params.min_park_rows > 0);
	float const sidewalk_width(get_sidewalk_width());

	for (auto i = plots.begin(); i != plots.end(); ++i) { // calculate num_x_plots and num_y_plots; these are used for determining edge power poles
		max_eq(num_x_plots, i->xpos+1U);
		max_eq(num_y_plots, i->ypos+1U);
	}
	for (auto i = plots.begin(); i != plots.end(); ++i) {
		tree_pos.clear();
		bcubes.clear();
		get_building_bcubes(*i, bcubes);
		size_t const plot_id(i - plots.begin()), buildings_end(bcubes.size());
		assert(plot_id < plot_colliders.size());
		vect_cube_t &colliders(plot_colliders[plot_id]); // used for pedestrians
		if (add_parking_lots && !i->is_park) {i->has_parking = gen_parking_lots_for_plot(*i, cars, city_id, plot_id, bcubes, colliders, rgen);}
		unsigned const driveways_start(driveways.size());
		if (is_residential) {add_house_driveways(*i, temp_cubes, detail_rgen, plot_id);}

		// driveways become blockers for other placed objects; make sure they extend into the road so that they intersect any placed streetlights or fire hydrants
		for (auto j = driveways.begin()+driveways_start; j != driveways.end(); ++j) {
			cube_t dw(*j);

			for (unsigned d = 0; d < 2; ++d) {
				if      (dw.d[d][0] == i->d[d][0]) {dw.d[d][0] -= sidewalk_width;}
				else if (dw.d[d][1] == i->d[d][1]) {dw.d[d][1] += sidewalk_width;}
			}
			bcubes.push_back(dw);
		} // for j
		if (city_params.assign_house_plots && plot_subdiv_sz > 0.0) {place_residential_plot_objects(*i, bcubes, colliders, detail_rgen);} // before placing trees
		place_trees_in_plot (*i, bcubes, colliders, tree_pos, detail_rgen, buildings_end);
		place_detail_objects(*i, bcubes, colliders, tree_pos, detail_rgen, is_residential, have_streetlights);
	} // for i
	connect_power_to_buildings(plots);
	if (have_cars) {add_cars_to_driveways(cars, plots, plot_colliders, city_id, rgen);}
	for (auto i = plot_colliders.begin(); i != plot_colliders.end(); ++i) {sort(i->begin(), i->end(), cube_by_x1());}
	bench_groups   .create_groups(benches,   all_objs_bcube);
	planter_groups .create_groups(planters,  all_objs_bcube);
	trashcan_groups.create_groups(trashcans, all_objs_bcube);
	fhydrant_groups.create_groups(fhydrants, all_objs_bcube);
	sstation_groups.create_groups(sstations, all_objs_bcube);
	divider_groups .create_groups(dividers,  all_objs_bcube);
	pool_groups    .create_groups(pools,     all_objs_bcube);
	ppole_groups   .create_groups(ppoles,    all_objs_bcube);

	if (0) { // debug info printing
		cout << TXT(benches.size()) << TXT(bench_groups.size()) << TXT(planters.size()) << TXT(planter_groups.size()) << TXT(trashcan_groups.size())
			 << TXT(fhydrants.size())  << TXT(fhydrant_groups.size()) << TXT(sstations.size()) << TXT(sstation_groups.size()) << TXT(dividers.size())
			 << TXT(divider_groups.size()) << TXT(pools.size()) << TXT(pool_groups.size()) << TXT(ppoles.size()) << TXT(ppole_groups.size()) << endl;
	}
	if (add_parking_lots) {
		cout << "parking lots: " << parking_lots.size() << ", spaces: " << num_spaces << ", filled: " << filled_spaces << ", benches: " << benches.size() << endl;
	}
}

/*static*/ bool city_obj_placer_t::subdivide_plot_for_residential(cube_t const &plot, float plot_subdiv_sz, unsigned parent_plot_ix, vect_city_zone_t &sub_plots) {
	if (min(plot.dx(), plot.dy()) < city_params.road_width) return 0; // plot is too small to divide
	assert(plot_subdiv_sz > 0.0);
	unsigned ndiv[2] = {0,0};
	float spacing[2] = {0,0};

	for (unsigned d = 0; d < 2; ++d) {
		float const plot_sz(plot.get_sz_dim(d));
		ndiv   [d] = max(1U, unsigned(round_fp(plot_sz/plot_subdiv_sz)));
		spacing[d] = plot_sz/ndiv[d];
	}
	if (ndiv[0] >= 100 || ndiv[1] >= 100) return 0; // too many plots? this shouldn't happen, but failing here is better than asserting or generating too many buildings
	unsigned const max_floors(0); // 0 is unlimited
	if (sub_plots.empty()) {sub_plots.reserve(2*(ndiv[0] + ndiv[1]) - 4);}

	for (unsigned y = 0; y < ndiv[1]; ++y) {
		float const y1(plot.y1() + spacing[1]*y), y2((y+1 == ndiv[1]) ? plot.y2() : (y1 + spacing[1])); // last sub-plot must end exactly at plot y2

		for (unsigned x = 0; x < ndiv[0]; ++x) {
			if (x > 0 && y > 0 && x+1 < ndiv[0] && y+1 < ndiv[1]) continue; // interior plot, no road access, skip
			float const x1(plot.x1() + spacing[0]*x), x2((x+1 == ndiv[0]) ? plot.x2() : (x1 + spacing[0])); // last sub-plot must end exactly at plot x2
			cube_t const c(x1, x2, y1, y2, plot.z1(), plot.z2());
			sub_plots.emplace_back(c, 0.0, 0, 1, get_street_dir(c, plot), 1, parent_plot_ix, max_floors); // cube, zval, park, res, sdir, capacity, ppix, nf; will favor x-dim for corner plots
		}
	} // for y
	return 1;
}

bool city_obj_placer_t::connect_power_to_point(point const &at_pos, bool near_power_pole) {
	float dmax_sq(0.0);

	for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts to connect a to a power pole without intersecting a building
		float dmin_sq(0.0);
		unsigned best_pole(0);
		point best_pos;

		for (auto p = ppoles.begin(); p != ppoles.end(); ++p) {
			point const cur_pos(p->get_nearest_connection_point(at_pos, near_power_pole));
			if (cur_pos == at_pos) continue; // bad point
			float const dsq(p2p_dist_sq(at_pos, cur_pos));
			if (dsq <= dmax_sq) continue; // this pole was previously flagged as bad
			if (dmin_sq == 0.0 || dsq < dmin_sq) {best_pos = cur_pos; dmin_sq = dsq; best_pole = (p - ppoles.begin());}
		} // for p
		if (dmin_sq == 0.0) return 0; // failed (no power poles?)
		if (ppoles[best_pole].add_wire(at_pos, best_pos, near_power_pole)) return 1; // add a wire pole for houses
		dmax_sq = dmin_sq; // prevent this pole from being used in the next iteration
	} // for n
	return 0; // failed
}
void city_obj_placer_t::connect_power_to_buildings(vector<road_plot_t> const &plots) {
	if (plots.empty() || ppoles.empty() || !have_buildings()) return;
	cube_t all_plots_bcube(plots.front());
	for (auto p = plots.begin()+1; p != plots.end(); ++p) {all_plots_bcube.union_with_cube(*p);} // query all buildings in the entire city rather than per-plot
	vector<point> ppts;
	get_building_power_points(all_plots_bcube, ppts);
	for (auto p = ppts.begin(); p != ppts.end(); ++p) {connect_power_to_point(*p, 1);} // near_power_pole=1
}

void city_obj_placer_t::move_to_not_intersect_driveway(point &pos, float radius, bool dim) const {
	cube_t test_cube;
	test_cube.set_from_sphere(pos, radius);

	// Note: this could be accelerated by iterating by plot, but this seems to already be fast enough (< 1ms)
	for (auto d = driveways.begin(); d != driveways.end(); ++d) {
		if (!d->intersects_xy(test_cube)) continue;
		bool const dir((d->d[dim][1] - pos[dim]) < (pos[dim] - d->d[dim][0]));
		pos[dim] = d->d[dim][dir] + (dir ? 1.0 : -1.0)*0.1*city_params.road_width;
		break; // maybe we should check for an adjacent driveway, but that would be rare and moving could result in oscillation
	}
}
void city_obj_placer_t::move_and_connect_streetlights(streetlights_t &sl) {
	for (auto s = sl.streetlights.begin(); s != sl.streetlights.end(); ++s) {
		if (!driveways.empty()) { // move to avoid driveways
			bool const dim(s->dir.y == 0.0); // direction to move in
			move_to_not_intersect_driveway(s->pos, 0.25*city_params.road_width, dim);
		}
		if (!ppoles.empty()) { // connect power
			point top(s->get_lpos());
			top.z += 1.05f*streetlight_ns::light_radius*city_params.road_width; // top of light
			connect_power_to_point(top, 0); // near_power_pole=0 because it may be too far away
		}
	} // for s
}

void city_obj_placer_t::draw_detail_objects(draw_state_t &dstate, bool shadow_only) {
	if (!dstate.check_cube_visible(all_objs_bcube, 1.0)) return; // check bcube, dist_scale=1.0
	draw_objects(benches,   bench_groups,    dstate, 0.16, shadow_only, 0); // dist_scale=0.16
	draw_objects(trashcans, trashcan_groups, dstate, 0.12, shadow_only, 0); // dist_scale=0.16
	draw_objects(fhydrants, fhydrant_groups, dstate, 0.07, shadow_only, 1); // dist_scale=0.07, has_immediate_draw=1
	draw_objects(sstations, sstation_groups, dstate, 0.15, shadow_only, 1); // dist_scale=0.15, has_immediate_draw=1
	draw_objects(ppoles,    ppole_groups,    dstate, 0.20, shadow_only, 0); // dist_scale=0.20
			
	if (!shadow_only) { // low profile, not drawn in shadow pass
		for (dstate.pass_ix = 0; dstate.pass_ix < 2; ++dstate.pass_ix) { // {dirt, stone}
			draw_objects(planters, planter_groups, dstate, 0.1, shadow_only, 0); // dist_scale=0.1
		}
	}
	for (dstate.pass_ix = 0; dstate.pass_ix < 4; ++dstate.pass_ix) { // {in-ground walls, in-ground water, above ground sides, above ground water}
		if (shadow_only && dstate.pass_ix <= 1) continue; // only above ground pools are drawn in the shadow pass; water surface is drawn to prevent light leaks, but maybe should extend z1 lower
		float const dist_scales[4] = {0.1, 0.5, 0.3, 0.5};
		draw_objects(pools, pool_groups, dstate, dist_scales[dstate.pass_ix], shadow_only, (dstate.pass_ix > 1)); // final 2 passes use immediate draw rather than qbd
	}
	// Note: not the most efficient solution, as it required processing blocks and binding shadow maps multiple times
	for (dstate.pass_ix = 0; dstate.pass_ix <= DIV_NUM_TYPES; ++dstate.pass_ix) { // {wall, fence, hedge, chainlink fence, chainlink fence posts}
		if (dstate.pass_ix == DIV_CHAINLINK && shadow_only) continue; // chainlink fence not drawn in the shadow pass
		draw_objects(dividers, divider_groups, dstate, 0.2, shadow_only, 0); // dist_scale=0.2
	}
	dstate.pass_ix = 0; // reset back to 0
}

template<typename T> bool proc_vector_sphere_coll(vector<T> const &objs, city_obj_groups_t const &groups, point &pos,
	point const &p_last, float radius, vector3d const &xlate, vector3d *cnorm)
{
	point const pos_bs(pos - xlate);
	unsigned start_ix(0);

	for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
		if (!sphere_cube_intersect(pos_bs, radius, *g)) continue;
		assert(start_ix <= g->ix && g->ix <= objs.size());

		for (auto i = objs.begin()+start_ix; i != objs.begin()+g->ix; ++i) {
			if (i->proc_sphere_coll(pos, p_last, radius, xlate, cnorm)) return 1;
		}
	} // for g
	return 0;
}
bool city_obj_placer_t::proc_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, vector3d *cnorm) const { // pos in in camera space
	if (!sphere_cube_intersect(pos, (radius + p2p_dist(pos,p_last)), (all_objs_bcube + xlate))) return 0;
	if (proc_vector_sphere_coll(benches,   bench_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(trashcans, trashcan_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(fhydrants, fhydrant_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(sstations, sstation_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(dividers,  divider_groups,  pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(pools,     pool_groups,     pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(ppoles,    ppole_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	// Note: no coll with tree_planters because the tree coll should take care of it
	return 0;
}

template<typename T> void check_vector_line_intersect(vector<T> const &objs, city_obj_groups_t const &groups, point const &p1, point const &p2, float &t, bool &ret) {
	unsigned start_ix(0);

	for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
		if (!check_line_clip(p1, p2, g->d)) continue;
		assert(start_ix <= g->ix && g->ix <= objs.size());
		for (auto i = objs.begin()+start_ix; i != objs.begin()+g->ix; ++i) {ret |= check_line_clip_update_t(p1, p2, t, i->bcube);}
	}
}
bool city_obj_placer_t::line_intersect(point const &p1, point const &p2, float &t) const { // Note: nothing to do for parking lots or tree_planters; p1/p2 in world space
	if (!all_objs_bcube.line_intersects(p1, p2)) return 0;
	bool ret(0);
	check_vector_line_intersect(benches,   bench_groups,    p1, p2, t, ret);
	check_vector_line_intersect(trashcans, trashcan_groups, p1, p2, t, ret);
	check_vector_line_intersect(fhydrants, fhydrant_groups, p1, p2, t, ret); // check bounding cube; cylinder intersection may be more accurate, but likely doesn't matter much
	check_vector_line_intersect(sstations, sstation_groups, p1, p2, t, ret);
	check_vector_line_intersect(dividers,  divider_groups,  p1, p2, t, ret);
	check_vector_line_intersect(pools,     pool_groups,     p1, p2, t, ret);
	check_vector_line_intersect(ppoles,    ppole_groups,    p1, p2, t, ret); // inaccurate, could be customized if needed
	return ret;
}

template<typename T> bool check_city_obj_bcube_pt_xy_contain(city_obj_groups_t const &groups, vector<T> const &objs, point const &pos, unsigned &obj_ix) {
	unsigned start_ix(0);

	for (auto i = groups.begin(); i != groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy(pos)) continue;
		assert(start_ix <= i->ix && i->ix <= objs.size());

		for (auto b = objs.begin()+start_ix; b != objs.begin()+i->ix; ++b) {
			if (pos.x < b->bcube.x1()) break; // objects are sorted by x1, none after this can match
			if (b->bcube.contains_pt_xy(pos)) {obj_ix = (b - objs.begin()); return 1;}
		}
	} // for i
	return 0;
}
bool city_obj_placer_t::get_color_at_xy(point const &pos, colorRGBA &color, bool skip_in_road) const {
	unsigned start_ix(0), obj_ix(0);
	if (check_city_obj_bcube_pt_xy_contain(bench_groups, benches, pos, obj_ix)) {color = texture_color(FENCE_TEX); return 1;}
	float const expand(0.15*city_params.road_width), x_test(pos.x + expand); // expand to approx tree diameter

	for (auto i = planter_groups.begin(); i != planter_groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy_exp(pos, expand)) continue;
		assert(start_ix <= i->ix && i->ix <= planters.size());

		for (auto p = planters.begin()+start_ix; p != planters.begin()+i->ix; ++p) {
			if (x_test < p->bcube.x1()) break; // planters are sorted by x1, none after this can match
			if (!p->bcube.contains_pt_xy_exp(pos, expand)) continue;
			// treat this as a tree rather than a planter by testing against a circle, since trees aren't otherwise included
			if (dist_xy_less_than(pos, p->pos, (p->radius + expand))) {color = DK_GREEN; return 1;}
		}
	} // for i
	start_ix = 0;

	if (!skip_in_road) { // fire hydrants are now placed on the edges of the road, so they're not inside plots and are skipped here
		for (auto i = fhydrant_groups.begin(); i != fhydrant_groups.end(); start_ix = i->ix, ++i) {
			if (!i->contains_pt_xy(pos)) continue;
			assert(start_ix <= i->ix && i->ix <= fhydrants.size());

			for (auto b = fhydrants.begin()+start_ix; b != fhydrants.begin()+i->ix; ++b) {
				if (pos.x < b->bcube.x1()) break; // fire_hydrants are sorted by x1, none after this can match
				if (dist_xy_less_than(pos, b->pos, b->radius)) {color = colorRGBA(1.0, 0.75, 0.0); return 1;} // orange/yellow color
			}
		} // for i
		start_ix = 0;
	}
	if (check_city_obj_bcube_pt_xy_contain(divider_groups, dividers, pos, obj_ix)) {
		assert(obj_ix < dividers.size());
		color = plot_divider_types[dividers[obj_ix].type].get_avg_color();
		return 1;
	}
	for (auto i = pool_groups.begin(); i != pool_groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy(pos)) continue;
		assert(start_ix <= i->ix && i->ix <= pools.size());

		for (auto b = pools.begin()+start_ix; b != pools.begin()+i->ix; ++b) {
			if (pos.x < b->bcube.x1()) break; // pools are sorted by x1, none after this can match
			if (!b->bcube.contains_pt_xy(pos)) continue;
			if (b->above_ground && !dist_xy_less_than(pos, point(b->bcube.xc(), b->bcube.yc(), b->bcube.z1()), b->get_radius())) continue; // circular in-ground pool
			color = b->wcolor; // return water color
			return 1;
		}
	} // for i
	start_ix = 0;
	if (check_city_obj_bcube_pt_xy_contain(sstation_groups, sstations, pos, obj_ix)) {color = colorRGBA(0.6, 0.8, 0.4, 1.0); return 1;} // light olive
	if (check_city_obj_bcube_pt_xy_contain(trashcan_groups, trashcans, pos, obj_ix)) {color = LT_BROWN; return 1;}
	// Note: ppoles are skipped for now
	return 0;
}

void city_obj_placer_t::get_occluders(pos_dir_up const &pdu, vect_cube_t &occluders) const {
	if (dividers.empty()) return; // dividers are currently the only occluders
	float const dmax(0.25f*(X_SCENE_SIZE + Y_SCENE_SIZE)); // set far clipping plane to 1/4 a tile (currently 2.0)
	unsigned start_ix(0);

	for (auto i = divider_groups.begin(); i != divider_groups.end(); start_ix = i->ix, ++i) {
		if (!dist_less_than(pdu.pos, i->closest_pt(pdu.pos), dmax) || !pdu.cube_visible(*i)) continue;
		assert(start_ix <= i->ix && i->ix <= dividers.size());

		for (auto d = dividers.begin()+start_ix; d != dividers.begin()+i->ix; ++d) {
			assert(d->type < DIV_NUM_TYPES);
			if (!plot_divider_types[d->type].is_occluder) continue; // skip
			if (d->bcube.z1() > pdu.pos.z || d->bcube.z2() < pdu.pos.z) continue; // z-range does not include the camera
			if (dist_less_than(pdu.pos, d->bcube.closest_pt(pdu.pos), dmax) && pdu.cube_visible(d->bcube)) {occluders.push_back(d->bcube);}
		}
	} // for i
}


