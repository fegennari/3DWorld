// 3D World - City Skyway and Moving Walkway
// by Frank Gennari
// 07/22/2024

#include "city_objects.h"
#include "lightmap.h" // for light_source

extern bool player_on_moving_ww, city_lights_custom_bcube;
extern int animate2, frame_counter;
extern float fticks;
extern vector<light_source> dl_sources;

int get_concrete_tid();
void set_tile_floor_texture();


class tile_drawer_t {
	uint64_t last_tile_id=0;
	bool tile_was_set=0;
public:
	void next_cube(cube_t const &c, draw_state_t &dstate, quad_batch_draw &qbd) {
		unsigned const tile_id(get_tile_id_containing_point_no_xyoff(c.get_cube_center()));
		if (tile_was_set && tile_id == last_tile_id) return;
		qbd.draw_and_clear();
		dstate.begin_tile(c.get_cube_center(), 1);
		last_tile_id = tile_id;
	}
	void end_draw(quad_batch_draw &qbd) {
		qbd.draw_and_clear();
		tile_was_set = 0;
	}
};

void draw_long_cube(cube_t const &c, colorRGBA const &color, draw_state_t &dstate, quad_batch_draw &qbd, tile_drawer_t &td, float dist_scale,
	bool shadow_only=0, bool skip_bottom=0, bool skip_top=0, float tscale=0.0, unsigned skip_dims=0, bool swap_tc_xy=0, float tex_tyoff=0.0)
{
	bool const dim(c.dx() < c.dy()); // longer dim; only supports splitting in one dim
	float const tile_sz(dim ? MESH_Y_SIZE*DY_VAL : MESH_X_SIZE*DX_VAL), length(c.get_sz_dim(dim));
	unsigned const num_segs(shadow_only ? 1U : unsigned(ceil(2.0*length/tile_sz)));
	float const step_len(length/num_segs);
	cube_t c2(c);

	for (unsigned s = 0; s < num_segs; ++s) { // split into segments, one per tile shadow map
		bool const beg(s == 0), end(s+1 == num_segs);
		c2.d[dim][1] = (end ? c.d[dim][1] : (c2.d[dim][0] + step_len));

		if (dstate.check_cube_visible(c2, dist_scale)) { // VFC
			unsigned skip_dims_seg(skip_dims);
			if (!(beg && dstate.camera_bs[dim] < c2.d[dim][0]) && !(end && dstate.camera_bs[dim] > c2.d[dim][1])) {skip_dims_seg |= (1 << unsigned(dim));} // skip int seg ends
			td.next_cube(c2, dstate, qbd);
			unsigned const verts_start(qbd.verts.size());
			dstate.draw_cube(qbd, c2, color, skip_bottom, tscale, skip_dims_seg, 0, 0, swap_tc_xy, 1.0, 1.0, 1.0, skip_top);

			if (tex_tyoff != 0.0) {
				tex_tyoff *= tscale;
				for (auto v = qbd.verts.begin()+verts_start; v != qbd.verts.end(); ++v) {v->t[1] += tex_tyoff;}
			}
		}
		c2.d[dim][0] += step_len;
	} // for s
}

// moving_walkway_t

moving_walkway_t::moving_walkway_t(cube_t const &c, bool dim_, bool dir_, float speed_) : cube_t(c), dim(dim_), dir(dir_), speed(speed_) {
	float const side_width(0.08*c.get_sz_dim(!dim));
	track = c;
	track.z2() = c.z1() + 0.05*c.dz();
	track.expand_in_dim(dim, -0.5*side_width); // slight shrink

	for (unsigned d = 0; d < 2; ++d) {
		sides[d] = c;
		sides[d].d[!dim][!d] = track.d[!dim][d] = c.d[!dim][d] + (d ? -1.0 : 1.0)*side_width;
	}
}
void moving_walkway_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, tile_drawer_t &td, bool shadow_only, bool draw_track, bool draw_sides) const {
	float const dist_scale = 0.5;

	if (draw_track) { // track
		float const tscale(1.0/track.get_sz_dim(!dim)); // scale to fit track width
		float const tex_tyoff(speed*move_time*(dir ? -1.0 : 1.0));
		draw_long_cube(track, WHITE, dstate, qbds.qbd, td, dist_scale, shadow_only, 1, 0, tscale, 0, dim, tex_tyoff); // skip_bottom=1, skip_top=0
		
		if (active && animate2 && !shadow_only) {
			move_time += fticks;
			if (move_time > 600.0*TICKS_PER_SECOND) {move_time = 0.0;} // reset every 10 min. to avoid precision problems
		}
	}
	if (draw_sides) { // sides and railings
		float const height(dz()), railing_thick(0.04*height), side_width(sides[0].get_sz_dim(!dim)), railing_shrink(0.15*side_width);

		for (unsigned d = 0; d < 2; ++d) { // sides + railings
			cube_t railing(sides[d]);
			railing.z1() += 0.8*height;
			railing.z2() += railing_thick;
			railing.expand_in_dim( dim,  railing_thick );
			railing.expand_in_dim(!dim, -railing_shrink);
			draw_long_cube(railing,  BKGRAY,  dstate, qbds.untex_qbd, td, 0.3*dist_scale, shadow_only, 0); // near black, skip_bottom=0
			draw_long_cube(sides[d], LT_GRAY, dstate, qbds.untex_qbd, td, 0.5*dist_scale, shadow_only, 1); // skip_bottom=1
		} // for d
	}
}
bool moving_walkway_t::proc_sphere_coll(point &pos, point const &p_last, float radius, point const &xlate, vector3d *cnorm) const {
	cube_t const bc(*this + xlate);
	if (!bc.contains_pt_xy_exp(pos, radius)) return 0; // not on/near walkway
	float const zval(max(pos.z, p_last.z));
	if (zval + radius < bc.z1()) return 0; // below walkway
	cube_t const track_xf(track + xlate);
	bool ret(0);

	if (track_xf.contains_pt_xy(pos)) { // on track
		pos.z = track_xf.z2() + radius;
	
		if (active && animate2) {
			if (frame_counter > last_update_frame) { // move at most once per frame
				pos[dim] += (dir ? 1.0 : -1.0)*fticks*speed;
				last_update_frame = frame_counter;
			}
			player_on_moving_ww = 1;
		}
		ret = 1;
	}
	for (unsigned d = 0; d < 2; ++d) {ret |= sphere_cube_int_update_pos(pos, radius, (sides[d] + xlate), p_last, 0, cnorm);} // check sides
	return ret;
}

// skyway_t

void skyway_t::init(cube_t const &c, bool dim_) {
	valid = 1;
	dim   = dim_;
	bcube = c;
	set_bsphere_from_bcube();
	float const height(bcube.dz()), width(bcube.get_sz_dim(!dim)), wall_width(0.05*width);
	cube_t center(bcube);
	center.expand_in_dim(!dim, -wall_width); // subtract off side walls
	bot = bcube;
	top = center;
	top.expand_in_dim(dim, -wall_width); // subtract off end walls
	cube_t nonbot(bcube);
	bot.z2() = nonbot.z1() = bcube.z1() + 0.085*height;
	top.z1() = bcube.z2() - 0.05*height;
	entrances.reserve(ww_conns.size());

	for (unsigned d = 0; d < 2; ++d) { // add sides, which will be cut by entrances
		cube_t side(nonbot);
		side.d[!dim][!d] = bcube.d[!dim][d] + (d ? -1.0 : 1.0)*wall_width;
		swap_cube_dims(side, !dim, 2); // swap so that subtract can be done in the XY plane
		sides.push_back(side);
	}
	for (skyway_conn_t const &conn : ww_conns) {
		cube_t entrance(conn);
		entrance.d[!dim][conn.dir] = center.d[!dim][!conn.dir];
		entrance.expand_in_dim( dim, -0.1*conn.get_sz_dim(dim)); // shrink sides
		entrance.expand_in_dim(2,    -0.55*FLOOR_THICK_VAL_OFFICE*conn.dz()); // lower the ceiling and raise the floor slightly
		entrances.push_back(entrance);
		cube_t sub(entrance);
		swap_cube_dims(sub, !dim, 2);
		subtract_cube_from_cubes(sub, sides);
		float const max_step_height(0.08*conn.dz()); // should be about the same for all conns
		float const dz(entrance.z1() - bot.z2());
		if (dz < max_step_height) continue; // no steps needed
		unsigned const num_steps(ceil(dz/max_step_height));
		bool const dir(!conn.dir); // dir relative to skyway, not building
		float const step_height(dz/(num_steps+1)), step_len((dir ? -1.0 : 1.0)*1.2*step_height), wall_pos(center.d[!dim][dir]);
		cube_t step(entrance);
		set_cube_zvals(step, bot.z2(), entrance.z1());
		step.d[!dim][ dir] = wall_pos; // starts at inner wall
		step.d[!dim][!dir] = wall_pos + step_len;

		for (unsigned n = 0; n < num_steps; ++n) {
			step.z2() -= step_height;
			assert(step.is_strictly_normalized());
			steps.push_back(step);
			step.translate_dim(!dim, step_len);
		}
	} // for conn
	// cut windows into long sides
	float const window_z1(center.z1() + 0.2*height), window_z2(window_z1 + 0.25*height);
	float const window_spacing(0.8*width), window_hwidth(0.3*window_spacing);
	vect_cube_t cut_sides, side_parts;

	for (cube_t &c : sides) {
		if (c.d[!dim][0] > window_z1 || c.d[!dim][1] < window_z2) {cut_sides.push_back(c); continue;} // not tall enough to split; note that dim is still swapped
		float const len(c.get_sz_dim(dim));
		if (len <= 2.0*window_spacing) {cut_sides.push_back(c); continue;} // not long enough to split
		unsigned const num_splits(len/window_spacing);
		float const split_len(len/num_splits);
		cube_t window(c); // !dim and Z are swapped
		window.d[!dim][0] = window_z1;
		window.d[!dim][1] = window_z2;
		side_parts.clear();
		side_parts.push_back(c);

		for (unsigned n = 1; n < num_splits; ++n) { // add a window between each split point
			set_wall_width(window, (c.d[dim][0] + n*split_len), window_hwidth, dim);
			subtract_cube_from_cubes(window, side_parts);
		}
		vector_add_to(side_parts, cut_sides);
	} // for c
	sides.swap(cut_sides);
	for (cube_t &c : sides) {swap_cube_dims(c, !dim, 2);} // swap dims back

	for (unsigned d = 0; d < 2; ++d) { // add ends
		cube_t end(center);
		end.z1() = nonbot.z1();
		end.d[dim][!d] = bcube.d[dim][d] + (d ? -1.0 : 1.0)*wall_width;
		sides.push_back(end);
	}
	sort(sides.begin(), sides.end(), cmp_by_tile()); // optimization for drawing
	// add moving walkways
	float const ww_hwidth(0.11*width), ww_height(0.8*ww_hwidth), ww_end_gap(0.75*width), centerline(bot.get_center_dim(!dim));
	float const speed = 0.005;
	cube_t ww_area(bot);
	set_cube_zvals(ww_area, bot.z2(), bot.z2()+ww_height);
	ww_area.expand_in_dim(dim, -ww_end_gap);
	// remove sections around walkway entrances
	vect_cube_t mww_segs;
	mww_segs.push_back(ww_area);

	for (skyway_conn_t const &conn : ww_conns) {
		cube_t sub(bcube);
		float const mww_gap(conn.building ? 2.0*conn.building->get_doorway_width() : 0.75*conn.get_sz_dim(dim));
		set_wall_width(sub, conn.get_center_dim(dim), 0.5*mww_gap, dim);
		subtract_cube_from_cubes(sub, mww_segs);
	}
	for (cube_t ww : mww_segs) {
		if (ww.get_sz_dim(dim) < 0.5*width) continue; // too short

		for (unsigned d = 0; d < 2; ++d) {
			set_wall_width(ww, (centerline + (d ? -1.0 : 1.0)*1.02*ww_hwidth), ww_hwidth, !dim);
			mwws.emplace_back(ww, dim, (bool(d) ^ dim), speed);
		}
	}
	// divide roof into panels and add lights
	float const length(top.get_sz_dim(dim)), max_panel_len(1.5*width);
	num_roof_panels = ceil(length/max_panel_len);
	float const light_dz(0.02*bcube.dz()), light_zval(top.z1() - light_dz), panel_len(length/num_roof_panels);

	for (unsigned n = 0; n < num_roof_panels; ++n) {
		point lpos(0.0, 0.0, light_zval);
		lpos[!dim] = top.get_center_dim(!dim);
		lpos[ dim] = top.d[dim][0] + (n + 0.5)*panel_len;
		lights.emplace_back(lpos);
	}
}

void skyway_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, bool shadow_only) const {
	float const dist_scale = 0.7;
	if (!valid || !dstate.check_cube_visible(bcube, dist_scale)) return; // VFC/distance culling
	//highres_timer_t timer("Draw Skyway"); // 0.03ms below, 0.08ms above, 0.09ms inside, 0.02ms at night (n shadows)
	tile_drawer_t td;
	if (!shadow_only) {select_texture(get_concrete_tid());}
	bind_default_flat_normal_map();
	colorRGBA const ext_color(LT_GRAY);
	bool const player_above_floor(dstate.camera_bs.z > bcube.z1());
	float const centerline(bcube.get_center_dim(!dim));
	float const side_shift(shadow_only ? 0.1*(top.x1() - bcube.x1()) : 0.0); // shrink by 10% of wall width in the shadow pass
	float const tscale = 16.0;
	
	for (cube_t c : sides) {
		bool const skip_bottom(c.z1() == bot.z2()), is_end(c.d[!dim][0] == top.d[!dim][0] && c.d[!dim][1] == top.d[!dim][1]);
		unsigned skip_dims(0);
		if (is_end) {skip_dims|= (1 << unsigned(!dim));} // don't need to draw sides of ends
		if (shadow_only) {c.expand_in_dim(!dim, (is_end ? 1.0 : -1.0)*side_shift);} // shrink walls slightly to prevent shadow acne
		draw_long_cube(c, ext_color, dstate, qbds.qbd, td, dist_scale, shadow_only, skip_bottom, 0, tscale, skip_dims);
	}
	if (player_above_floor) { // draw steps if player is above the floor
		for (cube_t const &c : steps) {
			if (!dstate.check_cube_visible(c, 0.25*dist_scale)) continue;
			td.next_cube(c, dstate, qbds.qbd);
			dstate.draw_cube(qbds.qbd, c, WHITE, 1, tscale); // skip_bottom=1
		}
	}
	if (shadow_only) {
		draw_long_cube(bot, ext_color, dstate, qbds.qbd, td, dist_scale, shadow_only, 0, 0, tscale); // draw all sides
	}
	else {
		draw_long_cube(bot, ext_color, dstate, qbds.qbd, td, dist_scale, shadow_only, 0, 1, tscale); // draw all sides except for top
		td.end_draw(qbds.qbd);
		set_tile_floor_texture();
		draw_long_cube(bot, GRAY, dstate, qbds.qbd, td, dist_scale, shadow_only, 1, 0, 4.0, 3); // draw top only
		bind_default_flat_normal_map();
	}
	td.end_draw(qbds.qbd);
	if (!player_above_floor) return; // skip the remainder of the non-visible interior items
	// draw moving walkway sides as untextured
	dstate.set_untextured_material();
	for (moving_walkway_t const &mww : mwws) {mww.draw(dstate, qbds, td, shadow_only, 0, 1);} // draw_track=0, draw_sides=1
	td.end_draw(qbds.untex_qbd);
	
	if (!shadow_only) { // draw moving walkway moving textured tracks; not needed for the shadow pass
		select_texture(get_texture_by_name("interiors/walkway_track.png", 1, 0, 1, 8.0));
		for (moving_walkway_t const &mww : mwws) {mww.draw(dstate, qbds, td, shadow_only, 1, 0);} // draw_track=1, draw_sides=0
		qbds.qbd.draw_and_clear();
		select_texture(WHITE_TEX);
	}
	// draw ramps, entryway frames, lights, and roof
	if (!shadow_only) {
		for (skyway_conn_t const &conn : ww_conns) { // ramps
			bool const dir(conn.dir);
			float const conn_zfloor(conn.z1() + 0.5*FLOOR_THICK_VAL_OFFICE*conn.dz());
			if (bot.z2() <= conn_zfloor) continue; // sloping down rather than up when entering skyway
			cube_t ramp(conn);
			ramp.d[!dim][!dir] += (dir ? -1.0 : 1.0)*0.1*conn.dz(); // extend into walkway
			ramp.expand_in_dim(dim, -0.1*conn.get_sz_dim(dim)); // shrink sides, same as entrances
			set_cube_zvals(ramp, min(conn_zfloor, bot.z2()), max(conn_zfloor, bot.z2()));
			if (!dstate.check_cube_visible(ramp, 0.1*dist_scale)) continue;
			point pts[4];
			for (unsigned d = 0; d < 2; ++d) {pts[2*d][dim] = pts[2*d+1][dim] = ramp.d[dim][d];} // sides of ramp
			pts[0][!dim] = pts[3][!dim] = ramp.d[!dim][!dir];
			pts[1][!dim] = pts[2][!dim] = ramp.d[!dim][ dir];
			pts[0].z = pts[3].z = ramp.z1();
			pts[1].z = pts[2].z = ramp.z2();
			vector3d const normal(((dim ^ dir) ? 1.0 : -1.0)*cross_product(pts[0]-pts[1], pts[3]-pts[0]).get_norm());
			qbds.untex_qbd.add_quad_pts(pts, GRAY, normal);
			td.next_cube(ramp, dstate, qbds.untex_qbd); // untextured
		} // for conn
	}
	for (cube_t const &entrance : entrances) { // entrance frames
		float const ecl(entrance.get_center_dim(!dim)), ewidth(entrance.get_sz_dim(!dim)), wwidth(0.25*ewidth), thickness(0.05*ewidth);
		bool const dir(ecl > centerline);
		cube_t e(entrance);
		if (bot.z2() < e.z1()) {e.d[!dim][dir] = ecl;} // raised walkway with stairs: frame ends at center of wall to avoid drawing its exterior under the bottom of the walkway
		else {e.d[!dim][dir] -= (dir ? -1.0 : 1.0)*thickness;} // flush stairs: extend to cover inside of wall (on walkway side)
		e.d[!dim][!dir] += (dir ? -1.0 : 1.0)*thickness; // extend outward
		e.expand_in_dim( dim, wwidth   );
		e.z2() += wwidth;
		min_eq(e.z1(), bot.z2()); // extends to the floor if needed
		if (!dstate.check_cube_visible(e, 0.1*dist_scale)) continue;
		cube_t etop(e), esides(e);
		etop.z1()  = esides.z2() = entrance.z2() - thickness;
		td.next_cube(e, dstate, qbds.untex_qbd);
		dstate.draw_cube(qbds.untex_qbd, etop, WHITE); // skip_bottom=0
		if (shadow_only) continue; // skip sides in shadow pass; top is needed to cover gap between shrunk walls

		for (unsigned d = 0; d < 2; ++d) { // sides
			cube_t eside(esides);
			eside.d[dim][!d] = eside.d[dim][d] + (d ? -1.0 : 1.0)*(wwidth + 2.0*thickness);
			dstate.draw_cube(qbds.untex_qbd, eside, WHITE, 1, 0.0, 4); // skip_bottom=1, skip dim Z
		}
	} // for e
	td.end_draw(qbds.untex_qbd);
		
	if (!lights.empty()) { // add ceiling lights
		bool const is_on(are_lights_on());
		dstate.s.set_cur_color(WHITE);
		if (is_on) {dstate.s.set_color_e(WHITE);}

		for (roof_light_t const &lpos : lights) { // draw light fixtures
			if (shadow_only && lpos.x == camera_pdu.pos.x && lpos.y == camera_pdu.pos.y) continue; // no self-shadows (zval differs)
			float const light_height(top.z1() - lpos.z), light_radius(4.0*light_height); // for the light fixture itself
			assert(light_height > 0.0);
			cube_t light_bc;
			light_bc.set_from_sphere(lpos, light_radius);
			set_cube_zvals(light_bc, lpos.z, top.z1()+0.1*light_height); // top slightly inside the glass to preventZ-fighting
			if (!dstate.check_cube_visible(light_bc, 0.12*dist_scale)) continue;
			if (!is_on && !shadow_only) {dstate.begin_tile(lpos, 1);}
			unsigned const ndiv(shadow_only ? 16 : max(4U, min(32U, unsigned(12.0f/p2p_dist(dstate.camera_bs, lpos)))));
			draw_fast_cylinder(point(lpos.x, lpos.y, light_bc.z1()), point(lpos.x, lpos.y, light_bc.z2()), light_radius, light_radius, ndiv, 0, 1); // untextured, draw ends
		} // for lpos
		if (is_on) {dstate.s.clear_color_e();}
	}
	// draw window frames/dividers; untextured black, so grouping into tiles for shadows isn't needed and they can all be batched together
	float const frame_width(0.75*top.dz());
	cube_t top_exp(top);
	top_exp.expand_in_dim(2, 0.25*frame_width); // extends outside the glass
	assert(qbds.untex_qbd.empty());

	for (unsigned dir = 0; dir < 2; ++dir) { // sides
		cube_t frame(top_exp);
		frame.d[!dim][!dir] = frame.d[!dim][dir] + (dir ? -1.0 : 1.0)*frame_width;
		if (shadow_only) {frame.d[!dim][dir] -= (dir ? -1.0 : 1.0)*side_shift;} // shift to cover the gap between the side
		dstate.draw_cube(qbds.untex_qbd, frame, BLACK); // draw all sides
	}
	float const panel_len((top.get_sz_dim(dim) - frame_width)/num_roof_panels);
	cube_t frame(top_exp);
	frame.d[dim][1] = frame.d[dim][0] + frame_width;
	frame.expand_in_dim(!dim, -frame_width); // meets the edge frame

	for (unsigned n = 0; n <= num_roof_panels; ++n) {
		dstate.draw_cube(qbds.untex_qbd, frame, BLACK, 0, 0.0, (1 << unsigned(!dim))); // skip ends
		frame.translate_dim(dim, panel_len);
	}
	qbds.untex_qbd.draw_and_clear();
		
	if (!shadow_only) { // draw transparent top glass panel; Z only; drawn last
		enable_blend();
		glDepthMask(GL_FALSE); // disable depth writing so that terrain and grass are drawn over the glass
		draw_long_cube(top, colorRGBA(1.0, 1.0, 1.0, 0.25), dstate, qbds.untex_qbd, td, dist_scale, shadow_only, 0, 0, 0.0, 3);
		td.end_draw(qbds.untex_qbd);
		glDepthMask(GL_TRUE);
		disable_blend();
	}
	dstate.unset_untextured_material();
}

bool skyway_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	if (!valid) return 0;
	cube_t const bc_cs(bcube + xlate);
	if (!bc_cs.contains_pt_xy_exp(pos_, radius_)) return 0; // expand by radius so that the player can step up from a walkway
	float const zval(max(pos_.z, p_last.z)), top_z2(top.z2() + xlate.z), bot_z2(bot.z2() + xlate.z);
	bool ret(0);

	if (zval > top_z2) { // above the skyway - walk on the top glass
		if (bc_cs.contains_pt_xy_exp(pos_, 0.1*radius_)) { // expand slightly for a bit of overhang so that player doesn't fall inside a walkway
			max_eq(pos_.z, (top_z2 + radius_));
			return 1;
		}
		return 0;
	}
	if (zval > bot_z2) { // above the bottom
		max_eq(pos_.z, (bot_z2 + radius_));
		for (moving_walkway_t const &mww : mwws) {mww.proc_sphere_coll(pos_, p_last, radius_, xlate, cnorm);}
		ret = 1;
	}
	for (cube_t const &side : sides) {
		cube_t const c(side + xlate);
		if (!c.contains_pt_xy_exp(pos_, radius_)) continue;

		if (zval > c.z2()) {
			max_eq(pos_.z, (c.z2() + radius_));
			ret = 1;
		}
		else {ret |= sphere_cube_int_update_pos(pos_, radius_, c, p_last, 0, cnorm);}
	} // for side
	for (cube_t const &s : steps) {
		cube_t const c(s + xlate);
		float const step_up_z(c.z2() + radius_);

		if (step_up_z + 0.1*radius_ > zval) { // can step up onto stair (with a bit of hysteresis)
			if (pos_[dim] > c.d[dim][0] && pos_[dim] < c.d[dim][1]) { // within stairs range
				if (pos_[!dim] > c.d[!dim][0]-radius_ && pos_[!dim] < c.d[!dim][1]+radius_) { // intersecting stair
					max_eq(pos_.z, step_up_z); // step up onto stair
					ret = 1;
					continue;
				}
			}
		}
		ret |= sphere_cube_int_update_pos(pos_, radius_, c, p_last, 0, cnorm);
	} // for s
	if (ret) return 1;
	return sphere_cube_int_update_pos(pos_, radius_, bc_cs, p_last, 0, cnorm); // exterior coll
}

void skyway_t::get_building_signs(vector<sign_t> &signs) const {
	float const sign_height(0.06*bcube.dz()), sign_hwidth(4.0*sign_height), sign_depth(0.2*sign_height);

	for (skyway_conn_t const &conn : ww_conns) { // ramps
		if (!conn.building || conn.building->name.empty()) continue; // no name, no sign
		bool const dir(!conn.dir); // dir relative to skyway, not building
		float const wall_pos(top.d[!dim][dir]); // inner wall of skyway
		float const sign_z1(conn.z2() + 0.5*sign_height);
		cube_t sign(conn);
		set_cube_zvals(sign, sign_z1, sign_z1+sign_height);
		sign.d[!dim][ dir] = wall_pos;
		sign.d[!dim][!dir] = wall_pos + (dir ? -1.0 : 1.0)*sign_depth;
		set_wall_width(sign, conn.get_center_dim(dim), sign_hwidth, dim);
		signs.emplace_back(sign, !dim, !dir, conn.building->name, WHITE, BLACK, 0, 0, 1, 0, 0, 1); // two_sided=0, emissive=0, small=1, scrolling=0, fs=0, in_skyway=1
		signs.back().draw_zmin = bcube.z1(); // skip draw if camera is below the bottom of the skyway
	} // for conn
}

bool skyway_t::are_lights_on() const {
	bool const always_on = 0;
	return (always_on || is_night(0.15)); // adjust to turn on earlier
}
void skyway_t::add_lights(vector3d const &xlate, cube_t &lights_bcube) const {
	if (lights.empty() || !are_lights_on())        return;
	if (get_camera_pos().z < bcube.z1() + xlate.z) return; // lights not visible from under skyway
	float const ldist(4.0*bcube.get_sz_dim(!dim)), light_dz(0.5*bcube.dz());
	colorRGBA const light_color(1.0, 0.95, 0.9, 1.0);

	for (roof_light_t const &light : lights) {
		if (!lights_bcube.contains_pt_xy(light)) continue; // not contained within the light volume
		point const lpos(light.x, light.y, light.z+light_dz); // make higher to produce better shadows
		if (!camera_pdu.sphere_visible_test((lpos + xlate), ldist)) continue; // VFC
		min_eq(lights_bcube.z1(), lpos.z-ldist);
		max_eq(lights_bcube.z2(), lpos.z); // pointed down
		dl_sources.emplace_back(ldist, lpos, lpos, light_color, 0, -plus_z, 0.35); // points down
		cube_t clip_cube(bcube);
		clip_cube.z2() = lpos.z; // extend to cover the top glass surface
		dl_sources.back().set_custom_bcube(clip_cube); // don't light surrounding buildings or walkways
		if (light.cached_smap) {dl_sources.back().assign_smap_id(uintptr_t(&light)/sizeof(void *));} // cache shadow map on second frame
		light.cached_smap = 1;
		city_lights_custom_bcube = 1;
	} // for lpos
}

