// 3D World - City Object Placement Header
// by Frank Gennari
// 08/14/21

#include "city_objects.h"
#include "tree_3dw.h" // for tree_placer_t

extern unsigned max_unique_trees;
extern tree_placer_t tree_placer;
extern city_params_t city_params;
extern object_model_loader_t building_obj_model_loader;

struct plot_divider_type_t {
	bool is_occluder;
	int tid;
	float wscale, hscale; // width and height scales
	colorRGBA color, map_color;
	string tex_name;

	plot_divider_type_t(string const &tn, float ws, float hs, bool ic, colorRGBA const &c, colorRGBA const &mc) :
		is_occluder(ic), tid(-1), wscale(ws), hscale(hs), color(c), map_color(mc), tex_name(tn) {}
	colorRGBA get_avg_color() const {return ((tid >= 0) ? texture_color(tid) : map_color).modulate_with(color);}
};
enum {DIV_WALL=0, DIV_FENCE, DIV_HEDGE, DIV_NUM_TYPES}; // types of plot dividers, with end terminator

plot_divider_type_t plot_divider_types[DIV_NUM_TYPES] = {
	plot_divider_type_t("cblock2.jpg", 0.50, 2.5, 1, WHITE, GRAY    ), // wall
	plot_divider_type_t("fence.jpg",   0.15, 2.0, 1, WHITE, LT_BROWN), // fence
	plot_divider_type_t("hedges.jpg",  1.00, 1.6, 0, GRAY,  GREEN   )  // hedge - too short to be an occluder
};

void add_house_driveways_for_plot(cube_t const &plot, vect_cube_t &driveways);
float get_sidewalk_width();
float get_inner_sidewalk_width();


bool city_obj_t::proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const {
	return sphere_cube_int_update_pos(pos_, radius_, (bcube + xlate), p_last, 1, 0, cnorm);
}

void bench_t::calc_bcube() {
	bcube.set_from_point(pos);
	bcube.expand_by(vector3d((dim ? 0.32 : 1.0), (dim ? 1.0 : 0.32), 0.0)*radius);
	bcube.z2() += 0.85*radius; // set bench height
}
/*static*/ void bench_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture(FENCE_TEX);} // normal map?
}
void bench_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const {
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

tree_planter_t::tree_planter_t(point const &pos_, float radius_, float height) : city_obj_t(pos_, radius_) {
	bcube.set_from_point(pos);
	bcube.expand_by_xy(radius);
	bcube.z2() += height;
}
/*static*/ void tree_planter_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture((dstate.pass_ix == 0) ? (int)DIRT_TEX : get_texture_by_name("roads/sidewalk.jpg"));}
}
void tree_planter_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const {
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

fire_hydrant_t::fire_hydrant_t(point const &pos_, float radius_, float height, vector3d const &orient_) : city_obj_t(pos_, radius_), cylin_radius(radius), orient(orient_) {
	bcube.set_from_sphere(*this);
	set_cube_zvals(bcube, pos.z, pos.z+height);
	pos.z += 0.5*height; // pos is bottom center point, make it the center
	max_eq(radius, 0.5f*height); // use a more accurate bounding sphere; Note: no cube root of (r*r + r*r + h*h)
}
/*static*/ void fire_hydrant_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {select_texture(WHITE_TEX);}
	if (!shadow_only) {dstate.s.set_cur_color(colorRGBA(1.0, 0.75, 0.0));}
}
/*static*/ void fire_hydrant_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {dstate.s.set_cur_color(WHITE);} // restore to default color
	city_obj_t::post_draw(dstate, shadow_only);
}
void fire_hydrant_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const { // Note: qbd is unused
	if (!dstate.check_cube_visible(bcube, dist_scale, shadow_only)) return;

	if (!shadow_only && building_obj_model_loader.is_model_valid(OBJ_MODEL_FHYDRANT)) {
		building_obj_model_loader.draw_model(dstate.s, pos, bcube, orient, WHITE, dstate.xlate, OBJ_MODEL_FHYDRANT, shadow_only);
	}
	else { // draw as a simple cylinder, untextured, top end only
		draw_fast_cylinder(point(pos.x, pos.y, bcube.z1()), point(pos.x, pos.y, bcube.z2()), 0.8*cylin_radius, 0.8*cylin_radius, (shadow_only ? 12 : N_CYL_SIDES), 0, 4);
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

/*static*/ void divider_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (shadow_only) return; // not textured
	assert(dstate.pass_ix < DIV_NUM_TYPES);
	plot_divider_type_t &pdt(plot_divider_types[dstate.pass_ix]);
	if (pdt.tid < 0) {pdt.tid = get_texture_by_name(pdt.tex_name);} // load/lookup texture if needed
	select_texture(pdt.tid);
}
void divider_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const {
	if (type != dstate.pass_ix) return; // this type not enabled in this pass
	if (!dstate.check_cube_visible(bcube, dist_scale, shadow_only)) return;
	assert(dstate.pass_ix < DIV_NUM_TYPES);
	dstate.draw_cube(qbd, bcube, color_wrapper(plot_divider_types[type].color), 1, 1.0/bcube.dz(), skip_dims); // skip bottom, scale texture to match the height
}

// passes: 0=in-ground walls, 1=in-ground water, 2=above ground sides, 3=above ground water
/*static*/ void swimming_pool_t::pre_draw(draw_state_t &dstate, bool shadow_only) {
	if (!shadow_only) {
		if      (dstate.pass_ix == 2) {select_texture(WHITE_TEX);} // sides/untextured
		else if (dstate.pass_ix == 0) {select_texture(get_texture_by_name("bathroom_tile.jpg"));} // walls
		else if (dstate.pass_ix == 1 || dstate.pass_ix == 3) {select_texture(get_texture_by_name("snow2.jpg"));} // water surface
		else {assert(0);}
	}
}
void swimming_pool_t::draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const {
	if ((dstate.pass_ix > 1) ^ above_ground) return; // not drawn in this pass
	if (!dstate.check_cube_visible(bcube, dist_scale, shadow_only)) return;

	if (above_ground) { // cylindrical; bcube should be square in XY
		point const camera_bs(camera_pdu.pos - dstate.xlate);
		float const radius(get_radius()), xc(bcube.xc()), yc(bcube.yc());
		unsigned const ndiv(shadow_only ? 24 : max(4U, min(32U, unsigned(4.0f*dist_scale*get_draw_tile_dist()/p2p_dist(camera_bs, bcube.closest_pt(camera_bs))))));

		if (dstate.pass_ix == 2) { // draw sides
			dstate.s.set_cur_color(color);
			draw_fast_cylinder(point(xc, yc, bcube.z1()), point(xc, yc, bcube.z2()), radius, radius, ndiv, 0, 0); // untextured, no ends
		}
		else if (dstate.pass_ix == 3) { // draw water surface
			dstate.s.set_cur_color(wcolor);
			draw_circle_normal(0.0, radius, ndiv, 0, point(xc, yc, (bcube.z2() - 0.1*bcube.dz()))); // shift slightly below the top
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
		float const radius(get_radius()), xc(bcube.xc() + xlate.x), yc(bcube.yc() + xlate.y), z1(bcube.z1() + xlate.z), z2(bcube.z2() + xlate.z);
		return sphere_vert_cylin_intersect(pos_, radius_, cylinder_3dw(point(xc, yc, z1), point(xc, yc, z2), radius, radius), cnorm); // checks sides
	}
	else {return city_obj_t::proc_sphere_coll(pos_, p_last, radius_, xlate, cnorm);}
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

// non-const because this sets driveway_t::car_ix
void city_obj_placer_t::add_cars_to_driveways(vector<car_t> &cars, vector<road_plot_t> const &plots, vector<vect_cube_t> &plot_colliders, unsigned city_id, rand_gen_t &rgen) {
	car_t car;
	car.park();
	car.cur_city = city_id;
	car.cur_road_type = TYPE_DRIVEWAY;

	for (auto i = driveways.begin(); i != driveways.end(); ++i) {
		if (rgen.rand_float() < 0.5) continue; // no car in this driveway 50% of the time
		car.cur_road = (unsigned short)i->plot_ix; // store plot_ix in road field
		car.cur_seg  = (unsigned short)(i - driveways.begin()); // store driveway index in cur_seg
		vector3d car_sz(city_params.get_nom_car_size()); // {length, width, height}
		cube_t const &plot(plots[i->plot_ix]);
		car.dim    = (i->y1() == plot.y1() || i->y2() == plot.y2()); // check which edge of the plot the driveway is connected to, which is more accurate than the aspect ratio
		if (i->get_sz_dim(car.dim) < 1.6*car_sz.x || i->get_sz_dim(!car.dim) < 1.25*car_sz.y) continue; // driveway is too small to fit this car
		car.dir    = rgen.rand_bool(); // randomly pulled in vs. backed in, since we don't know the direction to the house anyway
		car.height = car_sz.z;
		float const pad_l(0.75*car_sz.x), pad_w(0.6*car_sz.y); // needs to be a bit larger to fit trucks
		if (car.dim) {swap(car_sz.x, car_sz.y);}
		point cpos(0.0, 0.0, (i->z2() + 0.5*car_sz.z));
		cpos[ car.dim] = rgen.rand_uniform(i->d[ car.dim][0]+pad_l, i->d[ car.dim][1]-pad_l);
		cpos[!car.dim] = rgen.rand_uniform(i->d[!car.dim][0]+pad_w, i->d[!car.dim][1]-pad_w); // not quite centered
		car.bcube.set_from_point(cpos);
		car.bcube.expand_by(0.5*car_sz);
		// check if this car intersects another parked car; this can only happen if two driveways intersect, which should be rare
		bool intersects(0);

		for (auto c = cars.rbegin(); c != cars.rend(); ++c) {
			if (c->cur_road != i->plot_ix) break; // prev plot, done
			if (car.bcube.intersects(c->bcube)) {intersects = 1; break;}
		}
		if (intersects) continue; // skip
		i->add_car(cars.size());
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
template<typename T> void city_obj_groups_t::create_groups(vector<T> &objs) {
	vector<T> new_objs;
	new_objs.reserve(objs.size());
	reserve(by_tile.size()); // the number of actual groups

	for (auto g = by_tile.begin(); g != by_tile.end(); ++g) {
		unsigned const group_start(new_objs.size());
		cube_with_ix_t group;

		for (auto i = g->second.begin(); i != g->second.end(); ++i) {
			assert(*i < objs.size());
			group.assign_or_union_with_cube(objs[*i].bcube);
			new_objs.push_back(objs[*i]);
		}
		sort(new_objs.begin()+group_start, new_objs.end());
		group.ix = new_objs.size();
		push_back(group);
	} // for g
	objs.swap(new_objs);
	by_tile.clear(); // no longer needed
}

// Note: blockers are used for placement of objects within this plot; colliders are used for pedestrian AI
void city_obj_placer_t::place_detail_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders,
	vector<point> const &tree_pos, rand_gen_t &rgen, bool is_residential)
{
	float const car_length(city_params.get_nom_car_size().x); // used as a size reference for other objects

	// place fire_hydrants; don't add fire hydrants in parks
	if (!plot.is_park) {
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
				fire_hydrant_groups.add_obj(fire_hydrant, fire_hydrants);
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
}

void city_obj_placer_t::place_plot_dividers(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, rand_gen_t &rgen) {
	assert(plot_subdiv_sz > 0.0);
	sub_plots.clear();
	if (plot.is_park) return; // no dividers in parks
	subdivide_plot_for_residential(plot, plot_subdiv_sz, 0, sub_plots); // parent_plot_ix=0, not needed
	if (sub_plots.size() <= 1) return; // nothing to divide
	if (rgen.rand_bool()) {std::reverse(sub_plots.begin(), sub_plots.end());} // reverse half the time so that we don't prefer a divider in one side or the other
	unsigned const shrink_dim(rgen.rand_bool()); // mostly arbitrary, could maybe even make this a constant 0
	float const sz_scale(0.06*city_params.road_width);
	unsigned const dividers_start(dividers.size());

	for (auto i = sub_plots.begin(); i != sub_plots.end(); ++i) {
		unsigned const type(rgen.rand()%(DIV_NUM_TYPES + 1)); // use a consistent divider type for all sides of this plot
		if (type >= DIV_NUM_TYPES) continue; // no divider for this plot
		// should we remove or move houses fences for divided sub-plots? I'm not sure how that would actually be possible at this point; or maybe skip dividers if the house has a fence?
		plot_divider_type_t const &pdt(plot_divider_types[type]);
		float const hwidth(0.5*sz_scale*pdt.wscale), z2(i->z1() + sz_scale*pdt.hscale);
		float const shrink_border(1.5*get_inner_sidewalk_width()); // needed for pedestrians to move along the edge of the plot; slightly larger to prevent collisions
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
					c.translate_dim(dim, (dir ? -1.0 : 1.0)*hwidth); // move inside the plot so that edges line up
					// clip to the sides to remove overlap; may not line up with a neibhgoring divider of a different type/width, but hopefully okay
					for (unsigned d = 0; d < 2; ++d) {
						if (c.d[!dim][d] != plot.d[!dim][d]) {c.d[!dim][d] -= (d ? 1.0 : -1.0)*hwidth;}
					}
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
				divider_t divider(c, type, dim, dir, skip_dims);
				divider_groups.add_obj(divider, dividers);
				colliders.push_back(divider.bcube);
				blockers .push_back(divider.bcube);
			} // for dir
		} // for dim
	} // for i
}

void city_obj_placer_t::place_residential_plot_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, rand_gen_t &rgen) {
	// assumes place_plot_dividers() has been called first to populate sub_plots
	float const min_spacing_to_plot_edge(0.5*city_params.road_width);
	colorRGBA const side_colors[5] = {WHITE, WHITE, GRAY, LT_BROWN, LT_BLUE};

	for (auto i = sub_plots.begin(); i != sub_plots.end(); ++i) {
		if (!i->is_residential || i->is_park || i->street_dir == 0) continue; // not a residential plot along a road
		if (rgen.rand_bool()) continue; // only add pools 50% of the time
		bool const dim((i->street_dir-1)>>1), dir((i->street_dir-1)&1); // direction to the road
		bool const above_ground(rgen.rand_bool());
		cube_t pool_area(*i);
		pool_area.d[dim][dir] = pool_area.get_center_dim(dim); // limit the pool to the back yard
		float const dmin(min(pool_area.dx(), pool_area.dy())); // or should this be based on city_params.road_width?

		for (unsigned d = 0; d < 2; ++d) { // keep pools away from the edges of plots; applies to sub-plots on the corners
			max_eq(pool_area.d[d][0], plot.d[d][0]+min_spacing_to_plot_edge);
			min_eq(pool_area.d[d][1], plot.d[d][1]-min_spacing_to_plot_edge);
		}
		pool_area.expand_by_xy(-0.05*dmin); // small shrink to keep away from walls, fences, and hedges
		vector3d pool_sz;
		pool_sz.z = (above_ground ? rgen.rand_uniform(0.08, 0.12)*city_params.road_width : 0.01f*dmin);

		for (unsigned d = 0; d < 2; ++d) {
			pool_sz[d] = ((above_ground && d == 1) ? pool_sz[0] : rgen.rand_uniform(0.4, 0.7)*dmin); // above_ground_cylin pools have square bcubes
			pool_area.d[d][1] -= pool_sz[d]; // shrink so that pool_area is where (x1, x2) can be placed
		}
		if (!pool_area.is_normalized()) continue; // pool area is too small; this can only happen due to shrink at plot edges
		point pool_llc;
		pool_llc.z = i->z2();

		for (unsigned n = 0; n < 10; ++n) { // make some attempts to generate a valid pool location
			for (unsigned d = 0; d < 2; ++d) {pool_llc[d] = rgen.rand_uniform(pool_area.d[d][0], pool_area.d[d][1]);}
			cube_t pool(pool_llc, (pool_llc + pool_sz));
			if (has_bcube_int_xy(pool, blockers, 0.08*dmin)) continue; // intersects some other object
			float const grayscale(rgen.rand_uniform(0.7, 1.0));
			float const water_white_comp(rgen.rand_uniform(0.1, 0.3)), extra_green(rgen.rand_uniform(0.2, 0.5)), lightness(rgen.rand_uniform(0.5, 0.8));
			colorRGBA const color(above_ground ? side_colors[rgen.rand()%5]: colorRGBA(grayscale, grayscale, grayscale));
			colorRGBA const wcolor(lightness*water_white_comp, lightness*(water_white_comp + extra_green), lightness);
			pool_groups.add_obj(swimming_pool_t(pool, color, wcolor, above_ground), pools);
			pool.z2() += 0.1*city_params.road_width; // extend upward to make a better collider
			colliders.push_back(pool);
			blockers .push_back(pool);
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
	draw_state_t &dstate, float dist_scale, bool shadow_only, bool not_using_qbd)
{
	if (objs.empty()) return;
	T::pre_draw(dstate, shadow_only);
	unsigned start_ix(0);
	assert(qbd.empty());

	for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
		if (!dstate.check_cube_visible(*g, dist_scale, shadow_only)) continue; // VFC/distance culling for group
		if (not_using_qbd) {dstate.begin_tile(g->get_cube_center(), 1, 1);} // must setup shader and tile shadow map before drawing
		assert(start_ix <= g->ix && g->ix <= objs.size());

		for (unsigned i = start_ix; i < g->ix; ++i) {
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

void city_obj_placer_t::clear() {
	parking_lots.clear(); benches.clear(); planters.clear(); fire_hydrants.clear(); driveways.clear(); dividers.clear(); pools.clear();
	bench_groups.clear(); planter_groups.clear(); fire_hydrant_groups.clear(); divider_groups.clear(); pool_groups.clear();
	num_spaces = filled_spaces = 0;
}

void city_obj_placer_t::gen_parking_and_place_objects(vector<road_plot_t> &plots, vector<vect_cube_t> &plot_colliders, vector<car_t> &cars,
	unsigned city_id, bool have_cars, bool is_residential)
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
		if (city_params.assign_house_plots && plot_subdiv_sz > 0.0) {
			place_plot_dividers(*i, bcubes, colliders, detail_rgen); // before placing trees
			place_residential_plot_objects(*i, bcubes, colliders, detail_rgen);
		}
		place_trees_in_plot (*i, bcubes, colliders, tree_pos, detail_rgen, buildings_end);
		place_detail_objects(*i, bcubes, colliders, tree_pos, detail_rgen, is_residential);
	} // for i
	if (have_cars) {add_cars_to_driveways(cars, plots, plot_colliders, city_id, rgen);}
	for (auto i = plot_colliders.begin(); i != plot_colliders.end(); ++i) {sort(i->begin(), i->end(), cube_by_x1());}
	bench_groups.create_groups(benches);
	planter_groups.create_groups(planters);
	fire_hydrant_groups.create_groups(fire_hydrants);
	divider_groups.create_groups(dividers);
	pool_groups.create_groups(pools);

	if (0) { // debug info printing
		cout << TXT(benches.size()) << TXT(bench_groups.size()) << TXT(planters.size()) << TXT(planter_groups.size()) << TXT(fire_hydrants.size())
			 << TXT(fire_hydrant_groups.size()) << TXT(dividers.size()) << TXT(divider_groups.size()) << TXT(pools.size()) << TXT(pool_groups.size()) << endl;
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

void city_obj_placer_t::draw_detail_objects(draw_state_t &dstate, bool shadow_only) {
	draw_objects(benches,       bench_groups,        dstate, 0.16, shadow_only, 0); // dist_scale=0.16
	draw_objects(fire_hydrants, fire_hydrant_groups, dstate, 0.07, shadow_only, 1); // dist_scale=0.12, not_using_qbd=1
			
	if (!shadow_only) { // low profile, not drawn in shadow pass
		for (dstate.pass_ix = 0; dstate.pass_ix < 2; ++dstate.pass_ix) { // {dirt, stone}
			draw_objects(planters, planter_groups, dstate, 0.1, shadow_only, 0); // dist_scale=0.1
		}
	}
	for (dstate.pass_ix = 0; dstate.pass_ix < 4; ++dstate.pass_ix) { // {in-ground walls, in-ground water, above ground sides, above ground water}
		if (shadow_only && dstate.pass_ix <= 1) continue; // only above ground pools are drawn in the shadow pass; water surface is drawn to prevent light leaks, but maybe should extend z1 lower
		float const dist_scales[4] = {0.1, 0.5, 0.3, 0.5};
		draw_objects(pools, pool_groups, dstate, dist_scales[dstate.pass_ix], shadow_only, (dstate.pass_ix > 1)); // final 2 passes don't use qbd
	}
	// Note: not the most efficient solution, as it required processing blocks and binding shadow maps multiple times
	for (dstate.pass_ix = 0; dstate.pass_ix < DIV_NUM_TYPES; ++dstate.pass_ix) { // {wall, fence, hedge}
		draw_objects(dividers, divider_groups, dstate, 0.2, shadow_only, 0); // dist_scale=0.2
	}
	dstate.pass_ix = 0; // reset back to 0
}

template<typename T> bool proc_vector_sphere_coll(vector<T> const &objs, point &pos, point const &p_last, float radius, vector3d const &xlate, vector3d *cnorm) {
	for (auto i = objs.begin(); i != objs.end(); ++i) { // Note: could use bench_groups
		if (i->proc_sphere_coll(pos, p_last, radius, xlate, cnorm)) return 1;
	}
	return 0;
}
bool city_obj_placer_t::proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm) const {
	vector3d const xlate(get_camera_coord_space_xlate());
	if (proc_vector_sphere_coll(benches,       pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(fire_hydrants, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(dividers,      pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(pools,         pos, p_last, radius, xlate, cnorm)) return 1;
	// Note: no coll with tree_planters because the tree coll should take care of it
	return 0;
}

template<typename T> void check_vector_line_intersect(vector<T> const &objs, point const &p1, point const &p2, float &t, bool &ret) {
	for (auto i = objs.begin(); i != objs.end(); ++i) {ret |= check_line_clip_update_t(p1, p2, t, i->bcube);} // check bounding cube
}
bool city_obj_placer_t::line_intersect(point const &p1, point const &p2, float &t) const { // Note: nothing to do for parking lots or tree_planters
	bool ret(0);
	check_vector_line_intersect(benches,       p1, p2, t, ret);
	check_vector_line_intersect(fire_hydrants, p1, p2, t, ret); // check bounding cube; cylinder intersection may be more accurate, but likely doesn't matter much
	check_vector_line_intersect(dividers,      p1, p2, t, ret);
	check_vector_line_intersect(pools,         p1, p2, t, ret);
	return ret;
}

bool city_obj_placer_t::get_color_at_xy(point const &pos, colorRGBA &color, bool skip_in_road) const {
	unsigned start_ix(0);

	for (auto i = bench_groups.begin(); i != bench_groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy(pos)) continue;
		assert(start_ix <= i->ix && i->ix <= benches.size());
					
		for (auto b = benches.begin()+start_ix; b != benches.begin()+i->ix; ++b) {
			if (pos.x < b->bcube.x1()) break; // benches are sorted by x1, no bench after this can match
			if (b->bcube.contains_pt_xy(pos)) {color = texture_color(FENCE_TEX); return 1;}
		}
	} // for i
	float const expand(0.15*city_params.road_width), x_test(pos.x + expand); // expand to approx tree diameter
	start_ix = 0;

	for (auto i = planter_groups.begin(); i != planter_groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy_exp(pos, expand)) continue;
		assert(start_ix <= i->ix && i->ix <= planters.size());

		for (auto p = planters.begin()+start_ix; p != planters.begin()+i->ix; ++p) {
			if (x_test < p->bcube.x1()) break; // planters are sorted by x1, no planter after this can match
			if (!p->bcube.contains_pt_xy_exp(pos, expand)) continue;
			// treat this as a tree rather than a planter by testing against a circle, since trees aren't otherwise included
			if (dist_xy_less_than(pos, p->pos, (p->radius + expand))) {color = DK_GREEN; return 1;}
		}
	} // for i
	start_ix = 0;

	if (!skip_in_road) { // fire hydrants are now placed on the edges of the road, so they're not inside plots and are skipped here
		for (auto i = fire_hydrant_groups.begin(); i != fire_hydrant_groups.end(); start_ix = i->ix, ++i) {
			if (!i->contains_pt_xy(pos)) continue;
			assert(start_ix <= i->ix && i->ix <= fire_hydrants.size());

			for (auto b = fire_hydrants.begin()+start_ix; b != fire_hydrants.begin()+i->ix; ++b) {
				if (pos.x < b->bcube.x1()) break; // fire_hydrants are sorted by x1, no fire_hydrant after this can match
				if (dist_xy_less_than(pos, b->pos, b->radius)) {color = colorRGBA(1.0, 0.75, 0.0); return 1;} // orange/yellow color
			}
		} // for i
		start_ix = 0;
	}
	for (auto i = divider_groups.begin(); i != divider_groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy(pos)) continue;
		assert(start_ix <= i->ix && i->ix <= dividers.size());

		for (auto b = dividers.begin()+start_ix; b != dividers.begin()+i->ix; ++b) {
			if (pos.x < b->bcube.x1()) break; // dividers are sorted by x1, no divider after this can match
			
			if (b->bcube.contains_pt_xy(pos)) {
				assert(b->type < DIV_NUM_TYPES);
				color = plot_divider_types[b->type].get_avg_color(); return 1;
			}
		}
	} // for i
	start_ix = 0;

	for (auto i = pool_groups.begin(); i != pool_groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy(pos)) continue;
		assert(start_ix <= i->ix && i->ix <= pools.size());

		for (auto b = pools.begin()+start_ix; b != pools.begin()+i->ix; ++b) {
			if (pos.x < b->bcube.x1()) break; // pools are sorted by x1, no divider after this can match
			if (!b->bcube.contains_pt_xy(pos)) continue;
			if (b->above_ground && !dist_xy_less_than(pos, point(b->bcube.xc(), b->bcube.yc(), b->bcube.z1()), b->get_radius())) continue; // circular in-ground pool
			color = b->wcolor; // return water color
			return 1;
		}
	} // for i
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


