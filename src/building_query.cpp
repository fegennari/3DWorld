// 3D World - Building Query Functions
// by Frank Gennari 10/09/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"


extern bool draw_building_interiors, camera_in_building, player_near_toilet, player_in_unlit_room, building_has_open_ext_door, player_on_escalator, ctrl_key_pressed;
extern bool player_is_hiding, player_wait_respawn, had_building_interior_coll, player_in_int_elevator, player_in_walkway, player_on_house_stairs, player_on_moving_ww;
extern bool player_in_tunnel;
extern int camera_surf_collide, frame_counter, animate2, player_in_closet, player_in_elevator, player_in_basement, player_in_attic, player_in_water;
extern float CAMERA_RADIUS, C_STEP_HEIGHT, NEAR_CLIP, fticks, building_bcube_expand;
extern double camera_zh;
extern building_params_t global_building_params;
extern bldg_obj_type_t bldg_obj_types[];
extern building_t const *player_building;


float get_railing_height(room_object_t const &c);
cylinder_3dw get_railing_cylinder(room_object_t const &c);
bool sphere_vert_cylin_intersect_with_ends(point &center, float radius, cylinder_3dw const &c, vector3d *cnorm);
void register_in_closed_bathroom_stall();
pair<cube_t, colorRGBA> car_bcube_color_from_parking_space(room_object_t const &o);
void force_player_height(double height);
bool is_player_model_female();
bool get_sphere_poly_int_val(point const &sc, float sr, point const *const points, unsigned npoints, vector3d const &normal, float thickness, float &val, vector3d &cnorm);
float get_player_move_dist();
void get_shelf_brackets(room_object_t const &c, cube_t shelves[4], unsigned num_shelves, vect_cube_with_ix_t &brackets);


// assumes player is in this building; handles windows and exterior doors but not attics and basements
bool building_t::player_can_see_outside() const {
	vector3d const xlate(get_tiled_terrain_model_xlate());
	point const camera_pos(get_camera_pos()), camera_bs(camera_pos - xlate);
	float const floor_spacing(get_window_vspace());

	if (!has_int_windows()) { // no windows looking out
		// should we check if the door is visible?
		if (building_has_open_ext_door && !doors.empty()) { // maybe can see out a door
			// maybe can see out open door on first floor or a walkway
			if (have_walkway_ext_door || (camera_pos.z < (ground_floor_z1 + floor_spacing) && camera_pos.z > ground_floor_z1)) return 1;

			if (interior) { // check roof access stairs
				point camera_bot_bs(camera_bs);
				camera_bot_bs.z -= get_bldg_player_height();

				for (stairwell_t const &s : interior->stairwells) {
					if (s.roof_access && !s.is_u_shape() && s.contains_pt(camera_bot_bs)) return 1;
				}
			}
		}
		for (cube_t const &skylight : skylights) { // check skylights
			cube_t test_area(skylight);
			test_area.expand_by_xy(2.0*floor_spacing);
			test_area.z1() -= floor_spacing; // include the floor below
			if (test_area.contains_pt(camera_bs)) return 1; // outside may be visible through the skylight
		}
		if (player_can_see_out_mall_skylight(xlate)) return 1;
		return 0;
	}
	// maybe can see out a window
	if (!interior || !has_basement() || camera_pos.z > ground_floor_z1) return 1; // no interior, or not in the basement
	if (camera_pos.z < ground_floor_z1 - floor_spacing)     return 0; // on lower level of parking garage or extended basement
	if (point_in_extended_basement_not_basement(camera_bs)) return 0; // not visible from extended basement; not for rotated buildings

	for (auto const &s : interior->stairwells) { // check basement stairs
		if (s.z1() >= ground_floor_z1 || s.z2() <= ground_floor_z1) continue; // not basement stairs
		if (s.in_mall || s.in_ext_basement) continue; // mall and extended basement stairs don't count
		if (s.is_u_shape())                 continue; // can't see around the bend
		if (!is_rot_cube_visible(s, xlate)) continue; // VFC
		if (!s.contains_pt(camera_bs) && s.stairs_door_ix >= 0 && get_door(s.stairs_door_ix).open_amt == 0.0) continue; // closed stairs door, not visible
		if (!has_parking_garage) return 1; // assume exterior may be visible through normal stairs, but not parking garage stairs
		cube_t stairs_exp(s);
		stairs_exp.expand_by_xy(floor_spacing);
		if (!get_basement().contains_cube_xy(stairs_exp)) return 1; // stairs near the edge of the building basement - window may be visible
	} // for s
	if (has_pg_ramp() && !interior->ignore_ramp_placement) {} // what about ramp?
	return 0;
}
bool building_t::player_can_see_mall_skylight(vector3d const &xlate) const { // looking in or out
	for (cube_t const &skylight : interior->mall_info->skylights) {
		if (camera_pdu.cube_visible(skylight + xlate)) return 1;
	}
	return 0;
}
bool building_t::player_can_see_out_mall_skylight(vector3d const &xlate) const { // in mall looking out
	return (point_in_mall(camera_pos - xlate) && player_can_see_mall_skylight(xlate));
}
bool building_t::player_can_see_in_mall_skylight(vector3d const &xlate) const { // above mall looking in
	if (!has_mall_skylight()) return 0;
	point const cbs(camera_pos - xlate);
	return (cbs.z > ground_floor_z1 && point_in_mall(point(cbs.x, cbs.y, (ground_floor_z1 - 0.5*get_window_vspace()))) && player_can_see_mall_skylight(xlate));
}
bool player_in_windowless_building() {return (player_building != nullptr && !player_building->player_can_see_outside());}

bool player_cant_see_outside_building() {
	if (player_in_basement >= 3) { // player in extended basement; only outside view is through mall skylight
		return (!player_building || !player_building->player_can_see_out_mall_skylight(get_tiled_terrain_model_xlate()));
	}
	if (player_in_attic == 2) return 1; // player in windowless attic
	return player_in_windowless_building();
}
bool building_t::check_bcube_overlap_xy(building_t const &b, float expand_rel, float expand_abs) const {
	if (expand_rel == 0.0 && expand_abs == 0.0 && !bcube.intersects(b.bcube)) return 0;
	if (!is_rotated() && !b.is_rotated()) return 1; // above check is exact, top-level bcube check up to the caller
	if (b.bcube.contains_pt_xy(bcube.get_cube_center()) || bcube.contains_pt_xy(b.bcube.get_cube_center())) return 1; // slightly faster to include this check
	return (check_bcube_overlap_xy_one_dir(b, expand_rel, expand_abs) || b.check_bcube_overlap_xy_one_dir(*this, expand_rel, expand_abs));
}

// Note: only checks for point (x,y) value contained in one cube/N-gon/cylinder; assumes pt has already been rotated into local coordinate frame
bool building_t::check_part_contains_pt_xy(cube_t const &part, unsigned part_id, point const &pt) const {
	if (!part.contains_pt_xy(pt)) return 0; // check bounding cube
	if (is_cube())                return 1; // that's it
	vector<point> const &points(get_part_ext_verts(part_id));
	return point_in_polygon_2d(pt.x, pt.y, points.data(), points.size()); // 2D x/y containment
}
bool building_t::check_part_contains_cube_xy(cube_t const &part, unsigned part_id, cube_t const &c) const { // for placing roof objects
	for (unsigned i = 0; i < 4; ++i) { // check cylinder/ellipse
		point const pt(c.d[0][i&1], c.d[1][i>>1], 0.0); // XY only
		if (!check_part_contains_pt_xy(part, part_id, pt)) return 0;
	}
	return 1;
}
bool building_t::cube_int_parts_no_sec(cube_t const &c) const {
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		if (p->intersects_no_adj(c)) return 1;
	}
	return 0;
}
bool building_t::check_pt_in_retail_room(point const &p) const {
	if (!has_retail() || !interior || interior->rooms.empty()) return 0;
	return get_retail_room().contains_pt(p);
}
bool building_t::check_pt_in_or_near_walkway(point const &p, bool owned_only, bool inc_open_door, bool inc_conn_room) const {
	for (building_walkway_t const &w : walkways) {
		if (owned_only && !w.is_owner) continue;
		
		if (inc_open_door && w.has_skyway_conn()) { // test skyway walkway connection
			cube_t vis_area(w.skyway_conn);
			vis_area.expand_in_dim(w.dim, (player_in_walkway ? 2.0 : 0.1)*w.get_length()); // extend to include other nearby walkways
			if (vis_area.contains_pt(p)) return 1; // player in nearby skyway area or opposing walkway
		}
		if (p.z > w.bcube.z1() && p.z < w.bcube.z2()) { // check for Z overlap
			if ((inc_open_door ? w.get_bcube_inc_open_door() : w.bcube).contains_pt(p)) return 1;
			// check adjacent rooms if connected, walkway is visible through windows, and pos is off the ends of the walkways; forms a dog bone shape
			if (inc_conn_room && has_int_windows() && (p[w.dim] < w.bcube.d[w.dim][0] || p[w.dim] > w.bcube.d[w.dim][1]) && w.bcube_inc_rooms.contains_pt(p)) return 1;
		}
		if (!w.elevator_bcube.is_all_zeros() && w.elevator_bcube.contains_pt(p)) return 1; // check elevator connected to walkway
	} // for walkways
	return 0;
}
bool building_t::is_connected_with_walkway(building_t const &target, float zval) const {
	for (building_walkway_t const &w : walkways) {
		if (zval != 0.0 && (zval < w.bcube.z1() || zval > w.bcube.z2())) continue; // apply zval filter
		if (w.conn_bldg == &target) return 1;
	}
	return 0;
}

cube_t building_walkway_t::get_bcube_inc_open_door() const {
	cube_t bc(bcube);
	bc.expand_in_dim(dim, get_door_open_dist());
	return bc;
}

bool building_t::check_bcube_overlap_xy_one_dir(building_t const &b, float expand_rel, float expand_abs) const { // can be called before levels/splits are created
	// Note: easy cases are handled by check_bcube_overlap_xy() above
	point const center1(b.bcube.get_cube_center()), center2(bcube.get_cube_center());

	for (auto p1 = b.parts.begin(); p1 != b.parts.end(); ++p1) {
		point pts[9]; // {center, 00, 10, 01, 11, x0, x1, y0, y1}

		if (b.parts.size() == 1) {pts[0] = center1;} // single cube: we know we're rotating about its center
		else {
			pts[0] = p1->get_cube_center();
			b.do_xy_rotate(center1, pts[0]); // rotate into global space
		}
		cube_t c_exp(*p1);
		c_exp.expand_by_xy(expand_rel*p1->get_size() + vector3d(expand_abs, expand_abs, expand_abs));

		for (unsigned i = 0; i < 4; ++i) { // {00, 10, 01, 11}
			pts[i+1].assign(c_exp.d[0][i&1], c_exp.d[1][i>>1], 0.0); // XY only
			b.do_xy_rotate(center1, pts[i+1]); // rotate into global space
		}
		for (unsigned i = 0; i < 5; ++i) {do_xy_rotate_inv(center2, pts[i]);} // inverse rotate into local coord space - negate the sine term
		cube_t c_exp_rot(pts+1, 4); // use points 1-4
		pts[5] = 0.5*(pts[1] + pts[3]); // x0 edge center
		pts[6] = 0.5*(pts[2] + pts[4]); // x1 edge center
		pts[7] = 0.5*(pts[1] + pts[2]); // y0 edge center
		pts[8] = 0.5*(pts[3] + pts[4]); // y1 edge center

		for (auto p2 = parts.begin(); p2 != parts.end(); ++p2) {
			if (c_exp_rot.contains_pt_xy(p2->get_cube_center())) return 1; // quick and easy test for heavy overlap

			for (unsigned i = 0; i < 9; ++i) {
				if (check_part_contains_pt_xy(*p2, (p2 - parts.begin()), pts[i])) return 1; // Note: building geometry is likely not yet generated, this check should be sufficient
				//if (p2->contains_pt_xy(pts[i])) return 1;
			}
		}
	} // for p1
	return 0;
}

bool do_sphere_coll_polygon_sides(point &pos, cube_t const &part, float radius, bool interior_coll, vector<point> const &points, vector3d *cnorm) {
	point quad_pts[4]; // quads
	bool coll(0);

	for (unsigned S = 0; S < points.size(); ++S) { // generate vertex data quads
		for (unsigned d = 0, ix = 0; d < 2; ++d) {
			point const &p(points[(S+d)%points.size()]);
			for (unsigned e = 0; e < 2; ++e) {quad_pts[ix++].assign(p.x, p.y, part.d[2][d^e]);}
		}
		vector3d const normal((interior_coll ? -1.0 : 1.0)*get_poly_norm(quad_pts)); // invert to get interior normal
		float const rdist(dot_product_ptv(normal, pos, quad_pts[0]));
		if (rdist < 0.0 || rdist >= radius) continue; // too far or wrong side
		if (!sphere_poly_intersect(quad_pts, 4, pos, normal, rdist, radius)) continue;
		pos += normal*(radius - rdist);
		if (cnorm) {*cnorm = normal;}
		if (!interior_coll) return 1; // interior coll must test against all sides
		coll = 1;
	} // for S
	return coll;
}

// called for players and butterfiles
bool building_t::test_coll_with_sides(point &pos, point const &p_last, float radius, vector3d const &xlate, cube_t const &part, unsigned part_id, vector3d *cnorm) const {

	vect_point const &points(get_part_ext_verts(part_id));
	unsigned const num_steps(max(1U, min(100U, (unsigned)ceil(2.0*p2p_dist_xy(pos, p_last)/radius))));
	vector3d const step_delta((pos - p_last)/num_steps);
	pos = p_last - xlate;

	// if the player is moving too quickly, the intersection with a side polygon may be missed, which allows the player to travel through the building;
	// so we split the test into multiple smaller sphere collision steps
	for (unsigned step = 0; step < num_steps; ++step) {
		pos += step_delta;
		if (do_sphere_coll_polygon_sides(pos, part, radius, 0, points, cnorm)) {pos += xlate; return 1;} // interior_coll=0
	}
	if (max(pos.z, p_last.z) > part.z2() && point_in_polygon_2d(pos.x, pos.y, points.data(), num_sides)) { // test top plane (sphere on top of polygon?)
		pos.z = part.z2() + radius; // make sure it doesn't intersect the roof
		pos  += xlate;
		if (cnorm) {*cnorm = plus_z;}
		return 1;
	}
	pos += xlate;
	return 0;
}

cube_t building_t::get_interior_bcube(bool inc_ext_basement) const { // Note: called for indir lighting; could cache z2
	cube_t int_bcube(inc_ext_basement ? get_bcube_inc_extensions() : bcube);
	int_bcube.z2() = interior_z2;
	return int_bcube;
}
cube_t building_t::get_ext_vis_bcube() const {
	if (!has_mall_skylight()) return bcube;
	cube_t ext_vis_bc(bcube);
	for (cube_t const &s : interior->mall_info->skylights) {ext_vis_bc.union_with_cube(s);}
	return ext_vis_bc;
}
void building_t::union_with_coll_bcube(cube_t const &c) {
	coll_bcube.union_with_cube(c);
	for (unsigned d = 0; d < 2; ++d) {max_eq(building_bcube_expand, max((bcube.d[d][0] - c.d[d][0]), (c.d[d][1] - bcube.d[d][1])));}
}

void accumulate_shared_xy_area(cube_t const &c, cube_t const &sc, float &area) {
	if (c.intersects_xy(sc)) {area += (min(c.x2(), sc.x2()) - max(c.x1(), sc.x1()))*(min(c.y2(), sc.y2()) - max(c.y1(), sc.y1()));}
}

void apply_speed_factor(point &pos, point const &p_last, float speed_factor) {
	for (unsigned d = 0; d < 2; ++d) {pos[d] = speed_factor*pos[d] + (1.0 - speed_factor)*p_last[d];}
}

// Note: used for the player when check_interior=1; pos and p_last are in camera space
bool building_t::check_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only, vector3d *cnorm_ptr, bool check_interior) const {
	if (!is_valid()) return 0; // invalid building
	
	if (radius > 0.0) {
		// player can't be in multiple buildings at once; if they were in some other building last frame, they can't be in this building this frame
		if (check_interior && camera_in_building && player_building != this) return 0;

		// skip cube intersection test if point (camera) is in the mall elevator; use p_last for this test since we want to know where the player was the last frame
		if (!check_interior || !point_in_mall_elevator_entrance((p_last - xlate), 1)) {
			// skip zval check when the player is in the building because pos.z may be placed on the mesh,
			// which could be above the top of the building when the player is in the extended basement;
			// this should be legal because the player can't exit the building in the Z direction; maybe better to use p_last.z?
			bool const sc_xy_only(xy_only || (check_interior && camera_in_building && player_building == this));
			cube_t test_cube(coll_bcube);
			point const pos_bs(pos - xlate);
			if (check_interior && has_ext_basement()) {test_cube.union_with_cube(interior->basement_ext_bcube);} // include the extended basement for interior checks
			if (!(sc_xy_only ? sphere_cube_intersect_xy(pos_bs, radius, test_cube) : sphere_cube_intersect(pos_bs, radius, test_cube))) return 0;
		}
	}
	if (check_sphere_coll_inner(pos, p_last, xlate, radius, xy_only, cnorm_ptr, check_interior)) return 1;
	if (!check_interior) return 0;
	building_t const *const other_bldg(get_conn_bldg_for_pt((point(pos.x, pos.y, min(pos.z, p_last.z)) - xlate), radius)); // Note: not handling pos rotation
	if (other_bldg == nullptr) return 0; // no connected building at this location
	return other_bldg->check_sphere_coll_inner(pos, p_last, xlate, radius, xy_only, cnorm_ptr, check_interior); // check connected building
}
bool building_t::check_sphere_coll_inner(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only, vector3d *cnorm_ptr, bool check_interior) const {
	float const xy_radius(radius*global_building_params.player_coll_radius_scale);
	point pos2(pos), p_last2(p_last), center;
	bool had_coll(0), is_interior(0), is_in_attic(0), allow_outside_building(0), on_ext_stair(0), reduce_speed(0);
	float part_z2(bcube.z2());

	if (is_rotated()) {
		center = bcube.get_cube_center() + xlate;
		do_xy_rotate_inv(center, pos2); // inverse rotate - negate the sine term
		do_xy_rotate_inv(center, p_last2);
	}
	if (check_interior && draw_building_interiors && interior != nullptr) { // check for interior case first
		// this is the real zval for use in collsion detection, in building space
		float zval(((p_last2.z < ground_floor_z1) ? min(pos2.z, p_last2.z) : max(pos2.z, p_last2.z)) - xlate.z); // min/max away from the ground floor
		// if the player is in flight mode and enables collision detection at the ground floor or basement with their feet just below the floor,
		// clamp to Z1 to keep them inside the building rather than putting them on the roof
		if (zval < bcube.z1() && zval+camera_zh > bcube.z1() && bcube.contains_pt_xy(pos2 - xlate) && p_last == p_last2) {zval = bcube.z1();}
		point const pos2_bs(pos2 - xlate), query_pt(pos2_bs.x, pos2_bs.y, zval);

		// first check uses min of the two zvals to reject the basement, which is actually under the mesh
		if ((min(pos2.z, p_last2.z) + radius) > ground_floor_z1) {
			unsigned const floor_ix((zval - ground_floor_z1)/get_window_vspace());

			if (have_walkway_ext_door || (floor_ext_door_mask & (1 << floor_ix))) { // check for exterior doors on this floor
				for (auto d = doors.begin(); d != doors.end(); ++d) { // exterior doors
					if (d->type == tquad_with_ix_t::TYPE_RDOOR) continue; // doesn't apply to roof door
					cube_t bc(d->get_bcube());
					if (zval < bc.z1() || zval > bc.z2()) continue; // no Z overlap
					bool const door_dim(bc.dy() < bc.dx());
					bc.expand_in_dim( door_dim, 1.1*radius); // expand by radius plus some tolerance in door dim
					bc.expand_in_dim(!door_dim, -0.5*xy_radius); // shrink slightly in the other dim to prevent the player from clipping through the wall next to the door
					bc.z1() -= max(radius, (float)camera_zh); // account for player on steep slope up to door - require player head above doorframe bottom
					if (!bc.contains_pt(point(pos2_bs.x, pos2_bs.y, zval))) continue; // check if door can be used
					//if (floor_ix == 0) return 0; // ground floor, disable collsion to allow the player to walk through
					allow_outside_building = 1; // allow the player to exit, but for now still in the building
				} // for d
			}
		}
		// check for sphere contained in parts and maybe clamp to parts
		cube_t sc; sc.set_from_sphere(query_pt, radius); // sphere bounding cube; zvals are ignored

		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // include garages and sheds
			float cont_area(0.0);
			cube_t clamp_part(*i);

			if (is_basement(i)) {
				bool const in_ext_basement(point_in_extended_basement_not_basement(query_pt));
				if (in_ext_basement) {clamp_part = interior->basement_ext_bcube;}
				else if (!i->contains_pt(query_pt)) continue; // not in basement
				accumulate_shared_xy_area(*i, sc, cont_area);
				
				if (has_ext_basement() && zval < interior->basement_ext_bcube.z2()) { // below ext basement top (not in upper level of two level parking garage)
					// use the ext basement hallway if pos is in the basement, otherwise use the entire ext basement
					cube_t const &basement_cube(in_ext_basement ? interior->basement_ext_bcube : get_ext_basement_entrance());
					accumulate_shared_xy_area(basement_cube, sc, cont_area);
				}
			}
			else {
				if (!i->contains_pt(query_pt)) continue; // not interior to this part

				if (!allow_outside_building && !is_cube() && zval > i->z1() && zval < i->z2()) { // non-cube shaped building, in Z bounds, clamp_part is conservative
					//if (use_cylinder_coll()) {}
					vect_point const &points(get_part_ext_verts(i - parts.begin()));
					point const p_last2_bs(p_last2 - xlate);
					if (!point_in_polygon_2d(p_last2_bs.x, p_last2_bs.y, points.data(), points.size())) continue; // outside the building, even though inside part bcube
					pos2 -= xlate; // temporarily convert to building space
					bool const had_int(do_sphere_coll_polygon_sides(pos2, *i, radius, 1, points, cnorm_ptr));
					pos2 += xlate; // convert back to camera space
					if (had_int) {is_interior = had_coll = 1; break;} // interior_coll=1
				}
				for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
					if (zval < p->z1() || zval > p->z2()) continue; // wrong floor/part in stack
					accumulate_shared_xy_area(*p, sc, cont_area);
				}
			}
			if (!allow_outside_building && cont_area < 0.99*sc.get_area_xy()) { // sphere bounding cube not contained in union of parts - partially outside the building
				cube_t c(clamp_part + xlate); // convert to camera space
				c.expand_by_xy(-radius); // shrink part by sphere radius
				c.clamp_pt_xy(pos2); // force pos2 into interior of the cube to prevent the sphere from intersecting the part
				had_coll = 1;
			}
			is_interior = 1;
			break; // flag for interior collision detection
		} // for i
		if (point_in_attic(query_pt)) {is_interior = is_in_attic = 1;}
		is_interior |= point_in_mall_elevator_entrance(query_pt, 0); // inc_front_space=0

		if (!is_interior && has_mall() && interior->mall_info->city_elevator_ix >= 0) { // handle mall elevator exterior coll when doors are closed
			elevator_t const &e(get_elevator(interior->mall_info->city_elevator_ix));
			if (e.open_amt < 0.5) {had_coll |= sphere_cube_int_update_pos(pos2, radius, (e + xlate), p_last2, xy_only, cnorm_ptr);}
		}
		if (!is_interior) { // not interior to a part - check roof access
			float const floor_thickness(get_floor_thickness());

			for (auto i = interior->stairwells.begin(); i != interior->stairwells.end(); ++i) {
				if (!i->roof_access) continue;
				cube_t test_cube(*i);
				test_cube.expand_by_xy(-0.5f*radius);
				if (!test_cube.contains_pt_xy(pos2 - xlate)) continue; // pos not over stairs
				if (zval < i->z2() + radius + floor_thickness) {is_interior = 1; break;} // Note: don't have to check zval > i->z2() because we know that !is_interior
			}
		}
	}
	if (!xy_only && min(pos2.z, p_last2.z) > ground_floor_z1) { // only when above ground; skip if in basement
		// check for collision with exterior stairs, since they apply to both the interior and exterior case
		for (ext_step_t const &s : ext_steps) { // Note: includes balconies
			cube_t const c(s + xlate);
			if (s.enclosed && !c.contains_pt_xy(pos2))      continue; // don't clip through walls onto balconies
			if (!sphere_cube_intersect_xy(pos2, radius, c)) continue;
			float const zval(max(pos2.z, p_last2.z));
			// use 2x radius for stability: too low and we fall through the step some frames; too high and we get stuck on stairs below
			if (zval + radius < c.z1() || zval - 2.0*radius > c.z2()) continue; // no collision in Z

			if (!s.at_door && !s.is_base && (pos2[!s.dim] < c.d[!s.dim][0] || pos2[!s.dim] > c.d[!s.dim][1])) { // player to the side
				had_coll |= sphere_cube_int_update_pos(pos2, radius, c, p_last2, xy_only, cnorm_ptr); // check collision with sides
				continue;
			}
			float &wval(pos2[!s.dim]), &lval(pos2[s.dim]); // position relative to the {width, length) of the stairs

			if (c.contains_pt_xy(pos2)) {
				// handle collision with railing and building wall
				unsigned skip_dir(0);
				if (s.at_door) {skip_dir |= (1 << unsigned(  s.wall_dir));} // allow entering/exiting from house
				if (s.is_base) {skip_dir |= (1 << unsigned(1-s.wall_dir));} // allow entering/exiting from ground
				if (!(skip_dir & 1)) {max_eq(wval, (c.d[!s.dim][0] + xy_radius));}
				if (!(skip_dir & 2)) {min_eq(wval, (c.d[!s.dim][1] - xy_radius));}
			
				if (s.at_door || s.enclosed) { // handle collision with top end railing(s)
					if ( s.step_dir || s.enclosed) {max_eq(lval, (c.d[s.dim][0] + xy_radius));}
					if (!s.step_dir || s.enclosed) {min_eq(lval, (c.d[s.dim][1] - xy_radius));}
				}
			}
			if (zval + radius > c.z1()) { // step up
				if (s.step_up || c.contains_pt_xy(pos2)) {max_eq(pos2.z, (c.z2() + radius));} // only step up if on this stair, unless it's a wall we can climb
				if (wval > c.d[!s.dim][0] && wval < c.d[!s.dim][1]) {on_ext_stair = reduce_speed = 1;} // not entering or leaving from the sides
				had_coll = 1;
			}
			else {had_coll |= sphere_cube_int_update_pos(pos2, radius, c, p_last2, xy_only, cnorm_ptr);} // check collision with sides
		} // for i
		if (!ladder.is_all_zeros() && sphere_cube_intersect(pos2, radius, ladder)) {
			//had_coll = sphere_cube_int_update_pos(pos2, radius, ladder, p_last2, xy_only, cnorm_ptr); // no, don't want to collide when player is on the roof
			pos2.z   = p_last2.z + 0.5*get_player_move_dist()*cview_dir.z; // move up/down based on player vertical view (looking up vs. down)
			pos2.z   = min((ladder.z2() + radius), max((ladder.z1() + radius), pos2.z)); // clamp to ladder height range
			had_coll = 1;
		}
	}
	if (player_wait_respawn) {pos2.x = p_last2.x; pos2.y = p_last2.y;} // player can't move
	else if (reduce_speed) {apply_speed_factor(pos2, p_last2, 0.6);} // slow down to 60% when on exterior stairs or balconies

	if (on_ext_stair) { // only need to check for blockers at the bottom of stairs
		for (auto const &i : details) {
			if (i.type == DETAIL_OBJ_SHAD_ONLY) continue; // not a collider
			had_coll |= sphere_cube_int_update_pos(pos2, radius, (i + xlate), p_last2, xy_only, cnorm_ptr); // treat as cubes
		}
	}
	else if (is_interior) {
		point pos2_bs(pos2 - xlate);
		if (check_sphere_coll_interior(pos2_bs, (p_last2 - xlate), radius, is_in_attic, xy_only, cnorm_ptr)) {pos2 = pos2_bs + xlate; had_coll = 1;}
		had_building_interior_coll = 1; // required to skip city object collisions
	}
	else { // exterior
		if (allow_outside_building) return 0; // entering an exterior door - no collision with building exterior
		float const r_bias(0.5*get_window_vspace()); // enough to allow some falling but not pick up the floor below

		for (auto i = parts.begin(); i != parts.end(); ++i) {
			if (xy_only) {
				if (i->z1() < ground_floor_z1) continue; // skip basements, since they should be contained in the union of the ground floor (optimization)
				
				if (i->z1() > ground_floor_z1) { // only check ground floor
					if (has_complex_floorplan) {continue;} else {break;}
				}
			}
			// zval must include p_last2 for correct part selection in complex building floorplans (with flat roofs), but is not correct for walking on sloped roofs
			float const zval((roof_type == ROOF_TYPE_FLAT) ? max(pos2.z, p_last2.z) : pos2.z);
			cube_t const part_bc(*i + xlate);
			if (!xy_only && ((zval + radius < part_bc.z1()) || (zval - radius - r_bias > part_bc.z2()))) continue; // test z overlap
			if (radius == 0.0 && !(xy_only ? part_bc.contains_pt_xy(pos2) : part_bc.contains_pt(point(pos2.x, pos2.y, zval)))) continue; // no intersection; ignores p_last
			unsigned const part_id(i - parts.begin());
			bool part_coll(0);

			if (use_cylinder_coll()) {
				point const cc(part_bc.get_cube_center());
				float const crx(0.5*i->dx()), cry(0.5*i->dy()), r_sum(radius + max(crx, cry));
				if (!dist_xy_less_than(pos2, cc, r_sum)) continue; // no intersection

				if (fabs(crx - cry) < radius) { // close to a circle
					if (p_last2.z > part_bc.z2() && dist_xy_less_than(pos2, cc, max(crx, cry))) {
						pos2.z = part_bc.z2() + radius; // make sure it doesn't intersect the roof
						if (cnorm_ptr) {*cnorm_ptr = plus_z;}
					}
					else { // side coll
						vector2d const d((pos2.x - cc.x), (pos2.y - cc.y));
						float const mult(r_sum/d.mag());
						pos2.x = cc.x + mult*d.x;
						pos2.y = cc.y + mult*d.y;
						if (cnorm_ptr) {*cnorm_ptr = vector3d(d.x, d.y, 0.0).get_norm();} // no z-component
					}
					part_coll = 1;
				}
				else {
					part_coll |= test_coll_with_sides(pos2, p_last2, radius, xlate, part_bc, part_id, cnorm_ptr); // use polygon collision test
				}
			}
			else if (num_sides != 4) { // triangle, hexagon, octagon, etc.
				part_coll |= test_coll_with_sides(pos2, p_last2, radius, xlate, part_bc, part_id, cnorm_ptr);
			}
			else if (!xy_only && part_bc.contains_pt_xy_exp(pos2, radius) && p_last2.z > (i->z2() + xlate.z)) { // on top of building
				pos2.z = i->z2() + xlate.z + radius;
				if (cnorm_ptr) {*cnorm_ptr = plus_z;}
				part_coll = 1;
			}
			else if (sphere_cube_int_update_pos(pos2, radius, part_bc, p_last2, xy_only, cnorm_ptr)) { // cube
				part_coll = 1; // flag as colliding, continue to look for more collisions (inside corners)
			}
			if (part_coll && pos2.z < part_bc.z1()) {pos2.z = part_bc.z2() + radius;} // can't be under a building - make it on top of the building instead
			if (part_coll) {part_z2 = i->z2();}
			had_coll |= part_coll;
		} // for i
		if (has_chimney == 2) {had_coll |= sphere_cube_int_update_pos(pos2, radius, (get_fireplace() + xlate), p_last2, xy_only, cnorm_ptr);} // exterior fireplace
		for (cube_t const &fence : fences) {had_coll |= sphere_cube_int_update_pos(pos2, radius, (fence + xlate), p_last2, xy_only, cnorm_ptr);}

		// Note: driveways are handled elsewhere in the control flow
		if (!xy_only) { // don't need to check details and roof in xy_only mode because they're contained in the XY footprint of the parts (except balconies and stairs)
			for (auto const &i : details) {
				if (i.type == DETAIL_OBJ_SHAD_ONLY) continue; // not a collider

				if (i.type == ROOF_OBJ_WTOWER && pos2.z < i.z1() + 0.25*i.dz()) { // below custom water tower collision
					if (!sphere_cube_intersect(pos2, radius, (i + xlate))) continue;
					float const hwidth(0.25*(i.dx() + i.dy()));
					cube_t inner(i);
					inner.expand_by_xy(-0.1*hwidth);

					for (unsigned n = 0; n < 4; ++n) {
						bool const xd(n & 1), yd(n >> 1);
						cube_t leg(i);
						leg.d[0][!xd] = inner.d[0][xd];
						leg.d[1][!yd] = inner.d[1][yd];
						had_coll |= sphere_cube_int_update_pos(pos2, radius, (leg + xlate), p_last2, xy_only, cnorm_ptr);
					}
					cylinder_3dw const cylin(cube_bot_center(i)+xlate, cube_top_center(i)+xlate, 0.1*hwidth, 0.1*hwidth);
					had_coll |= sphere_vert_cylin_intersect(pos2, radius, cylin, cnorm_ptr);
				}
				else {had_coll |= sphere_cube_int_update_pos(pos2, radius, (i + xlate), p_last2, xy_only, cnorm_ptr);} // treat as cubes
			} // for i
			for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
				point const pos_xlate(pos2 - xlate);

				if (check_interior && had_coll && pos2.z - xlate.z > part_z2) { // player standing on top of a building with a sloped roof or roof access cover
					if (point_in_polygon_2d(pos_xlate.x, pos_xlate.y, i->pts, i->npts)) {
						vector3d const normal(i->get_norm());
						if (normal.z == 0.0) continue; // skip vertical sides as the player can't stand on them
						float const rdist(dot_product_ptv(normal, pos_xlate, i->pts[0]));

						if (draw_building_interiors && i->type == tquad_with_ix_t::TYPE_ROOF_ACC) { // don't allow walking on roof access tquads
							if (rdist < -radius*normal.z) continue; // player is below this tquad
							else {pos2.x = p_last2.x; pos2.y = p_last2.y; break;} // block the player from walking here (can only walk through raised opening)
						}
						pos2.z += (radius - rdist)/normal.z; // determine the distance we need to move vertically to achieve this diag separation
					}
				}
				else { // normal case for bouncing object, etc.
					vector3d const normal(i->get_norm());
					float const rdist(dot_product_ptv(normal, pos_xlate, i->pts[0]));

					if (fabs(rdist) < radius && sphere_poly_intersect(i->pts, i->npts, pos_xlate, normal, rdist, radius)) {
						pos2 += normal*(radius - rdist); // update current pos
						if (cnorm_ptr) {*cnorm_ptr = ((normal.z < 0.0) ? -1.0 : 1.0)*normal;} // make sure normal points up
						had_coll = 1; // flag as colliding
						break; // only use first colliding tquad
					}
				}
			} // for i
		}
		if (is_house && has_porch() && sphere_cube_int_update_pos(pos2, radius, (porch + xlate), p_last2, xy_only, cnorm_ptr)) {had_coll = 1;} // porch
	} // end !is_interior case
	if (!had_coll) return 0; // Note: no collisions with windows or doors, since they're colinear with walls

	if (is_rotated()) {
		do_xy_rotate(center, pos2); // rotate back around center
		if (cnorm_ptr) {do_xy_rotate_normal(*cnorm_ptr);} // rotate normal back
	}
	pos = pos2;
	return 1;
}

float room_object_t::get_radius() const {
	if (shape == SHAPE_CYLIN || shape == SHAPE_VERT_TORUS) {return 0.25f*(dx() + dy());} // vertical cylinder or torus: return average of x/y diameter
	if (shape == SHAPE_SPHERE) {return 0.5*dx();} // sphere, should be the same dx()/dy()/dz() value (but can't assert due to FP precision errors)
	std::cerr << "Error: get_radius() called on non-cylinder/cube object of type " << int(type) << endl;
	assert(0); // cubes don't have a radius
	return 0.0; // never gets here
}
cylinder_3dw room_object_t::get_cylinder() const {
	float const radius(get_radius());
	point const center(get_cube_center());
	return cylinder_3dw(point(center.x, center.y, z1()), point(center.x, center.y, z2()), radius, radius);
}

// Note: returns bit vector for each cube that collides; supports up to 32 cubes
unsigned check_cubes_collision(cube_t const *const cubes, unsigned num_cubes, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	unsigned coll_ret(0);

	for (unsigned n = 0; n < num_cubes; ++n) {
		if (!cubes[n].is_all_zeros() && sphere_cube_int_update_pos(pos, radius, cubes[n], p_last, 0, cnorm)) {coll_ret |= (1<<n);} // skip_z=0
	}
	return coll_ret;
}
unsigned get_closet_num_coll_cubes(room_object_t const &c) {
	// include closed closet door for large closets; skip collision check of open doors for large closets since this case is more complex
	return ((!c.is_open() && !c.is_small_closet()) ? 5U : 4U);
}
unsigned check_closet_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[5];
	get_closet_cubes(c, cubes, 1); // get cubes for walls and door; required to handle collision with closet interior; for_collision=1
	return check_cubes_collision(cubes, get_closet_num_coll_cubes(c), pos, p_last, radius, cnorm);
}
unsigned check_bed_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
	get_bed_cubes(c, cubes);
	unsigned num_to_check(5); // skip legs_bcube
	if (c.taken_level > 0) {--num_to_check;} // skip pillows
	if (c.taken_level > 2) {--num_to_check;} // skip mattress
	unsigned coll_ret(check_cubes_collision(cubes, num_to_check, pos, p_last, radius, cnorm));
	get_tc_leg_cubes(cubes[5], c, BED_HEAD_WIDTH, 0, cubes); // cubes[5] is not overwritten
	coll_ret |= (check_cubes_collision(cubes, 4, pos, p_last, radius, cnorm) << 5); // check legs
	return coll_ret;
}
bool can_use_table_coll(room_object_t const &c) {
	return (c.type == TYPE_DESK || c.type == TYPE_DRESSER || c.type == TYPE_NIGHTSTAND || c.type == TYPE_TABLE);
}
// actually applies to tables, desks, dressers, and nightstands
unsigned get_table_like_object_cubes(room_object_t const &c, cube_t cubes[7]) { // tables, desks, dressers, and nightstands
	if (c.type == TYPE_TABLE && c.item_flags > 0) {
		get_cubes_for_plastic_table(c, 0.12, cubes); // top_dz=0.12
		return 3; // {top, vert, base}
	}
	unsigned num(5);
	get_table_cubes(c, cubes); // top and 4 legs
	if (c.type == TYPE_DRESSER || c.type == TYPE_NIGHTSTAND) {cubes[num++] = get_dresser_middle(c);}
	else if (c.type == TYPE_DESK) {
		if (c.desk_has_drawers() ) {cubes[num++] = get_desk_drawers_part(c);}
		if (c.shape == SHAPE_TALL) {cubes[num++] = get_desk_top_back    (c);} // tall desk
	}
	return num;
}
unsigned check_table_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[7];
	unsigned const num(get_table_like_object_cubes(c, cubes));
	return check_cubes_collision(cubes, num, pos, p_last, radius, cnorm);
}
unsigned check_conf_table_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[2]; // {top, base}
	get_conf_table_cubes(c, cubes);
	return check_cubes_collision(cubes, 2, pos, p_last, radius, cnorm);
}
unsigned check_rdesk_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[3]; // {front, left, right}
	get_reception_desk_cubes(c, cubes);
	return check_cubes_collision(cubes, 3, pos, p_last, radius, cnorm);
}
unsigned check_chair_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[3], leg_cubes[4]; // seat, back, legs_bcube
	get_chair_cubes(c, cubes);
	get_tc_leg_cubes(cubes[2], c, c.get_chair_leg_width(), 1, leg_cubes);
	return (check_cubes_collision(cubes, 2, pos, p_last, radius, cnorm) | check_cubes_collision(leg_cubes, 4, pos, p_last, radius, cnorm));
}

float get_tub_water_level(room_object_t const &c);
void get_tub_cubes(room_object_t const &c, cube_t cubes[5], float radius_if_ball=0.0) { // bottom, front, back, 2 sides
	float const height(c.get_height()), width(c.get_width()), signed_depth((c.dir ? 1.0 : -1.0)*c.get_depth());
	cube_t bottom(c), front(c), back(c);
	float bot_z2(c.z1() + 0.07*height);
	// if the collider is a ball, assume density=0.5 and it floats centered at the water level; not using ball_type.density here
	if (radius_if_ball > 0.0) {max_eq(bot_z2, (c.z1() + min(get_tub_water_level(c), 1.0f)*c.dz() - radius_if_ball));}
	bottom.z2() = front.z1() = back.z1() = bot_z2;
	front.d[c.dim][ c.dir] -= 0.81*signed_depth;
	back .d[c.dim][!c.dir] += 0.86*signed_depth;
	cubes[0] = bottom;
	cubes[1] = front;
	cubes[2] = back;

	for (unsigned d = 0; d < 2; ++d) { // sides
		cube_t side(c);
		side.z1() = bottom.z2();
		side.d[!c.dim][d] -= (d ? 1.0 : -1.0)*0.85*width;
		side.d[c.dim][!c.dir] = front.d[c.dim][ c.dir];
		side.d[c.dim][ c.dir] = back .d[c.dim][!c.dir];
		cubes[d+3] = side;
	}
}
unsigned check_tub_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm, bool is_ball) {
	cube_t cubes[5]; // bottom, front, back, 2 sides
	get_tub_cubes(c, cubes, (is_ball ? radius : 0.0));
	return check_cubes_collision(cubes, 5, pos, p_last, radius, cnorm);
}
unsigned check_bookcase_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	if (c.is_on_floor()) return sphere_cube_int_update_pos(pos, radius, c, p_last, 0, cnorm); // fallen over
	cube_t top, middle, back, lr[2];
	get_bookcase_cubes(c, top, middle, back, lr);
	cube_t const cubes[3] = {middle, lr[0], lr[1]}; // ignore top and back since they should be contained in middle
	return check_cubes_collision(cubes, 3, pos, p_last, radius, cnorm);
}
unsigned check_bench_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[4]; // seat, lo side, hi side, [back]
	unsigned const num(get_bench_cubes(c, cubes));
	return check_cubes_collision(cubes, num, pos, p_last, radius, cnorm);
}
unsigned check_diving_board_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[2]; // board, base
	get_diving_board_cubes(c, cubes);
	return check_cubes_collision(cubes, 2, pos, p_last, radius, cnorm);
}

bool check_ramp_collision(room_object_t const &c, point &pos, float radius, vector3d *cnorm) {
	if (!sphere_cube_intersect(pos, radius, c)) return 0;
	bool const is_pg_ramp(c.type == TYPE_RAMP);
	float const thickness(is_pg_ramp ? RAMP_THICKNESS_SCALE*c.dz() : 0.0), half_thickness(0.5*thickness);
	tquad_t ramp(get_ramp_tquad(c));
	for (unsigned n = 0; n < ramp.npts; ++n) {ramp.pts[n].z -= half_thickness;} // shift to the centerline, which is what get_sphere_poly_int_val() requires
	vector3d const normal(ramp.get_norm());
	float dist(0.0); // distance from sphere to ramp surface
	vector3d coll_normal;
	if (!get_sphere_poly_int_val(pos, radius, ramp.pts, ramp.npts, normal, thickness, dist, coll_normal)) return 0;
	if (cnorm) {*cnorm = coll_normal;}
	pos += coll_normal*dist;
	pos += vector3d(coll_normal.x, coll_normal.y, 0.0)*0.1*radius; // move downhill (slowly when rolling)
	return 1;
}

unsigned get_stall_cubes(room_object_t const &c, cube_t sides[3]) { // {side, side, door}
	sides[0] = sides[1] = sides[2] = c;
	float const width(c.get_width());
	sides[0].d[!c.dim][1] -= 0.95*width;
	sides[1].d[!c.dim][0] += 0.95*width;
	if (!c.is_open()) {sides[2].d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*0.975*c.get_length();} // include closed door
	return (c.is_open() ? 2U : 3U);
}
// Note: these next two are intended to be called when maybe_inside_room_object() returns true
bool check_stall_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t sides[3];
	unsigned const num_cubes(get_stall_cubes(c, sides));
	return check_cubes_collision(sides, num_cubes, pos, p_last, radius, cnorm);
}

unsigned get_shower_cubes(room_object_t const &c, cube_t sides[2]) {
	sides[0] = sides[1] = c;
	bool const door_dim(c.dx() < c.dy()), door_dir(door_dim ? c.dir : c.dim), side_dir(door_dim ? c.dim : c.dir); // {c.dim, c.dir} => {dir_x, dir_y}
	sides[0].d[!door_dim][side_dir] -= (side_dir ? 1.0 : -1.0)*0.95*c.get_sz_dim(!door_dim); // shrink to just the outer glass wall of the shower
	if (!c.is_open()) {sides[1].d[door_dim][door_dir] += (door_dir ? -1.0 : 1.0)*0.95*c.get_sz_dim(door_dim);} // check collision with closed door
	return (c.is_open() ? 1U : 2U);
}
bool check_shower_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t sides[2];
	unsigned const num_cubes(get_shower_cubes(c, sides));
	return check_cubes_collision(sides, num_cubes, pos, p_last, radius, cnorm);
}
bool check_balcony_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[4];
	get_balcony_cubes(c, cubes);
	return check_cubes_collision(cubes, 4, pos, p_last, radius, cnorm);
}

void get_pool_table_cubes(room_object_t const &c, cube_t cubes[5]) { // body + 4 legs
	cube_t legs_bcube(c);
	legs_bcube.z2() -= 0.48*c.dz();
	legs_bcube.expand_by_xy(-0.13*c.get_width()); // legs are recessed a bit, but not quite centered; should be close enough
	get_tc_leg_cubes(legs_bcube, c, 0.11, 0, cubes+1);  // legs are cubes[1]:cubes[4]
	cubes[0] = c;
	cubes[0].z1() = cubes[1].z2(); // body starts at the top of the legs
}
unsigned check_pool_table_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[5]; // body + 4 legs
	get_pool_table_cubes(c, cubes);
	// if pos is over the interior, use the top surface zval; otherwise, use z2 and assume an edge collision; could also call subtract_cube_from_cube()
	cube_t const top(get_pool_table_top_surface(c));
	if (top.contains_pt_xy(pos)) {cubes[0].z2() = top.z1();}
	return check_cubes_collision(cubes, 5, pos, p_last, radius, cnorm);
}
cube_t get_shelves_no_bot_gap(room_object_t const &c) {
	cube_t bcube(c);
	bcube.z1() += c.dz()/(c.get_num_shelves() + 1); // raise z1 to the bottom of the first shelf
	return bcube;
}

unsigned get_all_shelf_rack_cubes(room_object_t const &c, cube_t cubes[9]) { // back, <= 5 shelves, top, <= 2 sides
	cube_t back, top, sides[2], shelves[5];
	unsigned const num_shelves(get_shelf_rack_cubes(c, back, top, sides, shelves));
	unsigned num(0);
	cubes[num++] = back;
	for (unsigned n = 0; n < num_shelves; ++n) {cubes[num++] = shelves[n];}
	if (!top.is_all_zeros()) {cubes[num++] = top;}
	
	if (!sides[0].is_all_zeros()) {
		for (unsigned d = 0; d < 2; ++d) {cubes[num++] = sides[d];}
	}
	return num;
}
unsigned check_shelf_rack_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[9];
	unsigned const num_cubes(get_all_shelf_rack_cubes(c, cubes));

	if (c.is_nonempty()) { // if there are objects on the shelf rack, expand the back to 80% of the full shelf width
		set_wall_width(cubes[0], c.get_center_dim(c.dim), 0.8*0.5*c.get_depth(), c.dim);
	}
	return check_cubes_collision(cubes, num_cubes, pos, p_last, radius, cnorm);
}

unsigned get_couch_cubes(room_object_t const &c, cube_t cubes[4]) { // bottom, back, arm, arm; excludes pillows and blanket
	bool const dir2(c.dim ^ c.dir);
	float const height(c.dz()), depth(c.get_depth()), width(c.get_width());
	cube_t c_clip(c);
	c_clip.d[!c.dim][dir2] -= (dir2 ? 1.0 : -1.0)*0.025*width; // subtract off the blanket on the right side
	cube_t bot(c_clip), arms(c_clip);
	bot.expand_in_dim(!c.dim, -0.1*width); // shrink off arms
	cube_t back(bot);
	bot .z2() -= 0.48*height;
	bot .z1() += 0.05*height;
	arms.z2() -= 0.32*height;
	back.z1()  = bot.z2();
	back.d[c.dim][c.dir] -= (c.dir ? 1.0 : -1.0)*0.65*depth;
	cubes[0] = bot;
	cubes[1] = back;

	for (unsigned d = 0; d < 2; ++d) {
		cube_t arm(arms);
		arm.d[!c.dim][!d] = back.d[!c.dim][d];
		cubes[d+2] = arm;
	}
	return 4;
}
unsigned check_couch_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[4];
	return check_cubes_collision(cubes, get_couch_cubes(c, cubes), pos, p_last, radius, cnorm);
}

unsigned get_cashreg_cubes(room_object_t const &c, cube_t cubes[2]) { // main, screen
	cube_t const bc(get_true_room_obj_bcube(c));
	bool const dir2(c.dim ^ c.dir);
	float const width(bc.get_sz_dim(!c.dim)), length(bc.get_sz_dim(c.dim));
	cube_t body(bc), screen(bc);
	body.z2() = screen.z1() = c.z1() + 0.7*c.dz();
	screen.d[ c.dim][ c.dir] -= (c.dir ? 1.0 : -1.0)*0.395*length; // black screen side
	screen.d[ c.dim][!c.dir] += (c.dir ? 1.0 : -1.0)*0.36*length;  // shorter rectangle side
	screen.d[!c.dim][ dir2 ] -= (dir2  ? 1.0 : -1.0)*0.47*width;   // subtract off space in front of screen
	screen.d[!c.dim][!dir2 ] += (dir2  ? 1.0 : -1.0)*0.26*width;   // subtract off space behind of screen
	cubes[0] = body;
	cubes[1] = screen; // conservative
	return 2;
}
unsigned check_cashreg_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[2];
	return check_cubes_collision(cubes, get_cashreg_cubes(c, cubes), pos, p_last, radius, cnorm);
}

bool maybe_inside_room_object(room_object_t const &obj, point const &pos, float radius) {
	return ((obj.is_open() && sphere_cube_intersect(pos, radius, obj)) || obj.contains_pt(pos));
}
cube_t get_true_room_obj_bcube(room_object_t const &c) { // for collisions, etc.
	if (c.type == TYPE_ATTIC_DOOR) {
		cube_t bcube(get_attic_access_door_cube(c));
		if (c.is_open()) {bcube.union_with_cube(get_ladder_bcube_from_open_attic_door(c, bcube));} // include ladder as well
		return bcube;
	}
	if (c.type == TYPE_RAILING && (c.flags & RO_FLAG_IN_POOL)) { // pool stairs railing
		cube_t c_ext(c);
		c_ext.z2() += 0.5*c_ext.dz(); // extend the top upward to block the player from walking onto the railing
		return c_ext;
	}
	if (c.type == TYPE_RAILING) { // only works for straight (non-sloped) top railings
		cube_t c_ext(c);
		c_ext.z2() = c.z1() + get_railing_height(c);
		return c_ext;
	}
	if (c.type == TYPE_CASHREG) {
		cube_t c_ext(c);
		c_ext.expand_in_dim(!c.dim, -0.25*c.get_width()); // object bcube is twice the size it should be; subtract this off
		return c_ext;
	}
	if (c.type == TYPE_PLANT) {
		cube_t c_pot(c);
		c_pot.expand_by_xy(-(1.0 - PLANT_POT_RADIUS)*c.get_radius()); // use XY radius of the pot; better for AI coll
		return c_pot;
	}
	if (c.type == TYPE_TREE)      {} // bcubes are not a good fit for trees; cube covers the pot and trunk, but not leaves; should not be colliding with leaves anyway
	if (c.type == TYPE_SHOWERTUB) {return get_shower_tub_wall   (c);} // only the end wall is a collider; the tub handles the bottom (what about curtains?)
	if (c.type == TYPE_SHELVES  ) {return get_shelves_no_bot_gap(c);}
	return c; // default cube case
}

bool room_object_t::is_player_collidable() const {
	// chairs are player collidable only when in attics or backrooms; trashcans are only player collidable in malls
	return (!no_coll() && (bldg_obj_types[type].player_coll || (type == TYPE_CHAIR && (in_attic() || (flags & RO_FLAG_BACKROOM))) ||
		(type == TYPE_TCAN && (in_mall() || in_hallway()))));
}

// Note: used for the player; pos and p_last are already in rotated coordinate space
// default player is actually too large to fit through doors and too tall to fit between the floor and celing,
// so player size/height must be reduced in the config file
bool building_t::check_sphere_coll_interior(point &pos, point const &p_last, float radius, bool is_in_attic, bool xy_only, vector3d *cnorm) const {
	pos.z = get_bcube_z1_inc_ext_basement(); // start at building z1 rather than the terrain height in case we're at the foot of a steep hill
	assert(interior);
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness()), attic_door_z_gap(0.2f*floor_thickness);
	float const xy_radius(radius*global_building_params.player_coll_radius_scale); // XY radius can be smaller to allow player to fit between furniture
	bool had_coll(0), on_stairs(0), on_attic_ladder(0);
	float obj_z(max(pos.z, p_last.z)); // use p_last to get orig zval
	cube_with_ix_t const &attic_access(interior->attic_access);
	unsigned reset_to_last_dims(0); // {x, y} bit flags, for attic roof collision
	double camera_height(get_player_height());
	static double prev_camera_height(0.0);
	cube_t tunnel_walk_area;

	if (!xy_only && 2.2f*radius < (floor_spacing - floor_thickness)) { // diameter is smaller than space between floor and ceiling
		// check Z collision with floors; no need to check ceilings; this will set pos.z correctly so that we can set skip_z=0 in later tests
		float const floor_test_zval(obj_z + floor_thickness); // move up by floor thickness to better handle steep stairs
		float closest_floor_zval(pos.z - radius); // start at the player's feet

		for (auto i = interior->floors.begin(); i != interior->floors.end(); ++i) { // Note: includes attic and basement floors
			if (!i->contains_pt_xy(pos)) continue; // sphere not in this floor
			float const z(i->z2());
			if (z <= floor_test_zval && z > closest_floor_zval) {closest_floor_zval = z; had_coll = 1;} // move up
		}
		if (interior->room_geom) { // check glass floors
			for (cube_t const &f : interior->room_geom->glass_floors) {
				if (!f.contains_pt_xy(pos)) continue; // sphere not on this floor
				if (point_in_elevator(pos) || point_in_stairwell(pos)) continue; // elevators and stairwells cut out glass floors, so ignore them in this case
				float const z(f.z2());
				if (z <= floor_test_zval && z > closest_floor_zval) {closest_floor_zval = z; had_coll = 1;} // move up
			}
		}
		// check tunnels
		if (!interior->tunnels.empty() && point_in_extended_basement_not_basement(pos)) {
			point const p_test(pos.x, pos.y, obj_z);
			bool found(0);

			// first pass checks bcube, second pass checks bcube_ext (which may overlap other tunnels)
			for (unsigned pass = 0; pass < 2 && !found; ++pass) {
				for (tunnel_seg_t const &t : interior->tunnels) {
					if (!((pass == 0) ? t.bcube : t.bcube_ext).contains_pt(p_test)) continue;
					float const z(t.bcube.z1());
					if (z > floor_test_zval || z < closest_floor_zval) continue; // wrong zval
					tunnel_walk_area = t.get_walk_area(pos, xy_radius);
					tunnel_walk_area.clamp_pt_xy(pos);
					closest_floor_zval = z; // move up
					had_coll = player_in_tunnel = found = 1;
				}
			} // for pass
		}
		if (had_coll) {pos.z = closest_floor_zval + radius; obj_z = max(pos.z, p_last.z);}

		if (has_attic() && attic_access.contains_pt_xy(pos)) {
			if (interior->attic_access_open && obj_z > (attic_access.z2() - floor_spacing)) { // on attic ladder - handle like a ramp
				apply_speed_factor(pos, p_last, 0.3); // slow down to 30% when climbing the ladder
				bool const dim(attic_access.ix >> 1), dir(attic_access.ix & 1);
				float const length(attic_access.get_sz_dim(dim)), t(CLIP_TO_01((pos[dim] - attic_access.d[dim][0])/length)), T(dir ? t : (1.0-t));
				pos.z = (attic_access.z2() + attic_door_z_gap) - floor_spacing*(1.0 - T) + radius;
				obj_z = max(pos.z, p_last.z);
				had_coll = on_attic_ladder = 1;
			}
			else if (!interior->attic_access_open && obj_z > attic_access.z2()) { // standing on closed attic access door
				pos.z = attic_access.z2() + radius; // move up
				obj_z = max(pos.z, p_last.z);
				had_coll = 1;
			}
		}
		if (is_in_attic) {
			vector3d roof_normal;
			float const beam_depth(get_attic_beam_depth());

			// check player's head against the roof/overhead beams to avoid clipping through it
			if (!point_in_attic(point(pos.x, pos.y, (pos.z + 1.1f*camera_height + beam_depth)), &roof_normal)) {
				for (unsigned d = 0; d < 2; ++d) { // reset pos X/Y if oriented toward the roof
					if (roof_normal[d] != 0.0 && ((roof_normal[d] < 0.0) ^ (pos[d] < p_last[d]))) {reset_to_last_dims |= (1<<d);}
				}
				if (cnorm) {*cnorm = roof_normal;}

				if (camera_height > prev_camera_height) { // less crouched than before
					camera_height = prev_camera_height;
					force_player_height(prev_camera_height); // force crouch
				}
			}
			// find the part the player is in and clamp our bsphere to it; currently attics are limited to a single part
			for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
				if (!i->contains_pt_xy(pos)) continue; // wrong part
				cube_t valid_area(*i);
				valid_area.expand_by_xy(-xy_radius);
				valid_area.clamp_pt_xy(pos);
				break;
			} // for i
			player_in_attic = (has_attic_window ? 1 : 2);
			had_coll = 1;
			obj_z    = max(pos.z, p_last.z);
		}
	}
	// *** Note ***: at this point pos.z is radius above the floors and *not* at the camera height, which is why we add camera_height to pos.z below
	// Note: this check must be after pos.z is set from interior->floors
	// pass in radius as wall_test_extra_z as a hack to allow player to step over a wall that's below the stairs connecting stacked parts
	had_coll |= interior->check_sphere_coll_walls_elevators_doors(*this, pos, p_last, xy_radius, radius, 1, cnorm); // is_player=1

	if (interior->room_geom) { // collision with room geometry
		vect_room_object_t const &objs(interior->room_geom->objs);
		float speed_factor(1.0);

		for (auto c = interior->room_geom->get_stairs_start(); c != objs.end(); ++c) { // check for and handle stairs first
			if (c->no_coll() || c->type != TYPE_STAIR) continue;
			if (!c->contains_pt_xy(pos))  continue; // sphere not on this stair
			if (obj_z < c->z1())          continue; // below the stair
			if (pos.z - radius > c->z2()) continue; // above the stair
			bool const is_u(c->shape == SHAPE_STAIRS_U), is_l(c->shape == SHAPE_STAIRS_L), is_ul(is_u || is_l);
			// U-shaped stairs have strange collisions, and it seems to work better to have the player stand on the z1 value rather than the z2;
			// but stairs that are tall (such as pool stairs) only work when standing on the z2 value
			float const stairs_zval(is_u ? c->z1() : c->z2());
			pos.z = stairs_zval + radius; // stand on the stair - this can happen for multiple stairs

			if (has_basement()) { // don't let the zval end where the camera is between the terrain mesh and the basement entrance plane or we'll see the terrain
				float const camera_z(pos.z + camera_height), mesh_z(get_basement().z2()), entrance_z(mesh_z + BASEMENT_ENTRANCE_SCALE*floor_thickness);
				if (camera_z > mesh_z && camera_z < entrance_z) {pos.z += 1.2*(entrance_z - camera_z);} // move up above the entrance plane
			}
			obj_z = max(pos.z, p_last.z);

			if (c->shape != SHAPE_STAIRS_FAN) { // force the sphere onto the stairs unless it's a fan (mall entrance)
				if (!is_ul || c->dir == 1) {max_eq(pos[!c->dim], (c->d[!c->dim][0] + xy_radius));}
				if (!is_ul || c->dir == 0) {min_eq(pos[!c->dim], (c->d[!c->dim][1] - xy_radius));}
			}
			had_coll     = on_stairs = 1;
			speed_factor = (is_u ? 0.875 : 0.75); // U-shaped stairs are a bit faster
		} // for c
		if (speed_factor < 1.0) {apply_speed_factor(pos, p_last, speed_factor);} // slow down when on stairs
		point pos_high(pos); // capture before apply object collisions; used for diving board
		max_eq(pos_high.z, p_last.z);

		for (auto c = objs.begin(); c != objs.end(); ++c) { // check for other objects to collide with (including stairs)
			if (!c->is_player_collidable()) continue;
			if (on_attic_ladder && c->type == TYPE_ATTIC_DOOR) continue; // collision with attic door/ladder is handled above

			if (c->type == TYPE_ELEVATOR) { // special handling for elevators
				cube_t car_inc_door(*c);
				car_inc_door.d[c->dim][c->dir] += get_wall_thickness()*(c->dir ? 1.0 : -1.0); // add extra space in front of the elevator for the door
				if (!car_inc_door.contains_pt_xy(pos)) continue;
				elevator_t const &elevator(get_elevator(c->room_id));
				bool ecoll(0);
				if      (obj_z >= elevator.z2()) {} // above the elevator shaft - not in the elevator
				else if (obj_z >= c->z2()) {max_eq(pos.z, (c->z2() + radius)); ecoll = 1;} // standing on the roof of the elevator
				else if (obj_z >= c->z1()) {max_eq(pos.z, (c->z1() + radius + get_fc_thickness())); ecoll = 1;} // inside the elevator
				
				if (ecoll) {
					// 3 if moving, 2 if doors are closed, otherwise 1; used to avoid drawing terrain as we pass through it when the elevator enters or leaves the basement
					player_in_elevator     = ((elevator.open_amt == 0.0) ? (elevator.is_moving ? 3 : 2) : 1);
					player_in_int_elevator = elevator.interior_room;
					had_coll = 1;
				}
				continue;
			}
			if ((c->type == TYPE_STAIR || (on_stairs && c->type != TYPE_RAILING)) && (obj_z + radius) > c->z2()) continue; // above the stair - allow it to be walked on

			if (c->is_sloped_ramp()) { // ramp, should be SHAPE_ANGLED
				// slight adjust so that player is above ramp when on the floor above, but below when falling through the ramp gap
				bool const is_falling(pos.z < obj_z);
				float const reff((is_falling ? 1.01 : 0.99)*radius); // effective radius
				if ((obj_z - reff) > c->z2()) continue; // above the top of the ramp - not walking on it
				if (!sphere_cube_intersect_xy(pos, xy_radius, *c)) continue; // optimization - not actually on the ramp
				float &pos_dim(pos[c->dim]);
				if (pos_dim < c->d[c->dim][0] || pos_dim > c->d[c->dim][1]) continue; // not yet over ramp
				float const length(c->get_length()), height(c->dz()), t(CLIP_TO_01((pos_dim - c->d[c->dim][0])/length));
				float const ztop(c->z1() + height*(c->dir ? t : (1.0-t))), zbot(ztop - RAMP_THICKNESS_SCALE*height);
				float const player_height(camera_height + NEAR_CLIP); // include near clip for collisions with the bottom of ramps
				float const player_bot_z(obj_z - radius), player_bot_z_step(player_bot_z + C_STEP_HEIGHT*radius), player_top_z(obj_z + player_height);
					
				if (ztop < player_bot_z_step) { // step onto or move along ramp top surface
					pos.z = ztop + reff; // use reff rather than radius to prevent jitter when stepping onto the upper edge of the ramp
					obj_z = max(pos.z, p_last.z);
					had_coll = 1;

					if (!(c->flags & RO_FLAG_OPEN)) { // top floor ramp can't be walked on if blocked by a ceiling
						for (auto i = interior->ceilings.begin(); i != interior->ceilings.end(); ++i) {
							if (!i->contains_pt_xy(pos)) continue; // sphere not in this ceiling
							if (player_bot_z_step > i->z2() || player_top_z < i->z1()) continue; // no Z overlap
							float const dz(player_top_z - i->z1()), delta(dz*length/height);
							pos_dim += (c->dir ? -1.0 : 1.0)*delta; // push player back down the ramp
							break; // only change zval once
						}
					}
					continue;
				}
				else if (c->type == TYPE_RAMP && zbot < player_top_z) { // colliding with sides or bottom of the ramp (not pool)
					cube_t ramp_ext(*c);
					ramp_ext.expand_in_dim(c->dim, 1.01*xy_radius); // extend to include the player radius at both ends

					if (ramp_ext.contains_pt_xy(pos)) { // colliding with the bottom
						float const dz(player_top_z - zbot), delta(dz*length/height);
						pos_dim += (c->dir ? 1.0 : -1.0)*delta;
						had_coll = 1;
						continue;
					} // else treat it as a cube side collision
				}
				else {continue;} // under the ramp: no collision
			} // end ramp case
			cube_t c_extended(get_true_room_obj_bcube(*c));
			// handle the player's head for stairs and attic doors; only applied to bottom floor of stairs to avoid getting stuck on steep stairs that span multiple floors
			if ((c->type == TYPE_STAIR && (c->flags & RO_FLAG_ADJ_BOT)) || c->type == TYPE_ATTIC_DOOR) {c_extended.z1() -= camera_height;}

			if (c->type == TYPE_DIV_BOARD && c_extended.contains_pt_xy(pos)) {
				// diving boards are special because the player can walk on them across different floor heights
				vector3d coll_norm;

				if (sphere_cube_int_update_pos(pos_high, xy_radius, c_extended, p_last, 0, &coll_norm)) {
					had_coll = 1;

					if (coll_norm.z == 1.0) { // update pos if on the top surface
						pos = pos_high;
						if (cnorm) {*cnorm = coll_norm;}
					}
					continue;
				}
			}
			if (c->type == TYPE_BAN_PEEL) {
				c_extended.z2() += radius; // increase height
				if (!sphere_cube_intersect(pos, xy_radius, c_extended)) continue;
				apply_speed_factor(pos, p_last, 2.0); // player slips; doesn't have much of an effect
				continue;
			}
			if (!sphere_cube_intersect(pos, xy_radius, c_extended)) continue; // optimization

			if (c->type == TYPE_RAILING && !(c->flags & RO_FLAG_IN_POOL)) { // stairs railing, not on a pool
				// only collide with railing at top of stairs, not when walking under stairs
				cylinder_3dw const railing(get_railing_cylinder(*c));
				float const t((pos[c->dim] - railing.p1[c->dim])/(railing.p2[c->dim] - railing.p1[c->dim]));
				float const railing_zval(railing.p1.z + CLIP_TO_01(t)*(railing.p2.z - railing.p1.z));
				if ((railing_zval - get_railing_height(*c)) > float(pos.z + camera_height) || railing_zval < (pos.z - radius)) continue; // no Z collision
			}
			if (c->type == TYPE_INT_LADDER && c->dz() > (camera_height + radius)) { // vertical ladder attached to or leaning against a wall
				if (c->in_factory() && (pos[c->dim] < c->get_center_dim(c->dim)) == c->dir) continue; // wrong side of vertical factory ladder
				pos.z = p_last.z + 0.25*get_player_move_dist()*cview_dir.z; // move up/down based on player vertical view (looking up vs. down)
				pos.z = min(float(c->z2() - camera_height), max(float(c->z1() + radius), pos.z)); // clamp to ladder height range
				obj_z = max(pos.z, p_last.z);
				had_coll = 1;
			}
			// Note: only vert pipes have player coll; ducts are not vert and are treated as cubes
			else if (c->is_vert_cylinder()) { // vertical cylinder
				cylinder_3dw cylin(c->get_cylinder());
				cylin.p2.z += radius; // extend upward by radius
				had_coll |= sphere_vert_cylin_intersect_with_ends(pos, xy_radius, cylin, cnorm);
			}
			else if (c->shape == SHAPE_SPHERE) { // sphere
				point const center(c->get_cube_center());
				if (!dist_less_than(pos, center, (c->get_radius() + radius))) continue; // xy_radius?
				if (cnorm) {*cnorm = (pos - center).get_norm();}
				had_coll = 1;
			}
			else if (c->type == TYPE_CLOSET) { // special case to handle closet interiors
				had_coll |= (bool)check_closet_collision(*c, pos, p_last, xy_radius, cnorm);
				
				if (c->contains_pt(pos)) {
					player_in_closet |= RO_FLAG_IN_CLOSET;
					if (interior->room_geom->closet_light_is_on(*c)) {player_in_closet |= RO_FLAG_LIT;}
					if (c->is_open()) {player_in_closet |= RO_FLAG_OPEN;} else {register_player_hiding(*c);} // player is hiding if the closet door is closed
				}
			}
			else if (c->type == TYPE_STALL && maybe_inside_room_object(*c, pos, xy_radius)) {
				// stall is open and intersecting player, or player is inside stall; perform collision test with sides only
				had_coll |= check_stall_collision(*c, pos, p_last, xy_radius, cnorm);
				if (!c->is_open() && c->contains_pt(pos)) {register_in_closed_bathroom_stall();}
			}
			else if (c->type == TYPE_SHOWER && maybe_inside_room_object(*c, pos, xy_radius)) {
				// shower is open and intersecting player, or player is inside shower; perform collision test with side only
				had_coll |= check_shower_collision(*c, pos, p_last, xy_radius, cnorm);
			}
			else if (c->type == TYPE_BALCONY) {
				had_coll |= check_balcony_collision(*c, pos, p_last, xy_radius, cnorm);
			}
			else if (c->type == TYPE_RDESK) {
				had_coll |= (check_rdesk_collision(*c, pos, p_last, xy_radius, cnorm) != 0);
			}
			else if (sphere_cube_int_update_pos(pos, xy_radius, c_extended, p_last, 0, cnorm)) { // assume it's a cube; skip_z=0
				if (c->type == TYPE_TOILET || (c->type == TYPE_URINAL && !is_player_model_female())) {player_near_toilet = 1;} // females can't use urinals
				had_coll = 1;
			}
			if ((c->type == TYPE_STALL || c->type == TYPE_SHOWER) && !c->is_open() && c->contains_pt(pos)) {register_player_hiding(*c);} // player is hiding in the stall/shower
		} // for c
	} // end interior->room_geom
	if (!player_wait_respawn) { // zombies push player, but not when player is dead
		for (person_t const &p : interior->people) {
			float const dist(p2p_dist(pos, p.pos)), r_sum(xy_radius + 0.5*p.get_width()); // ped_radius is a bit too large, multiply it by 0.5
			if (dist >= r_sum) continue; // no intersection
			vector3d const normal(vector3d(pos.x-p.pos.x, pos.y-p.pos.y, 0.0).get_norm()); // XY direction
			if (cnorm) {*cnorm = normal;}
			pos += normal*(r_sum - dist);
			if (!tunnel_walk_area.is_all_zeros()) {tunnel_walk_area.clamp_pt_xy(pos);} // re-clamp to avoid zombie pushing the player through a gate
			had_coll = 1;
		} // for i
	}
	for (unsigned d = 0; d < 2; ++d) { // apply attic roof pos reset at the end, to override other collisions that may move pos
		if (reset_to_last_dims & (1<<d)) {pos[d] = p_last[d];}
	}
	float const z_feet(pos.z - radius), z_head(pos.z + camera_height);
	handle_vert_cylin_tape_collision(pos, p_last, z_feet, z_head, xy_radius, 1); // is_player=1
	
	if (has_water() && get_water_cube().contains_pt_xy(pos) && z_feet < interior->water_zval) { // feet in the water
		float const submerge_amt(min(1.0f, ((interior->water_zval - z_feet) / (z_head - z_feet))));
		// player can hide in a pool if head is underwater, but this doesn't affect zombies because they don't enter the pool anyway
		player_is_hiding |= (has_pool() && submerge_amt == 1.0);
		apply_speed_factor(pos, p_last, (1.0 - 0.4*submerge_amt)); // up to 2x slower when more submerged
	}
	// not sure where this belongs, but the closet hiding logic is in this function, so I guess it goes here? player must be inside the building to see a windowless room anyway;
	// assume room is lit if the player is on the stairs; required for handling stairs cut into walls where the stairs pass through a room not marked as having stairs
	player_in_unlit_room   = (!on_stairs && check_pos_in_unlit_room(pos));
	player_on_house_stairs = (on_stairs && is_house);
	prev_camera_height     = camera_height; // update for this frame
	if (had_coll && !player_in_water && pos.z < p_last.z) {apply_building_fall_damage(p_last.z - pos.z);}
	return had_coll; // will generally always be true due to floors
}

bool building_t::check_cube_intersect_non_main_part(cube_t const &c) const {
	if (!c.intersects(bcube)) return 0;
	// check parts such as garage, shed, porch roof, porch support, and chimney for intersections
	for (auto p = get_real_parts_end(); p != parts.end(); ++p) {
		if (p->intersects(c)) return 1;
	}
	return 0;
}

bool building_t::point_in_elevator(point const &pos, bool check_elevator_car) const {
	if (!interior || (check_elevator_car && !has_room_geom())) return 0;

	for (elevator_t const &e : interior->elevators) {
		if (!e.contains_pt(pos)) continue;
		if (check_elevator_car && !interior->get_elevator_car(e).contains_pt(pos)) continue;
		return 1;
	}
	return 0;
}
bool building_interior_t::point_in_U_stairwell(point const &pos, float floor_spacing, bool mall_only) const {
	for (stairwell_t const &s : stairwells) {
		if (!s.is_u_shape())         continue;
		if (mall_only && !s.in_mall) continue;
		if (s.contains_pt(pos))      return 1;
		if (!s.not_an_exit_mask)     continue;
		// handle landing
		int const floor_ix(max(0, int((pos.z - s.z1())/floor_spacing)));
			
		if (s.not_an_exit_mask & (1 << floor_ix)) {
			cube_t stairs_ext(s);
			stairs_ext.d[s.dim][!s.dir] += (s.dir ? -1.0 : 1.0)*s.get_retail_landing_width(floor_spacing);
			if (stairs_ext.contains_pt(pos)) return 1;
		}
	} // for s
	return 0;
}

bool building_t::check_pos_in_unlit_room(point const &pos) const {
	if (!interior) return 0; // error?

	if (!interior->elevators_disabled) { // check if in elevator, if powered, since light is always on; assumes point is in the elevator car
		for (elevator_t const &e : interior->elevators) {
			if (e.get_bcube_padded(get_wall_thickness()).contains_pt(pos)) return 0; // include the wall in front
		}
	}
	set<unsigned> rooms_visited;
	return check_pos_in_unlit_room_recur(pos, rooms_visited);
}
bool building_t::check_pos_in_unlit_room_recur(point const &pos, set<unsigned> &rooms_visited, int known_room_id) const {
	int const room_id((known_room_id >= 0) ? known_room_id : get_room_containing_pt(pos));
	if (room_id < 0) return 0; // not in a room
	if (rooms_visited.find(room_id) != rooms_visited.end()) return 1; // already visited, return true
	room_t const &room(get_room(room_id));
	if (!room.interior && has_int_windows()) return 0; // room has windows and may be lit from outside; what about a room with an open exterior door?
	// Note: we can't easily check the connected building here, so assume it's lit; this doesn't fix the case where the conn room is in the other building
	if ( room.is_ext_basement_conn       ()) return 0; // be safe/conservative and treat rooms connected to two buildings as potentially lit
	if ( room.has_interior_window        ()) return 0; // light can go through the window
	float const floor_spacing(get_window_vspace());
	if (room.is_backrooms() && room.dz() > 1.5*floor_spacing) return 0; // multi-floor backrooms: assume lit as this is too difficult to determine
	unsigned const floor_ix(room.is_single_floor ? 0 : max(0.0f, (pos.z - room.z1()))/floor_spacing);
	if (room.get_has_skylight() && pos.z > (room.z2() - floor_spacing)) return 0; // top floor of room with a skylight
	if (room.has_elevator || (!room.is_ext_basement() && room.has_stairs_on_floor(floor_ix))) return 0; // assume light can come from stairs (not in ext basement) or open elevator

	if (pos.z > ground_floor_z1 && (pos.z < ground_floor_z1 + floor_spacing || !walkways.empty())) { // on the ground floor or have walkways
		if (is_room_adjacent_to_ext_door(room) && point_near_ext_door(pos, get_door_open_dist())) return 0; // handle windowless office building open exterior doors
	}
	// check if all lights are off
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if ((int)i->room_id != room_id)               continue; // wrong room
		if (!i->is_light_type() || !i->is_light_on()) continue; // not a light, or light not on
		//if (i->light_is_out()) continue; // broken light; continuing here doesn't work because sparks won't be drawn
		if (!room.is_single_floor && unsigned(max(0.0f, (i->z1() - room.z1()))/floor_spacing) != floor_ix) continue; // different floors
		return 0; // lit by a room light, including one in a closet (Note that closets are only in house bedrooms, which should always have windows anyway)
	}
	rooms_visited.insert(room_id); // mark this room as visited before making recursive calls

	// check for a light path through a series of open doors
	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
		if (i->get_bldg_conn()) {
			door_t const &door(interior->doors[i->first_door_ix]); // should be only one door
			if (door.open_amt > 0.0 && door.get_true_bcube().intersects(room)) return 0; // open ext basement conn door counts, assuming light is on or opposite door is open
			continue;
		}
		if ( i->not_a_room_separator())        continue; // skip non room separating doors
		if (!i->is_connected_to_room(room_id)) continue;
		assert(i->first_door_ix < interior->doors.size());

		for (unsigned dix = i->first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t const &door(interior->doors[dix]);
			if (!i->is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() > pos.z || door.z2() < pos.z) continue; // wrong floor
			if (door.open_amt == 0.0) continue; // fully closed
			// check adjacent room; since we're passing in room_id, we can reuse pos; Note: usually windowless rooms are connected to office building hallways
			if (!check_pos_in_unlit_room_recur(pos, rooms_visited, i->get_conn_room(room_id))) return 0; // if adjacent room is lit, return false
		} // for dix
	} // for i
	if (room.has_stairs && room.is_ext_basement()) {
		// extended basement room - check for stairs (should all be straight)
		for (stairwell_t const &s : interior->stairwells) {
			if (!s.intersects(room)) continue;
			bool room_is_above(room.z2() > s.z2());
			point const pos2(s.xc(), s.yc(), (room.zc() + (room_is_above ? -1.0 : 1.0)*floor_spacing));
			if (!check_pos_in_unlit_room_recur(pos2, rooms_visited)) return 0; // if adjacent room is lit, return false
		}
	}
	// check for light through connected office building hallways, open wall rooms, or adjacent/wall-less complex floorplan building rooms
	if (room.open_wall_mask || (has_complex_floorplan && !room.is_ext_basement())) {
		float const room_exp(0.5*get_wall_thickness());
		cube_t test_cube(room); // include adjacency, but don't expand enough to go through a wall

		for (unsigned d = 0; d < 2; ++d) {
			if (room.has_open_wall(d, 0)) {test_cube.d[d][0] -= room_exp;} // lo edge
			if (room.has_open_wall(d, 1)) {test_cube.d[d][1] += room_exp;} // hi edge
		}
		set_cube_zvals(test_cube, pos.z-0.1*floor_spacing, pos.z+floor_spacing); // clip to narrow Z-range
		unsigned const rooms_end((interior->ext_basement_hallway_room_id >= 0) ? interior->ext_basement_hallway_room_id : interior->rooms.size());

		for (unsigned r = 0; r < rooms_end; ++r) {
			if ((int)r == room_id) continue; // same room, skip
			room_t const &room2(interior->rooms[r]);
			if (!(room2.open_wall_mask || has_complex_floorplan) || !room2.intersects_no_adj(test_cube)) continue;
			point const center(room2.get_cube_center());
			if (!check_pos_in_unlit_room_recur(point(center.x, center.y, pos.z), rooms_visited, r)) return 0; // if adjacent room is lit, return false; room_id is known
		}
	}
	return 1;
}

bool building_t::is_room_windowless(room_t const &room) const { // or interior
	if (room.z2() <= ground_floor_z1)         return 1; // basement, no windows
	if (!bcube.contains_cube_xy_no_adj(room)) return 0; // adjacent to bcube, may have windows
	cube_t const &part(get_part_for_room(room));
	if (part.contains_cube_xy_no_adj(room))   return 1; // interior to its part, no windows
	return (room.ext_sides == 0); // maybe we should check this just after the basement test?
}

bool building_t::are_rooms_connected(room_t const &r1, room_t const &r2, float zval, bool check_door_open) const {
	float const expand(2.0*get_wall_thickness()); // expand so that adjacent rooms and doors overlap
	cube_t tc1(r1), tc2(r2);
	tc1.expand_by_xy(expand);
	if (!tc1.intersects(tc2)) return 0; // rooms not adjacent

	if (r1.open_wall_mask && r2.open_wall_mask) { // check open walls
		for (cube_t const &w : interior->open_walls) {
			if (w.intersects(r1) && w.intersects(r2)) return 1;
		}
	}
	tc2.expand_by_xy(expand);
	tc1.intersect_with_cube(tc2); // shared wall area

	for (auto const &ds : interior->door_stacks) {
		if (ds.not_a_room_separator()) continue;
		if (!tc1.contains_cube(ds))    continue;
		if (!check_door_open) return 1; // if there's a door stack, there must be a connecting door
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t const &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() < zval && door.z2() > zval && door.open_amt > 0.0 && tc1.contains_cube(door)) return 1;
		}
	} // for ds
	return 0;
}

bool building_t::all_room_int_doors_closed(unsigned room_ix, float zval) const {
	if (has_complex_floorplan) return 0; // not supported, as there may be missing walls
	room_t const &room(get_room(room_ix));
	if (room.is_sec_bldg)                         return 0; // garages and sheds don't count because objects can be seen through the windows
	if (room.open_wall_mask)                      return 0; // office hallways and open wall rooms can connect to other hallways with no doors
	if (room.is_parking() || room.is_backrooms()) return 0; // these cases are excluded because they have interior doors or ramps
	if (room.is_retail())                         return 0; // retail doesn't work because objects may be visible through stairs (similar to office hallway)
	if (room.is_mall_or_store())                  return 0; // malls and stores have gated doorways rather than doors that can be seen through
	if (room.has_interior_window())               return 0; // can see through the window (conference rooms, mall stores)
	min_eq(zval, room.z2()); max_eq(zval, room.z1()); // clamp to room bounds - is this needed?
	if (room.has_stairs_on_floor(room.get_floor_containing_zval(zval, get_window_vspace()))) return 0; // might be visible through stairs
	// if the room is a single floor, check doors within the full Z span; otherwise, only consider doors overlapping zval
	float const z1(room.is_single_floor ? room.z1() : zval), z2(room.is_single_floor ? room.z2() : zval);

	for (auto const &ds : interior->door_stacks) {
		if (!ds.is_connected_to_room(room_ix)) continue;
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t const &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() < z2 && door.z2() > z1 && door.open_amt > 0.0) return 0;
		}
	} // for ds
	return 1;
}

// Note: called on balls and particles
bool building_interior_t::check_sphere_coll(building_t const &building, point &pos, point const &p_last, float radius,
	vect_room_object_t::const_iterator self, vector3d &cnorm, float &hardness, int &obj_ix, bool is_ball) const
{
	bool had_coll(check_sphere_coll_walls_elevators_doors(building, pos, p_last, radius, 0.0, 0, &cnorm)); // is_player=0
	if (had_coll) {hardness = 1.0;}
	cube_t scube; scube.set_from_sphere(pos, radius);

	if (pos.z < p_last.z) { // check floor collision if falling
		max_eq(scube.z2(), (p_last.z - radius)); // extend up to touch the bottom of the last position to prevent it from going through the floor in one frame

		for (auto i = floors.begin(); i != floors.end(); ++i) {
			if (!i->intersects(scube) || !i->contains_pt_xy(pos)) continue; // no overlap
			pos.z = i->z2() + radius; // move to just touch the top of the floor
			cnorm = plus_z; // collision with top surface of floor
			had_coll = 1; hardness = 1.0;
		}
	}
	else { // check ceiling collision if rising
		min_eq(scube.z1(), (p_last.z + radius)); // extend down

		for (auto i = ceilings.begin(); i != ceilings.end(); ++i) {
			if (!i->intersects(scube)) continue; // overlap
			pos.z = i->z1() - radius; // move to just touch the top of the ceiling
			cnorm = -plus_z; // collision with top surface of ceiling
			had_coll = 1; hardness = 1.0;
		}
	}
	if (room_geom) { // check glass floors
		for (cube_t const &c : room_geom->glass_floors) {
			if (sphere_cube_int_update_pos(pos, radius, c, p_last, 0, &cnorm)) {had_coll = 1; hardness = 1.0;} // skip_z=0
		}
	}
	had_coll |= check_sphere_coll_room_objects(building, pos, p_last, radius, self, cnorm, hardness, obj_ix, is_ball);
	return had_coll;
}
bool building_interior_t::check_sphere_coll_room_objects(building_t const &building, point &pos, point const &p_last, float radius,
	vect_room_object_t::const_iterator self, vector3d &cnorm, float &hardness, int &obj_ix, bool is_ball) const
{
	if (!room_geom) return 0;
	bool had_coll(0);

	// Note: no collision check with expanded_objs
	for (auto c = room_geom->objs.begin(); c != room_geom->objs.end(); ++c) { // check for other objects to collide with
		room_object const type(c->type);
		// ignore blockers and railings, but allow more than c->no_coll()
		if (c == self || type == TYPE_BLOCKER || type == TYPE_PAPER || type == TYPE_PEN || type == TYPE_PENCIL || c->is_a_drink() || type == TYPE_FLOORING || type == TYPE_SIGN ||
			type == TYPE_WBOARD || type == TYPE_WALL_TRIM || type == TYPE_DRAIN || type == TYPE_CRACK || type == TYPE_SWITCH || type == TYPE_BREAKER || type == TYPE_OUTLET ||
			type == TYPE_VENT || type == TYPE_WIND_SILL || type == TYPE_TEESHIRT || type == TYPE_PANTS || type == TYPE_BLANKET || type == TYPE_FOLD_SHIRT ||
			type == TYPE_SPIWEB || type == TYPE_ELEC_WIRE) continue;
		if (type == TYPE_RAILING && (!(c->flags & RO_FLAG_TOS) || !c->is_open())) continue; // only railings at the top of stairs (non-sloped) with balusters have collisions
		if (type == TYPE_POOL_TILE && c->no_coll()) continue;
		cube_t const bc(get_true_room_obj_bcube(*c));
		if (!sphere_cube_intersect(pos, radius, bc)) continue; // no intersection (optimization)
		if (is_ball && type == TYPE_BOOK && c->dz() < 0.2*radius) continue; // large ball vs. small book on the floor: unstable, skip
		unsigned coll_ret(0);
		// add special handling for things like elevators and cubicles? right now these are only in office buildings, where there are no dynamic objects

		if (c->is_vert_cylinder()) { // vertical cylinder (including table)
			cylinder_3dw const cylin(c->get_cylinder());

			if (type == TYPE_TABLE) {
				cylinder_3dw top(cylin), base(cylin);
				top.p1.z  = base.p2.z = c->z2() - 0.12*c->dz(); // top shifted down by 0.12
				base.r1  *= 0.4; base.r2 *= 0.4; // vertical support has radius 0.08, legs have radius of 0.6, so use something in between
				coll_ret |= unsigned(sphere_vert_cylin_intersect_with_ends(pos, radius, top, &cnorm) || sphere_vert_cylin_intersect_with_ends(pos, radius, base, &cnorm));
			}
			else {coll_ret |= (unsigned)sphere_vert_cylin_intersect_with_ends(pos, radius, cylin, &cnorm);}
		}
		else if (c->shape == SHAPE_SPHERE) { // sphere
			point const center(c->get_cube_center());
			if (!dist_less_than(pos, center, (c->get_radius() + radius))) continue;
			cnorm = (pos - center).get_norm();
			coll_ret |= 1;
		}
		else if (c->shape == SHAPE_VERT_TORUS) { // pool float, etc.; treated as a vertical ring/cylinder with an inner hole
			float const ri(0.5*c->dz()), ro(c->get_radius() - ri);
			point p_int;

			if (sphere_torus_intersect(pos, radius, c->get_cube_center(), ri, ro, p_int, cnorm, 1)) { // calc_int=1
				pos = p_int; // move to not intersect
				coll_ret |= 1;
			}
		}
		else if (type == TYPE_SHELVES && c->in_mall()) { // mall shelves have fishtanks (pet stores) and flat items (clothing stores) that we can collide with
			cube_t shelves[4]; // max number of shelves
			unsigned const num_shelves(get_shelves_for_object(*c, shelves));
			coll_ret |= check_cubes_collision(shelves, num_shelves, pos, p_last, radius, &cnorm);
		}
		else { // assume it's a cube
			// some object types are special because they're common collision objects and they're not filled cubes
			if      (type == TYPE_CLOSET    ) {coll_ret |= check_closet_collision(*c, pos, p_last, radius, &cnorm);} // special case to handle closet interiors
			else if (type == TYPE_BED       ) {coll_ret |= check_bed_collision   (*c, pos, p_last, radius, &cnorm);}
			else if (can_use_table_coll(*c) ) {coll_ret |= check_table_collision (*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_CONF_TABLE) {coll_ret |= check_conf_table_collision(*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_RDESK     ) {coll_ret |= check_rdesk_collision (*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_CHAIR     ) {coll_ret |= check_chair_collision (*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_TUB       ) {coll_ret |= check_tub_collision   (*c, pos, p_last, radius, &cnorm, is_ball);}
			else if (type == TYPE_BCASE     ) {coll_ret |= check_bookcase_collision(*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_BENCH     ) {coll_ret |= check_bench_collision   (*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_DIV_BOARD ) {coll_ret |= check_diving_board_collision     (*c, pos, p_last, radius, &cnorm);}
			else if (c->is_sloped_ramp()    ) {coll_ret |= (unsigned)check_ramp_collision   (*c, pos,         radius, &cnorm);} // p_last is unused
			else if (type == TYPE_BALCONY   ) {had_coll |= (unsigned)check_balcony_collision(*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_POOL_TABLE) {coll_ret |= check_pool_table_collision(*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_SHELFRACK ) {coll_ret |= check_shelf_rack_collision(*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_COUCH     ) {coll_ret |= check_couch_collision     (*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_CASHREG   ) {coll_ret |= check_cashreg_collision   (*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_STALL  && maybe_inside_room_object(*c, pos, radius)) {coll_ret |= (unsigned)check_stall_collision (*c, pos, p_last, radius, &cnorm);}
			else if (type == TYPE_SHOWER && maybe_inside_room_object(*c, pos, radius)) {coll_ret |= (unsigned)check_shower_collision(*c, pos, p_last, radius, &cnorm);}
			else {coll_ret |= (unsigned)sphere_cube_int_update_pos(pos, radius, bc, p_last, 0, &cnorm);} // skip_z=0
		}
		if (coll_ret) { // collision with this object - set hardness
			if      (type == TYPE_RUG && cnorm != plus_z) {cnorm = plus_z; continue;} // rug collision can only be +z
			if      (type == TYPE_COUCH    ) {hardness = 0.6;} // couches are soft
			else if (type == TYPE_RUG      ) {hardness = 0.8;} // rug is somewhat soft, but placed on a hard floor
			else if (type == TYPE_BLINDS   ) {hardness = 0.6;} // blinds are soft
			else if (type == TYPE_PIZZA_BOX) {hardness = 0.7;} // somewhat soft
			else if (type == TYPE_TEESHIRT || type == TYPE_PANTS || type == TYPE_FOLD_SHIRT) {hardness = 0.7;} // somewhat soft, though collisions may be skipped
			else if (type == TYPE_BLANKET  ) {hardness = 0.5;} // soft
			else if (type == TYPE_BED && (coll_ret & 24)) {hardness = 0.5;} // pillow/mattress collision is very soft
			else if (type == TYPE_TUB && cnorm == plus_z) {hardness = 0.2;} // bottom of tub, possibly with water; will make a ball less likely to bounce out
			else {hardness = 1.0;}
			obj_ix   = (c - room_geom->objs.begin()); // may be overwritten, will be the last collided object
			had_coll = 1;
		}
	} // for c
	return had_coll;
}

room_object_t const &building_interior_t::get_elevator_car(elevator_t const &e) const {
	assert(room_geom);
	assert(e.car_obj_id < room_geom->objs.size());
	room_object_t const &obj(room_geom->objs[e.car_obj_id]); // elevator car for this elevator
	assert(obj.type == TYPE_ELEVATOR);
	return obj;
}

bool sphere_cube_int_zval(point const &pos, float zval, float radius, cube_t const &c) {
	return (zval > c.z1() && zval < c.z2() && sphere_cube_intersect_xy(pos, radius, c));
}
bool sphere_cube_int_update_pos_zval(point &pos, float zval, float radius, cube_t const &c, point const &p_last, vector3d *cnorm) {
	return (zval > c.z1() && zval < c.z2() && sphere_cube_int_update_pos(pos, radius, c, p_last, 1, cnorm)); // skip_z=1
}
bool check_door_coll(building_t const &building, door_t const &door, point &pos, point const &p_last, float radius, float obj_z, bool check_open_doors, vector3d *cnorm) {
	if (door.open) {
		if (!check_open_doors) return 0; // doors tend to block the player and other objects, don't collide with them unless they're closed
		cube_t door_bounds(door);
		door_bounds.expand_by_xy(door.get_width());
		if (!sphere_cube_intersect(pos, radius, door_bounds)) return 0; // check intersection with rough/conservative door bounds (optimization)
		tquad_with_ix_t const door_tq(building.set_interior_door_from_cube(door));
		vector3d normal(door_tq.get_norm());
		if (dot_product_ptv(normal, pos, door_tq.pts[0]) < 0.0) {normal.negate();} // use correct normal sign
		float rdist, thick;

		if (sphere_ext_poly_int_base(door_tq.pts[0], normal, pos, radius, door.get_thickness(), thick, rdist)) {
			if (sphere_poly_intersect(door_tq.pts, door_tq.npts, pos, normal, rdist, thick)) {
				pos += normal*(thick - rdist);
				if (cnorm) {*cnorm = normal;}
				return 1;
			}
		}
		return 0;
	}
	return sphere_cube_int_update_pos_zval(pos, obj_z, radius, door, p_last, cnorm);
}

cube_t escalator_t::get_ramp_bcube(bool exclude_sides) const {
	cube_t ramp(*this);
	if (exclude_sides) {ramp.expand_in_dim(!dim, -get_side_width());} // subtract off sides to get the belt/stairs
	cube_t lo_end(ramp), hi_end(ramp);
	ramp.expand_in_dim(dim, -end_ext); // subtract off ends
	ramp.z2() = z1() + delta_z;
	assert(ramp.is_strictly_normalized());
	return ramp;
}
void escalator_t::get_ends_bcube(cube_t &lo_end, cube_t &hi_end, bool exclude_sides) const {
	cube_t const ramp(get_ramp_bcube(exclude_sides));
	lo_end = hi_end = *this;
	for (unsigned d = 0; d < 2; ++d) {lo_end.d[!dim][d] = hi_end.d[!dim][d] = ramp.d[!dim][d];}
	float const side_height(get_side_height());
	lo_end.z2() = z1() + side_height;
	hi_end.z1() = ramp.z2();
	hi_end.z2() = hi_end.z1() + side_height;
	lo_end.d[dim][ dir] = ramp.d[dim][!dir];
	hi_end.d[dim][!dir] = ramp.d[dim][ dir];
}
cube_t escalator_t::get_side_for_end(cube_t const &end, bool lr) const {
	cube_t end_side(end);
	end_side.d[!dim][!lr] = end.d[!dim][lr] + (lr ? -1.0 : 1.0)*get_side_width();
	return end_side;
}
cube_t escalator_t::get_support_pillar() const {
	float const support_radius(0.15*get_width());
	cube_t support(*this);
	support.z2() = z1() + delta_z - get_upper_hang();
	set_wall_width(support, get_center_dim(!dim), support_radius, !dim);
	set_wall_width(support, (d[dim][dir] + (dir ? -1.0 : 1.0)*0.5*end_ext), support_radius, dim); // center of upper end
	return support;
}
void escalator_t::get_ramp_bottom_pts(cube_t const &ramp, point bot_pts[4]) const { // {lo-left, lo-right, hi-right, hi-left} (or reversed)
	bot_pts[0].z = bot_pts[1].z = ramp.z1() - bot_edge_shift; // shift below the steps for extra space and to hide the walking feet of people
	bot_pts[2].z = bot_pts[3].z = ramp.z2() - bot_edge_shift;
	bot_pts[0][ dim] = bot_pts[1][ dim] = ramp.d[dim][!dir];
	bot_pts[2][ dim] = bot_pts[3][ dim] = ramp.d[dim][ dir];
	bot_pts[0][!dim] = bot_pts[3][!dim] = ramp.d[!dim][0];
	bot_pts[1][!dim] = bot_pts[2][!dim] = ramp.d[!dim][1];
	if (dim ^ dir ^ 1) {std::reverse(bot_pts, bot_pts+4);} // use correct vertex winding order
}
void escalator_t::get_all_cubes(cube_t cubes[7]) const { // {lo left wall, lo right wall, lo floor, hi left wall, hi right wall, hi floor, pillar}
	cube_t ends[2];
	get_ends_bcube(ends[0], ends[1], 0); // exclude_sides=0
	unsigned ix(0);

	for (unsigned hi = 0; hi < 2; ++hi) {
		cube_t end(ends[hi]);
		for (unsigned lr = 0; lr < 2; ++lr) {cubes[ix++] = get_side_for_end(end, lr);}
		end.z2() = end.z1() + get_floor_thick();
		cubes[ix++] = end;
	}
	cubes[ix++] = get_support_pillar();
}

// Note: should be valid for player and other spherical objects; also checks escalators
bool building_interior_t::check_sphere_coll_walls_elevators_doors(building_t const &building, point &pos, point const &p_last, float radius,
	float wall_test_extra_z, bool is_player, vector3d *cnorm) const
{
	bool const check_open_doors(!is_player); // check_open_doors=0 for player to avoid getting the player stuck
	float obj_z(max(pos.z, p_last.z)); // use p_last to get orig zval
	float const wall_test_z(obj_z + wall_test_extra_z);
	bool had_coll(0);

	// Note: pos.z may be too small here and we should really use obj_z, so skip_z must be set to 1 in cube tests and obj_z tested explicitly instead
	for (unsigned d = 0; d < 2; ++d) { // check XY collision with walls
		for (auto i = walls[d].begin(); i != walls[d].end(); ++i) {
			if (wall_test_z < i->z1() || wall_test_z > i->z2()) continue; // wrong part/floor
			had_coll |= sphere_cube_int_update_pos(pos, radius, *i, p_last, 1, cnorm); // skip_z=1 (handled by zval test above)
		}
	}
	for (elevator_t const &e : elevators) {
		if (obj_z < e.z1() || obj_z > e.z2()) continue; // wrong part/floor

		if (room_geom && (e.open_amt > 0.75 || e.contains_pt(point(pos.x, pos.y, obj_z)))) { // elevator is mostly open, can enter || already inside elevator
			room_object_t const &obj(get_elevator_car(e)); // elevator car for this elevator

			if (obj_z > obj.z1() && obj_z < obj.z2()) { // same floor as elevator car - can enter it; otherwise can't enter elevator shaft
				cube_t cubes[5];
				unsigned const num_cubes(e.get_coll_cubes(cubes));
				for (unsigned n = 0; n < num_cubes; ++n) {had_coll |= sphere_cube_int_update_pos(pos, radius, cubes[n], p_last, 1, cnorm);} // skip_z=1
				continue; // done with elevator
			}
		}
		had_coll |= sphere_cube_int_update_pos(pos, radius, e, p_last, 1, cnorm); // handle as a single blocking cube; skip_z=1
	} // for e
	for (escalator_t const &e : escalators) {
		if (!sphere_cube_int_zval(pos, obj_z, radius, e)) continue; // Note: e.z2() extends up to the floor above, so this is conservative
		cube_t const ramp(e.get_ramp_bcube(0)); // exclude_sides=0
		cube_t ends[2];
		e.get_ends_bcube(ends[0], ends[1], 1); // exclude_sides=1
		if (obj_z > max(ramp.z2(), max(ends[0].z2(), ends[1].z2()))) continue; // above ramp and both ends
		bool handled(0);

		for (unsigned hi = 0; hi < 2; ++hi) { // check for standing on one of the ends
			if (!ends[hi].contains_pt(point(pos.x, pos.y, obj_z))) continue;
			float const floor_zmin(ends[hi].z1() + e.get_floor_thick() + radius);
			if (pos.z < floor_zmin) {pos.z = floor_zmin; max_eq(obj_z, pos.z); had_coll = 1;}
			if (is_player) {player_on_escalator = 1;}
			handled = 1;
		}
		if (!handled && sphere_cube_intersect_xy(pos, radius, ramp)) { // if entering or exiting, ignore the ramp
			point bot_pts[4]; // {lo-left, lo-right, hi-right, hi-left}
			e.get_ramp_bottom_pts(ramp, bot_pts);
			vector3d const normal(get_poly_norm(bot_pts));

			if (is_player) {
				cube_t const ramp_inner(e.get_ramp_bcube(1)); // exclude_sides=1

				if (ramp_inner.contains_pt_xy(pos)) { // player on escalator steps
					float const length(ramp_inner.get_sz_dim(e.dim)), height(ramp_inner.dz()), t(CLIP_TO_01((pos[e.dim] - ramp_inner.d[e.dim][0])/length));
					float const ztop(ramp_inner.z1() + height*(e.dir ? t : (1.0-t)));

					if (ztop < obj_z - radius + C_STEP_HEIGHT*radius) { // step onto or move along ramp top surface
						bool const is_falling(pos.z < obj_z);
						float const reff((is_falling ? 1.01 : 0.99)*radius); // effective radius
						pos.z = ztop + reff; // use reff rather than radius to prevent jitter when stepping onto the upper edge of the ramp
						apply_speed_factor(pos, p_last, 0.5); // slower when walking on escalator

						if (animate2) { // apply escalator movement
							float const delta((e.move_dir ? 1.0 : -1.0)*fticks*ESCALATOR_SPEED/NUM_STAIRS_PER_FLOOR_ESC);
							pos.z += delta*height;
							pos[e.dim] += (e.dir ? 1.0 : -1.0)*delta*length;
						}
						obj_z   = max(pos.z, p_last.z);
						handled = had_coll = player_on_moving_ww = player_on_escalator = 1; // set player_on_moving_ww=1 to suppress footsteps
						// clamp player to interior of escalator
						float const shrink(min(radius, 0.5f*ramp_inner.get_sz_dim(!e.dim))); // make sure it's normalized
						max_eq(pos[!e.dim], ramp_inner.d[!e.dim][0]+shrink);
						min_eq(pos[!e.dim], ramp_inner.d[!e.dim][1]-shrink);
					}
				}
			}
			if (!handled) {
				// treat the sloped part of the escalator as an extruded polygon; should be conservative for balls and correct for player exterior collisions with bottom
				float val(0.0);
				vector3d coll_norm;

				if (get_sphere_poly_int_val(pos, radius, bot_pts, 4, normal, e.get_side_deltaz(), val, coll_norm)) {
					if (is_player && coll_norm.z < 0.0) {coll_norm.z = 0.0; coll_norm.normalize();} // don't push player through floor when colliding with bottom
					if (cnorm) {*cnorm = coll_norm;}
					pos += coll_norm*val;
					had_coll = 1;
				}
			}
		}
		e.get_ends_bcube(ends[0], ends[1], 0); // exclude_sides=0

		for (unsigned hi = 0; hi < 2; ++hi) {
			if (!sphere_cube_int_zval(pos, obj_z, radius, ends[hi])) continue;
			for (unsigned lr = 0; lr < 2; ++lr) {had_coll |= sphere_cube_int_update_pos(pos, radius, e.get_side_for_end(ends[hi], lr), p_last, 1, cnorm);} // skip_z=1
		}
		if (!e.in_mall) { // check support pillar; not present in malls
			cube_t const pillar(e.get_support_pillar());
			had_coll |= sphere_cube_int_update_pos_zval(pos, obj_z, radius, pillar, p_last, cnorm);
		}
	} // for escalator e
	for (auto i = doors.begin(); i != doors.end(); ++i) {
		had_coll |= check_door_coll(building, *i, pos, p_last, radius, obj_z, check_open_doors, cnorm);
	}
	if (has_mall()) {
		for (store_doorway_t const &d : mall_info->store_doorways) {
			if (d.open_amt == 1.0) continue; // fully open
			cube_t gate(d);
			gate.z1() = d.get_gate_z1();
			had_coll |= sphere_cube_int_update_pos_zval(pos, obj_z, radius, gate, p_last, cnorm);
		}
	}
	if (conn_info) { // check for collision with closed door separating the adjacent building at the end of the connecting room
		door_t const *const conn_door(conn_info->get_door_to_conn_part(building, pos));
		if (conn_door != nullptr) {had_coll |= check_door_coll(building, *conn_door, pos, p_last, radius, obj_z, check_open_doors, cnorm);}
	}
	return had_coll;
}

bool get_line_clip_update_t(point const &p1, point const &p2, cube_t const &c, float &t) {
	float tmin(0.0), tmax(1.0);
	if (get_line_clip(p1, p2, c.d, tmin, tmax) && tmin < t) {t = tmin; return 1;}
	return 0;
}
template<typename T> bool get_line_clip_update_t(point const &p1, point const &p2, vector<T> const &cubes, float &t) { // cube_t, cube_with_ix_t, etc.
	bool had_coll(0);
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {had_coll |= get_line_clip_update_t(p1, p2, *i, t);}
	return had_coll;
}
bool building_interior_t::line_coll(building_t const &building, point const &p1, point const &p2, point &p_int, bool skip_room_geom) const { // for snakes and ball throwing
	bool had_coll(0);
	float t(1.0), tmin(1.0);
	had_coll |= get_line_clip_update_t(p1, p2, floors,   t);
	had_coll |= get_line_clip_update_t(p1, p2, ceilings, t);
	for (unsigned d = 0; d < 2; ++d) {had_coll |= get_line_clip_update_t(p1, p2, walls[d], t);}

	for (elevator_t const &e : elevators) {
		if (!check_line_clip(p1, p2, e.d)) continue;
		cube_t cubes[5];
		unsigned const num_cubes(e.get_coll_cubes(cubes));
		for (unsigned n = 0; n < num_cubes; ++n) {had_coll |= get_line_clip_update_t(p1, p2, cubes[n], t);}
	}
	for (escalator_t const &e : escalators) {
		if (!check_line_clip(p1, p2, e.d)) continue;
		cube_t cubes[7];
		e.get_all_cubes(cubes);
		for (unsigned n = 0; n < 7; ++n) {had_coll |= get_line_clip_update_t(p1, p2, cubes[n], t);}
		had_coll |= get_line_clip_update_t(p1, p2, e.get_ramp_bcube(0), t); // use conservative bcube of ramp
	}
	for (auto i = doors.begin(); i != doors.end(); ++i) {
		if (i->open) {
			cube_t door_bounds(*i);
			door_bounds.expand_by_xy(i->get_width());
			if (!check_line_clip(p1, p2, door_bounds.d)) continue; // check intersection with rough/conservative door bounds (optimization)
			tquad_with_ix_t const door(building.set_interior_door_from_cube(*i));
			vector3d const normal(door.get_norm()), v1(p2 - p1);
			point pts[2][4];
			gen_poly_planes(door.pts, door.npts, normal, i->get_thickness(), pts);
			float const dp(dot_product(v1, normal));
			bool const test_side(dp > 0.0);
			
			if (thick_poly_intersect(v1, p1, normal, pts, test_side, door.npts)) {
				tmin = dot_product_ptv(normal, pts[!test_side][0], p1)/dp; // calculate t, since thick_poly_intersect() doesn't return it; use the far edge of the door
				if (tmin > 0.0 && tmin < t) {t = tmin; had_coll = 1;}
			}
		}
		else {had_coll |= get_line_clip_update_t(p1, p2, i->get_true_bcube(), t);} // closed
	} // for i
	if (room_geom && !skip_room_geom) { // check room geometry
		for (auto c = room_geom->objs.begin(); c != room_geom->objs.end(); ++c) { // check for other objects to collide with (including stairs)
			if (c->no_coll() || c->type == TYPE_BLOCKER || c->type == TYPE_ELEVATOR) continue; // skip blockers and elevators

			if (c->type == TYPE_CLOSET) { // special case to handle closet interiors
				cube_t cubes[5];
				get_closet_cubes(*c, cubes, 1); // get cubes for walls and door; for_collision=1
				unsigned const n_end(get_closet_num_coll_cubes(*c));
				for (unsigned n = 0; n < n_end; ++n) {had_coll |= get_line_clip_update_t(p1, p2, cubes[n], t);}
			}
			else if (c->is_vert_cylinder()) { // vertical cylinder
				if (line_intersect_cylinder_with_t(p1, p2, c->get_cylinder(), 1, tmin) && tmin < t) {t = tmin; had_coll = 1;}
			}
			else if (c->shape == SHAPE_SPHERE) { // sphere
				float const radius(c->get_radius());
				if (sphere_test_comp(p1, c->get_cube_center(), (p1 - p2), radius*radius, tmin) && tmin < t) {t = tmin; had_coll = 1;}
			}
			//else if (c->type == TYPE_STALL || c->type == TYPE_SHOWER) {}
			else {had_coll |= get_line_clip_update_t(p1, p2, *c, t);} // assume it's a cube
		} // for c
		for (auto const &i : room_geom->trim_objs) {had_coll |= get_line_clip_update_t(p1, p2, i, t);} // include wall trim
	}
	if (had_coll) {p_int = p1 + t*(p2 - p1);}
	return had_coll;
}

void update_closest_pt(cube_t const &cube, point const &pos, point &closest, float pad_dist, float &dmin_sq) {
	cube_t cube_pad(cube);
	cube_pad.expand_by(pad_dist);
	point const cand(cube_pad.closest_pt(pos));
	float const dsq(p2p_dist_sq(pos, cand));
	if (dmin_sq < 0.0 || dsq < dmin_sq) {closest = cand; dmin_sq = dsq;}
}
template<typename T> void update_closest_pt(vector<T> const &cubes, point const &pos, point &closest, float pad_dist, float &dmin_sq) { // cube_t, cube_with_ix_t, etc.
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {update_closest_pt(*c, pos, closest, pad_dist, dmin_sq);}
}
point building_interior_t::find_closest_pt_on_obj_to_pos(building_t const &building, point const &pos, float pad_dist, bool no_ceil_floor) const { // for tape
	float dmin_sq(-1.0); // start at an invalid value
	point closest(pos); // start at pt - will keep this value if there are no objects

	if (!no_ceil_floor) {
		update_closest_pt(floors,   pos, closest, pad_dist, dmin_sq);
		update_closest_pt(ceilings, pos, closest, pad_dist, dmin_sq);
	}
	for (unsigned d = 0; d < 2; ++d) {update_closest_pt(walls[d], pos, closest, pad_dist, dmin_sq);}
	for (elevator_t const &e : elevators) {update_closest_pt(e, pos, closest, pad_dist, dmin_sq);} // ignores open elevator doors
	
	for (escalator_t const &e : escalators) { // ramps are ignored
		cube_t cubes[7];
		e.get_all_cubes(cubes);
		for (unsigned n = 0; n < 7; ++n) {update_closest_pt(cubes[n], pos, closest, pad_dist, dmin_sq);}
	}
	for (auto i = doors.begin(); i != doors.end(); ++i) {
		if (i->open) {} // handle open doors? - closest point on extruded polygon
		else {update_closest_pt(*i, pos, closest, pad_dist, dmin_sq);} // closed
	}
	if (room_geom) { // check room geometry
		for (auto c = room_geom->objs.begin(); c != room_geom->objs.end(); ++c) { // check for other objects to collide with (including stairs)
			if (c->no_coll() || c->type == TYPE_BLOCKER || c->type == TYPE_ELEVATOR) continue; // skip blockers and elevators
			if (c->type == TYPE_RAILING || c->type == TYPE_PLANT || c->type == TYPE_PLANT_MODEL || c->type == TYPE_TREE) continue; // complex shapes that are too hard to attach to

			if (c->type == TYPE_CLOSET) { // special case to handle closet interiors
				cube_t cubes[5];
				get_closet_cubes(*c, cubes, 1); // get cubes for walls and door; for_collision=1
				unsigned const n_end(get_closet_num_coll_cubes(*c));
				for (unsigned n = 0; n < n_end; ++n) {update_closest_pt(cubes[n], pos, closest, pad_dist, dmin_sq);}
			}
			if (c->is_vert_cylinder()) { // vertical cylinder
				point const center(c->get_cube_center());
				float const radius(c->get_radius() + pad_dist), dsq_xy(max(0.0f, (p2p_dist_xy_sq(pos, center) - radius*radius)));
				float const zval(max(c->z1()-pad_dist, min(c->z2()+pad_dist, pos.z))), dz(pos.z - zval), dsq(dsq_xy + dz*dz);
				if (dmin_sq < 0.0 || dsq < dmin_sq) {closest = point(center.x, center.y, zval) + radius*(point(pos.x, pos.y, 0.0) - point(center.x, center.y, 0.0)).get_norm(); dmin_sq = dsq;}
			}
			else if (c->shape == SHAPE_SPHERE) { // sphere
				point const center(c->get_cube_center());
				float const radius(c->get_radius() + pad_dist), dsq(max(0.0f, (p2p_dist_sq(pos, center) - radius*radius)));
				if (dmin_sq < 0.0 || dsq < dmin_sq) {closest = center + radius*(pos - center).get_norm(); dmin_sq = dsq;}
			}
			else {update_closest_pt(*c, pos, closest, pad_dist, dmin_sq);} // assume it's a cube
		} // for c
		for (auto const &i : room_geom->trim_objs) {update_closest_pt(i, pos, closest, pad_dist, dmin_sq);} // include wall trim
	}
	// what about exterior building walls?
	return closest;
}

bool line_quad_intersect(point const &p1, point const &p2, point const *const pts, float &t) {
	float tmin(0.0);
	if (line_poly_intersect(p1, p2, pts, 4, get_poly_norm(pts), tmin) && tmin < t) {t = tmin; return 1;}
	return 0;
}
bool line_int_polygon_sides(point const &p1, point const &p2, cube_t const &bcube, vect_point const &points, float &t) {
	point quad_pts[4]; // quads
	bool hit(0);

	for (unsigned S = 0; S < points.size(); ++S) { // generate vertex data quads
		for (unsigned d = 0, ix = 0; d < 2; ++d) {
			point const &p(points[(S+d)%points.size()]);
			for (unsigned e = 0; e < 2; ++e) {quad_pts[ix++].assign(p.x, p.y, bcube.d[2][d^e]);}
		}
		hit |= line_quad_intersect(p1, p2, quad_pts, t);
	} // for S
	return hit;
}
template<typename T> bool check_line_int_cubes(vector<T> const &cubes, point const &p1, point const &p2, float &t) {
	bool ret(0);
	for (cube_t const &c : cubes) {ret |= get_line_clip_update_t(p1, p2, c, t);}
	return ret;
}

// Note: p1/p2 are in building space
unsigned building_t::check_line_coll(point const &p1, point const &p2, float &t, bool occlusion_only, bool ret_any_pt, bool no_coll_pt, bool check_non_coll) const {
	point p1r(p1), p2r(p2); // copy before clipping
	if (!check_line_clip(p1r, p2r, bcube.d)) return BLDG_COLL_NONE; // no intersection (never returns here for vertical lines)
	
	if (is_rotated()) {
		point const center(bcube.get_cube_center());
		do_xy_rotate_inv(center, p1r); // inverse rotate - negate the sine term
		do_xy_rotate_inv(center, p2r);
	}
	float const pzmin(min(p1r.z, p2r.z)), pzmax(max(p1r.z, p2r.z));
	bool const vert(p1r.x == p2r.x && p1r.y == p2r.y);
	float tmin(0.0);
	unsigned coll(BLDG_COLL_NONE);
	if (check_non_coll && check_line_int_cubes(skylights, p1r, p2r, t)) return BLDG_COLL_SKYLIGHT; // nothing else is over the skylight

	for (auto i = parts.begin(); i != parts.end(); ++i) {
		if (pzmin > i->z2() || pzmax < i->z1()) continue; // no overlap in z
		bool hit(0);

		if (use_cylinder_coll()) { // vertical cylinder
			// Note: we know the line intersects the cylinder's bcube, and there's a good chance it intersects the cylinder, so we don't need any expensive early termination cases here
			point const cc(i->get_cube_center());
			vector3d const csz(i->get_size());

			if (vert) { // vertical line + vertical cylinder optimization + handling of ellipsoids
				if (!point_in_ellipse(p1r, cc, 0.5*csz.x, 0.5*csz.y)) continue; // no intersection (below test should return true as well)
				tmin = (i->z2() - p1r.z)/(p2r.z - p1r.z);
				if (tmin >= 0.0 && tmin < t) {t = tmin; hit = 1;}
			}
			else {
				float const radius(0.5*(occlusion_only ? min(csz.x, csz.y) : max(csz.x, csz.y))); // use conservative radius unless this is an occlusion query
				point const cp1(cc.x, cc.y, i->z1()), cp2(cc.x, cc.y, i->z2());
				if (p1r.z > i->z1() && p1r.z < i->z2() && dist_xy_less_than(cp1, p1r, radius)) {tmin = 0.0;} // p1 inside the cylinder
				else if (!line_int_cylinder(p1r, p2r, cp1, cp2, radius, radius, 1, tmin) || tmin >= t) continue; // conservative for non-occlusion rays

				if (!occlusion_only && csz.x != csz.y) { // ellipse; approximate
					vector3d const delta(p2r - p1r);
					float const rx_inv_sq(1.0/(0.25*csz.x*csz.x)), ry_inv_sq(1.0/(0.25*csz.y*csz.y));
					float t_step(0.1*max(csz.x, csz.y)/delta.mag());

					for (unsigned n = 0; n < 10; ++n) { // use an interative approach
						point const p(p1r + tmin*delta);
						if (point_in_ellipse_risq(p, cc, rx_inv_sq, ry_inv_sq)) {hit = 1; tmin -= t_step;} else {tmin += t_step;}
						if (hit) {t_step *= 0.5;} // converge on hit point
					}
					if (!hit) continue; // not actually a hit
				} // end ellipse case
				if (tmin < t) {t = tmin; hit = 1;} // re-check tmin in case it increased above t in the iterations above
			}
		}
		else if (!is_cube()) {
			vect_point const &points(get_part_ext_verts(i - parts.begin()));

			if (p1r.z != p2r.z) { // test vertical lines and top surface
				float const tz((i->z2() - p1r.z)/(p2r.z - p1r.z)); // t value at zval = top of cube

				if (tz >= 0.0 && tz < t) {
					float const xval(p1r.x + tz*(p2r.x - p1r.x)), yval(p1r.y + tz*(p2r.y - p1r.y));
					if (point_in_polygon_2d(xval, yval, points.data(), points.size())) {t = tz; hit = 1;} // XY plane test
				}
			}
			if (!vert && line_int_polygon_sides(p1r, p2r, *i, points, t)) {hit = 1;} // test building sides
		}
		else {hit |= get_line_clip_update_t(p1r, p2r, *i, t);} // cube

		if (hit) {
			if (occlusion_only) return BLDG_COLL_SIDE; // early exit
			if (vert) {coll = BLDG_COLL_ROOF;} // roof
			else {
				float const zval(p1r.z + t*(p2.z - p1r.z));
				coll = ((fabs(zval - i->z2()) < 0.0001*i->dz()) ? BLDG_COLL_ROOF : BLDG_COLL_SIDE); // test if clipped zval is close to the roof zval
			}
			if (ret_any_pt) return coll;
		}
	} // for i
	if (occlusion_only) return BLDG_COLL_NONE;
	if (check_line_int_cubes(details, p1r, p2r, t)) {coll = BLDG_COLL_DETAIL;} // details cube

	if (!no_coll_pt || !vert) { // vert line already tested building cylins/cubes, and marked coll roof, no need to test again unless we need correct coll_pt t-val
		for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
			if (line_poly_intersect(p1r, p2r, i->pts, i->npts, i->get_norm(), tmin) && tmin < t) {t = tmin; coll = BLDG_COLL_ROOF;} // roof quad
		}
	}
	if (!vert || check_non_coll) { // don't need to check fences for vertical collisions since they're horizontal blockers
		if (check_line_int_cubes(fences, p1r, p2r, t)) {coll = BLDG_COLL_FENCE;} // details cube
	}
	return coll; // Note: no collisions with windows or doors, since they're colinear with walls; no collision with interior for now
}

unsigned const OBJS_GRID_SZ = 4;

class building_color_query_geom_cache_t {
	struct cube_with_color_t : public cube_t {
		colorRGBA color;
		bool is_round;
		cube_with_color_t(cube_t const &cube, colorRGBA const &c, bool is_round_) : cube_t(cube), color(c), is_round(is_round_) {}
	};
	vector<cube_with_color_t> objs[OBJS_GRID_SZ*OBJS_GRID_SZ];
	cube_t bcube;
	float gxy_inv[2] = {};
	int cur_frame = 0;

	void add_obj(cube_with_color_t const &obj) {
		unsigned gix[2][2] = {}; // {x,y}, {lo,hi}

		for (unsigned d = 0; d < 2; ++d) {
			for (unsigned e = 0; e < 2; ++e) {gix[d][e] = min(OBJS_GRID_SZ-1U, unsigned(max(0.0f, (obj.d[d][e] - bcube.d[d][0])*gxy_inv[d])));}
		}
		for (unsigned y = gix[1][0]; y <= gix[1][1]; ++y) {
			for (unsigned x = gix[0][0]; x <= gix[0][1]; ++x) {objs[y*OBJS_GRID_SZ + x].push_back(obj);}
		}
	}
public:
	void gen_objs(building_t const &building, float z1, float z2) {
		if (cur_frame == frame_counter) return; // already handled by another thread
		assert(building.has_room_geom());
		for (unsigned n = 0; n < OBJS_GRID_SZ*OBJS_GRID_SZ; ++n) {objs[n].clear();}
		bcube = building.bcube;
		for (unsigned d = 0; d < 2; ++d) {gxy_inv[d] = OBJS_GRID_SZ/bcube.get_sz_dim(d);}
		float const sz_thresh(0.5*building.get_wall_thickness());

		for (unsigned d = 0; d < 2; ++d) { // add walls
			for (cube_t const &wall : building.interior->walls[d]) {
				if (z1 < wall.z2() && z2 > wall.z1()) {add_obj(cube_with_color_t(wall, WHITE, 0));} // is_round=0
			}
		}
		for (person_t const &p : building.interior->people) { // add people (drawn in red)
			cube_t const bcube(p.get_bcube());
			if (z1 < bcube.z2() && z2 > bcube.z1()) {add_obj(cube_with_color_t(bcube, RED, 1));} // is_round=1
		}
		if (building.has_room_geom()) {
			for (room_object_t const &obj : building.interior->room_geom->objs) { // add room objects
				if (obj.flags & RO_FLAG_INVIS) continue;
				room_object const type(obj.type);
				if (type != TYPE_PIPE && min(obj.dx(), obj.dy()) < sz_thresh) continue; // too small; exclude pipes
				if (type == TYPE_BOOK || type == TYPE_PLANT || type == TYPE_RAILING || obj.is_a_drink() || type == TYPE_PAPER || type == TYPE_PAINTCAN ||
					type == TYPE_WBOARD || type == TYPE_DRAIN || type == TYPE_PLATE || type == TYPE_LBASKET || type == TYPE_LAMP || type == TYPE_CUP || type == TYPE_LAPTOP ||
					is_ball_type(type) || type == TYPE_PAN || type == TYPE_PG_BEAM || type == TYPE_PLANT_MODEL || type == TYPE_PICTURE || type == TYPE_WINDOW ||
					type == TYPE_SIGN || type == TYPE_WALL_TRIM || type == TYPE_BLINDS || type == TYPE_PEN || type == TYPE_PENCIL || type == TYPE_HANGER_ROD ||
					type == TYPE_MONEY || type == TYPE_PHONE || type == TYPE_TPROLL || type == TYPE_SPRAYCAN || type == TYPE_MARKER || type == TYPE_BUTTON || type == TYPE_VALVE ||
					type == TYPE_CRACK || type == TYPE_SWITCH || type == TYPE_TAPE || type == TYPE_OUTLET || type == TYPE_CURB || type == TYPE_BRK_PANEL || type == TYPE_VENT ||
					type == TYPE_BREAKER || type == TYPE_DUCT || type == TYPE_TOY || type == TYPE_VASE || type == TYPE_URN || type == TYPE_STAPLER || type == TYPE_WIND_SILL ||
					type == TYPE_SPRINKLER || type == TYPE_FEXT_MOUNT || type == TYPE_FEXT_SIGN || type == TYPE_PIZZA_TOP || type == TYPE_TEESHIRT || type == TYPE_PANTS ||
					type == TYPE_DBG_SHAPE || type == TYPE_POOL_CUE || type == TYPE_WALL_MOUNT || type == TYPE_FALSE_DOOR || type == TYPE_FLASHLIGHT || type == TYPE_CANDLE ||
					type == TYPE_CAMERA || type == TYPE_CLOCK || type == TYPE_DOWNSPOUT || type == TYPE_CHIM_CAP || type == TYPE_FOOD_BOX || type == TYPE_LADDER ||
					type == TYPE_LAVALAMP || type == TYPE_POOL_LAD || type == TYPE_PADLOCK || type == TYPE_KEY || type == TYPE_HANGER || type == TYPE_CLOTHES ||
					type == TYPE_WALL_LAMP || type == TYPE_SILVER || type == TYPE_TOY_MODEL || type == TYPE_CEIL_FAN || type == TYPE_FOLD_SHIRT || type == TYPE_TRASH ||
					type == TYPE_INT_WINDOW || type == TYPE_INT_LADDER || type == TYPE_CONF_PHONE || type == TYPE_SPIWEB || type == TYPE_TREE || type == TYPE_ELEC_WIRE ||
					type == TYPE_ERASER || type == TYPE_SHOE) continue;
				if (z1 > obj.z2() || z2 < obj.z1()) continue; // zval test

				if (obj.type == TYPE_PARK_SPACE) {
					if (!obj.is_used()) continue; // no car in this space
					pair<cube_t, colorRGBA> const ret(car_bcube_color_from_parking_space(obj)); // Note: currently always white
					add_obj(cube_with_color_t(ret.first, ret.second, 0)); // is_round=0
				}
				else {add_obj(cube_with_color_t(obj, obj.get_color(), (obj.is_round() && obj.type != TYPE_PIPE)));}
			} // for obj
		}
		cur_frame = frame_counter;
	}
	bool query_objs(building_t const &building, point const &pos, float z1, float z2, colorRGBA &color) {
		if (frame_counter != cur_frame) {
			// the first thread to get here generates the objects; is this always thread safe?
#pragma omp critical(collect_building_objects)
			gen_objs(building, z1, z2);
		}
		float zmax(z1);
		unsigned gix[2] = {}; // {x,y}
		for (unsigned d = 0; d < 2; ++d) {gix[d] = min(OBJS_GRID_SZ-1U, unsigned(max(0.0f, (pos[d] - bcube.d[d][0])*gxy_inv[d])));}

		for (cube_with_color_t const &obj : objs[OBJS_GRID_SZ*gix[1] + gix[0]]) { // return topmost object; zval is already checked
			if (!obj.contains_pt_xy(pos) || obj.z2() <= zmax) continue;
			if (obj.is_round && !dist_xy_less_than(pos, obj.get_cube_center(), 0.5*obj.dx())) continue; // handle cylinders and spheres
			color = obj.color;
			zmax  = obj.z2();
		}
		return (zmax > z1);
	}
};

building_color_query_geom_cache_t building_color_query_geom_cache;

bool building_t::get_interior_color_at_xy(point const &pos_in, colorRGBA &color) const {
	if (!interior || !is_cube() || is_rotated()) return 0; // these cases aren't handled
	bool const player_in_this_building(camera_in_building && player_building == this);
	if (!player_in_this_building) return 0; // not currently handled; maybe could allow this if the zoom level is high enough
	point pos(pos_in); // set zval to the player's height if in this building, otherwise use the ground floor
	pos.z = (player_in_this_building ? camera_pos.z : (ground_floor_z1 + get_floor_thickness()));
	bool cont_in_part(0);

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
		if (i->contains_pt(pos)) {cont_in_part = 1; break;} // only check zval if player in building
	}
	// check the extended basement; note that this doesn't check for the player inside a room, but the player shouldn't really be underground within the bcube otherwise;
	// also note that we may not even get into this function if the ray doesn't intersect this building's above ground part
	if (!cont_in_part && has_ext_basement() && interior->basement_ext_bcube.contains_pt(pos)) {
		if (has_parking_garage) {cont_in_part = 1;} // parking garage + backrooms; fills the entire ext basement area
		else {cont_in_part = (get_room_containing_pt(pos) >= 0);} // check for actual room containment
	}
	if (!cont_in_part) return 0;
	float const z1(pos.z - 1.05*CAMERA_RADIUS), z2(z1 + get_floor_ceil_gap()); // approx span of one floor, including objects on the floor
	if (building_color_query_geom_cache.query_objs(*this, pos, z1, z2, color)) return 1;
	color = (is_house ? LT_BROWN : GRAY); // floor
	return 1;
}

bool building_t::point_in_mall_elevator_entrance(point const &pos, bool inc_front_space) const {
	if (!has_mall() || interior->mall_info->city_elevator_ix < 0) return 0; // no mall or above ground elevator
	if (pos.z < ground_floor_z1 - get_fc_thickness())             return 0; // below ground
	elevator_t const &e(get_elevator(interior->mall_info->city_elevator_ix));
	return ((inc_front_space ? e.get_bcube_padded(1.5*e.get_length()) : (cube_t)e).contains_pt(pos)); // in or standing in front of the mall underground elevator
}

void expand_convex_polygon_xy(vect_point &points, point const &center, float expand) {
	if (expand == 0.0) return;

	for (auto p = points.begin(); p != points.end(); ++p) {
		vector3d const dir((p->x - center.x), (p->y - center.y), 0.0); // only want XY component
		*p += dir*(expand/dir.mag());
	}
}

// Note: xy_radius == 0.0 is a point test; xy_radius > 0.0 is an intersection test;
// attic and basement queries only work with points; the xy_radius != 0 and coll_cube cases are only used with pedestrians
// return value: 0=no cont, 1=part, 2=attic, 3=ext basement, 4=roof access, 5=detail
int building_t::check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points,
	bool inc_attic, bool inc_ext_basement, bool inc_roof_acc, bool inc_details, cube_t *coll_cube) const
{
	if (coll_cube) {assert(!inc_attic && !inc_ext_basement && !inc_roof_acc);} // not supported
	point const pr(get_inv_rot_pos(pos));

	if (inc_ext_basement && point_in_extended_basement_not_basement(pr)) { // extended basement is not rotated
		// this check must be accurate; since the extended basement is sparse, we need to check every extended basement room
		// expand by wall thickness to avoid failing when crossing between two rooms that don't exactly line up, such as with conn room boundary with adj building
		return (interior->point_in_ext_basement_room(pr, get_window_vspace(), get_wall_thickness()) ? 3 : 0);
	}
	if (inc_ext_basement && has_pool() && interior->pool.contains_pt(pr)) return 3; // in the pool
	if (inc_ext_basement && point_in_mall_elevator_entrance(pr, 1))       return 3; // in or standing in front of the mall underground elevator; inc_front_space=1
	
	if (bcube.contains_pt_exp_xy_only(pos, xy_radius)) { // check parts
		if (inc_attic && point_in_attic(pr)) return 2;

		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (pr.z > i->z2() || pr.z < i->z1()) continue; // no overlap in z

			if (use_cylinder_coll()) { // vertical cylinder
				point const cc(i->get_cube_center());
				vector3d const csz(i->get_size());
				float const dx(cc.x - pr.x), dy(cc.y - pr.y), rx(0.5*csz.x + xy_radius), ry(0.5*csz.y + xy_radius);
				if (dx*dx/(rx*rx) + dy*dy/(ry*ry) > 1.0f) continue; // no intersection (below test should return true as well)
			}
			else if (num_sides != 4) {
				points = get_part_ext_verts(i - parts.begin()); // deep copy so that we can modify points below to apply the expand
				expand_convex_polygon_xy(points, i->get_cube_center(), xy_radius); // cylinder case: expand polygon by xy_radius; assumes a convex polygon
				if (!point_in_polygon_2d(pr.x, pr.y, points.data(), points.size())) continue; // XY plane test for top surface
			}
			else { // cube
				if (!i->contains_pt_xy_exp(pr, xy_radius)) continue;
			}
			if (coll_cube) {*coll_cube = *i;} // approximate for cylinder
			return 1;
		} // for i
	}
	if (has_porch()) { // check porch pillar
		assert(real_num_parts+1U < parts.size()); // must have {porch roof, porch support pillar}
		cube_t const &porch_pillar(parts[real_num_parts+1]);
		
		if (porch_pillar.contains_pt_xy_exp(pr, xy_radius)) {
			if (coll_cube) {*coll_cube = porch_pillar;}
			return 1;
		}
	}
	if (inc_roof_acc && !is_house && interior) { // check roof access stairs
		for (stairwell_t const &s : interior->stairwells) {
			if (s.roof_access && s.contains_pt_xy(pr) && pr.z > s.z1() && pr.z < s.z2() + get_window_vspace()) return 4;
		}
	}
	if (inc_details) {
		for (auto const &i : details) { // AC unit, balcony post, etc.
			if (i.contains_pt_exp_xy_only(pr, xy_radius)) {
				if (coll_cube) {*coll_cube = i;}
				return 5;
			}
		}
		for (cube_t const &fence : fences) {
			if (fence.contains_pt_exp_xy_only(pr, xy_radius)) {
				if (coll_cube) {*coll_cube = fence;}
				return 5;
			}
		}
	}
	return 0;
}

bool building_t::check_point_xy_in_part(point const &pos) const { // simpler/faster version of check_point_or_cylin_contained() with no z check
	if (!bcube.contains_pt_xy(pos)) return 0; // no intersection (bcube does not need to be rotated)
	point const pr(get_inv_rot_pos(pos));

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
		if (i->contains_pt_xy(pr)) return 1;
	}
	return 0;
}

template <typename T> bool has_cube_line_coll(point const &p1, point const &p2, vector<T> const &cubes) {
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {
		if (i->line_intersects(p1, p2)) return 1;
	}
	return 0;
}
bool building_t::check_for_wall_ceil_floor_int(point const &p1, point const &p2, bool inc_pg_br_walls) const { // and interior doors/windows/gates
	if (!interior) return 0;

	for (unsigned d = 0; d < 2; ++d) {
		if (has_cube_line_coll(p1, p2, interior->walls[d])) return 1;
	}
	if (p1.z != p2.z && has_cube_line_coll(p1, p2, interior->fc_occluders)) return 1; // skip for horizontal lines
	if (check_line_intersect_doors(p1, p2)) return 1;

	if (inc_pg_br_walls && has_parking_garage && has_room_geom() && min(p1.z, p2.z) < ground_floor_z1) {
		// Note: not using get_pgbr_wall_ix_for_pos() because p1 and p2 can be on different floors, or one in PG and the other in BR; also, only called a few times per frame
		for (unsigned d = 0; d < 2; ++d) {
			if (has_cube_line_coll(p1, p2, interior->room_geom->pgbr_walls[d])) return 1;
		}
	}
	if (check_line_int_interior_window(p1, p2)) return 1;

	if (has_mall()) { // check closed mall store gates
		for (store_doorway_t const &d : interior->mall_info->store_doorways) {
			if (d.open_amt == 1.0) continue; // fully open
			cube_t gate(d);
			gate.z1() = d.get_gate_z1();
			if (gate.line_intersects(p1, p2)) return 1;
		}
	}
	return 0;
}
bool building_t::line_intersect_stairs_or_ramp(point const &p1, point const &p2) const {
	if (!interior) return 0;
	if (has_pg_ramp() && interior->pg_ramp.line_intersects(p1, p2)) return 1;
	return has_cube_line_coll(p1, p2, interior->stairwells);
}

bool building_t::check_cube_on_or_near_stairs(cube_t const &c) const {
	float const doorway_width(get_doorway_width());

	for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
		cube_t tc(*s);
		tc.expand_in_dim(s->dim, doorway_width); // add extra space at both ends of stairs
		if (tc.intersects(c)) return 1;
	}
	return 0;
}

// Note: for procedural object placement; no expanded_objs, but includes blockers
bool building_t::overlaps_other_room_obj(cube_t const &c, unsigned objs_start, bool check_all, unsigned const *objs_end) const {
	assert(has_room_geom());
	vect_room_object_t &objs(interior->room_geom->objs);
	auto start(objs.begin()+objs_start), end((objs_end == nullptr) ? objs.end() : (objs.begin()+*objs_end));
	assert(objs_start <= objs.size() && start <= end);

	for (auto i = start; i != end; ++i) {
		if (i->type == TYPE_POOL_TILE) continue; // always excluded, since it's thin and objects can be mounted over it
		if (i->is_a_drink() && i->is_on_floor() && i->intersects(c)) return 1; // bottles and cans on the floor do count
		// Note: light switches/outlets/vents/pipes don't collide with the player or AI, but they collide with other placed objects to avoid blocking them;
		// however, it's okay to block outlets with furniture
		if ((check_all || !i->no_coll() || i->type == TYPE_SWITCH || i->type == TYPE_OUTLET || i->type == TYPE_VENT || i->type == TYPE_PIPE ||
			i->type == TYPE_FALSE_DOOR || i->is_pet_container()) && i->intersects(c)) return 1;
		if (i->type == TYPE_DESK && i->shape == SHAPE_TALL && i->intersects_xy_no_adj(c) && c.intersects_no_adj(get_desk_top_back(*i))) return 1; // check tall desk back
		if (i->type == TYPE_BOOK && (i->flags & RO_FLAG_ON_FLOOR) && i->intersects(c))   return 1; // books on floors count
	} // for i
	return 0;
}
// Note: for dynamic player object placement and book opening; includes expanded_objs but not blockers
bool building_t::overlaps_any_placed_obj(cube_t const &c) const { // Note: includes expanded_objs
	assert(has_room_geom());

	for (auto i = interior->room_geom->objs.begin(); i != interior->room_geom->objs.end(); ++i) { // include stairs and elevators
		if (i->type != TYPE_BLOCKER && i->intersects(c)) return 1; // only exclude blockers; maybe should use rat_coll flag?
	}
	for (auto i = interior->room_geom->expanded_objs.begin(); i != interior->room_geom->expanded_objs.end(); ++i) {
		if (i->type != TYPE_BLOCKER && i->intersects(c)) return 1; // maybe should use rat_coll flag?
	}
	return 0;
}

bool building_t::overlaps_or_adj_int_window(cube_t const &c) const {
	if (!interior || interior->int_windows.empty()) return 0;
	cube_t c_exp(c);
	c_exp.expand_by_xy(get_wall_thickness()); // account for frame thickness
	return has_bcube_int(c_exp, interior->int_windows);
}

bool building_t::check_skylight_intersection(cube_t const &c) const {
	if (skylights.empty()) return 0;
	cube_t test_cube(c);
	test_cube.expand_in_z(0.5*get_window_vspace()); // expand to include half a floor of space above and below the skylight
	return has_bcube_int(test_cube, skylights);
}

bool line_int_cube_exp(point const &p1, point const &p2, cube_t const &cube, vector3d const &expand) {
	cube_t tc(cube);
	tc.expand_by(expand);
	return tc.line_intersects(p1, p2);
}
bool line_int_cubes_exp(point const &p1, point const &p2, cube_t const *const cubes, unsigned num_cubes, vector3d const &expand) {
	for (unsigned n = 0; n < num_cubes; ++n) {
		if (line_int_cube_exp(p1, p2, cubes[n], expand)) return 1;
	}
	return 0;
}
// specialized for cube_t, line_bcube calculated internally
bool line_int_cubes_exp(point const &p1, point const &p2, vect_cube_t const &cubes, vector3d const &expand) {
	cube_t line_bcube(p1, p2);
	line_bcube.expand_by(expand);

	for (auto const &c : cubes) {
		if (line_bcube.intersects(c) && line_int_cube_exp(p1, p2, c, expand)) return 1;
	}
	return 0;
}
template<typename T> bool line_int_cubes_exp(point const &p1, point const &p2, vector<T> const &cubes, vector3d const &expand, cube_t const &line_bcube) {
	for (auto const &c : cubes) {
		if (line_bcube.intersects(c) && line_int_cube_exp(p1, p2, c, expand)) return 1;
	}
	return 0;
}
template bool line_int_cubes_exp(point const &p1, point const &p2, vector<cube_t> const &cubes, vector3d const &expand, cube_t const &line_bcube); // explicit instantiation

template<typename T> bool line_int_cubes(point const &p1, point const &p2, vector<T> const &cubes, cube_t const &line_bcube) {
	for (auto const &c : cubes) {
		if (line_bcube.intersects(c) && c.line_intersects(p1, p2)) return 1;
	}
	return 0;
}
template bool line_int_cubes(point const &p1, point const &p2, vector<cube_t> const &cubes, cube_t const &line_bcube); // explicit instantiation

unsigned get_ksink_cubes(room_object_t const &sink, cube_t cubes[3]) {
	assert(sink.type == TYPE_KSINK);
	cube_t dishwasher;
	if (!get_dishwasher_for_ksink(sink, dishwasher)) {cubes[0] = sink; return 1;} // no dishwahser - 1 cube
	room_object_t cabinet(sink);
	cubes[0] = split_cabinet_at_dishwasher(cabinet, dishwasher); // left part
	cubes[1] = cabinet; // right part
	cubes[2] = dishwasher;
	cubes[2].d[sink.dim][!sink.dir] = sink.d[sink.dim][!sink.dir]; // use the back of the cabinet, not the back of the dishwasher door
	return 3;
}

class cached_room_objs_t {
	vect_room_object_t rat_snake_objs, spider_objs; // Note: snakes use the same objs as rats because they're on the floor
	building_t const *building;
	int cur_frame;
public:
	cached_room_objs_t() : building(nullptr), cur_frame(0) {}

	void ensure_cached(building_t const &b) { // for now, we cache objects for both rats and spiders
		if (building == &b && cur_frame == frame_counter) return; // already cached
		building  = &b;
		cur_frame = frame_counter;
		rat_snake_objs.clear();
		spider_objs   .clear();
		if (global_building_params.num_rats_max    > 0 || global_building_params.num_snakes_max > 0) {b.get_objs_at_or_below_ground_floor(rat_snake_objs, 0);}
		if (global_building_params.num_spiders_max > 0) {b.get_objs_at_or_below_ground_floor(spider_objs, 1);}
	}
	vect_room_object_t const &get_objs(bool for_spider) const {return (for_spider ? spider_objs : rat_snake_objs);}
};
cached_room_objs_t cached_room_objs;

bool room_object_t::is_floor_collidable() const {return (bldg_obj_types[type].rat_coll || is_floor_clutter());}

bool room_object_t::is_spider_collidable() const { // include objects on the floor, walls, and ceilings
	if (!is_floor_collidable()) { // below are the exceptions: objects that are not floor collidable but spiders can walk on
		if (type != TYPE_LIGHT && type != TYPE_BRSINK && type != TYPE_MIRROR && type != TYPE_MWAVE && type != TYPE_HANGER_ROD &&
			type != TYPE_LAPTOP && type != TYPE_MONITOR && type != TYPE_CLOTHES && type != TYPE_TOASTER && type != TYPE_CABINET &&
			type != TYPE_INT_WINDOW && type != TYPE_IBEAM && !is_pet_container()) return 0;
	}
	// these objects are ignored by spiders and will be walked through; ceiling fan is probably okay because it connects to the ceiling with a small cylinder
	if (type == TYPE_CEIL_FAN || type == TYPE_BANANA || type == TYPE_BAN_PEEL) return 0; // not a cube
	if (type == TYPE_BOOK) return 0; // I guess books don't count, since they're too small to walk on?
	return 1;
}
bool room_object_t::is_vert_cylinder() const {
	if (shape != SHAPE_CYLIN) return 0; // not a cylinder
	if (type == TYPE_LIGHT  && (flags & RO_FLAG_ADJ_HI)) return 0; // horizontal cylinder wall light
	if (type == TYPE_PLATE  && get_length() < dz())      return 0; // vertical plate in dishwasher
	if (type == TYPE_TPROLL && !was_expanded())          return 0; // horizontal TP roll on holder
	if (type == TYPE_BUTTON || type == TYPE_HANGER_ROD || type == TYPE_CLOCK) return 0; // horizontal cylinders
	if (type != TYPE_DUCT   && type != TYPE_PIPE)        return 1; // only ducts and pipes can be horizontal cylinders
	return dir; // duct/pipe encoding for vertical is dim=x, dir=1
}

bool building_t::is_ground_floor_excluding_retail(float zval) const {
	if (zval == ground_floor_z1) return 1;
	if (!has_retail())           return 0;
	return (fabs(zval - (ground_floor_z1 + retail_floor_levels*get_window_vspace())) < get_floor_thickness()); // close to the second floor zval
}
float building_t::get_ground_floor_z_thresh(bool for_spider) const {
	return (ground_floor_z1 + (for_spider ? 1.0f : 0.25f)*get_window_vspace()); // rats are on the ground, while spiders can climb walls and insects can fly
}
bool building_t::get_begin_end_room_objs_on_ground_floor(float zval, bool for_spider, vect_room_object_t::const_iterator &b, vect_room_object_t::const_iterator &e) const {
	if (zval < get_ground_floor_z_thresh(for_spider)) { // optimized for the case of rats and spiders where most are on the ground floor or basement
		cached_room_objs.ensure_cached(*this);
		b = cached_room_objs.get_objs(for_spider).begin();
		e = cached_room_objs.get_objs(for_spider).end();
		return 1; // use cached objects
	}
	b = interior->room_geom->objs.begin();
	e = (for_spider ? interior->room_geom->objs.end() : interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators for rats, but include stairs for spiders
	return 0; // use standard objects
}
void building_t::get_objs_at_or_below_ground_floor(vect_room_object_t &ret, bool for_spider) const {
	float const z_thresh(get_ground_floor_z_thresh(for_spider));
	assert(has_room_geom());
	// skip buttons/stairs/elevators for rats, but include stairs for spiders
	auto objs_end(for_spider ? interior->room_geom->objs.end() : interior->room_geom->get_placed_objs_end());

	for (auto c = interior->room_geom->objs.begin(); c != objs_end; ++c) {
		if (c->z1() < z_thresh && c->is_collidable(for_spider)) {ret.push_back(*c);}
	}
	for (auto c = interior->room_geom->expanded_objs.begin(); c != interior->room_geom->expanded_objs.end(); ++c) {
		if (c->z1() < z_thresh && c->is_collidable(for_spider)) {ret.push_back(*c);}
	}
}

void get_approx_car_cubes(room_object_t const &cb, cube_t cubes[5]) {
	vector3d const sz(cb.get_size());
	float const bot_zval(cb.z1() + 0.12*sz.z);
	cubes[0] = cb; // body
	cubes[0].z1() = bot_zval;

	for (unsigned bf = 0; bf < 2; ++bf) { // {back, front}
		for (unsigned lr = 0; lr < 2; ++lr) { // {left, right, front}
			float const signed_side_len((bf ? -1.0 : 1.0)*sz[ cb.dim]);
			cube_t &c(cubes[2*bf + lr + 1]);
			c = cb;
			c.z2() = bot_zval;
			c.d[ cb.dim][ bf] += 0.1*signed_side_len; // shift in from ends
			c.d[ cb.dim][!bf] -= 0.7*signed_side_len; // middle of side is open
			c.d[!cb.dim][!lr] -= (lr ? -1.0 : 1.0)*0.85*sz[!cb.dim]; // middle 70% of ends is open
		} // for lr
	} // for bf
}

// for spiders; lg_cubes and sm_cubes are currently handled the same
void building_t::get_room_obj_cubes(room_object_t const &c, point const &pos, vect_cube_t &lg_cubes, vect_cube_t &sm_cubes, vect_cube_t &non_cubes) const {
	room_object const type(c.type);
	if (c.is_round() || c.is_floor_clutter()) {non_cubes.push_back(c);}
	else if (type == TYPE_RAILING || type == TYPE_RAMP || type == TYPE_BALCONY || type == TYPE_POOL_LAD || type == TYPE_OFF_CHAIR || type == TYPE_BAR_STOOL ||
		type == TYPE_LAVALAMP || type == TYPE_WFOUNTAIN || type == TYPE_INT_LADDER || type == TYPE_MACHINE || type == TYPE_CONF_PHONE || type == TYPE_XFORMER ||
		type == TYPE_US_FLAG || type == TYPE_BLDG_FOUNT || type == TYPE_THEFT_SENS)
	{
		non_cubes.push_back(c); // non-cubes; bar stools are close, should they be included?
		// allow walking on the floor above a parking garage ramp if there's no cutout; shrink ramp bcube to the ceiling of the top floor of the parking garage
		if (type == TYPE_RAMP && interior->ignore_ramp_placement && c.z2() >= ground_floor_z1) {non_cubes.back().z2() -= get_floor_thickness();}
	}
	else if (type == TYPE_SHELVES) {
		if (c.in_mall()) { // mall shelves have fishtanks (pet stores) and flat items (clothing stores), so we can handle walking on the shelves
			cube_t shelves[4]; // max number of shelves
			unsigned const num_shelves(get_shelves_for_object(c, shelves));
			lg_cubes.insert(lg_cubes.end(), shelves, shelves+num_shelves);
			// add brackets for small pet store spiders to walk on
			static vect_cube_with_ix_t brackets;
			get_shelf_brackets(c, shelves, num_shelves, brackets);
			vector_add_to(brackets, sm_cubes);
		}
		else { // storage room shelves have all sorts of items when expanded, which we don't collide with, so mark as non-cube
			non_cubes.push_back(get_shelves_no_bot_gap(c)); // allow spiders to crawl under shelves
		}
	}
	else if (type == TYPE_CLOSET && (c.is_open() || c.contains_pt(pos))) {
		cube_t cubes[5];
		get_closet_cubes(c, cubes, 1); // get cubes for walls and door; for_collision=1
		lg_cubes.insert(lg_cubes.end(), cubes, (cubes + get_closet_num_coll_cubes(c)));
	}
	else if (type == TYPE_BED) {
		cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
		get_bed_cubes(c, cubes);
		cubes[3].z1() = cubes[0].z1(); // extend mattress downward to include the frame
		lg_cubes.insert(lg_cubes.end(), cubes+1, cubes+5); // head, foot, mattress (+ frame), pillow (or should pillow be small?)
		get_tc_leg_cubes(cubes[5], c, BED_HEAD_WIDTH, 0, cubes); // cubes[5] is not overwritten
		sm_cubes.insert(sm_cubes.end(), cubes, cubes+4); // legs are small
	}
	else if (can_use_table_coll(c)) { // objects with legs
		cube_t cubes[7];
		unsigned const num(get_table_like_object_cubes(c, cubes));
		unsigned const legs_end(min(num, 5U)); // num_legs+1; typically 5, but can be less
		assert(num >= 2); // top and at least one leg-like cube
		lg_cubes.push_back(cubes[0]); // top of table, etc.
		sm_cubes.insert(sm_cubes.end(), cubes+1, cubes+legs_end); // skip table top; legs are small
		sm_cubes.insert(lg_cubes.end(), cubes+legs_end, cubes+num); // middle, drawers, and back; may be empty range
	}
	else if (type == TYPE_CONF_TABLE) {
		cube_t cubes[2]; // {top, base}
		get_conf_table_cubes(c, cubes);
		lg_cubes.insert(lg_cubes.end(), cubes, cubes+2);
	}
	else if (type == TYPE_RDESK) {
		cube_t cubes[3];
		get_reception_desk_cubes(c, cubes);
		lg_cubes.insert(lg_cubes.end(), cubes, cubes+3);
	}
	else if (type == TYPE_CHAIR) {
		cube_t cubes[3], leg_cubes[4]; // seat, back, legs_bcube
		get_chair_cubes(c, cubes);
		lg_cubes.insert(lg_cubes.end(), cubes, cubes+2); // seat, back
		get_tc_leg_cubes(cubes[2], c, c.get_chair_leg_width(), 1, leg_cubes);
		sm_cubes.insert(sm_cubes.end(), leg_cubes, leg_cubes+4); // legs are small
	}
	else if (type == TYPE_POOL_TABLE) {
		cube_t cubes[5]; // body + 4 legs
		get_pool_table_cubes(c, cubes);
		lg_cubes.insert(lg_cubes.end(), cubes, cubes+5); // are the legs lg_cubes or sm_cubes? they're larger than table legs
	}
	else if (type == TYPE_TUB) {
		cube_t cubes[5]; // bottom, front, back, 2 sides
		get_tub_cubes(c, cubes);
		lg_cubes.insert(lg_cubes.end(), cubes, cubes+5);
	}
	else if (type == TYPE_STOVE) {
		cube_t body(c), top(c);
		body.z2() = top.z1() = c.z2() - 0.22*c.dz();
		top.d[c.dim][c.dir] -= (c.dir ? 1.0 : -1.0)*0.94*c.get_depth();
		lg_cubes.push_back(body);
		lg_cubes.push_back(top); // or sm_cubes?
	}
	else if (type == TYPE_SHELFRACK) {
		cube_t cubes[9];
		unsigned const num_cubes(get_all_shelf_rack_cubes(c, cubes));

		for (unsigned n = 0; n < num_cubes; ++n) {
			cube_t &C(cubes[n]);
			// spiders can walk on empty shelves, top, and sides, but not shelves with items on them
			if (!c.is_nonempty() || C.z1() == c.z1() || C.z2() == c.z2()) {lg_cubes.push_back(C);}
			else {
				C.expand_in_dim(!c.dim, -0.05*c.dz()); // shrink ends of shelves and back slightly so that they don't interfere with shelf ends
				non_cubes.push_back(C);
			}
		}
	}
	else if (type == TYPE_COUCH) {
		cube_t cubes[4]; // bottom, back, arm, arm
		lg_cubes.insert(lg_cubes.end(), cubes, cubes+get_couch_cubes(c, cubes));
	}
	else if (type == TYPE_CASHREG) {
		cube_t cubes[2]; // body, screen
		get_cashreg_cubes(c, cubes);
		lg_cubes .push_back(cubes[0]); // walk on the body
		non_cubes.push_back(cubes[1]); // avoid the screen
	}
	else if (type == TYPE_CUBICLE) {
		//non_cubes.push_back(c); // using the colliders around the cubicles somewhat works
	}
	else if (type == TYPE_ATTIC_DOOR) {lg_cubes.push_back(get_true_room_obj_bcube(c));}
	// otherwise, treat as a large object; this includes: TYPE_BCASE, TYPE_KSINK (with dishwasher), TYPE_COUCH, TYPE_COLLIDER (cars)
	else {lg_cubes.push_back(c);}
}

// interior collision query used for rats, snakes, and insects: p1 and p2 are line end points; radius applies in X and Y, hheight is half height and applies in +/- z
// note that for_spider is used with insects, not spiders, but it's named this way because it's passed into nested calls that use this variable name
// return value: 0=no coll, 1=d0 wall, 2=d1 wall, 3=closed door d0, 4=closed door d1, 5=open door, 6=stairs, 7=elevator/escalator, 8=exterior wall, 9=room object, 10=closet, 11=cabinet
int building_t::check_line_coll_expand(point const &p1, point const &p2, float radius, float hheight, bool for_spider) const {
	assert(interior != nullptr);
	float const trim_thickness(get_trim_thickness()), zmin(min(p1.z, p2.z));
	float const obj_z1(min(p1.z, p2.z) - hheight), obj_z2(max(p1.z, p2.z) + hheight);
	vector3d const expand(radius, radius, hheight), expand_walls(expand + vector3d(trim_thickness, trim_thickness, 0.0)); // include the wall trim width
	cube_t line_bcube(p1, p2);
	line_bcube.expand_by(expand_walls); // use the larger/conservative expand

	// check interior walls, doors, stairwells, and elevators
	for (unsigned d = 0; d < 2; ++d) {
		if (line_int_cubes_exp(p1, p2, interior->walls[d], expand_walls, line_bcube)) return (d+1); // 2-3
	}
	for (auto const &ds : interior->door_stacks) {
		if (ds.z1() > obj_z2 || ds.z2() < obj_z1) continue; // wrong floor for this stack/part
		// calculate door bounds for bcube test, assuming it's open
		cube_t door_bounds(ds);
		door_bounds.expand_by_xy(ds.get_width());
		if (!line_bcube.intersects(door_bounds) || !line_int_cube_exp(p1, p2, door_bounds, expand)) continue; // optimization
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t const &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() > obj_z2 || door.z2() < obj_z1) continue; // wrong floor

			if (door.open) {
				if (line_int_cube_exp(p1, p2, get_door_bounding_cube(door), expand)) return 5; // approximate - uses bcube
			}
			else if (line_int_cube_exp(p1, p2, door, expand)) return (door.dim + 3); // 3-4
		} // for dix
	} // for door_stacks
	for (auto const &s : interior->stairwells) {
		if (!line_int_cube_exp(p1, p2, s, expand)) continue;
		if (s.shape == SHAPE_WALLED || s.is_u_shape()) return 6; // fully walled and U-shaped stairs always collide, even for spiders and insects
		if (!for_spider && s.z1() < zmin - 0.5f*get_window_vspace()) return 6; // not the ground floor - definitely a collision; but not for spiders or insects

		if (has_room_geom()) { // maybe we're under the stairs; check for individual stairs collisions; this condition should always be true
			for (auto c = interior->room_geom->get_stairs_start(); c != interior->room_geom->objs.end(); ++c) {
				if (c->no_coll()) continue;
				if ((c->type == TYPE_STAIR || (s.shape == SHAPE_WALLED_SIDES && c->type == TYPE_STAIR_WALL)) && line_int_cube_exp(p1, p2, *c, expand)) return 6;
			}
		}
	} // for s
	if (has_mall() && line_int_cube_exp(p1, p2, interior->mall_info->ent_stairs, expand)) return 6; // mall stairs
	if (line_int_cubes_exp(p1, p2, interior->elevators,  expand, line_bcube)) return 7; // collide with entire elevator
	
	for (escalator_t const &e : interior->escalators) {
		if (!line_int_cube_exp(p1, p2, e, expand)) continue;
		cube_t cubes[7];
		e.get_all_cubes(cubes);
		if (line_int_cubes_exp(p1, p2, cubes, 7, expand) || line_int_cube_exp(p1, p2, e.get_ramp_bcube(0), expand)) return 7; // test cubes and ramp bcube
	}
	// check exterior walls
	if (point_in_attic(p1) && point_in_attic(p2)) {} // both points in attic, no need to check exterior walls
	else if (real_num_parts > 1) { // only need to do this for multi-part buildings because the caller is assumed to check the building bcube
		// ray_cast_exterior_walls() doesn't really work here; instead we create a test cube and step it for every radius interval, and check if it ever exits the building
		unsigned const num_steps(min(100U, unsigned(ceil(p2p_dist(p1, p2)/radius))));

		if (num_steps > 1) {
			bool same_part(0);
			auto parts_end(get_real_parts_end_inc_sec());
			cube_t p2_cube(p1, p2);
			p2_cube.expand_by_xy(radius);

			for (auto i = parts.begin(); i != parts_end; ++i) { // includes garage/shed
				if (i->contains_pt(p1)) {same_part = i->contains_cube_xy(p2_cube);}
			}
			if (!same_part && point_in_extended_basement(p1)) {same_part = interior->cube_in_ext_basement_room(p2_cube, 1);} // xy_only=1

			if (!same_part) { // only need to check for going outside the part if p1 and p2 are in different parts (optimization)
				vector3d delta((p2 - p1)/num_steps);
				cube_t test_cube(p1, p1);
				test_cube.expand_by(expand_walls);

				for (unsigned step = 0; step < num_steps-1; ++step) { // take one less step (skip p2)
					test_cube += delta; // take one step first (skip p1)
					if (!is_cube_contained_in_parts(test_cube)) return 8; // extends outside the building walls
				}
			}
		}
	}
	if (has_pool()) { // check swimming pools
		cube_t pool(interior->pool);
		pool.z2() += 2.0*hheight; // raise the top edge of the pool to the height of the caller
		
		if (line_int_cube_exp(p1, p2, pool, expand)) {
			point const p(0.5*(p1 + p2)); // average
			bool const cdim(min(fabs(p.y - pool.y1()), fabs(p.y - pool.y2())) < min(fabs(p.x - pool.x1()), fabs(p.x - pool.x2())));
			return (cdim+1); // 2-3; treat this as a wall collision
		}
	}
	if (!has_room_geom()) return 0; // done (but really shouldn't get here)
	if (line_int_cubes_exp(p1, p2, interior->room_geom->glass_floors, expand_walls, line_bcube)) return 9; // counts as a room object?
	// check room objects and expanded objects (from closets)
	float t(0.0);
	vect_room_object_t::const_iterator b, e;
	bool const use_cached_objs(get_begin_end_room_objs_on_ground_floor(obj_z2, for_spider, b, e));

	for (unsigned vect_id = 0; vect_id < (use_cached_objs ? 1U : 2U); ++vect_id) {
		auto objs_beg((vect_id == 1) ? interior->room_geom->expanded_objs.begin() : b);
		auto objs_end((vect_id == 1) ? interior->room_geom->expanded_objs.end  () : e);

		for (auto c = objs_beg; c != objs_end; ++c) {
			if (c->z1() > obj_z2 || c->z2() < obj_z1) continue; // wrong floor
			// skip non-colliding objects except for balls and books (that the player can drop), computers under desks, and expanded objects from closets,
			// since rats must collide with these
			if (!(for_spider ? c->is_spider_collidable() : c->is_floor_collidable())) continue;
			if (!line_bcube.intersects(*c) || !line_int_cube_exp(p1, p2, get_true_room_obj_bcube(*c), expand)) continue;
			if (c->type == TYPE_RAMP && interior->ignore_ramp_placement && obj_z1 >= ground_floor_z1) continue;

			if (c->is_vert_cylinder()) { // vertical cylinder (should bathroom sink count?)
				if (!is_house && c->type == TYPE_WHEATER) return 9; // office building water heaters have pipes into the floor, more than a cylinder, so use their bcubes
				cylinder_3dw cylin(c->get_cylinder());
				cylin.p1.z -= hheight; cylin.p2.z += hheight; // extend top and bottom
				cylin.r1   += radius ; cylin.r2   += radius;
				if (line_intersect_cylinder_with_t(p1, p2, cylin, (p1.z != p2.z), t)) return 9; // check ends if zvals differ (such as when dropping a rat)
			}
			else if (c->shape == SHAPE_SPHERE) { // sphere (ball)
				float const rsum(radius + c->get_radius());
				if (sphere_test_comp(p1, c->get_cube_center(), (p1 - p2), rsum*rsum, t)) return 9; // approx; uses radius rather than height in z
			}
			else if (c->type == TYPE_CLOSET) {
				cube_t cubes[5];
				get_closet_cubes(*c, cubes, 1); // get cubes for walls and door; for_collision=1
				if (line_int_cubes_exp(p1, p2, cubes, get_closet_num_coll_cubes(*c), expand)) return 10; // closet
			}
			else if (c->type == TYPE_BED) {
				cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
				get_bed_cubes(*c, cubes);
				if (line_int_cube_exp (p1, p2, cubes[0], expand)) return 9; // check bed frame (in case p1.z is high enough)
				get_tc_leg_cubes(cubes[5], *c, BED_HEAD_WIDTH, 0, cubes); // cubes[5] is not overwritten
				if (line_int_cubes_exp(p1, p2, cubes, 4, expand)) return 9; // check legs
			}
			else if (can_use_table_coll(*c)) { // objects with legs
				cube_t cubes[7];
				unsigned const num(get_table_like_object_cubes(*c, cubes));
				if (line_int_cubes_exp(p1, p2, cubes, num, expand)) return 9;
			}
			else if (c->type == TYPE_CONF_TABLE) {
				cube_t cubes[2]; // {top, base}
				get_conf_table_cubes(*c, cubes);
				if (line_int_cubes_exp(p1, p2, cubes, 2, expand)) return 9;
			}
			else if (c->type == TYPE_RDESK) {
				cube_t cubes[3];
				get_reception_desk_cubes(*c, cubes);
				if (line_int_cubes_exp(p1, p2, cubes, 3, expand)) return 9;
			}
			else if (c->type == TYPE_CHAIR) {
				cube_t cubes[3], leg_cubes[4]; // seat, back, legs_bcube
				get_chair_cubes(*c, cubes);
				if (line_int_cube_exp(p1, p2, cubes[0], expand)) return 9; // check seat
				get_tc_leg_cubes(cubes[2], *c, c->get_chair_leg_width(), 1, leg_cubes);
				if (line_int_cubes_exp(p1, p2, leg_cubes, 4, expand)) return 9; // check legs
			}
			else if (c->type == TYPE_BCASE && !c->is_on_floor()) { // treat as a cube if fallen over
				cube_t top, middle, back, lr[2];
				get_bookcase_cubes(*c, top, middle, back, lr);
				cube_t const all_cubes[5] = {top, middle, back, lr[0], lr[1]}; // probably don't need to check the top and back, but okay to do so
				if (line_int_cubes_exp(p1, p2, all_cubes, 5, expand)) return 9;
			}
			else if (c->type == TYPE_KSINK) {
				cube_t cubes[3];
				if (line_int_cubes_exp(p1, p2, cubes, get_ksink_cubes(*c, cubes), expand)) return 11; // 1 or 3 sink parts; counts as a cabinet
			}
			else if (c->type == TYPE_COUCH) {
				cube_t couch_body(*c);
				couch_body.z1() += 0.06*c->dz(); // there's space under the couch
				if (line_int_cube_exp(p1, p2, couch_body, expand)) return 9;
			}
			else if (c->type == TYPE_POOL_TABLE) {
				cube_t cubes[5];
				get_pool_table_cubes(*c, cubes);
				if (line_int_cubes_exp(p1, p2, cubes, 5, expand)) return 11; // 1 or 3 sink parts; counts as a cabinet
			}
			else if (c->type == TYPE_SHELVES) {
				if (c->in_mall()) { // mall shelves have fishtanks (pet stores) and flat items (clothing stores)
					cube_t shelves[4]; // max number of shelves
					unsigned const num_shelves(get_shelves_for_object(*c, shelves));
					if (line_int_cubes_exp(p1, p2, shelves, num_shelves, expand)) return 9;
				}
				else { // storage room shelves have all sorts of items when expanded, which we don't test
					if (line_int_cube_exp(p1, p2, get_shelves_no_bot_gap(*c), expand)) return 9; // there's space under the shelves
				}
			}
			else if (c->type == TYPE_SHELFRACK) {
				cube_t cubes[9];
				unsigned const num(get_all_shelf_rack_cubes(*c, cubes));
				if (line_int_cubes_exp(p1, p2, cubes, num, expand)) return 9;
			}
			else if (c->type == TYPE_BENCH) {
				cube_t cubes[4]; // seat, lo side, hi side, [back]
				unsigned const num(get_bench_cubes(*c, cubes));
				if (line_int_cubes_exp(p1, p2, cubes, num, expand)) return 9;
			}
			else if (c->is_parked_car()) { // parked car
				cube_t cubes[5];
				get_approx_car_cubes(*c, cubes);
				if (line_int_cubes_exp(p1, p2, cubes, 5, expand)) return 9;
			}
			//else if (c->type == TYPE_STALL && maybe_inside_room_object(*c, p2, radius)) {} // is this useful? inside test only applied to end point
			else {
				return ((c->type == TYPE_COUNTER || c->type == TYPE_CABINET) ? 11 : 9); // intersection
			}
		} // for c
	} // for vect_id
	return 0;
}

// visibility query used for rats: ignores room objects
bool building_t::check_line_of_sight_large_objs(point const &p1, point const &p2) const {
	assert(interior != nullptr);
	cube_t const line_bcube(p1, p2);
	if (line_int_cubes(p1, p2, interior->fc_occluders, line_bcube)) return 0; // likely fastest test

	for (unsigned d = 0; d < 2; ++d) {
		if (line_int_cubes(p1, p2, interior->walls[d], line_bcube)) return 0;
	}
	for (auto const &ds : interior->door_stacks) {
		cube_t const ds_bcube(ds.get_true_bcube()); // expand to nonzero area
		if (!line_bcube.intersects(ds_bcube) || !ds_bcube.line_intersects(p1, p2)) continue; // no intersection with this stack
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t const &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			cube_t const door_bcube(door.get_true_bcube()); // expand to nonzero area
			if (!door.open && line_bcube.intersects(door_bcube) && door_bcube.line_intersects(p1, p2)) return 0;
		}
	} // for door_stacks
	if (line_int_cubes(p1, p2, interior->elevators, line_bcube)) return 0; // check elevators only because stairs and escalators aren't solid visual blockers

	// check exterior walls if there are multiple non-basement parts; skip for extended basement because this logic doesn't work and only interior walls are needed;
	// this is simpler than ray_cast_exterior_walls()
	if (real_num_parts > (1 + has_basement()) && !point_in_extended_basement(p1) && !point_in_extended_basement(p2)) {
		if (point_in_attic(p1) && point_in_attic(p2)) return 1; // always visible in attic
		float tot_len(0.0);
		auto parts_end(get_real_parts_end_inc_sec());
		bool contained(0);

		for (auto i = parts.begin(); i != parts_end; ++i) { // includes garage/shed
			if (i->contains_pt(p1) && i->contains_pt(p2)) {contained = 1; break;} // fully contained optimization
			point p1c(p1), p2c(p2);
			if (do_line_clip(p1c, p2c, i->d)) {tot_len += p2p_dist(p1c, p2c);}
		}
		// length of line segments clipped to parts < total line length => line not contained in sum of parts
		if (!contained && tot_len < 0.99*p2p_dist(p1, p2)) return 0;
	}
	return 1;
}

bool building_t::check_line_int_interior_window(point const &p1, point const &p2) const {
	if (!interior || interior->int_windows.empty()) return 0;
	return line_int_cubes(p1, p2, interior->int_windows, cube_t(p1, p2));
}

cube_t building_t::get_best_occluder(point const &camera_bs) const {
	cube_t ret;
	if (!interior) return ret;
	point const ray_end(camera_bs + 4.0*get_window_vspace()*cview_dir); // look up to 4 floor heights ahead
	float best_tmin(1.0);
	
	// check for walls in front of the player
	for (unsigned d = 0; d < 2; ++d) {
		for (cube_t const &w : interior->walls[d]) {
			float tmin(0.0), tmax(1.0);
			if (get_line_clip(camera_bs, ray_end, w.d, tmin, tmax) && tmin < best_tmin) {best_tmin = tmin; ret = w;}
		}
	} // for d
	if (!ret.is_all_zeros()) return ret;
	// if no good walls, find best floor occluder, or ceiling if on the ground floor
	float const above_penalty(get_window_vspace()); // prefer floors below
	float z_best(0.0);

	for (cube_t const &c : interior->fc_occluders) {
		if (!c.contains_pt_xy(camera_bs) || c.z1() < ground_floor_z1) continue; // skip basement cubes
		float const z(c.zc());
		float z_dist;
		if (z < camera_bs.z) {z_dist = camera_bs.z - z;} // floor below
		else {z_dist = z - camera_bs.z + above_penalty;} // floor above
		if (ret.is_all_zeros() || z_dist < z_best) {z_best = z_dist; ret = c;}
	} // for c
	return ret;
}

bool handle_vcylin_vcylin_int(point &p1, point const &p2, float rsum) {
	if (!dist_xy_less_than(p1, p2, rsum)) return 0; // no collision
	vector3d const delta(p1.x-p2.x, p1.y-p2.y, 0.0); // ignore zvals
	float const xy_dist(delta.mag());
	p1 += delta*((rsum - xy_dist)/xy_dist);
	return 1;
}
bool handle_dynamic_room_objs_coll(vect_room_object_t::const_iterator begin, vect_room_object_t::const_iterator end, point &pos, float radius, float z2) {
	for (auto c = begin; c != end; ++c) {
		if (c->no_coll() || !c->has_dstate()) continue; // Note: no test of player_coll flag
		if (pos.z > c->z2() || z2 < c->z1())  continue; // different floors
		assert(is_ball_type(c->type)); // currently, only balls have has_dstate()
		// treat the ball as a vertical cylinder because it's too complex to find the collision point of a vertical cylinder with a sphere
		point const center(c->get_cube_center());
		if (handle_vcylin_vcylin_int(pos, center, (radius + c->get_radius()))) return 1; // early exit on first hit
	} // for c
	return 0;
}

// Note: cur_obj_pos is the center of the snake while pos is the center of its head
template<typename T> void vect_animal_t<T>::update_delta_sum_for_animal_coll(point const &pos, point const &cur_obj_pos,
	float radius, float z1, float z2, float radius_scale, float &max_overlap, vector3d &delta_sum) const
{
	float const rsum_max(radius_scale*(radius + max_radius) + max_xmove), coll_x1(pos.x - rsum_max), coll_x2(pos.x + rsum_max);
	auto start(get_first_with_xv_gt(coll_x1)); // use a binary search to speed up iteration

	for (auto r = start; r != this->end(); ++r) {
		if (r->pos.x > coll_x2) break; // none after this can overlap - done
		if (r->pos == cur_obj_pos) continue; // skip ourself
		float coll_radius(r->get_xy_radius());
		if (!dist_xy_less_than(pos, r->pos, radius_scale*(radius + coll_radius))) continue; // no collision
		if (z2 < r->pos.z || z1 > (r->pos.z + r->get_height())) continue; // different floors; less likely to reject, so done last
		point r_hit_pos(r->pos);
		if (!r->detailed_sphere_coll(pos, radius, r_hit_pos, coll_radius)) continue; // check detailed collision for snakes
		float const overlap(radius_scale*(radius + coll_radius) - p2p_dist_xy(pos, r_hit_pos));
		vector3d const delta((pos - point(r_hit_pos.x, r_hit_pos.y, pos.z)).get_norm()*overlap);
		delta_sum += delta; // accumulate weighted delta across collisions
		max_eq(max_overlap, overlap);
	} // for r
}

// vertical cylinder collision detection with dynamic objects: balls, the player? people? other rats?
// only handles the first collision
bool building_t::check_and_handle_dynamic_obj_coll(point &pos, point const &cur_obj_pos, float radius,
	float z1, float z2, point const &camera_bs, bool for_spider, bool skip_player) const
{
	if (camera_surf_collide && !player_wait_respawn && !skip_player) { // check the player
		if (z1 < camera_bs.z && z2 > (camera_bs.z - get_bldg_player_height() - 0.1*CAMERA_RADIUS)) { // slighly below the player's feet
			if (handle_vcylin_vcylin_int(pos, camera_bs, (radius + get_scaled_player_radius()))) return 1;
		}
	}
	// check for collisions of people, at least for snakes (since rats will generally avoid people)
	for (person_t const &p : interior->people) {
		float const pz1(p.get_z1()), pz2(p.get_z2());
		if (z1 < pz2 && z2 > pz1 && handle_vcylin_vcylin_int(pos, p.pos, (radius + p.get_width()))) return 1;
	}
	assert(has_room_geom());
	// check dynamic objects such as balls; currently, we only check this for houses because they have balls on the floor;
	// in theory, the player can put a ball in an office building, but we don't handle that case because office buildings have tons of objects and this is too slow
	if (is_house) {
		vect_room_object_t::const_iterator b, e;
		get_begin_end_room_objs_on_ground_floor(z2, for_spider, b, e);
		handle_dynamic_room_objs_coll(b, e, pos, radius, z2);
	}
	float max_overlap(0.0);
	vector3d delta_sum;
	float const rat_radius_scale    = 0.7; // allow them to get a bit closer together, since radius is conservative
	float const spider_radius_scale = 0.75;
	float const snake_radius_scale  = 1.0;
	// no collisions with sewer_rats, pet_rats, sewer_spiders, or pet_snakes
	interior->room_geom->rats   .update_delta_sum_for_animal_coll(pos, cur_obj_pos, radius, z1, z2, rat_radius_scale,    max_overlap, delta_sum);
	interior->room_geom->spiders.update_delta_sum_for_animal_coll(pos, cur_obj_pos, radius, z1, z2, spider_radius_scale, max_overlap, delta_sum); // should we avoid squished spiders?
	interior->room_geom->snakes .update_delta_sum_for_animal_coll(pos, cur_obj_pos, radius, z1, z2, snake_radius_scale,  max_overlap, delta_sum);
	
	if (max_overlap > 0.0) { // we have at least one collision
		float const delta_mag(delta_sum.mag());
		if (delta_mag > max_overlap) {delta_sum *= max_overlap/delta_mag;} // clamp to max_overlap
		pos += 0.5*delta_sum; // dampen by 0.5 for stability
		return 1;
	}
	return 0;
}

void add_cube_int_volume(cube_t const &a, cube_t const &b, float &cont_vol) {
	if (a.intersects(b)) {cont_vol += (min(a.x2(), b.x2()) - max(a.x1(), b.x1()))*(min(a.y2(), b.y2()) - max(a.y1(), b.y1()))*(min(a.z2(), b.z2()) - max(a.z1(), b.z1()));}
}
bool building_t::is_cube_contained_in_parts(cube_t const &c) const {
	float cont_vol(0); // total shared volume
	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {add_cube_int_volume(c, *p, cont_vol);}
	
	if (has_ext_basement() && interior->basement_ext_bcube.intersects(c)) { // extended basement counts, even though it's not a part
		for (auto r = interior->ext_basement_rooms_start(); r != interior->rooms.end(); ++r) {add_cube_int_volume(c, *r, cont_vol);}
	}
	return (cont_vol > 0.99*c.get_volume()); // add a bit of tolerance
}

// returns: 0=not in basement, 1=on basement stairs, 2=fully in basement, 3=in extended basement
int building_t::check_player_in_basement(point const &pos) const {
	point const pos_rot(get_inv_rot_pos(pos)); // rotate pos into building space
	if (!is_pos_in_basement(pos_rot))                     return 0;
	if (point_in_extended_basement_not_basement(pos_rot)) return 3;
	
	if (interior && (pos_rot.z - get_bldg_player_height()) > (ground_floor_z1 - get_window_vspace())) { // only need to check if on the top floor of the basement
		if (point_in_stairwell(pos_rot)) return 1; // player on stairs, upper floor and windows/outside may be visible
		if (has_pg_ramp() && interior->pg_ramp.contains_pt(pos_rot)) return 1;
	}
	return 2; // player in basement but not on stairs or ramp
}

bool building_t::point_on_basement_stairs(point const &pos_bs) const {
	if (interior == nullptr) return 0;
	if (fabs(pos_bs.z - ground_floor_z1) > get_window_vspace()) return 0; // above first floor or below first floor of basement

	for (stairwell_t const &s : interior->stairwells) {
		if (s.z1() > ground_floor_z1 || s.z2() < ground_floor_z1) continue; // not extending to basement
		if (s.contains_pt(pos_bs)) return 1;
	}
	return 0;
}

bool building_t::is_obj_above_ramp(cube_t const &c) const {
	assert(interior);
	if (!has_pg_ramp()) return 0;
	return (interior->pg_ramp.contains_pt_xy(c.get_cube_center()) && c.z2() >= interior->pg_ramp.z2() && (c.z1() - interior->pg_ramp.z2()) <= get_window_vspace());
}
bool building_t::is_room_above_ramp(cube_t const &room, float zval) const {
	assert(interior);
	if (!has_pg_ramp()) return 0;
	return (interior->pg_ramp.intersects_xy(room) && zval >= interior->pg_ramp.z2() && (zval - interior->pg_ramp.z2()) <= get_window_vspace());
}
bool building_t::is_room_adjacent_to_ext_door(cube_t const &room, bool front_door_only) const {
	cube_t room_exp(room);
	room_exp.expand_by_xy(get_wall_thickness());

	for (tquad_with_ix_t const &d : doors) { // exterior doors
		if (!d.is_exterior_door() || d.type == tquad_with_ix_t::TYPE_RDOOR) continue;
		if (room_exp.contains_pt(d.get_bcube().get_cube_center())) return 1;
		if (front_door_only) return 0; // assumes the first door is the front door
	}
	return 0;
}
bool building_t::cube_int_ext_door(cube_t const &c) const {
	for (tquad_with_ix_t const &d : doors) { // exterior doors
		if (!d.is_exterior_door() || d.type == tquad_with_ix_t::TYPE_RDOOR) continue;
		if (c.intersects(d.get_bcube())) return 1;
	}
	return 0;
}

// for use with in dir lighting updates when interior doors are opened or closed
void building_t::get_rooms_for_door(unsigned door_ix, int room_ix[2]) const {
	door_t const &door(get_door(door_ix));
	cube_t dbc(door.get_true_bcube());
	dbc.expand_in_dim(door.dim, get_wall_thickness()); // make sure it overlaps the room
	room_ix[0] = room_ix[1] = -1; // reset to invalid rooms in case no rooms are found

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (!r->intersects(dbc)) continue;
		bool side(r->get_center_dim(door.dim) < door.get_center_dim(door.dim));
		
		if (room_ix[side] >= 0) { // must be basement door with room on top and bottom rather than each side
			assert(door.on_stairs);
			if (room_ix[!side] < 0) {side ^= 1;} // use the other side if unassigned
		}
		room_ix[side] = (r - interior->rooms.begin());
	} // for r
}
void building_t::get_lights_for_room_and_floor(unsigned room_ix, unsigned floor_ix, vector<unsigned> &light_ids) const {
	assert(has_room_geom());
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_type() || i->room_id != room_ix || get_floor_for_zval(i->z1()) != floor_ix) continue; // not a light, wrong room, or wrong floor
		light_ids.push_back(i - interior->room_geom->objs.begin());
	}
}
void building_t::get_lights_near_door(unsigned door_ix, vector<unsigned> &light_ids) const {
	light_ids.clear();
	int room_ix[2] = {};
	get_rooms_for_door(door_ix, room_ix);
	unsigned const floor_ix(get_floor_for_zval(get_door(door_ix).get_center_dim(2)));

	for (unsigned d = 0; d < 2; ++d) {
		if (room_ix[d] >= 0) {get_lights_for_room_and_floor((unsigned)room_ix[d], floor_ix, light_ids);}
	}
}

// it's up to the caller to check that the door is open or closed (depending on the query)
bool building_t::is_cube_visible_through_door(point const &viewer, cube_t const &c, door_t const &door) const {
	cube_t const door_bcube(door.get_true_bcube()); // must be nonzero area for this test to be correct
	if (door_bcube.contains_pt(viewer)) return 1; // viewer in doorway
	float const door_pos(door.get_center_dim(door.dim)), vpos(viewer[door.dim]);
	if ((c.get_center_dim(door.dim) < door_pos) == (vpos < door_pos)) return 0; // viewer and cube on same side of door - not visible through it (optimization)
	point pts[8];
	unsigned const npts(get_cube_corners(c.d, pts, viewer, 0)); // should return only the 6 visible corners
	cube_t proj_area;

	for (unsigned n = 0; n < npts; ++n) {
		float const t((door_pos - vpos)/(pts[n][door.dim] - vpos));
		point const proj(viewer + t*(pts[n] - viewer));
		proj_area.assign_or_union_with_pt(proj);
	}
	return door_bcube.intersects(proj_area);
}

void building_t::print_building_manifest() const { // Note: skips expanded_objs
	cout << TXT(parts.size()) << TXT(doors.size()) << TXT(details.size());

	if (interior) {
		unsigned const walls(interior->walls[0].size() + interior->walls[1].size()), rooms(interior->rooms.size());
		unsigned const floors(interior->floors.size()), ceilings(interior->ceilings.size()), door_stacks(interior->door_stacks.size());
		unsigned const doors(interior->doors.size()), stairs(interior->stairwells.size()), elevators(interior->elevators.size()), escalators(interior->escalators.size());
		cout << TXT(walls) << TXT(rooms) << TXT(floors) << TXT(ceilings) << TXT(door_stacks) << TXT(doors) << TXT(stairs) << TXT(elevators) << TXT(escalators);
	}
	if (has_room_geom()) {
		building_room_geom_t const &g(*interior->room_geom);
		unsigned const objects(g.objs.size()), models(g.obj_model_insts.size()), snakes(g.snakes.size());
		unsigned const rats(g.rats.size()), sewer_rats(g.sewer_rats.size()), pet_rats(g.pet_rats.size()), pet_birds(g.pet_birds.size());
		unsigned const spiders(g.spiders.size()), sewer_spiders(g.sewer_spiders.size()), pet_snakes(g.pet_snakes.size());
		cout << TXT(objects) << TXT(models) << TXT(rats) << TXT(sewer_rats) << TXT(pet_rats)
			 << TXT(spiders) << TXT(sewer_spiders) << TXT(snakes) << TXT(pet_snakes) << TXT(pet_birds) << endl;
		unsigned obj_counts[NUM_ROBJ_TYPES] = {};

		for (auto const &i : g.objs) {
			assert(i.type < NUM_ROBJ_TYPES);
			++obj_counts[i.type];
		}
		for (unsigned n = 0; n < NUM_ROBJ_TYPES; ++n) {
			if (obj_counts[n] == 0) continue;
			std::string const &name(bldg_obj_types[n].name);
			cout << name << ":" << std::string((18 - name.size()), ' ') << obj_counts[n] << endl;
		}
	}
	else {cout << endl;}
}
void building_t::print_building_stats() const {
	if (!interior) return;
	float const floor_thickness(get_floor_thickness()), floor_spacing(get_window_vspace());
	float const feet_per_unit(8.0/floor_spacing), feet_per_unit_sq(feet_per_unit*feet_per_unit);
	unsigned num_rooms(0), num_bedrooms(0), num_bathrooms(0);
	float square_footage(0.0);

	for (room_t const &room : interior->rooms) {
		if (room.is_sec_bldg)       continue; // garage/shed not included in count
		if (room.is_ext_basement()) continue; // extended basement not included in count
		if (room.is_hallway)        continue; // hallway not included in count
		unsigned const num_floors(calc_num_floors_room(room, floor_spacing, floor_thickness));
		assert(num_floors > 0);
		num_rooms += num_floors;
		if (room.z1() >= ground_floor_z1) {square_footage += num_floors*feet_per_unit_sq*room.get_area_xy();} // basement doesn't count

		for (unsigned f = 0; f < num_floors; ++f) {
			room_type const rtype(room.get_room_type(f));
			if (rtype == RTYPE_BED || rtype == RTYPE_MASTER_BED) {++num_bedrooms;}
			else if (is_bathroom(rtype)) {++num_bathrooms;}
		}
	} // for room
	cout << TXT(num_rooms) << TXT(num_bedrooms) << TXT(num_bathrooms) << TXT(square_footage) << endl;
}

