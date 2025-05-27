// 3D World - City Terrain and Road System
// by Frank Gennari
// 08/22/21

#include "city_terrain.h"

extern float water_plane_z;
extern city_params_t city_params;


// heightmap_query_t

float smooth_interp(float a, float b, float mix) {
	mix = mix * mix * (3.0 - 2.0 * mix); // cubic Hermite interoplation (smoothstep)
	return mix*a + (1.0 - mix)*b;
}

cube_t heightmap_query_t::get_cube_for_bounds(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float elevation) const {
	return cube_t(get_x_value(x1), get_x_value(x2), get_y_value(y1), get_y_value(y2), elevation, elevation);
}
cube_t heightmap_query_t::get_cube_for_cell(int x, int y) const {
	// Note: cube is actually twice the grid size in X/Y because the point can round in either direction and we need to guarantee any point that maps to (x,y) is contained by this cube
	float const xv(get_x_value(x)), yv(get_y_value(y)), z(get_height_clamped(x, y));
	return cube_t(xv-DX_VAL, xv+DX_VAL, yv-DY_VAL, yv+DY_VAL, z, z);
}

bool heightmap_query_t::any_underwater(unsigned x1, unsigned y1, unsigned x2, unsigned y2, bool check_border) const {
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

void heightmap_query_t::get_segment_end_pts(road_t const &r, unsigned six, unsigned eix, point &ps, point &pe) const {
	float const sv(r.dim ? get_y_value(six) : get_x_value(six)), ev(r.dim ? get_y_value(eix) : get_x_value(eix)); // start/end pos of bridge in road dim
	float const z1(r.get_start_z()), z2(r.get_end_z()), dz(z2 - z1), len(r.get_length()), v0(r.dim ? r.y1() : r.x1());
	ps[ r.dim] = sv;
	pe[ r.dim] = ev;
	ps[!r.dim] = pe[!r.dim] = r.get_center_dim(!r.dim);
	ps.z = z1 + dz*CLIP_TO_01((sv - v0)/len);
	pe.z = z1 + dz*CLIP_TO_01((ev - v0)/len);
}

void heightmap_query_t::flatten_region_to(cube_t const c, unsigned slope_width, bool decrease_only) {
	flatten_region_to(get_x_pos(c.x1()), get_y_pos(c.y1()), get_x_pos(c.x2()), get_y_pos(c.y2()), slope_width, (c.z1() - ROAD_HEIGHT), decrease_only);
}
void heightmap_query_t::flatten_region_to(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float elevation, bool decrease_only) {
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

float heightmap_query_t::flatten_sloped_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float z1, float z2, bool dim, unsigned border,
	unsigned skip_six, unsigned skip_eix, bool stats_only, bool decrease_only, bridge_t *bridge, tunnel_t *tunnel)
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

float heightmap_query_t::flatten_for_road(road_t const &road, unsigned border, bool stats_only, bool decrease_only, bridge_t *bridge, tunnel_t *tunnel) {
	float const z_adj(road.get_z_adj());
	unsigned const rx1(get_x_pos(road.x1())), ry1(get_y_pos(road.y1())), rx2(get_x_pos(road.x2())), ry2(get_y_pos(road.y2()));
	return flatten_sloped_region(rx1, ry1, rx2, ry2, road.d[2][road.slope]-z_adj, road.d[2][!road.slope]-z_adj, road.dim, border, 0, 0, stats_only, decrease_only, bridge, tunnel);
}


// city connector road path finding

bool city_road_connector_t::get_closer_dir(cube_t const &A, cube_t const &B, bool dim) {
	float const center(B.get_center_dim(dim));
	return (abs(A.d[dim][1] - center) < abs(A.d[dim][0] - center));
}
cube_t get_road_between_pts(point const &p1, point const &p2, float road_hwidth, bool expand_start, bool expand_end) {
	assert((p1.x == p2.x) != (p1.y == p2.y)); // roads must either be in X or Y
	bool const dim(p1.x == p2.x), dir(p1[dim] < p2[dim]);
	cube_t road(p1, p2);
	road.expand_in_dim(!dim, road_hwidth);
	if (expand_start) {road.d[dim][!dir] -= (dir ? 1.0 : -1.0)*road_hwidth;}
	if (expand_end  ) {road.d[dim][ dir] += (dir ? 1.0 : -1.0)*road_hwidth;}
	return road;
}
bool road_seg_valid(point const &p1, point const &p2, bool dim, vect_cube_t const &blockers, float road_hwidth, bool expand_start, bool expand_end) {
	float const length(fabs(p1[dim] - p2[dim]));
	if (length < 4.0*road_hwidth) return 0; // too short
	if (fabs(p1.z - p2.z)/(length - road_hwidth) > city_params.max_road_slope) return 0; // check slope
	return !has_bcube_int_xy_no_adj(get_road_between_pts(p1, p2, road_hwidth, expand_start, expand_end), blockers);
}
bool check_pt_valid(point const &pt, cube_t const exclude[2]) {
	return (pt.z > water_plane_z && !exclude[0].contains_pt_xy(pt) && !exclude[1].contains_pt_xy(pt));
}

float city_road_connector_t::calc_road_cost(point const &p1, point const &p2) {
	bool const dim(fabs(p2.x - p1.x) < fabs(p2.y - p1.y)), dir(p1[dim] < p2[dim]), slope((p1.z < p2.z) ^ dir);
	road_t const road(p1, p2, city_params.road_width, dim, slope, 0); // road_ix = 0
	if (!segment_road(road, 1)) return 0.0; // check_only=1
	float tot_dz(0.0);

	for (auto s = segments.begin(); s != segments.end(); ++s) {
		if (s->z2() < s->z1()) {swap(s->z2(), s->z1());} // swap zvals if needed
		assert(s->is_normalized());
		tot_dz += hq.flatten_for_road(*s, city_params.road_border, 1, 0, nullptr, nullptr); // check_only=1, no bridges or tunnels
	}
	return tot_dz;
}
float city_road_connector_t::calc_road_path_cost(vector<point> &pts) {
	if (pts.empty()) return 0.0; // failed case
	float cost(0.0);

	for (auto p = pts.begin(); p+1 != pts.end(); ++p) {
		float const seg_cost(calc_road_cost(*p, *(p+1)));
		if (seg_cost == 0.0) {pts.clear(); return 0.0;} // failed
		cost += seg_cost;
	}
	cost *= (pts.size() - 1); // add jog penalty
	return cost;
}

float city_road_connector_t::find_route_between_points(point const &p1, point const &p2, vect_cube_t const &blockers, vector<point> &pts,
	cube_t const &bcube1, cube_t const &bcube2, float road_hwidth, bool dim1, bool dir1, bool dim2, bool dir2)
{
	float const min_extend(4.0*road_hwidth), min_jog(4.0*road_hwidth);
	// can't route if the two endpoints are too close to add a jog/road
	if (min(fabs(p1.x - p2.x), fabs(p1.y - p2.y)) <= max(min_extend, 2.0f*HALF_DXY)) return 0.0;
	cube_t exclude[2] = {bcube1, bcube2};
	for (unsigned d = 0; d < 2; ++d) {exclude[d].expand_by_xy(min_jog);}

	if (dim1 == dim2) { // can add 2 points to create a jog
		assert(dir1 != dir2); // must be opposing
		float const gap_lo(min(p1[dim1], p2[dim1])), gap_hi(max(p1[dim1], p2[dim1]));

		if ((gap_hi - gap_lo) > 2.1f*min_jog) { // we have enough space for a jog
			point pt[2];
			pt[0][!dim1] = p1[!dim1]; pt[1][!dim1] = p2[!dim1];
			vector<point> cand_pts;
			float min_cost(0.0);

			for (unsigned n = 0; n < 10; ++n) { // choose best of 10 attempts
				pt[0][dim1] = pt[1][dim1] = rgen.rand_uniform((gap_lo + min_jog), (gap_hi - min_jog)); // choose a random jog pos
				for (unsigned d = 0; d < 2; ++d) {pt[d].z = hq.get_road_zval_at_pt(pt[d]);}
		
				if (check_pt_valid(pt[0], exclude) && check_pt_valid(pt[1], exclude) && road_seg_valid(p1, pt[0], dim1, blockers, road_hwidth, 0, 1) &&
					road_seg_valid(pt[0], pt[1], !dim1, blockers, road_hwidth, 1, 1) && road_seg_valid(pt[1], p2, dim1, blockers, road_hwidth, 1, 0))
				{
					cand_pts.clear();
					cand_pts.push_back(p1);
					for (unsigned d = 0; d < 2; ++d) {cand_pts.push_back(pt[d]);}
					cand_pts.push_back(p2);
					float const cost(calc_road_path_cost(cand_pts));
					if (cost > 0.0 && (min_cost == 0.0 || cost < min_cost)) {min_cost = cost; pts = cand_pts;}
				}
			} // for n
			return min_cost;
		}
	}
	// add one point to create a right angle bend
	point ipt;
	ipt[ dim1] = p2[ dim1];
	ipt[!dim1] = p1[!dim1];
	ipt.z = hq.get_road_zval_at_pt(ipt);

	if (check_pt_valid(ipt, exclude) && road_seg_valid(p1, ipt, dim1, blockers, road_hwidth, 0, 1) && road_seg_valid(ipt, p2, dim2, blockers, road_hwidth, 1, 0)) {
		pts.push_back(p1);
		pts.push_back(ipt);
		pts.push_back(p2);
		return calc_road_path_cost(pts); // Note: clears pts on failure
	}
	return 0.0; // failed
}

bool city_road_connector_t::segment_road(road_t const &road, bool check_only) {
	bool const dim(road.dim);
	float const road_len(road.get_length()), conn_pos(road.get_center_dim(!dim));
	assert(road_len > 0.0);
	unsigned const num_segs(ceil(road_len/city_params.conn_road_seg_len));
	assert(num_segs > 0 && num_segs < 1000); // sanity check
	segments.clear();
	if (num_segs == 1) {segments.push_back(road); return (fabs(road.get_slope_val()) < city_params.max_road_slope);} // single segment optimization
	float const seg_len(road_len/num_segs);
	float const min_zval(water_plane_z + 0.3*road.get_width()); // water_plane_z doesn't take into account tessellated waves; shift further up
	assert(seg_len <= city_params.conn_road_seg_len);
	road_t rs(road); // keep d[!dim][0], d[!dim][1], dim, and road_ix
	rs.z1() = road.d[2][road.slope];

	for (unsigned n = 0; n < num_segs; ++n) {
		rs.d[dim][1] = ((n+1 == num_segs) ? road.d[dim][1] : (rs.d[dim][0] + seg_len)); // make sure it ends exactly at the correct location
		point pos;
		pos[ dim] = rs.d[dim][1];
		pos[!dim] = conn_pos;
		rs.z2()   = max(hq.get_road_zval_at_pt(pos), min_zval); // terrain height at end of segment; don't go below water level
		rs.slope  = (rs.z2() < rs.z1());

		if (fabs(rs.get_slope_val()) > city_params.max_road_slope) { // slope is too high, clamp z2 to max allowed value
			if (n+1 == num_segs) {
				// Note: the height of the first/last segment may change slightly after placing the bend,
				// which can make this slope check fail when check_only=0 while it passed when check_only=1;
				// returning here will create a disconnected road segment and fail an assert later, so instead we allow the high slope
				if (check_only) return 0;
			}
			else {rs.z2() = rs.z1() + seg_len*city_params.max_road_slope*SIGN(rs.dz());}
		}
		segments.push_back(rs);
		rs.d[dim][0] = rs.d[dim][1]; rs.z1() = rs.z2(); // shift segment end point
	} // for n
	return 1;
}

// not roads, but uses similar routing logic; maybe this class should be renamed
bool city_road_connector_t::is_tline_seg_valid(point const &p1, point const &p2, float max_ground_clearance) const {
	vector3d const delta(p2 - p1);
	float const step_len(delta.xy_mag()), hmap_check_step(HALF_DXY);
	unsigned const num_hmap_steps(round_fp(step_len/hmap_check_step));
	vector3d const hmap_step_delta(delta/(num_hmap_steps + 1U)); // calculate here so that zvals are correct
	point test_pos(p1 + hmap_step_delta); // take one hmap step

	for (unsigned s = 0; s < num_hmap_steps; ++s) {
		float const hm_zval(hq.get_height_at(test_pos));
		if (test_pos.z < hm_zval + max_ground_clearance) return 0; // wire comes too close to terrain, fail
		test_pos += hmap_step_delta;
	}
	return 1;
}
bool check_and_add_tower_pt(point const &pos, float height, float clearance_radius, transmission_line_t &tline, vect_cube_t &blockers) {
	// Note: this doesn't check for two transmission line wires intersecting each other, but this should be unlikely, and I've never seen it happen
	cube_t tower_area(pos);
	tower_area.z1() -= height; // set base of tower
	if (tower_area.z1() < water_plane_z) return 0; // underwater, invalid tower location
	tower_area.expand_by_xy(clearance_radius);
	if (has_bcube_int_xy_no_adj(tower_area, blockers)) return 0;
	blockers.push_back(tower_area);
	tline.tower_pts.push_back(pos);
	return 1;
}
bool city_road_connector_t::route_transmission_line(transmission_line_t &tline, vect_cube_t &blockers, float road_width, float road_spacing) const {
	vector3d const vxy(tline.p2.x-tline.p1.x, tline.p2.y-tline.p1.y, 0.0);
	float const dist(vxy.xy_mag()), tower_spacing(1.0*road_spacing), max_ground_clearance(0.5*tline.tower_height), clearance_radius(0.5*road_width);
	unsigned const num_towers(2U + unsigned(floor(dist/tower_spacing))); // includes towers at the two end points
	vector3d const step_delta(vxy/(num_towers - 1U)); // divide distance by the number of spans
	point cur_pos(tline.p1); // first tower XY location
	cur_pos.z = hq.get_height_at(cur_pos) + tline.tower_height;
	if (!check_and_add_tower_pt(cur_pos, tline.tower_height, clearance_radius, tline, blockers)) return 0; // first tower point

	for (unsigned n = 0; n < num_towers-1; ++n) { // one tower was already added
		point next_pos(cur_pos + step_delta);
		next_pos.z = hq.get_height_at(next_pos) + tline.tower_height;

		if (!is_tline_seg_valid(cur_pos, next_pos, max_ground_clearance)) {
			// not valid, try adding an extra tower in the middle; if it's still invalid then fail
			point midpoint(0.5*(cur_pos + next_pos));
			midpoint.z = hq.get_height_at(midpoint) + tline.tower_height;
			if (!check_and_add_tower_pt(midpoint, tline.tower_height, clearance_radius, tline, blockers)) return 0;
			if (!is_tline_seg_valid(cur_pos, midpoint, max_ground_clearance) || !is_tline_seg_valid(midpoint, next_pos, max_ground_clearance)) return 0;
		}
		cur_pos = next_pos;
		if (!check_and_add_tower_pt(cur_pos, tline.tower_height, clearance_radius, tline, blockers)) return 0;
	} // for n
	return 1;
}

