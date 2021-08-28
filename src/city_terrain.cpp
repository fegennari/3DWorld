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
	float const xv(get_x_value(x)), yv(get_y_value(y)), z(get_height_clamped(x, y));
	return cube_t(xv, xv+DX_VAL, yv, yv+DY_VAL, z, z);
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

class road_router_t {
	struct xy_ix_t { // should we just use tile_xy_pair instead?
		int x, y;
		xy_ix_t() : x(0), y(0) {}
		xy_ix_t(int x_, int y_) : x(x_), y(y_) {}
		bool operator==(xy_ix_t const &v) const {return (x == v.x && y == v.y);}
		bool operator!=(xy_ix_t const &v) const {return !operator==(v);}
		bool operator< (xy_ix_t const &v) const {return (y == v.y) ? (x < v.x) : (y < v.y);} // compare y then x
	};
	struct cmp_xy_by_dist {
		xy_ix_t const &dest;
		cmp_xy_by_dist(xy_ix_t const &dest_) : dest(dest_) {}
		
		bool operator()(xy_ix_t const &a, xy_ix_t const &b) const { // comapre Euclidean distan squared
			int const dax(a.x - dest.x), day(a.y - dest.y), dbx(b.x - dest.x), dby(b.y - dest.y);
			return ((dax*dax + day*day) < (dbx*dbx + dby*dby));
		}
	};
	vect_cube_t blockers;
	heightmap_query_t const &hq;
public:
	road_router_t(vect_cube_t const &blockers_, heightmap_query_t const &hq_, cube_t const exclude[2]) : blockers(blockers_), hq(hq_) {
		// add the two exclude (city) cubes back into blockers; is it better to just pass blockers + active_blockers?
		for (unsigned d = 0; d < 2; ++d) {
			//blockers.push_back(exclude[d]);
			//blockers.back().expand_by(-HALF_DXY); // shrink to prevent initial collisions
		}
	}
private:
	xy_ix_t get_ix_for_pos(point const &pos) const {return xy_ix_t(hq.get_x_pos(pos.x), hq.get_y_pos(pos.y));}

	bool is_cell_valid(int x, int y, float zmin=0.0, float zmax=0.0) const {
		if (!hq.is_inside_terrain(x, y)) return 0;
		cube_t const c(hq.get_cube_for_cell(x, y));
		if (c.z1() < water_plane_z) return 0; // underwater
		if (zmin != zmax && (c.z1() < zmin || c.z2() > zmax)) return 0; // delta-z is too large
		return !has_bcube_int_xy_no_adj(c, blockers);
	}
	bool is_adj_cell_valid(int x, int y, float zval) const {
		float const max_dz(HALF_DXY*city_params.max_road_slope); // dz/dx < city_params.max_road_slope
		return is_cell_valid(x, y, zval-max_dz, zval+max_dz);
	}
	float get_valid_adj_cells(xy_ix_t const &cell, vector<xy_ix_t> &adj) const {
		assert(hq.is_inside_terrain(cell.x, cell.y));
		float const zval(hq.get_height(cell.x, cell.y));

		for (unsigned d = 0; d < 4; ++d) { // try +x, -x, +y, and -y
			xy_ix_t adj_cell(cell);
			((d>>1) ? adj_cell.y : adj_cell.x) += ((d&1) ? 1 : -1);
			if (is_adj_cell_valid(adj_cell.x, adj_cell.y, zval)) {adj.push_back(adj_cell);}
		}
		return zval;
	}
public:
	bool connect_endpoints(point const &p1, point const &p2, vector<point> &pts, bool dim1, bool dim2) { // do we need to pass in dir?
		// TODO: use A* path finding on the terrain heightmap, with heavy weight for introducing jogs
		xy_ix_t const start(get_ix_for_pos(p1)), end(get_ix_for_pos(p2));
		if (start == end) return 1; // success
		//deque<xy_ix_t> pend;
		vector<xy_ix_t> pend, next;
		pend.push_back(start);
		map<xy_ix_t, xy_ix_t> seen; // map to parent index; maybe could be a bit vector, but the heightmap could be very large (50M entries)
		seen.insert(make_pair(start, start));
		unsigned num_visited(0), path_len(0);
		cout << "connecting endpoints " << TXT(start.x) << TXT(start.y) << TXT(end.x) << TXT(end.y) << endl;

		while (!pend.empty()) {
			++num_visited;
			xy_ix_t cur(pend.back());
			pend.pop_back();
			auto it(seen.find(cur));
			assert(it != seen.end());
			xy_ix_t const prev(it->second);
			next.clear();
			float const zval(get_valid_adj_cells(cur, next));
			sort(next.begin(), next.end(), cmp_xy_by_dist(end)); // sort by distance to end

			for (auto i = next.begin(); i != next.end(); ++i) {
				if (*i == end) { // reached the end
					unsigned const pts_start(pts.size()); // typically 1
					bool cur_dim(dim2); // current road dim, starting from the end

					while (cur != start) { // reconstruct the path by walking backwards using seen
						auto it2(seen.find(cur));
						assert(it2 != seen.end());
						xy_ix_t const &parent(it2->second);
						assert(parent != cur); // check for cycles
						assert((cur.x == parent.x) != (cur.y == parent.y)); // must move in either x or y

						if ((cur.y == parent.y) != cur_dim) { // direction change, add a bend
							pts.push_back(hq.get_pos_for_cell(cur.x, cur.y));
							cur_dim ^= 1;
						}
						cur = parent;
						++path_len;
					} // for while()
					std::reverse(pts.begin()+pts_start, pts.end());
					unsigned const min_len(abs(end.x - start.x) + abs(end.y - start.y));
					cout << TXT(seen.size()) << TXT(pts.size()) << TXT(min_len) << TXT(path_len) << TXT(num_visited) << endl;
					return 1; // done
				}
				if (!seen.insert(make_pair(*i, cur)).second) continue; // already seen
				bool const is_turn(i->x != prev.x && i->y != prev.y);
				float const dz(fabs(zval - hq.get_height(i->x, i->y)));
				float const cost(dz + 1.0*is_turn*city_params.road_width);
				pend.push_back(*i);
			} // for i
		} // end while()
		cout << "failed" << endl;
		return 0;
	}
}; // end road_router_t

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
// Note: this is a simple fast cost estimate based on road length and elevation change at endpoints; it doesn't include elevation change along the length of the road
float get_road_cost(point const &p1, point const &p2) {return (p2p_dist(p1, p2) + 4.0f*fabs(p1.z - p2.z));} // 4x penalty for steepness

bool road_seg_valid(point const &p1, point const &p2, bool dim, vect_cube_t const &blockers, float road_hwidth, bool expand_start, bool expand_end) {
	float const length(fabs(p1[dim] - p2[dim]));
	if (length < 4.0*road_hwidth) return 0; // too short
	if (fabs(p1.z - p2.z)/(length - road_hwidth) > city_params.max_road_slope) return 0; // check slope
	return !has_bcube_int_xy_no_adj(get_road_between_pts(p1, p2, road_hwidth, expand_start, expand_end), blockers);
}
bool check_pt_valid(point const &pt, cube_t const exclude[2]) {
	return (pt.z > water_plane_z && !exclude[0].contains_pt_xy(pt) && !exclude[1].contains_pt_xy(pt));
}

float city_road_connector_t::find_route_between_points(point const &p1, point const &p2, vect_cube_t const &blockers, heightmap_query_t const &hq,
	vector<point> &pts, cube_t const exclude[2], float road_hwidth, bool dim1, bool dir1, bool dim2, bool dir2)
{
	pts.push_back(p1);
	road_router_t router(blockers, hq, exclude);
	bool success(0);
	
	if (router.connect_endpoints(p1, p2, pts, dim1, dim2)) {
		success = 1;
	}
	else if (dim1 == dim2) { // add 2 points to create a job
		assert(dir1 != dir2); // must be opposing
		point pt[2];
		pt[0][ dim1] = pt[1][ dim1] = 0.5f*(p1[dim1] + p2[dim1]); // halfway between the two end points
		pt[0][!dim1] = p1[!dim1]; pt[1][!dim1] = p2[!dim1];
		for (unsigned d = 0; d < 2; ++d) {pt[d].z = hq.get_road_zval_at_pt(pt[d]);}
		
		if (check_pt_valid(pt[0], exclude) && check_pt_valid(pt[1], exclude) && road_seg_valid(p1, pt[0], dim1, blockers, road_hwidth, 0, 1) &&
			road_seg_valid(pt[0], pt[1], !dim1, blockers, road_hwidth, 1, 1) && road_seg_valid(pt[1], p2, dim1, blockers, road_hwidth, 1, 0))
		{
			for (unsigned d = 0; d < 2; ++d) {pts.push_back(pt[d]);}
			success = 1;
		}
	}
	else { // add one point to create a right angle bend
		point ipt;
		ipt[ dim1] = p2[ dim1];
		ipt[!dim1] = p1[!dim1];
		ipt.z = hq.get_road_zval_at_pt(ipt);

		if (check_pt_valid(ipt, exclude) && road_seg_valid(p1, ipt, dim1, blockers, road_hwidth, 0, 1) && road_seg_valid(ipt, p2, dim2, blockers, road_hwidth, 1, 0)) {
			pts.push_back(ipt);
			success = 1;
		}
	}
	if (!success) {pts.clear(); return 0.0;} // failed
	pts.push_back(p2);
	// calculate cost
	float cost(0.0);
	for (auto p = pts.begin(); p+1 != pts.end(); ++p) {cost += get_road_cost(*p, *(p+1));}
	cost *= (pts.size() - 1); // add jog penalty
	return cost;
}

bool city_road_connector_t::segment_road(road_t const &road, heightmap_query_t const &hq, bool check_only) {
	bool const dim(road.dim);
	float const road_len(road.get_length()), conn_pos(road.get_center_dim(!dim));
	unsigned const num_segs(ceil(road_len/city_params.conn_road_seg_len));
	assert(num_segs > 0 && num_segs < 1000); // sanity check
	float const seg_len(road_len/num_segs);
	assert(seg_len <= city_params.conn_road_seg_len);
	road_t rs(road); // keep d[!dim][0], d[!dim][1], dim, and road_ix
	rs.z1() = road.d[2][road.slope];
	segments.clear();

	for (unsigned n = 0; n < num_segs; ++n) {
		rs.d[dim][1] = ((n+1 == num_segs) ? road.d[dim][1] : (rs.d[dim][0] + seg_len)); // make sure it ends exactly at the correct location
		point pos;
		pos[ dim] = rs.d[dim][1];
		pos[!dim] = conn_pos;
		rs.z2()   = hq.get_road_zval_at_pt(pos); // terrain height at end of segment
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

