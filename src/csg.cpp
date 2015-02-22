// 3D World - Shape primitive drawing
// by Frank Gennari
// 5/1/05
#include "3DWorld.h"
#include "mesh.h"
#include "csg.h"
#include "cobj_bsp_tree.h"


bool const VOXEL_MERGE        = 0;
bool const SUB_CUBE_MERGE     = 1;
bool const UNOVERLAP_COBJS    = 1;
bool const MERGE_COBJS        = 1;
bool const CHECK_COBJS        = 1;
bool const EFLAGS_STRICT      = 0;
bool const ONLY_SUB_PREV_NEG  = 1;
bool const CHECK_ADJACENCY    = 0; // doesn't seem to make any significant difference
int  const REMOVE_T_JUNCTIONS = 1; // fewer hole pixels and better ambient transitions but more cobjs and render time
float const REL_DMAX          = 0.2;


extern bool use_waypoints;


// *** RECT IMPLEMENTATION ***


rect::rect(float const r[3][2], unsigned d0, unsigned d1) { // projection from 3D => 2D

	assert(d0 >= 0 && d0 <= 3 && d1 >= 0 && d1 <= 3 && d0 != d1);
	d[0][0] = r[d0][0]; d[0][1] = r[d0][1]; d[1][0] = r[d1][0]; d[1][1] = r[d1][1];
	//assert(nonzero());
}


void rect::clip_to(float const c[2][2]) {

	for (unsigned i = 0; i < 2; ++i) {
		assert(c[i][0] < c[i][1]);
		d[i][0] = max(c[i][0], min(c[i][1], d[i][0]));
		d[i][1] = max(c[i][0], min(c[i][1], d[i][1]));
	}
}


// rr will be removed
void rect::subtract_from(rect const &rr, deque<rect> &new_rects) const { // subtract ourself from rr

	if (contains(rr.d)) return;
	unsigned i2[2], n[2];
	char vox[3][3]; // 3x3 pixels
	float vals[2][4];
	rect r2;

	// determine cutting planes
	for (unsigned i = 0; i < 2; ++i) {
		n[i]       = 0;
		vals[i][0] = rr.d[i][0];

		for (unsigned j = 0; j < 2; ++j) {
			if (d[i][j] > rr.d[i][0] && d[i][j] < rr.d[i][1]) vals[i][++n[i]] = d[i][j];
		}
		vals[i][++n[i]] = rr.d[i][1];
	}

	// build pixel table
	for (i2[0] = 0; i2[0] < n[0]; ++i2[0]) {
		for (i2[1] = 0; i2[1] < n[1]; ++i2[1]) {
			float pt[2];
			for (unsigned l = 0; l < 2; ++l) { // pt is the center of this pixel
				pt[l] = 0.5*(vals[l][i2[l]] + vals[l][i2[l]+1]);
			}
			vox[i2[0]][i2[1]] = !contains_pt(pt);
		}
	}

	// merge adjacent pixels (optional performance improvement)
	for (unsigned dim = 0; dim < 2; ++dim) { // rotate the cube slicing direction
		unsigned const d1(!dim);

		for (i2[dim] = 0; i2[dim] < n[dim]; ++i2[dim]) {
			bool all_in(1);

			for (i2[d1] = 0; i2[d1] < n[d1] && all_in; ++i2[d1]) {
				if (!vox[i2[0]][i2[1]]) all_in = 0;
			}
			if (all_in) {
				for (unsigned l = 0; l < 2; ++l) {
					for (unsigned m = 0; m < 2; ++m) {
						r2.d[l][m] = ((l == dim) ? vals[l][i2[l]+m] : rr.d[l][m]); // slice dimensions
					}
				}
				if (!r2.is_near_zero_area()) new_rects.push_back(r2);

				for (i2[d1] = 0; i2[d1] < n[d1]; ++i2[d1]) {
					vox[i2[0]][i2[1]] = 0; // remove pixel
				}
			}
		}
	}

	// generate output rects
	for (i2[0] = 0; i2[0] < n[0]; ++i2[0]) {
		for (i2[1] = 0; i2[1] < n[1]; ++i2[1]) {
			if (vox[i2[0]][i2[1]]) {
				for (unsigned l = 0; l < 2; ++l) {
					for (unsigned m = 0; m < 2; ++m) {
						r2.d[l][m] = vals[l][i2[l]+m];
					}
				}
				if (!r2.is_near_zero_area()) new_rects.push_back(r2);
			}
		}
	}
}


bool rect::merge_with(rect const &r) {

	for (unsigned j = 0; j < 2; ++j) {
		if (r.d[j][0] == d[j][0] && r.d[j][1] == d[j][1]) {
			for (unsigned k = 0; k < 2; ++k) {
				if (r.d[!j][!k] == d[!j][k]) {
					d[!j][k] = r.d[!j][k];
					return 1;
				}
			}
		}
	}
	return 0;
}


void rect::print() const {

	for (unsigned i = 0; i < 2; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			cout << d[i][j];
			if (j == 0) cout << " ";
		}
		if (i == 0) cout << ", ";
	}
}


// *** CUBE_T IMPLEMENTATION ***


void cube_t::set_from_points(point const *const pts, unsigned npts) {

	assert(npts > 0);
	UNROLL_3X(d[i_][0] = d[i_][1] = pts[0][i_];)
	for (unsigned i = 1; i < npts; ++i) {union_with_pt(pts[i]);} // get bounding xy rectangle
}


void cube_t::print() const {

	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			cout << d[i][j];
			if (i < 2 || j < 1) cout << ",";
		}
		cout << " ";
	}
}


bool cube_t::is_near_zero_area() const {
		
	UNROLL_3X(if (fabs(d[i_][0] - d[i_][1]) < TOLER) return 1;)
	return 0;
}


unsigned cube_t::get_split_dim(float &max_sz, float &sval, unsigned skip_dims) const {

	unsigned dim(0);
	max_sz = 0;

	for (unsigned i = 0; i < 3; ++i) {
		if (skip_dims & (1 << i)) continue;
		float const dim_sz(d[i][1] - d[i][0]);
		assert(dim_sz >= 0.0);
		
		if (dim_sz > max_sz) {
			max_sz = dim_sz;
			dim    = i;
		}
	}
	if (max_sz > 0.0) {
		sval = get_cube_center()[dim]; // center point (mean seems to work better than median)
		assert(!(skip_dims & (1 << dim)));
	}
	return dim;
}


bool cube_t::cube_intersection(const cube_t &cube, cube_t &res) const { // flags are not set
	
	for (unsigned i = 0; i < 3; ++i) {
		res.d[i][0] = max(d[i][0], cube.d[i][0]);
		res.d[i][1] = min(d[i][1], cube.d[i][1]);
		if (res.d[0] >= res.d[1]) return 0; // no intersection
	}
	return 1;
}


float cube_t::get_overlap_volume(const cube_t &cube) const {

	float volume(1.0);

	for (unsigned i = 0; i < 3; ++i) {
		float const val(min(d[i][1], cube.d[i][1]) - max(d[i][0], cube.d[i][0]));
		if (val <= 0.0) return 0.0; // no intersection
		volume *= val;
	}
	return volume;
}


vector3d cube_t::closest_side_dir(point const &pos) const { // for fragment velocity

	int dir(-1);
	float mdist(0.0);
	vector3d dv(zero_vector);

	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			float const dist(fabs(d[i][j] - pos[i]));

			if (dir == -1 || dist < mdist) {
				mdist = dist;
				dir   = (i << 1) + j;
			}
		}
	}
	dv[dir >> 1] = ((dir & 1) ? 1.0 : -1.0);
	return dv;
}


point cube_t::gen_rand_pt_in_cube() const {

	point pt;

	for (unsigned j = 0; j < 3; ++j) {
		pt[j] = rand_uniform(d[j][0], d[j][1]);
	}
	return pt;
}


int cube_t::closest_face(point const &pos) const {

	int face(0);
	float min_dist(0.0);

	for (unsigned dim = 0; dim < 3; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			float const dist(fabs(pos[dim] - d[dim][dir]));
				
			if (min_dist == 0 || dist < min_dist) {
				min_dist = dist;
				face     = (dim<<1) + dir;
			}
		}
	}
	return face;
}


bool cube_t::cube_merge(cube_t const &cube) { // simplified version of csg_cube::cube_merge() without the edge_flags

	unsigned compat[3], nc(0), ci(0);

	for (unsigned i = 0; i < 3; ++i) { // check compatability
		compat[i] = 1;
		for (unsigned j = 0; j < 2; ++j) {
			if (cube.d[i][j] != d[i][j]) compat[i] = 0;
		}
		if (compat[i]) ++nc; else ci = i;
	}
	if (nc >= 2) { // compatible
		if (nc == 3) return 1; // same cube, remove it

		for (unsigned i = 0; i < 2; ++i) { // only one iteration will merge something
			if (cube.d[ci][1-i] == d[ci][i]) { // adjacent
				d[ci][i] = cube.d[ci][i];
				return 1;
			}
		}
	}
	return 0;
}


void cube_t::get_points(point pts[8]) const {

	unsigned i[3];

	for (i[0] = 0; i[0] < 2; ++i[0]) {
		for (i[1] = 0; i[1] < 2; ++i[1]) {
			for (i[2] = 0; i[2] < 2; ++i[2]) {
				UNROLL_3X(pts[(((i[0]<<1)+i[1])<<1)+i[2]][i_] = d[i_][i[i_]];)
			}
		}
	}
}


// *** CSG_CUBE IMPLEMENTATION ***


// returns 1 if entire cube is removed
bool csg_cube::subtract_from_internal(const csg_cube &cube, vector<csg_cube> &output) const { // subtract ourself from cube

	if (contains_cube(cube)) return 1;
	unsigned const u[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}}; // unit vectors
	unsigned i3[3], n[3];
	char vox[3][3][3]; // 3x3x3 voxels
	float vals[3][4];

	// determine cutting planes
	for (unsigned i = 0; i < 3; ++i) {
		n[i]       = 0;
		vals[i][0] = cube.d[i][0];

		for (unsigned j = 0; j < 2; ++j) {
			if (d[i][j] > (cube.d[i][0] + TOLER) && d[i][j] < (cube.d[i][1] - TOLER)) vals[i][++n[i]] = d[i][j];
		}
		vals[i][++n[i]] = cube.d[i][1];
	}

	// build voxel table
	for (i3[0] = 0; i3[0] < n[0]; ++i3[0]) {
		for (i3[1] = 0; i3[1] < n[1]; ++i3[1]) {
			for (i3[2] = 0; i3[2] < n[2]; ++i3[2]) {
				point pt;
				for (unsigned l = 0; l < 3; ++l) { // pt is the center of this voxel
					pt[l] = 0.5*(vals[l][i3[l]] + vals[l][i3[l]+1]);
				}
				vox[i3[0]][i3[1]][i3[2]] = !contains_pt(pt);
			}
		}
	}

	// merge adjacent voxels
	if (VOXEL_MERGE) {
		for (unsigned dim = 0; dim < 3; ++dim) { // rotate the cube slicing direction
			unsigned const d1((dim+1)%3), d2((dim+2)%3);

			for (i3[dim] = 0; i3[dim] < n[dim]; ++i3[dim]) {
				bool all_in(1);

				for (i3[d1] = 0; i3[d1] < n[d1] && all_in; ++i3[d1]) {
					for (i3[d2] = 0; i3[d2] < n[d2] && all_in; ++i3[d2]) {
						if (!vox[i3[0]][i3[1]][i3[2]]) all_in = 0;
					}
				}
				if (all_in) { // full 2D slice
					unsigned char edgeflags(cube.eflags); // |= vs. &=~ ???
					if (i3[dim] != 0)        edgeflags &= ~EFLAGS[dim][0]; // first  side not on original cube edge
					if (i3[dim] != n[dim]-1) edgeflags &= ~EFLAGS[dim][1]; // second side not on original cube edge
					csg_cube ncube(edgeflags);

					for (unsigned l = 0; l < 3; ++l) {
						for (unsigned m = 0; m < 2; ++m) {
							ncube.d[l][m] = ((l == dim) ? vals[l][i3[l]+m] : cube.d[l][m]); // slice dimensions
						}
					}
					output.push_back(ncube);

					for (i3[d1] = 0; i3[d1] < n[d1]; ++i3[d1]) {
						for (i3[d2] = 0; i3[d2] < n[d2]; ++i3[d2]) {
							vox[i3[0]][i3[1]][i3[2]] = 0; // remove voxel
						}
					}
				}
			}
		}
	}

	// generate output cubes
	for (i3[0] = 0; i3[0] < n[0]; ++i3[0]) {
		for (i3[1] = 0; i3[1] < n[1]; ++i3[1]) {
			for (i3[2] = 0; i3[2] < n[2]; ++i3[2]) {
				if (!vox[i3[0]][i3[1]][i3[2]]) continue;
				unsigned char edgeflags(0);

				for (unsigned l = 0; l < 3; ++l) { // set edge flags to remove inside faces
					assert(n[l] > 0 && n[l] <= 3);
					if (i3[l] > 0) {
						if (vox[i3[0]-u[0][l]][i3[1]-u[1][l]][i3[2]-u[2][l]]) edgeflags |= EFLAGS[l][0];
					}
					else edgeflags |= (cube.eflags & EFLAGS[l][0]);
					
					if (i3[l] < n[l]-1) {
						if (vox[i3[0]+u[0][l]][i3[1]+u[1][l]][i3[2]+u[2][l]]) edgeflags |= EFLAGS[l][1];
					}
					else edgeflags |= (cube.eflags & EFLAGS[l][1]);
				}
				csg_cube ncube(edgeflags);
				
				for (unsigned l = 0; l < 3; ++l) {
					for (unsigned m = 0; m < 2; ++m) {
						ncube.d[l][m] = vals[l][i3[l]+m];
					}
				}
				if (SUB_CUBE_MERGE) {
					bool merged(0);
					for (unsigned c = 0; c < output.size() && !merged; ++c) {
						if (output[c].cube_merge(ncube)) merged = 1;
					}
					if (merged) continue;
				}
				if (!ncube.is_near_zero_area()) output.push_back(ncube);
			}
		}
	}
	return 0;
}


csg_cube::csg_cube(const coll_obj &cobj, bool use_bounding_cube) : eflags(cobj.cp.surfs) { // coll_obj constructor

	assert(use_bounding_cube || cobj.type == COLL_CUBE);
	copy_from(cobj);
	normalize();
}


inline void csg_cube::write_to_cobj(coll_obj &cobj) const {

	assert(cobj.type == COLL_CUBE);
	cobj.copy_from(*this);
	cobj.cp.surfs = eflags;
}


bool csg_cube::cube_intersection(const csg_cube &cube, csg_cube &res) const { // flags are not set

	res.eflags = 0; // fix later
	return cube_t::cube_intersection(cube, res);
}


// returns 1 if some work is done
bool csg_cube::subtract_from_cube(coll_obj_group &new_cobjs, coll_obj const &cobj) const { // subtract ourself from cobjs[index]

	assert(cobj.type == COLL_CUBE);
	if (!quick_intersect_test(cobj)) return 0; // no intersection
	csg_cube cube(cobj);
	if (!intersects(cube, TOLER))    return 0; // no intersection
	if (cube.is_zero_area())         return 1; // if zero area then remove it entirely
	vector<csg_cube> output;
	subtract_from_internal(cube, output);
	
	for (unsigned i = 0; i < output.size(); ++i) {
		new_cobjs.push_back(cobj); // keep properties of old coll cube
		output[i].write_to_cobj(new_cobjs.back());
		//assert(!intersects(new_cobjs.back()));
	}
	return 1;
}


float get_closest_pt_on_line_t(point const &pos, point const &l1, point const &l2) { // use for line lights as well?

	vector3d const L(l2 - l1);
	return CLIP_TO_01(dot_product((pos - l1), L)/L.mag_sq());
}


// subtract ourself from cobj; returns 1 if some work is done
bool csg_cube::subtract_from_cylinder(coll_obj_group &new_cobjs, coll_obj &cobj) const { // subtract ourself from cobjs[index]

	assert(cobj.is_cylinder());
	float const radius(max(cobj.radius, cobj.radius2)); // containment/intersection tests are conservative
	point const &p0(cobj.points[0]), &p1(cobj.points[1]);
	int const p(get_max_dim(p0 - p1));
	unsigned const m((p+1)%3), n((p+2)%3);
	
	if (p0[m] != p1[m] || p0[n] != p1[n]) { // not axis aligned, we approximate the cube as sphere bounding the min dim
		if (cobj.radius != cobj.radius2) return 0; // not handling non-axis aligned truncated cone
		float const radius(0.5*min_len());
		if (radius < cobj.radius) return 0; // cylinder is too thick
		point const center(get_cube_center());
		if (!pt_line_dist_less_than(center, p0, p1, radius)) return 0; // centerline is too far
		float const delta_t(radius/p2p_dist(p0, p1));
		float const t(get_closest_pt_on_line_t(center, p0, p1)), t0(t - delta_t), t1(t + delta_t);
		
		if (t0 > 0.0) { // add p0 segment
			new_cobjs.push_back(cobj); // keep properties of old coll cylinder
			new_cobjs.back().cp.surfs  = 0; // reset edge flags in case the ends become exposed
			new_cobjs.back().points[1] = p0 + t0*(p1 - p0); // set new p1
		}
		if (t1 < 1.0) { // add p1 segment
			new_cobjs.push_back(cobj); // keep properties of old coll cylinder
			new_cobjs.back().cp.surfs  = 0; // reset edge flags in case the ends become exposed
			new_cobjs.back().points[0] = p0 + t1*(p1 - p0); // set new p0
		}
		return 1;
	}
	assert(p0[p] != p1[p]);
		
	if (p0[p] > p1[p]) {
		swap(cobj.points[0], cobj.points[1]); // upside down (note that p0/p1 are also swapped)
		swap(cobj.radius,    cobj.radius2);
	}
	float const c[2][2] = {{p0[m]-radius, p0[m]+radius}, {p0[n]-radius, p0[n]+radius}};
	if (c[0][0] <  d[m][0] || c[0][1] >  d[m][1]) return 0; // no m-containment
	if (c[1][0] <  d[n][0] || c[1][1] >  d[n][1]) return 0; // no n-containment
	if (p0[p]   >= d[p][1] || p1[p]   <= d[p][0]) return 0; // no p-intersection
	if (p0[p]   >= d[p][0] && p1[p]   <= d[p][1]) return 1; // p-containment - remove
	csg_cube const cube(cobj, 1);
	if (cube.is_zero_area()) return 1; // if zero area then remove it entirely
	cobj.points[0][m] = cobj.points[1][m]; // update cobj.d?
	cobj.points[0][n] = cobj.points[1][n];
	unsigned nv(0);
	float vals[4];

	if (p0[p] < d[p][0]) {
		vals[nv++] = p0[p];   // A
		vals[nv++] = d[p][0]; // B
	}
	if (p1[p] > d[p][1]) {
		vals[nv++] = d[p][1]; // C
		vals[nv++] = p1[p];   // D
	}
	for (unsigned i = 0; i < nv; i += 2) {
		new_cobjs.push_back(cobj); // keep properties of old coll cylinder
		coll_obj &ncobj(new_cobjs.back());
		ncobj.cp.surfs = 0; // reset edge flags in case the ends become exposed
			
		for (unsigned j = 0; j < 2; ++j) {
			ncobj.points[j][p] = vals[i+j];
		}
		if (cobj.radius != cobj.radius2) { // calculate new radius values
			float rv[2];

			for (unsigned j = 0; j < 2; ++j) {
				rv[j] = cobj.radius + (cobj.radius2 - cobj.radius)*(vals[i+j] - p0[p])/(p1[p] - p0[p]);
			}
			ncobj.radius  = rv[0];
			ncobj.radius2 = rv[1];
		}
	} // for i
	return 1;
}


// returns 1 if some work is done
// see http://www.cs.fit.edu/~wds/classes/graphics/Clip/clip/clip.html
bool csg_cube::subtract_from_polygon(coll_obj_group &new_cobjs, coll_obj const &cobj) const { // subtract ourself from cobjs[index]

	// start by assuming *this intersects cobj.d (should have been tested already)
	assert(cobj.is_thin_poly()); // can't handle other cases yet
	if (contains_cube(cobj)) return 1; // contained - remove the entire cobj
	static vector<point> cur, next, new_poly;
	assert(cur.empty() && next.empty() && new_poly.empty());
	for (int i = 0; i < cobj.npoints; ++i) {cur.push_back(cobj.points[i]);}
	size_t const init_sz(new_cobjs.size());

	for (unsigned i = 0; i < 3 && !cur.empty(); ++i) {
		for (unsigned j = 0; j < 2 && !cur.empty(); ++j) {
			float const clip_val(d[i][j]); // clip cur polygon by this plane
			bool prev_outside(0);
			// put the outside part (tri/quad) (if any) in new_cobjs
			// put the inside part  (tri/quad) (if any) in cur
			for (unsigned p = 0; p <= cur.size(); ++p) {
				point const &pos(cur[(p == cur.size()) ? 0 : p]);
				bool const cur_outside(((pos[i] < clip_val) ^ j) != 0);
				bool write_int(0), write_cur(0);
				
				if (p == cur.size()) { // last point
					if (cur_outside != prev_outside) {write_int = 1;} // edge crossing
				}
				else if (p == 0 || prev_outside == cur_outside) { // first point or no edge crossing
					write_cur = 1;
				}
				else { // interior point, edge crossing
					write_int = write_cur = 1;
				}
				if (write_int) {
					vector3d const edge(pos - cur[p-1]);
					float const t((clip_val - cur[p-1][i])/edge[i]);
					point const p_int(cur[p-1] + edge*t);
					new_poly.push_back(p_int);
					next.push_back(p_int);
				}
				if (write_cur) {(cur_outside ? new_poly : next).push_back(pos);}
				prev_outside = cur_outside;
			} // for p
			if (!new_poly.empty()) {
				bool const split_quads(use_waypoints); // FIXME: waypoint issues with split polygons
				split_polygon_to_cobjs(cobj, new_cobjs, new_poly, split_quads);
				new_poly.clear();
			}
			cur.swap(next);
			next.clear();
		} // for j
	} // for i
	if (!cur.empty()) { // the remainder (cur) is the part to be removed
		cur.clear();
		return 1;
	}
	// else nothing removed
	assert(new_cobjs.size() >= init_sz); // can be equal if all pieces are tiny fragments that get removed
	new_cobjs.erase(new_cobjs.begin()+init_sz, new_cobjs.end()); // remove everything that was added
	return 0;
}

// same as above code, but ignores the outside part of the polygon
void clip_polygon_to_cube(cube_t const &cube, point const *const pts_in, unsigned npts_in, cube_t const &pts_bcube, vector<point> &pts_out) {

	if (!cube.intersects(pts_bcube)) {pts_out.clear(); return;}
	static vector<point> next;
	pts_out.resize(npts_in);
	for (unsigned i = 0; i < npts_in; ++i) {pts_out[i] = pts_in[i];}
	if (cube.contains_cube(pts_bcube)) return;

	for (unsigned i = 0; i < 3 && !pts_out.empty(); ++i) {
		for (unsigned j = 0; j < 2 && !pts_out.empty(); ++j) {
			float const clip_val(cube.d[i][j]); // clip cur polygon by this plane
			if (((pts_bcube.d[i][j] < clip_val) ^ j) == 0) continue; // don't need to clip in this dim
			bool prev_outside(0); // put the inside part (if any) in pts_out

			for (unsigned p = 0; p <= pts_out.size(); ++p) {
				point const &pos(pts_out[(p == pts_out.size()) ? 0 : p]);
				bool const cur_outside(((pos[i] < clip_val) ^ j) != 0);
				bool write_int(0), write_cur(0);
				
				if (p == pts_out.size()) { // last point
					if (cur_outside != prev_outside) {write_int = 1;} // edge crossing
				}
				else if (p == 0 || prev_outside == cur_outside) { // first point or no edge crossing
					write_cur = 1;
				}
				else { // interior point, edge crossing
					write_int = write_cur = 1;
				}
				if (write_int) {
					vector3d const edge(pos - pts_out[p-1]);
					float const t((clip_val - pts_out[p-1][i])/edge[i]);
					next.push_back(pts_out[p-1] + edge*t);
				}
				if (write_cur && !cur_outside) {next.push_back(pos);}
				prev_outside = cur_outside;
			} // for p
			pts_out.swap(next);
			next.clear();
		} // for j
	} // for i
}


bool csg_cube::subtract_from_thick_polygon(coll_obj_group &new_cobjs, coll_obj const &cobj) const { // subtract ourself from cobjs[index]

	assert(cobj.type == COLL_POLYGON && cobj.thickness > 0.0);
	assert(cobj.npoints == 3 || cobj.npoints == 4); // quad or triangle
	if (contains_cube(cobj)) return 1; // contained - remove the entire cobj

	if (is_axis_aligned(cobj.norm)) { // surprisingly not a common case, considering many axis aligned shapes can be represented with cubes
		// FIXME: special case where we can be more efficient or do a better job?
	}
	coll_obj polys[6]; // top, bot, and up to 4 sides
	vector<tquad_t> pts;
	thick_poly_to_sides(cobj.points, cobj.npoints, cobj.norm, cobj.thickness, pts); // split into 5-6 polygons, one per side
	assert(pts.size() <= 6);
	unsigned cur(0);

	for (unsigned i = 0; i < pts.size(); ++i) {
		polys[cur] = cobj;
		polys[cur].thickness = 0.1*MIN_POLY_THICK; // thin poly: 0 > val >= MIN_POLY_THICK
		polys[cur].norm      = pts[i].get_norm();
		polys[cur].set_from_pts(pts[i].pts, pts[i].npts);
		bool const removed(subtract_from_polygon(new_cobjs, polys[cur]));
		if (!removed) {++cur;} // if not removed, keep it
	}
	if (cur == pts.size()) return 0; // nothing was removed, no changes to make

	for (unsigned i = 0; i < cur; ++i) {
		new_cobjs.push_back(polys[i]);
	}
	// FIXME: add polygons to fill the gaps between the cube and the remaining polygons
	return 1;
}


float get_cube_dmax() {

	return REL_DMAX*(X_SCENE_SIZE + Y_SCENE_SIZE);
}


// could do this dynamically as cubes are split
bool csg_cube::cube_merge(csg_cube &cube) {

	unsigned compat[3], nc(0), ci(0);

	for (unsigned i = 0; i < 3; ++i) { // check compatability
		compat[i] = 1;
		for (unsigned j = 0; j < 2; ++j) {
			if (cube.d[i][j] != d[i][j]) compat[i] = 0;
			else if (EFLAGS_STRICT && ((cube.eflags ^ eflags) & EFLAGS[i][j])) compat[i] = 0;
		}
		if (compat[i]) ++nc; else ci = i;
	}
	if (nc >= 2) { // compatible
		if (nc == 3) { // same cube, remove it
			for (unsigned i = 0; i < 3; ++i) { // merge edge flags
				for (unsigned j = 0; j < 2; ++j) {
					if (!(cube.eflags & EFLAGS[i][j])) eflags &= ~EFLAGS[i][j];
				}
			}
			return 1;
		}
		float const dmax(get_cube_dmax());

		for (unsigned i = 0; i < 2; ++i) { // only one iteration will merge something
			if (cube.d[ci][1-i] == d[ci][i]) { // adjacent
				// only size test for real eflags cubes
				float const dval(fabs(cube.d[ci][i] - d[ci][1-i]));
				if (dval > dmax) return 0; // resulting cube will be too large
				d[ci][i] = cube.d[ci][i];
				eflags  &= ~EFLAGS[ci][i];
				eflags  |= (cube.eflags & EFLAGS[ci][i]);

				for (unsigned j = 0; j < 3; ++j) { // set shared edge flags
					if (j != ci) {
						for (unsigned k = 0; k < 2; ++k) {
							if (!(cube.eflags & EFLAGS[j][k])) eflags &= ~EFLAGS[j][k];
						}
					}
				}
				return 1;
			}
		}
	}
	if (CHECK_ADJACENCY /*&& nc == 1*/) { // no merge, determine contained faces (does gridbag need overlap to catch all of these cases?)
		bool adjacent(0);

		for (unsigned i = 0; i < 3 && !adjacent; ++i) { // could we use ci or something?
			for (unsigned j = 0; j < 2 && !adjacent; ++j) {
				if (cube.d[i][j] == d[i][!j]) { // potential adjacency
					adjacent = 1;
					bool contained[2] = {1, 1}; // cube,self = {c,s}
					for (unsigned k = 0; k < 3; ++k) {
						if (k == i) continue;
						for (unsigned l = 0; l < 2; ++l) {
							if (cube.d[k][l] < d[k][l]) contained[ l] = 0;
							if (cube.d[k][l] > d[k][l]) contained[!l] = 0;
						}
					}
					assert(!(contained[0] && contained[1])); // nc should be at least 2 if we get here
					
					if (contained[0]) {
						cube.eflags |= EFLAGS[i][j];
					}
					else if (contained[1]) {
						eflags |= EFLAGS[i][!j];
					}
				}
			}
		}
	}
	return 0;
}


void csg_cube::unset_adjacent_edge_flags(coll_obj &cobj) const {

	assert(cobj.type == COLL_CUBE);

	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			if (fabs(d[i][j] - cobj.d[i][!j]) < TOLER) { // shared opposing edge
				bool overlaps(1);

				for (unsigned k = 0; k < 3 && !overlaps; ++k) {
					if (k != i && (d[k][0] >= cobj.d[k][1] || d[k][1] <= cobj.d[k][0])) overlaps = 0;
				}
				if (overlaps) {
					cobj.cp.surfs &= ~EFLAGS[i][!j];
					return;
				}
			}
		}
	}
}


void csg_cube::unset_intersecting_edge_flags(coll_obj &cobj) const {

	assert(cobj.type == COLL_CUBE);

	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			if ((cobj.d[i][!j] >= d[i][0]-TOLER) && (cobj.d[i][!j] <= d[i][1]+TOLER)) { // cobj edge contained in cube in this dim
				bool overlaps(1);

				for (unsigned k = 0; k < 3 && !overlaps; ++k) {
					if (k != i && (d[k][0] >= cobj.d[k][1] || d[k][1] <= cobj.d[k][0])) overlaps = 0;
				}
				if (overlaps) {cobj.cp.surfs &= ~EFLAGS[i][!j];}
			}
		}
	}
}


// *** CSG ALGORITHM CODE ***


void coll_obj_group::remove_invalid_cobjs() {

	coll_obj_group cobjs2;

	for (size_t i = 0; i < size(); ++i) { // create new shapes vector with bad shapes removed
		if ((*this)[i].type != COLL_INVALID) cobjs2.push_back((*this)[i]);
	}
	this->swap(cobjs2);
}


void coll_obj_group::check_cubes() {

	if (!CHECK_COBJS) return;

	for (size_t i = 0; i < size(); ++i) {
		if ((*this)[i].type != COLL_CUBE) continue;
		csg_cube const cube((*this)[i]);
		
		if (cube.is_zero_area()) {
		  cout << "Zero area cube: "; cube.print(); cout << endl;
		  assert(0);
		}
	}
}


// Note: also sorts by alpha so that transparency works correctly
void coll_obj_group::merge_cubes() { // only merge compatible cubes

	if (!MERGE_COBJS) return;
	RESET_TIME;
	float const tolerance(-X_SCENE_SIZE*1.0E-6); // tiny negative tolerance to include adjacencies
	unsigned const ncobjs((unsigned)size());
	unsigned merged(0);
	cobj_bvh_tree cube_tree(this, 0, 0, 0, 1, 0); // cubes only
	cube_tree.add_cobjs(0);
	vector<unsigned> cids;

	for (unsigned i = 0; i < ncobjs; ++i) { // choose merge candidates
		if ((*this)[i].type != COLL_CUBE) continue;
		csg_cube cube((*this)[i]);
		if (cube.is_zero_area()) continue;
		cids.resize(0);
		cube_tree.get_intersecting_cobjs(cube, cids, i, tolerance, 0, -1);
		unsigned mi(0);

		for (vector<unsigned>::const_iterator it = cids.begin(); it != cids.end(); ++it) {
			unsigned const j(*it);
			assert(j < ncobjs && j != i);
			assert((*this)[j].type == COLL_CUBE);
			if (!(*this)[i].equal_params((*this)[j])) continue; // not compatible
			csg_cube cube2((*this)[j]);

			if (cube.cube_merge(cube2)) {
				(*this)[j].type = COLL_INVALID; // remove old coll obj
				++mi;
			}
		}
		if (mi > 0) { // cube has changed
			cube.write_to_cobj((*this)[i]);
			merged += mi;
			--i; // force this cube to be processed again
		}
	}
	if (merged > 0) remove_invalid_cobjs();
	cout << ncobjs << " => " << size() << endl;
	PRINT_TIME("Cube Merge");
}


void coll_obj_group::remove_overlapping_cubes(int min_split_destroy_thresh) { // objects specified later are the ones that are split/removed

	if (!UNOVERLAP_COBJS || empty()) return;
	RESET_TIME;
	unsigned const ncobjs((unsigned)size());
	vector<pair<unsigned, unsigned> > proc_order;
		
	for (unsigned i = 0; i < ncobjs; ++i) {
		if ((*this)[i].type == COLL_CUBE && (*this)[i].destroy >= min_split_destroy_thresh) {
			proc_order.push_back(make_pair((*this)[i].id, i));
		}
	}
	if (proc_order.empty()) return; // nothing to do
	sort(proc_order.begin(), proc_order.end());
	float const tolerance(X_SCENE_SIZE*1.0E-6); // tiny tolerance to prevent adjacencies
	cobj_bvh_tree cube_tree(this, 0, 0, 0, 1, 0); // cubes only
	cube_tree.add_cobjs(0);
	bool overlaps(0);
	coll_obj_group cur_cobjs, next_cobjs;
	vector<unsigned> cids;

	for (vector<pair<unsigned, unsigned> >::const_reverse_iterator it = proc_order.rbegin(); it != proc_order.rend(); ++it) {
		unsigned const i(it->second);
		csg_cube const cube((*this)[i]); // remove all other cobjs from cobjs[i] with lower id
		if (cube.is_zero_area()) continue;
		bool const neg((*this)[i].status == COLL_NEGATIVE);
		cids.resize(0);
		cube_tree.get_intersecting_cobjs(cube, cids, i, tolerance, 0, -1);
		if (cids.empty()) continue;
		cur_cobjs.resize(0);
		cur_cobjs.push_back((*this)[i]); // start with the current cobj
		bool was_removed(0);

		for (vector<unsigned>::const_iterator it = cids.begin(); it != cids.end(); ++it) {
			unsigned const j(*it);
			assert(j < size());
			assert((*this)[j].type == COLL_CUBE && j != i);
			if ((*this)[i].id < (*this)[j].id)              continue; // enforce ordering
			if (neg ^ ((*this)[j].status == COLL_NEGATIVE)) continue; // sign must be the same
			csg_cube sub_cube((*this)[j]);

			for (coll_obj_group::const_iterator c = cur_cobjs.begin(); c != cur_cobjs.end(); ++c) {
				if (sub_cube.subtract_from_cube(next_cobjs, *c)) {
					was_removed = overlaps = 1;
				}
				else { // didn't overlap
					next_cobjs.push_back(*c);
				}
			}
			cur_cobjs.clear();
			cur_cobjs.swap(next_cobjs);
		} // for it
		if (was_removed) {
			copy(cur_cobjs.begin(), cur_cobjs.end(), back_inserter(*this));
			(*this)[i].type = COLL_INVALID; // remove old coll obj
		}
		else {
			assert(cur_cobjs.size() == 1); // the original cobjs[i]
		}
	} // for i
	if (overlaps) remove_invalid_cobjs();
	cout << ncobjs << " => " << size() << endl;
	PRINT_TIME("Cube Overlap Removal");
}


// **********************************************


bool coll_obj::subtract_from_cobj(coll_obj_group &new_cobjs, csg_cube const &cube, bool include_polys) {

	bool removed(0);

	if (type == COLL_CUBE) {
		removed = cube.subtract_from_cube(new_cobjs, *this);
		if (!removed) cube.unset_adjacent_edge_flags(*this); // check adjacency and possibly remove some edge flags

		for (unsigned i = 0; i < new_cobjs.size(); ++i) {
			cube.unset_adjacent_edge_flags(new_cobjs[i]); // is this necessary?
		}
	}
	else if (cube.contains_cube(*this)) {
		removed = 1; // completely contained
	}
	else if (is_cylinder()) {
		removed = cube.subtract_from_cylinder(new_cobjs, *this);
	}
	else if (include_polys && type == COLL_POLYGON) {
		if (is_thin_poly()) {
			removed = cube.subtract_from_polygon(new_cobjs, *this);
		}
		else {
			removed = cube.subtract_from_thick_polygon(new_cobjs, *this); // implemented but not yet finished
		}
	}
	return removed;
}


void coll_obj_group::process_negative_shapes() { // negtive shapes should be non-overlapping

	RESET_TIME;
	unsigned const orig_ncobjs((unsigned)size());
	unsigned neg(0);
	coll_obj_group new_cobjs;

	for (unsigned i = 0; i < orig_ncobjs; ++i) { // find a negative cobj
		if ((*this)[i].status != COLL_NEGATIVE) continue;

		if ((*this)[i].type != COLL_CUBE) {
			std::cerr << "Only negative cube shapes are supported." << endl;
			exit(1);
		}
		unsigned ncobjs((unsigned)size()); // so as not to retest newly created subcubes
		csg_cube cube((*this)[i]); // the negative cube
		if (cube.is_zero_area()) continue;

		for (unsigned j = 0; j < ncobjs; ++j) { // find a positive cobj
			if (j != i && (*this)[j].status != COLL_NEGATIVE) {
				if (ONLY_SUB_PREV_NEG && (*this)[i].id < (*this)[j].id) continue; // positive cobj after negative cobj

				if ((*this)[j].subtract_from_cobj(new_cobjs, cube, 0)) {
					if (!new_cobjs.empty()) { // coll cube can be reused
						(*this)[j] = new_cobjs.back();
						new_cobjs.pop_back();
					}
					else {
						(*this)[j].type = COLL_INVALID; // remove old coll obj
					}
					copy(new_cobjs.begin(), new_cobjs.end(), back_inserter(*this)); // add in new fragments
					new_cobjs.clear();
				}
			}
		}
		(*this)[i].type = COLL_INVALID; // remove the negative cube since it is no longer needed
		++neg;
	}
	if (neg > 0) remove_invalid_cobjs();
	cout << orig_ncobjs << " => " << size() << endl;
	PRINT_TIME("Negative Shape Processing");
}


unsigned get_closest_val_index(float val, vector<double> const &sval) {

	for (unsigned i = 0; i < sval.size(); ++i) { // inefficient, assumes sval is small
		if (fabs(val - sval[i]) < TOLER) return i;
	}
	assert(0);
	return 0;
}


// split T-junctions of cubes in the same group
void coll_obj_group::subdiv_cubes() {

	if (!REMOVE_T_JUNCTIONS) return;
	RESET_TIME;
	unsigned ncobjs((unsigned)size()), num_remove(0);
	map<int, vector<unsigned> > id_map; // id to cobj indices map

	for (unsigned i = 0; i < ncobjs; ++i) {
		if ((*this)[i].type == COLL_INVALID || (*this)[i].type != COLL_CUBE) continue;
		if (REMOVE_T_JUNCTIONS == 1 && (*this)[i].counter != OBJ_CNT_REM_TJ) continue;
		id_map[(*this)[i].id].push_back(i);
	}
	for (map<int, vector<unsigned> >::const_iterator i = id_map.begin(); i != id_map.end(); ++i) {
		vector<unsigned> const &v(i->second);
		if (v.size() == 1) continue; // nothing to do
		set   <double> splits[3]; // x, y, z
		vector<double> svals [3]; // x, y, z

		for (unsigned j = 0; j < v.size(); ++j) {
			coll_obj const &c((*this)[v[j]]);

			for (unsigned d = 0; d < 3; ++d) {
				for (unsigned e = 0; e < 2; ++e) {
					splits[d].insert(c.d[d][e]);
				}
			}
		}
		for (unsigned d = 0; d < 3; ++d) {
			for (set<double>::const_iterator s = splits[d].begin(); s != splits[d].end(); ++s) {
				if (svals[d].empty() || (*s - svals[d].back()) > TOLER) svals[d].push_back(*s); // skip elements that are near equal
			}
			assert(svals[d].size() > 1);
		}
		for (unsigned j = 0; j < v.size(); ++j) {
			coll_obj const &c((*this)[v[j]]);
			unsigned bounds[3][2], tot_parts(1);

			for (unsigned d = 0; d < 3; ++d) {
				for (unsigned e = 0; e < 2; ++e) {
					bounds[d][e] = get_closest_val_index(c.d[d][e], svals[d]);
				}
				assert(bounds[d][0] < bounds[d][1] && bounds[d][1] < svals[d].size());
				tot_parts *= (bounds[d][1] - bounds[d][0]);
			}
			assert(tot_parts > 0);
			if (tot_parts == 1) continue; // no splits required
				
			for (unsigned x = bounds[0][0]; x < bounds[0][1]; ++x) {
				for (unsigned y = bounds[1][0]; y < bounds[1][1]; ++y) {
					for (unsigned z = bounds[2][0]; z < bounds[2][1]; ++z) {
						unsigned const xyz[3] = {x, y, z};
						push_back((*this)[v[j]]);

						for (unsigned d = 0; d < 3; ++d) {
							assert(xyz[d]+1 < svals[d].size());

							for (unsigned e = 0; e < 2; ++e) {
								back().d[d][e] = svals[d][xyz[d]+e];
							}
						}
					}
				}
			}
			(*this)[v[j]].type = COLL_INVALID;
			++num_remove;
		} // for j
	} // for i
	if (num_remove > 0) {remove_invalid_cobjs();}
	cout << ncobjs << " => " << size() << endl;
	PRINT_TIME("Subdiv Cubes");
}


bool comp_cobjs_by_draw_params(coll_obj const &a, coll_obj const &b) {
	if (a.cp.tid      < b.cp.tid)      return 1;
	if (b.cp.tid      < a.cp.tid)      return 0;
	if (a.group_id    < b.group_id)    return 1;
	if (b.group_id    < a.group_id)    return 0;
	if (a.cp.color    < b.cp.color)    return 1;
	if (b.cp.color    < a.cp.color)    return 0;
	if (a.cp.specular < b.cp.specular) return 1;
	if (b.cp.specular < a.cp.specular) return 0;
	if (a.cp.shine    < b.cp.shine)    return 1;
	if (b.cp.shine    < a.cp.shine)    return 0;
	if (a.type        < b.type)        return 1;
	if (b.type        < a.type)        return 0;
	return (a.points[0] < b.points[0]);
}

void coll_obj_group::sort_cobjs_for_rendering() {
	std::sort(begin(), end(), comp_cobjs_by_draw_params);
}


color_tid_vol::color_tid_vol(coll_obj const &cobj, float volume_, float thickness_, bool ua)
	: cid(cobj.id), tid(cobj.cp.tid), destroy(cobj.destroy), draw(cobj.cp.draw), unanchored(ua), is_2d(cobj.is_thin_poly()),
	volume(volume_), thickness(thickness_), tscale(cobj.cp.tscale), color(cobj.cp.color)
{
	if (cobj.type == COLL_CUBE && cobj.cp.light_atten > 0.0) {
		color.alpha += (1.0 - color.alpha)*(1.0 - exp(-cobj.cp.light_atten*thickness));
	}
	if (tscale == 0.0) { // calculate tscale from object size
		tscale = 3.0/(fabs(cobj.d[0][1] - cobj.d[0][0]) + fabs(cobj.d[1][1] - cobj.d[1][0]) + fabs(cobj.d[2][1] - cobj.d[2][0]));
	}
	copy_from(cobj);
}



