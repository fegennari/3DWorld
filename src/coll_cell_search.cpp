// 3D World - OpenGL CS184 Computer Graphics Project - collision detection coll cell iteration and search
// by Frank Gennari
// 2/6/06

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"


bool const CACHE_COBJ_LITES   = 0;
bool const CACHE_OCCLUDER     = 1;
unsigned const QLP_CACHE_SIZE = 10000;
float const GET_OCC_EXPAND    = 0.02;


int cobj_counter(0);

extern int coll_border, display_mode;
extern float czmin, czmax, zmin, zbottom, water_plane_z;
extern vector<coll_obj> coll_objects;


// *** Test BSP Tree of Static Collision Objects ***

class cobj_bsp_tree {
	struct bsp_node {
		float d[3][2]; // bbox
		unsigned b[3]; // {left, right, center} branches
		unsigned start, end; // index into cixs for leaves

		bsp_node() : start(0), end(0) {
			UNROLL_3X(d[i_][0] = d[i_][1] = 0.0;)
			UNROLL_3X(b[i_] = 0;) // kids of 0 means empty kids
		}
	};

	struct coll_state {
		point const &p1, p2;
		point &cpos;
		vector3d &cnorm;
		int &cindex;
		int ignore_cobj;
		float tmin, tmax;

		coll_state(point const &p1_, point const &p2_, point &cp, vector3d &cn, int &ci, int ic)
			: p1(p1_), p2(p2_), cpos(cp), cnorm(cn), cindex(ci), ignore_cobj(ic), tmin(0.0), tmax(1.0) {}
		void update_cpos() {cpos = p1 + (p2 - p1)*tmax;}
	};

	vector<coll_obj> const &cobjs;
	vector<bsp_node> nodes;
	vector<unsigned> cixs;
	unsigned root_node;

	inline coll_obj const &get_cobj(unsigned ix) const {
		assert(ix < cixs.size() && cixs[ix] < cobjs.size());
		return cobjs[cixs[ix]];
	}

	// Note: start by adding all cobjs to [start, end], then calc bbox, then split into branches and leaves
	void calc_node_bbox(bsp_node &n) const {
		assert(n.start < n.end);

		for (unsigned i = n.start; i < n.end; ++i) { // check contained
			if (i == 0) {
				copy_cube_d(get_cobj(i).d, n.d);
			}
			else {
				UNROLL_3X(n.d[i_][0] = min(n.d[i_][0], get_cobj(i).d[i_][0]);)
				UNROLL_3X(n.d[i_][1] = max(n.d[i_][1], get_cobj(i).d[i_][1]);)
			}
		}
	}

	void build_tree(unsigned nix) {
		assert(nix < nodes.size());
		bsp_node &n(nodes[nix]);
		assert(n.start < n.end);
		if ((n.end - n.start) == 1) return; // base case, one kid
		// write
	}

	bool test_cobj(point const &p1, point const &p2, coll_state &state, unsigned ix) const {
		float t(0.0);
		bool const ret(get_cobj(ix).line_int_exact(p1, p2, t, state.cnorm, state.tmin, state.tmax));
		if (ret) state.tmax = t;
		return ret;
	}

	bool search_tree(point p1, point p2, coll_state &state, unsigned nix) const {
		assert(nix < nodes.size());
		bsp_node const &n(nodes[nix]);
		assert(n.start < n.end);
		// write - clip p1 and p2
		bool ret(0);

		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			ret |= test_cobj(p1, p2, state, i);
		}

		// FIXME: interate in correct order based on (p2 - p1)
		for (unsigned i = 0; i < 3; ++i) { // check branches
			if (n.b[i] == 0) continue; // empty branch
			// write - bbox test
			ret |= search_tree(p1, p2, state, n.b[i]);
		}
		return ret;
	}

public:
	cobj_bsp_tree(vector<coll_obj> const &cobjs_) : cobjs(cobjs_), root_node(0) {}

	void clear() {
		cixs.clear();
		nodes.clear();
		root_node = 0;
	}

	void add_cobjs() {
		RESET_TIME;
		clear();
		cixs.reserve(cobjs.size()); // ???

		for (vector<coll_obj>::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
			if (i->status != COLL_STATIC) continue;
			//if (i->counter == cobj_counter) continue;
			cixs.push_back(i - cobjs.begin());
		}
		if (cixs.empty()) return; // nothing to be done
		bsp_node n;
		n.start = 0;
		n.end   = cixs.size();
		calc_node_bbox(n);
		root_node = nodes.size();
		nodes.push_back(n);
		build_tree(root_node);
		PRINT_TIME("Cobj BSP Tree Create");
	}

	bool check_coll_line_exact(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj) const {
		if (nodes.empty()) return 0;
		coll_state state(p1, p2, cpos, cnorm, cindex, ignore_cobj);
		bool const ret(search_tree(p1, p2, state, root_node));
		if (ret) state.update_cpos();
		return ret;
	}
};


// *** End Test Code ***



// false intersections are OK, used only for shadow calculations and such things
// cobjs == NULL runs in a special mode
bool check_vert_collision_sphere(point const &pos, float radius, int skip_dynamic, bool trans_test, vector<int> *cobjs) {

	assert(radius >= 0.0); // is radius == 0 OK?
	bool any_coll(0);
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos)) return 0; // object along edge
	coll_cell const &cell(v_collision_matrix[ypos][xpos]);
	if (cell.cvals.empty()) return 0;
	float const z1(pos.z-radius), z2(pos.z+radius); // Note that the comparison below adds in an addition radius to account for zmin/zmax inaccuracy
	//if (skip_dynamic && radius < HALF_DXY && (z1 > cell.zmax+radius || z2 < cell.zmin-radius)) return 0;
	float const x1(pos.x-radius), x2(pos.x+radius), y1(pos.y-radius), y2(pos.y+radius);

	for (int k = (int)cell.cvals.size()-1; k >= 0; --k) { // iterate backwards
		int const index(cell.cvals[k]);
		if (index < 0) continue;
		assert(unsigned(index) < coll_objects.size());
		coll_obj const &cobj(coll_objects[index]);
		if (cobjs != NULL && cobj.counter == cobj_counter) continue; // already seen

		if ((skip_dynamic && cobj.status == COLL_DYNAMIC) || cobj.no_collision() || (skip_dynamic >= 2 && !cobj.cp.draw) ||
			(skip_dynamic >= 3 && !cobj.fixed) || (trans_test && cobj.cp.color.alpha < 0.5))
		{
			coll_objects[index].counter = cobj_counter;
			continue;
		}
		if (z2 < cobj.d[2][0] || z1 > cobj.d[2][1]) continue;
		if (x2 < cobj.d[0][0] || x1 > cobj.d[0][1] || y2 < cobj.d[1][0] || y1 > cobj.d[1][1]) continue;
		int coll(0);

		switch (cobj.type) { // within bounding box of collision object
		case COLL_CUBE:
			if (skip_dynamic >= 2 && cobj.cp.surfs == EF_ALL) break; // all sides hidden
			if (!sphere_cube_intersect(pos, radius, cobj.d))  break;
			coll = 1;
			break;
		case COLL_CYLINDER:
			coll = dist_xy_less_than(pos, cobj.points[0], (cobj.radius+radius));
			break;
		case COLL_SPHERE:
			coll = dist_less_than(pos, cobj.points[0], (cobj.radius+radius));
			break;
		case COLL_CYLINDER_ROT:
			coll = sphere_intersect_cylinder(pos, radius, cobj.points[0], cobj.points[1], cobj.radius, cobj.radius2);
			//if (coll && z1 > cell.zmax+radius) {} // soemthing bad happened
			break;
		case COLL_POLYGON: // must be coplanar
			coll = sphere_ext_poly_intersect(cobj.points, cobj.npoints, cobj.norm, pos, radius, cobj.thickness, MIN_POLY_THICK2);
			break;
		} // switch
		if (coll) {
			any_coll = 1;
			
			if (cobjs != NULL) {
				coll_objects[index].counter = cobj_counter;
				cobjs->push_back(index);
			}
			else {
				return 1;
			}
		}
	} // for cl
	return any_coll;
}


// returns 1 if there is no intersection
bool coll_obj::cobj_plane_side_test(point const *pts, unsigned npts, point const &lpos) const {

	vector<vector<point> > const *ppts(NULL); // used as pointer to static variable

	for (unsigned i = 0; i < npts; ++i) {
		point const spts[3] = {pts[i], pts[(i+1)%npts], lpos};
		point const center(get_center(spts, 3));
		vector3d const pts_norm(get_poly_norm(spts));
		point pt;

		if (type == COLL_POLYGON) {
			if (thickness > MIN_POLY_THICK2) { // thick polygon
				if (!ppts) ppts = &thick_poly_to_sides(points, npoints, norm, thickness);

				for (unsigned i = 0; i < ppts->size(); ++i) {
					for (unsigned j = 0; j < (*ppts)[i].size(); ++j) {
						if (dot_product_ptv(pts_norm, center, (*ppts)[i][j]) > 0.0) return 0;
					}
				}
			}
			else { // thin polygon
				for (unsigned i = 0; i < unsigned(npoints); ++i) {
					if (dot_product_ptv(pts_norm, center, points[i]) > 0.0) return 0;
				}
			}
		}
		else { // use bounding cube (bad for pine tree cones)
			for (unsigned x = 0; x < 2; ++x) {
				pt[0] = d[0][x];
				for (unsigned y = 0; y < 2; ++y) {
					pt[1] = d[1][y];
					for (unsigned z = 0; z < 2; ++z) {
						pt[2] = d[2][z];
						if (dot_product_ptv(pts_norm, center, pt) > 0.0) return 0;
					}
				}
			}
		}
	}
	return 1;
}


// Note: used _only_ for static cobj shadows
bool coll_pt_vis_test_large(point pos, point pos2, vector<int> &cobjs, int cobj, float radius, int skip_dynamic,
							bool trans_test, point const *pts, unsigned npts, point const &lpos)
{
	assert(radius > 0.0);
	cobjs.resize(0);
	if (!do_line_clip_scene(pos, pos2, max(zbottom, czmin)-radius, czmax+radius)) return 0; // assumes pos is in the simulation region
	radius *= (2.0/SQRT3); // account for sphere projection step-by-step instead of true cylinder projection
	radius  = max(0.25f*HALF_DXY, radius); // lower bound so that nsteps isn't excessive
	vector3d vcf(pos2, pos);
	float const mag(1.1*vcf.mag()/radius); // stepsize == radius, multiply by 1.1 for larger overlap
	unsigned const nsteps(unsigned(mag + 0.5) + 1); // ceil? - should really subtract radius from pos
	assert(nsteps > 0);
	if (mag > TOLERANCE) vcf /= mag; // normalize
	++cobj_counter;
	assert(cobj < 0 || size_t(cobj) < coll_objects.size());
	if (cobj >= 0) coll_objects[cobj].counter = cobj_counter;
	
	for (unsigned i = 0; i < nsteps; ++i) {
		static vector<int> cand_cobjs;
		cand_cobjs.resize(0);
		check_vert_collision_sphere(pos, radius, skip_dynamic, trans_test, ((pts == NULL) ? &cobjs : &cand_cobjs));

		if (pts != NULL) {
			for (unsigned j = 0; j < cand_cobjs.size(); ++j) { // check candidate cobjs
				int const cindex(cand_cobjs[j]);
				assert(cindex >= 0 && cindex < (int)coll_objects.size() && cindex != cobj);
				coll_obj &c(coll_objects[cindex]);
				if (c.dynamic_shadows_only()) continue; // cobj has dynamic shadows only
				int region(~0);
				point center;
				float brad;
				c.bounding_sphere(center, brad);

				for (unsigned i = 0; i < npts && region != 0; ++i) {
					float const dist(p2p_dist(pts[i], center) + brad);
					vector3d const ldir((pts[i] - lpos).get_norm());
					region &= get_region( pts[i],              c.d);
					region &= get_region((pts[i] - ldir*dist), c.d);
				}
				if (region != 0 || c.cobj_plane_side_test(pts, npts, lpos)) continue; // all outside a plane of the view volume
				
				if (c.radius > radius) { // check for complete shadowing
					int vshadowed(1);

					for (unsigned i = 0; i < npts && vshadowed; ++i) {
						if (!c.line_intersect(pts[i], lpos)) vshadowed = 0;
					}
					if (vshadowed) { // completely contained in the shadow of this object
						c.shadow_depends.insert(cobj);
						return 1; // could continue to find other completely shadowing objects, but > 1 is unlikely
					}
				}
				cobjs.push_back(cindex);
			} // for j
		} // pts != NULL
		if (mag <= TOLERANCE) return 0;
		pos += vcf;
	} // for i
	return 0;
}


// false negatives are OK except when called from check_coll_line()
bool coll_obj::line_intersect(point const &p1, point const &p2) const {

	if (!check_line_clip(p1, p2, d)) return 0;
	
	switch (type) {
		case COLL_CUBE:
			return 1;
		case COLL_SPHERE:
			return line_sphere_intersect(p1, p2, points[0], radius);
		case COLL_CYLINDER:
		case COLL_CYLINDER_ROT:
			return line_intersect_cylinder(p1, p2, cylinder_3dw(points[0], points[1], radius, radius2), !(cp.surfs & 1));
		case COLL_POLYGON: // must be coplanar
			{
				assert(npoints >= 3);
				vector3d const v1(p2, p1);

				if (thickness > MIN_POLY_THICK2) { // test extruded (3D) polygon
					static vector<point> pts[2];
					gen_poly_planes(points, npoints, norm, thickness, pts);
					bool const test_side(dot_product(v1, norm) > 0.0);
					if (thick_poly_intersect(v1, p1, norm, pts, test_side, npoints)) return 1;
				}
				else { // test planar (2D) polygon
					if (line_poly_intersect(v1, p1, points, npoints, norm)) return 1;
				}
			}
			break;
	}
	return 0;
}


bool coll_obj::line_int_exact(point const &p1, point const &p2, float &t, vector3d &cnorm, float tmin, float tmax) const {

	float clip_tmin, clip_tmax;
	if (!get_line_clip(p1, p2, d, clip_tmin, clip_tmax) || clip_tmin > tmax || clip_tmax < tmin) return 0;
	vector3d v1(p2 - p1);
	
	switch (type) {
		case COLL_CUBE:
			{
				t = clip_tmin;
				if (t > tmax || t < tmin) return 0;
				get_closest_cube_norm(d, (p1 + v1*t), cnorm);
				return 1;
			}
		case COLL_SPHERE:
			{
				point coll_pos;
				v1.normalize();
				if (!line_sphere_int(v1, p1, points[0], radius, coll_pos, 0)) return 0;
				t = -1.0; // start at a bad value
				
				for (unsigned i = 0; i < 3; ++i) {
					if (fabs(p2[i] - p1[i]) > TOLERANCE) {
						t = (coll_pos[i] - p1[i])/(p2[i] - p1[i]);
						break;
					}
				}
				if (t > tmax || t < tmin) return 0;
				cnorm = (coll_pos - points[0]);
				if (!cnorm.normalize_test()) cnorm = plus_z; // arbitrary
				return 1;
			}
		case COLL_CYLINDER:
		case COLL_CYLINDER_ROT:
			{
				int const int_type(line_int_thick_cylinder(p1, p2, points[0], points[1], 0.0, 0.0, radius, radius2, 1, t));
				if (!int_type || t > tmax || t < tmin) return 0;
				
				if (int_type == 1) { // side intersection
					vector3d const cv(points[0] - points[1]);
					orthogonalize_dir((p1 - points[0]), cv, cnorm, 0);
					if (cnorm == zero_vector) orthogonalize_dir((p2 - points[0]), cv, cnorm, 0); // p1 is bad, so try p2

					if (radius != radius2) {
						if (!cnorm.normalize_test()) cnorm = plus_z; // arbitrary
						float const len(cv.mag());

						if (len > TOLERANCE) {
							float const dr(radius2 - radius), denom(sqrt(len*len + dr*dr));
							assert(denom > TOLERANCE);
							cnorm *= len/denom;
							cnorm += cv*(dr*len/denom);
						}
					}
				}
				else { // top/bottom intersection
					cnorm = (points[1] - points[0]);
					if ((int_type == 2) ^ (radius >= radius2)) cnorm.negate(); // r1 >= r2 => swap
				}
				if (!cnorm.normalize_test()) cnorm = plus_z; // arbitrary
				return 1;
			}
		case COLL_POLYGON: // must be coplanar
			{
				assert(npoints >= 3);

				if (thickness > MIN_POLY_THICK2) { // test extruded (3D) polygon
					t = 2.0; // start at a bad value
					float tval;
					static vector<point> pts[2];
					gen_poly_planes(points, npoints, norm, thickness, pts);
					unsigned const test_side(dot_product(v1, norm) > 0.0);
					point const *const points2(&(pts[test_side].front()));
					
					if (line_poly_intersect(p1, p2, points2, npoints, norm, tval) && (tval <= tmax && tval >= tmin)) {
						t     = tval;
						cnorm = get_poly_dir_norm(norm, p1, v1, t);
					}
					for (int j = 0; j < npoints; ++j) { // now test the <npoints> sides
						unsigned const jnext((j+1)%npoints);
						point const side_pts[4] = {pts[0][j], pts[0][jnext], pts[1][jnext], pts[1][j]};
						vector3d const side_norm(get_poly_norm(side_pts));
					
						if (line_poly_intersect(p1, p2, side_pts, 4, side_norm, tval)) {
							if (tval < t && (tval <= tmax && tval >= tmin)) {
								t     = tval;
								cnorm = get_poly_dir_norm(side_norm, p1, v1, t);
							}
						}
					}
					return (t <= tmax && t >= tmin);
				}
				if (!line_poly_intersect(p1, p2, points, npoints, norm, t)) return 0;
				if (t > tmax || t < tmin) return 0;
				cnorm = get_poly_dir_norm(norm, p1, v1, t);
				return 1;
			}
	}
	return 0;
}


struct base_intersector {

	point const &pos1, &pos2;
	float z1, z2;

	base_intersector(point const &pos1_, point const &pos2_)
		: pos1(pos1_), pos2(pos2_), z1(min(pos1.z, pos2.z)), z2(max(pos1.z, pos2.z)) {}
	bool test_cell_end(int xpos, int ypos) const {return 0;}
	bool had_intersection() const {return 0;}
};


// test_alpha: 0 = allow any alpha value, 1 = require alpha = 1.0, 2 = get interseced cobj with max alpha
class line_intersector : public base_intersector {

	int &cindex, skip_dynamic, test_alpha;
	float max_alpha;

public:
	line_intersector(point const &pos1_, point const &pos2_, int &cindex_, int skip_dynamic_, int test_alpha_)
		: base_intersector(pos1_, pos2_), cindex(cindex_), skip_dynamic(skip_dynamic_), test_alpha(test_alpha_), max_alpha(0.0) {}
	bool test_cell(coll_cell const &cell, int xpos, int ypos) const {return (!skip_dynamic || (z1 <= cell.zmax && z2 >= cell.zmin));}

	int proc_cobj(coll_cell const &cell, int const index) {
		coll_obj &cobj(coll_objects[index]);
		cobj.counter = cobj_counter;
		if (skip_dynamic && cobj.status == COLL_DYNAMIC)         return 2;
		if (skip_dynamic >= 2 && !cobj.cp.draw)                  return 2;
		if (test_alpha == 1 && cobj.is_semi_trans())             return 2; // semi-transparent, can see through
		if (test_alpha == 2 && cobj.cp.color.alpha <= max_alpha) return 2;
		if (cobj.status == COLL_STATIC && (z1 > cell.zmax || z2 < cell.zmin)) return 0;
		if (z1 > cobj.d[2][1] || z2 < cobj.d[2][0])              return 2; // clip this shape
		
		if (cobj.line_intersect(pos1, pos2)) {
			cindex    = index;
			if (test_alpha != 2) return 1;
			max_alpha = cobj.cp.color.alpha; // we need all intersections to find the max alpha
		}
		return 2;
	}
};


class line_intersector_exact : public base_intersector {

	point &cpos;
	vector3d &cnorm;
	int &cindex;
	bool test_alpha;
	int ignore_cobj;
	float splash_val, t, tmin;
	bool splash;
	point splash_pos;

public:
	line_intersector_exact(point const &pos1_, point const &pos2_, point &cpos_, vector3d &cnorm_,
		int &cindex_, float splash_val_, int ignore_cobj_, bool test_alpha_)
		: base_intersector(pos1_, pos2_), cpos(cpos_), cnorm(cnorm_), cindex(cindex_), test_alpha(test_alpha_),
		ignore_cobj(ignore_cobj_), splash_val(splash_val_), tmin(1.0), splash(0) {cindex = -1;}

	bool test_cell(coll_cell const &cell, int xpos, int ypos) { // always returns 1
		if (splash_val == 0.0 || splash) return 1;//(z1 <= cell.zmax && z2 >= cell.zmin);
		return test_cell_splash(xpos, ypos);
	}

	bool test_cell_splash(int xpos, int ypos) {
		if (!mesh_is_underwater(xpos, ypos)) return 1;
		float const wmz(water_matrix[ypos][xpos]);
		if ((pos1.z < wmz) ^ (pos2.z > wmz)) return 1;
		float const tz((wmz - pos1.z)/(pos2.z - pos1.z));
		splash_pos = pos1 + (pos2 - pos1)*tz;
		float const mx(get_xval(xpos)), my(get_yval(ypos));
		
		if (splash_pos.x > (mx-DX_VAL) && splash_pos.x < (mx+DX_VAL) && splash_pos.y > (my-DY_VAL) && splash_pos.y < (my+DY_VAL)) {
			add_splash(xpos, ypos, 25.0, 0.01); // dynamic water
			draw_splash(splash_pos.x, splash_pos.y, (wmz + 0.0001), splash_val);
			splash = 1;
		}
		return 1;
	}

	void finish() {
		if (splash) gen_line_of_bubbles(splash_pos, (had_intersection() ? cpos : pos2));
	}
	bool test_cell_end(int xpos, int ypos) const {return (cindex >= 0 && get_xpos(cpos.x) == xpos && get_ypos(cpos.y) == ypos);}
	bool had_intersection() const {return (cindex >= 0);}

	int proc_cobj(coll_cell const &cell, int const index) {
		if (index == ignore_cobj) return 2;
		coll_obj &cobj(coll_objects[index]);
		if (cobj.status == COLL_STATIC && (z1 > cell.zmax || z2 < cell.zmin)) return 0;
		cobj.counter = cobj_counter;
		if (z1 > cobj.d[2][1] || z2 < cobj.d[2][0]) return 2; // clip this shape
		if (test_alpha && cobj.is_semi_trans())     return 2; // semi-transparent, can see through
		if (cobj.cp.color.alpha == 0.0)             return 2; // transparent - what objects are these?
		vector3d temp_norm;

		if (had_intersection()) { // minor optimization (if z1/z2 have been updated, then don't need to check z)
			if (max(pos1.x, cpos.x) < cobj.d[0][0] || min(pos1.x, cpos.x) > cobj.d[0][1]) return 2;
			if (max(pos1.y, cpos.y) < cobj.d[1][0] || min(pos1.y, cpos.y) > cobj.d[1][1]) return 2;
		}
		if (cobj.line_int_exact(pos1, pos2, t, temp_norm, 0.0, tmin)) {
			cindex = index;
			cnorm  = temp_norm;
			tmin   = t;
			cpos   = pos1 + (pos2 - pos1)*tmin;
			z1     = min(pos1.z, cpos.z); // update z1/z2: minor optimization
			z2     = max(pos1.z, cpos.z);
		}
		return 2;
	}
};


bool is_contained(point const &pos, point const *const pts, unsigned npts, float const d[3][2]) {

	for (unsigned i = 0; i < npts; ++i) { // can almost skip two corners on a quad
		if (!check_line_clip(pos, pts[i], d)) return 0;
	}
	return 1;
}


class line_intersector_occlusion_cobjs : public base_intersector { // tests for polygon occlusion

	point const *const pts;
	unsigned npts;
	int coll_cobj;

public:
	line_intersector_occlusion_cobjs(point const &pos1_, point const &pos2_, const point *pts_, unsigned npts_)
		: base_intersector(pos1_, pos2_), pts(pts_), npts(npts_), coll_cobj(-1) {}
	bool test_cell(coll_cell const &cell, int xpos, int ypos) const {return (z1 <= cell.occ_zmax && z2 >= cell.occ_zmin);}
	int get_coll_cobj() const {return coll_cobj;}

	int proc_cobj(coll_cell const &cell, int const index) {
		coll_obj &cobj(coll_objects[index]);
		cobj.counter = cobj_counter;
		if (!cobj.is_occluder() || z1 > cobj.d[2][1] || z2 < cobj.d[2][0]) return 2; // clip this shape
		
		if (is_contained(pos1, pts, npts, cobj.d))  {
			coll_cobj = index;
			return 1;
		}
		return 2;
	}
};


class line_intersector_get_cobjs : public base_intersector { // gets potential occluders

	vector<int> &cobjs;

public:
	line_intersector_get_cobjs(point const &pos1_, point const &pos2_, vector<int> &cobjs_)
		: base_intersector(pos1_, pos2_), cobjs(cobjs_) {}
	bool test_cell(coll_cell const &cell, int xpos, int ypos) const {return (z1 <= cell.occ_zmax && z2 >= cell.occ_zmin);}

	int proc_cobj(coll_cell const &cell, int const index) {
		coll_obj &cobj(coll_objects[index]);
		cobj.counter = cobj_counter;
		if (!cobj.is_occluder() || !cobj.fixed || cobj.volume < 0.001) return 2; // not an occluder
		if (z1 > cobj.d[2][1] || z2 < cobj.d[2][0])                    return 2; // clipped
		if (check_line_clip_expand(pos1, pos2, cobj.d, GET_OCC_EXPAND)) cobjs.push_back(index);
		return 2;
	}
};


class line_intersector_cylinder : public base_intersector { // cylinder intersect cobj bbox - actually uses line intersects grown bbox, approximate

	vector<int> &cobjs;
	float radius;
	int skip_dynamic, c_obj;
	bool test_lighted;

public:
	line_intersector_cylinder(point const &pos1_, point const &pos2_, vector<int> &cobjs_, int c_obj_, float r, int sd, bool tl)
		: base_intersector(pos1_, pos2_), cobjs(cobjs_), radius(r), skip_dynamic(sd), c_obj(c_obj_), test_lighted(tl)
	{
		assert(radius >= 0.0);
	}
	bool test_cell(coll_cell const &cell, int xpos, int ypos) const { // what about 2*radius hack?
		return ((z1-radius) <= cell.zmax && (z2+radius) >= cell.zmin);
	}
	int proc_cobj(coll_cell const &cell, int const index) {
		coll_obj &cobj(coll_objects[index]);
		cobj.counter = cobj_counter;
		// *** careful, this is untested and may be incorrect for uses other than shadow filtering ***
		if (test_lighted  && cobj.status == COLL_STATIC && cobj.lighted == COBJ_LIT_FALSE)         return 2;
		if ((skip_dynamic && cobj.status == COLL_DYNAMIC) || (skip_dynamic >= 2 && !cobj.cp.draw)) return 2;
		if (cobj.cp.surfs == EF_ALL || (z1-radius) > cobj.d[2][1] || (z2+radius) < cobj.d[2][0])   return 2; // clip this shape

		if (check_line_clip_expand(pos1, pos2, cobj.d, radius)) { // line intersects expanded cube => cylinder intersects cube (conservative)
			cobjs.push_back(index);
		}
		return 2;
	}
};

template<typename T> class coll_cell_line_iterator {

	T &lint;
	bool skip, fast, skip_dynamic;
	int c_obj;

	bool skip_this_index(int index) const {

		if (index < 0 || index == c_obj)                 return 1; // bad or skipped index
		assert(unsigned(index) < coll_objects.size());
		if (coll_objects[index].counter == cobj_counter) return 1; // already seen
		if (coll_objects[index].no_collision())          return 1; // disabled
		return 0;
	}

	bool cobj_test(int xpos, int ypos) const {

		// we occasionally get here with duplicate x/y values when the line is short, but that should be ok
		coll_cell const &cell(v_collision_matrix[ypos][xpos]);
		if (!lint.test_cell(cell, xpos, ypos) || cell.cvals.empty()) return 0; // clip entire coll cell
		bool const subdiv(!cell.cvz.empty());

		if (!subdiv || !skip_dynamic) {
			for (int k = int(cell.cvals.size()-1); k >= 0; --k) { // iterate backwards
				int const index(cell.cvals[k]);
				assert(unsigned(index) < coll_objects.size());
				if (skip_this_index(index)) continue;
				if (subdiv && coll_objects[index].status == COLL_STATIC) break; // done with dynamic objects, move on to static
				int const val(lint.proc_cobj(cell, index));
				if (val == 0 || val == 1) return (val != 0);
			}
		}
		if (subdiv) {
			point cp1(lint.pos1), cp2(lint.pos2);
			float const xv(get_xval(xpos)), yv(get_yval(ypos));
			float const d[3][2] = {{(xv-DX_VAL), (xv+DX_VAL)}, {(yv-DY_VAL), (yv+DY_VAL)}, {cell.zmin, cell.zmax}};
			if (!do_line_clip(cp1, cp2, d)) return 0;
			unsigned const sz(cell.cvz.size());
			float const val(sz/(cell.zmax - cell.zmin));
			unsigned zs(min(sz-1, (unsigned)max(0, (int)floor(val*(min(cp1.z, cp2.z) - cell.zmin)))));
			unsigned ze(min(sz-1, (unsigned)max(0, (int)floor(val*(max(cp1.z, cp2.z) - cell.zmin)))));
			bool const dir(cp1.z < cp2.z);
			int const di(dir ? 1 : -1);
			if (dir == 0) swap(zs, ze);
			
			for (unsigned i = zs; ; i += di) {
				vector<int> const &cvz(cell.cvz[i]);

				for (unsigned z = 0; z < cvz.size(); ++z) {
					if (skip_this_index(cvz[z])) continue;
					int const val(lint.proc_cobj(cell, cvz[z]));
					if (val == 0 || val == 1) return (val != 0);
				}
				if (lint.had_intersection() || i == ze) break; // had intersection at a closer zval, so we're done
			}
		}
		return lint.test_cell_end(xpos, ypos);
	}

public:
	coll_cell_line_iterator(T &lint_, bool s, bool sd, int cobj, bool f=0) :
	  lint(lint_), skip(s), fast(s || f), skip_dynamic(sd), c_obj(cobj) {}

	bool do_iter(float radius=0.0) const {

		bool const do_bnd_test(radius > 2.0*HALF_DXY); // very large object
		int const xa(get_xpos(lint.pos1.x)), ya(get_ypos(lint.pos1.y)), xb(get_xpos(lint.pos2.x)), yb(get_ypos(lint.pos2.y));
		++cobj_counter;

		if (!do_bnd_test && xa == xb && ya == yb) {
			if (point_outside_mesh(xa, ya)) return 0;
			return cobj_test(xa, ya);
		}
		int const dx(xb - xa), dy(yb - ya), steps(max(1, ((abs(dx) > abs(dy)) ? abs(dx): abs(dy))/max(1, coll_border)));
		double const xinc(dx/(double)steps), yinc(dy/(double)steps);
		double x(xa), y(ya);
		int last_x(-1), last_y(-1), skipval(1);
		bool first(1);
		int const cb(max(1, coll_border)); // not sure if this is correct
		int bnds[2][2] = {{MESH_X_SIZE,0}, {MESH_Y_SIZE,0}}; // {x,y}{min,max}

		for (int k = 0; k <= steps; ++k) { // DDA algorithm
			int const xpos(int(x + 0.5)), ypos(int(y + 0.5));
			x += xinc;
			y += yinc;
			// might skip some coll cells during diag strides but having at least one cell border around cobjs should fix it
			if (skip && cb > 0) { // can miss some collisions with tree leaves, which have no coll_border
				if (--skipval) continue; // skip this cell
				skipval = (cb << 1);
			}
			unsigned num(1);
			int xv[2] = {xpos, xpos}, yv[2] = {ypos, ypos};

			if (!fast && !first && xpos != last_x && ypos != last_y) {
				bool const x_stays((fabs(last_x - x) + fabs(ypos - y)) < (fabs(xpos - x) + fabs(last_y - y)));
				if (x_stays) xv[1] = last_x; else yv[1] = last_y;
				++num;
			}
			for (unsigned n = 0; n < num; ++n) {
				if (do_bnd_test) {
					bnds[0][0] = min(bnds[0][0], xv[n]);
					bnds[0][1] = max(bnds[0][1], xv[n]);
					bnds[1][0] = min(bnds[1][0], yv[n]);
					bnds[1][1] = max(bnds[1][1], yv[n]);
				}
				else if (!point_outside_mesh(xv[n], yv[n]) && cobj_test(xv[n], yv[n])) return 1;
			}
			first  = 0;
			last_x = xpos;
			last_y = ypos;
		} // for k
		if (do_bnd_test) { // inefficient, but what can we do? take bounds of all cube corner points?
			int const x1(max(0,           (bnds[0][0] - int(radius/DX_VAL))));
			int const y1(max(0,           (bnds[1][0] - int(radius/DY_VAL))));
			int const x2(min(MESH_X_SIZE, (bnds[0][1] + int(radius/DX_VAL))));
			int const y2(min(MESH_Y_SIZE, (bnds[1][1] + int(radius/DY_VAL))));

			for (int y = y1; y <= y2; ++y) {
				for (int x = x1; x <= x2; ++x) {
					if (!point_outside_mesh(x, y)) cobj_test(x, y);
				}
			}
		}
		return 0;
	}
};


bool check_coll_line(point pos1, point pos2, int &cindex, int cobj, int skip_dynamic, int test_alpha) {

	cindex = -1;
	if (!do_line_clip_scene(pos1, pos2, czmin, czmax)) return 0;
	line_intersector lint(pos1, pos2, cindex, skip_dynamic, test_alpha);
	coll_cell_line_iterator<line_intersector> ccli(lint, 0, (skip_dynamic != 0), cobj);
	return ccli.do_iter();
}


bool check_coll_line_exact(point pos1, point pos2, point &cpos, vector3d &cnorm, int &cindex,
						   float splash_val, int ignore_cobj, bool fast, bool test_alpha, bool skip_dynamic) {

	cindex = -1;
	float z_lb(czmin), z_ub(czmax);
	
	if (splash_val > 0.0) { // handle water splashes
		z_lb = min(z_lb, zbottom);
		z_ub = max(z_ub, max(ztop, water_plane_z)); // max of dynamic and static water
	}
	if (!do_line_clip_scene(pos1, pos2, z_lb, z_ub)) return 0;
	line_intersector_exact lint(pos1, pos2, cpos, cnorm, cindex, splash_val, ignore_cobj, test_alpha);
	coll_cell_line_iterator<line_intersector_exact> ccli(lint, 0, skip_dynamic, -1, fast);
	ccli.do_iter();
	lint.finish();
	return (cindex >= 0);
}


bool cobj_contained(point pos1, const point *pts, unsigned npts, int cobj) {

	static int last_cobj(-1);

	if (CACHE_OCCLUDER && last_cobj >= 0 && last_cobj != cobj && !coll_objects[last_cobj].disabled()) {
		if (is_contained(pos1, pts, npts, coll_objects[last_cobj].d)) return 1;
	}
	assert(npts > 0);
	point pos2(pts[0]);
	if (!do_line_clip_scene(pos1, pos2, czmin, czmax)) return 0;
	line_intersector_occlusion_cobjs lint(pos1, pos2, pts, npts);
	coll_cell_line_iterator<line_intersector_occlusion_cobjs> ccli(lint, 1, 1, cobj);
	
	if (ccli.do_iter()) {
		last_cobj = lint.get_coll_cobj();
		assert(last_cobj >= 0);
		return 1;
	}
	return 0;
}


bool get_coll_line_cobjs(point pos1, point pos2, int cobj, vector<int> &cobjs) {

	cobjs.resize(0);
	if (!do_line_clip_scene(pos1, pos2, czmin, czmax)) return 0;
	line_intersector_get_cobjs lint(pos1, pos2, cobjs);
	coll_cell_line_iterator<line_intersector_get_cobjs> ccli(lint, 1, 1, cobj);
	ccli.do_iter();
	return (!cobjs.empty());
}


bool coll_pt_vis_test_large2(point pos1, point pos2, vector<int> &cobjs, int cobj, float radius, int skip_dynamic, bool tl) {

	assert(radius > 0.0);
	cobjs.resize(0);
	if (!do_line_clip_scene(pos1, pos2, max(zbottom, czmin), czmax)) return 0;
	line_intersector_cylinder lint(pos1, pos2, cobjs, cobj, radius, skip_dynamic, tl);
	coll_cell_line_iterator<line_intersector_cylinder> ccli(lint, 1, (skip_dynamic != 0), cobj);
	return ccli.do_iter(radius);
}


bool is_occluded(vector<int> const &occluders, point const *const pts, int npts, point const &camera) {

	unsigned const nocc(occluders.size());

	for (unsigned i = 0; i < nocc; ++i) { // cache last occluder?, promote to the front if occluded?
		coll_obj const &cobj(coll_objects[occluders[i]]);
		if (cobj.status == COLL_STATIC && is_contained(camera, pts, npts, cobj.d)) return 1;
	}
	return 0;
}


void get_occluders() { // 18M total, 380K unique

	RESET_TIME;
	if (!(display_mode & 0x08)) return;
	static unsigned startval(0), stopped_count(0);
	static bool first_run(1);
	unsigned const skipval(first_run ? 0 : 12); // spread update across many frames
	first_run = 0;
	if (++startval >= skipval) startval = 0;
	static point last_camera(FAR_CLIP, FAR_CLIP, FAR_CLIP);
	point const camera(get_camera_pos());
	
	if (p2p_dist(camera, last_camera) < 0.1*HALF_DXY) { // camera hasn't moved much
		if (skipval == 0 || stopped_count == skipval) return;
		++stopped_count;
	}
	else {
		stopped_count = 0;
		last_camera   = camera;
	}
	unsigned const ncobjs(coll_objects.size());

	for (unsigned i = startval; i < ncobjs; i += max(1U, skipval)) {
		coll_obj &cobj(coll_objects[i]);
		if (!cobj.fixed || cobj.status != COLL_STATIC || cobj.cp.surfs == EF_ALL) continue;
		get_coll_line_cobjs(camera, cube_center(cobj.d), i, coll_objects[i].occluders);
	}
	if (skipval <= 1) {PRINT_TIME("Occlusion Preprocessing");}
}



