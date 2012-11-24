// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 5/1/02

#include "3DWorld.h"
#include "mesh.h"
#include "transform_obj.h"
#include "player_state.h"
#include "physics_objects.h"
#include "openal_wrap.h"


bool const REMOVE_ALL_COLL   = 1;
bool const ALWAYS_ADD_TO_HCM = 0;
unsigned const CAMERA_STEPS  = 10;
unsigned const PURGE_THRESH  = 20;
float const CAMERA_MESH_DZ   = 0.1; // max dz on mesh


// Global Variables
bool camera_on_snow(0);
int camera_coll_id(-1);
float czmin(FAR_CLIP), czmax(-FAR_CLIP), coll_rmax(0.0), model_czmin(czmin), model_czmax(czmax);
point camera_last_pos(all_zeros); // not sure about this, need to reset sometimes
coll_obj_group coll_objects;

extern int camera_coll_smooth, game_mode, world_mode, xoff, yoff, camera_change, display_mode, scrolling, animate2;
extern int camera_in_air, mesh_scale_change, camera_invincible, camera_flight, do_run, num_smileys;
extern float TIMESTEP, temperature, zmin, base_gravity, ftick, tstep, zbottom, ztop, fticks;
extern double camera_zh;
extern dwobject def_objects[];
extern obj_type object_types[];
extern player_state *sstates;
extern platform_cont platforms;


void add_coll_point(int i, int j, int index, float zminv, float zmaxv, int add_to_hcm, int is_dynamic, int dhcm);
void free_all_coll_objects();



bool decal_obj::is_on_cobj(int cobj) const {

	if (cobj < 0) return 0;
	assert((unsigned)cobj < coll_objects.size()); // can this fail if the cobj was destroyed? coll_objects only increases in size
	coll_obj const &c(coll_objects[cobj]);
	// spheres and cylinders not supported - decals look bad on rounded objects
	if (c.status != COLL_STATIC || (c.type != COLL_CUBE && c.type != COLL_POLYGON)) return 0;
	point const center(ipos + get_platform_delta());
	if (!sphere_cube_intersect(center, SMALL_NUMBER, c)) return 0;
	if (c.type == COLL_CUBE) return 1;
	float t; // polygon case
	point const p1(center - orient*MIN_POLY_THICK), p2(center + orient*MIN_POLY_THICK);
	return line_poly_intersect(p1, p2, c.points, c.npoints, c.norm, t); // doesn't really work on extruded polygons
}


void decal_obj::check_cobj() {

	if (!status || cid < 0) return; // already disabled, or no bound cobj
	
	if (!is_on_cobj(cid)) { // try to find the cobj this is attached to (likely a split part of the original)
		int const xpos(get_xpos(ipos.x)), ypos(get_ypos(ipos.y));
		
		if (point_outside_mesh(xpos, ypos)) {
			status = 0;
			return;
		}
		vector<int> const &cvals(v_collision_matrix[ypos][xpos].cvals);
		cid = -1;

		for (unsigned i = 0; i < cvals.size(); ++i) {
			if (is_on_cobj(cvals[i])) {
				cid = cvals[i];
				break;
			}
		}
	}
	if (cid < 0) { // not found
		status = 0; // remove it
		return; // no longer on a cobj
	}
}


vector3d decal_obj::get_platform_delta() const {

	if (cid >= 0) {
		assert((unsigned)cid < coll_objects.size());
		int const pid(coll_objects[cid].platform_id);
		
		if (pid >= 0) {
			assert((unsigned)pid < platforms.size());
			return platforms[pid].get_delta();
		}
	}
	return all_zeros;
}


class cobj_manager_t {

	vector<int> index_stack;
	coll_obj_group &cobjs;
	unsigned index_top;

	void extend_index_stack(unsigned start, unsigned end) {
		index_stack.resize(end);

		for (size_t i = start; i < end; ++i) { // initialize
			index_stack[i] = (int)i; // put on the free list
		}
	}

public:
	unsigned cobjs_removed;

	cobj_manager_t(coll_obj_group &cobjs_) : cobjs(cobjs_), index_top(0), cobjs_removed(0) {
		extend_index_stack(0, cobjs.size());
	}

	void reserve_cobjs(size_t size) {
		size_t const old_size(cobjs.size());
		if (old_size >= size) return; // already large enough
		size_t const new_size(max(size, 2*old_size)); // prevent small incremental reserves
		cobjs.resize(new_size);
		extend_index_stack(old_size, new_size);
	}

	int get_next_avail_index() {
		size_t const old_size(cobjs.size());
		assert(index_stack.size() == old_size);
		if (index_top >= old_size) reserve_cobjs(2*index_top + 4); // approx double in size
		int const index(index_stack[index_top]);
		assert(size_t(index) < cobjs.size());
		assert(cobjs[index].status == COLL_UNUSED);
		cobjs[index].status = COLL_PENDING;
		index_stack[index_top++] = -1;
		return index;
	}

	void free_index(int index) {
		assert(cobjs[index].status != COLL_UNUSED);
		cobjs[index].status = COLL_UNUSED;

		if (!cobjs[index].fixed) {
			assert(index_top > 0);
			index_stack[--index_top] = index;
		}
	}

	bool swap_and_set_as_coll_objects(coll_obj_group &new_cobjs) {
		if (!cobjs.empty()) return 0;
		cobjs.swap(new_cobjs);
		unsigned const ncobjs(cobjs.size());
		cobjs.resize(cobjs.capacity()); // use up all available capacity
		extend_index_stack(0, cobjs.size());

		for (unsigned i = 0; i < ncobjs; ++i) {
			coll_obj temp_cobj(cobjs[i]);
			temp_cobj.add_as_fixed_cobj(); // don't need to remove it
			assert(cobjs[i].id == i);
		}
		return 1;
	}
};

cobj_manager_t cobj_manager(coll_objects);


void reserve_coll_objects(unsigned size) {
	cobj_manager.reserve_cobjs(size);
}


bool swap_and_set_as_coll_objects(coll_obj_group &new_cobjs) {
	return cobj_manager.swap_and_set_as_coll_objects(new_cobjs);
}


inline void get_params(int &x1, int &y1, int &x2, int &y2, const float d[3][2]) {

	x1 = max(0, get_xpos(d[0][0]));
	y1 = max(0, get_ypos(d[1][0]));
	x2 = min((MESH_X_SIZE-1), get_xpos(d[0][1]));
	y2 = min((MESH_Y_SIZE-1), get_ypos(d[1][1]));
}


void add_coll_cube_to_matrix(int index, int dhcm) {

	int x1, x2, y1, y2;
	coll_obj &cobj(coll_objects[index]);
	bool const is_dynamic(cobj.status == COLL_DYNAMIC);
	float ds[3][2];
	vector3d delta(zero_vector);

	// we adjust the size of the cube to account for all possible platform locations
	// if delta is not aligned with x/y/z axes then the boundary will be an over approximation, which is inefficient but ok
	if (cobj.platform_id >= 0) {
		assert(cobj.platform_id < (int)platforms.size());
		delta = platforms[cobj.platform_id].get_range();
	}
	for (unsigned j = 0; j < 3; ++j) {
		ds[j][0] = cobj.d[j][0] + min(delta[j], 0.0f);
		ds[j][1] = cobj.d[j][1] + max(delta[j], 0.0f);
	}
	get_params(x1, y1, x2, y2, ds);

	for (int i = y1; i <= y2; ++i) {
		for (int j = x1; j <= x2; ++j) {
			add_coll_point(i, j, index, ds[2][0], ds[2][1], 1, is_dynamic, dhcm);
		}
	}
}


int add_coll_cube(cube_t &cube, cobj_params const &cparams, int platform_id, int dhcm) {

	int const index(cobj_manager.get_next_avail_index());
	coll_obj &cobj(coll_objects[index]);
	cube.normalize();
	cobj.copy_from(cube);
	// cache the center point and radius
	cobj.points[0] = cobj.get_cube_center();
	coll_objects.set_coll_obj_props(index, COLL_CUBE, cobj.get_bsphere_radius(), 0.0, platform_id, cparams);
	add_coll_cube_to_matrix(index, dhcm);
	return index;
}


void add_coll_cylinder_to_matrix(int index, int dhcm) {

	int xx1, xx2, yy1, yy2;
	coll_obj &cobj(coll_objects[index]);
	float zminc(cobj.d[2][0]), zmaxc(cobj.d[2][1]), zmin0(zminc), zmax0(zmaxc);
	point const p1(cobj.points[0]), p2(cobj.points[1]);
	float const x1(p1.x), x2(p2.x), y1(p1.y), y2(p2.y), z1(p1.z), z2(p2.z);
	float const radius(cobj.radius), radius2(cobj.radius2), dr(radius2 - radius), rscale((z1-z2)/fabs(dr));
	float const length(p2p_dist(p1, p2)), dt(HALF_DXY/length), r_off(radius + dt*fabs(dr));
	int const radx(int(ceil(radius*DX_VAL_INV))+1), rady(int(ceil(radius*DY_VAL_INV))+1), rxry(radx*rady);
	int const xpos(get_xpos(x1)), ypos(get_ypos(y1));
	bool const is_dynamic(cobj.status == COLL_DYNAMIC);
	get_params(xx1, yy1, xx2, yy2, cobj.d);

	if (cobj.type == COLL_CYLINDER_ROT) {
		float xylen(sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)));
		float const rmin(min(radius, radius2)), rmax(max(radius, radius2));
		bool const vertical(x1 == x2 && y1 == y2), horizontal(fabs(z1 - z2) < TOLERANCE);
		bool const vert_trunc_cone(z1 != z2 && radius != radius2 && rmax > HALF_DXY);

		for (int i = yy1; i <= yy2; ++i) {
			float const yv(get_yval(i)), v2y(y1 - yv);

			for (int j = xx1; j <= xx2; ++j) {
				float xv(get_xval(j));

				if (vertical) { // vertical
					if (vert_trunc_cone) { // calculate zmin/zmax
						xv -= x1;
						float const rval(min(rmax, (float)sqrt(xv*xv + v2y*v2y)));
						zmaxc = ((rval > rmin) ? max(zmin0, min(zmax0, (rscale*(rval - rmin) + z2))) : zmax0);
					} // else near constant radius, so zminc/zmaxc are correct
				}
				else if (horizontal) {
					zminc = z1 - rmax;
					zmaxc = z1 + rmax;
				}
				else { // diagonal
					// too complex/slow to get this right, so just use the bbox zminc/zmaxc
				}
				int add_to_hcm(0);
				
				if (i >= yy1-1 && i <= yy2+1 && j >= xx1-1 && j <= xx2+1) {
					if (vertical) {
						add_to_hcm = ((i-ypos)*(i-ypos) + (j-xpos)*(j-xpos) <= rxry);
					}
					else {
						float const dist(fabs((x2-x1)*(y1-yv) - (x1-xv)*(y2-y1))/xylen - HALF_DXY);

						if (dist < rmax) {
							if (horizontal) {
								float const t(((x1-x2)*(x1-xv) + (y1-y2)*(y1-yv))/(xylen*xylen)); // location along cylinder axis
								add_to_hcm = (t >= -dt && t <= 1.0+dt && dist < min(rmax, (r_off + t*dr)));
							}
							else { // diagonal
								add_to_hcm = 1; // again, too complex/slow to get this right, so just be conservative
							}
						}
					}
				}
				add_coll_point(i, j, index, zminc, zmaxc, add_to_hcm, is_dynamic, dhcm);
			}
		}
	}
	else { // standard vertical constant-radius cylinder
		assert(cobj.type == COLL_CYLINDER);
		int const crsq(radx*rady);

		for (int i = yy1; i <= yy2; ++i) {
			for (int j = xx1; j <= xx2; ++j) {
				int const distsq((i - ypos)*(i - ypos) + (j - xpos)*(j - xpos));
				if (distsq <= crsq) add_coll_point(i, j, index, z1, z2, (distsq <= rxry), is_dynamic, dhcm);
			}
		}
	}
}


int add_coll_cylinder(float x1, float y1, float z1, float x2, float y2, float z2, float radius, float radius2,
					  cobj_params const &cparams, int platform_id, int dhcm)
{
	int type;
	int const index(cobj_manager.get_next_avail_index());
	coll_obj &cobj(coll_objects[index]);
	radius  = fabs(radius);
	radius2 = fabs(radius2);
	float const rav(max(radius, radius2));
	assert(radius > 0.0 || radius2 > 0.0);
	bool const nonvert(x1 != x2 || y1 != y2 || (fabs(radius - radius2)/max(radius, radius2)) > 0.2);

	if (nonvert) {
		type = COLL_CYLINDER_ROT;
	}
	else { // standard vertical constant-radius cylinder
		if (z2 < z1) swap(z2, z1);
		type = COLL_CYLINDER;
	}
	point *points = cobj.points;
	points[0].x = x1; points[0].y = y1; points[0].z = z1;
	points[1].x = x2; points[1].y = y2; points[1].z = z2;
	
	if (dist_less_than(points[0], points[1], TOLERANCE)) { // no near zero length cylinders
		cout << "pt0 = "; points[0].print(); cout << ", pt1 = "; points[1].print(); cout << endl;
		assert(0);
	}
	coll_objects.set_coll_obj_props(index, type, radius, radius2, platform_id, cparams);
	add_coll_cylinder_to_matrix(index, dhcm);
	return index;
}


void add_coll_sphere_to_matrix(int index, int dhcm) {

	int x1, x2, y1, y2;
	coll_obj &cobj(coll_objects[index]);
	point const pt(cobj.points[0]);
	bool const is_dynamic(cobj.status == COLL_DYNAMIC);
	float const radius(cobj.radius);
	int const radx(int(radius*DX_VAL_INV) + 1), rady(int(radius*DY_VAL_INV) + 1);
	int const xpos(get_xpos(pt.x)), ypos(get_ypos(pt.y));
	get_params(x1, y1, x2, y2, cobj.d);
	int const rxry(radx*rady), crsq(radx*rady);

	for (int i = y1; i <= y2; ++i) {
		for (int j = x1; j <= x2; ++j) {
			int const distsq((i - ypos)*(i - ypos) + (j - xpos)*(j - xpos));

			if (distsq <= crsq) { // nasty offset by HALF_DXY to account for discretization error
				float const dz(sqrt(max(0.0f, (radius*radius - max(0.0f, (distsq*dxdy - HALF_DXY))))));
				add_coll_point(i, j, index, (pt.z-dz), (pt.z+dz), (distsq <= rxry), is_dynamic, dhcm);
			}
		}
	}
}


// doesn't work for ellipses when X != Y
int add_coll_sphere(point const &pt, float radius, cobj_params const &cparams, int platform_id, int dhcm) {

	radius = fabs(radius);
	int const index(cobj_manager.get_next_avail_index());
	coll_objects[index].points[0] = pt;
	coll_objects.set_coll_obj_props(index, COLL_SPHERE, radius, radius, platform_id, cparams);
	add_coll_sphere_to_matrix(index, dhcm);
	return index;
}


void add_coll_polygon_to_matrix(int index, int dhcm) { // coll_obj member function?

	int x1, x2, y1, y2;
	coll_obj &cobj(coll_objects[index]);
	get_params(x1, y1, x2, y2, cobj.d);
	bool const is_dynamic(cobj.status == COLL_DYNAMIC);
	float const zminc(cobj.d[2][0]), zmaxc(cobj.d[2][1]); // thickness has already been added/subtracted

	if (cobj.thickness == 0.0 && (x2-x1) <= 1 && (y2-y1) <=1) { // small polygon
		for (int i = y1; i <= y2; ++i) {
			for (int j = x1; j <= x2; ++j) {
				add_coll_point(i, j, index, zminc, zmaxc, 1, is_dynamic, dhcm);
			}
		}
		return;
	}
	vector3d const norm(cobj.norm);
	float const dval(-dot_product(norm, cobj.points[0]));
	float const thick(0.5*cobj.thickness), dx(0.5*DX_VAL), dy(0.5*DY_VAL);
	float const dzx(norm.z == 0.0 ? 0.0 : DX_VAL*norm.x/norm.z), dzy(norm.z == 0.0 ? 0.0 : DY_VAL*norm.y/norm.z);
	float const delta_z(sqrt(dzx*dzx + dzy*dzy));
	vector<tquad_t> pts;
	if (cobj.thickness > MIN_POLY_THICK) thick_poly_to_sides(cobj.points, cobj.npoints, norm, cobj.thickness, pts);
	cube_t cube;
	cube.d[2][0] = zminc - SMALL_NUMBER;
	cube.d[2][1] = zmaxc + SMALL_NUMBER;

	for (int i = y1; i <= y2; ++i) {
		float const yv(get_yval(i));
		cube.d[1][0] = yv - dy - (i==y1)*thick;
		cube.d[1][1] = yv + dy + (i==y2)*thick;

		for (int j = x1; j <= x2; ++j) {
			float const xv(get_xval(j));
			float z1(zmaxc), z2(zminc);
			cube.d[0][0] = xv - dx - (j==x1)*thick;
			cube.d[0][1] = xv + dx + (j==x2)*thick;

			if (cobj.thickness > MIN_POLY_THICK) { // thick polygon	
				assert(!pts.empty());
				bool inside(0);

				for (unsigned k = 0; k < pts.size(); ++k) {
					point const *const p(pts[k].pts);
					vector3d const pn(get_poly_norm(p));

					if (fabs(pn.z) > 1.0E-3) { // ignore near-vertical polygon edges (for now)
						inside |= get_poly_zminmax(p, pts[k].npts, pn, -dot_product(pn, p[0]), cube, z1, z2);
					}
				}
				if (!inside) continue;
			}
			else if (!get_poly_zminmax(cobj.points, cobj.npoints, norm, dval, cube, z1, z2)) {
				continue;
			}
			// adjust z bounds so that they are for the entire cell x/y bounds, not a single point (conservative)
			z1 = max(zminc, (z1 - delta_z));
			z2 = min(zmaxc, (z2 + delta_z));
			add_coll_point(i, j, index, z1, z2, 1, is_dynamic, dhcm);
		} // for j
	} // for i
}


// must be planar, convex polygon with unique consecutive points
int add_coll_polygon(const point *points, int npoints, cobj_params const &cparams,
	float thickness, point const &xlate, int platform_id, int dhcm)
{
	assert(npoints >= 3 && points != NULL); // too strict?
	assert(npoints <= N_COLL_POLY_PTS);
	int const index(cobj_manager.get_next_avail_index());
	coll_obj &cobj(coll_objects[index]);
	//if (thickness == 0.0) thickness = MIN_POLY_THICK;
	cobj.norm = get_poly_norm(points);

	if (npoints == 4) { // average the norm from both triangles in case they're not coplanar
		point const p2[3] = {points[0], points[2], points[3]};
		vector3d const norm2(get_poly_norm(p2));
		cobj.norm += ((dot_product(cobj.norm, norm2) < 0.0) ? -norm2 : norm2);
		cobj.norm.normalize();
	}
	if (cobj.norm == zero_vector) {
		cout << "degenerate polygon created: points:" << endl;
		for (int i = 0; i < npoints; ++i) {points[i].print(); cout << endl;}
		cout << "valid: " << is_poly_valid(points) << endl;
		cobj.norm = plus_z; // FIXME: this shouldn't be possible, but FP accuracy/errors make this tough to prevent
	}
	for (int i = 0; i < npoints; ++i) {
		cobj.points[i] = points[i] + xlate;
	}
	cobj.npoints   = npoints;
	cobj.thickness = thickness;
	float brad;
	point center; // unused
	polygon_bounding_sphere(points, npoints, thickness, center, brad);
	coll_objects.set_coll_obj_props(index, COLL_POLYGON, brad, 0.0, platform_id, cparams);
	add_coll_polygon_to_matrix(index, dhcm);
	return index;
}


void coll_obj::add_as_fixed_cobj() {

	calc_size();
	fixed = 1;
	id    = add_coll_cobj();
}


int coll_obj::add_coll_cobj() {

	int cid(-1);
	cp.is_dynamic = (status == COLL_DYNAMIC);

	switch (type) {
	case COLL_CUBE:
		cid = add_coll_cube(*this, cp, platform_id);
		break;
	case COLL_SPHERE:
		cid = add_coll_sphere(points[0], radius, cp, platform_id);
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		cid = add_coll_cylinder(points[0], points[1], radius, radius2, cp, platform_id);
		break;
	case COLL_POLYGON:
		cid = add_coll_polygon(points, npoints, cp, thickness, all_zeros, platform_id);
		break;
	default:
		assert(0);
	}
	assert(cid >= 0 && size_t(cid) < coll_objects.size());
	coll_objects[cid].destroy  = destroy;
	coll_objects[cid].fixed    = fixed;
	coll_objects[cid].group_id = group_id;
	coll_objects[cid].v_fall   = v_fall;
	coll_objects[cid].texture_offset = texture_offset;
	return cid;
}


void coll_obj::re_add_coll_cobj(int index, int remove_old, int dhcm) {

	if (!fixed) return;
	assert(index >= 0);
	assert(id == -1 || id == (int)index);
	if (remove_old) remove_coll_object(id, 0); // might already have been removed

	switch (type) {
	case COLL_CUBE:
		add_coll_cube_to_matrix(index, dhcm);
		break;
	case COLL_SPHERE:
		add_coll_sphere_to_matrix(index, dhcm);
		break;
	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		add_coll_cylinder_to_matrix(index, dhcm);
		break;
	case COLL_POLYGON:
		add_coll_polygon_to_matrix(index, dhcm);
		break;
	default:
		assert(0);
	}
	cp.is_dynamic = 0;
	status        = COLL_STATIC;
	counter       = 0;
	id            = index;
}

void coll_cell::clear(bool clear_vectors) {

	if (clear_vectors) {
		if (cvals.capacity() > INIT_CCELL_SIZE) cvals.clear(); else cvals.resize(0);
	}
	zmin = occ_zmin =  FAR_CLIP;
	zmax = occ_zmax = -FAR_CLIP;
}


inline void coll_cell::update_zmm(float zmin_, float zmax_, coll_obj const &cobj) {
		
	assert(zmin_ <= zmax_);

	if (cobj.is_occluder()) {
		occ_zmin = min(zmin_, occ_zmin);
		occ_zmax = max(zmax_, occ_zmax);
	}
	zmin = min(zmin_, zmin);
	zmax = max(zmax_, zmax);
	assert(zmin <= occ_zmin);
	assert(zmax >= occ_zmax);
}


void cobj_stats() {

	unsigned ncv(0), nonempty(0), ncobj(0);
	unsigned const csize((unsigned)coll_objects.size());

	for (int y = 0; y < MESH_Y_SIZE; ++y) {
		for (int x = 0; x < MESH_X_SIZE; ++x) {
			unsigned const sz((unsigned)v_collision_matrix[y][x].cvals.size());
			ncv += sz;
			if (sz > 0) ++nonempty;
		}
	}
	for (unsigned i = 0; i < csize; ++i) {
		if (coll_objects[i].status == COLL_STATIC) ++ncobj;
	}
	if (ncobj > 0) {
		cout << "bins = " << XY_MULT_SIZE << ", ne = " << nonempty << ", cobjs = " << ncobj
			 << ", ent = " << ncv << ", per c = " << ncv/ncobj << ", per bin = " << ncv/XY_MULT_SIZE << endl;
	}
}


void add_coll_point(int i, int j, int index, float zminv, float zmaxv, int add_to_hcm, int is_dynamic, int dhcm) {

	assert(!point_outside_mesh(j, i));
	coll_cell &vcm(v_collision_matrix[i][j]);
	assert((unsigned)index < coll_objects.size());
	vcm.add_entry(index);
	coll_obj const &cobj(coll_objects[index]);
	unsigned const size((unsigned)vcm.cvals.size());

	if (size > 1 && cobj.status == COLL_STATIC && coll_objects[vcm.cvals[size-2]].status == COLL_DYNAMIC) {
		std::rotate(vcm.cvals.begin(), vcm.cvals.begin()+size-1, vcm.cvals.end()); // rotate last point to first point???
	}
	if (is_dynamic) return;

	// update the z values if this cobj is part of a vertically moving platform
	// if it's a cube then it's handled in add_coll_cube_to_matrix()
	if (cobj.type != COLL_CUBE && cobj.platform_id >= 0) {
		assert(cobj.platform_id < (int)platforms.size());
		vector3d const range(platforms[cobj.platform_id].get_range());

		if (range.x == 0.0 && range.y == 0.0) { // vertical platform
			if (range.z > 0.0) {
				zmaxv += range.z; // travels up
			}
			else {
				zminv += range.z; // travels down
			}
		}
	}
	if (dhcm == 0 && add_to_hcm && h_collision_matrix[i][j] < zmaxv && (mesh_height[i][j] + 2.0*object_types[SMILEY].radius) > zminv) {
		h_collision_matrix[i][j] = zmaxv;
	}
	if (add_to_hcm || ALWAYS_ADD_TO_HCM) {
		vcm.update_zmm(zminv, zmaxv, cobj);
		czmin = min(zminv, czmin);
		czmax = max(zmaxv, czmax);
	}
}


int remove_coll_object(int index, bool reset_draw) {

	if (index < 0) return 0;
	assert((size_t)index < coll_objects.size());
	coll_obj &c(coll_objects[index]);

	if (c.status == COLL_UNUSED) {
		assert(REMOVE_ALL_COLL);
		return 0;
	}
	if (c.status == COLL_FREED) return 0;
	coll_objects.remove_index_from_ids(index);
	if (reset_draw) c.cp.draw = 0;
	c.status   = COLL_FREED;
	c.waypt_id = -1; // is this necessary?
	
	if (c.status == COLL_STATIC) {
		//free_index(index); // can't do this here - object's collision id needs to be held until purge
		++cobj_manager.cobjs_removed;
		return 0;
	}
	int x1, y1, x2, y2;
	get_params(x1, y1, x2, y2, c.d);

	for (int i = y1; i <= y2; ++i) {
		for (int j = x1; j <= x2; ++j) {
			vector<int> &cvals(v_collision_matrix[i][j].cvals);
			
			for (unsigned k = 0; k < cvals.size() ; ++k) {
				if (cvals[k] == index) {
					cvals.erase(cvals.begin()+k); // can't change zmin or zmax (I think)
					break; // should only be in here once
				}
			}
		}
	}
	cobj_manager.free_index(index);
	return 1;
}


int remove_reset_coll_obj(int &index) {

	int const retval(remove_coll_object(index));
	index = -1;
	return retval;
}


void purge_coll_freed(bool force) {

	if (!force && cobj_manager.cobjs_removed < PURGE_THRESH) return;
	//RESET_TIME;

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			bool changed(0);
			coll_cell &vcm(v_collision_matrix[i][j]);
			unsigned const size((unsigned)vcm.cvals.size());

			for (unsigned k = 0; k < size && !changed; ++k) {
				if (coll_objects[vcm.cvals[k]].freed_unused()) changed = 1;
			}
			// Note: don't actually have to recalculate zmin/zmax unless a removed object was on the top or bottom of the coll cell
			if (!changed) continue;
			vcm.zmin = vcm.occ_zmin = mesh_height[i][j];
			vcm.zmax = vcm.occ_zmax = zmin;
			vector<int>::const_iterator in(vcm.cvals.begin());
			vector<int>::iterator o(vcm.cvals.begin());

			for (; in != vcm.cvals.end(); ++in) {
				coll_obj &cobj(coll_objects[*in]);

				if (!cobj.freed_unused()) {
					if (cobj.status == COLL_STATIC) vcm.update_zmm(cobj.d[2][0], cobj.d[2][1], cobj);
					*o++ = *in;
				}
			}
			vcm.cvals.erase(o, vcm.cvals.end()); // excess capacity?
			h_collision_matrix[i][j] = vcm.zmax; // need to think about add_to_hcm...
		}
	}
	unsigned const ncobjs((unsigned)coll_objects.size());

	for (unsigned i = 0; i < ncobjs; ++i) {
		if (coll_objects[i].status == COLL_FREED) cobj_manager.free_index(i);
	}
	cobj_manager.cobjs_removed = 0;
	//PRINT_TIME("Purge");
}


void remove_all_coll_obj() {

	camera_coll_id = -1; // camera is special - keeps state

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			v_collision_matrix[i][j].clear(1);
			h_collision_matrix[i][j] = mesh_height[i][j];
		}
	}
	for (unsigned i = 0; i < coll_objects.size(); ++i) {
		if (coll_objects[i].status != COLL_UNUSED) {
			coll_objects.remove_index_from_ids(i);
			cobj_manager.free_index(i);
		}
	}
	free_all_coll_objects();
}


void coll_obj_group::set_coll_obj_props(int index, int type, float radius, float radius2, int platform_id, cobj_params const &cparams) {
	
	coll_obj &cobj(at(index)); // Note: this is the *only* place a new cobj is allocated/created
	cobj.status      = (cparams.is_dynamic ? COLL_DYNAMIC : COLL_STATIC);
	cobj.texture_offset = zero_vector;
	cobj.cp          = cparams;
	cobj.id          = index;
	cobj.radius      = radius;
	cobj.radius2     = radius2;
	cobj.v_fall      = 0.0;
	cobj.type        = type;
	cobj.platform_id = platform_id;
	cobj.group_id    = -1;
	cobj.fixed       = 0;
	cobj.counter     = 0;
	cobj.destroy     = 0;
	cobj.coll_type   = 0;
	cobj.last_coll   = 0;
	cobj.is_billboard= 0;
	cobj.falling     = 0;
	cobj.calc_size();
	cobj.set_npoints();
	cobj.calc_bcube();
	if (cparams.is_dynamic) dynamic_ids.must_insert (index);
	if (cparams.draw      ) drawn_ids.must_insert   (index);
	if (platform_id >= 0  ) platform_ids.must_insert(index);
	if ((type == COLL_CUBE) && cparams.light_atten != 0.0) has_lt_atten = 1;
}


void coll_obj_group::remove_index_from_ids(int index) {

	if (index < 0)  return;
	coll_obj &cobj(at(index));
	if (cobj.fixed) return; // won't actually be freed
	if (cobj.status == COLL_DYNAMIC) coll_objects.dynamic_ids.must_erase (index);
	if (cobj.cp.draw               ) coll_objects.drawn_ids.must_erase   (index);
	if (cobj.platform_id >= 0      ) coll_objects.platform_ids.must_erase(index);
	cobj.cp.draw     = 0;
	cobj.platform_id = -1;
}


// only works for mesh collisions
int collision_detect_large_sphere(point &pos, float radius, unsigned flags) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	float const rdx(radius*DX_VAL_INV), rdy(radius*DY_VAL_INV);
	int const crdx((int)ceil(rdx)), crdy((int)ceil(rdy));
	int const x1(max(0, (xpos - (int)rdx))), x2(min(MESH_X_SIZE-1, (xpos + crdx)));
	int const y1(max(0, (ypos - (int)rdy))), y2(min(MESH_Y_SIZE-1, (ypos + crdy)));
	if (x1 > x2 || y1 > y2) return 0; // not sure if this can ever be false
	int const rsq(crdx*crdy);
	float const rsqf(radius*radius);
	bool const z_up(!(flags & (Z_STOPPED | FLOATING)));
	int coll(0);

	for (int i = y1; i <= y2; ++i) {
		float const yval(get_yval(i));
		int const y_dist((i - ypos)*(i - ypos));
		if (y_dist > rsq) continue; // can never satisfy condition below

		for (int j = x1; j <= x2; ++j) {
			if (y_dist + (j - xpos)*(j - xpos) > rsq) continue;
			point const mpt(get_xval(j), yval, mesh_height[i][j]);
			vector3d const v(pos, mpt);
			float const mag_sq(v.mag_sq());
			if (mag_sq >= rsqf) continue;
			float const old_z(pos.z), mag(sqrt(mag_sq));
			pos = mpt;

			if (mag < TOLERANCE) {
				pos.x += radius; // avoid divide by zero, choose x-direction (arbitrary)
			}
			else {
				pos += v*(radius/mag);
			}
			if (!z_up && pos.z >= old_z) pos.z = old_z;
			coll = 1; // return 1; ?
		}
	}
	assert(!is_nan(pos));
	return coll;
}


int check_legal_move(int x_new, int y_new, float zval, float radius, int &cindex) { // not dynamically updated

	if (point_outside_mesh(x_new, y_new)) return 0; // object out of simulation region
	coll_cell const &cell(v_collision_matrix[y_new][x_new]);
	if (cell.cvals.empty()) return 1;
	float const xval(get_xval(x_new)), yval(get_yval(y_new)), z1(zval - radius), z2(zval + radius);
	point const pval(xval, yval, zval);

	for (int k = (int)cell.cvals.size()-1; k >= 0; --k) { // iterate backwards
		int const index(cell.cvals[k]);
		if (index < 0) continue;
		assert(unsigned(index) < coll_objects.size());
		coll_obj &cobj(coll_objects[index]);
		if (cobj.no_collision()) continue;
		if (cobj.status == COLL_STATIC) {
			if (z1 > cell.zmax || z2 < cell.zmin) return 1; // should be OK here since this is approximate, not quite right with, but not quite right without
		}
		else continue; // smileys collision with dynamic objects can be handled by check_vert_collision()
		if (z1 > cobj.d[2][1] || z2 < cobj.d[2][0]) continue;
		bool coll(0);
		
		switch (cobj.type) {
		case COLL_CUBE:
			coll = ((pval.x + radius) >= cobj.d[0][0] && (pval.x - radius) <= cobj.d[0][1] &&
				    (pval.y + radius) >= cobj.d[1][0] && (pval.y - radius) <= cobj.d[1][1]);
			break;
		case COLL_SPHERE:
			coll = dist_less_than(pval, cobj.points[0], (cobj.radius + radius));
			break;
		case COLL_CYLINDER:
			coll = dist_xy_less_than(pval, cobj.points[0], (cobj.radius + radius));
			break;
		case COLL_CYLINDER_ROT:
			coll = (sphere_int_cylinder_sides(pval, radius, cobj.points[0], cobj.points[1], cobj.radius, cobj.radius2));
			break;
		case COLL_POLYGON: // must be coplanar
			{
				float thick, rdist;

				if (sphere_ext_poly_int_base(cobj.points[0], cobj.norm, pval, radius, cobj.thickness, thick, rdist)) {
					point const pos(pval - cobj.norm*rdist);
					assert(cobj.npoints > 0);
					coll = planar_contour_intersect(cobj.points, cobj.npoints, pos, cobj.norm);
				}
			}
			break;
		}
		if (coll) {
			cindex = index;
			return 0;
		}
	}
	return 1;
}


bool is_point_interior(point const &pos, float radius) { // is query point interior to mesh or cobjs

	if (!is_over_mesh(pos)) return 0; // off scene bounds, outside
	if (is_under_mesh(pos)) return 1; // under mesh, inside
	if (pos.z >= czmax) return 0; // above all cobjs (and mesh), outside
	int inside_count(0);
	float const scene_radius(get_scene_radius());

	for (int z = -1; z <= 1; ++z) {
		for (int y = -1; y <= 1; ++y) {
			for (int x = -1; x <= 1; ++x) {
				if (x == 0 && y == 0 && z == 0) continue;
				vector3d const dir(x, y, z);
				point const pos2(pos + dir*scene_radius);
				point cpos;
				vector3d cnorm;
				int cindex;

				if (check_coll_line_exact(pos, pos2, cpos, cnorm, cindex, 0.0, -1, 0, 0, 1)) {
					assert((unsigned)cindex < coll_objects.size());
					coll_obj const &cobj(coll_objects[cindex]);
					if (radius > 0.0 && cobj.line_intersect(pos, (pos + dir*radius))) return 1; // intersects a short line, so inside/close to a cobj

					if (cobj.type == COLL_CUBE) { // cube
						unsigned const dim(get_max_dim(cnorm)); // cube intersection dim
						bool const dir(cnorm[dim] > 0); // cube intersection dir
						inside_count += ((cobj.cp.surfs & EFLAGS[dim][dir]) ? 1 : -1); // inside if cube side is disabled
					}
					else if (cobj.type == COLL_POLYGON && cobj.thickness <= MIN_POLY_THICK) { // planar/thin polygon
						inside_count += ((dot_product(cnorm, cobj.norm) > 0.0) ? 1 : -1); // inside if hit the polygon back face
					}
					else {
						--inside_count; // sphere, cylinder, cone, or extruded polygon: assume outside
					}
				}
				else {
					--inside_count; // no collision, assume outside
				}
			}
		}
	}
	return (inside_count > 0);
}


// ************ begin vert_coll_detector ************


bool dwobject::proc_stuck(bool static_top_coll) {

	float const friction(object_types[type].friction_factor);
	if (friction < 2.0*STICK_THRESHOLD || friction < rand_uniform(2.0, 3.0)*STICK_THRESHOLD) return 0;
	flags |= (static_top_coll ? ALL_COLL_STOPPED : XYZ_STOPPED); // stuck in coll object
	status = 4;
	return 1;
}


bool vert_coll_detector::safe_norm_div(float rad, float radius, vector3d &norm) {

	if (fabs(rad) < 10.0*TOLERANCE) {
		obj.pos.x += radius; // arbitrary
		norm.assign(1.0, 0.0, 0.0);
		return 0;
	}
	return 1;
}


void vert_coll_detector::check_cobj(int index) {

	coll_obj const &cobj(coll_objects[index]);
	if (cobj.no_collision())                         return; // collisions are disabled for this cobj
	if (type == PROJC    && obj.source  == cobj.id)  return; // can't shoot yourself with a projectile
	if (player           && obj.coll_id == cobj.id)  return; // can't collide with yourself
	if (type == LANDMINE && obj.invalid_coll(cobj))  return;
	if (skip_dynamic && cobj.status == COLL_DYNAMIC) return;
	if (only_drawn   && !cobj.cp.might_be_drawn())   return;
	
	if (type == SMOKE) { // special case to allow smoke to pass through small spheres and cylinders
		switch (cobj.type) { // Note: cubes and polygons could be split into many small pieces, so we don't check their sizes
		case COLL_CYLINDER:
		case COLL_CYLINDER_ROT:
			if (cobj.radius2 < 0.25*object_types[type].radius) return;
		case COLL_SPHERE: // fallthrough from above
			if (cobj.radius  < 0.25*object_types[type].radius) return;
			break;
		}
	}
	if (z1 > cobj.d[2][1] || z2 < cobj.d[2][0]) return;
	if (pos.x < (cobj.d[0][0]-o_radius) || pos.x > (cobj.d[0][1]+o_radius)) return;
	if (pos.y < (cobj.d[1][0]-o_radius) || pos.y > (cobj.d[1][1]+o_radius)) return;
	bool const player_step(player && ((type == CAMERA && camera_change) || (cobj.d[2][1] - z1) <= o_radius*C_STEP_HEIGHT));
	check_cobj_intersect(index, 1, player_step);

	if (type == CAMERA && camera_zh > 0.0) {
		unsigned const nsteps((unsigned)ceil(camera_zh/o_radius));
		float const step_sz(camera_zh/nsteps), pz(pos.z), opz(obj.pos.z), poz(pold.z);

		for (unsigned i = 1; i <= nsteps; ++i) {
			float const step(i*step_sz);
			pos.z     = pz  + step;
			obj.pos.z = opz + step;
			pold.z    = poz + step;
			if (pos.z-o_radius > cobj.d[2][1] || pos.z+o_radius < cobj.d[2][0]) continue;
			check_cobj_intersect(index, 0, 0);
		}
		pos.z     = pz;
		obj.pos.z = opz;
		pold.z    = poz;
	}
}


void vert_coll_detector::check_cobj_intersect(int index, bool enable_cfs, bool player_step) {

	coll_obj const &cobj(coll_objects[index]);
	
	if (cobj.type == COLL_CUBE || cobj.type == COLL_CYLINDER) {
		if (o_radius > 0.9*LARGE_OBJ_RAD && !sphere_cube_intersect(pos, o_radius, cobj)) return;
	}
	vector3d norm(zero_vector), pvel(zero_vector);
	bool coll_top(0), coll_bot(0);
	
	if (cobj.platform_id >= 0) { // calculate platform velocity
		assert(cobj.platform_id < (int)platforms.size());
		pvel = platforms[cobj.platform_id].get_velocity();
	}
	vector3d const mdir(motion_dir - pvel*fticks); // not sure if this helps
	float zmaxc(cobj.d[2][1]), zminc(cobj.d[2][0]);

	switch (cobj.type) { // within bounding box of collision object
	case COLL_CUBE:
		{
			if (!sphere_cube_intersect(pos, o_radius, cobj, (pold - mdir), obj.pos, norm, cdir, 0)) break; // shouldn't get here much when this fails
			coll_top = (cdir == 5);
			coll_bot = (cdir == 4);
			lcoll    = 1;

			if (!coll_top && !coll_bot && player_step) {
				lcoll   = 0; // can step up onto the object
				obj.pos = pos; // reset pos
				norm    = zero_vector;
				break;
			}
			if (coll_top) { // +z collision
				if (cobj.contains_pt_xy(pos)) lcoll = 2;
				float const rdist(max(max(max((pos.x-(cobj.d[0][1]+o_radius)), ((cobj.d[0][0]-o_radius)-pos.x)),
					(pos.y-(cobj.d[1][1]+o_radius))), ((cobj.d[1][0]-o_radius)-pos.y)));
				
				if (rdist > 0.0) {
					obj.pos.z -= o_radius;
					if (o_radius > rdist) obj.pos.z += sqrt(o_radius*o_radius - rdist*rdist);
				}
				break;
			}
		}
		break;

	case COLL_SPHERE:
		{
			float const radius(cobj.radius + o_radius);
			float rad(p2p_dist_sq(pos, cobj.points[0])), reff(radius);

			if (player && cobj.cp.coll_func == landmine_collision) {
				reff += 1.5*object_types[type].radius; // landmine
			}
			if (type == LANDMINE && cobj.is_player()) {
				reff += 1.5*object_types[cobj.type].radius; // landmine
			}
			if (rad <= reff*reff) {
				lcoll = 1;
				rad   = sqrt(rad);
				if (!safe_norm_div(rad, radius, norm)) break;
				norm  = (pos - cobj.points[0])/rad;

				if (rad <= radius) {
					obj.pos = cobj.points[0] + norm*radius;
					assert(!is_nan(obj.pos));
				}
			}
		}
		break;

	case COLL_CYLINDER:
		{
			point const center(cobj.get_center_pt());
			float rad(p2p_dist_xy_sq(pos, center)), radius(cobj.radius); // rad is xy dist

			if (rad <= (radius + o_radius)*(radius + o_radius)) {
				rad    = sqrt(rad);
				lcoll  = 1;
				zmaxc += o_radius;
				zminc -= o_radius;
				float const pozm(pold.z - mdir.z);

				if (!(cobj.cp.surfs & 1) && pozm > (zmaxc - SMALL_NUMBER) && pos.z <= zmaxc) { // collision with top
					if (rad <= radius) lcoll = 2;
					norm.assign(0.0, 0.0, 1.0);
					float const rdist(rad - radius);
					obj.pos.z = zmaxc;
					coll_top  = 1;
					
					if (rdist > 0.0) {
						obj.pos.z -= o_radius;
						if (o_radius >= rdist) obj.pos.z += sqrt(o_radius*o_radius - rdist*rdist);
					}
				}
				else if (!(cobj.cp.surfs & 1) && pozm < (zminc + SMALL_NUMBER) && pos.z >= zminc) { // collision with bottom
					norm.assign(0.0, 0.0, -1.0);
					obj.pos.z = zminc - o_radius;
					coll_bot  = 1;
				}
				else { // collision with side
					if (player_step) {
						norm = plus_z;
						break; // OK, can step up onto cylinder
					}
					radius += o_radius;
					if (!safe_norm_div(rad, radius, norm)) break;
					float const objz(obj.pos.z);
					norm.assign((pos.x - center.x)/rad, (pos.y - center.y)/rad, 0.0);
					for (unsigned d = 0; d < 2; ++d) obj.pos[d] = center[d] + norm[d]*radius;

					/*if (objz > (zmaxc - o_radius) && objz < zmaxc) { // hit on the top edge
						obj.pos.x -= norm.x*o_radius;
						obj.pos.y -= norm.y*o_radius;
						norm.z += (objz - (zmaxc - o_radius))/o_radius; // denominator isn't exactly correct
						norm.normalize();
						obj.pos += norm*o_radius;
					}
					else if (objz < (zminc + o_radius) && objz > zminc) { // hit on the bottom edge
						obj.pos.x -= norm.x*o_radius;
						obj.pos.y -= norm.y*o_radius;
						norm.z -= ((zminc + o_radius) - objz)/o_radius; // denominator isn't exactly correct
						norm.normalize();
						obj.pos += norm*o_radius;
					}*/
				}
			}
		}
		break;

	case COLL_CYLINDER_ROT:
		if (sphere_intersect_cylinder_ipt(pos, o_radius, cobj.points[0], cobj.points[1],
			cobj.radius, cobj.radius2, !(cobj.cp.surfs & 1), obj.pos, norm, 1)) lcoll = 1;
		break;

	case COLL_POLYGON: // must be coplanar
		{
			float thick, rdist, val;
			norm = cobj.norm;
			if (dot_product_ptv(norm, (pold - mdir), cobj.points[0]) < 0.0) norm.negate(); // pos or cobj.points[0]?

			if (sphere_ext_poly_int_base(cobj.points[0], norm, pos, o_radius, cobj.thickness, thick, rdist)) {
				//if (rdist < 0) {rdist = -rdist; norm.negate();}

				if (sphere_poly_intersect(cobj.points, cobj.npoints, pos, norm, rdist, max(0.0f, (thick - MIN_POLY_THICK)))) {
					if (cobj.thickness > MIN_POLY_THICK) { // compute norm based on extruded sides
						vector<tquad_t> pts;
						thick_poly_to_sides(cobj.points, cobj.npoints, cobj.norm, cobj.thickness, pts);
						if (!sphere_intersect_poly_sides(pts, pos, o_radius, val, norm, 1)) break; // no collision
						bool intersects(0), inside(1);
						
						for (unsigned i = 0; i < pts.size(); ++i) { // inefficient and inexact, but closer to correct
							vector3d const norm2(pts[i].get_norm());
							float rdist2(dot_product_ptv(norm2, pos, cobj.points[0]));
							
							if (sphere_poly_intersect(pts[i].pts, pts[i].npts, pos, norm2, rdist2, o_radius)) {
								intersects = 1;
								break;
							}
							if (rdist2 > 0.0) inside = 0;
						}
						if (!intersects && !inside) break; // doesn't intersect any face and isn't completely inside
					}
					else {
						val = 1.01*(thick - rdist); // non-thick polygon
					}
					if (fabs(norm.z) < 0.5 && player_step) { // more horizontal than vertical edge
						norm = zero_vector;
						break; // can step up onto the object
					}
					assert(!is_nan(norm));
					obj.pos += norm*val; // calculate intersection point
					lcoll    = (norm.z > 0.99) ? 2 : 1; // top collision if normal is nearly vertical
				} // end sphere poly int
			} // end sphere int check
		} // end COLL_POLY scope
		break;
	} // switch
	if (!lcoll) return; // no collision
	assert(norm != zero_vector);
	assert(!is_nan(norm));
	bool is_moving(0);
	obj_type const &otype(object_types[type]);
	float const friction(otype.friction_factor);

	// collision with the top of a cube attached to a platform (on first iteration only)
	if (cobj.platform_id >= 0) {
		assert(cobj.platform_id < (int)platforms.size());
		platform const &pf(platforms[cobj.platform_id]);
		is_moving = (lcoll == 2 || friction >= STICK_THRESHOLD);

		if (animate2 && do_coll_funcs && enable_cfs && iter == 0) {
			if (is_moving) { // move with the platform (clip v if large -z?)
				obj.pos += pf.get_last_delta();
			}
			// the coll_top part isn't really right - we want to check for collsion with another object above
			else if ((coll_bot && pf.get_last_delta().z < 0.0) /*|| (coll_top && pf.get_last_delta().z > 0.0)*/) {
				if (player) {
					int const ix((type == CAMERA) ? CAMERA_ID : obj_index);
					smiley_collision(ix, -2, vector3d(0.0, 0.0, -1.0), pos, 2000.0, CRUSHED); // lots of damage
				} // other objects?
			}
		}
		// reset last pos (init_dir) if object is only moving on a platform
		bool const platform_moving(pf.is_moving());
		//if (type == BALL && platform_moving) obj.init_dir = obj.pos;
		if (platform_moving) obj.flags |= PLATFORM_COLL;
	}
	if (animate2 && !player && obj.health <= 0.1) obj.disable();
	vector3d v_old(zero_vector), v0(obj.velocity);
	bool const static_top_coll(lcoll == 2 && cobj.truly_static());

	if (is_moving || friction < STICK_THRESHOLD) {
		v_old = obj.velocity;

		if (otype.elasticity == 0.0 || cobj.cp.elastic == 0.0 || !obj.object_bounce(3, norm, cobj.cp.elastic, pos.z, 0.0, pvel)) {
			if (static_top_coll) {
				obj.flags |= STATIC_COBJ_COLL; // collision with top
				if (otype.flags & OBJ_IS_DROP) obj.velocity = zero_vector;
			}
			if (type != DYNAM_PART && obj.velocity != zero_vector) {
				assert(TIMESTEP > 0.0);
				if (friction > 0.0) obj.velocity *= (1.0 - min(1.0f, (tstep/TIMESTEP)*friction)); // apply kinetic friction
				//for (unsigned i = 0; i < 3; ++i) obj.velocity[i] *= (1.0 - fabs(norm[i])); // norm must be normalized
				orthogonalize_dir(obj.velocity, norm, obj.velocity, 0); // rolling friction model
			}
		}
		else if (already_bounced) {
			obj.velocity = v_old; // can only bounce once
		}
		else {
			already_bounced = 1;
			if (otype.flags & OBJ_IS_CYLIN) obj.init_dir.x += PI*signed_rand_float();
			
			if (type == BALL && cobj.status == COLL_STATIC) { // only static collisions to avoid camera/smiley bounce sounds
				float const vmag(obj.velocity.mag());
				if (vmag > 1.0) gen_sound(SOUND_BOING, obj.pos, min(1.0, 0.1*vmag));
			}
		}
	}
	else { // sticks
		if (cobj.status == COLL_STATIC) {
			if (!obj.proc_stuck(static_top_coll) && static_top_coll) obj.flags |= STATIC_COBJ_COLL; // coll with top
			obj.pos -= norm*(0.1*o_radius); // make sure it still intersects
		}
		obj.velocity = zero_vector; // I think this is correct
	}
	// only use cubes for now, because leaves colliding with tree leaves and branches and resetting the normals is too unstable
	if (cobj.type == COLL_CUBE && (otype.flags & OBJ_IS_FLAT)) obj.set_orient_for_coll(&norm);
		
	if ((otype.flags & OBJ_IS_CYLIN) && !already_bounced) {
		if (fabs(norm.z) == 1.0) { // z collision
			obj.set_orient_for_coll(&norm);
		}
		else { // roll in the direction of the slope with axis along z
			obj.orientation = vector3d(norm.x, norm.y, 0.0).get_norm();
			obj.init_dir.x  = 0.0;
			obj.angle       = 90.0;
		}
	}
	if (do_coll_funcs && enable_cfs && cobj.cp.coll_func != NULL) { // call collision function
		float energy_mult(1.0);
		if (type == PLASMA) energy_mult *= obj.init_dir.x*obj.init_dir.x; // size squared
		float const energy(get_coll_energy(v_old, obj.velocity, otype.mass));
			
		if (!cobj.cp.coll_func(cobj.cp.cf_index, obj_index, v_old, obj.pos, energy_mult*energy, type)) { // invalid collision - reset local collision
			lcoll = 0;
			obj   = temp;
			return;
		}
	}
	if (!(otype.flags & OBJ_IS_DROP) && type != LEAF && type != CHARRED && type != SHRAPNEL &&
		type != BEAM && type != LASER && type != FIRE && type != SMOKE && type != PARTICLE && type != WAYPOINT)
	{
		coll_objects[index].register_coll(TICKS_PER_SECOND, IMPACT);
	}
	obj.verify_data();
		
	if (!obj.disabled() && (otype.flags & EXPL_ON_COLL)) {
		if (cobj.type == COLL_CUBE && cobj.can_be_scorched()) {
			int const dir(cdir >> 1), ds((dir+1)%3), dt((dir+2)%3);
			float const sz(5.0*o_radius);
			float const dmin(min(min((cobj.d[ds][1] - obj.pos[ds]), (obj.pos[ds] - cobj.d[ds][0])),
					                min((cobj.d[dt][1] - obj.pos[dt]), (obj.pos[dt] - cobj.d[dt][0]))));
			if (dmin > sz) gen_decal((obj.pos - norm*o_radius), sz, norm, index, 0.75, BLACK);
		}
		obj.disable();
	}
	if (!obj.disabled() && (fabs(obj.velocity.z) > 1.0 || v0.z > 1.0) && !(obj.flags & STATIC_COBJ_COLL) &&
		((type == BLOOD && (rand()&1) == 0) || (type == CHUNK && !(obj.flags & TYPE_FLAG))))
	{
		gen_decal((obj.pos - norm*o_radius), 2.0*o_radius, norm, index, 1.0, BLOOD_C); // blood/bloody chunk on a non-bottom surface
	}
	deform_obj(obj, norm, v0);
	if (cnorm != NULL) *cnorm = norm;
	obj.flags |= OBJ_COLLIDED;
	coll      |= lcoll; // if not an invalid collision
	lcoll      = 0; // reset local collision
	init_reset_pos(); // reset local state
	if (friction >= STICK_THRESHOLD && (obj.flags & Z_STOPPED)) obj.pos.z = pos.z = z_old;
}


void vert_coll_detector::init_reset_pos() {

	temp = obj; // backup copy
	pos  = obj.pos; // reset local state
	z1   = pos.z - o_radius;
	z2   = pos.z + o_radius;
	if (type == CAMERA) z2 += camera_zh;
}


int vert_coll_detector::check_coll() {

	pold -= obj.velocity*tstep;
	assert(!is_nan(pold));
	assert(type >= 0 && type < NUM_TOT_OBJS);
	o_radius = obj.get_true_radius();
	init_reset_pos();

	if (only_cobj >= 0) {
		assert((unsigned)only_cobj < coll_objects.size());
		check_cobj(only_cobj);
		return coll;
	}
	for (int d = 0; d < 1+!skip_dynamic; ++d) {
		get_coll_sphere_cobjs_tree(obj.pos, o_radius, -1, *this, (d != 0));
	}
	return coll;
}


// ************ end vert_coll_detector ************


// 0 = non vert coll, 1 = X coll, 2 = Y coll, 3 = X + Y coll
int dwobject::check_vert_collision(int obj_index, int do_coll_funcs, int iter, vector3d *cnorm,
	vector3d const &mdir, bool skip_dynamic, bool only_drawn, int only_cobj)
{
	vert_coll_detector vcd(*this, obj_index, do_coll_funcs, iter, cnorm, mdir, skip_dynamic, only_drawn, only_cobj);
	return vcd.check_coll();
}


int dwobject::multistep_coll(point const &last_pos, int obj_index, unsigned nsteps) {

	int any_coll(0);
	vector3d cmove(pos, last_pos);
	float const dist(cmove.mag()); // 0.018

	if (dist < 1.0E-6 || nsteps == 1) {
		any_coll |= check_vert_collision(obj_index, 1, 0); // collision detection
	}
	else {
		float const step(dist/(float)nsteps);
		vector3d const dpos(cmove);
		cmove /= dist;
		pos    = last_pos; // Note: can get stuck in this position if forced off the mesh by a collision

		for (unsigned i = 0; i < nsteps && !disabled(); ++i) {
			point const lpos(pos);
			pos      += cmove*step;
			any_coll |= check_vert_collision(obj_index, (i==nsteps-1), 0, NULL, dpos); // collision detection

			if (type == CAMERA && !camera_change) {
				for (unsigned d = 0; d < 2; ++d) { // x,y
					if (dpos[d]*(pos[d] - lpos[d]) < 0.0) pos[d] = lpos[d]; // negative progress in this dimension, revert
				}
			}
		}
	}
	return any_coll;
}


void add_camera_cobj(point const &pos) {

	camera_coll_id = add_coll_sphere(pos, CAMERA_RADIUS,
		cobj_params(object_types[CAMERA].elasticity, object_types[CAMERA].color, 0, 1, camera_collision, 1));
}


void force_onto_surface_mesh(point &pos) { // for camera

	bool const cflight(game_mode && camera_flight);
	int coll(0);
	float const radius(CAMERA_RADIUS);
	dwobject camera_obj(def_objects[CAMERA]); // make a fresh copy

	if (!cflight) { // make sure camera is over simulation region
		camera_in_air = 0;
		player_clip_to_scene(pos);
	}
	remove_coll_object(camera_coll_id);
	camera_coll_id = -1;
	camera_obj.pos = pos;
	camera_obj.velocity.assign(0.0, 0.0, -1.0);
	//camera_obj.velocity += (pos - camera_last_pos)/tstep; // ???

	if (world_mode == WMODE_GROUND) {
		unsigned const nsteps(CAMERA_STEPS); // *** make the number of steps determined by fticks? ***
		coll  = camera_obj.multistep_coll(camera_last_pos, 0, nsteps);
		pos.x = camera_obj.pos.x;
		pos.y = camera_obj.pos.y;
		if (!cflight) player_clip_to_scene(pos);
	}
	else if (!cflight) {
		pos.z           = int_mesh_zval_pt_off(pos, 1, 0) + radius;
		camera_last_pos = pos;
		camera_change   = 0;
		return; // infinite terrain mode
	}
	if (cflight) {
		if (coll) pos.z = camera_obj.pos.z;
		float const mesh_z(int_mesh_zval_pt_off(pos, 1, 0));
		pos.z = min((camera_last_pos.z + float(C_STEP_HEIGHT*radius)), pos.z); // don't fall and don't rise too quickly
		if (pos.z + radius > zbottom) pos.z = max(pos.z, (mesh_z + radius)); // if not under the mesh
	}
	if (!cflight) {
		if (point_outside_mesh((get_xpos(pos.x) - xoff), (get_ypos(pos.y) - yoff))) {
			pos = camera_last_pos;
			camera_change = 0;
			return;
		}
		set_true_obj_height(pos, camera_last_pos, C_STEP_HEIGHT, sstates[CAMERA_ID].zvel, CAMERA, CAMERA_ID, cflight, camera_on_snow); // status return value is unused?
	}
	camera_on_snow = 0;

	if (display_mode & 0x10) { // walk on snow
		float zval;
		vector3d norm;
		
		if (get_snow_height(pos, radius, zval, norm, 1)) {
			pos.z = zval + radius;
			camera_on_snow = 1;
			camera_in_air  = 0;
		}
	}
	if (camera_coll_smooth) collision_detect_large_sphere(pos, radius, (unsigned char)0);
	if (temperature > W_FREEZE_POINT && is_underwater(pos, 1) && (rand()&1)) gen_bubble(pos);
	
	if (!cflight && camera_change == 0 && camera_last_pos.z != 0.0 && (pos.z - camera_last_pos.z) > CAMERA_MESH_DZ &&
		is_above_mesh(pos) && is_over_mesh(camera_last_pos))
	{
		pos = camera_last_pos; // mesh is too steep for camera to climb
	}
	else {
		camera_last_pos = pos;
	}
	if (camera_change == 1) {
		camera_last_pos = pos;
		camera_change   = 2;
	}
	else {
		camera_change = 0;
	}
	point pos2(pos);
	pos2.z += 0.5*camera_zh;
	add_camera_cobj(pos2);
}


// 0 = no change, 1 = moved up, 2 = falling (unused), 3 = stuck
int set_true_obj_height(point &pos, point const &lpos, float step_height, float &zvel, int type, int id,
	bool flight, bool on_snow, bool skip_dynamic, bool test_only)
{
	int const xpos(get_xpos(pos.x) - xoff), ypos(get_ypos(pos.y) - yoff);
	bool const is_camera(type == CAMERA), is_player(is_camera || (type == SMILEY && id >= 0));
	if (is_player) assert(id >= CAMERA_ID && id < num_smileys);

	if (point_outside_mesh(xpos, ypos)) {
		if (is_player) sstates[id].fall_counter = 0;
		zvel = 0.0;
		return 0;
	}
	float const radius(object_types[type].radius), step(step_height*radius), mh(int_mesh_zval_pt_off(pos, 1, 0)); // *** step height determined by fticks? ***
	pos.z = max(pos.z, (mh + radius));

	if ((display_mode & 0x10) && !test_only) { // walk on snow (smiley and camera, though doesn't actually set smiley z value correctly)
		float zval;
		vector3d norm;
		if (get_snow_height(pos, radius, zval, norm, 1)) pos.z = zval + radius;
	}
	float zmu(mh), z1(pos.z - radius), z2(pos.z + radius);
	if (is_camera /*|| type == WAYPOINT*/) z2 += camera_zh; // add camera height
	coll_cell const &cell(v_collision_matrix[ypos][xpos]);
	int any_coll(0), moved(0);
	float zceil, zfloor;

	for (int k = (int)cell.cvals.size()-1; k >= 0; --k) { // iterate backwards
		int const index(cell.cvals[k]);
		if (index < 0) continue;
		assert(unsigned(index) < coll_objects.size());
		coll_obj const &cobj(coll_objects[index]);
		if (cobj.d[2][0] > z2)         continue; // above the top of the object - can't affect it
		if (any_coll && cobj.d[2][1] < min(z1, min(zmu, zfloor))) continue; // top below a previously seen floor
		if (!cobj.contains_pt_xy(pos)) continue; // test bounding cube
		if (cobj.no_collision())       continue;
		if (skip_dynamic && cobj.status == COLL_DYNAMIC) continue;
		float zt, zb;
		int coll(0);
		
		switch (cobj.type) {
		case COLL_CUBE:
			zt   = cobj.d[2][1];
			zb   = cobj.d[2][0];
			coll = 1;
			break;

		case COLL_CYLINDER:
			if (dist_xy_less_than(pos, cobj.points[0], cobj.radius)) {
				zt   = cobj.d[2][1];
				zb   = cobj.d[2][0];
				coll = 1;
			}
			break;

		case COLL_SPHERE:
			{
				vector3d const vtemp(pos, cobj.points[0]);

				if (vtemp.mag_sq() <= cobj.radius*cobj.radius) {
					float const arg(cobj.radius*cobj.radius - vtemp.x*vtemp.x - vtemp.y*vtemp.y);
					assert(arg >= 0.0);
					float const sqrt_arg(sqrt(arg));
					zt   = cobj.points[0].z + sqrt_arg;
					zb   = cobj.points[0].z - sqrt_arg;
					coll = 1;
				}
			}
			break;

		case COLL_CYLINDER_ROT:
			{
				float t, rad;
				vector3d v1, v2;

				if (sphere_int_cylinder_pretest(pos, radius, cobj.points[0], cobj.points[1], cobj.radius, cobj.radius2, 0, v1, v2, t, rad)) {
					float const rdist(v2.mag());
					
					if (fabs(rdist) > TOLERANCE) {
						float const val(fabs((rad/rdist - 1.0)*v2.z));
						zt   = pos.z + val;
						zb   = pos.z - val;
						coll = 1;
					}
				}
			}
			break;

		case COLL_POLYGON: // must be coplanar, may not work correctly if vertical (but check_vert_collision() should take care of that case)
			{
				coll = 0;
				float const thick(0.5*cobj.thickness);
				bool const poly_z(fabs(cobj.norm.z) > 0.5); // hack to fix bouncy polygons and such - should use a better fix eventually

				if (cobj.thickness > MIN_POLY_THICK && !poly_z) {
					float val;
					vector3d norm;
					vector<tquad_t> pts;
					thick_poly_to_sides(cobj.points, cobj.npoints, cobj.norm, cobj.thickness, pts);
					if (!sphere_intersect_poly_sides(pts, pos, radius, val, norm, 0)) break; // no collision
					float const zminc(cobj.d[2][0]), zmaxc(cobj.d[2][1]);
					zb = zmaxc;
					zt = zminc;

					if (get_poly_zvals(pts, pos.x, pos.y, zb, zt)) {
						zb   = max(zminc, zb);
						zt   = min(zmaxc, zt);
						coll = (zb < zt);
					}
				}
				else if (point_in_polygon_2d(pos.x, pos.y, cobj.points, cobj.npoints, 0, 1)) {
					float const rdist(dot_product_ptv(cobj.norm, pos, cobj.points[0]));
					// works best if the polygon has a face oriented in +z or close
					// note that we don't care if the polygon is intersected in z
					zt   = pos.z + cobj.norm.z*(-rdist + thick);
					zb   = pos.z + cobj.norm.z*(-rdist - thick);
					coll = 1;
					if (zt < zb) swap(zt, zb);
				}
				// clamp to actual polygon bounds (for cases where the object intersects the polygon's plane but not the polygon)
				zt = max(cobj.d[2][0], min(cobj.d[2][1], zt));
				zb = max(cobj.d[2][0], min(cobj.d[2][1], zb));
			}
			break;
		} // end switch
		if (cobj.platform_id >= 0) zt -= 1.0E-6; // subtract a small value so that camera still collides with cobj

		if (coll) {
			if (zt < zb) {
				cout << "type = " << int(cobj.type) << ", zb = " << zb << ", zt = " << zt << ", pos.z = " << pos.z << endl;
				assert(0);
			}
			if (zt <= z1) zmu = max(zmu, zt);
			
			if (any_coll) {
				zceil  = max(zceil,  zt);
				zfloor = min(zfloor, zb);
			}
			else {
				zceil  = zt;
				zfloor = zb;
			}
			if (z2 > zb && z1 < zt) { // overlap: top of object above bottom of surface and bottom of object below top of surface
				if ((zt - z1) <= step) { // step up onto surface
					pos.z = max(pos.z, zt + radius);
					zmu   = max(zmu, zt);
					moved = 1;
				}
				else if (is_camera && camera_change) {
					pos.z = zt + radius;
					zmu   = max(zmu, zt);
				}
				else { // stuck against side of surface
					if (pos.z > zb) { // head inside the object
						if (is_player) sstates[id].fall_counter = 0;
						pos  = lpos; // reset to last known good position
						zvel = 0.0;
						return 3;
					}
					else { // fall down below zb - can recover
						pos.z = zb - radius;
					}
				}
			}
			any_coll = 1;
		} // if coll
	} // for k
	float const g_acc(base_gravity*GRAVITY*tstep*object_types[type].gravity), terminal_v(object_types[type].terminal_vel);
	bool falling(0);

	if (!any_coll || z2 < zfloor) {
		pos.z = mh;
		bool const on_ice(is_camera && (camera_coll_smooth || game_mode) && temperature <= W_FREEZE_POINT && is_underwater(pos));
		if (on_ice) pos.z = water_matrix[ypos][xpos]; // standing on ice
		pos.z += radius;
		if (!on_ice && type != WAYPOINT) modify_grass_at(pos, radius, (type != FIRE), (type == FIRE), 0, 0);
	}
	else {
		zceil = max(zceil, mh);

		if (z1 > zceil) { // bottom of object above all top surfaces
			pos.z = zceil + radius; // object falls to the floor
		}
		else {
			pos.z = zmu + radius;
		}
	}
	if ((is_camera && camera_change) || mesh_scale_change || on_snow) {
		zvel = 0.0;
	}
	else if ((pos.z - lpos.z) < -step) { // falling through the air
		falling = 1;
	}
	else {
		zvel = 0.0;
	}
	if (is_camera) {
		if (falling)        camera_in_air     = 1;
		if (!camera_in_air) camera_invincible = 0;
	}
	if (type == WAYPOINT) {
		// do nothing
	}
	else if (flight) {
		zvel = 0.0;
		if (is_player) sstates[id].fall_counter = 0;
	}
	else if (falling) {
		zvel  = max(-terminal_v, (zvel - g_acc));
		pos.z = max(pos.z, (lpos.z + tstep*zvel));
		
		if (is_player) {
			if (sstates[id].fall_counter == 0) sstates[id].last_dz = 0.0;
			++sstates[id].fall_counter;
			sstates[id].last_dz  += (pos.z - lpos.z);
			sstates[id].last_zvel = zvel;
		}
	}
	else if (is_player) { // falling for several frames continuously and finally stops
		if (sstates[id].fall_counter > 4 && sstates[id].last_dz < 0.0 && sstates[id].last_zvel < 0.0) player_fall(id);
		sstates[id].fall_counter = 0;
	}
	return moved;
}



