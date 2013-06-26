// 3D World - Voxel Code
// by Frank Gennari
// 2/25/12

#include "voxels.h"
#include "marching_cubes.h"
#include "upsurface.h" // for noise_gen_3d
#include "mesh.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "file_utils.h"
#include "openal_wrap.h"


bool const DEBUG_BLOCKS    = 0;
bool const PRE_ALLOC_COBJS = 1;
unsigned const NOISE_TSIZE = 64;

unsigned char const ON_EDGE_BIT    = 0x02;
unsigned char const ANCHORED_BIT   = 0x04;
unsigned char const UNDER_MESH_BIT = 0x08;


voxel_params_t global_voxel_params;
voxel_model_ground terrain_voxel_model;

extern bool group_back_face_cull;
extern int dynamic_mesh_scroll, rand_gen_index, scrolling, display_mode, display_framerate;
extern coll_obj_group coll_objects;


template class voxel_grid<float>; // explicit instantiation


// if size==old_size, we do nothing;  if size==0, we only free the texture
void noise_texture_manager_t::setup(unsigned size, int rseed, float mag, float freq, vector3d const &offset) {

	if (size == tsize) return; // nothing to do
	clear();
	tsize = size;
	if (size == 0) return; // nothing else to do
	voxels.clear();
	voxels.init(tsize, tsize, tsize, vector3d(1,1,1), all_zeros, 0.0, 1);
	voxels.create_procedural(mag, freq, offset, 1, rseed, 654+rand_gen_index);
	noise_tid = voxels.upload_to_3d_texture(GL_MIRRORED_REPEAT);
}


void noise_texture_manager_t::bind_texture(unsigned tu_id) const {

	set_3d_texture_as_current(noise_tid, tu_id);
}


void noise_texture_manager_t::clear() {

	free_texture(noise_tid);
	tsize = 0;
}


float noise_texture_manager_t::eval_at(point const &pos) const {

	unsigned ix;
	if (!voxels.get_ix(pos, ix)) return 0.0; // off the voxel grid
	return voxels[ix];
}


template<typename V> void voxel_grid<V>::init(unsigned nx_, unsigned ny_, unsigned nz_, vector3d const &vsz_,
	point const &center_, V default_val, unsigned num_blocks)
{
	vsz = vsz_;
	assert(vsz.x > 0.0 && vsz.y > 0.0 && vsz.z > 0.0);
	nx = nx_; ny = ny_; nz = nz_;
	xblocks = 1+(nx-1)/num_blocks; // ceil
	yblocks = 1+(ny-1)/num_blocks; // ceil
	unsigned const tot_size(nx * ny * nz);
	assert(tot_size > 0);
	clear();
	resize(tot_size, default_val);
	center = center_;
	lo_pos = center - 0.5*point((nx-1)*vsz.x, (ny-1)*vsz.y, (nz-1)*vsz.z);
}


template<typename V> bool voxel_grid<V>::read(FILE *fp) {

	// FIXME: what about nx, ny, nz, xblocks, yblocks, vsz, center, lo_pos
	assert(fp);
	unsigned sz(0);

	if (fread(&sz, sizeof(unsigned), 1, fp) != 1) {
		cerr << "Error reading voxel_grid size" << endl;
		return 0;
	}
	if (empty()) {
		resize(sz);
	}
	else if (sz != size()) {
		cerr << "Error reading voxel_grid size: expected " << size() << " but got " << sz << endl;
		return 0;
	}
	if (fread(&front(), sizeof(V), size(), fp) != size()) {
		cerr << "Error reading voxel_grid data" << endl;
		return 0;
	}
	return 1;
}


template<typename V> bool voxel_grid<V>::write(FILE *fp) const {

	assert(fp);
	unsigned const sz(size());
	
	if (fwrite(&sz, sizeof(unsigned), 1, fp) != 1) {
		cerr << "Error writing voxel_grid size" << endl;
		return 0;
	}
	if (fwrite(&front(), sizeof(V), size(), fp) != size()) {
		cerr << "Error writing voxel_grid data" << endl;
		return 0;
	}
	return 1;
}


bool voxel_model::from_file(string const &fn) {

	FILE *fp(fopen(fn.c_str(), "rb"));

	if (!fp) {
		cerr << "Error opening voxel file " << fn << " for read" << endl;
		return 0;
	}
	bool const success(read(fp) && outside.read(fp) && ao_lighting.read(fp)); // should ao_lighting be read or recalculated?
	fclose(fp);
	return success;
}


bool voxel_model::to_file(string const &fn) const {

	FILE *fp(fopen(fn.c_str(), "wb"));

	if (!fp) {
		cerr << "Error opening voxel file " << fn << " for write" << endl;
		return 0;
	}
	bool const success(write(fp) && outside.write(fp) && ao_lighting.write(fp)); // should ao_lighting be read or recalculated?
	fclose(fp);
	return success;
}


void voxel_manager::clear() {
	
	outside.clear();
	float_voxel_grid::clear();
}


void voxel_manager::create_procedural(float mag, float freq, vector3d const &offset, bool normalize_to_1, int rseed1, int rseed2) {

	noise_gen_3d ngen;
	ngen.set_rand_seeds(rseed1, rseed2);
	ngen.gen_sines(mag, freq); // create sine table
	unsigned const xyz_num[3] = {nx, ny, nz};
	vector<float> xyz_vals[3];
	ngen.gen_xyz_vals((lo_pos + offset), vsz, xyz_num, xyz_vals); // create xyz values
	float const zscale((params.invert ? -1.0 : 1.0)*params.z_gradient/(nz-1));

	#pragma omp parallel for schedule(static,1)
	for (int y = 0; y < (int)ny; ++y) { // generate voxel values
		for (unsigned x = 0; x < nx; ++x) {
			for (unsigned z = 0; z < nz; ++z) {
#if 1
				float val(ngen.get_val(x, y, z, xyz_vals) + z*zscale);
#else
				point pos(get_pt_at(x, y, z));
				pos += 20.0*fabs(ngen.get_val(0.01*pos))*vector3d(1,1,1); // warp
				float val(ngen.get_val(pos));
#endif
				if (normalize_to_1) val = CLIP_TO_pm1(val);
				set(x, y, z, val); // scale value?
			}
		}
	}
}


void voxel_manager::create_from_cobjs(coll_obj_group &cobjs, float filled_val) {

	for (coll_obj_group::iterator i = cobjs.begin(); i != cobjs.end(); ++i) { // doesn't seem to parallelize well
		if (i->cp.cobj_type != COBJ_TYPE_VOX_TERRAIN) continue; // skip it
		add_cobj_voxels(*i, filled_val);
	}
}


void voxel_manager::add_cobj_voxels(coll_obj &cobj, float filled_val) {

	if (cobj.cp.color.alpha < 0.5) return; // skip transparent objects (should they be here?)
	unsigned const num_test_pts = 4;
	cobj.calc_bcube(); // this is why cobj is passed by non-const reference
	cube_t const &bcube(cobj);
	int llc[3], urc[3];
	get_xyz(bcube.get_llc(), llc);
	get_xyz(bcube.get_urc(), urc);
	float const sphere_radius(0.5*vsz.mag()/num_test_pts); // FIXME: step the radius to generate grayscale intersection values?
	float const dv(1.0/(num_test_pts*num_test_pts*num_test_pts));

	// Note: if we can get the normals from the cobjs here we can use dual contouring for a more exact voxel model
	for (int y = max(0, llc[1]); y <= min(int(ny)-1, urc[1]); ++y) {
		for (int x = max(0, llc[0]); x <= min(int(nx)-1, urc[0]); ++x) {
			for (int z = max(0, llc[2]); z <= min(int(nz)-1, urc[2]); ++z) {
				float val(0.0);

				for (unsigned xx = 0; xx < num_test_pts; ++xx) {
					for (unsigned yy = 0; yy < num_test_pts; ++yy) {
						for (unsigned zz = 0; zz < num_test_pts; ++zz) {
							point const offset(xx*vsz.x/num_test_pts, yy*vsz.y/num_test_pts, zz*vsz.z/num_test_pts);
							if (cobj.sphere_intersects((get_pt_at(x, y, z) + offset), sphere_radius)) {val += dv;} // return value of 2 is considered an intersection
						}
					}
				}
				float &v(get_ref(x, y, z));
				if (val > 0.0) {v = max(v, val*filled_val);}
			}
		}
	}
}


void voxel_manager::atten_at_edges(float val) { // and top (5 edges)

	for (unsigned y = 0; y < ny; ++y) {
		float const vy(1.0 - 2.0*fabs(y - 0.5*ny)/float(ny)); // 0 at edges, 1 at center

		for (unsigned x = 0; x < nx; ++x) {
			float const vx(1.0 - 2.0*fabs(x - 0.5*nx)/float(nx)); // 0 at edges, 1 at center

			for (unsigned z = 0; z < nz; ++z) {
				float const vz(1.0 - 2.0*fabs(z - 0.5*nz)/float(nz)), v(0.25 - vx*vy*vz);
				if (v > 0.0) get_ref(x, y, z) += 8.0*val*v;
			}
		}
	}
}


void voxel_manager::atten_at_top_only(float val) {

	for (unsigned y = 0; y < ny; ++y) {
		for (unsigned x = 0; x < nx; ++x) {
			float top_atten_val(0.0);

			if (params.atten_top_mode == 1) { // atten to mesh
				point const pos(get_pt_at(x, y, 0));
				top_atten_val = interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1);
			}
			else if (params.atten_top_mode == 2) { // atten to random
				point const pos(get_pt_at(x, y, 0));
				top_atten_val = 2.0*eval_mesh_sin_terms(params.height_eval_freq*pos.x, params.height_eval_freq*pos.y);
			}
			for (unsigned z = 0; z < nz; ++z) {
				float &v(get_ref(x, y, z));
	
				if (params.atten_top_mode == 1) { // atten to mesh
					float const z_atten(((get_zv(z)) - top_atten_val)/(vsz.z*nz) - 0.5);
					if (z_atten > 0.0) v += val*z_atten;
				}
				else if (params.atten_top_mode == 2) { // atten to random
					v += top_atten_val + val*(z/float(nz) - 0.5);
				}
				else {
					float const z_atten(z/float(nz) - 0.75);
					if (z_atten > 0.0) v += val*z_atten;
				}
			}
		}
	}
}


void voxel_manager::atten_to_sphere(float val, float inner_radius, bool atten_inner, bool no_atten_zbot) {

	for (unsigned y = 0; y < ny; ++y) {
		float const vy(2.0*fabs(y - 0.5*ny)/float(ny)); // 1 at edges, 0 at center

		for (unsigned x = 0; x < nx; ++x) {
			float const vx(2.0*fabs(x - 0.5*nx)/float(nx)); // 1 at edges, 0 at center

			for (unsigned z = 0; z < nz; ++z) {
				float const deltaz(z - 0.5*nz), zval(no_atten_zbot ? max(0.0f, deltaz) : fabs(deltaz));
				float const vz(2.0*zval/float(nz)); // 1 at edges, 0 at center
				float const radius(sqrt(vx*vx + vy*vy + vz*vz));
				float adj(0.0);
				
				if (radius > inner_radius) {
					adj = (radius - inner_radius)/(1.0 - inner_radius);
				}
				else if (atten_inner) {
					adj = -(inner_radius - radius)/inner_radius;
				}
				get_ref(x, y, z) += val*adj;
			}
		}
	}
}


point voxel_manager::interpolate_pt(float isolevel, point const &pt1, point const &pt2, float const val1, float const val2) const {

	if (abs(isolevel - val1) < TOLERANCE) return pt1;
	if (abs(isolevel - val2) < TOLERANCE) return pt2;
	if (abs(val1     - val2) < TOLERANCE) return pt1;
	float const mu(CLIP_TO_01((isolevel - val1) / (val2 - val1)));
	point pt;
	UNROLL_3X(pt[i_] = pt1[i_] + mu * (pt2[i_] - pt1[i_]););
	return pt;
}


unsigned voxel_manager::add_triangles_for_voxel(tri_data_t::value_type &tri_verts, vertex_map_t<vertex_type_t> &vmap, unsigned x, unsigned y, unsigned z, bool count_only) const {

	float vals[8]; // 8 corner voxel values
	point pts[8]; // corner points
	unsigned cix(0);
	bool all_under_mesh(params.remove_under_mesh && (display_mode & 0x01)); // if mesh draw is enabled

	for (unsigned yhi = 0; yhi < 2; ++yhi) {
		for (unsigned xhi = 0; xhi < 2; ++xhi) {
			unsigned const xx(x+xhi), yy(y+yhi);
			float const xval(get_xv(xx)), yval(get_yv(yy));
			if (all_under_mesh) all_under_mesh = ((outside.get(xx, yy, z) & UNDER_MESH_BIT) != 0);
			
			for (unsigned zhi = 0; zhi < 2; ++zhi) {
				unsigned const zz(z+zhi), vix((xhi^yhi) + 2*yhi + 4*zhi);
				unsigned char const outside_val(outside.get(xx, yy, zz) & 7);
				if (outside_val) cix |= 1 << vix; // outside or on edge
				vals[vix] = (outside_val == ON_EDGE_BIT) ? params.isolevel : get(xx, yy, zz); // check on_edge status
				pts [vix] = point(xval, yval, get_zv(zz));
			}
		}
	}
	assert(cix < 256);
	if (all_under_mesh) return 0;
	unsigned const edge_val(voxel_detail::edge_table[cix]);
	if (edge_val == 0)  return 0; // no polygons
	int const *const tris(voxel_detail::tri_table[cix]);
	unsigned count(0);
	for (unsigned i = 0; tris[i] >= 0; i += 3) {++count;}
	if (count_only) {return count;}
	point vlist[12];

	for (unsigned i = 0; i < 12; ++i) {
		if (!(edge_val & (1 << i))) continue;
		unsigned const *eix = voxel_detail::edge_to_vals[i];
		vlist[i] = interpolate_pt(params.isolevel, pts[eix[0]], pts[eix[1]], vals[eix[0]], vals[eix[1]]);
	}
	for (unsigned i = 0; tris[i] >= 0; i += 3) {
		triangle const tri(vlist[tris[i]], vlist[tris[i+1]], vlist[tris[i+2]]);
		if (tri.pts[0] == tri.pts[1] || tri.pts[0] == tri.pts[2] || tri.pts[1] == tri.pts[2]) continue; // invalid triangle
		tri_verts.add_triangle(tri, vmap);
	}
	return count;
}


void voxel_manager::calc_outside_val(unsigned x, unsigned y, unsigned z, bool is_under_mesh) {

	bool const on_edge(params.make_closed_surface && ((x == 0 || x == nx-1) || (y == 0 || y == ny-1) || (z == 0 || z == nz-1)));
	unsigned char ival(on_edge? ON_EDGE_BIT : (((get(x, y, z) < params.isolevel) ^ params.invert) ? 1 : 0)); // on_edge is considered outside
	if (is_under_mesh) ival |= UNDER_MESH_BIT;
	outside.set(x, y, z, ival);
}


void voxel_manager::determine_voxels_outside() { // determine inside/outside points

	assert(!empty());
	assert(vsz.x > 0.0 && vsz.y > 0.0 && vsz.z > 0.0);
	outside.init(nx, ny, nz, vsz, center, 0, params.num_blocks);
	bool const sphere_mode(params.atten_sphere_mode());

	for (unsigned y = 0; y < ny; ++y) {
		for (unsigned x = 0; x < nx; ++x) {
			point const pos(get_pt_at(x, y, 0));
			int const xpos(get_xpos(pos.x)), ypos(get_xpos(pos.y));
			bool const no_zix(sphere_mode || !use_mesh || point_outside_mesh(xpos, ypos));
			unsigned const zix(no_zix ? 0 : max(0, int((z_min_matrix[ypos][xpos] - lo_pos.z)/vsz.z)));
			
			for (unsigned z = 0; z < nz; ++z) {
				calc_outside_val(x, y, z, (z < zix));
			}
		}
	}
}


void voxel_manager::remove_unconnected_outside() { // check for voxels connected to the mesh surface

	bool const keep_at_edge(params.keep_at_scene_edge == 1 || (params.keep_at_scene_edge == 2 && dynamic_mesh_scroll));
	remove_unconnected_outside_range(keep_at_edge, 0, 0, nx, ny, NULL, NULL);
}


struct block_group_t {
	unsigned v[2][2]; // {x,y} x {lo,hi}

	block_group_t() {v[0][0] = v[0][1] = v[1][0] = v[1][1] = 0;}
	unsigned area() const {return (v[0][1] - v[0][0])*(v[1][1] - v[1][0]);}

	void union_with_group(block_group_t const &g) { // not performance critical
		for (unsigned d = 0; d < 2; ++d) {
			v[d][0] = min(v[d][0], g.v[d][0]);
			v[d][1] = max(v[d][1], g.v[d][1]);
		}
	}
};


void voxel_model::remove_unconnected_outside_modified_blocks() {

	if (modified_blocks.empty()) return;
	int const pad = 1;
	vector<unsigned> to_proc(modified_blocks.begin(), modified_blocks.end());
	vector<unsigned> xy_updated;
	vector<point> updated_pts;
	vector<block_group_t> groups;
	unsigned const indiv_cost((2*pad+1)*(2*pad+1)), num_blocks(params.num_blocks);
	//unsigned group_work(0);

	for (unsigned i = 0; i < to_proc.size(); ++i) {
		int const xbix(to_proc[i] % num_blocks), ybix(to_proc[i] / num_blocks);
		block_group_t group;
		group.v[0][0] = max(0, xbix-pad); // x1
		group.v[1][0] = max(0, ybix-pad); // y1
		group.v[0][1] = min(num_blocks, xbix+pad+1U); // x2
		group.v[1][1] = min(num_blocks, ybix+pad+1U); // y2

		for (unsigned j = i+1; j < to_proc.size(); ++j) {
			int const xbix2(to_proc[j] % num_blocks), ybix2(to_proc[j] / num_blocks);
			block_group_t group2;
			group2.v[0][0] = min(group.v[0][0], (unsigned)max(0, xbix2-pad)); // x1
			group2.v[1][0] = min(group.v[1][0], (unsigned)max(0, ybix2-pad)); // y1
			group2.v[0][1] = max(group.v[0][1], min(num_blocks, xbix2+pad+1U)); // x2
			group2.v[1][1] = max(group.v[1][1], min(num_blocks, ybix2+pad+1U)); // y2

			if (group2.area() <= group.area() + indiv_cost) { // keep the merged group
				group = group2;
				std::swap(to_proc[j], to_proc.back());
				to_proc.pop_back();
				--j;
			}
		}
		groups.push_back(group);
	}
	for (vector<block_group_t>::const_iterator i = groups.begin(); i != groups.end(); ++i) {
		remove_unconnected_outside_range(1, i->v[0][0]*xblocks, i->v[1][0]*yblocks, min(nx, i->v[0][1]*xblocks), min(ny, i->v[1][1]*yblocks), &xy_updated, &updated_pts);
		//group_work += i->area();

		for (vector<unsigned>::const_iterator i = xy_updated.begin(); i != xy_updated.end(); ++i) {
			unsigned const x((*i)%nx), y((*i)/nx);
			assert(x < nx && y < ny);
			unsigned const bx1(max(0, (int)x-1)/xblocks), by1(max(0, (int)y-1)/yblocks);
			unsigned const bx2(min((int)nx-1, (int)x+1)/xblocks), by2(min((int)nx-1, (int)y+1)/yblocks);
		
			for (unsigned by = by1; by <= by2; ++by) {
				for (unsigned bx = bx1; bx <= bx2; ++bx) {
					unsigned const bix(by*num_blocks + bx);
					assert(bix < tri_data.size());
					modified_blocks.insert(bix);
				}
			}
		}
	}
	//cout << "blocks out " << modified_blocks.size() << " groups " << groups.size() << " group work " << group_work << " updated " << updated_pts.size() << " xy_up " << xy_updated.size() << endl;
	if (updated_pts.empty()) return;
	float const fragment_radius(0.5*vsz.mag());
	point center(all_zeros);

	for (vector<point>::const_iterator i = updated_pts.begin(); i != updated_pts.end(); ++i) {
		maybe_create_fragments(*i, fragment_radius, NO_SOURCE, 1);
		center += *i;
	}
	center /= updated_pts.size();
	gen_sound(SOUND_ROCK_FALL, center, CLIP_TO_01(0.05f*updated_pts.size()), 2.0);
}


// outside: 0=inside, 1=outside, 2=on_edge, 4-bit set=anchored, 8-bit set=under mesh
// NOTE: not thread safe due to class member <work>
void voxel_manager::remove_unconnected_outside_range(bool keep_at_edge, unsigned x1, unsigned y1, unsigned x2, unsigned y2,
	vector<unsigned> *xy_updated, vector<point> *updated_pts)
{
	assert(!outside.empty());
	vector<unsigned> &work(temp_work); // stack of voxels to process
	assert(work.empty());
	int const min_range[3] = {x1, y1, 0}, max_range[3] = {x2, y2, nz};

	if (params.atten_sphere_mode() || !use_mesh) { // sphere mode / not mesh mode
		unsigned const x(nx/2), y(ny/2); // add a single point at the center of the sphere (will only work for filled sphere center)

		if (x >= x1 && x <= x2 && y >= y1 && y <= y2) {
			unsigned const ix(outside.get_ix(x, y, nz/2));
			assert(outside[ix] != UNDER_MESH_BIT); // outside or above mesh
			work.push_back(ix); // inside, anchored to the mesh
			outside[ix] |= ANCHORED_BIT; // mark as anchored
		}
	}
	else { // add voxels along the mesh surface
		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				for (unsigned z = 0; z < nz; ++z) {
					unsigned const ix(outside.get_ix(x, y, z));
					if (outside[ix] != UNDER_MESH_BIT) continue; // outside or above mesh
					work.push_back(ix); // inside, anchored to the mesh
					outside[ix] |= ANCHORED_BIT; // mark as anchored
				}
			}
		}
	}

	// add voxels along the scene x/y boundary
	if (keep_at_edge) {
		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (x != x1 && x+1 != x2 && y != y1 && y+1 != y2) continue; // not on scene edge

				for (unsigned z = 0; z < nz; ++z) {
					unsigned const ix(outside.get_ix(x, y, z));
					if (outside[ix] == 1) continue; // outside
					work.push_back(ix); // inside, anchored to the mesh
					outside[ix] |= ANCHORED_BIT; // mark as anchored
				}
			}
		}
	}

	// flood fill anchored regions
	while (!work.empty()) {
		unsigned const cur(work.back());
		work.pop_back();
		assert(cur < outside.size());
		assert(outside[cur] & ANCHORED_BIT);
		unsigned const y(cur/(nz*nx)), cur_xz(cur - y*nz*nx), x(cur_xz/nz), z(cur_xz - x*nz);
		assert(outside.get_ix(x, y, z) == cur);

		for (unsigned dim = 0; dim < 3; ++dim) { // check neighbors
			for (unsigned dir = 0; dir < 2; ++dir) {
				int pos[3] = {x, y, z};
				pos[dim] += dir ? 1 : -1;
				if (pos[dim] < min_range[dim] || pos[dim] >= max_range[dim]) continue; // off the grid
				unsigned const ix(outside.get_ix(pos[0], pos[1], pos[2]));
				assert(ix < outside.size());
						
				if (outside[ix] == 0) { // inside
					work.push_back(ix);
					outside[ix] |= ANCHORED_BIT; // mark as anchored
				}
			}
		}
	} // while

	// if anchored or on_edge remove the anchored bit, else mark outside
	for (unsigned y = y1; y < y2; ++y) {
		for (unsigned x = x1; x < x2; ++x) {
			bool had_update(0);

			for (unsigned z = 0; z < nz; ++z) {
				unsigned const ix(outside.get_ix(x, y, z));

				if (outside[ix] > 1) { // anchored, on edge, or under mesh
					outside[ix] &= ~ANCHORED_BIT; // remove anchored bit
				}
				else if (outside[ix] != 1) { // inside and non-anchored
					if (updated_pts) {updated_pts->push_back(get_pt_at(x, y, z));}
					outside[ix] = 1; // make outside
					had_update  = 1;
					operator[](ix) = params.isolevel - (params.invert ? -TOLERANCE : TOLERANCE); // change voxel value to be outside
				}
			}
			if (had_update && xy_updated) {xy_updated->push_back(y*nx + x);}
		}
	}
}


bool voxel_manager::point_inside_volume(point const &pos) const {

	if (outside.empty()) return 0;
	unsigned ix(0);
	return (outside.get_ix(pos, ix) && !is_outside(ix));
}


bool voxel_manager::point_intersect(point const &center, point *int_pt) const {

	if (!point_inside_volume(center)) return 0;
	if (int_pt) {*int_pt = center;}
	return 1;
}


bool voxel_manager::sphere_intersect(point const &center, float radius, point *int_pt) const {

	if (point_intersect(center, int_pt))  return 1; // optimization
	if (radius == 0.0 || outside.empty()) return 0;
	cube_t bcube;
	bcube.set_from_sphere(center, radius);
	int llc[3], urc[3];
	get_xyz(bcube.get_llc(), llc);
	get_xyz(bcube.get_urc(), urc);
	int const num[3] = {nx, ny, nz};
	UNROLL_3X(llc[i_] = max(0, llc[i_]);)
	UNROLL_3X(urc[i_] = min(num[i_]-1, urc[i_]);)

	for (int y = llc[1]; y <= urc[1]; ++y) {
		for (int x = llc[0]; x <= urc[0]; ++x) {
			unsigned const xy_offset((x + y*nx)*nz);

			for (int z = llc[2]; z <= urc[2]; ++z) {
				if (is_outside(xy_offset+z)) continue;
				point const p(get_pt_at(x, y, z));
				if (!dist_less_than(p, center, radius)) continue;
				if (int_pt) {*int_pt = p;}
				return 1;
			}
		}
	}
	return 0;
}


bool voxel_manager::line_intersect(point const &p1, point const &p2, point *int_pt) const {

	if (outside.empty()) return 0;
	point pa(p1), pb(p2);
	if (!do_line_clip(pa, pb, get_raw_bbox().d)) return 0; // no bbox intersection
	if (point_intersect(pa, int_pt))             return 1; // first point intersects
	if (dist_less_than(pa, pb, TOLERANCE))       return 0; // near zero length line
	float const step0(0.5*min(vsz.x, min(vsz.y, vsz.z))), dist(p2p_dist(pa, pb));
	assert(step0 > 0);
	unsigned const num_steps(ceil(dist/step0));
	assert(num_steps > 0);
	vector3d const delta((pb - pa)/num_steps);
	point p(pa + delta); // first point has already been tested
	
	for (unsigned i = 0; i < num_steps; ++i) {
		if (point_intersect(p, int_pt)) return 1;
		p += delta;
	}
	return 0;
}


unsigned voxel_manager::upload_to_3d_texture(int wrap) const { // only works for float type

	vector<unsigned char> data;
	data.resize(size());

	for (unsigned i = 0; i < size(); ++i) {
		data[i] = (unsigned char)(255*CLIP_TO_01(fabs(operator[](i)))); // use fabs() to convert from [-1,1] to [0,1]
	}
	return create_3d_texture(nx, ny, nz, 1, data, GL_LINEAR, wrap);
}


float voxel_model::get_ao_lighting_val(point const &pos) const {

	if (ao_lighting.empty()) return 1.0;
	unsigned ix(0);
	if (!ao_lighting.get_ix(pos, ix)) return 1.0; // off the voxel grid
	assert(ix < ao_lighting.size());
	return ao_lighting[ix]/255.0;
}


sphere_t voxel_model::get_bsphere() const {

	sphere_t bsphere(center, 0.0);
	
	for (tri_data_t::const_iterator i = tri_data.begin(); i != tri_data.end(); ++i) {
		for (vector<vertex_type_t>::const_iterator v = i->begin(); v != i->end(); ++v) {
			bsphere.radius = max(bsphere.radius, p2p_dist_sq(center, v->v));
		}
	}
	bsphere.radius = sqrt(bsphere.radius);
	return bsphere;
}


bool voxel_model::has_triangles() const {

	for (tri_data_t::const_iterator i = tri_data.begin(); i != tri_data.end(); ++i) {
		if (!i->empty()) return 1;
	}
	return 0;
}


bool voxel_model::has_filled_at_edges() const {

	for (tri_data_t::const_iterator i = tri_data.begin(); i != tri_data.end(); ++i) {
		unsigned const num(i->num_verts());
		assert(!(num % 3));

		for (unsigned j = 0; j < num; j += 3) {
			point v[3];
			UNROLL_3X(v[i_] = i->get_vert(j+i_).v;)

			for (unsigned d = 0; d < 3; ++d) {
				if (v[0][d] == v[1][d] && v[0][d] == v[2][d]) return 1;
			}
		}
	}
	return 0;
}


unsigned voxel_model::get_block_ix(unsigned voxel_ix) const {

	assert(voxel_ix < size());
	unsigned const y(voxel_ix/(nz*nx)), vxz(voxel_ix - y*nz*nx), x(vxz/nz), bx(x/xblocks), by(y/yblocks);
	return by*params.num_blocks + bx;
}


void voxel_model::clear() {

	free_context();
	tri_data.clear();
	pt_to_ix.clear();
	modified_blocks.clear();
	ao_lighting.clear();
	boundary_vnmap.clear();
	voxel_manager::clear();
	volume_added = 0;
}


void voxel_model_ground::clear() {
	
	voxel_model::clear();
	for (unsigned i = 0; i < data_blocks.size(); ++i) {clear_block(i);} // unnecessary?
	data_blocks.clear();
	last_blocks_updated.clear();
}


// returns true if something was cleared
bool voxel_model::clear_block(unsigned block_ix) {

	if (tri_data.empty()) return 0; // tri_data was already cleared
	assert(block_ix < tri_data.size());
	bool const was_nonempty(!tri_data[block_ix].empty());
	tri_data[block_ix].clear();
	return was_nonempty;
}


bool voxel_model_ground::clear_block(unsigned block_ix) {

	bool const ret(voxel_model::clear_block(block_ix));

	if (add_cobjs) {
		assert(block_ix < data_blocks.size());

		for (vector<int>::const_iterator i = data_blocks[block_ix].cids.begin(); i != data_blocks[block_ix].cids.end(); ++i) {
			remove_coll_object(*i);
		}
		data_blocks[block_ix].clear();
	}
	return ret;
}


// returns the number of triangles created
unsigned voxel_model::create_block(unsigned block_ix, bool first_create, bool count_only) {

	assert(block_ix < tri_data.size());
	assert(tri_data[block_ix].empty());
	vertex_map_t<vertex_type_t> vmap(1);
	unsigned const xbix(block_ix%params.num_blocks), ybix(block_ix/params.num_blocks);
	unsigned count(0);
	
	for (unsigned y = ybix*yblocks; y < min(ny-1, (ybix+1)*yblocks); ++y) {
		for (unsigned x = xbix*xblocks; x < min(nx-1, (xbix+1)*xblocks); ++x) {
			for (unsigned z = 0; z < nz-1; ++z) {
				count += add_triangles_for_voxel(tri_data[block_ix], vmap, x, y, z, count_only);
			}
		}
	}
	if (!count_only) {
		if (first_create) { // after the first creation pt_to_ix is out of order
			pt_to_ix[block_ix].pt = (point((xbix+0.5)*xblocks, (ybix+0.5)*yblocks, nz/2)*vsz + lo_pos);
			pt_to_ix[block_ix].ix = block_ix;
		}
		create_block_hook(block_ix);
	}
	return count;
}


void voxel_model_ground::create_block_hook(unsigned block_ix) {

	if (!add_cobjs) return; // nothing to do
	cobj_params cparams(params.elasticity, params.base_color, 0, 0, NULL, 0, params.tids[0]);
	cparams.cobj_type = COBJ_TYPE_VOX_TERRAIN;
	assert(block_ix < data_blocks.size());
	assert(data_blocks[block_ix].cids.empty());
	tri_data_t::value_type const &td(tri_data[block_ix]);
	unsigned const num_verts(td.num_verts());
	assert((num_verts % 3) == 0);
	data_blocks[block_ix].cids.reserve(num_verts/3);

	#pragma omp critical(add_coll_polygon)
	for (unsigned v = 0; v < num_verts; v += 3) {
		point const pts[3] = {td.get_vert(v+0).v, td.get_vert(v+1).v, td.get_vert(v+2).v};
		int const cindex(add_coll_polygon(pts, 3, cparams, 0.0));
		assert(cindex >= 0 && (unsigned)cindex < coll_objects.size());
		if (add_as_fixed) {coll_objects[cindex].fixed = 1;} // mark as fixed so that lmap cells will be generated and cobjs will be re-added
		data_blocks[block_ix].cids.push_back(cindex);
	}
}


void voxel_model::calc_ao_dirs() {

	if (!ao_dirs.empty()) return; // already calculated

	for (int y = -1; y <= 1; ++y) {
		for (int x = -1; x <= 1; ++x) {
			for (int z = -1; z <= 1; ++z) {
				if (x == 0 && y == 0 && z == 0) continue;
				vector3d const delta(x*vsz.x, y*vsz.y, z*vsz.z);
				unsigned const nsteps(max(1, int(params.ao_radius/delta.mag())));
				int const dist(z + (x + y*(int)nx)*(int)nz);
				ao_dirs.push_back(step_dir_t(x, y, z, nsteps, dist));
			}
		}
	}
}


void voxel_model::calc_ao_lighting_for_block(unsigned block_ix, bool increase_only) {

	if (ao_lighting.empty()) return; // nothing to do
	float const norm(params.ao_weight_scale/ao_dirs.size());
	unsigned const xbix(block_ix%params.num_blocks), ybix(block_ix/params.num_blocks);
	unsigned char const end_ray_flags((display_mode & 0x01) ? UNDER_MESH_BIT : 0);
	unsigned const xstep(use_mesh ? max(1U, nx/MESH_X_SIZE ) : 1U);
	unsigned const ystep(use_mesh ? max(1U, ny/MESH_Y_SIZE ) : 1U);
	unsigned const zstep(use_mesh ? max(1U, nz/MESH_SIZE[2]) : 1U);
	unsigned const x_end(min(nx, (xbix+1)*xblocks)), y_end(min(ny, (ybix+1)*yblocks));
	int const voxel_sz[3] = {nx, ny, nz};
	
	#pragma omp parallel for schedule(static,1)
	for (int yi = ybix*yblocks; yi < (int)y_end; yi += ystep) {
		for (unsigned xi = xbix*xblocks; xi < x_end; xi += xstep) {
			if (xi == 0 || yi == 0 || xi >= nx-xstep || (unsigned)yi >= ny-ystep) continue; // at the mesh edges
			bool saw_inside(0);

			for (int zi = nz-2; zi >= 0; zi -= zstep) { // skip top zval
				unsigned const x(min(x_end-1, xi+xstep-1)), y(min(y_end-1, yi+ystep-1)), z(min(nz-1, zi+zstep-1));
				unsigned char const outside_val(outside.get(x, y, z));
				saw_inside |= (outside_val == 0 || (outside_val & end_ray_flags));
				if (!saw_inside) continue;
				if (increase_only && ao_lighting.get(x, y, z) == 255) continue;
				point const pos(ao_lighting.get_pt_at(x, y, z));
				if (use_mesh && !is_over_mesh(pos)) continue;
				float val(0.0);
				
				if (z+1 == nz || !(outside.get(x, y, z+1) & end_ray_flags)) { // above mesh
					for (vector<step_dir_t>::const_iterator i = ao_dirs.begin(); i != ao_dirs.end(); ++i) {
						float cur_val(1.0);
						int cur[3] = {x, y, z};
						// bias to pos side by 1 unit for positive steps to help compensate for grid point vs. grid center alignments
						unsigned max_steps(i->nsteps);
						UNROLL_3X(if (i->dir[i_] > 0) cur[i_] += 1;);
						UNROLL_3X(if (i->dir[i_]) max_steps = min(max_steps, (unsigned)max(0, ((i->dir[i_] < 0) ? cur[i_] : voxel_sz[i_]-cur[i_]-1))););
						int ix(outside.get_ix(cur[0], cur[1], cur[2]));

						for (unsigned s = 0; s < max_steps; ++s) { // take steps in this direction
							ix += i->dist_per_step; // increment first to skip the current voxel
						
							if (outside[ix] == 0 || (outside[ix] & end_ray_flags)) {
								cur_val = float(s)/float(i->nsteps);
								break; // voxel known to be inside the volume or under the mesh
							}
						}
						val += norm*cur_val;
						if (val >= 1.0) break;
					} // for i
					val = CLIP_TO_01(pow(val, params.ao_atten_power));
				}
				unsigned char const ao_val(255.0*val);

				for (unsigned yy = yi; yy < min(y_end, yi+ystep); ++yy) {
					for (unsigned xx = xi; xx < min(x_end, xi+xstep); ++xx) {
						for (unsigned zz = zi; zz < min(nz, zi+zstep); ++zz) {
							ao_lighting.set(xx, yy, zz, ao_val);
						}
					}
				}
			} // for z
		} // for x
	} // for y
}


void voxel_model_space::calc_ao_lighting_for_block(unsigned block_ix, bool increase_only) {

	voxel_model::calc_ao_lighting_for_block(block_ix, increase_only);
	free_ao_and_shadow_texture(); // will be recalculated if needed
}


void voxel_model::calc_ao_lighting() {

	if (empty() || scrolling) return; // too slow for scrolling
	if (params.ao_radius == 0.0 || params.ao_weight_scale == 0.0) return; // no AO lighting
	ao_lighting.init(nx, ny, nz, vsz, center, 255, params.num_blocks);
	calc_ao_dirs();

	for (unsigned block = 0; block < tri_data.size(); ++block) {
		calc_ao_lighting_for_block(block, 0);
	}
}


// returns true if something was updated
bool voxel_model::update_voxel_sphere_region(point const &center, float radius, float val_at_center,
	point *damage_pos, int shooter, unsigned num_fragments)
{
	assert(radius > 0.0);
	if (val_at_center == 0.0 || empty()) return 0;
	bool const material_removed(val_at_center < 0.0);
	if (params.invert) val_at_center *= -1.0; // is this correct?
	unsigned const num[3] = {nx, ny, nz};
	unsigned bounds[3][2]; // {x,y,z} x {lo,hi}
	std::set<unsigned> blocks_to_update;
	float const dist_adjust(0.5*vsz.mag()); // single voxel diagonal half-width
	bool saw_inside(0), saw_outside(0);

	for (unsigned d = 0; d < 3; ++d) {
		bounds[d][0] = max(0, min((int)num[d]-1, int(floor(((center[d] - radius) - lo_pos[d])/vsz[d]))));
		bounds[d][1] = max(0, min((int)num[d]-1, int(ceil (((center[d] + radius) - lo_pos[d])/vsz[d]))));
	}
	for (unsigned y = bounds[1][0]; y <= bounds[1][1]; ++y) {
		for (unsigned x = bounds[0][0]; x <= bounds[0][1]; ++x) {
			bool was_updated(0);

			for (unsigned z = bounds[2][0]; z <= bounds[2][1]; ++z) {
				point const pos(get_pt_at(x, y, z));
				if (use_mesh && is_under_mesh(pos)) continue;
				float const dist(max(0.0f, (p2p_dist(center, pos) - dist_adjust)));
				if (dist >= radius) continue; // too far
				// update voxel values, linear falloff with distance from center (ending at 0.0 at radius)
				float &val(get_ref(x, y, z));
				float const prev_val(val);
				val += val_at_center*(1.0 - dist/radius);
				if (params.normalize_to_1) val = CLIP_TO_pm1(val);
				if (val == prev_val) continue; // no change
				calc_outside_val(x, y, z, ((outside.get(x, y, z) & UNDER_MESH_BIT) != 0));
				was_updated = 1;
				((val      < params.isolevel) ? saw_outside : saw_inside) = 1;
				((prev_val < params.isolevel) ? saw_outside : saw_inside) = 1;
				if (damage_pos) {*damage_pos = pos;}
			}
			if (!was_updated) continue;
			// check adjacent voxels since we will need to update our neighbors at the boundaries
			unsigned const bx1(max((int)x-1, 0        )/xblocks), by1(max((int)y-1, 0        )/yblocks);
			unsigned const bx2(min((int)x+1, (int)nx-1)/xblocks), by2(min((int)y+1, (int)ny-1)/yblocks);

			for (unsigned by = by1; by <= by2; ++by) {
				for (unsigned bx = bx1; bx <= bx2; ++bx) {
					unsigned const block_ix(by*params.num_blocks + bx);
					assert(block_ix < tri_data.size());
					blocks_to_update.insert(block_ix);
				}
			}
		}
	}
	if (!saw_inside || !saw_outside) return 0; // nothing else to do
	std::copy(blocks_to_update.begin(), blocks_to_update.end(), inserter(modified_blocks, modified_blocks.begin()));

	if (material_removed) {
		maybe_create_fragments(center, radius, shooter, num_fragments);
		modify_grass_at(center, max(2.5*radius, vsz.mag()/sqrt(3.0)), 0, 0, 0, 0, 0, 0, 1); // remove any grass at this location
	}
	else {
		volume_added = 1;
	}
	return 1;
}


void voxel_model_ground::maybe_create_fragments(point const &center, float radius, int shooter, unsigned num_fragments) const {

	if (num_fragments == 0) return;
	float blend_val(fabs(eval_noise_texture_at(center)));
	blend_val = min(max(params.tex_mix_saturate*(blend_val - 0.5), -0.5), 0.5) + 0.5;
	colorRGBA color;
	blend_color(color, params.colors[1], params.colors[0], blend_val, 1);
	unsigned const tid(params.tids[blend_val > 0.5]);

	for (unsigned o = 0; o < num_fragments; ++o) {
		vector3d const velocity(signed_rand_vector(0.6));
		point fpos(center + signed_rand_vector_spherical(radius));
		gen_fragment(fpos, velocity, rand_uniform(1.0, 2.0), 0.5*rand_float(), color, tid, 1.0, shooter, 0);
	}
}


unsigned voxel_model::get_texture_at(point const &pos) const {
	return params.tids[fabs(eval_noise_texture_at(pos)) > 0.5];
}


void voxel_model::proc_pending_updates() {

	RESET_TIME;
	if (modified_blocks.empty()) return;
	bool something_removed(0);
	unsigned num_added(0);
	if (params.remove_unconnected == 2) {remove_unconnected_outside_modified_blocks();}
	vector<unsigned> blocks_to_update(modified_blocks.begin(), modified_blocks.end());
	
	// FIXME: can we only remove/add voxels within the modified region of each block?
	//        or, create the block first and only remove triangles that don't exist in the new block + add triangles that don't exist in the old block?
	for (unsigned i = 0; i < blocks_to_update.size(); ++i) {
		something_removed |= clear_block(blocks_to_update[i]);
	}
	if (something_removed) purge_coll_freed(0); // unecessary?

	#pragma omp parallel for schedule(static,1)
	for (int i = 0; i < (int)blocks_to_update.size(); ++i) {
		num_added += (create_block(blocks_to_update[i], 0, 0) > 0);
	}
	if (!(num_added > 0 || something_removed)) { // nothing updated
		modified_blocks.clear();
		return;
	}
	if (!boundary_vnmap.empty()) { // fix block boundary vertex normals
		for (unsigned i = 0; i < blocks_to_update.size(); ++i) {
			update_boundary_normals_for_block(blocks_to_update[i], 0);
		}
	}
	for (unsigned i = 0; i < blocks_to_update.size(); ++i) { // blocks will be sorted by y then x
		calc_ao_lighting_for_block(blocks_to_update[i], !volume_added); // update can only remove, so lighting can only increase
	}
	update_blocks_hook(blocks_to_update, num_added);
	modified_blocks.clear();
	volume_added = 0;
	PRINT_TIME("Process Voxel Updates");
}


void update_ao_texture(block_group_t const &group) {
	assert(group.area() > 0);
	update_smoke_indir_tex_range(group.v[0][0], group.v[0][1], group.v[1][0], group.v[1][1], 1);
}


void voxel_model_ground::update_blocks_hook(vector<unsigned> const &blocks_to_update, unsigned num_added) {

	if (add_cobjs) {
		vector<unsigned> cixs;
		cixs.reserve(num_added);
		bool can_undo_last_cboj_tree_add(!last_blocks_updated.empty());

		for (vector<unsigned>::const_iterator i = last_blocks_updated.begin(); i != last_blocks_updated.end(); ++i) {
			if (modified_blocks.find(*i) == modified_blocks.end()) {can_undo_last_cboj_tree_add = 0; break;}
		}
		for (unsigned i = 0; i < blocks_to_update.size(); ++i) {
			copy(data_blocks[blocks_to_update[i]].cids.begin(), data_blocks[blocks_to_update[i]].cids.end(), back_inserter(cixs));
		}
		if (can_undo_last_cboj_tree_add) {try_undo_last_add_to_cobj_tree(0);}
		last_blocks_updated = blocks_to_update;
		//sort(cixs.begin(), cixs.end()); // doesn't really help
		add_to_cobj_tree(cixs, 0);
	}
	if (!ao_lighting.empty()) {
		block_group_t cur_group;

		for (unsigned i = 0; i < blocks_to_update.size(); ++i) { // blocks will be sorted by y then x
			unsigned const xbix(blocks_to_update[i]%params.num_blocks), x1(xbix*xblocks), x2((xbix+1)*xblocks);
			unsigned const ybix(blocks_to_update[i]/params.num_blocks), y1(ybix*yblocks), y2((ybix+1)*yblocks);
			block_group_t group; // add a border of 1 to account for rounding errors
			group.v[0][0] = max(get_xpos(get_xv(x1))-1, 0); // x1
			group.v[1][0] = max(get_ypos(get_yv(y1))-1, 0); // y1
			group.v[0][1] = min(get_xpos(get_xv(x2))+1, MESH_X_SIZE); // x2
			group.v[1][1] = min(get_ypos(get_yv(y2))+1, MESH_Y_SIZE); // y2

			if (i == 0) {
				cur_group = group;
			}
			else if (group.v[1][0] <= cur_group.v[1][1]) { // y range overlap: new group begins at or before the old one ends
				cur_group.union_with_group(group);
			}
			else {
				update_ao_texture(cur_group);
				cur_group = group;
			}
		}
		if (cur_group.area() > 0) {update_ao_texture(cur_group);}
	}
}


void voxel_model::merge_vn_t::finalize() {

	assert(num > 0);
	if (num == 1) {normal = vn[0]->n; return;}
	normal = zero_vector;
	for (unsigned n = 0; n < num; ++n) {normal += vn[n]->n;} // average the vertex normals
	normal /= num;
	for (unsigned n = 0; n < num; ++n) {vn[n]->n = normal;} // make them equal
}


void voxel_model::update_boundary_normals_for_block(unsigned block_ix, bool calc_average) {

	unsigned const xbix(block_ix%params.num_blocks), ybix(block_ix/params.num_blocks);
	cube_t const bbox(get_xv(xbix*xblocks), get_xv(min(nx-1, (xbix+1)*xblocks)), get_yv(ybix*yblocks), get_yv(min(ny-1, (ybix+1)*yblocks)), 0.0, 0.0);

	for (tri_data_t::value_type::iterator i = tri_data[block_ix].begin(); i != tri_data[block_ix].end(); ++i) {
		if (i->v.x != bbox.d[0][0] && i->v.x != bbox.d[0][1] && i->v.y != bbox.d[1][0] && i->v.y != bbox.d[1][1]) continue; // not at a block boundary
		// Note: should be no duplicates within the same block
		if (calc_average) {boundary_vnmap[i->v].add(*i);} else {boundary_vnmap[i->v].update(*i);}
	}
}


void voxel_model::finalize_boundary_vmap() {

	for (vert_norm_map_t::iterator i = boundary_vnmap.begin(); i != boundary_vnmap.end(); ++i) {
		i->second.finalize();
	}
}


void voxel_model::build(bool verbose, bool do_ao_lighting) {

	RESET_TIME;
	float const atten_thresh((params.invert ? 1.0 : -1.0)*params.atten_thresh);

	switch (params.atten_at_edges) {
	case 0: break; // do nothing
	case 1: atten_at_top_only(atten_thresh); break;
	case 2: atten_at_edges   (atten_thresh); break;
	case 3: atten_to_sphere  (atten_thresh, params.radius_val, 0, 0); break;
	case 4: atten_to_sphere  (atten_thresh, params.radius_val, 1, 0); break;
	case 5: atten_to_sphere  (atten_thresh, params.radius_val, 1, 1); break;
	default: assert(0);
	}
	if (verbose) {PRINT_TIME("  Atten at Top/Edges");}
	determine_voxels_outside();
	if (verbose) {PRINT_TIME("  Determine Voxels Outside");}
	if (params.remove_unconnected) remove_unconnected_outside();
	if (verbose) {PRINT_TIME("  Remove Unconnected");}
	unsigned const tot_blocks(params.num_blocks*params.num_blocks);
	assert(pt_to_ix.empty() && tri_data.empty());
	pt_to_ix.resize(tot_blocks);
	tri_data.resize(tot_blocks, indexed_vntc_vect_t<vertex_type_t>(0));
	pre_build_hook();

	#pragma omp parallel for schedule(static,1)
	for (int block = 0; block < (int)tot_blocks; ++block) {
		create_block(block, 1, 0);
	}
	if (verbose) {PRINT_TIME("  Triangles to Model");}

	if (tot_blocks > 1) { // merge triangle vertices along block seams
		for (unsigned block_ix = 0; block_ix < tot_blocks; ++block_ix) {
			update_boundary_normals_for_block(block_ix, 1);
		}
		finalize_boundary_vmap();
		if (verbose) {PRINT_TIME("  Block Seam Merge");}
	}
	if (do_ao_lighting) {
		calc_ao_lighting();
		if (verbose) {PRINT_TIME("  Voxel AO Lighting");}
	}
}


void voxel_model_ground::pre_build_hook() {

	assert(data_blocks.empty());
	if (!add_cobjs)       return; // nothing to do
	data_blocks.resize(tri_data.size());
	if (!PRE_ALLOC_COBJS) return; // nothing to do
	unsigned num_triangles(0);
	#pragma omp parallel for schedule(static,1)
	for (int block = 0; block < (int)tri_data.size(); ++block) {
		num_triangles += create_block(block, 1, 1);
	}
	if (2*coll_objects.size() < num_triangles) {
		reserve_coll_objects(coll_objects.size() + 1.1*num_triangles); // reserve with 10% buffer
	}
}


void voxel_model_ground::build(bool add_cobjs_, bool add_as_fixed_, bool verbose) {

	add_cobjs    = add_cobjs_;
	add_as_fixed = add_as_fixed_;
	voxel_model::build(verbose);
}


void voxel_model::setup_tex_gen_for_rendering(shader_t &s) {

	assert(s.is_setup());
	noise_tex_gen.setup(NOISE_TSIZE, params.texture_rseed, 1.0, params.noise_freq);
	noise_tex_gen.bind_texture(5); // tu_id = 5
	const char *cnames[2] = {"color0", "color1"};
	unsigned const tu_ids[2] = {0,8};

	for (unsigned i = 0; i < 2; ++i) {
		set_active_texture(tu_ids[i]);
		select_texture(params.tids[i], 0);
		s.add_uniform_color(cnames[i], params.colors[i]);
	}
	set_active_texture(15); // TU_ID=15
	select_texture(params.tids[2], 0); // top texture
	set_active_texture(0);
}


void voxel_model_ground::setup_tex_gen_for_rendering(shader_t &s) {

	s.add_uniform_vector3d("tex_eval_offset", vector3d(DX_VAL*xoff2, DY_VAL*yoff2, 0.0));
	voxel_model::setup_tex_gen_for_rendering(s);
}


void voxel_model_space::calc_shadows(voxel_grid<unsigned char> &shadow_data) const {

	shadow_data.init(nx, ny, nz, vsz, center, 255, params.num_blocks); // default is unshadowed

	for (unsigned y = 0; y < ny; ++y) {
		for (unsigned x = 0; x < nx; ++x) {
			bool shadowed(0);
			
			for (unsigned z = 0; z < nz; ++z) {
				if (z+1 < nz) {shadow_data.set(x, y, z+1, (shadowed ? 0 : 255));} // z offset by 1, before shadowed is updated
				if ((outside.get(x, y, z) & 3) == 0) {shadowed = 1;} // inside
			}
		}
	}
}


void voxel_model_space::extract_shadow_edges(voxel_grid<unsigned char> const &shadow_data) {

	// use the z=nz-1 2D projection
	shadow_edge_tris.clear();
	unsigned const ndiv(max(nx, ny));
	float const step_delta(0.5*min(vsz.x, min(vsz.y, vsz.z)));
	point last(center);

	for (unsigned n = 0; n <= ndiv; ++n) { // first and last ray are the same
		float const angle(TWO_PI*((float)n/(float)ndiv)), dx(cosf(angle)), dy(sinf(angle));
		vector3d const step(step_delta*point(dx, dy, 0.0));
		point pos(center), shadow_edge_pos(all_zeros);
		bool last_unshadowed(0);

		while (1) {
			int i[3]; // x,y,z
			get_xyz(pos, i);
			if (!is_valid_range(i)) break; // off the grid
			bool const unshadowed(shadow_data.get(i[0], i[1], nz-1) != 0);
			if (unshadowed && !last_unshadowed) {shadow_edge_pos = pos;}
			last_unshadowed = unshadowed;
			pos += step;
		}
		if (shadow_edge_pos != all_zeros) {pos = shadow_edge_pos;} // found the edge
		if (n > 0) {shadow_edge_tris.push_back(triangle(pos, last, center));} // last is valid
		last = pos;
	}
}


void voxel_model_space::setup_tex_gen_for_rendering(shader_t &s) {

	voxel_model::setup_tex_gen_for_rendering(s);
	
	if (!ao_lighting.empty()) {
		if (ao_tid == 0) {ao_tid = create_3d_texture(nx, ny, nz, 1, ao_lighting, GL_LINEAR, GL_CLAMP_TO_EDGE);}
		set_3d_texture_as_current(ao_tid, 9);
	}
	if (!shadow_tid) {
		voxel_grid<unsigned char> shadow_data; // 0 == no light/in shadow, 255 = full light/no shadow
		calc_shadows(shadow_data);
		extract_shadow_edges(shadow_data);
		shadow_tid = create_3d_texture(nx, ny, nz, 1, shadow_data, GL_LINEAR, GL_CLAMP_TO_EDGE);
	}
	set_3d_texture_as_current(shadow_tid, 10);
}


void voxel_model::core_render(shader_t &s, bool is_shadow_pass, bool no_vfc) {

	for (vector<pt_ix_t>::const_iterator i = pt_to_ix.begin(); i != pt_to_ix.end(); ++i) {
		if (DEBUG_BLOCKS && s.is_setup()) {
			const char *cnames[2] = {"color0", "color1"};
			for (unsigned d = 0; d < 2; ++d) {
				s.add_uniform_color(cnames[d], ((((i->ix & 1) != 0) ^ ((i->ix & params.num_blocks) != 0)) ? RED : BLUE));
			}
		}
		assert(i->ix < tri_data.size());
		tri_data[i->ix].render(s, is_shadow_pass, GL_TRIANGLES, no_vfc);
	}
}


void voxel_model::render(bool is_shadow_pass) { // not const because of vbo caching, etc.

	if (empty()) return; // nothing to do
	shader_t s;
	set_fill_mode();
	glDisable(GL_LIGHTING); // custom lighting calculations from this point on
	
	if (!is_shadow_pass) {
		set_color_a(BLACK); // ambient will be set by indirect lighting in the shader
		float const min_alpha(0.0); // not needed (yet)
		bool const use_noise_tex(params.tids[0] != params.tids[1] || params.colors[0] != params.colors[1]);
		setup_procedural_shaders(s, min_alpha, 1, 1, 1, use_noise_tex, params.top_tex_used, params.tex_scale, params.noise_scale, params.tex_mix_saturate);
		setup_tex_gen_for_rendering(s);
	}
	BLACK.do_glColor();
	set_color_d(params.base_color);
	float const spec(0.0), shine(1.0);
	set_specular(spec, shine);
	if (group_back_face_cull) glEnable(GL_CULL_FACE);
	sort(pt_to_ix.begin(), pt_to_ix.end(), comp_by_dist(get_camera_pos())); // sort near to far
	core_render(s, is_shadow_pass);
	if (s.is_setup()) {s.end_shader();}
	if (group_back_face_cull) glDisable(GL_CULL_FACE);
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	set_specular(0.0, 1.0);
}


void voxel_model::free_context() {

	tri_data.free_vbos();
	noise_tex_gen.clear();
}


float voxel_model::eval_noise_texture_at(point const &pos) const {
	
	point const spos((pos + vector3d(DX_VAL*xoff2, DY_VAL*yoff2, 0.0))*params.noise_scale);
	point fpos;

	for (unsigned i = 0; i < 3; ++i) {
		int const ival((int)spos[i]);
		fpos[i] = spos[i] - ival; // create fractional part
		if (fpos[i] < 0.0) fpos[i] += 1.0; // make positive
		if ((spos[i] < 0.0) ^ (ival & 1)) fpos[i] = 1.0 - fpos[i]; // apply mirroring
		fpos[i] = fpos[i] - 0.5; // transform from [0,1] to [-0.5,0.5]
	}
	return noise_tex_gen.eval_at(fpos*NOISE_TSIZE);
}


void setup_voxel_landscape(voxel_params_t const &params, float default_val) {

	unsigned const nx((params.xsize > 0) ? params.xsize : MESH_X_SIZE);
	unsigned const ny((params.ysize > 0) ? params.ysize : MESH_Y_SIZE);
	unsigned const nz((params.zsize > 0) ? params.zsize : max((unsigned)MESH_Z_SIZE, (nx+ny)/4));
	//float const zlo(zbottom), zhi(max(ztop, zlo + Z_SCENE_SIZE)); // Note: does not include czmin/czmax range
	float const zlo(zbottom), zhi(min(czmin, zbottom) + Z_SCENE_SIZE); // Note: matches zhi of 3D volume textures/matrices for lighting
	// slightly smaller than 2.0 to avoid z-fighting issues at the edge of the mesh/water
	float const xsz((2.0*(1.0 - 0.05/MESH_X_SIZE)*X_SCENE_SIZE - DX_VAL)/(nx-1));
	float const ysz((2.0*(1.0 - 0.05/MESH_Y_SIZE)*Y_SCENE_SIZE - DY_VAL)/(ny-1));
	vector3d const vsz(xsz, ysz, (zhi - zlo)/nz);
	point const center(-0.5*DX_VAL, -0.5*DY_VAL, 0.5*(zlo + zhi));
	terrain_voxel_model.clear();
	terrain_voxel_model.set_params(params);
	terrain_voxel_model.init(nx, ny, nz, vsz, center, default_val, params.num_blocks);
}


void gen_voxel_landscape() {

	RESET_TIME;
	setup_voxel_landscape(global_voxel_params, 0.0);
	vector3d const gen_offset(DX_VAL*xoff2, DY_VAL*yoff2, 0.0);
	terrain_voxel_model.create_procedural(global_voxel_params.mag, global_voxel_params.freq, gen_offset,
		global_voxel_params.normalize_to_1, global_voxel_params.geom_rseed, 456+rand_gen_index);
	PRINT_TIME(" Voxel Gen");
	terrain_voxel_model.build(global_voxel_params.add_cobjs, 0, 1);
	PRINT_TIME(" Voxels to Triangles/Cobjs");
}


bool gen_voxels_from_cobjs(coll_obj_group &cobjs) {

	RESET_TIME;
	bool has_voxel_cobjs(0);

	for (coll_obj_group::const_iterator i = cobjs.begin(); i != cobjs.end() && !has_voxel_cobjs; ++i) {
		has_voxel_cobjs = (i->cp.cobj_type == COBJ_TYPE_VOX_TERRAIN);
	}
	if (!has_voxel_cobjs) return 0; // nothing to do
	voxel_params_t params(global_voxel_params);
	params.ao_atten_power = 0.7; // user-specified?
	params.ao_radius      = 0.5*(X_SCENE_SIZE + Y_SCENE_SIZE);
	setup_voxel_landscape(params, -1.0);
	terrain_voxel_model.create_from_cobjs(cobjs, 1.0);
	PRINT_TIME(" Cobjs Voxel Gen");
	terrain_voxel_model.build(params.add_cobjs, 1, 1);
	PRINT_TIME(" Cobjs Voxels to Triangles/Cobjs");
	return 1;
}


void gen_voxel_spherical(voxel_model &model, voxel_params_t &params, point const &center, float radius, unsigned size, int rseed) {

	params.remove_unconnected = 2; // always
	params.normalize_to_1 = 0;
	params.atten_at_edges = 4; // sphere (could be 3 or 4)
	params.freq           = 1.2;
	params.mag            = 1.2;
	params.radius_val     = 0.75; // seems to work well
	params.atten_thresh   = 3.0; // user-specified?
	params.tids[0]        = ROCK_TEX;
	params.tids[1]        = MOSSY_ROCK_TEX; // maybe change later
	float const vsz(2.0*radius/size);
	assert(model.empty());
	model.set_params(params);
	model.init(size, size, size, vector3d(vsz, vsz, vsz), center, -1.0, params.num_blocks);
	model.create_procedural(params.mag, params.freq, zero_vector, params.normalize_to_1, params.geom_rseed, rseed);
}


float gen_voxel_rock(voxel_model &model, point const &center, float radius, unsigned size, unsigned num_blocks, int rseed) {

	voxel_params_t params;
	params.num_blocks     = num_blocks; // in each of x and y - subdivision not needed? it produces seams
	params.ao_atten_power = 1.5; // user-specified? seems like a good value for asteroids
	params.ao_radius      = 1.0*radius;

	while (1) { // loop until we get a valid asteroid
		rseed = 27751*rseed + 123; // make unique for each iteration
		model.clear();
		gen_voxel_spherical(model, params, center, radius, size, rseed);
		model.build(0);
		if (model.has_filled_at_edges()) continue; // discard and recreate
		float const gen_radius(model.get_bsphere().radius);
		if (gen_radius > 0.0) return gen_radius; // nonempty
	}
	return 0.0; // never gets here
}


void voxel_file_err(string const &str, int &error) {
	cout << "Error reading voxel config option " << str << "." << endl;
	error = 1;
}


bool parse_voxel_option(FILE *fp) {

	int error(0);
	char strc[MAX_CHARS] = {0};

	if (!read_str(fp, strc)) return 0;
	string const str(strc);

	if (str == "xsize") {
		if (!read_nonzero_uint(fp, global_voxel_params.xsize)) voxel_file_err("xsize", error);
	}
	else if (str == "ysize") {
		if (!read_nonzero_uint(fp, global_voxel_params.ysize)) voxel_file_err("ysize", error);
	}
	else if (str == "zsize") {
		if (!read_nonzero_uint(fp, global_voxel_params.zsize)) voxel_file_err("zsize", error);
	}
	else if (str == "num_blocks") {
		if (!read_nonzero_uint(fp, global_voxel_params.num_blocks)) voxel_file_err("num_blocks", error);
	}
	else if (str == "add_cobjs") {
		if (!read_bool(fp, global_voxel_params.add_cobjs)) voxel_file_err("add_cobjs", error);
	}
	else if (str == "normalize_to_1") {
		if (!read_bool(fp, global_voxel_params.normalize_to_1)) voxel_file_err("normalize_to_1", error);
	}
	else if (str == "mag") {
		if (!read_float(fp, global_voxel_params.mag)) voxel_file_err("mag", error);
	}
	else if (str == "freq") {
		if (!read_float(fp, global_voxel_params.freq)) voxel_file_err("freq", error);
	}
	else if (str == "isolevel") {
		if (!read_float(fp, global_voxel_params.isolevel)) voxel_file_err("isolevel", error);
	}
	else if (str == "elasticity") {
		if (!read_float(fp, global_voxel_params.elasticity) || global_voxel_params.elasticity < 0.0) voxel_file_err("elasticity", error);
	}
	else if (str == "tex_scale") {
		if (!read_float(fp, global_voxel_params.tex_scale)) voxel_file_err("tex_scale", error);
	}
	else if (str == "noise_scale") {
		if (!read_float(fp, global_voxel_params.noise_scale)) voxel_file_err("noise_scale", error);
	}
	else if (str == "noise_freq") {
		if (!read_float(fp, global_voxel_params.noise_freq)) voxel_file_err("noise_freq", error);
	}
	else if (str == "tex_mix_saturate") {
		if (!read_float(fp, global_voxel_params.tex_mix_saturate)) voxel_file_err("tex_mix_saturate", error);
	}
	else if (str == "z_gradient") {
		if (!read_float(fp, global_voxel_params.z_gradient)) voxel_file_err("z_gradient", error);
	}
	else if (str == "radius_val") {
		if (!read_float(fp, global_voxel_params.radius_val) || global_voxel_params.radius_val < 0.0) voxel_file_err("radius_val", error);
	}
	else if (str == "height_eval_freq") {
		if (!read_float(fp, global_voxel_params.height_eval_freq)) voxel_file_err("height_eval_freq", error);
	}
	else if (str == "ao_radius") {
		if (!read_float(fp, global_voxel_params.ao_radius) || global_voxel_params.ao_radius < 0.0) voxel_file_err("ao_radius", error);
	}
	else if (str == "ao_weight_scale") {
		if (!read_float(fp, global_voxel_params.ao_weight_scale) || global_voxel_params.ao_weight_scale <= 0.0) voxel_file_err("ao_weight_scale", error);
	}
	else if (str == "ao_atten_power") {
		if (!read_float(fp, global_voxel_params.ao_atten_power) || global_voxel_params.ao_atten_power <= 0.0) voxel_file_err("ao_atten_power", error);
	}
	else if (str == "invert") {
		if (!read_bool(fp, global_voxel_params.invert)) voxel_file_err("invert", error);
	}
	else if (str == "make_closed_surface") {
		if (!read_bool(fp, global_voxel_params.make_closed_surface)) voxel_file_err("make_closed_surface", error);
	}
	else if (str == "remove_unconnected") {
		if (!read_uint(fp, global_voxel_params.remove_unconnected) || global_voxel_params.remove_unconnected > 2) voxel_file_err("remove_unconnected", error);
	}
	else if (str == "keep_at_scene_edge") {
		if (!read_uint(fp, global_voxel_params.keep_at_scene_edge) || global_voxel_params.keep_at_scene_edge > 2) voxel_file_err("keep_at_scene_edge", error);
	}
	else if (str == "remove_under_mesh") {
		if (!read_bool(fp, global_voxel_params.remove_under_mesh)) voxel_file_err("remove_under_mesh", error);
	}
	else if (str == "atten_top_mode") {
		if (!read_uint(fp, global_voxel_params.atten_top_mode) || global_voxel_params.atten_top_mode > 2) voxel_file_err("atten_top_mode", error);
	}
	else if (str == "atten_at_edges") {
		if (!read_uint(fp, global_voxel_params.atten_at_edges) || global_voxel_params.atten_at_edges > 5) voxel_file_err("atten_at_edges", error);
	}
	else if (str == "atten_thresh") {
		if (!read_float(fp, global_voxel_params.atten_thresh) || global_voxel_params.atten_thresh <= 0.0) voxel_file_err("atten_thresh", error);
	}
	else if (str == "geom_rseed") {
		if (!read_int(fp, global_voxel_params.geom_rseed)) voxel_file_err("geom_rseed", error);
	}
	else if (str == "texture_rseed") {
		if (!read_int(fp, global_voxel_params.texture_rseed)) voxel_file_err("texture_rseed", error);
	}
	else if (str == "tid1") {
		if (!read_str(fp, strc)) voxel_file_err("tid1", error);
		global_voxel_params.tids[0] = get_texture_by_name(std::string(strc));
	}
	else if (str == "tid2") {
		if (!read_str(fp, strc)) voxel_file_err("tid2", error);
		global_voxel_params.tids[1] = get_texture_by_name(std::string(strc));
	}
	else if (str == "tid_top") {
		if (!read_str(fp, strc)) voxel_file_err("tid_top", error);
		global_voxel_params.tids[2] = get_texture_by_name(std::string(strc));
		global_voxel_params.top_tex_used = 1;
	}
	else if (str == "base_color") {
		if (!read_color(fp, global_voxel_params.base_color)) voxel_file_err("base_color", error);
	}
	else if (str == "color1") {
		if (!read_color(fp, global_voxel_params.colors[0])) voxel_file_err("color1", error);
	}
	else if (str == "color2") {
		if (!read_color(fp, global_voxel_params.colors[1])) voxel_file_err("color2", error);
	}
	else {
		cout << "Unrecognized voxel keyword in input file: " << str << endl;
		error = 1;
	}
	return !error;
}


void render_voxel_data(bool shadow_pass) {
	terrain_voxel_model.render(shadow_pass);
}

void free_voxel_context() {
	terrain_voxel_model.free_context();
}

bool point_inside_voxel_terrain(point const &pos) {
	return terrain_voxel_model.point_inside_volume(pos);
}

float get_voxel_terrain_ao_lighting_val(point const &pos) {
	return terrain_voxel_model.get_ao_lighting_val(pos);
}

bool update_voxel_sphere_region(point const &center, float radius, float val_at_center, int shooter, unsigned num_fragments) {

	// optimization/hack to skip the update if the player didn't cause it and the camera can't see it
	if (shooter != CAMERA_ID && !camera_pdu.sphere_visible_test(center, radius)) return 0;
	return terrain_voxel_model.update_voxel_sphere_region(center, radius, val_at_center*(display_framerate ? 1.0 : -1.0), NULL, shooter, num_fragments);
}

void proc_voxel_updates() {
	terrain_voxel_model.proc_pending_updates();
}


