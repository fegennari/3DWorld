// 3D World - Grass Generation and Rendering Code
// by Frank Gennari
// 9/28/10

#include "3DWorld.h"
#include "grass.h"
#include "mesh.h"
#include "physics_objects.h"
#include "textures_3dw.h"
#include "lightmap.h"
#include "gl_ext_arb.h"
#include "shaders.h"


#define MOD_GEOM 0x01


bool grass_enabled(1);
unsigned grass_density(0);
float grass_length(0.02), grass_width(0.002);

extern bool has_dir_lights, has_snow, no_sun_lpos_update;
extern int island, default_ground_tex, read_landscape, display_mode, animate2, frame_counter;
extern unsigned create_voxel_landscape;
extern float vegetation, zmin, zmax, fticks, tfticks, h_sand[], h_dirt[], leaf_color_coherence, tree_deadness, relh_adj_tex;
extern colorRGBA leaf_base_color;
extern vector3d wind;
extern obj_type object_types[];
extern coll_obj_group coll_objects;

bool is_grass_enabled();


void grass_manager_t::grass_t::merge(grass_t const &g) {

	p   = (p + g.p)*0.5; // average locations
	float const dmag1(dir.mag()), dmag2(g.dir.mag());
	dir = (dir/dmag1 + g.dir/dmag2).get_norm() * (0.5*(dmag1 + dmag2)); // average directions and lengths independently
	n   = (n + g.n).get_norm(); // average normals
	//UNROLL_3X(c[i_] = (unsigned char)(unsigned(c[i_]) + unsigned(g.c[i_]))/2;) // don't average colors because they're used for the density filtering hash
	// keep original shadowed bit
	w  += g.w; // add widths to preserve surface area
}

void grass_manager_t::free_vbo() {
	delete_and_zero_vbo(vbo);
}

void grass_manager_t::clear() {
	free_vbo();
	grass.clear();
}

vector3d grass_manager_t::interpolate_mesh_normal(point const &pos) const {

	float const xp((pos.x + X_SCENE_SIZE)*DX_VAL_INV), yp((pos.y + Y_SCENE_SIZE)*DY_VAL_INV);
	int const x0((int)xp), y0((int)yp);
	if (point_outside_mesh(x0, y0) || point_outside_mesh(x0+1, y0+1)) return plus_z; // shouldn't get here
	float const xpi(fabs(xp - (float)x0)), ypi(fabs(yp - (float)y0)); // cubic_interpolate()?
	return vertex_normals[y0+0][x0+0]*((1.0 - xpi)*(1.0 - ypi))
			+ vertex_normals[y0+1][x0+0]*((1.0 - xpi)*ypi)
			+ vertex_normals[y0+0][x0+1]*(xpi*(1.0 - ypi))
		    + vertex_normals[y0+1][x0+1]*(xpi*ypi);
}

void grass_manager_t::add_grass_blade(point const &pos, float cscale, bool on_mesh) {

	vector3d const base_dir(on_mesh ? 0.5*(plus_z + interpolate_mesh_normal(pos)) : plus_z); // average mesh normal and +z for grass on mesh
	vector3d const dir((base_dir + rgen.signed_rand_vector(0.3)).get_norm());
	vector3d const norm(cross_product(dir, rgen.signed_rand_vector()).get_norm());
	float const ilch(1.0 - leaf_color_coherence), dead_scale(CLIP_TO_01(tree_deadness));
	float const base_color[3] = {0.25, 0.6, 0.08};
	float const mod_color [3] = {0.2,  0.2, 0.08};
	float const lbc_mult  [3] = {0.2,  0.4, 0.0 };
	float const dead_color[3] = {0.75, 0.6, 0.0 };
	unsigned char color[3];

	for (unsigned i = 0; i < 3; ++i) {
		float const ccomp(CLIP_TO_01(cscale*(base_color[i] + lbc_mult[i]*leaf_base_color[i] + ilch*mod_color[i]*rgen.rand_float())));
		color[i] = (unsigned char)(255.0*(dead_scale*dead_color[i] + (1.0 - dead_scale)*ccomp));
	}
	float const length(grass_length*rgen.rand_uniform(0.7, 1.3));
	float const width( grass_width *rgen.rand_uniform(0.7, 1.3));
	grass.push_back(grass_t(pos, dir*length, norm, color, width, on_mesh));
}

void grass_manager_t::create_new_vbo() {

	delete_vbo(vbo);
	vbo        = create_vbo();
	data_valid = 0;
}

void grass_manager_t::add_to_vbo_data(grass_t const &g, vector<grass_data_t> &data, unsigned &ix, vector3d &norm) const {

	point const p1(g.p), p2(p1 + g.dir + point(0.0, 0.0, 0.05*grass_length));
	vector3d const binorm(cross_product(g.dir, g.n).get_norm());
	vector3d const delta(binorm*(0.5*g.w));
	norm *= (g.shadowed ? 0.001 : 1.0);
	data[ix++].assign(p1-delta, norm, g.c);
	data[ix++].assign(p1+delta, norm, g.c);
	data[ix++].assign(p2,       norm, g.c);
	assert(ix <= data.size());
}

void grass_manager_t::begin_draw(float spec_weight) const {

	assert(vbo > 0);
	bind_vbo(vbo);
	grass_data_t::set_vbo_arrays();
	select_texture(GRASS_BLADE_TEX);
	set_specular(spec_weight, 20.0);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.99);
	glEnable(GL_COLOR_MATERIAL);
	glDisable(GL_NORMALIZE);
}

void grass_manager_t::end_draw() const {

	glDisable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	set_specular(0.0, 1.0);
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_TEXTURE_2D);
	bind_vbo(0);
	check_gl_error(40);
}


void grass_tile_manager_t::gen_block(unsigned bix) {

	for (unsigned y = 0; y < GRASS_BLOCK_SZ; ++y) {
		for (unsigned x = 0; x < GRASS_BLOCK_SZ; ++x) {
			float const xval(x*DX_VAL), yval(y*DY_VAL);

			for (unsigned n = 0; n < grass_density; ++n) {
				add_grass_blade(point(rgen.rand_uniform(xval, xval+DX_VAL), rgen.rand_uniform(yval, yval+DY_VAL), 0.0), TT_GRASS_COLOR_SCALE, 0); // no mesh normal
			}
		}
	}
	assert(bix+1 < vbo_offsets[0].size());
	vbo_offsets[0][bix+1] = grass.size(); // end of block/beginning of next block
}


void grass_tile_manager_t::gen_lod_block(unsigned bix, unsigned lod) {

	if (lod == 0) {
		gen_block(bix);
		return;
	}
	assert(bix+1 < vbo_offsets[lod].size());
	unsigned const search_dist(1*grass_density); // enough for one cell
	unsigned const start_ix(vbo_offsets[lod-1][bix]), end_ix(vbo_offsets[lod-1][bix+1]); // from previous LOD
	float const dmax(2.5*grass_width*(1 << lod));
	vector<unsigned char> used((end_ix - start_ix), 0); // initially all unused
	
	for (unsigned i = start_ix; i < end_ix; ++i) {
		if (used[i-start_ix]) continue; // already used
		grass.push_back(grass[i]); // seed with an existing grass blade
		float dmin_sq(dmax*dmax); // start at max allowed dist
		unsigned merge_ix(i); // start at ourself (invalid)
		unsigned const end_val(min(i+search_dist, end_ix));

		for (unsigned cur = i+1; cur < end_val; ++cur) {
			float const dist_sq(p2p_dist_sq(grass[i].p, grass[cur].p));
					
			if (dist_sq < dmin_sq) {
				dmin_sq  = dist_sq;
				merge_ix = cur;
			}
		}
		if (merge_ix > i) {
			assert(merge_ix < grass.size());
			assert(merge_ix-start_ix < used.size());
			grass.back().merge(grass[merge_ix]);
			used[merge_ix-start_ix] = 1;
		}
	} // for i
	//cout << "level " << lod << " num: " << (grass.size() - vbo_offsets[lod][bix]) << endl;
	vbo_offsets[lod][bix+1] = grass.size(); // end of current LOD block/beginning of next LOD block
}


void grass_tile_manager_t::clear() {

	grass_manager_t::clear();
	for (unsigned lod = 0; lod < NUM_GRASS_LODS; ++lod) {vbo_offsets[lod].clear();}
}


void grass_tile_manager_t::upload_data() {

	if (empty()) return;
	RESET_TIME;
	vector<grass_data_t> data(3*size()); // 3 vertices per grass blade

	for (unsigned i = 0, ix = 0; i < size(); ++i) {
		vector3d norm(plus_z); // use grass normal? 2-sided lighting? generate normals in vertex shader?
		//vector3d norm(grass[i].n);
		add_to_vbo_data(grass[i], data, ix, norm);
	}
	upload_to_vbo(vbo, data, 0, 1);
	data_valid = 1;
	PRINT_TIME("Grass Tile Upload");
}


void grass_tile_manager_t::gen_grass() {

	RESET_TIME;
	assert(NUM_GRASS_LODS > 0);
	assert((MESH_X_SIZE % GRASS_BLOCK_SZ) == 0 && (MESH_Y_SIZE % GRASS_BLOCK_SZ) == 0);

	for (unsigned lod = 0; lod < NUM_GRASS_LODS; ++lod) {
		vbo_offsets[lod].resize(NUM_RND_GRASS_BLOCKS+1);
		vbo_offsets[lod][0] = grass.size(); // start
		for (unsigned i = 0; i < NUM_RND_GRASS_BLOCKS; ++i) {gen_lod_block(i, lod);}
	}
	PRINT_TIME("Grass Tile Gen");
}


void grass_tile_manager_t::update() { // to be called once per frame

	if (!is_grass_enabled()) {
		clear();
		return;
	}
	if (empty()    ) {gen_grass();}
	if (vbo == 0   ) {create_new_vbo();}
	if (!data_valid) {upload_data();}
}


// Note: density won't work if grass is spatially sorted, which is currently the case for tiled terrain grass
void grass_tile_manager_t::render_block(unsigned block_ix, unsigned lod, float density, unsigned instance_count) {

	assert(density > 0.0 && density <= 1.0);
	assert(lod < NUM_GRASS_LODS);
	assert(block_ix+1 < vbo_offsets[lod].size());
	unsigned const start_ix(vbo_offsets[lod][block_ix]), end_ix(vbo_offsets[lod][block_ix+1]);
	assert(start_ix < end_ix && end_ix <= grass.size());
	unsigned const num_tris(ceil(density*(end_ix - start_ix)));
	if (num_tris == 0) return;
	
	if (instance_count > 0) {
		bind_vbo(vbo); // needed because incoming vbo is 0 (so that instance attrib array isn't bound to a vbo)
		glDrawArraysInstanced(GL_TRIANGLES, 3*start_ix, 3*num_tris, instance_count);
	}
	else {
		glDrawArrays(GL_TRIANGLES, 3*start_ix, 3*num_tris);
	}
}


class grass_manager_dynamic_t : public grass_manager_t {
	
	vector<unsigned> mesh_to_grass_map; // maps mesh x,y index to starting index in grass vector
	vector<unsigned char> modified; // only used for shadows
	vector<int> last_occluder;
	mutable vector<grass_data_t> vertex_data_buffer;
	bool shadows_valid, has_voxel_grass;
	int last_cobj;
	int last_light;
	point last_lpos;

	bool hcm_chk(int x, int y) const {
		return (!point_outside_mesh(x, y) && (mesh_height[y][x] + SMALL_NUMBER < h_collision_matrix[y][x]));
	}
	void check_and_update_grass(unsigned ix, unsigned min_up, unsigned max_up) {
		if (min_up > max_up) return; // nothing updated
		modified[ix] |= MOD_GEOM; // usually few duplicates each frame, except for cluster grenade explosions
		if (vbo > 0) {upload_data_to_vbo(min_up, max_up+1, 0);}
		//data_valid = 0;
	}

public:
	grass_manager_dynamic_t() : shadows_valid(0), has_voxel_grass(0), last_light(-1), last_lpos(all_zeros) {}
	void invalidate_shadows()  {shadows_valid = 0;}
	
	void clear() {
		grass_manager_t::clear();
		invalidate_shadows();
		mesh_to_grass_map.clear();
		modified.clear();
	}

	bool ao_lighting_too_low(point const &pos) {
		float const keep_prob(5.0*(get_voxel_terrain_ao_lighting_val(pos) - 0.8));
		return (keep_prob < 0.0 || (keep_prob < 1.0 && rgen.rand_float() > keep_prob)); // too dark
	}

	void gen_grass() {
		RESET_TIME;
		float const *h_tex(island ? h_sand     : h_dirt);
		ttex  const *lttex(island ? lttex_sand : lttex_dirt);
		int   const   NTEX(island ? NTEX_SAND  : NTEX_DIRT);
		float const dz_inv(1.0/(zmax - zmin));
		mesh_to_grass_map.resize(XY_MULT_SIZE+1, 0);
		modified.resize(XY_MULT_SIZE, 0);
		object_types[GRASS].radius = 0.0;
		rgen.pregen_floats(10000);
		unsigned num_voxel_polys(0), num_voxel_blades(0);
		
		for (int y = 0; y < MESH_Y_SIZE; ++y) {
			for (int x = 0; x < MESH_X_SIZE; ++x) {
				mesh_to_grass_map[y*MESH_X_SIZE+x] = (unsigned)grass.size();
				float const xval(get_xval(x)), yval(get_yval(y));

				//if (create_voxel_landscape) {
				if (coll_objects.has_voxel_cobjs) {
					float const blades_per_area(grass_density/dxdy);
					coll_cell const &cell(v_collision_matrix[y][x]);
					cube_t const test_cube(xval-0.5*DX_VAL, xval+0.5*DX_VAL, yval-0.5*DY_VAL, yval+0.5*DY_VAL, mesh_height[y][x], czmax+grass_length);
					float const nz_thresh = 0.4;

					for (unsigned k = 0; k < cell.cvals.size(); ++k) {
						int const index(cell.cvals[k]);
						if (index < 0) continue;
						assert(unsigned(index) < coll_objects.size());
						coll_obj const &cobj(coll_objects[index]);
						if (cobj.type != COLL_POLYGON || cobj.cp.cobj_type != COBJ_TYPE_VOX_TERRAIN) continue;
						if (cobj.norm.z < nz_thresh)     continue; // not oriented upward
						if (!cobj.intersects(test_cube)) continue;
						assert(cobj.npoints == 3 || cobj.npoints == 4); // triangles and quads
						float const density_scale((cobj.norm.z - nz_thresh)/(1.0 - nz_thresh)); // FIXME: better to use vertex normals and interpolate?
						unsigned const num_blades(blades_per_area*density_scale*polygon_area(cobj.points, cobj.npoints) + 0.5);
						++num_voxel_polys;

						for (unsigned n = 0; n < num_blades; ++n) {
							float const r1(rgen.rand_float()), r2(rgen.rand_float()), sqrt_r1(sqrt(r1));
							unsigned const ptix((cobj.npoints == 4 && n < num_blades/2) ? 3 : 1); // handle both triangles and quads
							point const pos((1 - sqrt_r1)*cobj.points[0] + (sqrt_r1*(1 - r2))*cobj.points[ptix] + (sqrt_r1*r2)*cobj.points[2]);
							if (!test_cube.contains_pt(pos)) continue; // bbox test
							if (ao_lighting_too_low(pos))    continue; // too dark
							add_grass_blade(pos, 0.8, 0); // use cobj.norm instead of mesh normal?
							++num_voxel_blades;
							has_voxel_grass = 1;
						}
					} // for k
				}

				// create mesh grass
				if (default_ground_tex >= 0 && default_ground_tex != GROUND_TEX) continue; // no grass
				if (x == MESH_X_SIZE-1 || y == MESH_Y_SIZE-1) continue; // mesh not drawn
				if (is_mesh_disabled(x, y) || is_mesh_disabled(x+1, y) || is_mesh_disabled(x, y+1) || is_mesh_disabled(x+1, y+1)) continue; // mesh disabled
				if (mesh_height[y][x] < water_matrix[y][x])   continue; // underwater (make this dynamically update?)
				bool const do_cobj_check(hcm_chk(x, y) || hcm_chk(x+1, y) || hcm_chk(x, y+1) || hcm_chk(x+1, y+1));
				float const vnz(vertex_normals[y][x].z);
				float const *const sti(sthresh[0][0]);
				float slope_scale(1.0);
				if (vnz < sti[1]) {slope_scale = CLIP_TO_01((vnz - sti[0])/(sti[1] - sti[0]));} // handle steep slopes (dirt/rock texture replaces grass texture)
				if (slope_scale == 0.0) continue; // no grass
				assert(vnz > 0.0);
				unsigned const tile_density(round_fp(grass_density/vnz)); // slightly more grass on steep slopes so that we have equal density over the surface, not just the XY projection

				for (unsigned n = 0; n < tile_density; ++n) {
					float const xv(rgen.rand_uniform(xval, xval + DX_VAL));
					float const yv(rgen.rand_uniform(yval, yval + DY_VAL));
					float const mh(interpolate_mesh_zval(xv, yv, 0.0, 0, 1));
					point const pos(xv, yv, mh);

					if (default_ground_tex < 0 && zmin < zmax) {
						float const relh(relh_adj_tex + (mh - zmin)*dz_inv);
						int k1, k2;
						float t(0.0);
						get_tids(relh, NTEX-1, h_tex, k1, k2, &t); // t==0 => use k1, t==1 => use k2
						int const id1(lttex[k1].id), id2(lttex[k2].id);
						if (id1 != GROUND_TEX && id2 != GROUND_TEX) continue; // not ground texture
						float density(1.0);
						if (id1 != GROUND_TEX) {density = t;}
						if (id2 != GROUND_TEX) {density = 1.0 - t;}
						if (rgen.rand_float() >= slope_scale*density) continue; // skip - density too low
					}
					// skip grass intersecting cobjs
					if (do_cobj_check && dwobject(GRASS, pos).check_vert_collision(0, 0, 0)) continue; // make a GRASS object for collision detection

					if (create_voxel_landscape) {
						if (point_inside_voxel_terrain(pos)) continue; // inside voxel volume
						if (ao_lighting_too_low(pos))        continue; // too dark
					}
					add_grass_blade(pos, 0.8, 1);
				} // for n
			} // for x
		} // for y
		if (has_voxel_grass) {cout << "voxel_polys: " << num_voxel_polys << ", voxel_blades: " << num_voxel_blades << endl;}
		mesh_to_grass_map[XY_MULT_SIZE] = (unsigned)grass.size();
		remove_excess_cap(grass);
		PRINT_TIME("Grass Generation");
	}

	unsigned char get_shadow_bits(int cid) const { // unused
		if (cid < 0) return MESH_SHADOW;
		assert((unsigned)cid < coll_objects.size());
		return ((coll_objects[cid].status == COLL_DYNAMIC) ? DYNAMIC_SHADOW : OBJECT_SHADOW);
	}

	bool is_pt_shadowed(point const &pos, bool skip_dynamic) {
		int const light(get_light());

		// determine if grass can be shadowed based on mesh shadow
		// Note: if mesh is shadowed this does *not* mean that grass is also shadowed
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
		bool unshadowed(1);

		for (int y = max(0, ypos-1); y <= min(MESH_Y_SIZE-1, ypos+1); ++y) { // test 3x3 window around the point
			for (int x = max(0, xpos-1); x <= min(MESH_X_SIZE-1, xpos+1); ++x) {
				if (x != xpos && y != ypos) continue; // no diagonals - faster, slightly less accurate
				if ((shadow_mask[light][y][x] & SHADOWED_ALL) != 0) unshadowed = 0;
			}
		}
		if (unshadowed) return 0; // no shadows on mesh, so no shadows on grass

		if (last_cobj >= 0) { // check to see if last cobj still intersects
			assert((unsigned)last_cobj < coll_objects.size());
			point lpos;
			if (get_light_pos(lpos, light) && coll_objects[last_cobj].line_intersect(pos, lpos)) return 1;
		}
		if (is_visible_to_light_cobj(pos, light, 0.0, -1, skip_dynamic, &last_cobj)) return 0;
		return 1;
	}

	void find_shadows() {
		if (empty()) return;
		RESET_TIME;
		last_cobj     = -1;
		shadows_valid = 1;
		data_valid    = 0;

		for (unsigned i = 0; i < grass.size(); ++i) {
			grass[i].shadowed = is_pt_shadowed((grass[i].p + grass[i].dir*0.5), 1); // per vertex shadows?
		}
		PRINT_TIME("Grass Find Shadows");
	}

	void upload_data_to_vbo(unsigned start, unsigned end, bool alloc_data) const {
		if (start == end) return; // nothing to update
		assert(start < end && end <= grass.size());
		unsigned const num_verts(3*(end - start)), block_size(3*4096); // must be a multiple of 3
		unsigned const vntc_sz(sizeof(grass_data_t));
		unsigned offset(3*start);
		vertex_data_buffer.resize(min(num_verts, block_size));
		bind_vbo(vbo);
		if (alloc_data) {upload_vbo_data(NULL, 3*grass.size()*vntc_sz);} // initial upload (setup, no data)
		
		for (unsigned i = start, ix = 0; i < end; ++i) {
			//vector3d norm(grass[i].n); // use grass normal? 2-sided lighting?
			//vector3d norm(surface_normals[get_ypos(p1.y)][get_xpos(p1.x)]);
			vector3d norm(grass[i].on_mesh ? interpolate_mesh_normal(grass[i].p) : plus_z); // use +z normal for voxels
			add_to_vbo_data(grass[i], vertex_data_buffer, ix, norm);

			if (ix == block_size || i+1 == end) { // filled block or last entry
				upload_vbo_sub_data(&vertex_data_buffer.front(), offset*vntc_sz, ix*vntc_sz); // upload part or all of the data
				offset += ix;
				ix = 0; // reset to the beginning of the buffer
			}
		}
		assert(offset == 3*end);
		bind_vbo(0);
	}

	float get_xy_bounds(point const &pos, float radius, int &x1, int &y1, int &x2, int &y2) const {
		if (empty() || !is_over_mesh(pos)) return 0.0;
		assert(radius > 0.0);

		if (!has_voxel_grass) { // determine radius at grass height
			float const mh(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1));
			if ((pos.z - radius) > (mh + grass_length)) return 0.0; // above grass
			if ((pos.z + radius) < mh)                  return 0.0; // below the mesh
			float const height(pos.z - (mh + grass_length));
			if (height > 0.0) {radius = sqrt(max(1.0E-6f, (radius*radius - height*height)));}
		}
		x1 = get_xpos(pos.x - radius);
		x2 = get_xpos(pos.x + radius);
		y1 = get_ypos(pos.y - radius);
		y2 = get_ypos(pos.y + radius);
		return radius;
	}

	unsigned get_start_and_end(int x, int y, unsigned &start, unsigned &end) const {
		unsigned const ix(y*MESH_X_SIZE + x);
		assert(ix+1 < mesh_to_grass_map.size());
		start = mesh_to_grass_map[ix];
		end   = mesh_to_grass_map[ix+1];
		assert(start <= end && end <= grass.size());
		return ix;
	}

	bool place_obj_on_grass(point &pos, float radius) const {
		int x1, y1, x2, y2;
		float const rad(get_xy_bounds(pos, radius, x1, y1, x2, y2));
		if (rad == 0.0) return 0;
		bool updated(0);

		for (int y = y1; y <= y2; ++y) {
			for (int x = x1; x <= x2; ++x) {
				if (point_outside_mesh(x, y)) continue;
				unsigned start, end;
				get_start_and_end(x, y, start, end);
				if (start == end) continue; // no grass at this location

				for (unsigned i = start; i < end; ++i) {
					if (grass[i].dir == zero_vector) continue; // removed
					float const dsq(p2p_dist_xy_sq(pos, grass[i].p));
					if (dsq > rad*rad) continue; // too far away
					pos.z   = max(pos.z, (grass[i].p.z + grass[i].dir.z + radius));
					updated = 1;
				}
			}
		}
		return updated;
	}

	float get_grass_density(point const &pos) {
		if (empty() || !is_over_mesh(pos)) return 0.0;
		int const x(get_xpos(pos.x)), y(get_ypos(pos.y));
		if (point_outside_mesh(x, y))      return 0.0;
		unsigned const ix(y*MESH_X_SIZE + x);
		assert(ix+1 < mesh_to_grass_map.size());
		unsigned const num_grass(mesh_to_grass_map[ix+1] - mesh_to_grass_map[ix]); // Note: ignores crushed/cut/removed grass
		return ((float)num_grass)/((float)grass_density);
	}

	void mesh_height_change(int x, int y, bool recalc_shadows) {
		assert(!point_outside_mesh(x, y));
		unsigned start, end;
		unsigned const ix(get_start_and_end(x, y, start, end));
		unsigned min_up(end+1), max_up(start);
		if (shadow_map_enabled()) {recalc_shadows = 0;}

		for (unsigned i = start; i < end; ++i) { // will do nothing if there's no grass here
			grass_t &g(grass[i]);
			if (!g.on_mesh || g.dir == zero_vector) continue; // not on mesh, or already "removed"
			float const mh(interpolate_mesh_zval(g.p.x, g.p.y, 0.0, 0, 1));

			if (fabs(g.p.z - mh) > 0.01*grass_width) { // is there any way we can check the ground texture to see if we sill have grass texture here?
				g.p.z  = mh;
				min_up = min(min_up, i);
				max_up = max(max_up, i);
				// not entirely correct, since the mesh shadow data won't necessarily have been updated when this is called
				// it may be more correct to defer this processing, or recalculate all grass lighting, but that's both complex and very expensive
				if (recalc_shadows) {g.shadowed = is_pt_shadowed((g.p + g.dir*0.5), 1);} 
			}
		} // for i
		check_and_update_grass(ix, min_up, max_up);
	}

	void modify_grass(point const &pos, float radius, bool crush, bool burn, bool cut, bool check_uw, bool add_color, bool remove, colorRGBA const &color) {
		if (burn && is_underwater(pos)) burn = 0;
		if (!burn && !crush && !cut && !check_uw && !add_color && !remove) return; // nothing left to do
		int x1, y1, x2, y2;
		float const rad(get_xy_bounds(pos, radius, x1, y1, x2, y2));
		if (rad == 0.0) return;
		color_wrapper cw;
		if (add_color) {cw.set_c3(color);}

		// modify grass within radius of pos
		for (int y = y1; y <= y2; ++y) {
			for (int x = x1; x <= x2; ++x) {
				if (point_outside_mesh(x, y)) continue;
				bool const maybe_underwater((burn || check_uw) && has_water(x, y) && mesh_height[y][x] <= water_matrix[y][x]);
				unsigned start, end;
				unsigned const ix(get_start_and_end(x, y, start, end));
				unsigned min_up(end+1), max_up(start);

				for (unsigned i = start; i < end; ++i) { // will do nothing if there's no grass here
					grass_t &g(grass[i]);
					if (g.dir == zero_vector) continue; // already "removed"
					float const dsq(p2p_dist_xy_sq(pos, g.p));
					if (dsq > rad*rad) continue; // too far away
					float const reld(sqrt(dsq)/rad);
					bool const underwater(maybe_underwater && g.on_mesh);
					bool updated(0);

					if (cut) {
						float const length(g.dir.mag());

						if (length > 0.25*grass_length) {
							g.dir  *= reld;
							updated = 1;
						}
					}
					if (crush) {
						vector3d const &sn(surface_normals[y][x]);
						float const length(g.dir.mag());

						if (fabs(dot_product(g.dir, sn)) > 0.1*length) { // update if not flat against the mesh
							float const dx(g.p.x - pos.x), dy(g.p.y - pos.y), atten_val(1.0 - (1.0 - reld)*(1.0 - reld));
							vector3d const new_dir(vector3d(dx, dy, -(sn.x*dx + sn.y*dy)/sn.z).get_norm()); // point away from crushing point

							if (dot_product(g.dir, new_dir) < 0.95*length) { // update if not already aligned
								g.dir   = (g.dir*(atten_val/length) + new_dir*(1.0 - atten_val)).get_norm()*length;
								g.n     = (g.n*atten_val + sn*(1.0 - atten_val)).get_norm();
								updated = 1;
							}
						}
					}
					if (add_color && !underwater) {
						UNROLL_3X(updated |= (g.c[0] != cw.c[0] || g.c[1] != cw.c[1] || g.c[2] != cw.c[2]);) // not already red
						
						if (updated) {
							float const atten_val(1.0 - color.alpha*(1.0 - reld)*(1.0 - reld));
							UNROLL_3X(g.c[i_] = (unsigned char)(atten_val*g.c[i_] + (1.0 - atten_val)*cw.c[i_]);)
						}
					}
					if (burn && !underwater) {
						float const atten_val(1.0 - (1.0 - reld)*(1.0 - reld));
						UNROLL_3X(updated |= (g.c[i_] > 0);)
						if (updated) {UNROLL_3X(g.c[i_] = (unsigned char)(atten_val*g.c[i_]);)}
					}
					if (check_uw && underwater && (g.p.z + g.dir.mag()) <= water_matrix[y][x]) {
						unsigned char uwc[3] = {120,  100, 50};
						UNROLL_3X(updated |= (g.c[i_] != uwc[i_]);)
						if (updated) {UNROLL_3X(g.c[i_] = (unsigned char)(0.9*g.c[i_] + 0.1*uwc[i_]);)}
					}
					if (remove) {
						// Note: if we're removing, it doesn't make sense to do any other operations since they won't have any effect
						g.dir   = zero_vector; // make zero length (can't actually remove it)
						updated = 1;
					}
					if (updated) {
						min_up = min(min_up, i);
						max_up = max(max_up, i);
					}
				} // for i
				check_and_update_grass(ix, min_up, max_up);
			} // for x
		} // for y
	}

	void upload_data(bool alloc_data) {
		if (empty()) return;
		RESET_TIME;
		upload_data_to_vbo(0, (unsigned)grass.size(), alloc_data);
		data_valid = 1;
		PRINT_TIME("Grass Upload VBO");
		cout << "mem used: " << grass.size()*sizeof(grass_t) << ", vmem used: " << 3*grass.size()*sizeof(grass_data_t) << endl;
	}

	void check_for_updates() {
		bool const vbo_invalid(vbo == 0);
		if (!shadows_valid && !shadow_map_enabled()) {find_shadows();}
		if (vbo_invalid) {create_new_vbo();}
		if (!data_valid) {upload_data(vbo_invalid);}
	}

	void draw_range(unsigned beg_ix, unsigned end_ix) const {
		assert(beg_ix <= end_ix && end_ix <= grass.size());
		if (beg_ix < end_ix) {glDrawArrays(GL_TRIANGLES, 3*beg_ix, 3*(end_ix - beg_ix));} // nonempty segment
	}

	static void setup_shaders(shader_t &s, bool distant) { // per-pixel dynamic lighting
		s.setup_enabled_lights(2, 2); // FS; L0-L1: static directional
		set_dlights_booleans(s, 1, 1); // FS
		s.check_for_fog_disabled();
		if (distant) {s.set_prefix("#define NO_GRASS_TEXTURE", 1);} // FS
		s.set_prefix("#define NO_DL_SPECULAR",   1); // FS ???
		s.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
		s.set_bool_prefix("use_shadow_map", shadow_map_enabled(), 1); // FS
		s.set_vert_shader("wind.part*+grass_texture.part+grass_pp_dl");
		s.set_frag_shader("linear_fog.part+dynamic_lighting.part*+ads_lighting.part*+shadow_map.part*+grass_with_dlights");
		//s.set_vert_shader("ads_lighting.part*+shadow_map.part*+wind.part*+grass_texture.part+grass");
		//s.set_frag_shader("linear_fog.part+textured_with_fog");
		s.begin_shader();
		s.setup_scene_bounds();
		setup_dlight_textures(s);
		s.add_uniform_color("color_scale", (distant ? texture_color(GRASS_BLADE_TEX)*1.06 : WHITE));
		if (shadow_map_enabled()) {set_smap_shader_for_all_lights(s);}
		setup_wind_for_shader(s, 1);
		s.add_uniform_int("tex0", 0);
		s.setup_fog_scale();
		s.add_uniform_float("height", grass_length);
	}

	// texture units used: 0: grass texture, 1: wind texture
	void draw() {
		if (empty()) return;
		//RESET_TIME;

		// determine if ligthing has changed and possibly calculate shadows/upload VBO data
		int const light(get_light());
		point const lpos(get_light_pos());

		if (light != last_light || (!no_sun_lpos_update && lpos != last_lpos)) {
			invalidate_shadows();
			last_light = light;
			last_lpos  = lpos;
		}
		check_for_updates();

		// check for dynamic light sources
		bool const use_grass_shader(!has_snow && (display_mode & 0x0100));
		shader_t s;
		if (use_grass_shader) {setup_shaders(s, 1);} // enables lighting and shadows as well
		begin_draw(0.2);

		// draw the grass
		bool last_visible(0);
		unsigned beg_ix(0);
		point const camera(get_camera_pos()), adj_camera(camera + point(0.0, 0.0, 2.0*grass_length));
		last_occluder.resize(XY_MULT_SIZE, -1);
		int last_occ_used(-1);
		vector<unsigned> nearby_ixs;

		for (int y = 0; y < MESH_Y_SIZE; ++y) {
			for (int x = 0; x < MESH_X_SIZE; ++x) {
				unsigned const ix(y*MESH_X_SIZE + x);
				if (mesh_to_grass_map[ix] == mesh_to_grass_map[ix+1]) continue; // empty section
				point const mpos(get_mesh_xyz_pos(x, y));
				float const grass_zmax((has_voxel_grass ? max(mpos.z, czmax) : mpos.z) + grass_length);
				bool visible(1);

				if (x+1 < MESH_X_SIZE && y+1 < MESH_Y_SIZE &&
					dot_product(surface_normals[y  ][x  ], (adj_camera - mpos                      )) < 0.0 &&
					dot_product(surface_normals[y  ][x+1], (adj_camera - get_mesh_xyz_pos(x+1, y  ))) < 0.0 &&
					dot_product(surface_normals[y+1][x+1], (adj_camera - get_mesh_xyz_pos(x+1, y+1))) < 0.0 &&
					dot_product(surface_normals[y+1][x  ], (adj_camera - get_mesh_xyz_pos(x,   y+1))) < 0.0) // back_facing
				{
					visible = 0;
				}
				else {
					cube_t const cube(mpos.x-grass_length, mpos.x+DX_VAL+grass_length,
									  mpos.y-grass_length, mpos.y+DY_VAL+grass_length, z_min_matrix[y][x], grass_zmax);
					visible = camera_pdu.cube_visible(cube); // could use camera_pdu.sphere_and_cube_visible_test()
				
					if (visible && (display_mode & 0x08)) {
						int &last_occ_cobj(last_occluder[y*MESH_X_SIZE + x]);

						if (last_occ_cobj >= 0 || ((frame_counter + y) & 7) == 0) { // only sometimes update if not previously occluded
							point pts[8];
							get_cube_points(cube.d, pts);

							if (x > 0 && last_occ_cobj >= 0 && last_occluder[y*MESH_X_SIZE + x-1] == last_occ_cobj &&
								!coll_objects[last_occ_cobj].disabled() && coll_objects[last_occ_cobj].intersects_all_pts(camera, (pts+4), 4))
							{
								visible = 0; // check last 4 points (first 4 should already be checked from the last grass block)
							}
							else {
								if (last_occ_cobj < 0) {last_occ_cobj = last_occ_used;}
								visible &= !cobj_contained_ref(camera, cube.get_cube_center(), pts, 8, -1, last_occ_cobj);
								if (visible) {last_occ_cobj = -1;}
							}
							if (last_occ_cobj >= 0) {last_occ_used = last_occ_cobj;}
						}
					}
				}
				if (use_grass_shader && visible && dist_less_than(camera, mpos, 500.0*grass_width)) { // nearby grass
					nearby_ixs.push_back(ix);
					visible = 0; // drawn in the second pass
				}
				if (visible && !last_visible) { // start a segment
					beg_ix = mesh_to_grass_map[ix];
				}
				else if (!visible && last_visible) { // end a segment
					draw_range(beg_ix, mesh_to_grass_map[ix]);
				}
				last_visible = visible;
			}
		}
		if (last_visible) {draw_range(beg_ix, (unsigned)grass.size());}

		if (!nearby_ixs.empty()) {
			s.end_shader();
			setup_shaders(s, 0);

			for (vector<unsigned>::const_iterator i = nearby_ixs.begin(); i != nearby_ixs.end(); ++i) {
				draw_range(mesh_to_grass_map[*i], mesh_to_grass_map[(*i)+1]);
			}
		}
		s.end_shader();
		end_draw();
		//PRINT_TIME("Draw Grass");
	}
};

grass_manager_dynamic_t grass_manager;


void setup_wind_for_shader(shader_t &s, unsigned tu_id) {

	static float time(0.0);
	if (animate2) time = tfticks;
	s.add_uniform_float("time", 0.5*time/TICKS_PER_SECOND);
	s.add_uniform_float("wind_x", wind.x);
	s.add_uniform_float("wind_y", wind.y);
	s.add_uniform_int("wind_noise_tex", tu_id);
	select_multitex(WIND_TEX, tu_id, 0);
}


bool no_grass() {
	return (grass_density == 0 || !grass_enabled || snow_enabled() || vegetation == 0.0 || read_landscape);
}


void gen_grass(bool full_regen) {

	if (!full_regen) { // update shadows only
		grass_manager.invalidate_shadows();
		return;
	}
	grass_manager.clear();
	if (no_grass()) return;
	grass_manager.gen_grass();
	cout << "grass: " << grass_manager.size() << " out of " << XY_MULT_SIZE*grass_density << endl;
}


void update_grass_vbos() {
	grass_manager.free_vbo();
}

void draw_grass() {
	if (!no_grass() && (display_mode & 0x02)) {grass_manager.draw();}
}

void modify_grass_at(point const &pos, float radius, bool crush, bool burn, bool cut, bool check_uw, bool add_color, bool remove, colorRGBA const &color) {
	if (!no_grass()) {grass_manager.modify_grass(pos, radius, crush, burn, cut, check_uw, add_color, remove, color);}
}

void grass_mesh_height_change(int xpos, int ypos, bool recalc_shadows) {
	if (!no_grass()) {grass_manager.mesh_height_change(xpos, ypos, recalc_shadows);}
}

bool place_obj_on_grass(point &pos, float radius) {
	return (!no_grass() && grass_manager.place_obj_on_grass(pos, radius));
}

float get_grass_density(point const &pos) {
	return (no_grass() ? 0.0 : grass_manager.get_grass_density(pos));
}



