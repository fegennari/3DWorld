// 3D World - Grass Generation and Rendering Code
// by Frank Gennari
// 9/28/10

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"
#include "textures_3dw.h"
#include "lightmap.h"
#include "gl_ext_arb.h"


#define MOD_GEOM   0x01
#define MOD_SHADOW 0x02

bool grass_enabled(1);
unsigned grass_density(0);
float grass_length(0.02), grass_width(0.002);

extern bool has_dir_lights, has_snow, disable_shaders;
extern int island, default_ground_tex, read_landscape, display_mode, animate2;
extern float vegetation, zmin, zmax, fticks, tfticks, h_sand[], h_dirt[], leaf_color_coherence, tree_deadness, relh_adj_tex;
extern colorRGBA leaf_base_color;
extern vector3d wind;
extern obj_type object_types[];
extern vector<coll_obj> coll_objects;
extern vector<light_source> dl_sources;
extern vector<shadow_sphere> shadow_objs;


bool snow_enabled();


class grass_manager_t {
	
	struct grass_t { // size = 48
		point p;
		vector3d dir, n;
		unsigned char c[3];
		unsigned char shadowed;
		float w;

		grass_t() {} // optimization
		grass_t(point const &p_, vector3d const &dir_, vector3d const &n_, unsigned char const *const c_, float w_)
			: p(p_), dir(dir_), n(n_), shadowed(0), w(w_) {c[0] = c_[0]; c[1] = c_[1]; c[2] = c_[2];}
	};

	vector<grass_t> grass;
	vector<unsigned> mesh_to_grass_map; // maps mesh x,y index to starting index in grass vector
	vector<unsigned char> modified; // only used for shadows
	vector<vector<unsigned> > sbins;
	unsigned vbo;
	bool vbo_valid, shadows_valid, data_valid;
	int last_cobj;
	int last_light;
	point last_lpos;

	bool hcm_chk(int x, int y) const {
		return (!point_outside_mesh(x, y) && (mesh_height[y][x] + SMALL_NUMBER < h_collision_matrix[y][x]));
	}

public:
	grass_manager_t() : vbo(0), vbo_valid(0), shadows_valid(0), data_valid(0), last_light(-1), last_lpos(all_zeros) {}
	~grass_manager_t() {clear();}
	size_t size() const {return grass.size() ;} // 2 points per grass blade
	bool empty()  const {return grass.empty();}
	void invalidate_vbo()     {vbo_valid     = 0;}
	void invalidate_shadows() {shadows_valid = 0;}
	
	void clear() {
		delete_vbo(vbo);
		vbo = 0;
		vbo_valid = shadows_valid = data_valid = 0;
		grass.clear();
		mesh_to_grass_map.clear();
		modified.clear();
	}

	void gen_grass() {
		RESET_TIME;
		float const *h_tex(island ? h_sand     : h_dirt);
		ttex  const *lttex(island ? lttex_sand : lttex_dirt);
		int   const   NTEX(island ? NTEX_SAND  : NTEX_DIRT);
		float const dz_inv(1.0/(zmax - zmin));
		mesh_to_grass_map.resize(XY_MULT_SIZE+1, 0);
		modified.resize(XY_MULT_SIZE, 0);
		
		for (int y = 0; y < MESH_Y_SIZE; ++y) {
			for (int x = 0; x < MESH_X_SIZE; ++x) {
				mesh_to_grass_map[y*MESH_X_SIZE+x] = grass.size();
				if (x == MESH_X_SIZE-1 || y == MESH_Y_SIZE-1) continue; // mesh not drawn
				if (is_mesh_disabled(x, y) || is_mesh_disabled(x+1, y) || is_mesh_disabled(x, y+1) || is_mesh_disabled(x+1, y+1)) continue; // mesh disabled
				if (mesh_height[y][x] < water_matrix[y][x]) continue; // underwater (make this dynamically update?)
				float const xval(get_xval(x)), yval(get_yval(y));

				for (unsigned n = 0; n < grass_density; ++n) {
					float const xv(rand_uniform2(xval, xval + DX_VAL));
					float const yv(rand_uniform2(yval, yval + DY_VAL));
					float const mh(interpolate_mesh_zval(xv, yv, 0.0, 0, 1));
					point const pos(xv, yv, mh);

					if (default_ground_tex >= 0 || zmax == zmin) {
						if (default_ground_tex >= 0 && default_ground_tex != GROUND_TEX) continue;
					}
					else {
						float const relh(relh_adj_tex + (mh - zmin)*dz_inv);
						int k1, k2;
						float t;
						get_tids(relh, NTEX-1, h_tex, k1, k2, t); // t==0 => use k1, t==1 => use k2
						int const id1(lttex[k1].id), id2(lttex[k2].id);
						if (id1 != GROUND_TEX && id2 != GROUND_TEX) continue; // not ground texture
						float density(1.0);
						if (id1 != GROUND_TEX) density = t;
						if (id2 != GROUND_TEX) density = 1.0 - t;
						if (rand_float() >= density) continue; // skip - density too low
					}
					if (hcm_chk(x, y) || hcm_chk(x+1, y) || hcm_chk(x, y+1) || hcm_chk(x+1, y+1)) { // skip grass intersecting cobjs
						dwobject obj(GRASS, pos); // make a GRASS object for collision detection
						object_types[GRASS].radius = 0.0;
						if (obj.check_vert_collision(0, 0, 0)) continue;
					}
					add_grass(pos);
				}
			}
		}
		mesh_to_grass_map[XY_MULT_SIZE] = grass.size();
		remove_excess_cap(grass);
		PRINT_TIME("Grass Generation");
	}

	void add_grass(point const &pos) {
		vector3d const base_dir(plus_z);
		//vector3d const base_dir(interpolate_mesh_normal(pos));
		vector3d const dir((base_dir + signed_rand_vector(0.3) + wind*0.3).get_norm()); // make dynamic based on local wind?
		vector3d const norm(cross_product(dir, signed_rand_vector()).get_norm());

		//(0.1, 0.35), (0.5, 0.75), (0.0, 0.1) // untextured white triangle
		//unsigned char color[3] = {75+rand()%50, 150+rand()%50, 25+rand()%20}; // (0.3,0.5), (0.6,0.8), (0.1,0.18)
		float const ilch(1.0 - leaf_color_coherence), dead_scale(CLIP_TO_01(tree_deadness));
		float const base_color[3] = {0.3,  0.6, 0.08};
		float const mod_color [3] = {0.2,  0.2, 0.08};
		float const lbc_mult  [3] = {0.2,  0.4, 0.0 };
		float const dead_color[3] = {0.75, 0.6, 0.0 };
		unsigned char color[3];

		for (unsigned i = 0; i < 3; ++i) {
			float const ccomp(CLIP_TO_01(base_color[i] + lbc_mult[i]*leaf_base_color[i] + ilch*mod_color[i]*rand_float()));
			color[i] = (unsigned char)(255.0*(dead_scale*dead_color[i] + (1.0 - dead_scale)*ccomp));
		}
		float const length(grass_length*rand_uniform2(0.7, 1.3));
		float const width( grass_width *rand_uniform2(0.7, 1.3));
		grass.push_back(grass_t(pos, dir*length, norm, color, width));
	}

	unsigned char get_shadow_bits(int cid) const {
		if (cid < 0) return MESH_SHADOW;
		assert((unsigned)cid < coll_objects.size());
		return ((coll_objects[cid].status == COLL_DYNAMIC) ? DYNAMIC_SHADOW : OBJECT_SHADOW);
	}

	unsigned char is_pt_shadowed(point const &pos, bool skip_dynamic) {
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
			if (get_light_pos(lpos, light) && coll_objects[last_cobj].line_intersect(pos, lpos)) return get_shadow_bits(last_cobj);
		}
		if (is_visible_to_light_cobj(pos, light, 0.0, -1, skip_dynamic, &last_cobj)) return 0;
		return get_shadow_bits(last_cobj);
	}

	void find_shadows() {
		if (empty()) return;
		RESET_TIME;
		last_cobj     = -1;
		shadows_valid = 1;
		data_valid    = 0;
		update_cobj_tree();

		for (unsigned i = 0; i < grass.size(); ++i) {
			grass[i].shadowed = is_pt_shadowed(grass[i].p + grass[i].dir*0.5, 1); // per vertex shadows?
		}
		PRINT_TIME("Grass Find Shadows");
	}

	void create_new_vbo() {
		delete_vbo(vbo);
		vbo        = create_vbo();
		vbo_valid  = 1;
		data_valid = 0;
	}

	vector3d interpolate_mesh_normal(point const &pos) const {
		float const xp((pos.x + X_SCENE_SIZE)*DX_VAL_INV), yp((pos.y + Y_SCENE_SIZE)*DY_VAL_INV);
		int const x0((int)xp), y0((int)yp);
		if (point_outside_mesh(x0, y0) || point_outside_mesh(x0+1, y0+1)) return plus_z; // shouldn't get here
		float const xpi(fabs(xp - (float)x0)), ypi(fabs(yp - (float)y0));
		return vertex_normals[y0+0][x0+0]*((1.0 - xpi)*(1.0 - ypi))
			 + vertex_normals[y0+1][x0+0]*((1.0 - xpi)*ypi)
			 + vertex_normals[y0+0][x0+1]*(xpi*(1.0 - ypi))
		     + vertex_normals[y0+1][x0+1]*(xpi*ypi);
	}

	void upload_data_to_vbo(unsigned start, unsigned end, bool create) const {
		if (start == end) return; // nothing to update
		assert(start < end && end <= grass.size());
		unsigned const num_verts(3*(end - start));
		unsigned const block_size(3*4096); // must be a multiple of 3
		unsigned const vntc_sz(sizeof(vert_norm_tc_color));
		unsigned offset(3*start);
		vector<vert_norm_tc_color> data(min(num_verts, block_size));
		bind_vbo(vbo);

		if (create) { // initial upload (setup, no data)
			upload_vbo_data(NULL, 3*grass.size()*vntc_sz);
		}
		for (unsigned i = start, ix = 0; i < end; ++i) {
			grass_t const &g(grass[i]);
			point const p1(g.p), p2(p1 + g.dir + point(0.0, 0.0, 0.05*grass_length));
			vector3d const binorm(cross_product(g.dir, g.n).get_norm());
			vector3d const delta(binorm*(0.5*g.w));
			float const nmag(g.shadowed ? 0.001 : 1.0);
			//vector3d const &norm(plus_z*nmag); // use grass normal? 2-sided lighting?
			//vector3d const &norm(g.n*nmag);
			//vector3d const &norm(surface_normals[get_ypos(p1.y)][get_xpos(p1.x)]*nmag);
			vector3d const &norm(interpolate_mesh_normal(p1)*nmag);
			float const tc_adj(0.1); // border around grass blade texture
			data[ix++].assign(p1-delta, norm, 1.0-tc_adj,     tc_adj, g.c);
			data[ix++].assign(p1+delta, norm, 1.0-tc_adj, 1.0-tc_adj, g.c);
			data[ix++].assign(p2,       norm,     tc_adj, 0.5,        g.c);
			assert(ix <= data.size());

			if (ix == block_size || i+1 == end) { // filled block or last entry
				upload_vbo_sub_data(&data.front(), offset*vntc_sz, ix*vntc_sz); // upload part or all of the data
				offset += ix;
				ix = 0; // reset to the beginning of the buffer
			}
		}
		assert(offset == 3*end);
		bind_vbo(0);
	}

	float get_xy_bounds(point const &pos, float radius, int &x1, int &y1, int &x2, int &y2) const {
		if (empty() || !is_over_mesh(pos)) return 0.0;

		// determine radius at grass height
		assert(radius > 0.0);
		float const mh(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1));
		if ((pos.z - radius) > (mh + grass_length)) return 0.0; // above grass
		if ((pos.z + radius) < mh)                  return 0.0; // below the mesh
		float const height(pos.z - (mh + grass_length));
		float const rad((height > 0.0) ? sqrt(max(1.0E-6f, (radius*radius - height*height))) : radius);
		x1 = get_xpos(pos.x - rad);
		x2 = get_xpos(pos.x + rad);
		y1 = get_ypos(pos.y - rad);
		y2 = get_ypos(pos.y + rad);
		return rad;
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
		unsigned const num_grass(mesh_to_grass_map[ix+1] - mesh_to_grass_map[ix]);
		return ((float)num_grass)/((float)grass_density);
	}

	void proc_dynamic_shadows(int light) {
		//RESET_TIME;
		point lpos;
		if (!get_light_pos(lpos, light)) return;
		sbins.resize(XY_MULT_SIZE); // move into constructor?
		last_cobj = -1;
		
		for (unsigned i = 0; i < shadow_objs.size(); ++i) {
			float const radius(shadow_objs[i].radius);
			if (radius < 2*LARGE_OBJ_RAD) continue; // for efficiency
			point const &pos(shadow_objs[i].pos);
			vector3d const dir((pos - lpos).get_norm());
			point pts[8];
			get_sphere_points(pos, radius, pts);
			int x1, y1, x2, y2, ret_val;
			get_shape_shadow_bb(pts, 8, light, 1, lpos, x1, y1, x2, y2, ret_val, DYNAMIC_SHADOW);

			for (int y = y1; y <= y2; ++y) {
				for (int x = x1; x <= x2; ++x) {
					assert(!point_outside_mesh(x, y));
					point const mesh_pos(point(get_xval(x), get_yval(y), mesh_height[y][x]));
					if (!line_sphere_intersect(lpos, mesh_pos, pos, (radius + HALF_DXY))) continue; // only helps slightly
					if (!camera_pdu.sphere_visible_test(mesh_pos, HALF_DXY))              continue; // only helps slightly
					sbins[y*MESH_X_SIZE + x].push_back(i);
				}
			}
		}
		for (int y = 0; y < MESH_Y_SIZE; ++y) { // check for dynamic shadows that need to be updated
			for (int x = 0; x < MESH_X_SIZE; ++x) {
				unsigned const ix(y*MESH_X_SIZE + x);
				bool const mod((modified[ix] & MOD_SHADOW) != 0);
				modified[ix] &= ~MOD_SHADOW;
				// check camera view frustum?
				if (mod || !sbins[ix].empty()) update_shadows(x, y, DYNAMIC_SHADOW, lpos, &sbins[ix]);
				sbins[ix].resize(0); // clear
			}
		}
		//PRINT_TIME("Grass Dynamic Shadows");
	}

	void update_shadows_for_cube(cube_t const &cube) {
		//RESET_TIME;
		int const light(get_light());
		point lpos;
		if (!get_light_pos(lpos, light)) return;
		point pts[8];
		get_cube_points(cube.d, pts);
		int x1, y1, x2, y2, ret_val;
		get_shape_shadow_bb(pts, 8, light, 1, lpos, x1, y1, x2, y2, ret_val, OBJECT_SHADOW);
		cube_t test_cube(cube);
		test_cube.expand_by(HALF_DXY);
		last_cobj = -1;

		for (int y = y1; y <= y2; ++y) {
			for (int x = x1; x <= x2; ++x) {
				assert(!point_outside_mesh(x, y));
				point const mesh_pos(point(get_xval(x), get_yval(y), mesh_height[y][x]));
				if (!check_line_clip(lpos, mesh_pos, test_cube.d)) continue;
				update_shadows(x, y, OBJECT_SHADOW, lpos);
			}
		}
		//PRINT_TIME("Grass Shadow Update");
	}

	void update_shadows(int x, int y, unsigned char shad_types, point const &lpos, vector<unsigned> const *const sixs=0) {
		if (point_outside_mesh(x, y)) return;
		unsigned start, end;
		unsigned const ix(get_start_and_end(x, y, start, end));
		if (start == end) return; // no grass at this location
		bool const REBUILD_TREE = 0;
		if (REBUILD_TREE) update_cobj_tree();
		unsigned min_up(end), max_up(start);
		bool const skip_dynamic((shad_types & DYNAMIC_SHADOW) == 0);

		for (unsigned i = start; i < end; ++i) {
			grass_t &g(grass[i]);
			unsigned char orig_shadowed(g.shadowed);
			g.shadowed &= ~shad_types; // unset these bits since they are now invalid
			if (g.shadowed) continue; // already shadowed with a non-checked type
			point const pos(g.p + g.dir*0.5);

			if (sixs) {
				for (vector<unsigned>::const_iterator j = sixs->begin(); j != sixs->end(); ++j) {				
					if (shadow_objs[*j].line_intersect(lpos, pos)) {
						g.shadowed |= shad_types;
						break;
					}
				}
			}
			else {
				g.shadowed = is_pt_shadowed(pos, skip_dynamic); // per vertex shadows?
			}
			if (g.shadowed & shad_types) modified[ix] |= MOD_SHADOW;
			
			if ((g.shadowed != 0) != (orig_shadowed != 0)) { // shadow nonzero-ness has changed, need to update
				min_up = min(min_up, i);
				max_up = max(max_up, i);
			}
		} // for i
		if (min_up > max_up) return; // nothing updated
		if (vbo_valid) upload_data_to_vbo(min_up, max_up+1, 0);
	}

	void modify_grass(point const &pos, float radius, bool crush, bool burn, bool cut, bool update_mh, bool check_uw) {
		if (burn && is_underwater(pos)) burn = 0;
		if (!burn && !crush && !cut && !update_mh && !check_uw) return; // nothing left to do
		int x1, y1, x2, y2;
		float const rad(get_xy_bounds(pos, radius, x1, y1, x2, y2));
		if (rad == 0.0) return;

		// modify grass within radius of pos
		for (int y = y1; y <= y2; ++y) {
			for (int x = x1; x <= x2; ++x) {
				if (point_outside_mesh(x, y)) continue;
				bool const underwater((burn || check_uw) && has_water(x, y) && mesh_height[y][x] <= water_matrix[y][x]);
				unsigned start, end;
				unsigned const ix(get_start_and_end(x, y, start, end));
				if (start == end) continue; // no grass at this location
				unsigned min_up(end), max_up(start);

				for (unsigned i = start; i < end; ++i) {
					grass_t &g(grass[i]);
					float const dsq(p2p_dist_xy_sq(pos, g.p));
					if (dsq > rad*rad) continue; // too far away
					float const reld(sqrt(dsq)/rad);
					bool updated(0);

					if (update_mh) {
						float const mh(interpolate_mesh_zval(g.p.x, g.p.y, 0.0, 0, 1));

						if (fabs(g.p.z - mh) > 0.01*grass_width) {
							g.p.z   = mh;
							updated = 1;
						}
					}
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
					if (updated) {
						min_up = min(min_up, i);
						max_up = max(max_up, i);
					}
				} // for i
				if (min_up > max_up) continue; // nothing updated
				modified[ix] |= MOD_GEOM; // usually few duplicates each frame, except for cluster grenade explosions
				if (vbo_valid) upload_data_to_vbo(min_up, max_up+1, 0);
				//data_valid = 0;
			} // for x
		} // for y
	}

	void heal_grass() {
		// WRITE: fix dir, length, and color for grass that has been crushed or burned using modified
	}

	void upload_data() {
		if (empty()) return;
		RESET_TIME;
		upload_data_to_vbo(0, grass.size(), 1);
		data_valid = 1;
		PRINT_TIME("Grass Upload VBO");
		cout << "mem used: " << grass.size()*sizeof(grass_t) << ", vmem used: " << 3*grass.size()*sizeof(vert_norm_tc_color) << endl;
	}

	void check_for_updates() {
		if (!shadows_valid) find_shadows();
		if (!vbo_valid    ) create_new_vbo();
		if (!data_valid   ) upload_data();
	}

	void draw() {
		if (empty()) return;

		// determine if ligthing has changed and possibly calculate shadows/upload VBO data
		int const light(get_light());
		point const lpos(get_light_pos());

		if (light != last_light || lpos != last_lpos) {
			invalidate_shadows();
			last_light = light;
			last_lpos  = lpos;
		}
		else if (!(display_mode & 0x10)) {
			proc_dynamic_shadows(light);
		}
		check_for_updates();

		// check for dynamic light sources
		bool const grass_wind(!disable_shaders && !has_snow && (display_mode & 0x0100) != 0);
		unsigned const num_dlights(enable_dynamic_lights());

		if (grass_wind) {
			set_shader_prefix("#define USE_LIGHT_COLORS",  0); // VS
			set_shader_prefix("#define USE_GOOD_SPECULAR", 0); // VS
#if 0 // per-pixel dynamic lighting - looks better, but slow
			setup_enabled_lights(2); // L0-L1: static directional
			set_dlights_booleans(1, 1); // FS
			unsigned const p(set_shader_prog("ads_lighting.part*+wind.part*+grass_pp_dl", "linear_fog.part+dynamic_lighting.part*+grass_with_dlights"));
			setup_scene_bounds(p);
			setup_dlight_textures(p);
#else // per-vertex dynamic lighting, limited to 6 lights - faster
			setup_enabled_lights(8); // L0-L1: static directional, L2-L7: dynamic point
			unsigned const p(set_shader_prog("ads_lighting.part*+wind.part*+grass", "linear_fog.part+simple_texture"));
#endif
			setup_wind_for_shader(p);
			setup_fog_scale(p);
			add_uniform_float(p, "height", grass_length);
		}

		// draw the grass
		set_array_client_state(1, 1, 1, 1);
		assert(vbo_valid && vbo > 0);
		bind_vbo(vbo);
		vert_norm_tc_color::set_vbo_arrays();
		select_multitex(GRASS_BLADE_TEX, 0);
		enable_blend();
		set_specular(0.2, 20.0);
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.75);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_NORMALIZE);
		glDrawArrays(GL_TRIANGLES, 0, 3*grass.size());
		glDisable(GL_COLOR_MATERIAL);
		glEnable(GL_NORMALIZE);
		set_specular(0.0, 1.0);
		disable_blend();
		glDisable(GL_ALPHA_TEST);
		if (grass_wind) unset_shader_prog();
		disable_multitex_a();
		disable_dynamic_lights(num_dlights);
		bind_vbo(0);
		check_gl_error(40);
	}
};

grass_manager_t grass_manager;


void setup_wind_for_shader(unsigned p) {

	static float time(0.0);
	if (animate2) time = tfticks;
	add_uniform_float(p, "time", 0.5*time/TICKS_PER_SECOND);
	add_uniform_float(p, "wind_x", wind.x);
	add_uniform_float(p, "wind_y", wind.y);
	add_uniform_int(p, "tex0", 0);
	add_uniform_int(p, "tex_noise", 1);
	select_multitex(CLOUD_RAW_TEX, 1, 0);
	set_multitex(0);
}


bool no_grass() {
	return (grass_density == 0 || !grass_enabled || snow_enabled() || vegetation == 0.0 || read_landscape);
}


void gen_grass(bool full_regen) {

	if (!full_regen) { // update shadows only
		grass_manager.invalidate_shadows();
		return;
	}
	bool const use_vbos(setup_gen_buffers());
		
	if (!use_vbos) {
		cout << "Warning: VBOs not supported, so grass cannot be enabled." << endl;
		grass_enabled = 0;
	}
	grass_manager.clear();
	if (no_grass()) return;
	grass_manager.gen_grass();
	cout << "grass: " << grass_manager.size() << " out of " << XY_MULT_SIZE*grass_density << endl;
}


void update_grass_vbos() {
	grass_manager.invalidate_vbo();
}

void draw_grass() {
	if (!no_grass() && (display_mode & 0x02)) grass_manager.draw();
}

void modify_grass_at(point const &pos, float radius, bool crush, bool burn, bool cut, bool update_mh, bool check_uw) {
	if (!no_grass()) grass_manager.modify_grass(pos, radius, crush, burn, cut, update_mh, check_uw);
}

void update_grass_shadows(int x, int y, unsigned char check_shad_types) {
	if (!no_grass()) grass_manager.update_shadows(x, y, check_shad_types, get_light_pos());
}

void update_grass_shadows_for_cube(cube_t const &cube) {
	if (!no_grass()) grass_manager.update_shadows_for_cube(cube);
}

bool place_obj_on_grass(point &pos, float radius) {
	return (!no_grass() && grass_manager.place_obj_on_grass(pos, radius));
}

float get_grass_density(point const &pos) {
	return (no_grass() ? 0.0 : grass_manager.get_grass_density(pos));
}



