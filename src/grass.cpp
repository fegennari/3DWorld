// 3D World - Grass Generation and Rendering Code
// by Frank Gennari
// 9/28/10

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"


float const GRASS_LENGTH = 0.02;
float const GRASS_WIDTH  = 0.002;


bool grass_enabled(1);
unsigned grass_density(0);

extern int island, default_ground_tex, read_landscape;
extern float vegetation, zmin, zmax, h_sand[], h_dirt[];
extern vector3d wind;
extern obj_type object_types[];
extern vector<coll_obj> coll_objects;


bool snow_enabled();


class grass_manager_t {
	
	struct grass_t { // size = 48
		point p;
		vector3d dir, n;
		unsigned char c[3];
		float w;

		grass_t(point const &p_, vector3d const &dir_, vector3d const &n_, unsigned char const *const c_, float w_)
			: p(p_), dir(dir_), n(n_), w(w_) {c[0] = c_[0]; c[1] = c_[1]; c[2] = c_[2];}
	};

	vector<grass_t> grass;
	vector<unsigned> mesh_to_grass_map; // maps mesh x,y index to starting index in grass vector
	vector<unsigned char> modified;
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
				if (x == MESH_X_SIZE-1 || y == MESH_Y_SIZE-1 || is_mesh_disabled(x, y)) continue; // mesh not drawn
				if (mesh_height[y][x] < water_matrix[y][x]) continue; // underwater (make this dynamically update?)
				float const xval(get_xval(x)), yval(get_yval(y));

				for (unsigned n = 0; n < grass_density; ++n) {
					float const xv(rand_uniform(xval, xval + DX_VAL));
					float const yv(rand_uniform(yval, yval + DY_VAL));
					float const mh(interpolate_mesh_zval(xv, yv, 0.0, 0, 1));
					point const pos(xv, yv, mh);

					if (default_ground_tex >= 0 || zmax == zmin) {
						if (default_ground_tex >= 0 && default_ground_tex != GROUND_TEX) continue;
					}
					else {
						float const relh((mh - zmin)*dz_inv);
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
		PRINT_TIME("Grass Generation");
	}

	void add_grass(point const &pos) {
		vector3d const dir((plus_z + signed_rand_vector(0.3) + wind*0.3).get_norm()); // FIXME: make dynamic? local wind?
		vector3d const norm(cross_product(dir, signed_rand_vector()).get_norm());
		// Vary color per vertex? Add precomputed lighting to color?
		//(0.1, 0.35), (0.5, 0.75), (0.0, 0.1) // untextured white triangle
		unsigned char const color[3] = {75+rand()%50, 150+rand()%50, 25+rand()%20};
		float const length(GRASS_LENGTH*rand_uniform(0.7, 1.3));
		float const width( GRASS_WIDTH *rand_uniform(0.7, 1.3));
		grass.push_back(grass_t(pos, dir*length, norm, color, width));
	}

	bool is_pt_shadowed(point const &pos) {
		int const light(get_light());
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
		bool shadowed(0), unshadowed(0);

		for (int y = max(0, ypos-1); y <= min(MESH_Y_SIZE-1, ypos+1); ++y) { // test 3x3 window around the point
			for (int x = max(0, xpos-1); x <= min(MESH_X_SIZE-1, xpos+1); ++x) {
				if (x != xpos && y != ypos) continue; // no diagonals - faster, slightly less accurate
				(((shadow_mask[light][y][x] & SHADOWED_ALL) != 0) ? shadowed : unshadowed) = 1;
			}
		}
		if (shadowed != unshadowed) return shadowed; // only one was set, so mesh is in agreement

		if (last_cobj >= 0) { // check to see if last cobj still intersects
			assert(last_cobj < (int)coll_objects.size());
			point lpos;
			if (get_light_pos(lpos, light) && coll_objects[last_cobj].line_intersect(pos, lpos)) return 1;
		}
		return !is_visible_to_light_cobj(pos, light, 0.0, -1, 1, &last_cobj); // neither (off the mesh) or both (conflict)
	}

	void find_shadows() {
		if (empty()) return;
		RESET_TIME;
		last_cobj     = -1;
		shadows_valid = 1;
		data_valid    = 0;

		for (unsigned i = 0; i < grass.size(); ++i) {
			point const p1(grass[i].p), p2(p1 + grass[i].dir);
			bool const shadowed(is_pt_shadowed((p1 + p2)*0.5)); // per vertex shadows?
			grass[i].n.normalize();
			grass[i].n *= (shadowed ? 0.001 : 1.0); // set normal near zero if shadowed
		}
		PRINT_TIME("Grass Find Shadows");
	}

	void create_new_vbo() {
		delete_vbo(vbo);
		vbo        = create_vbo();
		vbo_valid  = 1;
		data_valid = 0;
	}

	void upload_data_to_vbo(unsigned start, unsigned end) const {
		if (start == end) return; // nothing to update
		assert(start < end && end <= grass.size());
		vector<vert_norm_tc_color> data;
		data.reserve(3*(end - start));

		for (unsigned i = start; i < end; ++i) {
			point const p1(grass[i].p), p2(p1 + grass[i].dir + point(0.0, 0.0, 0.05*GRASS_LENGTH));
			vector3d const binorm(cross_product(grass[i].dir, grass[i].n).get_norm());
			vector3d const delta(binorm*(0.5*grass[i].w));
			//vector3d const norm(grass[i].n);
			vector3d const norm(plus_z*grass[i].n.mag());
			float const tc_adj(0.1); // border around grass blade texture
			data.push_back(vert_norm_tc_color(p1-delta, norm, 1.0-tc_adj,     tc_adj, grass[i].c));
			data.push_back(vert_norm_tc_color(p1+delta, norm, 1.0-tc_adj, 1.0-tc_adj, grass[i].c));
			data.push_back(vert_norm_tc_color(p2,       norm,     tc_adj, 0.5,        grass[i].c));
		}
		bind_vbo(vbo);

		if (start == 0 && end == grass.size()) { // full data, do full upload
			upload_vbo_data(&data.front(), data.size()*sizeof(vert_norm_tc_color));
		}
		else { // partial data, upload a subset
			upload_vbo_sub_data(&data.front(), 3*start*sizeof(vert_norm_tc_color), data.size()*sizeof(vert_norm_tc_color));
		}
		bind_vbo(0);
	}

	void modify_grass(point const &pos, float radius, bool crush, bool burn) {
		if (empty() || !is_over_mesh(pos)) return;

		// determine radius at grass height
		float const mh(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1));
		if ((pos.z - radius) > (mh + GRASS_LENGTH)) return; // above grass
		if ((pos.z + radius) < mh) return; // below the mesh
		float const height(pos.z - (mh + GRASS_LENGTH));
		float const rad((height > 0.0) ? sqrt(radius*radius - height*height) : radius), rad_sq(rad*rad);
		int const x1(get_xpos(pos.x - rad)), x2(get_xpos(pos.x + rad));
		int const y1(get_ypos(pos.y - rad)), y2(get_ypos(pos.y + rad));

		// modify grass within radius of pos
		for (int y = y1; y <= y2; ++y) {
			for (int x = x1; x <= x2; ++x) {
				if (point_outside_mesh(x, y)) continue;
				assert(vbo > 0);
				unsigned const ix(y*MESH_X_SIZE + x);
				assert(ix+1 < mesh_to_grass_map.size());
				unsigned const start(mesh_to_grass_map[ix]), end(mesh_to_grass_map[ix+1]);
				assert(start <= end && end <= grass.size());
				if (start == end) continue; // no grass at this location
				unsigned min_up(end), max_up(start);

				for (unsigned i = start; i < end; ++i) {
					grass_t &g(grass[i]);
					float const dsq(p2p_dist_xy_sq(pos, g.p));
					if (dsq > rad_sq) continue; // too far away
					float const reld(sqrt(dsq)/rad), atten_val(1.0 - (1.0 - reld)*(1.0 - reld));
					bool updated(0);

					if (crush) {
						vector3d const sn(surface_normals[y][x]);
						float const length(g.dir.mag()), dx(g.p.x - pos.x), dy(g.p.y - pos.y);

						if (fabs(dot_product(g.dir/length, sn)) > 0.1) { // update if not flat against the mesh
							vector3d const new_dir(vector3d(dx, dy, -(sn.x*dx + sn.y*dy)/sn.z).get_norm()); // point away from crusing point

							if (dot_product(g.dir/length, new_dir) < 0.95) { // update if not already aligned
								g.dir   = (g.dir*(atten_val/length) + new_dir*(1.0 - atten_val)).get_norm()*length;
								g.n     = (g.n.get_norm()*atten_val + sn*(1.0 - atten_val)).get_norm()*g.n.mag();
								updated = 1;
							}
						}
					}
					if (burn) {
						UNROLL_3X(updated |= (g.c[i_] > 0);)
						if (updated) {UNROLL_3X(g.c[i_] = (unsigned char)(atten_val*g.c[i_]);)}
					}
					if (updated) {
						min_up = min(min_up, i);
						max_up = max(max_up, i);
					}
				} // for i
				if (min_up > max_up) continue; // nothing updated
				// FIXME: only upload modified once per frame?
				modified[ix] = 1;
				if (vbo_valid) upload_data_to_vbo(min_up, max_up+1);
				//data_valid = 0;
			} // for x
		} // for y
	}

	void heal_grass() {
		// WRITE: fix dir and color for grass that has been crushed or burned using modified
	}

	void upload_data() {
		if (empty()) return;
		RESET_TIME;
		// remove excess capacity from grass?
		upload_data_to_vbo(0, grass.size());
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
		point lpos;
		get_light_pos(lpos, light);

		if (light != last_light || lpos != last_lpos) {
			invalidate_shadows();
			last_light = light;
			last_lpos  = lpos;
		}
		check_for_updates();

		// draw the grass
		assert(vbo_valid && vbo > 0);
		bind_vbo(vbo);
		vert_norm_tc_color::set_vbo_arrays();
		//set_lighted_sides(2);
		select_texture(GRASS_BLADE_TEX);
		set_specular(0.1, 10.0);
		enable_blend();
		//glEnable(GL_POLYGON_SMOOTH);
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.75);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_NORMALIZE);
		glDrawArrays(GL_TRIANGLES, 0, 3*grass.size());
		glDisable(GL_COLOR_MATERIAL);
		glEnable(GL_NORMALIZE);
		disable_blend();
		set_specular(0.0, 1.0);
		glDisable(GL_ALPHA_TEST);
		glDisable(GL_TEXTURE_2D);
		//set_lighted_sides(1);
		bind_vbo(0);
	}
};

grass_manager_t grass_manager;


bool no_grass() {
	return (grass_density == 0 || !grass_enabled || snow_enabled() || vegetation == 0.0 || read_landscape);
}


void gen_grass(bool full_regen) {

	if (!full_regen) { // update shadows only
		grass_manager.invalidate_shadows();
		return;
	}
	bool const use_vbos(setup_gen_buffers_arb());
		
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
	if (!no_grass()) grass_manager.draw();
}

void modify_grass_at(point const &pos, float radius, bool crush, bool burn) {
	if (!no_grass()) grass_manager.modify_grass(pos, radius, crush, burn);
}



