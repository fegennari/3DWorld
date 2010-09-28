// 3D World - Grass Generation and Rendering Code
// by Frank Gennari
// 9/28/10

#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"


float    const GRASS_LENGTH = 0.008;
unsigned const NUM_GRASS    = 100;


extern int island, default_ground_tex, read_landscape;
extern float vegetation, zmin, zmax, h_sand[], h_dirt[];
extern vector3d wind;


bool snow_enabled();


class grass_manager_t {
	
	vector<vert_norm_tc_color> data;
	unsigned vbo;
	bool vbo_valid;

public:
	grass_manager_t() : vbo(0), vbo_valid(0) {}
	~grass_manager_t() {clear();}
	size_t size() const {return (data.size() >> 1);} // 2 points per grass blade
	bool empty()  const {return data.empty();}
	void invalidate_vbo() {vbo_valid = 0;}
	
	void clear() {
		delete_vbo(vbo);
		vbo = 0;
		data.clear();
	}

	void add_grass(point const &p1) {
		vector3d const dir((plus_z + signed_rand_vector(0.25) + wind*0.25).get_norm()); // FIXME: make dynamic?
		point const p2(p1 + dir*(GRASS_LENGTH*rand_uniform(0.9, 1.1)));
		colorRGB const grass_color(rand_uniform(0.1, 0.35), rand_uniform(0.5, 0.75), rand_uniform(0.0, 0.1));
		// FIXME: use textured quad(s) with better normal instead of line
		data.push_back(vert_norm_tc_color(p1, plus_z, 0.0, 0.0, grass_color));
		data.push_back(vert_norm_tc_color(p2, plus_z, 0.1, 0.0, grass_color)); // second color different?
	}

	bool is_pt_shadowed(point const &pos) const {
		int const light(get_light());
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
		bool shadowed(0), unshadowed(0);

		for (int y = max(0, ypos-1); y <= min(MESH_Y_SIZE-1, ypos+1); ++y) { // test 3x3 window around the point
			for (int x = max(0, xpos-1); x <= min(MESH_X_SIZE-1, xpos+1); ++x) {
				(((shadow_mask[light][y][x] & SHADOWED_ALL) != 0) ? shadowed : unshadowed) = 1;
			}
		}
		if (shadowed != unshadowed) return shadowed; // only one was set, so mesh is in agreement
		return !is_visible_to_light_cobj(pos, light, 0.0, -1, 1); // neither (off the mesh) or both (conflict)
	}

	void gen_draw_data() {
		if (empty()) return;
		//cout << "mem used: " << data.size()*sizeof(vert_norm_tc_color) << endl;
		// remove excess capacity from data?

		for (vector<vert_norm_tc_color>::iterator i = data.begin(); i != data.end(); i += 2) {
			i->n = (i+1)->n = (is_pt_shadowed(i->v) ? zero_vector : ((i+1)->v - i->v).get_norm());
		}
		bool const use_vbos(setup_gen_buffers_arb());
		assert(use_vbos);
		delete_vbo(vbo);
		vbo = create_vbo();
		bind_vbo(vbo);
		upload_vbo_data(&data.front(), data.size()*sizeof(vert_norm_tc_color));
		bind_vbo(0);
		vbo_valid = 1;
	}

	void draw() {
		if (empty()) return;
		if (!vbo_valid) gen_draw_data();
		assert(vbo_valid && vbo > 0);
		bind_vbo(vbo);
		vert_norm_tc_color::set_vbo_arrays();
		//select_texture(GRASS_TEX);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_NORMALIZE);
		glDrawArrays(GL_LINES, 0, data.size());
		glDisable(GL_COLOR_MATERIAL);
		glEnable(GL_NORMALIZE);
		//glDisable(GL_TEXTURE_2D);
		bind_vbo(0);
	}
};

grass_manager_t grass_manager;


void gen_grass() {

	grass_manager.clear();
	if (NUM_GRASS == 0 || snow_enabled() || vegetation == 0.0 || read_landscape) return;
	float const *h_tex(island ? h_sand     : h_dirt);
	ttex  const *lttex(island ? lttex_sand : lttex_dirt);
	int   const   NTEX(island ? NTEX_SAND  : NTEX_DIRT);
	float const dz_inv(1.0/(zmax - zmin));
	
	for (int y = 0; y < MESH_Y_SIZE-1; ++y) {
		for (int x = 0; x < MESH_X_SIZE-1; ++x) {
			if (is_mesh_disabled(x, y)) continue; // mesh not drawn
			if (mesh_height[y][x] < water_matrix[y][x]) continue; // underwater
			float const xval(get_xval(x)), yval(get_yval(y));

			for (unsigned n = 0; n < NUM_GRASS; ++n) {
				float const xv(rand_uniform(xval, xval + DX_VAL));
				float const yv(rand_uniform(yval, yval + DY_VAL));
				float const mh(interpolate_mesh_zval(xv, yv, 0.0, 0, 1));

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
				grass_manager.add_grass(point(xv, yv, mh));
			}
		}
	}
	cout << "grass: " << grass_manager.size() << " out of " << XY_MULT_SIZE*NUM_GRASS << endl;
}


// called when light source moves, and to regen VBO(s)
void gen_grass_draw_data() {

	grass_manager.invalidate_vbo();
}


void draw_grass() {

	if (grass_manager.empty() || snow_enabled()) return;
	static int last_light(-1);
	static point last_lpos(all_zeros);
	int const light(get_light());
	point lpos;
	get_light_pos(lpos, light);

	if (light != last_light || lpos != last_lpos) {
		grass_manager.invalidate_vbo();
		gen_grass_draw_data();
		last_light = light;
		last_lpos  = lpos;
	}
	grass_manager.draw();
}



