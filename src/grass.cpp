// 3D World - Grass Generation and Rendering Code
// by Frank Gennari
// 9/28/10

#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"


bool has_grass(0);
pt_line_drawer grass_pld;


extern int island, default_ground_tex, read_landscape;
extern float vegetation, zmin, zmax, h_sand[], h_dirt[];
extern vector3d wind;


bool snow_enabled();


bool is_grass_shadowed(point const &pos, int light) {

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


void gen_grass() {

	has_grass = 0;
	grass_pld.clear();
	if (snow_enabled() || vegetation == 0.0 || read_landscape) return;
	int const light(get_light());
	float const *h_tex(island ? h_sand     : h_dirt);
	ttex  const *lttex(island ? lttex_sand : lttex_dirt);
	int   const   NTEX(island ? NTEX_SAND  : NTEX_DIRT);
	float const dz_inv(1.0/(zmax - zmin));
	float const grass_length(0.008);
	unsigned const num_grass(100);
	unsigned grass_count(0);
	
	for (int y = 0; y < MESH_Y_SIZE-1; ++y) {
		for (int x = 0; x < MESH_X_SIZE-1; ++x) {
			if (is_mesh_disabled(x, y)) continue; // mesh not drawn
			if (mesh_height[y][x] < water_matrix[y][x]) continue; // underwater
			float const xval(get_xval(x)), yval(get_yval(y));

			for (unsigned n = 0; n < num_grass; ++n) {
				float const xv(rand_uniform(xval, xval + DX_VAL));
				float const yv(rand_uniform(yval, yval + DY_VAL));
				float const mh(interpolate_mesh_zval(xv, yv, 0.0, 0, 1));

				if (default_ground_tex >= 0 || zmax == zmin) {
					if (default_ground_tex >= 0 && default_ground_tex != GROUND_TEX) continue;
				}
				else {
					// look for dominant green component in generated landscape texture instead?
					//float const relh(get_rel_height(mh, zmin, zmax));
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
				vector3d const dir((plus_z + signed_rand_vector(0.25) + wind*0.25).get_norm());
				point const p1(xv, yv, mh), p2(p1 + dir*grass_length);
				vector3d const normal(is_grass_shadowed(p1, light) ? zero_vector : dir);
				// FIXME: use textured quad(s) with better normal instead of line
				colorRGBA const grass_color(rand_uniform(0.1, 0.35), rand_uniform(0.5, 0.75), rand_uniform(0.0, 0.1), 1.0);
				grass_pld.add_line(p1, normal, grass_color, p2, normal, grass_color);
				++grass_count;
				has_grass = 1;
			}
		}
	}
	cout << "grass: " << grass_count << " out of " << XY_MULT_SIZE*num_grass << endl;
}


// called when light source moves, and to regen VBO(s)
void gen_grass_draw_data() {

	// write
}


void draw_grass() {

	if (!has_grass || snow_enabled()) return;
	static int last_light(-1);
	static point last_lpos(all_zeros);
	int const light(get_light());
	point lpos;
	get_light_pos(lpos, light);

	if (light != last_light || lpos != last_lpos) {
		gen_grass_draw_data();
		last_light = light;
		last_lpos  = lpos;
	}
	glDisable(GL_NORMALIZE);
	// FIXME: use VBO
	// FIXME: update shadows when light moves
	grass_pld.draw();
	glEnable(GL_NORMALIZE);
}



