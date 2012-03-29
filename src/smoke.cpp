// 3D World - dynamic volumetric smoke effects code
// by Frank Gennari
// 4/30/11
#include "3DWorld.h"
#include "mesh.h"
#include "lightmap.h"
#include "gl_ext_arb.h"


bool const DYNAMIC_SMOKE     = 1; // looks cool
int const SMOKE_SKIPVAL      = 6;
int const SMOKE_SEND_SKIP    = 8;
int const INDIR_LT_SEND_SKIP = 12;

float const SMOKE_DENSITY    = 1.0;
float const SMOKE_MAX_CELL   = 0.125;
float const SMOKE_MAX_VAL    = 100.0;
float const SMOKE_DIS_XY     = 0.05;
float const SMOKE_DIS_ZU     = 0.08;
float const SMOKE_DIS_ZD     = 0.03;


bool smoke_visible(0), smoke_exists(0), have_indir_smoke_tex(0);
unsigned smoke_tid(0), last_smoke_update(0);
colorRGB const_indir_color(BLACK);
cube_t cur_smoke_bb;
vector<unsigned char> smoke_tex_data; // several MB

extern bool disable_shaders, no_smoke_over_mesh, indir_lighting_updated, no_sun_lpos_update, create_voxel_landscape;
extern int animate2, display_mode;
extern float czmin0;
extern colorRGBA cur_ambient, cur_diffuse;
extern lmap_manager_t lmap_manager;



struct smoke_manager {
	bool enabled, smoke_vis;
	float tot_smoke;
	cube_t bbox;

	smoke_manager() {reset();}

	bool is_smoke_visible(point const &pos) const {
		return sphere_in_camera_view(pos, HALF_DXY, 0); // could use cube_visible()
	}
	void reset() {
		for (unsigned i = 0; i < 3; ++i) { // set backwards so that nothing intersects
			bbox.d[i][0] =  SCENE_SIZE[i];
			bbox.d[i][1] = -SCENE_SIZE[i];
		}
		tot_smoke = 0.0;
		enabled   = 0;
		smoke_vis = 0;
	}
	void add_smoke(int x, int y, int z, float smoke_amt) {
		point const pos(get_xval(x), get_yval(y), get_zval(z));

		if (is_smoke_visible(pos)) {
			bbox.union_with_pt(pos);
			cur_smoke_bb.union_with_pt(pos);
			smoke_vis = 1;
		}
		tot_smoke += smoke_amt;
		enabled    = 1;
	}
	void adj_bbox() {
		for (unsigned i = 0; i < 3; ++i) {
			float const dval(SCENE_SIZE[i]/MESH_SIZE[i]);
			bbox.d[i][0] -= dval;
			bbox.d[i][1] += dval;
		}
	}
};

smoke_manager smoke_man, next_smoke_man;


inline void adjust_smoke_val(float &val, float delta) {

	val = max(0.0f, min(SMOKE_MAX_VAL, (val + delta)));
}


void add_smoke(point const &pos, float val) {

	if (!DYNAMIC_SMOKE || (display_mode & 0x80) || val == 0.0 || pos.z >= czmax) return;
	lmcell *const lmc(lmap_manager.get_lmcell(pos));
	if (!lmc) return;
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos) || pos.z >= v_collision_matrix[ypos][xpos].zmax || pos.z < mesh_height[ypos][xpos]) return; // above all cobjs/outside
	if (no_smoke_over_mesh && !is_mesh_disabled(xpos, ypos)) return;
	//if (!check_coll_line(pos, point(pos.x, pos.y, czmax), cindex, -1, 1, 0)) return;
	adjust_smoke_val(lmc->smoke, SMOKE_DENSITY*val);
	if (smoke_man.is_smoke_visible(pos)) smoke_exists = 1;
}


void diffuse_smoke(int x, int y, int z, lmcell &adj, float pos_rate, float neg_rate, int dim, int dir) {

	float delta(0.0); // Note: not using fticks due to instability
	
	if (lmap_manager.is_valid_cell(x, y, z)) {
		lmcell &lmc(lmap_manager.vlmap[y][x][z]);
		unsigned char const flow(dir ? adj.pflow[dim] : lmc.pflow[dim]);
		if (flow == 0) return;
		float const cur_smoke(lmc.smoke);
		delta  = (flow/255.0)*(adj.smoke - cur_smoke); // diffusion out of current cell and into cell xyz (can be negative)
		delta *= ((delta < 0.0) ? neg_rate : pos_rate);
		adjust_smoke_val(lmc.smoke, delta);
		delta  = (lmc.smoke - cur_smoke); // actual change
	}
	else { // edge cell has infinite smoke capacity and zero total smoke
		delta = 0.5*(pos_rate + neg_rate);
	}
	adjust_smoke_val(adj.smoke, -delta);
}


void distribute_smoke_for_cell(int x, int y, int z) {

	if (!lmap_manager.is_valid_cell(x, y, z)) return;
	lmcell &lmc(lmap_manager.vlmap[y][x][z]);
	if (lmc.smoke == 0.0) return;
	if (lmc.smoke < 0.005f) {lmc.smoke = 0.0; return;}
	int const dx(rand()&1), dy(rand()&1); // randomize the processing order
	float const xy_rate(SMOKE_DIS_XY*SMOKE_SKIPVAL), z_rate[2] = {SMOKE_DIS_ZU, SMOKE_DIS_ZD};
	next_smoke_man.add_smoke(x, y, z, lmc.smoke);

	for (unsigned d = 0; d < 2; ++d) { // x/y
		diffuse_smoke(x+((d^dx) ? 1 : -1), y, z, lmc, xy_rate, xy_rate, 0, (d^dx));
		diffuse_smoke(x, y+((d^dy) ? 1 : -1), z, lmc, xy_rate, xy_rate, 1, (d^dy));
	}
	for (unsigned d = 0; d < 2; ++d) { // up, down
		diffuse_smoke(x, y, (z + (d ? 1 : -1)),  lmc, z_rate[!d], z_rate[d], 2, d);
	}
}


void distribute_smoke() { // called at most once per frame

	RESET_TIME;
	if (!DYNAMIC_SMOKE || !smoke_exists || !animate2) return;
	assert(SMOKE_SKIPVAL > 0);
	static int cur_skip(0);
	
	if (cur_skip == 0) {
		//cout << "tot_smoke: " << smoke_man.tot_smoke << ", enabled: " << smoke_exists << ", visible: " << smoke_visible << endl;
		smoke_man     = next_smoke_man;
		smoke_man.adj_bbox();
		smoke_visible = smoke_man.smoke_vis;
		smoke_exists  = smoke_man.enabled;
		cur_smoke_bb  = smoke_man.bbox; //cube_t(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, min(zbottom, czmin), max(ztop, czmax));
		next_smoke_man.reset();
	}
	
	//#pragma omp parallel for schedule(static,1) // only helps when skipval is small
	for (int y = cur_skip; y < MESH_Y_SIZE; y += SMOKE_SKIPVAL) { // split the computation across several frames
		for (int x = 0; x < MESH_X_SIZE; ++x) {
			lmcell *vldata(lmap_manager.vlmap[y][x]);
			if (vldata == NULL) continue;
			
			for (int z = 0; z < MESH_SIZE[2]; ++z) {
				distribute_smoke_for_cell(x, y, z);
			}
		}
	}
	cur_skip = (cur_skip+1) % SMOKE_SKIPVAL;
	//PRINT_TIME("Distribute Smoke");
}


float get_smoke_at_pos(point const &pos) {

	if (!DYNAMIC_SMOKE  || !smoke_exists)  return 0.0;
	if (pos.z <= czmin0 || pos.z >= czmax) return 0.0;
	int const x(get_xpos(pos.x)), y(get_ypos(pos.y)), z(get_zpos(pos.z));
	if (point_outside_mesh(x, y) || z < 0 || z >= MESH_SIZE[2]) return 0.0;
	lmcell const *const vldata(lmap_manager.vlmap[y][x]);
	return ((vldata == NULL) ? 0.0 : vldata[z].smoke);
}


void reset_smoke_tex_data() {

	smoke_tex_data.clear();
}


void update_smoke_row(vector<unsigned char> &data, lmcell const &default_lmc, unsigned x_start, unsigned x_end, unsigned y, bool update_lighting) {

	unsigned const zsize(MESH_SIZE[2]), ncomp(4);
	float const smoke_scale(1.0/SMOKE_MAX_CELL);
	colorRGB default_color;
	default_lmc.get_final_color(default_color, 1.0);

	for (unsigned x = x_start; x < x_end; ++x) {
		lmcell const *const vlm(lmap_manager.vlmap[y][x]);
		if (vlm == NULL && !update_lighting) continue; // x/y pairs that get into here should also be constant
		unsigned const off(zsize*(y*MESH_X_SIZE + x));
		bool const check_z_thresh((display_mode & 0x01) && !is_mesh_disabled(x, y));
		float const mh(mesh_height[y][x]);

		for (unsigned z = 0; z < zsize; ++z) {
			unsigned const off2(ncomp*(off + z));
			float const smoke_val((vlm == NULL) ? 0.0 : CLIP_TO_01(smoke_scale*vlm[z].smoke));
			data[off2+3] = (unsigned char)(255*smoke_val); // alpha: smoke
			if (!update_lighting && !indir_lighting_updated) continue; // lighting not needed
				
			if (check_z_thresh && get_zval(z+1) < mh) { // adjust by one because GPU will interpolate the texel
				UNROLL_3X(data[off2+i_] = 0;)
			}
			else {
				float const indir_scale(create_voxel_landscape ? get_voxel_terrain_ao_lighting_val(get_xyz_pos(x, y, z)) : 1.0);
				colorRGB color;

				if (vlm == NULL) {
					color  = default_color;
					color *= indir_scale;
				}
				else {
					vlm[z].get_final_color(color, 1.0, indir_scale);
				}
				UNROLL_3X(data[off2+i_] = (unsigned char)(255*CLIP_TO_01(color[i_]));) // lmc.pflow[i_]
			}
		} // for z
	} // for x
}


void update_smoke_indir_tex_range(unsigned x_start, unsigned x_end, unsigned y_start, unsigned y_end, bool update_lighting) {

	assert(y_start < y_end && y_end <= (unsigned)MESH_Y_SIZE);
	if (smoke_tex_data.empty()) return; // not allocated
	unsigned const ncomp(4);
	lmcell default_lmc;
	default_lmc.set_outside_colors();
	default_lmc.get_final_color(const_indir_color, 1.0);
	
	if (indir_lighting_updated) { // running with multiple threads, don't use openmp
		for (int y = y_start; y < (int)y_end; ++y) { // split the computation across several frames
			update_smoke_row(smoke_tex_data, default_lmc, x_start, x_end, y, update_lighting);
		}
	}
	else { // use openmp here
		#pragma omp parallel for schedule(static,1)
		for (int y = y_start; y < (int)y_end; ++y) { // split the computation across several frames
			update_smoke_row(smoke_tex_data, default_lmc, x_start, x_end, y, update_lighting);
		}
	}
	if (smoke_tid == 0) { // create texture
		cout << "Allocating " << MESH_SIZE[2] << " by " << MESH_X_SIZE << " by " << MESH_Y_SIZE << " smoke texture of " << smoke_tex_data.size() << " bytes." << endl;
		smoke_tid = create_3d_texture(MESH_SIZE[2], MESH_X_SIZE, MESH_Y_SIZE, ncomp, smoke_tex_data, GL_LINEAR, GL_CLAMP_TO_EDGE);
	}
	else { // update region/sync texture
		unsigned const off(ncomp*y_start*MESH_X_SIZE*MESH_SIZE[2]);
		assert(off < smoke_tex_data.size());
		update_3d_texture(smoke_tid, 0, 0, y_start, MESH_SIZE[2], MESH_X_SIZE, (y_end - y_start), ncomp, &smoke_tex_data[off]);
	}
}


bool upload_smoke_indir_texture() {

	//RESET_TIME;
	assert((MESH_Y_SIZE%SMOKE_SEND_SKIP) == 0);
	// ok when texture z size is not a power of 2
	unsigned const sz(MESH_X_SIZE*MESH_Y_SIZE*MESH_SIZE[2]), ncomp(4);

	if (disable_shaders || lmap_manager.vlmap == NULL || sz == 0) {
		have_indir_smoke_tex = 0;
		return 0;
	}
	if (smoke_tex_data.empty()) {
		free_texture(smoke_tid);
		smoke_tex_data.resize(ncomp*sz, 0);
	}
	else { // will recreate the texture
		assert(smoke_tex_data.size() == ncomp*sz); // sz should be constant (per config file/3DWorld session)
	}
	static colorRGBA last_cur_ambient(ALPHA0), last_cur_diffuse(ALPHA0);
	bool const lighting_changed(cur_ambient != last_cur_ambient || cur_diffuse != last_cur_diffuse);
	bool const full_update(smoke_tid == 0 || (!no_sun_lpos_update && lighting_changed));
	bool const could_have_smoke(smoke_exists || last_smoke_update > 0);
	if (!full_update && !could_have_smoke && !indir_lighting_updated && !lighting_changed) return 0; // return 1?

	if (full_update) {
		last_cur_ambient = cur_ambient;
		last_cur_diffuse = cur_diffuse;
	}
	if (smoke_exists) {
		last_smoke_update = SMOKE_SEND_SKIP;
	}
	else if (last_smoke_update > 0) {
		--last_smoke_update;
	}
	static int cur_block(0);
	unsigned const skipval(could_have_smoke ? SMOKE_SEND_SKIP : INDIR_LT_SEND_SKIP);
	unsigned const block_size(MESH_Y_SIZE/skipval);
	unsigned const y_start(full_update ? 0           :  cur_block*block_size);
	unsigned const y_end  (full_update ? MESH_Y_SIZE : (y_start + block_size));
	update_smoke_indir_tex_range(0, MESH_X_SIZE, y_start, y_end, full_update);
	cur_block = (full_update ? 0 : (cur_block+1) % skipval);
	if (cur_block == 0) indir_lighting_updated = 0; // only stop updating after we wrap around to the beginning again
	have_indir_smoke_tex = 1;
	//PRINT_TIME("Smoke Upload");
	return 1;
}

