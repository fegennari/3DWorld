// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/13/02

#include "3DWorld.h"
#include "mesh.h"
#include "openal_wrap.h"


int    const FORK_PROB            = 50;
float  const FORK_ATTEN           = 0.8;
float  const BASE_L_STRENGTH_MULT = 150.0;
float  const BASE_L_DAMAGE_MULT   = 80.0;
int    const L_FORK_END_PARAM     = 20;
float  const LIGHTNING_FREQ       = 100.0;
int    const DISCHARGE_RAD        = 20;
double const CHARGE_HALF_D        = 5.0;
float  const EVAP_AMOUNT          = 10.0;
float  const LITNING_ICE_EFF      = 0.5;
unsigned const MAX_LITN_FORKS     = 100;


int vmatrix_valid(0);
float L_STRENGTH_MULT, L_DAMAGE_MULT;
point litning_pos;
lightning l_strike;

extern int is_cloudy, game_mode, world_mode, ocean_set, iticks, DISABLE_WATER, animate2;
extern float zmin, temperature, lt_green_int, water_plane_z, ztop, fticks;
extern point ocean;
extern vector<valley> valleys;


void do_lightning_damage(point &pos, float damage, int hit_water);



inline float get_lit_h(int xpos, int ypos) {

	float h(h_collision_matrix[ypos][xpos]);
	if (!v_collision_matrix[ypos][xpos].cvals.empty()) h = max(h, v_collision_matrix[ypos][xpos].zmax);
	return h;
}


void compute_volume_matrix() {

	vmatrix_valid = 0; // just mark it as invalid
}


void compute_volume_matrix_if_invalid() {

	if (vmatrix_valid) return;
	short minval(0);
	float const dz((CLOUD_CEILING + ztop - zmin)/MESH_Z_SIZE);
	float zval(zmin + dz);
	vector<float> lh(XY_MULT_SIZE);
	vmatrix_valid = 1;

	for (int i = 0; i < MESH_Z_SIZE; ++i) { // clear volume matrix
		matrix_clear_2d(volume_matrix[i]);
	}
	for (int i = 0; i < MESH_Y_SIZE; ++i) { // compute lightning strike heights
		int const yoff(i*MESH_X_SIZE);
		
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			lh[yoff + j] = get_lit_h(j, i);
		}
	}
	for (int i = 1; i < MESH_Z_SIZE; ++i) { // calculate volume matrix shortest paths
		for (int j = 0; j < MESH_Y_SIZE; ++j) {
			int const yoff(j*MESH_X_SIZE);
			short **const vm = volume_matrix[i-1];

			for (int k = 0; k < MESH_X_SIZE; ++k) {
				if (lh[yoff + k] < zval) {
					minval = vm[j][k]+2;
					if (j > 0)             minval = min(minval, short(vm[j-1][k] + 3));
					if (j < MESH_Y_SIZE-1) minval = min(minval, short(vm[j+1][k] + 3));
					if (k > 0)             minval = min(minval, short(vm[j][k-1] + 3));
					if (k < MESH_X_SIZE-1) minval = min(minval, short(vm[j][k+1] + 3));
					volume_matrix[i][j][k] = minval;
				}
			}
		}
		zval += dz;
	}
}


void lightning::gen() {

	L_STRENGTH_MULT = BASE_L_STRENGTH_MULT*(MESH_Z_SIZE/50.0);
	L_DAMAGE_MULT   = BASE_L_DAMAGE_MULT  *(MESH_Z_SIZE/50.0);
	is_cloudy       = 1;

	if (enabled == -1) {
		enabled = 0;
		path.clear();
		return;
	}
	if (!animate2) {
		draw();
		return;
	}
	int const rnum(int(fticks*LIGHTNING_FREQ));

	if (rnum <= 1 || (rand()%rnum) != 0) {
		if (enabled == 1) {
			if (time < int(fticks*((float)LITNING_TIME))) {
				draw();
				if (!animate2) return;
				time += iticks;
			}
			else {
				enabled = 0;
			}
		}
		return;
	}
	compute_volume_matrix_if_invalid();
	enabled = 0;
	time    = 0;

	for (unsigned i = 0; i < path.size(); ++i) {
		path[i].destroy();
	}
	path.clear();
	int xpos(0), ypos(0);
	float max_e(0.0), charge(0);
	point pt(0, 0, 0);
	static int l_frame_counter(1);
	int const zpos(MESH_Z_SIZE - 1);

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			float const this_e(charge_dist[i][j]/((float)volume_matrix[zpos][i][j]));

			if (this_e > max_e) {
				max_e = this_e;
				xpos  = j;
				ypos  = i;
			}
		}
	}
	int const x0(max(0, (xpos - DISCHARGE_RAD))), x1(min(MESH_X_SIZE-1, (xpos + DISCHARGE_RAD)));
	int const y0(max(0, (ypos - DISCHARGE_RAD))), y1(min(MESH_Y_SIZE-1, (ypos + DISCHARGE_RAD)));

	for (int i = y0; i <= y1; ++i) {
		for (int j = x0; j <= x1; ++j) {
			float const distance(sqrt(float((xpos - j)*(xpos - j) + (ypos - i)*(ypos - i))));
			float const d_charge(charge_dist[i][j]/pow(2.0, double(distance/CHARGE_HALF_D)));
			charge            += d_charge;
			charge_dist[i][j] -= d_charge;
		}
	}
	float const d_charge(charge/XY_MULT_SIZE);

	for (int i = 0; i < MESH_Y_SIZE; ++i) { // redistribute charge
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			charge_dist[i][j] += d_charge;
		}
	}
	enabled = 1;
	cells_seen.clear();
	start.assign(get_xval(xpos), get_yval(ypos), (CLOUD_CEILING + ztop));
	gen_recur(pt, max_e, xpos, ypos, MESH_Z_SIZE-1, (CLOUD_CEILING + ztop), l_frame_counter);
	litning_pos    = end;
	litning_pos.z += 0.01;
	draw();
	++l_frame_counter;
	play_thunder(litning_pos, 4.0, 0.0);
}


void add_forks(int lforks[5][2], int &nforks, int val, int zpos, int x, int y) {

	int const val2(volume_matrix[zpos-1][y][x]);

	if (val2 < val) {
		nforks = 1;
		val    = val2;
		lforks[0][0] = x;
		lforks[0][1] = y;
	}
	else if (val2 == val) {
		lforks[nforks][0] = x;
		lforks[nforks][1] = y;
		++nforks;
	}
}


void lightning::gen_recur(point &start, float strength, int xpos, int ypos, int zpos, float zval, int l_frame_counter) {

	int i(0), hit_water(0), lforks[5][2];
	unsigned const path_id((unsigned)path.size());
	if (path_id >= MAX_LITN_FORKS) return;
	if (zpos < 0 || zpos >= MESH_Z_SIZE || point_outside_mesh(xpos, ypos)) return; // note this will return if MESH_Z_SIZE == 1
	path.push_back(line3d());
	int const max_points(MESH_X_SIZE + MESH_Y_SIZE + MESH_Z_SIZE);
	vector<point> points(max_points);
	int ptx(xpos), pty(ypos);
	float const dz((CLOUD_CEILING + ztop - zmin)/MESH_Z_SIZE);

	if (path_id > 0) {
		points[0] = start;
		++i;
	}
	for (; zpos > 0 && i < max_points-2; --zpos, ++i) {
		points[i].assign(get_xval(ptx), get_yval(pty), zval);
		int val(volume_matrix[zpos-1][pty][ptx]), val2(0), nforks(1), nforks2(0), x0(0);
		lforks[0][0] = ptx;
		lforks[0][1] = pty;
		if (pty > 0)             add_forks(lforks, nforks, val, zpos, ptx,   pty-1);
		if (pty < MESH_Y_SIZE-1) add_forks(lforks, nforks, val, zpos, ptx,   pty+1);
		if (ptx > 0)             add_forks(lforks, nforks, val, zpos, ptx-1, pty  );
		if (ptx < MESH_X_SIZE-1) add_forks(lforks, nforks, val, zpos, ptx+1, pty  );

		for (int j = 0; j < nforks; ++j) {
			cell_ix_t const cix(get_cell_ix(lforks[j][0], lforks[j][1], zpos));

			if (cells_seen.find(cix) == cells_seen.end()) {
				cells_seen.insert(cix);
				
				if (nforks2 < j) {
					for (unsigned d = 0; d < 2; ++d) lforks[nforks2][d] = lforks[j][d];
				}
				++nforks2;
			}
		}
		if (nforks2 == 0) {
			x0  = ((unsigned)rand())%nforks;
			ptx = lforks[x0][0];
			pty = lforks[x0][1];
			break;
		}
		x0  = ((unsigned)rand())%nforks2;
		ptx = lforks[x0][0];
		pty = lforks[x0][1];

		for (int j = 0; j < nforks2; ++j) {
			if (j != x0 && rand()%max(1, FORK_PROB) == 0) {
				gen_recur(points[i], FORK_ATTEN*strength, lforks[j][0], lforks[j][1], zpos, zval, l_frame_counter);
			}
		}
		if (val == 0) break;
		
		if (zval <= water_matrix[pty][ptx]) {
			hit_water = 1;
			break;
		}
		if (zval <= get_lit_h(ptx, pty))  break;
		if (ocean_set && zval < ocean.z)  break;
		if (rand()%L_FORK_END_PARAM == 0) break;
		zval -= dz;
	}
	++i;
	assert(i < max_points && path_id < path.size());
	points[i].assign(get_xval(ptx), get_yval(pty), (hit_water ? water_matrix[pty][ptx] : get_lit_h(ptx, pty)));
	float const dist_ratio(1.0/distance_to_camera(points[i]));
	points.resize(i+1);
	path[path_id].points = points;
	path[path_id].color  = LITN_C;
	path[path_id].width  = max(rand_uniform(1.0, 2.0), min(4.0f, L_STRENGTH_MULT*strength*dist_ratio));

	if (path_id == 0) {
		end = points[i]; // main branch
		path[0].width *= 2.0;
	}
	do_lightning_damage(points[i], L_DAMAGE_MULT*strength, hit_water);
}


void do_lightning_damage(point &pos, float damage, int hit_water) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos))         return; // object off edge
	if (pos.z - mesh_height[ypos][xpos] > 0.05) return; // didn't strike the mesh
	
	if (!DISABLE_WATER && wminside[ypos][xpos] == 1) { // interior water
		int const wsi(watershed_matrix[ypos][xpos].wsi);
		float ice_eff((temperature <= W_FREEZE_POINT) ? LITNING_ICE_EFF : 1.0);
		valleys[wsi].fvol += (MEL_EVAP_RATIO/SNOW_ACC)*ice_eff*accumulation_matrix[ypos][xpos]; // melt snow
		accumulation_matrix[ypos][xpos] = 0.0;

		if (hit_water) {
			draw_splash(pos.x, pos.y, (pos.z + 0.0001), 0.05);
			valleys[wsi].w_volume = max(0.0f, (valleys[wsi].w_volume - ice_eff*EVAP_AMOUNT));
		}
		else {
			surface_damage[ypos][xpos] += damage;
		}
	}
	else if (!DISABLE_WATER && wminside[ypos][xpos] == 2) { // exterior water
		if (hit_water) draw_splash(pos.x, pos.y, (pos.z + 0.0001), 0.05);
	}
	else {
		surface_damage[ypos][xpos] += damage;
		modify_grass_at(pos, HALF_DXY, 0, 1, 0, 0); // burn
	}
	if (ypos > 0)             surface_damage[ypos-1][xpos  ] += damage;
	if (xpos > 0)             surface_damage[ypos  ][xpos-1] += damage;
	if (xpos > 0 && ypos > 0) surface_damage[ypos-1][xpos-1] += damage;
	add_crater_to_landscape_texture(pos.x, pos.y, 0.5*damage);
}


void lightning::draw() const {

	float const lscale(LITNING_LINEAR_I);
	bool const smooth_lightning(1);
	if (smooth_lightning) {glEnable(GL_LINE_SMOOTH); enable_blend();}

	for (unsigned i = 0; i < path.size(); ++i) {
		assert(!path[i].points.empty());
		if (animate2) add_dynamic_light(0.6*path[i].width*lscale, path[i].points.back(), LITN_C);
		path[i].draw();
	}
	if (animate2) add_dynamic_light(7.8*lscale, litning_pos, LITN_C);
	if (smooth_lightning) {disable_blend(); glDisable(GL_LINE_SMOOTH);}
}




