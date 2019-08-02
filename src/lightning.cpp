// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/13/02

#include "3DWorld.h"
#include "mesh.h"
#include "openal_wrap.h"
#include "shaders.h"
#include "physics_objects.h" // for lightning_t


int    const FORK_PROB            = 50; // tuned for 64 z stacks
float  const FORK_ATTEN           = 0.8;
float  const BASE_L_STRENGTH_MULT = 150.0;
float  const BASE_L_DAMAGE_MULT   = 80.0;
float  const LIGHTNING_PERIOD     = 100.0; // average time between strikes in ticks
int    const DISCHARGE_RAD        = 20;
double const CHARGE_HALF_D        = 5.0;
float  const EVAP_AMOUNT          = 10.0;
float  const LITNING_ICE_EFF      = 0.5;
unsigned const MAX_LITN_FORKS     = 8;
bool const NEW_GEN_MODE           = 1;


int vmatrix_valid(0);
float L_STRENGTH_MULT, L_DAMAGE_MULT;
point litning_pos;
lightning_t l_strike;

extern int is_cloudy, game_mode, world_mode, iticks, DISABLE_WATER, animate2;
extern float zmin, temperature, lt_green_int, water_plane_z, ztop, fticks;
extern vector<valley> valleys;


void do_lightning_damage(point &pos, float damage, int hit_water);


inline float get_lit_h(int xpos, int ypos) {

	float h(h_collision_matrix[ypos][xpos]);
	if (!v_collision_matrix[ypos][xpos].cvals.empty()) {h = max(h, v_collision_matrix[ypos][xpos].zmax);}
	return h;
}


void compute_volume_matrix() {
	vmatrix_valid = 0; // just mark it as invalid
}


void compute_volume_matrix_if_invalid() {

	if (vmatrix_valid /*|| NEW_GEN_MODE*/) return;
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
		for (int j = 0; j < MESH_X_SIZE; ++j) {lh[yoff + j] = get_lit_h(j, i);}
	}
	for (int i = 1; i < MESH_Z_SIZE; ++i) { // calculate volume matrix shortest paths
		for (int j = 0; j < MESH_Y_SIZE; ++j) {
			int const yoff(j*MESH_X_SIZE);
			short **const vm = volume_matrix[i-1];

			for (int k = 0; k < MESH_X_SIZE; ++k) {
				if (lh[yoff + k] < zval) {
					minval = vm[j][k]+2;
					if (j > 0)             {min_eq(minval, short(vm[j-1][k] + 3));}
					if (j < MESH_Y_SIZE-1) {min_eq(minval, short(vm[j+1][k] + 3));}
					if (k > 0)             {min_eq(minval, short(vm[j][k-1] + 3));}
					if (k < MESH_X_SIZE-1) {min_eq(minval, short(vm[j][k+1] + 3));}
					volume_matrix[i][j][k] = minval;
				}
			}
		}
		zval += dz;
	}
}


void lightning_t::gen() {

	if (MESH_Z_SIZE == 0) return; // disabled
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
	int const rnum(int(LIGHTNING_PERIOD/fticks));

	if (rnum <= 1 || (rgen.rand()%rnum) != 0) {
		if (enabled == 1) {
			if (time < int(0.1f*LITNING_TIME/fticks)) {
				draw();
				if (!animate2) return;
				time += iticks;
			}
			else {enabled = 0;}
		}
		return;
	}
	compute_volume_matrix_if_invalid();
	enabled = 0;
	time    = 0;
	path.clear();
	int xpos(0), ypos(0);
	float max_e(0.0), charge(0);
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
		for (int j = 0; j < MESH_X_SIZE; ++j) {charge_dist[i][j] += d_charge;}
	}
	enabled = 1;
	cells_seen.clear();
	start.assign(get_xval(xpos), get_yval(ypos), (CLOUD_CEILING + ztop));
	if (NEW_GEN_MODE) {gen_recur_v2(start, max_e);}
	else {gen_recur(start, max_e, xpos, ypos, MESH_Z_SIZE-1, (CLOUD_CEILING + ztop));}
	litning_pos    = end;
	litning_pos.z += 0.01;
	draw();
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


void lightning_t::gen_recur(point const &start, float strength, int xpos, int ypos, int zpos, float zval) {

	int i(0), hit_water(0), full_path(1), lforks[5][2];
	unsigned const path_id((unsigned)path.size());
	if (path_id >= MAX_LITN_FORKS) return;
	if (zpos < 0 || zpos >= MESH_Z_SIZE || point_outside_mesh(xpos, ypos)) return; // note this will return if MESH_Z_SIZE == 1
	path.push_back(line3d());
	int const max_points(MESH_X_SIZE + MESH_Y_SIZE + MESH_Z_SIZE);
	vector<point> points(max_points);
	int ptx(xpos), pty(ypos);
	float const dz((CLOUD_CEILING + ztop - zmin)/MESH_Z_SIZE);
	points[0] = start;

	for (; zpos > 0 && i < max_points-2; --zpos, ++i) {
		if (i > 0) {points[i].assign(get_xval(ptx), get_yval(pty), zval);}
		int val(volume_matrix[zpos-1][pty][ptx]), nforks(1), nforks2(0), x0(0);
		lforks[0][0] = ptx;
		lforks[0][1] = pty;
		if (pty > 0)             {add_forks(lforks, nforks, val, zpos, ptx,   pty-1);}
		if (pty < MESH_Y_SIZE-1) {add_forks(lforks, nforks, val, zpos, ptx,   pty+1);}
		if (ptx > 0)             {add_forks(lforks, nforks, val, zpos, ptx-1, pty  );}
		if (ptx < MESH_X_SIZE-1) {add_forks(lforks, nforks, val, zpos, ptx+1, pty  );}

		for (int j = 0; j < nforks; ++j) {
			cell_ix_t const cix(get_cell_ix(lforks[j][0], lforks[j][1], zpos));

			if (cells_seen.insert(cix).second) { // if new	
				if (nforks2 < j) {
					for (unsigned d = 0; d < 2; ++d) {lforks[nforks2][d] = lforks[j][d];}
				}
				++nforks2;
			}
		}
		if (nforks2 == 0) {
			x0  = ((unsigned)rgen.rand())%nforks;
			ptx = lforks[x0][0];
			pty = lforks[x0][1];
			if (path_id != 0) {full_path = 0;} // truncate/merge non-main path fork
			break;
		}
		x0  = ((unsigned)rgen.rand())%nforks2;
		ptx = lforks[x0][0];
		pty = lforks[x0][1];

		for (int j = 0; j < nforks2; ++j) {
			if (j != x0 && i >= 2 && (rgen.rand()%max(1, FORK_PROB)) == 0) {
				gen_recur(points[i], FORK_ATTEN*strength, lforks[j][0], lforks[j][1], zpos, zval);
			}
		}
		if (val == 0) break;
		if (zval <= water_matrix[pty][ptx]) {hit_water = 1; break;}
		if (zval <= get_lit_h(ptx, pty))  break;
		if (rgen.rand()%MESH_Z_SIZE == 0) break;
		zval -= dz;
	} // for i
	assert(path_id < path.size());

	if (full_path) { // add lightning path end point
		points[++i].assign(get_xval(ptx), get_yval(pty), (hit_water ? water_matrix[pty][ptx] : get_lit_h(ptx, pty)));
	}
	else if (i < 2) { // small segment, discard
		path.pop_back();
		return;
	}
	points.resize(i+full_path);
	add_path(points, path_id, strength, full_path, hit_water);
}


void lightning_t::gen_recur_v2(point const &start, float strength) {

	int hit_water(0), full_path(1);
	unsigned const path_id((unsigned)path.size());
	if (path_id >= MAX_LITN_FORKS) return;
	if (path_id > 0 && !is_over_mesh(start)) return; // bad fork - skip
	path.push_back(line3d());
	vector<point> points;
	unsigned const max_steps(MESH_X_SIZE + MESH_Y_SIZE);
	float const step_sz(2.0*HALF_DXY);
	vector3d delta(0.0, 0.0, step_sz);
	point pos(start);

	for (unsigned step = 0; step < max_steps; ++step) {
		if (step > 1 && (rgen.rand()%max(1, 3*FORK_PROB/2)) == 0) {gen_recur_v2(pos, FORK_ATTEN*strength);}
		float depth(0.0);
		
		if (is_underwater(pos, 0, &depth)) { // terminate path - hit water
			pos.z += depth;
			points.push_back(pos);
			hit_water = 1;
			break;
		}
		int xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

		if (point_outside_mesh(xpos, ypos)) { // outside mesh
			if (path_id == 0) {xpos = max(0, min(MESH_X_SIZE-1, xpos)); ypos = max(0, min(MESH_Y_SIZE-1, ypos));} // primary path - apply clamp
			else {full_path = 0; break;} // terminate path
		}
		float const zval(get_lit_h(xpos, ypos));

		if (pos.z <= zval) { // terminate path - hit something
			if (!points.empty() && (points.back().z - zval) > 0.5*step_sz) { // hit top of an object, add a line segment down
				pos.z = zval;
				points.push_back(pos);
			} // else hit side of an object, end path here
			break;
		}
		points.push_back(pos);
		delta  += rgen.signed_rand_vector_spherical(0.75*step_sz); // add in random dir change
		delta  *= step_sz*rgen.rand_uniform(0.5, 1.0)/delta.mag(); // recompute length
		delta.z = -fabs(delta.z); // must be negative
		if (!is_over_mesh(pos + delta)) {delta.x *= -1.0; delta.y *= -1.0;} // if this delta takes us outside the mesh, invert the xy vals
		pos += delta;
	} // for step
	//cout << TXT(path_id) << TXT(points.size()) << TXT(hit_water) << TXT(full_path) << endl;
	assert(path_id < path.size());
	if (full_path) {points.push_back(pos);} // add lightning path end point
	if (points.size() < 2) {path.pop_back();} // small segment, discard
	else {add_path(points, path_id, strength, full_path, hit_water);}
}


void lightning_t::add_path(vector<point> &points, unsigned path_id, float strength, bool do_damage, bool hit_water) {

	float const dist_ratio(1.0/distance_to_camera(points.back()));
	path[path_id].color  = LITN_C;
	path[path_id].points = points;
	path[path_id].width  = max(rgen.rand_uniform(1.0, 2.0), min(4.0f, L_STRENGTH_MULT*strength*dist_ratio));

	if (path_id == 0) {
		end = points.back(); // main branch
		path[0].width *= 2.0;
	}
	if (do_damage) {do_lightning_damage(points.back(), L_DAMAGE_MULT*strength, hit_water);}
}


void do_lightning_damage(point &pos, float damage, int hit_water) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (point_outside_mesh(xpos, ypos))          return; // object off edge
	if (pos.z - mesh_height[ypos][xpos] > 0.05f) return; // didn't strike the mesh
	
	if (!DISABLE_WATER && wminside[ypos][xpos] == 1) { // interior water
		int const wsi(watershed_matrix[ypos][xpos].wsi);
		float ice_eff((temperature <= W_FREEZE_POINT) ? LITNING_ICE_EFF : 1.0);
		valleys[wsi].fvol += (MEL_EVAP_RATIO/SNOW_ACC)*ice_eff*accumulation_matrix[ypos][xpos]; // melt snow
		accumulation_matrix[ypos][xpos] = 0.0;

		if (hit_water) {
			draw_splash(pos.x, pos.y, pos.z, 0.05);
			valleys[wsi].w_volume = max(0.0f, (valleys[wsi].w_volume - ice_eff*EVAP_AMOUNT));
		}
		else {
			surface_damage[ypos][xpos] += damage;
		}
	}
	else if (!DISABLE_WATER && wminside[ypos][xpos] == 2) { // exterior water
		if (hit_water) {draw_splash(pos.x, pos.y, pos.z, 0.05);}
	}
	else {
		surface_damage[ypos][xpos] += damage;
		modify_grass_at(pos, HALF_DXY, 0, 1); // burn
	}
	if (ypos > 0)             {surface_damage[ypos-1][xpos  ] += damage;}
	if (xpos > 0)             {surface_damage[ypos  ][xpos-1] += damage;}
	if (xpos > 0 && ypos > 0) {surface_damage[ypos-1][xpos-1] += damage;}
	add_crater_to_landscape_texture(pos.x, pos.y, 0.5*damage);
}


void lightning_t::draw() const {

	float const lscale(LITNING_LINEAR_I);
	enable_blend();
	set_additive_blend_mode();
	shader_t s;
	s.begin_simple_textured_shader();

	for (unsigned i = 0; i < path.size(); ++i) {
		assert(!path[i].points.empty());
		if (animate2) {add_dynamic_light(0.6*path[i].width*lscale, path[i].points.back(), LITN_C);}
		path[i].draw_lines();
	}
	s.end_shader();
	set_std_blend_mode();
	disable_blend();
	if (animate2) {add_dynamic_light(7.8*lscale, litning_pos, LITN_C);}
}

