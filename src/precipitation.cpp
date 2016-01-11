// 3D World - Precipitation Physics and Rendering (Rain and Snow)
// by Frank Gennari
// 10/24/12
#include "3DWorld.h"
#include "physics_objects.h"
#include "draw_utils.h"
#include "shaders.h"
#include "mesh.h"


float const TT_PRECIP_DIST = 20.0;

extern int animate2, begin_motion, display_mode, camera_coll_id, precip_mode, DISABLE_WATER;
extern float temperature, fticks, zmin, water_plane_z, brightness, XY_SCENE_SIZE;
extern vector3d wind;
extern int coll_id[];
extern obj_group obj_groups[];


template <unsigned VERTS_PER_PRIM> class precip_manager_t {
protected:
	typedef vert_wrap_t vert_type_t;
	vector<vert_type_t> verts;
	rand_gen_t rgen;
	float prev_zmin, cur_zmin, prev_zmax, cur_zmax;

public:
	precip_manager_t() : prev_zmin(get_zmin()), cur_zmin(prev_zmin), prev_zmax(get_zmax()), cur_zmax(prev_zmax) {}
	void clear () {verts.clear();}
	bool empty () const {return verts.empty();}
	size_t size() const {return verts.size();}
	float get_zmin() const {return ((world_mode == WMODE_GROUND) ? zbottom : get_tiled_terrain_water_level());}
	float get_zmax() const {return get_cloud_zmax();}
	float get_precip_dist() const {return ((world_mode == WMODE_GROUND) ? XY_SCENE_SIZE : TT_PRECIP_DIST);}
	size_t get_num_precip() {return 700*get_precip_rate();} // similar to precip max objects
	bool in_range(point const &pos) const {return dist_xy_less_than(pos, get_camera_pos(), get_precip_dist());}
	vector3d get_velocity(float vz) const {return fticks*(0.02*wind + vector3d(0.0, 0.0, vz));}
	
	void pre_update() {
		cur_zmin = get_zmin();
		cur_zmax = get_zmax();

		// if zmin or zmax changes by more than some amount, then clear and regen point z-values so that rain/snow stays uniformly spaced in z
		if (fabs(prev_zmin - cur_zmin) > 0.05*(get_zmax() - cur_zmin) || fabs(prev_zmax - cur_zmax) > 0.25*(get_zmax() - cur_zmin)) {
			clear();
			prev_zmin = cur_zmin;
			prev_zmax = cur_zmax;
		}
		check_size();
	}
	point gen_pt(float zval) {
		point const camera(get_camera_pos());
		float const precip_dist(get_precip_dist());

		while (1) {
			vector3d const off(precip_dist*rgen.signed_rand_float(), precip_dist*rgen.signed_rand_float(), zval);
			if (off.x*off.x + off.y*off.y < precip_dist*precip_dist) {return (vector3d(camera.x, camera.y, 0.0) + off);}
		}
		return zero_vector; // never gets here
	}
	bool check_splash_dist(point const &pos) {
		point const camera(get_camera_pos());
		return (pos.z < camera.z && dist_less_than(camera, pos, 5.0)); // skip splashes above the camera (assuming the surface points up)
	}
	void maybe_add_rain_splash(point const &pos, point const &bot_pos, float z_int, deque<sphere_t> *splashes, int x, int y, bool in_water) {
		if (splashes == nullptr) return;
		float const t((z_int - pos.z)/(bot_pos.z - pos.z));
		point const cpos(pos + (bot_pos - pos)*t);
		if (!camera_pdu.point_visible_test(cpos)) return;
		if (check_splash_dist(cpos)) {splashes->push_back(sphere_t(cpos, 1.0));}
		if (in_water && (rgen.rand() & 1)) {add_splash(cpos, x, y, 0.5, 0.01, 0);} // 50% of the time
	}
	bool is_bot_pos_valid(point &pos, point const &bot_pos, deque<sphere_t> *splashes=nullptr) {
		if (world_mode != WMODE_GROUND) return 1;
		// check bottom of raindrop/snow below the mesh or top surface cobjs (even if just created)
		int const x(get_xpos(bot_pos.x)), y(get_ypos(bot_pos.y));
		if (point_outside_mesh(x, y)) return 1;
			
		if (!DISABLE_WATER && (display_mode & 0x04) && pos.z < water_matrix[y][x]) { // water collision
			if (rgen.rand() & 1) {maybe_add_rain_splash(pos, bot_pos, water_matrix[y][x], splashes, x, y, 1);} // 50% of the time
			return 0;
		}
		else if (pos.z < mesh_height[y][x]) { // mesh collision
			maybe_add_rain_splash(pos, bot_pos, mesh_height[y][x], splashes, x, y, 0); // line_intersect_mesh(pos, bot_pos, cpos);
			return 0;
		}
		else if (bot_pos.z < v_collision_matrix[y][x].zmax) { // possible cobj collision
			if (splashes != nullptr && check_splash_dist(bot_pos)) {
				point cpos;
				vector3d cnorm;
				int cindex;
				if (camera_pdu.point_visible_test(bot_pos) && check_coll_line_exact(pos, bot_pos, cpos, cnorm, cindex, 0.0, camera_coll_id)) {splashes->push_back(sphere_t(cpos, 1.0));}
			}
			return 0;
		}
		return 1;
	}
	void check_pos(point &pos, point const &bot_pos, deque<sphere_t> *splashes=nullptr) {
		if (pos == all_zeros) { // initial location
			vector3d const bot_delta(bot_pos - pos);
			for (unsigned attempt = 0; attempt < 16; ++attempt) { // make 16 attempts at choosing a valid starting z-value
				pos = gen_pt(rgen.rand_uniform(cur_zmin, cur_zmax));
				if (is_bot_pos_valid(pos, pos+bot_delta, nullptr)) break;
			}
		}
		else if (pos.z < cur_zmin)                          {pos = gen_pt(cur_zmax);} // start again near the top
		else if (!in_range(pos))                            {pos = gen_pt(pos.z   );} // move inside the range
		else if (!is_bot_pos_valid(pos, bot_pos, splashes)) {pos = gen_pt(cur_zmax);} // start again near the top
	}
	void check_size() {verts.resize(VERTS_PER_PRIM*get_num_precip(), all_zeros);}
};


class rain_manager_t : public precip_manager_t<2> {

	deque<sphere_t> splashes;

public:
	void update() {
		//RESET_TIME;
		pre_update();
		vector3d const v(get_velocity(-0.2)), vinc(v*(0.1/verts.size())), dir(0.1*v.get_norm()); // length is 0.1
		vector3d vcur(v);
		while (!splashes.empty() && splashes.front().radius > 4.0) {splashes.pop_front();} // remove old splashes from the front
		for (auto i = splashes.begin(); i != splashes.end(); ++i) {i->radius += 0.2*fticks;}

		//#pragma omp parallel for schedule(static,1) // not valid for splashes, and actually slower for light rain
		for (int i = 0; i < (int)verts.size(); i += 2) { // iterate in pairs
			check_pos(verts[i].v, verts[i+1].v, (begin_motion ? &splashes : nullptr));
			if (animate2) {verts[i].v += vcur; vcur += vinc;}
			verts[i+1].v = verts[i].v + dir;
		}
		//PRINT_TIME("Rain Update"); // 0.75ms for default rain intensity
	}
	void render() const { // partially transparent
		if (empty()) return;
		//RESET_TIME;
		colorRGBA color;
		get_avg_sky_color(color);
		color.alpha = 0.2;
		assert(!(size() & 1));
		enable_blend(); // split into point smooth and blend?
		shader_t s;
		s.begin_color_only_shader(color);
		draw_verts(verts, GL_LINES);
		line_tquad_draw_t drawer;
		point const camera(get_camera_pos());
		float const width = 0.002;

		for (unsigned i = 0; i < verts.size(); i += 2) { // iterate in pairs
			if (dist_less_than(verts[i].v, camera, 0.5) && (camera_pdu.point_visible_test(verts[i].v) || camera_pdu.point_visible_test(verts[i+1].v))) {
				drawer.add_line_as_tris(verts[i].v, verts[i+1].v, width, width, color, color);
			}
		}
		glDepthMask(GL_FALSE); // disable depth test
		drawer.draw(); // draw nearby raindrops as triangles
		s.end_shader();
		disable_blend();
		//PRINT_TIME("Rain Draw"); // similar to update time

		if (!splashes.empty()) {
			point const camera(get_camera_pos());
			float const size = 0.004; // 2x-8x rain line diameter
			ensure_filled_polygons();
			enable_blend();
			shader_t s;
			setup_smoke_shaders(s, 0.01, 0, 1, 1, 1, 1, 1, 0, 1);
			select_texture(FLARE2_TEX);
			quad_batch_draw qbd;
			
			for (auto i = splashes.begin(); i != splashes.end(); ++i) { // normal always faces up;
				float const sz(i->radius*size), alpha(0.75*(4.0 - i->radius)/3.0); // size increases with radius/time; alpha decreases with radius/time
				qbd.add_billboard(i->pos, camera, up_vector, colorRGBA(0.8, 0.9, 1.0, alpha), sz, sz, tex_range_t(), 0, &plus_z);
			}
			qbd.draw();
			s.end_shader();
			disable_blend();
			reset_fill_mode();
			//PRINT_TIME("Rain Draw + Splashes"); // 50% longer
		}
		glDepthMask(GL_TRUE);
	}
};


class snow_manager_t : public precip_manager_t<1> {
public:
	void update() {
		pre_update();
		vector3d const v(get_velocity(-0.02));
		float const vmult(0.1/verts.size());

		for (unsigned i = 0; i < verts.size(); ++i) {
			check_pos(verts[i].v, verts[i].v);
			if (animate2) {verts[i].v += (1.0 + i*vmult)*v;}
		}
	}
	void render() const {
		if (empty()) return;
		point_sprite_drawer psd;
		psd.reserve_pts(size());
		colorRGBA const color(WHITE*((world_mode == WMODE_GROUND) ? 1.5 : 1.0)*brightness); // constant
		for (vector<vert_type_t>::const_iterator i = verts.begin(); i != verts.end(); ++i) {psd.add_pt(vert_color(i->v, color));}
		psd.draw(WHITE_TEX, 1.0); // unblended pixels
	}
};


rain_manager_t rain_manager;
snow_manager_t snow_manager;


void draw_local_precipitation() {

	if (!(precip_mode & 1)) return;

	if (temperature <= W_FREEZE_POINT) { // draw snow
		rain_manager.clear();
		//if (world_mode != WMODE_INF_TERRAIN) return; // only drawn snow for tiled terrain?
		snow_manager.update();
		snow_manager.render();
	}
	else { // draw rain
		snow_manager.clear();
		rain_manager.update();
		rain_manager.render();
	}
}

