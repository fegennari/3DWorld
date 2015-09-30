// 3D World - Precipitation Physics and Rendering (Rain and Snow)
// by Frank Gennari
// 10/24/12
#include "3DWorld.h"
#include "physics_objects.h"
#include "draw_utils.h"
#include "shaders.h"
#include "mesh.h"


float const TT_PRECIP_DIST = 20.0;

extern int animate2, display_mode;
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
	size_t get_num_precip() {return obj_groups[coll_id[PRECIP]].max_objects();}
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
	void check_pos(point &pos, point const &bot_pos) {
		if (0) {}
		else if (pos == all_zeros) {pos = gen_pt(rgen.rand_uniform(cur_zmin, cur_zmax));} // initial location
		else if (pos.z < cur_zmin) {pos = gen_pt(cur_zmax);} // start again near the top
		else if (!in_range(pos)  ) {pos = gen_pt(pos.z);} // move inside the range
		else if (world_mode == WMODE_GROUND) { // check bottom of raindrop below the mesh or top surface cobjs
			int const x(get_xpos(bot_pos.x)), y(get_ypos(bot_pos.y));
			if (!point_outside_mesh(x, y) && (bot_pos.z < mesh_height[y][x] || bot_pos.z < v_collision_matrix[y][x].zmax)) {pos = gen_pt(cur_zmax);} // start again near the top
		}
	}
	void check_size() {verts.resize(VERTS_PER_PRIM*get_num_precip(), all_zeros);}
};


class rain_manager_t : public precip_manager_t<2> {
public:
	void update() {
		//RESET_TIME;
		pre_update();
		vector3d const v(get_velocity(-0.2)), vinc(v*(0.1/verts.size())), dir(0.1*v.get_norm()); // length is 0.1
		vector3d vcur(v);

		//#pragma omp parallel for schedule(static,1)
		for (unsigned i = 0; i < verts.size(); i += 2) { // iterate in pairs
			check_pos(verts[i].v, verts[i+1].v);
			if (animate2) {verts[i].v += vcur; vcur += vinc;}
			verts[i+1].v = verts[i].v + dir;
		}
		//PRINT_TIME("Rain Update");
	}
	void render() const { // partially transparent
		if (empty()) return;
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
			if (dist_less_than(verts[i].v, camera, 0.5) && camera_pdu.point_visible_test(verts[i].v)) {
				drawer.add_line_as_tris(verts[i].v, verts[i+1].v, width, width, color, color);
			}
		}
		drawer.draw(); // draw nearby raindrops as triangles
		s.end_shader();
		disable_blend();
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

	if (!is_precip_enabled()) return;

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

