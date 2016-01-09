// 3D World - Object Group Drawing Code
// by Frank Gennari
// 4/27/12
#include "3DWorld.h"
#include "mesh.h"
#include "transform_obj.h"
#include "player_state.h"
#include "tree_leaf.h"
#include "textures_3dw.h"
#include "dynamic_particle.h"
#include "physics_objects.h"
#include "shaders.h"
#include "draw_utils.h"
#include "gl_ext_arb.h"


bool const DEBUG_COLORCODE   = 0;
bool const DEBUG_COLOR_COLLS = 0;
bool const SHOW_DRAW_TIME    = 0;
float const NDIV_SCALE       = 1.6;


// Global Variables
quad_batch_draw puddle_qbd;
pt_line_drawer obj_pld, snow_pld;


extern bool underwater, smoke_exists;
extern int display_mode, num_groups, teams, begin_motion, UNLIMITED_WEAPONS;
extern int window_width, window_height, game_mode, draw_model, animate2;
extern float fticks, TIMESTEP, base_gravity, brightness, indir_vert_offset, cobj_z_bias;
extern point star_pts[];
extern vector3d up_norm;
extern vector<spark_t> sparks;
extern obj_group obj_groups[];
extern obj_type object_types[];
extern player_state *sstates;
extern int coll_id[];


void draw_group(obj_group &objg, shader_t &s, lt_atten_manager_t &lt_atten_manager);
void draw_sized_point(dwobject &obj, float radius, float cd_scale, const colorRGBA &color, const colorRGBA &tcolor,
					  bool do_texture, shader_t &shader, int is_chunky=0);
void draw_ammo(obj_group &objg, float radius, const colorRGBA &color, int ndiv, int j, shader_t &shader, lt_atten_manager_t &lt_atten_manager);
void draw_smiley_part(point const &pos, point const &pos0, vector3d const &orient, int type,
					  int use_orient, int ndiv, shader_t &shader, float scale=1.0);
void draw_smiley(point const &pos, vector3d const &orient, float radius, int ndiv, int time,
				 float health, int id, mesh2d const *const mesh, shader_t &shader);
void draw_powerup(point const &pos, float radius, int ndiv, int type, const colorRGBA &color, shader_t &shader, lt_atten_manager_t &lt_atten_manager);
void draw_rolling_obj(point const &pos, point &lpos, float radius, int status, int ndiv, bool on_platform, int tid, xform_matrix *matrix, shader_t &shader);
void draw_skull(point const &pos, vector3d const &orient, float radius, int status, int ndiv, int time, shader_t &shader, bool burned);
void draw_rocket(point const &pos, vector3d const &orient, float radius, int type, int ndiv, int time, shader_t &shader);
void draw_seekd(point const &pos, vector3d const &orient, float radius, int type, int ndiv, shader_t &shader);
void draw_landmine(point pos, float radius, int ndiv, int time, int source, bool in_ammo, shader_t &shader);
void draw_plasma(point const &pos, point const &part_pos, float radius, float size, int ndiv, bool gen_parts, bool add_halo, int time, shader_t &shader);
void draw_chunk(point const &pos, float radius, vector3d const &v, vector3d const &vdeform, int charred, int ndiv, shader_t &shader);
void draw_grenade(point const &pos, vector3d const &orient, float radius, int ndiv, int time, bool in_ammo, bool is_cgrenade, shader_t &shader);
void draw_star(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate);
void draw_sawblade(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate, int ndiv, bool bloody);
void draw_shell_casing(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius,
					   float angle, float cd_scale, unsigned char type, shader_t &shader);
colorRGBA get_glow_color(dwobject const &obj, bool shrapnel_cscale);

int same_team(int source, int target); // gameplay
void draw_one_star(colorRGBA const &colorA, colorRGBA const &colorB, point const &pos, float radius, int ndiv, bool add_halo); // universe
void setup_shield_shader(shader_t &shader, int noise_tu_id); // ship.cpp



void set_obj_specular(unsigned flags, float specular_brightness, shader_t &shader) {

	if (flags & SPECULAR) {
		shader.set_specular(specular_brightness, 50.0);
	}
	else if (flags & LOW_SPECULAR) {
		shader.set_specular(0.4*specular_brightness, 10.0);
	}
}


void check_drawing_flags(unsigned flags, int init_draw, shader_t &shader) {

	if (init_draw) {
		set_obj_specular(flags, 0.5*brightness, shader);
		if (flags & BLEND) enable_blend();
	}
	else {
		if (flags & BLEND) disable_blend();
		if (flags & (SPECULAR | LOW_SPECULAR)) {shader.clear_specular();}
	}
}


void set_emissive_only(colorRGBA const &color, shader_t &shader) {

	shader.set_color_e(color);
	shader.set_cur_color(colorRGBA(BLACK, color.alpha));
}


void set_color_by_status(int status, shader_t &shader) {

	colorRGBA const colors[6] = {BLACK, RED, WHITE, YELLOW, BLUE, GRAY};
	assert(status >= 0 && status < 6);
	shader.set_cur_color(colors[status]);
}


void set_color_v2(const colorRGBA &color, int status, shader_t &shader) {
	if (DEBUG_COLORCODE) {set_color_by_status(status, shader);} else {shader.set_cur_color(color);}
}

inline bool is_droplet(int type) {
	return ((object_types[type].flags & OBJ_IS_DROP) || type == HAIL || type == CHARRED);
}

inline bool get_cull_face(int type, colorRGBA const &color) {
	return (color.alpha < 1.0 && type != ROCKET && type != STAR5 && type != GRENADE && type != LANDMINE);
}

void select_no_texture() {
	select_texture(WHITE_TEX);
}

void scale_color_uw(colorRGBA &color, point const &pos) {
	water_color_atten_at_pos(color, pos); // ???
	if (underwater) {color.R *= 0.45; color.G *= 0.45; color.B *= 0.85;}
}


// angle is in degrees
vector3d get_rotation_dirs_and_normal(vector3d const &o, float angle, vector3d &v1, vector3d &v2) {
	/*
	tXX  + c	tXY + sZ	tXZ - sY	0
	tXY - sZ	tYY + c		tYZ + sX	0
	tXZ + sY	tYZ - sX	tZZ + c		0
	0			0			0			1
		c = cos (theta), s = sin (theta), t = 1-cos (theta), and <X,Y,Z>
	x [-r, 0, 0, 1]^-1 => [tX^2 + c,  tXY + sZ, tXZ  - sY, 1]
	x [ 0, 0, q, 1]^-1 => [tXZ  + sY, tYZ - sX, tZ^2 + c,  1]
	*/
	angle *= TO_RADIANS;
	float const c(cos(angle)), s(sin(angle)), t(1.0 - c);
	v1 = vector3d((t*o.x*o.x + c),     (t*o.x*o.y - s*o.z), (t*o.x*o.z - s*o.y));
	v2 = vector3d((t*o.x*o.z + s*o.y), (t*o.y*o.z - s*o.x), (t*o.z*o.z + c));
	return cross_product(v2, v1).get_norm();
}


void add_rotated_triangle(point const &pos, vector3d const &o, float radius, float angle, colorRGBA const &color, vector<vert_norm_color> &verts) {

	point p1, p2;
	vector3d const normal(get_rotation_dirs_and_normal(o, angle, p1, p2));
	verts.push_back(vert_norm_color((pos + 1.5*radius*p1), normal, color));
	verts.push_back(vert_norm_color((pos - 1.5*radius*p1), normal, color));
	verts.push_back(vert_norm_color((pos + 3.0*radius*p2), normal, color));
}


void add_rotated_textured_triangle(point const &pos, vector3d const &o, float radius, float angle, float tscale, colorRGBA const &color, vector<vert_norm_tc_color> &verts) {

	point p1, p2;
	vector3d const normal(get_rotation_dirs_and_normal(o, angle, p1, p2));
	float const ts(123.456*radius), tt(654.321*radius); // pseudo-random
	verts.push_back(vert_norm_tc_color((pos + 1.5*radius*p1), normal, ts, tt, color));
	verts.push_back(vert_norm_tc_color((pos - 1.5*radius*p1), normal, ts+2*tscale*radius, tt, color));
	verts.push_back(vert_norm_tc_color((pos + 3.0*radius*p2), normal, ts, tt+2*tscale*radius, color));
}


void draw_polygon_side(point const *points, int npoints, vector3d const &normal, vector<vert_norm_color> &verts, colorRGBA const &color) {

	for (int i = 0; i < ((npoints == 3) ? 3 : 6); ++i) {verts.push_back(vert_norm_color(points[quad_to_tris_ixs[i]], normal, color));} // 1-2 triangles
}

void get_sorted_thick_poly_faces(point pts[2][4], pair<int, unsigned> faces[6], point const *points,
	unsigned npoints, vector3d const &norm, float thick, bool bfc);

void add_thick_triangle(point const &pos, vector3d const &o, float radius, float angle, float tscale,
	vector<vert_norm_color> &verts, float thickness, colorRGBA const &color)
{
	//tscale = tscale*radius
	unsigned const npoints(3);
	point p1, p2;
	vector3d const norm(get_rotation_dirs_and_normal(o, angle, p1, p2));
	point points[3] = {(pos + 1.5*radius*p1), (pos - 1.5*radius*p1), (pos + 3.0*radius*p2)};
	point pts[2][4];
	pair<int, unsigned> faces[6];
	get_sorted_thick_poly_faces(pts, faces, points, npoints, norm, thickness, 0);
	unsigned const nsides(unsigned(npoints)+2);
	
	for (unsigned fi = 0; fi < nsides; ++fi) { // draw back to front
		unsigned const s(faces[fi].second);

		if (s < 2) { // draw front and back
			if (!s) {std::reverse(pts[s], pts[s]+npoints);}
			draw_polygon_side(pts[s], npoints, (s ? norm : -norm), verts, color); // draw bottom surface
			if (!s) {std::reverse(pts[s], pts[s]+npoints);}
		}
		else { // draw sides
			unsigned const i(s-2), ii((i+1)%npoints);
			point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
			draw_polygon_side(side_pts, 4, get_poly_norm(side_pts), verts, color);
		}
	}
}


void draw_solid_object_groups() {

	draw_waypoints();
	draw_select_groups(1);
	if (display_mode & 0x0200) {d_part_sys.draw();}
}


void draw_transparent_object_groups() {
	draw_select_groups(0);
}


void draw_select_groups(int solid) {

	if (!begin_motion) return;
	shader_t s;
	s.set_prefix("#define USE_WINDING_RULE_FOR_NORMAL", 1); // FS
	bool const force_tsl = 1;
	int const lt_atten(solid ? 0 : 2); // sphere light atten
	float const burn_tex_scale = 0.5;
	setup_smoke_shaders(s, 0.01, 0, 1, 1, 1, 1, 1, lt_atten, 1, 0, 0, 1, force_tsl, burn_tex_scale);
	if (cobj_z_bias < 0.002)     {s.add_uniform_float("z_bias", 0.002);} // reset larger
	if (indir_vert_offset > 0.1) {s.add_uniform_float("indir_vert_offset", 0.1);} // reset smaller
	lt_atten_manager_t lt_atten_manager(s);
	if (!solid) {lt_atten_manager.enable();}
	select_no_texture();

	for (int i = 0; i < num_groups; ++i) {
		obj_group &objg(obj_groups[i]);

		if (objg.enabled && objg.temperature_ok() && objg.end_id > 0) {
			if (!(object_types[objg.type].flags & SEMI_TRANSPARENT) == solid) {
				draw_group(objg, s, lt_atten_manager);
			}
		}
	}
	if (!puddle_qbd.empty() || !obj_pld.empty()) { // use the same shader
		enable_blend();

		if (!puddle_qbd.verts.empty()) { // draw puddles
			glDepthMask(GL_FALSE);
			select_texture(BLUR_TEX);
			puddle_qbd.draw_and_clear();
			glDepthMask(GL_TRUE);
		}
		select_no_texture();
		obj_pld.draw_and_clear();
		disable_blend();
	}
	s.end_shader();

	if (!snow_pld.empty()) { // draw snowflakes from points in a custom geometry shader
		select_texture(object_types[SNOW].tid);
		glDepthMask(GL_FALSE);
		shader_t s;
		s.setup_enabled_lights(2, 1); // VS
		s.set_vert_shader("ads_lighting.part*+two_lights_no_xform");
		s.set_frag_shader("simple_texture");
		s.set_geom_shader("pt_billboard_tri"); // point => 1 triangle
		//s.set_geom_shader("output_textured_quad.part+pt_billboard"); // point => 1 quad (2 triangles)
		s.begin_shader();
		s.add_uniform_float("size", 2.0*object_types[SNOW].radius);
		s.add_uniform_int("tex0", 0);
		s.add_uniform_float("min_alpha",   0.0);
		s.add_uniform_float("color_scale", 2.0);
		check_drawing_flags(object_types[SNOW].flags, 1, s);
		s.clear_specular();
		snow_pld.draw_and_clear();
		check_drawing_flags(object_types[SNOW].flags, 0, s);
		s.add_uniform_float("color_scale", 1.0);
		s.end_shader();
		glDepthMask(GL_TRUE);
	}
}


struct wap_obj {

	int id, ndiv;
	wap_obj(int id_, int ndiv_) : id(id_), ndiv(ndiv_) {}
};


void draw_obj(obj_group &objg, vector<wap_obj> *wap_vis_objs, int type, float radius,
	const colorRGBA &color, int ndiv, int j, bool in_ammo, shader_t &shader, lt_atten_manager_t &lt_atten_manager)
{
	dwobject const &obj(objg.get_obj(j));
	point const &pos(obj.pos);
	bool const cull_face(get_cull_face(type, color));
	if (cull_face) {glEnable(GL_CULL_FACE);}
	
	switch (type) {
	case SMILEY:
		if (!(obj.flags & CAMERA_VIEW)) {
			draw_smiley(pos, obj.orientation, radius, ndiv, obj.time, obj.health, j, (in_ammo ? NULL : &objg.get_td()->get_mesh(j)), shader);
		}
		break;
	case SFPART:
		draw_smiley_part(pos, pos, obj.orientation, obj.direction, 1, ndiv, shader);
		break;
	case CHUNK:
		draw_chunk(pos, radius, obj.init_dir, obj.vdeform, (obj.flags & TYPE_FLAG), ndiv, shader);
		break;
	case SKULL:
		draw_skull(pos, obj.orientation, radius, obj.status, ndiv, obj.time, shader, (obj.direction == 1));
		break;
	case ROCKET:
		draw_rocket(pos, obj.init_dir, radius, obj.type, ndiv, obj.time, shader);
		break;
	case SEEK_D:
		draw_seekd(pos, obj.init_dir, radius, obj.type, ndiv, shader);
		break;
	case LANDMINE:
		draw_landmine(pos, radius, ndiv, obj.time, obj.source, in_ammo, shader);
		break;
	case PLASMA:
		draw_plasma(pos, pos, radius, (in_ammo ? 1.0 : obj.init_dir.x), ndiv, !in_ammo, 1, obj.time, shader);
		break;
	case GRENADE:
		draw_grenade(pos, obj.init_dir, radius, ndiv, (in_ammo ? 0 : obj.time), in_ammo, 0, shader);
		break;
	case CGRENADE:
		draw_grenade(pos, obj.init_dir, radius, ndiv, (in_ammo ? 0 : obj.time), in_ammo, 1, shader);
		break;
	case BALL:
		// Note: this is the only place where drawing an object modifies its physics state, but it's difficult to move the code
		draw_rolling_obj(pos, objg.get_obj(j).init_dir, radius, obj.status, ndiv, ((obj.flags & PLATFORM_COLL) != 0),
			dodgeball_tids[(game_mode == 2) ? (j%NUM_DB_TIDS) : 0], (in_ammo ? NULL : &objg.get_td()->get_matrix(j)), shader);
		break;
	case POWERUP:
	case HEALTH:
	case SHIELD:
		draw_powerup(pos, radius, ndiv, ((type == POWERUP) ? (int)obj.direction : -1), color, shader, lt_atten_manager);
		break;
	case WA_PACK:
		wap_vis_objs[!wid_need_weapon((int)obj.direction)].push_back(wap_obj(j, ndiv));
		break;
	case WEAPON:
		wap_vis_objs[0].push_back(wap_obj(j, ndiv));
		break;
	case AMMO:
		wap_vis_objs[1].push_back(wap_obj(j, ndiv));
		break;
	case SAWBLADE:
		draw_sawblade(pos, obj.orientation, obj.init_dir, radius, obj.angle, 1, ndiv, (obj.direction != 0));
		break;
	default:
		if (obj.vdeform != all_ones) {
			fgPushMatrix();
			translate_to(pos);
			scale_by(obj.vdeform);
			rotate_into_plus_z(obj.orientation); // might be unnesessary
			draw_sphere_vbo(zero_vector, radius, ndiv, 0);
			fgPopMatrix();
		}
		else {
			draw_sphere_vbo(pos, radius, ndiv, 0); // (object_types[type].tid >= 0)
		}
	}
	if (cull_face) {glDisable(GL_CULL_FACE);}
}


// Note: incorrect if there is both a sun and a moon
bool is_object_shadowed(dwobject &obj, float cd_scale, float radius) { // only used for snowflakes

	bool is_shadowed((obj.flags & SHADOWED) != 0); // previous value
	float const pt_size(cd_scale/distance_to_camera(obj.pos)); // approx pixel size
	int const skipval(min(20, int(8.0/pt_size)));

	if (skipval <= 1 || (obj.time % skipval) == 0) {
		is_shadowed = !is_visible_to_light_cobj(obj.pos, get_specular_light(), radius, obj.coll_id, 0);
		set_bit_flag_to(obj.flags, SHADOWED, is_shadowed);
	}
	return is_shadowed;
}


struct tid_color_to_ix_t {
	int tid;
	colorRGBA c;
	unsigned ix;

	tid_color_to_ix_t(int tid_, colorRGBA const &c_, unsigned ix_) : tid(tid_), c(c_), ix(ix_) {}

	bool operator<(tid_color_to_ix_t const &v) const {
		if (tid != v.tid) {return (tid < v.tid);}
		if (c != v.c) {return (c < v.c);}
		return (ix < v.ix);
	}
};


colorRGBA get_textured_color(int tid, colorRGBA const &color) {

	if (tid < 0) {return color;}
	return colorRGBA(texture_color(tid), 1.0).modulate_with(color); // don't use texture alpha
}


void draw_and_clear_tris(vector<vert_norm_color> &vn, vector<vert_norm_tc_color> &vntc) {

	draw_and_clear_verts(vn,   GL_TRIANGLES);
	draw_and_clear_verts(vntc, GL_TRIANGLES);
}


void draw_group(obj_group &objg, shader_t &s, lt_atten_manager_t &lt_atten_manager) {

	RESET_TIME;
	s.clear_specular();
	int const type(objg.get_ptype());
	obj_type const &otype(object_types[type]);
	int tid(otype.tid);
	float const radius(otype.radius), cd_scale(NDIV_SCALE*radius*window_width);
	unsigned const flags(otype.flags);
	bool do_texture(select_texture(tid));
	colorRGBA color(otype.color);
	s.set_cur_color(color);
	check_drawing_flags(flags, 1, s);
	int const clip_level((type == SMILEY || type == LANDMINE || type == ROCKET || type == BALL) ? 2 : 0);
	unsigned num_drawn(0);

	if (type == LEAF) { // leaves
		static vector<pair<unsigned, unsigned> > ordering;
		ordering.resize(0);
		ordering.reserve(objg.end_id);

		for (unsigned j = 0; j < objg.end_id; ++j) {
			dwobject const &obj(objg.get_obj(j));
			if (obj.disabled()) continue;
			float const leaf_scale(obj.init_dir.z);
			assert(leaf_scale > 0.0);
			if (!sphere_in_camera_view(obj.pos, leaf_scale, 0)) continue;
			int const tree_type(obj.source);
			assert(tree_type >= 0 && tree_type < NUM_TREE_TYPES);
			ordering.push_back(make_pair(tree_type, j));
		}
		if (!ordering.empty()) {
			num_drawn += (unsigned)ordering.size();
			sort(ordering.begin(), ordering.end()); // sort by texture id
			int last_tid(-1);
			if (s.is_setup()) {s.disable();}
			shader_t ls;
			setup_smoke_shaders(ls, 0.99, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1); // TSL=1
			ls.set_specular(0.75, 25.0);
			static quad_batch_draw qbd;

			for (unsigned j = 0; j < ordering.size(); ++j) {
				dwobject const &obj(objg.get_obj(ordering[j].second));
				int const tree_type(ordering[j].first), tid(tree_types[tree_type].leaf_tex);
				float const leaf_scale(obj.init_dir.z), leaf_x_ar(tree_types[tree_type].leaf_x_ar);
				assert(tid >= 0);
			
				if (draw_model == 0 && tid != last_tid) {
					qbd.draw_and_clear();
					select_texture(tid);
					last_tid = tid;
				}
				point pos(obj.pos);
				if (place_obj_on_grass(pos, leaf_scale)) {pos.z = 0.5*(obj.pos.z + pos.z-leaf_scale);} // leaf is partially on grass
				float const t(((float)obj.time)/((float)otype.lifetime));
				colorRGBA const dry_color(1.0, 0.7, 0.1); // this is the final color, even for partially burnt leaves - oh well
				colorRGBA leaf_color(WHITE);
				UNROLL_3X(leaf_color[i_] *= obj.vdeform[i_];) // vdeform is the color
				if (leaf_color.get_luminance() > 0.1) {blend_color(leaf_color, dry_color, leaf_color, t, 0);} // not mostly burned
				vector3d dirs[2] = {(leaf_points[3] - leaf_points[0]), (leaf_points[1] - leaf_points[0])};
				
				for (unsigned d = 0; d < 2; ++d) {
					dirs[d]   *= 0.5*leaf_scale;
					dirs[d].x *= leaf_x_ar;
					rotate_vector3d(plus_z, -obj.init_dir.x, dirs[d]);
					rotate_vector3d(obj.orientation, -obj.angle/TO_DEG, dirs[d]);
				}
				qbd.add_quad_dirs((pos + dirs[1]), dirs[0], -dirs[1], leaf_color, cross_product(dirs[0], dirs[1]).get_norm());
			} // for j
			qbd.draw_and_clear();
			ls.clear_specular();
			ls.end_shader();
			if (s.is_setup()) {s.enable();} // back to the original shader
		}
	} // leaf
	else if (objg.large_radius()) { // large objects
		vector<wap_obj> wap_vis_objs[2];
		bool const gm_smiley(game_mode && type == SMILEY);
		float const radius_ext(gm_smiley ? 2.0*radius : radius); // double smiley radius to account for weapons
		vector<unsigned> smiley_weapons_to_draw;

		for (unsigned j = 0; j < objg.end_id; ++j) {
			dwobject const &obj(objg.get_obj(j));
			if (obj.disabled() || ((obj.flags & CAMERA_VIEW) && type != SMILEY)) continue;
			point const &pos(obj.pos);
			if (!sphere_in_camera_view(pos, radius_ext, clip_level)) continue;
			if (type == SMILEY) smiley_weapons_to_draw.push_back(j);

			if (DEBUG_COLORCODE) {
				set_color_by_status(obj.status, s);
			}
			else if (type != SMILEY && type != SFPART && type != ROCKET && type != CHUNK &&
				type != LANDMINE && type != PLASMA && type != POWERUP && type != HEALTH && type != SHIELD)
			{
				colorRGBA color2(color);
				scale_color_uw(color2, pos);
				s.set_cur_color(color2);
			}
			++num_drawn;
			float const pt_size(cd_scale/distance_to_camera(pos));
			int const ndiv(min(N_SPHERE_DIV, max(3, int(min(pt_size, 3.0f*sqrt(pt_size))))));
			draw_obj(objg, wap_vis_objs, type, radius, color, ndiv, j, 0, s, lt_atten_manager);
		} // for j
		for (unsigned k = 0; k < wap_vis_objs[0].size(); ++k) { // draw weapons
			dwobject const &obj(objg.get_obj(wap_vis_objs[0][k].id));
			int const wid((int)obj.direction);
			vector3d v(obj.pos);
			if (wid != W_BLADE) {v.x += 0.7*radius;}
			draw_weapon_simple(v, -plus_x, 0.5*radius, obj.coll_id, wid, 0.25, s);
		}
		for (unsigned k = 0; k < wap_vis_objs[1].size(); ++k) { // draw ammo
			unsigned const j(wap_vis_objs[1][k].id);
			assert(j < objg.max_objects());
			draw_ammo(objg, radius, color, wap_vis_objs[1][k].ndiv, j, s, lt_atten_manager);
		}
		if (!wap_vis_objs[0].empty() || !wap_vis_objs[1].empty()) {
			check_drawing_flags(otype.flags, 1, s);
			select_texture(tid);
			begin_sphere_draw(0);
		}
		for (unsigned j = 0; j < 2; ++j) {
			for (unsigned k = 0; k < wap_vis_objs[j].size(); ++k) {
				wap_obj const &wa(wap_vis_objs[j][k]);
				set_obj_specular(flags, 0.5*brightness, s);
				dwobject const &obj(objg.get_obj(wa.id));
				s.set_cur_color(color);
				draw_sphere_vbo_back_to_front(obj.pos, radius, wa.ndiv, 0);
			}
		}
		if (!wap_vis_objs[0].empty() || !wap_vis_objs[1].empty()) {end_sphere_draw();}

		if (!smiley_weapons_to_draw.empty()) {
			for (vector<unsigned>::const_iterator i = smiley_weapons_to_draw.begin(); i != smiley_weapons_to_draw.end(); ++i) {
				draw_weapon_in_hand(*i, s); // Note: view culling doesn't use correct bounding sphere for all weapons
			}
		}
	} // large objects
	else { // small objects
		colorRGBA const &base_color(otype.color);
		vector<tid_color_to_ix_t> tri_fragments, sphere_fragments;
		vector<vert_norm_color> shrapnel_verts;
		vector<pair<float, unsigned> > particles_to_draw;
		int selected_particle(-1);
		if (type == PARTICLE && (rand()&3) == 0) {selected_particle = rand()%objg.end_id;}

		for (unsigned j = 0; j < objg.end_id; ++j) {
			dwobject &obj(objg.get_obj(j));
			point const &pos(obj.pos);
			if (obj.disabled() || (obj.flags & CAMERA_VIEW))      continue;
			float const tradius(obj.get_true_radius()); // differs from radius for fragments
			if (!sphere_in_camera_view(pos, tradius, clip_level)) continue;
			colorRGBA color2;
			++num_drawn;

			if (type == FRAGMENT) {
				tid = -obj.coll_id - 2; // should we sort fragments by texture id?
				UNROLL_3X(color2[i_] = obj.init_dir[i_];)
				color2.alpha = obj.vdeform.y;
			}
			else if (do_texture && type != SNOW) {
				color2 = WHITE;
			}
			else {
				color2 = base_color;
			}
			if (type != SHRAPNEL && type != PARTICLE) {
				if (type == DROPLET) {select_liquid_color(color2, pos);}
				scale_color_uw(color2, pos);
			}
			switch (type) {
			case SHELLC:
				set_color_v2(color2, obj.status, s);
				draw_shell_casing(pos, obj.orientation, obj.init_dir, tradius, obj.angle, cd_scale, obj.direction, s);
				break;
			case STAR5:
				set_color_v2(color2, obj.status, s);
				draw_star(pos, obj.orientation, obj.init_dir, tradius, obj.angle, 1);
				break;
			case SHRAPNEL:
				add_rotated_triangle(obj.pos, obj.orientation, tradius, obj.angle, get_glow_color(obj, 1), shrapnel_verts);
				break;
			case PARTICLE:
				particles_to_draw.push_back(make_pair(-distance_to_camera_sq(pos), j));
				if (animate2 && j == selected_particle && obj.time < otype.lifetime/3) {gen_smoke(pos, 0.25, 0.1);} // max one particle per frame
				break;
			case SAND:
			case DIRT:
			case ROCK:
				{
					colorRGBA tcolor(get_textured_color(tid, color2));
					color2 *= obj.orientation.y;
					if (do_texture) {tcolor *= obj.orientation.y;}
					draw_sized_point(obj, tradius, obj.orientation.x*cd_scale, color2, tcolor, do_texture, s, (type == DIRT || type == ROCK));
					break;
				}
			case FRAGMENT: // draw_fragment()?
				((obj.flags & TYPE_FLAG) ? tri_fragments : sphere_fragments).push_back(tid_color_to_ix_t(tid, color2, j)); // if shatterable, use triangle
				break;
			default:
				if (DEBUG_COLOR_COLLS) {
					int cindex;
					float const time(TIMESTEP*fticks);
					point const pos2(pos + obj.velocity*time - point(0.0, 0.0, -base_gravity*GRAVITY*time*time*otype.gravity));
					s.set_cur_color(check_coll_line(pos, pos2, cindex, -1, 0, 0) ? RED : GREEN);
					vert_wrap_t const lines[2] = {pos, pos2};
					draw_verts(lines, 2, GL_LINES);
				}
				draw_sized_point(obj, tradius, cd_scale, color2, get_textured_color(tid, color2), do_texture, s, 0);
			} // switch (type)
		} // for j
		sort(tri_fragments.begin(), tri_fragments.end()); // sort by tid
		vector<vert_norm_color> fragment_vn;
		vector<vert_norm_tc_color> fragment_vntc;
		int last_tid(-1), emission_loc(-1);

		for (vector<tid_color_to_ix_t>::const_iterator i = tri_fragments.begin(); i != tri_fragments.end(); ++i) { // Note: needs 2-sided lighting
			dwobject const &obj(objg.get_obj(i->ix));
			bool const is_emissive(obj.direction > 0); // hot object, add color
			
			if (is_emissive || i->tid != last_tid) { // emit on state change
				draw_and_clear_tris(fragment_vn, fragment_vntc);
				if (i->tid != last_tid) {select_texture(i->tid); last_tid = i->tid;}
			}
			float const tradius(obj.get_true_radius());

			if (i->tid < 0) { // not textured, use thick triangle
				add_thick_triangle(obj.pos, obj.orientation, tradius, obj.angle, 0.0, fragment_vn, 0.2*tradius, i->c);
			}
			else {
				add_rotated_textured_triangle(obj.pos, obj.orientation, tradius, obj.angle, obj.vdeform.z, i->c, fragment_vntc); // obj.vdeform.z = tscale
			}
			if (is_emissive) { // hot object, add color
				colorRGBA emissive_color(get_glow_color(obj, 1));
				emissive_color *= 2.0*(obj.direction/255.0); // can be greater than 1.0
				emissive_color.alpha = i->c.alpha;
				s.ensure_uniform_loc(emission_loc, "emission");
				s.set_uniform_color(emission_loc, emissive_color);
				draw_and_clear_tris(fragment_vn, fragment_vntc); // emit immediately
				s.set_uniform_color(emission_loc, BLACK); // clear
			}
		}
		draw_and_clear_tris(fragment_vn, fragment_vntc); // draw any remaining triangles
		sort(sphere_fragments.begin(), sphere_fragments.end()); // sort by tid

		for (vector<tid_color_to_ix_t>::const_iterator i = sphere_fragments.begin(); i != sphere_fragments.end(); ++i) {
			dwobject &obj(objg.get_obj(i->ix));
			select_texture(i->tid);
			draw_sized_point(obj, obj.get_true_radius(), cd_scale, i->c, get_textured_color(tid, i->c), (i->tid >= 0), s, 2);
		}
		if (!particles_to_draw.empty() || !shrapnel_verts.empty()) { // draw particles and shrapnel as emissive
			s.add_uniform_float("emissive_scale", 1.0); // make colors emissive

			if (!particles_to_draw.empty()) {
				sort(particles_to_draw.begin(), particles_to_draw.end()); // sort back to front
				quad_batch_draw particle_qbd;
				point const camera(get_camera_pos());

				for (auto i = particles_to_draw.begin(); i != particles_to_draw.end(); ++i) {
					dwobject const &obj(objg.get_obj(i->second));
					float const tradius(obj.get_true_radius()); // == radius in all cases?
					colorRGBA const glow_color(get_glow_color(obj, 0));
					if (glow_color.alpha > 0.0) {particle_qbd.add_billboard(obj.pos, camera, up_vector, glow_color, 1.2*tradius, 1.2*tradius);}
				}
				bool const can_use_additive_blend(0 && smoke_exists); // too much of a transition when smoke appears/disappears, so disable for now
				if (can_use_additive_blend) {set_additive_blend_mode();}
				particle_qbd.draw();
				if (can_use_additive_blend) {set_std_blend_mode();}
			}
			draw_verts(shrapnel_verts, GL_TRIANGLES);
			s.add_uniform_float("emissive_scale", 0.0); // reset
		}
	} // small object
	check_drawing_flags(flags, 0, s);
	select_no_texture();

	if (SHOW_DRAW_TIME) {
		cout << "type = " << objg.type << ", num = " << objg.end_id << ", drawn = " << num_drawn << " ";
		PRINT_TIME("Group");
	}
}


void draw_low_res_sphere_pair(point const &pos, float radius, vector3d const &v, vector3d const &vdeform,
	colorRGBA const *c1, colorRGBA const *c2, int ndiv, bool do_texture, shader_t &shader)
{
	ndiv = min(ndiv, max(3, int(3 + 1.5*(v.x + v.y + v.z))));
	translate_to(pos);
	vector3d scale((0.8+0.5*fabs(v.x)), (0.8+0.5*fabs(v.y)), (0.8+0.5*fabs(v.z)));
	scale *= radius;
	
	if (vdeform != all_ones) { // apply deformation
		float vdmin(1.0);
		UNROLL_3X(scale[i_] *= vdeform[i_]; vdmin = min(vdmin, vdeform[i_]);)
		if (vdmin < 1.0) {scale *= pow(1.0/vdmin, 1.0/3.0);}
	}
	scale_by(scale);
	fgRotate(360.0*(v.x - v.y), v.x, v.y, (v.z+0.01));
	if (c1) {shader.set_cur_color(*c1);}
	bind_draw_sphere_vbo(do_texture, 1);
	draw_sphere_vbo_pre_bound(ndiv, do_texture);
	if (c2) {shader.set_cur_color(*c2);}
	fgTranslate(0.1*(v.x-v.y), 0.1*(v.y-v.z), 0.1*(v.x-v.z));
	fgRotate(360.0*(v.z - v.x), v.y, v.z, (v.x+0.01));
	draw_sphere_vbo_pre_bound(ndiv, do_texture);
	bind_vbo(0);
}


void draw_sized_point(dwobject &obj, float radius, float cd_scale, const colorRGBA &color, const colorRGBA &tcolor,
					  bool do_texture, shader_t &shader, int is_chunky)
{
	point pos(obj.pos);
	point const camera(get_camera_pos());
	float point_dia(cd_scale/p2p_dist(camera, pos));
	if (do_zoom) point_dia *= ZOOM_FACTOR;
	int const type(obj.type);
	bool const draw_large(point_dia >= 2.5);
	bool const draw_snowflake(draw_large && type == SNOW);
	bool const tail_type((object_types[type].flags & TAIL_WHEN_FALL) != 0);
	bool const tail(tail_type && obj.status == 1 && obj.velocity.z < RAIN_TAIL_MIN_V && !(obj.flags & OBJ_COLLIDED));

	if (tail && !draw_large) { // draw rain with lines
		point pos2(pos);
		pos2.z -= 2.0*fticks*TIMESTEP*obj.velocity.z;
		colorRGBA color2(color);
		color2.alpha *= min(1.0f, 0.5f*point_dia);
		obj_pld.add_line(pos2, (camera - pos2), ALPHA0, pos, (camera - pos), color2);
		return;
	}
	bool const precip((object_types[type].flags & IS_PRECIP) != 0);
	bool const is_drop(type != SNOW && (object_types[type].flags & OBJ_IS_DROP)); // can get snow when changing temp
	if (!precip && point_dia < 0.4) return; // clip it

	if (is_drop && obj.status == 4 && (obj.flags & STATIC_COBJ_COLL)) { // draw as puddle
		if (!draw_large || is_underwater(pos)) return; // don't draw
		assert(!do_texture);
		colorRGBA color2(color);
		if (type == RAIN) color2.alpha *= 0.5; // rain is mostly transparent when small
		puddle_qbd.add_billboard(pos, (pos + plus_z), plus_x, color2, 5.0*radius, 5.0*radius);
		return;
	}
	if (draw_snowflake) { // draw as a point to be converted to a billboard by the geometry shader
		bool const is_shadowed(is_object_shadowed(obj, cd_scale, radius));
		// Note: color is scaled by 0.5 here (and 2.0 in the shader to cancel) to allow for blue > 1.0
		snow_pld.add_pt(pos, (is_shadowed ? zero_vector : (get_light_pos() - pos)), (do_texture ? tcolor : color)*0.5);
		return;
	}
	if (!draw_large) { // draw as a point
		bool const scatters(type == RAIN || type == SNOW);
		vector3d const n((scatters ? get_light_pos() : camera) - pos);
		obj_pld.add_pt(pos, n, (do_texture ? tcolor : color));
		return;
	}
	set_color_v2(color, obj.status, shader);
	bool const cull_face(get_cull_face(type, color));
	fgPushMatrix();

	if (cull_face) {
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
	}

	// draw as a sphere
	if (is_chunky) {
		assert(!tail);
		vector3d const v((is_chunky == 2) ? obj.orientation : obj.init_dir);
		draw_low_res_sphere_pair(pos, radius, v, all_ones, NULL, NULL, 16, do_texture, shader);
	}
	else {
		int ndiv(int(4.0*sqrt(point_dia)));

		if (is_droplet(type)) {
			ndiv = min(ndiv/2, N_SPHERE_DIV/2);
		}
		else if (type == ROCK || type == SAND || type == DIRT || type == FRAGMENT) {
			ndiv /= 2;
		}
		ndiv = max(4, min(ndiv, N_SPHERE_DIV));
	
		if (ndiv > 3 && tail) { // cone on the tail of the raindrop
			draw_fast_cylinder(pos, pos+point(0.0, 0.0, 2.5*radius), radius, 0.0, (ndiv>>1), 0);
			translate_to(pos + vector3d(0.0, 0.0, -0.6*radius));
			fgScale(1.0, 1.0, 2.0);
			pos = all_zeros;
		}
		draw_sphere_vbo(pos, radius, ndiv, do_texture);
	}
	if (cull_face) glDisable(GL_CULL_FACE);
	fgPopMatrix();
}


void draw_ammo(obj_group &objg, float radius, const colorRGBA &color, int ndiv, int j, shader_t &shader, lt_atten_manager_t &lt_atten_manager) {

	dwobject const &obj(objg.get_obj(j));
	point pos(obj.pos);
	vector<wap_obj> wap_vis_objs[2]; // not actually used
	int const atype(get_ammo_or_obj((int)obj.direction));
	if (atype < 0) return; // can this happen?
	obj_type const &otype(object_types[atype]);
	check_drawing_flags(otype.flags, 1, shader);
	if (otype.tid >= 0) {select_texture(otype.tid);}
	bool const cull_face(get_cull_face(atype, color));
	if (cull_face) {glEnable(GL_CULL_FACE);}

	switch (atype) {
	case SHELLC: // M16
		pos.z -= 0.5*radius;
		set_brass_material(shader); // looks gray through the blue ammo shell
		draw_cylinder_at(pos, 1.0*radius, 0.2*radius, 0.2*radius, ndiv, 1);
		pos.z += radius;
		draw_sphere_vbo(pos, 0.2*radius, ndiv, 0);
		break;
	case PROJECTILE: // shotgun
		for (unsigned n = 0; n < 2; ++n) { // two shells in one ammo
			point pos2(pos);
			pos2.x += (1.0 - 2.0*n)*0.3*radius;
			shader.set_cur_color(RED);
			pos2.z -= 0.5*radius;
			draw_cylinder_at(pos2, 1.2*radius, 0.3*radius, 0.3*radius, ndiv, 1);
			set_brass_material(shader); // looks gray through the blue ammo shell
			pos2.z -= 0.2*radius;
			draw_cylinder_at(pos2, 0.4*radius, 0.32*radius, 0.32*radius, ndiv, 1);
			check_drawing_flags(otype.flags, 1, shader);
		}
		break;
	case BEAM: // laser
		set_emissive_only(RED, shader);
		pos.z -= 0.5*radius;
		draw_cylinder_at(pos, 1.0*radius, 0.1*radius, 0.1*radius, ndiv, 1);
		shader.clear_color_e();
		break;
	case STAR5: // throwing star
		shader.set_cur_color(otype.color);
		draw_star(pos, obj.orientation, obj.init_dir, 0.4*radius, obj.angle, 0);
		break;
	case GASSED:
		shader.set_cur_color(otype.color);
		draw_sphere_vbo(pos, 0.6*radius, ndiv, 1);
		break;
	default:
		shader.set_cur_color(otype.color);
		draw_obj(objg, wap_vis_objs, atype, 0.4*radius, color, ndiv, j, 1, shader, lt_atten_manager);
	}
	if (cull_face) {glDisable(GL_CULL_FACE);}
	if (otype.tid >= 0) {select_no_texture();}
}


colorRGBA get_powerup_color(int powerup) {

	switch (powerup) {
	case PU_NONE:         return BLACK;
	case PU_DAMAGE:       return CYAN;
	case PU_REGEN:        return BLUE;
	case PU_SHIELD:       return GREEN;
	case PU_SPEED:        return WHITE;
	case PU_FLIGHT:       return PURPLE;
	case PU_INVISIBILITY: return BLACK; // not useful for smileys
	}
	return BLACK;
}


inline void rotate_to_dir(vector3d const &dir) { // normalized to +y (for smileys)

	fgRotate(atan2(dir.y, dir.x)*TO_DEG-90.0, 0.0, 0.0, 1.0);
	fgRotate(safe_acosf(-dir.z)*TO_DEG-90.0, 1.0, 0.0, 0.0);
}


void draw_smiley_part(point const &pos, point const &pos0, vector3d const &orient, int type, int use_orient, int ndiv, shader_t &shader, float scale) {

	assert(type < NUM_SMILEY_PARTS);
	float const radius(scale*object_types[SFPART].radius);
	colorRGBA const sf_color[NUM_SMILEY_PARTS] = {BLACK, RED, PINK};
	shader.set_cur_color(sf_color[type]);

	switch (type) {
	case SF_EYE:
		draw_sphere_vbo(pos, 1.0*radius, ndiv, 0);
		break;
	case SF_NOSE:
		draw_sphere_vbo(pos, 1.2*radius, ndiv, 0);
		break;
	case SF_TONGUE:
		fgPushMatrix();
		translate_to(pos);
		if (use_orient) rotate_to_dir(orient);
		fgScale(1.5*radius, 3.0*radius, 0.375*radius);
		draw_sphere_vbo(all_zeros, 1.0, ndiv, 0);
		fgPopMatrix();
		break;
	default: assert(0);
	}
}


colorRGBA mult_alpha(colorRGBA const &c, float alpha) {
	return colorRGBA(c.R, c.G, c.B, c.A*alpha);
}


void draw_smiley(point const &pos, vector3d const &orient, float radius, int ndiv, int time,
				 float health, int id, mesh2d const *const mesh, shader_t &shader)
{
	colorRGBA color;
	int const powerup(sstates[id].powerup), ndiv2(max(3, (ndiv>>1)));
	fgPushMatrix();
	translate_to(pos);
	rotate_to_dir(orient);
	point pos2(-0.4*radius, 0.85*radius, 0.3*radius);

	// draw eyes
	for (unsigned i = 0; i < 2; ++i) {
		if (health > 10.0) {
			float const scale((powerup == PU_SPEED) ? 1.5 : 1.0);
			point const pos3 ((powerup == PU_SPEED) ? (pos2 + point(0.0, 0.2*radius, 0.0)) : pos2);
			draw_smiley_part(pos3, pos, orient, SF_EYE, 0, ndiv2, shader, scale); // eyes
		}
		else {
			shader.set_cur_color(BLACK);

			for (unsigned l = 0; l < 2; ++l) {
				point pts[2];
				for (unsigned p = 0; p < 2; ++p) {pts[p] = pos2 + radius*point((p ? -0.12 : 0.12), 0.1, ((p^l) ? -0.12 : 0.12));}
				draw_fast_cylinder(pts[0], pts[1], 0.02*radius, 0.02*radius, ndiv2, 0);
			}
		}
		pos2.x *= -1.0;
	}

	// draw nose
	if (powerup != PU_INVISIBILITY || same_team(id, -1)) { // show nose even if invisible if same team as player
		point pos3(0.0, 1.1*radius, 0.0);
		draw_smiley_part(pos3, pos, orient, SF_NOSE, 0, ndiv2, shader); // nose
	}
	float alpha(1.0);

	switch (powerup) {
		case PU_DAMAGE: // devil horns
			shader.set_cur_color(RED);
			draw_cylinder_at(point( 0.3*radius, 0.7*radius, 0.6*radius), 0.6*radius, 0.1*radius, 0.0, ndiv2, 0);
			draw_cylinder_at(point(-0.3*radius, 0.7*radius, 0.6*radius), 0.6*radius, 0.1*radius, 0.0, ndiv2, 0);
			break;

		case PU_REGEN: // raindrops
			if (animate2) {
				int const cid(coll_id[WDROPLET]), k(obj_groups[cid].choose_object());
				obj_groups[cid].create_object_at(k, (pos + point(0.0, 0.0, 4.0*radius)));
				vadd_rand(obj_groups[cid].get_obj(k).pos, 0.04);
			}
			break;

		case PU_SHIELD: // shield + gas mask
			// gas generated in gameplay.cpp::update_game_frame()
			break;

		case PU_SPEED:
			if (animate2 && !(rand()&3)) gen_smoke(pos);
			break;

		case PU_FLIGHT: // propeller or wings?
			shader.set_cur_color(BLACK);
			draw_cylinder(0.5*radius, 0.05*radius, 0.05*radius, ndiv2, 0, 0, 0, 0.9*radius);
			fgPushMatrix();
			fgTranslate(0.0, 0.0, 1.4*radius);
			fgRotate(float((30*time)%360), 0.0, 0.0, 1.0);
			fgScale(1.0, 0.25, 0.05);
			draw_sphere_vbo(point( 0.5*radius, 0.0, 0.0), 0.5*radius, ndiv, 0); // propeller
			draw_sphere_vbo(point(-0.5*radius, 0.0, 0.0), 0.5*radius, ndiv, 0); // propeller
			fgPopMatrix();
			break;

		case PU_INVISIBILITY: {
				float const put(float(sstates[id].powerup_time)/TICKS_PER_SECOND), init_put(float(POWERUP_TIME)/TICKS_PER_SECOND);
				if ((init_put - put) < 1.0) {alpha = (1.0 - (init_put - put));} // fading out
				else if (put < 1.0)         {alpha = (1.0 - put);} // fading in
				else { // fully invisible
					fgPopMatrix();
					return;
				}
				if (alpha < 1.0) enable_blend();
				break;
			}
	} // switch (powerup)
	if (powerup != PU_INVISIBILITY && time%10 < 5 && powerup >= 0) { // powerup
		shader.set_cur_color(get_powerup_color(powerup));
	}
	else if (health >= 50.0) {
		shader.set_cur_color(mult_alpha(YELLOW, alpha));
	}
	else {
		shader.set_cur_color(colorRGBA(1.0, (0.25 + 0.015*health), 0.0, alpha));
	}
	if (game_mode == 2) { // dodgeball
		select_texture(CAMOFLAGE_TEX);
	}
	else {
		select_smiley_texture(id);
	}
	if (mesh) { // main body
		mesh->draw_perturbed_sphere(all_zeros, radius, ndiv, 1);
		//if (mesh->size > 0) ndiv = mesh->size;
	}
	else {
		draw_sphere_vbo(all_zeros, radius, ndiv, 1);
	}
	select_no_texture();

	if (teams > 1) { // draw team headband
		shader.set_cur_color(mult_alpha(get_smiley_team_color(id), alpha));
		fgPushMatrix();
		fgScale(1.0, 1.0, 0.5);
		draw_sphere_vbo(point(0.0, 0.0, 0.9*radius), 0.94*radius, ndiv, 0);
		fgPopMatrix();
	}

	// draw unique identifier
	shader.set_cur_color(mult_alpha(get_smiley_team_color(id+1, 1), alpha)); // ignore teams and use max_colors
	fgPushMatrix();
	fgScale(1.0, 1.0, 0.3);
	draw_sphere_vbo(point(0.0, 0.0, (0.8/0.3)*radius), 0.65*radius, ndiv, 0);
	fgPopMatrix();
	
	// draw mouth
	float const hval(0.004*(100.0 - min(160.0f, health)));
	shader.set_cur_color(mult_alpha(BLACK, alpha));
	point const pts[4] = {point(-0.5, 0.95, -0.2-hval), point(-0.15, 0.95, -0.4), point( 0.15, 0.95, -0.4), point( 0.5, 0.95, -0.2-hval)};

	for (unsigned i = 0; i < 3; ++i) {
		draw_fast_cylinder(radius*pts[i], radius*pts[i+1], 0.04*radius, 0.04*radius, ndiv2, 0);
	}

	// draw tongue
	if (sstates[id].kill_time < int(2*TICKS_PER_SECOND) || powerup == PU_DAMAGE) { // stick your tongue out at a dead enemy
		point pos4(0.0, 0.8*radius, -0.4*radius);
		draw_smiley_part(pos4, pos, orient, SF_TONGUE, 0, ndiv2, shader);
	}
	if (game_mode == 2 && (sstates[id].p_ammo[W_BALL] > 0 || UNLIMITED_WEAPONS)) { // dodgeball
		select_texture(select_dodgeball_texture(id));
		shader.set_cur_color(mult_alpha(object_types[BALL].color, alpha));
		draw_cube_mapped_sphere(point(0.0, 1.3*radius, 0.0), 0.8*object_types[BALL].radius, ndiv/2, 1);
		select_no_texture();
	}
	fgPopMatrix();
	vector3d hit_dir;
	int const hit(get_smiley_hit(hit_dir, id));

	if (hit > 0) { // hit - draw damage or shields
		enable_blend();
		fgPushMatrix();
		translate_to(pos);
		rotate_sphere_tex_to_dir(hit_dir);
		select_texture(SBLUR_TEX);

		if (powerup == PU_SHIELD || sstates[id].shields > 0.01) { // new shields mode
			shader.disable();
			shader_t shield_shader;
			setup_shield_shader(shield_shader, 11);
			shield_shader.set_cur_color(colorRGBA(GREEN, alpha*hit/HIT_TIME));
			set_additive_blend_mode();
			draw_sphere_vbo(all_zeros, 1.015*radius, ndiv, 1);
			set_std_blend_mode();
			shield_shader.end_shader();
			shader.enable();
		}
		else {
			shader.set_cur_color(colorRGBA(BLOOD_C, alpha*hit/HIT_TIME)); // black color for burns?
			shader.add_uniform_float("min_alpha", 0.05);
			draw_sphere_vbo(all_zeros, 1.015*radius, ndiv, 1);
			shader.add_uniform_float("min_alpha", 0.01);
		}
		fgPopMatrix();

		if (powerup == PU_SHIELD) {
			select_texture(PLASMA_TEX);
			shader.set_cur_color(colorRGBA(PURPLE, alpha)); // not scaled
			draw_sphere_vbo(pos, 1.05*radius, ndiv, 1);
		}
		select_no_texture();
		disable_blend();
	}
	if (alpha < 1.0) {disable_blend();}
}


void draw_powerup(point const &pos, float radius, int ndiv, int type, const colorRGBA &color, shader_t &shader, lt_atten_manager_t &lt_atten_manager) {

	ndiv = 3*ndiv/2; // increase ndiv for better transparency effect
	colorRGBA const ecolor(((type == -1) ? color : get_powerup_color(type)), 0.02); // low alpha
	set_emissive_only(ecolor, shader);
	lt_atten_manager.next_sphere(40.0, 1.0, pos, 0.7*radius); // dense gas
	draw_sphere_vbo(pos, 0.7*radius, ndiv, 0); // draw flare/billboard?
	shader.clear_color_e();
	shader.set_cur_color(colorRGBA(color, 0.05)); // low alpha
	lt_atten_manager.next_sphere(20.0, 1.6, pos, radius); // light glass
	draw_sphere_vbo(pos, radius, ndiv, 0);
}


void draw_rolling_obj(point const &pos, point &lpos, float radius, int status, int ndiv, bool on_platform, int tid, xform_matrix *matrix, shader_t &shader) {

	select_texture(tid);
	fgPushMatrix();
	translate_to(pos);
	
	if (matrix) {
		if (on_platform) {lpos = pos;} // reset now so there's no rotation
		apply_obj_mesh_roll(*matrix, pos, lpos, radius, ((status == 1) ? 0.01 : 0.0), ((status == 1) ? 0.2 : 1.0));
	}
	//draw_sphere_vbo(all_zeros, radius, 2*ndiv, 1);
	draw_cube_mapped_sphere(all_zeros, radius, ndiv, 1);
	fgPopMatrix();
	lpos = pos;
}


void draw_skull(point const &pos, vector3d const &orient, float radius, int status, int ndiv, int time, shader_t &shader, bool burned) {

	float const burn_val(burned ? 0.5*(2.0*(float)time/((float)object_types[SKULL].lifetime) - 1.0) : -1.0);
	if (burn_val > -1.0) {shader.add_uniform_float("burn_offset", burn_val);}
	shader.add_uniform_float("min_alpha", 0.9);
	fgPushMatrix();
	translate_to(pos);
	rotate_from_v2v(orient, vector3d(0.0, -1.0, 0.0));
	fgRotate(180.0, 0.0, 1.0, 0.0);
	draw_sphere_vbo(all_zeros, radius, 2*ndiv, 1);
	fgPopMatrix();
	shader.add_uniform_float("min_alpha", 0.01);
	if (burn_val > -1.0) {shader.add_uniform_float("burn_offset", -1.0);} // reset
}


void draw_rocket(point const &pos, vector3d const &orient, float radius, int type, int ndiv, int time, shader_t &shader) {

	fgPushMatrix();
	translate_to(pos);
	rotate_by_vector(orient);
	shader.set_cur_color(RED);
	uniform_scale(radius);
	vert_wrap_t const verts[6] = {point(0.0, 0.0, 0.0), point(1.8, 0.0, -2.0), point(-1.8, 0.0, -2.0), point(0.0, 0.0, 0.0), point( 0.0, 1.8, -2.0), point(0.0, -1.8, -2.0)};
	draw_verts(verts, 6, GL_TRIANGLES);
	shader.set_cur_color(object_types[ROCKET].color);
	fgScale(1.0, 1.0, 2.0);
	draw_sphere_vbo_raw(ndiv, 0);
	draw_cylinder(-1.1, 1.0, 1.0, ndiv);
	fgPopMatrix();
	if (type == ROCKET) {gen_rocket_smoke(pos, orient, radius);}
}


void draw_seekd(point const &pos, vector3d const &orient, float radius, int type, int ndiv, shader_t &shader) {

	fgPushMatrix();
	translate_to(pos);
	rotate_by_vector(orient);
	uniform_scale(radius);
	draw_fast_cylinder(point(0.0, 0.0, -2.0 ), point(0.0, 0.0, -1.0), 1.0, 0.0, ndiv, 0);
	shader.set_cur_color(BLACK);
	draw_fast_cylinder(point(0.0, 0.0, -2.25), point(0.0, 0.0, -2.0), 1.0, 1.0, ndiv, 0);
	fgScale(1.0, 1.0, 1.5);
	fgRotate(90.0, -1.0, 0.0, 0.0);
	fgRotate(90.0,  0.0, 1.0, 0.0);
	shader.set_cur_color(WHITE);
	select_texture(SKULL_TEX);
	draw_sphere_vbo_raw(ndiv, 1);
	select_no_texture();
	fgPopMatrix();
	if (type == SEEK_D) {gen_rocket_smoke(pos, orient, radius);}
}


colorRGBA const &get_landmine_light_color(int time) {

	return ((time < 40) ? GREEN : (((time/6)&1) ? RED : BLUE));
}


float get_landmine_sensor_height(float radius, int time) {

	return ((time <= 6) ? 0 : ((time > 16) ? 1.5*radius : (radius + 0.05*radius*(time - 6))));
}


void draw_landmine(point pos, float radius, int ndiv, int time, int source, bool in_ammo, shader_t &shader) {

	assert(radius > 0.0 && ndiv > 0);

	if (!in_ammo) {
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

		if (!point_outside_mesh(xpos, ypos) && !is_mesh_disabled(xpos, ypos) &&
			pos.z < (interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1) + 1.05*radius))
		{
			pos.z -= 0.8*radius; // appears to sink into the ground
		}
	}
	if (!DEBUG_COLORCODE) {shader.set_cur_color(WHITE);}
	draw_sphere_vbo(pos, radius, ndiv, 1); // main body
	select_no_texture();
	fgPushMatrix();
	translate_to(pos);
	float const val(get_landmine_sensor_height(radius, time));

	if (time > 6) {
		shader.set_cur_color(GRAY);
		draw_cylinder(val, 0.05*radius, 0.05*radius, ndiv/2); // sensor pole
	}
	pos.z += val;

	if (time > 26) {
		float val;
		fgTranslate(0.0, 0.0, 1.4*radius);
		fgRotate(90.0, 1.0, 0.0, 0.0);
		
		if (time > 36) {
			fgRotate(10.0*(time%36), 0.0, 1.0, 0.0);
			val = 0.42*radius;
		}
		else {
			val = 0.04*radius*(time - 26) + 0.02*radius;
		}
		if (teams > 1) {shader.set_cur_color(get_smiley_team_color(source));} // use team color
		draw_cylinder(0.4*radius, 0.0, val, ndiv); // sensor
	}
	fgPopMatrix();

	if (time > 5) {
		pos.z += 0.15*radius;
		shader.set_color_e(get_landmine_light_color(time));
		draw_sphere_vbo(pos, 0.15*radius, max(3, ndiv/2), 0); // warning light
		shader.clear_color_e();
	}
	select_texture(object_types[LANDMINE].tid);
}


colorRGBA get_plasma_color(float size) {
	return colorRGBA(1.0, size/5.0, max(0.0f, 0.5f*(size-3.0f)), 0.9);
}


void draw_plasma(point const &pos, point const &part_pos, float radius, float size, int ndiv, bool gen_parts, bool add_halo, int time, shader_t &shader) {

	if (animate2) {radius *= rand_uniform(0.99, 1.01) + 0.1*(0.5 + 0.1*(abs((time % 20) - 10)));}

	if (1) { // slower, but looks much nicer
		select_texture(WHITE_TEX); // texture is procedural
		draw_one_star(RED, YELLOW, pos, size*radius, ndiv, add_halo);
		shader.make_current();
	}
	else {
		colorRGBA const color(get_plasma_color(size + 0.5*(0.5 + 0.16*abs((time % 12) - 6))));
		set_emissive_only(color, shader);
		//draw_sphere_vbo(pos, size*radius, ndiv, 1);
		draw_cube_mapped_sphere(pos, size*radius, ndiv/2, 1);
		shader.clear_color_e();
	}
	if (gen_parts && animate2 && !is_underwater(part_pos, 1) && (rand()&15) == 0) {gen_particles(part_pos, 1);}
}


void draw_chunk(point const &pos, float radius, vector3d const &v, vector3d const &vdeform, int charred, int ndiv, shader_t &shader) {

	fgPushMatrix();
	draw_low_res_sphere_pair(pos, radius*(0.5 + fabs(v.x)), v, vdeform, &(charred ? BLACK : YELLOW), &(charred ? DK_GRAY : BLOOD_C), ndiv, 0, shader);
	fgPopMatrix();
}


void draw_grenade(point const &pos, vector3d const &orient, float radius, int ndiv, int time, bool in_ammo, bool is_cgrenade, shader_t &shader) {

	fgPushMatrix();
	translate_to(pos);
	uniform_scale(radius);
	fgPushMatrix();
	if (!is_cgrenade) {fgScale(0.8, 0.8, 1.2);} // rotate also?
	//shader.set_cur_color(BLACK);
	(is_cgrenade ? set_gold_material(shader) : set_copper_material(shader));
	draw_sphere_vbo_raw(ndiv, 0);
	set_obj_specular(object_types[GRENADE].flags, brightness, shader);
	fgPopMatrix();

	float const stime(1.0 - float(time)/float(object_types[is_cgrenade ? CGRENADE : GRENADE].lifetime)), sval(0.2 + 0.8*stime);
	vector3d const vr((orient.x == 0.0 && orient.y == 0.0) ? vector3d(1.0, 0.0, 0.0) : vector3d(orient.x, orient.y, 0.0));
	vector3d vd(plus_z);
	rotate_vector3d_norm(vr, -0.25*PI, vd);
	rotate_about(45.0, vr);
	draw_fast_cylinder(point(0.0, 0.0, 0.7), point(0.0, 0.0, 1.2), 0.3, 0.3, max(3, ndiv/2), 0);
	shader.set_cur_color(GRAY);
	draw_fast_cylinder(point(0.0, 0.0, 1.0), point(0.0, 0.0, 1.0+sval), 0.05, 0.05, max(3, ndiv/2), 0); // fuse
	fgPopMatrix();

	if (!animate2) return;
	point const spos(pos + vd*((1.0 + sval)*radius));
	colorRGBA scolor;
	blend_color(scolor, YELLOW, ORANGE, rand_uniform(0.3, 0.7), 1);
	float const size(radius*rand_uniform(0.5, 0.7));
	sparks.push_back(spark_t(spos, scolor, size));
	add_dynamic_light(0.15, spos, scolor); // out of sync by a frame?
	if (!in_ammo && (rand()&15) == 0) {gen_particles(spos, 1, 0.5, 1);}
}


void draw_star(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate) {
	
	fgPushMatrix();
	translate_to(pos);

	if (rotate) {
		rotate_by_vector(init_dir, -90.0);
		rotate_about(angle, orient);
	}
	vert_norm points[3*N_STAR_POINTS];

	for (int i = N_STAR_POINTS-1, ix = 0; i >= 0; --i) { // Note: needs 2-sided lighting
		points[ix++] = vert_norm(2.0*radius*star_pts[(i == 0) ? (N_STAR_POINTS<<1)-1 : (i<<1)-1], orient);
		points[ix++] = vert_norm(2.0*radius*star_pts[i<<1], orient);
		points[ix++] = vert_norm(2.0*radius*star_pts[(i<<1)+1], orient);
	}
	draw_verts(points, 3*N_STAR_POINTS, GL_TRIANGLES);
	fgPopMatrix();
}


void draw_sawblade(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate, int ndiv, bool bloody) {
	
	fgPushMatrix();
	translate_to(pos);

	if (rotate) {
		rotate_by_vector(init_dir, -90.0);
		rotate_about(angle, orient);
	}
	select_texture(bloody ? SAW_B_TEX : SAW_TEX);
	draw_circle_normal(0.0, radius, ndiv, 0);
	fgPopMatrix();
}


void draw_shell_casing(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius,
					   float angle, float cd_scale, unsigned char type, shader_t &shader)
{
	assert(type == 0 || type == 1);
	float const point_size(cd_scale/distance_to_camera(pos));
	int const ndiv(max(3, min(N_SPHERE_DIV/2, int(point_size))));
	fgPushMatrix();
	translate_to(pos);
	fgRotate(TO_DEG*init_dir.x, 0.0, 0.0, 1.0);
	//rotate_by_vector(init_dir);
	rotate_about(angle, orient);
	uniform_scale(radius); // Note: needs 2-sided lighting

	if (type == 0) { // M16 shell casing
		set_brass_material(shader);
		draw_cylinder(4.0, 1.0, 1.0, ndiv);
	}
	else { // shotgun shell casing
		shader.set_cur_color(RED);
		draw_fast_cylinder(point(0.0, 0.0, -2.0), point(0.0, 0.0,  2.8), 1.2,  1.2,  ndiv, 0);
		set_brass_material(shader);
		draw_fast_cylinder(point(0.0, 0.0, -2.8), point(0.0, 0.0, -1.2), 1.28, 1.28, ndiv, 0);
	}
	set_obj_specular(object_types[SHELLC].flags, 0.5*brightness, shader); // reset
	if (point_size > 1.0) {draw_circle_normal(0, ((type == 0) ? 1.0 : 1.28), ndiv, 0, ((type == 0) ? 0.0 : -2.8));}
	fgPopMatrix();
}


colorRGBA get_glowing_obj_color(point const &pos, int time, int lifetime, float &stime, bool shrapnel_cscale, bool fade) {

	stime = ((float)time)/((float)lifetime)*(shrapnel_cscale ? 1.0 : 0.3);
	if (is_underwater(pos)) stime *= 5.0;
	stime = min(1.0f, 5.0f*stime);

	if (fade) {
		return colorRGBA(1.0, min(1.0, (2.0 - 2.0*stime)), max(0.0f, (1.0f - 1.5f*stime)), min(1.0, (4.0 - 4.0*stime)));
	}
	else {
		return colorRGBA((1.0 - 0.9*stime), max(0.0f, (0.9f - 2.0f*stime)), max(0.0f, (0.6f - 4.0f*stime)), 1.0);
	}
}


colorRGBA get_glow_color(dwobject const &obj, bool shrapnel_cscale) {

	float stime;
	colorRGBA color(get_glowing_obj_color(obj.pos, obj.time, object_types[obj.type].lifetime, stime, shrapnel_cscale, ((obj.flags & TYPE_FLAG) != 0)));
	if (shrapnel_cscale) color *= CLIP_TO_01(1.0f - stime);
	return color;
}


void update_precip_rate(float val) {
	obj_groups[coll_id[PRECIP]].update_app_rate(val, 2, 1000);
	if (val < 1.0) {obj_pld.free_mem();}
}

unsigned get_precip_rate() {return obj_groups[coll_id[PRECIP]].app_rate;}
float get_rain_intensity() {return (is_rain_enabled() ? min(get_precip_rate()/100.0, 1.0) : 0.0);}

