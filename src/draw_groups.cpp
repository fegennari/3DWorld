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


bool const DEBUG_COLORCODE   = 0;
bool const DEBUG_COLOR_COLLS = 0;
bool const SHOW_DRAW_TIME    = 0;
float const NDIV_SCALE       = 1.6;


struct puddle_t {
	point pos;
	float radius;
	colorRGBA color;

	puddle_t() {}
	puddle_t(point const &pos_, float radius_, colorRGBA const &color_) : pos(pos_), radius(radius_), color(color_) {}
};


// Global Variables
vector<puddle_t> puddles;
pt_line_drawer obj_pld;
pt_line_drawer_hdr snow_pld;


extern bool underwater;
extern int display_mode, num_groups, teams, begin_motion, UNLIMITED_WEAPONS;
extern int window_width, window_height, game_mode, draw_model, animate2;
extern float fticks, TIMESTEP, base_gravity, brightness, indir_vert_offset, cobj_z_bias;
extern point star_pts[];
extern vector3d up_norm;
extern GLUquadricObj* quadric;
extern vector<spark_t> sparks;
extern obj_group obj_groups[];
extern obj_type object_types[];
extern player_state *sstates;
extern int coll_id[];



void draw_group(obj_group &objg, shader_t &s);
void draw_sized_point(dwobject &obj, float radius, float cd_scale, const colorRGBA &color, const colorRGBA &tcolor,
					  bool do_texture, int is_chunky=0);
void draw_weapon2(dwobject const &obj, float radius);
void draw_ammo(obj_group &objg, float radius, const colorRGBA &color, int ndiv, int j);
void draw_smiley_part(point const &pos, point const &pos0, vector3d const &orient, int type,
					  int use_orient, int ndiv, float scale=1.0);
void draw_smiley(point const &pos, vector3d const &orient, float radius, int ndiv, int time,
				 float health, int id, mesh2d const *const mesh);
void draw_powerup(point const &pos, float radius, int ndiv, int type, const colorRGBA &color);
void draw_rolling_obj(point const &pos, point &lpos, float radius, int status, int ndiv, bool on_platform, int tid, xform_matrix *matrix);
void draw_skull(point const &pos, vector3d const &orient, float radius, int status, int ndiv);
void draw_rocket(point const &pos, vector3d const &orient, float radius, int type, int ndiv, int time);
void draw_seekd(point const &pos, vector3d const &orient, float radius, int type, int ndiv);
void draw_landmine(point pos, float radius, int ndiv, int time, int source, bool in_ammo);
void draw_plasma(point const &pos, point const &part_pos, float radius, float size, int ndiv, int shpere_tex, bool gen_parts, int time);
void draw_chunk(point const &pos, float radius, vector3d const &v, vector3d const &vdeform, int charred, int ndiv);
void draw_grenade(point const &pos, vector3d const &orient, float radius, int ndiv, int time, bool in_ammo, bool is_cgrenade);
void draw_star(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate);
void draw_shell_casing(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius,
					   float angle, float cd_scale, unsigned char type);
colorRGBA get_glow_color(dwobject const &obj, bool shrapnel_cscale);

int same_team(int source, int target); // gameplay



void set_obj_specular(unsigned flags, float specular_brightness) {

	if (flags & SPECULAR) {
		set_specular(specular_brightness, 50.0);
	}
	else if (flags & LOW_SPECULAR) {
		set_specular(0.4*specular_brightness, 10.0);
	}
}


void check_drawing_flags(unsigned flags, int init_draw) {

	if (init_draw) {
		set_obj_specular(flags, 0.5*brightness);
		if (flags & BLEND) enable_blend();
	}
	else {
		if (flags & BLEND) disable_blend();
		if (flags & (SPECULAR | LOW_SPECULAR)) set_specular(0.0, 1.0);
	}
}


void set_color_by_status(int status) {

	colorRGBA const colors[6] = {BLACK, RED, WHITE, YELLOW, BLUE, GRAY};
	assert(status >= 0 && status < 6);
	set_color_alpha(colors[status]);
}


void set_color_v2(const colorRGBA &color, int status) {

	if (DEBUG_COLORCODE) {
		set_color_by_status(status);
	}
	else {
		set_color_alpha(color);
	}
}


void set_emissive_color_obj(colorRGBA const &color) {
	color.do_glColor();
	//set_emissive_color(color);
	//set_color(color);
}


inline bool is_droplet(int type) {
	return ((object_types[type].flags & OBJ_IS_DROP) || type == HAIL || type == CHARRED);
}

inline bool get_cull_face(int type, colorRGBA const &color) {
	return (color.alpha < 1.0 && type != ROCKET && type != STAR5);
}

void draw_unit_sphere(int ndiv, bool do_texture) {
	draw_sphere_dlist_raw(ndiv, do_texture);
}

void select_no_texture() {
	//glDisable(GL_TEXTURE_2D);
	select_texture(WHITE_TEX, 0);
}


void scale_color_uw(colorRGBA &color, point const &pos) {

	water_color_atten_at_pos(color, pos); // ???
	if (!underwater) return;
	color.R *= 0.45;
	color.G *= 0.45;
	color.B *= 0.85;
}


void draw_rotated_triangle(point const &pos, vector3d const &o, float radius, float angle, float tscale,
	float thickness=0.0, int tid=-1, colorRGBA const &color=WHITE, bool thick_poly=0)
{
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
	float const r(1.5*radius), q(3.0*radius); // o must be normalized
	float const c(cos(angle)), s(sin(angle)), t(1.0 - c);
	point const p1(r*(t*o.x*o.x + c),     r*(t*o.x*o.y - s*o.z), r*(t*o.x*o.z - s*o.y));
	point const p2(q*(t*o.x*o.z + s*o.y), q*(t*o.y*o.z - s*o.x), q*(t*o.z*o.z + c));

	if (thick_poly) {
		coll_obj cobj;
		cobj.type      = COLL_POLYGON;
		cobj.thickness = thickness;
		cobj.npoints   = 3;
		cobj.cp.color  = color;
		cobj.cp.tid    = tid;
		cobj.cp.tscale = tscale*radius;
		cobj.points[0] = (pos + p1);
		cobj.points[1] = (pos - p1);
		cobj.points[2] = (pos + p2);
		cobj.norm      = get_poly_norm(cobj.points);
		cobj.draw_extruded_polygon(tid, NULL);
	}
	else {
		cross_product(p2, p1).get_norm().do_glNormal();
		float const ts(123.456*radius), tt(654.321*radius);
		if (tscale != 0.0) glTexCoord2f(ts, tt);
		(pos + p1).do_glVertex();
		if (tscale != 0.0) glTexCoord2f(ts+2*tscale*radius, tt);
		(pos - p1).do_glVertex();
		if (tscale != 0.0) glTexCoord2f(ts, tt+2*tscale*radius);
		(pos + p2).do_glVertex();
	}
}


void draw_solid_object_groups() {

	draw_waypoints();
	draw_select_groups(1);
	if (display_mode & 0x0200) d_part_sys.draw();
}


void draw_transparent_object_groups() {

	draw_select_groups(0);
}


void draw_select_groups(int solid) {

	if (!begin_motion) return;
	float const orig_ivo(indir_vert_offset), orig_czb(cobj_z_bias); // store original variable values FIXME: pass into setup_smoke_shaders?
	shader_t s;
	bool const force_tsl(1);
	indir_vert_offset = min(0.1f, indir_vert_offset); // smaller
	cobj_z_bias       = max(0.002f, cobj_z_bias); // larger
	colorRGBA const orig_fog_color(setup_smoke_shaders(s, 0.01, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, force_tsl));
	select_no_texture();
	BLACK.do_glColor();

	for (int i = 0; i < num_groups; ++i) {
		obj_group &objg(obj_groups[i]);

		if (objg.enabled && objg.temperature_ok() && objg.end_id > 0) {
			if ((objg.large_radius() && !(object_types[objg.type].flags & SEMI_TRANSPARENT)) == solid) {
				draw_group(objg, s);
			}
		}
	}
	glDisable(GL_TEXTURE_2D);

	if (s.is_setup()) {
		end_smoke_shaders(s, orig_fog_color);
		indir_vert_offset  = orig_ivo; // restore original variable values
		cobj_z_bias        = orig_czb;
	}
	if (!snow_pld.empty()) { // draw snowflakes from points in a custom geometry shader
		set_specular(0.0, 1.0); // disable
		select_texture(object_types[SNOW].tid, 1, 1);
		check_drawing_flags(object_types[SNOW].flags, 1);
		glDepthMask(GL_FALSE);
		shader_t s;
		s.setup_enabled_lights(2, 1); // VS
		s.set_prefix("#define USE_LIGHT_COLORS", 0); // VS
		s.set_vert_shader("ads_lighting.part*+two_lights_no_xform");
		s.set_frag_shader("simple_texture");
		s.set_geom_shader("pt_billboard_tri", GL_POINTS, GL_TRIANGLE_STRIP, 3);
		//s.set_geom_shader("output_textured_quad.part+pt_billboard", GL_POINTS, GL_TRIANGLE_STRIP, 4);
		s.begin_shader();
		s.add_uniform_float("size", 2.0*object_types[SNOW].radius);
		s.add_uniform_int("tex0", 0);
		s.add_uniform_float("min_alpha", 0.0);
		snow_pld.draw_and_clear();
		s.end_shader();
		glDepthMask(GL_TRUE);
		glDisable(GL_TEXTURE_2D);
		check_drawing_flags(object_types[SNOW].flags, 0);
	}
}



struct wap_obj {

	int id, ndiv;
	wap_obj(int id_, int ndiv_) : id(id_), ndiv(ndiv_) {}
};


void draw_obj(obj_group &objg, vector<wap_obj> *wap_vis_objs, int type, float radius, const colorRGBA &color, int ndiv, int j, bool in_ammo) {

	dwobject const &obj(objg.get_obj(j));
	point const &pos(obj.pos);
	bool const cull_face(get_cull_face(type, color));
	if (cull_face) glEnable(GL_CULL_FACE);
	
	switch (type) {
	case SMILEY:
		if (!(obj.flags & CAMERA_VIEW)) {
			draw_smiley(pos, obj.orientation, radius, ndiv, obj.time, obj.health, j,
				(in_ammo ? NULL : &objg.get_td()->get_mesh(j)));
		}
		break;
	case SFPART:
		draw_smiley_part(pos, pos, obj.orientation, obj.direction, 1, ndiv);
		break;
	case CHUNK:
		draw_chunk(pos, radius, obj.init_dir, obj.vdeform, (obj.flags & TYPE_FLAG), ndiv);
		break;
	case SKULL:
		draw_skull(pos, obj.orientation, radius, obj.status, ndiv);
		break;
	case ROCKET:
		draw_rocket(pos, obj.init_dir, radius, obj.type, ndiv, obj.time);
		break;
	case SEEK_D:
		draw_seekd(pos, obj.init_dir, radius, obj.type, ndiv);
		break;
	case LANDMINE:
		draw_landmine(pos, radius, ndiv, obj.time, obj.source, in_ammo);
		break;
	case PLASMA:
		draw_plasma(pos, pos, radius, (in_ammo ? 1.0 : obj.init_dir.x), ndiv, 1, !in_ammo, obj.time);
		break;
	case GRENADE:
		draw_grenade(pos, obj.init_dir, radius, ndiv, (in_ammo ? 0 : obj.time), in_ammo, 0);
		break;
	case CGRENADE:
		draw_grenade(pos, obj.init_dir, radius, ndiv, (in_ammo ? 0 : obj.time), in_ammo, 1);
		break;
	case BALL:
		// FIXME: this is the only place where drawing an object modifies its physics state, but it's difficult to move the code
		draw_rolling_obj(pos, objg.get_obj(j).init_dir, radius, obj.status, ndiv, ((obj.flags & PLATFORM_COLL) != 0),
			dodgeball_tids[(game_mode == 2) ? (j%NUM_DB_TIDS) : 0], (in_ammo ? NULL : &objg.get_td()->get_matrix(j)));
		break;
	case POWERUP:
	case HEALTH:
	case SHIELD:
		draw_powerup(pos, radius, ndiv, ((type == POWERUP) ? (int)obj.direction : -1), color);
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
	default:
		if (obj.vdeform != all_ones) {
			glPushMatrix();
			translate_to(pos);
			scale_by(obj.vdeform);
			rotate_into_plus_z(obj.orientation); // might be unnesessary
			gluSphere(quadric, radius, ndiv, ndiv);
			glPopMatrix();
		}
		else {
			draw_sphere_dlist(pos, radius, ndiv, 0);
		}
	}
	if (cull_face) glDisable(GL_CULL_FACE);
}


bool is_object_shadowed(dwobject &obj, float cd_scale, float radius) {

	bool is_shadowed((obj.flags & SHADOWED) != 0); // previous value
	float const pt_size(cd_scale/distance_to_camera(obj.pos)); // approx pixel size
	int const skipval(min(20, int(8.0/pt_size)));

	if (skipval <= 1 || (obj.time % skipval) == 0) {
		is_shadowed = pt_is_shadowed(obj.pos, get_specular_light(), radius, obj.coll_id, 0, (pt_size < 2.0));
		if (is_shadowed) obj.flags |= SHADOWED; else obj.flags &= ~SHADOWED;
	}
	return is_shadowed;
}


void draw_group(obj_group &objg, shader_t &s) {

	RESET_TIME;
	set_specular(0.0, 1.0); // disable
	colorRGBA color2;
	set_fill_mode();
	glEnable(GL_NORMALIZE);
	int const type(objg.get_ptype());
	obj_type const &otype(object_types[type]);
	int tid(otype.tid);
	float const radius(otype.radius), cd_scale(NDIV_SCALE*radius*window_width);
	unsigned const flags(otype.flags);
	bool do_texture(select_texture(tid, 1, 1));
	colorRGBA color(otype.color), tcolor(color);
	set_color_alpha(color);
	gluQuadricTexture(quadric, do_texture);
	check_drawing_flags(flags, 1);
	int const clip_level((type == SMILEY || type == LANDMINE || type == ROCKET || type == BALL) ? 2 : 0);
	unsigned num_drawn(0);

	if (type == LEAF) { // leaves
		static vector<pair<unsigned, unsigned> > ordering;
		ordering.resize(0);
		ordering.reserve(objg.end_id);
		float const leaf_size(get_leaf_size());

		for (unsigned j = 0; j < objg.end_id; ++j) {
			dwobject const &obj(objg.get_obj(j));
			if (obj.disabled()) continue;
			float const leaf_scale(2.0*leaf_size*obj.init_dir.z);
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
			set_specular(0.1, 10.0);
			if (s.is_setup()) {s.disable();}
			shader_t ls;
			ls.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
			colorRGBA const orig_fog_color(setup_smoke_shaders(ls, 0.0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1)); // TSL=1
			ls.add_uniform_float("base_color_scale", 0.0); // hack to force usage of material properties instead of color
			ls.add_uniform_float("ambient_scale",    0.0);
			quad_batch_draw qbd;
			set_color(BLACK);

			for (unsigned j = 0; j < ordering.size(); ++j) {
				dwobject const &obj(objg.get_obj(ordering[j].second));
				int const tree_type(ordering[j].first), tid(tree_types[tree_type].leaf_tex);
				float const leaf_scale(2.0*leaf_size*obj.init_dir.z), leaf_x_ar(tree_types[tree_type].leaf_x_ar);
				assert(tid >= 0);
			
				if (draw_model == 0 && tid != last_tid) {
					qbd.draw_and_clear();
					select_texture(tid, 0, 1);
					last_tid = tid;
				}
				point pos(obj.pos);
				if (place_obj_on_grass(pos, leaf_scale)) {pos.z = 0.5*(obj.pos.z + pos.z-leaf_scale);} // leaf is partially on grass
				float const t(((float)obj.time)/((float)otype.lifetime));
				colorRGBA const dry_color(1.0, 0.7, 0.1); // this is the final color, even for partially burnt leaves - oh well
				colorRGBA leaf_color(WHITE);
				UNROLL_3X(leaf_color[i_] *= obj.vdeform[i_];) // vdeform.x is color_scale
				if (leaf_color != BLACK) {blend_color(leaf_color, dry_color, leaf_color, t, 0);}
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
			ls.add_uniform_float("base_color_scale", 1.0);
			ls.add_uniform_float("ambient_scale",    1.0);
			end_smoke_shaders(ls, orig_fog_color);
			if (s.is_setup()) {s.enable();} // back to the original shader
			set_specular(0.0, 1.0);
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
				set_color_by_status(obj.status);
			}
			else if (type != SMILEY && type != SFPART && type != ROCKET && type != CHUNK &&
				type != LANDMINE && type != PLASMA && type != POWERUP && type != HEALTH && type != SHIELD)
			{
				colorRGBA color2(color);
				scale_color_uw(color2, pos);
				set_color_alpha(color2);
			}
			++num_drawn;
			float const pt_size(cd_scale/distance_to_camera(pos));
			int const ndiv(min(N_SPHERE_DIV, max(3, int(min(pt_size, 3.0f*sqrt(pt_size))))));
			draw_obj(objg, wap_vis_objs, type, radius, color, ndiv, j, 0);
		} // for j
		for (unsigned k = 0; k < wap_vis_objs[0].size(); ++k) { // draw weapons
			draw_weapon2(objg.get_obj(wap_vis_objs[0][k].id), radius);
		}
		for (unsigned k = 0; k < wap_vis_objs[1].size(); ++k) { // draw ammo
			unsigned const j(wap_vis_objs[1][k].id);
			assert(j < objg.max_objects());
			draw_ammo(objg, radius, color, wap_vis_objs[1][k].ndiv, j);
		}
		if (!wap_vis_objs[0].empty() || !wap_vis_objs[1].empty()) {
			check_drawing_flags(otype.flags, 1);
			select_texture(tid, 1, 1);
			gluQuadricTexture(quadric, do_texture);
		}
		for (unsigned j = 0; j < 2; ++j) {
			for (unsigned k = 0; k < wap_vis_objs[j].size(); ++k) {
				wap_obj const &wa(wap_vis_objs[j][k]);
				set_obj_specular(flags, 0.5*brightness);
				dwobject const &obj(objg.get_obj(wa.id));
				set_color_alpha(color);
				draw_subdiv_sphere(obj.pos, radius, wa.ndiv, 0, 0);
			}
		}
		if (!smiley_weapons_to_draw.empty()) {
			for (vector<unsigned>::const_iterator i = smiley_weapons_to_draw.begin(); i != smiley_weapons_to_draw.end(); ++i) {
				draw_weapon_in_hand(*i); // Note: view culling doesn't use correct bounding sphere for all weapons
			}
		}
	} // large objects
	else { // small objects
		quad_batch_draw particle_qbd;

		if (type == SHRAPNEL) {
			glBegin(GL_TRIANGLES);
		}
		else if (type == PARTICLE) {
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.01);
		}
		for (unsigned j = 0; j < objg.end_id; ++j) {
			dwobject &obj(objg.get_obj(j));
			point const &pos(obj.pos);
			if (obj.disabled() || (obj.flags & CAMERA_VIEW))      continue;
			float const tradius(obj.get_true_radius()); // differs from radius for fragments
			if (!sphere_in_camera_view(pos, tradius, clip_level)) continue;
			++num_drawn;

			if (type == FRAGMENT) {
				tid = -obj.coll_id - 2; // should we sort fragments by texture id?
				do_texture = select_texture(tid, 1, 1);
				UNROLL_3X(color2[i_] = obj.init_dir[i_];)
				color2.alpha = obj.vdeform.y;
			}
			else if (do_texture && type != SNOW) {
				color2 = WHITE;
			}
			else {
				color2 = object_types[obj.type].color;
			}
			if (type != SHRAPNEL && type != PARTICLE) {
				if (type == DROPLET) select_liquid_color(color2, pos);
				scale_color_uw(color2, pos);

				if (do_texture) {
					assert(tid >= 0);
					tcolor = texture_color(tid);
					tcolor.alpha = 1.0; // don't use texture alpha
					UNROLL_3X(tcolor[i_] *= color2[i_];)
				}
				else {
					tcolor = color2;
				}
			}
			switch (type) {
			case SHELLC:
				set_color_v2(color2, obj.status);
				draw_shell_casing(pos, obj.orientation, obj.init_dir, tradius, obj.angle, cd_scale, obj.direction);
				break;
			case STAR5:
				set_color_v2(color2, obj.status);
				draw_star(pos, obj.orientation, obj.init_dir, tradius, obj.angle, 1);
				break;
			case SHRAPNEL:
				set_emissive_color_obj(get_glow_color(obj, 1));
				draw_rotated_triangle(obj.pos, obj.orientation, tradius, obj.angle, 0.0);
				break;
			case PARTICLE:
				{
					colorRGBA const glow_color(get_glow_color(obj, 0));
					if (glow_color.alpha > 0.0) {particle_qbd.add_billboard(obj.pos, get_camera_pos(), up_vector, glow_color, 1.2*tradius, 1.2*tradius);}
				}
				break;

			case SAND:
			case DIRT:
			case ROCK:
				color2 *= obj.orientation.y;
				if (do_texture) tcolor *= obj.orientation.y;
				draw_sized_point(obj, tradius, obj.orientation.x*cd_scale, color2, tcolor, do_texture, (type == DIRT || type == ROCK));
				break;

			case FRAGMENT: // draw_fragment()?
				if (obj.vdeform.z > 0.0) { // shatterable - use triangle
					set_color_v2(color2, obj.status);
					bool const use_thick(tid < 0); // when not textured
					if (!use_thick) glBegin(GL_TRIANGLES); // Note: needs 2-sided lighting
					draw_rotated_triangle(pos, obj.orientation, tradius, obj.angle, (do_texture ? obj.vdeform.z : 0.0), 0.2*tradius, tid, color2, use_thick); // obj.vdeform.z = tscale
					if (!use_thick) glEnd();
					break;
				}
				draw_sized_point(obj, tradius, cd_scale, color2, tcolor, do_texture, 2);
				break;

			default:
				if (DEBUG_COLOR_COLLS) {
					int cindex;
					float const time(TIMESTEP*fticks);
					point const pos2(pos + obj.velocity*time - point(0.0, 0.0, -base_gravity*GRAVITY*time*time*otype.gravity));
					set_color(check_coll_line(pos, pos2, cindex, -1, 0, 0) ? RED : GREEN);
					draw_line(pos, pos2);
				}
				draw_sized_point(obj, tradius, cd_scale, color2, tcolor, do_texture, 0);
			} // switch (type)
		} // for j
		if (!puddles.empty()) { // draw puddles
			glDepthMask(GL_FALSE);
			select_texture(BLUR_TEX);
			plus_z.do_glNormal();
			glBegin(GL_TRIANGLES);

			for (vector<puddle_t>::const_iterator p = puddles.begin(); p != puddles.end(); ++p) {
				set_color_alpha(p->color);
				draw_billboard(p->pos, (p->pos + plus_z), plus_x, 5.0*p->radius, 5.0*p->radius);
			}
			glEnd();
			glDepthMask(GL_TRUE);
			select_no_texture();
			puddles.resize(0);
		}
		if (type == SHRAPNEL) {
			clear_emissive_color();
			glEnd();
		}
		else if (type == PARTICLE) {
			particle_qbd.draw();
			glDisable(GL_ALPHA_TEST);
		}
		if (!obj_pld.empty()) {
			glEnable(GL_COLOR_MATERIAL); // unnecessary?
			select_no_texture();
			if (s.is_setup()) {s.add_uniform_float("base_color_scale", 0.0);} // hack to force usage of material properties instead of color
			obj_pld.draw_and_clear();
			if (s.is_setup()) {s.add_uniform_float("base_color_scale", 1.0);}
			glDisable(GL_COLOR_MATERIAL);
		}
	} // small object
	check_drawing_flags(flags, 0);
	gluQuadricTexture(quadric, GL_FALSE);
	select_no_texture();

	if (SHOW_DRAW_TIME) {
		cout << "type = " << objg.type << ", num = " << objg.end_id << ", drawn = " << num_drawn << " ";
		PRINT_TIME("Group");
	}
}


void draw_sized_point(dwobject &obj, float radius, float cd_scale, const colorRGBA &color, const colorRGBA &tcolor,
					  bool do_texture, int is_chunky)
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
		puddles.push_back(puddle_t(pos, radius, color2));
		return;
	}
	if (draw_snowflake) { // draw as a point to be converted to a billboard by the geometry shader
		bool const is_shadowed(is_object_shadowed(obj, cd_scale, radius));
		snow_pld.add_pt(pos, (is_shadowed ? zero_vector : (get_light_pos() - pos)), (do_texture ? tcolor : color));
		return;
	}
	if (!draw_large) { // draw as a point
		bool const scatters(type == RAIN || type == SNOW);
		vector3d const n((scatters ? get_light_pos() : camera) - pos);
		obj_pld.add_pt(pos, n, (do_texture ? tcolor : color));
		return;
	}
	set_color_v2(color, obj.status);
	bool const cull_face(get_cull_face(type, color));
	glPushMatrix();

	if (cull_face) {
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
	}

	// draw as a sphere
	if (is_chunky) {
		assert(!tail);
		vector3d const v((is_chunky == 2) ? obj.orientation : obj.init_dir);
		int const ndiv(max(3, int(3 + 1.5*(v.x + v.y + v.z))));
		translate_to(pos);
		vector3d const scale((0.8+0.5*fabs(v.x)), (0.8+0.5*fabs(v.y)), (0.8+0.5*fabs(v.z)));
		scale_by(scale*radius);
		glRotatef(360.0*(v.x - v.y), v.x, v.y, (v.z+0.01));
		draw_unit_sphere(ndiv, do_texture);
		glTranslatef(0.1*(v.x-v.y), 0.1*(v.y-v.z), 0.1*(v.x-v.z));
		glRotatef(360.0*(v.z - v.x), v.y, v.z, (v.x+0.01));
		draw_unit_sphere(ndiv, do_texture);
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
	
		if (quadric != 0 && ndiv > 3 && tail) { // cone on the tail of the raindrop
			translate_to(pos);
			gluCylinder(quadric, radius, 0.0, 2.5*radius, (ndiv>>1), 1);
			glTranslatef(0.0, 0.0, -0.6*radius);
			glScalef(1.0, 1.0, 2.0);
			pos = all_zeros;
		}
		draw_sphere_dlist(pos, radius, ndiv, do_texture);
	}
	if (cull_face) glDisable(GL_CULL_FACE);
	glPopMatrix();
}


void draw_weapon2(dwobject const &obj, float radius) {

	int const wid((int)obj.direction);
	vector3d d(-1.0, 0.0, 0.0), v(obj.pos);
	if (wid != W_BLADE) v.x += 0.7*radius;
	draw_weapon_simple(v, d, 0.5*radius, obj.coll_id, wid, 0.25);
}


void draw_ammo(obj_group &objg, float radius, const colorRGBA &color, int ndiv, int j) {

	dwobject const &obj(objg.get_obj(j));
	point pos(obj.pos);
	vector<wap_obj> wap_vis_objs[2]; // not actually used
	int const atype(get_ammo_or_obj((int)obj.direction));

	if (atype >= 0) {
		check_drawing_flags(object_types[atype].flags, 1);
		int const textured(select_texture(object_types[atype].tid, 1, 1));
		gluQuadricTexture(quadric, textured);
		set_color_alpha(object_types[atype].color);
		bool const cull_face(get_cull_face(atype, color));
		if (cull_face) glEnable(GL_CULL_FACE);

		switch (atype) {
		case SHELLC: // M16
			pos.z -= 0.5*radius;
			draw_cylinder(pos, 1.0*radius, 0.2*radius, 0.2*radius, ndiv, 1, 1);
			pos.z += radius;
			draw_subdiv_sphere(pos, 0.2*radius, ndiv, 0, 0);
			break;
		case PROJECTILE: // shotgun
			for (unsigned n = 0; n < 2; ++n) { // two shells in one ammo
				point pos2(pos);
				pos2.x += (1.0 - 2.0*n)*0.3*radius;
				set_color_alpha(RED);
				pos2.z -= 0.5*radius;
				draw_cylinder(pos2, 1.2*radius, 0.3*radius, 0.3*radius, ndiv, 1, 1);
				set_color_alpha(GOLD);
				pos2.z -= 0.2*radius;
				draw_cylinder(pos2, 0.4*radius, 0.32*radius, 0.32*radius, ndiv, 1, 1);
			}
			break;
		case BEAM: // laser
			set_color_alpha(RED);
			pos.z -= 0.5*radius;
			draw_cylinder(pos, 1.0*radius, 0.1*radius, 0.1*radius, ndiv, 1, 1);
			break;
		case STAR5: // throwing star
			draw_star(pos, obj.orientation, obj.init_dir, 0.4*radius, obj.angle, 0);
			break;
		case GASSED:
			draw_subdiv_sphere(pos, 0.6*radius, ndiv, 1, 0);
			break;
		default:
			draw_obj(objg, wap_vis_objs, atype, 0.4*radius, color, ndiv, j, 1);
		}
		if (cull_face) glDisable(GL_CULL_FACE);
	}
}


colorRGBA get_powerup_color(int powerup) {

	switch (powerup) {
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

	glRotatef(atan2(dir.y, dir.x)*TO_DEG-90.0, 0.0, 0.0, 1.0);
	glRotatef(safe_acosf(-dir.z)*TO_DEG-90.0, 1.0, 0.0, 0.0);
}


void draw_smiley_part(point const &pos, point const &pos0, vector3d const &orient, int type, int use_orient, int ndiv, float scale) {

	assert(type < NUM_SMILEY_PARTS);
	float const radius(scale*object_types[SFPART].radius);
	colorRGBA const sf_color[NUM_SMILEY_PARTS] = {BLACK, RED, PINK};
	set_color_alpha(sf_color[type]);

	switch (type) {
	case SF_EYE:
		draw_sphere_dlist(pos, 1.0*radius, ndiv, 0);
		break;
	case SF_NOSE:
		draw_sphere_dlist(pos, 1.2*radius, ndiv, 0);
		break;
	case SF_TONGUE:
		glPushMatrix();
		translate_to(pos);
		if (use_orient) rotate_to_dir(orient);
		glScalef(1.0, 2.0, 0.25);
		draw_sphere_dlist(all_zeros, 1.5*radius, ndiv, 0);
		glPopMatrix();
		break;
	default: assert(0);
	}
}


colorRGBA mult_alpha(colorRGBA const &c, float alpha) {
	return colorRGBA(c.R, c.G, c.B, c.A*alpha);
}


void draw_smiley(point const &pos, vector3d const &orient, float radius, int ndiv, int time,
				 float health, int id, mesh2d const *const mesh)
{
	colorRGBA color;
	int const powerup(sstates[id].powerup), ndiv2(max(3, (ndiv>>1)));
	float const dist(distance_to_camera(pos));
	glPushMatrix();
	translate_to(pos);
	rotate_to_dir(orient);
	//uniform_scale(radius);
	point pos2(-0.4*radius, 0.85*radius, 0.3*radius);

	// draw eyes
	for (unsigned i = 0; i < 2; ++i) {
		if (health > 10.0) {
			float const scale((powerup == PU_SPEED) ? 1.5 : 1.0);
			point const pos3 ((powerup == PU_SPEED) ? (pos2 + point(0.0, 0.2*radius, 0.0)) : pos2);
			draw_smiley_part(pos3, pos, orient, SF_EYE, 0, ndiv2, scale); // eyes
		}
		else {
			set_color_alpha(BLACK);
			enable_blend();
			glLineWidth(min(8.0f, max(1.0f, 6.0f/dist)));
			glPushMatrix();
			translate_to(pos2);
			uniform_scale(radius);
			glBegin(GL_LINES);

			for (unsigned l = 0; l < 2; ++l) {
				for (unsigned p = 0; p < 2; ++p) {
					glVertex3f((p ? -0.12 : 0.12), 0.1, ((p^l) ? -0.12 : 0.12));
				}
			}
			glEnd();
			glPopMatrix();
			glLineWidth(1.0);
			disable_blend();
		}
		pos2.x *= -1.0;
	}

	// draw nose
	if (powerup != PU_INVISIBILITY || same_team(id, -1)) { // show nose even if invisible if same team as player
		point pos3(0.0, 1.1*radius, 0.0);
		draw_smiley_part(pos3, pos, orient, SF_NOSE, 0, ndiv2); // nose
	}
	float alpha(1.0);

	switch (powerup) {
		case PU_DAMAGE: // devil horns
			set_color_alpha(RED);
			glPushMatrix();
			glTranslatef( 0.3*radius, 0.7*radius, 0.6*radius);
			draw_cylinder(0.6*radius, 0.1*radius, 0.0, ndiv2, 1, 0);
			glTranslatef(-0.6*radius, 0.0, 0.0);
			draw_cylinder(0.6*radius, 0.1*radius, 0.0, ndiv2, 1, 0);
			glPopMatrix();
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
			set_color_alpha(BLACK);
			glPushMatrix();
			glTranslatef(0.0, 0.0, 0.9*radius);
			draw_cylinder(0.5*radius, 0.05*radius, 0.05*radius, ndiv2, 1, 0);
			glTranslatef(0.0, 0.0, 0.5*radius);
			glRotatef(float((30*time)%360), 0.0, 0.0, 1.0);
			glScalef(1.0, 0.25, 0.05);
			glTranslatef( 0.5*radius, 0.0, 0.0);
			draw_sphere_dlist(all_zeros, 0.5*radius, ndiv, 0); // propeller
			glTranslatef(-1.0*radius, 0.0, 0.0);
			draw_sphere_dlist(all_zeros, 0.5*radius, ndiv, 0); // propeller
			glPopMatrix();
			break;

		case PU_INVISIBILITY:
			{
				float const put(float(sstates[id].powerup_time)/TICKS_PER_SECOND), init_put(float(POWERUP_TIME)/TICKS_PER_SECOND);

				if ((init_put - put) < 1.0) { // fading out
					alpha = (1.0 - (init_put - put));
				}
				else if (put < 1.0) { // fading in
					alpha = (1.0 - put);
				}
				else { // fully invisible
					glPopMatrix();
					return;
				}
				if (alpha < 1.0) enable_blend();
				break;
			}
	} // switch (powerup)
	if (powerup != PU_INVISIBILITY && time%10 < 5 && powerup >= 0) { // powerup
		set_color_alpha(get_powerup_color(powerup));
	}
	else if (health >= 50.0) {
		set_color_alpha(mult_alpha(YELLOW, alpha));
	}
	else {
		set_color_alpha(colorRGBA(1.0, (0.25 + 0.015*health), 0.0, alpha));
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
		draw_sphere_dlist(all_zeros, radius, ndiv, 1);
	}
	select_no_texture();

	if (teams > 1) { // draw team headband
		set_color_alpha(mult_alpha(get_smiley_team_color(id), alpha));
		glPushMatrix();
		glTranslatef(0.0, 0.0, 0.45*radius);
		glScalef(1.0, 1.0, 0.5);
		gluSphere(quadric, 0.94*radius, ndiv, ndiv2);
		glPopMatrix();
	}
	// draw unique identifier
	int temp(teams);
	teams = 10; // 10 default colors
	set_color_alpha(mult_alpha(get_smiley_team_color(id+1), alpha));
	teams = temp;
	glPushMatrix();
	glTranslatef(0.0, 0.0, 0.8*radius);
	glScalef(1.0, 1.0, 0.3);
	gluSphere(quadric, 0.65*radius, ndiv, ndiv2);
	glPopMatrix();
	
	// draw mouth
	float const hval(0.004*(100.0 - min(160.0f, health)));
	set_color_alpha(mult_alpha(BLACK, alpha));
	enable_blend();
	glLineWidth(min(8.0f, max(1.0f, 5.0f/dist)));
	glPushMatrix();
	uniform_scale(radius);
	glBegin(GL_LINE_STRIP);
	glVertex3f(-0.5,  0.95, -0.2-hval);
	glVertex3f(-0.15, 0.95, -0.4);
	glVertex3f( 0.15, 0.95, -0.4);
	glVertex3f( 0.5,  0.95, -0.2-hval);
	glEnd();
	glPopMatrix();
	glLineWidth(1.0);
	disable_blend();

	// draw tongue
	if (sstates[id].kill_time < int(2*TICKS_PER_SECOND) || powerup == PU_DAMAGE) { // stick your tongue out at a dead enemy
		point pos4(0.0, 0.8*radius, -0.4*radius);
		draw_smiley_part(pos4, pos, orient, SF_TONGUE, 0, ndiv2);
	}
	if (game_mode == 2 && (sstates[id].p_ammo[W_BALL] > 0 || UNLIMITED_WEAPONS)) { // dodgeball
		select_texture(select_dodgeball_texture(id), 1, 1);
		set_color_alpha(mult_alpha(object_types[BALL].color, alpha));
		draw_sphere_dlist(point(0.0, 1.3*radius, 0.0), 0.8*object_types[BALL].radius, ndiv, 1);
		select_no_texture();
	}
	glPopMatrix();
	vector3d hit_dir;
	int const hit(get_smiley_hit(hit_dir, id));

	if (hit > 0) { // hit - draw damage or shields
		select_texture(SBLUR_TEX);
		enable_blend();
		colorRGBA color2((sstates[id].shields < 0.01) ? BLOOD_C : GREEN);
		color2.alpha = alpha*hit/6.0;
		set_color_alpha(color2);
		glPushMatrix();
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.05);
		translate_to(pos);
		rotate_sphere_tex_to_dir(hit_dir);
		draw_sphere_dlist(all_zeros, 1.015*radius, ndiv, 1);
		glDisable(GL_ALPHA_TEST);
		glPopMatrix();

		if (powerup == PU_SHIELD) {
			select_texture(PLASMA_TEX);
			set_color_alpha(PURPLE, alpha); // not scaled
			draw_sphere_dlist(pos, 1.05*radius, ndiv, 1);
		}
		disable_blend();
		select_no_texture();
	}
	if (alpha < 1.0) disable_blend();
}


void draw_powerup(point const &pos, float radius, int ndiv, int type, const colorRGBA &color) {

	set_emissive_color_obj((type == -1) ? color : get_powerup_color(type));
	draw_subdiv_sphere(pos, 0.7*radius, ndiv, 0, 0); // draw flare/billboard?
	clear_emissive_color();
	set_color_alpha(color);
	draw_subdiv_sphere(pos, radius, ndiv, 0, 0);
}


void draw_rolling_obj(point const &pos, point &lpos, float radius, int status, int ndiv, bool on_platform, int tid, xform_matrix *matrix) {

	select_texture(tid, 1, 1);
	glPushMatrix();
	translate_to(pos);
	
	if (matrix) {
		if (on_platform) lpos = pos; // reset now so there's no rotation
		apply_obj_mesh_roll(*matrix, pos, lpos, radius, ((status == 1) ? 0.01 : 0.0), ((status == 1) ? 0.2 : 1.0));
	}
	draw_sphere_dlist(all_zeros, radius, 2*ndiv, 1);
	glPopMatrix();
	lpos = pos;
}


void draw_skull(point const &pos, vector3d const &orient, float radius, int status, int ndiv) {

	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.9);
	glPushMatrix();
	translate_to(pos);
	rotate_from_v2v(orient, vector3d(0.0, -1.0, 0.0));
	glRotatef(180.0, 0.0, 1.0, 0.0);
	draw_sphere_dlist_back_to_front(all_zeros, radius, 2*ndiv, 1);
	glPopMatrix();
	glDisable(GL_ALPHA_TEST);
}


void draw_rocket(point const &pos, vector3d const &orient, float radius, int type, int ndiv, int time) {

	glPushMatrix();
	translate_to(pos);
	rotate_by_vector(orient, 0.0);
	set_color_alpha(RED);
	uniform_scale(radius);
	glBegin(GL_TRIANGLES);
	glVertex3f( 0.0,  0.0,  0.0);
	glVertex3f( 1.8,  0.0, -2.0);
	glVertex3f(-1.8,  0.0, -2.0);
	glVertex3f( 0.0,  0.0,  0.0);
	glVertex3f( 0.0,  1.8, -2.0);
	glVertex3f( 0.0, -1.8, -2.0);
	glEnd();
	set_color_alpha(object_types[ROCKET].color);
	glScalef(1.0, 1.0, -2.0);
	draw_unit_sphere(ndiv, 0);
	gluCylinder(quadric, 1.0, 1.0, 1.1, ndiv, 1);
	glPopMatrix();
	if (type == ROCKET) gen_rocket_smoke(pos, orient, radius);
}


void draw_seekd(point const &pos, vector3d const &orient, float radius, int type, int ndiv) {

	assert(quadric);
	glPushMatrix();
	translate_to(pos);
	rotate_by_vector(orient, 0.0);
	uniform_scale(radius);
	glPushMatrix();
	glTranslatef(0.0, 0.0, -2.0);
	set_color_alpha(BLACK);
	gluCylinder(quadric, 1.0, 0.0, 1.0, ndiv, 1);
	glTranslatef(0.0, 0.0, -0.25);
	gluCylinder(quadric, 1.0, 1.0, 0.25, ndiv, 1);
	glPopMatrix();

	glScalef(1.0, 1.0, 1.5);
	glRotatef(90.0, -1.0, 0.0, 0.0);
	glRotatef(90.0,  0.0, 1.0, 0.0);
	set_color_alpha(WHITE);
	select_texture(SKULL_TEX);
	draw_unit_sphere(ndiv, 1);
	select_no_texture();

	glPopMatrix();
	if (type == SEEK_D) gen_rocket_smoke(pos, orient, radius);
}


colorRGBA const &get_landmine_light_color(int time) {

	return ((time < 40) ? GREEN : (((time/6)&1) ? RED : BLUE));
}


float get_landmine_sensor_height(float radius, int time) {

	return ((time <= 6) ? 0 : ((time > 16) ? 1.5*radius : (radius + 0.05*radius*(time - 6))));
}


void draw_landmine(point pos, float radius, int ndiv, int time, int source, bool in_ammo) {

	assert(radius > 0.0 && ndiv > 0);

	if (!in_ammo) {
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

		if (!point_outside_mesh(xpos, ypos) && !is_mesh_disabled(xpos, ypos) &&
			pos.z < (interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 1) + 1.05*radius))
		{
			pos.z -= 0.8*radius; // appears to sink into the ground
		}
	}
	if (!DEBUG_COLORCODE) set_color_alpha(WHITE);
	draw_subdiv_sphere(pos, radius, ndiv, 1, 0); // main body
	select_no_texture();
	glPushMatrix();
	translate_to(pos);
	float const val(get_landmine_sensor_height(radius, time));

	if (time > 6) {
		set_color_alpha(GRAY);
		gluCylinder(quadric, 0.05*radius, 0.05*radius, val, ndiv, 1);
		if (teams > 1) set_color_alpha(get_smiley_team_color(source)); // use team color
		gluDisk(quadric, 0, 0.05*radius, ndiv, 1); // sensor
	}
	pos.z += val;

	if (time > 20) {
		float val;
		glTranslatef(0.0, 0.0, 1.4*radius);
		glRotatef(((time > 26) ? 90.0 : 15.0*(time - 20)), 1.0, 0.0, 0.0);
		
		if (time > 36) {
			glRotatef(10.0*(time%36), 0.0, 1.0, 0.0);
			val = 0.42*radius;
		}
		else if (time > 26){
			val = 0.04*radius*(time - 26) + 0.02*radius;
		}
		if (time > 26) gluCylinder(quadric, 0.0, val, 0.4*radius, ndiv, 1); // sensor pole
	}
	glPopMatrix();

	if (time > 5) {
		pos.z += 0.15*radius;
		set_emissive_color_obj(get_landmine_light_color(time));
		draw_subdiv_sphere(pos, 0.15*radius, ndiv/2, 0, 0); // warning light
		clear_emissive_color();
	}
	select_texture(object_types[LANDMINE].tid, 1, 1);
}


colorRGBA get_plasma_color(float size) {

	return colorRGBA(1.0, size/5.0, max(0.0f, 0.5f*(size-3.0f)), 0.9);
}


void draw_plasma(point const &pos, point const &part_pos, float radius, float size, int ndiv, int shpere_tex, bool gen_parts, int time) {

	int const tmode(shpere_tex ? GL_SPHERE_MAP : GL_EYE_LINEAR);
	colorRGBA const color(get_plasma_color(size + 0.5*(0.5 + 0.16*abs((time % 12) - 6))));

	if (animate2) {
		setup_texgen(0.2*rand_uniform(0.95, 1.05)/radius, 0.2*rand_uniform(0.95, 1.05)/radius, rand_float(), rand_float(), 0.0, tmode);
	}
	else {
		setup_texgen(0.2/radius, 0.2/radius, 0.0, 0.0, 0.0, tmode);
	}
	set_emissive_color_obj(color);
	if (animate2) radius *= rand_uniform(0.99, 1.01) + 0.1*(0.5 + 0.1*(abs((time % 20) - 10)));
	draw_sphere_dlist(pos, size*radius, ndiv, 1);
	clear_emissive_color();
	disable_texgen();
	if (gen_parts && animate2 && !is_underwater(part_pos, 1) && (rand()&15) == 0) gen_particles(part_pos, 1);
}


void draw_chunk(point const &pos, float radius, vector3d const &v, vector3d const &vdeform, int charred, int ndiv) {

	ndiv    = min(ndiv, max(3, int(3 + 1.5*(v.x + v.y + v.z))));
	radius *= (0.5 + fabs(v.x));
	glPushMatrix();
	translate_to(pos);
	vector3d scale((0.8+0.5*fabs(v.x)), (0.8+0.5*fabs(v.y)), (0.8+0.5*fabs(v.z)));
	scale *= radius;
	
	if (vdeform != all_ones) { // apply deformation
		float vdmin(1.0);
		UNROLL_3X(scale[i_] *= vdeform[i_]; vdmin = min(vdmin, vdeform[i_]);)
		if (vdmin < 1.0) scale *= pow(1.0/vdmin, 1.0/3.0);
	}
	scale_by(scale);
	glRotatef(360.0*(v.x - v.y), v.x, v.y, (v.z+0.01));
	set_color_alpha((charred ? BLACK : YELLOW));
	draw_unit_sphere(ndiv, 0);
	set_color_alpha((charred ? DK_GRAY : BLOOD_C));
	glTranslatef(0.1*(v.x-v.y), 0.1*(v.y-v.z), 0.1*(v.x-v.z));
	glRotatef(360.0*(v.z - v.x), v.y, v.z, (v.x+0.01));
	draw_unit_sphere(ndiv, 0);
	glPopMatrix();
}


void draw_grenade(point const &pos, vector3d const &orient, float radius, int ndiv, int time, bool in_ammo, bool is_cgrenade) {

	assert(quadric);
	glPushMatrix();
	translate_to(pos);
	uniform_scale(radius);
	glPushMatrix();
	if (!is_cgrenade) glScalef(0.8, 0.8, 1.2); // rotate also?
	set_color_alpha(BLACK);
	draw_unit_sphere(ndiv, 0);
	glPopMatrix();

	float const stime(1.0 - float(time)/float(object_types[is_cgrenade ? CGRENADE : GRENADE].lifetime)), sval(0.2 + 0.8*stime);
	vector3d const vr((orient.x == 0.0 && orient.y == 0.0) ? vector3d(1.0, 0.0, 0.0) : vector3d(orient.x, orient.y, 0.0));
	vector3d vd(plus_z);
	rotate_vector3d_norm(vr, -0.25*PI, vd);
	rotate_about(45.0, vr);
	glTranslatef(0.0, 0.0, 0.7);
	glDisable(GL_CULL_FACE);
	gluCylinder(quadric, 0.3, 0.3, 0.5, max(3, ndiv/2), 1);
	set_color_alpha(GRAY);
	glTranslatef(0.0, 0.0, 0.3);
	gluCylinder(quadric, 0.05, 0.05, sval, max(3, ndiv/4), 1); // fuse
	glPopMatrix();

	if (!animate2) return;
	point const spos(pos + vd*((1.0 + sval)*radius));
	colorRGBA scolor;
	blend_color(scolor, YELLOW, ORANGE, rand_uniform(0.3, 0.7), 1);
	float const size(radius*rand_uniform(0.5, 0.7));
	sparks.push_back(spark_t(spos, scolor, size));
	add_dynamic_light(0.15, spos, scolor); // out of sync by a frame?
	if (!in_ammo && (rand()&15) == 0) gen_particles(spos, 1, 0.5, 1);
}


void draw_star(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate) { // not all variables used
	
	glPushMatrix();
	translate_to(pos);
	uniform_scale(2.0*radius);
	up_norm.do_glNormal();

	if (rotate) {
		rotate_by_vector(init_dir, -90.0);
		if (angle != 0.0) rotate_about(angle, orient);
	}
	glBegin(GL_TRIANGLES); // Note: needs 2-sided lighting

	for (int i = N_STAR_POINTS-1; i >= 0; --i) {
		int ii((i == 0) ? (N_STAR_POINTS<<1)-1 : (i<<1)-1);
		star_pts[ii].do_glVertex();
		ii = (i << 1);
		star_pts[ii].do_glVertex();
		++ii;
		star_pts[ii].do_glVertex();
	}
	glEnd();
	glPopMatrix();
}


void draw_shell_casing(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius,
					   float angle, float cd_scale, unsigned char type)
{
	float const point_size(cd_scale/distance_to_camera(pos));
	int const ndiv(max(3, min(N_SPHERE_DIV/2, int(point_size))));
	//glDepthMask(1);
	glPushMatrix();
	translate_to(pos);
	glRotatef(TO_DEG*init_dir.x, 0.0, 0.0, 1.0);
	//rotate_by_vector(init_dir, 0.0);
	rotate_about(angle, orient);
	uniform_scale(radius); // Note: needs 2-sided lighting

	if (type == 0) { // M16 shell casing
		gluCylinder(quadric, 1.0, 1.0, 4.0, ndiv, 1);
		if (point_size > 1.0) gluDisk(quadric, 0, 1.0, ndiv, 1);
	}
	else if (type == 1) { // shotgun shell casing
		set_color_alpha(RED);
		glTranslatef(0.0, 0.0, -2.0);
		gluCylinder(quadric, 1.2, 1.2, 4.8, ndiv, 1);
		set_color_alpha(GOLD);
		glTranslatef(0.0, 0.0, -0.8);
		gluCylinder(quadric, 1.28, 1.28, 1.6, ndiv, 1);
		if (point_size > 1.0) gluDisk(quadric, 0, 1.28, ndiv, 1);
	}
	else {
		assert(0);
	}
	glPopMatrix();
	//glDepthMask(0);
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
