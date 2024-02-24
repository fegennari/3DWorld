// 3D World - Shadow Mapping using Shaders
// by Frank Gennari
// 1/21/11
#include "3DWorld.h"
#include "collision_detect.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "model3d.h"
#include "shadow_map.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>


unsigned const NUM_CSM_CASCADES = 4; // must agree with the values in csm_layers.geom
unsigned const MAT4x4_SIZE      = sizeof(glm::mat4);

// texture storage datatypes for local shadow maps
#if 1 // integer
//int const SHADOW_MAP_DATATYPE = GL_UNSIGNED_BYTE; // 8-bit shadow maps
int const SHADOW_MAP_DATATYPE = GL_UNSIGNED_SHORT; // 16-bit shadow maps
//int const SHADOW_MAP_DATATYPE = GL_UNSIGNED_INT; // 32-bit shadow maps (overkill)
int const SMAP_INTERNAL_FORMAT  = GL_DEPTH_COMPONENT;
unsigned const smap_bytes_per_pixel = sizeof(unsigned short);
#else // floating-point
int const SHADOW_MAP_DATATYPE  = GL_FLOAT; // floating-point shadow maps
int const SMAP_INTERNAL_FORMAT = GL_DEPTH_COMPONENT32F;
unsigned const smap_bytes_per_pixel = sizeof(float);
#endif

bool voxel_shadows_updated(0);
unsigned shadow_map_sz(0), scene_smap_vbo_invalid(0), empty_smap_tid(0);
pos_dir_up orig_camera_pdu;
ubo_wrap_t shadow_matrix_ubo; // reused for dlights shadow map matrices

extern bool snow_shadows, enable_depth_clamp, flashlight_on, interior_shadow_maps, enable_ground_csm;
extern int window_width, window_height, animate2, display_mode, tree_mode, ground_effects_level, num_trees, camera_coll_id;
extern unsigned enabled_lights;
extern float NEAR_CLIP, tree_deadness, vegetation;
extern vector<shadow_sphere> shadow_objs;
extern set<unsigned> moving_cobjs;
extern coll_obj_group coll_objects;
extern cobj_draw_groups cdraw_groups;
extern platform_cont platforms;

void draw_trees(bool shadow_only=0, bool reflection_pass=0);
void free_light_source_gl_state();
void set_shadow_tex_params(unsigned &tid, bool is_array, bool is_csm, bool use_white_border=0);


cube_t get_scene_bounds() {
	return cube_t(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, min(zbottom, czmin), max(ztop, czmax));
}

float approx_pixel_width(unsigned smap_sz) {
	return 0.5*sqrt(X_SCENE_SIZE*X_SCENE_SIZE + Y_SCENE_SIZE*Y_SCENE_SIZE) / smap_sz;
}
int get_smap_ndiv(float radius, unsigned smap_sz) {
	// dynamic based on distance(camera, line(lpos, scene_center))?
	return min(N_SPHERE_DIV, max(4, int(0.5*radius/approx_pixel_width(smap_sz))));
}
int get_def_smap_ndiv(float radius) {return get_smap_ndiv(radius, shadow_map_sz);}


struct ground_mode_smap_data_t : public cached_dynamic_smap_data_t { // used for ground mode sun and moon directional lights

	ground_mode_smap_data_t(unsigned tu_id_) : cached_dynamic_smap_data_t(tu_id_, shadow_map_sz) {}
	virtual void render_scene_shadow_pass(point const &lpos);
	virtual bool needs_update(point const &lpos);
};

vector<ground_mode_smap_data_t> smap_data;

void ensure_smap_data() { // sun/moon

	if (smap_data.empty()) {
		for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
			ground_mode_smap_data_t smd(GLOBAL_SMAP_START_TU_ID+l); // tu_ids 6 and 7
			smd.is_csm = enable_ground_csm;
			smap_data.push_back(smd);
		}
	}
	assert(smap_data.size() == NUM_LIGHT_SRC);
}


class smap_vertex_cache_t : public vbo_wrap_t { // used for sun and moon shadows of dynamic objects in ground mode

	unsigned num_verts1, num_verts2;
	vector<vert_wrap_t> dverts;
	vector<unsigned> movable_cids;

	void end_block1(unsigned size) {num_verts1 = size;}

	void upload(vector<vert_wrap_t> const &verts) {
		if (!verts.empty()) {
			create_and_upload(verts);
			assert(vbo_valid());
		}
		num_verts2 = verts.size();
	}
	int get_ndiv(coll_obj const &c, unsigned smap_sz, unsigned fixed_ndiv) {
		if (fixed_ndiv) return fixed_ndiv;
		if (c.type == COLL_SPHERE) {return get_smap_ndiv(c.radius, smap_sz);}
		if (c.type == COLL_CYLINDER || c.type == COLL_CYLINDER_ROT || c.type == COLL_CAPSULE || c.type == COLL_TORUS) {return get_smap_ndiv(max(c.radius, c.radius2), smap_sz);}
		return 0; // ndiv is unused
	}
public:
	void render() const {
		if (num_verts2 == 0) return; // empty
		shader_t s;
		s.begin_shadow_map_shader();
		assert(num_verts1 <= num_verts2);
		pre_render();
		if (num_verts1 > 0) {draw_verts<vert_wrap_t>(NULL, num_verts1, GL_TRIANGLES);}

		if (num_verts2 > num_verts1) {
			glEnable(GL_CULL_FACE);
			//glCullFace(GL_FRONT); // faster, but artifacts at surface intersections such as stairs, especially with more z bias
			draw_verts<vert_wrap_t>(NULL, (num_verts2 - num_verts1), GL_TRIANGLES, num_verts1);
			//glCullFace(GL_BACK);
			glDisable(GL_CULL_FACE);
		}
		post_render();
		s.end_shader();
	}
	void render_dynamic() {
		draw_and_clear_verts(dverts, GL_TRIANGLES);
	}
	smap_vertex_cache_t() : num_verts1(0), num_verts2(0) {}

	void free() {
		clear_vbo();
		num_verts1 = num_verts2 = 0;
		dverts.clear();
		movable_cids.clear();
	}

	void add_cobjs(unsigned smap_sz, unsigned fixed_ndiv, bool enable_vfc) {
		if (coll_objects.drawn_ids.empty()) return; // do nothing
		if (vbo_valid()) return; // already valid
		// only valid if drawing trees, small trees, and scenery separately
		vector<vert_wrap_t> verts;
		vector<pair<float, unsigned> > z_sorted;
		movable_cids.clear();

		for (cobj_id_set_t::const_iterator i = coll_objects.drawn_ids.begin(); i != coll_objects.drawn_ids.end(); ++i) {
			coll_obj const &c(coll_objects.get_cobj(*i));
			assert(c.cp.draw);
			if (c.no_shadow_map()) continue;
			// Note: since these are static/drawn once, we can't do any VFC for the camera, but we can do VFC for the light frustum of local light sources 
			if (enable_vfc && !c.is_cobj_visible()) continue; // Note: assumes camera_du == light_pdu
			if (!enable_vfc && c.is_movable()) {movable_cids.push_back(*i); continue;}

			if (c.dgroup_id >= 0) {
				vector<unsigned> const &group_cids(cdraw_groups.get_draw_group(c.dgroup_id, c));
				
				for (auto j = group_cids.begin(); j != group_cids.end(); ++j) {
					coll_obj const &c2(cdraw_groups.get_cobj(*j));
					c2.get_shadow_triangle_verts(verts, get_ndiv(c2, smap_sz, fixed_ndiv)); // no cube optimization here
				}
			}
			else if (c.type == COLL_CUBE) {z_sorted.push_back(make_pair(-c.d[2][1], *i));}
			else {c.get_shadow_triangle_verts(verts, get_ndiv(c, smap_sz, fixed_ndiv));}
		}
		end_block1(verts.size());
		sort(z_sorted.begin(), z_sorted.end());

		for (vector<pair<float, unsigned> >::const_iterator i = z_sorted.begin(); i != z_sorted.end(); ++i) {
			coll_objects[i->second].get_shadow_triangle_verts(verts, 1); // no view_dir
		}
		upload(verts);
	}

	void draw_shadow_sphere(point const &pos, float radius, int shader_loc, unsigned smap_sz, unsigned fixed_ndiv) {
		// it might be faster to draw circles, but probably not because VBOs can't be used,
		// since circles aren't rotationally symmetric and would need to be rotated to align with the projection vector
		shader_t::set_uniform_vector4d(shader_loc, vector4d(pos, radius));
		draw_sphere_vbo_pre_bound((fixed_ndiv ? fixed_ndiv : get_smap_ndiv(radius, smap_sz)), 0);
	}

	void draw_shadow_cobj(coll_obj const &c, unsigned smap_sz, unsigned fixed_ndiv, int shader_loc, unsigned char eflags) { // handle spheres specially
		if (c.type == COLL_TORUS) { // special case optimized for torus
			bind_vbo(0); // unbind sphere VBO
			unsigned const ndiv(get_smap_ndiv(c.radius, smap_sz));
			shader_t::set_uniform_vector4d(shader_loc, vector4d(all_zeros, 1.0));
			draw_rot_torus(c.points[0], c.norm, c.radius2, c.radius, ndiv, ndiv);
			bind_draw_sphere_vbo(0, 0); // no tex coords or normals
			return;
		}
		if (c.type == COLL_CAPSULE || c.type == COLL_SPHERE) {draw_shadow_sphere(c.points[0], c.radius, shader_loc, smap_sz, fixed_ndiv);}
		if (c.type == COLL_CAPSULE) {draw_shadow_sphere(c.points[1], c.radius2, shader_loc, smap_sz, fixed_ndiv);}
		if (c.type != COLL_SPHERE) {c.get_shadow_triangle_verts(dverts, get_ndiv(c, smap_sz, fixed_ndiv), 1, eflags);} // skip_spheres=1
	}

	void add_draw_dynamic(pos_dir_up const &pdu, unsigned smap_sz, unsigned fixed_ndiv, point const &camera_pos, vector3d const &light_dir, float back_face_thresh) {
		if (shadow_objs.empty() && movable_cids.empty()) return; // no dynamic objects
		//timer_t timer("Add Draw Dynamic");
		shader_t shader;
		shader.begin_shadow_map_shader(0, 1); // use_alpha_mask=0, enable_xlate_scale=1
		int const shader_loc(shader.get_uniform_loc("xlate_scale"));
		assert(shader_loc >= 0);
		bind_draw_sphere_vbo(0, 0); // no tex coords or normals
		bool const is_camera(dist_less_than(pdu.pos, camera_pos, 0.25*CAMERA_RADIUS));
		unsigned char eflags(0);
		UNROLL_3X(if (fabs(light_dir[i_]) > back_face_thresh) {eflags |= EFLAGS[i_][light_dir[i_] > 0.0];});

		for (auto i = movable_cids.begin(); i != movable_cids.end(); ++i) {
			coll_obj const &c(coll_objects.get_cobj(*i));
			if (c.no_shadow_map() || !c.is_movable()) continue; // should we remove it from the list in this case?
			if (!c.check_pdu_visible(pdu)) continue;
			if (dist_less_than(c.get_center_pt(), pdu.pos, pdu.near_) && c.contains_point(pdu.pos)) continue; // this cobj must be casting the light
			
			if (c.dgroup_id >= 0) {
				vector<unsigned> const &group_cids(cdraw_groups.get_draw_group(c.dgroup_id, c));
				
				for (auto j = group_cids.begin(); j != group_cids.end(); ++j) {
					draw_shadow_cobj(cdraw_groups.get_cobj(*j), smap_sz, fixed_ndiv, shader_loc, eflags);
				}
			}
			else {draw_shadow_cobj(c, smap_sz, fixed_ndiv, shader_loc, eflags);}
		}
		for (auto i = shadow_objs.begin(); i != shadow_objs.end(); ++i) {
			if (!pdu.sphere_visible_test(i->pos, i->radius)) continue; // VFC against light volume (may be culled earlier)
			if (is_camera && flashlight_on && i->is_player)  continue; // skip the camera shadow for flashlight
			if (i->pos == pdu.pos) continue; // this sphere must be casting the light
			//if (i->contains_point(pdu.pos)) continue; too strong

			if (i->ctype != COLL_SPHERE) {
				coll_obj const &c(coll_objects.get_cobj(i->cid));
				c.get_shadow_triangle_verts(dverts, get_ndiv(c, smap_sz, fixed_ndiv));
			}
			else {
				draw_shadow_sphere(i->pos, i->radius, shader_loc, smap_sz, fixed_ndiv);
			}
		}
		bind_vbo(0); // clear any bound sphere VBOs
		shader.set_uniform_vector4d(shader_loc, vector4d(all_zeros, 1.0)); // reset to identity transform
		render_dynamic();
		shader.end_shader();
	}

	void register_movable_cobj(unsigned cid) {
		coll_obj const &c(coll_objects.get_cobj(cid));
		if (!c.no_shadow_map() && c.is_movable()) {movable_cids.push_back(cid);}
	}
}; // smap_vertex_cache_t

smap_vertex_cache_t smap_vertex_cache;


bool shadow_map_enabled() {return (shadow_map_sz > 0);}
void free_smap_vbo() {smap_vertex_cache.free();}
void register_movable_cobj_shadow(unsigned cid) {smap_vertex_cache.register_movable_cobj(cid);}


xform_matrix get_texture_matrix(xform_matrix const &camera_mv_matrix) {
	
	// This matrix transforms every coordinate {x,y,z} to {x,y,z}* 0.5 + 0.5 
	// Moving from unit cube [-1,1] to [0,1]  
	glm::mat4 const bias(	
		0.5, 0.0, 0.0, 0.0, 
		0.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.5, 0.0,
		0.5, 0.5, 0.5, 1.0);
	return (bias * fgGetPJM() * fgGetMVM() * glm::affineInverse((glm::mat4)camera_mv_matrix));
}


void smap_texture_array_t::reserve_num_layers(unsigned num) {
	if (num <= num_layers) return; // have enough layers
	free_gl_state();
	num_layers = num;
}

unsigned smap_texture_array_t::new_layer() {
	assert(num_layers_used <= num_layers);
	if (num_layers_used == num_layers) {reserve_num_layers(max(2U*num_layers, 1U));} // double in size
	assert(num_layers_used < num_layers);
	return num_layers_used++;
}

void smap_texture_array_t::free_gl_state() {free_texture(tid); gpu_mem = 0;}

void smap_data_state_t::free_gl_state() {
	if (is_arrayed()) {gen_id = 0;} else {free_texture(local_tid);}
	free_fbo(fbo_id);
}

void smap_data_state_t::bind_tex_array(smap_texture_array_t *tex_arr_) { // must be unbound and bound to non-null
	assert(tex_arr_); assert(!tex_arr);
	tex_arr  = tex_arr_;
	layer_id = tex_arr->new_layer();
}

float const cascade_plane_dscales[NUM_CSM_CASCADES] = {0.03, 0.1, 0.3, 1.0};

float get_cascade_plane_dscale(unsigned cascade_ix) {
	assert(cascade_ix < NUM_CSM_CASCADES);
	float dscale(cascade_plane_dscales[cascade_ix]);
	if (world_mode == WMODE_GROUND) {dscale *= min(1.0f, (X_SCENE_SIZE + Y_SCENE_SIZE)/camera_pdu.far_);}
	return dscale;
}
float get_cascade_plane_dscale_last_at_1(unsigned cascade_ix) {
	return ((cascade_ix == NUM_CSM_CASCADES-1) ? 1.0 : get_cascade_plane_dscale(cascade_ix));
}
void set_csm_near_far(pos_dir_up &frustum, unsigned cascade_ix) {
	if (cascade_ix > 0) {frustum.near_ = frustum.far_*get_cascade_plane_dscale_last_at_1(cascade_ix-1);}
	frustum.far_ *= get_cascade_plane_dscale_last_at_1(cascade_ix);
}

bool smap_data_t::set_smap_shader_for_light(shader_t &s, int light, xform_matrix const *const mvm) const {

	if (!shadow_map_enabled() || !is_light_enabled(light)) return 0;
	point lpos; // unused
	bool const light_is_valid(light_valid(0xFF, light, lpos));
	string sm_tex_str("sm_tex"), sm_scale_str("sm_scale");
	s.add_uniform_int  (append_ix(sm_tex_str,   light, 0), tu_id);
	s.add_uniform_float(append_ix(sm_scale_str, light, 0), (light_is_valid ? 1.0 : 0.0));

	if (is_csm) {
		// cascade_plane_distances only needs to be done only once per light, but we have at most two lights, so it should be okay
		assert(cascade_matrices.size() == NUM_CSM_CASCADES);

		for (unsigned i = 0; i < NUM_CSM_CASCADES; ++i) {
			string cpd_str("cascade_plane_distances"), smap_matrix_str("smap_matrix"), matrix_str(append_ix(smap_matrix_str, light, 0));
			s.add_uniform_float(append_ix(cpd_str, i, 1), camera_pdu.far_*get_cascade_plane_dscale(i));
			s.add_uniform_matrix_4x4(append_ix(matrix_str, i, 1), cascade_matrices[i].get_ptr(), 0);
		}
	}
	else {
		xform_matrix tm(texture_matrix);
		if (mvm) {tm *= (*mvm) * glm::affineInverse((glm::mat4)fgGetMVM());} // Note: works for translate, but not scale?
		string smap_matrix_str("smap_matrix");
		s.add_uniform_matrix_4x4(append_ix(smap_matrix_str, light, 0), tm.get_ptr(), 0);
	}
	bind_smap_texture(light_is_valid);
	return 1;
}

void setup_and_bind_shadow_matrix_ubo() {
	shadow_matrix_ubo.allocate_with_size(MAX_DLIGHT_SMAPS*MAT4x4_SIZE, 0); // static draw?
	shadow_matrix_ubo.bind_base(0); // U_smap_matrix_dl is at binding point 0
	shadow_matrix_ubo.pre_render();
}

bool local_smap_data_t::set_smap_shader_for_light(shader_t &s, bool &arr_tex_set) const { // point/spot lights
	if (!shadow_map_enabled()) return 0;
	assert(tu_id >= LOCAL_SMAP_START_TU_ID);
	assert(is_arrayed());
	assert(!is_csm); // not supported for point/spot lights

	if (!arr_tex_set) { // Note: assumes all lights use the same texture array
		arr_tex_set = bind_smap_texture(); // not setup until we bind a valid texture
		bool const tex_ret(s.add_uniform_int("smap_tex_arr_dl", tu_id));
		if (!tex_ret) {cerr << "Error: unable to set shader uniform 'smap_tex_arr_dl'." << endl;}
		assert(tex_ret); // Note: we can assert this returns true, though it makes shader debugging harder
	}
	assert(layer_id < MAX_DLIGHT_SMAPS);
	upload_ubo_sub_data(texture_matrix.get_ptr(), layer_id*MAT4x4_SIZE, MAT4x4_SIZE); // shadow_matrix_ubo should be bound
	return 1;
}

// default empty shadow map for use in shaders when the shadow map hasn't been generated yet; for tiled terrain and disabled sun/moon
unsigned get_empty_smap_tid(bool is_csm) {
	if (empty_smap_tid == 0) {
		set_shadow_tex_params(empty_smap_tid, 0, is_csm); // is_array=0
		char const zero_data[16] = {0}; // large enough for any type
		if (is_csm) {glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, SMAP_INTERNAL_FORMAT, 1, 1, NUM_CSM_CASCADES, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_DATATYPE, zero_data);} // 1x1 texel
		else        {glTexImage2D(GL_TEXTURE_2D,       0, SMAP_INTERNAL_FORMAT, 1, 1,                   0, GL_DEPTH_COMPONENT, SHADOW_MAP_DATATYPE, zero_data);} // 1x1 texel
	}
	return empty_smap_tid;
}

void bind_default_sun_moon_smap_textures() { // for tiled terrain and building exteriors; uses different TUs compared to tiled terrain disable_shadow_maps()
	// ensure shadow map TU_IDs are valid in case we try to use them without binding to a tile;
	// maybe this fixes the occasional shader undefined behavior warning we get with the debug callback?
	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid_and_enabled(l)) continue;
		bind_texture_tu(get_empty_smap_tid(), GLOBAL_SMAP_START_TU_ID+l); // bind empty shadow map
	}
}

bool smap_data_t::bind_smap_texture(bool light_valid) const {
	// Note: the is_allocated() check shouldn't be required, but can happen in tiled terrain mode when switching between combined_gu mode (safe and conservative)
	// due to some disagreement between the update pass and draw pass during reflection drawing
	bool const use_tid(light_valid && is_allocated()); // otherwise, we know that sm_scale will be 0.0 and we won't do the lookup
	bind_texture_tu((use_tid ? get_tid() : get_empty_smap_tid(is_csm)), tu_id);
	return use_tid;
}


void upload_shadow_data_to_shader(shader_t &s) {
	for (unsigned l = 0; l < smap_data.size(); ++l) {smap_data[l].set_smap_shader_for_light(s, l);} // {sun, moon}
}

void set_smap_shader_for_all_lights(shader_t &s, float z_bias) {
	s.add_uniform_float("z_bias", z_bias);
	upload_shadow_data_to_shader(s);
}

void pre_bind_smap_tus(shader_t &s) {
	if (is_light_enabled(0)) {s.add_uniform_int("sm_tex0", GLOBAL_SMAP_START_TU_ID+0);}
	if (is_light_enabled(1)) {s.add_uniform_int("sm_tex1", GLOBAL_SMAP_START_TU_ID+1);}
}


// should this be a pos_dir_up member function?
void set_smap_mvm_pjm(point const &eye, point const &center, vector3d const &up_dir, float angle, float aspect, float near_clip, float far_clip) {

	fgMatrixMode(FG_PROJECTION);
	fgLoadIdentity();
	fgPerspective(2.0*angle/TO_RADIANS, aspect, near_clip, far_clip);
	fgMatrixMode(FG_MODELVIEW);
	fgLoadIdentity();
	fgLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up_dir.x, up_dir.y, up_dir.z);
}

pos_dir_up get_pt_cube_frustum_pdu(point const &pos_, cube_t const &bounds) {

	assert(bounds.is_strictly_normalized());
	point const center(bounds.get_cube_center());
	point pos(pos_);
	
	if (world_mode == WMODE_INF_TERRAIN) { // in TT mode, the tile shadow maps must be more stable so that they don't jitter when the scene shifts
		// making light_dir constant causes perf problems due to detection of dir changes causing spurious updates
		//pos += center; // this works, but makes shadows misalign at tile boundaries
		pos *= 10.0; // move the light further away (closer to infinite) to reduce the problem at the cost of some precision loss
		// maybe a better solution is to use a real directional light + ortho projection here?
	}
	float const dist(p2p_dist(pos, center));
	vector3d const light_dir((center - pos)/dist); // almost equal to lpos (point light)
	vector3d up_dir(zero_vector);
	up_dir[get_min_dim(light_dir)] = 1.0;
	point corners[8];
	get_cube_corners(bounds.d, corners);
	float const radius(bounds.get_bsphere_radius());

	// tighter bounds / higher quality / slower
	vector3d dirs[2]; // x, y (up)
	orthogonalize_dir(up_dir, light_dir, dirs[1], 1);
	cross_product(light_dir, dirs[1], dirs[0]);
	float rx(0.0), ry(0.0);

	for (unsigned i = 0; i < 8; ++i) {
		vector3d const delta((corners[i] - pos).get_norm());
		max_eq(rx, fabs(dot_product(dirs[0], delta)));
		max_eq(ry, fabs(dot_product(dirs[1], delta)));
	}
	float const frustum_skew_val(1.0 + 0.5*bounds.dz()/dist);
	float const angle(atan2(frustum_skew_val*ry, 1.0f)), aspect(rx/ry);
	return pos_dir_up(pos, light_dir, up_dir, angle, max(NEAR_CLIP, dist-radius), max(NEAR_CLIP, dist)+radius, aspect); // force near_clip < far_clip
}


void draw_scene_bounds_and_light_frustum(point const &lpos) {

	shader_t s;
	s.begin_color_only_shader();

	if (enable_ground_csm) {
		static pos_dir_up capture_pdu;
		if (display_mode & 0x10) {capture_pdu = camera_pdu;} // capture key = '5'

		if (capture_pdu.valid) {
			ensure_outlined_polygons();
			colorRGBA const frustum_colors[NUM_CSM_CASCADES] = {RED, GREEN, BLUE, YELLOW};

			for (unsigned i = 0; i < NUM_CSM_CASCADES; ++i) {
				s.set_cur_color(frustum_colors[i]);
				pos_dir_up frustum(capture_pdu);
				set_csm_near_far(frustum, i);
				frustum.draw_frustum();
			}
			set_fill_mode();
		}
	}
	else {
		enable_blend();
		// draw scene bounds
		s.set_cur_color(colorRGBA(1.0, 1.0, 1.0, 0.25)); // white
		draw_simple_cube(get_scene_bounds(), 0);
		// draw light frustum
		s.set_cur_color(colorRGBA(1.0, 1.0, 0.0, 0.25)); // yellow
		get_pt_cube_frustum_pdu(lpos, get_scene_bounds()).draw_frustum();
		disable_blend();
	}
	s.end_shader();
}


void set_shadow_tex_params(unsigned &tid, bool is_array, bool is_csm, bool use_white_border) {

	bool const nearest(0); // nearest filter: sharper shadow edges, but needs more biasing
	setup_texture(tid, 0, 0, 0, 0, 0, nearest, 1.0, is_array);
	// This is to allow usage of textureProj function in the shader
	int const target(get_2d_texture_target(is_array));
	glTexParameteri(target, GL_TEXTURE_COMPARE_MODE, (is_csm ? GL_NONE : GL_COMPARE_R_TO_TEXTURE));
	glTexParameteri(target, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);

	if (use_white_border) {
		// clamps area outside the shadow volume to unshadowed rather than extending the value at the border;
		// for buildings, this removes the incorrect shadow that extends above the player's head when near a wall, but causes some light leaking above walls
		glTexParameteri (target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri (target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		glTexParameterfv(target, GL_TEXTURE_BORDER_COLOR, &WHITE.R);
	}
}


bool no_sparse_smap_update() {

	if (world_mode != WMODE_GROUND) return 0;
	bool const leaf_wind(num_trees > 0 && (display_mode & 0x0100) != 0 && (tree_mode & 1) != 0 && tree_deadness < 1.0 && vegetation > 0.0);
	if (leaf_wind || !shadow_objs.empty() || platforms.any_active()) return 1;
	//return !coll_objects.drawn_ids.empty();
	//return !coll_objects.dynamic_ids.empty();
	return 0;
}


bool smap_data_t::needs_update(point const &lpos) {

	bool const ret(lpos != last_lpos);
	last_lpos = lpos;
	return (ret || !is_allocated());
}

bool ground_mode_smap_data_t::needs_update(point const &lpos) {

	bool const has_dynamic(!is_allocated() || scene_smap_vbo_invalid || no_sparse_smap_update()); // Note: force two frames of updates the first time the smap is created by setting has_dynamic
	bool const ret(smap_data_t::needs_update(lpos) || has_dynamic || last_has_dynamic || voxel_shadows_updated); // Note: see view clipping in indexed_vntc_vect_t<T>::render()
	last_has_dynamic = has_dynamic;
	return ret;
}


void smap_texture_array_t::ensure_tid(unsigned xsize, unsigned ysize) {

	assert(xsize > 0 && ysize > 0 && num_layers > 0);
	if (tid) return; // already allocated
	++gen_id;
	set_shadow_tex_params(tid, 1, 0); // is_array=1, is_csm=0
	glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, SMAP_INTERNAL_FORMAT, xsize, ysize, num_layers, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_DATATYPE, NULL);
	check_gl_error(630);
	gpu_mem += xsize*ysize*num_layers*smap_bytes_per_pixel;
}


glm::mat4 perspective_from_frustum(pos_dir_up const &pdu) {
	return glm::perspective(glm::radians(pdu.angle), (float)pdu.A, pdu.near_, pdu.far_);
}
glm::mat4 setup_csm_matrix(pos_dir_up const &frustum, point const &lpos) {
	point pts[9];
	frustum.get_frustum_corners(pts);

	if (world_mode == WMODE_GROUND) { // clamp to scene bounds for a tighter fit
		cube_t scene_bounds(get_scene_bounds());
		scene_bounds.expand_by(DX_VAL + DY_VAL); // expand slightly so that shadow map covers the edges
		for (unsigned i = 0; i < 8; ++i) {scene_bounds.clamp_pt(pts[i]);}
	}
	point const center(get_center_arb(pts, 8));
	vector3d const light_dir(lpos.get_norm()); // pointing toward the light
	unsigned npts(8);
	
	if (world_mode == WMODE_GROUND) { // project center in the light direction to zmax plane to make sure we capture all geometry
		float const plane_z(max(ztop, czmax)), dz(plane_z - center.z);
		if (dz > 0.0 && light_dir.z > 0.0) {pts[npts++] = center + (dz/light_dir.z)*light_dir;}
	}
	auto const light_view(glm::lookAt(vec3_from_vector3d(center + 0.1*light_dir), vec3_from_vector3d(center), vec3_from_vector3d(plus_y)));
	for (unsigned i = 0; i < npts; ++i) {pts[i] = vector3d_from_vec3(light_view * glm::vec4(vec3_from_vector3d(pts[i]), 1.0));}
	cube_t bcube(pts, npts);
	float const z_mult = 10.0f; // tune this parameter according to the scene
	if (bcube.z1() < 0.0) {bcube.z1() *= z_mult;} else {bcube.z1() /= z_mult;}
	if (bcube.z2() < 0.0) {bcube.z2() /= z_mult;} else {bcube.z2() *= z_mult;}
	glm::mat4 const light_projection(glm::ortho(bcube.x1(), bcube.x2(), bcube.y1(), bcube.y2(), bcube.z1(), bcube.z2()));
	return light_projection * light_view;
}

// hack to send a vector of matrices from CSM setup to rendering without passing it down the entire call tree
smap_data_t const *active_smap_data = nullptr;

bool is_csm_active() {return (active_smap_data != nullptr && active_smap_data->is_csm);}

void shader_csm_render_setup(shader_t &s) {
	assert(active_smap_data != nullptr);
	active_smap_data->set_csm_matrices(s);
}
void smap_data_t::set_csm_matrices(shader_t &s) const { // send matrices to the geometry shader
	assert(is_csm);
	assert(cascade_matrices.size() == NUM_CSM_CASCADES);

	for (unsigned i = 0; i < NUM_CSM_CASCADES; ++i) {
		string matrix_str("light_space_matrices");
		s.add_uniform_matrix_4x4(append_ix(matrix_str, i, 1), cascade_matrices[i].get_ptr(), 0);
	}
}

// if bounds is passed in, calculate pdu from it; otherwise, assume the user has alreay caclulated pdu;
// note that we store the light pos and other light-specific data in the shadow map, which means we can't reuse it across lights, but we can reuse the texture
void smap_data_t::create_shadow_map_for_light(point const &lpos, cube_t const *const bounds, bool use_world_space, bool no_update, bool force_update) {

	// setup render state
	assert(smap_sz > 0);
	bool do_update(!no_update && (is_csm || force_update || needs_update(lpos))); // must be called first, because this may indirectly update bounds
	xform_matrix camera_mv_matrix; // starts as identity matrix
	if (!use_world_space) {camera_mv_matrix = fgGetMVM();} // cache the camera modelview matrix before we change it
	fgPushMatrix();
	fgMatrixMode(FG_PROJECTION);
	fgPushMatrix();
	fgMatrixMode(FG_MODELVIEW);
	if (bounds && (do_update || !pdu.valid)) {pdu = get_pt_cube_frustum_pdu(lpos, *bounds);} // else pdu should have been set by the caller

	if (is_csm) {
		cascade_matrices.resize(NUM_CSM_CASCADES);
		
		for (unsigned i = 0; i < NUM_CSM_CASCADES; ++i) {
			pos_dir_up frustum(camera_pdu);
			set_csm_near_far(frustum, i);
			cascade_matrices[i] = setup_csm_matrix(frustum, lpos);
		}
	}
	else { // treat as perspective projection, even for directional lights
		set_smap_mvm_pjm(pdu.pos, (pdu.pos + pdu.dir), pdu.upv, pdu.angle, pdu.A, pdu.near_, pdu.far_);
		texture_matrix = get_texture_matrix(camera_mv_matrix);
	}
	check_gl_error(201);

	if (do_update) {
		// setup textures and framebuffer
		if (!is_allocated()) {
			free_fbo(fbo_id); // free existing fbo so that it can be recreated and bound to the new texture

			if (is_arrayed()) {
				tex_arr->ensure_tid(smap_sz, smap_sz); // create texture array if needed; point local to array texture so that we know it's bound
				gen_id = tex_arr->gen_id; // tag with current generation
			}
			else if (is_csm) {
				set_shadow_tex_params(local_tid, 1, 1, 1); // is_array=1, is_csm=1, use_white_border=1
				glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, SMAP_INTERNAL_FORMAT, smap_sz, smap_sz, NUM_CSM_CASCADES, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_DATATYPE, nullptr);
				check_gl_error(631);
			}
			else { // non-arrayed
				set_shadow_tex_params(local_tid, 0, 0); // is_array=0, is_csm=0
				glTexImage2D(GL_TEXTURE_2D, 0, SMAP_INTERNAL_FORMAT, smap_sz, smap_sz, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_DATATYPE, nullptr);
			}
		}
		assert(is_allocated());
		// render from the light POV to a FBO, store depth values only
		enable_fbo(fbo_id, get_tid(), 1, 0, is_csm, get_layer()); // is_depth_fbo=1, multisample=0
		glViewport(0, 0, smap_sz, smap_sz);
		glClear(GL_DEPTH_BUFFER_BIT);
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color rendering, we only want to write to the Z-Buffer
		// save state and update variables for fast rendering with correct clipping
		unsigned const orig_enabled_lights(enabled_lights);
		pos_dir_up const camera_pdu_(camera_pdu);
		camera_pdu     = pdu;
		enabled_lights = 0; // disable lighting so that shaders that auto-detect enabled lights don't try to do lighting
		ensure_filled_polygons();
		// render the scene
		check_gl_error(202);
		active_smap_data = this;
		render_scene_shadow_pass(lpos);
		active_smap_data = nullptr;
		// restore state variables
		camera_pdu     = camera_pdu_;
		enabled_lights = orig_enabled_lights;
		disable_fbo();
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		set_standard_viewport();
		set_fill_mode();
	}
	restore_prev_mvm_pjm_state();
	// Now rendering from the camera POV, using the FBO to generate shadows
	check_gl_error(203);
}

unsigned get_smap_bytes_per_pixel() {return smap_bytes_per_pixel;}

unsigned smap_data_t::get_gpu_mem() const {
	if (!is_allocated()) return 0;
	if (tex_arr)         return 0; // not tex_arr->gpu_mem, as this is shared across shadow maps and added once in local_smap_manager_t::get_gpu_mem()
	return smap_bytes_per_pixel*smap_sz*smap_sz*(is_csm ? NUM_CSM_CASCADES : 1);
}

void draw_mesh_shadow_pass(point const &lpos, unsigned smap_sz) {

	if (!(display_mode & 0x01) || ground_effects_level == 0) return;
	fgPushMatrix();
	float const val(1.0/dot_product(lpos.get_norm(), plus_z));
	fgTranslate(0.0, 0.0, -val*approx_pixel_width(smap_sz)); // translate down slightly to reduce shadow aliasing problems
	display_mesh(1);
	fgPopMatrix();
}

void draw_outdoor_shadow_pass(point const &lpos, unsigned smap_sz) {

	render_voxel_data(1);
	if (snow_shadows) {draw_snow(1);} // slow
	draw_trees(1); // shadow_only=1
	draw_scenery(1); // shadow_only=1
	draw_mesh_shadow_pass(lpos, smap_sz);
}

void ground_mode_smap_data_t::render_scene_shadow_pass(point const &lpos) {

	point const camera_pos_(camera_pos);
	camera_pos = lpos;
	vector3d const light_dir(-pdu.pos.get_norm()); // approximate as directional light; should be close enough for culling cube faces
	// add static objects
	smap_vertex_cache.add_cobjs(smap_sz, 0, 0); // no VFC for static cobjs
	smap_vertex_cache.render();
	render_models(1, 0);
	smap_vertex_cache.add_draw_dynamic(pdu, smap_sz, 0, camera_pos_, light_dir, 0.01);
	draw_outdoor_shadow_pass(lpos, smap_sz); // add snow, trees, scenery, and mesh
	voxel_shadows_updated = 0;
	camera_pos = camera_pos_;
}


// Note: not meant to shadow voxel terrain, snow, trees, scenery, mesh, etc. - basically designed to shadow cobjs and dynamic objects
void local_smap_data_t::render_scene_shadow_pass(point const &lpos) {
	
	point const camera_pos_(camera_pos);
	camera_pos = lpos;
	if (enable_depth_clamp) {glDisable(GL_DEPTH_CLAMP);} // no depth clamping (due to light fixtures in front of the near clip plane)

	if (world_mode == WMODE_GROUND) {
		// Note: using the global cobjs here may be less efficient for smap generation since we can't VFC (due to caching/sharing with other shadow maps)
		// however, it's simpler and more efficient for memory usage since there is only one buffer shared across all smaps, plus dlight smaps aren't generated each frame
		// Note: don't use fixed_ndiv in this call: since this may be shared with the global smap, both control flow paths should generate the same cobjs geometry;
		// it likely doesn't matter, since cobj geometry will generally be cached during the global smap pass (which is run first)
		smap_vertex_cache.add_cobjs(smap_sz, 0, 0); // no VFC for static cobjs
		smap_vertex_cache.render();
		render_models(2, 0);
		unsigned const fixed_ndiv = 24;
		// high back_face_thresh of 0.75 to avoid shadow artifacts for close cubes (should be at least 1/sqrt(2) for 90 deg FOV, and < 1.0)
		smap_vertex_cache.add_draw_dynamic(pdu, smap_sz, fixed_ndiv, camera_pos_, pdu.dir, 0.75);
		if (outdoor_shadows) {draw_outdoor_shadow_pass(lpos, smap_sz);}
	}
	else if (world_mode == WMODE_INF_TERRAIN) { // Note: not really a clean case split; should pass this in somehow, or use a different class in tiled terrain mode (cities)
		if (interior_shadow_maps) {
			draw_buildings(2, 0, zero_vector); // only need to draw buildings
		}
		else {
			render_models(2, 0, 1); // opaque only
		
			if (!interior_shadow_maps) { // all of this is here to draw tree shadows in tiled terrain mode, which is not needed for building interiors
				vector3d const xlate(get_tiled_terrain_model_xlate());
				camera_pdu.pos += xlate;
				fgPushMatrix();
				translate_to(-xlate);
				draw_tiled_terrain_decid_tree_shadows();
				fgPopMatrix();
				camera_pdu.pos -= xlate;
			}
		}
	}
	else {assert(0);} // not supported in universe mode
	if (enable_depth_clamp) {glEnable(GL_DEPTH_CLAMP);}
	camera_pos = camera_pos_;
}

bool local_smap_data_t::needs_update(point const &lpos) {
	
	// Note: scene_smap_vbo_invalid is reset at the end of the global create_shadow_map() call, so this call must be done before that
	bool has_dynamic(!get_tid() || scene_smap_vbo_invalid); // Note: scene_smap_vbo_invalid test is conservative
	
	for (auto i = shadow_objs.begin(); i != shadow_objs.end() && !has_dynamic; ++i) { // test dynamic objects
		has_dynamic |= pdu.sphere_visible_test(i->pos, i->radius);
	}
	if (scene_smap_vbo_invalid) { // maybe invalid due to moving cobjs (normally can't get here due to conservative test above)
		for (auto i = moving_cobjs.begin(); i != moving_cobjs.end() && !has_dynamic; ++i) {
			has_dynamic |= pdu.cube_visible(coll_objects.get_cobj(*i));
		}
	}
	if (!has_dynamic) {has_dynamic |= platforms.any_moving_platforms_in_view(pdu);} // test platforms
	// Note: maybe should check for moving objects as well - but they can only move if pushed by a dynamic shadow object (player) which is probably also in the light's view
	bool const ret(smap_data_t::needs_update(lpos) || has_dynamic || last_has_dynamic);
	last_has_dynamic = has_dynamic;
	return ret;
}


void create_shadow_map_inner(bool no_update) {

	cube_t const bounds(get_scene_bounds());
	
	for (unsigned l = 0; l < smap_data.size(); ++l) { // {sun, moon}
		point lpos;
		if (light_valid_and_enabled(l, lpos)) {smap_data[l].create_shadow_map_for_light(lpos, &bounds, 0, no_update);}
	}
}

void create_shadow_map() {

	if (!shadow_map_enabled()) return; // disabled
	//RESET_TIME;

	// save state
	int const do_zoom_(do_zoom), animate2_(animate2), display_mode_(display_mode);
	orig_camera_pdu = camera_pdu;

	// set to shadow map state
	do_zoom  = 0;
	animate2 = 0; // disable any animations or generated effects
	display_mode &= ~0x08; // disable occlusion culling

	// check VBO
	if (scene_smap_vbo_invalid == 2) {free_smap_vbo();} // force rebuild of VBO

	// render shadow maps to textures
	check_gl_error(198);
	ensure_smap_data();
	check_gl_error(199);
	create_shadow_map_inner(0); // no_update=0
	scene_smap_vbo_invalid = 0; // needs to be after dlights update

	// restore old state
	check_gl_error(200);
	do_zoom      = do_zoom_;
	animate2     = animate2_;
	display_mode = display_mode_;
	//PRINT_TIME("Shadow Map Creation");
}

void update_shadow_matrices() {
	if (!shadow_map_enabled() || smap_data.empty()) return; // disabled
	assert(scene_smap_vbo_invalid != 2);
	create_shadow_map_inner(1); // no_update=1
}


void free_shadow_map_textures() {
	for (unsigned l = 0; l < smap_data.size(); ++l) {smap_data[l].free_gl_state();}
	free_smap_vbo();
	free_light_source_gl_state(); // free any shadow maps within light sources
}




