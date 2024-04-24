// 3D World - light source implementation
// by Frank Gennari
// 5/17/15
#include "mesh.h"
#include "lightmap.h"
#include "shaders.h"
#include "sinf.h"
#include "shadow_map.h"


float const MAX_SMAP_FOV = 0.8; // acos(1.0 - 0.8) = 78 degrees

cube_t smap_light_clip_cube;

extern bool dl_smap_enabled, skip_light_vis_test, enable_dlight_shadows;
extern int display_mode, camera_coll_id, max_tius;
extern unsigned shadow_map_sz;
extern float fticks, CAMERA_RADIUS;
extern vector<light_source> light_sources_a;
extern vector<light_source_trig> light_sources_d;
extern coll_obj_group coll_objects;


bool bind_point_t::is_valid() { // used with placed dlights

	if (disabled) return 0;
	if (!bound  ) return 1; // if no binding point, always valid
	if (!valid  ) return 0; // already determined to be invalid

	if (bind_cobj < 0) { // cobj not yet found
		if (!check_point_contained_tree(bind_pos, bind_cobj, dynamic_cobj)) {valid = 0; return 0;}
		bind_pos = get_updated_bind_pos(); // update bind_pos to point to center of cobj, if movable
		return 1;
	}
	coll_obj const &cobj(coll_objects.get_cobj(bind_cobj));
	if (cobj.status != COLL_STATIC || !cobj.contains_point(get_updated_bind_pos())) {valid = 0; return 0;} // check status and also containment, in case coll id was reused
	return 1;
}

point bind_point_t::get_updated_bind_pos() const {
	// if bound to a movable cobj, return the cobj center so that this object moves with the cobj
	// Note: movable cobjs can only be created during scene loading, so we don't need to worry about the cobj slot being reused after the original cobj is destroyed
	if (!bound || !valid || disabled || bind_cobj < 0) return bind_pos; // no change
	coll_obj const &cobj(coll_objects.get_cobj(bind_cobj));
	return (cobj.is_movable() ? cobj.get_center_pt() : bind_pos);
}


// *** LIGHT_SOURCE IMPLEMENTATION ***


// radius == 0.0 is really radius == infinity (no attenuation)
light_source::light_source(float sz, point const &p, point const &p2, colorRGBA const &c, bool id, vector3d const &d, float bw, float ri, bool icf, float nc, float fc) :
	dynamic(id), enabled(1), is_cube_face(icf), radius(sz), radius_inv((radius == 0.0) ? 0.0 : 1.0/radius),
	r_inner(ri), bwidth(bw), near_clip(nc), far_clip(fc), pos(p), pos2(p2), dir(d.get_norm()), color(c)
{
	assert(bw > 0.0 && bw <= 1.0);
	assert(r_inner <= radius);
	assert(!(is_directional() && is_line_light())); // can't be both
}


void light_source::add_color(colorRGBA const &c) {

	color = color*color.alpha + c*c.alpha;
	color.alpha = 1.0;
}


float light_source::get_intensity_at(point const &p, point &updated_lpos) const {

	if (radius == 0.0) return color[3]; // no falloff
	updated_lpos = (is_cube_light ? cube_t(pos, pos2).closest_pt(p) : pos);

	if (is_line_light()) {
		vector3d const L(pos2 - pos);
		updated_lpos += L*CLIP_TO_01(dot_product((p - pos), L)/L.mag_sq());
	}
	if (fabs(p.z - updated_lpos.z) > radius) return 0.0; // fast test
	float const dist_sq(p2p_dist_sq(updated_lpos, p));
	if (dist_sq > radius*radius) return 0.0;
	float const rscale((radius - sqrt(dist_sq))*radius_inv);
	return rscale*rscale*color[3]; // quadratic 1/r^2 attenuation
}


float light_source::get_dir_intensity(vector3d const &obj_dir) const {

	if (!is_directional()) return 1.0;
	float const dp(dot_product(obj_dir, dir));
	if (dp >= 0.0f && (bwidth + LT_DIR_FALLOFF) < 0.5f) return 0.0;
	float const dp_norm(0.5f*(-dp*InvSqrt(obj_dir.mag_sq()) + 1.0f)); // dp = -1.0 to 1.0, bw = 0.0 to 1.0
	return CLIP_TO_01(2.0f*(dp_norm + bwidth + LT_DIR_FALLOFF - 1.0f)*LT_DIR_FALLOFF_INV);
}


cube_t light_source::calc_bcube(bool add_pad, float sqrt_thresh, bool clip_to_scene_bcube, float falloff) const {

	assert(radius > 0.0);
	assert(sqrt_thresh < 1.0);
	cube_t bcube(pos, pos2);
	bcube.expand_by(radius*(1.0 - sqrt_thresh));

	if (is_very_directional()) {
		cube_t bcube2;
		calc_bounding_cylin(sqrt_thresh, clip_to_scene_bcube, falloff).calc_bcube(bcube2);
		if (add_pad) {bcube2.expand_by(vector3d(DX_VAL, DY_VAL, DZ_VAL));} // add one grid unit
		bcube.intersect_with_cube(bcube2);
	}
	if (!custom_bcube.is_all_zeros()) {
		//assert(bcube.contains_cube(custom_bcube)); // too strong?
		assert(bcube.intersects(custom_bcube));
		bcube.intersect_with_cube(custom_bcube);
	}
	return bcube;
}

void light_source::get_bounds(cube_t &bcube, int bnds[3][2], float sqrt_thresh, bool clip_to_scene_bcube, vector3d const &bounds_offset) const {

	if (radius == 0.0) { // global light source
		for (unsigned d = 0; d < 3; ++d) {
			bcube.d[d][0] = -SCENE_SIZE[d];
			bcube.d[d][1] =  SCENE_SIZE[d];
			bnds[d][0]    = 0;
			bnds[d][1]    = MESH_SIZE[d]-1;
		}
	}
	else { // local point/spot/line light source
		bcube = calc_bcube(1, sqrt_thresh, clip_to_scene_bcube); // padded

		for (unsigned d = 0; d < 3; ++d) {
			UNROLL_2X(bnds[d][i_] = max(0, min(MESH_SIZE[d]-1, get_dim_pos((bcube.d[d][i_] + bounds_offset[d]), d)));)
		}
	}
}

float light_source::calc_cylin_end_radius(float falloff) const {
	float const d(1.0f - 2.0f*(bwidth + ((falloff > 0.0) ? falloff : LT_DIR_FALLOFF))); // use default falloff if zero
	return radius*sqrt(1.0f/(d*d) - 1.0f);
}
cylinder_3dw light_source::calc_bounding_cylin(float sqrt_thresh, bool clip_to_scene_bcube, float falloff) const {

	float const rad(radius*(1.0 - sqrt_thresh));
	if (is_line_light()) {return cylinder_3dw(pos, pos2, rad, rad);}
	assert(is_very_directional()); // not for use with point lights or spotlights larger than a hemisphere
	point pos2(pos + dir*rad);
	float end_radius((1.0 - sqrt_thresh)*calc_cylin_end_radius(falloff));

	if (clip_to_scene_bcube) { // Note: not correct in general, but okay for bcube calculation for large light sources
		point pos1(pos);
		cube_t const scene_bounds(get_scene_bounds());
		float const len1(p2p_dist(pos, pos2)), max_len(scene_bounds.furthest_dist_to_pt(pos));

		if (max_len < len1) { // spotlight cone projects outside the scene - clip the length to produce a smaller footprint
			float t(max_len/len1);
			end_radius *= t; // scale end_radius by the clipped line length
			pos2 = pos1 + t*(pos2 - pos1); // clip the line
			cube_t cylin_bcube;
			cylinder_3dw(pos, pos2, 0.0, end_radius).calc_bcube(cylin_bcube); // compute bounding cube of clipped cylinder
			cylin_bcube.intersect_with_cube(scene_bounds); // clip to scene bounds
			float const new_max_len(cylin_bcube.furthest_dist_to_pt(pos)); // compute max dist to new clipped scene, which should be a smaller (less conservative) value
			t = new_max_len/max_len;
			end_radius *= t; // scale end_radius by the clipped line length again
			pos2 = pos1 + t*(pos2 - pos1); // clip the line again
			//cout << TXT(len1) << TXT(max_len) << TXT(new_max_len) << endl;
		}
	}
	return cylinder_3dw(pos, pos2, 0.0, end_radius);
}


bool light_source::is_visible() const {

	if (!enabled) return 0;
	if (radius == 0.0) return 1;
	bool const line_light(is_line_light());
	point const lpos(line_light ? 0.5*(pos + pos2) : pos);
	float const lradius(line_light ? (radius + 0.5*p2p_dist(pos, pos2)) : radius); // use capsule bounding sphere for line light
	if (!camera_pdu.sphere_visible_test(lpos, lradius)) return 0; // view frustum culling

	if (line_light || is_very_directional()) { // use tighter bounding cube
		cube_t light_bc(calc_bcube());
		// since dlights are only applied within the scene bounds, we can ignore any part of the bcube (light volume) that's outside of the scene
		light_bc.intersect_with_cube(get_scene_bounds()); // optional optimization
		if (!camera_pdu.cube_visible(light_bc)) return 0;
	}
	if (!line_light) {
		if (radius < 0.5) return 1; // don't do anything more expensive for small light sources
		if (sphere_cobj_occluded(get_camera_pos(), pos, max(0.5f*radius, r_inner))) return 0; // approximate occlusion culling, can miss lights but rarely happens
	}
	if (skip_light_vis_test) return 1;
	if (dynamic || radius < 0.65 || !(display_mode & 0x08)) return 1; // dynamic/platform lights (common case), small/medium lights, or occlusion culling disabled
	unsigned const num_rays = 100;
	static rand_gen_t rgen;
	static vector<vector3d> dirs;
	static map<pair<point, point>, point> ray_map;
	//shader_t shader;
	//shader.begin_color_only_shader(RED);
	//RESET_TIME;
	point const camera(get_camera_pos());
	int prev_cindex(-1);
	if (!check_coll_line_tree(pos, camera, prev_cindex, camera_coll_id, 0, 1, 1, 0, 1)) return 1; // light center is visible
	unsigned cur_dir(0);
	bool const directional(is_directional()), very_dir(is_very_directional());
	vector3d vortho[2];
	if (very_dir) {get_ortho_vectors(dir, vortho);}
	float const cylin_end_radius(very_dir ? calc_cylin_end_radius() : 0.0), radius_adj(radius*(1.0 - SQRT_CTHRESH));

	if (dirs.empty()) { // start with 26 uniformly distributed directions
		for (int x = -1; x <= 1; ++x) {
			for (int y = -1; y <= 1; ++y) {
				for (int z = -1; z <= 1; ++z) {
					if (x == 0 && y == 0 && z == 0) continue;
					dirs.push_back(vector3d(x, y, z).get_norm());
				}
			}
		}
	}
	for (unsigned n = 0; n < num_rays; ++n) { // for static scene lights we do ray queries
		vector3d ray_dir;
		point start_pos(pos);
		
		if (very_dir && n < num_rays/4) { // uniformly spaced around the cylinder perimeter
			float const theta(TWO_PI*float(n)/float(num_rays/4));
			ray_dir = radius*dir + cylin_end_radius*(SINF(theta)*vortho[0] + COSF(theta)*vortho[1]);
		}
		else if (directional) { // randomly spaced within cylinder volume
			while (1) {
				if (cur_dir >= dirs.size()) {dirs.push_back(rgen.signed_rand_vector_norm());}
				ray_dir = dirs[cur_dir++];
				if ((bwidth + LT_DIR_FALLOFF) < 0.5f && dot_product(ray_dir, dir) < 0.0f) {ray_dir = -ray_dir;} // backwards
				if (get_dir_intensity(-ray_dir) > 0.0f) break;
			}
		}
		else { // randomly spaced around the unit sphere
			if (cur_dir >= dirs.size()) {dirs.push_back(rgen.signed_rand_vector_norm());}
			ray_dir = dirs[cur_dir++];
			if (cur_dir > 26 && dir != zero_vector && dot_product(dir, ray_dir) < 0.0) {ray_dir = -ray_dir;} // invert direction
		}
		if (line_light) {start_pos += (float(n)/float(num_rays-1))*(pos2 - pos);} // fixed spacing along the length of the line
		point const end_pos(start_pos + radius_adj*ray_dir);
		pair<point, point> const key(start_pos, end_pos);
		auto it(ray_map.find(key));
		point cpos;
		int cindex(-1);
		
		if (it != ray_map.end()) {cpos = it->second;} // intersection point is cached
		else { // not found in cache, computer intersection point and add it
			vector3d cnorm; // unused
			if (check_coll_line_exact_tree(start_pos, end_pos, cpos, cnorm, cindex, camera_coll_id, 0, 1, 1, 0, 1)) {cpos -= SMALL_NUMBER*ray_dir;} // move away from coll pos
			else {cpos = end_pos;} // clamp to end_pos if no int
			if (cindex < 0 || coll_objects[cindex].truly_static()) {ray_map[key] = cpos;}
		}
		//draw_subdiv_sphere(cpos, 0.01, N_SPHERE_DIV/2, 0, 0);
		if (!camera_pdu.sphere_visible_test(cpos, 0.1*radius)) continue; // point not visible
		if ((prev_cindex < 0 || !coll_objects[prev_cindex].line_intersect(cpos, camera)) && // doesn't intersect the previous cobj
			!check_coll_line_tree(cpos, camera, cindex, camera_coll_id, 0, 1, 1, 0, 0)) return 1; // visible
		prev_cindex = cindex;
	}
	//shader.end_shader();
	//PRINT_TIME("Light Source Vis");
	return 0; // not visible
}


void light_source::combine_with(light_source const &l) { // Note: unused

	assert(radius > 0.0);
	float const w1(radius*radius*radius), w2(l.radius*l.radius*l.radius), wsum(w1 + w2), wa(w1/wsum), wb(w2/wsum);
	radius     = pow(wsum, (1.0f/3.0f));
	radius_inv = 1.0/radius;
	pos       *= wa;
	pos       += l.pos*wb; // weighted average
	blend_color(color, color, l.color, wa, 1);
}

bool light_source::try_merge_into(light_source &ls) const {

	if (ls.radius < radius) return 0; // shouldn't get here because of radius sort
	if (!dist_less_than(pos, ls.pos, 0.2*min(HALF_DXY, radius))) return 0;
	if (ls.bwidth != bwidth || ls.r_inner != r_inner || ls.dynamic != dynamic) return 0;
	if (is_directional() && dot_product(dir, ls.dir) < 0.95) return 0;
	if (is_line_light() || ls.is_line_light()) return 0; // don't merge line lights
	if (is_neg_light () != ls.is_neg_light ()) return 0; // don't merge neg lights (looks bad)
	colorRGBA lcolor(color);
	float const rr(radius/ls.radius);
	lcolor.alpha *= rr*rr; // scale by radius ratio squared
	ls.add_color(lcolor);
	return 1;
}

void light_source::pack_to_floatv(float *data) const {

	// store light_source as: {pos.xyz, radius} {color.rgba} {dir.xyz|pos2.xyz, bwidth} [smap_index]
	// Note: we don't really need to store the z-component of dir because we can calculate it from sqrt(1 - x*x - y*y),
	//       but doing this won't save us any texture data so it's not worth the trouble
	assert(data);
	UNROLL_3X(*(data++) = pos[i_];)
	*(data++) = radius;
	UNROLL_3X(*(data++) = 0.5*(1.0 + color[i_]);) // map [-1,1] => [0,1] for negative light support
	*(data++) = color[3];

	if (is_line_light()) {
		UNROLL_3X(*(data++) = pos2[i_];)
		*(data++) = 0.0; // pack bwidth as 0 to indicate a line light
	}
	else {
		UNROLL_3X(*(data++) = 0.5*(1.0 + dir[i_]);) // map [-1,1] to [0,1]
		*(data++) = bwidth; // [0,1]
	}
	if (smap_enabled()) {
		assert(dl_smap_enabled);
		// the smap index is stored as a float in [0,1] and converted to an 8-bit int by multiplying by 255;
		// the int contains 7 index bits for up to 127 shadow maps + the 8th bit stores is_cube_face
		*(data++) = float(smap_index + (is_cube_face ? 128 : 0))/255.0f;
	}
	else {*(data++) = 0.0f;} // shadow map disabled
}

void light_source_trig::advance_timestep() {

	if (enabled && rot_rate != 0.0) {rotate_vector3d(rot_axis, rot_rate*fticks/TICKS_PER_SECOND, dir); dir.normalize();}

	if (bind_point_t::valid && bound) { // shift light if bound to a movable cobj
		point const new_bind_pos(get_updated_bind_pos());
		if (new_bind_pos != bind_pos) {shift_by(new_bind_pos - bind_pos); dynamic = 1;} // if the light moves, flag it as dynamic
	}
	if (!bind_point_t::valid || bind_point_t::disabled) {release_smap();} // free shadow map if invalid as an optimization
	if (!triggers.is_active() && !sensor.enabled()) return; // triggers and sensors not active
	if (active_time == 0.0 && triggers.check_for_activate_this_frame()) {register_activate(0);} // not active - check if activated by a non-player object
	enabled = (active_time > 0.0); // light on by default
	// if both triggers and sensor are enabled, use sensor to determine enabled state (but still check auto on/off time for triggers in case sensor is deactivated later)
	if (/*triggers.empty() &&*/ sensor.enabled()) {enabled = sensor.check_active();}
	
	if (enabled) {
		if (triggers.get_auto_off_time() > 0.0) {active_time = max(0.0f, (active_time - fticks));} // decrease active time in auto off mode
	}
	else {
		if (triggers.get_auto_on_time()  > 0.0) {inactive_time += fticks;} // increase inactive time in auto on mode
	}
}

void light_source_trig::shift_by(vector3d const &vd) {
	
	light_source::shift_by(vd);
	bind_point_t::shift_by(vd);
	triggers.shift_by(vd);
}

bool light_source_trig::need_update_indir() {

	if (!dynamic_indir) return 0;
	if (pos == last_pos && dir == last_dir) return 0;
	last_pos = pos;
	last_dir = dir;
	return 1;
}

void light_source_trig::set_rotate(vector3d const &axis, float rotate) {
	
	assert(axis != zero_vector);
	assert(dir != zero_vector);
	rot_rate = rotate;
	rot_axis = axis.get_norm();
}

bool light_source_trig::check_activate(point const &p, float radius, int activator) {

	//if (active_time > 0.0) return 1; // already activated, don't reset timing
	float const auto_on_time(triggers.get_auto_on_time());
	unsigned trigger_mode(0); // 0 = not triggered, 1 bit = proximity, 2 bit = action, 4 bit = auto on
	if (auto_on_time > 0.0 && inactive_time > TICKS_PER_SECOND*auto_on_time) {inactive_time = 0.0; trigger_mode = 4;} // turn on, reset inactive_time
	trigger_mode |= triggers.register_activator_pos(p, radius, activator, 1);
	if (trigger_mode == 0) return 0; // not yet triggered
	register_activate((trigger_mode & 2) != 0);
	return 1;
}

void light_source_trig::register_activate(bool player_triggered) {

	float const auto_off_time(triggers.get_auto_off_time());
	bool const is_off(active_time == 0.0);
	//if (auto_off_time == 0.0 || trigger.requires_action) {active_time = (is_off ? ((auto_off_time == 0.0) ? 1.0 : auto_off_time) : 0.0);} // toggle mode
	if (auto_off_time == 0.0)  {active_time = (is_off ? 1.0 : 0.0);} // toggle mode
	else if (player_triggered) {active_time = (is_off ? auto_off_time : 0.0);} // toggle mode from user action with auto off
	else {active_time = auto_off_time;} // reset active time (on duration)
	active_time *= TICKS_PER_SECOND; // convert from seconds to ticks
}


// ************ SHADOW MAPS ***********

// local shadow maps are for point/spot lights used for cases such as buildings and streetlights and are capped at 64 max to agree with the shader;
// this cap can be increased, but some GPUs may not have enough uniform slots for more, and it would take more GPU memory;
// when the player rotates, different light sources may enter the view frustum as lights with cached shadows leave it, which can thrash the cache;
// it may make sense to reserve some extra slots to help with this, though only if the max_shadow_maps config option is set less than the max;
// however, this doesn't work because the shader wants to index texture layer and shadow matrix by the same index,
// and we may not have enough uniforms for caching
class local_smap_manager_t {

	bool use_tex_array;
	vector<local_smap_data_t> smap_data;
	vector<unsigned> free_list;
	smap_texture_array_t smap_tex_arr;

	local_smap_data_t &get_smap(unsigned index) {assert(index < smap_data.size()); return smap_data[index];} // does bounds checking

	void search_free_list_slot_for_id_and_make_last(unsigned id) {
		for (unsigned i = 0; i < free_list.size(); ++i) {
			if (get_smap(free_list[i]).user_smap_id != id) continue; // wrong ID
			swap(free_list[i], free_list.back());
			return; // done
		}
	}
public:
	local_smap_manager_t(bool use_tex_array_=1) : use_tex_array(use_tex_array_) {}

	unsigned new_smap(unsigned size, unsigned user_smap_id, bool &matched_smap_id) {
		unsigned index(0);
		
		if (free_list.empty()) { // allocate a new smap
			index = smap_data.size();
			if (index >= MAX_DLIGHT_SMAPS) return 0; // not enough shader uniforms
			unsigned tu_id(LOCAL_SMAP_START_TU_ID);
			if (!use_tex_array) {tu_id += index;} // if not using texture arrays, we need to allocate a unique tu_id for each shadow map
			if ((int)tu_id >= max_tius) return 0; // not enough TIU's (none for texture array) - fail
			local_smap_data_t smd(tu_id);
			if (use_tex_array) {smd.bind_tex_array(&smap_tex_arr); assert(*smd.get_layer() == index);} // Note: must be <= GL_MAX_ARRAY_TEXTURE_LAYERS (which is 2048)
			smap_data.push_back(smd);
		}
		else { // use free list element
			if (user_smap_id > 0 || // caller requested a specific shadow map - see if it's available
				get_smap(free_list.back()).user_smap_id > 0) // caller requested an unbound smap, and the last free list element is used - try to find an unused smap
			{
				search_free_list_slot_for_id_and_make_last(user_smap_id); // Note: only one ID should match
			}
			index = free_list.back(); // use the most recently used free list entry, which may have been moved to the back by the loop above
			free_list.pop_back();
		}
		local_smap_data_t &smd(get_smap(index));
		assert(!smd.used);
		matched_smap_id      = (user_smap_id > 0 && smd.user_smap_id == user_smap_id);
		smd.used             = 1; // mark as used (for error checking)
		smd.last_has_dynamic = 1; // force recreation
		smd.outdoor_shadows  = 0; // reset to default
		smd.user_smap_id     = user_smap_id; // tag this shadow map with the caller's ID (if provided); if nonzero, this should be a unique value
		smd.last_lpos        = zero_vector;

		if (size > 0 && smd.smap_sz != size) { // size change - free and reallocate
			smd.free_gl_state();
			smd.smap_sz = size;
		}
		if (use_tex_array) {assert(*smd.get_layer() == index);}
		return index + 1; // offset by 1
	}
	void free_smap(unsigned index) {
		assert(index > 0 && index <= smap_data.size());
		assert(smap_data[index-1].used);
		smap_data[index-1].used = 0;
		free_list.push_back(index-1);
	}
	local_smap_data_t &get(unsigned index) { // Note: index is offset by 1; index 0 is invalid
		assert(index > 0 && index <= smap_data.size());
		assert(smap_data[index-1].used);
		return smap_data[index-1];
	}
	void invalidate_cached_smap_id(unsigned smap_id) {
		if (smap_id == 0) return; // invalid/unset smap_id
		
		for (auto &i : smap_data) { // there should be at most one match, but it shouldn't hurt to iterate until the end just in case
			if (i.user_smap_id == smap_id) {i.user_smap_id = 0;} // invalidate if the ID matches
		}
	}
	void free_gl_state() {
		for (auto &i : smap_data) {i.free_gl_state();}
		smap_tex_arr.free_gl_state();
	}
	void clear() {
		free_gl_state();
		smap_data.clear();
		free_list.clear();
	}
	unsigned get_gpu_mem() const {
		unsigned mem(0);
		for (auto &i : smap_data) {mem += i.get_gpu_mem();}
		return mem + smap_tex_arr.gpu_mem;
	}
};

unsigned const NUM_SMAP_MGRS = 2;
local_smap_manager_t local_smap_manager[NUM_SMAP_MGRS]; // {normal/city, building_interiors}

void free_light_source_gl_state() { // free shadow maps
	for (unsigned i = 0; i < NUM_SMAP_MGRS; ++i) {local_smap_manager[i].free_gl_state();}
}
unsigned get_dlights_smap_gpu_mem() {
	unsigned mem(0);
	for (unsigned i = 0; i < NUM_SMAP_MGRS; ++i) {mem += local_smap_manager[i].get_gpu_mem();}
	return mem;
}

local_smap_manager_t &light_source::get_smap_mgr() const {
	assert(smap_mgr_id < NUM_SMAP_MGRS);
	return local_smap_manager[smap_mgr_id];
}

void light_source::setup_and_bind_smap_texture(shader_t &s, bool &arr_tex_set) const {
	if (smap_index > 0) {get_smap_mgr().get(smap_index).set_smap_shader_for_light(s, arr_tex_set);}
}

pos_dir_up light_source::calc_pdu(bool dynamic_cobj, bool is_cube_face, float falloff) const {

	float const cos_theta(1.0 - min((is_cube_face ? 0.35f : MAX_SMAP_FOV), 2.0f*(bwidth + falloff))); // Note: cube face is 49.5 degrees (must be > 45)
	float const angle(acosf(cos_theta)); // half FOV
	int const dim(get_min_dim(dir));
	vector3d temp(zero_vector), up_dir;
	temp[dim] = 1.0; // choose up axis
	orthogonalize_dir(temp, dir, up_dir, 1);
	int cindex(-1);
	float t(0.0);
	vector3d cnorm; // unused
	float nclip(0.001*radius); // min value

	if (near_clip > 0.0) {nclip = max(nclip, near_clip);}
	else if (world_mode == WMODE_GROUND && check_point_contained_tree(pos, cindex, dynamic_cobj)) {
		// if light is inside a light fixture, move the near clip plane so that the light fixture cobj is outside the view frustum
		assert(cindex >= 0);
		point const start_pos(pos + dir*radius);
		if (coll_objects[cindex].line_int_exact(start_pos, pos, t, cnorm)) {nclip += (1.0 - t)*radius;}
	}
	float const fclip((far_clip == 0.0) ? radius : far_clip);
	return pos_dir_up(pos, dir, up_dir, angle, nclip, max(fclip, nclip+0.01f*radius), 1.0, 1); // force near_clip < far_clip; AR=1.0
}

void light_source::draw_light_cone(shader_t &shader, float alpha) const {

	if (!is_enabled_spotlight()) return; // disabled or not a spotlight
	//if (!is_visible()) return; // too slow?
	if (!camera_pdu.sphere_visible_test(pos, radius)) return; // view frustum culling
	if (dist_less_than(pos, camera_pos, 0.25*CAMERA_RADIUS)) return; // skip player flashlight
	unsigned const ndiv    = 40;
	unsigned const nstacks = 16;
	cylinder_3dw const cylin(calc_bounding_cylin());
	vector3d const dir(cylin.p2 - cylin.p1);
	vector<vert_norm_color> verts(2*(ndiv+1)); // triangle strip
	float const stack_delta(1.0/nstacks);

	for (unsigned n = 0; n < nstacks; ++n) { // stacks
		float const v1(n*stack_delta), v2(v1 + stack_delta);
		point const ce[2] = {(cylin.p1 + v2*dir), (cylin.p1 + v1*dir)}; // swap points to make second point the tip/zero radius
		vector3d v12;
		vector_point_norm const &vpn(gen_cylinder_data(ce, v2*cylin.r2, v1*cylin.r2, ndiv, v12));
		color_wrapper cw_tip(colorRGBA(color, (1.0 - v1)*alpha)), cw_end(colorRGBA(color, (1.0 - v2)*alpha));

		for (unsigned S = 0; S <= ndiv; ++S) {
			unsigned const s(S%ndiv), vix(2*S);
			vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]);
			verts[vix+0].assign(vpn.p[(s<<1)+0], normal, cw_end);
			verts[vix+1].assign(vpn.p[(s<<1)+1], normal, cw_tip);
		}
		draw_verts(verts, GL_TRIANGLE_STRIP); // untextured triangle strip with normals and color gradient from tip to end
	}
}


bool light_source_trig::is_shadow_map_enabled() const {
	if (!use_smap || no_shadows || shadow_map_sz == 0 || !enable_dlight_shadows) return 0;
	if (is_line_light())    return 0; // line lights don't support shadow maps
	if (dir == zero_vector) return 0; // point light: need cube map, skip for now
	//if (!is_enabled())      return 0; // disabled or destroyed
	if (is_directional() && !is_cube_face) {} // directional vs. hemisphere: use 2D shadow map for both
	return 1;
}

bool light_source_trig::check_shadow_map() {
	if (!is_shadow_map_enabled()) return 0;
	if (!is_enabled())            return 0; // disabled or destroyed
	bool const force_update(rot_rate != 0.0); // force shadow map update if rotating
	return setup_shadow_map(LT_DIR_FALLOFF, dynamic_cobj, outdoor_shadows, force_update, sm_size);
}

bool light_source::alloc_shadow_map(bool &matched_smap_id, unsigned sm_size) {
	if (smap_index == 0) {
		smap_index = get_smap_mgr().new_smap(sm_size, user_smap_id, matched_smap_id);
		if (smap_index == 0) return 0; // allocation failed (at max)
	}
	return 1;
}
void light_source::update_shadow_map(bool matched_smap_id, float falloff, bool dynamic_cobj, bool outdoor_shadows, bool force_update) {
	local_smap_data_t &smap(get_smap_mgr().get(smap_index));
	smap.pdu = calc_pdu(dynamic_cobj, is_cube_face, falloff); // Note: could cache this in the light source for static lights
	smap.outdoor_shadows = outdoor_shadows;

	if (0 && (display_mode & 0x10)) { // draw light/shadow frustum for debugging
		shader_t shader;
		shader.begin_color_only_shader(RED);
		smap.pdu.draw_frustum();
		shader.end_shader();
	}
	smap_light_clip_cube = custom_bcube;
	// if matched_smap_id==1, we can skip the shadow map update
	smap.create_shadow_map_for_light(pos, nullptr, 1, matched_smap_id, force_update); // no bcube, in world space, no texture array (layer=nullptr)
	smap_light_clip_cube.set_to_zeros();
}
bool light_source::setup_shadow_map(float falloff, bool dynamic_cobj, bool outdoor_shadows, bool force_update, unsigned sm_size) {
	bool matched_smap_id(0);
	if (!alloc_shadow_map(matched_smap_id, sm_size)) return 0;
	update_shadow_map(matched_smap_id, falloff, dynamic_cobj, outdoor_shadows, force_update);
	return 1;
}

void light_source::release_smap() {
	if (smap_index > 0) {get_smap_mgr().free_smap(smap_index); smap_index = 0;}
}

void light_source::invalidate_cached_smap_id(unsigned smap_id) const {
	get_smap_mgr().invalidate_cached_smap_id(smap_id);
}


template<typename T> void shift_ls_vect(T &v, vector3d const &vd) {
	for (auto i = v.begin(); i != v.end(); ++i) {i->shift_by(vd);}
}
void shift_light_sources(vector3d const &vd) {
	shift_ls_vect(light_sources_a, vd);
	shift_ls_vect(light_sources_d, vd);
}

void draw_spotlight_cones() {

	bool has_spotlight(0);
	for (auto i = light_sources_d.begin(); i != light_sources_d.end() && !has_spotlight; ++i) {has_spotlight |= i->is_enabled_spotlight();}
	if (!has_spotlight) return;
	float const alpha = 0.5;
	unsigned depth_tid(0);
	shader_t s;
	s.set_vert_shader("vert_plus_normal");
	s.set_frag_shader("depth_utils.part+spotlight_volume");
	s.set_prefix("#define USE_DEPTH_TRANSPARENCY", 1); // FS
	s.begin_shader();
	setup_depth_trans_texture(s, depth_tid);
	glDepthMask(GL_FALSE); // no depth writing
	enable_blend();
	for (auto i = light_sources_d.begin(); i != light_sources_d.end(); ++i) {i->draw_light_cone(s, alpha);}
	disable_blend();
	glDepthMask(GL_TRUE);
	s.end_shader();
	free_texture(depth_tid);
}

