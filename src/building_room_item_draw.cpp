// 3D World - Building Interior Room Item Drawing
// by Frank Gennari 4/17/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#include "profiler.h"
#include "openal_wrap.h"


bool const DEBUG_AI_COLLIDERS  = 0;
bool const ADD_WORKER_HARDHATS = 0; // doesn't look corret yet

unsigned room_geom_mem(0);
quad_batch_draw candle_qbd;
vect_room_object_t pending_objs;
object_model_loader_t building_obj_model_loader;

extern bool camera_in_building, player_in_tunnel, player_in_mall, building_alarm_active, is_cube_map_reflection, building_has_open_ext_door;
extern int display_mode, frame_counter, animate2, player_in_basement, player_in_elevator;
extern unsigned room_mirror_ref_tid;
extern float fticks, office_chair_rot_rate, building_ambient_scale;
extern point actual_player_pos, player_candle_pos, pre_reflect_camera_pos_bs;
extern vector4d clip_plane;
extern colorRGB cur_diffuse, cur_ambient;
extern cube_t smap_light_clip_cube, reflection_clip_cube;
extern pos_dir_up camera_pdu;
extern building_t const *player_building;
extern carried_item_t player_held_object;
extern room_object_t cur_room_mirror;
extern building_params_t global_building_params;

unsigned get_num_screenshot_tids();
tid_nm_pair_t get_phone_tex(room_object_t const &c);
template< typename T > void gen_quad_ixs(vector<T> &ixs, unsigned size, unsigned ix_offset);
void draw_emissive_billboards(quad_batch_draw &qbd, int tid);
void draw_car_in_pspace(car_t &car, shader_t &s, vector3d const &xlate, bool shadow_only, unsigned btype);
void set_car_model_color(car_t &car, unsigned btype);
bldg_obj_type_t get_taken_obj_type(room_object_t const &obj);
int get_toilet_paper_nm_id();
void setup_monitor_screen_draw(room_object_t const &monitor, rgeom_mat_t &mat, std::string &onscreen_text);
void add_tv_or_monitor_screen(room_object_t const &c, rgeom_mat_t &mat, std::string const &onscreen_text, rgeom_mat_t *text_mat);
void get_shower_head_pos_dir(room_object_t const &c, point &head_pos, vector3d &head_dir);
bool check_clock_time();
bool have_fish_model();
void register_fishtank(room_object_t const &obj, bool is_visible);
void end_fish_draw(shader_t &s, bool inc_pools_and_fb);
void calc_cur_ambient_diffuse();
void reset_interior_lighting_and_end_shader(shader_t &s);
bool has_cars_enabled();

bool has_key_3d_model      () {return building_obj_model_loader.is_model_valid(OBJ_MODEL_KEY);}
bool has_office_chair_model() {return building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR);}
bool is_flashing_light_on  () {return (tid_nm_pair_t(RED_TEX).get_emissive_val() > 0.5);}

colorRGBA room_object_t::get_model_color() const {return building_obj_model_loader.get_avg_color(get_model_id());}


class rgeom_alloc_t {
	deque<rgeom_storage_t> free_list; // one per unique texture ID/material
public:
	void alloc_safe(rgeom_storage_t &s) {
#pragma omp critical(rgeom_alloc)
		alloc(s);
	}
	void alloc(rgeom_storage_t &s) { // attempt to use free_list entry to reuse existing capacity
		if (free_list.empty()) return; // no pre-alloc
		//cout << TXT(free_list.size()) << TXT(free_list.back().get_tot_vert_capacity()) << endl; // total mem usage is 913/1045

		// try to find a free list element with the same tex so that we balance out material memory usage/capacity better
		for (unsigned i = 0; i < free_list.size(); ++i) {
			if (!free_list[i].tex.is_compatible(s.tex)) continue;
			s.swap_vectors(free_list[i]); // transfer existing capacity from free list
			free_list[i].swap(free_list.back());
			free_list.pop_back();
			return; // done
		}
	}
	void free(rgeom_storage_t &s) {
		s.clear(); // in case the caller didn't clear it
		if (s.get_mem_usage() == 0) return; // no memory allocated, no point in adding to the free list
		free_list.push_back(rgeom_storage_t(s.tex)); // record tex of incoming element
		s.swap_vectors(free_list.back()); // transfer existing capacity to free list; clear capacity from s
	}
	size_t get_mem_usage() const {
		size_t mem(free_list.size()*sizeof(rgeom_storage_t));
		for (auto i = free_list.begin(); i != free_list.end(); ++i) {
			//cout << i->tex.tid << "\t" << i->tex.shadowed << "\t" << i->get_tot_vert_capacity() << "\t" << i->get_mem_usage() << endl; // TESTING
			mem += i->get_mem_usage();
		}
		return mem;
	}
	size_t size() const {return free_list.size();}
	void clear_free_list() {free_list.clear();}
};
rgeom_alloc_t rgeom_alloc; // static allocator with free list, shared across all buildings; not thread safe

void print_building_rgeom_stats() {
	size_t const size(rgeom_alloc.size());
	if (size > 0) {cout << "rgeom_alloc: size: " << size << " mem (MB): " << in_mb(rgeom_alloc.get_mem_usage()) << endl;} // start=462MB
}
void clear_building_rgeom_free_list() {rgeom_alloc.clear_free_list();} // unused, but may be useful for testing


vbo_cache_t::vbo_cache_entry_t vbo_cache_t::alloc(unsigned size, bool is_index) {
	assert(size > 0); // not required, but a good sanity check
	auto &e(entries[is_index]);
	unsigned const max_size(size + size/5); // no more than 20% wasted cap
	auto best_fit(e.end());
	unsigned target_sz(size);

	for (auto i = e.begin(); i != e.end(); ++i) {
		target_sz += 8; // increase with every iteration
		if (i->size < size || i->size > max_size) continue; // too small or too large
		if (best_fit != e.end() && i->size >= best_fit->size) continue; // not the best fit
		best_fit = i;
		if (i->size < target_sz) break; // close fit, done
	}
	if (best_fit != e.end()) { // found a compatible VBO to reuse
		vbo_cache_entry_t const ret(*best_fit); // deep copy
		swap(*best_fit, e.back());
		e.pop_back();
		++v_reuse; s_reuse += ret.size; ++v_used; s_used += ret.size;
		assert(v_free > 0); assert(s_free >= size);
		--v_free; s_free -= size;
		return ret; // done
	} // for i
	++v_alloc; s_alloc += size; ++v_used; s_used += size;
	return vbo_cache_entry_t(create_vbo(), 0); // create a new VBO
}
void vbo_cache_t::free(unsigned &vbo, unsigned size, bool is_index) {
	if (!vbo) return; // nothing allocated
	if (size == 0) {delete_and_zero_vbo(vbo); return;} // shouldn't get here?
	assert(v_used > 0); assert(s_used >= size);
	--v_used; s_used -= size; ++v_free; s_free += size;
	entries[is_index].emplace_back(vbo, size);
	vbo = 0;
}
void vbo_cache_t::clear() { // unused
	for (unsigned d = 0; d < 2; ++d) {
		for (vbo_cache_entry_t &entry : entries[d]) {
			delete_vbo(entry.vbo);
			room_geom_mem -= entry.size;
		}
		entries[d].clear();
	}
}
void vbo_cache_t::print_stats() const {
	// v_alloc / s_alloc: number of VBOs / size allocated
	// v_used  / s_used : number of VBOs / size currently in use
	// v_reuse / s_reuse: number of VBOs / size reused, cumulative
	// v_free  / s_free : number of VBOs / size in free list
	cout << "VBOs: A " << v_alloc << " U " << v_used << " R " << v_reuse << " F " << v_free
		 << " | SZ MB: A " << (s_alloc>>20) << " U " << (s_used>>20) << " R " << (s_reuse>>20) << " F " << (s_free>>20) << endl; // in MB
}

/*static*/ vbo_cache_t rgeom_mat_t::vbo_cache;

/*static*/ void rgeom_mat_t::print_vbo_cache_stats() {
	if (vbo_cache.in_use()) {vbo_cache.print_stats();}
}
void print_rgeom_vbo_cache_stats() {rgeom_mat_t::print_vbo_cache_stats();}

void rgeom_storage_t::clear(bool free_memory) {
	if (free_memory) {clear_container(quad_verts);} else {quad_verts.clear();}
	if (free_memory) {clear_container(itri_verts);} else {itri_verts.clear();}
	if (free_memory) {clear_container(indices   );} else {indices   .clear();}
}
void rgeom_storage_t::swap_vectors(rgeom_storage_t &s) { // Note: doesn't swap tex
	quad_verts.swap(s.quad_verts);
	itri_verts.swap(s.itri_verts);
	indices.swap(s.indices);
}
void rgeom_storage_t::swap(rgeom_storage_t &s) {
	swap_vectors(s);
	std::swap(tex, s.tex);
}

void rgeom_mat_t::clear() {
	clear_vbos();
	clear_vectors();
	num_verts = num_ixs = 0;
}
void rgeom_mat_t::clear_vbos() {
	vbo_cache.free(vao_mgr.vbo,  vert_vbo_sz, 0);
	vbo_cache.free(vao_mgr.ivbo, ixs_vbo_sz,  1);
	vao_mgr.clear_vaos(); // Note: VAOs not reused because they generally won't be used with the same {vbo, ivbo} pair
	vert_vbo_sz = ixs_vbo_sz = 0;
}

void rotate_verts(vector<rgeom_mat_t::vertex_t> &verts, building_t const &building) {
	point const center(building.bcube.get_cube_center());

	for (auto &v : verts) {
		building.do_xy_rotate(center, v.v);
		vector3d n(v.get_norm());
		building.do_xy_rotate_normal(n);
		v.set_norm(n);
	}
}
void rgeom_mat_t::create_vbo(building_t const &building) {
	if (building.is_rotated()) { // rotate all vertices to match the building rotation
		rotate_verts(quad_verts, building);
		rotate_verts(itri_verts, building);
	}
	create_vbo_inner();
	rgeom_alloc.free(*this); // vertex and index data is no longer needed and can be returned to the free list
}
void rgeom_mat_t::create_vbo_inner() {
	assert(itri_verts.empty() == indices.empty());
	unsigned const qsz(quad_verts.size()*sizeof(vertex_t)), itsz(itri_verts.size()*sizeof(vertex_t)), tot_verts_sz(qsz + itsz);
	// in most cases when num_verts starts out nonzero there is no actual update for this material, but accurately skipping the VBO update is difficult;
	// hashing the vertex data is too slow, and simply summing the verts is inaccurate for things like light switch rotations and buildings far from the origin
	num_verts = quad_verts.size() + itri_verts.size();
	if (num_verts == 0) return; // nothing to do
	gen_quad_ixs(indices, 6*(quad_verts.size()/4), itri_verts.size()); // append indices for quad_verts
	num_ixs = indices.size();
	unsigned const ix_data_sz(num_ixs*sizeof(unsigned));

	if (vao_mgr.vbo && tot_verts_sz <= vert_vbo_sz && vao_mgr.ivbo && ix_data_sz <= ixs_vbo_sz) { // reuse previous VBOs
		update_indices(vao_mgr.ivbo, indices);
		check_bind_vbo(vao_mgr.vbo);
	}
	else { // create a new VBO
		clear_vbos(); // free any existing VBO memory
		auto vret(vbo_cache.alloc(tot_verts_sz, 0)); // verts
		auto iret(vbo_cache.alloc(ix_data_sz,   1)); // indices
		vao_mgr.vbo  = vret.vbo;
		vao_mgr.ivbo = iret.vbo;
		check_bind_vbo(vao_mgr.vbo);

		if (vret.size == 0) { // newly created
			vert_vbo_sz = tot_verts_sz;
			upload_vbo_data(nullptr, tot_verts_sz);
			room_geom_mem += tot_verts_sz;
		}
		else { // existing
			vert_vbo_sz = vret.size;
			assert(tot_verts_sz <= vert_vbo_sz);
		}
		if (iret.size == 0) { // newly created
			ixs_vbo_sz = ix_data_sz;
			upload_to_vbo(vao_mgr.ivbo, indices, 1, 1);
			room_geom_mem += indices.size()*sizeof(unsigned);
		}
		else { // existing
			ixs_vbo_sz = iret.size;
			assert(ix_data_sz <= ixs_vbo_sz );
			update_indices(vao_mgr.ivbo, indices, ix_data_sz);
		}
	}
	if (itsz > 0) {upload_vbo_sub_data(itri_verts.data(), 0,    itsz);}
	if (qsz  > 0) {upload_vbo_sub_data(quad_verts.data(), itsz, qsz );}
	bind_vbo(0);
	check_gl_error(475);

	if (num_verts >= 32) {dir_mask = 63;} // too many verts, assume all orients
	else {
		dir_mask = 0;
		for (unsigned n = 0; n < 2; ++n) {
			for (auto const &v : (n ? itri_verts : quad_verts)) {
				for (unsigned d = 0; d < 3; ++d) {
					if (v.n[d] < 0) {dir_mask |= 1<<(2*d);} else if (v.n[d] > 0) {dir_mask |= 1<<(2*d+1);}
				}
			}
		} // for n
	}
	// calculate bcube, for use in VFC for the shadow pass;
	// only enable for small blocks to reduce runtime overhead, plus they're more likely to be occluded
	if (num_verts >= 256) {bcube.set_to_zeros();}
	else {
		bcube.set_from_point(itri_verts.empty() ? quad_verts.front().v : itri_verts.front().v);
		for (auto const &v : itri_verts) {bcube.union_with_pt(v.v);}
		for (auto const &v : quad_verts) {bcube.union_with_pt(v.v);}
	}
}

bool brg_batch_draw_t::has_ext_geom() const {
	for (tile_block_t const &tb : ext_by_tile) {
		for (mat_entry_t const &e : tb.to_draw) {
			if (!e.mats.empty()) return 1;
		}
	}
	return 0;
}
void brg_batch_draw_t::clear() {
	to_draw.clear();
	ext_by_tile.clear();
	tid_to_first_mat_map.clear();
	cur_tile_slot = 0;
}
void brg_batch_draw_t::set_camera_dir_mask(point const &camera_bs, cube_t const &bcube) {
	camera_dir_mask = 0;
	for (unsigned d = 0; d < 3; ++d) {
		if (camera_bs[d] < bcube.d[d][1]) {camera_dir_mask |= 1<<(2*d  );}
		if (camera_bs[d] > bcube.d[d][0]) {camera_dir_mask |= 1<<(2*d+1);}
	}
}
void brg_batch_draw_t::next_tile(cube_t const &bcube) {
	for (unsigned i = 0; i < ext_by_tile.size(); ++i) { // try to find an unused slot
		if (!ext_by_tile[i].bcube.is_all_zeros()) continue; // already used
		ext_by_tile[i].bcube = bcube;
		cur_tile_slot = i;
		return;
	}
	cur_tile_slot = ext_by_tile.size();
	ext_by_tile.emplace_back(bcube); // create a new slot
}
void brg_batch_draw_t::add_material(rgeom_mat_t const &m, bool is_ext_tile) {
	if (is_ext_tile) {assert(cur_tile_slot < ext_by_tile.size());}
	vector<mat_entry_t>& dest(is_ext_tile ? ext_by_tile[cur_tile_slot].to_draw : to_draw);

	for (auto &i : dest) { // check all existing materials for a matching texture, etc.
		if (i.tex.is_compat_ignore_shadowed(m.tex)) {i.mats.push_back(&m); return;} // found existing material
	}
	dest.emplace_back(m); // add a new material entry
}
void brg_batch_draw_t::draw_and_clear_batch(vector<mat_entry_t> &batch, tid_nm_pair_dstate_t &state) {
	for (auto &i : batch) {
		if (i.mats.empty()) continue; // empty slot
		i.tex.set_gl(state);
		for (auto const &m : i.mats) {m->draw_inner(0);} // shadow_only=0
		i.tex.unset_gl(state);
		i.mats.clear(); // clear mats but not batch
	}
}
void brg_batch_draw_t::draw_and_clear(shader_t &s) {
	if (to_draw.empty()) return;
	tid_nm_pair_dstate_t state(s);
	enable_blend(); // needed for rugs, book, and sign text
	draw_and_clear_batch(to_draw, state);
	disable_blend();
	indexed_vao_manager_with_shadow_t::post_render();
}
void brg_batch_draw_t::draw_and_clear_ext_tiles(shader_t &s, vector3d const &xlate) {
	if (ext_by_tile.empty()) return;
	tid_nm_pair_dstate_t state(s);
	enable_blend(); // needed for sign text

	for (tile_block_t &tb : ext_by_tile) {
		if (tb.to_draw.empty()) continue; // skip empty batches
		try_bind_tile_smap_at_point((tb.bcube.get_cube_center() + xlate), s);
		draw_and_clear_batch(tb.to_draw, state);
	}
	disable_blend();
	indexed_vao_manager_with_shadow_t::post_render();
}
void brg_batch_draw_t::clear_ext_tiles() {
	for (tile_block_t &tb : ext_by_tile) {tb.bcube = cube_t();} // reset for next frame
}

// shadow_only: 0=non-shadow pass, 1=shadow pass, 2=shadow pass with alpha mask texture
void rgeom_mat_t::draw(tid_nm_pair_dstate_t &state, brg_batch_draw_t *bbd, int shadow_only, int reflection_pass, bool exterior_geom) {
	assert(!(exterior_geom && shadow_only)); // exterior geom shadows are not yet supported
	if (shadow_only && !en_shadows)         return; // shadows not enabled for this material (picture, whiteboard, rug, etc.)
	if (shadow_only && tex.emissive == 1.0) return; // assume this is a light source and shouldn't produce shadows
	//if (reflection_pass == 2 && !en_shadows && tex.tid < 0)  return; // disable small non-shadow casting untextured objects in cube map reflections (optimization)
	if (reflection_pass == 2 && tex.no_reflect) return; // disable for cube map reflections (optimization)
	if (reflection_pass && tex.tid == REFLECTION_TEXTURE_ID) return; // don't draw reflections of mirrors as this doesn't work correctly
	if (bbd != nullptr  && tex.tid == REFLECTION_TEXTURE_ID) return; // only draw mirror reflections for player building (which has a null bbd)
	if (num_verts == 0) return; // Note: should only happen when reusing materials and all objects using this material were removed
	// VFC test for shadow pass on sparse materials that have their bcubes calculated; only really helps with backrooms;
	// here we don't add xlate to bcube because it's the location of a light source that's already in building space, not camera space
	if (shadow_only && !bcube.is_all_zeros() && !camera_pdu.cube_visible(bcube)) return;
	vao_setup(shadow_only);

	// Note: the shadow pass doesn't normally bind textures and set uniforms, so we don't need to combine those calls into batches
	if (bbd != nullptr && !shadow_only) { // add to batch draw (optimization)
		if (dir_mask > 0 && bbd->camera_dir_mask > 0 && (dir_mask & bbd->camera_dir_mask) == 0) return; // check for visible surfaces
		bbd->add_material(*this, exterior_geom);
	}
	else { // draw this material now
		if (shadow_only != 1) {tex.set_gl  (state);} // ignores texture scale for now; enable alpha texture for shadow pass
		draw_inner(shadow_only);
		if (shadow_only != 1) {tex.unset_gl(state);}
	}
}
void rgeom_mat_t::pre_draw(int shadow_only) const {
	vao_mgr.pre_render(shadow_only != 0);
}
void rgeom_mat_t::draw_geom() const {
	glDrawRangeElements(GL_TRIANGLES, 0, num_verts, num_ixs, GL_UNSIGNED_INT, nullptr);
	++num_frame_draw_calls;
}
void rgeom_mat_t::draw_inner(int shadow_only) const {
	pre_draw(shadow_only);
	draw_geom();
}
void rgeom_mat_t::vao_setup(bool shadow_only) {
	vao_mgr.create_and_upload(vector<vertex_t>(), vector<unsigned>(), shadow_only, 0, 1); // pass empty vectors because data is already uploaded; dynamic_level=0, setup_pointers=1
}
void rgeom_mat_t::upload_draw_and_clear(tid_nm_pair_dstate_t &state) { // Note: called by draw_interactive_player_obj() and water_draw_t
	if (empty()) return; // nothing to do; can this happen?
	create_vbo_inner();
	draw(state, nullptr, 0, 0, 0); // no brg_batch_draw_t, shadow=reflection=exterior=0
	clear();
	indexed_vao_manager_with_shadow_t::post_render();
}

void building_materials_t::clear() {
	invalidate();
	for (iterator m = begin(); m != end(); ++m) {m->clear();}
	vector<rgeom_mat_t>::clear();
}
unsigned building_materials_t::count_all_verts(bool shadow_only, bool reflect_only) const {
	unsigned num_verts(0);
	
	for (const_iterator m = begin(); m != end(); ++m) {
		if (shadow_only  && !m->en_shadows    ) continue;
		if (reflect_only &&  m->tex.no_reflect) continue;
		num_verts += m->num_verts;
	}
	return num_verts;
}
rgeom_mat_t &building_materials_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows) {
	// for now we do a simple linear search because there shouldn't be too many unique materials
	for (iterator m = begin(); m != end(); ++m) {
		if (!m->tex.is_compatible(tex)) continue;
		if (inc_shadows) {m->enable_shadows();} // Note: m->en_shadows should already be set
		// tscale diffs don't make new materials; copy tscales from incoming tex; this field may be used locally by the caller, but isn't used for drawing
		m->tex.tscale_x = tex.tscale_x; m->tex.tscale_y = tex.tscale_y;
		if (m->get_tot_vert_capacity() == 0) {rgeom_alloc.alloc_safe(*m);} // existing but empty entry, allocate capacity from the allocator free list
		return *m;
	}
	emplace_back(tex); // not found, add a new material
	if (inc_shadows) {back().enable_shadows();}
	rgeom_alloc.alloc_safe(back());
	return back();
}
void building_materials_t::create_vbos(building_t const &building) { // up to ~100 materials and ~2M verts
	for (iterator m = begin(); m != end(); ++m) {m->create_vbo(building);}
	valid = 1;
}
void building_materials_t::draw(brg_batch_draw_t *bbd, shader_t &s, int shadow_only, int reflection_pass, bool exterior_geom) {
	if (!valid) return; // pending generation of data, don't draw yet
	//highres_timer_t timer("Draw Materials"); // 0.0168
	tid_nm_pair_dstate_t state(s);
	for (iterator m = begin(); m != end(); ++m) {m->draw(state, bbd, shadow_only, reflection_pass, exterior_geom);}
}
void building_materials_t::upload_draw_and_clear(shader_t &s) {
	tid_nm_pair_dstate_t state(s);
	for (iterator m = begin(); m != end(); ++m) {m->upload_draw_and_clear(state);}
}

void building_room_geom_t::add_tquad(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex,
	colorRGBA const &color, bool invert_tc_x, bool exclude_frame, bool no_tc)
{
	assert(tquad.npts == 4); // quads only, for doors
	add_tquad_to_verts(bg, tquad, bcube, tex, color, mats_doors.get_material(tex, 1).quad_verts, invert_tc_x, exclude_frame, no_tc, 1); // inc_shadows=1, no_rotate=1
}

void building_room_geom_t::clear() {
	clear_materials();
	objs.clear();
	light_bcubes.clear();
}
void building_room_geom_t::clear_materials() { // clears material VBOs
	mats_static .clear();
	mats_alpha  .clear();
	mats_dynamic.clear();
	mats_doors  .clear();
	mats_lights .clear();
	mats_exterior.clear();
	mats_ext_detail.clear();
	obj_model_insts.clear(); // these are associated with static VBOs
	door_handles   .clear();
	for (unsigned d = 0; d < 2; ++d) {mats_glass[d].clear();}
	clear_small_materials();
}
void building_room_geom_t::clear_small_materials() { // clears small material VBOs
	mats_small   .clear();
	mats_text    .clear();
	mats_amask   .clear();
	mats_detail  .clear();
	mats_alpha_sm.clear();
}
// Note: used for room lighting changes; detail object changes are not supported
void building_room_geom_t::check_invalid_draw_data() {
	if (invalidate_mats_mask & (1 << MAT_TYPE_SMALL  )) { // small objects
		mats_small   .invalidate();
		mats_amask   .invalidate();
		mats_text    .invalidate(); // Note: for now text is assigned to type MAT_TYPE_SMALL since it's always drawn with small objects
		mats_alpha_sm.invalidate();
	}
	if (invalidate_mats_mask & (1 << MAT_TYPE_STATIC )) { // large objects and 3D models
		mats_static  .invalidate(); // obj_model_insts will also be recreated
		mats_alpha   .invalidate();
		mats_exterior.invalidate(); // not needed since this is immutable?
	}
	//if (invalidate_mats_mask & (1 << MAT_TYPE_TEXT  )) {mats_text    .invalidate();} // text objects
	if (invalidate_mats_mask & (1 << MAT_TYPE_DYNAMIC )) {mats_dynamic .invalidate();} // dynamic objects
	if (invalidate_mats_mask & (1 << MAT_TYPE_DOORS   )) {mats_doors   .invalidate();} // door_handles will also be created
	if (invalidate_mats_mask & (1 << MAT_TYPE_LIGHTS  )) {mats_lights  .invalidate();}
	if (invalidate_mats_mask & (1 << MAT_TYPE_DETAIL  )) {mats_detail  .invalidate();}
	invalidate_mats_mask = 0; // reset for next frame
}
void building_room_geom_t::invalidate_draw_data_for_obj(room_object_t const &obj, bool was_taken) {
	if (obj.is_dynamic() || (obj.type == TYPE_BUTTON && obj.in_elevator())) { // elevator buttons are drawn as dynamic objects
		update_dynamic_draw_data();
		return;
	}
	bldg_obj_type_t const type(was_taken ? get_taken_obj_type(obj) : get_room_obj_type(obj));
	if (type.lg_sm & 2)            {invalidate_small_geom ();} // small objects
	if (type.lg_sm & 1)            {invalidate_static_geom();} // large objects and 3D models
	if (type.is_model )            {invalidate_model_geom ();} // model
	else if (type.lg_sm == 0)      {invalidate_detail_geom();} // detail object
	if (obj.type == TYPE_MACHINE)  {invalidate_detail_geom();} // machine details (valves, gauges, vents)
	if (obj.type == TYPE_CEIL_FAN) {invalidate_lights_geom();} // invalidate the light on the fan as well
	if (obj.type == TYPE_CLOCK)    {update_dynamic_draw_data();}
	if (obj.type == TYPE_CUP && obj.is_nonempty()) {invalidate_small_geom();} // cup with coffee
}
// Note: called when adding, removing, or moving objects
void building_room_geom_t::update_draw_state_for_room_object(room_object_t const &obj, building_t &building, bool was_taken) {
	invalidate_draw_data_for_obj(obj, was_taken);
	//if (type.ai_coll) {building.invalidate_nav_graph();} // removing this object should not affect the AI navigation graph
	modified_by_player = 1; // flag so that we avoid re-generating room geom if the player leaves and comes back
}
unsigned building_room_geom_t::get_num_verts() const {
	return mats_static .count_all_verts() +
		mats_small     .count_all_verts() +
		mats_text      .count_all_verts() +
		mats_detail    .count_all_verts() +
		mats_dynamic   .count_all_verts() +
		mats_lights    .count_all_verts() +
		mats_amask     .count_all_verts() +
		mats_alpha     .count_all_verts() +
		mats_alpha_sm  .count_all_verts() +
		mats_doors     .count_all_verts() +
		mats_exterior  .count_all_verts() +
		mats_ext_detail.count_all_verts();
}

building_materials_t &building_room_geom_t::get_building_mat(tid_nm_pair_t const &tex, bool dynamic, unsigned small, bool transparent, bool exterior) {
	assert(!(dynamic && exterior));
	assert(small <= 3); // 0=mats_static, 1=mats_small, 2=mats_detail, 3=doors
	if (transparent) {assert(!dynamic && !exterior);} // transparent objects must be static and can't be exterior
	if (transparent) return (small ? mats_alpha_sm   : mats_alpha);
	if (exterior)    return (small ? mats_ext_detail : mats_exterior);
	if (dynamic)     return mats_dynamic;
	if (small == 2)  return mats_detail;
	if (small == 3)  return mats_doors; // treated as a door
	if (small)       return ((tex.tid == FONT_TEXTURE_ID) ? mats_text : mats_small);
	return mats_static;
}
rgeom_mat_t &building_room_geom_t::get_metal_material(bool inc_shadows, bool dynamic, unsigned small, bool exterior, bool no_reflect,
	colorRGBA const &spec_color, float mag, float shine)
{
	tid_nm_pair_t tex(-1, 1.0f, inc_shadows, 0, no_reflect);
	tex.set_metal_specular(spec_color, mag, shine);
	tex.metalness = 1.0;
	return get_material(tex, inc_shadows, dynamic, small, 0, exterior);
}
tid_nm_pair_t get_scratched_metal_tex(float tscale, bool inc_shadows) {
	tid_nm_pair_t tex(get_texture_by_name("metals/60_scratch_metal.jpg"), tscale, inc_shadows);
	tex.set_metal_specular(WHITE, 0.5, 40.0); // less shiny
	tex.metalness = 0.5; // half reflective
	return tex;
}
rgeom_mat_t &building_room_geom_t::get_scratched_metal_material(float tscale, bool inc_shadows, bool dynamic, unsigned small, bool exterior) {
	return get_material(get_scratched_metal_tex(tscale, inc_shadows), inc_shadows, dynamic, small, 0, exterior);
}

void room_object_t::set_as_bottle(unsigned rand_id, unsigned max_type, bool no_empty, unsigned exclude_mask, bool make_empty, bool allow_transparent) {
	assert(max_type > 0 && max_type < NUM_BOTTLE_TYPES);
	obj_id = (uint16_t)rand_id;
	// cycle with a prime number until a valid type is selected; it's up to the caller to not exclude everything and make this infinite loop
	while (get_bottle_type() > max_type || ((1 << get_bottle_type()) & exclude_mask)) {obj_id += 13;}
	if      (no_empty  ) {obj_id &= 127;} // strip off second empty bit
	else if (make_empty) {obj_id |= BOTTLE_EMPTY_MASK;}
	bottle_params_t const &bp(bottle_params[get_bottle_type()]);
	color = (is_bottle_empty() ? bp.empty_color : bp.glass_color);
	if (allow_transparent && bp.transparent && is_bottle_empty()) {color.A = 0.33;} // make transparent
}

void building_room_geom_t::create_static_vbos(building_t const &building) {
	float const tscale(2.0/obj_scale);
	tid_nm_pair_t const &wall_tex(building.get_material().wall_tex); // interior wall texture
	static vect_room_object_t rugs;
	rugs.clear();

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible() || i->is_dynamic()) continue; // skip invisible and dynamic objects
		if (!i->is_strictly_normalized()) {cout << "Denormalized object of type " << unsigned(i->type) << ": " << i->str() << endl; assert(0);}
		assert(i->type < NUM_ROBJ_TYPES);

		switch (i->type) {
		case TYPE_TABLE:   add_table   (*i, tscale, 0.1, 0.08); break; // top_dz=10% of height, leg_width=8% of height
		case TYPE_CHAIR:   add_chair   (*i, tscale); break;
		case TYPE_STAIR:   add_stair   (*i, tscale, tex_origin, 0); break; // is_small_pass=0
		case TYPE_STAIR_WALL: add_stairs_wall(*i, tex_origin, wall_tex); break;
		case TYPE_RUG:     rugs.push_back(*i); break; // must be drawn last - save until later
		case TYPE_PICTURE: add_picture (*i); break;
		case TYPE_WBOARD:  add_picture (*i); break;
		case TYPE_BOOK:    add_book    (*i, 1, 0, 0); break; // lg
		case TYPE_BCASE:   add_bookcase(*i, 1, 0, 0, tscale, 0); break; // lg
		case TYPE_WINE_RACK: add_wine_rack(*i, 1, 0, tscale); break;
		case TYPE_DESK:    add_desk    (*i, tscale, 1, 0); break;
		case TYPE_RDESK:   add_reception_desk(*i, tscale); break;
		case TYPE_BED:     add_bed     (*i, 1, 0, tscale); break;
		case TYPE_WINDOW:  add_window  (*i, tscale); break;
		case TYPE_TUB:     add_tub_outer (*i); break;
		case TYPE_SINK:    add_sink_water(*i); break;
		case TYPE_TV: case TYPE_MONITOR: add_tv_picture(*i); break;
		case TYPE_CUBICLE: add_cubicle (*i, tscale); break;
		case TYPE_STALL:   add_br_stall(*i, 1, 0); break; // inc_lg=1, inc_sm=0
		case TYPE_COUNTER: add_counter (*i, tscale, 1, 0); break; // lg
		case TYPE_CABINET: add_cabinet (*i, *i, tscale, 1, 0); break; // parent is self; lg
		case TYPE_KSINK:   add_counter (*i, tscale, 1, 0); break; // counter with kitchen  sink; lg
		case TYPE_BRSINK:  add_counter (*i, tscale, 1, 0); break; // counter with bathroom sink; lg
		case TYPE_VANITY:  add_counter (*i, tscale, 1, 0); break; // counter with bathroom sink; lg
		case TYPE_PLANT:   add_potted_plant(*i, 1, 0); break; // pot only
		case TYPE_TREE:    add_tree(*i, 1, 0); break; // pot only
		case TYPE_DRESSER: case TYPE_NIGHTSTAND: add_dresser(*i, tscale, 1, 0); break;
		case TYPE_DRESS_MIR: add_dresser_mirror(*i, tscale); break;
		case TYPE_FLOORING:add_flooring(*i, tscale); break;
		case TYPE_CLOSET:  add_closet  (*i, wall_tex, building.get_trim_color(), 1, 0); break; // inc_lg=1, inc_sm=0
		case TYPE_MIRROR:  add_mirror  (*i); break;
		case TYPE_MED_CAB: add_med_cab (*i); break;
		case TYPE_SHOWER:  add_shower  (*i, tscale, 1, 0); break; // inc_lg=1, inc_sm=0
		case TYPE_SHOWERTUB: add_shower_tub(*i, wall_tex, building.get_trim_color(), tscale, 1, 0); break; // inc_lg=1, inc_sm=0
		case TYPE_BLINDS:  add_blinds  (*i); break;
		case TYPE_FPLACE:  add_fireplace(*i, tscale); break;
		case TYPE_FCABINET: add_filing_cabinet(*i, 1, 0); break; // lg
		case TYPE_SIGN:    add_sign    (*i, 1, 1, 1); break; // lg, exterior_only=1
		case TYPE_PIPE:    add_pipe(*i, 1); break; // add_exterior=1
		case TYPE_WIND_SILL: add_window_sill  (*i); break;
		case TYPE_EXT_STEP:  add_exterior_step(*i); break;
		case TYPE_BALCONY:   add_balcony(*i, building.ground_floor_z1, building.is_in_city); break;
		case TYPE_FALSE_DOOR: add_false_door(*i); break;
		case TYPE_RAILING:   if (i->is_exterior()) {add_railing(*i);}  break; // exterior only
		case TYPE_DOWNSPOUT: add_downspout(*i); break;
		case TYPE_SHELFRACK: add_rack(*i, 1, 0); break; // add_rack=1, add_objs=0
		case TYPE_CHIM_CAP:  add_chimney_cap(*i); break;
		case TYPE_LADDER:    add_ext_ladder(*i); break;
		case TYPE_CHECKOUT:  add_checkout(*i, tscale); break;
		case TYPE_FISHTANK:  add_fishtank(*i); break;
		case TYPE_OFF_PILLAR:add_wall_or_pillar(*i, tex_origin, wall_tex); break;
		case TYPE_SHELF_WALL:add_shelf_wall(*i, wall_tex); break;
		case TYPE_CONF_TABLE:add_conference_table(*i, tscale); break;
		case TYPE_INT_WINDOW:add_int_window(*i); break;
		case TYPE_BUCKET:    add_bucket(*i, 0, 1); break; // draw_metal=0, draw_liquid=1
		case TYPE_DWASHER:   add_dishwasher(*i); break;
		case TYPE_IBEAM:     add_ibeam(*i); break;
		case TYPE_CHEM_TANK: add_chem_tank(*i, 0); break; // draw_label=0
		case TYPE_HVAC_UNIT: add_hvac_unit(*i); break;
		case TYPE_VENT_FAN:  add_vent_fan_frame(*i); break;
		case TYPE_VENDING:   add_vending_machine(*i); break;
		case TYPE_LOCKER:    add_locker(*i); break;
		case TYPE_PARK_GATE: add_parking_gate (*i); break;
		case TYPE_CONV_BELT: add_conveyor_belt(*i, 0); break; // draw_dynamic=0
		case TYPE_COAT_RACK: add_coat_rack(*i, tscale); break;
		//case TYPE_FRIDGE: if (i->is_open()) {} break; // draw open fridge?
		case TYPE_ELEVATOR: break; // not handled here
		case TYPE_BLOCKER:  break; // not drawn
		case TYPE_COLLIDER: break; // not drawn
		default: break;
		} // end switch
	} // for i
	for (escalator_t const &e : building.interior->escalators) {add_escalator(e, building.get_window_vspace(), 0);} // draw_dynamic=0
	add_skylights_details(building);
	for (room_object_t &rug : rugs) {add_rug(rug);} // rugs are added last so that alpha blending of their edges works
	// Note: verts are temporary, but cubes are needed for things such as collision detection with the player and ray queries for indir lighting
	//highres_timer_t timer2("Gen Room Geom VBOs"); // < 2ms
	mats_static  .create_vbos(building);
	mats_alpha   .create_vbos(building);
	mats_exterior.create_vbos(building); // Note: ideally we want to include window dividers from trim_objs, but that may not have been created yet
}

void building_room_geom_t::create_small_static_vbos(building_t const &building, bool skip_mall_objs) {
	//highres_timer_t timer("Gen Room Geom Small"); // up to 36ms on new computer for buildings with large retail areas
	float const floor_ceil_gap(building.get_floor_ceil_gap());
	colorRGBA const &trim_color(building.get_trim_color());
	bldg_industrial_info_t const *ind_info(building.interior->ind_info.get());
	model_objs.clear(); // currently model_objs are only created for small objects in drawers, so we clear this here
	add_small_static_objs_to_verts(expanded_objs, 0, 0, trim_color, 0, floor_ceil_gap, ind_info); // inc_text=0

	if (skip_mall_objs && first_mall_obj_ix > 0) { // draw pre-mall objects and then elevator buttons, etc.
		add_small_static_objs_to_verts(objs, 0, first_mall_obj_ix,          trim_color, 0, floor_ceil_gap, ind_info); // inc_text=0
		add_small_static_objs_to_verts(objs, last_mall_obj_ix, objs.size(), trim_color, 0, floor_ceil_gap, ind_info); // inc_text=0
	}
	else { // draw full range
		add_small_static_objs_to_verts(objs, 0, 0, trim_color, 0, floor_ceil_gap, ind_info); // inc_text=0
	}
	add_attic_interior_and_rafters(building, 2.0/obj_scale, 0); // only if there's an attic; detail_pass=0
	for (tunnel_seg_t const &t : building.interior->tunnels) {add_tunnel(t);}
}

void building_room_geom_t::add_nested_objs_to_verts(vect_room_object_t const &objs_to_add) {
	vector_add_to(objs_to_add, pending_objs); // nested objects are added at the end so that small and text materials are thread safe
}
void building_room_geom_t::add_small_static_objs_to_verts(vect_room_object_t const &objs_to_add, unsigned six, unsigned eix, colorRGBA const &trim_color,
	bool inc_text, float floor_ceil_gap, bldg_industrial_info_t const *ind_info)
{
	if (objs_to_add.empty()) return; // don't add untextured material, otherwise we may fail the (num_verts > 0) assert
	float const tscale(2.0/obj_scale);
	if (eix == 0) {eix = objs_to_add.size();}

	for (unsigned i = six; i < eix; ++i) { // Note: iterating with indices to avoid invalid ref when add_nested_objs_to_verts() is called
		room_object_t const &c(objs_to_add[i]);
		if (!c.is_visible() || c.is_dynamic()) continue; // skip invisible and dynamic objects
		if (!c.is_strictly_normalized()) {cerr << "Denormalized object of type " << int(c.type) << ": " << c.str() << endl; assert(0);}
		assert(c.type < NUM_ROBJ_TYPES);

		switch (c.type) {
		case TYPE_STAIR:     add_stair    (c, tscale, tex_origin, 1); break; // is_small_pass=1
		case TYPE_BOOK:      add_book     (c, 0, 1, inc_text); break; // sm, maybe text
		case TYPE_BCASE:     add_bookcase (c, 0, 1, inc_text, tscale, 0); break; // sm, maybe text
		case TYPE_BED:       add_bed      (c, 0, 1, tscale); break;
		case TYPE_DESK:      add_desk     (c, tscale, 0, 1); break;
		case TYPE_DRESSER: case TYPE_NIGHTSTAND: add_dresser(c, tscale, 0, 1); break;
		case TYPE_TCAN:      add_trashcan (c); break;
		case TYPE_BUCKET:    add_bucket   (c, 1, 0); break; // draw_metal=1, draw_liquid=0
		case TYPE_SIGN:      add_sign     (c, 1, inc_text); break; // sm, maybe text
		case TYPE_CLOSET:    add_closet   (c, tid_nm_pair_t(), trim_color, 0, 1); break; // add closet wall trim and interior objects, don't need wall_tex; inc_lg=0, inc_sm=1
		case TYPE_SHOWER:    add_shower   (c, tscale, 0, 1); break; // inc_lg=0, inc_sm=1
		case TYPE_O_SHOWER:  add_shower   (c, tscale, 0, 1); break; // inc_lg=0, inc_sm=1
		case TYPE_SHOWERTUB: add_shower_tub(c, tid_nm_pair_t(), trim_color, tscale, 0, 1); break; // don't need wall_tex; inc_lg=0, inc_sm=1
		case TYPE_RAILING:   if (!c.is_exterior()) {add_railing(c);}  break; // interior only
		case TYPE_PLANT:     add_potted_plant(c, 0, 1); break; // plant only
		case TYPE_TREE:      add_tree(c, 0, 1); break; // tree only
		case TYPE_CRATE:     add_crate    (c); break; // not small but only added to windowless rooms
		case TYPE_BOX:       add_box      (c); break; // not small but only added to windowless rooms
		case TYPE_SHELVES:   add_shelves  (c, tscale); break; // not small but only added to windowless rooms
		case TYPE_SHELFRACK: add_rack(c, 0, 1, 0); break; // add_rack=0, add_objs=1, obj_text_pass=0
		case TYPE_MWAVE:     add_mwave    (c); break;
		case TYPE_COMPUTER:  add_computer (c); break;
		case TYPE_KEYBOARD:  add_keyboard (c); break;
		case TYPE_WINE_RACK: add_wine_rack(c, 0, 1, tscale); break;
		case TYPE_BOTTLE:    add_bottle   (c); break;
		case TYPE_DRINK_CAN: add_drink_can(c); break;
		case TYPE_VASE:      add_vase     (c); break;
		case TYPE_URN:       add_vase     (c); break;
		case TYPE_PAPER:     add_paper    (c); break;
		case TYPE_PAINTCAN:  add_paint_can(c); break;
		case TYPE_PEN: case TYPE_PENCIL: case TYPE_MARKER: add_pen_pencil_marker(c); break;
		case TYPE_TESTTUBE:  add_test_tube (c); break;
		case TYPE_LG_BALL:   add_lg_ball   (c); break;
		case TYPE_HANGER_ROD:add_hanger_rod(c); break;
		case TYPE_DRAIN:     add_drain_pipe(c); break;
		case TYPE_ELEC_WIRE: add_electrical_wire_pair(c); break;
		case TYPE_KEY:       if (has_key_3d_model()) {model_objs.push_back(c);} else {add_key(c);} break; // draw or add as 3D model
		case TYPE_SILVER:    if (c.was_expanded()  ) {model_objs.push_back(c);} break; // only draw here if expanded
		case TYPE_FOLD_SHIRT:if (c.was_expanded()  ) {model_objs.push_back(c);} break; // only draw here if expanded
		case TYPE_MONEY:     add_money   (c); break;
		case TYPE_PHONE:     add_phone   (c); break;
		case TYPE_TPROLL:    add_tproll  (c); break;
		case TYPE_TAPE:      add_tape    (c); break;
		case TYPE_SPRAYCAN:  add_spraycan(c); break;
		case TYPE_CRACK:     add_crack   (c); break;
		case TYPE_SWITCH:    add_switch  (c, 0); break; // draw_detail_pass=0
		case TYPE_BREAKER:   add_breaker (c); break;
		case TYPE_PLATE:     add_plate   (c); break;
		case TYPE_LAPTOP:    add_laptop  (c); break;
		case TYPE_PIZZA_BOX: add_pizza_box(c); break;
		case TYPE_PIZZA_TOP: add_pizza_top(c); break;
		case TYPE_BUTTON:    if (!c.in_elevator()) {add_button(c, 1, 0);} break; // skip buttons inside elevators, which are drawn as dynamic objects; inc_geom=1, inc_text=0
		case TYPE_LBASKET:   add_laundry_basket(c); break;
		case TYPE_TOASTER:   add_toaster_proxy (c); break;
		case TYPE_WHEATER:   add_water_heater  (c); break;
		case TYPE_FURNACE:   add_furnace       (c); break;
		case TYPE_BRK_PANEL: add_breaker_panel (c); break; // only added to basements
		case TYPE_ATTIC_DOOR:add_attic_door(c, tscale); break;
		case TYPE_RAMP:      add_pg_ramp   (c, tscale); break;
		case TYPE_TOY:       add_toy(c); break;
		case TYPE_PAN:       add_pan(c); break;
		case TYPE_COUNTER: add_counter (c, tscale, 0, 1); break; // sm
		case TYPE_KSINK:   add_counter (c, tscale, 0, 1); break; // sm
		case TYPE_CABINET: add_cabinet (c, c, tscale, 0, 1); break; // parent is self; sm
		case TYPE_VANITY:  add_counter (c, tscale, 0, 1); break; // sm
		case TYPE_FCABINET: add_filing_cabinet(c,  0, 1); break; // sm
		case TYPE_STAPLER: add_stapler(c); break;
		case TYPE_ERASER:  add_eraser (c); break;
		case TYPE_BAR_SOAP:   add_bar_soap(c); break;
		case TYPE_FEXT_MOUNT: add_fire_ext_mount(c); break;
		case TYPE_FEXT_SIGN:  add_fire_ext_sign (c); break;
		case TYPE_TEESHIRT:   add_teeshirt(c); break;
		case TYPE_PANTS:      add_pants   (c); break;
		case TYPE_BLANKET:    add_blanket (c); break;
		case TYPE_SERVER:     add_server  (c); break;
		case TYPE_POOL_BALL:  add_pool_ball(c); break;
		case TYPE_POOL_CUE:   add_pool_cue (c); break;
		case TYPE_WALL_MOUNT: add_wall_mount (c); break;
		case TYPE_POOL_TILE:  add_pool_tile  (c, tscale); break;
		case TYPE_POOL_FLOAT: add_pool_float(c); break;
		case TYPE_BENCH:      add_bench(c); break;
		case TYPE_DIV_BOARD:  add_diving_board(c); break;
		case TYPE_FLASHLIGHT: add_flashlight(c); break;
		case TYPE_CANDLE:     add_candle(c); break;
		case TYPE_CAMERA:     add_camera(c); break;
		case TYPE_CLOCK:      add_clock (c, 0); break; // add_dynamic=0
		case TYPE_FOOD_BOX:   add_food_box(c); break;
		case TYPE_SAFE:       add_safe(c); break;
		case TYPE_LAVALAMP:   add_lava_lamp(c); break;
		case TYPE_TRASH:      add_trash(c); break;
		case TYPE_METAL_BAR:  add_metal_bar(c); break;
		case TYPE_THEFT_SENS: add_theft_sensor(c); break;
		case TYPE_INT_LADDER: add_int_ladder(c); break;
		case TYPE_MACHINE:    add_machine(c, floor_ceil_gap, ind_info); break;
		case TYPE_SPIWEB:     add_spider_web(c); break;
		case TYPE_PET_CAGE:   add_pet_cage(c); break;
		case TYPE_CATWALK:    add_catwalk(c); break;
		case TYPE_DUCT:       add_duct(c); break;
		case TYPE_WARN_LIGHT: add_warning_light(c); break;
		case TYPE_PALLET:     add_pallet(c); break;
		case TYPE_CHEM_TANK:  add_chem_tank(c, 1); break; // draw_label=1
		case TYPE_CUP:        add_cup_liquid(c); break;
		case TYPE_HARDHAT:    add_hard_hat(c); break;
		case TYPE_TOPHAT:     add_tophat(c); break;
		case TYPE_COMP_MOUSE: add_comp_mouse(c); break;
		case TYPE_JAIL_BARS:  add_jail_bars(c); break;
		case TYPE_GUN:        add_gun(c); break;
		case TYPE_STICK_NOTE: add_sticky_note(c); break;
		case TYPE_GYM_WEIGHT: add_gym_weight (c); break;
		case TYPE_FOOD_TRAY:  add_food_tray  (c); break;
		case TYPE_STALL:      add_br_stall(c, 0, 1); break; // inc_lg=0, inc_sm=1
		case TYPE_DBG_SHAPE:  add_debug_shape(c); break;
		default: break;
		} // end switch
	} // for i
}

void building_room_geom_t::create_text_vbos(bool skip_mall_objs) {
	//highres_timer_t timer("Gen Room Geom Text");
	if (skip_mall_objs && first_mall_obj_ix > 0) { // draw pre-mall objects and then elevator buttons, etc.
		add_text_objs_to_verts(objs, 0, first_mall_obj_ix);
		add_text_objs_to_verts(objs, last_mall_obj_ix, objs.size());
	}
	else { // draw full range
		add_text_objs_to_verts(objs);
	}
	add_text_objs_to_verts(expanded_objs);
}
void building_room_geom_t::add_text_objs_to_verts(vect_room_object_t const &objs_to_add, unsigned six, unsigned eix) {
	if (eix == 0) {eix = objs_to_add.size();}

	for (unsigned i = six; i < eix; ++i) { // Note: iterating with indices to avoid invalid ref when add_nested_objs_to_verts() is called
		room_object_t const &c(objs_to_add[i]);
		if (!c.is_visible()) continue; // skip invisible objects
		switch (c.type) {
		case TYPE_BOOK:   add_book    (c, 0, 0, 1); break; // text only
		case TYPE_BCASE:  add_bookcase(c, 0, 0, 1, 1.0, 0); break; // text only
		case TYPE_SIGN:   add_sign    (c, 0, 1); break; // text only
		case TYPE_BUTTON: add_button  (c, 0, 1); break; // inc_geom=0, inc_text=1
		case TYPE_SHELFRACK:  add_rack(c, 0, 1, 1); break; // add_rack=0, add_objs=1, obj_text_pass=1
		default: break;
		} // end switch
	} // for i
}

void building_room_geom_t::create_detail_vbos(building_t const &building) {
	// currently only small objects that are non-interactive and can't be taken; TYPE_SWITCH almost counts; also, anything in the basement not seen from outside the building
	auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators
	tid_nm_pair_t const &wall_tex(building.get_material().wall_tex); // interior wall texture

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_visible()) continue;

		switch (i->type) {
		case TYPE_OUTLET:     add_outlet(*i); break;
		case TYPE_VENT:       add_vent  (*i); break;
		case TYPE_SWITCH:     add_switch(*i, 1); break; // draw_detail_pass=0
		case TYPE_PG_WALL:    add_wall_or_pillar (*i, tex_origin, wall_tex); break;
		case TYPE_PG_PILLAR:  add_basement_pillar(*i, wall_tex); break;
		case TYPE_PG_BEAM:    add_basement_beam  (*i, wall_tex); break;
		case TYPE_PARK_SPACE: add_parking_space  (*i, wall_tex.tscale_x); break;
		case TYPE_PIPE:       add_pipe(*i, 0); break; // add_exterior=0
		case TYPE_SPRINKLER:  add_sprinkler(*i); break;
		case TYPE_VALVE:      add_valve(*i); break;
		case TYPE_GAUGE:      add_gauge(*i); break;
		case TYPE_CURB:       add_curb (*i); break;
		case TYPE_CHIMNEY:    add_chimney(*i, building.get_material().side_tex); break; // uses exterior wall texture
		default: break;
		} // end switch
	} // for i
	for (auto const &i : trim_objs) {
		assert(i.type == TYPE_WALL_TRIM);
		add_wall_trim(i);
	}
	add_attic_interior_and_rafters(building, 2.0/obj_scale, 1); // only if there's an attic; detail_pass=1
	mats_detail    .create_vbos(building);
	mats_ext_detail.create_vbos(building);
}

vector3d get_obj_model_rotated_dir(room_object_t const &obj, building_t const *const building) {
	vector3d dir(obj.get_dir());
	if (!obj.rotates()) return dir;
	float const angle(123.4*obj.x1() + 456.7*obj.y1() + 567.8*obj.z1()); // random rotation angle based on position
	vector3d const rand_dir(vector3d(sin(angle), cos(angle), 0.0).get_norm());
	dir = ((dot_product(rand_dir, dir) < 0.0) ? -rand_dir : rand_dir); // random, but facing in the correct general direction

	if (obj.type == TYPE_RCHAIR && building) { // rotate to face the center of the room
		vector3d const center_dir(building->get_room(obj.room_id).get_cube_center() - obj.get_cube_center());
		if (SIGN(dir.x) != SIGN(center_dir.x)) {dir.x *= -1.0;}
		if (SIGN(dir.y) != SIGN(center_dir.y)) {dir.y *= -1.0;}
	}
	return dir;
}

void building_room_geom_t::create_obj_model_insts(building_t const &building) { // handle drawing of 3D models
	//highres_timer_t timer("Gen Room Model Insts");
	map<unsigned, vector3d> saved_office_chair_dirs;

	for (obj_model_inst_t const &i : obj_model_insts) { // save any office chair rotations from chairs the player spun
		if (i.dir.x == 0.0 || i.dir.y == 0.0) continue; // axis aligned, not random rotation
		if (get_room_object_by_index(i.obj_id).type == TYPE_OFF_CHAIR) {saved_office_chair_dirs[i.obj_id] = i.dir;}
	}
	obj_model_insts.clear();
	bool const is_residential(building.is_residential());

	for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
		auto const &obj_vect((vect_id == 1) ? expanded_objs : objs);
		unsigned const obj_id_offset((vect_id == 1) ? objs.size() : 0);
		auto objs_end((vect_id == 1) ? expanded_objs.end() : get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = obj_vect.begin(); i != objs_end; ++i) {
			if (!i->is_visible() || !i->is_obj_model_type()) continue;
			if (i->type == TYPE_KEY || (i->type == TYPE_SILVER && i->was_expanded())) continue; // drawn as small object model
			vector3d dir(get_obj_model_rotated_dir(*i, &building));
			if (building.is_rotated()) {building.do_xy_rotate_normal(dir);}
			unsigned const obj_id(i - obj_vect.begin() + obj_id_offset);
			// don't draw objects in interior rooms if the player is outside the building (useful for office bathrooms);
			// draw interior objects in residential buildings, such as plumbing fixtures visible through open interior bathroom doors, but not if in basement
			bool const int_vis_only(i->is_interior() && (!is_residential || i->z1() < building.ground_floor_z1));
			obj_model_insts.emplace_back(obj_id, dir, int_vis_only);

			if (i->type == TYPE_OFF_CHAIR) { // apply saved office chair rotations
				auto it(saved_office_chair_dirs.find(obj_id));
				if (it != saved_office_chair_dirs.end()) {obj_model_insts.back().dir = it->second;}
			}
			//get_untextured_material().add_cube_to_verts_untextured(get_true_room_obj_bcube(*i), WHITE); // for debugging of model bcubes
		} // for i
	} // for vect_id
}

void building_room_geom_t::create_lights_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Light"); // 0.75ms
	float const tscale(2.0/obj_scale);
	auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->is_visible() && i->type == TYPE_LIGHT) {add_light(*i, tscale);}
		else if (i->type == TYPE_CEIL_FAN) {add_ceiling_fan_light(*i, get_room_object_by_index(i->obj_id));} // pass in both the fan and the light
	}
	mats_lights.create_vbos(building);
}

void building_room_geom_t::create_dynamic_vbos(building_t const &building, point const &camera_bs, vector3d const &xlate, bool play_clock_tick) {
	//highres_timer_t timer(string("Gen Room Geom Dynamic ") + (building.is_house ? "house" : "office"));
	float const clock_sound_dist(10.0*CAMERA_RADIUS), floor_spacing(building.get_window_vspace());
	
	// is it better to have a rgeom type just for clocks that gets updated when the second/minute changes?
	// unclear if this would help, since we need to iterate over objs in either case, and that may be more expensive than drawing (and would be shared the current way);
	// plus when we get here we often want to update both dynamic objects and clocks anyway
	if (!obj_dstate.empty() || have_clock || building_alarm_active || have_conv_belt) { // we have an object with dynamic state, a dynamic clock, an alarm, or a conveyor belt
		auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (!i->is_visible()) continue; // only visible objects; can't do VFC because this is not updated every frame
			
			if (i->type == TYPE_CLOCK) {
				add_clock(*i, 1); // add_dynamic=1
				if (!play_clock_tick || (i->item_flags & 1)) continue; // skip for digital clocks
				point const pos(i->get_cube_center());
				if (fabs(pos.z - camera_bs.z) > 0.75*floor_spacing) continue; // different floors
				float const dist(p2p_dist(camera_bs, pos));
				if (dist > clock_sound_dist) continue;
				if (!building.get_room(i->room_id).contains_pt(camera_bs)) continue; // only play ticking when the player and clock are in the same room (ignoring rotation)
				gen_sound_thread_safe(SOUND_CLICK, (pos + xlate), 0.5*(1.0 - dist/clock_sound_dist));
				continue;
			}
			if (i->type == TYPE_THEFT_SENS) {
				if (i->is_active()) {add_theft_sensor(*i, 1);} // alarm_mode=1
				continue;
			}
			if (i->type == TYPE_CONV_BELT) {
				add_conveyor_belt(*i, 1); // draw_dynamic=1
				continue;
			}
			if (!i->is_dynamic()) continue; // only dynamic objects
			switch (i->type) {
			case TYPE_LG_BALL  : add_lg_ball  (*i); break;
			case TYPE_POOL_BALL: add_pool_ball(*i); break;
			default: assert(0); // not a supported dynamic object type
			}
		} // for i
	}
	for (elevator_t const &e : building.interior->elevators) {
		assert(e.car_obj_id < objs.size());
		assert(e.button_id_start < e.button_id_end && e.button_id_end <= objs.size());
		float const fc_thick_scale(building.get_elevator_fc_thick_scale());
		add_elevator_doors(e, fc_thick_scale); // add dynamic elevator doors
		// draw elevator car for this elevator
		float const e_floor_spacing(building.get_elevator_floor_spacing(e));
		unsigned const floor_offset(building.calc_floor_offset(e.z1(), e_floor_spacing));
		add_elevator(objs[e.car_obj_id], e, 2.0/obj_scale, fc_thick_scale, floor_offset, e_floor_spacing,
			floor_spacing, building.has_parking_garage, !building.interior->elevators_disabled);

		for (auto j = objs.begin() + e.button_id_start; j != objs.begin() + e.button_id_end; ++j) {
			if (j->type == TYPE_BLOCKER || j->type == TYPE_ELEC_WIRE) continue; // button was removed?
			assert(j->type == TYPE_BUTTON);
			if (j->in_elevator()) {add_button(*j, 1, 0);} // add button as a dynamic object if it's inside the elevator; inc_geom=1, inc_text=0
		}
	} // for e
	for (escalator_t  const &e : building.interior->escalators) {add_escalator(e, floor_spacing, 1);} // draw_dynamic=1
	for (tunnel_seg_t const &t : building.interior->tunnels   ) {add_tunnel_water(t);}
	mats_dynamic.create_vbos(building);
}

colorRGBA building_t::get_door_handle_color() const {
	colorRGBA const house_handle_colors [3] = {LT_GRAY, GRAY, BRASS_C};
	colorRGBA const office_handle_colors[3] = {GRAY, DK_GRAY, BLACK};
	return (is_house ? house_handle_colors : office_handle_colors)[(interior->doors.size() + interior->rooms.size() + mat_ix) % 3];
}
void building_room_geom_t::create_door_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Doors"); // 0.1ms
	vect_door_t const &doors(building.interior->doors);
	bool const have_door_handle_model(global_building_params.add_door_handles && building_obj_model_loader.is_model_valid(OBJ_MODEL_DOOR_HANDLE));
	colorRGBA const handle_color(global_building_params.add_door_handles ? building.get_door_handle_color() : WHITE);
	bool const residential(building.is_residential()), has_br_tex(building.has_backrooms_texture());
	door_handles.clear();

	for (door_t const &d : doors) { // interior doors; opens_out=0, exterior=0
		door_rotation_t drot;

		if (d.for_jail) { // jail cell metal door
			add_jail_cell_door(d, building, drot);
		}
		else { // normal interior door
			bool const house_door(building.is_residential() ||
				(d.z1() < building.ground_floor_z1 && has_br_tex && !building.get_basement().contains_cube(d))); // backrooms uses house door
			uint8_t const door_type(house_door ? (uint8_t)tquad_with_ix_t::TYPE_HDOOR : (uint8_t)tquad_with_ix_t::TYPE_ODOOR);
			building.add_door_verts(d, *this, drot, door_type, d.dim, d.open_dir, d.open_amt, 0, 0,
				d.on_stairs, d.hinge_side, d.use_min_open_amt(), d.get_mult_floor()); // opens_out=0, exterior=0
		}
		maybe_add_door_sign(d, drot);
		if (!global_building_params.add_door_handles) continue;
		if (d.on_stairs) continue; // skip basement stairs doors since they're not drawn when open anyway
		
		if (d.for_jail) {} // handle drawn inside add_jail_cell_door
		else if (have_door_handle_model) { // add model to door_handles
			bool const handle_side(d.get_handle_side());
			float const handle_height(0.04*d.dz());
			point const door_center(d.get_cube_center());
			vector3d handle_dir;
			handle_dir[!d.dim] = (handle_side ? -1.0 : 1.0);
			point handle_center(door_center);
			handle_center.z += (residential ? 0.15 : -0.18)*handle_height;
			handle_center   -= 0.393*d.get_width()*handle_dir;
			float sin_term(0.0), cos_term(0.0);
			point pivot;
			
			if (d.open_amt > 0.0) { // rotate handle_dir and handle_center
				// similar to rotate_and_shift_door()
				float const rot_angle(-float(drot.angle)*TO_RADIANS*(d.hinge_side ? -1.0 : 1.0));
				sin_term = sin(rot_angle);
				cos_term = cos(rot_angle);
				pivot    = d.get_cube_center();
				pivot[!d.dim] = d.d[!d.dim][!handle_side];
				do_xy_rotate_normal(sin_term, cos_term, handle_dir);
			}
			for (unsigned side = 0; side < 2; ++side) {
				//if (d.for_jail == 1 && bool(side) != d.open_dir) continue; // no handle on the inside/cell side of the door
				point side_pos(handle_center);
				side_pos[d.dim] += (side ? 1.0 : -1.0)*0.68*handle_height;

				if (d.open_amt > 0.0) {
					do_xy_rotate(sin_term, cos_term, pivot, side_pos);
					side_pos[d.dim] += drot.shift;
				}
				bool const mirror(bool(side) ^ d.open_dir ^ d.hinge_side ^ 1);
				door_handles.emplace_back(side_pos, handle_height, d.dim, bool(side), mirror, (d.open_amt == 0.0), handle_dir);
			}
		}
		else {add_door_handle(d, drot, handle_color, residential);}
	} // for d
	if (building.has_mall()) {
		float const window_vspace(building.get_window_vspace()), wall_thickness(building.get_wall_thickness()), trim_thick(building.get_trim_thickness());

		for (store_doorway_t const &d : building.interior->mall_info->store_doorways) {
			if (d.open_amt == 1.0) continue; // open gate
			cube_t gate(d);
			gate.z2() = d.z2() - 0.1*window_vspace - trim_thick;
			gate.z1() = d.get_gate_z1() + trim_thick;
			if (gate.dz() <= 0.0) continue; // none visible
			gate.expand_in_dim( d.dim, -0.4*d.get_sz_dim(d.dim)); // shrink
			gate.expand_in_dim(!d.dim, -0.5*wall_thickness); // shrink (same as cbox)
			add_store_gate(gate, d.dim, d.open_amt);
		} // for d
	}
	mats_doors.create_vbos(building);
}

bool ceiling_fan_is_on(room_object_t &obj, vect_room_object_t const &objs) {
	if (!obj.is_powered() || !(obj.flags & RO_FLAG_ROTATING)) return 0; // only enabled for some fans
	if (obj.in_factory ()) return 1; // always on if rotating
	assert(obj.obj_id < objs.size());
	return (objs[obj.obj_id].type == TYPE_LIGHT && objs[obj.obj_id].is_light_on()); // fan is on if light is on
}
void rotate_dir_about_z(vector3d &dir, float angle) { // Note: assumes dir is normalized
	if (angle == 0.0) return;
	assert(dir.z == 0.0); // dir must be in XY plane
	float const new_angle(atan2(dir.y, dir.x) + angle);
	dir.assign(cosf(new_angle), sinf(new_angle), 0.0);
}
void apply_room_obj_rotate(room_object_t &obj, obj_model_inst_t &inst, vect_room_object_t const &objs) {
	if (!animate2 || !(obj.flags & RO_FLAG_ROTATING)) return;

	if (obj.type == TYPE_OFF_CHAIR) {
		if (office_chair_rot_rate == 0.0) {obj.flags &= ~RO_FLAG_ROTATING; return;} // if no longer rotating, clear rotation bit
		rotate_dir_about_z(inst.dir, office_chair_rot_rate*fticks);
	}
	else if (obj.type == TYPE_HANGER || obj.type == TYPE_CLOTHES || obj.type == TYPE_JUMPSUIT) {
		inst.dir = obj.get_dir(); // reset before applying rotate
		float const angle(((obj.flags & RO_FLAG_ADJ_LO) ? -1.0 : 1.0)*0.08*TWO_PI);
		rotate_dir_about_z(inst.dir, -angle); // limited rotation angle
	}
	else if (obj.type == TYPE_CEIL_FAN) { // rotate the entire model rather than just the fan blades
		if (ceiling_fan_is_on(obj, objs)) {rotate_dir_about_z(inst.dir, 0.2*fticks);} // fan is on if light is on
	}
	else {
		cerr << "Error: apply_room_obj_rotate() on unsupported object type " << unsigned(obj.type) << endl;
		assert(0); // unsupported object type
	}
}
float get_camera_z_rotate() {return -atan2(cview_dir.y, cview_dir.x);}

void building_room_geom_t::draw_interactive_player_obj(carried_item_t const &c, shader_t &s, vector3d const &xlate) { // held by the player
	static rgeom_mat_t mat; // allocated memory is reused across frames; VBO is recreated every time; untextured
	bool needs_blend(0), reset_mat_nm_tid(0);

	if (c.type == TYPE_SPRAYCAN || c.type == TYPE_MARKER) {
		room_object_t c_rot(c);
		c_rot.dir = 0; // facing up
		unsigned const dim(get_max_dim(c.get_size()));

		if (dim != 2) { // if not oriented in Z
			UNROLL_2X(swap(c_rot.d[dim][i_], c_rot.d[2][i_]);); // rotate into Z dir
			c_rot.translate(c.get_cube_center() - c_rot.get_cube_center()); // translate it back to the correct location
		}
		if (c.type == TYPE_SPRAYCAN) {add_spraycan_to_material(c_rot, mat, 1);} // draw_bottom=1
		else {add_pen_pencil_marker_to_material(c_rot, mat);}
	}
	else if (c.type == TYPE_TPROLL || c.type == TYPE_TAPE) { // apply get_player_cview_rot_matrix()?
		if (c.type == TYPE_TPROLL) {mat.tex.nm_tid = get_toilet_paper_nm_id(); reset_mat_nm_tid = 1;}
		add_vert_roll_to_material(c, mat, c.get_remaining_capacity_ratio(), 1); // player_held=1; unfortunately, we don't support texturing the TP roll here
		needs_blend = 1;
	}
	else if (c.type == TYPE_BOOK) {
		static building_room_geom_t tmp_rgeom;
		float const z_rot_angle(get_camera_z_rotate() - PI_TWO);
		tmp_rgeom.add_book(c, 1, 1, 1, 0.0, 0, 0, z_rot_angle); // draw lg/sm/text
		enable_blend(); // needed for book text
		tmp_rgeom.mats_small.upload_draw_and_clear(s);
		tmp_rgeom.mats_text .upload_draw_and_clear(s);
		disable_blend();
		return;
	}
	else if (c.type == TYPE_PHONE) {
		float const z_rot_angle(get_camera_z_rotate());
		cube_t const screen(add_phone_frame_and_return_screen_if_on(c, mat, 1)); // in_hand=1
		rotate_verts(mat.itri_verts, plus_z, z_rot_angle, c.get_cube_center(), 0); // rotate all quad verts about Z axis
		
		if (!screen.is_all_zeros()) {
			static rgeom_mat_t screen_mat;
			screen_mat.tex = get_phone_tex(c);
			screen_mat.add_cube_to_verts(screen, WHITE, all_zeros, ~EF_Z2, 0, 1); // mirror_x=1
			rotate_verts(screen_mat.quad_verts, plus_z, z_rot_angle, c.get_cube_center(), 0); // rotate all quad verts about Z axis
			tid_nm_pair_dstate_t state(s);
			screen_mat.upload_draw_and_clear(state);
		}
	}
	else if (c.type == TYPE_RAT) { // draw the rat facing away from the player
		bool const is_dead(c.is_broken()); // upside down if dead; shadow_pass=0; not animated
		building_obj_model_loader.draw_model(s, c.get_cube_center(), c, cview_dir, rat_color, xlate, OBJ_MODEL_RAT, 0, 0, nullptr, 0, 0, 0, is_dead);
		check_mvm_update();
		return; // don't need to run the code below
	}
	else if (c.is_medicine()) {
		static building_room_geom_t tmp_rgeom;
		room_object_t bottle(c);
		vector3d const sz(bottle.get_size());

		for (unsigned d = 0; d < 2; ++d) { // if bottle was on its side, draw it as upright
			if (sz[d] < sz.z) continue;
			float const val(0.5*(sz[d] - sz.z));
			bottle.expand_in_dim(2,  val);
			bottle.expand_in_dim(d, -val);
			break;
		}
		tmp_rgeom.add_bottle(bottle, 1, get_camera_z_rotate()); // add_bottom=1
		tmp_rgeom.mats_small.upload_draw_and_clear(s);
	}
	else if (c.type == TYPE_FIRE_EXT) {
		building_obj_model_loader.draw_model(s, c.get_cube_center(), c, cview_dir, WHITE, xlate, OBJ_MODEL_FIRE_EXT);
		check_mvm_update();
		return; // don't need to run the code below
	}
	else if (c.type == TYPE_CANDLE) {
		static building_room_geom_t tmp_rgeom;
		room_object_t obj(c);
		obj.z2() -= (1.0 - c.get_remaining_capacity_ratio())*0.9*c.dz(); // slowly burn down shorter over time
		tmp_rgeom.add_candle(obj);
		tmp_rgeom.mats_small.upload_draw_and_clear(s);
		
		if (c.is_lit()) { // add flame; will be drawn after building interior geom
			player_candle_pos = cube_top_center(obj);
			float const radius(0.8*c.get_radius()), height(4.0*radius);
			point const camera_bs(get_camera_pos() - xlate);
			point center(player_candle_pos + 0.2*height*plus_z);
			center += 0.1*radius*(camera_bs - center).get_norm(); // move slightly in front of the wick
			candle_qbd.add_animated_billboard(center, camera_bs, up_vector, WHITE, radius, 0.5*height, fract(2.0f*tfticks/TICKS_PER_SECOND));
			
			if (animate2) { // add smoke particles about once every 20 frames
				static rand_gen_t smoke_rgen;

				if ((smoke_rgen.rand() % 20) == 0) {
					particle_manager.add_particle((center + 0.5*height*plus_z), 0.0001*plus_z, colorRGBA(GRAY, 0.25), 1.2*radius, PART_EFFECT_SMOKE);
				}
			}
		}
	}
	else if (c.type == TYPE_ERASER) { // Note: could make erasers rounded if we had support for non-vertical rounded rectangles
		float const sz_delta(c.get_width() - c.dz()); // rotate vertically
		cube_t c_rot(c);
		
		if (c.dim) { // swap XY
			float const resize_amt(0.5*(c.dx() - c.dy()));
			c_rot.expand_in_x(-resize_amt);
			c_rot.expand_in_y( resize_amt);
		}
		c_rot.z2() += sz_delta;
		c_rot.expand_in_y(-0.5*sz_delta);
		mat.add_cube_to_verts_untextured(c_rot, c.color); // simple untextured cube; all sides drawn
		rotate_verts(mat.quad_verts, plus_z, (get_camera_z_rotate() + PI_TWO), c.get_cube_center(), 0); // rotate all itri verts about Z axis
	}
	else if (c.type == TYPE_FLASHLIGHT) {
		unsigned const dim(get_max_dim(c.get_size()));
		room_object_t c_rot(c);
		c_rot.dir = 1; // points outward
		c_rot.swap_dims(0, dim); // swap so that dim==X
		c_rot.translate(c.get_cube_center() - c_rot.get_cube_center()); // translate it back to the correct location
		mat.tex.set_metal_specular(); // shiny
		add_flashlight_to_material(c_rot, mat, N_CYL_SIDES);
		rotate_verts(mat.itri_verts, plus_z, get_camera_z_rotate(), c.get_cube_center(), 0); // rotate all itri verts about Z axis
	}
	else {assert(0);}
	if (needs_blend) {enable_blend();}
	tid_nm_pair_dstate_t state(s);
	mat.upload_draw_and_clear(state);
	if (reset_mat_nm_tid) {mat.tex.nm_tid = -1;}
	if (needs_blend) {disable_blend();}
	mat.tex.set_specular(0.0, 0.0); // clear specular
}

void draw_candle_flames() {
	draw_emissive_billboards(candle_qbd, FIRE_TEX);
}

class water_draw_t {
	rgeom_mat_t mat;
	float tex_off;
	rand_gen_t rgen;

	void apply_vert_texture(unsigned verts_start, float tscale, float speed) {
		speed *= tex_off;
		for (auto i = mat.itri_verts.begin() + verts_start; i != mat.itri_verts.end(); ++i) {i->t[1] *= tscale; i->t[1] += speed;}
	}
public:
	water_draw_t() : mat(rgeom_mat_t(tid_nm_pair_t(FOAM_TEX))), tex_off(0.0) {}

	void add_water_for_sink(room_object_t const &obj) {
		if (!obj.is_active()) return;
		bool const is_cube(obj.type == TYPE_KSINK || obj.type == TYPE_BRSINK || obj.type == TYPE_VANITY);
		float const dz(obj.dz());
		cube_t c;
		point pos(obj.get_cube_center());
		pos[obj.dim] += (obj.dir ? -1.0 : 1.0)*(is_cube ? 0.25 : 0.095)*obj.get_length(); // move toward the back of the sink
		c.set_from_sphere(pos, (is_cube ? 0.02 : 0.0055)*dz);
		set_cube_zvals(c, (obj.z1() + (is_cube ? 0.7 : 0.6)*dz), (obj.z1() + (is_cube ? 1.3 : 0.925)*dz));
		unsigned const verts_start(mat.itri_verts.size());
		mat.add_vcylin_to_verts(c, colorRGBA(WHITE, 0.5), 0, 0, 0, 0, 1.0, 1.0, 0.2);
		apply_vert_texture(verts_start, 1.2, 1.0);
	}
	void add_water_for_shower(room_object_t const &obj, unsigned obj_ix, particle_manager_t &particle_manager) {
		if (!obj.is_active()) return;
		point head_pos;
		vector3d head_dir;
		get_shower_head_pos_dir(obj, head_pos, head_dir);
		colorRGBA const color(WHITE, 0.4);
		unsigned const verts_start(mat.itri_verts.size()), num_streams(12);
		float const delta_angle(TWO_PI/num_streams), height(obj.dz()), max_dist(2.0*height), water_height(0.005*height);
		float const width(0.5*(obj.dx() + obj.dy())), head_radius(0.07*width), spread_radius(0.75*head_radius);
		float const spray_radius_top(0.04*head_radius), spray_radius_bot(0.12*head_radius), spread_radius_bot(7.0*spread_radius);
		vector3d d_ref(vector_from_dim_dir(obj.dim, obj.dir));
		vector3d const d1(cross_product(d_ref, head_dir).get_norm()), d2(cross_product(d1, head_dir).get_norm());
		float const spash_prob(min(1.0, 0.8/fticks)); // 32 per second
		unsigned const splash_stream_ix(rgen.rand() % round_fp(num_streams/spash_prob)); // may select an index larger than num_streams

		for (unsigned n = 0; n < num_streams; ++n) {
			float const theta(n*delta_angle);
			vector3d const offset(sin(theta)*d1 + cos(theta)*d2);
			point water_start(head_pos + spread_radius*offset), water_end(head_pos + spread_radius_bot*offset + max_dist*head_dir);
			do_line_clip(water_start, water_end, obj.d); // should only clip water_end
			point const ray_end(water_end + 2.0*spray_radius_bot*head_dir); // extend a bit further into the floor
			mat.add_cylin_to_verts(ray_end, water_start, spray_radius_bot, spray_radius_top, color, 0, 0); // no ends
			if (!animate2 || n != splash_stream_ix) continue; // no splash for this stream
			water_end.z += water_height; // place above the water plane
			particle_manager.add_particle(water_end, zero_vector, WHITE, 2.5*spray_radius_bot, PART_EFFECT_SPLASH);
		} // for n
		apply_vert_texture(verts_start, 4.0, 3.0);
		// draw a puddle on the floor
		point puddle_center(obj.xc(), obj.yc(), (obj.z1() + water_height));
		puddle_center[obj.dim] += (obj.dir ? -1.0 : 1.0)*0.15*width; // move further from the wall, closer to where the water hits
		float const puddle_radius((0.4 + 0.01*sin(TWO_PI*tex_off + 100.0*(puddle_center.x + puddle_center.y)))*width); // slightly changes over time, different for each shower
		mat.add_disk_to_verts(puddle_center, puddle_radius, plus_z, color);
		// add steam
		float const steam_prob(0.05/fticks);

		if (animate2 && rgen.rand_float() < steam_prob) {
			float const theta(rgen.rand_float()*TWO_PI), head_dist(rgen.rand_uniform(0.1, 0.9));
			point water_end(head_pos + max_dist*head_dir);
			do_line_clip(head_pos, water_end, obj.d); // should only clip water_end
			point steam_pos(head_pos + head_dist*(water_end - head_pos)); // somewhere along the water path
			steam_pos += (spread_radius + head_dist*(spread_radius_bot - spread_radius))*(sin(theta)*d1 + cos(theta)*d2); // random shift to the side
			particle_manager.add_particle(steam_pos, 0.0008*plus_z, colorRGBA(WHITE, 0.15), 0.2*width, PART_EFFECT_STEAM, obj_ix, 1.5); // coll_radius=1.5
		}
	}
	void draw_and_clear(shader_t &s) {
		if (mat.empty()) return;
		glDepthMask(GL_FALSE); // disable depth writing - fixes sky visible through exterior wall, but then not drawn in front of exterior wall
		enable_blend();
		tid_nm_pair_dstate_t state(s);
		mat.upload_draw_and_clear(state);
		disable_blend();
		glDepthMask(GL_TRUE); // re-enable depth writing
		if (animate2) {tex_off += 0.02*fticks;} // animate the texture
		if (tex_off > 1.0) {tex_off -= 1.0;}
	}
};

water_draw_t water_draw;

class lava_lamp_draw_t {
	bool closest_is_on=0;
	float dmin_sq=0.0, cur_time=0.0;
	quad_batch_draw qbd;
public:
	void add_lava_lamp(room_object_t const &obj, point const &camera_bs, building_t const &building) {
		// Note: must match code in building_room_geom_t::add_lava_lamp()
		float const radius(obj.get_radius()), r_top(0.5*radius), height(obj.dz());
		point center(obj.get_cube_center());
		if (building.is_rotated()) {building.do_xy_rotate(building.bcube.get_cube_center(), center);}
		vector3d const vdir((camera_bs - center).get_norm()), vside(cross_product(vdir, plus_z).get_norm());
		point pts[4] = {(center - radius*vside), (center + radius*vside), (center + r_top*vside), (center - r_top*vside)};
		pts[0].z = pts[1].z = obj.z1() + 0.42 *height; // z1
		pts[2].z = pts[3].z = obj.z2() - 0.167*height; // z2
		bool const is_on(obj.is_light_on());
		color_wrapper const cw(is_on ? WHITE : GRAY); // we can't do proper lighting here, so make the color gray if it's off
		qbd.add_quad_pts(pts, cw, vdir, tex_range_t());
		float const dsq(p2p_dist_sq(camera_bs, center));
		
		if (dmin_sq == 0.0 || dsq < dmin_sq) {
			dmin_sq = dsq;
			closest_is_on = is_on;
		}
	}
	void next_frame() {
		if (closest_is_on) {cur_time += fticks;} // update even if empty/not visible
		dmin_sq = 0.0; // reset for next frame
	}
	void draw_and_clear(shader_t &s) { // Note: not drawn in the shadow pass, so the interior part doesn't cast a shadow
		if (qbd.empty()) {
			if (cur_time > 10000) {cur_time = 0.0;} // reset time when not visible to avoid FP accuracy issues
			return;
		}
		shader_t lls;
		lls.set_vert_shader("no_lighting_tex_coord");
		lls.set_frag_shader("lava_lamp");
		lls.begin_shader();
		lls.add_uniform_float("time", cur_time/TICKS_PER_SECOND);
		qbd.draw_and_clear();
		s.make_current(); // switch back to the normal shader
	}
}; // lava_lamp_draw_t

lava_lamp_draw_t lava_lamp_draw;

int room_object_t::get_model_id() const { // Note: first 8 bits is model ID, last 8 bits is sub-model ID
	assert(type >= TYPE_TOILET);
	if (type == TYPE_MONITOR)    return OBJ_MODEL_TV        ; // monitor has same model as TV
	if (type == TYPE_GBIKE  )    return OBJ_MODEL_BICYCLE   ; // same model as city bicycle
	if (type == TYPE_XFORMER)    return OBJ_MODEL_SUBSTATION; // same model as city substation
	if (type == TYPE_US_FLAG)    return OBJ_MODEL_FLAG      ; // same model as city flag
	if (type == TYPE_BLDG_FOUNT) return OBJ_MODEL_FOUNTAIN + pack_sub_model_id(item_flags); // same models as city fountains; select a sub_model_id
	int id((int)type + OBJ_MODEL_TOILET - TYPE_TOILET);
	// choose a sub_model_id for these types using bits 8-15
	if (type == TYPE_HANGER || type == TYPE_CLOTHES || type == TYPE_PLANT_MODEL || type == TYPE_SHOE || type == TYPE_HOSP_BED || type == TYPE_APPLE || type == TYPE_EX_MACHINE) {
		id += pack_sub_model_id(item_flags);
	}
	return id;
}

void building_t::draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, shader_t &amask_shader, occlusion_checker_noncity_t &oc, vector3d const &xlate,
	unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building, bool mall_visible)
{
	if (!has_room_geom()) return;

	if (0 && (display_mode & 0x20) && !shadow_only && player_in_building && bcube.contains_pt(camera_pdu.pos - xlate)) { // debug visualization of light occluders
		vect_colored_cube_t cc;
		gather_interior_cubes(cc, get_bcube_inc_extensions());
		select_texture(WHITE_TEX);
		for (colored_cube_t const &c : cc) {s.set_cur_color(c.color); draw_simple_cube(c);}
		return;
	}
	interior->room_geom->draw(bbd, s, amask_shader, *this, oc, xlate, building_ix, shadow_only, reflection_pass, inc_small, player_in_building, mall_visible);
	//enable_blend();
	//glDisable(GL_CULL_FACE);
	//s.set_cur_color(colorRGBA(1.0, 0.0, 0.0, 0.5)); // for use with debug visualization
}
void building_t::gen_and_draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, shader_t &amask_shader, occlusion_checker_noncity_t &oc, vector3d const &xlate,
	unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building, bool ext_basement_conn_visible, bool mall_visible)
{
	if (!interior) return;
	if (!global_building_params.enable_rotated_room_geom && is_rotated()) return; // rotated buildings: need to fix texture coords, room object collisions, mirrors, etc.

	if (!shadow_only && !player_in_building && !mall_visible && !camera_pdu.point_visible_test(bcube.get_cube_center() + xlate)) {
		// skip if none of the building parts are visible to the camera; this is rare, so it may not help
		bool any_part_visible(0);

		for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
			if (camera_pdu.cube_visible(*p + xlate)) {any_part_visible = 1; break;}
		}
		// check the extended basement bcube if visible through a connecting room
		any_part_visible |= (ext_basement_conn_visible && has_ext_basement() && camera_pdu.cube_visible(interior->basement_ext_bcube + xlate));
		any_part_visible |= point_in_attic(camera_pdu.pos - xlate); // check the attic

		if (!any_part_visible) {
			// exterior geometry such as roof vents and chimney caps may still be visible even if no parts are visible, so draw them
			if (is_house && !reflection_pass && has_room_geom()) {interior->room_geom->mats_exterior.draw(bbd, s, 0, 0, 1);}
			return;
		}
	}
	if (!has_room_geom()) {
		interior->room_geom.reset(new building_room_geom_t(bcube.get_llc()));
		// capture state before generating backrooms, which may add more doors
		interior->room_geom->init_num_doors   = interior->doors      .size();
		interior->room_geom->init_num_dstacks = interior->door_stacks.size();
		rand_gen_t rgen;
		rgen.set_state(building_ix, (parts.size() + 17*interior->rgen_seed_ix)); // set to something canonical per building
		interior->room_geom->decal_manager.rgen = rgen; // copy rgen for use with decals
		gen_room_details(rgen, building_ix); // generate so that we can draw it
		assert(has_room_geom());
	}
	if (has_room_geom() && (inc_small == 2 || inc_small == 3)) {add_wall_and_door_trim_if_needed();} // gen trim (exterior and interior) when close to the player
	draw_room_geom(bbd, s, amask_shader, oc, xlate, building_ix, shadow_only, reflection_pass, inc_small, player_in_building, mall_visible);
}
void building_t::clear_room_geom(bool even_if_player_modified) {
	clear_clothing_textures();
	if (!has_room_geom()) return;

	if (interior->room_geom->modified_by_player && !even_if_player_modified) { // keep the player's modifications and don't delete the room geom
		interior->room_geom->clear_materials(); // but we can still clear the materials
		return;
	}
	// restore pre-room_geom door state by removing any doors added to backrooms
	assert(interior->room_geom->init_num_doors   <= interior->doors      .size());
	assert(interior->room_geom->init_num_dstacks <= interior->door_stacks.size());
	interior->doors      .resize(interior->room_geom->init_num_doors  );
	interior->door_stacks.resize(interior->room_geom->init_num_dstacks);
	interior->room_geom->clear(); // free VBO data before deleting the room_geom object
	interior->room_geom.reset();
	invalidate_nav_graph(); // required since interior doors may be removed
}
void building_t::clear_small_room_geom_vbos() {
	if (this == player_building) return; // not for the player building
	if (has_room_geom()) {interior->room_geom->clear_small_materials();}
}
void building_t::clear_and_regen_new_seed() {
	clear_room_geom(1); // even_if_player_modified=1
	if (interior) {++interior->rgen_seed_ix;}
}

void get_stove_burner_locs(room_object_t const &stove, point locs[4]) {
	vector3d const sz(stove.get_size());
	bool const dim(stove.dim), dir(stove.dir);
	float const zval(stove.z2() - 0.23*sz.z), dsign(dir ? -1.0 : 1.0);

	for (unsigned w = 0; w < 2; ++w) { // width dim
		float const wval(stove.d[!dim][0] + ((bool(w) ^ dim ^ dir ^ 1) ? 0.72 : 0.28)*sz[!dim]); // left/right burners

		for (unsigned d = 0; d < 2; ++d) { // depth dim
			float const dval(stove.d[dim][dir] + dsign*(d ? 0.66 : 0.375)*sz[dim]); // front/back burners
			point pos(dval, wval, zval);
			if (dim) {swap(pos.x, pos.y);}
			locs[2*w + d] = pos;
		} // for d
	} // for w
}
void draw_stove_flames(room_object_t const &stove, point const &camera_bs, shader_t &s) {
	if (stove.item_flags == 0) return; // no burners on
	static quad_batch_draw flame_qbd; // reused across frames
	point locs[4];
	get_stove_burner_locs(stove, locs);
	float const dsign(stove.dir ? -1.0 : 1.0);

	for (unsigned w = 0; w < 2; ++w) { // width dim
		float const radius((w ? 0.09 : 0.07)*stove.dz());

		for (unsigned d = 0; d < 2; ++d) { // depth dim
			if (!(stove.item_flags & (1U<<(2U*w + d)))) continue; // burner not on
			flame_qbd.add_quad_dirs(locs[2*w + d], dsign*radius*plus_x, -dsign*radius*plus_y, WHITE); // use a negative Y to get the proper CW order; flip with dsign for symmetry
		}
	} // for w
	if (flame_qbd.empty()) return;
	select_texture(get_texture_by_name("interiors/gas_burner.png"));
	s.set_color_e(WHITE); // emissive
	enable_blend();
	set_additive_blend_mode();
	glDepthMask(GL_FALSE);
	flame_qbd.draw_and_clear();
	glDepthMask(GL_TRUE);
	set_std_blend_mode();
	disable_blend();
	s.set_color_e(BLACK);
}

void draw_obj_model(obj_model_inst_t const &i, room_object_t const &obj, shader_t &s, vector3d const &xlate, point const &obj_center, bool shadow_only,
	int mirror_dim=3, bool using_custom_tid=0)
{
	bool const emissive_first_mat(!shadow_only && obj.type == TYPE_LAMP      && obj.is_light_on());
	bool const emissive_body_mat (!shadow_only && obj.type == TYPE_WALL_LAMP && obj.is_light_on());
	bool const use_low_z_bias(obj.type == TYPE_CUP && !shadow_only);
	bool const untextured(obj.flags & RO_FLAG_UNTEXTURED);
	bool const upside_down((obj.type == TYPE_RAT || obj.type == TYPE_ROACH || obj.type == TYPE_INSECT) && obj.is_broken());
	if (emissive_first_mat) {s.set_color_e(LAMP_COLOR*0.4);}
	if (use_low_z_bias    ) {s.add_uniform_float("norm_bias_scale", 0.5*DEF_NORM_BIAS_SCALE);} // half the default value
	// Note: lamps are the most common and therefore most expensive models to draw
	int const model_id(obj.get_model_id()); // first 8 bits is model ID, last 8 bits is sub-model ID
	unsigned rot_only_mat_mask(0);
	vector3d dir(i.dir);

	if (dir != obj.get_dir()) { // handle models that have rotating parts; similar to car_draw_state_t::draw_helicopter()
		if      (obj.type == TYPE_CEIL_FAN ) {rot_only_mat_mask =  1;} // only the first material (fan blades) rotate
		else if (obj.type == TYPE_OFF_CHAIR) {rot_only_mat_mask = ~1;} // all but the first material (base) rotates

		if (rot_only_mat_mask > 0) { // draw the rotated part
			building_obj_model_loader.draw_model(s, obj_center, obj, dir, obj.color, xlate, model_id, shadow_only,
				0, nullptr, ~rot_only_mat_mask, untextured, 0, upside_down, emissive_body_mat, 0, 3, using_custom_tid);
			dir = obj.get_dir(); // base model rotation based on object dim/dir orient and not instance rotation vector
		}
	}
	// disable the leg and rubber feet for hanging monitors; materials are {glass screen, plastic body, logo + metal + leg, object, rubber feet}
	if (obj.is_tv_or_monitor() && obj.is_hanging()) {rot_only_mat_mask |= 20;}

	building_obj_model_loader.draw_model(s, obj_center, obj, dir, obj.color, xlate, model_id, shadow_only,
		0, nullptr, rot_only_mat_mask, untextured, 0, upside_down, emissive_body_mat, 0, mirror_dim, using_custom_tid);
	if (!shadow_only && obj.type == TYPE_STOVE) {draw_stove_flames(obj, (camera_pdu.pos - xlate), s);} // draw blue burner flame
	if (use_low_z_bias    ) {s.add_uniform_float("norm_bias_scale", DEF_NORM_BIAS_SCALE);} // restore to the defaults
	if (emissive_first_mat) {s.set_color_e(BLACK);}
}

void brg_batch_draw_t::draw_obj_models(shader_t &s, vector3d const &xlate, bool shadow_only) const {
	for (obj_model_inst_with_obj_t const &i : models_to_draw) {
		point const obj_center(i.obj.get_cube_center());
		try_bind_tile_smap_at_point(obj_center, s);
		draw_obj_model(i, i.obj, s, xlate, obj_center, shadow_only);
	}
	if (!models_to_draw.empty()) {check_mvm_update();}
}

float get_ao_shadow(room_object_t const &c, bool enable_indir) {
	room_object const type(c.type);
	// include types that don't contribute to indir lighting; these always contribute AO shadows
	if (type == TYPE_BAR_STOOL || type == TYPE_SHELVES || type == TYPE_CONV_BELT || type == TYPE_EX_MACHINE) return 0.25; // light shadow
	if (type == TYPE_OFF_CHAIR || type == TYPE_BENCH || type == TYPE_RCHAIR || type == TYPE_CASHREG || type == TYPE_CHEM_TANK || type == TYPE_HOSP_BED) return 0.5; // medium shadow
	if (type == TYPE_WHEELCHAIR || type == TYPE_OP_TABLE || type == TYPE_TROLLEY || type == TYPE_EX_MACHINE) return 0.5; // medium shadow
	if (type == TYPE_PARK_SPACE && c.is_used()) return 0.75; // parked car; dense shadow
	//TYPE_PALLET? or is that too expensive?
	if (enable_indir) return 0.0; // skip objects below because they're already handled by indir lighting

	if (type == TYPE_TABLE) {
		if (c.is_glass_table()) return 0.0; // no shadow
		return ((c.shape == SHAPE_TALL) ? 0.25 : 0.5); // small/medium shadow
	}
	if (type == TYPE_BED) {return ((c.taken_level > 2) ? 0.35 : 0.5);} // reduced AO when the mattress has been taken and light gets through the slats
	if (type == TYPE_BCASE || type == TYPE_DRESSER || type == TYPE_NIGHTSTAND || type == TYPE_COUCH || type == TYPE_CONF_TABLE) return 0.75; // dense shadow
	if (type == TYPE_SINK || type == TYPE_TOILET || type == TYPE_STALL) return 0.25; // light shadow
	if (type == TYPE_CHAIR || type == TYPE_DESK || type == TYPE_RDESK || type == TYPE_POOL_TABLE || type == TYPE_MACHINE || type == TYPE_XFORMER) return 0.5; // med shadow
	return 0.0; // no shadow
}

void draw_and_clear_flares(quad_batch_draw &qbd, shader_t &s, colorRGBA const &color) {
	if (qbd.empty()) return;
	s.set_color_e(RED);
	draw_and_clear_blur_qbd(qbd);
	s.clear_color_e();
}
quad_batch_draw flare_qbd; // resed across/between frames

cube_t union_of_two_cubes(cube_t const &a, cube_t const &b) {
	cube_t bc(a);
	bc.union_with_cube(b);
	return bc;
}

bool check_custom_occluder_cull(cube_t const &c, point const &viewer, vect_cube_t const &occluders, bool occluder_dim);

// Note: non-const because it creates the VBO; inc_small: 0=large only, 1=large+small, 2=large+small+ext detail, 3=large+small+ext detail+int detail, 4=ext only
void building_room_geom_t::draw(brg_batch_draw_t *bbd, shader_t &s, shader_t &amask_shader, building_t const &building, occlusion_checker_noncity_t &oc,
	vector3d const &xlate, unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building, bool mall_visible)
{
	if (empty()) return; // no geom
	unsigned const num_screenshot_tids(get_num_screenshot_tids());
	static int last_frame(0);
	static unsigned num_geom_this_frame(0); // used to limit per-frame geom gen time; doesn't apply to shadow pass, in case shadows are cached
	if (frame_counter < 100 || frame_counter > last_frame) {num_geom_this_frame = 0; last_frame = frame_counter;} // unlimited for the first 100 frames
	point const camera_bs(camera_pdu.pos - xlate);
	float const floor_spacing(building.get_window_vspace()), ground_floor_z1(building.ground_floor_z1);
	bool const draw_ext_only(inc_small == 4), check_occlusion(display_mode & 0x08), is_industrial(building.is_industrial());
	if (draw_ext_only) {inc_small = 0;}
	// don't draw ceiling lights when player is above the building unless there's a light placed on a skylight
	bool const draw_lights(!draw_ext_only && (camera_bs.z < building.bcube.z2() + (building.has_skylight_light ? 20.0*floor_spacing : 0.0)));
	// only industrial, parking garages, backrooms, and attics have detail objects that cast shadows
	bool const draw_detail_objs(inc_small >= 2 && (!shadow_only || is_industrial || building.point_in_attic(camera_bs) || building.is_pos_in_pg_or_backrooms(camera_bs)));
	bool const draw_int_detail_objs(inc_small >= 3 && !shadow_only);
	bool player_in_bldg_normal_pass(player_in_building && !shadow_only && !reflection_pass);
	// update clocks if moved to next second; only applies to the player's building
	bool const update_clocks(player_in_bldg_normal_pass && inc_small >= 2 && have_clock && check_clock_time());
	bool const player_in_doorway(building.point_near_ext_door(camera_bs, get_door_open_dist()));
	bool const player_in_building_or_doorway(player_in_building || player_in_doorway);
	if (bbd != nullptr) {bbd->set_camera_dir_mask(camera_bs, ((camera_bs.z < ground_floor_z1) ? building.get_bcube_inc_extensions() : building.bcube));}
	brg_batch_draw_t *const bbd_in(bbd); // capture bbd for instance drawing before setting to null if player_in_building
	if (player_in_building_or_doorway) {bbd = nullptr;} // use immediate drawing when player is in the building because draw order matters for alpha blending
	bool const update_tunnel_water(animate2 && player_in_building && !shadow_only && (player_in_tunnel || building.interior->point_near_tunnel_entrance(camera_bs)));
	bool enable_indir(0), update_escalators(0);

	if (animate2 && player_in_bldg_normal_pass) { // maybe update escalators
		for (escalator_t const &e : building.interior->escalators) {
			cube_t const bc(e.get_ramp_bcube(1)); // exclude_sides=1
			if (!building.is_rot_cube_visible(bc, xlate)) continue; // VFC
			if (check_occlusion && building.check_obj_occluded(bc, camera_bs, oc, reflection_pass)) continue;
			update_escalators = 1;
			break;
		}
	}
	if (player_in_bldg_normal_pass) { // indir lighting auto update logic
		static bool last_enable_indir(0);
		enable_indir = enable_building_indir_lighting_no_cib();
		if (enable_indir != last_enable_indir) {invalidate_mats_mask |= 0xFF;} // update all geom when material lighting changes due to indir
		last_enable_indir = enable_indir;
	}
	if (has_pictures && num_pic_tids != num_screenshot_tids) {
		invalidate_static_geom(); // user created a new screenshot texture, and this building has pictures - recreate room geom
		num_pic_tids = num_screenshot_tids;
	}
	if (animate2 && (update_clocks || update_escalators || update_tunnel_water || (player_in_bldg_normal_pass && have_conv_belt))) {update_dynamic_draw_data();}
	check_invalid_draw_data();

	// generate vertex data in the shadow pass or if we haven't hit our generation limit unless this is the first frame; must be consistent for static and small geom
	// Note that the distance cutoff for mats_static and mats_small is different, so we generally won't be creating them both
	// unless the player just appeared by this building, or we need to update the geometry; in either case this is higher priority and we want to update both
	if (shadow_only || frame_counter <= 1 || num_geom_this_frame < max(global_building_params.max_room_geom_gen_per_frame, 1U)) {
		if (!mats_static.valid) { // create static materials if needed
			//highres_timer_t timer("Create Static VBOs");
			create_obj_model_insts(building);
			create_static_vbos    (building);
			if (!shadow_only) {++num_geom_this_frame;}
		}
		// as an optimization, we can skip drawing mall objects if not visible and the mall has not yet been drawn
		bool const has_mall(building.has_mall());
		bool const skip_mall_objs(!mall_geom_drawn && has_mall && (shadow_only || !(mall_visible || building.player_can_see_inside_mall(xlate))));

		if (has_mall && !skip_mall_objs && !mall_geom_drawn) { // first draw of mall geom requires clearing existing mats_small and mats_text
			mats_small.clear();
			mats_text .clear();
			mall_geom_drawn = 1; // once we've drawn the mall geom we must always draw it so that shadows are updated and objects can be moved out of the mall
		}
		bool const create_small(inc_small && !mats_small.valid), create_text(draw_int_detail_objs && !mats_text.valid);
		//highres_timer_t timer("Create Small + Text VBOs", (create_small || create_text));

		// Note: shelf rack book text is drawn in the text pass to make it thread safe
		if (create_small && create_text) { // MT case
#pragma omp parallel num_threads(2)
			if (omp_get_thread_num_3dw() == 0) {create_small_static_vbos(building, skip_mall_objs);} else {create_text_vbos(skip_mall_objs);}
		}
		else { // serial case
			if (create_small) {create_small_static_vbos(building, skip_mall_objs);}
			if (create_text ) {create_text_vbos(skip_mall_objs);}
		}
		add_small_static_objs_to_verts(pending_objs, 0, 0, building.get_trim_color(), create_text, building.get_floor_ceil_gap(), building.interior->ind_info.get());
		pending_objs.clear();

		// upload VBO data serially
		if (create_small) {
			mats_small   .create_vbos(building);
			mats_amask   .create_vbos(building);
			mats_alpha_sm.create_vbos(building);
		}
		if (create_text) {mats_text.create_vbos(building);}
		if (!shadow_only) {num_geom_this_frame += (unsigned(create_small) + unsigned(create_text));}

		// Note: not created on the shadow pass unless trim_objs has been created so that we don't miss including it;
		// the trim_objs test is needed to handle parking garage and attic objects, which are also drawn as details
		if (draw_detail_objs && (!shadow_only || !trim_objs.empty()) && !mats_detail.valid) { // create detail materials if needed (mats_detail and mats_ext_detail)
			create_detail_vbos(building);
			if (!shadow_only) {++num_geom_this_frame;}
		}
	}
	if (draw_lights && !mats_lights .valid) {create_lights_vbos (building);} // create lights  materials if needed (no limit)
	if (inc_small   && !mats_dynamic.valid) {create_dynamic_vbos(building, camera_bs, xlate, update_clocks);} // create dynamic materials if needed (no limit)
	if (!mats_doors.valid) {create_door_vbos(building);} // create door materials if needed (no limit)
	if (!shadow_only) {enable_blend();} // needed for rugs and book text
	bool const cube_map_ref(reflection_pass && is_cube_map_reflection);
	int const ref_pass(reflection_pass ? (cube_map_ref ? 2 : 1) : 0); // set ref_pass=1 for cube maps to disable some small objects
	assert(s.is_setup());
	if (!draw_ext_only) {mats_static .draw(bbd, s, shadow_only, ref_pass);}
	if (draw_lights)    {mats_lights .draw(bbd, s, shadow_only, ref_pass);}
	if (inc_small  )    {mats_dynamic.draw(bbd, s, shadow_only, ref_pass);}
	if (draw_detail_objs && inc_small >= 3) {mats_detail.draw(bbd, s, shadow_only, ref_pass);} // now included in the shadow pass
	
	// draw exterior geom; shadows not supported; always use bbd;
	// skip in reflection pass because that control flow doesn't work and is probably not needed (except for L-shaped house?)
	if (!shadow_only && !reflection_pass && player_in_basement < 2 && !cube_map_ref) { // skip for player fully in the basement
		// is there a way to incrementally blend/dither this geometry in? it's not drawn in the correct order for alpha blending to work
		mats_exterior.draw(bbd_in, s, shadow_only, ref_pass, 1); // exterior_geom=1
		if (draw_detail_objs) {mats_ext_detail.draw(bbd_in, s, shadow_only, ref_pass, 1);} // exterior_geom=1
	}
	if (!draw_ext_only) {mats_doors.draw(bbd, s, shadow_only, ref_pass);}
	if (inc_small     ) {mats_small.draw(bbd, s, shadow_only, ref_pass);}

	if (0 && player_in_building && !shadow_only && !reflection_pass && (display_mode & 0x20)) {
		cout << "static: " << mats_static.count_all_verts() << endl;
		cout << "small : " << mats_small.count_all_verts() << endl;
		cout << "sm_sh : " << mats_small.count_all_verts(1) << endl;
		cout << "sm_ref: " << mats_small.count_all_verts(0, 1) << endl;
		cout << "text  : " << mats_text.count_all_verts() << endl;
		cout << "txt_r : " << mats_text.count_all_verts(0, 1) << endl;
		cout << "detail: " << mats_detail.count_all_verts() << endl;
		cout << "det_r : " << mats_detail.count_all_verts(0, 1) << endl;
		cout << "dynam : " << mats_dynamic.count_all_verts() << endl;
		cout << "lights: " << mats_lights.count_all_verts() << endl;
		cout << "amask : " << mats_amask.count_all_verts() << endl;
		cout << "alpha : " << mats_alpha.count_all_verts() << endl;
		cout << "a_sm  : " << mats_alpha_sm.count_all_verts() << endl;
		cout << "doors : " << mats_doors.count_all_verts() << endl;
		cout << "ext   : " << mats_exterior.count_all_verts() << endl;
		cout << "ext_dt: " << mats_ext_detail.count_all_verts() << endl;
		cout << endl;
	}
	if (!mats_amask.empty() && !draw_ext_only) { // draw plant leaves, spider webs, etc. using alpha masks in the detail pass
		if (shadow_only) {
			if (!amask_shader.is_setup()) {amask_shader.begin_simple_textured_shader(0.9);} // need to use texture with alpha test
			else {amask_shader.make_current();} // min_alpha should be left at 0.9 from the previous call
			mats_amask.draw(nullptr, amask_shader, 2, 0); // shadow pass with alpha mask; no brg_batch_draw
			s.make_current(); // switch back to the normal shader
		}
		else if (reflection_pass) {
			mats_amask.draw(nullptr, s, 0, 1); // no brg_batch_draw
		}
		// this is expensive: only enable for the main draw pass and skip for buildings the player isn't in, except for industrial building metal grates
		else if (inc_small >= 2 && (player_in_building || is_industrial)) {
			// without the special shader these won't look correct when drawn through windows
			// used for both plant/tree leaves and spider webs;
			// plants are above ground and in malls with high min_alpha; metal stairs are above ground and in basements with high min_alpha;
			// spider webs are in ext basement with low min_alpha
			float const min_alpha((player_in_basement >= 3 && !building.has_mall()) ? 0.1 : 0.9);

			if (!amask_shader.is_setup()) {
				setup_building_draw_shader(amask_shader, min_alpha, 1, 1, 0); // enable_indir=1, force_tsl=1, use_texgen=0, water_damage=0.0
			}
			else {
				amask_shader.make_current();
				amask_shader.add_uniform_float("min_alpha", min_alpha); // set min_alpha in case it changed
			}
			mats_amask.draw(nullptr, amask_shader, 0, 0); // no brg_batch_draw
			s.make_current(); // switch back to the normal shader
		}
	}
	if (draw_int_detail_objs) {mats_text.draw(bbd, s, shadow_only, ref_pass);} // text must be drawn last; drawn as interior detail objects
	if (!shadow_only) {disable_blend();}
	indexed_vao_manager_with_shadow_t::post_render();
	if (draw_ext_only) return; // done
	bool const disable_cull_face(0); // better but slower?
	if (disable_cull_face) {glDisable(GL_CULL_FACE);}
	point const building_center(building.bcube.get_cube_center());
	oc.set_exclude_bix(building_ix);
	water_sound_manager_t water_sound_manager(camera_bs);
	rgeom_mat_t monitor_screens_mat, onscreen_text_mat=rgeom_mat_t(tid_nm_pair_t(FONT_TEXTURE_ID));
	string onscreen_text;
	bool const is_rotated(building.is_rotated()), is_player_building(&building == player_building), has_pri_hall(building.has_pri_hall());
	bool check_clip_cube(shadow_only && !is_rotated && !smap_light_clip_cube.is_all_zeros()); // check clip cube for shadow pass; not implemented for rotated buildings
	check_clip_cube |= (reflection_pass && !reflection_clip_cube.is_all_zeros()); // handle reflection clip cube as well
	bool const skip_interior_objs(!player_in_building_or_doorway && !shadow_only), has_windows(building.has_windows());
	bool const player_in_industrial(building.point_in_industrial(camera_bs));
	bool const player_in_this_basement(player_in_building && player_in_basement >= 2), player_above_this_basement(player_in_building && player_in_basement == 0);
	float const one_floor_above(camera_bs.z + floor_spacing);
	float two_floors_below(camera_bs.z - 2.0*floor_spacing);
	if (player_in_industrial) {min_eq(two_floors_below, ground_floor_z1);} // industrial lights reach more than 2 floors
	cube_t const clip_cube_bs(shadow_only ? (smap_light_clip_cube - xlate) : reflection_clip_cube);
	// skip for rotated buildings and planar reflection pass, since reflected pos may be in a different room; should we use actual_player_pos for shadow_only mode?
	int const camera_room((is_rotated || (reflection_pass && !cube_map_ref)) ? -1 : building.get_room_containing_camera(camera_bs));
	int camera_part(-1), cull_room_ix(-1);
	unsigned last_room_ix(building.interior->rooms.size()), last_culled_room_ix(last_room_ix), last_floor_ix(0); // start at an invalid value
	bool camera_in_closed_room(0), last_room_closed(0), obj_drawn(0), no_cull_room(0);
	auto model_to_cull(obj_model_insts.end()), model_to_not_cull(obj_model_insts.end());
	vect_cube_t const *extra_occluders(nullptr);
	bool extra_occluders_dim(0);

	if (camera_room >= 0) {
		// check for stairs in case an object is visible through an open door on a floor above or below
		room_t const &room(building.get_room(camera_room));
		unsigned const camera_floor(room.get_floor_containing_zval(max(camera_bs.z, room.z1()), floor_spacing)); // clamp zval to room range
		camera_in_closed_room = (!room.has_stairs_on_floor(camera_floor) && building.all_room_int_doors_closed(camera_room, camera_bs.z));
		camera_part           = room.part_id;

		if (check_occlusion && room.is_store()) { // add mall store occluders for the room the player is in
			for (store_info_t const &s : building.interior->mall_info->stores) {
				if ((int)s.room_id != camera_room) continue;
				extra_occluders     = &s.occluders;
				extra_occluders_dim = !s.dim; // currently, all occluder walls are perpendicular to the store dim
				break;
			}
			no_cull_room = 1; // don't need to cull objects in the same store as the player
		}
	}
	// draw object models; ~0.4ms average for malls
	for (auto i = obj_model_insts.begin(); i != obj_model_insts.end(); ++i) {
		if (i == model_to_cull) continue;
		if (skip_interior_objs && i->int_vis_only) continue; // interior, not visible
		room_object_t &obj(get_room_object_by_index(i->obj_id));
		if ((int)obj.room_id == cull_room_ix)                 continue; // cull all objects in this room
		if (check_clip_cube && !clip_cube_bs.intersects(obj)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
		// optimization: only draw models in the same room as the mirror for malls; applies to furniture stores, clothing stores, and bathrooms
		if (reflection_pass && player_in_mall && cur_room_mirror.room_id > 0 && obj.room_id != cur_room_mirror.room_id) continue;
		bool const in_camera_room((int)obj.room_id == camera_room);
		room_object const type(obj.type);

		if (shadow_only) {
			if (type == TYPE_CEIL_FAN) continue; // not shadow casting; would shadow its own light
			if (obj.is_exterior())     continue; // outdoors; no indoor shadow
			if (player_in_mall && type == TYPE_HANGER) continue; // too small to cast a shadow for lights high above
			if (player_in_mall && type == TYPE_CLOTHES && !in_camera_room) continue; // optimization for expensive clothing
			if (type == TYPE_KEY || type == TYPE_SILVER || type == TYPE_FOLD_SHIRT)                 continue; // small
			if (obj.z1() > camera_bs.z || (obj.z2() < two_floors_below && type != TYPE_BLDG_FOUNT)) continue; // above or more than 2 floors below light, except mall fountains
		}
		else if (cube_map_ref) {
			if (type == TYPE_KEY || type == TYPE_PADLOCK || type == TYPE_SILVER || type == TYPE_CUP || type == TYPE_HANGER ||
				type == TYPE_BANANA || type == TYPE_BAN_PEEL || type == TYPE_APPLE || type == TYPE_TOY_MODEL) continue; // small
			if (!in_camera_room) { // cull expensive objects in different room from the player
				if (type == TYPE_CONF_PHONE || type == TYPE_VIS_PHONE || type == TYPE_SHOE || type == TYPE_SHOEBOX || type == TYPE_CLOTHES ||
					type == TYPE_FOLD_SHIRT || type == TYPE_FIRE_EXT || type == TYPE_TOASTER || type == TYPE_LAMP) continue;
			}
		}
		else if ((type == TYPE_SILVER /*|| type == TYPE_FOLD_SHIRT*/) && camera_bs.z < obj.z1()) continue; // not visible from below (except on glass table?)
		point obj_center(obj.get_cube_center());

		if (!shadow_only && !building.is_house && !has_windows && !building.point_in_mall(obj_center)) { // windowless building
			if (obj.z1() > one_floor_above || obj.z2() < two_floors_below) continue; // more than one floor of difference
		}
		if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
		
		// distance culling
		if (!shadow_only) {
			float cull_dist(32.0*(obj.dx() + obj.dy() + obj.dz()));

			if (player_in_building_or_doorway) { // increase culling distance for player building
				if (has_pri_hall && building.pri_hall.contains_pt(obj_center)) {cull_dist *= 10.0;} // primary hallway objects visible down a long hallway
				else if (type == TYPE_FIRE_EXT || type == TYPE_PADLOCK)        {cull_dist *= 3.0 ;} // padlocks are small but can be seen from far away
				else if (building.check_pt_in_retail_room(obj_center))         {cull_dist *= 2.5 ;} // increased culling distance for retail areas
				else if (building.point_in_mall          (obj_center))         {cull_dist *= ((type == TYPE_HANGER) ? 1.2 : 2.0);} // increase dist for malls; less for hangers
			}
			else { // player not in this building
				if (camera_in_building) {cull_dist *= 0.7;} // decrease culling distance for player in a different building (optimization)
				if (obj.is_interior())  {cull_dist *= 0.7;} // interior room objects are unlikely to be seen by the player
			}
			if (!dist_less_than(camera_bs, obj_center, cull_dist)) continue; // too far
		}
		bool cull(0);
		cull |= (player_in_this_basement    && obj_center.z > ground_floor_z1); // player in basement, obj not
		cull |= (player_above_this_basement && obj_center.z < ground_floor_z1); // obj in basement, player not

		if (cull) { // basement separation; check for primary stairs visibility (for reception desk chair, etc.)
			vect_stairwell_t const &sw(building.interior->stairwells);
			if (sw.empty() || sw.front().is_u_shape() || !sw.front().line_intersects(camera_bs, obj_center)) continue;
		}
		//if (shadow_only && type == TYPE_SHOEBOX) {} // draw as a cube? maybe TYPE_TOASTER and TYPE_FOLD_SHIRT as well?
		room_t const &room(building.get_room(obj.room_id));

		if (player_in_building && obj.room_id != last_culled_room_ix && !obj.is_exterior()) { // new room; apply room-based VFC + occlusion culling
			last_culled_room_ix = obj.room_id;

			if (!in_camera_room) { // camera not in this room
				cube_t c(room);
				c.expand_in_z (-building.get_fc_thickness  ());
				c.expand_by_xy(-building.get_wall_thickness());
				if (!(is_rotated ? building.is_rot_cube_visible(c, xlate) : camera_pdu.cube_visible(c + xlate)) ||
					(check_occlusion && building.check_obj_occluded(c, camera_bs, oc, reflection_pass))) {cull_room_ix = obj.room_id; continue;} // cull entire room
			}
		}
		// mall shoe shadows are expensive; only enable if player is in the room
		if (shadow_only && (type == TYPE_SHOE || type == TYPE_SHOEBOX) && !room.contains_pt(actual_player_pos) && building.point_in_mall(obj_center)) continue;
		if (!(is_rotated ? building.is_rot_cube_visible(obj, xlate) : camera_pdu.cube_visible(obj + xlate))) continue; // VFC
		// check for parking garage vs. mall/backrooms separation
		if (!shadow_only && player_in_this_basement && room.is_ext_basement() != (player_in_basement >= 3) && !building.is_cube_visible_through_extb_door(camera_bs, obj)) continue;
		//highres_timer_t timer(string("Draw ") + get_room_obj_type(obj).name + (shadow_only ? " Shadow" : (reflection_pass ? " Ref" : "")), (display_mode & 0x20));

		if (check_occlusion && i != model_to_not_cull) { // occlusion culling
			if (!(no_cull_room && in_camera_room)) { // skip culling of objects in the same room as the player/light, except for extra_occluders
				bool no_check_occlude(0);

				if (type == TYPE_HANGER && obj.is_hanging() && i+1 != obj_model_insts.end()) { // hanger
					room_object_t &obj2(get_room_object_by_index((i+1)->obj_id));

					if (obj2.type == TYPE_CLOTHES || obj2.type == TYPE_JUMPSUIT) { // cull hanger and clothing/jumpsuit together
						if (building.check_obj_occluded(union_of_two_cubes(obj, obj2), camera_bs, oc, reflection_pass)) {model_to_cull = i+1; continue;}
						model_to_not_cull = i+1;
					}
				}
				if (type == TYPE_SHOE && (obj.flags & RO_FLAG_ADJ_LO) && i+1 != obj_model_insts.end()) { // first shoe in pair
					room_object_t &obj2(get_room_object_by_index((i+1)->obj_id));

					if (obj2.type == TYPE_SHOE && (obj.flags & RO_FLAG_ADJ_HI)) { // second shoe in pair
						if (building.check_obj_occluded(union_of_two_cubes(obj, obj2), camera_bs, oc, reflection_pass)) {model_to_cull = i+1; continue;}
						model_to_not_cull = i+1;
						no_check_occlude  = 1;
					}
				}
				if (!no_check_occlude && building.check_obj_occluded(obj, camera_bs, oc, reflection_pass)) continue;
			}
		}
		if (extra_occluders && check_custom_occluder_cull(obj, camera_bs, *extra_occluders, extra_occluders_dim)) continue;

		if (camera_room >= 0) {
			unsigned const floor_ix(room.get_floor_containing_zval(obj.zc(), floor_spacing));

			if (obj.room_id != last_room_ix || floor_ix != last_floor_ix) { // new room or new floor
				last_room_closed = building.all_room_int_doors_closed(obj.room_id, obj.zc());
				last_room_ix     = obj.room_id;
				last_floor_ix    = floor_ix;
			}
			// if either the camera or the object are in different rooms with closed doors,
			// on the same floor (not separated by stairs) of the same part (not visible across windows), then the object isn't visible
			if ((last_room_closed || camera_in_closed_room) && !in_camera_room && (room.part_id == camera_part || !has_windows)) continue;
		}
		apply_room_obj_rotate(obj, *i, objs); // Note: may modify obj by clearing flags and inst by updating dir
		
		if (bbd_in && !shadow_only && !is_rotated && (type == TYPE_WALL_LAMP || obj.is_exterior())) { // draw exterior objects later; not for rotated buildings
			// wall lamp has transparent glass and must be drawn last; fire escape and wall lamp use outdoor lighting
			bbd_in->models_to_draw.emplace_back(*i, obj);
		}
		else { // draw now
			//draw_simple_cube(obj); // TESTING
			int mirror_dim(3); // 3=none
			bool const using_custom_tid(building.bind_custom_clothing_texure(obj));
			if (type == TYPE_SHOE && (obj.flags & RO_FLAG_ADJ_TOP)) {mirror_dim = 1;} // shoes may be mirrored in !obj.dim (Y in model space)
			bool const is_metal(type == TYPE_SILVER || type == TYPE_KEY);
			if (is_metal) {s.add_uniform_float("metalness", 1.0);}
			draw_obj_model(*i, obj, s, xlate, obj_center, shadow_only, mirror_dim, using_custom_tid);
			if (is_metal) {s.add_uniform_float("metalness", 0.0);}
			obj_drawn = 1;
		}
		// check for security camera monitor if player is in this building; must be on and active
		if (player_in_building && type == TYPE_MONITOR && obj.is_tv_monitor_on() && obj.is_active()) {
			onscreen_text.clear();
			setup_monitor_screen_draw(obj, monitor_screens_mat, onscreen_text);
			add_tv_or_monitor_screen (obj, monitor_screens_mat, onscreen_text, &onscreen_text_mat);
			s.set_color_e(WHITE); // emissive
			tid_nm_pair_dstate_t screen_state(s, 1), text_state(s, 0); // no_set_texture=1/0
			monitor_screens_mat.upload_draw_and_clear(screen_state);

			if (!onscreen_text_mat.empty()) {
				enable_blend();
				onscreen_text_mat.upload_draw_and_clear(text_state);
				disable_blend();
			}
			s.set_color_e(BLACK);
		}
		if (player_in_building && !shadow_only && type == TYPE_SINK) { // sink
			if (in_camera_room) {water_sound_manager.register_running_water(obj, building);}
			water_draw.add_water_for_sink(obj);
		}
	} // for i
	if (!skip_interior_objs && !cube_map_ref && !door_handles.empty()) { // optimization: skip door handles for player outside building
		colorRGBA const handle_color(building.get_door_handle_color());
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_DOOR_HANDLE)); // L, W, H
		vector3d const exp_val((0.5/sz.z)*sz);
		s.add_uniform_float("metalness", 1.0); // door handles are metal

		for (door_handle_t const &h : door_handles) {
			cube_t bc(h.center);
			bc.expand_by(h.height*exp_val);
			if (check_clip_cube && !clip_cube_bs.intersects(bc)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first

			if (shadow_only) {
				if (bc.z1() > camera_bs.z || bc.z2() < two_floors_below) continue; // above or more than two floors below the light
			}
			else if (!building.is_house && !has_windows) { // windowless building
				if (bc.z1() > one_floor_above || bc.z2() < two_floors_below) continue; // more than one floor of difference
			}
			point obj_center(h.center);
			if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
			if (h.closed && (camera_bs[h.dim] < obj_center[h.dim]) == h.hdir) continue; // opposite side of closed door
			if (player_in_this_basement    && obj_center.z > ground_floor_z1) continue; // player in basement, obj not
			if (player_above_this_basement && obj_center.z < ground_floor_z1) continue; // obj in basement, player not
			if (!(is_rotated ? building.is_rot_cube_visible(bc, xlate) : camera_pdu.cube_visible(bc + xlate))) continue; // VFC

			if (!h.closed) { // check if blocked by open door
				float const dp(dot_product_ptv(vector3d(-h.dir.y, h.dir.x, 0.0), camera_bs, obj_center)*(h.mirror ? 1.0 : -1.0));
				if (dp < -0.5*p2p_dist(camera_bs, obj_center)) continue; // rotate h.dir by 90 degrees
			}
			if (check_occlusion && building.check_obj_occluded(bc, camera_bs, oc, reflection_pass)) continue;
			// can door handles (and models in general) be drawn with instancing? requires custom shader setup and many control flow changes
			building_obj_model_loader.draw_model(s, obj_center, bc, h.dir, handle_color, xlate, OBJ_MODEL_DOOR_HANDLE, shadow_only, 0, nullptr, 0, 0, 0, h.mirror);
			obj_drawn = 1;
		} // for h
		s.add_uniform_float("metalness", 0.0); // reset
	}
	if (player_in_bldg_normal_pass) { // only drawn for the player building
		// key/silverware/folded shirt: not drawn in the shadow or reflection passes; no emissive or rotated objects; no animations
		for (room_object_t const &c : model_objs) {
			float const cull_dist(((c.type == TYPE_KEY) ? 200.0 : 100.0)*c.max_len());
			point obj_center(c.get_cube_center());
			if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
			if (!shadow_only && !dist_less_than(camera_bs, obj_center, cull_dist))                           continue; // too far away
			if (player_in_this_basement    && obj_center.z > ground_floor_z1)                                continue; // player in basement, obj not
			if (player_above_this_basement && obj_center.z < ground_floor_z1)                                continue; // obj in basement, player not
			if (!(is_rotated ? building.is_rot_cube_visible(c, xlate) : camera_pdu.cube_visible(c + xlate))) continue; // VFC
			if (check_occlusion && building.check_obj_occluded(c, camera_bs, oc, reflection_pass))           continue;
			vector3d dir(c.get_dir());
			if (is_rotated) {building.do_xy_rotate_normal(dir);}
			int swap_xy_mode(0);

			if (c.type == TYPE_KEY && c.is_hanging()) { // make hanging key vertical
				if (c.dim) {swap_xy_mode = 1;}
				dir = plus_z;
			}
			bool const is_metal(c.type == TYPE_SILVER || c.type == TYPE_KEY);
			if (is_metal) {s.add_uniform_float("metalness", 1.0);}
			building_obj_model_loader.draw_model(s, obj_center, c, dir, c.color, xlate, c.get_model_id(), shadow_only, 0, nullptr, 0, 0, 0, 0, 0, 0, 3, 0, swap_xy_mode);
			if (is_metal) {s.add_uniform_float("metalness", 0.0);}
			obj_drawn = 1;
		} // for c
	}
	if (player_in_building_or_doorway && !cube_map_ref) {
		// draw animals; skip animal shadows for industrial areas, retail rooms, and malls/stores as an optimization (must agree with lighting code)
		if (shadow_only && (player_in_industrial || (building.has_retail() && building.get_retail_part().contains_pt(camera_bs)) || building.point_in_mall(camera_bs))) {}
		else {draw_animals(s, building, oc, xlate, camera_bs, shadow_only, reflection_pass, check_clip_cube);}
	}
	if (ADD_WORKER_HARDHATS && player_in_building_or_doorway && is_industrial && !building.interior->people.empty()) {
		// draw hard hats on industrial building people; this doesn't look correct because the hat clips through the models when animated
		static rgeom_mat_t mat;
		mat.tex.set_specular(0.6, 70.0);

		for (person_t const &p : building.interior->people) {
			if (p.is_zombie) continue; // doesn't look correct on zombies
			if (p.pos.z < ground_floor_z1) continue; // hard hat not needed in parking garage, mall, backrooms, or extended basement
			bool const dim(0), dir(1);
			float const hh_depth(0.42*p.radius), hh_width(0.8*hh_depth), hh_height(0.56*hh_depth), top_of_heat(p.get_z2());
			point const pos(p.pos + 0.15*hh_depth*p.dir); // center of hardhat is slighly forward
			cube_t hhat;
			set_cube_zvals(hhat, top_of_heat-0.75*hh_height, top_of_heat+0.25*hh_height);
			set_wall_width(hhat, pos[ dim], 0.5*hh_depth,  dim); // set thickness
			set_wall_width(hhat, pos[!dim], 0.5*hh_width, !dim); // set width
			room_object_t const hh(hhat, TYPE_HARDHAT, 0, dim, dir, 0, 1.0, SHAPE_CUBE, hardhat_colors[p.ssn%NUM_HARDHAT_COLORS]);
			unsigned const verts_start(mat.itri_verts.size());
			add_hard_hat_to_material(hh, mat);
			rotate_verts(mat.itri_verts, plus_z, -atan2(p.dir.y, p.dir.x), pos, verts_start);
		} // for p
		tid_nm_pair_dstate_t state(s);
		mat.upload_draw_and_clear(state);
	}
	if (disable_cull_face) {glEnable(GL_CULL_FACE);} // re-enable face culling
	if (obj_drawn) {check_mvm_update();} // needed after popping model transform matrix

	// draw water for sinks that are turned on, lava lamps, fish in fishtanks, and AO shadows; these aren't visible when the player is outside looking in through a window
	if (player_in_building_or_doorway && !shadow_only && !cube_map_ref) {
		// only update and draw fish in the player building (since two extended basement bcubes can overlap); skip in reflection pass
		bool const draw_fish(!reflection_pass && have_fish_model() && (player_in_doorway || (is_player_building && building.point_in_building_or_basement_bcube(camera_bs))));
		float const ao_z_off(1.1*building.get_flooring_thick()); // slightly above rugs and flooring
		float const ao_zmin(camera_bs.z - 2.0*floor_spacing), ao_max_ext(0.5*floor_spacing);
		static quad_batch_draw ao_qbd;
		bool const inc_pools_and_fb(draw_fish && building.begin_fish_draw());
		auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators
		point pts[4];

		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->type == TYPE_KSINK || i->type == TYPE_BRSINK || i->type == TYPE_VANITY) { // TYPE_SINK is handled above
				if (i->room_id == camera_room) {water_sound_manager.register_running_water(*i, building);}
				if (!i->is_active()) continue; // not turned on
				if (!(is_rotated ? building.is_rot_cube_visible(*i, xlate) : camera_pdu.cube_visible(*i + xlate))) continue; // VFC
				water_draw.add_water_for_sink(*i);
			}
			else if (i->type == TYPE_LAVALAMP) {
				if (!(is_rotated ? building.is_rot_cube_visible(*i, xlate) : camera_pdu.cube_visible(*i + xlate))) continue; // VFC
				lava_lamp_draw.add_lava_lamp(*i, camera_bs, building);
			}
			else if (i->type == TYPE_FISHTANK) {
				if (draw_fish && i->item_flags == TYPE_FISH) { // fishtank with fish
					bool visible(is_rotated ? building.is_rot_cube_visible(*i, xlate) : camera_pdu.cube_visible(*i + xlate)); // VFC
					if (visible && check_occlusion && building.check_obj_occluded(*i, camera_bs, oc, reflection_pass)) {visible = 0;}
					register_fishtank(*i, visible);
				}
			}
			else if (i->type == TYPE_WARN_LIGHT) {
				if (i->is_powered() && is_flashing_light_on()) {
					float const radius(3.0*i->get_radius());
					point const light_center(get_warning_light_src_pos(*i));
					flare_qbd.add_billboard(light_center, camera_bs, plus_x, RED, radius, radius);
				}
			}
			else if (!reflection_pass && i->is_shower()) { // Note: only closest shower plays water sound
				if (i->room_id == camera_room) {water_sound_manager.register_running_water(*i, building);}
				water_draw.add_water_for_shower(*i, (i - objs.begin()), particle_manager);
			}
			if (!reflection_pass && i->z1() < camera_bs.z && i->z1() > ao_zmin - max(0.0f, (i->dz() - floor_spacing))) {
				// camera not below or too far above this object; handle tall objects
				float const ao_shadow(get_ao_shadow(*i, enable_indir));

				if (ao_shadow > 0.0) { // add AO shadow quad on the floor below the object
					if (!is_rotated && !camera_pdu.cube_visible(*i + xlate)) continue; // VFC - may not help much
					float rscale(0.5 + 0.5*(1.0 - ao_shadow)); // 0.5 will be the size of the object; dense shadow is sharper/smaller radius
					if (i->type == TYPE_CASHREG || i->type == TYPE_PARK_SPACE) {rscale *= 0.75;} // bcube is larger than it should be for cash registers and parked cars
					float const dx(i->dx()), dy(i->dy()), rx(min(rscale*dx, 0.5f*dx+ao_max_ext)), ry(min(rscale*dy, 0.5f*dy+ao_max_ext));

					set_z_plane_rect_pts(point(i->xc(), i->yc(), (i->z1() + ao_z_off)), rx, ry, pts);

					if (is_rotated) {
						for (unsigned n = 0; n < 4; ++n) {building.do_xy_rotate(building_center, pts[n]);}
					}
					ao_qbd.add_quad_pts(pts, colorRGBA(0, 0, 0, ao_shadow), plus_z);
				}
			}
		} // for i
		for (person_t const &p : building.interior->people) {
			if (p.is_on_stairs) continue; // AO shadows may look wrong on stairs
			cube_t const pbc(p.get_bcube());
			if (pbc.z1() > camera_bs.z || pbc.z1() < ao_zmin) continue; // camera below or too far above, skip
			if (!is_rotated && !camera_pdu.cube_visible(pbc + xlate)) continue; // VFC - may not help much
			set_z_plane_rect_pts(point(p.pos.x, p.pos.y, (pbc.z1() + ao_z_off)), 0.4*p.radius, 0.4*p.radius, pts);
			ao_qbd.add_quad_pts(pts, colorRGBA(0, 0, 0, 0.4), plus_z);
		} // for p
		// Note: animals are generally too small to have AO shadows
		lava_lamp_draw.draw_and_clear(s);
		if (!reflection_pass) {lava_lamp_draw.next_frame();}
		if (draw_fish) {end_fish_draw(s, inc_pools_and_fb);}
		draw_and_clear_blur_qbd(ao_qbd);
		if (!building.is_factory()) {draw_and_clear_flares(flare_qbd, s, RED);} // factory flares are drawn later
	}
	if (!shadow_only && !reflection_pass) {water_sound_manager.finalize();}
	water_draw.draw_and_clear(s);

	if (player_in_building && !shadow_only && player_held_object.is_valid() && &building == player_building) {
		// draw the item the player is holding; actual_player_pos should be the correct position for reflections
		point const obj_pos((reflection_pass ? actual_player_pos : camera_bs) + CAMERA_RADIUS*cview_dir - vector3d(0.0, 0.0, 0.5*CAMERA_RADIUS));
		player_held_object.translate(obj_pos - player_held_object.get_cube_center());
		//unsigned room_id(building.get_room_containing_pt(obj_pos));
		//if (room_id >= 0) {player_held_object.room_id = room_id;} // is this necessary?
		if (is_ball_type(player_held_object.type)) {draw_ball_in_building(player_held_object, s);} // the only supported dynamic object type
		else if (player_held_object.can_use()) {draw_interactive_player_obj(player_held_object, s, xlate);}
		else {assert(0);}
	}
	// alpha blended, should be drawn near last
	decal_manager.draw_building_interior_decals(s, player_in_building_or_doorway, shadow_only); // draw decals in this building
	
	if (player_in_building && !shadow_only) { // ideally should be drawn after all buildings, but the shaders won't be setup correctly
		if (!building.is_factory()) {particle_manager.draw(s, xlate);} // factory smoke is drawn later
		fire_manager.draw(s, xlate);
	}
	if (!shadow_only && (!mats_alpha.empty() || (inc_small && !mats_alpha_sm.empty()))) { // draw last; not shadow casters; for shower glass, etc.
		enable_blend();
		glDepthMask(GL_FALSE); // disable depth writing
		s.add_uniform_float("refract_ix", 1.6); // refractive glass windows, tables, etc.
		mats_alpha.draw(bbd, s, shadow_only, reflection_pass);
		s.add_uniform_float("refract_ix", 1.0); // reset
		if (inc_small && !cube_map_ref) {mats_alpha_sm.draw(bbd, s, shadow_only, reflection_pass);} // bottles, etc.
		glDepthMask(GL_TRUE);
		disable_blend();
		indexed_vao_manager_with_shadow_t::post_render();
	}
}

void building_t::subtract_stairs_and_elevators_from_cube(cube_t const &c, vect_cube_t &cube_parts, bool inc_stairs, bool inc_elevators) const {
	cube_parts.clear();
	cube_parts.push_back(c);
	if (!interior) return; // error?

	if (inc_elevators) {
		for (elevator_t const &e : interior->elevators) { // clip out holes for elevators
			if (e.intersects(c)) {subtract_cube_from_cubes(e, cube_parts);}
		}
	}
	if (inc_stairs) {
		for (stairwell_t const &s : interior->stairwells) { // clip out holes for stairs
			if (s.intersects(c)) {subtract_cube_from_cubes(s, cube_parts);}
		}
	}
}
bool building_t::glass_floor_visible(vector3d const &xlate, bool from_outside_building) const {
	if (!has_glass_floor())      return 0;
	if (player_in_elevator >= 2) return 0; // not visible from within an elevator with the doors closed
	if (!from_outside_building && !get_retail_room().contains_pt(get_camera_pos() - xlate)) return 0; // wrong room (should always have U-shaped stairs that block visibility)
	return is_rot_cube_visible(get_bcubes_union(interior->room_geom->glass_floors), xlate); // VFC
}
void building_t::draw_glass_surfaces(vector3d const &xlate) const {
	// currently there are only glass floors, but this could be used for drawing showers and fishtanks as well
	if (!has_room_geom()) return;
	vect_cube_t const &glass_floors(interior->room_geom->glass_floors);
	if (glass_floors.empty()) return;
	point const camera_bs(get_camera_pos() - xlate);
	bool const player_is_above(camera_bs.z > glass_floors.front().zc());
	rgeom_mat_t &mat(interior->room_geom->mats_glass[player_is_above]);

	if (mat.empty()) { // create geometry
		float const wall_thickness(get_wall_thickness());

		for (cube_t const &c : glass_floors) {
			vect_cube_t floor_parts;
			subtract_stairs_and_elevators_from_cube(c, floor_parts, 1, 0); // inc_stairs=1, inc_elevators=0 (since we check for player in elevator)
			bool const was_split(floor_parts.size() > 1);
			interior->room_geom->glass_floor_split |= was_split;
			colorRGBA color(GLASS_COLOR, global_building_params.glass_floor_alpha);
			if (was_split) {color.A *= 2.0;} // if split, back face culling is enabled, bottom side is not drawn, so double the alpha value
			min_eq(color.A, 1.0f); // clamp, in case it was specified > 0.5 in the config file

			for (cube_t const &f : floor_parts) {
				// skip faces along building exterior walls; assumes glass floor is in a retail room that spans the entire building bcube;
				// also skip interior faces that are not part of the original floor cube's edge
				unsigned skip_faces(0);

				for (unsigned dim = 0; dim < 2; ++dim) {
					for (unsigned dir = 0; dir < 2; ++dir) {
						if (f.d[dim][dir] != c.d[dim][dir] || fabs(f.d[dim][dir] - bcube.d[dim][dir]) < wall_thickness) {skip_faces |= ~get_face_mask(dim, dir);}
					}
				}
				mat.add_cube_to_verts_untextured(f, color, skip_faces);
			} // for f
		} // for c
		if (!player_is_above) {reverse(mat.quad_verts.begin(), mat.quad_verts.end());} // reverse so that top surface is drawn before bottom surface for correct alpha blending
		mat.create_vbo_inner();
	}
	calc_cur_ambient_diffuse();
	enable_blend();
	// cull back faces if floor was split due to stairs or an elevator, since the bottom won't alpha blend properly
	if (interior->room_geom->glass_floor_split) {glEnable(GL_CULL_FACE);}
	colorRGBA const indoor_light_color(get_retail_light_color());
	// draw mirror reflection if player is on the top surface of the glass floor
	bool const enable_reflection(ENABLE_GLASS_FLOOR_REF && room_mirror_ref_tid > 0 && point_over_glass_floor(camera_bs, 1)); // inc_escalator=1
	shader_t s;
	if (enable_reflection) {s.set_prefix("#define ENABLE_REFLECTION", 1);} // FS
	s.set_vert_shader("glass_surface");
	s.set_frag_shader("fresnel.part*+glass_surface");
	s.begin_shader();
	s.add_uniform_color("light_color", colorRGB(indoor_light_color*0.5 + cur_diffuse*0.2 + cur_ambient*0.5));
	
	if (enable_reflection) {
		s.add_uniform_int("reflection_tex", 0);
		bind_2d_texture(room_mirror_ref_tid);
	}
	mat.vao_setup (0); // shadow_only=0
	mat.draw_inner(0); // shadow_only=0
	indexed_vao_manager_with_shadow_t::post_render();
	s.end_shader();
	if (interior->room_geom->glass_floor_split) {glDisable(GL_CULL_FACE);}
	disable_blend();
}

void building_t::draw_factory_alpha(vector3d const &xlate) const { // smoke and light flares
	if (!has_room_geom() || !is_factory()) return;
	shader_t s;
	s.begin_simple_textured_shader();
	interior->room_geom->particle_manager.draw(s, xlate);
	draw_and_clear_flares(flare_qbd, s, RED);
}

void building_t::draw_z_prepass(point const &camera_bs) const {
	if (!has_room_geom()) return;
	if (camera_bs.z < ground_floor_z1) return; // underground, no special occluders
	bool check_shelfracks(0);
	if      (has_retail  () && get_retail_part().contains_pt(camera_bs)) {check_shelfracks = 1;} // player in retail
	else if (is_warehouse() && point_in_industrial          (camera_bs)) {check_shelfracks = 1;} // player in the warehouse
	if (!check_shelfracks) return;
	quad_batch_draw qbd;
	point const rot_center(bcube.get_cube_center());

	for (cube_t const &c : interior->room_geom->shelf_rack_occluders[0]) { // use shelfrack backs
		bool const dim(c.dy() < c.dx()); // wall dim
		float const center(c.get_center_dim(dim));
		point pts[4] = {point(c.x1(), c.y1(), c.z1()), point(c.x2(), c.y2(), c.z1()), point(c.x2(), c.y2(), c.z2()), point(c.x1(), c.y1(), c.z2())};
		
		for (unsigned n = 0; n < 4; ++n) {
			pts[n][dim] = center; // shrink to zero area in dim
			if (is_rotated()) {do_xy_rotate(rot_center, pts[n]);}
		}
		qbd.add_quad_pts(pts, WHITE);
	} // for c
	qbd.draw();
}

void draw_billboards(quad_batch_draw &qbd, int tid, bool no_depth_write=1, bool do_blend=1) {
	if (qbd.empty()) return;
	if (no_depth_write) {glDepthMask(GL_FALSE);} // disable depth write
	if (do_blend) {enable_blend();}
	select_texture(tid);
	bind_default_flat_normal_map(); // no normal map
	qbd.draw_and_clear();
	if (do_blend) {disable_blend();}
	if (no_depth_write) {glDepthMask(GL_TRUE);}
}
void draw_emissive_billboards(quad_batch_draw &qbd, int tid) {
	if (qbd.empty()) return;
	shader_t s;
	s.begin_simple_textured_shader(); // unlit/emissive and textured
	set_additive_blend_mode();
	draw_billboards(qbd, tid);
	set_std_blend_mode();
}

class particle_texture_manager_t {
	// none, sparks, clouds, smoke, splash, bubble, droplet, steam
	string const fns[NUM_PART_EFFECTS] = {"", "", "", "", "water_splash.png", "white_circle.png", "white_circle.png"};
	int         tids[NUM_PART_EFFECTS] = {-1, BLUR_CENT_TEX, BLUR_TEX, SMOKE_PUFF_TEX, -1, -1, -1, SMOKE_PUFF_TEX};
public:
	int get_tid(unsigned effect) {
		assert(effect < NUM_PART_EFFECTS);
		int &tid(tids[effect]);
		if (tid < 0 && !fns[effect].empty()) {tid = get_texture_by_name(fns[effect]);} // load if needed
		return tid;
	}
};
particle_texture_manager_t particle_texture_manager;

void particle_manager_t::draw(shader_t &s, vector3d const &xlate) { // non-const because qbd is modified
	if (particles.empty() || !camera_pdu.cube_visible(get_bcube() + xlate)) return; // no particles are visible
	point const viewer_bs(camera_pdu.pos - xlate);
	vector<sphere_t> bubbles;

	for (particle_t const &p : particles) {
		if (!camera_pdu.sphere_visible_test((p.pos + xlate), p.radius)) continue; // VFC
		vector3d const vdir(viewer_bs - p.pos), up_dir((p.vel == zero_vector) ? plus_z : p.vel);
		vector3d v1(cross_product(vdir, up_dir).get_norm()*p.radius);
		vector3d v2(cross_product(v1,   vdir  ).get_norm()*p.radius);
		assert(p.effect < NUM_PART_EFFECTS);
		if      (p.effect == PART_EFFECT_SPARK  ) {v2 *= (1.0 + 1500.0*p.vel.mag());} // stretch in velocity dir
		else if (p.effect == PART_EFFECT_DROPLET) {v2 *= 3.0;} // stretch vertically
		if (p.effect == PART_EFFECT_BUBBLE) {bubbles.emplace_back(p.pos, p.radius);}
		else {qbds[p.effect].add_quad_dirs(p.pos, v1, v2, p.color, plus_z);} // use +z form the normal
	} // for p
	for (unsigned i = 0; i < NUM_PART_EFFECTS; ++i) {
		quad_batch_draw &qbd(qbds[i]);
		if (qbd.empty()) continue;
		int const tid(particle_texture_manager.get_tid(i));

		if (i == PART_EFFECT_SPARK) { // draw emissive particles with a custom shader
			draw_emissive_billboards(qbd, tid); // smooth alpha blended edges
			s.make_current();
		}
		else {draw_billboards(qbd, tid);} // no depth write, blend
	} // for i
	if (!bubbles.empty()) {
		// draw bubbles as spheres since billboards don't work well with the underwater postprocessing shader (due to depth and blend issues);
		// we can't make them transparent though because the underwater effect is run later and won't alpha blend properly
		// Note: can probably use instanced drawing here
		s.add_uniform_float("ambient_scale", 0.5);
		select_texture(WHITE_TEX);
		bind_default_flat_normal_map(); // no normal map
		s.set_cur_color(colorRGBA(0.6, 0.8, 1.0)); // blue-green tinted
		begin_sphere_draw(0); // untextured
		for (sphere_t const &b : bubbles) {draw_sphere_vbo(b.pos, b.radius, N_SPHERE_DIV, 0);} // textured=0
		end_sphere_draw();
		s.add_uniform_float("ambient_scale", building_ambient_scale); // reset
		check_mvm_update(); // needed after sphere drawing applies transforms
	}
}

void fire_manager_t::draw(shader_t &s, vector3d const &xlate) {
	if (fires.empty() || !camera_pdu.cube_visible(get_bcube() + xlate)) return; // no particles are visible
	point const viewer_bs(camera_pdu.pos - xlate);

	for (fire_t const &f : fires) {
		float const height(f.get_height());
		point const center(f.get_center());
		if (!camera_pdu.sphere_visible_test((center + xlate), max(f.radius, 0.5f*height))) continue; // VFC
		qbd.add_animated_billboard(center, viewer_bs, up_vector, WHITE, 1.5*f.radius, 0.5*height, fract(2.0f*f.time/TICKS_PER_SECOND));
	}
	draw_emissive_billboards(qbd, FIRE_TEX);
	s.make_current();
}

template<bool check_sz, typename T> bool are_pts_occluded_by_any_cubes(point const &pt, point const *const pts, unsigned npts,
	cube_t const &occ_area, T begin, T end, unsigned dim=3, float min_sz=0.0, float max_sep_dist=0.0)
{
	assert(npts > 0);

	for (auto c = begin; c != end; ++c) {
		if (check_sz && c->get_sz_dim(!dim) < min_sz) break; // too small an occluder; since cubes are sorted by size in this dim, we can exit the loop here
		if (dim <= 2 && (pt[dim] < c->d[dim][0]) == (pts[0][dim] < c->d[dim][0])) continue; // skip if cube face does not separate pt from the first point (dim > 2 disables)
		if (max_sep_dist > 0.0 && fabs(pt[dim] - c->get_center_dim(dim)) > max_sep_dist) continue; // check only one floor below/ceiling above
		if (!c->intersects(occ_area)) continue; // not between the object and viewer
		if (!check_line_clip(pt, pts[0], c->d)) continue; // first point does not intersect
		bool not_occluded(0);

		for (unsigned p = 1; p < npts; ++p) { // skip first point
			if (!check_line_clip(pt, pts[p], c->d)) {not_occluded = 1; break;}
		}
		if (!not_occluded) return 1;
	} // for c
	return 0;
}
template<bool check_sz, typename T> bool are_pts_occluded_by_any_cubes(point const &pt, point const *const pts, unsigned npts,
	cube_t const &occ_area, vector<T> const &cubes, unsigned dim=3, float min_sz=0.0, float max_sep_dist=0.0)
{
	return are_pts_occluded_by_any_cubes<check_sz>(pt, pts, npts, occ_area, cubes.begin(), cubes.end(), dim, min_sz, max_sep_dist);
}

car_t car_from_parking_space(room_object_t const &o) {
	rand_gen_t rgen;
	rgen.set_state(333*o.obj_id, o.obj_id+3);
	rgen.rand_mix();
	car_t car;
	car.dim     = o.dim;
	car.dir     = o.dir; // or random?
	car.cur_seg = o.obj_id; // store the random seed in car.cur_seg
	point center(o.get_cube_center());
	center[ o.dim] += 0.03*o.get_length()*rgen.signed_rand_float(); // small random misalign front/back
	center[!o.dim] += 0.05*o.get_width ()*rgen.signed_rand_float(); // small random misalign side
	car.set_bcube(point(center.x, center.y, o.z1()), get_nom_car_size());
	return car;
}
pair<cube_t, colorRGBA> car_bcube_color_from_parking_space(room_object_t const &o, unsigned btype) {
	car_t car(car_from_parking_space(o));
	set_car_model_color(car, btype);
	return make_pair(car.bcube, car.get_color());
}
bool check_cube_occluded(cube_t const &cube, vect_cube_t const &occluders, point const &viewer) {
	if (occluders.empty()) return 0;
	point pts[8];
	unsigned const npts(get_cube_corners(cube.d, pts, viewer, 0)); // should return only the 6 visible corners
	cube_t occ_area(cube);
	occ_area.union_with_pt(viewer); // any occluder must intersect this cube
	return are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, occluders); // no dim because cubes are of mixed dim and we can't use that optimization
}
struct comp_car_by_dist {
	vector3d const &viewer;
	comp_car_by_dist(vector3d const &viewer_) : viewer(viewer_) {}
	bool operator()(car_t const &c1, car_t const &c2) const {
		return (p2p_dist_xy_sq(c1.bcube.get_cube_center(), viewer) > p2p_dist_xy_sq(c2.bcube.get_cube_center(), viewer));
	}
};

bool building_t::has_cars_to_draw(bool player_in_building) const {
	if (!has_cars_enabled()) return 0;
	if (!has_room_geom())    return 0;
	if (is_parking())        return 1; // parking structure
	if (player_in_building && has_parking_garage) return 1; // parking garage cars are drawn if the player is in the building
	if (interior->room_geom->has_garage_car)      return 1; // have car in a garage
	return 0;
}
void building_t::draw_cars_in_building(shader_t &s, vector3d const &xlate, bool player_in_this_building, bool shadow_only) const {
	if (!has_room_geom()) return; // can get here in rare cases, maybe only shadow_only pass
	point viewer(camera_pdu.pos - xlate); // building space
	bool const check_occlusion(display_mode & 0x08), is_parking_str(is_parking());
	float const floor_spacing(get_window_vspace());
	vect_room_object_t const &objs(interior->room_geom->objs);
	unsigned const pg_wall_start(interior->room_geom->wall_ps_start);
	assert(pg_wall_start < objs.size());
	static vector<car_t> cars_to_draw; // reused across frames
	cars_to_draw.clear();

	if (interior->room_geom->has_garage_car && player_in_basement < 2) { // car in a house garage, and player not in the basement
		room_object_t const &obj(objs[pg_wall_start]);
		assert(obj.type == TYPE_PARK_SPACE); // must be a parking space

		// skip if viewer is in this building and on a different floor
		if (obj.is_used() && (!player_in_this_building || !check_occlusion || !has_int_garage ||
			int((viewer.z - bcube.z1())/floor_spacing) == int((obj.z2() - bcube.z1())/floor_spacing)))
		{
			car_t const car(car_from_parking_space(obj));
			if (camera_pdu.cube_visible(car.bcube + xlate)) {cars_to_draw.push_back(car);}
		}
	}
	else if ((player_in_this_building && has_parking_garage) || is_parking_str) { // cars in parking garage that the player is in, or in parking structure
		int const camera_floor(get_floor_for_zval(viewer.z));
		float max_vis_zval(ground_floor_z1);

		if (is_parking_str) {max_vis_zval = bcube.z2();}
		else { // underground parking garage
			// only draw if the player or light is in the basement, or the player is on the first floor where a car may be visible through the stairs
			if (!shadow_only) {max_vis_zval += floor_spacing;} // player on first floor?
			if (viewer.z > max_vis_zval) return;
		}
		viewer = get_inv_rot_pos(viewer); // not needed because there are no cars in rotated buildings?
		auto objs_ps_end(objs.begin() + interior->room_geom->wall_ps_end);
		assert(objs_ps_end <= objs.end());
		vect_cube_t occluders; // should this be split out per PG level?

		// start at walls, since parking spaces are added after those, and continue across all parking garage levels
		for (auto i = objs.begin()+pg_wall_start; i != objs_ps_end; ++i) {
			if (check_occlusion && (i->type == TYPE_PG_WALL || (is_parking_str && i->type == TYPE_STAIR_WALL))) {occluders.push_back(*i);}
			if (i->type != TYPE_PARK_SPACE || !i->is_used())                       continue; // not a space, or no car in this space
			if (player_in_this_building && i->z2() < viewer.z - 2.0*floor_spacing) continue; // more than a floor below - skip
			car_t car(car_from_parking_space(*i));
			int const car_floor(get_floor_for_zval(car.bcube.zc())), num_floors_apart(abs(camera_floor - car_floor));

			if (shadow_only) {
				if (num_floors_apart > 0) continue; // different floors
			}
			else if (check_occlusion && player_in_this_building) {
				if (num_floors_apart > 1) continue; // more than one floor apart
				if (num_floors_apart > 0 && !player_in_mall && !line_intersect_stairs_or_ramp(viewer, car.get_center(), 1)) continue; // skip_u_stairs=1
			}
			if (!shadow_only && player_in_basement == 3 && !is_cube_visible_through_extb_door(viewer, car.bcube)) continue; // player in extb; cars only visible through door
			if (camera_pdu.cube_visible(car.bcube + xlate)) {cars_to_draw.push_back(car);}
		} // for i
		if (cars_to_draw.empty()) return;

		if (check_occlusion) {
			// gather occluders from parking garage ceilings and floors (below ground floor or in parking structure)
			for (auto const &ceiling : interior->fc_occluders) {
				if (ceiling.z1() <= max_vis_zval) {occluders.push_back(ceiling);}
			}
			auto in(cars_to_draw.begin()), out(in);

			for (; in != cars_to_draw.end(); ++in) { // filter out occluded cars
				if (!check_cube_occluded(in->bcube, occluders, viewer)) {*(out++) = *in;}
			}
			cars_to_draw.erase(out, cars_to_draw.end());
		} // end check_occlusion
		std::sort(cars_to_draw.begin(), cars_to_draw.end(), comp_car_by_dist(viewer)); // required for correct window alpha blending
	}
	if (cars_to_draw.empty()) return;
	
	if (!s.is_setup()) { // caller didn't set up the shader
		if (shadow_only) {s.begin_shadow_map_shader();} // this path should be unused
		else {setup_building_draw_shader(s, 0.0, 1, 0, 0);} // min_alpha=0.0, enable_indir=1, force_tsl=0, use_texgen=0, water_damage=0.0
	}
	for (auto &car : cars_to_draw) {draw_car_in_pspace(car, s, xlate, shadow_only, btype);}
	check_mvm_update(); // needed after popping model transform matrix
}

void append_line_pt(vector<vert_wrap_t> &line_pts, point const &pos) {
	if (line_pts.size() > 1) {line_pts.emplace_back(line_pts.back());} // duplicate point to create a line segment
	line_pts.emplace_back(pos);
}
void building_t::debug_people_in_building(shader_t &s, point const &camera_bs) const {
	if (!has_people() && !DEBUG_AI_COLLIDERS) return;
	shader_t color_shader;
	color_shader.begin_color_only_shader(YELLOW);
	vector<vert_wrap_t> line_pts;
	unsigned const ndiv = 16;

	for (person_t const &p : interior->people) {
		// use different colors if the path uses the nav grid or was shortened
		colorRGBA const &path_color(p.path.uses_nav_grid ? (p.path.is_shortened ? LT_BLUE : GREEN) : (p.path.is_shortened ? ORANGE : YELLOW));
		color_shader.set_cur_color(path_color);
		for (point const &v : p.path) {append_line_pt(line_pts, v);}
		if (p.target_valid ())        {append_line_pt(line_pts, p.target_pos);} // next target - not dest
		if (!line_pts.empty())        {append_line_pt(line_pts, p.pos);} // add starting point if there's a valid path
		if (line_pts.size() > 1)      {draw_verts    (line_pts, GL_LINES);}
		line_pts.clear();
		float const sradius(0.25*p.radius);
		point from;

		for (path_pt_t const &v : p.path) { // Note: backwards
			if (v == p.path.front()) {from = v; continue;} // skip first point
			color_shader.set_cur_color(path_color*(v.fixed ? 0.5 : 1.0));
			draw_sphere_vbo(v, sradius, ndiv, 0);
			point const arrow_start(v + (2.0*sradius)*(from - v).get_norm());
			draw_fast_cylinder(v, arrow_start, 0.5*sradius, 0.0, ndiv, 0);
			from = v;
		}
		assert(p.goal_type < NUM_GOAL_TYPES);
		colorRGBA const goal_colors[NUM_GOAL_TYPES] = {BLACK, BLUE, PINK, MAGENTA, RED, ORANGE, PURPLE}; // NONE, ROOM, ELEVATOR, ESCALATOR, PLAYER, PLAYER_LAST_POS, SOUND
		color_shader.set_cur_color(goal_colors[p.goal_type]);
		if (!p.path.empty ()) {draw_sphere_vbo(p.path.front(), sradius, ndiv, 0);} // draw last point/dest
		if (p.target_valid()) {draw_sphere_vbo(p.target_pos,   sradius, ndiv, 0);} // draw target pos
	} // for p
	if (DEBUG_AI_COLLIDERS && (frame_counter & 1)) { // debug avoid cubes on alternating frames
		vect_cube_t avoid;
		interior->get_avoid_cubes(avoid, (camera_bs.z - get_bldg_player_height()), (camera_bs.z + CAMERA_RADIUS),
			CAMERA_RADIUS*global_building_params.player_coll_radius_scale, get_floor_thickness(), get_floor_ceil_gap(), 1, 0); // same_as_player=1, skip_stairs=0
		color_shader.set_cur_color(RED);
		//for (cube_t &c : avoid) {c.expand_by_xy(0.75*interior->people.front().radius);} // expand by ref person's radius
		for (cube_t const &c : avoid) {draw_simple_cube(c);}
	}
	color_shader.end_shader();
	s.make_current();
}

bool check_custom_occluder_cull(cube_t const &c, point const &viewer, vect_cube_t const &occluders, bool occluder_dim) {
	// check any extra occluders from special rooms such as mall stores
	// Note: ignores building rotate since malls should not be in a rotated building
	point pts[8];
	unsigned const npts(get_cube_corners(c.d, pts, viewer, 0)); // should return only the 6 visible corners
	cube_t occ_area(c);
	occ_area.union_with_pt(viewer); // any occluder must intersect this cube
	return are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, occluders, occluder_dim);
}

// Note: c is in local building space and viewer_in is in non-rotated building space
bool building_t::check_obj_occluded(cube_t const &c, point const &viewer_in, occlusion_checker_noncity_t const &oc,
	bool reflection_pass, bool c_is_building_part, bool skip_basement_check) const
{
	if (!interior) return 0; // could probably make this an assert
	//highres_timer_t timer("Check Object Occlusion"); // 0.001ms
	point const viewer(get_inv_rot_pos(viewer_in)); // rotate viewer pos into building space
	bool const player_in_building(reflection_pass || point_in_building_or_basement_bcube(viewer)); // if reflection pass, assume the player is in this building
	bool const targ_in_basement(c.z2() <= ground_floor_z1);
	float const floor_spacing(get_window_vspace()), ground_floor_ceiling(ground_floor_z1 + floor_spacing);
	bool checked_conn_ret(0);
	
	if (targ_in_basement && !skip_basement_check) { // fully inside basement
		if (viewer.z > ground_floor_ceiling) return 1; // viewer not on first floor
		
		if (!player_in_building) { // player not in this building
			checked_conn_ret = (camera_in_building && player_in_basement && interior_visible_from_other_building_ext_basement(oc.get_xlate(), oc.query_is_for_light));
			if (!checked_conn_ret) return 1;
		}
	}
	point pts[8];
	unsigned const npts(get_cube_corners(c.d, pts, viewer, 0)); // should return only the 6 visible corners

	if (reflection_pass && clip_plane.w != 0.0 && !is_rotated()) { // check for all visible corners of object behind the clip plane
		bool any_not_behind(0);

		for (unsigned n = 0; n < npts; ++n) {
			if (dot_product(pts[n], clip_plane) + clip_plane.w > 0.0) {any_not_behind = 1; break;}
		}
		if (!any_not_behind) return 1;
	}
	cube_t occ_area(c);
	occ_area.union_with_pt(viewer); // any occluder must intersect this cube
	point const center(c.get_cube_center());

	if (!c_is_building_part && viewer.z > ground_floor_z1) { // viewer above ground
		if (has_retail()) {
			cube_t const &retail_part(get_retail_part());

			if (viewer.z < retail_part.z2() && retail_part.contains_pt(center)) {
				// both the object and the viewer are in the ground floor of a retail building - check shelf rack backs as occluders
				if (reflection_pass) { // use actual camera, not camera reflected in glass floor; should be correct for vertical reflection
					if (pre_reflect_camera_pos_bs != zero_vector && !is_rotated()) { // have pre-reflect pos; rotated building should not get into this case
						cube_t occ_area2(c);
						occ_area2.union_with_pt(pre_reflect_camera_pos_bs);
						if (check_shelfrack_occlusion(pre_reflect_camera_pos_bs, pts, get_cube_corners(c.d, pts, pre_reflect_camera_pos_bs, 0), occ_area2)) return 1;
					}
				}
				else {
					if (check_shelfrack_occlusion(viewer, pts, npts, occ_area)) return 1;
				}
				if (retail_part.contains_pt(viewer)) return 0; // no walls/ceilings/floors in retail area; done
			}
		}
		else if (is_warehouse() && point_in_industrial(viewer) && point_in_industrial(center)) { // both the object and the viewer are in the warehouse
			if (check_warehouse_shelf_occlusion(viewer, pts, npts, occ_area)) return 1; // blocked by a shelf
		}
	}
	vector3d const dir(viewer - center);
	bool const pri_dim(fabs(dir.x) < fabs(dir.y));
	
	if (!c_is_building_part && !reflection_pass) {
		// check walls of this building; not valid for reflections because the reflected camera may be on the other side of a wall/mirror
		if (targ_in_basement && is_pos_in_pg_or_backrooms(viewer) && is_pos_in_pg_or_backrooms(center)) { // object and pos are both in the parking garage or backrooms
			if (interior->has_backrooms) {
				// check for occlusion from the wall segments on either side of the extended basement door that separates it from the basement
				bool const dim(interior->extb_wall_dim), dir(interior->extb_wall_dir);
				door_t const &door(interior->get_ext_basement_door());
				cube_t wall(get_basement());
				min_eq(wall.z1(), interior->basement_ext_bcube.z1()); // cover both basement and ext basement in Z
				float const wall_pos(wall.d[dim][dir]);
				wall.d[dim][!dir] = wall_pos;
				wall.d[dim][ dir] = wall_pos + (dir ? 1.0 : -1.0)*get_wall_thickness(); // extend slightly into ext basement
				assert(wall.is_strictly_normalized());
				cube_t walls[2] = {wall, wall}; // lo, hi
				for (unsigned d = 0; d < 2; ++d) {walls[d].d[!dim][!d] = door.d[!dim][d];} // exclude the door
				if (are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, walls, walls+2, dim, 0.0)) return 1; // no size check
			}
			if (check_pg_br_wall_occlusion(viewer, pts, npts, occ_area, dir)) return 1;
		}
		else { // regular walls case, with size check (helps with light bcubes)
			bool const in_ext_basement(point_in_extended_basement_not_basement(viewer) || point_in_extended_basement_not_basement(center));

			for (unsigned D = 0; D < 2; ++D) {
				bool const d(bool(D) ^ pri_dim); // try primary dim first
				float const min_sz(min(floor_spacing, 0.5f*c.get_sz_dim(!d))); // account for perspective; min with floor_spacing to allow for room sized queries
				// exclude mall back hallway walls, since they are a mostly either exterior walls and duplicated with store walls
				unsigned const extb_walls_start(interior->extb_walls_start[d]), extb_walls_end(interior->mall_hall_walls_start[d]);
				vect_cube_t const &walls(interior->walls[d]);
				assert(extb_walls_start <= walls.size());
				if (are_pts_occluded_by_any_cubes<1>(viewer, pts, npts, occ_area, walls.begin()+(in_ext_basement ? extb_walls_start : 0),
					walls.begin()+(in_ext_basement ? extb_walls_end : extb_walls_start), d, min_sz)) return 1;
			} // for D
		}
	}
	if (!c_is_building_part && player_in_building) {
		// viewer inside this building; includes shadow_only case and reflection_pass (even if reflected camera is outside the building);
		// okay for retail glass floor reflections because the reflection won't cross an opaque floor or ceiling in fc_occluders
		// check floors/ceilings of this building
		if (fabs(viewer.z - c.zc()) > (reflection_pass ? 1.0 : 0.5)*floor_spacing) { // on different floors
			float max_sep_dist(floor_spacing);
			
			if (point_in_mall(viewer) || point_in_mall(center)) { // taller ceilings in malls
				max_sep_dist = get_mall_floor_spacing();
			}
			else if (has_tall_retail()) { // handle ceilings more than one part tall
				cube_t const &retail_part(get_retail_part());
				if (retail_part.contains_pt(viewer) || retail_part.contains_pt(center)) {max_sep_dist *= retail_floor_levels;}
			}
			if (are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, interior->fc_occluders, 2, 0.0, max_sep_dist)) return 1;
		}
	}
	else if (camera_in_building) { // player in some other building
		if (checked_conn_ret  ) return 0; // player in ext basement connected to this building; skip below checks in this case
		if (player_in_basement) return 1; // if player is in the basement of a different building, they probably can't see an object in this building
		if (player_in_windowless_building()) return 1; // player inside another windowless office building, objects in this building not visible
		if (is_rotated()) return 0; // not implemented yet - need to rotate viewer and pts into coordinate space of player_building

		if (player_building != nullptr && player_building->interior) { // check walls of the building the player is in
			if (player_building != this) { // otherwise player_in_this_building should be true; note that we can get here from building_t::add_room_lights()
				if (player_building->has_ind_info() && viewer.z > player_building->ground_floor_z1) {
					cube_with_ix_t const *nww(player_building->interior->ind_info->non_window_walls);

					for (unsigned n = 0; n < 3; ++n) {
						cube_with_ix_t const &c(nww[n]);
						if (building_has_open_ext_door && (c.ix & 4)) continue; // open exterior door, and this wall has an exterior door; skip it
						bool const dim(c.ix & 2);
						if ((viewer[dim] < c.d[dim][0]) == (pts[0][dim] < c.d[dim][0]) || !c.intersects(occ_area) || !check_line_clip(viewer, pts[0], c.d)) continue;
						bool not_occluded(0);

						for (unsigned p = 1; p < npts; ++p) { // skip first point
							if (!check_line_clip(viewer, pts[p], c.d)) {not_occluded = 1; break;}
						}
						if (!not_occluded) return 1;
					} // for n
				}
				for (unsigned D = 0; D < 2; ++D) {
					bool const d(bool(D) ^ pri_dim); // try primary dim first
					unsigned const extb_walls_start(player_building->interior->extb_walls_start[d]);
					vect_cube_t const &walls(player_building->interior->walls[d]); // skip extended basement walls, since no other buildings will be visible from there
					assert(extb_walls_start <= walls.size());
					if (are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, walls.begin(), walls.begin()+extb_walls_start, d)) return 1;
				} // for D
				if (fabs(viewer.z - c.zc()) > 0.5*floor_spacing) { // check floors and ceilings of the building the player is in
					if (are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, player_building->interior->fc_occluders, 2, 0.0, floor_spacing)) return 1;
				}
			}
		}
	}
	else if (viewer.z < bcube.z2()) { // player not in a building and not above this building
		if (is_rotated())                       return 0; // not implemented yet - c is not an axis aligned cube in global coordinate space
		if (has_windows() && oc.is_occluded(c)) return 1; // check other buildings; not needed for windowless buildings since they shouldn't be drawn (and wrong for walkways)
	}
	if (!c_is_building_part && viewer.z > ground_floor_ceiling && is_cube()) {
		// player above first floor of this building; check if object is occluded by a roof; we don't check bcube.z2() becase a lower part roof may be an occluder
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			float const roof_z(p->z2());
			if (viewer.z < roof_z) continue; // viewer below the roof
			if (c.z2()   > roof_z) continue; // object above the roof
			if (p->contains_pt_xy(viewer) && p->contains_pt_xy(c.get_cube_center())) continue; // skip stacked parts case; should be handled by floor occluders
			cube_t roof(*p);
			roof.z1() = roof_z - get_fc_thickness();
			
			// check if on top floor of roof with a skylight, or on the stairs of the floor below; could split the roof in this case, but that may not make much of a difference
			if ((c.z2() > (roof_z - floor_spacing) || (c.z2() > (roof_z - 2.0f*floor_spacing) && check_cube_on_or_near_stairs(c))) && check_skylight_intersection(roof)) {
				continue;
			}
			bool not_occluded(0);

			for (unsigned p = 0; p < npts; ++p) {
				if (!check_line_clip(viewer, pts[p], roof.d)) {not_occluded = 1; break;}
			}
			if (!not_occluded) return 1;
		} // for p
	}
	return 0;
}
bool building_t::check_pg_br_wall_occlusion(point const &viewer, point const *const pts, unsigned npts, cube_t const &occ_area, vector3d const &view_dir) const {
	if (!has_room_geom()) return 0;
	bool const pri_dim(fabs(view_dir.x) < fabs(view_dir.y));
	index_pair_t start, end;
	get_pgbr_wall_ix_for_pos(viewer, start, end);
	// in cases where we clipped under the building the range will be empty and there will be no occlusion

	for (unsigned D = 0; D < 2; ++D) {
		bool const d(bool(D) ^ pri_dim); // try primary dim first
		vect_cube_t const &walls(interior->room_geom->pgbr_walls[d]);
		if (are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, walls.begin()+start.ix[d], walls.begin()+end.ix[d], d, 0.0)) return 1; // no size check
	}
	return 0;
}
bool building_t::check_shelfrack_occlusion(point const &viewer, point const *const pts, unsigned npts, cube_t const &occ_area) const {
	vect_cube_t const &back_occ(interior->room_geom->shelf_rack_occluders[0]), &top_acc(interior->room_geom->shelf_rack_occluders[1]);
	if (!has_room_geom() || back_occ.empty()) return 0;
	// if viewer (maybe light) is above shelf racks (first one should be on the ground floor), use tops as occluders
	if (viewer.z > back_occ.front().z2() && are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, top_acc, 2)) return 1;
	bool const long_dim(get_retail_long_dim());
	return are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, back_occ, !long_dim); // use backs as occluders
}
bool building_t::check_warehouse_shelf_occlusion(point const &viewer, point const *const pts, unsigned npts, cube_t const &occ_area) const {
	if (!has_room_geom() || !interior->ind_info || interior->room_geom->shelf_rack_occluders[0].empty()) return 0;
	return are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occ_area, interior->room_geom->shelf_rack_occluders[0], !interior->ind_info->entrance_dim);
}

bool building_t::is_entire_building_occluded(point const &viewer, occlusion_checker_noncity_t const &oc) const {
	if (is_rotated()) return 0; // not handled (optimization)

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
		if (has_basement() && (i - parts.begin()) == (int)basement_part_ix) continue; // skip the basement, which isn't visible from outside the building
		if (!check_obj_occluded(*i, viewer, oc, 0, 1)) return 0; // c_is_building_part=1
	}
	if (interior_visible_from_other_building_ext_basement(oc.get_xlate(), oc.query_is_for_light)) return 0;
	return 1; // all parts occluded
}

bool paint_draw_t::have_any_sp() const {
	for (unsigned i = 0; i <= NUM_SP_EMISSIVE_COLORS; ++i) {
		if (!sp_qbd[i].empty()) return 1;
	}
	return 0;
}
quad_batch_draw &paint_draw_t::get_paint_qbd(bool is_marker, unsigned emissive_color_id) {
	if (is_marker) return m_qbd;
	return sp_qbd[(emissive_color_id == 0) ? 0 : (1 + ((emissive_color_id-1)%NUM_SP_EMISSIVE_COLORS))]; // choose spraypaint emissive color
}
void paint_draw_t::draw_paint(shader_t &s) const {
	bool const have_sp(have_any_sp());
	if (!have_sp && m_qbd.empty()) return; // nothing to do
	glDepthMask(GL_FALSE); // disable depth write
	enable_blend();
	bind_default_flat_normal_map(); // no normal map

	if (have_sp) {
		select_texture(BLUR_CENT_TEX); // spraypaint - smooth alpha blended edges

		for (unsigned i = 0; i <= NUM_SP_EMISSIVE_COLORS; ++i) {
			quad_batch_draw const &qbd(sp_qbd[i]);
			if (qbd.empty()) continue;
			if (i > 0) {s.set_color_e(sp_emissive_colors[i-1]);}
			qbd.draw();
			if (i > 0) {s.clear_color_e();}
		} // for i
	}
	if (!m_qbd.empty()) {
		select_texture(get_texture_by_name("circle.png", 0, 0, 1, 0.0, 1, 1, 1)); // markers - sharp edges, used as alpha mask with white background color
		m_qbd.draw();
	}
	disable_blend();
	glDepthMask(GL_TRUE);
}
void paint_draw_t::clear() {
	for (unsigned i = 0; i <= NUM_SP_EMISSIVE_COLORS; ++i) {sp_qbd[i].clear();}
	m_qbd.clear();
}

void building_decal_manager_t::commit_pend_tape_qbd() {
	pend_tape_qbd.add_quads(tape_qbd);
	pend_tape_qbd.clear();
}
void building_decal_manager_t::add_burn_spot(point const &pos, float radius) {
	// if there are too many existing spots, remove the first (oldest) one; this can be slow, but shouldn't happen very often
	unsigned const max_spots = 100;
	unsigned const num_spots(burn_qbd.verts.size()/6); // 6 verts/2 triangles per quad
	if (num_spots >= max_spots) {burn_qbd.verts.erase(burn_qbd.verts.begin(), (burn_qbd.verts.begin() + 6*(num_spots - max_spots + 1)));}
	burn_qbd.add_quad_dirs(pos, -plus_x*radius, plus_y*radius, BLACK); // -x!
}
void building_decal_manager_t::add_blood_or_stain(point const &pos, float radius, colorRGBA const &color, bool is_blood, unsigned dim, bool dir) {
	assert(dim <= 2);
	tex_range_t tex_range(tex_range_t::from_atlas((rgen.rand()&1), (rgen.rand()&1), 2, 2)); // 2x2 texture atlas
	tex_range.swap_xy = rgen.rand_bool();
	vector3d ds, dt, dn;
	dn[dim] = (dir ? 1.0 : -1.0);
	ds[(dim+1)%3] = -1.0;
	dt[(dim+2)%3] =  1.0;
	if (!dir) {swap(ds, dt);} // use correct winding order
	blood_qbd[!is_blood].add_quad_dirs(pos, ds*radius, dt*radius, color, dn, tex_range); // -x!
}
void building_decal_manager_t::draw_building_interior_decals(shader_t &s, bool player_in_building, bool shadow_only) const {
	if (shadow_only) { // shadow pass, draw tape only
		if (player_in_building) {
			tape_qbd.draw(); // somewhat inefficient, since we have to send all the data for every light source
			pend_tape_qbd.draw();
		}
		return;
	}
	paint_draw[1].draw_paint(s); // draw exterior paint always - this will show up on windows (even when looking outside into another part of the same building)
	if (!player_in_building) return;
	paint_draw[0].draw_paint(s); // draw interior paint

	if (!tp_qbd.empty()) { // toilet paper squares: double sided, lit from top
		glDisable(GL_CULL_FACE); // draw both sides
		select_texture(WHITE_TEX);
		select_texture(get_toilet_paper_nm_id(), 5); // apply normal map
		tp_qbd.draw(); // use a VBO for this if the player leaves the building and then comes back?
		bind_default_flat_normal_map(); // no normal map
		glEnable(GL_CULL_FACE);
	}
	if (!tape_qbd.empty() || !pend_tape_qbd.empty()) { // tape lines: single sided so that lighting works, both sides drawn independently
		select_texture(WHITE_TEX);
		tape_qbd.draw();
		pend_tape_qbd.draw();
	}
	if (!blood_qbd[0].empty() || !blood_qbd[1].empty() || !glass_qbd.empty() || !burn_qbd.empty()) { // draw alpha blended decals
		glDepthMask(GL_FALSE); // disable depth write
		enable_blend();
		int const blood_tids[2] = {BLOOD_SPLAT_TEX, get_texture_by_name("atlas/blood_white.png")};

		for (unsigned i = 0; i < 2; ++i) {
			if (blood_qbd[i].empty()) continue;
			select_texture(blood_tids[i]);
			blood_qbd[i].draw();
		}
		if (!glass_qbd.empty()) {
			select_texture(get_texture_by_name("interiors/broken_glass.png"));
			glass_qbd.draw();
		}
		if (!burn_qbd.empty()) {
			select_texture(BLUR_CENT_TEX);
			burn_qbd.draw();
		}
		disable_blend();
		glDepthMask(GL_TRUE);
	}
}


