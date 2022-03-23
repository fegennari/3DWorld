// 3D World - Building Interior Room Item Drawing
// by Frank Gennari 4/17/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#include "subdiv.h" // for sd_sphere_d
#include "profiler.h"

unsigned const MAX_ROOM_GEOM_GEN_PER_FRAME = 1;
colorRGBA const rat_color(GRAY); // make the rat's fur darker

object_model_loader_t building_obj_model_loader;

extern bool camera_in_building;
extern int display_mode, frame_counter, animate2, player_in_basement;
extern float office_chair_rot_rate, cur_dlight_pcf_offset;
extern point pre_smap_player_pos;
extern cube_t smap_light_clip_cube;
extern pos_dir_up camera_pdu;
extern building_t const *player_building;
extern carried_item_t player_held_object;
extern building_params_t global_building_params;

unsigned get_num_screenshot_tids();
tid_nm_pair_t get_phone_tex(room_object_t const &c);
template< typename T > void gen_quad_ixs(vector<T> &ixs, unsigned size, unsigned ix_offset);
void draw_car_in_pspace(car_t &car, shader_t &s, vector3d const &xlate, bool shadow_only);

bool has_key_3d_model() {return building_obj_model_loader.is_model_valid(OBJ_MODEL_KEY);}

colorRGBA room_object_t::get_model_color() const {return building_obj_model_loader.get_avg_color(get_model_id());}

// skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2 to match CSG cube flags
void rgeom_mat_t::add_cube_to_verts(cube_t const &c, colorRGBA const &color, point const &tex_origin,
	unsigned skip_faces, bool swap_tex_st, bool mirror_x, bool mirror_y, bool inverted)
{
	//assert(c.is_normalized()); // no, bathroom window is denormalized
	vertex_t v;
	v.set_c4(color);

	// Note: stolen from draw_cube() with tex coord logic, back face culling, etc. removed
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions, drawn as {Z, X, Y}
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (skip_faces & (1 << (2*(2-n) + j))) continue; // skip this face
			v.set_ortho_norm(n, (bool(j) ^ inverted));
			v.v[n] = c.d[n][j];

			for (unsigned s1 = 0; s1 < 2; ++s1) {
				v.v[d[1]] = c.d[d[1]][s1];
				v.t[swap_tex_st] = ((tex.tscale_x == 0.0) ? float(s1) : tex.tscale_x*(v.v[d[1]] - tex_origin[d[1]])); // tscale==0.0 => fit texture to cube

				for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
					bool const s2(bool(k^j^s1)^inverted^1); // need to orient the vertices differently for each side
					v.v[d[0]] = c.d[d[0]][s2];
					v.t[!swap_tex_st] = ((tex.tscale_y == 0.0) ? float(s2) : tex.tscale_y*(v.v[d[0]] - tex_origin[d[0]]));
					quad_verts.push_back(v);
					if (mirror_x) {quad_verts.back().t[0] = 1.0 - v.t[0];} // use for pictures and books
					if (mirror_y) {quad_verts.back().t[1] = 1.0 - v.t[1];} // used for books
				} // for k
			} // for s1
		} // for j
	} // for i
}
// untextured version of the above function
void rgeom_mat_t::add_cube_to_verts_untextured(cube_t const &c, colorRGBA const &color, unsigned skip_faces) { // add an inverted flag?
	vertex_t v;
	v.set_c4(color);

	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (skip_faces & (1 << (2*(2-n) + j))) continue; // skip this face
			v.set_ortho_norm(n, j);
			v.v[n] = c.d[n][j];

			for (unsigned s1 = 0; s1 < 2; ++s1) {
				v.v[d[1]] = c.d[d[1]][s1];
				v.v[d[0]] = c.d[d[0]][!(j^s1)]; quad_verts.push_back(v);
				v.v[d[0]] = c.d[d[0]][ (j^s1)]; quad_verts.push_back(v);
			} // for s1
		} // for j
	} // for i
}

template<typename T> void add_inverted_triangles(T &verts, vector<unsigned> &indices, unsigned verts_start, unsigned ixs_start) {
	unsigned const verts_end(verts.size()), numv(verts_end - verts_start);
	verts.resize(verts_end + numv);

	for (unsigned i = verts_start; i < verts_end; ++i) {
		verts[i+numv] = verts[i];
		verts[i+numv].invert_normal();
	}
	unsigned const ixs_end(indices.size()), numi(ixs_end - ixs_start);
	indices.resize(ixs_end + numi);
	for (unsigned i = 0; i < numi; ++i) {indices[ixs_end + i] = (indices[ixs_end - i - 1] + numv);} // copy in reverse order
}

void swap_cube_z_xy(cube_t &c, bool dim) {
	swap(c.z1(), c.d[dim][0]);
	swap(c.z2(), c.d[dim][1]);
}

void rgeom_mat_t::add_ortho_cylin_to_verts(cube_t const &c, colorRGBA const &color, int dim, bool draw_bot, bool draw_top, bool two_sided, bool inv_tb,
	float rs_bot, float rs_top, float side_tscale, float end_tscale, bool skip_sides, unsigned ndiv, float side_tscale_add, bool swap_txy)
{
	if (dim == 2) { // Z: this is our standard v_cylinder
		add_vcylin_to_verts(c, color, draw_bot, draw_top, two_sided, inv_tb, rs_bot, rs_top, side_tscale, end_tscale, skip_sides, ndiv, side_tscale_add, swap_txy);
		return;
	}
	cube_t c_rot(c);
	swap_cube_z_xy(c_rot, dim);
	unsigned const itri_verts_start_ix(itri_verts.size()), ixs_start_ix(indices.size());
	add_vcylin_to_verts(c_rot, color, draw_bot, draw_top, two_sided, inv_tb, rs_bot, rs_top, side_tscale, end_tscale, skip_sides, ndiv, side_tscale_add, swap_txy);
	
	for (auto v = itri_verts.begin()+itri_verts_start_ix; v != itri_verts.end(); ++v) { // swap triangle vertices and normals
		std::swap(v->v[2], v->v[dim]);
		std::swap(v->n[2], v->n[dim]);
	}
	std::reverse(indices.begin()+ixs_start_ix, indices.end()); // fix winding order
}
void rgeom_mat_t::add_vcylin_to_verts(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top, bool two_sided, bool inv_tb,
	float rs_bot, float rs_top, float side_tscale, float end_tscale, bool skip_sides, unsigned ndiv, float side_tscale_add, bool swap_txy)
{
	point const center(c.get_cube_center());
	float const radius(0.5*min(c.dx(), c.dy())); // cube X/Y size should be equal/square
	add_cylin_to_verts(point(center.x, center.y, c.z1()), point(center.x, center.y, c.z2()), radius*rs_bot, radius*rs_top,
		color, draw_bot, draw_top, two_sided, inv_tb, side_tscale, end_tscale, skip_sides, ndiv, side_tscale_add, swap_txy);
}
void rgeom_mat_t::add_cylin_to_verts(point const &bot, point const &top, float bot_radius, float top_radius, colorRGBA const &color, bool draw_bot,
	bool draw_top, bool two_sided, bool inv_tb, float side_tscale, float end_tscale, bool skip_sides, unsigned ndiv, float side_tscale_add, bool swap_txy)
{
	assert((!skip_sides) || draw_bot || draw_top); // must draw something
	point const ce[2] = {bot, top};
	float const ndiv_inv(1.0/ndiv), half_end_tscale(0.5*end_tscale);
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, bot_radius, top_radius, ndiv, v12));
	color_wrapper const cw(color);
	unsigned itris_start(itri_verts.size()), ixs_start(indices.size()), itix(itris_start), iix(ixs_start);

	if (!skip_sides) {
		itri_verts.resize(itris_start + 2*(ndiv+1));
		indices.resize(ixs_start + 6*ndiv);
		unsigned const ixs_off[6] = {1,2,0, 3,2,1}; // 1 quad = 2 triangles

		for (unsigned i = 0; i <= ndiv; ++i) { // vertex data
			unsigned const s(i%ndiv);
			float const ts(side_tscale*(1.0f - i*ndiv_inv) + side_tscale_add);
			norm_comp const normal(0.5*(vpn.n[s] + vpn.n[(i+ndiv-1)%ndiv])); // normalize?
			itri_verts[itix++].assign(vpn.p[(s<<1)+0], normal, (swap_txy ? 0.0 : ts), (swap_txy ? ts : 0.0), cw);
			itri_verts[itix++].assign(vpn.p[(s<<1)+1], normal, (swap_txy ? 1.0 : ts), (swap_txy ? ts : 1.0), cw);
		}
		for (unsigned i = 0; i < ndiv; ++i) { // index data
			unsigned const ix0(itris_start + 2*i);
			for (unsigned j = 0; j < 6; ++j) {indices[iix++] = ix0 + ixs_off[j];}
		}
		// room object drawing uses back face culling and single sided lighting; to make lighting two sided, need to add verts with inverted normals/winding dirs
		if (two_sided) {add_inverted_triangles(itri_verts, indices, itris_start, ixs_start);}
	}
	// maybe add top and bottom end cap using triangles, currently using all TCs=0.0
	unsigned const num_ends((unsigned)draw_top + (unsigned)draw_bot);
	itris_start = itix = itri_verts.size();
	ixs_start   = iix  = indices.size();
	itri_verts.resize(itris_start + (ndiv + 1)*num_ends);
	indices.resize(ixs_start + 3*ndiv*num_ends);

	for (unsigned bt = 0; bt < 2; ++bt) {
		if (!(bt ? draw_top : draw_bot)) continue; // this disk not drawn
		norm_comp const normal((bool(bt) ^ inv_tb) ? v12 : -v12);
		unsigned const center_ix(itix);
		itri_verts[itix++].assign(ce[bt], normal, half_end_tscale, half_end_tscale, cw); // center

		for (unsigned I = 0; I < ndiv; ++I) {
			unsigned const i(bt ? ndiv-I-1 : I); // invert winding order for top face
			vector3d const &side_normal(vpn.n[i]);
			itri_verts[itix++].assign(vpn.p[(i<<1) + bt], normal, half_end_tscale*(side_normal.x + 1.0), half_end_tscale*(side_normal.y + 1.0), cw); // assign tcs from side normal
			indices[iix++] = center_ix; // center
			indices[iix++] = center_ix + i + 1;
			indices[iix++] = center_ix + ((i+1)%ndiv) + 1;
		}
	} // for bt
	if (inv_tb) {std::reverse(indices.begin()+ixs_start, indices.end());} // reverse the order to swap triangle winding order
	if (two_sided) {add_inverted_triangles(itri_verts, indices, itris_start, ixs_start);}
}

void rgeom_mat_t::add_disk_to_verts(point const &pos, float radius, bool normal_z_neg, colorRGBA const &color) {
	assert(radius > 0.0);
	color_wrapper const cw(color);
	norm_comp const nc(normal_z_neg ? -plus_z : plus_z);
	unsigned const ndiv(N_CYL_SIDES), itris_start(itri_verts.size());
	float const css(-1.0*TWO_PI/(float)ndiv), sin_ds(sin(css)), cos_ds(cos(css));
	float sin_s(0.0), cos_s(1.0);
	itri_verts.emplace_back(pos, nc, 0.5, 0.5, cw);

	for (unsigned i = 0; i < ndiv; ++i) {
		float const s(sin_s), c(cos_s);
		itri_verts.emplace_back((pos + point(radius*s, radius*c, 0.0)), nc, 0.5*(1.0 + s), 0.5*(1.0 + c), cw);
		indices.push_back(itris_start); // center
		indices.push_back(itris_start + i + 1);
		indices.push_back(itris_start + ((i+1)%ndiv) + 1);
		sin_s = s*cos_ds + c*sin_ds;
		cos_s = c*cos_ds - s*sin_ds;
	}
}

// Note: size can be nonuniform in X/Y/Z
void rgeom_mat_t::add_sphere_to_verts(point const &center, vector3d const &size, colorRGBA const &color, bool low_detail,
	vector3d const &skip_hemi_dir, xform_matrix const *const matrix)
{
	static vector<vert_norm_tc> cached_verts[2]; // high/low detail, reused across all calls
	static vector<vert_norm_comp_tc> cached_vncs[2];
	static vector<unsigned> cached_ixs[2];
	vector<vert_norm_tc> &verts(cached_verts[low_detail]);
	vector<vert_norm_comp_tc> &vncs(cached_vncs[low_detail]);
	vector<unsigned> &ixs(cached_ixs[low_detail]);
	bool const draw_hemisphere(skip_hemi_dir != zero_vector);

	if (verts.empty()) { // not yet created, create and cache verts
		sd_sphere_d sd(all_zeros, 1.0, (low_detail ? N_SPHERE_DIV/2 : N_SPHERE_DIV));
		sphere_point_norm spn;
		sd.gen_points_norms(spn);
		sd.get_quad_points(verts, &ixs);
		assert((ixs.size()&3) == 0); // must be a multiple of 4
		vncs.resize(verts.size());
		for (unsigned i = 0; i < verts.size(); ++i) {vncs[i] = vert_norm_comp_tc(verts[i].v, verts[i].n, verts[i].t[0], verts[i].t[1]);}
	}
	color_wrapper const cw(color);
	unsigned const ioff(itri_verts.size());

	if (matrix) { // must apply matrix transform to verts and normals and reconstruct norm_comps
		for (auto i = verts.begin(); i != verts.end(); ++i) {
			point pt(i->v*size);
			vector3d normal(i->n);
			matrix->apply_to_vector3d(pt); matrix->apply_to_vector3d(normal);
			itri_verts.emplace_back((pt + center), normal, i->t[0], i->t[1], cw);
		}
	}
	else { // can use vncs (norm_comps)
		for (auto i = vncs.begin(); i != vncs.end(); ++i) {itri_verts.emplace_back((i->v*size + center), *i, i->t[0], i->t[1], cw);}
	}
	for (auto i = ixs.begin(); i != ixs.end(); i += 4) { // indices are for quads, but we want triangles, so expand them
		if (draw_hemisphere) { // only draw one hemisphere; can drop some verts as well, but that's difficult
			vector3d const face_normal(verts[*(i+0)].n + verts[*(i+1)].n + verts[*(i+2)].n); // use one triangle, no need to normalize
			if (dot_product(face_normal, skip_hemi_dir) > 0.0) continue; // skip this face/quad (optimization)
		}
		indices.push_back(*(i+0) + ioff); indices.push_back(*(i+1) + ioff); indices.push_back(*(i+2) + ioff);
		indices.push_back(*(i+3) + ioff); indices.push_back(*(i+0) + ioff); indices.push_back(*(i+2) + ioff);
	} // for i
	assert(indices.back() < itri_verts.size());
}

void rgeom_mat_t::add_triangle_to_verts(point const v[3], colorRGBA const &color, bool two_sided, float tscale) {
	color_wrapper cw(color);
	norm_comp normal(get_poly_norm(v));
	float const ts[3] = {0.0, 0.0, tscale}, tt[3] = {0.0, tscale, 0.0}; // hard-coded for now, maybe pass in?

	for (unsigned side = 0; side < 2; ++side) {
		for (unsigned n = 0; n < 3; ++n) {
			indices.push_back(itri_verts.size()); // since we only support indexed triangles, we have to assign each vertex its own index
			unsigned const ix(side ? 2-n : n); // reverse order for side=1
			itri_verts.emplace_back(v[ix], normal, ts[ix], tt[ix], cw);
		}
		if (side == 0) {
			if (!two_sided) break; // single sided
			normal.invert_normal();
		}
	} // for side
}

class rgeom_alloc_t {
	deque<rgeom_storage_t> free_list; // one per unique texture ID/material
public:
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
		free_list.push_back(rgeom_storage_t(s.tex)); // record tex of incoming element
		s.swap_vectors(free_list.back()); // transfer existing capacity to free list; clear capacity from s
	}
	unsigned get_mem_usage() const {
		unsigned mem(free_list.size()*sizeof(rgeom_storage_t));
		for (auto i = free_list.begin(); i != free_list.end(); ++i) {mem += i->get_mem_usage();}
		return mem;
	}
	unsigned size() const {return free_list.size();}
};

rgeom_alloc_t rgeom_alloc; // static allocator with free list, shared across all buildings; not thread safe

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
	vao_mgr.clear_vbos();
	clear_vectors();
	num_verts = num_ixs = 0;
}

void rotate_verts(vector<rgeom_mat_t::vertex_t> &verts, building_t const &building) {
	point const center(building.bcube.get_cube_center());

	for (auto i = verts.begin(); i != verts.end(); ++i) {
		building.do_xy_rotate(center, i->v);
		vector3d n(i->get_norm());
		building.do_xy_rotate_normal(n);
		i->set_norm(n);
	}
}

void rgeom_mat_t::create_vbo(building_t const &building) {
	if (building.is_rotated()) { // rotate all vertices to match the building rotation
		rotate_verts(quad_verts, building);
		rotate_verts(itri_verts, building);
	}
	create_vbo_inner();
	rgeom_alloc.free(*this); // vertex and index data is no longer needed and can be cleared
}
void rgeom_mat_t::create_vbo_inner() {
	assert(itri_verts.empty() == indices.empty());
	unsigned qsz(quad_verts.size()*sizeof(vertex_t)), itsz(itri_verts.size()*sizeof(vertex_t));
	num_verts   = quad_verts.size() + itri_verts.size();
	assert(num_verts > 0); // too strong?
	vao_mgr.vbo = ::create_vbo();
	check_bind_vbo(vao_mgr.vbo);
	upload_vbo_data(nullptr, num_verts*sizeof(vertex_t));
	if (itsz > 0) {upload_vbo_sub_data(itri_verts.data(), 0,    itsz);}
	if (qsz  > 0) {upload_vbo_sub_data(quad_verts.data(), itsz, qsz );}
	bind_vbo(0);
	gen_quad_ixs(indices, 6*(quad_verts.size()/4), itri_verts.size()); // append indices for quad_verts
	create_vbo_and_upload(vao_mgr.ivbo, indices, 1, 1); // indices should always be nonempty
	num_ixs = indices.size();

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
}

void brg_batch_draw_t::set_camera_dir_mask(point const &camera_bs, cube_t const &bcube) {
	camera_dir_mask = 0;
	for (unsigned d = 0; d < 3; ++d) {
		if (camera_bs[d] < bcube.d[d][1]) {camera_dir_mask |=  1<<(2*d);}
		if (camera_bs[d] > bcube.d[d][0]) {camera_dir_mask |= 1<<(2*d+1);}
	}
}
void brg_batch_draw_t::add_material(rgeom_mat_t const &m) {
#if 0 // doesn't seem like this helps for ~100 materials, but it may help when more materials are added
	unsigned const tid((m.tex.tid < 0) ? WHITE_TEX : m.tex.tid);
	if (tid >= tid_to_first_mat_map.size()) {tid_to_first_mat_map.resize(tid+1, -1);}
	int &tex_ix(tid_to_first_mat_map[tid]);
	unsigned const start_ix((tex_ix >= 0) ? tex_ix : 0); // start searching at the first material using this texture
	assert(start_ix <= to_draw.size());
		
	for (auto i = to_draw.begin() + start_ix; i != to_draw.end(); ++i) {
		if (i->tex.is_compat_ignore_shadowed(m.tex)) {i->mats.push_back(&m); return;} // found existing material
	}
	if (tex_ix < 0) {tex_ix = to_draw.size();} // cache this material index for later calls
#else
	for (auto &i : to_draw) { // check all existing materials for a matching texture, etc.
		if (i.tex.is_compat_ignore_shadowed(m.tex)) {i.mats.push_back(&m); return;} // found existing material
	}
#endif
	to_draw.emplace_back(m); // add a new material entry
}
void brg_batch_draw_t::draw_and_clear(shader_t &s) {
	tid_nm_pair_dstate_t state(s);
	enable_blend(); // needed for rugs and book text

	for (auto &i : to_draw) {
		if (i.mats.empty()) continue; // empty slot
		i.tex.set_gl(state);
		for (auto const &m : i.mats) {m->draw_inner(0);} // shadow_only=0
		i.tex.unset_gl(state);
		i.mats.clear(); // clear mats but not to_draw
	}
	disable_blend();
	indexed_vao_manager_with_shadow_t::post_render();
}

// shadow_only: 0=non-shadow pass, 1=shadow pass, 2=shadow pass with alpha mask texture
void rgeom_mat_t::draw(tid_nm_pair_dstate_t &state, brg_batch_draw_t *bbd, int shadow_only, bool reflection_pass) {
	if (shadow_only && !en_shadows)  return; // shadows not enabled for this material (picture, whiteboard, rug, etc.)
	if (shadow_only && tex.emissive) return; // assume this is a light source and shouldn't produce shadows (also applies to bathroom windows, which don't produce shadows)
	if (reflection_pass && tex.tid == REFLECTION_TEXTURE_ID) return; // don't draw reflections of mirrors as this doesn't work correctly
	assert(num_verts > 0); // too strong? should be okay to remove this check
	if (num_verts == 0) return;
	vao_setup(shadow_only);

	// Note: the shadow pass doesn't normally bind textures and set uniforms, so we don't need to combine those calls into batches
	if (bbd != nullptr && !shadow_only) { // add to batch draw (optimization)
		if (dir_mask > 0 && bbd->camera_dir_mask > 0 && (dir_mask & bbd->camera_dir_mask) == 0) return; // check for visible surfaces
		bbd->add_material(*this);
	}
	else { // draw this material now
		if (shadow_only != 1) {tex.set_gl  (state);} // ignores texture scale for now; enable alpha texture for shadow pass
		draw_inner(shadow_only);
		if (shadow_only != 1) {tex.unset_gl(state);}
	}
}
void rgeom_mat_t::draw_inner(int shadow_only) const {
	vao_mgr.pre_render(shadow_only != 0);
	glDrawRangeElements(GL_TRIANGLES, 0, num_verts, num_ixs, GL_UNSIGNED_INT, nullptr);
}
void rgeom_mat_t::vao_setup(bool shadow_only) {
	vao_mgr.create_and_upload(vector<vertex_t>(), vector<unsigned>(), shadow_only, 0, 1); // pass empty vectors because data is already uploaded; dynamic_level=0, setup_pointers=1
}
void rgeom_mat_t::upload_draw_and_clear(tid_nm_pair_dstate_t &state) { // Note: called by draw_interactive_player_obj() and water_draw_t
	if (empty()) return; // nothing to do; can this happen?
	create_vbo_inner();
	draw(state, nullptr, 0, 0); // no brg_batch_draw_t
	clear();
}

void building_materials_t::clear() {
	for (iterator m = begin(); m != end(); ++m) {m->clear();}
	vector<rgeom_mat_t>::clear();
}
unsigned building_materials_t::count_all_verts() const {
	unsigned num_verts(0);
	for (const_iterator m = begin(); m != end(); ++m) {num_verts += m->num_verts;}
	return num_verts;
}
rgeom_mat_t &building_materials_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows) {
	// for now we do a simple linear search because there shouldn't be too many unique materials
	for (iterator m = begin(); m != end(); ++m) {
		if (!m->tex.is_compatible(tex)) continue;
		if (inc_shadows) {m->enable_shadows();}
		// tscale diffs don't make new materials; copy tscales from incoming tex; this field may be used locally by the caller, but isn't used for drawing
		m->tex.tscale_x = tex.tscale_x; m->tex.tscale_y = tex.tscale_y;
		return *m;
	}
	emplace_back(tex); // not found, add a new material
	if (inc_shadows) {back().enable_shadows();}
	rgeom_alloc.alloc(back());
	return back();
}
void building_materials_t::create_vbos(building_t const &building) {
	for (iterator m = begin(); m != end(); ++m) {m->create_vbo(building);}
}
void building_materials_t::draw(brg_batch_draw_t *bbd, shader_t &s, int shadow_only, bool reflection_pass) {
	//highres_timer_t timer("Draw Materials"); // 0.0168
	static vector<iterator> text_mats;
	text_mats.clear();
	tid_nm_pair_dstate_t state(s);

	// first pass, draw regular materials (excluding text)
	for (iterator m = begin(); m != end(); ++m) {
		if (m->tex.tid == FONT_TEXTURE_ID) {text_mats.push_back(m);} // skip in this pass
		else {m->draw(state, bbd, shadow_only, reflection_pass);}
	}
	// second pass, draw text (if it exists) so that alpha blending works; really only needed for the building the player is in
	for (auto m = text_mats.begin(); m != text_mats.end(); ++m) {(*m)->draw(state, bbd, shadow_only, reflection_pass);}
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
	has_elevators = 0;
}
void building_room_geom_t::clear_materials() { // clears all materials
	clear_lit_materials();
	mats_lights.clear();
	mats_detail.clear();
}
void building_room_geom_t::clear_lit_materials() { // clears materials that have color based on room lighting
	// Note: mats_detail and mats_lights are excluded
	clear_static_vbos();
	clear_static_small_vbos();
	mats_dynamic.clear();
	mats_doors  .clear();
}
void building_room_geom_t::clear_static_vbos() { // used to clear pictures
	mats_static.clear();
	obj_model_insts.clear(); // these are associated with static VBOs
	mats_alpha.clear();
}
void building_room_geom_t::clear_static_small_vbos() {
	mats_small.clear();
	mats_amask.clear();
}

rgeom_mat_t &building_room_geom_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows, bool dynamic, unsigned small, bool transparent) {
	// small: 0=mats_static, 1=mats_small, 2=mats_detail
	return (dynamic ? mats_dynamic : (small ? ((small == 2) ? mats_detail : mats_small) : (transparent ? mats_alpha : mats_static))).get_material(tex, inc_shadows);
}
rgeom_mat_t &building_room_geom_t::get_metal_material(bool inc_shadows, bool dynamic, unsigned small) {
	tid_nm_pair_t tex(-1, 1.0, inc_shadows);
	tex.set_specular(0.8, 60.0);
	return get_material(tex, inc_shadows, dynamic, small);
}

void room_object_t::set_as_bottle(unsigned rand_id, unsigned max_type, bool no_empty) {
	assert(max_type > 0 && max_type < NUM_BOTTLE_TYPES);
	obj_id = (uint16_t)rand_id;
	while (get_bottle_type() > max_type) {obj_id += 13;} // cycle with a prime number until a valid type is selected
	if (no_empty) {obj_id &= 127;} // strip off second empty bit
	color  = bottle_params[get_bottle_type()].color;
}

void building_room_geom_t::create_static_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom"); // 2.35ms
	float const tscale(2.0/obj_scale);
	mats_static.clear();
	mats_alpha .clear();
	tid_nm_pair_t const &wall_tex(building.get_material().wall_tex);

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible() || i->is_dynamic()) continue; // skip invisible and dynamic objects
		assert(i->is_strictly_normalized());
		assert(i->type < NUM_ROBJ_TYPES);

		switch (i->type) {
		case TYPE_TABLE:   add_table   (*i, tscale, 0.12, 0.08); break; // top_dz=12% of height, leg_width=8% of height
		case TYPE_CHAIR:   add_chair   (*i, tscale); break;
		case TYPE_STAIR:   add_stair   (*i, tscale, tex_origin); break;
		case TYPE_STAIR_WALL: add_stairs_wall(*i, tex_origin, wall_tex); break;
		case TYPE_RUG:     add_rug     (*i); break;
		case TYPE_PICTURE: add_picture (*i); break;
		case TYPE_WBOARD:  add_picture (*i); break;
		case TYPE_BOOK:    add_book    (*i, 1, 0); break;
		case TYPE_BCASE:   add_bookcase(*i, 1, 0, tscale, 0); break;
		case TYPE_WINE_RACK: add_wine_rack(*i, 1, 0, tscale); break;
		case TYPE_DESK:    add_desk    (*i, tscale, 1, 0); break;
		case TYPE_RDESK:   add_reception_desk(*i, tscale); break;
		case TYPE_BED:     add_bed     (*i, 1, 0, tscale); break;
		case TYPE_WINDOW:  add_window  (*i, tscale); break;
		case TYPE_TUB:     add_tub_outer(*i); break;
		case TYPE_TV: case TYPE_MONITOR: add_tv_picture(*i); break;
		case TYPE_CUBICLE: add_cubicle (*i, tscale); break;
		case TYPE_STALL:   add_br_stall(*i); break;
		case TYPE_SIGN:    add_sign    (*i, 1, 0); break;
		case TYPE_COUNTER: add_counter (*i, tscale); break;
		case TYPE_KSINK:   add_counter (*i, tscale); break; // counter with kitchen sink
		case TYPE_BRSINK:  add_counter (*i, tscale); break; // counter with bathroom sink
		case TYPE_CABINET: add_cabinet (*i, tscale); break;
		case TYPE_PLANT:   add_potted_plant(*i, 1, 0); break; // pot only
		case TYPE_DRESSER: case TYPE_NIGHTSTAND: add_dresser(*i, tscale, 1, 0); break;
		case TYPE_FLOORING:add_flooring(*i, tscale); break;
		case TYPE_CLOSET:  add_closet  (*i, wall_tex, 1, 0); break;
		case TYPE_MIRROR:  add_mirror  (*i); break;
		case TYPE_SHOWER:  add_shower  (*i, tscale); break;
		case TYPE_MWAVE:   add_mwave   (*i); break;
		case TYPE_BLINDS:  add_blinds  (*i); break;
		case TYPE_FPLACE:  add_fireplace(*i, tscale); break;
		case TYPE_ELEVATOR: break; // not handled here
		case TYPE_BLOCKER:  break; // not drawn
		case TYPE_COLLIDER: break; // not drawn
		default: break;
		} // end switch
	} // for i
	// Note: verts are temporary, but cubes are needed for things such as collision detection with the player and ray queries for indir lighting
	//timer_t timer2("Create VBOs"); // < 2ms
	mats_static.create_vbos(building);
	mats_alpha .create_vbos(building);
	//cout << "static: size: " << rgeom_alloc.size() << " mem: " << rgeom_alloc.get_mem_usage() << endl; // start=78MB, peak=193MB
}

void building_room_geom_t::create_small_static_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Small"); // 7.8ms, slow building at 26,16
	mats_small.clear();
	mats_amask.clear();
	model_objs.clear(); // currently model_objs are only created for small objects in drawers, so we clear this here
	add_small_static_objs_to_verts(expanded_objs);
	add_small_static_objs_to_verts(objs);
	mats_small.create_vbos(building);
	mats_amask.create_vbos(building);
}

void building_room_geom_t::add_small_static_objs_to_verts(vect_room_object_t const &objs_to_add) {
	if (objs_to_add.empty()) return; // don't add untextured material, otherwise we may fail the (num_verts > 0) assert
	float const tscale(2.0/obj_scale);
	// Note: we no longer have to add an untextured material here for book covers/text because the mats_detail should add that first

	for (auto i = objs_to_add.begin(); i != objs_to_add.end(); ++i) {
		if (!i->is_visible() || i->is_dynamic()) continue; // skip invisible and dynamic objects
		assert(i->is_strictly_normalized());
		assert(i->type < NUM_ROBJ_TYPES);

		switch (i->type) {
		case TYPE_BOOK:      add_book     (*i, 0, 1); break;
		case TYPE_BCASE:     add_bookcase (*i, 0, 1, tscale, 0); break;
		case TYPE_BED:       add_bed      (*i, 0, 1, tscale); break;
		case TYPE_DESK:      add_desk     (*i, tscale, 0, 1); break;
		case TYPE_DRESSER: case TYPE_NIGHTSTAND: add_dresser(*i, tscale, 0, 1); break;
		case TYPE_TCAN:      add_trashcan (*i); break;
		case TYPE_SIGN:      add_sign     (*i, 0, 1); break;
		case TYPE_CLOSET:    add_closet   (*i, tid_nm_pair_t(), 0, 1); break; // add closet wall trim and interior objects, don't need wall_tex
		case TYPE_RAILING:   add_railing  (*i); break;
		case TYPE_PLANT:     add_potted_plant(*i, 0, 1); break; // plant only
		case TYPE_CRATE:     add_crate    (*i); break; // not small but only added to windowless rooms
		case TYPE_BOX:       add_box      (*i); break; // not small but only added to windowless rooms
		case TYPE_SHELVES:   add_shelves  (*i, tscale); break; // not small but only added to windowless rooms
		case TYPE_COMPUTER:  add_computer (*i); break;
		case TYPE_KEYBOARD:  add_keyboard (*i); break;
		case TYPE_WINE_RACK: add_wine_rack(*i, 0, 1, tscale); break;
		case TYPE_BOTTLE:    add_bottle   (*i); break;
		case TYPE_PAPER:     add_paper    (*i); break;
		case TYPE_PAINTCAN:  add_paint_can(*i); break;
		case TYPE_PEN: case TYPE_PENCIL: case TYPE_MARKER: add_pen_pencil_marker(*i); break;
		case TYPE_LG_BALL:   add_lg_ball  (*i); break;
		case TYPE_HANGER_ROD:add_hanger_rod(*i); break;
		case TYPE_DRAIN:     add_drain_pipe(*i); break;
		case TYPE_KEY:       if (has_key_3d_model()) {model_objs.push_back(*i);} else {add_key(*i);} break; // draw or add as 3D model
		case TYPE_MONEY:     add_money (*i); break;
		case TYPE_PHONE:     add_phone (*i); break;
		case TYPE_TPROLL:    add_tproll(*i); break;
		case TYPE_TAPE:      add_tape  (*i); break;
		case TYPE_SPRAYCAN:  add_spraycan(*i); break;
		case TYPE_CRACK:     add_crack (*i); break;
		case TYPE_SWITCH:    add_switch(*i, 0); break; // draw_detail_pass=0
		case TYPE_PLATE:     add_plate (*i); break;
		case TYPE_LAPTOP:    add_laptop(*i); break;
		case TYPE_BUTTON:    if (!(i->flags & RO_FLAG_IN_ELEV)) {add_button(*i);} break; // skip buttons inside elevators, which are drawn as dynamic objects
		case TYPE_LBASKET:   add_laundry_basket(*i); break;
		case TYPE_WHEATER:   add_water_heater  (*i); break; // small since this object is only added to basements
		case TYPE_TOASTER:   add_toaster_proxy(*i);  break;
		default: break;
		} // end switch
	} // for i
}

void building_room_geom_t::create_detail_vbos(building_t const &building) {
	mats_detail.clear();
	// currently only small objects that are non-interactive and can't be taken; TYPE_SWITCH almost counts
	auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators
	tid_nm_pair_t const &wall_tex(building.get_material().wall_tex);

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_visible()) continue;

		switch (i->type) {
		case TYPE_OUTLET:     add_outlet(*i); break;
		case TYPE_SWITCH:     add_switch(*i, 1); break; // draw_detail_pass=0
		case TYPE_PG_WALL:    add_parking_garage_wall(*i, tex_origin, wall_tex); break;
		case TYPE_PARK_SPACE: add_parking_space(*i, tex_origin, wall_tex.tscale_x); break;
		case TYPE_RAMP:       add_pg_ramp(*i, tex_origin, wall_tex.tscale_x); break;
		case TYPE_PIPE:       add_pipe(*i); break;
		case TYPE_CURB:       add_curb(*i); break;
		default: break;
		} // end switch
	} // for i
	for (auto const &i : trim_objs) {
		assert(i.type == TYPE_WALL_TRIM);
		add_wall_trim(i);
	}
	mats_detail.create_vbos(building);
}

void building_room_geom_t::create_obj_model_insts(building_t const &building) { // handle drawing of 3D models
	//highres_timer_t timer("Gen Room Model Insts");
	obj_model_insts.clear();

	for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
		auto const &obj_vect((vect_id == 1) ? expanded_objs : objs);
		unsigned const obj_id_offset((vect_id == 1) ? objs.size() : 0);
		auto objs_end((vect_id == 1) ? expanded_objs.end() : get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = obj_vect.begin(); i != objs_end; ++i) {
			if (!i->is_visible() || !i->is_obj_model_type()) continue;
			vector3d dir(i->get_dir());

			if (i->flags & RO_FLAG_RAND_ROT) {
				float const angle(123.4*i->x1() + 456.7*i->y1() + 567.8*i->z1()); // random rotation angle based on position
				vector3d const rand_dir(vector3d(sin(angle), cos(angle), 0.0).get_norm());
				dir = ((dot_product(rand_dir, dir) < 0.0) ? -rand_dir : rand_dir); // random, but facing in the correct general direction
			}
			if (building.is_rotated()) {building.do_xy_rotate_normal(dir);}
			obj_model_insts.emplace_back((i - obj_vect.begin() + obj_id_offset), dir);
			//get_untextured_material().add_cube_to_verts_untextured(*i, WHITE); // for debugging of model bcubes
		} // for i
	} // for vect_id
}

void building_room_geom_t::create_lights_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Light"); // 0.75ms
	float const tscale(2.0/obj_scale);
	auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->is_visible() && i->type == TYPE_LIGHT) {add_light(*i, tscale);}
	}
	mats_lights.create_vbos(building);
}

void building_room_geom_t::create_dynamic_vbos(building_t const &building) {
	//highres_timer_t timer(string("Gen Room Geom Dynamic ") + (building.is_house ? "house" : "office"));
	
	if (!obj_dstate.empty()) { // we have an object with dynamic state
		auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (!i->is_dynamic() || !i->is_visible()) continue; // only visible + dynamic objects; can't do VFC because this is not updated every frame
			switch (i->type) {
			case TYPE_LG_BALL: add_lg_ball(*i); break;
			default: assert(0); // not a supported dynamic object type
			}
		} // for i
	}
	for (auto e = building.interior->elevators.begin(); e != building.interior->elevators.end(); ++e) {
		assert(e->car_obj_id < objs.size());
		assert(e->button_id_start < e->button_id_end && e->button_id_end <= objs.size());
		float const fc_thick_scale(building.get_elevator_fc_thick_scale());
		add_elevator_doors(*e, fc_thick_scale); // add dynamic elevator doors
		add_elevator(objs[e->car_obj_id], 2.0/obj_scale, fc_thick_scale, building.calc_floor_offset(e->z1()), building.has_parking_garage); // draw elevator car for this elevator

		for (auto j = objs.begin() + e->button_id_start; j != objs.begin() + e->button_id_end; ++j) {
			if (j->type == TYPE_BLOCKER) continue; // button was removed?
			assert(j->type == TYPE_BUTTON);
			if (j->flags & RO_FLAG_IN_ELEV) {add_button(*j);} // add button as a dynamic object if it's inside the elevator
		}
	} // for e
	mats_dynamic.create_vbos(building);
}

void building_room_geom_t::create_door_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Doors"); // 0.1ms
	vector<door_t> const &doors(building.interior->doors);
	uint8_t const door_type(building.is_house ? (uint8_t)tquad_with_ix_t::TYPE_HDOOR : (uint8_t)tquad_with_ix_t::TYPE_ODOOR);

	for (auto i = doors.begin(); i != doors.end(); ++i) {
		building.add_door_verts(*i, *this, door_type, i->dim, i->open_dir, i->open, 0, 0, i->on_stairs, i->hinge_side); // opens_out=0, exterior=0
	}
	mats_doors.create_vbos(building);
}

void rotate_dir_about_z(vector3d &dir, float angle) { // Note: assumes dir is normalized
	if (angle == 0.0) return;
	assert(dir.z == 0.0); // dir must be in XY plane
	float const new_angle(atan2(dir.y, dir.x) + angle);
	dir.assign(cosf(new_angle), sinf(new_angle), 0.0);
}
void apply_room_obj_rotate(room_object_t &obj, obj_model_inst_t &inst) {
	if (!(obj.flags & RO_FLAG_ROTATING)) return;

	if (obj.type == TYPE_OFF_CHAIR) {
		if (office_chair_rot_rate == 0.0) {obj.flags &= ~RO_FLAG_ROTATING; return;} // if no longer rotating, clear rotation bit
		rotate_dir_about_z(inst.dir, office_chair_rot_rate*fticks);
	}
	else if (obj.type == TYPE_HANGER || obj.type == TYPE_CLOTHES) {
		inst.dir = obj.get_dir(); // reset before applying rotate
		float const angle(((obj.flags & RO_FLAG_ADJ_LO) ? -1.0 : 1.0)*0.08*TWO_PI);
		rotate_dir_about_z(inst.dir, -angle); // limited rotation angle
	}
	else {
		assert(0); // unsupported object type
	}
}

/*static*/ void building_room_geom_t::draw_interactive_player_obj(carried_item_t const &c, shader_t &s, vector3d const &xlate) {
	static rgeom_mat_t mat; // allocated memory is reused across frames; VBO is recreated every time
	bool needs_blend(0);

	if (c.type == TYPE_SPRAYCAN || c.type == TYPE_MARKER) {
		room_object_t c_rot(c);
		c_rot.dir = 0; // facing up
		unsigned const dim(get_max_dim(c.get_size()));

		if (dim != 2) { // if not oriented in Z
			UNROLL_2X(swap(c_rot.d[dim][i_], c_rot.d[2][i_]);); // rotate into Z dir
			c_rot.translate(c.get_cube_center() - c_rot.get_cube_center()); // translate it back to the correct location
		}
		if (c.type == TYPE_SPRAYCAN) {add_spraycan_to_material(c_rot, mat);}
		else {add_pen_pencil_marker_to_material(c_rot, mat);}
	}
	else if (c.type == TYPE_TPROLL || c.type == TYPE_TAPE) { // apply get_player_cview_rot_matrix()?
		add_vert_roll_to_material(c, mat, c.get_remaining_capacity_ratio(), 1); // player_held=1
		needs_blend = 1;
	}
	else if (c.type == TYPE_BOOK) {
		static building_room_geom_t tmp_rgeom;
		float const z_rot_angle(-(atan2(cview_dir.y, cview_dir.x) + PI_TWO));
		tmp_rgeom.add_book(c, 1, 1, 0.0, 0, 0, z_rot_angle);
		enable_blend(); // needed for book text
		tmp_rgeom.mats_small.upload_draw_and_clear(s);
		disable_blend();
		return;
	}
	else if (c.type == TYPE_PHONE) {
		float const z_rot_angle(-atan2(cview_dir.y, cview_dir.x));

		if (c.flags & (RO_FLAG_EMISSIVE | RO_FLAG_OPEN)) { // phone is ringing or locked screen
			static rgeom_mat_t screen_mat;
			screen_mat.tex = get_phone_tex(c);
			screen_mat.add_cube_to_verts(c, WHITE, all_zeros, ~EF_Z2, 0, 1); // mirror_x=1
			rotate_verts(screen_mat.quad_verts, plus_z, z_rot_angle, c.get_cube_center(), 0); // rotate all quad verts about Z axis
			tid_nm_pair_dstate_t state(s);
			screen_mat.upload_draw_and_clear(state);
		}
		else {mat.add_cube_to_verts(c, BLACK, all_zeros, ~EF_Z2);} // screen drawn as black
		mat.add_cube_to_verts(c, c.color, all_zeros, EF_Z2);
		rotate_verts(mat.quad_verts, plus_z, z_rot_angle, c.get_cube_center(), 0); // rotate all quad verts about Z axis
	}
	else if (c.type == TYPE_RAT) { // draw the rat facing the player
		building_obj_model_loader.draw_model(s, c.get_cube_center(), c, cview_dir, rat_color, xlate, OBJ_MODEL_RAT, 0); // facing away from the player; shadow_pass=0
		check_mvm_update();
		return; // don't need to run the code below
	}
	else {assert(0);}
	if (needs_blend) {enable_blend();}
	tid_nm_pair_dstate_t state(s);
	mat.upload_draw_and_clear(state);
	if (needs_blend) {disable_blend();}
}

class water_draw_t {
	rgeom_mat_t mat;
	float tex_off;
public:
	water_draw_t() : mat(rgeom_mat_t(tid_nm_pair_t(FOAM_TEX))), tex_off(0.0) {}

	void add_water_for_sink(room_object_t const &obj) {
		if (!obj.is_active()) return; // not turned on
		bool const is_cube(obj.type == TYPE_KSINK || obj.type == TYPE_BRSINK);
		float const dz(obj.dz());
		cube_t c;
		point pos(obj.get_cube_center());
		pos[obj.dim] += (obj.dir ? -1.0 : 1.0)*(is_cube ? 0.25 : 0.095)*obj.get_length(); // move toward the back of the sink
		c.set_from_sphere(pos, (is_cube ? 0.02 : 0.0055)*dz);
		set_cube_zvals(c, (obj.z1() + (is_cube ? 0.7 : 0.6)*dz), (obj.z1() + (is_cube ? 1.3 : 0.925)*dz));
		unsigned const verts_start(mat.itri_verts.size());
		mat.add_vcylin_to_verts(c, colorRGBA(WHITE, 0.5), 0, 0, 0, 0, 1.0, 1.0, 0.2);
		for (auto i = mat.itri_verts.begin() + verts_start; i != mat.itri_verts.end(); ++i) {i->t[1] *= 1.2; i->t[1] += tex_off;}
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

int room_object_t::get_model_id() const {
	assert(type >= TYPE_TOILET);
	if (type == TYPE_MONITOR) return OBJ_MODEL_TV; // monitor has same model as TV
	int id((int)type + OBJ_MODEL_TOILET - TYPE_TOILET);
	if (type == TYPE_HANGER || type == TYPE_CLOTHES) {id += ((int)item_flags << 8);} // choose a sub_model_id for these types using bits 8-15
	return id;
}

void building_t::draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, occlusion_checker_noncity_t &oc, vector3d const &xlate, unsigned building_ix,
	bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building)
{
	if (!interior || !interior->room_geom) return;
	if (ENABLE_MIRROR_REFLECTIONS && !shadow_only && !reflection_pass && player_in_building) {find_mirror_needing_reflection(xlate);}
	interior->room_geom->draw(bbd, s, *this, oc, xlate, building_ix, shadow_only, reflection_pass, inc_small, player_in_building);
}
void building_t::gen_and_draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, occlusion_checker_noncity_t &oc, vector3d const &xlate, vect_cube_t &ped_bcubes,
	unsigned building_ix, int ped_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building)
{
	if (!interior) return;
	if (!global_building_params.enable_rotated_room_geom && is_rotated()) return; // rotated buildings: need to fix texture coords, room object collisions, mirrors, etc.

	if (!shadow_only && !camera_pdu.point_visible_test(bcube.get_cube_center() + xlate)) {
		// skip if none of the building parts are visible to the camera; this is rare, so it may not help
		bool any_part_visible(0);

		for (auto p = parts.begin(); p != parts.end(); ++p) {
			if (camera_pdu.cube_visible(*p + xlate)) {any_part_visible = 1; break;}
		}
		if (!any_part_visible) return;
	}
	if (!has_room_geom()) {
		rand_gen_t rgen;
		rgen.set_state(building_ix, parts.size()); // set to something canonical per building
		ped_bcubes.clear();
		if (ped_ix >= 0) {get_ped_bcubes_for_building(ped_ix, ped_bcubes);}
		gen_room_details(rgen, ped_bcubes, building_ix); // generate so that we can draw it
		assert(has_room_geom());
	}
	if (has_room_geom() && inc_small == 2) {add_wall_and_door_trim_if_needed();} // gen trim when close to the player
	draw_room_geom(bbd, s, oc, xlate, building_ix, shadow_only, reflection_pass, inc_small, player_in_building);
}
void building_t::clear_room_geom(bool force) {
	if (!has_room_geom()) return;
	if (interior->room_geom->modified_by_player) return; // keep the player's modifications and don't delete the room geom
	interior->room_geom->clear(); // free VBO data before deleting the room_geom object
	interior->room_geom.reset();
}

void draw_stove_flames(room_object_t const &stove, point const &camera_bs, shader_t &s) {
	if (stove.item_flags == 0) return; // no burners on
	static quad_batch_draw flame_qbd; // reused across frames
	vector3d const sz(stove.get_size());
	bool const dim(stove.dim), dir(stove.dir);
	float const zval(stove.z2() - 0.23*sz.z), dsign(dir ? -1.0 : 1.0);
	
	for (unsigned w = 0; w < 2; ++w) { // width dim
		float const radius((w ? 0.09 : 0.07)*sz.z), wval(stove.d[!dim][0] + ((bool(w) ^ dim ^ dir ^ 1) ? 0.72 : 0.28)*sz[!dim]); // left/right burners

		for (unsigned d = 0; d < 2; ++d) { // depth dim
			if (!(stove.item_flags & (1U<<(2U*w + d)))) continue; // burner not on
			float const dval(stove.d[dim][dir] + dsign*(d ? 0.66 : 0.375)*sz[dim]); // front/back burners
			point pos(dval, wval, zval);
			if (dim) {swap(pos.x, pos.y);}
			flame_qbd.add_quad_dirs(pos, dsign*radius*plus_x, -dsign*radius*plus_y, WHITE); // use a negative Y to get the proper CW order; flip with dsign for symmetry
		} // for d
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

class spider_draw_t {
	rgeom_mat_t mat;
	unsigned cur_vert_pos;
	bool is_setup;

	void add_eye(point const &pos, float radius) {
		cube_t eye_bc(pos);
		eye_bc.expand_by(radius);
		mat.add_sphere_to_verts(eye_bc, DK_RED); // eye
	}
	void init() {
		// generate spider geometry; centered at (0,0,0) with radius=1.0; head is in +X
		colorRGBA const color(BLACK);
		float const body_zval(-0.3);
		cube_t abdomen(point(-0.8, 0.0, body_zval)), body(point(0.0, 0.0, body_zval));
		abdomen.expand_by(vector3d(0.50, 0.35, 0.35));
		body   .expand_by(vector3d(0.45, 0.30, 0.20));
		mat.add_sphere_to_verts(abdomen, color);
		mat.add_sphere_to_verts(body,    color);
		assign_tc_range(0.0, 0.0, 0.0); // head and body aren't animated
		float const leg_radius(0.03);

		for (unsigned d = 0; d < 2; ++d) { // {left, right}
			float const d_sign(d ? -1.0 : 1.0);
			add_eye(point(0.30, 0.080*d_sign, body_zval+0.14), 0.026);
			add_eye(point(0.40, 0.045*d_sign, body_zval+0.08), 0.028);
			add_eye(point(0.44, 0.020*d_sign, body_zval+0.04), 0.016);
			add_eye(point(0.43, 0.055*d_sign, body_zval+0.03), 0.015);
			float const fang_radius(0.05);
			point const fang_top(0.44, 0.05*d_sign, body_zval-0.04), fang_bot(fang_top - vector3d(0.0, 0.0, 0.2));
			mat.add_sphere_to_verts(fang_top, vector3d(fang_radius, fang_radius, fang_radius), color); // top of fang
			mat.add_cylin_to_verts (fang_bot, fang_top, 0.0, fang_radius, color, 0, 0); // fang
			// hourglass shape? colorRGBA(0.7, 0.2, 0.0)
			assign_tc_range(0.0, 0.0, 0.0); // not animated

			// add legs
			for (unsigned n = 0; n < 4; ++n) {
				float const ts(n/4.0);
				point const joint(0.12*(n - 1.5), 0.26*d_sign, body_zval);
				point const knee (2.0*joint.x, 2.0*joint.y,  0.5);
				point const ankle(2.8*knee .x, 2.8*knee .y,  0.0);
				point const foot (3.5*knee .x, 3.5*knee .y, -1.0);
				vector3d const sphere_radius(leg_radius, leg_radius, leg_radius);
				float const joint_tt(0.0*d_sign), knee_tt(0.3*d_sign), ankle_tt(0.7*d_sign), foot_tt(1.0*d_sign);
				mat.add_sphere_to_verts(joint, sphere_radius, color); // round body joint
				assign_tc_range(ts, joint_tt, joint_tt);
				mat.add_cylin_to_verts(joint, knee, leg_radius, leg_radius, color, 0, 0);
				assign_tc_range(ts, joint_tt, knee_tt);
				mat.add_sphere_to_verts(knee, sphere_radius, color); // round knee joint
				assign_tc_range(ts, knee_tt, knee_tt);
				mat.add_cylin_to_verts(ankle, knee, leg_radius, leg_radius, color, 0, 0);
				assign_tc_range(ts, ankle_tt, knee_tt);
				mat.add_sphere_to_verts(ankle, sphere_radius, color); // round ankle joint
				assign_tc_range(ts, ankle_tt, ankle_tt);
				mat.add_cylin_to_verts(foot,  ankle, 0.1*leg_radius, leg_radius, color, 0, 0);
				assign_tc_range(ts, foot_tt, ankle_tt);
			} // for n
		} // for d
		mat.create_vbo_inner();
		mat.clear_vectors(1); // free_memory=1: vector data no longer needed
		is_setup = 1;
	}
	void assign_tc_range(float ts, float tt_lo, float tt_hi) {
		// tcs are used for animation:
		// ts is the position of the leg from front to back: {0.0, 0.25, 0.5, 0.75}
		// tt is the joint index: negative for left, positive for right; 0.0 for body joint, 0.5 for knee, 1.0 for foot
		for (auto i = mat.itri_verts.begin()+cur_vert_pos; i != mat.itri_verts.end(); ++i) {
			i->t[0] = ts;
			i->t[1] = i->t[1]*tt_hi + (1.0 - i->t[1])*tt_lo; // existing t[1] is 0.0 for low end and 1.0 for high end
		}
		cur_vert_pos = mat.itri_verts.size();
	}
public:
	spider_draw_t() : cur_vert_pos(0), is_setup(0) {}
	void clear() {mat.clear(); is_setup = 0;}

	void draw(vect_spider_t const &spiders, shader_t &s, bool shadow_only) {
		if (spiders.empty()) return; // nothing to draw
		if (!is_setup) {init();}
		mat.vao_setup(shadow_only);
		s.set_specular(0.2, 80.0); // FIXME: building interior lighting specular looks wrong for rotated models and with sunlight
		select_texture(WHITE_TEX);
		int const animation_id = 8; // custom spider animation
		s.add_uniform_int("animation_id", animation_id);

		for (spider_t const &S : spiders) { // future work: use instancing
			s.add_uniform_float("animation_time", S.anim_time);
			fgPushMatrix();
			fgTranslate(S.pos.x, S.pos.y, (S.pos.z + S.radius)); // shift up by radius when drawing
			city_model_loader_t::rotate_model_from_plus_x_to_dir(S.dir);
			uniform_scale(S.radius);
			check_mvm_update();
			mat.draw_inner(shadow_only);
			fgPopMatrix();
		} // for S
		s.add_uniform_float("animation_time", 0.0); // reset animation time
		s.add_uniform_int("animation_id", 0); // clear animation
		s.clear_specular();
		indexed_vao_manager_with_shadow_t::post_render();
	}
};
spider_draw_t spider_draw;

// Note: non-const because it creates the VBO; inc_small: 0=large only, 1=large+small, 2=large+small+detail
void building_room_geom_t::draw(brg_batch_draw_t *bbd, shader_t &s, building_t const &building, occlusion_checker_noncity_t &oc, vector3d const &xlate,
	unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building)
{
	if (empty()) return; // no geom
	unsigned const num_screenshot_tids(get_num_screenshot_tids());
	static int last_frame(0);
	static unsigned num_geom_this_frame(0); // used to limit per-frame geom gen time; doesn't apply to shadow pass, in case shadows are cached
	if (frame_counter > last_frame) {num_geom_this_frame = 0; last_frame = frame_counter;}
	point const camera_bs(camera_pdu.pos - xlate);
	bool const draw_lights(camera_bs.z < building.bcube.z2()); // don't draw ceiling lights when player is above the building
	if (player_in_building) {bbd = nullptr;} // use immediate drawing when player is in the building because draw order matters for alpha blending
	if (bbd != nullptr) {bbd->set_camera_dir_mask(camera_bs, building.bcube);}

	if (lighting_invalid) { // set in set_obj_lit_state_to()
		clear_lit_materials();
		lighting_invalid = 0;
	}
	if (lights_changed) {
		mats_lights.clear();
		lights_changed = 0;
	}
	if (has_pictures && num_pic_tids != num_screenshot_tids) {
		clear_static_vbos(); // user created a new screenshot texture, and this building has pictures - recreate room geom
		num_pic_tids = num_screenshot_tids;
	}
	// generate vertex data in the shadow pass or if we haven't hit our generation limit; must be consistent for static and small geom
	if (shadow_only || num_geom_this_frame < MAX_ROOM_GEOM_GEN_PER_FRAME) {
		if (mats_static.empty()) { // create static materials if needed
			create_obj_model_insts(building);
			create_static_vbos(building);
			++num_geom_this_frame;
		}
		if (inc_small && mats_small.empty()) { // create small materials if needed
			create_small_static_vbos(building);
			++num_geom_this_frame;
		}
		// Note: not created on the shadow pass, because trim_objs may not have been created yet and we would miss including it
		if (inc_small == 2 && !shadow_only && mats_detail.empty()) { // create detail materials if needed
			create_detail_vbos(building);
			++num_geom_this_frame;
		}
	}
	if (draw_lights && mats_lights .empty()) {create_lights_vbos (building);} // create lights  materials if needed (no limit)
	if (inc_small   && mats_dynamic.empty()) {create_dynamic_vbos(building);} // create dynamic materials if needed (no limit); drawn with small objects
	if (mats_doors.empty()) {create_door_vbos(building);} // create door materials if needed (no limit)
	enable_blend(); // needed for rugs and book text
	assert(s.is_setup());
	mats_static.draw(bbd, s, shadow_only, reflection_pass); // this is the slowest call
	if (draw_lights)    {mats_lights .draw(bbd, s, shadow_only, reflection_pass);}
	if (inc_small  )    {mats_dynamic.draw(bbd, s, shadow_only, reflection_pass);}
	if (inc_small == 2) {mats_detail .draw(bbd, s, shadow_only, reflection_pass);}
	mats_doors.draw(bbd, s, shadow_only, reflection_pass);
	
	if (inc_small) {
		mats_small.draw(bbd, s, shadow_only, reflection_pass);

		if (player_in_building) { // if we're not in the building, don't draw alpha mask materials at all; without the special shader they won't look correct when drawn through windows
			if (shadow_only) {
				shader_t amask_shader;
				amask_shader.begin_simple_textured_shader(0.9); // need to use texture with alpha test
				mats_amask.draw(nullptr, amask_shader, 2, 0); // shadow pass with alpha mask; no brg_batch_draw
				s.make_current(); // switch back to the normal shader
			}
			else if (reflection_pass) {
				mats_amask.draw(nullptr, s, 0, 1); // no brg_batch_draw
			}
			else { // this is expensive: only enable for the current building and the main draw pass
				shader_t amask_shader;
				setup_building_draw_shader(amask_shader, 0.9, 1, 1, 0); // min_alpha=0.9, enable_indir=1, force_tsl=1, use_texgen=0
				mats_amask.draw(nullptr, amask_shader, 0, 0); // no brg_batch_draw
				s.make_current(); // switch back to the normal shader
			}
		}
	}
	disable_blend();
	indexed_vao_manager_with_shadow_t::post_render();
	bool const disable_cull_face(0); // better but slower?
	if (disable_cull_face) {glDisable(GL_CULL_FACE);}
	point const building_center(building.bcube.get_cube_center());
	bool const is_rotated(building.is_rotated());
	oc.set_exclude_bix(building_ix);
	bool obj_drawn(0);
	water_sound_manager_t water_sound_manager(camera_bs);
	bool const check_clip_cube(shadow_only && !is_rotated && !smap_light_clip_cube.is_all_zeros()); // check clip cube for shadow pass; not implemented for rotated buildings
	bool const check_occlusion(display_mode & 0x08);

	// draw object models
	for (auto i = obj_model_insts.begin(); i != obj_model_insts.end(); ++i) {
		room_object_t &obj(get_room_object_by_index(i->obj_id));
		if (!player_in_building && !shadow_only && obj.is_interior()) continue; // don't draw objects in interior rooms if the player is outside the building (useful for office bathrooms)
		if (check_clip_cube && !smap_light_clip_cube.intersects(obj + xlate)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
		point obj_center(obj.get_cube_center());
		if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
		if (!shadow_only && !dist_less_than(camera_bs, obj_center, 100.0*obj.dz())) continue; // too far away (obj.max_len()?)
		bool const is_sink(player_in_building && !shadow_only && obj.type == TYPE_SINK);
		if (is_sink) {water_sound_manager.register_running_water(obj, building);}
		if (!(is_rotated ? building.is_rot_cube_visible(obj, xlate) : camera_pdu.cube_visible(obj + xlate))) continue; // VFC
		if (check_occlusion && building.check_obj_occluded(obj, camera_bs, oc, reflection_pass)) continue;
		bool const is_emissive(!shadow_only && obj.type == TYPE_LAMP && obj.is_lit());
		if (is_emissive) {s.set_color_e(LAMP_COLOR*0.4);}
		apply_room_obj_rotate(obj, *i); // Note: may modify obj by clearing flags
		bool const use_low_z_bias(obj.type == TYPE_CUP && !shadow_only);
		bool const untextured(obj.flags & RO_FLAG_UNTEXTURED);
		
		if (use_low_z_bias) {
			s.add_uniform_float("norm_bias_scale",   5.0); // half the default value
			s.add_uniform_float("dlight_pcf_offset", 0.5*cur_dlight_pcf_offset);
		}
		// Note: lamps are the most common and therefore most expensive models to draw
		building_obj_model_loader.draw_model(s, obj_center, obj, i->dir, obj.color, xlate, obj.get_model_id(), shadow_only, 0, 0, 0, untextured);
		if (!shadow_only && obj.type == TYPE_STOVE) {draw_stove_flames(obj, camera_bs, s);} // draw blue burner flame
		
		if (use_low_z_bias) { // restore to the defaults
			s.add_uniform_float("norm_bias_scale",   10.0);
			s.add_uniform_float("dlight_pcf_offset", cur_dlight_pcf_offset);
		}
		if (is_emissive) {s.set_color_e(BLACK);}
		if (is_sink) {water_draw.add_water_for_sink(obj);}
		obj_drawn = 1;
	} // for i
	if (player_in_building) { // only drawn for the player building
		if (!shadow_only && !reflection_pass) { // these models aren't drawn in the shadow or reflection passes; no emissive or rotated objects
			for (auto i = model_objs.begin(); i != model_objs.end(); ++i) {
				point obj_center(i->get_cube_center());
				if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
				if (!shadow_only && !dist_less_than(camera_bs, obj_center, 100.0*i->max_len())) continue; // too far away
				if (!(is_rotated ? building.is_rot_cube_visible(*i, xlate) : camera_pdu.cube_visible(*i + xlate))) continue; // VFC
				if (check_occlusion && building.check_obj_occluded(*i, camera_bs, oc, reflection_pass)) continue;
				vector3d dir(i->get_dir());
				if (is_rotated) {building.do_xy_rotate_normal(dir);}
				building_obj_model_loader.draw_model(s, obj_center, *i, dir, i->color, xlate, i->get_model_id(), shadow_only, 0, 0);
				obj_drawn = 1;
			} // for model_objs
		}
		if (!rats.empty()) {
			if (!shadow_only) {
				int const animation_id = 7; // custom rat animation
				s.add_uniform_int("animation_id", animation_id);
			}
			//spiders.clear(); // FIXME

			for (rat_t &rat : rats) {
				cube_t const bcube(rat.get_bcube());
				if (check_clip_cube && !smap_light_clip_cube.intersects(bcube + xlate)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
				if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
				if (check_occlusion && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass)) continue;
	#if 0 // FIXME: hack to draw rats as spiders
				float const radius(0.5*rat.radius);
				spiders.push_back(spider_t((rat.pos + vector3d(0.0, 0.0, radius)), radius, rat.dir));
				spiders.back().anim_time = 20.0*rat.anim_time;
				continue;
	#endif
				point const pos(bcube.get_cube_center());
				bool const animate(rat.anim_time > 0.0 && !shadow_only); // can't see the animation in the shadow pass anyway
				if (!shadow_only) {s.add_uniform_float("animation_time", rat.anim_time);}
				colorRGBA const color(rat_color); // make the rat's fur darker
				//colorRGBA const color(blend_color(RED, WHITE, rat.fear, 0)); // used for debugging fear
				//colorRGBA const color(blend_color(RED, WHITE, rat.attacking, 0));
				cube_t const rat_bcube(rat.get_bcube_with_dir());
				building_obj_model_loader.draw_model(s, pos, rat_bcube, rat.dir, color, xlate, OBJ_MODEL_RAT, shadow_only, 0, animate);

				if (rat.attacking) { // draw red glowing eyes
					s.set_color_e(colorRGBA(0.5, 0.0, 0.0, 1.0)); // light emissive red
					s.set_cur_color(RED);
					select_texture(WHITE_TEX);
					s.add_uniform_float("animation_time", 0.0); // clear animations
					point eyes_center(pos + vector3d(0.0, 0.0, 0.09*rat.height) + 0.85*rat.get_hlength()*rat.dir);
					vector3d const eye_sep_dir(0.21*rat.hwidth*cross_product(rat.dir, plus_z).get_norm());

					for (unsigned d = 0; d < 2; ++d) { // draw left and right eye, untextured
						draw_sphere_vbo((eyes_center + (d ? 1.0 : -1.0)*eye_sep_dir), 0.05*rat.height, 16, 0);
					}
					s.set_color_e(BLACK);
				}
				obj_drawn = 1;
			} // for rat
			if (!shadow_only) {s.add_uniform_int("animation_id", 0);} // reset
		} // end rats drawing
		spider_draw.draw(spiders, s, shadow_only);
	}
	if (disable_cull_face) {glEnable(GL_CULL_FACE);}
	if (obj_drawn) {check_mvm_update();} // needed after popping model transform matrix

	if (player_in_building && !shadow_only) { // draw water for sinks that are turned on
		auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->type != TYPE_KSINK && i->type != TYPE_BRSINK) continue; // TYPE_SINK is handled above
			water_sound_manager.register_running_water(*i, building);
			water_draw.add_water_for_sink(*i);
		}
	}
	water_sound_manager.finalize();
	water_draw.draw_and_clear(s);

	if (player_in_building && !shadow_only && player_held_object.is_valid()) {
		// draw the item the player is holding; pre_smap_player_pos should be the correct position for reflections
		point const obj_pos((reflection_pass ? pre_smap_player_pos : camera_bs) + CAMERA_RADIUS*cview_dir - vector3d(0.0, 0.0, 0.5*CAMERA_RADIUS));
		player_held_object.translate(obj_pos - player_held_object.get_cube_center());
		if (player_held_object.type == TYPE_LG_BALL) {draw_lg_ball_in_building(player_held_object, s);} // the only supported dynamic object type
		else if (player_held_object.can_use()) {draw_interactive_player_obj(player_held_object, s, xlate);}
		else {assert(0);}
	}
	// alpha blended, should be drawn near last
	decal_manager.draw_building_interior_decals(player_in_building, shadow_only); // draw decals in this building
	
	if (!shadow_only && !mats_alpha.empty()) { // draw last; not shadow casters; for shower glass, etc.
		enable_blend();
		glDepthMask(GL_FALSE); // disable depth writing
		mats_alpha.draw(bbd, s, shadow_only, reflection_pass);
		glDepthMask(GL_TRUE);
		disable_blend();
		indexed_vao_manager_with_shadow_t::post_render();
	}
}

template<bool check_sz> bool are_pts_occluded_by_any_cubes(point const &pt, point const *const pts, unsigned npts, vect_cube_t const &cubes, unsigned dim, float min_sz=0.0) {
	assert(npts > 0);

	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (check_sz && c->get_sz_dim(!dim) < min_sz) break; // too small an occluder; since cubes are sorted by size in this dim, we can exit the loop here
		if (dim <= 2 && (pt[dim] < c->d[dim][0]) == (pts[0][dim] < c->d[dim][0])) continue; // skip if cube face does not separate pt from the first point (dim > 2 disables)
		if (!check_line_clip(pt, pts[0], c->d)) continue; // first point does not intersect
		bool not_occluded(0);

		for (unsigned p = 1; p < npts; ++p) { // skip first point
			if (!check_line_clip(pt, pts[p], c->d)) {not_occluded = 1; break;}
		}
		if (!not_occluded) return 1;
	} // for c
	return 0;
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
bool check_cube_occluded(cube_t const &cube, vect_cube_t const &occluders, point const &viewer) {
	if (occluders.empty()) return 0;
	point pts[8];
	unsigned const npts(get_cube_corners(cube.d, pts, viewer, 0)); // should return only the 6 visible corners
	return are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, occluders, 3); // set invalid dim of 3 because cubes are of mixed dim and we can't use that optimization
}
struct comp_car_by_dist {
	vector3d const &viewer;
	comp_car_by_dist(vector3d const &viewer_) : viewer(viewer_) {}
	bool operator()(car_t const &c1, car_t const &c2) const {
		return (p2p_dist_xy_sq(c1.bcube.get_cube_center(), viewer) > p2p_dist_xy_sq(c2.bcube.get_cube_center(), viewer));
	}
};

bool building_t::has_cars_to_draw(bool player_in_building) const {
	if (!has_room_geom()) return 0;
	if (player_in_building && has_parking_garage) return 1; // parking garage cars are drawn if the player is in the building
	if (interior->room_geom->has_garage_car)      return 1; // have car in a garage
	return 0;
}
void building_t::draw_cars_in_building(shader_t &s, vector3d const &xlate, bool player_in_building, bool shadow_only) const {
	assert(has_room_geom());
	point viewer(camera_pdu.pos - xlate);
	bool const player_in_this_building(player_in_building && bcube.contains_pt(viewer)); // check before rotating
	bool const check_occlusion(display_mode & 0x08);
	float const floor_spacing(get_window_vspace());
	vect_room_object_t const &objs(interior->room_geom->objs);
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	unsigned const pg_wall_start(interior->room_geom->wall_ps_start);
	assert(pg_wall_start < objs.size());
	static vector<car_t> cars_to_draw; // reused across frames
	cars_to_draw.clear();

	if (interior->room_geom->has_garage_car && player_in_basement != 2) { // car in a house garage
		room_object_t const &obj(objs[pg_wall_start]);
		assert(obj.type == TYPE_PARK_SPACE); // must be a parking space

		// skip if viewer is in this building and on a different floor
		if ((obj.flags & RO_FLAG_USED) && (!player_in_this_building || !check_occlusion || !has_int_garage ||
			int((viewer.z - bcube.z1())/floor_spacing) == int((obj.z2() - bcube.z1())/floor_spacing)))
		{
			car_t car(car_from_parking_space(obj));
			if (camera_pdu.cube_visible(car.bcube + xlate)) {cars_to_draw.push_back(car);}
		}
	}
	else if (player_in_this_building && has_parking_garage) { // cars in parking garages
		// only draw if the player or light is in the basement, or the player is on the first floor where a car may be visible through the stairs
		float max_vis_zval(ground_floor_z1);
		if (!shadow_only) {max_vis_zval += floor_spacing;} // player on first floor?
		if (viewer.z > max_vis_zval) return;
		maybe_inv_rotate_point(viewer); // not needed because there are no cars in rotated buildings?

		// start at walls, since parking spaces are added after those
		for (auto i = (objs.begin() + pg_wall_start); i != objs_end; ++i) {
			if (i->type != TYPE_PARK_SPACE) continue;
			if (!(i->flags & RO_FLAG_USED)) continue; // no car in this space
			if (i->z2() < viewer.z - 2.0*floor_spacing) continue; // move than a floor below - skip
			car_t car(car_from_parking_space(*i));
			if (!shadow_only && check_occlusion && viewer.z > ground_floor_z1 && !line_intersect_stairs_or_ramp(viewer, car.get_center())) continue;
			if (camera_pdu.cube_visible(car.bcube + xlate)) {cars_to_draw.push_back(car);}
		}
		if (cars_to_draw.empty()) return;

		if (check_occlusion) {
			vect_cube_t occluders; // should this be split out per PG level?

			for (auto i = (objs.begin() + pg_wall_start); i != objs_end; ++i) {
				if (i->type != TYPE_PG_WALL || i->item_flags != 0) continue; // not parking garage wall (breaking is incorrect for multiple PG levels)
				occluders.push_back(*i);
			}
			// gather occluders from parking garage ceilings and floors (below ground floor)
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
		if (shadow_only) {s.begin_color_only_shader();} // this path should be unused
		else {setup_building_draw_shader(s, 0.0, 1, 0, 0);} // min_alpha=0.0, enable_indir=1, force_tsl=0, use_texgen=0
	}
	for (auto &car : cars_to_draw) {draw_car_in_pspace(car, s, xlate, shadow_only);}
}

// Note: c is in local building space and viewer_in is in non-rotated building space
bool building_t::check_obj_occluded(cube_t const &c, point const &viewer_in, occlusion_checker_noncity_t &oc, bool reflection_pass, bool c_is_building_part) const {
	if (!interior) return 0; // could probably make this an assert
	//highres_timer_t timer("Check Object Occlusion"); // 0.001ms
	point viewer(viewer_in);
	maybe_inv_rotate_point(viewer); // rotate viewer pos into building space
	// if fully inside basement, and viewer outside building or not on first floor, will be occluded
	if (c.z2() < ground_floor_z1 && (viewer.z > (ground_floor_z1 + get_window_vspace()) || !bcube.contains_pt(viewer))) return 1;
	point pts[8];
	unsigned const npts(get_cube_corners(c.d, pts, viewer, 0)); // should return only the 6 visible corners
	
	if (!reflection_pass && !c_is_building_part) {
		// check walls of this building; not valid for reflections because the reflected camera may be on the other side of a wall/mirror
		for (unsigned d = 0; d < 2; ++d) {
			if (are_pts_occluded_by_any_cubes<1>(viewer, pts, npts, interior->walls[d], d, c.get_sz_dim(!d))) return 1; // with size check (helps with light bcubes)
		}
	}
	if (!c_is_building_part && (reflection_pass || bcube.contains_pt(viewer))) {
		// viewer inside this building; includes shadow_only case and reflection_pass (even if reflected camera is outside the building);
		// check floors/ceilings of this building
		if (fabs(viewer.z - c.zc()) > (reflection_pass ? 1.0 : 0.5)*get_window_vspace()) { // on different floors
			if (are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, interior->fc_occluders, 2)) return 1;
		}
	}
	else if (camera_in_building) { // player in some other building
		if (is_rotated()) return 0; // not implemented yet - need to rotate viewer and pts into coordinate space of player_building

		if (player_building != nullptr && player_building->interior) { // check walls of the building the player is in
			if (player_building != this) { // otherwise player_in_this_building should be true; note that we can get here from building_t::add_room_lights()
				for (unsigned d = 0; d < 2; ++d) { // check walls of the building the player is in; can't use min_sz due to perspective effect of walls near the camera
					if (are_pts_occluded_by_any_cubes<0>(viewer, pts, npts, player_building->interior->walls[d], d)) return 1;
				}
			}
		}
	}
	else if (viewer.z < bcube.z2()) { // player not in a building and not above this building
		if (is_rotated()) return 0; // not implemented yet - c is not an axis aligned cube in global coordinate space
		if (oc.is_occluded(c)) return 1; // check other buildings
	}
	else if (!c_is_building_part && is_simple_cube()) { // player above this building; check if object is occluded by the roof
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			cube_t roof(*p);
			roof.z1() = roof.z2() - get_fc_thickness();
			bool not_occluded(0);

			for (unsigned p = 0; p < npts; ++p) {
				if (!check_line_clip(viewer, pts[p], roof.d)) {not_occluded = 1; break;}
			}
			if (!not_occluded) return 1;
		} // for p
		return 0;
	}
	return 0;
}

bool building_t::is_entire_building_occluded(point const &viewer, occlusion_checker_noncity_t &oc) const {
	if (is_rotated()) return 0; // not handled (optimization)

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
		if (has_basement() && (i - parts.begin()) == (int)basement_part_ix) continue; // skip the basement, which isn't visible from outside the building
		if (!check_obj_occluded(*i, viewer, oc, 0, 1)) return 0; // c_is_building_part=1
	}
	return 1; // all parts occluded
}

void paint_draw_t::draw_paint() const {
	if (qbd[0].empty() && qbd[1].empty()) return; // nothing to do
	glDepthMask(GL_FALSE); // disable depth write
	enable_blend();

	if (!qbd[0].empty()) {
		select_texture(BLUR_CENT_TEX); // spraypaint - smooth alpha blended edges
		qbd[0].draw();
	}
	if (!qbd[1].empty()) {
		select_texture(get_texture_by_name("circle.png", 0, 0, 1, 0.0, 1, 1, 1)); // markers - sharp edges, used as alpha mask with white background color
		qbd[1].draw();
	}
	disable_blend();
	glDepthMask(GL_TRUE);
}

void building_decal_manager_t::commit_pend_tape_qbd() {
	pend_tape_qbd.add_quads(tape_qbd);
	pend_tape_qbd.clear();
}
void building_decal_manager_t::draw_building_interior_decals(bool player_in_building, bool shadow_only) const {
	if (shadow_only) { // shadow pass, draw tape only
		if (player_in_building) {
			tape_qbd.draw(); // somewhat inefficient, since we have to send all the data for every light source
			pend_tape_qbd.draw();
		}
		return;
	}
	paint_draw[1].draw_paint(); // draw exterior paint always - this will show up on windows (even when looking outside into another part of the same building)
	if (!player_in_building) return;
	paint_draw[0].draw_paint(); // draw interior paint

	if (!tp_qbd.empty()) { // toilet paper squares: double sided, lit from top
		glDisable(GL_CULL_FACE); // draw both sides
		select_texture(WHITE_TEX);
		tp_qbd.draw(); // use a VBO for this if the player leaves the building and then comes back?
		glEnable(GL_CULL_FACE);
	}
	if (!tape_qbd.empty() || !pend_tape_qbd.empty()) { // tape lines: single sided so that lighting works, both sides drawn independently
		select_texture(WHITE_TEX);
		tape_qbd.draw();
		pend_tape_qbd.draw();
	}
	if (!blood_qbd.empty()) {
		select_texture(BLOOD_SPLAT_TEX);
		glDepthMask(GL_FALSE); // disable depth write
		enable_blend();
		blood_qbd.draw();
		disable_blend();
		glDepthMask(GL_TRUE);
	}
}


