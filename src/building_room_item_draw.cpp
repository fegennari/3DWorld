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

vect_room_object_t pending_objs;
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
void set_car_model_color(car_t &car);
bldg_obj_type_t get_taken_obj_type(room_object_t const &obj);

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
		unsigned const ixs_off[6] = {1,2,0, 3,2,1}; // 1 quad = 2 triangles
		bool const flat_sides(ndiv <= 6 && side_tscale == 0.0); // hack to draw bolts untextured with flat sides, since no other cylinders have only 6 sides
		unsigned const num_side_verts(flat_sides ? 4*ndiv : 2*(ndiv+1)), unique_verts_per_side(flat_sides ? 4 : 2);
		itri_verts.resize(itris_start + num_side_verts);
		indices.resize(ixs_start + 6*ndiv);

		if (flat_sides) {
			for (unsigned i = 0; i < ndiv; ++i) { // vertex data
				unsigned const in((i+1)%ndiv);
				point const pts[4] = {vpn.p[(i<<1)+0], vpn.p[(i<<1)+1], vpn.p[(in<<1)+0], vpn.p[(in<<1)+1]};
				norm_comp const normal(get_poly_norm(pts));
				for (unsigned n = 0; n < 4; ++n) {itri_verts[itix++].assign(pts[n], normal, 0.0, 0.0, cw);} // all tcs=0
			}
		}
		else {
			for (unsigned i = 0; i <= ndiv; ++i) { // vertex data
				unsigned const s(i%ndiv);
				float const ts(side_tscale*(1.0f - i*ndiv_inv) + side_tscale_add);
				norm_comp const normal(0.5*(vpn.n[s] + vpn.n[(i+ndiv-1)%ndiv])); // normalize?
				itri_verts[itix++].assign(vpn.p[(s<<1)+0], normal, (swap_txy ? 0.0 : ts), (swap_txy ? ts : 0.0), cw);
				itri_verts[itix++].assign(vpn.p[(s<<1)+1], normal, (swap_txy ? 1.0 : ts), (swap_txy ? ts : 1.0), cw);
			}
		}
		for (unsigned i = 0; i < ndiv; ++i) { // index data
			unsigned const ix0(itris_start + unique_verts_per_side*i);
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
	unsigned get_mem_usage() const {
		unsigned mem(free_list.size()*sizeof(rgeom_storage_t));
		for (auto i = free_list.begin(); i != free_list.end(); ++i) {
			//cout << i->tex.tid << "\t" << i->tex.shadowed << "\t" << (i->quad_verts.capacity() + i->itri_verts.capacity()) << "\t" << i->get_mem_usage() << endl; // TESTING
			mem += i->get_mem_usage();
		}
		return mem;
	}
	unsigned size() const {return free_list.size();}
};
rgeom_alloc_t rgeom_alloc; // static allocator with free list, shared across all buildings; not thread safe


vbo_cache_t::vbo_cache_entry_t vbo_cache_t::alloc(unsigned size, bool is_index) {
	assert(size > 0); // not required, but a good sanity check
	auto &e(entries[is_index]);
	//if ((v_used % 1000) == 0) {print_stats();} // TESTING
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
	if (best_fit != e.end()) {
		vbo_cache_entry_t const ret(*best_fit); // deep copy
		swap(*best_fit, e.back());
		e.pop_back();
		++v_reuse; s_reuse += ret.size; ++v_used; s_used += ret.size; // Note: s_reuse can overflow
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
		for (vbo_cache_entry_t &entry : entries[d]) {delete_vbo(entry.vbo);}
		entries[d].clear();
	}
}
void vbo_cache_t::print_stats() const {
	cout << "VBOs: A " << v_alloc << " U " << v_used << " R " << v_reuse << " F " << v_free
		 << "  SZ: A " << (s_alloc>>20) << " U " << (s_used>>20) << " R " << (s_reuse>>20) << " F " << (s_free>>20) << endl; // in MB
}

/*static*/ vbo_cache_t rgeom_mat_t::vbo_cache;

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
	unsigned const qsz(quad_verts.size()*sizeof(vertex_t)), itsz(itri_verts.size()*sizeof(vertex_t)), tot_verts_sz(qsz + itsz);
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
		}
		else { // existing
			vert_vbo_sz = vret.size;
			assert(tot_verts_sz <= vert_vbo_sz);
		}
		if (iret.size == 0) { // newly created
			ixs_vbo_sz = ix_data_sz;
			upload_to_vbo(vao_mgr.ivbo, indices, 1, 1);
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
}

void brg_batch_draw_t::set_camera_dir_mask(point const &camera_bs, cube_t const &bcube) {
	camera_dir_mask = 0;
	for (unsigned d = 0; d < 3; ++d) {
		if (camera_bs[d] < bcube.d[d][1]) {camera_dir_mask |=  1<<(2*d);}
		if (camera_bs[d] > bcube.d[d][0]) {camera_dir_mask |= 1<<(2*d+1);}
	}
}
void brg_batch_draw_t::add_material(rgeom_mat_t const &m) {
	for (auto &i : to_draw) { // check all existing materials for a matching texture, etc.
		if (i.tex.is_compat_ignore_shadowed(m.tex)) {i.mats.push_back(&m); return;} // found existing material
	}
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
	if (num_verts == 0) return; // Note: should only happen when reusing materials and all objects using this material were removed
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
void rgeom_mat_t::pre_draw(int shadow_only) const {
	vao_mgr.pre_render(shadow_only != 0);
}
void rgeom_mat_t::draw_geom() const {
	glDrawRangeElements(GL_TRIANGLES, 0, num_verts, num_ixs, GL_UNSIGNED_INT, nullptr);
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
	draw(state, nullptr, 0, 0); // no brg_batch_draw_t
	clear();
	indexed_vao_manager_with_shadow_t::post_render();
}

void building_materials_t::clear() {
	invalidate();
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
		if (m->get_tot_vert_capacity() == 0) {rgeom_alloc.alloc_safe(*m);} // existing but empty entry, allocate capacity from the allocator free list
		return *m;
	}
	emplace_back(tex); // not found, add a new material
	if (inc_shadows) {back().enable_shadows();}
	rgeom_alloc.alloc_safe(back());
	return back();
}
void building_materials_t::create_vbos(building_t const &building) {
	for (iterator m = begin(); m != end(); ++m) {m->create_vbo(building);}
	valid = 1;
}
void building_materials_t::draw(brg_batch_draw_t *bbd, shader_t &s, int shadow_only, bool reflection_pass) {
	if (!valid) return; // pending generation of data, don't draw yet
	//highres_timer_t timer("Draw Materials"); // 0.0168
	tid_nm_pair_dstate_t state(s);
	for (iterator m = begin(); m != end(); ++m) {m->draw(state, bbd, shadow_only, reflection_pass);}
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
	mats_static .clear();
	mats_alpha  .clear();
	mats_small  .clear();
	mats_text   .clear();
	mats_amask  .clear();
	mats_dynamic.clear();
	mats_doors  .clear();
	mats_lights .clear();
	mats_detail .clear();
	obj_model_insts.clear(); // these are associated with static VBOs
}
// Note: used for room lighting changes; detail object changes are not supported
void building_room_geom_t::check_invalid_draw_data() {
	if (invalidate_mats_mask & (1 << MAT_TYPE_SMALL  )) { // small objects
		mats_small.invalidate();
		mats_amask.invalidate();
		mats_text .invalidate(); // Note: for now text is assigned to type MAT_TYPE_SMALL since it's always drawn with small objects
	}
	if (invalidate_mats_mask & (1 << MAT_TYPE_STATIC )) { // large objects and 3D models
		mats_static.invalidate();
		mats_alpha .invalidate();
		obj_model_insts.clear();
	}
	//if (invalidate_mats_mask & (1 << MAT_TYPE_TEXT   )) {mats_text   .invalidate();} // text objects
	if (invalidate_mats_mask & (1 << MAT_TYPE_DYNAMIC)) {mats_dynamic.invalidate();} // dynamic objects
	if (invalidate_mats_mask & (1 << MAT_TYPE_DOORS  )) {mats_doors  .invalidate();}
	if (invalidate_mats_mask & (1 << MAT_TYPE_LIGHTS )) {mats_lights .invalidate();}
	invalidate_mats_mask = 0; // reset for next frame
}
void building_room_geom_t::invalidate_draw_data_for_obj(room_object_t const &obj, bool was_taken) {
	if (obj.is_dynamic() || (obj.type == TYPE_BUTTON && (obj.flags & RO_FLAG_IN_ELEV))) { // elevator buttons are drawn as dynamic objects
		update_dynamic_draw_data();
		return;
	}
	bldg_obj_type_t const type(was_taken ? get_taken_obj_type(obj) : get_room_obj_type(obj));
	if ( type.lg_sm & 2 )                  {invalidate_small_geom ();} // small objects
	if ((type.lg_sm & 1) || type.is_model) {invalidate_static_geom();} // large objects and 3D models
}
// Note: called when adding, removing, or moving objects
void building_room_geom_t::update_draw_state_for_room_object(room_object_t const &obj, building_t &building, bool was_taken) {
	invalidate_draw_data_for_obj(obj, was_taken);
	//if (type.ai_coll) {building.invalidate_nav_graph();} // removing this object should not affect the AI navigation graph
	modified_by_player = 1; // flag so that we avoid re-generating room geom if the player leaves and comes back
}

rgeom_mat_t &building_room_geom_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows, bool dynamic, unsigned small, bool transparent) {
	// small: 0=mats_static, 1=mats_small, 2=mats_detail
	if (small == 1 && tex.tid == FONT_TEXTURE_ID) {return mats_text.get_material(tex, inc_shadows);} // inc_shadows should be 0
	return (dynamic ? mats_dynamic : (small ? ((small == 2) ? mats_detail : mats_small) : (transparent ? mats_alpha : mats_static))).get_material(tex, inc_shadows);
}
rgeom_mat_t &building_room_geom_t::get_metal_material(bool inc_shadows, bool dynamic, unsigned small, colorRGBA const &spec_color) {
	tid_nm_pair_t tex(-1, 1.0, inc_shadows);
	tex.set_specular_color(spec_color, 0.8, 60.0);
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
		case TYPE_BOOK:    add_book    (*i, 1, 0, 0); break; // lg
		case TYPE_BCASE:   add_bookcase(*i, 1, 0, 0, tscale, 0); break; // lg
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
	//highres_timer_t timer2("Gen Room Geom VBOs"); // < 2ms
	mats_static.create_vbos(building);
	mats_alpha .create_vbos(building);
	//cout << "static: size: " << rgeom_alloc.size() << " mem: " << rgeom_alloc.get_mem_usage() << endl; // start=50MB, peak=76MB
}

void building_room_geom_t::create_small_static_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Small"); // 7.8ms, slow building at 26,16
	model_objs.clear(); // currently model_objs are only created for small objects in drawers, so we clear this here
	add_small_static_objs_to_verts(expanded_objs);
	add_small_static_objs_to_verts(objs);
}

void building_room_geom_t::add_nested_objs_to_verts(vect_room_object_t const &objs_to_add) {
	vector_add_to(objs_to_add, pending_objs); // nested objects are added at the end so that small and text materials are thread safe
}
void building_room_geom_t::add_small_static_objs_to_verts(vect_room_object_t const &objs_to_add, bool inc_text) {
	if (objs_to_add.empty()) return; // don't add untextured material, otherwise we may fail the (num_verts > 0) assert
	float const tscale(2.0/obj_scale);

	for (unsigned i = 0; i < objs_to_add.size(); ++i) { // Note: iterating with indices to avoid invalid ref when add_nested_objs_to_verts() is called
		room_object_t const &c(objs_to_add[i]);
		if (!c.is_visible() || c.is_dynamic()) continue; // skip invisible and dynamic objects
		assert(c.is_strictly_normalized());
		assert(c.type < NUM_ROBJ_TYPES);

		switch (c.type) {
		case TYPE_BOOK:      add_book     (c, 0, 1, inc_text); break; // sm, maybe text
		case TYPE_BCASE:     add_bookcase (c, 0, 1, inc_text, tscale, 0); break; // sm, maybe text
		case TYPE_BED:       add_bed      (c, 0, 1, tscale); break;
		case TYPE_DESK:      add_desk     (c, tscale, 0, 1); break;
		case TYPE_DRESSER: case TYPE_NIGHTSTAND: add_dresser(c, tscale, 0, 1); break;
		case TYPE_TCAN:      add_trashcan (c); break;
		case TYPE_SIGN:      add_sign     (c, 0, inc_text); break; // sm, maybe text
		case TYPE_CLOSET:    add_closet   (c, tid_nm_pair_t(), 0, 1); break; // add closet wall trim and interior objects, don't need wall_tex
		case TYPE_RAILING:   add_railing  (c); break;
		case TYPE_PLANT:     add_potted_plant(c, 0, 1); break; // plant only
		case TYPE_CRATE:     add_crate    (c); break; // not small but only added to windowless rooms
		case TYPE_BOX:       add_box      (c); break; // not small but only added to windowless rooms
		case TYPE_SHELVES:   add_shelves  (c, tscale); break; // not small but only added to windowless rooms
		case TYPE_COMPUTER:  add_computer (c); break;
		case TYPE_KEYBOARD:  add_keyboard (c); break;
		case TYPE_WINE_RACK: add_wine_rack(c, 0, 1, tscale); break;
		case TYPE_BOTTLE:    add_bottle   (c); break;
		case TYPE_PAPER:     add_paper    (c); break;
		case TYPE_PAINTCAN:  add_paint_can(c); break;
		case TYPE_PEN: case TYPE_PENCIL: case TYPE_MARKER: add_pen_pencil_marker(c); break;
		case TYPE_LG_BALL:   add_lg_ball   (c); break;
		case TYPE_HANGER_ROD:add_hanger_rod(c); break;
		case TYPE_DRAIN:     add_drain_pipe(c); break;
		case TYPE_KEY:       if (has_key_3d_model()) {model_objs.push_back(c);} else {add_key(c);} break; // draw or add as 3D model
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
		case TYPE_BUTTON:    if (!(c.flags & RO_FLAG_IN_ELEV)) {add_button(c);} break; // skip buttons inside elevators, which are drawn as dynamic objects
		case TYPE_LBASKET:   add_laundry_basket(c); break;
		case TYPE_TOASTER:   add_toaster_proxy (c); break;
		case TYPE_WHEATER:   add_water_heater  (c); break;
		case TYPE_BRK_PANEL: add_breaker_panel (c); break; // only added to basements
		default: break;
		} // end switch
	} // for i
}

void building_room_geom_t::create_text_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Text");

	for (unsigned d = 0; d < 2; ++d) { // {objs, expanded_objs}
		vect_room_object_t const &v(d ? expanded_objs : objs);

		for (auto i = v.begin(); i != v.end(); ++i) {
			if (!i->is_visible()) continue; // skip invisible objects
			switch (i->type) {
			case TYPE_BOOK:  add_book    (*i, 0, 0, 1); break; // text only
			case TYPE_BCASE: add_bookcase(*i, 0, 0, 1, 1.0, 0); break; // text only
			case TYPE_SIGN:  add_sign    (*i, 0, 1); break; // text only
			default: break;
			} // end switch
		} // for i
	} // for d
}

void building_room_geom_t::create_detail_vbos(building_t const &building) {
	// currently only small objects that are non-interactive and can't be taken; TYPE_SWITCH almost counts; also, anything in the basement not seen from outside the building
	auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators
	tid_nm_pair_t const &wall_tex(building.get_material().wall_tex);

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_visible()) continue;

		switch (i->type) {
		case TYPE_OUTLET:     add_outlet(*i); break;
		case TYPE_VENT:       add_vent(*i); break;
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
		tmp_rgeom.add_book(c, 1, 1, 1, 0.0, 0, 0, z_rot_angle); // draw lg/sm/text
		enable_blend(); // needed for book text
		tmp_rgeom.mats_small.upload_draw_and_clear(s);
		tmp_rgeom.mats_text .upload_draw_and_clear(s);
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
	if (!has_room_geom()) return;
	if (ENABLE_MIRROR_REFLECTIONS && !shadow_only && !reflection_pass && player_in_building) {find_mirror_needing_reflection(xlate);}
	interior->room_geom->draw(bbd, s, *this, oc, xlate, building_ix, shadow_only, reflection_pass, inc_small, player_in_building);
}
void building_t::gen_and_draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, occlusion_checker_noncity_t &oc, vector3d const &xlate,
	unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building)
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
		gen_room_details(rgen, building_ix); // generate so that we can draw it
		assert(has_room_geom());
	}
	if (has_room_geom() && inc_small == 2) {add_wall_and_door_trim_if_needed();} // gen trim when close to the player
	draw_room_geom(bbd, s, oc, xlate, building_ix, shadow_only, reflection_pass, inc_small, player_in_building);
}
void building_t::clear_room_geom() {
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
	// could make mats[1] a subset of mats[0], but index logic is more complex
	rgeom_mat_t mats[2], web_mat; // mats: {high detail/67K verts, low detail shadow pass/4.2K verts}
	bool is_setup;

	void add_eye(rgeom_mat_t &mat, point const &pos, float radius) {
		cube_t eye_bc(pos);
		eye_bc.expand_by(radius);
		mat.add_sphere_to_verts(eye_bc, DK_RED, 1); // low_detail=1
	}
	void init() {
		// generate spider geometry; centered at (0,0,0) with radius=1.0; head is in +X
		colorRGBA const color(BLACK);
		float const body_zval(-0.3), leg_radius(0.03);
		point const abdomen_center(-0.8, 0.0, body_zval);
		cube_t abdomen(abdomen_center), body(point(0.0, 0.0, body_zval));
		abdomen.expand_by(vector3d(0.50, 0.35, 0.35));
		body   .expand_by(vector3d(0.45, 0.30, 0.20));

		for (unsigned n = 0; n < 2; ++n) { // {high detail, low detail shadow pass}
			bool const low_detail(n == 1);
			unsigned const ndiv(low_detail ? N_SPHERE_DIV/4 : N_SPHERE_DIV/2);
			unsigned cur_vert_pos(0);
			rgeom_mat_t &mat(mats[n]);
			mat.add_sphere_to_verts(abdomen, color, low_detail);
			mat.add_sphere_to_verts(body,    color, low_detail);

			if (!low_detail) { // red markings on back
				cube_t marking(abdomen_center + vector3d(-0.02, 0.0, 0.01));
				marking.expand_by(vector3d(0.48, 0.3, 0.35));
				mat.add_sphere_to_verts(marking, colorRGBA(0.5, 0.05, 0.0));
			}
			assign_tc_range(mat, cur_vert_pos, 0.0, 0.0, 0.0); // head and body aren't animated

			for (unsigned d = 0; d < 2; ++d) { // {left, right}
				float const d_sign(d ? -1.0 : 1.0);

				if (!low_detail) { // eyes and fangs are high detail only
					float const fang_radius(0.05);
					point const fang_top(0.44, 0.05*d_sign, body_zval-0.04), fang_bot(fang_top - vector3d(0.0, 0.0, 0.2));
					add_eye(mat, point(0.30, 0.080*d_sign, body_zval+0.14), 0.026);
					add_eye(mat, point(0.40, 0.045*d_sign, body_zval+0.08), 0.028);
					add_eye(mat, point(0.44, 0.020*d_sign, body_zval+0.04), 0.016);
					add_eye(mat, point(0.43, 0.055*d_sign, body_zval+0.03), 0.015);
					mat.add_sphere_to_verts(fang_top, vector3d(fang_radius, fang_radius, fang_radius), color, 1); // top of fang; low_detail=1
					mat.add_cylin_to_verts(fang_bot, fang_top, 0.0, fang_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv); // fang
					assign_tc_range(mat, cur_vert_pos, 0.0, 0.0, 0.0); // not animated
				}
				// add legs
				for (unsigned l = 0; l < 4; ++l) {
					float const ts(l/4.0);
					point const joint(0.12*(l - 1.5), 0.26*d_sign, body_zval);
					point const knee (2.0*joint.x, 2.0*joint.y,  0.5);
					point const ankle(2.8*knee .x, 2.8*knee .y,  0.0);
					point const foot (3.5*knee .x, 3.5*knee .y, -1.0);
					vector3d const sphere_radius(leg_radius, leg_radius, leg_radius);
					float const joint_tt(0.0*d_sign), knee_tt(0.3*d_sign), ankle_tt(0.7*d_sign), foot_tt(1.0*d_sign);
					if (!low_detail) {mat.add_sphere_to_verts(joint, sphere_radius, color, 1);} // round body joint; high detail only; low_detail=1
					assign_tc_range(mat, cur_vert_pos, ts, joint_tt, joint_tt);
					mat.add_cylin_to_verts(joint, knee, leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
					assign_tc_range(mat, cur_vert_pos, ts, joint_tt, knee_tt);
					if (!low_detail) {mat.add_sphere_to_verts(knee, sphere_radius, color, 1);} // round knee joint; high detail only; low_detail=1
					assign_tc_range(mat, cur_vert_pos, ts, knee_tt, knee_tt);
					mat.add_cylin_to_verts(ankle, knee, leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
					assign_tc_range(mat, cur_vert_pos, ts, ankle_tt, knee_tt);
					if (!low_detail) {mat.add_sphere_to_verts(ankle, sphere_radius, color, 1);} // round ankle joint; high detail only; low_detail=1
					assign_tc_range(mat, cur_vert_pos, ts, ankle_tt, ankle_tt);
					mat.add_cylin_to_verts(foot,  ankle, 0.1*leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
					assign_tc_range(mat, cur_vert_pos, ts, foot_tt, ankle_tt);
				} // for n
			} // for d
			mat.create_vbo_inner();
			mat.clear_vectors(1); // free_memory=1: vector data no longer needed
		} // for n
		is_setup = 1;
	}
	void assign_tc_range(rgeom_mat_t &mat, unsigned &cur_vert_pos, float ts, float tt_lo, float tt_hi) {
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
	spider_draw_t() : is_setup(0) {}
	void clear() {mats[0].clear(); mats[1].clear(); is_setup = 0;}

	void draw(vect_spider_t const &spiders, shader_t &s, building_t const &building, occlusion_checker_noncity_t &oc,
		vector3d const &xlate, bool shadow_only, bool reflection_pass, bool check_clip_cube)
	{
		if (spiders.empty()) return; // nothing to draw
		int anim_time_loc(-1);
		point const camera_bs(camera_pdu.pos - xlate);
		bool const check_occlusion(display_mode & 0x08);
		bool any_drawn(0);
		rgeom_mat_t &mat(mats[shadow_only]);

		for (spider_t const &S : spiders) { // future work: use instancing
			if (!shadow_only && S.on_web) { // draw spider webs
				cube_t strand(S.pos - 0.3*S.radius*S.upv); // end of abdomen (matches body_zval)
				set_cube_zvals(strand, (S.pos.z + 1.3*S.radius), S.web_start_zval);
				strand.expand_by_xy(0.02*S.radius);

				if (strand.z1() < strand.z2() && camera_pdu.cube_visible(strand + xlate)) {
					web_mat.add_vcylin_to_verts(strand, WHITE, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 16); // ndiv=16
				}
			}
			cube_t const bcube(S.get_bcube());
			if (check_clip_cube && !smap_light_clip_cube.intersects(bcube + xlate)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
			if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
			if (check_occlusion && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass)) continue;

			if (!any_drawn) { // setup shaders
				if (!is_setup) {init();}
				mat.vao_setup(shadow_only);
				s.set_specular(0.5, 80.0);
				select_texture(WHITE_TEX);
				int const animation_id = 8; // custom spider animation
				s.add_uniform_int("animation_id", animation_id);
				s.add_uniform_float("animation_scale",    1.0); // not using a model, nominal size is 1.0
				s.add_uniform_float("model_delta_height", 1.0); // not using a model, nominal size is 1.0
				s.add_uniform_float("bump_map_mag",       0.0);
				mat.pre_draw(shadow_only);
				anim_time_loc = s.get_uniform_loc("animation_time");
				any_drawn = 1;
			}
			s.set_uniform_float(anim_time_loc, S.anim_time);
			fgPushMatrix();
			translate_to(S.pos);
			vector3d const dir(S.dir.get_norm()), upv(S.upv.get_norm()); // normalize, just in case
			vector3d const right(cross_product(S.upv, S.dir).get_norm());
			xform_matrix xm;
			float *m(xm.get_ptr());
			m[0] = dir.x; // x
			m[1] = dir.y;
			m[2] = dir.z;
			m[4] = right.x; // y
			m[5] = right.y;
			m[6] = right.z;
			m[8] = upv.x; // z
			m[9] = upv.y;
			m[10]= upv.z;
			fgMultMatrix(xm);
			uniform_scale(S.radius);
			check_mvm_update();
			// Note: we could use hardware instancing here if the MVM is used as a uniform in the shader, but we also need to set per-instance anim_time
			mat.draw_geom();
			fgPopMatrix();
		} // for S
		if (any_drawn) { // reset state
			check_mvm_update(); // make sure to reset MVM
			s.add_uniform_float("bump_map_mag",   1.0);
			s.add_uniform_float("animation_time", 0.0); // reset animation time
			s.add_uniform_int  ("animation_id",   0); // clear animation
			s.clear_specular();
			indexed_vao_manager_with_shadow_t::post_render();
		}
		if (!web_mat.empty()) {
			select_texture(WHITE_TEX);
			tid_nm_pair_dstate_t state(s);
			web_mat.upload_draw_and_clear(state);
		}
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
	if (frame_counter < 100 || frame_counter > last_frame) {num_geom_this_frame = 0; last_frame = frame_counter;} // unlimited for the first 100 frames
	point const camera_bs(camera_pdu.pos - xlate);
	bool const draw_lights(camera_bs.z < building.bcube.z2()); // don't draw ceiling lights when player is above the building
	bool const draw_detail_objs(inc_small == 2 && (!shadow_only || building.has_parking_garage)); // only parking garages have detail objects that cast shadows
	if (player_in_building) {bbd = nullptr;} // use immediate drawing when player is in the building because draw order matters for alpha blending
	if (bbd != nullptr) {bbd->set_camera_dir_mask(camera_bs, building.bcube);}

	if (has_pictures && num_pic_tids != num_screenshot_tids) {
		invalidate_static_geom(); // user created a new screenshot texture, and this building has pictures - recreate room geom
		num_pic_tids = num_screenshot_tids;
	}
	check_invalid_draw_data();
	// generate vertex data in the shadow pass or if we haven't hit our generation limit; must be consistent for static and small geom
	// Note that the distance cutoff for mats_static and mats_small is different, so we generally won't be creating them both
	// unless the player just appeared by this building, or we need to update the geometry; in either case this is higher priority and we want to update both
	if (shadow_only || num_geom_this_frame < MAX_ROOM_GEOM_GEN_PER_FRAME) {
		if (!mats_static.valid) { // create static materials if needed
			create_obj_model_insts(building);
			create_static_vbos(building);
			if (!shadow_only) {++num_geom_this_frame;}
		}
		bool const create_small(inc_small && !mats_small.valid), create_text(inc_small == 2 && !mats_text.valid);
		//highres_timer_t timer("Create Small + Text VBOs", (create_small || create_text));

		if (create_small && create_text) {
#pragma omp parallel num_threads(2)
			if (omp_get_thread_num_3dw() == 0) {create_small_static_vbos(building);} else {create_text_vbos(building);}
		}
		else if (create_small) {create_small_static_vbos(building);}
		else if (create_text ) {create_text_vbos        (building);}
		add_small_static_objs_to_verts(pending_objs, create_text);
		pending_objs.clear();

		// upload VBO data serially
		if (create_small) {
			mats_small.create_vbos(building);
			mats_amask.create_vbos(building);
		}
		if (create_text) {mats_text.create_vbos(building);}
		if (!shadow_only) {num_geom_this_frame += (unsigned(create_small) + unsigned(create_text));}

		// Note: not created on the shadow pass, because trim_objs may not have been created yet and we would miss including it
		if (draw_detail_objs && !shadow_only && !mats_detail.valid) { // create detail materials if needed
			create_detail_vbos(building);
			if (!shadow_only) {++num_geom_this_frame;}
		}
	}
	if (draw_lights && !mats_lights .valid) {create_lights_vbos (building);} // create lights  materials if needed (no limit)
	if (inc_small   && !mats_dynamic.valid) {create_dynamic_vbos(building);} // create dynamic materials if needed (no limit); drawn with small objects
	if (!mats_doors.valid) {create_door_vbos(building);} // create door materials if needed (no limit)
	enable_blend(); // needed for rugs and book text
	assert(s.is_setup());
	mats_static.draw(bbd, s, shadow_only, reflection_pass); // this is the slowest call
	if (draw_lights)      {mats_lights .draw(bbd, s, shadow_only, reflection_pass);}
	if (inc_small  )      {mats_dynamic.draw(bbd, s, shadow_only, reflection_pass);}
	if (draw_detail_objs) {mats_detail .draw(bbd, s, shadow_only, reflection_pass);}
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
	if (inc_small == 2 && !shadow_only) {mats_text.draw(bbd, s, shadow_only, reflection_pass);} // text must be drawn last; drawn as detail objects
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
		bool const is_emissive(!shadow_only && obj.type == TYPE_LAMP && obj.is_light_on());
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
			for (rat_t &rat : rats) {
				cube_t const bcube(rat.get_bcube());
				if (check_clip_cube && !smap_light_clip_cube.intersects(bcube + xlate)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
				if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
				if (check_occlusion && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass)) continue;
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
			model3d::bind_default_flat_normal_map();
		} // end rats drawing
		spider_draw.draw(spiders, s, building, oc, xlate, shadow_only, reflection_pass, check_clip_cube);
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
	decal_manager.draw_building_interior_decals(s, player_in_building, shadow_only); // draw decals in this building
	
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
pair<cube_t, colorRGBA> car_bcube_color_from_parking_space(room_object_t const &o) {
	car_t car(car_from_parking_space(o));
	set_car_model_color(car);
	return make_pair(car.bcube, car.get_color());
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
	if (!has_room_geom()) return; // can get here in rare cases, maybe only shadow_only pass
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
		if (obj.is_used() && (!player_in_this_building || !check_occlusion || !has_int_garage ||
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
			if (i->type != TYPE_PARK_SPACE || !i->is_used()) continue; // not a space, or no car in this space
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


