// 3D World - Building Interior Room Item Drawing
// by Frank Gennari 4/17/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#include "subdiv.h" // for sd_sphere_d
#include "profiler.h"

unsigned const MAX_ROOM_GEOM_GEN_PER_FRAME = 1;

object_model_loader_t building_obj_model_loader;

extern bool camera_in_building;
extern int display_mode, frame_counter, animate2;
extern float office_chair_rot_rate;
extern point pre_smap_player_pos;
extern pos_dir_up camera_pdu;
extern building_t const *player_building;
extern room_object_t player_held_object;

unsigned get_num_screenshot_tids();

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
			v.set_ortho_norm(n, j);
			if (inverted) {v.invert_normal();}
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
void rgeom_mat_t::add_cube_to_verts_untextured(cube_t const &c, colorRGBA const &color, unsigned skip_faces) {
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

void rgeom_mat_t::add_sphere_to_verts(cube_t const &c, colorRGBA const &color, bool low_detail, vector3d const &skip_hemi_dir, xform_matrix const *const matrix) {
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
	point const center(c.get_cube_center()), size(0.5*c.get_size());
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

class rgeom_alloc_t {
	deque<rgeom_storage_t> free_list; // one per unique texture ID/material
public:
	void alloc(rgeom_storage_t &s) { // attempt to use free_list entry to reuse existing capacity
		if (free_list.empty()) return; // no pre-alloc
		//cout << TXT(free_list.size()) << TXT(free_list.back().get_tot_vert_capacity()) << endl;

		// try to find a free list element with the same tex so that we balance out material memory usage/capacity better
		for (unsigned i = 0; i < free_list.size(); ++i) {
			if (free_list[i].tex.tid != s.tex.tid) continue;
			s.swap_vectors(free_list[i]); // transfer existing capacity from free list
			free_list[i].swap(free_list.back());
			free_list.pop_back();
			return; // done
		}
		//s.swap(free_list.back());
		//free_list.pop_back();
	}
	void free(rgeom_storage_t &s) {
		s.clear(); // in case the caller didn't clear it
		free_list.push_back(rgeom_storage_t(s.tex)); // record tex of incoming element
		s.swap_vectors(free_list.back()); // transfer existing capacity to free list; clear capacity from s
	}
};

rgeom_alloc_t rgeom_alloc; // static allocator with free list, shared across all buildings; not thread safe

void rgeom_storage_t::clear() {
	quad_verts.clear();
	itri_verts.clear();
	indices.clear();
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
	vbo.clear();
	delete_and_zero_vbo(ivbo);
	rgeom_storage_t::clear();
	num_qverts = num_itverts = num_ixs = 0;
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
	num_qverts  = quad_verts.size();
	num_itverts = itri_verts.size();
	num_ixs     = indices.size();
	unsigned qsz(num_qverts*sizeof(vertex_t)), itsz(num_itverts*sizeof(vertex_t));
	vbo.vbo = ::create_vbo();
	check_bind_vbo(vbo.vbo);
	upload_vbo_data(nullptr, get_tot_vert_count()*sizeof(vertex_t));
	upload_vbo_sub_data(quad_verts.data(), 0, qsz);
	upload_vbo_sub_data(itri_verts.data(), qsz, itsz);
	bind_vbo(0);

	if (!indices.empty()) { // we have some indexed quads
		for (auto i = indices.begin(); i != indices.end(); ++i) {*i += num_qverts;} // shift indices to match the new vertex location
		create_vbo_and_upload(ivbo, indices, 1, 1);
	}
}

void rgeom_mat_t::draw(shader_t &s, bool shadow_only, bool reflection_pass) {
	if (shadow_only && !en_shadows)  return; // shadows not enabled for this material (picture, whiteboard, rug, etc.)
	if (shadow_only && tex.emissive) return; // assume this is a light source and shouldn't produce shadows (also applies to bathroom windows, which don't produce shadows)
	if (reflection_pass && tex.tid == REFLECTION_TEXTURE_ID) return; // don't draw reflections of mirrors as this doesn't work correctly
	assert(vbo.vbo_valid());
	assert(num_qverts > 0 || num_itverts > 0);
	if (!shadow_only) {tex.set_gl(s);} // ignores texture scale for now
	vbo.pre_render();
	vertex_t::set_vbo_arrays();
	if (num_qverts > 0) {draw_quads_as_tris(num_qverts);}

	if (num_itverts > 0) { // index quads, used for cylinders and spheres
		assert(ivbo > 0);
		bind_vbo(ivbo, 1);
		//glDisable(GL_CULL_FACE); // two sided lighting requires fewer verts (no duplicates), but must be set in the shader
		glDrawRangeElements(GL_TRIANGLES, num_qverts, (num_qverts + num_itverts), num_ixs, GL_UNSIGNED_INT, nullptr);
		//glEnable(GL_CULL_FACE);
		bind_vbo(0, 1);
	}
	if (!shadow_only) {tex.unset_gl(s);}
}

void rgeom_mat_t::upload_draw_and_clear(shader_t &s) {
	if (empty()) return; // nothing to do
	create_vbo_inner();
	draw(s, 0, 0);
	clear();
}

void building_materials_t::clear() {
	for (iterator m = begin(); m != end(); ++m) {m->clear();}
	vector<rgeom_mat_t>::clear();
}
unsigned building_materials_t::count_all_verts() const {
	unsigned num_verts(0);
	for (const_iterator m = begin(); m != end(); ++m) {num_verts += m->get_tot_vert_count();}
	return num_verts;
}
rgeom_mat_t &building_materials_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows) {
	// for now we do a simple linear search because there shouldn't be too many unique materials
	for (iterator m = begin(); m != end(); ++m) {
		if (m->tex != tex) continue;
		if (inc_shadows) {m->enable_shadows();}
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
void building_materials_t::draw(shader_t &s, bool shadow_only, bool reflection_pass) {
	static vector<iterator> text_mats;
	text_mats.clear();

	// first pass, draw regular materials (excluding text)
	for (iterator m = begin(); m != end(); ++m) {
		if (m->tex.tid == FONT_TEXTURE_ID) {text_mats.push_back(m);} // skip in this pass
		else {m->draw(s, shadow_only, reflection_pass);}
	}
	// second pass, draw text (if it exists) so that alpha blending works
	for (auto m = text_mats.begin(); m != text_mats.end(); ++m) {(*m)->draw(s, shadow_only, reflection_pass);}
}
void building_materials_t::upload_draw_and_clear(shader_t &s) {
	for (iterator m = begin(); m != end(); ++m) {m->upload_draw_and_clear(s);}
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
void building_room_geom_t::clear_materials() { // can be called to update textures, lighting state, etc.
	clear_static_vbos();
	clear_static_small_vbos();
	mats_dynamic.clear();
	mats_lights.clear();
	mats_doors.clear();
}
void building_room_geom_t::clear_static_vbos() { // used to clear pictures
	mats_static.clear();
	obj_model_insts.clear(); // these are associated with static VBOs
	mats_alpha.clear();
}
void building_room_geom_t::clear_static_small_vbos() {
	mats_small.clear();
	mats_plants.clear();
}

rgeom_mat_t &building_room_geom_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows, bool dynamic, bool small, bool transparent) {
	return (dynamic ? mats_dynamic : (small ? mats_small : (transparent ? mats_alpha : mats_static))).get_material(tex, inc_shadows);
}
rgeom_mat_t &building_room_geom_t::get_metal_material(bool inc_shadows, bool dynamic, bool small) {
	tid_nm_pair_t tex(-1, 1.0, inc_shadows);
	tex.set_specular(0.8, 60.0);
	return get_material(tex, inc_shadows, dynamic, small);
}

void room_object_t::set_as_bottle(unsigned rand_id, unsigned max_type) {
	assert(max_type > 0 && max_type < NUM_BOTTLE_TYPES);
	obj_id = (uint16_t)rand_id;
	while (get_bottle_type() > max_type) {obj_id += 13;} // cycle with a prime number until a valid type is selected
	color  = bottle_params[get_bottle_type()].color;
}

void building_room_geom_t::create_static_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom"); // 2.35ms
	float const tscale(2.0/obj_scale);
	mats_static.clear();
	mats_alpha .clear();

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible() || i->is_dynamic()) continue; // skip invisible and dynamic objects
		assert(i->is_strictly_normalized());
		assert(i->type < NUM_ROBJ_TYPES);

		switch (i->type) {
		case TYPE_TABLE:   add_table   (*i, tscale, 0.12, 0.08); break; // top_dz=12% of height, leg_width=8% of height
		case TYPE_CHAIR:   add_chair   (*i, tscale); break;
		case TYPE_STAIR:   add_stair   (*i, tscale, tex_origin); break;
		case TYPE_STAIR_WALL: add_stairs_wall(*i, tex_origin, building.get_material().wall_tex); break;
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
		case TYPE_CLOSET:  add_closet  (*i, building.get_material().wall_tex, 1, 0); break;
		case TYPE_MIRROR:  add_mirror  (*i); break;
		case TYPE_SHOWER:  add_shower  (*i, tscale); break;
		case TYPE_MWAVE:   add_mwave   (*i); break;
		case TYPE_BLINDS:  add_blinds  (*i); break;
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
}

void building_room_geom_t::create_small_static_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Small"); // 7.8ms, slow building at 26,16
	mats_small .clear();
	mats_plants.clear();
	model_objs .clear(); // currently model_objs are only created for small objects in drawers, so we clear this here
	add_small_static_objs_to_verts(expanded_objs);
	add_small_static_objs_to_verts(objs);
	mats_small .create_vbos(building);
	mats_plants.create_vbos(building);
}

void building_room_geom_t::add_small_static_objs_to_verts(vector<room_object_t> const &objs_to_add) {
	float const tscale(2.0/obj_scale);
	get_untextured_material(0, 0, 1); // must ensure book covers are drawn before the text so that alpha blending works properly, so add the untextured cover material first

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
		case TYPE_TCAN:      add_trashcan(*i); break;
		case TYPE_SIGN:      add_sign     (*i, 0, 1); break;
		case TYPE_WALL_TRIM: add_wall_trim(*i); break;
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
		case TYPE_KEY: if (has_key_3d_model()) {model_objs.push_back(*i);} else {add_key(*i);} break; // draw or add as 3D model
		case TYPE_MONEY:     add_money(*i); break;
		case TYPE_PHONE:     add_phone(*i); break;
		case TYPE_TPROLL:    add_tproll(*i); break;
		case TYPE_SPRAYCAN:  add_spraycan(*i); break;
		case TYPE_CRACK:     add_crack(*i); break;
		case TYPE_BUTTON:    if (!(i->flags & RO_FLAG_IN_ELEV)) {add_button(*i);} break; // skip buttons inside elevators, which are drawn as dynamic objects
		default: break;
		}
	} // for i
}

void building_room_geom_t::create_obj_model_insts(building_t const &building) { // handle drawing of 3D models
	//highres_timer_t timer("Gen Room Model Insts");
	obj_model_insts.clear();

	for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
		auto const &obj_vect((vect_id == 1) ? expanded_objs : objs);
		unsigned const obj_id_offset((vect_id == 1) ? objs.size() : 0);
		auto objs_end((vect_id == 1) ? expanded_objs.end() : get_std_objs_end()); // skip buttons/stairs/elevators

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
	//highres_timer_t timer("Gen Room Geom Light"); // 0.3ms
	float const tscale(2.0/obj_scale);
	auto objs_end(get_std_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->is_visible() && i->type == TYPE_LIGHT) {add_light(*i, tscale);}
	}
	mats_lights.create_vbos(building);
}

void building_room_geom_t::create_dynamic_vbos(building_t const &building) {
	//highres_timer_t timer(string("Gen Room Geom Dynamic ") + (building.is_house ? "house" : "office"));
	
	if (building.is_house) { // currently, only houses have dynamic objects (balls)
		auto objs_end(get_std_objs_end()); // skip buttons/stairs/elevators

		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (!i->is_visible() || !i->is_dynamic()) continue; // only visible + dynamic objects; can't do VFC because this is not updated every frame
			switch (i->type) {
			case TYPE_LG_BALL: add_lg_ball(*i); break;
			default: assert(0); // not a supported dynamic object type
			}
		} // for i
	}
	for (auto e = building.interior->elevators.begin(); e != building.interior->elevators.end(); ++e) {
		assert(e->car_obj_id < objs.size());
		assert(e->button_id_start < e->button_id_end && e->button_id_end <= objs.size());
		add_elevator_doors(*e); // add dynamic elevator doors
		add_elevator(objs[e->car_obj_id], 2.0/obj_scale); // draw elevator car for this elevator

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

	for (auto i = doors.begin(); i != doors.end(); ++i) {
		building.add_door_verts(*i, *this, tquad_with_ix_t::TYPE_HDOOR, i->dim, i->open_dir, i->open, 0, 0, i->on_stairs); // opens_out=0, exterior=0
	}
	mats_doors.create_vbos(building);
}

void building_room_geom_t::expand_object(room_object_t &c) {
	if (c.flags & RO_FLAG_EXPANDED) return; // already expanded
	switch (c.type) {
	case TYPE_CLOSET:    expand_closet   (c); break;
	case TYPE_SHELVES:   expand_shelves  (c); break;
	case TYPE_WINE_RACK: expand_wine_rack(c); break;
	default: assert(0); // not a supported expand type
	}
	c.flags |= RO_FLAG_EXPANDED; // flag as expanded
}

void rotate_dir_about_z(vector3d &dir, float rate) { // Note: assumes dir is normalized
	if (rate == 0.0) return;
	assert(dir.z == 0.0); // dir must be in XY plane
	float const new_angle(atan2(dir.y, dir.x) + rate*fticks);
	dir.assign(cosf(new_angle), sinf(new_angle), 0.0);
}
void apply_room_obj_rotate(room_object_t &obj, obj_model_inst_t &inst) {
	if (!(obj.flags & RO_FLAG_ROTATING)) return;
	if (office_chair_rot_rate == 0.0) {obj.flags &= ~RO_FLAG_ROTATING; return;} // if no longer rotating, clear rotation bit
	assert(obj.type == TYPE_OFF_CHAIR); // only office chairs are supported for now
	rotate_dir_about_z(inst.dir, office_chair_rot_rate);
}

/*static*/ void building_room_geom_t::draw_interactive_player_obj(room_object_t const &c, shader_t &s) {
	static rgeom_mat_t mat = rgeom_mat_t(tid_nm_pair_t()); // allocated memory is reused across frames; VBO is recreated every time

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
	else if (c.type == TYPE_TPROLL) { // apply get_player_cview_rot_matrix()?
		add_vert_tproll_to_material(c, mat);
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
	else {assert(0);}
	if (c.type == TYPE_TPROLL) {enable_blend();}
	mat.upload_draw_and_clear(s);
	if (c.type == TYPE_TPROLL) {disable_blend();}
}

class water_draw_t {
	rgeom_mat_t mat;
	float tex_off;
public:
	water_draw_t() : mat(rgeom_mat_t(tid_nm_pair_t(FOAM_TEX))), tex_off(0.0) {}

	void add_water_for_sink(room_object_t const &obj) {
		if (!(obj.flags & RO_FLAG_IS_ACTIVE)) return; // not turned on
		bool const is_cube(obj.type == TYPE_KSINK || obj.type == TYPE_BRSINK);
		float const dz(obj.dz());
		cube_t c;
		point pos(obj.get_cube_center());
		pos[obj.dim] += (obj.dir ? -1.0 : 1.0)*(is_cube ? 0.25 : 0.095)*obj.get_sz_dim(obj.dim); // move toward the back of the sink
		c.set_from_sphere(pos, (is_cube ? 0.02 : 0.0055)*dz);
		set_cube_zvals(c, (obj.z1() + (is_cube ? 0.7 : 0.6)*dz), (obj.z1() + (is_cube ? 1.3 : 0.925)*dz));
		unsigned const verts_start(mat.itri_verts.size());
		mat.add_vcylin_to_verts(c, colorRGBA(WHITE, 0.5), 0, 0, 0, 0, 1.0, 1.0, 0.2);
		for (auto i = mat.itri_verts.begin() + verts_start; i != mat.itri_verts.end(); ++i) {i->t[1] *= 1.2; i->t[1] += tex_off;}
	}
	void draw_and_clear(shader_t &s) {
		if (mat.empty()) return;
		enable_blend();
		mat.upload_draw_and_clear(s);
		disable_blend();
		if (animate2) {tex_off += 0.02*fticks;} // animate the texture
		if (tex_off > 1.0) {tex_off -= 1.0;}
	}
};

water_draw_t water_draw;

// Note: non-const because it creates the VBO
void building_room_geom_t::draw(shader_t &s, building_t const &building, occlusion_checker_noncity_t &oc, vector3d const &xlate,
	unsigned building_ix, bool shadow_only, bool reflection_pass, bool inc_small, bool player_in_building)
{
	if (empty()) return; // no geom
	unsigned const num_screenshot_tids(get_num_screenshot_tids());
	static int last_frame(0);
	static unsigned num_geom_this_frame(0); // used to limit per-frame geom gen time; doesn't apply to shadow pass, in case shadows are cached
	if (frame_counter > last_frame) {num_geom_this_frame = 0; last_frame = frame_counter;}
	bool const can_update_geom(shadow_only || num_geom_this_frame < MAX_ROOM_GEOM_GEN_PER_FRAME); // must be consistent for static and small geom

	if (materials_invalid) { // set in set_obj_lit_state_to()
		clear_materials();
		materials_invalid = 0;
	}
	if (lights_changed) {
		mats_lights.clear();
		lights_changed = 0;
	}
	if (has_pictures && num_pic_tids != num_screenshot_tids) {
		clear_static_vbos(); // user created a new screenshot texture, and this building has pictures - recreate room geom
		num_pic_tids = num_screenshot_tids;
	}
	if (mats_static.empty() && can_update_geom) { // create static materials if needed
		create_obj_model_insts(building);
		create_static_vbos(building);
		++num_geom_this_frame;
	}
	if (inc_small && mats_small.empty() && can_update_geom) { // create small materials if needed
		create_small_static_vbos(building);
		++num_geom_this_frame;
	}
	if (mats_lights .empty()) {create_lights_vbos (building);} // create lights  materials if needed (no limit)
	if (mats_dynamic.empty()) {create_dynamic_vbos(building);} // create dynamic materials if needed (no limit)
	if (mats_doors  .empty()) {create_door_vbos   (building);} // create door    materials if needed (no limit)
	enable_blend(); // needed for rugs and book text
	assert(s.is_setup());
	mats_static .draw(s, shadow_only, reflection_pass);
	mats_lights .draw(s, shadow_only, reflection_pass);
	mats_dynamic.draw(s, shadow_only, reflection_pass);
	mats_doors  .draw(s, shadow_only, reflection_pass);

	if (inc_small) {
		mats_small.draw(s, shadow_only, reflection_pass);

		if (player_in_building) { // if we're not in the building, don't draw plants at all; without the special shader they won't look correct when drawn through windows
			if (shadow_only) {
				shader_t plant_shader;
				plant_shader.begin_simple_textured_shader(0.9); // need to use texture with alpha test
				mats_plants.draw(s, 0, 0);
				s.make_current(); // switch back to the normal shader
			}
			else if (reflection_pass) {
				mats_plants.draw(s, 0, 1);
			}
			else { // this is expensive: only enable for the current building and the main draw pass
				shader_t plant_shader;
				setup_building_draw_shader(plant_shader, 0.9, 1, 1, 0); // min_alpha=0.5, enable_indir=1, force_tsl=1, use_texgen=1
				mats_plants.draw(plant_shader, 0, 0);
				s.make_current(); // switch back to the normal shader
			}
		}
	}
	disable_blend();
	vbo_wrap_t::post_render();
	bool const disable_cull_face(0); // better but slower?
	if (disable_cull_face) {glDisable(GL_CULL_FACE);}
	point const camera_bs(camera_pdu.pos - xlate), building_center(building.bcube.get_cube_center());
	bool const is_rotated(building.is_rotated());
	oc.set_exclude_bix(building_ix);
	bool obj_drawn(0);
	water_sound_manager_t water_sound_manager(camera_bs);

	// draw object models
	for (auto i = obj_model_insts.begin(); i != obj_model_insts.end(); ++i) {
		room_object_t &obj(get_room_object_by_index(i->obj_id));
		if (!player_in_building && obj.is_interior()) continue; // don't draw objects in interior rooms if the player is outside the building (useful for office bathrooms)
		point obj_center(obj.get_cube_center());
		if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
		if (!shadow_only && !dist_less_than(camera_bs, obj_center, 100.0*obj.dz())) continue; // too far away (obj.max_len()?)
		bool const is_sink(player_in_building && !shadow_only && obj.type == TYPE_SINK);
		if (is_sink) {water_sound_manager.register_running_water(obj, building);}
		if (!(is_rotated ? building.is_rot_cube_visible(obj, xlate) : camera_pdu.cube_visible(obj + xlate))) continue; // VFC
		if ((display_mode & 0x08) && building.check_obj_occluded(obj, camera_bs, oc, reflection_pass)) continue;
		bool const is_emissive(obj.type == TYPE_LAMP && obj.is_lit());
		if (is_emissive) {s.set_color_e(LAMP_COLOR*0.4);}
		apply_room_obj_rotate(obj, *i); // Note: may modify obj by clearing flags
		building_obj_model_loader.draw_model(s, obj_center, obj, i->dir, obj.color, xlate, obj.get_model_id(), shadow_only, 0, 0);
		if (is_emissive) {s.set_color_e(BLACK);}
		if (is_sink) {water_draw.add_water_for_sink(obj);}
		obj_drawn = 1;
	} // for i
	if (player_in_building && !shadow_only && !reflection_pass) { // these models aren't drawn in the shadow or reflection passes; no emissive or rotated objects
		for (auto i = model_objs.begin(); i != model_objs.end(); ++i) {
			point obj_center(i->get_cube_center());
			if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
			if (!shadow_only && !dist_less_than(camera_bs, obj_center, 100.0*i->max_len())) continue; // too far away
			if (!(is_rotated ? building.is_rot_cube_visible(*i, xlate) : camera_pdu.cube_visible(*i + xlate))) continue; // VFC
			if ((display_mode & 0x08) && building.check_obj_occluded(*i, camera_bs, oc, reflection_pass)) continue;
			vector3d dir(i->get_dir());
			if (is_rotated) {building.do_xy_rotate_normal(dir);}
			building_obj_model_loader.draw_model(s, obj_center, *i, dir, i->color, xlate, i->get_model_id(), shadow_only, 0, 0);
			obj_drawn = 1;
		} // for i
	}
	if (disable_cull_face) {glEnable(GL_CULL_FACE);}
	if (obj_drawn) {check_mvm_update();} // needed after popping model transform matrix

	if (player_in_building && !shadow_only) { // draw water for sinks that are turned on
		auto objs_end(get_std_objs_end()); // skip buttons/stairs/elevators

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
		else if (player_held_object.can_use()) {draw_interactive_player_obj(player_held_object, s);}
		else {assert(0);}
	}
	// alpha blended, should be drawn near last
	if (!shadow_only) {draw_building_interior_paint(3, &building);} // draw both interior and exterior
	if (player_in_building && !shadow_only) {draw_building_interior_decals(&building);} // for blood decals, only drawn for current building
	
	if (!shadow_only && !mats_alpha.empty()) { // draw last; not shadow casters; for shower glass, etc.
		enable_blend();
		glDepthMask(GL_FALSE); // disable depth writing
		mats_alpha.draw(s, shadow_only, reflection_pass);
		glDepthMask(GL_TRUE);
		disable_blend();
		vbo_wrap_t::post_render();
	}
}

bool are_pts_occluded_by_any_cubes(point const &pt, point const *const pts, unsigned npts, vect_cube_t const &cubes, unsigned dim) {
	assert(npts > 0 && dim <= 2);

	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if ((pt[dim] < c->d[dim][0]) == (pts[0][dim] < c->d[dim][0])) continue; // skip if cube face does not separate pt from the first point
		if (!check_line_clip(pt, pts[0], c->d)) continue; // first point does not intersect
		bool not_occluded(0);

		for (unsigned p = 1; p < npts; ++p) { // skip first point
			if (!check_line_clip(pt, pts[p], c->d)) {not_occluded = 1; break;}
		}
		if (!not_occluded) return 1;
	} // for c
	return 0;
}

bool building_t::check_obj_occluded(cube_t const &c, point const &viewer_in, occlusion_checker_noncity_t &oc, bool reflection_pass) const {
	if (!interior) return 0; // could probably make this an assert
	//highres_timer_t timer("Check Object Occlusion"); // 0.001ms
	point viewer(viewer_in);
	maybe_inv_rotate_point(viewer); // rotate viewer pos into building space
	if (c.z2() < ground_floor_z1 && !bcube.contains_pt(viewer)) return 1; // fully inside basement, and viewer outside building: will be occluded
	point pts[8];
	unsigned const npts(get_cube_corners(c.d, pts, viewer, 0)); // should return only the 6 visible corners
	
	if (!reflection_pass) { // check walls of this building; not valid for reflections because the reflected camera may be on the other side of a wall/mirror
		for (unsigned d = 0; d < 2; ++d) {
			if (are_pts_occluded_by_any_cubes(viewer, pts, npts, interior->walls[d], d)) return 1;
		}
	}
	if (reflection_pass || bcube.contains_pt(viewer)) { // viewer inside this building; includes shadow_only case and reflection_pass (even if reflected camera is outside the building)
		// check floors of this building (and technically also ceilings)
		if (fabs(viewer.z - c.zc()) > (reflection_pass ? 1.0 : 0.5)*get_window_vspace()) { // on different floors
			if (are_pts_occluded_by_any_cubes(viewer, pts, npts, interior->floors, 2)) return 1;
		}
	}
	else if (camera_in_building) { // player in some other building
		if (is_rotated()) return 0; // not implemented yet - need to rotate viewer and pts into coordinate space of player_building

		if (player_building != nullptr && player_building->interior) { // check walls of the building the player is in
			assert(player_building != this); // otherwise player_in_this_building should be true
			
			for (unsigned d = 0; d < 2; ++d) { // check walls of the building the player is in
				if (are_pts_occluded_by_any_cubes(viewer, pts, npts, player_building->interior->walls[d], d)) return 1;
			}
		}
	}
	else { // player not in a building
		if (is_rotated()) return 0; // not implemented yet - c is not an axis aligned cube in global coordinate space
		if (oc.is_occluded(c)) return 1; // check other buildings
	}
	return 0;
}

