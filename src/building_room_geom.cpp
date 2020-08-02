// 3D World - Building Interior Room Geometry Drawing
// by Frank Gennari 7/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#include "subdiv.h" // for sd_sphere_d
#include "profiler.h"
#pragma warning(disable : 26812) // prefer enum class over enum

bool const ADD_BOOK_COVERS = 1;
bool const ADD_BOOK_TITLES = 1;
unsigned const MAX_ROOM_GEOM_GEN_PER_FRAME = 1;
colorRGBA const WOOD_COLOR(0.9, 0.7, 0.5); // light brown, multiplies wood texture color

object_model_loader_t building_obj_model_loader;

extern int display_mode, frame_counter;
extern pos_dir_up camera_pdu;

int get_rand_screenshot_texture(unsigned rand_ix);
unsigned get_num_screenshot_tids();

void gen_text_verts(vector<vert_tc_t> &verts, point const &pos, string const &text, float tsize, vector3d const &column_dir, vector3d const &line_dir, bool use_quads=0);
string const &gen_book_title(unsigned rand_id, string *author, unsigned split_len);


unsigned get_face_mask(unsigned dim, bool dir) {return ~(1 << (2*(2-dim) + dir));} // skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2
unsigned get_skip_mask_for_xy(bool dim) {return (dim ? EF_Y12 : EF_X12);}

// skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2 to match CSG cube flags
void rgeom_mat_t::add_cube_to_verts(cube_t const &c, colorRGBA const &color, vector3d const &tex_origin,
	unsigned skip_faces, bool swap_tex_st, bool mirror_x, bool mirror_y, bool inverted)
{
	//assert(c.is_normalized()); // no, bathroom window is denormalized
	vertex_t v;
	v.set_c4(color);

	// Note: stolen from draw_cube() with tex coord logic, back face culling, etc. removed
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
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

void rgeom_mat_t::add_vcylin_to_verts(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top, bool two_sided, bool ts_tb, bool inv_tb, float rs_bot, float rs_top) {
	assert(!(ts_tb && inv_tb));
	point const center(c.get_cube_center());
	point const ce[2] = {point(center.x, center.y, c.z1()), point(center.x, center.y, c.z2())};
	unsigned const ndiv(N_CYL_SIDES), num_ends((unsigned)draw_top + (unsigned)draw_bot);
	float const radius(0.5*min(c.dx(), c.dy())), ndiv_inv(1.0/ndiv); // cube X/Y size should be equal/square
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius*rs_bot, radius*rs_top, ndiv, v12));
	color_wrapper const cw(color);
	unsigned itris_start(itri_verts.size()), ixs_start(indices.size()), itix(itris_start), iix(ixs_start);
	itri_verts.resize(itris_start + 2*(ndiv+1));
	indices.resize(ixs_start + 6*ndiv);
	unsigned const ixs_off[6] = {1,2,0, 3,2,1}; // 1 quad = 2 triangles

	for (unsigned i = 0; i <= ndiv; ++i) { // vertex data
		unsigned const s(i%ndiv);
		float const ts(1.0f - i*ndiv_inv);
		norm_comp const normal(0.5*(vpn.n[s] + vpn.n[(i+ndiv-1)%ndiv])); // normalize?
		itri_verts[itix++].assign(vpn.p[(s<<1)+0], normal, ts, 0.0, cw.c);
		itri_verts[itix++].assign(vpn.p[(s<<1)+1], normal, ts, 1.0, cw.c);
	}
	for (unsigned i = 0; i < ndiv; ++i) { // index data
		unsigned const ix0(itris_start + 2*i);
		for (unsigned j = 0; j < 6; ++j) {indices[iix++] = ix0 + ixs_off[j];}
	}
	// room object drawing uses back face culling and single sided lighting; to make lighting two sided, need to add verts with inverted normals/winding dirs
	if (two_sided) {add_inverted_triangles(itri_verts, indices, itris_start, ixs_start);}
	// maybe add top and bottom end cap using triangles, currently using all TCs=0.0
	itris_start = itix = itri_verts.size();
	ixs_start   = iix  = indices.size();
	itri_verts.resize(itris_start + (ndiv + 1)*num_ends);
	indices.resize(ixs_start + 3*ndiv*num_ends);

	for (unsigned bt = 0; bt < 2; ++bt) {
		if (!(bt ? draw_top : draw_bot)) continue; // this disk not drawn
		norm_comp const normal((bool(bt) ^ inv_tb) ? plus_z : -plus_z);
		unsigned const center_ix(itix);
		itri_verts[itix++].assign(ce[bt], normal, 0.0, 0.0, cw.c); // center

		for (unsigned i = 0; i < ndiv; ++i) {
			vector3d const &side_normal(vpn.n[i]);
			itri_verts[itix++].assign(vpn.p[(i<<1) + bt], normal, side_normal.x, side_normal.y, cw.c); // assign tcs based on side normal
			indices[iix++] = center_ix; // center
			indices[iix++] = center_ix + i + 1;
			indices[iix++] = center_ix + ((i+1)%ndiv) + 1;
		}
	} // for bt
	if (inv_tb) {std::reverse(indices.begin()+ixs_start, indices.end());} // reverse the order to swap triangle winding order
	if (two_sided) {add_inverted_triangles(itri_verts, indices, itris_start, ixs_start);}
}

void rgeom_mat_t::add_sphere_to_verts(cube_t const &c, colorRGBA const &color) {
	static vector<vert_norm_tc> verts;
	verts.clear();
	static sd_sphere_d sd(all_zeros, 1.0, N_SPHERE_DIV); // resed across all calls
	static sphere_point_norm spn;
	if (!spn.get_points()) {sd.gen_points_norms(spn);} // calculate once and reuse
	sd.get_quad_points(verts); // could use indexed triangles, but this only returns indexed quads
	color_wrapper cw;
	cw.set_c4(color);
	point const center(c.get_cube_center()), size(0.5*c.get_size());
	for (auto i = verts.begin(); i != verts.end(); ++i) {quad_verts.emplace_back(vert_norm_comp_tc((i->v*size + center), i->n, i->t[0], i->t[1]), cw);}
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

void rgeom_mat_t::create_vbo() {
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
	rgeom_alloc.free(*this); // vertex and index data is no longer needed and can be cleared
}

void rgeom_mat_t::draw(shader_t &s, bool shadow_only) {
	if (shadow_only && !en_shadows)  return; // shadows not enabled for this material (picture, whiteboard, rug, etc.)
	if (shadow_only && tex.emissive) return; // assume this is a light source and shouldn't produce shadows
	//if (no_small_features && tex.tid == FONT_TEXTURE_ID) return; // book text is always small
	assert(vbo.vbo_valid());
	assert(num_qverts > 0 || num_itverts > 0);
	if (!shadow_only) {tex.set_gl(s);} // ignores texture scale for now
	vbo.pre_render();
	vertex_t::set_vbo_arrays();
	if (num_qverts > 0) {draw_quads_as_tris(num_qverts);}

	if (num_itverts > 0) { // index quads, used for cylinders
		assert(ivbo > 0);
		bind_vbo(ivbo, 1);
		//glDisable(GL_CULL_FACE); // two sided lighting requires fewer verts (no duplicates), but must be set in the shader
		glDrawRangeElements(GL_TRIANGLES, num_qverts, (num_qverts + num_itverts), num_ixs, GL_UNSIGNED_INT, nullptr);
		//glEnable(GL_CULL_FACE);
		bind_vbo(0, 1);
	}
	if (!shadow_only) {tex.unset_gl(s);}
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
void building_materials_t::create_vbos() {
	for (iterator m = begin(); m != end(); ++m) {m->create_vbo();}
}
void building_materials_t::draw(shader_t &s, bool shadow_only) {
	for (iterator m = begin(); m != end(); ++m) {m->draw(s, shadow_only);}
}

void get_tc_leg_cubes(cube_t const &c, float width, cube_t cubes[4]) {
	float const leg_width(0.5f*width*(c.dx() + c.dy())); // make legs square

	for (unsigned y = 0; y < 2; ++y) {
		for (unsigned x = 0; x < 2; ++x) {
			cube_t leg(c);
			leg.d[0][x] += (x ? -1.0f : 1.0f)*(c.dx() - leg_width);
			leg.d[1][y] += (y ? -1.0f : 1.0f)*(c.dy() - leg_width);
			cubes[2*y+x] = leg;
		}
	}
}
void building_room_geom_t::add_tc_legs(cube_t const &c, colorRGBA const &color, float width, float tscale) {
	rgeom_mat_t &mat(get_wood_material(tscale));
	cube_t cubes[4];
	get_tc_leg_cubes(c, width, cubes);
	for (unsigned i = 0; i < 4; ++i) {mat.add_cube_to_verts(cubes[i], color, c.get_llc(), EF_Z12);} // skip top and bottom faces
}

colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c) {
	if (display_mode & 0x10) return c; // disable this when using indir lighting
	return c * (0.5f + 0.5f*min(sqrt(o.light_amt), 1.5f)); // use c.light_amt as an approximation for ambient lighting due to sun/moon
}
colorRGBA apply_light_color(room_object_t const &o) {return apply_light_color(o, o.color);} // use object color

tid_nm_pair_t const untex_shad_mat(-1, 2.0); // make sure it's different from default tid_nm_pair_t so that it's not grouped with shadowed materials

void building_room_geom_t::add_table(room_object_t const &c, float tscale) { // 6 quads for top + 4 quads per leg = 22 quads = 88 verts
	cube_t top(c), legs_bcube(c);
	top.z1() += 0.85*c.dz(); // 15% of height
	legs_bcube.z2() = top.z1();
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(top, color, c.get_llc()); // all faces drawn
	add_tc_legs(legs_bcube, color, 0.08, tscale);
}

void building_room_geom_t::add_chair(room_object_t const &c, float tscale) { // 6 quads for seat + 5 quads for back + 4 quads per leg = 27 quads = 108 verts
	float const height(c.dz()*((c.type == TYPE_SM_CHAIR) ? 1.333 : 1.0)); // effective height if the chair wasn't short
	cube_t seat(c), back(c), legs_bcube(c);
	seat.z1() += 0.32*height;
	seat.z2()  = back.z1() = seat.z1() + 0.07*height;
	legs_bcube.z2() = seat.z1();
	back.d[c.dim][c.dir] += 0.88f*(c.dir ? -1.0f : 1.0f)*c.get_sz_dim(c.dim);
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.2*tscale), 1).add_cube_to_verts(seat, apply_light_color(c), c.get_llc()); // all faces drawn
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(back, color, c.get_llc(), EF_Z1); // skip bottom face
	add_tc_legs(legs_bcube, color, 0.15, tscale);
}

void building_room_geom_t::add_stair(room_object_t const &c, float tscale, vector3d const &tex_origin) {
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.5*tscale), 1).add_cube_to_verts(c, colorRGBA(0.85, 0.85, 0.85), tex_origin); // all faces drawn
}

tid_nm_pair_t get_tex_auto_nm(int tid, float tscale) {
	return tid_nm_pair_t(tid, get_normal_map_for_bldg_tid(tid), tscale, tscale);
}

void building_room_geom_t::add_elevator(room_object_t const &c, float tscale) {
	// elevator car, all materials are dynamic
	float const thickness(0.051*c.dz());
	cube_t floor(c), ceil(c), back(c);
	floor.z2() = floor.z1() + thickness;
	ceil. z1() = ceil. z2() - thickness;
	floor.expand_by_xy(-0.5f*thickness);
	ceil .expand_by_xy(-0.5f*thickness);
	back.d[c.dim][c.dir] = back.d[c.dim][!c.dir] + (c.dir ? 1.0 : -1.0)*thickness;
	vector3d const tex_origin(c.get_llc());
	unsigned const front_face_mask(get_face_mask(c.dim, c.dir)), floor_ceil_face_mask(front_face_mask & 60); // +Z faces
	tid_nm_pair_t const paneling(get_tex_auto_nm(PANELING_TEX, 2.0f*tscale));
	get_material(get_tex_auto_nm(TILE_TEX, tscale), 1, 1).add_cube_to_verts(floor, WHITE, tex_origin, floor_ceil_face_mask);
	get_material(get_tex_auto_nm(get_rect_panel_tid(), tscale), 1, 1).add_cube_to_verts(ceil, WHITE, tex_origin, floor_ceil_face_mask);
	get_material(paneling, 1, 1).add_cube_to_verts(back, WHITE, tex_origin, front_face_mask, !c.dim);

	for (unsigned d = 0; d < 2; ++d) { // side walls
		cube_t side(c);
		side.d[!c.dim][!d] = side.d[!c.dim][d] + (d ? -1.0 : 1.0)*thickness;
		get_material(paneling, 1, 1).add_cube_to_verts(side, WHITE, tex_origin, get_face_mask(!c.dim, !d), c.dim);
	}
}

void building_room_geom_t::add_light(room_object_t const &c, float tscale) {
	// Note: need to use a different texture (or -1) for is_on because emissive flag alone does not cause a material change
	bool const is_on(c.is_lit());
	tid_nm_pair_t tp(((is_on || c.shape == SHAPE_SPHERE) ? (int)WHITE_TEX : (int)PLASTER_TEX), tscale);
	tp.emissive = is_on;
	rgeom_mat_t &mat(mats_lights.get_material(tp, 0)); // no shadows
	if      (c.shape == SHAPE_CUBE  ) {mat.add_cube_to_verts  (c, c.color, c.get_llc(), EF_Z2);} // untextured, skip top face
	else if (c.shape == SHAPE_CYLIN ) {mat.add_vcylin_to_verts(c, c.color, 1, 0);} // bottom only
	else if (c.shape == SHAPE_SPHERE) {mat.add_sphere_to_verts(c, c.color);}
	else {assert(0);}
}

void building_room_geom_t::add_rug(room_object_t const &c) {
	bool const swap_tex_st(c.dy() < c.dx()); // rug textures are oriented with the long side in X, so swap the coordinates (rotate 90 degrees) if our rug is oriented the other way
	get_material(tid_nm_pair_t(c.get_rug_tid(), 0.0)).add_cube_to_verts(c, WHITE, c.get_llc(), 61, swap_tex_st); // only draw top/+z face
}

void building_room_geom_t::add_picture(room_object_t const &c) { // also whiteboards
	bool const whiteboard(c.type == TYPE_WBOARD);
	int picture_tid(WHITE_TEX);

	if (!whiteboard) { // picture
		int const user_tid(get_rand_screenshot_texture(c.obj_id));
		picture_tid  = ((user_tid >= 0) ? (unsigned)user_tid : c.get_picture_tid()); // if user texture is valid, use that instead
		num_pic_tids = get_num_screenshot_tids();
		has_pictures = 1;
	}
	unsigned skip_faces(get_face_mask(c.dim, c.dir)); // only the face oriented outward
	bool const mirror_x(!whiteboard && !(c.dim ^ c.dir));
	vector3d const tex_origin(c.get_llc());
	get_material(tid_nm_pair_t(picture_tid, 0.0)).add_cube_to_verts(c, WHITE, tex_origin, skip_faces, !c.dim, mirror_x);
	// add a frame
	cube_t frame(c);
	vector3d exp;
	exp.z = exp[!c.dim] = (whiteboard ? 0.04 : 0.06)*c.dz(); // frame width
	exp[c.dim] = (whiteboard ? -0.1 : -0.25)*c.get_sz_dim(c.dim); // shrink in this dim
	frame.expand_by(exp);
	get_material(tid_nm_pair_t()).add_cube_to_verts(frame, (whiteboard ? GRAY : BLACK), tex_origin, skip_faces, 0);
	
	if (whiteboard) { // add a marker ledge
		cube_t ledge(c);
		ledge.z2() = ledge.z1() + 0.016*c.dz(); // along the bottom edge
		ledge.d[c.dim][c.dir] += (c.dir ? 1.5 : -1.5)*c.get_sz_dim(c.dim); // extrude outward
		get_material(untex_shad_mat, 1).add_cube_to_verts(ledge, GRAY, tex_origin, (1 << (2*(2-c.dim) + !c.dir)), 0); // shadowed
	}
}

void building_room_geom_t::add_book_title(string const &title, cube_t const &title_area, rgeom_mat_t &mat, colorRGBA const &color,
	unsigned hdim, unsigned tdim, unsigned wdim, bool cdir, bool ldir, bool wdir)
{
	vector3d column_dir(zero_vector), line_dir(zero_vector), normal(zero_vector);
	column_dir[hdim] = (cdir ? -1.0 : 1.0); // along book height
	line_dir  [tdim] = (ldir ? -1.0 : 1.0); // along book thickness
	normal    [wdim] = (wdir ? -1.0 : 1.0); // along book width
	static vector<vert_tc_t> verts;
	verts.clear();
	gen_text_verts(verts, all_zeros, title, 1.0, column_dir, line_dir, 1); // use_quads=1 (could cache this for c.obj_id + dim/dir bits)
	assert(!verts.empty());
	cube_t text_bcube(verts[0].v);
	for (auto i = verts.begin()+2; i != verts.end(); i += 2) {text_bcube.union_with_pt(i->v);} // only need to include opposite corners
	float const wscale(title_area.get_sz_dim(hdim)/text_bcube.get_sz_dim(hdim)), hscale(title_area.get_sz_dim(tdim)/text_bcube.get_sz_dim(tdim));
	float width_scale(wscale), height_scale(hscale);
	min_eq(width_scale,  1.5f*height_scale); // use a reasonable aspect ratio
	min_eq(height_scale, 1.5f*width_scale );
	float const title_start_hdim(title_area.d[hdim][cdir] + column_dir[hdim]*0.5*title_area.get_sz_dim(hdim)*(1.0 -  width_scale/wscale)); // centered
	float const title_start_tdim(title_area.d[tdim][ldir] + line_dir  [tdim]*0.5*title_area.get_sz_dim(tdim)*(1.0 - height_scale/hscale)); // centered
	if (dot_product(normal, cross_product((verts[1].v - verts[0].v), (verts[2].v - verts[1].v))) < 0.0) {std::reverse(verts.begin(), verts.end());} // swap vertex winding order
	color_wrapper const cw(color);
	norm_comp const nc(normal);

	for (auto i = verts.begin(); i != verts.end(); ++i) {
		i->v[wdim] = title_area.d[wdim][!wdir]; // spine pos
		i->v[hdim] = (i->v[hdim] - text_bcube.d[hdim][cdir])*width_scale  + title_start_hdim;
		i->v[tdim] = (i->v[tdim] - text_bcube.d[tdim][ldir])*height_scale + title_start_tdim;
		mat.quad_verts.emplace_back(vert_norm_comp_tc(i->v, nc, i->t[0], i->t[1]), cw);
	} // for i
}

void building_room_geom_t::add_book(room_object_t const &c, bool inc_lg, bool inc_sm, float tilt_angle, unsigned extra_skip_faces, bool no_title) {
	bool const upright(c.get_sz_dim(!c.dim) < c.dz());
	bool const tdir(upright ? (c.dim ^ c.dir ^ bool(c.obj_id%7)) : 1); // sometimes upside down when upright
	bool const ldir(!tdir), cdir(c.dim ^ c.dir ^ upright ^ ldir); // colum and line directions (left/right/top/bot) + mirror flags for front cover
	unsigned const tdim(upright ? !c.dim : 2), hdim(upright ? 2 : !c.dim); // thickness dim, height dim (c.dim is width dim)
	float const thickness(c.get_sz_dim(tdim)), width(c.get_sz_dim(c.dim)), cov_thickness(0.125*thickness), indent(0.02*width);
	cube_t bot(c), top(c), spine(c), pages(c), cover(c);
	bot.d[tdim][1] = c.d[tdim][0] + cov_thickness;
	top.d[tdim][0] = c.d[tdim][1] - cov_thickness;
	pages.d[tdim][0] = spine.d[tdim][0] = bot.d[tdim][1];
	pages.d[tdim][1] = spine.d[tdim][1] = top.d[tdim][0];
	vector3d shrink(zero_vector);
	shrink[c.dim] = shrink[upright ? 2 : !c.dim] = -indent;
	pages.expand_by(shrink);
	spine.d[c.dim][c.dir] = pages.d[c.dim][!c.dir];
	vector3d const tex_origin(c.get_llc());
	vector3d axis, about(c.get_urc());
	axis[c.dim] = 1.0; // along book width
	tilt_angle *= (c.dim ? -1.0 : 1.0);
	bool has_cover(0);

	if (inc_lg) { // add book geom
		colorRGBA const color(apply_light_color(c));
		// skip top face, bottom face if not tilted, thickness dim if upright
		unsigned const skip_faces(extra_skip_faces | ((tilt_angle == 0.0) ? EF_Z1 : 0) | (upright ? get_skip_mask_for_xy(tdim) : EF_Z2));
		rgeom_mat_t &mat(get_material(tid_nm_pair_t(), 0)); // unshadowed, since shadows are too small to have much effect
		unsigned const qv_start(mat.quad_verts.size());
		mat.add_cube_to_verts(bot,   color, tex_origin, (extra_skip_faces | EF_Z1)); // untextured, skip bottom face
		mat.add_cube_to_verts(top,   color, tex_origin, (extra_skip_faces | (upright ? EF_Z1 : 0))); // untextured, skip bottom face if upright
		mat.add_cube_to_verts(spine, color, tex_origin, skip_faces); // untextured
		mat.add_cube_to_verts(pages, apply_light_color(c, WHITE), tex_origin, (skip_faces | ~get_face_mask(c.dim, !c.dir))); // untextured
		rotate_verts(mat.quad_verts, axis, tilt_angle, about, qv_start);
	}
	if (ADD_BOOK_COVERS && inc_sm && c.enable_pictures() && (upright || (c.obj_id&2))) { // add picture to book cover
		vector3d expand;
		float const height(c.get_sz_dim(hdim)), img_width(0.9*width), img_height(min(0.9f*height, 0.67f*img_width)); // use correct aspect ratio
		expand[ hdim] = -0.5f*(height - img_height);
		expand[c.dim] = -0.5f*(width  - img_width);
		expand[ tdim] = 0.1*indent; // expand outward, other dims expand inward
		cover.expand_by(expand);
		int const picture_tid(c.get_picture_tid()); // not using user screenshot images
		bool const swap_xy(upright ^ (!c.dim));
		rgeom_mat_t &cover_mat(get_material(tid_nm_pair_t(picture_tid, 0.0), 0, 0, 1));
		unsigned const qv_start(cover_mat.quad_verts.size());
		cover_mat.add_cube_to_verts(cover, WHITE, tex_origin, get_face_mask(tdim, tdir), swap_xy, ldir, !cdir); // no shadows, small=1
		rotate_verts(cover_mat.quad_verts, axis, tilt_angle, about, qv_start);
		has_cover = 1;
	} // end cover image
	bool const add_spine_title(c.obj_id & 7); // 7/8 of the time

	if (ADD_BOOK_TITLES && inc_sm && !no_title && (!upright || add_spine_title)) {
		unsigned const SPLIT_LINE_SZ = 24;
		string const &title(gen_book_title(c.obj_id, nullptr, SPLIT_LINE_SZ)); // select our title text
		if (title.empty()) return; // no title
		colorRGBA text_color(BLACK);
		for (unsigned i = 0; i < 3; ++i) {text_color[i] = ((c.color[i] > 0.5) ? 0.0 : 1.0);} // invert + saturate to contrast with book cover
		text_color = apply_light_color(c, text_color);
		rgeom_mat_t &mat(get_material(tid_nm_pair_t(FONT_TEXTURE_ID), 0, 0, 1)); // no shadows, small=1
		unsigned const qv_start(mat.quad_verts.size());

		if (add_spine_title) { // add title along spine
			cube_t title_area(c);
			vector3d expand;
			expand[ hdim] = -4.0*indent; // shrink
			expand[ tdim] = -1.0*indent; // shrink
			expand[c.dim] =  0.2*indent; // expand outward
			title_area.expand_by(expand);
			add_book_title(title, title_area, mat, text_color, hdim, tdim, c.dim, cdir, ldir, c.dir);
		}
		if (!upright && (!add_spine_title || (c.obj_id%3))) { // add title to front cover if upright
			cube_t title_area_fc(c);
			title_area_fc.z1()  = title_area_fc.z2();
			title_area_fc.z2() += 0.2*indent;
			title_area_fc.expand_in_dim(c.dim, -4.0*indent);
			bool const top_dir(c.dim ^ c.dir);
			
			if (has_cover) { // place above cover; else, place in center
				title_area_fc.d[!c.dim][!top_dir] = cover.d[!c.dim][top_dir];
				title_area_fc.expand_in_dim(!c.dim, -1.0*indent);
			}
			add_book_title(title, title_area_fc, mat, text_color, c.dim, !c.dim, 2, !c.dir, !top_dir, 0); // {columns, lines, normal}
		}
		rotate_verts(mat.quad_verts, axis, tilt_angle, about, qv_start);
	} // end pages
}

void building_room_geom_t::add_bookcase(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale, bool no_shelves, float sides_scale) {
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	unsigned const skip_faces(~get_face_mask(c.dim, !c.dir)); // skip back face
	unsigned const skip_faces_shelves(skip_faces | get_skip_mask_for_xy(!c.dim)); // skip back face and sides
	float const width(c.get_sz_dim(!c.dim)), depth((c.dir ? -1.0 : 1.0)*c.get_sz_dim(c.dim)); // signed depth
	float const side_thickness(0.06*sides_scale*width);
	vector3d const tex_origin(c.get_llc());
	cube_t middle(c);

	for (unsigned d = 0; d < 2; ++d) { // left/right sides
		cube_t lr(c);
		lr.d[!c.dim][d] += (d ? -1.0f : 1.0f)*(width - side_thickness);
		if (inc_lg) {get_wood_material(tscale).add_cube_to_verts(lr, color, tex_origin, (skip_faces | EF_Z1));} // side
		middle.d[!c.dim][!d] = lr.d[!c.dim][d];
	}
	cube_t top(middle);
	top.z1()   += c.dz() - side_thickness; // make same width as sides
	middle.z2() = top.z1();
	if (inc_lg) {get_wood_material(tscale).add_cube_to_verts(top, color, tex_origin, skip_faces_shelves);} // top
	cube_t back(middle);
	back.d[c.dim] [c.dir]  += 0.94*depth;
	middle.d[c.dim][!c.dir] = back.d[c.dim][c.dir];
	if (inc_lg) {get_wood_material(tscale).add_cube_to_verts(back, color, tex_origin, get_face_mask(c.dim, c.dir));} // back - only face oriented outward
	if (no_shelves) return;
	// add shelves
	rand_gen_t rgen;
	rgen.set_state(c.obj_id+1, c.room_id+1);
	unsigned const num_shelves(3 + ((rgen.rand() + c.obj_id)%3)); // 3-5
	float const shelf_dz(middle.dz()/(num_shelves+0.25)), shelf_thick(0.03*c.dz());
	cube_t shelves[5];
	
	for (unsigned i = 0; i < num_shelves; ++i) {
		cube_t &shelf(shelves[i]);
		shelf = middle; // copy XY parts
		shelf.z1() += (i+0.25)*shelf_dz;
		shelf.z2()  = shelf.z1() + shelf_thick;
		if (inc_lg) {get_wood_material(tscale).add_cube_to_verts(shelf, color, tex_origin, skip_faces_shelves);} // Note: mat reference may be invalidated by adding books
	}
	// add books
	for (unsigned i = 0; i < num_shelves; ++i) {
		if (rgen.rand_float() < 0.2) continue; // no books on this shelf
		cube_t const &shelf(shelves[i]);
		unsigned const num_spaces(22 + (rgen.rand()%11)); // 22-32 books per shelf
		float const book_space(shelf.get_sz_dim(!c.dim)/num_spaces);
		float pos(shelf.d[!c.dim][0]), shelf_end(shelf.d[!c.dim][1]), last_book_pos(pos), min_height(0.0);
		unsigned skip_mask(0);
		bool prev_tilted(0);

		for (unsigned n = 0; n < num_spaces; ++n) {
			if (rgen.rand_float() < 0.12) {
				unsigned const skip_end(n + (rgen.rand()%8) + 1); // skip 1-8 books
				for (; n < skip_end; ++n) {skip_mask |= (1<<n);}
			}
		}
		for (unsigned n = 0; n < num_spaces; ++n) {
			if ((pos + 0.7*book_space) > shelf_end) break; // not enough space for another book
			float const width(book_space*rgen.rand_uniform(0.7, 1.3));
			if (!prev_tilted && (skip_mask & (1<<n))) {pos += width; continue;} // skip this book, and don't tilt the next one
			float const height(max((shelf_dz - shelf_thick)*rgen.rand_uniform(0.6, 0.98), min_height));
			float const right_pos(min((pos + width), shelf_end)), avail_space(right_pos - last_book_pos);
			float tilt_angle(0.0);
			cube_t book;
			book.z1() = shelf.z2();
			book.d[c.dim][ c.dir] = shelf.d[c.dim][ c.dir] + depth*rgen.rand_uniform(0.0, 0.25); // facing out
			book.d[c.dim][!c.dir] = shelf.d[c.dim][!c.dir]; // facing in
			min_height = 0.0;

			if (avail_space > 1.1f*height && rgen.rand_float() < 0.5) { // book has space to fall over 50% of the time
				book.d[!c.dim][0] = last_book_pos + rgen.rand_uniform(0.0, (right_pos - last_book_pos - height)); // shift a random amount within the gap
				book.d[!c.dim][1] = book.d[!c.dim][0] + height;
				book.z2() = shelf.z2() + width;
			}
			else { // upright
				if (!prev_tilted && avail_space > 2.0*width && (right_pos + book_space) < shelf_end && n+1 < num_spaces) { // rotates about the URC
					float const lean_width(min((avail_space - width), rgen.rand_uniform(0.1, 0.6)*height)); // use part of the availabe space to lean
					tilt_angle = asinf(lean_width/height);
					float const delta_z(height - sqrt(height*height - lean_width*lean_width)); // move down to touch the bottom of the bookshelf when rotated
					book.z1() -= delta_z;
					min_height = rgen.rand_uniform(0.95, 1.05)*(height - delta_z); // make sure the book this book is leaning on is tall enough
				}
				book.d[!c.dim][0] = pos;
				book.d[!c.dim][1] = right_pos; // clamp to edge of bookcase interior
				book.z2() = book.z1() + height;
				assert(pos < right_pos);
			}
			assert(book.is_strictly_normalized());
			colorRGBA const &book_color(book_colors[rgen.rand() % NUM_BOOK_COLORS]);
			bool const backwards((rgen.rand()&3) == 0), book_dir(c.dir ^ backwards ^ 1); // spine facing out 75% of the time
			room_object_t obj(book, TYPE_BOOK, c.room_id, c.dim, book_dir, c.flags, c.light_amt, room_obj_shape::SHAPE_CUBE, book_color);
			obj.obj_id = c.obj_id + 123*i + 1367*n;
			add_book(obj, inc_lg, inc_sm, tilt_angle, skip_faces, backwards); // detailed book, no title if backwards
			pos += width;
			last_book_pos = pos;
			prev_tilted   = (tilt_angle != 0.0); // don't tilt two books in a row
		} // for n
	} // for i
}

void building_room_geom_t::add_desk(room_object_t const &c, float tscale) {
	// desk top and legs, similar to add_table()
	cube_t top(c), legs_bcube(c);
	top.z1() += 0.8*c.dz();
	legs_bcube.z2() = top.z1();
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(top, color, c.get_llc()); // all faces drawn
	add_tc_legs(legs_bcube, color, 0.06, tscale);

	if (c.shape == SHAPE_TALL) { // add top/back section of desk; this part is outside the bcube
		room_object_t c_top_back(c);
		c_top_back.z1() = top.z2();
		c_top_back.z2() = top.z2() + 1.8*c.dz();
		c_top_back.d[c.dim][c.dir] += 0.75*(c.dir ? -1.0 : 1.0)*c.get_sz_dim(c.dim);
		add_bookcase(c_top_back, 1, 1, tscale, 1, 0.4); // no_shelves=1, side_width=0.4, both large and small
	}
}

void add_pillow(cube_t const &c, rgeom_mat_t &mat, colorRGBA const &color, vector3d const &tex_origin) {
	unsigned const ndiv = 24; // number of quads in X and Y
	float const ndiv_inv(1.0/ndiv), dx_inv(1.0/c.dx()), dy_inv(1.0/c.dy());
	color_wrapper cw(color);
	auto &verts(mat.itri_verts); // Note: could cache verts
	unsigned const start(verts.size()), stride(ndiv + 1);
	float dists[ndiv+1];
	norm_comp const nc(plus_z);

	for (unsigned x = 0; x <= ndiv; ++x) {
		float const v(2.0f*x*ndiv_inv - 1.0f); // centered on 0 in range [-1, 1]
		dists[x] = 0.5*SIGN(v)*sqrt(abs(v)) + 0.5; // nonlinear spacing, closer near the edges, convert back to [0, 1] range
	}
	for (unsigned y = 0; y <= ndiv; ++y) {
		float const yval(c.y1() + dists[y]*c.dy()), ey(2.0f*max(0.0f, min((yval - c.y1()), (c.y2() - yval)))*dy_inv);

		for (unsigned x = 0; x <= ndiv; ++x) {
			float const xval(c.x1() + dists[x]*c.dx()), ex(2.0f*max(0.0f, min((xval - c.x1()), (c.x2() - xval)))*dx_inv), zval(c.z1() + c.dz()*pow(ex*ey, 0.2f));
			verts.emplace_back(vert_norm_comp_tc(point(xval, yval, zval), nc, mat.tex.tscale_x*(xval - tex_origin.x), mat.tex.tscale_y*(yval - tex_origin.y)), cw);
		} // for x
	} // for y
	for (unsigned y = 0; y <= ndiv; ++y) {
		for (unsigned x = 0; x <= ndiv; ++x) {
			unsigned const off(start + y*stride + x);
			vector3d const &v(verts[off].v);
			vector3d normal(zero_vector);
			if (x > 0    && y >    0) {normal += cross_product((v - verts[off-stride].v), (verts[off-1].v - v));} // LL
			if (x < ndiv && y >    0) {normal += cross_product((v - verts[off+1].v), (verts[off-stride].v - v));} // LR
			if (x < ndiv && y < ndiv) {normal += cross_product((v - verts[off+stride].v), (verts[off+1].v - v));} // UR
			if (x > 0    && y < ndiv) {normal += cross_product((v - verts[off-1].v), (verts[off+stride].v - v));} // UL
			verts[off].set_norm(normal.get_norm()); // this is the slowest line
		} // for x
	} // for y
	for (unsigned y = 0; y < ndiv; ++y) {
		for (unsigned x = 0; x < ndiv; ++x) {
			unsigned const off(start + y*stride + x);
			mat.indices.push_back(off + 0); // T1
			mat.indices.push_back(off + 1);
			mat.indices.push_back(off + stride+1);
			mat.indices.push_back(off + 0); // T2
			mat.indices.push_back(off + stride+1);
			mat.indices.push_back(off + stride);
		} // for x
	} // for y
}

void building_room_geom_t::add_bed(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale) {
	float const height(c.dz()), length(c.get_sz_dim(c.dim)), width(c.get_sz_dim(!c.dim));
	bool const is_wide(width > 0.7*length);
	cube_t frame(c), head(c), foot(c), mattress(c), legs_bcube(c), pillow(c);
	head.d[c.dim][!c.dir] += (c.dir ? 1.0 : -1.0)*0.96*length;
	foot.d[c.dim][ c.dir] -= (c.dir ? 1.0 : -1.0)*0.97*length;
	mattress.d[c.dim][ c.dir] = head.d[c.dim][!c.dir];
	mattress.d[c.dim][!c.dir] = foot.d[c.dim][ c.dir];
	frame.z1() += 0.3*height;
	frame.z2() -= 0.65*height;
	foot.z2()  -= 0.15*height;
	mattress.z1()   = head.z1()   = foot.z1() = frame.z2();
	mattress.z2()   = pillow.z1() = mattress.z1() + 0.2*height;
	pillow.z2()     = pillow.z1() + 0.13*height;
	legs_bcube.z2() = frame.z1();
	float const pillow_space((is_wide ? 0.08 : 0.23)*width);
	pillow.expand_in_dim(!c.dim, -pillow_space);
	pillow.d[c.dim][ c.dir] = mattress.d[c.dim][ c.dir] + (c.dir ? -1.0 : 1.0)*0.02*length; // head
	pillow.d[c.dim][!c.dir] = pillow  .d[c.dim][ c.dir] + (c.dir ? -1.0 : 1.0)*(is_wide ? 0.25 : 0.6)*pillow.get_sz_dim(!c.dim);
	mattress.expand_in_dim(!c.dim, -0.02*width);
	colorRGBA const sheet_color(apply_light_color(c));
	tid_nm_pair_t const sheet_tex(c.get_sheet_tid(), tscale);

	if (inc_lg) {
		colorRGBA const color(apply_light_color(c, WOOD_COLOR));
		add_tc_legs(legs_bcube, color, 0.04, tscale);
		rgeom_mat_t &wood_mat(get_wood_material(tscale));
		vector3d const tex_origin(c.get_llc());
		wood_mat.add_cube_to_verts(frame, color, tex_origin);
		wood_mat.add_cube_to_verts(head, color, tex_origin, EF_Z1);
		wood_mat.add_cube_to_verts(foot, color, tex_origin, EF_Z1);
		unsigned const mattress_skip_faces(EF_Z1 | get_skip_mask_for_xy(c.dim));
		rgeom_mat_t &sheet_mat(get_material(sheet_tex, 1));
		sheet_mat.add_cube_to_verts(mattress, sheet_color, tex_origin, mattress_skip_faces);
	}
	if (inc_sm) {
		rgeom_mat_t &pillow_mat(get_material(sheet_tex, 1, 0, 1)); // small=1

		if (is_wide) { // two pillows
			for (unsigned d = 0; d < 2; ++d) {
				cube_t p(pillow);
				p.d[!c.dim][d] += (d ? -1.0 : 1.0)*0.55*pillow.get_sz_dim(!c.dim);
				add_pillow(p, pillow_mat, sheet_color, tex_origin);
			}
		}
		else {add_pillow(pillow, pillow_mat, sheet_color, tex_origin);} // one pillow
	}
}

void building_room_geom_t::add_trashcan(room_object_t const &c) {
	rgeom_mat_t &mat(get_material(untex_shad_mat, 1));
	colorRGBA const color(apply_light_color(c));

	if (c.shape == room_obj_shape::SHAPE_CYLIN) {
		mat.add_vcylin_to_verts(c, color, 1, 0, 1, 0, 1, 0.7, 1.0); // untextured, bottom only, two_sided cylinder with inverted bottom normal
	}
	else { // sloped cube; this shape is rather unique, so is drawn inline; untextured
		cube_t base(c);
		base.expand_by_xy(vector3d(-0.2*c.dx(), -0.2*c.dy(), 0.0)); // shrink base by 40%
		auto &verts(mat.quad_verts);
		rgeom_mat_t::vertex_t v;
		v.set_c4(color);
		v.set_ortho_norm(2, 1); // +z
		
		for (unsigned i = 0; i < 4; ++i) { // bottom
			bool const xp(i==0||i==1), yp(i==1||i==2);
			v.v.assign(base.d[0][xp], base.d[1][yp], base.z1());
			v.t[0] = float(xp); v.t[1] = float(yp); // required for normal mapping ddx/ddy on texture coordinate
			verts.push_back(v);
		}
		for (unsigned dim = 0; dim < 2; ++dim) { // x,y
			for (unsigned dir = 0; dir < 2; ++dir) {
				unsigned const six(verts.size());

				for (unsigned i = 0; i < 4; ++i) {
					bool const tb(i==1||i==2), lohi(i==0||i==1);
					v.v[ dim] = (tb ? (cube_t)c : base).d[ dim][dir];
					v.v[!dim] = (tb ? (cube_t)c : base).d[!dim][lohi];
					v.v.z  = c.d[2][tb];
					//v.t[0] = float(tb); v.t[1] = float(lohi); // causes a seam between triangles due to TBN basis change, so leave at 0.0
					verts.push_back(v);
				}
				for (unsigned i = 0; i < 4; ++i) {verts.push_back(verts[six+3-i]);} // add reversed quad for opposing face
				norm_comp n(cross_product((verts[six].v - verts[six+1].v), (verts[six].v - verts[six+2].v)).get_norm());
				for (unsigned i = 0; i < 4; ++i) {verts[six+i].set_norm(n);} // front face
				n.invert_normal();
				for (unsigned i = 4; i < 8; ++i) {verts[six+i].set_norm(n);} // back face
			} // for dir
		} // for dim
	}
}

void building_room_geom_t::add_br_stall(room_object_t const &c) {
	rgeom_mat_t &mat(get_material(untex_shad_mat, 1));
	colorRGBA const color(apply_light_color(c));
	point const tex_origin(c.get_llc()); // doesn't really need to be set, since stall is untextured
	float const dz(c.dz()), wall_thick(0.0125*dz), frame_thick(2.0*wall_thick), door_gap(0.2*wall_thick);
	cube_t sides(c), front(c);
	sides.z2() -= 0.35*dz;
	sides.z1() += 0.15*dz;
	sides.d[c.dim][!c.dir] += (c.dir ? 1.0 : -1.0)*wall_thick; // shorten for door
	front.d[c.dim][ c.dir] = sides.d[c.dim][!c.dir];
	cube_t side1(sides), side2(sides), front1(front), front2(front), door(front);
	door.z2() -= 0.38*dz;
	door.z1() += 0.18*dz;
	side1.d[!c.dim][1] = side1.d[!c.dim][0] + wall_thick;
	side2.d[!c.dim][0] = side2.d[!c.dim][1] - wall_thick;
	door.expand_in_dim(!c.dim, -frame_thick);
	front1.d[!c.dim][1] = door.d[!c.dim][0];
	front2.d[!c.dim][0] = door.d[!c.dim][1];
	door.expand_in_dim(!c.dim, -door_gap);
	unsigned const side_skip_mask(get_skip_mask_for_xy(c.dim));
	mat.add_cube_to_verts(side1,  color, tex_origin, side_skip_mask);
	mat.add_cube_to_verts(side2,  color, tex_origin, side_skip_mask);
	mat.add_cube_to_verts(front1, color, tex_origin, EF_Z12);
	mat.add_cube_to_verts(front2, color, tex_origin, EF_Z12);
	mat.add_cube_to_verts(door,   color, tex_origin);
}

void building_room_geom_t::add_cubicle(room_object_t const &c, float tscale) {
	int const tid(get_texture_by_name((c.obj_id & 1) ? "carpet/carpet1.jpg" : "carpet/carpet2.jpg")); // select from one of 2 textures
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(tid, tscale), 1));
	colorRGBA const color(apply_light_color(c));
	point const tex_origin(c.get_llc());
	float const wall_thick(0.07*c.dz()), frame_thick(5.0*wall_thick);
	bool const is_short(c.shape == SHAPE_SHORT);
	cube_t sides(c), front(c), back(c);
	if (is_short) {back.z2() -= 0.4*c.dz();}
	sides.d[c.dim][!c.dir] += (c.dir ? 1.0 : -1.0)*wall_thick; // front
	sides.d[c.dim][ c.dir] -= (c.dir ? 1.0 : -1.0)*wall_thick; // back
	front.d[c.dim][ c.dir] = sides.d[c.dim][!c.dir];
	back .d[c.dim][!c.dir] = sides.d[c.dim][ c.dir];
	cube_t side1(sides), side2(sides), front1(front), front2(front);
	side1 .d[!c.dim][1] = side1.d[!c.dim][0] + wall_thick;
	side2 .d[!c.dim][0] = side2.d[!c.dim][1] - wall_thick;
	front1.d[!c.dim][1] = front.d[!c.dim][0] + frame_thick;
	front2.d[!c.dim][0] = front.d[!c.dim][1] - frame_thick;
	unsigned const side_skip_mask (EF_Z12 | get_skip_mask_for_xy( c.dim));
	unsigned const front_skip_mask(EF_Z12 | get_skip_mask_for_xy(!c.dim));
	mat.add_cube_to_verts(side1,  color, tex_origin, side_skip_mask);
	mat.add_cube_to_verts(side2,  color, tex_origin, side_skip_mask);
	mat.add_cube_to_verts(front1, color, tex_origin, front_skip_mask);
	mat.add_cube_to_verts(front2, color, tex_origin, front_skip_mask);
	mat.add_cube_to_verts(back,   color, tex_origin, EF_Z12);
	// black edges on walls
	rgeom_mat_t &edge_mat(get_material(untex_shad_mat, 1)); // shadowed?
	unsigned const side_edge_skip_mask (~(EF_Z2 | (is_short ? ~get_face_mask(c.dim, c.dir) : 0)));
	unsigned const front_edge_skip_mask(~(EF_Z2 | get_skip_mask_for_xy(!c.dim)));
	edge_mat.add_cube_to_verts(side1,  BKGRAY, tex_origin, side_edge_skip_mask);
	edge_mat.add_cube_to_verts(side2,  BKGRAY, tex_origin, side_edge_skip_mask);
	edge_mat.add_cube_to_verts(front1, BKGRAY, tex_origin, front_edge_skip_mask);
	edge_mat.add_cube_to_verts(front2, BKGRAY, tex_origin, front_edge_skip_mask);
	edge_mat.add_cube_to_verts(back,   BKGRAY, tex_origin, ~EF_Z2);
}

class sign_helper_t {
	map<string, unsigned> txt_to_id;
	vector<string> text;
public:
	unsigned register_text(string const &t) {
		auto it(txt_to_id.find(t));
		if (it != txt_to_id.end()) return it->second; // found
		unsigned const id(text.size());
		txt_to_id[t] = id; // new text, insert it
		text.push_back(t);
		assert(text.size() == txt_to_id.size());
		return id;
	}
	string const &get_text(unsigned id) const {
		assert(id < text.size());
		return text[id];
	}
};

sign_helper_t sign_helper;

unsigned register_sign_text(string const &text) {return sign_helper.register_text(text);}

void building_room_geom_t::add_sign(room_object_t const &c, bool inc_back, bool inc_text) {
	if (inc_back) {
		unsigned const skip_faces(~get_face_mask(c.dim, !c.dir)); // skip back face
		get_material(tid_nm_pair_t(), 0).add_cube_to_verts(c, WHITE, zero_vector, skip_faces); // back of the sign, always white (for now)
	}
	if (!inc_text) return;
	// add sign text
	cube_t ct(c); // text area is slightly smaller than full cube
	ct.expand_in_dim(!c.dim, -0.1*c.get_sz_dim(!c.dim));
	ct.expand_in_dim(2, 0.1*c.dz());
	vector3d col_dir(zero_vector), normal(zero_vector);
	bool const ldir(c.dim ^ c.dir);
	col_dir[!c.dim] = (ldir  ? 1.0 : -1.0);
	normal [ c.dim] = (c.dir ? 1.0 : -1.0);
	static vector<vert_tc_t> verts;
	verts.clear();
	string const &text(sign_helper.get_text(c.obj_id));
	assert(!text.empty());
	point pos;
	pos[c.dim] = ct.d[c.dim][c.dir] + (c.dir ? 1.0 : -1.0)*0.1*ct.get_sz_dim(c.dim); // normal
	gen_text_verts(verts, pos, text, 1.0, col_dir, plus_z, 1); // use_quads=1
	assert(!verts.empty());
	cube_t text_bcube(verts[0].v);
	for (auto i = verts.begin()+2; i != verts.end(); i += 2) {text_bcube.union_with_pt(i->v);} // only need to include opposite corners
	float const width_scale(ct.get_sz_dim(!c.dim)/text_bcube.get_sz_dim(!c.dim)), height_scale(ct.dz()/text_bcube.dz());
	if (dot_product(normal, cross_product((verts[1].v - verts[0].v), (verts[2].v - verts[1].v))) < 0.0) {std::reverse(verts.begin(), verts.end());} // swap vertex winding order
	tid_nm_pair_t tex(FONT_TEXTURE_ID);
	if (c.color.A == 0.0) {tex.emissive = 1;}
	rgeom_mat_t &mat(get_material(tex, 0, 0, 1));
	color_wrapper const cw(colorRGBA(apply_light_color(c), 1.0)); // set alpha=1.0
	norm_comp const nc(normal);

	for (auto i = verts.begin(); i != verts.end(); ++i) {
		i->v[!c.dim] = i->v[!c.dim]*width_scale + ct.d[!c.dim][!ldir]; // line
		i->v.z       = i->v.z*height_scale + ct.z1(); // column
		mat.quad_verts.emplace_back(vert_norm_comp_tc(i->v, nc, i->t[0], i->t[1]), cw);
	}
}

void building_room_geom_t::add_window(room_object_t const &c, float tscale) {
	unsigned const skip_faces(get_skip_mask_for_xy(!c.dim) | EF_Z12); // only enable faces in dim
	cube_t window(c);
	swap(window.d[c.dim][0], window.d[c.dim][1]); // denormalized
	get_material(tid_nm_pair_t(get_bath_wind_tid(), tscale), 0).add_cube_to_verts(window, c.color, c.get_llc(), skip_faces); // no apply_light_color()
}

void building_room_geom_t::add_tub_outer(room_object_t const &c) {
	rgeom_mat_t &mat(get_material(untex_shad_mat, 1));
	colorRGBA const color(apply_light_color(c));
	mat.add_cube_to_verts(c, color, zero_vector, EF_Z12); // shadowed, no top/bottom faces
}

void building_room_geom_t::add_tv_picture(room_object_t const &c) {
	if (c.obj_id & 1) return; // TV is off half the time
	cube_t screen(c);
	screen.d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*0.35*c.get_sz_dim(c.dim);
	screen.expand_in_dim(!c.dim, -0.03*c.get_sz_dim(!c.dim)); // shrink the sides in
	screen.z1() += 0.09*c.dz();
	screen.z2() -= 0.04*c.dz();
	unsigned skip_faces(get_face_mask(c.dim, c.dir)); // only the face oriented outward
	tid_nm_pair_t tex(c.get_picture_tid(), 0.0);
	tex.emissive = 1;
	get_material(tex).add_cube_to_verts(screen, WHITE, c.get_llc(), skip_faces, !c.dim, !(c.dim ^ c.dir));
}

void building_room_geom_t::clear() {
	clear_materials();
	objs.clear();
	light_bcubes.clear();
	has_elevators = 0;
}
void building_room_geom_t::clear_materials() { // can be called to update textures, lighting state, etc.
	clear_static_vbos();
	mats_small.clear();
	mats_dynamic.clear();
	mats_lights.clear();
}
void building_room_geom_t::clear_static_vbos() { // used to clear pictures
	mats_static.clear();
	obj_model_insts.clear(); // these are associated with static VBOs
}

rgeom_mat_t &building_room_geom_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows, bool dynamic, bool small) {
	return (dynamic ? mats_dynamic : (small ? mats_small : mats_static)).get_material(tex, inc_shadows);
}
rgeom_mat_t &building_room_geom_t::get_wood_material(float tscale) {
	return get_material(get_tex_auto_nm(WOOD2_TEX, tscale), 1); // hard-coded for common material
}
colorRGBA get_textured_wood_color() {return WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX));}

colorRGBA room_object_t::get_color() const {
	switch (type) {
	case TYPE_TABLE:    return get_textured_wood_color();
	case TYPE_CHAIR: case TYPE_SM_CHAIR: return (color + get_textured_wood_color())*0.5; // 50% seat color / 50% wood legs color
	case TYPE_STAIR:    return LT_GRAY; // close enough
	case TYPE_ELEVATOR: return LT_BROWN; // ???
	case TYPE_RUG:      return texture_color(get_rug_tid());
	case TYPE_PICTURE:  return texture_color(get_picture_tid());
	case TYPE_WBOARD:   return WHITE;
	case TYPE_BCASE:    return get_textured_wood_color();
	case TYPE_DESK:     return get_textured_wood_color();
	case TYPE_BED:      return (color.modulate_with(texture_color(get_sheet_tid())) + get_textured_wood_color())*0.5; // half wood and half cloth
	default: return color; // TYPE_LIGHT, TYPE_TCAN, TYPE_BOOK, TYPE_BED
	}
	return color; // Note: probably should always set color so that we can return it here
}

void building_room_geom_t::create_static_vbos() {
	//highres_timer_t timer("Gen Room Geom"); // 2.1ms
	float const tscale(2.0/obj_scale);
	obj_model_insts.clear();

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible()) continue;
		assert(i->is_strictly_normalized());
		assert(i->type < NUM_TYPES);

		switch (i->type) {
		case TYPE_TABLE:   add_table   (*i, tscale); break;
		case TYPE_CHAIR:   add_chair   (*i, tscale); break;
		case TYPE_SM_CHAIR:add_chair   (*i, tscale); break;
		case TYPE_STAIR:   add_stair   (*i, tscale, tex_origin); break;
		case TYPE_RUG:     add_rug     (*i); break;
		case TYPE_PICTURE: add_picture (*i); break;
		case TYPE_WBOARD:  add_picture (*i); break;
		case TYPE_BOOK:    add_book    (*i, 1, 0); break;
		case TYPE_BCASE:   add_bookcase(*i, 1, 0, tscale, 0); break;
		case TYPE_DESK:    add_desk    (*i, tscale); break;
		case TYPE_TCAN:    add_trashcan(*i); break;
		case TYPE_BED:     add_bed     (*i, 1, 0, tscale); break;
		case TYPE_WINDOW:  add_window  (*i, tscale); break;
		case TYPE_TUB:     add_tub_outer(*i); break;
		case TYPE_TV:      add_tv_picture(*i); break;
		case TYPE_CUBICLE: add_cubicle (*i, tscale); break;
		case TYPE_STALL:   add_br_stall(*i); break;
		case TYPE_SIGN:    add_sign    (*i, 1, 0); break;
		case TYPE_PLANT:    break; // TODO
		case TYPE_ELEVATOR: break; // not handled here
		case TYPE_BLOCKER:  break; // not drawn
		case TYPE_COLLIDER: break; // not drawn
		default: break;
		} // end switch
		if (i->type >= TYPE_TOILET && i->type <= TYPE_COUCH) { // handle drawing of 3D models
			obj_model_insts.emplace_back((i - objs.begin()), (i->type + OBJ_MODEL_TOILET - TYPE_TOILET), i->color);
		}
	} // for i
	// Note: verts are temporary, but cubes are needed for things such as collision detection with the player and ray queries for indir lighting
	//timer_t timer2("Create VBOs"); // < 2ms
	mats_static.create_vbos();
}
void building_room_geom_t::create_small_static_vbos() {
	//highres_timer_t timer("Gen Room Geom Small"); // 1.3ms
	float const tscale(2.0/obj_scale);

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible()) continue;
		assert(i->is_strictly_normalized());
		assert(i->type < NUM_TYPES);

		switch (i->type) {
		case TYPE_BOOK:  add_book    (*i, 0, 1); break;
		case TYPE_BCASE: add_bookcase(*i, 0, 1, tscale, 0); break;
		case TYPE_BED:   add_bed     (*i, 0, 1, tscale); break;
		case TYPE_SIGN:  add_sign    (*i, 0, 1); break;
		default: break;
		}
	} // for i
	mats_small.create_vbos();
}
void building_room_geom_t::create_lights_vbos() {
	//highres_timer_t timer("Gen Room Geom Light"); // 0.3ms
	float const tscale(2.0/obj_scale);

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (i->is_visible() && i->type == TYPE_LIGHT) {add_light(*i, tscale);}
	}
	mats_lights.create_vbos();
}
void building_room_geom_t::create_dynamic_vbos() {
	if (!has_elevators) return; // currently only elevators are dynamic, can skip this step if there are no elevators

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible() || i->type != TYPE_ELEVATOR) continue; // only elevators for now
		add_elevator(*i, 2.0/obj_scale);
	}
	mats_dynamic.create_vbos();
}


void building_room_geom_t::draw(shader_t &s, vector3d const &xlate, bool shadow_only, bool inc_small, bool player_in_building) { // non-const because it creates the VBO
	if (empty()) return; // no geom
	unsigned const num_screenshot_tids(get_num_screenshot_tids());
	static int last_frame(0);
	static unsigned num_geom_this_frame(0); // used to limit per-frame geom gen time; doesn't apply to shadow pass, in case shadows are cached
	if (frame_counter > last_frame) {num_geom_this_frame = 0; last_frame = frame_counter;}

	if (lights_changed) {
		mats_lights.clear();
		lights_changed = 0;
	}
	if (has_pictures && num_pic_tids != num_screenshot_tids) {
		clear_static_vbos(); // user created a new screenshot texture, and this building has pictures - recreate room geom
		num_pic_tids = num_screenshot_tids;
	}
	if (mats_static.empty() && (shadow_only || num_geom_this_frame < MAX_ROOM_GEOM_GEN_PER_FRAME)) { // create static materials if needed
		create_static_vbos();
		++num_geom_this_frame;
	}
	if (inc_small && mats_small.empty() && (shadow_only || num_geom_this_frame < MAX_ROOM_GEOM_GEN_PER_FRAME)) { // create small materials if needed
		create_small_static_vbos();
		++num_geom_this_frame;
	}
	if (mats_lights .empty()) {create_lights_vbos ();} // create lights  materials if needed (no limit)
	if (mats_dynamic.empty()) {create_dynamic_vbos();} // create dynamic materials if needed (no limit)
	enable_blend(); // needed for rugs and book text
	mats_static .draw(s, shadow_only);
	mats_lights .draw(s, shadow_only);
	mats_dynamic.draw(s, shadow_only);
	if (inc_small) {mats_small.draw(s, shadow_only);}
	disable_blend();
	vbo_wrap_t::post_render();
	bool obj_drawn(0);

	// draw object models
	for (auto i = obj_model_insts.begin(); i != obj_model_insts.end(); ++i) {
		assert(i->obj_id < objs.size());
		auto const &obj(objs[i->obj_id]);
		if (!player_in_building && obj.is_interior()) continue; // don't draw objects in interior rooms if the player is outside the building (useful for office bathrooms)
		if (!shadow_only && !dist_less_than((camera_pdu.pos - xlate), obj.get_llc(), 100.0*obj.dz())) continue; // too far away
		if (!camera_pdu.cube_visible(obj + xlate)) continue; // VFC
		vector3d dir(zero_vector);
		dir[obj.dim] = (obj.dir ? 1.0 : -1.0);
		building_obj_model_loader.draw_model(s, obj.get_cube_center(), obj, dir, i->color, xlate, i->model_id, shadow_only, 0, 0);
		obj_drawn = 1;
	}
	if (obj_drawn) {check_mvm_update();} // needed after popping model transform matrix
}


