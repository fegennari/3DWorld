// 3D World - Building Interior Room Geometry Drawing
// by Frank Gennari 7/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#include "subdiv.h" // for sd_sphere_d
#include "profiler.h"
#include "scenery.h" // for s_plant

bool const ADD_BOOK_COVERS = 1;
bool const ADD_BOOK_TITLES = 1;
unsigned const MAX_ROOM_GEOM_GEN_PER_FRAME = 1;
colorRGBA const STAIRS_COLOR_TOP(0.7, 0.7, 0.7);
colorRGBA const STAIRS_COLOR_BOT(0.9, 0.9, 0.9);

vect_cube_t temp_cubes;
vector<room_object_t> temp_objects;
object_model_loader_t building_obj_model_loader;

extern bool camera_in_building;
extern int display_mode, frame_counter, player_in_closet;
extern float office_chair_rot_rate;
extern pos_dir_up camera_pdu;
extern building_t const *player_building;
extern room_object_t player_held_object;

int get_rand_screenshot_texture(unsigned rand_ix);
unsigned get_num_screenshot_tids();

void gen_text_verts(vector<vert_tc_t> &verts, point const &pos, string const &text, float tsize, vector3d const &column_dir, vector3d const &line_dir, bool use_quads=0);
string const &gen_book_title(unsigned rand_id, string *author, unsigned split_len);


unsigned get_face_mask(unsigned dim, bool dir) {return ~(1 << (2*(2-dim) + dir));} // skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2
unsigned get_skip_mask_for_xy(bool dim) {return (dim ? EF_Y12 : EF_X12);}
tid_nm_pair_t get_tex_auto_nm(int tid, float tscale=1.0) {return tid_nm_pair_t(tid, get_normal_map_for_bldg_tid(tid), tscale, tscale);}
int get_counter_tid    () {return get_texture_by_name("marble2.jpg");}
int get_paneling_nm_tid() {return get_texture_by_name("normal_maps/paneling_NRM.jpg", 1);}
int get_blinds_tid     () {return get_texture_by_name("interiors/blinds.jpg", 1, 0, 1, 8.0);} // use high aniso

vect_cube_t get_temp_cubes() {temp_cubes.clear(); return temp_cubes;}
vector<room_object_t> &get_temp_objects() {temp_objects.clear(); return temp_objects;}

// skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2 to match CSG cube flags
void rgeom_mat_t::add_cube_to_verts(cube_t const &c, colorRGBA const &color, point const &tex_origin,
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

void rgeom_mat_t::add_sphere_to_verts(cube_t const &c, colorRGBA const &color, bool low_detail, xform_matrix const *const matrix) {
	static vector<vert_norm_tc> cached_verts[2]; // high/low detail, reused across all calls
	vector<vert_norm_tc> &verts(cached_verts[low_detail]);

	if (verts.empty()) { // not yet created, create and cache verts
		sd_sphere_d sd(all_zeros, 1.0, (low_detail ? N_SPHERE_DIV/2 : N_SPHERE_DIV));
		sphere_point_norm spn;
		sd.gen_points_norms(spn);
		sd.get_quad_points(verts);
	}
	color_wrapper const cw(color);
	point const center(c.get_cube_center()), size(0.5*c.get_size());

	for (auto i = verts.begin(); i != verts.end(); ++i) {
		point pt(i->v*size);
		vector3d normal(i->n);
		if (matrix) {matrix->apply_to_vector3d(pt); matrix->apply_to_vector3d(normal);}
		quad_verts.emplace_back((pt + center), normal, i->t[0], i->t[1], cw);
	}
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

void rotate_vertst(vector<rgeom_mat_t::vertex_t> &verts, building_t const &building) {
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
		rotate_vertst(quad_verts, building);
		rotate_vertst(itri_verts, building);
	}
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
void building_materials_t::create_vbos(building_t const &building) {
	for (iterator m = begin(); m != end(); ++m) {m->create_vbo(building);}
}
void building_materials_t::draw(shader_t &s, bool shadow_only, bool reflection_pass) {
	for (iterator m = begin(); m != end(); ++m) {m->draw(s, shadow_only, reflection_pass);}
}

void get_tc_leg_cubes_abs_width(cube_t const &c, float leg_width, cube_t cubes[4]) {
	for (unsigned y = 0; y < 2; ++y) {
		for (unsigned x = 0; x < 2; ++x) {
			cube_t leg(c);
			leg.d[0][x] += (x ? -1.0f : 1.0f)*(c.dx() - leg_width);
			leg.d[1][y] += (y ? -1.0f : 1.0f)*(c.dy() - leg_width);
			cubes[2*y+x] = leg;
		}
	}
}
void get_tc_leg_cubes(cube_t const &c, float width, cube_t cubes[4]) {
	get_tc_leg_cubes_abs_width(c, get_tc_leg_width(c, width), cubes);
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
colorRGBA building_room_geom_t::apply_wood_light_color(room_object_t const &o) const {return apply_light_color(o, wood_color);}
colorRGBA apply_light_color(room_object_t const &o) {return apply_light_color(o, o.color);} // use object color

tid_nm_pair_t const untex_shad_mat(-1, 2.0); // make sure it's different from default tid_nm_pair_t so that it's not grouped with shadowed materials

void building_room_geom_t::add_table(room_object_t const &c, float tscale, float top_dz, float leg_width) { // 6 quads for top + 4 quads per leg = 22 quads = 88 verts
	cube_t top(c), legs_bcube(c);
	top.z1() += (1.0 - top_dz)*c.dz();
	legs_bcube.z2() = top.z1();
	colorRGBA const color(apply_wood_light_color(c));
	rgeom_mat_t &mat(get_wood_material(tscale));

	if (c.shape == SHAPE_CYLIN) { // cylindrical table
		vector3d const size(c.get_size());
		legs_bcube.expand_by_xy(-0.46*size);
		mat.add_vcylin_to_verts(top, color, 1, 1, 0, 0, 1.0, 1.0, 16.0, 2.0); // draw top and bottom with scaled side texture coords
		mat.add_vcylin_to_verts(legs_bcube, color, 1, 1, 0, 0, 1.0, 1.0, 1.0); // support
		cube_t feet(c);
		feet.z2() = c.z1() + 0.1*c.dz();
		feet.expand_by_xy(-0.2*size);

		for (unsigned d = 0; d < 2; ++d) { // add crossed feet
			cube_t foot(feet);
			foot.expand_in_dim(d, -0.27*size[d]);
			mat.add_cube_to_verts(foot, color, tex_origin, EF_Z1); // skip bottom surface
		}
	}
	else { // cube or short table
		assert(c.shape == SHAPE_CUBE || c.shape == SHAPE_SHORT);
		mat.add_cube_to_verts(top, color, c.get_llc()); // all faces drawn
		add_tc_legs(legs_bcube, color, leg_width, tscale);
	}
}

void get_chair_cubes(room_object_t const &c, cube_t cubes[3]) {
	float const height(c.dz()*((c.shape == SHAPE_SHORT) ? 1.333 : 1.0)); // effective height if the chair wasn't short
	cube_t seat(c), back(c), legs_bcube(c);
	seat.z1() += 0.32*height;
	seat.z2()  = back.z1() = seat.z1() + 0.07*height;
	legs_bcube.z2() = seat.z1();
	back.d[c.dim][c.dir] += 0.88f*(c.dir ? -1.0f : 1.0f)*c.get_sz_dim(c.dim);
	cubes[0] = seat; cubes[1] = back; cubes[2] = legs_bcube;
}
void building_room_geom_t::add_chair(room_object_t const &c, float tscale) { // 6 quads for seat + 5 quads for back + 4 quads per leg = 27 quads = 108 verts
	cube_t cubes[3]; // seat, back, legs_bcube
	get_chair_cubes(c, cubes);
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.2*tscale), 1).add_cube_to_verts(cubes[0], apply_light_color(c), c.get_llc()); // seat; all faces drawn
	colorRGBA const color(apply_wood_light_color(c));
	get_wood_material(tscale).add_cube_to_verts(cubes[1], color, c.get_llc(), EF_Z1); // back; skip bottom face
	add_tc_legs(cubes[2], color, 0.15, tscale); // legs
}

room_object_t get_dresser_middle(room_object_t const &c) {
	room_object_t middle(c);
	middle.z1() += 0.12*c.dz();
	middle.z2() -= 0.06*c.dz(); // at bottom of top surface
	middle.expand_by_xy(-0.5*get_tc_leg_width(c, 0.10)); // shrink by half leg width
	return middle;
}
void building_room_geom_t::add_dresser(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm) { // or nightstand
	if (inc_lg) {
		add_table(c, tscale, 0.06, 0.10);
		get_wood_material(tscale).add_cube_to_verts(get_dresser_middle(c), apply_wood_light_color(c), c.get_llc()); // all faces drawn
	}
	if (inc_sm) { // add drawers
		room_object_t middle(get_dresser_middle(c));
		middle.expand_in_dim(!c.dim, -0.5*get_tc_leg_width(c, 0.10));
		add_dresser_drawers(middle, tscale);
	}
}

void get_drawer_cubes(room_object_t const &c, vect_cube_t &drawers, bool front_only) {
	drawers.clear();
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	rgen.rand_mix();
	float const width(c.get_sz_dim(!c.dim)), depth(c.get_sz_dim(c.dim)), height(c.dz());
	bool is_lg(width > 2.0*height);
	unsigned const num_rows(2 + (rgen.rand() & 1)); // 2-3
	float const row_spacing(height/num_rows), door_thick(0.05*height), border(0.1*row_spacing), dir_sign(c.dir ? 1.0 : -1.0);
	float const drawer_extend(((c.type == TYPE_DESK) ? 0.5 : 0.8)*dir_sign*depth);
	cube_t d_row(c);
	d_row.d[ c.dim][!c.dir]  = c.d[c.dim][c.dir];
	d_row.d[ c.dim][ c.dir] += dir_sign*door_thick; // expand out a bit
	unsigned num_cols(1); // 1 for nightstand
	float vpos(c.z1());

	for (unsigned n = 0; n < num_rows; ++n) { // at most 12 drawers
		if (is_lg && (num_cols == 1 || rgen.rand_bool())) {num_cols = 2 + (rgen.rand() % 3);} // 2-4, 50% of the time keep same as prev row
		float const col_spacing(width/num_cols);
		float hpos(c.d[!c.dim][0]);
		set_cube_zvals(d_row, (vpos + border), (vpos + row_spacing - border));

		for (unsigned m = 0; m < num_cols; ++m) {
			cube_t drawer(d_row); // front part of the drawer
			drawer.d[!c.dim][0] = hpos + border;
			drawer.d[!c.dim][1] = hpos + col_spacing - border;

			if (c.drawer_flags & (1 << drawers.size())) { // make a drawer open
				drawer.d[c.dim][c.dir] += drawer_extend;
				if (front_only) {drawer.d[c.dim][!c.dir] += drawer_extend;} // translate the other side as well
			}
			drawers.push_back(drawer);
			hpos += col_spacing;
		} // for m
		vpos += row_spacing;
	} // for n
}

/*static*/ room_object_t building_room_geom_t::get_item_in_drawer(room_object_t const &c, cube_t const &drawer, unsigned drawer_ix) {
	assert(drawer_ix < 16);
	if (c.item_flags & (1U << drawer_ix)) {return room_object_t();} // item has been taken
	vector3d const sz(drawer.get_size());
	cube_t inside(drawer);
	inside.z1() += (0.2/0.9)*sz.z; // place on bottom surface
	rand_gen_t rgen;
	rgen.set_state((123*drawer_ix + 1), (456*c.room_id + 777*c.obj_id + 1));
	unsigned const type_ix(rgen.rand() % 10); // 0-9
	room_object_t obj; // starts as no item

	// TODO: add TYPE_PAPER, TYPE_PEN, TYPE_PENCIL, TYPE_BOOK, TYPE_BOTTLE, cell phone, wallet, money, etc.
	switch (type_ix) {
	case 0: case 2: case 3: case 4: // box - common
		obj = room_object_t(inside, TYPE_BOX, c.room_id, rgen.rand_bool(), rgen.rand_bool());
		
		for (unsigned d = 0; d < 2; ++d) {
			obj.d[d][0] += rgen.rand_uniform(0.1, 0.4)*sz[d];
			obj.d[d][1] -= rgen.rand_uniform(0.1, 0.4)*sz[d];
		}
		obj.z2() -= rgen.rand_uniform(0.1, 0.6)*sz.z;
		break;
	case 5: // key - rare
		obj = room_object_t(inside, TYPE_KEY, c.room_id, rgen.rand_bool(), rgen.rand_bool());
		obj.expand_in_dim( obj.dim, -0.40*sz[ obj.dim]); // long  dim
		obj.expand_in_dim(!obj.dim, -0.46*sz[!obj.dim]); // short dim
		for (unsigned d = 0; d < 2; ++d) {obj.translate_dim(d, 0.4*sz[d]*rgen.rand_uniform(-1.0, 1.0));}
		obj.z2() = obj.z1() + 0.05*sz.z;
		break;
	}
	obj.flags |= RO_FLAG_WAS_EXP;
	return obj;
}

void building_room_geom_t::add_dresser_drawers(room_object_t const &c, float tscale) { // or nightstand
	vect_cube_t &drawers(get_temp_cubes());
	get_drawer_cubes(c, drawers, 1); // front_only=1
	assert(drawers.size() <= 16); // we only have 16 bits to store drawer flags
	float const depth(c.get_sz_dim(c.dim)), height(c.dz()), drawer_thick(0.05*height), handle_thick(0.75*drawer_thick), dir_sign(c.dir ? 1.0 : -1.0), handle_width(0.07*height);
	get_metal_material(0, 0, 1); // ensure material exists so that door_mat reference is not invalidated
	rgeom_mat_t &drawer_mat(get_material(get_tex_auto_nm(WOOD2_TEX, 2.0*tscale), 1, 0, 1)); // shadowed, small=1
	rgeom_mat_t &handle_mat(get_metal_material(0, 0, 1)); // untextured, unshadowed, small=1
	colorRGBA const drawer_color(apply_light_color(c, WHITE)); // lighter color than dresser
	colorRGBA const handle_color(apply_light_color(c, GRAY_BLACK));
	unsigned const door_skip_faces(~get_face_mask(c.dim, !c.dir));
	vector<room_object_t> &objects(get_temp_objects());

	for (auto i = drawers.begin(); i != drawers.end(); ++i) {
		float const dwidth(i->get_sz_dim(!c.dim)), handle_shrink(0.5*dwidth - handle_width);
		unsigned door_skip_faces_mod(door_skip_faces);

		if (i->d[c.dim][!c.dir] != c.d[c.dim][c.dir]) { // drawer is not flush with front face, so it's open
			float const dheight(i->dz());
			cube_t drawer_body(*i);
			drawer_body.d[c.dim][!c.dir] = c. d[c.dim][ c.dir]; // flush with front
			drawer_body.d[c.dim][ c.dir] = i->d[c.dim][!c.dir]; // inside of drawer face
			drawer_body.expand_in_dim(!c.dim, -0.05*dwidth);
			drawer_body.expand_in_dim(2,      -0.05*dheight);
			cube_t bottom(drawer_body), left(drawer_body), right(drawer_body), back(drawer_body);
			left.z1() = right.z1() = bottom.z2() = drawer_body.z2() - 0.8*dheight;
			left.z2() = right.z2() = drawer_body.z2() - 0.1*dheight; // sides slightly shorter than the front and back
			left .d[!c.dim][1] -= 0.87*dwidth;
			right.d[!c.dim][0] += 0.87*dwidth;
			back.d[c.dim][ c.dir] = c.d[c.dim][c.dir] + 0.25f*dir_sign*drawer_thick; // flush with front face and narrow
			unsigned const skip_mask_front_back(get_skip_mask_for_xy(c.dim));
			colorRGBA const blr_color(drawer_color*0.4 + apply_wood_light_color(c)*0.4); // halfway between base and drawer colors, but slightly darker
			// swap the texture orientation of drawers to make them stand out more
			drawer_mat.add_cube_to_verts(bottom, blr_color, tex_origin,  skip_mask_front_back, 1);
			drawer_mat.add_cube_to_verts(left,   blr_color, tex_origin, (skip_mask_front_back | EF_Z1), 1);
			drawer_mat.add_cube_to_verts(right,  blr_color, tex_origin, (skip_mask_front_back | EF_Z1), 1);
			// draw inside face of back of drawer;
			// normally this wouldn't be visible here, but it's easier to drawn than holes for the drawers and it doesn't look as bad as doing nothing;
			// it would be better to cut a hole into the front of the desk for the drawer to slide into, but that seems to be difficult
			drawer_mat.add_cube_to_verts(back, drawer_color, tex_origin, get_face_mask(c.dim,c.dir), 1);
			door_skip_faces_mod = 0; // need to draw interior face
			room_object_t const obj(get_item_in_drawer(c, drawer_body, (i - drawers.begin())));
			if (obj.type != TYPE_NONE) {objects.push_back(obj);}
		}
		drawer_mat.add_cube_to_verts(*i, drawer_color, tex_origin, door_skip_faces_mod, 1); // swap the texture orientation of drawers to make them stand out more
		// add door handle
		cube_t handle(*i);
		handle.d[c.dim][!c.dir]  = i->d[c.dim][c.dir];
		handle.d[c.dim][ c.dir] += dir_sign*handle_thick; // expand out a bit
		handle.d[!c.dim][0] = i->d[!c.dim][0] + handle_shrink;
		handle.d[!c.dim][1] = i->d[!c.dim][1] - handle_shrink;
		handle.z1() = i->z1() + 0.8*i->dz();
		handle.z2() = handle.z1() + 0.1*i->dz();
		handle_mat.add_cube_to_verts(handle, handle_color, tex_origin, door_skip_faces); // same skip_faces
	} // for i
	add_small_static_objs_to_verts(objects); // add any objects that were found in open drawers; must be small static objects
}

tid_nm_pair_t get_scaled_wall_tex(tid_nm_pair_t const &wall_tex) {
	tid_nm_pair_t wall_tex_scaled(wall_tex);
	wall_tex_scaled.tscale_x *= 2.0;
	wall_tex_scaled.tscale_y *= 2.0;
	return wall_tex_scaled;
}

// cubes: front left, left side, front right, right side, door
void get_closet_cubes(room_object_t const &c, cube_t cubes[5]) {
	float const width(c.get_sz_dim(!c.dim)), depth(c.get_sz_dim(c.dim)), height(c.dz());
	bool const use_small_door(c.is_small_closet());
	float const wall_width(use_small_door ? 0.5*(width - 0.5*height) : 0.05*width), wall_shift(width - wall_width), wall_thick(WALL_THICK_VAL*(1.0f - FLOOR_THICK_VAL_HOUSE)*height);
	assert(wall_shift > 0.0);
	cube_t doors(c), walls[2] = {c, c}; // left, right
	walls[0].d[!c.dim][1] -= wall_shift;
	walls[1].d[!c.dim][0] += wall_shift;

	for (unsigned d = 0; d < 2; ++d) {
		cube_t front(walls[d]), side(walls[d]);
		front.d[c.dim][!c.dir] = side.d[c.dim][c.dir] = front.d[c.dim][c.dir] - (c.dir ? 1.0 : -1.0)*wall_thick;

		if (!(c.flags & (d ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO))) { // only need to draw side wall when not adjacent to room wall
			side.d[!c.dim][!d] += (d ? 1.0f : -1.0f)*(wall_width - wall_thick);
			cubes[2*d+1] = side;
		}
		cubes[2*d] = front;
		doors.d[!c.dim][d] = walls[d].d[!c.dim][!d]; // clip door to space between walls
	} // for d
	doors.d[c.dim][ c.dir] -= (c.dir ? 1.0 : -1.0)*0.2*wall_thick; // shift in slightly
	doors.d[c.dim][!c.dir] += (c.dir ? 1.0 : -1.0)*(depth - 0.8*wall_thick); // make it narrow
	cubes[4] = doors; // Note: this is for closed door; caller must handle open door
}
cube_t get_closet_interior_space(room_object_t const &c, cube_t const cubes[5]) {
	cube_t interior(c);
	if (!cubes[1].is_all_zeros()) {interior.d[!c.dim][0] = cubes[1].d[!c.dim][1];} // left  side (if wall exists)
	if (!cubes[3].is_all_zeros()) {interior.d[!c.dim][1] = cubes[3].d[!c.dim][0];} // right side (if wall exists)
	interior.d[ c.dim][c.dir] = cubes[2].d[c.dim][!c.dir]; // set to inside of front wall
	assert(interior.is_strictly_normalized());
	return interior;
}
void add_closet_objects(room_object_t const &c, vector<room_object_t> &objects) {
	cube_t ccubes[5]; // unused
	get_closet_cubes(c, ccubes);
	cube_t const interior(get_closet_interior_space(c, ccubes));
	float const depth(interior.get_sz_dim(c.dim)), box_sz(0.25*depth), window_vspacing(c.dz()*(1.0 + FLOOR_THICK_VAL_HOUSE));
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	unsigned const num_boxes((rgen.rand()%3) + (rgen.rand()%4)); // 0-5
	vect_cube_t &cubes(get_temp_cubes());
	room_object_t C(c);
	C.flags |= RO_FLAG_WAS_EXP;
	C.type   = TYPE_BOX;
	vector3d sz;
	point center;

	for (unsigned n = 0; n < num_boxes; ++n) {
		for (unsigned d = 0; d < 2; ++d) {
			sz    [d] = box_sz*rgen.rand_uniform(0.5, 1.0); // x,y half width
			center[d] = rgen.rand_uniform(interior.d[d][0]+sz[d], interior.d[d][1]-sz[d]); // randomly placed within the bounds of the shelf
		}
		C.set_from_point(center);
		set_cube_zvals(C, interior.z1(), (interior.z1() + box_sz*rgen.rand_uniform(0.8, 1.5)));
		C.expand_by_xy(sz);
		if (has_bcube_int(C, cubes)) continue; // intersects - just skip it, don't try another placement
		C.color = colorRGBA(rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0)); // add minor color variation
		C.dim   = rgen.rand_bool();
		C.dir   = rgen.rand_bool();
		objects.push_back(C);
		cubes.push_back(C);
	} // for n
	// add hanger rod
	float const hr_radius(0.015*window_vspacing);
	room_object_t hanger_rod(interior, TYPE_HANGER_ROD, c.room_id, c.dim, c.dir, (RO_FLAG_NOCOLL | RO_FLAG_INTERIOR));
	hanger_rod.z1() = c.z1() + 0.8*window_vspacing;
	hanger_rod.z2() = hanger_rod.z1() + 2.0*hr_radius;
	set_wall_width(hanger_rod, c.get_center_dim(c.dim), hr_radius, c.dim);
	objects.push_back(hanger_rod);
}

void building_room_geom_t::add_closet(room_object_t const &c, tid_nm_pair_t const &wall_tex, bool inc_lg, bool inc_sm) { // no lighting scale, houses only
	bool const use_small_door(c.is_small_closet()), draw_interior(c.is_open() || player_in_closet);
	cube_t cubes[5];
	get_closet_cubes(c, cubes);

	if (inc_lg) { // draw closet walls and doors
		rgeom_mat_t &wall_mat(get_material(get_scaled_wall_tex(wall_tex), 1));
		unsigned const skip_faces(~get_face_mask(c.dim, !c.dir) | EF_Z12); // skip top, bottom, and face that's against the wall

		for (unsigned d = 0; d < 2; ++d) {
			bool const adj_room_wall(c.flags & (d ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));
			unsigned const extra_skip_faces(adj_room_wall ? ~get_face_mask(!c.dim, d) : 0); // adjacent to room wall, skip that face
			// only need to draw side wall when not adjacent to room wall; skip front face of side wall
			if (!adj_room_wall) {wall_mat.add_cube_to_verts(cubes[2*d+1], c.color, tex_origin, (skip_faces | extra_skip_faces | ~get_face_mask(c.dim, c.dir)));}
			unsigned const front_wall_skip_flags((draw_interior ? EF_Z12 : skip_faces) | extra_skip_faces);
			wall_mat.add_cube_to_verts(cubes[2*d], c.color, tex_origin, front_wall_skip_flags); // Note: c.color should be wall color
		} // for d
		cube_t const &doors(cubes[4]);
		point const llc(doors.get_llc());

		if (use_small_door) { // small house closet door
			tid_nm_pair_t const tp(get_int_door_tid(), 0.0);

			if (c.is_open()) {
				float const door_width(doors.get_sz_dim(!c.dim)), door_thickness(doors.get_sz_dim(c.dim));
				cube_t door(doors);
				door.d[ c.dim][c.dir] += (c.dir ? 1.0 : -1.0)*door_width;
				door.d[!c.dim][1    ] -= (door_width - door_thickness);
				get_material(tp, 1)         .add_cube_to_verts(door, WHITE, llc, ~get_skip_mask_for_xy(!c.dim), c.dim, !c.dir); // draw front and back faces
				get_material(untex_shad_mat).add_cube_to_verts(door, WHITE, llc, ~get_skip_mask_for_xy( c.dim)); // draw edges untextured
			}
			else {
				unsigned const door_skip_faces(get_face_mask(c.dim, c.dir) & (player_in_closet ? get_face_mask(c.dim, !c.dir) : 0xFF));
				get_material(tp, 1).add_cube_to_verts(doors, WHITE, llc, door_skip_faces, !c.dim); // draw only front face, back face if player in closet
			}
		}
		else { // 4 panel folding door
			assert(!c.is_open()); // not yet handled
			float const doors_width(doors.get_sz_dim(!c.dim)), door_spacing(0.25*doors_width), door_gap(0.01*door_spacing);
			int const tid(get_rect_panel_tid());
			float tx(1.0/doors_width), ty(0.25/doors.dz());
			if (!c.dim) {swap(tx, ty);} // swap so that ty is always in Z
			tid_nm_pair_t const door_tex(tid, get_normal_map_for_bldg_tid(tid), tx, ty); // 4x1 panels
			rgeom_mat_t &door_mat(get_material(door_tex, 1));

			for (unsigned n = 0; n < 4; ++n) { // draw closet door in 4 parts
				cube_t door(doors);
				door.d[!c.dim][0] = doors.d[!c.dim][0] + n    *door_spacing + door_gap; // left edge
				door.d[!c.dim][1] = doors.d[!c.dim][0] + (n+1)*door_spacing - door_gap; // right edge
				door_mat.add_cube_to_verts(door, WHITE, llc, skip_faces);
			}
		}
	} // end inc_lg
	if (inc_sm) { // add wall trim for each side of the closet door
		float const height(c.dz()), window_vspacing(height*(1.0 + FLOOR_THICK_VAL_HOUSE));
		float const trim_height(0.04*window_vspacing), trim_thickness(0.1*WALL_THICK_VAL*window_vspacing);
		float const wall_thick(WALL_THICK_VAL*(1.0f - FLOOR_THICK_VAL_HOUSE)*height), trim_plus_wall_thick(trim_thickness + wall_thick);
		colorRGBA const trim_color(WHITE); // assume trim is white
		bool const draw_interior_trim(1 || draw_interior); // always enable so that we don't have to regenerate small geom when closet doors are opened or closed

		for (unsigned is_side = 0; is_side < 2; ++is_side) { // front wall, side wall
			for (unsigned d = 0; d < 2; ++d) {
				unsigned skip_faces((draw_interior_trim ? 0 : ~get_face_mask(c.dim, !c.dir)) | EF_Z1);
				cube_t trim;
				
				if (is_side) { // sides of closet
					if (c.flags & (d ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO)) continue; // adjacent to room wall, skip that wall
					trim = c;
					trim.d[!c.dim][!d]     = trim.d[!c.dim][d];
					trim.d[!c.dim][ d]    += (d     ? 1.0 : -1.0)*trim_thickness; // expand away from wall
					trim.d[ c.dim][c.dir] += (c.dir ? 1.0 : -1.0)*trim_thickness; // expand to cover the outside corner gap; doesn't really line up properly for angled ceiling trim though
				}
				else { // front of closet on sides of door
					trim = cubes[2*d]; // start with front wall on this side
					trim.d[!c.dim][!d    ] -= (d ? -1.0 : 1.0)*0.1*trim_thickness; // shrink slightly to avoid z-fighting with closet wall/door frame when the door is open
					trim.d[ c.dim][!c.dir]  = trim.d[c.dim][c.dir];
					trim.d[ c.dim][ c.dir] += (c.dir ? 1.0 : -1.0)*trim_thickness; // expand away from wall
					// expand to cover the outside corner gap if not along room wall, otherwise hide the face; doesn't really line up properly for angled ceiling trim though
					if (c.flags & (d ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO)) {skip_faces |= ~get_face_mask(!c.dim, d);} // disable hidden faces
					else {trim.d[!c.dim][d] += (d ? 1.0 : -1.0)*trim_thickness;}
				}
				bool const trim_dim(c.dim ^ bool(is_side)), trim_dir(is_side ? d : c.dir);
				cube_t btrim(trim); // bottom trim
				if (is_side) {btrim.d[!c.dim][!d    ] -= (d     ? 1.0 : -1.0)*trim_plus_wall_thick;} // expand on the inside of the closet
				else         {btrim.d[ c.dim][!c.dir] -= (c.dir ? 1.0 : -1.0)*trim_plus_wall_thick;} // expand on the inside of the closet
				btrim.z2() = c.z1() + trim_height;
				get_material(tid_nm_pair_t(), 0, 0, 1).add_cube_to_verts(btrim, trim_color, tex_origin, skip_faces); // is_small, untextured, no shadows; both interior and exterior
				set_cube_zvals(trim, c.z2()-trim_height, c.z2());
				add_wall_trim(room_object_t(trim, TYPE_WALL_TRIM, c.room_id, trim_dim, trim_dir, 0, 1.0, SHAPE_ANGLED, trim_color)); // ceiling trim, missing end caps; exterior only
			} // for d
		} // for is_side
	} // end inc_sm
	if (inc_sm /*&& draw_interior*/) { // Note: now always drawn to avoid recreating all small objects when the player opens/closes a closet door
		// Note: only for the closets the player has manually opened
		if (!(c.flags & RO_FLAG_EXPANDED)) { // add boxes if not expanded
			vector<room_object_t> &objects(get_temp_objects());
			add_closet_objects(c, objects);
			add_small_static_objs_to_verts(objects);
		}
	}
}
void building_room_geom_t::expand_closet(room_object_t const &c) {add_closet_objects(c, expanded_objs);}

void building_room_geom_t::add_hanger_rod(room_object_t const &c) { // is_small=1
	get_wood_material(1.0, 1, 0, 1).add_ortho_cylin_to_verts(c, LT_GRAY, !c.dim, 0, 0, 0, 0, 1.0, 1.0, 0.25, 1.0, 0, 16, 0.0, 1); // 16 sides, swap_txy=1
}

void building_room_geom_t::add_drain_pipe(room_object_t const &c) { // is_small=1
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(), 0, 0, 1)); // unshadowed, small
	mat.add_vcylin_to_verts(c, apply_light_color(c), 0, 0); // draw sides only
	mat.add_disk_to_verts(point(c.xc(), c.yc(), c.z2()), 0.5*c.dx(), 0, BLACK); // draw top as black
}

void building_room_geom_t::add_key(room_object_t const &c) { // is_small=1
	rgeom_mat_t &key_mat(get_metal_material(0, 0, 1)); // untextured, unshadowed, small=1
	key_mat.add_cube_to_verts(c, apply_light_color(c)); // TODO: make thisa key shape
}

int get_box_tid() {return get_texture_by_name("interiors/box.jpg");}
int get_crate_tid(room_object_t const &c) {return get_texture_by_name((c.obj_id & 1) ? "interiors/crate2.jpg" : "interiors/crate.jpg");}

void building_room_geom_t::add_crate(room_object_t const &c) { // is_small=1
	// Note: draw as "small", not because crates are small, but because they're only added to windowless rooms and can't be easily seen from outside a building
	get_material(tid_nm_pair_t(get_crate_tid(c), 0.0), 1, 0, 1).add_cube_to_verts(c, apply_light_color(c), zero_vector, EF_Z1); // skip bottom face (even for stacked crate?)
}

void building_room_geom_t::add_box(room_object_t const &c) { // is_small=1
	// Note: draw as "small", not because boxes are small, but because they're only added to windowless rooms and can't be easily seen from outside a building
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_box_tid(), get_texture_by_name("interiors/box_normal.jpg", 1), 0.0, 0.0), 1, 0, 1)); // is_small=1
	unsigned const verts_start(mat.quad_verts.size());
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, EF_Z1); // skip bottom face (even for stacked box?)
	assert(mat.quad_verts.size() == verts_start + 20); // there should be 5 quads (+z -x +x -y +y) / 20 verts (no -z)
	float const sz(2048), x1(12/sz), x2(576/sz), x3(1458/sz), y1(1-1667/sz), y2(1-1263/sz), y3(1-535/sz); //, x4(2032/sz), y4(1-128/sz); // Note: we don't use all parts of the texture
	mat.quad_verts[verts_start+0].set_tc(x1, y2); // z (top)
	mat.quad_verts[verts_start+1].set_tc(x2, y2);
	mat.quad_verts[verts_start+2].set_tc(x2, y3);
	mat.quad_verts[verts_start+3].set_tc(x1, y3);

	for (unsigned d = 0; d < 2; ++d) { // for each end
		unsigned const ix_shift((1 + 2*d)*c.dim); // needed to make sure the up icon actually faces up
		unsigned ix(verts_start + 4*d + 4);

		for (unsigned e = 0; e < 2; ++e) { // x, y
			bool const f(c.dim ^ bool(e));
			mat.quad_verts[ix+((0+ix_shift)&3)].set_tc(x2, (f ? y1 : y2));
			mat.quad_verts[ix+((1+ix_shift)&3)].set_tc(x3, (f ? y1 : y2));
			mat.quad_verts[ix+((2+ix_shift)&3)].set_tc(x3, (f ? y2 : y3));
			mat.quad_verts[ix+((3+ix_shift)&3)].set_tc(x2, (f ? y2 : y3));
			ix += 8; // skip the other face
		} // for e
	} // for d
}

void building_room_geom_t::add_paint_can(room_object_t const &c) {
	float const side_tscale_add(fract(11111*c.x1() + 22222*c.y1() + 33333*c.z1())); // somewhat random
	rgeom_mat_t &side_mat(get_material(tid_nm_pair_t(get_texture_by_name("interiors/paint_can_label.png")), 1, 0, 1)); // shadows, small
	side_mat.add_vcylin_to_verts(c, apply_light_color(c), 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 24, side_tscale_add); // draw sides only; random texture rotation
	point top(c.get_cube_center());
	top.z = c.z2();
	get_metal_material(1, 0, 1).add_disk_to_verts(top, 0.5*min(c.dx(), c.dy()), 0, apply_light_color(c, LT_GRAY)); // shadowed, specular metal; small=1
}

unsigned get_shelves_for_object(room_object_t const &c, cube_t shelves[4]) {
	unsigned const num_shelves(2 + (c.room_id % 3)); // 2-4 shelves
	float const thickness(0.02*c.dz()), bracket_thickness(0.8*thickness), z_step(c.dz()/(num_shelves + 1)); // include a space at the bottom
	cube_t shelf(c);
	shelf.z2() = shelf.z1() + thickness; // set shelf thickness
	shelf.d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*bracket_thickness; // leave space behind the shelf for brackets

	for (unsigned s = 0; s < num_shelves; ++s) {
		shelf.translate_dim(2, z_step); // move up one step
		shelves[s] = shelf; // record for later use
	}
	return num_shelves;
}
void add_if_not_intersecting(room_object_t const &obj, vector<room_object_t> &objects, vect_cube_t &cubes) {
	if (!has_bcube_int(obj, cubes)) {objects.push_back(obj); cubes.push_back(obj);}
}
void get_shelf_objects(room_object_t const &c_in, cube_t const shelves[4], unsigned num_shelves, vector<room_object_t> &objects) {
	room_object_t c(c_in);
	c.flags |= RO_FLAG_WAS_EXP;
	bool const is_house(c.is_house());
	vector3d const c_sz(c.get_size());
	float const dz(c_sz.z), width(c_sz[c.dim]), thickness(0.02*dz), bracket_thickness(0.75*thickness);
	float const z_step(dz/(num_shelves + 1)), shelf_clearance(z_step - thickness - bracket_thickness), sz_scale(is_house ? 0.5 : 1.0);
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	vect_cube_t &cubes(get_temp_cubes());

	for (unsigned s = 0; s < num_shelves; ++s) {
		cube_t const &S(shelves[s]);
		room_object_t C(c);
		vector3d sz;
		point center;
		cubes.clear();
		// add crates/boxes
		unsigned const num_boxes(rgen.rand() % (is_house ? 8 : 13)); // 0-12

		for (unsigned n = 0; n < num_boxes; ++n) {
			for (unsigned d = 0; d < 2; ++d) {
				sz[d] = 0.5*width*sz_scale*(is_house ? 1.5 : 1.0)*rgen.rand_uniform(0.45, 0.8); // x,y half width
				center[d] = rgen.rand_uniform(S.d[d][0]+sz[d], S.d[d][1]-sz[d]); // randomly placed within the bounds of the shelf
			}
			C.set_from_point(center);
			set_cube_zvals(C, S.z2(), (S.z2() + shelf_clearance*sz_scale*rgen.rand_uniform(0.4, 0.98)));
			C.expand_by_xy(sz);
			if (has_bcube_int(C, cubes)) continue; // intersects - just skip it, don't try another placement
			C.color = colorRGBA(rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0)); // add minor color variation
			C.dim   = c.dim ^ bool(rgen.rand()&3) ^ 1; // make the box label face outside 75% of the time

			if (is_house || rgen.rand_bool()) {C.type = TYPE_BOX; objects.push_back(C);}
			else {
				C.obj_id = rgen.rand(); // used to select crate texture
				C.type   = TYPE_CRATE;
				objects.push_back(C);
			}
			cubes.push_back(C);
		} // for n
		// add computers; what about monitors?
		bool const top_shelf(s+1 == num_shelves);
		unsigned const num_comps(rgen.rand() % (is_house ? 3 : 6)); // 0-5
		float const h_val(0.21*1.1*dz), cheight(0.75*h_val), cwidth(0.44*cheight), cdepth(0.9*cheight); // fixed AR=0.44 to match the texture
		sz[ c.dim] = 0.5*cdepth;
		sz[!c.dim] = 0.5*cwidth;
		C.dim  = c.dim; C.dir = !c.dir; // reset dim, flip dir
		C.type = TYPE_COMPUTER;

		if (2.1*sz.x < c_sz.x && 2.1*sz.y < c_sz.y && (top_shelf || cheight < shelf_clearance)) { // if it fits in all dims
			for (unsigned n = 0; n < num_comps; ++n) {
				for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(S.d[d][0]+sz[d], S.d[d][1]-sz[d]);} // randomly placed within the bounds of the shelf
				C.set_from_point(center);
				set_cube_zvals(C, S.z2(), S.z2()+cheight);
				C.expand_by_xy(sz);
				add_if_not_intersecting(C, objects, cubes);
			} // for n
		}
		// add keyboards
		unsigned const num_kbds(rgen.rand() % 5); // 0-4
		float const kbd_hwidth(0.7*0.5*1.1*2.0*h_val), kbd_depth(0.6*kbd_hwidth), kbd_height(0.06*kbd_hwidth);
		sz[ c.dim] = 0.5*kbd_depth;
		sz[!c.dim] = kbd_hwidth;
		C.type = TYPE_KEYBOARD;

		if (2.1*sz.x < c_sz.x && 2.1*sz.y < c_sz.y && (top_shelf || kbd_height < shelf_clearance)) { // if it fits in all dims
			for (unsigned n = 0; n < num_kbds; ++n) {
				for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(S.d[d][0]+sz[d], S.d[d][1]-sz[d]);} // randomly placed within the bounds of the shelf
				C.set_from_point(center);
				set_cube_zvals(C, S.z2(), S.z2()+kbd_height);
				C.expand_by_xy(sz);
				C.dir = rgen.rand_bool();
				add_if_not_intersecting(C, objects, cubes);
			} // for n
		}
		// add bottles
		unsigned const num_bottles(((rgen.rand()&3) == 0) ? 0 : (rgen.rand() % 9)); // 0-8, 75% chance
		C.dir  = C.dim = 0;
		C.type = TYPE_BOTTLE;

		for (unsigned n = 0; n < num_bottles; ++n) {
			// same as building_t::place_bottle_on_obj()
			float const bottle_height(z_step*rgen.rand_uniform(0.4, 0.7)), bottle_radius(z_step*rgen.rand_uniform(0.07, 0.11));
			if (min(c_sz.x, c_sz.y) < 5.0*bottle_radius) continue; // shelf not wide/deep enough to add this bottle
			for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform((S.d[d][0] + 2.0*bottle_radius), (S.d[d][1] - 2.0*bottle_radius));} // place at least 2*radius from edge
			C.set_from_sphere(center, bottle_radius);
			set_cube_zvals(C, S.z2(), S.z2()+bottle_height);
			C.set_as_bottle(rgen.rand() & 127); // no empty bottles - don't set 7th bit
			add_if_not_intersecting(C, objects, cubes);
		} // for n
		// add paint cans
		float const pc_height(0.64*z_step), pc_radius(0.28*z_step), pc_edge_spacing(1.1*pc_radius);

		if (2.1*pc_edge_spacing < min(c_sz.x, c_sz.y)) { // shelf is wide/deep enough for paint cans
			unsigned const num_pcans(((rgen.rand()&3) == 0) ? 0 : (rgen.rand() % 7)); // 0-6, 75% chance
			C.color = WHITE;
			C.type  = TYPE_PAINTCAN;
			C.shape = SHAPE_CYLIN;

			for (unsigned n = 0; n < num_pcans; ++n) {
				for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform((S.d[d][0] + pc_edge_spacing), (S.d[d][1] - pc_edge_spacing));} // place at least pc_edge_spacing from edge
				C.set_from_sphere(center, pc_radius);
				set_cube_zvals(C, S.z2(), S.z2()+pc_height);
				add_if_not_intersecting(C, objects, cubes);
			} // for n
			C.shape = SHAPE_CUBE; // reset for next object type
		}
		// add large balls to houses
		float const ball_radius(0.048*1.1*dz); // 4.7 inches

		if (is_house && 2.1*ball_radius < min(c_sz.x, c_sz.y)) { // shelf is wide/deep enough for paint cans
			unsigned const num_balls(rgen.rand() % 3); // 0-2
			C.color = WHITE;
			C.type  = TYPE_LG_BALL;
			C.shape = SHAPE_SPHERE;
			center.z = S.z2() + ball_radius;

			for (unsigned n = 0; n < num_balls; ++n) {
				for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform((S.d[d][0] + ball_radius), (S.d[d][1] - ball_radius));} // place at least ball_radius from edge
				C.set_from_sphere(center, ball_radius);
				C.item_flags = rgen.rand_bool(); // random type
				add_if_not_intersecting(C, objects, cubes);
			} // for n
			C.shape      = SHAPE_CUBE; // reset for next object type
			C.item_flags = 0;          // reset for next object type
		}
	} // for s
}
void building_room_geom_t::add_shelves(room_object_t const &c, float tscale) {
	// Note: draw as "small", not because shelves are small, but because they're only added to windowless rooms and can't be easily seen from outside a building
	// draw back in case it's against a window, even though that shouldn't happen
	bool const is_house(c.is_house());
	cube_t shelves[4]; // max number of shelves
	unsigned const num_shelves(get_shelves_for_object(c, shelves));
	// add wooden shelves
	unsigned const skip_faces(is_house ? 0 : ~get_face_mask(c.dim, c.dir)); // skip back face at wall if it's a house because it could be against a window (though it really shouldn't be)
	rgeom_mat_t &wood_mat(get_wood_material(tscale, 1, 0, 1)); // inc_shadows=1, dynamic=0, small=1
	colorRGBA const shelf_color(apply_light_color(c));

	for (unsigned s = 0; s < num_shelves; ++s) {
		wood_mat.add_cube_to_verts(shelves[s], shelf_color, c.get_llc(), skip_faces, !c.dim); // make wood grain horizontal
	}
	// add support brackets
	vector3d const c_sz(c.get_size());
	float const dz(c_sz.z), length(c_sz[!c.dim]), width(c_sz[c.dim]), thickness(0.02*dz), bracket_thickness(0.8*thickness);
	unsigned const num_brackets(2 + round_fp(0.5*length/dz));
	float const b_offset(0.05*dz), b_step((length - 2*b_offset)/(num_brackets-1)), bracket_width(1.8*thickness);
	rgeom_mat_t &metal_mat(get_metal_material(1, 0, 1)); // shadowed, specular metal; small=1
	colorRGBA const bracket_color(apply_light_color(c, LT_GRAY));

	for (unsigned s = 0; s < num_shelves; ++s) {
		cube_t bracket(shelves[s]);
		bracket.z2()  = bracket.z1(); // below the shelf
		bracket.z1() -= bracket_thickness;
		bracket.d[c.dim][!c.dir] -= (c.dir ? -1.0 : 1.0)*0.1*width; // shorten slightly
		bracket.d[!c.dim][1] = bracket.d[!c.dim][0] + bracket_width; // set width
		bracket.translate_dim(!c.dim, b_offset);

		for (unsigned b = 0; b < num_brackets; ++b) {
			metal_mat.add_cube_to_verts(bracket, bracket_color, zero_vector, (skip_faces | EF_Z2)); // skip top faces, maybe back

			if (s == 0) { // add vertical brackets on first shelf
				cube_t vbracket(bracket);
				set_cube_zvals(vbracket, c.z1(), c.z2());
				vbracket.d[c.dim][ c.dir] = c         .d[c.dim][c.dir] + (c.dir ? -1.0 : 1.0)*0.1*bracket_thickness; // nearly against the wall
				vbracket.d[c.dim][!c.dir] = shelves[s].d[c.dim][c.dir]; // against the shelf
				metal_mat.add_cube_to_verts(vbracket, bracket_color, zero_vector, (skip_faces | EF_Z12)); // skip top/bottom faces, maybe back
			}
			bracket.translate_dim(!c.dim, b_step);
		} // for b
	} // for s
	// add objects to the shelves
	if (c.flags & RO_FLAG_EXPANDED) return; // shelves have already been expanded, don't need to create contained objects below
	vector<room_object_t> &objects(get_temp_objects());
	get_shelf_objects(c, shelves, num_shelves, objects);
	add_small_static_objs_to_verts(objects);
}
void building_room_geom_t::expand_shelves(room_object_t const &c) {
	cube_t shelves[4]; // max number of shelves
	unsigned const num_shelves(get_shelves_for_object(c, shelves));
	get_shelf_objects(c, shelves, num_shelves, expanded_objs);
}

void building_room_geom_t::add_keyboard(room_object_t const &c) {
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_texture_by_name("interiors/keyboard.jpg"), 0.0), 1, 0, 1)); // shadows, small
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, ~EF_Z2, c.dim, (c.dim ^ c.dir ^ 1), c.dir); // top face only
	get_material(tid_nm_pair_t(), 0, 0, 1).add_cube_to_verts(c, apply_light_color(c, BKGRAY), zero_vector, EF_Z12); // sides, no shadows, small
}

void building_room_geom_t::add_obj_with_front_texture(room_object_t const &c, string const &texture_name, colorRGBA const &sides_color, bool is_small) {
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_texture_by_name(texture_name), 0.0), 1, 0, is_small)); // shadows
	unsigned const front_mask(get_face_mask(c.dim, c.dir));
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, front_mask, !c.dim, (c.dim ^ c.dir ^ 1), 0); // front face only
	get_material(tid_nm_pair_t(), 1, 0, is_small).add_cube_to_verts(c, apply_light_color(c, sides_color), zero_vector, ~front_mask); // sides, shadows
}
void building_room_geom_t::add_computer(room_object_t const &c) {add_obj_with_front_texture(c, "interiors/computer.jpg",  BKGRAY, 1);} // is_small=1
void building_room_geom_t::add_mwave   (room_object_t const &c) {add_obj_with_front_texture(c, "interiors/microwave.jpg", GRAY,   0);} // is_small=0

void building_room_geom_t::add_mirror(room_object_t const &c) {
	tid_nm_pair_t tp(REFLECTION_TEXTURE_ID, 0.0);
	if (ENABLE_MIRROR_REFLECTIONS) {tp.emissive = 1.0;}
	get_material(tp, 0).add_cube_to_verts(c, c.color, zero_vector, get_face_mask(c.dim, c.dir), !c.dim); // draw only the front face
	get_material(tid_nm_pair_t(), 0).add_cube_to_verts(c, apply_light_color(c), zero_vector, get_skip_mask_for_xy(c.dim)); // draw only the sides untextured
}

void building_room_geom_t::add_shower(room_object_t const &c, float tscale) {
	bool const xdir(c.dim), ydir(c.dir), dirs[2] = {xdir, ydir}; // placed in this corner
	vector3d const sz(c.get_size());
	float const signs[2] = {(xdir ? -1.0f : 1.0f), (ydir ? -1.0f : 1.0f)};
	colorRGBA const color(apply_light_color(c));
	// add tile material along walls and floor
	int const skip_faces[2] = {(EF_Z1 | (xdir ? EF_X2 : EF_X1)), (EF_Z1 | (ydir ? EF_Y2 : EF_Y1))};
	rgeom_mat_t &tile_mat(get_material(tid_nm_pair_t(get_texture_by_name("bathroom_tile.jpg"), 2.5*tscale), 0)); // no shadows
	cube_t bottom(c), sides[2] = {c, c};
	bottom.z2() = c.z1() + 0.025*sz.z;
	tile_mat.add_cube_to_verts(bottom,   color, zero_vector, (skip_faces[0] | skip_faces[1]));

	for (unsigned d = 0; d < 2; ++d) {
		sides[d].d[d][!dirs[d]] -= signs[d]*0.98*sz[d];
		sides[d].z1() = bottom.z2();
		tile_mat.add_cube_to_verts(sides[d], color, zero_vector, skip_faces[d]);
	}
	// add metal frame around glass
	colorRGBA const metal_color(apply_light_color(c, GRAY));
	rgeom_mat_t &metal_mat(get_metal_material(1)); // shadowed, specular metal
	cube_t fc(c); // corner frame
	set_cube_zvals(fc, bottom.z2(), (c.z2() - 0.05*sz.z)); // slightly shorter than tile
	cube_t fxy[2] = {fc, fc};
	float const glass_bot(fc.z1() + 0.02*sz.z), glass_top(fc.z2() - 0.02*sz.z);

	for (unsigned d = 0; d < 2; ++d) {
		cube_t &f(fxy[d]);
		f.d[ d][ dirs[ d]]  = sides[d].d[d][!dirs[d]];
		f.d[ d][!dirs[ d]]  = sides[d].d[d][!dirs[d]] + signs[d]*0.04*sz[d];
		f.d[!d][ dirs[!d]] += signs[!d]*0.94*sz[!d];
		f.d[!d][!dirs[!d]] -= signs[!d]*0.02*sz[!d];
		metal_mat.add_cube_to_verts(f, metal_color, zero_vector, skip_faces[d]);
		fc.d[!d][0] = f.d[!d][0]; fc.d[!d][1] = f.d[!d][1];
	}
	metal_mat.add_cube_to_verts(fc, metal_color, zero_vector, EF_Z1);

	for (unsigned d = 0; d < 2; ++d) { // add top and bottom bars; these overlap with vertical frame cubes, but it should be okay and simplifies the math
		unsigned const bar_skip_faces(get_skip_mask_for_xy(d));
		cube_t tb_bars(fxy[d]);
		tb_bars.union_with_cube(fc);
		cube_t bot_bar(tb_bars), top_bar(tb_bars);
		bot_bar.z2() = glass_bot;
		top_bar.z1() = glass_top;
		metal_mat.add_cube_to_verts(bot_bar, metal_color, zero_vector, bar_skip_faces); // the track
		metal_mat.add_cube_to_verts(top_bar, metal_color, zero_vector, bar_skip_faces);
	}
	// add shower head
	float const radius(0.5f*(sz.x + sz.y));
	bool const head_dim(sz.y < sz.x);
	point start_pos, end_pos, base_pos, head_pos;
	start_pos.z = c.z1() + 0.75*sz.z;
	start_pos[ head_dim] = sides[head_dim].d[head_dim][!dirs[head_dim]];
	start_pos[!head_dim] = c.get_center_dim(!head_dim);
	base_pos = start_pos;
	base_pos[head_dim] += signs[head_dim]*0.06*sz[head_dim]; // move out from the wall
	end_pos  = base_pos;
	end_pos[head_dim] += signs[head_dim]*0.02*sz[head_dim]; // move out from the wall a bit more
	head_pos = base_pos;
	head_pos.z -= 0.05*sz.z;
	head_pos[ head_dim] += signs[head_dim]*0.09*sz[head_dim];
	metal_mat.add_cylin_to_verts(base_pos,  head_pos, 0.02*radius, 0.10*radius, metal_color, 0, 1); // draw top/wide end only
	metal_mat.add_cylin_to_verts(start_pos, end_pos,  0.02*radius, 0.02*radius, metal_color, 0, 0); // no ends
	// add door handle
	bool const hdim(c.dx() < c.dy()), hdir(dirs[hdim]); // large dim
	float const frame_width(fc.dx()), door_width(sz[!hdim]);
	cube_t handle(c);
	handle.d[ hdim][ hdir] = c.d[hdim][!hdir]     - (hdir ? -1.0 : 1.0)*(0.19*frame_width + 0.02*sz[hdir]); // place on the glass but slightly offset
	handle.d[ hdim][!hdir] = handle.d[hdim][hdir] + (hdir ? -1.0 : 1.0)*1.0*frame_width;
	handle.d[!hdim][0] += 0.20*door_width;
	handle.d[!hdim][1] -= 0.77*door_width;
	handle.z1() += 0.48*sz.z;
	handle.z2() -= 0.42*sz.z;
	metal_mat.add_cube_to_verts(handle, metal_color, zero_vector, 0); // draw all faces
	// add drain
	cube_t drain;
	drain.set_from_point(bottom.get_cube_center());
	set_cube_zvals(drain, bottom.z2(), (bottom.z2() + 0.05*bottom.dz())); // very small height
	drain.expand_by_xy(0.06*radius); // set radius
	get_material(tid_nm_pair_t(MANHOLE_TEX, 0.0), 0).add_vcylin_to_verts(drain, metal_color, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 1); // draw top only, not sides, unshadowed
	// add transparent glass
	colorRGBA const glass_color(apply_light_color(c, colorRGBA(1.0, 1.0, 1.0, 0.25)));
	rgeom_mat_t &glass_mat(get_material(tid_nm_pair_t(), 0, 0, 0, 1)); // no shadows; transparent=1

	for (unsigned d = 0; d < 2; ++d) {
		cube_t glass(fc);
		glass.z1()  = glass_bot;
		glass.z2()  = glass_top;
		glass.z1() += 0.002*sz.z; // to prevent z-fighting
		glass.z2() -= 0.002*sz.z; // to prevent z-fighting
		glass.d[d][!dirs[d]] = glass. d[d][ dirs[d]]; // remove overlap with frame
		glass.d[d][ dirs[d]] = fxy[d].d[d][!dirs[d]];
		glass.expand_in_dim( d, -0.01*frame_width); // to prevent z-fighting
		glass.expand_in_dim(!d, -0.20*frame_width);
		glass_mat.add_cube_to_verts(glass, glass_color, zero_vector, 0, 0, 0, 0, 1); // inside surface, inverted
		glass_mat.add_cube_to_verts(glass, glass_color, zero_vector, (EF_Z1 | (d ? EF_Y12 : EF_X12))); // outside surface
	}
}

void building_room_geom_t::add_bottle(room_object_t const &c) {
	// obj_id: bits 1-2 for type, bits 6-7 for emptiness, bit 6 for cap color
	// for now, no texture, but could use a bottle label texture for the central cylinder
	unsigned const bottle_ndiv = 16; // use smaller ndiv to reduce vertex count
	tid_nm_pair_t tex;
	tex.set_specular(0.5, 80.0);
	rgeom_mat_t &mat(get_material(tex, 1, 0, 1)); // inc_shadows=1, dynamic=0, small=1
	colorRGBA const color(apply_light_color(c)), cap_colors[2] = {LT_GRAY, GOLD};
	vector3d const sz(c.get_size());
	unsigned const dim(get_max_dim(sz)), dim1((dim+1)%3), dim2((dim+2)%3);
	bool const is_empty(c.is_bottle_empty()), add_bottom(max(sz.x, sz.y) > sz.z); // add bottom if bottle is on its side
	float const dir_sign(c.dir ? -1.0 : 1.0), radius(0.25f*(sz[dim1] + sz[dim2])); // base should be square (default/avg radius is 0.15*height)
	cube_t sphere(c), main_cylin(c), top_cylin(c);
	sphere.d[dim][ c.dir] = c.d[dim][c.dir] + dir_sign*0.5*sz[dim];
	sphere.d[dim][!c.dir] = sphere.d[dim][c.dir] + dir_sign*0.3*sz.z;
	main_cylin.d[dim][!c.dir] = sphere.d[dim][c.dir] + dir_sign*0.15*sz.z;
	top_cylin .d[dim][ c.dir] = main_cylin.d[dim][!c.dir]; // there will be some intersection, but that should be okay
	top_cylin.expand_in_dim(dim1, -0.32*sz[dim1]); // smaller radius
	top_cylin.expand_in_dim(dim2, -0.32*sz[dim2]); // smaller radius
	cube_t cap(top_cylin);
	top_cylin.d[dim][!c.dir] = cap.d[dim][c.dir] = c.d[dim][!c.dir] - dir_sign*0.08*sz[dim]; // set cap thickness
	cap.expand_in_dim(dim1, 0.034*sz[dim1]); // slightly larger radius
	cap.expand_in_dim(dim2, 0.034*sz[dim2]); // slightly larger radius
	mat.add_sphere_to_verts(sphere, color, 1); // low_detail=1
	mat.add_ortho_cylin_to_verts(main_cylin, color, dim, (add_bottom && !c.dir), (add_bottom && c.dir), 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bottle_ndiv);
	mat.add_ortho_cylin_to_verts(top_cylin,  color, dim, 0, is_empty, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bottle_ndiv); // draw neck of bottle; draw top if empty

	if (!is_empty) { // draw cap if nonempty
		mat.add_ortho_cylin_to_verts(cap, apply_light_color(c, cap_colors[bool(c.obj_id & 64)]), dim, c.dir, !c.dir, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bottle_ndiv);
	}
	// Note: we could add a bottom sphere to make it a capsule, then translate below the surface in -z to flatten the bottom
	main_cylin.expand_in_dim(dim1, 0.02*radius); // expand slightly in radius
	main_cylin.expand_in_dim(dim2, 0.02*radius); // expand slightly in radius
	main_cylin.d[dim][c.dir] += dir_sign*0.3*sz.z; main_cylin.d[dim][!c.dir] -= dir_sign*0.15*sz.z; // shrink in length
	get_material(tid_nm_pair_t(), 0, 0, 1).add_ortho_cylin_to_verts(main_cylin, apply_light_color(c, WHITE), dim, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bottle_ndiv); // white label
}

void building_room_geom_t::add_paper(room_object_t const &c) {
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(c.get_paper_tid(), 0.0), 0, 0, 1)); // map texture to quad
	unsigned const qv_start(mat.quad_verts.size());
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, ~EF_Z2, c.dim, (c.dim ^ c.dir ^ 1), c.dir); // unshadowed, top face only, with proper orient
	
	if (c.flags & RO_FLAG_RAND_ROT) { // add slight rotation to misalign the paper
		float const angle((PI/8.0)*(fract(123.456*c.obj_id) - 0.5));
		rotate_verts(mat.quad_verts, plus_z, angle, c.get_cube_center(), qv_start);
	}
}

void building_room_geom_t::add_pen_pencil(room_object_t const &c_) {
	room_object_t c(c_);
	if (!c.dir) {swap(c.d[c.dim][0], c.d[c.dim][1]); c.dir = 1;} // put in canonical orientation; okay if denormalized
	colorRGBA const color(apply_light_color(c));
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(), 0, 0, 1)); // unshadowed, small
	bool const is_pen(c.type == TYPE_PEN);
	float const length(c.get_sz_dim(c.dim));
	cube_t body(c), point(c);
	body.d[c.dim][1] = point.d[c.dim][0] = c.d[c.dim][1] - (is_pen ? 0.08 : 0.10)*length;
	// non-AA rotation/direction?

	if (is_pen) { // point is at the top
		mat.add_ortho_cylin_to_verts(body,  color, c.dim, 1, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 16); // 16-sided cylinder
		mat.add_ortho_cylin_to_verts(point, color, c.dim, 0, 1, 0, 0, 1.0, 0.2, 1.0, 1.0, 0, 16); // 16-sided truncated cone
	}
	else { // eraser is at the bottom and point is at the top
		colorRGBA const lt_wood(1.0, 0.8, 0.6);
		cube_t end_part(body), lead(point);
		point.d[c.dim][1] = lead    .d[c.dim][0] = c.d[c.dim][1] - 0.03*length;
		body .d[c.dim][0] = end_part.d[c.dim][1] = c.d[c.dim][0] + 0.09*length;
		cube_t eraser(end_part), metal(end_part);
		metal.d[c.dim][0] = eraser.d[c.dim][1] = c.d[c.dim][0] + 0.04*length;
		mat.add_ortho_cylin_to_verts(body,   color,                          c.dim, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 6); // 6-sided cylinder, should really be made flat sided
		mat.add_ortho_cylin_to_verts(point,  apply_light_color(c, lt_wood),  c.dim, 0, 0, 0, 0, 1.0, 0.3, 1.0, 1.0, 0, 12); // 12-sided truncated cone
		mat.add_ortho_cylin_to_verts(lead,   apply_light_color(c, DK_GRAY),  c.dim, 0, 0, 0, 0, 0.3, 0.0, 1.0, 1.0, 0, 12); // 12-sided cone
		mat.add_ortho_cylin_to_verts(metal,  apply_light_color(c, DK_GREEN), c.dim, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 12); // 12-sided cylinder
		mat.add_ortho_cylin_to_verts(eraser, apply_light_color(c, PINK),     c.dim, 1, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 12); // 12-sided cylinder
	}
}

void building_room_geom_t::add_flooring(room_object_t const &c, float tscale) {
	get_material(tid_nm_pair_t(MARBLE_TEX, 0.8*tscale), 1).add_cube_to_verts(c, apply_light_color(c), tex_origin, ~EF_Z2); // top face only
}

void building_room_geom_t::add_wall_trim(room_object_t const &c) {
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(), 0, 0, 1)); // inc_shadows=0, dynamic=0, small=1

	if (c.shape == SHAPE_ANGLED) { // single quad
		point pts[4];
		pts[0][!c.dim] = pts[1][!c.dim] = c.d[!c.dim][0];
		pts[2][!c.dim] = pts[3][!c.dim] = c.d[!c.dim][1];
		pts[0][ c.dim] = pts[3][ c.dim] = c.d[ c.dim][!c.dir];
		pts[1][ c.dim] = pts[2][ c.dim] = c.d[ c.dim][ c.dir];
		pts[0].z = pts[3].z = c.z1();
		pts[1].z = pts[2].z = c.z2();
		if (c.dir ^ c.dim) {swap(pts[0], pts[3]); swap(pts[1], pts[2]);} // change winding order/normal sign
		rgeom_mat_t::vertex_t v;
		v.set_norm(get_poly_norm(pts));
		v.set_c4(c.color);
		float const tcs[2][4] = {{0,0,1,1}, {0,1,1,0}};

		for (unsigned n = 0; n < 4; ++n) {
			v.v = pts[n];
			v.t[0] = tcs[0][n];
			v.t[1] = tcs[1][n];
			mat.quad_verts.push_back(v);
		}
	}
	else { // cube
		unsigned skip_faces(0);
		if      (c.shape == SHAPE_TALL ) {skip_faces = 0;} // door/window side trim
		else if (c.shape == SHAPE_SHORT) {skip_faces = get_skip_mask_for_xy(!c.dim);} // door top trim: skip ends
		else                             {skip_faces = get_skip_mask_for_xy(!c.dim) | EF_Z1;} // wall trim: skip bottom surface and short sides
		if (c.flags & RO_FLAG_ADJ_LO) {skip_faces |= ~get_face_mask(c.dim, 0);}
		if (c.flags & RO_FLAG_ADJ_HI) {skip_faces |= ~get_face_mask(c.dim, 1);}
		skip_faces |= ((c.flags & RO_FLAG_ADJ_BOT) ? EF_Z1 : 0) | ((c.flags & RO_FLAG_ADJ_TOP) ? EF_Z2 : 0);
		mat.add_cube_to_verts(c, c.color, tex_origin, skip_faces); // is_small, untextured, no shadows, not light scale
	}
}

void building_room_geom_t::add_blinds(room_object_t const &c) {
	bool const vertical(!(c.flags & RO_FLAG_HANGING));
	colorRGBA const color(c.color); // room color not applied as it looks wrong when viewed from outside the building
	// fit the texture to the cube; blinds have a fixed number of slats that compress when they are shortened
	// should these be partially transparent/backlit like bathroom windows? I guess not, most blinds are plastic or wood rather than material
	int const nm_tid(get_texture_by_name("interiors/blinds_hn.jpg", 1, 0, 1, 8.0)); // use high aniso
	float tx(vertical ? 1.0/c.dz() : 0.0), ty(vertical ? 0.5/c.get_sz_dim(!c.dim) : 0.0);
	if (c.dim) {swap(tx, ty);}
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_blinds_tid(), nm_tid, tx, ty), 1));
	unsigned df1(~get_skip_mask_for_xy(!c.dim)), df2(~EF_Z12);
	if (vertical) {swap(df1, df2);} // swap sides vs. top/bottom
	vector3d const llc(c.get_llc());
	mat.add_cube_to_verts(c, color, llc, ~get_skip_mask_for_xy(c.dim), (c.dim ^ vertical ^ 1)); // draw front and back faces
	mat.add_cube_to_verts(c, color, llc, df1, (c.dim ^ vertical)); // draw sides / top and bottom
	get_material(untex_shad_mat, 1).add_cube_to_verts(c, texture_color(get_blinds_tid()).modulate_with(color), llc, df2); // draw top and bottom / front and back untextured
}

float get_railing_height(room_object_t const &c) {
	bool const is_u_stairs(c.flags & (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI));
	return (is_u_stairs ? 0.70 : 0.35)*c.dz(); // use a larger relative height for lo/hi railings on U-shaped stairs
}
cylinder_3dw get_railing_cylinder(room_object_t const &c) {
	float const radius(0.5*c.get_sz_dim(!c.dim)), center(c.get_center_dim(!c.dim)), height(get_railing_height(c));
	point p[2];

	for (unsigned d = 0; d < 2; ++d) {
		p[d].z = c.d[2][d] + height;
		p[d][!c.dim] = center;
		p[d][ c.dim] = c.d[c.dim][c.dir^bool(d)^1];
	}
	return cylinder_3dw(p[0], p[1], radius, radius);
}
void building_room_geom_t::add_railing(room_object_t const &c) {
	cylinder_3dw const railing(get_railing_cylinder(c));
	bool const is_u_stairs(c.flags & (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI));
	float const pole_radius(0.75*railing.r1), length(c.get_sz_dim(c.dim)), height(get_railing_height(c));
	tid_nm_pair_t tex;
	tex.set_specular(0.7, 70.0);
	rgeom_mat_t &mat(get_material(tex, 1, 0, 1)); // inc_shadows=1, dynamic=0, small=1
	mat.add_cylin_to_verts(railing.p1, railing.p2, railing.r1, railing.r2, c.color, 1, 1); // draw sloped railing with both ends

	if (!is_u_stairs) {
		for (unsigned d = 0; d < 2; ++d) { // add the two vertical poles
			point pt(d ? railing.p2 : railing.p1);
			pt[c.dim] += (c.dir^bool(d) ? 1.0 : -1.0)*0.01*length; // shift slightly inward toward the center
			float const hscale(d ? 1.25 : 1.0); // shorten for lower end, which rests on the step
			point const p1(pt - vector3d(0, 0, hscale*height)), p2(pt - vector3d(0, 0, 0.02*(d ? 1.0 : -1.0)*height));
			mat.add_cylin_to_verts(p1, p2, pole_radius, pole_radius, c.color, 0, 0); // no top or bottom
		}
		if (c.is_open()) { // add balusters
			unsigned const num(NUM_STAIRS_PER_FLOOR - 1);
			float const step_sz(1.0/(num+1)), radius(0.75*pole_radius), bot_radius(0.85*pole_radius);
			vector3d const delta(0, 0, -height);

			for (unsigned n = 0; n < num; ++n) {
				float const t((n+1)*step_sz);
				point const pt(t*railing.p1 + (1.0 - t)*railing.p2);
				mat.add_cylin_to_verts((pt + delta), pt, radius, radius, c.color, 0, 0, 0, 0, 1.0, 1.0, 0, 16); // only 16 sides, no top or bottom
			}
			mat.add_cylin_to_verts((railing.p1 + delta), (railing.p2 + delta), bot_radius, bot_radius, c.color, 1, 1); // bottom bar with both ends
		}
	}
}

void building_room_geom_t::add_stair(room_object_t const &c, float tscale, vector3d const &tex_origin) { // Note: no room lighting color atten
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(MARBLE_TEX, 1.5*tscale), 1));
	float const width(c.get_sz_dim(!c.dim)); // use width as a size reference because this is constant for a set of stairs and in a relative small range
	cube_t top(c), bot(c);
	bot.z2() = top.z1() = c.z2() - min(0.025*width, 0.25*c.dz()); // set top thickness
	top.d[c.dim][!c.dir] += (c.dir ? -1.0 : 1.0)*0.0125*width; // extension
	top.expand_in_dim(!c.dim, 0.01*width); // make slightly wider
	mat.add_cube_to_verts(top, STAIRS_COLOR_TOP, tex_origin); // all faces drawn
	mat.add_cube_to_verts(bot, STAIRS_COLOR_BOT, tex_origin, EF_Z2); // skip top face
}

void building_room_geom_t::add_stairs_wall(room_object_t const &c, vector3d const &tex_origin, tid_nm_pair_t const &wall_tex) { // Note: no room lighting color atten
	get_material(get_scaled_wall_tex(wall_tex), 1).add_cube_to_verts(c, c.color, tex_origin); // all faces drawn
}

void building_room_geom_t::add_elevator(room_object_t const &c, float tscale) {
	// elevator car, all materials are dynamic; no lighting scale
	float const thickness(1.005*0.5*FLOOR_THICK_VAL_OFFICE*c.dz());
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
	rgeom_mat_t &paneling_mat(get_material(paneling, 1, 1));
	paneling_mat.add_cube_to_verts(back, WHITE, tex_origin, front_face_mask, !c.dim);

	for (unsigned d = 0; d < 2; ++d) { // side walls
		cube_t side(c);
		side.d[!c.dim][!d] = side.d[!c.dim][d] + (d ? -1.0 : 1.0)*thickness;
		paneling_mat.add_cube_to_verts(side, WHITE, tex_origin, get_face_mask(!c.dim, !d), c.dim);
	}
}

void building_room_geom_t::add_light(room_object_t const &c, float tscale) {
	// Note: need to use a different texture (or -1) for is_on because emissive flag alone does not cause a material change
	bool const is_on(c.is_lit());
	tid_nm_pair_t tp(((is_on || c.shape == SHAPE_SPHERE) ? (int)WHITE_TEX : (int)PLASTER_TEX), tscale);
	tp.emissive = (is_on ? 1.0 : 0.0);
	rgeom_mat_t &mat(mats_lights.get_material(tp, 0)); // no shadows
	if      (c.shape == SHAPE_CUBE  ) {mat.add_cube_to_verts  (c, c.color, c.get_llc(), EF_Z2);} // untextured, skip top face
	else if (c.shape == SHAPE_CYLIN ) {mat.add_vcylin_to_verts(c, c.color, 1, 0);} // bottom only
	else if (c.shape == SHAPE_SPHERE) {mat.add_sphere_to_verts(c, c.color);}
	else {assert(0);}
}

void building_room_geom_t::add_rug(room_object_t const &c) {
	bool const swap_tex_st(c.dy() < c.dx()); // rug textures are oriented with the long side in X, so swap the coordinates (rotate 90 degrees) if our rug is oriented the other way
	get_material(tid_nm_pair_t(c.get_rug_tid(), 0.0)).add_cube_to_verts(c, c.color, c.get_llc(), 61, swap_tex_st); // only draw top/+z face
}

void building_room_geom_t::add_picture(room_object_t const &c) { // also whiteboards
	bool const whiteboard(c.type == TYPE_WBOARD);
	int picture_tid(WHITE_TEX);

	if (!whiteboard && !(c.flags & RO_FLAG_TAKEN1)) { // picture, not taken/frame only
		int const user_tid(get_rand_screenshot_texture(c.obj_id));
		picture_tid  = ((user_tid >= 0) ? (unsigned)user_tid : c.get_picture_tid()); // if user texture is valid, use that instead
		num_pic_tids = get_num_screenshot_tids();
		has_pictures = 1;
	}
	unsigned skip_faces(get_face_mask(c.dim, c.dir)); // only the face oriented outward
	bool const mirror_x(!whiteboard && !(c.dim ^ c.dir));
	vector3d const tex_origin(c.get_llc());
	get_material(tid_nm_pair_t(picture_tid, 0.0)).add_cube_to_verts(c, c.color, tex_origin, skip_faces, !c.dim, mirror_x);
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
		mat.quad_verts.emplace_back(i->v, nc, i->t[0], i->t[1], cw);
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

void building_room_geom_t::add_bookcase(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale, bool no_shelves, float sides_scale, point const *const use_this_tex_origin) {
	colorRGBA const color(apply_wood_light_color(c));
	unsigned const skip_faces(~get_face_mask(c.dim, !c.dir)); // skip back face
	unsigned const skip_faces_shelves(skip_faces | get_skip_mask_for_xy(!c.dim)); // skip back face and sides
	float const width(c.get_sz_dim(!c.dim)), depth((c.dir ? -1.0 : 1.0)*c.get_sz_dim(c.dim)); // signed depth
	float const side_thickness(0.06*sides_scale*width);
	point const tex_origin(use_this_tex_origin ? *use_this_tex_origin : c.get_llc());
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
	c.set_rand_gen_state(rgen);
	unsigned const num_shelves(3 + ((17*c.room_id + int(1000.0*fabs(c.z1())))%3)); // 3-5, randomly selected by room ID and floor
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
			bool const backwards((rgen.rand()%10) == 0), book_dir(c.dir ^ backwards ^ 1); // spine facing out 90% of the time
			room_object_t obj(book, TYPE_BOOK, c.room_id, c.dim, book_dir, c.flags, c.light_amt, SHAPE_CUBE, book_color);
			obj.obj_id = c.obj_id + 123*i + 1367*n;
			add_book(obj, inc_lg, inc_sm, tilt_angle, skip_faces, backwards); // detailed book, no title if backwards
			pos += width;
			last_book_pos = pos;
			prev_tilted   = (tilt_angle != 0.0); // don't tilt two books in a row
		} // for n
	} // for i
}
void building_room_geom_t::expand_bookcase(room_object_t const &c) {
	// WRITE, if needed
}

void add_wine_rack_bottles(room_object_t const &c, vector<room_object_t> &objects) {
	float const height(c.dz()), width(c.get_sz_dim(!c.dim)), depth(c.get_sz_dim(c.dim)), shelf_thick(0.1*depth);
	unsigned const num_rows(max(1, round_fp(2.0*height/depth))), num_cols(max(1, round_fp(2.0*width/depth)));
	float const row_step((height - shelf_thick)/num_rows), col_step((width - shelf_thick)/num_cols);
	float const space_w(col_step - shelf_thick), space_h(row_step - shelf_thick), diameter(min(space_w, space_h) - shelf_thick);
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	room_object_t bottle(c);
	bottle.type   = TYPE_BOTTLE;
	bottle.dir   ^= 1;
	bottle.flags |= (RO_FLAG_WAS_EXP | RO_FLAG_NO_CONS); // wine bottles on wine racks are not consumable
	set_wall_width(bottle, (c.d[!c.dim][0] + 0.5f*(col_step + shelf_thick)), 0.5*diameter, !c.dim); // center in this dim

	for (unsigned i = 0; i < num_cols; ++i) { // columns/vertical
		bottle.z1() = c.z1() + shelf_thick; // rest on the top of the shelf
		bottle.z2() = bottle.z1() + diameter;

		for (unsigned j = 0; j < num_rows; ++j) { // rows/horizontal
			if (rgen.rand()%3) { // add a bottle 67% of the time; add_bottom=1
				bottle.set_as_bottle(3 + 64*rgen.rand_bool()); // always wine, but mixed cap color
				objects.push_back(bottle);
			}
			bottle.translate_dim(2, row_step); // translate in Z
		}
		bottle.translate_dim(!c.dim, col_step); // translate in !dim
	} // for i
}
void building_room_geom_t::add_wine_rack(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale) {
	if (inc_lg) { // add wooden frame
		float const height(c.dz()), width(c.get_sz_dim(!c.dim)), depth(c.get_sz_dim(c.dim)), shelf_thick(0.1*depth);
		unsigned const num_rows(max(1, round_fp(2.0*height/depth))), num_cols(max(1, round_fp(2.0*width/depth)));
		float const row_step((height - shelf_thick)/num_rows), col_step((width - shelf_thick)/num_cols);
		colorRGBA const color(apply_wood_light_color(c));
		rgeom_mat_t &wood_mat(get_wood_material(tscale));
		cube_t frame(c);
		frame.d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*0.09*c.get_sz_dim(c.dim); // slightly less depth so that bottles stick out a bit

		// create rows and columns of cubbies by intersecting horizontal and vertical cubes
		for (unsigned i = 0; i <= num_rows; ++i) { // rows/horizontal
			cube_t hc(frame);
			hc.z1() = c .z1() + row_step*i;
			hc.z2() = hc.z1() + shelf_thick;
			wood_mat.add_cube_to_verts(hc, color, tex_origin, 0); // draw all faces, even the back, in case it's visible through the window
		}
		for (unsigned i = 0; i <= num_cols; ++i) { // columns/vertical
			cube_t vc(frame);
			vc.d[!c.dim][0] = c .d[!c.dim][0] + col_step*i;
			vc.d[!c.dim][1] = vc.d[!c.dim][0] + shelf_thick;
			wood_mat.add_cube_to_verts(vc, color, tex_origin, 0); // draw all faces, even the back, in case it's visible through the window
		}
	}
	if (inc_sm && !(c.flags & RO_FLAG_EXPANDED)) { // add wine bottles if not expanded
		vector<room_object_t> &objects(get_temp_objects());
		add_wine_rack_bottles(c, objects);
		add_small_static_objs_to_verts(objects);
	}
}
void building_room_geom_t::expand_wine_rack(room_object_t const &c) {
	add_wine_rack_bottles(c, expanded_objs);
}

room_object_t get_desk_drawers_part(room_object_t const &c) {
	bool const side(c.obj_id & 1);
	float const desk_width(c.get_sz_dim(!c.dim)), height(c.dz());
	room_object_t drawers(c);
	drawers.z1() += 0.05*height; // shift up a bit from the floor
	drawers.z2()  = c.z1() + 0.85*height;
	drawers.expand_by_xy(-0.15*get_tc_leg_width(c, 0.06));
	drawers.d[!c.dim][!side] += (side ? 1.0 : -1.0)*0.75*desk_width; // put the drawers off to one side
	return drawers;
}
void building_room_geom_t::add_desk(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm) {
	// desk top and legs, similar to add_table()
	float const height(c.dz());
	cube_t top(c);
	top.z1() += 0.85*height;
	vector3d const tex_origin(c.get_llc());
	colorRGBA const color(apply_wood_light_color(c));

	if (inc_lg) {
		cube_t legs_bcube(c);
		legs_bcube.z2() = top.z1();
		get_wood_material(tscale).add_cube_to_verts(top, color, tex_origin); // all faces drawn
		add_tc_legs(legs_bcube, color, 0.06, tscale);
	}
	if (c.room_id & 3) { // add drawers 75% of the time
		room_object_t drawers(get_desk_drawers_part(c));
		if (inc_lg) {get_wood_material(tscale).add_cube_to_verts(drawers, color, tex_origin);} // all faces drawn

		if (inc_sm) {
			bool const side(c.obj_id & 1);
			drawers.d[!c.dim][side] -= (side ? 1.0 : -1.0)*0.85*get_tc_leg_width(c, 0.06); // make sure the drawers can pull out without hitting the desk legs
			add_dresser_drawers(drawers, tscale);
		}
	}
	if (inc_lg && c.shape == SHAPE_TALL) { // add top/back section of desk; this part is outside the bcube
		room_object_t c_top_back(c);
		set_cube_zvals(c_top_back, top.z2(), (top.z2() + 1.8*height));
		c_top_back.d[c.dim][c.dir] += 0.75*(c.dir ? -1.0 : 1.0)*c.get_sz_dim(c.dim);
		add_bookcase(c_top_back, 1, 1, tscale, 1, 0.4, &tex_origin); // no_shelves=1, side_width=0.4, both large and small, use same tex origin
	}
}

void building_room_geom_t::add_reception_desk(room_object_t const &c, float tscale) {
	vector3d const sz(c.get_size());
	float const top_z1(c.z1() + 0.94*sz.z), depth(sz[c.dim]), width(sz[!c.dim]), overhang(0.04*depth), lr_width(0.2*width), cutlen(depth - lr_width);
	assert(width > depth && cutlen > 0.0);
	colorRGBA const color(apply_light_color(c));
	// wood paneling sides
	rgeom_mat_t &side_mat(get_material(tid_nm_pair_t(PANELING_TEX, get_paneling_nm_tid(), 4.0*tscale, 4.0*tscale), 1)); // with shadows
	vector3d const tex_origin(c.get_llc());
	unsigned const lr_dim_mask(~get_face_mask(c.dim, c.dir));
	cube_t base(c);
	base.z2() = top_z1;
	base.expand_by_xy(-overhang);
	cube_t front(base), left(base), right(base);
	front.d[ c.dim][!c.dir] -= (c.dir ? -1.0 : 1.0)*cutlen;
	left .d[!c.dim][1] -= (width - lr_width);
	right.d[!c.dim][0] += (width - lr_width);
	left .d[ c.dim][c.dir] = right.d[ c.dim][c.dir] = front.d[ c.dim][!c.dir];
	side_mat.add_cube_to_verts(front, color, tex_origin, EF_Z2);
	side_mat.add_cube_to_verts(left,  color, tex_origin, (EF_Z2 | lr_dim_mask)); // skip top face
	side_mat.add_cube_to_verts(right, color, tex_origin, (EF_Z2 | lr_dim_mask)); // skip top face
	// shiny marble top
	// Note: I wanted to add cylinders to the left and right top to round the corners like in the mapx lobby, but it's not easy to get the textures to line up here
	tid_nm_pair_t top_tex(get_counter_tid(), 2.5*tscale);
	top_tex.set_specular(0.5, 80.0);
	rgeom_mat_t &top_mat(get_material(top_tex, 1)); // with shadows
	cube_t top_front(front), top_left(left), top_right(right);
	top_front.z1() = top_left.z1() = top_right.z1() = top_z1;
	top_front.z2() = top_left.z2() = top_right.z2() = c.z2();
	top_front.expand_by_xy(overhang);
	top_left .expand_by_xy(overhang);
	top_right.expand_by_xy(overhang);
	top_left.d[c.dim][c.dir] = top_right.d[c.dim][c.dir] = top_front.d[c.dim][!c.dir]; // remove the overlap
	top_mat.add_cube_to_verts(top_front, color, tex_origin, 0); // all faces drawn
	top_mat.add_cube_to_verts(top_left,  color, tex_origin, lr_dim_mask);
	top_mat.add_cube_to_verts(top_right, color, tex_origin, lr_dim_mask);
}

void add_pillow(cube_t const &c, rgeom_mat_t &mat, colorRGBA const &color, vector3d const &tex_origin) {
	unsigned const ndiv = 24; // number of quads in X and Y
	float const ndiv_inv(1.0/ndiv), dx_inv(1.0/c.dx()), dy_inv(1.0/c.dy());
	color_wrapper cw(color);
	auto &verts(mat.itri_verts); // Note: could cache verts
	unsigned const start(verts.size()), stride(ndiv + 1);
	float dists[ndiv+1] = {};
	norm_comp const nc(plus_z);

	for (unsigned x = 0; x <= ndiv; ++x) {
		float const v(2.0f*x*ndiv_inv - 1.0f); // centered on 0 in range [-1, 1]
		dists[x] = 0.5*SIGN(v)*sqrt(abs(v)) + 0.5; // nonlinear spacing, closer near the edges, convert back to [0, 1] range
	}
	for (unsigned y = 0; y <= ndiv; ++y) {
		float const yval(c.y1() + dists[y]*c.dy()), ey(2.0f*max(0.0f, min((yval - c.y1()), (c.y2() - yval)))*dy_inv);

		for (unsigned x = 0; x <= ndiv; ++x) {
			float const xval(c.x1() + dists[x]*c.dx()), ex(2.0f*max(0.0f, min((xval - c.x1()), (c.x2() - xval)))*dx_inv), zval(c.z1() + c.dz()*pow(ex*ey, 0.2f));
			verts.emplace_back(point(xval, yval, zval), nc, mat.tex.tscale_x*(xval - tex_origin.x), mat.tex.tscale_y*(yval - tex_origin.y), cw);
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

void get_bed_cubes(room_object_t const &c, cube_t cubes[6]) {
	float const height(c.dz()), length(c.get_sz_dim(c.dim)), width(c.get_sz_dim(!c.dim));
	bool const is_wide(width > 0.7*length), add_posts(is_wide && (c.obj_id & 1)); // no posts for twin beds
	float const head_width(0.04), foot_width(add_posts ? head_width : 0.03f); // relative to length
	cube_t frame(c), head(c), foot(c), mattress(c), legs_bcube(c), pillow(c);
	head.d[c.dim][!c.dir] += (c.dir ? 1.0 : -1.0)*(1.0 - head_width)*length;
	foot.d[c.dim][ c.dir] -= (c.dir ? 1.0 : -1.0)*(1.0 - foot_width)*length;
	mattress.d[c.dim][ c.dir] = head.d[c.dim][!c.dir];
	mattress.d[c.dim][!c.dir] = foot.d[c.dim][ c.dir];
	frame.z1() += 0.3*height;
	frame.z2() -= 0.65*height;
	foot.z2()  -= 0.2*height;
	mattress.z1()   = head.z1()   = foot.z1() = frame.z2();
	mattress.z2()   = pillow.z1() = mattress.z1() + 0.2*height;
	pillow.z2()     = pillow.z1() + 0.13*height;
	legs_bcube.z2() = frame.z1();
	mattress.expand_in_dim(!c.dim, -0.02*width); // shrink slightly
	float const pillow_space((is_wide ? 0.08 : 0.23)*width);
	pillow.expand_in_dim(!c.dim, -pillow_space);
	pillow.d[c.dim][ c.dir] = mattress.d[c.dim][ c.dir] + (c.dir ? -1.0 : 1.0)*0.02*length; // head
	pillow.d[c.dim][!c.dir] = pillow  .d[c.dim][ c.dir] + (c.dir ? -1.0 : 1.0)*(is_wide ? 0.25 : 0.6)*pillow.get_sz_dim(!c.dim);
	cubes[0] = frame; cubes[1] = head; cubes[2] = foot; cubes[3] = mattress; cubes[4] = pillow; cubes[5] = legs_bcube;
}
void building_room_geom_t::add_bed(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale) {
	float const height(c.dz()), length(c.get_sz_dim(c.dim)), width(c.get_sz_dim(!c.dim));
	bool const is_wide(width > 0.7*length), add_posts(is_wide && (c.obj_id & 1)), add_canopy(add_posts && (c.obj_id & 2)); // no posts for twin beds
	float const head_width(0.04), foot_width(add_posts ? head_width : 0.03f); // relative to length
	cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
	get_bed_cubes(c, cubes);
	cube_t const &frame(cubes[0]), &head(cubes[1]), &foot(cubes[2]), &mattress(cubes[3]), &pillow(cubes[4]), &legs_bcube(cubes[5]);
	colorRGBA const sheet_color(apply_light_color(c));
	tid_nm_pair_t const sheet_tex(c.get_sheet_tid(), tscale);

	if (inc_lg) {
		bool const no_mattress(c.flags & RO_FLAG_TAKEN3);
		colorRGBA const color(apply_wood_light_color(c));
		add_tc_legs(legs_bcube, color, max(head_width, foot_width), tscale);
		if (no_mattress) {get_wood_material(4.0*tscale);} // pre-allocate slats material if needed
		rgeom_mat_t &wood_mat(get_wood_material(tscale));
		vector3d const tex_origin(c.get_llc());

		if (no_mattress) { // mattress is gone, draw the slats on the bottom of the bed
			unsigned const num_slats = 12;
			unsigned const slat_skip_faces(get_skip_mask_for_xy(!c.dim));
			float const side_width(0.08*width), slat_spacing(length/num_slats), slat_width(0.45*slat_spacing), slat_gap(0.5*(slat_spacing - slat_width));
			cube_t sides[2] = {frame, frame}, slat(frame);
			sides[0].d[!c.dim][1] -= (width - side_width);
			sides[1].d[!c.dim][0] += (width - side_width);
			for (unsigned d = 0; d < 2; ++d) {wood_mat.add_cube_to_verts(sides[d], color, tex_origin);}
			slat.expand_in_dim(!c.dim, -side_width); // flush with sides
			cube_t ends[2] = {slat, slat};
			ends[0].d[c.dim][1] = frame.d[c.dim][0] + slat_gap;
			ends[1].d[c.dim][0] = frame.d[c.dim][1] - slat_gap;
			for (unsigned d = 0; d < 2; ++d) {wood_mat.add_cube_to_verts(ends[d], color, tex_origin, slat_skip_faces);}
			slat.d[c.dim][1] = slat.d[c.dim][0] + slat_spacing;
			slat.expand_in_dim(c.dim, -slat_gap); // add gap between slats to each side
			slat.expand_in_dim(2, -0.25*frame.dz()); // make thinner in Z
			rgeom_mat_t &slat_mat(get_wood_material(4.0*tscale));
			colorRGBA const slat_color(color*1.5); // make them lighter in color

			for (unsigned n = 0; n < num_slats; ++n) {
				slat_mat.add_cube_to_verts(slat, slat_color, tex_origin, (slat_skip_faces | EF_Z1));
				slat.translate_dim(c.dim, slat_spacing);
			}
		}
		else {
			wood_mat.add_cube_to_verts(frame, color, tex_origin);
		}
		wood_mat.add_cube_to_verts(head, color, tex_origin, EF_Z1);
		wood_mat.add_cube_to_verts(foot, color, tex_origin, EF_Z1);
		
		if (add_posts) { // maybe add bed posts and canopy; these extend outside of the bed bcube, but that probably doesn't matter
			float const post_width(min(head_width, foot_width)*length);
			cube_t posts_area(c);
			posts_area.z1() = foot.z2(); // start at the foot
			posts_area.z2() = posts_area.z1() + (add_canopy ? 1.4 : 0.6)*height; // higher posts for canopy bed
			cube_t posts[4];
			get_tc_leg_cubes_abs_width(posts_area, post_width, posts);
			bool const use_cylinders(!add_canopy && (c.obj_id & 4));

			for (unsigned i = 0; i < 4; ++i) {
				if (!add_canopy && posts[i].d[c.dim][!c.dir] == c.d[c.dim][!c.dir]) {posts[i].translate_dim(2, -(head.z2() - foot.z2()));} // make footboard posts shorter
				if (use_cylinders) {wood_mat.add_vcylin_to_verts(posts[i], color, 0, 1);}
				else {wood_mat.add_cube_to_verts(posts[i], color, tex_origin, EF_Z1);} // skip bottom face
			}
			if (add_canopy) {
				for (unsigned i = 0; i < 4; ++i) { // add 4 horizontal cube bars along the top of the bed connecting each adjacent pair of posts
					cube_t top(posts[i]);
					unsigned const next_ix[4] = {1, 3, 0, 2};
					top.union_with_cube(posts[next_ix[i]]); // next post
					top.z1() = top.z2() - post_width; // height = width
					bool const dim(top.dx() < top.dy());
					top.expand_in_dim(dim, -post_width); // remove overlaps with the post
					wood_mat.add_cube_to_verts(top, color, tex_origin, get_skip_mask_for_xy(dim));
				}
				// TODO: add material to the top?
			}
		}
		if (!no_mattress) {
			unsigned const mattress_skip_faces(EF_Z1 | get_skip_mask_for_xy(c.dim));
			if (c.flags & RO_FLAG_TAKEN2) {get_material(untex_shad_mat, 1).add_cube_to_verts(mattress, sheet_color, tex_origin, mattress_skip_faces);} // sheets taken, bare mattress
			else {get_material(sheet_tex, 1).add_cube_to_verts(mattress, sheet_color, tex_origin, mattress_skip_faces);} // draw matterss with sheets
		}
	}
	if (inc_sm && !(c.flags & RO_FLAG_TAKEN1)) { // draw pillows if not taken
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

	if (c.shape == SHAPE_CYLIN) {
		mat.add_vcylin_to_verts(c, color, 1, 0, 1, 1, 0.7, 1.0); // untextured, bottom only, two_sided cylinder with inverted bottom normal
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

	if (c.shape == SHAPE_SHORT) { // wall separating urinals, drawn as a single cube
		mat.add_cube_to_verts(c, color, tex_origin, ~get_face_mask(c.dim, c.dir));
		return;
	}
	float const dz(c.dz()), wall_thick(0.0125*dz), frame_thick(2.0*wall_thick), door_gap(0.3*wall_thick);
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
	unsigned const side_skip_mask(get_skip_mask_for_xy(c.dim));
	mat.add_cube_to_verts(side1,  color, tex_origin, side_skip_mask);
	mat.add_cube_to_verts(side2,  color, tex_origin, side_skip_mask);
	mat.add_cube_to_verts(front1, color, tex_origin, EF_Z12);
	mat.add_cube_to_verts(front2, color, tex_origin, EF_Z12);

	if (c.is_open()) {
		// draw open door?
	}
	else {
		door.expand_in_dim(!c.dim, -door_gap);
		mat.add_cube_to_verts(door, color, tex_origin);
	}
}

int get_cubicle_tid(room_object_t const &c) {return get_texture_by_name((c.obj_id & 1) ? "carpet/carpet1.jpg" : "carpet/carpet2.jpg");} // select from one of 2 textures

void building_room_geom_t::add_cubicle(room_object_t const &c, float tscale) {
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_cubicle_tid(c), tscale), 1));
	colorRGBA const color(apply_light_color(c));
	point const tex_origin(c.get_llc());
	float const dz(c.dz()), wall_thick(0.07*dz), frame_thick(8.0*wall_thick), dir_sign(c.dir ? 1.0 : -1.0);
	bool const is_short(c.shape == SHAPE_SHORT);
	cube_t sides(c), front(c), back(c);
	if (is_short) {back.z2() -= 0.4*dz;}
	sides.d[c.dim][!c.dir] += dir_sign*wall_thick; // front
	sides.d[c.dim][ c.dir] -= dir_sign*wall_thick; // back
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
	rgeom_mat_t &edge_mat(get_material(tid_nm_pair_t(), 0)); // unshadowed
	unsigned const side_edge_skip_mask (~(EF_Z2 | (is_short ? ~get_face_mask(c.dim, c.dir) : 0)));
	unsigned const front_edge_skip_mask(~(EF_Z2 | get_skip_mask_for_xy(!c.dim)));
	colorRGBA const edge_color(apply_light_color(c, BKGRAY));
	edge_mat.add_cube_to_verts(side1,  edge_color, tex_origin, side_edge_skip_mask);
	edge_mat.add_cube_to_verts(side2,  edge_color, tex_origin, side_edge_skip_mask);
	edge_mat.add_cube_to_verts(front1, edge_color, tex_origin, front_edge_skip_mask);
	edge_mat.add_cube_to_verts(front2, edge_color, tex_origin, front_edge_skip_mask);
	edge_mat.add_cube_to_verts(back,   edge_color, tex_origin, ~EF_Z2);
	// desk surface
	rgeom_mat_t &surf_mat(get_material(tid_nm_pair_t(MARBLE_TEX, 4.0*tscale), 1));
	colorRGBA const surf_color(apply_light_color(c, LT_GRAY));
	cube_t surface(sides);
	set_cube_zvals(surface, (c.z1() + 0.45*dz), (c.z1() + 0.50*dz));
	cube_t surf1(surface), surf2(surface), surf3(surface); // left, right, back
	surf1.d[!c.dim][0] = side1.d[!c.dim][1];
	surf1.d[!c.dim][1] = surf3.d[!c.dim][0] = front1.d[!c.dim][1];
	surf2.d[!c.dim][0] = surf3.d[!c.dim][1] = front2.d[!c.dim][0];
	surf2.d[!c.dim][1] = side2.d[!c.dim][0];
	surf3.d[ c.dim][!c.dir] = surface.d[c.dim][c.dir] - dir_sign*frame_thick;
	surf_mat.add_cube_to_verts(surf1, surf_color, tex_origin, get_skip_mask_for_xy( c.dim));
	surf_mat.add_cube_to_verts(surf2, surf_color, tex_origin, get_skip_mask_for_xy( c.dim));
	surf_mat.add_cube_to_verts(surf3, surf_color, tex_origin, get_skip_mask_for_xy(!c.dim));
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
		unsigned const skip_faces((c.flags & RO_FLAG_HANGING) ? 0 : ~get_face_mask(c.dim, !c.dir)); // skip back face, unless hanging
		get_material(tid_nm_pair_t(), 0).add_cube_to_verts(c, WHITE, zero_vector, skip_faces); // back of the sign, always white (for now)
	}
	if (!inc_text) return;
	// add sign text
	cube_t ct(c); // text area is slightly smaller than full cube
	ct.expand_in_dim(!c.dim, -0.1*c.get_sz_dim(!c.dim));
	ct.expand_in_dim(2, -0.05*c.dz());
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
	if (c.flags & RO_FLAG_EMISSIVE) {tex.emissive = 1.0;}
	rgeom_mat_t &mat(get_material(tex, 0, 0, 1));
	color_wrapper const cw(apply_light_color(c)); // set alpha=1.0
	norm_comp const nc(normal);

	for (auto i = verts.begin(); i != verts.end(); ++i) {
		i->v[!c.dim] = i->v[!c.dim]*width_scale + ct.d[!c.dim][!ldir]; // line
		i->v.z       = i->v.z*height_scale + ct.z1(); // column
		mat.quad_verts.emplace_back(i->v, nc, i->t[0], i->t[1], cw);
	}
}

void building_room_geom_t::add_counter(room_object_t const &c, float tscale) { // for kitchens
	float const dz(c.dz()), depth(c.get_sz_dim(c.dim)), dir_sign(c.dir ? 1.0 : -1.0);
	cube_t top(c), cabinet_gap;
	top.z1() += 0.95*dz;
	tid_nm_pair_t const marble_tex(get_counter_tid(), 2.5*tscale);
	rgeom_mat_t &top_mat(get_material(marble_tex, 1));
	colorRGBA const top_color(apply_light_color(c, WHITE));

	if (c.type == TYPE_KSINK || c.type == TYPE_BRSINK) { // counter with kitchen or bathroom sink
		float const sdepth(0.8*depth), width(c.get_sz_dim(!c.dim)), swidth(min(1.4f*sdepth, 0.75f*width));
		vector3d const center(c.get_cube_center());
		vector3d faucet_pos(center);
		faucet_pos[c.dim] -= dir_sign*0.56*sdepth;
		cube_t sink(center, center), faucet1(faucet_pos, faucet_pos);
		set_cube_zvals(sink,   (top.z1() - 0.25*dz), top.z1());
		set_cube_zvals(faucet1, top.z2(), (top.z2() + 0.30*dz));
		sink.expand_in_dim( c.dim, 0.5*sdepth);
		sink.expand_in_dim(!c.dim, 0.5*swidth);
		faucet1.expand_in_dim( c.dim, 0.04*sdepth);
		faucet1.expand_in_dim(!c.dim, 0.05*swidth);
		cube_t faucet2(faucet1);
		faucet2.z1()  = faucet1.z2();
		faucet2.z2() += 0.035*dz;
		faucet2.d[c.dim][c.dir] += dir_sign*0.28*sdepth;
		vect_cube_t &cubes(get_temp_cubes());
		subtract_cube_from_cube(top, sink, cubes);
		for (auto i = cubes.begin(); i != cubes.end(); ++i) {top_mat.add_cube_to_verts(*i, top_color, tex_origin);} // should always be 4 cubes
		colorRGBA const sink_color(apply_light_color(c, GRAY));
		get_metal_material(0).add_cube_to_verts(sink, sink_color, tex_origin, EF_Z2, 0, 0, 0, 1); // basin: inverted, skip top face, unshadowed
		rgeom_mat_t &metal_mat(get_metal_material(1)); // shadowed, specular metal (specular doesn't do much because it's flat, but may make more of a diff using a cylinder later)
		metal_mat.add_cube_to_verts(faucet1, sink_color, tex_origin, EF_Z12); // vertical part of faucet, skip top and bottom faces
		metal_mat.add_cube_to_verts(faucet2, sink_color, tex_origin, 0); // horizontal part of faucet, draw all faces

		if (c.type == TYPE_BRSINK) { // bathroom sink
			metal_mat.add_cube_to_verts(sink, sink_color, tex_origin, EF_Z2); // outside of basin, no top surface, shadowed
			cube_t front(c);
			front.z2() = top.z1();
			front.z1() = sink.z1() - 0.1*dz; // slightly below the sink basin
			front.d[c.dim][!c.dir] += dir_sign*0.94*depth;
			get_material(marble_tex, 1).add_cube_to_verts(front, top_color, tex_origin, EF_Z2); // front surface, no top face; same as top_mat
		}
		if (c.type == TYPE_KSINK && width > 3.5*depth) { // kitchen sink - add dishwasher
			bool const side((c.flags & RO_FLAG_ADJ_LO) ? 1 : ((c.flags & RO_FLAG_ADJ_HI) ? 0 : (c.obj_id & 1))); // left/right of the sink
			unsigned const dw_skip_faces(~get_face_mask(c.dim, !c.dir));
			colorRGBA const dw_color(apply_light_color(c, LT_GRAY));
			cube_t dishwasher(c);
			dishwasher.z1() += 0.05*dz;
			dishwasher.z2() -= 0.05*dz;
			dishwasher.d[ c.dim][!c.dir]  = c.d[c.dim][c.dir] - dir_sign*0.1*depth;
			dishwasher.d[ c.dim][ c.dir] += dir_sign*0.05*depth; // front
			dishwasher.d[!c.dim][!side ]  = sink.d[!c.dim][side] + (side ? 1.0 : -1.0)*0.1*depth;
			dishwasher.d[!c.dim][ side ]  = dishwasher.d[!c.dim][!side] + (side ? 1.0 : -1.0)*1.05*depth;
			metal_mat.add_cube_to_verts(dishwasher, dw_color, tex_origin, dw_skip_faces);
			cube_t dishwasher_back(dishwasher), handle(dishwasher);
			dishwasher_back.d[c.dim][!c.dir] = c.d[c.dim][!c.dir]; // flush with the cabinet
			metal_mat.add_cube_to_verts(dishwasher_back, dw_color, tex_origin, ~dw_skip_faces); // draw only the back face, in case it's visible through a window
			handle.z1() += 0.77*dz;
			handle.z2() -= 0.10*dz;
			handle.expand_in_dim(!c.dim, -0.1*depth);
			handle.d[c.dim][ c.dir]  = handle.d[c.dim][!c.dir] = dishwasher.d[c.dim][c.dir];
			handle.d[c.dim][ c.dir] += dir_sign*0.04*depth; // front
			metal_mat.add_ortho_cylin_to_verts(handle, sink_color, !c.dim, 1, 1); // add handle as a cylinder in the proper dim with both ends
			cabinet_gap = dishwasher;
		}
	}
	else { // regular counter top
		top_mat.add_cube_to_verts(top, top_color, tex_origin); // top surface, all faces
	}
	if (c.type != TYPE_BRSINK) { // add wood sides of counter/cabinet
		float const overhang(0.05*depth);
		room_object_t cabinet(c);
		cabinet.z2() = top.z1();
		//cabinet.expand_in_dim(!c.dim, -overhang); // add side overhang: disable to allow cabinets to be flush with objects
		cabinet.d[c.dim][c.dir] -= dir_sign*overhang; // add front overhang

		if (!cabinet_gap.is_all_zeros()) { // split cabinet into two parts to avoid the dishwasher
			room_object_t left_part(cabinet);
			left_part.d[!c.dim][1] = cabinet_gap.d[!c.dim][0];
			cabinet  .d[!c.dim][0] = cabinet_gap.d[!c.dim][1];
			left_part.flags &= ~RO_FLAG_ADJ_HI;
			cabinet  .flags &= ~RO_FLAG_ADJ_LO;
			add_cabinet(left_part, tscale);
		}
		add_cabinet(cabinet, tscale); // draw the wood part
	}
	if (c.item_flags) { // add backsplash, 50% chance of tile vs. matching marble
		tid_nm_pair_t const bs_tex((c.room_id & 1) ? marble_tex : tid_nm_pair_t(get_texture_by_name("bathroom_tile.jpg"), 2.5*tscale));
		rgeom_mat_t &bs_mat(get_material(bs_tex, 0)); // no shadows
		cube_t bsz(c);
		bsz.z1()  = c.z2();
		bsz.z2() += 0.33*c.dz();

		if (c.item_flags & 1) { // back
			cube_t bs(bsz);
			bs.d[c.dim][c.dir] -= (c.dir ? 1.0 : -1.0)*0.99*depth;
			bs_mat.add_cube_to_verts(bs, top_color, zero_vector, (EF_Z1 | ~get_face_mask(c.dim, !c.dir)));
		}
		for (unsigned d = 0; d < 2; ++d) { // handle the other dim
			if (!(c.item_flags & (1<<(d+1)))) continue; // not adjacent in this dir
			cube_t bs(bsz);
			bs.d[!c.dim][!d] -= (d ? -1.0 : 1.0)*(c.get_sz_dim(!c.dim) - 0.01*depth);
			bs_mat.add_cube_to_verts(bs, top_color, zero_vector, (EF_Z1 | ~get_face_mask(!c.dim, d)));
		}
	}
}

void building_room_geom_t::add_cabinet(room_object_t const &c, float tscale) { // for kitchens
	assert(c.is_strictly_normalized());
	unsigned const skip_faces((c.type == TYPE_COUNTER) ? EF_Z12 : EF_Z2); // skip top face (can't skip back in case it's against a window)
	get_wood_material(tscale).add_cube_to_verts(c, apply_wood_light_color(c), tex_origin, skip_faces);
	// add cabinet doors
	float const cab_depth(c.get_sz_dim(c.dim)), door_height(0.8*c.dz()), door_width(0.75*door_height), door_thick(0.05*door_height);
	cube_t front(c);
	if (c.flags & RO_FLAG_ADJ_LO) {front.d[!c.dim][0] += cab_depth;} // exclude L-joins of cabinets from having doors; assumes all cabinets are the same depth
	if (c.flags & RO_FLAG_ADJ_HI) {front.d[!c.dim][1] -= cab_depth;}
	float const handle_thick(0.75*door_thick), cab_width(front.get_sz_dim(!c.dim));
	if (cab_width < 0.0) return; // this seems to happen on occasion; maybe it's a bug, or maybe the random size parameters can lead to bad values; either way, skip it
	float door_spacing(1.2*door_width);
	unsigned const num_doors(floor(cab_width/door_spacing));
	if (num_doors == 0) return; // is this possible?
	assert(num_doors < 1000); // sanity check
	door_spacing = cab_width/num_doors;
	float const tb_border(0.5f*(c.dz() - door_height)), side_border(0.16*door_width), dir_sign(c.dir ? 1.0 : -1.0);
	float lo(front.d[!c.dim][0]);
	get_metal_material(0); // ensure material exists so that door_mat reference is not invalidated
	rgeom_mat_t &door_mat(get_material(get_tex_auto_nm(WOOD2_TEX, 2.0*tscale), 0)); // unshadowed
	rgeom_mat_t &handle_mat(get_metal_material(0)); // untextured, unshadowed
	colorRGBA const door_color(apply_light_color(c, WHITE)); // lighter color than cabinet
	colorRGBA const handle_color(apply_light_color(c, GRAY_BLACK));
	unsigned const door_skip_faces(~get_face_mask(c.dim, !c.dir));
	cube_t door(c);
	door.d[ c.dim][!c.dir]  = door.d[c.dim][c.dir];
	door.d[ c.dim][ c.dir] += dir_sign*door_thick; // expand out a bit
	door.expand_in_dim(2, -tb_border  ); // shrink in Z
	cube_t handle(door);
	handle.d[ c.dim][!c.dir]  = door.d[c.dim][c.dir];
	handle.d[ c.dim][ c.dir] += dir_sign*handle_thick; // expand out a bit
	handle.expand_in_dim(2, -0.4*door.dz()); // shrink in Z
	
	for (unsigned n = 0; n < num_doors; ++n) {
		float const hi(lo + door_spacing);
		door.d[!c.dim][0] = lo;
		door.d[!c.dim][1] = hi;
		door.expand_in_dim(!c.dim, -side_border); // shrink in XY
		door_mat.add_cube_to_verts(door, door_color, tex_origin, door_skip_faces);
		lo = hi;
		// add door handle
		float const dwidth(door.get_sz_dim(!c.dim)), hwidth(0.04*door.dz()); // the actual door and handle widths
		float const near_side(0.1*dwidth), far_side(dwidth - near_side - hwidth);
		bool const side(n & 1); // alternate
		handle.d[!c.dim][0] = door.d[!c.dim][0] + (side ? near_side : far_side);
		handle.d[!c.dim][1] = door.d[!c.dim][1] - (side ? far_side : near_side);
		handle_mat.add_cube_to_verts(handle, handle_color, tex_origin, door_skip_faces); // same skip_faces
	} // for n
}

void building_room_geom_t::add_window(room_object_t const &c, float tscale) { // frosted window blocks
	// Maybe windows should be refractive + blurred to simulate frosted glass?
	// - Using a separate drawing pass like reflections could be slow because there can be multiple windows in one bathroom,
	//   and they can be seen when the player is outside the bathroom.
	// - What about post-processing blur?
	//   - Can't use the stencil buffer because that's used (and cleared) by the main drawing code after drawing windows.
	//   - Drawing windows last won't properly alpha blend with other windows, showers, water, etc., and the depth buffer may be wrong.
	//   - Drawing windows in world space on the CPU (like blast effects) doesn't correctly handle occlusion by fragments of lower depth (such as interior walls).
	unsigned const skip_faces(get_skip_mask_for_xy(!c.dim) | EF_Z12); // only enable faces in dim
	cube_t window(c);
	tid_nm_pair_t tex(get_bath_wind_tid(), 0.0); // fit texture to the window front/back faces
	if (c.is_lit()) {tex.emissive = 0.33;} // one third emissive
	get_material(tex, 0).add_cube_to_verts(window, c.color, c.get_llc(), skip_faces); // no apply_light_color()
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
	tid_nm_pair_t tex(((c.shape == SHAPE_SHORT) ? c.get_comp_monitor_tid() : c.get_picture_tid()), 0.0); // computer monitor vs. TV
	tex.emissive = 1.0;
	get_material(tex).add_cube_to_verts(screen, WHITE, c.get_llc(), skip_faces, !c.dim, !(c.dim ^ c.dir));
}

void building_room_geom_t::add_potted_plant(room_object_t const &c, bool inc_pot, bool inc_plant) {
	float const plant_diameter(2.0*c.get_radius()), stem_radius(0.04*plant_diameter);
	float const pot_height(max(0.6*plant_diameter, 0.3*c.dz())), pot_top(c.z1() + pot_height), dirt_level(pot_top - 0.15*pot_height);
	float const cx(c.get_center_dim(0)), cy(c.get_center_dim(1));
	point const base_pos(cx, cy, dirt_level); // base of plant trunk, center of dirt disk

	if (inc_pot) {
		// draw the pot, tapered with narrower bottom; draw the bottom of the pot if there's no dirt
		bool const no_dirt(c.flags & RO_FLAG_TAKEN2);
		float const pot_radius(0.4*plant_diameter);
		get_material(untex_shad_mat, 1).add_cylin_to_verts(point(cx, cy, c.z1()), point(cx, cy, pot_top), 0.65*pot_radius, pot_radius, apply_light_color(c), no_dirt, 0, 1, 0);
		
		if (!no_dirt) { // draw dirt in the pot as a disk if not taken
			rgeom_mat_t &dirt_mat(get_material(tid_nm_pair_t(get_texture_by_name("rock2.png")), 1)); // use dirt texture
			dirt_mat.add_disk_to_verts(base_pos, 0.947*pot_radius, 0, apply_light_color(c, WHITE));
		}
	}
	if (inc_plant && !(c.flags & RO_FLAG_TAKEN1)) { // plant not taken
		// draw plant leaves
		s_plant plant;
		plant.create_no_verts(base_pos, (c.z2() - base_pos.z), stem_radius, c.obj_id, 0, 1); // land_plants_only=1
		static vector<vert_norm_comp> points;
		points.clear();
		plant.create_leaf_points(points, 10.0, 1.5, 4); // plant_scale=10.0 seems to work well; more levels and rings
		auto &leaf_verts(mats_plants.get_material(tid_nm_pair_t(plant.get_leaf_tid()), 1).quad_verts);
		color_wrapper const leaf_cw(apply_light_color(c, plant.get_leaf_color()));
		float const ts[4] = {0,1,1,0}, tt[4] = {0,0,1,1};
		for (unsigned i = 0; i < points.size(); ++i) {leaf_verts.emplace_back(vert_norm_comp_tc(points[i], ts[i&3], tt[i&3]), leaf_cw);}
		// draw plant stem
		colorRGBA const stem_color(plant.get_stem_color());
		mats_plants.get_material(get_tex_auto_nm(WOOD2_TEX), 1).add_cylin_to_verts(point(cx, cy, base_pos.z), point(cx, cy, c.z2()), stem_radius, 0.0f, stem_color, 0, 0); // stem
	}
}

int get_lg_ball_tid   (room_object_t const &c) {return get_texture_by_name((c.item_flags & 1) ? "interiors/basketball.png" : "interiors/soccer_ball_diffuse.png");}
int get_lg_ball_nm_tid(room_object_t const &c) {return ((c.item_flags & 1) ? -1 : get_texture_by_name("interiors/soccer_ball_normal.png"));}

void building_room_geom_t::add_lg_ball(room_object_t const &c) { // is_small=1
	bool const dynamic(c.is_dynamic()); // either small or dynamic
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_lg_ball_tid(c), get_lg_ball_nm_tid(c), 0.0, 0.0), 1, dynamic, !dynamic));
	// rotate the texture coords when the ball is rolling
	mat.add_sphere_to_verts(c, apply_light_color(c), 0, (c.has_dstate() ? &get_dstate(c).rot_matrix : nullptr)); // low_detail=0
}
/*static*/ void building_room_geom_t::draw_lg_ball_in_building(room_object_t const &c, shader_t &s) {
	float const angle(atan2(cview_dir.y, cview_dir.x)); // angle of camera view in XY plane, for rotating about Z
	xform_matrix rot_matrix(get_rotation_matrix(plus_z, angle));
	rgeom_mat_t mat(tid_nm_pair_t(get_lg_ball_tid(c), get_lg_ball_nm_tid(c), 0.0, 0.0));
	mat.add_sphere_to_verts(c, apply_light_color(c), 0, &rot_matrix); // low_detail=0
	mat.tex.set_gl(s);
	draw_quad_verts_as_tris(mat.quad_verts);
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
rgeom_mat_t &building_room_geom_t::get_wood_material(float tscale, bool inc_shadows, bool dynamic, bool small) {
	return get_material(get_tex_auto_nm(WOOD2_TEX, tscale), inc_shadows, dynamic, small); // hard-coded for common material
}
rgeom_mat_t &building_room_geom_t::get_metal_material(bool inc_shadows, bool dynamic, bool small) {
	tid_nm_pair_t tex;
	tex.set_specular(0.8, 60.0);
	return get_material(tex, inc_shadows, dynamic, small);
}
colorRGBA get_textured_wood_color() {return WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX));} // Note: uses default WOOD_COLOR, not the per-building random variant
colorRGBA get_counter_color      () {return (get_textured_wood_color()*0.75 + texture_color(get_counter_tid())*0.25);}

// Note: we could add colorRGBA(0.8, 0.9, 1.0, 0.4) for water bottles, but transparent objects require removing interior faces such as half of the sphere
colorRGBA const bottle_colors[NUM_BOTTLE_COLORS] = {colorRGBA(0.4, 0.7, 1.0), colorRGBA(0.2, 0.1, 0.05), colorRGBA(0.1, 0.4, 0.1), BLACK}; // water, Coke, green beer bottle, wine

void room_object_t::set_as_bottle(unsigned rand_id) {
	obj_id = rand_id;
	color  = bottle_colors[obj_id%NUM_BOTTLE_COLORS];
}

colorRGBA room_object_t::get_color() const {
	switch (type) {
	case TYPE_TABLE:    return get_textured_wood_color();
	case TYPE_CHAIR:    return (color + get_textured_wood_color())*0.5; // 50% seat color / 50% wood legs color
	case TYPE_STAIR:    return (STAIRS_COLOR_TOP*0.5 + STAIRS_COLOR_BOT*0.5).modulate_with(texture_color(MARBLE_TEX));
	case TYPE_STAIR_WALL: return texture_color(STUCCO_TEX);
	case TYPE_ELEVATOR: return LT_BROWN; // ???
	case TYPE_RUG:      return texture_color(get_rug_tid());
	case TYPE_PICTURE:  return texture_color(get_picture_tid());
	case TYPE_BCASE:    return get_textured_wood_color();
	case TYPE_WINE_RACK:return get_textured_wood_color();
	case TYPE_DESK:     return get_textured_wood_color();
	case TYPE_RDESK:    return (texture_color(PANELING_TEX)*0.5 + texture_color(get_counter_tid())*0.5);
	case TYPE_BED:      return (color.modulate_with(texture_color(get_sheet_tid())) + get_textured_wood_color())*0.5; // half wood and half cloth
	case TYPE_COUNTER:  return get_counter_color();
	case TYPE_KSINK:    return (get_counter_color()*0.9 + GRAY*0.1); // counter, with a bit of gray mixed in from the sink
	case TYPE_BRSINK:   return texture_color(get_counter_tid()).modulate_with(color);
	case TYPE_CABINET:  return get_textured_wood_color();
	case TYPE_PLANT:    return (color*0.75 + blend_color(GREEN, BROWN, 0.5, 0)*0.25); // halfway between green and brown, as a guess; mix in 75% of pot color
	case TYPE_CLOSET:   return (color*0.5 + WHITE*0.5); // half white door and half wall color
	case TYPE_DRESSER:  return  get_textured_wood_color();
	case TYPE_NIGHTSTAND:return get_textured_wood_color();
	case TYPE_FLOORING: return texture_color(MARBLE_TEX).modulate_with(color);
	case TYPE_CRATE:    return texture_color(get_crate_tid(*this)).modulate_with(color);
	case TYPE_BOX:      return texture_color(get_box_tid()).modulate_with(color);
	case TYPE_CUBICLE:  return texture_color(get_cubicle_tid(*this));
	case TYPE_SHELVES:  return (WHITE*0.75 + get_textured_wood_color()*0.25); // mostly white walls (sparse), with some wood mixed in
	case TYPE_KEYBOARD: return BKGRAY;
	case TYPE_COMPUTER: return BKGRAY;
	case TYPE_MWAVE:    return GRAY;
	case TYPE_SHOWER:   return colorRGBA(WHITE, 0.25); // partially transparent - does this actually work?
	case TYPE_BLINDS:   return texture_color(get_blinds_tid()).modulate_with(color);
	case TYPE_LG_BALL:  return texture_color(get_lg_ball_tid(*this));
	case TYPE_HANGER_ROD:return get_textured_wood_color();
	default: return color; // TYPE_LIGHT, TYPE_TCAN, TYPE_BOOK, TYPE_BOTTLE, TYPE_PEN_PENCIL, etc.
	}
	if (is_obj_model_type()) {return color.modulate_with(building_obj_model_loader.get_avg_color(get_model_id()));} // handle models
	return color; // Note: probably should always set color so that we can return it here
}

void building_room_geom_t::create_static_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom"); // 3.3ms
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
		case TYPE_TCAN:    add_trashcan(*i); break;
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
	//highres_timer_t timer("Gen Room Geom Small"); // 5.6ms
	mats_small .clear();
	mats_plants.clear();
	add_small_static_objs_to_verts(objs);
	add_small_static_objs_to_verts(expanded_objs);
	mats_small .create_vbos(building);
	mats_plants.create_vbos(building);
}
void building_room_geom_t::add_small_static_objs_to_verts(vector<room_object_t> const &objs_to_add) {
	float const tscale(2.0/obj_scale);

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
		case TYPE_PEN: case TYPE_PENCIL: add_pen_pencil(*i); break;
		case TYPE_LG_BALL:   add_lg_ball  (*i); break;
		case TYPE_HANGER_ROD:add_hanger_rod(*i); break;
		case TYPE_DRAIN:     add_drain_pipe(*i); break;
		case TYPE_KEY:       add_key       (*i); break;
		default: break;
		}
	} // for i
}
void building_room_geom_t::create_obj_model_insts(building_t const &building) { // handle drawing of 3D models
	obj_model_insts.clear();
	auto objs_end(objs.begin() + stairs_start); // skip stairs and elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_visible() || !i->is_obj_model_type()) continue;
		vector3d dir(zero_vector);
		dir[i->dim] = (i->dir ? 1.0 : -1.0);

		if (i->flags & RO_FLAG_RAND_ROT) {
			float const angle(123.4*i->x1() + 456.7*i->y1() + 567.8*i->z1()); // random rotation angle based on position
			vector3d const rand_dir(vector3d(sin(angle), cos(angle), 0.0).get_norm());
			dir = ((dot_product(rand_dir, dir) < 0.0) ? -rand_dir : rand_dir); // random, but facing in the correct general direction
		}
		if (building.is_rotated()) {building.do_xy_rotate_normal(dir);}
		obj_model_insts.emplace_back((i - objs.begin()), i->get_model_id(), dir);
		//get_material(tid_nm_pair_t()).add_cube_to_verts(*i, WHITE, tex_origin); // for debugging of model bcubes
	} // for i
}
void building_room_geom_t::create_lights_vbos(building_t const &building) {
	//highres_timer_t timer("Gen Room Geom Light"); // 0.3ms
	float const tscale(2.0/obj_scale);
	auto objs_end(objs.begin() + stairs_start); // skip stairs and elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->is_visible() && i->type == TYPE_LIGHT) {add_light(*i, tscale);}
	}
	mats_lights.create_vbos(building);
}
void building_room_geom_t::create_dynamic_vbos(building_t const &building) {
	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible() || !i->is_dynamic()) continue; // only visible + dynamic objects; can't do VFC because this is not updated every frame
		switch (i->type) {
		case TYPE_ELEVATOR: add_elevator(*i, 2.0/obj_scale); break;
		case TYPE_LG_BALL:  add_lg_ball (*i); break;
		default: assert(0); // not a supported dynamic object type
		}
	} // for i
	mats_dynamic.create_vbos(building);
}

void building_room_geom_t::expand_object(room_object_t &c) {
	if (c.flags & RO_FLAG_EXPANDED) return; // already expanded
	switch (c.type) {
	case TYPE_BCASE:     expand_bookcase (c); break;
	case TYPE_CLOSET:    expand_closet   (c); break;
	case TYPE_SHELVES:   expand_shelves  (c); break;
	case TYPE_WINE_RACK: expand_wine_rack(c); break;
	default: assert(0); // not a supported expand type
	}
	c.flags |= RO_FLAG_EXPANDED; // flag as expanded
}

struct occlusion_stats_t {
	unsigned nobj, next, nnear, nvis, ndraw;
	int last_frame_counter;
	occlusion_stats_t() : last_frame_counter(0) {reset();}
	void reset() {nobj = next = nnear = nvis = ndraw = 0;}

	void update() {
		if (frame_counter == last_frame_counter) return; // same frame
		last_frame_counter = frame_counter;
		cout << TXT(nobj) << TXT(next) << TXT(nnear) << TXT(nvis) << TXT(ndraw) << endl;
		reset();
	}
};

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

occlusion_stats_t occlusion_stats;

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
	enable_blend(); // needed for rugs and book text
	mats_static .draw(s, shadow_only, reflection_pass);
	mats_lights .draw(s, shadow_only, reflection_pass);
	mats_dynamic.draw(s, shadow_only, reflection_pass);

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
	//if (!obj_model_insts.empty()) {glDisable(GL_CULL_FACE);} // better but slower?
	point const camera_bs(camera_pdu.pos - xlate), building_center(building.bcube.get_cube_center());
	bool const is_rotated(building.is_rotated());
	oc.set_exclude_bix(building_ix);
	bool obj_drawn(0);

	// draw object models
	for (auto i = obj_model_insts.begin(); i != obj_model_insts.end(); ++i) {
		auto &obj(get_room_object_by_index(i->obj_id));
		//++occlusion_stats.nobj;
		if (!player_in_building && obj.is_interior()) continue; // don't draw objects in interior rooms if the player is outside the building (useful for office bathrooms)
		//++occlusion_stats.next;
		point obj_center(obj.get_cube_center());
		if (is_rotated) {building.do_xy_rotate(building_center, obj_center);}
		if (!shadow_only && !dist_less_than(camera_bs, obj_center, 100.0*obj.dz())) continue; // too far away
		//++occlusion_stats.nnear;
		if (!(is_rotated ? building.is_rot_cube_visible(obj, xlate) : camera_pdu.cube_visible(obj + xlate))) continue; // VFC
		//++occlusion_stats.nvis;
		if ((display_mode & 0x08) && building.check_obj_occluded(obj, camera_bs, oc, reflection_pass)) continue;
		//++occlusion_stats.ndraw;
		bool const is_emissive(i->model_id == OBJ_MODEL_LAMP && obj.is_lit());
		if (is_emissive) {s.set_color_e(LAMP_COLOR*0.4);}
		apply_room_obj_rotate(obj, *i); // Note: may modify obj by clearing flags
		building_obj_model_loader.draw_model(s, obj_center, obj, i->dir, obj.color, xlate, i->model_id, shadow_only, 0, 0);
		if (is_emissive) {s.set_color_e(BLACK);}
		obj_drawn = 1;
	} // for i
	//occlusion_stats.update();
	//if (!obj_model_insts.empty()) {glEnable(GL_CULL_FACE);}
	if (obj_drawn) {check_mvm_update();} // needed after popping model transform matrix

	if (player_in_building && !shadow_only && player_held_object.is_valid()) { // draw the item the player is holding
		assert(player_held_object.type == TYPE_LG_BALL); // this is currently the only supported object type
		player_held_object.translate((camera_bs + CAMERA_RADIUS*cview_dir - vector3d(0.0, 0.0, 0.5*CAMERA_RADIUS)) - player_held_object.get_cube_center());
		bind_vbo(0); // not using a VBO here
		draw_lg_ball_in_building(player_held_object, s);
	}
	if (!shadow_only && !mats_alpha.empty()) { // draw last; not shadow casters
		enable_blend();
		mats_alpha.draw(s, shadow_only, reflection_pass);
		disable_blend();
		vbo_wrap_t::post_render();
	}
}

bool are_pts_occluded_by_any_cubes(point const &pt, point const *const pts, unsigned npts, vect_cube_t const &cubes, unsigned dim) {
	assert(npts > 0 && dim <= 2);

	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		// use dim as an optimization?
		if (!check_line_clip(pt, pts[0], c->d)) continue; // first point does not intersect
		bool not_occluded(0);

		for (unsigned p = 1; p < npts; ++p) { // skip first point
			if (!check_line_clip(pt, pts[p], c->d)) {not_occluded = 1; break;}
		}
		if (!not_occluded) return 1;
	} // for c
	return 0;
}

bool building_t::check_obj_occluded(cube_t const &c, point const &viewer, occlusion_checker_noncity_t &oc, bool reflection_pass) const {
	if (!interior)    return 0; // could probably make this an assert
	if (is_rotated()) return 0; // TODO: implement rotated building occlusion culling; cubes are not actually cubes; seems messy
	//highres_timer_t timer("Check Object Occlusion");
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
		if (player_building != nullptr && player_building->interior) { // check walls of the building the player is in
			assert(player_building != this); // otherwise player_in_this_building should be true
			
			for (unsigned d = 0; d < 2; ++d) { // check walls of the building the player is in
				if (are_pts_occluded_by_any_cubes(viewer, pts, npts, player_building->interior->walls[d], d)) return 1;
			}
		}
	}
	else { // player not in a building
		if (oc.is_occluded(c)) return 1; // check other buildings
	}
	return 0;
}


