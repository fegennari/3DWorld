// 3D World - Building Interior Room Geometry
// by Frank Gennari 4/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#pragma warning(disable : 26812) // prefer enum class over enum


extern int display_mode;

int get_rand_screenshot_texture(unsigned rand_ix);
unsigned get_num_screenshot_tids();


colorRGBA const WOOD_COLOR(0.9, 0.7, 0.5); // light brown, multiplies wood texture color

// skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2 to match CSG cube flags
void rgeom_mat_t::add_cube_to_verts(cube_t const &c, colorRGBA const &color, unsigned skip_faces, bool swap_tex_st, bool mirror_x) {
	vertex_t v;
	v.set_c4(color);

	// Note: stolen from draw_cube() with tex coord logic, back face culling, etc. removed
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (skip_faces & (1 << (2*(2-n) + j))) continue; // skip this face
			v.set_ortho_norm(n, j);
			v.v[n] = c.d[n][j];

			for (unsigned s1 = 0; s1 < 2; ++s1) {
				v.v[d[1]] = c.d[d[1]][s1];
				v.t[swap_tex_st] = ((tex.tscale_x == 0.0) ? float(s1) : tex.tscale_x*v.v[d[1]]); // tscale==0.0 => fit texture to cube

				for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
					bool const s2(k^j^s1^1); // need to orient the vertices differently for each side
					v.v[d[0]] = c.d[d[0]][s2];
					v.t[!swap_tex_st] = ((tex.tscale_y == 0.0) ? float(s2) : tex.tscale_y*v.v[d[0]]);
					quad_verts.push_back(v);
					if (mirror_x) {quad_verts.back().t[0] = 1.0 - v.t[0];} // use for pictures
				} // for k
			} // for s1
		} // for j
	} // for i
}

void rgeom_mat_t::add_vcylin_to_verts(cube_t const &c, colorRGBA const &color) {
	float const radius(0.5*min(c.dx(), c.dy())); // should be equal/square
	point const center(c.get_cube_center());
	point const ce[2] = {point(center.x, center.y, c.z1()), point(center.x, center.y, c.z2())};
	unsigned const ndiv(N_CYL_SIDES);
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius, radius, ndiv, v12));
	float const ndiv_inv(1.0/ndiv);
	unsigned qix(quad_verts.size()), tix(tri_verts.size());
	quad_verts.resize(qix + 4*ndiv);
	color_wrapper cw(color);

	for (unsigned i = 0; i < ndiv; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			unsigned const S(i + j), s(S%ndiv);
			float const ts(1.0f - S*ndiv_inv);
			vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]); // normalize?
			quad_verts[qix++].assign(vpn.p[(s<<1)+ j], normal, ts, 1.0*( j), cw.c);
			quad_verts[qix++].assign(vpn.p[(s<<1)+!j], normal, ts, 1.0*(!j), cw.c);
		}
	} // for i
	assert(qix == quad_verts.size());
	// add bottom end cap using triangles, currently using all TCs=0.0
	tri_verts.resize(tix + 3*ndiv);

	for (unsigned i = 0; i < ndiv; ++i) {
		for (unsigned j = 0; j < 2; ++j) {tri_verts[tix++].assign(vpn.p[((i + j)%ndiv)<<1], -plus_z, 0.0, 0.0, cw.c);}
		tri_verts[tix++].assign(ce[0], -plus_z, 0.0, 0.0, cw.c); // center
	}
	assert(tix == tri_verts.size());
}

void rgeom_mat_t::create_vbo() {
	num_tverts = tri_verts.size();
	num_qverts = quad_verts.size();
	vector_add_to(tri_verts, quad_verts);
	vbo.create_and_upload(quad_verts);
	clear_container(tri_verts);  // no longer needed
	clear_container(quad_verts); // no longer needed
}

void rgeom_mat_t::draw(shader_t &s, bool shadow_only) {
	if (shadow_only && tex.emissive) return; // assume this is a light source and shouldn't produce shadows
	assert(vbo.vbo_valid());
	assert(num_tverts > 0 || num_qverts > 0);
	if (!shadow_only) {tex.set_gl(s);} // ignores texture scale for now
	vbo.pre_render();
	vertex_t::set_vbo_arrays();
	if (num_qverts > 0) {draw_quads_as_tris(num_qverts);}
	if (num_tverts > 0) {glDrawArrays(GL_TRIANGLES, num_qverts, num_tverts);}
	tex.unset_gl(s);
}

void building_materials_t::clear() {
	for (iterator m = begin(); m != end(); ++m) {m->clear();}
	vector<rgeom_mat_t>::clear();
}
unsigned building_materials_t::count_all_verts() const {
	unsigned num_verts(0);
	for (const_iterator m = begin(); m != end(); ++m) {num_verts += m->num_qverts + m->num_tverts;}
	return num_verts;
}
rgeom_mat_t &building_materials_t::get_material(tid_nm_pair_t const &tex) {
	// for now we do a simple linear search because there shouldn't be too many unique materials
	for (iterator m = begin(); m != end(); ++m) {
		if (m->tex == tex) return *m;
	}
	emplace_back(tex); // not found, add a new material
	return back();
}
void building_materials_t::create_vbos() {
	for (iterator m = begin(); m != end(); ++m) {m->create_vbo();}
}
void building_materials_t::draw(shader_t &s, bool shadow_only) {
	for (iterator m = begin(); m != end(); ++m) {m->draw(s, shadow_only);}
}

void get_tc_leg_cubes(cube_t const &c, float width, cube_t cubes[4]) {
	for (unsigned y = 0; y < 2; ++y) {
		for (unsigned x = 0; x < 2; ++x) {
			cube_t leg(c);
			leg.d[0][x] += (1.0f - width)*(x ? -1.0f : 1.0f)*c.dx();
			leg.d[1][y] += (1.0f - width)*(y ? -1.0f : 1.0f)*c.dy();
			cubes[2*y+x] = leg;
		}
	}
}
void building_room_geom_t::add_tc_legs(cube_t const &c, colorRGBA const &color, float width, float tscale) {
	rgeom_mat_t &mat(get_wood_material(tscale));
	cube_t cubes[4];
	get_tc_leg_cubes(c, width, cubes);
	for (unsigned i = 0; i < 4; ++i) {mat.add_cube_to_verts(cubes[i], color, (EF_Z1 | EF_Z2));} // skip top and bottom faces
}

colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c) {
	if (display_mode & 0x10) return c; // disable this when using indir lighting
	return c * (0.5f + 0.5f*min(sqrt(o.light_amt), 1.5f)); // use c.light_amt as an approximation for ambient lighting due to sun/moon
}

void building_room_geom_t::add_table(room_object_t const &c, float tscale) { // 6 quads for top + 4 quads per leg = 22 quads = 88 verts
	cube_t top(c), legs_bcube(c);
	top.z1() += 0.85*c.dz(); // 15% of height
	legs_bcube.z2() = top.z1();
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(top, color); // all faces drawn
	add_tc_legs(legs_bcube, color, 0.08, tscale);
}

void building_room_geom_t::add_chair(room_object_t const &c, float tscale) { // 6 quads for seat + 5 quads for back + 4 quads per leg = 27 quads = 108 verts
	float const height(c.dz());
	cube_t seat(c), back(c), legs_bcube(c);
	seat.z1() += 0.32*height;
	seat.z2()  = back.z1() = seat.z1() + 0.07*height;
	legs_bcube.z2() = seat.z1();
	back.d[c.dim][c.dir] += 0.88f*(c.dir ? -1.0f : 1.0f)*c.get_sz_dim(c.dim);
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.2*tscale)).add_cube_to_verts(seat, apply_light_color(c, c.color)); // all faces drawn
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(back, color, EF_Z1); // skip bottom face
	add_tc_legs(legs_bcube, color, 0.15, tscale);
}

void building_room_geom_t::add_stair(room_object_t const &c, float tscale) {
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.5*tscale)).add_cube_to_verts(c, colorRGBA(0.85, 0.85, 0.85)); // all faces drawn
}

unsigned get_face_mask(bool dim, bool dir) {return ~(1 << (2*(2-dim) + dir));} // skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2

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
	unsigned const front_face_mask(get_face_mask(c.dim, c.dir)), floor_ceil_face_mask(front_face_mask & 60); // +Z faces
	tid_nm_pair_t const paneling(get_tex_auto_nm(PANELING_TEX, 2.0f*tscale));
	get_material(get_tex_auto_nm(TILE_TEX, tscale), 1).add_cube_to_verts(floor, WHITE, floor_ceil_face_mask);
	get_material(get_tex_auto_nm(get_rect_panel_tid(), tscale), 1).add_cube_to_verts(ceil, WHITE, floor_ceil_face_mask);
	get_material(paneling, 1).add_cube_to_verts(back, WHITE, front_face_mask, !c.dim);

	for (unsigned d = 0; d < 2; ++d) { // side walls
		cube_t side(c);
		side.d[!c.dim][!d] = side.d[!c.dim][d] + (d ? -1.0 : 1.0)*thickness;
		get_material(paneling, 1).add_cube_to_verts(side, WHITE, get_face_mask(!c.dim, !d), c.dim);
	}
}

void building_room_geom_t::add_light(room_object_t const &c, float tscale) {
	// Note: need to use a different texture (or -1) for is_on because emissive flag alone does not cause a material change
	bool const is_on(c.is_lit());
	tid_nm_pair_t tp((is_on ? (int)WHITE_TEX : (int)PLASTER_TEX), tscale);
	tp.emissive = is_on;
	rgeom_mat_t &mat(get_material(tp));
	if      (c.shape == SHAPE_CUBE ) {mat.add_cube_to_verts  (c, c.color, EF_Z2);} //untextured, skip top face
	else if (c.shape == SHAPE_CYLIN) {mat.add_vcylin_to_verts(c, c.color);}
	else {assert(0);}
}

void building_room_geom_t::add_rug(room_object_t const &c) {
	bool const swap_tex_st(c.dy() < c.dx()); // rug textures are oriented with the long side in X, so swap the coordinates (rotate 90 degrees) if our rug is oriented the other way
	get_material(tid_nm_pair_t(c.get_rug_tid(), 0.0)).add_cube_to_verts(c, WHITE, 61, swap_tex_st); // only draw top/+z face
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
	unsigned skip_faces(~(1 << (2*(2-c.dim) + c.dir))); // only the face oriented outward
	bool const mirror_x(!whiteboard && !(c.dim ^ c.dir));
	get_material(tid_nm_pair_t(picture_tid, 0.0)).add_cube_to_verts(c, WHITE, skip_faces, !c.dim, mirror_x);
	// add a frame
	cube_t frame(c);
	vector3d exp;
	exp.z = exp[!c.dim] = (whiteboard ? 0.04 : 0.06)*c.dz(); // frame width
	exp[c.dim] = -0.1*c.get_sz_dim(c.dim); // shrink in this dim
	frame.expand_by(exp);
	get_material(tid_nm_pair_t()).add_cube_to_verts(frame, (whiteboard ? GRAY : BLACK), skip_faces, 0);
	
	if (whiteboard) { // add a marker ledge
		cube_t ledge(c);
		ledge.z2() = ledge.z1() + 0.016*c.dz(); // along the bottom edge
		ledge.d[c.dim][c.dir] += (c.dir ? 1.5 : -1.5)*c.get_sz_dim(c.dim); // extrude outward
		get_material(tid_nm_pair_t()).add_cube_to_verts(ledge, GRAY, (1 << (2*(2-c.dim) + !c.dir)), 0);
	}
}

void building_room_geom_t::add_book(room_object_t const &c) {
	// TODO - WRITE
}
void building_room_geom_t::add_bookcase(room_object_t const &c) {
	// TODO - WRITE
}
void building_room_geom_t::add_desk(room_object_t const &c) {
	// TODO - WRITE
}
void building_room_geom_t::add_trashcan(room_object_t const &c) {
	// TODO - WRITE
}

void building_room_geom_t::clear() {
	clear_materials();
	objs.clear();
	light_bcubes.clear();
	has_elevators = 0;
}
void building_room_geom_t::clear_materials() { // can be called to update textures, lighting state, etc.
	materials_s.clear();
	materials_d.clear();
}

rgeom_mat_t &building_room_geom_t::get_wood_material(float tscale) {
	return get_material(get_tex_auto_nm(WOOD2_TEX, tscale)); // hard-coded for common material
}

colorRGBA room_object_t::get_color() const {
	switch (type) {
	case TYPE_TABLE:    return WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX));
	case TYPE_CHAIR:    return (color + WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX)))*0.5; // 50% seat color / 50% wood legs color
	case TYPE_STAIR:    return LT_GRAY; // close enough
	case TYPE_ELEVATOR: return LT_BROWN; // ???
	case TYPE_RUG:      return texture_color(get_rug_tid());
	case TYPE_PICTURE:  return texture_color(get_picture_tid());
	case TYPE_WBOARD:   return WHITE;
	case TYPE_BCASE:    return WOOD_COLOR;
	case TYPE_DESK:     return WOOD_COLOR;
	case TYPE_TCAN:     return BLACK;
	default: return color;
	}
	return color; // Note: probably should always set color so that we can return it here
}

void building_room_geom_t::create_static_vbos() {
	float const tscale(2.0/obj_scale);

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible()) continue;
		assert(i->is_strictly_normalized());
		switch (i->type) {
		case TYPE_NONE:  assert(0); // not supported
		case TYPE_TABLE:   add_table   (*i, tscale); break;
		case TYPE_CHAIR:   add_chair   (*i, tscale); break;
		case TYPE_STAIR:   add_stair   (*i, tscale); break;
		case TYPE_LIGHT:   add_light   (*i, tscale); break; // light fixture
		case TYPE_RUG:     add_rug     (*i); break;
		case TYPE_PICTURE: add_picture (*i); break;
		case TYPE_WBOARD:  add_picture (*i); break;
		case TYPE_BOOK:    add_book    (*i); break;
		case TYPE_BCASE:   add_bookcase(*i); break;
		case TYPE_DESK:    add_desk    (*i); break;
		case TYPE_TCAN:    add_trashcan(*i); break;
		case TYPE_ELEVATOR: break; // not handled here
		default: assert(0); // undefined type
		}
	} // for i
	// Note: verts are temporary, but cubes are needed for things such as collision detection with the player and ray queries for indir lighting
	materials_s.create_vbos();
}

void building_room_geom_t::create_dynamic_vbos() {
	if (!has_elevators) return; // currently only elevators are dynamic, can skip this step if there are no elevators

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible() || i->type != TYPE_ELEVATOR) continue; // only elevators for now
		add_elevator(*i, 2.0/obj_scale);
	}
	materials_d.create_vbos();
}

void building_room_geom_t::draw(shader_t &s, bool shadow_only) { // non-const because it creates the VBO
	if (empty()) return; // no geom
	unsigned const num_screenshot_tids(get_num_screenshot_tids());

	if (has_pictures && num_pic_tids != num_screenshot_tids) {
		clear_materials(); // user created a new screenshot texture, and this building has pictures - recreate room geom
		num_pic_tids = num_screenshot_tids;
	}
	if (materials_s.empty()) {create_static_vbos ();} // create static  materials if needed
	if (materials_d.empty()) {create_dynamic_vbos();} // create dynamic materials if needed
	enable_blend(); // needed for rugs
	materials_s.draw(s, shadow_only);
	materials_d.draw(s, shadow_only);
	disable_blend();
	vbo_wrap_t::post_render();
}


