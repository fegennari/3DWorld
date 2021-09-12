// 3D World - Building Interior Room Geometry Drawing
// by Frank Gennari 7/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "scenery.h" // for s_plant
#include "shaders.h"

bool const ADD_BOOK_COVERS = 1;
bool const ADD_BOOK_TITLES = 1;
colorRGBA const STAIRS_COLOR_TOP(0.7, 0.7, 0.7);
colorRGBA const STAIRS_COLOR_BOT(0.9, 0.9, 0.9);

vect_cube_t temp_cubes;
vector<room_object_t> temp_objects;
vect_cube_t &get_temp_cubes() {temp_cubes.clear(); return temp_cubes;}
vector<room_object_t> &get_temp_objects() {temp_objects.clear(); return temp_objects;}

extern int display_mode, player_in_closet;

int get_rand_screenshot_texture(unsigned rand_ix);
unsigned get_num_screenshot_tids();

void gen_text_verts(vector<vert_tc_t> &verts, point const &pos, string const &text, float tsize, vector3d const &column_dir, vector3d const &line_dir, bool use_quads=0);
string const &gen_book_title(unsigned rand_id, string *author, unsigned split_len);

unsigned get_face_mask(unsigned dim, bool dir) {return ~(1 << (2*(2-dim) + dir));} // skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2
unsigned get_skip_mask_for_xy(bool dim) {return (dim ? EF_Y12 : EF_X12);}
tid_nm_pair_t get_tex_auto_nm(int tid, float tscale=1.0, bool shadowed=1) {return tid_nm_pair_t(tid, get_normal_map_for_bldg_tid(tid), tscale, tscale, 0.0, 0.0, shadowed);}
int get_counter_tid    () {return get_texture_by_name("marble2.jpg");}
int get_paneling_nm_tid() {return get_texture_by_name("normal_maps/paneling_NRM.jpg", 1);}
int get_blinds_tid     () {return get_texture_by_name("interiors/blinds.jpg", 1, 0, 1, 8.0);} // use high aniso
int get_money_tid      () {return get_texture_by_name("interiors/dollar20.jpg");}
int get_crack_tid(room_object_t const &obj) {return get_texture_by_name(((5*obj.obj_id + 7*obj.room_id) & 1) ? "interiors/cracked_glass2.jpg" : "interiors/cracked_glass.jpg");}

colorRGBA get_textured_wood_color() {return WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX));} // Note: uses default WOOD_COLOR, not the per-building random variant
colorRGBA get_counter_color      () {return (get_textured_wood_color()*0.75 + texture_color(get_counter_tid())*0.25);}

rgeom_mat_t &building_room_geom_t::get_wood_material(float tscale, bool inc_shadows, bool dynamic, bool small) {
	return get_material(get_tex_auto_nm(WOOD2_TEX, tscale, inc_shadows), inc_shadows, dynamic, small); // hard-coded for common material
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

void get_table_cubes(room_object_t const &c, cube_t cubes[5], bool is_desk) {
	assert(c.shape != SHAPE_CYLIN); // can't call this on cylindrical table
	cube_t top(c), legs_bcube(c);
	top.z1() += (is_desk ? 0.85 : 0.88)*c.dz(); // Note: default table with top_dz=0.12, leg_width=0.08; desk is 0.15/0.06
	legs_bcube.z2() = top.z1();
	cubes[0] = top;
	get_tc_leg_cubes(legs_bcube, (is_desk ? 0.06 : 0.08), (cubes+1));
}
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

void building_room_geom_t::add_dresser_drawers(room_object_t const &c, float tscale) { // or nightstand
	vect_cube_t &drawers(get_temp_cubes());
	get_drawer_cubes(c, drawers, 1); // front_only=1
	assert(drawers.size() <= 16); // we only have 16 bits to store drawer flags
	float const height(c.dz()), drawer_thick(0.05*height), handle_thick(0.75*drawer_thick), dir_sign(c.dir ? 1.0 : -1.0), handle_width(0.07*height);
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
			left .d[!c.dim][1]    -= 0.87*dwidth; // set width of left  side
			right.d[!c.dim][0]    += 0.87*dwidth; // set width of right side
			back.d[c.dim][ c.dir]  = c.d[c.dim][c.dir] + 0.25f*dir_sign*drawer_thick; // flush with front face and narrow
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
			cube_t interior(drawer_body);
			interior.z1() = bottom.z2();
			interior.d[!c.dim][0] = left .d[!c.dim][1];
			interior.d[!c.dim][1] = right.d[!c.dim][0];
			interior.d[ c.dim][!c.dir] = back.d[c.dim][c.dir];
			room_object_t const obj(get_item_in_drawer(c, interior, (i - drawers.begin())));
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
void get_closet_cubes(room_object_t const &c, cube_t cubes[5], bool for_collision) {
	float const width(c.get_sz_dim(!c.dim)), depth(c.get_sz_dim(c.dim)), height(c.dz());
	bool const use_small_door(c.is_small_closet()), doors_fold(!use_small_door && (c.flags & RO_FLAG_HANGING));
	// small closets: door does not collide when open; large closets: edges of door still collide even when open
	float const wall_width(use_small_door ? 0.5*(width - 0.5*height) : ((for_collision && c.is_open()) ? (doors_fold ? 0.2 : 0.25) : 0.05)*width);
	float const wall_thick(WALL_THICK_VAL*(1.0f - FLOOR_THICK_VAL_HOUSE)*height), wall_shift(width - wall_width);
	assert(wall_shift > 0.0);
	cube_t doors(c), walls[2] = {c, c}; // left, right
	walls[0].d[!c.dim][1] -= wall_shift;
	walls[1].d[!c.dim][0] += wall_shift;

	for (unsigned d = 0; d < 2; ++d) {
		cube_t front(walls[d]), side(walls[d]);
		front.d[c.dim][!c.dir] = side.d[c.dim][c.dir] = front.d[c.dim][c.dir] - (c.dir ? 1.0 : -1.0)*wall_thick; // set front wall/door thickness

		if (!(c.flags & (d ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO))) { // only need to draw side wall when not adjacent to room wall
			side.d[!c.dim][!d] += (d ? 1.0f : -1.0f)*(wall_width - wall_thick);
			cubes[2*d+1] = side; // left or right side
		}
		cubes[2*d] = front; // front left or front right
		doors.d[!c.dim][d] = walls[d].d[!c.dim][!d]; // clip door to space between walls
	} // for d
	doors.d[c.dim][ c.dir] -= (c.dir ? 1.0 : -1.0)*0.2*wall_thick; // shift in slightly
	doors.d[c.dim][!c.dir] += (c.dir ? 1.0 : -1.0)*(depth - 0.8*wall_thick); // make it narrow
	cubes[4] = doors; // Note: this is for closed door; caller must handle open door
}

void add_quad_to_mat(rgeom_mat_t &mat, point const pts[4], float const ts[4], float const tt[4], color_wrapper const &cw) {
	norm_comp normal(get_poly_norm(pts));
	for (unsigned n = 0; n < 4; ++n) {mat.quad_verts.emplace_back(pts[n], normal, ts[n], tt[n], cw);}
}

void building_room_geom_t::add_closet(room_object_t const &c, tid_nm_pair_t const &wall_tex, bool inc_lg, bool inc_sm) { // no lighting scale, houses only
	bool const open(c.is_open()), use_small_door(c.is_small_closet()), draw_interior(open || player_in_closet);
	cube_t cubes[5];
	get_closet_cubes(c, cubes);

	if (inc_lg) { // draw closet walls and doors
		rgeom_mat_t &wall_mat(get_material(get_scaled_wall_tex(wall_tex), 1));
		// need to draw the face that's against the wall for the shadow pass if the closet light is on, if the player is in the closet, or if the doors are open
		unsigned const skip_faces(EF_Z12); // skip top and bottom

		for (unsigned d = 0; d < 2; ++d) {
			bool const adj_room_wall(c.flags & (d ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));
			unsigned const extra_skip_faces(adj_room_wall ? ~get_face_mask(!c.dim, d) : 0); // adjacent to room wall, skip that face
			// only need to draw side wall when not adjacent to room wall; skip front face of side wall
			if (!adj_room_wall) {wall_mat.add_cube_to_verts(cubes[2*d+1], c.color, tex_origin, (skip_faces | extra_skip_faces | get_skip_mask_for_xy(c.dim)));}
			unsigned const front_wall_skip_flags((draw_interior ? EF_Z12 : skip_faces) | extra_skip_faces);
			wall_mat.add_cube_to_verts(cubes[2*d], c.color, tex_origin, front_wall_skip_flags); // Note: c.color should be wall color
		} // for d
		cube_t const &doors(cubes[4]);
		point const llc(doors.get_llc());
		float const out_sign(c.dir ? 1.0 : -1.0);

		if (use_small_door) { // small house closet door
			tid_nm_pair_t const tp(get_int_door_tid(), 0.0);

			if (open) {
				float const door_width(doors.get_sz_dim(!c.dim)), door_thickness(doors.get_sz_dim(c.dim));
				cube_t door(doors);
				door.d[ c.dim][c.dir] += out_sign*door_width;
				door.d[!c.dim][1    ] -= (door_width - door_thickness);
				get_material(tp, 1)       .add_cube_to_verts(door, WHITE, llc, ~get_skip_mask_for_xy(!c.dim), c.dim, !c.dir); // draw front and back faces
				get_untextured_material(1).add_cube_to_verts(door, WHITE, llc, ~get_skip_mask_for_xy( c.dim)); // draw edges untextured
			}
			else {
				unsigned const door_skip_faces(get_face_mask(c.dim, c.dir) & (player_in_closet ? get_face_mask(c.dim, !c.dir) : 0xFF));
				get_material(tp, 1).add_cube_to_verts(doors, WHITE, llc, door_skip_faces, !c.dim); // draw only front face, back face if player in closet
			}
		}
		else { // 4 panel folding door
			float const doors_width(doors.get_sz_dim(!c.dim)), door_thickness(doors.get_sz_dim(c.dim));
			float const door_spacing(0.25*doors_width), door_gap(0.01*door_spacing);
			int const tid(get_rect_panel_tid());
			float tx(1.0/doors_width), ty(0.25/doors.dz());
			if (!c.dim) {swap(tx, ty);} // swap so that ty is always in Z
			tid_nm_pair_t const door_tex(tid, get_normal_map_for_bldg_tid(tid), tx, ty); // 4x1 panels
			rgeom_mat_t &door_mat(get_material(door_tex, 1));
			bool const doors_fold(c.flags & RO_FLAG_HANGING); // else doors slide

			if (doors_fold && open) { // draw open bifold doors open on both
				// Note: this doesn't always look correct because doors can intersect other objects such as lights and dressers, and they have no edge quads
				float const panel_len(0.2*doors.get_sz_dim(!c.dim) - 2.0*door_gap), open_amt(0.5*panel_len), extend(sqrt(panel_len*panel_len - open_amt*open_amt));
				float const nom_pos(doors.d[c.dim][!c.dir]), front_pos(nom_pos + out_sign*extend), z1(doors.z1()), z2(doors.z2());
				float const ts[4] = {0.0, 0.25, 0.25, 0.0}, tt[4] = {0.0, 0.0, 0.25, 0.25};
				color_wrapper const cw(WHITE);
				point side_pt, out_pt, inner_pt; // left side door points in this order from left to right, forming a V-shape pointing outward
				side_pt[c.dim] = inner_pt[c.dim] = nom_pos;
				out_pt [c.dim] = front_pos;

				for (unsigned side = 0; side < 2; ++side) {
					float const open_sign(side ? -1.0 : 1.0), side_pos(doors.d[!c.dim][side]);
					side_pt [!c.dim] = side_pos;
					out_pt  [!c.dim] = side_pos + open_sign*open_amt;
					inner_pt[!c.dim] = side_pos + 2*open_sign*open_amt;
					// outside faces
					point pts1o[4] = {point(side_pt.x, side_pt.y, z1), point(out_pt.x,   out_pt.y,   z1), point(out_pt.x,   out_pt.y,   z2), point(side_pt.x, side_pt.y, z2)};
					point pts2o[4] = {point(out_pt.x,  out_pt.y,  z1), point(inner_pt.x, inner_pt.y, z1), point(inner_pt.x, inner_pt.y, z2), point(out_pt.x,  out_pt.y,  z2)};
					if (!c.dim) {std::reverse(pts1o, pts1o+4); std::reverse(pts2o, pts2o+4);} // reverse the winding order
					add_quad_to_mat(door_mat, pts1o, ts, tt, cw);
					add_quad_to_mat(door_mat, pts2o, ts, tt, cw);
					// inside faces
					point pts1i[4], pts2i[4];
					
					for (unsigned n = 0; n < 4; ++n) { // create inside surfaces of doors with inverted winding order and normal
						pts1i[n] = pts1o[3-n]; pts1i[n][c.dim] += open_sign*door_thickness;
						pts2i[n] = pts2o[3-n]; pts2i[n][c.dim] += open_sign*door_thickness;
					}
					add_quad_to_mat(door_mat, pts1i, ts, tt, cw);
					add_quad_to_mat(door_mat, pts2i, ts, tt, cw);
					point edge[4] = {pts2i[2], pts2i[1], pts2o[2], pts2o[1]}; // the 4 inner points for the two sides of the door; independent of reverse
					float const tc_zeros[4] = {0,0,0,0};
					add_quad_to_mat(door_mat, edge, tc_zeros, tc_zeros, cw); // Note: edge is at an odd angle, not perpendicular to the door; is this good enough?
				} // for side
			}
			else { // draw closet door in 4 cube panels
				float const mid_gap(doors_fold ? door_gap : 0.0), gaps[5] = {door_gap, mid_gap, door_gap, mid_gap, door_gap};
				unsigned const open_n[4] = {0, 0, 3, 3};

				for (unsigned n = 0; n < 4; ++n) {
					unsigned const N(open ? open_n[n] : n);
					cube_t door(doors);
					door.d[!c.dim][0] = doors.d[!c.dim][0] +  N   *door_spacing + gaps[N  ]; // left  edge
					door.d[!c.dim][1] = doors.d[!c.dim][0] + (N+1)*door_spacing - gaps[N+1]; // right edge
					if (!doors_fold && (n == 1 || n == 2)) {door.translate_dim(c.dim, -1.1*out_sign*door_thickness);} // inset the inner sliding doors
					door_mat.add_cube_to_verts(door, WHITE, llc, skip_faces);
				}
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
				get_untextured_material(0, 0, 1).add_cube_to_verts_untextured(btrim, trim_color, skip_faces); // is_small, untextured, no shadows; both interior and exterior
				set_cube_zvals(trim, c.z2()-trim_height, c.z2());
				add_wall_trim(room_object_t(trim, TYPE_WALL_TRIM, c.room_id, trim_dim, trim_dir, 0, 1.0, SHAPE_ANGLED, trim_color)); // ceiling trim, missing end caps; exterior only
			} // for d
		} // for is_side
		// Note: always drawn to avoid recreating all small objects when the player opens/closes a closet door, and so that objects can be seen through the cracks in the doors
		if (!(c.flags & RO_FLAG_EXPANDED)) { // add boxes if not expanded
			vector<room_object_t> &objects(get_temp_objects());
			add_closet_objects(c, objects);
			add_small_static_objs_to_verts(objects);
		}
	} // end inc_sm
}

void building_room_geom_t::add_hanger_rod(room_object_t const &c) { // is_small=1
	get_wood_material(1.0, 1, 0, 1).add_ortho_cylin_to_verts(c, LT_GRAY, !c.dim, 0, 0, 0, 0, 1.0, 1.0, 0.25, 1.0, 0, 16, 0.0, 1); // 16 sides, swap_txy=1
}

void building_room_geom_t::add_drain_pipe(room_object_t const &c) { // is_small=1
	rgeom_mat_t &mat(get_untextured_material(0, 0, 1)); // unshadowed, small
	mat.add_vcylin_to_verts(c, apply_light_color(c), 0, 0); // draw sides only
	mat.add_disk_to_verts(point(c.xc(), c.yc(), c.z2()), 0.5*c.dx(), 0, BLACK); // draw top as black
}

void building_room_geom_t::add_key(room_object_t const &c) { // is_small=1
	rgeom_mat_t &key_mat(get_metal_material(0, 0, 1)); // untextured, unshadowed, small=1
	key_mat.add_cube_to_verts(c, apply_light_color(c)); // placeholder for case where there's no key 3D model
}

void building_room_geom_t::add_money(room_object_t const &c) { // is_small=1
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_money_tid(), 0.0), 0, 0, 1));
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, ~EF_Z2, c.dim, (c.dim ^ c.dir ^ 1), c.dir); // top face only, no shadows
	unsigned const verts_start(mat.quad_verts.size());
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, EF_Z12); // sides, no shadows
	for (auto i = mat.quad_verts.begin() + verts_start; i != mat.quad_verts.end(); ++i) {i->t[0] = i->t[1] = 0.0;} // set tex coords to 0 for sides to use border texture color
}

tid_nm_pair_t get_phone_tex(room_object_t const &c) {
	bool const is_ringing(c.flags & RO_FLAG_EMISSIVE), is_locked(c.flags & RO_FLAG_OPEN);
	if (!is_ringing && !is_locked) {return tid_nm_pair_t();} // off - untextured black
	// disable texture compression since it looks bad
	tid_nm_pair_t tp(get_texture_by_name((is_ringing ? "interiors/phone_ring_screen.jpg" : "interiors/phone_lock_screen.jpg"), 0, 0, 0, 0.0, 0), 0.0);
	tp.emissive = 1.0; // lit
	return tp;
}
void building_room_geom_t::add_phone(room_object_t const &c) { // is_small=1
	rgeom_mat_t &mat(get_untextured_material(0, 0, 1));
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, EF_Z12); // sides, no shadows

	if (c.flags & (RO_FLAG_EMISSIVE | RO_FLAG_OPEN)) { // ringing or locked screen: top, no shadows, lit
		get_material(get_phone_tex(c), 0, 0, 1).add_cube_to_verts(c, WHITE, zero_vector, ~EF_Z2, c.dim, (c.dim ^ c.dir ^ 1), c.dir);
	}
	else {mat.add_cube_to_verts(c, apply_light_color(c, BLACK), zero_vector, ~EF_Z2);} // top, no shadows, unlit
}

void add_tproll_to_material(room_object_t const &c, rgeom_mat_t &mat) {
	colorRGBA const tp_color(apply_light_color(c));
	float const radius(0.5*c.dz()), rod_shrink(-0.7*radius), roll_shrink(0.75*rod_shrink*fract(123.456*c.obj_id)); // randomly partially empty (25-100%)
	cube_t roll(c);
	roll.expand_in_dim(c.dim, roll_shrink);
	roll.expand_in_dim(2,     roll_shrink); // z
	mat.add_ortho_cylin_to_verts(roll, tp_color, !c.dim, 1, 1); // c.dim applies to the wall; the roll is oriented perpendicular
	cube_t square(roll); // hanging square of TP
	set_cube_zvals(square, c.z1(), c.zc()); // starts at the centerline (tangent) and extends to the bottom
	if (c.flags & RO_FLAG_HANGING) {square.z1() -= 3.0*c.dz();} // player has pulled it down lower
	square.d[c.dim][c.dir] = square.d[c.dim][!c.dir]; // shrink to zero thickness at outer edge
	mat.add_cube_to_verts(square, tp_color, zero_vector, ~get_skip_mask_for_xy(c.dim)); // only draw front/back faces
}
void building_room_geom_t::add_vert_tproll_to_material(room_object_t const &c, rgeom_mat_t &mat, float sz_ratio) {
	cube_t hole(c);
	hole.expand_by_xy(-0.3*c.dx());
	cube_t tube(hole);
	mat.add_vcylin_to_verts(tube, apply_light_color(c, LT_BROWN), 0, 0, 1); // tube, sides only, two sided (only need inside)
	if (sz_ratio == 0.0) return; // empty, tube only, don't need to draw the rest of the roll
	cube_t roll(c);
	if (sz_ratio < 1.0) {roll.expand_by_xy(-0.3*(1.0 - sz_ratio)*c.dx());} // partially used
	hole.expand_in_dim(2, 0.0025*c.dz()); // expand slightly to avoid z-fighting
	mat.add_vcylin_to_verts(hole, ALPHA0, 1, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 1); // hole - top/bottom surface only to mask off the outer part of the roll
	mat.add_vcylin_to_verts(roll, apply_light_color(c),  1, 1); // paper roll
}
void building_room_geom_t::add_tproll(room_object_t const &c) { // is_small=1
	if (c.flags & RO_FLAG_WAS_EXP) { // bare TP roll from a box
		add_vert_tproll_to_material(c, get_untextured_material(1, 0, 1)); // shadowed, small
		return;
	}
	float const radius(0.5*c.dz()), rod_shrink(-0.7*radius), length(c.get_sz_dim(!c.dim));
	if (!(c.flags & RO_FLAG_TAKEN1)) {add_tproll_to_material(c, get_untextured_material(1, 0, 1));} // draw the roll if not taken
	// draw the holder attached to the wall
	rgeom_mat_t &holder_mat(get_metal_material(1, 0, 1)); // untextured, shadowed, small=1
	colorRGBA const holder_color(apply_light_color(c, GRAY));
	cube_t rod(c), plate(c);
	rod.expand_in_dim( c.dim, rod_shrink);
	rod.expand_in_dim( 2,     rod_shrink); // z
	rod.expand_in_dim(!c.dim, 0.05*length); // will go slightly outside the bounds of c
	float const rod_end(rod.d[!c.dim][0]); // arbitrarily choose lower end
	plate.expand_in_dim(2, -0.65*radius); // z
	plate.d[c.dim][!c.dir] -= (c.dir ? -1.0 : 1.0)*0.6*radius;
	plate.d[!c.dim][0] = rod_end - 0.08*length; // set thickness; will go slightly outside the bounds of c
	plate.d[!c.dim][1] = rod_end;
	holder_mat.add_ortho_cylin_to_verts(rod, holder_color, !c.dim, 0, 1);
	holder_mat.add_cube_to_verts(plate, holder_color, zero_vector, ~get_face_mask(c.dim, c.dir)); // skip the face against the wall
}

void building_room_geom_t::add_spraycan_to_material(room_object_t const &c, rgeom_mat_t &mat) {
	unsigned const dim(get_max_dim(c.get_size()));
	bool const add_bottom(dim != 2); // if on its side
	cube_t can(c), cap(c);
	can.d[dim][!c.dir] = cap.d[dim][c.dir] = (c.d[dim][c.dir] + 0.7*c.get_sz_dim(dim)*(c.dir ? -1.0 : 1.0)); // point between top of can and bottom of cap
	mat.add_ortho_cylin_to_verts(can, apply_light_color(c, DK_GRAY), dim, (add_bottom && !c.dir), (add_bottom && c.dir)); // sides + bottom (if on side)
	mat.add_ortho_cylin_to_verts(cap, apply_light_color(c), dim, c.dir, !c.dir); // sides + top
}
void building_room_geom_t::add_spraycan(room_object_t const &c) { // is_small=1
	add_spraycan_to_material(c, get_untextured_material(1, 0, 1));
}

void building_room_geom_t::add_button(room_object_t const &c) {
	bool const in_elevator(c.flags & RO_FLAG_IN_ELEV);
	tid_nm_pair_t tp; // unshadowed
	if (c.flags & RO_FLAG_IS_ACTIVE) {tp.emissive = 1.0;} // make it lit when active
	colorRGBA const color(apply_light_color(c));
	rgeom_mat_t &mat(get_material(tp, 0, in_elevator, !in_elevator)); // (in_elevator ? dynamic : small)
	if      (c.shape == SHAPE_CUBE ) {mat.add_cube_to_verts(c, color, all_zeros, ~get_face_mask(c.dim, !c.dir));} // square button
	else if (c.shape == SHAPE_CYLIN) {mat.add_ortho_cylin_to_verts(c, color, c.dim, !c.dir, c.dir);} // round button
	else {assert(0);}
	
	if (!in_elevator) { // add the frame for better color contrast
		float const expand(0.7*c.dz());
		cube_t frame(c);
		frame.d[c.dim][c.dir] -= 0.9*(c.dir ? 1.0 : -1.0)*c.get_sz_dim(c.dim); // shorten to a sliver against the elevator wall
		frame.expand_in_dim(!c.dim, expand);
		frame.expand_in_dim(2,      expand); // Z
		get_untextured_material(0, 0, 1).add_cube_to_verts(frame, apply_light_color(c, DK_GRAY), all_zeros, ~get_face_mask(c.dim, !c.dir)); // small
	}
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
	float const sz(2048), x1(12/sz), x2(576/sz), x3(1458/sz), y1(1-1667/sz), y2(1-1263/sz), y3(1-535/sz); //, x4(2032/sz), y4(1-128/sz); // Note: we don't use all parts of the texture
	unsigned verts_start(mat.quad_verts.size());
	colorRGBA const color(apply_light_color(c));
	mat.add_cube_to_verts(c, color, zero_vector, (c.is_open() ? EF_Z2 : EF_Z1)); // skip bottom face (even for stacked box?)
	assert(mat.quad_verts.size() == verts_start + 20); // there should be 5 quads (+z -x +x -y +y) / 20 verts (no -z)
	mat.quad_verts[verts_start+0].set_tc(x1, y2); // z (top or inside bottom)
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
	if (c.is_open()) {
		// draw the inside of the box
		verts_start = mat.quad_verts.size(); // update
		cube_t box(c);
		box.expand_by_xy(-0.001*c.get_size()); // slight shrink of inside of box to prevent z-fighting
		box.z1() += 0.01*c.dz(); // shrink a bit more in Z
		mat.add_cube_to_verts(box, color, zero_vector, EF_Z2, 0, 0, 0, 1); // skip top face; draw inverted
		assert(mat.quad_verts.size() == verts_start + 20); // there should be 5 quads (+z -x +x -y +y) / 20 verts (no +z)
		float const ts[4] = {x2, x3, x3, x2}, tt[4] = {y1, y1, y2, y2};

		for (unsigned side = 0; side < 5; ++side) { // make all sides use a subset of the texture that has no markings
			unsigned ix(verts_start + 4*side);
			for (unsigned n = 0; n < 4; ++n) {mat.quad_verts[ix+n].set_tc(ts[n], tt[n]);}
		}
		if (!(c.flags & RO_FLAG_WAS_EXP)) { // draw open box flaps, but not if box is in drawer/shelf/closet because there may not be space for flaps
			vector3d const box_sz(c.get_size());
			float const flap_len(0.485*min(box_sz.x, box_sz.y)); // same length in both dims; slightly less than half-width because this is the base of the triangle
			color_wrapper const cw(color);
			unsigned const up_verts[2][4] = {{0,1,0,2}, {3,2,1,3}}; // vertex indices on upward pointing outside flap edges

			for (unsigned d = 0; d < 2; ++d) { // x/y
				for (unsigned e = 0; e < 2; ++e) { // side dir
					unsigned const side_ix(2*d+e);
					bool const against_wall(c.flags & (RO_FLAG_ADJ_LO << side_ix)); // encoded in adj flags
					cube_t C(box);
					C.d[d][!e] = C.d[d][e];
					C.d[d][ e] = C.d[d][e] + (e ? 1.0 : -1.0)*(against_wall ? 0.05 : 1.0)*flap_len;
					float const zbot(C.z2()), dz(against_wall ? flap_len : 0.25*min(flap_len, box_sz.z)); // tilted somewhat upward; pointing up if against wall
					point const pts[4] = {point(C.x1(), C.y1(), zbot), point(C.x2(), C.y1(), zbot), point(C.x2(), C.y2(), zbot), point(C.x1(), C.y2(), zbot)};
					unsigned const ix(mat.quad_verts.size());
					add_quad_to_mat(mat, pts, ts, tt, cw);
					for (unsigned n = 0; n < 2; ++n) {mat.quad_verts[ix + up_verts[n][side_ix]].v.z += dz;}
					
					// add bottom surface with inverted normal in reverse order
					for (unsigned n = 0; n < 4; ++n) {
						mat.quad_verts.push_back(mat.quad_verts[ix+3-n]);
						mat.quad_verts.back().invert_normal();
					}
				} // for e
			} // for d
		}
	}
}

void building_room_geom_t::add_paint_can(room_object_t const &c) {
	float const side_tscale_add(fract(11111*c.x1() + 22222*c.y1() + 33333*c.z1())); // somewhat random
	rgeom_mat_t &side_mat(get_material(tid_nm_pair_t(get_texture_by_name("interiors/paint_can_label.png")), 1, 0, 1)); // shadows, small
	side_mat.add_vcylin_to_verts(c, apply_light_color(c), 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 24, side_tscale_add); // draw sides only; random texture rotation
	point top(c.get_cube_center());
	top.z = c.z2();
	get_metal_material(1, 0, 1).add_disk_to_verts(top, 0.5*min(c.dx(), c.dy()), 0, apply_light_color(c, LT_GRAY)); // shadowed, specular metal; small=1
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
	if (c.flags & RO_FLAG_INTERIOR) { // add support brackets to interior shelves; skip them if against an exterior wall in case they intersect a window
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
				metal_mat.add_cube_to_verts_untextured(bracket, bracket_color, (skip_faces | EF_Z2)); // skip top faces, maybe back

				if (s == 0) { // add vertical brackets on first shelf
					cube_t vbracket(bracket);
					set_cube_zvals(vbracket, c.z1(), c.z2());
					vbracket.d[c.dim][ c.dir] = c         .d[c.dim][c.dir] + (c.dir ? -1.0 : 1.0)*0.1*bracket_thickness; // nearly against the wall
					vbracket.d[c.dim][!c.dir] = shelves[s].d[c.dim][c.dir]; // against the shelf
					metal_mat.add_cube_to_verts_untextured(vbracket, bracket_color, (skip_faces | EF_Z12)); // skip top/bottom faces, maybe back
				}
				bracket.translate_dim(!c.dim, b_step);
			} // for b
		} // for s
	}
	// add objects to the shelves
	if (c.flags & RO_FLAG_EXPANDED) return; // shelves have already been expanded, don't need to create contained objects below
	vector<room_object_t> &objects(get_temp_objects());
	get_shelf_objects(c, shelves, num_shelves, objects);
	add_small_static_objs_to_verts(objects);
}

void building_room_geom_t::add_obj_with_top_texture(room_object_t const &c, string const &texture_name, colorRGBA const &sides_color, bool is_small) {
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_texture_by_name(texture_name), 0.0), 1, 0, is_small)); // shadows
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, ~EF_Z2, c.dim, (c.dim ^ c.dir ^ 1), c.dir); // top face only
	get_untextured_material(1, 0, is_small).add_cube_to_verts(c, apply_light_color(c, sides_color), zero_vector, EF_Z12); // sides, no shadows, small
}
void building_room_geom_t::add_obj_with_front_texture(room_object_t const &c, string const &texture_name, colorRGBA const &sides_color, bool is_small) {
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_texture_by_name(texture_name), 0.0), 1, 0, is_small)); // shadows
	unsigned const front_mask(get_face_mask(c.dim, c.dir));
	mat.add_cube_to_verts(c, apply_light_color(c), zero_vector, front_mask, !c.dim, (c.dim ^ c.dir ^ 1), 0); // front face only
	get_untextured_material(1, 0, is_small).add_cube_to_verts(c, apply_light_color(c, sides_color), zero_vector, ~front_mask); // sides, shadows
}

void building_room_geom_t::add_keyboard(room_object_t const &c) {add_obj_with_top_texture  (c, "interiors/keyboard.jpg",  BKGRAY, 1);} // is_small=1
void building_room_geom_t::add_laptop  (room_object_t const &c) {add_obj_with_top_texture  (c, "interiors/laptop.jpg",    BKGRAY, 1);} // is_small=1
void building_room_geom_t::add_computer(room_object_t const &c) {add_obj_with_front_texture(c, "interiors/computer.jpg",  BKGRAY, 1);} // is_small=1
void building_room_geom_t::add_mwave   (room_object_t const &c) {add_obj_with_front_texture(c, "interiors/microwave.jpg", GRAY,   0);} // is_small=0

void building_room_geom_t::add_mirror(room_object_t const &c) {
	tid_nm_pair_t tp(REFLECTION_TEXTURE_ID, 0.0);
	if (ENABLE_MIRROR_REFLECTIONS) {tp.emissive = 1.0;}
	get_material(tp, 0).add_cube_to_verts(c, c.color, zero_vector, get_face_mask(c.dim, c.dir), !c.dim); // draw only the front face
	get_untextured_material(0).add_cube_to_verts(c, apply_light_color(c), zero_vector, get_skip_mask_for_xy(c.dim)); // draw only the sides untextured
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
	bool const hdim(c.dx() < c.dy()), hdir(dirs[hdim]); // larger dim
	float const frame_width(fc.dx()), door_width(sz[!hdim]);

	if (1 || !c.is_open()) { // only draw handle if the door is closed; the math to figure out where the handle goes on the open door is complex
		bool const hside(!dirs[!hdim]);
		float const wall_dist(0.77*door_width), handle_thickness(0.8*frame_width), hdir_sign(hdir ? -1.0 : 1.0), wall_pos(c.d[hdim][!hdir]);
		cube_t handle(c);
		handle.z1() += 0.48*sz.z;
		handle.z2() -= 0.42*sz.z;

		if (c.is_open()) { // move it into the open position; the math for this is pretty complex, so it's somewhat of a trial-and error with the constants
			bool const odir(dirs[hdim]);
			float const odir_sign(odir ? -1.0 : 1.0), hside_sign(hside ? 1.0 : -1.0), inner_extend(wall_pos + odir_sign*(wall_dist - 2.0*frame_width));
			float const open_glass_pos(c.d[!hdim][!hside] + hside_sign*0.39*frame_width);
			handle.d[!hdim][ hside] = open_glass_pos;
			handle.d[!hdim][!hside] = open_glass_pos + hside_sign*handle_thickness; // outer edge
			handle.d[ hdim][!odir ] = inner_extend;
			handle.d[ hdim][ odir ] = inner_extend + odir_sign*0.03*door_width;
		}
		else { // closed
			float const hside_sign(hside ? -1.0 : 1.0), glass_pos(wall_pos - hdir_sign*(0.19*frame_width + 0.02*sz[hdir])); // place on the glass but slightly offset
			handle.d[ hdim][ hdir ]  = glass_pos;
			handle.d[ hdim][!hdir ]  = glass_pos + hdir_sign*handle_thickness; // outer edge
			handle.d[!hdim][ hside] += hside_sign*0.20*door_width; // distance to outer wall/corner
			handle.d[!hdim][!hside] -= hside_sign*wall_dist;
		}
		metal_mat.add_cube_to_verts(handle, metal_color); // draw all faces
	}
	// add drain
	cube_t drain;
	drain.set_from_point(bottom.get_cube_center());
	set_cube_zvals(drain, bottom.z2(), (bottom.z2() + 0.05*bottom.dz())); // very small height
	drain.expand_by_xy(0.06*radius); // set radius
	get_material(tid_nm_pair_t(MANHOLE_TEX, 0.0), 0).add_vcylin_to_verts(drain, metal_color, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 1); // draw top only, not sides, unshadowed
	// add transparent glass
	colorRGBA const glass_color(apply_light_color(c, colorRGBA(1.0, 1.0, 1.0, 0.25)));
	rgeom_mat_t &glass_mat(get_untextured_material(0, 0, 0, 1)); // no shadows; transparent=1

	for (unsigned d = 0; d < 2; ++d) { // for each dim
		bool const dir(dirs[d]);
		cube_t glass(fc); // start from the frame at the corner
		glass.z1()  = glass_bot;
		glass.z2()  = glass_top;
		glass.z1() += 0.002*sz.z; // to prevent z-fighting
		glass.z2() -= 0.002*sz.z; // to prevent z-fighting
		glass.d[d][!dir] = glass. d[d][ dir]; // corner point; remove overlap with frame
		glass.d[d][ dir] = fxy[d].d[d][!dir]; // edge near the wall
		glass.expand_in_dim( d, -0.01*frame_width); // to prevent z-fighting
		glass.expand_in_dim(!d, -0.20*frame_width); // set thickness

		if (bool(d) != hdim && c.is_open()) { // draw open door
			bool const odir(dirs[!d]);
			float const width(glass.get_sz_dim(d)), thickness(glass.get_sz_dim(!d)), delta(width - thickness);
			glass.d[ d][! dir] -= ( dir ? -1.0 : 1.0)*delta; // shrink width to thickness
			glass.d[!d][!odir] += (odir ? -1.0 : 1.0)*delta; // expand thickness to width
			// draw frame part of door
			float const door_frame_width(0.4*frame_width);
			cube_t top(glass), bot(glass), side(glass);
			set_cube_zvals(top, (glass.z2() + 0.01*door_frame_width), (glass.z2() + door_frame_width)); // prevent z-fighting
			set_cube_zvals(bot, (glass.z1() - door_frame_width), (glass.z1() - 0.01*door_frame_width));
			side.d[!d][!odir] = top.d[!d][!odir] = bot.d[!d][!odir] = glass.d[!d][!odir] + (odir ? -1.0 : 1.0)*door_frame_width;
			side.d[!d][ odir] = glass.d[!d][!odir]; // flush with glass on this end
			rgeom_mat_t &metal_mat2(get_metal_material(1)); // get the metal material again, in case the reference was invaldiated
			metal_mat2.add_cube_to_verts(top,  metal_color); // draw all faces
			metal_mat2.add_cube_to_verts(bot,  metal_color); // draw all faces
			metal_mat2.add_cube_to_verts(side, metal_color, all_zeros, EF_Z12); // skip top and bottom faces
		}
		glass_mat.add_cube_to_verts(glass, glass_color, zero_vector, 0, 0, 0, 0, 1); // inside surface, inverted
		glass_mat.add_cube_to_verts(glass, glass_color, zero_vector, (EF_Z1 | (d ? EF_Y12 : EF_X12))); // outside surface
	} // for d
}

void building_room_geom_t::add_bottle(room_object_t const &c) {
	// obj_id: bits 1-3 for type, bits 6-7 for emptiness, bit 6 for cap color
	// for now, no texture, but could use a bottle label texture for the central cylinder
	unsigned const bottle_ndiv = 16; // use smaller ndiv to reduce vertex count
	tid_nm_pair_t tex(-1, 1.0, 1); // shadowed
	tex.set_specular(0.5, 80.0);
	rgeom_mat_t &mat(get_material(tex, 1, 0, 1)); // inc_shadows=1, dynamic=0, small=1
	colorRGBA const color(apply_light_color(c)), cap_colors[2] = {LT_GRAY, GOLD};
	vector3d const sz(c.get_size());
	unsigned const dim(get_max_dim(sz)), dim1((dim+1)%3), dim2((dim+2)%3);
	bool const is_empty(c.is_bottle_empty()), add_bottom(dim != 2); // add bottom if bottle is on its side
	float const dir_sign(c.dir ? -1.0 : 1.0), radius(0.25f*(sz[dim1] + sz[dim2])); // base should be square (default/avg radius is 0.15*height)
	float const length(sz[dim]); // AKA height, if standing up
	cube_t sphere(c), body(c), neck(c);
	sphere.d[dim][ c.dir] = c.d[dim][c.dir] + dir_sign*0.5*length;
	sphere.d[dim][!c.dir] = sphere.d[dim][c.dir] + dir_sign*0.3*length;
	body  .d[dim][!c.dir] = sphere.d[dim][c.dir] + dir_sign*0.15*length;
	neck  .d[dim][ c.dir] = body.d[dim][!c.dir]; // there will be some intersection, but that should be okay
	neck.expand_in_dim(dim1, -0.29*sz[dim1]); // smaller radius
	neck.expand_in_dim(dim2, -0.29*sz[dim2]); // smaller radius
	cube_t cap(neck);
	neck.d[dim][!c.dir] = cap.d[dim][c.dir] = c.d[dim][!c.dir] - dir_sign*0.08*length; // set cap thickness
	cap.expand_in_dim(dim1, -0.006*sz[dim1]); // slightly larger radius than narrow end of neck
	cap.expand_in_dim(dim2, -0.006*sz[dim2]); // slightly larger radius than narrow end of neck
	// draw as a sphere
	vector3d skip_hemi_dir(zero_vector);
	skip_hemi_dir[dim] = -dir_sign;
	mat.add_sphere_to_verts(sphere, color, 1, skip_hemi_dir); // low_detail=1
	mat.add_ortho_cylin_to_verts(body, color, dim, (add_bottom && !c.dir), (add_bottom && c.dir), 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bottle_ndiv);
	// draw neck of bottle as a truncated cone; draw top if empty
	mat.add_ortho_cylin_to_verts(neck, color, dim, (is_empty && c.dir), (is_empty && !c.dir), 0, 0, (c.dir ? 0.85 : 1.0), (c.dir ? 1.0 : 0.85), 1.0, 1.0, 0, bottle_ndiv);

	if (!is_empty) { // draw cap if nonempty
		bool const draw_bot(c.flags & RO_FLAG_WAS_EXP);
		mat.add_ortho_cylin_to_verts(cap, apply_light_color(c, cap_colors[bool(c.obj_id & 64)]), dim, (draw_bot || c.dir), (draw_bot || !c.dir), 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bottle_ndiv);
	}
	// add the label
	// Note: we could add a bottom sphere to make it a capsule, then translate below the surface in -z to flatten the bottom
	body.expand_in_dim(dim1, 0.03*radius); // expand slightly in radius
	body.expand_in_dim(dim2, 0.03*radius); // expand slightly in radius
	body.d[dim][c.dir] += dir_sign*0.24*length; body.d[dim][!c.dir] -= dir_sign*0.12*length; // shrink in length
	string const &texture_fn(bottle_params[c.get_bottle_type()].texture_fn);
	rgeom_mat_t &label_mat(get_material(tid_nm_pair_t(texture_fn.empty() ? -1 : get_texture_by_name(texture_fn)), 0, 0, 1));
	label_mat.add_ortho_cylin_to_verts(body, apply_light_color(c, WHITE), dim, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bottle_ndiv); // white label
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

void building_room_geom_t::add_pen_pencil_marker_to_material(room_object_t const &c_, rgeom_mat_t &mat) {
	room_object_t c(c_);
	unsigned const dim(get_max_dim(c.get_size()));
	if (dim != 2 && !c.dir) {swap(c.d[dim][0], c.d[dim][1]); c.dir = 1;} // put in canonical orientation if not vertical; okay if denormalized
	colorRGBA const color(apply_light_color(c));
	bool const is_pen(c.type == TYPE_PEN), is_marker(c.type == TYPE_MARKER);
	float const length(c.get_sz_dim(dim));
	cube_t body(c), point(c);
	body.d[dim][1] = point.d[dim][0] = c.d[dim][1] - (is_marker ? 0.25 : (is_pen ? 0.08 : 0.10))*length;
	// non-AA rotation/direction?

	if (is_marker) {
		mat.add_ortho_cylin_to_verts(body,  apply_light_color(c, WHITE), dim, 1, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 16); // 16-sided cylinder, always white
		mat.add_ortho_cylin_to_verts(point, color, dim, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 16); // 16-sided cylinder (cap)
	}
	else if (is_pen) { // point is at the top
		mat.add_ortho_cylin_to_verts(body,  color, dim, 1, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 16); // 16-sided cylinder
		mat.add_ortho_cylin_to_verts(point, color, dim, 0, 0, 0, 0, 1.0, 0.2, 1.0, 1.0, 0, 16); // 16-sided truncated cone
	}
	else { // eraser is at the bottom and point is at the top
		colorRGBA const lt_wood(1.0, 0.8, 0.6);
		cube_t end_part(body), lead(point);
		point.d[dim][1] = lead    .d[dim][0] = c.d[dim][1] - 0.03*length;
		body .d[dim][0] = end_part.d[dim][1] = c.d[dim][0] + 0.09*length;
		cube_t eraser(end_part), metal(end_part);
		metal.d[dim][0] = eraser.d[dim][1] = c.d[dim][0] + 0.04*length;
		mat.add_ortho_cylin_to_verts(body,   color,                          dim, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 6); // 6-sided cylinder, should really be made flat sided
		mat.add_ortho_cylin_to_verts(point,  apply_light_color(c, lt_wood),  dim, 0, 0, 0, 0, 1.0, 0.3, 1.0, 1.0, 0, 12); // 12-sided truncated cone
		mat.add_ortho_cylin_to_verts(lead,   apply_light_color(c, DK_GRAY),  dim, 0, 0, 0, 0, 0.3, 0.0, 1.0, 1.0, 0, 12); // 12-sided cone
		mat.add_ortho_cylin_to_verts(metal,  apply_light_color(c, DK_GREEN), dim, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 12); // 12-sided cylinder
		mat.add_ortho_cylin_to_verts(eraser, apply_light_color(c, PINK),     dim, 1, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 12); // 12-sided cylinder
	}
}
void building_room_geom_t::add_pen_pencil_marker(room_object_t const &c) {
	add_pen_pencil_marker_to_material(c, get_untextured_material(0, 0, 1)); // unshadowed, small
}

void building_room_geom_t::add_flooring(room_object_t const &c, float tscale) {
	get_material(tid_nm_pair_t(MARBLE_TEX, 0.8*tscale)).add_cube_to_verts(c, apply_light_color(c), tex_origin, ~EF_Z2); // top face only, unshadowed
}

void building_room_geom_t::add_wall_trim(room_object_t const &c) {
	rgeom_mat_t &mat(get_untextured_material(0, 0, 1)); // inc_shadows=0, dynamic=0, small=1

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
		mat.add_cube_to_verts_untextured(c, c.color, skip_faces); // is_small, untextured, no shadows
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
	get_untextured_material(1).add_cube_to_verts(c, texture_color(get_blinds_tid()).modulate_with(color), llc, df2); // draw top and bottom / front and back untextured
}

void building_room_geom_t::add_fireplace(room_object_t const &c, float tscale) {
	float const dir_sign(c.dir ? -1.0 : 1.0), depth(c.get_sz_dim(c.dim)), width(c.get_sz_dim(!c.dim)), botz(c.z1() + 0.1*c.dz());
	float const face_pos(c.d[c.dim][!c.dir] - 0.4*dir_sign*depth); // front face pos
	cube_t base(c), front(c), bot(c), top(c);
	base .d[c.dim][!c.dir] = front.d[c.dim][c.dir] = face_pos;
	top  .d[c.dim][!c.dir] = face_pos + dir_sign*0.02*width;
	front.d[c.dim][!c.dir] = face_pos + 0.025*dir_sign*depth; // set front thickness
	front.expand_in_dim(!c.dim, -0.1*width); // shrink
	front.z1()  = bot.z2() = botz;
	front.z2() -= 0.2*c.dz();
	base.z2()   = top.z1() = c.z2() - 0.04*c.dz();
	bot.expand_in_dim(!c.dim, 0.01*width); // expand slightly
	top.expand_in_dim(!c.dim, 0.02*width); // expand slightly
	colorRGBA const color(apply_light_color(c));
	point const tex_origin(c.get_llc());
	int const fp_tid(get_texture_by_name("interiors/fireplace.jpg"));
	unsigned const skip_back_face(~get_face_mask(c.dim, c.dir));
	rgeom_mat_t &brick_mat(get_material(tid_nm_pair_t(BRICK2_TEX, 1.0*tscale, 1), 1)); // shadowed
	brick_mat.add_cube_to_verts(base, color, tex_origin, (skip_back_face | EF_Z1), !c.dim); // skip back and bottom faces
	rgeom_mat_t &front_mat(get_material(tid_nm_pair_t(fp_tid, 0.0, 0), 0)); // unshadowed
	front_mat.add_cube_to_verts(front, color, tex_origin, get_face_mask(c.dim, !c.dir), !c.dim, (c.dim ^ c.dir)); // front face only
	rgeom_mat_t &fside_mat(get_material(tid_nm_pair_t(fp_tid, 0.01, 0), 0)); // unshadowed, small TC so that we get only black area
	fside_mat.add_cube_to_verts(front, color, tex_origin, (get_skip_mask_for_xy(c.dim) | EF_Z1)); // sides only
	rgeom_mat_t &marble_mat(get_material(tid_nm_pair_t(MARBLE_TEX, 1.5*tscale, 1), 1));
	marble_mat.add_cube_to_verts(bot, color, tex_origin, (skip_back_face | EF_Z1)); // skip back and bottom faces
	marble_mat.add_cube_to_verts(top, color, tex_origin,  skip_back_face); // skip back face
}

float get_railing_height(room_object_t const &c) {
	bool const is_u_stairs(c.flags & (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI));
	return (is_u_stairs ? 0.70 : 0.35)*c.dz(); // use a larger relative height for lo/hi railings on U-shaped stairs
}
cylinder_3dw get_railing_cylinder(room_object_t const &c) {
	float const radius(0.5*c.get_sz_dim(!c.dim)), center(c.get_center_dim(!c.dim)), height(get_railing_height(c));
	point p[2];

	for (unsigned d = 0; d < 2; ++d) {
		p[d].z = ((c.flags & RO_FLAG_TOS) ? c.z1() : c.d[2][d]) + height; // top railing is level, otherwise sloped
		p[d][!c.dim] = center;
		p[d][ c.dim] = c.d[c.dim][c.dir^bool(d)^1];
	}
	return cylinder_3dw(p[0], p[1], radius, radius);
}
void building_room_geom_t::add_railing(room_object_t const &c) {
	cylinder_3dw const railing(get_railing_cylinder(c));
	bool const is_u_stairs(c.flags & (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI)), is_top_railing(c.flags & RO_FLAG_TOS), draw_ends(!(c.flags & RO_FLAG_INTERIOR));
	float const pole_radius(0.75*railing.r1), length(c.get_sz_dim(c.dim)), height(get_railing_height(c));
	tid_nm_pair_t tex(-1, 1.0, 1); // shadowed
	tex.set_specular(0.7, 70.0);
	rgeom_mat_t &mat(get_material(tex, 1, 0, 1)); // inc_shadows=1, dynamic=0, small=1
	mat.add_cylin_to_verts(railing.p1, railing.p2, railing.r1, railing.r2, c.color, draw_ends, draw_ends); // draw sloped railing

	if (!is_u_stairs && !(c.flags & RO_FLAG_ADJ_TOP)) {
		for (unsigned d = 0; d < 2; ++d) { // add the two vertical poles
			point pt(d ? railing.p2 : railing.p1);
			if (!is_top_railing) {pt[c.dim] += ((c.dir ^ bool(d)) ? 1.0 : -1.0)*0.01*length;} // shift slightly inward toward the center
			float const hscale((d && !is_top_railing) ? 1.25 : 1.0); // shorten for lower end, which rests on the step (unless top railing)
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

// Note: there is a lot duplicated with building_room_geom_t::add_elevator(), but we need a separate function for adding interior elevator buttons
cube_t get_elevator_car_panel(room_object_t const &c) {
	float const dz(c.dz()), thickness(elevator_fc_thick_scale*dz), signed_thickness((c.dir ? 1.0 : -1.0)*thickness);
	float const width(c.get_sz_dim(!c.dim)), frame_width(0.2*width), panel_width(min(0.9f*frame_width, 0.25f*dz)), front_face(c.d[c.dim][c.dir] - signed_thickness);
	cube_t panel(c);
	panel.d[c.dim][ c.dir] = front_face; // flush front inner wall
	panel.d[c.dim][!c.dir] = front_face - 0.1*signed_thickness; // set panel thickness
	panel.d[!c.dim][0] = c.d[!c.dim][0] + 0.5*(frame_width - panel_width) + thickness; // edge near the wall
	panel.d[!c.dim][1] = panel.d[!c.dim][0] + panel_width - thickness; // edge near door
	panel.z1() += 0.28*dz; panel.z2() -= 0.28*dz;
	return panel;
}
void building_room_geom_t::add_elevator(room_object_t const &c, float tscale) { // elevator car; dynamic=1
	// elevator car, all materials are dynamic; no lighting scale
	float const dz(c.dz()), thickness(elevator_fc_thick_scale*dz), dir_sign(c.dir ? 1.0 : -1.0), signed_thickness(dir_sign*thickness);
	cube_t floor(c), ceil(c), back(c);
	floor.z2() = floor.z1() + thickness;
	ceil. z1() = ceil. z2() - thickness;
	floor.expand_by_xy(-0.5f*thickness);
	ceil .expand_by_xy(-0.5f*thickness);
	back.d[c.dim][c.dir] = c.d[c.dim][!c.dir] + signed_thickness;
	vector3d const tex_origin(c.get_llc());
	unsigned const front_face_mask(get_face_mask(c.dim, c.dir)), floor_ceil_face_mask(front_face_mask & (EF_X12 | EF_Y12)); // +Z faces
	tid_nm_pair_t const paneling(get_tex_auto_nm(PANELING_TEX, 2.0f*tscale));
	get_material(get_tex_auto_nm(TILE_TEX, tscale), 1, 1).add_cube_to_verts(floor, WHITE, tex_origin, floor_ceil_face_mask);
	get_material(get_tex_auto_nm(get_rect_panel_tid(), tscale), 1, 1).add_cube_to_verts(ceil, WHITE, tex_origin, floor_ceil_face_mask);
	rgeom_mat_t &paneling_mat(get_material(paneling, 1, 1));
	paneling_mat.add_cube_to_verts(back, WHITE, tex_origin, front_face_mask, !c.dim);
	float const width(c.get_sz_dim(!c.dim)), frame_width(0.2*width), spacing(0.02*width), front_face(c.d[c.dim][c.dir] - signed_thickness);
	cube_t front(c);
	front.d[c.dim][ c.dir] -= (c.dir ? 1.0 : -1.0)*spacing; // slight gap with elevator doors
	front.d[c.dim][!c.dir]  = front_face;

	for (unsigned d = 0; d < 2; ++d) {
		// side walls
		unsigned const side_skip_faces(get_face_mask(!c.dim, !d));
		cube_t side(c);
		side.d[!c.dim][!d] = c.d[!c.dim][d] + (d ? -1.0 : 1.0)*thickness;
		paneling_mat.add_cube_to_verts(side,  WHITE, tex_origin, side_skip_faces, c.dim);
		// front sides of doors
		front.d[!c.dim][ d] = side.d[!c.dim][!d];
		front.d[!c.dim][!d] = c   .d[!c.dim][ d] + (d ? -1.0 : 1.0)*frame_width;
		paneling_mat.add_cube_to_verts(front, WHITE, tex_origin, (get_face_mask(c.dim, !c.dir) & side_skip_faces), !c.dim); // draw front and inside side
	}
	// add button panel
	cube_t const panel(get_elevator_car_panel(c));
	get_untextured_material(0, 1).add_cube_to_verts(panel, DK_GRAY, all_zeros, ~get_face_mask(c.dim, c.dir));
	// add floor numbers to either the panel (or the buttons themselves?)
	unsigned const num_floors(c.drawer_flags), cur_floor(c.item_flags);
	assert(num_floors > 1);
	float const button_spacing(panel.dz()/(num_floors + 1)); // add extra spacing on bottom and top of panel
	float const panel_width(panel.get_sz_dim(!c.dim)), text_height(min(0.5f*panel_width, 0.8f*button_spacing));
	float const inner_button_radius(min(0.6f*thickness, min(0.35f*button_spacing, 0.25f*panel.get_sz_dim(!c.dim)))); // approx match to elevator
	bool const ldir(c.dim ^ c.dir);
	point text_pos;
	text_pos[ c.dim] = panel.d[c.dim][!c.dir] - 0.01*signed_thickness; // slightly in front of the panel
	text_pos[!c.dim] = max((panel.d[!c.dim][0] + 0.05f*panel_width), (panel.get_center_dim(!c.dim) - (1.5f*inner_button_radius + text_height))) +
		0.6f*ldir*text_height; // shift by approx width of font character(s) because we're aligning to the right side rather than the left
	vector3d col_dir(zero_vector), normal(zero_vector);
	col_dir[!c.dim] = (ldir ? -1.0 : 1.0);
	normal [ c.dim] = -dir_sign; // opposite dir from front of elevator
	static vector<vert_tc_t> verts;
	static ostringstream oss; // reused across buttons
	tid_nm_pair_t tp(FONT_TEXTURE_ID), lit_tp(tp);
	lit_tp.emissive = 1.0;
	get_material(lit_tp, 0, 1); // make sure it's allocated
	rgeom_mat_t &mat(get_material(tp, 0, 1)); // dynamic=1
	color_wrapper const cw(BLACK), lit_cw(colorRGBA(1.0, 0.9, 0.5));
	norm_comp const nc(normal);

	for (unsigned f = 0; f < num_floors; ++f) {
		bool const is_lit(f == cur_floor);
		text_pos.z = panel.z1() + (f + 1)*button_spacing - 0.5*text_height;
		verts.clear();
		oss.str("");
		oss << (f+1); // floor index
		gen_text_verts(verts, text_pos, oss.str(), 1000.0*text_height, col_dir, plus_z, 1); // use_quads=1
		assert(!verts.empty());
		if (dot_product(normal, cross_product((verts[1].v - verts[0].v), (verts[2].v - verts[1].v))) < 0.0) {std::reverse(verts.begin(), verts.end());} // swap vertex winding order
		rgeom_mat_t &cur_mat(is_lit ? get_material(lit_tp, 0, 1) : mat);
		for (auto i = verts.begin(); i != verts.end(); ++i) {cur_mat.quad_verts.emplace_back(i->v, nc, i->t[0], i->t[1], (is_lit ? lit_cw : cw));}
	} // for f
}

void building_room_geom_t::add_elevator_doors(elevator_t const &e) {
	float const spacing(e.get_wall_thickness()), closed_door_width(0.99*0.5*e.get_sz_dim(!e.dim)), open_door_width(1.12*e.get_frame_width());
	rgeom_mat_t &mat(get_untextured_material(1, 1));
	float open_z1(e.z1()), open_z2(e.z2());

	if (e.open_amt > 0.0) { // only draw the doors as open for the floor the elevator car is on
		assert(e.car_obj_id < objs.size());
		room_object_t const &car(objs[e.car_obj_id]); // elevator car for this elevator
		float const z_shrink(0.8*elevator_fc_thick_scale*car.dz()); // shrink slightly to avoid clipping through the ceiling and floor
		open_z1 = car.z1() + z_shrink; open_z2 = car.z2() - z_shrink;
		assert(open_z1 < open_z2);
	}
	for (unsigned d = 0; d < 2; ++d) { // left/right doors, untextured for now
		unsigned const skip_faces(((e.open_amt > 0.0) ? 0 : EF_Z12) | ~get_face_mask(!e.dim, !d)); // skip top and bottom if fully closed
		cube_t door(e);
		door.d[e.dim][!e.dir] = door.d[e.dim][e.dir] + (e.dir ? -1.0f : 1.0f)*spacing; // set correct thickness
		door.expand_in_dim(e.dim, -0.2*spacing); // shrink slightly to make thinner
		door.d[!e.dim][d] = e.d[!e.dim][!d] + (d ? 1.0f : -1.0f)*closed_door_width; // this section is always closed

		if (e.open_amt > 0.0) { // only draw the doors as open for the floor the elevator car is on
			cube_t open_door(door);
			open_door.d[!e.dim][d] += (d ? 1.0f : -1.0f)*(open_door_width - closed_door_width)*e.open_amt;
			open_door.z1() = door.z2() = open_z1; open_door.z2() = open_z2;
			mat.add_cube_to_verts(open_door, GRAY, all_zeros, skip_faces); // open part
			if (door.dz() > 0.0) {mat.add_cube_to_verts(door, GRAY, all_zeros, skip_faces);} // bottom part
			door.z1() = open_z2; door.z2() = e.z2(); // top part
		}
		if (door.dz() > 0.0) {mat.add_cube_to_verts(door, GRAY, all_zeros, skip_faces);} // all or top part
	} // for d
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
	get_untextured_material(); // ensure frame material is valid
	rgeom_mat_t &picture_mat(get_material(tid_nm_pair_t(picture_tid, 0.0)));
	unsigned const picture_qv_start(picture_mat.quad_verts.size());
	picture_mat.add_cube_to_verts(c, c.color, tex_origin, skip_faces, !c.dim, mirror_x);
	// add a frame
	cube_t frame(c);
	vector3d exp;
	exp.z = exp[!c.dim] = (whiteboard ? 0.04 : 0.06)*c.dz(); // frame width
	exp[c.dim] = (whiteboard ? -0.1 : -0.25)*c.get_sz_dim(c.dim); // shrink in this dim
	frame.expand_by(exp);
	rgeom_mat_t &frame_mat(get_untextured_material());
	unsigned const frame_qv_start(frame_mat.quad_verts.size());
	frame_mat.add_cube_to_verts_untextured(frame, (whiteboard ? GRAY : BLACK), skip_faces);
	
	if (whiteboard) { // add a marker ledge
		cube_t ledge(c);
		ledge.z2() = ledge.z1() + 0.016*c.dz(); // along the bottom edge
		ledge.d[c.dim][c.dir] += (c.dir ? 1.5 : -1.5)*c.get_sz_dim(c.dim); // extrude outward
		get_untextured_material(1).add_cube_to_verts_untextured(ledge, GRAY, (1 << (2*(2-c.dim) + !c.dir))); // shadowed
	}
	else if (c.flags & RO_FLAG_RAND_ROT) { // apply a random rotation
		float const angle(0.2*(fract(PI*c.obj_id + 1.61803*c.item_flags) - 0.5)); // random rotation based on obj_id and item flags
		point rotate_pt(c.get_cube_center());
		rotate_pt.z += 0.45*c.dz(); // rotate about a point near the top of the picture
		vector3d normal(zero_vector);
		normal[c.dim] = (c.dir ? -1.0 : 1.0);
		rotate_verts(picture_mat.quad_verts, normal, angle, rotate_pt, picture_qv_start);
		rotate_verts(frame_mat  .quad_verts, normal, angle, rotate_pt, frame_qv_start  );
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

void building_room_geom_t::add_book(room_object_t const &c, bool inc_lg, bool inc_sm, float tilt_angle, unsigned extra_skip_faces, bool no_title, float z_rot_angle) {
	bool const is_held(z_rot_angle != 0.0); // held by the player, and need to draw the bottom
	bool const draw_cover_as_small((c.flags & RO_FLAG_WAS_EXP) || is_held); // books in drawers, held, or dropped are always drawn as small objects
	if (draw_cover_as_small && !inc_sm) return; // nothing to draw
	bool const upright(c.get_sz_dim(!c.dim) < c.dz()); // on a bookshelf
	bool const tdir(upright ? (c.dim ^ c.dir ^ bool(c.obj_id%7)) : 1); // sometimes upside down when upright
	bool const ldir(!tdir), cdir(c.dim ^ c.dir ^ upright ^ ldir); // colum and line directions (left/right/top/bot) + mirror flags for front cover
	bool const was_dropped(c.flags & RO_FLAG_TAKEN1); // or held
	bool const shadowed(was_dropped && !is_held); // only shadowed if dropped by the player, since otherwise shadows are too small to have much effect; skip held objects (don't work)
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
	vector3d axis, tilt_about(c.get_urc()), zrot_about(c.get_cube_center());
	axis[c.dim] = 1.0; // along book width
	tilt_angle *= (c.dim ? -1.0 : 1.0);
	bool has_cover(0);
	colorRGBA const color(apply_light_color(c));
	// skip top face, bottom face if not tilted, thickness dim if upright
	unsigned const sides_mask(upright ? get_skip_mask_for_xy(tdim) : (is_held ? EF_Z12 : EF_Z2)), spine_mask(~get_face_mask(c.dim, !c.dir)); // masks of faces to draw
	unsigned const skip_faces(extra_skip_faces | ((tilt_angle == 0.0) ? EF_Z1 : 0) | sides_mask);

	if (z_rot_angle == 0.0 && (c.flags & RO_FLAG_RAND_ROT) && (c.obj_id%3) == 0) { // books placed on tables/desks are sometimes randomly rotated a bit
		z_rot_angle = (PI/12.0)*(fract(123.456*c.obj_id) - 0.5);
	}
	if (c.is_open()) {
		assert(!upright && !is_held);
		// draw book as open?
	}
	if (draw_cover_as_small || inc_lg) { // draw large faces: outside faces of covers and spine
		rgeom_mat_t &mat(get_untextured_material(shadowed, 0, draw_cover_as_small));
		unsigned const qv_start(mat.quad_verts.size());
		mat.add_cube_to_verts_untextured(c, color, (extra_skip_faces | ~(sides_mask | spine_mask))); // untextured
		rotate_verts(mat.quad_verts, axis,   tilt_angle,  tilt_about, qv_start);
		rotate_verts(mat.quad_verts, plus_z, z_rot_angle, zrot_about, qv_start);
	}
	if (draw_cover_as_small || inc_sm) { // draw small faces: insides of covers, edges, and pages
		rgeom_mat_t &mat(get_untextured_material(shadowed, 0, 1));
		unsigned const qv_start(mat.quad_verts.size());
		mat.add_cube_to_verts_untextured(bot,   color, (extra_skip_faces | (was_dropped ? 0 : (EF_Z1 | ~get_face_mask(tdim, 0))))); // untextured, skip bottom face if not dropped
		mat.add_cube_to_verts_untextured(top,   color, (extra_skip_faces | (upright ? EF_Z1 : 0) | ~get_face_mask(tdim, 1))); // untextured, skip top face, skip bottom face if upright
		mat.add_cube_to_verts_untextured(spine, color, (skip_faces | spine_mask)); // untextured, skip back of spine (drawn as lg geom)
		mat.add_cube_to_verts_untextured(pages, apply_light_color(c, WHITE), (skip_faces | spine_mask)); // untextured
		rotate_verts(mat.quad_verts, axis,   tilt_angle,  tilt_about, qv_start);
		rotate_verts(mat.quad_verts, plus_z, z_rot_angle, zrot_about, qv_start);
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
		cover_mat.add_cube_to_verts(cover, WHITE, zero_vector, get_face_mask(tdim, tdir), swap_xy, ldir, !cdir); // no shadows, small=1
		rotate_verts(cover_mat.quad_verts, axis,   tilt_angle,  tilt_about, qv_start);
		rotate_verts(cover_mat.quad_verts, plus_z, z_rot_angle, zrot_about, qv_start);
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
		rotate_verts(mat.quad_verts, axis,   tilt_angle,  tilt_about, qv_start);
		rotate_verts(mat.quad_verts, plus_z, z_rot_angle, zrot_about, qv_start);
	} // end titles
}

void building_room_geom_t::add_bookcase(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale, bool no_shelves, float sides_scale,
	point const *const use_this_tex_origin, vector<room_object_t> *books)
{
	colorRGBA const color(apply_wood_light_color(c));
	unsigned const skip_faces(c.was_moved() ? 0 : ~get_face_mask(c.dim, !c.dir)); // skip back face, unless moved by the player and no longer against the wall
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
	
	if (inc_lg) {
		unsigned const back_skip_faces(c.was_moved() ? ~get_skip_mask_for_xy(c.dim) : get_face_mask(c.dim, c.dir)); // back - only face oriented outward
		get_wood_material(tscale).add_cube_to_verts(back, color, tex_origin, back_skip_faces);
	}
	if (no_shelves) return;
	// add shelves
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	unsigned const num_shelves(3 + ((17*c.room_id + int(1000.0*fabs(c.z1())))%3)); // 3-5, randomly selected by room ID and floor
	float const shelf_dz(middle.dz()/(num_shelves+0.25)), shelf_thick(0.03*c.dz());
	unsigned const skip_book_flags(c.get_combined_flags());
	cube_t shelves[5];
	
	for (unsigned i = 0; i < num_shelves; ++i) {
		cube_t &shelf(shelves[i]);
		shelf = middle; // copy XY parts
		shelf.z1() += (i+0.25)*shelf_dz;
		shelf.z2()  = shelf.z1() + shelf_thick;
		if (inc_lg) {get_wood_material(tscale).add_cube_to_verts(shelf, color, tex_origin, skip_faces_shelves);} // Note: mat reference may be invalidated by adding books
	}
	// add books
	for (unsigned i = 0, book_ix = 0; i < num_shelves; ++i) {
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
			colorRGBA const &book_color(book_colors[rgen.rand() % NUM_BOOK_COLORS]);
			bool const backwards((rgen.rand()%10) == 0), book_dir(c.dir ^ backwards ^ 1); // spine facing out 90% of the time

			if (!(skip_book_flags & (1<<(book_ix&31)))) { // may have more than 32 books, and will wrap in that case
				assert(book.is_strictly_normalized());
				room_object_t obj(book, TYPE_BOOK, c.room_id, c.dim, book_dir, c.flags, c.light_amt, SHAPE_CUBE, book_color);
				obj.obj_id     = c.obj_id + 123*i + 1367*n;
				obj.item_flags = (uint16_t)book_ix;
				if (inc_lg || inc_sm) {add_book(obj, inc_lg, inc_sm, tilt_angle, skip_faces, backwards);} // detailed book, no title if backwards
				if (books) {books->push_back(obj);}
			}
			++book_ix;
			pos += width;
			last_book_pos = pos;
			prev_tilted   = (tilt_angle != 0.0); // don't tilt two books in a row
		} // for n
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
	vector3d const tex_origin(c.get_llc());

	if (inc_lg) {
		bool const no_mattress(c.flags & RO_FLAG_TAKEN3);
		colorRGBA const color(apply_wood_light_color(c));
		add_tc_legs(legs_bcube, color, max(head_width, foot_width), tscale);
		if (no_mattress) {get_wood_material(4.0*tscale);} // pre-allocate slats material if needed
		rgeom_mat_t &wood_mat(get_wood_material(tscale));

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
			if (c.flags & RO_FLAG_TAKEN2) {get_untextured_material(1).add_cube_to_verts(mattress, sheet_color, tex_origin, mattress_skip_faces);} // sheets taken, bare mattress
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
	rgeom_mat_t &mat(get_untextured_material(1, 0, 1)); // inc_shadows=1, dynamic=0, small=1
	colorRGBA const color(apply_light_color(c));

	if (c.shape == SHAPE_CYLIN) {
		mat.add_vcylin_to_verts(c, color, 1, 0, 1, 1, 0.7, 1.0); // untextured, bottom only, two_sided truncated cone with inverted bottom normal
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

void building_room_geom_t::add_laundry_basket(room_object_t const &c) {
	// Note: no alpha test is enabled in the shader when drawing this, so the holes in the material may not be drawn correctly against objects such as exterior walls
	rgeom_mat_t &tex_mat(get_material(tid_nm_pair_t(get_texture_by_name("interiors/plastic_mesh.png")), 1, 0, 1)); // inc_shadows=1, dynamic=0, small=1
	cube_t bot(c), mid(c), top(c);
	bot.z2() = mid.z1() = c.z1() + 0.12*c.dz();
	mid.z2() = top.z1() = c.z2() - 0.12*c.dz();
	colorRGBA const color(apply_light_color(c));
	tex_mat  .add_vcylin_to_verts(mid, color, 0, 0, 1, 1); // two_sided cylinder
	rgeom_mat_t &solid_mat(get_untextured_material(0, 0, 1)); // inc_shadows=0, dynamic=0, small=1
	solid_mat.add_vcylin_to_verts(bot, color, 1, 0, 1, 1); // two_sided cylinder with bottom
	solid_mat.add_vcylin_to_verts(top, color, 0, 0, 1, 1); // two_sided cylinder
}

void building_room_geom_t::add_br_stall(room_object_t const &c) {
	rgeom_mat_t &mat(get_untextured_material(1));
	colorRGBA const color(apply_light_color(c));
	point const tex_origin(c.get_llc()); // doesn't really need to be set, since stall is untextured

	if (c.shape == SHAPE_SHORT) { // wall separating urinals, drawn as a single cube
		mat.add_cube_to_verts_untextured(c, color, ~get_face_mask(c.dim, c.dir));
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
	mat.add_cube_to_verts_untextured(side1,  color, side_skip_mask);
	mat.add_cube_to_verts_untextured(side2,  color, side_skip_mask);
	mat.add_cube_to_verts_untextured(front1, color, EF_Z12);
	mat.add_cube_to_verts_untextured(front2, color, EF_Z12);

	if (c.is_open()) {
		// draw open door?
	}
	else {
		door.expand_in_dim(!c.dim, -door_gap);
		mat.add_cube_to_verts_untextured(door, color);
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
	rgeom_mat_t &edge_mat(get_untextured_material(0)); // unshadowed
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
		get_untextured_material(0).add_cube_to_verts_untextured(c, WHITE, skip_faces); // back of the sign, always white (for now)
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

bool get_dishwasher_for_ksink(room_object_t const &c, cube_t &dishwasher) {
	if (c.type != TYPE_KSINK) return 0; // error?
	float const dz(c.dz()), depth(c.get_sz_dim(c.dim)), width(c.get_sz_dim(!c.dim)), dir_sign(c.dir ? 1.0 : -1.0);
	if (width <= 3.5*depth) return 0; // too small to fit a dishwasher
	dishwasher = c;
	bool const side((c.flags & RO_FLAG_ADJ_LO) ? 1 : ((c.flags & RO_FLAG_ADJ_HI) ? 0 : (c.obj_id & 1))); // left/right of the sink
	dishwasher.z1() += 0.05*dz;
	dishwasher.z2() -= 0.05*dz;
	dishwasher.d[ c.dim][!c.dir]  = c.d[c.dim][c.dir] - dir_sign*0.1*depth;
	dishwasher.d[ c.dim][ c.dir] += dir_sign*0.05*depth; // front
	dishwasher.d[!c.dim][!side ]  = c.get_center_dim(!c.dim) + (side ? 1.0 : -1.0)*0.66*depth;
	dishwasher.d[!c.dim][ side ]  = dishwasher.d[!c.dim][!side] + (side ? 1.0 : -1.0)*1.05*depth;
	return 1;
}
room_object_t split_cabinet_at_dishwasher(room_object_t &cabinet, cube_t const &dishwasher) { // Note: modifies cabinet
	room_object_t left_part(cabinet);
	left_part.d[!cabinet.dim][1] = dishwasher.d[!cabinet.dim][0];
	cabinet  .d[!cabinet.dim][0] = dishwasher.d[!cabinet.dim][1];
	left_part.flags &= ~RO_FLAG_ADJ_HI;
	cabinet  .flags &= ~RO_FLAG_ADJ_LO;
	cabinet.drawer_flags >>= 8; // left half has first 8 door bits, right half has last 8 door bits
	return left_part;
}
cube_t get_sink_cube(room_object_t const &c) {
	assert(c.type == TYPE_KSINK || c.type == TYPE_BRSINK);
	float const dz(c.dz()), sdepth(0.8*c.get_sz_dim(c.dim)), swidth(min(1.4f*sdepth, 0.75f*c.get_sz_dim(!c.dim)));
	vector3d const center(c.get_cube_center());
	cube_t sink(center, center);
	set_cube_zvals(sink, (c.z2() - 0.3*dz), (c.z2() - 0.05*dz));
	sink.expand_in_dim( c.dim, 0.5*sdepth);
	sink.expand_in_dim(!c.dim, 0.5*swidth);
	return sink;
}

void building_room_geom_t::add_counter(room_object_t const &c, float tscale) { // for kitchens
	float const dz(c.dz()), depth(c.get_sz_dim(c.dim)), dir_sign(c.dir ? 1.0 : -1.0);
	cube_t top(c), dishwasher;
	top.z1() += 0.95*dz;
	tid_nm_pair_t const marble_tex(get_counter_tid(), 2.5*tscale);
	rgeom_mat_t &top_mat(get_material(marble_tex, 1));
	colorRGBA const top_color(apply_light_color(c, WHITE));

	if (c.type == TYPE_KSINK || c.type == TYPE_BRSINK) { // counter with kitchen or bathroom sink
		float const sdepth(0.8*depth);
		vector3d faucet_pos(c.get_cube_center());
		faucet_pos[c.dim] -= dir_sign*0.56*sdepth;
		cube_t const sink(get_sink_cube(c));
		cube_t faucet1(faucet_pos, faucet_pos);
		set_cube_zvals(faucet1, top.z2(), (top.z2() + 0.30*dz));
		faucet1.expand_in_dim( c.dim, 0.04*sdepth);
		faucet1.expand_in_dim(!c.dim, 0.07*sdepth);
		cube_t faucet2(faucet1);
		faucet2.z1()  = faucet1.z2();
		faucet2.z2() += 0.035*dz;
		faucet2.d[c.dim][c.dir] += dir_sign*0.28*sdepth;
		vect_cube_t &cubes(get_temp_cubes());
		subtract_cube_from_cube(top, sink, cubes);
		for (auto i = cubes.begin(); i != cubes.end(); ++i) {top_mat.add_cube_to_verts(*i, top_color, tex_origin);} // should always be 4 cubes
		colorRGBA const sink_color(apply_light_color(c, GRAY));
		rgeom_mat_t &basin_mat(get_metal_material(0));
		basin_mat.add_cube_to_verts(sink, sink_color, tex_origin, EF_Z2, 0, 0, 0, 1); // basin: inverted, skip top face, unshadowed
		
		if (c.drawer_flags > 0) { // draw outside of sink basin if any drawers are open
			cube_t sink_outer(sink);
			sink_outer.expand_by_xy(0.01*dz); // expand by sink basin thickness
			sink_outer.z1() -= 0.1*dz;
			basin_mat.add_cube_to_verts(sink_outer, sink_color, tex_origin, (~get_face_mask(c.dim, !c.dir) | EF_Z2)); // skip back and top
		}
		rgeom_mat_t &metal_mat(get_metal_material(1)); // shadowed, specular metal (specular doesn't do much because it's flat, but may make more of a diff using a cylinder later)
		metal_mat.add_cube_to_verts_untextured(faucet1, sink_color, EF_Z12); // vertical part of faucet, skip top and bottom faces
		metal_mat.add_cube_to_verts_untextured(faucet2, sink_color, 0); // horizontal part of faucet, draw all faces

		if (c.type == TYPE_BRSINK) { // bathroom sink
			metal_mat.add_cube_to_verts_untextured(sink, sink_color, EF_Z2); // outside of basin, no top surface, shadowed
			cube_t front(c);
			front.z2() = top.z1();
			front.z1() = sink.z1() - 0.1*dz; // slightly below the sink basin
			front.d[c.dim][!c.dir] += dir_sign*0.94*depth;
			get_material(marble_tex, 1).add_cube_to_verts(front, top_color, tex_origin, EF_Z2); // front surface, no top face; same as top_mat
		}
		else if (c.type == TYPE_KSINK && get_dishwasher_for_ksink(c, dishwasher)) { // kitchen sink - add dishwasher if wide enough
			unsigned const dw_skip_faces(~get_face_mask(c.dim, !c.dir));
			colorRGBA const dw_color(apply_light_color(c, LT_GRAY));
			metal_mat.add_cube_to_verts_untextured(dishwasher, dw_color, dw_skip_faces);
			cube_t dishwasher_back(dishwasher), handle(dishwasher);
			dishwasher_back.d[c.dim][!c.dir] = c.d[c.dim][!c.dir]; // flush with the cabinet
			metal_mat.add_cube_to_verts_untextured(dishwasher_back, dw_color, ~dw_skip_faces); // draw only the back face, in case it's visible through a window
			handle.z1() += 0.77*dz;
			handle.z2() -= 0.10*dz;
			handle.expand_in_dim(!c.dim, -0.1*depth);
			handle.d[c.dim][ c.dir]  = handle.d[c.dim][!c.dir] = dishwasher.d[c.dim][c.dir];
			handle.d[c.dim][ c.dir] += dir_sign*0.04*depth; // front
			metal_mat.add_ortho_cylin_to_verts(handle, sink_color, !c.dim, 1, 1); // add handle as a cylinder in the proper dim with both ends
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
		if (!dishwasher.is_all_zeros()) {add_cabinet(split_cabinet_at_dishwasher(cabinet, dishwasher), tscale);}
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

float get_cabinet_doors(room_object_t const &c, vect_cube_t &doors) {
	float const cab_depth(c.get_sz_dim(c.dim)), door_height(0.8*c.dz()), door_thick(0.05*door_height);
	cube_t front(c);
	if (c.flags & RO_FLAG_ADJ_LO) {front.d[!c.dim][0] += cab_depth;} // exclude L-joins of cabinets from having doors; assumes all cabinets are the same depth
	if (c.flags & RO_FLAG_ADJ_HI) {front.d[!c.dim][1] -= cab_depth;}
	float const cab_width(front.get_sz_dim(!c.dim));
	if (cab_width < 0.0) return 0.0; // this seems to happen on occasion; maybe it's a bug, or maybe the random size parameters can lead to bad values; either way, skip it
	float door_width(0.75*door_height), door_spacing(1.2*door_width);
	unsigned const num_doors(floor(cab_width/door_spacing));
	if (num_doors == 0) return 0.0; // is this possible?
	assert(num_doors < 1000); // sanity check
	door_spacing = cab_width/num_doors;
	float const tb_border(0.5f*(c.dz() - door_height)), side_border(0.16*door_width), dir_sign(c.dir ? 1.0 : -1.0);
	door_width = (door_spacing - 2.0*side_border); // recalculate actual value
	float lo(front.d[!c.dim][0]);
	cube_t door0(c);
	door0.d[ c.dim][!c.dir]  = door0.d[c.dim][c.dir];
	door0.d[ c.dim][ c.dir] += dir_sign*door_thick; // expand out a bit
	door0.expand_in_dim(2, -tb_border); // shrink in Z

	for (unsigned n = 0; n < num_doors; ++n) {
		cube_t door(door0);
		float const hi(lo + door_spacing);
		door.d[!c.dim][0] = lo;
		door.d[!c.dim][1] = hi;
		door.expand_in_dim(!c.dim, -side_border); // shrink in XY
		doors.push_back(door);
		lo = hi; // advance to next door
	} // for n
	return door_width;
}
void get_cabinet_or_counter_doors(room_object_t const &c, vect_cube_t &doors) {
	doors.clear();
	if (c.type == TYPE_CABINET) {get_cabinet_doors(c, doors); return;}
	room_object_t cabinet(c), dishwasher; // start with counter
	cabinet.z2() -= 0.05*c.dz(); // remove counter top
	cabinet.d[c.dim][c.dir] -= (c.dir ? 1.0 : -1.0)*0.05*c.get_sz_dim(c.dim); // add front overhang
	
	if (c.type == TYPE_KSINK && get_dishwasher_for_ksink(c, dishwasher)) {
		get_cabinet_doors(split_cabinet_at_dishwasher(cabinet, dishwasher), doors);
		doors.resize(8); // hack to force right side cabinet doors to use the correct set of second 8 drawer_flags bits
	}
	get_cabinet_doors(cabinet, doors);
}

void building_room_geom_t::add_cabinet(room_object_t const &c, float tscale) { // for kitchens
	assert(c.is_strictly_normalized());
	static vect_cube_t doors;
	doors.clear();
	float const door_width(get_cabinet_doors(c, doors)), dir_sign(c.dir ? 1.0 : -1.0);
	rgeom_mat_t &wood_mat(get_wood_material(tscale));
	bool const any_doors_open(c.drawer_flags > 0);
	colorRGBA const cabinet_color(apply_wood_light_color(c));
	unsigned skip_faces((c.type == TYPE_COUNTER) ? EF_Z12 : EF_Z2); // skip top face (can't skip back in case it's against a window)

	if (any_doors_open) {
		unsigned const skip_front_face(~get_face_mask(c.dim, c.dir));
		float const wall_thickness(0.04*c.dz());
		vect_cube_t &cubes(get_temp_cubes());
		cubes.push_back(c); // start with entire cabinet
		cube_t interior(c);
		interior.expand_by(-wall_thickness);
		wood_mat.add_cube_to_verts(interior, cabinet_color*0.5, tex_origin, skip_front_face, 0, 0, 0, 1); // darker interior; skip front face; inverted

		for (unsigned n = 0; n < doors.size(); ++n) { // draw open doors as holes
			if (!(c.drawer_flags & (1 << n))) continue; // not open
			cube_t hole(doors[n]), frame(hole);
			hole.expand_in_dim(c.dim, c.get_sz_dim(c.dim)); // expand so that it cuts entirely through the cabinet
			frame.d[c.dim][ c.dir]  = frame.d[c.dim][!c.dir];
			frame.d[c.dim][!c.dir] -= dir_sign*wall_thickness; // move inward by door thickness
			wood_mat.add_cube_to_verts(frame, cabinet_color, tex_origin, get_skip_mask_for_xy(c.dim), 0, 0, 0, 1); // skip front/back face; inverted
			subtract_cube_from_cubes(hole, cubes, nullptr, 1); // clip_in_z=1
		} // for n
		// draw front faces with holes cut in them for open doors
		for (auto i = cubes.begin(); i != cubes.end(); ++i) {wood_mat.add_cube_to_verts(*i, cabinet_color, tex_origin, ~skip_front_face);}
		skip_faces |= skip_front_face; // front face drawn above, don't draw it again below
	}
	wood_mat.add_cube_to_verts(c, cabinet_color, tex_origin, skip_faces); // draw wood exterior

	// add cabinet doors; maybe these should be small objects, but there are at most a few cabinets per house and none in office buildings
	if (doors.empty()) return; // no doors
	get_metal_material(0); // ensure material exists so that door_mat reference is not invalidated
	rgeom_mat_t &door_mat(get_material(get_tex_auto_nm(WOOD2_TEX, 2.0*tscale, any_doors_open), any_doors_open)); // only shadowed if a door is open
	rgeom_mat_t &handle_mat(get_metal_material(0)); // untextured, unshadowed
	colorRGBA const door_color(apply_light_color(c, WHITE)); // lighter color than cabinet
	colorRGBA const handle_color(apply_light_color(c, GRAY_BLACK));
	unsigned const door_skip_faces(~get_face_mask(c.dim, !c.dir));
	float const door_thick(doors[0].get_sz_dim(c.dim)), handle_thick(0.75*door_thick);
	float const hwidth(0.04*doors[0].dz()), near_side(0.1*door_width), far_side(door_width - near_side - hwidth);

	for (unsigned n = 0; n < doors.size(); ++n) {
		bool const is_open(c.drawer_flags & (1 << n)), handle_side(n & 1); // alternate handle side
		cube_t &door(doors[n]);

		if (is_open) { // make this door open
			door.d[ c.dim][c.dir] += dir_sign*(door_width - door_thick); // expand out to full width
			door.d[!c.dim][!handle_side] -= (handle_side ? -1.0 : 1.0)*(door_width - door_thick); // shrink to correct thickness
		}
		door_mat.add_cube_to_verts(door, door_color, tex_origin, door_skip_faces);
		// add door handle
		cube_t handle(door);
		handle.d[c.dim][!c.dir]  = door.d[c.dim][c.dir];
		handle.d[c.dim][ c.dir] += dir_sign*handle_thick; // expand out a bit
		handle.expand_in_dim(2, -0.4*door.dz()); // shrink in Z

		if (is_open) { // rotate 90 degrees
			handle.d[!c.dim][!handle_side] = door.d[!c.dim][handle_side];
			handle.d[!c.dim][ handle_side] = door.d[!c.dim][handle_side] + (handle_side ? 1.0 : -1.0)*handle_thick; // expand out a bit
			handle.d[ c.dim][0] = door.d[c.dim][0] + (!c.dir ? near_side : far_side);
			handle.d[ c.dim][1] = door.d[c.dim][1] - (!c.dir ? far_side : near_side);
		}
		else {
			handle.d[!c.dim][0] = door.d[!c.dim][0] + (handle_side ? near_side : far_side);
			handle.d[!c.dim][1] = door.d[!c.dim][1] - (handle_side ? far_side : near_side);
		}
		handle_mat.add_cube_to_verts_untextured(handle, handle_color, door_skip_faces); // same skip_faces
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

void building_room_geom_t::add_switch(room_object_t const &c) { // light switch, etc.
	unsigned const skip_faces(~get_face_mask(c.dim, c.dir)), front_face_mask(get_face_mask(c.dim, !c.dir)); // skip face that's against the wall
	vector3d const sz(c.get_size());
	cube_t plate(c), rocker(c);
	plate.d[c.dim][!c.dir] -= (c.dir ? -1.0 : 1.0)*0.70*sz[c.dim]; // front face of plate
	set_wall_width(rocker, plate.d[c.dim][!c.dir], 0.15*sz[c.dim], c.dim);
	rocker.expand_in_dim(!c.dim, -0.27*sz[!c.dim]); // shrink horizontally
	rocker.expand_in_dim(2,      -0.39*sz[!c.dim]); // shrink vertically
	rgeom_mat_t &front_mat(get_material(tid_nm_pair_t(get_texture_by_name("interiors/light_switch.jpg"), 0.0, 0), 0, 0, 1));
	front_mat.add_cube_to_verts(plate, c.color, zero_vector, front_face_mask, !c.dim); // textured front face
	rgeom_mat_t &mat(get_untextured_material(0, 0, 1)); // unshadowed, small
	mat.add_cube_to_verts(plate,  c.color, zero_vector, (skip_faces | ~front_face_mask)); // skip front face; always fully lit to match wall
	unsigned const qv_start(mat.quad_verts.size());
	mat.add_cube_to_verts(rocker, c.color, zero_vector, (skip_faces | EF_Z1)); // skip bottom face
	vector3d rot_axis(zero_vector);
	rot_axis[!c.dim] = ((c.dir ^ c.is_open()) ? 1.0 : -1.0);
	rotate_verts(mat.quad_verts, rot_axis, 0.015*PI, plate.get_cube_center(), qv_start); // rotate rocker slightly about base plate center; could be optimized by caching
}

void building_room_geom_t::add_plate(room_object_t const &c) { // is_small=1
	// select plate texture based on room and a property of this building; plates in the same room will match
	unsigned const NUM_PLATE_TEXTURES = 6;
	string const plate_textures[NUM_PLATE_TEXTURES] = {"plates/plate1.png", "plates/plate2.jpg", "plates/plate3.jpg", "plates/plate4.jpg", "plates/plate5.jpg", "plates/plate6.jpg"};
	int const tid(get_texture_by_name(plate_textures[(c.room_id + stairs_start) % NUM_PLATE_TEXTURES]));
	rgeom_mat_t &top_mat(get_material(tid_nm_pair_t(tid, 0.0, 0), 0, 0, 1)); // unshadowed, small
	colorRGBA color(apply_light_color(c));
	UNROLL_3X(min_eq(color[i_], 0.9f);); // clamp color to 90% max to avoid over saturation
	top_mat.add_vcylin_to_verts(c, color, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 1); // top surface, skip sides
	get_untextured_material(0, 0, 1).add_vcylin_to_verts(c, color, 0, 0, 0, 0, 0.8, 1.0); // bottom: truncated cone, untextured, unshadowed, sloped sides only
}

void building_room_geom_t::add_tub_outer(room_object_t const &c) {
	rgeom_mat_t &mat(get_untextured_material(1));
	colorRGBA const color(apply_light_color(c));
	mat.add_cube_to_verts_untextured(c, color, EF_Z12); // shadowed, no top/bottom faces
}

void building_room_geom_t::add_crack(room_object_t const &c) { // in window? (TV and computer monitor cracks are drawn below)
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_crack_tid(c), 0.0, 0), 0, 0, 1)); // unshadowed, small
	mat.add_cube_to_verts(c, apply_light_color(c), all_zeros, get_face_mask(c.dim, c.dir), !c.dim, (c.obj_id&1), (c.obj_id&2)); // X/Y mirror based on obj_id
}

void building_room_geom_t::add_tv_picture(room_object_t const &c) {
	bool const is_broken(c.flags & RO_FLAG_BROKEN);
	if ((c.obj_id & 1) && !is_broken) return; // TV is off half the time
	cube_t screen(c);
	screen.d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*0.35*c.get_sz_dim(c.dim);
	screen.expand_in_dim(!c.dim, -0.03*c.get_sz_dim(!c.dim)); // shrink the sides in
	screen.z1() += 0.09*c.dz();
	screen.z2() -= 0.04*c.dz();
	unsigned skip_faces(get_face_mask(c.dim, c.dir)); // only the face oriented outward

	if (is_broken) {
		rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_crack_tid(c), 0.0)));
		mat.add_cube_to_verts(screen, apply_light_color(c, WHITE), c.get_llc(), skip_faces, !c.dim, (c.obj_id&1), (c.obj_id&2)); // X/Y mirror based on obj_id
	}
	else {
		tid_nm_pair_t tex(((c.shape == SHAPE_SHORT) ? c.get_comp_monitor_tid() : c.get_tv_tid()), 0.0); // computer monitor vs. TV
		tex.emissive = 1.0;
		get_material(tex).add_cube_to_verts(screen, WHITE, c.get_llc(), skip_faces, !c.dim, !(c.dim ^ c.dir));
	}
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
		get_untextured_material(1).add_cylin_to_verts(point(cx, cy, c.z1()), point(cx, cy, pot_top), 0.65*pot_radius, pot_radius, apply_light_color(c), no_dirt, 0, 1, 0);
		
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

xform_matrix get_player_cview_rot_matrix() {
	float const angle(atan2(cview_dir.y, cview_dir.x)); // angle of camera view in XY plane, for rotating about Z
	return get_rotation_matrix(plus_z, angle);
}

void building_room_geom_t::add_lg_ball(room_object_t const &c) { // is_small=1
	bool const dynamic(c.is_dynamic()); // either small or dynamic
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_lg_ball_tid(c), get_lg_ball_nm_tid(c), 0.0, 0.0), 1, dynamic, !dynamic));
	// rotate the texture coords when the ball is rolling
	mat.add_sphere_to_verts(c, apply_light_color(c), 0, zero_vector, (c.has_dstate() ? &get_dstate(c).rot_matrix : nullptr)); // low_detail=0
}
/*static*/ void building_room_geom_t::draw_lg_ball_in_building(room_object_t const &c, shader_t &s) {
	//highres_timer_t timer("Draw Ball"); // 0.105ms
	xform_matrix const rot_matrix(get_player_cview_rot_matrix());
	// Note: since we're using indexed triangles now, we can't simply call draw_quad_verts_as_tris(); instead we create temp VBO/IBO; not the most efficient solution, but it should work
	static rgeom_mat_t mat = rgeom_mat_t(tid_nm_pair_t()); // allocated memory is reused across frames; VBO/IBO are recreated every time
	mat.tex = tid_nm_pair_t(get_lg_ball_tid(c), get_lg_ball_nm_tid(c), 0.0, 0.0);
	mat.add_sphere_to_verts(c, apply_light_color(c), 0, zero_vector, &rot_matrix); // low_detail=0
	mat.upload_draw_and_clear(s);
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
	case TYPE_MONEY:    return texture_color(get_money_tid());
	case TYPE_PHONE:    return color*0.5; // 50% case color, 50% black
	case TYPE_LAPTOP:   return BKGRAY; // black-gray case, ignore logo colors
	case TYPE_TPROLL:   return (WHITE*0.75  + GRAY*0.25);
	case TYPE_SPRAYCAN: return (DK_GRAY*0.5 + color*0.5);
	case TYPE_CRACK:    return ALPHA0; // transparent
	case TYPE_FPLACE:   return texture_color(BRICK2_TEX).modulate_with(color);
	default: return color; // TYPE_LIGHT, TYPE_TCAN, TYPE_BOOK, TYPE_BOTTLE, TYPE_PEN_PENCIL, etc.
	}
	if (is_obj_model_type()) {return color.modulate_with(get_model_color());} // handle models
	return color; // Note: probably should always set color so that we can return it here
}

