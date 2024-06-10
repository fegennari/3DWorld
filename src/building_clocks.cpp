// 3D World - Building Analog and Digital Wall Clocks
// by Frank Gennari 11/30/2023

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "scenery.h" // for s_plant
#include "shaders.h"
#include <ctime>

extern int frame_counter;

unsigned get_face_mask(unsigned dim, bool dir);
colorRGBA apply_light_color(room_object_t const &o);


struct clock_time_t {
	unsigned hours=0, mins=0, secs=0;
	int last_update_frame=-1;

	void update(bool use_12_hours=1) {
		if (last_update_frame == frame_counter) return; // already valid, since we're only using second resolution and the frame time should be less than a second
		last_update_frame = frame_counter;
		std::time_t const cur_time(std::time(nullptr)); // get the current time point
		std::tm const cal_time(*std::localtime(std::addressof(cur_time)));
		hours = cal_time.tm_hour;
		mins  = cal_time.tm_min;
		secs  = cal_time.tm_sec;
		if (use_12_hours) {hours = hours % 12;} // convert 24h => 12h
		if (hours == 0  ) {hours = 12;}
	}
	bool operator==(clock_time_t const &c) const {return (secs == c.secs && mins == c.mins && hours == c.hours);}
	void print_time() const {cout << hours << ":" << mins << ":" << secs << endl;}
};

clock_time_t cur_clock_time;

bool check_clock_time() { // returns true if time has moved at least a second since the last call (should be last frame)
	clock_time_t const prev_clock_time(cur_clock_time);
	cur_clock_time.update();
	return !(cur_clock_time == prev_clock_time);
}

/*
  a
f   b
  g
e   c
  d
*/
void add_display_digit(rgeom_mat_t &mat, cube_t const &face, colorRGBA const &on_color, colorRGBA const &off_color, unsigned number, bool dim, bool dir, bool ddir) {
	float const height(face.dz()), swidth((ddir ? -1.0 : 1.0)*face.get_sz_dim(!dim)); // signed width
	float const bar_to_gap(4.0), gap_dz(height/(2*bar_to_gap + 3)), gap_dx(swidth/(bar_to_gap + 2)), bar_dz(bar_to_gap*gap_dz), bar_dx(bar_to_gap*gap_dx);
	float const z2(face.z1() + gap_dz), z3(z2 + bar_dz), z4(z3 + gap_dz), z5(z4 + bar_dz), x2(face.d[!dim][ddir] + gap_dx), x3(x2 + bar_dx);
	cube_t a(face), b(face), c(face), d(face), e(face), f(face), g(face);
	d.z2() = e.z1() = c.z1() = z2;
	e.z2() = c.z2() = g.z1() = z3;
	g.z2() = f.z1() = b.z1() = z4;
	f.z2() = b.z2() = a.z1() = z5;
	e.d[!dim][!ddir] = f.d[!dim][!ddir] = a.d[!dim][ ddir] = d.d[!dim][ ddir] = g.d[!dim][ ddir] = x2;
	b.d[!dim][ ddir] = c.d[!dim][ ddir] = a.d[!dim][!ddir] = d.d[!dim][!ddir] = g.d[!dim][!ddir] = x3;
	cube_t const segs[7] = {a, b, c, d, e, f, g};
	assert(number <= 9);
	unsigned const num_to_segs[10] = {0x3F, 0x06, 0x5B, 0x4F, 0x66, 0x6D, 0x7D, 0x07, 0x7F, 0x6F};
	unsigned const seg_mask(num_to_segs[number]), skip_faces(get_face_mask(dim, dir));

	for (unsigned n = 0; n < 7; ++n) {
		colorRGBA const &color((seg_mask & (1 << n)) ? on_color : off_color);
		mat.add_cube_to_verts_untextured(segs[n], color, skip_faces); // only draw front face
	}
}
void add_display_digit_pair(rgeom_mat_t &mat, cube_t const &face, colorRGBA const &on_color, colorRGBA const &off_color,
	unsigned number, bool dim, bool dir, bool ddir, bool skip_leading_zero)
{
	float const width(face.get_sz_dim(!dim)), space(0.2*width), digit_width(0.5*(width - space));
	cube_t digits[2] = {face, face};
	digits[0].d[!dim][1] = face.d[!dim][0] + digit_width; // left /tens digit
	digits[1].d[!dim][0] = face.d[!dim][1] - digit_width; // right/ones digit
	assert(number < 100); // must fit in 2 digits
	unsigned const nums[2] = {(number / 10), (number % 10)}; // tens, ones
	
	for (unsigned d = 0; d < 2; ++d) {
		bool const skip(skip_leading_zero && d == 0 && nums[d] == 0); // if skipping leading zero, use off color for all digits
		add_display_digit(mat, digits[bool(d) ^ ddir], (skip ? off_color : on_color), off_color, nums[d], dim, dir, ddir);
	}
}
void add_display_colon(rgeom_mat_t &mat, cube_t const &face, colorRGBA const &color, bool dim, bool dir) {
	float const width(face.get_sz_dim(!dim)), shrink(0.35*width), inner_width(width - 2*shrink); // colon occupies the middle 50%
	cube_t face_inner(face);
	face_inner.expand_in_dim(2,    -shrink);
	face_inner.expand_in_dim(!dim, -shrink);
	cube_t dots[2] = {face_inner, face_inner};
	dots[0].z2() = face_inner.z1() + inner_width; // bottom dot
	dots[1].z1() = face_inner.z2() - inner_width; // top    dot
	unsigned const skip_faces(get_face_mask(dim, dir));
	for (unsigned d = 0; d < 2; ++d) {mat.add_cube_to_verts_untextured(dots[d], color, skip_faces);} // only draw front face
}

// hand_pos is the position around the clock from 0.0 at 12:00 to 1.0 clockwise
void add_clock_hand(rgeom_mat_t &mat, point const &center, float length, float stub_len, float width, float hand_pos, bool dim, bool dir) {
	// start with a vertical quad, then rotate into place
	float const bot_hwidth(0.5*width), top_hwidth(0.2*bot_hwidth), angle(TWO_PI*hand_pos);
	point pts[6]; // start at the origin
	pts[0].z = pts[1].z =  -stub_len; // zbot
	pts[2].z = pts[5].z = 0.2*length; // zmid
	pts[3].z = pts[4].z =     length; // ztop
	pts[0][!dim] = pts[5][!dim] = -bot_hwidth;
	pts[1][!dim] = pts[2][!dim] =  bot_hwidth;
	pts[3][!dim] =  top_hwidth;
	pts[4][!dim] = -top_hwidth;
	vector3d rot_axis;
	rot_axis[dim] = (dir ? 1.0 : -1.0);
	rotate_vector3d_multi(rot_axis, angle, pts, 6);
	norm_comp const normal(rot_axis);
	color_wrapper const cw(BLACK);
	float const ts[6] = {0.0, 1.0, 1.0, 1.0, 0.0, 0.0}, tt[6] = {0.0, 0.0, 0.5, 1.0, 1.0, 0.5};
	unsigned const vix(mat.itri_verts.size());
	for (unsigned n = 0; n < 6; ++n) {mat.itri_verts.emplace_back((pts[n] + center), normal, ts[n], tt[n], cw);}
	unsigned ixs[12] = {0,1,2, 0,2,5, 5,2,3, 5,3,4}; // 4 triangles
	if (dim ^ dir ^ 1) {reverse(ixs, ixs+12);} // use correct winding order
	for (unsigned n = 0; n < 12; ++n) {mat.indices.push_back(vix + ixs[n]);}
}

void building_room_geom_t::add_clock(room_object_t const &c, bool add_dynamic) {
	if (add_dynamic) {check_clock_time();} // may be needed for the first frame a clock is visible when the seconds haven't yet changed

	if (c.item_flags & 1) { // digital clock
		if (add_dynamic) {
			colorRGBA const on_color(RED), off_color(on_color*0.05);
			cube_t face(c);
			face.expand_in_dim(!c.dim, -0.1*c.get_width ());
			face.expand_in_dim(2,      -0.1*c.get_height());
			face.d[c.dim][0] = face.d[c.dim][1] = c.d[c.dim][c.dir] + (c.dir ? 1.0 : -1.0)*0.1*c.get_depth(); // move out a bit from the front
			tid_nm_pair_t tp; // unshadowed
			tp.emissive = 1.0;
			rgeom_mat_t &mat(get_material(tp, 0, 1)); // unshadowed, dynamic
			// format is HH:MM:SS ; colon_sz includes the spaces between digits and colons
			bool const ddir(c.dim ^ c.dir ^ 1);
			float const digit_to_colons_sz(3.5), colon_sz((ddir ? -1.0 : 1.0)*face.get_sz_dim(!c.dim)/(3*digit_to_colons_sz + 2)), digit_sz(digit_to_colons_sz*colon_sz);
			float const he(face.d[!c.dim][ddir] + digit_sz), ms(he + colon_sz), me(ms + digit_sz), ss(me + colon_sz);
			cube_t dh(face), dm(face), ds(face), chm(face), cms(face); // digits, colons
			dh .d[!c.dim][!ddir] = chm.d[!c.dim][ddir] = he;
			chm.d[!c.dim][!ddir] = dm .d[!c.dim][ddir] = ms;
			dm .d[!c.dim][!ddir] = cms.d[!c.dim][ddir] = me;
			cms.d[!c.dim][!ddir] = ds .d[!c.dim][ddir] = ss;
			bool const show_colons(cur_clock_time.secs & 1); // alternate each second
			colorRGBA const &colon_color(show_colons ? on_color : off_color);
			add_display_digit_pair(mat, dh, on_color, off_color, cur_clock_time.hours, c.dim, c.dir, ddir, 1);
			add_display_colon(mat, chm, colon_color, c.dim, c.dir);
			add_display_digit_pair(mat, dm, on_color, off_color, cur_clock_time.mins,  c.dim, c.dir, ddir, 0);
			add_display_colon(mat, cms, colon_color, c.dim, c.dir);
			add_display_digit_pair(mat, ds, on_color, off_color, cur_clock_time.secs,  c.dim, c.dir, ddir, 0);
		}
		else {
			get_untextured_material(1, 0, 1).add_cube_to_verts_untextured(c, apply_light_color(c), ~get_face_mask(c.dim, !c.dir)); // shadowed, small; skip back face
		}
	}
	else { // analog clock
		float const radius(0.5*c.dz());
		point center(c.get_cube_center());
		center[c.dim] = c.d[c.dim][c.dir];

		if (add_dynamic) {
			float const second_pos(cur_clock_time.secs/60.0), minute_pos((cur_clock_time.mins + second_pos)/60.0), hour_pos((cur_clock_time.hours + minute_pos)/12.0); // [0.0, 1.0]
			float const step_dist((c.dir ? 1.0 : -1.0)*0.1*c.get_depth());
			rgeom_mat_t& mat(get_untextured_material(0, 1)); // unshadowed, dynamic
			center[c.dim] += step_dist;
			add_clock_hand(mat, center, 0.56*radius, 0.1*radius, 0.06*radius, hour_pos,   c.dim, c.dir); // hour   hand
			center[c.dim] += step_dist;
			add_clock_hand(mat, center, 0.86*radius, 0.1*radius, 0.06*radius, minute_pos, c.dim, c.dir); // minute hand
			center[c.dim] += step_dist;
			add_clock_hand(mat, center, 0.88*radius, 0.2*radius, 0.03*radius, second_pos, c.dim, c.dir); // second hand
		}
		else {
			get_untextured_material(1, 0, 1).add_ortho_cylin_to_verts(c, apply_light_color(c), c.dim, 0, 0); // shadowed, small; draw sides only
			rgeom_mat_t& face_mat(get_material(tid_nm_pair_t(get_texture_by_name("interiors/clock_face.png")), 1, 0, 1)); // shadows, small
			vector3d const face_dir(vector_from_dim_dir(c.dim, c.dir));
			bool const swap_txy(1), inv_ts(c.dir ^ c.dim), inv_tt(c.dir ^ c.dim ^ 1);
			face_mat.add_disk_to_verts(center, radius, face_dir, WHITE, swap_txy, inv_ts, inv_tt); // always white
		}
	}
}

void building_t::add_clock(cube_t const &clock, unsigned room_id, float tot_light_amt, bool dim, bool dir, bool digital) {
	assert(has_room_geom());
	room_object_t obj(clock, TYPE_CLOCK, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, (digital ? SHAPE_CUBE : SHAPE_CYLIN), (digital ? BLACK : WHITE));
	if (digital) {obj.item_flags = 1;}
	interior->room_geom->objs.push_back(obj);
	interior->room_geom->have_clock = 1; // flag so that we know to update the draw state
}


