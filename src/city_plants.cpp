// 3D World - City Vegetation/Plants
// by Frank Gennari
// 03/31/26

#include "city_objects.h"


// hedge_draw_t

void add_leaf_verts(point const &pos, vector3d const &normal, float angle, float leaf_sz, vector<vert_norm_comp_tc_comp> &verts) {
	vector3d tangent;
	rotate_vector3d(cross_product(normal, ((fabs(normal.x) < fabs(normal.y)) ? plus_x : plus_y)), normal, angle, tangent);
	vector3d const binormal(cross_product(normal, tangent)), dx(leaf_sz*tangent), dy(leaf_sz*binormal);
	point const pts[4] = {(pos - dx - dy), (pos + dx - dy), (pos + dx + dy), (pos - dx + dy)};
	float const ts [4] = {0.0, 1.0, 1.0, 0.0}, tt[4] = {0.0, 0.0, 1.0, 1.0};
	for (unsigned i = 0; i < 4; ++i) {verts.emplace_back(pts[i], normal, ts[i], tt[i]);}
}

void hedge_draw_t::create(cube_t const &bc) {
	bcube = bc - bc.get_cube_center(); // centered on the origin
	unsigned const target_num_leaves(40000);
	vector3d const sz(bcube.get_size());
	float const leaf_sz(0.05*sz.z), surf_area(sz.x*sz.y + 2.0f*sz.z*(sz.x + sz.y));
	float const side_areas[5] = {sz.y*sz.z, sz.y*sz.z, sz.x*sz.z, sz.x*sz.z, sz.x*sz.y};
	rand_gen_t rgen;
	vector<vert_norm_comp_tc_comp> verts;
	verts.reserve(4*target_num_leaves);

	for (unsigned n = 0; n < 5; ++n) { // {+X, -X, +Y, -Y, +Z} sides
		unsigned const dim(n>>1), dir(n&1), d1((dim+1)%3), d2((dim+2)%3);
		unsigned const num_this_face(target_num_leaves*side_areas[n]/surf_area);
		point pos;
		pos[dim] = bcube.d[dim][!dir];

		for (unsigned n = 0; n < num_this_face; ++n) {
			pos[d1] = rgen.rand_uniform(bcube.d[d1][0], bcube.d[d1][1]);
			pos[d2] = rgen.rand_uniform(bcube.d[d2][0], bcube.d[d2][1]);
			vector3d const normal(rgen.signed_rand_vector_spherical_norm());
			float const angle(TWO_PI*rgen.rand_float());
			add_leaf_verts(pos, normal, angle, leaf_sz, verts);
		}
	} // for s
	num_verts = verts.size();
	create_and_upload(verts, 0, 1);
}

void begin_leaf_draw(shader_t &s, int tid) {
	select_texture(tid);
	bind_default_flat_normal_map();
	enable_blend(); // slightly smoother, but a bit of background shows through
	s.add_uniform_float("min_alpha", 0.9);
	s.add_uniform_int("two_sided_lighting", 1);
	s.set_cur_color(WHITE);
}
void end_leaf_draw(shader_t &s) {
	drawable_t::post_render();
	s.add_uniform_int("two_sided_lighting", 0); // reset
	s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA); // restore to the default
	disable_blend();
}

void hedge_draw_t::draw_and_clear(shader_t &s) {
	if (empty()) return;
	if (!vbo_valid()) {create(to_draw.front());}
	begin_leaf_draw(s, get_texture_by_name("pine2.jpg"));
	pre_render();
	vector3d const sz_mult(bcube.get_size().inverse());
	// we can almost use an instance_render here, but that requires changes to the shader or a custom shader

	for (auto c = to_draw.begin(); c != to_draw.end(); ++c) {
		bool const swap_dims((c->dx() < c->dy()) ^ (bcube.dx() < bcube.dy())); // wrong dim, rotate 90 degrees
		vector3d sz(c->get_size());
		if (swap_dims) {swap(sz.x, sz.y);}
		fgPushMatrix();
		translate_to(c->get_cube_center()); // align the center
		if (swap_dims) {fgRotate(90.0, 0.0, 0.0, 1.0);}
		scale_by(sz_mult*sz); // scale to match the size
		s.upload_mvm();
		draw_quads_as_tris(num_verts);
		fgPopMatrix();
	} // for c
	end_leaf_draw(s);
	to_draw.clear();
}

// ivy_manager_t

void ivy_manager_t::ivy_wall_t::gen(cube_t const &bcube, unsigned face_mask, rand_gen_t &rgen) {
	leaves.bcube = branches.bcube = bcube; // both the same
	float const leaf_sz(0.05*bcube.dz()), branch_radius(0.1*leaf_sz);
	vector<vert_norm_comp_tc_comp> leaf_verts, branch_verts;

	for (unsigned n = 0; n < 4; ++n) { // {+X, -X, +Y, -Y} sides
		if (!(face_mask & (1<<n)))   continue; // face not enabled
		if (rgen.rand_float() < 0.3) continue; // no ivy on this face
		unsigned const dim(n>>1), dir(n&1), d1(1-dim), d2(2);
		unsigned const num(50 + (rgen.rand() % 50)); // temp placeholder
		float const wall_face(bcube.d[dim][dir]), shift_dist(0.1*leaf_sz*(dir ? 1.0 : -1.0));
		point pos;
		pos[dim] = wall_face + shift_dist; // move slightly away from the surface
		vector3d const normal(vector_from_dim_dir(dim, dir));
		vert_norm_comp_tc_comp bv;
		bv.set_norm(normal);
		bv.v[dim] = wall_face + 0.5*shift_dist; // in front of the wall but behind the leaf

		for (unsigned n = 0; n < num; ++n) { // TODO: handle Z-fighting of overlapping leaves
			pos[d1] = rgen.rand_uniform(bcube.d[d1][0], bcube.d[d1][1]);
			pos[d2] = rgen.rand_uniform(bcube.d[d2][0], bcube.d[d2][1]);
			float const angle(TWO_PI*rgen.rand_float());
			add_leaf_verts(pos, normal, angle, leaf_sz, leaf_verts);
			// add placeholder vertical branch quad
			float const p1a(pos[d1] - branch_radius), p1b(pos[d1] + branch_radius);
			float const p2a(bcube.z1()), p2b(pos[d2] - 0.5*leaf_sz), tscale1(0.5), tscale2(1.0); // Note: tscale must be in [0.0, 1.0]
			float const vd1[4] = {p1a, p1b, p1b, p1a}, vd2[4] = {p2a, p2a, p2b, p2b};
			float const tcs[4] = {0.0, 1.0, 1.0, 0.0}, tct[4] = {0.0, 0.0, 1.0, 1.0};
			
			for (unsigned i = 0; i < 4; ++i) {
				bv.v[d1] = vd1[i];
				bv.v[d2] = vd2[i];
				bv.set_tcs(tscale1*tcs[i], tscale2*tct[i]);
				branch_verts.push_back(bv);
			}
		} // for n
	} // for s
	assert(leaf_verts.empty() == branch_verts.empty());
	if (leaf_verts.empty()) return; // empty
	leaves  .num_verts = leaf_verts  .size();
	branches.num_verts = branch_verts.size();
	leaves  .create_and_upload(leaf_verts,   0, 1);
	branches.create_and_upload(branch_verts, 0, 1);
}

size_t ivy_manager_t::get_gpu_mem() const {
	size_t mem(0);
	for (auto const &kv : ivy_walls) {mem += kv.second.get_gpu_mem();}
	return mem;
}
void ivy_manager_t::clear() {
	//assert(to_draw.empty()); // can't clear mid-draw; but maybe this can happen if two cities are placed close together and the player is between them, so we allow it
	to_draw.clear();
	for (auto &kv : ivy_walls) {kv.second.clear();}
	ivy_walls.clear();
}
void ivy_manager_t::add_wall(cube_t const &wall, unsigned face_mask, unsigned wall_ix, unsigned plot_ix, unsigned city_ix) {
	if (city_ix != cur_city_ix) { // city change
		clear();
		cur_city_ix = city_ix;
	}
	if (((13*plot_ix) % 5) == 0) return; // some plots have no ivy
	if (((17*wall_ix) % 5) == 0) return; // some walls have no ivy
	ivy_wall_t &w(ivy_walls[wall_ix]);

	if (w.leaves.bcube.is_all_zeros()) { // new wall
		rand_gen_t rgen;
		rgen.set_state(wall_ix+1, plot_ix+1);
		w.gen(wall, face_mask, rgen);
	}
	else { // existing wall
		assert(w.leaves.bcube == wall && w.branches.bcube == wall);
	}
	if (!w.empty()) {to_draw.push_back(wall_ix);} // not checked for duplicates, but there shouldn't be any
}

void drawable_t::draw() const {
	if (num_verts == 0) return;
	pre_render();
	draw_quads_as_tris(num_verts);
}
void ivy_manager_t::draw_and_clear(shader_t &s) {
	if (to_draw.empty()) return;
	//cout << TXT(to_draw.size()) << TXT(ivy_walls.size()) << endl;
	// draw branches first
	select_texture(BARK1_TEX); // TODO: another (custom) texture?
	bind_default_flat_normal_map();
	s.set_cur_color(LT_BROWN); // darker than the texture
	
	for (uint32_t wix : to_draw) {
		auto it(ivy_walls.find(wix));
		assert(it != ivy_walls.end()); // must be found
		it->second.branches.draw();
	}
	// draw leaves second
	begin_leaf_draw(s, PLANT1_TEX);
	for (uint32_t wix : to_draw) {ivy_walls.find(wix)->second.leaves.draw();}
	end_leaf_draw(s);
	to_draw.clear();
}


