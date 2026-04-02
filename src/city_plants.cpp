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

void ivy_wall_t::gen(cube_t const &wall, unsigned face_mask, rand_gen_t &rgen) {
	leaves.bcube = branches.bcube = wall; // both the same
	vector<vertex_t> leaf_verts, branch_verts;

	for (unsigned n = 0; n < 4; ++n) { // {+X, -X, +Y, -Y} sides
		if (!(face_mask & (1<<n)))   continue; // face not enabled
		if (rgen.rand_float() < 0.3) continue; // no ivy on this face
		place_on_wall_face(wall, (n>>1), (n&1), leaf_verts, branch_verts, rgen);
	}
	assert(leaf_verts.empty() == branch_verts.empty());
	if (leaf_verts.empty()) return; // empty
	leaves  .num_verts = leaf_verts  .size();
	branches.num_verts = branch_verts.size();
	leaves  .create_and_upload(leaf_verts,   0, 1);
	branches.create_and_upload(branch_verts, 0, 1);
}

class ivy_builder_t {
	typedef ivy_wall_t::vertex_t vertex_t;
	bool dim, dir, first_side=0;
	float leaf_sz;
	cube_t const &wall;
	rand_gen_t &rgen;
	vector<cylinder_3dw> cylins; // for branches
	vector<tquad_t> leaves; // 4 points

	bool check_contained_on_wall(cube_t const &c) const {
		if (c.d[!dim][0] < wall.d[!dim][0] || c.d[!dim][1] > wall.d[!dim][1]) return 0; // off the wall horizontally
		if (c.z1() < wall.z1() || c.z2() > wall.z2()) return 0; // off the wall vertically
		return 1;
	}
public:
	ivy_builder_t(float leaf_sz_, cube_t const &wall_, bool dim_, bool dir_, rand_gen_t &rgen_) :
		dim(dim_), dir(dir_), leaf_sz(leaf_sz_), wall(wall_), rgen(rgen_) {}

	void next_plant() { // aka clear()
		cylins.clear();
		leaves.clear();
		first_side = rgen.rand_bool();
	}
	bool add_leaf(point const &pos, vector3d const &branch_dir, vector3d const &side_dir, float lsz, unsigned cur_branch_leaves_start) {
		float const radius(0.5*lsz);
		tquad_t leaf(4);

		for (unsigned i = 0; i < 4; ++i) {
			// should leaves be angled away from the wall to add variety and reduce the chance of Z-fighting?
			leaf.pts[i] = pos + radius*(((i>>1) ? 1.0 : -1.0)*branch_dir + ((bool(i&1) ^ bool(i&2)) ? 1.0 : -1.0)*side_dir);
		}
		cube_t const leaf_bc(leaf.get_bcube());
		if (!check_contained_on_wall(leaf_bc)) return 0; // off the wall

		for (auto i = leaves.begin(); i < leaves.begin()+cur_branch_leaves_start; ++i) {
			if (i->get_bcube().intersects(leaf_bc)) return 0; // intersects a leaf from a different branch; too restrictive?
		}
		leaves.push_back(leaf);
		return 1;
	}
	void add_leaves_to_branch(cylinder_3dw const &c) {
		float const blen(c.get_length());
		unsigned const nleaves(unsigned(0.75 * blen / leaf_sz) + rgen.rand_bool());
		if (nleaves == 0) return; // short branch, no leaves
		unsigned const cur_branch_leaves_start(leaves.size());
		vector3d const branch_delta(c.p2 - c.p1), branch_dir(branch_delta.get_norm()), wall_normal(vector_from_dim_dir(dim, dir));
		vector3d const side_dir(cross_product(wall_normal, branch_dir)); // should be normalized
		vector3d const pos_step(c.p1 + branch_delta/nleaves);
		float const r_step((c.r2 - c.r1)/nleaves);
		point cur_pt(c.p1); // first leaf is at starting point
		float radius(c.r1);

		for (unsigned n = 0; n < nleaves; ++n) {
			float const lsz(rgen.rand_uniform(0.8, 1.2)*leaf_sz); // leaf size +/- 20%
			bool const side((leaves.size() & 1) ^ first_side);
			point leaf_pt(cur_pt);
			leaf_pt += (radius + 0.5*leaf_sz)*(side ? 1.0 : -1.0)*side_dir; // offset to the side of the branch, alternating sides
			add_leaf(leaf_pt, branch_dir, side_dir, lsz, cur_branch_leaves_start);
			cur_pt  += pos_step;
			radius  += r_step;
		} // for n
	}
	bool add_branch(point const &p1, point const &p2, float r1, float r2) {
		cube_t pt_bc;
		pt_bc.set_from_sphere(p1, r1);
		if (!check_contained_on_wall(pt_bc)) return 0; // p1 off the wall
		pt_bc.set_from_sphere(p2, r2);
		if (!check_contained_on_wall(pt_bc)) return 0; // p2 off the wall

		for (cylinder_3dw const &c : cylins) {
			if (p1 == c.p1 || p1 == c.p2) continue; // skip the cylinder we're connected to
			// TODO: cylinder-cylinder intersection test
		}
		cylins.emplace_back(p1, p2, r1, r2);
		return 1;
	}
	void add_leaves() {
		for (cylinder_3dw const &c : cylins) {add_leaves_to_branch(c);}
	}
	void create_branch_verts(vector<vertex_t> &verts) const {
		unsigned const ndiv = 12; // should this vary by branch radius?
		point prev_p2;

		for (cylinder_3dw const &c : cylins) {
			bool const is_join(c.p1 == prev_p2); // this cylinder is joined to the previous cylinder (part of the same branch)
			// TODO: see tree and snake code that joins cylinders
			prev_p2 = c.p2;
		} // for c
	}
	void create_leaf_verts(vector<vertex_t> &verts) const {
		float const tcs[4] = {0.0, 1.0, 1.0, 0.0}, tct[4] = {0.0, 0.0, 1.0, 1.0};

		for (tquad_t const &l : leaves) {
			vector3d const leaf_normal(l.get_norm());
			for (unsigned i = 0; i < 4; ++i) {verts.emplace_back(l.pts[i], leaf_normal, tcs[i], tct[i]);}
		}
	}
}; // end ivy_builder_t

void ivy_wall_t::place_on_wall_face(cube_t const &wall, bool dim, bool dir, vector<vertex_t> &lverts, vector<vertex_t> &bverts, rand_gen_t &rgen) {
	// generation steps:
	// * select random start points at the base of the wall
	// * create branch with upward direction and random curve
	// * clamp to or stop at wall bounds
	// * starting at random points along a branch, create a split point that becomes a new branch and widens the branch segments below
	// * place leaves along branches with stem touching the branch
	// * give leaves random orients, but not enough that they clip through walls
	// * check that branches don't intersect
	// * check that leaves don't intersect -or- handle Z-fighting of overlapping leaves
#if 1
	unsigned const num_roots(4 + (rgen.rand() % 5)); // 4-8
	float const leaf_sz(0.05*wall.dz()), branch_radius(0.1*leaf_sz);
	ivy_builder_t builder(leaf_sz, wall, dim, dir, rgen);

	for (unsigned n = 0; n < num_roots; ++n) {
		builder.next_plant();
		// TODO: add branches
		builder.add_leaves();
		builder.create_branch_verts(bverts);
		builder.create_leaf_verts  (lverts);
	} // for n
#else
	unsigned const d1(1-dim), d2(2), num(50 + (rgen.rand() % 50));
	float const wall_face(wall.d[dim][dir]), shift_dist(0.1*leaf_sz*(dir ? 1.0 : -1.0));
	point pos;
	pos[dim] = wall_face + shift_dist; // move slightly away from the surface
	vector3d const normal(vector_from_dim_dir(dim, dir));
	vertex_t bv;
	bv.set_norm(normal);
	bv.v[dim] = wall_face + 0.5*shift_dist; // in front of the wall but behind the leaf

	for (unsigned n = 0; n < num; ++n) {
		pos[d1] = rgen.rand_uniform(wall.d[d1][0], wall.d[d1][1]);
		pos[d2] = rgen.rand_uniform(wall.d[d2][0], wall.d[d2][1]);
		float const angle(TWO_PI*rgen.rand_float());
		add_leaf_verts(pos, normal, angle, leaf_sz, lverts);
		// add placeholder vertical branch quad
		float const p1a(pos[d1] - branch_radius), p1b(pos[d1] + branch_radius);
		float const p2a(wall.z1()), p2b(pos[d2] - 0.5*leaf_sz), tscale1(0.5), tscale2(1.0); // Note: tscale must be in [0.0, 1.0]
		float const vd1[4] = {p1a, p1b, p1b, p1a}, vd2[4] = {p2a, p2a, p2b, p2b};
		float const tcs[4] = {0.0, 1.0, 1.0, 0.0}, tct[4] = {0.0, 0.0, 1.0, 1.0};

		for (unsigned i = 0; i < 4; ++i) {
			bv.v[d1] = vd1[i];
			bv.v[d2] = vd2[i];
			bv.set_tcs(tscale1*tcs[i], tscale2*tct[i]);
			bverts.push_back(bv);
		}
	} // for n
#endif
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


