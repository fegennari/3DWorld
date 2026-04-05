// 3D World - City Vegetation/Plants
// by Frank Gennari
// 03/31/26

#include "city_objects.h"


void add_cylin_indices_tris(vector<unsigned> &idata, unsigned ndiv, unsigned ix_start); // from animal_draw.cpp

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
	create_and_upload(verts, vector<unsigned>(), 0, 1); // no indices
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
	vector<unsigned> branch_ixs;

	for (unsigned n = 0; n < 4; ++n) { // {+X, -X, +Y, -Y} sides
		if (!(face_mask & (1<<n)))   continue; // face not enabled
		if (rgen.rand_float() < 0.3) continue; // no ivy on this face
		place_on_wall_face(wall, (n>>1), (n&1), leaf_verts, branch_verts, branch_ixs, rgen);
	}
	assert(leaf_verts.empty() == branch_verts.empty());
	if (leaf_verts.empty()) return; // empty
	leaves  .num_verts = leaf_verts  .size();
	branches.num_verts = branch_verts.size();
	branches.num_ixs   = branch_ixs  .size();
	leaves  .create_and_upload(leaf_verts,   vector<unsigned>(), 0, 1); // no indices
	branches.create_and_upload(branch_verts, branch_ixs,         0, 1);
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
	point swap_not_dim_z_to_xy(point const &p) const {
		return point(p[!dim], p.z, p[dim]);
	}
	bool cylins_intersect(cylinder_3dw const &c1, cylinder_3dw const &c2) const { // conservative
		float const r_sum(max(c1.r1, c1.r2) + max(c2.r1, c2.r2)); // works best for fixed radius cylinders
		// line_seg_line_seg_dist_2d() expects lines to be in the same XY plane, so swap dims to make that happen; this won't change the distance
		return (line_seg_line_seg_dist_2d(swap_not_dim_z_to_xy(c1.p1), swap_not_dim_z_to_xy(c1.p2), swap_not_dim_z_to_xy(c2.p1), swap_not_dim_z_to_xy(c2.p2)) < r_sum);
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
			// TODO: should leaves be angled away from the wall to add variety and reduce the chance of Z-fighting?
			leaf.pts[i] = pos + radius*(((bool(i&1) ^ bool(i&2)) ? 1.0 : -1.0)*branch_dir + ((i>>1) ? -1.0 : 1.0)*side_dir);
		}
		cube_t const leaf_bc(leaf.get_bcube());
		if (!check_contained_on_wall(leaf_bc)) return 0; // off the wall
		point const center(leaf_bc.get_cube_center());

		// check for overlaps with leaves previously added for this plant; allow some amount of overlap; shouldn't Z-fight because dist from wall is random
		for (auto i = leaves.begin(); i < leaves.begin()+cur_branch_leaves_start; ++i) {
			cube_t const ibc(i->get_bcube());
			if (center.z > ibc.z1() && center.z < ibc.z2() && center[!dim] > ibc.d[!dim][0] && center[!dim] < ibc.d[!dim][1]) return 0; // too much overlap
		}
		leaves.push_back(leaf);
		return 1;
	}
	void add_leaves_to_branch(cylinder_3dw const &c) {
		float const blen(c.get_length());
		unsigned nleaves(unsigned(1.0 * blen / leaf_sz) + rgen.rand_bool());
		if (leaves.empty()) {max_eq(nleaves, 1U);} // seed branch must have at least one leaf so that leaves is never empty
		if (nleaves == 0) return; // short branch, no leaves
		unsigned const cur_branch_leaves_start(leaves.size());
		vector3d const branch_delta(c.p2 - c.p1), branch_dir(branch_delta.get_norm()), wall_normal(vector_from_dim_dir(dim, dir));
		vector3d const side_dir(cross_product(wall_normal, branch_dir)); // should be normalized
		vector3d const pos_step(branch_delta/(nleaves + 1)); // skip last slot because it may be the tip or leaf loc of next segment
		float const r_step((c.r2 - c.r1)/nleaves);
		point cur_pt(c.p1); // first leaf is at starting point
		float radius(c.r1);

		for (unsigned n = 0; n < nleaves; ++n) {
			float const lsz(rgen.rand_uniform(0.8, 1.2)*leaf_sz); // leaf size +/- 20%
			bool const side((leaves.size() & 1) ^ first_side);
			vector3d const side_dir_leaf((side ? 1.0 : -1.0)*side_dir); // in correct direction
			point leaf_pt(cur_pt);
			leaf_pt += 0.1*rgen.signed_rand_float()*pos_step; // add some random position jitter
			leaf_pt += (radius + 0.5*leaf_sz)*side_dir_leaf; // offset to the side of the branch, alternating sides
			add_leaf(leaf_pt, branch_dir, side_dir_leaf, lsz, cur_branch_leaves_start);
			cur_pt  += pos_step;
			radius  += r_step;
		} // for n
	}
	bool add_branch_seg(point const &p1, point const &p2, float r1, float r2, bool is_new_branch) {
		cylinder_3dw const cand(p1, p2, r1, r2);

		// skip placement check for root cylinder; it will be < wall.z2(), the caller should guarantee it's valid, and it's illegal to have empty cylins
		if (!cylins.empty()) {
			cube_t pt_bc;
			pt_bc.set_from_sphere(p1, r1);
			if (!check_contained_on_wall(pt_bc)) return 0; // p1 off the wall
			pt_bc.set_from_sphere(p2, r2);
			if (!check_contained_on_wall(pt_bc)) return 0; // p2 off the wall

			for (cylinder_3dw const &c : cylins) { // does this check is_new_branch?
				if (p1 == c.p1 || p1 == c.p2) continue; // skip the cylinder we're connected to
				if (cylins_intersect(cand, c)) return 0; // intersects an existing cylinder
			}
		}
		cylins.push_back(cand);
		return 1;
	}
	void end_branch() {
		assert(!cylins.empty());
		cylins.back().r2 = 0.0; // taper the end of the branch
	}
	cylinder_3dw const &select_random_cylin() const {
		assert(!cylins.empty());
		return cylins[rgen.rand() % cylins.size()];
	}
	void add_leaves() {
		assert(!cylins.empty());
		for (cylinder_3dw const &c : cylins) {add_leaves_to_branch(c);}
		assert(!leaves.empty());
	}
	void create_branch_verts(vector<vertex_t> &verts, vector<unsigned> &ixs) const {
		unsigned const ndiv = 8; // or could go up to 12, but 8 seems like enough
		float const ndiv_inv(1.0/ndiv);
		unsigned data_pos(verts.size()), cylin_ix(0);
		point prev_p2;
		assert(!cylins.empty());

		for (cylinder_3dw const &c : cylins) {
			bool const is_join(c.p1 == prev_p2); // this cylinder is joined to the previous cylinder (part of the same branch)
			point const ce[2] = {c.p1, c.p2};
			vector3d v12;
			vector_point_norm const &vpn(gen_cylinder_data(ce, c.r1, c.r2, ndiv, v12));
			if (!is_join) {cylin_ix = 0; data_pos = verts.size();} // reset for new branch

			for (unsigned j = is_join; j < 2; ++j) {
				float const ty(((cylin_ix + j) & 1) ? 1.0 : 0.0); // alternates between texture ends (mirrors); tc_comp is limited to [0.0, 1.0] range

				for (unsigned S = 0; S < ndiv; ++S) {
					float const tx(fabs(S*ndiv_inv - 0.5f)); // [0.0, 1.0]
					vector3d const n(0.5f*(vpn.n[S] + vpn.n[(S+ndiv-1)%ndiv])); // average face normals to get vert normals, don't need to normalize
					verts.emplace_back(vpn.p[(S<<1)+j], n, tx, ty);
				}
			} // for j
			add_cylin_indices_tris(ixs, ndiv, data_pos); // create index data
			data_pos += ndiv;
			prev_p2   = c.p2;
			++cylin_ix;
		} // for c
	}
	void create_leaf_verts(vector<vertex_t> &verts) const {
		float const tcs[4] = {0.0, 1.0, 1.0, 0.0}, tct[4] = {0.0, 0.0, 1.0, 1.0};
		assert(!leaves.empty());

		for (tquad_t const &l : leaves) {
			vector3d const leaf_normal(l.get_norm());
			for (unsigned i = 0; i < 4; ++i) {verts.emplace_back(l.pts[i], leaf_normal, tcs[i], tct[i]);}
		}
	}
}; // end ivy_builder_t

void ivy_wall_t::place_on_wall_face(cube_t const &wall, bool dim, bool dir, vector<vertex_t> &lverts, vector<vertex_t> &bverts, vector<unsigned> &bixs, rand_gen_t &rgen) {
	// generation steps:
	// * select random start points at the base of the wall
	// * create branch with upward direction and random curve
	// * starting at random points along a branch, create a split point that becomes a new branch and widens the branch segments below
	// * place leaves along branches with stem touching the branch
	// * give leaves random orients, but not enough that they clip through walls
	// * check that leaves don't overlap too much with other leaves
	float const leaf_sz(0.05*wall.dz());
	float const wall_len(wall.get_sz_dim(!dim)), wall_edge_space(4.0*leaf_sz), root_spacing(4.0*leaf_sz);
	if (wall_len <= 2.0*wall_edge_space) return; // wall too short; shouldn't happen
	float const pos_lo(wall.d[!dim][0] + wall_edge_space), pos_hi(wall.d[!dim][1] - wall_edge_space), wall_face(wall.d[dim][dir]);
	vector3d const wall_normal(vector_from_dim_dir(dim, dir));
	unsigned const num_roots(4 + (rgen.rand() % 5)); // 4-8
	ivy_builder_t builder(leaf_sz, wall, dim, dir, rgen);
	vector<float> root_vals;

	for (unsigned n = 0; n < num_roots; ++n) {
		// find root pos at base of wall that's not too close to a previous root
		bool root_valid(0);
		point root;

		for (unsigned N = 0; N < 100 && !root_valid; ++N) {
			root[!dim] = rgen.rand_uniform(pos_lo, pos_hi);
			root_valid = 1;

			for (float const &v : root_vals) {
				if (fabs(root[!dim] - v) < root_spacing) {root_valid = 0; break;}
			}
		} // for N
		if (!root_valid) break; // no more roots can be placed
		// determine branch size
		float const branch_radius(rgen.rand_uniform(0.08, 0.12)*leaf_sz), seg_len(12.0*branch_radius);
		root.z    = wall.z1(); // assume this is the ground
		root[dim] = wall_face + branch_radius*(dir ? 1.0 : -1.0); // move slightly away from the surface
		// add branches; main branch in cylinder segments, then connect secondary branches
		unsigned const num_branches(4 + (rgen.rand() % 3)); // 4-6
		point pos(root);
		float radius(branch_radius);
		builder.next_plant();

		for (unsigned B = 0; B < num_branches; ++B) {
			unsigned const num_segs(8 + (rgen.rand() % 5)); // 8-12
			vector3d prev_dir;
			
			if (B > 0) {
				cylinder_3dw const &c(builder.select_random_cylin());
				
				if (c.r2 == 0.0) { // if top is a point/branch end, split at the bottom of the cylinder
					pos    = c.p1; // splits at bottom of cylinder
					radius = c.r1;
				}
				else {
					pos    = c.p2; // splits at top of cylinder
					radius = c.r2; // same radius as cylinder
				}
				// TODO: smaller radius, and closer to wall? does this widen the previous branch?
				prev_dir = (c.p2 - c.p1).get_norm();
			}
			for (unsigned S = 0; S < num_segs; ++S) {
				bool placed(0);

				for (unsigned N = 0; N < 10; ++N) { // 10 tries to place a branch segment
					point pos2(pos);

					if (S == 0) { // new branch
						if (B == 0) { // root; nearly vertical
							pos2.z     += rgen.rand_uniform(1.6, 2.4)*seg_len; // increase height ~2x seg len
							pos2[!dim] += rgen.signed_rand_float()*0.25*seg_len;
						}
						else { // branching point
							vector3d new_dir(prev_dir);
							rotate_vector3d(wall_normal, TO_RADIANS*rgen.rand_uniform(30.0, 60.0)*(rgen.rand_bool() ? 1.0 : -1.0), new_dir); // 30-60 deg from prev branch
							pos2 += (rgen.rand_uniform(0.8, 1.2)*seg_len)*new_dir; // start in the new direction
						}
					}
					else { // continuation
						// TODO: branch should form a curve, not be completely random
						vector3d new_dir(prev_dir);
						rotate_vector3d(wall_normal, 10.0*TO_RADIANS*rgen.signed_rand_float(), new_dir); // +/- 10 deg from prev branch
						new_dir.z = fabs(new_dir.z); // must move upward
						pos2 += (rgen.rand_uniform(0.8, 1.2)*seg_len)*new_dir;
					}
					if (!builder.add_branch_seg(pos, pos2, radius, radius, (S == 0))) continue; // constant radius; reject if branch can't be placed
					prev_dir = (pos2 - pos).get_norm();
					pos      = pos2;
					placed   = 1;
					break; // done
				} // for N
				if (!placed) break; // can't place any more segments on this branch
			} // for S
			builder.end_branch();
		} // for B
		builder.add_leaves();
		builder.create_branch_verts(bverts, bixs);
		builder.create_leaf_verts  (lverts);
	} // for n
	ivy_faces |= (1 << dir); // mark face as having ivy
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
void ivy_manager_t::add_wall(cube_t const &wall, bool dim, unsigned wall_ix, unsigned plot_ix, unsigned city_ix, point const &camera_bs) {
	if (city_ix != cur_city_ix) { // city change
		clear();
		cur_city_ix = city_ix;
	}
	if (((13*plot_ix) % 5) == 0) return; // some plots have no ivy
	if (((17*wall_ix) % 5) == 0) return; // some walls have no ivy
	ivy_wall_t &w(ivy_walls[wall_ix]);

	if (w.leaves.bcube.is_all_zeros()) { // new wall
		unsigned face_mask(dim ? 12 : 3); // either both X sides or both Y sides
		rand_gen_t rgen;
		rgen.set_state(wall_ix+1, plot_ix+1);
		w.gen(wall, face_mask, rgen);
	}
	else { // existing wall
		assert(w.leaves.bcube == wall && w.branches.bcube == wall);
	}
	if (w.empty()) return; // no ivy
	assert(w.ivy_faces);
	if (w.ivy_faces == 1 && camera_bs[dim] > wall.d[dim][1]) return; // facing opposite side
	if (w.ivy_faces == 2 && camera_bs[dim] < wall.d[dim][0]) return; // facing opposite side
	to_draw.push_back(wall_ix); // not checked for duplicates, but there shouldn't be any
}

void drawable_t::draw() const {
	if (num_verts == 0) return;
	pre_render();
	if (num_ixs == 0) {draw_quads_as_tris(num_verts);} // quads
	else {draw_indexed_tri_verts(num_verts, num_ixs, GL_TRIANGLES);} // indexed triangles
}
void ivy_manager_t::draw_and_clear(shader_t &s) {
	if (to_draw.empty()) return;
	// draw branches first
	select_texture(BARK1_TEX);
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


