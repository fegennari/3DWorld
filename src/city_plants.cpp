// 3D World - City Vegetation/Plants
// by Frank Gennari
// 03/31/26

#include "city_objects.h"
#include "subdiv.h" // for sd_sphere_d


void add_cylin_indices_tris(vector<unsigned> &idata, unsigned ndiv, unsigned ix_start); // from animal_draw.cpp
void rotate_verts(point *verts, unsigned num_verts, vector3d const &axis, float angle, vector3d const &about);

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

void ivy_wall_t::gen(cube_t const &wall, vect_cube_t const &avoid, float leaf_sz, unsigned face_mask, bool rand_select, bool is_house, rand_gen_t &rgen) {
	//highres_timer_t timer("Gen Ivy"); // 127 147.263 3.2178 1.15955
	// maintained_amt really should be per-yard, not per-wall, but the mapping from yards to walls isn't tracked, and both faces of a wall and processed together
	maintained_amt = CLIP_TO_01(rgen.rand_uniform(-1.0, 3.0)); // more often maintained
	leaves.bcube   = branches.bcube = wall; // both the same
	vector<vertex_t> leaf_verts, branch_verts;
	vector<unsigned> branch_ixs;

	for (unsigned n = 0; n < 4; ++n) { // {-X, +X, -Y, +Y} sides
		if (!(face_mask & (1<<n))) continue; // face not enabled
		if (rand_select && rgen.rand_float() < 0.3) continue; // no ivy on this face
		place_on_wall_face(wall, avoid, (n>>1), (n&1), is_house, leaf_sz, leaf_verts, branch_verts, branch_ixs, rgen);
	}
	assert(leaf_verts.empty() == branch_verts.empty());
	if (leaf_verts.empty()) return; // empty
	leaves  .num_verts = leaf_verts  .size();
	branches.num_verts = branch_verts.size();
	branches.num_ixs   = branch_ixs  .size();
	leaves  .create_and_upload(leaf_verts,   vector<unsigned>(), 0, 1); // no indices
	branches.create_and_upload(branch_verts, branch_ixs,         0, 1);
	indexed_vao_manager_t::post_render();
}

class ivy_builder_t {
	typedef ivy_wall_t::vertex_t vertex_t;
	bool dim, dir, first_side=0, has_horizontal=0, is_house=0;
	unsigned main_cylins_end=0;
	float leaf_sz;
	cube_t const &wall;
	rand_gen_t &rgen;
	vector<cylinder_3dw> cylins; // for branches
	vector<tquad_t> leaves; // 4 points
	vector<sphere_t> branch_bends; // vertical <=> horizontal transitions
	sphere_point_norm spn;
	vector<unsigned> sphere_ixs;
	vector<vert_norm_tc> sphere_verts;

	bool check_contained_on_wall_xy(cube_t const &c) const {
		return (c.d[!dim][0] > wall.d[!dim][0] && c.d[!dim][1] < wall.d[!dim][1]); // check if off the wall horizontally
	}
	bool check_contained_on_wall(cube_t const &c) const {
		if (!check_contained_on_wall_xy(c)) return 0; // off the wall horizontally
		if (!is_house && !has_horizontal && c.z2() > wall.z2()) return 0; // above the wall; allow if we already have a horizontal branch on the wall
		if (c.z1() < wall.z1()) return 0; // below the wall; unlikely
		return 1;
	}
	point swap_not_dim_z_to_xy(point const &p) const {
		return point(p[!dim], p.z, p[dim]);
	}
	static cube_t cylin_bcube_conservative(cylinder_3dw const &c) {
		cube_t bc(c.p1, c.p2);
		bc.expand_by(max(c.r1, c.r2));
		return bc;
	}
	bool cylins_intersect(cylinder_3dw const &c1, cylinder_3dw const &c2, cube_t const &c1_bc) const { // conservative
		if (!c1_bc.intersects(cylin_bcube_conservative(c2))) return 0;
		float const r_sum(max(c1.r1, c1.r2) + max(c2.r1, c2.r2)); // works best for fixed radius cylinders
		// line_seg_line_seg_dist_2d() expects lines to be in the same XY plane, so swap dims to make that happen; this won't change the distance
		return (line_seg_line_seg_dist_2d(swap_not_dim_z_to_xy(c1.p1), swap_not_dim_z_to_xy(c1.p2), swap_not_dim_z_to_xy(c2.p1), swap_not_dim_z_to_xy(c2.p2)) < r_sum);
	}
public:
	ivy_builder_t(float leaf_sz_, cube_t const &wall_, bool dim_, bool dir_, bool is_house_, rand_gen_t &rgen_) :
		dim(dim_), dir(dir_), is_house(is_house_), leaf_sz(leaf_sz_), wall(wall_), rgen(rgen_) {}

	void next_plant() { // aka clear()
		cylins.clear();
		leaves.clear();
		first_side = rgen.rand_bool();
		main_cylins_end = 0;
	}
	bool add_leaf(point const &pos, vector3d const &branch_dir, vector3d const &side_dir, bool side, float lsz, unsigned cur_branch_leaves_start, bool is_horizontal) {
		//bool const is_hanging(!is_horizontal && fabs(pos[dim] - wall.d[dim][dir]) > leaf_sz);
		float const radius(0.5*lsz), r_test(3.0*radius); // conservative for r_test
		tquad_t leaf(4);

		for (unsigned i = 0; i < 4; ++i) {
			leaf.pts[i] = pos + (radius*((bool(i&1) ^ bool(i&2)) ? 1.0 : -1.0))*branch_dir; // set width
			if (!(i>>1)) {leaf.pts[i] += lsz*side_dir;} // extend away from branch
		}
		rotate_verts(leaf.pts, 2, branch_dir, 0.75*PI_TWO*rgen.rand_float()*(side ? -1.0 : 1.0), pos); // rotate tip away from wall
		cube_t const leaf_bc(leaf.get_bcube());
		if (!is_horizontal && !check_contained_on_wall(leaf_bc)) return 0; // off the wall; allow if horizontal
		point const center(leaf_bc.get_cube_center());

		// check for overlaps with leaves previously added for this plant; allow some amount of overlap; shouldn't Z-fight because dist from wall is random
		for (auto i = leaves.begin(); i < leaves.begin()+cur_branch_leaves_start; ++i) {
			if (dist_less_than(center, i->pts[0], r_test) && dist_less_than(center, i->get_bcube().get_cube_center(), radius)) return 0; // too much overlap
		}
		leaves.push_back(leaf);
		return 1;
	}
	void add_leaves_to_branch(cylinder_3dw const &c) {
		bool const at_branch_end(c.r2 == 0.0), is_horizontal(c.p1.z == c.p2.z);
		float const blen(c.get_length());
		unsigned nleaves(unsigned(1.5 * blen / leaf_sz) + at_branch_end + rgen.rand_bool());
		if (leaves.empty()) {max_eq(nleaves, 1U);} // seed branch must have at least one leaf so that leaves is never empty
		if (nleaves == 0) return; // short branch, no leaves
		unsigned const cur_branch_leaves_start(leaves.size());
		vector3d const branch_delta(c.p2 - c.p1), branch_dir(branch_delta.get_norm());
		vector3d const wall_normal(is_horizontal ? plus_z : vector_from_dim_dir(dim, dir));
		vector3d const side_dir(cross_product(wall_normal, branch_dir).get_norm()); // must normalize for hanging vines
		vector3d const pos_step(branch_delta/(nleaves + !at_branch_end)); // skip last slot if not at branch end because it may be the tip or leaf loc of next segment
		float const r_step((c.r2 - c.r1)/nleaves), rmax(max(c.r1, c.r2));
		point cur_pt(c.p1); // first leaf is at starting point
		float radius(c.r1);

		for (unsigned n = 0; n < nleaves; ++n) {
			float const lsz_scale(max(0.1f, radius/rmax)); // smaller leaves near ends of branches
			float const lsz(rgen.rand_uniform(0.8, 1.4)*lsz_scale*leaf_sz); // leaf size -20% to +40%
			bool const side((leaves.size() & 1) ^ first_side);
			vector3d const side_dir_leaf((side ? 1.0 : -1.0)*side_dir); // in correct direction
			point leaf_pt(cur_pt);
			leaf_pt += 0.1*rgen.signed_rand_float()*pos_step; // add some random position jitter
			leaf_pt += radius*side_dir_leaf; // offset to the side of the branch, alternating sides
			add_leaf(leaf_pt, branch_dir, side_dir_leaf, side, lsz, cur_branch_leaves_start, is_horizontal);
			cur_pt  += pos_step;
			radius  += r_step;
		} // for n
	}
	bool add_branch_seg(point const &p1, point &p2, float r1, float r2, bool is_new_branch, bool &is_horizontal, bool not_on_wall=0) {
		if (is_horizontal) {
			if (not_on_wall) { // check wall length only
				if (p2[!dim] < wall.d[!dim][0] || p2[!dim] > wall.d[!dim][1]) return 0;
			}
			else if (!wall.contains_pt_xy(p2)) return 0; // p2 off the top or end of the wall; allow extending over but check the centerline; p1 is assumed to be valid
		}
		else { // vertical
			cube_t pt_bc;
			pt_bc.set_from_sphere(p2, r2);
			if (!check_contained_on_wall_xy(pt_bc)) return 0; // p2 off the wall; p1 is assumed to be valid
		}
		assert(p1 != p2);
		assert(r1 > 0.0 && r2 > 0.0);
		if (is_horizontal && not_on_wall) {max_eq(p2.z, wall.z2());} // don't fall below the top of the wall
		cylinder_3dw cand(p1, p2, r1, r2);
		cube_t const cand_bc(cylin_bcube_conservative(cand));

		for (cylinder_3dw const &c : cylins) { // does this check is_new_branch?
			if (p1 == c.p1 || p1 == c.p2)           continue; // skip the cylinder we're connected to
			if (cylins_intersect(cand, c, cand_bc)) return 0; // intersects an existing cylinder
		}
		if (!not_on_wall && !is_horizontal && p2.z > wall.z2()) { // off the top of the wall; make horizontal
			if (wall.get_sz_dim(dim) == 0.0) return 0; // zero width wall, can't place ivy on top of it
			// add right angle bend and continue along the top of the wall
			p2.z = cand.p2.z = wall.z2() + r2;
			is_horizontal = has_horizontal = 1;
			branch_bends.emplace_back(p2, r2);
		}
		cylins.push_back(cand);
		return 1;
	}
	void end_branch() {
		assert(!cylins.empty());
		cylinder_3dw &tip(cylins.back());
		tip.r2 = 0.0; // taper the end of the branch
		
		if (!branch_bends.empty() &&tip.p2 == branch_bends.back().pos) { // branch ends at the bend
			branch_bends.pop_back(); // remove bend
			if (tip.get_length() < 4.0*tip.r1) {tip.p2.z += 2.0*tip.r1;} // extend up slightly so that tip isn't too blunt
		}
	}
	bool end_main_branches() {
		main_cylins_end = cylins.size();
		return !cylins.empty();
	}
	cylinder_3dw const &select_random_cylin() const {
		unsigned const cylins_end(main_cylins_end ? main_cylins_end : cylins.size());
		assert(cylins_end > 0 && cylins_end <= cylins.size());
		return cylins[rgen.rand() % cylins_end];
	}
	static void select_split_pos_radius(cylinder_3dw const &c, point &pos, float &radius) {
		if (c.r2 == 0.0) { // if top is a point/branch end, split at the bottom of the cylinder
			pos    = c.p1; // splits at bottom of cylinder
			radius = c.r1;
		}
		else {
			pos    = c.p2; // splits at top of cylinder
			radius = c.r2; // same radius as cylinder
		}
	}
	bool is_at_bend_pt(point const &p) const {
		for (sphere_t const &s : branch_bends) {
			if (p == s.pos) return 1;
		}
		return 0;
	}
	void add_leaves() {
		assert(!cylins.empty());
		for (cylinder_3dw const &c : cylins) {add_leaves_to_branch(c);}
		assert(!leaves.empty());
	}
	void create_branch_verts(vector<vertex_t> &verts, vector<unsigned> &ixs) {
		unsigned const ndiv = 8; // or could go up to 12, but 8 seems like enough
		float const ndiv_inv(1.0/ndiv);
		unsigned data_pos(verts.size()), cylin_ix(0);
		bool prev_horizontal(0);
		point prev_p2;
		assert(!cylins.empty());
		ixs.reserve(ixs.size() + 6*ndiv*(cylins.size() + ndiv*branch_bends.size())); // exact; verts size depends on joins

		for (cylinder_3dw const &c : cylins) {
			bool const is_horizontal(c.p1.z == c.p2.z);
			bool is_join(c.p1 == prev_p2); // this cylinder is joined to the previous cylinder (part of the same branch)
			if (is_horizontal && !prev_horizontal) {is_join = 0;} // vertical => horizontal transition is drawn as a sphere and doesn't join
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
			prev_horizontal = is_horizontal;
			++cylin_ix;
		} // for c
		for (sphere_t const &s : branch_bends) {
			if (sphere_verts.empty()) { // generate unit sphere data for branch bends if we haven't already
				sd_sphere_d sd(all_zeros, 1.0, ndiv);
				sd.gen_points_norms(spn);	
				sd.get_itri_points(sphere_verts, sphere_ixs);
				assert(!(sphere_ixs.size() % 6)); // must be a multiple of 6 (triangle pairs)
			}
			unsigned const verts_start(verts.size());
			for (unsigned   ix : sphere_ixs  ) {ixs.push_back(verts_start + ix);}
			for (auto const &v : sphere_verts) {verts.emplace_back((s.radius*v.v + s.pos), v.n, v.t[0], v.t[1]);}
		} // for s
	}
	void create_leaf_verts(vector<vertex_t> &verts) const {
		// Note: no reserve(), since we're adding to existing ivy plant leaves
		float const tcs[4] = {0.0, 1.0, 1.0, 0.0}, tct[4] = {0.0, 0.0, 1.0, 1.0};
		assert(!leaves.empty());
		verts.reserve(verts.size() + 4*leaves.size()); // exact

		for (tquad_t const &l : leaves) {
			vector3d const leaf_normal(l.get_norm());
			for (unsigned i = 0; i < 4; ++i) {verts.emplace_back(l.pts[i], leaf_normal, tcs[i], tct[i]);}
		}
	}
}; // end ivy_builder_t

void ivy_wall_t::place_on_wall_face(cube_t const &wall, vect_cube_t const &avoid, bool dim, bool dir, bool is_house, float leaf_sz,
	vector<vertex_t> &lverts, vector<vertex_t> &bverts, vector<unsigned> &bixs, rand_gen_t &rgen)
{
	// generation steps:
	// * select random start points at the base of the wall
	// * create branch with upward direction and random curve
	// * starting at random points along a branch, create a split point that becomes a new branch and widens the branch segments below
	// * place leaves along branches with stem touching the branch
	// * give leaves random orients, but not enough that they clip through walls
	// * check that leaves don't overlap too much with other leaves
	float const wall_len(wall.get_sz_dim(!dim)), wall_edge_space(4.0*leaf_sz), root_spacing(8.0*leaf_sz);
	if (wall_len <= 2.0*wall_edge_space) return; // wall too short; shouldn't happen
	float const pos_lo(wall.d[!dim][0] + wall_edge_space), pos_hi(wall.d[!dim][1] - wall_edge_space);
	float const wall_thick(wall.get_sz_dim(dim)), wall_face(wall.d[dim][dir]), dsign(dir ? 1.0 : -1.0);
	vector3d const wall_normal(vector_from_dim_dir(dim, dir));
	unsigned const num_roots(4 + (rgen.rand() % 9)); // 4-12
	ivy_builder_t builder(leaf_sz, wall, dim, dir, is_house, rgen);
	vector<float> root_vals;

	for (unsigned n = 0; n < num_roots; ++n) {
		float const branch_radius(rgen.rand_uniform(0.08, 0.12)*leaf_sz), seg_len(12.0*branch_radius);
		// find root pos at base of wall that's not too close to a previous root
		bool root_valid(0);
		point root;
		root.z    = wall.z1(); // assume this is the ground
		root[dim] = wall_face + branch_radius*dsign; // move slightly away from the surface

		for (unsigned N = 0; N < 100 && !root_valid; ++N) {
			root[!dim] = rgen.rand_uniform(pos_lo, pos_hi);
			if (check_vect_cube_contains_pt(avoid, root)) continue; // blocked
			root_valid = 1;

			for (float const &v : root_vals) {
				if (fabs(root[!dim] - v) < root_spacing) {root_valid = 0; break;}
			}
			root_vals.push_back(root[!dim]);
		} // for N
		if (!root_valid) break; // no more roots can be placed
		// add branches; main branch in cylinder segments, then connect secondary branches
		unsigned const num_branches(6 + (rgen.rand() % 3)); // 6-8
		point pos(root);
		float radius(branch_radius);
		builder.next_plant();

		for (unsigned B = 0; B < num_branches; ++B) {
			unsigned num_segs(10 + (rgen.rand() % 6)); // 10-15
			vector3d prev_dir;
			bool is_horizontal(0);
			
			if (B > 0) { // add secondary branch
				cylinder_3dw const &c(builder.select_random_cylin());
				builder.select_split_pos_radius(c, pos, radius);
				is_horizontal = (c.p1.z == c.p2.z || builder.is_at_bend_pt(c.p1)); // starting at bend point will be horizontal since it can't ascend any more
				prev_dir = (c.p2 - c.p1).get_norm();
				radius   = max(0.65f*branch_radius, rgen.rand_uniform(0.7, 0.9)*radius); // smaller radius
				if (is_horizontal) {pos.z = wall.z2() + branch_radius;} // move closer to top of wall
				else {pos[dim] = wall_face + branch_radius*dsign;} // move closer to side of wall
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
							vector3d const vrot(is_horizontal ? plus_z : wall_normal);
							vector3d new_dir(prev_dir);
							rotate_vector3d(vrot, TO_RADIANS*rgen.rand_uniform(30.0, 60.0)*rgen.rand_sign(), new_dir); // 30-60 deg from prev branch
							pos2 += (rgen.rand_uniform(0.8, 1.2)*seg_len)*new_dir; // start in the new direction
						}
					}
					else { // continuation; should this be a curve rather than random?
						vector3d new_dir(prev_dir);
						float const rot_angle(10.0*TO_RADIANS*rgen.signed_rand_float()); // +/- 10 deg from prev branch

						if (is_horizontal) {
							rotate_vector3d(plus_z, rot_angle, new_dir);
						}
						else { // vertical
							rotate_vector3d(wall_normal, rot_angle, new_dir);

							if (new_dir.z < 0.0 && new_dir.z < prev_dir.z) { // should move upward
								// can't clamp because angle may be too sharp from prev_dir, so rotate in the other dir if we pointed more downward
								new_dir = prev_dir;
								rotate_vector3d(wall_normal, -rot_angle, new_dir);
							}
						}
						pos2 += (rgen.rand_uniform(0.8, 1.2)*seg_len)*new_dir;
					}
					if (check_vect_cube_contains_pt(avoid, pos2)) continue; // blocked
					bool const was_horizontal(is_horizontal);
					if (!builder.add_branch_seg(pos, pos2, radius, radius, (S == 0), is_horizontal)) continue; // constant radius; reject if branch can't be placed
					if (was_horizontal) {ivy_faces |= 3;} // mark both faces as having ivy since the top of the wall will be visible from both sides
					
					if (is_horizontal && !was_horizontal) { // crossing the top of the wall; choose a new random XY dir that won't cross the opposite side of the wall
						prev_dir[dim] = rgen.rand_uniform(0.1, 1.0)*min(wall_thick, seg_len); // move toward the opposite side of the wall
						prev_dir.z    = 0.0; // prev_dir[!dim] remains unchanged
						prev_dir.normalize();
						if (S+1 == num_segs) {++num_segs;} // don't end on the bend point
					}
					else {prev_dir = (pos2 - pos).get_norm();} // vertical
					pos    = pos2;
					placed = 1;
					break; // done
				} // for N
				if (!placed) break; // can't place any more segments on this branch
			} // for S
			builder.end_branch();
		} // for B
		// add branches sticking out into the air on the sides and top if not well maintained
		if (maintained_amt < 1.0) { // maybe add upward curve
			bool const has_branches(builder.end_main_branches());
			unsigned const num_extra_branches(has_branches ? round_fp(8.0*(1.0 - maintained_amt)*rgen.rand_float()) : 0); // 0-8

			for (unsigned B = 0; B < num_extra_branches; ++B) {
				unsigned num_segs(6 + (rgen.rand() % 5)); // 6-10
				cylinder_3dw const &c(builder.select_random_cylin());
				builder.select_split_pos_radius(c, pos, radius);
				bool is_horizontal(c.p1.z == c.p2.z || builder.is_at_bend_pt(c.p1));
				radius = max(0.65f*branch_radius, rgen.rand_uniform(0.7, 0.9)*radius); // smaller radius
				vector3d cur_dir;
				// Note: no need to move closer to wall if radius is reduced

				for (unsigned S = 0; S < num_segs; ++S) {
					bool placed(0);

					for (unsigned N = 0; N < 10; ++N) { // 10 tries to place a branch segment
						vector3d new_dir(cur_dir), adj_dir(rgen.signed_rand_vector_spherical_norm());
						if (!is_horizontal && dot_product(adj_dir, wall_normal) < 0.0) {adj_dir.negate();} // face away from the wall
						adj_dir.z = fabs(adj_dir.z); // must point upward
						if (S == 0) {new_dir = adj_dir;} // new branch
						else {new_dir = (0.9*new_dir + 0.1*adj_dir).get_norm();} // continuation
						point pos2(pos + (rgen.rand_uniform(0.7, 1.0)*seg_len)*new_dir); // start in the new direction
						if (check_vect_cube_contains_pt(avoid, pos2)) continue; // blocked
						if (!builder.add_branch_seg(pos, pos2, radius, radius, (S == 0), is_horizontal, 1)) continue; // not_on_wall=1
						new_dir.z -= 0.2; // sags with the weight of each new segment
						cur_dir = new_dir;
						pos     = pos2;
						placed  = 1;
						break; // done
					} // for N
					if (!placed) break; // can't place any more segments on this branch
				} // for S
				builder.end_branch();
			} // for B
		}
		builder.add_leaves();
		builder.create_branch_verts(bverts, bixs);
		builder.create_leaf_verts  (lverts);
	} // for n
	ivy_faces |= (1 << (unsigned)dir); // mark face as having ivy
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
void ivy_manager_t::add_wall(cube_t const &wall, vect_cube_t const &avoid, bool dim, bool is_house, unsigned skip_dirs,
	unsigned wall_ix, unsigned plot_ix, unsigned city_ix, float leaf_sz, point const &camera_bs)
{
	if (city_ix != cur_city_ix) { // city change
		cur_city_ix = city_ix;
		clear();
	}
	bool const rand_select(skip_dirs == 0); // house walls have skip_dirs set and all enabled have ivy

	if (rand_select) { // apply filtering for plot divider wall
		if (((13*plot_ix) % 5) == 0) return; // some plots have no ivy
		if (((17*wall_ix) % 5) == 0) return; // some walls have no ivy
	}
	ivy_wall_t &w(ivy_walls[wall_ix]);

	if (w.leaves.bcube.is_all_zeros()) { // new wall
		unsigned face_mask(dim ? 12 : 3); // either both X sides or both Y sides
		if (skip_dirs & 1) {face_mask &=  5;} // skip lo dirs
		if (skip_dirs & 2) {face_mask &= 10;} // skip hi dirs
		rand_gen_t rgen;
		rgen.set_state(wall_ix+1, plot_ix+1);
		rgen.rand_mix();
		w.gen(wall, avoid, leaf_sz, face_mask, rand_select, is_house, rgen);
	}
	else { // existing wall
		assert(w.leaves.bcube == wall && w.branches.bcube == wall);
	}
	if (w.empty()) return; // no ivy
	assert(w.ivy_faces);

	if (skip_dirs == 0) { // back face culling is only for walls, since ivy may be visible through house windows
		if (w.ivy_faces == 1 && camera_bs[dim] > wall.d[dim][1]) return; // facing opposite side
		if (w.ivy_faces == 2 && camera_bs[dim] < wall.d[dim][0]) return; // facing opposite side
	}
	to_draw.push_back(wall_ix); // not checked for duplicates, but there shouldn't be any
}

void drawable_t::draw() const {
	if (num_verts == 0) return;
	bool const indexed(num_ixs > 0);
	pre_render(indexed, indexed); // do_bind_vbo=indexed
	if (indexed) {draw_indexed_tri_verts(num_verts, num_ixs, GL_TRIANGLES);} // indexed triangles
	else {draw_quads_as_tris(num_verts);} // quads
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
	bind_vbo(0, 1); // unbind index VBO
	// draw leaves second; ideally we want these to have two sided leaf lighting, but that's not available in this shader
	begin_leaf_draw(s, PLANT1_TEX);
	for (uint32_t wix : to_draw) {ivy_walls.find(wix)->second.leaves.draw();}
	end_leaf_draw(s);
	to_draw.clear();
}

// ponds with lily pads and cat tails

bool has_circle_overlap(sphere_t const &circle, vector<sphere_t> const &circles) {
	for (sphere_t const &c : circles) {
		if (dist_xy_less_than(circle.pos, c.pos, (circle.radius + c.radius))) return 1;
	}
	return 0;
}
bool point_in_ellipse(point const &p, cube_t const &c);

// untextured; can be used with quad_batch_draw
void add_vert_circle_verts(point const &pos, float radius, unsigned ndiv, color_wrapper const &cw, vector<vert_norm_tc_color> &verts) {
	vector_point_norm const &vpn(gen_cylinder_data(pos, pos, radius, radius, ndiv));
	unsigned const center_ix(verts.size());
	verts.emplace_back(pos, plus_z, 0.0, 0.0, cw); // center

	for (unsigned S = 0; S <= ndiv; ++S) {
		if (S > 1) { // tri_fan_push
			verts.push_back(verts[center_ix]);
			verts.push_back(verts[verts.size()-2]);
		}
		verts.emplace_back(vpn.p[(S%ndiv)<<1], plus_z, 0.0, 0.0, cw);
	}
}

void pond_t::gen_vegetation(park_heightmap_t const &hmap) {
	// add lily pads
	rand_gen_t rgen;
	rgen.set_state(rseed, 3*rseed+1);
	unsigned const num_lp(10 + (rgen.rand() % 31)); // 10-40
	lily_pads.reserve(num_lp);

	for (unsigned n = 0; n < num_lp; ++n) {
		sphere_t lpad;
		lpad.pos.z  = (rgen.rand() % 8); // stores orient
		lpad.radius = rgen.rand_uniform(0.5, 1.0)*0.03*radius;
		cube_t lp_area(bcube), lp_avoid(bcube);
		lp_area.expand_by_xy(-2.0*lpad.radius); // don't place too close to the pond edge
		lp_avoid.expand_by(-vector3d(0.4*bcube.dx(), 0.4*bcube.dy(), 0.0)); // avoid the center
		bool success(0);

		for (unsigned N = 0; N < 100; ++N) {
			for (unsigned d = 0; d < 2; ++d) {lpad.pos[d] = rgen.rand_uniform(bcube.d[d][0], bcube.d[d][1]);}
			if (point_in_ellipse(lpad.pos, lp_area) && !point_in_ellipse(lpad.pos, lp_avoid) && !has_circle_overlap(lpad, lily_pads)) {success = 1; break;}
		}
		if (!success) break; // shouldn't happen
		lily_pads.push_back(lpad);
	} // for n
	// add cattails
	float const depth(bcube.dz()), ct_max_height(1.0*depth), ct_max_radius(0.1*ct_max_height); // depth is slightly larger than actual pond depth
	float const ct_min_z1(get_water_zval() - 0.25*depth), ct_max_z1(get_water_zval() - 0.02*depth), clump_radius(0.02*(bcube.dx() + bcube.dy()));
	unsigned const num_ct(30 + (rgen.rand() % 21)); // 30-50
	cat_tails.reserve(num_ct);
	point ct_pos;

	for (unsigned n = 0; n < num_ct; ++n) {
		// place around perimeter in shallow water, clustered into groups
		for (unsigned N = 0; N < 100; ++N) { // 100 random tries
			if (N < 2 && ct_pos != all_zeros) { // first 2 iterations: place close to prev cat tail to produce clumps
				ct_pos += clump_radius*rgen.signed_rand_vector_spherical_xy_norm();
				if (!bcube.contains_pt_xy(ct_pos)) continue;
			}
			else {
				gen_xy_pos_in_cube(ct_pos, bcube, rgen);
			}
			ct_pos.z = hmap.get_zval_at_pos(cube_bot_center(ct_pos)); // use zval from the bottom of the pond
			if (ct_pos.z < ct_min_z1 || ct_pos.z > ct_max_z1) continue; // too deep or too shallow
			cube_t ctail(ct_pos);
			ctail.z2() = ct_pos.z + rgen.rand_uniform(0.75, 1.0)*ct_max_height;
			float const ct_radius(rgen.rand_uniform(0.75, 1.0)*ct_max_radius);
			ctail.expand_by_xy(ct_radius);
			if (has_bcube_int_xy(ctail, cat_tails)) continue; // too close to another cat tail
			if (has_circle_overlap(sphere_t(ct_pos, ct_radius), lily_pads)) continue; // too close to a lily pad
			cat_tails.push_back(ctail);
			break;
		} // for N
	} // for n
}

void pond_t::draw_lily_pads(draw_state_t &dstate, city_draw_qbds_t &qbds, bool shadow_only, float dist) const {
	bool const draw_bot(dstate.camera_bs.z < get_water_zval());
	float const dz_off((draw_bot ? -1.0 : 1.0)*max(0.0001f*bcube.dz(), 0.00025f*dist)), z1(get_water_zval() + 2.0*dz_off), z2(z1 + dz_off);
	color_wrapper const cw(WHITE);
	if (!shadow_only) {select_texture(get_texture_by_name("lilypad.png"));} // set in case we drew cat tails previously

	for (sphere_t const &lp : lily_pads) { // draw lily pads
		if (shadow_only) { // circular
			add_vert_circle_verts(point(lp.pos.x, lp.pos.y, z1), lp.radius, 16, cw, qbds.qbd.verts); // ndiv=16
		}
		else { // textured quad
			unsigned const orient(round_fp(lp.pos.z));
			bool const mx(orient & 1), my(orient & 2), swap_xy(orient & 4);
			cube_t lpad;
			lpad.set_from_sphere(lp);
			set_cube_zvals(lpad, z1, z2);
			dstate.draw_cube(qbds.qbd, lpad, cw, !draw_bot, 0.0, 3, mx, my, swap_xy, 1.0, 1.0, 1.0, draw_bot);
		}
	} // for lp
}

struct plant_bender_t {
	float bend_amt;
	point bot;

	plant_bender_t(point const &bot_, point const &top, float flex_val) : bot(bot_) {
		assert(bot.z < top.z);
		vector3d const delta(top - bot);
		bend_amt = flex_val*min(2.0f, delta.xy_mag()/delta.z)/delta.mag();
	}
	void bend(point &pos) {
		float const length(p2p_dist(pos, bot));
		pos.z = max(bot.z, (pos.z - bend_amt*length*length)); // quadratic bend function; can't go below the base
		vector3d const new_delta(pos - bot);
		pos = bot + new_delta*(length/new_delta.mag()); // adjust to preserve original segment length but keep the new direction
	}
};

void add_back_side_verts(vector<vert_norm_tc> &verts, unsigned start_ix) {
	unsigned const end_ix(verts.size());
	assert(start_ix <= end_ix);

	for (unsigned i = start_ix; i < end_ix; ++i) {
		verts.push_back(verts[i]);
		verts.back().invert_normal();
	}
	reverse(verts.begin()+end_ix, verts.end()); // reverse winding order
}

void pond_t::draw_cat_tails(draw_state_t &dstate, bool shadow_only) const {
	//if (shadow_only) return; // stem and leaves are too thin and shadows cause artifacts
	//highres_timer_t timer("draw_cat_tails"); // 6.2ms => 3.4ms => 2.4ms => 1.9ms with 1000 cat tails
	ctdd.create_verts(cat_tails, rseed);
	ctdd.draw(dstate, shadow_only);
}

void pond_t::cat_tail_draw_data_t::create_verts(vect_cube_t const &cat_tails, unsigned rseed) {
	if (!leaf_qverts.empty()) return; // already setup
	unsigned const ndiv(16), nstacks(8);
	float const nstacks_inv(1.0/nstacks);
	rand_gen_t rgen;
	rgen.set_state(rseed, 3*rseed+1);

	// generate sphere verts once and reuse with correct pos and radius
	vector<vert_norm_tc> sphere_verts;
	vector<unsigned> sphere_ixs;
	{ // open a scope for sd/spn
		sd_sphere_d sd(all_zeros, 1.0, ndiv);
		sphere_point_norm spn;
		sd.gen_points_norms(spn);
		sd.get_itri_points(sphere_verts, sphere_ixs);
	}
	for (cube_t const &ct : cat_tails) {
		// stem
		float const ct_height(ct.dz()), ct_radius(0.5*ct.dx()), lean_amt(rgen.rand_uniform(0.0, 2.0));
		float const stem_radius(0.07*ct_radius), stem_radius_step(nstacks_inv*stem_radius);
		point const bot(cube_bot_center(ct));
		point const top(cube_top_center(ct) + lean_amt*ct_radius*rgen.signed_rand_vector_spherical_xy_norm()); // random lean
		vector3d const stem_delta(top - bot), stack_step(nstacks_inv*stem_delta);
		float cur_radius(stem_radius), cur_ts(0.0);
		point cur_stem(bot), prev_stem_draw(cur_stem);
		plant_bender_t stem_bender(bot, top, 2.8);

		for (unsigned s = 0; s < nstacks; ++s) { // split into nstacks cylinder + cone segments and bend them
			point const next_stem(cur_stem + stack_step);
			point next_stem_draw(next_stem);
			stem_bender.bend(next_stem_draw); // apply bend; cylinder ends won't line up exactly right, but it's not too noticeable

			if (s+1 < nstacks) { // truncated cone (quads)
				float const t((s+1)*nstacks_inv), next_radius((1.0 - t*t)*stem_radius); // slow sqrt taper
				assert(next_radius > 0.0);
				gen_cylinder_quads(leaf_qverts, gen_cylinder_data(prev_stem_draw, next_stem_draw, cur_radius, next_radius, ndiv), 0, nstacks_inv, cur_ts);
				// can we join the verts and normals between this cylinder and the previous one?
				cur_radius = next_radius;
				cur_ts    += nstacks_inv;
				cur_stem   = next_stem;
				prev_stem_draw = next_stem_draw;
			}
			else { // cone for last stack (triangles)
				gen_cone_triangles(leaf_tverts, gen_cylinder_data(prev_stem_draw, next_stem_draw, cur_radius, 0.0, ndiv), 0, cur_ts, 1.0);
			}
		} // for s
		// leaves
		unsigned const nleaves(5 + (rgen.rand() % 4)); // 5-8

		for (unsigned n = 0; n < nleaves; ++n) {
			unsigned const qv_start(leaf_qverts.size()), tv_start(leaf_tverts.size());
			float const leaf_len(rgen.rand_uniform(0.42, 0.75)*ct_height), leaf_base_hwidth(rgen.rand_uniform(0.75, 1.0)*0.02*leaf_len), bend_amt(rgen.rand_uniform(0.6, 0.9));
			vector3d const dir(rgen.signed_rand_vector_spherical_xy_norm()), side_dir(cross_product(dir, plus_z)), side_delta(leaf_base_hwidth*side_dir);
			point const base(bot + dir*ct_radius*rgen.rand_uniform(0.1, 0.2)), tip(base + dir*ct_radius*bend_amt + leaf_len*plus_z);
			vector3d const delta(tip - base);
			plant_bender_t leaf_bender(base, tip, 4.0);
			point pa(bot), p1a(pa - side_delta), p2a(pa + side_delta);
			vector3d normal1(-dir);

			for (unsigned s = 0; s < nstacks; ++s) {
				float const t1(s*nstacks_inv), t2(t1 + nstacks_inv);
				point pb(base + t2*delta);
				leaf_bender.bend(pb);
				vector3d const normal2(-cross_product((pb - pa), side_dir).get_norm()); // normal of segment

				if (s+1 < nstacks) { // quad
					vector3d const side_vb((1.0 - t2*t2*t2*t2)*side_delta); // slow taper
					point const p1b(pb - side_vb), p2b(pb + side_vb);
					leaf_qverts.emplace_back(p1a, normal1, 0.0, t1);
					leaf_qverts.emplace_back(p2a, normal1, 1.0, t1);
					leaf_qverts.emplace_back(p2b, normal2, 1.0, t2);
					leaf_qverts.emplace_back(p1b, normal2, 0.0, t2);
					pa = pb; p1a = p1b; p2a = p2b; normal1 = normal2;
				}
				else { // triangle tip
					leaf_tverts.emplace_back(p1a, normal1, 0.0, t1 );
					leaf_tverts.emplace_back(p2a, normal1, 1.0, t1 );
					leaf_tverts.emplace_back(pb,  normal2, 0.5, 1.0);
				}
			} // for s
			add_back_side_verts(leaf_qverts, qv_start);
			add_back_side_verts(leaf_tverts, tv_start);
		} // for n
		// capsule (flower head)
		float const cap_t1(rgen.rand_uniform(0.75, 0.8)), cap_t2(cap_t1 + rgen.rand_uniform(0.08, 0.12)), cap_radius(0.2*ct_radius); // parametric position along bot=>top
		point cap_pos[2];

		for (unsigned d = 0; d < 2; ++d) {
			cap_pos[d] = bot + (d ? cap_t2 : cap_t1)*stem_delta;
			stem_bender.bend(cap_pos[d]);
		}
		// stretch a sphere into a capsule shape
		unsigned const num_sv(sphere_verts.size()), ix_off(cap_verts.size());
		vector3d const cap_delta(cap_pos[1] - cap_pos[0]);
		point const cap_center(0.5*(cap_pos[1] + cap_pos[0]));
		float const cap_dz(0.5*cap_delta.mag()/cap_radius);
		vector<point> verts(num_sv), norms(num_sv);

		for (unsigned i = 0; i < num_sv; ++i) {
			point &v(verts[i]);
			v = sphere_verts[i].v;
			float const t(sqrt(fabs(v.z))); // 0 at center, 1 at ends
			v.z += SIGN(v.z)*t*cap_dz; // stretch out ends
			norms[i] = sphere_verts[i].n; // not quite correct for points along the sphere center/capsule sides
			norms[i].z *= t; // FIXME: not quite correct
			norms[i].normalize();
		} // for i
		rotate_vector3d_by_vr_multi(plus_z, cap_delta, verts.data(), num_sv);
		rotate_vector3d_by_vr_multi(plus_z, cap_delta, norms.data(), num_sv);
		for (unsigned ix : sphere_ixs) {cap_ixs.push_back(ix + ix_off);}
		for (unsigned i = 0; i < num_sv; ++i) {cap_verts.emplace_back((cap_radius*verts[i] + cap_center), norms[i], sphere_verts[i].t);}
	} // for ct
}

void pond_t::cat_tail_draw_data_t::draw(draw_state_t &dstate, bool shadow_only) const {
	glEnable(GL_CULL_FACE); // so that correct face of leaves is drawn
	select_no_texture();
	dstate.s.set_cur_color(BROWN);
	set_ptr_state(cap_verts.data(), cap_verts.size(), 0, 1);
	draw_indexed_tri_verts(cap_verts.size(), cap_ixs.size(), GL_TRIANGLES, (void *)cap_ixs.data());
	++num_frame_draw_calls;
	unset_ptr_state(cap_verts.data());
	dstate.s.set_cur_color(colorRGBA(0.4, 0.6, 0.2)); // make the green even darker
	if (!shadow_only) {select_texture(GRASS_BLADE_TEX);}
	draw_verts(leaf_tverts, GL_TRIANGLES);
	draw_quad_verts_as_tris(leaf_qverts);
	glDisable(GL_CULL_FACE);
}


