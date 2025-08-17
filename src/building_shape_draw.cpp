// 3D World - Building Room Geometry Primitive Shape Drawing
// by Frank Gennari 8/10/25

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "subdiv.h" // for sd_sphere_d


// skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2 to match CSG cube flags
void rgeom_mat_t::add_cube_to_verts(cube_t const &c, colorRGBA const &color, point const &tex_origin, unsigned skip_faces,
	bool swap_tex_st, bool mirror_x, bool mirror_y, bool inverted, bool z_dim_uses_ty, float tx_add, float ty_add)
{
	//assert(c.is_normalized()); // no, bathroom window is denormalized
	vertex_t v;
	v.set_c4(color);
	tx_add += tex.txoff;
	ty_add += tex.tyoff;

	// Note: stolen from draw_cube() with tex coord logic, back face culling, etc. removed
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions, drawn as {Z, X, Y}
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);
		bool const tex_st(swap_tex_st ^ (z_dim_uses_ty && (d[1] == 2)));

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (skip_faces & (1 << (2*(2-n) + j))) continue; // skip this face
			v.set_ortho_norm(n, (bool(j) ^ inverted));
			v.v[n] = c.d[n][j];

			for (unsigned s1 = 0; s1 < 2; ++s1) {
				v.v[d[1]] = c.d[d[1]][s1];
				v.t[tex_st] = ((tex.tscale_x == 0.0) ? float(s1) : (tex.tscale_x*(v.v[d[1]] - tex_origin[d[1]]) + tx_add)); // tscale==0.0 => fit texture to cube

				for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
					bool const s2(bool(k^j^s1)^inverted^1); // need to orient the vertices differently for each side
					v.v[d[0]] = c.d[d[0]][s2];
					v.t[!tex_st] = ((tex.tscale_y == 0.0) ? float(s2) : (tex.tscale_y*(v.v[d[0]] - tex_origin[d[0]]) + ty_add));
					quad_verts.push_back(v);
					if (mirror_x) {quad_verts.back().t[0] = 1.0 - v.t[0];} // use for pictures and books
					if (mirror_y) {quad_verts.back().t[1] = 1.0 - v.t[1];} // used for books
				} // for k
			} // for s1
		} // for j
	} // for i
}
void rgeom_mat_t::add_cube_to_verts_two_sided(cube_t const &c, colorRGBA const &color, point const &tex_origin,
unsigned skip_faces, bool swap_tex_st, bool mirror_x, bool mirror_y)
{
	for (unsigned s = 0; s < 2; ++s) {add_cube_to_verts(c, color, tex_origin, skip_faces, swap_tex_st, mirror_x, mirror_y, bool(s));}
}
// untextured version of the above function
void rgeom_mat_t::add_cube_to_verts_untextured(cube_t const &c, colorRGBA const &color, unsigned skip_faces) { // add an inverted flag?
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

void apply_half_or_quarter(int half_or_quarter, unsigned &s_end) { // 0=full circle, 1=half circle, 2=quarter circle, 3=half a full circle in the other dim, 4=eighth circle
	if      (half_or_quarter == 0) {} // full
	else if (half_or_quarter == 1) {s_end /= 2;} // half
	else if (half_or_quarter == 2) {s_end /= 4;} // quarter
	else if (half_or_quarter == 3) {} // half in other dim - not handled here
	else if (half_or_quarter == 4) {s_end /= 8;} // eighth
	else {assert(0);}
}

void rgeom_mat_t::add_ortho_cylin_to_verts(cube_t const &c, colorRGBA const &color, int dim, bool draw_bot, bool draw_top, bool two_sided,
	bool inv_tb, float rs_bot, float rs_top, float side_tscale, float end_tscale, bool skip_sides, unsigned ndiv, float side_tscale_add,
	bool swap_txy, float len_tc2, float len_tc1, int half_or_quarter, int swap_end_txy)
{
	if (dim == 2) { // Z: this is our standard v_cylinder
		add_vcylin_to_verts(c, color, draw_bot, draw_top, two_sided, inv_tb, rs_bot, rs_top, side_tscale, end_tscale, skip_sides,
			ndiv, side_tscale_add, swap_txy, len_tc2, len_tc1, half_or_quarter, swap_end_txy);
		return;
	}
	cube_t c_rot(c);
	c_rot.swap_dims(2, dim);
	unsigned const itri_verts_start_ix(itri_verts.size()), ixs_start_ix(indices.size());
	add_vcylin_to_verts(c_rot, color, draw_bot, draw_top, two_sided, inv_tb, rs_bot, rs_top, side_tscale, end_tscale, skip_sides,
		ndiv, side_tscale_add, swap_txy, len_tc2, len_tc1, half_or_quarter, swap_end_txy);
	for (auto v = itri_verts.begin()+itri_verts_start_ix; v != itri_verts.end(); ++v) {v->swap_dims(2, dim);} // swap triangle vertices and normals
	std::reverse(indices.begin()+ixs_start_ix, indices.end()); // fix winding order
}
void rgeom_mat_t::add_vcylin_to_verts(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top, bool two_sided,
	bool inv_tb, float rs_bot, float rs_top, float side_tscale, float end_tscale, bool skip_sides, unsigned ndiv, float side_tscale_add,
	bool swap_txy, float len_tc2, float len_tc1, int half_or_quarter, int swap_end_txy)
{
	point const center(c.get_cube_center());
	float const radius(0.5*min(c.dx(), c.dy())); // cube X/Y size should be equal/square
	add_cylin_to_verts(point(center.x, center.y, c.z1()), point(center.x, center.y, c.z2()), radius*rs_bot, radius*rs_top, color, draw_bot, draw_top,
		two_sided, inv_tb, side_tscale, end_tscale, skip_sides, ndiv, side_tscale_add, swap_txy, len_tc2, len_tc1, half_or_quarter, swap_end_txy);
}
void rgeom_mat_t::add_vcylin_to_verts_tscale(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top) {
	vector3d const sz(c.get_size());
	// make side_tscale an exact multiple of 1.0 so that there are no seams
	float const side_tscale(round_fp(max(1.0, tex.tscale_x*0.5*PI*(sz.x + sz.y)))), len_tscale(tex.tscale_y*sz.z), end_tscale(0.25*(tex.tscale_x*sz.x + tex.tscale_y*sz.y));
	add_vcylin_to_verts(c, color, draw_bot, draw_top, 0, 0, 1.0, 1.0, side_tscale, end_tscale, 0, N_CYL_SIDES, 0.0, 0, len_tscale);
}
void rgeom_mat_t::add_cylin_to_verts(point const &bot, point const &top, float bot_radius, float top_radius, colorRGBA const &color,
	bool draw_bot, bool draw_top, bool two_sided, bool inv_tb, float side_tscale, float end_tscale, bool skip_sides, unsigned ndiv,
	float side_tscale_add, bool swap_txy, float len_tc2, float len_tc1, int half_or_quarter, int swap_end_txy)
{
	assert((!skip_sides) || draw_bot || draw_top); // must draw something
	point const ce[2] = {bot, top};
	float const ndiv_inv(1.0/ndiv), half_end_tscale(0.5*end_tscale);
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, bot_radius, top_radius, ndiv, v12));
	color_wrapper const cw(color);
	unsigned itris_start(itri_verts.size()), ixs_start(indices.size()), itix(itris_start), iix(ixs_start);

	if (!skip_sides) {
		unsigned const ixs_off[6] = {1,2,0, 3,2,1}; // 1 quad = 2 triangles
		bool const flat_sides(ndiv <= 6 && side_tscale == 0.0); // hack to draw bolts untextured with flat sides, since no other cylinders have only 6 sides
		unsigned ndiv_draw(ndiv);
		assert(half_or_quarter != 3); // no half-in-other-dim
		apply_half_or_quarter(half_or_quarter, ndiv_draw);
		unsigned const num_side_verts(flat_sides ? 4*ndiv_draw : 2*(ndiv_draw+1)), unique_verts_per_side(flat_sides ? 4 : 2);
		itri_verts.resize(itris_start + num_side_verts);
		indices.resize(ixs_start + 6*ndiv_draw);

		if (flat_sides) {
			for (unsigned i = 0; i < ndiv_draw; ++i) { // vertex data
				unsigned const in((i+1)%ndiv);
				point const pts[4] = {vpn.p[(i<<1)+0], vpn.p[(i<<1)+1], vpn.p[(in<<1)+0], vpn.p[(in<<1)+1]};
				norm_comp const normal(get_poly_norm(pts));
				for (unsigned n = 0; n < 4; ++n) {itri_verts[itix++].assign(pts[n], normal, 0.0, 0.0, cw);} // all tcs=0
			}
		}
		else {
			for (unsigned i = 0; i <= ndiv_draw; ++i) { // vertex data
				unsigned const s(i%ndiv);
				float const ts(side_tscale*(1.0f - i*ndiv_inv) + side_tscale_add);
				norm_comp const normal(0.5*(vpn.n[s] + vpn.n[(i+ndiv-1)%ndiv])); // normalize?
				itri_verts[itix++].assign(vpn.p[(s<<1)+0], normal, (swap_txy ? len_tc1 : ts), (swap_txy ? ts : len_tc1), cw);
				itri_verts[itix++].assign(vpn.p[(s<<1)+1], normal, (swap_txy ? len_tc2 : ts), (swap_txy ? ts : len_tc2), cw);
			}
		}
		for (unsigned i = 0; i < ndiv_draw; ++i) { // index data
			unsigned const ix0(itris_start + unique_verts_per_side*i);
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
		assert(half_or_quarter == 0); // half and quarter disk are not supported
		norm_comp const normal((bool(bt) ^ inv_tb) ? v12 : -v12);
		unsigned const center_ix(itix);
		itri_verts[itix++].assign(ce[bt], normal, half_end_tscale, half_end_tscale, cw); // center

		for (unsigned I = 0; I < ndiv; ++I) {
			unsigned const i(bt ? ndiv-I-1 : I); // invert winding order for top face
			vector3d const side_normal(0.5*(vpn.n[i] + vpn.n[(i+ndiv-1)%ndiv])); // normalize?
			float const ts(half_end_tscale*(side_normal.x*((swap_end_txy & 1) ? -1.0 : 1.0) + 1.0));
			float const tt(half_end_tscale*(side_normal.y*((swap_end_txy & 2) ? -1.0 : 1.0) + 1.0));
			itri_verts[itix++].assign(vpn.p[(i<<1) + bt], normal, ts, tt, cw); // assign tcs from side normal
			indices[iix++] = center_ix; // center
			indices[iix++] = center_ix + i + 1;
			indices[iix++] = center_ix + ((i+1)%ndiv) + 1;
		}
	} // for bt
	if (inv_tb) {std::reverse(indices.begin()+ixs_start, indices.end());} // reverse the order to swap triangle winding order
	if (two_sided) {add_inverted_triangles(itri_verts, indices, itris_start, ixs_start);}
}

void rgeom_mat_t::add_disk_to_verts(point const &pos, float radius, vector3d const &dir, colorRGBA const &color, bool swap_txy, bool inv_ts, bool inv_tt, unsigned ndiv) {
	assert(radius > 0.0);
	color_wrapper const cw(color);
	norm_comp const nc(dir);
	unsigned const itris_start(itri_verts.size());
	float const css(-1.0*TWO_PI/(float)ndiv), sin_ds(sin(css)), cos_ds(cos(css));
	float sin_s(0.0), cos_s(1.0);
	vector3d const v1(cross_product(dir, (fabs(dir.x) > fabs(dir.y) ? plus_y : plus_x)).get_norm()), v2(cross_product(dir, v1).get_norm());
	itri_verts.emplace_back(pos, nc, 0.5, 0.5, cw);

	for (unsigned i = 0; i < ndiv; ++i) {
		float const s(sin_s), c(cos_s), ts(0.5*(1.0 + (swap_txy ? c : s))), tt(0.5*(1.0 + (swap_txy ? s : c)));
		itri_verts.emplace_back((pos + (radius*s)*v1 + (radius*c)*v2), nc, (inv_ts ? 1.0-ts : ts), (inv_tt ? 1.0-tt : tt), cw);
		indices.push_back(itris_start); // center
		indices.push_back(itris_start + i + 1);
		indices.push_back(itris_start + ((i+1)%ndiv) + 1);
		sin_s = s*cos_ds + c*sin_ds;
		cos_s = c*cos_ds - s*sin_ds;
	}
}

// Note: intended for untextured materials
void rgeom_mat_t::add_round_rect_to_verts(cube_t const &c, float corner_radius, colorRGBA const &color, bool draw_top, bool draw_bot, bool skip_sides,
	bool two_sided, float rs_bot, float rs_top, float side_tscale, float end_tscale, unsigned ndiv)
{
	assert(corner_radius > 0.0);
	cube_t c_inner(c);
	c_inner.expand_by_xy(-corner_radius);
	assert(c_inner.is_strictly_normalized());
	vector2d const sz(c_inner.get_size_xy());
	// start with a vertical cylinder at the cube LLC
	point const llc(c_inner.get_llc());
	unsigned itris_start(itri_verts.size());
	cube_t corner(llc);
	corner.z2() = c.z2();
	corner.expand_by_xy(corner_radius);
	add_vcylin_to_verts(corner, color, draw_bot, draw_top, two_sided, 0, rs_bot, rs_top, side_tscale, end_tscale, skip_sides);
	
	// stretch the curved sides over the cube shape; normals remain unchanged
	for (auto v = itri_verts.begin()+itris_start; v != itri_verts.end(); ++v) {
		for (unsigned d = 0; d < 2; ++d) {
			if      (v->v[d] >  llc[d]) {v->v[d] +=     sz[d];} // stretch this vertex
			else if (v->v[d] == llc[d]) {v->v[d] += 0.5*sz[d];} // center point; half stretch
		}
	}
}

// Note: size can be nonuniform in X/Y/Z
void rgeom_mat_t::add_sphere_to_verts(point const &center, vector3d const &size, colorRGBA const &color, bool low_detail,
	vector3d const &skip_hemi_dir, tex_range_t const &tr, xform_matrix const *const matrix, float ts_add, float tt_add)
{
	static vector<vert_norm_tc>      cached_verts[2]; // high/low detail, reused across all calls
	static vector<vert_norm_comp_tc> cached_vncs [2];
	static vector<unsigned>          cached_ixs  [2];
	vector<vert_norm_tc>      &verts(cached_verts[low_detail]);
	vector<vert_norm_comp_tc> &vncs (cached_vncs [low_detail]);
	vector<unsigned>          &ixs  (cached_ixs  [low_detail]);
	unsigned const ndiv(get_rgeom_sphere_ndiv(low_detail));

	if (verts.empty()) { // not yet created, create and cache verts
		sd_sphere_d sd(all_zeros, 1.0, ndiv);
		sphere_point_norm spn;
		sd.gen_points_norms(spn);
		sd.get_quad_points(verts, &ixs);
		assert((ixs.size()&3) == 0); // must be a multiple of 4
		vncs.resize(verts.size());
		for (unsigned i = 0; i < verts.size(); ++i) {vncs[i] = vert_norm_comp_tc(verts[i].v, verts[i].n, verts[i].t[0], verts[i].t[1]);} // vntc => vnctc
	}
	color_wrapper const cw(color);
	unsigned const ioff(itri_verts.size());
	float const tscale[2] = {(tr.x2 - tr.x1), (tr.y2 - tr.y1)}; // scale existing [0.0, 1.0] texture coords into the specified range

	if (matrix) { // must apply matrix transform to verts and normals and reconstruct norm_comps
		for (auto i = verts.begin(); i != verts.end(); ++i) {
			point pt(i->v*size);
			vector3d normal(i->n);
			matrix->apply_to_vector3d(pt); matrix->apply_to_vector3d(normal);
			itri_verts.emplace_back((pt + center), normal, (tr.x1 + i->t[0]*tscale[0] + ts_add), (tr.y1 + i->t[1]*tscale[1] + tt_add), cw);
		}
	}
	else if (skip_hemi_dir != zero_vector) { // only draw one hemisphere; assumes skip_hemi_dir is along a primary {x, y, z} axis
		// spheres are generated in a circle in the XY plane in the outer loop, and from top to bottom in Z in the inner loop;
		// to draw the top half we need to draw/include the "equator" plus the first/top half of each circular band
		unsigned const stride(ndiv+1), t_end(ndiv/2), out_stride(t_end + 1), indices_start(indices.size());
		bool const inv_hemi(skip_hemi_dir == plus_x || skip_hemi_dir == plus_y || skip_hemi_dir == plus_z);
		bool const swap_x(skip_hemi_dir == plus_x || skip_hemi_dir == -plus_x), swap_y(skip_hemi_dir == plus_y || skip_hemi_dir == -plus_y);
		assert(vncs.size() == stride*stride);

		for (unsigned s = 0; s <= ndiv; ++s) { // XY circular band
			for (unsigned t = 0; t <= t_end; ++t) { // vertical slices in Z; start in the middle; assumes ndiv is an even number
				auto v(vncs[s*stride + t]); // deep copy
				unsigned const ix(itri_verts.size());
				// rotate into correct orientation by swapping coordinates and mirroring; must apply to both the vertex and the normal
				if (inv_hemi) {v.v.z = -v.v.z; v.n[2] = 255 - v.n[2];}
				if (swap_x  ) {std::swap(v.v.x, v.v.z); std::swap(v.n[0], v.n[2]);}
				if (swap_y  ) {std::swap(v.v.y, v.v.z); std::swap(v.n[1], v.n[2]);}
				if (swap_x || swap_y) {v.v.z = -v.v.z; v.n[2] = 255 - v.n[2];}
				itri_verts.emplace_back((v.v*size + center), v, (tr.x1 + v.t[0]*tscale[0] + ts_add), (tr.y1 + v.t[1]*tscale[1] + tt_add), cw);
				if (t == t_end || s == ndiv) continue; // no indices added for last s or t values
				unsigned const qixs[4] = {ix, ix+out_stride, ix+out_stride+1, ix+1};
				for (unsigned i = 0; i < 6; ++i) {indices.push_back(qixs[quad_to_tris_ixs[i]]);} // quads (2 triangles)
			} // for t
		} // for s
		assert(indices.back() < itri_verts.size());
		if (inv_hemi) {std::reverse(indices.begin()+indices_start, indices.end());} // set correct winding order
		return;
	}
	else { // can use vncs (norm_comps)
		for (auto i = vncs.begin(); i != vncs.end(); ++i) {
			itri_verts.emplace_back((i->v*size + center), *i, (tr.x1 + i->t[0]*tscale[0] + ts_add), (tr.y1 + i->t[1]*tscale[1] + tt_add), cw);
		}
	}
	for (auto i = ixs.begin(); i != ixs.end(); i += 4) { // indices are for quads, but we want triangles, so expand them
		indices.push_back(*(i+0) + ioff); indices.push_back(*(i+1) + ioff); indices.push_back(*(i+2) + ioff);
		indices.push_back(*(i+3) + ioff); indices.push_back(*(i+0) + ioff); indices.push_back(*(i+2) + ioff);
	}
	assert(indices.back() < itri_verts.size());
}

void rgeom_mat_t::add_vert_torus_to_verts(point const &center, float r_inner, float r_outer, colorRGBA const &color,
	float tscale, bool low_detail, int half_or_quarter, float s_offset, unsigned ndivo, unsigned ndivi, float spiral_offset)
{
	unsigned const def_ndiv(get_rgeom_sphere_ndiv(low_detail)); // calculate ndiv if not set
	if (ndivo == 0) {ndivo = def_ndiv;}
	if (ndivi == 0) {ndivi = def_ndiv;}
	unsigned s_end(ndivo), t_end(ndivi), sin_cos_off(0);
	apply_half_or_quarter(half_or_quarter, s_end);
	if (half_or_quarter == 3) {t_end /= 2; sin_cos_off += 3*ndivi/4;} // half of a full circle (+z half)
	bool const is_offset(spiral_offset != 0.0);
	float const ts_tt(tscale/ndivi), ds(TWO_PI/ndivo), cds(cos(ds)), sds(sin(ds));
	vector<float> const &sin_cos(gen_torus_sin_cos_vals(ndivi));
	color_wrapper const cw(color);
	float zval(0.0);
	s_offset *= TWO_PI;
	if (is_offset) {spiral_offset /= (ndivo*r_outer);}

	for (unsigned s = 0; s < s_end; ++s) { // outer
		float const theta(s*ds + s_offset), ct(cos(theta)), st(sin(theta)), ct2(ct*cds - st*sds), st2(st*cds + ct*sds);
		point const pos [2] = {point(ct, st, zval), point(ct2, st2, (zval + spiral_offset))};
		point const vpos[2] = {(center + pos[0]*r_outer), (center + pos[1]*r_outer)};
		unsigned const tri_ix_start(itri_verts.size()), ixs_start(indices.size());

		// Note: drawn as one triangle strip
		for (unsigned t = 0; t <= t_end; ++t) { // inner
			unsigned const t_((t + sin_cos_off) % ndivi);
			float const cp(sin_cos[(t_<<1)+0]), sp(sin_cos[(t_<<1)+1]);

			for (unsigned i = 0; i < 2; ++i) {
				vector3d delta(pos[1-i]*sp); // normal
				delta.z += cp;
				if (is_offset) {delta.normalize();}
				itri_verts.emplace_back((vpos[1-i] + delta*r_inner), delta, ts_tt*(s+1-i), ts_tt*t, cw);
			}
		} // for t
		zval += spiral_offset;
		for (unsigned n = 0; n < 3; ++n) {indices.push_back(tri_ix_start + n);} // first triangle

		for (unsigned n = tri_ix_start+3; n < itri_verts.size(); ++n) { // each vertex after this creates a new triangle
			unsigned const ix1(indices[indices.size()-2]), ix2(indices.back()); // two previous indices
			indices.push_back(ix1);
			indices.push_back(ix2);
			indices.push_back(n); // new triangle index
		}
		// swap the winding order of every other triangle, stepping in triangle pairs
		for (unsigned i = ixs_start; i < indices.size(); i += 6) {std::swap(indices[i+4], indices[i+5]);}
	} // for s
}
void rgeom_mat_t::add_contained_vert_torus_to_verts(cube_t const &c, colorRGBA const &color, float tscale, bool low_detail) { // unused
	float const r_inner(0.5*c.dz()), r_outer(0.25*(c.dx() + c.dy()) - r_inner);
	assert(r_inner > 0.0 && r_outer > 0.0); // cube must be wider than it is tall
	add_vert_torus_to_verts(c.get_cube_center(), r_inner, r_outer, color, tscale, low_detail);
}
void rgeom_mat_t::add_ortho_torus_to_verts(point const &center, float r_inner, float r_outer, unsigned dim, colorRGBA const &color,
	float tscale, bool low_detail, int half_or_quarter, float s_offset, unsigned ndivo, unsigned ndivi, float spiral_offset)
{
	assert(dim < 3);
	unsigned const verts_start(itri_verts.size()), ixs_start(indices.size());
	add_vert_torus_to_verts(all_zeros, r_inner, r_outer, color, tscale, low_detail, half_or_quarter, s_offset, ndivo, ndivi, spiral_offset);
	
	if (dim < 2) { // swap X or Y with Z
		for (auto i = itri_verts.begin()+verts_start; i != itri_verts.end(); ++i) {
			std::swap(i->v[dim], i->v[2]);
			std::swap(i->n[dim], i->n[2]);
		}
		reverse(indices.begin()+ixs_start, indices.end()); // reverse winding order
	}
	for (auto i = itri_verts.begin()+verts_start; i != itri_verts.end(); ++i) {i->v += center;}
}

void rgeom_mat_t::add_triangle_to_verts(point const v[3], colorRGBA const &color, bool two_sided, float tscale) {
	color_wrapper cw(color);
	norm_comp normal(get_poly_norm(v));
	float const ts[3] = {0.0, 0.0, tscale}, tt[3] = {0.0, tscale, 0.0}; // hard-coded for now, maybe pass in?

	for (unsigned side = 0; side < 2; ++side) {
		for (unsigned n = 0; n < 3; ++n) {
			indices.push_back(itri_verts.size()); // since we only support indexed triangles, we have to assign each vertex its own index
			unsigned const ix(side ? 2-n : n); // reverse order for side=1
			itri_verts.emplace_back(v[ix], normal, ts[ix], tt[ix], cw);
		}
		if (side == 0) {
			if (!two_sided) break; // single sided
			normal.invert_normal();
		}
	} // for side
}
void rgeom_mat_t::add_quad_to_verts(point const v[4], colorRGBA const &color, float tscale) { // 4 points must be planar
	color_wrapper cw(color);
	norm_comp normal(get_poly_norm(v));
	unsigned const vix(itri_verts.size());
	float const ts[4] = {0.0, tscale, tscale, 0.0}, tt[4] = {0.0, 0.0, tscale, tscale}; // hard-coded for now, maybe pass in?
	for (unsigned n = 0; n < 4; ++n) {itri_verts.emplace_back(v[n], normal, ts[n], tt[n], cw);}
	for (unsigned n = 0; n < 6; ++n) {indices.push_back(vix + quad_to_tris_ixs[n]);}
}

