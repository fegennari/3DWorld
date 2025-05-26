// 3D World - Building Basement Machine Creation and Drawing
// by Frank Gennari 10/15/2024

#include "function_registry.h"
#include "buildings.h"


int get_cylin_duct_tid();
int get_cube_duct_tid ();
int get_ac_unit_tid   (unsigned ix);
int get_metal_texture (unsigned id);
tid_nm_pair_t get_metal_plate_tex(float tscale, bool shadowed);
colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c);
float get_merged_pipe_radius(float r1, float r2, float exponent);


colorRGBA choose_pipe_color(rand_gen_t &rgen) {
	unsigned const NCOLORS = 6;
	colorRGBA const colors[NCOLORS] = {DARK_BRASS_C, COPPER_C, BRONZE_C, DK_BROWN, BROWN, LT_GRAY};
	return colors[rgen.rand() % NCOLORS];
}
colorRGBA choose_machine_part_color(rand_gen_t &rgen, bool is_textured) { // shade of gray
	float const lum(is_textured ? rgen.rand_uniform(0.5, 1.0) : rgen.rand_uniform(0.1, 0.6));
	return colorRGBA(lum, lum, lum);
}

void draw_metal_handle_wheel(cube_t const &c, unsigned dim, colorRGBA const &color, colorRGBA const &shaft_color, rgeom_mat_t &mat, rgeom_mat_t &shaft_mat) {
	float const radius(0.5*c.get_sz_dim((dim+1)%3));
	// draw the outer handle
	float const r_inner(0.12*radius), r_outer(radius - r_inner), r_bar(0.075*radius), r_shaft(0.1*radius);
	point const center(c.get_cube_center());
	mat.add_ortho_torus_to_verts(center, r_inner, r_outer, dim, color);
	// draw inner sphere for handle
	mat.add_sphere_to_verts(center, 0.2*radius, color, 1); // low_detail=1
	// draw horizontal and vertical bars
	unsigned const dims[2] = {(dim+1)%3, (dim+2)%3};
	unsigned const verts_start(mat.itri_verts.size());
	cube_t bar;
	set_wall_width(bar, center[dim], r_bar, dim);

	for (unsigned d = 0; d < 2; ++d) {
		set_wall_width(bar, center[dims[ d]], r_bar,   dims[ d]);
		set_wall_width(bar, center[dims[!d]], r_outer, dims[!d]);
		mat.add_ortho_cylin_to_verts(bar, color, dims[!d], 0, 0); // draw sides only
	}
	// rotate a random-ish amount
	float const rot_angle((c.x1() + c.y1() + c.z1())/radius);
	rotate_verts(mat.itri_verts, vector_from_dim_dir(dim, 1), rot_angle, center, verts_start);
	// draw the shaft
	cube_t shaft(c);
	for (unsigned d = 0; d < 2; ++d) {set_wall_width(shaft, center[dims[d]], r_shaft, dims[d]);}
	shaft_mat.add_ortho_cylin_to_verts(shaft, shaft_color, dim, 1, 1); // draw sides and ends
}

void select_pipe_location(point &p1, point &p2, cube_t const &region, float radius, unsigned dim, rand_gen_t &rgen) {
	float const edge_space(1.5*radius);
	p1[dim] = region.d[dim][0];
	p2[dim] = region.d[dim][1];

	for (unsigned n = 0; n < 2; ++n) {
		unsigned const d((dim+n+1)%3);
		p1[d] = p2[d] = rgen.rand_uniform(region.d[d][0]+edge_space, region.d[d][1]-edge_space);
	}
}

float get_cylin_radius(vector3d const &sz, unsigned dim) {
	return 0.25*(sz[(dim+1)%3] + sz[(dim+2)%3]);
}
void building_room_geom_t::add_machine_pipe_in_region(room_object_t const &c, cube_t const &region, float rmax, unsigned dim,
	vect_sphere_t &pipe_ends, rand_gen_t &rgen, bool allow_valve, bool valve_on_cylin, bool add_coil, bool is_high)
{
	assert(region.is_strictly_normalized());
	min_eq(rmax, get_cylin_radius(region.get_size(), dim));
	float const radius(rgen.rand_uniform((add_coil ? 0.4 : 0.25), 1.0)*rmax);
	point p1, p2;
	select_pipe_location(p1, p2, region, radius, dim, rgen);

	for (sphere_t const &pe : pipe_ends) {
		if (dist_less_than(p1, pe.pos, (radius + pe.radius))) return; // too close to a previously placed pipe
	}
	pipe_ends.emplace_back(p1, radius);
	colorRGBA const color(choose_pipe_color(rgen)), spec_color(get_specular_color(color)); // special case metals

	if (add_coil) {
		bool const sparse(c.in_factory()); // larger coil gap for factories
		float const coil_radius(rgen.rand_uniform(0.05, 0.15)*radius);
		float const coil_gap(rgen.rand_uniform((sparse ? 2.5 : 1.5), (sparse ? 5.0 : 4.0))*coil_radius);
		p1[dim] -= coil_radius; p2[dim] += coil_radius; // must start and end inside the object
		float const length(p2[dim] - p1[dim]);
		add_spring(p1, radius, coil_radius, length, coil_gap, dim, color, spec_color, sparse);
		return;
	}
	rgeom_mat_t &pipe_mat(get_metal_material(1, 0, 1, 0, spec_color));
	pipe_mat.add_cylin_to_verts(p1, p2, radius, radius, apply_light_color(c, color), 0, 0, 0, 0, 1.0, 1.0, 0, 16); // shadowed, small

	// maybe add a valve or gauge
	if (!allow_valve || rgen.rand_float() > 0.6) return; // no valve or gauge; added 60% of the time
	// valve: vertical facing up for horizontal pipes, horizontal facing away from the machine for vertical pipes; gauge: facing X/Y
	bool const is_gauge(rgen.rand_bool());
	float const pipe_len(p2[dim] - p1[dim]), vg_radius(radius*(is_gauge ? rgen.rand_uniform(1.0, 2.4) : rgen.rand_uniform(2.4, 3.6)));
	assert(pipe_len > 0.0);
	if (pipe_len < 4.0*vg_radius) return; // too short to add a valve or gauge
	point attach_pt(p1);
	float const end_pad(max((valve_on_cylin ? 0.4 : 0.25)*pipe_len, 1.5*vg_radius)); // place closer to the middle if part is a cylinder
	float const vg_pos(rgen.rand_uniform((p1[dim] + end_pad), (p2[dim] - end_pad)));
	float const fitting_radius(1.25*radius), fitting_hlen(1.6*radius);
	unsigned vdim(0), odim(0);
	bool vdir(0);

	if (dim == 2) { // vertical pipe
		vector3d const pipe_dir(p1 - c.get_cube_center());
		vdim = (fabs(pipe_dir.x) < fabs(pipe_dir.y));
		vdir = ((c.dim == bool(vdim)) ? c.dir : (pipe_dir[vdim] > 0.0)); // toward the front, or the closer side
		odim = 1 - vdim;
	}
	else { // horizontal pipe
		if (!is_gauge && rgen.rand_bool()) { // vertical valve
			vdim = 2;
			vdir = !is_high; // if horizontal, place below high pipes and above low pipes
			odim = (1 - dim);
		}
		else { // horizontal valve or gauge
			vdim = (1 - dim);
			vdir = ((c.dim == bool(vdim)) ? c.dir : rgen.rand_bool()); // toward the front, or a random side
			odim = 2;
		}
	}
	bool const is_top_gauge(is_gauge && dim < 2 && rgen.rand_bool()); // horizontal pipes only
	if (is_top_gauge) {attach_pt.z += (fitting_radius + 1.5*vg_radius);} // shift up
	else {attach_pt[vdim] += (vdir ? 1.0 : -1.0)*radius;} // shift toward the front

	for (auto pe = pipe_ends.begin(); pe+1 != pipe_ends.end(); ++pe) { // skip pipe end added above
		if (dist_less_than(attach_pt, pe->pos, (vg_radius + pe->radius))) return; // too close to a previously placed pipe
	}
	pipe_ends.emplace_back(attach_pt, vg_radius);
	cube_t vg;
	set_wall_width(vg, vg_pos,          vg_radius,  dim);
	set_wall_width(vg, attach_pt[odim], vg_radius, odim); // pipe centerline
	colorRGBA const fitting_color(is_known_metal_color(color) ? BRASS_C : color); // brass if metal, otherwise keep the pipe color

	if (is_gauge) { // gauge
		set_wall_width(vg, attach_pt[vdim], max(1.1f*(fitting_radius - radius), vg_radius*rgen.rand_uniform(0.2, 0.25)), vdim); // set depth
		unsigned flags(RO_FLAG_NOCOLL | (is_top_gauge ? RO_FLAG_ADJ_BOT : 0));
		add_gauge(room_object_t(vg, TYPE_GAUGE, c.room_id, vdim, vdir, flags, c.light_amt, SHAPE_CYLIN, fitting_color)); // same color as fitting
	}
	else { // valve
		vg.d[vdim][!vdir] = attach_pt[vdim];
		vg.d[vdim][ vdir] = attach_pt[vdim] + (vdir ? 1.0 : -1.0)*vg_radius*rgen.rand_uniform(0.4, 0.5); // set depth
		colorRGBA const handle_color(handle_colors[rgen.rand() % NUM_HANDLE_COLORS]), shaft_color(choose_pipe_color(rgen));
		rgeom_mat_t &handle_mat(get_metal_material(1, 0, 1)); // shadowed
		draw_metal_handle_wheel(vg, vdim, apply_light_color(c, handle_color), apply_light_color(c, shaft_color), handle_mat, handle_mat);
	}
	// draw the fitting
	rgeom_mat_t &fitting_mat(get_metal_material(0, 0, 1, 0, get_specular_color(fitting_color))); // unshadowed
	p1[dim] = vg_pos - fitting_hlen;
	p2[dim] = vg_pos + fitting_hlen;
	fitting_mat.add_cylin_to_verts(p1, p2, fitting_radius, fitting_radius, apply_light_color(c, fitting_color), 1, 1, 0, 0, 1.0, 1.0, 0, N_CYL_SIDES); // draw ends
}

void clip_cylin_to_square(cube_t &c, unsigned dim, bool clip_to_z2=0) { // about the center
	unsigned const d1((dim+1)%3), d2((dim+2)%3); // the non-cylin long dims
	float const da(c.get_sz_dim(d1)), db(c.get_sz_dim(d2)), orig_z2(c.z2());
	if (da < db) {c.expand_in_dim(d2, -0.5*(db - da));} else {c.expand_in_dim(d1, -0.5*(da - db));} // shrink
	if (clip_to_z2 && dim != 2) {c.translate_dim(2, (orig_z2 - c.z2()));} // shift upward to original zval
}

cube_t place_obj_on_cube_side(cube_t const &c, bool dim, bool dir, float hheight, float hwidth, float depth, float pad, rand_gen_t &rgen) {
	float const wpad(pad*hwidth), hpad(pad*hheight);
	cube_t obj;
	obj.d[dim][!dir] = c.d[dim][dir]; // edge of part
	obj.d[dim][ dir] = c.d[dim][dir] + (dir ? 1.0 : -1.0)*depth; // set depth
	set_wall_width(obj, rgen.rand_uniform(c.d[!dim][0]+wpad, c.d[!dim][1]-wpad), hwidth, !dim);
	set_wall_width(obj, rgen.rand_uniform(c.d[2   ][0]+hpad, c.d[2   ][1]-hpad), hheight, 2  );
	return obj;
}

tid_nm_pair_t get_machine_part_texture(bool is_cylin, vector3d const sz, float &tscale, bool untextured, rand_gen_t &rgen) {
	if (untextured) { // for hospitals
		tid_nm_pair_t tex(WHITE_TEX, 0.0, 1); // shadowed
		tex.set_specular(0.5, 40.0);
		return tex;
	}
	if (!is_cylin) {tscale /= sz.get_max_val();} // scale to fit the cube largest dim; cylinder tscale is unused in drawing
	int tid(-1), nm_tid(-1);
	
	if (rgen.rand_float() < 0.1) { // built-in textures; could maybe use shiphull.jpg and get_cube_duct_tid()
		tid     = get_met_plate_tid();
		nm_tid  = get_mplate_nm_tid();
		tscale *= 2.0;
	}
	else {tid = get_metal_texture(rgen.rand());}
	tid_nm_pair_t tex(tid, nm_tid, tscale, tscale, 0.0, 0.0, 1); // shadowed
	tex.set_specular(0.1, 20.0);
	return tex;
}
vector3d calc_cylin_tscales(float tscale, vector3d const &sz, unsigned dim) { // returns {side_tscale, end_tscale, len_tscale}
	float const length(sz[dim]), radius(get_cylin_radius(sz, dim)), diameter(2.0*radius), circumference(PI*radius);
	assert(length > 0.0 && radius > 0.0);
	float const len_tscale(tscale), side_tscale(len_tscale*circumference/length), end_tscale(len_tscale*diameter/length);
	return vector3d(max(1, round_fp(side_tscale)), end_tscale, len_tscale); // round side_tscale to an exact multiple to ensure it tiles
}

void add_pipe_with_bend(rgeom_mat_t &mat, colorRGBA const &color, point const &bot_pt, point const &top_pt, point const &bend, unsigned ndiv, float radius, bool draw_ends) {
	mat.add_sphere_to_verts(bend, vector3d(radius, radius, radius), color, (ndiv == 16), -plus_z); // round part, low detail if ndiv==16, top hemisphere (always bends down)
	mat.add_cylin_to_verts (bot_pt, bend, radius, radius, color, draw_ends, 0, 0, 0, 1.0, 1.0, 0, ndiv); // vertical
	mat.add_cylin_to_verts (top_pt, bend, radius, radius, color, draw_ends, 0, 0, 0, 1.0, 1.0, 0, ndiv); // horizontal
}

void building_room_geom_t::add_spring(point pos, float radius, float r_wire, float length, float coil_gap, unsigned dim,
	colorRGBA const &color, colorRGBA const &spec_color, bool sparse)
{
	assert(dim < 3);
	unsigned num_coils(max(1, round_fp(length/(2.0*r_wire + coil_gap))));
	unsigned const ndivo(N_CYL_SIDES / (sparse ? 2 : 1)), ndivi(N_CYL_SIDES/2);
	float const coil_spacing(length/num_coils);
	rgeom_mat_t &mat(get_metal_material(1, 0, 1, 0, spec_color)); // shadowed, small
	// create one coil, then copy and translate to get the others
	unsigned const verts_start(mat.itri_verts.size()), indices_start(mat.indices.size());
	mat.add_ortho_torus_to_verts(pos, r_wire, radius, dim, color, 1.0, 0, 0, 0.0, ndivo, ndivi, coil_spacing);
	unsigned const verts_end(mat.itri_verts.size()), indices_end(mat.indices.size());

	for (unsigned n = 1; n < num_coils; ++n) {
		unsigned const ix_off(mat.itri_verts.size() - verts_start);
		float const xlate(n*coil_spacing);

		for (unsigned v = verts_start; v != verts_end; ++v) {
			mat.itri_verts.push_back(mat.itri_verts[v]);
			mat.itri_verts.back().v[dim] += xlate;
		}
		for (unsigned i = indices_start; i != indices_end; ++i) {mat.indices.push_back(mat.indices[i] + ix_off);}
	} // for n
}

rand_gen_t get_machine_info(room_object_t const &c, float floor_ceil_gap, cube_t &base, cube_t parts[2], cube_t &support,
	bool is_cylins[2], unsigned cylin_dims[2], unsigned &num_parts, bool &parts_swapped)
{
	rand_gen_t rgen(c.create_rgen());
	float const height(c.dz()), width(c.get_width());
	bool const two_part(width > rgen.rand_uniform(1.2, 1.6)*c.get_depth());
	cube_t main(c);
	num_parts = (two_part ? 2U : 1U);
	base      = c;
	base.z2() = main.z1() = c.z1() + rgen.rand_uniform(0.04, 0.1)*min(height, floor_ceil_gap);
	parts[0] = parts[1] = main;

	if (two_part) { // Note: two_part is not consistent between template and machine
		float const split_pos(main.d[!c.dim][0] +     rgen.rand_uniform( 0.4, 0.6) *width);
		parts[0].d[!c.dim][1] = split_pos - max(0.0f, rgen.rand_uniform(-0.1, 0.1))*width; // maybe add a gap
		parts[1].d[!c.dim][0] = split_pos + max(0.0f, rgen.rand_uniform(-0.1, 0.1))*width; // maybe add a gap
		if (rgen.rand_bool()) {swap(parts[0], parts[1]); parts_swapped = 1;} // remove any bias toward the left/right
	}
	// calculate size and shape of each part
	for (unsigned n = 0; n < num_parts; ++n) {
		cube_t &part(parts[n]);
		vector3d const psz(part.get_size());
		// fewer cylinders in factories; don't place two cylinders
		bool const is_cylin(!is_cylins[0] && rgen.rand_float() < (c.in_factory() ? 0.2 : 0.5) && max(psz.x, psz.y) < 1.5*min(psz.x, psz.y));
		bool const rb1(rgen.rand_bool()), rb2(rgen.rand_bool()); // make rgen independent of size

		if (is_cylin) { // 50% chance vert/Z, 25% chance X, 25% chance Y
			bool const must_be_vert(psz.z > 2.0*min(psz.x, psz.y)); // tall thin cylinders must be vertical
			cylin_dims[n] = ((must_be_vert || rb1) ? 2 : rb2);
		}
		is_cylins[n] = is_cylin;
		float const max_shrink_val(is_cylin ? 0.15 : 0.25);
		for (unsigned d = 0; d < 2; ++d) {part.expand_in_dim(d, -rgen.rand_uniform(0.1, 1.0)*max_shrink_val*part.get_sz_dim(d));} // shrink in both dims independently
		part.z2() -= rgen.rand_uniform(0.05, 0.2)*height; // make a bit shorter

		if (is_cylin) { // make it square
			clip_cylin_to_square(part, cylin_dims[n], 1); // clip_to_z2=1

			if (two_part && part.intersects(parts[1-n])) { // if adjacent, shrink to force a gap; maybe can't happen now?
				part.expand_by_xy(-0.1*min(part.dx(), part.dy()));
				clip_cylin_to_square(part, cylin_dims[n]);
			}
			if (cylin_dims[n] < 2 && part.zc() > main.z1()) { // horizontal cylinder - add a supporting cube (should be always)
				support = part;
				set_cube_zvals(support, main.z1(), part.zc());
				for (unsigned d = 0; d < 2; ++d) {support.expand_in_dim(d, -rgen.rand_uniform(0.05, 0.45)*part.get_sz_dim(d));} // shrink in both dims
			}
		}
	} // for n
	return rgen; // return for use in drawing
}
unsigned get_machine_part_cubes(room_object_t const &c, float floor_ceil_gap, cube_t cubes[4]) { // base part0 [part1] [support]
	cube_t support;
	bool parts_swapped(0), is_cylins [2] = {0, 0};
	unsigned num_parts(0), cylin_dims[2] = {2, 2}; // defaults to Z
	get_machine_info(c, floor_ceil_gap, cubes[0], cubes+1, support, is_cylins, cylin_dims, num_parts, parts_swapped);
	unsigned ncubes(num_parts + 1);
	if (!support.is_all_zeros()) {cubes[ncubes++] = support;}
	return ncubes;
}

void building_room_geom_t::add_machine(room_object_t const &c, float floor_ceil_gap, bldg_industrial_info_t const *ind_info) { // components are shadowed and small
	// can use AC Unit, metal plate texture, buttons, lights, etc.
	float const height(c.dz()), width(c.get_width()), depth(c.get_depth()), orig_floor_ceil_gap(floor_ceil_gap);
	cube_t base, parts[2], support;
	bool parts_swapped(0), is_cylins[2] = {0, 0};
	unsigned num_parts(0), cylin_dims[2] = {2, 2}; // defaults to Z
	rand_gen_t rgen(get_machine_info(c, floor_ceil_gap, base, parts, support, is_cylins, cylin_dims, num_parts, parts_swapped));
	bool const dim(c.dim), dir(c.dir), two_part(num_parts == 2), in_factory(c.in_factory()), in_box(c.was_expanded());
	float const pipe_rmax(0.033*min(height, min(width, depth))), back_wall_pos(c.d[dim][!dir]);
	bool const metal_base( in_factory && (c.flags & RO_FLAG_FROM_SET  )); // factory machine grid
	bool const untextured(!in_factory && (c.flags & RO_FLAG_UNTEXTURED)); // for hospitals
	unsigned base_skip_faces(EF_Z1); // skip bottom
	if (!in_factory) {base_skip_faces |= ~get_face_mask(c.dim, !c.dir);} // skip back face of base against wall for non-factory machines
	int const base_tid(metal_base ? get_texture_by_name("metals/249_iron_metal_plate.jpg") : get_concrete_tid());
	rgeom_mat_t &base_mat(get_material(tid_nm_pair_t(base_tid, 3.0/(width + depth)), 1, 0, 1)); // shadowed, small
	base_mat.add_cube_to_verts(base, apply_light_color(c, (metal_base ? WHITE : c.color)), all_zeros, base_skip_faces);
	vect_cube_t avoid, shapes;
	vect_sphere_t pipe_ends;
	if (in_factory && ind_info) {floor_ceil_gap = ind_info->floor_space.dz();}
	cube_t main(c);
	main.z1() = base.z2();

	// add/draw shapes/objects for each part
	for (unsigned n = 0; n < num_parts; ++n) {
		cube_t const &part(parts[n]);
		vector3d const part_sz(part.get_size());
		float const size_scale(min(1.0f, orig_floor_ceil_gap/part_sz.z)); // not too large for factory machines
		bool const is_cylin(is_cylins[n]);
		unsigned const cylin_dim(cylin_dims[n]);
		avoid .clear();
		shapes.clear();

		// draw main part and horizontal cylinder support
		if (1) { // create a new scope for local variables
			float tscale(1.0);
			tid_nm_pair_t const part_tex(get_machine_part_texture(is_cylin, part_sz, tscale, untextured, rgen));
			rgeom_mat_t &part_mat(get_material(part_tex, 1, 0, 1)); // shadowed, small
			colorRGBA const part_color(apply_light_color(c, choose_machine_part_color(rgen, (part_tex.tid >= 0))));
			
			if (is_cylin) {
				vector3d const ts(calc_cylin_tscales(tscale, part_sz, cylin_dim)); // {side_tscale, end_tscale, len_tscale}
				part_mat.add_ortho_cylin_to_verts(part, part_color, cylin_dim, (cylin_dim < 2), 1, 0, 0, 1.0, 1.0, ts.x, ts.y, 0, 32, 0.0, 0, ts.z); // draw sides and top

				if (!support.is_all_zeros()) { // add/draw cylinder support; should only be for the first part
					tid_nm_pair_t const tex(get_machine_part_texture(0, support.get_size(), tscale, untextured, rgen)); // is_cylin=0
					rgeom_mat_t &support_mat(get_material(tex, 1, 0, 1)); // shadowed, small
					colorRGBA const support_color(apply_light_color(c, choose_machine_part_color(rgen, (part_tex.tid >= 0))));
					support_mat.add_cube_to_verts(support, support_color, part.get_llc(), EF_Z12); // skip top and bottom
				}
			}
			else {part_mat.add_cube_to_verts(part, part_color, part.get_llc(), EF_Z1);} // skip bottom
		}
		// maybe add a breaker panel if a cube
		if (!is_cylin && rgen.rand_float() < 0.6) {
			unsigned const flags((!in_factory && !in_box && rgen.rand_float() < 0.2) ? RO_FLAG_OPEN : 0); // open 20% of the time; not for factory or boxed machines
			float panel_hheight(size_scale*part_sz.z*rgen.rand_uniform(0.15, 0.22)), panel_hwidth(size_scale*part_sz[!dim]*rgen.rand_uniform(0.15, 0.22));
			min_eq(panel_hheight, 1.5f*panel_hwidth ); // set reasonable aspect ratio
			min_eq(panel_hwidth,  1.5f*panel_hheight); // set reasonable aspect ratio
			float const panel_depth(min(panel_hheight, panel_hwidth)*rgen.rand_uniform(0.3, 0.4));
			cube_t panel(place_obj_on_cube_side(part, dim, dir, panel_hheight, panel_hwidth, panel_depth, 1.0, rgen)); // Note: may extend a bit outside c in the front
			assert(panel.is_strictly_normalized());

			if (!has_bcube_int(panel, avoid)) { // Note: checking avoid here in case we add some object before the panel later
				colorRGBA const panel_color(0.3, 0.4, 0.4); // hard-coded, for now, darker than normal breaker panels
				room_object_t const bp(panel, TYPE_BRK_PANEL, c.room_id, dim, !dir, flags, c.light_amt, SHAPE_CUBE, panel_color);
				add_breaker_panel(bp);
				avoid.push_back(panel);

				if (bp.is_open()) { // draw breakers inside as a texured quad
					avoid.back().expand_in_dim(!dim, 0.1*panel_hwidth); // add space for the open door
					rgeom_mat_t &breaker_mat(get_material(tid_nm_pair_t(get_texture_by_name("interiors/breaker_panel.jpg"), 0.0), 0, 0, 1)); // unshadowed, small
					float const border(0.1*min(panel_hwidth, panel_hheight)), new_hheight(panel_hheight - border), new_hwidth(panel_hwidth - border);
					panel.expand_in_dim(!dim, -border); // shrink
					panel.expand_in_dim(2,    -border);
					float const target_ar(1.25), ar(new_hwidth/new_hheight);
					if (ar < target_ar) {set_wall_width(panel, panel.zc(),  new_hwidth /target_ar,  2  );} // too tall
					else {set_wall_width(panel, panel.get_center_dim(!dim), new_hheight*target_ar, !dim);} // too wide
					breaker_mat.add_cube_to_verts(panel, apply_light_color(c, WHITE), part.get_llc(), get_face_mask(dim, dir));
				}
			}
		}
		// maybe add a valve handle to the front if a cube
		if (!is_cylin && rgen.rand_float() < 0.75) {
			float const valve_radius(size_scale*min(5.0f*pipe_rmax, 0.4f*min(part_sz[!dim], part_sz.z))*rgen.rand_uniform(0.75, 1.0));
			float const valve_depth(valve_radius*rgen.rand_uniform(0.4, 0.5));
			cube_t lower_part(part);
			if (in_factory) {min_eq(lower_part.z2(), (part.z1() + max(2.0f*valve_radius, 0.15f*floor_ceil_gap)));} // not too high for factory machines
			cube_t valve(place_obj_on_cube_side(lower_part, dim, dir, valve_radius, valve_radius, valve_depth, 1.1, rgen)); // Note: may extend a bit outside c in the front
			assert(valve.is_strictly_normalized());

			if (!has_bcube_int(valve, avoid)) {
				unsigned const NUM_HANDLE_COLORS = 5;
				colorRGBA const handle_colors[NUM_HANDLE_COLORS] = {colorRGBA(0.5, 0.0, 0.0), DK_RED, colorRGBA(0.5, 0.25, 0.0), LT_GRAY, BKGRAY};
				colorRGBA const handle_color(handle_colors[rgen.rand() % NUM_HANDLE_COLORS]), shaft_color(choose_pipe_color(rgen));
				rgeom_mat_t &mat(get_metal_material(1, 0, 1)); // shadowed, small, specular metal
				// use the same material for the handle and the shaft
				draw_metal_handle_wheel(valve, dim, apply_light_color(c, handle_color), apply_light_color(c, shaft_color), mat, mat);
				avoid.push_back(valve);
			}
		}
		// calculate exterior orients, used for placing vents, AC units, and pipes
		unsigned orients[3]={};
		orients[0] = 2*dim + dir; // front side

		if (two_part) { // side not facing the other part
			orients[1] = 2*(!dim) + (bool(n)^parts_swapped);
		}
		else { // both sides
			for (unsigned d = 0; d < 2; ++d) {orients[d+1] = 2*(!dim) + d;}
		}
		// add vents or AC unit fans
		if (!is_cylin && rgen.rand_float() < 0.9) {
			for (unsigned m = 0; m < (two_part ? 2U : 3U); ++m) {
				if (rgen.rand_float() > (in_factory ? 0.75 : 0.5)) continue;
				bool const vdim(orients[m] >> 1), vdir(orients[m] & 1);
				bool const add_ac_unit(rgen.rand_float() < ((vdim == dim) ? 0.75 : 0.25)), mirror_y(add_ac_unit);
				float const pwidth(part_sz[!vdim]), pheight(part_sz.z);

				for (unsigned N = 0; N < 10; ++N) { // make 10 attempts to place it
					float const vhwidth(sqrt(size_scale)*(add_ac_unit ? rgen.rand_uniform(0.36, 0.4) : rgen.rand_uniform(0.24, 0.32))*min(pwidth, 1.0f*pheight));
					float const vhheight((add_ac_unit ? 0.7 : rgen.rand_uniform(0.4, 0.6))*vhwidth), vdepth(rgen.rand_uniform(0.02, 0.4)*vhwidth);
					cube_t vent(place_obj_on_cube_side(part, vdim, vdir, vhheight, vhwidth, vdepth, (add_ac_unit ? 1.02 : 1.2), rgen));
					vent.intersect_with_cube(main); // can't extend outside or into the base
					assert(vent.is_strictly_normalized());
					if (has_bcube_int(vent, avoid)) continue;
					int const tid(add_ac_unit ? get_ac_unit_tid(rgen.rand()) : get_texture_by_name("interiors/vent.jpg"));
					room_object_t const vent_obj(vent, TYPE_VENT, c.room_id, vdim, !vdir, 0, c.light_amt);
					add_flat_textured_detail_wall_object(vent_obj, c.color, tid, 0, 0, 0, mirror_y); // side vent; skip_z1_face=0, draw_all_faces=0, detail=0 (small)
					avoid.push_back(vent);
					break; // done/success
				} // for N
			} // for n
		}
		// add pipes if gap is large enough
		float const cylin_radius(0.25*(part_sz[(cylin_dim+1)%3] + part_sz[(cylin_dim+2)%3])), pipe_rmax(0.05*part_sz.get_min_val());
		bool const has_gap(!in_factory && fabs(part.d[dim][dir] - c.d[dim][dir]) > pipe_rmax); // only for interior (basement) machines, not in factories next to windows
		cube_t region(part);
		region.d[dim][ dir] = (is_cylin ? part.get_center_dim(dim) : part.d[dim][!dir]); // center plane of cylinder or back of cube
		region.d[dim][!dir] = back_wall_pos; // wall behind the machine
		assert(region.is_strictly_normalized());

		if (in_factory && ind_info) { // connect to neighbors on low X/Y edges
			float const min_height(get_player_eye_height() + 2.8*pipe_rmax); // add space for pipe/coil

			if (part_sz.z > 1.25*min_height) {
				bool const coil_dim(rgen.rand_bool());

				for (unsigned ndim = 0; ndim < 2; ++ndim) {
					if (!(c.flags & (ndim ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO))) continue; // no neighbor in this dim
					float const row_spacing(ind_info->machine_row_spacing[ndim]);
					cube_t region2(part);
					max_eq(region2.z1(), (c.z1() + min_height));
					region2.d[ndim][1] = (is_cylin ? part.get_center_dim(ndim) : part.d[ndim][0]); // center plane of cylinder or low edge of cube
					region2.d[ndim][0] = (is_cylin ? part.get_center_dim(ndim) : part.d[ndim][1]) - row_spacing; // opposite edge of this part on adjacet machine
					assert(region2.is_strictly_normalized());
					rand_gen_t rgen2(rgen); // copy rgen so that we can make pipe connections differ across machine instances
					rgen2.rseed1 += (ndim ? c.state_flags : c.drawer_flags); // make it different per row or column
					//rgen2.rseed1 += 17*c.state_flags + 31*c.drawer_flags; // make it different per row and column
					unsigned const num_pipes((rgen2.rand() % 4) + 3); // 3-6
					unsigned num_coils(0);
					pipe_ends.clear();

					for (unsigned n = 0; n < num_pipes; ++n) {
						bool const is_coil(!is_cylin && bool(ndim) == coil_dim && num_coils == 0 && rgen2.rand_bool()); // max of 1 coil across both dims; cubes only
						num_coils += is_coil;
						add_machine_pipe_in_region(c, region2, (is_coil ? 2.0 : 1.25)*pipe_rmax, ndim, pipe_ends, rgen2, 1, is_cylin, is_coil, 1); // allow_valve=1, is_high=1
					}
				} // for n
			}
		}
		else if ((is_cylin || has_gap) && cylin_dim != unsigned(dim)) { // if there's a gap between the machine and the wall; not for cylinders facing the wall
			// add pipe(s) connecting to back wall in {dim, !dir} or ceiling
			unsigned const num_pipes((rgen.rand() % 4) + 1); // 1-4
			bool const allow_valve(!is_cylin); // add valves to cubes only so that they don't intersect the parts
			pipe_ends.clear();
			for (unsigned n = 0; n < num_pipes; ++n) {add_machine_pipe_in_region(c, region, pipe_rmax, dim, pipe_ends, rgen, allow_valve);}
		}
		// if there's space and floor_ceil_gap was specified; skip for factories to avoid clipping through lights and beams; skip for boxed machines
		if (!in_factory && !in_box && height < floor_ceil_gap) {
			float const ceil_zval(c.z1() + floor_ceil_gap);
			// add pipes up to the ceiling
			unsigned const num_pipes(rgen.rand() % 4); // 0-3

			if (num_pipes > 0) {
				bool allow_valve(1);
				cube_t region(part);
				set_cube_zvals(region, part.z2(), ceil_zval);
				
				if (is_cylin) {
					if (cylin_dim == 2) {region.expand_by_xy(-(1.0 - SQRTOFTWOINV)*cylin_radius);} // vertical; use square inscribed in circle
					else {region.z1() = part.zc(); allow_valve = 0;} // horizontal; extend to center of cylinder; no valves
				}
				pipe_ends.clear();
				for (unsigned n = 0; n < num_pipes; ++n) {add_machine_pipe_in_region(c, region, pipe_rmax, 2, pipe_ends, rgen, allow_valve);}
			}
			// maybe add a spider web between the machine and the wall; only for cubes and vertical cylinders, and basement machines
			if (!in_factory && (!is_cylin || cylin_dim == 2) && rgen.rand_float() < 0.75) {
				float const web_pos(is_cylin ? part.get_center_dim(!dim) : rgen.rand_uniform(part.d[!dim][0], part.d[!dim][1]));
				room_object_t web(c, TYPE_SPIWEB, c.room_id, !dim, !dir, 0, c.light_amt);
				set_cube_zvals(web, (part.z2() - 0.1*part_sz.z), (ceil_zval + 0.02*floor_ceil_gap));
				set_wall_width(web, web_pos, 0.0, !dim); // zero thickness
				web.d[dim][!dir] = back_wall_pos;
				web.d[dim][ dir] = max(part.d[dim][0], min(part.d[dim][1], (back_wall_pos + (dir ? 1.0f : -1.0f)*web.dz()))); // clamp to part range
				add_spider_web(web);
			}
		}
		// add a duct to the wall
		if (!is_cylin && has_gap && rgen.rand_float() < 0.65) {
			bool const two_sided(in_box); // since end is visible
			float const radius(rgen.rand_uniform(2.0, 4.0)*pipe_rmax); // wider than a pipe
			point p1, p2;
			select_pipe_location(p1, p2, region, radius, dim, rgen);
			rgeom_mat_t &duct_mat(get_material(tid_nm_pair_t(get_cylin_duct_tid(), 1.0, 1), 1, 0, 1)); // shadowed, small
			duct_mat.add_cylin_to_verts(p1, p2, radius, radius, apply_light_color(c, LT_GRAY), 0, 0, two_sided, 0, 1.0, 1.0, 0, 32, 1.0, 1); // ndiv=32, swap_txy=1
			cube_t vent_bc(p1, p2);
			vent_bc.expand_by(radius); // close enough
			avoid.push_back(vent_bc);
		}
		// add smaller shapes to this one
		unsigned num_shapes(rgen.rand() % (is_cylin ? 5 : 9) + (two_part ? 0 : 1)); // 0-4 for cylin, 0-8 for cube, +1 for single part
		if (in_factory) {num_shapes = 2*num_shapes + rgen.rand_bool();} // 2x for factories
		float const max_shape_sz(0.25*part_sz.get_min_val());
		point const part_center(part.get_cube_center());

		for (unsigned N = 0; N < num_shapes; ++N) {
			bool const add_cylin(rgen.rand_float() < (is_cylin ? 0.5 : 0.25));

			for (unsigned m = 0; m < 10; ++m) { // 10 attempts to place this shape
				unsigned const face(rgen.rand() & 3); // {top, front, left, right}
				unsigned const face_dims[4] = {2, dim, !dim, !dim}, face_dirs[4] = {1, dir, 0, 1};
				unsigned const fdim(face_dims[face]), fdir(face_dirs[face]);
				unsigned cdim(2); // defaults to Z
				
				if (add_cylin) { // cylin on part side can be oriented in X, Y, or Z
					if (is_cylin) {
						if (fdim == cylin_dim) {cdim = cylin_dim;} // end attachment is always the same orient
						else {cdim = (rgen.rand_bool() ? fdim : cylin_dim);} // side attachment can be horizontal or vertical
					}
					else {cdim = fdim;} // cylinder always follows surface orient for cubes
				}
				vector3d half_sz;
				cube_t shape;
				for (unsigned d = 0; d < 3; ++d) {half_sz[d] = rgen.rand_uniform(0.4, 1.0)*max_shape_sz;}
				bool const centered(add_cylin && is_cylin && cdim == cylin_dim && rgen.rand_float() < 0.75); // centered cylinders in same dim 75% of the time

				if (add_cylin) { // make sure cylinder ends are square; could use clip_cylin_to_square(), but this always makes shapes smaller
					unsigned const d1((cdim+1)%3), d2((cdim+2)%3); // the non-cylin long dims
					half_sz[d1] = half_sz[d2] = 0.5*(half_sz[d1] + half_sz[d2]);
				}
				for (unsigned d = 0; d < 3; ++d) {
					float val(0.0);
					if (d == fdim) {val = part.d[fdim][fdir];} // attach to surface
					else if (centered || (d != cylin_dim && is_cylin && !(add_cylin && cdim == cylin_dim))) {val = part.get_center_dim(d);} // centered on the side
					else {val = rgen.rand_uniform(part.d[d][0], part.d[d][1]);} // chose a random attachment point
					set_wall_width(shape, val, half_sz[d], d);
				}
				shape.intersect_with_cube(main);
				if (add_cylin) {clip_cylin_to_square(shape, cdim);}
				assert(shape.is_strictly_normalized());
				if (!part.intersects(shape) || part.contains_cube(shape)) continue; // shouldn't happen?
				if (two_part && parts[1-n].intersects(shape)) continue; // safer to check this, though it may be okay to allow
				
				if (add_cylin) { // cylinder case is more difficult
					float const shape_radius(0.25*(shape.dx() + shape.dy()));
					point const shape_center(shape.get_cube_center());
					if (!circle_rect_intersect(shape_center, shape_radius, part, cdim)) continue; // new cylinder vs. part bcube
					
					if (is_cylin && cdim == cylin_dim) { // cylins same dim
						vector3d delta(part_center - shape.get_cube_center());
						delta[cdim] = 0.0; // ignore this dim
						if (delta.mag_sq() >= (cylin_radius + shape_radius)*(cylin_radius + shape_radius)) continue;
					}
				}
				if (is_cylin && !circle_rect_intersect(part_center, cylin_radius, shape, cylin_dim)) continue; // new bcube vs. part cylinder
				if (has_bcube_int(shape, avoid)) continue;
				bool bad_place(0);

				for (cube_t const &s : shapes) { // check previously placed shapes for Z-fighting at clip boundaries
					if (!s.intersects(shape)) continue;
					if (s.x1()==shape.x1() || s.x2()==shape.x2() || s.y1()==shape.y1() || s.y2()==shape.y2() || s.z1()==shape.z1() || s.z2()==shape.z2()) {bad_place = 1; break;}
				}
				if (bad_place) continue;
				float tscale(1.0);
				tid_nm_pair_t const tex(get_machine_part_texture(add_cylin, 2.0*half_sz, tscale, untextured, rgen));
				rgeom_mat_t &mat(get_material(tex, 1, 0, 1)); // shadowed, small
				colorRGBA const pcolor(choose_machine_part_color(rgen, (tex.tid >= 0))), spec_color(get_specular_color(pcolor)), lcolor(apply_light_color(c, pcolor));
				
				if (add_cylin) {
					vector3d const ts(calc_cylin_tscales(tscale, 2.0*half_sz, cdim)); // {side_tscale, end_tscale, len_tscale}
					mat.add_ortho_cylin_to_verts(shape, lcolor, cdim, 1, 1, 0, 0, 1.0, 1.0, ts.x, ts.y, 0, 32, 0.0, 0, ts.z); // draw top/bot
				}
				else {mat.add_cube_to_verts(shape, lcolor, shape.get_llc(), 0);} // draw all faces since we don't track which are visible
				shapes.push_back(shape);
				break; // success/done
			} // for m
		} // for N
		// add pipes down to the floor
		unsigned const num_floor_pipes((rgen.rand() % (in_factory ? 6 : 5)) + in_factory); // 0-4, 1-6 for factories

		if (num_floor_pipes > 0) {
			colorRGBA const color(choose_pipe_color(rgen)), spec_color(get_specular_color(color)), pcolor(apply_light_color(c, color));
			rgeom_mat_t &mat(get_metal_material(1, 0, 1, 0, spec_color));
			bool used_sides[3] = {};

			for (unsigned n = 0; n < num_floor_pipes; ++n) {
				for (unsigned N = 0; N < 4; ++N) { // make 4 placement attempts
					unsigned const orient_ix(rgen.rand() % (two_part ? 2 : 3)); // pick a random exterior side

					if (is_cylin) { // cylinder can have at most one pipe per orient since they're at the same horizontal position
						if (used_sides[orient_ix]) continue;
						used_sides[orient_ix] = 1;
					}
					bool const fdim(orients[orient_ix] >> 1), fdir(orients[orient_ix] & 1);
					float const radius(rgen.rand_uniform(0.25, 1.0)*pipe_rmax), edge_space(1.5*radius), face_edge(part.d[fdim][fdir]);
					cube_t base_valid(base);
					base_valid.expand_by_xy(-edge_space); // shrink
					point top_pos;

					for (unsigned d = 0; d < 3; ++d) {
						if (d == unsigned(fdim)) {top_pos[d] = (is_cylin ? part.get_center_dim(fdim) : face_edge);} // attach to surface, or centerline for cylinder
						else if (is_cylin && (/*d != cylin_dim*/unsigned(fdim) == cylin_dim && d == 2)) {top_pos[d] = part.get_center_dim(d);} // centered on side of cylinder part
						else {top_pos[d] = rgen.rand_uniform(part.d[d][0]+edge_space, part.d[d][1]-edge_space);} // chose a random attachment point
					}
					point bend_pos(top_pos);
					bend_pos[fdim] = face_edge + rgen.rand_uniform(2.0, 5.0)*radius*(fdir ? 1.0 : -1.0); // extend a random amount; may go outside the base
					base_valid.clamp_pt_xy(bend_pos);
					if (fabs(bend_pos[fdim] - face_edge) < 2.0*radius) continue; // too short/too close to the part face
					point bot_pos(bend_pos);
					bot_pos.z = base.z2(); // top of base platform
					cube_t pipe_bcube(top_pos, bot_pos);
					pipe_bcube.expand_by(radius);
					if (has_bcube_int(pipe_bcube, avoid) || has_bcube_int(pipe_bcube, shapes)) continue; // skip if any other object is intersected
					add_pipe_with_bend(mat, pcolor, bot_pos, top_pos, bend_pos, 16, radius, 0); // ndiv=16, draw_ends=0
					avoid.push_back(pipe_bcube); // is this needed?
					break; // done/success
				} // for N
			} // for n
		}
	} // for n (parts)
	if (two_part) {
		// connect the two parts with pipe and coils
		unsigned const num_pipes((rgen.rand() % (in_factory ? 6 : 4)) + 1); // 1-4, 1-6 for factories
		cube_t region(parts[0]), region2(parts[1]);
		float side_pos[2] = {};

		for (unsigned d = 0; d < 2; ++d) {
			bool const end_on_cylin(is_cylins[d] && cylin_dims[d] == unsigned(!dim));
			side_pos[d] = ((is_cylins[d] && !end_on_cylin) ? parts[d].get_center_dim(!dim) : parts[d].d[!dim][bool(d)^parts_swapped^1]);
			if (!end_on_cylin) continue; // not end-on horizontal cylinder
			cube_t &r(d ? region2 : region);
			float const cylin_radius(0.5*r.dz()), exp_val(-(1.0 - SQRTOFTWOINV)*cylin_radius);
			region.expand_in_z  (     exp_val); // clip to inscribed cube
			region.expand_in_dim(dim, exp_val);
		} // for d
		min_eq(region.z2(), region2.z2()); // shared Z range
		max_eq(region.d[dim][0], region2.d[dim][0]); // shared range
		min_eq(region.d[dim][1], region2.d[dim][1]);

		if (region.is_strictly_normalized() && side_pos[0] != side_pos[1]) {
			bool const either_cylin(is_cylins[0] || is_cylins[1]), allow_valve(!either_cylin);
			float const coil_prob(either_cylin ? 0.2 : 0.5); // less likely for cylinders since some of the coils are hidden inside the cylider
			region.d[!dim][0] = min(side_pos[0], side_pos[1]);
			region.d[!dim][1] = max(side_pos[0], side_pos[1]);
			assert(region.is_strictly_normalized());
			unsigned num_coils(0);
			pipe_ends.clear();
			
			for (unsigned n = 0; n < num_pipes; ++n) {
				bool const is_coil(num_coils < (either_cylin ? 1U : 2U) && rgen.rand_float() < coil_prob); // max of 2 coils
				num_coils += is_coil;
				add_machine_pipe_in_region(c, region, (is_coil ? 2.0 : 1.0)*pipe_rmax, !dim, pipe_ends, rgen, allow_valve, either_cylin, is_coil);
			}
			// add flanges if large and connecting to a cube?
		}
	}
}

vector2d get_machine_max_sz(cube_t const &place_area, float min_gap, float max_place_sz, float sz_cap_mult) {
	vector2d const place_sz(place_area.dx(), place_area.dy());
	vector2d avail_sz, max_sz;
	for (unsigned d = 0; d < 2; ++d) {avail_sz[d] = min(max_place_sz, min(0.4f*place_sz[d], sz_cap_mult*(place_sz[d] - min_gap)));}
	for (unsigned d = 0; d < 2; ++d) {max_sz  [d] = min(avail_sz[d], 2.0f*avail_sz[!d]);} // keep aspect ratio <= 2:1
	return max_sz;
}
bool building_t::add_machines_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool less_clearance) {
	float const floor_spacing(get_window_vspace()), fc_gap(get_floor_ceil_gap());
	float const clearance_factor(less_clearance ? 0.2 : 1.0), min_place_sz((less_clearance ? 0.25 : 0.5)*floor_spacing), max_place_sz(1.5*floor_spacing);
	float const min_clearance(clearance_factor*get_min_front_clearance_inc_people()), min_gap(max(clearance_factor*get_doorway_width(), min_clearance));
	unsigned const num_machines((rgen.rand() % 5) + 1); // 1-5
	cube_t const place_area(get_walkable_room_bounds(room)); // ignore trim?
	cube_t avoid;
	avoid.set_from_sphere(point(place_area.xc(), place_area.yc(), zval), min_clearance);
	bool any_placed(0), no_opposite_sides(0), used_orients[4] = {};
	vect_room_object_t &objs(interior->room_geom->objs);
	vector2d max_sz(get_machine_max_sz(place_area, min_gap, max_place_sz, 0.5)); // start with a larger gap that allows two opposing machines

	if (min(max_sz.x, max_sz.y) < min_place_sz) { // not enough space
		max_sz = get_machine_max_sz(place_area, min_gap, max_place_sz, 1.0); // try again with larger size cap, which means we can't place machines on opposite walls
		if (min(max_sz.x, max_sz.y) < min_place_sz) return 0; // still too small of a room to place a machine
		no_opposite_sides = 1;
	}
	for (unsigned n = 0; n < num_machines; ++n) {
		float const height(fc_gap*rgen.rand_uniform(0.6, 0.9));
		// similar to place_obj_along_wall(), but with custom size logic that depends on dim
		cube_t c;
		set_cube_zvals(c, zval, zval+height);

		for (unsigned i = 0; i < 25; ++i) { // make 25 attempts to place the object
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
			if (no_opposite_sides && used_orients[unsigned(dim) + unsigned(!dir)]) continue; // can't place opposing object
			float const hwidth(0.5*max_sz[!dim]*rgen.rand_uniform(0.75, 1.0));
			float const depth(min(2.0f*hwidth, max_sz[dim])*rgen.rand_uniform(0.75, 1.0));
			float center(rgen.rand_uniform(place_area.d[!dim][0]+hwidth, place_area.d[!dim][1]-hwidth)); // random position
			c.d[dim][ dir] = place_area.d[dim][dir];
			c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
			set_wall_width(c, center, hwidth, !dim);
			if (c.intersects(avoid) || overlaps_other_room_obj(c, objs_start) || is_obj_placement_blocked(c, room, 1)) continue; // inc_open_doors=1
			objs.emplace_back(c, TYPE_MACHINE, room_id, dim, !dir, 0, tot_light_amt, SHAPE_CUBE, LT_GRAY, rgen.rand()); // add more randomness
			set_obj_id(objs);
			if (no_opposite_sides) {used_orients[unsigned(dim) + unsigned(dir)] = 1;} // mark this orient as used
			any_placed = 1;
			break; // done
		} // for i
	} // for n
	return any_placed;
}

void building_t::add_machines_to_factory(rand_gen_t rgen, room_t const &room, cube_t const &place_area, float zval,
	unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned objs_start_inc_beams, cube_t const &ladder)
{
	assert(is_factory() && interior->ind_info);
	bool const edim(interior->ind_info->entrance_dim), edir(interior->ind_info->entrance_dir);
	float const floor_spacing(get_window_vspace()), fc_gap(room.dz()), max_place_sz(1.0*floor_spacing), max_height(fc_gap - floor_spacing);
	float const doorway_width(get_doorway_width()), min_gap(max(doorway_width, get_min_front_clearance_inc_people())), ceil_zval(room.z2() - get_fc_thickness());
	if (max_height <= 0.0) return; // should never happen
	bool const check_ladder(!ladder.is_all_zeros());
	vect_room_object_t &objs(interior->room_geom->objs);
	vector2d max_sz(get_machine_max_sz(place_area, min_gap, max_place_sz, 0.5));
	cube_t avoid(interior->ind_info->entrance_area);
	avoid.expand_in_dim(edim, doorway_width); // don't block entrance area
	
	if (1) { // add a 2D grid of machines to the center of the place area
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // used for machines and tanks
		float const height(max_height*rgen.rand_uniform(0.6, 0.8)), aisle_spacing(1.25*doorway_width);
		vector2d const machine_sz(max_place_sz*vector2d(rgen.rand_uniform(0.8, 1.0), rgen.rand_uniform(0.8, 1.0)));
		float const tank_radius(0.5*machine_sz.get_min_val()), pipe_radius(0.04*tank_radius);
		float const tank_height(max(2.25*tank_radius, 0.75*height)), tank_z2(zval + tank_height);
		unsigned const item_flags(rgen.rand()), rand_seed(rgen.rand());
		cube_t center_area(place_area);
		center_area.expand_by_xy(-(max_place_sz + 0.5*doorway_width)); // space for machines along the wall
		vector2d const center_sz(center_area.get_size_xy());
		bool const tank_dim(rgen.rand_bool()), tank_dir(rgen.rand_bool()); // side of grid the tank is on
		colorRGBA const pipe_color(COPPER_C), fitting_color(BRASS_C);
		vector2d &spacing(interior->ind_info->machine_row_spacing);
		unsigned num_xy[2]={};

		for (unsigned d = 0; d < 2; ++d) {
			num_xy [d] = center_sz[d]/(machine_sz[d] + aisle_spacing);
			spacing[d] = center_sz[d]/num_xy[d];
		}
		unsigned machine_range[2][2] = {{0, num_xy[0]-1}, {0, num_xy[1]-1}};
		if (tank_dir) {--machine_range[tank_dim][1];} else {++machine_range[tank_dim][0];} // clip off tank row/column
		vector<room_object> obj_grid(num_xy[0]*num_xy[1], TYPE_NONE);
		vector<point> tank_conn_pts;
		float merged_pipe_radius(0.0);
		cube_t tank_conn_pipe;

		for (unsigned ny = 0; ny < num_xy[1]; ++ny) {
			for (unsigned nx = 0; nx < num_xy[0]; ++nx) {
				bool const is_tank((tank_dim ? ny : nx) == (tank_dir ? num_xy[tank_dim]-1 : 0));
				point const center((center_area.x1() + (nx + 0.5)*spacing.x), (center_area.y1() + (ny + 0.5)*spacing.y), zval);
				cube_t c(center);
				c.z2() += height;
				if (is_tank) {c.expand_by_xy(tank_radius);} // tank is square
				else         {c.expand_by_xy(0.5*machine_sz);}
				if (c.intersects(avoid) || (check_ladder && c.intersects(ladder)) || overlaps_obj_or_placement_blocked(c, room, objs_start)) continue;
				unsigned const gix(ny*num_xy[0] + nx);

				if (is_tank) { // make it a chemical tank; the tank itself is smaller to make room for pipes
					cube_t tank(c);
					tank.z2() = tank_z2;
					tank.expand_by_xy(-0.1*tank_radius);
					room_object_t const tank_obj(tank, TYPE_CHEM_TANK, room_id, dim, dir, RO_FLAG_IN_FACTORY, tot_light_amt, SHAPE_CYLIN, WHITE, item_flags);
					objs.push_back(tank_obj);
					// add a gauge to the side of the tank
					float const gauge_radius(0.08*tank_radius), gauge_height(1.25*tank_radius);
					add_chem_tank_gauge(tank_obj, gauge_radius, gauge_height);
					// add a pipe to the top; the tank itself has a pipe going down to the floor
					cube_t pipe(center);
					set_cube_zvals(pipe, (tank_z2 - pipe_radius), (tank_z2 + 0.15*height - pipe_radius));
					pipe.expand_by_xy(pipe_radius);
					objs.emplace_back(pipe, TYPE_PIPE, room_id, 0, 1, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, pipe_color); // vertical
					tank_conn_pipe.assign_or_union_with_cube(pipe);
					merged_pipe_radius = get_merged_pipe_radius(merged_pipe_radius, pipe_radius, 3.0);
					tank_conn_pts.emplace_back(center.x, center.y, pipe.z2());
				}
				else { // make it a machine
					unsigned flags(RO_FLAG_IN_FACTORY | RO_FLAG_FROM_SET); // tag as part of a group
					cube_t space_x(c), space_y(c);
					space_x.x1() -= spacing.x;
					space_y.y1() -= spacing.y;
					
					if (nx > 0 && obj_grid[gix - 1        ] == TYPE_MACHINE) { // has prev X neighbor	
						if (!check_ladder || !space_x.intersects(ladder)) {flags |= RO_FLAG_ADJ_LO;} // add pipes if not blocked by a ladder
					}
					if (ny > 0 && obj_grid[gix - num_xy[0]] == TYPE_MACHINE) { // has prev Y neighbor
						if (!check_ladder || !space_y.intersects(ladder)) {flags |= RO_FLAG_ADJ_HI;} // add pipes if not blocked by a ladder
					}
					objs.emplace_back(c, TYPE_MACHINE, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, LT_GRAY, item_flags);
					room_object_t &machine(objs.back());
					machine.obj_id       = rand_seed;
					machine.state_flags  = nx;
					machine.drawer_flags = ny;
					bool added_light(0);

					// add corner machine warning lights
					for (unsigned ydir = 0; ydir < 2 && !added_light; ++ydir) {
						if (ny != machine_range[1][ydir]) continue; // not edge row

						for (unsigned xdir = 0; xdir < 2 && !added_light; ++xdir) {
							if (nx != machine_range[0][xdir]) continue; // not edge row
							cube_t cubes[4];
							get_machine_part_cubes(machine, get_floor_ceil_gap(), cubes);
							cube_t const &part(cubes[1]); // first (should be only) part
							assert(part.is_strictly_normalized());
							float const light_height(0.14*height), light_radius(0.12*light_height);
							point center(0.0, 0.0, (part.z1() + min(floor_spacing, 0.4f*part.dz())));

							for (unsigned D = 0; D < 2; ++D) { // extend dim
								for (unsigned d = 0; d < 2; ++d) {
									bool const dir(d ? ydir : xdir), extend_dim(d != D);
									// contained in part range for dim D, outside part in !D; may extend outside machine bcube, but shouldn't intersect anything
									center[d] = (extend_dim ? c : part).d[d][dir] + (dir ? -1.0 : 1.0)*(extend_dim ? -1.5 : 3.0)*light_radius;
								}
								cube_t light(center);
								light.z2() += light_height;
								light.expand_by_xy(light_radius);
								objs.emplace_back(light, TYPE_WARN_LIGHT, room_id, 0, 0, (RO_FLAG_IN_FACTORY | RO_FLAG_NOCOLL | RO_FLAG_LIT), tot_light_amt, SHAPE_CYLIN);
								// add a horizontal rod connecting the light to the machine part
								bool const bar_dir(D ? xdir : ydir);
								cube_t bar(center);
								bar.expand_by(0.2*light_radius);
								bar.d[!D][!bar_dir] = /*part.d[!D][bar_dir]*/part.get_center_dim(!D); // use center in case part is a cylinder
								objs.emplace_back(bar, TYPE_METAL_BAR, room_id, !D, bar_dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, GRAY);
							} // for D
							added_light = 1;
						} // for xdir
					} // for ydir
				}
				obj_grid[gix] = (is_tank ? TYPE_CHEM_TANK : TYPE_MACHINE);
			} // for nx
		} // for ny
		if (merged_pipe_radius > 0.0) { // create horizontal pipe connecting tanks
			float const conn_pipe_center(tank_conn_pipe.get_center_dim(tank_dim));
			set_wall_width(tank_conn_pipe, tank_conn_pipe.z2(), merged_pipe_radius, 2); // set zvals
			set_wall_width(tank_conn_pipe, conn_pipe_center,    merged_pipe_radius, tank_dim);
			// add a bend and another segment to the ceiling, floor, or wall
			bool const pri_ext_dir(rgen.rand_bool());
			float const tank_gap(tank_radius + merged_pipe_radius), wall_gap(merged_pipe_radius + 0.26*floor_spacing); // leave wall gap for duct width
			cube_t h_pipe(tank_conn_pipe), v_pipe;
			unsigned h_pipe_flags(RO_FLAG_NOCOLL | RO_FLAG_LIT);
			bool connected(0), conn_dir(0), conn_up(0);

			for (unsigned n = 0; n < 50 && !connected; ++n) { // 50 attempts to extend the horizontal pipe and connect it vertically
				for (unsigned D = 0; D < 2 && !connected; ++D) { // try to extend both sides
					bool const d(bool(D) ^ pri_ext_dir);
					float const pipe_end(tank_conn_pipe.d[!tank_dim][d]), wall_pos(place_area.d[!tank_dim][d]), dsign(d ? 1.0 : -1.0);
					float cand_pos(0.0);
					if (d == 0) {cand_pos = rgen.rand_uniform((wall_pos + wall_gap), (pipe_end - tank_gap));}
					else        {cand_pos = rgen.rand_uniform((pipe_end + tank_gap), (wall_pos - wall_gap));}
					cube_t ext_pipe(tank_conn_pipe);
					ext_pipe.d[!tank_dim][!d] = pipe_end + dsign*tank_gap; // shift so as not to intersect the tank
					ext_pipe.d[!tank_dim][ d] = cand_pos + dsign*merged_pipe_radius; // extend to outer edge
					if (ext_pipe.intersects(avoid) || overlaps_obj_or_placement_blocked(ext_pipe, room, objs_start)) continue;
					ext_pipe.d[!tank_dim][ d] = cand_pos; // clip off overlapping end so that we can add a bend
					// horizontal pipe is okay, check vertical pipe up and down
					v_pipe = ext_pipe;
					set_wall_width(v_pipe, cand_pos, merged_pipe_radius, !tank_dim);
					bool init_up(rgen.rand_bool());

					for (unsigned e = 0; e < 2 && !connected; ++e) {
						conn_up = (init_up ^ bool(e));
						if (conn_up) {set_cube_zvals(v_pipe, ext_pipe.zc(), ceil_zval);} // up
						else         {set_cube_zvals(v_pipe, zval, ext_pipe.zc());} // down
						if (v_pipe.intersects(avoid) || overlaps_obj_or_placement_blocked(v_pipe, room, (conn_up ? objs_start_inc_beams : objs_start), conn_up)) continue;
						h_pipe.union_with_cube(ext_pipe); // extend pipe
						unsigned const flags(RO_FLAG_LIT | (conn_up ? RO_FLAG_NOCOLL : 0)); // cast shadows; only collide if going down
						objs.emplace_back(v_pipe, TYPE_PIPE,    room_id, 0, 1, flags, tot_light_amt, SHAPE_CYLIN, pipe_color);
						objs.emplace_back(v_pipe, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_NOCOLL, 1.0); // add a blocker over the pipe so that we don't intersect it with sprinklers
						connected = 1; conn_dir = d;
					} // for e
				} // for D
			} // for n
			objs.emplace_back(h_pipe, TYPE_PIPE, room_id, !tank_dim, 0, h_pipe_flags, tot_light_amt, SHAPE_CYLIN, pipe_color); // horizontal
			// add pipe fittings
			unsigned const fittings_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING);
			float const fitting_exp(0.2*merged_pipe_radius), fitting_hlen(1.2*merged_pipe_radius);
			cube_t h_fitting(h_pipe);
			h_fitting.expand_in_dim(tank_dim, fitting_exp);
			h_fitting.expand_in_z(fitting_exp);

			if (connected) {
				// add v_pipe fitting into ceiling or floor
				unsigned vfflags(fittings_flags);
				cube_t fitting(v_pipe);
				fitting.expand_by_xy(fitting_exp);
				if (conn_up) {fitting.z1() = v_pipe.z2() - fitting_hlen; vfflags |= RO_FLAG_ADJ_LO;} // up
				else         {fitting.z2() = v_pipe.z1() + fitting_hlen; vfflags |= RO_FLAG_ADJ_HI;} // down
				objs.emplace_back(fitting, TYPE_PIPE, room_id, 0, 1, vfflags, tot_light_amt, SHAPE_CYLIN, fitting_color); // vertical
				// add a bend
				float const end_pos(h_pipe.d[!tank_dim][conn_dir]), hpipe_zc(h_pipe.zc()), dsign(conn_dir ? 1.0 : -1.0);
				// make this section rounded and shadow casting, since it forms the elbow of the pipe
				unsigned const hfflags(fittings_flags | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI | RO_FLAG_LIT | (conn_dir ? RO_FLAG_ADJ_TOP : RO_FLAG_ADJ_BOT));
				fitting = h_fitting;
				fitting.d[!tank_dim][ conn_dir] = end_pos; // flush with pipe end
				fitting.d[!tank_dim][!conn_dir] = end_pos - dsign*2.0*fitting_hlen; // extend into pipe
				objs.emplace_back(fitting, TYPE_PIPE, room_id, !tank_dim, 0, hfflags, tot_light_amt, SHAPE_CYLIN, fitting_color); // horizontal
				set_wall_width(fitting, end_pos, (merged_pipe_radius + fitting_exp), !tank_dim);
				set_cube_zvals(fitting, hpipe_zc, hpipe_zc);
				if (conn_up) {fitting.z2() += merged_pipe_radius + fitting_hlen;} // up
				else         {fitting.z1() -= merged_pipe_radius + fitting_hlen;} // down
				vfflags ^= (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI); // swap end bits
				objs.emplace_back(fitting, TYPE_PIPE, room_id, 0, 1, vfflags, tot_light_amt, SHAPE_CYLIN, fitting_color); // vertical
				// add valve with fitting along v_pipe
				float const valve_radius(2.5*merged_pipe_radius), valve_depth(0.45*valve_radius), outside_edge(end_pos + dsign*merged_pipe_radius);
				float const valve_height(conn_up ? (hpipe_zc + 2.0*valve_radius) : (zval + 0.8*get_player_eye_height()));
				set_wall_width(fitting, valve_height, 1.5*fitting_hlen, 2); // set zval
				vfflags |= (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI); // set both top and bottom end bits
				objs.emplace_back(fitting, TYPE_PIPE, room_id, 0, 1, vfflags, tot_light_amt, SHAPE_CYLIN, fitting_color); // vertical
				cube_t valve;
				set_wall_width(valve, valve_height,     valve_radius, 2); // set zval
				set_wall_width(valve, conn_pipe_center, valve_radius, tank_dim);
				valve.d[!tank_dim][!conn_dir] = outside_edge; // outside edge of pipe
				valve.d[!tank_dim][ conn_dir] = outside_edge + dsign*valve_depth;
				objs.emplace_back(valve, TYPE_VALVE, room_id, !tank_dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, RED);
			}
			for (point const &pos : tank_conn_pts) {
				// horizontal fitting to connector pipe
				cube_t fitting(h_fitting);
				set_wall_width(fitting, pos[!tank_dim], fitting_hlen, !tank_dim); // Note: extends off the ends of the pipe
				objs.emplace_back(fitting, TYPE_PIPE, room_id, !tank_dim, 0, (fittings_flags | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, fitting_color);
				// vertical fittings for top and bottom of pipe
				fitting.set_from_point(pos);
				fitting.expand_by_xy(1.2*pipe_radius);
				fitting.z1() -= merged_pipe_radius + 1.5*pipe_radius;
				objs.emplace_back(fitting, TYPE_PIPE, room_id, 0, 1, (fittings_flags | RO_FLAG_ADJ_LO), tot_light_amt, SHAPE_CYLIN, fitting_color); // top
				set_cube_zvals(fitting, (tank_z2 - pipe_radius), (tank_z2 + 1.5*pipe_radius));
				objs.emplace_back(fitting, TYPE_PIPE, room_id, 0, 1, (fittings_flags | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, fitting_color); // bottom
			} // for pos
		}
	}
	// add machines along the walls
	unsigned const num_machines((rgen.rand() % 11) + 10); // 10-20
	avoid.d[edim][!edir] = room.d[edim][!edir]; // extend to opposite end of the factory so that we have space to place a ladder and sprinkler pipe

	for (unsigned n = 0; n < num_machines; ++n) {
		float const height(min(fc_gap*rgen.rand_uniform(0.3, 0.7), max_height));
		cube_t c;
		set_cube_zvals(c, zval, zval+height);

		for (unsigned i = 0; i < 25; ++i) { // make 25 attempts to place the object
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
			float const hwidth(0.5*max_sz[!dim]*rgen.rand_uniform(0.75, 1.0));
			float const depth(min(2.0f*hwidth, max_sz[dim])*rgen.rand_uniform(0.75, 1.0));
			float center(rgen.rand_uniform(place_area.d[!dim][0]+hwidth, place_area.d[!dim][1]-hwidth)); // random position
			c.d[dim][ dir] = place_area.d[dim][dir];
			c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
			set_wall_width(c, center, hwidth, !dim);
			if (c.intersects(avoid) || overlaps_obj_or_placement_blocked(c, room, objs_start)) continue;
			objs.emplace_back(c, TYPE_MACHINE, room_id, dim, !dir, RO_FLAG_IN_FACTORY, tot_light_amt, SHAPE_CUBE, LT_GRAY, rgen.rand()); // add more randomness
			set_obj_id(objs);
			break; // done
		} // for i
	} // for n
}

void building_t::add_chem_tank_gauge(room_object_t const &tank, float radius, float height) {
	bool gdim(tank.dim), gdir(tank.dir);
	if (gdim == 1 && gdir == 1 && height < 0.5*tank.dz()) {gdir ^= 1;} // swap gauge if in the bottom half to avoid placing over label
	float const gauge_depth(0.25*radius), tank_side(tank.d[gdim][gdir]);
	cube_t gauge;
	set_wall_width(gauge, (tank.z1() + height), radius, 2); // Z
	set_wall_width(gauge, tank.get_center_dim(!gdim), radius, !gdim);
	gauge.d[gdim][!gdir] = tank_side;
	gauge.d[gdim][ gdir] = tank_side + (gdir ? 1.0 : -1.0)*gauge_depth;
	interior->room_geom->objs.emplace_back(gauge, TYPE_GAUGE, tank.room_id, gdim, gdir, RO_FLAG_NOCOLL, tank.light_amt, SHAPE_CYLIN, BRASS_C);
}

// Note: assumes building is a factory, but has no effect if it's not
void building_room_geom_t::set_factory_machine_seed(unsigned rseed1, unsigned rseed2) {
	for (room_object_t& obj : objs) {
		if (obj.type != TYPE_MACHINE || !obj.in_factory())    continue;
		if (obj.item_flags == rseed1 && obj.obj_id == rseed2) continue; // already set
		obj.item_flags = rseed1;
		obj.obj_id     = rseed2;
		invalidate_draw_data_for_obj(obj);
	} // for obj
}

