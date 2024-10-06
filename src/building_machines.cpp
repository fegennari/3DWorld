// 3D World - Building Basement Machine Creation and Drawing
// by Frank Gennari 10/15/2024

#include "function_registry.h"
#include "buildings.h"


int get_cylin_duct_tid();
int get_cube_duct_tid ();
int get_ac_unit_tid   (unsigned ix);
tid_nm_pair_t get_metal_plate_tex(float tscale, bool shadowed);
colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c);


colorRGBA choose_pipe_color(rand_gen_t &rgen) {
	unsigned const NCOLORS = 6;
	colorRGBA const colors[NCOLORS] = {BRASS_C, COPPER_C, BRONZE_C, DK_BROWN, BROWN, LT_GRAY};
	return colors[rgen.rand() % NCOLORS];
}
colorRGBA choose_machine_part_color(rand_gen_t &rgen, bool is_small) {
	if (!is_small || rgen.rand_float() < 0.6) { // shade of gray
		float const lum(rgen.rand_uniform(0.1, 0.6));
		return colorRGBA(lum, lum, lum);
	}
	return choose_pipe_color(rgen);
}

void select_pipe_location(point &p1, point &p2, cube_t const &region, float radius, bool dim, rand_gen_t &rgen) {
	float const edge_space(1.5*radius);
	p1[dim] = region.d[dim][0];
	p2[dim] = region.d[dim][1];

	for (unsigned n = 0; n < 2; ++n) {
		unsigned const d(n ? 2 : !dim);
		p1[d] = p2[d] = rgen.rand_uniform(region.d[d][0]+edge_space, region.d[d][1]-edge_space);
	}
}

void building_room_geom_t::add_machine_pipe_in_region(room_object_t const &c, cube_t const &region, float rmax, bool dim, rand_gen_t &rgen) {
	assert(region.is_strictly_normalized());
	min_eq(rmax, 0.25f*min(region.dz(), region.get_sz_dim(!dim)));
	float const radius(rgen.rand_uniform(0.25, 1.0)*rmax);
	point p1, p2;
	select_pipe_location(p1, p2, region, radius, dim, rgen);
	colorRGBA const color(choose_pipe_color(rgen)), spec_color(get_specular_color(color)); // special case metals
	get_metal_material(1, 0, 1, 0, spec_color).add_cylin_to_verts(p1, p2, radius, radius, apply_light_color(c, color), 0, 0, 0, 0, 1.0, 1.0, 0, 16); // shadowed, small
}

void clip_cylin_to_square(cube_t &c, unsigned dim) { // about the center
	unsigned const d1((dim+1)%3), d2((dim+2)%3); // the non-cylin long dims
	float const da(c.get_sz_dim(d1)), db(c.get_sz_dim(d2));
	if (da < db) {c.expand_in_dim(d2, -0.5*(db - da));} else {c.expand_in_dim(d1, -0.5*(da - db));} // shrink
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

void building_room_geom_t::add_machine(room_object_t const &c) { // components are shadowed and small
	// can use AC Unit, metal plate texture, buttons, lights, etc.
	rand_gen_t rgen(c.create_rgen());
	float const height(c.dz()), width(c.get_width()), depth(c.get_depth());
	float const pipe_rmax(0.033*min(height, min(width, depth)));
	bool const dim(c.dim), dir(c.dir), two_part(width > rgen.rand_uniform(1.5, 2.2)*depth);
	unsigned const num_parts(two_part ? 2U : 1U);
	colorRGBA const base_color(apply_light_color(c, c.color));
	cube_t base(c), main(c);
	base.z2() = main.z1() = c.z1() + rgen.rand_uniform(0.04, 0.1)*height;
	cube_t parts[2] = {main, main};
	get_metal_material(1, 0, 1).add_cube_to_verts_untextured(base, base_color, EF_Z1); // skip bottom
	bool parts_swapped(0), is_cylins[2] = {0, 0};
	static vect_cube_t avoid, shapes;

	if (two_part) {
		float const split_pos(main.d[!dim][0] + rgen.rand_uniform(0.4, 0.6)*width);
		parts[0].d[!dim][1] = split_pos - max(0.0f, rgen.rand_uniform(-0.1, 0.1))*width; // maybe add a gap
		parts[1].d[!dim][0] = split_pos + max(0.0f, rgen.rand_uniform(-0.1, 0.1))*width; // maybe add a gap
		if (rgen.rand_bool()) {swap(parts[0], parts[1]); parts_swapped = 1;} // remove any bias toward the left/right
	}
	// calculate size and shape of each part
	for (unsigned n = 0; n < num_parts; ++n) {
		cube_t &part(parts[n]);
		bool const is_cylin(!is_cylins[0] && rgen.rand_float() < 0.4 && max(part.dx(), part.dy()) < 1.5*min(part.dx(), part.dy())); // don't place two cylinders
		is_cylins[n] = is_cylin;
		float const max_shrink_val(is_cylin ? 0.15 : 0.25);
		for (unsigned d = 0; d < 2; ++d) {part.expand_in_dim(d, -rgen.rand_uniform(0.1, 1.0)*max_shrink_val*part.get_sz_dim(d));} // shrink in both dims independently
		part.z2() -= rgen.rand_uniform(0.05, 0.2)*height; // make a bit shorter

		if (is_cylin) { // make it square
			clip_cylin_to_square(part, 2); // Z
			if (two_part && part.intersects(parts[1-n])) {part.expand_by_xy(-0.1*min(part.dx(), part.dy()));} // if adjacent, shrink to force a gap
		}
	}
	// add/draw shapes/objects for each part
	for (unsigned n = 0; n < num_parts; ++n) {
		cube_t const &part(parts[n]);
		bool const is_cylin(is_cylins[n]);
		vector3d const part_sz(part.get_size());
		colorRGBA const color(apply_light_color(c, choose_machine_part_color(rgen, 0))); // is_small=0
		avoid .clear();
		shapes.clear();

		if (is_cylin) {
			if (rgen.rand_bool()) { // make it textured
				tid_nm_pair_t const mplate_tex(get_metal_plate_tex(1.0, 1)); // shadowed
				rgeom_mat_t &mp_mat(get_material(mplate_tex, 1, 0, 1)); // shadowed, small
				mp_mat.add_vcylin_to_verts(part, color, 0, 1, 0, 0, 1.0, 1.0, 2.0, 1.0, 0, 32, 0.0, 0, 2.0); // draw sides and top; tscale=2.0
			}
			else {get_metal_material(1, 0, 1).add_vcylin_to_verts(part, color, 0, 1);} // draw sides and top
		}
		else {get_metal_material(1, 0, 1).add_cube_to_verts_untextured(part, color, EF_Z1);} // skip bottom

		if (!is_cylin && rgen.rand_float() < 0.6) { // maybe add a breaker panel if a cube
			//unsigned const flags(rgen.rand_bool() ? RO_FLAG_OPEN : 0); // open 50% of the time
			unsigned const flags(0); // not open, because the breakers aren't currently drawn
			float panel_hheight(part_sz.z*rgen.rand_uniform(0.15, 0.22)), panel_hwidth(part_sz[!dim]*rgen.rand_uniform(0.15, 0.22));
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
			}
		}
		if (!is_cylin && rgen.rand_float() < 0.75) { // maybe add a valve handle to the front if a cube
			float const valve_radius(min(5.0f*pipe_rmax, 0.4f*min(part_sz[!dim], part_sz.z))*rgen.rand_uniform(0.75, 1.0));
			float const valve_depth(valve_radius*rgen.rand_uniform(0.4, 0.5));
			cube_t valve(place_obj_on_cube_side(part, dim, dir, valve_radius, valve_radius, valve_depth, 1.1, rgen)); // Note: may extend a bit outside c in the front
			assert(valve.is_strictly_normalized());

			if (!has_bcube_int(valve, avoid)) {
				unsigned const NUM_HANDLE_COLORS = 5;
				colorRGBA const handle_colors[NUM_HANDLE_COLORS] = {colorRGBA(0.5, 0.0, 0.0), colorRGBA(0.7, 0.0, 0.0), colorRGBA(0.5, 0.25, 0.0), LT_GRAY, BKGRAY};
				colorRGBA const handle_color(handle_colors[rgen.rand() % NUM_HANDLE_COLORS]), shaft_color(choose_pipe_color(rgen));
				rgeom_mat_t &mat(get_metal_material(1, 0, 1)); // shadowed, small, specular metal
				// use the same material for the handle and the shaft
				draw_metal_handle_wheel(valve, dim, apply_light_color(c, handle_color), apply_light_color(c, shaft_color), mat, mat);
				avoid.push_back(valve);
			}
		}
		// add vents or AC unit fans
		if (!is_cylin && rgen.rand_float() < 0.9) {
			unsigned orients[3]={};
			orients[0] = 2*dim + dir; // front side

			if (two_part) { // side not facing the other part
				orients[1] = 2*(!dim) + (bool(n)^parts_swapped);
			}
			else { // both sides
				for (unsigned d = 0; d < 2; ++d) {orients[d+1] = 2*(!dim) + d;}
			}
			for (unsigned m = 0; m < (two_part ? 2 : 3); ++m) {
				if (rgen.rand_bool()) continue;
				bool const vdim(orients[m] >> 1), vdir(orients[m] & 1);
				bool const add_ac_unit(rgen.rand_float() < ((vdim == dim) ? 0.75 : 0.25)), mirror_y(add_ac_unit);
				float const pwidth(part_sz[!vdim]), pheight(part_sz.z);

				for (unsigned N = 0; N < 10; ++N) { // make 10 attempts to place it
					float const vhwidth((add_ac_unit ? rgen.rand_uniform(0.36, 0.4) : rgen.rand_uniform(0.24, 0.32))*min(pwidth, 1.0f*pheight));
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
		float const pipe_rmax(0.05*part_sz.get_min_val());
		bool const has_gap(fabs(part.d[dim][dir] - c.d[dim][dir]) > pipe_rmax); // add pipes if gap is large enough
		cube_t region(part);
		region.d[dim][ dir] = (is_cylin ? part.get_center_dim(dim) : part.d[dim][!dir]); // center plane of cylinder or back of cube
		region.d[dim][!dir] = c.d[dim][!dir]; // wall behind the machine
		assert(region.is_strictly_normalized());

		if (is_cylin || has_gap) { // if there's a gap between the machine and the wall
			// add pipe(s) connecting to back wall in {dim, !dir} or ceiling; there's no check for intersecting pipes
			unsigned const num_pipes(rgen.rand() % 4); // 0-3
			for (unsigned n = 0; n < num_pipes; ++n) {add_machine_pipe_in_region(c, region, pipe_rmax, dim, rgen);}
		}
		if (!is_cylin && has_gap && rgen.rand_float() < 0.65) { // add a duct to the wall
			float const radius(rgen.rand_uniform(2.0, 4.0)*pipe_rmax); // wider than a pipe
			point p1, p2;
			select_pipe_location(p1, p2, region, radius, dim, rgen);
			rgeom_mat_t &duct_mat(get_material(tid_nm_pair_t(get_cylin_duct_tid(), 1.0, 1), 1, 0, 1)); // shadowed, small
			duct_mat.add_cylin_to_verts(p1, p2, radius, radius, apply_light_color(c, LT_GRAY), 0, 0, 0, 0, 1.0, 1.0, 0, 32, 1.0, 1); // ndiv=32, swap_txy=1
			cube_t vent_bc(p1, p2);
			vent_bc.expand_by(radius); // close enough
			avoid.push_back(vent_bc);
		}
		// add smaller shapes to this one
		unsigned const num_shapes(rgen.rand() % (is_cylin ? 5 : 9) + (two_part ? 0 : 1)); // 0-4 for cylin, 0-8 for cube, +1 for single part
		float const max_shape_sz(0.25*part_sz.get_min_val()), cylin_radius(0.25*(part_sz.x + part_sz.y));
		point const part_center(part.get_cube_center());

		for (unsigned N = 0; N < num_shapes; ++N) {
			bool const add_cylin(rgen.rand_float() < 0.33);

			for (unsigned m = 0; m < 10; ++m) { // 10 attempts to place this shape
				unsigned const face(rgen.rand() & 3); // {top, front, left, right}
				unsigned const face_dims[4] = {2, dim, !dim, !dim}, face_dirs[4] = {1, dir, 0, 1};
				unsigned const fdim(face_dims[face]), fdir(face_dirs[face]);
				unsigned cdim(2); // defaults to Z
				if (add_cylin && fdim < 2) {cdim = (rgen.rand() % 3);} // cylin on part side can be oriented in X, Y, or Z
				vector3d half_sz;
				cube_t shape;
				for (unsigned d = 0; d < 3; ++d) {half_sz[d] = rgen.rand_uniform(0.4, 1.0)*max_shape_sz;}

				if (add_cylin) { // make sure cylinder ends are square; could use clip_cylin_to_square(), but this always makes shapes smaller
					unsigned const d1((cdim+1)%3), d2((cdim+2)%3); // the non-cylin long dims
					half_sz[d1] = half_sz[d2] = 0.5*(half_sz[d1] + half_sz[d2]);
				}
				for (unsigned d = 0; d < 3; ++d) {
					float const val((d == fdim) ? part.d[fdim][fdir] : rgen.rand_uniform(part.d[d][0], part.d[d][1]));
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
					if (is_cylin && cdim == 2 && !dist_xy_less_than(part_center, shape.get_cube_center(), (cylin_radius + shape_radius))) continue; // vcylin vs. vcylin
				}
				if (is_cylin && !circle_rect_intersect(part_center, cylin_radius, shape, 2)) continue; // new bcube vs. part cylinder
				if (has_bcube_int(shape, avoid)) continue;
				bool bad_place(0);

				for (cube_t const &s : shapes) { // check previously placed shapes for Z-fighting at clip boundaries
					if (!s.intersects(shape)) continue;
					if (s.x1()==shape.x1() || s.x2()==shape.x2() || s.y1()==shape.y1() || s.y2()==shape.y2() || s.z1()==shape.z1() || s.z2()==shape.z2()) {bad_place = 1; break;}
				}
				if (bad_place) continue;
				// TODO: textures
				colorRGBA const pcolor(choose_machine_part_color(rgen, 1)), spec_color(get_specular_color(pcolor)), lcolor(apply_light_color(c, pcolor));
				rgeom_mat_t &mat(get_metal_material(1, 0, 1, 0, spec_color)); // shadowed, small, specular metal
				if (add_cylin) {mat.add_ortho_cylin_to_verts(shape, lcolor, cdim, 1, 1);} // draw top and bottom in case they're visible
				else           {mat.add_cube_to_verts_untextured(shape, lcolor, 0);} // draw all faces since we don't track which are visible
				shapes.push_back(shape);
				break; // success/done
			} // for m
		} // for n
	} // for n
	if (two_part) { // connect the two parts with pipe(s); there's no check for intersecting pipes
		unsigned const num_pipes(rgen.rand() % 4); // 0-3

		if (num_pipes > 0) {
			cube_t region(parts[0]);
			min_eq(region.z2(), parts[1].z2()); // shared Z range
			max_eq(region.d[dim][0], parts[1].d[dim][0]); // shared range
			min_eq(region.d[dim][1], parts[1].d[dim][1]);
			float side_pos[2] = {};
			for (unsigned d = 0; d < 2; ++d) {side_pos[d] = (is_cylins[d] ? parts[d].get_center_dim(!dim) : parts[d].d[!dim][bool(d)^parts_swapped]);}

			if (side_pos[0] != side_pos[1]) {
				region.d[!dim][0] = min(side_pos[0], side_pos[1]);
				region.d[!dim][1] = max(side_pos[0], side_pos[1]);
				assert(region.is_strictly_normalized());
				// add flanges if large and connecting to a cube?
				for (unsigned n = 0; n < num_pipes; ++n) {add_machine_pipe_in_region(c, region, pipe_rmax, !dim, rgen);}
			}
		}
	}
}

bool building_t::add_machines_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const floor_spacing(get_window_vspace()), fc_gap(get_floor_ceil_gap());
	float const min_clearance(get_min_front_clearance_inc_people()), min_gap(max(get_doorway_width(), min_clearance));
	cube_t const place_area(get_walkable_room_bounds(room)); // ignore trim?
	cube_t avoid;
	avoid.set_from_sphere(point(place_area.xc(), place_area.yc(), zval), min_clearance);
	unsigned const flags(0), num_machines((rgen.rand() % 4) + 1); // 1-4
	float const sz_cap_mult((num_machines > 0) ? 0.5 : 1.0); // need to guarantee gap if two machines are placed on opposite walls of a narrow room
	vector2d const place_sz(place_area.dx(), place_area.dy());
	vector2d avail_sz, max_sz;
	for (unsigned d = 0; d < 2; ++d) {avail_sz[d] = min(0.4f*place_sz[d], sz_cap_mult*(place_sz[d] - min_gap));}
	for (unsigned d = 0; d < 2; ++d) {max_sz  [d] = min(avail_sz[d], 2.0f*avail_sz[!d]);} // keep aspect ratio <= 2:1
	if (min(max_sz.x, max_sz.y) < 0.5*floor_spacing) return 0; // too small of a room to place a machine
	vect_room_object_t &objs(interior->room_geom->objs);
	bool any_placed(0);

	for (unsigned n = 0; n < num_machines; ++n) {
		float const height(fc_gap*rgen.rand_uniform(0.6, 0.9));
		// similar to place_obj_along_wall(), but with custom size logic that depends on dim
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
			if (c.intersects(avoid) || overlaps_other_room_obj(c, objs_start)) continue;
			if (interior->is_blocked_by_stairs_or_elevator(c) || is_cube_close_to_doorway(c, room, 0.0, 1)) continue;
			objs.emplace_back(c, TYPE_MACHINE, room_id, dim, !dir, flags, tot_light_amt, SHAPE_CUBE, LT_GRAY);
			objs.back().item_flags = rgen.rand(); // add more randomness
			set_obj_id(objs);
			any_placed = 1;
			break; // done
		} // for i
	} // for n
	return any_placed;
}
