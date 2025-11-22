// 3D World - Building Retail Room Setup and Object Placement
// by Frank Gennari 2/08/2025

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t


extern object_model_loader_t building_obj_model_loader;

unsigned get_srack_num_shelves(room_object_t const &c);


bool check_for_overlap(cube_t const &c, vect_room_object_t const &objs, unsigned objs_start, float xy_spacing) {
	cube_t test_cube(c);
	test_cube.expand_by_xy(xy_spacing);
	return has_bcube_int(test_cube, objs, objs_start);
}
bool building_t::get_retail_long_dim() const {
	cube_t const &part(get_retail_part());
	return (part.dx() < part.dy());
}
bool building_t::add_retail_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, light_ix_assign_t &light_ix_assign) {
	// Note: this room should occupy the entire floor, so walkable room bounds == room == part
	assert(has_room_geom());
	for (unsigned d = 0; d < 2; ++d) {interior->room_geom->shelf_rack_occluders[d].clear();}
	float const floor_spacing(get_window_vspace()), dx(room.dx()), dy(room.dy()), spacing(0.7);
	float const door_width(get_doorway_width()), se_pad(0.8*door_width), nom_aisle_width(1.5*door_width), rack_height(SHELF_RACK_HEIGHT_FS*floor_spacing);
	unsigned const nx(max(1U, unsigned(spacing*dx/floor_spacing))), ny(max(1U, unsigned(spacing*dy/floor_spacing))); // same spacing as room lights
	bool const dim(dx < dy); // long dim (could also use get_retail_long_dim())
	float const length(dim ? dy : dx), width(dim ? dx : dy), max_rack_width(0.5*floor_spacing);
	if (width < 4.0*nom_aisle_width) return 0; // too small for shelf racks
	unsigned const nrows((dim ? nx : ny)-1), nracks(max(2U, (dim ? ny : nx)/4));
	if (nrows == 0) return 0; // too small for shelf racks
	float row_aisle_width(nom_aisle_width), aisle_spacing((width - row_aisle_width)/nrows), rack_width(aisle_spacing - row_aisle_width);
	assert(rack_width > 0.0);
	
	if (rack_width > max_rack_width) { // rack is too wide; widen the aisle instead
		row_aisle_width += 0.5*(rack_width - max_rack_width); // first aisle width increases by half the rack width; aisle_spacing is unchanged
		rack_width       = max_rack_width;
	}
	float const wall_thickness(get_wall_thickness()), pillar_width(2.0*wall_thickness), fc_gap(get_floor_ceil_gap()), clearance(get_min_front_clearance_inc_people());
	float const col_aisle_width(nom_aisle_width + pillar_width), rack_spacing((length - col_aisle_width)/nracks), rack_length(rack_spacing - col_aisle_width);
	float const e_width(0.9*door_width), pair_width(2.1*e_width), end_pad(1.2*door_width); // for escalators
	assert(rack_length > 0.0);
	vect_room_object_t &objs(interior->room_geom->objs);
	vect_cube_t pillars, rack_bcubes;
	cube_t rack, pillar, srack_avoid;
	set_cube_zvals(rack,   zval, (zval + rack_height));
	set_cube_zvals(pillar, zval, room.z2()-get_fc_thickness()); // up to the ceiling
	unsigned const objs_start(objs.size()), style_id(rgen.rand()); // same style for each rack
	unsigned rack_id(0);
	bool const skip_middle_row(nrows & 1), tall_retail(has_tall_retail());
	for (unsigned d = 0; d < 2; ++d) {interior->room_geom->shelf_rack_occluders[d].reserve((nrows - skip_middle_row)*nracks);}

	if (tall_retail) {
		srack_avoid = room;
		set_wall_width(srack_avoid, room.get_center_dim(!dim), (door_width + 0.5*pair_width), !dim);
	}
	//vect_room_object_t temp_objs;
	
	for (unsigned n = 0; n < nrows; ++n) { // n+1 aisles
		if (skip_middle_row && n == nrows/2) continue; // skip middle aisle rack to allow direct access to central stairs and elevator/escalator
		float const rack_lo(room.d[!dim][0] + row_aisle_width + n*aisle_spacing);
		rack.d[!dim][0] = rack_lo;
		rack.d[!dim][1] = rack_lo + rack_width;
		set_wall_width(pillar, (rack_lo + 0.5*rack_width), 0.5*pillar_width, !dim); // centered on the rack

		for (unsigned r = 0; r < nracks; ++r) {
			float const start(room.d[dim][0] + col_aisle_width + r*rack_spacing);
			rack.d[dim][0] = start;
			rack.d[dim][1] = start + rack_length;
			bool was_shortened(0);

			for (unsigned n = 0; n < 5; ++n) { // try to trim ends to make it fit
				cube_t cand(rack);
				if      (n == 1) {cand.d[dim][0] += 0.25*rack_length;} // 25% length reduction
				else if (n == 2) {cand.d[dim][1] -= 0.25*rack_length;}
				else if (n == 3) {cand.d[dim][0] += 0.50*rack_length;} // 50% length reduction
				else if (n == 4) {cand.d[dim][1] -= 0.50*rack_length;}
				// no other items have been placed in this room yet (other than lights), and there are no interior doors,
				// and all exterior doors have aisle_spacing around them, so we only need to check for stairs and elevators
				cube_t test_cube(cand);
				test_cube.expand_by_xy(se_pad); // add extra padding
				if (interior->is_blocked_by_stairs_or_elevator(test_cube)) {was_shortened = 1; continue;} // blocked
				rack_bcubes.push_back(cand); // can add to the upper floor even if blocked by an escalator because it can be clipped in length
				if (tall_retail && cand.intersects_xy(srack_avoid)) break; // may block potential escalator; shortening won't help
				add_shelf_rack(cand, dim, style_id, rack_id, room_id, 0, 0, 1, rgen); // extra_flags=0, item_category=0, add_occluers=1
				//for (unsigned n = 0; n < 2; ++n) {interior->room_geom->get_shelfrack_objects(objs.back(), temp_objs, n);}
				break; // done
			} // for n
			if (!was_shortened && r > 0) { // place a pillar at the end of the rack
				pillar.d[dim][0] = start - pillar_width;
				pillar.d[dim][1] = start;
				add_retail_pillar(pillar, zval, room_id, tall_retail);
				if (tall_retail) {pillars.push_back(pillar);} // needed for upper glass floors
			}
		} // for r
	} // for n
	if (tall_retail) { // maybe add a pair of escalators and a partial upper level glass floor
		bool const add_glass_floor = 1;
		float const floor_thickness(get_floor_thickness()), e_height(room.dz() - floor_thickness), delta_z(e_height - get_floor_ceil_gap());
		float const e_length(1.0*delta_z + 2.0*door_width); // upward at 45 degree angle + entrance/exit
		float const door_extra_pad(0.5*door_width + (add_glass_floor ? 2.0*floor_spacing : 0.0)); // add extra spaacing so that the glass floor at the top is wide enough
		cube_t centered;
		set_cube_zvals(centered, zval, (zval + e_height));
		set_wall_width(centered, room.get_center_dim(!dim), 0.5*pair_width, !dim);
		set_wall_width(centered, room.get_center_dim( dim), (0.5*e_length + end_pad), dim);
		// try to find a valid escalator placement by starting at the center and translating to each side
		float const gap(centered.d[dim][0] - room.d[dim][0] - door_extra_pad); // should be the same at each end; room has exterior walls, so no place area shrink

		if (gap > 0.0) { // we have the space to add an escalator; should always be true
			unsigned const num_steps = 10;
			float const step_len(gap/num_steps);
			bool const first_dir(rgen.rand_bool());
			bool success(0);

			for (unsigned step = 0; step <= num_steps && !success; ++step) { // take N steps to each side
				for (unsigned D = 0; D < (step ? 2U : 1U); ++D) { // first step is length 0 and has no dir
					bool const dir(bool(D) ^ first_dir);
					cube_t cand(centered);
					cand.translate_dim(dim, (dir ? 1.0 : -1.0)*step*step_len);
					// Note: we don't need to pad stairs and elevators because we've already padded the ends of the escalator, and double padding on each side isn't needed
					if (has_bcube_int(cand, interior->elevators)) continue; // no need to check for other escalators (should be none)
					bool has_int(0);

					for (stairwell_t const &s : interior->stairwells) { // check stairs coll
						if (s.z2() < room.z2() || !s.intersects(room)) continue; // wrong stairs
						cube_t stairs_ext(s);
						stairs_ext.d[s.dim][!s.dir] += (s.dir ? -1.0 : 1.0)*s.get_retail_landing_width(floor_spacing);
						if (stairs_ext.intersects(cand)) {has_int = 1; break;}
					}
					if (has_int) continue;
					cube_t upper_conn; // union of high ends of both escalators

					for (unsigned s = 0; s < 2; ++s) { // place two side-by-side escalators with opposite directions
						// extend 90% of floor thickness below; enough to hide building people animated feet, but not enough to clip through the ceiling below
						escalator_t e(cand, dim, dir, (bool(s) ^ dim ^ dir), door_width, delta_z, 0.9*floor_thickness, room_id, 0); // in_mall=0
						e.expand_in_dim(dim, -end_pad);
						e.d[!dim][!s] = cand.d[!dim][s] + (s ? -1.0 : 1.0)*e_width;
						interior->escalators.push_back(e);
						cube_t lo_end, hi_end;
						e.get_ends_bcube(lo_end, hi_end, 0); // exclude_sides=0
						upper_conn.assign_or_union_with_cube(hi_end);
					}
					success = 1;

					// add an upper floor connected to the escalator; it's too late to add to interior->floors, since the VBO has already been created
					if (!add_glass_floor) break; // done
					float const floor_z1(upper_conn.z1() - 0.5*floor_thickness), floor_z2(upper_conn.z1());
					cube_t upper_floor(room);
					set_cube_zvals(upper_floor, floor_z1, floor_z2);
					upper_floor.d[dim][!dir] = upper_conn.d[dim][!dir] + (dir ? 1.0 : -1.0)*0.1*upper_conn.get_sz_dim(dim); // connect to upper end of escalator
					interior->room_geom->glass_floors.push_back(upper_floor);
					// add upper floor railing to each side of the escalator
					colorRGBA const railing_color(railing_colors[rgen.rand()%3]);
					unsigned const railing_flags(RO_FLAG_TOS | RO_FLAG_OPEN | RO_FLAG_ADJ_TOP);
					float const railing_thickness(0.5*wall_thickness), railing_z1(floor_z2 + 0.25*railing_thickness);
					cube_t railing(upper_floor);
					set_cube_zvals(railing, railing_z1, railing_z1+fc_gap);
					railing.d[dim][dir] = upper_floor.d[dim][!dir] + (dir ? 1.0 : -1.0)*railing_thickness;
					railing.expand_in_dim(!dim, -get_trim_thickness()); // shrink slightly to prevent Z-fighting with building exterior
					vect_cube_t obj_parts;
					
					for (unsigned d = 0; d < 2; ++d) {
						cube_t r(railing);
						r.d[!dim][!d] = upper_conn.d[!dim][d]; // flush with the escalator side
						subtract_stairs_and_elevators_from_cube(r, obj_parts);
						for (cube_t const &c : obj_parts) {objs.emplace_back(c, TYPE_RAILING, room_id, !dim, d, railing_flags, 1.0, SHAPE_CUBE, railing_color);}
					}
					// add support beams connected to pillars
					float const beam_hwidth(0.36*pillar_width), beam_thickness(0.3*pillar_width), beam_frame_depth(0.25*beam_hwidth);
					float const beam_z2(upper_floor.z1() + 0.1*beam_thickness), beam_z1(beam_z2 - beam_thickness); // slightly overlapping with bottom of upper floor
					float const metal_frame_z2(upper_floor.z2() + 0.1*beam_thickness);
					vector<float> pillar_xy[2];

					for (cube_t const &pillar : pillars) {
						if (!pillar.intersects_xy(upper_floor)) continue;
						pillar_xy[0].push_back(pillar.xc());
						pillar_xy[1].push_back(pillar.yc());
					}
					for (unsigned d = 0; d < 2; ++d) {
						sort_and_unique(pillar_xy[d]);
						cube_t beam_area(upper_floor);
						set_cube_zvals(beam_area, beam_z1, beam_z2); // below the floor
						unsigned const skip_faces(get_skip_mask_for_xy(!d)); // skip ends

						for (float v : pillar_xy[d]) {
							set_wall_width(beam_area, v, beam_hwidth, d);
							subtract_stairs_and_elevators_from_cube(beam_area, obj_parts);

							for (cube_t const &c : obj_parts) {
								objs.emplace_back(c, TYPE_METAL_BAR, room_id, !d, 0, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, BLACK, skip_faces); // item_flags=skip_faces
							}
						} // for v
					} // for d
					bool const add_wall_frame = 1;

					if (add_wall_frame || !pillar_xy[0].empty()) { // wall frame, or at least one pillar
						// add an additional metal bar under the railing to cover the exposed end
						cube_t beam(upper_floor);
						set_cube_zvals(beam, beam_z1, metal_frame_z2); // was previously upper_floor.zc()
						set_wall_width(beam, upper_floor.d[dim][!dir], 0.5*beam_thickness, dim);
						unsigned const item_flags(get_skip_mask_for_xy(!dim)); // skip ends
						objs.emplace_back(beam, TYPE_METAL_BAR, room_id, !dim, !dir, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, BLACK, item_flags);
					}
					if (add_wall_frame) { // add framing around the exterior wall
						for (unsigned fdim = 0; fdim < 2; ++fdim) {
							for (unsigned fdir = 0; fdir < 2; ++fdir) {
								if ((bool)fdim == dim && (bool)fdir == !dir) continue; // skip railing side
								cube_t beam(upper_floor);
								set_cube_zvals(beam, beam_z1, metal_frame_z2); // below the floor spanning to slightly above
								beam.d[fdim][!fdir] = room.d[fdim][fdir] + (fdir ? -1.0 : 1.0)*beam_frame_depth;
								unsigned const item_flags(get_skip_mask_for_xy(!fdim) | ~get_face_mask(fdim, fdir)); // skip ends and side facing wall
								objs.emplace_back(beam, TYPE_METAL_BAR, room_id, fdim, !fdir, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, BLACK, item_flags);
							} // for fdir
						} // for fdim
					}
					if (1) { // add frame around escalator
						cube_t beam(upper_conn);
						set_cube_zvals(beam, beam_z1, metal_frame_z2); // below the floor spanning to slightly above
						beam.expand_by_xy(beam_frame_depth);
						beam.d[dim][!dir] = upper_floor.d[dim][!dir]; // clip to upper floor edge
						unsigned const item_flags(~get_face_mask(dim, !dir)); // skip side facing railing
						objs.emplace_back(beam, TYPE_METAL_BAR, room_id, dim, dir, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, BLACK, item_flags);
					}
					// add another level of shelf racks on this floor, but no shopping carts;
					// stairs/elevators/escalators/pillars extend through both floors, so we don't need to recheck intersections, but we may need to increase their height
					unsigned const objs_end(objs.size());
					float const min_rack_len(min(0.5f*rack_length, 1.0f*floor_spacing)), min_pillar_z2(floor_z2 + 0.85*fc_gap);
					cube_t upper_place_area(upper_floor), avoid_area(upper_conn);
					upper_place_area.expand_by_xy(-nom_aisle_width); // add space around the edges for aisles
					avoid_area.expand_by_xy(1.25*door_width);

					for (cube_t const &rack : rack_bcubes) {
						if (!upper_place_area.intersects_xy(rack)) continue;
						cube_t cand(rack);
						cand.intersect_with_cube_xy(upper_place_area); // clip to fit in upper floor area
						if (cand.get_sz_dim(dim) < min_rack_len) continue; // too short
						cand.translate_dim(2, (floor_z2 - rack.z1())); // move to the floor above
						if (cand.intersects_xy(avoid_area)) {cand.d[dim][!dir] = avoid_area.d[dim][dir];} // too close to escalator, shorten
						add_shelf_rack(cand, dim, style_id, rack_id, room_id, RO_FLAG_ON_FLOOR, 0, 1, rgen); // flag so that bot surf is drawn; item_category=0, add_occluders=0
					} // for rack
					if (!pillars.empty()) { // raise office pillar outer cubes above upper level
						for (unsigned i = objs_start; i < objs_end; ++i) {
							room_object_t &obj(objs[i]);
							if (obj.type == TYPE_OFF_PILLAR && upper_place_area.intersects_xy(obj)) {max_eq(obj.z2(), min_pillar_z2);}
						}
					}
					// re-enable this floor on any elevator passing through it; add short beams under elevator entrances
					for (elevator_t &e : interior->elevators) {
						if (e.skip_floors_mask == 0 || e.z2() < floor_z2 || e.z1() > floor_z1 || !upper_floor.contains_cube_xy(e)) continue;
						// check for opening clearance
						float const front_clearance(ELEVATOR_STAND_DIST*floor_spacing + 0.5*clearance);
						cube_t const bc_pad(e.get_bcube_padded(front_clearance));
						if (!upper_floor.contains_cube_xy(bc_pad)) continue; // too close to edge of floor; can this fail?
						if (upper_conn.intersects(bc_pad))         continue; // too close to escalator
						unsigned const floor_ix(round_fp((floor_z1 - e.z1())/floor_spacing));
						e.skip_floors_mask &= ~(1ULL << floor_ix); // unset this floor

						// add beams around the edges of the elevator
						for (unsigned fdim = 0; fdim < 2; ++fdim) {
							for (unsigned fdir = 0; fdir < 2; ++fdir) {
								float const edge(e.d[fdim][fdir]);
								cube_t beam(e);
								set_cube_zvals(beam, beam_z1, metal_frame_z2); // below the floor spanning to slightly above
								beam.d[fdim][!fdir] = edge - (fdir ? 1.0 : -1.0)*0.4*wall_thickness; // interior - cover most of the door gap at the floor
								beam.d[fdim][ fdir] = edge + (fdir ? 1.0 : -1.0)*beam_frame_depth; // exterior
								beam.expand_in_dim(!fdim, beam_frame_depth); // cover the corners
								unsigned const item_flags(get_skip_mask_for_xy(!fdim) | ~get_face_mask(fdim, !fdir)); // skip ends and side facing elevator
								objs.emplace_back(beam, TYPE_METAL_BAR, room_id, fdim, fdir, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, BLACK, item_flags);
							} // for fdir
						} // for fdim
					} // for e
					break; // done
				} // for dir
			} // for step
		}
	}
	if (!doors.empty()) { // add checkout counters
		// find union of primary stairs and central elevators, which forms the other bounds of our checkout counter
		cube_t blocked(room.get_cube_center()); // block off the center

		for (auto const &e : interior->elevators) {
			if (e.intersects(room)) {blocked.union_with_cube(e);}
		}
		for (auto const &e : interior->escalators) {
			if (e.intersects(room)) {blocked.union_with_cube(e);}
		}
		for (auto const &s : interior->stairwells) {
			if (s.intersects(room)) {blocked.union_with_cube(s);}
		}
		for (unsigned dir = 0; dir < 2; ++dir) { // each end of room/each door
			cube_t checkout(room);
			checkout.d[dim][!dir] = blocked.d[dim][dir];
			checkout.expand_in_dim(dim, -1.6*nom_aisle_width); // shrink
			unsigned const check_objs_start(skip_middle_row ? objs.size() : objs_start); // can skip objs check if middle row is skipped
			add_checkout_objs(checkout, zval, room_id, 1.0, check_objs_start, dim, dir, rgen.rand_bool()); // tot_light_amt=1.0
		}
	}
	add_shopping_carts_to_room(rgen, room, zval, room_id, 1.0, objs_start, (nrows*nracks/4)); // light_amt=1.0
	add_cameras_to_room(rgen, room, zval, room_id, 1.0, objs.size()); // tot_light_amt=1.0
	//cout << TXT(temp_objs.size()) << endl;

	if (has_tall_retail()) { // add a small wall light at the top of the stairs since the stairwell is extra tall
		// find central stairs (should be first stairs), but might not be the vertical segment in the retail area
		for (stairwell_t const &s : interior->stairwells) {
			if (s.z2() < room.z2() || !s.intersects(room)) continue; // wrong stairs
			add_U_stair_landing_lights(s, room_id, light_ix_assign.get_next_ix(), (room.z2() - floor_spacing));
			break;
		}
	}
	return 1;
}

void building_t::add_shopping_carts_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float light_amt, unsigned objs_start, unsigned max_carts) {
	if (max_carts == 0 || !building_obj_model_loader.is_model_valid(OBJ_MODEL_SHOP_CART)) return;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SHOP_CART)); // D, W, H
	float const height(0.425*get_window_vspace()), hlen(0.5*height*sz.x/sz.z), hwidth(0.5*height*sz.y/sz.z);
	float const clearance(get_min_front_clearance_inc_people());
	unsigned const num_carts(1 + (rgen.rand() % max_carts));
	cube_t cart_area(room);
	cart_area.expand_by_xy(-(max(hlen, hwidth) + get_wall_thickness()));
	cube_t cart;
	set_cube_zvals(cart, zval, zval+height);

	if (cart_area.is_strictly_normalized()) { // should always be true
		for (unsigned n = 0; n < num_carts; ++n) {
			bool const dim(rgen.rand_bool());
			vector2d const sz((dim ? hwidth : hlen), (dim ? hlen : hwidth));

			for (unsigned N = 0; N < 10; ++N) { // 10 random placement tries
				for (unsigned d = 0; d < 2; ++d) {set_wall_width(cart, rgen.rand_uniform(cart_area.d[d][0], cart_area.d[d][1]), sz[d], d);}
				cube_t c_exp(cart);
				c_exp.expand_in_dim(dim, clearance); // make sure there's clearance in front and behind to push; should align carts to aisles
				if (overlaps_other_room_obj(c_exp, objs_start) || interior->is_blocked_by_stairs_or_elevator(c_exp)) continue;
				interior->room_geom->objs.emplace_back(cart, TYPE_SHOP_CART, room_id, dim, rgen.rand_bool(), 0, 1.0);
				break; // success/done
			}
		} // for n
	}
}

bool building_t::add_small_retail_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float light_amt) { // for prisons, etc.
	// Note: simplified version of mall retail stores from building_t::add_mall_store_objs()
	bool const dim(room.dx() < room.dy()); // long dim
	float const window_vspace(get_window_vspace()), door_width(get_doorway_width());
	cube_t const place_area(get_walkable_room_bounds(room));
	float const dx(place_area.dx()), dy(place_area.dy()), spacing(0.8), nom_aisle_width(1.2*door_width);
	unsigned const nx(max(1U, unsigned(spacing*dx/window_vspace))), ny(max(1U, unsigned(spacing*dy/window_vspace)));
	float const length(dim ? dy : dx), width(dim ? dx : dy), max_rack_width(0.45*window_vspace);
	unsigned const nrows((dim ? nx : ny)-1), nracks(max(2U, (dim ? ny : nx)/4));
	if (width < 4.0*nom_aisle_width || nrows < 2) return 0; // can't fit at least two rows
	float row_aisle_width(nom_aisle_width), aisle_spacing((width - row_aisle_width)/nrows), rack_width(aisle_spacing - row_aisle_width);
	assert(rack_width > 0.0);

	if (rack_width > max_rack_width) { // rack is too wide; widen the aisle instead
		row_aisle_width += 0.5*(rack_width - max_rack_width); // first aisle width increases by half the rack width; aisle_spacing is unchanged
		rack_width       = max_rack_width;
	}
	float const rack_spacing((length - nom_aisle_width)/nracks), rack_length(rack_spacing - nom_aisle_width);
	assert(rack_length > 0.0);
	cube_t rack;
	set_cube_zvals(rack, zval, (zval + 0.65*window_vspace));
	unsigned const style_id(rgen.rand() & 1); // same style for each rack: 3-4 shelves, no top or sides
	unsigned rack_id(0);

	for (unsigned n = 0; n < nrows; ++n) { // n+1 aisles
		float const rack_lo(place_area.d[!dim][0] + row_aisle_width + n*aisle_spacing);
		rack.d[!dim][0] = rack_lo;
		rack.d[!dim][1] = rack_lo + rack_width;

		for (unsigned r = 0; r < nracks; ++r) {
			float const start(place_area.d[dim][0] + nom_aisle_width + r*rack_spacing);
			rack.d[dim][0] = start;
			rack.d[dim][1] = start + rack_length;
			cube_t test_cube(rack);
			test_cube.expand_by_xy(0.6*door_width); // add extra padding
			if (is_obj_placement_blocked(test_cube, room, 1, 1)) continue; // inc doors and check open_dir
			add_shelf_rack(rack, dim, style_id, rack_id, room_id, 0, RETAIL_FOOD+1, 0, rgen); // add_occluders=0
		} // for r
	} // for n
	if (rack_id == 0) return 0; // no racks were added
	// add cash register/checkout counter?
	add_door_sign("Store", room, zval, room_id);
	return 1;
}

void building_t::add_retail_pillar(cube_t const &pillar, float zval, unsigned room_id, bool is_tall) {
	assert(pillar.is_strictly_normalized());
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(pillar, TYPE_OFF_PILLAR, room_id, 0, 0, 0);
	if (!is_tall) return;
	// add bare top part of pillars
	objs.back().flags |= RO_FLAG_ADJ_TOP; // must draw the top surface
	cube_t pillar_top(pillar);
	objs.back().z2() = pillar_top.z1() = zval + get_floor_ceil_gap(); // prev was bottom, next is top
	pillar_top.expand_by_xy(-0.1*get_wall_thickness()); // slight shrink
	assert(pillar_top.is_strictly_normalized());
	objs.emplace_back(pillar_top, TYPE_OFF_PILLAR, room_id, 0, 0, RO_FLAG_ADJ_HI); // flag as upper/concrete
}

void building_t::add_U_stair_landing_lights(stairwell_t const &s, unsigned room_id, unsigned light_ix, float floor_zval) {
	float const floor_spacing(get_window_vspace());
	float const landing_width(s.get_retail_landing_width(floor_spacing));
	float const radius(0.15*landing_width), dsign(s.dir ? -1.0 : 1.0), wall_pos(s.d[s.dim][!s.dir] + dsign*landing_width);
	cube_t light;
	light.d[s.dim][!s.dir] = wall_pos;
	light.d[s.dim][ s.dir] = wall_pos - dsign*0.01*floor_spacing; // extend out from wall
	set_wall_width(light, s.get_center_dim(!s.dim), radius, !s.dim); // width
	set_wall_width(light, (floor_zval - get_fc_thickness() + 0.75*floor_spacing), radius, 2); // height
	colorRGBA const light_color(get_light_color_temp(0.6)); // slightly blue-ish white
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(light, TYPE_LIGHT, room_id, s.dim, s.dir, (RO_FLAG_NOCOLL | RO_FLAG_LIT | RO_FLAG_ADJ_HI | RO_FLAG_RSTAIRS), 0.0, SHAPE_CYLIN, light_color);
	objs.back().obj_id = light_ix;
}

void building_t::add_checkout_objs(cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool dim, bool dir, bool cr_dir) {
	float const floor_spacing(get_window_vspace());
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t checkout(place_area);
	float const centerline(place_area.get_center_dim(!dim)), cr_space(1.5*get_doorway_width());
	float checkout_len(checkout.get_sz_dim(dim));

	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_CASHREG)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_CASHREG)); // D, W, H
		float const height(0.5*floor_spacing), hlen(0.5*height*sz.x/sz.z), hwidth(0.5*height*sz.y/sz.z);
		if (2.0*hlen > checkout_len) return; // not enough space to fit
		cube_t cashreg;
		vect_cube_t to_add;
		set_cube_zvals(cashreg, zval, (zval + height));
		set_wall_width(cashreg, centerline, hwidth, !dim);
		unsigned const num(max(1U, min(4U, unsigned(0.25*checkout_len/hlen))));

		for (unsigned n = 0; n < num; ++n) {
			set_wall_width(cashreg, (checkout.d[dim][0] + checkout_len*(n+1)/(num+1)), hlen, dim);
			// check shelfracks; shouldn't need to check stairs or elevators
			if (!check_for_overlap(cashreg, objs, objs_start, cr_space)) {to_add.push_back(cashreg);} // don't add directly as it may collide with the next cashreg
		}
		for (cube_t const &c : to_add) {objs.emplace_back(c, TYPE_CASHREG, room_id, dim, cr_dir, 0);}
	}
	else if (checkout_len > 0.4*floor_spacing) { // add if large enough
		checkout.expand_in_dim(dim, -0.1*(checkout_len - 0.4*floor_spacing)); // shrink slightly more if long
		set_wall_width(checkout, centerline, 0.2f*floor_spacing, !dim); // set width
		set_cube_zvals(checkout, zval, (zval + 0.4*floor_spacing));
		checkout_len = checkout.get_sz_dim(dim); // update after shrink
		unsigned const num_segs(max(1, round_fp(checkout_len/(1.33*floor_spacing))));
		float const seg_len(checkout_len/num_segs), seg_half_gap((num_segs == 1) ? 0.0 : 0.15*seg_len);

		for (unsigned n = 0; n < num_segs; ++n) {
			float const lo_end(checkout.d[dim][0] + n*seg_len);
			cube_t seg(checkout);
			seg.d[dim][0] = lo_end + seg_half_gap;
			seg.d[dim][1] = lo_end + seg_len - seg_half_gap;
			// check shelfracks; shouldn't need to check stairs or elevators
			if (!check_for_overlap(seg, objs, objs_start, cr_space)) {objs.emplace_back(seg, TYPE_CHECKOUT, room_id, dim, dir, 0);}
		}
	}
}

cube_t building_t::add_shelf_rack(cube_t const &c, bool dim, unsigned style_id, unsigned &rack_id, unsigned room_id,
	unsigned extra_flags, unsigned item_category, bool add_occluders, rand_gen_t &rgen)
{
	bool const is_empty(rgen.rand_float() < 0.05); // 5% empty
	unsigned flags(extra_flags | (is_empty ? 0 : RO_FLAG_NONEMPTY));
	unsigned const max_num_shelves((item_category == RETAIL_ELECTRONICS+1) ? 3 : 5); // Note: mall item_category has +1 added
	if (is_school() || is_prison()) {flags |= RO_FLAG_ADJ_HI;} // flag as no_alcohol
	room_object_t srack(c, TYPE_SHELFRACK, room_id, !dim, 0, flags, 1.0, SHAPE_CUBE, WHITE); // tot_light_amt=1.0
	srack.obj_id       = style_id;  // common for all racks
	srack.item_flags   = rack_id++; // unique per rack
	srack.drawer_flags = item_category;
	while (get_srack_num_shelves(srack) > max_num_shelves) {++srack.obj_id;}
	interior->room_geom->objs.push_back(srack);
	cube_t back, top, sides[2], shelves[5];
	unsigned const num_shelves(get_shelf_rack_cubes(srack, back, top, sides, shelves));
	interior->room_geom->get_shelfrack_objects(srack, interior->room_geom->objs, 1); // add_models_mode=1
	
	if (add_occluders) {
		interior->room_geom->shelf_rack_occluders[0].push_back(back);
		interior->room_geom->shelf_rack_occluders[1].push_back(top.is_all_zeros() ? shelves[num_shelves-1] : top); // top, or top shelf
	}
	return back; // return back for mall stores so that they can be used as occluders
}

point building_t::get_retail_upper_stairs_landing_center() const {
	assert(has_tall_retail());
	float const floor_spacing(get_window_vspace());
	room_t const &room(get_retail_room());

	for (stairwell_t const &s : interior->stairwells) {
		if (s.z2() < room.z2() || !s.intersects(room)) continue; // wrong stairs
		point pos;
		pos[!s.dim] = s.get_center_dim(!s.dim);
		pos[ s.dim] = s.d[s.dim][!s.dir] + (s.dir ? -1.0 : 1.0)*0.5*s.get_retail_landing_width(floor_spacing); // halfway out
		pos.z = room.z2() - floor_spacing;
		return pos;
	} // for s
	assert(0); // should never get here
	return all_zeros;
}

bool building_t::point_over_glass_floor(point const &pos, bool inc_escalator) const {
	if (!has_glass_floor()) return 0;
	float const floor_spacing(get_window_vspace());
	bool is_above_floor(0);

	for (cube_t const &f : interior->room_geom->glass_floors) {
		if (pos.z < f.z2() || pos.z > f.z2() + floor_spacing) continue; // below, or more than one floor above
		if (f.contains_pt_xy(pos)) return 1;
		is_above_floor = 1;
	}
	return (inc_escalator && is_above_floor && point_in_escalator(pos));
}

colorRGBA building_t::get_retail_light_color() const {
	assert(has_room_geom());
	vect_room_object_t const &objs(interior->room_geom->objs);
	assert(interior->room_geom->retail_start < objs.size());

	// find the first light placed in the retail area, starting with the first retail room object (which should be a light)
	for (auto i = objs.begin()+interior->room_geom->retail_start; i != objs.end(); ++i) {
		if (i->type == TYPE_LIGHT) {return (i->is_light_on() ? i->color : BLACK);}
	}
	assert(0);
	return BLACK; // not found?
}

cube_t escalator_t::get_ramp_bcube(bool exclude_sides) const {
	cube_t ramp(*this);
	if (exclude_sides) {ramp.expand_in_dim(!dim, -get_side_width());} // subtract off sides to get the belt/stairs
	cube_t lo_end(ramp), hi_end(ramp);
	ramp.expand_in_dim(dim, -end_ext); // subtract off ends
	ramp.z2() = z1() + delta_z;
	assert(ramp.is_strictly_normalized());
	return ramp;
}
void escalator_t::get_ends_bcube(cube_t &lo_end, cube_t &hi_end, bool exclude_sides) const {
	cube_t const ramp(get_ramp_bcube(exclude_sides));
	lo_end = hi_end = *this;
	for (unsigned d = 0; d < 2; ++d) {lo_end.d[!dim][d] = hi_end.d[!dim][d] = ramp.d[!dim][d];}
	float const side_height(get_side_height());
	lo_end.z2() = z1() + side_height;
	hi_end.z1() = ramp.z2();
	hi_end.z2() = hi_end.z1() + side_height;
	lo_end.d[dim][ dir] = ramp.d[dim][!dir];
	hi_end.d[dim][!dir] = ramp.d[dim][ dir];
}
cube_t escalator_t::get_side_for_end(cube_t const &end, bool lr) const {
	cube_t end_side(end);
	end_side.d[!dim][!lr] = end.d[!dim][lr] + (lr ? -1.0 : 1.0)*get_side_width();
	return end_side;
}
cube_t escalator_t::get_support_pillar() const {
	assert(!in_mall);
	float const support_radius(0.15*get_width());
	cube_t support(*this);
	support.z2() = z1() + delta_z - get_upper_hang();
	set_wall_width(support, get_center_dim(!dim), support_radius, !dim);
	set_wall_width(support, (d[dim][dir] + (dir ? -1.0 : 1.0)*0.5*end_ext), support_radius, dim); // center of upper end
	return support;
}
void escalator_t::get_ramp_bottom_pts(cube_t const &ramp, point bot_pts[4]) const { // {lo-left, lo-right, hi-right, hi-left} (or reversed)
	bot_pts[0].z = bot_pts[1].z = ramp.z1() - bot_edge_shift; // shift below the steps for extra space and to hide the walking feet of people
	bot_pts[2].z = bot_pts[3].z = ramp.z2() - bot_edge_shift;
	bot_pts[0][ dim] = bot_pts[1][ dim] = ramp.d[dim][!dir];
	bot_pts[2][ dim] = bot_pts[3][ dim] = ramp.d[dim][ dir];
	bot_pts[0][!dim] = bot_pts[3][!dim] = ramp.d[!dim][0];
	bot_pts[1][!dim] = bot_pts[2][!dim] = ramp.d[!dim][1];
	if (dim ^ dir ^ 1) {std::reverse(bot_pts, bot_pts+4);} // use correct vertex winding order
}
unsigned escalator_t::get_all_cubes(cube_t cubes[7]) const { // {lo left wall, lo right wall, lo floor, hi left wall, hi right wall, hi floor, [pillar]}
	cube_t ends[2];
	get_ends_bcube(ends[0], ends[1], 0); // exclude_sides=0
	unsigned ix(0);

	for (unsigned hi = 0; hi < 2; ++hi) {
		cube_t end(ends[hi]);
		for (unsigned lr = 0; lr < 2; ++lr) {cubes[ix++] = get_side_for_end(end, lr);}
		end.z2() = end.z1() + get_floor_thick();
		cubes[ix++] = end;
	}
	if (in_mall) return 6; // no pillar
	cubes[ix++] = get_support_pillar();
	return 7;
}

