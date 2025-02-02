// 3D World - City Object Placement
// by Frank Gennari
// 05/21/23

#include "city_objects.h"
#include "tree_3dw.h" // for tree_placer_t
//#include "profiler.h"

float pond_max_depth(0.0);

extern bool enable_model3d_custom_mipmaps, player_in_walkway, player_in_skyway;
extern int display_mode, animate2, player_in_basement;
extern unsigned max_unique_trees;
extern float fticks;
extern colorRGBA sun_color;
extern tree_placer_t tree_placer;
extern city_params_t city_params;
extern object_model_loader_t building_obj_model_loader;
extern plot_divider_type_t plot_divider_types[];
extern textured_mat_t pool_deck_mats[];

city_flag_t create_flag(bool dim, bool dir, point const &base_pt, float height, float length, int flag_id=-1);
void get_building_ext_basement_bcubes(cube_t const &city_bcube, vect_cube_t &bcubes);
void get_walkways_for_city(cube_t const &city_bcube, vect_bldg_walkway_t &walkway_cands);
void add_house_driveways_for_plot(cube_t const &plot, vect_cube_t &driveways);
bool connect_buildings_to_skyway(cube_t &m_bcube, bool m_dim, cube_t const &city_bcube, vector<skyway_conn_t> &ww_conns);
void get_city_building_walkways(cube_t const &city_bcube, vector<building_walkway_t *> &bwws);
float get_inner_sidewalk_width();
cube_t get_plot_coll_region(cube_t const &plot_bcube);
void play_hum_sound(point const &pos, float gain, float pitch);
bool enable_instanced_pine_trees();

bool are_birds_enabled() {return building_obj_model_loader.is_model_valid(OBJ_MODEL_BIRD_ANIM);}


bool city_obj_placer_t::gen_parking_lots_for_plot(cube_t const &full_plot, vector<car_t> &cars, unsigned city_id, unsigned plot_ix,
	vect_cube_t &bcubes, vect_cube_t &colliders, vect_cube_t const &plot_cuts, rand_gen_t &rgen, bool add_cars)
{
	vector3d const nom_car_size(city_params.get_nom_car_size()); // {length, width, height}
	float const space_width(PARK_SPACE_WIDTH *nom_car_size.y); // add 50% extra space between cars
	float const space_len  (PARK_SPACE_LENGTH*nom_car_size.x); // space for car + gap for cars to drive through
	float const pad_dist   (max(1.0f*nom_car_size.x, get_min_obj_spacing())); // one car length or min building spacing
	float const sidewalk_width(get_sidewalk_width());
	cube_t plot(full_plot); // plot shrunk by padding
	plot.expand_by_xy(-pad_dist);
	if (bcubes.empty()) return 0; // shouldn't happen, unless buildings are disabled; skip to avoid perf problems with an entire plot of parking lot
	unsigned const buildings_end(bcubes.size()), first_corner(rgen.rand()&3); // 0-3
	bool const car_dim(rgen.rand() & 1), car_dir(rgen.rand() & 1); // car_dim: 0=cars face in X; 1=cars face in Y
	float const xsz(car_dim ? space_width : space_len), ysz(car_dim ? space_len : space_width);
	bool has_parking(0);
	vector_add_to(plot_cuts, bcubes); // avoid plot cuts from skylights and underground elevators to blocking bcubes
	unsigned const bcubes_coll_end(bcubes.size());
	vector<hcap_with_dist_t> hcap_cands;
	car_t car;
	car.park();
	car.cur_city = city_id;
	car.cur_road = plot_ix; // store plot_ix in road field
	car.cur_road_type = TYPE_PLOT;

	for (unsigned c = 0; c < 4; ++c) { // generate 0-4 parking lots per plot, starting at the corners, in random order
		unsigned const cix((first_corner + c) & 3), xdir(cix & 1), ydir(cix >> 1), wdir(car_dim ? xdir : ydir), rdir(car_dim ? ydir : xdir);
		float const dx(xdir ? -xsz : xsz), dy(ydir ? -ysz : ysz), dw(car_dim ? dx : dy), dr(car_dim ? dy : dx); // delta-width and delta-row
		point const corner_pos(plot.d[0][xdir], plot.d[1][ydir], (plot.z1() + 0.1*ROAD_HEIGHT)); // shift up slightly to avoid z-fighting
		assert(dw != 0.0 && dr != 0.0);
		unsigned const parking_lot_ix(parking_lots.size());
		// start as min size at the corner
		parking_lot_t cand(cube_t(corner_pos, corner_pos), car_dim, car_dir, rdir, city_params.min_park_spaces, city_params.min_park_rows, parking_lot_ix);
		cand.d[!car_dim][!wdir] += cand.row_sz  *dw;
		cand.d[ car_dim][!rdir] += cand.num_rows*dr;
		if (!plot.contains_cube_xy(cand)) {continue;} // can't fit a min size parking lot in this plot, so skip it (shouldn't happen)
		if (has_bcube_int_xy(cand, bcubes, pad_dist)) continue; // intersects a building - skip (can't fit min size parking lot)
		cand.z2() += plot.dz(); // probably unnecessary
		parking_lot_t park(cand);

		// try to add more parking spaces in a row
		for (; plot.contains_cube_xy(cand); ++cand.row_sz, cand.d[!car_dim][!wdir] += dw) {
			if (has_bcube_int_xy(cand, bcubes, pad_dist)) break; // intersects a building - done
			park = cand; // success: increase parking lot to this size
		}
		cand = park;
		// try to add more rows of parking spaces
		// the space between rows is too small for cars to drive, but this is really limited by the tiling of the parking lot texture
		for (; plot.contains_cube_xy(cand); ++cand.num_rows, cand.d[car_dim][!rdir] += dr) {
			if (has_bcube_int_xy(cand, bcubes, pad_dist)) break; // intersects a building - done
			park = cand; // success: increase parking lot to this size
		}
		assert(park.row_sz >= city_params.min_park_spaces && park.num_rows >= city_params.min_park_rows);
		assert(park.dx() > 0.0 && park.dy() > 0.0);
		bool &cdir(park.dir);

		if (1) { // add driveways connecting to parking lot
			bool const dw_dir(plot.get_center_dim(!car_dim) < park.get_center_dim(!car_dim)); // connect to road on closer side
			float const dw_pad_dist(0.25*pad_dist);

			for (unsigned d = 0; d < 2; ++d) { // try both dirs
				float const back_of_lot(park.d[car_dim][!cdir]); // driveway connects to this side
				cube_t driveway(park);
				driveway.d[ car_dim][ cdir] = back_of_lot;
				driveway.d[ car_dim][!cdir] = back_of_lot + (cdir ? -1.0 : 1.0)*1.25*space_width;
				driveway.d[!car_dim][  dw_dir] = full_plot.d[!car_dim][dw_dir] + (dw_dir ? 1.0 : -1.0)*sidewalk_width; // extend to the road, with a small gap for the curb

				if (has_bcube_int_xy(driveway, bcubes, dw_pad_dist)) { // too close to an object such as a building (or its door); try other dir
					cdir ^= 1;
				}
				else { // add the driveway
					// try to extend the driveway in the other dim if there's space, so that cars have room to back out of parking spaces
					float const max_ext(space_len - 0.5*space_width), dw_end(driveway.d[!car_dim][!dw_dir]);
					cube_t dw_ext(driveway);
					dw_ext.d[!car_dim][dw_dir] = dw_end;
					float const ext_factor[5] = {1.0, 0.9, 0.8, 0.7, 0.6};

					for (unsigned n = 0; n < 5; ++n) { // 5 possible extend dists
						float const extend_to(dw_end + (dw_dir ? -1.0 : 1.0)*ext_factor[n]*max_ext);
						dw_ext.d[!car_dim][!dw_dir] = extend_to;
						
						if (full_plot.contains_cube_xy(dw_ext) && !has_bcube_int_xy(dw_ext, bcubes, dw_pad_dist)) {
							driveway.d[!car_dim][!dw_dir] = extend_to;
							break;
						}
					} // for n
					driveways.emplace_back(driveway, !car_dim, dw_dir, plot_ix, parking_lot_ix);
					bcubes.push_back(driveway); // add to list of blocker bcubes
					break;
				}
			} // for d
		}
		if (rgen.rand_float() < 0.4) { // add solar roofs over parking lots 40% of the time
			bool blocked(0);

			if (skyway.valid) { // check for skyway blocking solar panels
				cube_t blocked_area(skyway.bcube);
				blocked_area.expand_by_xy(0.5*city_params.road_width);
				blocked = blocked_area.intersects_xy(park);
			}
			if (!blocked) {
				cube_t roof_bc(park);
				roof_bc.z1()  = plot.z2();
				roof_bc.z2() += 5.0*nom_car_size.z;
				roof_bc.expand_by_xy(0.06*roof_bc.dz()); // legs are outside of the parking area
				bool const panel_dir(car_dim ? 1 : rgen.rand_bool()); // if north/south, face south (northern hemisphere); if east/west, choose a random dir
				parking_solar_t const ps(roof_bc, car_dim, panel_dir, park.row_sz, park.num_rows);
				p_solar_groups.add_obj(ps, p_solars);
				vector_add_to(ps.get_legs(), colliders); // add legs to colliders but not blockers
				cube_t blocker(roof_bc);
				blocker.z1() += 0.5*roof_bc.dz(); // top half
				blocker.expand_by_xy(0.5*nom_car_size.x); // add tree clearance
				bcubes.push_back(blocker); // required for trees
			}
		}
		//car.cur_seg = (unsigned short)parking_lot_ix; // store parking lot index in cur_seg; or rather don't, because sort invalidates the index
		parking_lots.push_back(park);
		bcubes.push_back(park); // add to list of blocker bcubes so that no later parking lots or other city objects overlap this one
		//parking_lots.back().expand_by_xy(0.5*pad_dist); // re-add half the padding for drawing (breaks texture coord alignment)
		unsigned const nspaces(park.row_sz*park.num_rows);
		num_spaces += nspaces;

		// create parking spaces, fill the parking lot with cars, and assign handicap spaces
		vector<unsigned char> &used_spaces(parking_lots.back().used_spaces);
		used_spaces.resize(nspaces, 0); // start empty
		car.dim = car_dim; car.dir = cdir;
		point pos(corner_pos.x, corner_pos.y, plot.z2());
		pos[car_dim] += 0.5*dr + (cdir ? 0.15 : -0.15)*fabs(dr); // offset for centerline, biased toward the front of the parking space
		float const car_density(add_cars ? rgen.rand_uniform(city_params.min_park_density, city_params.max_park_density) : 0.0);

		for (unsigned row = 0; row < park.num_rows; ++row) { // car lengths
			pos[!car_dim] = corner_pos[!car_dim] + 0.5*dw; // start at the low end; add half offset for centerline
			bool prev_was_bad(0);

			for (unsigned col = 0; col < park.row_sz; ++col) { // car widths
				point const center(pos.x, pos.y, park.z2());
				parking_space_t pspace(center, car_dim, cdir, parking_lot_ix, row, col);

				if (prev_was_bad) { // previous car did a bad parking job, leave this space empty
					prev_was_bad   = 0;
					pspace.blocked = 1;
				}
				else if (add_cars && rgen.rand_float() < car_density) { // only half the spaces are filled on average
					point cpos(pos);
					cpos[ car_dim] += 0.05*dr*rgen.rand_uniform(-1.0, 1.0); // randomness of front amount
					cpos[!car_dim] += 0.12*dw*rgen.rand_uniform(-1.0, 1.0); // randomness of side  amount

					if (col+1 != park.row_sz && (rgen.rand()&15) == 0) { // occasional bad parking job
						cpos[!car_dim] += dw*rgen.rand_uniform(0.3, 0.35)*(rgen.rand_bool() ? 1.0 : -1.0);
						prev_was_bad    = 1;
					}
					car.set_bcube(pos, nom_car_size);
					cars.push_back(car);
					if ((rgen.rand()&7) == 0) {cars.back().dir ^= 1;} // pack backwards 1/8 of the time
					used_spaces[row*park.row_sz + col] = 1;
					has_parking     = 1;
					pspace.occupied = 1;
					++filled_spaces;
				}
				hcap_cands.emplace_back(hcap_space_t(center, 0.25*space_width, car_dim, cdir, pspaces.size()), plot, bcubes, buildings_end);
				pspaces.push_back(pspace);
				pos[!car_dim] += dw;
			} // for col
			pos[car_dim] += dr;
		} // for row
		if (!add_cars) continue;
		// generate colliders for each group of used parking space columns
		cube_t cur_cube(park); // set zvals, etc.
		unsigned row_min(park.num_rows), row_max(0);
		bool inside(0);

		for (unsigned col = 0; col <= park.row_sz; ++col) { // one extra iteration
			// mark space as blocked if any spaces in the row are blocked;
			// avoids creating diagonally adjacent colliders that cause dead ends and confuse path finding;
			// doesn't affect parking spaces or cars
			bool blocked(0);

			if (col < park.row_sz) { // not a valid column
				for (unsigned row = 0; row < park.num_rows; ++row) {
					if (used_spaces[row*park.row_sz + col] == 0) continue; // no car
					blocked = 1;
					min_eq(row_min, row);
					max_eq(row_max, row);
				}
			}
			if (!inside && blocked) { // start a new segment
				cur_cube.d[!car_dim][0] = corner_pos[!car_dim] + col*dw;
				inside = 1;
			}
			else if (inside && !blocked) { // end the current segment
				cur_cube.d[!car_dim][1] = corner_pos[!car_dim] + col*dw;
				cur_cube.d[ car_dim][0] = corner_pos[ car_dim] +  row_min   *dr; // set row span for range of cars
				cur_cube.d[ car_dim][1] = corner_pos[ car_dim] + (row_max+1)*dr;
				cur_cube.normalize();
				cur_cube.d[ car_dim][!cdir] -= (cdir ? -1.0 : 1.0)*0.25*fabs(dr); // remove 25% of back of parking space (by DW) since cars don't usually block this
				colliders.push_back(cur_cube);
				row_min = park.num_rows; row_max = 0; inside = 0;
			}
		} // for col
	} // for c
	bcubes.erase(bcubes.begin()+buildings_end, bcubes.begin()+bcubes_coll_end); // erase colliders that were added above from bcubes
	// assign handicap spots
	unsigned const num_hcap_spots((hcap_cands.size() + 10)/20); // 5% of total spots, rounded to the center

	if (num_hcap_spots > 0) {
		sort(hcap_cands.begin(), hcap_cands.end());

		for (unsigned n = 0; n < num_hcap_spots; ++n) {
			hcap_groups.add_obj(hcap_space_t(hcap_cands[n]), hcaps);
			assert(hcap_cands[n].pspace_ix < pspaces.size());
			pspaces[hcap_cands[n].pspace_ix].is_hcap = 1;
		}
	}
	return has_parking;
}

// non-const because this sets driveway_t::car_ix through add_car()
void city_obj_placer_t::add_cars_to_driveways(vector<car_t> &cars, vector<road_plot_t> const &plots, vector<vect_cube_t> &plot_colliders, unsigned city_id, rand_gen_t &rgen) {
	car_t car;
	car.park();
	car.cur_city = city_id;
	car.cur_road_type = TYPE_DRIVEWAY;
	vector3d const nom_car_size(city_params.get_nom_car_size()); // {length, width, height}

	for (auto i = driveways.begin(); i != driveways.end(); ++i) {
		if (rgen.rand_float() < 0.5) continue; // no car in this driveway 50% of the time
		car.cur_road = (unsigned short)i->plot_ix; // store plot_ix in road field
		car.cur_seg  = (unsigned short)(i - driveways.begin()); // store driveway index in cur_seg
		cube_t const &plot(plots[i->plot_ix]);
		car.dim = (i->y1() == plot.y1() || i->y2() == plot.y2()); // check which edge of the plot the driveway is connected to, which is more accurate than the aspect ratio
		if (i->get_sz_dim(car.dim) < 1.6*nom_car_size.x || i->get_sz_dim(!car.dim) < 1.25*nom_car_size.y) continue; // driveway is too small to fit this car
		car.dir = rgen.rand_bool(); // randomly pulled in vs. backed in, since we don't know the direction to the house anyway
		float const pad_l(0.75*nom_car_size.x), pad_w(0.6*nom_car_size.y); // needs to be a bit larger to fit trucks
		point cpos(0.0, 0.0, i->z2());
		cpos[ car.dim] = rgen.rand_uniform(i->d[ car.dim][0]+pad_l, i->d[ car.dim][1]-pad_l);
		cpos[!car.dim] = rgen.rand_uniform(i->d[!car.dim][0]+pad_w, i->d[!car.dim][1]-pad_w); // not quite centered
		car.set_bcube(cpos, nom_car_size);
		// check if this car intersects another parked car; this can only happen if two driveways intersect, which should be rare
		bool intersects(0);

		for (auto c = cars.rbegin(); c != cars.rend(); ++c) {
			if (c->cur_road != i->plot_ix) break; // prev plot, done
			if (car.bcube.intersects(c->bcube)) {intersects = 1; break;}
		}
		if (intersects) continue; // skip
		i->in_use = 2; // permanently in use
		cars.push_back(car);
		plot_colliders[i->plot_ix].push_back(car.bcube); // prevent pedestrians from walking through this parked car
	} // for i
}

int city_obj_placer_t::select_dest_parking_space(unsigned driveway_ix, bool allow_hcap, bool reserve_spot, float car_len, rand_gen_t &rgen) const {
	if (pspaces.empty()) return -1; // no parking spaces; error?
	assert(driveway_ix < driveways.size());
	driveway_t const &driveway(driveways[driveway_ix]);
	int const park_lot_ix(driveway.park_lot_ix);
	if (park_lot_ix < 0) return -1; // error?
	parking_lot_t const &parking_lot(get_parking_lot(park_lot_ix));
	bool const rdir(parking_lot.dir ^ parking_lot.row_dir);
	unsigned const entrance_row_ix(rdir ? 0 : parking_lot.num_rows-1);
	int const next_space_inc((rdir ? 1 : -1) * parking_lot.row_sz); // skip by columns
	float const min_dw_pad(1.1*car_len); // for backing up out of parking space
	// select an available spot on the row next to the driveway
	vector<unsigned> avail_spaces;

	for (unsigned i = 0; i < pspaces.size(); ++i) {
		parking_space_t const &pi(pspaces[i]);
		if ((int)pi.p_lot_ix != park_lot_ix || !pi.is_avail()) continue;
		if (!allow_hcap && pi.is_hcap)    continue;
		if (pi.row_ix != entrance_row_ix) continue; // can only enter along the first row, which is connected to the driveway
		// check driveway padding around space for pulling in and backing up
		float const centerline(pi.center[driveway.dim]);
		if ((driveway.d[driveway.dim][1] - centerline) < min_dw_pad || (centerline - driveway.d[driveway.dim][0]) < min_dw_pad) continue;
		unsigned psix(i);
		bool col_has_car(0);

		// pull all the way forward
		for (int j = int(i)+next_space_inc; j >= 0 && j < (int)pspaces.size(); j += next_space_inc) {
			parking_space_t const &pj(pspaces[j]);
			if ((int)pj.p_lot_ix != park_lot_ix) break; // end of parking lot, done
			assert(pj.col_ix == pi.col_ix);
			if (pj.has_active_car) {col_has_car = 1; break;}
			if (!pj.is_avail()) break; // blocked, done
			psix = j; // move up one space
		}
		if (col_has_car) continue; // don't block another car in; this is only needed because parking lots don't have space between columns
		avail_spaces.push_back(psix);
	} // for i
	if (avail_spaces.empty()) return -1; // no space found
	unsigned const sel_space(avail_spaces[rgen.rand() % avail_spaces.size()]); // select a random space
	if (reserve_spot) {pspaces[sel_space].add_car();} // mark space as occupied; this is non-const
	return sel_space;
}

bool check_pt_and_place_blocker(point const &pos, vect_cube_t &blockers, float radius, float blocker_spacing, bool add_blocker=1) {
	cube_t bc(pos);
	if (has_bcube_int_xy(bc, blockers, radius)) return 0; // intersects a building, parking lot, or other object - skip
	if (!add_blocker) return 1;
	bc.expand_by_xy(blocker_spacing);
	blockers.push_back(bc); // prevent trees and benches from being too close to each other
	return 1;
}
bool try_place_obj(cube_t const &plot, vect_cube_t &blockers, rand_gen_t &rgen, float radius, float blocker_spacing, unsigned num_tries, point &pos, bool add_blocker) {
	for (unsigned t = 0; t < num_tries; ++t) {
		pos = rand_xy_pt_in_cube(plot, radius, rgen);
		if (check_pt_and_place_blocker(pos, blockers, radius, blocker_spacing, add_blocker)) return 1; // success
	}
	return 0;
}
void place_tree(point const &pos, float radius, int ttype, vect_cube_t &colliders, vector<point> *tree_pos,
	bool allow_bush, bool add_bush, bool is_sm_tree, bool has_planter, float custom_size=0.0, float pine_xy_sz=1.0)
{
	tree_placer.add(pos, custom_size, ttype, allow_bush, add_bush, is_sm_tree, pine_xy_sz); // use same tree type
	// use 15% of the placement radius for collision (trunk + planter), smaller if no planter
	cube_t bcube;
	bcube.set_from_sphere(pos, (has_planter ? 0.15 : 0.05)*radius);
	bcube.z2() += max(radius, 0.25f*city_params.road_width); // increase cube height; make sure it's taller than people
	colliders.push_back(bcube);
	if (tree_pos != nullptr) {tree_pos->push_back(pos);}
}
void resize_blockers_for_trees(vect_cube_t &blockers, unsigned six, unsigned eix, float resize_amt) {
	float const height_thresh(1.5*city_params.get_nom_car_size().z);

	for (auto i = blockers.begin()+six; i != blockers.begin()+eix; ++i) {
		float const resize_scale((i->dz() < height_thresh) ? 1.0 : 0.5); // reduced shrink for tall objects such as fences and walls since trees may clip through them
		i->expand_by_xy(resize_scale*resize_amt);
	}
}

template<typename T> bool intersects_city_obj(cube_t const &c, vector<T> const &objs, cube_t const &exclude=cube_t()) {
	for (T const &i : objs) {if (i.bcube != exclude && i.bcube.intersects(c)) return 1;}
	return 0;
}
template<typename T> bool intersects_city_obj_xy(cube_t const &c, vector<T> const &objs) {
	for (T const &i : objs) {if (i.bcube.intersects_xy(c)) return 1;}
	return 0;
}
bool city_obj_placer_t::check_walkway_coll_xy(point const &pos, float radius) const {
	cube_t test_cube(pos);
	test_cube.expand_by_xy(radius);
	return (intersects_city_obj_xy(test_cube, walkways) || intersects_city_obj_xy(test_cube, elevators));
}

void city_obj_placer_t::place_trees_in_plot(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders,
	vector<point> &tree_pos, vect_cube_t const &plot_cuts, rand_gen_t &rgen, unsigned buildings_end)
{
	if (city_params.max_trees_per_plot == 0) return;
	float const radius(city_params.tree_spacing*city_params.get_nom_car_size().x); // in multiples of car length
	float const spacing(max(radius, get_min_obj_spacing())), radius_exp(2.0*spacing);
	vector3d const plot_sz(plot.get_size());
	if (min(plot_sz.x, plot_sz.y) < 2.0*radius_exp) return; // plot is too small for trees of this size
	unsigned num_trees(city_params.max_trees_per_plot);
	if (plot.is_park) {num_trees += (rgen.rand() % city_params.max_trees_per_plot);} // allow up to twice as many trees in parks
	assert(buildings_end <= blockers.size());
	// shrink non-building blockers (parking lots, driveways, fences, walls, hedges) to allow trees to hang over them; okay if they become denormalized
	unsigned const input_blockers_end(blockers.size());
	float const non_buildings_overlap(0.7*radius);
	resize_blockers_for_trees(blockers, buildings_end, input_blockers_end, -non_buildings_overlap);
	bool const has_planter(plot.is_commercial()); // only commercial trees

	for (unsigned n = 0; n < num_trees; ++n) {
		bool const is_sm_tree((rgen.rand()%3) == 0); // 33% of the time is a pine/palm tree
		int ttype(-1); // Note: okay to leave at -1; also, don't have to set to a valid tree type
		if (is_sm_tree) {ttype = (plot.is_park ? (rgen.rand()&1) : 2);} // pine/short pine in parks, palm in city blocks
		else {ttype = rgen.rand()%100;} // random type
		bool const is_palm(is_sm_tree && ttype == 2);
		bool const allow_bush(plot.is_park && max_unique_trees == 0); // can't place bushes if tree instances are enabled (generally true) because bushes may be instanced in non-parks
		bool const add_bush(0); // not yet supported
		float const bldg_extra_radius(is_palm ? 0.5f*radius : 0.0f); // palm trees are larger and must be kept away from buildings, but can overlap with other trees
		float const pine_xy_sz((is_sm_tree && plot.is_park) ? rgen.rand_uniform(0.5, 0.8) : 1.0); // randomly narrower
		float const coll_radius(spacing + bldg_extra_radius);
		point pos;
		if (!try_place_obj(plot, blockers, rgen, coll_radius, (radius - bldg_extra_radius), 10, pos, 1)) continue; // 10 tries per tree, extra spacing for palm trees
		if (point_in_cubes_xy_exp(plot_cuts, pos, radius_exp))  continue; // no retry
		// check walkways; waklway elevators haven't been placed yet for this plot, so add extra padding
		if (check_walkway_coll_xy(pos, (coll_radius + radius))) continue; // no retry
		// size is randomly selected by the tree generator using default values; allow bushes in parks
		place_tree(pos, radius, ttype, colliders, &tree_pos, allow_bush, add_bush, is_sm_tree, has_planter, 0.0, pine_xy_sz);
		if (plot.is_park) continue; // skip row logic and just place trees randomly throughout the park
		// now that we're here, try to place more trees at this same distance from the road in a row
		bool const dim(min((pos.x - plot.x1()), (plot.x2() - pos.x)) < min((pos.y - plot.y1()), (plot.y2() - pos.y)));
		bool const dir((pos[dim] - plot.d[dim][0]) < (plot.d[dim][1] - pos[dim]));
		float const step(1.25*radius_exp*(dir ? 1.0 : -1.0)); // positive or negative (must be > 2x radius spacing)
					
		for (; n < city_params.max_trees_per_plot; ++n) {
			pos[dim] += step;
			if (pos[dim] < plot.d[dim][0]+radius || pos[dim] > plot.d[dim][1]-radius) break; // outside place area
			if (!check_pt_and_place_blocker(pos, blockers, coll_radius, (spacing - bldg_extra_radius))) continue; // placement failed
			if (point_in_cubes_xy_exp(plot_cuts, pos, radius_exp)) continue;
			if (check_walkway_coll_xy(pos, coll_radius))           continue; // hit walkway
			place_tree(pos, radius, ttype, colliders, &tree_pos, allow_bush, add_bush, is_sm_tree, has_planter); // use same tree type
		} // for n
	} // for n
	resize_blockers_for_trees(blockers, buildings_end, input_blockers_end, non_buildings_overlap); // undo initial expand
}

void city_obj_groups_t::clear() {
	vector<cube_with_ix_t>::clear();
	by_tile.clear();
	bcube.set_to_zeros();
}
void city_obj_groups_t::insert_obj_ix(cube_t const &c, unsigned ix) {
	by_tile[get_tile_id_for_cube(c)].push_back(ix);
}
template<typename T> void city_obj_groups_t::add_obj(T const &obj, vector<T> &objs) {
	insert_obj_ix(obj.bcube, objs.size());
	objs.push_back(obj);
}
template<typename T> void city_obj_groups_t::create_groups(vector<T> &objs, cube_t &all_objs_bcube) {
	vector<cube_with_ix_t>::clear();
	bcube.set_to_zeros();
	vector<T> new_objs;
	new_objs.reserve(objs.size());
	reserve(by_tile.size()); // the number of actual groups

	for (auto g = by_tile.begin(); g != by_tile.end(); ++g) {
		unsigned const group_start(new_objs.size());
		cube_with_ix_t group;

		for (auto i = g->second.begin(); i != g->second.end(); ++i) {
			assert(*i < objs.size());
			group.assign_or_union_with_cube(objs[*i].get_outer_bcube());
			new_objs.push_back(objs[*i]);
		}
		sort(new_objs.begin()+group_start, new_objs.end());
		group.ix = new_objs.size();
		push_back(group);
		bcube.assign_or_union_with_cube(group);
	} // for g
	all_objs_bcube.assign_or_union_with_cube(bcube);
	objs.swap(new_objs);
	by_tile.clear(); // no longer needed
}

void add_cube_to_colliders_and_blockers(cube_t const &cube, vect_cube_t &blockers, vect_cube_t &colliders) {
	colliders.push_back(cube);
	blockers .push_back(cube);
}
vect_bird_place_t *select_bird_loc_dest(bool add_pigeons, bool add_birds, vect_bird_place_t &pigeon_locs, vect_bird_place_t &bird_locs, rand_gen_t &rgen) {
	if (!add_pigeons && !add_birds) return nullptr; // error?
	if (!add_pigeons) return &bird_locs; // always a bird
	if (rgen.rand_float() < 0.25) return &pigeon_locs; // add pigeon 25% of the time
	if (add_birds) return &bird_locs;
	return nullptr;
}
template<typename T> void add_bird_loc(T const &obj, vect_bird_place_t &dest, rand_gen_t &rgen) {
	dest.add_placement_top_center(obj.get_bird_bcube(), rgen);
}
template<> void add_bird_loc(mailbox_t const &obj, vect_bird_place_t &dest, rand_gen_t &rgen) {
	// for mailboxes, start the bird facing toward the road so that it doesn't fly into a house
	dest.emplace_back(cube_top_center(obj.get_bird_bcube()), obj.dim, obj.dir, 0); // use_orient=0
}
template<typename T> void add_objs_top_center(T const &objs, unsigned start_ix, bool add_pigeons, bool add_birds,
	vect_bird_place_t &pigeon_locs, vect_bird_place_t &bird_locs, rand_gen_t &rgen)
{
	for (auto i = objs.begin()+start_ix; i != objs.end(); ++i) {
		vect_bird_place_t *const dest(select_bird_loc_dest(add_pigeons, add_birds, pigeon_locs, bird_locs, rgen));
		if (dest != nullptr) {add_bird_loc(*i, *dest, rgen);}
	}
}

void choose_edge_pos(cube_t const &region, float border, bool dim, bool dir, point &pos, rand_gen_t &rgen) {
	pos[ dim] = region.d[dim][dir];
	pos[!dim] = rgen.rand_uniform(region.d[!dim][0]+border, region.d[!dim][1]-border);
}
bool check_path_coll_xy(cube_t const &c, vector<park_path_t> const &paths, unsigned paths_start) {
	for (auto p = paths.begin()+paths_start; p != paths.end(); ++p) {
		if (p->check_cube_coll_xy(c)) return 1;
	}
	return 0;
}
bool check_path_tree_coll(park_path_t const &path, vector<point> const &tree_pos) {
	for (point const &pos : tree_pos) { // check for collisions with tree trunks
		cube_t bc; bc.set_from_sphere(pos, 0.1*path.hwidth); // small size
		if (path.check_cube_coll_xy(bc)) return 1;
	}
	return 0;
}
float get_power_pole_height() {return 0.9*city_params.road_width;}

// this version is for checking against blockers placed in a previous step
bool is_placement_blocked(cube_t const &cube, vect_cube_t const &blockers, cube_t const &exclude, unsigned prev_blockers_end, float expand=0.0, bool exp_dim=0) {
	cube_t query_cube(cube);
	query_cube.expand_in_dim(exp_dim, expand);

	for (auto b = blockers.begin(); b != blockers.begin()+prev_blockers_end; ++b) {
		if (*b != exclude && b->intersects_xy_no_adj(query_cube)) return 1;
	}
	return 0;
}
// this version is for checking against blockers placed in the current step
bool is_placement_blocked_recent(cube_t const &cube, vect_cube_t const &blockers, unsigned blockers_start) {
	assert(blockers_start <= blockers.size());

	for (auto b = blockers.begin()+blockers_start; b != blockers.end(); ++b) {
		if (b->intersects_xy_no_adj(cube)) return 1;
	}
	return 0;
}
bool check_close_to_door(point const &pos, float min_dist, unsigned building_ix) {
	point door_pos;
	return (get_building_door_pos_closest_to(building_ix, pos, door_pos, 1) && dist_xy_less_than(pos, door_pos, min_dist)); // close to door; inc_garage_door=1
}

// Note: blockers are used for placement of objects within this plot; colliders are used for pedestrian AI
void city_obj_placer_t::place_detail_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders,
	vector<point> const &tree_pos, vect_cube_t const &pond_blockers, rand_gen_t &rgen, bool have_streetlights)
{
	bool const is_residential(plot.is_residential), is_park(plot.is_park);
	float const car_length(city_params.get_nom_car_size().x); // used as a size reference for other objects
	float const min_obj_spacing(get_min_obj_spacing()), sidewalk_width(get_sidewalk_width());
	unsigned const benches_start(benches.size()), trashcans_start(trashcans.size()), substations_start(sstations.size());
	unsigned const fountains_start(fountains.size()), ppoles_start(ppoles.size()), paths_start(ppaths.size());

	// place paths in parks
	if (is_park) {
		//highres_timer_t timer("Add Park Paths"); // < 1ms
		unsigned const num_paths(1 + (rgen.rand_float() < 0.67)); // 1-2
		unsigned const num_path_segs = 100;
		float const path_hwidth(0.08*city_params.road_width), path_height(0.01*path_hwidth);
		float const plot_min_edge(min(plot.dx(), plot.dy())), edge_border(max(2.0f*path_hwidth, 0.2f*plot_min_edge));

		if (plot_min_edge > 3.0*edge_border) { // park has enough space for a path; should always be true
			point start, end;
			start.z = end.z = plot.z2() + path_height;
			bool dim(rgen.rand_bool());

			for (unsigned n = 0; n < num_paths; ++n, dim ^= 1) { // alternate dims for each path
				for (unsigned N = 0; N < 100; ++N) { // make 100 tries
					choose_edge_pos(plot, edge_border, dim, 0, start, rgen); // choose starting point
					choose_edge_pos(plot, edge_border, dim, 1, end,   rgen); // choose ending point on the opposite edge
					float const path_curve(rgen.rand_float()*20.0*path_hwidth);
					float const sine_mult((1 + (rgen.rand()%3))*TWO_PI); // 1-3 cycles
					vector3d const seg_delta((end - start)/num_path_segs);
					cube_t valid_region(plot);
					valid_region.expand_in_dim(!dim, -2.0f*path_hwidth); // shrink to keep path inside the park on the opposite edges
					park_path_t path(path_hwidth, GRAY, plot);
					point cur(start);
					start[dim] -= path_hwidth; // extend outside the plot to avoid a gap; will be clipped during drawing
					end  [dim] += path_hwidth;
					path.pts.push_back(start);

					for (unsigned s = 0; s+1 < num_path_segs; ++s) {
						cur += seg_delta;
						point p(cur);
						float const t(float(s+1)/num_path_segs), tc(2.0*fabs(t - 0.5)); // t is the position along the path in [0.0, 1.0]
						p[!dim] += path_curve * sin(sine_mult*t) * (1.0 - tc*tc*tc); // sine wave in cubic envelope
						valid_region.clamp_pt_xy(p);
						path.pts.push_back(p);
					} // for s
					path.pts.push_back(end);
					path.calc_bcube_bsphere();
					if (check_path_tree_coll(path, tree_pos)) continue;
					ppath_groups.add_obj(path, ppaths);
					break; // success
				} // for N
			} // for n
		}
		if (1) { // try to place pond(s)
			float const pond_border(max(sidewalk_width, path_hwidth));
			vect_cube_t active_pond_blockers;

			for (cube_t const &pb : pond_blockers) {
				if (pb.intersects_xy(plot)) {active_pond_blockers.push_back(pb);}
			}
			for (unsigned n = 0; n < 100; ++n) { // 100 tries to place a pond
				float const sz(city_params.road_width*rgen.rand_uniform(0.5, 1.0));
				vector3d pond_sz; // radius
				for (unsigned d = 0; d < 2; ++d) {pond_sz[d] = sz*rgen.rand_uniform(1.0, 1.5);}
				cube_t pond_area(plot);
				pond_area.expand_by(-pond_sz);
				pond_area.expand_by_xy(-pond_border);
				if (pond_area.dx() <= 0.0 || pond_area.dy() <= 0.0) continue; // not enough space; shouldn't fail for reasonable road vs. plot sizes
				point center(0.0, 0.0, plot.z2());
				for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(pond_area.d[d][0], pond_area.d[d][1]);}
				cube_t pond(center);
				pond.expand_by(pond_sz);
				bool blocked(0);

				for (auto p = ppaths.begin()+paths_start; p != ppaths.end(); ++p) { // check paths
					if (p->check_cube_coll_xy(pond)) {blocked = 1; break;}
				}
				if (blocked) continue;
				if (has_bcube_int_xy(pond, active_pond_blockers)) continue; // check underground basement rooms
				float const depth(city_params.road_width*rgen.rand_uniform(0.1, 0.5));
				pond_t const pond_obj(center, pond_sz.x, pond_sz.y, depth);

				for (point const &p : tree_pos) { // check trees; trees on the edge are okay
					if (pond.contains_pt_xy(p) && pond_obj.point_contains_xy(p)) {blocked = 1; break;}
				}
				if (blocked) continue;
				pond_groups.add_obj(pond_obj, ponds);
				add_cube_to_colliders_and_blockers(pond, colliders, blockers);
				max_eq(pond_max_depth, depth);
				break; // success
			} // for n
		}
		// place picnic tables
		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_PICNIC)) {
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_PICNIC)); // W, D, H
			float const pt_len(0.42*car_length), pt_width(pt_len*sz.y/sz.x), pt_height(pt_len*sz.z/sz.x);

			for (unsigned n = 0; n < city_params.max_benches_per_plot; ++n) { // use the same max count as benches
				point pos;
				if (!try_place_obj(plot, blockers, rgen, 0.5*max(pt_len, pt_width), 0.0, 1, pos, 0)) continue; // 1 try
				picnic_t const pt(pos, pt_height, rgen.rand_bool(), 0); // random dim; dir=0 since picnic tables are symmetric
				if (check_path_coll_xy(pt.bcube, ppaths, paths_start)) continue; // check park path collision
				picnic_groups.add_obj(pt, picnics);
				add_cube_to_colliders_and_blockers(pt.bcube, colliders, blockers);
				blockers.back().expand_by_xy(0.5*pt_width); // add extra padding to the sides and ends
			} // for n
		}
	} // end is_park
	if (!walkways.empty()) {
		// place vertical pillars supporting walkways connecting buildings, and elevators leading up to walkways
		cube_t pillar_area(plot);
		pillar_area.expand_by_xy(-2.0*sidewalk_width); // don't block sidewalks or nearby areas

		for (walkway_t &w : walkways) {
			if (!w.bcube.intersects_xy(plot)) continue;
			float const length(w.get_length()), width(w.get_width()), height(w.bcube.z1() - plot.z2());
			bool const to_skyway(w.open_ends[0] || w.open_ends[1]);
			if (!to_skyway && (length < 8.0*width || length < 1.6*height)) continue; // too short or high for a support pillar/elevator; skyway walkways are always supported
			bool const is_concrete(height < 0.25*length); // concrete when shorter, steel when taller
			float const pillar_hlen((is_concrete ? 0.1 : 0.15)*width), pillar_hwid((is_concrete ? 0.3 : 0.15)*width);
			cube_t pillar;
			set_wall_width(pillar, w.bcube.get_center_dim(!w.dim), pillar_hwid, !w.dim);
			float const pos_offset[3] = {0.5, 0.333, 0.667}; // try center and one third from each end

			for (unsigned n = 0; n < 3; ++n) {
				float len_pos(w.bcube.d[w.dim][0] + pos_offset[n]*length);

				if (n == 0 && len_pos > plot.d[w.dim][0] && len_pos < plot.d[w.dim][1]) { // center point is in the plot; clamp to pillar_area
					max_eq(len_pos, pillar_area.d[w.dim][0]+pillar_hlen);
					min_eq(len_pos, pillar_area.d[w.dim][1]-pillar_hlen);
				}
				set_wall_width(pillar, len_pos, pillar_hlen, w.dim);
				if (!pillar_area.contains_cube_xy(pillar)) continue; // not contained in plot interior
				set_cube_zvals(pillar, plot.z2(), w.bcube.z1());
				cube_t pillar_exp(pillar);
				pillar_exp.expand_by_xy(min_obj_spacing);
				if (has_bcube_int_no_adj(pillar_exp, blockers)) continue; // skip the walkway; what else can this intersect, only parking lots?
				if (intersects_city_obj(pillar_exp, elevators) || intersects_city_obj(pillar_exp, walkways, w.bcube)) continue; // exclude ourself
				pillar_groups.add_obj(pillar_t(pillar, is_concrete), pillars);
				add_cube_to_colliders_and_blockers(pillar, colliders, blockers);

				// maybe add walkway elevator
				float const e_width(2.0*sidewalk_width), e_depth(0.9*e_width), e_clearance(0.9*e_depth);
				bool const first_side(rgen.rand_bool());

				for (unsigned side = 0; side < 2; ++side) {
					bool const dir(bool(side) ^ first_side ^ 1);
					cube_t ebc;
					set_cube_zvals(ebc, plot.z2(), w.bcube.z2());
					set_wall_width(ebc, len_pos, 0.5*e_width, w.dim); // set width, parallel to pillar
					ebc.d[!w.dim][0] = ebc.d[!w.dim][1] = w.bcube.d[!w.dim][dir]; // edge of walkway
					ebc.d[!w.dim][dir] += (dir ? 1.0 : -1.0)*e_depth; // set depth
					cube_t ebc_exp(ebc);
					ebc_exp.d[!w.dim][dir] += (dir ? 1.0 : -1.0)*e_clearance; // add clearance in front of door
					if (!pillar_area.contains_cube_xy(ebc_exp))  continue; // not contained in plot interior
					if (has_bcube_int_no_adj(ebc_exp, blockers)) continue;
					if (intersects_city_obj(ebc_exp, elevators) || intersects_city_obj(ebc_exp, walkways, w.bcube)) continue; // exclude ourself
					ww_elevator_t const elevator(ebc, !w.dim, dir, w.floor_spacing, w.bcube);
					w.attach_elevator(elevator);
					wwe_groups.add_obj(elevator, elevators);
					add_cube_to_colliders_and_blockers(ebc, colliders, blockers);
					break; // only one side needed
				} // for side
				break; // success
			} // for n
		} // for w
	}
	// place fire_hydrants if the model has been loaded; don't add fire hydrants in parks
	if (!is_park && building_obj_model_loader.is_model_valid(OBJ_MODEL_FHYDRANT)) {
		// we want the fire hydrant on the edge of the sidewalk next to the road, not next to the plot; this makes it outside the plot itself
		float const radius(0.04*car_length), height(0.18*car_length), dist_from_road(-0.5*radius - sidewalk_width);
		point pos(0.0, 0.0, plot.z2()); // XY will be assigned below

		for (unsigned dim = 0; dim < 2; ++dim) {
			pos[!dim] = plot.get_center_dim(!dim);

			for (unsigned dir = 0; dir < 2; ++dir) {
				pos[dim] = plot.d[dim][dir] - (dir ? 1.0 : -1.0)*dist_from_road; // move into the sidewalk along the road
				// Note: will skip placement if too close to a previously placed tree, but that should be okay as it is relatively uncommon
				if (!check_pt_and_place_blocker(pos, blockers, radius, 2.0*radius)) continue; // bad placement, skip
				vector3d orient(zero_vector);
				orient[!dim] = (dir ? 1.0 : -1.0); // oriented perpendicular to the road
				fire_hydrant_t const fire_hydrant(pos, radius, height, orient);
				fhydrant_groups.add_obj(fire_hydrant, fhydrants);
				colliders.push_back(fire_hydrant.bcube);
			} // for dir
		} // for dim
	}
	// place dumpsters in city blocks
	if (!plot.is_residential_not_park() && have_buildings() && building_obj_model_loader.is_model_valid(OBJ_MODEL_DUMPSTER)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_DUMPSTER)); // W, D, H
		float const height(0.3*car_length), width(height*sz.y/sz.z), depth(height*sz.x/sz.z);
		point const center(plot.get_cube_center());
		vect_cube_with_ix_t bcubes; // we need the building index for the get_building_door_pos_closest_to() call
		get_building_bcubes(plot, bcubes);
		
		for (cube_with_ix_t const &b : bcubes) {
			if (rgen.rand_float() < 0.1) continue; // skip this building 10% of the time; placements often fail anyway, and usually fail for non-cube buildings
			// find the bcube side furthest from the road/closest to the plot center
			bool dim(0), dir(0);
			float dmin_sq(0.0);

			for (unsigned d = 0; d < 2; ++d) {
				point edge_center;
				edge_center[!d] = b.get_center_dim(!d);

				for (unsigned e = 0; e < 2; ++e) {
					edge_center[d] = b.d[d][e];
					float const dsq(p2p_dist_xy_sq(edge_center, center));
					if (dmin_sq == 0.0 || dsq < dmin_sq) {dim = d; dir = e; dmin_sq = dsq;}
				}
			} // for d
			float const lo(b.d[!dim][0] + 0.7*width), hi(b.d[!dim][1] - 0.7*width);
			if (lo >= hi) continue; // wall too short to place a dumpster
			float const wall_pos(b.d[dim][dir]);
			point pos;
			pos.z     = plot.z1();
			pos[ dim] = wall_pos + (dir ? 1.0 : -1.0)*0.65*depth; // move away from the wall
			pos[!dim] = rgen.rand_uniform(lo, hi);
			dumpster_t const dumpster(pos, height, dim, dir);
			cube_t check_cube(dumpster.bcube);
			check_cube.d[dim][dir] += (dir ? 1.0 : -1.0)*2.0*depth; // add extra clearance in front
			if (is_placement_blocked(check_cube, blockers, b, blockers.size())) continue; // check blockers; no expand
			point const center(dumpster.bcube.get_cube_center()); // pos shifted up
			point p_include(center);
			p_include[dim] = wall_pos - (dir ? 1.0 : -1.0)*1.0*depth; // inside the wall
			if (!check_sphere_coll_building(p_include, 0.0, 0, b.ix)) continue; // *must* be next to an exterior wall; radius=0.0, xy_only=0
			if (check_close_to_door(pos, 2.0*width, b.ix))            continue;
			dumpster_groups.add_obj(dumpster, dumpsters);
			add_cube_to_colliders_and_blockers(dumpster.bcube, colliders, blockers);
			blockers.back().d[dim][dir] = check_cube.d[dim][dir]; // include clearance
		} // for b
	}
	// place fountains in parks and 25% of the time in city blocks
	if ((is_park || (!is_residential && (rgen.rand() & 3) == 0)) && building_obj_model_loader.is_model_valid(OBJ_MODEL_FOUNTAIN)) {
		float const radius(0.35 * car_length), spacing(max(1.5f*radius, min_obj_spacing));
		cube_t place_area(plot);
		place_area.expand_by_xy(-sidewalk_width); // not too close to sidewalks
		point pos;

		if (try_place_obj(place_area, blockers, rgen, (is_park ? 1.0 : 1.5)*radius, spacing, 10, pos, 0)) { // 10 tries
			float const height(2.0*radius); // seems about right for the current fountain models
			fountain_t const fountain(pos, height, rgen.rand()); // random model_select

			if (!is_park || !check_path_coll_xy(fountain.bcube, ppaths, paths_start)) { // check park path collision
				fountain_groups.add_obj(fountain, fountains);
				// don't place blockers until we've added benches, since we want to allow benches closer to the fountains
			}
		}
	}
	// place benches in parks and non-residential areas, and next to fountains
	if (!plot.is_residential_not_park()) {
		float const bench_radius(0.3 * car_length), bench_spacing(max(bench_radius, 1.5f*min_obj_spacing)); // add a bit of extra space

		for (unsigned n = 0; n < city_params.max_benches_per_plot; ++n) {
			point pos;
			if (!try_place_obj(plot, blockers, rgen, bench_spacing, 0.0, 1, pos, 0)) continue; // 1 try
			bool bench_dim(0), bench_dir(0);
			float dmin(0.0);

			for (unsigned dim = 0; dim < 2; ++dim) {
				for (unsigned dir = 0; dir < 2; ++dir) {
					float const dist(fabs(pos[dim] - plot.d[dim][dir])); // find closest distance to road (plot edge) and orient bench that way
					if (dmin == 0.0 || dist < dmin) {bench_dim = !dim; bench_dir = !dir; dmin = dist;}
				}
			}
			bench_t const bench(pos, bench_radius, bench_dim, bench_dir);
			if (is_park && check_path_coll_xy(bench.bcube, ppaths, paths_start)) continue; // check park path collision
			bench_groups.add_obj(bench, benches);
			add_cube_to_colliders_and_blockers(bench.bcube, colliders, blockers);
			blockers.back().expand_in_dim(!bench.dim, 0.25*bench.radius); // add extra padding in front (for seat access) and back (which extends outside bcube)
		} // for n
		cube_t place_area(plot);
		place_area.expand_by_xy(-(bench_radius + sidewalk_width));

		for (auto f = fountains.begin()+fountains_start; f != fountains.end(); ++f) {
			bool const dim(rgen.rand_bool()); // add benches in a random dir around the fountain

			for (unsigned dir = 0; dir < 2; ++dir) {
				point pos;
				pos.z     = plot.z1();
				pos[ dim] = f->bcube.d[dim][dir] + (dir ? 1.0 : -1.0)*1.5*bench_radius;
				pos[!dim] = f->bcube.get_center_dim(!dim);
				if (!place_area.contains_pt_xy(pos) || has_bcube_int_xy(cube_t(pos), blockers, bench_spacing)) continue; // excludes fountains
				bench_t const bench(pos, bench_radius, !dim, dir); // face toward the fountain
				if (is_park && check_path_coll_xy(bench.bcube, ppaths, paths_start)) continue; // check park path collision
				bench_groups.add_obj(bench, benches);
				add_cube_to_colliders_and_blockers(bench.bcube, colliders, blockers);
				blockers.back().expand_in_dim(!bench.dim, 0.25*bench.radius); // add extra padding in front (for seat access) and back (which extends outside bcube)
			} // for dir
		} // for f
	}
	for (auto f = fountains.begin()+fountains_start; f != fountains.end(); ++f) { // now add the fountain blockers
		add_cube_to_colliders_and_blockers(f->bcube, colliders, blockers);
	}
	// place planters; don't add planters in parks or residential areas
	if (plot.is_commercial()) {
		float const planter_height(0.05*car_length), planter_radius(0.25*car_length);

		for (auto i = tree_pos.begin(); i != tree_pos.end(); ++i) {
			planter_groups.add_obj(tree_planter_t(*i, planter_radius, planter_height), planters); // no colliders for planters; pedestrians avoid the trees instead
		}
	}
	// place power poles if there are houses or streetlights
	point corner_pole_pos(all_zeros);

	if ((is_residential && have_city_buildings()) || have_streetlights) {
		float const road_width(city_params.road_width), pole_radius(0.015*road_width), height(get_power_pole_height());
		float const xspace(plot.dx() + road_width), yspace(plot.dy() + road_width); // == city_params.road_spacing?
		float const xyspace[2] = {0.5f*xspace, 0.5f*yspace};
		// we can move in toward the plot center so that they don't block pedestrians, but then they can block driveways;
		// if we move them into the road, then they block traffic light crosswalks;
		// so we move them toward the road in an assymetic way and allow the pole to be not centered with the wires
		float const offset(0.075*road_width), extra_offset(get_power_pole_offset()); // assymmetric offset to match the aspect ratio of stoplights
		float const pp_x(plot.x2() + offset + extra_offset), pp_y(plot.y2() + offset);
		point pts[3]; // one on the corner and two on each side: {corner, x, y}
		for (unsigned i = 0; i < 3; ++i) {pts[i].assign(pp_x, pp_y, plot.z2());} // start at plot upper corner
		pts[1].x -= 0.5*xspace;
		pts[2].y -= 0.5*yspace;
		unsigned const dims[3] = {3, 1, 2};

		for (unsigned i = 0; i < 3; ++i) {
			point pos(pts[i]);
			float wires_offset(0.0);

			if (i > 0 && !driveways.empty()) { // check driveways for non-corner pole
				bool const dim(i == 2);
				float const prev_val(pos[dim]);
				move_to_not_intersect_driveway(pos, (2.0*pole_radius + sidewalk_width), dim);
				wires_offset = prev_val - pos[dim];
			}
			point base(pos);
			if (i == 1) {base.y += extra_offset;} // shift the pole off the sidewalk and off toward the road to keep it out of the way of pedestrians
			bool const at_line_end[2] = {0, 0};
			bool const at_grid_edge(plot.xpos+1U == num_x_plots || plot.ypos+1U == num_y_plots);
			ppole_groups.add_obj(power_pole_t(base, pos, pole_radius, height, wires_offset, xyspace, dims[i], at_grid_edge, at_line_end, is_residential), ppoles);
			if (i == 0) {corner_pole_pos = base;}
		} // for i
		if (plot.xpos == 0) { // no -x neighbor plot, but need to add the power poles there
			unsigned const pole_ixs[2] = {0, 2};

			for (unsigned i = 0; i < 2; ++i) {
				point pt(pts[pole_ixs[i]]);
				pt.x -= xspace;
				bool const at_line_end[2] = {1, 0};
				ppole_groups.add_obj(power_pole_t(pt, pt, pole_radius, height, 0.0, xyspace, dims[pole_ixs[i]], 1, at_line_end, is_residential), ppoles);
			}
		}
		if (plot.ypos == 0) { // no -y neighbor plot, but need to add the power poles there
			unsigned const pole_ixs[2] = {0, 1};

			for (unsigned i = 0; i < 2; ++i) {
				point pt(pts[pole_ixs[i]]);
				pt.y -= yspace;
				point base(pt);
				if (i == 1) {base.y += extra_offset;}
				bool const at_line_end[2] = {0, 1};
				ppole_groups.add_obj(power_pole_t(base, pt, pole_radius, height, 0.0, xyspace, dims[pole_ixs[i]], 1, at_line_end, is_residential), ppoles);
			}
		}
		if (plot.xpos == 0 && plot.ypos == 0) { // pole at the corner of the grid
			point pt(pts[0]);
			pt.x -= xspace;
			pt.y -= yspace;
			bool const at_line_end[2] = {1, 1};
			ppole_groups.add_obj(power_pole_t(pt, pt, pole_radius, height, 0.0, xyspace, dims[0], 1, at_line_end, is_residential), ppoles);
		}
		for (auto i = ppoles.begin()+ppoles_start; i != ppoles.end(); ++i) {colliders.push_back(i->get_ped_occluder());}
	}
	// place substations in commercial cities, near the corner pole that routes power into the ground, if the model has been loaded
	if (!is_residential && corner_pole_pos != all_zeros && building_obj_model_loader.is_model_valid(OBJ_MODEL_SUBSTATION)) {
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
		float const ss_height(0.08*city_params.road_width), dist_from_corner(0.12); // distance from corner relative to plot size
		vector3d ss_center((1.0 - dist_from_corner)*corner_pole_pos + dist_from_corner*plot.get_cube_center());
		ss_center.z -= 0.03*ss_height; // shift down slightly because the wire extends a bit below the base
		substation_t const ss(ss_center, ss_height, dim, dir);

		if (!has_bcube_int_xy(ss.bcube, blockers, max(0.2f*ss_height, min_obj_spacing))) { // skip if intersects a building or parking lot
			sstation_groups.add_obj(ss, sstations);
			add_cube_to_colliders_and_blockers(ss.bcube, colliders, blockers);
		}
	}
	// place trashcans next to sidewalks in commercial cities and parks; after substations so that we don't block them
	if (!plot.is_residential_not_park()) {
		float const tc_height(0.18*car_length), tc_radius(0.4*tc_height), dist_from_corner(0.06); // dist from corner relative to plot size

		for (unsigned d = 0; d < 4; ++d) { // try all 4 corners
			vector3d const tc_center((1.0 - dist_from_corner)*point(plot.d[0][d&1], plot.d[1][d>>1], plot.z2()) + dist_from_corner*plot.get_cube_center());
			trashcan_t const trashcan(tc_center, tc_radius, tc_height, 0); // is_cylin=0

			if (!has_bcube_int_xy(trashcan.bcube, blockers, max(1.5f*tc_radius, min_obj_spacing))) { // skip if intersects a building or parking lot, with some padding
				trashcan_groups.add_obj(trashcan, trashcans);
				add_cube_to_colliders_and_blockers(trashcan.bcube, colliders, blockers);
			}
		} // for d
	}
	// place newsracks along non-residential city streets
	if (!is_residential) {
		unsigned const NUM_NR_COLORS = 8;
		colorRGBA const nr_colors[NUM_NR_COLORS] = {WHITE, DK_BLUE, BLUE, ORANGE, RED, DK_GREEN, YELLOW, GRAY_BLACK};
		float const dist_from_corner(-0.03); // dist from corner relative to plot size; negative is outside the plot in the street area
		point pos(0, 0, plot.z2());

		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				if (rgen.rand_float() < 0.65) continue; // add to about a third of edges
				unsigned const num_this_side(1 + (rgen.rand() % 5)); // 1-5
				bool const place_together(num_this_side > 1 && rgen.rand_float() < 0.75); // 75% of the time
				unsigned place_region(0);
				pos[dim] = (1.0 - dist_from_corner)*plot.d[dim][dir] + dist_from_corner*plot.get_center_dim(dim); // set distance from the plot edge

				for (unsigned n = 0; n < num_this_side; ++n) {
					float const nr_height(0.28*car_length*rgen.rand_uniform(0.9, 1.11));
					float const nr_width(0.44*nr_height*rgen.rand_uniform(0.8, 1.25)), nr_depth(0.44*nr_height*rgen.rand_uniform(0.8, 1.25));
					float road_pos(0.0);
					if (n == 0 || !place_together) {place_region = (rgen.rand()&3);} // select a new placement region
					// streetlights are at 0.25 and 0.75, and telephone poles are at 0.5, so skip those ranges
					switch (place_region) {
					case 0: road_pos = rgen.rand_uniform(0.10, 0.20); break;
					case 1: road_pos = rgen.rand_uniform(0.30, 0.45); break;
					case 2: road_pos = rgen.rand_uniform(0.55, 0.70); break;
					case 3: road_pos = rgen.rand_uniform(0.80, 0.90); break;
					}
					pos[!dim] = plot.d[!dim][0] + road_pos*plot.get_sz_dim(!dim); // random pos along plot
					newsrack_t const newsrack(pos, nr_height, nr_width, nr_depth, dim, !dir, rgen.rand(), nr_colors[rgen.rand() % NUM_NR_COLORS]); // random style
					cube_t test_cube(newsrack.bcube);
					test_cube.d[dim][!dir] += (dir ? -1.0 : 1.0)*nr_depth; // add front clearance; faces inward from the plot
					if (has_bcube_int_xy(test_cube, blockers, max(0.1f*nr_width, min_obj_spacing))) continue; // skip if intersects a building or parking lot, with padding
					bool has_ppole_int(0);

					for (auto i = ppoles.begin()+ppoles_start; i != ppoles.end(); ++i) {
						if (sphere_cube_intersect_xy(i->get_base(), i->get_pole_radius(), test_cube)) {has_ppole_int = 1; break;}
					}
					if (has_ppole_int) continue; // should be rare
					nrack_groups.add_obj(newsrack, newsracks);
					colliders.push_back(newsrack.bcube); // no clearance
					blockers .push_back(test_cube); // includes clearance
				} // for n
			} // for dir
		} // for dim
	}
	// place manholes in adjacent roads, only on +x and +y sides; this will leave one edge road in each dim with no manholes
	if (1) {
		float const radius(0.125*car_length), dist_from_plot_edge(0.35*city_params.road_width); // centered in one lane of the road
		point pos(0.0, 0.0, plot.z2()); // XY will be assigned below
		bool const dir(0); // hard-coded for now

		for (unsigned dim = 0; dim < 2; ++dim) {
			pos[!dim] = plot.d[!dim][0] + 0.35*plot.get_sz_dim(!dim); // not centered
			pos[ dim] = plot.d[ dim][dir] + (dir ? 1.0 : -1.0)*dist_from_plot_edge; // move into the center of the road
			manhole_groups.add_obj(manhole_t(pos, radius), manholes); // Note: colliders not needed
		}
	}
	// maybe place a flag in a city commercial plot or a park
	if ((!is_residential && rgen.rand_float() < 0.3) || (is_park && rgen.rand_float() < 0.75)) { // 30% of the time for commerical plots, 75% of the time for parks
		float const length(0.25*city_params.road_width*rgen.rand_uniform(0.8, 1.25)), pradius(0.05*length);
		float const height((is_park ? 0.8 : 1.0)*city_params.road_width*rgen.rand_uniform(0.8, 1.25));
		float const spacing(is_park ? 2.0*pradius : 1.25*length); // parks have low items, so we only need to avoid colliding with the pole; otherwise need to check tall buildings
		cube_t place_area(plot);
		place_area.expand_by_xy(-0.5*city_params.road_width); // shrink slightly to keep flags away from power lines
		point base_pt;

		if (try_place_obj(place_area, blockers, rgen, spacing, 0.0, (is_park ? 5 : 20), base_pt, 1)) { // make up to 5/20 tries
			base_pt.z = plot.z2();
			cube_t pole(base_pt);
			pole.expand_by_xy(pradius);
			pole.z2() += height;

			if (check_walkway_coll_xy(base_pt, length)) {} // skip if collides with or is under a walkway
			else if (!is_park || !check_path_coll_xy(pole, ppaths, paths_start)) { // check park path collision
				bool dim(0), dir(0); // facing dir
				get_closest_dim_dir_xy(pole, plot, dim, dir); // face the closest plot edge
				flag_groups.add_obj(create_flag(dim, dir, base_pt, height, length), flags);
				colliders.push_back(pole); // only the pole itself is a collider
			}
		}
	}
	bool const add_pigeons(!is_residential && building_obj_model_loader.is_model_valid(OBJ_MODEL_PIGEON)); // only in cities with office buildings
	bool const add_birds(are_birds_enabled());

	if (add_pigeons || add_birds) {
		float const base_height(min(0.06f*car_length, 0.03f*min(plot.dx(), plot.dy()))), place_radius(4.0*base_height), obj_edge_spacing(0.25*base_height);
		// find all bird placements
		vect_bird_place_t pigeon_locs;

		// maybe place on benches, trashcans, substations, and fountains
		for (auto i = benches.begin()+benches_start; i != benches.end(); ++i) {
			if (i->bcube.get_sz_dim(!i->dim) <= 2.0*obj_edge_spacing) continue;
			vect_bird_place_t *const dest(select_bird_loc_dest(add_pigeons, add_birds, pigeon_locs, bird_locs, rgen));
			if (dest == nullptr) continue;
			dest->add_placement(i->get_bird_bcube(), !i->dim, i->dir, rgen.rand_bool(), obj_edge_spacing, rgen); // random orient_dir
		}
		for (auto i = trashcans.begin()+trashcans_start; i != trashcans.end(); ++i) {
			cube_t const bcube(i->get_bird_bcube());
			if (min(bcube.dx(), bcube.dy()) <= 2.0*obj_edge_spacing) continue;
			vect_bird_place_t *const dest(select_bird_loc_dest(add_pigeons, add_birds, pigeon_locs, bird_locs, rgen));
			if (dest == nullptr) continue;
			cube_t top_place(bcube);
			top_place.expand_by_xy(-1.5*obj_edge_spacing); // small shrink
			dest->add_placement_rand_dim_dir(top_place, obj_edge_spacing, rgen); // facing outward on a random side of the rim
		}
		add_objs_top_center(sstations, substations_start, add_pigeons, add_birds, pigeon_locs, bird_locs, rgen);
		add_objs_top_center(fountains, fountains_start,   add_pigeons, add_birds, pigeon_locs, bird_locs, rgen);

		// place pigeons
		if (add_pigeons) {
			// place some random pigeons; use place_radius because pigeon radius hasn't been calculated yet
			unsigned const count_mod(is_park ? 9 : 5), num_pigeons(rgen.rand() % count_mod); // 0-4, 0-8 for parks
			for (unsigned n = 0; n < num_pigeons; ++n) {pigeon_locs.emplace_back(rand_xy_pt_in_cube(plot, place_radius, rgen), rgen);}

			for (unsigned i = 0; i < pigeon_locs.size(); ++i) {
				bird_place_t p(pigeon_locs[i]);
				float const height(base_height*rgen.rand_uniform(0.8, 1.2));
				// the current model's tail extends below its feet; move down slightly so that feet are on the object, though the tail may clip through the object;
				// the feet gap isn't really visible when placed on the ground since there are no shadows, and it looks better than having the tail clip through the ground
				if (!p.on_ground) {p.pos.z -= 0.15*height;}
				pigeon_t const pigeon(p.pos, height, p.orient);
				if (p.on_ground && has_bcube_int_xy(pigeon.bcube, blockers, 2.0*pigeon.radius)) continue; // placed on the ground - check for collisions
				if (p.on_ground) {blockers.push_back(pigeon.bcube);} // not needed? don't need to add to pedestrian colliders
				pigeon_groups.add_obj(pigeon, pigeons);

				if (i == 0 && p.on_ground) { // place a group of pigeons on the ground for the first pigeon
					float const place_range(0.5*car_length);
					unsigned const group_size(rgen.rand() % count_mod); // 0-4 more, 0-8 for parks
					cube_t valid_range(plot);
					valid_range.expand_by_xy(-place_radius);

					for (unsigned N = 0; N < group_size; ++N) {
						bird_place_t p2(p);
						p2.pos += rgen.signed_rand_vector_spherical_xy(place_range);
						if (!valid_range.contains_pt_xy(p2.pos) || pigeon.bcube.contains_pt_xy_exp(p2.pos, base_height)) continue;
						p2.set_rand_orient(rgen);
						pigeon_locs.push_back(p2);
					}
				}
			} // for i
		}
	}
}

// dim=narrow dimension of fence; dir=house front/back dir in dimension dim; side=house left/right side in dimension !dim
float extend_fence_to_house(cube_t &fence, cube_t const &house, float fence_hwidth, float fence_height, bool dim, bool dir, bool side) {
	float &fence_end(fence.d[!dim][!side]);
	fence_end = house.d[!dim][side]; // adjacent to the house
	set_wall_width(fence, house.d[dim][dir], fence_hwidth, dim);
	// try to expand to the wall edge of two part houses by doing a line intersection query
	point p1, p2;
	p1.z     = p2.z    = fence.z1() + 0.25*fence_height; // slightly up from the bottom edge of the fence
	p1[ dim] = p2[dim] = fence.d[dim][!dir]; // use the side that overlaps the house bcube
	p1[!dim] = fence_end - (side ? -1.0 : 1.0)*fence_hwidth; // pull back slightly so that the start point isn't exactly at the house edge
	p2[!dim] = house.d[!dim][!side]; // end point is the opposite side of the house
	point p_int;
	if (!check_city_building_line_coll_bs(p1, p2, p_int)) return 0.0; // if this fails, house bcube must be wrong; should this be asserted?
	float const dist(fabs(fence_end - p_int[!dim]));
	fence_end = p_int[!dim];
	assert(fence.is_strictly_normalized());
	return dist;
}
bool check_valid_house_obj_place(point const &pos, float height, float radius, float wall_pos, bool dim, bool dir, cube_t const &bcube,
	cube_with_ix_t const &house, vect_cube_t const &blockers, unsigned prev_blockers_end, unsigned yard_blockers_start)
{
	point const center(pos + 0.5*height*plus_z);
	if (is_placement_blocked(bcube, blockers, house, prev_blockers_end))   return 0; // check blockers from prev step; no expand
	if (is_placement_blocked_recent(bcube, blockers, yard_blockers_start)) return 0; // check prev blockers in this yard
	if (check_sphere_coll_building(center, radius, 0, house.ix))           return 0; // xy_only=0
	point p_include(center);
	p_include[dim] = wall_pos - (dir ? 1.0 : -1.0)*0.25*radius; // slightly inside the wall
	if (!check_sphere_coll_building(p_include, 0.0, 0, house.ix)) return 0; // *must* be next to an exterior wall or fence; radius=0.0, xy_only=0
	if (check_close_to_door(pos, 4.0*radius, house.ix))           return 0;
	return 1;
}

void city_obj_placer_t::place_residential_plot_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<road_t> const &roads,
	vect_cube_t const &pool_blockers, unsigned driveways_start, unsigned city_ix, rand_gen_t &rgen)
{
	assert(plot_subdiv_sz > 0.0);
	sub_plots.clear();
	if (plot.is_park) return; // no dividers in parks
	subdivide_plot_for_residential(plot, roads, plot_subdiv_sz, 0, city_ix, sub_plots); // parent_plot_ix=0, not needed
	if (sub_plots.size() <= 1) return; // nothing to divide
	has_residential_plots = 1;
	if (rgen.rand_bool()) {std::reverse(sub_plots.begin(), sub_plots.end());} // reverse half the time so that we don't prefer a divider in one side or the other
	unsigned const shrink_dim(rgen.rand_bool()); // mostly arbitrary, could maybe even make this a constant 0
	float const sz_scale(0.06*city_params.road_width); // about 3 feet
	unsigned const dividers_start(dividers.size()), prev_blockers_end(blockers.size()); // prev_blockers_end is the end of blockers placed by previous steps for this plot
	vect_cube_with_ix_t bcubes; // we need the building index for the get_building_door_pos_closest_to() call

	for (auto i = sub_plots.begin(); i != sub_plots.end(); ++i) { // populate each yard
		unsigned const yard_blockers_start(blockers.size());
		// place plot dividers
		float hwidth(0.0), translate_dist[2] = {0.0, 0.0};
		unsigned const type(rgen.rand()%DIV_NUM_TYPES); // use a consistent divider type for all sides of this plot
		// chain link fence is not a primary divider; also, can't place a swimming pool here because it's not enclosed
		bool const add_divider(type != DIV_CHAINLINK);

		if (add_divider) {
			// should we remove or move house fences for divided sub-plots? I'm not sure how that would actually be possible at this point;
			// or maybe skip dividers if the house has a fence?
			plot_divider_type_t const &pdt(plot_divider_types[type]);
			float const z2(i->z1() + sz_scale*pdt.hscale);
			float const shrink_border(1.5*get_inner_sidewalk_width()); // needed for pedestrians to move along the edge of the plot; slightly larger to prevent collisions
			unsigned const prev_dividers_end(dividers.size());
			cube_t place_area(plot);
			place_area.expand_by_xy(-shrink_border);
			hwidth = 0.5*sz_scale*pdt.wscale;

			for (unsigned dim = 0; dim < 2; ++dim) {
				for (unsigned dir = 0; dir < 2; ++dir) {
					float const div_pos(i->d[dim][dir]);
					if (div_pos == plot.d[dim][dir]) continue; // sub-plot is against the plot border, don't need to add a divider
					bool const back_of_plot(i->d[!dim][0] != plot.d[!dim][0] && i->d[!dim][1] != plot.d[!dim][1]); // back of the plot, opposite the street
					unsigned const skip_dims(0); // can't make this (back_of_plot ? (1<<(1-dim)) : 0) because the edge may be showing at borders of different divider types
					cube_t c(*i);
					c.intersect_with_cube_xy(place_area);
					c.z2() = z2;
					set_wall_width(c, div_pos, hwidth, dim); // centered on the edge of the plot

					if (dim == shrink_dim) {
						translate_dist[dir] = (dir ? -1.0 : 1.0)*hwidth;
						c.translate_dim(dim, translate_dist[dir]); // move inside the plot so that edges line up
						// clip to the sides to remove overlap; may not line up with a neighboring divider of a different type/width, but hopefully okay
						for (unsigned d = 0; d < 2; ++d) {
							if (c.d[!dim][d] != plot.d[!dim][d]) {c.d[!dim][d] -= (d ? 1.0 : -1.0)*hwidth;}
						}
					}
					else {
						c.expand_in_dim(!dim, -0.001*hwidth); // fix for z-fighting
					}
					if (!back_of_plot) { // check for overlap of other plot dividers to the left and right
						cube_t test_cube(c);
						test_cube.expand_by_xy(4.0*hwidth); // expand so that adjacency counts as intersection
						bool overlaps(0);

						for (auto d = (dividers.begin()+dividers_start); d != (dividers.begin()+prev_dividers_end) && !overlaps; ++d) {
							overlaps |= (d->dim == bool(dim) && test_cube.contains_pt_xy(d->bcube.get_cube_center()));
						}
						if (overlaps) continue; // overlaps a previous divider, skip this one
					}
					divider_groups.add_obj(divider_t(c, type, dim, dir, skip_dims), dividers);
					add_cube_to_colliders_and_blockers(c, colliders, blockers);
				} // for dir
			} // for dim
		} // end dividers

		// place yard objects
		// Note: can't check for collisions with fire escapes and balcony support pillars here because they haven't been placed yet
		if (!i->is_residential || i->is_park || i->street_dir == 0) continue; // not a residential plot along a road
		bcubes.clear();
		if (have_buildings()) {get_building_bcubes(*i, bcubes);}
		if (bcubes.empty()) continue; // no house, skip adding other yard objects
		assert(bcubes.size() == 1); // there should be exactly one building/house in this sub-plot
		cube_with_ix_t const &house(bcubes.front());
		bool const sdim((i->street_dir-1)>>1), sdir((i->street_dir-1)&1); // direction to the road
		float const plot_z(i->z2());

		// attempt place swimming pool; often unsuccessful
		bool const placed_pool(add_divider && place_swimming_pool(plot, *i, house, sdim, sdir, shrink_dim, prev_blockers_end, hwidth,
			translate_dist, pool_blockers, blockers, colliders, rgen));
		bool const add_umbrella(building_obj_model_loader.is_model_valid(OBJ_MODEL_BIG_UMBRELLA) && rgen.rand_float() < 0.25); // 25% of the time
		bool placed_obj(placed_pool);

		if (!placed_obj || add_umbrella) { // place other objects
			cube_t place_area(*i);
			place_area.d[sdim][sdir] = house.d[sdim][!sdir]; // limit to the back yard
			place_area.expand_by_xy(-(0.1*city_params.road_width + hwidth)); // add some spacing, including divider width
			cube_t center_area(place_area);
			center_area.expand_by_xy(-0.2*city_params.road_width); // shrink the center placement area
			float const dx(center_area.dx()), dy(center_area.dy());

			if (dx > 0.0 && dy > 0.0) { // we have enough space
				if (!placed_obj && building_obj_model_loader.is_model_valid(OBJ_MODEL_SWINGSET) && rgen.rand_float() < 0.7) { // 70% of the time if ther was no pool
					float const ss_height(0.2*city_params.road_width);

					for (unsigned n = 0; n < 10; ++n) { // make some attempts to generate a valid swingset location
						point const ss_pos(rgen.gen_rand_cube_point_xy(center_area, plot_z));
						bool dim(0);
						if      (dx < 0.5*dy) {dim = 0;}
						else if (dy < 0.5*dx) {dim = 1;}
						else                  {dim = rgen.rand_bool();}
						swingset_t const ss(ss_pos, ss_height, dim, rgen.rand_bool()); // random dir, though dir doesn't really matter
						if (!place_area.contains_cube_xy(ss.bcube)) continue; // too close to back yard edge
						if (is_placement_blocked(ss.bcube, blockers, cube_t(), prev_blockers_end)) continue; // intersects some other object
						swing_groups.add_obj(ss, swings);
						add_cube_to_colliders_and_blockers(ss.bcube, colliders, blockers);
						placed_obj = 1;
						break; // success
					} // for n
				}
				if (!placed_obj && building_obj_model_loader.is_model_valid(OBJ_MODEL_TRAMPOLINE)) { // add a trampoline if there was no pool or swingset
					float const tramp_height(0.18*city_params.road_width);

					for (unsigned n = 0; n < 10; ++n) { // make some attempts to generate a valid trampoline location
						point const ss_pos(rgen.gen_rand_cube_point_xy(center_area, plot_z));
						trampoline_t const tramp(ss_pos, tramp_height, rgen);
						if (!place_area.contains_cube_xy(tramp.bcube)) continue; // too close to back yard edge
						if (is_placement_blocked(tramp.bcube, blockers, cube_t(), prev_blockers_end)) continue; // intersects some other object
						tramp_groups.add_obj(tramp, tramps);
						add_cube_to_colliders_and_blockers(tramp.bcube, colliders, blockers);
						break; // success
					} // for n
				}
				if (add_umbrella) { // maybe place a beach umbrella
					float const umbrella_height(0.25*city_params.road_width);

					for (unsigned n = 0; n < 10; ++n) { // make some attempts to generate a valid umbrella location
						point const ss_pos(rgen.gen_rand_cube_point_xy(center_area, plot_z));
						umbrella_t const umbrella(ss_pos, umbrella_height, 0, 0); // dim=0, dir=0
						if (!place_area.contains_cube_xy(umbrella.bcube)) continue; // too close to back yard edge
						//if (placed_pool && umbrella.bcube.intersects_xy(pools.back().bcube))             continue; // intersects the pool
						if (is_placement_blocked(umbrella.bcube, blockers, cube_t(), blockers.size())) continue; // intersects some other object; include all objects
						umbrella_groups.add_obj(umbrella, umbrellas);
						add_cube_to_colliders_and_blockers(umbrella.bcube, colliders, blockers);
						break; // success
					} // for n
				}
				// maybe place clothes line
				if (!placed_pool && rgen.rand_float() < 0.75) { // 75% of the time (but often fails)
					float const height(0.14*city_params.road_width), cl_z2(plot_z + height);

					for (unsigned n = 0; n < 40; ++n) { // make some attempts to generate a valid pair of points
						point const p1(rgen.gen_rand_cube_point_xy(center_area, cl_z2));
						bool const cdim(rgen.rand_bool()), cdir(rgen.rand_bool()); // should we always choose the furthest edge of center_area?
						point p2(p1);
						p2[cdim] += (cdir ? 1.0 : -1.0)*rgen.rand_uniform(0.4, 0.8)*city_params.road_width;
						if (!center_area.contains_pt_xy(p2)) continue;
						clothesline_t const cline(p1, p2, height, rgen);
						cube_t test_cube(cline.bcube);
						test_cube.expand_in_dim(!cdim, 0.75*height); // add extra padding to the sides - for example, for swing sets
						if (is_placement_blocked(test_cube, blockers, cube_t(), blockers.size())) continue; // intersects some other object; include all objects
						cline_groups.add_obj(cline, clines);
						blockers.push_back(cline.bcube);

						// add ped colliders for each pole, but not for the line itself
						for (unsigned d = 0; d < 2; ++d) {
							cube_t c(cline.ends[d]);
							c.z1() = plot_z;
							c.expand_by_xy(cline.pradius);
							colliders.push_back(c);
						}
						break; // done
					} // for n
				}
			} // end dx/dy check
		} // end place other objects
		
		// place short pine trees by the front
		if (!enable_instanced_pine_trees()) {
			unsigned const num_trees(rgen.rand() % 5); // 0-4
			float const pine_xy_sz = 0.5; // narrow

			for (unsigned n = 0; n < num_trees; ++n) {
				bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random house wall dim/dir
				float const wall_pos(house.d[dim][dir]);
				float const tree_scale(rgen.rand_uniform(0.25, 0.3)), radius(3.0*sz_scale*tree_scale), height(3.0*radius); // height is just a guess
				point pos;
				pos.z = plot_z;
				pos[ dim] = wall_pos + (dir ? 1.0 : -1.0)*1.5*radius; // place near this wall of the house
				pos[!dim] = rgen.rand_uniform(house.d[!dim][0], house.d[!dim][1]); // on the corner is okay
				cube_t const tree_bc(get_cube_height_radius(pos, radius, height));
				if (!check_valid_house_obj_place(pos, radius, radius, wall_pos, dim, dir, tree_bc, house, blockers, prev_blockers_end, yard_blockers_start)) continue;
				// Note: we can't test objects such as balconies and fire escapes, so we might end up with a pine tree intersecting them
				int const ttype(0); // 0=pine, 1=short pine, 2=palm
				place_tree(pos, radius, ttype, colliders, nullptr, 0, 0, 1, 0, tree_scale, pine_xy_sz); // tree_pos=nullptr, allow_bush=0, add_bush=0, is_sm_tree=1, has_planter=0
				blockers.push_back(tree_bc); // includes branches
			} // for n
		}
		// maybe place a trashcan next to the house
		if (1) {
			float const tc_height(0.18*city_params.get_nom_car_size().x), radius(0.3*tc_height);

			for (unsigned n = 0; n < 4; ++n) { // make some attempts to generate a valid trashcan location
				bool dim(0), dir(0);
				unsigned const house_side(rgen.rand() % 3); // any side but the front
				if (house_side == 0) {dim = sdim; dir = !sdir;} // put it in the back yard
				else {dim = !sdim; dir = rgen.rand_bool();} // put it on a random side of the house
				if (radius > 0.25*house.get_sz_dim(!dim)) continue; // house wall is too short - shouldn't happen in the normal case
				float const wall_pos(house.d[dim][dir]);
				point pos;
				pos.z = plot_z;
				pos[ dim] = wall_pos + (dir ? 1.0 : -1.0)*1.75*radius; // place near this wall of the house
				pos[!dim] = rgen.rand_uniform((house.d[!dim][0] + radius), (house.d[!dim][1] - radius));
				trashcan_t const trashcan(pos, radius, tc_height, 1); // is_cylin=1
				if (!check_valid_house_obj_place(pos, tc_height, radius, wall_pos, dim, dir, trashcan.bcube, house, blockers, prev_blockers_end, yard_blockers_start)) continue;
				trashcan_groups.add_obj(trashcan, trashcans);
				add_cube_to_colliders_and_blockers(trashcan.bcube, colliders, blockers);
				break; // success
			} // for n
		}
		// maybe place a bike next to the house wall or fence;
		// it would be nice if we could place a bike lying down, but the current model drawing code only supports rotations about the Z axis in the XY plane
		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_BICYCLE) && rgen.rand_float() < 0.6) { // 60% of the time (not always successful)
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_BICYCLE)); // L, W, H
			float const bike_height(1.3*sz_scale), bike_len(bike_height*sz.x/sz.z), bike_width(bike_height*sz.y/sz.z), wall_extend(0.5*bike_len);

			for (unsigned n = 0; n < 4; ++n) { // make some attempts to generate a valid bike location
				bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random house wall dim/dir
				if (bike_len > 0.5*house.get_sz_dim(!dim)) continue; // house wall is too short - shouldn't happen in the normal case
				float const wall_pos(house.d[dim][dir]);
				point pos;
				pos.z = plot_z;
				pos[ dim] = wall_pos + (dir ? 1.0 : -1.0)*0.6*bike_width; // place near this wall of the house (with a small gap for the window sill and FP error)
				pos[!dim] = rgen.rand_uniform((house.d[!dim][0] + wall_extend), (house.d[!dim][1] - wall_extend));
				bicycle_t const bike(pos, bike_height, !dim, rgen.rand_bool()); // random dir, against the wall
				if (is_placement_blocked(bike.bcube, blockers, house, prev_blockers_end))   continue; // check blockers from prev step; no expand
				if (is_placement_blocked_recent(bike.bcube, blockers, yard_blockers_start)) continue; // check prev blockers in this yard
				point const center(bike.bcube.get_cube_center()); // pos shifted up
				point p_include(center);
				p_include[dim] = wall_pos - (dir ? 1.0 : -1.0)*0.2*bike_width; // slightly inside the wall
				if (!check_sphere_coll_building(p_include, 0.0, 0, house.ix)) continue; // *must* be next to an exterior wall or fence; radius=0.0, xy_only=0
				// bikes aren't sphere-ish, so can't call check_valid_house_obj_place(); instead we check three points in the front, middle, and back
				point p_exclude[3] = {center, center, center}; // {middle, front, back}
				for (unsigned d = 0; d < 2; ++d) {p_exclude[d+1][!dim] += (d ? 1.0 : -1.0)*0.5*(bike_len - bike_width);} // front and back
				bool blocked(0);
				for (unsigned d = 0; d < 3; ++d) {blocked |= check_sphere_coll_building(p_exclude[d], 0.5*bike_width, 0, house.ix);}
				if (blocked) continue;
				if (check_close_to_door(pos, 1.0*bike_len, house.ix)) continue;
				bike_groups.add_obj(bike, bikes);
				add_cube_to_colliders_and_blockers(bike.bcube, colliders, blockers);
				break; // success
			} // for n
		}
		// maybe place potted plants next to the house
		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_PLANT)) {
			unsigned const num_plants(rgen.rand() % 9); // 0-8

			for (unsigned n = 0; n < num_plants; ++n) { // make one attempt per plant
				bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random house wall dim/dir
				float const wall_pos(house.d[dim][dir]), plant_height(sz_scale*rgen.rand_uniform(0.8, 1.25));
				point pos;
				pos.z = plot_z;
				pos[ dim] = wall_pos; // place at the house wall - will move away from the wall once we know the radius
				pos[!dim] = rgen.rand_uniform(house.d[!dim][0], house.d[!dim][1]); // on the corner is okay
				potted_plant_t plant(pos, plant_height, rgen.rand_bool(), rgen.rand_bool(), rgen.rand()); // random dim/dir/model
				float const radius(0.5*plant.bcube.get_sz_dim(dim)); // only care about size in the dim facing the wall
				plant.translate_dim(dim, (dir ? 1.0 : -1.0)*1.2*radius);
				if (!check_valid_house_obj_place(plant.pos, 0.0, radius, wall_pos, dim, dir, plant.bcube, house, blockers, prev_blockers_end, yard_blockers_start)) continue;
				plant_groups.add_obj(plant, plants);
				add_cube_to_colliders_and_blockers(plant.bcube, colliders, blockers);
			} // for n
		}
		// maybe place some flowers around the house
		// hack to enable custom alpha mipmaps for the flower model
		bool const orig_enable_model3d_custom_mipmaps(enable_model3d_custom_mipmaps);
		enable_model3d_custom_mipmaps = 1;
		bool const have_flowers(building_obj_model_loader.is_model_valid(OBJ_MODEL_FLOWER));
		enable_model3d_custom_mipmaps = orig_enable_model3d_custom_mipmaps;

		if (have_flowers) {
			unsigned const num_flower_groups(rgen.rand() % 4); // 0-3

			for (unsigned n = 0; n < num_flower_groups; ++n) { // make one attempt per flower
				bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random house wall dim/dir for all flowers in this group
				float const wall_pos(house.d[dim][dir]), flower_height(sz_scale*rgen.rand_uniform(0.75, 1.2));
				point pos;
				pos.z = plot_z;
				pos[ dim] = wall_pos; // place at the house wall - will move away from the wall once we know the radius
				pos[!dim] = rgen.rand_uniform(house.d[!dim][0], house.d[!dim][1]); // on the corner is okay
				unsigned const model_select(rgen.rand());
				flower_t flower(pos, flower_height, rgen.rand_bool(), rgen.rand_bool(), model_select); // random dim/dir
				float const radius(0.5*flower.bcube.get_sz_dim(dim)); // only care about size in the dim facing the wall
				float const row_spacing(1.25*flower.bcube.get_sz_dim(!dim));
				flower.translate_dim(dim, (dir ? 1.0 : -1.0)*1.25*radius);
				if (!check_valid_house_obj_place(flower.pos, 0.0, radius, wall_pos, dim, dir, flower.bcube, house, blockers, prev_blockers_end, yard_blockers_start)) continue;
				flower_groups.add_obj(flower, flowers);
				add_cube_to_colliders_and_blockers(flower.bcube, colliders, blockers);
				if (rgen.rand_bool()) continue; // single flower

				// create a row of flowers extending to the left and right of this one
				for (unsigned pdir = 0; pdir < 2; ++pdir) {
					for (unsigned m = 0; m < 4; ++m) {
						point pos2(flower.pos.x, flower.pos.y, pos.z);
						pos2[!dim] += (pdir ? 1.0 : -1.0)*(m+1)*row_spacing;
						flower_t const flower2(pos2, flower_height, rgen.rand_bool(), rgen.rand_bool(), model_select); // random dim/dir
						if (!check_valid_house_obj_place(flower2.pos, 0.0, radius, wall_pos, dim, dir, flower2.bcube, house, blockers, prev_blockers_end, yard_blockers_start)) break;
						flower_groups.add_obj(flower2, flowers);
						add_cube_to_colliders_and_blockers(flower2.bcube, colliders, blockers);
					}
				} // for pdir
			} // for n
		}
	} // for i (sub_plots)
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_MAILBOX)) {
		// place mailboxes on residential streets
		assert(driveways_start <= driveways.size());
		float const mbox_height(1.1*sz_scale);

		for (auto dw = driveways.begin()+driveways_start; dw != driveways.end(); ++dw) {
			if (rgen.rand_bool()) continue; // only 50% of houses have mailboxes along the road
			bool const dim(dw->dim), dir(dw->dir);
			unsigned pref_side(2); // starts invalid
			point pos(dw->get_cube_center()); // start at the driveway center

			for (auto i = sub_plots.begin(); i != sub_plots.end(); ++i) { // find subplot for this driveway
				if (!i->contains_pt_xy(pos)) continue; // wrong subplot
				pref_side = (pos[!dim] < i->get_center_dim(!dim)); // place mailbox on the side of the driveway closer to the center of the plot
				break; // done
			}
			if (pref_side == 2) continue; // no subplot found? error, or just skip the mailbox?
			// place at end of driveway at the road, but far enough back to leave space for peds
			pos[dim] = dw->d[dim][dir] - (dir ? 1.0 : -1.0)*(get_inner_sidewalk_width() + 0.5*mbox_height);

			for (unsigned n = 0; n < 2; ++n) {
				unsigned const side(pref_side ^ n);
				pos[!dim] = dw->d[!dim][side] + (side ? 1.0 : -1.0)*0.5*mbox_height; // off to the side of the driveway
				mailbox_t const mbox(pos, mbox_height, dim, dir);
				if (is_placement_blocked(mbox.bcube, blockers, *dw, prev_blockers_end, 0.0, 0)) continue; // try the other side
				mbox_groups.add_obj(mbox, mboxes);
				add_cube_to_colliders_and_blockers(mbox.bcube, colliders, blockers);
				break; // done
			} // for n
		} // for dw
	}
}

bool city_obj_placer_t::place_swimming_pool(road_plot_t const &plot, city_zone_t const &yard, cube_with_ix_t const &house,
	bool dim, bool dir, bool shrink_dim, unsigned prev_blockers_end, float divider_hwidth, float const translate_dist[2],
	vect_cube_t const &pool_blockers, vect_cube_t &blockers, vect_cube_t &colliders, rand_gen_t &rgen)
{
	if (rgen.rand_float() < 0.05) return 0; // add pools 95% of the time (since placing them often fails anyway)
	cube_t pool_area(yard);
	pool_area.d[dim][dir] = house.d[dim][!dir]; // limit the pool to the back yard

	for (unsigned d = 0; d < 2; ++d) {
		if (yard.d[!dim][d] == plot.d[!dim][d]) {pool_area.d[!dim][d] = house.d[!dim][d];} // adjacent to road - constrain to house projection so that side fence can be placed
	}
	float const dmin(min(pool_area.dx(), pool_area.dy())); // or should this be based on city_params.road_width?
	if (dmin < 0.75f*city_params.road_width) return 0; // back yard is too small to add a pool
	float const inground_pool_depth(0.17*city_params.road_width);
	cube_t inground_pool_area(pool_area);
	inground_pool_area.z1() -= inground_pool_depth;
	bool const has_blockers(has_bcube_int(inground_pool_area, pool_blockers)); // underground rooms below this yard
	bool const above_ground(rgen.rand_float() < (has_blockers ? 0.75 : 0.25)); // prefer above ground if there are underground rooms; otherwise, prefer in-ground
	float const min_pool_spacing_to_plot_edge(0.5*city_params.road_width), sz_scale(0.06*city_params.road_width); // about 3 feet

	for (unsigned d = 0; d < 2; ++d) { // keep pools away from the edges of plots; applies to sub-plots on the corners
		max_eq(pool_area.d[d][0], plot.d[d][0]+min_pool_spacing_to_plot_edge);
		min_eq(pool_area.d[d][1], plot.d[d][1]-min_pool_spacing_to_plot_edge);
	}
	pool_area.expand_by_xy(-(0.05*dmin + 0.05*city_params.road_width)); // small shrink to keep away from walls, fences, and hedges
	vector3d pool_sz;
	pool_sz.z = (above_ground ? rgen.rand_uniform(0.08, 0.12)*city_params.road_width : 0.01f*dmin);

	for (unsigned d = 0; d < 2; ++d) {
		pool_sz[d] = ((above_ground && d == 1) ? pool_sz[0] : rgen.rand_uniform(0.5, 0.7)*dmin); // above_ground_cylin pools have square bcubes
		pool_area.d[d][1] -= pool_sz[d]; // shrink so that pool_area is where (x1, x2) can be placed
	}
	if (!pool_area.is_normalized()) return 0; // pool area is too small; this can only happen due to shrink at plot edges
	colorRGBA const pool_side_colors[5] = {WHITE, WHITE, GRAY, LT_BROWN, LT_BLUE};

	for (unsigned n = 0; n < 20; ++n) { // make some attempts to generate a valid pool location
		point const pool_llc(rgen.gen_rand_cube_point_xy(pool_area, yard.z2()));
		cube_t pool(pool_llc, (pool_llc + pool_sz));
		cube_t tc(pool);
		tc.expand_by_xy(0.08*dmin);
		if (!above_ground) {tc.z1() -= inground_pool_depth;}
		if (is_placement_blocked(tc, blockers, cube_t(), prev_blockers_end))   continue; // intersects some other object
		if (!above_ground && has_blockers && has_bcube_int(tc, pool_blockers)) continue; // blocked by an extended basement room
		float const grayscale(rgen.rand_uniform(0.7, 1.0)), dirtyness(CLIP_TO_01(rgen.rand_uniform(-3.0, 1.2)));
		float const water_white_comp(rgen.rand_uniform(0.1, 0.3)), extra_green(rgen.rand_uniform(0.2, 0.5)), lightness(rgen.rand_uniform(0.5, 0.8));
		colorRGBA const color(above_ground ? pool_side_colors[rgen.rand()%5]: colorRGBA(grayscale, grayscale, grayscale));
		colorRGBA wcolor(lightness*water_white_comp, lightness*(water_white_comp + extra_green), lightness, 0.5); // A=0.5
		blend_color(wcolor, colorRGBA(0.1, 0.3, 0.15, 1.0), wcolor, dirtyness, 1); // blend in dark green dirty water; calc_alpha=1

		// add fences along the sides of the house to separate the back yard from the front yard; if fences can't be added, then don't add the pool either
		plot_divider_type_t const &fence_pdt(plot_divider_types[DIV_CHAINLINK]);
		float const fence_hwidth(0.5*sz_scale*fence_pdt.wscale), fence_height(sz_scale*fence_pdt.hscale), fence_z2(yard.z1() + fence_height);
		float const expand(-1.5*sz_scale*plot_divider_types[DIV_HEDGE].wscale); // shrink by widest divider to avoid false intersection with orthogonal dividers
		cube_t subplot_shrunk(yard);
		// translate so that fences line up with dividers; inexact if different width dividers are on each side
		for (unsigned d = 0; d < 2; ++d) {subplot_shrunk.d[shrink_dim][d] += translate_dist[d];}
		subplot_shrunk.expand_by_xy(-divider_hwidth); // shrink by half width of surrounding dividers
		divider_t fences[2];
		bool bad_fence_place(0);

		for (unsigned side = 0; side < 2; ++side) { // left/right
			bool fence_dim(dim);
			cube_t fence(subplot_shrunk);
			fence.z2() = fence_z2;

			if (yard.d[!dim][side] == plot.d[!dim][side]) { // at the edge of the plot, wrap the fence around in the back yard instead
				fence_dim ^= 1;
				extend_fence_to_house(fence, house, fence_hwidth, fence_height, !dim, side, !dir);

				if (is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {
					// Note: can't safely move the fence to the middle of the house if it intersects the pool because it may intersect or block a door
					if ((house.d[!dim][!side] > pool.d[!dim][side]) ^ side) {bad_fence_place = 1; break;} // fence at back of house does not contain the pool
					extend_fence_to_house(fence, house, fence_hwidth, fence_height, !dim, !side, !dir); // try the back edge of the house
					if (is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {bad_fence_place = 1; break;} // blocked by a driveway, etc.
				}
			}
			else { // at the front of the house
				float const ext_dist(extend_fence_to_house(fence, house, fence_hwidth, fence_height, dim, dir, side));

				// check if front fence position is bad, or fence extension is too long (may block off the porch)
				if (ext_dist > 0.33*house.get_sz_dim(!dim) || is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {
					extend_fence_to_house(fence, house, fence_hwidth, fence_height, dim, !dir, side); // try the back edge of the house
					if (is_placement_blocked(fence, blockers, house, prev_blockers_end, expand, !fence_dim)) {bad_fence_place = 1; break;} // blocked by a driveway, etc.
				}
			}
			fences[side] = divider_t(fence, DIV_CHAINLINK, fence_dim, dir, 0); // Note: dir is unused in divider_t so doesn't have to be set correctly
		} // for side
		if (bad_fence_place) continue; // failed to fence off the pool, don't place it here
		cube_t pool_full_height(pool);
		if (!above_ground) {pool_full_height.z1() -= inground_pool_depth;} // actual pool extends below the ground
		bool const sloped(!above_ground && rgen.rand_bool()); // 50% of in-ground pools have sloped bottoms
		swimming_pool_t const pool_obj(pool_full_height, color, wcolor, above_ground, sloped, dim, dir);
		pool_groups.add_obj(pool_obj, pools);
		unsigned const pre_pool_blockers_end(blockers.size());
		cube_t pool_collider(pool), ladder;
		pool_collider.z2() += 0.1*city_params.road_width; // extend upward to make a better collider
		if (above_ground) {pool_collider.d[dim][dir] += (dir ? 1.0 : -1.0)*pool_obj.get_ladder_depth();} // include external ladder
		add_cube_to_colliders_and_blockers(pool_collider, colliders, blockers);
		bool added_deck(0);

		for (unsigned side = 0; side < 2; ++side) {
			divider_t const &fence(fences[side]);
			assert(fence.bcube.is_strictly_normalized());
			divider_groups.add_obj(fence, dividers);
			add_cube_to_colliders_and_blockers(fence.bcube, colliders, blockers);
		}
		if (!above_ground && rgen.rand_float() < 0.8) { // in-ground pool, add pool deck 80% of the time
			bool had_cand(0);

			for (unsigned d = 0; d < 2 && !had_cand; ++d) { // check for projections in both dims
				float const deck_height(0.4*pool.dz()), deck_shrink(0.8*deck_height);
				float const plo(pool.d[d][0]), phi(pool.d[d][1]), hlo(house.d[d][0]), hhi(house.d[d][1]);
				if (max(plo, hlo) > min(phi, hhi)) continue; // no shared projection in this dim
				// randomly choose to extend deck to cover both objects or restrict to the shared length
				float const lo((rgen.rand_bool() ? min(plo, hlo) : max(plo, hlo)) + deck_shrink), hi((rgen.rand_bool() ? max(phi, hhi) : min(phi, hhi)) - deck_shrink);
				if (hi - lo < 0.25*pool.get_sz_dim(d)) continue; // not enough shared area
				bool const dir(house.get_center_dim(!d) < pool.get_center_dim(!d)); // dir from house to pool
				cube_t deck;
				set_cube_zvals(deck, pool.z1(), (pool.z1() + deck_height)); // lower height than the pool
				deck.d[ d][0   ] = lo;
				deck.d[ d][1   ] = hi;
				deck.d[!d][ dir] = pool .d[!d][!dir];
				deck.d[!d][!dir] = house.d[!d][ dir];
				assert(deck.is_strictly_normalized());
				had_cand = 1; // we have a candidate, even if it was blocked, so there are no more cands to check
				// check for partial intersections of objects such as trashcans, etc. and skip the deck in that case; allow contained objects to rest on the deck
				bool valid(1);

				for (auto b = blockers.begin()+prev_blockers_end; b != blockers.begin()+pre_pool_blockers_end; ++b) { // skip pool, pool fence, and house
					if (*b != house && deck.intersects_xy_no_adj(*b) && !deck.contains_cube_xy(*b)) {valid = 0; break;}
				}
				if (!valid) continue;
				float const roof_height(0.181*city_params.road_width);
				cube_t roof;
				
				if (deck.get_sz_dim(!d) > 0.8*roof_height) { // add if wide enough
					cube_t deck_with_roof(deck);
					deck_with_roof.z2() += roof_height; // extend upward
					cube_t const part_bounds(register_deck_and_get_part_bounds(house.ix, deck_with_roof));
					
					if (part_bounds.is_strictly_normalized()) { // valid house part edge; this only fails for one house
						roof = deck;
						roof.z2() = deck_with_roof.z2();
						roof.z1() = roof.z2() - 0.035*roof.dz();
						max_eq(roof.d[d][0], part_bounds.d[d][0]); // clip to range of house part(s)
						min_eq(roof.d[d][1], part_bounds.d[d][1]);
						assert(roof.is_strictly_normalized());
						if (roof.get_sz_dim(d) < 2.0*roof_height) {roof.set_to_zeros();} // disable if roof it too short
					}
				}
				pdeck_groups.add_obj(pool_deck_t(deck, roof, rgen.rand(), !d, dir), pdecks);// choose a random material
				blockers.push_back(pdecks.back().bcube); // blocker for other objects, but not a collider for people or the player
				added_deck = 1;
				break; // only needed in one dim
			} // for d
		}
		if (!above_ground && building_obj_model_loader.is_model_valid(OBJ_MODEL_POOL_LAD)) { // in-ground pool, add pool ladder
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_POOL_LAD)); // D, W, H
			float const height(1.6*sz_scale), hdepth(0.5*height*sz.x/sz.z), hwidth(0.5*height*sz.y/sz.z);

			if (hdepth < 0.2*min(pool.dx(), pool.dy())) { // add if pool is large enough; should be true
				ladder.z1() = pool.z2() - 0.425*height;
				ladder.z2() = ladder.z1() + height;
				set_wall_width(ladder, rgen.rand_uniform((pool.d[!dim][0] + 2.0*hdepth), (pool.d[!dim][1] - 2.0*hdepth)), hdepth, !dim); // pos along pool length
				set_wall_width(ladder, (pool.d[dim][dir] - (dir ? 1.0 : -1.0)*1.15*hwidth), hwidth, dim); // at pool edge
				plad_groups.add_obj(pool_ladder_t(ladder, dim, !dir), pladders); // inside the pool, facing the street/house - no placement check or blockers added
			}
		}
		if (added_deck) {
			pool_deck_t &pdeck(pdecks.back());
			pdeck.calc_pillars(ladder); // pass in the ladder pos to avoid placing pillars there
			vector_add_to(pdeck.pillars, colliders); // include pillars for AI collisions
			
			if (building_obj_model_loader.is_model_valid(OBJ_MODEL_DECK_CHAIR)) { // add deck chair
				vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_DECK_CHAIR)); // D, W, H
				float const height(1.2*sz_scale), hdepth(0.5*height*sz.x/sz.z), hwidth(0.5*height*sz.y/sz.z), radius(max(hwidth, hdepth));

				if (radius < 0.2*pdeck.get_width() && radius < 0.4*pdeck.get_depth()) { // add if deck is large enough; uses radius rather than hwidth/hdepth
					cube_t place_area(pdeck.base), chair;
					place_area.expand_in_dim( dim, -0.25*hdepth);
					place_area.expand_in_dim(!dim, -0.10*hwidth);
					point pos, door_pos;
					// use custom blockers that includes the ladder and pillars but not the deck
					vect_cube_t chair_blockers(blockers.begin(), blockers.begin()+pre_pool_blockers_end);
					chair_blockers.push_back(ladder);
					vector_add_to(pdeck.pillars, chair_blockers);
					
					if (get_building_door_pos_closest_to(house.ix, place_area.get_cube_center(), door_pos, 1)) { // check for house back door at pool deck
						cube_t door_blocker;
						door_blocker.set_from_sphere(door_pos, 1.5*radius);
						chair_blockers.push_back(door_blocker);
					}
					for (unsigned n = 0; n < 4; ++n) { // 4 outer layer placement iterations
						if (!try_place_obj(place_area, chair_blockers, rgen, radius, radius, 10, pos, 0)) break; // failed
						if (check_sphere_coll_building(pos, radius, 0, house.ix)) continue; // collided with part of the house (AC unit, etc.)
						pos.z = place_area.z2();
						chair.set_from_point(pos);
						chair.z2() += height;
						chair.expand_in_dim( dim, hdepth);
						chair.expand_in_dim(!dim, hwidth);
						chair_groups.add_obj(deck_chair_t(chair, dim, !dir), chairs);
						add_cube_to_colliders_and_blockers(chair, colliders, blockers);
						break; // success
					} // for n
				}
			}
		}
		// maybe place pool float in the pool
		cube_t float_bcube;
		point center;

		if (rgen.rand_float() < 0.6) {
			float const radius(0.05*city_params.road_width);
			colorRGBA const color(pfloat_colors[rgen.rand() % NUM_PFLOAT_COLORS]);
			
			if (rgen.rand_float() < 0.5) { // place inside the pool
				if (pool_obj.place_obj_on_water(center, radius, 0.1*radius, rgen)) {
					pool_float_t const pfloat(center, radius, color);

					if (!pfloat.bcube.intersects(ladder)) {
						pfloat_groups.add_obj(pfloat, pfloats);
						float_bcube = pfloat.bcube; // record so that ball placement avoids intersecting the pool float
					}
				}
			}
			else { // place next to the pool
				for (unsigned n = 0; n < 10; ++n) {
					pool_float_t const pfloat(rgen.gen_rand_cube_point_xy(pool_area, plot.z2()), radius, color);
					if (has_bcube_int(pfloat.bcube, blockers)) continue; // Note: may intersect the in-ground pool ladder
					pfloat_groups.add_obj(pfloat, pfloats);
					float_bcube = pfloat.bcube; // record so that ball placement avoids intersecting the pool float
					add_cube_to_colliders_and_blockers(pfloat.bcube, colliders, blockers);
					break;
				}
			}
		}
		// maybe place beach ball in the pool or next to it
		if (rgen.rand_float() < 0.7) {
			float const radius(ball_types[BALL_TYPE_BEACH].radius*0.18*city_params.road_width/(12*8));
			vector3d const orient(rgen.signed_rand_vector_norm()); // random orient

			if (rgen.rand_float() < 0.5) { // place inside the pool
				if (pool_obj.place_obj_on_water(center, radius, 0.1*radius, rgen)) {
					beach_ball_t const ball(center, radius, orient);
					if (!ball.bcube.intersects(ladder) && !ball.bcube.intersects(float_bcube)) {bball_groups.add_obj(ball, bballs);}
				}
			}
			else { // place next to the pool
				for (unsigned n = 0; n < 10; ++n) {
					beach_ball_t const ball(rgen.gen_rand_cube_point_xy(pool_area, plot.z2()), radius, orient);
					if (has_bcube_int(ball.bcube, blockers) || ball.bcube.intersects(float_bcube)) continue; // Note: may intersect the in-ground pool ladder
					bball_groups.add_obj(ball, bballs);
					add_cube_to_colliders_and_blockers(ball.bcube, colliders, blockers);
					break;
				}
			}
		}
		return 1; // success - done with pool
	} // for n
	return 0; // placement failed
}

void city_obj_placer_t::place_birds(cube_t const &city_bcube, rand_gen_t &rgen) {
	city_zval = city_bcube.z2();
	if (!are_birds_enabled()) return;
	bird_poop_manager.init(city_bcube);

	for (power_pole_t const &pp : ppoles) { // must check for walkway clearance
		if (check_bird_walkway_clearance(pp.bcube)) {add_bird_loc(pp, bird_locs, rgen);}
		// add locs to top power lines
		point top_wires[2][3][2];
		pp.get_top_wire_end_pts(top_wires);
		
		for (unsigned d = 0; d < 2; ++d) {
			unsigned const n(rgen.rand() % 3); // select a random wire
			point const &A(top_wires[d][n][0]), &B(top_wires[d][n][1]);
			if (A == B) continue; // no wire
			float const t(rgen.rand_uniform(0.2, 0.8));
			point const pos(t*B + (1.0 - t)*A);
			cube_t c; c.set_from_sphere(pos, pp.get_pole_radius());
			if (check_bird_walkway_clearance(c)) {bird_locs.emplace_back(pos, !d, rgen.rand_bool());}
		}
	} // for pp
	vect_bird_place_t unused;
	add_objs_top_center(newsracks, 0, 0, 1, unused, bird_locs, rgen);
	add_objs_top_center(mboxes,    0, 0, 1, unused, bird_locs, rgen);
	add_objs_top_center(stopsigns, 0, 0, 1, unused, bird_locs, rgen);
	add_objs_top_center(swings,    0, 0, 1, unused, bird_locs, rgen);
	add_objs_top_center(picnics,   0, 0, 1, unused, bird_locs, rgen);
	add_objs_top_center(clines,    0, 0, 1, unused, bird_locs, rgen);

	for (sign_t const &sign : signs) {
		if (sign.small) continue; // skip small signs above doors
		if (!sign.free_standing) continue; // only free standing signs, since signs on building roofs are ver far from the player and too close to buildings
		vect_bird_place_t *const dest(select_bird_loc_dest(0, 1, unused, bird_locs, rgen));
		if (dest != nullptr) {add_bird_loc(sign, *dest, rgen);}
	}
	// include houses/office buildings and streetlights?
	for (auto i = dividers.begin(); i != dividers.end(); ++i) {
		if (rgen.rand() & 3) continue; // only add one in 4, since there are so many
		bird_locs.add_placement_centerline(i->get_bird_bcube(), i->dim, rgen.rand_bool(), rgen); // place somewhere along the divider with random dir
	}
	// place initial birds; some bird_locs may have been added by place_residential_plot_objects(), which is called first
	unsigned const num_locs(bird_locs.size()), num_place(min(200U, num_locs/5U)), num_tries(2*num_place); // 2 tries on average per bird
	float const car_length(city_params.get_nom_car_size().x), bird_height(0.075f*car_length);
	unsigned num_added(0);

	for (unsigned n = 0; n < num_tries; ++n) {
		unsigned const loc_ix(rgen.rand()%num_locs);
		bird_place_t &p(bird_locs[loc_ix]);
		if (p.in_use) continue;
		// birds are taller than pigeons because their wingspan is included in their bcube, which means their height is shorter relative to their radius
		float const height(bird_height*rgen.rand_uniform(0.8, 1.2));
		point bird_pos(p.pos);
		bird_pos.z += BIRD_ZVAL_ADJ*height;
		colorRGBA color(WHITE*rgen.rand_uniform(0.1, 1.0)); // random grayscale color
		bird_groups.add_obj(city_bird_t(bird_pos, height, p.orient, color, loc_ix, rgen), birds);
		p.in_use = 1;
		++num_added;
		if (num_added == num_place) break; // done
	} // for n
}

void city_obj_placer_t::add_house_driveways(road_plot_t const &plot, vect_cube_t &temp_cubes, rand_gen_t &rgen, unsigned plot_ix) {
	cube_t plot_z(plot);
	plot_z.z1() = plot_z.z2() = plot.z2() + 0.0002*city_params.road_width; // shift slightly up to avoid Z-fighting
	temp_cubes.clear();
	add_house_driveways_for_plot(plot_z, temp_cubes);

	for (auto i = temp_cubes.begin(); i != temp_cubes.end(); ++i) {
		bool dim(0), dir(0);
		get_closest_dim_dir_xy(*i, plot, dim, dir);
		driveways.emplace_back(*i, dim, dir, plot_ix);
	}
}

void city_obj_placer_t::place_stopsigns_in_isec(road_isec_t &isec) {
	if (isec.has_stoplight) return; // can't have both a stoplight and a stopsign
	if (isec.num_conn == 2) return; // skip for 2-way intersections (bends)
	float const height(isec.get_stop_sign_height()), width(0.3*height);

	for (unsigned n = 0; n < 4; ++n) { // place stop signs on each connector
		if (!(isec.conn & (1 << n))) continue; // no road in this dir
		bool const dim((n>>1) != 0), dir((n&1) == 0); // Note: dir is inverted here to represent car dir
		stopsign_t const ssign(isec.get_stop_sign_pos(n), height, width, dim, !dir, isec.num_conn);
		stopsign_groups.add_obj(ssign, stopsigns);
		isec.make_stop_signs();
	} // for n
}

void city_obj_placer_t::place_objects_in_isec(road_isec_t &isec, bool is_residential, vector<point> const &hospital_signs, rand_gen_t &rgen) {
	if (/*!is_residential &&*/ isec.num_conn == 2) { // bend in road at city corner
		// place some traffic cones
		float const car_length(city_params.get_nom_car_size().x); // used as a size reference for other objects
		float const cone_radius(0.075*car_length), isec_radius(0.25*(isec.dx() + isec.dy())), place_radius(2.0*isec_radius);
		point center(cube_top_center(isec));
		// set center of curvature
		center.x += ((isec.conn & 1) ? -1.0 : 1.0)*isec_radius;
		center.y += ((isec.conn & 4) ? -1.0 : 1.0)*isec_radius;

		for (unsigned A = 0; A < 360; A += 9) {
			float const angle(TO_RADIANS*A), dx(sin(angle)), dy(cos(angle));
			// filter out positions blocked by roads connected to this intersection, including cones along axes
			if (dx <  0.2 && (isec.conn & 1)) continue; // -x/W
			if (dx > -0.2 && (isec.conn & 2)) continue; // +x/E
			if (dy <  0.2 && (isec.conn & 4)) continue; // -y/S
			if (dy > -0.2 && (isec.conn & 8)) continue; // +y/N
			point pos(center);
			pos.x += place_radius*dx;
			pos.y += place_radius*dy;
			for (unsigned d = 0; d < 2; ++d) {pos[d] += 0.5*cone_radius*rgen.signed_rand_float();} // add a small amount of random jitter to the position
			tcone_groups.add_obj(traffic_cone_t(pos, cone_radius), tcones);
			// Note: colliders not needed since people normally don't walk here (through zombies can chase the player here);
			// also, colliders don't work anyway because these cones are in the intersection, which counts as the road, not a plot
		} // for angle
	}
	if (!is_residential) { // maybe add hospital signs; only consider the hospital closest to the center of the intersection
		point const center(isec.get_cube_center());
		float const dmin(2.0*city_params.road_spacing);
		float dmin_sq(dmin*dmin);
		point closest_hospital;

		for (point const &pos : hospital_signs) {
			float const dxa(fabs(center.x - pos.x)), dya(fabs(center.y - pos.y));
			if (dxa > city_params.road_spacing && dya > city_params.road_spacing) continue; // not along either road
			float const dsq(dxa*dxa + dya*dya);
			if (dsq < dmin_sq) {closest_hospital = pos; dmin_sq = dsq;}
		}
		if (closest_hospital != all_zeros) {
			float const dx(center.x - closest_hospital.x), dy(center.y - closest_hospital.y), dxa(fabs(dx)), dya(fabs(dy));
			bool const dim(dxa < dya), dir(dim ? (dy > 0.0) : (dx > 0.0)); // turn direction is in the dim that's further
			// dim=1, dir=1: hospital to N; W right arrow, E left  arrow
			// dim=1, dir=0: hospital to S; W left  arrow, E right arrow
			// dim=0, dir=1: hospital to E; N right arrow, S left  arrow
			// dim=0, dir=0: hospital to W; N left  arrow, S right arrow
			isec.hospital_dir |= ((1<<unsigned(dir)) << (dim ? 0 : 6)) | ((1<<unsigned(!dir)) << (dim ? 2 : 4));
		}
	}
}

void add_to_plot_colliders(cube_t const &bcube, float bcube_expand, vector<road_plot_t> const &plots, vector<vect_cube_t> &plot_colliders) {
	cube_t bcube_ext(bcube);
	bcube_ext.expand_by_xy(bcube_expand);

	for (unsigned i = 0; i < plots.size(); ++i) { // linear iteration; stop signs are added to isecs, not plots; seems to be only 0.05ms per call
		if (plots[i].intersects_xy(bcube_ext)) {plot_colliders[i].push_back(bcube); return;}
	}
}
void city_obj_placer_t::add_ssign_and_slight_plot_colliders(vector<road_plot_t> const &plots, vector<road_isec_t> const isecs[3], vector<vect_cube_t> &plot_colliders) const {
	assert(plots.size() == plot_colliders.size());
	float const bcube_expand(2.0*get_sidewalk_width()); // include sidewalk stop signs in their associated plots
	for (stopsign_t const &s : stopsigns) {add_to_plot_colliders(s.bcube, bcube_expand, plots, plot_colliders);}

	for (unsigned n = 0; n < 3; ++n) { // {2-way, 3-way, 4-way}
		for (road_isec_t const &isec : isecs[n]) {
			for (unsigned i = 0; i < plots.size(); ++i) { // another linear iteration, for the same reason
				isec.add_stoplight_bcubes_in_region(get_plot_coll_region(plots[i]), plot_colliders[i]);
			}
		}
	}
}

template<typename T> void city_obj_placer_t::draw_objects(vector<T> const &objs, city_obj_groups_t const &groups, draw_state_t &dstate,
	float dist_scale, bool shadow_only, bool has_immediate_draw, bool draw_qbd_as_quads, float specular, float shininess)
{
	if (groups.empty() || !dstate.check_cube_visible(groups.get_bcube(), dist_scale)) return;
	T::pre_draw(dstate, shadow_only);
	unsigned start_ix(0);
	assert(city_draw_qbds_t::empty());

	for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
		if (!dstate.check_cube_visible(*g, dist_scale)) continue; // VFC/distance culling for group
		if (has_immediate_draw) {dstate.begin_tile(g->get_cube_center(), 1, 1);} // must setup shader and tile shadow map before drawing
		assert(start_ix <= g->ix && g->ix <= objs.size());

		for (unsigned i = start_ix; i < g->ix; ++i) {
			T const &obj(objs[i]);
			if (dstate.check_sphere_visible(obj.pos, obj.get_bsphere_radius(shadow_only))) {obj.draw(dstate, *this, dist_scale, shadow_only);}
		}
		if (!city_draw_qbds_t::empty() || !dstate.hedge_draw.empty()) { // we have something to draw
			if (!has_immediate_draw) {dstate.begin_tile(g->get_cube_center(), 1, 1);} // will_emit_now=1, ensure_active=1
			bool must_restore_state(!dstate.hedge_draw.empty());
			dstate.hedge_draw.draw_and_clear(dstate.s);

			if (has_untex_verts()) { // draw untextured verts before qbd so that textures such as text can alpha blend on top
				dstate.set_untextured_material();
				untex_qbd.draw_and_clear(); // matte

				if (!untex_spec_qbd.empty()) { // specular
					dstate.s.set_specular(specular, shininess); // shiny; values are per-object type
					untex_spec_qbd.draw_and_clear();
					dstate.s.clear_specular();
				}
				dstate.unset_untextured_material();
				must_restore_state = 1;
			}
			if (!shadow_only && must_restore_state) {T::pre_draw(dstate, shadow_only);} // re-setup for below draw call and/or next tile
			if (draw_qbd_as_quads) {qbd.draw_and_clear_quads();} else {qbd.draw_and_clear();} // draw this group with current smap

			if (!emissive_qbd.empty()) { // draw emissive materials
				dstate.s.add_uniform_float("emissive_scale", 1.0); // 100% emissive
				if (draw_qbd_as_quads) {emissive_qbd.draw_and_clear_quads();} else {emissive_qbd.draw_and_clear();}
				dstate.s.add_uniform_float("emissive_scale", 0.0); // reset
			}
		}
	} // for g
	T::post_draw(dstate, shadow_only);
}

void city_obj_placer_t::gen_parking_and_place_objects(vector<road_plot_t> &plots, vector<vect_cube_t> &plot_colliders, vector<car_t> &cars, vector<road_t> const &roads,
	vector<road_isec_t> isecs[3], cube_t const &city_bcube, vect_cube_t const &plot_cuts, unsigned city_id, bool have_cars, bool is_residential, bool have_streetlights)
{
	// Note: fills in plots.has_parking
	//highres_timer_t timer("Gen Parking Lots and Place Objects");
	vect_cube_t blockers, temp_cubes, underground_blockers; // blockers, driveways, extended basement rooms
	vector<point> tree_pos, hospital_signs;
	rand_gen_t rgen, detail_rgen;
	rgen.set_state(city_id, 123);
	detail_rgen.set_state(3145739*(city_id+1), 1572869*(city_id+1));
	if (city_params.max_trees_per_plot > 0) {tree_placer.begin_block(0); tree_placer.begin_block(1);} // both small and large trees
	bool const add_parking_lots(/*have_cars &&*/ !is_residential && city_params.min_park_spaces > 0 && city_params.min_park_rows > 0);
	float const sidewalk_width(get_sidewalk_width());
	get_building_ext_basement_bcubes(city_bcube, underground_blockers); // used for inground swimming pools and ponds in parks

	for (auto i = plots.begin(); i != plots.end(); ++i) { // calculate num_x_plots and num_y_plots; these are used for determining edge power poles
		max_eq(num_x_plots, i->xpos+1U);
		max_eq(num_y_plots, i->ypos+1U);
	}
	if (!is_residential) { // commercial city office buildings; add walkways and skyway
		vect_bldg_walkway_t walkway_cands;
		get_walkways_for_city(city_bcube, walkway_cands);
		// Note: not added to colliders since walkways are above pedestrians; not added to blockers since walkways are above most objects
		
		if (add_skyway(city_bcube, walkway_cands, rgen)) {
			walkway_cands.clear();
			get_walkways_for_city(city_bcube, walkway_cands); // query walkways again, this time including skyways
		}
		for (bldg_walkway_t const &w : walkway_cands) {walkway_groups.add_obj(walkway_t(w), walkways);}
	}
	for (auto i = plots.begin(); i != plots.end(); ++i) {
		tree_pos.clear();
		blockers.clear();
		get_building_bcubes(*i, blockers);
		size_t const plot_id(i - plots.begin()), buildings_end(blockers.size());
		assert(plot_id < plot_colliders.size());
		vect_cube_t &colliders(plot_colliders[plot_id]); // used for pedestrians
		if (add_parking_lots && !i->is_park) {i->has_parking = gen_parking_lots_for_plot(*i, cars, city_id, plot_id, blockers, colliders, plot_cuts, rgen, have_cars);}
		unsigned const driveways_start(driveways.size());
		if (is_residential) {add_house_driveways(*i, temp_cubes, detail_rgen, plot_id);}

		// driveways become blockers for other placed objects; make sure they extend into the road so that they intersect any placed streetlights or fire hydrants
		for (auto j = driveways.begin()+driveways_start; j != driveways.end(); ++j) {
			cube_t dw(*j);

			for (unsigned d = 0; d < 2; ++d) {
				if      (dw.d[d][0] == i->d[d][0]) {dw.d[d][0] -= sidewalk_width;}
				else if (dw.d[d][1] == i->d[d][1]) {dw.d[d][1] += sidewalk_width;}
			}
			blockers.push_back(dw);
		} // for j
		for (cube_t const &c : plot_cuts) { // add underground elevators and skylights to blockers since they have already been placed
			if (i->intersects_xy(c)) {blockers.push_back(c);}
		}
		if (city_params.assign_house_plots && plot_subdiv_sz > 0.0) {
			place_residential_plot_objects(*i, blockers, colliders, roads, underground_blockers, driveways_start, city_id, detail_rgen); // before placing trees
		}
		place_trees_in_plot  (*i, blockers, colliders, tree_pos, plot_cuts, detail_rgen, buildings_end);
		place_detail_objects (*i, blockers, colliders, tree_pos, underground_blockers, detail_rgen, have_streetlights);
		add_objs_on_buildings(*i, blockers, colliders, hospital_signs);
	} // for i (plot)
	for (unsigned n = 0; n < 3; ++n) {
		for (road_isec_t &isec : isecs[n]) {
			place_stopsigns_in_isec(isec); // Note: not a plot, can't use plot colliders
			place_objects_in_isec(isec, is_residential, hospital_signs, rgen);
		}
	}
	bind_elevators_to_building_walkways(city_bcube); // after all plots are populate and elevators are added
	add_ssign_and_slight_plot_colliders(plots, isecs, plot_colliders);
	connect_power_to_buildings(plots);
	if (have_cars && is_residential) {add_cars_to_driveways(cars, plots, plot_colliders, city_id, rgen);}
	place_birds(city_bcube, rgen); // after placing other objects
	bench_groups   .create_groups(benches,   all_objs_bcube);
	planter_groups .create_groups(planters,  all_objs_bcube);
	trashcan_groups.create_groups(trashcans, all_objs_bcube);
	fhydrant_groups.create_groups(fhydrants, all_objs_bcube);
	sstation_groups.create_groups(sstations, all_objs_bcube);
	fountain_groups.create_groups(fountains, all_objs_bcube);
	divider_groups .create_groups(dividers,  all_objs_bcube);
	pool_groups    .create_groups(pools,     all_objs_bcube);
	plad_groups    .create_groups(pladders,  all_objs_bcube);
	chair_groups   .create_groups(chairs,    all_objs_bcube);
	pdeck_groups   .create_groups(pdecks,    all_objs_bcube);
	ppole_groups   .create_groups(ppoles,    all_objs_bcube);
	hcap_groups    .create_groups(hcaps,     all_objs_bcube);
	manhole_groups .create_groups(manholes,  all_objs_bcube);
	mbox_groups    .create_groups(mboxes,    all_objs_bcube);
	tcone_groups   .create_groups(tcones,    all_objs_bcube);
	pigeon_groups  .create_groups(pigeons,   all_objs_bcube);
	bird_groups    .create_groups(birds,     all_objs_bcube);
	sign_groups    .create_groups(signs,     all_objs_bcube);
	stopsign_groups.create_groups(stopsigns, all_objs_bcube);
	flag_groups    .create_groups(flags,     all_objs_bcube);
	nrack_groups   .create_groups(newsracks, all_objs_bcube);
	cline_groups   .create_groups(clines,    all_objs_bcube);
	ppath_groups   .create_groups(ppaths,    all_objs_bcube);
	swing_groups   .create_groups(swings,    all_objs_bcube);
	tramp_groups   .create_groups(tramps,    all_objs_bcube);
	umbrella_groups.create_groups(umbrellas, all_objs_bcube);
	bike_groups    .create_groups(bikes,     all_objs_bcube);
	dumpster_groups.create_groups(dumpsters, all_objs_bcube);
	plant_groups   .create_groups(plants,    all_objs_bcube);
	flower_groups  .create_groups(flowers,   all_objs_bcube);
	picnic_groups  .create_groups(picnics,   all_objs_bcube);
	pond_groups    .create_groups(ponds,     all_objs_bcube);
	walkway_groups .create_groups(walkways,  all_objs_bcube);
	pillar_groups  .create_groups(pillars,   all_objs_bcube);
	wwe_groups     .create_groups(elevators, all_objs_bcube);
	uge_groups     .create_groups(ug_elevs,  all_objs_bcube);
	p_solar_groups .create_groups(p_solars,  all_objs_bcube);
	bball_groups   .create_groups(bballs,    all_objs_bcube);
	pfloat_groups  .create_groups(pfloats,   all_objs_bcube);
	if (skyway.valid) {all_objs_bcube.assign_or_union_with_cube(skyway.bcube);}
	if (add_parking_lots) {cout << "parking lots: " << parking_lots.size() << ", spaces: " << num_spaces << ", filled: " << filled_spaces << endl;}
}

void city_obj_placer_t::remap_parking_lot_ixs() {
	// sort of parking lots has invalidated the indices in parking spaces and driveways; generate the permutation vector and update the indices
	unsigned const num_pl(parking_lots.size());
	vector<unsigned> ix_map(num_pl, 0);
	
	for (unsigned i = 0; i < num_pl; ++i) {
		unsigned const orig_ix(parking_lots[i].orig_ix);
		assert(parking_lots[i].orig_ix < num_pl);
		assert(ix_map[orig_ix] == 0); // not yet assigned
		ix_map[orig_ix] = i;
	}
	for (driveway_t &d : driveways) {
		if (d.park_lot_ix >= 0) {assert(d.park_lot_ix < (int)num_pl); d.park_lot_ix = ix_map[d.park_lot_ix];}
	}
	for (parking_space_t &p : pspaces) {assert(p.p_lot_ix < num_pl); p.p_lot_ix = ix_map[p.p_lot_ix];}
}

// also adds signs on the ground near building doors
void city_obj_placer_t::add_objs_on_buildings(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> &hospital_signs) {
	// add signs and flags attached to buildings; note that some signs and flags may have already been added at this point
	vector<sign_t     > signs_to_add;
	vector<city_flag_t> flags_to_add;
	add_city_building_signs(plot, signs_to_add);
	add_city_building_flags(plot, flags_to_add);
	for (city_flag_t const &flag : flags_to_add) {flag_groups.add_obj(flag, flags);}
	cube_t plot_inner(plot);
	plot_inner.expand_by_xy(-get_sidewalk_width());

	for (sign_t const &sign : signs_to_add) {
		if (sign.sign_id >= 0 && !signs.empty() && sign.sign_id == signs.back().sign_id) continue; // already added a sign for this group

		if (sign.free_standing) { // sign on the ground, not on the building
			cube_t bcube_ext(sign.bcube);
			bcube_ext.expand_in_dim(sign.dim, 0.5*sign.bcube.get_sz_dim(!sign.dim)); // expand by half the sign width
			if (has_bcube_int(bcube_ext, blockers))       continue; // blocked, skip
			if (!plot_inner.contains_cube_xy(sign.bcube)) continue; // must stay inside the plot center
			add_cube_to_colliders_and_blockers(sign.bcube, colliders, blockers);
		}
		sign_groups.add_obj(sign, signs);
		if (sign.text == "Hospital" || sign.text == "Emergency") {hospital_signs.push_back(sign.pos);}
	} // for sign
}

/*static*/ bool city_obj_placer_t::subdivide_plot_for_residential(cube_t const &plot, vector<road_t> const &roads,
	float plot_subdiv_sz, unsigned parent_plot_ix, unsigned city_ix, vect_city_zone_t &sub_plots)
{
	if (min(plot.dx(), plot.dy()) < city_params.road_width) return 0; // plot is too small to divide
	assert(plot_subdiv_sz > 0.0);
	unsigned ndiv[2] = {0,0};
	float spacing[2] = {0,0};

	for (unsigned d = 0; d < 2; ++d) {
		float const plot_sz(plot.get_sz_dim(d));
		ndiv   [d] = max(1U, unsigned(round_fp(plot_sz/plot_subdiv_sz)));
		spacing[d] = plot_sz/ndiv[d];
	}
	if (ndiv[0] >= 100 || ndiv[1] >= 100) return 0; // too many plots? this shouldn't happen, but failing here is better than asserting or generating too many buildings
	unsigned const max_floors(0); // 0 is unlimited
	if (sub_plots.empty()) {sub_plots.reserve(2*(ndiv[0] + ndiv[1]) - 4);}

	for (unsigned y = 0; y < ndiv[1]; ++y) {
		float const y1(plot.y1() + spacing[1]*y), y2((y+1 == ndiv[1]) ? plot.y2() : (y1 + spacing[1])); // last sub-plot must end exactly at plot y2

		for (unsigned x = 0; x < ndiv[0]; ++x) {
			if (x > 0 && y > 0 && x+1 < ndiv[0] && y+1 < ndiv[1]) continue; // interior plot, no road access, skip
			float const x1(plot.x1() + spacing[0]*x), x2((x+1 == ndiv[0]) ? plot.x2() : (x1 + spacing[0])); // last sub-plot must end exactly at plot x2
			cube_t const c(x1, x2, y1, y2, plot.z1(), plot.z2());
			unsigned const street_dir(get_street_dir(c, plot)); // will favor x-dim for corner plots
			sub_plots.emplace_back(c, 0.0, 0, 1, street_dir, 1, parent_plot_ix, city_ix, max_floors); // cube, zval, park, res, sdir, capacity, ppix, cix, nf

			if (!roads.empty()) { // find the road name and determine the street address
				bool const sdim((street_dir - 1) >> 1), sdir((street_dir - 1) & 1);
				unsigned road_ix(0);
				float dmin(FLT_MAX), road_pos(0.0);

				for (auto r = roads.begin(); r != roads.end(); ++r) {
					if (r->dim == sdim) continue;
					float const dist(fabs(r->d[sdim][!sdir] - c.d[sdim][sdir]));
					if (dist < dmin) {dmin = dist; road_ix = (r - roads.begin()); road_pos = (c.get_center_dim(!sdim) - r->d[!sdim][0]);}
				}
				string const &road_name(roads[road_ix].get_name(city_ix));
				unsigned street_number(100 + 100*((7*road_ix + 13*road_name.size())%20)); // start at 100-2000 randomly based on road name and index
				street_number += 2*round_fp(25.0*road_pos/plot.get_sz_dim(!sdim)); // generate an even street number; should it differ by more than 1?
				if (sdir) {++street_number;} // make it an odd number if on this side of the road
				sub_plots.back().address    = std::to_string(street_number) + " " + road_name;
				sub_plots.back().street_num = street_number;
			}
		} // for x
	} // for y
	return 1;
}

bool city_obj_placer_t::connect_power_to_point(point const &at_pos, bool near_power_pole) {
	float dmax_sq(0.0);

	for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts to connect a to a power pole without intersecting a building
		float dmin_sq(0.0);
		unsigned best_pole(0);
		point best_pos;

		for (auto p = ppoles.begin(); p != ppoles.end(); ++p) {
			point const cur_pos(p->get_nearest_connection_point(at_pos, near_power_pole));
			if (cur_pos == at_pos) continue; // bad point
			float const dsq(p2p_dist_sq(at_pos, cur_pos));
			if (dsq <= dmax_sq) continue; // this pole was previously flagged as bad
			if (dmin_sq == 0.0 || dsq < dmin_sq) {best_pos = cur_pos; dmin_sq = dsq; best_pole = (p - ppoles.begin());}
		} // for p
		if (dmin_sq == 0.0) return 0; // failed (no power poles?)
		if (ppoles[best_pole].add_wire(at_pos, best_pos, near_power_pole)) return 1; // add a wire pole for houses
		dmax_sq = dmin_sq; // prevent this pole from being used in the next iteration
	} // for n
	return 0; // failed
}
void city_obj_placer_t::connect_power_to_buildings(vector<road_plot_t> const &plots) {
	if (plots.empty() || ppoles.empty() || !have_buildings()) return;
	cube_t all_plots_bcube(plots.front());
	for (auto p = plots.begin()+1; p != plots.end(); ++p) {all_plots_bcube.union_with_cube(*p);} // query all buildings in the entire city rather than per-plot
	vector<point> ppts;
	get_building_power_points(all_plots_bcube, ppts); // get points on house roofs
	for (point const &p : ppts) {connect_power_to_point(p, 1);} // near_power_pole=1
}

bool city_obj_placer_t::move_to_not_intersect_driveway(point &pos, float radius, bool dim) const {
	cube_t test_cube;
	test_cube.set_from_sphere(pos, radius);

	// Note: this could be accelerated by iterating by plot, but this seems to already be fast enough (< 1ms)
	for (auto d = driveways.begin(); d != driveways.end(); ++d) {
		if (!d->intersects_xy(test_cube)) continue;
		bool const dir((d->d[dim][1] - pos[dim]) < (pos[dim] - d->d[dim][0]));
		pos[dim] = d->d[dim][dir] + (dir ? 1.0 : -1.0)*0.1*city_params.road_width;
		return 1; // maybe we should check for an adjacent driveway, but that would be rare and moving could result in oscillation
	}
	return 0;
}
void city_obj_placer_t::finalize_streetlights_and_power(streetlights_t &sl, vector<vect_cube_t> &plot_colliders) {
	bool was_moved(0);

	for (auto s = sl.streetlights.begin(); s != sl.streetlights.end(); ++s) {
		if (!driveways.empty()) { // move to avoid driveways
			bool const dim(s->dir.y == 0.0); // direction to move in
			was_moved |= move_to_not_intersect_driveway(s->pos, 0.25*city_params.road_width, dim);
		}
		if (!ppoles.empty()) { // connect power
			point top(s->get_lpos());
			top.z += 1.05f*streetlight_ns::light_radius*city_params.road_width; // top of light
			connect_power_to_point(top, 0); // near_power_pole=0 because it may be too far away
		}
		if (s->on_bridge_or_tunnel) continue; // not inside a plot
		assert(s->plot_ix < plot_colliders.size());
		cube_t collider;
		collider.set_from_point(s->pos);
		collider.expand_by_xy(streetlight_ns::get_streetlight_pole_radius());
		collider.z2() += streetlight_ns::get_streetlight_height();
		plot_colliders[s->plot_ix].push_back(collider);
	} // for s
	if (was_moved) {sl.sort_streetlights_by_yx();} // must re-sort if a streetlight was moved
}

void city_obj_placer_t::add_manhole(point const &pos, float radius, bool is_over_road) {
	for (manhole_t const &m : manholes) { // first check to see if this is a duplicate
		if (m.pos.x == pos.x && m.pos.y == pos.y) return; // duplicate
	}
	point pos2(pos);
	if (!is_over_road && proc_sphere_coll(pos2, pos, zero_vector, radius, nullptr)) return; // check for blocker
	manhole_groups.add_obj(manhole_t(pos, radius), manholes);
	manhole_groups.rebuild(manholes, all_objs_bcube); // re-sort by tile
}

void city_obj_placer_t::add_city_ug_elevator_entrances(vect_ug_elev_info_t const &uges) {
	for (ug_elev_info_t const &e : uges) {uge_groups.add_obj(ug_elevator_t(e), ug_elevs);}
}

bool city_obj_placer_t::add_skyway(cube_t const &city_bcube, vect_bldg_walkway_t const &walkway_cands, rand_gen_t rgen) {
	if (!city_params.add_skyways) return 0;
	if (walkway_cands.empty())    return 0; // only add a skyway if this city has walkways
	bool const dim(city_bcube.dx() < city_bcube.dy()); // use longer dim
	unsigned const num_plots_wide(dim ? num_x_plots : num_y_plots); // num plots in !dim
	float centerline(city_bcube.get_center_dim(!dim)); // skyway is centered over the road by default
	float const road_width(city_params.road_width);

	if (num_plots_wide & 1) { // odd number of plots - skyway is off-center
		centerline += (rgen.rand_bool() ? 1.0 : -1.0)*0.5*(city_bcube.get_sz_dim(!dim) - road_width)/num_plots_wide; // shift half a plot to a random side
	}
	cube_t skyway_bc(city_bcube);
	skyway_bc.expand_in_dim(dim, -road_width); // shrink ends
	set_wall_width(skyway_bc, centerline, 0.4*road_width, !dim);
	cube_t clearance_area(skyway_bc);
	clearance_area.expand_by_xy(0.5*road_width);
	float ww_z1(city_bcube.z2() + 1.5*get_power_pole_height());
	
	for (bldg_walkway_t const &w : walkway_cands) {
		if (w.intersects_xy(clearance_area)) {max_eq(ww_z1, w.z2());} // make sure skyway is above all walkways
	}
	ww_z1 += 0.5*road_width;
	set_cube_zvals(skyway_bc, ww_z1, (ww_z1 + 0.4*road_width));
	if (!connect_buildings_to_skyway(skyway_bc, dim, city_bcube, skyway.ww_conns)) return 0;
	skyway.init(skyway_bc, dim); // only add skyway if it can connect to buildings
	vector<sign_t> skyway_signs;
	skyway.get_building_signs(skyway_signs);
	for (sign_t const &sign : skyway_signs) {sign_groups.add_obj(sign, signs);}
	return 1;
}

void city_obj_placer_t::bind_elevators_to_building_walkways(cube_t const &city_bcube) const {
	if (elevators.empty()) return;
	vector<building_walkway_t *> bwws;
	get_city_building_walkways(city_bcube, bwws);

	for (building_walkway_t *bww : bwws) { // walkways and elevators should be small enough that we can iterate over the cross product
		for (ww_elevator_t const &e : elevators) {
			if (bww->bcube.intersects(e.bcube)) {bww->attach_elevator(e.bcube); break;} // should be adjacent, and only one
		}
	}
}

void city_obj_placer_t::next_frame() {
	if (!animate2) return;
	float const fticks_stable(min(fticks, 1.0f)); // cap to 1/40s to improve stability
	point const camera_bs(get_camera_building_space());

	if (all_objs_bcube.contains_pt_xy(camera_bs)) { // player in city (approximate, since all_objs_bcube doesn't cover the entire city)
		for (swingset_t    &s : swings   ) {s.next_frame(camera_bs, fticks_stable);}
		for (ww_elevator_t &e : elevators) {e.next_frame(camera_bs, fticks_stable);}
	}
	next_frame_birds(camera_bs, fticks_stable);
}

void city_obj_placer_t::draw_detail_objects(draw_state_t &dstate, bool shadow_only) {
	float const dist_scale((player_in_basement >= 2) ? 0.1 : 1.0); // small distance scale for player in mall since only cur city is visible through skylight
	if (!dstate.check_cube_visible(all_objs_bcube, dist_scale)) return; // check bcube
	dstate.pass_ix = 0;
	draw_objects(benches,   bench_groups,    dstate, 0.16, shadow_only, 0); // dist_scale=0.16, has_immediate_draw=0
	draw_objects(fhydrants, fhydrant_groups, dstate, 0.06, shadow_only, 1);
	draw_objects(sstations, sstation_groups, dstate, 0.15, shadow_only, 1);
	draw_objects(fountains, fountain_groups, dstate, 0.20, shadow_only, 1);
	draw_objects(mboxes,    mbox_groups,     dstate, 0.04, shadow_only, 1);
	draw_objects(ppoles,    ppole_groups,    dstate, 0.20, shadow_only, 0);
	draw_objects(signs,     sign_groups,     dstate, 0.25, shadow_only, 1, 1); // draw_qbd_as_quads=1
	draw_objects(flags,     flag_groups,     dstate, 0.18, shadow_only, 1);
	draw_objects(newsracks, nrack_groups,    dstate, 0.10, shadow_only, 0);
	draw_objects(tcones,    tcone_groups,    dstate, 0.08, shadow_only, 1);
	draw_objects(swings,    swing_groups,    dstate, 0.06, shadow_only, 1);
	draw_objects(tramps,    tramp_groups,    dstate, 0.10, shadow_only, 1);
	draw_objects(umbrellas, umbrella_groups, dstate, 0.18, shadow_only, 1);
	draw_objects(bikes,     bike_groups,     dstate, 0.025,shadow_only, 1);
	draw_objects(dumpsters, dumpster_groups, dstate, 0.15, shadow_only, 1);
	draw_objects(plants,    plant_groups,    dstate, 0.04, shadow_only, 1);
	draw_objects(flowers,   flower_groups,   dstate, 0.06, shadow_only, 1);
	draw_objects(picnics,   picnic_groups,   dstate, 0.14, shadow_only, 1);
	draw_objects(chairs,    chair_groups,    dstate, 0.10, shadow_only, 1);
	draw_objects(walkways,  walkway_groups,  dstate, 0.25, shadow_only, 1);
	draw_objects(p_solars,  p_solar_groups,  dstate, 0.40, shadow_only, 0);
	draw_objects(elevators, wwe_groups,      dstate, 0.15, shadow_only, 0); // draw first pass opaque geometry
	draw_objects(ug_elevs,  uge_groups,      dstate, 0.20, shadow_only, 0);
	draw_objects(bballs,    bball_groups,    dstate, 0.12, shadow_only, 1);
	draw_objects(pfloats,   pfloat_groups,   dstate, 0.15, shadow_only, 1);
	
	if (!shadow_only) { // non shadow casting objects
		draw_objects(hcaps,    hcap_groups,    dstate, 0.12, shadow_only, 0);
		draw_objects(manholes, manhole_groups, dstate, 0.07, shadow_only, 1);
		draw_objects(pigeons,  pigeon_groups,  dstate, 0.03, shadow_only, 1);
		draw_objects(birds,    bird_groups,    dstate, 0.03, shadow_only, 1);
		draw_objects(pladders, plad_groups,    dstate, 0.06, shadow_only, 1);
		draw_objects(ppaths,   ppath_groups,   dstate, 0.25, shadow_only, 0, 1); // draw_qbd_as_quads=1

		for (dstate.pass_ix = 0; dstate.pass_ix < 2; ++dstate.pass_ix) { // {dirt bottom, dark blur}
			draw_objects(ponds, pond_groups, dstate, 0.30, shadow_only, 1); // dist_scale=0.30, has_immediate_draw=1
		}
	}
	for (dstate.pass_ix = 0; dstate.pass_ix < 3; ++dstate.pass_ix) { // {line, poles, clothes}
		if (shadow_only && dstate.pass_ix == 0) continue; // skip line in the shadow pass because its too narrow to cast a good shadow
		draw_objects(clines, cline_groups, dstate, 0.09, shadow_only, 1); // has_immediate_draw=1
	}
	dstate.s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA); // reset back to default after drawing 3D models such as fire hydrants and substations
	
	for (dstate.pass_ix = 0; dstate.pass_ix < 2; ++dstate.pass_ix) { // {concrete cube, metal cylinder}
		bool const is_cylin(dstate.pass_ix > 0);
		draw_objects(pillars, pillar_groups, dstate, 0.20, shadow_only, is_cylin); // dist_scale=0.25, has_immediate_draw=cylinder
	}
	for (dstate.pass_ix = 0; dstate.pass_ix < 2; ++dstate.pass_ix) { // {cube/city, cylinder/residential}
		bool const is_cylin(dstate.pass_ix > 0);
		draw_objects(trashcans, trashcan_groups, dstate, (is_cylin ? 0.08 : 0.10), shadow_only, is_cylin); // has_immediate_draw=cylinder
	}
	if (!shadow_only) { // low profile, not drawn in shadow pass
		for (dstate.pass_ix = 0; dstate.pass_ix < 2; ++dstate.pass_ix) { // {dirt, stone}
			draw_objects(planters, planter_groups, dstate, 0.1, shadow_only, 0); // dist_scale=0.1
		}
	}
	for (dstate.pass_ix = 0; dstate.pass_ix < 4; ++dstate.pass_ix) { // {0=in-ground walls, 1=in-ground water, 2=above ground sides, 3=above ground water}
		if (shadow_only && dstate.pass_ix != 2) continue; // only above ground pools are drawn in the shadow pass
		float const dist_scales[4] = {0.1, 0.5, 0.3, 0.5};
		draw_objects(pools, pool_groups, dstate, dist_scales[dstate.pass_ix], shadow_only, (dstate.pass_ix > 1)); // final 2 passes use immediate draw rather than qbd
	}
	if (!shadow_only && light_factor > 0.5 && !pools.empty()) { // 4=pool underwater caustics pass
		cube_t draw_region(all_objs_bcube);
		draw_region.expand_by(0.1*dstate.draw_tile_dist);

		if (draw_region.contains_pt(dstate.camera_bs)) { // camera in or near this city
			draw_state_t dstate2; // create a new draw state for a new shader
			dstate2.copy_from(dstate);
			dstate2.pass_ix = 4;
			shader_t &s(dstate2.s);
			// Note: this shader doesn't support fog, but caustics should only be drawn when close to the player
			s.setup_enabled_lights(1, (1 << SHADER_TYPE_FRAG)); // sun only, fragment shader
			s.set_prefix(make_shader_bool_prefix("use_shadow_map", 1), SHADER_TYPE_FRAG); // shadow maps always enabled, use fragment shader
			s.set_vert_shader("per_pixel_lighting");
			s.set_frag_shader("ads_lighting.part*+shadow_map.part*+caustics_overlay");
			s.begin_shader();
			select_texture(WATER_CAUSTIC_TEX);
			s.add_uniform_int  ("caustic_tex", 0);
			s.add_uniform_float("time", tfticks);
			s.add_uniform_float("shad_bias_scale", CITY_BIAS_SCALE);
			s.add_uniform_color("color_scale", sun_color);
			s.set_cur_color(WHITE); // not needed?
			setup_tile_shader_shadow_map(s);
			glPolygonOffset(-1.0, -1.0); // useful for avoiding z-fighting
			glEnable(GL_POLYGON_OFFSET_FILL);
			enable_blend();
			set_additive_blend_mode();
			set_std_depth_func_with_eq(); // <=
			glDepthMask(GL_FALSE); // disable depth writing
			draw_objects(pools, pool_groups, dstate2, 0.075, shadow_only, 1); // has_immediate_draw=1
			glDepthMask(GL_TRUE);
			set_std_depth_func(); // <
			set_std_blend_mode();
			disable_blend();
			glDisable(GL_POLYGON_OFFSET_FILL);
			if (dstate.s.is_setup()) {dstate.s.enable();} // enable(), not make_current(), because we may need to update the MVM for some reason
		}
	}
	// Note: not the most efficient solution, as it required processing blocks and binding shadow maps multiple times
	for (dstate.pass_ix = 0; dstate.pass_ix <= DIV_NUM_TYPES; ++dstate.pass_ix) { // {wall, fence, hedge, chainlink fence, chainlink fence posts}
		if (dstate.pass_ix == DIV_CHAINLINK && shadow_only) continue; // chainlink fence not drawn in the shadow pass
		draw_objects(dividers, divider_groups, dstate, 0.2, shadow_only, 0); // dist_scale=0.2
	}
	for (dstate.pass_ix = 0; dstate.pass_ix < NUM_POOL_DECK_PASSES; ++dstate.pass_ix) { // {wood, concrete, roof, pillars}
		if (shadow_only && dstate.pass_ix < NUM_POOL_DECK_TYPES) continue; // decks don't cast shadows, but the roof and pillars do
		draw_objects(pdecks, pdeck_groups, dstate, 0.26, shadow_only, 0); // dist_scale=0.3
	}
	for (dstate.pass_ix = 0; dstate.pass_ix < 3; ++dstate.pass_ix) { // {sign front, sign back + pole, 4-way sign}
		draw_objects(stopsigns, stopsign_groups, dstate, 0.1, shadow_only, 0); // dist_scale=0.1
	}
	dstate.pass_ix = 0; // reset back to 0
	if (!shadow_only) {bird_poop_manager.draw(dstate.s, dstate.xlate);}
	skyway.draw(dstate, *this, shadow_only); // must be last due to transparent roof
}
void city_obj_placer_t::draw_transparent_objects(draw_state_t &dstate) {
	if (!dstate.check_cube_visible(all_objs_bcube, 1.0)) return; // check bcube, dist_scale=1.0
	dstate.pass_ix = 2; // water surface
	draw_objects(ponds, pond_groups, dstate, 0.30, 0, 1); // dist_scale=0.30, shadow_only=0, has_immediate_draw=1
	dstate.pass_ix = 1; // transparent glass surfaces
	draw_objects(elevators, wwe_groups, dstate, 0.20, 0, 0); // dist_scale=0.20, shadow_only=0, has_immediate_draw=0
	dstate.pass_ix = 0; // reset back to 0
	skyway.draw_glass_surfaces(dstate, *this);
}

void city_obj_placer_t::add_lights(vector3d const &xlate, cube_t &lights_bcube) const {
	skyway.add_lights(xlate, lights_bcube);
}

template<typename T> bool proc_vector_sphere_coll(vector<T> const &objs, city_obj_groups_t const &groups, point &pos,
	point const &p_last, float radius, vector3d const &xlate, vector3d *cnorm)
{
	if (groups.empty() || !sphere_cube_intersect((pos - xlate), radius, groups.get_bcube())) return 0;
	point const pos_bs(pos - xlate);
	unsigned start_ix(0);

	for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
		if (!sphere_cube_intersect(pos_bs, radius, *g)) continue;
		assert(start_ix <= g->ix && g->ix <= objs.size());

		for (auto i = objs.begin()+start_ix; i != objs.begin()+g->ix; ++i) {
			if (i->proc_sphere_coll(pos, p_last, radius, xlate, cnorm)) return 1;
		}
	} // for g
	return 0;
}
bool city_obj_placer_t::proc_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, vector3d *cnorm) const { // pos in in camera space
	if (!sphere_cube_intersect(pos, (radius + p2p_dist(pos, p_last)), (all_objs_bcube + xlate))) return 0;
	bool const skyway_coll(skyway.proc_sphere_coll(pos, p_last, radius, xlate, cnorm)); // must be before walkways
	player_in_skyway |= skyway_coll;

	// special handling for player walking on/in walkways; we need to handle collisions with the top surface when above, so proc_vector_sphere_coll() can't be used here
	for (walkway_t const &w : walkways) {
		if (w.proc_sphere_coll(pos, p_last, radius, xlate, cnorm)) return 1;
	}
	if (skyway_coll) return 1;
	if (proc_vector_sphere_coll(benches,   bench_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(trashcans, trashcan_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(fhydrants, fhydrant_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(sstations, sstation_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(fountains, fountain_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(dividers,  divider_groups,  pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(pools,     pool_groups,     pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(ppoles,    ppole_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(mboxes,    mbox_groups,     pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(signs,     sign_groups,     pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(stopsigns, stopsign_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(flags,     flag_groups,     pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(newsracks, nrack_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(swings,    swing_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(tramps,    tramp_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(umbrellas, umbrella_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(bikes,     bike_groups,     pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(dumpsters, dumpster_groups, pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(plants,    plant_groups,    pos, p_last, radius, xlate, cnorm)) return 1; // optional?
	if (proc_vector_sphere_coll(chairs,    chair_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(picnics,   picnic_groups,   pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(ponds,     pond_groups,     pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(pillars,   pillar_groups,   pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(elevators, wwe_groups,      pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(ug_elevs,  uge_groups,      pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(pdecks,    pdeck_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(p_solars,  p_solar_groups,  pos, p_last, radius, xlate, cnorm)) return 1;
	if (proc_vector_sphere_coll(clines,    cline_groups,    pos, p_last, radius, xlate, cnorm)) return 1;
	// Note: no coll with tree_planters because the tree coll should take care of it;
	// no coll with hcaps, manholes, tcones, flowers, pladders, bballs, pfloats, pigeons, ppaths, or birds
	return 0;
}

template<typename T> void check_vector_line_intersect(vector<T> const &objs, city_obj_groups_t const &groups, point const &p1, point const &p2, float &t, bool &ret) {
	if (groups.empty() || !check_line_clip(p1, p2, groups.get_bcube().d)) return;
	unsigned start_ix(0);

	for (auto g = groups.begin(); g != groups.end(); start_ix = g->ix, ++g) {
		if (!check_line_clip(p1, p2, g->d)) continue;
		assert(start_ix <= g->ix && g->ix <= objs.size());
		for (auto i = objs.begin()+start_ix; i != objs.begin()+g->ix; ++i) {ret |= check_line_clip_update_t(p1, p2, t, i->bcube);}
	}
}
bool city_obj_placer_t::line_intersect(point const &p1, point const &p2, float &t) const { // p1/p2 in world space; uses object bcubes
	if (!all_objs_bcube.line_intersects(p1, p2)) return 0;
	bool ret(0);
	check_vector_line_intersect(benches,   bench_groups,    p1, p2, t, ret);
	check_vector_line_intersect(trashcans, trashcan_groups, p1, p2, t, ret);
	check_vector_line_intersect(fhydrants, fhydrant_groups, p1, p2, t, ret); // check bounding cube; cylinder intersection may be more accurate, but likely doesn't matter much
	check_vector_line_intersect(sstations, sstation_groups, p1, p2, t, ret);
	check_vector_line_intersect(fountains, fountain_groups, p1, p2, t, ret);
	check_vector_line_intersect(dividers,  divider_groups,  p1, p2, t, ret);
	check_vector_line_intersect(pools,     pool_groups,     p1, p2, t, ret);
	check_vector_line_intersect(ppoles,    ppole_groups,    p1, p2, t, ret); // inaccurate, could be customized if needed
	check_vector_line_intersect(signs,     sign_groups,     p1, p2, t, ret);
	check_vector_line_intersect(stopsigns, stopsign_groups, p1, p2, t, ret);
	check_vector_line_intersect(flags,     flag_groups,     p1, p2, t, ret);
	check_vector_line_intersect(newsracks, nrack_groups,    p1, p2, t, ret);
	check_vector_line_intersect(walkways,  walkway_groups,  p1, p2, t, ret);
	check_vector_line_intersect(pillars,   pillar_groups,   p1, p2, t, ret);
	check_vector_line_intersect(elevators, wwe_groups,      p1, p2, t, ret);
	check_vector_line_intersect(ug_elevs,  uge_groups,      p1, p2, t, ret);
	check_vector_line_intersect(dumpsters, dumpster_groups, p1, p2, t, ret);
	check_vector_line_intersect(picnics,   picnic_groups,   p1, p2, t, ret);
	// Note: nothing to do for parking lots, tree_planters, hcaps, manholes, tcones, flowers, pladders, chairs, pdecks, bballs, pfloats, clines, pigeons, ppaths, or birds;
	// mboxes, swings, tramps, umbrellas, bikes, plants, ponds, p_solars, and momorail are ignored because they're small or not simple shapes
	return ret;
}

bool city_obj_placer_t::intersects_parking_lot(cube_t const &c) const { // Note: currently called before parking lots and driveways are added
	return (has_bcube_int_xy(c, parking_lots) || has_bcube_int_xy(c, driveways)); // no acceleration structure for these, so do a linear iteration
}

template<typename T> bool check_city_obj_pt_xy_contains(city_obj_groups_t const &groups, vector<T> const &objs, point const &pos, unsigned &obj_ix, bool is_cylin=0) {
	if (groups.empty() || !groups.get_bcube().contains_pt_xy(pos)) return 0;
	unsigned start_ix(0);

	for (auto i = groups.begin(); i != groups.end(); start_ix = i->ix, ++i) {
		if (!i->contains_pt_xy(pos)) continue;
		assert(start_ix <= i->ix && i->ix <= objs.size());

		for (auto b = objs.begin()+start_ix; b != objs.begin()+i->ix; ++b) {
			if (pos.x < b->bcube.x1()) break; // objects are sorted by x1, none after this can match
			if (!b->check_point_contains_xy(pos)) continue;
			if (is_cylin && !dist_xy_less_than(pos, b->pos, b->get_overlay_radius())) continue; // cylinder case
			obj_ix = (b - objs.begin());
			return 1;
		}
	} // for i
	return 0;
}
bool city_obj_placer_t::get_color_at_xy_pre_road(point const &pos, colorRGBA &color) const { // check walkways because they can be over both roads and plots
	if (skyway.valid && skyway.bcube.contains_pt_xy(pos)) {color = WHITE; return 1;}
	unsigned obj_ix(0);
	if (check_city_obj_pt_xy_contains(walkway_groups, walkways, pos, obj_ix, 0)) {color = walkways[obj_ix].map_mode_color; return 1;} // is_cylin=0
	if (check_city_obj_pt_xy_contains(p_solar_groups, p_solars, pos, obj_ix, 0)) {color = colorRGBA(0.4, 0.4, 1.0); return 1;} // placed over parking lots; light blue
	return 0;
}
bool city_obj_placer_t::get_color_at_xy(point const &pos, vect_cube_t const &plot_cuts, colorRGBA &color, bool skip_in_road) const {
	unsigned obj_ix(0);
	if (check_city_obj_pt_xy_contains(bench_groups, benches, pos, obj_ix, 0)) {color = texture_color(FENCE_TEX); return 1;} // is_cylin=0
	float const expand(0.15*city_params.road_width), x_test(pos.x + expand); // expand to approx tree diameter

	if (!planter_groups.empty() && planter_groups.get_bcube().contains_pt_xy(pos)) {
		unsigned start_ix(0);

		for (auto i = planter_groups.begin(); i != planter_groups.end(); start_ix = i->ix, ++i) {
			if (!i->contains_pt_xy_exp(pos, expand)) continue;
			assert(start_ix <= i->ix && i->ix <= planters.size());

			for (auto p = planters.begin()+start_ix; p != planters.begin()+i->ix; ++p) {
				if (x_test < p->bcube.x1()) break; // planters are sorted by x1, none after this can match
				if (!p->bcube.contains_pt_xy_exp(pos, expand)) continue;
				// treat this as a tree rather than a planter by testing against a circle, since trees aren't otherwise included
				if (dist_xy_less_than(pos, p->pos, (p->radius + expand))) {color = DK_GREEN; return 1;}
			}
		} // for i
	}
	// fire hydrants are now placed on the edges of the road, so they're not inside plots and are skipped here
	if (!skip_in_road && check_city_obj_pt_xy_contains(fhydrant_groups, fhydrants, pos, obj_ix, 1)) {color = colorRGBA(1.0, 0.75, 0.0); return 1;} // orange/yellow; is_cylin=1
	
	if (check_city_obj_pt_xy_contains(divider_groups, dividers, pos, obj_ix)) {
		assert(obj_ix < dividers.size());
		color = plot_divider_types[dividers[obj_ix].type].get_avg_color();
		return 1;
	}
	if (!pool_groups.empty() && pool_groups.get_bcube().contains_pt_xy(pos)) {
		unsigned start_ix(0);

		for (auto i = pool_groups.begin(); i != pool_groups.end(); start_ix = i->ix, ++i) {
			if (!i->contains_pt_xy(pos)) continue;
			assert(start_ix <= i->ix && i->ix <= pools.size());

			for (auto b = pools.begin()+start_ix; b != pools.begin()+i->ix; ++b) {
				if (pos.x < b->bcube.x1()) break; // pools are sorted by x1, none after this can match
				if (!b->bcube.contains_pt_xy(pos)) continue;
				if (b->above_ground && !dist_xy_less_than(pos, cube_bot_center(b->bcube), b->get_radius())) continue; // circular in-ground pool
				color = colorRGBA(b->wcolor, 1.0); // return water color with alpha=1.0
				return 1;
			}
		} // for i
	}
	if (check_city_obj_pt_xy_contains(pdeck_groups, pdecks, pos, obj_ix)) {
		assert(obj_ix < pdecks.size());
		if (pdecks[obj_ix].has_roof()) {color = pool_deck_mats[NUM_POOL_DECK_TYPES  ].get_avg_color();} // roof visible
		else                           {color = pool_deck_mats[pdecks[obj_ix].mat_id].get_avg_color();} // deck visible
		return 1;
	}
	if (!skip_in_road && check_city_obj_pt_xy_contains(nrack_groups, newsracks, pos, obj_ix)) { // now placed in roads
		assert(obj_ix < newsracks.size());
		color = newsracks[obj_ix].color;
		return 1;
	}
	if (!pond_groups.empty() && pond_groups.get_bcube().contains_pt_xy(pos)) {
		for (pond_t const &pond : ponds) { // sparse enough to iterate over
			float const xv((pos.x - pond.pos.x)/(0.5*pond.bcube.dx())), yv((pos.y - pond.pos.y)/(0.5*pond.bcube.dy()));
			if ((xv*xv + yv*yv) < 1.0) {color = BLUE; return 1;} // ellipse collision
		}
	}
	if (check_city_obj_pt_xy_contains(sstation_groups, sstations, pos, obj_ix, 0)) {color = colorRGBA(0.6, 0.8, 0.4, 1.0); return 1;} // light olive
	if (check_city_obj_pt_xy_contains(trashcan_groups, trashcans, pos, obj_ix, 0)) {color = colorRGBA(0.8, 0.6, 0.3, 1.0); return 1;} // tan
	if (check_city_obj_pt_xy_contains(ppath_groups,    ppaths,    pos, obj_ix, 0)) {color = GRAY ; return 1;} // can/should we restrict this to only run when inside a park?
	if (check_city_obj_pt_xy_contains(fountain_groups, fountains, pos, obj_ix, 1)) {color = GRAY ; return 1;} // is_cylin=1
	if (check_city_obj_pt_xy_contains(tramp_groups,    tramps,    pos, obj_ix, 1)) {color = (BKGRAY*0.75 + tramps[obj_ix].color*0.25); return 1;} // is_cylin=1
	if (check_city_obj_pt_xy_contains(dumpster_groups, dumpsters, pos, obj_ix, 0)) {color = colorRGBA(0.1, 0.4, 0.1, 1.0); return 1;} // dark green
	if (check_city_obj_pt_xy_contains(umbrella_groups, umbrellas, pos, obj_ix, 1)) {color = WHITE; return 1;} // is_cylin=1
	if (check_city_obj_pt_xy_contains(picnic_groups,   picnics,   pos, obj_ix, 0)) {color = BROWN; return 1;}
	if (check_city_obj_pt_xy_contains(wwe_groups,      elevators, pos, obj_ix, 0)) {color = colorRGBA(0.8, 1.0, 0.8, 1.0); return 1;} // slightly blue-green glass; transparent?
	if (check_city_obj_pt_xy_contains(uge_groups,      ug_elevs,  pos, obj_ix, 0)) {color = LT_GRAY; return 1;}
	if (check_vect_cube_contains_pt_xy(plot_cuts, pos)) {color = colorRGBA(0.7, 0.7, 1.0); return 1;} // mall skylight; very light blue
	// Note: ppoles, hcaps, manholes, mboxes, tcones, flowers, pladders, chairs, stopsigns, flags, clines, pigeons, birds, swings, umbrellas, bikes, and plants are skipped;
	// pillars aren't visible under walkways;
	// free standing signs can be added, but they're small and expensive to iterate over and won't contribute much
	return 0;
}

void city_obj_placer_t::get_occluders(pos_dir_up const &pdu, vector3d const &xlate, vect_cube_t &occluders) const {
	if (player_in_skyway && skyway.valid) {
		if ((skyway.bcube + xlate).contains_pt(pdu.pos)) {occluders.push_back(skyway.get_floor_occluder() + xlate);} // should always be valid?
	}
	if (player_in_walkway) {
		for (walkway_t const &w : walkways) {
			if ((w.bcube + xlate).contains_pt(pdu.pos)) {occluders.push_back(w.get_floor_occluder() + xlate); break;} // can be only one
		}
	}
	//if (player_in_elevator) {}
	if (dividers.empty()) return; // dividers are currently the only other occluders
	float const dmax(0.25f*(X_SCENE_SIZE + Y_SCENE_SIZE)); // set far clipping plane to 1/4 a tile (currently 2.0)
	unsigned start_ix(0);

	for (auto i = divider_groups.begin(); i != divider_groups.end(); start_ix = i->ix, ++i) { // no divider_groups.get_bcube() test here?
		cube_t const gbc(*i + xlate);
		if (!dist_less_than(pdu.pos, gbc.closest_pt(pdu.pos), dmax) || !pdu.cube_visible(gbc)) continue;
		assert(start_ix <= i->ix && i->ix <= dividers.size());

		for (auto d = dividers.begin()+start_ix; d != dividers.begin()+i->ix; ++d) {
			assert(d->type < DIV_NUM_TYPES);
			if (!plot_divider_types[d->type].is_occluder) continue; // skip
			cube_t const bc(d->bcube + xlate);
			if (bc.z1() > pdu.pos.z || bc.z2() < pdu.pos.z) continue; // z-range does not include the camera
			if (dist_less_than(pdu.pos, bc.closest_pt(pdu.pos), dmax) && pdu.cube_visible(bc)) {occluders.push_back(bc);}
		}
	} // for i
}

void city_obj_placer_t::get_plot_cuts(cube_t const &plot, vect_cube_t &cuts) const {
	for (swimming_pool_t const &p : pools) {
		if (!p.above_ground && plot.intersects_xy(p.bcube)) {cuts.push_back(p.bcube);} // zvals are ignored
	}
}
bool city_obj_placer_t::cube_int_underground_obj(cube_t const &c) const { // Note: not useful for generating buildings (ext basements) because pools are added later
	for (swimming_pool_t const &p : pools) {
		if (!p.above_ground && c.intersects(p.bcube)) return 1; // zvals are checked
	}
	return 0;
}
void city_obj_placer_t::get_ponds_in_xy_range(cube_t const &range, vect_cube_t &pond_bcs) const {
	if (!range.intersects_xy(pond_groups.get_bcube())) return;

	for (pond_t const &p : ponds) {
		if (range.intersects_xy(p.bcube)) {pond_bcs.push_back(p.bcube);}
	}
}

bool city_obj_placer_t::update_depth_if_underwater(point const &pos, float &depth) const {
	if (!all_objs_bcube.contains_pt(pos)) return 0;

	// both pools and ponds are non-overlapping in XY, so if we find an intersection, we can stop and return the depth
	for (swimming_pool_t const &p : pools) { // not many pools, don't need to use pool_groups
		if (p.update_depth_if_underwater(pos, depth)) return 1;
	}
	for (pond_t const &p : ponds) {
		if (p.update_depth_if_underwater(pos, depth)) return 1;
	}
	return 0;
}

void city_obj_placer_t::play_sounds() {
	if (ppole_groups.empty()) return;
	point const camera_bs(get_camera_building_space());
	if (!all_objs_bcube.contains_pt(camera_bs)) return; // player not in this city
	float const sound_dist(0.65*city_params.road_width);
	cube_t test_cube;
	test_cube.set_from_sphere(camera_bs, sound_dist);
	unsigned start_ix(0);
	
	for (auto g = ppole_groups.begin(); g != ppole_groups.end(); start_ix = g->ix, ++g) {
		if (!g->intersects(test_cube)) continue;
		assert(start_ix <= g->ix && g->ix <= ppoles.size());
		
		for (auto i = ppoles.begin()+start_ix; i != ppoles.begin()+g->ix; ++i) {
			if (!i->has_transformer() || !i->bcube.intersects(test_cube)) continue;
			point const transformer_pos(i->get_transformer_center());
			float const dist(p2p_dist(camera_bs, transformer_pos));
			if (dist > sound_dist) continue;
			play_hum_sound(transformer_pos, (1.0 - dist/sound_dist), 0.6); // 60Hz
		}
	} // for g
}


