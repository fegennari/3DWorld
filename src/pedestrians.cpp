// 3D World - Pedestrians for Procedural Cities
// by Frank Gennari
// 12/6/18
#include "city.h"

float const PED_WIDTH_SCALE  = 0.5;
float const PED_HEIGHT_SCALE = 2.5;

extern int animate2;
extern city_params_t city_params;


string pedestrian_t::str() const {
	std::ostringstream oss;
	oss << TXT(vel.mag()) << TXT(radius) << TXT(city) << TXT(plot) << TXTi(stuck_count) << TXT(collided); // Note: pos, vel, dir not printed
	return oss.str();
}

bool pedestrian_t::is_valid_pos(cube_t const &plot_cube, vector<cube_t> const &colliders) const {
	if (!plot_cube.contains_pt_xy(pos)) return 0; // outside the plot (not yet allowing crossing roads at crosswalks)
	if (check_buildings_ped_coll(pos, radius, plot)) return 0;
	float const xmin(pos.x - radius), xmax(pos.x + radius);

	for (auto i = colliders.begin(); i != colliders.end(); ++i) {
		if (i->x2() < xmin) continue; // to the left
		if (i->x1() > xmax) break; // to the right - sorted from left to right, so no more colliders can intersect - done
		if (sphere_cube_intersect(pos, radius, *i)) return 0;
	}
	return 1;
}

bool pedestrian_t::check_ped_ped_coll(vector<pedestrian_t> &peds, unsigned pid) const {
	assert(pid < peds.size());

	for (auto i = peds.begin()+pid+1; i != peds.end(); ++i) { // check every ped after this one
		if (i->plot != plot || i->city != city) break; // moved to a new plot or city, no collision, done
		if (dist_xy_less_than(pos, i->pos, (radius + i->radius))) {i->collided  = 1; return 1;} // collision
	}
	return 0;
}

bool pedestrian_t::try_place_in_plot(cube_t const &plot_cube, vector<cube_t> const &colliders, unsigned plot_id, rand_gen_t &rgen) {
	pos    = rand_xy_pt_in_cube(plot_cube, radius, rgen);
	pos.z += radius; // place on top of the plot
	plot   = plot_id;
	if (!is_valid_pos(plot_cube, colliders)) return 0; // failed
	return 1; // success
}

void pedestrian_t::next_frame(cube_t const &plot_cube, vector<cube_t> const &colliders, vector<pedestrian_t> &peds, unsigned pid, rand_gen_t &rgen, float delta_dir) {
	if (vel == zero_vector) return; // not moving
	point const prev_pos(pos); // assume this ped starts out not colliding
	move();

	if (!collided && is_valid_pos(plot_cube, colliders) && !check_ped_ped_coll(peds, pid)) {stuck_count = 0;} // no collisions, nothing else to do
	else {
		collided = 0; // reset for next frame
		pos = prev_pos; // restore to previous valid pos
		if (++stuck_count > 8) {pos += rgen.signed_rand_vector_spherical_xy()*(0.1*radius);} // shift randomly by 10% radius to get unstuck
		vector3d new_vel(rgen.signed_rand_vector_spherical_xy()); // try a random new direction
		if (dot_product(vel, new_vel) > 0.0) {new_vel *= -1.0;} // negate if pointing in the same dir
		vel = new_vel * (vel.mag()/new_vel.mag()); // normalize to original velocity
	}
	dir = (delta_dir/vel.mag())*vel + (1.0 - delta_dir)*dir; // merge velocity into dir gradually for smooth turning
	dir.normalize();
}


unsigned ped_model_loader_t::num_models() const {return city_params.ped_model_files.size();}

city_model_t const &ped_model_loader_t::get_model(unsigned id) const {
	assert(id < num_models());
	return city_params.ped_model_files[id];
}


float ped_manager_t::get_ped_radius() const {return 0.05*city_params.road_width;} // or should this be relative to player/camera radius?

void ped_manager_t::expand_cube_for_ped(cube_t &cube) const {
	float const radius(get_ped_radius());
	cube.expand_by_xy(radius); // PED_WIDTH_SCALE*radius for models?
	cube.z2() += PED_HEIGHT_SCALE*radius;
}

void ped_manager_t::init(unsigned num) {
	if (num == 0) return;
	timer_t timer("Gen Peds");
	peds.reserve(num);
	float const radius(get_ped_radius()); // currently, all pedestrians are the same size
	unsigned const num_models(ped_model_loader.num_models());

	for (unsigned n = 0; n < num; ++n) {
		pedestrian_t ped(radius); // Note: constant radius for now, but can vary this to create peds of different sizes

		if (gen_ped_pos(ped)) {
			if (city_params.ped_speed > 0.0) {ped.vel = rgen.signed_rand_vector_spherical_xy(city_params.ped_speed);}
			ped.model_id = ((num_models > 0) ? (rgen.rand()%num_models) : 0);
			peds.push_back(ped);
		}
	}
	cout << "Pedestrians: " << peds.size() << endl; // testing
	if (peds.empty()) return;
	sort(peds.begin(), peds.end());

	// construct by_city
	unsigned const max_city(peds.back().city), max_plot(peds.back().plot);
	by_city.resize(max_city + 2); // one per city + terminator

	for (unsigned city = 0, pix = 0; city <= max_city; ++city) {
		while (pix < peds.size() && peds[pix].city == city) {++pix;}
		unsigned const cur_plot((pix < peds.size()) ? peds[pix].plot : max_plot+1);
		by_city[city+1].assign(pix, cur_plot); // next city begins here
	}
	// construct by_plot
	by_plot.resize((max_plot + 2), 0); // one per by_plot + terminator

	for (unsigned plot = 0, pix = 0; plot <= max_plot; ++plot) {
		while (pix < peds.size() && peds[pix].plot == plot) {++pix;}
		by_plot[plot+1] = pix; // next plot begins here
	}
}

bool ped_manager_t::proc_sphere_coll(point &pos, float radius, vector3d *cnorm) const { // Note: no p_last; for potential use with ped/ped collisions
	float const rsum(get_ped_radius() + radius);

	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		cube_t const city_bcube(get_expanded_city_bcube_for_peds(city));
		if (pos.z > city_bcube.z2() + rsum) continue; // above the peds
		if (!sphere_cube_intersect_xy(pos, radius, city_bcube)) continue;
		unsigned const plot_start(by_city[city].plot_ix), plot_end(by_city[city+1].plot_ix);

		for (unsigned plot = plot_start; plot < plot_end; ++plot) {
			cube_t const plot_bcube(get_expanded_city_plot_bcube_for_peds(city, plot));
			if (!sphere_cube_intersect_xy(pos, radius, plot_bcube)) continue;
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				if (!dist_less_than(pos, peds[i].pos, rsum)) continue;
				if (cnorm) {*cnorm = (pos - peds[i].pos).get_norm();}
				return 1; // return on first coll
			}
		} // for plot
	} // for city
	return 0;
}

bool ped_manager_t::line_intersect_peds(point const &p1, point const &p2, float &t) const {
	bool ret(0);

	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		if (!get_expanded_city_bcube_for_peds(city).line_intersects(p1, p2)) continue;
		unsigned const plot_start(by_city[city].plot_ix), plot_end(by_city[city+1].plot_ix);

		for (unsigned plot = plot_start; plot < plot_end; ++plot) {
			if (!get_expanded_city_plot_bcube_for_peds(city, plot).line_intersects(p1, p2)) continue;
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				float tmin(0.0);
				if (line_sphere_int_closest_pt_t(p1, p2, peds[i].pos, peds[i].radius, tmin) && tmin < t) {t = tmin; ret = 1;}
			}
		} // for plot
	} // for city
	return ret;
}

void ped_manager_t::next_frame() {
	if (!animate2) return;
	//timer_t timer("Ped Update"); // ~1.6ms for 10K peds
	float const delta_dir(1.0 - pow(0.7f, fticks)); // controls pedestrian turning rate

	for (auto i = peds.begin(); i != peds.end(); ++i) {
		// FIXME_PEDS: navigation with destination
		// TODO: road_gen.get_city().get_isec().mark_crosswalk_in_use(dim, dir)
		cube_t const plot(get_city_plot_bcube_for_peds(i->city, i->plot));
		auto const &colliders(get_colliders_for_plot(i->city, i->plot));
		i->next_frame(plot, colliders, peds, (i - peds.begin()), rgen, delta_dir);
	} // for i
}

void being_sphere_draw(shader_t &s, bool &in_sphere_draw, bool textured) {
	if (in_sphere_draw) return;
	if (!textured) {select_texture(WHITE_TEX);} // currently not textured
	s.set_cur_color(YELLOW); // smileys
	begin_sphere_draw(textured);
	in_sphere_draw = 1;
}
void end_sphere_draw(bool &in_sphere_draw) {
	if (!in_sphere_draw) return;
	end_sphere_draw();
	in_sphere_draw = 0;
}

void ped_manager_t::draw(vector3d const &xlate, bool use_dlights, bool shadow_only, bool is_dlight_shadows) {
	if (empty()) return;
	if (is_dlight_shadows && !city_params.car_shadows) return; // use car_shadows as ped_shadows
	if (shadow_only && !is_dlight_shadows) return; // don't add to precomputed shadow map
	//timer_t timer("Ped Draw"); // ~1ms
	bool const use_models(ped_model_loader.num_models() > 0);
	float const draw_dist(is_dlight_shadows ? 0.8*camera_pdu.far_ : (use_models ? 500.0 : 2000.0)*get_ped_radius()); // smaller view dist for models
	pos_dir_up pdu(camera_pdu); // decrease the far clipping plane for pedestrians
	pdu.far_ = draw_dist;
	pdu.pos -= xlate; // adjust for local translate
	dstate.xlate = xlate;
	fgPushMatrix();
	translate_to(xlate);
	dstate.pre_draw(xlate, use_dlights, shadow_only, 1); // always_setup_shader=1
	bool const textured(shadow_only && 0); // disabled for now
	bool in_sphere_draw(0);

	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		if (!pdu.cube_visible(get_expanded_city_bcube_for_peds(city))) continue; // city not visible - skip
		unsigned const plot_start(by_city[city].plot_ix), plot_end(by_city[city+1].plot_ix);
		assert(plot_start <= plot_end);

		for (unsigned plot = plot_start; plot < plot_end; ++plot) {
			assert(plot < by_plot.size());
			cube_t const plot_bcube(get_expanded_city_plot_bcube_for_peds(city, plot));
			if (!pdu.cube_visible(plot_bcube)) continue; // plot not visible - skip
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);
			assert(ped_start <= ped_end);
			if (ped_start == ped_end) continue; // no peds on this plot
			dstate.ensure_shader_active(); // needed for use_smap=0 case
			if (!shadow_only) {dstate.begin_tile(plot_bcube.get_cube_center(), 1);} // use the plot's tile's shadow map

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				pedestrian_t const &ped(peds[i]);
				assert(ped.city == city && ped.plot == plot);
				if (!dist_less_than(pdu.pos, ped.pos, draw_dist)) continue; // too far - skip

				if (!use_models) { // or distant?
					if (!pdu.sphere_visible_test(ped.pos, ped.radius)) continue; // not visible - skip
					being_sphere_draw(dstate.s, in_sphere_draw, textured);
					int const ndiv = 16; // currently hard-coded
					draw_sphere_vbo(ped.pos, ped.radius, ndiv, textured);
				}
				else {
					float const width(PED_WIDTH_SCALE*ped.radius), height(PED_HEIGHT_SCALE*ped.radius);
					cube_t bcube;
					bcube.set_from_sphere(ped.pos, width);
					bcube.z1() = ped.pos.z - ped.radius;
					bcube.z2() = bcube.z1() + height;
					if (!pdu.sphere_visible_test(bcube.get_cube_center(), 0.5*height)) continue; // not visible - skip
					end_sphere_draw(in_sphere_draw);
					bool const low_detail(!shadow_only && !dist_less_than(pdu.pos, ped.pos, 0.5*draw_dist)); // low detail for non-shadow pass at half draw dist
					ped_model_loader.draw_model(dstate.s, ped.pos, bcube, ped.dir, ALPHA0, xlate, ped.model_id, shadow_only, low_detail);
				}
			} // for i
		} // for plot
	} // for city
	end_sphere_draw(in_sphere_draw);
	dstate.end_draw();
	fgPopMatrix();
}

