// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "lightmap.h" // for light_source
#include "cobj_bsp_tree.h"
#include <thread>

bool const USE_BKG_THREAD = 1;

extern int MESH_Z_SIZE, display_mode, display_framerate, camera_surf_collide, animate2;
extern unsigned LOCAL_RAYS, MAX_RAY_BOUNCES, NUM_THREADS;
extern float indir_light_exp;
extern std::string lighting_update_text;
extern vector<light_source> dl_sources;

bool enable_building_people_ai();


bool ray_cast_cube(point const &p1, point const &p2, cube_t const &c, vector3d &cnorm, float &t) {
	float tmin(0.0), tmax(1.0);
	if (!get_line_clip(p1, p2, c.d, tmin, tmax) || tmin >= t) return 0;
	t = tmin;
	get_closest_cube_norm(c.d, (p1 + (p2 - p1)*t), cnorm);
	return 1;
}

class cube_bvh_t : public cobj_tree_simple_type_t<colored_cube_t> {
	virtual void calc_node_bbox(tree_node &n) const {
		assert(n.start < n.end);
		for (unsigned i = n.start; i < n.end; ++i) {n.assign_or_union_with_cube(objects[i]);} // bcube union
	}
public:
	vect_colored_cube_t &get_objs() {return objects;}

	bool ray_cast(point const &p1, point const &p2, vector3d &cnorm, colorRGBA &ccolor, float &t) const {
		if (nodes.empty()) return 0;
		bool ret(0);
		node_ix_mgr nixm(nodes, p1, p2);
		unsigned const num_nodes((unsigned)nodes.size());

		for (unsigned nix = 0; nix < num_nodes;) {
			tree_node const &n(nodes[nix]);
			if (!nixm.check_node(nix)) continue; // Note: modifies nix

			for (unsigned i = n.start; i < n.end; ++i) { // check leaves
				if (ray_cast_cube(p1, p2, objects[i], cnorm, t)) {ccolor = objects[i].color; ret = 1;}
			}
		}
		return ret;
	}
};

bool follow_ray_through_cubes_recur(point const &p1, point const &p2, point const &start, vect_cube_t const &cubes,
	vect_cube_t::const_iterator end, vect_cube_t::const_iterator parent, vector3d &cnorm, float &t)
{
	for (auto c = cubes.begin(); c != end; ++c) {
		if (c == parent || !c->contains_pt(start)) continue;
		float tmin(0.0), tmax(1.0);
		if (!get_line_clip(p1, p2, c->d, tmin, tmax) || tmax >= t) continue;
		point const cpos(p1 + (p2 - p1)*(tmax + 0.001)); // move slightly beyond the hit point
		if ((end - cubes.begin()) > 1 && follow_ray_through_cubes_recur(p1, p2, cpos, cubes, end, c, cnorm, t)) return 1;
		t = tmax;
		get_closest_cube_norm(c->d, (p1 + (p2 - p1)*t), cnorm); // this is the final point, update cnorm
		cnorm.negate(); // reverse hit dir
		return 1;
	} // for c
	return 0;
}

// Note: static objects only; excludes people; pos in building space
bool building_t::ray_cast_interior(point const &pos, vector3d const &dir, cube_bvh_t const &bvh, point &cpos, vector3d &cnorm, colorRGBA &ccolor) const {

	if (!interior || is_rotated() || !is_simple_cube()) return 0; // these cases are not yet supported
	float const extent(bcube.get_max_extent());
	cube_t clip_cube(bcube);
	clip_cube.expand_by(0.01*extent); // expand slightly so that collisions with objects on the edge are still considered interior
	point p1(pos), p2(pos + dir*(2.0*extent));
	if (!do_line_clip(p1, p2, clip_cube.d)) return 0; // ray does not intersect building bcube
	building_mat_t const &mat(get_material());
	float t(1.0); // start at p2
	bool hit(0);
	cpos = p2; // use far clip point for clip cube if there is no hit

	// check parts (exterior walls); should chimneys and porch roofs be included?
	if (follow_ray_through_cubes_recur(p1, p2, p1, parts, get_real_parts_end_inc_sec(), parts.end(), cnorm, t)) { // interior ray - find furthest exit point
		ccolor = mat.wall_color.modulate_with(mat.wall_tex.get_avg_color());
		p2  = p1 + (p2 - p1)*t; t = 1.0; // clip p2 to t (minor optimization)
		hit = 1;
	}
	else { // check for exterior rays
		bool hit(0);
		for (auto p = parts.begin(); p != parts.end(); ++p) {hit |= ray_cast_cube(p1, p2, *p, cnorm, t);} // find closest entrance point
		
		if (hit) { // exterior hit - don't need to check interior geometry
			ccolor = side_color.modulate_with(mat.side_tex.get_avg_color());
			cpos   = p1 + (p2 - p1)*t;
			return 1;
		}
	}
	//for (auto r = roof_tquads.begin(); r != roof_tquads.end(); ++r) {} // WRITE; use roof_color/mat.roof_tex; only needed for exterior rays?
	bvh.ray_cast(p1, p2, cnorm, ccolor, t);
	if (t == 1.0) {cpos = p2; return hit;} // no intersection with bvh
	cpos = p1 + (p2 - p1)*t;
	return 1;
}

unsigned elevator_t::get_coll_cubes(cube_t cubes[5]) const {
	if (!open) {cubes[0] = *this; return 1;} // closed, entire elevator is 1 cube
	float const wwidth(get_wall_thickness()), fwidth(get_frame_width());
	for (unsigned i = 0; i < 5; ++i) {cubes[i] = *this;} // start with full cube
	cubes[0].d[dim][dir] = cubes[0].d[dim][!dir] + (dir ? 1.0 : -1.0)*wwidth; // back
	
	for (unsigned d = 0; d < 2; ++d) {
		cube_t &side(cubes[d+1]), &front(cubes[d+3]);
		side.d[!dim][d] = side.d[!dim][!d] + (d ? 1.0 : -1.0)*wwidth; // sides
		front.d[dim][!dir] = front.d[dim][dir] + (dir ? -1.0 : 1.0)*wwidth; // front side parts
		front.d[!dim][d] = front.d[!dim][!d] + (d ? 1.0 : -1.0)*fwidth;
	}
	return 5; // back + 2 sides + 2 front parts
}

template<typename T> void add_colored_cubes(vector<T> const &cubes, colorRGBA const &color, vect_colored_cube_t &cc) {
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {cc.emplace_back(*c, color);}
}
void building_t::gather_interior_cubes(vect_colored_cube_t &cc) const {

	if (!interior) return; // nothing to do
	building_mat_t const &mat(get_material());
	colorRGBA const wall_color(mat.wall_color.modulate_with(mat.wall_tex.get_avg_color()));
	for (unsigned d = 0; d < 2; ++d) {add_colored_cubes(interior->walls[d], wall_color, cc);}

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		cube_t cubes[5];
		unsigned const num_cubes(e->get_coll_cubes(cubes));
		// for now elevators are treated the same as walls with the same color, even though the inside of open elevators is wood
		for (unsigned n = 0; n < num_cubes; ++n) {cc.emplace_back(cubes[n], wall_color);} // can only assign the same color to all sides of the cube
	}
	add_colored_cubes(interior->ceilings, mat.ceil_color .modulate_with(mat.ceil_tex .get_avg_color()), cc);
	add_colored_cubes(interior->floors,   mat.floor_color.modulate_with(mat.floor_tex.get_avg_color()), cc);
	add_colored_cubes(details,            detail_color.   modulate_with(mat.roof_tex. get_avg_color()), cc); // should this be included?
	if (!has_room_geom()) return; // nothing else to add
	vector<room_object_t> const &objs(interior->room_geom->objs);
	cc.reserve(cc.size() + objs.size());
		
	for (auto c = objs.begin(); c != objs.end(); ++c) {
		if (c->shape == SHAPE_CYLIN  ) continue; // cylinders (lights) are not cubes
		if (c->type  == TYPE_ELEVATOR) continue; // elevator cars/internals can move so should not contribute to lighting
		if (c->type  == TYPE_BLOCKER ) continue; // blockers are not drawn
		// to be more accurate, we could use the actual cubes of tables and chairs, but this adds a lot of complexity, increases lighting time, and makes little difference
		cc.emplace_back(*c, c->get_color());
		if (c->type == TYPE_TABLE) {cc.back().z1() += 0.85*c->dz();} // at least be a bit more accurate for tables by using only the top
	}
}


class building_indir_light_mgr_t {
	bool is_running, is_done, kill_thread, lighting_updated, needs_to_join;
	int cur_bix, cur_light;
	unsigned cur_tid;
	vector<unsigned char> tex_data;
	vector<unsigned> light_ids;
	set<unsigned> lights_complete;
	cube_bvh_t bvh;
	lmap_manager_t lmgr;
	std::thread rt_thread;

	void init_lmgr(bool clear_lighting) {
		if (clear_lighting) {lmgr.reset_all();}
		if (lmgr.is_allocated()) return; // already setup
		unsigned const tot_sz(XY_MULT_SIZE*MESH_SIZE[2]); // Note: MESH_SIZE[2], not MESH_Z_SIZE; want clipped size that lmap uses rather than user-specified size
		lmcell init_lmcell;
		lmgr.alloc(tot_sz, MESH_X_SIZE, MESH_Y_SIZE, MESH_SIZE[2], (unsigned char **)nullptr, init_lmcell);
	}
	void start_lighting_compute(building_t const &b) {
		assert(cur_light >= 0);
		init_lmgr(0); // clear_lighting=0
		is_running = 1;
		lighting_updated = 1;

		if (USE_BKG_THREAD) { // start a thread to compute cur_light for building b
			rt_thread = std::thread(&building_indir_light_mgr_t::cast_light_ray, this, b);
			needs_to_join = 1;
		}
		else {
			timer_t timer("Ray Cast Building Light");
			cast_light_ray(b);
		}
	}
	void calc_reflect_ray(point &pos, point const &cpos, vector3d &dir, vector3d const &cnorm, rand_gen_t &rgen, float tolerance) const {
		vector3d v_ref;
		calc_reflection_angle(dir, v_ref, cnorm);
		v_ref.normalize();
		vector3d const rand_dir(rgen.signed_rand_vector().get_norm());
		dir = (v_ref + rand_dir).get_norm(); // diffuse reflection: new dir is mix 50% specular with 50% random
		if (dot_product(dir, cnorm) < 0.0) {dir.negate();} // make sure it points away from the surface (is this needed?)
		pos = cpos + tolerance*dir; // move slightly away from the surface
	}
	void cast_light_ray(building_t const &b) {
		// Note: modifies lmgr, but otherwise thread safe
		unsigned const num_rt_threads(NUM_THREADS - (USE_BKG_THREAD ? 1 : 0)); // reserve a thread for the main thread if running in the background
		vector<room_object_t> const &objs(b.interior->room_geom->objs);
		assert((unsigned)cur_light < objs.size());
		room_object_t const &ro(objs[cur_light]);
		colorRGBA const lcolor(ro.get_color());
		cube_t const scene_bounds(get_scene_bounds_bcube()); // expected by lmap update code
		point const ray_scale(scene_bounds.get_size()/b.bcube.get_size()), llc_shift(scene_bounds.get_llc() - b.bcube.get_llc()*ray_scale);
		float const tolerance(1.0E-5*b.bcube.get_max_extent()), light_zval(ro.z1() - 0.01*ro.dz()); // set slightly below bottom of light
		float const surface_area(ro.dx()*ro.dy() + 2.0f*(ro.dx() + ro.dy())*ro.dz()); // bottom + 4 sides (top is occluded), 0.0003 for houses
		float weight(100.0f*(surface_area/0.0003f)/LOCAL_RAYS); // normalize to the number of rays
		if (b.has_pri_hall()) {weight *= 0.8;} // floorplan is open and well lit, indir lighting value seems too high
		if (b.is_house) {weight *= 2.0;} // houses have dimmer lights and seem to work better with more indir
		unsigned const NUM_PRI_SPLITS = 16;
		int const num_rays(LOCAL_RAYS/NUM_PRI_SPLITS);

#pragma omp parallel for schedule(dynamic) num_threads(num_rt_threads)
		for (int n = 0; n < num_rays; ++n) {
			if (kill_thread) continue;
			rand_gen_t rgen;
			rgen.set_state(n+1, cur_light);
			vector3d pri_dir(rgen.signed_rand_vector_spherical(1.0).get_norm());
			pri_dir.z = -fabs(pri_dir.z); // make sure dir points down
			point origin, init_cpos, cpos;
			vector3d init_cnorm, cnorm;
			colorRGBA ccolor(WHITE);
			// select a random point on the light cube (close enough for (ro.shape == SHAPE_CYLIN))
			for (unsigned d = 0; d < 2; ++d) {origin[d] = rgen.rand_uniform(ro.d[d][0], ro.d[d][1]);}
			origin.z = light_zval;
			init_cpos = origin; // init value
			if (!b.ray_cast_interior(origin, pri_dir, bvh, init_cpos, init_cnorm, ccolor)) continue;
			colorRGBA const init_color(lcolor.modulate_with(ccolor));
			if (init_color.get_luminance() < 0.1) continue; // done

			for (unsigned splits = 0; splits < NUM_PRI_SPLITS; ++splits) {
				point pos(origin);
				vector3d dir(pri_dir);
				colorRGBA cur_color(init_color);
				calc_reflect_ray(pos, init_cpos, dir, init_cnorm, rgen, tolerance);

				for (unsigned bounce = 1; bounce < MAX_RAY_BOUNCES; ++bounce) { // allow up to MAX_RAY_BOUNCES bounces
					cpos = pos; // init value
					bool const hit(b.ray_cast_interior(pos, dir, bvh, cpos, cnorm, ccolor));

					if (cpos != pos) { // accumulate light along the ray from pos to cpos (which is always valid) with color cur_color
						point const p1(pos*ray_scale + llc_shift), p2(cpos*ray_scale + llc_shift); // transform building space to global scene space
						add_path_to_lmcs(&lmgr, nullptr, p1, p2, weight, cur_color, LIGHTING_LOCAL, 0); // local light, no bcube
					}
					if (!hit) break; // done
					cur_color = cur_color.modulate_with(ccolor);
					if (cur_color.get_luminance() < 0.1) break; // done
					calc_reflect_ray(pos, cpos, dir, cnorm, rgen, tolerance);
				} // for bounce
			} // for splits
		} // for n
		is_running = 0;
	}
	void wait_for_finish(bool force_kill) {
		// Note: for now the time taken to process a light should be pretty fast so we just block until finished; set kill_thread=1 to be faster
		if (force_kill) {kill_thread = 1;}
		while (is_running) {alut_sleep(0.01);}
		kill_thread = 0;
	}
	void update_volume_light_texture() { // full update, 6.6ms for z=128
		//timer_t timer("Lighting Tex Create");
		indir_light_tex_from_lmap(cur_tid, lmgr, tex_data, MESH_X_SIZE, MESH_Y_SIZE, MESH_SIZE[2], indir_light_exp, 1); // local_only=1
	}
	void maybe_join_thread() {
		if (needs_to_join) {rt_thread.join(); needs_to_join = 0;}
	}
public:
	building_indir_light_mgr_t() : is_running(0), is_done(0), kill_thread(0), lighting_updated(0), needs_to_join(0), cur_bix(-1), cur_light(-1), cur_tid(0) {}

	void clear() {
		is_done = lighting_updated = 0;
		cur_bix = cur_light = -1;
		tex_data.clear();
		light_ids.clear();
		lights_complete.clear();
		end_rt_job();
		lmgr.reset_all(); // clear lighting values back to 0
		bvh.clear();
	}
	void end_rt_job() {
		wait_for_finish(1); // force_kill=1
		maybe_join_thread();
	}
	void free_indir_texture() {free_texture(cur_tid);}

	void register_cur_building(building_t const &b, unsigned bix, point const &target, unsigned &tid) { // target is in building space
		if ((int)bix != cur_bix) { // change to a different building
			clear();
			cur_bix = bix;
			assert(!is_running);
			build_bvh(b);
		}
		if (cur_tid > 0 && is_done) return; // nothing else to do

		if (display_framerate && (is_running || lighting_updated)) { // show progress to the user
			std::ostringstream oss;
			oss << "Lights: " << lights_complete.size() << " / " << light_ids.size();
			lighting_update_text = oss.str();
		}
		if (is_running) return; // still running, let it continue

		if (lighting_updated) { // update lighting texture based on incremental progress
			maybe_join_thread();
			update_volume_light_texture();
			lighting_updated = 0;
		}
		// nothing is running and there is more work to do, find the nearest light to the target and process it
		if (cur_light >= 0) {lights_complete.insert(cur_light); cur_light = -1;} // mark the most recent light as complete
		b.order_lights_by_priority(target, light_ids);

		for (auto i = light_ids.begin(); i != light_ids.end(); ++i) {
			if (lights_complete.find(*i) == lights_complete.end()) {cur_light = *i; break;} // find an incomplete light
		}
		if (cur_light >= 0) {start_lighting_compute(b);} // this light is next
		else {is_done = 1;} // no more lights to process
		//cout << "Process light " << lights_complete.size() << " of " << light_ids.size() << endl;
		tid = cur_tid;
	}
	void build_bvh(building_t const &b) {
		bvh.clear();
		b.gather_interior_cubes(bvh.get_objs());
		bvh.build_tree_top(0); // verbose=0
	}
	cube_bvh_t const &get_bvh() const {return bvh;}
};

building_indir_light_mgr_t building_indir_light_mgr;

void free_building_indir_texture() {building_indir_light_mgr.free_indir_texture();}
void end_building_rt_job() {building_indir_light_mgr.end_rt_job();}

void building_t::create_building_volume_light_texture(unsigned bix, point const &target, unsigned &tid) const {
	if (!has_room_geom()) return; // error?
	building_indir_light_mgr.register_cur_building(*this, bix, target, tid);
}

bool building_t::ray_cast_camera_dir(point const &camera_bs, point &cpos, colorRGBA &ccolor) const {
	assert(!USE_BKG_THREAD); // not legal to call when running lighting in a background thread
	building_indir_light_mgr.build_bvh(*this);
	vector3d cnorm; // unused
	return ray_cast_interior(camera_bs, cview_dir, building_indir_light_mgr.get_bvh(), cpos, cnorm, ccolor);
}

void building_t::order_lights_by_priority(point const &target, vector<unsigned> &light_ids) const {

	light_ids.clear();
	if (!has_room_geom()) return; // error?
	vector<room_object_t> const &objs(interior->room_geom->objs);
	vector<pair<float, unsigned>> to_sort;
	float const window_vspacing(get_window_vspace());
	float const diag_dist_sq(bcube.dx()*bcube.dx() + bcube.dy()*bcube.dy()), other_floor_penalty(0.25*diag_dist_sq);

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (i->type != TYPE_LIGHT || !i->is_lit()) continue; // not a light, or light not on
		float dist_sq(p2p_dist_sq(i->get_cube_center(), target));
		dist_sq *= 0.005f*window_vspacing/(i->dx()*i->dy()); // account for the size of the light, larger lights smaller/higher priority

		if (i->z1() < target.z || i->z2() > (target.z + window_vspacing)) { // penalty if on a different floor
			dist_sq += (i->has_stairs() ? 0.25 : 1.0)*other_floor_penalty; // less penalty for lights on stairs
		}
		// reduce distance for lights visible to target?
		to_sort.emplace_back(dist_sq, (i - objs.begin()));
	} // for i
	sort(to_sort.begin(), to_sort.end()); // sort by increasing distance
	for (auto i = to_sort.begin(); i != to_sort.end(); ++i) {light_ids.push_back(i->second);}
}

bool line_int_cubes(point const &p1, point const &p2, vect_cube_t const &cubes) {
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (c->line_intersects(p1, p2)) return 1;
	}
	return 0;
}
bool building_t::is_light_occluded(point const &lpos, point const &camera_bs) const {
	// Note: assumes the light is inside the building
	// exterior walls have windows and don't generally occlude lights; room objects and doors are too small to occlude; elevators are too sparse to occlude
	for (unsigned d = 0; d < 2; ++d) {
		if (line_int_cubes(lpos, camera_bs, interior->walls[d])) return 1;
	}
	if (line_int_cubes(lpos, camera_bs, interior->floors  )) return 1;
	if (line_int_cubes(lpos, camera_bs, interior->ceilings)) return 1;
	return 0;
}
void building_t::clip_ray_to_walls(point const &p1, point &p2) const { // Note: assumes p1.z == p2.z
	float t(1.0), tmin(0.0), tmax(1.0);

	for (unsigned d = 0; d < 2; ++d) {
		for (auto c = interior->walls[d].begin(); c != interior->walls[d].end(); ++c) {
			if (p1.z < c->z1() || p1.z > c->z2()) continue; // no z overlap
			if (get_line_clip_xy(p1, p2, c->d, tmin, tmax) && tmin < t) {t = tmin;} // optimization: only clip p2?
		}
	}
	if (t < 1.0) {p2 = p1 + (p2 - p1)*t;} // clip to t
}

void building_t::refine_light_bcube(point const &lpos, float light_radius, cube_t &light_bcube) const {
	// base: 173613 / bcube: 163942 / clipped bcube: 161455 / tight: 159005 / rays: 101205 / no ls bcube expand: 74538
	// starts with building bcube clipped to light bcube
	//timer_t timer("refine_light_bcube");
	cube_t tight_bcube;

	// first determine the union of all intersections with parts; ignore zvals here so that we get the same result for every floor
	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		//if (lpos.z < p->z1() || lpos.z > p->z2()) continue; // light zval not included (doesn't seem to work)
		if (!light_bcube.intersects_xy(*p)) continue;
		cube_t c(light_bcube);
		c.intersect_with_cube_xy(*p);
		if (tight_bcube.is_all_zeros()) {tight_bcube = c;} else {tight_bcube.union_with_cube_xy(c);}
	}
	tight_bcube.z1() = light_bcube.z1();
	tight_bcube.z2() = light_bcube.z2();
	// next cast a number of horizontal rays in a circle around the light to see how far they reach; any walls hit occlude the light and reduce the bcube
	unsigned const NUM_RAYS = 180; // every 2 degrees
	if (NUM_RAYS == 0) {light_bcube = tight_bcube; return;}
	cube_t rays_bcube(lpos, lpos);

	for (unsigned n = 0; n < NUM_RAYS; ++n) {
		float const angle(TWO_PI*n/NUM_RAYS), dx(light_radius*sin(angle)), dy(light_radius*cos(angle));
		point p1(lpos), p2(lpos + point(dx, dy, 0.0));
		bool const ret(do_line_clip(p1, p2, tight_bcube.d));
		assert(ret);
		assert(p1 == lpos);
		clip_ray_to_walls(p1, p2);
		rays_bcube.union_with_pt(p2);
	} // for n
	light_bcube = rays_bcube;
}

// Note: non const because this caches light_bcubes
void building_t::add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building, int ped_ix, vect_cube_t &ped_bcubes, cube_t &lights_bcube) {

	if (!has_room_geom()) return; // error?
	vector<room_object_t> &objs(interior->room_geom->objs); // non-const, light flags are updated
	vect_cube_t &light_bcubes(interior->room_geom->light_bcubes);
	point const camera_bs(camera_pdu.pos - xlate); // camera in building space
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), camera_z(camera_bs.z);
	assert(interior->room_geom->stairs_start <= objs.size());
	auto objs_end(objs.begin() + interior->room_geom->stairs_start); // skip stairs and elevators
	unsigned camera_part(parts.size()); // start at an invalid value
	bool camera_by_stairs(0), camera_near_building(camera_in_building);
	int camera_room(-1);
	ped_bcubes.clear();

	if (camera_in_building) {
		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (i->contains_pt(camera_bs)) {camera_part = (i - parts.begin()); break;}
		}
		for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) { // conservative but less efficient
			cube_t tc(*r);
			tc.expand_by_xy(wall_thickness); // to include camera in doorway

			// Note: stairs that connect stacked parts aren't flagged with has_stairs because stairs only connect to the bottom floor, but they're partially handled below
			if (tc.contains_pt(camera_bs)) {
				camera_room = (r - interior->rooms.begin());
				camera_by_stairs = r->has_stairs;
				assert(r->rtype < NUM_RTYPES);
				if (display_mode & 0x20) {lighting_update_text = room_names[r->rtype];} // debugging, key '6'
				break;
			}
		}
	}
	else {
		cube_t bcube_exp(bcube);
		bcube_exp.expand_by_xy(2.0*window_vspacing);
		camera_near_building = bcube_exp.contains_pt(camera_bs);
	}
	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_LIGHT || !i->is_lit()) continue; // not a light, or light not on
		point const lpos(i->get_cube_center()); // centered in the light fixture
		if (!lights_bcube.contains_pt_xy(lpos)) continue; // not contained within the light volume
		//if (is_light_occluded(lpos, camera_bs)) continue; // too strong a test in general, but may be useful for selecting high importance lights
		float const floor_z(i->z2() - window_vspacing), ceil_z(i->z2());
		bool const floor_is_above(camera_z < floor_z), floor_is_below(camera_z > ceil_z);
		assert(i->room_id < interior->rooms.size());
		room_t const &room(interior->rooms[i->room_id]);
		// less culling if either the light or the camera is by stairs and light is on the floor above or below
		bool stairs_light(0);

		if ((camera_z > floor_z-window_vspacing) && (camera_z < ceil_z+window_vspacing)) { // light is on the floor above or below the camera
			stairs_light = (i->has_stairs() || room.has_stairs || camera_by_stairs); // either the light or the camera is by the stairs

			if (!stairs_light /*&& floor_is_below*/ && camera_room >= 0) { // what about camera in room adjacent to one with stairs?
				cube_t cr(interior->rooms[camera_room]);
				cr.expand_by_xy(2.0*wall_thickness);

				for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
					if (r->has_stairs && r->intersects_no_adj(cr)) {stairs_light = 1; break;}
				}
			}
		}
		if (floor_is_above || floor_is_below) { // light is on a different floor from the camera
			if (camera_in_building && (room.part_id == camera_part ||
				(room.contains_pt_xy(camera_bs) && camera_z < ceil_z+window_vspacing && camera_z > floor_z-window_vspacing)))
			{
				// player is on a different floor of the same building part, or more than one floor away in a part stack, and can't see a light from the floor above/below
				if (!stairs_light) continue; // camera in building and on wrong floor, don't add light
				if (camera_z < (i->z2() - 2.0*window_vspacing) || camera_z > (i->z2() + window_vspacing)) continue; // light is on the stairs, add if one floor above/below
			}
			else { // camera outside the building (or the part that contains this light)
				float const xy_dist(p2p_dist_xy(camera_bs, lpos));
				if (!stairs_light && ((camera_z - lpos.z) > 2.0f*xy_dist || (lpos.z - camera_z) > 1.0f*xy_dist)) continue; // light viewed at too high an angle

				if (camera_in_building) { // camera and light are in different buildings/parts
					if (camera_part >= real_num_parts) continue; // camera in garage or shed
					assert(camera_part < parts.size());
					cube_t const &cpart(parts[camera_part]);
					if (cpart.z2() <= room.z1() || cpart.z1() >= room.z2()) continue; // light in a different vertical stack than the camera
					// is it better to check if light half sphere is occluded by the floor above/below?
					assert(room.part_id < parts.size());
					cube_t const &part(parts[room.part_id]);
					bool visible[2] = {0};

					for (unsigned d = 0; d < 2; ++d) { // for each dim
						bool const dir(camera_bs[d] > lpos[d]);
						if ((camera_bs[d] > part.d[d][dir]) ^ dir) continue; // camera not on the outside face of the part containing this room, so can't see through any windows
						visible[d] = (room.ext_sides & (1 << (2*d + dir)));
					}
					if (!visible[0] && !visible[1]) continue; // room is not on the exterior of the building on either side facing the camera
				}
			}
		} // end camera on different floor case
		float const light_radius(6.0f*(i->dx() + i->dy())), cull_radius(0.95*light_radius);
		if (!camera_pdu.sphere_visible_test((lpos + xlate), cull_radius)) continue; // VFC
		// check visibility of bcube of light sphere clipped to building bcube; this excludes lights behind the camera and improves shadow map assignment quality
		cube_t sphere_bc; // in building space
		sphere_bc.set_from_sphere(lpos, cull_radius);
		cube_t clipped_bc(sphere_bc);
		clipped_bc.intersect_with_cube(bcube);
		if (!stairs_light) {clipped_bc.z1() = floor_z; clipped_bc.z2() = ceil_z;} // clip zval to current floor if light not in a room with stairs or elevator
		if (!camera_pdu.cube_visible(clipped_bc + xlate)) continue; // VFC
		// update lights_bcube and add light(s)
		min_eq(lights_bcube.z1(), (lpos.z - light_radius));
		max_eq(lights_bcube.z2(), (lpos.z + 0.1f*light_radius)); // pointed down - don't extend as far up
		float const bwidth = 0.25; // as close to 180 degree FOV as we can get without shadow clipping
		colorRGBA const color(i->get_color()*1.1); // make it extra bright
		assert(i->obj_id < light_bcubes.size());
		cube_t &light_bcube(light_bcubes[i->obj_id]);

		if (light_bcube.is_all_zeros()) { // not yet calculated - calculate and cache
			light_bcube = clipped_bc;
			refine_light_bcube(lpos, light_radius, light_bcube);
		}
		clipped_bc.x1() = light_bcube.x1(); clipped_bc.x2() = light_bcube.x2(); // copy X/Y but keep orig zvals
		clipped_bc.y1() = light_bcube.y1(); clipped_bc.y2() = light_bcube.y2();
		clipped_bc.expand_by_xy(wall_thickness); // expand by wall thickness so that offset exterior doors are properly handled
		clipped_bc.intersect_with_cube(sphere_bc); // clip to original light sphere, which still applies (only need to expand at building exterior)
		if (!camera_pdu.cube_visible(clipped_bc + xlate)) continue; // VFC - post clip
		dl_sources.emplace_back(light_radius, lpos, lpos, color, 0, -plus_z, bwidth); // points down, white for now
		dl_sources.back().set_custom_bcube(clipped_bc);
		bool dynamic_shadows(0);

		if (camera_near_building) {
			if (camera_surf_collide && camera_in_building && lpos.z > camera_bs.z && clipped_bc.contains_pt(camera_bs)) {dynamic_shadows = 1;} // camera shadow
			else if (animate2 && enable_building_people_ai()) { // check moving people
				if (ped_ix >= 0 && ped_bcubes.empty()) {get_ped_bcubes_for_building(ped_ix, building_id, ped_bcubes);} // get cubes on first light

				for (auto c = ped_bcubes.begin(); c != ped_bcubes.end(); ++c) {
					if (lpos.z > c->z2() && c->intersects(clipped_bc)) {dynamic_shadows = 1; break;}
				}
			}
		}
		// if there are no dynamic shadows, we can reuse the previous frame's shadow map;
		// need to handle the case where a shadow caster moves out of the light's influence and leaves a shadow behind;
		// also need to handle the case where the light is added on the frame the room geom is generated when the shadow map is not yet created;
		// requiring two consecutive frames of no dynamic objects should fix both problems
		bool const cache_shadows(!dynamic_shadows && (i->flags & RO_FLAG_NODYNAM)); // no dynamic object on this frame or the last frame
		if (cache_shadows) {dl_sources.back().assign_smap_id(uintptr_t(&(*i))/sizeof(void *));} // use memory address as a unique ID
		if (dynamic_shadows) {i->flags &= ~RO_FLAG_NODYNAM;} else {i->flags |= RO_FLAG_NODYNAM;}
		dl_sources.back().assign_smap_mgr_id(1); // use a different smap manager than the city (cars + streetlights) so that they don't interfere with each other
		
		if (camera_near_building && lpos.z > camera_bs.z) { // only when the player is near/inside a building and can't see the light bleeding through the floor
			// add a smaller unshadowed light with 360 deg FOV to illuminate the ceiling and other areas as cheap indirect lighting
			// since we're not enabling shadows for this light, it could incorrectly illuminate objects on the floor above, so we must make it small, unless it's on the top floor
			bool const at_top(lpos.z > room.z2() - 0.5f*get_window_vspace());
			float const radius_scale(at_top ? 0.55 : 0.5);
			point const lpos_up(lpos - vector3d(0.0, 0.0, 2.0*i->dz()));
			dl_sources.emplace_back(radius_scale*((room.is_hallway ? 0.3 : room.is_office ? 0.35 : 0.5))*light_radius, lpos_up, lpos_up, color);
			dl_sources.back().set_custom_bcube(clipped_bc); // Note: could reduce clipped_bc further if needed
			dl_sources.back().disable_shadows();
		}
	} // for i
}

bool building_t::toggle_room_light(point const &closest_to) { // Note: called by the player; closest_to is in building space, not camera space
	if (!has_room_geom()) return 0; // error?
	vector<room_object_t> &objs(interior->room_geom->objs);
	auto objs_end(objs.begin() + interior->room_geom->stairs_start); // skip stairs and elevators
	float const window_vspacing(get_window_vspace());
	float closest_dist_sq(0.0);
	unsigned closest_light(0);

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_LIGHT) continue; // not a light
		if (i->z1() < closest_to.z || i->z1() > (closest_to.z + window_vspacing)) continue; // light is on the wrong floor
		float const dist_sq(p2p_dist_sq(closest_to, i->get_cube_center()));
		if (closest_dist_sq == 0.0 || dist_sq < closest_dist_sq) {closest_dist_sq = dist_sq; closest_light = (i - objs.begin());}
	} // for i
	if (closest_dist_sq == 0.0) return 0; // no light found
	assert(closest_light < objs.size());
	objs[closest_light].toggle_lit_state(); // Note: doesn't update indir lighting or room light value
	interior->room_geom->clear_and_recreate_lights(); // recreate light geom with correct emissive properties
	return 1;
}

bool building_t::set_room_light_state_to(room_t const &room, float zval, bool make_on) { // called by AI people
	if (!has_room_geom()) return 0; // error?
	if (room.is_hallway)  return 0; // don't toggle lights for hallways, which can have more than one light
	vector<room_object_t> &objs(interior->room_geom->objs);
	auto objs_end(objs.begin() + interior->room_geom->stairs_start); // skip stairs and elevators
	float const window_vspacing(get_window_vspace());
	bool updated(0);

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_LIGHT) continue; // not a light
		if (i->z1() < zval || i->z1() > (zval + window_vspacing) || !room.contains_cube_xy(*i)) continue; // light is on the wrong floor or in the wrong room
		if (i->is_lit() != make_on) {i->toggle_lit_state(); updated = 1;} // Note: doesn't update indir lighting or room light value
	} // for i
	if (updated) {interior->room_geom->clear_and_recreate_lights();} // recreate light geom with correct emissive properties
	return updated;
}

float room_t::get_light_amt() const { // Note: not normalized to 1.0
	float ext_perim(0.0);

	for (unsigned d = 0; d < 4; ++d) {
		if (ext_sides & (1<<d)) {ext_perim += get_sz_dim(d>>1);} // add length of each exterior side, assuming it has windows
	}
	return ext_perim/get_area_xy(); // light per square meter = exterior perimeter over area
}

