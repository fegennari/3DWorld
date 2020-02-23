// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "lightmap.h" // for light_source
#include "cobj_bsp_tree.h"
#include <atomic>

extern int MESH_Z_SIZE, display_mode;
extern unsigned LOCAL_RAYS, NUM_THREADS;
extern float indir_light_exp;
extern vector<light_source> dl_sources;


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
	if (follow_ray_through_cubes_recur(p1, p2, p1, parts, get_real_parts_end(), parts.end(), cnorm, t)) { // interior ray - find furthest exit point
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

template<typename T> void add_colored_cubes(vector<T> const &cubes, colorRGBA const &color, vect_colored_cube_t &cc) {
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {cc.emplace_back(*c, color);}
}
void building_t::gather_interior_cubes(vect_colored_cube_t &cc) const {

	if (!interior) return; // nothing to do
	building_mat_t const &mat(get_material());
	colorRGBA const wall_color(mat.wall_color.modulate_with(mat.wall_tex.get_avg_color()));
	for (unsigned d = 0; d < 2; ++d) {add_colored_cubes(interior->walls[d], wall_color, cc);}
	add_colored_cubes(interior->elevators, wall_color, cc); // for now elevators are treated the same as walls
	add_colored_cubes(interior->ceilings, mat.ceil_color .modulate_with(mat.ceil_tex .get_avg_color()), cc);
	add_colored_cubes(interior->floors,   mat.floor_color.modulate_with(mat.floor_tex.get_avg_color()), cc);
	add_colored_cubes(details,            detail_color.   modulate_with(mat.roof_tex. get_avg_color()), cc); // should this be included?

	if (has_room_geom()) {
		vector<room_object_t> const &objs(interior->room_geom->objs);
		cc.reserve(cc.size() + objs.size());
		
		for (auto c = objs.begin(); c != objs.end(); ++c) {
			if (c->shape == SHAPE_CYLIN) continue; // cylinders (lights) are not cubes
			cc.emplace_back(*c, c->get_color()); // TODO: to be more accurate, we should use the actual cubes of tables and chairs (which adds a lot of complexity)
			if (c->type == TYPE_TABLE) {cc.back().z1() += 0.85*c->dz();} // at least be a bit more accurate for tables by using only the top
		}
	}
	//cout << TXT(interior->room_geom->objs.size()) << TXT(interior->room_geom->stairs_start) << TXT(cc.size()) << endl;
}

void building_t::ray_cast_room_light(point const &lpos, colorRGBA const &lcolor, cube_bvh_t const &bvh, rand_gen_t &rgen, lmap_manager_t *lmgr, float weight) const {

	// see ray_trace_local_light_source()
	float const tolerance(1.0E-5*bcube.get_max_extent());
	cube_t const scene_bounds(get_scene_bounds_bcube()); // expected by lmap update code
	point const ray_scale(scene_bounds.get_size()/bcube.get_size()), llc_shift(scene_bounds.get_llc() - bcube.get_llc()*ray_scale);

	for (unsigned n = 0; n < LOCAL_RAYS; ++n) {
		vector3d dir(rgen.signed_rand_vector_spherical(1.0).get_norm());
		dir.z = -fabs(dir.z); // make sure dir points down
		point pos(lpos), cpos;
		vector3d cnorm, v_ref;
		colorRGBA cur_color(lcolor), ccolor(WHITE);

		for (unsigned bounce = 0; bounce < 4; ++bounce) { // allow up to 4 bounces
			cpos = pos; // init value
			bool const hit(ray_cast_interior(pos, dir, bvh, cpos, cnorm, ccolor));
			
			if (lmgr != nullptr && cpos != pos) { // accumulate light along the ray from pos to cpos (which is always valid) with color cur_color
				point const p1(pos*ray_scale + llc_shift), p2(cpos*ray_scale + llc_shift); // transform building space to global scene space
				add_path_to_lmcs(lmgr, nullptr, p1, p2, weight, cur_color, LIGHTING_LOCAL, (bounce == 0)); // local light, no bcube
			}
			if (!hit) break; // done
			cur_color = cur_color.modulate_with(ccolor);
			if (cur_color.get_luminance() < 0.05) break; // done
			calc_reflection_angle(dir, v_ref, cnorm);
			v_ref.normalize();
			vector3d const rand_dir(rgen.signed_rand_vector().get_norm());
			dir = (v_ref + rand_dir).get_norm(); // diffuse reflection: new dir is mix 50% specular with 50% random
			if (dot_product(dir, cnorm) < 0.0) {dir.negate();} // make sure it points away from the surface (is this needed?)
			pos = cpos + tolerance*dir; // move slightly away from the surface
		} // for bounce
	} // for n
}

void building_t::ray_cast_building(lmap_manager_t *lmgr, float weight) const {

	if (!has_room_geom()) return; // error?
	timer_t timer("Ray Cast Building");
	vector<room_object_t> &objs(interior->room_geom->objs);
	unsigned const objs_size(interior->room_geom->stairs_start); // skip stairs
	assert(objs_size <= objs.size());
	std::atomic<unsigned> count(0);
	cube_bvh_t bvh;
	gather_interior_cubes(bvh.get_objs());
	bvh.build_tree_top(0); // verbose=0
	weight *= 100.0f/LOCAL_RAYS; // normalize to the number of rays

#pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
	for (int i = 0; i < (int)objs_size; ++i) {
		room_object_t const &ro(objs[i]);
		if (ro.type != TYPE_LIGHT || !ro.is_lit()) continue; // not a light, or light not on
		point lpos(ro.get_cube_center());
		lpos.z = ro.z1() - 0.01*ro.dz(); // set slightly below bottom of light
		rand_gen_t rgen;
		rgen.set_state(i, 123);
		ray_cast_room_light(lpos, ro.get_color(), bvh, rgen, lmgr, weight);
		++count;
	} // for i
	cout << "Lights: " << count << endl;
}

unsigned building_t::create_building_volume_light_texture() const {

	if (!has_room_geom()) return 0; // error?
	lmap_manager_t local_lmap_manager; // cache and reuse this?
	unsigned const zsize(MESH_SIZE[2]); // not MESH_Z_SIZE; we want the clipped size that lmap uses rather than the user-specified size

	if (!local_lmap_manager.is_allocated()) { // first time setup
		unsigned const tot_sz(XY_MULT_SIZE*zsize);
		lmcell init_lmcell;
		local_lmap_manager.alloc(tot_sz, MESH_X_SIZE, MESH_Y_SIZE, zsize, (unsigned char **)nullptr, init_lmcell);
	}
	float const weight(1.0);
	ray_cast_building(&local_lmap_manager, weight);
	return indir_light_tex_from_lmap(local_lmap_manager, MESH_X_SIZE, MESH_Y_SIZE, zsize, indir_light_exp); // indir_light_exp applies to local lighting
}

bool building_t::ray_cast_camera_dir(vector3d const &xlate, point &cpos, colorRGBA &ccolor) const {
	cube_bvh_t bvh;
	gather_interior_cubes(bvh.get_objs());
	bvh.build_tree_top(0); // verbose=0
	vector3d cnorm; // unused
	return ray_cast_interior((get_camera_pos() - xlate), cview_dir, bvh, cpos, cnorm, ccolor);
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
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
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
void building_t::add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building, cube_t &lights_bcube) {

	if (!has_room_geom()) return; // error?
	vector<room_object_t> const &objs(interior->room_geom->objs);
	vect_cube_t &light_bcubes(interior->room_geom->light_bcubes);
	point const camera_bs(camera_pdu.pos - xlate); // camera in building space
	float const window_vspacing(get_window_vspace()), camera_z(camera_bs.z);
	assert(interior->room_geom->stairs_start <= objs.size());
	auto objs_end(objs.begin() + interior->room_geom->stairs_start); // skip stairs
	unsigned camera_part(parts.size()); // start at an invalid value
	bool camera_by_stairs(0), camera_near_building(camera_in_building);

	if (camera_in_building) {
		for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
			if (i->contains_pt(camera_bs)) {camera_part = (i - parts.begin()); break;}
		}
		for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) { // conservative but less efficient
			// Note: stairs that connect stacked parts aren't flagged with has_stairs because stairs only connect to the bottom floor, but they're partially handled below
			if (r->contains_pt(camera_bs)) {camera_by_stairs = r->has_stairs; break;}
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
		float const floor_z(i->z2() - window_vspacing), ceil_z(i->z2());
		bool const floor_is_above(camera_z < floor_z), floor_is_below(camera_z > ceil_z);
		assert(i->room_id < interior->rooms.size());
		room_t const &room(interior->rooms[i->room_id]);
		// less culling if either the light or the camera is by stairs and light is on the floor above or below
		bool const stairs_light((i->has_stairs() || room.has_stairs || camera_by_stairs) && (camera_z > floor_z-window_vspacing) && (camera_z < ceil_z+window_vspacing));
		//if (is_light_occluded(lpos, camera_bs)) continue; // too strong a test in general, but may be useful for selecting high importance lights

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
					assert(camera_part < real_num_parts);
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
		float const light_radius(7.0f*(i->dx() + i->dy())), cull_radius(0.95*light_radius);
		if (!camera_pdu.sphere_visible_test((lpos + xlate), cull_radius)) continue; // VFC
		// check visibility of bcube of light sphere clipped to building bcube; this excludes lights behind the camera and improves shadow map assignment quality
		cube_t clipped_bc; // in building space
		clipped_bc.set_from_sphere(lpos, cull_radius);
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
		if (!camera_pdu.cube_visible(clipped_bc + xlate)) continue; // VFC - post clip
		dl_sources.emplace_back(light_radius, lpos, lpos, color, 0, -plus_z, bwidth); // points down, white for now
		dl_sources.back().set_custom_bcube(clipped_bc);

		if (camera_near_building) { // only when the player is near/inside a building and can't see the light bleeding through the floor
			// add a smaller unshadowed light with 360 deg FOV to illuminate the ceiling and other areas as cheap indirect lighting
			point const lpos_up(lpos - vector3d(0.0, 0.0, 2.0*i->dz()));
			dl_sources.emplace_back(0.5*((room.is_hallway ? 0.3 : room.is_office ? 0.3 : 0.5))*light_radius, lpos_up, lpos_up, color);
			dl_sources.back().set_custom_bcube(clipped_bc); // Note: could reduce clipped_bc further if needed
			dl_sources.back().disable_shadows();
		}
	} // for i
}

float room_t::get_light_amt() const { // Note: not normalized to 1.0
	float ext_perim(0.0);

	for (unsigned d = 0; d < 4; ++d) {
		if (ext_sides & (1<<d)) {ext_perim += get_sz_dim(d>>1);} // add length of each exterior side, assuming it has windows
	}
	return ext_perim/get_area_xy(); // light per square meter = exterior perimeter over area
}

