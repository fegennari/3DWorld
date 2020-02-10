// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "lightmap.h" // for light_source

extern vector<light_source> dl_sources;


bool ray_cast_cube(point const &p1, point const &p2, cube_t const &c, vector3d &cnorm, float &t) {
	float tmin(0.0), tmax(1.0);
	if (!get_line_clip(p1, p2, c.d, tmin, tmax) || tmin >= t) return 0;
	t = tmin;
	get_closest_cube_norm(c.d, (p1 + (p2 - p1)*t), cnorm);
	return 1;
}
bool ray_cast_vect_cube(point const &p1, point const &p2, vect_cube_t const &cubes, vector3d &cnorm, float &t) {
	bool ret(0);
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {ret |= ray_cast_cube(p1, p2, *c, cnorm, t);}
	return ret;
}

bool building_t::ray_cast_interior(point const &pos, vector3d const &dir, point &cpos, vector3d &cnorm, colorRGBA &ccolor) const { // pos in building space
	if (!interior || is_rotated() || !is_simple_cube()) return 0; // these cases are not yet supported
	float const extent(bcube.get_max_extent());
	cube_t clip_cube(bcube);
	clip_cube.expand_by(0.01*extent); // expand slightly so that collisions with objects on the edge are still considered interior
	point p1(pos), p2(pos + dir*(2.0*extent));
	if (!do_line_clip(p1, p2, clip_cube.d)) return 0; // ray does not intersect building bcube
	building_mat_t const &mat(get_material());
	float t(1.0); // at p2
	bool hit_side(0);

	// check parts (exterior walls)
	for (auto p = parts.begin(); p != parts.end(); ++p) { // should chimneys and porch roofs be included?
		if (p->contains_pt(p1)) { // interior ray - find exit point
			float tmin(0.0), tmax(1.0);
			if (!get_line_clip(p1, p2, p->d, tmin, tmax) || tmax >= t) continue;
			t = tmax;
			get_closest_cube_norm(p->d, (p1 + (p2 - p1)*t), cnorm);
			cnorm.negate(); // reverse hit dir
		}
		else { // exterior ray - find entrance point
			if (!ray_cast_cube(p1, p2, *p, cnorm, t)) continue;
		}
		hit_side = 1;
	} // for p
	if (hit_side) {ccolor = side_color.modulate_with(mat.side_tex.get_avg_color());}

	for (auto r = roof_tquads.begin(); r != roof_tquads.end(); ++r) {
		// WRITE; use roof_color/mat.roof_tex
	}
	// check walls, floors, and ceilings
	bool const hit_wall(ray_cast_vect_cube(p1, p2, interior->walls[0], cnorm, t) || ray_cast_vect_cube(p1, p2, interior->walls[1], cnorm, t));
	if (hit_wall) {ccolor = mat.wall_color.modulate_with(mat.wall_tex.get_avg_color());}
	if (ray_cast_vect_cube(p1, p2, interior->ceilings, cnorm, t)) {ccolor = mat.ceil_color .modulate_with(mat.ceil_tex .get_avg_color());}
	if (ray_cast_vect_cube(p1, p2, interior->floors,   cnorm, t)) {ccolor = mat.floor_color.modulate_with(mat.floor_tex.get_avg_color());}
	if (ray_cast_vect_cube(p1, p2, details,            cnorm, t)) {ccolor = detail_color.   modulate_with(mat.roof_tex. get_avg_color());} // should this be included?

	if (has_room_geom()) {
		vector<room_object_t> const &objs(interior->room_geom->objs);
		point const cur_p2(p1 + (p2 - p1)*t); // clip to reduce the chance of a stairwell intersection
		bool hit_stairs(0);

		for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
			if (s->line_intersects(p1, cur_p2)) {hit_stairs = 1; break;}
		}
		auto objs_end(hit_stairs ? objs.end() : (objs.begin() + interior->room_geom->stairs_start)); // stairs optimization

		for (auto c = objs.begin(); c != objs.end(); ++c) {
			if (ray_cast_cube(p1, p2, *c, cnorm, t)) {ccolor = c->get_color();}
		}
	}
	if (t == 1.0) return 0; // no intersection
	cpos = p1 + (p2 - p1)*t;
	return 1;
}

void building_t::add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building, cube_t &lights_bcube) const {

	if (!has_room_geom()) return; // error?
	vector<room_object_t> &objs(interior->room_geom->objs);
	point const camera_bs(camera_pdu.pos - xlate); // camera in building space
	float const window_vspacing(get_window_vspace()), camera_z(camera_bs.z);
	assert(interior->room_geom->stairs_start <= objs.size());
	auto objs_end(objs.begin() + interior->room_geom->stairs_start); // skip stairs
	unsigned camera_part(parts.size()); // start at an invalid value
	bool camera_by_stairs(0);

	if (camera_in_building) {
		for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
			if (i->contains_pt(camera_bs)) {camera_part = (i - parts.begin()); break;}
		}
		for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) { // conservative but less efficient
			// Note: stairs that connect stacked parts aren't flagged with has_stairs because stairs only connect to the bottom floor, but they're partially handled below
			if (r->contains_pt(camera_bs)) {camera_by_stairs = r->has_stairs; break;}
		}
	}
	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_LIGHT || !i->is_lit()) continue; // not a light, or light not on
		point const lpos(i->get_cube_center()); // centered in the light fixture
		if (!lights_bcube.contains_pt_xy(lpos)) continue; // not contained within the light volume
		float const floor_z(i->z2() - window_vspacing), ceil_z(i->z2());
		bool const floor_is_above(camera_z < floor_z), floor_is_below(camera_z > ceil_z);
		// less culling if either the light or the camera is by stairs and light is on the floor above or below
		bool const stairs_light((i->has_stairs() || camera_by_stairs) && (camera_z > floor_z-window_vspacing) && (camera_z < ceil_z+window_vspacing));
		assert(i->room_id < interior->rooms.size());
		room_t const &room(interior->rooms[i->room_id]);

		if (floor_is_above || floor_is_below) { // light is on a different floor from the camera
			if (camera_in_building && room.part_id == camera_part ||
				(room.contains_pt_xy(camera_bs) && camera_z < ceil_z+window_vspacing && camera_z > floor_z-window_vspacing))
			{
				// player is on a different floor of the same building part, or more than one floor away in a part stack, and can't see a light from the floor above/below
				if (!stairs_light) continue; // camera in building and on wrong floor, don't add light
				if (camera_z < (i->z2() - 2.0*window_vspacing) || camera_z > (i->z2() + window_vspacing)) continue; // light is on the stairs, add if one floor above/below
			}
			else { // camera outside the building (or the part that contains this light)
				float const xy_dist(p2p_dist_xy(camera_bs, lpos));
				if (!stairs_light && ((camera_z - lpos.z) > 2.0f*xy_dist || (lpos.z - camera_z) > 0.5f*xy_dist)) continue; // light viewed at too high an angle

				if (camera_in_building) { // camera and light are in different buildings/parts
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
		colorRGBA color;
		if (is_house) {color = colorRGBA(1.0, 1.0, 0.8);} // house - yellowish
		else if (room.is_hallway || room.is_office) {color = colorRGBA(0.8, 0.8, 1.0);} // office building - blueish
		else {color = colorRGBA(1.0, 1.0, 1.0);} // white - small office
		dl_sources.emplace_back(light_radius, lpos, lpos, color, 0, -plus_z, bwidth); // points down, white for now
		dl_sources.back().set_building_id(building_id);

		if (camera_in_building) { // only when the player is inside a building and can't see the light bleeding through the floor
			// add a smaller unshadowed light with 360 deg FOV to illuminate the ceiling and other areas as cheap indirect lighting
			point const lpos_up(lpos - vector3d(0.0, 0.0, 2.0*i->dz()));
			dl_sources.emplace_back(0.5*((room.is_hallway ? 0.3 : room.is_office ? 0.3 : 0.5))*light_radius, lpos_up, lpos_up, color);
			dl_sources.back().set_building_id(building_id);
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

