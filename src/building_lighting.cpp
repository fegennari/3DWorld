// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "lightmap.h" // for light_source
#include "cobj_bsp_tree.h"
#include "profiler.h"
#include <thread>

bool  const USE_BKG_THREAD      = 1;
bool  const INDIR_BASEMENT_EN   = 1;
bool  const INDIR_ATTIC_ENABLE  = 1;
bool  const INDIR_BLDG_ENABLE   = 1;
unsigned INDIR_LIGHT_FLOOR_SPAN = 5; // in number of floors, generally an odd number to represent current floor and floors above/below; 0 is unlimited
float const ATTIC_LIGHT_RADIUS_SCALE = 2.0; // larger radius in attic, since space is larger

vector<point> enabled_bldg_lights;

extern bool camera_in_building, player_in_walkway, some_person_has_idle_animation;
extern int MESH_Z_SIZE, display_mode, display_framerate, camera_surf_collide, animate2, frame_counter, building_action_key, player_in_basement, player_in_elevator, player_in_attic;
extern unsigned LOCAL_RAYS, MAX_RAY_BOUNCES, NUM_THREADS;
extern float indir_light_exp, fticks;
extern double tfticks;
extern colorRGB cur_ambient, cur_diffuse;
extern std::string lighting_update_text;
extern vector<light_source> dl_sources;
extern building_t const *player_building;

bool enable_building_people_ai();
bool check_cube_occluded(cube_t const &cube, vect_cube_t const &occluders, point const &viewer);
bool cube_visible_in_building_mirror_reflection(cube_t const &c);
void calc_cur_ambient_diffuse();
colorRGBA get_textured_wood_color();
bool bed_has_canopy_mat(room_object_t const &c);
int get_canopy_texture();
colorRGBA get_canopy_base_color(room_object_t const &c);
void get_water_heater_cubes(room_object_t const &wh, cube_t cubes[2]);
bool line_int_polygon_sides(point const &p1, point const &p2, cube_t const &bcube, vect_point const &points, float &t);
void get_pool_table_cubes(room_object_t const &c, cube_t cubes[5]);
unsigned get_couch_cubes(room_object_t const &c, cube_t cubes[4]);
unsigned get_cashreg_cubes(room_object_t const &c, cube_t cubes[2]);;

bool check_indir_enabled(bool in_basement, bool in_attic) {
	if (in_basement) return INDIR_BASEMENT_EN;
	if (in_attic   ) return INDIR_ATTIC_ENABLE;
	return INDIR_BLDG_ENABLE; // not basement or attic
}
bool enable_building_indir_lighting_no_cib() {
	if (!(display_mode & 0x10)) return 0; // key 5
	if (!check_indir_enabled(player_in_basement, player_in_attic)) return 0;
	if (MESH_SIZE[2] == 0) return 0; // no volume texture allocated
	return 1;
}
bool enable_building_indir_lighting() {return (camera_in_building && enable_building_indir_lighting_no_cib());}

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
	vect_tquad_with_ix_t roof_tquads; // these aren't cubes, so they don't go into the BVH; should only be 2-4 of these

	vect_colored_cube_t &get_objs() {return objects;}

	void clear() {
		cobj_tree_simple_type_t<colored_cube_t>::clear();
		roof_tquads.clear();
	}
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
	vect_cube_t::const_iterator end, vect_cube_t::const_iterator parent, unsigned depth, vector3d &cnorm, float &t)
{
	if (depth > 32) return 0; // occasional ray trace query can cause a stack overflow, likely due to FP error, so limit recursion to 32

	for (auto c = cubes.begin(); c != end; ++c) {
		if (c == parent || !c->contains_pt(start)) continue;
		float tmin(0.0), tmax(1.0);
		if (!get_line_clip(p1, p2, c->d, tmin, tmax) || tmax >= t) continue;
		point const cpos(p1 + (p2 - p1)*(tmax + 0.001)); // move slightly beyond the hit point
		if ((end - cubes.begin()) > 1 && follow_ray_through_cubes_recur(p1, p2, cpos, cubes, end, c, depth+1, cnorm, t)) return 1;
		t = tmax;
		get_closest_cube_norm(c->d, (p1 + (p2 - p1)*t), cnorm); // this is the final point, update cnorm
		cnorm.negate(); // reverse hit dir
		return 1;
	} // for c
	return 0;
}
bool building_t::ray_cast_exterior_walls(point const &p1, point const &p2, vector3d &cnorm, float &t) const {
	return follow_ray_through_cubes_recur(p1, p2, p1, parts, get_real_parts_end_inc_sec(), parts.end(), 0, cnorm, t);
}

colorRGBA building_interior_t::get_attic_ceiling_color() const {
	switch (attic_type) {
	case ATTIC_TYPE_RAFTERS   : return texture_color(FENCE_TEX );
	case ATTIC_TYPE_WOOD      : return texture_color(get_plywood_tid());
	case ATTIC_TYPE_PLASTER   : return texture_color(STUCCO_TEX);
	case ATTIC_TYPE_FIBERGLASS: return texture_color(get_insulation_tid()).modulate_with(colorRGBA(1.0, 0.7, 0.6));
	default: assert(0);
	}
	return WHITE; // never gets here
}

bool line_contained_in_cube(point const &p1, point const &p2, vect_cube_t const &cubes, vect_cube_t::const_iterator end) {
	for (auto i = cubes.begin(); i != end; ++i) {
		if (i->contains_pt(p1) && i->contains_pt(p2)) return 1;
	}
	return 0;
}

void building_t::set_building_colors(building_colors_t &bcolors) const {
	building_mat_t const &mat(get_material());
	bcolors.side_color = side_color.modulate_with(mat.side_tex.get_avg_color());
	bcolors.wall_color = mat.wall_color.modulate_with(mat.wall_tex.get_avg_color());
	if (has_basement()) {bcolors.basement_wall_color = get_basement_wall_texture().get_avg_color();} // this one is particularly expensive
	if (has_attic   ()) {bcolors.attic_color         = interior->get_attic_ceiling_color();}
}

// Note: static objects only; excludes people; pos in building space
bool building_t::ray_cast_interior(point const &pos, vector3d const &dir, cube_t const &valid_area, cube_bvh_t const &bvh, bool in_attic, bool in_ext_basement,
	building_colors_t const &bcolors, point &cpos, vector3d &cnorm, colorRGBA &ccolor, rand_gen_t *rgen) const
{
	if (!interior || is_rotated()) return 0; // these cases are not yet supported
	float const extent(valid_area.get_max_dim_sz());
	cube_t clip_cube(valid_area);
	clip_cube.expand_by(0.01*extent); // expand slightly so that collisions with objects on the edge are still considered interior
	point p1(pos), p2(pos + dir*(2.0*extent));
	if (!do_line_clip(p1, p2, clip_cube.d)) return 0; // ray does not intersect clip cube
	float t(1.0); // start at p2
	bool hit(0);
	cpos = p2; // use far clip point for clip cube if there is no hit

	if (in_attic) {
		for (tquad_with_ix_t const &tq : bvh.roof_tquads) {
			vector3d const normal(-tq.get_norm()); // negate because we're testing the inside of the roof
			float t_tq(t);
			if (!line_poly_intersect(p1, p2, tq.pts, tq.npts, normal, t_tq) || t_tq >= t) continue;
			t = t_tq; cnorm = normal; hit = 1; ccolor = bcolors.attic_color;
		}
		if (hit) {p2  = p1 + (p2 - p1)*t; t = 1.0;} // clip p2 to t (minor optimization)
	}
	else if (in_ext_basement) {} // no exterior walls to check
	else if (line_contained_in_cube(p1, p2, parts, get_real_parts_end_inc_sec())) {} // both points in same part; no exterior wall hit
	// check parts (exterior walls); should chimneys and porch roofs be included?
	else if (ray_cast_exterior_walls(p1, p2, cnorm, t)) { // interior ray (common case) - find furthest exit point
		p2  = p1 + (p2 - p1)*t; t = 1.0; // clip p2 to t (minor optimization)
		hit = 1;
		if (p2.z < ground_floor_z1) {ccolor = bcolors.basement_wall_color;} // basement wall
		else {ccolor = bcolors.wall_color;} // non-basement wall
	}
	else { // check for exterior rays (uncommon case)
		bool hit(0);

		if (!is_cube()) {
			int const part_ix(get_part_ix_containing_pt(p2));
			if (part_ix >= 0 && line_int_polygon_sides(p1, p2, parts[part_ix], get_part_ext_verts(part_ix), t)) {hit = 1;}
		}
		if (!hit) {
			auto const parts_end(get_real_parts_end_inc_sec());
			for (auto p = parts.begin(); p != parts_end; ++p) {hit |= ray_cast_cube(p1, p2, *p, cnorm, t);} // find closest entrance point
		}
		if (hit) { // exterior hit (ray outside building) - don't need to check interior geometry
			cpos   = p1 + (p2 - p1)*t;
			ccolor = bcolors.side_color;
			return 1;
		}
	}
	bvh.ray_cast(p1, p2, cnorm, ccolor, t);

	if (t == 1.0) { // no intersection with bvh
		cpos = p2;
		if (!hit) return 0;
		if (rgen && p2.z > ground_floor_z1 && has_int_windows() && rgen->rand_bool()) return 0; // 50% chance of exiting through a window
		return 1;
	}
	cpos = p1 + (p2 - p1)*t;
	return 1;
}

unsigned elevator_t::get_coll_cubes(cube_t cubes[5]) const {
	//if (!is_open) {cubes[0] = *this; return 1;} // closed, entire elevator is 1 cube
	float const wwidth(get_wall_thickness()), fwidth(get_frame_width());
	for (unsigned i = 0; i < 5; ++i) {cubes[i] = *this;} // start with full cube
	cubes[0].d[dim][dir] = cubes[0].d[dim][!dir] + (dir ? 1.0 : -1.0)*wwidth; // back
	
	for (unsigned d = 0; d < 2; ++d) {
		cube_t &side(cubes[d+1]);
		side.d[!dim][d] = side.d[!dim][!d] + (d ? 1.0 : -1.0)*wwidth; // sides
	}
	if (open_amt > 0.0) { // open at least partially
		for (unsigned d = 0; d < 2; ++d) {
			cube_t &front(cubes[d+3]);
			front.d[ dim][!dir] = front.d[dim][dir] + (dir ? -1.0 : 1.0)*wwidth; // front side parts
			front.d[!dim][d]    = front.d[!dim][!d] + (d   ? 1.0 : -1.0)*fwidth;
		}
		return 5; // back + 2 sides + 2 front parts
	}
	else { // closed
		cube_t &front(cubes[3]);
		front.d[dim][!dir] = front.d[dim][dir] + (dir ? -1.0 : 1.0)*wwidth; // front (sides + door)
		return 4; // back + 2 sides + front
	}
}

template<typename T> void add_colored_cubes(vector<T> const &cubes, colorRGBA const &color, cube_t const &ext_bcube, vect_colored_cube_t &cc) {
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (c->intersects(ext_bcube)) {cc.emplace_back(*c, color);}
	}
}
void add_colored_cubes(cube_t const *const cubes, unsigned num_cubes, colorRGBA const &color, vect_colored_cube_t &cc) { // Note: no ext_bcube intersection check
	for (unsigned n = 0; n < num_cubes; ++n) {
		if (!cubes[n].is_all_zeros()) {cc.emplace_back(cubes[n], color);}
	}
}
void building_t::gather_interior_cubes(vect_colored_cube_t &cc, cube_t const &ext_bcube) const {
	if (!interior) return; // nothing to do
	building_mat_t const &mat(get_material());
	colorRGBA const wall_color(wall_color.modulate_with(mat.wall_tex.get_avg_color()));
	float const floor_spacing(get_window_vspace()), z1(ext_bcube.z1()), z2(ext_bcube.z2());
	float const stairs_z1(z1 - floor_spacing), stairs_z2(z2 + floor_spacing); // stairs extend an extra floor up and down to block rays in stairwells
	for (unsigned d = 0; d < 2; ++d) {add_colored_cubes(interior->walls[d], wall_color, ext_bcube, cc);}

	for (elevator_t const &e : interior->elevators) {
		if (!e.intersects(ext_bcube)) continue;
		cube_t cubes[5];
		unsigned const num_cubes(e.get_coll_cubes(cubes));
		// for now elevators are treated the same as walls with the same color, even though the inside of open elevators is wood
		add_colored_cubes(cubes, num_cubes, wall_color, cc); // can only assign the same color to all sides of the cube
	}
	for (escalator_t const &e : interior->escalators) {
		cube_t cubes[7];
		e.get_all_cubes(cubes);
		add_colored_cubes(cubes, 6, LT_GRAY, cc); // walls and floors are light gray
		cc.emplace_back(cubes[6], WHITE); // pillar is white
		//e.get_ramp_bcube(0) - ignore ramp for now because it's not a cube
	}
	for (door_t const &d : interior->doors) {
		if (d.open || !d.intersects(ext_bcube)) continue; // add only closed doors
		cc.emplace_back(d.get_true_bcube(), WHITE);
	}
	for (cube_t const &i : interior->floors) {
		if (!i.intersects(ext_bcube)) continue;
		tid_nm_pair_t tex;
		cc.emplace_back(i, get_floor_tex_and_color(i, tex).modulate_with(tex.get_avg_color()));
	}
	for (cube_t const &i : interior->ceilings) {
		if (!i.intersects(ext_bcube)) continue;
		tid_nm_pair_t tex;
		cc.emplace_back(i, get_ceil_tex_and_color(i, tex).modulate_with(tex.get_avg_color()));
	}
	// should glass floors be included? maybe not, since we want refractions rather than reflections; plus they're mostly transparent and near white
	//for (cube_t const &i : interior->room_geom->glass_floors) {if (i.intersects(ext_bcube)) {cc.emplace_back(i, GLASS_COLOR);}}
	add_colored_cubes(details, detail_color.modulate_with(mat.roof_tex.get_avg_color()), ext_bcube, cc); // should this be included?
	if (!has_room_geom()) return; // nothing else to add
	vect_room_object_t const &objs(interior->room_geom->objs);
	static vect_cube_t temp; // used across calls for subtracting holes
		
	for (auto c = objs.begin(); c != objs.end(); ++c) { // Note: ignores expanded objects (including shelf rack objects)
		room_object const type(c->type);
		if (!c->is_visible())         continue;
		if (type  == TYPE_ELEVATOR)   continue; // elevator cars/internals can move so should not contribute to lighting
		if (type  == TYPE_SHOWER  )   continue; // transparent
		if (type  == TYPE_BLOCKER || type == TYPE_COLLIDER) continue; // blockers and colliders are not drawn
		if (c->shape == SHAPE_ANGLED) continue; // sloped objects such as parking garage and pool ramps aren't cubes
		if (c->is_exterior())         continue; // interior objects only
		// skip other object types that are too small, not cube shaped, or not interior
		if (type == TYPE_WALL_TRIM || type == TYPE_BOOK || type == TYPE_CRACK || type == TYPE_PLANT || type == TYPE_RAILING || type == TYPE_SHELVES || type == TYPE_BOTTLE ||
			type == TYPE_PEN || type == TYPE_PENCIL || is_ball_type(type) || type == TYPE_HANGER_ROD || type == TYPE_DRAIN || type == TYPE_MONEY || type == TYPE_PHONE ||
			type == TYPE_TPROLL || type == TYPE_SPRAYCAN || type == TYPE_MARKER || type == TYPE_BUTTON || type == TYPE_SWITCH || type == TYPE_TAPE || type == TYPE_OUTLET ||
			type == TYPE_PARK_SPACE || type == TYPE_RAMP || type == TYPE_PIPE || type == TYPE_VENT || type == TYPE_BREAKER || type == TYPE_KEY || type == TYPE_HANGER ||
			type == TYPE_FESCAPE || type == TYPE_CUP || type == TYPE_CLOTHES || type == TYPE_LAMP || type == TYPE_OFF_CHAIR || type == TYPE_LIGHT || type == TYPE_SIGN ||
			type == TYPE_PAPER || type == TYPE_WALL_LAMP || type == TYPE_RCHAIR || type == TYPE_SILVER || type == TYPE_STAPLER || type == TYPE_WIND_SILL || type == TYPE_BALCONY ||
			type == TYPE_TOY_MODEL || type == TYPE_CEIL_FAN || type == TYPE_PLANT_MODEL || type == TYPE_POOL_FLOAT || type == TYPE_BENCH || type == TYPE_DIV_BOARD ||
			type == TYPE_POOL_LAD || type == TYPE_FLASHLIGHT || type == TYPE_CANDLE || type == TYPE_CAMERA || type == TYPE_CLOCK || type == TYPE_BAR_STOOL || type == TYPE_PADLOCK ||
			type == TYPE_WFOUNTAIN || type == TYPE_BANANA || type == TYPE_BAN_PEEL || type == TYPE_VALVE) continue;
		bool const is_stairs(type == TYPE_STAIR || type == TYPE_STAIR_WALL);
		if (c->z1() > (is_stairs ? stairs_z2 : z2) || c->z2() < (is_stairs ? stairs_z1 : z1)) continue;
		if (!c->intersects_xy(ext_bcube)) continue;
		colorRGBA const color(c->get_color());
		
		if (c->is_round()) {
			cube_t inner_cube(*c);

			if (type == TYPE_WHEATER) { // {tank, pipes}
				cube_t cubes[2];
				get_water_heater_cubes(*c, cubes);
				cc.emplace_back(cubes[1], color); // add pipes directly
				inner_cube = cubes[0]; // tank
			}
			float const shrink(c->get_radius()*(1.0 - 1.0/SQRT2));
			inner_cube.expand_by_xy(-shrink); // shrink to inscribed cube in XY
			if (c->shape == SHAPE_SPHERE) {inner_cube.expand_in_dim(2, -shrink);} // shrink in Z as well
			if (   type  == TYPE_TABLE  ) {inner_cube.z1() += 0.88*c->dz();} // top of table
			cc.emplace_back(inner_cube, color);
		}
		else if (type == TYPE_CLOSET) { // Note: lighting cubes and indir lighting are *not* updated when closet doors are opened and closed
			cube_t cubes[5];
			get_closet_cubes(*c, cubes, 1); // for_collision=1
			if (!c->is_small_closet() && c->is_open()) {cubes[4].set_to_zeros();} // ignore open doors of large closets
			add_colored_cubes(cubes, 5, color, cc); // include door, whether closed or open
		}
		else if (type == TYPE_BED) { // Note: posts are not included
			colorRGBA const wood_color(get_textured_wood_color()), sheets_color(c->color.modulate_with(texture_color(c->get_sheet_tid())));
			cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
			get_bed_cubes(*c, cubes);
			add_colored_cubes(cubes,   3, wood_color,   cc); // frame, head, foot
			add_colored_cubes(cubes+3, 2, sheets_color, cc); // mattress, pillow

			if (bed_has_canopy_mat(*c)) { // add canopy (but not posts); only the top surface, not the corner triangles
				cube_t canopy(*c);
				canopy.z1() = cubes[2].z2() + 1.4*c->dz();
				canopy.z2() = canopy.z1() + 0.001*c->dz(); // almost zero thickness
				cc.emplace_back(canopy, get_canopy_base_color(*c).modulate_with(texture_color(get_canopy_texture())));
			}
			get_tc_leg_cubes(cubes[5], *c, BED_HEAD_WIDTH, 0, cubes); // cubes[5] is not overwritten
			add_colored_cubes(cubes, 4, wood_color, cc); // legs
		}
		else if (type == TYPE_DESK || type == TYPE_DRESSER || type == TYPE_NIGHTSTAND || type == TYPE_TABLE) { // objects with legs
			if (c->is_glass_table()) continue; // skip glass table (transparent with thin legs)
			cube_t cubes[7];
			unsigned const num(get_table_like_object_cubes(*c, cubes));
			add_colored_cubes(cubes, num, color, cc);
		}
		else if (type == TYPE_RDESK) {
			cube_t cubes[3];
			get_reception_desk_cubes(*c, cubes);
			add_colored_cubes(cubes, 3, color, cc);
		}
		else if (type == TYPE_CHAIR) {
			colorRGBA const wood_color(get_textured_wood_color());
			cube_t cubes[3], leg_cubes[4]; // seat, back, legs_bcube
			get_chair_cubes(*c, cubes);
			cc.emplace_back(cubes[0], color     ); // seat
			cc.emplace_back(cubes[1], wood_color); // back
			get_tc_leg_cubes(cubes[2], *c, CHAIR_LEG_WIDTH, 1, leg_cubes);
			add_colored_cubes(leg_cubes, 4, wood_color, cc);
		}
		else if (type == TYPE_CUBICLE || (type == TYPE_STALL && c->shape != SHAPE_SHORT)) { // cubicle or bathroom stall - hollow
			bool const is_stall(type != TYPE_CUBICLE);
			cube_t sides(*c);
			float const dz(c->dz());
			
			if (is_stall) { // set halfway between sides and door height for simplicity; cubicle ls already the correct height
				sides.z2() -= 0.365*dz;
				sides.z1() += 0.165*dz;
			}
			cube_t inside(sides);
			inside.expand_by_xy(-(is_stall ? 0.0125 : 0.07)*dz); // shrink by wall thickness
			temp.clear();
			subtract_cube_from_cube(sides, inside, temp);
			assert(temp.size() == 4); // -y, +y, -x, +x
			unsigned const front_ix(3 - (2*c->dim + c->dir)); // dim|dir:front_ix: 00:3, 01:2, 10:1, 11:0

			for (unsigned n = 0; n < 4; ++n) { // front at dim,!dir
				if ((!is_stall || c->is_open()) && n == front_ix) continue; // open front, skip
				cc.emplace_back(temp[n], color);
			}
		}
		else if (type == TYPE_POOL_TABLE) {
			cube_t const top(get_pool_table_top_surface(*c));
			cube_t cubes[5]; // body + 4 legs
			get_pool_table_cubes(*c, cubes);
			cubes[0].z2() = top.z1(); // clip top of body to bottom of top
			for (unsigned n = 0; n < 5; ++n) {cc.emplace_back(cubes[n], BROWN);} // body and legs are brown
			cc.emplace_back(top, GREEN); // top surface is green
		}
		else if (type == TYPE_SHELFRACK) {
			cube_t back, top, sides[2], shelves[5];
			unsigned const num_shelves(get_shelf_rack_cubes(*c, back, top, sides, shelves));
			cc.emplace_back(back, color*0.67);
			add_colored_cubes(shelves, num_shelves, color, cc);
			if (!top     .is_all_zeros()) {cc.emplace_back(top, color);}
			if (!sides[0].is_all_zeros()) {add_colored_cubes(sides, 2, color, cc);}
		}
		else if (type == TYPE_COUCH) {
			cube_t cubes[4]; // bottom, back, arm, arm
			unsigned const num_cubes(get_couch_cubes(*c, cubes));
			for (unsigned n = 0; n < num_cubes; ++n) {cc.emplace_back(cubes[n], color);}
		}
		else if (type == TYPE_CASHREG) {
			cube_t cubes[2]; // body, screen
			get_cashreg_cubes(*c, cubes);
			cc.emplace_back(cubes[0], color); // only add the body
		}
		else { // single cube
			cube_t bc(*c); // handle 3D models that don't fill the entire cube
			bool const dim(c->dim), dir(c->dir);

			if (type == TYPE_TUB) {
				bc.z2() -= 0.9*bc.dz(); // bottom
				cube_t top(*c);
				top.z1() = bc.z2();
				cube_t inside(top);
				inside.expand_by_xy(-0.1*c->get_sz_dim(dim)); // shrink
				temp.clear();
				subtract_cube_from_cube(top, inside, temp);
				assert(temp.size() == 4); // -y, +y, -x, +x
				add_colored_cubes(temp, color, ext_bcube, cc);
			}
			else if (type == TYPE_HOOD) {
				bc.expand_in_dim(!dim, 0.25*c->get_sz_dim(!dim)); // shrink width
				bc.d[dim][dir] -= (dir ? 1.0 : -1.0)*0.25*c->get_sz_dim(dim); // shift front in
				bc.z2() += 0.25*bc.dz(); // shift bottom up
			}
			else if (type == TYPE_STOVE ) {bc.z2() -= 0.20*bc.dz();}
			else if (type == TYPE_TOILET) {bc.z2() -= 0.33*bc.dz();}
			else if (type == TYPE_SINK  ) {bc.z2() -= 0.20*bc.dz(); bc.z1() += 0.65*bc.dz();}
			else if (type == TYPE_MONITOR || type == TYPE_TV) {bc.expand_in_dim(dim, -0.3*bc.get_sz_dim(dim));} // reduce thickness
			else if (type == TYPE_BRSINK) {bc.z1() += 0.60*bc.dz();}
			else if (type == TYPE_ATTIC_DOOR) {bc = get_attic_access_door_cube(*c, 0);} // inc_ladder=0: includes door but not ladder
			else if (type == TYPE_SHOWERTUB ) {bc = get_shower_tub_wall(*c);}
			cc.emplace_back(bc, color);
		}
	} // for c
}


colorRGBA get_outdoor_light_color() {
	calc_cur_ambient_diffuse(); // needed for correct outdoor color
	return (cur_ambient + cur_diffuse); // sum of each
}

unsigned const IS_WINDOW_BIT = (1<<24); // if this bit is set, the light is from a window; if not, it's from a light room object

class building_indir_light_mgr_t {
	bool is_running, kill_thread, lighting_updated, needs_to_join, need_bvh_rebuild, update_windows, is_negative_light, in_ext_basement;
	int cur_bix, cur_light, cur_floor;
	unsigned cur_tid;
	colorRGBA outdoor_color;
	cube_t valid_area, light_bounds;
	vector<unsigned char> tex_data;
	vector<unsigned> light_ids;
	vector<pair<float, unsigned>> lights_to_sort;
	deque<unsigned> remove_queue;
	set<unsigned> lights_complete, lights_seen;
	vect_cube_with_ix_t windows;
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
			rt_thread = std::thread(&building_indir_light_mgr_t::cast_light_rays, this, b);
			needs_to_join = 1;
		}
		else {
			// per-light time for large office building: orig: 194ms, per-floor BVH: 96ms, clip rays to floor: 44ms, now 37ms
			highres_timer_t timer("Ray Cast Building Light");
			cast_light_rays(b);
		}
	}
	vector3d get_reflect_dir(vector3d const &dir, vector3d const &cnorm) {
		vector3d v_ref;
		calc_reflection_angle(dir, v_ref, cnorm);
		v_ref.normalize();
		return v_ref;
	}
	void calc_reflect_ray(point &pos, point const &cpos, vector3d &dir, vector3d const &cnorm, vector3d const &v_ref, rand_gen_t &rgen, float tolerance) const {
		vector3d const rand_dir(rgen.signed_rand_vector().get_norm());
		dir = (v_ref + rand_dir).get_norm(); // diffuse reflection: new dir is mix 50% specular with 50% random
		if (dot_product(dir, cnorm) < 0.0) {dir.negate();} // make sure it points away from the surface (is this needed?)
		pos = cpos + tolerance*dir; // move slightly away from the surface
	}
	void cast_light_rays(building_t const &b) {
		// Note: modifies lmgr, but otherwise thread safe
		unsigned const num_rt_threads(max(1U, (NUM_THREADS - (USE_BKG_THREAD ? 1 : 0)))); // reserve a thread for the main thread if running in the background
		unsigned base_num_rays(LOCAL_RAYS), dim(2), dir(0); // default dim is z; dir=2 is omnidirectional
		cube_t const scene_bounds(get_scene_bounds_bcube()); // expected by lmap update code
		vector3d const ray_scale(scene_bounds.get_size()/light_bounds.get_size()), llc_shift(scene_bounds.get_llc() - light_bounds.get_llc()*ray_scale);
		float const tolerance(1.0E-5*valid_area.get_max_dim_sz());
		bool const is_window(cur_light & IS_WINDOW_BIT);
		bool in_attic(0), in_ext_basement(0), is_skylight(0);
		float weight(100.0), light_radius(0.0);
		point light_center;
		cube_t light_cube;
		colorRGBA lcolor, pri_lcolor;
		vector3d light_dir; // points toward the light
		assert(cur_light >= 0);

		if (is_window) { // window
			unsigned const window_ix(cur_light & ~IS_WINDOW_BIT);
			assert(window_ix < windows.size());
			cube_with_ix_t const &window(windows[window_ix]);
			float surface_area(0.0);
			light_cube = window;

			if (window.dz() < min(window.dx(), window.dy())) { // skylight; we could encode skylights as a different ix, but testing aspect ratio is easier
				is_skylight    = 1;
				surface_area   = window.dx()*window.dy();
				base_num_rays *= 8; // more rays, since skylights are larger and can cover multiple rooms
				weight        *= 10.0; // stronger due to direct sun/moon/cloud lighting and reduced occlusion from buildings and terrain
				light_cube.translate_dim(2, -b.get_fc_thickness()); // shift slightly down into the building to avoid collision with the roof/ceiling
				// select primary light rays oriented away from the sun/moon; doesn't work well due to reduced ray scattering
				light_dir     = get_light_pos().get_norm(); // more accurate, but requires indir to be recomputed when sun/moon pos changes
				//light_dir     = plus_z; // make it vertical so that it doesn't need to be updated when the sun/moon pos changes
				lcolor        = cur_ambient*2.0; // split rays into two groups for ambient and diffuse
				pri_lcolor    = cur_diffuse;
				dir           = 1; // pointed up
			}
			else { // normal window
				assert(window.ix < 4); // encodes 2*dim + dir
				dim =  bool(window.ix >> 1);
				dir = !bool(window.ix &  1); // cast toward the interior
				surface_area = window.dz()*window.get_sz_dim(!bool(dim));
				light_cube.translate_dim(dim, (dir ? 1.0 : -1.0)*0.5*b.get_wall_thickness()); // shift slightly inside the building to avoid collision with the exterior wall
				lcolor = outdoor_color;
			}
			// light intensity scales with surface area, since incoming light is a constant per unit area (large windows = more light)
			weight *= surface_area/0.0016f; // a fraction the surface area weight of lights
		}
		else { // room light or lamp, pointing downward
			vect_room_object_t const &objs(b.interior->room_geom->objs);
			assert((unsigned)cur_light < objs.size());
			room_object_t const &ro(objs[cur_light]);
			bool const light_in_basement(ro.z1() < b.ground_floor_z1), is_lamp(ro.type == TYPE_LAMP);
			light_cube      = ro;
			light_cube.z1() = light_cube.z2() = (ro.z1() - 0.01*ro.dz()); // set slightly below bottom of light
			light_center    = light_cube.get_cube_center();
			in_attic        = ro.in_attic();
			in_ext_basement = (light_in_basement && b.point_in_extended_basement_not_basement(light_center));
			if (in_attic) {base_num_rays *= 4;} // more rays in attic, since light is large and there are only 1-2 of them
			if (is_lamp ) {base_num_rays /= 2;} // half the rays for lamps
			if (is_lamp ) {dir = 2;} // onmidirectional; dim stays at 2/Z
			float const surface_area(ro.dx()*ro.dy() + 2.0f*(ro.dx() + ro.dy())*ro.dz()); // bottom + 4 sides (top is occluded), 0.0003 for houses
			lcolor  = (is_lamp ? LAMP_COLOR : ro.get_color());
			weight *= surface_area/0.0003f;
			if (b.has_pri_hall())     {weight *= 0.70;} // floorplan is open and well lit, indir lighting value seems too high
			if (ro.type == TYPE_LAMP) {weight *= 0.33;} // lamps are less bright
			if (light_in_basement)    {weight *= ((b.has_parking_garage && !in_ext_basement) ? 0.25 : 0.5);} // basement is darker, parking garages are even darker
			if (in_attic)             {weight *= ATTIC_LIGHT_RADIUS_SCALE*ATTIC_LIGHT_RADIUS_SCALE;} // based on surface area rather than radius
			if (ro.is_round())        {light_radius = ro.get_radius();}
		}
		if (b.is_house)        {weight *=  2.0;} // houses have dimmer lights and seem to work better with more indir
		if (is_negative_light) {weight *= -1.0;}
		weight /= base_num_rays; // normalize to the number of rays
		unsigned NUM_PRI_SPLITS(is_window ? 4 : 16); // we're counting primary rays for windows, use fewer primary splits to reduce noise at the cost of increased time
		max_eq(base_num_rays, NUM_PRI_SPLITS);
		int const num_rays(base_num_rays/NUM_PRI_SPLITS);
		building_colors_t bcolors;
		b.set_building_colors(bcolors);
		
		// Note: dynamic scheduling is faster, and using blocks doesn't help
#pragma omp parallel for schedule(dynamic) num_threads(num_rt_threads)
		for (int n = 0; n < num_rays; ++n) {
			if (kill_thread) continue;
			rand_gen_t rgen;
			rgen.set_state(n+1, cur_light); // should be deterministic, though add_path_to_lmcs() is not (due to thread races)
			vector3d pri_dir;
			colorRGBA ray_lcolor(lcolor), ccolor(WHITE);
			bool const is_skylight_dir(is_skylight && (n&1)); // alternate between sky ambient and sun/moon directional
			
			if (is_skylight_dir) { // skylight directional diffuse
				pri_dir    = light_dir;
				ray_lcolor = pri_lcolor;
			}
			else { // omidirectional or sky ambient from windows
				pri_dir = rgen.signed_rand_vector_spherical().get_norm(); // should this be cosine weighted for windows?
				if (is_window && ((pri_dir[dim] > 0.0) ^ dir)) {pri_dir[dim] *= -1.0;} // reflect light if needed about window plane to ensure it enters the room
				//if (!is_window && dim == 2 && dir == 2 && pri_dir.z > 0.0) {pri_dir.z = -pri_dir.z;} // must point down
			}
			float const lum_thresh(0.1*ray_lcolor.get_luminance());
			point origin, init_cpos, cpos;
			vector3d init_cnorm, cnorm;

			// select a random point on the light cube
			for (unsigned N = 0; N < 10; ++N) { // 10 attempts to find a point within the light shape
				for (unsigned d = 0; d < 3; ++d) {
					float const lo(light_cube.d[d][0]), hi(light_cube.d[d][1]);
					origin[d] = ((lo == hi) ? lo : rgen.rand_uniform(lo, hi));
				}
				if (light_radius == 0.0 || dist_xy_less_than(origin, light_center, light_radius)) break; // done/success
			} // for N
			init_cpos = origin; // init value
			bool const hit(b.ray_cast_interior(origin, pri_dir, valid_area, bvh, in_attic, in_ext_basement, bcolors, init_cpos, init_cnorm, ccolor, &rgen));

			// room lights already contribute direct lighting, so we skip this ray; however, windows don't, so we add their primary ray contribution
			if (is_window && /*!is_skylight_dir*/!is_skylight && init_cpos != origin) {
				point const p1(origin*ray_scale + llc_shift), p2(init_cpos*ray_scale + llc_shift); // transform building space to global scene space
				add_path_to_lmcs(&lmgr, nullptr, p1, p2, weight, ray_lcolor*NUM_PRI_SPLITS, LIGHTING_LOCAL, 0); // local light, no bcube; scale color based on splits
			}
			if (!hit) continue; // done
			colorRGBA const init_color(ray_lcolor.modulate_with(ccolor));
			if (init_color.get_luminance() < lum_thresh) continue; // done (Note: get_weighted_luminance() will discard too much blue light)
			vector3d const v_ref(get_reflect_dir(pri_dir, init_cnorm));

			for (unsigned splits = 0; splits < NUM_PRI_SPLITS; ++splits) {
				point pos(origin);
				vector3d dir(pri_dir);
				colorRGBA cur_color(init_color);
				calc_reflect_ray(pos, init_cpos, dir, init_cnorm, v_ref, rgen, tolerance);

				for (unsigned bounce = 1; bounce < MAX_RAY_BOUNCES; ++bounce) { // allow up to MAX_RAY_BOUNCES bounces
					cpos = pos; // init value
					bool const hit(b.ray_cast_interior(pos, dir, valid_area, bvh, in_attic, in_ext_basement, bcolors, cpos, cnorm, ccolor, &rgen));

					if (cpos != pos) { // accumulate light along the ray from pos to cpos (which is always valid) with color cur_color
						point const p1(pos*ray_scale + llc_shift), p2(cpos*ray_scale + llc_shift); // transform building space to global scene space
						add_path_to_lmcs(&lmgr, nullptr, p1, p2, weight, cur_color, LIGHTING_LOCAL, 0); // local light, no bcube
					}
					if (!hit) break; // done
					cur_color = cur_color.modulate_with(ccolor);
					if (cur_color.get_luminance() < lum_thresh) break; // done
					calc_reflect_ray(pos, cpos, dir, cnorm, get_reflect_dir(dir, cnorm), rgen, tolerance);
				} // for bounce
			} // for splits
		} // for n
		is_running = 0; // flag as done
	}
	void wait_for_finish(bool force_kill) {
		// Note: for now the time taken to process a light should be pretty fast so we just block until finished; set kill_thread=1 to be faster
		if (force_kill) {kill_thread = 1;}
		while (is_running) {sleep_for_ms(10);}
		kill_thread = 0;
	}
	void update_volume_light_texture() { // full update, 6.6ms for z=128
		init_lmgr(0); // init on first call; clear_lighting=0
		//highres_timer_t timer("Lighting Tex Create");
		indir_light_tex_from_lmap(cur_tid, lmgr, tex_data, MESH_X_SIZE, MESH_Y_SIZE, MESH_SIZE[2], indir_light_exp, 1); // local_only=1
	}
	void maybe_join_thread() {
		if (needs_to_join) {rt_thread.join(); needs_to_join = 0;}
	}
	void add_to_remove_queue(unsigned light_ix) {
		for (unsigned v : remove_queue) {
			if (v == light_ix) return; // already in the queue
		}
		remove_queue.push_back(light_ix);
	}
	void add_window_lights(building_t const &b, point const &target) {
		float const window_vspacing(b.get_window_vspace());
		int const target_room(b.get_room_containing_pt(target)); // generally always should be >= 0

		for (auto i = windows.begin(); i != windows.end(); ++i) {
			if (i->ix & 4) continue; // blocked bit is set, skip
			point const center(i->get_cube_center());
			if (center.z < valid_area.z1() || center.z > valid_area.z2()) continue; // wrong floor
			float dist_sq(p2p_dist_sq(center, target));
			float const surface_area(i->dx()*i->dy() + i->dx()*i->dz() + i->dy()*i->dz());
			dist_sq *= 0.05f*window_vspacing/surface_area; // account for the size of the window, larger window smaller/higher priority
			if (target_room >= 0 && b.get_room_containing_pt(center) == target_room) {dist_sq *= 0.1;} // prioritize the room the player is in
			lights_to_sort.emplace_back(dist_sq, ((i - windows.begin()) | IS_WINDOW_BIT));
		} // for i
	}
	void get_windows(building_t const &b) {
		b.get_all_windows(windows);
		update_windows = 0;
	}
	void sort_lights_by_priority() {
		light_ids.clear();
		sort(lights_to_sort.begin(), lights_to_sort.end()); // sort by increasing distance
		for (auto const &light : lights_to_sort) {light_ids.push_back(light.second);}
		lights_to_sort.clear();
	}
public:
	building_indir_light_mgr_t() : is_running(0), kill_thread(0), lighting_updated(0), needs_to_join(0), need_bvh_rebuild(0),
		update_windows(0), is_negative_light(0), in_ext_basement(0), cur_bix(-1), cur_light(-1), cur_floor(-1), cur_tid(0) {}

	cube_t get_light_bounds() const {return light_bounds;}

	void invalidate_lighting() {
		is_negative_light = in_ext_basement = 0;
		cur_light = -1;
		remove_queue.clear();
		lights_complete.clear();
		lights_seen.clear();
		end_rt_job();
		lmgr.reset_all(); // clear lighting values back to 0
	}
	void clear() {
		lighting_updated = need_bvh_rebuild = update_windows = 0;
		cur_bix = cur_floor = -1;
		invalidate_lighting();
		tex_data.clear();
		light_ids.clear();
		lights_to_sort.clear();
		windows.clear();
		update_volume_light_texture(); // reset lighting from prev building, or reset to dark when entering first building
		bvh.clear();
	}
	void end_rt_job() {
		wait_for_finish(1); // force_kill=1
		maybe_join_thread();
	}
	void free_indir_texture() {free_texture(cur_tid);}

	void register_cur_building(building_t const &b, unsigned bix, point const &target, unsigned &tid) { // target is in building space
		bool floor_change(0), need_rebuild(0);

		if ((int)bix != cur_bix) { // change to a different building
			clear();
			cur_bix = bix;
			assert(!is_running);
			build_bvh(b, target);
			get_windows(b);
		}
		else {
			if (update_windows) {get_windows(b);}
			unsigned const new_floor(b.get_floor_for_zval(target.z));
			bool const new_in_ext_basement(b.point_in_extended_basement_not_basement(target));
			floor_change  = (cur_floor >= 0 && cur_floor != (int)new_floor);
			floor_change |= (new_in_ext_basement != in_ext_basement); // treat extended basement threshold cross as a floor change
			need_rebuild  = (floor_change && !light_bounds.contains_cube(get_valid_area(b, target, new_floor)));
		}
		if (need_rebuild) { // handle floor change
			invalidate_lighting();
			build_bvh(b, target);
		}
		calc_cur_ambient_diffuse(); // needed for correct outdoor color
		colorRGBA const cur_outdoor_color(get_outdoor_light_color());

		if (!windows.empty() && cur_outdoor_color != outdoor_color) {
			// outdoor color change, need to update lighting
			// Note1: we could remove and re-add window lights, but that may be more work than clearing and re-adding both lights and windows
			// Note2: we could store separate volumes for indoor vs. outdoor lights, but that would add some overhead per-update, while sun/moon changes would be uncommon
			invalidate_lighting();
			outdoor_color = cur_outdoor_color;
		}
		if (display_framerate && (is_running || lighting_updated)) { // show progress to the user
			std::ostringstream oss;
			oss << "Lights: " << lights_complete.size() << " / " << (max(light_ids.size(), lights_seen.size()) + remove_queue.size());
			lighting_update_text = oss.str();
		}
		if (is_running) return; // still running, let it continue

		if (lighting_updated) { // update lighting texture based on incremental progress
			maybe_join_thread();
			update_volume_light_texture();
			lighting_updated = 0;
		}
		// nothing is running and there is more work to do, find the nearest light to the target and process it
		if (!need_rebuild) {need_bvh_rebuild |= floor_change;} // rebuild on player floor change if not rebuilt above
		if (need_bvh_rebuild) {build_bvh(b, target);}
		
		if (cur_light >= 0) {
			if (!is_negative_light) {lights_complete.insert(cur_light);} // mark the most recent light as complete if not a light removal
			cur_light = -1;
		}
		if (!remove_queue.empty()) { // remove an existing light; must run even when player_in_elevator>=2 (doors closed/moving) to remove elevator light at old pos
			cur_light = remove_queue.front();
			remove_queue.pop_front();
			is_negative_light = 1;
		}
		else { // find a new light to add
			if (player_in_elevator >= 2) return; // pause updates for player in closed elevator since lighting is not visible
			is_negative_light = 0; // back to normal positive lights
			b.get_lights_with_priorities(target, valid_area, lights_to_sort);
			add_window_lights(b, target);
			sort_lights_by_priority();

			for (auto i = light_ids.begin(); i != light_ids.end(); ++i) {
				lights_seen.insert(*i); // must track lights across all floors seen for correct progress update
				if (cur_light < 0 && lights_complete.find(*i) == lights_complete.end()) {cur_light = *i;} // find an incomplete light
			}
		}
		if (cur_light >= 0) {start_lighting_compute(b);} // this light is next
		tid = cur_tid;
	}
	void register_light_state_change(unsigned light_ix, bool light_is_on, bool in_elevator, bool geom_changed) {
		if (geom_changed) {
			// TODO: have to first remove the lighting with the old geom, and then re-add it - but the geometry has already been updated
			return;
		}
		unsigned const num_erased(lights_complete.erase(light_ix)); // light is no longer completed; erase its state
		bool const is_cur_light((int)light_ix == cur_light && is_running);
		// Note: we can't just stop in the middle, because that will leave cur_light in an invalid/incomplete state
		// Note: if door state changed since this light was turned on, removing it may leave some light
		if ((geom_changed || !light_is_on || (in_elevator && num_erased)) && (num_erased || is_cur_light)) {add_to_remove_queue(light_ix);} // must remove the light instead
	}
	static cube_t get_valid_area(building_t const &b, point const &target, unsigned target_floor) {
		cube_t VA;

		// clip per light source to current floor; note that this will exclude stairs going up or down
		if (b.point_in_attic(target)) {
			VA = b.get_attic_part();
			set_cube_zvals(VA, b.interior->attic_access.z1(), b.interior_z2);
		}
		else {
			bool const in_ext_basement(b.point_in_extended_basement_not_basement(target));
			float const floor_spacing(b.get_window_vspace()), building_z1(b.get_bcube_z1_inc_ext_basement());
			VA = (in_ext_basement ? b.interior->basement_ext_bcube : b.bcube);
			VA.z1() = building_z1 + target_floor*floor_spacing;
			VA.z2() = VA.z1() + floor_spacing;

			// handle multi-floor rooms; this is only called once when enabling indir lighting,
			// so it will handle starting in a tall room, but it won't work if the player walks into a tall room later;
			// but, since the front door generally opens into the living room, which is often the tall room, it will work in many cases
			if (in_ext_basement) {
				if (b.has_pool()) { // check if light source room has a pool, and include the pool in our valid area
					room_t const &pool_room(b.get_room(b.interior->pool.room_ix));

					if (pool_room.contains_pt(target) || b.interior->pool.contains_pt(target)) {
						min_eq(VA.z1(), b.interior->pool.orig_z1); // bottom of the pool
						max_eq(VA.z2(), pool_room.z2()); // ceiling above the pool
					}
				}
			}
			else if (b.has_tall_retail()) { // handle lights on tall retail ceilings
				cube_t const &retail_room(b.get_retail_part());
				if (retail_room.contains_pt(target)) {max_eq(VA.z2(), retail_room.z2());} // use ceiling of retail part
			}
			else if (target.z > b.ground_floor_z1) { // not in the basement - check for tall rooms
				int const target_room(b.get_room_containing_pt(target));
				
				if (target_room >= 0) {
					room_t const &room(b.get_room(target_room));
					if (room.is_single_floor) {max_eq(VA.z2(), room.z2());} // include room height
				}
			}
		}
		return VA;
	}
	void build_bvh(building_t const &b, point const &target) {
		//highres_timer_t timer("Build BVH");
		cur_floor  = b.get_floor_for_zval(target.z);
		valid_area = get_valid_area(b, target, cur_floor);
		in_ext_basement = b.point_in_extended_basement_not_basement(target);

		if (in_ext_basement) { // extended basement
			light_bounds = valid_area;
			set_cube_zvals(light_bounds, b.interior->basement_ext_bcube.z1(), b.interior->basement_ext_bcube.z2()); // cover the entire extended basement range
			// very high aspect ratio cubes cause banding artifacts in lighting, so increase the height if needed to avoid this;
			// note that there's significant loss of vertical resolution, though lighting should be a bit faster
			max_eq(light_bounds.z2(), (light_bounds.z1() + 0.25f*max(light_bounds.dx(), light_bounds.dy())));
		}
		else if (INDIR_LIGHT_FLOOR_SPAN == 0) {light_bounds = b.get_interior_bcube(0);} // unlimited floor span/entire building; inc_ext_basement=0
		else if (INDIR_LIGHT_FLOOR_SPAN == 1) {light_bounds = valid_area;} // single floor only
		else if (light_bounds.contains_cube(valid_area)) {} // valid area already contained, no update needed (will fail if light_bounds is not set)
		else { // limited number of floors
			float const dist_above_below(0.5f*(INDIR_LIGHT_FLOOR_SPAN - 1)*b.get_window_vspace());
			cube_t const interior_bcube(b.get_interior_bcube(0)); // inc_ext_basement=0
			light_bounds = valid_area;
			light_bounds.expand_in_dim(2, dist_above_below); // expand in Z
			// clamp Z range to interior_bcube and add any excess height to the other end of the range
			if      (light_bounds.z1() < interior_bcube.z1()) {light_bounds.z2() += (interior_bcube.z1() - light_bounds.z1());}
			else if (light_bounds.z2() > interior_bcube.z2()) {light_bounds.z1() -= (light_bounds.z2() - interior_bcube.z2());}
			max_eq(light_bounds.z1(), interior_bcube.z1());
			min_eq(light_bounds.z2(), interior_bcube.z2());
			// Note: if we get here, we assume lighting has already been invalidated
		}
		bvh.clear();
		b.gather_interior_cubes(bvh.get_objs(), valid_area);
		if (b.point_in_attic(target)) {b.get_attic_roof_tquads(bvh.roof_tquads);}
		bvh.build_tree_top(0); // verbose=0
		need_bvh_rebuild = 0;
	}
	void invalidate_bvh    () {need_bvh_rebuild = 1;} // Note: can't directly clear bvh because a thread may be using it
	void invalidate_windows() {update_windows = 1;}
	cube_bvh_t const &get_bvh() const {return bvh;}
};

building_indir_light_mgr_t building_indir_light_mgr;

void free_building_indir_texture() {building_indir_light_mgr.free_indir_texture();}
void end_building_rt_job() {building_indir_light_mgr.end_rt_job();}
cube_t get_building_indir_light_bounds() {return building_indir_light_mgr.get_light_bounds();}

void building_t::create_building_volume_light_texture(unsigned bix, point const &target, unsigned &tid) const {
	if (!has_room_geom()) return; // error?
	building_indir_light_mgr.register_cur_building(*this, bix, target, tid);
}

bool building_t::ray_cast_camera_dir(point const &camera_bs, point &cpos, colorRGBA &ccolor) const { // unused - for debugging; excludes attic and extended basement
	assert(!USE_BKG_THREAD); // not legal to call when running lighting in a background thread
	building_indir_light_mgr.build_bvh(*this, camera_bs);
	building_colors_t bcolors;
	set_building_colors(bcolors);
	vector3d cnorm; // unused
	return ray_cast_interior(camera_bs, cview_dir, bcube, building_indir_light_mgr.get_bvh(), point_in_attic(camera_bs),
		point_in_extended_basement_not_basement(camera_bs), bcolors, cpos, cnorm, ccolor);
}

// Note: target is building space camera
void building_t::get_lights_with_priorities(point const &target, cube_t const &valid_area, vector<pair<float, unsigned>> &lights_to_sort) const {
	if (!has_room_geom()) return; // error?
	//if (is_rotated()) {} // do we need to handle this case?
	vect_room_object_t const &objs(interior->room_geom->objs);
	float const window_vspacing(get_window_vspace()), diag_dist_sq(bcube.dx()*bcube.dx() + bcube.dy()*bcube.dy()), other_floor_penalty(0.25*diag_dist_sq);
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	int const target_room(get_room_containing_pt(target)); // generally always should be >= 0

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_type() || !i->is_light_on() || i->is_broken2()) continue; // not a light, or light not on, or fully broken
		bool const light_in_basement(i->z1() < ground_floor_z1);
		if (!check_indir_enabled(light_in_basement, i->in_attic())) continue;
		if (!valid_area.contains_cube(*i)) continue; // outside valid area

		if (i->in_elevator()) { // elevator light
			if (get_elevator(i->obj_id).may_be_moving()) continue; // possibly moving elevator or elevator doors, don't update light yet
		}
		point const center(i->get_cube_center());
		float dist_sq(p2p_dist_sq(center, target));
		dist_sq *= 0.005f*window_vspacing/(i->dx()*i->dy()); // account for the size of the light, larger lights smaller/higher priority

		if (i->z1() < target.z || i->z2() > (target.z + window_vspacing)) { // penalty if on a different floor
			dist_sq += (i->has_stairs() ? 0.25 : 1.0)*other_floor_penalty; // less penalty for lights on stairs
		}
		if (target_room >= 0 && get_room_containing_pt(center) == target_room) {dist_sq *= 0.1;} // prioritize the room the player is in
		// reduce distance for lights visible to target?
		lights_to_sort.emplace_back(dist_sq, (i - objs.begin()));
	} // for i
}

bool get_wall_quad_window_area(vect_vnctcc_t const &wall_quad_verts, unsigned i, cube_t &c, float &tx1, float &tx2, float &tz1, float &tz2) {
	auto const &v0(wall_quad_verts[i]);
	c = cube_t(v0.v);
	tx1 = tx2 = v0.t[0]; tz1 = tz2 = v0.t[1]; // tex coord ranges (xy, z); should generally be whole integers

	for (unsigned j = 1; j < 4; ++j) {
		auto const &vj(wall_quad_verts[i + j]);
		c.union_with_pt(vj.v);
		min_eq(tx1, vj.t[0]);
		max_eq(tx2, vj.t[0]);
		min_eq(tz1, vj.t[1]);
		max_eq(tz2, vj.t[1]);
	}
	if (tx1 == tx2 || tz1 == tz2) return 0; // wall is too small to contain a window
	assert(tx2 - tx1 < 1000.0f && tz2 - tz1 < 1000.0f); // sanity check - less than 1000 windows in each dim
	assert(c.dz() > 0.0);
	return 1;
}

// get windows used as indir lighting sources
void building_t::get_all_windows(vect_cube_with_ix_t &windows) const { // Note: ix encodes 2*dim+dir
	windows.clear();
	if (!has_int_windows() || is_rotated() || !is_cube()) return; // no windows; rotated and non-cube buildings are not handled
	float const window_h_border(WINDOW_BORDER_MULT*get_window_h_border()), window_v_border(WINDOW_BORDER_MULT*get_window_v_border()); // (0, 1) range
	float const wall_thickness(get_wall_thickness());
	vect_room_object_t blinds;

	if (is_house && has_room_geom()) { // find all bedroom blinds and use them to partially block windows
		for (room_object_t const &c : interior->room_geom->objs) {
			if (c.type != TYPE_BLINDS) continue;
			blinds.push_back(c);
			blinds.back().expand_in_dim(c.dim, wall_thickness); // make sure the intersect the windows
		}
	}
	static vect_vnctcc_t wall_quad_verts;
	wall_quad_verts.clear();
	get_all_drawn_window_verts_as_quads(wall_quad_verts);

	for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
		cube_t c;
		float tx1, tx2, tz1, tz2;
		if (!get_wall_quad_window_area(wall_quad_verts, i, c, tx1, tx2, tz1, tz2)) continue;
		bool const dim(c.dy() < c.dx()), dir(wall_quad_verts[i].get_norm()[dim] > 0.0);
		assert(c.get_sz_dim(dim) == 0.0); // must be zero size in one dim (X or Y oriented); could also use the vertex normal
		float const d_tx_inv(1.0f/(tx2 - tx1)), d_tz_inv(1.0f/(tz2 - tz1));
		float const window_width(c.get_sz_dim(!dim)*d_tx_inv), window_height(c.dz()*d_tz_inv); // window_height should be equal to window_vspacing
		float const border_xy(window_width*window_h_border), border_z(window_height*window_v_border);
		cube_t window(c); // copy dim <dim>

		for (float z = tz1; z < tz2; z += 1.0) { // each floor
			float const bot_edge(c.z1() + (z - tz1)*window_height);
			set_cube_zvals(window, bot_edge+border_z, bot_edge+window_height-border_z);

			for (float xy = tx1; xy < tx2; xy += 1.0) { // windows along each wall
				float const low_edge(c.d[!dim][0] + (xy - tx1)*window_width);
				window.d[!dim][0] = low_edge + border_xy;
				window.d[!dim][1] = low_edge + window_width - border_xy;

				for (room_object_t const &b : blinds) {
					if (!b.intersects(window)) continue;
					
					if (b.is_hanging()) { // horizontal (comes from the top)
						min_eq(window.z2(), b.z1());
					}
					else { // vertical (comes from the sides)
						bool const side(b.get_center_dim(!dim) < window.get_center_dim(!dim));
						if (side) {max_eq(window.d[!dim][0], b.d[!dim][1]);} // left  side
						else      {min_eq(window.d[!dim][1], b.d[!dim][0]);} // right side
					}
				} // for b
				bool is_blocked(window.dz() <= 0.0 || window.get_sz_dim(!dim) <= 0.0);
				cube_t window_clipped(window);

				if (!is_blocked && !walkways.empty()) {
					// check for windows blocked by walkways and clip them, even though some light could come through if not blocked by a door
					float const wz(window.zc()), wpos(window_clipped.d[dim][dir]);
					float &wlo(window_clipped.d[!dim][0]), &whi(window_clipped.d[!dim][1]);
					assert(wlo < whi);

					for (building_walkway_t const &w : walkways) {
						if (wz < w.bcube.z1() || wz > w.bcube.z2()) continue; // no Z overlap
						if (wpos < w.bcube.d[dim][0] - wall_thickness || wpos > w.bcube.d[dim][1] + wall_thickness) continue; // wrong wall
						float const blo(w.bcube.d[!dim][0]), bhi(w.bcube.d[!dim][1]);
						if (whi <= blo || wlo >= bhi) continue; // no overlap
						if (wlo >= blo && whi <= bhi) {is_blocked = 1; break;} // fully contained
						if      (wlo > blo && wlo < bhi) {wlo = bhi;} // clip low  edge
						else if (whi > blo && whi < bhi) {whi = blo;} // clip high edge
						assert(wlo < whi);
					} // for w
				}
				windows.emplace_back(window_clipped, (4*is_blocked + 2*dim + dir)); // Note: must include blocked window for seen cache to work
			} // for xy
		} // for z
	} // for i
	if (get_light_pos().z > 0.0) { // if primary light (sun/moon) is above the horizon, add skylight indir lighting
		for (cube_t const &skylight : skylights) {windows.emplace_back(skylight, 0);} // add skylights as vertical windows with ix=0
	}
}

bool building_t::register_indir_lighting_state_change(unsigned light_ix, bool is_door_change) const {
	if (!enable_building_indir_lighting()) return 0; // no update needed
	assert(has_room_geom());
	assert(light_ix < interior->room_geom->objs.size());
	room_object_t const &obj(interior->room_geom->objs[light_ix]);
	assert(obj.is_light_type());
	if (is_door_change && !obj.is_light_on()) return 0; // light off, no state change
	building_indir_light_mgr.register_light_state_change(light_ix, obj.is_light_on(), obj.in_elevator(), is_door_change);
	return 1;
}
void building_t::register_indir_lighting_geom_change() const {
	building_indir_light_mgr.invalidate_bvh();
}
void building_t::register_blinds_state_change() const {
	building_indir_light_mgr.invalidate_windows();
	register_indir_lighting_geom_change();
}

bool building_t::is_light_occluded(point const &lpos, point const &camera_bs) const {
	assert(interior);
	// Note: assumes the light is inside the building
	// exterior walls have windows and don't generally occlude lights; room objects and doors are too small to occlude; elevators are too sparse to occlude
	if (line_intersect_walls(lpos, camera_bs)) return 1; // check interior walls
	if (line_int_cubes(lpos, camera_bs, interior->fc_occluders, cube_t(lpos, camera_bs))) return 1; // check floors and ceilings
	return 0;
}
void building_t::clip_ray_to_walls(point const &p1, point &p2, vect_cube_t const walls[2]) const { // Note: assumes p1.z == p2.z
	for (unsigned d = 0; d < 2; ++d) {
		bool const dir(p2[d] < p1[d]);
		float mult((p2[!d] - p1[!d])/(p2[d] - p1[d]));

		for (auto c = walls[d].begin(); c != walls[d].end(); ++c) {
			float const v(c->d[d][dir]);
			if ((p1[d] < v) == (p2[d] < v)) continue; // no crossing
			if (p1.z < c->z1() || p1.z > c->z2()) continue; // no z overlap
			float const lo(c->d[!d][0]), hi(c->d[!d][1]); // expand by wall_edge_spacing
			if ((p1[!d] < lo && p2[!d] < lo) || (p1[!d] > hi && p2[!d] > hi)) continue; // both to one side of the cube
			float const cross_pos(p1[!d] + mult*(v - p1[d]));
			if (cross_pos < lo || cross_pos > hi) continue; // doesn't cross within the cube bounds
			p2[ d] = v;
			p2[!d] = cross_pos;
			mult   = (p2[!d] - p1[!d])/(p2[d] - p1[d]);
		} // for c
	} // for d
}

bool do_line_clip_xy_p2(point const &p1, point &p2, cube_t const &c) {
	float tmin(0.0), tmax(1.0);
	if (!get_line_clip_xy(p1, p2, c.d, tmin, tmax)) return 0;
	p2 = p1 + tmax*(p2 - p1);
	return 1;
}

void building_t::refine_light_bcube(point const &lpos, float light_radius, room_t const &room, cube_t &light_bcube, bool is_parking_garage) const {
	// base: 173613 / bcube: 163942 / clipped bcube: 161455 / tight: 159005 / rays: 101205 / no ls bcube expand: 74538
	// starts with building bcube clipped to light bcube
	//highres_timer_t timer("refine_light_bcube"); // 0.035ms average
	assert(has_room_geom());
	cube_t tight_bcube, part;
	static vect_cube_t other_parts, walls[2];
	other_parts.clear();

	// first determine the union of all intersections with parts; ignore zvals here so that we get the same result for every floor
	if (lpos.z < ground_floor_z1) { // light in basement, including parking garage
		assert(has_basement());
		cube_t const basement(get_basement());

		if (!is_parking_garage && point_in_extended_basement_not_basement(lpos)) { // light in extended basement
			tight_bcube = part = interior->basement_ext_bcube;
			if (light_bcube.intersects_xy(basement)) {tight_bcube.union_with_cube(basement);}
			other_parts.push_back(basement); // include the basement itself
		}
		else {
			tight_bcube = part = basement;

			if (has_ext_basement()) { // include extended basement
				if (light_bcube.intersects_xy(interior->basement_ext_bcube)) {tight_bcube.union_with_cube(interior->basement_ext_bcube);}
				other_parts.push_back(interior->basement_ext_bcube);
			}
		}
	}
	else { // light above ground
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // secondary buildings aren't handled here
			//if (lpos.z < p->z1() || lpos.z > p->z2()) continue; // light zval not included (doesn't seem to work)
			if (!light_bcube.intersects_xy(*p)) continue;
			cube_t c(light_bcube);
			c.intersect_with_cube_xy(*p);
			if (tight_bcube.is_all_zeros()) {tight_bcube = c;} else {tight_bcube.union_with_cube_xy(c);}
			if (p->contains_pt(lpos)) {part = *p;} else {other_parts.push_back(*p);}
		} // for p
	}
	assert(!part.is_all_zeros());
	tight_bcube.z1() = light_bcube.z1();
	tight_bcube.z2() = light_bcube.z2();
	// next cast a number of horizontal rays in a circle around the light to see how far they reach; any walls hit occlude the light and reduce the bcube;
	// it's important that we don't use any tests against parts/rooms that depend on the light zval because the result must be valid for all lights in the stack
	unsigned const NUM_RAYS = 180; // every 2 degrees
	if (NUM_RAYS == 0) {light_bcube = tight_bcube; return;}
	cube_t rays_bcube(lpos, lpos), room_exp(get_walkable_room_bounds(room));
	float const wall_thickness(get_wall_thickness()), tolerance(0.01*wall_thickness);
	room_exp.expand_by_xy(wall_thickness + tolerance); // to include points on the border + some FP error
	// pre-compute the nearby walls we will use for clipping
	for (unsigned d = 0; d < 2; ++d) {walls[d].clear();}
	tight_bcube.z1() = tight_bcube.z2() - get_floor_ceil_gap(); // limit to a single floor to exclude walls on the floor below (for backrooms)

	if (is_pos_in_pg_or_backrooms(lpos)) {
		index_pair_t start, end;
		get_pgbr_wall_ix_for_pos(lpos, start, end);

		for (unsigned d = 0; d < 2; ++d) {
			vect_cube_t const &pbgr_walls(interior->room_geom->pgbr_walls[d]);

			for (auto w = pbgr_walls.begin()+start.ix[d]; w != pbgr_walls.begin()+end.ix[d]; ++w) {
				if (tight_bcube.intersects(*w)) {walls[d].push_back(*w);}
			}
		}
	}
	if (!is_parking_garage || interior->has_backrooms) { // still need to check for backrooms to handle wall adjacent to parking garage
		for (unsigned d = 0; d < 2; ++d) {
			for (cube_t const &c : interior->walls[d]) {
				if (tight_bcube.intersects(c)) {walls[d].push_back(c);}
			}
		}
	}
	point ray_origin(lpos);
	// if this is a tall room, lower the ray origin to the ground floor to include light paths through ground floor doorways
	bool use_first_floor(room.is_single_floor);
	float const floor_spacing(get_window_vspace());

	if (!use_first_floor && lpos.z - room.z1() > 1.5*floor_spacing) { // light above the ground floor
		// check if adjacent to a tall room, in which case we also need to lower the ray origin
		for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
			if (!r->is_single_floor || !are_rooms_connected(room, *r, (room.z1() + 0.5*floor_spacing), 0)) continue; // check_door_open=0
			use_first_floor = 1;
			break;
		}
	}
	if (use_first_floor) {min_eq(ray_origin.z, (room.z1() + get_floor_ceil_gap()));}

	for (unsigned n = 0; n < NUM_RAYS; ++n) {
		float const angle(TWO_PI*n/NUM_RAYS);
		point p2(ray_origin + point(light_radius*sin(angle), light_radius*cos(angle), 0.0));
		// test for bad rays; this can fail on rotated buildings if lights aren't rotated properly, and in other cases when I'm experimenting, so it's allowed
		if (!do_line_clip_xy_p2(ray_origin, p2, tight_bcube)) continue; // bad ray, skip

		if (other_parts.empty() || part.contains_pt_xy(p2)) {
			clip_ray_to_walls(ray_origin, p2, walls); // the simple case where we don't need to handle part boundaries (optimization)
		}
		else {
			// find the point where this ray exits the building by following it through all parts; parts should be exactly adjacent to each other in X or Y
			point cur_pt(p2);
			bool const ret(do_line_clip_xy_p2(ray_origin, cur_pt, part)); // exit point of the starting part
			assert(ret);
			auto prev_part(other_parts.end());

			for (unsigned n = 0; n < other_parts.size(); ++n) { // could be while(1), but this is safer in case we run into an infinite loop due to FP errors
				bool found(0);

				for (auto p = other_parts.begin(); p != other_parts.end(); ++p) {
					if (p == prev_part || !p->contains_pt_xy_exp(cur_pt, tolerance)) continue; // ray does not continue into this new part
					point new_pt(p2);
					if (do_line_clip_xy_p2(cur_pt, new_pt, *p)) {cur_pt = new_pt; prev_part = p; found = 1; break;} // ray continues into this part
				}
				if (!found || cur_pt == p2) break; // ray has exited the building, or we've reached it's end point; done
			} // for n
			p2 = cur_pt;
			if (room_exp.contains_pt_xy_inclusive(p2)) {} // ray ends in this room; it may have exited the building from this room
			else {clip_ray_to_walls(ray_origin, p2, walls);} // ray ends in another room, need to clip it to the building walls
		}
		rays_bcube.union_with_pt(p2);
	} // for n
	light_bcube = rays_bcube;
}

void assign_light_for_building_interior(light_source &ls, void const *obj, cube_t const &light_bcube, bool cache_shadows, bool is_lamp=0) {
	unsigned const smap_id(uintptr_t(obj)/sizeof(void *) + is_lamp); // use memory address as a unique ID; add 1 for lmaps to keep them separate
	ls.assign_smap_mgr_id(1); // use a different smap manager than the city (cars + streetlights) so that they don't interfere with each other
	if (!light_bcube.is_all_zeros()) {ls.set_custom_bcube(light_bcube);}
	ls.assign_smap_id(cache_shadows ? smap_id : 0); // if cache_shadows, mark so that shadow map can be reused in later frames
	if (!cache_shadows) {ls.invalidate_cached_smap_id(smap_id);}
}
bool setup_light_for_building_interior(light_source &ls, room_object_t &obj, cube_t const &light_bcube, bool force_smap_update, unsigned shadow_caster_hash) {
	// If there are no dynamic shadows, we can reuse the previous frame's shadow map;
	// hashing object positions should handle the case where a shadow caster moves out of the light's influence and leaves a shadow behind;
	// also need to handle the case where the light is added on the frame the room geom is generated when the shadow map is not yet created;
	// requiring two consecutive frames of no dynamic objects should fix this
	uint16_t const sc_hash16(shadow_caster_hash ^ (shadow_caster_hash >> 16)); // combine upper and lower 16 bits into a single 16-bit value
	// cache if no objects moved (based on position hashing) this frame or last frame, and we're not forced to do an update
	bool const shadow_update(obj.item_flags != sc_hash16), cache_shadows(!shadow_update && !force_smap_update && (obj.flags & RO_FLAG_NODYNAM));
	assign_light_for_building_interior(ls, &obj, light_bcube, cache_shadows);
	if (shadow_update) {obj.flags &= ~RO_FLAG_NODYNAM;} else {obj.flags |= RO_FLAG_NODYNAM;} // store prev update state in object flag
	obj.item_flags = sc_hash16; // store current object hash in item flags
	return shadow_update;
}

cube_t building_t::get_rotated_bcube(cube_t const &c, bool inv_rotate) const {
	if (!is_rotated()) return c;
	point const center(bcube.get_cube_center());
	float const z(c.z2()); // top edge
	point corners[4] = {point(c.x1(), c.y1(), z), point(c.x2(), c.y1(), z), point(c.x2(), c.y2(), z), point(c.x1(), c.y2(), z)};
	
	for (unsigned n = 0; n < 4; ++n) {
		if (inv_rotate) {do_xy_rotate_inv(center, corners[n]);} else {do_xy_rotate(center, corners[n]);}
	}
	cube_t ret;
	ret.set_from_points(corners, 4);
	ret.z1() = c.z1();
	return ret;
}
bool building_t::is_rot_cube_visible(cube_t const &c, vector3d const &xlate, bool inc_mirror_reflections) const {
	cube_t const c_xf((is_rotated() ? get_rotated_bcube(c) : c) + xlate);
	if (camera_pdu.cube_visible(c_xf)) return 1;
	if (inc_mirror_reflections && cube_visible_in_building_mirror_reflection(c)) return 1;
	return 0;
}

float get_radius_for_room_light(room_object_t const &obj) {
	float radius(6.0f*(obj.dx() + obj.dy()));
	//if (obj.type == TYPE_LAMP) {radius *= 1.0;}
	if (obj.flags & RO_FLAG_ADJ_HI) {radius *= 2.0;} // wall lights have a larger radius since they're not as centered in the room
	if (obj.in_attic()) {radius *= ATTIC_LIGHT_RADIUS_SCALE;}
	return radius;
}

bool check_for_shadow_caster(vect_cube_with_ix_t const &cubes, cube_t const &light_bcube, point const &lpos,
	float dmax, bool has_stairs, vector3d const &xlate, unsigned &shadow_caster_hash)
{
	bool ret(0);

	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (lpos.z < c->z1()) continue; // light is below the object's bottom; assumes lights are spotlights pointed downward
		if (!c->intersects(light_bcube)) continue; // object not within light area of effect
		point const center(c->get_cube_center());
		if (dmax > 0.0 && !dist_less_than(lpos, center, dmax)) continue; // too far from light to cast a visible shadow
		
		if (!has_stairs) {
			// check for camera visibility of the union of the cube and the intersection points of the light rays
			// projected through the top 4 corners of the cube to the floor (cube z1 value)
			float const z(c->z2()); // top edge of the cube
			float const floor_z(light_bcube.z1()); // assume light z1 is on the floor; maybe could also use c->z1(), which at least works for people
			point const corners[4] = {point(c->x1(), c->y1(), z), point(c->x2(), c->y1(), z), point(c->x2(), c->y2(), z), point(c->x1(), c->y2(), z)};
			cube_t cube_ext(*c);

			for (unsigned n = 0; n < 4; ++n) {
				vector3d const delta(corners[n] - lpos);
				float const t((floor_z - lpos.z)/delta.z);
				cube_ext.union_with_pt(lpos + delta*t); // union with floor hit pos
			}
			if (!camera_pdu.cube_visible(cube_ext + xlate)) { // VFC
				hash_mix_point(c->get_size(), shadow_caster_hash); // hash size rather than position to include the object's presence but not its position (sort of)
				continue;
			}
		}
		shadow_caster_hash += c->ix;
		hash_mix_point(center, shadow_caster_hash);
		ret = 1;
	} // for c
	return ret;
}
void check_for_dynamic_shadow_casters(vector<person_t> const &people, vect_cube_with_ix_t &ped_bcubes, vect_cube_with_ix_t const &moving_objs,
	cube_t const &light_bcube, point const &lpos, float dmax, bool has_stairs, vector3d const &xlate, bool check_people, unsigned &shadow_caster_hash)
{
	if (check_people && animate2) { // update shadow_caster_hash for moving people, but not for lamps, because their light points toward the floor
		if (ped_bcubes.empty()) { // get all cubes on first light
			for (person_t const &p : people) {
				// if this person is waiting and their location isn't changing,
				// assume they have an idle animation playing and use the frame counter to make sure their shadows are updated each frame
				unsigned const ix((some_person_has_idle_animation && p.waiting_start > 0) ? frame_counter : 0);
				ped_bcubes.emplace_back(p.get_bcube(), ix);
			}
		}
		check_for_shadow_caster(ped_bcubes, light_bcube, lpos, dmax, has_stairs, xlate, shadow_caster_hash);
	}
	// update shadow_caster_hash for moving objects
	check_for_shadow_caster(moving_objs, light_bcube, lpos, dmax, has_stairs, xlate, shadow_caster_hash);
}
// rat_t/spider_t/snake_t
template<typename T> void get_animal_shadow_casters(vector<T> &animals, vect_cube_with_ix_t &moving_objs, vector3d const &xlate, unsigned base_ix) {
	for (auto i = animals.begin(); i != animals.end(); ++i) {
		if (!i->is_moving()) continue;
		// calculate which animals are visible to the player and consider those for shadow casters;
		// the shadow_non_visible flag will also be used to cull them in the shadow pass;
		// assumes non-vertical/hanging animals (such as spiders) may cast a shadow far below them
		cube_t const bcube(i->get_bcube());
		i->shadow_non_visible = (i->get_upv().z > 0.99 && !camera_pdu.cube_visible(bcube + xlate));
		if (!i->shadow_non_visible) {moving_objs.emplace_back(bcube, (i - animals.begin() + base_ix));} // only add if shadow visible
	} // for i
}
void expand_cube_zvals(cube_t &c, float z1, float z2) {
	min_eq(c.z1(), z1);
	max_eq(c.z2(), z2);
}

bool check_cube_visible_through_cut(vect_cube_t const &cuts, cube_t const &light_bounds, point const &lpos, point const &camera_bs, float light_radius, bool floor_is_above) {
	float const light_dist(p2p_dist(camera_bs, lpos)); // upper bound on line length

	for (cube_t const &s : cuts) { // check visible stairs/ramps
		if (!light_bounds.intersects_xy(s)) continue; // light does not reach the cut
		float const z(floor_is_above ? s.z2() : s.z1());
		point const pts[4] = {point(s.x1(), s.y1(), z), point(s.x2(), s.y1(), z), point(s.x2(), s.y2(), z), point(s.x1(), s.y2(), z)};

		for (unsigned n = 0; n < 4; ++n) { // check if the ray through any corner of the gap hits the light bcube
			vector3d const delta(pts[n] - camera_bs);
			float const corner_dist(delta.mag());
			vector3d const dir(delta/corner_dist);
			if (!line_intersect_sphere(camera_bs, dir, lpos, light_radius)) continue; // test bsphere
			float const ray_len(max(light_dist, 1.1f*corner_dist)); // make sure the ray projects past the corner of the cut to the floor above or below
			if (light_bounds.line_intersects(camera_bs, (camera_bs + dir*ray_len))) return 1; // test bounds
		}
	} // for s
	return 0;
}

bool add_dlight_if_visible(point const &pos, float radius, colorRGBA const &color, vector3d const &xlate, cube_t &lights_bcube, vector3d const &dir=zero_vector, float bwidth=1.0) {
	if (!lights_bcube.contains_pt_xy(pos)) return 0;
	if (!camera_pdu.sphere_visible_test((pos + xlate), radius)) return 0; // VFC
	dl_sources.emplace_back(radius, pos, pos, color, 0, dir, bwidth); // is_dynamic=0
	dl_sources.back().disable_shadows();
	expand_cube_zvals(lights_bcube, (pos.z - radius), (pos.z + radius));
	return 1;
}

// Note: non const because this caches light_bcubes
void building_t::add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building, bool sec_camera_mode,
	occlusion_checker_noncity_t &oc, vect_cube_with_ix_t &ped_bcubes, cube_t &lights_bcube)
{
	if (!has_room_geom()) return; // error?
	point const camera_bs(camera_pdu.pos - xlate), building_center(bcube.get_cube_center()); // camera in building space
	bool const player_by_ext_door(!camera_in_building && point_near_ext_door(camera_bs, get_door_open_dist()));
	bool walkway_only(0), same_floor_only(0), same_or_adj_floor_only(0);

	if (!camera_in_building && !has_windows()) { // can't see interior through windows
		bool above_skylight(0);

		for (cube_with_ix_t &sl : skylights) {
			if (camera_bs.z > sl.z2()) {above_skylight = 1; break;}
		}
		if (!above_skylight) {
			if (!player_by_ext_door) { // interior lights not visible
				if (check_pt_in_or_near_walkway(camera_bs, 1, 1, 1)) {walkway_only = 1;} // player in or near walkway
				else if (has_int_windows() && player_building != nullptr && is_connected_with_walkway(*player_building, camera_bs.z)) {walkway_only = 1;}
				else return; // no lights visible
			}
			// ground floor door may have stairs visible and includes floor above and below; walkway only includes the current floor
			((camera_bs.z < ground_floor_z1 + get_window_vspace()) ? same_or_adj_floor_only : same_floor_only) = 1;
		}
	}
	else if ((display_mode & 0x08) && !camera_in_building && !bcube.contains_pt_xy(camera_bs) && is_entire_building_occluded(camera_bs, oc)) return;
	// Note: camera_bs is used to test against bcube, lpos_rot, and anything else in global space; camera_rot is used to test against building interior objects
	point const camera_rot(get_inv_rot_pos(camera_bs)); // rotate camera into building space; use this pos below except with building bcube, occlusion checks, or lpos_rot
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), fc_thick(get_fc_thickness()), fc_gap(get_floor_ceil_gap());
	float const room_xy_expand(0.75*wall_thickness), player_feet_zval(camera_bs.z - get_bldg_player_height()), ground_floor_z2(ground_floor_z1 + window_vspacing);
	bool const check_building_people(enable_building_people_ai()), check_attic(camera_in_building && has_attic() && interior->attic_access_open);
	bool const camera_in_basement(camera_bs.z < ground_floor_z1), camera_in_ext_basement(camera_in_building && point_in_extended_basement_not_basement(camera_rot));
	bool const show_room_name(display_mode & 0x20); // debugging, key '6'
	cube_t const &attic_access(interior->attic_access);
	vect_cube_t &light_bcubes(interior->room_geom->light_bcubes);
	vect_room_object_t &objs(interior->room_geom->objs); // non-const, light flags are updated
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	// Note: we check (camera_bs.z > attic_access.z1()) rather than player_in_attic because we want to include the case where the player is in the access door
	bool const player_on_attic_stairs(check_attic && camera_bs.z > attic_access.z1() && attic_access.contains_pt_xy(camera_rot));
	bool const player_in_pool(camera_in_building && has_pool() && interior->pool.contains_pt(camera_rot));
	unsigned camera_part(parts.size()); // start at an invalid value
	bool camera_by_stairs(0), camera_on_stairs(0), camera_by_L_stairs(0), camera_somewhat_by_stairs(0), camera_in_hallway(0), camera_can_see_ext_basement(0);
	bool camera_near_building(camera_in_building), check_ramp(0), stairs_or_ramp_visible(0), camera_room_tall(0), camera_in_closed_room(0);
	float camera_z(camera_bs.z), up_light_zmin(camera_z);
	// if player is in the pool, increase camera zval to the top of the pool so that lights in the room above are within a single floor and not culled
	if (player_in_pool) {camera_z = interior->pool.z2();}
	int camera_room(-1), L_stairs_room(-1);
	vect_cube_t cuts_above, cuts_below, cuts_above_nonvis; // only used when player is in the building
	vect_cube_with_ix_t moving_objs;
	cube_t floor_above_region, floor_below_region; // filters for lights on the floors above/below based on stairs
	bool saw_open_stairs(0);
	ped_bcubes.clear();
	bool const track_lights(0 && camera_in_building && !sec_camera_mode && animate2); // used for debugging
	bool const camera_above_ground_floor(camera_z > ground_floor_z2), camera_feet_above_basement(player_feet_zval >= ground_floor_z1);
	if (track_lights) {enabled_bldg_lights.clear();}

	if (camera_in_building) {
		run_light_motion_detect_logic(camera_rot);

		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (i->contains_pt(camera_rot)) {camera_part = (i - parts.begin()); break;}
		}
		int const room_ix(get_room_containing_camera(camera_rot));

		if (room_ix >= 0) {
			// Note: stairs connecting stacked parts aren't flagged with has_stairs because stairs only connect to the bottom floor, but they're partially handled below
			room_t const &room(get_room(room_ix));
			unsigned const camera_floor(room.get_floor_containing_zval(max(camera_bs.z, room.z1()), window_vspacing)); // clamp zval to room range
			unsigned const room_type(room.get_room_type(camera_floor));
			assert(room_type < NUM_RTYPES);
			camera_room       = room_ix;
			camera_by_stairs  = camera_somewhat_by_stairs = room.has_stairs_on_floor(camera_floor);
			camera_in_hallway = room.is_hallway;
			check_ramp        = (has_pg_ramp() && !interior->ignore_ramp_placement);
			camera_room_tall  = (room.is_single_floor && camera_bs.z > room.z1() + window_vspacing);
			if (room.is_single_floor) {up_light_zmin = max(camera_bs.z-window_vspacing, room.z1());} // player can see upward lights on walls when in a tall ceiling room
			if (show_room_name) {lighting_update_text = room_names[room_type];}

			// stairs and ramps only allow light to pass if visible to the player
			float const zval(room.z1() + camera_floor*window_vspacing), ceil_below_z(zval - fc_thick), floor_below_z(zval + fc_thick);
			float const ceil_above_z(ceil_below_z + window_vspacing), floor_above_z(floor_below_z + window_vspacing);
			vect_cube_with_ix_t stair_ramp_cuts; // ix enables cuts {below, above}

			for (stairwell_t const &s : interior->stairwells) { // check stairs
				if (s.z1() > camera_z || s.z2() < camera_z) continue; // wrong floor
				
				if (s.is_l_shape()) {
					L_stairs_room = get_room_containing_pt(s.get_cube_center()); // there should be only one
					if (s.intersects(room)) {camera_by_L_stairs = 1;} // same room
					else if (L_stairs_room >= 0) { // adjacent connected rooms
						if (are_rooms_connected(room, get_room(L_stairs_room), (zval + 0.5*window_vspacing), 1)) {camera_by_L_stairs = 1;} // check_door_open=1
					}
				}
				unsigned cut_mask(3); // both dirs enabled by default

				if (!s.contains_pt(camera_rot)) { // disable these optimizations if the player is on the stairs
					float const center_val(s.get_center_dim(s.dim));
					bool const dir_val((camera_rot[s.dim] < center_val) ^ s.dir);
					
					if (s.is_u_shape()) { // light may be visible on the back wall of the stairs from the light above or on the edges of the stairs from the light below
						if (dir_val) continue; // back facing - light not visible
						// floors above and below are only visible from one side of stairs
						cube_t vis_region(bcube); // start with full building bcube
						vis_region.d[s.dim][s.dir] = center_val;
						floor_above_region.assign_or_union_with_cube(vis_region);
						floor_below_region.assign_or_union_with_cube(vis_region);
					}
					else if (s.has_walled_sides()) {
						if (!dir_val) { // facing updards stairs
							cut_mask = 2; // floor below not visible
							cube_t vis_region(bcube); // start with full building bcube
							vis_region.d[s.dim][!s.dir] = center_val;
							floor_above_region.assign_or_union_with_cube(vis_region); // floor above may be visible
						}
					}
					else {saw_open_stairs = 1;}
				}
				stair_ramp_cuts.emplace_back(s, cut_mask);
			} // for s
			// some parking garages may have open stairs connecting to the floor above; a light may be visible through these even if there are also walled or U-shaped stairs
			if (saw_open_stairs) {floor_above_region.set_to_zeros(); floor_below_region.set_to_zeros();}

			if (has_pg_ramp()) { // check ramp if player is in the parking garage or the backrooms doorway
				if (room.contains_pt(interior->pg_ramp.get_cube_center()) ||
					(has_ext_basement() && interior->get_ext_basement_door().get_clearance_bcube().contains_pt(camera_rot)))
				{
					stair_ramp_cuts.emplace_back(interior->pg_ramp, 3); // both dirs
					cube_t &region((camera_floor == 0) ? floor_above_region : floor_below_region); // lights on PG level above or below may be visible through ramp opening
					//region.set_to_zeros(); // conservative
					
					if (!region.is_all_zeros()) {
						bool const rdim(interior->pg_ramp.ix >> 1);
						cube_t ramp_ext(interior->pg_ramp);
						ramp_ext.expand_in_dim(rdim, interior->pg_ramp.get_sz_dim(rdim)); // expand to capture lights off the ends of the ramp
						region.union_with_cube_xy(ramp_ext); // close enough?
					}
				}
			}
			for (cube_with_ix_t const &s : stair_ramp_cuts) {
				if (s.contains_pt(camera_rot)) {camera_on_stairs = stairs_or_ramp_visible = camera_by_stairs = camera_somewhat_by_stairs = 1;} // player on this stairs or ramp

				if ((s.ix & 2) && s.z2() >= floor_above_z) { // cut above
					cube_t cut(s);
					set_cube_zvals(cut, ceil_above_z, floor_above_z);
					bool const visible(is_rot_cube_visible(cut, xlate) && !((display_mode & 0x08) && check_obj_occluded(cut, camera_bs, oc, 0))); // VFC + occlusion culling
					(visible ? cuts_above : cuts_above_nonvis).push_back(cut);
					stairs_or_ramp_visible |= visible;
				}
				if ((s.ix & 1) && s.z1() <= ceil_below_z) { // cut below
					cube_t cut(s);
					set_cube_zvals(cut, ceil_below_z, floor_below_z);
					bool const visible(is_rot_cube_visible(cut, xlate) && !((display_mode & 0x08) && check_obj_occluded(cut, camera_bs, oc, 0))); // VFC + occlusion culling
					if (visible) {cuts_below.push_back(cut);}
					stairs_or_ramp_visible |= visible;
				}
			} // for s
			// set camera_somewhat_by_stairs when camera is in room with stairs, or adjacent to one with stairs
			if (stairs_or_ramp_visible) {camera_somewhat_by_stairs |= bool(room_or_adj_room_has_stairs(camera_room, camera_bs.z, 1, 1));} // inc_adj_rooms=1, check_door_open=1
			// if player is by the stairs in a room with all closed doors, it's still possible to see a light shining through a door of the floor above or below
			camera_in_closed_room = (!camera_by_stairs && !player_on_attic_stairs && all_room_int_doors_closed(room_ix, camera_z));
		}
		else if (point_in_attic(camera_rot)) {
			if (show_room_name) {lighting_update_text = room_names[RTYPE_ATTIC];}
		}
		//lighting_update_text = ((is_sphere_lit(camera_rot, get_scaled_player_radius()) || is_sphere_lit((camera_rot - vector3d(0.0, 0.0, get_player_height())), get_scaled_player_radius())) ? "Lit" : "Unlit");
	}
	else {
		cube_t bcube_exp(bcube);
		bcube_exp.expand_by_xy(2.0*window_vspacing);
		camera_can_see_ext_basement = interior_visible_from_other_building_ext_basement(xlate, 1); // expand_for_light=1
		camera_near_building        = bcube_exp.contains_pt(camera_bs) || camera_can_see_ext_basement;
	}
	if (camera_near_building) { // build moving objects vector
		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->is_moving() && i->is_visible()) {moving_objs.emplace_back(*i, (i - objs.begin() + 1));}
		}
	}
	if (has_room_geom() && frame_counter <= (int)interior->room_geom->last_animal_update_frame+1) { // animals were updated this frame or the previous frame
		if (has_retail() && get_retail_part().contains_pt(camera_rot)) {} // optimization: no dynamic animal shadows in retail area
		else { // add a base index to each animal group to make all moving objects unique
			get_animal_shadow_casters(interior->room_geom->rats,    moving_objs, xlate, 10000);
			get_animal_shadow_casters(interior->room_geom->spiders, moving_objs, xlate, 20000);
			get_animal_shadow_casters(interior->room_geom->snakes,  moving_objs, xlate, 30000);
		}
	}
	//highres_timer_t timer("Lighting", camera_in_building); // 13.8ms => 13.1ms => 12.7ms => 3.6ms => 3.2ms => 0.81
	//unsigned num_add(0);
	unsigned last_room_ix(interior->rooms.size()); // start at an invalid value
	bool last_room_closed(0);

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (camera_near_building && !walkway_only) { // handle light emitting objects in the player's building
			bool was_added(0);

			if (i->type == TYPE_LAVALAMP && i->is_light_on()) { // should we do an occlusion query?
				if (camera_room_tall && camera_room == (int)i->room_id) {} // special case for player in a tall room - skip the continue case below
				else if (int((camera_bs.z - ground_floor_z1)/window_vspacing) != int((i->zc() - ground_floor_z1)/window_vspacing)) continue; // different floor
				if (!add_dlight_if_visible(i->get_cube_center(), 10.0*i->get_radius(), colorRGBA(1.0, 0.75, 0.25), xlate, lights_bcube)) continue;
				was_added = 1;
			}
			else if (i->type == TYPE_FISHTANK && i->is_light_on()) { // should we do an occlusion query?
				if (camera_room_tall && camera_room == (int)i->room_id) {} // special case for player in a tall room - skip the continue case below
				else if (int((camera_bs.z - ground_floor_z1)/window_vspacing) != int((i->zc() - ground_floor_z1)/window_vspacing)) continue; // different floor
				if (!add_dlight_if_visible(cube_top_center(*i), 1.25*(i->dx() + i->dy()), colorRGBA(0.8, 0.9, 1.0), xlate, lights_bcube, -plus_z, 0.3)) continue; // pointed downward
				was_added = 1;
			}
			if (was_added) {
				cube_t room(get_room(i->room_id)); // should we clip to the current floor in Z?
				assert(room.contains_pt(dl_sources.back().get_pos()));
				room.expand_by_xy(wall_thickness);
				dl_sources.back().set_custom_bcube(room);
				continue;
			}
		}
		if (!i->is_light_on() || !i->is_light_type()) continue; // light not on, or not a light or lamp
		point lpos(i->get_cube_center()); // centered in the light fixture
		min_eq(lpos.z, (i->z2() - 0.0125f*window_vspacing)); // make sure the light isn't too close to the ceiling (if shifted up to avoid a door intersection)
		point lpos_rot(lpos); // lpos in global space
		if (is_rotated()) {do_xy_rotate(building_center, lpos_rot);}
		if (!lights_bcube.contains_pt_xy(lpos_rot)) continue; // not contained within the light volume
		bool const is_in_attic(i->in_attic()), is_in_windowless_attic(is_in_attic && !has_attic_window), is_exterior(i->is_exterior());
		if (walkway_only && !is_exterior) continue;
		bool const light_in_basement(lpos.z < ground_floor_z1), is_in_elevator(i->in_elevator()), is_in_closet(i->in_closet());
		// basement, attic, and elevator lights are only visible when player is in the building;
		// elevator test is questionable because it can be open on the ground floor of a room with windows in a small office building, but should be good enough
		if (!camera_in_building && ((light_in_basement && !camera_can_see_ext_basement) || is_in_windowless_attic || is_in_elevator)) continue;
		if ((is_in_elevator || is_in_closet) && camera_z > lpos.z) continue; // elevator or closet light on the floor below the player
		if (light_in_basement && camera_above_ground_floor)        continue; // basement lights only visible if player is on basement or ground floor
		room_t const &room(get_room(i->room_id));
		bool const in_ext_basement(room.is_ext_basement()), in_retail_room(room.is_retail());
		bool light_vis_from_basement(light_in_basement);

		if (camera_in_ext_basement && !light_in_basement && lpos.z < ground_floor_z2) { // check if light by basement stairs; can happen with office ground floor hall lights
			for (stairwell_t const &s : interior->stairwells) {
				if (s.z1() >= ground_floor_z1 || s.z2() < ground_floor_z1) continue; // not basement stairs
				if (s.is_u_shape()) continue; // light from above U-shaped stairs can't be seen from the bottom
				cube_t stairs_exp(s);
				stairs_exp.expand_in_dim(s.dim, 2.0*window_vspacing);
				if (stairs_exp.intersects_xy(*i)) {light_vis_from_basement = 1; break;}
			}
		}
		if (in_ext_basement && camera_feet_above_basement)      continue; // light  in extended basement, and camera not in basement or on basement stairs
		if (camera_in_ext_basement && !light_vis_from_basement) continue; // camera in extended basement, and light  not in basement
		
		// if the player is in the parking garage, retail lights may not be visible through stairs; this check may not be valid if ramps extend to the first floor in the future
		if (camera_in_basement && !camera_on_stairs && in_retail_room && !interior->stairwells.empty()) {
			stairwell_t const &s(interior->stairwells.front()); // assume first stairs placed are the ones connecting retail to parking garage, and no other stairs
			if (s.is_u_shape()) continue; // U-shaped stairs - light is never visible
			float const entrance_plane(s.d[s.dim][!s.dir]);
			if ((camera_bs[s.dim] < entrance_plane) != s.dir) continue; // back facing
			cube_t opening(s);
			set_wall_width(opening, entrance_plane, wall_thickness, s.dim); // shrink to a narrow strip
			set_cube_zvals(opening, room.z1()-window_vspacing, room.z1()); // limit to upper floor of parking garage
			if (!is_rot_cube_visible(opening, xlate)) continue; // stairs opening not visible; can't test upper floor opening because walls may be lit
		}
		//if (is_light_occluded(lpos_rot, camera_bs))  continue; // too strong a test in general, but may be useful for selecting high importance lights
		//if (!camera_in_building && i->is_interior()) continue; // skip interior lights when camera is outside the building: makes little difference, not worth the trouble
		bool const is_lamp(i->type == TYPE_LAMP), is_single_floor(room.is_single_floor || is_in_elevator), wall_light(i->flags & RO_FLAG_ADJ_HI);
		int const cur_floor(is_single_floor ? 0 : (i->z1() - room.z1())/window_vspacing); // garages and sheds are all one floor
		float const level_z(is_in_attic ? interior->attic_access.z1() : (room.z1() + cur_floor*window_vspacing)), floor_z(level_z + fc_thick);
		float ceil_z(0.0);
		if      (is_in_attic         ) {ceil_z = interior_z2;} // top of interior/attic
		else if (is_in_elevator      ) {ceil_z = get_elevator(i->obj_id).z2();} // top of elevator shaft
		else if (room.is_single_floor) {ceil_z = room.z2();} // top of current room/part (garage shed, etc.)
		else                           {ceil_z = (level_z + window_vspacing - fc_thick);} // normal room light
		float const floor_below_zval(floor_z - window_vspacing), ceil_above_zval(ceil_z + window_vspacing);
		// Note: we use level_z rather than floor_z for floor_is_above test so that it agrees with the threshold logic for player_in_basement
		bool const floor_is_above((camera_z < level_z) && !is_single_floor), floor_is_below(camera_z > ceil_z+fc_thick); // check floor_ix transition points
		if (same_floor_only && (floor_is_above || floor_is_below)) continue;
		if (same_or_adj_floor_only && (camera_z < level_z-window_vspacing || camera_z > ceil_z+fc_thick+window_vspacing)) continue;
		if (!is_house && floor_is_below &&  in_ext_basement && !camera_in_ext_basement && lpos.z   < get_basement().z1()) continue; // light  in lower level extb, player not in extb
		if (!is_house && floor_is_above && !in_ext_basement &&  camera_in_ext_basement && camera_z < get_basement().z1()) continue; // player in lower level extb, light  not in extb
		if (!camera_on_stairs && floor_is_below && !floor_below_region.is_all_zeros() && !floor_below_region.contains_pt_xy(lpos)) continue; // check floor_below_region
		if (!camera_on_stairs && floor_is_above && !floor_above_region.is_all_zeros() && !floor_above_region.contains_pt_xy(lpos) && !(i->flags & RO_FLAG_TOS)) continue; // check floor_below_region
		if (is_exterior && (floor_is_below || floor_is_above)) continue; // different floor of walkway - not visible
		if (player_by_ext_door && camera_z < level_z - window_vspacing) continue; // can't see more than one floor above current floor through open door
		cube_t const &room_part(get_part_for_room(room));
		bool const camera_in_room_part_xy(room_part.contains_pt_xy(camera_rot)), in_camera_room((int)i->room_id == camera_room);
		bool const camera_room_same_part(room.part_id == camera_part || (is_house && camera_in_room_part_xy)); // treat stacked house parts as the same
		bool const has_stairs_this_floor(!is_in_attic && room.has_stairs_on_floor(cur_floor));
		bool const light_room_has_stairs_or_ramp(i->has_stairs() || has_stairs_this_floor || (check_ramp && is_room_above_ramp(room, i->z1())));
		bool const is_over_pool(has_pool() && (int)i->room_id == interior->pool.room_ix);
		bool const light_and_camera_by_L_stairs((in_camera_room || (int)i->room_id == L_stairs_room) && has_stairs_this_floor && camera_by_L_stairs);
		//bool const light_room_is_tall(room.is_single_floor && lpos.z > room.z1() + window_vspacing);
		// special case for light shining down from above stairs or ramp when the player is below
		bool const light_above_stairs(lpos.z > camera_z && light_room_has_stairs_or_ramp); // or ramp
		bool stairs_light(0), camera_in_elevator(0), cull_if_not_by_stairs(0), in_camera_walkway(0), in_walkway_near_camera(0);
		//if (!light_above_stairs && camera_in_basement && !light_in_basement) {cull_if_not_by_stairs = 1;} // light may not be visible from basement; check not needed?

		if (is_in_elevator) {
			elevator_t const &e(get_elevator(i->obj_id));
			room_object_t const &car(interior->get_elevator_car(e)); // elevator car for this elevator
			
			if (car.contains_pt(camera_rot)) {camera_in_elevator = 1;} // player inside elevator
			else if (e.open_amt == 0.0) { // closed elevator
				if (floor_is_above || floor_is_below)          continue; // viewed from a different floor
				if ((lpos[e.dim] < camera_rot[e.dim]) ^ e.dir) continue; // camera facing the back of the elevator: light not visible
			}
		}
		if (camera_in_building && !is_in_elevator && !is_in_attic && !floor_is_above && !floor_is_below) {
			if (i->room_id != last_room_ix) { // new room
				last_room_closed = all_room_int_doors_closed(i->room_id, lpos.z);
				last_room_ix     = i->room_id;
			}
			// if either the camera or the light are in different rooms with closed doors,
			// on the same floor (not separated by stairs) of the same part (not visible across windows), then the light isn't visible
			if ((last_room_closed || camera_in_closed_room) && i->room_id != camera_room && (room.part_id == camera_part || !has_windows())) continue;
		}
		if (!camera_in_elevator) { // none of the below culling applies when the player is in the elevator
			// if the light is in the basement and the camera isn't, it's not visible unless the player is by the stairs
			if ( light_in_basement && player_in_basement == 0 && !camera_somewhat_by_stairs) continue;
			//if (!light_vis_from_basement && player_in_basement >= 2 && !light_room_has_stairs_or_ramp) continue; // player is fully in basement but light isn't; too strong
			if (!light_vis_from_basement && player_in_basement >= 3) continue; // player is in the extended basement but the light isn't in the basement
			bool const check_stairs(stairs_or_ramp_visible || light_above_stairs);
			bool const camera_within_one_floor(camera_z > floor_below_zval && camera_z < ceil_above_zval);
			
			// less culling if either the light or the camera is by stairs and light is on the floor above or below
			if (camera_within_one_floor && check_stairs) { // light is on the floor above or below the camera
				if (camera_in_hallway && camera_by_stairs && camera_room_same_part) {
					// special case for player in an office building primary hallway with stairs; only handle the case where the light is in the hallway above or below;
					// if camera is on the stairs or a ramp this also counts because this may be connecting two rooms in two different parts
					stairs_light = (in_camera_room || camera_on_stairs);
				}
				else {
					stairs_light = (light_room_has_stairs_or_ramp || camera_somewhat_by_stairs); // either the light or the camera is by the stairs
				}
			}
			// include lights above or near the parking garage ramp
			stairs_light |= (light_in_basement && has_pg_ramp() && interior->pg_ramp.contains_pt_xy_exp(i->get_cube_center(), window_vspacing));

			if (player_in_pool && is_over_pool) {
				// camera and light both in pool room - keep it
			}
			else if (check_attic && floor_is_below && camera_bs.z > attic_access.z1() && room.contains_cube_xy(attic_access)) {
				// camera in attic, possibly looking down through attic access door, and light is in the room below - keep it
			}
			else if (check_attic && floor_is_above && lpos.z > attic_access.z2() && camera_room >= 0 && get_room(camera_room).contains_cube_xy(attic_access)) {
				// light in attic, and camera in room with attic access
			}
			else if (light_and_camera_by_L_stairs) {
				// lights may be visible multiple rooms above or below through the hole in the center of L-shaped stairs
			}
			else if (player_on_attic_stairs && floor_is_below) {
				// lights on floor below may be visible through attic opening
			}
			else if (floor_is_above || (floor_is_below && !camera_room_tall)) { // light is on a different floor from the camera
				bool const parts_are_stacked(camera_part < real_num_parts && (parts[camera_part].z2() <= room.z1() || parts[camera_part].z1() >= room.z2()));

				// the basement is a different part, but it's still the same vertical stack; consider this the same effective part if the camera is in the basement above the room's part
				if (camera_in_ext_basement && !in_camera_room && has_stairs_this_floor && in_ext_basement &&
					(camera_by_stairs || (light_room_has_stairs_or_ramp && camera_somewhat_by_stairs)))
				{
					// camera and light are on different floors of different rooms of the extended basement in two rooms connected by stairs
				}
				else if (camera_in_building && (camera_room_same_part // camera and light in same part; can't see through a window, can only see through ceiling/floor/stairs
					|| (player_in_basement || light_in_basement) // basement only visible through ceiling/floor/stairs
					|| (parts_are_stacked && camera_within_one_floor))) // stacked parts maybe connected with stairs; shouldn't be visible through windows
				{
					// player is on a different floor of the same building part, or more than one floor away in a part stack, and can't see a light from the floor above/below
					if (!stairs_light) continue; // camera in building and on wrong floor, don't add light; will always return if more than one floor away
					cull_if_not_by_stairs = 1;
				}
				else { // camera outside the building (or the part that contains this light)
					float const xy_dist(p2p_dist_xy(camera_bs, lpos_rot));
					if (!stairs_light && ((camera_z - lpos_rot.z) > 2.0f*xy_dist || (lpos_rot.z - camera_z) > 1.0f*xy_dist)) continue; // light viewed at too high an angle

					if (camera_in_building) { // camera and light are in different buildings/parts
						if (camera_part >= real_num_parts) continue; // camera in garage or shed
						if (parts_are_stacked)             continue; // light in a different vertical stack than the camera; not visible through windows

						if (!is_rotated()) { // check exterior wall visibility; this part doesn't work for rotated buildings
							// is it better to check if light half sphere is occluded by the floor above/below?
							bool visible[2] = {0};

							for (unsigned d = 0; d < 2; ++d) { // for each dim
								bool const dir(camera_bs[d] > lpos_rot[d]);
								// if we're not on the outside face of the part containing this room, we can't see through any windows
								if ((camera_rot[d] > room_part.d[d][dir]) ^ dir) continue;
								visible[d] = (room.ext_sides & (1 << (2*d + dir)));
								// if we're by the stairs, a light in the other part may be visible through a door in the wall separating the rooms/parts
								if (camera_somewhat_by_stairs) {visible[d] |= (fabs(room.d[d][dir] - parts[camera_part].d[d][!dir]) < 2.0*wall_thickness);}
							}
							if (!visible[0] && !visible[1]) continue; // room is not on the exterior of the building on either side facing the camera
						}
					}
				}
			} // end camera on different floor case
		} // end !camera_in_elevator
		float const light_radius(get_radius_for_room_light(*i)), cull_radius(0.95*light_radius);
		// note that the same lights are used for the reflection pass, so a light behind the player won't be active in a mirror reflection
		if (!camera_pdu.sphere_visible_test((lpos_rot + xlate), cull_radius)) continue; // VFC
		// ext basement connector room must include the other building's ext basement, and it's simplest to just expand it by the max length of that room plus approx hallway width
		bool const is_ext_conn_light(i->is_exterior());
		float const light_bcube_expand(is_ext_conn_light ? (EXT_BASEMENT_JOIN_DIST*window_vspacing + get_doorway_width()) : 0.0);
		// check visibility of bcube of light sphere clipped to building bcube; this excludes lights behind the camera and improves shadow map assignment quality
		cube_t sphere_bc, light_clip_cube; // in building space, unrotated
		sphere_bc.set_from_sphere(lpos, cull_radius);
		bool light_in_walkway(0);

		if (light_in_basement) { // clip to basement + ext basement
			light_clip_cube = get_basement();
			if (has_ext_basement()) {light_clip_cube.union_with_cube(interior->basement_ext_bcube);}
			if (is_in_elevator)     {light_clip_cube.z2() = bcube.z2();} // extends up the elevator shaft into floors above the basement
			assert(light_clip_cube.contains_pt(lpos)); // Note: may not be contained in building bcube
			light_clip_cube.expand_by_xy(light_bcube_expand);
		}
		else if (is_exterior) { // exterior lights are only in walkways and are only visible when the player is inside or in a connected room
			assert(!walkways.empty());

			for (building_walkway_t const &w : walkways) { // Note: walkways shouldn't be added to rotated buildings
				if (!w.is_owner || !w.bcube.contains_pt(lpos)) continue;
				if (w.bcube_inc_rooms.contains_pt(camera_rot)) {in_camera_walkway = in_walkway_near_camera = 1;} // camera in or near the walkway
				else if (w.has_skyway_conn() && w.skyway_conn.contains_pt(camera_rot)) {in_walkway_near_camera = 1;} // camera in skyway near the walkway
				light_clip_cube  = w.bcube;
				light_in_walkway = 1;
				break;
			}
			assert(light_in_walkway);
		}
		else if (room.is_sec_bldg) { // secondary buildings only light their single room and the exterior door
			light_clip_cube = room;
			light_clip_cube.expand_by_xy(wall_thickness); // needed for garage door
		}
		else { // clip to bcube
			if (is_rotated()) {light_clip_cube = get_rotated_bcube(bcube, 1);} // inv_rotate=1
			else {light_clip_cube = bcube;}
			assert(light_clip_cube.contains_pt(lpos));
		}
		cube_t clipped_bc(sphere_bc);
		clipped_bc.intersect_with_cube(light_clip_cube);
		// clip zval to current floor if light not in a room with stairs, elevator, pool, tall ceilings, or L-shaped stairs
		if (!stairs_light && !is_single_floor && !is_over_pool && !light_and_camera_by_L_stairs) {max_eq(clipped_bc.z1(), (floor_z - fc_thick));}
		float const z2_adj((wall_light ? 1.1 : 1.0)*fc_thick); // prevent Z-fighting on stairs lights
		min_eq(clipped_bc.z2(), (ceil_z + z2_adj)); // ceiling is always valid, since lights point downward

		if (!clipped_bc.is_strictly_normalized()) {
			cout << "Error: Invalid light bcube: " << TXT(clipped_bc.str()) << TXT(sphere_bc.str()) << TXT(light_clip_cube.str())
				 << TXT(lpos.str()) << TXT(floor_z) << TXT(ceil_z) << TXT(is_lamp) << TXT(is_in_elevator) << TXT(light_in_basement) << endl;
			assert(0);
		}
		if (!is_rot_cube_visible(clipped_bc, xlate, 1)) continue; // VFC; inc_mirror_reflections=1
		bool recheck_coll(0);

		if (cull_if_not_by_stairs) { // test light visibility through stairs and ramp cuts
			vect_cube_t const &cuts(floor_is_above ? cuts_above : cuts_below);
			point test_pt(lpos);
			if (floor_is_below) {test_pt.z = floor_z;} // if light is below us, pick a point on the floor below it as this is more likely to be visible
			bool maybe_visible(0);

			for (cube_t const &s : cuts) { // check visible stairs/ramps
				if (s.line_intersects(camera_rot, test_pt)) {maybe_visible = 1; break;} // center visible through gap
			}
			if (!maybe_visible && check_cube_visible_through_cut(cuts, clipped_bc, lpos, camera_rot, cull_radius, floor_is_above)) {maybe_visible = recheck_coll = 1;}

			if (!maybe_visible && light_above_stairs) { // check for light shining down the stairs behind or out of view of the player
				for (cube_t const &s : cuts_above_nonvis) { // compute bcube of light rays through gap to the floor below
					if (!dist_less_than(s.get_cube_center(), lpos, cull_radius)) continue; // light not near this cut
					float const z(s.z2()), t((lpos.z - floor_below_zval)/(lpos.z - floor_z)); // slightly greater than 2.0
					point const pts[4] = {point(s.x1(), s.y1(), z), point(s.x2(), s.y1(), z), point(s.x2(), s.y2(), z), point(s.x1(), s.y2(), z)};
					cube_t vis_bcube(s);
					for (unsigned n = 0; n < 4; ++n) {vis_bcube.union_with_pt(lpos + t*(pts[n] - lpos));}
					if (is_rot_cube_visible(vis_bcube, xlate) && !((display_mode & 0x08) && check_obj_occluded(vis_bcube, camera_bs, oc, 0))) {maybe_visible = 1; break;}
				} // for s
			}
			if (!maybe_visible) continue;
		}
		//if (line_intersect_walls(lpos, camera_rot)) continue; // straight line visibility test - for debugging, or maybe future use in assigning priorities
		//if (check_cube_occluded(clipped_bc, interior->fc_occluders, camera_rot)) continue; // legal, but may not help much
		bool const is_fully_broken(i->is_broken2());
		
		// run flicker logic for broken lights; this is done later in the control flow because updating light geometry can be expensive
		if (i->is_broken() || is_fully_broken) {
			static rand_gen_t rgen;

			if (animate2 && tfticks > i->light_amt) { // flickering; time for state transition
				float const flicker_time_secs(rgen.rand_uniform(0.1, 1.0));
				float const on_amt(is_fully_broken ? 0.02 : 1.0), off_amt(is_fully_broken ? 1.0 : 0.1);
				float const delay_mult(i->is_open() ? off_amt : on_amt);
				i->light_amt = tfticks + flicker_time_secs*TICKS_PER_SECOND*delay_mult; // schedule time for next transition
				i->flags    ^= RO_FLAG_OPEN;
				// regenerate lights geometry (can be somewhat slow); only update if player is below the level of the light and the light itself is visible
				if (camera_bs.z < i->z2() && is_rot_cube_visible(*i, xlate) && !((display_mode & 0x08) && check_obj_occluded(*i, camera_bs, oc, 0))) {
					interior->room_geom->invalidate_lights_geom();
				}
			}
			// emit sparks only if player in building, which should be true since we shouldn't get here otherwise
			if (animate2 && is_fully_broken && camera_in_building && rgen.rand_float() < 22*fticks/TICKS_PER_SECOND) { // 22/s
				cube_t sparks_area(*i);
				sparks_area.z1() = floor_z;
				sparks_area.expand_by_xy(0.75*window_vspacing);

				if (is_rot_cube_visible(sparks_area, xlate)) {
					interior->room_geom->particle_manager.add_for_obj(*i, 0.16*i->dz(), -plus_z, 0.001, 1, 1, PART_EFFECT_SPARK, (i - objs.begin()));
				}
			}
			if (!i->is_open()) continue; // not currently on
		}
		// update lights_bcube and add light(s)
		if (is_lamp) { // lamps are generally against a wall and not in a room with stairs and only illuminate that one room
			expand_cube_zvals(lights_bcube, (lpos_rot.z - min(window_vspacing, light_radius)), (lpos_rot.z + min(window_vspacing, light_radius)));
		}
		else {
			float const plus_z_factor(wall_light ? 1.0 : 0.1); // if not a wall light it's pointed down - don't extend as far up
			expand_cube_zvals(lights_bcube, (lpos_rot.z - light_radius), (lpos_rot.z + plus_z_factor*light_radius));
		}
		colorRGBA color;
		unsigned shadow_caster_hash(0);
		bool lamp_was_moved(0), was_refined(0);

		if (is_lamp) { // no light refinement, since lamps are not aligned between floors; refinement doesn't help as much with houses anyway
			if (i->obj_id == 0) { // this lamp has not yet been assigned a light bcube (ID 0 will never be valid because the bedroom will have a light assigned first)
				assert(!light_bcubes.empty());
				i->obj_id = (uint16_t)light_bcubes.size();
				light_bcubes.push_back(cube_t()); // allocate a new entry
			}
			color = LAMP_COLOR; // soft white
		}
		else {
			color = i->get_color()*1.1; // make it extra bright
		}
		if (is_in_elevator) {
			elevator_t const &e(get_elevator(i->obj_id));
			room_object_t const &car(interior->get_elevator_car(e)); // elevator car for this elevator
			assert(car.contains_pt(lpos));
			cube_t clip_cube(car); // light is constrained to the elevator car
			clip_cube.expand_in_dim(!e.dim, 0.1*room_xy_expand); // expand sides to include walls adjacent to elevator (enough to account for FP error)
			// allow light to extend outside open elevator door; full light radius if closed, a small amount to lit the interior edge of each floor when closed
			float const light_extend((e.open_amt > 0.0) ? light_radius : 0.1*room_xy_expand);
			clip_cube.d[e.dim][e.dir] += (e.dir ? 1.0 : -1.0)*light_extend;
			clipped_bc.intersect_with_cube(clip_cube); // Note: clipped_bc is likely contained in clip_cube and could be replaced with it
			if (e.may_be_moving()) {hash_mix_point(e.get_llc(), shadow_caster_hash);} // make sure to update shadows if elevator or its doors are potentially moving
			if (e.open_amt > 0.0 && e.open_amt < 1.0) {shadow_caster_hash += hash_by_bytes<float>()(e.open_amt);} // update shadows if door is opening or closing
		}
		else {
			if (is_in_attic || is_exterior || room.is_sec_bldg) {} // nothing else to do
			else {
				assert(i->obj_id < light_bcubes.size());
				cube_t &light_bcube(light_bcubes[i->obj_id]);

				if (is_lamp && i->was_moved()) {
					i->flags &= ~RO_FLAG_MOVED; // clear moved flag, since we noticed it was moved
					light_bcube.set_to_zeros(); // will be recalculated below
					lamp_was_moved = 1; // force update of upward pointing light shadows
				}
				if (light_bcube.is_all_zeros()) { // not yet calculated - calculate and cache
					light_bcube = clipped_bc;
					bool const is_parking_garage(light_in_basement && has_parking_garage && !in_ext_basement);
					refine_light_bcube(lpos, light_radius, room, light_bcube, is_parking_garage); // incorrect for rotated buildings?
				}
				was_refined     = (light_bcube.get_area_xy() < clipped_bc.get_area_xy());
				clipped_bc.x1() = light_bcube.x1(); clipped_bc.x2() = light_bcube.x2(); // copy X/Y but keep orig zvals
				clipped_bc.y1() = light_bcube.y1(); clipped_bc.y2() = light_bcube.y2();
				clipped_bc.expand_by_xy(light_bcube_expand);

				if (!is_over_pool) {
					// clip to the current floor zval if there are no floor stairs or ramp cutouts for the light to pass through;
					// unclear if this helps with framerate or light culling/shadow map issues, but it seems like a good idea; possibly redundant with z1 clip above
					float clip_z1(is_single_floor ? room.z1() : (lpos.z - fc_gap));
					cube_t test_cube(clipped_bc);
					test_cube.z1() = clip_z1;

					for (stairwell_t const &s : interior->stairwells) {
						if (s.intersects(test_cube)) {min_eq(clip_z1, s.z1());}
					}
					if (has_pg_ramp() && interior->pg_ramp.intersects(test_cube)) {min_eq(clip_z1, interior->pg_ramp.z1());}
					max_eq(clipped_bc.z1(), (clip_z1 - fc_thick)); // slightly lower to include the floor
				}
			}
			// expand so that offset exterior doors are properly handled, but less for walkway lights
			bool const is_upper_floor(!room.is_single_floor && lpos.z > ground_floor_z1 + window_vspacing);
			clipped_bc.expand_by_xy((light_in_walkway ? 0.1 : (is_upper_floor ? 0.65 : 1.0))*room_xy_expand);
			clipped_bc.intersect_with_cube(sphere_bc); // clip to original light sphere, which still applies (only need to expand at building exterior)
		}
		if (!clipped_bc.contains_pt(lpos)) {
			static bool had_invalid_light_bcube_warning = 0;

			if (!had_invalid_light_bcube_warning) { // only print once
				cout << "Error: Invalid light bcube: " << TXT(clipped_bc.str()) << TXT(lpos.str()) << TXT(room.str()) << TXT(bcube.str()) << TXT(is_lamp) << TXT(is_in_elevator) << endl;
				had_invalid_light_bcube_warning = 1;
			}
			continue; // can fail in rare cases when very far from the origin, likely due to FP error, so skip light in this case
		}
		if (recheck_coll && was_refined) { // test on refined bcube
			vect_cube_t const &cuts(floor_is_above ? cuts_above : cuts_below);
			if (!check_cube_visible_through_cut(cuts, clipped_bc, lpos_rot, camera_bs, cull_radius, floor_is_above)) continue;
		}
		bool maybe_walkway(0);

		// handle light hitting open office building doors by expanding outward
		if (!is_house && !light_in_basement && !light_in_walkway && !doors.empty()) {
			// if this light is in a room connected to a walkway door, use the part containing the room rather than the building bcube;
			// that way a walkway connnecting to a recessed door (part edge inside bcube) will have a door that's properly lit
			bool const ground_floor(lpos.z < ground_floor_z1 + window_vspacing);
			maybe_walkway = (!ground_floor && check_pt_in_or_near_walkway(lpos, 0, 0, 1) && is_room_adjacent_to_ext_door(room));
			cube_t const &test_cube(maybe_walkway ? parts[room.part_id] : bcube);

			for (unsigned d = 0; d < 2; ++d) {
				if (clipped_bc.d[d][0] < test_cube.d[d][0]) {clipped_bc.d[d][0] = sphere_bc.d[d][0];}
				if (clipped_bc.d[d][1] > test_cube.d[d][1]) {clipped_bc.d[d][1] = sphere_bc.d[d][1];}
			}
		}
		if (!is_rot_cube_visible(clipped_bc, xlate, 1)) continue; // VFC - post clip; inc_mirror_reflections=1
		if ((display_mode & 0x08) && !clipped_bc.contains_pt(camera_rot) && check_obj_occluded(clipped_bc, camera_bs, oc, 0)) continue; // occlusion culling (expensive)
		//float const bwidth(wall_light ? 1.0 : 0.25); // wall light omnidirectional; ceiling light as close to 180 degree FOV as can get without shadow clipping; shadows are wrong
		float const bwidth(0.25); // as close to 180 degree FOV as we can get without shadow clipping
		// should bwidth be set smaller for (in_retail_room && has_tall_retail())?
		vector3d dir;
		if (wall_light) {dir[i->dim] = (i->dir ? 1.0 : -1.0);} else {dir = -plus_z;} // points down, unless it's a wall light
		dl_sources.emplace_back(light_radius, lpos_rot, lpos_rot, color, 0, dir, bwidth);
		if (track_lights) {enabled_bldg_lights.push_back(lpos_rot);}
		//++num_add;
		// use smaller shadow radius for retail rooms, since there are so many lights (meaning shadows are less visible and perf is more important)
		float const light_radius_shadow((in_retail_room ? RETAIL_SMAP_DSCALE : 1.0)*light_radius);
		bool force_smap_update(0);

		// check for dynamic shadows; check the player first; use full light_radius_shadow
		if (camera_surf_collide && (camera_in_building || in_camera_walkway || (player_in_walkway && maybe_walkway) || camera_can_see_ext_basement) &&
			dist_less_than(lpos_rot, camera_bs, light_radius_shadow))
		{
			bool player_in_this_room(player_on_attic_stairs && (is_in_attic || room.intersects_xy(interior->attic_access))); // ladder case
			player_in_this_room |= (player_in_pool && is_over_pool); // pool case
			player_in_this_room |= in_camera_walkway; // walkway case

			if (clipped_bc.contains_pt(camera_rot) || clipped_bc.contains_pt(point(camera_rot.x, camera_rot.y, player_feet_zval)) || player_in_this_room) {
				// must update shadow maps for the room above if the player is on the stairs or in the same room when there are stairs
				// must update even if light is visible (meaning shadows aren't due to < 90 degree FOV) because the old shadow may become visible when the player moves
				bool const check_floor_above(camera_on_stairs || (camera_by_stairs && in_camera_room) || is_single_floor);

				if (is_lamp || player_in_this_room || (player_in_attic && is_in_attic) ||
					(lpos_rot.z > player_feet_zval && (check_floor_above || lpos_rot.z < (camera_bs.z + window_vspacing))))
				{
					//if (dl_sources.back().calc_pdu(0, 0, 0.0).sphere_visible_test(camera_bs, CAMERA_RADIUS)) {} // not really correct, and doesn't help much
					// player shadow, based on head to feet Z-range; includes lamps (with no zval test)
					force_smap_update   = 1; // always update, even if stationary; required to get correct shadows when player stands still and takes/moves objects
					shadow_caster_hash ^= 0xdeadbeef; // update hash when player enters or leaves the light's area
				}
			}
		}
		if (!force_smap_update) {
			bool check_dynamic_shadows(camera_near_building);
			// handle people visible through skylights when player is above and light is on the top floor of a room with a skylight
			check_dynamic_shadows |= (room.get_has_skylight() && camera_bs.z > lpos.z && lpos.z > (room.z2() - 0.5*window_vspacing));
			
			if (check_dynamic_shadows) {
				// use full light radius for attics since they're more open
				float const dshadow_radius((is_in_attic ? 1.0 : (in_retail_room ? RETAIL_SMAP_DSCALE : PERSON_INT_SMAP_DSCALE))*light_radius);
				if (building_action_key) {force_smap_update = 1;} // toggling a door state or interacting with objects invalidates shadows in the building for that frame
				check_for_dynamic_shadow_casters(interior->people, ped_bcubes, moving_objs, clipped_bc, lpos_rot,
					dshadow_radius, stairs_light, xlate, (check_building_people && !is_lamp), shadow_caster_hash); // no people shadows for lam[s
			}
		}
		if (!force_smap_update && interior->last_active_door_ix >= 0) { // check for door opening or closing; since this is player controlled, there should be at most one
			door_t const &door(interior->get_door(interior->last_active_door_ix));
			force_smap_update |= (door.is_partially_open() && clipped_bc.intersects(door));
		}
		// end dynamic shadows check
		cube_t const clipped_bc_rot(is_rotated() ? get_rotated_bcube(clipped_bc) : clipped_bc);
		setup_light_for_building_interior(dl_sources.back(), *i, clipped_bc_rot, force_smap_update, shadow_caster_hash);
		
		// add upward pointing light (sideways for wall lights); only when player is near/inside a building (optimization); not for lights hanging on ceiling fans
		if ((camera_near_building || in_walkway_near_camera) && (is_lamp || wall_light || lpos_rot.z > up_light_zmin) && !i->is_hanging()) {
			cube_t light_bc2(clipped_bc);

			if (is_in_elevator) {
				light_bc2.intersect_with_cube(get_elevator(i->obj_id)); // clip to elevator to avoid light leaking onto walls outside but near the elevator
			}
			else if (!is_in_attic && !is_exterior) {
				// expand slightly so that points exactly on the room bounds and exterior doors are included; not for backrooms because it already contains the wall width
				cube_t room_exp(get_walkable_room_bounds(room));
				room_exp.expand_by((room.is_backrooms() ? 0.1 : 1.0)*room_xy_expand); // smaller expand for backrooms

				if (room.open_wall_mask && !room.is_hallway) { // don't clamp on open wall sides, except for hallways
					for (unsigned d = 0; d < 2; ++d) {
						if (!room.has_open_wall(d, 0)) {max_eq(light_bc2.d[d][0], room_exp.d[d][0]);} // clamp lo
						if (!room.has_open_wall(d, 1)) {min_eq(light_bc2.d[d][1], room_exp.d[d][1]);} // clamp hi
					}
					max_eq(light_bc2.z1(), room_exp.z1());
					min_eq(light_bc2.z2(), room_exp.z2()); // not needed?
				}
				else {
					light_bc2.intersect_with_cube(room_exp); // upward facing light is for this room only
				}
				min_eq(light_bc2.z2(), (ceil_z + fc_thick)); // doesn't reach higher than the ceiling of this room
			}
			if (is_rotated()) {light_bc2 = get_rotated_bcube(light_bc2);}
			float sec_light_radius(0.0);
			point sec_lpos(lpos_rot);

			if (is_lamp) { // add a second shadowed light source pointing up
				bool const cache_shadows(!lamp_was_moved);
				dl_sources.emplace_back(light_radius, lpos_rot, lpos_rot, color, 0, plus_z, 0.5*bwidth); // points up
				// lamps are static and have no dynamic shadows, so always cache their shadow maps
				assign_light_for_building_interior(dl_sources.back(), &(*i), light_bc2, cache_shadows, 1); // is_lamp=1
				sec_light_radius = 0.15*light_radius;
				dl_sources.emplace_back(sec_light_radius, sec_lpos, sec_lpos, color); // add an additional small unshadowed light for ambient effect
			}
			else { // add a second, smaller unshadowed light for the upper hemisphere or omnidirectional for wall lights
				// the secondary light is unshadowed and won't pick up shadows from any stairs in the room, so reduce the radius
				float const rscale((room.is_hallway ? 0.25 : (room.is_office ? 0.45 : 0.5))*(has_stairs_this_floor ? 0.67 : 1.0));
				sec_light_radius = rscale*light_radius;

				if (wall_light) {
					sec_lpos += (2.0*i->get_sz_dim(i->dim))*dir; // shift to the side
					dl_sources.emplace_back(sec_light_radius, sec_lpos, sec_lpos, color, 0, -dir, 1.0); // omnidirectional
					// clip to the centerline of the wall itself; needed for stairs wall lights
					light_bc2.d[i->dim][!i->dir] = i->d[i->dim][!i->dir] + (i->dir ? -1.0 : 1.0)*0.5*wall_thickness;
				}
				else { // ceiling light
					sec_lpos.z -= wall_thickness; // shift down somewhat; use a constant that works with recessed lights at doorways
					dl_sources.emplace_back(sec_light_radius, sec_lpos, sec_lpos, color, 0, plus_z, 0.5); // hemisphere that points up
				}
			}
			if (!light_bc2.is_all_zeros()) {
				cube_t sphere_bc;
				sphere_bc.set_from_sphere(sec_lpos, sec_light_radius);
				assert(light_bc2.intersects(sphere_bc));
				light_bc2.intersect_with_cube(sphere_bc);
				dl_sources.back().set_custom_bcube(light_bc2);
			}
			dl_sources.back().disable_shadows();
		} // end upward light
	} // for i (objs)
	//if (camera_in_building) {cout << name << ": " << num_add << endl;} // TESTING

	// add skylight lights
	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
		if (skylights.empty() || camera_in_basement) continue;
		point sun_moon_pos;
		if (!light_valid_and_enabled(l) || !get_light_pos(sun_moon_pos, l)) continue;
		float const camera_zval_check((camera_somewhat_by_stairs ? 2.0 : 1.0)*window_vspacing); // allow two floors below if player can see skylight from stairs
		bool const dir_always_vert  = 1;
		bool const add_sky_lighting = 1;

		for (cube_with_ix_t &sl : skylights) { // add virtual spotlights for each skylight to simulate sun light
			if (!lights_bcube.intersects_xy(sl))        continue; // not contained within the light volume
			if (camera_z < sl.z1() - camera_zval_check) continue; // player below the floor with the skylight or the one below; invalid when flying outside the building?
			cube_t lit_area(sl);
			lit_area.z2() += fc_thick; // include the tops of the skylight
			point lpos(sl.get_cube_center());
			if (is_rotated()) {do_xy_rotate(building_center, lpos);} // ???
			vector3d const light_dir((sun_moon_pos - lpos).get_norm());
			float const light_dist(3.0*window_vspacing); // larger is more physically correct (directional), but produces lower shadow resolution due to wasted texels
			lpos += light_dist*light_dir;
			bcube.clamp_pt_xy(lpos); // must be within the XY bounds of the bcube to pick up shadows from this building
			bool room_has_stairs(0);

			for (room_t const &room : interior->rooms) {
				if (room.get_has_skylight() && room.intersects(sl)) {
					lit_area.union_with_cube_xy(room);
					unsigned const num_floors(calc_num_floors_room(room, window_vspacing, 2.0*fc_thick));
					room_has_stairs |= room.has_stairs_on_floor(num_floors - 1); // check for stairs on top floor
				}
			}
			if (camera_somewhat_by_stairs && !room_has_stairs && camera_z < sl.z1() - window_vspacing) continue; // player below the floor with the skylight
			lit_area.z1() -= (room_has_stairs ? 2.0 : 1.0)*window_vspacing; // floor below the skylight; two floors if a room has stairs
			lit_area.expand_by_xy(room_xy_expand); // include walls
			if (!is_rot_cube_visible(lit_area, xlate, 1)) continue; // VFC; inc_mirror_reflections=1
			// further constrain bcube based on projected light rays through the corners of the skylight to prevent unshadowed light from passing through exterior walls
			cube_t proj_bcube(sl);
			float const z(sl.z2()); // top edge
			point const corners[4] = {point(sl.x1(), sl.y1(), z), point(sl.x2(), sl.y1(), z), point(sl.x2(), sl.y2(), z), point(sl.x1(), sl.y2(), z)};
			float const extend_amt(window_vspacing/light_dist); // length of light ray reaching the floor is this much longer than the length of the ray to the skylight
			for (unsigned n = 0; n < 4; ++n) {proj_bcube.union_with_pt(corners[n] + extend_amt*(corners[n] - lpos));}
			cube_t clipped_area(lit_area);
			clipped_area.intersect_with_cube_xy(proj_bcube);
			colorRGBA const outdoor_color(get_outdoor_light_color()); // required to calculate cur_ambient and cur_diffuse
			float const dx(sl.dx()), dy(sl.dy()), diag_sz(sqrt(dx*dx + dy*dy));
			float const light_radius(1.2*diag_sz + 1.5*light_dist); // determined experimentally
			unsigned const dlights_start(dl_sources.size());

			if (!is_rot_cube_visible(clipped_area, xlate, 1)) {} // VFC, again; inc_mirror_reflections=1
			else if ((display_mode & 0x08) && !clipped_area.contains_pt(camera_rot) && check_obj_occluded(clipped_area, camera_bs, oc, 0)) {} // occlusion culling
			else { // add the light
				colorRGBA const color(add_sky_lighting ? colorRGBA(cur_diffuse) : outdoor_color);
				float corner_horiz_dist(0.0);

				if (dir_always_vert) { // make the light dir vertical/Z to avoid aliasing artifacts, though this is less physically correct
					for (unsigned n = 0; n < 4; ++n) {max_eq(corner_horiz_dist, p2p_dist_xy(lpos, corners[n]));} // find the furthest corner
				}
				else { // change the light direction correctly
					corner_horiz_dist = 0.5*diag_sz; // calculate the radius
				}
				float const dp(light_dist/sqrt(light_dist*light_dist + corner_horiz_dist*corner_horiz_dist)), bwidth(0.5*(1.0 - dp));
				// check for dynamic shadow casters
				bool force_smap_update(building_action_key); // update if a door is open or closed
				unsigned shadow_caster_hash(0);

				if (camera_surf_collide && camera_in_building && clipped_area.contains_pt(camera_rot)) {
					force_smap_update   = 1;
					shadow_caster_hash ^= 0xdeadbeef; // update hash when player enters or leaves the light's area
				}
				if (!force_smap_update && ((camera_near_building && camera_bs.z > clipped_area.z1()) || camera_bs.z > z)) {
					check_for_dynamic_shadow_casters(interior->people, ped_bcubes, moving_objs, clipped_area, lpos,
						0.0, room_has_stairs, xlate, check_building_people, shadow_caster_hash); // dmax=0
				}
				hash_mix_point(lpos, shadow_caster_hash); // update when light (sun/moon) pos changes
				bool const cache_shadows(!force_smap_update && sl.ix == shadow_caster_hash);
				sl.ix = shadow_caster_hash; // store new hashval in the skylight for next frame
				vector3d const spotlight_dir(dir_always_vert ? -plus_z : -light_dir); // points either downward or away from the sun/moon
				dl_sources.emplace_back(light_radius, lpos, lpos, color, 0, spotlight_dir, bwidth);
				assign_light_for_building_interior(dl_sources.back(), &sl, clipped_area, cache_shadows);
			}
			if (add_sky_lighting) { // add a weaker unshadowed vertical light using cur_ambient
				point lpos2(sl.get_cube_center());
				lpos2.z += light_dist; // directly above the skylight
				cube_t proj_bcube2(sl);
				for (unsigned n = 0; n < 4; ++n) {proj_bcube2.union_with_pt(corners[n] + extend_amt*(corners[n] - lpos2));}
				cube_t clipped_area2(lit_area);
				clipped_area2.intersect_with_cube_xy(proj_bcube2);

				if (!is_rot_cube_visible(clipped_area2, xlate, 1)) {} // VFC; inc_mirror_reflections=1
				else if ((display_mode & 0x08) && !clipped_area2.contains_pt(camera_rot) && check_obj_occluded(clipped_area2, camera_bs, oc, 0)) {} // occlusion culling
				else { // add the light
					dl_sources.emplace_back(light_radius, lpos2, lpos2, cur_ambient, 0);
					dl_sources.back().set_custom_bcube(clipped_area2);
					dl_sources.back().disable_shadows();
				}
			}
			if (dl_sources.size() > dlights_start) { // a light was added, update lights_bcube
				expand_cube_zvals(lights_bcube, lit_area.z1(), lpos.z); // must include the light
			}
		} // for skylight sl
	} // for l
	if (camera_in_building) {
		interior->room_geom->particle_manager.add_lights(xlate, *this, oc, lights_bcube);
		interior->room_geom->fire_manager    .add_lights(xlate, *this, oc, lights_bcube);
	}
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

bool check_bcube_visible_for_building(cube_t const &bcube, vector3d const &xlate, building_t const &building, occlusion_checker_noncity_t &oc, cube_t &lights_bcube) {
	if (!lights_bcube.intersects_xy(bcube)) return 0;
	cube_t const bcube_cs(bcube + xlate); // maybe we should check each parent object separately? but it could be rare if there's more than one
	if (!camera_pdu.cube_visible(bcube_cs)) return 0; // no particles are visible
	if ((display_mode & 0x08) && building.check_obj_occluded(bcube_cs, get_camera_pos(), oc)) return 0;
	return 1;
}
void particle_manager_t::add_lights(vector3d const &xlate, building_t const &building, occlusion_checker_noncity_t &oc, cube_t &lights_bcube) const {
	if (particles.empty()) return;
	if (!check_bcube_visible_for_building(get_bcube(), xlate, building, oc, lights_bcube)) return;

	for (particle_t const &p : particles) { // only sparks create light
		if (p.effect == PART_EFFECT_SPARK) {add_dlight_if_visible(p.pos, 20.0*p.radius, p.color, xlate, lights_bcube);}
	}
}
void fire_manager_t::add_lights(vector3d const &xlate, building_t const &building, occlusion_checker_noncity_t &oc, cube_t &lights_bcube) const {
	if (fires.empty()) return;
	if (!check_bcube_visible_for_building(get_bcube(), xlate, building, oc, lights_bcube)) return;
	for (fire_t const &f : fires) {add_dlight_if_visible(f.get_center(), 5.0*f.radius, ORANGE, xlate, lights_bcube);}
}

float room_t::get_light_amt() const { // Note: not normalized to 1.0
	float ext_perim(0.0);

	// add length of each exterior side, assuming it has windows; this is approximate because it treats partially exterior walls and fully exterior
	for (unsigned d = 0; d < 4; ++d) {
		if (ext_sides & (1<<d)) {ext_perim += get_sz_dim(d>>1);}
	}
	return ext_perim/get_area_xy(); // light per square meter = exterior perimeter over area
}

