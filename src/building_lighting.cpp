// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "lightmap.h" // for light_source
#include "cobj_bsp_tree.h"
#include "profiler.h"
#include <thread>

bool const USE_BKG_THREAD      = 1;
bool const INDIR_BASEMENT_ONLY = 0;
bool const INDIR_VOL_PER_FLOOR = 0;

extern bool camera_in_building, player_in_attic;
extern int MESH_Z_SIZE, display_mode, display_framerate, camera_surf_collide, animate2, frame_counter, building_action_key, player_in_basement, player_in_elevator;
extern unsigned LOCAL_RAYS, MAX_RAY_BOUNCES, NUM_THREADS;
extern float indir_light_exp;
extern double camera_zh, tfticks;
extern colorRGB cur_ambient, cur_diffuse;
extern std::string lighting_update_text;
extern vector<light_source> dl_sources;

bool enable_building_people_ai();
bool check_cube_occluded(cube_t const &cube, vect_cube_t const &occluders, point const &viewer);
void calc_cur_ambient_diffuse();
colorRGBA get_textured_wood_color();
bool bed_has_canopy_mat(room_object_t const &c);
int get_canopy_texture();
colorRGBA get_canopy_base_color(room_object_t const &c);
void get_water_heater_cubes(room_object_t const &wh, cube_t cubes[2]);

bool enable_building_indir_lighting_no_cib() {
	if (!(display_mode & 0x10)) return 0; // key 5
	if (INDIR_BASEMENT_ONLY && !player_in_basement) return 0;
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
// Note: static objects only; excludes people; pos in building space
bool building_t::ray_cast_interior(point const &pos, vector3d const &dir, cube_t const &valid_area, cube_bvh_t const &bvh,
	point &cpos, vector3d &cnorm, colorRGBA &ccolor, rand_gen_t *rgen) const
{
	if (!interior || is_rotated() || !is_simple_cube()) return 0; // these cases are not yet supported
	float const extent(valid_area.get_max_extent());
	cube_t clip_cube(valid_area);
	clip_cube.expand_by(0.01*extent); // expand slightly so that collisions with objects on the edge are still considered interior
	point p1(pos), p2(pos + dir*(2.0*extent));
	if (!do_line_clip(p1, p2, clip_cube.d)) return 0; // ray does not intersect clip cube
	building_mat_t const &mat(get_material());
	float t(1.0); // start at p2
	bool hit(0);
	cpos = p2; // use far clip point for clip cube if there is no hit

	// check parts (exterior walls); should chimneys and porch roofs be included?
	if (ray_cast_exterior_walls(p1, p2, cnorm, t)) { // interior ray (common case) - find furthest exit point
		p2  = p1 + (p2 - p1)*t; t = 1.0; // clip p2 to t (minor optimization)
		hit = 1;
		if (p2.z < ground_floor_z1) {ccolor = get_basement_wall_texture().get_avg_color();} // basement wall
		else {ccolor = mat.wall_color.modulate_with(mat.wall_tex.get_avg_color());} // non-basement wall
	}
	else { // check for exterior rays (uncommon case)
		bool hit(0);
		for (auto p = parts.begin(); p != parts.end(); ++p) {hit |= ray_cast_cube(p1, p2, *p, cnorm, t);} // find closest entrance point
		
		if (hit) { // exterior hit (ray outside building) - don't need to check interior geometry
			cpos   = p1 + (p2 - p1)*t;
			ccolor = side_color.modulate_with(mat.side_tex.get_avg_color());
			return 1;
		}
	}
	bvh.ray_cast(p1, p2, cnorm, ccolor, t);

	if (t == 1.0) { // no intersection with bvh
		cpos = p2;
		if (!hit) return 0;
		if (rgen && p2.z > ground_floor_z1 && has_windows() && rgen->rand_bool()) return 0; // 50% chance of exiting through a window
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

template<typename T> void add_colored_cubes(vector<T> const &cubes, colorRGBA const &color, float z1, float z2, vect_colored_cube_t &cc) {
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (c->z1() < z2 && c->z2() > z1) {cc.emplace_back(*c, color);}
	}
}
void add_colored_cubes(cube_t const *const cubes, unsigned num_cubes, colorRGBA const &color, vect_colored_cube_t &cc) {
	for (unsigned n = 0; n < num_cubes; ++n) {
		if (!cubes[n].is_all_zeros()) {cc.emplace_back(cubes[n], color);}
	}
}
void building_t::gather_interior_cubes(vect_colored_cube_t &cc, int only_this_floor) const {
	if (!interior) return; // nothing to do
	building_mat_t const &mat(get_material());
	colorRGBA const wall_color(wall_color.modulate_with(mat.wall_tex.get_avg_color()));
	float z1(bcube.z1()), z2(bcube.z2()), stairs_z1(z1), stairs_z2(z2); // start with full bcube Z range

	if (only_this_floor >= 0) { // clip per light source to current floor
		float const floor_spacing(get_window_vspace());
		z1 += only_this_floor*floor_spacing;
		z2  = z1 + floor_spacing;
		stairs_z1 = z1 - floor_spacing; // stairs extend an extra floor up and down to block rays in stairwells
		stairs_z2 = z2 + floor_spacing;
	}
	for (unsigned d = 0; d < 2; ++d) {add_colored_cubes(interior->walls[d], wall_color, z1, z2, cc);}

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		if (e->z1() > z2 || e->z2() < z1) continue;
		cube_t cubes[5];
		unsigned const num_cubes(e->get_coll_cubes(cubes));
		// for now elevators are treated the same as walls with the same color, even though the inside of open elevators is wood
		add_colored_cubes(cubes, num_cubes, wall_color, cc); // can only assign the same color to all sides of the cube
	}
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
		if (i->open) continue; // add only closed doors
		if (i->z1() > z2 || i->z2() < z1) continue;
		cc.emplace_back(i->get_true_bcube(), WHITE);
	}
	colorRGBA const ccolor(is_house ? mat.house_ceil_color .modulate_with(mat.house_ceil_tex .get_avg_color()) : mat.ceil_color .modulate_with(mat.ceil_tex .get_avg_color()));
	colorRGBA const fcolor(is_house ? mat.house_floor_color.modulate_with(mat.house_floor_tex.get_avg_color()) : mat.floor_color.modulate_with(mat.floor_tex.get_avg_color()));
	colorRGBA const dcolor(detail_color.modulate_with(mat.roof_tex. get_avg_color()));
	add_colored_cubes(interior->ceilings, ccolor, z1, z2, cc);
	add_colored_cubes(interior->floors,   fcolor, z1, z2, cc);
	add_colored_cubes(details,            dcolor, z1, z2, cc); // should this be included?
	if (!has_room_geom()) return; // nothing else to add
	vect_room_object_t const &objs(interior->room_geom->objs);
	static vect_cube_t temp; // used across calls for subtracting holes
		
	for (auto c = objs.begin(); c != objs.end(); ++c) {
		if (!c->is_visible()) continue;
		if (c->type  == TYPE_ELEVATOR) continue; // elevator cars/internals can move so should not contribute to lighting
		if (c->type  == TYPE_SHOWER  ) continue; // transparent
		if (c->type  == TYPE_BLOCKER || c->type == TYPE_COLLIDER) continue; // blockers and colliders are not drawn
		// skip other object types that are too small, not cube shaped, or not interior
		if (c->type == TYPE_WALL_TRIM || c->type == TYPE_BOOK || c->type == TYPE_CRACK || c->type == TYPE_PLANT || c->type == TYPE_RAILING || c->type == TYPE_SHELVES ||
			c->type == TYPE_BOTTLE || c->type == TYPE_PEN || c->type == TYPE_PENCIL || c->type == TYPE_LG_BALL || c->type == TYPE_HANGER_ROD || c->type == TYPE_DRAIN ||
			c->type == TYPE_MONEY || c->type == TYPE_PHONE || c->type == TYPE_TPROLL || c->type == TYPE_SPRAYCAN || c->type == TYPE_MARKER || c->type == TYPE_BUTTON ||
			c->type == TYPE_SWITCH || c->type == TYPE_TAPE || c->type == TYPE_OUTLET || c->type == TYPE_PARK_SPACE || c->type == TYPE_RAMP || c->type == TYPE_PIPE ||
			c->type == TYPE_VENT || c->type == TYPE_BREAKER || c->type == TYPE_KEY || c->type == TYPE_HANGER || c->type == TYPE_FESCAPE || c->type == TYPE_CUP ||
			c->type == TYPE_CLOTHES || c->type == TYPE_LAMP || c->type == TYPE_OFF_CHAIR || c->type == TYPE_LIGHT || c->type == TYPE_SIGN || c->type == TYPE_PAPER ||
			c->type == TYPE_PLANT) continue;
		bool const is_stairs(c->type == TYPE_STAIR || c->type == TYPE_STAIR_WALL);
		if (c->z1() > (is_stairs ? stairs_z2 : z2) || c->z2() < (is_stairs ? stairs_z1 : z1)) continue;
		colorRGBA const color(c->get_color());
		
		if (c->shape == SHAPE_CYLIN || c->shape == SHAPE_SPHERE) {
			cube_t inner_cube(*c);

			if (c->type == TYPE_WHEATER) { // {tank, pipes}
				cube_t cubes[2];
				get_water_heater_cubes(*c, cubes);
				cc.emplace_back(cubes[1], color); // add pipes directly
				inner_cube = cubes[0]; // tank
			}
			float const shrink(c->get_radius()*(1.0 - 1.0/SQRT2));
			inner_cube.expand_by_xy(-shrink); // shrink to inscribed cube in XY
			if (c->shape == SHAPE_SPHERE) {inner_cube.expand_in_dim(2, -shrink);} // shrink in Z as well
			if (c->type  == TYPE_TABLE  ) {inner_cube.z1() += 0.88*c->dz();} // top of table
			cc.emplace_back(inner_cube, color);
		}
		else if (c->type == TYPE_CLOSET) {
			cube_t cubes[5];
			get_closet_cubes(*c, cubes, 1); // for_collision=1
			add_colored_cubes(cubes, 5, color, cc); // include door, whether closed or open
		}
		else if (c->type == TYPE_BED) { // Note: posts are not included
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
			get_tc_leg_cubes(cubes[5], 0.04, cubes); // head_width=0.04; cubes[5] is not overwritten
			add_colored_cubes(cubes,   4, wood_color,   cc); // legs
		}
		else if (c->type == TYPE_DESK || c->type == TYPE_DRESSER || c->type == TYPE_NIGHTSTAND || c->type == TYPE_TABLE) { // objects with legs
			if (c->is_glass_table()) continue; // skip glass table (transparent with thin legs)
			cube_t cubes[5];
			get_table_cubes(*c, cubes); // top and 4 legs
			add_colored_cubes(cubes, 5, color, cc);
			if      (c->type == TYPE_DRESSER || c->type == TYPE_NIGHTSTAND) {cc.emplace_back(get_dresser_middle(*c), color);}
			else if (c->type == TYPE_DESK) {
				if (c->desk_has_drawers() ) {cc.emplace_back(get_desk_drawers_part(*c), color);}
				if (c->shape == SHAPE_TALL) {cc.emplace_back(get_desk_top_back(*c, cubes[0]), color);} // tall desk
			}
		}
		else if (c->type == TYPE_CHAIR) {
			colorRGBA const wood_color(get_textured_wood_color());
			cube_t cubes[3], leg_cubes[4]; // seat, back, legs_bcube
			get_chair_cubes(*c, cubes);
			cc.emplace_back(cubes[0], color     ); // seat
			cc.emplace_back(cubes[1], wood_color); // back
			get_tc_leg_cubes(cubes[2], 0.15, leg_cubes); // width=0.15
			add_colored_cubes(leg_cubes, 4, wood_color, cc);
		}
		else if (c->type == TYPE_CUBICLE || (c->type == TYPE_STALL && c->shape != SHAPE_SHORT)) { // cubicle or bathroom stall - hollow
			bool const is_stall(c->type != TYPE_CUBICLE);
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
		else { // single cube
			cube_t bc(*c); // handle 3D models that don't fill the entire cube

			if (c->type == TYPE_COUCH) {
				bc.z2() -= 0.5*bc.dz(); // bottom
				cube_t top(*c);
				top.z1() = bc.z2();
				top.d[c->dim][c->dir] += (c->dir ? -1.0 : 1.0)*0.6*bc.get_sz_dim(c->dim); // reduce thickness
				top.expand_in_dim(!c->dim, -0.1*c->get_sz_dim(!c->dim));
				cc.emplace_back(top, color);
			}
			else if (c->type == TYPE_TUB) {
				bc.z2() -= 0.9*bc.dz(); // bottom
				cube_t top(*c);
				top.z1() = bc.z2();
				cube_t inside(top);
				inside.expand_by_xy(-0.1*c->get_sz_dim(c->dim)); // shrink
				temp.clear();
				subtract_cube_from_cube(top, inside, temp);
				assert(temp.size() == 4); // -y, +y, -x, +x
				add_colored_cubes(temp, color, z1, z2, cc);
			}
			else if (c->type == TYPE_STOVE ) {bc.z2() -= 0.20*bc.dz();}
			else if (c->type == TYPE_TOILET) {bc.z2() -= 0.33*bc.dz();}
			else if (c->type == TYPE_SINK  ) {bc.z2() -= 0.20*bc.dz(); bc.z1() += 0.65*bc.dz();}
			else if (c->type == TYPE_MONITOR || c->type == TYPE_TV) {bc.expand_in_dim(c->dim, -0.3*bc.get_sz_dim(c->dim));} // reduce thickness
			else if (c->type == TYPE_BRSINK) {bc.z1() += 0.60*bc.dz();}
			else if (c->type == TYPE_ATTIC_DOOR) {bc = get_attic_access_door_cube(*c);} // Note: includes door but not ladder
			cc.emplace_back(bc, color);
		}
	} // for c
}


unsigned const IS_WINDOW_BIT = (1<<24); // if this bit is set, the light is from a window; if not, it's from a light room object

class building_indir_light_mgr_t {
	bool is_running, kill_thread, lighting_updated, needs_to_join, need_bvh_rebuild, update_windows, is_negative_light;
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
		// TODO: Some type of blur to remove noise that doesn't blur across walls
		// Note: modifies lmgr, but otherwise thread safe
		unsigned const num_rt_threads(max(1U, (NUM_THREADS - (USE_BKG_THREAD ? 1 : 0)))); // reserve a thread for the main thread if running in the background
		unsigned base_num_rays(LOCAL_RAYS), dim(2), dir(0); // default dim is z; dir=2 is omnidirectional
		cube_t const scene_bounds(get_scene_bounds_bcube()); // expected by lmap update code
		point const ray_scale(scene_bounds.get_size()/light_bounds.get_size()), llc_shift(scene_bounds.get_llc() - light_bounds.get_llc()*ray_scale);
		float const tolerance(1.0E-5*valid_area.get_max_extent());
		bool const is_window(cur_light & IS_WINDOW_BIT);
		float weight(100.0);
		cube_t light_cube;
		colorRGBA lcolor;
		assert(cur_light >= 0);

		if (is_window) { // window
			unsigned const window_ix(cur_light & ~IS_WINDOW_BIT);
			assert(window_ix < windows.size());
			cube_with_ix_t const &window(windows[window_ix]);
			assert(window.ix < 4); // encodes 2*dim + dir
			dim =  bool(window.ix >> 1);
			dir = !bool(window.ix &  1); // cast toward the interior
			// light intensity scales with surface area, since incoming light is a constant per unit area (large windows = more light)
			float const surface_area(window.dz()*window.get_sz_dim(!bool(dim)));
			lcolor     = outdoor_color;
			weight    *= surface_area/0.0012f; // 1/4 the surface area weight of lights
			light_cube = window;
			light_cube.translate_dim(dim, (dir ? 1.0 : -1.0)*0.5*b.get_wall_thickness()); // shift slightly inside the building to avoid collision with the exterior wall
		}
		else { // room light or lamp, pointing downward
			vect_room_object_t const &objs(b.interior->room_geom->objs);
			assert((unsigned)cur_light < objs.size());
			room_object_t const &ro(objs[cur_light]);
			light_cube      = ro;
			light_cube.z1() = light_cube.z2() = (ro.z1() - 0.01*ro.dz()); // set slightly below bottom of light
			bool const light_in_basement(ro.z1() < b.ground_floor_z1), is_lamp(ro.type == TYPE_LAMP);
			if (is_lamp) {base_num_rays /= 2;} // half the rays for lamps
			if (is_lamp) {dir = 2;} // onmidirectional; dim stays at 2/Z
			float const surface_area(ro.dx()*ro.dy() + 2.0f*(ro.dx() + ro.dy())*ro.dz()); // bottom + 4 sides (top is occluded), 0.0003 for houses
			lcolor  = (is_lamp ? LAMP_COLOR : ro.get_color());
			weight *= surface_area/0.0003f;
			if (b.has_pri_hall())     {weight *= 0.70;} // floorplan is open and well lit, indir lighting value seems too high
			if (ro.type == TYPE_LAMP) {weight *= 0.33;} // lamps are less bright
			if (light_in_basement)    {weight *= (b.has_parking_garage ? 0.25 : 0.5);} // basement is darker, parking garages are even darker
		}
		if (b.is_house) {weight *= 2.0 ;} // houses have dimmer lights and seem to work better with more indir
		if (is_negative_light) {weight *= -1.0;}
		weight /= base_num_rays; // normalize to the number of rays
		unsigned NUM_PRI_SPLITS(is_window ? 4 : 16); // we're counting primary rays for windows, use fewer primary splits to reduce noise at the cost of increased time
		max_eq(base_num_rays, NUM_PRI_SPLITS);
		int const num_rays(base_num_rays/NUM_PRI_SPLITS);
		
		// Note: dynamic scheduling is faster, and using blocks doesn't help
#pragma omp parallel for schedule(dynamic) num_threads(num_rt_threads)
		for (int n = 0; n < num_rays; ++n) {
			if (kill_thread) continue;
			rand_gen_t rgen;
			rgen.set_state(n+1, cur_light); // should be deterministic, though add_path_to_lmcs() is not (due to thread races)
			vector3d pri_dir(rgen.signed_rand_vector_spherical(1.0).get_norm()); // should this be cosine weighted for windows?
			if (dir < 2 && ((pri_dir[dim] > 0.0) ^ dir)) {pri_dir[dim] *= -1.0;}
			point origin, init_cpos, cpos;
			vector3d init_cnorm, cnorm;
			colorRGBA ccolor(WHITE);

			// select a random point on the light cube (close enough for (ro.shape == SHAPE_CYLIN))
			for (unsigned d = 0; d < 3; ++d) {
				float const lo(light_cube.d[d][0]), hi(light_cube.d[d][1]);
				origin[d] = ((lo == hi) ? lo : rgen.rand_uniform(lo, hi));
			}
			init_cpos = origin; // init value
			bool const hit(b.ray_cast_interior(origin, pri_dir, valid_area, bvh, init_cpos, init_cnorm, ccolor, &rgen));

			// room lights lights already contribute direct lighting, so we skip this ray; however, windows don't, so we add their primary ray contribution
			if (is_window && init_cpos != origin) {
				point const p1(origin*ray_scale + llc_shift), p2(init_cpos*ray_scale + llc_shift); // transform building space to global scene space
				add_path_to_lmcs(&lmgr, nullptr, p1, p2, weight, lcolor*NUM_PRI_SPLITS, LIGHTING_LOCAL, 0); // local light, no bcube; scale color based on splits
			}
			if (!hit) continue; // done
			colorRGBA const init_color(lcolor.modulate_with(ccolor));
			if (init_color.get_weighted_luminance() < 0.1) continue; // done
			vector3d const v_ref(get_reflect_dir(pri_dir, init_cnorm));

			for (unsigned splits = 0; splits < NUM_PRI_SPLITS; ++splits) {
				point pos(origin);
				vector3d dir(pri_dir);
				colorRGBA cur_color(init_color);
				calc_reflect_ray(pos, init_cpos, dir, init_cnorm, v_ref, rgen, tolerance);

				for (unsigned bounce = 1; bounce < MAX_RAY_BOUNCES; ++bounce) { // allow up to MAX_RAY_BOUNCES bounces
					cpos = pos; // init value
					bool const hit(b.ray_cast_interior(pos, dir, valid_area, bvh, cpos, cnorm, ccolor, &rgen));

					if (cpos != pos) { // accumulate light along the ray from pos to cpos (which is always valid) with color cur_color
						point const p1(pos*ray_scale + llc_shift), p2(cpos*ray_scale + llc_shift); // transform building space to global scene space
						add_path_to_lmcs(&lmgr, nullptr, p1, p2, weight, cur_color, LIGHTING_LOCAL, 0); // local light, no bcube
					}
					if (!hit) break; // done
					cur_color = cur_color.modulate_with(ccolor);
					if (cur_color.get_weighted_luminance() < 0.1) break; // done
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
			dist_sq *= 0.05f*window_vspacing/(i->dz()*(i->dx() + i->dy())); // account for the size of the window, larger window smaller/higher priority
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
	building_indir_light_mgr_t() : is_running(0), kill_thread(0), lighting_updated(0), needs_to_join(0),
		need_bvh_rebuild(0), update_windows(0), is_negative_light(0), cur_bix(-1), cur_light(-1), cur_floor(-1), cur_tid(0) {}

	cube_t get_light_bounds() const {return light_bounds;}

	void invalidate_lighting() {
		is_negative_light = 0;
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
		bool floor_change(0);

		if ((int)bix != cur_bix) { // change to a different building
			clear();
			cur_bix = bix;
			assert(!is_running);
			build_bvh(b, target);
			get_windows(b);
		}
		else {
			if (update_windows) {get_windows(b);}
			floor_change = (cur_floor >= 0 && cur_floor != (int)b.get_floor_for_zval(target.z));
		}
		if (INDIR_VOL_PER_FLOOR && floor_change) { // handle floor change
			invalidate_lighting();
			build_bvh(b, target);
		}
		calc_cur_ambient_diffuse(); // needed for correct outdoor color
		colorRGBA const cur_outdoor_color(blend_color(cur_ambient, cur_diffuse, 0.5, 0)); // a mix of each

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
		if (!INDIR_VOL_PER_FLOOR) {need_bvh_rebuild |= floor_change;} // rebuild on player floor change
		if (need_bvh_rebuild) {build_bvh(b, target);}
		
		if (cur_light >= 0) {
			if (!is_negative_light) {lights_complete.insert(cur_light);} // mark the most recent light as complete if not a light removal
			cur_light = -1;
		}
		if (!remove_queue.empty()) { // remove an existing light; must run even when player_in_elevator==2 to remove elevator light at old pos
			cur_light = remove_queue.front();
			remove_queue.pop_front();
			is_negative_light = 1;
		}
		else { // find a new light to add
			if (player_in_elevator == 2) return; // pause updates for player in closed elevator since lighting is not visible
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
	void build_bvh(building_t const &b, point const &target) {
		//highres_timer_t timer("Build BVH");
		cur_floor  = b.get_floor_for_zval(target.z);
		valid_area = b.bcube;
		
		// clip per light source to current floor; note that this will exclude stairs going up or down
		if (b.point_in_attic(target)) {
			set_cube_zvals(valid_area, b.interior->attic_access.z1(), b.interior_z2);
		}
		else {
			float const floor_spacing(b.get_window_vspace());
			valid_area.z1() += cur_floor*floor_spacing;
			valid_area.z2()  = valid_area.z1() + floor_spacing;
		}
		light_bounds = (INDIR_VOL_PER_FLOOR ? valid_area : b.get_interior_bcube());
		bvh.clear();
		b.gather_interior_cubes(bvh.get_objs(), cur_floor);
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

bool building_t::ray_cast_camera_dir(point const &camera_bs, point &cpos, colorRGBA &ccolor) const {
	assert(!USE_BKG_THREAD); // not legal to call when running lighting in a background thread
	building_indir_light_mgr.build_bvh(*this, camera_bs);
	vector3d cnorm; // unused
	return ray_cast_interior(camera_bs, cview_dir, bcube, building_indir_light_mgr.get_bvh(), cpos, cnorm, ccolor);
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
		if (!i->is_light_type() || !i->is_light_on())  continue; // not a light, or light not on
		bool const light_in_basement(i->z1() < ground_floor_z1);
		if (INDIR_BASEMENT_ONLY && !light_in_basement) continue; // not a basement light
		if (!valid_area.contains_cube(*i))             continue; // outside valid area

		if (i->flags & RO_FLAG_IN_ELEV) { // elevator light
			elevator_t const &e(interior->elevators[i->obj_id]);
			if (e.was_called) continue; // moving elevator, don't update light yet
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
	tx1 = v0.t[0]; tx2 = tx1; tz1 = v0.t[1]; tz2 = tz1; // tex coord ranges (xy, z); should generally be whole integers

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

void building_t::get_all_windows(vect_cube_with_ix_t &windows) const { // Note: ix encodes 2*dim+dir
	windows.clear();
	if (!has_windows() || is_rotated()) return; // no windows; rotated buildings not handled
	float const border_mult(0.94); // account for the frame part of the window texture, which is included in the interior cutout of the window
	float const window_h_border(border_mult*get_window_h_border()), window_v_border(border_mult*get_window_v_border()); // (0, 1) range
	vect_room_object_t blinds;

	if (is_house && has_room_geom()) { // find all bedroom blinds and use them to partially block windows
		float const width_expand(get_wall_thickness());

		for (room_object_t const &c : interior->room_geom->objs) {
			if (c.type != TYPE_BLINDS) continue;
			blinds.push_back(c);
			blinds.back().expand_in_dim(c.dim, width_expand); // make sure the intersect the windows
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
					
					if (b.flags & RO_FLAG_HANGING) { // horizontal (comes from the top)
						min_eq(window.z2(), b.z1());
					}
					else { // vertical (comes from the sides)
						bool const side(b.get_center_dim(!dim) < window.get_center_dim(!dim));
						if (side) {max_eq(window.d[!dim][0], b.d[!dim][1]);} // left  side
						else      {min_eq(window.d[!dim][1], b.d[!dim][0]);} // right side
					}
				} // for b
				bool const is_blocked(window.dz() <= 0.0 || window.get_sz_dim(!dim) <= 0.0);
				windows.emplace_back(window, (4*is_blocked + 2*dim + dir)); // Note: must include blocked window for seen cache to work
			}
		} // for z
	} // for i
}

bool building_t::register_indir_lighting_state_change(unsigned light_ix, bool is_door_change) const {
	if (!enable_building_indir_lighting()) return 0; // no update needed
	assert(has_room_geom());
	assert(light_ix < interior->room_geom->objs.size());
	room_object_t const &obj(interior->room_geom->objs[light_ix]);
	assert(obj.is_light_type());
	if (is_door_change && !obj.is_light_on()) return 0; // light off, no state change
	building_indir_light_mgr.register_light_state_change(light_ix, obj.is_light_on(), (obj.flags & RO_FLAG_IN_ELEV), is_door_change);
	return 1;
}
void building_t::register_indir_lighting_geom_change() const {
	building_indir_light_mgr.invalidate_bvh();
}
void building_t::register_blinds_state_change() const {
	building_indir_light_mgr.invalidate_windows();
	register_indir_lighting_geom_change();
}

bool line_int_cubes(point const &p1, point const &p2, vect_cube_t const &cubes) {
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (c->line_intersects(p1, p2)) return 1;
	}
	return 0;
}
bool building_t::is_light_occluded(point const &lpos, point const &camera_bs) const {
	assert(interior);
	// Note: assumes the light is inside the building
	// exterior walls have windows and don't generally occlude lights; room objects and doors are too small to occlude; elevators are too sparse to occlude
	for (unsigned d = 0; d < 2; ++d) {
		if (line_int_cubes(lpos, camera_bs, interior->walls[d])) return 1;
	}
	if (line_int_cubes(lpos, camera_bs, interior->fc_occluders)) return 1;
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

void building_t::refine_light_bcube(point const &lpos, float light_radius, cube_t const &room, cube_t &light_bcube, bool is_parking_garage) const {
	// base: 173613 / bcube: 163942 / clipped bcube: 161455 / tight: 159005 / rays: 101205 / no ls bcube expand: 74538
	// starts with building bcube clipped to light bcube
	//highres_timer_t timer("refine_light_bcube"); // 0.035ms average
	assert(interior);
	cube_t tight_bcube, part;
	static vect_cube_t other_parts, walls[2];
	other_parts.clear();

	// first determine the union of all intersections with parts; ignore zvals here so that we get the same result for every floor
	if (is_parking_garage) {tight_bcube = part = room;} // parking garage is the entire room; other_parts remains empty
	else {
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
	cube_t rays_bcube(lpos, lpos), room_exp(room);
	float const wall_thickness(get_wall_thickness()), tolerance(0.01*wall_thickness);
	room_exp.expand_by_xy(wall_thickness + tolerance); // to include points on the border + some FP error
	// pre-compute the nearby walls we will use for clipping
	for (unsigned d = 0; d < 2; ++d) {walls[d].clear();}

	if (is_parking_garage) {
		auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
		unsigned const pg_wall_start(interior->room_geom->wall_ps_start);
		assert(pg_wall_start < interior->room_geom->objs.size());

		for (auto i = (interior->room_geom->objs.begin() + pg_wall_start); i != objs_end; ++i) {
			if (i->type != TYPE_PG_WALL || i->item_flags != 0) continue; // not parking garage wall (breaking is incorrect for multiple PG levels)
			if (tight_bcube.intersects(*i)) {walls[i->dim].push_back(*i);}
		}
	}
	else {
		for (unsigned d = 0; d < 2; ++d) {
			for (cube_t const &c : interior->walls[d]) {
				if (tight_bcube.intersects(c)) {walls[d].push_back(c);}
			}
		}
	}
	for (unsigned n = 0; n < NUM_RAYS; ++n) {
		float const angle(TWO_PI*n/NUM_RAYS), dx(light_radius*sin(angle)), dy(light_radius*cos(angle));
		point p2(lpos + point(dx, dy, 0.0));
		// test for bad rays; this can fail on rotated buildings if lights aren't rotated properly, and in other cases when I'm experimenting, so it's allowed
		if (!do_line_clip_xy_p2(lpos, p2, tight_bcube)) continue; // bad ray, skip

		if (other_parts.empty() || part.contains_pt_xy(p2)) {
			clip_ray_to_walls(lpos, p2, walls); // the simple case where we don't need to handle part boundaries (optimization)
		}
		else {
			// find the point where this ray exits the building by following it through all parts; parts should be exactly adjacent to each other in X or Y
			point cur_pt(p2);
			bool const ret(do_line_clip_xy_p2(lpos, cur_pt, part)); // exit point of the starting part
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
			else {clip_ray_to_walls(lpos, p2, walls);} // ray ends in another room, need to clip it to the building walls
		}
		rays_bcube.union_with_pt(p2);
	} // for n
	light_bcube = rays_bcube;
}

void assign_light_for_building_interior(light_source &ls, room_object_t const &obj, cube_t const &light_bcube, bool cache_shadows, bool is_lamp) {
	unsigned const smap_id(uintptr_t(&obj)/sizeof(void *) + is_lamp); // use memory address as a unique ID; add 1 for lmaps to keep them separate
	ls.assign_smap_mgr_id(1); // use a different smap manager than the city (cars + streetlights) so that they don't interfere with each other
	if (!light_bcube.is_all_zeros()) {ls.set_custom_bcube(light_bcube);}
	ls.assign_smap_id(cache_shadows ? smap_id : 0); // if cache_shadows, mark so that shadow map can be reused in later frames
	if (!cache_shadows) {ls.invalidate_cached_smap_id(smap_id);}
}
void setup_light_for_building_interior(light_source &ls, room_object_t &obj, cube_t const &light_bcube, bool force_smap_update, unsigned shadow_caster_hash) {
	// If there are no dynamic shadows, we can reuse the previous frame's shadow map;
	// hashing object positions should handle the case where a shadow caster moves out of the light's influence and leaves a shadow behind;
	// also need to handle the case where the light is added on the frame the room geom is generated when the shadow map is not yet created;
	// requiring two consecutive frames of no dynamic objects should fix this
	uint16_t const sc_hash16(shadow_caster_hash ^ (shadow_caster_hash >> 16)); // combine upper and lower 16 bits into a single 16-bit value
	// cache if no objects moved (based on position hashing) this frame or last frame, and we're not forced to do an update
	bool const shadow_update(obj.item_flags != sc_hash16), cache_shadows(!shadow_update && !force_smap_update && (obj.flags & RO_FLAG_NODYNAM));
	assign_light_for_building_interior(ls, obj, light_bcube, cache_shadows, 0); // is_lamp=0
	if (shadow_update) {obj.flags &= ~RO_FLAG_NODYNAM;} else {obj.flags |= RO_FLAG_NODYNAM;} // store prev update state in object flag
	obj.item_flags = sc_hash16; // store current object hash in item flags
}

cube_t building_t::get_rotated_bcube(cube_t const &c) const {
	if (!is_rotated()) return c;
	point const center(bcube.get_cube_center());
	float const z(c.z2()); // top edge
	point corners[4] = {point(c.x1(), c.y1(), z), point(c.x2(), c.y1(), z), point(c.x2(), c.y2(), z), point(c.x1(), c.y2(), z)};
	for (unsigned n = 0; n < 4; ++n) {do_xy_rotate(center, corners[n]);}
	cube_t ret;
	ret.set_from_points(corners, 4);
	ret.z1() = c.z1();
	return ret;
}
bool building_t::is_rot_cube_visible(cube_t const &c, vector3d const &xlate) const {
	return camera_pdu.cube_visible((is_rotated() ? get_rotated_bcube(c) : c) + xlate);
}

float get_radius_for_room_light(room_object_t const &obj) {return 6.0f*(obj.dx() + obj.dy());}

bool check_for_shadow_caster(vect_cube_t const &cubes, cube_t const &light_bcube, point const &lpos,
	float dmax, bool has_stairs, vector3d const &xlate, unsigned &shadow_caster_hash)
{
	bool ret(0);

	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (lpos.z < c->z2()) continue; // light is below the object; assumes lights are spotlights pointed downward
		if (!c->intersects(light_bcube)) continue; // object not within light area of effect
		point const center(c->get_cube_center());
		if (!dist_less_than(lpos, center, dmax)) continue; // too far from light to cast a visible shadow
		
		if (!has_stairs) { // the below check is inaccurate, skip and just return 1
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
			//if ((display_mode & 0x08) && check_obj_occluded((cube_ext + xlate), get_camera_pos(), oc, 0)) continue; // occlusion culling - is this useful?
			if (!camera_pdu.cube_visible(cube_ext + xlate)) { // VFC
				shadow_caster_hash += hash_point(c->get_size()); // hash size rather than position to include the object's presence but not its position (sort of)
				continue;
			}
		}
		shadow_caster_hash += hash_point(center);
		ret = 1;
	} // for c
	return ret;
}

// Note: non const because this caches light_bcubes
void building_t::add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building,
	occlusion_checker_noncity_t &oc, vect_cube_t &ped_bcubes, cube_t &lights_bcube)
{
	if (!has_room_geom()) return; // error?
	point const camera_bs(camera_pdu.pos - xlate), building_center(bcube.get_cube_center()); // camera in building space
	if ((display_mode & 0x08) && !bcube.contains_pt_xy(camera_bs) && is_entire_building_occluded(camera_bs, oc)) return; // check occlusion (optimization)
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), fc_thick(get_fc_thickness());
	float const camera_z(camera_bs.z), room_xy_expand(0.75*wall_thickness);
	bool const check_building_people(enable_building_people_ai()), check_attic(camera_in_building && has_attic() && interior->attic_access_open);
	bool const show_room_name(display_mode & 0x20); // debugging, key '6'
	cube_t const &attic_access(interior->attic_access);
	vect_cube_t &light_bcubes(interior->room_geom->light_bcubes);
	vect_room_object_t &objs(interior->room_geom->objs); // non-const, light flags are updated
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	point camera_rot(camera_bs);
	maybe_inv_rotate_point(camera_rot); // rotate camera pos into building space; should use this pos below except with building bcube, occlusion checks, or lpos_rot
	bool const player_on_attic_stairs(player_in_attic && interior->attic_access_open && interior->attic_access.contains_pt_xy(camera_rot));
	unsigned camera_part(parts.size()); // start at an invalid value
	unsigned camera_floor(0);
	bool camera_by_stairs(0), camera_on_stairs(0), camera_somewhat_by_stairs(0), camera_in_hallway(0), camera_near_building(camera_in_building);
	int camera_room(-1);
	vect_cube_t moving_objs;
	ped_bcubes.clear();

	if (camera_in_building) {
		run_light_motion_detect_logic(camera_bs);

		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (i->contains_pt(camera_rot)) {camera_part = (i - parts.begin()); break;}
		}
		int const room_ix(get_room_containing_pt(camera_rot));

		if (room_ix >= 0) {
			// Note: stairs connecting stacked parts aren't flagged with has_stairs because stairs only connect to the bottom floor, but they're partially handled below
			room_t const &room(get_room(room_ix));
			camera_floor      = max(0.0f, (camera_rot.z - room.z1()))/window_vspacing;
			camera_room       = room_ix;
			camera_by_stairs  = camera_somewhat_by_stairs = room.has_stairs_on_floor(camera_floor);
			camera_in_hallway = room.is_hallway;
			unsigned const room_type(room.get_room_type(camera_floor));
			assert(room_type < NUM_RTYPES);
			if (show_room_name) {lighting_update_text = room_names[room_type];}
			register_player_in_building(camera_bs, building_id); // required for AI following logic

			if (camera_by_stairs) { // by stairs - check if we're actually on the stairs
				for (auto i = interior->stairwells.begin(); i != interior->stairwells.end(); ++i) {
					if (i->contains_pt(camera_rot)) {camera_on_stairs = 1; break;}
				}
			}
			if (!camera_on_stairs && has_pg_ramp() && !interior->ignore_ramp_placement &&
				interior->pg_ramp.contains_pt(camera_rot - vector3d(0.0, 0.0, (CAMERA_RADIUS + camera_zh)))) // what about on a ramp?
			{
				camera_on_stairs = camera_by_stairs = camera_somewhat_by_stairs = 1; // ramp counts as stairs
			}
			if (!camera_somewhat_by_stairs) { // what about camera in room adjacent to one with stairs? maybe set camera_somewhat_by_stairs
				cube_t cr(get_room(camera_room));
				cr.expand_by_xy(2.0*wall_thickness);

				for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
					if (r->has_stairs_on_floor(camera_floor) && r->intersects_no_adj(cr)) {camera_somewhat_by_stairs = 1; break;}
				}
			}
		}
		else if (point_in_attic(camera_rot)) {
			if (show_room_name) {lighting_update_text = room_names[RTYPE_ATTIC];}
		}
		//lighting_update_text = ((is_sphere_lit(camera_rot, get_scaled_player_radius()) || is_sphere_lit((camera_rot - vector3d(0.0, 0.0, camera_zh)), get_scaled_player_radius())) ? "Lit" : "Unlit");
	}
	else {
		cube_t bcube_exp(bcube);
		bcube_exp.expand_by_xy(2.0*window_vspacing);
		camera_near_building = bcube_exp.contains_pt(camera_bs);
	}
	if (camera_near_building) { // build moving objects vector
		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->is_visible() && i->is_moving()) {moving_objs.push_back(*i);}
		}
	}
	if (bcube.contains_pt(camera_bs)) { // camera in building bcube; matches rat update logic
		for (rat_t const &rat : interior->room_geom->rats) {
			if (rat.is_moving()) {moving_objs.push_back(rat.get_bcube());}
		}
		for (spider_t const &spider : interior->room_geom->spiders) {
			if (spider.is_moving()) {moving_objs.push_back(spider.get_bcube());}
		}
	}
	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_on() || !i->is_light_type()) continue; // light not on, or not a light or lamp
		point lpos(i->get_cube_center()); // centered in the light fixture
		min_eq(lpos.z, (i->z2() - 0.0125f*window_vspacing)); // make sure the light isn't too close to the ceiling (if shifted up to avoid a door intersection)
		point lpos_rot(lpos);
		if (is_rotated()) {do_xy_rotate(building_center, lpos_rot);}
		if (!lights_bcube.contains_pt_xy(lpos_rot)) continue; // not contained within the light volume
		// basement lights are only visible if the player is inside the building on the basement or ground floor
		bool const light_in_basement(lpos.z < ground_floor_z1), is_in_elevator(i->flags & RO_FLAG_IN_ELEV);
		bool const is_in_closet(i->flags & RO_FLAG_IN_CLOSET), is_in_attic(i->flags & RO_FLAG_IN_ATTIC);
		if ((is_in_elevator || is_in_closet) && camera_z > lpos.z) continue; // elevator or closet light on the floor below the player
		if (light_in_basement && (camera_z > (ground_floor_z1 + window_vspacing) || !bcube.contains_pt(camera_bs))) continue;
		//if (is_light_occluded(lpos_rot, camera_bs))  continue; // too strong a test in general, but may be useful for selecting high importance lights
		//if (!camera_in_building && i->is_interior()) continue; // skip interior lights when camera is outside the building: makes little difference, not worth the trouble
		room_t const &room(get_room(i->room_id));
		bool const is_lamp(i->type == TYPE_LAMP), is_single_floor(room.is_sec_bldg || is_in_elevator);
		int const cur_floor(is_single_floor ? 0 : (i->z1() - room.z1())/window_vspacing); // garages and sheds are all one floor
		float const level_z(is_in_attic ? interior->attic_access.z1() : (room.z1() + cur_floor*window_vspacing));
		float const floor_z(level_z + fc_thick), ceil_z(is_in_attic ? interior_z2 : (is_single_floor ? room.z2() : (level_z + window_vspacing - fc_thick)));
		float const floor_below_zval(floor_z - window_vspacing), ceil_above_zval(ceil_z + window_vspacing);
		// Note: we use level_z rather than floor_z for floor_is_above test so that it agrees with the threshold logic for player_in_basement
		bool const floor_is_above((camera_z < level_z) && !is_single_floor), floor_is_below(camera_z > ceil_z);
		bool const camera_room_same_part(room.part_id == camera_part), has_stairs_this_floor(!is_in_attic && room.has_stairs_on_floor(cur_floor));
		bool const has_ramp(!interior->ignore_ramp_placement && is_room_above_ramp(room, i->z1()));
		bool const light_room_has_stairs_or_ramp(i->has_stairs() || has_stairs_this_floor || has_ramp);
		bool stairs_light(0), player_in_elevator(0);

		if (is_in_elevator) {
			assert(i->obj_id < interior->elevators.size());
			elevator_t const &e(interior->elevators[i->obj_id]);
			assert(e.car_obj_id < objs.size());
			room_object_t const &car(objs[e.car_obj_id]); // elevator car for this elevator
			if (car.contains_pt(camera_rot)) {player_in_elevator = 1;} // player inside elevator
			else if (e.open_amt == 0.0 && (floor_is_above || floor_is_below)) continue; // closed elevator viewed from a different floor
		}
		if (!player_in_elevator) { // none of the below culling applies when the player is in the elevator
			// if the light is in the basement and the camera isn't, it's not visible unless the player is by the stairs
			if ( light_in_basement && player_in_basement == 0 && !camera_somewhat_by_stairs) continue;
			
			if (!light_in_basement && player_in_basement == 2) { // the player is fully in the basement but the light isn't
				if (!light_room_has_stairs_or_ramp) continue; // it's not visible unless the room with the light has stairs or a ramp up to it (parking garages)

				if (!has_ramp) { // no ramp, but we know there are stairs - check for basement connector stairs
					bool has_basement_stairs(0);

					for (auto const &s : interior->stairwells) {
						if (s.z1() < ground_floor_z1 && s.intersects(room)) {has_basement_stairs = 1; break;}
					}
					if (!has_basement_stairs) continue; // room has stairs, but not basement stairs
				}
			}
			// less culling if either the light or the camera is by stairs and light is on the floor above or below
			if (camera_z > floor_below_zval && camera_z < ceil_above_zval) { // light is on the floor above or below the camera
				if (camera_in_hallway && camera_by_stairs && camera_room_same_part) {
					// special case for player in an office building primary hallway with stairs; only handle the case where the light is in the hallway above or below;
					// if camera is on the stairs or a ramp this also counts because this may be connecting two rooms in two different parts
					stairs_light = ((int)i->room_id == camera_room || camera_on_stairs);
				}
				else {
					stairs_light = (light_room_has_stairs_or_ramp || camera_somewhat_by_stairs); // either the light or the camera is by the stairs
				}
			}
			if (check_attic && floor_is_below && camera_bs.z > attic_access.z1() && room.contains_cube_xy(attic_access)) {
				// camera in attic, possibly looking down through attic access door, and light is in the room below - keep it
			}
			else if (check_attic && floor_is_above && lpos.z > attic_access.z2() && camera_room >= 0 && get_room(camera_room).contains_cube_xy(attic_access)) {
				// light in attic, and camera in room with attic access
			}
			else if (floor_is_above || floor_is_below) { // light is on a different floor from the camera
				// the basement is a different part, but it's still the same vertical stack; consider this the same effective part if the camera is in the basement above the room's part
				if (camera_in_building && (camera_room_same_part || ((player_in_basement || light_in_basement) && parts[room.part_id].contains_pt_xy(camera_rot)) ||
					(room.contains_pt_xy(camera_rot) && camera_z < ceil_above_zval && camera_z > floor_below_zval)))
				{
					// player is on a different floor of the same building part, or more than one floor away in a part stack, and can't see a light from the floor above/below
					if (!stairs_light) continue; // camera in building and on wrong floor, don't add light
				}
				else { // camera outside the building (or the part that contains this light)
					float const xy_dist(p2p_dist_xy(camera_bs, lpos_rot));
					if (!stairs_light && ((camera_z - lpos_rot.z) > 2.0f*xy_dist || (lpos_rot.z - camera_z) > 1.0f*xy_dist)) continue; // light viewed at too high an angle

					if (camera_in_building) { // camera and light are in different buildings/parts
						if (camera_part >= real_num_parts) continue; // camera in garage or shed
						assert(camera_part < parts.size());
						cube_t const &cpart(parts[camera_part]);
						if (cpart.z2() <= room.z1() || cpart.z1() >= room.z2()) continue; // light in a different vertical stack than the camera

						if (!is_rotated()) { // check exterior wall visibility; this part doesn't work for rotated buildings
							// is it better to check if light half sphere is occluded by the floor above/below?
							cube_t const &part(get_part_for_room(room));
							bool visible[2] = {0};

							for (unsigned d = 0; d < 2; ++d) { // for each dim
								bool const dir(camera_bs[d] > lpos_rot[d]);
								if ((camera_rot[d] > part.d[d][dir]) ^ dir) continue; // camera not on the outside face of the part containing this room, so can't see through any windows
								visible[d] = (room.ext_sides & (1 << (2*d + dir)));
							}
							if (!visible[0] && !visible[1]) continue; // room is not on the exterior of the building on either side facing the camera
						}
					}
				}
			} // end camera on different floor case
		} // end !player_in_elevator
		float const light_radius(get_radius_for_room_light(*i)), cull_radius(0.95*light_radius), dshadow_radius(0.8*light_radius); // what about light_radius for lamps?
		if (!camera_pdu.sphere_visible_test((lpos_rot + xlate), cull_radius)) continue; // VFC
		// check visibility of bcube of light sphere clipped to building bcube; this excludes lights behind the camera and improves shadow map assignment quality
		cube_t sphere_bc; // in building space
		sphere_bc.set_from_sphere(lpos, cull_radius);
		cube_t clipped_bc(sphere_bc);
		clipped_bc.intersect_with_cube(bcube);

		if (!stairs_light && !is_in_elevator) { // clip zval to current floor if light not in a room with stairs or elevator
			max_eq(clipped_bc.z1(), (floor_z - fc_thick));
		}
		min_eq(clipped_bc.z2(), (ceil_z + fc_thick)); // ceiling is always valid, since lights point downward
		assert(clipped_bc.is_strictly_normalized());
		if (!is_rot_cube_visible(clipped_bc, xlate)) continue; // VFC
		//if (line_intersect_walls(lpos, camera_rot)) continue; // straight line visibility test - for debugging, or maybe future use in assigning priorities
		//if (check_cube_occluded(clipped_bc, interior->fc_occluders, camera_rot)) continue; // legal, but may not help much
		
		// run flicker logic for broken lights; this is done later in the control flow because updating light gemetry can be expensive
		if (animate2 && i->is_broken()) {
			static rand_gen_t rgen;

			if (tfticks > i->light_amt) { // time for state transition
				float const delay_mult(i->is_open() ? 0.1 : 1.0);
				i->light_amt = tfticks + rgen.rand_uniform(0.1, 1.0)*TICKS_PER_SECOND*delay_mult; // schedule time for next transition
				i->flags    ^= RO_FLAG_OPEN;
				// regenerate lights geometry (can be somewhat slow); only update if player is below the level of the light
				if (camera_bs.z < i->z2()) {interior->room_geom->invalidate_lights_geom();}
			}
			if (!i->is_open()) continue; // not currently on
		}
		// update lights_bcube and add light(s)

		if (is_lamp) { // lamps are generally against a wall and not in a room with stairs and only illuminate that one room
			min_eq(lights_bcube.z1(), (lpos_rot.z - min(window_vspacing, light_radius)));
			max_eq(lights_bcube.z2(), (lpos_rot.z + min(window_vspacing, light_radius)));
		}
		else {
			min_eq(lights_bcube.z1(), (lpos_rot.z - light_radius));
			max_eq(lights_bcube.z2(), (lpos_rot.z + 0.1f*light_radius)); // pointed down - don't extend as far up
		}
		float const bwidth = 0.25; // as close to 180 degree FOV as we can get without shadow clipping
		colorRGBA color;
		unsigned shadow_caster_hash(0);

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
			elevator_t const &e(interior->elevators[i->obj_id]);
			room_object_t const &car(objs[e.car_obj_id]); // elevator car for this elevator
			assert(car.contains_pt(lpos));
			cube_t clip_cube(car); // light is constrained to the elevator car
			clip_cube.expand_in_dim(!e.dim, 0.1*room_xy_expand); // expand sides to include walls adjacent to elevator (enough to account for FP error)
			if (e.open_amt > 0.0) {clip_cube.d[e.dim][e.dir] += (e.dir ? 1.0 : -1.0)*light_radius;} // allow light to extend outside open elevator door
			clipped_bc.intersect_with_cube(clip_cube); // Note: clipped_bc is likely contained in clip_cube and could be replaced with it
			if (e.was_called) {shadow_caster_hash += hash_point(e.get_llc());} // make sure to update shadows if elevator is moving
			if (e.open_amt > 0.0 && e.open_amt < 1.0) {shadow_caster_hash += hash_by_bytes<float>()(e.open_amt);} // update shadows if door is opening or closing
		}
		else {
			if (is_in_attic) {} // nothing else to do?
			else if (room.is_sec_bldg) {clipped_bc.intersect_with_cube(room);} // secondary buildings only light their single room
			else {
				assert(i->obj_id < light_bcubes.size());
				cube_t &light_bcube(light_bcubes[i->obj_id]);

				if (is_lamp && (i->flags & RO_FLAG_MOVED)) {
					i->flags &= ~RO_FLAG_MOVED; // clear moved flag, since we noticed it was moved
					light_bcube.set_to_zeros(); // will be recalculated below
				}
				if (light_bcube.is_all_zeros()) { // not yet calculated - calculate and cache
					light_bcube = clipped_bc;
					refine_light_bcube(lpos, light_radius, room, light_bcube, (light_in_basement && has_parking_garage));
				}
				clipped_bc.x1() = light_bcube.x1(); clipped_bc.x2() = light_bcube.x2(); // copy X/Y but keep orig zvals
				clipped_bc.y1() = light_bcube.y1(); clipped_bc.y2() = light_bcube.y2();
			}
			clipped_bc.expand_by_xy(room_xy_expand); // expand so that offset exterior doors are properly handled
			clipped_bc.intersect_with_cube(sphere_bc); // clip to original light sphere, which still applies (only need to expand at building exterior)
		}
		if (!clipped_bc.contains_pt(lpos)) {
			cout << "Error: Invalid light bcube: " << TXT(clipped_bc.str()) << TXT(lpos.str()) << TXT(room.str()) << TXT(bcube.str()) << TXT(is_lamp) << endl;
			continue; // can fail in rare cases when very far from the origin, likely due to FP error, so skip light in this case
		}
		if (!is_rot_cube_visible(clipped_bc, xlate)) continue; // VFC - post clip
		if ((display_mode & 0x08) && !clipped_bc.contains_pt(camera_rot) && check_obj_occluded(clipped_bc, camera_bs, oc, 0)) continue; // occlusion culling
		dl_sources.emplace_back(light_radius, lpos_rot, lpos_rot, color, 0, -plus_z, bwidth); // points down
		bool force_smap_update(0);

		// check for dynamic shadows
		if (camera_surf_collide && camera_in_building && dist_less_than(lpos_rot, camera_bs, dshadow_radius)) {
			bool player_on_ladder_this_room(player_on_attic_stairs && (is_in_attic || room.intersects_xy(interior->attic_access)));

			if (clipped_bc.contains_pt(camera_rot) || player_on_ladder_this_room) {
				// must update shadow maps for the room above if the player is on the stairs or in the same room when there are stairs
				bool const check_floor_above(camera_on_stairs || (camera_by_stairs && camera_room == i->room_id));

				if (is_lamp || player_on_ladder_this_room || (player_in_attic && is_in_attic) ||
					(lpos_rot.z > camera_bs.z && (check_floor_above || lpos_rot.z < (camera_bs.z + window_vspacing))))
				{
					// player shadow; includes lamps (with no zval test)
					force_smap_update   = 1; // always update, even if stationary; required to get correct shadows when player stands still and takes/moves objects
					shadow_caster_hash ^= 0xdeadbeef; // update hash when player enters or leaves the light's area
				}
			}
		}
		if (!force_smap_update && camera_near_building) {
			if (building_action_key) {
				force_smap_update = 1; // toggling a door state or interacting with objects will generally invalidate shadows in the building for that frame
			}
			if (check_building_people && !is_lamp) { // update shadow_caster_hash for moving people, but not for lamps, because their light points toward the floor
				if (ped_bcubes.empty()) { // get all cubes on first light
					for (person_t const &p : interior->people) {ped_bcubes.push_back(p.get_bcube());}
				}
				check_for_shadow_caster(ped_bcubes, clipped_bc, lpos_rot, dshadow_radius, stairs_light, xlate, shadow_caster_hash);
			}
			// update shadow_caster_hash for moving objects
			check_for_shadow_caster(moving_objs, clipped_bc, lpos_rot, dshadow_radius, stairs_light, xlate, shadow_caster_hash);
		}
		// end dynamic shadows check
		cube_t const clipped_bc_rot(is_rotated() ? get_rotated_bcube(clipped_bc) : clipped_bc);
		setup_light_for_building_interior(dl_sources.back(), *i, clipped_bc_rot, force_smap_update, shadow_caster_hash);
		
		if (camera_near_building && (is_lamp || lpos_rot.z > camera_bs.z)) { // only when the player is near/inside a building (optimization)
			cube_t light_bc2(clipped_bc);

			if (is_in_elevator) {
				light_bc2.intersect_with_cube(interior->elevators[i->obj_id]); // clip to elevator to avoid light leaking onto walls outside but near the elevator
			}
			else if (!is_in_attic) {
				cube_t room_exp(room);
				room_exp.expand_by(room_xy_expand); // expand slightly so that points exactly on the room bounds and exterior doors are included
				light_bc2.intersect_with_cube(room_exp); // upward facing light is for this room only
				min_eq(light_bc2.z2(), (ceil_z + fc_thick)); // doesn't reach higher than the ceiling of this room
			}
			if (is_rotated()) {light_bc2 = get_rotated_bcube(light_bc2);}

			if (is_lamp) { // add a second shadowed light source pointing up
				dl_sources.emplace_back(light_radius, lpos_rot, lpos_rot, color, 0, plus_z, 0.5*bwidth); // points up
				// lamps are static and have no dynamic shadows, so always cache their shadow maps
				assign_light_for_building_interior(dl_sources.back(), *i, light_bc2, 1, 1); // cache_shadows=1, is_lamp=1
				dl_sources.emplace_back(0.15*light_radius, lpos_rot, lpos_rot, color); // add an additional small unshadowed light for ambient effect
				dl_sources.back().set_custom_bcube(light_bc2); // not sure if this is helpful, but should be okay
			}
			else { // add a second, smaller unshadowed light for the upper hemisphere
				point const lpos_up(lpos_rot - vector3d(0.0, 0.0, 2.0*i->dz()));
				// the upward pointing light is unshadowed and won't pick up shadows from any stairs in the room, so reduce the radius
				float const rscale((room.is_hallway ? 0.25 : room.is_office ? 0.45 : 0.5)*(has_stairs_this_floor ? 0.67 : 1.0));
				dl_sources.emplace_back(rscale*light_radius, lpos_up, lpos_up, color, 0, plus_z, 0.5);
			}
			if (!light_bc2.is_all_zeros()) {dl_sources.back().set_custom_bcube(light_bc2);}
			dl_sources.back().disable_shadows();
		}
	} // for i
}

float room_t::get_light_amt() const { // Note: not normalized to 1.0
	float ext_perim(0.0);

	// add length of each exterior side, assuming it has windows; this is approximate because it treats partially exterior walls and fully exterior
	for (unsigned d = 0; d < 4; ++d) {
		if (ext_sides & (1<<d)) {ext_perim += get_sz_dim(d>>1);}
	}
	return ext_perim/get_area_xy(); // light per square meter = exterior perimeter over area
}

