// 3D World - Scenery Classes Implementations (plants, rocks, logs stumps)
// by Frank Gennari
// 12/5/02

#include "scenery.h"
#include "mesh.h"
#include "shaders.h"
#include "gl_ext_arb.h"


bool     const USE_VOXEL_ROCKS = 0;
unsigned const ROCK_NDIV       = 24;
unsigned const ROCK_VOX_SZ     = 32;
float    const SHADOW_VAL      = 0.5;
float    const PT_LINE_THRESH  = 800.0;


colorRGBA const stem_c(0.4, 0.6, 0.2, 1.0);
colorRGBA const leaf_c(0.7, 0.7, 0.7, 1.0);

// tid, stemc, leafc
plant_type const pltype[NUM_PLANT_TYPES] = {
	plant_type(MJ_LEAF_TEX, stem_c,   leaf_c, ALPHA0),
	plant_type(PLANT1_TEX,  stem_c,   leaf_c, ALPHA0),
	plant_type(PLANT2_TEX,  stem_c,   leaf_c, PURPLE),
	plant_type(PLANT3_TEX,  stem_c,   leaf_c, colorRGBA(0.9, 0.1, 0.05)),
	plant_type(PLANT4_TEX,  stem_c,   leaf_c, ALPHA0),
	plant_type(COFFEE_TEX,  LT_BROWN, WHITE,  ALPHA0)
};


int DISABLE_SCENERY(0), has_scenery(0), has_scenery2(0);


extern int num_trees, xoff2, yoff2, rand_gen_index, window_width, do_zoom, display_mode, draw_model, DISABLE_WATER;
extern float zmin, zmax_est, water_plane_z, tree_scale, vegetation, fticks, ocean_wave_height;
extern pt_line_drawer tree_scenery_pld; // we can use this for plant trunks
extern rand_gen_t global_rand_gen;


int get_bark_tex_for_tree_type(int type);


inline float get_pt_line_thresh() {return PT_LINE_THRESH*(do_zoom ? ZOOM_FACTOR : 1.0);}


bool skip_uw_draw(point const &pos, float radius) {

	// used in tiled terrain mode to skip underwater rocks - otherwise, the fog calculation is incorrect (needs special air/water fog transition handling)
	if (world_mode != WMODE_INF_TERRAIN || DISABLE_WATER || !(display_mode & 0x04)) return 0;
	return (get_camera_pos().z > water_plane_z && (pos.z + radius) < (get_water_z_height() - ocean_wave_height)); // water_plane_z
}


// ************ SCENERY OBJECT CLASSES ************


bool scenery_obj::check_sphere_coll(point &center, float sphere_radius) const { // sphere-sphere intersection

	float const rsum(radius + sphere_radius);
	if (!dist_less_than(center, pos, rsum)) return 0;
	vector3d const normal((center - pos).get_norm());
	center = pos + normal*(1.001*rsum); // move center out along sphere normal so it doesn't intersect
	return 1;
}


void scenery_obj::shift_by(vector3d const &vd) {
		
	pos   += vd;
	pos.z -= dz; // reset to original z value
	dz     = 0.0;
}


void scenery_obj::gen_spos(int x, int y, int use_xy) {

	if (use_xy) {
		float const px(get_xval(x) + 0.5*DX_VAL*rand2d()), py(get_yval(y) + 0.5*DY_VAL*rand2d());
		pos.assign(px, py, interpolate_mesh_zval(px, py, 0.0, 1, 1));
		return;
	}
	do {
		pos.x = (X_SCENE_SIZE-2.0*DX_VAL)*signed_rand_float2();
		pos.y = (Y_SCENE_SIZE-2.0*DY_VAL)*signed_rand_float2();
		pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
	} while (pos.z <= zmin);
}


float get_new_zval(point const &pos, int x1, int y1, int x2, int y2) {

	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
	if (xpos < x1 || ypos < y1 || xpos > x2 || ypos > y2 || point_outside_mesh(xpos, ypos)) return pos.z;
	return interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
}


bool scenery_obj::update_zvals(int x1, int y1, int x2, int y2) {

	float const delta_z(get_new_zval(pos, x1, y1, x2, y2) - pos.z);
	if (fabs(delta_z) < 1.0E-6) return 0; // no significant change
	dz     = delta_z;
	pos.z += dz;
	return 1;
}


bool scenery_obj::in_camera_view(float brad, vector3d const &xlate) const {

	return sphere_in_camera_view(pos+xlate, ((brad == 0.0) ? radius : brad), 0);
}


bool scenery_obj::is_visible(bool shadow_only, float bradius, vector3d const &xlate) const {

	if (shadow_only ? !is_over_mesh(pos+xlate) : !in_camera_view(bradius, xlate)) return 0;
	return (shadow_only || !skip_uw_draw(pos+xlate, radius));
}


float scenery_obj::get_shadowed_color(point const &p, float eff_radius) const { // not used on rock_shapes

	if (world_mode != WMODE_GROUND) return 1.0;
	return (is_visible_to_light_cobj(p, get_light(), eff_radius, coll_id, 0) ? 1.0 : SHADOW_VAL);
}


float scenery_obj::get_size_scale(float dist_to_camera, float scale_val, float scale_exp) const {
	return ((scale_val == 0.0) ? 1.0 : min(1.0f, pow(scale_val/dist_to_camera, scale_exp)));
}


colorRGBA scenery_obj::get_atten_color(colorRGBA c, vector3d const &xlate) const {

	water_color_atten_at_pos(c, (pos + xlate + point(0.0, 0.0, radius)));
	return c;
}


void scenery_obj::remove_cobjs() {

	remove_reset_coll_obj(coll_id);
}


void rock_shape3d::create(int x, int y, bool use_xy) {

	int rs_rock(rand2());
	gen_spos(x, y, use_xy);
	gen_rock(48, 0.05/tree_scale, rs_rock, (rand2()&1));
	radius = 0.0;
	
	for (unsigned i = 0; i < points.size(); ++i) { // calculate radius
		radius = max(radius, points[i].mag_sq()); // check if correct - might need pos or average/center calculation
	}
	radius = sqrt(radius);
	pos.z += 0.1*radius;
}


struct edge {

	unsigned face;
	bool dir;
	edge(unsigned f=0, bool d=0) : face(f), dir(d) {}
};


class edge_seen_set {

	set<pair<unsigned, unsigned> > seen;

public:
	bool find(unsigned a, unsigned b) const {
		assert(a != b);
		return (seen.find(make_pair(min(a, b), max(a,b))) != seen.end());
	}
	void insert(unsigned a, unsigned b) {
		assert(a != b);
		seen.insert(make_pair(min(a, b), max(a,b)));
	}
};


void rock_shape3d::gen_rock(unsigned nverts, float size, int rand_seed, int type) {

	set_rand2_state(rand_seed, 10423232);
	nverts    = max(nverts, 4U);
	color     = LT_BROWN;
	tex_scale = 5.0;
	tid       = DARK_ROCK_TEX;

	if (type == 0) { // this doesn't really look like a rock
		alloc_shape(nverts, (2*nverts - 4), 0);
		unsigned face_counter(4);

		// create default 4x4 prism shape
		// vertices
		for (unsigned i = 0; i < 4; ++i) {
			points[i] = signed_rand_vector2(size);
		}

		// face vertices
		unsigned const fv[4][3] = {{2,1,0}, {2,0,3}, {1,2,3}, {0,1,3}};

		for (unsigned i = 0; i < 4; ++i) {
			for (unsigned j = 0; j < 3; ++j) {
				faces[i].v[j] = fv[i][j];
			}
		}
		for (unsigned i = 4; i < points.size(); ++i) {
			point center;
			unsigned const face_id(rand2()%face_counter);
			get_face_normal(face_id);
			get_triangle_center(center, face_id);

			for (unsigned j = 0; j < 3; ++j) {
				points[i][j] = center[j] + faces[face_id].norm[j]*size*rand2d();
			}
			add_vertex(i, face_id, face_counter);
		}
	}
	else if (type == 1) {
		alloc_shape(nverts, 10*nverts, 0); // not sure how many faces yet

		for (unsigned i = 0; i < nverts; ++i) {
			points[i] = gen_rand_vector2(size);
		}
		unsigned face(0);
		vector<unsigned char> used(nverts, 0);
		edge_seen_set edges_seen;
		deque<edge> edges; // incomplete faces

		for (unsigned cv = 0; cv < nverts; ++cv) { // is this outer loop necessary?
			if (used[cv]) continue; // finished with this vertex
			unsigned imin(0);
			float dmin(0.0);

			for (unsigned i = 0; i < nverts; ++i) { // find closest point to cv
				if (i == cv) continue;
				float const d(p2p_dist_sq(points[cv], points[i]));
				if (dmin == 0.0 || d < dmin) {dmin = d; imin = i;}
			}
			assert(dmin != 0.0);
			assert(!edges_seen.find(cv, imin));
			used[cv]   = 1;
			used[imin] = 1;

			for (unsigned d = 0; d < 2; ++d) { // two faces on every new edge
				assert(face < faces.size());
				faces[face].v[0] = cv; // start a new face
				faces[face].v[1] = imin;
				edges.push_back(edge(face++, (d != 0)));
			}
			while (!edges.empty()) {
				edge e(edges.back());
				edges.pop_back();
				assert(e.face < faces.size());
				unsigned *v(faces[e.face].v);
				float dmin(0.0);

				for (unsigned i = 0; i < nverts; ++i) {
					if (i == v[0] || i == v[1]) continue;

					if (dmin > 0.0) { // not quite right
						// points[i] must lie on the edir side of the line (points[v[1]] - points[v[0]])
						float const dp(dot_product((points[i] - points[v[0]]), (points[v[0]] - points[v[1]])));
						if ((dp < 0.0)^e.dir) continue;
					}
					vector3d const A(points[v[0]], points[i]), B(points[v[1]], points[i]);
					point2d<float> const dist(A.mag(), B.mag());
					float const d(dist.mag_sq() - 0.05*size*fabs(dot_product(points[v[1]].get_norm(), cross_product(A, B).get_norm())));
					if (dmin == 0.0 || d < dmin) {dmin = d; v[2] = i;}
				}
				assert(dmin != 0.0);
				used[v[2]] = 1;
				
				for (unsigned d = 0; d < 2; ++d) {
					if (!edges_seen.find(v[d], v[2])) { // new edge that has no opposite face
						edges_seen.insert(v[d], v[2]);
						assert(face < faces.size());
						faces[face].v[0] = v[d]; // start a new face
						faces[face].v[1] = v[2];
						float const dp(dot_product((points[v[!d]] - points[v[d]]), (points[v[d]] - points[v[2]])));
						edges.push_back(edge(face++, (dp > 0.0)));
					}
				}
			}
		}
		faces.resize(face);
		//cout << "nverts = " << nverts << ", nfaces = " << faces.size() << endl;
	}
	else {
		assert(0);
	}
	for (unsigned i = 0; i < faces.size(); ++i) {
		faces[i].color_id = 0;
	}
	gen_face_normals();
}

void rock_shape3d::add_cobjs() {

	coll_id = add_coll_sphere(pos, 0.5*radius, cobj_params(0.9, BROWN, 0, 0, rock_collision, 0, DARK_ROCK_TEX));
}

bool rock_shape3d::do_impact_damage(point const &pos_, float radius_) {

	if (radius < 0.02 || p2p_dist(pos_, pos) > (0.75*radius_ + 0.5*radius)) return 0;
	float const size_scale(pow(0.99f, fticks));
	radius *= size_scale; // chip off parts of the rock to make it smaller
	scale  *= size_scale;
	remove_cobjs();
	add_cobjs();
	return 1;
}

void rock_shape3d::draw(bool shadow_only, vector3d const &xlate) const { // Note: assumes texture is already setup

	if (!is_visible(shadow_only, 0.0, xlate)) return;
	(shadow_only ? WHITE : get_atten_color(color*get_shadowed_color(pos, 0.5*radius), xlate)).set_for_cur_shader();
	draw_using_vbo();
}

void rock_shape3d::draw_using_vbo() const {

	unsigned const vert_size(3*faces.size());

	if (vbo == 0) {
		vector<vert_norm_tc> verts;
		verts.reserve(vert_size);
		get_triangle_verts(verts);
		create_vbo_and_upload(vbo, verts, 0, 0);
	}
	else {
		bind_vbo(vbo);
	}
	draw_verts<vert_norm_tc>(NULL, vert_size, GL_TRIANGLES);
	bind_vbo(0);
}


bool rock_shape3d::update_zvals(int x1, int y1, int x2, int y2) {

	if (!scenery_obj::update_zvals(x1, y1, x2, y2)) return 0;
	clear_vbo(); // clear and recreate points if pos changes
	return 1;
}


void rock_shape3d::clear_vbo() {
	delete_and_zero_vbo(vbo);
}


upsurface *surface_cache::get_surface(bool fixed_sz_rock_cache) {

	seed_pair sp;
	
	if (fixed_sz_rock_cache) { // cache size is 256
		sp = seed_pair((global_rand_gen.rseed1 & 15), (global_rand_gen.rseed2 & 15));
	}
	else {
		sp = seed_pair(global_rand_gen.rseed1, global_rand_gen.rseed2);
	}
	surface_map::const_iterator it(scache.find(sp));
		
	if (it != scache.end()) {
		assert(it->second);
		it->second->inc_ref();
		return it->second;
	}
	upsurface *surface(new upsurface);
	scache[sp] = surface;
	surface->inc_ref();
	return surface;
}

void surface_cache::clear() {

	for (surface_map::const_iterator i = scache.begin(); i != scache.end(); ++i) {
		assert(i->second);
		delete i->second;
	}
	scache.clear();
}

void surface_cache::clear_unref() {

	vector<seed_pair> to_erase;

	for (surface_map::iterator i = scache.begin(); i != scache.end(); ++i) {
		assert(i->second);
			
		if (i->second->unrefed()) {
			delete i->second;
			to_erase.push_back(i->first);
		}
	}
	for (vector<seed_pair>::const_iterator i = to_erase.begin(); i != to_erase.end(); ++i) {
		scache.erase(*i);
	}
}


surface_cache surface_rock_cache;


void surface_rock::create(int x, int y, int use_xy, vbo_vnt_block_manager_t &vbo_manager, bool fixed_sz_rock_cache) {

	gen_spos(x, y, use_xy);
	radius  = rand_uniform2(0.1, 0.2)*rand_float2()/tree_scale;
	dir     = signed_rand_vector2_norm();
	surface = surface_rock_cache.get_surface(fixed_sz_rock_cache);
	assert(surface);

	if (surface->ssize == 0) { // not inited
		surface->gen(0.5, rand_uniform2(0.5, 5.0), 10, rand_uniform2(0.5, 2.0));
		surface->setup(ROCK_NDIV, 0.0, 0);
		surface->setup_draw_sphere(all_zeros, rand_uniform2(0.25, 1.0), 0.0, ROCK_NDIV, NULL);
		surface->calc_rmax();
	}
	scale = radius/surface->rmax;
	surface->sd.get_quad_points(vbo_manager.get_pts_vector_for_adding());
	vbo_mgr_ix = vbo_manager.get_offset_for_last_points_added();
}

void surface_rock::add_cobjs() {
	coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, BROWN, 0, 0, rock_collision, 1, ROCK_SPHERE_TEX));
}

void surface_rock::draw(float sscale, bool shadow_only, vector3d const &xlate, float scale_val, vbo_vnt_block_manager_t &vbo_manager) const {

	assert(surface);
	if (!is_visible(shadow_only, 0.0, xlate)) return;
	colorRGBA const color(shadow_only ? WHITE : get_atten_color(WHITE, xlate)*get_shadowed_color(pos+xlate, radius));
	float const dist(distance_to_camera(pos+xlate));

	if (!shadow_only && 2*get_pt_line_thresh()*radius < dist) { // draw as point
		tree_scenery_pld.add_textured_pt(pos+xlate, color, ROCK_SPHERE_TEX);
		return;
	}
	color.set_for_cur_shader();
	fgPushMatrix();
	translate_to(pos);
	uniform_scale(scale*get_size_scale(dist, scale_val));
	rotate_into_plus_z(dir);
	assert(vbo_mgr_ix >= 0);
	vbo_manager.render_range(vbo_mgr_ix, vbo_mgr_ix+1);
	fgPopMatrix();
}

void surface_rock::update_points_vbo(vbo_vnt_block_manager_t &vbo_manager) {

	assert(vbo_mgr_ix >= 0);
	vector<vert_norm_tc> points;
	surface->sd.get_quad_points(points);
	vbo_manager.update_range(points, WHITE, vbo_mgr_ix, vbo_mgr_ix+1); // color is unused
}

bool surface_rock::update_zvals(int x1, int y1, int x2, int y2, vbo_vnt_block_manager_t &vbo_manager) {

	if (!scenery_obj::update_zvals(x1, y1, x2, y2)) return 0;
	if (vbo_mgr_ix >= 0) {update_points_vbo(vbo_manager);}
	remove_cobjs();
	add_cobjs();
	return 1;
}

void surface_rock::destroy() {

	//delete surface;
	if (surface) surface->dec_ref();
	vbo_mgr_ix = -1;
	surface    = NULL;
	scenery_obj::destroy();
}


void s_rock::create(int x, int y, int use_xy) {

	for (unsigned i = 0; i < 3; ++i) {
		scale[i] = rand_uniform2(0.8, 1.3);
	}
	gen_spos(x, y, use_xy);
	size   = 0.02*rand_uniform2(0.2, 0.8)/tree_scale;
	if ((rand2()&3) == 0) size *= rand_uniform2(1.2, 8.0);
	dir    = signed_rand_vector2_norm();
	angle  = rand_uniform2(0.0, 360.0);
	radius = size*(scale.x + scale.y + scale.z)/3.0;
	pos.z += radius*rand_uniform2(-0.1, 0.25);
}

void s_rock::add_cobjs() {
	coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, BROWN, 0, 0, rock_collision, 1, ROCK_SPHERE_TEX));
}

void s_rock::draw(float sscale, bool shadow_only, vector3d const &xlate, float scale_val) const {

	float const rmax(1.3*radius);
	if (!is_visible(shadow_only, rmax, xlate)) return;
	colorRGBA const color(shadow_only ? WHITE : get_atten_color(WHITE, xlate)*get_shadowed_color(pos+xlate, rmax));
	float const dist(distance_to_camera(pos+xlate));

	if (!shadow_only && 2*get_pt_line_thresh()*radius < dist) { // draw as point
		tree_scenery_pld.add_textured_pt(pos+xlate, color, ROCK_SPHERE_TEX);
		return;
	}
	color.set_for_cur_shader();
	int const ndiv(max(4, min(N_SPHERE_DIV, (shadow_only ? get_smap_ndiv(radius) : int(sscale*radius/dist)))));
	fgPushMatrix();
	translate_to(pos);
	rotate_about(angle, dir);
	scale_by(size*get_size_scale(dist, scale_val)*scale);
	draw_sphere_vbo_raw(ndiv, 1);
	fgPopMatrix();
}


void voxel_rock::create(int x, int y, int use_xy) {

	gen_spos(x, y, use_xy);
	radius = 0.2*rand_uniform2(0.5, 1.0)*rand_float2()/tree_scale;
	rseed  = rand2();
}

void voxel_rock::build_model() {

	float const gen_radius(gen_voxel_rock(model, all_zeros, 1.0, ROCK_VOX_SZ, 1, rseed)); // will be translated to pos and scaled by radius during rendering
	assert(gen_radius > 0.0);
	radius /= gen_radius;
}

void voxel_rock::add_cobjs() {
	coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, LT_GRAY, 0, 0, rock_collision, 1, get_tid()));
}

void voxel_rock::draw(float sscale, bool shadow_only, vector3d const &xlate, float scale_val, shader_t &s, bool use_model_texgen) {

	assert(radius > 0.0);
	if (!is_visible(shadow_only, radius, xlate)) return;
	unsigned const lod_level = 0;
	colorRGBA const color(shadow_only ? WHITE : get_atten_color(WHITE, xlate)*get_shadowed_color(pos+xlate, radius));
	color.set_for_cur_shader();
	fgPushMatrix();
	translate_to(pos);
	uniform_scale(radius*get_size_scale(distance_to_camera(pos+xlate), scale_val));
	
	if (use_model_texgen) {
		model.setup_tex_gen_for_rendering(s);
	}
	else {
		select_texture(get_tid());
	}
	model.core_render(s, lod_level, shadow_only, 1); // disable view frustum culling because it's incorrect (due to transform matrices)
	fgPopMatrix();
}

void voxel_rock::destroy() {

	model.clear();
	scenery_obj::destroy();
}


void s_log::shift_by(vector3d const &vd) {

	scenery_obj::shift_by(vd);
	pt2   += vd;
	pt2.z -= dz;
	dz     = 0.0;
}

int s_log::create(int x, int y, int use_xy, float minz) {

	gen_spos(x, y, use_xy);
	radius  = rand_uniform2(0.003, 0.008)/tree_scale;
	radius2 = rand_uniform2(0.9*radius, 1.1*radius);
	length  = rand_uniform2(max(0.03/tree_scale, 4.0*radius), min(0.15/tree_scale, 20.0*radius));
	dir     = signed_rand_vector2_norm();
	dir.x  *= length;
	dir.y  *= length;
	pt2.x   = pos.x + dir.x;
	pt2.y   = pos.y + dir.y;

	if (world_mode == WMODE_GROUND && (pt2.x > X_SCENE_SIZE-DX_VAL || pt2.x < -X_SCENE_SIZE+DX_VAL || pt2.y > Y_SCENE_SIZE-DY_VAL || pt2.y < -Y_SCENE_SIZE+DY_VAL)) {
		return 0; // off the end of the mesh
	}
	pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1) + rand_uniform2(0.7, 0.99)*radius;
	pt2.z = interpolate_mesh_zval(pt2.x, pt2.y, 0.0, 1, 1) + rand_uniform2(0.7, 0.99)*radius2;
	if (max(pos.z, pt2.z) < minz)       return 0;
	if (pos.z <= zmin || pt2.z <= zmin) return 0; // bad z value
	dir.z  = (pt2.z - pos.z);
	length = dir.mag(); // recalculate
	dir   /= -length; // something is backwards
	type   = (char)get_tree_type_from_height(max(pos.z, pt2.z), global_rand_gen);
	return (type >= 0);
}

void s_log::add_cobjs() {
	coll_id = add_coll_cylinder(pos, pt2, radius, radius2, cobj_params(0.8, BROWN, 0, 0, NULL, 0, get_tid()));
}

void s_log::draw(float sscale, bool shadow_only, vector3d const &xlate, float scale_val) const {

	float const sz(max(length, max(radius, radius2)));
	point const center((pos + pt2)*0.5 + xlate);
	if (type < 0 || (shadow_only ? !is_over_mesh(center) : !in_camera_view(sz, xlate))) return;
	colorRGBA const color(shadow_only ? WHITE : get_tree_trunk_color(type, 0)*get_shadowed_color(center, sz));
	float const dist(distance_to_camera(center));

	if (!shadow_only && get_pt_line_thresh()*(radius + radius2) < dist) { // draw as line
		tree_scenery_pld.add_textured_line(pos+xlate, pt2+xlate, color, get_tid());
		return;
	}
	color.set_for_cur_shader();
	int const ndiv(max(3, min(N_CYL_SIDES, (shadow_only ? get_smap_ndiv(2.0*radius) : int(2.0*sscale*radius/dist)))));
	fgPushMatrix();
	translate_to(pos);
	rotate_by_vector(dir);
	if (!shadow_only) select_texture(TREE_END_TEX);
	draw_circle_normal(0.0, radius,  ndiv, 1, 0.0);
	draw_circle_normal(0.0, radius2, ndiv, 0, length);
	if (!shadow_only) select_texture(get_tid());
	draw_cylin_fast(radius, radius2, length, ndiv, 1);
	fgPopMatrix();
}

bool s_log::update_zvals(int x1, int y1, int x2, int y2) {

	float const orig_pz(pos.z);
	if (!scenery_obj::update_zvals(x1, y1, x2, y2)) return 0;
	pt2.z += (pos.z - orig_pz); // apply the same offset to pt2 even though it might be at a different mesh location
	return 1;
}


int s_stump::create(int x, int y, int use_xy, float minz) {

	gen_spos(x, y, use_xy);
	radius  = rand_uniform2(0.005, 0.01)/tree_scale;
	radius2 = rand_uniform2(0.8*radius, radius);
	pos.z  -= 2.0*radius;
	if (pos.z < minz) return 0;
	height  = rand_uniform2(0.01/tree_scale, min(0.05/tree_scale, 4.0*radius)) + 0.015;
		
	if ((rand2()&3) == 0) {
		height  *= rand_uniform2(1.0, 5.0); // larger stump = upright dead tree
		radius  *= 1.5;
		radius2 *= 1.3;
	}
	type = (char)get_tree_type_from_height(pos.z, global_rand_gen);
	return (type >= 0);
}

void s_stump::add_cobjs() {
	coll_id = add_coll_cylinder(pos-point(0.0, 0.0, 0.2*height), pos+point(0.0, 0.0, height), radius, radius2, cobj_params(0.8, BROWN, 0, 0, NULL, 0, get_tid()));
}

bool s_stump::check_sphere_coll(point &center, float sphere_radius) const {
	return sphere_vert_cylin_intersect(center, sphere_radius, cylinder_3dw(pos-point(0.0, 0.0, 0.2*height), pos+point(0.0, 0.0, height), radius, radius));
}

void s_stump::draw(float sscale, bool shadow_only, vector3d const &xlate, float scale_val) const {

	float const sz(max(height, max(radius, radius2)));
	point const center(pos + point(0.0, 0.0, 0.5*height) + xlate);
	if (type < 0 || (shadow_only ? !is_over_mesh(center) : !in_camera_view(sz, xlate))) return;
	colorRGBA const color(shadow_only ? WHITE : get_tree_trunk_color(type, 0)*get_shadowed_color(center, sz));
	float const dist(distance_to_camera(center));

	if (!shadow_only && get_pt_line_thresh()*(radius + radius2) < dist) { // draw as line
		tree_scenery_pld.add_textured_line((pos+xlate - point(0.0, 0.0, 0.2*height)), (pos+xlate + point(0.0, 0.0, height)), color, get_tid());
		return;
	}
	color.set_for_cur_shader();
	int const ndiv(max(3, min(N_CYL_SIDES, (shadow_only ? get_smap_ndiv(2.2*radius) : int(2.2*sscale*radius/dist)))));
	fgPushMatrix();
	translate_to(pos - point(0.0, 0.0, 0.2*height));
	if (!shadow_only) select_texture(TREE_END_TEX);
	draw_circle_normal(0.0, radius2, ndiv, 0, 1.2*height);
	if (!shadow_only) select_texture(get_tid());
	draw_cylin_fast(radius, radius2, 1.2*height, ndiv, 1);
	fgPopMatrix();
}


int s_plant::create(int x, int y, int use_xy, float minz, vbo_vnc_block_manager_t &vbo_manager) {

	vbo_mgr_ix = -1;
	type   = rand2()%NUM_PLANT_TYPES;
	gen_spos(x, y, use_xy);
	if (pos.z < minz) return 0;
	if (get_rel_height(pos.z, -zmax_est, zmax_est) > 0.62) return 0; // altitude too high for plants
	radius = rand_uniform2(0.0025, 0.0045)/tree_scale;
	height = rand_uniform2(0.2, 0.4)/tree_scale + 0.025;
	gen_points(vbo_manager);
	return 1;
}

void s_plant::create2(point const &pos_, float height_, float radius_, int type_, int calc_z, vbo_vnc_block_manager_t &vbo_manager) {

	vbo_mgr_ix = -1;
	type   = abs(type_)%NUM_PLANT_TYPES;
	pos    = pos_;
	radius = radius_;
	height = height_;
	if (calc_z) pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
	gen_points(vbo_manager);
}

void s_plant::add_cobjs() {

	point cpos(pos), cpos2(pos), bpos(pos);
	float const wscale(radius*tree_scale/0.004), r2(radius+0.07*wscale*(height+0.03));
	cpos.z  += height;
	cpos2.z += 3.0*height/(36.0*height + 4.0);
	bpos.z  -= 0.1*height;
	coll_id  = add_coll_cylinder(bpos, cpos, radius, 0.0, cobj_params(0.4, pltype[type].stemc, 0, 0, NULL, 0, WOOD_TEX)); // trunk
	if (!no_leaves) {coll_id2 = add_coll_cylinder(cpos2, cpos, r2, radius, cobj_params(0.4, pltype[type].leafc, 0, 0, NULL, 0, pltype[type].tid));} // leaves
}

bool s_plant::check_sphere_coll(point &center, float sphere_radius) const { // used in tiled terrain mode, ignores no_leaves
	return sphere_vert_cylin_intersect(center, sphere_radius, cylinder_3dw(pos-point(0.0, 0.0, 0.1*height), pos+point(0.0, 0.0, height), radius, radius));
}

void s_plant::create_leaf_points(vector<vert_norm> &points) const {

	// Note: could scale leaves for different plant types differently in x vs. y to allow for non-square textures (tighter bounds = lower fillrate)
	float const wscale(250.0*radius*tree_scale), theta0((int(1.0E6*height)%360)*TO_RADIANS);
	unsigned const nlevels(unsigned(36.0*height*tree_scale)), nrings(3);
	float rdeg(30.0);
	points.resize(4*nlevels*nrings);

	for (unsigned j = 0, ix = 0; j < nlevels; ++j) { // could do the same optimizations as the high detail pine tree
		float const sz(0.07*(height + 0.03/tree_scale)*((nlevels - j + 3.0)/(float)nlevels));
		float const z((j + 3.0)*height/(nlevels + 4.0));

		for (unsigned k = 0; k < nrings; ++k) {
			float const theta(TWO_PI*(3.3*j + 0.2*k) + theta0);
			int const val(int(((int(1.0E6*height))*(5463*j + 537879*k))%301));
			rdeg += 0.01*(val - 150);
			add_rotated_quad_pts(&points.front(), ix, theta, z, pos, sz*wscale, sz*rdeg/45.0);
		}
	}
}

void s_plant::gen_points(vbo_vnc_block_manager_t &vbo_manager) {

	if (vbo_mgr_ix >= 0) return; // already generated
	vector<vert_norm> &pts(vbo_manager.temp_points);
	create_leaf_points(pts);
	vbo_mgr_ix = vbo_manager.add_points_with_offset(pts, pltype[type].leafc);
	no_leaves  = 0;

	if (pltype[type].berryc.A != 0.0) { // create berries
		rand_gen_t rgen;

		for (unsigned i = 0; i < pts.size(); i += 4) { // iterate over each leaf quad
			point const start(0.5*(pts[i+0].v + pts[i+1].v)), end(0.5*(pts[i+2].v + pts[i+3].v)); // along centerline
			unsigned const num((rgen.rand()&3) + 7);
			float const t0(rgen.rand_uniform(0.1, 0.2)), dt((rgen.rand_uniform(0.8, 0.9) - t0)/(num-1));
			float const err(0.1*dt*p2p_dist(start, end));
			float t(t0);

			for (unsigned n = 0; n < num; ++n, t += dt) {
				berries.push_back(start + t*(end - start) + rgen.signed_rand_vector(err));
			}
		}
	}
}

// to be called when the plant is translated or zval changes
void s_plant::update_points_vbo(vbo_vnc_block_manager_t &vbo_manager) {

	assert(vbo_mgr_ix >= 0);
	create_leaf_points(vbo_manager.temp_points);
	vbo_manager.update_range(vbo_manager.temp_points, pltype[type].leafc, vbo_mgr_ix, vbo_mgr_ix+1);
}

bool s_plant::update_zvals(int x1, int y1, int x2, int y2, vbo_vnc_block_manager_t &vbo_manager) {

	float orig_z(pos.z);
	if (!scenery_obj::update_zvals(x1, y1, x2, y2)) return 0;
	for (vector<vert_wrap_t>::iterator i = berries.begin(); i != berries.end(); ++i) {i->v.z += (pos.z - orig_z);}

	if (vbo_mgr_ix >= 0) { // leaf pos changed - VBO data needs to be updated
		//disable_leaves(); // remove the leaves (simpler and more efficient)
		update_points_vbo(vbo_manager); // regenerate leaf points and re-upload VBO sub-data (slower)
	}
	remove_cobjs();
	add_cobjs();
	return 1;
}

bool s_plant::is_shadowed() const {

	if (world_mode != WMODE_GROUND) return 0;
	int const light(get_light());

	for (unsigned i = 0; i < 3; ++i) {
		point p(pos);
		p.z += 0.5*i*height;
		if (is_visible_to_light_cobj(p, light, (radius + height), coll_id, 0)) return 0;
	}
	return 1;
}

void s_plant::draw_stem(float sscale, bool shadow_only, vector3d const &xlate) const {

	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	if (shadow_only ? !is_over_mesh(pos+xlate) : !sphere_in_camera_view(pos2, (height + radius), 0)) return;
	bool const shadowed(shadow_only ? 0 : is_shadowed());
	colorRGBA color(pltype[type].stemc*(shadowed ? SHADOW_VAL : 1.0));
	float const dist(distance_to_camera(pos2));

	if (!shadow_only && 2*get_pt_line_thresh()*radius < dist) { // draw as line
		tree_scenery_pld.add_textured_line((pos+xlate - point(0.0, 0.0, 0.1*height)), (pos+xlate + point(0.0, 0.0, 0.75*height)), color, WOOD_TEX);
	}
	else {
		int const ndiv(max(3, min(N_CYL_SIDES, (shadow_only ? get_smap_ndiv(2.0*radius) : int(2.0*sscale*radius/dist)))));
		color.set_for_cur_shader();
		draw_fast_cylinder((pos - point(0.0, 0.0, 0.1*height)), (pos + point(0.0, 0.0, height)), radius, 0.0, ndiv, 1);
	}
}

void s_plant::draw_leaves(shader_t &s, vbo_vnc_block_manager_t &vbo_manager, bool shadow_only, vector3d const &xlate) const {

	if (no_leaves) return;
	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	if (shadow_only ? !is_over_mesh(pos2) : !sphere_in_camera_view(pos2, 0.5*(height + radius), 0)) return;
	bool const shadowed(shadow_only ? 0 : is_shadowed());
	if (shadowed) {s.add_uniform_float("normal_scale", 0.0);}
	select_texture((draw_model == 0) ? pltype[type].tid : WHITE_TEX); // could pre-bind textures and select using shader int, but probably won't improve performance
	assert(vbo_mgr_ix >= 0);
	vbo_manager.render_range(vbo_mgr_ix, vbo_mgr_ix+1);
	if (shadowed) {s.add_uniform_float("normal_scale", 1.0);}
}

void s_plant::draw_berries(shader_t &s, vector3d const &xlate) const {

	if (berries.empty()) return;
	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	if (!sphere_in_camera_view(pos2, 0.5*(height + radius), 0)) return;
	float const size_scale(get_pt_line_thresh()*radius/distance_to_camera(pos2));
	if (size_scale < 1.2) return; // too small/far away
	int const ndiv(max(4, min(16, int(2.0*size_scale))));
	s.set_cur_color(pltype[type].berryc);
	fgPushMatrix();
	uniform_scale(0.25*radius);
	s.add_uniform_float("vertex_offset_scale", 1.0/(0.25*radius)); // enable
	int const loc(s.get_attrib_loc("vertex_offset", 0));
	shader_float_matrix_uploader<3,1>::enable(loc, 1, &berries.front().v.x); // FIXME_VBO
	draw_sphere_vbo_raw(ndiv, 0, 0, berries.size()); // ndiv=16, untextured
	shader_float_matrix_uploader<3,1>::disable(loc);
	s.add_uniform_float("vertex_offset_scale", 0.0); // disable
	fgPopMatrix();
}

void s_plant::remove_cobjs() {

	remove_reset_coll_obj(coll_id2);
	scenery_obj::remove_cobjs();
}

void s_plant::destroy() {

	berries.clear();
	remove_cobjs();
	scenery_obj::destroy(); // will remove coll_id twice, which is OK
}


// ************ SCENERY OBJECT INTERFACE/WRAPPERS/DRIVERS ************


template<typename T> void draw_scenery_vector(vector<T> &v, float sscale, bool shadow_only, vector3d const &xlate, float scale_val) {
	for (unsigned i = 0; i < v.size(); ++i) {v[i].draw(sscale, shadow_only, xlate, scale_val);}
}

template<typename T> void add_scenery_vector_cobjs(vector<T> &v) {
	for (unsigned i = 0; i < v.size(); ++i) {
		if (is_over_mesh(v[i].get_pos())) {v[i].add_cobjs();}
	}
}

template<typename T> bool check_scenery_vector_sphere_coll(vector<T> const &v, point &center, float radius) {
	bool coll(0);
	for (unsigned i = 0; i < v.size(); ++i) {coll |= v[i].check_sphere_coll(center, radius);}
	return coll;
}

template<typename T> void shift_scenery_vector(vector<T> &v, vector3d const &vd) {
	for (unsigned i = 0; i < v.size(); ++i) {v[i].shift_by(vd);}
}

template<typename T> void free_scenery_vector(vector<T> &v) {
	for (unsigned i = 0; i < v.size(); ++i) {v[i].destroy();}
}

template<typename T> void update_scenery_zvals_vector(vector<T> &v, int x1, int y1, int x2, int y2, bool &updated) {
	
	for (unsigned i = 0; i < v.size(); ++i) { // zval has change, remove and re-add cobjs
		if (v[i].update_zvals(x1, y1, x2, y2)) {
			v[i].remove_cobjs();
			v[i].add_cobjs();
			updated = 1;
		}
	}
}


void scenery_group::clear_vbos() {
	
	for (unsigned i = 0; i < rock_shapes.size(); ++i) {
		rock_shapes[i].clear_vbo();
	}
	for (unsigned i = 0; i < voxel_rocks.size(); ++i) {
		voxel_rocks[i].free_context();
	}
	plant_vbo_manager.clear_vbo();
	rock_vbo_manager.clear_vbo();
}

void scenery_group::clear() {

	free_scenery();
	rock_shapes.clear();
	surface_rocks.clear();
	voxel_rocks.clear();
	rocks.clear();
	logs.clear();
	stumps.clear();
	plants.clear();
	clear_vbos();
	generated = 0;
}

void scenery_group::free_scenery() {

	free_scenery_vector(rock_shapes);
	free_scenery_vector(surface_rocks);
	free_scenery_vector(voxel_rocks);
	free_scenery_vector(rocks);
	free_scenery_vector(logs);
	free_scenery_vector(stumps);
	free_scenery_vector(plants);
}

void scenery_group::add_cobjs() {

	add_scenery_vector_cobjs(rock_shapes);
	add_scenery_vector_cobjs(surface_rocks);
	add_scenery_vector_cobjs(voxel_rocks);
	add_scenery_vector_cobjs(rocks);
	add_scenery_vector_cobjs(logs);
	add_scenery_vector_cobjs(stumps);
	add_scenery_vector_cobjs(plants);
}


bool scenery_group::check_sphere_coll(point &center, float radius) const {

	bool coll(0);
	coll |= check_scenery_vector_sphere_coll(rock_shapes,   center, radius);
	coll |= check_scenery_vector_sphere_coll(surface_rocks, center, radius);
	coll |= check_scenery_vector_sphere_coll(voxel_rocks,   center, radius);
	coll |= check_scenery_vector_sphere_coll(rocks,         center, radius);
	coll |= check_scenery_vector_sphere_coll(logs,          center, radius);
	coll |= check_scenery_vector_sphere_coll(stumps,        center, radius);
	coll |= check_scenery_vector_sphere_coll(plants,        center, radius);
	return coll;
}

void scenery_group::shift(vector3d const &vd) {

	shift_scenery_vector(rock_shapes,   vd);
	shift_scenery_vector(surface_rocks, vd);
	shift_scenery_vector(voxel_rocks,   vd);
	shift_scenery_vector(rocks,         vd);
	shift_scenery_vector(logs,          vd);
	shift_scenery_vector(stumps,        vd);
	shift_scenery_vector(plants,        vd);
}

// update region is inclusive: [x1,x2]x[y1,y2]
bool scenery_group::update_zvals(int x1, int y1, int x2, int y2) { // inefficient, should use spatial subdivision

	assert(x1 <= x2 && y1 <= y2);
	bool updated(0);
	update_scenery_zvals_vector(rock_shapes, x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(voxel_rocks, x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(rocks,       x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(logs,        x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(stumps,      x1, y1, x2, y2, updated);

	for (unsigned i = 0; i < plants.size(); ++i) { // zval has change, remove and re-add cobjs
		updated |= plants[i].update_zvals(x1, y1, x2, y2, plant_vbo_manager); // different signature (takes plant_vbo_manager)
	}
	for (unsigned i = 0; i < surface_rocks.size(); ++i) { // zval has change, remove and re-add cobjs
		updated |= surface_rocks[i].update_zvals(x1, y1, x2, y2, rock_vbo_manager); // different signature (takes rock_vbo_manager)
	}
	return updated;
}

void scenery_group::do_rock_damage(point const &pos, float radius, float damage) {

	for (unsigned i = 0; i < rock_shapes.size(); ++i) {
		if (rock_shapes[i].do_impact_damage(pos, radius)) rock_collision(0, -1, zero_vector, pos, damage, IMPACT);
	}
}

void scenery_group::add_plant(point const &pos, float height, float radius, int type, int calc_z) {

	assert(height > 0.0 && radius > 0.0);
	plants.push_back(s_plant());
	plants.back().create2(pos, height, radius, type, calc_z, plant_vbo_manager);
}

void scenery_group::gen(int x1, int y1, int x2, int y2, float vegetation_, bool fixed_sz_rock_cache) {

	//RESET_TIME;
	unsigned const num_voxel_rock_lod_levels = 1;
	unsigned const smod(max(200U, unsigned(3.321*XY_MULT_SIZE/(tree_scale+1))));
	float const min_stump_z(water_plane_z + 0.010*zmax_est);
	float const min_plant_z(water_plane_z + 0.016*zmax_est);
	float const min_log_z  (water_plane_z - 0.040*zmax_est);
	generated = 1;

	for (int i = y1; i < y2; ++i) {
		for (int j = x1; j < x2; ++j) {
			global_rand_gen.rseed1 = 786433* (i + yoff2) + 196613 *rand_gen_index;
			global_rand_gen.rseed2 = 6291469*(j + xoff2) + 1572869*rand_gen_index;
			int const val(rand2_seed_mix()%smod);
			if (val > 100) continue;
			rand2_mix();
			bool const veg((global_rand_gen.rseed1&127)/128.0 < vegetation_);
			
			if (veg && rand2()%100 < 35) { // Note: numbers below were based on 30% plants
				plants.push_back(s_plant()); // 35%
				if (!plants.back().create(j, i, 1, min_plant_z, plant_vbo_manager)) plants.pop_back();
			}
			else if (val < 5) { // 3.5%
				rock_shapes.push_back(rock_shape3d());
				rock_shapes.back().create(j, i, 1);
			}
			else if (val < 15) { // 7%
				surface_rocks.push_back(surface_rock());
				surface_rocks.back().create(j, i, 1, rock_vbo_manager, fixed_sz_rock_cache);
			}
			else if (USE_VOXEL_ROCKS && val < 35) { // FIXME: too slow, and need special shaders for texturing
				voxel_rocks.push_back(voxel_rock(&voxel_rock_ntg, num_voxel_rock_lod_levels));
				voxel_rocks.back().create(j, i, 1);
			}
			else if (val < 50) { // 24.5%
				rocks.push_back(s_rock());
				rocks.back().create(j, i, 1);
			}
			else if (veg && val < 85) { // 24.5%
				logs.push_back(s_log());
				if (!logs.back().create(j, i, 1, min_log_z)) logs.pop_back();
			}
			else if (veg) {
				stumps.push_back(s_stump()); // 10.5%
				if (!stumps.back().create(j, i, 1, min_stump_z)) stumps.pop_back();
			}
		}
	}
	if (!fixed_sz_rock_cache) {surface_rock_cache.clear_unref();}
	sort(plants.begin(), plants.end()); // sort by type

	if (!voxel_rocks.empty()) {
		RESET_TIME;

		#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < (int)voxel_rocks.size(); ++i) {
			voxel_rocks[i].build_model();
		}
		PRINT_TIME("Gen Voxel Rocks");
	}
	//PRINT_TIME("Gen Scenery");
}

void scenery_group::draw_plant_leaves(shader_t &s, bool shadow_only, vector3d const &xlate) {

	s.set_specular(0.25, 20.0); // a small amount of specular
	plant_vbo_manager.upload();
	plant_vbo_manager.begin_render();

	for (unsigned i = 0; i < plants.size(); ++i) {
		plants[i].draw_leaves(s, plant_vbo_manager, shadow_only, xlate);
	}
	plant_vbo_manager.end_render();
	s.set_specular(0.0, 1.0);
}


void scenery_group::draw_opaque_objects(shader_t &s, bool shadow_only, vector3d const &xlate, bool draw_pld, float scale_val) {

	set_fill_mode();
	select_texture(DARK_ROCK_TEX);

	for (unsigned i = 0; i < rock_shapes.size(); ++i) {
		rock_shapes[i].draw(shadow_only, xlate);
	}
	int const sscale(int((do_zoom ? ZOOM_FACTOR : 1.0)*window_width));
	rock_vbo_manager.upload();
	rock_vbo_manager.begin_render();
	if (!shadow_only) {select_texture(ROCK_SPHERE_TEX);}

	for (unsigned i = 0; i < surface_rocks.size(); ++i) {
		surface_rocks[i].draw(sscale, shadow_only, xlate, scale_val, rock_vbo_manager);
	}
	rock_vbo_manager.end_render();
	draw_scenery_vector(rocks, sscale, shadow_only, xlate, scale_val);

	for (unsigned i = 0; i < voxel_rocks.size(); ++i) {
		voxel_rocks[i].draw(sscale, shadow_only, xlate, scale_val, s, 0); // Note: no model texgen
	}
	draw_scenery_vector(logs,   sscale, shadow_only, xlate, scale_val);
	draw_scenery_vector(stumps, sscale, shadow_only, xlate, scale_val);
	if (!shadow_only) {select_texture(WOOD_TEX);} // plant stems use wood texture

	for (unsigned i = 0; i < plants.size(); ++i) {
		plants[i].draw_stem(sscale, shadow_only, xlate);
	}
	if (!shadow_only) { // no berry shadows
		select_texture(WHITE_TEX); // berries are untextured
		s.set_specular(0.9, 80.0);
		for (unsigned i = 0; i < plants.size(); ++i) {plants[i].draw_berries(s, xlate);}
		s.set_specular(0.0, 1.0);
	}
	if (draw_pld) {tree_scenery_pld.draw_and_clear();}
}


void scenery_group::draw(bool draw_opaque, bool draw_transparent, bool shadow_only, vector3d const &xlate) {

	if (draw_opaque) { // draw stems, rocks, logs, and stumps
		shader_t s;
		s.begin_simple_textured_shader(0.0, !shadow_only); // with lighting, unless shadow_only
		draw_opaque_objects(s, shadow_only, xlate, 1);
		s.end_shader();
	}
	if (draw_transparent && !plants.empty()) { // draw leaves
		shader_t s;
		set_leaf_shader(s, 0.9, 1, 0, 0, shadow_only);
		draw_plant_leaves(s, shadow_only, xlate);
		s.end_shader();
	}
}


scenery_group all_scenery;


void gen_scenery() {

	if (has_scenery2) return; // don't generate scenery if some has already been added
	all_scenery.clear();
	all_scenery = scenery_group(); // really force a clear
	has_scenery = 0;
	if (DISABLE_SCENERY) return;
	has_scenery = 1;
	all_scenery.gen(1, 1, MESH_X_SIZE-1, MESH_Y_SIZE-1, vegetation, 0);
	all_scenery.add_cobjs();
}


void add_plant(point const &pos, float height, float radius, int type, int calc_z) {
	all_scenery.add_plant(pos, height, radius, type, calc_z);
	has_scenery2 = 1;
}

void draw_scenery(bool draw_opaque, bool draw_transparent, bool shadow_only) {
	if (!has_scenery && !has_scenery2) return;
	set_fill_mode();
	all_scenery.draw(draw_opaque, draw_transparent, shadow_only);
}

void add_scenery_cobjs() {
	all_scenery.add_cobjs();
}

void shift_scenery(vector3d const &vd) {
	if (!has_scenery2) return; // dynamically created, not placed
	all_scenery.shift(vd);
}

// update region is inclusive: [x1,x2]x[y1,y2]
bool update_scenery_zvals(int x1, int y1, int x2, int y2) {
	return all_scenery.update_zvals(x1, y1, x2, y2);
}

void free_scenery() {
	all_scenery.free_scenery();
}

void clear_scenery_vbos() {
	all_scenery.clear_vbos();
}

void do_rock_damage(point const &pos, float radius, float damage) {
	all_scenery.do_rock_damage(pos, radius, damage);
}



