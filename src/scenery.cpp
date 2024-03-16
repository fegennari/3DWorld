// 3D World - Scenery Classes Implementations (plants, rocks, logs stumps)
// by Frank Gennari
// 12/5/02

#include "scenery.h"
#include "mesh.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "tree_leaf.h" // for tree_type
#include "tree_3dw.h" // for get_closest_tree_type
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>


bool const ENABLE_PLANT_SHADOWS = 1;
unsigned const ROCK_NDIV        = 24;
unsigned const ROCK_VOX_SZ      = 32;
unsigned const VOX_ROCK_NUM_LOD = 1;
unsigned const NUM_VROCK_MODELS = 100;
float    const SHADOW_VAL       = 0.5;
float    const PT_LINE_THRESH   = 800.0;
float    const LEAFY_PLANT_WIND = 0.25; // less wind than regular plants


colorRGBA const stem_c(0.4, 0.6, 0.2, 1.0);
colorRGBA const leaf_c(0.7, 0.7, 0.7, 1.0);

// tid, stemc, leafc, berryc, leaf_length, leaf_width_base, leaf_width_end
plant_type const pltype[NUM_PLANT_TYPES] = {
	plant_type(MJ_LEAF_TEX, stem_c,   leaf_c, ALPHA0),
	plant_type(PLANT1_TEX,  stem_c,   leaf_c, ALPHA0),
	plant_type(PLANT2_TEX,  stem_c,   leaf_c, PURPLE),
	plant_type(PLANT3_TEX,  stem_c,   leaf_c, colorRGBA(0.9, 0.1, 0.05)),
	plant_type(PLANT4_TEX,  stem_c,   leaf_c, ALPHA0),
	plant_type(COFFEE_TEX,  LT_BROWN, WHITE,  ALPHA0),
	plant_type(GRASS_BLADE_TEX, GREEN, DK_GREEN, ALPHA0, 2.0, 0.2, 0.0) // seaweed
};


int DISABLE_SCENERY(0), has_scenery(0), has_scenery2(0);


extern bool underwater, has_snow;
extern int num_trees, xoff2, yoff2, rand_gen_index, window_width, do_zoom, display_mode, tree_mode, draw_model, DISABLE_WATER, animate2, frame_counter, use_voxel_rocks;
extern float zmin, zmax_est, water_plane_z, tree_scale, vegetation, fticks, ocean_wave_height;
extern pt_line_drawer tree_scenery_pld; // we can use this for plant trunks
extern voxel_params_t global_voxel_params;
extern tree_type tree_types[];


int get_bark_tex_for_tree_type(int type); // for small trees
bool is_pine_tree_type(int type);


inline float get_pt_line_thresh   () {return PT_LINE_THRESH*(do_zoom ? ZOOM_FACTOR : 1.0);}
inline float get_min_water_plane_z() {return (get_water_z_height() - ocean_wave_height);}


bool skip_uw_draw(point const &pos, float radius) {

	// used in tiled terrain mode to skip underwater rocks - otherwise, the fog calculation is incorrect (needs special air/water fog transition handling)
	if (world_mode != WMODE_INF_TERRAIN || DISABLE_WATER || !(display_mode & 0x04)) return 0;
	return (get_camera_pos().z > water_plane_z && (pos.z + radius) < get_min_water_plane_z());
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

bool scenery_obj::check_visible(bool shadow_only, float bradius, point const &p, int level) const {
	// Note: VFC now seems to improve perf for tiled terrain mode plants
	if (world_mode == WMODE_INF_TERRAIN) {return camera_pdu.sphere_visible_test(p, ((bradius == 0.0) ? radius : bradius));}
	return (shadow_only ? is_over_mesh(p) : sphere_in_camera_view(p, ((bradius == 0.0) ? radius : bradius), level));
}
bool scenery_obj::is_visible(bool shadow_only, float bradius, vector3d const &xlate) const {
	if (!check_visible(shadow_only, bradius, pos+xlate)) return 0;
	return (shadow_only || !skip_uw_draw(pos+xlate, radius));
}

float scenery_obj::get_size_scale(float dist_to_camera, float scale_val, float scale_exp) const {
	assert(scale_val >= 0.0);
	return ((scale_val == 0.0) ? 1.0 : min(1.0f, pow(scale_val/dist_to_camera, scale_exp)));
}

colorRGBA scenery_obj::get_atten_color(colorRGBA c, vector3d const &xlate) const {
	water_color_atten_at_pos(c, (pos + xlate + point(0.0, 0.0, radius)));
	return c;
}

void scenery_obj::remove_cobjs() {remove_reset_coll_obj(coll_id);}


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

	vector<uint64_t> seen;

public:
	edge_seen_set(unsigned const sz) {assert(sz <= 64); seen.resize(sz, 0);}
	bool find(unsigned a, unsigned b) const {
		assert(a != b);
		return ((seen[min(a,b)] & ((uint64_t)1 << max(a,b))) != 0);
	}
	void insert(unsigned a, unsigned b) {
		assert(a != b);
		seen[min(a,b)] |= ((uint64_t)1 << max(a,b));
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
		for (unsigned i = 0; i < 4; ++i) {points[i] = signed_rand_vector2(size);}

		// face vertices
		unsigned const fv[4][3] = {{2,1,0}, {2,0,3}, {1,2,3}, {0,1,3}};

		for (unsigned i = 0; i < 4; ++i) {
			for (unsigned j = 0; j < 3; ++j) {faces[i].v[j] = fv[i][j];}
		}
		for (unsigned i = 4; i < points.size(); ++i) {
			point center;
			unsigned const face_id(rand2()%face_counter);
			get_face_normal(face_id);
			get_triangle_center(center, face_id);
			for (unsigned j = 0; j < 3; ++j) {points[i][j] = center[j] + faces[face_id].norm[j]*size*(float)rand2d();}
			add_vertex(i, face_id, face_counter);
		}
	}
	else if (type == 1) {
		alloc_shape(nverts, 10*nverts, 0); // not sure how many faces yet
		for (unsigned i = 0; i < nverts; ++i) {points[i] = gen_rand_vector2(size);}
		unsigned face(0);
		assert(nverts <= 64);
		uint64_t used(0);
		edge_seen_set edges_seen(nverts);
		deque<edge> edges; // incomplete faces

		for (unsigned cv = 0; cv < nverts; ++cv) { // is this outer loop necessary?
			if (used & (1ULL<<cv)) continue; // finished with this vertex
			unsigned imin(0);
			float dmin(0.0);

			for (unsigned i = 0; i < nverts; ++i) { // find closest point to cv
				if (i == cv) continue;
				float const d(p2p_dist_sq(points[cv], points[i]));
				if (dmin == 0.0 || d < dmin) {dmin = d; imin = i;}
			}
			assert(dmin != 0.0);
			assert(!edges_seen.find(cv, imin));
			used |= (1ULL<<cv);
			used |= (1ULL<<imin);

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
				point const &pv0(points[v[0]]), &pv1(points[v[1]]);
				float const pv1_mag(pv1.mag());
				float dmin(0.0);

				for (unsigned i = 0; i < nverts; ++i) {
					if (i == v[0] || i == v[1]) continue;

					if (dmin > 0.0) { // not quite right
						// points[i] must lie on the edir side of the line (points[v[1]] - points[v[0]])
						float const dp(dot_product((points[i] - pv0), (pv0 - pv1)));
						if ((dp < 0.0) ^ e.dir) continue;
					}
					vector3d const A(pv0, points[i]), B(pv1, points[i]), cp(cross_product(A, B));
					float const d((A.mag_sq() + B.mag_sq()) - 0.05f*size*fabs(dot_product(pv1, cp))/(pv1_mag*cp.mag()));
					if (dmin == 0.0 || d < dmin) {dmin = d; v[2] = i;}
				}
				assert(dmin != 0.0);
				used |= (1ULL<<v[2]);
				
				for (unsigned d = 0; d < 2; ++d) {
					if (edges_seen.find(v[d], v[2])) continue; // not a new edge - has opposite face
					edges_seen.insert(v[d], v[2]); // not a new edge - has opposite face
					assert(face < faces.size());
					faces[face].v[0] = v[d]; // start a new face
					faces[face].v[1] = v[2];
					float const dp(dot_product((points[v[!d]] - points[v[d]]), (points[v[d]] - points[v[2]])));
					edges.push_back(edge(face++, (dp > 0.0)));
				}
			} // while
		} // for cv
		faces.resize(face);
		//cout << "nverts = " << nverts << ", nfaces = " << faces.size() << endl;
	}
	else {
		assert(0);
	}
	for (unsigned i = 0; i < faces.size(); ++i) {faces[i].color_id = 0;}
	gen_face_normals();
}

void rock_shape3d::add_cobjs() {
	coll_id = add_coll_sphere(pos, 0.5*radius, cobj_params(0.9, color, 0, 0, rock_collision, 0, DARK_ROCK_TEX));
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

void rock_shape3d::draw(bool shadow_only, bool reflection_pass, vector3d const &xlate) const { // Note: assumes texture is already setup

	if (!is_visible(shadow_only, 0.0, xlate))     return;
	if (reflection_pass && pos.z < water_plane_z) return;
	(shadow_only ? WHITE : get_atten_color(color, xlate)).set_for_cur_shader();
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


p_upsurface surface_cache::get_surface(bool fixed_sz_rock_cache) {

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
	p_upsurface surface(new upsurface);
	scache[sp] = surface;
	surface->inc_ref();
	return surface;
}

void surface_cache::clear_unref() {

	for (surface_map::iterator i = scache.begin(); i != scache.end(); ) { // Note: no ++i
		assert(i->second);
		if (i->second->unrefed()) {scache.erase(i++);} else {++i;}
	}
}


surface_cache surface_rock_cache;


void surface_rock::create(int x, int y, int use_xy, bool fixed_sz_rock_cache) {

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
}
void surface_rock::gen_points(vbo_vnt_block_manager_t &vbo_manager) {
	surface->sd.get_quad_points(vbo_manager.get_pts_vector_for_adding()); // use_tri_strip=0 (quads mode)
	vbo_mgr_ix = vbo_manager.get_offset_for_last_points_added();
}
unsigned surface_rock::get_num_verts() const {return 4*ROCK_NDIV*ROCK_NDIV;} // one quad/4 verts per face

void surface_rock::add_cobjs() {
	coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, WHITE, 0, 0, rock_collision, 1, ROCK_SPHERE_TEX));
}

void surface_rock::draw(float sscale, bool shadow_only, bool reflection_pass,
	vector3d const &xlate, float scale_val, vbo_vnt_block_manager_t &vbo_manager) const
{
	assert(surface);
	if (!is_visible(shadow_only, 0.0, xlate))     return;
	if (reflection_pass && pos.z < water_plane_z) return;
	colorRGBA const color(shadow_only ? WHITE : get_atten_color(WHITE, xlate));
	float const dist(distance_to_camera(pos+xlate));

	if (!shadow_only && 2*get_pt_line_thresh()*radius < dist) { // draw as point
		tree_scenery_pld.add_textured_pt(pos+xlate, color, ROCK_SPHERE_TEX);
		return;
	}
	// Note: rock size scales with distance so that it shrinks to zero area rather than popping in and out;
	// this means that we can't store the final points in the VBO vertex data and must apply a transform here
	color.set_for_cur_shader();
	fgPushMatrix();
	translate_to(pos);
	uniform_scale(scale*get_size_scale(dist, scale_val));
	rotate_into_plus_z(dir);
	assert(vbo_mgr_ix >= 0);
	vbo_manager.render_single(vbo_mgr_ix);
	fgPopMatrix();
}

void surface_rock::destroy() {
	if (surface) {surface->dec_ref();}
	vbo_mgr_ix = -1;
	surface    = NULL;
}


void s_rock::create(int x, int y, int use_xy) {

	UNROLL_3X(scale[i_] = rand_uniform2(0.8, 1.3);)
	gen_spos(x, y, use_xy);
	size   = 0.02*rand_uniform2(0.2, 0.8)/tree_scale;
	if ((rand2()&3) == 0) size *= rand_uniform2(1.2, 8.0);
	dir    = signed_rand_vector2_norm();
	angle  = rand_uniform2(0.0, 360.0);
	radius = size*(scale.x + scale.y + scale.z)/3.0f;
	pos.z += radius*rand_uniform2(-0.1, 0.25);
}

void s_rock::add_cobjs() {
	coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, WHITE, 0, 0, rock_collision, 1, ROCK_SPHERE_TEX));
}

void s_rock::draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val) const {

	if (!is_visible(shadow_only, 1.3*radius, xlate)) return;
	if (reflection_pass && pos.z < water_plane_z)    return;
	colorRGBA const color(shadow_only ? WHITE : get_atten_color(WHITE, xlate));
	float const dist(distance_to_camera(pos+xlate));

	if (!shadow_only && 2*get_pt_line_thresh()*radius < dist) { // draw as point
		tree_scenery_pld.add_textured_pt(pos+xlate, color, ROCK_SPHERE_TEX);
		return;
	}
	color.set_for_cur_shader();
	int const ndiv(max(4, min(N_SPHERE_DIV, (shadow_only ? get_def_smap_ndiv(radius) : int(sscale*radius/dist)))));
	fgPushMatrix();
	translate_to(pos);
	rotate_about(angle, dir);
	scale_by(size*get_size_scale(dist, scale_val)*scale);
	draw_sphere_vbo_raw(ndiv, 1);
	fgPopMatrix();
}



unsigned voxel_rock_manager_t::gen_model_ix(int rseed) {

	models.resize(NUM_VROCK_MODELS);
	unsigned const ix(rseed % models.size());
	if (!models[ix]) {to_gen.insert(ix);} // mark new model for building
	return ix;
}

void voxel_rock_manager_t::build_models(unsigned num_lod_levels) {

	if (to_gen.empty()) return;
	//timer_t timer("Gen Voxel Rocks");
	vector<unsigned> const gen_queue(to_gen.begin(), to_gen.end());

//#pragma omp parallel for schedule(dynamic,1) // not needed, gen_voxel_rock() already uses openmp
	for (int i = 0; i < (int)gen_queue.size(); ++i) { // build scheduled models
		unsigned const ix(gen_queue[i]);
		assert(!models[ix]);
		models[ix].reset(new voxel_model_rock(&ntg, num_lod_levels));
		gen_voxel_rock(*models[ix], all_zeros, 1.0, ROCK_VOX_SZ, 1, ix);
	}
	to_gen.clear();
}

void voxel_rock_manager_t::free_context() {
	for (auto i = models.begin(); i != models.end(); ++i) {if (*i) {(*i)->free_context();}}
	ntg.clear();
}

voxel_rock_manager_t voxel_rock_manager;


void voxel_rock::create(int x, int y, int use_xy) {

	gen_spos(x, y, use_xy);
	radius   = 0.2*rand_uniform2(0.5, 1.0)*rand_float2()/tree_scale;
	rseed    = rand2();
	model_ix = voxel_rock_manager.gen_model_ix(rseed); // will be translated to pos and scaled by radius during rendering
}

void voxel_rock::build_model() {

	float const gen_radius(voxel_rock_manager.get_model(model_ix).get_bsphere().radius);
	assert(gen_radius > 0.0);
	radius /= gen_radius;
}

unsigned voxel_rock::get_tid() const {
	return ((global_voxel_params.tids[0] > 0) ? global_voxel_params.tids[0] : voxel_rock_manager.get_model(model_ix).get_params().tids[0]);
}

void voxel_rock::add_cobjs() {
	coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, WHITE, 0, 0, rock_collision, 1, get_tid()));
}

void voxel_rock::draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val, shader_t &s, bool use_model_texgen) {

	assert(radius > 0.0);
	if (!is_visible(shadow_only, radius, xlate))  return;
	if (reflection_pass && pos.z < water_plane_z) return;
	unsigned const lod_level = 0; // Note: max is VOX_ROCK_NUM_LOD, but makes little difference in runtime
	colorRGBA const color(shadow_only ? WHITE : get_atten_color(WHITE, xlate));
	color.set_for_cur_shader();
	fgPushMatrix();
	translate_to(pos);
	uniform_scale(radius*get_size_scale(distance_to_camera(pos+xlate), scale_val));
	voxel_model_rock &model(voxel_rock_manager.get_model(model_ix));
	if (use_model_texgen) {model.setup_tex_gen_for_rendering(s);} else {select_texture(get_tid());} // set for first rock and reuse?
	model.core_render(s, lod_level, shadow_only, 1); // disable view frustum culling because it's incorrect (due to transform matrices)
	fgPopMatrix();
}


void burnable_scenery_obj::next_frame() {
	if (!animate2 || world_mode != WMODE_GROUND || burn_amt == 1.0) return;
	fire_amt = get_ground_fire_intensity(pos, 2.0*radius);
	if (fire_amt > 0.0) {burn_amt = min(1.0, (burn_amt + 0.003*fticks*fire_amt));}
}
void burnable_scenery_obj::draw_fire(fire_drawer_t &fire_drawer, float rscale, unsigned ix) const {
	if (fire_amt == 0.0 || burn_amt >= 1.0) return; // no fire, or all burned out
	float const fire_radius(rscale*get_bsphere_radius());
	fire_drawer.add_fire((get_center() + vector3d(0.0, 0.0, 0.75*fire_radius)), fire_radius, ((animate2 ? frame_counter : 0) + ix)); // one flame
}


void wood_scenery_obj::calc_type() {type = (char)get_tree_type_from_height(pos.z, global_rand_gen, 1);}

bool wood_scenery_obj::is_from_large_trees() const {
	if (tree_mode == 0) return 0; // no trees enabled: any solution is acceptable, so use the simplest one
	if (tree_mode == 2) return 0; // small trees only: correct/simple
	if (tree_mode == 3) return !is_pine_tree_type(type); // both large and small trees: choose based on pine vs. decid tree type
	return 1;
}
void wood_scenery_obj::cache_closest_tree_type(tree_cont_t const &trees) {
	if (closest_tree_type < 0) {closest_tree_type = trees.get_closest_tree_type(pos);} // compute once and cache
}
int wood_scenery_obj::get_tid() const {
	if (!is_from_large_trees()) {return get_bark_tex_for_tree_type(type);} // small tree
	// else large trees only, or large (non-pine) trees at this height
	return tree_types[closest_tree_type].bark_tex;
}
colorRGBA wood_scenery_obj::get_bark_color(vector3d const &xlate) const {
	colorRGBA const base_color(is_from_large_trees() ? tree_types[closest_tree_type].barkc : get_tree_trunk_color(type, 0));
	return get_atten_color(blend_color(BLACK, base_color, burn_amt, 0), xlate);
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
	calc_type();
	return (type >= 0);
}

void s_log::add_cobjs() {
	coll_id = add_coll_cylinder(pos, pt2, radius, radius2, cobj_params(0.8, get_tree_trunk_color(type, 1), 0, 0, NULL, 0, get_tid()));
}

void s_log::draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val) const {

	if (type < 0) return;
	point const center((pos + pt2)*0.5 + xlate);
	if (!check_visible(shadow_only, get_bsphere_radius(), center)) return;
	if (reflection_pass && 0.5f*(pos.z + pt2.z) < water_plane_z) return;
	colorRGBA const color(shadow_only ? WHITE : get_bark_color(xlate));
	float const dist(distance_to_camera(center));

	if (!shadow_only && get_pt_line_thresh()*(radius + radius2) < dist) { // draw as line
		tree_scenery_pld.add_textured_line(pos+xlate, pt2+xlate, color, get_tid());
		return;
	}
	color.set_for_cur_shader();
	int const ndiv(max(3, min(N_CYL_SIDES, (shadow_only ? get_def_smap_ndiv(2.0*radius) : int(2.0*sscale*radius/dist)))));
	fgPushMatrix();
	translate_to(pos);
	rotate_by_vector(dir);
	if (!shadow_only) {select_texture(TREE_END_TEX);}
	if (dot_product_ptv(dir, get_camera_pos(), center) > 0.0) {draw_circle_normal(0.0, radius,  ndiv, 1, 0.0);}
	else                                                      {draw_circle_normal(0.0, radius2, ndiv, 0, length);}
	if (!shadow_only) {select_texture(get_tid());}
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
	if (pos.z < minz) return 0;
	radius  = rand_uniform2(0.005, 0.01)/tree_scale;
	radius2 = rand_uniform2(0.8*radius, radius);
	pos.z  -= 2.0*radius;
	height  = rand_uniform2(0.01/tree_scale, min(0.05/tree_scale, 4.0*radius)) + 0.015;
		
	if ((rand2()&3) == 0) {
		height  *= rand_uniform2(1.0, 5.0); // larger stump = upright dead tree
		radius  *= 1.5;
		radius2 *= 1.3;
	}
	calc_type();
	return (type >= 0);
}

void s_stump::add_cobjs() {
	coll_id = add_coll_cylinder(pos-point(0.0, 0.0, 0.2*height), pos+point(0.0, 0.0, height),
		radius, radius2, cobj_params(0.8, get_tree_trunk_color(type, 1), 0, 0, NULL, 0, get_tid()));
}

bool s_stump::check_sphere_coll(point &center, float sphere_radius) const {
	if (!dist_less_than(center, pos, (max(height, max(radius, radius2)) + sphere_radius))) return 0; // sphere-sphere coll optimization
	return sphere_vert_cylin_intersect(center, sphere_radius, cylinder_3dw(pos-point(0.0, 0.0, 0.2*height), pos+point(0.0, 0.0, height), radius, radius));
}

void s_stump::draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val) const {

	if (type < 0) return;
	point const center(pos + point(0.0, 0.0, 0.5*height) + xlate);
	if (!check_visible(shadow_only, get_bsphere_radius(), center)) return;
	if (reflection_pass && pos.z < water_plane_z) return;
	colorRGBA const color(shadow_only ? WHITE : get_bark_color(xlate));
	float const dist(distance_to_camera(center));

	if (!shadow_only && get_pt_line_thresh()*(radius + radius2) < dist) { // draw as line
		tree_scenery_pld.add_textured_line((pos+xlate - point(0.0, 0.0, 0.2*height)), (pos+xlate + point(0.0, 0.0, height)), color, get_tid());
		return;
	}
	color.set_for_cur_shader();
	int const ndiv(max(3, min(N_CYL_SIDES, (shadow_only ? get_def_smap_ndiv(2.2*radius) : int(2.2*sscale*radius/dist)))));

	if (get_camera_pos().z > pos.z + height) { // only draw top if visible
		if (!shadow_only) {select_texture(TREE_END_TEX);}
		draw_circle_normal(0.0, radius2, ndiv, 0, pos+vector3d(0.0, 0.0, height));
	}
	if (!shadow_only) {select_texture(get_tid());}
	draw_fast_cylinder(pos-vector3d(0.0, 0.0, 0.2*height), pos+vector3d(0.0, 0.0, height), radius, radius2, ndiv, 1);
}


// returns 0 for no plants, 1 for land plants, and 2 for water plants
int plant_base::create(int x, int y, int use_xy, float minz) {

	gen_spos(x, y, use_xy);

	if (pos.z < minz) {
		if (pos.z + (0.4/tree_scale + 0.025) > get_min_water_plane_z()) return 0; // max plant height above min water plane zval
		if (NUM_WATER_PLANT_TYPES == 0) return 0; // no water plant types
		return 2; // water plant
	}
	if (get_rel_height(pos.z, -zmax_est, zmax_est) > 0.62) return 0; // altitude too high for plants
	if (NUM_LAND_PLANT_TYPES == 0) return 0; // no land plant types
	return 1; // land plant
}

void plant_base::next_frame() {
	burnable_scenery_obj::next_frame();
	if (burn_amt == 1.0) {remove_cobjs();} // remove leaves, but keep stem for plant
}
colorRGBA plant_base::get_plant_color(vector3d const &xlate) const {
	return get_atten_color(blend_color(BLACK, WHITE, burn_amt, 0), xlate);
}


int s_plant::create(int x, int y, int use_xy, float minz) {

	vbo_mgr_ix = -1;
	int const ret(plant_base::create(x, y, use_xy, minz));
	if (ret == 0) return 0;
	if (ret == 2) {type = NUM_LAND_PLANT_TYPES + (rand2() % NUM_WATER_PLANT_TYPES);} // water plant
	else          {type = rand2() % NUM_LAND_PLANT_TYPES;} // land plant
	radius = rand_uniform2(0.0025, 0.0045)/tree_scale;
	height = rand_uniform2(0.2, 0.4)/tree_scale + 0.025;
	return 1;
}

void s_plant::create2(point const &pos_, float height_, float radius_, int type_, int calc_z) {
	create_no_verts(pos_, height_, radius_, type_, calc_z);
}

void s_plant::create_no_verts(point const &pos_, float height_, float radius_, int type_, int calc_z, bool land_plants_only) {

	vbo_mgr_ix = -1;
	type   = abs(type_) % (land_plants_only ? (int)NUM_LAND_PLANT_TYPES : (int)NUM_PLANT_TYPES);
	pos    = pos_;
	radius = radius_;
	height = height_;
	if (calc_z) {pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);}
}

colorRGBA const &s_plant::get_leaf_color() const {return pltype[type].leafc;}
colorRGBA const &s_plant::get_stem_color() const {return pltype[type].stemc;}
int s_plant::get_leaf_tid() const {return pltype[type].tid;}

void s_plant::add_cobjs() {

	point cpos(pos), cpos2(pos), bpos(pos);
	float const wscale(radius*tree_scale/0.004f), r2(radius+0.07f*wscale*(height+0.03f));
	cpos.z  += height;
	cpos2.z += 3.0*height/(36.0*height + 4.0);
	bpos.z  -= 0.1*height;
	coll_id  = add_coll_cylinder(cpos2, cpos, r2, radius, cobj_params(0.4, get_leaf_color(), 0, 0, NULL, 0, get_leaf_tid())); // leaves
	coll_id2 = add_coll_cylinder(bpos, cpos, radius, 0.0, cobj_params(0.4, get_stem_color(), 0, 0, NULL, 0, WOOD_TEX      )); // trunk
}

bool s_plant::check_sphere_coll(point &center, float sphere_radius) const { // used in tiled terrain mode
	if (!dist_less_than(center, pos, (0.5f*(height + radius) + sphere_radius))) return 0; // sphere-sphere coll optimization
	return sphere_vert_cylin_intersect(center, sphere_radius, cylinder_3dw(pos-point(0.0, 0.0, 0.1*height), get_top_pt(), radius, radius));
}

void s_plant::create_leaf_points(vector<vert_norm_comp> &points, float plant_scale, float nlevels_scale, unsigned nrings) const {
	
	// Note: could scale leaves for different plant types differently in x vs. y to allow for non-square textures (tighter bounds = lower fillrate)
	float const wscale(250.0*radius*plant_scale), theta0((int(1.0E6*height)%360)*TO_RADIANS);
	float const llen(wscale*pltype[type].leaf_length), blwidth(wscale*pltype[type].leaf_width_base), elwidth(wscale*pltype[type].leaf_width_end);
	unsigned const nlevels(max(1U, unsigned(36.0*height*plant_scale*nlevels_scale)));
	float rdeg(30.0);
	points.resize(4*nlevels*nrings);

	for (unsigned j = 0, ix = 0; j < nlevels; ++j) { // could do the same optimizations as the high detail pine tree
		float const sz(0.07*(height + 0.03/plant_scale)*((nlevels - j + 3.0)/(float)nlevels));
		float const z((j + 3.0)*height/(nlevels + 4.0));

		for (unsigned k = 0; k < nrings; ++k) {
			float const theta(TWO_PI*(3.3*j + 0.2*k) + theta0);
			int const val(int(((int(1.0E6f*height))*(5463*j + 537879*k))%301));
			rdeg += 0.01f*(val - 150);
			add_rotated_quad_pts(&points.front(), ix, theta, z, pos, sz*blwidth, sz*elwidth, sz*llen, sz*rdeg/45.0f);
		}
	}
}

void s_plant::gen_points(vbo_vnc_block_manager_t &vbo_manager, vector<vert_norm_comp> &pts) {

	if (vbo_mgr_ix >= 0) return; // already generated
	create_leaf_points(pts, tree_scale);
	assert(!pts.empty());
	vbo_mgr_ix = vbo_manager.add_points_with_offset(pts, get_leaf_color());

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
		} // for i
	}
}

// to be called when the plant is translated or zval changes
void s_plant::update_points_vbo(vbo_vnc_block_manager_t &vbo_manager) {

	assert(vbo_mgr_ix >= 0);
	static vector<vert_norm_comp> pts; // reused across calls
	create_leaf_points(pts, tree_scale);
	vbo_manager.update_range(pts, get_leaf_color(), vbo_mgr_ix, vbo_mgr_ix+1);
}

bool s_plant::update_zvals(int x1, int y1, int x2, int y2, vbo_vnc_block_manager_t &vbo_manager) {

	float orig_z(pos.z);
	if (!scenery_obj::update_zvals(x1, y1, x2, y2)) return 0;
	for (vector<vert_wrap_t>::iterator i = berries.begin(); i != berries.end(); ++i) {i->v.z += (pos.z - orig_z);}
	if (vbo_mgr_ix >= 0) {update_points_vbo(vbo_manager);} // leaf pos changed, regenerate leaf points and re-upload VBO sub-data
	return 1;
}

bool s_plant::is_shadowed() const {

	if (world_mode != WMODE_GROUND) return 0;
	int const light(get_light());

	for (unsigned i = 0; i < 3; ++i) {
		point p(pos);
		p.z += 0.5*i*height;
		if (is_visible_to_light_cobj(p, light, (radius + height), coll_id2, 0)) return 0;
	}
	return 1;
}

bool s_plant::is_water_plant() const {return (type >= (int)NUM_LAND_PLANT_TYPES);}

void s_plant::draw_stem(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate) const {

	if (world_mode == WMODE_INF_TERRAIN && is_water_plant() && (reflection_pass || (!shadow_only && pos.z < water_plane_z && get_camera_pos().z > water_plane_z))) return; // underwater, skip
	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	if (!check_visible(shadow_only, (height + radius), pos2)) return;
	bool const shadowed(shadow_only ? 0 : is_shadowed());
	colorRGBA color(get_stem_color()*(shadowed ? SHADOW_VAL : 1.0));
	if (is_water_plant() && !shadow_only) {water_color_atten_at_pos(color, pos+xlate);}
	float const dist(distance_to_camera(pos2));

	if (!shadow_only && 2*get_pt_line_thresh()*radius < dist) { // draw as line
		tree_scenery_pld.add_textured_line((pos+xlate - point(0.0, 0.0, 0.1*height)), (pos+xlate + point(0.0, 0.0, 0.75*height)), color, WOOD_TEX);
	}
	else {
		int const ndiv(max(3, min(N_CYL_SIDES, (shadow_only ? get_def_smap_ndiv(2.0*radius) : int(2.0*sscale*radius/dist)))));
		if (!shadow_only) {color.set_for_cur_shader();}
		draw_fast_cylinder((pos - point(0.0, 0.0, 0.1*height)), get_top_pt(), radius, 0.0, ndiv, !shadow_only, 0, 0, nullptr, 6.0);
	}
}

void s_plant::shader_state_t::set_color_scale(shader_t &s, colorRGBA const &color) {
	s.ensure_uniform_loc(color_scale_loc, "color_scale");
	s.set_uniform_color(color_scale_loc, color);
}
void s_plant::shader_state_t::set_normal_scale(shader_t &s, float normal_scale) {
	s.ensure_uniform_loc(normal_scale_loc, "normal_scale");
	s.set_uniform_float(normal_scale_loc, normal_scale);
}
void s_plant::shader_state_t::set_wind_scale(shader_t &s, float wscale) {
	if (wscale == wind_scale) return;
	s.ensure_uniform_loc(wind_scale_loc, "wind_scale");
	s.set_uniform_float(wind_scale_loc, wscale);
	wind_scale = wscale;
}
void s_plant::shader_state_t::set_wind_add(shader_t &s, float w_add) {
	s.ensure_uniform_loc(wind_add_loc, "wind_add");
	s.set_uniform_float(wind_add_loc, w_add);
}

void s_plant::draw_leaves(shader_t &s, vbo_vnc_block_manager_t &vbo_manager, bool shadow_only, bool reflection_pass, vector3d const &xlate, shader_state_t &state) const {

	if (burn_amt == 1.0) return;
	if (world_mode == WMODE_INF_TERRAIN && is_water_plant() && (reflection_pass || (!shadow_only && pos.z < water_plane_z && get_camera_pos().z > water_plane_z))) return; // underwater, skip
	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	if (!check_visible(shadow_only, 0.5f*(height + radius), pos2)) return;
	bool const shadowed((shadow_only || (ENABLE_PLANT_SHADOWS && shadow_map_enabled())) ? 0 : is_shadowed());
	float const wind_scale(berries.empty() ? (is_water_plant() ? 5.0 : 1.0) : 0.0); // no wind if this plant type has berries
	bool const set_color(!shadow_only && (is_water_plant() || burn_amt > 0.0));
	if (set_color) {state.set_color_scale(s, get_plant_color(xlate));}
	if (shadowed ) {state.set_normal_scale(s, 0.0);}
	state.set_wind_scale(s, wind_scale);
	select_texture((draw_model == 0) ? get_leaf_tid() : WHITE_TEX); // could pre-bind textures and select using shader int, but probably won't improve performance
	assert(vbo_mgr_ix >= 0);
	vbo_manager.render_single(vbo_mgr_ix);
	if (set_color) {state.set_color_scale(s, WHITE);}
	if (shadowed ) {state.set_normal_scale(s, 1.0);}
}

void s_plant::draw_berries(shader_t &s, vector3d const &xlate) const { // drawn using instancing

	if (berries.empty() || burn_amt > 0.5) return;
	point const pos2(pos + xlate + point(0.0, 0.0, 0.5*height));
	if (!sphere_in_camera_view(pos2, 0.5f*(height + radius), 2)) return;
	float const size_scale(get_pt_line_thresh()*radius/distance_to_camera(pos2));
	if (size_scale < 1.2) return; // too small/far away
	int const ndiv(max(4, min(16, int(2.0*size_scale))));
	s.set_cur_color(pltype[type].berryc);
	fgPushMatrix();
	uniform_scale(0.25*radius);
	int const vo_loc(s.get_attrib_loc("vertex_offset", 0)), vos_loc(s.get_uniform_loc("vertex_offset_scale"));
	s.set_uniform_float(vos_loc, 1.0/(0.25*radius)); // enable
	shader_float_matrix_uploader<3,1>::enable(vo_loc, 1, (float const *)get_dynamic_vbo_ptr(&berries.front().v.x, berries.size()*sizeof(point)));
	draw_sphere_vbo_raw(ndiv, 0, 0, berries.size()); // untextured
	shader_float_matrix_uploader<3,1>::disable(vo_loc);
	s.set_uniform_float(vos_loc, 0.0); // disable
	fgPopMatrix();
}

void s_plant::remove_cobjs() {
	remove_reset_coll_obj(coll_id2);
	scenery_obj::remove_cobjs();
}


int leafy_plant::create(int x, int y, int use_xy, float minz, unsigned plant_ix_) {
	
	vbo_mgr_ix = -1;
	plant_ix   = plant_ix_;
	int const ret(plant_base::create(x, y, use_xy, minz));
	if (ret == 0) return 0;
	if (ret == 2) {type = LEAFY_PLANT_UW;} // underwater
	else {
		float const relh(get_rel_height(pos.z, -zmax_est, zmax_est));
		if      (relh < 0.46) {type = LEAFY_PLANT_DIRT;} // dirt/sand
		else if (relh < 0.60) {type = LEAFY_PLANT_GRASS;} // grass
		else if (relh < 0.75) {type = LEAFY_PLANT_ROCK;} // rock
		else return 0; // snow
	}
	radius = rand_uniform2(0.06, 0.12)/tree_scale;
	gen_leaves();
	return 1;
}

void leafy_plant::create2(point const &pos_, float radius_, int type_, int calc_z, unsigned plant_ix_) {

	vbo_mgr_ix = -1;
	plant_ix   = plant_ix_;
	type   = abs(type_)%NUM_LEAFY_PLANT_TYPES;
	pos    = pos_;
	radius = radius_;
	if (calc_z) {pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);}
	gen_leaves();
}

void leafy_plant::gen_leaves() {

	rand_gen_t rgen;
	rgen.set_state(rand2(), 123);
	leaves.resize(rgen.rand_uniform_uint(4, 8));
	float const delta_angle(TWO_PI/leaves.size());

	for (auto i = leaves.begin(); i != leaves.end(); ++i) { // for each leaf
		float const dxy(rgen.rand_uniform(0.7, 1.3)), dz(rgen.rand_uniform(-0.1, 0.4));
		float const angle(delta_angle*((i - leaves.begin()) + 0.5*rgen.rand_float()));
		float const rscale(rgen.rand_uniform(0.5, 1.0));
		vector3d const delta((1.2*rscale*radius)*point(-dxy*cos(angle), -dxy*sin(angle), dz));
		i->m = glm::translate(i->m, vec3_from_vector3d(pos + delta));
		i->m = glm::scale(i->m, vec3_from_vector3d(vector3d(1.0, 1.0, 0.75)*(rscale*radius)));
		i->m = glm::rotate(i->m, TO_RADIANS*135.0f, glm::vec3(0, 1, 0));
		i->m = glm::rotate(i->m, angle, glm::vec3(-1, 0, -1));
	}
}

void leafy_plant::gen_points(vbo_vnt_block_manager_t &vbo_manager, vector<vert_norm_tc> const &sphere_verts) {

	vector<vert_norm_tc> &verts(vbo_manager.get_pts_vector_for_adding());

	for (auto i = leaves.begin(); i != leaves.end(); ++i) { // for each leaf
		glm::mat3 const nm(glm::inverseTranspose(glm::mat3(i->m)));

		for (auto v = sphere_verts.begin(); v != sphere_verts.end(); ++v) {
			verts.emplace_back(v->v, vector3d_from_vec3(nm * vec3_from_vector3d(v->n)), v->t); // Note: n is normalized in the shader
			i->m.apply_to_vector3d(verts.back().v);
		}
	}
	vbo_mgr_ix = vbo_manager.get_offset_for_last_points_added();
}

void leafy_plant::add_cobjs() {
	coll_id = add_coll_sphere(pos, radius, cobj_params(0.5, WHITE, 0, 0, leafy_plant_collision, plant_ix, get_tid()));
}

void leafy_plant::next_frame() {

	plant_base::next_frame();
	float const delta_energy(abs(cur_motion_energy - prev_motion_energy));
	if (delta_energy > 1.0) {motion_amt = min(1.0, (motion_amt + 0.05*delta_energy));} // increase wind effect to make leaves move
	prev_motion_energy = cur_motion_energy;
	cur_motion_energy  = 0.0;
	if (motion_amt == 0.0) return;
	motion_amt *= pow(2.0, -0.3*fticks); // exponential decay
	if (motion_amt < 0.01) {motion_amt = 0.0;} // done
}

bool leafy_plant::update_zvals(int x1, int y1, int x2, int y2, vbo_vnt_block_manager_t &vbo_manager) {
	
	if (!scenery_obj::update_zvals(x1, y1, x2, y2)) return 0;
	//if (vbo_mgr_ix >= 0) {update_points_vbo(vbo_manager);}
	delta_z += dz;
	return 1;
}

int leafy_plant::get_tid() const {
	unsigned const tids[NUM_LEAFY_PLANT_TYPES] = {LEAF2_TEX, PLANT3_TEX, LEAF_TEX, PAPAYA_TEX}; // LEAF3_TEX is okay but has artifacts at a distance; PALM_FROND_TEX needs clipping
	return tids[type];
}

void leafy_plant::draw_leaves(shader_t &s, bool shadow_only, bool reflection_pass, vector3d const &xlate, s_plant::shader_state_t &state, vbo_vnt_block_manager_t &vbo_manager) const {
	
	if (burn_amt == 1.0) return;
	if (!is_visible(shadow_only, radius, xlate))  return;
	if (reflection_pass && pos.z < water_plane_z) return;
	(shadow_only ? WHITE : get_plant_color(xlate)).set_for_cur_shader(); // no underwater case yet
	bool const is_underwater(pos.z < water_plane_z);
	select_texture(get_tid());
	assert(vbo_mgr_ix >= 0);
	if (delta_z != 0.0) {fgPushMatrix(); fgTranslate(0, 0, delta_z);} // not the cleanest or most efficient solution, but much simpler than updating the VBO data
	if (motion_amt > 0.0) {state.set_wind_add(s, 0.005*motion_amt);}
	if (is_underwater) {state.set_wind_scale(s, 5.0*LEAFY_PLANT_WIND);}
	vbo_manager.render_single(vbo_mgr_ix);
	if (motion_amt > 0.0) {state.set_wind_add(s, 0.0);} // restore orig value
	if (is_underwater) {state.set_wind_scale(s, LEAFY_PLANT_WIND);}
	if (delta_z != 0.0) {fgPopMatrix();}
}


// ************ SCENERY OBJECT INTERFACE/WRAPPERS/DRIVERS ************


template<typename T> void draw_scenery_vector(vector<T> &v, float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val) {
	for (unsigned i = 0; i < v.size(); ++i) {v[i].draw(sscale, shadow_only, reflection_pass, xlate, scale_val);}
}

template<typename T> void draw_scenery_vector_fires(vector<T> const &v, fire_drawer_t &fire_drawer, float rscale) {
	for (unsigned i = 0; i < v.size(); ++i) {v[i].draw_fire(fire_drawer, rscale, i);}
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

template<typename T> void free_scenery_vector_cobjs(vector<T> &v) {
	for (unsigned i = 0; i < v.size(); ++i) {v[i].remove_cobjs();}
}

template<typename T> void update_bcube(vector<T> &v, cube_t &bcube) {
	for (unsigned i = 0; i < v.size(); ++i) {v[i].add_bounds_to_bcube(bcube);}
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
template<typename T, typename ARG> void update_scenery_zvals_vector(vector<T> &v, int x1, int y1, int x2, int y2, bool &updated, ARG &arg) {

	for (unsigned i = 0; i < v.size(); ++i) { // zval has change, remove and re-add cobjs
		if (v[i].update_zvals(x1, y1, x2, y2, arg)) {
			v[i].remove_cobjs();
			v[i].add_cobjs();
			updated = 1;
		}
	}
}


void scenery_group::clear_vbos() {
	
	for (unsigned i = 0; i < rock_shapes.size(); ++i) {rock_shapes[i].clear_vbo();}
	plant_vbo_manager.clear_vbo();
	rock_vbo_manager .clear_vbo();
	leafy_vbo_manager.clear_vbo();
}

void scenery_group::clear() {

	for (auto i = surface_rocks.begin(); i != surface_rocks.end(); ++i) {i->destroy();}
	clear_vbos();
	free_cobjs();
	rock_shapes.clear();
	surface_rocks.clear();
	voxel_rocks.clear();
	rocks.clear();
	logs.clear();
	stumps.clear();
	plants.clear();
	leafy_plants.clear();
	plant_vbo_manager.clear();
	rock_vbo_manager .clear();
	leafy_vbo_manager.clear();
	generated = 0;
}

void scenery_group::free_cobjs() {

	free_scenery_vector_cobjs(rock_shapes);
	free_scenery_vector_cobjs(surface_rocks);
	free_scenery_vector_cobjs(voxel_rocks);
	free_scenery_vector_cobjs(rocks);
	free_scenery_vector_cobjs(logs);
	free_scenery_vector_cobjs(stumps);
	free_scenery_vector_cobjs(plants);
	free_scenery_vector_cobjs(leafy_plants);
}

void scenery_group::add_cobjs() {

	add_scenery_vector_cobjs(rock_shapes);
	add_scenery_vector_cobjs(surface_rocks);
	add_scenery_vector_cobjs(voxel_rocks);
	add_scenery_vector_cobjs(rocks);
	add_scenery_vector_cobjs(logs);
	add_scenery_vector_cobjs(stumps);
	add_scenery_vector_cobjs(plants);
	add_scenery_vector_cobjs(leafy_plants);
	calc_bcube(); // okay to put here, since add_cobjs() is always called
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
	coll |= check_scenery_vector_sphere_coll(leafy_plants,  center, radius);
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
	shift_scenery_vector(leafy_plants,  vd);
}

void scenery_group::calc_bcube() {

	all_bcube.set_to_zeros();
	update_bcube(rock_shapes,   all_bcube);
	update_bcube(surface_rocks, all_bcube);
	update_bcube(voxel_rocks,   all_bcube);
	update_bcube(rocks,         all_bcube);
	update_bcube(logs,          all_bcube);
	update_bcube(stumps,        all_bcube);
	update_bcube(plants,        all_bcube);
	update_bcube(leafy_plants,  all_bcube);
}

// update region is inclusive: [x1,x2]x[y1,y2]
bool scenery_group::update_zvals(int x1, int y1, int x2, int y2) { // inefficient, should use spatial subdivision

	assert(x1 <= x2 && y1 <= y2);
	bool updated(0);
	update_scenery_zvals_vector(rock_shapes,   x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(voxel_rocks,   x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(rocks,         x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(logs,          x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(stumps,        x1, y1, x2, y2, updated);
	update_scenery_zvals_vector(plants,        x1, y1, x2, y2, updated, plant_vbo_manager);
	update_scenery_zvals_vector(leafy_plants,  x1, y1, x2, y2, updated, leafy_vbo_manager);
	update_scenery_zvals_vector(surface_rocks, x1, y1, x2, y2, updated); // rock_vbo_manager not updated (pos is applied to the MVM), so not passed into the call
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
	plants.back().create2(pos, height, radius, type, calc_z);
}
void scenery_group::add_leafy_plant(point const &pos, float radius, int type, int calc_z) {
	assert(radius > 0.0);
	unsigned const plant_ix(leafy_plants.size());
	leafy_plants.push_back(leafy_plant());
	leafy_plants.back().create2(pos, radius, type, calc_z, plant_ix);
}

bool check_valid_scenery_pos(scenery_obj const &obj) {return check_valid_scenery_pos(obj.get_pos(), obj.get_radius());}

void scenery_group::gen(int x1, int y1, int x2, int y2, float vegetation_, bool fixed_sz_rock_cache, tree_cont_t const &trees) { // called in tiled terrain mode

	//RESET_TIME;
	unsigned const smod(max(200U, unsigned(3.321f*XY_MULT_SIZE/(tree_scale+1))));
	float const min_stump_z(water_plane_z + 0.010*zmax_est);
	float const min_plant_z(water_plane_z + 0.016*zmax_est);
	float const min_log_z  (water_plane_z - 0.040*zmax_est);
	generated = 1;

	for (int i = y1; i < y2; ++i) {
		for (int j = x1; j < x2; ++j) {
			global_rand_gen.rseed1 = 786433* (i + yoff2) + 196613 *rand_gen_index;
			global_rand_gen.rseed2 = 6291469*(j + xoff2) + 1572869*rand_gen_index;
			int const val(rand2_seed_mix()%smod);
			if (val >= 150) continue;
			rand2_mix();
			bool const veg((global_rand_gen.rseed1&127)/128.0 < vegetation_);
			
			if (val >= 100) { // +50% leafy plants
				leafy_plant plant;
				if (veg && plant.create(j, i, 1, min_plant_z, leafy_plants.size())) {
					if (!check_valid_scenery_pos(plant.get_pos(), plant.get_radius())) continue;
					leafy_plants.push_back(plant);
				}
			}
			else if (veg && rand2()%100 < 35) { // Note: numbers below were based on 30% plants but we now have 35% plants
				s_plant plant; // 35%
				if (plant.create(j, i, 1, min_plant_z)) {
					if (!check_valid_scenery_pos(plant)) continue;
					plants.push_back(plant);
					plant.add_bounds_to_bcube(all_bcube);
				}
			}
			else if (val < 5) { // 3.5%
				rock_shapes.push_back(rock_shape3d());
				rock_shapes.back().create(j, i, 1);
				if (!check_valid_scenery_pos(rock_shapes.back())) {rock_shapes.pop_back(); continue;}
				rock_shapes.back().add_bounds_to_bcube(all_bcube);
			}
			else if (val < 15) { // 7%
				surface_rocks.push_back(surface_rock());
				surface_rocks.back().create(j, i, 1, fixed_sz_rock_cache);
				if (!check_valid_scenery_pos(surface_rocks.back())) {surface_rocks.pop_back(); continue;}
				surface_rocks.back().add_bounds_to_bcube(all_bcube);
			}
			else if ((use_voxel_rocks == 1 || (use_voxel_rocks >= 2 && vegetation == 0.0)) && val < 35) { // 0=never, 1=always, 2=only when no vegetation
				voxel_rocks.push_back(voxel_rock());
				voxel_rocks.back().create(j, i, 1);
				if (!check_valid_scenery_pos(voxel_rocks.back())) {voxel_rocks.pop_back(); continue;}
			}
			else if (val < 50) { // 24.5%
				rocks.push_back(s_rock());
				rocks.back().create(j, i, 1);
				if (!check_valid_scenery_pos(rocks.back())) {rocks.pop_back(); continue;}
				rocks.back().add_bounds_to_bcube(all_bcube);
			}
			else if (veg && val < 85) { // 24.5%
				s_log log;
				if (log.create(j, i, 1, min_log_z)) {
					if (!check_valid_scenery_pos(log)) continue;
					logs.push_back(log);
					log.add_bounds_to_bcube(all_bcube);
				}
			}
			else if (veg) { // 10.5%
				s_stump stump;
				if (stump.create(j, i, 1, min_stump_z)) {
					if (!check_valid_scenery_pos(stump)) continue;
					stumps.push_back(stump);
					stump.add_bounds_to_bcube(all_bcube);
				}
			}
		} // for j
	} // for i
	if (!fixed_sz_rock_cache) {surface_rock_cache.clear_unref();}
	post_gen_setup(trees);
}

void scenery_group::post_gen_setup(tree_cont_t const &trees) {

	if (!leafy_plants.empty()) {
		bool const use_tri_strip = 1;
		vector<vert_norm_tc> sphere_verts;
		leafy_vbo_manager.clear();
		add_sphere_quads(sphere_verts, nullptr, all_zeros, 1.0, 16, use_tri_strip,  0.5, 1.0, 0.125, 1.0); // only emit the textured top part of the sphere + the 'stem'
		if (use_tri_strip) {leafy_vbo_manager.set_prim_type(GL_TRIANGLE_STRIP);}
		sort(plants.begin(), plants.end()); // sort by type, before creating VBO data
		unsigned num_lp_leaves(0);
		for (auto const &p : leafy_plants) {num_lp_leaves += p.num_leaves();}
		leafy_vbo_manager.reserve_pts(num_lp_leaves*sphere_verts.size());
		for (auto &p : leafy_plants) {p.gen_points(leafy_vbo_manager, sphere_verts);}
	}
	if (!surface_rocks.empty()) {
		unsigned tot_num_verts(0);
		for (auto const &r : surface_rocks) {tot_num_verts += r.get_num_verts();}
		rock_vbo_manager.reserve_pts(tot_num_verts);
		for (auto &r : surface_rocks) {r.gen_points(rock_vbo_manager);}
	}
	sort(plants.begin(), plants.end()); // sort by type, before creating VBO data
	for (s_plant &plant : plants) {plant.gen_points(plant_vbo_manager, temp_pts);}
	if (!voxel_rocks.empty()) {voxel_rock_manager.build_models(VOX_ROCK_NUM_LOD);}
		
	for (auto &r : voxel_rocks) {
		r.build_model();
		r.add_bounds_to_bcube(all_bcube);
	}
	for (auto &log   : logs  ) {log  .cache_closest_tree_type(trees);}
	for (auto &stump : stumps) {stump.cache_closest_tree_type(trees);}
	//PRINT_TIME("Gen Scenery");
}

void scenery_group::draw_plant_leaves(shader_t &s, bool shadow_only, vector3d const &xlate, bool reflection_pass) {

	if (plants.empty() && leafy_plants.empty()) return; // nothing to draw
	bool const do_update(!shadow_only && !reflection_pass && world_mode == WMODE_GROUND);
	s.set_specular(0.25, 20.0); // a small amount of specular
	s.add_uniform_float("wind_zscale", 2.0);
	s_plant::shader_state_t state;

	if (!plants.empty()) {
		plant_vbo_manager.upload();
		plant_vbo_manager.begin_render();

		for (unsigned i = 0; i < plants.size(); ++i) {
			if (do_update) {plants[i].next_frame();}
			plants[i].draw_leaves(s, plant_vbo_manager, shadow_only, reflection_pass, xlate, state);
		}
		plant_vbo_manager.end_render();
	}
	if (!leafy_plants.empty()) {
		state.set_wind_scale(s, LEAFY_PLANT_WIND);
		s.add_uniform_float("tex_coord_weight", 2.0); // using tex coords, not texgen from vert ID
		leafy_vbo_manager.upload();
		leafy_vbo_manager.begin_render();

		for (unsigned i = 0; i < leafy_plants.size(); ++i) {
			if (do_update) {leafy_plants[i].next_frame();}
			leafy_plants[i].draw_leaves(s, shadow_only, reflection_pass, xlate, state, leafy_vbo_manager);
		}
		leafy_vbo_manager.end_render();
		s.add_uniform_float("tex_coord_weight", 0.0); // reset
	}
	s.clear_specular();
}


void scenery_group::draw_opaque_objects(shader_t &s, shader_t &vrs, bool shadow_only, vector3d const &xlate, bool draw_pld, float scale_val, bool reflection_pass) {

	if (!shadow_only && !reflection_pass && world_mode == WMODE_GROUND) { // apply burn
		for (auto i = logs.begin()  ; i != logs.end()  ; ++i) {i->next_frame();}
		for (auto i = stumps.begin(); i != stumps.end(); ++i) {i->next_frame();}
	}
	select_texture(DARK_ROCK_TEX);
	for (unsigned i = 0; i < rock_shapes.size(); ++i) {rock_shapes[i].draw(shadow_only, reflection_pass, xlate);}
	int const sscale(int((do_zoom ? ZOOM_FACTOR : 1.0)*window_width));
	rock_vbo_manager.upload();
	rock_vbo_manager.begin_render();
	if (!shadow_only) {select_texture(ROCK_SPHERE_TEX);}

	for (unsigned i = 0; i < surface_rocks.size(); ++i) {
		surface_rocks[i].draw(sscale, shadow_only, reflection_pass, xlate, scale_val, rock_vbo_manager);
	}
	rock_vbo_manager.end_render();

	if (!rocks.empty()) {
		begin_sphere_draw(1); // textured=1
		draw_scenery_vector(rocks, sscale, shadow_only, reflection_pass, xlate, scale_val);
		end_sphere_draw();
	}
	draw_scenery_vector(logs,   sscale, shadow_only, reflection_pass, xlate, scale_val);
	draw_scenery_vector(stumps, sscale, shadow_only, reflection_pass, xlate, scale_val);
	if (!shadow_only) {select_texture(WOOD_TEX);} // plant stems use wood texture
	for (unsigned i = 0; i < plants.size(); ++i) {plants[i].draw_stem(sscale, shadow_only, reflection_pass, xlate);}

	if (!shadow_only && !plants.empty()) { // no berry shadows
		select_texture(WHITE_TEX); // berries are untextured
		s.set_specular(0.9, 80.0);
		begin_sphere_draw(0); // textured=0
		for (unsigned i = 0; i < plants.size(); ++i) {plants[i].draw_berries(s, xlate);}
		end_sphere_draw();
		s.clear_specular();
	}
	if (draw_pld) {tree_scenery_pld.draw_and_clear();}

	if (!voxel_rocks.empty()) {
		bool const need_restore(setup_voxel_rocks_shader(vrs, shadow_only)); // uses a custom shader
		
		for (unsigned i = 0; i < voxel_rocks.size(); ++i) {
			voxel_rocks[i].draw(sscale, shadow_only, reflection_pass, xlate, scale_val, (shadow_only ? s : vrs), 0); // Note: no model texgen
		}
		if (need_restore) {s.make_current();} // restore original shader
	}
}

bool scenery_group::setup_voxel_rocks_shader(shader_t &vrs, bool shadow_only) const {

	if (voxel_rocks.empty() || shadow_only)  return 0; // setup not needed
	if (vrs.is_setup()) {vrs.make_current(); return 1;} // already setup
	bool const v(world_mode == WMODE_GROUND), use_noise_tex(0), use_bmap(0), use_smap(v); // FIXME: no TT shadow maps, fog not setup (not needed?)
	setup_procedural_shaders(vrs, 0.0, v, v, use_smap, use_bmap, use_noise_tex,
		global_voxel_params.top_tex_used, global_voxel_params.tex_scale, global_voxel_params.noise_scale, global_voxel_params.tex_mix_saturate);
	vrs.set_cur_color(WHITE);
	return 1;
}


void scenery_group::draw(bool shadow_only, vector3d const &xlate) {

	if (all_bcube.is_zero_area() || !camera_pdu.cube_visible(all_bcube)) return; // empty, or no scenery is visible
	// draw stems, rocks, logs, and stumps
	shader_t s, vrs;

	if (!shadow_only) {
		setup_smoke_shaders(s, 0.0, 0, 1, 0, 1, 1, 0, 0, 2, 0, 0, 1, 0, 0.0, 0.0, 0, 0, 1); // direct lighting + dlights + shadow map, is_outside=1
	}
	else {
		s.begin_simple_textured_shader(0.0, !shadow_only); // with lighting, unless shadow_only
	}
	draw_opaque_objects(s, vrs, shadow_only, xlate, 1);
	s.end_shader();

	if (!(plants.empty() && leafy_plants.empty())) { // draw leaves
		set_leaf_shader(s, 0.9, 0, 0, shadow_only, get_plant_leaf_wind_mag(shadow_only), underwater, ENABLE_PLANT_SHADOWS, (ENABLE_PLANT_SHADOWS && !shadow_only), 1, shadow_only);
		draw_plant_leaves(s, shadow_only, xlate);
		s.end_shader();
	}
}

void scenery_group::draw_fires(shader_t &s) const {
	
	fire_drawer_t fire_drawer;
	draw_scenery_vector_fires(logs,   fire_drawer, 1.8);
	draw_scenery_vector_fires(stumps, fire_drawer, 1.8);
	draw_scenery_vector_fires(plants, fire_drawer, 1.5);
	draw_scenery_vector_fires(leafy_plants, fire_drawer, 1.8);
	fire_drawer.draw(s);
}

void scenery_group::leafy_plant_coll(unsigned plant_ix, float energy) {
	assert(plant_ix < leafy_plants.size());
	leafy_plants[plant_ix].obj_collision(energy);
}

bool scenery_group::choose_butterfly_dest(point &dest, sphere_t &plant_bsphere, rand_gen_t &rgen) const {
	unsigned const tot_plants(plants.size() + leafy_plants.size());
	if (tot_plants == 0) return 0; // no plants
	unsigned const plant_ix(rgen.rand() % tot_plants);

	if (plant_ix < plants.size()) { // choose a plant
		s_plant const &plant(plants[plant_ix]);
		if (plant.is_water_plant() || plant.get_pos().z < water_plane_z) return 0; // butterfly can't land on an underwater plant, try again next time
		dest = plant.get_top_pt();
		plant_bsphere = sphere_t(plant.get_pos(), plant.get_bsphere_radius());
	}
	else { // choose a leafy plant
		leafy_plant const &plant(leafy_plants[plant_ix - plants.size()]);
		if (plant.get_pos().z < water_plane_z) return 0; // butterfly can't land on an underwater plant, try again next time
		dest = plant.get_top_pt();
		plant_bsphere = sphere_t(plant.get_pos(), plant.get_bsphere_radius());
	}
	return 1;
}


scenery_group all_scenery;


void gen_scenery(tree_cont_t const &trees) { // called in ground mode

	if (has_scenery2) {all_scenery.post_gen_setup(trees); return;} // don't generate scenery if some has already been added
	all_scenery.clear();
	all_scenery = scenery_group(); // really force a clear
	has_scenery = 0;
	if (DISABLE_SCENERY) return;
	has_scenery = 1;
	all_scenery.gen(1, 1, MESH_X_SIZE-1, MESH_Y_SIZE-1, vegetation, 0, trees);
	all_scenery.add_cobjs();
}


void add_plant(point const &pos, float height, float radius, int type, int calc_z) {
	all_scenery.add_plant(pos, height, radius, type, calc_z);
	has_scenery2 = 1;
}
void add_leafy_plant(point const &pos, float radius, int type, int calc_z) {
	all_scenery.add_leafy_plant(pos, radius, type, calc_z);
	has_scenery2 = 1;
}

void draw_scenery(bool shadow_only) {
	if (has_scenery || has_scenery2) {all_scenery.draw(shadow_only);}
}
void draw_scenery_fires(shader_t &s) {
	if (has_scenery || has_scenery2) {all_scenery.draw_fires(s);}
}

void add_scenery_cobjs() {all_scenery.add_cobjs();}

void shift_scenery(vector3d const &vd) {
	if (!has_scenery2) return; // dynamically created, not placed
	all_scenery.shift(vd);
}

// update region is inclusive: [x1,x2]x[y1,y2]
bool update_scenery_zvals(int x1, int y1, int x2, int y2) {
	return all_scenery.update_zvals(x1, y1, x2, y2);
}

void free_scenery_cobjs() {all_scenery.free_cobjs();}

void clear_scenery_vbos() {
	all_scenery.clear_vbos();
	voxel_rock_manager.free_context();
}

void do_rock_damage(point const &pos, float radius, float damage) {
	all_scenery.do_rock_damage(pos, radius, damage);
}

bool leafy_plant_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type) {
	all_scenery.leafy_plant_coll(index, energy);
	return 1;
}


void s_plant::write_to_cobj_file(std::ostream &out) const {
	// 'G': // place plant: xpos ypos height radius type [zpos], type: PLANT_MJ = 0, PLANT1, PLANT2, PLANT3, PLANT4
	//fscanf(fp, "%f%f%f%f%i%f", &pos.x, &pos.y, &fvals[0], &fvals[1], &ivals[0], &pos.z)
	//add_plant(pos, xf.scale*fvals[0], xf.scale*fvals[1], ivals[0], !use_z);
	out << "G " << pos.x << " " << pos.y << " " << height << " " << radius << " " << type << " " << pos.z << endl;
}
void scenery_group::write_plants_to_cobj_file(ostream &out) const { // Note: only plants are written because only plants can be placed in the cobj file
	for (auto i = plants.begin(); i != plants.end(); ++i) {i->write_to_cobj_file(out);}
}
void write_plants_to_cobj_file(ostream &out) {all_scenery.write_plants_to_cobj_file(out);}


