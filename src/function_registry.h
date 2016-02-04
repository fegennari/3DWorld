// 3D World - Global Function Registry Header
// by Frank Gennari
// 9/8/12

#ifndef _FUNCTION_REGISTRY_H_
#define _FUNCTION_REGISTRY_H_

#include "3DWorld.h"


// function prototypes - main (3DWorld.cpp, etc.)
bool check_gl_error(unsigned loc_id);
void enable_blend();
void disable_blend();
void set_std_blend_mode();
void set_additive_blend_mode();
void set_array_client_state(bool va, bool tca, bool na, bool ca, bool actually_set_state=1);
void reset_fog();
void set_perspective_near_far(float near_clip, float far_clip);
void set_perspective(float fovy, float nc_scale=1.0);
float get_star_alpha(bool obscured_by_clouds=0);
float get_moon_light_factor();
void setup_basic_fog();
void set_multisample(bool enable);
void check_zoom();
void reset_camera_pos();
void move_camera_pos_xy(vector3d const &v, float dist);
void move_camera_pos(vector3d const &v, float dist);
void update_cpos();
void advance_camera(int dir);
bool open_file(FILE *&fp, char const *const fn, std::string const &file_type, char const *const mode="r");
void fire_weapon();
bool has_extension(std::string const &ext);
bool is_shift_key_pressed();
bool is_ctrl_key_pressed();
bool is_alt_key_pressed();

// function prototypes - visibility
void calc_mesh_shadows(unsigned l, point const &lpos, float const *const mh, unsigned char *smask, int xsize, int ysize,
					   float const *sh_in_x=NULL, float const *sh_in_y=NULL, float *sh_out_x=NULL, float *sh_out_y=NULL);
void calc_visibility(unsigned light_sources);
bool is_visible_to_light_cobj(point const &pos, int light, float radius, int cobj, int skip_dynamic, int *cobj_ix=NULL);
bool coll_pt_vis_test(point pos, point pos2, float dist, int &index, int cobj, int skip_dynamic, int test_alpha);
void set_camera_pdu();
bool sphere_cobj_occluded(point const &viewer, point const &sc, float radius);
bool cube_cobj_occluded(point const &viewer, cube_t const &cube);
bool sphere_in_view(pos_dir_up const &pdu, point const &pos, float radius, int max_level, bool no_frustum_test=0);
int  get_light_pos(point &lpos, int light);
void update_sun_shadows();
void update_sun_and_moon();
bool light_valid(unsigned light_sources, int l, point &lpos);
bool light_valid_and_enabled(int l, point &lpos);

// function prototypes - mesh_intersect
bool sphere_visible_to_pt(point const &pt, point const &center, float radius);
bool is_visible_from_light(point const &pos, point const &lpos, int fast);
bool line_intersect_surface_cached(point const &v1, point const &v2, int &xpos, int &ypos, float &zval, int fast=0);
bool line_intersect_mesh(point const &v1, point const &v2, int &xpos, int &ypos, float &zval, int fast=0, bool cached=0);
bool line_intersect_mesh(point const &v1, point const &v2, int fast=0);
bool line_intersect_mesh(point const &v1, point const &v2, point &cpos, int fast=0, bool cached=0);
void gen_mesh_bsp_tree();

// function prototypes - build_world
void create_object_groups();
bool is_rain_enabled();
bool is_snow_enabled();
bool is_ground_wet();
bool is_ground_snowy();
void shift_all_objs(vector3d const &vd);
void process_platforms_falling_moving_and_light_triggers();
void set_global_state();
void process_groups();
void gen_scene(int generate_mesh, int gen_trees, int keep_sin_table, int update_zvals, int rgt_only);
void init_models();
void free_models();

// function prototypes - display_world
void glClearColor_rgba(const colorRGBA &color);
void set_standard_viewport();
point get_sun_pos();
point get_moon_pos();
colorRGBA get_bkg_color(point const &p1, vector3d const &v12);
void draw_scene_from_custom_frustum(pos_dir_up const &pdu, bool reflection_pass, bool include_mesh, bool disable_occ_cull);

// function prototypes - draw_world
void set_fill_mode();
void ensure_filled_polygons();
void reset_fill_mode();
int get_universe_ambient_light(bool for_universe_draw);
void set_gl_light_pos(int light, point const &pos, float w, shader_t *shader=NULL);
void set_light_ds_color(int light, colorRGBA const &diffuse, shader_t *shader=NULL);
void set_light_a_color(int light, colorRGBA const &ambient, shader_t *shader=NULL);
void set_light_colors(int light, colorRGBA const &ambient, colorRGBA const &diffuse, shader_t *shader=NULL);
void set_colors_and_enable_light(int light, colorRGBA const &ambient, colorRGBA const &diffuse, shader_t *shader=NULL);
void clear_colors_and_disable_light(int light, shader_t *shader=NULL);
void setup_gl_light_atten(int light, float c_a, float l_a, float q_a, shader_t *shader=NULL);
int get_light();
void draw_camera_weapon(bool want_has_trans);
void draw_solid_object_groups();
void draw_transparent_object_groups();
void draw_select_groups(int solid);
colorRGBA get_powerup_color(int powerup);
void update_precip_rate(float val);
unsigned get_precip_rate();
float get_rain_intensity();
float get_snow_intensity();
bool is_light_enabled(int l);
void enable_light    (int l);
void disable_light   (int l);
colorRGBA get_glowing_obj_color(point const &pos, int time, int lifetime, float &stime, bool shrapnel_cscale, bool fade);
colorRGBA const &get_landmine_light_color(int time);
float get_landmine_sensor_height(float radius, int time);
colorRGBA get_plasma_color(float size);
bool set_dlights_booleans(shader_t &s, bool enable, int shader_type, bool no_dl_smap=0);
float setup_underwater_fog(shader_t &s, int shader_type);
unsigned get_sky_zval_texture();
void invalidate_snow_coverage();
void setup_smoke_shaders(shader_t &s, float min_alpha, int use_texgen, bool keep_alpha, bool indir_lighting, bool direct_lighting, bool dlights, bool smoke_en,
	int has_lt_atten=0, bool use_smap=0, int use_bmap=0, bool use_spec_map=0, bool use_mvm=0, bool force_tsl=0, float burn_tex_scale=0.0,
	float triplanar_texture_scale=0.0, bool use_depth_trans=0, bool enable_reflections=0, bool is_outside=0);
void set_tree_branch_shader(shader_t &s, bool direct_lighting, bool dlights, bool use_smap);
void setup_procedural_shaders(shader_t &s, float min_alpha, bool indir_lighting, bool dlights, bool use_smap, bool use_bmap, bool use_noise_tex,
	bool z_top_test, float tex_scale=1.0, float noise_scale=1.0, float tex_mix_saturate=1.0);
void setup_object_render_data();
void end_group(int &last_group_id);
void draw_coll_surfaces(bool draw_trans, bool reflection_pass);
void draw_stars(float alpha);
void draw_sun();
void draw_moon();
void draw_earth();
void apply_red_sky(colorRGBA &color);
colorRGBA get_cloud_color();
void get_avg_sky_color(colorRGBA &avg_color);
float get_cloud_density(point const &pt, vector3d const &dir);
void free_cloud_textures();
void draw_puffy_clouds(int order, bool no_update=0);
float get_cloud_zmax();
void set_cloud_uniforms(shader_t &s, unsigned tu_id);
void draw_cloud_planes(float terrain_zmin, bool reflection_pass, bool draw_ceil, bool draw_floor);
void draw_sky(bool camera_side, bool no_update=0);
void compute_brightness();
void setup_water_plane_texgen(float s_scale, float t_scale, shader_t &shader, int mode);
void draw_water_plane(float zval, float terrain_zmin, unsigned reflection_tid);
void draw_splashes();
void draw_bubbles();
void draw_cracks_and_decals();
void draw_smoke_and_fires();
void add_camera_filter(colorRGBA const &color, unsigned time, int tid, unsigned ix, bool fades=0);
void draw_camera_filters(vector<camera_filter> &cfs);
point world_space_to_screen_space(point const &pos);
void restore_prev_mvm_pjm_state();
bool is_sun_flare_visible();
void draw_projectile_effects();
void draw_splash(float x, float y, float z, float size, colorRGBA color=WATER_C);
void draw_framerate(float val);
void draw_compass_and_alt();
void exec_universe_text(std::string const &text);
void set_silver_material(shader_t &shader, float alpha=1.0);
void set_gold_material  (shader_t &shader, float alpha=1.0);
void set_copper_material(shader_t &shader, float alpha=1.0);
void set_brass_material (shader_t &shader, float alpha=1.0);

// function prototypes - draw shapes
bool is_above_mesh(point const &pos);
bool check_face_containment(cube_t const &cube, int dim, int dir, int cobj);
float get_mesh_zmax(point const *const pts, unsigned npts);
void add_shadow_obj(point const &pos, float radius, int coll_id);
void add_coll_shadow_objs();
void get_occluders();

// function prototypes - draw primitives
void get_ortho_vectors(vector3d const &v12, vector3d *vab, int force_dim=-1);
vector_point_norm const &gen_cylinder_data(point const ce[2], float radius1, float radius2, unsigned ndiv, vector3d &v12,
										   float const *const perturb_map=NULL, float s_beg=0.0, float s_end=1.0, int force_dim=-1);
void draw_cylinder(float length, float radius1, float radius2, int ndiv, bool draw_ends=0, bool first_end_only=0, bool last_end_only=0, float z_offset=0.0);
void draw_cylinder_at(point const &p1, float length, float radius1, float radius2, int ndiv, bool draw_ends=0, bool first_end_only=0, bool last_end_only=0);
void draw_circle_normal(float r_inner, float r_outer, int ndiv, int invert_normals, point const &pos);
void draw_circle_normal(float r_inner, float r_outer, int ndiv, int invert_normals, float zval=0.0);
void begin_cylin_vertex_buffering();
void flush_cylin_vertex_buffer();
void gen_cone_triangles(vector<vert_norm_tc> &verts, vector_point_norm const &vpn, bool two_sided_lighting=0, float tc_t0=0.0, float tc_t1=1.0, vector3d const &xlate=zero_vector);
void gen_cylinder_triangle_strip(vector<vert_norm_tc> &verts, vector_point_norm const &vpn, bool two_sided_lighting=0, float tc_t0=0.0, float tc_t1=1.0, vector3d const &xlate=zero_vector);
void draw_fast_cylinder(point const &p1, point const &p2, float radius1, float radius2, int ndiv, bool texture,
	int draw_sides_ends=0, bool two_sided_lighting=0, float const *const perturb_map=NULL, float tex_scale_len=1.0, float tex_t_start=0.0, point const *inst_pos=NULL, unsigned num_insts=0);
void draw_cylindrical_section(float length, float r_inner, float r_outer, int ndiv, bool texture=0, float tex_scale_len=1.0, float z_offset=0.0);
void draw_cube_mapped_sphere(point const &center, float radius, unsigned ndiv, bool texture=0);
void get_sphere_triangles(vector<vert_wrap_t> &verts, point const &pos, float radius, int ndiv);
void draw_subdiv_sphere(point const &pos, float radius, int ndiv, point const &vfrom, float const *perturb_map,
						int texture, bool disable_bfc, unsigned char const *const render_map=NULL, float const *const exp_map=NULL,
						point const *const pt_shift=NULL, float expand=0.0, float s_beg=0.0, float s_end=1.0, float t_beg=0.0, float t_end=1.0);
void draw_subdiv_sphere(point const &pos, float radius, int ndiv, int texture, bool disable_bfc);
void draw_subdiv_sphere_section(point const &pos, float radius, int ndiv, int texture,
								float s_beg, float s_end, float t_beg, float t_end);
void rotate_sphere_tex_to_dir(vector3d const &dir);
void draw_single_colored_sphere(point const &pos, float radius, int ndiv, colorRGBA const &color);
void draw_torus(float ri, float ro, unsigned ndivi, unsigned ndivo, float tex_scale_i=1.0, float tex_scale_o=1.0);
void rotate_towards_camera(point const &pos);
void enable_flares(int tid);
void disable_flares();
void draw_tquad(float xsize, float ysize, float z, int prim_type=GL_TRIANGLE_FAN);
void draw_one_tquad(float x1, float y1, float x2, float y2, float z, int prim_type=GL_TRIANGLE_FAN);
int get_line_as_quad_pts(point const &p1, point const &p2, float w1, float w2, point pts[4]);
void draw_simple_cube(cube_t const &c, bool texture);
void draw_cube(point const &pos, float sx, float sy, float sz, bool texture, bool scale_ndiv=0,
			   float texture_scale=1.0, bool proportional_texture=0, vector3d const *const view_dir=NULL);
void gen_quad_tex_coords(float *tdata, unsigned num, unsigned stride);
void gen_quad_tri_tex_coords(float *tdata, unsigned num, unsigned stride);
void free_sphere_vbos();
void setup_sphere_vbos();
void draw_cylin_fast(float r1, float r2, float l, int ndiv, bool texture, float tex_scale_len=1.0, float z_offset=0.0);
void begin_sphere_draw(bool textured);
void end_sphere_draw();
void bind_draw_sphere_vbo(bool textured, bool normals=1);
void draw_sphere_vbo_pre_bound(int ndiv, bool textured, bool half=0, unsigned num_instances=1);
void draw_sphere_vbo_raw(int ndiv, bool textured, bool half=0, unsigned num_instances=1);
void draw_sphere_vbo(point const &pos, float radius, int ndiv, bool textured, bool half=0, bool bfc=0);
void draw_sphere_vbo_back_to_front(point const &pos, float radius, int ndiv, bool textured, bool enable_front=1, bool enable_back=1);

// function prototypes - draw mesh
float integrate_water_dist(point const &targ_pos, point const &src_pos, float const water_z);
void water_color_atten_pt(float *c, int x, int y, point const &pos, point const &p1, point const &p2);
void set_landscape_texgen(float tex_scale, int xoffset, int yoffset, int xsize, int ysize, shader_t &shader, unsigned detail_tu_id);
void display_mesh(bool shadow_pass=0, bool no_update=0);
void draw_water_sides(shader_t &shader, int check_zvals);
float get_tt_fog_top();
float get_tt_fog_bot();
float get_tt_cloud_level();
float get_inf_terrain_fog_dist();
float get_tt_fog_based_far_clip(float min_camera_dist);

// function prototypes - tiled mesh
bool using_tiled_terrain_hmap_tex();
float get_tiled_terrain_height_tex(float xval, float yval);
vector3d get_tiled_terrain_height_tex_norm(int x, int y);
bool write_default_hmap_modmap();
float update_tiled_terrain(float &min_camera_dist);
void pre_draw_tiled_terrain();
void draw_tiled_terrain(bool reflection_pass);
void draw_tiled_terrain_lightning(bool reflection_pass);
void end_tiled_terrain_lightning();
void draw_tiled_terrain_clouds(bool reflection_pass);
void clear_tiled_terrain();
void reset_tiled_terrain_state();
void clear_tiled_terrain_shaders();
float get_tiled_terrain_water_level();
void update_tiled_terrain_grass_vbos();
void draw_tiled_terrain_water(shader_t &s, float zval);
bool check_player_tiled_terrain_collision();
bool line_intersect_tiled_mesh(point const &v1, point const &v2, point &p_int);
void change_inf_terrain_fire_mode(int val);
void inf_terrain_fire_weapon();
void inf_terrain_undo_hmap_mod();
void setup_tt_fog_pre(shader_t &s);
void setup_tt_fog_post(shader_t &s);

// function prototypes - precipitation
void draw_local_precipitation();

// function prototypes - map_view
void draw_overhead_map();

// function prototypes - gen_obj
void gen_and_draw_stars(float alpha, bool half_sphere=0, bool no_update=0);
void gen_star(star &star1, int half_sphere=0);
void rand_xy_point(float zval, point &pt, unsigned flags);
void gen_object_pos(point &position, unsigned flags);
void gen_bubble(point const &pos, float r=0.0, colorRGBA const &c=WATER_C);
void gen_line_of_bubbles(point const &p1, point const &p2, float r=0.0, colorRGBA const &c=WATER_C);
bool gen_arb_smoke(point const &pos, colorRGBA const &bc, vector3d const &iv, float r, float den, float dark, float dam,
	int src, int dt, bool as, float spread=1.0);
void gen_smoke(point const &pos, float zvel_scale=1.0, float radius_scale=1.0);
bool gen_fire(point const &pos, float size, int source, bool allow_close=0, bool is_static=0, float light_bwidth=1.0, float intensity=1.0);
void gen_decal(point const &pos, float radius, vector3d const &orient, int tid, int cid=-1, colorRGBA const &color=BLACK,
	bool is_glass=0, bool rand_angle=0, int lifetime=60*TICKS_PER_SECOND, tex_range_t const &tr=tex_range_t());
void gen_particles(point const &pos, unsigned num, float lt_scale=1.0, bool fade=0);
int gen_fragment(point const &pos, vector3d const &velocity, float size_mult, float time_mult,
	colorRGBA const &color, int tid, float tscale, int source, bool tri_fragment, float hotness=0.0);
void gen_leaf_at(point const *const points, vector3d const &normal, int type, colorRGB const &color);
void add_water_particles(point const &pos, vector3d const &vadd, float vmag, float gen_radius, float mud_mix, float blood_mix, unsigned num);
void add_explosion_particles(point const &pos, vector3d const &vadd, float vmag, float gen_radius, colorRGBA const &color, unsigned num);
void gen_gauss_rand_arr();

// function prototypes - mesh_gen
bool bmp_to_chars(char const *const fname, unsigned char **&data);
void gen_mesh(int surface_type, int keep_sin_table, int update_zvals);
float do_glaciate_exp(float value);
float get_rel_wpz();
void init_terrain_mesh();
float eval_mesh_sin_terms(float xv, float yv);
float get_exact_zval(float xval, float yval);
void reset_offsets();
float get_median_height(float distribution_pos);
float get_water_z_height();
float get_cur_temperature();
void update_mesh(float dms, bool do_regen_trees);
bool is_under_mesh(point const &p);
bool read_mesh(const char *filename, float zmm=0.0);
bool write_mesh(const char *filename);
bool load_state(const char *filename);
bool save_state(const char *filename);

// function prototypes - physics
float get_max_t(int obj_type);
void init_objects();
void set_coll_rmax(float rmax);
void change_timestep(float mult_factor);
vector3d get_local_wind(int xpos, int ypos, float zval, bool no_use_mesh=0);
vector3d get_local_wind(point const &pt, bool no_use_mesh=0);
void reanimate_objects();
void seed_water_on_mesh(float amount);
void accumulate_object(point const &pos, int type, float amount);
void shift_other_objs(vector3d const &vd);
void advance_physics_objects();
void reset_other_objects_status();
void auto_advance_time();

// function prototypes - ai
void advance_smiley(dwobject &obj, int smiley_id);
void shift_player_state(vector3d const &vd, int smiley_id);
void player_clip_to_scene(point &pos);
void smiley_action(int smiley_id);

// function prototypes - matrix
void set_scene_constants();
void alloc_matrices();
void delete_matrices();
void compute_matrices();
void update_matrix_element(int xpos, int ypos);
void update_mesh_height(int xpos, int ypos, int rad, float scale, float offset, int mode, bool is_large_change);
vector3d get_matrix_surf_norm(float **matrix, unsigned char **enabled, int xsize, int ysize, int x, int y);
void calc_matrix_normal_at(float **matrix, vector3d **vn, vector3d **sn, unsigned char **enabled, int xsize, int ysize, int xpos, int ypos);
void calc_matrix_normals(float **matrix, vector3d **vn, vector3d **sn, unsigned char **enabled, int xsize, int ysize);
void get_matrix_point(int xpos, int ypos, point &pt);
int  is_in_ice(int xpos, int ypos);
float interpolate_mesh_zval(float xval, float yval, float rad, int use_real_equation, int ignore_ice, bool clamp_xy=0);
float int_mesh_zval_pt_off(point const &pos, int use_real_equation, int ignore_ice);
void calc_motion_direction();
float lowest_mesh_point(point const &pt, float radius);
float highest_mesh_point(point const &pt, float radius);

// function prototypes - collision detection
void reserve_coll_objects(unsigned size);
bool swap_and_set_as_coll_objects(coll_obj_group &new_cobjs);
int  add_coll_cube(cube_t &cube, cobj_params const &cparams, int platform_id=-1, int dhcm=0);
int  add_coll_cylinder(point const &p1, point const &p2, float radius, float radius2, cobj_params const &cparams, int platform_id=-1, int dhcm=0);
int  add_coll_capsule (point const &p1, point const &p2, float radius, float radius2, cobj_params const &cparams, int platform_id=-1, int dhcm=0);
int  add_coll_sphere(point const &pt, float radius, cobj_params const &cparams, int platform_id=-1, int dhcm=0);
int  add_coll_polygon(const point *points, int npoints, cobj_params const &cparams, float thickness, int platform_id=-1, int dhcm=0);
int  add_simple_coll_polygon(const point *points, int npoints, cobj_params const &cparams, vector3d const &normal, int dhcm=0);
int  remove_coll_object(int index, bool reset_draw=1);
int  remove_reset_coll_obj(int &index);
void purge_coll_freed(bool force);
void remove_all_coll_obj();
void cobj_stats();
int  collision_detect_large_sphere(point &pos, float radius, unsigned flags);
int  check_legal_move(int x_new, int y_new, float zval, float radius, int &cindex);
bool is_point_interior(point const &pos, float radius);
bool decal_contained_in_cobj(coll_obj const &cobj, point const &pos, vector3d const &norm, float radius, int dir);
void gen_explosion_decal(point const &pos, float radius, vector3d const &coll_norm, coll_obj const &cobj, int dir);
void proc_moving_cobjs();

// function prototypes - objects
void free_cobj_draw_group_vbos();

// function prototypes - coll_cell_search
void build_static_moving_cobj_tree();
void build_cobj_tree(bool dynamic=0, bool verbose=1);
bool check_coll_line_exact_tree(point const &p1, point const &p2, point &cpos, vector3d &cnorm,
	int &cindex, int ignore_cobj, bool dynamic=0, int test_alpha=0, bool skip_non_drawn=0, bool include_voxels=1, bool skip_init_colls=0);
bool check_coll_line_tree(point const &p1, point const &p2, int &cindex, int ignore_cobj, bool dynamic=0, int test_alpha=0,
	bool skip_non_drawn=0, bool include_voxels=1, bool skip_init_colls=0);
bool cobj_contained_tree(point const &p1, point const &p2, point const &viewer, point const *const pts, unsigned npts,
	int ignore_cobj, int &cobj);
void get_coll_line_cobjs_tree(point const &pos1, point const &pos2, int ignore_cobj,
	vector<int> *cobjs, cobj_query_callback *cqc, bool dynamic, bool occlude);
void get_coll_sphere_cobjs_tree(point const &center, float radius, int cobj, vert_coll_detector &vcd, bool dynamic);
bool check_point_contained_tree(point const &p, int &cindex, bool dynamic);
bool have_occluders();
void get_intersecting_cobjs_tree(cube_t const &cube, vector<unsigned> &cobjs, int ignore_cobj, float toler,
	bool dynamic, bool check_ccounter, int id_for_cobj_int);
bool check_coll_line(point const &pos1, point const &pos2, int &cindex, int c_obj, int skip_dynamic, int test_alpha, bool include_voxels=1, bool skip_init_colls=0);
bool check_coll_line_exact(point pos1, point pos2, point &cpos, vector3d &coll_norm, int &cindex, float splash_val=0.0,
						   int ignore_cobj=-1, bool fast=0, bool test_alpha=0, bool skip_dynamic=0, bool include_voxels=1, bool skip_init_colls=0);
bool cobj_contained_ref(point const &pos1, point const &center, const point *pts, unsigned npts, int cobj, int &last_cobj);
bool cobj_contained(point const &pos1, point const &center, const point *pts, unsigned npts, int cobj);
bool is_cube_occluded(cube_t const &cube, point const &viewer);
bool is_occluded(vector<int> const &occluders, point const *const pts0, int npts, point const &camera);
void add_camera_cobj(point const &pos);
void force_onto_surface_mesh(point &pos);
int  set_true_obj_height(point &pos, point const &lpos, float step_height, float &zvel, int type, int id,
	bool flight, bool on_snow, bool skip_dynamic=0, bool test_only=0);

// function prototypes - math3d
float fix_angle(float angle);
void get_face_normal(shape3d &shape, int face_id);
void calc_reflection_angle(vector3d const &v_inc, vector3d &v_ref, vector3d const &norm);
bool calc_refraction_angle(vector3d const &v_inc, vector3d &v_ref, vector3d const &norm, float n1, float n2);
float get_fresnel_reflection(vector3d const &v_inc, vector3d const &norm, float n1, float n2);
float get_reflected_weight(float fresnel_ref, float alpha);
float get_coll_energy(vector3d const &v1, vector3d const &v2, float mass);
point triangle_centroid(point const &p1, point const &p2, point const &p3);
float triangle_area(point const &p1, point const &p2, point const &p3);
float polygon_area(point const *const points, unsigned npoints);
float get_closest_pt_on_line_t(point const &pos, point const &l1, point const &l2);
point get_closest_pt_on_line(point const &pos, point const &l1, point const &l2);
bool planar_contour_intersect(const point *points, unsigned npoints, point const &pos, vector3d const &norm);
bool point_in_polygon_2d(float xval, float yval, const point *points, int npts, int dx=0, int dy=1);
bool point_in_convex_planar_polygon(vector<point> const &pts, point const &normal, point const &pt);
bool get_poly_zminmax(point const *const pts, unsigned npts, vector3d const &norm, float dval,
					  cube_t const &cube, float &z1, float &z2);
bool get_poly_zvals(vector<tquad_t> const &pts, float xv, float yv, float &z1, float &z2);
void gen_poly_planes(point const *const points, unsigned npoints, vector3d const &norm, float thick, point pts[2][4]);
void thick_poly_to_sides(point const *const points, unsigned npoints, vector3d const &norm, float thick, vector<tquad_t> &sides);
bool line_int_plane(point const &p1, point const &p2, point const &pp0, vector3d const &norm, point &p_int, float &t, bool ignore_t);
bool line_int_plane_test_only(point const &p1, point const &p2, point const &pp0, vector3d const &norm);
bool thick_poly_intersect(vector3d const &v1, point const &p1, vector3d const &norm,
						  point const pts[2][4], bool test_side, unsigned npoints);
bool sphere_intersect_poly_sides(vector<tquad_t> const &pts, point const &center, float radius, float &dist, vector3d &norm, bool strict);
bool pt_line_seg_dist_less_than(point const &P, point const &L1, point const &L2, float dist);
bool sphere_poly_intersect(const point *points, unsigned npoints, point const &pos, vector3d const &norm, float rdist, float radius);
bool sphere_ext_poly_int_base(point const &pt, vector3d const &norm, point const &pos, float radius,
							  float thickness, float &thick, float &rdist);
bool sphere_ext_poly_intersect(point const *const points, unsigned npoints, vector3d const &norm,
							   point const &pos, float radius, float thickness, float t_adj);
bool sphere_test_comp(point const &p2, point const &p1, vector3d const &v1, float r2sq, float &t);
bool line_sphere_intersect(point const &p1, point const &p2, point const &c, float r);
bool circle_test_comp(point const &p2, point const &p1, vector3d const &v1, vector3d norm, float r2sq, float &t);
void dir_to_sphere_s_t(vector3d const &dir, vector3d const &sdir, double &s, double &t);
bool line_sphere_intersect_s_t(point const &p1, point const &p2, point const &sc, float radius,
							   vector3d const &sdir, double &s, double &t);
bool line_sphere_int(vector3d const &v1, point const &p1, point const &center, float radius, point &lsint, bool test_neg_t);
bool line_intersect_sphere(point const &p1, vector3d const &v12, point const &sc, float radius, float &rad, float &dist, float &t);
bool sphere_vert_cylin_intersect(point &center, float radius, cylinder_3dw const &c);
void get_sphere_border_pts(point *qp, point const &pos, point const &viewed_from, float radius, unsigned num_pts);
void get_sphere_points(point const &pos, float radius, point *pts, unsigned npts, vector3d const &dir);
bool line_torus_intersect(point const &p1, point const &p2, point const &tc, float ri, float ro, float &t);
bool sphere_torus_intersect(point const &sc, float sr, point const &tc, float ri, float ro, point &p_int, vector3d &norm, bool calc_int);
bool circle_rect_intersect(point const &pos, float radius, cube_t const &cube, int dim);
bool sphere_cube_intersect(point const &pos, float radius, cube_t const &cube);
bool sphere_cube_intersect(point const &pos, float radius, cube_t const &cube, point const &p_last,
						   point &p_int, vector3d &norm, unsigned &cdir, bool check_int, bool skip_z=0);
bool approx_poly_cylin_int(point const *const pts, unsigned npts, cylinder_3dw const &cylin);
bool do_line_clip(point &v1, point &v2, float const d[3][2]);
bool get_line_clip(point const &v1, point const &v2, float const d[3][2], float &tmin, float &tmax);
bool check_line_clip_expand(point const &v1, point const &v2, float const d[3][2], float expand);
float line_line_dist(point const &p1a, point const &p1b, point const &p2a, point const &p2b);
float get_cylinder_params(point const &cp1, point const &cp2, point const &pos, vector3d &v1, vector3d &v2);
int  line_intersect_trunc_cone(point const &p1, point const &p2, point const &cp1, point const &cp2,
							   float r1, float r2, bool check_ends, float &t, bool swap_ends=0);
bool line_intersect_cylinder(point const &p1, point const &p2, cylinder_3dw const &c, bool check_ends);
int  line_int_thick_cylinder(point const &p1, point const &p2, point const &cp1, point const &cp2,
							 float ri1, float ri2, float ro1, float ro2, bool check_ends, float &t);
bool cylin_proj_circle_z_SAT_test(point const &cc, float cr, point const &cp1, point const &cp2, float r1, float r2);
bool sphere_int_cylinder_pretest(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2,
								 bool check_ends, vector3d &v1, vector3d &v2, float &t, float &rad);
bool sphere_intersect_cylinder_ipt(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2,
							   bool check_ends, point &p_int, vector3d &norm, bool calc_int);
void cylinder_quad_projection(point *pts, cylinder_3dw const &c, vector3d const &v1, vector3d &v2, int &npts);
template<typename T> pointT<T> get_center_arb(pointT<T> const *const pts, int npts);
unsigned get_cube_corners(float const d[3][2], point corners[8], point const &viewed_from=all_zeros, bool all_corners=1);
void get_closest_cube_norm(float const d[3][2], point const &p, vector3d &norm);
void cylinder_bounding_sphere(point const *const pts, float r1, float r2, point &center, float &radius);
void polygon_bounding_sphere(const point *pts, int npts, float thick, point &center, float &radius);
void add_rotated_quad_pts(vert_norm *points, unsigned &ix, float theta, float z, point const &pos, float xscale1, float xscale2, float yscale, float zscale);
void vproj_plane(vector3d const &vin, vector3d const &n, vector3d &vout);
template<typename T> void rotate_vector3d(pointT<T> vin, pointT<T> const &vrot, double angle, pointT<T> &vout);
template<typename T> void rotate_vector3d_multi(pointT<T> const &vrot, double angle, pointT<T> *vout, unsigned nv);
void rotate_vector3d_x2(point const &vrot, double angle, point &vout1, point &vout2);
void rotate_vector3d_by_vr(vector3d v1, vector3d v2, vector3d &vout);
cube_t rotate_cube(cube_t const &cube, vector3d const &axis, float angle_in_radians);
void mirror_about_plane(vector3d const &norm, point const &pt);
float angle_of_projected_vectors(vector3d const &v1, vector3d const &v2, vector3d n);
vector3d rtp_to_xyz(float radius, double theta, double phi);
vector3d gen_rand_vector_uniform(float mag);
vector3d gen_rand_vector(float mag, float zscale=1.0, float phi_term=PI);
vector3d gen_rand_vector2(float mag, float zscale=1.0, float phi_term=PI);
vector3d lead_target(point const &ps, point const &pt, vector3d const &vs, vector3d const &vt, float vweap);
vector3d get_firing_dir(vector3d const &src, vector3d const &dest, float fvel, float gravity_scale);

// function prototypes - water
bool get_water_enabled(int x, int y);
bool has_water(int x, int y);
bool mesh_is_underwater(int x, int y);
void water_color_atten_at_pos(colorRGBA &c, point const &pos);
void select_water_ice_texture(shader_t &shader, colorRGBA &color);
void set_tt_water_specular(shader_t &shader);
colorRGBA get_tt_water_color();
void draw_water(bool no_update=0);
void add_splash(point const &pos, int xpos, int ypos, float energy, float radius, bool add_sound, vector3d const &vadd=zero_vector);
bool add_water_section(float x1, float y1, float x2, float y2, float zval, float wvol);
void float_downstream(point &pos, float radius);
void change_water_level(float water_level);
void calc_watershed();
bool is_underwater(point const &pos, int check_bottom=0, float *depth=NULL);
void select_liquid_color(colorRGBA &color, int xpos, int ypos);
void select_liquid_color(colorRGBA &color, point const &pos);
float get_blood_mix(point const &pos);
void add_water_spring(point const &pos, vector3d const &vel, float rate, float diff, int calc_z, int gen_vel);
void shift_water_springs(vector3d const &vd);
void update_water_zval(int x, int y, float old_mh);

// function prototypes - lightning
void compute_volume_matrix();

// function prototypes - textures
void load_texture_names();
void load_textures();
int texture_lookup(std::string const &name);
int get_texture_by_name(std::string const &name, bool is_normal_map=0, bool invert_y=0);
bool select_texture(int id);
void update_player_bbb_texture(float extra_blood, bool recreate);
float get_tex_ar(int id);
void bind_1d_texture(unsigned tid, bool is_array=0);
void bind_2d_texture(unsigned tid, bool is_array=0, bool multisample=0);
void setup_texture(unsigned &tid, bool mipmap, bool wrap_s, bool wrap_t,
	bool mirror_s=0, bool mirror_t=0, bool nearest=0, float anisotropy=1.0, bool is_array=0, bool multisample=0);
void setup_1d_texture(unsigned &tid, bool mipmap, bool wrap, bool mirror, bool nearest);
void depth_buffer_to_texture(unsigned &tid);
void frame_buffer_RGB_to_texture(unsigned &tid);
void free_textures();
void reset_textures();
void free_texture(unsigned &tid);
void setup_landscape_tex_colors(colorRGBA const &c1, colorRGBA const &c2);
colorRGBA texture_color(int tid);
unsigned get_texture_size(int tid, bool dim);
void get_lum_alpha(colorRGBA const &color, int tid, float &luminance, float &alpha);
std::string get_file_extension(std::string const &filename, unsigned level, bool make_lower);
unsigned get_noise_tex_3d(unsigned tsize, unsigned ncomp);
colorRGBA get_landscape_texture_color(int xpos, int ypos);
void update_lttex_ix(int &ix);
void get_tids(float relh, int &k1, int &k2, float *t=NULL);
void create_landscape_texture();
float add_crater_to_landscape_texture(float xval, float yval, float radius);
void add_color_to_landscape_texture(colorRGBA const &color, float xval, float yval, float radius);
void add_snow_to_landscape_texture(point const &pos, float acc);
void update_landscape_texture();
void gen_tex_height_tables();
void setup_texgen_full(float sx, float sy, float sz, float sw, float tx, float ty, float tz, float tw, shader_t &shader, int mode);
void setup_texgen(float xscale, float yscale, float tx, float ty, float z_off, shader_t &shader, int mode);
void get_poly_texgen_dirs(vector3d const &norm, vector3d v[2]);
void setup_polygon_texgen(vector3d const &norm, float const scale[2], float const xlate[2], vector3d const &offset, bool swap_txy, shader_t &shader, int mode);
void get_tex_coord(vector3d const &dir, vector3d const &sdir, unsigned txsize, unsigned tysize, int &tx, int &ty, bool invert);
float get_texture_component(unsigned tid, float u, float v, int comp);
colorRGBA get_texture_color(unsigned tid, float u, float v);
vector2d get_billboard_texture_uv(point const *const points, point const &pos);
bool is_billboard_texture_transparent(point const *const points, point const &pos, int tid);

// function prototypes - sun flares
void DoFlares(point const &from, point const &at, point const &light, float near_clip, float size, float intensity, int start_ix=0);
void load_flare_textures();
void free_flare_textures();

// function prototypes - gameplay/ai
bool camera_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool smiley_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool landmine_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool health_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool shield_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool powerup_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool weapon_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool ammo_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool pack_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool rock_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool sball_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool dodgeball_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool skull_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
bool sawblade_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);

void gen_rocket_smoke(point const &pos, vector3d const &orient, float radius);
void gen_landmine_scorch(point const &pos);
int get_smiley_hit(vector3d &hdir, int index);
void blast_radius(point const &pos, int type, int obj_index, int shooter, int chain_level);
void create_explosion(point const &pos, int shooter, int chain_level, float damage, float size, int type, bool cview);
void do_area_effect_damage(point const &pos, float effect_radius, float damage, int index, int source, int type);
void switch_player_weapon(int val);
void draw_beams();
void show_blood_on_camera();
void update_weapon_cobjs();
int select_dodgeball_texture(int shooter);
sphere_t get_weapon_bsphere(int weapon);
void draw_weapon_simple(point const &pos, vector3d const &dir, float radius, int cid, int wid, float scale, shader_t &shader);
void draw_weapon_in_hand(int shooter, shader_t &shader);
bool weap_has_transparent(int shooter);
void draw_scheduled_weapons();
void add_weapon_lights(int shooter);
void show_crosshair(colorRGBA const &color, int do_zoom);
void show_user_stats();
void show_other_messages();
void print_text_onscreen(std::string const &text, colorRGBA const &color, float size, int time, int priority=0);
void print_weapon(int weapon_id);
bool check_underwater(int who, float &depth);
void player_fall(int id);
void update_camera_velocity(vector3d const &v);
void init_game_state();
void gamemode_rand_appear();
bool has_invisibility(int id);
void init_smileys();
void init_game_mode();
void update_game_frame();
void change_game_mode();
void free_dodgeballs(bool camera, bool smileys);
int gen_smiley_or_player_pos(point &pos, int index);
colorRGBA get_smiley_team_color(int smiley_id, bool ignore_teams=0);
void select_smiley_texture(int smiley_id);
void free_smiley_textures();
void clear_cached_waypoints();
int get_ammo_or_obj(int wid);
int wid_need_weapon(int wid);
void draw_teleporters();
bool maybe_teleport_object(point &opos, float oradius, int player_id);
void player_teleported(point const &pos, int player_id);

// function prototypes - explosion
void update_blasts();
void draw_blasts();

// function prototypes - scenery
void gen_scenery();
void draw_scenery(bool shadow_only=0);
bool update_scenery_zvals(int x1, int y1, int x2, int y2);
void free_scenery();
void do_rock_damage(point const &pos, float radius, float damage);
void add_scenery_cobjs();
void shift_scenery(vector3d const &vd);
void add_plant(point const &pos, float height, float radius, int type, int calc_z);

// function prototypes - grass
void setup_wind_for_shader(shader_t &s, unsigned tu_id);
bool no_grass();
void gen_grass();
void update_grass_vbos();
void draw_grass();
void modify_grass_at(point const &pos, float radius, bool crush=0, bool burn=0, bool cut=0, bool check_uw=0, bool add_color=0, bool remove=0, colorRGBA const &color=BLACK);
void grass_mesh_height_change(int xpos, int ypos);
void flower_mesh_height_change(int xpos, int ypos, int rad);
bool place_obj_on_grass(point &pos, float radius);
float get_grass_density(point const &pos);

// function prototypes - draw mech
void build_hmv_shape();
void delete_hmv_shape();
void add_shape_coll_objs();
void shift_hmv(vector3d const &vd);

// function prototypes - tree + sm_tree (see also tree_3dw.h)
colorRGBA get_tree_trunk_color(int type, bool modulate_with_texture);
int get_tree_class_from_height(float zpos, bool pine_trees_only);
int get_tree_type_from_height(float zpos, rand_gen_t &rgen);
float get_plant_leaf_wind_mag(bool shadow_only);
void setup_leaf_wind(shader_t &s, float wind_mag, bool underwater);
void set_leaf_shader(shader_t &s, float min_alpha, unsigned tc_start_ix=0, bool enable_opacity=0, bool no_dlights=0, float wind_mag=0.0, bool underwater=0, bool use_fs_smap=0, bool enable_smap=1);
bool update_decid_tree_zvals(int x1, int y1, int x2, int y2);
bool update_small_tree_zvals(int x1, int y1, int x2, int y2);
void exp_damage_trees(point const &epos, float damage, float bradius, int type);

// function prototypes - ship
upos_point_type const &get_player_pos();
vector3d const &get_player_dir();
vector3d const &get_player_up();
vector3d const &get_player_velocity();
float get_player_radius();
void set_player_pos(point const &pos_);
void set_player_dir(vector3d const &dir_);
void set_player_up(vector3d const &upv_);
void stop_player_ship();
void init_universe_display();
void set_univ_pdu();
void setup_current_system(float sun_intensity=1.0);
void apply_univ_physics();
void draw_universe(bool static_only=0, bool skip_closest=0, int no_distant=0, bool gen_only=0);
void draw_universe_stats();
void clear_univ_obj_contexts();
void clear_cached_shaders();

// function prototypes - lightmap
void update_flow_for_voxels(vector<cube_t> const &cubes);
void shift_light_sources(vector3d const &vd);
void regen_lightmap();
void clear_lightmap();
void build_lightmap(bool verbose);
void add_smoke(point const &pos, float val);
void distribute_smoke();
float get_smoke_at_pos(point const &pos);
void add_line_light(point const &p1, point const &p2, colorRGBA const &color, float size, float intensity=1.0);
void add_dynamic_light(float sz, point const &p, colorRGBA const &c=WHITE, vector3d const &d=plus_z, float bw=1.0, point *line_end_pos=nullptr, bool is_static_pos=0);
colorRGBA gen_fire_color(float &cval, float &inten, float rate=1.0);
void clear_dynamic_lights();
void add_camera_flashlight();
void add_camera_candlelight();
void add_dynamic_lights_ground();
void get_xyz(point const &p, int v[3]);
void get_xyz_v2(point &p, int const v[3]);
void update_smoke_indir_tex_range(unsigned x_start, unsigned x_end, unsigned y_start, unsigned y_end, bool update_lighting);
bool upload_smoke_indir_texture();
void upload_dlights_textures(cube_t const &bounds);
void setup_dlight_textures(shader_t &s, bool enable_dlights_smap=1);
bool is_in_darkness(point const &pos, float radius, int cobj);
void get_indir_light(colorRGBA &a, point const &p);

// function prototypes - tessellate
void split_polygon_to_cobjs(coll_obj const &cobj, coll_obj_group &split_polygons, vector<point> const &poly_pt, bool split_quads);

// function prototypes - shaders
char const *append_ix(std::string &s, unsigned i, bool as_array);
bool setup_shaders();
void clear_shaders();
void reload_all_shaders();
void check_mvm_update();
void upload_mvm_to_shader(shader_t &s, char const *const var_name);
void set_point_sprite_mode(bool enabled);

// function prototypes - snow
bool snow_enabled();
void gen_snow_coverage();
void reset_snow_vbos();
void draw_snow(bool shadow_only=0);
bool get_snow_height(point const &p, float radius, float &zval, vector3d &norm, bool crush_snow=0);
bool crush_snow_at_pt(point const &p, float radius);

// function prototypes - waypoints
void create_waypoints(vector<user_waypt_t> const &user_waypoints);
void shift_waypoints(vector3d const &vd);
void draw_waypoints();

// function prototypes - destroy_cobj
void destroy_coll_objs(point const &pos, float damage, int shooter, int damage_type, float force_radius=0.0);
void check_falling_cobjs();

// function prototypes - shadow_map
cube_t get_scene_bounds();
bool shadow_map_enabled();
int get_def_smap_ndiv(float radius);
void set_smap_shader_for_all_lights(shader_t &s, float z_bias=DEF_Z_BIAS);
pos_dir_up get_pt_cube_frustum_pdu(point const &pos, cube_t const &bounds);
void draw_scene_bounds_and_light_frustum(point const &lpos);
void create_shadow_map();
void free_shadow_map_textures();

// function prototypes - raytrace
float get_scene_radius();
void kill_current_raytrace_threads();
void check_update_global_lighting(unsigned lights);

// function prototypes - voxels
void gen_voxel_landscape();
bool gen_voxels_from_cobjs(coll_obj_group &cobjs);
float gen_voxel_rock(voxel_model &model, point const &center, float radius, unsigned size, unsigned num_blocks=1, int rseed=456);
bool parse_voxel_option(FILE *fp);
void render_voxel_data(bool shadow_pass);
void free_voxel_context();
bool point_inside_voxel_terrain(point const &pos);
float get_voxel_terrain_ao_lighting_val(point const &pos);
bool update_voxel_sphere_region(point const &center, float radius, float val_at_center, int shooter, unsigned num_fragments=0);
void proc_voxel_updates();
bool check_voxel_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj, bool exact);
void get_voxel_coll_sphere_cobjs(point const &center, float radius, int ignore_cobj, vert_coll_detector &vcd);
bool write_voxel_brushes();
void change_voxel_editing_mode(int val);
void undo_voxel_brush();
void modify_voxels();

// function prototypes - screenshot
void read_depth_buffer(unsigned window_width, unsigned window_height, vector<float> &depth, bool normalize=0);
void read_pixels(unsigned window_width, unsigned window_height, vector<unsigned char> &buf);
int screenshot(unsigned window_width, unsigned window_height, char *file_path, bool write_bmp);
int write_jpeg(unsigned window_width, unsigned window_height, char *file_path);

// function prototypes - spray paint
void toggle_spraypaint_mode();
void change_spraypaint_color(int val);
void draw_spraypaint_crosshair();
void spray_paint(bool mode);

// function prototypes - edit_ui
void next_selected_menu_ix();
bool ui_intercept_keyboard(unsigned char key, bool is_special);
bool ui_intercept_mouse(int button, int state, int x, int y, bool is_up_down);
void draw_enabled_ui_menus();

// function prototypes - transform_obj
void fgMatrixMode(int val);
void fgPushMatrix();
void fgPopMatrix();
void fgLoadIdentity();
void fgPushIdentityMatrix();
void fgTranslate(float x, float y, float z);
void fgScale(float x, float y, float z);
void fgScale(float s);
void fgRotate(float angle, float x, float y, float z);
void fgRotateRadians(float angle, float x, float y, float z);
void fgPerspective(float fov_y, float aspect, float near_clip, float far_clip);
void fgOrtho(float left, float right, float bottom, float top, float zNear, float zFar);
void fgLookAt(float eyex, float eyey, float eyez, float centerx, float centery, float centerz, float upx, float upy, float upz);
void deform_obj(dwobject &obj, vector3d const &norm, vector3d const &v0);
void update_deformation(dwobject &obj);

// function prototypes - draw_text
void load_font_texture_atlas(std::string const &fn="");
void free_font_texture_atlas();
void draw_text(colorRGBA const &color, float x, float y, float z, char const *text, float tsize=1.0);

// function prototypes - postproc_effects
void set_xy_step(shader_t &s);
void setup_depth_tex(shader_t &s, int tu_id);

// function prototypes - video_capture
void start_video_capture(std::string const &fn);
void end_video_capture();
void toggle_video_capture();
void video_capture_end_frame();

#include "inlines.h"


#endif // _FUNCTION_REGISTRY_H_

