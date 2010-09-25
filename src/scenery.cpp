// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 12/5/02

#include "3DWorld.h"
#include "mesh.h"
#include "shape_line3d.h"
#include "upsurface.h"

using namespace std;

bool     const NO_ISLAND_SCENERY = 1;
unsigned const ROCK_NDIV         = 24;
float    const SHADOW_VAL        = 0.5;
float    const PT_LINE_THRESH    = 800.0;


enum {PLANT_MJ = 0, PLANT1, PLANT2, PLANT3, PLANT4, COFFEE, NUM_PLANT_TYPES};

struct plant_type {

	int tid;
	colorRGBA stemc, leafc;

	plant_type(int tid_, colorRGBA const &sc, colorRGBA const &lc) : tid(tid_), stemc(sc), leafc(lc) {}
};

colorRGBA const stem_c(0.4, 0.6, 0.2, 1.0);
colorRGBA const leaf_c(0.7, 0.7, 0.7, 1.0);

// tid, stemc, leafc
plant_type const pltype[NUM_PLANT_TYPES] = {
	plant_type(MJ_LEAF_TEX, stem_c,   leaf_c),
	plant_type(PLANT1_TEX,  stem_c,   leaf_c),
	plant_type(PLANT2_TEX,  stem_c,   leaf_c),
	plant_type(PLANT3_TEX,  stem_c,   leaf_c),
	plant_type(PLANT4_TEX,  stem_c,   leaf_c),
	plant_type(COFFEE_TEX,  LT_BROWN, WHITE)
};

colorRGBA const log_colors[4] = {BLACK, PTREE_C, TREE_C, TREE_C};


int DISABLE_SCENERY(0), has_scenery(0), has_scenery2(0);
float msms2(1.0);


extern int num_trees, xoff2, yoff2, rand_gen_index, island, window_width, do_zoom, display_mode, DISABLE_WATER;
extern long rseed1, rseed2;
extern float zmin, zmax_est, water_plane_z, mesh_scale, mesh_scale2, vegetation, max_water_height;
extern GLUquadricObj* quadric;
extern pt_line_drawer tree_scenery_pld; // we can use this for plant trunks


void gen_scenery_deterministic();

inline float get_pt_line_thresh() {return PT_LINE_THRESH*(do_zoom ? ZOOM_FACTOR : 1.0);}


// ************ SCENERY OBJECT CLASSES ************


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


bool scenery_obj::in_camera_view(float brad) const {

	return sphere_in_camera_view(pos, ((brad == 0.0) ? radius : brad), 0);
}


float scenery_obj::get_shadowed_color(point const &p, float eff_radius) const { // not used on rock_shapes

	return (is_visible_to_light_cobj(p, get_light(), eff_radius, coll_id, 0) ? 1.0 : SHADOW_VAL);
}


colorRGBA scenery_obj::get_atten_color(colorRGBA c) const {

	if (!DISABLE_WATER && (pos.z + radius) < max_water_height) {
		int const x(get_xpos(pos.x)), y(get_ypos(pos.y));
		
		if (!point_outside_mesh(x, y) && (pos.z + radius) < water_matrix[y][x]) {
			water_color_atten(((float *)&(c.red)), x, y, (pos + point(0.0, 0.0, radius)));
		}
	}
	return c;
}


void scenery_obj::remove_cobjs() {

	remove_reset_coll_obj(coll_id);
}


void rock_shape3d::create(int x, int y, bool use_xy) {

	int rs_rock(rand2());
	gen_spos(x, y, use_xy);
	gen_rock(48, 0.05/msms2, rs_rock, (rand2()&1));
	radius = 0.0;
	
	for (unsigned i = 0; i < points.size(); ++i) { // calculate radius
		float const dist(points[i].mag()); // check if correct - might need pos or average/center calculation
		radius = max(radius, dist);
	}
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


void rock_shape3d::gen_rock(unsigned nverts, float size, int &rand_seed, int type) {

	rseed1    = rand_seed;
	rseed2    = 10423232; // whatever
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
			get_triangle_center(center, face_id, 1);

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
		unsigned face(0), tot_used(0);
		vector<bool> used(nverts, 0);
		edge_seen_set edges_seen;
		deque<edge> edges; // incomplete faces

		for (unsigned cv = 0; cv < nverts; ++cv) { // is this outer loop necessary?
			if (used[cv]) continue; // finished with this vertex
			unsigned imin(0);
			float dmin(0.0);

			for (unsigned i = 0; i < nverts; ++i) {
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
	radius *= 0.99; // chip off parts of the rock to make it smaller
	scale  *= 0.99;
	remove_cobjs();
	add_cobjs();
	return 1;
}


void rock_shape3d::draw() const {

	if (in_camera_view()) {
		set_color(get_atten_color(color*get_shadowed_color(pos, 0.5*radius)));
		shape3d::draw(1);
	}
}


struct surface_cache {

	typedef pair<long, long> seed_pair;
	typedef map<seed_pair, upsurface *> surface_map;
	surface_map scache;

	upsurface *get_surface() {
		seed_pair const sp(rseed1, rseed2);
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

	void clear() {
		for (surface_map::const_iterator i = scache.begin(); i != scache.end(); ++i) {
			assert(i->second);
			delete i->second;
		}
		scache.clear();
	}

	void clear_unref() {
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
};

surface_cache surface_rock_cache;


class surface_rock : public scenery_obj { // size = 1456+

	float scale;
	vector3d dir;
	upsurface *surface;

public:
	surface_rock() : surface(NULL) {}

	void create(int x, int y, int use_xy) {
		gen_spos(x, y, use_xy);
		radius  = 0.5*rand_uniform2(0.2/msms2, 0.8/msms2)*rand_float2();
		dir     = signed_rand_vector2_norm();
		surface = surface_rock_cache.get_surface();
		assert(surface);

		if (surface->ssize == 0) { // not inited
			surface->gen(0.5, rand_uniform2(0.5, 5.0), 10, rand_uniform2(0.5, 2.0));
			surface->setup(ROCK_NDIV, 0.0, 0);
			surface->setup_draw_sphere(all_zeros, rand_uniform2(0.25, 1.0), 0.0, ROCK_NDIV, NULL);
			surface->calc_rmax();
		}
		scale   = radius/surface->rmax;
	}

	void add_cobjs() {
		coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, BROWN, 0, 0, rock_collision, 1, ROCK_SPHERE_TEX));
	}

	void draw(float sscale) const {
		assert(surface);
		if (!in_camera_view(0.0)) return;
		colorRGBA const color(get_atten_color(WHITE)*get_shadowed_color(pos, radius));
		float const dist(distance_to_camera(pos));

		if (2*get_pt_line_thresh()*radius < dist) { // draw as point
			tree_scenery_pld.add_textured_pt(pos, color, ROCK_SPHERE_TEX);
			return;
		}
		set_color(color);
		select_texture(ROCK_SPHERE_TEX);
		glPushMatrix();
		translate_to(pos);
		uniform_scale(scale);
		rotate_into_plus_z(dir);
		surface->sd.draw_ndiv_pow2(sscale*radius/dist);
		glPopMatrix();
	}

	void destroy() {
		//delete surface;
		if (surface) surface->dec_ref();
		surface = NULL;
		scenery_obj::destroy();
	}
};


class s_rock : public scenery_obj { // size = 48

	float size, angle;
	vector3d scale, dir;

public:
	void create(int x, int y, int use_xy) {
		for (unsigned i = 0; i < 3; ++i) {
			scale[i] = rand_uniform2(0.8, 1.3);
		}
		gen_spos(x, y, use_xy);
		size   = 0.02*rand_uniform2(0.2/msms2, 0.8/msms2);
		if ((rand2()&3) == 0) size *= rand_uniform2(1.2, 8.0);
		dir    = signed_rand_vector2_norm();
		angle  = rand_uniform2(0.0, 360.0);
		radius = size*(scale.x + scale.y + scale.z)/3.0;
		pos.z += radius*rand_uniform2(-0.1, 0.25);
	}

	void add_cobjs() {
		coll_id = add_coll_sphere(pos, radius, cobj_params(0.95, BROWN, 0, 0, rock_collision, 1, ROCK_SPHERE_TEX));
	}

	void draw(float sscale) const {
		float const rmax(1.3*radius);
		if (!in_camera_view(rmax)) return;
		colorRGBA const color(get_atten_color(WHITE)*get_shadowed_color(pos, rmax));
		float const dist(distance_to_camera(pos));

		if (2*get_pt_line_thresh()*radius < dist) { // draw as point
			tree_scenery_pld.add_textured_pt(pos, color, ROCK_SPHERE_TEX);
			return;
		}
		set_color(color);
		int const ndiv(max(4, min(N_SPHERE_DIV, int(sscale*radius/dist))));
		select_texture(ROCK_SPHERE_TEX);
		glPushMatrix();
		translate_to(pos);
		rotate_about(angle, dir);
		scale_by(scale);
		draw_sphere_dlist(all_zeros, size, ndiv, 1);
		glPopMatrix();
	}
};


class s_log : public scenery_obj { // size = 57 (60)

	float length, radius2;
	point pt2;
	vector3d dir;

public:
	void shift_by(vector3d const &vd) {
		scenery_obj::shift_by(vd);
		pt2   += vd;
		pt2.z -= dz;
		dz     = 0.0;
	}

	int create(int x, int y, int use_xy, float minz) {
		gen_spos(x, y, use_xy);
		radius  = rand_uniform2(0.003/msms2, 0.008/msms2);
		radius2 = rand_uniform2(0.9*radius, 1.1*radius);
		length  = rand_uniform2(max(0.03/msms2, 4.0*radius), min(0.15/msms2, 20.0*radius));
		dir     = signed_rand_vector2_norm();
		dir.x  *= length;
		dir.y  *= length;
		pt2.x   = pos.x + dir.x;
		pt2.y   = pos.y + dir.y;

		if (pt2.x > X_SCENE_SIZE-DX_VAL || pt2.x < -X_SCENE_SIZE+DX_VAL || pt2.y > Y_SCENE_SIZE-DY_VAL || pt2.y < -Y_SCENE_SIZE+DY_VAL) {
			return 0; // off the end of the mesh
		}
		pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1) + rand_uniform2(0.7, 0.99)*radius;
		pt2.z = interpolate_mesh_zval(pt2.x, pt2.y, 0.0, 1, 1) + rand_uniform2(0.7, 0.99)*radius2;
		if (max(pos.z, pt2.z) < minz)       return 0;
		if (pos.z <= zmin || pt2.z <= zmin) return 0; // bad z value
		dir.z  = (pt2.z - pos.z);
		length = dir.mag(); // recalculate
		dir   /= -length; // something is backwards
		type   = (char)get_tree_type_from_height(max(pos.z, pt2.z));
		return 1;
	}

	void add_cobjs() {
		coll_id = add_coll_cylinder(pos, pt2, radius, radius2, cobj_params(0.8, BROWN, 0, 0, NULL, 0, WOOD_TEX));
	}

	void draw(float sscale) const {
		float const sz(max(length, max(radius, radius2)));
		if (type == 0 || !in_camera_view(sz)) return;
		colorRGBA const color(log_colors[type]*get_shadowed_color((pos + pt2)*0.5, sz));
		float const dist(distance_to_camera(pos));

		if (get_pt_line_thresh()*(radius + radius2) < dist) { // draw as line
			tree_scenery_pld.add_textured_line(pos, pt2, color, WOOD_TEX);
			return;
		}
		set_color(color);
		int const ndiv(max(3, min(N_CYL_SIDES, int(2.0*sscale*radius/dist))));
		select_texture(WOOD_TEX);
		glPushMatrix();
		translate_to(pos);
		rotate_by_vector(dir, 0.0);
		gluCylinder(quadric, radius, radius2, length, ndiv, 1);
		select_texture(TREE_END_TEX);
		draw_circle_normal(0.0, radius, ndiv, 1);
		glTranslatef(0.0, 0.0, length);
		draw_circle_normal(0.0, radius2, ndiv, 0);
		glPopMatrix();
	}

	bool update_zvals(int x1, int y1, int x2, int y2) {
		float const orig_pz(pos.z);
		if (!scenery_obj::update_zvals(x1, y1, x2, y2)) return 0;
		pt2.z += (pos.z - orig_pz); // apply the same offset to pt2 even though it might be at a different mesh location
		return 1;
	}
};


class s_stump : public scenery_obj { // size = 29 (32)

	float radius2, height;

public:
	int create(int x, int y, int use_xy, float minz) {
		gen_spos(x, y, use_xy);
		radius  = rand_uniform2(0.005/msms2, 0.01/msms2);
		radius2 = rand_uniform2(0.8*radius, radius);
		pos.z  -= 2.0*radius;
		if (pos.z < minz) return 0;
		height  = rand_uniform2(0.01/msms2, min(0.05/msms2, 4.0*radius)) + 0.015;
		
		if ((rand2()&3) == 0) {
			height  *= rand_uniform2(1.0, 5.0); // larger stump = upright dead tree
			radius  *= 1.5;
			radius2 *= 1.3;
		}
		type = (char)get_tree_type_from_height(pos.z);
		return 1;
	}

	void add_cobjs() {
		coll_id = add_coll_cylinder(pos, point(pos.x, pos.y, (pos.z + height)), radius, radius2,
			cobj_params(0.8, BROWN, 0, 0, NULL, 0, WOOD_TEX));
	}

	void draw(float sscale) const {
		float const sz(max(height, max(radius, radius2)));
		if (type == 0 || !in_camera_view(sz)) return;
		colorRGBA const color(log_colors[type]*get_shadowed_color(point(pos.x, pos.y, (pos.z + 0.5*height)), sz));
		float const dist(distance_to_camera(pos));

		if (get_pt_line_thresh()*(radius + radius2) < dist) { // draw as line
			tree_scenery_pld.add_textured_line(pos, (pos + point(0.0, 0.0, height)), color, WOOD_TEX);
			return;
		}
		set_color(color);
		int const ndiv(max(3, min(N_CYL_SIDES, int(2.2*sscale*radius/dist))));
		select_texture(WOOD_TEX);
		glPushMatrix();
		translate_to(pos);
		gluCylinder(quadric, radius, radius2, height, ndiv, 1);
		select_texture(TREE_END_TEX);
		glTranslatef(0.0, 0.0, height);
		gluDisk(quadric, 0.0, radius2, ndiv, 1);
		glPopMatrix();
	}
};


class s_plant : public scenery_obj { // size = 40

	int coll_id2;
	float height;
	vector<point> points;

public:
	s_plant() : coll_id2(-1), height(1.0) {}
	bool operator<(s_plant const &p) const {return (type < p.type);}

	int create(int x, int y, int use_xy, float minz) {
		points.clear();
		type   = rand2()%NUM_PLANT_TYPES;
		gen_spos(x, y, use_xy);
		if (pos.z < minz) return 0;
		radius = rand_uniform2(0.0025/msms2, 0.0045/msms2);
		height = rand_uniform2(0.2/msms2, 0.4/msms2) + 0.025;
		return 1;
	}

	void create2(point const &pos_, float height_, float radius_, int type_, int calc_z) {
		points.clear();
		type   = abs(type_)%NUM_PLANT_TYPES;
		pos    = pos_;
		radius = radius_;
		height = height_;
		if (calc_z) pos.z = interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1);
	}

	void add_cobjs() {
		point cpos(pos), cpos2(pos);
		float const wscale(radius*msms2/0.004);
		float const r2(radius+0.07*wscale*(height+0.03));
		cobj_params const cp(0.4, pltype[type].leafc.modulate_with(texture_color(pltype[type].tid)), 0, 0, NULL, 0, -1, 1.0, 0, 1);
		cpos.z  += height;
		cpos2.z += 3.0*height/(36.0*height + 4.0);
		coll_id  = add_coll_cylinder(pos,   cpos, radius, radius, cp); // trunk
		coll_id2 = add_coll_cylinder(cpos2, cpos, r2,     radius, cp); // leaves
	}

	void gen_points() {
		float const wscale(250.0*radius*msms2);
		float const ms(mesh_scale*mesh_scale2), theta0((int(1.0E6*height)%360)*TO_RADIANS);
		unsigned const nlevels(unsigned(36.0*height*ms)), nrings(3);
		float rdeg(30.0);

		for (unsigned j = 0; j < nlevels; ++j) { // could do the same optimizations as the high detail pine tree
			for (unsigned k = 0; k < nrings; ++k) {
				float const sz(0.07*(height + 0.03/ms)*((nlevels - j + 3.0)/(float)nlevels));
				float const theta(TWO_PI*(3.3*j + k/5.0) + theta0);
				float const z((j + 3.0)*height/(nlevels + 4.0));
				int const val(int(((int(1.0E6*height))*(5463*j + 537879*k))%301));
				rdeg += 0.01*(val - 150);
				add_rotated_quad_pts(points, theta, rdeg/45.0, z, pos, vector3d(sz*wscale, sz*wscale, sz));
			}
		}
	}

	void draw(float sscale, int mode) { // modifies points, so non-const
		if (!sphere_in_camera_view(pos, (height + radius), 2)) return;
		int const light(get_light());
		float color_scale(SHADOW_VAL);

		for (unsigned i = 0; i < 3; ++i) {
			point p(pos);
			p.z += 0.5*i*height;

			if (is_visible_to_light_cobj(p, light, (radius + height), coll_id, 0)) {
				color_scale = 1.0;
				break;
			}
		}
		if (mode & 1) {
			colorRGBA color(pltype[type].stemc*color_scale);
			float const dist(distance_to_camera(pos));

			if (2*get_pt_line_thresh()*radius < dist) { // draw as line
				tree_scenery_pld.add_textured_line(pos, (pos + point(0.0, 0.0, 0.75*height)), color, WOOD_TEX);
			}
			else {
				int const ndiv(max(3, min(N_CYL_SIDES, int(5.0*sscale*radius/dist))));
				select_texture(WOOD_TEX);
				set_color(color);
				draw_fast_cylinder(pos, (pos + point(0.0, 0.0, height)), radius, 0.0, ndiv, 1);
			}
		}
		if (mode & 2) {
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.75);
			set_color(pltype[type].leafc*color_scale);
			select_texture(pltype[type].tid);
			glNormal3f(0.0, 0.0, 1.0);
			set_lighted_sides(2);
			if (points.empty()) gen_points();
			draw_quads_from_pts(points);
			glDisable(GL_ALPHA_TEST);
			set_lighted_sides(1);
		}
	}

	void remove_cobjs() {
		remove_reset_coll_obj(coll_id2);
		scenery_obj::remove_cobjs();
	}

	void destroy() {
		remove_cobjs();
		scenery_obj::destroy(); // will remove coll_id twice, which is OK
	}
};


// ************ SCENERY OBJECT INTERFACE/WRAPPERS/DRIVERS ************


vector<rock_shape3d> rock_shapes;
vector<surface_rock> surface_rocks;
vector<s_rock>       rocks;
vector<s_log>        logs;
vector<s_stump>      stumps;
vector<s_plant>      plants;


void clear_scenery_objs() {

	free_scenery();
	rock_shapes.clear();
	surface_rocks.clear();
	rocks.clear();
	logs.clear();
	stumps.clear();
	plants.clear();
}


void gen_scenery() {

	if (has_scenery2) return; // don't generate scenery if some has already been added
	clear_scenery_objs();
	has_scenery = 0;
	if (DISABLE_SCENERY || (NO_ISLAND_SCENERY && island)) return;
	has_scenery = 1;
	msms2       = mesh_scale*mesh_scale2;
	gen_scenery_deterministic();
	add_scenery_cobjs();
}


void add_plant(point const &pos, float height, float radius, int type, int calc_z) {

	assert(height > 0.0 && radius > 0.0);
	plants.push_back(s_plant());
	plants.back().create2(pos, height, radius, type, calc_z);
	has_scenery2 = 1;
}


void gen_scenery_deterministic() {

	unsigned const smod(max(200U, unsigned(3.321*XY_MULT_SIZE/(mesh_scale*mesh_scale2+1))));
	float const min_stump_z(water_plane_z + 0.010*zmax_est);
	float const min_plant_z(water_plane_z + 0.016*zmax_est);
	float const min_log_z  (water_plane_z - 0.040*zmax_est);
	clear_scenery_objs();

	for (int i = get_ext_y1(); i < get_ext_y2(); ++i) {
		for (int j = get_ext_x1(); j < get_ext_x2(); ++j) {
			rseed1 = 786433* (i + yoff2) + 196613 *rand_gen_index;
			rseed2 = 6291469*(j + xoff2) + 1572869*rand_gen_index;
			int const val((rand2() + rand2())%smod);
			if (val > 100) continue;
			bool const veg((rseed1&127)/128.0 < vegetation);
			
			if (veg && rand2()%100 < 30) {
				plants.push_back(s_plant());
				if (!plants.back().create(j, i, 1, min_plant_z)) plants.pop_back();
			}
			else if (val < 5) {
				rock_shapes.push_back(rock_shape3d());
				rock_shapes.back().create(j, i, 1);
			}
			else if (val < 15) {
				surface_rocks.push_back(surface_rock());
				surface_rocks.back().create(j, i, 1);
			}
			else if (val < 50) {
				rocks.push_back(s_rock());
				rocks.back().create(j, i, 1);
			}
			else if (veg && val < 85) {
				logs.push_back(s_log());
				if (!logs.back().create(j, i, 1, min_log_z)) logs.pop_back();
			}
			else if (veg) {
				stumps.push_back(s_stump());
				if (!stumps.back().create(j, i, 1, min_stump_z)) stumps.pop_back();
			}
		}
	}
	surface_rock_cache.clear_unref();
	sort(plants.begin(), plants.end()); // sort by type
}


template<typename T> void draw_scenery_vector(vector<T> &v, float sscale) {
	for (unsigned i = 0; i < v.size(); ++i) v[i].draw(sscale);
}

template<typename T> void add_scenery_vector_cobjs(vector<T> &v) {
	for (unsigned i = 0; i < v.size(); ++i) {
		if (is_over_mesh(v[i].get_pos())) v[i].add_cobjs();
	}
}

template<typename T> void shift_scenery_vector(vector<T> &v, vector3d const &vd) {
	for (unsigned i = 0; i < v.size(); ++i) v[i].shift_by(vd);
}

template<typename T> void free_scenery_vector(vector<T> &v) {
	for (unsigned i = 0; i < v.size(); ++i) v[i].destroy();
}

template<typename T> void update_scenery_zvals_vector(vector<T> &v, int x1, int y1, int x2, int y2) {
	
	for (unsigned i = 0; i < v.size(); ++i) { // zval has change, remove and re-add cobjs
		if (v[i].update_zvals(x1, y1, x2, y2)) {
			v[i].remove_cobjs();
			v[i].add_cobjs();
		}
	}
}


void draw_scenery(bool draw_opaque, bool draw_transparent) {

	if (!has_scenery && !has_scenery2) return;
	set_fill_mode();
	assert(quadric != 0);
	int const sscale(int((do_zoom ? ZOOM_FACTOR : 1.0)*window_width));

	if (draw_opaque) {
		for (unsigned i = 0; i < rock_shapes.size(); ++i) { // draw rock shapes
			rock_shapes[i].draw();
		}
		draw_scenery_vector(surface_rocks, sscale);
		draw_scenery_vector(rocks,  sscale); // can unset gluQuadricTexture
		gluQuadricTexture(quadric, GL_TRUE);
		draw_scenery_vector(logs,   sscale);
		draw_scenery_vector(stumps, sscale);
		gluQuadricTexture(quadric, GL_FALSE);

		for (unsigned i = 0; i < plants.size(); ++i) {
			plants[i].draw(sscale, 1); // draw stem
		}
		glDisable(GL_TEXTURE_2D);
		tree_scenery_pld.draw_and_clear();
	}
	if (draw_transparent) {
		enable_blend();

		for (unsigned i = 0; i < plants.size(); ++i) {
			plants[i].draw(sscale, 2); // draw leaves
		}
		disable_blend();
		glDisable(GL_TEXTURE_2D);
	}
	// *** draw grass? ***
}


void add_scenery_cobjs() {

	add_scenery_vector_cobjs(rock_shapes);
	add_scenery_vector_cobjs(surface_rocks);
	add_scenery_vector_cobjs(rocks);
	add_scenery_vector_cobjs(logs);
	add_scenery_vector_cobjs(stumps);
	add_scenery_vector_cobjs(plants);
}


void shift_scenery(vector3d const &vd) {

	if (!has_scenery2) return; // dynamically created, not placed
	shift_scenery_vector(rock_shapes,   vd);
	shift_scenery_vector(surface_rocks, vd);
	shift_scenery_vector(rocks,         vd);
	shift_scenery_vector(logs,          vd);
	shift_scenery_vector(stumps,        vd);
	shift_scenery_vector(plants,        vd);
}


// update region is inclusive: [x1,x2]x[y1,y2]
void update_scenery_zvals(int x1, int y1, int x2, int y2) { // inefficient, should use spatial subdivision

	assert(x1 <= x2 && y1 <= y2);
	// test if there are any cobjs within this region?
	update_scenery_zvals_vector(rock_shapes,   x1, y1, x2, y2);
	update_scenery_zvals_vector(surface_rocks, x1, y1, x2, y2);
	update_scenery_zvals_vector(rocks,         x1, y1, x2, y2);
	update_scenery_zvals_vector(logs,          x1, y1, x2, y2);
	update_scenery_zvals_vector(stumps,        x1, y1, x2, y2);
	update_scenery_zvals_vector(plants,        x1, y1, x2, y2);
}


void free_scenery() {

	free_scenery_vector(rock_shapes);
	free_scenery_vector(surface_rocks);
	free_scenery_vector(rocks);
	free_scenery_vector(logs);
	free_scenery_vector(stumps);
	free_scenery_vector(plants);
}



