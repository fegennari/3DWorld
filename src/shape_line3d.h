// 3D World - shape3d, line3d, and scenery objects
// by Frank Gennari
// 4/5/07

#ifndef _SHAPE_LINE3D_H_
#define _SHAPE_LINE3D_H_


#include "3DWorld.h"


class scenery_obj { // size = 28

protected:
	int coll_id, type;
	float radius, dz;
	point pos;

public:
	scenery_obj() : coll_id(-1), type(-1), radius(0.0), dz(0.0), pos(all_zeros) {}
	void shift_by(vector3d const &vd);
	void gen_spos(int x, int y, int use_xy);
	bool update_zvals(int x1, int y1, int x2, int y2);
	bool in_camera_view(float brad=0.0) const;
	float get_shadowed_color(point const &p, float eff_radius) const;
	void set_atten_color(colorRGBA c) const;
	void remove_cobjs();
	void destroy() {remove_cobjs();}
	point get_pos() const {return pos;}
};


struct scolor { // size = 28

	int tid;
	float spec1, spec2;
	colorRGBA c;
};


struct face3d { // size = 28

	unsigned v[3], color_id;
	vector3d norm;
};


class shape3d : public scenery_obj { // size = 72

protected:
	int tid;
	float tex_scale, scale;
	colorRGBA color;
	vector<point>  points;
	vector<face3d> faces;
	vector<scolor> colors;

public:
	shape3d() : tid(0), tex_scale(1.0), scale(1.0) {}
	bool alloc_shape(unsigned npoints, unsigned nfaces, unsigned ncolors);
	bool read_from_file(char *filename);
	void set_texture(int tid_, float tscale=1.0) {tid = tid_; tex_scale = tscale;}
	void move_to(point const &pos_)          {pos   = pos_;}
	void set_scale(float scale_)             {scale = scale_;}
	void set_shape_color(colorRGBA const &c) {color = c;}
	void translate(vector3d const &vd)       {shift_by(vd);}
	unsigned get_num_faces() const           {return faces.size();}
	void gen_face_normals();
	void get_face_normal(unsigned face_id);
	void get_triangle_center(point &center, unsigned face_id, unsigned quality);
	void add_vertex(unsigned vertex, unsigned face_id, unsigned &face_counter);
	void draw(bool skip_color_set=0) const;
	void add_cobjs(vector<int> &cids, bool draw);
	void destroy();
	bool is_allocated() const {return (!points.empty() && !faces.empty());}
	colorRGBA get_color() const {return color;}
};


class rock_shape3d : public shape3d { // size = 72

public:
	void create(int x, int y, bool use_xy);
	void gen_rock(unsigned nverts, float size, int &rand_seed, int type);
	void add_cobjs();
	bool do_impact_damage(point const &pos_, float radius_);
	void draw() const;
};


#endif // _SHAPE_LINE3D_H_

