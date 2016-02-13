// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/27/02

#include "3DWorld.h"
#include "textures_3dw.h"
#include "shape_line3d.h"
#include "collision_detect.h"


float const wheelr = 0.118;

float hmv_scale(1.0);
point hmv_pos(all_zeros);
shape3d hmv_shape;
std::vector<int> hmv_coll_obj;

extern int load_hmv;


void add_hmv_coll_objs(point &pos, float scale);


void build_hmv_shape() {

	if (!load_hmv) return;
	delete_hmv_shape();
	hmv_shape.set_shape_color(WHITE);
	cout << "Reading mesh hmv.mesh." << endl;
	if (!hmv_shape.read_from_file("hmv.mesh")) {load_hmv = 0;}
	else {hmv_shape.gen_face_normals();}
	hmv_shape.set_scale(hmv_scale);
	hmv_shape.move_to(hmv_pos);
}


void delete_hmv_shape() {hmv_shape.destroy();}

// well, not actually a draw function...
void add_shape_coll_objs() {
	if (load_hmv) add_hmv_coll_objs(hmv_pos, hmv_scale);
}


void add_hmv_coll_objs(point &pos, float scale) {

	float x(pos.x), y(pos.y), z(pos.z);

	if (hmv_coll_obj.empty()) {
		hmv_coll_obj.reserve(hmv_shape.get_num_faces() + 4);
	}
	else {
		for (unsigned i = 0; i < hmv_coll_obj.size(); ++i) {remove_reset_coll_obj(hmv_coll_obj[i]);}
		purge_coll_freed(1);
	}
	// wheels
	cobj_params const cp(0.9, BKGRAY, 1, 0, nullptr, 0, -1, 1.0, 0, 0.5, 20.0); // drawn
	point const pts[4] = {point(x+0.19*scale, y+0.04*scale, z), point(x+0.19*scale, y+0.46*scale, z), point(x+0.81*scale, y+0.04*scale, z), point(x+0.81*scale, y+0.46*scale, z)};
	for (unsigned i = 0; i < 4; ++i) {hmv_coll_obj.push_back(add_coll_torus(pts[i], plus_y, 0.6*wheelr, 0.4*wheelr, cp));}
	// body
	hmv_shape.add_cobjs(hmv_coll_obj, 1); // drawn
}


void shift_hmv(vector3d const &vd) {
	hmv_pos += vd;
	hmv_shape.translate(vd);
}



