// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/27/02

#include "3DWorld.h"
#include "textures_3dw.h"
#include "shape_line3d.h"
#include "collision_detect.h"


int const DRAW_HMV_COBJS = 1;
float const wheelr       = 0.118;


float hmv_scale(1.0);
point hmv_pos(all_zeros);
shape3d hmv_shape;
std::vector<int> hmv_coll_obj;

extern int load_hmv;
extern GLUquadricObj* quadric;


void add_hmv_coll_objs(point &pos, float scale);



void draw_hmv() {

	if (DRAW_HMV_COBJS || !load_hmv) return;
	set_fill_mode();

	// draw body
	hmv_shape.draw(); // CAMOFLAGE_TEX

	// draw wheels
	glPushMatrix();
	translate_to(hmv_pos);
	uniform_scale(hmv_scale);
	set_color(BKGRAY);
	glRotatef(90.0, 1.0, 0.0, 0.0); // x = x, y = z, z = -y
	unsigned const nsides((3*N_CYL_SIDES)/2);

	for (unsigned i = 0; i < 4; ++i) {
		point const pos((0.19 + (float(i>1))*.62), 0.0, (-(float(i&1))*0.42 - 0.08));
		draw_fast_cylinder(pos, pos+vector3d(0.0, 0.0, 0.08), wheelr, wheelr, nsides, 0, 1);
	}
	glPopMatrix();
}


void build_hmv_shape() {

	if (!load_hmv) return;
	delete_hmv_shape();
	hmv_shape.set_shape_color(WHITE);
	cout << "Reading mesh hmv.mesh." << endl;

	if (!hmv_shape.read_from_file("hmv.mesh")) {
		load_hmv = 0;
	}
	else {
		hmv_shape.gen_face_normals();
	}
	hmv_shape.set_scale(hmv_scale);
	hmv_shape.move_to(hmv_pos);
}


void delete_hmv_shape() {

	hmv_shape.destroy();
}


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
		for (unsigned i = 0; i < hmv_coll_obj.size(); ++i) {
			remove_reset_coll_obj(hmv_coll_obj[i]);
		}
		purge_coll_freed(1);
	}
	// wheels
	cobj_params const cp(0.9, BKGRAY, DRAW_HMV_COBJS, 0);
	hmv_coll_obj.push_back(add_coll_cylinder(x+0.19*scale, y+0.08*scale, z, x+0.19*scale, y,            z, wheelr, wheelr, cp));
	hmv_coll_obj.push_back(add_coll_cylinder(x+0.19*scale, y+0.5*scale,  z, x+0.19*scale, y+0.42*scale, z, wheelr, wheelr, cp));
	hmv_coll_obj.push_back(add_coll_cylinder(x+0.81*scale, y+0.08*scale, z, x+0.81*scale, y,            z, wheelr, wheelr, cp));
	hmv_coll_obj.push_back(add_coll_cylinder(x+0.81*scale, y+0.5*scale,  z, x+0.81*scale, y+0.42*scale, z, wheelr, wheelr, cp));

	// body
	hmv_shape.add_cobjs(hmv_coll_obj, DRAW_HMV_COBJS);
}


void shift_hmv(vector3d const &vd) {

	hmv_pos += vd;
	hmv_shape.translate(vd);
}



