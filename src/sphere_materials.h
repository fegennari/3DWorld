// 3D World - Throwable Sphere Materials header
// by Frank Gennari
// 9/9/16
#include "3DWorld.h"

#pragma once

struct sphere_mat_t {
	bool shadows=0, emissive=0, reflective=0;
	int destroyable=0, tid=-1, nm_tid=-1;
	float radius_scale=1.0, alpha=1.0, metal=1.0, spec_mag=0.0, shine=1.0, hardness=0.8, density=1.0, light_atten=0.0, refract_ix=1.0, light_radius=0.0;
	colorRGB diff_c=WHITE, spec_c=WHITE;
	std::string name;

	std::string get_name() const;
};

sphere_mat_t &get_cur_sphere_mat();

