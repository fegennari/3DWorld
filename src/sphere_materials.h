// 3D World - Throwable Sphere Materials header
// by Frank Gennari
// 9/9/16
#include "3DWorld.h"

#pragma once

struct sphere_mat_t {
	bool shadows, emissive, reflective;
	int destroy_thresh;
	float alpha, metal, spec_mag, shine, hardness, density, light_atten, refract_ix, light_radius;
	colorRGB diff_c, spec_c;
	std::string name;

	sphere_mat_t() : shadows(0), emissive(0), reflective(0), destroy_thresh(0), alpha(1.0), metal(1.0), spec_mag(0.0), shine(1.0),
		hardness(0.8), density(1.0), light_atten(0.0), refract_ix(1.0), light_radius(0.0), diff_c(WHITE), spec_c(WHITE) {}
};

sphere_mat_t &get_cur_sphere_mat();

