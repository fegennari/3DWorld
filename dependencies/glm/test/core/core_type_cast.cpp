///////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL Mathematics Copyright (c) 2005 - 2014 G-Truc Creation (www.g-truc.net)
///////////////////////////////////////////////////////////////////////////////////////////////////
// Created : 2013-05-06
// Updated : 2013-05-06
// Licence : This source is under MIT License
// File    : test/core/type_cast.cpp
///////////////////////////////////////////////////////////////////////////////////////////////////

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>

struct my_vec2
{
	operator glm::vec2() { return glm::vec2(x, y); }
	float x, y;
};

int test_vec2_cast()
{
	glm::vec2 A(1.0f, 2.0f);
	glm::lowp_vec2 B(A);
	glm::mediump_vec2 C(A);
	glm::highp_vec2 D(A);
	
	glm::vec2 E = static_cast<glm::vec2>(A);
	glm::lowp_vec2 F = static_cast<glm::lowp_vec2>(A);
	glm::mediump_vec2 G = static_cast<glm::mediump_vec2>(A);
	glm::highp_vec2 H = static_cast<glm::highp_vec2>(A);
	
	my_vec2 I;
	glm::vec2 J = static_cast<glm::vec2>(I);
	glm::vec2 K(7.8f);

	int Error(0);
	
	Error += glm::all(glm::equal(A, E)) ? 0 : 1;
	Error += glm::all(glm::equal(B, F)) ? 0 : 1;
	Error += glm::all(glm::equal(C, G)) ? 0 : 1;
	Error += glm::all(glm::equal(D, H)) ? 0 : 1;
	
	return Error;
}

int test_vec3_cast()
{
	glm::vec3 A(1.0f, 2.0f, 3.0f);
	glm::lowp_vec3 B(A);
	glm::mediump_vec3 C(A);
	glm::highp_vec3 D(A);
	
	glm::vec3 E = static_cast<glm::vec3>(A);
	glm::lowp_vec3 F = static_cast<glm::lowp_vec3>(A);
	glm::mediump_vec3 G = static_cast<glm::mediump_vec3>(A);
	glm::highp_vec3 H = static_cast<glm::highp_vec3>(A);
	
	int Error(0);
	
	Error += glm::all(glm::equal(A, E)) ? 0 : 1;
	Error += glm::all(glm::equal(B, F)) ? 0 : 1;
	Error += glm::all(glm::equal(C, G)) ? 0 : 1;
	Error += glm::all(glm::equal(D, H)) ? 0 : 1;
	
	return Error;
}

int test_vec4_cast()
{
	glm::vec4 A(1.0f, 2.0f, 3.0f, 4.0f);
	glm::lowp_vec4 B(A);
	glm::mediump_vec4 C(A);
	glm::highp_vec4 D(A);
	
	glm::vec4 E = static_cast<glm::vec4>(A);
	glm::lowp_vec4 F = static_cast<glm::lowp_vec4>(A);
	glm::mediump_vec4 G = static_cast<glm::mediump_vec4>(A);
	glm::highp_vec4 H = static_cast<glm::highp_vec4>(A);
	
	int Error(0);
	
	Error += glm::all(glm::equal(A, E)) ? 0 : 1;
	Error += glm::all(glm::equal(B, F)) ? 0 : 1;
	Error += glm::all(glm::equal(C, G)) ? 0 : 1;
	Error += glm::all(glm::equal(D, H)) ? 0 : 1;
	
	return Error;
}

int main()
{
	int Error = 0;

	Error += test_vec2_cast();
	Error += test_vec3_cast();
	Error += test_vec4_cast();

	return Error;
}
