///////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL Mathematics Copyright (c) 2005 - 2014 G-Truc Creation (www.g-truc.net)
///////////////////////////////////////////////////////////////////////////////////////////////////
// Created : 2013-11-22
// Updated : 2013-11-22
// Licence : This source is under MIT licence
// File    : test/gtx/io.cpp
///////////////////////////////////////////////////////////////////////////////////////////////////

#define GLM_FORCE_RADIANS
#include <glm/gtc/type_precision.hpp>
#include <glm/gtx/io.hpp>
#include <iostream>
#include <typeinfo>

namespace {
  
  template <typename CTy, typename CTr>
  std::basic_ostream<CTy,CTr>&
  operator<<(std::basic_ostream<CTy,CTr>& os, glm::precision const& a)
  {
    typename std::basic_ostream<CTy,CTr>::sentry const cerberus(os);

    if (cerberus) {
      switch (a) {
      case glm::highp:   os << "hi"; break;
      case glm::mediump: os << "md"; break;
      case glm::lowp:    os << "lo"; break;
      }
    }

    return os;
  }

} // namespace {

template <typename T, glm::precision P, typename OS>
int test_io_vec(OS& os)
{
  os << '\n'
     << typeid(OS).name()
     << '\n';
  
  glm::detail::tvec2<T,P> const v2(0, 1);
  glm::detail::tvec3<T,P> const v3(2, 3, 4);
  glm::detail::tvec4<T,P> const v4(5, 6, 7, 8);

  os << "vec2<" << typeid(T).name() << ',' << P << ">: " << v2 << '\n'
     << "vec3<" << typeid(T).name() << ',' << P << ">: " << v3 << '\n'
     << "vec4<" << typeid(T).name() << ',' << P << ">: " << v4 << '\n';

  glm::io::precision_guard const iopg;
  
  glm::io::precision()   = 2;
  glm::io::value_width() = 1 + 2 + 1 + glm::io::precision();
  
  os << "vec2<" << typeid(T).name() << ',' << P << ">: " << v2 << '\n'
     << "vec3<" << typeid(T).name() << ',' << P << ">: " << v3 << '\n'
     << "vec4<" << typeid(T).name() << ',' << P << ">: " << v4 << '\n';

  return 0;
}

template <typename T, glm::precision P, typename OS>
int test_io_mat(OS& os)
{
  os << '\n'
     << typeid(OS).name()
     << '\n';

  glm::detail::tvec2<T,P> const v2_1( 0,  1);
  glm::detail::tvec2<T,P> const v2_2( 2,  3);
  glm::detail::tvec2<T,P> const v2_3( 4,  5);
  glm::detail::tvec2<T,P> const v2_4( 6,  7);
  glm::detail::tvec3<T,P> const v3_1( 8,  9, 10);
  glm::detail::tvec3<T,P> const v3_2(11, 12, 13);
  glm::detail::tvec3<T,P> const v3_3(14, 15, 16);
  glm::detail::tvec3<T,P> const v3_4(17, 18, 19);
  glm::detail::tvec4<T,P> const v4_1(20, 21, 22, 23);
  glm::detail::tvec4<T,P> const v4_2(24, 25, 26, 27);
  glm::detail::tvec4<T,P> const v4_3(28, 29, 30, 31);
  glm::detail::tvec4<T,P> const v4_4(32, 33, 34, 35);

#if 0
  os << "mat2x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x2<T,P>(v2_1, v2_2) << '\n'
     << "mat2x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x3<T,P>(v3_1, v3_2) << '\n'
     << "mat2x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x4<T,P>(v4_1, v4_2) << '\n'
     << "mat3x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x2<T,P>(v2_1, v2_2, v2_3) << '\n'
     << "mat3x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x3<T,P>(v3_1, v3_2, v3_3) << '\n'
     << "mat3x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x4<T,P>(v4_1, v4_2, v4_3) << '\n'
     << "mat4x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x2<T,P>(v2_1, v2_2, v2_3, v2_4) << '\n'
     << "mat4x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x3<T,P>(v3_1, v3_2, v3_3, v3_4) << '\n'
     << "mat4x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x4<T,P>(v4_1, v4_2, v4_3, v4_4) << '\n';
#endif
  
  glm::io::precision_guard const iopg;
  
  glm::io::precision()   = 2;
  glm::io::value_width() = 1 + 2 + 1 + glm::io::precision();
  
  os << "mat2x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x2<T,P>(v2_1, v2_2) << '\n'
     << "mat2x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x3<T,P>(v3_1, v3_2) << '\n'
     << "mat2x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x4<T,P>(v4_1, v4_2) << '\n'
     << "mat3x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x2<T,P>(v2_1, v2_2, v2_3) << '\n'
     << "mat3x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x3<T,P>(v3_1, v3_2, v3_3) << '\n'
     << "mat3x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x4<T,P>(v4_1, v4_2, v4_3) << '\n'
     << "mat4x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x2<T,P>(v2_1, v2_2, v2_3, v2_4) << '\n'
     << "mat4x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x3<T,P>(v3_1, v3_2, v3_3, v3_4) << '\n'
     << "mat4x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x4<T,P>(v4_1, v4_2, v4_3, v4_4) << '\n';
  
  os << glm::io::column_major
     << "mat2x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x2<T,P>(v2_1, v2_2) << '\n'
     << "mat2x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x3<T,P>(v3_1, v3_2) << '\n'
     << "mat2x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat2x4<T,P>(v4_1, v4_2) << '\n'
     << "mat3x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x2<T,P>(v2_1, v2_2, v2_3) << '\n'
     << "mat3x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x3<T,P>(v3_1, v3_2, v3_3) << '\n'
     << "mat3x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat3x4<T,P>(v4_1, v4_2, v4_3) << '\n'
     << "mat4x2<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x2<T,P>(v2_1, v2_2, v2_3, v2_4) << '\n'
     << "mat4x3<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x3<T,P>(v3_1, v3_2, v3_3, v3_4) << '\n'
     << "mat4x4<" << typeid(T).name() << ',' << P << ">: " << glm::detail::tmat4x4<T,P>(v4_1, v4_2, v4_3, v4_4) << '\n';
  
  return 0;
}

int main()
{
	int Error(0);

	Error += test_io_vec<float, glm::highp>(std::cout);
	Error += test_io_vec<float, glm::highp>(std::wcout);
	Error += test_io_vec<int, glm::mediump>(std::cout);
	Error += test_io_vec<int, glm::mediump>(std::wcout);
	Error += test_io_vec<glm::uint, glm::lowp>(std::cout);
	Error += test_io_vec<glm::uint, glm::lowp>(std::wcout);

	Error += test_io_mat<float, glm::highp>(std::cout);
	Error += test_io_mat<float, glm::lowp>(std::wcout);

	return Error;
}
