/// @ref core
/// @file glm/dvec3.hpp

#pragma once
#include "detail/type_vec3.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_DOUBLE))
		typedef vec<3, double, lowp>		dvec3;
#	elif(defined(GLM_PRECISION_MEDIUMP_DOUBLE))
		typedef vec<3, double, mediump>		dvec3;
#	else //defined(GLM_PRECISION_HIGHP_DOUBLE)
		/// 3 components vector of double-precision floating-point numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<3, double, highp>		dvec3;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
