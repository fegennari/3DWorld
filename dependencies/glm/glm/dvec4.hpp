/// @ref core
/// @file glm/dvec4.hpp

#pragma once
#include "detail/type_vec4.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_DOUBLE))
		typedef vec<4, double, lowp>		dvec4;
#	elif(defined(GLM_PRECISION_MEDIUMP_DOUBLE))
		typedef vec<4, double, mediump>		dvec4;
#	else //defined(GLM_PRECISION_HIGHP_DOUBLE)
		/// 4 components vector of double-precision floating-point numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<4, double, highp>		dvec4;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
