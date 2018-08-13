/// @ref core
/// @file glm/dvec2.hpp

#pragma once
#include "detail/type_vec2.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_DOUBLE))
		typedef vec<2, double, lowp>		dvec2;
#	elif(defined(GLM_PRECISION_MEDIUMP_DOUBLE))
		typedef vec<2, double, mediump>		dvec2;
#	else //defined(GLM_PRECISION_HIGHP_DOUBLE)
		/// 2 components vector of double-precision floating-point numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<2, double, highp>		dvec2;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
