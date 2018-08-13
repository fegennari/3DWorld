/// @ref core
/// @file glm/fvec2.hpp

#pragma once
#include "detail/type_vec2.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_FLOAT))
		typedef vec<2, float, lowp>		vec2;
#	elif(defined(GLM_PRECISION_MEDIUMP_FLOAT))
		typedef vec<2, float, mediump>	vec2;
#	else //defined(GLM_PRECISION_HIGHP_FLOAT)
		/// 2 components vector of single-precision floating-point numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<2, float, highp>	vec2;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
