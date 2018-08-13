/// @ref core
/// @file glm/ivec2.hpp

#pragma once
#include "detail/type_vec2.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_INT))
		typedef vec<2, int, lowp>		ivec2;
#	elif(defined(GLM_PRECISION_MEDIUMP_INT))
		typedef vec<2, int, mediump>	ivec2;
#	else //defined(GLM_PRECISION_HIGHP_INT)
		/// 2 components vector of signed integer numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<2, int, highp>		ivec2;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
