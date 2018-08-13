/// @ref core
/// @file glm/ivec4.hpp

#pragma once
#include "detail/type_vec4.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_INT))
		typedef vec<4, int, lowp>		ivec4;
#	elif(defined(GLM_PRECISION_MEDIUMP_INT))
		typedef vec<4, int, mediump>	ivec4;
#	else //defined(GLM_PRECISION_HIGHP_INT)
		/// 4 components vector of signed integer numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<4, int, highp>		ivec4;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
