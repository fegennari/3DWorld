/// @ref core
/// @file glm/ivec3.hpp

#pragma once
#include "detail/type_vec3.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_INT))
		typedef vec<3, int, lowp>		ivec3;
#	elif(defined(GLM_PRECISION_MEDIUMP_INT))
		typedef vec<3, int, mediump>	ivec3;
#	else //defined(GLM_PRECISION_HIGHP_INT)
		/// 3 components vector of signed integer numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<3, int, highp>		ivec3;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
