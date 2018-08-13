/// @ref core
/// @file glm/fvec3.hpp

#pragma once
#include "detail/type_vec3.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_FLOAT))
		typedef vec<3, float, lowp>		vec3;
#	elif(defined(GLM_PRECISION_MEDIUMP_FLOAT))
		typedef vec<3, float, mediump>	vec3;
#	else //defined(GLM_PRECISION_HIGHP_FLOAT)
		/// 3 components vector of single-precision floating-point numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<3, float, highp>	vec3;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
