/// @ref core
/// @file glm/uvec3.hpp

#pragma once
#include "detail/type_vec3.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_UINT))
		typedef vec<3, uint, lowp>		uvec3;
#	elif(defined(GLM_PRECISION_MEDIUMP_UINT))
		typedef vec<3, uint, mediump>	uvec3;
#	else //defined(GLM_PRECISION_HIGHP_UINT)
		/// 3 components vector of unsigned integer numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<3, uint, highp>		uvec3;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
