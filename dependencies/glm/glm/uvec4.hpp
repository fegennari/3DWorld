/// @ref core
/// @file glm/uvec4.hpp

#pragma once
#include "detail/type_vec4.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_UINT))
		typedef vec<4, uint, lowp>		uvec4;
#	elif(defined(GLM_PRECISION_MEDIUMP_UINT))
		typedef vec<4, uint, mediump>	uvec4;
#	else //defined(GLM_PRECISION_HIGHP_UINT)
		/// 4 components vector of unsigned integer numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<4, uint, highp>		uvec4;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
