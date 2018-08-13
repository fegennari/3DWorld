/// @ref core
/// @file glm/uvec2.hpp

#pragma once
#include "detail/type_vec2.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_UINT))
		typedef vec<2, uint, lowp>		uvec2;
#	elif(defined(GLM_PRECISION_MEDIUMP_UINT))
		typedef vec<2, uint, mediump>	uvec2;
#	else //defined(GLM_PRECISION_HIGHP_UINT)
		/// 2 components vector of unsigned integer numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<2, uint, highp>		uvec2;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
