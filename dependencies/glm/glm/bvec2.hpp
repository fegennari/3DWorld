/// @ref core
/// @file glm/bvec2.hpp

#pragma once
#include "detail/type_vec2.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_BOOL))
		typedef vec<2, bool, lowp>		bvec2;
#	elif(defined(GLM_PRECISION_MEDIUMP_BOOL))
		typedef vec<2, bool, mediump>	bvec2;
#	else //defined(GLM_PRECISION_HIGHP_BOOL)
		/// 2 components vector of boolean.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<2, bool, highp>		bvec2;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
