/// @ref core
/// @file glm/bvec3.hpp

#pragma once
#include "detail/type_vec3.hpp"

namespace glm
{
	/// @addtogroup core
	/// @{

#	if(defined(GLM_PRECISION_LOWP_BOOL))
		typedef vec<3, bool, lowp>		bvec3;
#	elif(defined(GLM_PRECISION_MEDIUMP_BOOL))
		typedef vec<3, bool, mediump>	bvec3;
#	else //defined(GLM_PRECISION_HIGHP_BOOL)
		/// 3 components vector of boolean.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<3, bool, highp>		bvec3;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
