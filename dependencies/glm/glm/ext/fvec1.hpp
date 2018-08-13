/// @ref ext_vec1
/// @file glm/ext/vec1.hpp
///
/// @see core (dependence)
///
/// @defgroup ext_vec1 GLM_EXT_fvec1
/// @ingroup ext
///
/// Include <glm/ext/vec1.hpp> to use the features of this extension.
///
/// Expose vec1 vector type.

#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_fvec1 extension included")
#endif

namespace glm
{
	/// @addtogroup ext_fvec1
	/// @{

#	if(defined(GLM_PRECISION_LOWP_FLOAT))
		typedef vec<1, float, lowp>			vec1;
#	elif(defined(GLM_PRECISION_MEDIUMP_FLOAT))
		typedef vec<1, float, mediump>		vec1;
#	else //defined(GLM_PRECISION_HIGHP_FLOAT)
		/// 1 components vector of single-precision floating-point numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<1, float, highp>		vec1;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
