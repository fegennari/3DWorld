/// @ref ext_dvec1
/// @file glm/ext/dvec1.hpp
///
/// @see core (dependence)
///
/// @defgroup ext_dvec1 GLM_EXT_dvec1
/// @ingroup ext
///
/// Include <glm/ext/dvec1.hpp> to use the features of this extension.
///
/// Expose dvec1 vector type.

#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_dvec1 extension included")
#endif

namespace glm
{
	/// @addtogroup ext_dvec1
	/// @{

#	if(defined(GLM_PRECISION_LOWP_DOUBLE))
		typedef vec<1, double, lowp>		dvec1;
#	elif(defined(GLM_PRECISION_MEDIUMP_DOUBLE))
		typedef vec<1, double, mediump>		dvec1;
#	else //defined(GLM_PRECISION_HIGHP_DOUBLE)
		/// 1 components vector of double-precision floating-point numbers.
		///
		/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
		typedef vec<1, double, highp>		dvec1;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
