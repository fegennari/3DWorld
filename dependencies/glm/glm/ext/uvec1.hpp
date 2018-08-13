/// @ref ext_uvec1
/// @file glm/ext/uvec1.hpp
///
/// @see core (dependence)
///
/// @defgroup ext_uvec1 GLM_EXT_uvec1
/// @ingroup ext
///
/// Include <glm/ext/uvec1.hpp> to use the features of this extension.
///
/// Exposes uvec1 vector type.

#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_uvec1 extension included")
#endif

namespace glm
{
	/// @addtogroup ext_uvec1
	/// @{

#	if(defined(GLM_PRECISION_LOWP_INT))
	typedef vec<1, uint, lowp>			uvec1;
#	elif(defined(GLM_PRECISION_MEDIUMP_INT))
	typedef vec<1, uint, mediump>		uvec1;
#else //defined(GLM_PRECISION_HIGHP_INT)
	/// 1 component vector of unsigned integer numbers.
	///
	/// @see ext_uvec1 extension.
	typedef vec<1, uint, highp>			uvec1;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm

