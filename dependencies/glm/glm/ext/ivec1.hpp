/// @ref ext_ivec1
/// @file glm/ext/ivec1.hpp
///
/// @see core (dependence)
///
/// @defgroup ext_ivec1 GLM_EXT_ivec1
/// @ingroup ext
///
/// Include <glm/ext/ivec1.hpp> to use the features of this extension.
///
/// Exposes ivec1 vector type.

#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_ivec1 extension included")
#endif

namespace glm
{
	/// @addtogroup ext_ivec1
	/// @{

#	if(defined(GLM_PRECISION_LOWP_INT))
		typedef vec<1, int, lowp>			ivec1;
#	elif(defined(GLM_PRECISION_MEDIUMP_INT))
		typedef vec<1, int, mediump>		ivec1;
	#else //defined(GLM_PRECISION_HIGHP_INT)
		/// 1 component vector of signed integer numbers.
		///
		/// @see ext_ivec1 extension.
		typedef vec<1, int, highp>			ivec1;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm

