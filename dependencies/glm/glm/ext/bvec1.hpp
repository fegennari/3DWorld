/// @ref ext_bvec1
/// @file glm/ext/bvec1.hpp
///
/// @see core (dependence)
///
/// @defgroup ext_bvec1 GLM_EXT_bvec1
/// @ingroup ext
///
/// Include <glm/ext/bvec1.hpp> to use the features of this extension.
///
/// Expose bvec1 type.

#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_bvec1 extension included")
#endif

namespace glm
{
	/// @addtogroup ext_bvec1
	/// @{

#	if(defined(GLM_PRECISION_LOWP_BOOL))
		typedef vec<1, bool, lowp>		bvec1;
#	elif(defined(GLM_PRECISION_MEDIUMP_BOOL))
		typedef vec<1, bool, mediump>	bvec1;
#	else //defined(GLM_PRECISION_HIGHP_BOOL)
		/// 1 components vector of boolean.
		///
		/// @see ext_bvec1 extension.
		typedef vec<1, bool, highp>		bvec1;
#	endif//GLM_PRECISION

	/// @}
}//namespace glm
