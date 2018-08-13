/// @ref ext_uvec1_precision
/// @file glm/ext/uvec1_precision.hpp
///
/// @see core (dependence)
///
/// @defgroup ext_uvec1_precision GLM_EXT_uvec1_precision
/// @ingroup ext
///
/// Include <glm/ext/uvec1_precision.hpp> to use the features of this extension.
///
/// Add Exposes highp_uvec1, mediump_uvec1 and lowp_uvec1 types.

#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_uvec1_precision extension included")
#endif

namespace glm
{
	/// @addtogroup ext_uvec1_precision
	/// @{

	/// 1 component vector of unsigned integer values.
	///
	/// @see ext_vec1
	typedef vec<1, uint, highp>			highp_uvec1;

	/// 1 component vector of unsigned integer values.
	///
	/// @see ext_vec1
	typedef vec<1, uint, mediump>		mediump_uvec1;

	/// 1 component vector of unsigned integer values.
	///
	/// @see ext_vec1
	typedef vec<1, uint, lowp>			lowp_uvec1;

	/// @}
}//namespace glm
