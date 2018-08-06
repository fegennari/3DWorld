/// @ref core
/// @file glm/vec4.hpp

#pragma once
#include "bvec4.hpp"
#include "bvec4_precision.hpp"
#include "fvec4.hpp"
#include "fvec4_precision.hpp"
#include "dvec4.hpp"
#include "dvec4_precision.hpp"
#include "ivec4.hpp"
#include "ivec4_precision.hpp"
#include "uvec4.hpp"
#include "uvec4_precision.hpp"

namespace glm
{
#	if GLM_HAS_TEMPLATE_ALIASES
		template <typename T, qualifier Q = defaultp> using tvec4 = vec<4, T, Q>;
#	endif
}//namespace glm

