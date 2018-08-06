/// @ref core
/// @file glm/vec2.hpp

#pragma once
#include "bvec2.hpp"
#include "bvec2_precision.hpp"
#include "fvec2.hpp"
#include "fvec2_precision.hpp"
#include "dvec2.hpp"
#include "dvec2_precision.hpp"
#include "ivec2.hpp"
#include "ivec2_precision.hpp"
#include "uvec2.hpp"
#include "uvec2_precision.hpp"

namespace glm
{
#	if GLM_HAS_TEMPLATE_ALIASES
		template <typename T, qualifier Q = defaultp> using tvec2 = vec<2, T, Q>;
#	endif
}//namespace glm
