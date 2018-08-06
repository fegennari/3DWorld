/// @ref core
/// @file glm/vec3.hpp

#pragma once
#include "bvec3.hpp"
#include "bvec3_precision.hpp"
#include "fvec3.hpp"
#include "fvec3_precision.hpp"
#include "dvec3.hpp"
#include "dvec3_precision.hpp"
#include "ivec3.hpp"
#include "ivec3_precision.hpp"
#include "uvec3.hpp"
#include "uvec3_precision.hpp"

namespace glm
{
#	if GLM_HAS_TEMPLATE_ALIASES
		template <typename T, qualifier Q = defaultp> using tvec3 = vec<3, T, Q>;
#	endif
}//namespace glm

