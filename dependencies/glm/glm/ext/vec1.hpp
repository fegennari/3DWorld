/// @ref core
/// @file glm/vec1.hpp

#pragma once
#include "bvec1.hpp"
#include "bvec1_precision.hpp"
#include "fvec1.hpp"
#include "fvec1_precision.hpp"
#include "dvec1.hpp"
#include "dvec1_precision.hpp"
#include "ivec1.hpp"
#include "ivec1_precision.hpp"
#include "uvec1.hpp"
#include "uvec1_precision.hpp"

namespace glm
{
#	if GLM_HAS_TEMPLATE_ALIASES
	template <typename T, qualifier Q = defaultp> using tvec1 = vec<1, T, Q>;
#	endif
}//namespace glm
