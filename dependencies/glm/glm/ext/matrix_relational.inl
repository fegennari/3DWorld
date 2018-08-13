/// @ref ext_vector_relational
/// @file glm/ext/vector_relational.inl

// Dependency:
#include "../vector_relational.hpp"
#include "../common.hpp"
#include "../detail/type_vec.hpp"

namespace glm
{
	template<typename genType>
	GLM_FUNC_QUALIFIER bool equal(genType const& x, genType const& y, genType const& epsilon)
	{
		return abs(x - y) <= epsilon;
	}

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, bool, Q> equal(vec<L, T, Q> const& x, vec<L, T, Q> const& y, T epsilon)
	{
		return equal(x, y, vec<L, T, Q>(epsilon));
	}

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, bool, Q> equal(vec<L, T, Q> const& x, vec<L, T, Q> const& y, vec<L, T, Q> const& epsilon)
	{
		return lessThanEqual(abs(x - y), epsilon);
	}

	template<length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<C, bool, Q> equal(mat<C, R, T, Q> const& a, mat<C, R, T, Q> const& b)
	{
		return equal(a, b, static_cast<T>(0));
	}

	template<length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<C, bool, Q> equal(mat<C, R, T, Q> const& a, mat<C, R, T, Q> const& b, T epsilon)
	{
		return equal(a, b, vec<C, T, Q>(epsilon));
	}

	template<length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<C, bool, Q> equal(mat<C, R, T, Q> const& a, mat<C, R, T, Q> const& b, vec<C, T, Q> const& epsilon)
	{
		vec<C, bool, Q> Result;
		for(length_t i = 0, n = C; i < n; ++i)
			Result[i] = all(equal(a[i], b[i], epsilon[i]));
		return Result;
	}

	template<typename genType>
	GLM_FUNC_QUALIFIER bool notEqual(genType const& x, genType const& y, genType const& epsilon)
	{
		return abs(x - y) > epsilon;
	}

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, bool, Q> notEqual(vec<L, T, Q> const& x, vec<L, T, Q> const& y, T epsilon)
	{
		return notEqual(x, y, vec<L, T, Q>(epsilon));
	}

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, bool, Q> notEqual(vec<L, T, Q> const& x, vec<L, T, Q> const& y, vec<L, T, Q> const& epsilon)
	{
		return greaterThan(abs(x - y), epsilon);
	}

	template<length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<C, bool, Q> notEqual(mat<C, R, T, Q> const& x, mat<C, R, T, Q> const& y)
	{
		return notEqual(x, y, static_cast<T>(0));
	}

	template<length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<C, bool, Q> notEqual(mat<C, R, T, Q> const& x, mat<C, R, T, Q> const& y, T epsilon)
	{
		return notEqual(x, y, vec<C, T, Q>(epsilon));
	}

	template<length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<C, bool, Q> notEqual(mat<C, R, T, Q> const& a, mat<C, R, T, Q> const& b, vec<C, T, Q> const& epsilon)
	{
		vec<C, bool, Q> Result;
		for(length_t i = 0, n = C; i < n; ++i)
			Result[i] = any(notEqual(a[i], b[i], epsilon[i]));
		return Result;
	}
}//namespace glm
