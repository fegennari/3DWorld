///////////////////////////////////////////////////////////////////////////////////
/// OpenGL Image (gli.g-truc.net)
///
/// Copyright (c) 2008 - 2013 G-Truc Creation (www.g-truc.net)
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
/// 
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
/// 
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.
///
/// @ref core
/// @file gli/core/comparison.hpp
/// @date 2013-02-04 / 2013-02-04
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#ifndef GLI_CORE_COMPARISON_INCLUDED
#define GLI_CORE_COMPARISON_INCLUDED

#include "image.hpp"
#include "texture1d.hpp"
#include "texture1d_array.hpp"
#include "texture2d.hpp"
#include "texture2d_array.hpp"
#include "texture3d.hpp"
#include "texture_cube.hpp"
#include "texture_cube_array.hpp"

namespace gli
{
	bool operator==(image const & ImageA, image const & ImageB);
	bool operator!=(image const & ImageA, image const & ImageB);

	bool operator==(gli::texture1D const & A, gli::texture1D const & B);
	bool operator==(gli::texture1DArray const & A, gli::texture1DArray const & B);
	bool operator==(gli::texture2D const & A, gli::texture2D const & B);
	bool operator==(gli::texture2DArray const & A, gli::texture2DArray const & B);
	bool operator==(gli::texture3D const & A, gli::texture3D const & B);
	bool operator==(gli::textureCube const & A, gli::textureCube const & B);
	bool operator==(gli::textureCubeArray const & A, gli::textureCubeArray const & B);
	bool operator!=(gli::texture1D const & A, gli::texture1D const & B);
	bool operator!=(gli::texture1DArray const & A, gli::texture1DArray const & B);
	bool operator!=(gli::texture2D const & A, gli::texture2D const & B);
	bool operator!=(gli::texture2DArray const & A, gli::texture2DArray const & B);
	bool operator!=(gli::texture3D const & A, gli::texture3D const & B);
	bool operator!=(gli::textureCube const & A, gli::textureCube const & B);
	bool operator!=(gli::textureCubeArray const & A, gli::textureCubeArray const & B);
}//namespace gli

#include "comparison.inl"

#endif//GLI_CORE_COMPARISON_INCLUDED
