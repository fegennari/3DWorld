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
/// @file gli/core/view.hpp
/// @date 2013-02-07 / 2013-02-07
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#ifndef GLI_CORE_VIEW_INCLUDED
#define GLI_CORE_VIEW_INCLUDED

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
	image view(image const & Image);

	template <typename texture>
	texture view(texture const & Texture);

	texture1D view(
		texture1D const & Texture,
		texture1D::size_type const & BaseLevel,
		texture1D::size_type const & MaxLevel);

	texture2D view(
		texture2D const & Texture,
		texture2D::size_type const & BaseLevel,
		texture2D::size_type const & MaxLevel);

	texture3D view(
		texture3D const & Texture,
		texture3D::size_type const & BaseLevel,
		texture3D::size_type const & MaxLevel);

	// texture can be texture1DArray, texture2DArray
	template <typename texture>
	texture view(
		texture const & Texture,
		typename texture::size_type const & BaseLayer,
		typename texture::size_type const & MaxLayer,
		typename texture::size_type const & BaseLevel,
		typename texture::size_type const & MaxLevel);

	textureCube view(
		textureCube const & Texture,
		textureCube::size_type const & BaseFace,
		textureCube::size_type const & MaxFace,
		textureCube::size_type const & BaseLevel,
		textureCube::size_type const & MaxLevel);

	textureCubeArray view(
		textureCubeArray const & Texture,
		textureCubeArray::size_type const & BaseLayer,
		textureCubeArray::size_type const & MaxLayer,
		textureCubeArray::size_type const & BaseFace,
		textureCubeArray::size_type const & MaxFace,
		textureCubeArray::size_type const & BaseLevel,
		textureCubeArray::size_type const & MaxLevel);
}//namespace gli

#include "view.inl"

#endif//GLI_CORE_VIEW_INCLUDED
