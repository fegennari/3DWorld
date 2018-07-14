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
/// @file gli/core/view.inl
/// @date 2013-02-07 / 2013-02-07
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli
{
	inline image view(image const & Image)
	{
		return Image;
	}

	template <typename texture>
	inline texture view(texture const & Texture)
	{
		return Texture;
	}

	inline texture1D view
	(
		texture1D const & Texture,
		texture1D::size_type const & BaseLevel,
		texture1D::size_type const & MaxLevel
	)
	{
		assert(BaseLevel <= MaxLevel);
		assert(BaseLevel < Texture.levels());
		assert(MaxLevel < Texture.levels());

		return texture1D(Texture,
			BaseLevel, MaxLevel);
	}

	inline texture2D view
	(
		texture2D const & Texture,
		texture2D::size_type const & BaseLevel,
		texture2D::size_type const & MaxLevel
	)
	{
		assert(BaseLevel <= MaxLevel);
		assert(BaseLevel < Texture.levels());
		assert(MaxLevel < Texture.levels());

		return texture2D(Texture,
			BaseLevel, MaxLevel);
	}

	inline texture3D view
	(
		texture3D const & Texture,
		texture3D::size_type const & BaseLevel,
		texture3D::size_type const & MaxLevel
	)
	{
		assert(BaseLevel <= MaxLevel);
		assert(BaseLevel < Texture.levels());
		assert(MaxLevel < Texture.levels());

		return texture3D(Texture,
			BaseLevel, MaxLevel);
	}

	// texture can be texture1DArray, texture2DArray
	template <typename texture>
	inline texture view
	(
		texture const & Texture,
		typename texture::size_type const & BaseLayer,
		typename texture::size_type const & MaxLayer,
		typename texture::size_type const & BaseLevel,
		typename texture::size_type const & MaxLevel
	)
	{
		assert(BaseLevel <= MaxLevel);
		assert(BaseLevel < Texture.levels());
		assert(MaxLevel < Texture.levels());
		assert(BaseLayer <= MaxLayer);
		assert(BaseLayer < Texture.layers());
		assert(MaxLayer < Texture.layers());

		texture View(Texture,
			BaseLayer, MaxLayer,
			BaseLevel, MaxLevel);

		return View;
	}

	inline textureCube view
	(
		textureCube const & Texture,
		textureCube::size_type const & BaseFace,
		textureCube::size_type const & MaxFace,
		textureCube::size_type const & BaseLevel,
		textureCube::size_type const & MaxLevel
	)
	{
		assert(BaseLevel <= MaxLevel);
		assert(BaseLevel < Texture.levels());
		assert(MaxLevel < Texture.levels());
		assert(BaseFace <= MaxFace);
		assert(BaseFace < Texture.faces());
		assert(MaxFace < Texture.faces());

		textureCube View(Texture,
			BaseFace, MaxFace,
			BaseLevel, MaxLevel);

		return View;
	}

	inline textureCubeArray view
	(
		textureCubeArray const & Texture,
		textureCubeArray::size_type const & BaseLayer,
		textureCubeArray::size_type const & MaxLayer,
		textureCubeArray::size_type const & BaseFace,
		textureCubeArray::size_type const & MaxFace,
		textureCubeArray::size_type const & BaseLevel,
		textureCubeArray::size_type const & MaxLevel
	)
	{
		assert(BaseLevel <= MaxLevel);
		assert(BaseLevel < Texture.levels());
		assert(MaxLevel < Texture.levels());
		assert(BaseFace <= MaxFace);
		assert(BaseFace < Texture.faces());
		assert(MaxFace < Texture.faces());
		assert(BaseLayer <= MaxLayer);
		assert(BaseLayer < Texture.layers());
		assert(MaxLayer < Texture.layers());

		textureCubeArray View(Texture,
			BaseLayer, MaxLayer,
			BaseFace, MaxFace,
			BaseLevel, MaxLevel);

		return View;
	}

}//namespace gli
