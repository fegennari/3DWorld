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
/// @file gli/core/addressing.inl
/// @date 2012-11-19 / 2012-11-19
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli{
namespace detail
{
	inline size_type imageAddressing
	(
		storage const & Storage,
		size_type const & LayerOffset, 
		size_type const & FaceOffset, 
		size_type const & LevelOffset
	)
	{
		assert(LayerOffset < Storage.layers());
		assert(FaceOffset < Storage.faces());
		assert(LevelOffset < Storage.levels());

		size_type LayerSize = Storage.layerSize(0, Storage.faces() - 1, 0, Storage.levels() - 1);
		size_type FaceSize = Storage.faceSize(0, Storage.levels() - 1);
		size_type BaseOffset = LayerSize * LayerOffset + FaceSize * FaceOffset; 

		for(size_type Level(0); Level < LevelOffset; ++Level)
			BaseOffset += Storage.levelSize(Level);

		return BaseOffset;
	}

	inline size_type texelLinearAdressing
	(
		dimensions1_type const & Dimensions,
		dimensions1_type const & TexelCoord
	)
	{
		assert(TexelCoord < Dimensions);

		return TexelCoord;
	}

	inline size_type texelLinearAdressing
	(
		dimensions2_type const & Dimensions,
		dimensions2_type const & TexelCoord
	)
	{
		assert(TexelCoord.x < Dimensions.x);
		assert(TexelCoord.y < Dimensions.y);

		return TexelCoord.x + Dimensions.x * TexelCoord.y;
	}

	inline size_type texelLinearAdressing
	(
		dimensions3_type const & Dimensions,
		dimensions3_type const & TexelCoord
	)
	{
		assert(TexelCoord.x < Dimensions.x);
		assert(TexelCoord.y < Dimensions.y);
		assert(TexelCoord.z < Dimensions.z);

		return TexelCoord.x + Dimensions.x * (TexelCoord.y + Dimensions.y * TexelCoord.z);
	}

	inline size_type texelMortonAdressing
	(
		dimensions1_type const & Dimensions,
		dimensions1_type const & TexelCoord
	)
	{
		assert(TexelCoord < Dimensions);

		return TexelCoord;
	}

	inline size_type texelMortonAdressing
	(
		dimensions2_type const & Dimensions,
		dimensions2_type const & TexelCoord
	)
	{
		assert(TexelCoord.x < Dimensions.x);
		assert(TexelCoord.y < Dimensions.y);

		return glm::bitfieldInterleave(TexelCoord.x, TexelCoord.y);
	}

	inline size_type texelMortonAdressing
	(
		dimensions3_type const & Dimensions,
		dimensions3_type const & TexelCoord
	)
	{
		assert(TexelCoord.x < Dimensions.x);
		assert(TexelCoord.y < Dimensions.y);
		assert(TexelCoord.z < Dimensions.z);

		return glm::bitfieldInterleave(TexelCoord.x, TexelCoord.y, TexelCoord.z);
	}
}//namespace detail
}//namespace gli
