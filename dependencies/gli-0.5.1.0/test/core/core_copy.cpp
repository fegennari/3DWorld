//////////////////////////////////////////////////////////////////////////////////
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
/// @file gli/core/copy.cpp
/// @date 2013-02-04 / 2013-11-25
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <gli/gli.hpp>

int test_texture1D
(
	std::vector<gli::format> const & Formats, 
	gli::texture1D::size_type const & TextureSize
)
{
	int Error(0);

	for(std::size_t i = 0; i < Formats.size(); ++i)
	{
		gli::texture1D TextureA(
			gli::texture1D::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::texture1D::dimensions_type(TextureSize));

		gli::texture1D TextureB = gli::copy(TextureA);

		Error += TextureA == TextureB ? 0 : 1;

		gli::texture1D TextureC(TextureA, 
			gli::texture1D::size_type(1), gli::texture1D::size_type(2));

		Error += TextureA[1] == TextureC[0] ? 0 : 1;
		Error += TextureA[2] == TextureC[1] ? 0 : 1;

		gli::texture1D TextureD = gli::copy(TextureC);

		Error += TextureC == TextureD ? 0 : 1;

		gli::texture1D TextureG = gli::copy(TextureA, 0, TextureA.levels() - 1);
		Error += TextureA == TextureG ? 0 : 1;

		gli::texture1D TextureE = gli::copy(TextureA, 1, TextureA.levels() - 2);
		Error += TextureA[1] == TextureE[0] ? 0 : 1;

		gli::texture1D TextureF(TextureA, 1, TextureA.levels() - 2); 

		Error += TextureE == TextureF ? 0 : 1;
	}

	return Error;
}

int test_texture1DArray
(
	std::vector<gli::format> const & Formats, 
	gli::texture1D::size_type const & TextureSize
)
{
	int Error(0);

	for(std::size_t i = 0; i < Formats.size(); ++i)
	{
		gli::texture1DArray TextureA(
			gli::texture1DArray::size_type(4),
			gli::texture1DArray::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::texture1DArray::dimensions_type(TextureSize));

		gli::texture1DArray TextureB = gli::copy(TextureA);

		Error += TextureA == TextureB ? 0 : 1;

		gli::texture1DArray TextureC(TextureA, 
			gli::texture1DArray::size_type(0), TextureA.layers() - 1,
			gli::texture1DArray::size_type(1), gli::texture1DArray::size_type(2));

		Error += TextureA[0][1] == TextureC[0][0] ? 0 : 1;
		Error += TextureA[0][2] == TextureC[0][1] ? 0 : 1;
		Error += TextureA[1][1] == TextureC[1][0] ? 0 : 1;
		Error += TextureA[1][2] == TextureC[1][1] ? 0 : 1;

		gli::texture1DArray TextureD = gli::copy(TextureC);

		Error += TextureC == TextureD ? 0 : 1;

		gli::texture1DArray TextureG = gli::copy(
			TextureA, 
			0, TextureA.layers() - 1, 
			0, TextureA.levels() - 1);
		Error += TextureA == TextureG ? 0 : 1;

		gli::texture1DArray TextureE = gli::copy(
			TextureA, 
			1, TextureA.layers() - 1, 
			0, TextureA.levels() - 1);
		Error += TextureA[1] == TextureE[0] ? 0 : 1;

		gli::texture1DArray TextureF(
			TextureA, 
			1, TextureA.layers() - 1, 
			0, TextureA.levels() - 1); 

		Error += TextureE == TextureF ? 0 : 1;

		gli::texture1DArray TextureK(
			gli::texture1DArray::size_type(4),
			gli::texture1DArray::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::texture1DArray::dimensions_type(TextureSize));

		gli::texture1DArray TextureH(TextureK, 1, 2, 1, 2);
		gli::texture1DArray TextureI = gli::copy(TextureH);

		Error += TextureH == TextureI ? 0 : 1;

		gli::texture1DArray TextureJ = gli::copy(TextureK, 1, 2, 1, 2);
		Error += TextureH == TextureJ ? 0 : 1;
		Error += TextureI == TextureJ ? 0 : 1;
	}

	return Error;
}

int test_texture2D
(
	std::vector<gli::format> const & Formats, 
	gli::texture2D::size_type const & TextureSize
)
{
	int Error(0);

	for(std::size_t i = 0; i < Formats.size(); ++i)
	{
		gli::texture2D TextureA(
			gli::texture2D::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::texture2D::dimensions_type(TextureSize));

		gli::texture2D TextureB = gli::copy(TextureA);
		Error += TextureA == TextureB ? 0 : 1;

		gli::texture2D TextureC(TextureA, 
			gli::texture2D::size_type(1), gli::texture2D::size_type(2));

		Error += TextureA[1] == TextureC[0] ? 0 : 1;
		Error += TextureA[2] == TextureC[1] ? 0 : 1;

		gli::texture2D TextureD = gli::copy(TextureC);

		Error += TextureC == TextureD ? 0 : 1;

		gli::texture2D TextureG = gli::copy(TextureA, 0, TextureA.levels() - 1);
		Error += TextureA == TextureG ? 0 : 1;

		gli::texture2D TextureE = gli::copy(TextureA, 1, TextureA.levels() - 1);
		Error += TextureA[1] == TextureE[0] ? 0 : 1;

		gli::texture2D TextureF(TextureA, 1, TextureA.levels() - 1); 

		Error += TextureE == TextureF ? 0 : 1;
	}

	return Error;
}

int test_texture2DArray
(
	std::vector<gli::format> const & Formats, 
	gli::texture1D::size_type const & TextureSize
)
{
	int Error(0);

	for(std::size_t i = 0; i < Formats.size(); ++i)
	{
		gli::texture2DArray TextureA(
			gli::texture2DArray::size_type(4),
			gli::texture2DArray::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::texture2DArray::dimensions_type(TextureSize));

		gli::texture2DArray TextureB = gli::copy(TextureA);

		Error += TextureA == TextureB ? 0 : 1;

		gli::texture2DArray TextureC(TextureA, 
			gli::texture2DArray::size_type(0), TextureA.layers() - 1,
			gli::texture2DArray::size_type(1), gli::texture2DArray::size_type(2));

		Error += TextureA[0][1] == TextureC[0][0] ? 0 : 1;
		Error += TextureA[0][2] == TextureC[0][1] ? 0 : 1;
		Error += TextureA[1][1] == TextureC[1][0] ? 0 : 1;
		Error += TextureA[1][2] == TextureC[1][1] ? 0 : 1;

		gli::texture2DArray TextureD = gli::copy(TextureC);

		Error += TextureC == TextureD ? 0 : 1;

		gli::texture2DArray TextureG = gli::copy(
			TextureA, 
			0, TextureA.layers() - 1, 
			0, TextureA.levels() - 1);
		Error += TextureA == TextureG ? 0 : 1;

		gli::texture2DArray TextureE = gli::copy(
			TextureA, 
			1, TextureA.layers() - 1, 
			0, TextureA.levels() - 1);
		Error += TextureA[1] == TextureE[0] ? 0 : 1;

		gli::texture2DArray TextureF(
			TextureA, 
			1, TextureA.layers() - 1, 
			0, TextureA.levels() - 1); 

		Error += TextureE == TextureF ? 0 : 1;

		gli::texture2DArray TextureK(
			gli::texture2DArray::size_type(4),
			gli::texture2DArray::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::texture2DArray::dimensions_type(TextureSize));

		gli::texture2DArray TextureH(TextureK, 1, 2, 1, 2);
		gli::texture2DArray TextureI = gli::copy(TextureH);

		Error += TextureH == TextureI ? 0 : 1;

		gli::texture2DArray TextureJ = gli::copy(TextureK, 1, 2, 1, 2);
		Error += TextureH == TextureJ ? 0 : 1;
		Error += TextureI == TextureJ ? 0 : 1;
	}

	return Error;
}

int test_texture3D
(
	std::vector<gli::format> const & Formats, 
	gli::texture3D::size_type const & TextureSize
)
{
	int Error(0);

	for(std::size_t i = 0; i < Formats.size(); ++i)
	{
		gli::texture3D TextureA(
			gli::texture3D::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::texture3D::dimensions_type(TextureSize));

		gli::texture3D TextureB = gli::copy(TextureA);

		Error += TextureA == TextureB ? 0 : 1;

		gli::texture3D TextureC(TextureA, 
			gli::texture3D::size_type(1), gli::texture3D::size_type(2));

		Error += TextureA[1] == TextureC[0] ? 0 : 1;
		Error += TextureA[2] == TextureC[1] ? 0 : 1;

		gli::texture3D TextureD = gli::copy(TextureC);

		Error += TextureC == TextureD ? 0 : 1;

		gli::texture3D TextureG = gli::copy(TextureA, 0, TextureA.levels() - 1);
		Error += TextureA == TextureG ? 0 : 1;

		gli::texture3D TextureE = gli::copy(TextureA, 1, TextureA.levels() - 1);
		Error += TextureA[1] == TextureE[0] ? 0 : 1;

		gli::texture3D TextureF(TextureA, 1, TextureA.levels() - 1); 

		Error += TextureE == TextureF ? 0 : 1;
	}

	return Error;
}

int test_textureCube
(
	std::vector<gli::format> const & Formats, 
	gli::texture1D::size_type const & TextureSize
)
{
	int Error(0);

	for(std::size_t i = 0; i < Formats.size(); ++i)
	{
		gli::textureCube TextureA(
			gli::textureCube::size_type(6),
			gli::textureCube::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::textureCube::dimensions_type(TextureSize));

		gli::textureCube TextureB = gli::copy(TextureA);

		Error += TextureA == TextureB ? 0 : 1;

		gli::textureCube TextureC(TextureA, 
			gli::textureCube::size_type(0), TextureA.faces() - 1,
			gli::textureCube::size_type(1), gli::textureCube::size_type(2));

		Error += TextureA[0][1] == TextureC[0][0] ? 0 : 1;
		Error += TextureA[0][2] == TextureC[0][1] ? 0 : 1;
		Error += TextureA[1][1] == TextureC[1][0] ? 0 : 1;
		Error += TextureA[1][2] == TextureC[1][1] ? 0 : 1;

		gli::textureCube TextureD = gli::copy(TextureC);

		Error += TextureC == TextureD ? 0 : 1;

		gli::textureCube TextureG = gli::copy(
			TextureA, 
			0, TextureA.faces() - 1, 
			0, TextureA.levels() - 1);
		Error += TextureA == TextureG ? 0 : 1;

		gli::textureCube TextureE = gli::copy(
			TextureA, 
			1, TextureA.faces() - 1, 
			0, TextureA.levels() - 1);
		Error += TextureA[1] == TextureE[0] ? 0 : 1;

		gli::textureCube TextureF(
			TextureA, 
			1, TextureA.faces() - 1, 
			0, TextureA.levels() - 1); 

		Error += TextureE == TextureF ? 0 : 1;

		gli::textureCube TextureK(
			gli::textureCube::size_type(6),
			gli::textureCube::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::textureCube::dimensions_type(TextureSize));

		gli::textureCube TextureH(TextureK, 1, 2, 1, 2);
		gli::textureCube TextureI = gli::copy(TextureH);

		Error += TextureH == TextureI ? 0 : 1;

		gli::textureCube TextureJ = gli::copy(TextureK, 1, 2, 1, 2);
		Error += TextureH == TextureJ ? 0 : 1;
		Error += TextureI == TextureJ ? 0 : 1;
	}

	return Error;
}

int test_textureCubeArray
(
	std::vector<gli::format> const & Formats, 
	gli::texture1D::size_type const & TextureSize
)
{
	int Error(0);

	for(std::size_t i = 0; i < Formats.size(); ++i)
	{
		gli::textureCubeArray TextureA(
			gli::textureCubeArray::size_type(4),
			gli::textureCubeArray::size_type(6),
			gli::textureCubeArray::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::textureCubeArray::dimensions_type(TextureSize));

		gli::textureCubeArray TextureB = gli::copy(TextureA);

		Error += TextureA == TextureB ? 0 : 1;

		gli::textureCubeArray TextureC(TextureA, 
			gli::textureCubeArray::size_type(0), TextureA.layers() - 1,
			gli::textureCubeArray::size_type(0), TextureA.faces() - 1,
			gli::textureCubeArray::size_type(1), gli::textureCubeArray::size_type(2));

		Error += TextureA[0][0][1] == TextureC[0][0][0] ? 0 : 1;
		Error += TextureA[0][0][2] == TextureC[0][0][1] ? 0 : 1;
		Error += TextureA[0][1][1] == TextureC[0][1][0] ? 0 : 1;
		Error += TextureA[0][1][2] == TextureC[0][1][1] ? 0 : 1;

		gli::textureCubeArray TextureD = gli::copy(TextureC);

		Error += TextureC == TextureD ? 0 : 1;

		gli::textureCubeArray TextureG = gli::copy(
			TextureA, 
			gli::textureCubeArray::size_type(0), TextureA.layers() - 1,
			gli::textureCubeArray::size_type(0), TextureA.faces() - 1, 
			gli::textureCubeArray::size_type(0), TextureA.levels() - 1);
		Error += TextureA == TextureG ? 0 : 1;

		gli::textureCubeArray TextureE = gli::copy(
			TextureA,
			gli::textureCubeArray::size_type(0), TextureA.layers() - 1,
			gli::textureCubeArray::size_type(1), TextureA.faces() - 1,
			gli::textureCubeArray::size_type(0), TextureA.levels() - 1);
		Error += TextureA[0][1] == TextureE[0][0] ? 0 : 1;

		gli::textureCubeArray TextureF(
			TextureA, 
			gli::textureCubeArray::size_type(0), TextureA.layers() - 1,
			gli::textureCubeArray::size_type(1), TextureA.faces() - 1, 
			gli::textureCubeArray::size_type(0), TextureA.levels() - 1); 

		Error += TextureE == TextureF ? 0 : 1;

		gli::textureCubeArray TextureK(
			gli::textureCubeArray::size_type(4),
			gli::textureCubeArray::size_type(6),
			gli::textureCubeArray::size_type(glm::log2(float(TextureSize))) + 1,
			Formats[i],
			gli::textureCubeArray::dimensions_type(TextureSize));

		gli::textureCubeArray TextureH(TextureK, 1, 2, 1, 2, 1, 2);
		gli::textureCubeArray TextureI = gli::copy(TextureH);

		Error += TextureH == TextureI ? 0 : 1;

		gli::textureCubeArray TextureJ = gli::copy(TextureK, 1, 2, 1, 2, 1, 2);
		Error += TextureH == TextureJ ? 0 : 1;
		Error += TextureI == TextureJ ? 0 : 1;
	}

	return Error;
}

int main()
{
	int Error(0);

	std::vector<gli::format> Formats;
	Formats.push_back(gli::RGBA8_UNORM);
	Formats.push_back(gli::RGB8_UNORM);
	Formats.push_back(gli::RGB_DXT1);
	Formats.push_back(gli::RGB_BP_UNORM);
	Formats.push_back(gli::RGBA32F);
	std::size_t const TextureSize = 32;

	Error += test_textureCubeArray(Formats, TextureSize);
	Error += test_texture1D(Formats, TextureSize);
	Error += test_texture1DArray(Formats, TextureSize);
	Error += test_texture2D(Formats, TextureSize);
	Error += test_texture2DArray(Formats, TextureSize);
	Error += test_texture3D(Formats, TextureSize);
	Error += test_textureCube(Formats, TextureSize);

	return Error;
}
