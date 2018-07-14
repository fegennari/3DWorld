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
/// @file gli/core/core_comparison.cpp
/// @date 2013-02-04 / 2013-02-04
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <gli/gli.hpp>
#include <map>

int test_texture1D()
{
	int Error(0);

	std::vector<glm::u8vec4> Color(6);
	Color.push_back(glm::u8vec4(255,   0,   0, 255));
	Color.push_back(glm::u8vec4(255, 127,   0, 255));
	Color.push_back(glm::u8vec4(255, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255, 255, 255));
	Color.push_back(glm::u8vec4(  0,   0, 255, 255));

	gli::texture1D TextureA(
		glm::log2(32) + 1, 
		gli::RGBA8_UNORM, 
		gli::texture1D::dimensions_type(32));

	{
		gli::texture1D TextureB(
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture1D::dimensions_type(32));

		Error += TextureA == TextureB ? 0 : 1;
		assert(!Error);
		Error += TextureA != TextureB ? 1 : 0;
		assert(!Error);
	}

	{
		gli::texture1D TextureC(TextureA);

		Error += TextureA == TextureC ? 0 : 1;
		assert(!Error);
		Error += TextureA != TextureC ? 1 : 0;
		assert(!Error);
	}

	{
		gli::texture1D TextureD(TextureA, 
			TextureA.baseLevel(), TextureA.maxLevel());

		Error += TextureA == TextureD ? 0 : 1;
		assert(!Error);
		Error += TextureA != TextureD ? 1 : 0;
		assert(!Error);
	}

	{
		gli::texture1D TextureE(
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture1D::dimensions_type(32));

		*TextureE[TextureE.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureE ? 0 : 1;
		assert(!Error);
		Error += TextureA == TextureE ? 1 : 0;
		assert(!Error);
	}

	{
		gli::texture1D TextureB(
			1, 
			gli::RGBA8_UNORM, 
			gli::texture1D::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		assert(!Error);
		Error += TextureA == TextureB ? 1 : 0;
		assert(!Error);
	}

	{
		gli::texture1D TextureB(
			glm::log2(32) + 1, 
			gli::RGBA8_SNORM, 
			gli::texture1D::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		assert(!Error);
		Error += TextureA == TextureB ? 1 : 0;
		assert(!Error);
	}

	{
		gli::texture1D TextureB(
			glm::log2(64) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture1D::dimensions_type(64));

		gli::texture1D TextureC(TextureB, 
			TextureB.baseLevel() + 1, TextureB.maxLevel());

		Error += TextureA == TextureC ? 0 : 1;
		assert(!Error);
		Error += TextureA != TextureC ? 1 : 0;
		assert(!Error);
	}

	return Error;
}

int test_texture1DArray()
{
	int Error(0);

	std::vector<glm::u8vec4> Color(6);
	Color.push_back(glm::u8vec4(255,   0,   0, 255));
	Color.push_back(glm::u8vec4(255, 127,   0, 255));
	Color.push_back(glm::u8vec4(255, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255, 255, 255));
	Color.push_back(glm::u8vec4(  0,   0, 255, 255));

	gli::texture1DArray TextureA(
		gli::texture1D::size_type(1),
		glm::log2(32) + 1, 
		gli::RGBA8_UNORM, 
		gli::texture1D::dimensions_type(32));

	{
		gli::texture1DArray TextureB(
			gli::texture1D::size_type(1),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture1D::dimensions_type(32));

		Error += TextureA == TextureB ? 0 : 1;
		Error += TextureA != TextureB ? 1 : 0;
	}

	{
		gli::texture1DArray TextureC(TextureA);

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	{
		gli::texture1DArray TextureD(TextureA, 
			TextureA.baseLayer(),
			TextureA.maxLayer(),
			TextureA.baseLevel(),
			TextureA.maxLevel());

		Error += TextureA == TextureD ? 0 : 1;
		Error += TextureA != TextureD ? 1 : 0;
	}

	{
		gli::texture1DArray TextureE(
			gli::texture1D::size_type(1),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture1DArray::dimensions_type(32));

		*TextureE[0][TextureE.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureE ? 0 : 1;
		Error += TextureA == TextureE ? 1 : 0;
	}

	{
		gli::texture1DArray TextureB(
			gli::texture1D::size_type(1),
			1, 
			gli::RGBA8_UNORM, 
			gli::texture1DArray::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture1DArray TextureB(
			gli::texture1D::size_type(1),
			glm::log2(32) + 1, 
			gli::RGBA8_SNORM, 
			gli::texture1DArray::dimensions_type(32));

		*TextureB[0][TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture1DArray TextureB(
			gli::texture1D::size_type(1),
			glm::log2(64) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture1DArray::dimensions_type(64));

		gli::texture1DArray TextureC(TextureB, 
			TextureB.baseLayer(), TextureB.maxLayer(),
			TextureB.baseLevel() + 1, TextureB.maxLevel());

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	return Error;
}

int test_texture2D()
{
	int Error(0);

	std::vector<glm::u8vec4> Color(6);
	Color.push_back(glm::u8vec4(255,   0,   0, 255));
	Color.push_back(glm::u8vec4(255, 127,   0, 255));
	Color.push_back(glm::u8vec4(255, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255, 255, 255));
	Color.push_back(glm::u8vec4(  0,   0, 255, 255));

	gli::texture2D TextureA(
		glm::log2(32) + 1, 
		gli::RGBA8_UNORM, 
		gli::texture2D::dimensions_type(32));

	{
		gli::texture2D TextureB(
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture2D::dimensions_type(32));

		Error += TextureA == TextureB ? 0 : 1;
		Error += TextureA != TextureB ? 1 : 0;
	}

	{
		gli::texture2D TextureC(TextureA);

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	{
		gli::texture2D TextureD(TextureA, 
			TextureA.baseLevel(), TextureA.maxLevel());

		Error += TextureA == TextureD ? 0 : 1;
		Error += TextureA != TextureD ? 1 : 0;
	}

	{
		gli::texture2D TextureE(
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture2D::dimensions_type(32));

		*TextureE[TextureE.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureE ? 0 : 1;
		Error += TextureA == TextureE ? 1 : 0;
	}

	{
		gli::texture2D TextureB(
			1, 
			gli::RGBA8_UNORM, 
			gli::texture2D::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture2D TextureB(
			glm::log2(32) + 1, 
			gli::RGBA8_SNORM, 
			gli::texture2D::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture2D TextureB(
			glm::log2(64) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture2D::dimensions_type(64));

		gli::texture2D TextureC(TextureB, 
			TextureB.baseLevel() + 1, TextureB.maxLevel());

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	return Error;
}

int test_texture2DArray()
{
	int Error(0);

	std::vector<glm::u8vec4> Color(6);
	Color.push_back(glm::u8vec4(255,   0,   0, 255));
	Color.push_back(glm::u8vec4(255, 127,   0, 255));
	Color.push_back(glm::u8vec4(255, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255, 255, 255));
	Color.push_back(glm::u8vec4(  0,   0, 255, 255));

	gli::texture2DArray TextureA(
		gli::texture2DArray::size_type(1),
		glm::log2(32) + 1, 
		gli::RGBA8_UNORM, 
		gli::texture2DArray::dimensions_type(32));

	{
		gli::texture2DArray TextureB(
			gli::texture2DArray::size_type(1),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture2DArray::dimensions_type(32));

		Error += TextureA == TextureB ? 0 : 1;
		Error += TextureA != TextureB ? 1 : 0;
	}

	{
		gli::texture2DArray TextureC(TextureA);

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	{
		gli::texture2DArray TextureD(TextureA, 
			TextureA.baseLayer(), TextureA.maxLayer(),
			TextureA.baseLevel(), TextureA.maxLevel());

		Error += TextureA == TextureD ? 0 : 1;
		Error += TextureA != TextureD ? 1 : 0;
	}

	{
		gli::texture2DArray TextureE(
			gli::texture2DArray::size_type(1),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture2DArray::dimensions_type(32));

		*TextureE[0][TextureE.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureE ? 0 : 1;
		Error += TextureA == TextureE ? 1 : 0;
	}

	{
		gli::texture2DArray TextureB(
			gli::texture2DArray::size_type(1),
			gli::texture2DArray::size_type(1),
			gli::RGBA8_UNORM, 
			gli::texture2DArray::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture2DArray TextureB(
			gli::texture2DArray::size_type(1),
			glm::log2(32) + 1, 
			gli::RGBA8_SNORM, 
			gli::texture2DArray::dimensions_type(32));

		*TextureB[0][TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture2DArray TextureB(
			gli::texture2DArray::size_type(1),
			glm::log2(64) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture2DArray::dimensions_type(64));

		gli::texture2DArray TextureC(TextureB, 
			TextureB.baseLayer(), TextureB.maxLayer(),
			TextureB.baseLevel() + 1, TextureB.maxLevel());

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	return Error;
}

int test_texture3D()
{
	int Error(0);

	std::vector<glm::u8vec4> Color(6);
	Color.push_back(glm::u8vec4(255,   0,   0, 255));
	Color.push_back(glm::u8vec4(255, 127,   0, 255));
	Color.push_back(glm::u8vec4(255, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255, 255, 255));
	Color.push_back(glm::u8vec4(  0,   0, 255, 255));

	gli::texture3D TextureA(
		glm::log2(32) + 1, 
		gli::RGBA8_UNORM, 
		gli::texture3D::dimensions_type(32));

	{
		gli::texture3D TextureB(
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture3D::dimensions_type(32));

		Error += TextureA == TextureB ? 0 : 1;
		Error += TextureA != TextureB ? 1 : 0;
	}

	{
		gli::texture3D TextureC(TextureA);

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	{
		gli::texture3D TextureD(TextureA, 
			TextureA.baseLevel(), TextureA.maxLevel());

		Error += TextureA == TextureD ? 0 : 1;
		Error += TextureA != TextureD ? 1 : 0;
	}

	{
		gli::texture3D TextureE(
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture3D::dimensions_type(32));

		*TextureE[TextureE.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureE ? 0 : 1;
		Error += TextureA == TextureE ? 1 : 0;
	}

	{
		gli::texture3D TextureB(
			1, 
			gli::RGBA8_UNORM, 
			gli::texture3D::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture3D TextureB(
			glm::log2(32) + 1, 
			gli::RGBA8_SNORM, 
			gli::texture3D::dimensions_type(32));

		*TextureB[TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::texture3D TextureB(
			glm::log2(64) + 1, 
			gli::RGBA8_UNORM, 
			gli::texture3D::dimensions_type(64));

		gli::texture3D TextureC(TextureB, 
			TextureB.baseLevel() + 1, TextureB.maxLevel());

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	return Error;
}

int test_textureCube()
{
	int Error(0);

	std::vector<glm::u8vec4> Color(6);
	Color.push_back(glm::u8vec4(255,   0,   0, 255));
	Color.push_back(glm::u8vec4(255, 127,   0, 255));
	Color.push_back(glm::u8vec4(255, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255, 255, 255));
	Color.push_back(glm::u8vec4(  0,   0, 255, 255));

	gli::textureCube TextureA(
		gli::textureCube::size_type(6),
		glm::log2(32) + 1, 
		gli::RGBA8_UNORM, 
		gli::textureCube::dimensions_type(32));

	{
		gli::textureCube TextureB(
			gli::textureCube::size_type(6),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::textureCube::dimensions_type(32));

		Error += TextureA == TextureB ? 0 : 1;
		Error += TextureA != TextureB ? 1 : 0;
	}

	{
		gli::textureCube TextureC(TextureA);

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	{
		gli::textureCube TextureD(TextureA, 
			TextureA.baseFace(), TextureA.maxFace(),
			TextureA.baseLevel(), TextureA.maxLevel());

		Error += TextureA == TextureD ? 0 : 1;
		Error += TextureA != TextureD ? 1 : 0;
	}

	{
		gli::textureCube TextureE(
			gli::textureCube::size_type(6),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::textureCube::dimensions_type(32));

		*TextureE[TextureE.faces() - 1][TextureE.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureE ? 0 : 1;
		Error += TextureA == TextureE ? 1 : 0;
	}

	{
		gli::textureCube TextureB(
			gli::textureCube::size_type(6),
			gli::textureCube::size_type(1),
			gli::RGBA8_UNORM, 
			gli::textureCube::dimensions_type(32));

		*TextureB[TextureB.faces() - 1][TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::textureCube TextureB(
			gli::textureCube::size_type(6),
			glm::log2(32) + 1, 
			gli::RGBA8_SNORM, 
			gli::textureCube::dimensions_type(32));

		*TextureB[TextureB.faces() - 1][TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::textureCube TextureB(
			gli::textureCube::size_type(6),
			glm::log2(64) + 1, 
			gli::RGBA8_UNORM, 
			gli::textureCube::dimensions_type(64));

		gli::textureCube TextureC(TextureB, 
			TextureB.baseFace(), TextureB.maxFace(),
			TextureB.baseLevel() + 1, TextureB.maxLevel());

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	return Error;
}

int test_textureCubeArray()
{
	int Error(0);

	std::vector<glm::u8vec4> Color(6);
	Color.push_back(glm::u8vec4(255,   0,   0, 255));
	Color.push_back(glm::u8vec4(255, 127,   0, 255));
	Color.push_back(glm::u8vec4(255, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255,   0, 255));
	Color.push_back(glm::u8vec4(  0, 255, 255, 255));
	Color.push_back(glm::u8vec4(  0,   0, 255, 255));

	gli::textureCubeArray TextureA(
		gli::textureCubeArray::size_type(1),
		gli::textureCubeArray::size_type(6),
		glm::log2(32) + 1, 
		gli::RGBA8_UNORM, 
		gli::textureCubeArray::dimensions_type(32));

	{
		gli::textureCubeArray TextureB(
			gli::textureCubeArray::size_type(1),
			gli::textureCubeArray::size_type(6),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::textureCubeArray::dimensions_type(32));

		Error += TextureA == TextureB ? 0 : 1;
		Error += TextureA != TextureB ? 1 : 0;
	}

	{
		gli::textureCubeArray TextureC(TextureA);

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	{
		gli::textureCubeArray TextureD(TextureA, 
			TextureA.baseLayer(), TextureA.maxLayer(),
			TextureA.baseFace(), TextureA.maxFace(),
			TextureA.baseLevel(), TextureA.maxLevel());

		Error += TextureA == TextureD ? 0 : 1;
		Error += TextureA != TextureD ? 1 : 0;
	}

	{
		gli::textureCubeArray TextureE(
			gli::textureCubeArray::size_type(1),
			gli::textureCubeArray::size_type(6),
			glm::log2(32) + 1, 
			gli::RGBA8_UNORM, 
			gli::textureCubeArray::dimensions_type(32));

		*TextureE[0][TextureE.faces() - 1][TextureE.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureE ? 0 : 1;
		Error += TextureA == TextureE ? 1 : 0;
	}

	{
		gli::textureCubeArray TextureB(
			gli::textureCubeArray::size_type(1),
			gli::textureCubeArray::size_type(6),
			gli::textureCubeArray::size_type(1),
			gli::RGBA8_UNORM, 
			gli::textureCubeArray::dimensions_type(32));

		*TextureB[0][TextureB.faces() - 1][TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::textureCubeArray TextureB(
			gli::textureCubeArray::size_type(1),
			gli::textureCubeArray::size_type(6),
			glm::log2(32) + 1, 
			gli::RGBA8_SNORM, 
			gli::textureCubeArray::dimensions_type(32));

		*TextureB[0][TextureB.faces() - 1][TextureB.levels() - 1].data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		Error += TextureA != TextureB ? 0 : 1;
		Error += TextureA == TextureB ? 1 : 0;
	}

	{
		gli::textureCubeArray TextureB(
			gli::textureCubeArray::size_type(1),
			gli::textureCubeArray::size_type(6),
			glm::log2(64) + 1, 
			gli::RGBA8_UNORM, 
			gli::textureCubeArray::dimensions_type(64));

		gli::textureCubeArray TextureC(TextureB, 
			TextureB.baseLayer(), TextureB.maxLayer(),
			TextureB.baseFace(), TextureB.maxFace(),
			TextureB.baseLevel() + 1, TextureB.maxLevel());

		Error += TextureA == TextureC ? 0 : 1;
		Error += TextureA != TextureC ? 1 : 0;
	}

	return Error;
}

class C
{};

int test_map()
{
	int Error(0);

	gli::texture2D TextureA(
		gli::RGBA8_UNORM, 
		gli::texture2D::dimensions_type(32));

	gli::texture2D TextureB(
		gli::RGBA8_UNORM, 
		gli::texture2D::dimensions_type(64));

	
	std::map<int, gli::texture2D> Map;

	Map.insert(std::make_pair(0, TextureA));
	Map.insert(std::make_pair(0, TextureB));
	
	std::map<int, C> MapC;
	MapC.insert(std::make_pair(0, C()));
	MapC.insert(std::make_pair(0, C()));

	return Error;
}

int main()
{
	int Error(0);

	Error += test_texture1D();
	Error += test_texture1DArray();
	Error += test_texture2D();
	Error += test_texture2DArray();
	Error += test_texture3D();
	Error += test_textureCube();
	Error += test_textureCubeArray();
	Error += test_map();
		
	return Error;
}
