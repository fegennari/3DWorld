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
/// @file gli/core/texture_2d_array.cpp
/// @date 2011-10-11 / 2013-11-25
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <gli/gli.hpp>

int test_alloc()
{
	int Error(0);

	std::vector<gli::format> Formats;
	Formats.push_back(gli::RGBA8_UNORM);
	Formats.push_back(gli::RGB8_UNORM);
	Formats.push_back(gli::R8_SNORM);
	Formats.push_back(gli::RGB_DXT1);
	Formats.push_back(gli::RGB_BP_UNORM);
	Formats.push_back(gli::RGBA32F);

	std::vector<std::size_t> Sizes;
	Sizes.push_back(16);
	Sizes.push_back(32);
	Sizes.push_back(15);
	Sizes.push_back(17);
	Sizes.push_back(1);

	for(std::size_t FormatIndex = 0; FormatIndex < Formats.size(); ++FormatIndex)
	for(std::size_t SizeIndex = 0; SizeIndex < Sizes.size(); ++SizeIndex)
	{
		gli::texture2DArray TextureA(
			1,
			gli::texture2DArray::size_type(glm::log2(int(Sizes[SizeIndex])) + 1),
			Formats[FormatIndex],
			gli::texture2DArray::dimensions_type(Sizes[SizeIndex]));

		gli::texture2DArray TextureB(
			1,
			Formats[FormatIndex],
			gli::texture2DArray::dimensions_type(Sizes[SizeIndex]));

		Error += TextureA == TextureB ? 0 : 1;
	}

	return Error;
}

int test_texture2DArray_clear()
{
	int Error(0);

	glm::u8vec4 const Orange(255, 127, 0, 255);

	gli::texture2DArray Texture(
		2,
		gli::texture2DArray::size_type(glm::log2(16u) + 1),
		gli::RGBA8U,
		gli::texture2DArray::dimensions_type(16));

	Texture.clear<glm::u8vec4>(Orange);

	return Error;
}

int test_texture2DArray_query()
{
	int Error(0);

	gli::texture2DArray::size_type Layers(2);
	gli::texture2DArray::size_type Levels(2);

	gli::texture2DArray Texture(
		Layers, 
		Levels,
		gli::RGBA8U,
		gli::texture2DArray::dimensions_type(2));

	gli::texture2DArray::size_type Size = Texture.size();

	Error += Size == sizeof(glm::u8vec4) * 5 * Layers ? 0 : 1;
	Error += Texture.format() == gli::RGBA8U ? 0 : 1;
	Error += Texture.layers() == Layers ? 0 : 1;
	Error += Texture.levels() == Levels ? 0 : 1;
	Error += !Texture.empty() ? 0 : 1;
	Error += Texture.dimensions().x == 2 ? 0 : 1;
	Error += Texture.dimensions().y == 2 ? 0 : 1;

	return Error;
}

int test_texture2DArray_access()
{
	int Error(0);

	{
		gli::texture2DArray Texture2DArray(
			gli::texture2DArray::size_type(2),
			gli::texture2DArray::size_type(1),
			gli::RGBA8U,
			gli::texture2DArray::dimensions_type(2));
		assert(!Texture2DArray.empty());

		gli::texture2D Texture0 = Texture2DArray[0];
		gli::texture2D Texture1 = Texture2DArray[1];
		
		std::size_t Size0 = Texture0.size();
		std::size_t Size1 = Texture1.size();

		Error += Size0 == sizeof(glm::u8vec4) * 4 ? 0 : 1;
		Error += Size1 == sizeof(glm::u8vec4) * 4 ? 0 : 1;

		for(std::size_t i = 0; i < 4; ++i)
		{
			*(Texture0.data<glm::u8vec4>() + i) = glm::u8vec4(255, 127, 0, 255);
			*(Texture1.data<glm::u8vec4>() + i) = glm::u8vec4(0, 127, 255, 255);
		}

		glm::u8vec4 * PointerA = Texture0.data<glm::u8vec4>();
		glm::u8vec4 * PointerB = Texture1.data<glm::u8vec4>();

		glm::u8vec4 * Pointer0 = Texture2DArray.data<glm::u8vec4>() + 0;
		glm::u8vec4 * Pointer1 = Texture2DArray.data<glm::u8vec4>() + 4;

		Error += PointerA == Pointer0 ? 0 : 1;
		Error += PointerB == Pointer1 ? 0 : 1;

		glm::u8vec4 ColorA = *Texture0.data<glm::u8vec4>();
		glm::u8vec4 ColorB = *Texture1.data<glm::u8vec4>();

		glm::u8vec4 Color0 = *Pointer0;
		glm::u8vec4 Color1 = *Pointer1;

		Error += glm::all(glm::equal(Color0, glm::u8vec4(255, 127, 0, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Color1, glm::u8vec4(0, 127, 255, 255))) ? 0 : 1;
	}

	return Error;
}

struct test
{
	test(
		gli::format const & Format,
		gli::texture2DArray::dimensions_type const & Dimensions,
		gli::texture2DArray::size_type const & Size) :
		Format(Format),
		Dimensions(Dimensions),
		Size(Size)
	{}

	gli::format Format;
	gli::texture2DArray::dimensions_type Dimensions;
	gli::texture2DArray::size_type Size;
};

int test_texture2DArray_size()
{
	int Error(0);

	std::vector<test> Tests;
	Tests.push_back(test(gli::RGBA8U, gli::texture2DArray::dimensions_type(4), 64 * 2));
	Tests.push_back(test(gli::R8U, gli::texture2DArray::dimensions_type(4), 16 * 2));
	Tests.push_back(test(gli::RGB_DXT1, gli::texture2DArray::dimensions_type(4), 8 * 2));
	Tests.push_back(test(gli::RGB_DXT1, gli::texture2DArray::dimensions_type(2), 8 * 2));
	Tests.push_back(test(gli::RGB_DXT1, gli::texture2DArray::dimensions_type(1), 8 * 2));
	Tests.push_back(test(gli::RGBA_DXT5, gli::texture2DArray::dimensions_type(4), 16 * 2));

	for(std::size_t i = 0; i < Tests.size(); ++i)
	{
		gli::texture2DArray Texture2DArray(
			gli::texture2DArray::size_type(2),
			gli::texture2DArray::size_type(1),
			Tests[i].Format,
			gli::texture2DArray::dimensions_type(4));

		Error += Texture2DArray.size() == Tests[i].Size ? 0 : 1;
		assert(!Error);
	}

	for(std::size_t i = 0; i < Tests.size(); ++i)
	{
		gli::texture2DArray Texture2DArray(
			gli::texture2DArray::size_type(2),
			gli::texture2DArray::size_type(1),
			Tests[i].Format,
			gli::texture2DArray::dimensions_type(4));

		gli::texture2D Texture2D = Texture2DArray[0];

		Error += Texture2DArray.size() == Tests[i].Size ? 0 : 1;
		assert(!Error);
	}

	return Error;
}

int main()
{
	int Error(0);

	Error += test_alloc();
	Error += test_texture2DArray_size();
	Error += test_texture2DArray_query();
	Error += test_texture2DArray_clear();
	Error += test_texture2DArray_access();

	return Error;
}

