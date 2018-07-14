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
/// @file gli/core/texture_2d.cpp
/// @date 2011-10-11 / 2013-11-25
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <gli/gli.hpp>

inline gli::texture2D radial
(
	gli::texture2D::dimensions_type const & Size, 
	gli::texture2D::texcoord_type const & Center,
	float const & Radius,
	gli::texture2D::texcoord_type const & Focal
)
{
	gli::texture2D Result(1, gli::RGB8U, gli::texture2D::dimensions_type(Size));
	glm::u8vec3 * DstData = (glm::u8vec3*)Result.data();

	for(std::size_t y = 0; y < Result.dimensions().y; ++y)
	for(std::size_t x = 0; x < Result.dimensions().x; ++x)
	{
		float Value = glm::radialGradient(
			Center * glm::vec2(Size), 
			Radius, 
			Focal * glm::vec2(Size),
			glm::vec2(x, y));

		std::size_t Index = x + y * Result.dimensions().x;

		*(DstData + Index) = glm::u8vec3(glm::u8(glm::clamp(Value * 255.f, 0.f, 255.f)));
	}

	return Result;
}

inline gli::texture2D linear
(
	gli::texture2D::dimensions_type const & Size, 
	gli::texture2D::texcoord_type const & Point0, 
	gli::texture2D::texcoord_type const & Point1
)
{
	gli::texture2D Result(1, gli::RGB8U, gli::texture2D::dimensions_type(Size));
	glm::u8vec3 * DstData = (glm::u8vec3*)Result.data();

	for(std::size_t y = 0; y < Result.dimensions().y; ++y)
	for(std::size_t x = 0; x < Result.dimensions().x; ++x)
	{
		float Value = glm::linearGradient(
			Point0 * glm::vec2(Size), 
			Point1 * glm::vec2(Size),
			gli::texture2D::texcoord_type(x, y));

		std::size_t Index = x + y * Result.dimensions().x;

		*(DstData + Index) = glm::u8vec3(glm::u8(glm::clamp(Value * 255.f, 0.f, 255.f)));
	}

	return Result;
}

int test_create()
{
	int Error(0);

	gli::texture2D TextureA = radial(
		gli::texture2D::dimensions_type(128), gli::texture2D::texcoord_type(0.5), 16.f, gli::texture2D::texcoord_type(0.7));

	gli::texture2D TextureB = linear(
		gli::texture2D::dimensions_type(128), gli::texture2D::texcoord_type(0.5), gli::texture2D::texcoord_type(0.7));

	Error += TextureA != TextureB ? 0 : 1;

	return Error;
}

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
		gli::texture2D TextureA(
			gli::texture2D::size_type(glm::log2(int(Sizes[SizeIndex])) + 1),
			Formats[FormatIndex],
			gli::texture2D::dimensions_type(Sizes[SizeIndex]));

		gli::texture2D TextureB(
			Formats[FormatIndex],
			gli::texture2D::dimensions_type(Sizes[SizeIndex]));

		Error += TextureA == TextureB ? 0 : 1;
	}

	return Error;
}

int test_texture2d_clear()
{
	int Error(0);

	glm::u8vec4 const Orange(255, 127, 0, 255);

	gli::texture2D Texture(
		gli::texture2D::size_type(glm::log2(16u) + 1),
		gli::RGBA8U,
		gli::texture2D::dimensions_type(16));

	//Texture.clear<glm::u8vec4>(Orange);
	Texture.clear(Orange);

	return Error;
}

int test_texture2d_query()
{
	int Error(0);

	gli::texture2D Texture(
		gli::texture2D::size_type(2),
		gli::RGBA8U,
		gli::texture2D::dimensions_type(2));

	Error += Texture.size() == sizeof(glm::u8vec4) * 5 ? 0 : 1;
	Error += Texture.format() == gli::RGBA8U ? 0 : 1;
	Error += Texture.levels() == 2 ? 0 : 1;
	Error += !Texture.empty() ? 0 : 1;
	Error += Texture.dimensions().x == 2 ? 0 : 1;
	Error += Texture.dimensions().y == 2 ? 0 : 1;

	return Error;
}

int test_texture2d_image_access()
{
	int Error(0);

	{
		glm::u8vec4 const Orange(255, 127, 0, 255);

		gli::image Image0(gli::image::dimensions_type(2, 2, 1), sizeof(glm::u8vec4), gli::image::dimensions_type(1));
		for(std::size_t i = 0; i < Image0.size(); ++i)
			*(Image0.data<glm::byte>() + i) = glm::byte(i);

		gli::image Image1(gli::image::dimensions_type(1, 1, 1), sizeof(glm::u8vec4), gli::image::dimensions_type(1));
		for(std::size_t i = 0; i < Image1.size(); ++i)
			*(Image1.data<glm::byte>() + i) = glm::byte(i + 100);

		gli::texture2D Texture(
			gli::texture2D::size_type(2),
			gli::RGBA8U,
			gli::texture2D::dimensions_type(2));

		/// TODO copy function
		/// Texture[0] = Image0;
		/// Texture[1] = Image1;

		/// Error += Texture[0] == Image0 ? 0 : 1;
		/// Error += Texture[1] == Image1 ? 0 : 1;
	}

	{
		gli::texture2D Texture(
			gli::texture2D::size_type(2),
			gli::RGBA8U,
			gli::texture2D::dimensions_type(2));
		assert(!Texture.empty());

		gli::image Image0 = Texture[0];
		gli::image Image1 = Texture[1];
		
		std::size_t Size0 = Image0.size();
		std::size_t Size1 = Image1.size();

		Error += Size0 == sizeof(glm::u8vec4) * 4 ? 0 : 1;
		Error += Size1 == sizeof(glm::u8vec4) * 1 ? 0 : 1;

		*Image0.data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);
		*Image1.data<glm::u8vec4>() = glm::u8vec4(0, 127, 255, 255);

		glm::u8vec4 * PointerA = Image0.data<glm::u8vec4>();
		glm::u8vec4 * PointerB = Image1.data<glm::u8vec4>();

		glm::u8vec4 * Pointer0 = Texture.data<glm::u8vec4>() + 0;
		glm::u8vec4 * Pointer1 = Texture.data<glm::u8vec4>() + 4;

		Error += PointerA == Pointer0 ? 0 : 1;
		Error += PointerB == Pointer1 ? 0 : 1;

		glm::u8vec4 ColorA = *Image0.data<glm::u8vec4>();
		glm::u8vec4 ColorB = *Image1.data<glm::u8vec4>();

		glm::u8vec4 Color0 = *Pointer0;
		glm::u8vec4 Color1 = *Pointer1;

		Error += glm::all(glm::equal(Color0, glm::u8vec4(255, 127, 0, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Color1, glm::u8vec4(0, 127, 255, 255))) ? 0 : 1;
	}

	{
		gli::texture2D Texture(
			gli::texture2D::size_type(1),
			gli::RGBA8U,
			gli::texture2D::dimensions_type(2));

		std::size_t SizeA = Texture.size();
		Error += SizeA == sizeof(glm::u8vec4) * 4 ? 0 : 1;

		gli::image Image0 = Texture[0];
		
		std::size_t Size0 = Image0.size();
		Error += Size0 == sizeof(glm::u8vec4) * 4 ? 0 : 1;

		*Image0.data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);

		glm::u8vec4 * PointerA = Image0.data<glm::u8vec4>();
		glm::u8vec4 * Pointer0 = Texture.data<glm::u8vec4>() + 0;
		Error += PointerA == Pointer0 ? 0 : 1;

		glm::u8vec4 ColorA = *PointerA;
		glm::u8vec4 Color0 = *Pointer0;
		Error += glm::all(glm::equal(Color0, glm::u8vec4(255, 127, 0, 255))) ? 0 : 1;
	}

	return Error;
}

struct test
{
	test(
		gli::format const & Format,
		gli::texture2D::dimensions_type const & Dimensions,
		gli::texture2D::size_type const & Size) :
		Format(Format),
		Dimensions(Dimensions),
		Size(Size)
	{}

	gli::format Format;
	gli::texture2D::dimensions_type Dimensions;
	gli::texture2D::size_type Size;
};

int test_texture2d_image_size()
{
	int Error(0);

	std::vector<test> Tests;
	Tests.push_back(test(gli::RGBA8U, gli::texture2D::dimensions_type(4), 64));
	Tests.push_back(test(gli::R8U, gli::texture2D::dimensions_type(4), 16));
	Tests.push_back(test(gli::RGBA_DXT1, gli::texture2D::dimensions_type(4), 8));
	Tests.push_back(test(gli::RGBA_DXT1, gli::texture2D::dimensions_type(2), 8));
	Tests.push_back(test(gli::RGBA_DXT1, gli::texture2D::dimensions_type(1), 8));
	Tests.push_back(test(gli::RGBA_DXT5, gli::texture2D::dimensions_type(4), 16));

	for(std::size_t i = 0; i < Tests.size(); ++i)
	{
		gli::texture2D Texture(
			gli::texture2D::size_type(1),
			Tests[i].Format,
			gli::texture2D::dimensions_type(4));

		Error += Texture.size() == Tests[i].Size ? 0 : 1;
		assert(!Error);
	}

	for(std::size_t i = 0; i < Tests.size(); ++i)
	{
		gli::texture2D Texture(
			gli::texture2D::size_type(1),
			Tests[i].Format,
			gli::texture2D::dimensions_type(4));

		gli::image Image = Texture[0];

		Error += Image.size() == Tests[i].Size ? 0 : 1;
		assert(!Error);
	}

	return Error;
}

namespace fetch
{
	int test()
	{
		int Error(0);

		gli::texture2D Texture(1, gli::RGBA8U, gli::texture2D::dimensions_type(4, 2));
		*(Texture.data<glm::u8vec4>() + 0) = glm::u8vec4(255,   0,   0, 255);
		*(Texture.data<glm::u8vec4>() + 1) = glm::u8vec4(255, 128,   0, 255);
		*(Texture.data<glm::u8vec4>() + 2) = glm::u8vec4(255, 255,   0, 255);
		*(Texture.data<glm::u8vec4>() + 3) = glm::u8vec4(128, 255,   0, 255);
		*(Texture.data<glm::u8vec4>() + 4) = glm::u8vec4(  0, 255,   0, 255);
		*(Texture.data<glm::u8vec4>() + 5) = glm::u8vec4(  0, 255, 255, 255);
		*(Texture.data<glm::u8vec4>() + 6) = glm::u8vec4(  0,   0, 255, 255);
		*(Texture.data<glm::u8vec4>() + 7) = glm::u8vec4(255,   0, 255, 255);

		glm::u8vec4 Data0 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(0, 0), 0);
		glm::u8vec4 Data1 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(1, 0), 0);
		glm::u8vec4 Data2 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(2, 0), 0);
		glm::u8vec4 Data3 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(3, 0), 0);
		glm::u8vec4 Data4 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(0, 1), 0);
		glm::u8vec4 Data5 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(1, 1), 0);
		glm::u8vec4 Data6 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(2, 1), 0);
		glm::u8vec4 Data7 = Texture.fetch<glm::u8vec4>(gli::texture2D::dimensions_type(3, 1), 0);

		return Error;
	}
}//namespace fetch

int main()
{
	int Error(0);

	Error += test_alloc();
	Error += test_texture2d_image_size();
	Error += test_texture2d_query();
	Error += test_texture2d_clear();
	Error += test_texture2d_image_access();
	Error += test_create();
	Error += fetch::test();

	return Error;
}

