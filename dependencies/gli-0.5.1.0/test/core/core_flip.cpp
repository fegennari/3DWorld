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

template <typename texture, typename genType>
int test_texture
(
	texture const & Texture,
	genType const & ClearColor,
	genType const & FirstColor
)
{
	int Error(0);

	texture TextureA = gli::copy(Texture);
	TextureA.template clear<genType>(ClearColor);
	*TextureA.template data<genType>() = FirstColor;

	texture TextureB = gli::flip(TextureA);
	Error += TextureA != TextureB ? 0 : 1;

	texture TextureC = gli::flip(TextureB);
	Error += TextureC == TextureA ? 0 : 1;

	return Error;
}

int main()
{
	int Error(0);

	std::size_t const TextureSize = 32;
	std::size_t const Levels = gli::level_count(TextureSize);

	Error += test_texture(
		gli::texture2D(Levels, gli::R8_UNORM, gli::texture2D::dimensions_type(TextureSize)),
		glm::uint8(255), glm::uint8(0));

	Error += test_texture(
		gli::texture2D(Levels, gli::RGB8_UNORM, gli::texture2D::dimensions_type(TextureSize)),
		glm::u8vec3(255, 128, 0), glm::u8vec3(0, 128, 255));

	Error += test_texture(
		gli::texture2D(Levels, gli::RGBA8_UNORM, gli::texture2D::dimensions_type(TextureSize)),
		glm::u8vec4(255, 128, 0, 255), glm::u8vec4(0, 128, 255, 255));

	Error += test_texture(
		gli::texture2D(Levels, gli::RGBA32F, gli::texture2D::dimensions_type(TextureSize)),
		glm::f32vec4(1.0, 0.5, 0.0, 1.0), glm::f32vec4(0.0, 0.5, 1.0, 1.0));

	Error += test_texture(
		gli::texture2DArray(4, Levels, gli::RGBA8_UNORM, gli::texture2DArray::dimensions_type(TextureSize)),
		glm::u8vec4(255, 128, 0, 255), glm::u8vec4(0, 128, 255, 255));

	Error += test_texture(
		gli::texture2DArray(4, Levels, gli::RGBA32F, gli::texture2DArray::dimensions_type(TextureSize)),
		glm::f32vec4(1.0, 0.5, 0.0, 1.0), glm::f32vec4(0.0, 0.5, 1.0, 1.0));

	return Error;
}
