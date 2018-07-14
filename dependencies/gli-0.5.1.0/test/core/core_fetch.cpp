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
/// @file gli/core/fetch.cpp
/// @date 2013-11-25 / 2013-11-25
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <gli/gli.hpp>
#include <glm/gtc/epsilon.hpp>

namespace textureLod
{
	int test()
	{
		int Error(0);

		gli::texture2D Texture(gli::load_dds("rgba4.dds"));
		glm::u8vec3 TexelA = gli::textureLod<glm::u8vec3>(Texture, gli::texture2D::texcoord_type(0.0f, 0.0f), 0);

		return Error;
	}
}//namespace fetch

namespace fetch_rgba8_unorm
{
	int test()
	{
		int Error(0);

		{
			gli::texture2D Texture(gli::RGBA8_UNORM, gli::texture2D::dimensions_type(2, 2));
			gli::texelWrite<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(0, 0), 0, glm::u8vec4(255, 0, 0, 255));
			gli::texelWrite<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(1, 0), 0, glm::u8vec4(255, 255, 0, 255));
			gli::texelWrite<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(1, 1), 0, glm::u8vec4(0, 255, 0, 255));
			gli::texelWrite<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(0, 1), 0, glm::u8vec4(0, 0, 255, 255));
			gli::texelWrite<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(0, 0), 1, glm::u8vec4(255, 128, 0, 255));
			gli::save_dds(Texture, "rgba8_unorm_4pixels.dds");
		}

		{
			gli::texture2D Texture(gli::load_dds("rgba8_unorm_4pixels.dds"));
			glm::u8vec4 A = gli::texelFetch<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(0, 0), 0);
			Error += glm::all(glm::equal(A, glm::u8vec4(255, 0, 0, 255)));
			glm::u8vec4 B = gli::texelFetch<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(1, 0), 0);
			Error += glm::all(glm::equal(B, glm::u8vec4(255, 255, 0, 255)));
			glm::u8vec4 C = gli::texelFetch<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(1, 1), 0);
			Error += glm::all(glm::equal(C, glm::u8vec4(0, 255, 0, 255)));
			glm::u8vec4 D = gli::texelFetch<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(0, 1), 0);
			Error += glm::all(glm::equal(D, glm::u8vec4(0, 0, 255, 255)));
			glm::u8vec4 E = gli::texelFetch<glm::u8vec4>(Texture, gli::texture2D::dimensions_type(0, 0), 1);
			Error += glm::all(glm::equal(E, glm::u8vec4(255, 128, 0, 255)));
		}

		return Error;
	}
}//namespace fetch

namespace fetch_rgb32f
{
	int test()
	{
		int Error(0);

		{
			gli::texture2D Texture(gli::RGB32F, gli::texture2D::dimensions_type(2, 2));
			gli::texelWrite<glm::vec3>(Texture, gli::texture2D::dimensions_type(0, 0), 0, glm::vec3(1.0f, 0.0f, 0.0f));
			gli::texelWrite<glm::vec3>(Texture, gli::texture2D::dimensions_type(1, 0), 0, glm::vec3(1.0f, 1.0f, 0.0f));
			gli::texelWrite<glm::vec3>(Texture, gli::texture2D::dimensions_type(1, 1), 0, glm::vec3(0.0f, 1.0f, 0.0f));
			gli::texelWrite<glm::vec3>(Texture, gli::texture2D::dimensions_type(0, 1), 0, glm::vec3(0.0f, 0.0f, 1.0f));
			gli::texelWrite<glm::vec3>(Texture, gli::texture2D::dimensions_type(0, 0), 1, glm::vec3(1.0f, 0.5f, 0.0f));
			gli::save_dds(Texture, "rgba32f_4pixels.dds");
		}

		{
			gli::texture2D Texture(gli::load_dds("rgba32f_4pixels.dds"));
			glm::vec3 A = gli::texelFetch<glm::vec3>(Texture, gli::texture2D::dimensions_type(0, 0), 0);
			Error += glm::all(glm::epsilonEqual(A, glm::vec3(1.0f, 0.0f, 0.0f), 0.001f)) ? 0 : 1;
			glm::vec3 B = gli::texelFetch<glm::vec3>(Texture, gli::texture2D::dimensions_type(1, 0), 0);
			Error += glm::all(glm::epsilonEqual(B, glm::vec3(1.0f, 1.0f, 0.0f), 0.001f)) ? 0 : 1;
			glm::vec3 C = gli::texelFetch<glm::vec3>(Texture, gli::texture2D::dimensions_type(1, 1), 0);
			Error += glm::all(glm::epsilonEqual(C, glm::vec3(0.0f, 1.0f, 0.0f), 0.001f)) ? 0 : 1;
			glm::vec3 D = gli::texelFetch<glm::vec3>(Texture, gli::texture2D::dimensions_type(0, 1), 0);
			Error += glm::all(glm::epsilonEqual(D, glm::vec3(0.0f, 0.0f, 1.0f), 0.001f)) ? 0 : 1;
			glm::vec3 E = gli::texelFetch<glm::vec3>(Texture, gli::texture2D::dimensions_type(0, 0), 1);
			Error += glm::all(glm::epsilonEqual(E, glm::vec3(1.0f, 0.5f, 0.0f), 0.001f)) ? 0 : 1;
		}

		return Error;
	}
}//namespace fetch

int main()
{
	int Error(0);
	
	Error += fetch_rgba8_unorm::test();
	Error += fetch_rgb32f::test();
	Error += textureLod::test();

	return Error;
}
