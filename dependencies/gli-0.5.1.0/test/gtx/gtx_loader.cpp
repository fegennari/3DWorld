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
/// @file gli/test/gtx_loader.cpp
/// @date 2009-11-14 / 2013-11-25
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#define GLM_SWIZZLE_XYZW
#include <glm/glm.hpp>
#include <glm/gtc/type_precision.hpp>

#include <gli/gli.hpp>
#include <gli/core/save_dds.hpp>
#include <gli/core/fetch.hpp>

#include <vector>
 
//#include <boost/format.hpp>

int test_image_wip()
{
	int Error(0);
	//gli::wip::texture2D<glm::u8vec3, gli::wip::plain> Texture;
	//gli::wip::texture2D<glm::u8vec3, gli::wip::plain>::image Mipmap = Texture[0];

	//glm::vec2 Texcoord(0);
	//Texture[0](Texcoord);

	//gli::wip::plain<glm::u8vec3> Surface;
	//gli::wip::fetch(Surface);
	//gli::wip::fetch(Texture);

	return Error;
}

int test_image_export()
{
	int Error(0);
	//gli::texture2D Texture = gli::load<gli::TGA>("../test_rgb8.tga");
	//gli::texture2D TextureMipmaped = gli::generateMipmaps(Texture, 0);

	//gli::save(TextureMipmaped, 0, "../test0.tga");
	//gli::save(TextureMipmaped, 1, "../test1.tga");
	//gli::save(TextureMipmaped, 2, "../test2.tga");
	//gli::save(TextureMipmaped, 3, "../test3.tga");

	return Error;
}

int test_image_export_dds()
{
	int Error(0);

	{
		gli::texture2D Texture(gli::load_dds("../../data/test_rgb8.dds"));
		assert(!Texture.empty());
		gli::save_dds(Texture, "../../data/test_rgb8_saved.dds");
	}
	{
		gli::texture2D Texture(gli::load_dds("../../data/test_rgb8_saved.dds"));
		assert(!Texture.empty());
		gli::save_dds(Texture, "../../data/test_rgb8_reloaded.dds");
	}
	{
		gli::texture2D Texture(gli::load_dds("../../data/test_dxt1.dds"));
		assert(!Texture.empty());
		gli::save_dds(Texture, "../../data/test_dxt1_saved.dds");
	}
	{
		gli::texture2D Texture(gli::load_dds("../../data/test_bc1.dds"));
		assert(!Texture.empty());
		gli::save_dds(Texture, "../../data/test_bc1_saved.dds");
	}

	return Error;
}

/*
int test_texture1d()
{
	int Error(0);

	{
		gli::texture1D Texture(
			gli::texture1D::size_type(1), 
			gli::RGBA8U, 
			gli::texture1D::dimensions_type(256));

		for(gli::texture1D::size_type TexelIndex = 0; TexelIndex < Texture.size<glm::u8vec4>(); ++TexelIndex)
			*(Texture.data<glm::u8vec4>() + TexelIndex) = glm::u8vec4(255, 128, 0, 255);

		gli::saveDDS10(Texture, "../../data/texture1D_1level_rgba8u_256.dds");
	}

	{
		gli::texture1D Texture(
			gli::texture1D::size_type(glm::log2(256.f)), 
			gli::RGBA8U, 
			gli::texture1D::dimensions_type(256));

		for(gli::texture1D::size_type TexelIndex = 0; TexelIndex < Texture.size<glm::u8vec4>(); ++TexelIndex)
			*(Texture.data<glm::u8vec4>() + TexelIndex) = glm::u8vec4(255, 128, 0, 255);

		gli::saveDDS10(Texture, "../../data/texture1D_mipmaps_rgba8u_256.dds");
	}

	return 0;
}
*/

int test_texture2d_load()
{
	int Error(0);

	{
		gli::texture2D Texture(gli::load_dds("../../data/test_rgb8.dds"));

		gli::texture2D::format_type Format = Texture.format();

	}

	return Error;
}

int test_texture2d()
{
	int Error(0);

	{
		gli::texture2D Texture(
			gli::texture2D::size_type(1), 
			gli::RGBA8_UNORM, 
			gli::texture2D::dimensions_type(256));

		for(gli::texture2D::size_type TexelIndex = 0; TexelIndex < Texture.size<glm::u8vec4>(); ++TexelIndex)
			*(Texture.data<glm::u8vec4>() + TexelIndex) = glm::u8vec4(255, 128, 0, 255);

		gli::save_dds(Texture, "../../data/texture2D_1level_rgba8_unorm_256.dds");
	}

	{
		gli::texture2D Texture(
			gli::texture2D::size_type(1), 
			gli::RGBA16_UNORM, 
			gli::texture2D::dimensions_type(256));

		for(gli::texture2D::size_type TexelIndex = 0; TexelIndex < Texture.size<glm::u16vec4>(); ++TexelIndex)
			*(Texture.data<glm::u16vec4>() + TexelIndex) = glm::u16vec4(65535, 32768, 0, 65535);

		gli::save_dds(Texture, "../../data/texture2D_1level_rgba16_unorm_256.dds");
	}

	{
		gli::texture2D Texture(
			gli::texture2D::size_type(glm::log2(256.f)), 
			gli::RGBA8_UNORM, 
			gli::texture2D::dimensions_type(256));

		for(gli::texture2D::size_type TexelIndex = 0; TexelIndex < Texture.size<glm::u8vec4>(); ++TexelIndex)
			*(Texture.data<glm::u8vec4>() + TexelIndex) = glm::u8vec4(255, 128, 0, 255);

		gli::save_dds(Texture, "../../data/texture2D_mipmaps_rgba8_unorm_256.dds");
	}

	{
		gli::texture2D Texture(
			gli::texture2D::size_type(glm::log2(256.f)), 
			gli::RGBA16_UNORM, 
			gli::texture2D::dimensions_type(256));

		for(gli::texture2D::size_type TexelIndex = 0; TexelIndex < Texture.size<glm::u16vec4>(); ++TexelIndex)
			*(Texture.data<glm::u16vec4>() + TexelIndex) = glm::u16vec4(65535, 32768, 0, 65535);

		gli::save_dds(Texture, "../../data/texture2D_mipmaps_rgba16_unorm_256.dds");
	}

	return 0;
}

int test_texture2DArray()
{
	int Error(0);

	{
		gli::texture2DArray Texture(
			gli::texture2D::size_type(2), 
			gli::texture2D::size_type(1), 
			gli::RGBA8_UNORM, 
			gli::texture2D::dimensions_type(8));

		std::vector<glm::u8vec4> Color;
		Color.push_back(glm::u8vec4(255, 128,   0, 255));
		Color.push_back(glm::u8vec4(  0, 128, 255, 255));

		for(gli::texture2DArray::size_type LayerIndex = 0; LayerIndex < Texture.layers(); ++LayerIndex)
		for(gli::texture2D::size_type TexelIndex = 0; TexelIndex < Texture[LayerIndex].size<glm::u8vec4>(); ++TexelIndex)
			*(Texture[LayerIndex].data<glm::u8vec4>() + TexelIndex) = Color[LayerIndex];

		gli::save_dds(Texture, "../../data/texture2DArray_rgba8u_256.dds");
	}

	{
		gli::texture2D TextureA(gli::load_dds("../../data/test-rgb8-256a.dds"));
		gli::texture2D TextureB(gli::load_dds("../../data/test-rgb8-256b.dds"));

		assert(TextureA.format() == TextureB.format());
		assert(glm::all(glm::equal(TextureA.dimensions(), TextureB.dimensions())));

		gli::texture2DArray Texture(
			gli::texture2D::size_type(2), 
			TextureA.levels(), 
			TextureA.format(), 
			TextureA.dimensions());

		memcpy(Texture[0].data(), TextureA.data(), TextureA.size());
		memcpy(Texture[1].data(), TextureB.data(), TextureB.size());

		gli::texture2D Texture0 = gli::copy(Texture[0]);
		gli::texture2D Texture1 = gli::copy(Texture[1]);

		gli::save_dds(Texture0, "../../data/test-rgb8-256-layer0-saved.dds");
		gli::save_dds(Texture1, "../../data/test-rgb8-256-layer1-saved.dds");
		gli::save_dds(Texture, "../../data/test-rgb8-256-array-saved.dds");
	}

	{
		gli::texture2DArray Texture(gli::load_dds("../../data/test-rgb8-256-array-saved.dds"));
		assert(!Texture.empty());
	}

	return 0;
}

int test_textureCube()
{
	int Error(0);

	{
		gli::textureCube Texture(
			gli::textureCube::size_type(6),
			gli::textureCube::size_type(1), 
			gli::RGBA8_UNORM, 
			gli::textureCube::dimensions_type(256));

		std::vector<glm::u8vec4> Color;
		Color.push_back(glm::u8vec4(255,   0,   0, 255));
		Color.push_back(glm::u8vec4(255, 128,   0, 255));
		Color.push_back(glm::u8vec4(255, 255,   0, 255));
		Color.push_back(glm::u8vec4(  0, 255,   0, 255));
		Color.push_back(glm::u8vec4(  0, 128, 255, 255));
		Color.push_back(glm::u8vec4(  0,   0, 255, 255));

		for(gli::textureCube::size_type FaceIndex = 0; FaceIndex < Texture.faces(); ++FaceIndex)
		for(gli::texture2D::size_type TexelIndex = 0; TexelIndex < Texture[FaceIndex].size<glm::u8vec4>(); ++TexelIndex)
			*(Texture[FaceIndex].data<glm::u8vec4>() + TexelIndex) = Color[FaceIndex];

		gli::save_dds(Texture, "../../data/textureCube_rgba8u_256.dds");
	}

	{
		gli::textureCube Texture(gli::load_dds("../../data/cube.dds"));
		gli::save_dds(Texture, "../../data/cube_saved.dds");
	}

	{
		gli::textureCube Texture(gli::load_dds("../../data/cube_saved.dds"));
	}

	return 0;
}

int main()
{
	int Error(0);

	Error += test_texture2DArray();
	Error += test_textureCube();
	Error += test_texture2d_load();
	Error += test_texture2d();
	Error += test_image_wip();
	Error += test_image_export_dds();

	return Error;
}
