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
/// @file gli/core/image.cpp
/// @date 2011-10-07 / 2013-11-25
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <gli/gli.hpp>

int test_image_ctor()
{
	int Error(0);

	gli::image ImageA(
		gli::RGBA8U,
		gli::image::dimensions_type(4, 4, 1));
	gli::image ImageB(
		gli::image::dimensions_type(4, 4, 1), 
		4, 
		gli::image::dimensions_type(1));
	gli::image ImageC = ImageA;

	Error += ImageA == ImageB ? 0 : 1;
	Error += ImageC == ImageB ? 0 : 1;

	return Error;
}

int test_image_data()
{
	int Error(0);

	gli::image ImageA;
	Error += ImageA.empty() ? 0 : 1;
	assert(!Error);

	gli::image ImageB(
		gli::image::dimensions_type(1, 1, 1), 
		sizeof(glm::u8vec4), 
		gli::image::dimensions_type(1, 1, 1));
	Error += ImageB.size() == sizeof(glm::u8vec4) ? 0 : 1;

	*ImageB.data<glm::u8vec4>() = glm::u8vec4(255, 127, 0, 255);
	Error += !ImageB.empty() ? 0 : 1;
	assert(!Error);

	gli::storage Storage(
		2, 1, 1, 
		gli::storage::dimensions3_type(1), 
		gli::FORMAT_NULL,
		sizeof(glm::u8vec4),
		gli::storage::dimensions3_type(1));

	std::vector<glm::u8vec4> Data(2);
	Data[0] = glm::u8vec4(  0, 127, 255, 255);
	Data[1] = glm::u8vec4(255, 127,   0, 255);

	memcpy(Storage.data(), &Data[0][0], Data.size() * sizeof(glm::u8vec4));

	gli::image ImageC(Storage, 0, 0, 0, 0, 0, 1);
	Error += !ImageC.empty() ? 0 : 1;
	assert(!Error);

	glm::u8vec4 ValueB = *(ImageB.data<glm::u8vec4>() + 0);
	glm::u8vec4 ValueC = *(ImageC.data<glm::u8vec4>() + 1);
	glm::u8vec4 ValueD = *(ImageC.data<glm::u8vec4>() + 0);

	Error += glm::all(glm::equal(ValueB, ValueC)) ? 0 : 1;
	assert(!Error);

	return Error;
}

int test_image_query()
{
	int Error(0);

	gli::image Image(gli::image::dimensions_type(1, 1, 1), sizeof(glm::u8vec4), gli::image::dimensions_type(1, 1, 1));

	Error += Image.size() == sizeof(glm::u8vec4) ? 0 : 1;
	Error += !Image.empty() ? 0 : 1;
	Error += Image.dimensions().x == 1 ? 0 : 1;
	Error += Image.dimensions().y == 1 ? 0 : 1;
	Error += Image.dimensions().z == 1 ? 0 : 1;

	return Error;
}

namespace fetch
{
	int test()
	{
		int Error(0);

		gli::image Image(gli::RGBA8U, gli::image::dimensions_type(4, 2, 1));
		*(Image.data<glm::u8vec4>() + 0) = glm::u8vec4(255,   0,   0, 255);
		*(Image.data<glm::u8vec4>() + 1) = glm::u8vec4(255, 128,   0, 255);
		*(Image.data<glm::u8vec4>() + 2) = glm::u8vec4(255, 255,   0, 255);
		*(Image.data<glm::u8vec4>() + 3) = glm::u8vec4(128, 255,   0, 255);
		*(Image.data<glm::u8vec4>() + 4) = glm::u8vec4(  0, 255,   0, 255);
		*(Image.data<glm::u8vec4>() + 5) = glm::u8vec4(  0, 255, 255, 255);
		*(Image.data<glm::u8vec4>() + 6) = glm::u8vec4(  0,   0, 255, 255);
		*(Image.data<glm::u8vec4>() + 7) = glm::u8vec4(255,   0, 255, 255);

		glm::u8vec4 Data0 = Image.load<glm::u8vec4>(gli::image::dimensions_type(0, 0, 0));
		glm::u8vec4 Data1 = Image.load<glm::u8vec4>(gli::image::dimensions_type(1, 0, 0));
		glm::u8vec4 Data2 = Image.load<glm::u8vec4>(gli::image::dimensions_type(2, 0, 0));
		glm::u8vec4 Data3 = Image.load<glm::u8vec4>(gli::image::dimensions_type(3, 0, 0));
		glm::u8vec4 Data4 = Image.load<glm::u8vec4>(gli::image::dimensions_type(0, 1, 0));
		glm::u8vec4 Data5 = Image.load<glm::u8vec4>(gli::image::dimensions_type(1, 1, 0));
		glm::u8vec4 Data6 = Image.load<glm::u8vec4>(gli::image::dimensions_type(2, 1, 0));
		glm::u8vec4 Data7 = Image.load<glm::u8vec4>(gli::image::dimensions_type(3, 1, 0));

		Error += glm::all(glm::equal(Data0, glm::u8vec4(255,   0,   0, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Data1, glm::u8vec4(255, 128,   0, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Data2, glm::u8vec4(255, 255,   0, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Data3, glm::u8vec4(128, 255,   0, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Data4, glm::u8vec4(  0, 255,   0, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Data5, glm::u8vec4(  0, 255, 255, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Data6, glm::u8vec4(  0,   0, 255, 255))) ? 0 : 1;
		Error += glm::all(glm::equal(Data7, glm::u8vec4(255,   0, 255, 255))) ? 0 : 1;

		return Error;
	}
}//namespace fetch

int main()
{
	int Error(0);

	Error += test_image_ctor();
	Error += test_image_data();
	Error += test_image_query();
	Error += fetch::test();
		
	return Error;
}

