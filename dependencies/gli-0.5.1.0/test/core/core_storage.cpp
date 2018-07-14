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
/// @file gli/core/storage.cpp
/// @date 2012-06-25 / 2013-01-13
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <gli/core/storage.hpp>

int test_storage_layer_size()
{
	int Error(0);

	gli::storage Storage(
		2, 1, 1,
		gli::storage::dimensions3_type(2, 2, 1),
		gli::FORMAT_NULL,
		sizeof(glm::u8vec4),
		gli::storage::dimensions3_type(1));

	std::vector<glm::u8vec4> Data(8, glm::u8vec4(0));
	for(std::size_t i = 0; i < 4; ++i)
		Data[i + 0] = glm::u8vec4(255, 127, 0, 255);
	for(std::size_t i = 0; i < 4; ++i)
		Data[i + 4] = glm::u8vec4(0, 127, 255, 255);

	memcpy(Storage.data(), &Data[0][0], Data.size() * sizeof(glm::u8vec4));

	Error += Storage.blockSize() == sizeof(glm::u8vec4) ? 0 : 1;
	Error += Storage.levelSize(0) == sizeof(glm::u8vec4) * 2 * 2 ? 0 : 1;
	Error += Storage.faceSize(0, Storage.levels() - 1) == sizeof(glm::u8vec4) * 2 * 2 ? 0 : 1;
	Error += Storage.layerSize(0, Storage.faces() - 1, 0, Storage.levels() - 1) == sizeof(glm::u8vec4) * 2 * 2 ? 0 : 1;
	Error += Storage.size() == sizeof(glm::u8vec4) * 2 * 2 * 2 ? 0 : 1;

	return Error;
}

int test_storage_face_size()
{
	int Error(0);

	gli::storage Storage(
		1, 6, 1,
		gli::storage::dimensions3_type(2, 2, 1),
		gli::FORMAT_NULL,
		sizeof(glm::u8vec4),
		gli::storage::dimensions3_type(1));

	gli::storage::size_type BlockSize = Storage.blockSize();
	Error += BlockSize == sizeof(glm::u8vec4) ? 0 : 1;

	gli::storage::size_type LevelSize = Storage.levelSize(0);
	Error += LevelSize == sizeof(glm::u8vec4) * 2 * 2 ? 0 : 1;

	gli::storage::size_type FaceSize = Storage.faceSize(0, Storage.levels() - 1);
	Error += FaceSize == sizeof(glm::u8vec4) * 2 * 2 ? 0 : 1;

	gli::storage::size_type LayerSize = Storage.layerSize(0, Storage.faces() - 1, 0, Storage.levels() - 1);
	Error += LayerSize == sizeof(glm::u8vec4) * 2 * 2 * 6 ? 0 : 1;

	gli::storage::size_type Size = Storage.size();
	Error += Size == sizeof(glm::u8vec4) * 2 * 2 * 6 ? 0 : 1;

	return Error;
}

int main()
{
	int Error(0);

	Error += test_storage_layer_size();
    Error += test_storage_face_size();

	assert(!Error);

	return Error;
}
