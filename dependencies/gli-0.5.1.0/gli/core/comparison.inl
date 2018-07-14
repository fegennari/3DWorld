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
/// @file gli/core/comparison.inl
/// @date 2013-02-04 / 2013-02-04
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli{
namespace detail
{
	template <typename texture>
	inline bool equalDataValues(texture const & TextureA, texture const & TextureB)
	{
		// Compare the pointer
		if(TextureA.data() == TextureB.data())
			return true;

		// Compare the actual data
        for(typename texture::size_type Level = 0; Level < TextureA.levels(); ++Level)
        for(typename texture::size_type i = 0; i < TextureA[Level].template size<glm::byte>(); ++i)
        {
            glm::byte A = *(TextureA[Level].template data<glm::byte>() + i);
            glm::byte B = *(TextureB[Level].template data<glm::byte>() + i);
            if(A != B)
                return false;
        }

		return true;
	}

	inline bool equalData(texture1D const & TextureA, texture1D const & TextureB)
	{
		return equalDataValues(TextureA, TextureB);
	}

	inline bool equalData(texture1DArray const & TextureA, texture1DArray const & TextureB)
	{
		for(std::size_t Layer = 0; Layer < TextureA.layers(); ++Layer)
			if(!equalDataValues(TextureA[Layer], TextureB[Layer]))
				return false;
		return true;
	}

	inline bool equalData(texture2D const & TextureA, texture2D const & TextureB)
	{
		return equalDataValues(TextureA, TextureB);
	}

	inline bool equalData(texture2DArray const & TextureA, texture2DArray const & TextureB)
	{
		for(std::size_t Layer = 0; Layer < TextureA.layers(); ++Layer)
			if(!equalDataValues(TextureA[Layer], TextureB[Layer]))
				return false;
		return true;
	}

	inline bool equalData(texture3D const & TextureA, texture3D const & TextureB)
	{
		return equalDataValues(TextureA, TextureB);
	}

	inline bool equalData(textureCube const & TextureA, textureCube const & TextureB)
	{
		for(std::size_t Face = 0; Face < TextureA.faces(); ++Face)
			if(!equalDataValues(TextureA[Face], TextureB[Face]))
				return false;
		return true;
	}

	inline bool equalData(textureCubeArray const & TextureA, textureCubeArray const & TextureB)
	{
		for(std::size_t Layer = 0; Layer < TextureA.layers(); ++Layer)
			for(std::size_t Face = 0; Face < TextureA[Layer].faces(); ++Face)
				if(!equalDataValues(TextureA[Layer][Face], TextureB[Layer][Face]))
					return false;
		return true;
	}

	inline bool notEqualData(texture1D const & TextureA, texture1D const & TextureB)
	{
		return !equalDataValues(TextureA, TextureB);
	}

	inline bool notEqualData(texture1DArray const & TextureA, texture1DArray const & TextureB)
	{
		for(std::size_t Layer = 0; Layer < TextureA.layers(); ++Layer)
			if(!equalDataValues(TextureA[Layer], TextureB[Layer]))
				return true;
		return false;
	}

	inline bool notEqualData(texture2D const & TextureA, texture2D const & TextureB)
	{
		return !equalDataValues(TextureA, TextureB);
	}

	inline bool notEqualData(texture2DArray const & TextureA, texture2DArray const & TextureB)
	{
		for(std::size_t Layer = 0; Layer < TextureA.layers(); ++Layer)
			if(!equalDataValues(TextureA[Layer], TextureB[Layer]))
				return true;
		return false;
	}

	inline bool notEqualData(texture3D const & TextureA, texture3D const & TextureB)
	{
		return !equalDataValues(TextureA, TextureB);
	}

	inline bool notEqualData(textureCube const & TextureA, textureCube const & TextureB)
	{
		for(std::size_t Face = 0; Face < TextureA.faces(); ++Face)
			if(!equalDataValues(TextureA[Face], TextureB[Face]))
				return true;
		return false;
	}

	inline bool notEqualData(textureCubeArray const & TextureA, textureCubeArray const & TextureB)
	{
		for(std::size_t Layer = 0; Layer < TextureA.layers(); ++Layer)
			for(std::size_t Face = 0; Face < TextureA[Layer].faces(); ++Face)
				if(!equalDataValues(TextureA[Layer][Face], TextureB[Layer][Face]))
					return true;
		return false;
	}

}//namespace detail

	inline bool operator==(image const & ImageA, image const & ImageB)
	{
		if(!glm::all(glm::equal(ImageA.dimensions(), ImageB.dimensions())))
			return false;
		if(ImageA.size() != ImageB.size())
			return false;

		for(std::size_t i = 0; i < ImageA.size<glm::byte>(); ++i)
			if(*(ImageA.data<glm::byte>() + i) != *(ImageB.data<glm::byte>() + i))
				return false;

		return true;
	}

	inline bool operator!=(image const & ImageA, image const & ImageB)
	{
		if(!glm::all(glm::equal(ImageA.dimensions(), ImageB.dimensions())))
			return true;
		if(ImageA.size() != ImageB.size())
			return true;

		for(std::size_t i = 0; i < ImageA.size<glm::byte>(); ++i)
			if(*(ImageA.data<glm::byte>() + i) != *(ImageB.data<glm::byte>() + i))
				return true;

		return false;
	}

	template <typename texture>
	inline bool equal(texture const & TextureA, texture const & TextureB)
	{
		if(TextureA.empty() && TextureB.empty())
			return true;
		if(TextureA.empty() != TextureB.empty())
			return false;
		if(!glm::all(glm::equal(TextureA.dimensions(), TextureB.dimensions())))
			return false;
		if(TextureA.layers() != TextureB.layers())
			return false;
		if(TextureA.faces() != TextureB.faces())
			return false;
		if(TextureA.levels() != TextureB.levels())
			return false;
		if(TextureA.format() != TextureB.format())
			return false;
		if(TextureA.size() != TextureB.size())
			return false;

		return detail::equalData(TextureA, TextureB);
	}

	template <typename texture>
	inline bool notEqual(texture const & TextureA, texture const & TextureB)
	{
		if(TextureA.empty() && TextureB.empty())
			return false;
		if(TextureA.empty() != TextureB.empty())
			return true;
		if(!glm::all(glm::equal(TextureA.dimensions(), TextureB.dimensions())))
			return true;
		if(TextureA.layers() != TextureB.layers())
			return true;
		if(TextureA.faces() != TextureB.faces())
			return true;
		if(TextureA.levels() != TextureB.levels())
			return true;
		if(TextureA.format() != TextureB.format())
			return true;
		if(TextureA.size() != TextureB.size())
			return true;

		return detail::notEqualData(TextureA, TextureB);
	}

	inline bool operator==(gli::texture1D const & A, gli::texture1D const & B)
	{
		return gli::equal(A, B);
	}

	inline bool operator==(gli::texture1DArray const & A, gli::texture1DArray const & B)
	{
		return gli::equal(A, B);
	}

	inline bool operator==(gli::texture2D const & A, gli::texture2D const & B)
	{
		return gli::equal(A, B);
	}

	inline bool operator==(gli::texture2DArray const & A, gli::texture2DArray const & B)
	{
		return gli::equal(A, B);
	}

	inline bool operator==(gli::texture3D const & A, gli::texture3D const & B)
	{
		return gli::equal(A, B);
	}

	inline bool operator==(gli::textureCube const & A, gli::textureCube const & B)
	{
		return gli::equal(A, B);
	}

	inline bool operator==(gli::textureCubeArray const & A, gli::textureCubeArray const & B)
	{
		return gli::equal(A, B);
	}

	inline bool operator!=(gli::texture1D const & A, gli::texture1D const & B)
	{
		return gli::notEqual(A, B);
	}

	inline bool operator!=(gli::texture1DArray const & A, gli::texture1DArray const & B)
	{
		return gli::notEqual(A, B);
	}

	inline bool operator!=(gli::texture2D const & A, gli::texture2D const & B)
	{
		return gli::notEqual(A, B);
	}

	inline bool operator!=(gli::texture2DArray const & A, gli::texture2DArray const & B)
	{
		return gli::notEqual(A, B);
	}

	inline bool operator!=(gli::texture3D const & A, gli::texture3D const & B)
	{
		return gli::notEqual(A, B);
	}

	inline bool operator!=(gli::textureCube const & A, gli::textureCube const & B)
	{
		return gli::notEqual(A, B);
	}

	inline bool operator!=(gli::textureCubeArray const & A, gli::textureCubeArray const & B)
	{
		return gli::notEqual(A, B);
	}
}//namespace gli
