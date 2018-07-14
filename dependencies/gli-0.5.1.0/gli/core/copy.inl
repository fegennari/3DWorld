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
/// @file gli/core/copy.inl
/// @date 2013-01-23 / 2013-02-03
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli{

inline storage copy(storage const & Storage)
{
	storage Result(
		Storage.layers(),
		Storage.faces(),
		Storage.levels(),
		Storage.format(),
		Storage.dimensions(0));

	memcpy(Result.data(), Storage.data(), Storage.size());
		
	return Result;
}

/*
inline image copy(image const & Image)
{
	image Result;

	memcpy(Result.data<glm::byte>(), Image.data<glm::byte>(), Image.size());
		
	return Result;
}
*/

template <>
inline texture1D copy(texture1D const & Texture)
{
	texture1D Copy(
		Texture.levels(), 
		Texture.format(), 
		Texture.dimensions());

	memcpy(
		Copy.data<glm::byte>(), 
		Texture.data<glm::byte>(), 
		Copy.size<glm::byte>());
		
	return Copy;
}

inline texture1D copy
(
	texture1D const & Texture,
	texture1D::size_type const & BaseLevel,
	texture1D::size_type const & MaxLevel
)
{
	assert(BaseLevel <= MaxLevel);
	assert(BaseLevel < Texture.levels());
	assert(MaxLevel < Texture.levels());
	
	texture1D Copy(
		MaxLevel - BaseLevel + 1, 
		Texture.format(), 
		texture1D::dimensions_type(Texture[BaseLevel].dimensions().x));

	memcpy(
		Copy.data<glm::byte>(), 
		Texture[BaseLevel].data<glm::byte>(), 
		Copy.size<glm::byte>());
		
	return Copy;
}

template <>
inline texture1DArray copy(texture1DArray const & Texture)
{
	texture1DArray Copy(
		Texture.layers(),
		Texture.levels(), 
		Texture.format(), 
		Texture.dimensions());

	for(texture1DArray::size_type Layer = 0; Layer < Copy.layers(); ++Layer)
	{
		memcpy(
			Copy[Layer].data<glm::byte>(), 
			Texture[Layer].data<glm::byte>(), 
			Copy[Layer].size<glm::byte>());
	}
		
	return Copy;
}

inline texture1DArray copy
(
	texture1DArray const & Texture,
	texture1DArray::size_type const & BaseLayer,
	texture1DArray::size_type const & MaxMayer,
	texture1DArray::size_type const & BaseLevel,
	texture1DArray::size_type const & MaxLevel
)
{
	assert(BaseLevel <= MaxLevel);
	assert(BaseLevel < Texture.levels());
	assert(MaxLevel < Texture.levels());
	assert(BaseLayer <= MaxMayer);
	assert(BaseLayer < Texture.layers());
	assert(MaxMayer < Texture.layers());

	texture1DArray Copy(
		MaxMayer - BaseLayer + 1, 
		MaxLevel - BaseLevel + 1, 
		Texture.format(), 
		texture1DArray::dimensions_type(Texture[BaseLayer][BaseLevel].dimensions().x));

	for(texture1DArray::size_type Layer = 0; Layer < Copy.layers(); ++Layer)
	{
		memcpy(
			Copy[Layer].data<glm::byte>(), 
			Texture[Layer + BaseLayer][BaseLevel].data<glm::byte>(), 
			Copy[Layer].size<glm::byte>());
	}

	return Copy;
}

template <>
inline texture2D copy(texture2D const & Texture)
{
	texture2D Copy(
		Texture.levels(), 
		Texture.format(), 
		Texture.dimensions());

	memcpy(
		Copy.data<glm::byte>(), 
		Texture.data<glm::byte>(), 
		Copy.size<glm::byte>());
		
	return Copy;
}

inline texture2D copy
(
	texture2D const & Texture,
	texture2D::size_type const & BaseLevel,
	texture2D::size_type const & MaxLevel
)
{
	assert(BaseLevel <= MaxLevel);
	assert(BaseLevel < Texture.levels());
	assert(MaxLevel < Texture.levels());
	
	texture2D Copy(
		MaxLevel - BaseLevel + 1, 
		Texture.format(), 
		texture2D::dimensions_type(Texture[BaseLevel].dimensions().x));

	memcpy(
		Copy.data<glm::byte>(), 
		Texture[BaseLevel].data<glm::byte>(), 
		Copy.size<glm::byte>());
		
	return Copy;
}

template <>
inline texture2DArray copy(texture2DArray const & Texture)
{
	texture2DArray Copy(
		Texture.layers(),
		Texture.levels(), 
		Texture.format(), 
		Texture.dimensions());

	for(texture2DArray::size_type Layer = 0; Layer < Copy.layers(); ++Layer)
	{
		memcpy(
			Copy[Layer].data<glm::byte>(), 
			Texture[Layer].data<glm::byte>(), 
			Copy[Layer].size<glm::byte>());
	}
		
	return Copy;
}

inline texture2DArray copy
(
	texture2DArray const & Texture,
	texture2DArray::size_type const & BaseLayer,
	texture2DArray::size_type const & MaxMayer,
	texture2DArray::size_type const & BaseLevel,
	texture2DArray::size_type const & MaxLevel
)
{
	assert(BaseLevel <= MaxLevel);
	assert(BaseLevel < Texture.levels());
	assert(MaxLevel < Texture.levels());
	assert(BaseLayer <= MaxMayer);
	assert(BaseLayer < Texture.layers());
	assert(MaxMayer < Texture.layers());

	texture2DArray Copy(
		MaxMayer - BaseLayer + 1, 
		MaxLevel - BaseLevel + 1, 
		Texture.format(), 
		texture2DArray::dimensions_type(Texture[BaseLayer][BaseLevel].dimensions()));

	for(texture2DArray::size_type Layer = 0; Layer < Copy.layers(); ++Layer)
	{
		memcpy(
			Copy[Layer].data<glm::byte>(), 
			Texture[Layer + BaseLayer][BaseLevel].data<glm::byte>(), 
			Copy[Layer].size<glm::byte>());
	}
		
	return Copy;
}

template <>
inline texture3D copy(texture3D const & Texture)
{
	texture3D Copy(
		Texture.levels(), 
		Texture.format(), 
		Texture.dimensions());

	memcpy(
		Copy.data<glm::byte>(), 
		Texture.data<glm::byte>(), 
		Texture.size<glm::byte>());
		
	return Copy;
}

inline texture3D copy
(
	texture3D const & Texture,
	texture3D::size_type const & BaseLevel,
	texture3D::size_type const & MaxLevel
)
{
	assert(BaseLevel <= MaxLevel);
	assert(BaseLevel < Texture.levels());
	assert(MaxLevel < Texture.levels());
	
	texture3D Copy(
		MaxLevel - BaseLevel + 1, 
		Texture.format(), 
		texture3D::dimensions_type(Texture[BaseLevel].dimensions()));

	memcpy(
		Copy.data<glm::byte>(), 
		Texture[BaseLevel].data<glm::byte>(), 
		Copy.size<glm::byte>());
		
	return Copy;
}

template <>
inline textureCube copy(textureCube const & Texture)
{
	textureCube Copy(
		Texture.faces(),
		Texture.levels(), 
		Texture.format(), 
		Texture.dimensions());

	for(textureCube::size_type Face = 0; Face < Copy.faces(); ++Face)
	{
		memcpy(
			Copy[Face].data<glm::byte>(), 
			Texture[Face].data<glm::byte>(), 
			Texture[Face].size<glm::byte>());
	}
		
	return Copy;
}

inline textureCube copy
(
	textureCube const & Texture,
	textureCube::size_type const & BaseFace,
	textureCube::size_type const & MaxFace,
	textureCube::size_type const & BaseLevel,
	textureCube::size_type const & MaxLevel
)
{
	assert(BaseLevel <= MaxLevel);
	assert(BaseLevel < Texture.levels());
	assert(MaxLevel < Texture.levels());
	assert(BaseFace <= MaxFace);
	assert(BaseFace < Texture.faces());
	assert(MaxFace < Texture.faces());

	textureCube Copy(
		MaxFace - BaseFace + 1, 
		MaxLevel - BaseLevel + 1, 
		Texture.format(), 
		textureCube::dimensions_type(Texture[BaseFace][BaseLevel].dimensions()));

	for(textureCube::size_type Face = 0; Face < Copy.faces(); ++Face)
	{
		memcpy(
			Copy[Face].data<glm::byte>(), 
			Texture[Face + BaseFace][BaseLevel].data<glm::byte>(), 
			Copy[Face].size<glm::byte>());
	}
		
	return Copy;
}

template <>
inline textureCubeArray copy(textureCubeArray const & Texture)
{
	textureCubeArray Copy(
		Texture.layers(),
		Texture.faces(),
		Texture.levels(), 
		Texture.format(), 
		Texture.dimensions());

	for(textureCubeArray::size_type Layer = 0; Layer < Copy.layers(); ++Layer)
	for(textureCubeArray::size_type Face = 0; Face < Copy[Layer].faces(); ++Face)
	{
		memcpy(
			Copy[Layer][Face].data<glm::byte>(), 
			Texture[Layer][Face].data<glm::byte>(), 
			Copy[Layer][Face].size<glm::byte>());
	}

	return Copy;
}

inline textureCubeArray copy
(
	textureCubeArray const & Texture,
	textureCubeArray::size_type const & BaseLayer,
	textureCubeArray::size_type const & MaxLayer,
	textureCubeArray::size_type const & BaseFace,
	textureCubeArray::size_type const & MaxFace,
	textureCubeArray::size_type const & BaseLevel,
	textureCubeArray::size_type const & MaxLevel
)
{
	assert(BaseLevel <= MaxLevel);
	assert(BaseLevel < Texture.levels());
	assert(MaxLevel < Texture.levels());
	assert(BaseFace <= MaxFace);
	assert(BaseFace < Texture.faces());
	assert(MaxFace < Texture.faces());
	assert(BaseLayer <= MaxLayer);
	assert(BaseLayer < Texture.layers());
	assert(MaxLayer < Texture.layers());

	textureCubeArray Copy(
		MaxLayer - BaseLayer + 1, 
		MaxFace - BaseFace + 1, 
		MaxLevel - BaseLevel + 1, 
		Texture.format(), 
		textureCube::dimensions_type(Texture[BaseLayer][BaseFace][BaseLevel].dimensions()));

	for(textureCubeArray::size_type Layer = 0; Layer < Copy.layers(); ++Layer)
	for(textureCubeArray::size_type Face = 0; Face < Copy[Layer].faces(); ++Face)
	{
		memcpy(
			Copy[Layer][Face].data<glm::byte>(), 
			Texture[Layer + BaseLayer][Face + BaseFace][BaseLevel].data<glm::byte>(), 
			Copy[Layer][Face].size<glm::byte>());
	}

	return Copy;
}

}//namespace gli
