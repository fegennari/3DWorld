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
/// @file gli/core/texture_cube_array.inl
/// @date 2013-01-10 / 2013-01-13
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli
{
	inline textureCubeArray::textureCubeArray() :
		BaseLayer(0),
		MaxLayer(0),
		BaseFace(0),
		MaxFace(0),
		BaseLevel(0),
		MaxLevel(0),
		Format(FORMAT_NULL)
	{}

	inline textureCubeArray::textureCubeArray
	(
		size_type const & Layers,
		size_type const & Faces,
		size_type const & Levels,
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			Layers,
			Faces,
			Levels,
			Format,
			storage::dimensions_type(Dimensions, 1)),
		BaseLayer(0), 
		MaxLayer(Layers - 1), 
		BaseFace(0), 
		MaxFace(Faces - 1), 
		BaseLevel(0), 
		MaxLevel(Levels - 1),
		Format(Format)
	{}

	inline textureCubeArray::textureCubeArray
	(
		size_type const & Layers,
		size_type const & Faces,
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			Layers,
			Faces,
			size_type(glm::log2(int(glm::max(Dimensions.x, Dimensions.y))) + 1),
			Format,
			storage::dimensions_type(Dimensions, 1)),
		BaseLayer(0), 
		MaxLayer(Layers - 1), 
		BaseFace(0), 
		MaxFace(Faces - 1), 
		BaseLevel(0), 
		MaxLevel(glm::log2(int(glm::max(Dimensions.x, Dimensions.y)))),
		Format(Format)
	{}

	inline textureCubeArray::textureCubeArray
	(
		storage const & Storage
	) :
		Storage(Storage),
		BaseLayer(0), 
		MaxLayer(Storage.layers() - 1), 
		BaseFace(0), 
		MaxFace(Storage.faces() - 1), 
		BaseLevel(0), 
		MaxLevel(Storage.levels() - 1),
		Format(Storage.format())
	{}

	inline textureCubeArray::textureCubeArray
	(
		storage const & Storage,
		format_type const & Format,
		size_type BaseLayer,
		size_type MaxLayer,
		size_type BaseFace,
		size_type MaxFace,
		size_type BaseLevel,
		size_type MaxLevel
	) :
		Storage(Storage),
		BaseLayer(BaseLayer),
		MaxLayer(MaxLayer),
		BaseFace(BaseFace),
		MaxFace(MaxFace),
		BaseLevel(BaseLevel),
		MaxLevel(MaxLevel),
		Format(Format)
	{}

	inline textureCubeArray::textureCubeArray
	(
		textureCubeArray const & Texture,
		size_type const & BaseLayer,
		size_type const & MaxLayer,
		size_type const & BaseFace,
		size_type const & MaxFace,
		size_type const & BaseLevel,
		size_type const & MaxLevel
	) :
		Storage(Texture.Storage),
		BaseLayer(Texture.baseLayer() + BaseLayer),
		MaxLayer(Texture.baseLayer() + MaxLayer),
		BaseFace(Texture.baseFace() + BaseFace),
		MaxFace(Texture.baseFace() + MaxFace),
		BaseLevel(Texture.baseLevel() + BaseLevel),
		MaxLevel(Texture.baseLevel() + MaxLevel),
		Format(Texture.format())
	{}

	inline textureCubeArray::textureCubeArray
	(
		textureCube const & Texture,
		size_type const & BaseFace,
		size_type const & MaxFace,
		size_type const & BaseLevel,
		size_type const & MaxLevel
	) :
		Storage(Texture),
		BaseLayer(Texture.baseLayer()),
		MaxLayer(Texture.maxLayer()),
		BaseFace(Texture.baseFace() + BaseFace),
		MaxFace(Texture.baseFace() + MaxFace),
		BaseLevel(Texture.baseLevel() + BaseLevel),
		MaxLevel(Texture.baseLevel() + MaxLevel),
		Format(Texture.format())
	{}

	inline textureCubeArray::textureCubeArray
	(
		texture2D const & Texture,
		size_type const & BaseLevel,
		size_type const & MaxLevel
	) :
		Storage(Texture),
		BaseLayer(Texture.baseLayer()),
		MaxLayer(Texture.maxLayer()),
		BaseFace(Texture.baseFace()),
		MaxFace(Texture.maxFace()),
		BaseLevel(Texture.baseLevel() + BaseLevel),
		MaxLevel(Texture.baseLevel() + MaxLevel),
		Format(Texture.format())
	{}

	inline textureCubeArray::operator storage() const
	{
		return this->Storage;
	}

	inline textureCube textureCubeArray::operator[] 
	(
		size_type const & Layer
	) const
	{
		assert(Layer < this->layers());

		return textureCube(
			this->Storage, this->format(),
			this->baseLayer() + Layer, this->baseLayer() + Layer, 
			this->baseFace(), this->maxFace(),
			this->baseLevel(), this->maxLevel());
	}

	inline bool textureCubeArray::empty() const
	{
		return this->Storage.empty();
	}

	inline textureCubeArray::dimensions_type textureCubeArray::dimensions() const
	{
		return textureCubeArray::dimensions_type(this->Storage.dimensions(this->baseLevel()));
	}

	inline textureCube::format_type textureCubeArray::format() const
	{
		return this->Format;
	}

	inline textureCubeArray::size_type textureCubeArray::layers() const
	{
		return this->maxLayer() - this->baseLayer() + 1;
	}
	
	inline textureCubeArray::size_type textureCubeArray::faces() const
	{
		return this->maxFace() - this->baseFace() + 1;
	}
	
	inline textureCubeArray::size_type textureCubeArray::levels() const
	{
		return this->maxLevel() - this->baseLevel() + 1;
	}

	inline textureCubeArray::size_type textureCubeArray::size() const
	{
		assert(!this->empty());

		return this->Storage.layerSize(
			this->baseFace(), this->maxFace(),
			this->baseLevel(), this->maxLevel()) * this->layers();
	}

	inline void * textureCubeArray::data()
	{
		assert(!this->empty());

		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	inline void const * textureCubeArray::data() const
	{
		assert(!this->empty());
		
		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	template <typename genType>
	inline textureCubeArray::size_type textureCubeArray::size() const
	{
		assert(sizeof(genType) <= this->Storage.blockSize());
		return this->size() / sizeof(genType);
	}

	template <typename genType>
	inline genType * textureCubeArray::data()
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<genType *>(this->data());
	}

	template <typename genType>
	inline genType const * textureCubeArray::data() const
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<genType const *>(this->data());
	}

	inline void textureCubeArray::clear()
	{
		for(size_type Layer = 0; Layer < this->layers(); ++Layer)
		for(size_type Face = 0; Face < this->faces(); ++Face)
			(*this)[Layer][Face].clear();
	}

	template <typename genType>
	inline void textureCubeArray::clear(genType const & Texel)
	{
		assert(this->Storage.blockSize() == sizeof(genType));

		for(size_type Layer = 0; Layer < this->layers(); ++Layer)
		for(size_type Face = 0; Face < this->faces(); ++Face)
			(*this)[Layer][Face].clear<genType>(Texel);
	}

	inline textureCubeArray::size_type textureCubeArray::baseLayer() const
	{
		return this->BaseLayer;
	}

	inline textureCubeArray::size_type textureCubeArray::maxLayer() const
	{
		return this->MaxLayer;
	}

	inline textureCubeArray::size_type textureCubeArray::baseFace() const
	{
		return this->BaseFace;
	}

	inline textureCubeArray::size_type textureCubeArray::maxFace() const
	{
		return this->MaxFace;
	}

	inline textureCubeArray::size_type textureCubeArray::baseLevel() const
	{
		return this->BaseLevel;
	}

	inline textureCubeArray::size_type textureCubeArray::maxLevel() const
	{
		return this->MaxLevel;
	}
}//namespace gli
