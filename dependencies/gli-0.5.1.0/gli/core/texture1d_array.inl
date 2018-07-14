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
/// @file gli/core/texture1d_array.inl
/// @date 2012-06-27 / 2013-01-12
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli
{
	inline texture1DArray::texture1DArray() :
		BaseLayer(0), 
		MaxLayer(0), 
		BaseFace(0), 
		MaxFace(0), 
		BaseLevel(0), 
		MaxLevel(0),
		Format(FORMAT_NULL)
	{}

	inline texture1DArray::texture1DArray
	(
		size_type const & Layers,
		size_type const & Levels,
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			Layers, 
			1, 
			Levels,
			Format,
			storage::dimensions_type(Dimensions, 1, 1)),
		BaseLayer(0), 
		MaxLayer(Layers - 1), 
		BaseFace(0), 
		MaxFace(0), 
		BaseLevel(0), 
		MaxLevel(Levels - 1),
		Format(Format)
	{}

	inline texture1DArray::texture1DArray
	(
		size_type const & Layers,
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			Layers,
			1,
			size_type(glm::log2(int(Dimensions)) + 1),
			Format,
			storage::dimensions_type(Dimensions, 1, 1)),
		BaseLayer(0), 
		MaxLayer(Layers - 1), 
		BaseFace(0), 
		MaxFace(0), 
		BaseLevel(0), 
		MaxLevel(glm::log2(int(Dimensions))),
		Format(Format)
	{}

	inline texture1DArray::texture1DArray
	(
		storage const & Storage
	) :
		Storage(Storage),
		BaseLayer(0), 
		MaxLayer(Storage.layers() - 1), 
		BaseFace(0), 
		MaxFace(0), 
		BaseLevel(0), 
		MaxLevel(Storage.levels() - 1),
		Format(Storage.format())
	{}

	inline texture1DArray::texture1DArray
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

	inline texture1DArray::texture1DArray
	(
		texture1DArray const & Texture,
		size_type const & BaseLayer,
		size_type const & MaxLayer,
		size_type const & BaseLevel,
		size_type const & MaxLevel
	) :
		Storage(Texture.Storage),
		BaseLayer(Texture.baseLayer() + BaseLayer),
		MaxLayer(Texture.baseLayer() + MaxLayer),
		BaseFace(Texture.baseFace()),
		MaxFace(Texture.maxFace()),
		BaseLevel(Texture.baseLevel() + BaseLevel),
		MaxLevel(Texture.baseLevel() + MaxLevel),
		Format(Texture.format())
	{}

	inline texture1DArray::texture1DArray
	(
		texture1D const & Texture,
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

	inline texture1DArray::operator storage() const
	{
		return this->Storage;
	}

	inline texture1D texture1DArray::operator[]
	(
		size_type const & Layer
	) const
	{
		assert(!this->empty());
		assert(Layer < this->layers());

		return texture1D(
			this->Storage, this->format(),
			this->baseLayer() + Layer, this->baseLayer() + Layer, 
			this->baseFace(), 	this->maxFace(),
			this->baseLevel(), this->maxLevel());
	}

	inline bool texture1DArray::empty() const
	{
		return this->Storage.empty();
	}

	inline texture1DArray::dimensions_type texture1DArray::dimensions() const
	{
		assert(!this->empty());

		return texture1DArray::dimensions_type(this->Storage.dimensions(this->baseLevel()).x);
	}

	inline texture1DArray::format_type texture1DArray::format() const
	{
		return this->Format;
	}

	inline texture1DArray::size_type texture1DArray::layers() const
	{
		return this->maxLayer() - this->baseLayer() + 1;
	}

	inline texture1DArray::size_type texture1DArray::faces() const
	{
		return this->maxFace() - this->baseFace() + 1;
	}

	inline texture1DArray::size_type texture1DArray::levels() const
	{
		return this->maxLevel() - this->baseLevel() + 1;
	}

	inline texture1DArray::size_type texture1DArray::size() const
	{
		assert(!this->empty());

		return this->Storage.layerSize(
			this->baseFace(), this->maxFace(),
			this->baseLevel(), this->maxLevel()) * this->layers();
	}

	inline void * texture1DArray::data()
	{
		assert(!this->empty());

		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	inline void const * texture1DArray::data() const
	{
		assert(!this->empty());

		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	template <typename genType>
	inline texture1DArray::size_type texture1DArray::size() const
	{
		assert(!this->empty());
		assert(sizeof(genType) <= this->Storage.blockSize());

		return this->size() / sizeof(genType);
	}

	template <typename genType>
	inline genType * texture1DArray::data()
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<genType *>(this->data());
	}

	template <typename genType>
	inline genType const * texture1DArray::data() const
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<const genType *>(this->data());
	}

	inline void texture1DArray::clear()
	{
		for(size_type Layer = 0; Layer < this->layers(); ++Layer)
			(*this)[Layer].clear();
	}

	template <typename genType>
	inline void texture1DArray::clear(genType const & Texel)
	{
		assert(this->Storage.blockSize() == sizeof(genType));

		for(size_type Layer = 0; Layer < this->layers(); ++Layer)
			(*this)[Layer].clear<genType>(Texel);
	}

	inline texture1DArray::size_type texture1DArray::baseLayer() const
	{
		return this->BaseLayer;
	}

	inline texture1DArray::size_type texture1DArray::maxLayer() const
	{
		return this->MaxLayer;
	}

	inline texture1DArray::size_type texture1DArray::baseFace() const
	{
		return this->BaseFace;
	}

	inline texture1DArray::size_type texture1DArray::maxFace() const
	{
		return this->MaxFace;
	}

	inline texture1DArray::size_type texture1DArray::baseLevel() const
	{
		return this->BaseLevel;
	}

	inline texture1DArray::size_type texture1DArray::maxLevel() const
	{
		return this->MaxLevel;
	}
}//namespace gli
