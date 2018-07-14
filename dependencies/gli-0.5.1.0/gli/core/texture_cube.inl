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
/// @file gli/core/texture_cube.inl
/// @date 2011-04-06 / 2012-12-12
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli
{
	inline textureCube::textureCube() :
		BaseLayer(0),
		MaxLayer(0),
		BaseFace(0),
		MaxFace(0),
		BaseLevel(0),
		MaxLevel(0),
		Format(FORMAT_NULL)
	{}

	inline textureCube::textureCube
	(
		size_type const & Faces,
		size_type const & Levels,
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			1,
			Faces,
			Levels,
			Format,
			storage::dimensions_type(Dimensions, 1)),
		BaseLayer(0), 
		MaxLayer(0), 
		BaseFace(0), 
		MaxFace(Faces - 1), 
		BaseLevel(0), 
		MaxLevel(Levels - 1),
		Format(Format)
	{}

	inline textureCube::textureCube
	(
		size_type const & Faces,
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			1,
			Faces,
			size_type(glm::log2(int(glm::max(Dimensions.x, Dimensions.y))) + 1),
			Format,
			storage::dimensions_type(Dimensions, 1)),
		BaseLayer(0), 
		MaxLayer(0), 
		BaseFace(0), 
		MaxFace(Faces - 1), 
		BaseLevel(0), 
		MaxLevel(glm::log2(int(glm::max(Dimensions.x, Dimensions.y)))),
		Format(Format)
	{}

	inline textureCube::textureCube
	(
		storage const & Storage
	) :
		Storage(Storage),
		BaseLayer(0), 
		MaxLayer(0), 
		BaseFace(0), 
		MaxFace(Storage.faces() - 1), 
		BaseLevel(0), 
		MaxLevel(Storage.levels() - 1),
		Format(Storage.format())
	{}

	inline textureCube::textureCube
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

	inline textureCube::textureCube
	(
		textureCube const & Texture,
		size_type const & BaseFace,
		size_type const & MaxFace,
		size_type const & BaseLevel,
		size_type const & MaxLevel
	) :
		Storage(Texture.Storage),
		BaseLayer(Texture.baseLayer()),
		MaxLayer(Texture.maxLayer()),
		BaseFace(Texture.baseFace() + BaseFace),
		MaxFace(Texture.baseFace() + MaxFace),
		BaseLevel(Texture.baseLevel() + BaseLevel),
		MaxLevel(Texture.baseLevel() + MaxLevel),
		Format(Texture.format())
	{}

	inline textureCube::textureCube
	(
		textureCubeArray const & Texture,
		size_type const & BaseLayer,
		size_type const & BaseFace,
		size_type const & MaxFace,
		size_type const & BaseLevel,
		size_type const & MaxLevel
	) :
		Storage(Texture),
		BaseLayer(Texture.baseLayer() + BaseLayer),
		MaxLayer(Texture.baseLayer() + BaseLayer),
		BaseFace(Texture.baseFace() + BaseFace),
		MaxFace(Texture.baseFace() + MaxFace),
		BaseLevel(Texture.baseLevel() + BaseLevel),
		MaxLevel(Texture.baseLevel() + MaxLevel),
		Format(Texture.format())
	{}

	inline textureCube::textureCube
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

	inline textureCube::operator storage() const
	{
		return this->Storage;
	}

	inline texture2D textureCube::operator[] (size_type const & Face) const
	{
		assert(Face < this->faces());

		return texture2D(
			this->Storage, this->format(),
			this->baseLayer(), this->maxLayer(), 
			this->baseFace() + Face, 	this->baseFace() + Face,
			this->baseLevel(), this->maxLevel());
	}

	inline bool textureCube::empty() const
	{
		return this->Storage.empty();
	}

	inline textureCube::dimensions_type textureCube::dimensions() const
	{
		return textureCube::dimensions_type(this->Storage.dimensions(this->baseLevel()));
	}

	inline textureCube::format_type textureCube::format() const
	{
		return this->Format;
	}

	inline textureCube::size_type textureCube::layers() const
	{
		assert(this->maxLayer() - this->baseLayer() + 1 == 1);
		return 1;
	}
	
	inline textureCube::size_type textureCube::faces() const
	{
		return this->maxFace() - this->baseFace() + 1;
	}
	
	inline textureCube::size_type textureCube::levels() const
	{
		return this->maxLevel() - this->baseLevel() + 1;
	}

	inline textureCube::size_type textureCube::size() const
	{
		assert(!this->empty());

		return this->Storage.layerSize(
			this->baseFace(), this->maxFace(),
			this->baseLevel(), this->maxLevel());
	}

	inline void * textureCube::data()
	{
		assert(!this->empty());

		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	inline void const * textureCube::data() const
	{
		assert(!this->empty());
		
		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	template <typename genType>
	inline textureCube::size_type textureCube::size() const
	{
		assert(!this->empty());
		assert(sizeof(genType) <= this->Storage.blockSize());

		return this->size() / sizeof(genType);
	}

	template <typename genType>
	inline genType * textureCube::data()
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<genType *>(this->data());
	}

	template <typename genType>
	inline genType const * textureCube::data() const
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<genType const *>(this->data());
	}

	inline void textureCube::clear()
	{
		for(size_type Face = 0; Face < this->faces(); ++Face)
			(*this)[Face].clear();
	}

	template <typename genType>
	inline void textureCube::clear(genType const & Texel)
	{
		assert(this->Storage.blockSize() == sizeof(genType));

		for(size_type Face = 0; Face < this->faces(); ++Face)
			(*this)[Face].clear<genType>(Texel);
	}

	inline textureCube::size_type textureCube::baseLayer() const
	{
		return this->BaseLayer;
	}

	inline textureCube::size_type textureCube::maxLayer() const
	{
		return this->MaxLayer;
	}

	inline textureCube::size_type textureCube::baseFace() const
	{
		return this->BaseFace;
	}

	inline textureCube::size_type textureCube::maxFace() const
	{
		return this->MaxFace;
	}

	inline textureCube::size_type textureCube::baseLevel() const
	{
		return this->BaseLevel;
	}

	inline textureCube::size_type textureCube::maxLevel() const
	{
		return this->MaxLevel;
	}
}//namespace gli
