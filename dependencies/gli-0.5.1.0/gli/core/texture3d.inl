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
/// @file gli/core/texture3d.inl
/// @date 2013-01-13 / 2013-01-13
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli
{
	inline texture3D::texture3D() :
		BaseLayer(0),
		MaxLayer(0),
		BaseFace(0),
		MaxFace(0),
		BaseLevel(0),
		MaxLevel(0),
		Format(FORMAT_NULL)
	{}

	inline texture3D::texture3D
	(
		size_type const & Levels,
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			1,
			1,
			Levels,
			Format,
			storage::dimensions_type(Dimensions)),
		BaseLayer(0), 
		MaxLayer(0), 
		BaseFace(0), 
		MaxFace(0), 
		BaseLevel(0), 
		MaxLevel(Levels - 1),
		Format(Format)
	{}

	inline texture3D::texture3D
	(
		format_type const & Format,
		dimensions_type const & Dimensions
	) :
		Storage(
			1,
			1,
			size_type(glm::log2(int(glm::max(glm::max(Dimensions.x, Dimensions.y), Dimensions.z))) + 1),
			Format,
			storage::dimensions_type(Dimensions)),
		BaseLayer(0),
		MaxLayer(0),
		BaseFace(0),
		MaxFace(0),
		BaseLevel(0),
		MaxLevel(glm::log2(int(glm::max(glm::max(Dimensions.x, Dimensions.y), Dimensions.z)))),
		Format(Format)
	{}

	inline texture3D::texture3D
	(
		storage const & Storage
	) :
		Storage(Storage),
		BaseLayer(0), 
		MaxLayer(0), 
		BaseFace(0), 
		MaxFace(0), 
		BaseLevel(0), 
		MaxLevel(Storage.levels() - 1),
		Format(Storage.format())
	{}

	inline texture3D::texture3D
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
 
	inline texture3D::texture3D
	(
		texture3D const & Texture,
		size_type const & BaseLevel,
		size_type const & MaxLevel
	) :
		Storage(Texture.Storage),
		BaseLayer(Texture.baseLayer()),
		MaxLayer(Texture.maxLayer()),
		BaseFace(Texture.baseFace()),
		MaxFace(Texture.maxFace()),
		BaseLevel(Texture.baseLevel() + BaseLevel),
		MaxLevel(Texture.baseLevel() + MaxLevel),
		Format(Texture.format())
	{}

	inline texture3D::operator storage() const
	{
		return this->Storage;
	}

	inline image texture3D::operator[]
	(
		texture3D::size_type const & Level
	) const
	{
		assert(Level < this->levels());

		return image(
			this->Storage,
			this->baseLayer(),
			this->maxLayer(),
			this->baseFace(),
			this->maxFace(),
			this->baseLevel() + Level,
			this->baseLevel() + Level);
	}

	inline bool texture3D::empty() const
	{
		return this->Storage.empty();
	}

	inline texture3D::dimensions_type texture3D::dimensions() const
	{
		return texture3D::dimensions_type(this->Storage.dimensions(this->baseLevel()));
	}

	inline texture3D::format_type texture3D::format() const
	{
		return this->Format;
	}

	inline texture3D::size_type texture3D::layers() const
	{
		return 1;
	}

	inline texture3D::size_type texture3D::faces() const
	{
		return 1;
	}

	inline texture3D::size_type texture3D::levels() const
	{
		return this->maxLevel() - this->baseLevel() + 1;
	}

	inline texture3D::size_type texture3D::size() const
	{
		assert(!this->empty());

		return this->Storage.faceSize(this->baseLevel(), this->maxLevel());
	}

	inline void * texture3D::data()
	{
		assert(!this->empty());

		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	inline void const * texture3D::data() const
	{
		assert(!this->empty());
		
		size_type const offset = detail::imageAddressing(
			this->Storage, this->baseLayer(), this->baseFace(), this->baseLevel());

		return this->Storage.data() + offset;
	}

	template <typename genType>
	inline texture3D::size_type texture3D::size() const
	{
		assert(!this->empty());
		assert(sizeof(genType) <= this->Storage.blockSize());

		return this->size() / sizeof(genType);
	}

	template <typename genType>
	inline genType * texture3D::data()
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<genType *>(this->data());
	}

	template <typename genType>
	inline genType const * texture3D::data() const
	{
		assert(!this->empty());
		assert(this->Storage.blockSize() >= sizeof(genType));

		return reinterpret_cast<genType const *>(this->data());
	}

	inline void texture3D::clear()
	{
		memset(this->data<glm::byte>(), 0, this->size<glm::byte>());
	}

	template <typename genType>
	inline void texture3D::clear(genType const & Texel)
	{
		assert(this->Storage.blockSize() == sizeof(genType));

		for(size_type TexelIndex = 0; TexelIndex < this->size<genType>(); ++TexelIndex)
			*(this->data<genType>() + TexelIndex) = Texel;
	}

	inline texture3D::size_type texture3D::baseLayer() const
	{
		return this->BaseLayer;
	}

	inline texture3D::size_type texture3D::maxLayer() const
	{
		return this->MaxLayer;
	}

	inline texture3D::size_type texture3D::baseFace() const
	{
		return this->BaseFace;
	}

	inline texture3D::size_type texture3D::maxFace() const
	{
		return this->MaxFace;
	}

	inline texture3D::size_type texture3D::baseLevel() const
	{
		return this->BaseLevel;
	}

	inline texture3D::size_type texture3D::maxLevel() const
	{
		return this->MaxLevel;
	}
}//namespace gli
