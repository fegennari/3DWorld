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
/// @file gli/core/storage.hpp
/// @date 2012-06-21 / 2013-01-12
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#ifndef GLI_CORE_STORAGE_INCLUDED
#define GLI_CORE_STORAGE_INCLUDED

// STD
#include <vector>
#include <queue>
#include <string>
#include <cassert>
#include <cmath>
#include <cstring>

// GLM
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/gradient_paint.hpp>
#include <glm/gtx/integer.hpp>
#include <glm/gtx/bit.hpp>
#include <glm/gtx/multiple.hpp>
#include <glm/gtx/number_precision.hpp>
#include <glm/gtx/raw_data.hpp>
#include <glm/gtx/scalar_relational.hpp>

#include "shared_ptr.hpp"
#include "header.hpp"
#include "format.hpp"

namespace gli{
namespace detail
{
	typedef std::size_t size_type;
	typedef glm::uint dimensions1_type;
	typedef glm::uvec2 dimensions2_type;
	typedef glm::uvec3 dimensions3_type;
	typedef glm::uvec4 dimensions4_type;
}//namespace detail

	class storage
	{
	public:
		typedef glm::uint dimensions1_type;
		typedef glm::uvec2 dimensions2_type;
		typedef glm::uvec3 dimensions3_type;
		typedef glm::uvec4 dimensions4_type;
		typedef dimensions3_type dimensions_type;
		typedef float texcoord1_type;
		typedef glm::vec2 texcoord2_type;
		typedef glm::vec3 texcoord3_type;
		typedef glm::vec4 texcoord4_type;
		typedef std::size_t size_type;
		typedef gli::format format_type;

	public:
		storage();

		storage(
			size_type const & Layers,
			size_type const & Faces,
			size_type const & Levels,
			format_type const & Format,
			dimensions_type const & Dimensions);

		storage(
			size_type const & Layers,
			size_type const & Faces,
			size_type const & Levels,
			dimensions_type const & Dimensions,
			format_type const & Format,
			size_type const & BlockSize,
			dimensions_type const & BlockDimensions);

		bool empty() const;
		size_type size() const; // Express is bytes
		format_type format() const;
		size_type layers() const;
		size_type faces() const;
		size_type levels() const;

		size_type blockSize() const; // Express is bytes
		dimensions_type blockDimensions() const; // Express is bytes
		dimensions_type dimensions(size_type const & Level) const;

		glm::byte * data();
		glm::byte const * data() const;

		size_type levelSize(
			size_type const & Level) const;
		size_type faceSize(
			size_type const & BaseLevel,
			size_type const & MaxLevel) const;
		size_type layerSize(
			size_type const & BaseFace,
			size_type const & MaxFace,
			size_type const & BaseLevel,
			size_type const & MaxLevel) const;

	private:
		struct impl
		{
			impl();

			explicit impl(
				size_type const & Layers, 
				size_type const & Faces,
				size_type const & Levels,
				format_type const & Format,
				dimensions_type const & Dimensions,
				size_type const & BlockSize,
				dimensions_type const & BlockDimensions);

			size_type const Layers; 
			size_type const Faces;
			size_type const Levels;
			format_type const Format;
			dimensions_type const Dimensions;
			size_type const BlockSize;
			dimensions_type const BlockDimensions;
			std::vector<glm::byte> Data;
		};

		shared_ptr<impl> Impl;
	};

/*
	storage extractLayers(
		storage const & Storage, 
		storage::size_type const & Offset, 
		storage::size_type const & Size);
*/
/*
	storage extractFace(
		storage const & Storage, 
		face const & Face);
*/
/*
	storage extractLevels(
		storage const & Storage, 
		storage::size_type const & Offset, 
		storage::size_type const & Size);
*/
/*
	void copy_layers(
		storage const & SourceStorage, 
		storage::size_type const & SourceLayerOffset,
		storage::size_type const & SourceLayerSize,
		storage & DestinationStorage, 
		storage::size_type const & DestinationLayerOffset);

	void copy_faces(
		storage const & SourceStorage, 
		storage::size_type const & SourceLayerOffset,
		storage::size_type const & SourceFaceOffset,
		storage::size_type const & SourceLayerSize,
		storage & DestinationStorage, 
		storage::size_type const & DestinationLayerOffset,
		storage::size_type const & DestinationFaceOffset);

	void copy_levels(
		storage const & SourceStorage, 
		storage::size_type const & SourceLayerOffset,
		storage::size_type const & SourceFaceOffset,
		storage::size_type const & SourceLevelOffset,
		storage::size_type const & SourceLayerSize,
		storage & DestinationStorage, 
		storage::size_type const & DestinationLayerOffset,
		storage::size_type const & DestinationFaceOffset,
		storage::size_type const & DestinationlevelOffset);
*/

	std::size_t block_size(format const & Format);
	glm::uvec3 block_dimensions(format const & Format);
	std::size_t component_count(format const & Format);
	bool is_compressed(format const & Format);

	std::size_t level_count(storage::dimensions1_type const & Dimensions);
	std::size_t level_count(storage::dimensions2_type const & Dimensions);
	std::size_t level_count(storage::dimensions3_type const & Dimensions);

}//namespace gli

#include "storage.inl"

#endif//GLI_CORE_STORAGE_INCLUDED
