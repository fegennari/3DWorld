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
/// @file gli/core/addressing.hpp
/// @date 2012-11-19 / 2012-11-19
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#ifndef GLI_CORE_ADDRESSING_INCLUDED
#define GLI_CORE_ADDRESSING_INCLUDED

#include "storage.hpp"

namespace gli{
namespace detail
{
	size_type imageAddressing(
		storage const & Storage,
		size_type const & LayerOffset, 
		size_type const & FaceOffset, 
		size_type const & LevelOffset);

	size_type texelLinearAdressing(
		dimensions1_type const & Dimensions,
		dimensions1_type const & TexelCoord);

	size_type texelLinearAdressing(
		dimensions2_type const & Dimensions,
		dimensions2_type const & TexelCoord);

	size_type texelLinearAdressing(
		dimensions3_type const & Dimensions,
		dimensions3_type const & TexelCoord);

	size_type texelMortonAdressing(
		dimensions1_type const & Dimensions,
		dimensions1_type const & TexelCoord);

	size_type texelMortonAdressing(
		dimensions2_type const & Dimensions,
		dimensions2_type const & TexelCoord);

	size_type texelMortonAdressing(
		dimensions3_type const & Dimensions,
		dimensions3_type const & TexelCoord);

}//namespace detail
/*
	template <typename texture>
	class sampler
	{
		enum wrap
		{
			MIRROR,
			CLAMP_TO_EDGE
		};

		enum filter
		{
			NEAREST,
			LINEAR,
			CUBIC
		};

	public:
		sampler(texture const & Texture)
		{
			switch(Texture.format())
			{
			default:
				break;
			case RGB_DXT1:
				Impl.reset(new impl_fetch_dxt1());
				break;
			}
		}

		template <typename genType>
		genType operator() (texture::dimension_type const & Dimensions) const

	private:
		struct impl
		{
			virtual ~impl(){}

			template <typename genType>
			virtual genType operator() (texture::dimension_type const & Dimensions) const = 0;
		};

		class impl_fetch_color : public impl
		{
		public: 
			impl_fetch_color(){}

			template <typename genType>
			virtual genType operator() (
				texture const & Texture,
				texture::texcoord_type const & Texcoord) const;

		private:

		};

		class impl_fetch_dxt1 : public impl
		{
		public:
			impl_fetch_dxt1(){}

			template <typename genType>
			virtual genType operator() (texture::dimension_type const & Dimensions) const;

		private:

		};

		shared_ptr<impl> Impl;
	};
*/
}//namespace gli

#include "addressing.inl"

#endif//GLI_CORE_ADDRESSING_INCLUDED
