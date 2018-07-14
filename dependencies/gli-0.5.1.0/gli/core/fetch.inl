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
/// @file gli/gtx/fetch.inl
/// @date 2008-12-19 / 2013-01-13
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli
{
	template <typename genType>
	inline genType texelFetch
	(
		texture2D const & Texture, 
		texture2D::dimensions_type const & TexCoord,
		texture2D::size_type const & Level
	)
	{
		assert(!Texture.empty());
		assert(!is_compressed(Texture.format()));

		image::dimensions_type Dimensions = Texture[Level].dimensions();
		genType const * const Data = reinterpret_cast<genType const * const >(Texture[Level].data());

		return reinterpret_cast<genType const * const>(Data)[TexCoord.x + TexCoord.y * Dimensions.x];
	}

	template <typename genType>
	void texelWrite
	(
		texture2D & Texture,
		texture2D::dimensions_type const & Texcoord,
		texture2D::size_type const & Level,
		genType const & Color
	)
	{
		assert(!Texture.empty());
		assert(!is_compressed(Texture.format()));

		genType * Data = Texture[Level].data<genType>();
		std::size_t Index = Texcoord.x + Texcoord.y * Texture[Level].dimensions().x;
		
		std::size_t Capacity = Texture[Level].size();
		assert(Index < Capacity);

		*(Data + Index) = Color;
	}

	template <typename genType>
	inline genType textureLod
	(
		texture2D const & Texture, 
		texture2D::texcoord_type const & TexCoord, 
		texture2D::size_type const & Level
	)
	{
		//assert(Texture.format() == R8U || Texture.format() == RG8U || Texture.format() == RGB8U || Texture.format() == RGBA8U);

		image::dimensions_type Dimensions = Texture[Level].dimensions(); 
		genType const * const Data = reinterpret_cast<genType const * const>(Texture[Level].data());

		std::size_t s_below = std::size_t(glm::floor(TexCoord.s * float(Dimensions.x - 1)));
		std::size_t s_above = std::size_t(glm::ceil( TexCoord.s * float(Dimensions.x - 1)));
		std::size_t t_below = std::size_t(glm::floor(TexCoord.t * float(Dimensions.y - 1)));
		std::size_t t_above = std::size_t(glm::ceil( TexCoord.t * float(Dimensions.y - 1)));

		float s_step = 1.0f / float(Dimensions.x);
		float t_step = 1.0f / float(Dimensions.y);

		float s_below_normalized = s_below / float(Dimensions.x);
		float s_above_normalized = s_above / float(Dimensions.x);
		float t_below_normalized = t_below / float(Dimensions.y);
		float t_above_normalized = t_above / float(Dimensions.y);

		genType Value1 = reinterpret_cast<genType const * const>(Data)[s_below + t_below * Dimensions.x];
		genType Value2 = reinterpret_cast<genType const * const>(Data)[s_above + t_below * Dimensions.x];
		genType Value3 = reinterpret_cast<genType const * const>(Data)[s_above + t_above * Dimensions.x];
		genType Value4 = reinterpret_cast<genType const * const>(Data)[s_below + t_above * Dimensions.x];

		float BlendA = float(TexCoord.s - s_below_normalized) * float(Dimensions.x - 1);
		float BlendB = float(TexCoord.s - s_below_normalized) * float(Dimensions.x - 1);
		float BlendC = float(TexCoord.t - t_below_normalized) * float(Dimensions.y - 1);

		genType ValueA(glm::mix(Value1, Value2, BlendA));
		genType ValueB(glm::mix(Value4, Value3, BlendB));

		return genType(glm::mix(ValueA, ValueB, BlendC));
	}

}//namespace gli
