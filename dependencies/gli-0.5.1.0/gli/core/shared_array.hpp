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
/// @file gli/core/shared_array.hpp
/// @date 2008-12-19 / 2013-01-12
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#ifndef GLI_SHARED_ARRAY_INCLUDED
#define GLI_SHARED_ARRAY_INCLUDED

namespace gli
{
	template <typename T>
	class shared_array
	{
	public:
		shared_array();
		shared_array(shared_array const & SharedArray);
		shared_array(T * Pointer);
		virtual ~shared_array();

		void reset();
		void reset(T * Pointer);

		T & operator*();
		T * operator->();
		T const & operator*() const;
		T const * const operator->() const;

		T * get();
		T const * const get() const;

		shared_array & operator=(shared_array const & SharedArray);
		bool operator==(shared_array const & SharedArray) const;
		bool operator!=(shared_array const & SharedArray) const;

	private:
		long * Counter;
		T * Pointer;
	};
}//namespace gli

#include "shared_array.inl"

#endif //GLI_SHARED_ARRAY_INCLUDED
