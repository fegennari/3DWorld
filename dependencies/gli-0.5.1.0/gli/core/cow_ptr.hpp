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
/// @file gli/core/cow_ptr.hpp
/// @date 2012-09-24 / 2013-01-12
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#ifndef GLI_COW_PTR_INCLUDED
#define GLI_COW_PTR_INCLUDED

#include "shared_ptr.hpp"

namespace gli
{
	template <typename T>
	class cow_ptr
	{
	public:
		cow_ptr(T* Pointer) :
			SharedPtr(Pointer)
		{}

		cow_ptr(shared_ptr<T> const & SharedPtr) :
			SharedPtr(SharedPtr)
		{}

		cow_ptr(cow_ptr<T> const & CowPtr) :
			SharedPtr(CowPtr.Pointer)
		{}

		cow_ptr & operator=(cow_ptr<T> const & CowPtr)
		{
			// no need to check for self-assignment with shared_ptr

			this->SharedPtr = CowPtr.SharedPtr; 
			return *this;
		}

		T const & operator*() const
		{
			return *this->SharedPtr;
		}

		T & operator*()
		{
			this->detach();
			return *this->SharedPtr;
		}

		T const * operator->() const
		{
			return this->SharedPtr.operator->();
		}

		T * operator->()
		{
			this->detach();
			return this->SharedPtr.operator->();
		}

	private:
		void detach()
		{
			T* Tmp = this->SharedPtr.get();
			if(!(Tmp == 0 || this->SharedPtr.unique()))
				this->SharedPtr = cow_ptr<T>(new T(*Tmp));
		}

		shared_ptr<T> SharedPtr;
	};
}//namespace gli

#endif//GLI_COW_PTR_INCLUDED
