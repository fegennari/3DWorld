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
/// @file gli/core/shared_array.inl
/// @date 2008-12-19 / 2013-01-12
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

namespace gli
{
	template <typename T>
	shared_array<T>::shared_array() :
		Counter(0),
		Pointer(0)
	{}

	template <typename T>
	shared_array<T>::shared_array
	(
		shared_array<T> const & SharedArray
	) :
		Counter(SharedArray.Counter),
		Pointer(SharedArray.Pointer)
	{
		if(this->Counter)
			(*this->Counter)++;
	}

	template <typename T>
	shared_array<T>::shared_array
	(
		T * Pointer
	) :
		Pointer(0),
		Counter(0)
	{
		if(Pointer)
			this->reset(Pointer);
	}

	template <typename T>
	shared_array<T>::~shared_array()
	{
		this->reset();
	}

	template <typename T>
	void shared_array<T>::reset()
	{
		if(!this->Pointer)
			return;

		(*this->Counter)--;
		if(*this->Counter <= 0)
		{
			delete this->Counter;
			delete[] this->Pointer;
		}

		this->Counter = 0;
		this->Pointer = 0;
	}

	template <typename T>
	void shared_array<T>::reset(T * Pointer)
	{
		this->reset();
	
		this->Counter = new long(1);
		this->Pointer = Pointer;
	}

	template <typename T>
	shared_array<T>& shared_array<T>::operator=
	(
		shared_array<T> const & SharedArray
	)
	{
		this->reset();

		this->Counter = SharedArray.Counter;
		this->Pointer = SharedArray.Pointer;
		(*this->Counter)++;

		return *this;
	}

	template <typename T>
	bool shared_array<T>::operator==(shared_array<T> const & SharedArray) const
	{
		return this->Pointer == SharedArray.Pointer;
	}

	template <typename T>
	bool shared_array<T>::operator!=(shared_array<T> const & SharedArray) const
	{
		return this->Pointer != SharedArray.Pointer;
	}

	template <typename T>
	T & shared_array<T>::operator*()
	{
		return *this->Pointer;
	}

	template <typename T>
	T * shared_array<T>::operator->()
	{
		return this->Pointer;
	}

	template <typename T>
	T const & shared_array<T>::operator*() const
	{
		return * this->Pointer;
	}

	template <typename T>
	T const * const shared_array<T>::operator->() const
	{
		return this->Pointer;
	}

	template <typename T>
	T * shared_array<T>::get()
	{
		return this->Pointer;
	}

	template <typename T>
	T const * const shared_array<T>::get() const
	{
		return this->Pointer;
	}
}//namespace gli
