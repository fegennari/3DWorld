// 3D World - OpenGL Virtual Environment
// by Frank Gennari
// 6/22/05

#ifndef _MEEMORY_ALLOC_H_
#define _MEEMORY_ALLOC_H_

#include <vector>
#include <assert.h>


template <typename T> class reusable_mem { // can be allocated and freed but only hits the system memory allocation when the requested size increases

	bool alloced;
	vector<T> data;

public:
	reusable_mem() : alloced(0) {}

	void reusable_malloc(T *&alloced_data, size_t alloc_size) {
		assert(alloc_size > 0 && !alloced);
		alloced = 1;
		data.resize(alloc_size);
		alloced_data = &data.front();
	}

	void reusable_free(T *&alloced_data) {
		assert(alloced && alloced_data == &data.front());
		alloced_data = NULL;
		alloced      = 0;
	}
};



#endif

