// 3D World
// by Frank Gennari
// 3/19/05
#ifndef _ALLOCATORS_H_
#define _ALLOCATORS_H_

template <typename T>
class single_free_list_allocator {
	vector<T *> free_list;
public:
	typedef T value_type;
	typedef T* pointer;
	typedef const T* const_pointer;

	template <typename O> struct rebind {typedef single_free_list_allocator<O> other;};
	void operator=(single_free_list_allocator const &a) {assert(free_list.empty() && a.free_list.empty());} // free list must be empty during copy
	single_free_list_allocator(single_free_list_allocator const &a) {assert(a.free_list.empty());}
	single_free_list_allocator() {}
	~single_free_list_allocator() {
		for (auto i = free_list.begin(); i != free_list.end(); ++i) {free(*i);}
	}
	T* allocate(std::size_t n) {
		assert(n == 1);
		if (free_list.empty()) {return static_cast<pointer>(malloc(sizeof(T)));}
		T *ptr(free_list.back());
		free_list.pop_back();
		return ptr;
	}
	void deallocate(T* ptr, std::size_t n) {
		assert(n == 1);
		free_list.push_back(ptr);
	}
	void construct(pointer p, const T& val) {new(static_cast<void*>(p)) T(val);}
	void construct(pointer p) {new(static_cast<void*>(p)) T();}
	void destroy(pointer p) {p->~T();}
};

#endif // _ALLOCATORS_H_

