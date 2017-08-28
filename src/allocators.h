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
	template <typename O> struct rebind {typedef single_free_list_allocator<O> other;};
	~single_free_list_allocator() {
		for (auto i = free_list.begin(); i != free_list.end(); ++i) {delete *i;}
	}
	T* allocate(std::size_t n) {
		assert(n == 1);
		if (free_list.empty()) {return new T;}
		T *ptr(free_list.back());
		free_list.pop_back();
		return ptr;
	}
	void deallocate(T* ptr, std::size_t n) {
		assert(n == 1);
		free_list.push_back(ptr);
	}
};

#endif // _ALLOCATORS_H_

