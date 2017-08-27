// 3D World
// by Frank Gennari
// 3/19/05
#ifndef _ALLOCATORS_H_
#define _ALLOCATORS_H_

template <typename T, std::size_t pool_size = 1024>
class pool_allocator {
private:
	std::vector<T*> d_pools;
	T*              d_next;
	T*              d_end;
public:
	template <typename O> struct rebind {typedef pool_allocator<O, pool_size> other;};
	pool_allocator() : d_next(), d_end() {}
	~pool_allocator() {
		std::for_each(this->d_pools.rbegin(), this->d_pools.rend(), [](T* memory){ operator delete(memory); });
	}
	typedef T value_type;
	T* allocate(std::size_t n) {
		if (std::size_t(this->d_end - this->d_next) < n) {
			if (pool_size < n) {
				// custom allocation for bigger number of objects
				this->d_pools.push_back(static_cast<T*>(operator new(sizeof(T) * n)));
				return this->d_pools.back();
			}
			this->d_pools.push_back(static_cast<T*>(operator new(sizeof(T) * pool_size)));
			this->d_next = this->d_pools.back();
			this->d_end  = this->d_next + pool_size;
		}
		T* rc(this->d_next);
		this->d_next += n;
		return rc;
	}
	void deallocate(T*, std::size_t) {
		// this could try to recycle buffers
	}
};

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

