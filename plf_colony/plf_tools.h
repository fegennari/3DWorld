// Copyright (c) 2026, Matthew Bentley (mattreecebentley@gmail.com) www.plflib.org

// zLib license (https://www.zlib.net/zlib_license.html):
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not
// 	claim that you wrote the original software. If you use this software
// 	in a product, an acknowledgement in the product documentation would be
// 	appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
// 	misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.



#if !defined(PLF_COMPILER_DEFINES)
	#define PLF_COMPILER_DEFINES // Persistent define, so one container hold another plf container as elements doesn't cause re-definitions
	// Compiler-specific defines:

	#include <memory> // std::allocator/allocator_traits. This inclusion is also to pull in at least one standard library header so that any necessary library macros are usable below (eg. _LIBCPP_VERSION)

	// Define default cases before possibly redefining:
	#define PLF_NOEXCEPT throw()
	#define PLF_NOEXCEPT_ALLOCATOR
	#define PLF_CONSTEXPR
	#define PLF_CONSTFUNC


	#define PLF_EXCEPTIONS_SUPPORT

	#if !(defined(__cpp_exceptions) || defined(__EXCEPTIONS) || defined(_CPPUNWIND))
		#undef PLF_EXCEPTIONS_SUPPORT
		#include <exception> // std::terminate
	#endif


	#if defined(_MSC_VER) && !defined(__clang__) && !defined(__GNUC__)
		#if _MSC_VER >= 1600
			#define PLF_MOVE_SEMANTICS_SUPPORT
		#endif
		#if _MSC_VER >= 1700
			#define PLF_TYPE_TRAITS_SUPPORT
			#define PLF_ALLOCATOR_TRAITS_SUPPORT
		#endif
		#if _MSC_VER >= 1800
			#define PLF_VARIADICS_SUPPORT // Variadics, in this context, means both variadic templates and variadic macros are supported
			#define PLF_DEFAULT_SUPPORT // Both support for default template arguments and defaulted class functions
			#define PLF_INITIALIZER_LIST_SUPPORT
		#endif
		#if _MSC_VER >= 1900
			#define PLF_CPP11_SUPPORT
			#define PLF_ALIGNMENT_SUPPORT
			#undef PLF_NOEXCEPT
			#undef PLF_NOEXCEPT_ALLOCATOR
			#define PLF_NOEXCEPT noexcept
			#define PLF_NOEXCEPT_ALLOCATOR noexcept(noexcept(allocator_type()))
			#define PLF_IS_ALWAYS_EQUAL_SUPPORT
			#define PLF_VOIDT_SUPPORT
		#endif

		#if defined(_MSVC_LANG)
			#if _MSVC_LANG >= 201703L
				#undef PLF_CONSTEXPR
				#define PLF_CONSTEXPR constexpr
			#endif

			#if _MSVC_LANG >= 202002L && _MSC_VER >= 1929
				#define PLF_CPP20_SUPPORT
				#undef PLF_CONSTFUNC
				#define PLF_CONSTFUNC constexpr

				#if _MSC_VER < 1951
					// Suppress incorrect (MSVC older version bug) warnings re: constant expressions in constexpr-if statements
					#pragma warning ( push )
					#pragma warning ( disable : 4127 )
				#endif
			#endif

			#if _MSVC_LANG >= 202302L && _MSC_VER >= 1944
				#define PLF_CONSTEVAL_SUPPORT
			#endif
		#endif
	#elif defined(__cplusplus) && __cplusplus >= 201103L // C++11 support, at least
		#if defined(__GNUC__) && defined(__GNUC_MINOR__) && !defined(__clang__) // If compiler is GCC/G++
			#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 3) || __GNUC__ > 4
				#define PLF_MOVE_SEMANTICS_SUPPORT
				#define PLF_VARIADICS_SUPPORT
			#endif
			#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 4) || __GNUC__ > 4
				#define PLF_DEFAULT_SUPPORT
				#define PLF_INITIALIZER_LIST_SUPPORT
			#endif
			#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 6) || __GNUC__ > 4
				#undef PLF_NOEXCEPT
				#undef PLF_NOEXCEPT_ALLOCATOR
				#define PLF_NOEXCEPT noexcept
				#define PLF_NOEXCEPT_ALLOCATOR noexcept(noexcept(allocator_type()))
			#endif
			#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 7) || __GNUC__ > 4
				#define PLF_ALLOCATOR_TRAITS_SUPPORT
			#endif
			#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 8) || __GNUC__ > 4
				#define PLF_ALIGNMENT_SUPPORT
			#endif
			#if __GNUC__ >= 5 // GCC v4.9 and below do not support std::is_trivially_copyable
				#define PLF_TYPE_TRAITS_SUPPORT
			#endif
			#if __GNUC__ > 6
				#define PLF_CPP11_SUPPORT
				#define PLF_IS_ALWAYS_EQUAL_SUPPORT
			#endif
		#elif defined(__clang__) && !defined(__GLIBCXX__) && !defined(_LIBCPP_CXX03_LANG) && __clang_major__ >= 3
			#define PLF_DEFAULT_SUPPORT
			#define PLF_ALLOCATOR_TRAITS_SUPPORT
			#define PLF_TYPE_TRAITS_SUPPORT

			#if __has_feature(cxx_alignas) && __has_feature(cxx_alignof)
				#define PLF_ALIGNMENT_SUPPORT
			#endif
			#if __has_feature(cxx_noexcept)
				#define PLF_CPP11_SUPPORT
				#undef PLF_NOEXCEPT
				#undef PLF_NOEXCEPT_ALLOCATOR
				#define PLF_NOEXCEPT noexcept
				#define PLF_NOEXCEPT_ALLOCATOR noexcept(noexcept(allocator_type()))
			#endif
			#if __clang_major__ >= 4 || __clang_minor__ >= 8
				#define PLF_IS_ALWAYS_EQUAL_SUPPORT
			#endif
			#if __has_feature(cxx_rvalue_references) && !defined(_LIBCPP_HAS_NO_RVALUE_REFERENCES)
				#define PLF_MOVE_SEMANTICS_SUPPORT
			#endif
			#if __has_feature(cxx_variadic_templates) && !defined(_LIBCPP_HAS_NO_VARIADICS)
				#define PLF_VARIADICS_SUPPORT
			#endif
			#if (__clang_major__ == 3 && __clang_minor__ >= 1) || __clang_major__ > 3
				#define PLF_INITIALIZER_LIST_SUPPORT
			#endif
		#elif defined(__GLIBCXX__) // Using another compiler type with libstdc++ - we are assuming full c++11 compliance for compiler - which may not be true
			#define PLF_DEFAULT_SUPPORT

			#if __GLIBCXX__ >= 20080606
				#define PLF_MOVE_SEMANTICS_SUPPORT
				#define PLF_VARIADICS_SUPPORT
			#endif
			#if __GLIBCXX__ >= 20090421
				#define PLF_INITIALIZER_LIST_SUPPORT
			#endif
			#if __GLIBCXX__ >= 20120322
				#define PLF_ALLOCATOR_TRAITS_SUPPORT
				#undef PLF_NOEXCEPT
				#undef PLF_NOEXCEPT_ALLOCATOR
				#define PLF_NOEXCEPT noexcept
				#define PLF_NOEXCEPT_ALLOCATOR noexcept(noexcept(allocator_type()))
			#endif
			#if __GLIBCXX__ >= 20130322
				#define PLF_ALIGNMENT_SUPPORT
			#endif
			#if __GLIBCXX__ >= 20150422 // libstdc++ v4.9 and below do not support std::is_trivially_copyable
				#define PLF_TYPE_TRAITS_SUPPORT
			#endif
			#if __GLIBCXX__ >= 20160111
				#define PLF_IS_ALWAYS_EQUAL_SUPPORT
			#endif
		#elif defined(_LIBCPP_CXX03_LANG) || defined(_LIBCPP_HAS_NO_RVALUE_REFERENCES) // Special case for checking C++11 support with libCPP
			#if !defined(_LIBCPP_HAS_NO_VARIADICS)
				#define PLF_VARIADICS_SUPPORT
			#endif
		#else // Assume type traits and initializer support for other compilers and standard library implementations
			#define PLF_DEFAULT_SUPPORT
			#define PLF_MOVE_SEMANTICS_SUPPORT
			#define PLF_VARIADICS_SUPPORT
			#define PLF_TYPE_TRAITS_SUPPORT
			#define PLF_ALLOCATOR_TRAITS_SUPPORT
			#define PLF_ALIGNMENT_SUPPORT
			#define PLF_INITIALIZER_LIST_SUPPORT
			#undef PLF_NOEXCEPT
			#undef PLF_NOEXCEPT_ALLOCATOR
			#define PLF_NOEXCEPT noexcept
			#define PLF_NOEXCEPT_ALLOCATOR noexcept(noexcept(allocator_type()))
			#define PLF_IS_ALWAYS_EQUAL_SUPPORT
			#define PLF_CPP11_SUPPORT
		#endif

		#if __cplusplus >= 201703L && ((defined(__clang__) && ((__clang_major__ == 3 && __clang_minor__ == 9) || __clang_major__ > 3)) || (defined(__GNUC__) && __GNUC__ >= 7) || (!defined(__clang__) && !defined(__GNUC__))) // assume correct C++17 implementation for non-gcc/clang compilers
			#undef PLF_CONSTEXPR
			#define PLF_CONSTEXPR constexpr
			#define PLF_VOIDT_SUPPORT
		#endif

		#if __cplusplus >= 202001L && ((((defined(__clang__) && __clang_major__ >= 15) || (defined(__GNUC__) && (__GNUC__ >= 12))) && ((defined(_LIBCPP_VERSION) && _LIBCPP_VERSION >= 15) || (defined(__GLIBCXX__) &&	_GLIBCXX_RELEASE >= 12))) || (!defined(__clang__) && !defined(__GNUC__))) // detecting the library version is required so that, for example, std::to_address can be used when PLF_CPP20_SUPPORT is defined
			#define PLF_CPP20_SUPPORT
			#undef PLF_CONSTFUNC
			#define PLF_CONSTFUNC constexpr
		#endif

		#if __cplusplus >= 202302L && (((defined(__clang__) && __clang_major__ >= 14) || (defined(__GNUC__) && (__GNUC__ >= 12))) || (!defined(__clang__) && !defined(__GNUC__)))
			#define PLF_CONSTEVAL_SUPPORT
		#endif
	#endif


	#if defined(PLF_IS_ALWAYS_EQUAL_SUPPORT) && defined(PLF_MOVE_SEMANTICS_SUPPORT) && defined(PLF_ALLOCATOR_TRAITS_SUPPORT) && (__cplusplus >= 201703L || (defined(_MSVC_LANG) && (_MSVC_LANG >= 201703L)))
		#define PLF_NOEXCEPT_MOVE_ASSIGN(the_allocator) noexcept(std::allocator_traits<the_allocator>::propagate_on_container_move_assignment::value || std::allocator_traits<the_allocator>::is_always_equal::value)
		#define PLF_NOEXCEPT_SWAP(the_allocator) noexcept(std::allocator_traits<the_allocator>::propagate_on_container_swap::value || std::allocator_traits<the_allocator>::is_always_equal::value)
	#else
		#define PLF_NOEXCEPT_MOVE_ASSIGN(the_allocator)
		#define PLF_NOEXCEPT_SWAP(the_allocator)
	#endif


	#ifdef PLF_ALLOCATOR_TRAITS_SUPPORT
		#ifdef PLF_VARIADICS_SUPPORT
			#define PLF_CONSTRUCT(the_allocator, allocator_instance, location, ...) std::allocator_traits<the_allocator>::construct(allocator_instance, location, __VA_ARGS__)
		#else
			#define PLF_CONSTRUCT(the_allocator, allocator_instance, location, data)	std::allocator_traits<the_allocator>::construct(allocator_instance, location, data)
		#endif

		#define PLF_DESTROY(the_allocator, allocator_instance, location)				std::allocator_traits<the_allocator>::destroy(allocator_instance, location)
		#define PLF_ALLOCATE(the_allocator, allocator_instance, size, hint)			std::allocator_traits<the_allocator>::allocate(allocator_instance, size, hint)
		#define PLF_DEALLOCATE(the_allocator, allocator_instance, location, size)	std::allocator_traits<the_allocator>::deallocate(allocator_instance, location, size)
	#else
		#ifdef PLF_VARIADICS_SUPPORT
			#define PLF_CONSTRUCT(the_allocator, allocator_instance, location, ...) 	(allocator_instance).construct(location, __VA_ARGS__)
		#else
			#define PLF_CONSTRUCT(the_allocator, allocator_instance, location, data)	(allocator_instance).construct(location, data)
		#endif

		#define PLF_DESTROY(the_allocator, allocator_instance, location)				(allocator_instance).destroy(location)
		#define PLF_ALLOCATE(the_allocator, allocator_instance, size, hint)			(allocator_instance).allocate(size, hint)
		#define PLF_DEALLOCATE(the_allocator, allocator_instance, location, size)	(allocator_instance).deallocate(location, size)
	#endif


	#define PLF_CONSTRUCT_ELEMENT(location, element) PLF_CONSTRUCT(allocator_type, *this, pointer_cast<pointer>(location), element)

#endif // PLF_COMPILER_DEFINES




#if (defined(PLF_INCLUDE_TOOLS) || defined(PLF_INCLUDE_UNINITIALIZED_TOOLS)) && !defined(PLF_TOOLS)  // uninitialized_tools uses plf::make_move_iterator, hence needs these tools
	#define PLF_TOOLS

	#include <memory> // to_address

	#ifdef PLF_CPP20_SUPPORT
		#include <ranges>
	#endif

	#if defined(PLF_MOVE_SEMANTICS_SUPPORT)
		#include <iterator> // move_iterator
	#endif


	namespace plf
	{

		// std:: tool replacements for C++03/98/11 support:
		template <bool condition, class T = void>
		struct enable_if
		{
			typedef T type;
		};

		template <class T>
		struct enable_if<false, T>
		{};



		template <bool flag, class is_true, class is_false> struct conditional;

		template <class is_true, class is_false> struct conditional<true, is_true, is_false>
		{
			typedef is_true type;
		};

		template <class is_true, class is_false> struct conditional<false, is_true, is_false>
		{
			typedef is_false type;
		};



		template <class element_type>
		struct less
		{
			bool operator() (const element_type &a, const element_type &b) const PLF_NOEXCEPT
			{
				return a < b;
			}
		};



		template<class element_type>
		struct equal_to
		{
			const element_type &value;

			explicit equal_to(const element_type &store_value) PLF_NOEXCEPT:
				value(store_value)
			{}

			bool operator() (const element_type &compare_value) const PLF_NOEXCEPT
			{
				return value == compare_value;
			}
		};



		#ifdef PLF_CPP20_SUPPORT
			#define PLF_TO_ADDRESS(pointer) std::to_address(pointer)
		#else
			#define PLF_TO_ADDRESS(pointer) &*(pointer)
		#endif



		// For converting the underlying skipfield storage * type to void * when the allocator supplies non-trivial pointers.
		// The void * conversion is technically unnecessary since it will be implicitly converted when required, but it's more straightforward than having to identify the underlying storage type:
		template <class source_pointer_type>
		static PLF_CONSTFUNC void * void_cast(const source_pointer_type source_pointer) PLF_NOEXCEPT
		{
			return static_cast<void *>(PLF_TO_ADDRESS(source_pointer));
		}



		#ifdef PLF_MOVE_SEMANTICS_SUPPORT
			template <class iterator_type>
			static PLF_CONSTFUNC std::move_iterator<iterator_type> make_move_iterator(iterator_type it)
			{
				return std::move_iterator<iterator_type>(std::move(it));
			}
		#endif



		enum priority { performance = 1, memory_use = 4};



		#ifdef PLF_CPP20_SUPPORT
			// For matching ranges which return input_iterator's and which match the container's element type:
			template <typename range_type, class element_type>
			concept compatible_range = std::ranges::input_range<range_type> && std::convertible_to<std::ranges::range_reference_t<range_type>, element_type>;

			// For non-C++23-compliant libraries which do not include std::ranges::from_range_t - this is so the rangesv3 constructor overloads will work unambiguously:
			namespace ranges
			{
				struct from_range_t {};
				constexpr from_range_t from_range;
			}
		#endif



		// To simplify conversion when allocator supplies non-raw pointers:
		template <class destination_pointer_type, class source_pointer_type>
		static PLF_CONSTFUNC destination_pointer_type pointer_cast(const source_pointer_type source_pointer) PLF_NOEXCEPT
		{
			#if defined(PLF_TYPE_TRAITS_SUPPORT) && defined(PLF_CPP20_SUPPORT) // constexpr necessary to avoid a branch for every call
				if constexpr (std::is_trivially_constructible<destination_pointer_type>::value)
				{
					if constexpr (std::is_trivially_constructible<source_pointer_type>::value)
					{
						return reinterpret_cast<destination_pointer_type>(source_pointer);
					}
					else
					{
						return reinterpret_cast<destination_pointer_type>(std::to_address(source_pointer));
					}
				}
				else
				{
					return destination_pointer_type(std::to_address(source_pointer));
				}
			#else
				return destination_pointer_type(&*source_pointer);
			#endif
		}

	} // plf namepsace

#endif // PLF_TOOLS




#if defined(PLF_INCLUDE_BIT_TOOLS) && !defined(PLF_BIT_TOOLS)
	#define PLF_BIT_TOOLS // Persistent define so same functions don't get included twice

	#include <cstddef> // std::size_t

	#ifdef PLF_CPP20_SUPPORT
		#include <bit> // std::popcount, std::countr/etc
	#endif

	#if defined(_MSC_VER)
		#include <intrin.h> // _BitScanForward
	#endif


	namespace plf
	{
		// A popcount which works pre-c++20:
		template <typename storage_type>
		static PLF_CONSTFUNC std::size_t popcount(storage_type value)
		{
			#ifdef PLF_CPP20_SUPPORT
				return std::popcount(value);
			#else
				std::size_t total = 0;
				for (; value; ++total) value &= value - 1; // Kernighan's algorithm
				return total;
			#endif
		}



		// These countr/countl implementations work in pre-C++20 modes, but also skip zero-checks in >= C++20 if the architecture was going to end up using BSR instead of some other instruction set.
		// Hence if you use them, you must make sure value != 0 (or for countr_one/countl_one, that value != std::numeric_limits<storage_type>::max()). These implementations do not work with types > unsigned long long.

		template<typename storage_type>
		static PLF_CONSTFUNC std::size_t countr_zero(const storage_type value)
		{
			#ifdef PLF_CPP20_SUPPORT
				[[assume(value != 0)]];
			#endif

			#ifdef PLF_CONSTEVAL_SUPPORT
				if consteval { return std::countr_zero(value); } // Use whatever the library's constexpr-friendly version is
			#endif

			#if defined(_MSC_VER) && (!defined(PLF_CPP20_SUPPORT) || (!defined(__AVX2__) && !defined(_M_CEE_PURE) && ((defined(_M_IX86) && !defined(_M_HYBRID_X86_ARM64)) || (defined(_M_X64) && !defined(_M_ARM64EC))))) // Matches MSVC and clang under MSVC
				unsigned long result;

				#if defined(_M_X64) && _MSC_VER >= 1900 // support for _BitScanForward64
					if PLF_CONSTEXPR (sizeof(storage_type) > 4)
					{
						_BitScanForward64(&result, value);
					}
					else
				#endif
				{ // Note: the ~ appears to invert the ::max() value here (ie. interpreting -1 as an unsigned type), but actually, since the type may be a smaller number of bits than unsigned int, this is potentially only the inversion of the lower bits.
					_BitScanForward(&result, static_cast<unsigned int>(~static_cast<storage_type>(-1) | value));
				}

				return static_cast<std::size_t>(result);
			#elif ((defined(__GNUC__) && __GNUC__ >= 4) || defined(__clang__)) && !defined(PLF_CPP20_SUPPORT)
				if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(unsigned))
				{
					return static_cast<std::size_t>(__builtin_ctz(value));
				}
				else if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(unsigned long))
				{
					return static_cast<std::size_t>(__builtin_ctzl(value));
				}
				#ifdef PLF_CPP11_SUPPORT
					else if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(unsigned long long))
				#else
					else if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(std::size_t) && sizeof(std::size_t) >= 8)
				#endif
				{
					return static_cast<std::size_t>(__builtin_ctzll(value));
				}
			#elif defined(PLF_CPP20_SUPPORT)
				return std::countr_zero(value);
			#else
				for (storage_type bit_index = 0; ; ++bit_index)
				{
					if (value & (storage_type(1) << bit_index)) return static_cast<std::size_t>(bit_index);
				}
			#endif
		}



		template<typename storage_type>
		static PLF_CONSTFUNC std::size_t countr_one(const storage_type value)
		{
			return plf::countr_zero(~value);
		}




		template<typename storage_type>
		static PLF_CONSTFUNC std::size_t countl_zero(const storage_type value)
		{
			#ifdef PLF_CPP20_SUPPORT
				[[assume(value != 0)]];
			#endif

			#ifdef PLF_CONSTEVAL_SUPPORT
				if consteval { return std::countl_zero(value); }
			#endif

			#if defined(_MSC_VER) && (!defined(PLF_CPP20_SUPPORT) || (!defined(__AVX2__) && !defined(_M_CEE_PURE) && ((defined(_M_IX86) && !defined(_M_HYBRID_X86_ARM64)) || (defined(_M_X64) && !defined(_M_ARM64EC)))))
				unsigned long result;

				#if defined(_M_X64) && _MSC_VER >= 1900
					if PLF_CONSTEXPR (sizeof(storage_type) > 4)
					{
						_BitScanReverse64(&result, value);
					}
					else
				#endif
				{
					_BitScanReverse(&result, value);
				}

				return static_cast<std::size_t>(((sizeof(storage_type) * 8) - 1) - result);
			#elif ((defined(__GNUC__) && __GNUC__ >= 4) || defined(__clang__)) && !defined(PLF_CPP20_SUPPORT)
				if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(unsigned))
				{
					return static_cast<std::size_t>(__builtin_clz(value) - (sizeof(unsigned) - sizeof(storage_type)));
				}
				else if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(unsigned long))
				{
					return static_cast<std::size_t>(__builtin_clzl(value) - (sizeof(unsigned long) - sizeof(storage_type)));
				}
				#ifdef PLF_CPP11_SUPPORT
					else if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(unsigned long long))
					{
						return static_cast<std::size_t>(__builtin_clzll(value) - (sizeof(unsigned long long) - sizeof(storage_type)));
					}
				#else
					else if PLF_CONSTEXPR (sizeof(storage_type) <= sizeof(std::size_t) && sizeof(std::size_t) >= 8)
					{
						return static_cast<std::size_t>(__builtin_clzll(value) - (sizeof(std::size_t) - sizeof(storage_type)));
					}
				#endif
			#elif defined(PLF_CPP20_SUPPORT)
				return std::countl_zero(value);
			#else
				for (storage_type bit_index = (sizeof(storage_type) * 8) - 1; ; --bit_index)
				{
					if (value & (storage_type(1) << bit_index)) return ((sizeof(storage_type) * 8) - 1) - bit_index;
				}
			#endif
		}



		template<typename storage_type>
		static PLF_CONSTFUNC std::size_t countl_one(const storage_type value)
		{
			return plf::countl_zero(~value);
		}

	} // plf namespace

#endif // PLF_INCLUDE_BIT_TOOLS




#if defined(PLF_INCLUDE_UNINITIALIZED_TOOLS) && !defined(PLF_UNINITIALIZED_TOOLS)
	#define PLF_UNINITIALIZED_TOOLS

	#if defined(PLF_TYPE_TRAITS_SUPPORT) && defined(PLF_VOIDT_SUPPORT)
		#include <type_traits> // void_t, false_type
		#include <utility> // declval
	#endif

	// Note: plf::uninitialized_move relies on PLF_TOOLS for plf::make_move_iterator, hence <memory> is already included at this point

	namespace plf
	{

		// Allocator-aware uninitialized_copy/move/fill_n:

		// Template to check whether an allocator has a custom construct function, or just relies on allocator_traits (eg. std::allocator since C++20):
		#ifdef PLF_VOIDT_SUPPORT
			template<typename allocator_type, typename = std::void_t<>>
			struct allocator_has_construct : std::false_type {};

			// Tests for dummy type int:
			template<typename allocator_type>
			struct allocator_has_construct< allocator_type, std::void_t<decltype(std::declval<allocator_type&>().construct(std::declval<int*>(), std::declval<int>()))> > : std::true_type {};
		#endif



		template <class allocator_type, class iterator_type, class iterator_type2>
		void uninitialized_copy(iterator_type begin, const iterator_type end, iterator_type2 destination,
		#ifdef PLF_CPP20_SUPPORT
			[[maybe_unused]]
		#endif
			allocator_type &alloc)
		{
			#ifdef PLF_VOIDT_SUPPORT
				if PLF_CONSTEXPR (!allocator_has_construct<allocator_type>::value) // If allocator has no construct method, we can take advantage of optimized routines for POD types
				{
					std::uninitialized_copy(begin, end, destination);
				}
				else
			#endif
			{
				for (; begin != end; ++begin, ++destination)
				{
					PLF_CONSTRUCT(allocator_type, alloc, PLF_TO_ADDRESS(destination), *begin);
				}
			}
		}



		#ifdef PLF_MOVE_SEMANTICS_SUPPORT
			template <class allocator_type, class iterator_type, class iterator_type2>
			void uninitialized_move(iterator_type begin, const iterator_type end, iterator_type2 destination, allocator_type &alloc)
			{
				uninitialized_copy(plf::make_move_iterator(begin), plf::make_move_iterator(end), destination, alloc);
			}
		#endif



		template <class allocator_type, class iterator_type, class element_type>
		void uninitialized_fill_n(iterator_type begin, std::size_t size, const element_type &element,
		#ifdef PLF_CPP20_SUPPORT
			[[maybe_unused]]
		#endif
			allocator_type &alloc)
		{
			#ifdef PLF_VOIDT_SUPPORT
				if PLF_CONSTEXPR (!allocator_has_construct<allocator_type>::value)
				{
					std::uninitialized_fill_n(begin, size, element);
				}
				else
			#endif
			{
				for (; size != 0; ++begin, --size)
				{
					PLF_CONSTRUCT(allocator_type, alloc, PLF_TO_ADDRESS(begin), element);
				}
			}
		}

	} // plf namespace
#endif // PLF_UNINITIALIZED_TOOLS



#undef PLF_INCLUDE_BIT_TOOLS
#undef PLF_INCLUDE_UNINITIALIZED_TOOLS
#undef PLF_INCLUDE_TOOLS
