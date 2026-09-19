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


#ifndef PLF_COLONY_H
#define PLF_COLONY_H

#ifndef PLF_COMPILER_DEFINES
	#define PLF_COLONY_DEFINES // ie. No encapsulating unit/class has previously defined the compiler feature macros in plf_tools.h below, so allow this header to undefine them at it's end.
#endif

#define PLF_INCLUDE_UNINITIALIZED_TOOLS
#define PLF_INCLUDE_TOOLS
#include "plf_tools.h"


#include <algorithm> // std::fill_n, std::sort, std::swap
#include <cassert>	// assert
#include <cstring>	// memset, memcpy, size_t
#include <limits>  // std::numeric_limits
#include <memory> // std::allocator, std::allocator_traits, std::to_address
#include <iterator> // std::bidirectional_iterator_tag, std::move_iterator, std::distance for range insert
#include <stdexcept> // std::length_error


#ifdef PLF_TYPE_TRAITS_SUPPORT
	#include <cstddef> // offsetof, used in blank()
	#include <type_traits> // std::is_trivially_destructible, type_identity_t, etc
#endif

#ifdef PLF_MOVE_SEMANTICS_SUPPORT
	#include <utility> // std::move
#endif

#ifdef PLF_INITIALIZER_LIST_SUPPORT
	#include <initializer_list>
#endif


#ifdef PLF_CPP20_SUPPORT
	#include <concepts>
	#include <compare> // std::strong_ordering
	#include <ranges>


	namespace plf
	{
		// For getting std:: overload for reverse_iterator to match colony iterators specifically (see bottom of header):
		template <class T>
		concept colony_iterator_concept = requires { typename T::colony_iterator_tag; };
	}
#endif



#ifndef PLF_SORT_FUNCTION
	#define PLF_SORT_FUNCTION std::sort
	#define PLF_SORT_FUNCTION_DEFINED
#endif



namespace plf
{


struct limits
{
	size_t min, max;
	PLF_CONSTFUNC limits(const size_t minimum, const size_t maximum) PLF_NOEXCEPT : min(minimum), max(maximum) {}
};




template <class element_type, class allocator_type = std::allocator<element_type>, plf::priority priority = performance>
class colony : private allocator_type // Empty base class optimisation - inheriting allocator functions
{
	#ifdef PLF_ALIGNMENT_SUPPORT
		typedef typename plf::conditional<(priority == performance && (sizeof(element_type) > 10 || alignof(element_type) > 10)), unsigned short, unsigned char>::type		skipfield_type; // Note: unsigned short is equivalent to uint_least16_t ie. Using 16-bit unsigned integer in best-case scenario, greater-than-16-bit unsigned integer where platform doesn't support 16-bit types. unsigned char is always == 1 byte, as opposed to uint_8, which may not be.
	#else
		typedef typename plf::conditional<(priority == performance && sizeof(element_type) > 10), unsigned short, unsigned char>::type		skipfield_type;
	#endif

public:
	// Standard container typedefs:

	#ifdef PLF_ALLOCATOR_TRAITS_SUPPORT
		typedef typename std::allocator_traits<allocator_type>::size_type 			size_type;
		typedef typename std::allocator_traits<allocator_type>::difference_type 	difference_type;
		typedef typename std::allocator_traits<allocator_type>::pointer				pointer;
		typedef typename std::allocator_traits<allocator_type>::const_pointer		const_pointer;
	#else
		typedef typename allocator_type::size_type			size_type;
		typedef typename allocator_type::difference_type	difference_type;
		typedef typename allocator_type::pointer				pointer;
		typedef typename allocator_type::const_pointer		const_pointer;
	#endif

	typedef element_type value_type;
	typedef element_type &			reference;
	typedef const element_type &	const_reference;


	// Iterator forward declarations:
	template <bool is_const> class			colony_iterator;
	typedef colony_iterator<false>			iterator;
	typedef colony_iterator<true> 			const_iterator;
	friend class colony_iterator<false>;
	friend class colony_iterator<true>;

	template <bool is_const_r> class			colony_reverse_iterator;
	typedef colony_reverse_iterator<false>	reverse_iterator;
	typedef colony_reverse_iterator<true>	const_reverse_iterator;
	friend class colony_reverse_iterator<false>;
	friend class colony_reverse_iterator<true>;


	#ifdef PLF_ALIGNMENT_SUPPORT
		// The element as allocated in memory needs to be at-least 2*skipfield_type width in order to support free list indexes in erased element memory space, so:
		// make the size of this struct the larger of alignof(T), sizeof(T) or 2*skipfield_type (the latter is only relevant for type char/uchar), and
		// make the alignment alignof(T).
		// This type is used mainly for correct pointer arithmetic while iterating over elements in memory.
		struct alignas(alignof(element_type)) aligned_element_struct
		{
			 // Using char as sizeof is always guaranteed to be 1 byte regardless of the number of bits in a byte on given computer, whereas for example, uint8_t would fail on machines where there are more than 8 bits in a byte eg. Texas Instruments C54x DSPs.
			char data[
			(sizeof(element_type) < (sizeof(skipfield_type) * 2)) ?
			((sizeof(skipfield_type) * 2) < alignof(element_type) ? alignof(element_type) : (sizeof(skipfield_type) * 2)) :
			((sizeof(element_type) < alignof(element_type)) ? alignof(element_type) : sizeof(element_type))
			];
		};


		// We combine the allocation of elements and skipfield into one allocation to save performance. This memory must be allocated as an aligned type with the same alignment as T in order for the elements to align with memory boundaries correctly (which won't happen if we allocate as char or uint_8). But the larger the sizeof in the type we use for allocation, the greater the chance of creating a lot of unused memory in the skipfield portion of the allocated block. So we create a type that is sizeof(alignof(T)), as in most cases alignof(T) < sizeof(T). If alignof(t) >= sizeof(t) this makes no difference.
		struct alignas(alignof(element_type)) aligned_allocation_struct
		{
		  char data[alignof(element_type)];
		};
	#else
		struct aligned_element_struct
		{
			char data[(sizeof(element_type) < (sizeof(skipfield_type) * 2)) ? (sizeof(skipfield_type) * 2) : sizeof(element_type)];
		};

		struct aligned_allocation_struct
		{
		  char data[sizeof(skipfield_type)];
		};
	#endif


private:

	// Calculate the capacity of a group's elements+skipfield memory block when expressed in multiples of the value_type's alignment (rounding up).
	static size_type get_aligned_block_capacity(const skipfield_type elements_per_group)
	{
		return ((elements_per_group * (sizeof(aligned_element_struct) + sizeof(skipfield_type))) + sizeof(skipfield_type) + sizeof(aligned_allocation_struct) - 1) / sizeof(aligned_allocation_struct);
	}


	// forward declarations for typedefs below
	struct group;
	struct item_index_tuple; // for use in sort()

	// These two need to be raw pointers as instances of this type have pointer arithmetic done on them:
	typedef aligned_element_struct *	aligned_pointer_type; // pointer to the (potentially overaligned) element type, not the original element type
	typedef skipfield_type *			skipfield_pointer_type;

	#ifdef PLF_ALLOCATOR_TRAITS_SUPPORT
		typedef typename std::allocator_traits<allocator_type>::template rebind_alloc<group>							group_allocator_type;
		typedef typename std::allocator_traits<allocator_type>::template rebind_alloc<skipfield_type>				skipfield_allocator_type;
		typedef typename std::allocator_traits<allocator_type>::template rebind_alloc<aligned_allocation_struct> aligned_struct_allocator_type;
		typedef typename std::allocator_traits<allocator_type>::template rebind_alloc<item_index_tuple> 			tuple_allocator_type;
		typedef typename std::allocator_traits<allocator_type>::template rebind_alloc<unsigned char> 				uchar_allocator_type;

		typedef typename std::allocator_traits<group_allocator_type>::pointer				group_pointer_type;
		typedef typename std::allocator_traits<aligned_struct_allocator_type>::pointer	aligned_struct_pointer_type;
		typedef typename std::allocator_traits<tuple_allocator_type>::pointer				tuple_pointer_type;
		typedef typename std::allocator_traits<uchar_allocator_type>::pointer				uchar_pointer_type;
	#else
		typedef typename allocator_type::template rebind<group>::other 						group_allocator_type;
		typedef typename allocator_type::template rebind<skipfield_type>::other 			skipfield_allocator_type;
		typedef typename allocator_type::template rebind<aligned_allocation_struct>::other	aligned_struct_allocator_type;
		typedef typename allocator_type::template rebind<item_index_tuple>::other			tuple_allocator_type;
		typedef typename allocator_type::template rebind<unsigned char>::other				uchar_allocator_type;

		typedef typename group_allocator_type::pointer				group_pointer_type;
		typedef typename aligned_struct_allocator_type::pointer	aligned_struct_pointer_type;
		typedef typename tuple_allocator_type::pointer				tuple_pointer_type;
		typedef typename uchar_allocator_type::pointer				uchar_pointer_type;
	#endif



	// group == element memory block + skipfield + block metadata

	// Skipfield implementation notes:
	// (1) Follows low-complexity jump-counting pattern rules as described here: archive.org/details/matt_bentley_-_the_low_complexity_jump-counting_pattern
	// (2) Initialized to 0 by-default, which means 'non-erased' ie. either no element has ever been constructed there, or an element has been constructed there. Whereas non-zero means an element has been constructed there and subsequently erased. The value of the first and last non-zero nodes in a run of non-zero nodes determines jump length. See the paper for details.
	// (3) This definition means we can bulk-initialize each group's skipfield to 0, rather than bulk-initialize to non-zero then subsequently change individual skipfield nodes to 0 upon insertion - which is obviously slower. Defining unconstructed elements as 0 has no impact on iteration since they're after end(). Note that this definition is violated slightly during splice: unconstructed element nodes at the end of the destination colony's back block will have their corresponding skipfield nodes flipped to 'erased' in order to make iteration work, because they will not longer be after end() once the source colony's blocks are appended.
	// (4) There will always be one additional skipfield node allocated compared to the group's number of elements. This ensures a faster ++ iterator operation (fewer checks are required when it is present). The extra node is unused and always 0, but checked, and not having it will result in out-of-bounds memory errors.


	struct group
	{
		skipfield_pointer_type					skipfield;			// Skipfield storage. The element and skipfield arrays are allocated contiguously, in a single allocation, in this implementation, hence the skipfield pointer also functions as a 'one-past-end' pointer for the elements array. This is present before elements in the group struct as it is referenced constantly by the ++ operator, hence having it first results in a minor performance increase.
		group_pointer_type						next_group;			// Next group in the linked list of all groups. NULL if no following group. 2nd in struct because it is so frequently used during iteration.
		const aligned_struct_pointer_type 	elements;			// Element storage. Allocated as a block of chars, this memory is then divided between elements & skipfield
		group_pointer_type						previous_group;	// Previous group in the linked list of all groups. NULL if no preceding group.
		skipfield_type 							free_list_head;	// The index of the last erased element in the group. The last erased element will, in turn, contain the number of the index of the next erased element, and so on. If this is == maximum skipfield_type value then free_list is empty ie. no erasures have occurred in the group (or if they have, the erased locations have subsequently been reused via insert/emplace/assign).
		const skipfield_type 					capacity;			// The element capacity of this particular group - can also be calculated from reinterpret_cast<aligned_pointer_type>(group->skipfield) - group->elements, however this space is effectively free due to struct padding and the sizeof(skipfield_type), and calculating it once is faster in benchmarking.
		skipfield_type 							size; 				// The total number of active elements in group - changes with insert and erase commands - used to check for empty group in erase function, as an indication to remove the group. Also used in combination with capacity to check if group is full, which is used in the next/previous/advance/distance overloads, and range-erase.
		group_pointer_type						erasures_list_next_group, erasures_list_previous_group; // The next and previous groups in the list of groups with erasures ie. with active erased-element free lists. NULL if no next or previous group.
		size_type									group_number;		// Used for comparison (> < >= <= <=>) iterator operators (used by distance function and user).




		#ifdef PLF_VARIADICS_SUPPORT
			group(aligned_struct_allocator_type &aligned_struct_allocator, const skipfield_type elements_per_group, const group_pointer_type previous):
				next_group(NULL),
				elements(PLF_ALLOCATE(aligned_struct_allocator_type, aligned_struct_allocator, get_aligned_block_capacity(elements_per_group), (previous == NULL) ? NULL : previous->elements)),
				previous_group(previous),
				free_list_head(std::numeric_limits<skipfield_type>::max()),
				capacity(elements_per_group),
				size(1),
				erasures_list_next_group(NULL),
				erasures_list_previous_group(NULL),
				group_number((previous == NULL) ? 0 : previous->group_number + 1u)
			{
				skipfield = pointer_cast<skipfield_pointer_type>(front() + elements_per_group);
				std::memset(plf::void_cast(skipfield), 0, sizeof(skipfield_type) * (static_cast<size_type>(elements_per_group) + 1u));
			}
		#else
			// This is a hack around the fact that allocator_type::construct only supports copy construction in C++03 and copy elision does not occur on the vast majority of compilers in this circumstance. So to avoid running out of memory (and losing performance) from allocating the same block twice, we 'move' in the 'copy' constructor.
			group(aligned_struct_allocator_type &aligned_struct_allocator, const skipfield_type elements_per_group, const group_pointer_type previous) PLF_NOEXCEPT:
				elements(PLF_ALLOCATE(aligned_struct_allocator_type, aligned_struct_allocator, get_aligned_block_capacity(elements_per_group), (previous == NULL) ? 0 : previous->elements)),
				previous_group(previous),
				capacity(elements_per_group)
			{}



			// Not a real copy constructor ie. actually a move constructor. Only used for allocator.construct in C++03 for reasons stated above:
			group(const group &source):
				skipfield(pointer_cast<skipfield_pointer_type>(source.front() + source.capacity)),
				next_group(NULL),
				elements(source.elements),
				previous_group(source.previous_group),
				free_list_head(std::numeric_limits<skipfield_type>::max()),
				capacity(source.capacity),
				size(1),
				erasures_list_next_group(NULL),
				erasures_list_previous_group(NULL),
				group_number((source.previous_group == NULL) ? 0 : source.previous_group->group_number + 1u)
			{
				std::memset(plf::void_cast(skipfield), 0, sizeof(skipfield_type) * (static_cast<size_type>(capacity) + 1u));
			}
		#endif



		void reset(const skipfield_type increment, const group_pointer_type next, const group_pointer_type previous, const size_type group_num) PLF_NOEXCEPT
		{
			next_group = next;
			free_list_head = std::numeric_limits<skipfield_type>::max();
			previous_group = previous;
			size = increment;
			erasures_list_next_group = NULL;
			erasures_list_previous_group = NULL;
			group_number = group_num;

			std::memset(plf::void_cast(skipfield), 0, sizeof(skipfield_type) * static_cast<size_type>(capacity)); // capacity + 1 is not necessary here as the final skipfield node is never written to after initialization
		}



		bool no_erasures() const PLF_NOEXCEPT
		{
			return free_list_head == std::numeric_limits<skipfield_type>::max();
		}



		aligned_pointer_type first_element() const PLF_NOEXCEPT
		{
			return pointer_cast<aligned_pointer_type>(elements) + *skipfield;
		}



		aligned_pointer_type front() const PLF_NOEXCEPT
		{
			return pointer_cast<aligned_pointer_type>(elements);
		}



		aligned_pointer_type past_back() const PLF_NOEXCEPT
		{
			return pointer_cast<aligned_pointer_type>(skipfield);
		}
	};



	// colony member variables:

	iterator 				end_iterator, begin_iterator;
	group_pointer_type	erasure_groups_head,	// Head of doubly-linked list of groups which have erased-element memory locations available for re-use
								unused_groups_head;	// Head of singly-linked list of reserved groups retained by erase()/clear() or created by reserve()
	size_type				total_size, total_capacity;
	skipfield_type 		min_block_capacity, max_block_capacity;

	// Under most compilers, the following will be 1 byte each - since the skipfield types are 2 bytes each at max, the following add nothing to the class sizeof in a 64-bit build, due to padding
	group_allocator_type group_allocator;
	aligned_struct_allocator_type aligned_struct_allocator;
	skipfield_allocator_type skipfield_allocator;
	tuple_allocator_type tuple_allocator;



	void check_capacities_conformance(const plf::limits capacities) const
	{
		PLF_CONSTFUNC plf::limits hard_capacities = block_capacity_hard_limits();

		if (capacities.min < hard_capacities.min || capacities.min > capacities.max || capacities.max > hard_capacities.max)
		{
			#ifdef PLF_EXCEPTIONS_SUPPORT
				throw std::length_error("Supplied memory block capacity limits are either invalid or outside of block_capacity_hard_limits()");
			#else
				std::terminate();
			#endif
		}
	}



	void blank() PLF_NOEXCEPT
	{
		#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT // allocator_traits always available when is_always_equal is available
			if PLF_CONSTEXPR (std::is_standard_layout<colony>::value && std::allocator_traits<allocator_type>::is_always_equal::value && std::is_trivially_destructible<group_pointer_type>::value)
			{ // If all pointer types are trivial, we can just nuke the member variables from orbit with memset (NULL is always 0):
				std::memset(static_cast<void *>(this), 0, offsetof(colony, min_block_capacity));
			}
			else
		#endif
		{
			end_iterator.group_pointer = NULL;
			end_iterator.element_pointer = NULL;
			end_iterator.skipfield_pointer = NULL;
			begin_iterator.group_pointer = NULL;
			begin_iterator.element_pointer = NULL;
			begin_iterator.skipfield_pointer = NULL;
			erasure_groups_head = NULL;
			unused_groups_head = NULL;
			total_size = 0;
			total_capacity = 0;
		}
	}



	static PLF_CONSTFUNC size_t max_size_static() PLF_NOEXCEPT
	{
		#ifdef PLF_ALLOCATOR_TRAITS_SUPPORT
			return static_cast<size_t>(std::allocator_traits<allocator_type>::max_size(allocator_type()));
		#else
			return static_cast<size_t>(std::numeric_limits<size_type>::max() / sizeof(element_type)); // Substituting because allocator::max_size is not static
		#endif
	}



	void reserve_and_fill(const size_type size, const element_type &element)
	{
		if (size != 0)
		{
			reserve(size);
			end_iterator.group_pointer->next_group = unused_groups_head;
			fill_unused_groups(size, element, 0, NULL, begin_iterator.group_pointer);
		}
	}



	template <class iterator_type>
	void reserve_and_range_fill(const size_type size, const iterator_type &it)
	{
		if (size != 0)
		{
			reserve(size);
			end_iterator.group_pointer->next_group = unused_groups_head;
			range_fill_unused_groups(size, it, 0, NULL, begin_iterator.group_pointer);
		}
	}



public:

	// Adaptive minimum based around aligned size, sizeof(group) and sizeof(colony):
	static PLF_CONSTFUNC skipfield_type block_capacity_default_min() PLF_NOEXCEPT
	{
		const skipfield_type adaptive_size = static_cast<skipfield_type>(((sizeof(colony) + sizeof(group)) * 2) / sizeof(aligned_element_struct));
		const skipfield_type max_block_capacity = block_capacity_default_max(); // Necessary to check against in situations with > 64bit pointer sizes and small sizeof(T)
		return std::max(static_cast<skipfield_type>(8), std::min(adaptive_size, max_block_capacity));
	}



	// Adaptive maximum based on numeric_limits and best outcome from multiple benchmark's (on balance) in terms of memory usage and performance:
	static PLF_CONSTFUNC skipfield_type block_capacity_default_max() PLF_NOEXCEPT
	{
		return static_cast<skipfield_type>(std::min(std::min(static_cast<size_t>(std::numeric_limits<skipfield_type>::max()), static_cast<size_t>(8192u)), max_size_static()));
	}



	static PLF_CONSTFUNC plf::limits block_capacity_default_limits() PLF_NOEXCEPT
	{
		return plf::limits(static_cast<size_t>(block_capacity_default_min()), static_cast<size_t>(block_capacity_default_max()));
	}



	// Default constructors:

	PLF_CONSTFUNC explicit colony(const allocator_type &alloc) PLF_NOEXCEPT:
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(block_capacity_default_min()),
		max_block_capacity(block_capacity_default_max()),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{}



	PLF_CONSTFUNC colony() PLF_NOEXCEPT_ALLOCATOR:
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(block_capacity_default_min()),
		max_block_capacity(block_capacity_default_max()),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{}



	PLF_CONSTFUNC colony(const plf::limits block_limits, const allocator_type &alloc):
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(static_cast<skipfield_type>(block_limits.min)),
		max_block_capacity(static_cast<skipfield_type>(block_limits.max)),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		check_capacities_conformance(block_limits);
	}



	PLF_CONSTFUNC explicit colony(const plf::limits block_limits):
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(static_cast<skipfield_type>(block_limits.min)),
		max_block_capacity(static_cast<skipfield_type>(block_limits.max)),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		check_capacities_conformance(block_limits);
	}



	// Copy constructors:
	#ifdef PLF_CPP20_SUPPORT
		colony(const colony &source, const std::type_identity_t<allocator_type> &alloc):
	#else
		colony(const colony &source, const allocator_type &alloc):
	#endif
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(std::max(source.min_block_capacity, static_cast<skipfield_type>(std::min(source.total_size, static_cast<size_type>(source.max_block_capacity))))), // min group size is set to value closest to total number of elements in source colony, in order to not create unnecessary small groups in the range-insert below, then reverts to the original min group size afterwards. This effectively saves a call to reserve.
		max_block_capacity(source.max_block_capacity),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{ // can skip checking for skipfield conformance here as source will have already checked theirs. Same applies for other copy and move constructors below
		reserve_and_range_fill(source.total_size, source.begin_iterator);
		min_block_capacity = source.min_block_capacity; // reset to correct value for future operations
	}



	colony(const colony &source):
		#if (defined(__cplusplus) && __cplusplus >= 201103L) || _MSC_VER >= 1700
			allocator_type(std::allocator_traits<allocator_type>::select_on_container_copy_construction(source)),
		#else
			allocator_type(source),
		#endif
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(std::max(source.min_block_capacity, static_cast<skipfield_type>(std::min(source.total_size, static_cast<size_type>(source.max_block_capacity))))),
		max_block_capacity(source.max_block_capacity),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		reserve_and_range_fill(source.total_size, source.begin_iterator);
		min_block_capacity = source.min_block_capacity; // reset to correct value for future operations
	}



	#ifdef PLF_MOVE_SEMANTICS_SUPPORT
		// Move constructors:

		#ifdef PLF_CPP20_SUPPORT
			colony(colony &&source, const std::type_identity_t<allocator_type> &alloc):
		#else
			colony(colony &&source, const allocator_type &alloc):
		#endif
			allocator_type(alloc),
			end_iterator(source.end_iterator),
			begin_iterator(source.begin_iterator),
			erasure_groups_head(source.erasure_groups_head),
			unused_groups_head(source.unused_groups_head),
			total_size(source.total_size),
			total_capacity(source.total_capacity),
			min_block_capacity(source.min_block_capacity),
			max_block_capacity(source.max_block_capacity),
			group_allocator(alloc),
			aligned_struct_allocator(alloc),
			skipfield_allocator(alloc),
			tuple_allocator(alloc)
		{
			#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT
				if PLF_CONSTEXPR (!std::allocator_traits<allocator_type>::is_always_equal::value)
			#endif
			{
				if (alloc != static_cast<allocator_type &>(source))
				{
					blank();
					reserve_and_range_fill(source.total_size, plf::make_move_iterator(source.begin_iterator));
					source.destroy_all_data();
				}
			}

			source.blank();
		}



		colony(colony &&source) PLF_NOEXCEPT:
			allocator_type(static_cast<allocator_type &>(source)),
			end_iterator(std::move(source.end_iterator)),
			begin_iterator(std::move(source.begin_iterator)),
			erasure_groups_head(std::move(source.erasure_groups_head)),
			unused_groups_head(std::move(source.unused_groups_head)),
			total_size(source.total_size),
			total_capacity(source.total_capacity),
			min_block_capacity(source.min_block_capacity),
			max_block_capacity(source.max_block_capacity),
			group_allocator(*this),
			aligned_struct_allocator(*this),
			skipfield_allocator(*this),
			tuple_allocator(*this)
		{
			source.blank();
		}
	#endif



	// Fill constructors:

	colony(const size_type fill_number, const element_type &element, const plf::limits block_limits, const allocator_type &alloc = allocator_type()):
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(static_cast<skipfield_type>(block_limits.min)),
		max_block_capacity(static_cast<skipfield_type>(block_limits.max)),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		check_capacities_conformance(block_limits);
		reserve_and_fill(fill_number, element);
	}



	colony(const size_type fill_number, const element_type &element, const allocator_type &alloc = allocator_type()) :
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(block_capacity_default_min()),
		max_block_capacity(block_capacity_default_max()),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		reserve_and_fill(fill_number, element);
	}



	// Default-value fill constructors:

	colony(const size_type fill_number, const plf::limits block_limits, const allocator_type &alloc = allocator_type()):
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(static_cast<skipfield_type>(block_limits.min)),
		max_block_capacity(static_cast<skipfield_type>(block_limits.max)),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		check_capacities_conformance(block_limits);
		reserve_and_fill(fill_number, element_type());
	}



	colony(const size_type fill_number, const allocator_type &alloc = allocator_type()):
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(block_capacity_default_min()),
		max_block_capacity(block_capacity_default_max()),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		reserve_and_fill(fill_number, element_type());
	}



	// Range constructors:

	template<typename iterator_type>
	colony(const typename plf::enable_if<!std::numeric_limits<iterator_type>::is_integer, iterator_type>::type &first, const iterator_type &last, const plf::limits block_limits, const allocator_type &alloc = allocator_type()):
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(static_cast<skipfield_type>(block_limits.min)),
		max_block_capacity(static_cast<skipfield_type>(block_limits.max)),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		check_capacities_conformance(block_limits);
		assign<iterator_type>(first, last);
	}



	template<typename iterator_type>
	colony(const typename plf::enable_if<!std::numeric_limits<iterator_type>::is_integer, iterator_type>::type &first, const iterator_type &last, const allocator_type &alloc = allocator_type()):
		allocator_type(alloc),
		erasure_groups_head(NULL),
		unused_groups_head(NULL),
		total_size(0),
		total_capacity(0),
		min_block_capacity(block_capacity_default_min()),
		max_block_capacity(block_capacity_default_max()),
		group_allocator(*this),
		aligned_struct_allocator(*this),
		skipfield_allocator(*this),
		tuple_allocator(*this)
	{
		assign<iterator_type>(first, last);
	}



	#ifdef PLF_INITIALIZER_LIST_SUPPORT
		// Initializer-list constructors:

		colony(const std::initializer_list<element_type> &element_list, const plf::limits block_limits, const allocator_type &alloc = allocator_type()):
			allocator_type(alloc),
			erasure_groups_head(NULL),
			unused_groups_head(NULL),
			total_size(0),
			total_capacity(0),
			min_block_capacity(static_cast<skipfield_type>(block_limits.min)),
			max_block_capacity(static_cast<skipfield_type>(block_limits.max)),
			group_allocator(*this),
			aligned_struct_allocator(*this),
			skipfield_allocator(*this),
			tuple_allocator(*this)
		{
			check_capacities_conformance(block_limits);
			reserve_and_range_fill(static_cast<size_type>(element_list.size()), element_list.begin());
		}



		colony(const std::initializer_list<element_type> &element_list, const allocator_type &alloc = allocator_type()):
			colony(element_list, block_capacity_default_limits(), alloc)
		{}
	#endif



	#ifdef PLF_CPP20_SUPPORT
		// Ranges v3 constructors:

		template<plf::compatible_range<element_type> range_type>
		colony(ranges::from_range_t, range_type &&rg, const plf::limits block_limits, const allocator_type &alloc = allocator_type()):
			allocator_type(alloc),
			erasure_groups_head(NULL),
			unused_groups_head(NULL),
			total_size(0),
			total_capacity(0),
			min_block_capacity(static_cast<skipfield_type>(block_limits.min)),
			max_block_capacity(static_cast<skipfield_type>(block_limits.max)),
			group_allocator(*this),
			aligned_struct_allocator(*this),
			skipfield_allocator(*this),
			tuple_allocator(*this)
		{
			check_capacities_conformance(block_limits);
			reserve_and_range_fill(static_cast<size_type>(std::ranges::distance(rg)), std::ranges::begin(rg));
		}



		template<plf::compatible_range<element_type> range_type>
		colony(plf::ranges::from_range_t, range_type &&rg, const allocator_type &alloc = allocator_type()):
			colony(plf::ranges::from_range, std::move(rg), block_capacity_default_limits(), alloc)
		{}
	#endif



	// Everything else:

	iterator begin() PLF_NOEXCEPT
	{
		return begin_iterator;
	}



	const_iterator begin() const PLF_NOEXCEPT
	{
		return begin_iterator;
	}



	iterator end() PLF_NOEXCEPT
	{
		return end_iterator;
	}



	const_iterator end() const PLF_NOEXCEPT
	{
		return end_iterator;
	}



	const_iterator cbegin() const PLF_NOEXCEPT
	{
		return begin_iterator;
	}



	const_iterator cend() const PLF_NOEXCEPT
	{
		return end_iterator;
	}



	reverse_iterator rbegin() PLF_NOEXCEPT
	{
		return reverse_iterator(end_iterator);
	}



	const_reverse_iterator rbegin() const PLF_NOEXCEPT
	{
		return crbegin();
	}



	reverse_iterator rend() PLF_NOEXCEPT
	{
		return reverse_iterator(begin_iterator);
	}



	const_reverse_iterator rend() const PLF_NOEXCEPT
	{
		return crend();
	}



	const_reverse_iterator crbegin() const PLF_NOEXCEPT
	{
		return const_reverse_iterator(end_iterator);
	}



	const_reverse_iterator crend() const PLF_NOEXCEPT
	{
		return const_reverse_iterator(begin_iterator);
	}



	~colony() PLF_NOEXCEPT
	{
		destroy_all_data();
	}




private:

	group_pointer_type allocate_new_group(const skipfield_type elements_per_group, const group_pointer_type previous = NULL)
	{
		if (max_size() - total_capacity < elements_per_group) // Just in case max_size is a lower amount than the actual memory available (uncommon platform). Comparison avoids overflow.
		{
			#ifdef PLF_EXCEPTIONS_SUPPORT
				throw std::length_error("New block allocation would create capacity greater than max_size()");
			#else
				std::terminate();
			#endif
		}

		const group_pointer_type new_group = PLF_ALLOCATE(group_allocator_type, group_allocator, 1, end_iterator.group_pointer);
		total_capacity += elements_per_group; // I don't know why GCC creates better/smaller codegen when this is placed here rather than at bottom of function, But it does.

		#ifdef PLF_EXCEPTIONS_SUPPORT
			try
			{
				#ifdef PLF_VARIADICS_SUPPORT
					PLF_CONSTRUCT(group_allocator_type, group_allocator, new_group, aligned_struct_allocator, elements_per_group, previous);
				#else
					PLF_CONSTRUCT(group_allocator_type, group_allocator, new_group, group(aligned_struct_allocator, elements_per_group, previous));
				#endif
			}
			catch (...)
			{
				PLF_DEALLOCATE(group_allocator_type, group_allocator, new_group, 1);
				total_capacity -= elements_per_group;
				throw;
			}
		#else
			#ifdef PLF_VARIADICS_SUPPORT
				PLF_CONSTRUCT(group_allocator_type, group_allocator, new_group, aligned_struct_allocator, elements_per_group, previous);
			#else
				PLF_CONSTRUCT(group_allocator_type, group_allocator, new_group, group(aligned_struct_allocator, elements_per_group, previous));
			#endif
		#endif

		return new_group;
	}



	void deallocate_group(const group_pointer_type the_group) PLF_NOEXCEPT
	{
		PLF_DEALLOCATE(aligned_struct_allocator_type, aligned_struct_allocator, the_group->elements, get_aligned_block_capacity(the_group->capacity));
		PLF_DEALLOCATE(group_allocator_type, group_allocator, the_group, 1);
	}



	void deallocate_group_remove_capacity(const group_pointer_type the_group) PLF_NOEXCEPT
	{
		total_capacity -= the_group->capacity;
		deallocate_group(the_group);
	}



	PLF_CONSTFUNC void destroy_element(
	#ifdef PLF_CPP_20_SUPPORT
		[[maybe_unused]]
	#endif
		const aligned_pointer_type element) PLF_NOEXCEPT
	{
		#if defined(PLF_TYPE_TRAITS_SUPPORT) && defined(PLF_CPP20_SUPPORT)  // To avoid codegen for this function with trivially-destructible types. CPP20 because we don't want to trigger a branch for every destruction
			if constexpr (!std::is_trivially_destructible<element_type>::value)
		#endif
		PLF_DESTROY(allocator_type, *this, pointer_cast<pointer>(element));
	}



	void destroy_remainder(
	#ifdef PLF_CPP_20_SUPPORT
		[[maybe_unused]]
	#endif
		const_iterator it) PLF_NOEXCEPT
	{
		#if defined(PLF_TYPE_TRAITS_SUPPORT)
			if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
		#endif
		{
			while (it != end_iterator) destroy_element(it++.element_pointer);
		}
	}



	#ifdef PLF_CPP_20_SUPPORT
		void destroy_group([[maybe_unused]] const_iterator current, [[maybe_unused]] const aligned_pointer_type end) PLF_NOEXCEPT
	#else
		void destroy_group(const_iterator current, const aligned_pointer_type end) PLF_NOEXCEPT
	#endif
	{
		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
		#endif
		{
			do
			{
				destroy_element(current.element_pointer);
				current.element_pointer += static_cast<size_type>(*++current.skipfield_pointer) + 1u;
				current.skipfield_pointer += *current.skipfield_pointer;
			} while(current.element_pointer != end);
		}
	}



	void destroy_dealloc_begin_group(const aligned_pointer_type end) PLF_NOEXCEPT
	{
		destroy_group(begin_iterator, end);
		deallocate_group(begin_iterator.group_pointer);
	}



	void destroy_all_data() PLF_NOEXCEPT
	{
		if (begin_iterator.group_pointer != NULL)
		{
			end_iterator.group_pointer->next_group = unused_groups_head; // Link used and unused_group lists together

			#ifdef PLF_TYPE_TRAITS_SUPPORT
				if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
			#endif
			{
				if (total_size != 0)
				{
					while (begin_iterator.group_pointer != end_iterator.group_pointer) // Erase elements without bothering to update skipfield - much faster:
					{
						const group_pointer_type next_group = begin_iterator.group_pointer->next_group;
						destroy_dealloc_begin_group(begin_iterator.group_pointer->past_back());
						begin_iterator.group_pointer = next_group;
						begin_iterator.set_to_first_element_in_group();
					}

					destroy_dealloc_begin_group(end_iterator.element_pointer);
					begin_iterator.group_pointer = unused_groups_head;
				}
			}

			while (begin_iterator.group_pointer != NULL)
			{
				const group_pointer_type next_group = begin_iterator.group_pointer->next_group;
				deallocate_group(begin_iterator.group_pointer);
				begin_iterator.group_pointer = next_group;
			}
		}
	}



	void initialize(const skipfield_type first_group_size)
	{
		end_iterator.group_pointer = begin_iterator.group_pointer = allocate_new_group(first_group_size);
		end_iterator.element_pointer = begin_iterator.element_pointer = begin_iterator.group_pointer->front();
		end_iterator.skipfield_pointer = begin_iterator.skipfield_pointer = begin_iterator.group_pointer->skipfield;
	}



	void edit_free_list(const skipfield_pointer_type location, const skipfield_type value) PLF_NOEXCEPT
	{
		PLF_DESTROY(skipfield_allocator_type, skipfield_allocator, location);
		PLF_CONSTRUCT(skipfield_allocator_type, skipfield_allocator, location, value);
	}



	void edit_free_list_prev(const aligned_pointer_type location, const skipfield_type value) PLF_NOEXCEPT // Write to the 'previous erased element' index in the erased element memory location
	{
		edit_free_list(pointer_cast<skipfield_pointer_type>(location), value);
	}



	void edit_free_list_next(const aligned_pointer_type location, const skipfield_type value) PLF_NOEXCEPT // Ditto 'next'
	{
		edit_free_list(pointer_cast<skipfield_pointer_type>(location) + 1, value);
	}



	void edit_free_list_head(const aligned_pointer_type location, const skipfield_type value) PLF_NOEXCEPT
	{
		const skipfield_pointer_type converted_location = pointer_cast<skipfield_pointer_type>(location);
		edit_free_list(converted_location, value);
		edit_free_list(converted_location + 1, std::numeric_limits<skipfield_type>::max());
	}



	void update_skipblock(const iterator &new_location, const skipfield_type prev_free_list_index) PLF_NOEXCEPT
	{
		const skipfield_type new_value = static_cast<skipfield_type>(*(new_location.skipfield_pointer) - 1);

		if (new_value != 0) // ie. skipfield was not originally length 1, hence we need to truncate it
		{
			// set (new) start and (original) end of skipblock to new value:
			*(new_location.skipfield_pointer + new_value) = *(new_location.skipfield_pointer + 1) = new_value;

			// transfer free list node to new start node:
			++(erasure_groups_head->free_list_head);

			if (prev_free_list_index != std::numeric_limits<skipfield_type>::max()) // ie. not the tail free list node
			{
				edit_free_list_next(new_location.group_pointer->front() + prev_free_list_index, erasure_groups_head->free_list_head);
			}

			edit_free_list_head(new_location.element_pointer + 1, prev_free_list_index);
		}
		else // single-node skipblock, remove skipblock
		{
			erasure_groups_head->free_list_head = prev_free_list_index;

			if (prev_free_list_index != std::numeric_limits<skipfield_type>::max()) // ie. not the last free list node
			{
				edit_free_list_next(new_location.group_pointer->front() + prev_free_list_index, std::numeric_limits<skipfield_type>::max());
			}
			else // remove this group from the list of groups with erasures
			{
				erasure_groups_head = erasure_groups_head->erasures_list_next_group; // update_skipblock is only used within insert/emplace, where the head group is being used, so no need for additional checks here
			}
		}

		*(new_location.skipfield_pointer) = 0;
		++(new_location.group_pointer->size);

		if (new_location.group_pointer == begin_iterator.group_pointer && new_location.element_pointer < begin_iterator.element_pointer)
		{ /* ie. begin_iterator was moved forwards as the result of an erasure at some point, this erased element is before the current begin, hence, set current begin iterator to this element */
			begin_iterator = new_location;
		}

		++total_size;
	}



	void update_subsequent_group_numbers(size_type current_group_number, group_pointer_type update_group) PLF_NOEXCEPT
	{
		do
		{
			update_group->group_number = current_group_number++;
			update_group = update_group->next_group;
		} while (update_group != NULL);
	}



	void reset_group_numbers() PLF_NOEXCEPT
	{
		update_subsequent_group_numbers(0, begin_iterator.group_pointer);
	}



	void reset_group_numbers_if_necessary() PLF_NOEXCEPT
	{
		if (end_iterator.group_pointer->group_number == std::numeric_limits<size_type>::max())
		#ifdef PLF_CPP20_SUPPORT
			[[unlikely]]
		#endif
		{
			reset_group_numbers();
		}
	}



	group_pointer_type reuse_unused_group() PLF_NOEXCEPT
	{
		const group_pointer_type reused_group = unused_groups_head;
		unused_groups_head = reused_group->next_group;
		reused_group->reset(1, NULL, end_iterator.group_pointer, end_iterator.group_pointer->group_number + 1u);
		return reused_group;
	}



public:


	void reset() PLF_NOEXCEPT
	{
		destroy_all_data();
		blank();
	}



	iterator insert(const element_type &element)
	{
		if (end_iterator.element_pointer != NULL)
		{
			if (erasure_groups_head == NULL) // ie. there are no erased elements
			{
				if (end_iterator.element_pointer != end_iterator.group_pointer->past_back())
				{
					PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, element);

					const iterator return_iterator = end_iterator;
					++end_iterator.element_pointer;
					++end_iterator.skipfield_pointer;
					++(end_iterator.group_pointer->size);
					++total_size;
					return return_iterator;
				}
				else
				{
					reset_group_numbers_if_necessary();
					group_pointer_type next_group;

					if (unused_groups_head == NULL)
					{
						next_group = allocate_new_group(static_cast<skipfield_type>(std::min(total_size, static_cast<size_type>(max_block_capacity))), end_iterator.group_pointer);

						#ifndef PLF_EXCEPTIONS_SUPPORT
							PLF_CONSTRUCT_ELEMENT(next_group->elements, element);
						#else
							#ifdef PLF_TYPE_TRAITS_SUPPORT
								if PLF_CONSTEXPR (std::is_nothrow_copy_constructible<element_type>::value)
								{
									PLF_CONSTRUCT_ELEMENT(next_group->elements, element);
								}
								else
							#endif
							{
								try
								{
									PLF_CONSTRUCT_ELEMENT(next_group->elements, element);
								}
								catch (...)
								{
									deallocate_group_remove_capacity(next_group);
									throw;
								}
							}
						#endif
					}
					else
					{
						PLF_CONSTRUCT_ELEMENT(unused_groups_head->elements, element);
						next_group = reuse_unused_group();
					}

					end_iterator.group_pointer->next_group = next_group;
					end_iterator.group_pointer = next_group;
					end_iterator.element_pointer = next_group->front() + 1;
					end_iterator.skipfield_pointer = next_group->skipfield + 1;
					++total_size;

					return iterator(next_group, next_group->front(), next_group->skipfield);
				}
			}
			else // there are erased elements, reuse those memory locations
			{
				const iterator new_location(erasure_groups_head, erasure_groups_head->front() + erasure_groups_head->free_list_head, erasure_groups_head->skipfield + erasure_groups_head->free_list_head);

				// We always reuse the element at the start of the skipblock, this is also where the free-list information for that skipblock is stored. Get the previous free-list node's index from this memory space, before we write to our element to it. 'Next' index is always the free_list_head (as represented by the maximum value of the skipfield type) here so we don't need to get it:
				const skipfield_type prev_free_list_index = *pointer_cast<skipfield_pointer_type>(new_location.element_pointer);
				PLF_CONSTRUCT_ELEMENT(new_location.element_pointer, element);
				update_skipblock(new_location, prev_free_list_index);

				return new_location;
			}
		}
		else // ie. newly-constructed colony, no insertions yet and no groups
		{
			initialize(min_block_capacity);

			#ifndef PLF_EXCEPTIONS_SUPPORT
				PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, element);
			#else
				#ifdef PLF_TYPE_TRAITS_SUPPORT
					if PLF_CONSTEXPR (std::is_nothrow_copy_constructible<element_type>::value)
					{
						PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, element);
					}
					else
				#endif
				{
					try
					{
						PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, element);
					}
					catch (...)
					{
						reset();
						throw;
					}
				}
			#endif

			++end_iterator.skipfield_pointer;
			total_size = 1;
			return begin_iterator;
		}
	}



	#ifdef PLF_CPP20_SUPPORT
		iterator insert([[maybe_unused]] const_iterator &hint, const element_type &element) // Note: hint is ignored, purely to serve other standard library functions
		{
			return insert(element);
		}
	#endif



	#ifdef PLF_MOVE_SEMANTICS_SUPPORT
		iterator insert(element_type &&element) // The move-insert function is near-identical to the regular insert function, with the exception of the element construction method and is_nothrow tests.
		{
			if (end_iterator.element_pointer != NULL)
			{
				if (erasure_groups_head == NULL)
				{
					if (end_iterator.element_pointer != end_iterator.group_pointer->past_back())
					{
						PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, std::move(element));

						const iterator return_iterator = end_iterator;
						++end_iterator.element_pointer;
						++end_iterator.skipfield_pointer;
						++(end_iterator.group_pointer->size);
						++total_size;

						return return_iterator;
					}
					else
					{
						reset_group_numbers_if_necessary();
						group_pointer_type next_group;

						if (unused_groups_head == NULL)
						{
							next_group = allocate_new_group(static_cast<skipfield_type>(std::min(total_size, static_cast<size_type>(max_block_capacity))), end_iterator.group_pointer);

							#ifndef PLF_EXCEPTIONS_SUPPORT
								PLF_CONSTRUCT_ELEMENT(next_group->elements, std::move(element));
							#else
								#ifdef PLF_TYPE_TRAITS_SUPPORT
									if PLF_CONSTEXPR (std::is_nothrow_move_constructible<element_type>::value)
									{
										PLF_CONSTRUCT_ELEMENT(next_group->elements, std::move(element));
									}
									else
								#endif
								{
									try
									{
										PLF_CONSTRUCT_ELEMENT(next_group->elements, std::move(element));
									}
									catch (...)
									{
										deallocate_group_remove_capacity(next_group);
										throw;
									}
								}
							#endif
						}
						else
						{
							PLF_CONSTRUCT_ELEMENT(unused_groups_head->elements, std::move(element));
							next_group = reuse_unused_group();
						}

						end_iterator.group_pointer->next_group = next_group;
						end_iterator.group_pointer = next_group;
						end_iterator.element_pointer = next_group->front() + 1;
						end_iterator.skipfield_pointer = next_group->skipfield + 1;
						++total_size;

						return iterator(next_group, next_group->front(), next_group->skipfield);
					}
				}
				else
				{
					const iterator new_location(erasure_groups_head, erasure_groups_head->front() + erasure_groups_head->free_list_head, erasure_groups_head->skipfield + erasure_groups_head->free_list_head);

					const skipfield_type prev_free_list_index = *pointer_cast<skipfield_pointer_type>(new_location.element_pointer);
					PLF_CONSTRUCT_ELEMENT(new_location.element_pointer, std::move(element));
					update_skipblock(new_location, prev_free_list_index);

					return new_location;
				}
			}
			else
			{
				initialize(min_block_capacity);

				#ifndef PLF_EXCEPTIONS_SUPPORT
					PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, std::move(element));
				#else
					#ifdef PLF_TYPE_TRAITS_SUPPORT
						if PLF_CONSTEXPR (std::is_nothrow_move_constructible<element_type>::value)
						{
							PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, std::move(element));
						}
						else
					#endif
					{
						try
						{
							PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, std::move(element));
						}
						catch (...)
						{
							reset();
							throw;
						}
					}
				#endif

				++end_iterator.skipfield_pointer;
				total_size = 1;
				return begin_iterator;
			}
		}



		#ifdef PLF_CPP20_SUPPORT
			iterator insert([[maybe_unused]] const_iterator &hint, element_type &&element)
			{
				return insert(std::forward<element_type &&>(element));
			}
		#endif
	#endif



	#ifdef PLF_VARIADICS_SUPPORT
		template<typename... arguments>
		iterator emplace(arguments &&... parameters) // The emplace function is near-identical to the regular insert function, with the exception of the element construction method, and change to is_nothrow tests.
		{
			if (end_iterator.element_pointer != NULL)
			{
				if (erasure_groups_head == NULL)
				{
					if (end_iterator.element_pointer != end_iterator.group_pointer->past_back())
					{
						PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, std::forward<arguments>(parameters) ...);

						const iterator return_iterator = end_iterator;
						++end_iterator.element_pointer;
						++end_iterator.skipfield_pointer;
						++(end_iterator.group_pointer->size);
						++total_size;
						return return_iterator;
					}

					reset_group_numbers_if_necessary();
					group_pointer_type next_group;

					if (unused_groups_head == NULL)
					{
						next_group = allocate_new_group(static_cast<skipfield_type>(std::min(total_size, static_cast<size_type>(max_block_capacity))), end_iterator.group_pointer);

						#ifndef PLF_EXCEPTIONS_SUPPORT
							PLF_CONSTRUCT_ELEMENT(next_group->elements, std::forward<arguments>(parameters) ...);
						#else
							#ifdef PLF_TYPE_TRAITS_SUPPORT
								if PLF_CONSTEXPR (std::is_nothrow_constructible<element_type, arguments...>::value)
								{
									PLF_CONSTRUCT_ELEMENT(next_group->elements, std::forward<arguments>(parameters) ...);
								}
								else
							#endif
							{
								try
								{
									PLF_CONSTRUCT_ELEMENT(next_group->elements, std::forward<arguments>(parameters) ...);
								}
								catch (...)
								{
									deallocate_group_remove_capacity(next_group);
									throw;
								}
							}
						#endif
					}
					else
					{
						PLF_CONSTRUCT_ELEMENT(unused_groups_head->elements, std::forward<arguments>(parameters) ...);
						next_group = reuse_unused_group();
					}

					end_iterator.group_pointer->next_group = next_group;
					end_iterator.group_pointer = next_group;
					end_iterator.element_pointer = next_group->front() + 1;
					end_iterator.skipfield_pointer = next_group->skipfield + 1;
					++total_size;

					return iterator(next_group, next_group->front(), next_group->skipfield);
				}
				else
				{
					const iterator new_location(erasure_groups_head, erasure_groups_head->front() + erasure_groups_head->free_list_head, erasure_groups_head->skipfield + erasure_groups_head->free_list_head);

					const skipfield_type prev_free_list_index = *pointer_cast<skipfield_pointer_type>(new_location.element_pointer);
					PLF_CONSTRUCT_ELEMENT(new_location.element_pointer, std::forward<arguments>(parameters) ...);
					update_skipblock(new_location, prev_free_list_index);

					return new_location;
				}
			}
			else
			{
				initialize(min_block_capacity);

				#ifndef PLF_EXCEPTIONS_SUPPORT
					PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, std::forward<arguments>(parameters) ...);
				#else
					#ifdef PLF_TYPE_TRAITS_SUPPORT
						if PLF_CONSTEXPR (std::is_nothrow_constructible<element_type, arguments...>::value)
						{
							PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, std::forward<arguments>(parameters) ...);
						}
						else
					#endif
					{
						try
						{
							PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer++, std::forward<arguments>(parameters) ...);
						}
						catch (...)
						{
							reset();
							throw;
						}
					}
				#endif

				++end_iterator.skipfield_pointer;
				total_size = 1;
				return begin_iterator;
			}
		}



		#ifdef PLF_CPP20_SUPPORT
			template<typename... arguments>
			iterator emplace_hint([[maybe_unused]] const_iterator &hint, arguments &&... parameters)
			{
				return emplace(std::forward<arguments>(parameters) ...);
			}
		#endif
	#endif




private:

	// For catch blocks in fill() and range_fill()
	void recover_from_partial_fill()
	{
		#ifdef PLF_EXCEPTIONS_SUPPORT
			#ifdef PLF_TYPE_TRAITS_SUPPORT
				if PLF_CONSTEXPR ((!std::is_copy_constructible<element_type>::value && !std::is_nothrow_move_constructible<element_type>::value) || !std::is_nothrow_copy_constructible<element_type>::value) // to avoid unnecessary codegen, since this function will never be called if this line isn't true
			#endif
			{
				const skipfield_type elements_constructed_before_exception = static_cast<skipfield_type>(end_iterator.element_pointer - end_iterator.group_pointer->front());
				end_iterator.group_pointer->size = elements_constructed_before_exception;
				end_iterator.skipfield_pointer = end_iterator.group_pointer->skipfield + elements_constructed_before_exception;
				total_size += elements_constructed_before_exception;
				unused_groups_head = end_iterator.group_pointer->next_group;
				end_iterator.group_pointer->next_group = NULL;
			}
		#endif
	}



	void fill(const element_type &element, const skipfield_type size)
	{
		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (std::is_nothrow_copy_constructible<element_type>::value)
			{
				if PLF_CONSTEXPR (sizeof(aligned_element_struct) != sizeof(element_type))
				{
					alignas (alignof(aligned_element_struct)) element_type aligned_copy = element; // to avoid potentially violating memory boundaries in line below, create an initial object copy of same (but aligned) type
					plf::uninitialized_fill_n(end_iterator.element_pointer, size, *pointer_cast<aligned_pointer_type>(&aligned_copy), static_cast<allocator_type &>(*this));
				}
				else
				{
					plf::uninitialized_fill_n(pointer_cast<pointer>(end_iterator.element_pointer), size, element, static_cast<allocator_type &>(*this));
				}

				end_iterator.element_pointer += size;
			}
			else
		#endif
		{
			const aligned_pointer_type fill_end = end_iterator.element_pointer + size;

			do
			{
				#ifdef PLF_EXCEPTIONS_SUPPORT
					try
					{
						PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, element);
					}
					catch (...)
					{
						recover_from_partial_fill();
						throw;
					}
				#else
					PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, element);
				#endif
			} while (++end_iterator.element_pointer != fill_end);
		}

		total_size += size;
	}



	// For catch blocks in range_fill_skipblock and fill_skipblock - update existing skipblock and free-list indexes to reflect partially-reused skipblock:
	#ifdef PLF_CPP20_SUPPORT
		void recover_from_partial_skipblock_fill([[maybe_unused]] const aligned_pointer_type location, [[maybe_unused]] const aligned_pointer_type current_location, [[maybe_unused]] const skipfield_pointer_type skipfield_pointer, [[maybe_unused]] const skipfield_type prev_free_list_node)
	#else
		void recover_from_partial_skipblock_fill(const aligned_pointer_type location, const aligned_pointer_type current_location, const skipfield_pointer_type skipfield_pointer, const skipfield_type prev_free_list_node)
	#endif
	{
		#ifdef PLF_EXCEPTIONS_SUPPORT // to avoid unnecessary codegen
			#ifdef PLF_TYPE_TRAITS_SUPPORT
				if PLF_CONSTEXPR ((!std::is_copy_constructible<element_type>::value && !std::is_nothrow_move_constructible<element_type>::value) || !std::is_nothrow_copy_constructible<element_type>::value)
			#endif
			{
				const skipfield_type elements_constructed_before_exception = static_cast<skipfield_type>(current_location - location);
				erasure_groups_head->size += elements_constructed_before_exception;
				total_size += elements_constructed_before_exception;

				// Update skipblock:
				const skipfield_type new_start_node_value = *skipfield_pointer - elements_constructed_before_exception;
				const skipfield_pointer_type new_start_node = skipfield_pointer + elements_constructed_before_exception;
				std::memset(plf::void_cast(skipfield_pointer), 0, elements_constructed_before_exception * sizeof(skipfield_type)); // Reset skipfield for elements written before exception
				*new_start_node = *(new_start_node + new_start_node_value - 1) = new_start_node_value; // Create new skipblock for unused elements

				// Update free list of erased elements:
				edit_free_list_head(location + elements_constructed_before_exception, prev_free_list_node);

				const skipfield_type new_skipblock_head_index = static_cast<skipfield_type>(current_location - erasure_groups_head->front());
				erasure_groups_head->free_list_head = new_skipblock_head_index;

				if (prev_free_list_node != std::numeric_limits<skipfield_type>::max())
				{
					edit_free_list_next(erasure_groups_head->front() + prev_free_list_node, new_skipblock_head_index);
				}
			}
		#endif
	}



	void fill_skipblock(const element_type &element, const aligned_pointer_type location, const skipfield_pointer_type skipfield_pointer, const skipfield_type size)
	{
		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (std::is_nothrow_copy_constructible<element_type>::value)
			{
				if PLF_CONSTEXPR (sizeof(aligned_element_struct) != sizeof(element_type))
				{
					alignas (alignof(aligned_element_struct)) element_type aligned_copy = element;
					plf::uninitialized_fill_n(location, size, *pointer_cast<aligned_pointer_type>(&aligned_copy), static_cast<allocator_type &>(*this));
				}
				else
				{
					plf::uninitialized_fill_n(pointer_cast<pointer>(location), size, element, static_cast<allocator_type &>(*this));
				}
			}
			else
		#endif
		{
			const aligned_pointer_type fill_end = location + size;
			#ifdef PLF_EXCEPTIONS_SUPPORT
				const skipfield_type prev_free_list_node = *pointer_cast<skipfield_pointer_type>(location); // in case of exception, grabbing indexes before free_list node is reused
			#endif

			for (aligned_pointer_type current_location = location; current_location != fill_end; ++current_location)
			{
				#ifdef PLF_EXCEPTIONS_SUPPORT
					try
					{
						PLF_CONSTRUCT_ELEMENT(current_location, element);
					}
					catch (...)
					{
						recover_from_partial_skipblock_fill(location, current_location, skipfield_pointer, prev_free_list_node);
						throw;
					}
				#else
					PLF_CONSTRUCT_ELEMENT(current_location, element);
				#endif
			}
		}

		std::memset(plf::void_cast(skipfield_pointer), 0, size * sizeof(skipfield_type)); // reset skipfield nodes within skipblock to 0
		erasure_groups_head->size += size;
		total_size += size;
	}



	void fill_unused_groups(size_type size, const element_type &element, size_type group_number, group_pointer_type previous_group, const group_pointer_type current_group)
	{
		for (end_iterator.group_pointer = current_group; end_iterator.group_pointer->capacity < size; end_iterator.group_pointer = end_iterator.group_pointer->next_group)
		{
			const skipfield_type capacity = end_iterator.group_pointer->capacity;
			end_iterator.group_pointer->reset(capacity, end_iterator.group_pointer->next_group, previous_group, group_number++);
			previous_group = end_iterator.group_pointer;
			size -= static_cast<size_type>(capacity);
			end_iterator.element_pointer = end_iterator.group_pointer->front();
			fill(element, capacity);
		}

		// Deal with final group (partial fill)
		unused_groups_head = end_iterator.group_pointer->next_group;
		end_iterator.group_pointer->reset(static_cast<skipfield_type>(size), NULL, previous_group, group_number);
		end_iterator.element_pointer = end_iterator.group_pointer->front();
		end_iterator.skipfield_pointer = end_iterator.group_pointer->skipfield + size;
		fill(element, static_cast<skipfield_type>(size));
	}



public:

	// Fill insert

	void insert(size_type size, const element_type &element)
	{
		if (size == 0)
		{
			return;
		}
		else if (size == 1)
		{
			insert(element);
			return;
		}

		if (total_size == 0)
		{
			prepare_groups_for_assign(size);
			fill_unused_groups(size, element, 0, NULL, begin_iterator.group_pointer);
			return;
		}

		reserve(total_size + size);

		// Use up erased locations if available:
		while(erasure_groups_head != NULL) // skipblock loop: breaks when colony is exhausted of reusable skipblocks, or returns if size == 0
		{
			const aligned_pointer_type element_pointer = erasure_groups_head->front() + erasure_groups_head->free_list_head;
			const skipfield_pointer_type skipfield_pointer = erasure_groups_head->skipfield + erasure_groups_head->free_list_head;
			const skipfield_type skipblock_size = *skipfield_pointer;

			if (erasure_groups_head == begin_iterator.group_pointer && element_pointer < begin_iterator.element_pointer)
			{
				begin_iterator.element_pointer = element_pointer;
				begin_iterator.skipfield_pointer = skipfield_pointer;
			}

			if (skipblock_size <= size)
			{
				erasure_groups_head->free_list_head = *pointer_cast<skipfield_pointer_type>(element_pointer); // set free list head to previous free list node
				fill_skipblock(element, element_pointer, skipfield_pointer, skipblock_size);
				size -= skipblock_size;

				if (!erasure_groups_head->no_erasures()) // ie. there are more skipblocks to be filled in this group
				{
					edit_free_list_next(erasure_groups_head->front() + erasure_groups_head->free_list_head, std::numeric_limits<skipfield_type>::max()); // set 'next' index of new free list head to 'end' (numeric max)
				}
				else
				{
					erasure_groups_head = erasure_groups_head->erasures_list_next_group; // change groups
				}

				if (size == 0) return;
			}
			else // skipblock is larger than remaining number of elements
			{
				const skipfield_type prev_index = *pointer_cast<skipfield_pointer_type>(element_pointer); // save before element location is overwritten
				fill_skipblock(element, element_pointer, skipfield_pointer, static_cast<skipfield_type>(size));
				const skipfield_type new_skipblock_size = static_cast<skipfield_type>(skipblock_size - size);

				// Update skipfield (earlier nodes already memset'd in fill_skipblock function):
				*(skipfield_pointer + size) = new_skipblock_size;
				*(skipfield_pointer + skipblock_size - 1) = new_skipblock_size;
				erasure_groups_head->free_list_head += static_cast<skipfield_type>(size); // set free list head to new start node

				// Update free list with new head:
				edit_free_list_head(element_pointer + size, prev_index);

				if (prev_index != std::numeric_limits<skipfield_type>::max())
				{
					edit_free_list_next(erasure_groups_head->front() + prev_index,  erasure_groups_head->free_list_head); // set 'next' index of previous skipblock to new start of skipblock
				}

				return;
			}
		}


		// Use up remaining available element locations in end group:
		// This variable is either the remaining capacity of the group or the number of elements yet to be inserted, whichever is smaller:
		const skipfield_type group_remainder = static_cast<skipfield_type>(std::min(static_cast<size_type>(end_iterator.group_pointer->past_back() - end_iterator.element_pointer), size));

		if (group_remainder != 0)
		{
			fill(element, group_remainder);
			end_iterator.group_pointer->size += group_remainder;

			if (size == group_remainder) // ie. remaining capacity was >= remaining elements to be filled
			{
				end_iterator.skipfield_pointer = end_iterator.group_pointer->skipfield + end_iterator.group_pointer->size;
				return;
			}

			size -= group_remainder;
		}


		// Use unused groups:
		end_iterator.group_pointer->next_group = unused_groups_head;

		if ((std::numeric_limits<size_type>::max() - end_iterator.group_pointer->group_number) < size)
		#ifdef PLF_CPP20_SUPPORT
			[[unlikely]]
		#endif
		{
			reset_group_numbers();
		}

		fill_unused_groups(size, element, end_iterator.group_pointer->group_number + 1u, end_iterator.group_pointer, unused_groups_head);
	}



private:

	template <class iterator_type>
	void range_fill(iterator_type &it, const skipfield_type size)
	{
		const aligned_pointer_type fill_end = end_iterator.element_pointer + size;

		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (std::is_nothrow_copy_constructible<element_type>::value)
			{
				do
				{
					PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, *it++);
				} while (++end_iterator.element_pointer != fill_end);
			}
			#ifdef PLF_MOVE_SEMANTICS_SUPPORT
				else if PLF_CONSTEXPR (std::is_nothrow_move_constructible<element_type>::value && !std::is_copy_constructible<element_type>::value)
				{
					do
					{
						PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, std::move(*it++));
					} while (++end_iterator.element_pointer != fill_end);
				}
			#endif
			else
		#endif
		{
			do
			{
				#ifdef PLF_EXCEPTIONS_SUPPORT
					try
					{
				#endif
					#ifdef PLF_TYPE_TRAITS_SUPPORT
						if PLF_CONSTEXPR (!std::is_copy_constructible<element_type>::value)
						{
							PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, std::move(*it++));
						}
						else
					#endif
						PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, *it++);
				#ifdef PLF_EXCEPTIONS_SUPPORT
					}
					catch (...)
					{
						recover_from_partial_fill();
						throw;
					}
				#endif
			} while (++end_iterator.element_pointer != fill_end);
		}

		total_size += size;
	}



	template <class iterator_type>
	void range_fill_skipblock(iterator_type &it, const aligned_pointer_type location, const skipfield_pointer_type skipfield_pointer, const skipfield_type size)
	{
		const aligned_pointer_type fill_end = location + size;

		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (std::is_nothrow_copy_constructible<element_type>::value)
			{
				for (aligned_pointer_type current_location = location; current_location != fill_end; ++current_location)
				{
					PLF_CONSTRUCT_ELEMENT(current_location, *it++);
				}
			}
			#ifdef PLF_MOVE_SEMANTICS_SUPPORT
				else if PLF_CONSTEXPR (std::is_nothrow_move_constructible<element_type>::value && !std::is_copy_constructible<element_type>::value)
				{
					for (aligned_pointer_type current_location = location; current_location != fill_end; ++current_location)
					{
						PLF_CONSTRUCT_ELEMENT(current_location, std::move(*it++));
					}
				}
			#endif
			else
		#endif
		{
			#ifdef PLF_EXCEPTIONS_SUPPORT
				const skipfield_type prev_free_list_node = *pointer_cast<skipfield_pointer_type>(location); // in case of exception, grabbing indexes before free_list node is reused
			#endif

			for (aligned_pointer_type current_location = location; current_location != fill_end; ++current_location)
			{
				#ifdef PLF_EXCEPTIONS_SUPPORT
					try
					{
				#endif
					#ifdef PLF_TYPE_TRAITS_SUPPORT
						if PLF_CONSTEXPR (!std::is_copy_constructible<element_type>::value)
						{
							PLF_CONSTRUCT_ELEMENT(current_location, std::move(*it++));
						}
						else
					#endif
					PLF_CONSTRUCT_ELEMENT(current_location, *it++);
				#ifdef PLF_EXCEPTIONS_SUPPORT
					}
					catch (...)
					{
						recover_from_partial_skipblock_fill(location, current_location, skipfield_pointer, prev_free_list_node);
						throw;
					}
				#endif
			}
		}

		std::memset(plf::void_cast(skipfield_pointer), 0, size * sizeof(skipfield_type)); // reset skipfield nodes within skipblock to 0
		erasure_groups_head->size += size;
		total_size += size;
	}



	template <class iterator_type>
	void range_fill_unused_groups(size_type size, iterator_type it, size_type group_number, group_pointer_type previous_group, const group_pointer_type current_group)
	{
		for (end_iterator.group_pointer = current_group; end_iterator.group_pointer->capacity < size; end_iterator.group_pointer = end_iterator.group_pointer->next_group)
		{
			const skipfield_type capacity = end_iterator.group_pointer->capacity;
			end_iterator.group_pointer->reset(capacity, end_iterator.group_pointer->next_group, previous_group, group_number++);
			previous_group = end_iterator.group_pointer;
			size -= static_cast<size_type>(capacity);
			end_iterator.element_pointer = end_iterator.group_pointer->front();
			range_fill(it, capacity);
		}

		// Deal with final group (partial fill)
		unused_groups_head = end_iterator.group_pointer->next_group;
		end_iterator.group_pointer->reset(static_cast<skipfield_type>(size), NULL, previous_group, group_number);
		end_iterator.element_pointer = end_iterator.group_pointer->front();
		end_iterator.skipfield_pointer = end_iterator.group_pointer->skipfield + size;
		range_fill(it, static_cast<skipfield_type>(size));
	}



	template <class iterator_type>
	void range_insert (iterator_type it, size_type size) // this is near-identical to the fill insert, with the only alteration being incrementing an iterator for construction, rather than using a const element. And the fill etc function calls are changed to range_fill to match this pattern. See fill insert for code explanations
	{
		if (size == 0)
		{
			return;
		}
		else if (size == 1)
		{
			insert(*it);
			return;
		}

		if (total_size == 0)
		{
			prepare_groups_for_assign(size);
			range_fill_unused_groups(size, it, 0, NULL, begin_iterator.group_pointer);
			return;
		}

		reserve(total_size + size);

		while(erasure_groups_head != NULL)
		{
			const aligned_pointer_type element_pointer = erasure_groups_head->front() + erasure_groups_head->free_list_head;
			const skipfield_pointer_type skipfield_pointer = erasure_groups_head->skipfield + erasure_groups_head->free_list_head;
			const skipfield_type skipblock_size = *skipfield_pointer;

			if (erasure_groups_head == begin_iterator.group_pointer && element_pointer < begin_iterator.element_pointer)
			{
				begin_iterator.element_pointer = element_pointer;
				begin_iterator.skipfield_pointer = skipfield_pointer;
			}

			if (skipblock_size <= size)
			{
				erasure_groups_head->free_list_head = *pointer_cast<skipfield_pointer_type>(element_pointer);
				range_fill_skipblock(it, element_pointer, skipfield_pointer, skipblock_size);
				size -= skipblock_size;

				if (!erasure_groups_head->no_erasures())
				{
					edit_free_list_next(erasure_groups_head->front() + erasure_groups_head->free_list_head, std::numeric_limits<skipfield_type>::max());
				}
				else
				{
					erasure_groups_head = erasure_groups_head->erasures_list_next_group;
				}

				if (size == 0) return;
			}
			else
			{
				const skipfield_type prev_index = *pointer_cast<skipfield_pointer_type>(element_pointer);
				range_fill_skipblock(it, element_pointer, skipfield_pointer, static_cast<skipfield_type>(size));
				const skipfield_type new_skipblock_size = static_cast<skipfield_type>(skipblock_size - size);

				*(skipfield_pointer + size) = new_skipblock_size;
				*(skipfield_pointer + skipblock_size - 1) = new_skipblock_size;
				erasure_groups_head->free_list_head += static_cast<skipfield_type>(size);
				edit_free_list_head(element_pointer + size, prev_index);

				if (prev_index != std::numeric_limits<skipfield_type>::max())
				{
					edit_free_list_next(erasure_groups_head->front() + prev_index, erasure_groups_head->free_list_head);
				}

				return;
			}
		}

		const skipfield_type group_remainder = static_cast<skipfield_type>(std::min(static_cast<size_type>(end_iterator.group_pointer->past_back() - end_iterator.element_pointer), size));

		if (group_remainder != 0)
		{
			range_fill(it, group_remainder);
			end_iterator.group_pointer->size += group_remainder;

			if (size == group_remainder)
			{
				end_iterator.skipfield_pointer = end_iterator.group_pointer->skipfield + end_iterator.group_pointer->size;
				return;
			}

			size -= group_remainder;
		}


		end_iterator.group_pointer->next_group = unused_groups_head;

		if ((std::numeric_limits<size_type>::max() - end_iterator.group_pointer->group_number) < size)
		#ifdef PLF_CPP20_SUPPORT
			[[unlikely]]
		#endif
		{
			reset_group_numbers();
		}

		range_fill_unused_groups(size, it, end_iterator.group_pointer->group_number + 1u, end_iterator.group_pointer, unused_groups_head);
	}



public:

	// Range insert:

	template <class iterator_type>
	void insert (const typename plf::enable_if<!std::numeric_limits<iterator_type>::is_integer, iterator_type>::type first, const iterator_type last)
	{
		range_insert(first, static_cast<size_type>(std::distance(first, last)));
	}



	template <class iterator_type>
	void insert (const iterator &first, const iterator &last)
	{
		range_insert(first, static_cast<size_type>(first.distance(last)));
	}



	template <class iterator_type>
	void insert (const const_iterator &first, const const_iterator &last)
	{
		range_insert(first, static_cast<size_type>(first.distance(last)));
	}



	template <class iterator_type>
	void insert (const reverse_iterator &first, const reverse_iterator &last)
	{
		range_insert(first, static_cast<size_type>(first.distance(last)));
	}



	template <class iterator_type>
	void insert (const const_reverse_iterator &first, const const_reverse_iterator &last)
	{
		range_insert(first, static_cast<size_type>(first.distance(last)));
	}



	#ifdef PLF_MOVE_SEMANTICS_SUPPORT
		// Range insert, move_iterator overload:

		template <class iterator_type>
		void insert (const std::move_iterator<iterator_type> first, const std::move_iterator<iterator_type> last)
		{
			range_insert(first, static_cast<size_type>(std::distance(first.base(), last.base())));
		}
	#endif



	#ifdef PLF_INITIALIZER_LIST_SUPPORT
		// Initializer-list insert:

		void insert (const std::initializer_list<element_type> &element_list)
		{
			range_insert(element_list.begin(), static_cast<size_type>(element_list.size()));
		}
	#endif



	#ifdef PLF_CPP20_SUPPORT
		template<plf::compatible_range<element_type> range_type>
		void insert_range(range_type &&the_range)
		{
			range_insert(std::ranges::begin(the_range), static_cast<size_type>(std::ranges::distance(the_range)));
		}
	#endif



private:


	void add_to_groups_with_erasures_list(const group_pointer_type group_to_add) PLF_NOEXCEPT
	{
		group_to_add->erasures_list_next_group = erasure_groups_head;

		if (erasure_groups_head != NULL)
		{
			erasure_groups_head->erasures_list_previous_group = group_to_add;
		}

		erasure_groups_head = group_to_add;
	}



	void remove_from_groups_with_erasures_list(const group_pointer_type group_pointer) PLF_NOEXCEPT
	{
		if (group_pointer != erasure_groups_head)
		{
			group_pointer->erasures_list_previous_group->erasures_list_next_group = group_pointer->erasures_list_next_group;

			if (group_pointer->erasures_list_next_group != NULL)
			{
				group_pointer->erasures_list_next_group->erasures_list_previous_group = group_pointer->erasures_list_previous_group;
			}
		}
		else
		{
			erasure_groups_head = erasure_groups_head->erasures_list_next_group;
		}
	}



	void remove_from_groups_with_erasures_list_if_necessary(const group_pointer_type group_pointer) PLF_NOEXCEPT
	{
		if (!group_pointer->no_erasures()) remove_from_groups_with_erasures_list(group_pointer);
	}



	void reset_only_group_left(const group_pointer_type group_pointer) PLF_NOEXCEPT
	{
		erasure_groups_head = NULL;
		group_pointer->reset(0, NULL, NULL, 0);

		// Reset begin and end iterators:
		end_iterator.element_pointer = begin_iterator.element_pointer = group_pointer->front();
		end_iterator.skipfield_pointer = begin_iterator.skipfield_pointer = group_pointer->skipfield;
	}



	void add_to_unused_groups_list(const group_pointer_type group_pointer) PLF_NOEXCEPT
	{
		group_pointer->next_group = unused_groups_head;
		unused_groups_head = group_pointer;
	}




public:

	iterator erase(const const_iterator &it) // if uninitialized/invalid iterator supplied, function could generate an exception
	{
		assert(total_size != 0);
		assert(it.group_pointer != NULL); // ie. not uninitialized iterator
		assert(it.element_pointer != end_iterator.element_pointer); // ie. != end()
		assert(*(it.skipfield_pointer) == 0); // ie. element pointed to by iterator has not been erased previously

		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value) // Avoid the function call if possible
		#endif
		destroy_element(it.element_pointer);
		--total_size;

		if (--(it.group_pointer->size) != 0) // ie. non-empty group at this point in time, don't consolidate
		{
			// Code logic for following section:
			// ---------------------------------
			// If current skipfield node has no skipblock on either side, create new skipblock of size 1
			// If node only has skipblock on left, set current node and start node of the skipblock to left node value + 1.
			// If node only has skipblock on right, make this node the start node of the skipblock and update end node
			// If node has skipblocks on left and right, set start node of left skipblock and end node of right skipblock to the values of the left + right nodes + 1

			// Optimization explanation:
			// The contextual logic below is the same as that in the insert() functions but in this case the value of the current skipfield node will always be
			// zero (since it is not yet erased), meaning no additional manipulations are necessary for the previous skipfield node comparison - we only have to check against zero
			const skipfield_type prev_skipfield = *(it.skipfield_pointer - (it.skipfield_pointer != it.group_pointer->skipfield)); // true if previous node is erased or this node is at beginning of skipfield
			const skipfield_type after_skipfield = *(it.skipfield_pointer + 1);	// NOTE: boundary test (checking against end-of-elements) is able to be skipped due to the extra skipfield node (compared to element field) - which is present to enable faster iterator operator ++ operations

			if ((prev_skipfield == 0) & (after_skipfield == 0)) // no consecutive erased elements
			{
				*it.skipfield_pointer = 1; // solo skipped node
				const skipfield_type index = static_cast<skipfield_type>(it.element_pointer - it.group_pointer->front());

				if (!it.group_pointer->no_erasures())
				{
					edit_free_list_next(it.group_pointer->front() + it.group_pointer->free_list_head, index); // set prev free list head's 'next index' number to the index of the current element
				}
				else
				{
					add_to_groups_with_erasures_list(it.group_pointer);
				}

				edit_free_list_head(it.element_pointer, it.group_pointer->free_list_head);
				it.group_pointer->free_list_head = index;
			}
			else if ((prev_skipfield != 0) & (after_skipfield == 0)) // previous erased consecutive elements, none following
			{
				*(it.skipfield_pointer - prev_skipfield) = *it.skipfield_pointer = static_cast<skipfield_type>(prev_skipfield + 1);
			}
			else if ((prev_skipfield == 0) & (after_skipfield != 0)) // following erased consecutive elements, none preceding
			{
				*(it.skipfield_pointer + after_skipfield) = *(it.skipfield_pointer) = after_skipfield + 1;

				const skipfield_type following_previous = *(pointer_cast<skipfield_pointer_type>(it.element_pointer + 1));
				const skipfield_type following_next = *(pointer_cast<skipfield_pointer_type>(it.element_pointer + 1) + 1);
				edit_free_list_prev(it.element_pointer, following_previous);
				edit_free_list_next(it.element_pointer, following_next);

				const skipfield_type index = static_cast<skipfield_type>(it.element_pointer - it.group_pointer->front());

				if (following_previous != std::numeric_limits<skipfield_type>::max())
				{
					edit_free_list_next(it.group_pointer->front() + following_previous, index); // Set next index of previous free list node to this node's 'next' index
				}

				if (following_next != std::numeric_limits<skipfield_type>::max())
				{
					edit_free_list_prev(it.group_pointer->front() + following_next, index);	// Set previous index of next free list node to this node's 'previous' index
				}
				else
				{
					it.group_pointer->free_list_head = index;
				}
			}
			else // both preceding and following consecutive erased elements - erased element is between two skipblocks
			{
				*(it.skipfield_pointer) = 1; // This line necessary in order for get_iterator() to work - ensures that erased element skipfield nodes are always non-zero

				// Join the skipblocks
				*(it.skipfield_pointer - prev_skipfield) = *(it.skipfield_pointer + after_skipfield) = static_cast<skipfield_type>(prev_skipfield + after_skipfield + 1);

				// Remove the following skipblock's entry from the free list
				const skipfield_type following_previous = *(pointer_cast<skipfield_pointer_type>(it.element_pointer + 1));
				const skipfield_type following_next = *(pointer_cast<skipfield_pointer_type>(it.element_pointer + 1) + 1);

				if (following_previous != std::numeric_limits<skipfield_type>::max())
				{
					edit_free_list_next(it.group_pointer->front() + following_previous, following_next); // Set next index of previous free list node to this node's 'next' index
				}

				if (following_next != std::numeric_limits<skipfield_type>::max())
				{
					edit_free_list_prev(it.group_pointer->front() + following_next, following_previous); // Set previous index of next free list node to this node's 'previous' index
				}
				else
				{
					it.group_pointer->free_list_head = following_previous;
				}
			}

			iterator return_iterator(it.group_pointer, it.element_pointer + after_skipfield + 1, it.skipfield_pointer + after_skipfield + 1);

			if (return_iterator.element_pointer == it.group_pointer->past_back() && it.group_pointer != end_iterator.group_pointer)
			{
				return_iterator.group_pointer = it.group_pointer->next_group;
				return_iterator.set_to_first_element_in_group();
			}

			if (it.element_pointer == begin_iterator.element_pointer) begin_iterator = return_iterator; // If original iterator was first element in hive, update it's value with the next non-erased element:

			return return_iterator;
		}

		// else: group is empty, consolidate groups
		const bool in_back_block = (it.group_pointer == end_iterator.group_pointer), in_front_block = (it.group_pointer == begin_iterator.group_pointer);

		if (in_back_block & in_front_block) // ie. only group in colony
		{
			// Reset skipfield and free list rather than clearing - leads to fewer allocations/deallocations:
			reset_only_group_left(it.group_pointer);
			return end_iterator;
		}
		else if ((!in_back_block) & in_front_block) // ie. Remove first group, change first group to next group
		{
			begin_iterator.group_pointer = it.group_pointer->next_group; // Make the next group the first group
			begin_iterator.set_to_first_element_in_group();
			begin_iterator.group_pointer->previous_group = NULL; // Cut off this group from the chain
			// note: end iterator only needs to be changed if the deleted group was the final group in the chain ie. not in this case

			remove_from_groups_with_erasures_list(it.group_pointer); // Not a back group, so will always have prior erasures at this point.
			deallocate_group_remove_capacity(it.group_pointer);
			return begin_iterator;
		}
		else if (!(in_back_block | in_front_block)) // this is a non-first group but not final group in chain: delete the group, then link previous group to the next group in the chain:
		{
			it.group_pointer->next_group->previous_group = it.group_pointer->previous_group;
			const group_pointer_type return_group = it.group_pointer->previous_group->next_group = it.group_pointer->next_group; // close the chain, removing this group from it
			remove_from_groups_with_erasures_list(it.group_pointer); // As above, not a back group, hence will always have prior erasures at this point.

			if PLF_CONSTEXPR (priority == performance)
			{
				if (it.group_pointer->next_group != end_iterator.group_pointer)
				{
					deallocate_group_remove_capacity(it.group_pointer);
				}
				else
				{ // ie. second to last block in iterative sequence
					add_to_unused_groups_list(it.group_pointer);
				}
			}
			else
			{
				deallocate_group_remove_capacity(it.group_pointer);
			}

			// Return next group's first non-erased element:
			return iterator(return_group, return_group->front() + *(return_group->skipfield), return_group->skipfield + *(return_group->skipfield));
		}
		else // this is a non-first group and the final group in the chain
		{
			remove_from_groups_with_erasures_list_if_necessary(it.group_pointer); // Will not have prior erasures if only one element was ever inserted into this group and that was the element we erased.
			it.group_pointer->previous_group->next_group = NULL;
			end_iterator.group_pointer = it.group_pointer->previous_group; // end iterator needs to be changed as element supplied was the back element of the colony
			end_iterator.set_to_past_back_of_group();

			if PLF_CONSTEXPR (priority == performance)
			{
				add_to_unused_groups_list(it.group_pointer);
			}
			else
			{
				if (unused_groups_head != NULL) // priority == memory_use, if there are other reserved blocks already, get rid of this one
				{
					deallocate_group_remove_capacity(it.group_pointer);
				}
				else // otherwise, retain - prevents unnecessary allocations/deallocations with stack-like usage
				{
					add_to_unused_groups_list(it.group_pointer);
				}
			}

			return end_iterator;
		}
	}



private:


	void partially_erase_group(const const_iterator &start, const aligned_pointer_type end)
	{
		// For the partial block erasures, we have to remove the existing skipblocks within the range from the intra-block free list of skipblocks. However if there're no erasures in the block, we can avoid doing so.
		const_iterator current = start;
		skipfield_type erasure_count = 0;

		// First erase all elements until end of block & remove all skipblocks post-initial position from the free_list. Then, either update preceding skipblock or create new one:

		if (start.group_pointer->no_erasures())
		{
			erasure_count += static_cast<skipfield_type>(end - start.element_pointer);
			add_to_groups_with_erasures_list(start.group_pointer);

			#ifdef PLF_TYPE_TRAITS_SUPPORT
				if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
			#endif
			{
				do // Avoid checking skipfield as there are no erased elements
				{
					destroy_element(current.element_pointer);
				} while (++current.element_pointer != end);
			}
		}
		else
		{
			do // we already know here that first element is not erased
			{
				#ifdef PLF_TYPE_TRAITS_SUPPORT
					if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
				#endif
				destroy_element(current.element_pointer);

				++erasure_count;
				++current.element_pointer;
				++current.skipfield_pointer;

				if (*current.skipfield_pointer != 0)
				{
					// jump over and remove skipblock
					const skipfield_type prev_free_list_index = *(pointer_cast<skipfield_pointer_type>(current.element_pointer));
					const skipfield_type next_free_list_index = *(pointer_cast<skipfield_pointer_type>(current.element_pointer) + 1);

					current.element_pointer += *(current.skipfield_pointer);
					current.skipfield_pointer += *(current.skipfield_pointer);

					if (next_free_list_index == std::numeric_limits<skipfield_type>::max() && prev_free_list_index == std::numeric_limits<skipfield_type>::max()) // if this is the last skipblock in the free list
					{
						current.group_pointer->free_list_head = std::numeric_limits<skipfield_type>::max();
						erasure_count += static_cast<skipfield_type>(end - current.element_pointer);

						#ifdef PLF_TYPE_TRAITS_SUPPORT
							if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
						#endif
						{
							while (current.element_pointer != end) destroy_element(current.element_pointer++); // Avoids checking skipfield, as there are no erased elements left in block
						}

						break; // end overall while loop
					}
					else if (next_free_list_index == std::numeric_limits<skipfield_type>::max()) // if this is the head of the free list
					{
						current.group_pointer->free_list_head = prev_free_list_index; // make free list head equal to next free list node
						edit_free_list_next(current.group_pointer->front() + prev_free_list_index, std::numeric_limits<skipfield_type>::max());
					}
					else // either a tail or middle free list node
					{
						edit_free_list_prev(current.group_pointer->front() + next_free_list_index, prev_free_list_index);

						if (prev_free_list_index != std::numeric_limits<skipfield_type>::max()) // ie. not the tail free list node
						{
							edit_free_list_next(current.group_pointer->front() + prev_free_list_index, next_free_list_index);
						}
					}
				}
			} while (current.element_pointer != end);
		}


		// Update jump-counting skipfield:
		const size_type distance_to_end = static_cast<skipfield_type>(end - start.element_pointer);
		const skipfield_type start_index = static_cast<skipfield_type>(start.element_pointer - start.group_pointer->front());
		const size_type previous_node_value = (start_index == 0) ? 0 : *(start.skipfield_pointer - 1);

		if (previous_node_value == 0) // start element is either at start of block, or previous element is non-erased so no adjacent skipblock
		{
			*(start.skipfield_pointer) = *(start.skipfield_pointer + distance_to_end - 1) = static_cast<skipfield_type>(distance_to_end); // set start and end node of skipblock

			if (!start.group_pointer->no_erasures())
			{
				edit_free_list_next(start.group_pointer->front() + start.group_pointer->free_list_head, start_index);
			}

			edit_free_list_head(start.element_pointer, start.group_pointer->free_list_head);
			start.group_pointer->free_list_head = start_index;
		}
		else
		{
			// Just update existing skipblock, no need to create new free list node:
			*(start.skipfield_pointer - previous_node_value) = *(start.skipfield_pointer + distance_to_end - 1) = static_cast<skipfield_type>(previous_node_value + distance_to_end);
		}

		if (distance_to_end > 2) // if the skipblock is longer than 2 nodes, fill in the middle nodes with non-zero values so that get_iterator() will work
		{
			std::memset(plf::void_cast(start.skipfield_pointer + 1), 1, sizeof(skipfield_type) * (distance_to_end - 2));
		}

		// Update group and hive size:
		start.group_pointer->size -= erasure_count;
		total_size -= erasure_count;
	}



public:


	// Range erase:

	iterator erase(const const_iterator &iterator1, const const_iterator &iterator2)	// if uninitialized/invalid iterators supplied, function could generate an exception. If iterator1 > iterator2, behaviour is undefined.
	{
		// General code logic: if iterator1 and iterator2 point to elements in the same block, we skip to code section 3 (final block).
		// If they aren't and iterator1 isn't the first non-erased element in first block, we erase part of that block and update accordingly in code Section 1.
		// If it is the first non-erased element, it gets handled in code section 2.
		// In code Section 2 we fully erase and remove all intermediate blocks which aren't the final block. This's optimal as we can just discard the blocks and not do any skipfield or free list updating.
		// In code Section 3 we either partially or fully erase (if iterator2 == end()) the final block in the supplied sequence. If iterator2 was the first non-erased element in it's block or iterator1 == iterator2, this is caught at this point and no action is taken.

		assert(iterator1 <= iterator2);

		const_iterator current = iterator1;

		if (iterator1.group_pointer != iterator2.group_pointer)
		{
			// Section 1: process first block, if partial block erasure
			// ========================================================
			if (current.element_pointer != current.group_pointer->first_element()) // if iterator1 is not the first non-erased element in it's block - most common case
			{
				partially_erase_group(iterator1, iterator1.group_pointer->past_back());
				current.group_pointer = current.group_pointer->next_group;
			}


			// Section 2: remove all intermediate blocks before final block (including first block if it's a full block erasure rather than partial)
			// ====================================================================================================================================
			const group_pointer_type previous_group = current.group_pointer->previous_group;

			while (current.group_pointer != iterator2.group_pointer)
			{
				#ifdef PLF_TYPE_TRAITS_SUPPORT
					if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
				#endif
				{
					current.set_to_first_element_in_group();
					destroy_group(current, current.group_pointer->past_back());
				}

				remove_from_groups_with_erasures_list_if_necessary(current.group_pointer);
				total_size -= current.group_pointer->size;
				const group_pointer_type current_group = current.group_pointer;
				current.group_pointer = current.group_pointer->next_group;

				if (current_group != end_iterator.group_pointer && current_group->next_group != end_iterator.group_pointer)
				{
					deallocate_group_remove_capacity(current_group);
				}
				else
				{
					add_to_unused_groups_list(current_group);
				}
			}

			current.set_to_first_element_in_group();
			current.group_pointer->previous_group = previous_group; // Join this group to the previous non-removed group

			if (previous_group != NULL)
			{
				previous_group->next_group = current.group_pointer;
			}
			else
			{
				begin_iterator = iterator(iterator2.group_pointer, iterator2.element_pointer, iterator2.skipfield_pointer); // This line is included primarily to avoid a secondary if statement within the if block below - it is not needed otherwise
			}
		}


		// Section 3: final block
		// =======================================================================
		// Code logic:
		// If not erasing entire final block, 1. Destruct elements (if non-trivial destructor), 2. add skipblock location to group's free list of skipblocks, and 3. update skipfield.
		// If erasing entire block, 1. Destruct elements (if non-trivial destructor), 2. if no elements left in hive reset the group, otherwise 3. reset end_iterator and remove group from groups-with-erasures list (if prior erasures are present).
		// Note: only way that entire block can be erased is if iterator2 == end() in this case, hence why we reset end_iterator.

		if (current.element_pointer != iterator2.element_pointer) // in case iterator2 was at beginning of it's block - also covers empty range case (first == last)
		{
			if (iterator2.element_pointer != end_iterator.element_pointer || current.element_pointer != current.group_pointer->first_element()) // ie. not erasing entire block. Second condition can only (potentially) be true if iterator1 & iterator2 are in same block.
			{
				partially_erase_group(current, iterator2.element_pointer);
				// Note: there is no need to check for an immediately-adjacent skipblock following the range, as that would imply that iterator2 pointed to an erased element,
				// making it an invalid iterator and unable to be supplied to erase().

				if (iterator1.element_pointer == begin_iterator.element_pointer)
				{
					begin_iterator = iterator(iterator2.group_pointer, iterator2.element_pointer, iterator2.skipfield_pointer);
				}
			}
			else // ie. full block erasure - this can only happen if iterator2 = end() - as it's not valid for an iterator to be pointing one-past the end of a block unless that is the current end().
			{
				#ifdef PLF_TYPE_TRAITS_SUPPORT
					if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
				#endif
				{
					destroy_group(current, iterator2.element_pointer);
				}

				if ((total_size -= current.group_pointer->size) != 0) // ie. hive is not empty
				{
					remove_from_groups_with_erasures_list_if_necessary(current.group_pointer);
					current.group_pointer->previous_group->next_group = current.group_pointer->next_group;

					end_iterator.group_pointer = current.group_pointer->previous_group;
					end_iterator.set_to_past_back_of_group();
					add_to_unused_groups_list(current.group_pointer);
				}
				else // ie. hive is now empty
				{
					// Reset skipfield and free list rather than clearing - leads to fewer allocations/deallocations:
					reset_only_group_left(current.group_pointer);
				}

				return end_iterator;
			}
		}

		return iterator(iterator2.group_pointer, iterator2.element_pointer, iterator2.skipfield_pointer);
	}



private:

	void prepare_groups_for_assign(const size_type size)
	{
		// Destroy all elements if non-trivial:
		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
		#endif
		{
			destroy_remainder(begin_iterator);
		}

		reserve(size);

		// Join all unused_groups to main chain:
		end_iterator.group_pointer->next_group = unused_groups_head;

		begin_iterator.set_to_front_of_group();
		erasure_groups_head = NULL;
		total_size = 0;
	}



public:


	// Fill assign:

	void assign(size_type size, const element_type &element)
	{
		if (size == 0)
		{
			reset();
			return;
		}

		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR ((std::is_trivially_destructible<element_type>::value && std::is_trivially_constructible<element_type>::value && std::is_trivially_copy_assignable<element_type>::value) || !std::is_copy_assignable<element_type>::value) // ie. If there is no benefit nor difference to assigning vs constructing, or if we can't assign, use faster method:
			{
				prepare_groups_for_assign(size);
				fill_unused_groups(size, element, 0, NULL, begin_iterator.group_pointer);
			}
			else
		#endif
		{
			if (total_size == 0)
			{
				prepare_groups_for_assign(size);
				fill_unused_groups(size, element, 0, NULL, begin_iterator.group_pointer);
			}
			else if (size < total_size)
			{
				iterator current = begin_iterator;

				do
				{
					*current++ = element;
				} while (--size != 0);

				erase(current, end_iterator);
			}
			else
			{
				iterator current = begin_iterator;

				do
				{
					*current = element;
				} while (++current != end_iterator);

				insert(size - total_size, element);
			}
		}
	}



private:


	void reset_group_range_assign(iterator &it) PLF_NOEXCEPT
	{
		std::memset(plf::void_cast(it.group_pointer->skipfield), 0, it.group_pointer->capacity * sizeof(skipfield_type));
		it.group_pointer->size = static_cast<skipfield_type>(it.element_pointer - it.group_pointer->front());
	}



	void make_back_group_range_assign(iterator &it) PLF_NOEXCEPT
	{
		// Add all subsequent active groups to unused_groups list:
 		if (it.group_pointer != end_iterator.group_pointer)
		{
			end_iterator.group_pointer->next_group = unused_groups_head;
			unused_groups_head = it.group_pointer->next_group;
		}

		end_iterator = it;
		it.group_pointer->next_group = NULL;
	}



	void finish_range_assign(iterator &it) PLF_NOEXCEPT
	{
		reset_group_range_assign(it);
		make_back_group_range_assign(it);
	}



	void check_iterator_end_of_block(const_iterator &it) PLF_NOEXCEPT
	{
		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
		#endif
		{
			if (it.element_pointer == it.group_pointer->past_back())
			{
				it.group_pointer = it.group_pointer->next_group;
				it.set_to_first_element_in_group();
			}
		}
	}



	// Range assign core:

	template <class iterator_type>
	void range_assign(iterator_type it, size_type size)
	{
		if (size == 0)
		{
			reset();
			return;
		}

		if (total_size == 0)
		{
			prepare_groups_for_assign(size);
			range_fill_unused_groups(size, it, 0, NULL, begin_iterator.group_pointer);
		}
		else
		{
			erasure_groups_head = NULL;
			total_size = 0;
			begin_iterator.set_to_front_of_group();

			for (iterator current(begin_iterator); current.group_pointer != NULL; current.group_pointer = current.group_pointer->next_group)
			{
				current.set_to_front_of_group();
				current.group_pointer->free_list_head = std::numeric_limits<skipfield_type>::max();

				for (const aligned_pointer_type end = (current.group_pointer == end_iterator.group_pointer) ? end_iterator.element_pointer : current.group_pointer->past_back(); current.element_pointer != end;)
				{
					if (*(current.skipfield_pointer) != 0)
					{
						skipfield_type skipblock_length = *(current.skipfield_pointer);
						const_iterator next_element(current.group_pointer, current.element_pointer + skipblock_length, current.skipfield_pointer + skipblock_length);
						skipblock_length = (skipblock_length > size) ? static_cast<skipfield_type>(size) : skipblock_length;

						for (const aligned_pointer_type fill_end = current.element_pointer + skipblock_length; current.element_pointer != fill_end; ++current.element_pointer)
						{
							#ifdef PLF_EXCEPTIONS_SUPPORT
								#ifdef PLF_TYPE_TRAITS_SUPPORT
									if PLF_CONSTEXPR (std::is_nothrow_copy_constructible<element_type>::value)
									{
										PLF_CONSTRUCT_ELEMENT(current.element_pointer, *it++);
									}
									else
								#endif
								{
									try
									{
										PLF_CONSTRUCT_ELEMENT(current.element_pointer, *it++);
									}
									catch (...)
									{
										#ifdef PLF_TYPE_TRAITS_SUPPORT
											if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
										#endif
										{
											check_iterator_end_of_block(next_element);
											destroy_remainder(next_element);
										}

										current.skipfield_pointer = current.group_pointer->skipfield + (current.element_pointer - current.group_pointer->front());
										finish_range_assign(current);
										throw;
									}
								}
							#else
								PLF_CONSTRUCT_ELEMENT(current.element_pointer, *it++);
							#endif

							++total_size;
						}

						current.skipfield_pointer += skipblock_length;

						if ((size -= skipblock_length) == 0)
						{
							#ifdef PLF_TYPE_TRAITS_SUPPORT
								if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
							#endif
							{
								check_iterator_end_of_block(next_element);
								destroy_remainder(next_element);
							}

							finish_range_assign(current);
							return;
						}
					}
					else
					{
						#ifdef PLF_EXCEPTIONS_SUPPORT
							#ifdef PLF_TYPE_TRAITS_SUPPORT
								if PLF_CONSTEXPR (std::is_nothrow_copy_assignable<element_type>::value)
								{
									*pointer_cast<pointer>(current.element_pointer) = *it++;
								}
								else
							#endif
							{
								try
								{
									*pointer_cast<pointer>(current.element_pointer) = *it++;
								}
								catch (...)
								{
									#ifdef PLF_TYPE_TRAITS_SUPPORT
										if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
									#endif
									destroy_remainder(current);

									finish_range_assign(current);
									throw;
								}
							}
						#else
							*pointer_cast<pointer>(current.element_pointer) = *it++;
						#endif

						++total_size;

						if (--size == 0)
						{
							#ifdef PLF_TYPE_TRAITS_SUPPORT
								if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
							#endif
							destroy_remainder(++iterator(current)); // ++ to potentially allow for skipping over a skipblock
							++current.element_pointer; // As opposed to just incrementing, as we do here
							++current.skipfield_pointer;
							finish_range_assign(current);
							return;
						}

						++current.element_pointer;
						++current.skipfield_pointer;
					}
				}

				reset_group_range_assign(current);
			}

			// Use up any remaining space at end of end block (would not be correctly identified above because the skipfield in unused nodes is 0)
			for (const aligned_pointer_type end = end_iterator.group_pointer->past_back(); end_iterator.element_pointer != end;)
			{
				PLF_CONSTRUCT_ELEMENT(end_iterator.element_pointer, *it++);
				++total_size;
				++end_iterator.group_pointer->size;
				++end_iterator.element_pointer;
				++end_iterator.skipfield_pointer;

				if (--size == 0) return;
			}

			// Indicates we've reached the end of existing groups with elements, now can only reused unused groups or create new ones:
			range_insert(it, size);
		}
	}



public:

	// Range assign:

	template <class iterator_type>
	void assign(const typename plf::enable_if<!std::numeric_limits<iterator_type>::is_integer, iterator_type>::type &first, const iterator_type &last)
	{
		range_assign(first, static_cast<size_type>(std::distance(first, last)));
	}



	// Overloads for colony iterators, since std::distance overloads for these are not possible without C++20 concepts and we must stick with ADL:
	template <class iterator_type, bool is_const>
	void assign(const colony_iterator<is_const> &first, const colony_iterator<is_const> &last)
	{
		range_assign(first, static_cast<size_type>(first.distance(last)));
	}



	template <class iterator_type, bool is_const>
	void assign(const colony_reverse_iterator<is_const> &first, const colony_reverse_iterator<is_const> &last)
	{
		range_assign(first, static_cast<size_type>(first.distance(last)));
	}



	#ifdef PLF_MOVE_SEMANTICS_SUPPORT
		// Range assign, move_iterator overload:

		template <class iterator_type>
		void assign (const std::move_iterator<iterator_type> &first, const std::move_iterator<iterator_type> &last)
		{
			range_assign(first, static_cast<size_type>(std::distance(first.base(), last.base())));
		}
	#endif



	#ifdef PLF_INITIALIZER_LIST_SUPPORT
		// Initializer-list assign:

		void assign(const std::initializer_list<element_type> &element_list)
		{
			range_assign(element_list.begin(), static_cast<size_type>(element_list.size()));
		}
	#endif



	#ifdef PLF_CPP20_SUPPORT
		template<plf::compatible_range<element_type> range_type>
		void assign_range(range_type &&the_range)
		{
			range_assign(std::ranges::begin(the_range), static_cast<size_type>(std::ranges::distance(the_range)));
		}



		[[nodiscard]]
	#endif
	bool empty() const PLF_NOEXCEPT
	{
		return total_size == 0;
	}



	size_type size() const PLF_NOEXCEPT
	{
		return total_size;
	}



	PLF_CONSTFUNC size_type max_size() const PLF_NOEXCEPT
	{
		#ifdef PLF_ALLOCATOR_TRAITS_SUPPORT
			return std::allocator_traits<allocator_type>::max_size(*this);
		#else
			return allocator_type::max_size();
		#endif
	}



	size_type capacity() const PLF_NOEXCEPT
	{
		return total_capacity;
	}



	size_type memory() const PLF_NOEXCEPT
	{
		size_type memory_use = sizeof(*this); // sizeof colony structure
		end_iterator.group_pointer->next_group = unused_groups_head; // temporarily link the active & reserved (unused) groups in order to only have one loop below instead of two

		for(group_pointer_type current = begin_iterator.group_pointer; current != NULL; current = current->next_group)
		{
			memory_use += sizeof(group) + (get_aligned_block_capacity(current->capacity) * sizeof(aligned_allocation_struct)); // add element/skipfield memory block sizes + size of the group struct
		}

		end_iterator.group_pointer->next_group = NULL; // unlink active & reserved groups
		return memory_use;
	}



	static PLF_CONSTFUNC size_type block_metadata_memory(const size_type block_capacity) PLF_NOEXCEPT
	{
		return sizeof(group) + ((get_aligned_block_capacity(static_cast<skipfield_type>(block_capacity)) - block_capacity) * sizeof(aligned_allocation_struct));
	}



	static PLF_CONSTFUNC size_type block_allocation_amount(size_type block_capacity) PLF_NOEXCEPT
	{
		if (block_capacity > std::numeric_limits<skipfield_type>::max()) block_capacity = std::numeric_limits<skipfield_type>::max();

		return sizeof(aligned_allocation_struct) * get_aligned_block_capacity(static_cast<skipfield_type>(block_capacity));
	}



	static PLF_CONSTFUNC size_type max_elements_per_allocation(const size_type allocation_amount) PLF_NOEXCEPT
	{
		// Get a rough approximation of the number of elements + skipfield units we can fit in the amount expressed:
		PLF_CONSTFUNC size_type num_units = allocation_amount / (sizeof(aligned_element_struct) + sizeof(skipfield_type));
		PLF_CONSTFUNC plf::limits hard_capacities = block_capacity_hard_limits();

		// Truncate the amount to the implementation's hard block capacity max limit:
		if (num_units > hard_capacities.max) num_units = hard_capacities.max;

		// Adjust num_units downward based on (a) the additional skipfield node necessary per-block in this implementation and
		// (b) any additional memory waste required in order to allocate the skipfield in multiples of the element type's alignof:
		if ((	/* Explanation: elements and skipfield are allocated in a single allocation to save performance.
				In order for the elements to be correctly aligned in memory, this single allocation is aligned to the alignof
				the element type, so the first line below is the allocation amount in bytes required for the skipfield
				when allocated in multiples of the element type's alignof. The + sizeof(skipfield_type) adds the additional skipfield node
				as mentioned, and the (num_units + 1) minus 1 byte rounds up the integer division: */
			(((((num_units + 1) * sizeof(aligned_allocation_struct)) - 1) + sizeof(skipfield_type)) / sizeof(aligned_allocation_struct))
				/* the second line is the amount of memory in bytes necessary for the elements themselves: */
			+ (num_units * sizeof(aligned_element_struct)))
				/* then we compare against the desired allocation amount: */
			> allocation_amount)
		{
			--num_units; // In this implementation it is not possible for the necessary adjustment to be greater than 1 element+skipfield sizeof
		}

		if (num_units < hard_capacities.min) num_units = 0;

		return num_units;
	}



	#ifdef PLF_COLONY_TEST_DEBUG // used for debugging during internal testing only:
		size_type group_size_sum() const PLF_NOEXCEPT
		{
			size_type temp = 0;

			for (group_pointer_type current = begin_iterator.group_pointer; current != NULL; current = current->next_group)
			{
				temp += current->size;
			}

			return temp;
		}
	#endif



private:

	// get all elements contiguous in memory and shrink to fit, remove erasures and free lists. Invalidates all iterators and pointers to elements.
	void consolidate(const skipfield_type new_min, const skipfield_type new_max)
	{
		colony temp(plf::limits(new_min, new_max));
		temp.reserve(total_size);
		temp.end_iterator.group_pointer->next_group = temp.unused_groups_head;

		#if defined(PLF_MOVE_SEMANTICS_SUPPORT) && defined(PLF_TYPE_TRAITS_SUPPORT)
			if PLF_CONSTEXPR (!std::is_trivially_copyable<element_type>::value && std::is_nothrow_move_constructible<element_type>::value)
			{
				temp.range_fill_unused_groups(total_size, plf::make_move_iterator(begin_iterator), 0, NULL, temp.begin_iterator.group_pointer);
			}
			else
		#endif
		{
			temp.range_fill_unused_groups(total_size, begin_iterator, 0, NULL, temp.begin_iterator.group_pointer);
		}

		#ifdef PLF_MOVE_SEMANTICS_SUPPORT
			*this = std::move(temp); // Avoid generating 2nd temporary
		#else
			swap(temp);
		#endif
	}



public:


	void reshape(const plf::limits block_limits)
	{
		check_capacities_conformance(block_limits);
		const skipfield_type new_min = static_cast<skipfield_type>(block_limits.min), new_max = static_cast<skipfield_type>(block_limits.max);

		if (total_capacity != 0)
		{
			if (total_size != 0)
			{
				if (min_block_capacity > new_max || max_block_capacity < new_min) // If none of the original blocks could potentially fit within the new limits, skip checking of blocks and just consolidate:
				{
					consolidate(new_min, new_max);
					return;
				}

				if (min_block_capacity < new_min || max_block_capacity > new_max) // ie. If existing blocks could be outside of the new limits
				{
					// Otherwise need to check all group sizes here (not just back one, which is most likely largest), because splice might append smaller blocks after a larger block:
					for (group_pointer_type current_group = begin_iterator.group_pointer; current_group != NULL; current_group = current_group->next_group)
					{
						if (current_group->capacity < new_min || current_group->capacity > new_max)
						{
							consolidate(new_min, new_max);
							return;
						}
					}
				}
			}
			else // include first group to be checked in the loop below
			{
				begin_iterator.group_pointer->next_group = unused_groups_head;
				unused_groups_head = begin_iterator.group_pointer;
			}

			// If a consolidation or throw has not occured, process reserved/unused groups and deallocate where they don't fit the new limits:

			for (group_pointer_type current_group = unused_groups_head, previous_group = NULL; current_group != NULL;)
			{
				const group_pointer_type next_group = current_group->next_group;

				if (current_group->capacity < new_min || current_group->capacity > new_max)
				{
					deallocate_group_remove_capacity(current_group);

					if (previous_group == NULL)
					{
						unused_groups_head = next_group;
					}
					else
					{
						previous_group->next_group = next_group;
					}
				}
				else
				{
					previous_group = current_group;
				}

				current_group = next_group;
			}

			if (total_size == 0)
			{
				if (unused_groups_head == NULL)
				{
					blank();
				}
				else
				{
					begin_iterator.group_pointer = unused_groups_head;
					unused_groups_head = begin_iterator.group_pointer->next_group;
					begin_iterator.group_pointer->next_group = NULL;
				}
			}
		}

		min_block_capacity = new_min;
		max_block_capacity = new_max;
	}



	PLF_CONSTFUNC plf::limits block_capacity_limits() const PLF_NOEXCEPT
	{
		return plf::limits(static_cast<size_t>(min_block_capacity), static_cast<size_t>(max_block_capacity));
	}



	static PLF_CONSTFUNC plf::limits block_capacity_hard_limits() PLF_NOEXCEPT
	{
		return plf::limits(3, std::min(static_cast<size_t>(std::numeric_limits<skipfield_type>::max()), max_size_static()));
	}



	void clear() PLF_NOEXCEPT
	{
		if (total_size == 0) return;

		// Destroy all elements if element type is non-trivial:
		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR (!std::is_trivially_destructible<element_type>::value)
		#endif
		{
			destroy_remainder(begin_iterator);
		}

		if (begin_iterator.group_pointer != end_iterator.group_pointer)
		{ // Move all other groups onto the unused_groups list
			end_iterator.group_pointer->next_group = unused_groups_head;
			unused_groups_head = begin_iterator.group_pointer->next_group;
			end_iterator.group_pointer = begin_iterator.group_pointer; // other parts of iterator reset in the function below
		}

		reset_only_group_left(begin_iterator.group_pointer);
		erasure_groups_head = NULL;
		total_size = 0;
	}



	colony & operator = (const colony &source)
	{
		assert(&source != this);

		#ifdef PLF_ALLOCATOR_TRAITS_SUPPORT
			if PLF_CONSTEXPR (std::allocator_traits<allocator_type>::propagate_on_container_copy_assignment::value)
		#endif
		{
			#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT
				if PLF_CONSTEXPR (!std::allocator_traits<allocator_type>::is_always_equal::value)
			#endif
			{
				if (static_cast<allocator_type &>(*this) != static_cast<const allocator_type &>(source))
				{ // Deallocate existing blocks as source allocator is not necessarily able to do so
					reset();
				}
			}

			static_cast<allocator_type &>(*this) = static_cast<const allocator_type &>(source);
			// Reconstruct rebinds:
			group_allocator = group_allocator_type(*this);
			aligned_struct_allocator = aligned_struct_allocator_type(*this);
			skipfield_allocator = skipfield_allocator_type(*this);
			tuple_allocator = tuple_allocator_type(*this);
		}

		range_assign(source.begin_iterator, source.total_size);
		return *this;
	}



	#ifdef PLF_MOVE_SEMANTICS_SUPPORT
	private:

		void move_assign(colony &&source) PLF_NOEXCEPT
		{
			#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT
				if PLF_CONSTEXPR ((std::is_trivially_copyable<allocator_type>::value || std::allocator_traits<allocator_type>::is_always_equal::value) && std::is_trivially_copyable<group_pointer_type>::value)
				{
					std::memcpy(static_cast<void *>(this), static_cast<void *>(&source), sizeof(colony));
				}
				else
			#endif
			{
				end_iterator = std::move(source.end_iterator);
				begin_iterator = std::move(source.begin_iterator);
				erasure_groups_head = std::move(source.erasure_groups_head);
				unused_groups_head =  std::move(source.unused_groups_head);
				total_size = source.total_size;
				total_capacity = source.total_capacity;
				min_block_capacity = source.min_block_capacity;
				max_block_capacity = source.max_block_capacity;

				#ifdef PLF_ALLOCATOR_TRAITS_SUPPORT
					if PLF_CONSTEXPR(std::allocator_traits<allocator_type>::propagate_on_container_move_assignment::value)
				#endif
				{
					static_cast<allocator_type &>(*this) = static_cast<allocator_type &>(source);
					// Reconstruct rebinds:
					group_allocator = group_allocator_type(*this);
					aligned_struct_allocator = aligned_struct_allocator_type(*this);
					skipfield_allocator = skipfield_allocator_type(*this);
					tuple_allocator = tuple_allocator_type(*this);
				}
			}
		}



	public:

		// Move assignment
		#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT
			colony & operator = (colony &&source) PLF_NOEXCEPT(std::allocator_traits<allocator_type>::propagate_on_container_move_assignment::value || std::allocator_traits<allocator_type>::is_always_equal::value)
		#else
			colony & operator = (colony &&source)
		#endif
		{
			assert(&source != this);
			destroy_all_data();

			#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT
				if PLF_CONSTEXPR (std::allocator_traits<allocator_type>::propagate_on_container_move_assignment::value || std::allocator_traits<allocator_type>::is_always_equal::value)
				{ // Note: we need this to be constexpr to avoid warning errors on the potentially-throwing section below
					move_assign(std::move(source));
				}
				else
			#endif
			if (static_cast<allocator_type &>(*this) == static_cast<allocator_type &>(source))
			{
				move_assign(std::move(source));
			}
			else // Allocator isn't propagatable so move elements from source and deallocate the source's blocks. Could throw here:
			{
				#ifdef PLF_TYPE_TRAITS_SUPPORT
					if PLF_CONSTEXPR (!(std::is_move_constructible<element_type>::value && std::is_move_assignable<element_type>::value))
					{
						range_assign(source.begin_iterator, source.total_size);
					}
					else
				#endif
				{
					range_assign(plf::make_move_iterator(source.begin_iterator), source.total_size);
				}

				source.destroy_all_data();
			}

			source.blank();
			return *this;
		}
	#endif



	#ifdef PLF_INITIALIZER_LIST_SUPPORT
		colony & operator = (const std::initializer_list<element_type> &element_list)
		{
			range_assign(element_list.begin(), static_cast<size_type>(element_list.size()));
			return *this;
		}
	#endif



	friend bool operator == (const colony &lh, const colony &rh) PLF_NOEXCEPT
	{
		if (lh.total_size != rh.total_size) return false;

		for (const_iterator lh_iterator = lh.begin_iterator, rh_iterator = rh.begin_iterator; lh_iterator != lh.end_iterator; ++lh_iterator, ++rh_iterator)
		{
			if (*lh_iterator != *rh_iterator) return false;
		}

		return true;
	}



	friend bool operator != (const colony &lh, const colony &rh) PLF_NOEXCEPT
	{
		return !(lh == rh);
	}



	#ifdef PLF_CPP20_SUPPORT
		friend constexpr std::strong_ordering operator <=> (const colony &lh, const colony &rh)
		{
			return std::lexicographical_compare_three_way(lh.begin(), lh.end(), rh.begin(), rh.end());
		}
	#endif



	void shrink_to_fit()
	{
		if (total_size == total_capacity)
		{
			return;
		}
		else if (total_size == 0)
		{
			reset();
			return;
		}

		consolidate(min_block_capacity, max_block_capacity);
	}



	void trim_capacity() PLF_NOEXCEPT
	{
		if (end_iterator.element_pointer == NULL) return; // empty colony

		while(unused_groups_head != NULL)
		{
			const group_pointer_type next_group = unused_groups_head->next_group;
			deallocate_group_remove_capacity(unused_groups_head);
			unused_groups_head = next_group;
		}

		if (begin_iterator.element_pointer == end_iterator.element_pointer) // ie. clear() has been called prior
		{
			deallocate_group(begin_iterator.group_pointer);
			blank();
		}
	}



	void trim_capacity(const size_type capacity_retain) PLF_NOEXCEPT
	{
		const size_type capacity_difference = total_capacity - capacity_retain;

		if (end_iterator.element_pointer == NULL || total_capacity <= capacity_retain || total_size >= capacity_retain || capacity_difference < min_block_capacity) return;

		size_type capacity_to_remove_remaining = capacity_difference;

		for (group_pointer_type current_group = unused_groups_head, previous_group = NULL; current_group != NULL;)
		{
			const group_pointer_type next_group = current_group->next_group;

			if (capacity_to_remove_remaining >= current_group->capacity)
			{
				capacity_to_remove_remaining -= current_group->capacity;
				deallocate_group(current_group);

				if (previous_group == NULL)
				{
					unused_groups_head = next_group;
				}
				else
				{
					previous_group->next_group = next_group;
				}

				if (capacity_to_remove_remaining < min_block_capacity) break;
			}
			else
			{
				previous_group = current_group;
			}

			current_group = next_group;
		}


		if (begin_iterator.element_pointer == end_iterator.element_pointer) // ie. clear() has been called prior
		{
			if (capacity_to_remove_remaining >= begin_iterator.group_pointer->capacity)
			{
				capacity_to_remove_remaining -= begin_iterator.group_pointer->capacity;
				deallocate_group(begin_iterator.group_pointer);

				if (unused_groups_head != NULL) // some of the reserved blocks were not removed as they were too large, so use one of these to make the new begin group
				{
					end_iterator.group_pointer = begin_iterator.group_pointer = unused_groups_head;
					end_iterator.element_pointer = begin_iterator.element_pointer = unused_groups_head->front();
					end_iterator.skipfield_pointer = begin_iterator.skipfield_pointer = unused_groups_head->skipfield;

					unused_groups_head = unused_groups_head->next_group;
					begin_iterator.group_pointer->next_group = NULL;
				}
				else
				{
					blank();
					return;
				}
			}
		}

		total_capacity -= capacity_difference - capacity_to_remove_remaining;
	}



	void reserve(size_type new_capacity)
	{
		if (new_capacity == 0 || new_capacity <= total_capacity) return; // ie. We already have enough space allocated

		if (new_capacity > max_size())
		{
			#ifdef PLF_EXCEPTIONS_SUPPORT
				throw std::length_error("Capacity requested via reserve() greater than max_size()");
			#else
				std::terminate();
			#endif
		}

		new_capacity -= total_capacity;

		size_type number_of_max_groups = new_capacity / max_block_capacity;
		skipfield_type remainder = static_cast<skipfield_type>(new_capacity - (number_of_max_groups * max_block_capacity)), negative_remainder = 0;
		group_pointer_type deallocatable_group = NULL;

		if (remainder == 0)
		{
			remainder = max_block_capacity;
			--number_of_max_groups;
		}
		else
		{
			// Here we try to increase iteration performance by deallocating a small unused group and allocating one larger group.
			// This also means that if remainder < min_block_capacity we don't have to allocate a min capacity group and spread the difference over subsequent groups (see subsequent if block).
			// The smaller group is not deallocated immediately so that, in the event that an exception is triggered when allocating the larger group, we don't end up with lower capacity than before reserve().

			if (unused_groups_head != NULL && max_block_capacity - remainder >= min_block_capacity)
			{
				deallocatable_group = unused_groups_head;
				group_pointer_type prev_unused_group = NULL;

				do
				{
					const skipfield_type current_capacity = deallocatable_group->capacity;

					// If there exists an unused group which's of low-enough capacity, deallocate that later and add it's capacity to the remainder group:
					if (std::numeric_limits<skipfield_type>::max() - current_capacity > remainder && /* <- to make sure we don't overflow in next line */
						max_block_capacity >= current_capacity + remainder)
					{
						remainder += current_capacity;
						const group_pointer_type next_group = deallocatable_group->next_group;

						if (prev_unused_group != NULL)
						{
							prev_unused_group->next_group = next_group;
						}
						else
						{
						  	unused_groups_head = next_group;
						}

						break;
					}

					prev_unused_group = deallocatable_group;
					deallocatable_group = deallocatable_group->next_group;
				} while (deallocatable_group != NULL);
			}


			if (remainder < min_block_capacity) // Implies we were unable to consolidate remainder with an existing unused group, in the if-block above
			{
				// Note: negative_remainder is used to take the difference between the minimum block capacity limit and the actual remainder, and spread this negative difference over subsequent blocks which are in the usual case at max capacity.
				negative_remainder = min_block_capacity - remainder;
				remainder = min_block_capacity;

	  			// This line checks to see - if we have to reduce the size of the max-capacity blocks to spread the negative_remainder out - whether even reducing the max blocks to min capacity will be enough to keep the capacity under max_size(). We add 1 for the initial (remainder) block. This guards against situations where, for example, the min/max limits are very similar so spreading the negative remainder out is less doable:
				if (max_size() - total_capacity < ((number_of_max_groups + 1) * min_block_capacity))
				{
					#ifdef PLF_EXCEPTIONS_SUPPORT
						throw std::length_error("Reserve cannot increase capacity to >= n without being > max_size() due to current capacity() and block capacity limits");
					#else
						std::terminate();
					#endif
				}
			}
		}


		group_pointer_type current_group, first_unused_group;

		if (begin_iterator.group_pointer == NULL) // Most common scenario - empty hive
		{
			initialize(remainder);
			begin_iterator.group_pointer->size = 0; // Note: this is set to 1 by default in the initialize function (which is optimised for insert())

			if (number_of_max_groups == 0) return;

			// Make the first allocated unused group:
			const skipfield_type new_block_capacity = (max_block_capacity - negative_remainder < min_block_capacity) ? min_block_capacity : max_block_capacity - negative_remainder;
			negative_remainder -= max_block_capacity - new_block_capacity;
			first_unused_group = current_group = allocate_new_group(new_block_capacity, begin_iterator.group_pointer);
			--number_of_max_groups;
		}
		else // Non-empty hive, add first new unused group:
		{
			#ifdef PLF_EXCEPTIONS_SUPPORT
				try
				{
					first_unused_group = current_group = allocate_new_group(remainder, end_iterator.group_pointer);
				}
				catch (...)
				{
					if (deallocatable_group != NULL) // roll back group removal
					{
						add_to_unused_groups_list(deallocatable_group);
					}

					throw;
				}
			#else
				first_unused_group = current_group = allocate_new_group(remainder, end_iterator.group_pointer);
			#endif

			// We've now successfully allocated another group which is guaranteed to be larger than this group, so capacity is larger than it was before reserve() was called even if the other allocations below trigger an exception, and we can deallocate the group:
			if (deallocatable_group != NULL) deallocate_group_remove_capacity(deallocatable_group);
		}


		while (number_of_max_groups != 0)
		{
			const skipfield_type new_block_capacity = (max_block_capacity - negative_remainder < min_block_capacity) ? min_block_capacity : max_block_capacity - negative_remainder;
			negative_remainder -= max_block_capacity - new_block_capacity;

			#ifdef PLF_EXCEPTIONS_SUPPORT
				try
				{
					current_group->next_group = allocate_new_group(new_block_capacity, current_group);
				}
				catch (...)
				{
					current_group->next_group = unused_groups_head;
					unused_groups_head = first_unused_group;
					throw;
				}
			#else
				current_group->next_group = allocate_new_group(new_block_capacity, current_group);
			#endif

			current_group = current_group->next_group;
			--number_of_max_groups;
		}

		current_group->next_group = unused_groups_head;
		unused_groups_head = first_unused_group;
	}



private:

	template <bool is_const>
	colony_iterator<is_const> get_it(const pointer element_pointer) const PLF_NOEXCEPT
	{
		if (end_iterator.group_pointer != NULL)
		{
			const aligned_pointer_type aligned_element_pointer = pointer_cast<aligned_pointer_type>(element_pointer);
			// Note: we start with checking the back group first, as it will be the largest group in most cases, so there's a statistically-higher chance of the element being within it.

			// Special case for back group in case the element was in a group which became empty and got moved to the unused_groups list or was deallocated, and then that memory was re-used (ie. it became the current back group). The following prevents the function from mistakenly returning an iterator which is beyond the back element of the colony:
			if (aligned_element_pointer >= end_iterator.group_pointer->front() && aligned_element_pointer < end_iterator.element_pointer)
			{
				const skipfield_pointer_type skipfield_pointer = end_iterator.group_pointer->skipfield + (aligned_element_pointer - end_iterator.group_pointer->front());
				// The first test below checks to see whether the element is live or erased. The second test checks to see if the pointer points to an erased element in a block of memory which got deallocated by the container, then another block got allocated by the container which contained that memory space, but the new block is not exactly aligned with the old block. eg. colony<structure> where structure is struct {int x, y} and the new memory block is offset by 1 int such that 'element_pointer' points at y.
				return (*skipfield_pointer == 0 && ((reinterpret_cast<char *>(element_pointer) - reinterpret_cast<char *>(end_iterator.group_pointer->elements)) % sizeof(aligned_element_struct) == 0)) ?
					colony_iterator<is_const>(end_iterator.group_pointer, aligned_element_pointer, skipfield_pointer) : colony_iterator<is_const>(end_iterator);
			}

			// All other groups, if any exist:
			for (group_pointer_type current_group = end_iterator.group_pointer->previous_group; current_group != NULL; current_group = current_group->previous_group)
			{
				if (aligned_element_pointer >= current_group->front() && aligned_element_pointer < current_group->past_back())
				{
					const skipfield_pointer_type skipfield_pointer = current_group->skipfield + (aligned_element_pointer - current_group->front());
					return (*skipfield_pointer == 0 && ((reinterpret_cast<char *>(element_pointer) - reinterpret_cast<char *>(current_group->elements)) % sizeof(aligned_element_struct) == 0)) ?
						colony_iterator<is_const>(current_group, aligned_element_pointer, skipfield_pointer) : colony_iterator<is_const>(end_iterator);
				}
			}
		}

		return colony_iterator<is_const>(end_iterator);
	}



public:

	iterator get_iterator(const pointer element_pointer) PLF_NOEXCEPT
	{
		return get_it<false>(element_pointer);
	}



	const_iterator get_iterator(const const_pointer element_pointer) const PLF_NOEXCEPT
	{
		return get_it<true>(const_cast<pointer>(element_pointer));
	}



	bool is_active(const const_iterator &it) const PLF_NOEXCEPT
	{
		if (end_iterator.group_pointer != NULL)
		{
			// Schema: check (a) that the group the iterator belongs to is still active and not deallocated or in the unused_groups list, then (b) that the element is not erased. (a) prevents an out-of-bounds memory access if the group is deallocated. Same reasoning as get_iterator for loop conditions

			// Special case for back group, same reasoning as in get_it():
			if (it.group_pointer == end_iterator.group_pointer && it.element_pointer >= end_iterator.group_pointer->front() && it.element_pointer < end_iterator.element_pointer)
			{
				return (*it.skipfield_pointer == 0);
			}

			for (group_pointer_type current_group = end_iterator.group_pointer->previous_group; current_group != NULL; current_group = current_group->previous_group)
			{
				if (it.group_pointer == current_group && it.element_pointer >= current_group->front() && it.element_pointer < current_group->past_back()) // 2nd 2 conditions necessary in case the group contained the element which the iterator points to, has been deallocated from the colony previously, but then the same pointer address is re-supplied via an allocator for a subsequent group allocation (in which case the group's element block memory location may be different)
				{
					return (*it.skipfield_pointer == 0);
				}
			}
		}

		return false;
	}



	allocator_type get_allocator() const PLF_NOEXCEPT
	{
		return static_cast<allocator_type>(*this);
	}



private:

	void source_blocks_incompatible()
	{
		#ifdef PLF_EXCEPTIONS_SUPPORT
			throw std::length_error("A source memory block capacity is outside of the destination's minimum or maximum memory block capacity limits - please change either the source or the destination's min/max block capacity limits using reshape() before calling splice() in this case");
		#else
			std::terminate();
		#endif
	}



public:

	void splice(colony &source)
	{
		// Process: if there are unused memory spaces at the end of the current back group of the chain, convert them
		// to skipped elements and add the locations to the group's free list.
		// Then link the destination's groups to the source's groups and nullify the source.
		// If the source has more unused memory spaces in the back group than the destination, swap them before processing to reduce the number of locations added to a free list and also subsequent jumps during iteration.

		assert(&source != this);

		if (source.total_size == 0) return;

		// Throw if incompatible block capacities found in source:
		if (source.min_block_capacity > max_block_capacity || source.max_block_capacity < min_block_capacity) // ie. source blocks cannot possibly fit within *this's block capacity limits
		{
			source_blocks_incompatible();
		}
		else if (source.min_block_capacity < min_block_capacity || source.max_block_capacity > max_block_capacity) // ie. source blocks may or may not fit
		{
			for (group_pointer_type current_group = source.begin_iterator.group_pointer; current_group != NULL; current_group = current_group->next_group)
			{
				if (current_group->capacity < min_block_capacity || current_group->capacity > max_block_capacity) source_blocks_incompatible();
			}
		}


		if (total_size != 0)
		{
			// If there's more unused element locations in back memory block of destination than in back memory block of source, swap with source to reduce number of skipped elements during iteration:
			if ((end_iterator.group_pointer->past_back() - end_iterator.element_pointer) > (source.end_iterator.group_pointer->past_back() - source.end_iterator.element_pointer))
			{
				swap(source);
				// Swap back unused groups list and block capacity limits so that source and *this retain their original ones:
				std::swap(source.unused_groups_head, unused_groups_head);
				std::swap(source.min_block_capacity, min_block_capacity);
				std::swap(source.max_block_capacity, max_block_capacity);
			}


			// Add source list of groups-with-erasures to destination list of groups-with-erasures:
			if (source.erasure_groups_head != NULL)
			{
				if (erasure_groups_head != NULL)
				{
					group_pointer_type tail_group = erasure_groups_head;

					while (tail_group->erasures_list_next_group != NULL)
					{
						tail_group = tail_group->erasures_list_next_group;
					}

					tail_group->erasures_list_next_group = source.erasure_groups_head;
					source.erasure_groups_head->erasures_list_previous_group = tail_group;
				}
				else
				{
					erasure_groups_head = source.erasure_groups_head;
				}
			}


			const skipfield_type distance_to_end = static_cast<skipfield_type>(end_iterator.group_pointer->past_back() - end_iterator.element_pointer);

			if (distance_to_end != 0) // 0 == edge case
			{	 // Mark unused element memory locations from back group as skipped/erased:
				// Update skipfield:
				const skipfield_type previous_node_value = *(end_iterator.skipfield_pointer - 1);

				if (previous_node_value == 0) // no previous skipblock
				{
					*end_iterator.skipfield_pointer = distance_to_end;
					*(end_iterator.skipfield_pointer + distance_to_end - 1) = distance_to_end;

					if (distance_to_end > 2) // make erased middle nodes non-zero for get_iterator and is_active
					{
						std::memset(plf::void_cast(end_iterator.skipfield_pointer + 1), 1, sizeof(skipfield_type) * (distance_to_end - 2));
					}

					const skipfield_type index = static_cast<skipfield_type>(end_iterator.element_pointer - end_iterator.group_pointer->front());

					if (!end_iterator.group_pointer->no_erasures())
					{
						edit_free_list_next(end_iterator.group_pointer->front() + end_iterator.group_pointer->free_list_head, index); // set prev free list head's 'next index' number to the index of the current element
					}
					else
					{
						add_to_groups_with_erasures_list(end_iterator.group_pointer);
					}

					edit_free_list_head(end_iterator.element_pointer, end_iterator.group_pointer->free_list_head);
					end_iterator.group_pointer->free_list_head = index;
				}
				else
				{ // update previous skipblock, no need to update free list:
					*(end_iterator.skipfield_pointer - previous_node_value) = *(end_iterator.skipfield_pointer + distance_to_end - 1) = static_cast<skipfield_type>(previous_node_value + distance_to_end);

					if (distance_to_end > 1) // make erased middle nodes non-zero for get_iterator and is_active
					{
						std::memset(plf::void_cast(end_iterator.skipfield_pointer), 1, sizeof(skipfield_type) * (distance_to_end - 1));
					}
				}
			}


			// Join the destination and source group chains:
			end_iterator.group_pointer->next_group = source.begin_iterator.group_pointer;
			source.begin_iterator.group_pointer->previous_group = end_iterator.group_pointer;

			// Update group numbers if necessary:
			if (source.begin_iterator.group_pointer->group_number <= end_iterator.group_pointer->group_number)
			{
				size_type source_group_count = 0;

				for (group_pointer_type current_group = source.begin_iterator.group_pointer; current_group != NULL; current_group = current_group->next_group, ++source_group_count) {}

				if ((std::numeric_limits<size_type>::max() - end_iterator.group_pointer->group_number) >= source_group_count)
				{
					update_subsequent_group_numbers(end_iterator.group_pointer->group_number + 1u, source.begin_iterator.group_pointer);
				}
				else
				#ifdef PLF_CPP20_SUPPORT
					[[unlikely]]
				#endif
				{
					reset_group_numbers();
				}
			}

			end_iterator = source.end_iterator;
			total_size += source.total_size;
			total_capacity += source.total_capacity;
		}
		else // If *this is empty():
		{
			// Preserve unused_groups_head and de-link so that destroy_all_data doesn't remove them:
			const group_pointer_type original_unused_groups = unused_groups_head;
			unused_groups_head = NULL;
			destroy_all_data();
			unused_groups_head = original_unused_groups;

			// Move source data to *this:
			end_iterator = source.end_iterator;
			begin_iterator = source.begin_iterator;
			erasure_groups_head = source.erasure_groups_head;
			total_size = source.total_size;
			total_capacity = source.total_capacity;

			// Add capacity for unused groups back into *this:
			for (group_pointer_type current = original_unused_groups; current != NULL; current = current->next_group)
			{
				total_capacity += current->capacity;
			}
		}


		// Reset source values:
		const group_pointer_type original_unused_groups_head = source.unused_groups_head; // grab value before it gets wiped
		source.blank(); // blank source before adding capacity from unused groups back in

		if (original_unused_groups_head != NULL) // If there were unused groups in source, re-link them and remove their capacity count from *this while adding it to source:
		{
			size_type source_unused_groups_capacity = 0;

			// Count capacity in source unused_groups:
			for (group_pointer_type current = original_unused_groups_head; current != NULL; current = current->next_group)
			{
				source_unused_groups_capacity += current->capacity;
			}

			total_capacity -= source_unused_groups_capacity;
			source.total_capacity = source_unused_groups_capacity;

			// Establish first group from source unused_groups as first active group in source, link rest as reserved groups:
			source.unused_groups_head = original_unused_groups_head->next_group;
			source.begin_iterator.group_pointer = original_unused_groups_head;
			source.begin_iterator.set_to_front_of_group();
			source.end_iterator = source.begin_iterator;
			original_unused_groups_head->reset(0, NULL, NULL, 0);
		}
	}



	#ifdef PLF_MOVE_SEMANTICS_SUPPORT
		void splice(colony &&source)
		{
			splice(source);
		}
	#endif



private:


	struct item_index_tuple
	{
		pointer original_location;
		size_type original_index;

		item_index_tuple(const pointer _item, const size_type _index) PLF_NOEXCEPT:
			original_location(_item),
			original_index(_index)
		{}
	};



	template <class comparison_function>
	struct sort_dereferencer
	{
		comparison_function stored_instance;

		explicit sort_dereferencer(const comparison_function &function_instance):
			stored_instance(function_instance)
		{}

		bool operator() (const item_index_tuple first, const item_index_tuple second)
		{
			return stored_instance(*(first.original_location), *(second.original_location));
		}
	};




	// Try and find space in the unused blocks or the back block instead of allocating for sort:
	template <class the_type>
	aligned_pointer_type get_free_space() const PLF_NOEXCEPT
	{
		const size_type number_of_elements_needed = ((total_size * sizeof(the_type)) + sizeof(aligned_element_struct) - 1) / sizeof(aligned_element_struct); // rounding up

		if (number_of_elements_needed < max_block_capacity)
		{
			if (static_cast<size_type>(end_iterator.group_pointer->past_back() - end_iterator.element_pointer) >= number_of_elements_needed)
			{ // there is enough space at the back of the back block
				return end_iterator.element_pointer;
			}

			for (group_pointer_type current = unused_groups_head; current != NULL; current = current->next_group)
			{
				if (current->capacity >= number_of_elements_needed) return current->front(); // ie. there is enough space in one of the unused blocks
			}
		}

		return NULL;
	}



public:

	template <class comparison_function>
	void sort(comparison_function compare)
	{
		if (total_size < 2) return;

  		#ifdef PLF_TYPE_TRAITS_SUPPORT
			if PLF_CONSTEXPR ((std::is_trivially_copyable<element_type>::value || std::is_move_assignable<element_type>::value) && sizeof(element_type) <= sizeof(element_type *) * 2) // If element is <= 2 pointers, just copy to an array and sort that then copy back - consumes less memory
		#else
			if PLF_CONSTEXPR (sizeof(element_type) <= sizeof(element_type *) * 2)
		#endif
		{
			pointer sort_array = pointer_cast<pointer>(get_free_space<element_type>());
			const bool need_to_allocate = (sort_array == NULL);

			if (need_to_allocate)
			{
				sort_array = PLF_ALLOCATE(allocator_type, *this, total_size, end_iterator.skipfield_pointer);
			}

			const pointer end = sort_array + total_size;

			#if defined(PLF_TYPE_TRAITS_SUPPORT) && defined(PLF_MOVE_SEMANTICS_SUPPORT)
				if PLF_CONSTEXPR (!std::is_trivially_copy_constructible<element_type>::value && std::is_nothrow_move_constructible<element_type>::value)
				{
					plf::uninitialized_move(begin_iterator, end_iterator, sort_array, static_cast<allocator_type &>(*this));
				}
				else
			#endif
			{
				plf::uninitialized_copy(begin_iterator, end_iterator, sort_array, static_cast<allocator_type &>(*this));
			}

			PLF_SORT_FUNCTION(sort_array, end, compare);

			#if defined(PLF_TYPE_TRAITS_SUPPORT) && defined(PLF_MOVE_SEMANTICS_SUPPORT)
				if PLF_CONSTEXPR (!std::is_trivially_copy_assignable<element_type>::value && std::is_nothrow_move_assignable<element_type>::value)
				{
					std::copy(plf::make_move_iterator(sort_array), plf::make_move_iterator(end), begin_iterator);
				}
				else
			#endif
			{
				std::copy(sort_array, end, begin_iterator);

				#ifdef PLF_TYPE_TRAITS_SUPPORT
					if (!std::is_trivially_destructible<element_type>::value)
				#endif
				{
					for (element_type *current = sort_array; current != end; ++current)
					{
						PLF_DESTROY(allocator_type, *this, current);
					}
				}
			}

			if (need_to_allocate)
			{
				PLF_DEALLOCATE(allocator_type, *this, sort_array, total_size);
			}
		}
 		else
 		{
			tuple_pointer_type sort_array = pointer_cast<tuple_pointer_type>(get_free_space<item_index_tuple>());
			const bool need_to_allocate = (sort_array == NULL);

			if (need_to_allocate)
			{
				sort_array = PLF_ALLOCATE(tuple_allocator_type, tuple_allocator, total_size, end_iterator.skipfield_pointer);
			}

			tuple_pointer_type tuple_pointer = sort_array;

			// Construct pointers to all elements in the sequence:
			size_type index = 0;

			for (iterator current_element = begin_iterator; current_element != end_iterator; ++current_element, ++tuple_pointer, ++index)
			{
				#ifdef PLF_VARIADICS_SUPPORT
					PLF_CONSTRUCT(tuple_allocator_type, tuple_allocator, tuple_pointer, &*current_element, index);
				#else
					PLF_CONSTRUCT(tuple_allocator_type, tuple_allocator, tuple_pointer, item_index_tuple(&*current_element, index));
				#endif
			}

			// Now, sort the pointers by the values they point to:
			PLF_SORT_FUNCTION(sort_array, tuple_pointer, sort_dereferencer<comparison_function>(compare));

			// Sort the actual elements via the tuple array:
			index = 0;

			for (tuple_pointer_type current_tuple = sort_array; current_tuple != tuple_pointer; ++current_tuple, ++index)
			{
				if (current_tuple->original_index != index)
				{
					#ifdef PLF_MOVE_SEMANTICS_SUPPORT
						element_type end_value = std::move(*(current_tuple->original_location));
					#else
						element_type end_value = *(current_tuple->original_location);
					#endif
					size_type destination_index = index;
					size_type source_index = current_tuple->original_index;

					do
					{
						#ifdef PLF_MOVE_SEMANTICS_SUPPORT
							*(sort_array[destination_index].original_location) = std::move(*(sort_array[source_index].original_location));
						#else
							*(sort_array[destination_index].original_location) = *(sort_array[source_index].original_location);
						#endif
						destination_index = source_index;
						source_index = sort_array[destination_index].original_index;
						sort_array[destination_index].original_index = destination_index;
					} while (source_index != index);

					#ifdef PLF_MOVE_SEMANTICS_SUPPORT
						*(sort_array[destination_index].original_location) = std::move(end_value);
					#else
						*(sort_array[destination_index].original_location) = end_value;
					#endif
				}
			}

			if (need_to_allocate)
			{
				PLF_DEALLOCATE(tuple_allocator_type, tuple_allocator, sort_array, total_size);
			}
		}
	}



	void sort()
	{
		sort(plf::less<element_type>());
	}



	template <class comparison_function>
	size_type unique(comparison_function compare)
	{
		if (total_size < 2) return 0;

		size_type count = 0;
		const const_iterator end = end_iterator;

		for(const_iterator current = begin_iterator, previous = begin_iterator; ++current != end; previous = current)
		{
			if (compare(*current, *previous))
			{
				const size_type original_count = ++count;
				const_iterator last = current;

				while(++last != end && compare(*last, *previous))
				{
					++count;
				}

				if (count != original_count)
				{
					current = erase(current, last); // optimised range-erase
				}
				else
				{
					current = erase(current);
				}

				if (last == end) break;
			}
		}

		return count;
	}



	size_type unique()
	{
		return unique(std::equal_to<element_type>());
	}



	void swap(colony &source) PLF_NOEXCEPT_SWAP(allocator_type)
	{
		assert(&source != this);

		#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT
			if PLF_CONSTEXPR (std::allocator_traits<allocator_type>::is_always_equal::value && std::is_trivially_copyable<group_pointer_type>::value) // if all pointer types are trivial we can just copy using memcpy - avoids constructors/destructors etc and is faster
			{
				char temp[sizeof(colony)];
				std::memcpy(static_cast<void *>(&temp), static_cast<void *>(this), sizeof(colony));
				std::memcpy(static_cast<void *>(this), static_cast<void *>(&source), sizeof(colony));
				std::memcpy(static_cast<void *>(&source), static_cast<void *>(&temp), sizeof(colony));
			}
			#ifdef PLF_MOVE_SEMANTICS_SUPPORT
				else if PLF_CONSTEXPR (std::is_move_assignable<group_pointer_type>::value && std::is_move_constructible<group_pointer_type>::value)
				{
					colony temp(std::move(source));
					source = std::move(*this);
					*this = std::move(temp);
				}
				else
			#endif
		#endif
		{
			// Otherwise, make the reads/writes as contiguous in memory as-possible (yes, it is faster than using std::swap with the individual variables):
			const iterator 				swap_end_iterator = end_iterator, swap_begin_iterator = begin_iterator;
			const group_pointer_type	swap_erasure_groups_head = erasure_groups_head, swap_unused_groups_head = unused_groups_head;
			const size_type				swap_total_size = total_size, swap_total_capacity = total_capacity;
			const skipfield_type 		swap_min_block_capacity = min_block_capacity, swap_max_block_capacity = max_block_capacity;

			end_iterator = source.end_iterator;
			begin_iterator = source.begin_iterator;
			erasure_groups_head = source.erasure_groups_head;
			unused_groups_head = source.unused_groups_head;
			total_size = source.total_size;
			total_capacity = source.total_capacity;
			min_block_capacity = source.min_block_capacity;
			max_block_capacity = source.max_block_capacity;

			source.end_iterator = swap_end_iterator;
			source.begin_iterator = swap_begin_iterator;
			source.erasure_groups_head = swap_erasure_groups_head;
			source.unused_groups_head = swap_unused_groups_head;
			source.total_size = swap_total_size;
			source.total_capacity = swap_total_capacity;
			source.min_block_capacity = swap_min_block_capacity;
			source.max_block_capacity = swap_max_block_capacity;

			#ifdef PLF_IS_ALWAYS_EQUAL_SUPPORT
				if PLF_CONSTEXPR (std::allocator_traits<allocator_type>::propagate_on_container_swap::value && !std::allocator_traits<allocator_type>::is_always_equal::value)
			#endif
			{
				std::swap(static_cast<allocator_type &>(source), static_cast<allocator_type &>(*this));

				// Reconstruct rebinds for swapped allocators:
				group_allocator = group_allocator_type(*this);
				aligned_struct_allocator = aligned_struct_allocator_type(*this);
				skipfield_allocator = skipfield_allocator_type(*this);
				tuple_allocator = tuple_allocator_type(*this);
				source.group_allocator = group_allocator_type(source);
				source.aligned_struct_allocator = aligned_struct_allocator_type(source);
				source.skipfield_allocator = skipfield_allocator_type(source);
				source.tuple_allocator = tuple_allocator_type(source);
			} // else: undefined behaviour, as per standard
		}
	}



	// Because it's going to take a lot of defining of very specific allocator and pointer types to make this work with smart pointers - this class and functions associated with it will not work if the allocator supplied to the container returns smart pointers.
	struct colony_data : public uchar_allocator_type
	{
		aligned_pointer_type * const block_pointers; 	// array of pointers to element memory blocks
		unsigned char * * const bitfield_pointers;		// array of pointers to bitfields in the form of unsigned char arrays representing whether an element is erased or not (0 for erased).
		size_t * const block_capacities; 				// array of the number of elements in each memory block
		const size_t number_of_blocks;					// size of each of the arrays above


		colony_data(const typename colony::size_type size) :
			block_pointers(pointer_cast<aligned_pointer_type *>(PLF_ALLOCATE(uchar_allocator_type, *this, size * sizeof(aligned_pointer_type), NULL))),
			bitfield_pointers(pointer_cast<unsigned char **>(PLF_ALLOCATE(uchar_allocator_type, *this, size * sizeof(unsigned char *), NULL))),
			block_capacities(pointer_cast<size_t *>(PLF_ALLOCATE(uchar_allocator_type, *this, size * sizeof(size_t), NULL))),
			number_of_blocks(size)
		{}


		~colony_data()
		{
			for (size_t index = 0; index != number_of_blocks; ++index)
			{
				PLF_DEALLOCATE(uchar_allocator_type, *this, bitfield_pointers[index], (block_capacities[index] + 7) / 8);
			}

			PLF_DEALLOCATE(uchar_allocator_type, *this, pointer_cast<unsigned char *>(block_pointers), number_of_blocks * sizeof(aligned_pointer_type));
			PLF_DEALLOCATE(uchar_allocator_type, *this, pointer_cast<unsigned char *>(bitfield_pointers), number_of_blocks * sizeof(unsigned char *));
			PLF_DEALLOCATE(uchar_allocator_type, *this, pointer_cast<unsigned char *>(block_capacities), number_of_blocks * sizeof(size_t));
		}
	};



private:

	void setup_data_cell(colony_data *data, const group_pointer_type current_group, const size_t capacity, const size_t group_number)
	{
		const size_t bitfield_capacity = (capacity + 7) / 8; // round up

		data->block_pointers[group_number] = current_group->front();
		unsigned char *bitfield_location = data->bitfield_pointers[group_number] = PLF_ALLOCATE(uchar_allocator_type, (*data), bitfield_capacity, NULL);
		data->block_capacities[group_number] = capacity;
		std::memset(bitfield_location, 0, bitfield_capacity);

		skipfield_pointer_type skipfield_pointer = current_group->skipfield;
		const unsigned char * const end = bitfield_location + bitfield_capacity;

		for (size_t index = 0; bitfield_location != end; ++bitfield_location)
		{
			for (unsigned char offset = 0; offset != 8 && index != capacity; ++index, ++offset, ++skipfield_pointer)
			{
				*bitfield_location |= static_cast<unsigned char>(static_cast<int>(!*skipfield_pointer) << offset);
			}
		}
	}



public:

	colony_data * data()
	{
		colony_data *data = new colony_data(end_iterator.group_pointer->group_number + 1);
		size_t group_number = 0;

		for (group_pointer_type current_group = begin_iterator.group_pointer; current_group != end_iterator.group_pointer; current_group = current_group->next_group, ++group_number)
		{
			setup_data_cell(data, current_group, current_group->capacity, group_number);
		}

		// Special case for end group:
		setup_data_cell(data, end_iterator.group_pointer, static_cast<size_t>(end_iterator.element_pointer - end_iterator.group_pointer->front()), group_number);

		return data;
	}




	// Iterators:
	template <bool is_const>
	class colony_iterator
	{
	private:
		typedef typename colony::group_pointer_type 		group_pointer_type;
		typedef typename colony::aligned_pointer_type 	aligned_pointer_type;
		typedef typename colony::skipfield_pointer_type skipfield_pointer_type;

		#ifdef PLF_DEFAULT_SUPPORT
			group_pointer_type		group_pointer{NULL};
			aligned_pointer_type 	element_pointer{NULL};
			skipfield_pointer_type	skipfield_pointer{NULL};
		#else
			group_pointer_type		group_pointer;
			aligned_pointer_type 	element_pointer;
			skipfield_pointer_type	skipfield_pointer;
		#endif


		void set_to_front_of_group() PLF_NOEXCEPT
		{
			element_pointer = group_pointer->front();
			skipfield_pointer = group_pointer->skipfield;
		}



		void set_to_first_element_in_group() PLF_NOEXCEPT
		{
			element_pointer = group_pointer->front() + *(group_pointer->skipfield);
			skipfield_pointer = group_pointer->skipfield + *(group_pointer->skipfield);
		}



		void set_to_past_back_of_group() PLF_NOEXCEPT
		{
			element_pointer = group_pointer->past_back();
			skipfield_pointer = group_pointer->skipfield + group_pointer->capacity;
		}


	public:
		struct colony_iterator_tag {};
		typedef std::bidirectional_iterator_tag	iterator_category;
		typedef std::bidirectional_iterator_tag	iterator_concept;
		typedef typename colony::value_type 			value_type;
		typedef typename colony::difference_type		difference_type;
		typedef colony_reverse_iterator<is_const> 	reverse_type;
		typedef typename plf::conditional<is_const, typename colony::const_pointer, typename colony::pointer>::type		pointer;
		typedef typename plf::conditional<is_const, typename colony::const_reference, typename colony::reference>::type	reference;

		friend class colony;
		friend class colony_reverse_iterator<false>;
		friend class colony_reverse_iterator<true>;

		// Friend functions:

		template <class distance_type>
		friend void advance(colony_iterator &it, distance_type distance)
		{
			it.advance(static_cast<difference_type>(distance));
		}



		friend colony_iterator next(const colony_iterator &it, const difference_type distance)
		{
			colony_iterator return_iterator(it);
			return_iterator.advance(static_cast<difference_type>(distance));
			return return_iterator;
		}



		friend colony_iterator prev(const colony_iterator &it, const difference_type distance)
		{
			colony_iterator return_iterator(it);
			return_iterator.advance(-(static_cast<difference_type>(distance)));
			return return_iterator;
		}



		friend typename colony_iterator::difference_type distance(const colony_iterator &first, const colony_iterator &last)
		{
			return first.distance(last);
		}



		colony_iterator() PLF_NOEXCEPT
		#ifdef PLF_DEFAULT_SUPPORT
			= default;
		#else
			: group_pointer(NULL),
			element_pointer(NULL),
			skipfield_pointer(NULL)
			{}
		#endif



		colony_iterator (const colony_iterator &source) PLF_NOEXCEPT
		#ifdef PLF_DEFAULT_SUPPORT
			= default;
		#else
			: group_pointer(source.group_pointer),
			element_pointer(source.element_pointer),
			skipfield_pointer(source.skipfield_pointer)
			{}
		#endif



		#ifdef PLF_DEFAULT_SUPPORT
			template <bool is_const_it = is_const, class = typename plf::enable_if<is_const_it>::type >
			colony_iterator(const colony_iterator<false> &source) PLF_NOEXCEPT:
		#else
			colony_iterator(const colony_iterator<!is_const> &source) PLF_NOEXCEPT:
		#endif
			group_pointer(source.group_pointer),
			element_pointer(source.element_pointer),
			skipfield_pointer(source.skipfield_pointer)
		{}



		#ifdef PLF_MOVE_SEMANTICS_SUPPORT
			colony_iterator(colony_iterator &&source) PLF_NOEXCEPT
			#ifdef PLF_DEFAULT_SUPPORT
				= default;
			#else
				: group_pointer(std::move(source.group_pointer)),
				element_pointer(std::move(source.element_pointer)),
				skipfield_pointer(std::move(source.skipfield_pointer))
				{}
			#endif



			#ifdef PLF_DEFAULT_SUPPORT
				template <bool is_const_it = is_const, class = typename plf::enable_if<is_const_it>::type >
				colony_iterator(colony_iterator<false> &&source) PLF_NOEXCEPT:
			#else
				colony_iterator(colony_iterator<!is_const> &&source) PLF_NOEXCEPT:
			#endif
				group_pointer(std::move(source.group_pointer)),
				element_pointer(std::move(source.element_pointer)),
				skipfield_pointer(std::move(source.skipfield_pointer))
			{}
		#endif



		colony_iterator & operator = (const colony_iterator &source) PLF_NOEXCEPT
		#ifdef PLF_DEFAULT_SUPPORT
			= default;
		#else
			{
				group_pointer = source.group_pointer;
				element_pointer = source.element_pointer;
				skipfield_pointer = source.skipfield_pointer;
				return *this;
			}
		#endif



		#ifdef PLF_DEFAULT_SUPPORT
			template <bool is_const_it = is_const, class = typename plf::enable_if<is_const_it>::type >
			colony_iterator & operator = (const colony_iterator<false> &source) PLF_NOEXCEPT
		#else
			colony_iterator & operator = (const colony_iterator<!is_const> &source) PLF_NOEXCEPT
		#endif
		{
			group_pointer = source.group_pointer;
			element_pointer = source.element_pointer;
			skipfield_pointer = source.skipfield_pointer;
			return *this;
		}



		#ifdef PLF_MOVE_SEMANTICS_SUPPORT
			colony_iterator & operator = (colony_iterator &&source) PLF_NOEXCEPT
			#ifdef PLF_DEFAULT_SUPPORT
				= default;
			#else
				{
					assert(&source != this);
					group_pointer = std::move(source.group_pointer);
					element_pointer = std::move(source.element_pointer);
					skipfield_pointer = std::move(source.skipfield_pointer);
					return *this;
				}
			#endif



			#ifdef PLF_DEFAULT_SUPPORT
				template <bool is_const_it = is_const, class = typename plf::enable_if<is_const_it>::type >
				colony_iterator & operator = (colony_iterator<false> &&source) PLF_NOEXCEPT
			#else
				colony_iterator & operator = (colony_iterator<!is_const> &&source) PLF_NOEXCEPT
			#endif
			{
				group_pointer = std::move(source.group_pointer);
				element_pointer = std::move(source.element_pointer);
				skipfield_pointer = std::move(source.skipfield_pointer);
				return *this;
			}
		#endif



		bool operator == (const colony_iterator &rh) const PLF_NOEXCEPT
		{
			return (element_pointer == rh.element_pointer);
		}



		bool operator == (const colony_iterator<!is_const> &rh) const PLF_NOEXCEPT
		{
			return (element_pointer == rh.element_pointer);
		}



		bool operator != (const colony_iterator &rh) const PLF_NOEXCEPT
		{
			return (element_pointer != rh.element_pointer);
		}



		bool operator != (const colony_iterator<!is_const> &rh) const PLF_NOEXCEPT
		{
			return (element_pointer != rh.element_pointer);
		}



		reference operator * () const // may cause exception with uninitialized iterator
		{
			return *pointer_cast<pointer>(element_pointer);
		}



		pointer operator -> () const
		{
			return pointer_cast<pointer>(element_pointer);
		}



		colony_iterator & operator ++ ()
		{
			assert(group_pointer != NULL); // covers uninitialised colony_iterator

			#if defined(_MSC_VER) && !defined(__clang__) && !defined(__GNUC__) // This version creates significantly faster release code under MSVC in some scenarios, but not in any other compilers
				element_pointer += *(++skipfield_pointer) + 1;
				skipfield_pointer += *skipfield_pointer;

				if (element_pointer == group_pointer->past_back() && group_pointer->next_group != NULL) // Second condition allows iterator to reach end(), which may be 1 past end of block, if block has been fully used
				{
					group_pointer = group_pointer->next_group;
					set_to_front_of_group();
					element_pointer += *skipfield_pointer;
					skipfield_pointer += *skipfield_pointer;
				}
			#else
				skipfield_type skip = *(++skipfield_pointer);

				if ((element_pointer += static_cast<size_type>(skip) + 1u) == group_pointer->past_back() && group_pointer->next_group != NULL)
				{
					group_pointer = group_pointer->next_group;
					const aligned_pointer_type elements = group_pointer->front();
					const skipfield_pointer_type skipfield = group_pointer->skipfield;
					skip = *skipfield;
					element_pointer = elements + skip;
					skipfield_pointer = skipfield;
				}

				skipfield_pointer += skip;
			#endif

			return *this;
		}



		colony_iterator operator ++(int)
		{
			const colony_iterator copy(*this);
			++*this;
			return copy;
		}



		colony_iterator & operator -- ()
		{
			assert(group_pointer != NULL);

			if (--skipfield_pointer >= group_pointer->skipfield) // ie. not already at beginning of group prior to decrementation
			{
				element_pointer -= static_cast<size_type>(*skipfield_pointer) + 1u;
				if ((skipfield_pointer -= *skipfield_pointer) >= group_pointer->skipfield) return *this; // ie. skipfield jump value does not takes us beyond beginning of group
			}

			group_pointer = group_pointer->previous_group;
			const skipfield_pointer_type skipfield = group_pointer->skipfield + group_pointer->capacity - 1;
			const skipfield_type skip = *skipfield;
			element_pointer = (group_pointer->past_back() - 1) - skip;
			skipfield_pointer = skipfield - skip;
			return *this;
		}



		colony_iterator operator -- (int)
		{
			const colony_iterator copy(*this);
			--*this;
			return copy;
		}



		// Less-than etc operators retained as GCC codegen synthesis from <=> is slower and bulkier for same operations:
		template <bool is_const_it>
		bool operator > (const colony_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return ((group_pointer == rh.group_pointer) & (element_pointer > rh.element_pointer)) ||
				(group_pointer != rh.group_pointer && group_pointer->group_number > rh.group_pointer->group_number);
		}



		template <bool is_const_it>
		bool operator < (const colony_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return rh > *this;
		}



		template <bool is_const_it>
		bool operator >= (const colony_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return !(rh > *this);
		}



		template <bool is_const_it>
		bool operator <= (const colony_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return !(*this > rh);
		}



		#ifdef PLF_CPP20_SUPPORT
			template <bool is_const_it>
			std::strong_ordering operator <=> (const colony_iterator<is_const_it> &rh) const noexcept
			{
				return (element_pointer == rh.element_pointer) ? std::strong_ordering::equal : ((*this > rh) ? std::strong_ordering::greater : std::strong_ordering::less);
			}
		#endif



	private:
		// Used by cend(), erase() etc:
		colony_iterator(const group_pointer_type group_p, const aligned_pointer_type element_p, const skipfield_pointer_type skipfield_p) PLF_NOEXCEPT:
			group_pointer(group_p),
			element_pointer(element_p),
			skipfield_pointer(skipfield_p)
		{}



		// These 3 functions are used by advance:

		void set_element_location_from_skipfield() PLF_NOEXCEPT
		{
			element_pointer = group_pointer->front() + (skipfield_pointer - group_pointer->skipfield);
		}



		void advance_from_group_beginning(skipfield_type distance) PLF_NOEXCEPT
		{
			if (group_pointer->no_erasures())
			{
				element_pointer = group_pointer->front() + distance;
				skipfield_pointer = group_pointer->skipfield + distance;
			}
			else
			{
				skipfield_pointer = group_pointer->skipfield + *(group_pointer->skipfield); // in case first element(s) in block are erased

				do
				{
					++skipfield_pointer;
					skipfield_pointer += *skipfield_pointer;
				} while(--distance != 0);

				set_element_location_from_skipfield();
			}
		}



		aligned_pointer_type get_end() const PLF_NOEXCEPT
		{
			#ifndef PLF_COLONY_DISABLE_EXPENSIVE_BOUNDING
				if (group_pointer->next_group == NULL)
				{
					colony_iterator current(group_pointer, NULL, 0);
					current.advance_from_group_beginning(group_pointer->size);
					return current.element_pointer;
				}
				else
			#endif
			{
				return group_pointer->past_back();
			}
		}



		// Advance implementation:

		void advance(difference_type distance) // Cannot be noexcept due to the possibility of an uninitialized iterator
		{
			assert(group_pointer != NULL); // covers uninitialized colony_iterator && empty group

			// Now, run code based on the nature of the distance type - negative, positive or zero:
			if (distance > 0) // ie. +=
			{
				// Code explanation:
				// For the initial state of the iterator, we don't know which elements have been erased before that element in that group.
				// So for the first group, we follow the following logic:
				// 1. If no elements have been erased in the group, we do simple pointer addition to progress, either to within the group (if the distance is small enough) or the end of the group and subtract from distance accordingly.
				// 2. If any of the first group's elements have been erased, we manually iterate, as we don't know whether the erased elements occur before or after the initial iterator position, and we subtract 1 from the distance amount each time we iterate. Iteration continues until either distance becomes zero, or we reach the end of the group.

				// For all subsequent groups, we follow this logic:
				// 1. If distance is larger than the total number of non-erased elements in a group, we skip that group and subtract the number of elements in that group from distance.
				// 2. If distance is smaller than the total number of non-erased elements in a group, then:
				//   a. If there are no erased elements in the group we simply add distance to group->elements to find the new location for the iterator.
				//   b. If there are erased elements in the group, we manually iterate and subtract 1 from distance on each iteration, until the new iterator location is found ie. distance = 0.

				// Note: incrementing element_pointer is avoided until necessary to avoid needless calculations.

				aligned_pointer_type end = get_end();
				if (element_pointer >= end) return; // Bound to end()

				// Special case for initial element pointer and initial group (we don't know how far into the group the element pointer is)
				if (element_pointer != group_pointer->first_element()) // Otherwise we skip this section and treat the first block as we would an intermediary block
				{
					const difference_type distance_from_end = end - element_pointer;

					if (group_pointer->no_erasures())
					{
						if (distance < distance_from_end)
						{
							element_pointer += distance;
							skipfield_pointer += distance;
							return;
						} // distance >= distance_from_end
						else if (group_pointer->next_group == NULL) // either we're going to reach end() or go beyond it, so bound to back of block
						{
							element_pointer = end;
							skipfield_pointer += distance_from_end;
							return;
						}
						else
						{
							distance -= distance_from_end;
						}
					}
					else
					{
						const skipfield_pointer_type endpoint = skipfield_pointer + distance_from_end;

						while(true)
						{
							++skipfield_pointer;
							skipfield_pointer += *skipfield_pointer;
							--distance;

							if (skipfield_pointer == endpoint)
							{
								break;
							}
							else if (distance == 0)
							{
								set_element_location_from_skipfield();
								return;
							}
						}

						if (group_pointer->next_group == NULL)
						{
							element_pointer = end;
							return;
						}
					}

					group_pointer = group_pointer->next_group;

					if (distance == 0)
					{
						set_to_first_element_in_group();
						return;
					}
				}


				// Intermediary groups - at the start of this code block and the subsequent block, the position of the iterator is assumed to be the first non-erased element in the current group:
				while (static_cast<difference_type>(group_pointer->size) <= distance)
				{
					if (group_pointer->next_group == NULL) // either we would reach end() or go beyond it, so bound to end() - necessary to reach rbegin()
					{
						advance_from_group_beginning(group_pointer->size);
						return;
					}

					distance -= group_pointer->size;
					group_pointer = group_pointer->next_group;

					if (distance == 0)
					{
						set_to_first_element_in_group();
						return;
					}
				}


				// Final group (if not already reached):
				// We already know group size > distance due to the intermediary group checks above - so it's safe to ignore endpoint check condition while incrementing here.
				advance_from_group_beginning(static_cast<skipfield_type>(distance));
			}
			else if (distance < 0)
			{
				if (group_pointer->previous_group == NULL && element_pointer == group_pointer->first_element()) return; // if we are at begin(), bound to that

				distance = -distance;

				// Special case for initial element pointer and initial group (we don't know how far into the group the element pointer is)
 				if (element_pointer != group_pointer->past_back()) // Optimization: if end() is cheaply calculable (edge case, is one past end of block), and iterator == end(), skip this first section and treat current group like an intermediary group. Note: only end() can be 1-past end of block, no other iterator.
				{
					if (group_pointer->no_erasures())
					{
						const difference_type distance_from_beginning = static_cast<difference_type>(element_pointer - group_pointer->front());

						if (distance <= distance_from_beginning)
						{
							element_pointer -= distance;
							skipfield_pointer -= distance;
							return;
						}
						else if (group_pointer->previous_group == NULL) // we've gone before begin(), so bound to that
						{
							element_pointer = group_pointer->front();
							skipfield_pointer = group_pointer->skipfield;
							return;
						}
						else
						{
							distance -= distance_from_beginning;
						}
					}
					else
					{
						for (const skipfield_pointer_type begin = group_pointer->skipfield + *(group_pointer->skipfield); skipfield_pointer != begin;)
						{
							--skipfield_pointer;
							skipfield_pointer -= *skipfield_pointer;

							if (--distance == 0)
							{
								set_element_location_from_skipfield();
								return;
							}
						}

						if (group_pointer->previous_group == NULL)
						{
							set_to_first_element_in_group(); // this is first group, so bound to begin()
							return;
						}
					}

					group_pointer = group_pointer->previous_group;
				}


				// Intermediary groups - at the start of this code block and the subsequent block, the position of the iterator is assumed to be either the first non-erased element in the next group over, or end():
				while(static_cast<difference_type>(group_pointer->size) < distance)
				{
					if (group_pointer->previous_group == NULL) // we've gone beyond begin(), so bound to it
					{
						set_to_first_element_in_group();
						return;
					}

					distance -= group_pointer->size;
					group_pointer = group_pointer->previous_group;
				}


				// Final group (if not already reached above):
				if (static_cast<difference_type>(group_pointer->size) == distance) // go to front of group
				{
					set_to_first_element_in_group();
				}
				else if (group_pointer->no_erasures())
				{
					element_pointer = group_pointer->past_back() - distance;
					skipfield_pointer = (group_pointer->skipfield + group_pointer->size) - distance;
				}
				else // ie. no more groups to traverse but there are erased elements in this group
				{
					skipfield_pointer = group_pointer->skipfield + group_pointer->capacity;

					do
					{
						--skipfield_pointer;
						skipfield_pointer -= *skipfield_pointer;
					} while(--distance != 0);

					set_element_location_from_skipfield();
				}
			}
		}



		// distance implementation:

		difference_type distance(const colony_iterator &last) const
		{
			// Code logic:
			// If iterators are the same, return 0
			// Otherwise, find which iterator is later in colony, copy that to iterator2. Copy the lower to iterator1.
			// If they are not pointing to elements in the same group, process the intermediate groups and add distances,
			// skipping manual incrementation in all but the initial and final groups.
			// In the initial and final groups, manual incrementation must be used to calculate distance, if there have been any erasures in those groups.
			// If there are no prior erasures in either of those groups we can use pointer arithmetic to calculate the distances.

			assert(!(group_pointer == NULL) && !(last.group_pointer == NULL));  // Check that they are both initialized

			difference_type distance = 0;
			colony_iterator iterator1 = *this, iterator2 = last;
			const bool swap_iterators = iterator1 > iterator2;

			if (swap_iterators)
			{
				iterator1 = last;
				iterator2 = *this;
			}

			if (iterator1.group_pointer != iterator2.group_pointer) // if not in same group, process intermediate groups
			{
				// Process initial group:
				if (iterator1.group_pointer->no_erasures()) // ie. we can do simple addition
				{
					distance += static_cast<difference_type>(iterator1.group_pointer->past_back() - iterator1.element_pointer);
				}
				else if (iterator1.element_pointer == iterator1.group_pointer->first_element()) // edge case
				{
					distance += static_cast<difference_type>(iterator1.group_pointer->size);
				}
				else // Manually iterate to find distance to end of group:
				{
					for (const skipfield_pointer_type end = iterator1.skipfield_pointer + (iterator1.group_pointer->past_back() - iterator1.element_pointer); iterator1.skipfield_pointer != end;)
					{
						++iterator1.skipfield_pointer;
						iterator1.skipfield_pointer += *iterator1.skipfield_pointer;
						++distance;
					}
				}

				// Process all other intermediate groups:
				iterator1.group_pointer = iterator1.group_pointer->next_group;

				while (iterator1.group_pointer != iterator2.group_pointer)
				{
					distance += static_cast<difference_type>(iterator1.group_pointer->size);
					iterator1.group_pointer = iterator1.group_pointer->next_group;
				}

				iterator1.skipfield_pointer = iterator1.group_pointer->skipfield + *(iterator1.group_pointer->skipfield);
			}


			if (iterator2.group_pointer->no_erasures()) // ie. simple subtraction is possible
			{
				distance += iterator2.skipfield_pointer - iterator1.skipfield_pointer;
			}
 			else if (iterator1.element_pointer == iterator2.group_pointer->first_element() && iterator2.element_pointer + 1 + *(iterator2.skipfield_pointer + 1) == iterator2.group_pointer->past_back()) // ie. if iterator1 is at beginning of block (have to check this in case first and last are in the same block to begin with) and iterator2 is last element in the block. This check won't work for back block (as end() may not be past end of block), but that's fine, it's just an optimization.
			{
				distance += static_cast<difference_type>(iterator2.group_pointer->size) - 1;
			}
			else
			{
				while (iterator1.skipfield_pointer != iterator2.skipfield_pointer)
				{
					++iterator1.skipfield_pointer;
					iterator1.skipfield_pointer += *iterator1.skipfield_pointer;
					++distance;
				}
			}


			if (swap_iterators) distance = -distance;

			return distance;
		}
	}; // colony_iterator




	// Reverse iterators:
	// Note: despite iterator being trivially_copyable, defining reverse_iterator via std::reverse_iterator<iterator> results in it not being trivially_copyable here, at least under libstdc++.
	// Hence we have a custom implementation.

	template <bool is_const_r>
	class colony_reverse_iterator
	{
	protected:
		iterator current;

	public:
		struct colony_iterator_tag {};
		typedef std::bidirectional_iterator_tag 	iterator_category;
		typedef std::bidirectional_iterator_tag 	iterator_concept;
		typedef iterator 									iterator_type;
		typedef typename colony::value_type 		value_type;
		typedef typename colony::difference_type	difference_type;
		typedef typename plf::conditional<is_const_r, typename colony::const_pointer, typename colony::pointer>::type		pointer;
		typedef typename plf::conditional<is_const_r, typename colony::const_reference, typename colony::reference>::type	reference;

		friend class colony;


		template <class distance_type>
		friend void advance(colony_reverse_iterator &it, const distance_type distance)
		{
			it.advance(static_cast<difference_type>(distance));
		}



		friend colony_reverse_iterator next(const colony_reverse_iterator &it, const difference_type distance)
		{
			colony_reverse_iterator return_iterator(it);
			return_iterator.advance(static_cast<difference_type>(distance));
			return return_iterator;
		}



		friend colony_reverse_iterator prev(const colony_reverse_iterator &it, const difference_type distance)
		{
			colony_reverse_iterator return_iterator(it);
			return_iterator.advance(static_cast<difference_type>(-distance));
			return return_iterator;
		}



		friend typename colony_reverse_iterator::difference_type distance(const colony_reverse_iterator &first, const colony_reverse_iterator &last)
		{
			return first.distance(last);
		}



		colony_reverse_iterator (const colony_reverse_iterator &source) PLF_NOEXCEPT
		#ifdef PLF_DEFAULT_SUPPORT
			= default;
		#else
			: current(source.current) {}
		#endif


		#ifdef PLF_DEFAULT_SUPPORT
			template <bool is_const_rit = is_const_r, class = typename plf::enable_if<is_const_rit>::type >
			colony_reverse_iterator (const colony_reverse_iterator<false> &source) PLF_NOEXCEPT:
		#else
			colony_reverse_iterator (const colony_reverse_iterator<!is_const_r> &source) PLF_NOEXCEPT:
		#endif
			current(source.current)
		{}


		colony_reverse_iterator (const colony_iterator<is_const_r> &source) PLF_NOEXCEPT:
			current(source)
		{}


		#ifdef PLF_DEFAULT_SUPPORT
			template <bool is_const_rit = is_const_r, class = typename plf::enable_if<is_const_rit>::type >
			colony_reverse_iterator (const colony_iterator<false> &source) PLF_NOEXCEPT:
		#else
			colony_reverse_iterator (const colony_iterator<!is_const_r> &source) PLF_NOEXCEPT:
		#endif
			current(source)
		{}


		#ifdef PLF_MOVE_SEMANTICS_SUPPORT
			colony_reverse_iterator (colony_reverse_iterator &&source) PLF_NOEXCEPT
			#ifdef PLF_DEFAULT_SUPPORT
				= default;
			#else
				: current(std::move(source.current)) {}
			#endif


			#ifdef PLF_DEFAULT_SUPPORT
				template <bool is_const_rit = is_const_r, class = typename plf::enable_if<is_const_rit>::type >
				colony_reverse_iterator (colony_reverse_iterator<false> &&source) PLF_NOEXCEPT:
			#else
				colony_reverse_iterator (colony_iterator<!is_const_r> &&source) PLF_NOEXCEPT:
			#endif
				current(std::move(source.current))
			{}
		#endif


		colony_reverse_iterator& operator = (const colony_iterator<is_const_r> &source) PLF_NOEXCEPT
		{
			current = source;
			++current;
			return *this;
		}


		#ifdef PLF_DEFAULT_SUPPORT
			template <bool is_const_rit = is_const_r, class = typename plf::enable_if<is_const_rit>::type >
			colony_reverse_iterator& operator = (const colony_iterator<false> &source) PLF_NOEXCEPT
		#else
			colony_reverse_iterator& operator = (const colony_iterator<!is_const_r> &source) PLF_NOEXCEPT
		#endif
		{
			current = source;
			++current;
			return *this;
		}


		colony_reverse_iterator& operator = (const colony_reverse_iterator &source) PLF_NOEXCEPT
		#ifdef PLF_DEFAULT_SUPPORT
			= default;
		#else
			{
				current = source.current;
				return *this;
			}
		#endif


		#ifdef PLF_DEFAULT_SUPPORT
			template <bool is_const_rit = is_const_r, class = typename plf::enable_if<is_const_rit>::type >
			colony_reverse_iterator& operator = (const colony_reverse_iterator<false> &source) PLF_NOEXCEPT
		#else
			colony_reverse_iterator& operator = (const colony_reverse_iterator<!is_const_r> &source) PLF_NOEXCEPT
		#endif
		{
			current = source.current;
			return *this;
		}


		#ifdef PLF_MOVE_SEMANTICS_SUPPORT
			colony_reverse_iterator& operator = (colony_reverse_iterator &&source) PLF_NOEXCEPT
			#ifdef PLF_DEFAULT_SUPPORT
				= default;
			#else
				{
					assert(&source != this);
					current = std::move(source.current);
					return *this;
				}
			#endif


			#ifdef PLF_DEFAULT_SUPPORT
				template <bool is_const_rit = is_const_r, class = typename plf::enable_if<is_const_rit>::type >
				colony_reverse_iterator& operator = (colony_reverse_iterator<false> &&source) PLF_NOEXCEPT
			#else
				colony_reverse_iterator& operator = (colony_reverse_iterator<!is_const_r> &&source) PLF_NOEXCEPT
			#endif
			{
				assert(&source != this);
				current = std::move(source.current);
				return *this;
			}
		#endif



		bool operator == (const colony_reverse_iterator &rh) const PLF_NOEXCEPT
		{
			return (current == rh.current);
		}



		bool operator == (const colony_reverse_iterator<!is_const_r> &rh) const PLF_NOEXCEPT
		{
			return (current == rh.current);
		}



		bool operator != (const colony_reverse_iterator &rh) const PLF_NOEXCEPT
		{
			return (current != rh.current);
		}



		bool operator != (const colony_reverse_iterator<!is_const_r> &rh) const PLF_NOEXCEPT
		{
			return (current != rh.current);
		}



		reference operator * () const PLF_NOEXCEPT
		{
			return *--iterator(current);
		}



		pointer operator -> () const PLF_NOEXCEPT
		{
			return pointer_cast<pointer>(--iterator(current).element_pointer);
		}



		colony_reverse_iterator & operator ++ ()
		{
			--current;
			return *this;
		}



		colony_reverse_iterator operator ++ (int)
		{
			const colony_reverse_iterator copy(*this);
			++*this;
			return copy;
		}



		colony_reverse_iterator & operator -- ()
		{
			++current;
			return *this;
		}



		colony_reverse_iterator operator -- (int)
		{
			const colony_reverse_iterator copy(*this);
			++current;
			return copy;
		}



		colony_iterator<is_const_r> base() const PLF_NOEXCEPT
		{
			return colony_iterator<is_const_r>(current);
		}



		template <bool is_const_it>
		bool operator > (const colony_reverse_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return (rh.current > current);
		}



		template <bool is_const_it>
		bool operator < (const colony_reverse_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return (current > rh.current);
		}



		template <bool is_const_it>
		bool operator >= (const colony_reverse_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return !(current > rh.current);
		}



		template <bool is_const_it>
		bool operator <= (const colony_reverse_iterator<is_const_it> &rh) const PLF_NOEXCEPT
		{
			return !(rh.current > current);
		}



		#ifdef PLF_CPP20_SUPPORT
			template <bool is_const_it>
			std::strong_ordering operator <=> (const colony_reverse_iterator<is_const_it> &rh) const noexcept
			{
				return (rh.current <=> current);
			}
		#endif



	private:

		void advance(const difference_type distance)
		{
 			current.advance(-distance);
		}



 		difference_type distance(const colony_reverse_iterator &last) const
 		{
 			return last.current.distance(current);
 		}
	}; // colony_reverse_iterator


}; // colony



} // plf namespace




namespace std
{


template <class element_type, class allocator_type>
void swap (plf::colony<element_type, allocator_type> &a, plf::colony<element_type, allocator_type> &b) PLF_NOEXCEPT_SWAP(allocator_type)
{
	a.swap(b);
}



template <class element_type, class allocator_type, class predicate_function>
typename plf::colony<element_type, allocator_type>::size_type erase_if(plf::colony<element_type, allocator_type> &container, predicate_function predicate)
{
	typedef typename plf::colony<element_type, allocator_type> colony;
	typedef typename colony::const_iterator 	const_iterator;
	typedef typename colony::size_type 			size_type;
	size_type count = 0;
	const const_iterator end = container.cend();

	for (const_iterator current = container.cbegin(); current != end; ++current)
	{
		if (predicate(*current))
		{
			const size_type original_count = ++count;
			const_iterator last = current;

			while(++last != end && predicate(*last))
			{
				++count;
			}

			if (count != original_count)
			{
				current = container.erase(current, last); // optimised range-erase
			}
			else
			{
				current = container.erase(current);
			}

			if (last == end) break;
		}
	}

	return count;
}



template <class element_type, class allocator_type>
typename plf::colony<element_type, allocator_type>::size_type erase(plf::colony<element_type, allocator_type> &container, const element_type &value)
{
	return erase_if(container, plf::equal_to<element_type>(value));
}



#ifdef PLF_CPP20_SUPPORT
	// std::reverse_iterator overload, to allow use of colony with ranges and make_reverse_iterator primarily:
	template <plf::colony_iterator_concept it_type>
	class reverse_iterator<it_type> : public it_type::reverse_type
	{
	public:
		typedef typename it_type::reverse_type rit;
		using rit::rit;
	};
#endif

} // namespace std



#ifdef PLF_COLONY_DEFINES
	#include "plf_tools_undef.h"
#endif


#ifdef PLF_SORT_FUNCTION_DEFINED
	#undef PLF_SORT_FUNCTION
	#undef PLF_SORT_FUNCTION_DEFINED
#endif

#endif // PLF_COLONY_H
