// NOTE: This test suite throws exceptions to test catch handling. You might see this show up in the debug output when you run it in MSVC. It doesn't mean something has gone wrong.

#define PLF_COLONY_TEST_DEBUG

#include "plf_tools.h"

#include <numeric> // std::accumulate
#include <functional> // std::greater, std::bind2nd
#include <vector> // range-insert testing
#include <algorithm> // std::find
#include <cstdio> // log redirection, printf
#include <cstdlib> // abort

#ifdef PLF_MOVE_SEMANTICS_SUPPORT
	#include <utility> // std::move
#endif

#ifdef PLF_CPP20_SUPPORT
	#include <ranges>
#endif


#include "plf_colony.h"



void message(const char *message_text)
{
	printf("%s\n", message_text);
}


void title1(const char *title_text)
{
	printf("\n\n\n*** %s ***\n", title_text);
	printf("===========================================\n\n\n");
}

void title2(const char *title_text)
{
	printf("\n\n--- %s ---\n\n", title_text);
}


void failpass(const char *test_type, bool condition)
{
	printf("%s: ", test_type);

	if (condition)
	{
		printf("Pass\n");
	}
	else
	{
		printf("Fail\n");
		getchar();
		abort();
	}
}



#ifdef PLF_VARIADICS_SUPPORT
	struct perfect_forwarding_test
	{
		const bool success;

		perfect_forwarding_test(int&& /*perfect1*/, int& perfect2)
			: success(true)
		{
			perfect2 = 1;
		}

		template <typename T, typename U>
		perfect_forwarding_test(T&& /*imperfect1*/, U&& /*imperfect2*/)
			: success(false)
		{}
	};



	class non_copyable_type
	{
	private:
		int i;
	public:
		non_copyable_type(const non_copyable_type&); // non construction-copyable
		non_copyable_type& operator=(const non_copyable_type&); // non copyable
		non_copyable_type(int a) : i(a) {}
	};

#endif


struct small_struct
{
	double *empty_field_1;
	double unused_number;
	unsigned int empty_field2;
	double *empty_field_3;
	int number;
	unsigned int empty_field4;

	int operator * () const { return number; };
	bool operator == (const small_struct &source) const { return source.number == number; };
	bool operator != (const small_struct &source) const { return source.number != number; };
	bool operator > (const small_struct &source) const { return number > source.number; };
	bool operator < (const small_struct &source) const { return number < source.number; };
	bool operator >= (const small_struct &source) const { return number >= source.number; };
	bool operator <= (const small_struct &source) const { return number <= source.number; };

	small_struct(const unsigned int num): number(num) {};
};



struct exceptions_test
{
	int num;

	exceptions_test(const int num2): num(num2) {};

	exceptions_test(const exceptions_test &et)
	{
		if (et.num == 4) throw 123;
		num = et.num;
	};
};





int global_counter = 0;

struct small_struct_non_trivial
{
	double *empty_field_1;
	double unused_number;
	unsigned int empty_field2;
	double *empty_field_3;
	int number;
	unsigned int empty_field4;

	small_struct_non_trivial(const int num) : number(num) {}
	~small_struct_non_trivial() { ++global_counter; }
};





int main()
{
	freopen("error.log","w", stderr); // For catching assertion failure info when run outside of a command line prompt

	using namespace std;
	using namespace plf;


	#ifdef PLF_TYPE_TRAITS_SUPPORT
	{
		title2("Trivially copyable iterators tests");

		{
			typedef plf::colony<small_struct_non_trivial> colony_type;
			failpass("Trivially-copyable iterators for non-trivial types", std::is_trivially_copyable<colony_type::iterator>::value);
			failpass("Trivially-copyable const_iterators for non-trivial types", std::is_trivially_copyable<colony_type::const_iterator>::value);

			failpass("Trivially-copyable reverse_iterators for non-trivial types ", std::is_trivially_copyable<colony_type::reverse_iterator>::value);
			failpass("Trivially-copyable reverse_const_iterators for non-trivial types", std::is_trivially_copyable<colony_type::const_reverse_iterator>::value);
		}

		{
			typedef plf::colony<int> colony_type;
			failpass("Trivially-copyable iterators for trivial types", std::is_trivially_copyable<colony_type::iterator>::value);
			failpass("Trivially-copyable const_iterators for trivial types", std::is_trivially_copyable<colony_type::const_iterator>::value);

			failpass("Trivially-copyable reverse_iterators for trivial types", std::is_trivially_copyable<colony_type::reverse_iterator>::value);
			failpass("Trivially-copyable reverse_const_iterators for trivial types", std::is_trivially_copyable<colony_type::const_reverse_iterator>::value);
		}
	}
	#endif


	printf("\nExceptions will be thrown (and caught) by this suite in order to test sub-function exception recovery. Tests loop many times in order to catch edge cases (via random number generation). Press enter to continue.\n");
	getchar();


	for (unsigned int looper = 0; looper != 100; ++looper)
	{
		{
			title1("Colony");
			title2("Test Basics");

			colony<int *> p_colony;

			failpass("Colony empty", p_colony.empty());

			int ten = 10;
			p_colony.insert(&ten);

			failpass("Colony not-empty", !p_colony.empty());

			title2("Iterator tests");

			failpass("Begin() working", **p_colony.begin() == 10);
			failpass("End() working", p_colony.begin() != p_colony.end());


			p_colony.clear();

			failpass("Begin = End after clear", p_colony.begin() == p_colony.end());

			int twenty = 20;

			for (unsigned int temp = 0; temp != 200; ++temp)
			{
				p_colony.insert(&ten);
				p_colony.insert(&twenty);
			}

			int total = 0, numtotal = 0;

			for(colony<int *>::iterator the_iterator = p_colony.begin(); the_iterator != p_colony.end(); ++the_iterator)
			{
				++total;
				numtotal += **the_iterator;
			}

			failpass("Iteration count test", total == 400);
			failpass("Iterator access test", numtotal == 6000);

			colony<int *>::iterator plus_twenty = p_colony.begin();
			advance(plus_twenty, 20);
			colony<int *>::iterator plus_two_hundred = p_colony.begin();
			advance(plus_two_hundred, 200);

			failpass("Iterator + distance test", distance(p_colony.begin(), plus_twenty) == 20);
			failpass("Iterator - distance test", distance(plus_two_hundred, p_colony.begin()) == -200);

			{
				colony<int> d_colony(1000, 1, plf::limits(20, 20));

				for (colony<int>::iterator current = d_colony.begin(), end = d_colony.end(); current != end;)
				{
					if ((rand() & 7) == 0)
					{
						current = d_colony.erase(current);
					}
					else
					{
						++current;
					}
				}

				int d_size = static_cast<int>(d_colony.size());

				for (int counter = 0; counter != 10000; ++counter)
				{
					const int dist1 = rand() % (d_size - 2), dist2 = rand() % ((d_size - 2) - dist1);
					colony<int>::iterator first = d_colony.begin(), last;
					advance(first, dist1);
					last = first;
					advance(last, dist2);

					const int dist = static_cast<int>(distance(first, last));

					if (dist != dist2)
					{
						printf("positive distance overload fuzz-test failed, real distance = %d, reported distance = %d, counter = %d, suite loop = %d", dist2, dist, counter, looper);
						getchar();
						abort();
					}
				}

				failpass("Positive distance overload fuzz-test", true);


				for (int counter = 0; counter != 10000; ++counter)
				{
					const int dist1 = rand() % (d_size - 2), dist2 = rand() % ((d_size - 2) - dist1);
					colony<int>::iterator first = d_colony.end(), last;
					advance(first, -dist1);
					last = first;
					advance(last, -dist2);

					const int dist = static_cast<int>(distance(last, first));

					if (dist != dist2)
					{
						printf("Negative distance overload fuzz-test failed, real distance = %d, reported distance = %d, counter = %d, suite loop = %d", dist2, dist, counter, looper);
						getchar();
						abort();
					}
				}

				failpass("Negative distance overload fuzz-test", true);



				for (colony<int>::iterator current = d_colony.begin(), end = d_colony.end(); current!= end;)
				{
					if ((rand() & 3) == 0)
					{
						current = d_colony.erase(current);
					}
					else
					{
						++current;
					}
				}

				d_size = static_cast<int>(d_colony.size());


				for (int counter = 0; counter != 10000; ++counter)
				{
					const int dist1 = rand() % (d_size - 2), dist2 = rand() % (d_size - 2);
					colony<int>::iterator first = d_colony.begin(), last = d_colony.begin();
					advance(first, dist1);
					advance(last, dist2);
					if (last > d_colony.end()) last = d_colony.end();

					const int dist = static_cast<int>(distance(first, last));

					if (dist != dist2 - dist1)
					{
						printf("positive/negative distance overload fuzz-test failed, real distance = %d, reported distance = %d, counter = %d, suite loop = %d", dist2, dist, counter, looper);
						getchar();
						abort();
					}
				}

				failpass("Positive/negative distance overload fuzz-test", true);
			}



			#ifdef PLF_CPP20_SUPPORT
				colony<int *>::const_iterator plus_two_hundred_c = plus_two_hundred;
				colony<int *> colony_copy(plus_twenty, plus_two_hundred_c);

				total = 0;

				for(colony<int *>::iterator the_iterator = colony_copy.begin(); the_iterator != colony_copy.end(); ++the_iterator)
				{
					++total;
				}

				failpass("Range-constructor with differing iterator types test", total == 180);
			#endif


			colony<int *>::iterator next_iterator = next(p_colony.begin(), 5);
			colony<int *>::const_iterator prev_iterator = prev(p_colony.cend(), 300);

			failpass("Iterator next test", distance(p_colony.begin(), next_iterator) == 5);
			failpass("Const iterator prev test", distance(p_colony.cend(), prev_iterator) == -300);

			#if defined(__cplusplus) && __cplusplus >= 201402L
				colony<int *>::iterator prev_iterator2 = prev(p_colony.end(), 300);
				failpass("Iterator/Const iterator equality operator test", prev_iterator == prev_iterator2);
			#endif

			prev_iterator = p_colony.begin();
			advance(prev_iterator, 5);
			failpass("Iterator/Const iterator equality operator test 2", prev_iterator == next_iterator);

			colony<int *> p_colony2;
			p_colony2 = p_colony;
			colony<int *> p_colony3(p_colony);
			colony<int *> p_colony4(p_colony2, p_colony2.get_allocator());

			colony<int *>::iterator it1 = p_colony.begin();
			colony<int *>::const_iterator cit(it1);

			failpass("Copy test", p_colony2.size() == 400);
			failpass("Copy construct test", p_colony3.size() == 400);
			failpass("Allocator-extended copy construct test", p_colony4.size() == 400);


			failpass("Equality operator test", p_colony == p_colony2);
			failpass("Equality operator test 2", p_colony2 == p_colony3);

			p_colony2.insert(&ten);

			failpass("Inequality operator test", p_colony2 != p_colony3);

			#ifdef PLF_CPP20_SUPPORT
				failpass("Spaceship operator test", (p_colony2 <=> p_colony3) != 0);
			#endif

			numtotal = 0;
			total = 0;

			for (colony<int *>::reverse_iterator the_iterator = p_colony.rbegin(); the_iterator != p_colony.rend(); ++the_iterator)
			{
				++total;
				numtotal += **the_iterator;
			}

			failpass("Reverse iteration count test", total == 400);
			failpass("Reverse iterator access test", numtotal == 6000);

			{
				colony<int> numbers;
				numbers.insert(1);
				colony<int>::reverse_iterator rit = numbers.rend(), rit2 = prev(rit, 1), rit3 = numbers.rbegin();
				failpass("Reverse iteration prev/rbegin/rend test", rit2 == rit3);
			}


			colony<int *>::reverse_iterator r_iterator = p_colony.rbegin();
			advance(r_iterator, 50);

			failpass("Reverse iterator advance and distance test", distance(p_colony.rbegin(), r_iterator) == 50);

			colony<int *>::reverse_iterator r_iterator2 = next(r_iterator, 2);

			failpass("Reverse iterator next and distance test", distance(p_colony.rbegin(), r_iterator2) == 52);

			numtotal = 0;
			total = 0;

			for(colony<int *>::iterator the_iterator = p_colony.begin(); the_iterator < p_colony.end(); advance(the_iterator, 2))
			{
				++total;
				numtotal += **the_iterator;
			}

			failpass("Multiple iteration test", total == 200);
			failpass("Multiple iteration access test", numtotal == 2000);

			numtotal = 0;
			total = 0;

			for(colony<int *>::const_iterator the_iterator = p_colony.cbegin(); the_iterator != p_colony.cend(); ++the_iterator)
			{
				++total;
				numtotal += **the_iterator;
			}

			failpass("Const_iterator test", total == 400);
			failpass("Const_iterator access test", numtotal == 6000);


			numtotal = 0;
			total = 0;

			for(colony<int *>::const_reverse_iterator the_iterator = --colony<int *>::const_reverse_iterator(p_colony.crend()); the_iterator != colony<int *>::const_reverse_iterator(p_colony.crbegin()); --the_iterator)
			{
				++total;
				numtotal += **the_iterator;
			}

			failpass("Const_reverse_iterator -- test", total == 399);
			failpass("Const_reverse_iterator -- access test", numtotal == 5980);

			total = 0;

			for(colony<int *>::iterator the_iterator = ++colony<int *>::iterator(p_colony.begin()); the_iterator < p_colony.end(); ++the_iterator)
			{
				++total;
				the_iterator = p_colony.erase(the_iterator);
			}

			failpass("Partial erase iteration test", total == 200);
			failpass("Post-erase size test", p_colony.size() == 200);

			#ifdef PLF_INITIALIZER_LIST_SUPPORT
				{
					colony<int> trim_colony(2000, 10, {200, 200});
					trim_colony.reserve(4000);
					trim_colony.trim_capacity(3000);
					failpass("trim_capacity(n) test", trim_colony.capacity() == 3000);
				}
			#endif

			const unsigned int temp_capacity = static_cast<unsigned int>(p_colony.capacity());
			p_colony.shrink_to_fit();
			failpass("Shrink_to_fit test", p_colony.capacity() < temp_capacity);
			failpass("Shrink_to_fit test 2", p_colony.capacity() == 200);

			total = 0;

			for(colony<int *>::reverse_iterator the_iterator = p_colony.rbegin(); the_iterator != p_colony.rend();)
			{
				colony<int *>::iterator it = the_iterator.base();
				the_iterator = colony<int *>::reverse_iterator(p_colony.erase(--it));
				++total;
			}

			failpass("Full erase reverse iteration test", total == 200);
			failpass("Post-erase size test", p_colony.size() == 0);

			for (unsigned int temp = 0; temp != 200; ++temp)
			{
				p_colony.insert(&ten);
				p_colony.insert(&twenty);
			}

			total = 0;

			for(colony<int *>::iterator the_iterator = --colony<int *>::iterator(p_colony.end()); the_iterator != p_colony.begin(); --the_iterator)
			{
				++total;
			}

			failpass("Negative iteration test", total == 399);


			total = 0;

			for(colony<int *>::iterator the_iterator = --(colony<int *>::iterator(p_colony.end())); the_iterator > p_colony.begin(); advance(the_iterator, -2))
			{
				++total;
			}

			failpass("Negative multiple iteration test", total == 200);

			#ifdef PLF_MOVE_SEMANTICS_SUPPORT
				p_colony2 = std::move(p_colony);
				failpass("Move test", p_colony2.size() == 400);

				p_colony.insert(&ten);

				failpass("Insert to post-moved-colony test", p_colony.size() == 1);

				colony<int *> p_colony5(p_colony2);
				colony<int *> p_colony6(std::move(p_colony5), p_colony2.get_allocator());

				failpass("Allocator-extended move construct test", p_colony6.size() == 400);
			#else
				p_colony2 = p_colony;
			#endif

			p_colony3 = p_colony2;

			failpass("Copy test 2", p_colony3.size() == 400);

			p_colony2.insert(&ten);

			p_colony2.swap(p_colony3);

			failpass("Swap test", p_colony2.size() == p_colony3.size() - 1);

			std::swap(p_colony2, p_colony3);

			failpass("Swap test 2", p_colony3.size() == p_colony2.size() - 1);

			failpass("max_size() test", p_colony2.max_size() > p_colony2.size());

		}


		{
			title2("Iterator comparison tests");

			colony<int> i_colony;

			for (int temp = 0; temp != 10; ++temp)
			{
				i_colony.insert(temp);
			}

			colony<int>::iterator it1 = i_colony.begin(), it2 = i_colony.begin();

			++it2;
			++it2;
			++it2;

			failpass("Iterator ++ test", *it2 == 3);

			failpass("Iterator > test", it2 > it1);

			failpass("Iterator >= test", it2 >= it1);

			failpass("Iterator < test", it1 < it2);

			failpass("Iterator <= test", it1 <= it2);

			failpass("Iterator != test", it2 != it1);

			#ifdef PLF_CPP20_SUPPORT
				failpass("Iterator <=> test 1", (it2 <=> it1) == std::strong_ordering::greater);

				failpass("Iterator <=> test 2", (it1 <=> it2) == std::strong_ordering::less);

				it1 = it2;

				failpass("Iterator <=> test 3", (it1 <=> it2) == std::strong_ordering::equal);
			#endif
		}


		{
			title2("Insert and Erase tests");

			colony<int> i_colony;

			for (int temp = 0; temp != 500000; ++temp)
			{
				i_colony.insert(temp);
			}


			failpass("Size after insert test", i_colony.size() == 500000);


			colony<int>::iterator found_item = std::find(i_colony.begin(), i_colony.end(), 5000);;

			failpass("std::find iterator test", *found_item == 5000);


			colony<int>::reverse_iterator found_item2 = std::find(i_colony.rbegin(), i_colony.rend(), 5000);;

			failpass("std::find reverse_iterator test", *found_item2 == 5000);


			for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end(); ++the_iterator)
			{
				the_iterator = i_colony.erase(the_iterator);
			}

			failpass("Erase alternating test", i_colony.size() == 250000);

			do
			{
				for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end();)
				{
					if ((rand() & 7) == 0)
					{
						the_iterator = i_colony.erase(the_iterator);
					}
					else
					{
						++the_iterator;
					}
				}
			} while (!i_colony.empty());

			failpass("Erase randomly till-empty test", i_colony.size() == 0);


			i_colony.reset();
			i_colony.reshape(plf::limits(100, i_colony.block_capacity_limits().max));

			i_colony.insert(30000, 1); // fill-insert 30000 elements

			failpass("Size after reinitialize + fill-insert test", i_colony.size() == 30000);

			unsigned short count2 = 0;

			do
			{
				for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end();)
				{
					if ((rand() & 7) == 0)
					{
						the_iterator = i_colony.erase(the_iterator);
						++count2;
					}
					else
					{
						++the_iterator;
					}
				}

			} while (count2 < 15000);

			failpass("Erase randomly till half-empty test", i_colony.size() == 30000u - count2);

			i_colony.insert(count2, 1);

			failpass("Size after reinsert test", i_colony.size() == 30000);




			unsigned int sum = 0;

			for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end();)
			{
				if (++sum == 3)
				{
					sum = 0;
					the_iterator = i_colony.erase(the_iterator);
				}
				else
				{
					i_colony.insert(1);
					++the_iterator;
				}
			}

			failpass("Alternating insert/erase test", i_colony.size() == 45001);


			do
			{
				for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end();)
				{
					if ((rand() & 3) == 0)
					{
						++the_iterator;
						i_colony.insert(1);
					}
					else
					{
						the_iterator = i_colony.erase(the_iterator);
					}
				}
			} while (!i_colony.empty());;

			failpass("Random insert/erase till empty test", i_colony.size() == 0);


			i_colony.insert(500000, 10);

			failpass("Insert post-erase test", i_colony.size() == 500000);
			colony<int>::iterator it2 = i_colony.begin();
			advance(it2, 250000);


			for (; it2 != i_colony.end();)
			{
				it2 = i_colony.erase(it2);
			}

			failpass("Large multi-increment iterator test", i_colony.size() == 250000);

			i_colony.insert(250000, 10);

			colony<int>::iterator end_iterator = i_colony.end();
			colony<int>::iterator end_iterator2 = i_colony.end();
			advance(end_iterator, -250000);

			for (unsigned int count = 0; count != 250000; ++count, --end_iterator2){}
			failpass("Large multi-decrement iterator test 1", end_iterator == end_iterator2);


			for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != end_iterator;)
			{
				the_iterator = i_colony.erase(the_iterator);
			}

			failpass("Large multi-decrement iterator test", i_colony.size() == 250000);



			i_colony.insert(250000, 10);
			int total = 0;

			for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end(); ++the_iterator)
			{
				total += *the_iterator;
			}

			failpass("Re-insert post-heavy-erasure test", total == 5000000);


			end_iterator = i_colony.end();
			advance(end_iterator, -50001);
			colony<int>::iterator begin_iterator = i_colony.begin();
			advance(begin_iterator, 300000);

			for (colony<int>::iterator the_iterator = begin_iterator; the_iterator != end_iterator;)
			{
				the_iterator = i_colony.erase(the_iterator);
			}

			failpass("Non-end decrement + erase test", i_colony.size() == 350001);


			i_colony.insert(100000, 10);

			begin_iterator = i_colony.begin();
			advance(begin_iterator, 300001);


			for (colony<int>::iterator the_iterator = begin_iterator; the_iterator != i_colony.end();)
			{
				the_iterator = i_colony.erase(the_iterator);
			}

			failpass("Non-beginning increment + erase test", i_colony.size() == 300001);

			colony<int>::iterator temp_iterator = i_colony.begin();
			advance(temp_iterator, 20); // Advance test 1

			unsigned int index = static_cast<unsigned int>(distance(i_colony.begin(), temp_iterator));
			failpass("Advance + iterator-to-index test", index == 20);

			i_colony.erase(temp_iterator);

			failpass("is_active test 1", i_colony.is_active(temp_iterator) == false);
			failpass("get_iterator test 1", i_colony.get_iterator(&(*temp_iterator)) == i_colony.end());

			temp_iterator = i_colony.begin(); // Check edge-case with advance when erasures present in initial group

			failpass("is_active test 2", i_colony.is_active(temp_iterator) == true);
			failpass("get_iterator test 2", i_colony.get_iterator(&(*temp_iterator)) == temp_iterator);

			advance(temp_iterator, 500);

			index = static_cast<unsigned int>(distance(i_colony.begin(), temp_iterator));

			failpass("Advance + iterator-to-index test", index == 500);

			colony<int>::iterator temp2 = i_colony.get_iterator(&(*temp_iterator));
			colony<int>::const_iterator temp3 = i_colony.get_iterator(const_cast<colony<int>::const_pointer>(&(*temp_iterator)));

			failpass("Pointer-to-iterator test", temp2 != i_colony.end());
			failpass("Const_pointer-to-const_iterator test", temp3 != i_colony.end());

			temp2 = i_colony.begin();
			advance(temp2, 500);

			failpass("Index-to-iterator test", temp2 == temp_iterator);


			for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end();)
			{
				the_iterator = i_colony.erase(the_iterator);
			}

			failpass("Total erase test", i_colony.empty());

			// Test get_iterator etc on an empty colony with empty memory blocks retained:
			failpass("is_active test 3", i_colony.is_active(temp_iterator) == false);
  			failpass("get_iterator test 3", i_colony.get_iterator(&(*temp_iterator)) == i_colony.end());

			i_colony.trim_capacity();

			// Test get_iterator etc on a colony with no blocks:
			failpass("is_active test 4", i_colony.is_active(temp_iterator) == false);
  			failpass("get_iterator test 4", i_colony.get_iterator(&(*temp_iterator)) == i_colony.end());

			i_colony.reshape(plf::limits(3, i_colony.block_capacity_limits().max));

			const unsigned int temp_capacity2 = static_cast<unsigned int>(i_colony.capacity());
			i_colony.reserve(100000);
			failpass("Post-reset reserve test", temp_capacity2 != i_colony.capacity());

			i_colony.insert(110000, 1);

			failpass("Post-reserve insert test", i_colony.size() == 110000);

			unsigned int count = 110000;

			for (unsigned int loop1 = 0; loop1 != 50000; ++loop1)
			{
				for (unsigned int loop = 0; loop != 10; ++loop)
				{
					if ((rand() & 7) == 0)
					{
						i_colony.insert(1);
						++count;
					}
				}


				for (colony<int>::iterator the_iterator = i_colony.begin(); the_iterator != i_colony.end();)
				{
					if ((rand() & 7) == 0)
					{
						the_iterator = i_colony.erase(the_iterator);
						--count;
					}
					else
					{
						++the_iterator;
					}
				}
			}

			failpass("Multiple sequential small insert/erase commands test", count == i_colony.size());
		}


		{
			title2("Range-erase tests");

			colony<int> i_colony;

			int counter = 0;

			for (; counter != 1000; ++counter)
			{
				i_colony.insert(counter);
			}


			colony<int>::iterator it1 = i_colony.begin(), it2 = i_colony.begin();

			advance(it1, 500);
			advance(it2, 800);

			i_colony.erase(it1, it2);

			counter = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				++counter;
			}

			failpass("Simple range-erase test 1", counter == 700 && i_colony.size() == 700);


			it1 = it2 = i_colony.begin();

			advance(it1, 400);
			advance(it2, 500); // This should put it2 past the point of previous erasures

			i_colony.erase(it1, it2);

			counter = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				++counter;
			}

			failpass("Simple range-erase test 2", counter == 600 && i_colony.size() == 600);



			it2 = it1 = i_colony.begin();

			advance(it1, 4);
			advance(it2, 9); // This should put it2 past the point of previous erasures

			i_colony.erase(it1, it2);

			counter = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				++counter;
			}

			failpass("Simple range-erase test 3", counter == 595 && i_colony.size() == 595);




			it2 = it1 = i_colony.begin();

			advance(it2, 50);

			i_colony.erase(it1, it2);

			counter = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				++counter;
			}

			failpass("Range-erase from begin()", counter == 545 && i_colony.size() == 545);




			it1 = i_colony.begin();
			it2 = i_colony.end();

			advance(it1, 345); // Test erasing and validity when it removes the final group in colony
			i_colony.erase(it1, it2);

			counter = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				++counter;
			}

			failpass("Range-erase to end()", counter == 345 && i_colony.size() == 345);



			i_colony.clear();

			for (counter = 0; counter != 3000; ++counter)
			{
				i_colony.insert(counter);
			}

			for (colony<int>::iterator it = i_colony.begin(); it < i_colony.end(); ++it)
			{
				it = i_colony.erase(it);
			}

			it2 = it1 = i_colony.begin();

			advance(it1, 4);
			advance(it2, 600);
			i_colony.erase(it1, it2);

			counter = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				++counter;
			}

			failpass("Range-erase with colony already half-erased, alternating erasures", counter == 904 && i_colony.size() == 904);



			i_colony.clear();

			for (counter = 0; counter != 3000; ++counter)
			{
				i_colony.insert(counter);
			}

			for (colony<int>::iterator it = i_colony.begin(); it < i_colony.end(); ++it)
			{
				if ((rand() & 1) == 0)
				{
					it = i_colony.erase(it);
				}
			}

			if (i_colony.size() < 400)
			{
				for (counter = 0; counter != 400; ++counter)
				{
					i_colony.insert(counter);
				}
			}

			it1 = i_colony.begin();
			it2 = i_colony.end();

			advance(it1, 400);
			i_colony.erase(it1, it2);

			counter = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				++counter;
			}

			failpass("Range-erase with colony already third-erased, randomized erasures", counter == 400 && i_colony.size() == 400);



			unsigned int size, range1, range2, internal_loop_counter;

			for (unsigned int loop_counter = 0; loop_counter != 50; ++loop_counter)
			{
				i_colony.clear();

				for (counter = 0; counter != 1000; ++counter)
				{
					i_colony.insert(counter);
				}

				internal_loop_counter = 0;

				while (!i_colony.empty())
				{
					it2 = it1 = i_colony.begin();

					size = static_cast<unsigned int>(i_colony.size());
					range1 = rand() % size;
					range2 = range1 + 1 + (rand() % (size - range1));
					advance(it1, static_cast<int>(range1));
					advance(it2, static_cast<int>(range2));

					if (it2 > i_colony.end()) it2 = i_colony.end();

					i_colony.erase(it1, it2);

					counter = 0;

					for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
					{
						++counter;
					}

					if (i_colony.size() != static_cast<unsigned int>(counter))
					{
						printf("Fuzz-test range-erase randomly Fail: loop counter: %u, internal_loop_counter: %u.\n", loop_counter, internal_loop_counter);
						getchar();
						abort();
					}

					if (i_colony.size() != i_colony.group_size_sum())
					{
						printf("Fuzz-test range-erase randomly Fail - group_size_sum failure: loop counter: %u, internal_loop_counter: %u, size: %u, group_size_sum: %u.\n", loop_counter, internal_loop_counter, static_cast<unsigned int>(i_colony.size()), static_cast<unsigned int>(i_colony.group_size_sum()));
						getchar();
						abort();
					}

					if (i_colony.size() > 2)
					{ // Test to make sure our stored erased_locations are valid
						i_colony.insert(1);
						i_colony.insert(10);
					}

					++internal_loop_counter;
				}
			}

			failpass("Fuzz-test range-erase randomly until empty", i_colony.size() == 0);



			for (unsigned int loop_counter = 0; loop_counter != 50; ++loop_counter)
			{
				i_colony.clear();
				internal_loop_counter = 0;

				i_colony.insert(10000, 10);

				while (!i_colony.empty())
				{
					it2 = it1 = i_colony.begin();

					size = static_cast<unsigned int>(i_colony.size());
					range1 = rand() % size;
					range2 = range1 + 1 + (rand() % (size - range1));
					advance(it1, static_cast<int>(range1));
					advance(it2, static_cast<int>(range2));
					if (it2 > i_colony.end()) it2 = i_colony.end();

					i_colony.erase(it1, it2);

					counter = 0;

					for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
					{
						++counter;
					}

					if (i_colony.size() != i_colony.group_size_sum())
					{
						printf("Fuzz-test range-erase + fill-insert randomly Fails during erase - group_size_sum failure: loop counter: %u, internal_loop_counter: %u, size: %u, group_size_sum: %u.\n", loop_counter, internal_loop_counter, static_cast<unsigned int>(i_colony.size()), static_cast<unsigned int>(i_colony.group_size_sum()));
						getchar();
						abort();
					}

					if (i_colony.size() != static_cast<unsigned int>(counter))
					{
						printf("Fuzz-test range-erase + fill-insert randomly Fails during erase: loop counter: %u, internal_loop_counter: %u.\n", loop_counter, internal_loop_counter);
						getchar();
						abort();
					}

					if (i_colony.size() > 100)
					{ // Test to make sure our stored erased_locations are valid & fill-insert is functioning properly in these scenarios
						const unsigned int extra_size = rand() & 127;
						i_colony.insert(extra_size, 5);

						if (i_colony.size() != i_colony.group_size_sum())
						{
							printf("Fuzz-test range-erase + fill-insert randomly Fails during insert - group_size_sum failure: loop counter: %u, internal_loop_counter: %u, size: %u, group_size_sum: %u.\n", loop_counter, internal_loop_counter, static_cast<unsigned int>(i_colony.size()), static_cast<unsigned int>(i_colony.group_size_sum()));
							getchar();
							abort();
						}

						if (i_colony.size() != static_cast<unsigned int>(counter) + extra_size)
						{
							printf("Fuzz-test range-erase + fill-insert randomly Fails during fill-insert: loop counter: %u, internal_loop_counter: %u.\n", loop_counter, internal_loop_counter);
							getchar();
							abort();
						}

						counter = 0;

						for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
						{
							++counter;
						}

						if (i_colony.size() != static_cast<unsigned int>(counter))
						{
							printf("Fuzz-test range-erase + fill-insert randomly Fails during counter-test fill-insert: loop counter: %u, internal_loop_counter: %u.\n", loop_counter, internal_loop_counter);
							getchar();
							abort();
						}
					}

					++internal_loop_counter;
				}
			}

			failpass("Fuzz-test range-erase + fill-insert randomly until empty", i_colony.size() == 0);

			i_colony.erase(i_colony.begin(), i_colony.end());

			failpass("Range-erase when colony is empty test (crash test)", i_colony.size() == 0);

			i_colony.insert(10, 1);

			i_colony.erase(i_colony.begin(), i_colony.begin());

			failpass("Range-erase when range is empty test (crash test)", i_colony.size() == 10);
		}



		{
			title1("Non-trivial type tests");

			colony<small_struct_non_trivial> ss_nt;
			colony<small_struct_non_trivial>::iterator ss_it1, ss_it2;

			small_struct_non_trivial ss(5);

			unsigned int size, range1 = 0, range2 = 0, internal_loop_counter;
			int counter;

			ss_nt.insert(10000, ss);

			failpass("Non-trivial type insert test", ss_nt.size() == 10000);


			for (colony<small_struct_non_trivial>::iterator ss_it = ss_nt.begin(); ss_it != ss_nt.end(); ++ss_it)
			{
				ss_it = ss_nt.erase(ss_it);
				++range1;
			}

			failpass("Non-trivial type erase half of all elements", ss_nt.size() == 5000);


			colony<small_struct_non_trivial> ss_nt2;
			ss_nt2.insert(5, ss);

			ss_nt = ss_nt2;

			failpass("Non-trivial type range-assign", ss_nt.size() == 5);



			for (unsigned int loop_counter = 0; loop_counter != 50; ++loop_counter)
			{
				ss_nt.clear();

				for (counter = 0; counter != 1000; ++counter)
				{
					ss_nt.insert(counter);
				}

				internal_loop_counter = 0;

				while (!ss_nt.empty())
				{
					ss_it2 = ss_it1 = ss_nt.begin();

					size = static_cast<unsigned int>(ss_nt.size());
					range1 = rand() % size;
					range2 = range1 + 1 + (rand() % (size - range1));
					advance(ss_it1, static_cast<int>(range1));
					advance(ss_it2, static_cast<int>(range2));
					if (ss_it2 > ss_nt.end()) ss_it2 = ss_nt.end();

					ss_nt.erase(ss_it1, ss_it2);

					counter = 0;

					for (colony<small_struct_non_trivial>::iterator ss_it = ss_nt.begin(); ss_it != ss_nt.end(); ++ss_it)
					{
						++counter;
					}

					if (ss_nt.size() != static_cast<unsigned int>(counter))
					{
						printf("Fuzz-test range-erase randomly Fail: loop counter: %u, internal_loop_counter: %u.\n", loop_counter, internal_loop_counter);
						getchar();
						abort();
					}

					if (ss_nt.size() != ss_nt.group_size_sum())
					{
						printf("Fuzz-test range-erase randomly Fail - group_size_sum failure: loop counter: %u, internal_loop_counter: %u, size: %u, group_size_sum: %u.\n", loop_counter, internal_loop_counter, static_cast<unsigned int>(ss_nt.size()), static_cast<unsigned int>(ss_nt.group_size_sum()));
						getchar();
						abort();
					}

					if (ss_nt.size() > 2)
					{ // Test to make sure our stored erased_locations are valid
						ss_nt.insert(1);
						ss_nt.insert(10);
					}

					++internal_loop_counter;
				}
			}

			failpass("Non-trivial type fuzz-test range-erase randomly until empty", ss_nt.size() == 0);
		}




		{
			title2("Sort tests");

			colony<int> i_colony;

			i_colony.reserve(50000);

			for (unsigned int temp = 0; temp != 50000; ++temp)
			{
				i_colony.insert(rand() & 65535);
			}

			i_colony.sort();

			bool sorted = true;
			int previous = 0;

			for (colony<int>::iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				if (previous > *current)
				{
					sorted = false;
					break;
				}

				previous = *current;
			}

			failpass("Less-than sort test", sorted);

			i_colony.unique();


			bool unique = true;
			previous = -1;

			for (colony<int>::iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				if (previous == *current)
				{
					unique = false;
					break;
				}

				previous = *current;
			}

			failpass("Unique test", unique);

			i_colony.sort(std::greater<int>());

			previous = 65536;

			for (colony<int>::iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				if (previous < *current)
				{
					sorted = false;
					break;
				}

				previous = *current;
			}

			failpass("Greater-than sort test", sorted);


			i_colony.clear();
			i_colony.trim_capacity();
			i_colony.reserve(512);

			for (unsigned int temp = 0; temp != 200; ++temp)
			{
				i_colony.insert(rand() & 65535);
			}

			i_colony.sort();

			previous = 0;

			for (colony<int>::iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				if (previous > *current)
				{
					sorted = false;
					break;
				}

				previous = *current;
			}

			failpass("Less-than sort with allocation optimization 1 test - small types", sorted);


			colony<small_struct> ss_colony;
			ss_colony.reserve(512);

			for (unsigned int temp = 0; temp != 200; ++temp)
			{
				ss_colony.insert(rand() & 65535);
			}

			ss_colony.sort();

			small_struct ss_previous(0);

			for (colony<small_struct>::iterator current = ss_colony.begin(); current != ss_colony.end(); ++current)
			{
				if (ss_previous > *current)
				{
					sorted = false;
					break;
				}

				ss_previous = *current;
			}

			failpass("Less-than sort with allocation optimization 2 test - large type 1", sorted);


			ss_colony.clear();
			ss_colony.reshape(plf::limits(100, 100));
			ss_colony.reserve(200);

			for (unsigned int temp = 0; temp != 75; ++temp)
			{
				ss_colony.insert(rand() & 65535);
			}

			ss_colony.sort();

			ss_previous = 0;

			for (colony<small_struct>::iterator current = ss_colony.begin(); current != ss_colony.end(); ++current)
			{
				if (ss_previous > *current)
				{
					sorted = false;
					break;
				}

				ss_previous = *current;
			}

			failpass("Less-than sort with allocation optimization 3 test - large type 2", sorted);
		}




		{
			title2("Different insertion-style tests");

			#ifdef PLF_INITIALIZER_LIST_SUPPORT
				colony<int> i_colony({1, 2, 3});

				failpass("Initializer-list constructor test", i_colony.size() == 3);
			#else
				colony<int> i_colony(3, 1);
			#endif

			colony<int> i_colony2(i_colony.begin(), i_colony.end());

			failpass("Range constructor test", i_colony2.size() == 3);

			#ifdef PLF_CPP20_SUPPORT
				std::ranges::take_view<std::ranges::ref_view<plf::colony<int>>> rng = i_colony2 | std::ranges::views::take(2);

				colony<int> i_colony_range(plf::ranges::from_range, rng);

				failpass("Rangesv3 constructor test", i_colony_range.size() == 2);
			#endif

			colony<int> i_colony3(5000, 2, plf::limits(100, 250));

			failpass("Fill construction test", i_colony3.size() == 5000);

			i_colony2.insert(500000, 5);

			failpass("Fill insertion test", i_colony2.size() == 500003);

			std::vector<int> some_ints(500, 2);

			i_colony2.insert(some_ints.begin(), some_ints.end());

			failpass("Range insertion test", i_colony2.size() == 500503);

			#ifdef PLF_CPP20_SUPPORT
				i_colony2.insert(some_ints.begin(), some_ints.cend());

				failpass("Range insertion with differing iterators test", i_colony2.size() == 501003);
			#else
				i_colony2.insert(some_ints.begin(), some_ints.end());
			#endif

			#ifdef PLF_MOVE_SEMANTICS_SUPPORT
	 			i_colony3.insert(plf::make_move_iterator(i_colony2.begin()), plf::make_move_iterator(i_colony2.end()));

	 			failpass("Range move-insertion test", i_colony3.size() == 506003);

				i_colony3.clear();
			#endif

			i_colony2.reset();
			i_colony2.reserve(50000);
			i_colony2.insert(60000, 1);

			int total = 0;

			for (colony<int>::iterator it = i_colony2.begin(); it != i_colony2.end(); ++it)
			{
				total += *it;
			}

			failpass("Reserve + fill insert test", i_colony2.size() == 60000 && total == 60000);


			i_colony2.clear();
			i_colony2.reserve(5000);
			i_colony2.insert(60, 1);

			total = 0;

			for (colony<int>::iterator it = i_colony2.begin(); it != i_colony2.end(); ++it)
			{
				total += *it;
			}

			failpass("Reserve + fill insert test 2", i_colony2.size() == 60 && total == 60);

			i_colony2.insert(6000, 1);

			total = 0;

			for (colony<int>::iterator it = i_colony2.begin(); it != i_colony2.end(); ++it)
			{
				total += *it;
			}

			failpass("Reserve + fill + fill test", i_colony2.size() == 6060 && total == 6060);

			i_colony2.reserve(18000);
			i_colony2.insert(6000, 1);

			total = 0;

			for (colony<int>::iterator it = i_colony2.begin(); it != i_colony2.end(); ++it)
			{
				total += *it;
			}

			failpass("Reserve + fill + fill + reserve + fill test", i_colony2.size() == 12060 && total == 12060);

			i_colony2.clear();
			i_colony2.insert(6000, 2);

			failpass("Clear + fill test", i_colony2.size() == 6000 && *(i_colony2.begin()) == 2);

			i_colony.insert(i_colony2.begin(), i_colony2.end());

			failpass("Range insert when not empty test", i_colony.size() == 6003);

			i_colony2.reset();
			i_colony2.reshape(plf::limits(50, 70));
			i_colony2.reserve(150);

			failpass("Improved reserve accuracy test", i_colony2.capacity() == 150);
		}


		{
			title2("Assign tests");

			colony<int> i_colony(50, 2);

			i_colony.assign(50, 1);

			int total = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				total += *it;
			}

			failpass("Equal capacity assign test", i_colony.size() == 50 && total == 50);


			i_colony.assign(10, 2);

			total = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				total += *it;
			}

			failpass("Lesser capacity assign test", i_colony.size() == 10 && total == 20);


			i_colony.assign(2000, 20);

			total = 0;

			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				total += *it;
			}

			failpass("Greater capacity assign test", i_colony.size() == 2000 && total == 40000);

			i_colony.clear();


			for (unsigned int internal_loop_counter = 0; internal_loop_counter != 10; ++internal_loop_counter)
			{
				const unsigned int capacity = rand() & 65535;
				i_colony.assign(capacity, 1);

				total = 0;

				for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
				{
					total += *it;
				}

				if (i_colony.size() != capacity)
				{
					printf("Fuzz-test assign capacity Fail: global loop counter: %u, internal loop counter: %u.\n", looper, internal_loop_counter);
					getchar();
					abort();
				}

				if (i_colony.size() != static_cast<unsigned int>(total))
				{
					printf("Fuzz-test assign sum Fail: global loop counter: %u, internal loop counter: %u.\n", looper, internal_loop_counter);
					getchar();
					abort();
				}
			}

			message("Fuzz-test assign passed.");


			i_colony.clear();

			std::vector<int> i_vector;

			for (int counter = 1; counter != 11; ++counter)
			{
				i_vector.push_back(counter);
			}

			i_colony.assign(i_vector.begin(), i_vector.end());

			colony<int>::iterator it = i_colony.begin();
			bool fail = false;

			for (int counter = 1; counter != 11; ++counter, ++it)
			{
				if (*it != counter)
				{
					fail = true;
					break;
				}
			}

			failpass("Range assign test", i_colony.size() == 10 && !fail);


			i_colony.clear();


			for (unsigned int internal_loop_counter = 0; internal_loop_counter != 10; ++internal_loop_counter)
			{
				const unsigned int capacity = rand() & 65535;
				i_vector.assign(capacity, 1);
				i_colony.assign(i_vector.begin(), i_vector.end());

				total = 0;

				for (colony<int>::iterator it3 = i_colony.begin(); it3 != i_colony.end(); ++it3)
				{
					total += *it3;
				}

				if (i_colony.size() != capacity)
				{
					printf("Fuzz-test range assign capacity Fail: global loop counter: %u, internal loop counter: %u.\n", looper, internal_loop_counter);
					getchar();
					abort();
				}

				if (i_colony.size() != static_cast<unsigned int>(total))
				{
					printf("Fuzz-test range assign sum Fail: global loop counter: %u, internal loop counter: %u.\n", looper, internal_loop_counter);
					getchar();
					abort();
				}
			}

			message("Fuzz-test range assign passed.");


			i_colony.clear();


			for (unsigned int internal_loop_counter = 0; internal_loop_counter != 10; ++internal_loop_counter)
			{
				const unsigned int capacity = rand() & 65535;
				i_vector.assign(capacity, 1);
				i_colony.assign(i_vector.begin(), i_vector.end());

				total = 0;
				unsigned int subtract = 0;

				for (colony<int>::iterator it3 = i_colony.begin(); it3 != i_colony.end(); )
				{
					if ((rand() & 7) == 0)
					{
						it3 = i_colony.erase(it3);
						++subtract;
					}
					else
					{
						total += *it3++;
					}
				}

				if (i_colony.size() != capacity - subtract)
				{
					printf("Fuzz-test range assign + erase capacity Fail: global loop counter: %u, internal loop counter: %u.\n", looper, internal_loop_counter);
					getchar();
					abort();
				}

				if (i_colony.size() != static_cast<unsigned int>(total))
				{
					printf("Fuzz-test range assign + erase sum Fail: global loop counter: %u, internal loop counter: %u.\n", looper, internal_loop_counter);
					getchar();
					abort();
				}
			}

			message("Fuzz-test range assign + erase passed.");


			i_colony.clear();

			#ifdef PLF_INITIALIZER_LIST_SUPPORT
				i_colony.assign({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
				it = i_colony.begin();

				for (int counter = 1; counter != 11; ++counter, ++it)
				{
					if (*it != counter)
					{
						fail = true;
						break;
					}
				}

				failpass("Initializer_list assign test", i_colony.size() == 10 && !fail);
			#endif
		}


		#ifdef PLF_VARIADICS_SUPPORT
			{
				title2("Perfect Forwarding tests");

				colony<perfect_forwarding_test> pf_colony;

				int lvalue = 0;
				int &lvalueref = lvalue;

				pf_colony.emplace(7, lvalueref);

				failpass("Perfect forwarding test", (*pf_colony.begin()).success);
				failpass("Perfect forwarding test 2", lvalueref == 1);
			}


			{
				title2("Basic emplace test");

				colony<small_struct> ss_colony;
				int total1 = 0, total2 = 0;

				for (int counter = 0; counter != 100; ++counter)
				{
					ss_colony.emplace(counter);
					total1 += counter;
				}

				for (colony<small_struct>::iterator it = ss_colony.begin(); it != ss_colony.end(); ++it)
				{
					total2 += it->number;
				}

				failpass("Basic emplace test", total1 == total2);
				failpass("Basic emplace test 2", ss_colony.size() == 100);
			}


			{
				title2("Non-copyable type test");

				colony<non_copyable_type> temp;

				temp.emplace(1);
				temp.emplace(2);

				failpass("Non-copyable size test", temp.size() == 2);
			}

		#endif


		{
			title2("Misc function tests");

			colony<int> colony1;
			colony1.reshape(plf::limits(50, 100));

			colony1.insert(27);

			failpass("Change_group_sizes min-size test", colony1.capacity() == 50);

			for (int counter = 0; counter != 100; ++counter)
			{
				colony1.insert(counter);
			}

			failpass("Change_group_sizes max-size test", colony1.capacity() == 200);

			colony1.clear();
			colony1.reshape(plf::limits(200, 255));

			colony1.insert(27);

			failpass("Reshape min-size test", colony1.capacity() == 200);

			plf::limits temp_limits = colony1.block_capacity_limits();

			failpass("get_block_capacity_limits test", temp_limits.min == 200 && temp_limits.max == 255);

			temp_limits = plf::colony<int>::block_capacity_hard_limits();

			failpass("get_block_capacity_limits test", temp_limits.min == 3 && temp_limits.max == 255);

			for (int counter = 0; counter != 3300; ++counter)
			{
				colony1.insert(counter);
			}

			failpass("Reshape max-size test", colony1.capacity() == 3460);

			colony1.reshape(plf::limits(150, 150));

			failpass("Reshape test 3", colony1.capacity() == 3450);

			colony1.reshape(plf::limits(200, 200));

			failpass("Reshape test 4", colony1.capacity() == 3400);

		}

		{
			title2("Splice tests");

			{
				colony<int> colony1, colony2;

				for (int number = 0; number != 20; ++number)
				{
					colony1.insert(number);
					colony2.insert(number + 20);
				}

				colony1.splice(colony2);

				int check_number = 0;
				bool fail = false;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					if (check_number++ != *current)
					{
						fail = true;
					}
				}

				failpass("Small splice test 1", fail == false);
			}


			{
				colony<int> colony1, colony2;

				for(int number = 0; number != 100; ++number)
				{
					colony1.insert(number);
					colony2.insert(number + 100);
				}

				colony1.splice(colony2);

				int check_number = 0;
				bool fail = false;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					if (check_number++ != *current)
					{
						fail = true;
					}
				}

				failpass("Small splice test 2", fail == false);
			}


			{
				colony<int> colony1, colony2;

				for(int number = 0; number != 100000; ++number)
				{
					colony1.insert(number);
					colony2.insert(number + 100000);
				}

				colony1.splice(colony2);

				int check_number = 0;
				bool fail = false;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					if (check_number++ != *current)
					{
						fail = true;
					}
				}

				failpass("Large splice test 1", fail == false);
			}


			{
				colony<int> colony1, colony2;

				for(int number = 0; number != 100; ++number)
				{
					colony1.insert(number);
					colony2.insert(number + 100);
				}


				for (colony<int>::iterator current = colony2.begin(); current != colony2.end();)
				{
					if ((rand() & 7) == 0)
					{
						current = colony2.erase(current);
					}
					else
					{
						++current;
					}
				}


				colony1.splice(colony2);

				int check_number = -1;
				bool fail = false;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					if (check_number >= *current)
					{
						fail = true;
					}

					check_number = *current;
				}

				failpass("Erase + splice test 1", fail == false);
			}



			{
				colony<int> colony1, colony2;

				for(int number = 0; number != 100; ++number)
				{
					colony1.insert(number);
					colony2.insert(number + 100);
				}



				for (colony<int>::iterator current = colony2.begin(); current != colony2.end();)
				{
					if ((rand() & 3) == 0)
					{
						current = colony2.erase(current);
					}
					else
					{
						++current;
					}
				}


				for (colony<int>::iterator current = colony1.begin(); current != colony1.end();)
				{
					if ((rand() & 1) == 0)
					{
						current = colony1.erase(current);
					}
					else
					{
						++current;
					}
				}


				colony1.splice(colony2);

				int check_number = -1;
				bool fail = false;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					if (check_number >= *current)
					{
						fail = true;
					}

					check_number = *current;
				}

				failpass("Erase + splice test 2", fail == false);
			}



			{
				colony<int> colony1, colony2;

				colony1.reshape(plf::limits(200, 200));
				colony2.reshape(plf::limits(200, 200));

				for(int number = 0; number != 100; ++number)
				{
					colony1.insert(number + 150);
				}


				for(int number = 0; number != 150; ++number)
				{
					colony2.insert(number);
				}


				colony1.splice(colony2);

				int check_number = -1;
				bool fail = false;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					if (check_number >= *current)
					{
						fail = true;
					}

					check_number = *current;
				}

				failpass("Unequal size splice test 1", fail == false);
			}



			{
				colony<int> colony1(plf::limits(200, 200)), colony2(plf::limits(200, 200));

				for(int number = 0; number != 100; ++number)
				{
					colony1.insert(100 - number);
				}


				for(int number = 0; number != 150; ++number)
				{
					colony2.insert(250 - number);
				}


				colony1.splice(colony2);

				int check_number = 255;
				bool fail = false;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					if (check_number < *current)
					{
						fail = true;
					}

					check_number = *current;
				}

				failpass("Unequal size splice test 2", fail == false);
			}



			{
				colony<int> colony1, colony2;

				for(int number = 0; number != 100000; ++number)
				{
					colony1.insert(number + 200000);
				}


				for(int number = 0; number != 200000; ++number)
				{
					colony2.insert(number);
				}

				for (colony<int>::iterator current = colony2.begin(); current != colony2.end();)
				{
					if ((rand() & 1) == 0)
					{
						current = colony2.erase(current);
					}
					else
					{
						++current;
					}
				}


				for (colony<int>::iterator current = colony1.begin(); current != colony1.end();)
				{
					if ((rand() & 1) == 0)
					{
						current = colony1.erase(current);
					}
					else
					{
						++current;
					}
				}


				colony1.erase(--(colony1.end()));
				colony2.erase(--(colony2.end()));

				const unsigned int total = static_cast<unsigned int>(colony1.size() + colony2.size());

				colony1.splice(colony2);

				unsigned int check_number = 0;

				for (colony<int>::iterator current = colony1.begin(); current != colony1.end(); ++current)
				{
					++check_number;
				}

				failpass("Large unequal size + erase splice test 1", check_number == total);


				do
				{
					for (colony<int>::iterator current = colony1.begin(); current != colony1.end();)
					{
						if ((rand() & 3) == 0)
						{
							current = colony1.erase(current);
						}
						else if ((rand() & 7) == 0)
						{
							colony1.insert(433);
							++current;
						}
						else
						{
							++current;
						}
					}

				} while (!colony1.empty());

				failpass("Post-splice insert-and-erase randomly till-empty test", colony1.size() == 0);
			}
		}

		{
			title2("erase_if tests");

			colony<int> i_colony(100, 100);

			i_colony.insert(100, 200);
			colony<int> i_colony2 = i_colony;

			std::erase(i_colony, 100);

			failpass("non-member erase test 1", i_colony.size() == 100);

			erase(i_colony2, 200);


			failpass("non-member erase test 2", i_colony2.size() == 100);


			i_colony.clear();

			for(int count = 0; count != 1000; ++count)
			{
				i_colony.insert((rand() & 1));
			}

			i_colony2 = i_colony;

			const int count0 = static_cast<int>(std::count(i_colony.begin(), i_colony.end(), 0));
			const int count1 = 1000 - count0;

			erase(i_colony, 0);
			failpass("random non-member erase test 1", static_cast<int>(i_colony.size()) == count1);

			erase(i_colony2, 1);
			failpass("random non-member erase test 2",  static_cast<int>(i_colony2.size()) == count0);


			i_colony.clear();

			for(int count = 0; count != 1000; ++count)
			{
				i_colony.insert(count);
			}

			#ifdef PLF_MOVE_SEMANTICS_SUPPORT // approximating checking for C++11 here
				erase_if(i_colony, std::bind(std::greater<int>(), std::placeholders::_1, 499));
			#else // C++03 or lower
				erase_if(i_colony, std::bind2nd(std::greater<int>(), 499));
			#endif

			failpass("erase_if test",	static_cast<int>(i_colony.size()) == 500);

		}

		#ifdef PLF_CPP20_SUPPORT
		{
			title2("range/fill-insertion when skipblocks are present testing");

			colony<int> i_colony(100, 1);

			i_colony.erase(next(i_colony.begin(), 80), i_colony.end());
			i_colony.insert(10, 2);

			unsigned int total1s = 0, total2s = 0, total3s = 0;

			for (colony<int>::const_iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				total1s += *current == 1;
				total2s += *current == 2;
			}

			failpass("fill-insert smaller than pre-existing skipblock test", static_cast<int>(i_colony.size()) == 90 && total1s == 80 && total2s == 10);

			i_colony.insert(40, 3);
			total2s = 0;

			for (colony<int>::const_iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				total2s += *current == 2;
				total3s += *current == 3;
			}

			failpass("fill-insert larger than pre-existing skipblock test", static_cast<int>(i_colony.size()) == 130 && total2s == 10 && total3s == 40);



			i_colony.erase(next(i_colony.begin(), 120), i_colony.end());
			i_colony.insert({4, 4, 4, 4, 4});

			unsigned int total4s = 0;
			total3s = 0;

			for (colony<int>::const_iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				total3s += *current == 3;
				total4s += *current == 4;
			}

			failpass("range-insert smaller than pre-existing skipblock test", static_cast<int>(i_colony.size()) == 125 && total3s == 30 && total4s == 5);

			i_colony.insert({5, 5, 5, 5, 5, 5, 5, 5, 5, 5});

			total4s = 0;
			unsigned int total5s = 0;

			for (colony<int>::const_iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				total4s += *current == 4;
				total5s += *current == 5;
			}

			failpass("range-insert larger than pre-existing skipblock test", static_cast<int>(i_colony.size()) == 135 && total4s == 5 && total5s == 10);

		}
		#endif

		{
			title2("data() tests");

			colony<int> i_colony(10000, 5);

			int sum1 = 0, sum2 = 0, range1 = 0, range2 = 0;


			// Erase half of all elements and sum the rest:
			for (colony<int>::iterator it = i_colony.begin(); it != i_colony.end(); ++it)
			{
				it = i_colony.erase(it);
				sum1 += *it;
				++range1;
			}


			colony<int>::colony_data *data = i_colony.data();

			// Manually sum using raw memory blocks:
			for (unsigned int block_num = 0; block_num != data->number_of_blocks; ++block_num)
			{
				colony<int>::aligned_element_struct *current_element = data->block_pointers[block_num];
				const unsigned char *bitfield_location = data->bitfield_pointers[block_num];
				const size_t capacity = data->block_capacities[block_num];

				size_t char_index = 0;
				unsigned char offset = 0;

				for (size_t index = 0; index != capacity; ++index, ++current_element)
				{
					if ((bitfield_location[char_index] >> offset) & 1)
					{
						sum2 += *(reinterpret_cast<int *>(current_element));
						++range2;
					}

					if (++offset == 8)
					{
						offset = 0;
						++char_index;
					}
				}
			}

			delete data;

			failpass("Manual summing pass over elements obtained from data()", (sum1 == sum2) && (range1 == range2));
		}


		{
			title2("range/fill insert partial recovery tests");

			colony<exceptions_test> i_colony;
			exceptions_test input_data[10] = {6, 6, 6, 6, 4, 6, 6, 6, 6, 1};

			try
			{
				i_colony.insert(input_data, input_data + 10);
			}
			catch(...)
			{} // do nothing

			int accumulator = 0;
			for (colony<exceptions_test>::iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				accumulator += current->num;
			}

			failpass("fill-insert initializer-list exceptions test", static_cast<int>(i_colony.size()) == 4 && accumulator == 24);

			i_colony.clear();
			i_colony.reserve(20);
			i_colony.insert(20, input_data[9]); // 20 1's
			colony<exceptions_test>::iterator start = next(i_colony.begin(), 5), end = next(i_colony.begin(), 15);
			i_colony.erase(start, end); // Create skip-block in middle of colony sequence

			try
			{
				i_colony.insert(input_data, input_data + 10); // Insert directly into the skipblock and trigger an exception
			}
			catch(...)
			{} // do nothing

			accumulator = 0;
			for (colony<exceptions_test>::iterator current = i_colony.begin(); current != i_colony.end(); ++current)
			{
				accumulator += current->num;
			}

			failpass("fill-insert into skip-block exceptions test", static_cast<int>(i_colony.size()) == 14 && accumulator == 34);
		}
	}


	title1("Test Suite PASS - Press ENTER to Exit");
	getchar();

	return 0;
}
