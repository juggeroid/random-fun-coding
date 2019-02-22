#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <random>
#include <string_view>
#include <utility>
#include <vector>
#include <cmath>
#include <chrono>
#include <atomic>
#include <fstream>
#include <iterator>
#include <numeric>

template <typename T> std::ostream& operator<<(std::ostream& o, std::vector<T> const& vector)
{ for (auto const& element: vector) o << element << " "; return o }

using u64 = std::uint64_t; using i32 = std::int32_t;
using u32 = std::uint32_t; using i64 = std::int64_t;
using std::size_t;
using namespace std::literals;
using tuple_statistics_t = std::tuple<u64, std::vector<u64>>;

namespace utility {
	template
	<typename T, typename = std::enable_if_t<std::is_unsigned_v<T>, T>>
	std::uint32_t vector_median(std::vector<T> const& vector) {
		const auto size = vector.size();
		std::accumulate(std::begin(vector), std::end(vector), 0) / size;
	}
	template <typename T> void vector_output_to_file(std::vector<T> const& vector, std::string const& filename) {
		std::ofstream output_file {filename, std::ios::app};
		std::ostream_iterator<T> output_iterator(output_file, "\n");
		std::copy(std::begin(vector), std::end(vector), output_iterator);
		return;
	}
}

namespace benchmark_timer {
	template <typename C = std::chrono::high_resolution_clock> class timer_c {
	const typename C::time_point start_point;
	public: timer_c() : start_point(C::now()) {}
		template <typename R = typename C::duration::rep, typename U = typename C::duration>
		R elapsed() const {
			std::atomic_thread_fence(std::memory_order_relaxed);
			auto counted_time = std::chrono::duration_cast<U>(C::now() - start_point).count();
			std::atomic_thread_fence(std::memory_order_relaxed);
			return static_cast<R>(counted_time);
		}
	};
	using precise_stopwatch   = timer_c<>;
	using system_stopwatch    = timer_c<std::chrono::system_clock>;
	using monotonic_stopwatch = timer_c<std::chrono::steady_clock>;
}


namespace string_generator_utility {

	std::string generate_string(std::string::size_type length,
		std::vector<char> const& alphabet = { 'a', 'b', 'c', 'd' }) {

		const auto alphabet_size = alphabet.size();
		if (!alphabet_size) throw std::invalid_argument("alphabet cannot be empty");

			      std::string result;
		static std::random_device rd;
		static    std::mt19937_64 alphabet_index_generator{ rd() };
		const                auto supremum = alphabet.size() - 1;
		std::uniform_int_distribution
		<std::mt19937_64::result_type> range {0, supremum};

		result.reserve(length);

		for (size_t index = 0; index < length; ++index) {
			const auto letter = range(alphabet_index_generator);
			result += alphabet[letter];
		}
		return result;
	}

	std::vector<std::string> generate_n_strings(u32 vsize, u32 string_length,
		std::vector<char> const& alphabet = { 'a', 'b', 'c', 'd' }) {
		std::vector<std::string> generated;
		generated.reserve(vsize);
		std::generate_n(std::back_inserter(generated), vsize, [&]() { return generate_string(string_length, alphabet); });
		return generated;
	}

} // namespace string_generator_utility

namespace algo {

	tuple_statistics_t naive_substring(std::string_view haystack, std::string_view needle) {

		const auto haystack_size = haystack.size();
		const auto needle_size   = needle.size();
		       u64 comparisons   = 0;

		std::vector<u64> matches;
		for (size_t index = 0; index < haystack_size - needle_size + 1; ++index) {
			size_t offset = 0;
			for (; offset < needle_size; ++offset, ++comparisons) {
				if (haystack[index + offset] != needle[offset]) break;
			}
			if (offset == needle_size) matches.emplace_back(index);
		}
		return { comparisons, matches };
	}

	tuple_statistics_t rabin_karp(std::string_view haystack, std::string_view needle, i32 alphabet_size = 4) {

		const           auto  needle_size   = needle.size();
		const           auto  haystack_size = haystack.size();
		static constexpr i32  prime         = 101;
		const            i32  h             = std::pow(alphabet_size, needle.size() - 1);
		std::vector<u64> matches;

		size_t  comparisons   = 0;
		i64     needle_hash   = 0,
			haystack_hash = 0;

		for (size_t i = 0; i < needle_size; ++i) {
			needle_hash = (alphabet_size * needle_hash + needle[i]) % prime;
			haystack_hash = (alphabet_size * haystack_hash + haystack[i]) % prime;
		}

		for (size_t i = 0; i <= haystack_size - needle_size; ++i, ++comparisons) {
			if (haystack_hash == needle_hash) {
				size_t p = 0;
				for (size_t j = i; j < i + needle.size(); ++j, ++p, ++comparisons) {
					if (needle[p] != haystack[j]) break;
				}
				if (p == needle_size) matches.emplace_back(i);
			}
			haystack_hash = (alphabet_size * (haystack_hash - haystack[i] * h) + haystack[i + needle.size()]) % prime;
			if (haystack_hash < 0) haystack_hash += prime;
		}
		return { comparisons, matches };
	}

	tuple_statistics_t rabin_karp_hash_improved(std::string_view haystack, std::string_view needle) {

		const auto haystack_size = haystack.size();
		const auto needle_size   = needle.size();
		u64 comparisons          = 0;

		std::vector<u64> matches;
		static const auto hash_function = [](std::string_view::iterator begin, std::string_view::iterator end) {
			i64 hash = 5381;
			for (; begin != end; ++begin) hash = ((hash << 5) + hash) + *begin;
			return hash;
		};

		const auto needle_hashed = hash_function(std::begin(needle), std::end(needle));
		for (std::size_t index = 0; index < haystack_size - needle_size + 1; ++index, ++comparisons) {
			const auto substring_hash = hash_function
			(std::begin(haystack) + index, std::begin(haystack) + index + needle_size);
			if (substring_hash == needle_hashed) matches.emplace_back(index);
		}
		return { comparisons, matches };
	}

	tuple_statistics_t boyer_moore_horspool(std::string_view haystack, std::string_view needle) {

		const      auto      haystack_size = haystack.size();
		const      auto      needle_size   = needle.size();
		            u64      comparisons   = 0;
		std::vector<u64>     matches;
		u32 bad_match_table[256];

		for (size_t i = 0; i < 256; ++i)             bad_match_table[i] = needle_size;
		for (size_t i = 0; i < needle_size - 1; ++i) bad_match_table[needle[i]] = needle_size - 1 - i;
		auto skip = 0;

		while (haystack_size - skip >= needle_size) {
			auto i = needle_size - 1;
			while (haystack[skip + i] == needle[i]) {
				++comparisons;
				if (!i) matches.emplace_back(skip);
				--i;
			}
			skip = skip + bad_match_table[haystack[skip + needle_size - 1]];
		}
		return { comparisons, matches };
	}

	tuple_statistics_t knuth_morris_pratt(std::string_view haystack, std::string_view needle) {

		const      auto    haystack_size = haystack.size();
		const      auto    needle_size = needle.size();
		            u64    comparisons = 0;
		std::vector<u64>   matches;
		    std::size_t    index = 0,
				   offset = 0;
		std::vector<char> preproccessing_array;
		preproccessing_array.reserve(needle_size);

		const auto longest_proper_suffix = [&preproccessing_array, &needle_size, &needle]() {
			auto length = 0, i = 1;
			preproccessing_array[0] = 0;
			while (i < needle_size) {
				if (needle[i] == needle[length]) {
					++length; preproccessing_array[i] = length; ++i;
				}
				else {
					if (length != 0) length = preproccessing_array[length - 1];
					else             preproccessing_array[i] = 0; ++i;
				}
			}
		};

		longest_proper_suffix();
		while (index < haystack_size) {
			if (needle[offset] == haystack[index]) { ++offset; ++index; ++comparisons; }
			if (offset == needle_size) { matches.emplace_back(index - offset); offset = preproccessing_array[offset - 1]; }
			else if (index < haystack_size&& needle[offset] != haystack[index]) {
				if (offset != 0) offset = preproccessing_array[offset - 1];
				else ++index;
			}
		}
		return { comparisons, matches };
	}

} // namespace algo

void test_framework(std::size_t vector_size, std::size_t string_size, std::size_t substring_size) {

	auto v = string_generator_utility::generate_n_strings(vector_size, string_size);
	auto substr = string_generator_utility::generate_string(substring_size);
	std::vector<std::size_t> comparisons_total;
	comparisons_total.reserve(vector_size);

	std::cout << "[substr] success: finished generating the data set!\n\n";
	//! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated
	benchmark_timer::precise_stopwatch timer;
	comparisons_total.clear();
	for (std::size_t index = 0; index < v.size(); ++index)
	{
		auto [naive_comparisons, naive_matches] = algo::naive_substring(v[index], substr);
		comparisons_total.push_back(naive_comparisons);
	}
	auto elapsed = timer.elapsed() / 1000 / 1000;
	std::cout << "[substr] naive_substring elapsed time:\t\t" << elapsed << " ms\n";
	std::cout << "[substr] naive_substring comp_avg:\t\t" << utility::vector_median(comparisons_total) << "\n";
	utility::vector_output_to_file(std::vector<long long>{elapsed}, "naive_time.txt");

	//! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated
	benchmark_timer::precise_stopwatch timer_2;
	comparisons_total.clear();
	for (std::size_t index = 0; index < v.size(); ++index)
	{
		auto [kmp_comparisons, kmp_matches] = algo::knuth_morris_pratt(v[index], substr);
		comparisons_total.push_back(kmp_comparisons);
	}
	auto elapsed_2 = timer_2.elapsed() / 1000 / 1000;
	std::cout << "[substr] knuth_morris_pratt elapsed time:\t" << elapsed_2 << " ms\n";
	std::cout << "[substr] knuth_morris_pratt comp_avg:\t\t" << utility::vector_median(comparisons_total) << "\n";
	utility::vector_output_to_file(std::vector<long long>{elapsed_2}, "kmp_time.txt");


	//! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated
	benchmark_timer::precise_stopwatch timer_3;
	comparisons_total.clear();
	for (std::size_t index = 0; index < v.size(); ++index)
	{
		auto [rk_comparisons, rk_matches] = algo::rabin_karp(v[index], substr);
		comparisons_total.push_back(rk_comparisons);
	}
	auto elapsed_3 = timer_3.elapsed() / 1000 / 1000;
	std::cout << "[substr] rabin_karp elapsed time:\t\t" << elapsed_3 << " ms\n";
	std::cout << "[substr] rabin_karp comp_avg:\t\t\t" << utility::vector_median(comparisons_total) << "\n";
	utility::vector_output_to_file(std::vector<long long>{elapsed_3}, "rk_time.txt");

	//! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated
	/*
	benchmark_timer::precise_stopwatch timer_4;
	comparisons_total.clear();
	for (std::size_t index = 0; index < v.size(); ++index)
	{
		auto [rki_comparisons,     rki_matches] = algo::rabin_karp_hash_improved (v[index], substr);
		comparisons_total.push_back(rki_comparisons);
	}
	std::cout << "[substr] rabin_karp_improved elapsed time:\t" << timer_4.elapsed() / 1000 / 1000 << " ms\n";
	std::cout << "[substr] rabin_karp_improved comp_avg:\t\t"   << utility::vector_median(comparisons_total) << "\n\n";
	*/

	//! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated //! @deprecated
	benchmark_timer::precise_stopwatch timer_5;
	comparisons_total.clear();
	for (std::size_t index = 0; index < v.size(); ++index)
	{
		auto [bmh_comparisons, bmh_matches] = algo::boyer_moore_horspool(v[index], substr);
		comparisons_total.push_back(bmh_comparisons);
	}
	auto elapsed_5 = timer_5.elapsed() / 1000 / 1000;
	std::cout << "[substr] boyer_moore_horspool elapsed time:\t" << elapsed_5 << " ms\n";
	std::cout << "[substr] boyer_moore_horspool comp_avg:\t\t" << utility::vector_median(comparisons_total) << "\n";
	utility::vector_output_to_file(std::vector<long long>{elapsed_5}, "bmh_time.txt");

	benchmark_timer::precise_stopwatch timer_6;
	comparisons_total.clear();
	for (std::size_t index = 0; index < v.size(); ++index)
	{
		auto d = std::search(std::begin(v[index]), std::end(v[index]), std::boyer_moore_horspool_searcher(substr.begin(), substr.end()));
	}
	auto elapsed_6 = timer_6.elapsed() / 1000 / 1000;
	std::cout << "[substr] std::search (BMH) elapsed time:\t" << elapsed_6 << " ms\n";
}

int main() {
	std::cout << "[substr] ok ready\n";
	test_framework(2000, 100000, 100);
	return 0;
}
