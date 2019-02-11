#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <random>
#include <string_view>
#include <utility>
#include <vector>
#include <cmath>

template <typename T>
std::ostream &operator<<(std::ostream &o, std::vector<T> const &vector)
{ for (auto const &element: vector) o << element << " "; return o; }

using u64 = std::uint64_t; using i32 = std::int32_t;
using u32 = std::uint32_t; using i64 = std::int64_t;
using std::size_t;

namespace benchmark_timer {
    template <typename C = std::chrono::high_resolution_clock> class timer_c {
        const typename C::time_point start_point;
        public: timer_c(): start_point {C::now()} {}
        template <typename R = typename C::duration::rep, typename U = typename C::duration>
        R elapsed() const {
            std::atomic_thread_fence(std::memory_order_relaxed);
            auto counted_time = std::chrono::duration_cast<U>(C::now() - start_point).count();
            std::atomic_thread_fence(std::memory_order_relaxed);
            return static_cast<R>(counted_time);
        }
    };
    using precise_stopwatch     = timer_c<>;
    using system_stopwatch      = timer_c<std::chrono::system_clock>;
    using monotonic_stopwatch   = timer_c<std::chrono::steady_clock>;
} // namespace benchmark_timer

namespace string_generator_utility
{
    std::string generate_string(std::string::size_type length,
                                std::vector<char>&& alphabet = {'a', 'b', 'c', 'd'}) {

                      std::string result;
        static std::random_device rd;
        static    std::mt19937_64 mersenne_twister_generator {rd ()};
        const                auto supremum = alphabet.size() - 1;
        std::uniform_int_distribution
        <std::mt19937_64::result_type> range {0, supremum};

        for (size_t index = 0; index < length; ++index) {
            const auto letter = range(mersenne_twister_generator);
            result += alphabet[letter];
        }
        return result;
    }

    std::vector<std::string> generate_n_strings(u32 vsize, u32 string_length,
                                                std::vector<char>&& alphabet = {'a', 'b', 'c', 'd'}) {
        std::vector<std::string> generated;
        generated.reserve(vsize);
        std::generate_n(std::back_inserter(generated), vsize,
                        [&]() { return generate_string(string_length, std::move(alphabet)); });
        return generated;
    }

} // namespace string_generator_utility

namespace algo
{

  using tuple_statistics_t = std::tuple<u64, std::vector<u64>>;

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
    return {comparisons, matches};
  }

  tuple_statistics_t rabin_karp(std::string_view haystack, std::string_view needle, u64 prime = 101) {

      const auto   needle_size   = needle.size();
      const auto   haystack_size = haystack.size();
      const u64    alphabet_size = 4;
      const u64    h             = std::pow(alphabet_size, needle.size() - 1);

            size_t comparisons   = 0,
                   calculations  = 0;
            u64    needle_hash   = 0,
                   haystack_hash = 0;

      std::vector<u64> matches;

      for (size_t i = 0; i < needle_size; ++i) {
          needle_hash   = (alphabet_size * needle_hash   + needle[i])   % prime;
          haystack_hash = (alphabet_size * haystack_hash + haystack[i]) % prime;
      }

      for (size_t i = 0; i <= haystack_size - needle_size; ++i, ++calculations) {
          if (haystack_hash == needle_hash) {
              size_t p = 0;
              for (size_t j = i; j < i + needle.size(); ++j, ++p, ++comparisons) {
                  if (needle[p] != haystack[j]) break;
              }
              if (p == needle.size()) matches.emplace_back(i);
          }
          haystack_hash = (alphabet_size * (haystack_hash - haystack[i] * h) + haystack[i + needle.size()]) % prime;
          if (haystack_hash < 0) haystack_hash += prime;
      }
      return {comparisons, matches};
  }

  tuple_statistics_t rabin_karp_hash_improved(std::string_view haystack, std::string_view needle) {

    const auto haystack_size = haystack.size();
    const auto needle_size   = needle.size();
           u64 comparisons   = 0;

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
    return {comparisons, matches};
  }

  tuple_statistics_t boyer_moore_horspool(std::string_view haystack, std::string_view needle) {
      
      const      auto      haystack_size = haystack.size();
      const      auto      needle_size   = needle.size();
                  u64      comparisons   = 0;
      std::vector<u64>     matches;
      std::array<u32, 256> bad_match_table;

      for (size_t i = 0; i < 256; ++i)             bad_match_table[i] = needle_size;
      for (size_t i = 0; i < needle_size - 1; ++i) bad_match_table[needle[i]] = needle_size - 1 - i;
      auto skip = 0;

      while (haystack_size - skip >= needle_size) {
          auto i = needle_size - 1;
          while (haystack[skip + i] == needle[i]) {
              if (!i) matches.emplace_back(skip);
              --i;
          }
          skip = skip + bad_match_table[haystack[skip + needle_size - 1]];
      }
      return {comparisons, matches};
  }

} // namespace algo

int main() {
  auto v = string_generator_utility::generate_n_strings(10, 20);

  for (std::size_t index = 0; index < v.size(); ++index) {
    auto [indices_naive, vector_naive] = algo::naive_substring(v[index], "a");
    std::cout << "naive: \t\t" << v[index] << ": " << indices_naive << " \t~ " << vector_naive << "\n";
    auto [indices_rk, vector_rk] = algo::rabin_karp(v[index], "a");
    std::cout << "rk: \t\t" << v[index] << ": " << indices_rk << " \t~ " << vector_rk << "\n";
    auto [indices_rki, vector_rki] = algo::rabin_karp_hash_improved(v[index], "a");
    std::cout << "rki: \t\t" << v[index] << ": " << indices_rki << " \t~ " << vector_rki << "\n\n";
  }

  return 0;
}
