#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <random>
#include <string_view>
#include <utility>
#include <vector>

template <typename T>
std::ostream &operator<<(std::ostream &o, std::vector<T> const &vector)
{ for (auto const &element: vector) o << element << " "; return o; }

namespace string_generator_utility {
    std::string generate_string(std::string::size_type length,
                                std::vector<char>&& alphabet = {'a', 'b', 'c', 'd'}) {
        std::string result;
        static std::random_device rd;
        static std::mt19937_64 mersenne_twister_generator {rd ()};
        const auto supremum = alphabet.size() - 1;
        std::uniform_int_distribution<std::size_t> range {0, supremum};
        for (std::size_t index = 0; index < length; ++index) {
            const auto letter = range(mersenne_twister_generator);
            result += alphabet[letter];
        }
        return result;
    }

    std::vector<std::string> generate_n_strings(std::uint32_t vector_size, std::uint32_t string_length,
                                                std::vector<char>&& alphabet = {'a', 'b', 'c', 'd'}) {
        std::vector<std::string> generated;
        generated.reserve(vector_size);
        std::generate_n(std::back_inserter(generated), vector_size,
                        [&]() { return generate_string(string_length, std::move(alphabet)); });
        return generated;
    }

} // namespace string_generator_utility

namespace algo {

  std::vector<std::int64_t> naive_substring(std::string_view haystack, std::string_view needle) {

    const auto haystack_size = haystack.size();
    const auto needle_size = needle.size();
    assert(haystack_size >= needle_size);

    std::vector<std::int64_t> result;
    for (std::size_t index = 0; index < haystack_size - needle_size + 1; ++index) {
        std::size_t offset = 0;
        for (; offset < needle_size; ++offset) {
            if (haystack[index + offset] != needle[offset])
                break;
        }
        if (offset == needle_size) 
            result.push_back(index);
    }
    return result;
  }

  std::vector<std::int64_t> rabin_karp_hash(std::string_view haystack, std::string_view needle) {

    const auto haystack_size = haystack.size();
    const auto needle_size = needle.size();
    assert(haystack_size >= needle_size);

    std::vector<std::int64_t> matches;
    static const auto hash_function = [](std::string_view::iterator begin, std::string_view::iterator end) {
        std::int64_t hash = 5381;
        for (; begin != end; ++begin)
            hash = ((hash << 5) + hash) + *begin;
        return hash;
    };

    const auto needle_hashed = hash_function(std::begin(needle), std::end(needle));
    for (std::size_t index = 0; index < haystack_size - needle_size + 1; ++index) {
        const auto substring_hash = hash_function
        (
            std::begin(haystack) + index, std::begin(haystack) + index + needle_size
        );
        if (substring_hash == needle_hashed) 
            matches.push_back(index);
    }
    return matches;
  }
} // namespace algo

int main() {
  auto vector = string_generator_utility::generate_n_strings(25, 75);
  std::cout << "naive substring:\n";
  for (std::size_t index = 0; index < vector.size(); ++index)
  {
      std::cout << vector[index] << ": ";
      auto shift = algo::naive_substring(vector[index], "bad");
      std::cout << shift << "\n";
  }
  std::cout << "\nrabin-karp-substring:\n";
  for (std::size_t index = 0; index < vector.size(); ++index)
  {
      std::cout << vector[index] << ": ";
      auto shift = algo::rabin_karp_hash(vector[index], "bad");
      std::cout << shift << "\n";
  }
  return 0;
}
