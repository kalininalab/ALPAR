/**
 * MIT License
 *
 * Copyright (c) 2017 Thibaut Goetghebuer-Planchon <tessil@gmx.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef TSL_SPARSE_GROWTH_POLICY_H
#define TSL_SPARSE_GROWTH_POLICY_H

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <ratio>
#include <stdexcept>

namespace tsl {
 namespace sh {

    /**
     * Grow the hash table by a factor of GrowthFactor keeping the bucket count to a
     * power of two. It allows the table to use a mask operation instead of a modulo
     * operation to map a hash to a bucket.
     *
     * GrowthFactor must be a power of two >= 2.
     */
    template <std::size_t GrowthFactor>
    class power_of_two_growth_policy {
     public:
        /**
         * Called on the hash table creation and on rehash. The number of buckets for
         * the table is passed in parameter. This number is a minimum, the policy may
         * update this value with a higher value if needed (but not lower).
         *
         * If 0 is given, min_bucket_count_in_out must still be 0 after the policy
         * creation and bucket_for_hash must always return 0 in this case.
         */
        explicit power_of_two_growth_policy(std::size_t &min_bucket_count_in_out) {
            if (min_bucket_count_in_out > max_bucket_count()) {
                throw std::length_error("The hash table exceeds its maximum size.");
            }

            if (min_bucket_count_in_out > 0) {
                min_bucket_count_in_out =
                        round_up_to_power_of_two(min_bucket_count_in_out);
                m_mask = min_bucket_count_in_out - 1;
            } else {
                m_mask = 0;
            }
        }

        /**
         * Return the bucket [0, bucket_count()) to which the hash belongs.
         * If bucket_count() is 0, it must always return 0.
         */
        std::size_t bucket_for_hash(std::size_t hash) const noexcept {
            return hash & m_mask;
        }

        /**
         * Return the number of buckets that should be used on next growth.
         */
        std::size_t next_bucket_count() const {
            if ((m_mask + 1) > max_bucket_count() / GrowthFactor) {
                throw std::length_error("The hash table exceeds its maximum size.");
            }

            return (m_mask + 1) * GrowthFactor;
        }

        /**
         * Return the maximum number of buckets supported by the policy.
         */
        std::size_t max_bucket_count() const {
            // Largest power of two.
            return (std::numeric_limits<std::size_t>::max() / 2) + 1;
        }

        /**
         * Reset the growth policy as if it was created with a bucket count of 0.
         * After a clear, the policy must always return 0 when bucket_for_hash is
         * called.
         */
        void clear() noexcept { m_mask = 0; }

     private:
        static std::size_t round_up_to_power_of_two(std::size_t value) {
            if (is_power_of_two(value)) {
                return value;
            }

            if (value == 0) {
                return 1;
            }

            --value;
            for (std::size_t i = 1; i < sizeof(std::size_t) * CHAR_BIT; i *= 2) {
                value |= value >> i;
            }

            return value + 1;
        }

        static constexpr bool is_power_of_two(std::size_t value) {
            return value != 0 && (value & (value - 1)) == 0;
        }

     protected:
        static_assert(is_power_of_two(GrowthFactor) && GrowthFactor >= 2,
                      "GrowthFactor must be a power of two >= 2.");

        std::size_t m_mask;
    };

    /**
     * Grow the hash table by GrowthFactor::num / GrowthFactor::den and use a modulo
     * to map a hash to a bucket. Slower but it can be useful if you want a slower
     * growth.
     */
    template <class GrowthFactor = std::ratio<3, 2>>
    class mod_growth_policy {
     public:
        explicit mod_growth_policy(std::size_t &min_bucket_count_in_out) {
            if (min_bucket_count_in_out > max_bucket_count()) {
                throw std::length_error("The hash table exceeds its maximum size.");
            }

            if (min_bucket_count_in_out > 0) {
                m_mod = min_bucket_count_in_out;
            } else {
                m_mod = 1;
            }
        }

        std::size_t bucket_for_hash(std::size_t hash) const noexcept {
            return hash % m_mod;
        }

        std::size_t next_bucket_count() const {
            if (m_mod == max_bucket_count()) {
                throw std::length_error("The hash table exceeds its maximum size.");
            }

            const double next_bucket_count =
                    std::ceil(double(m_mod) * REHASH_SIZE_MULTIPLICATION_FACTOR);
            if (!std::isnormal(next_bucket_count)) {
                throw std::length_error("The hash table exceeds its maximum size.");
            }

            if (next_bucket_count > double(max_bucket_count())) {
                return max_bucket_count();
            } else {
                return std::size_t(next_bucket_count);
            }
        }

        std::size_t max_bucket_count() const { return MAX_BUCKET_COUNT; }

        void clear() noexcept { m_mod = 1; }

     private:
        static constexpr double REHASH_SIZE_MULTIPLICATION_FACTOR =
                1.0 * GrowthFactor::num / GrowthFactor::den;
        static const std::size_t MAX_BUCKET_COUNT =
                std::size_t(double(std::numeric_limits<std::size_t>::max() /
                                   REHASH_SIZE_MULTIPLICATION_FACTOR));

        static_assert(REHASH_SIZE_MULTIPLICATION_FACTOR >= 1.1,
                      "Growth factor should be >= 1.1.");

        std::size_t m_mod;
    };

    /**
     * Grow the hash table by using prime numbers as bucket count. Slower than
     * tsl::sh::power_of_two_growth_policy in general but will probably distribute
     * the values around better in the buckets with a poor hash function.
     *
     * To allow the compiler to optimize the modulo operation, a lookup table is
     * used with constant primes numbers.
     *
     * With a switch the code would look like:
     * \code
     * switch(iprime) { // iprime is the current prime of the hash table
     *     case 0: hash % 5ul;
     *             break;
     *     case 1: hash % 17ul;
     *             break;
     *     case 2: hash % 29ul;
     *             break;
     *     ...
     * }
     * \endcode
     *
     * Due to the constant variable in the modulo the compiler is able to optimize
     * the operation by a series of multiplications, substractions and shifts.
     *
     * The 'hash % 5' could become something like 'hash - (hash * 0xCCCCCCCD) >> 34)
     * * 5' in a 64 bits environment.
     */
    class prime_growth_policy {
     public:
        explicit prime_growth_policy(std::size_t &min_bucket_count_in_out) {
            auto it_prime = std::lower_bound(primes().begin(), primes().end(),
                                             min_bucket_count_in_out);
            if (it_prime == primes().end()) {
                throw std::length_error("The hash table exceeds its maximum size.");
            }

            m_iprime =
                    static_cast<unsigned int>(std::distance(primes().begin(), it_prime));
            if (min_bucket_count_in_out > 0) {
                min_bucket_count_in_out = *it_prime;
            } else {
                min_bucket_count_in_out = 0;
            }
        }

        std::size_t bucket_for_hash(std::size_t hash) const noexcept {
            return mod_prime()[m_iprime](hash);
        }

        std::size_t next_bucket_count() const {
            if (m_iprime + 1 >= primes().size()) {
                throw std::length_error("The hash table exceeds its maximum size.");
            }

            return primes()[m_iprime + 1];
        }

        std::size_t max_bucket_count() const { return primes().back(); }

        void clear() noexcept { m_iprime = 0; }

     private:
      #if SIZE_MAX <= UINT32_MAX
         static const std::array<std::size_t, 40> &primes() {
            static const std::array<std::size_t, 40> PRIMES = {{
               1ul, 5ul, 17ul, 29ul, 37ul, 53ul, 67ul, 79ul, 97ul, 131ul, 193ul, 257ul, 389ul, 521ul,
               769ul, 1031ul, 1543ul, 2053ul, 3079ul, 6151ul, 12289ul, 24593ul, 49157ul, 98317ul, 196613ul,
               393241ul, 786433ul, 1572869ul, 3145739ul, 6291469ul, 12582917ul, 25165843ul, 50331653ul,
               100663319ul, 201326611ul, 402653189ul, 805306457ul, 1610612741ul, 3221225473ul, 4294967291ul
            }};
            return PRIMES;
         }
      #else /* SIZE_MAX <= UINT64_MAX */
         static const std::array<std::size_t, 71> &primes() {
            static const std::array<std::size_t, 71> PRIMES = {{
               1ull, 5ull, 17ull, 29ull, 41ull, 53ull, 73ull, 97ull, 131ull, 193ull, 257ull, 389ull, 521ull, 769ull, 1031ull, 1543ull, 2053ull, 3079ull,
               6151ull, 12289ull, 24593ull, 49157ull, 98317ull, 196613ull, 393241ull, 786433ull, 1572869ull, 3145739ull, 6291469ull, 12582917ull,
               25165843ull, 50331653ull, 100663319ull, 201326611ull, 402653189ull, 805306457ull, 1610612741ull, 3221225473ull, 6442450939ull,
               12884901877ull, 25769803751ull, 51539607503ull, 103079215007ull, 206158430017ull, 412316860027ull, 824633720047ull, 1649267440087ull,
               3298534880167ull, 6597069760331ull, 13194139520663ull, 26388279041353ull, 52776558082697ull, 105553116165401ull, 211106232330821ull,
               422212464661667ull, 844424929323337ull, 1688849858646677ull, 3377699717293357ull, 6755399434586687ull, 13510798869173377ull,
               27021597738346709ull, 54043195476693431ull, 108086390953386863ull, 216172781906773733ull, 432345563813547511ull, 864691127627095019ull,
               1729382255254190023ull, 3458764510508380061ull, 6917529021016760191ull, 13835058042033520309ull, 18446744073709551557ull
            }};
            return PRIMES;
         }
      #endif

     // MOD_PRIME[iprime](hash) returns hash % PRIMES[iprime]. This table allows
     // for faster modulo as the compiler can optimize the modulo code better
     // with a constant known at the compilation.
      #if SIZE_MAX <= UINT32_MAX
         static const std::array<std::size_t (*)(std::size_t), 40> &mod_prime() {
            static const std::array<std::size_t (*)(std::size_t), 40> MOD_PRIME = {{
               &mod<0>,  &mod<1>,  &mod<2>,  &mod<3>,  &mod<4>,  &mod<5>,  &mod<6>,  &mod<7>,  &mod<8>,  &mod<9>,
               &mod<10>, &mod<11>, &mod<12>, &mod<13>, &mod<14>, &mod<15>, &mod<16>, &mod<17>, &mod<18>, &mod<19>,
               &mod<20>, &mod<21>, &mod<22>, &mod<23>, &mod<24>, &mod<25>, &mod<26>, &mod<27>, &mod<28>, &mod<29>,
               &mod<30>, &mod<31>, &mod<32>, &mod<33>, &mod<34>, &mod<35>, &mod<36>, &mod<37>, &mod<38>, &mod<39>
            }};
            return MOD_PRIME;
         }
      #else /* SIZE_MAX <= UINT64_MAX */
         static const std::array<std::size_t (*)(std::size_t), 71> &mod_prime() {
            static const std::array<std::size_t (*)(std::size_t), 71> MOD_PRIME = {{
               &mod<0>, &mod<1>, &mod<2>, &mod<3>, &mod<4>, &mod<5>, &mod<6>, &mod<7>, &mod<8>, &mod<9>, &mod<10>,
               &mod<11>, &mod<12>, &mod<13>, &mod<14>, &mod<15>, &mod<16>, &mod<17>, &mod<18>, &mod<19>, &mod<20>,
               &mod<21>, &mod<22>, &mod<23>, &mod<24>, &mod<25>, &mod<26>, &mod<27>, &mod<28>, &mod<29>, &mod<30>,
               &mod<31>, &mod<32>, &mod<33>, &mod<34>, &mod<35>, &mod<36>, &mod<37>, &mod<38>, &mod<39>, &mod<40>,
               &mod<41>, &mod<42>, &mod<43>, &mod<44>, &mod<45>, &mod<46>, &mod<47>, &mod<48>, &mod<49>, &mod<50>,
               &mod<51>, &mod<52>, &mod<53>, &mod<54>, &mod<55>, &mod<56>, &mod<57>, &mod<58>, &mod<59>, &mod<60>,
               &mod<61>, &mod<62>, &mod<63>, &mod<64>, &mod<65>, &mod<66>, &mod<67>, &mod<68>, &mod<69>, &mod<70>
            }};
            return MOD_PRIME;
         }
      #endif

        template <unsigned int IPrime>
        static std::size_t mod(std::size_t hash) {
            return hash % primes()[IPrime];
        }

     private:
        unsigned int m_iprime;
    };

 }  // namespace sh
}  // namespace tsl

#endif
