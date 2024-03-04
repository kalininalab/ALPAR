#include "kmerAmino.h"

/**
 * This is the length of a k-mer.
 */
size5K_t kmerAmino::k;              // length of a k-mer
kmerAmino_t kmerAmino::mask;      // bit-mask to earase all bits that exceed the k-mer length

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param kmer_length k-mer length
 */
void kmerAmino::init(const size5K_t& kmer_length) {
    k = kmer_length; mask = 0b0u;
    for (uint64_t i = 0; i < 5*k; ++i) {
        mask <<= 01u;    // fill all bits within the k-mer length with ones
        mask |= 01u;    // the remaining zero bits can be used to mask bits
    }
}

/**
 * This function shifts a k-mer adding a new character to the left.
 *
 * @param kmer bit sequence
 * @param chr left character
 * @return right character
 */
void kmerAmino::shift_left(kmerAmino_t& kmer, char& chr) {
    uint64_t left = util::amino_char_to_bits(chr);    // new leftmost character
    kmer >>= 05u;    // shift all current bits to the right by five positions
    kmer |= left << (5*k-05u);    // encode the new character within the leftmost five bits
    kmer &= mask;    // set all bits to zero that exceed the k-mer length
}

/**
 * This function shifts a k-mer adding a new character to the right.
 *
 * @param kmer bit sequence
 * @param chr right character
 * @return left character
 */
void kmerAmino::shift_right(kmerAmino_t& kmer, char& chr) {
    uint64_t right = util::amino_char_to_bits(chr);    // new rightmost character
    kmer <<= 05u;    // shift all current bits to the left by five positions
    kmer |= right;    // encode the new character within the rightmost five bits
    kmer &= mask;    // set all bits to zero that exceed the k-mer length
}
