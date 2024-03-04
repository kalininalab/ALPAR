#include "kmer.h"

/*
 * This class contains functions for working with k-mer types.
 */
size2K_t kmer::k;      // length of a k-mer (including gap positions)
kmer_t   kmer::mask;   // bit-mask to erase all bits that exceed the k-mer length

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param length k-mer length
 */
void kmer::init(const size2K_t& length) {
    k = length; mask = 0b0u;
    for (size2K_t i = 0; i < k; ++i)  // fill all bits within the k-mer length with ones
        (mask <<= 02u) |= 0b11u;     // the remaining zero bits can be used to mask bits
}

/**
 * This function shifts a k-mer appending a new character to the right.
 *
 * @param kmer bit sequence
 * @param uint_fast8_t right character in binary-code
 */
void kmer::shift(kmer_t& kmer, uint_fast8_t& right) {
    kmer <<= 02u;    // shift all current bits to the left by two positions
    kmer |= right;    // encode the new rightmost character
    kmer &= mask;    // set all bits to zero that exceed the k-mer length
}

/**
 * This function shifts a k-mer appending a new character to the right.
 *
 * @param kmer bit sequence
 * @param chr right character
 */
void kmer::shift(kmer_t& kmer, char& c_right) {
    kmer <<= 02u;    // shift all current bits to the left by two positions
    kmer |= util::char_to_bits(c_right);    // encode the new rightmost character
    kmer &= mask;    // set all bits to zero that exceed the k-mer length
}



/**
 * This function unshifts a k-mer returning the character on the right.
 *
 * @param kmer bit sequence
 * @param chr right character
 */
void kmer::unshift(kmer_t& kmer, char& chr) {
    kmer >>= 02u;    // shift all current bits to the right by two positions
}

/**
 * This function constructs the reverse complement of a given k-mer.
 *
 * @param kmer bit sequence
 */
void kmer::reverse_complement(kmer_t& kmer) {
    kmer_t bits = ~kmer;    // flip the original k-mer
    kmer_t rcmp;
    for (size2K_t i = 0; i < k; ++i) {
        rcmp <<= 02u;    // shift in the first base
        rcmp |= bits & 0b11u;    // transfer the character
        bits >>= 02u;    // shift out the last base
    }
    kmer = rcmp;
}

/**
 * This function constructs the r.c. representative of a given k-mer.
 *
 * @param kmer bit sequence
 * @return 1 if inverted, 0 otherwise
 */
bool kmer::reverse_represent(kmer_t& kmer) {
    kmer_t bits = ~kmer;    // flip the original k-mer
    kmer_t rcmp;
    for (size2K_t i = 0; i < k; ++i) {
        rcmp <<= 02u;    // shift in the first base
        rcmp |= bits & 0b11u;    // transfer the character
        bits >>= 02u;    // shift out the last base
    }
    // return the lexicographically smaller
    if  (kmer < rcmp) return false;    // not reversed
    else kmer = rcmp; return true;    // reversed
}

