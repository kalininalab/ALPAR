#include <iostream>
using namespace std;

#ifndef maxK     // max. k-mer length defined
#define maxK 12  // as preprocessor directive
#endif

#define CLASS_NAME   kmerAmino_t
#define STORAGE_TYPE uint5K_t
#define INDEX_TYPE   size5K_t
#define BIT_LENGTH   (5*maxK)
#define LEX_INTEGER_COMPARATORS
#include "byte.h"

#include "util.h"

/**
 * This class contains functions for working with k-mer types.
 */
class kmerAmino {

 private:

    /**
     * This is a bit-mask to erase all bits that exceed the k-mer length.
     */
    static kmerAmino_t mask;

 public:

    /**
     * This is the length of a k-mer (including gap positions).
     */
    static size5K_t k;


    /**
     * This function initializes the k-mer length and bit-mask.
     *
     * @param length k-mer length
     */
    static void init(const size5K_t& length);

    /**
     * This function shifts a k-mer appending a new character to the right.
     *
     * @param kmer bit sequence
     * @param chr right character
     */
    static void shift_left(kmerAmino_t& kmer, char& chr);

    /**
     * This function shifts a k-mer adding a new character to the right.
     *
     * @param kmer bit sequence
     * @param c right character
     */
    static void shift_right(kmerAmino_t& kmer, char& c);
};
