#include <iostream>
using namespace std;

#ifndef maxN     // max. color number defined
#define maxN 64  // as preprocessor directive
#endif

#define CLASS_NAME   color_t
#define STORAGE_TYPE uint1N_t
#define INDEX_TYPE   size1N_t
#define BIT_LENGTH   (1*maxN)
#define SET_ELEMENT_COMPARATORS
#include "byte.h"

/**
 * This class contains functions for working with color types.
 */
class color {

 private:

    /**
     * This is a bit-mask to erase all bits that exceed the color number.
     */
    static color_t mask;

 public:

    /**
     * This is the number of colors.
     */
    static size1N_t n;

    /**
     * This function initializes the color number and bit-mask.
     *
     * @param number color number
     */
    static void init(const size1N_t& number);

    /**
     * This function shifts a color appending a new bit char to the right.
     *
     * @param color bit sequence
     * @param chr right color bit
     */
    static void shift(color_t& color, const char& chr);

    /**
     * This function unshifts a color returning the bit char on the right.
     *
     * @param color bit sequence
     * @param chr right color bit
     */
    static void unshift(color_t& color, char& chr);

    /**
     * This function constructs the bit complement of a given color set.
     *
     * @param color bit sequence
     */
    static void complement(color_t& color);

    /**
     * This function constructs the representative of a given color set.
     *
     * @param color bit sequence
     * @return 1 if inverted, 0 otherwise
     */
    static bool represent(color_t& color);

    /**
     * This function tests if two splits of colors are compatible.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @return true, if compatible
     */
    static bool is_compatible(const color_t& c1, const color_t& c2);

    /**
     * This function tests if three splits of colors are weakly compatible.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @param c3 bit sequence
     * @return true, if weakly compatible
     */
    static bool is_weakly_compatible(const color_t& c1, const color_t& c2, const color_t& c3);

 protected:

};
