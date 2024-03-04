#if !defined(CLASS_NAME) // must be defined in the including header
    #error "CLASS_NAME is not defined (byte.h)"
#endif
#if !defined(STORAGE_TYPE) // must be defined in the including header
    #error "STORAGE_TYPE is not defined (byte.h)"
#endif
#if !defined(INDEX_TYPE) // must be defined in the including header
    #error "INDEX_TYPE is not defined (byte.h)"
#endif
#if !defined(BIT_LENGTH) // must be defined in the including header
    #error "BIT_LENGTH is not defined (byte.h)"
#endif

#if __GNUC__ <= 8 // workaround for call to non-constexpr function error
    #define _constexpr 
#else
    #define _constexpr constexpr
#endif

#if BIT_LENGTH <= 8
    #define STORAGE_BITS 8
    typedef uint_least8_t STORAGE_TYPE;
#elif BIT_LENGTH <= 16
    #define STORAGE_BITS 16
    typedef uint_least16_t STORAGE_TYPE;
#elif BIT_LENGTH <= 32
    #define STORAGE_BITS 32
    typedef uint_least32_t STORAGE_TYPE;
#else // _LENGTH <= 64
    #define STORAGE_BITS 64
    typedef uint_least64_t STORAGE_TYPE;
#endif

#if BIT_LENGTH <= 255
    typedef uint_fast8_t INDEX_TYPE;
#elif BIT_LENGTH <= 65535
    typedef uint_fast16_t INDEX_TYPE;
#elif BIT_LENGTH <= 4294967295
    typedef uint_fast32_t INDEX_TYPE;
#else // _LENGTH <= 18446744073709551615
    typedef uint_fast64_t INDEX_TYPE;
#endif

#define MAX_STORAGE_BITS 64 // threshold to switch from single to array representation
#define ARRAY_LENGTH ((BIT_LENGTH / STORAGE_BITS) + (bool)(BIT_LENGTH % STORAGE_BITS))

#if defined(__has_include)
    #if __has_include(<immintrin.h>)
        #include <immintrin.h>
        #if defined(__POPCNT__)
            #if STORAGE_BITS <= 32
                #define _popcnt(_X) _mm_popcnt_u32(_X)
            #elif STORAGE_BITS <= 64
                #define _popcnt(_X) _mm_popcnt_u64(_X)
            #endif
        #endif
        #if defined(__BMI__)
            #if STORAGE_BITS <= 32
                #define _tzcnt(_X) (_X? _tzcnt_u32(_X): STORAGE_BITS)
            #elif STORAGE_BITS <= 64
                #define _tzcnt(_X) (_X? _tzcnt_u64(_X): STORAGE_BITS)
            #endif
        #endif
        #if defined(__BMI2__)
            #if STORAGE_BITS <= 32
                #define _pext(_X,_Y) _pext_u32(_X,_Y)
                #define _pdep(_X,_Y) _pdep_u32(_X,_Y)
            #elif STORAGE_BITS <= 64
                #define _pext(_X,_Y) _pext_u64(_X,_Y)
                #define _pdep(_X,_Y) _pdep_u64(_X,_Y)
            #endif
        #endif
    #endif
#endif

#if !defined(_popcnt)
    #define _popcnt(_X)\
    [] (const STORAGE_TYPE& X) -> INDEX_TYPE {\
    if (!X) return 0; if (!~X) return STORAGE_BITS;\
        INDEX_TYPE _count; STORAGE_TYPE _temp = X;\
        for (_count = 0; _temp; ++_count) {\
           _temp &= _temp-1;\
        } return _count; } (_X)
#endif
#if !defined(_tzcnt)
    #define _tzcnt(_X)\
    [] (const STORAGE_TYPE& X) -> INDEX_TYPE {\
    if (!X) return STORAGE_BITS; if (!~X) return 0;\
        INDEX_TYPE _count; STORAGE_TYPE _temp = X^(X-1);\
        for (_count = -1; _temp; ++_count) {\
           _temp >>= 1;\
        } return _count; } (_X)
#endif
#if !defined(_pext)
    #define _pext(_X,_Y)\
    [] (const STORAGE_TYPE& X, const STORAGE_TYPE& Y) -> STORAGE_TYPE {\
    if (!X || !Y) return 0; if (!~Y) return X;\
        STORAGE_TYPE _byte = 0, _mask = Y;\
        for (STORAGE_TYPE _pos = 1; _mask; _pos <<= 1) {\
            if (X & _mask & -_mask) _byte |= _pos;\
           _mask &= _mask-1;\
        } return _byte; } (_X,_Y)
#endif
#if !defined(_pdep)
    #define _pdep(_X,_Y)\
    [] (const STORAGE_TYPE& X, const STORAGE_TYPE& Y) -> STORAGE_TYPE {\
    if (!X || !Y) return 0; if (!~Y) return X;\
        STORAGE_TYPE _byte = 0, _mask = Y;\
        for (STORAGE_TYPE _pos = 1; _mask; _pos <<= 1) {\
            if (X & _pos) _byte |= _mask & -_mask;\
           _mask &= _mask-1;\
        } return _byte; } (_X,_Y)
#endif

// ############################ BEGIN CLASS DEFINITION ############################ //

class CLASS_NAME {
 private:
   #if BIT_LENGTH <= MAX_STORAGE_BITS
     STORAGE_TYPE byte{};
   #else
     STORAGE_TYPE byte[ARRAY_LENGTH]{};
   #endif
    friend struct std::hash<CLASS_NAME>;

 public:
    constexpr CLASS_NAME() noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
       #else // empty default constructor
       #endif
    }
    constexpr CLASS_NAME(const STORAGE_TYPE& value) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte = value;
       #else
         byte[0] = value;
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             byte[i] = 0;
       #endif
    }
    constexpr CLASS_NAME(const CLASS_NAME& other) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte = other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             byte[i] = other.byte[i];
       #endif
    }
    constexpr CLASS_NAME& operator=(const STORAGE_TYPE& value) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte = value;
       #else
         byte[0] = value;
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             byte[i] = 0;
       #endif
        return *this;
    }
    constexpr CLASS_NAME& operator=(const CLASS_NAME& other) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte = other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             byte[i] = other.byte[i];
       #endif
        return *this;
    }

 // ########################## BIT& ASSIGNMENT OPERATORS ########################## //

    constexpr CLASS_NAME& operator&=(const STORAGE_TYPE& value) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte &= value;
       #else
         byte[0] &= value;
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             byte[i] = 0;
       #endif
        return *this;
    }
    constexpr CLASS_NAME& operator&=(const CLASS_NAME& other) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte &= other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             byte[i] &= other.byte[i];
       #endif
        return *this;
    }
    constexpr CLASS_NAME& operator|=(const STORAGE_TYPE& value) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte |= value;
       #else
         byte[0] |= value;
      // for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
      //     byte[i] = ...
       #endif
        return *this;
    }
    constexpr CLASS_NAME& operator|=(const CLASS_NAME& other) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte |= other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             byte[i] |= other.byte[i];
       #endif
        return *this;
    }
    constexpr CLASS_NAME& operator^=(const STORAGE_TYPE& value) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte ^= value;
       #else
         byte[0] ^= value;
      // for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
      //     byte[i] = ...
       #endif
        return *this;
    }
    constexpr CLASS_NAME& operator^=(const CLASS_NAME& other) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte ^= other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             byte[i] ^= other.byte[i];
       #endif
        return *this;
    }

 // ########################## BIT& ARITHMETIC OPERATORS ########################## //

    constexpr STORAGE_TYPE operator&(const STORAGE_TYPE& value) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte & value;
       #else
         return byte[0] & value;
      // for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
      //     return 0;
       #endif
    }
    constexpr CLASS_NAME operator&(const CLASS_NAME& other) const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = byte & other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             result.byte[i] = byte[i] & other.byte[i];
       #endif
        return result;
    }
    constexpr CLASS_NAME operator|(const STORAGE_TYPE& value) const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = byte | value;
       #else
         result.byte[0] = byte[0] | value;
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             result.byte[i] = byte[i];
       #endif
        return result;
    }
    constexpr CLASS_NAME operator|(const CLASS_NAME& other) const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = byte | other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             result.byte[i] = byte[i] | other.byte[i];
       #endif
        return result;
    }
    constexpr CLASS_NAME operator^(const STORAGE_TYPE& value) const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = byte ^ value;
       #else
         result.byte[0] = byte[0] ^ value;
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             result.byte[i] = byte[i];
       #endif
        return result;
    }
    constexpr CLASS_NAME operator^(const CLASS_NAME& other) const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = byte ^ other.byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             result.byte[i] = byte[i] ^ other.byte[i];
       #endif
        return result;
    }

 // ########################### BITWISE SHIFT OPERATORS ########################### //

    constexpr CLASS_NAME& operator<<=(const INDEX_TYPE& pos) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte <<= pos;
       #else
       const INDEX_TYPE &shl = pos, &shr = STORAGE_BITS - pos;
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
             byte[i] = (byte[i] << shl) | (byte[i-1] >> shr);
         byte[0] <<= shl;
       #endif
        return *this;
    }
    constexpr CLASS_NAME& operator>>=(const INDEX_TYPE& pos) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte >>= pos;
       #else
       const INDEX_TYPE &shr = pos, &shl = STORAGE_BITS - pos;
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH-1; ++i)
             byte[i] = (byte[i] >> shr) | (byte[i+1] << shl);
         byte[ARRAY_LENGTH-1] >>= shr;
       #endif
        return *this;
    }
    constexpr CLASS_NAME operator<<(const INDEX_TYPE& pos) const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = byte << pos;
       #else
       const INDEX_TYPE &shl = pos, &shr = STORAGE_BITS - pos;
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
             result.byte[i] = (byte[i] << shl) | (byte[i-1] >> shr);
         result.byte[0] = byte[0] << shl;
       #endif
        return result;
    }
    constexpr CLASS_NAME operator>>(const INDEX_TYPE& pos) const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = byte >> pos;
       #else
       const INDEX_TYPE &shr = pos, &shl = STORAGE_BITS - pos;
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH-1; ++i)
             result.byte[i] = (byte[i] >> shr) | (byte[i+1] << shl);
         result.byte[ARRAY_LENGTH-1] = byte[ARRAY_LENGTH-1] >> shr;
       #endif
        return result;
    }
    constexpr CLASS_NAME operator~() const noexcept {
        CLASS_NAME result;
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         result.byte = ~byte;
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
             result.byte[i] = ~byte[i];
       #endif
        return result;
    }

 // ########################## ELEMENT ACCESS& MODIFIERS ########################## //

    constexpr CLASS_NAME& set(const INDEX_TYPE& pos) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte |= ((STORAGE_TYPE) 1 << pos);
       #else
         byte[pos / STORAGE_BITS] |= ((STORAGE_TYPE) 1 << (pos % STORAGE_BITS));
       #endif
        return *this;
    }
    constexpr CLASS_NAME& reset(const INDEX_TYPE& pos) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte &= ~((STORAGE_TYPE) 1 << pos);
       #else
         byte[pos / STORAGE_BITS] &= ~((STORAGE_TYPE) 1 << (pos % STORAGE_BITS));
       #endif
        return *this;
    }
    constexpr bool test(const INDEX_TYPE& pos) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte & ((STORAGE_TYPE) 1 << pos);
       #else
         return byte[pos / STORAGE_BITS] & ((STORAGE_TYPE) 1 << (pos % STORAGE_BITS));
       #endif
    }
    INDEX_TYPE popcnt() const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return _popcnt(byte);
       #else
         INDEX_TYPE count = _popcnt(byte[0]);
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             count += _popcnt(byte[i]);
         return count;
       #endif
    }
    INDEX_TYPE tzcnt() const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return _tzcnt(byte);
       #else
         INDEX_TYPE count = _tzcnt(byte[0]);
         for (INDEX_TYPE i = 1; !byte[i-1] && i != ARRAY_LENGTH; ++i)
             count += _tzcnt(byte[i]);
         return count;
       #endif
    }
    CLASS_NAME& pext(const CLASS_NAME& other) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte = _pext(byte, other.byte);
       #else
         INDEX_TYPE shl = 0, shr = STORAGE_BITS, pos = 0; STORAGE_TYPE temp;
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
         {   temp = _pext(byte[i], other.byte[i]); byte[i] = 0;
             if (!shl) byte[pos]  = temp;
                 else  byte[pos] |= temp << shl, byte[pos+1] = temp >> shr;
             shl += _popcnt(other.byte[i]); pos += shl / STORAGE_BITS;
             shl %=  STORAGE_BITS;          shr  = STORAGE_BITS - shl;    }
       #endif
        return *this;
    }
    CLASS_NAME& pdep(const CLASS_NAME& other) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         byte = _pdep(byte, other.byte);
       #else
         INDEX_TYPE shr = 0, shl = STORAGE_BITS, pos = 0; CLASS_NAME temp;
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i)
         {   if (!shr) temp.byte[i] = byte[pos];
                 else  temp.byte[i] = byte[pos] >> shr | byte[pos+1] << shl;
             shr += _popcnt(other.byte[i]); pos += shr / STORAGE_BITS;
             shr %=  STORAGE_BITS;          shl  = STORAGE_BITS - shr;    }
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH; ++i) /* ------- */
             byte[i] = _pdep(temp.byte[i], other.byte[i]);
       #endif
        return *this;
    }

 // ########################## BIT& COMPARISON OPERATORS ########################## //

    constexpr bool operator==(const STORAGE_TYPE& value) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte == value;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
             if (byte[i] != 0) return false;
         return byte[0] == value;
       #endif
    }
    constexpr bool operator==(const CLASS_NAME& other) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte == other.byte;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
             if (byte[i] != other.byte[i]) return false;
         return byte[0] == other.byte[0];
       #endif
    }
    constexpr bool operator!=(const STORAGE_TYPE& value) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte != value;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
             if (byte[i] != 0) return true;
         return byte[0] != value;
       #endif
    }
    constexpr bool operator!=(const CLASS_NAME& other) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte != other.byte;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
             if (byte[i] != other.byte[i]) return true;
         return byte[0] != other.byte[0];
       #endif
    }

 #if defined(LEX_INTEGER_COMPARATORS)
    constexpr bool operator<(const CLASS_NAME& other) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte < other.byte;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
                if (byte[i] < other.byte[i]) return true;
           else if (byte[i] > other.byte[i]) return false;
         return byte[0] < other.byte[0];
       #endif
    }
    constexpr bool operator<=(const CLASS_NAME& other) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte <= other.byte;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
                if (byte[i] < other.byte[i]) return true;
           else if (byte[i] > other.byte[i]) return false;
         return byte[0] <= other.byte[0];
       #endif
    }
    constexpr bool operator>(const CLASS_NAME& other) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte > other.byte;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
                if (byte[i] > other.byte[i]) return true;
           else if (byte[i] < other.byte[i]) return false;
         return byte[0] > other.byte[0];
       #endif
    }
    constexpr bool operator>=(const CLASS_NAME& other) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte >= other.byte;
       #else
         for (INDEX_TYPE i = ARRAY_LENGTH-1; i != 0; --i)
                if (byte[i] > other.byte[i]) return true;
           else if (byte[i] < other.byte[i]) return false;
         return byte[0] >= other.byte[0];
       #endif
    }
    #undef LEX_INTEGER_COMPARATORS
 #endif

 // ########################## SET& COMPARISON OPERATORS ########################## //

 #if defined(SET_ELEMENT_COMPARATORS)
   _constexpr bool operator<(const CLASS_NAME& other) const noexcept {
        if (*this == other) return false;
       { INDEX_TYPE x = (*this).popcnt(); INDEX_TYPE y = other.popcnt(); if (x != y)  return x < y; }
       { INDEX_TYPE x = (*this).tzcnt();  INDEX_TYPE y = other.tzcnt();  if (x != y)  return x < y;
         CLASS_NAME X = (*this);          CLASS_NAME Y = other;
               do { X.reset(x);    x = X.tzcnt();    Y.reset(y);    y = Y.tzcnt(); }
                                                                      while (x == y); return x < y; }
    }
   _constexpr bool operator<=(const CLASS_NAME& other) const noexcept {
        if (*this == other) return true;
       { INDEX_TYPE x = (*this).popcnt(); INDEX_TYPE y = other.popcnt(); if (x != y)  return x < y; }
       { INDEX_TYPE x = (*this).tzcnt();  INDEX_TYPE y = other.tzcnt();  if (x != y)  return x < y;
         CLASS_NAME X = (*this);          CLASS_NAME Y = other;
               do { X.reset(x);    x = X.tzcnt();    Y.reset(y);    y = Y.tzcnt(); }
                                                                      while (x == y); return x < y; }
    }
   _constexpr bool operator>(const CLASS_NAME& other) const noexcept {
        if (*this == other) return false;
       { INDEX_TYPE x = (*this).popcnt(); INDEX_TYPE y = other.popcnt(); if (x != y)  return x > y; }
       { INDEX_TYPE x = (*this).tzcnt();  INDEX_TYPE y = other.tzcnt();  if (x != y)  return x > y;
         CLASS_NAME X = (*this);          CLASS_NAME Y = other;
               do { X.reset(x);    x = X.tzcnt();    Y.reset(y);    y = Y.tzcnt(); }
                                                                      while (x == y); return x > y; }
    }
   _constexpr bool operator>=(const CLASS_NAME& other) const noexcept {
        if (*this == other) return true;
       { INDEX_TYPE x = (*this).popcnt(); INDEX_TYPE y = other.popcnt(); if (x != y)  return x > y; }
       { INDEX_TYPE x = (*this).tzcnt();  INDEX_TYPE y = other.tzcnt();  if (x != y)  return x > y;
         CLASS_NAME X = (*this);          CLASS_NAME Y = other;
               do { X.reset(x);    x = X.tzcnt();    Y.reset(y);    y = Y.tzcnt(); }
                                                                      while (x == y); return x > y; }
    }
    static constexpr bool disjoint(const CLASS_NAME& x, const CLASS_NAME& y) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return !(x.byte & y.byte);
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH-1; ++i)
             if (x.byte[i] & y.byte[i]) return false;
         return !(x.byte[ARRAY_LENGTH-1] & y.byte[ARRAY_LENGTH-1]);
       #endif
    }
    static constexpr bool disjoint(const CLASS_NAME& x, const CLASS_NAME& y, const CLASS_NAME& z) noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return !(x.byte & y.byte & z.byte);
       #else
         for (INDEX_TYPE i = 0; i != ARRAY_LENGTH-1; ++i)
             if (x.byte[i] & y.byte[i] & z.byte[i]) return false;
         return !(x.byte[ARRAY_LENGTH-1] & y.byte[ARRAY_LENGTH-1] & z.byte[ARRAY_LENGTH-1]);
       #endif
    }
    #undef SET_ELEMENT_COMPARATORS
 #endif

// ############################# END CLASS DEFINITION ############################# //

constexpr STORAGE_TYPE operator%(const STORAGE_TYPE& value) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return byte % value;
       #else
         STORAGE_TYPE hash = byte[0];
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             hash ^= byte[i];
         return hash % value;
       #endif
    }
};

template<> struct std::hash<CLASS_NAME> {
    constexpr STORAGE_TYPE operator()(const CLASS_NAME& obj) const noexcept {
       #if BIT_LENGTH <= MAX_STORAGE_BITS
         return obj.byte;
       #else
         STORAGE_TYPE hash = obj.byte[0];
         for (INDEX_TYPE i = 1; i != ARRAY_LENGTH; ++i)
             hash ^= obj.byte[i];
         return hash;
       #endif
    }
};

#undef CLASS_NAME
#undef STORAGE_TYPE
#undef INDEX_TYPE
#undef BIT_LENGTH

#undef STORAGE_BITS
#undef MAX_STORAGE_BITS
#undef ARRAY_LENGTH

#undef _popcnt
#undef _tzcnt
#undef _pext
#undef _pdep
