/*************

ans_fast.cpp

by cbloom

v2 02-21-2014

-------------------------

fast implementation of tANS and rANS

intended to demonstrate how fast versions of ANS look
and measure performance

-------------------------

please link to this :

http://cbloomrants.blogspot.com/2014/02/02-18-14-understanding-ans-conclusion.html

*******************/

#include "cblib/FileEnum.h"
#include "cblib/FileUtil.h"
#include "cblib/Log.h"
#include "cblib/LogUtil.h"
#include "cblib/Rand.h"
#include "cblib/Util.h"
#include "cblib/Win32Util.h"
#include "cblib/inc.h"
#include "cblib/vector.h"
#include "cblib/vector_s.h"

#define MAX_L 4096 // to avoid some allocations
#define MAX_ALPHABET 256

//-----------------------------------------------------

#define NOINLINE __declspec(noinline)

extern "C" unsigned long __cdecl _byteswap_ulong(unsigned long _Long);
#pragma intrinsic(_byteswap_ulong)

extern "C" unsigned __int64 __cdecl _byteswap_uint64(unsigned __int64 val);
#pragma intrinsic(_byteswap_uint64)

extern "C" unsigned long __cdecl _lrotl(unsigned long, int);
#pragma intrinsic(_lrotl)

extern "C" unsigned __int64 __cdecl _rotl64(unsigned __int64 _Val, int _Shift);
#pragma intrinsic(_rotl64)

extern "C" unsigned long __cdecl _lrotr(unsigned long, int);
#pragma intrinsic(_lrotr)

extern "C" unsigned __int64 __cdecl _rotr64(unsigned __int64 _Val, int _Shift);
#pragma intrinsic(_rotr64)

extern "C" unsigned __int64 __umulh(unsigned __int64 a, unsigned __int64 b);
#pragma intrinsic(__umulh)

extern "C" unsigned char _BitScanReverse(
    unsigned long* Index, unsigned long Mask);
#pragma intrinsic(_BitScanReverse)

static inline int bsr(unsigned long val)
{
    ASSERT(val != 0);
    unsigned long b = 0;
    _BitScanReverse(&b, val);
    return (int)b;
}

// ilog2ceil = ceil(log2(val))
static inline int ilog2ceil(unsigned long val)
{
    ASSERT(val != 0);
    if (val == 1)
        return 0;
    else
        return bsr(val - 1) + 1;
}

static inline int ilog2floor(unsigned long val)
{
    if (cb::IsPow2(val))
        return ilog2ceil(val);
    else
        return ilog2ceil(val) - 1;
}

//-----------------------------------------------------

USE_CB

struct sort_sym {
    int sym;
    float rank;
    sort_sym() {}
    sort_sym(int s, float r)
        : sym(s)
        , rank(r)
    {
    }
    bool operator<(const sort_sym& rhs) const { return rank < rhs.rank; }
};

static NOINLINE void normalize_counts(
    uint32* to, const uint32* from, int alphabet, int to_sum_desired)
{
    uint32 from_sum = cb::sum(from, from + alphabet);
    ASSERT_RELEASE(from_sum > 0);

    double scale = (double)to_sum_desired / from_sum;

    uint32 to_sum = 0;
        for
            LOOP(i, alphabet)
            {
                if (from[i] == 0) {
                    to[i] = 0;
                } else {
                    double from_scaled = from[i] * scale;
                    uint32 down = (uint32)(from_scaled);

                    to[i] = (from_scaled * from_scaled < down * (down + 1))
                        ? down
                        : down + 1;

                    ASSERT(to[i] > 0);
                    to_sum += to[i];
                }
            }

        int32 correction = to_sum_desired - to_sum;
        if (correction != 0) {
            // lprintfvar(correction);
            int32 correction_sign = (correction > 0) ? 1 : -1;

            vector_s<sort_sym, MAX_ALPHABET> heap;
            heap.reserve(alphabet);

                for
                    LOOP(i, alphabet)
                    {
                        if (from[i] == 0)
                            continue;
                        ASSERT(to[i] != 0);
                        if (to[i] > 1 || correction_sign == 1) {
                            float change
                                = (float)logf(
                                      1.0f + correction_sign * (1.0f / to[i]))
                                * from[i];

                            heap.push_back(sort_sym(i, change));
                        }
                    }

                std::make_heap(heap.begin(), heap.end());

                while (correction != 0) {
                    ASSERT_RELEASE(!heap.empty());
                    std::pop_heap(heap.begin(), heap.end());
                    sort_sym ss = heap.back();
                    heap.pop_back();

                    int i = ss.sym;
                    ASSERT(from[i] != 0);

                    to[i] += correction_sign;
                    correction -= correction_sign;
                    ASSERT(to[i] != 0);

                    if (to[i] > 1 || correction_sign == 1) {
                        float change = (float)logf(1.0f
                                           + correction_sign * (1.0f / to[i]))
                            * from[i];

                        heap.push_back(sort_sym(i, change));
                        std::push_heap(heap.begin(), heap.end());
                    }
                }
        }

        ASSERT(cb::sum(to, to + alphabet) == (uint32)to_sum_desired);
}

static double entropy(const uint32* counts, int alphabet)
{
    double H = 0;
    uint32 countsum = 0;
        for
            LOOP(a, alphabet)
            {
                uint32 c = counts[a];
                if (c == 0)
                    continue;
                countsum += c;
                H -= c * log2(c);
            }
        H = log2(countsum) + H / countsum;
        return H;
}

//=========================

struct tans_encode_table {
    int L; // renorm is to [L,2L-1]
    int state_count, state_bits, alphabet;

    vector_s<uint16, MAX_L> encode_table_packed; // encode_table_packed[ x +
    // packed_table_offset[s] ]

    // 5 byte encode entry
    /*
    #pragma pack(push)
    #pragma pack(1) // surprisingly does not hurt
    struct encode_entry { int16 packed_table_offset; uint16 max_state_thresh;
    uint8 max_state_numbits; };
    #pragma pack(pop)
    */
    // 8-byte encode entry with packed_table_ptr
    struct encode_entry {
        uint16* packed_table_ptr;
        uint16 max_state_thresh;
        uint8 max_state_numbits;
    };
    vector_s<encode_entry, MAX_ALPHABET> encode_sym_table; // index by symbol
};

struct tans_decode_table {
    int L; // renorm is to [L,2L-1]
    int state_count, state_bits, alphabet;

// 4 byte decode_entry :
// speed the same no matter what
// struct decode_entry { uint16 next_state; uint8 num_bits; uint8 sym; };
#pragma pack(push)
#pragma pack(1)
    struct decode_entry {
        uint8 sym;
        uint8 num_bits;
        uint16 next_state;
    };
#pragma pack(pop)

    vector_s<decode_entry, MAX_L> decode_table_packed; // index by [state - L]
};

//===================================

// radix rank is a fixed point with radix_fpbits bits of precision

/*

// more precise radix
// 1791466.75
#define RADIX_BITS(L) L
enum { radix_fpbits = 0 };

/*/

// faster
// 1791581.13
#define RADIX_BITS(L) 8
enum { radix_fpbits = 8 };

/**/

//===================================

static NOINLINE void tans_make_encode_table(tans_encode_table* tables, int L,
    int alphabet, const uint32* normalized_counts)
{
    ASSERT(
        cb::sum(normalized_counts, normalized_counts + alphabet) == (uint32)L);

    tables->L = L;
    tables->state_count = 2 * L;
    tables->alphabet = alphabet;

    int L_bits = ilog2ceil(L);
    // current construction requires L a power of 2 :
    //  in general that is not required
    ASSERT_RELEASE((1 << L_bits) == L);

    tables->state_bits = L_bits + 1;

    tables->encode_table_packed.resize(L, 0);
    tables->encode_sym_table.resize(alphabet);

    //===================================================
    // radix sort

    ASSERT((int)cb::sum(normalized_counts, normalized_counts + alphabet) == L);

    int radix_bits = RADIX_BITS(L_bits);
    int radix_size = 1 << radix_bits;

    uint32 radix_numer = ((radix_size - 1) << radix_fpbits);
    // uint32 radix_numer = (radix_size<<radix_fpbits)-1;

    vector_s<uint32, MAX_L + 1> radix_histo;
    radix_histo.resize(radix_size, 0);

    // count ranks to histo :
        for
            LOOP(sym, alphabet)
            {
                uint32 count = normalized_counts[sym];
                if (count == 0)
                    continue;

                uint32 invp = radix_numer / count;
                ASSERT_RELEASE(invp > 0 || count == (uint32)L);

                uint32 rank = invp;

                for
                    LOOP(c, (int)count)
                    {
                        radix_histo[(rank >> radix_fpbits)]++;
                        rank += invp;
                    }
            }

        // sum radix histo to cumulative prob :
        int radix_cum = 0;
        for
            LOOP(r, radix_size)
            {
                uint32 cur = radix_histo[r];
                radix_histo[r] = radix_cum;
                radix_cum += cur;
            }

        ASSERT_RELEASE(radix_cum == L);

        // radix_histo[r] now tells you where that radix goes in the array

        // read out of radix and build encode table at the same time :
        uint32 cumulative_count = 0;
        for
            LOOP(sym, alphabet)
            {
                uint32 count = normalized_counts[sym];
                if (count == 0) {
                    tables->encode_sym_table[sym].packed_table_ptr = NULL;
                } else if (count
                    == 1) // helps speed a tiny bit to special case 1
                {
                    int num_bits = L_bits;

                    tans_encode_table::encode_entry& ee
                        = tables->encode_sym_table[sym];

                    ee.max_state_numbits = check_value_cast<uint8>(num_bits);
                    ee.max_state_thresh
                        = check_value_cast<uint16>(2 << num_bits);

                    ee.packed_table_ptr = tables->encode_table_packed.data()
                        + check_value_cast<int16>(cumulative_count - count);

                    cumulative_count += 1;

                    uint32 index = radix_numer >> radix_fpbits;

                    uint32 to = radix_histo[index];
                    radix_histo[index]++;

                    int to_state = L + to;

                    // int from_state = count + c;
                    // packed_table_ptr[ from_state ] =
                    // check_value_cast<uint16>( to_state );

                    ee.packed_table_ptr[count]
                        = check_value_cast<uint16>(to_state);
                } else {
                    int num_bits = L_bits - ilog2ceil(count);
                    ASSERT(num_bits == ilog2floor(L / count));
                    ASSERT(num_bits == L_bits - bsr(2 * count - 1));

                    tans_encode_table::encode_entry& ee
                        = tables->encode_sym_table[sym];

                    ee.max_state_numbits = check_value_cast<uint8>(num_bits);
                    ee.max_state_thresh
                        = check_value_cast<uint16>((count + count) << num_bits);

                    // + cumulative_count to pack us in the array
                    // - nc because that is the lowest state that will index
                    // this array
                    // tables->encode_sym_table[a].packed_table_offset =
                    // check_value_cast<int16>(cumulative_count - nc);
                    ee.packed_table_ptr = tables->encode_table_packed.data()
                        + check_value_cast<int16>(cumulative_count - count);

                    cumulative_count += count;

                    uint32 invp = radix_numer / count;
                    uint32 rank = invp;

                    uint16* packed_table_ptr = ee.packed_table_ptr;
                    packed_table_ptr += count;

                        for
                            LOOP(c, (int)count)
                            {
                                uint32 index = rank >> radix_fpbits;
                                rank += invp;

                                uint32 to = radix_histo[index];
                                radix_histo[index]++;

                                int to_state = L + to;

                                // int from_state = count + c;
                                // packed_table_ptr[ from_state ] =
                                // check_value_cast<uint16>( to_state );

                                *packed_table_ptr++
                                    = check_value_cast<uint16>(to_state);
                            }
                }
            }

        ASSERT_RELEASE(cumulative_count == (uint32)L);
}

static NOINLINE void tans_make_decode_table(tans_decode_table* tables, int L,
    int alphabet, const uint32* normalized_counts)
{
    ASSERT_RELEASE(alphabet <= MAX_ALPHABET);
    ASSERT(
        cb::sum(normalized_counts, normalized_counts + alphabet) == (uint32)L);

    tables->L = L;
    tables->state_count = 2 * L;
    tables->alphabet = alphabet;

    int L_bits = ilog2ceil(L);
    // current construction requires L a power of 2 :
    ASSERT_RELEASE((1 << L_bits) == L);

    tables->state_bits = L_bits + 1;

    tables->decode_table_packed.resize(L);

    //=======================================================================

    int radix_bits = RADIX_BITS(L_bits);
    int radix_size = 1 << radix_bits;

    uint32 radix_numer = ((radix_size - 1) << radix_fpbits);
    // uint32 radix_numer = (radix_size<<radix_fpbits)-1;

    vector_s<uint32, MAX_L + 1> radix_histo;
    radix_histo.resize(radix_size, 0);

    // count ranks to histo :
        for
            LOOP(sym, alphabet)
            {
                uint32 count = normalized_counts[sym];
                if (count == 0)
                    continue;

                uint32 invp = radix_numer / count;
                ASSERT_RELEASE(invp > 0 || count == (uint32)L);

                uint32 rank = invp;

                for
                    LOOP(c, (int)count)
                    {
                        radix_histo[(rank >> radix_fpbits)]++;
                        rank += invp;
                    }
            }

        // sum radix histo to cumulative prob :
        int radix_cum = 0;
        for
            LOOP(r, radix_size)
            {
                uint32 cur = radix_histo[r];
                radix_histo[r] = radix_cum;
                radix_cum += cur;
            }

        ASSERT_RELEASE(radix_cum == L);

        // radix_histo[r] now tells you where that radix goes in the array

        tans_decode_table::decode_entry* detable
            = tables->decode_table_packed.data();

        // count ranks to histo :
        for
            LOOP(sym, alphabet)
            {
                uint32 count = normalized_counts[sym];
                if (count == 0) {
                    // nop
                } else if (count == 1) // helps speed a little to special case 1
                {
                    uint32 rank = radix_numer;

                    uint32 index = rank >> radix_fpbits;
                    uint32 to = radix_histo[index];
                    radix_histo[index]++;

                    // int from_state = 1;

                    // decoder numbits :
                    ASSERT(bsr(1) == 0);
                    // int num_bits = L_bits - bsr(from_state);
                    int num_bits = L_bits;

                    ASSERT(detable[to].next_state == 0);
                    detable[to].next_state
                        = check_value_cast<uint16>(1 << num_bits);
                    detable[to].num_bits = check_value_cast<uint8>(num_bits);
                    detable[to].sym = check_value_cast<uint8>(sym);
                } else {
                    uint32 invp = radix_numer / count;
                    uint32 rank = invp;

                    // simpler, faster !

                        for
                            LOOP(c, (int)count)
                            {
                                uint32 index = rank >> radix_fpbits;
                                uint32 to = radix_histo[index];
                                radix_histo[index]++;

                                int from_state = count + c;

                                // decoder numbits :
                                int num_bits = L_bits - bsr(from_state);

                                ASSERT(detable[to].next_state == 0);
                                detable[to].sym = check_value_cast<uint8>(sym);
                                detable[to].num_bits
                                    = check_value_cast<uint8>(num_bits);
                                detable[to].next_state
                                    = check_value_cast<uint16>(
                                        from_state << num_bits);

                                rank += invp;
                            }
                }
            }
}

//===============================================================================

#define BITOUT_VARS(bout_bits, bout_numbits, bout_ptr)                         \
    uint64 bout_bits;                                                          \
    int64 bout_numbits;                                                        \
    uint8* bout_ptr;

#define BITOUT_START(bout_bits, bout_numbits, bout_ptr, end)                   \
    do {                                                                       \
        bout_bits = 0;                                                         \
        bout_numbits = 0;                                                      \
        bout_ptr = (uint8*)end;                                                \
    } while (0)

// to left justify the new val :
// this way does the masking for you :
// bout_bits |= ((uint64)val) << (64 - nb);

// this way is faster but leaves crud in the low bits
// @@ this is only okay if rotated val does not run into numbits
// bout_bits |= _rotr64((uint64)val,nb); \

#define BITOUT_PUT(bout_bits, bout_numbits, val, nb)                           \
    do {                                                                       \
        ASSERT((bout_numbits + nb) <= 64);                                     \
        bout_bits >>= nb;                                                      \
        bout_bits |= _rotr64((uint64)val, (int)nb);                            \
        bout_numbits += nb;                                                    \
    } while (0)

#define BITOUT_FLUSH(bout_bits, bout_numbits, bout_ptr)                        \
    do {                                                                       \
        *((uint64*)(bout_ptr - 8))                                             \
            = _byteswap_uint64(bout_bits >> (64 - bout_numbits));              \
        bout_ptr -= bout_numbits >> 3;                                         \
        bout_numbits &= 7;                                                     \
    } while (0)

#define BITOUT_END(bout_bits, bout_numbits, bout_ptr)                          \
    do {                                                                       \
        while (bout_numbits > 8) {                                             \
            bout_ptr--;                                                        \
            *bout_ptr = (uint8)(bout_bits >> (64 - bout_numbits));             \
            bout_numbits -= 8;                                                 \
        }                                                                      \
        if (bout_numbits == 0) {                                               \
            bout_numbits = 8;                                                  \
        } else {                                                               \
            ASSERT(bout_numbits >= 1 && bout_numbits <= 8);                    \
            bout_ptr--;                                                        \
            *bout_ptr = (uint8)(bout_bits >> 56);                              \
        }                                                                      \
        ASSERT(bout_numbits >= 1 && bout_numbits <= 8);                        \
    } while (0)

//===============================================================
/*

tans_encode

encodes a buffer
writes compressed data *back* from "to_end"
returns last pointer and last # of bits

*/

static NOINLINE std::pair<void*, int> tans_encode(const uint8* buffer,
    ptrdiff_t length, void* to_end, const tans_encode_table* table)
{
    const tans_encode_table::encode_entry* eetable
        = table->encode_sym_table.data();

    uint64 state = table->L;

    BITOUT_VARS(bout_bits, bout_numbits, bout_ptr);

    BITOUT_START(bout_bits, bout_numbits, bout_ptr, to_end);

#define ENCODE_ONE(state)                                                      \
    do {                                                                       \
        sym = *bufptr--;                                                       \
        ee = eetable + sym;                                                    \
        msnb = ee->max_state_numbits;                                          \
        msnb += (state >= ee->max_state_thresh);                               \
        BITOUT_PUT(bout_bits, bout_numbits, state, msnb);                      \
        state = ee->packed_table_ptr[state >> msnb];                           \
    } while (0)

    const uint8* bufptr = buffer + length - 1;

    for (ptrdiff_t i = 0; i < (length % 4); i++) {
        int sym;
        const tans_encode_table::encode_entry* ee;
        int64 msnb;

        ENCODE_ONE(state);
    }

    BITOUT_FLUSH(bout_bits, bout_numbits, bout_ptr);

    for (ptrdiff_t i = 0; i < (length / 4); i++) {
        int sym;
        const tans_encode_table::encode_entry* ee;
        int64 msnb;

        // can unroll four encodes per flush :

        ENCODE_ONE(state);
        ENCODE_ONE(state);
        ENCODE_ONE(state);
        ENCODE_ONE(state);

        BITOUT_FLUSH(bout_bits, bout_numbits, bout_ptr);
    }

    BITOUT_PUT(bout_bits, bout_numbits, state, table->state_bits);
    BITOUT_END(bout_bits, bout_numbits, bout_ptr);

    return std::pair<void*, int>(bout_ptr, (int)bout_numbits);
}

static NOINLINE std::pair<void*, int> tans_encode_int2(const uint8* buffer,
    ptrdiff_t length, void* to_end, const tans_encode_table* table)
{
    const tans_encode_table::encode_entry* eetable
        = table->encode_sym_table.data();

    uint64 state1 = table->L;
    uint64 state2 = table->L;

    BITOUT_VARS(bout_bits, bout_numbits, bout_ptr);

    BITOUT_START(bout_bits, bout_numbits, bout_ptr, to_end);

    const uint8* bufptr = buffer + length - 1;

    for (ptrdiff_t i = 0; i < (length % 4); i++) {
        int sym;
        const tans_encode_table::encode_entry* ee;
        int64 msnb;

        ENCODE_ONE(state1);
    }

    BITOUT_FLUSH(bout_bits, bout_numbits, bout_ptr);

    for (ptrdiff_t i = 0; i < (length / 4); i++) {
        int sym;
        const tans_encode_table::encode_entry* ee;
        int64 msnb;

        // can unroll four encodes per flush :

        ENCODE_ONE(state1);
        ENCODE_ONE(state2);
        ENCODE_ONE(state1);
        ENCODE_ONE(state2);

        BITOUT_FLUSH(bout_bits, bout_numbits, bout_ptr);
    }

    BITOUT_PUT(bout_bits, bout_numbits, state1, table->state_bits);
    BITOUT_PUT(bout_bits, bout_numbits, state2, table->state_bits);

    BITOUT_END(bout_bits, bout_numbits, bout_ptr);

    return std::pair<void*, int>(bout_ptr, (int)bout_numbits);
}

#define BITIN_VARS(bitin_bits, bitin_numbits, bitin_ptr)                       \
    uint64 bitin_bits;                                                         \
    int64 bitin_numbits;                                                       \
    uint8* bitin_ptr;

#define BITIN_START(bitin_bits, bitin_numbits, bitin_ptr, begin_ptr, lastbits) \
    do {                                                                       \
        bitin_ptr = (uint8*)begin_ptr;                                         \
        bitin_bits = *bitin_ptr++;                                             \
        bitin_bits <<= 56;                                                     \
        ASSERT(lastbits >= 1 && lastbits <= 8);                                \
        bitin_numbits = lastbits;                                              \
    } while (0)

// bitin_numbits == 64 is rare but possible
//	it messes up the shift >>bitin_numbits
// can either check for it or handle it
// checking seems to be faster
// bitin_bits |= next8 >> bitin_numbits;
/*
#define BITIN_REFILL(bitin_bits,bitin_numbits,bitin_ptr) do { \
                ASSERT( bitin_numbits > 0 && bitin_numbits <= 64 ); \
                int64 bytesToGet = (64 - bitin_numbits)>>3; \
                uint64 next8 = _byteswap_uint64( *( (uint64 *)bitin_ptr ) ); \
                bitin_ptr += bytesToGet; \
                bitin_bits |= (next8 >> 1) >> (bitin_numbits-1); \
                bitin_numbits += bytesToGet<<3; \
                ASSERT( bitin_numbits >= 56 && bitin_numbits <= 64 ); \
        } while(0)
/*/
#define BITIN_REFILL(bitin_bits, bitin_numbits, bitin_ptr)                     \
    do {                                                                       \
        if (bitin_numbits < 64) {                                              \
            ASSERT(bitin_numbits > 0 && bitin_numbits <= 64);                  \
            uint64 next8 = _byteswap_uint64(*((uint64*)bitin_ptr));            \
            int64 bytesToGet = (64 - bitin_numbits) >> 3;                      \
            bitin_ptr += bytesToGet;                                           \
            bitin_bits |= next8 >> bitin_numbits;                              \
            bitin_numbits += bytesToGet << 3;                                  \
            ASSERT(bitin_numbits >= 56 && bitin_numbits <= 64);                \
        }                                                                      \
    } while (0)
/**/

// nb==0 is possible but very rare :
// most files don't generate any nb==0 transitions at all
// which is faster depends on the file
// files like "pic" where nb == 0 is common, the if is slower cuz its
// unpredictable
// files like "book1" where the if is super predictable, it's faster (379 vs 368
// mbps)
//*
#define BITIN_OR(bitin_bits, bitin_numbits, nb, ret)                           \
    do {                                                                       \
        ASSERT(nb <= bitin_numbits);                                           \
        ret |= (bitin_bits >> 1) >> (63 - nb);                                 \
        bitin_bits <<= nb;                                                     \
        bitin_numbits -= nb;                                                   \
    } while (0)
/*/
#define BITIN_OR(bitin_bits,bitin_numbits,nb,ret) if ( nb > 0 ) { \
                ASSERT( nb <= bitin_numbits ); \
                ret |= bitin_bits >> (64 - nb); \
                bitin_bits <<= nb; \
                bitin_numbits -= nb; \
        }
/**/

static NOINLINE void tans_decode(void* comp_ptr, int last_nbits, uint8* tobuf,
    ptrdiff_t tobuflen, const tans_decode_table* table)
{
    BITIN_VARS(bitin_bits, bitin_numbits, bitin_ptr);
    BITIN_START(bitin_bits, bitin_numbits, bitin_ptr, comp_ptr, last_nbits);
    BITIN_REFILL(bitin_bits, bitin_numbits, bitin_ptr);

    uint64 state = 0;
    BITIN_OR(bitin_bits, bitin_numbits, table->state_bits, state);

    const tans_decode_table::decode_entry* detable
        = table->decode_table_packed.data() - table->L;

    uint8* toptr = tobuf;

#define DECODE_ONE(state)                                                      \
    do {                                                                       \
        de = detable + state;                                                  \
        nb = de->num_bits;                                                     \
        state = de->next_state;                                                \
        BITIN_OR(bitin_bits, bitin_numbits, nb, state);                        \
        *toptr++ = (uint8)de->sym;                                             \
    } while (0)

    for (ptrdiff_t i = 0; i < (tobuflen % 4); i++) {
        ASSERT(state >= table->L && state < 2 * table->L);

        const tans_decode_table::decode_entry* de;
        int64 nb;

        DECODE_ONE(state);
    }

    BITIN_REFILL(bitin_bits, bitin_numbits, bitin_ptr);

    for (ptrdiff_t i = 0; i < (tobuflen / 4); i++) {
        ASSERT(state >= table->L && state < 2 * table->L);

        const tans_decode_table::decode_entry* de;
        int64 nb;

        DECODE_ONE(state);
        DECODE_ONE(state);
        DECODE_ONE(state);
        DECODE_ONE(state);

        BITIN_REFILL(bitin_bits, bitin_numbits, bitin_ptr);
    }

    ASSERT(state == table->L);
}

static NOINLINE void tans_decode_int2(void* comp_ptr, int last_nbits,
    uint8* tobuf, ptrdiff_t tobuflen, const tans_decode_table* table)
{
    BITIN_VARS(bitin_bits, bitin_numbits, bitin_ptr);
    BITIN_START(bitin_bits, bitin_numbits, bitin_ptr, comp_ptr, last_nbits);
    BITIN_REFILL(bitin_bits, bitin_numbits, bitin_ptr);

    uint64 state1 = 0;
    uint64 state2 = 0;
    BITIN_OR(bitin_bits, bitin_numbits, table->state_bits, state2);
    BITIN_OR(bitin_bits, bitin_numbits, table->state_bits, state1);

    const tans_decode_table::decode_entry* detable
        = table->decode_table_packed.data() - table->L;

    uint8* toptr = tobuf;

    BITIN_REFILL(bitin_bits, bitin_numbits, bitin_ptr);

    for (ptrdiff_t i = 0; i < (tobuflen / 4); i++) {
        const tans_decode_table::decode_entry* de;
        int64 nb;

        DECODE_ONE(state2);
        DECODE_ONE(state1);
        DECODE_ONE(state2);
        DECODE_ONE(state1);

        BITIN_REFILL(bitin_bits, bitin_numbits, bitin_ptr);
    }

    for (ptrdiff_t i = 0; i < (tobuflen % 4); i++) {
        const tans_decode_table::decode_entry* de;
        int64 nb;

        DECODE_ONE(state1);
    }

    ASSERT(state1 == table->L);
    ASSERT(state2 == table->L);
}

//================================================================

void log_tables(const tans_encode_table* entable,
    const tans_decode_table* detable, const uint32* normalized_counts)
{
    int L = entable->L;

        for
            LOOP(k, L)
            {
                const tans_decode_table::decode_entry& de
                    = detable->decode_table_packed[k];
                lprintf("%c", de.sym);
            }
        lprintf("\n");

        lprintf("decode:\n");
        for
            LOOP(k, L)
            {
                const tans_decode_table::decode_entry& de
                    = detable->decode_table_packed[k];
                lprintf("%2d: %c -> %2d + (%d bits) = [%2d,%2d]\n", k, de.sym,
                    de.next_state - L, de.num_bits, de.next_state - L,
                    de.next_state - L + (1 << de.num_bits) - 1);
            }

        lprintf("encode:\n");
        for
            LOOP(a, detable->alphabet)
            {
                if (normalized_counts[a] == 0)
                    continue;
                const tans_encode_table::encode_entry& ee
                    = entable->encode_sym_table[a];
                int nb = ee.max_state_numbits;
                int thresh = ee.max_state_thresh - L;
                int count = normalized_counts[a];

                lprintf("%c : b=%d+(t>=%d) : {", a, nb, thresh);
                for
                    LOOP(i, count)
                    {
                        int s = ee.packed_table_ptr[count + i] - L;
                        lprintf("%d", s);
                        if (i < count - 1)
                            lprintf(",");
                    }
                lprintf("}\n");
            }
}

//================================================================
/*

rANS for comparison

no attempt to make encoding fast

64-bit state
32-bit renormalization

*/

#define RANS64_RENORM_LO_SHIFT (31)
#define RANS64_RENORM_LO (1ULL << RANS64_RENORM_LO_SHIFT)
#define RANS64_RENORM_HI_SHIFT (63)
#define RANS64_RENORM_HI (1ULL << RANS64_RENORM_HI_SHIFT)

struct rans_encode_table {
    struct entry {
        uint64 rcp_freq; // Fixed-point reciprocal frequency
        uint16 bias; // Bias
        uint16 M_minus_freq; // M - freq
        uint16 freq;
        uint16 rcp_shift; // Reciprocal shift
    };

    entry table[MAX_ALPHABET];

    int cumprobtot_bits;
};

void rans_make_encode_table(rans_encode_table* entable,
    const uint32* normalized_counts, int M_bits, int alphabet)
{
    ASSERT(cb::sum(normalized_counts, normalized_counts + alphabet)
        == (1UL << M_bits));

    entable->cumprobtot_bits = M_bits;

    uint32 M = 1U << M_bits;

    uint32 low = 0;
        for
            LOOP(a, alphabet)
            {
                uint32 freq = normalized_counts[a];
                if (freq == 0)
                    continue;

                /*

                reciprocals from multiplication

                https://github.com/rygorous/ryg_rans

                */

                rans_encode_table::entry& e = entable->table[a];

                e.freq = cb::check_value_cast<uint16>(freq);
                e.M_minus_freq = cb::check_value_cast<uint16>(M - freq);
                if (freq < 2) {
                    // freq=0 symbols are never valid to encode, so it doesn't
                    // matter what
                    // we set our values to.
                    //
                    // freq=1 is tricky, since the reciprocal of 1 is 1;
                    // unfortunately,
                    // our fixed-point reciprocal approximation can only
                    // multiply by values
                    // smaller than 1.
                    //
                    // So we use the "next best thing": rcp_freq=~0,
                    // rcp_shift=0.
                    // This gives:
                    //   q = mul_hi(x, rcp_freq) >> rcp_shift
                    //     = mul_hi(x, (1<<64) - 1)) >> 0
                    //     = floor(x - x/(2^64))
                    //     = x - 1 if 1 <= x < 2^64
                    // and we know that x>0 (x=0 is never in a valid
                    // normalization interval).
                    //
                    // So we now need to choose the other parameters such that
                    //   x_new = x*M + start
                    // plug it in:
                    //     x*M + start                   (desired result)
                    //   = bias + x + q*M_minus_freq        (*)
                    //   = bias + x + (x - 1)*(M - 1)    (plug in q=x-1,
                    //   M_minus_freq)
                    //   = bias + 1 + (x - 1)*M
                    //   = x*M + (bias + 1 - M)
                    //
                    // so we have start = bias + 1 - M, or equivalently
                    //   bias = start + M - 1.
                    e.rcp_freq = ~0ull;
                    e.rcp_shift = 0;
                    e.bias = cb::check_value_cast<uint16>(low + M - 1);
                } else {
                    // Alverson, "Integer Division using reciprocals"
                    // shift=ceil(log2(freq))
                    uint32 shift = ilog2ceil(freq);
                    uint64 x0, x1, t0, t1;

                    // long divide ((uint128) (1 << (shift + 63)) + freq-1) /
                    // freq
                    // by splitting it into two 64:64 bit divides (this works
                    // because
                    // the dividend has a simple form.)
                    x0 = freq - 1;
                    x1 = 1ull << (shift + 31);

                    t1 = x1 / freq;
                    x0 += (x1 % freq) << 32;
                    t0 = x0 / freq;

                    e.rcp_freq = t0 + (t1 << 32);
                    ASSERT(shift >= 1);
                    e.rcp_shift = cb::check_value_cast<uint8>(shift - 1);

                    // With these values, 'q' is the correct quotient, so we
                    // have bias=start.
                    e.bias = cb::check_value_cast<uint16>(low);
                }

                low += freq;
            }
        ASSERT_RELEASE(low == (1U << M_bits));
}

static NOINLINE uint8* rans_encode(const uint8* buffer, ptrdiff_t buflen,
    uint8* comp_end, const rans_encode_table* entable)
{
    const uint64 min_x = RANS64_RENORM_LO;
    uint64 x = min_x;
    uint8* comp_ptr = comp_end;

    int cumprobtot_bits = entable->cumprobtot_bits;

    const rans_encode_table::entry* table = entable->table;

    for (ptrdiff_t i = buflen - 1; i >= 0; i--) {
        int sym = buffer[i];

        const rans_encode_table::entry& e = table[sym];

        uint64 freq = e.freq;

        uint64 x_max = freq << (RANS64_RENORM_HI_SHIFT - cumprobtot_bits);
        if (x >= x_max) {
            comp_ptr -= 4;
            *((uint32*)comp_ptr) = (uint32)x;
            x >>= 32;
        }

        // DURING_ASSERT( int low = entable->low[sym] );
        // DURING_ASSERT( uint64 x_check = ((x / freq) << cumprobtot_bits) + (x
        // % freq) + low );

        // x = C(s,x)
        x = x + e.bias
            + (__umulh(x, e.rcp_freq) >> e.rcp_shift) * e.M_minus_freq;

        // ASSERT( x == x_check );

        ASSERT(x >= min_x);
    }

    ASSERT(x >= min_x);

    comp_ptr -= 8;
    *((uint64*)comp_ptr) = x;

    return comp_ptr;
}

// rans_encode_int2
// encode two interleaved states

static NOINLINE uint8* rans_encode_int2(const uint8* buffer, ptrdiff_t buflen,
    uint8* comp_end, const rans_encode_table* entable)
{
    const uint64 min_x = RANS64_RENORM_LO;
    uint64 x1 = min_x;
    uint64 x2 = min_x;
    uint8* comp_ptr = comp_end;

    int cumprobtot_bits = entable->cumprobtot_bits;

    const rans_encode_table::entry* table = entable->table;

    for (ptrdiff_t i = buflen - 2; i >= 0; i -= 2) {
        int sym1 = buffer[i + 1];
        int sym2 = buffer[i];

        const rans_encode_table::entry& e1 = table[sym1];
        const rans_encode_table::entry& e2 = table[sym2];

        uint64 x1_max = ((uint64)e1.freq)
            << (RANS64_RENORM_HI_SHIFT - cumprobtot_bits);
        if (x1 >= x1_max) {
            comp_ptr -= 4;
            *((uint32*)comp_ptr) = (uint32)x1;
            x1 >>= 32;
        }

        uint64 x2_max = ((uint64)e2.freq)
            << (RANS64_RENORM_HI_SHIFT - cumprobtot_bits);
        if (x2 >= x2_max) {
            comp_ptr -= 4;
            *((uint32*)comp_ptr) = (uint32)x2;
            x2 >>= 32;
        }

        x1 = x1 + e1.bias
            + (__umulh(x1, e1.rcp_freq) >> e1.rcp_shift) * e1.M_minus_freq;

        ASSERT(x1 >= min_x);

        x2 = x2 + e2.bias
            + (__umulh(x2, e2.rcp_freq) >> e2.rcp_shift) * e2.M_minus_freq;

        ASSERT(x2 >= min_x);
    }

    if ((buflen & 1)) {
        int sym1 = buffer[0];

        const rans_encode_table::entry& e1 = table[sym1];

        uint64 x1_max = ((uint64)e1.freq)
            << (RANS64_RENORM_HI_SHIFT - cumprobtot_bits);
        if (x1 >= x1_max) {
            comp_ptr -= 4;
            *((uint32*)comp_ptr) = (uint32)x1;
            x1 >>= 32;
        }

        x1 = x1 + e1.bias
            + (__umulh(x1, e1.rcp_freq) >> e1.rcp_shift) * e1.M_minus_freq;
    }

    comp_ptr -= 8;
    *((uint64*)comp_ptr) = x1;

    comp_ptr -= 8;
    *((uint64*)comp_ptr) = x2;

    return comp_ptr;
}

struct rans_decode_table {
    // 8-byte decode entry
    // 350 mb/s
    struct entry {
        uint16 freq;
        uint16 xm_minus_low;
        uint8 sym;
        uint16 pad;
    };
    /*
    // 4-byte decode entry
    // 300 mb/s
    #pragma pack(push)
    #pragma pack(1)
    struct entry { uint16 freq; uint16 xm_minus_low; uint8 sym; };
    #pragma pack(pop)
    */

    // you could also do just a 1-byte cum2sym table
    // and then have freq & low as per-sym tables

    entry table[4096];
    int cumprobtot_bits;
};

static NOINLINE void rans_make_decode_table(rans_decode_table* detable,
    const uint32* normalized_counts, int M_bits, int alphabet)
{
    ASSERT(cb::sum(normalized_counts, normalized_counts + alphabet)
        == (1UL << M_bits));

    detable->cumprobtot_bits = M_bits;

    int low = 0;
        for
            LOOP(a, alphabet)
            {
                int count = normalized_counts[a];
                if (count == 0)
                    continue;

                /*
                // simple way :
                rans_decode_table::entry e;
                e.freq = cb::check_value_cast<uint16>( count );
                //e.low = cb::check_value_cast<uint16>( low );
                e.sym = cb::check_value_cast<uint8>( a );
                rans_decode_table::entry * table_ptr = detable->table + low;
                for LOOP(c,count)
                {
                        e.xm_minus_low = cb::check_value_cast<uint16>( c );
                        table_ptr[c] = e;
                }
                /*/
                // way faster :
                // fill first entry :
                rans_decode_table::entry* table_ptr = detable->table + low;
                rans_decode_table::entry& e = table_ptr[0];
                e.freq = cb::check_value_cast<uint16>(count);
                e.sym = cb::check_value_cast<uint8>(a);
                e.xm_minus_low = 0;
                if (count > 1) {
                    // copy first entry to all slots as a U64
                    COMPILER_ASSERT(sizeof(e) == sizeof(uint64));
                    uint64* table_u64 = (uint64*)table_ptr;
                    uint64 te = table_u64[0];
                    for (int c = 1; c < count; c++) {
                        table_u64[c] = te;
                        // fix xm_minus_low :
                        ((rans_decode_table::entry*)(table_u64 + c))
                            ->xm_minus_low
                            = cb::check_value_cast<uint16>(c);
                    }
                }
                /**/

                low += count;
            }

        ASSERT_RELEASE(low == (1 << M_bits));
}

static NOINLINE const uint8* rans_decode(const uint8* comp_start, uint8* buffer,
    ptrdiff_t buflen, const rans_decode_table* detable)
{
    const uint64 min_x = RANS64_RENORM_LO;
    const uint8* comp_ptr = comp_start;
    uint64 x = *((uint64*)comp_ptr);
    comp_ptr += 8;

    int cumprobtot_bits = detable->cumprobtot_bits;
    uint64 mask = (1 << cumprobtot_bits) - 1;

    const rans_decode_table::entry* table = detable->table;

    for (ptrdiff_t i = 0; i < buflen; i++) {
        ASSERT(x >= min_x);

        uint64 xm = x & mask;
        const rans_decode_table::entry& e = table[xm];

        // x = e.freq * (x >> cumprobtot_bits) + xm - e.low;
        x = e.freq * (x >> cumprobtot_bits) + e.xm_minus_low;

        buffer[i] = (uint8)e.sym;

        if (x < min_x) {
            x <<= 32;
            x |= *((uint32*)comp_ptr);
            comp_ptr += 4;
        }
    }

    return comp_ptr;
}

static NOINLINE const uint8* rans_decode_int2(const uint8* comp_start,
    uint8* buffer, ptrdiff_t buflen, const rans_decode_table* detable)
{
    const uint64 min_x = RANS64_RENORM_LO;
    const uint8* comp_ptr = comp_start;

    uint64 x2 = *((uint64*)comp_ptr);
    comp_ptr += 8;

    uint64 x1 = *((uint64*)comp_ptr);
    comp_ptr += 8;

    int cumprobtot_bits = detable->cumprobtot_bits;
    uint64 mask = (1 << cumprobtot_bits) - 1;

    const rans_decode_table::entry* table = detable->table;

    ptrdiff_t i = 0;

    if (buflen & 1) {
        ASSERT(x1 >= min_x);

        uint64 xm1 = x1 & mask;
        const rans_decode_table::entry& e1 = table[xm1];

        x1 = e1.freq * (x1 >> cumprobtot_bits) + e1.xm_minus_low;

        buffer[i] = (uint8)e1.sym;

        if (x1 < min_x) {
            x1 <<= 32;
            x1 |= *((uint32*)comp_ptr);
            comp_ptr += 4;
        }

        i++;
    }

    for (; i < buflen; i += 2) {
        ASSERT(x2 >= min_x);

        uint64 xm2 = x2 & mask;
        const rans_decode_table::entry& e2 = table[xm2];

        x2 = e2.freq * (x2 >> cumprobtot_bits) + e2.xm_minus_low;

        buffer[i] = (uint8)e2.sym;

        if (x2 < min_x) {
            x2 <<= 32;
            x2 |= *((uint32*)comp_ptr);
            comp_ptr += 4;
        }

        ASSERT(x1 >= min_x);

        uint64 xm1 = x1 & mask;
        const rans_decode_table::entry& e1 = table[xm1];

        x1 = e1.freq * (x1 >> cumprobtot_bits) + e1.xm_minus_low;

        buffer[i + 1] = (uint8)e1.sym;

        if (x1 < min_x) {
            x1 <<= 32;
            x1 |= *((uint32*)comp_ptr);
            comp_ptr += 4;
        }
    }

    return comp_ptr;
}

//================================================================

#ifdef _DEBUG
int timing_repeats = 1;
#else
int timing_repeats = 10;
#endif

// int L = 1024;
int L = 2048; // 11 bit
// int L = 4096; //

//===================================================================

int64 TestFile(const char* filename)
{
    int64 filelen;
    uint8* filebuf = (uint8*)cb::ReadWholeFile(filename, &filelen);
    if (filebuf == NULL) {
        lprintf("failed to load file\n");
        return 0;
    }

    uint8* decodebuf = (uint8*)CBALLOC(filelen);

    lprintfvar(filelen);

    uint32 unnormalized_counts[256] = { 0 };

        for
            LOOP(i, filelen) { unnormalized_counts[filebuf[i]]++; }

        int alphabet = 256;

        while (alphabet > 1 && unnormalized_counts[alphabet - 1] == 0)
            alphabet--;

        lprintfvar(alphabet);

        double H = entropy(unnormalized_counts, 256);
        lprintfvar(H);

        vector<uint32> normalized_counts;
        normalized_counts.resize(alphabet, (uint32)0);

        normalize_counts(
            normalized_counts.data(), unnormalized_counts, alphabet, L);

        ASSERT_RELEASE(
            cb::sum(normalized_counts.begin(), normalized_counts.end())
            == (uint32)L);

        int L_bits = ilog2ceil(L);
        ASSERT_RELEASE(L == (1 << L_bits));
        lprintfvar(L_bits);

        //---------------------------------------------------

        if (0) // log the tables
        {
            tans_encode_table entables;
            tans_make_encode_table(
                &entables, L, alphabet, normalized_counts.data());

            tans_decode_table detables;
            tans_make_decode_table(
                &detables, L, alphabet, normalized_counts.data());

            log_tables(&entables, &detables, normalized_counts.data());
        }

        //---------------------------------------------------

        vector<uint8> comp;
        comp.resize((int)filelen + 4096);
        uint8* comp_start = comp.data();
        uint8* comp_end = comp_start + comp.size() - 8;
        *comp_end = 0xCC;

        //---------------------------------------------------

        if (1) // rans
        {
            lprintf("rANS:\n");

            //---------------------------------------------

            uint8* rans_end;
            Timer::tsc_type encode_time_withtable = (uint64)-1;
            Timer::tsc_type encode_time = (uint64)-1;
            Timer::tsc_type decode_time_withtable = (uint64)-1;
            Timer::tsc_type decode_time = (uint64)-1;

            TIMING_LOOP_PRE(timing_repeats)
            {

                {
                    Timer::tsc_type t1 = Timer::rdtsc();

                    rans_encode_table entables;
                    rans_make_encode_table(
                        &entables, normalized_counts.data(), L_bits, alphabet);

                    Timer::tsc_type t2 = Timer::rdtsc();

                    rans_end = rans_encode_int2(
                        filebuf, (ptrdiff_t)filelen, comp_end, &entables);
                    // rans_end =
                    // rans_encode(filebuf,(ptrdiff_t)filelen,comp_end,&entables);

                    Timer::tsc_type t3 = Timer::rdtsc();
                    encode_time_withtable
                        = MIN(encode_time_withtable, (t3 - t1));
                    encode_time = MIN(encode_time, (t3 - t2));
                }

                TrashTheCache();

                {
                    Timer::tsc_type t1 = Timer::rdtsc();

                    rans_decode_table detables;
                    rans_make_decode_table(
                        &detables, normalized_counts.data(), L_bits, alphabet);

                    Timer::tsc_type t2 = Timer::rdtsc();

                    rans_decode_int2(
                        rans_end, decodebuf, (ptrdiff_t)filelen, &detables);
                    // rans_decode(rans_end,decodebuf,(ptrdiff_t)filelen,&detables);

                    Timer::tsc_type t3 = Timer::rdtsc();
                    decode_time_withtable
                        = MIN(decode_time_withtable, (t3 - t1));
                    decode_time = MIN(decode_time, (t3 - t2));
                }
            }
            TIMING_LOOP_POST();

            lprintf("rans wrote : %.6f bpb , %d bytes\n",
                (comp_end - rans_end) * 8.0 / filelen,
                (int)(comp_end - rans_end));

            lprintf("ticks to encode: %.2f decode: %.2f\n",
                encode_time / (double)filelen, decode_time / (double)filelen);

            lprintf("mbps encode: %.2f decode: %.2f\n", (double)filelen
                    / (1000 * 1000 * Timer::ConvertTicksToSeconds(encode_time)),
                (double)filelen / (1000 * 1000 * Timer::ConvertTicksToSeconds(
                                                     decode_time)));

            lprintf("withtable ticks to encode: %.2f decode: %.2f\n",
                encode_time_withtable / (double)filelen,
                decode_time_withtable / (double)filelen);

            lprintf("withtable mbps encode: %.2f decode: %.2f\n",
                (double)filelen / (1000 * 1000 * Timer::ConvertTicksToSeconds(
                                                     encode_time_withtable)),
                (double)filelen / (1000 * 1000 * Timer::ConvertTicksToSeconds(
                                                     decode_time_withtable)));

            lprintf("build table ticks encode: %a decode: %a\n",
                (encode_time_withtable - encode_time),
                (decode_time_withtable - decode_time));

            //===============================================
            // check it :

                for
                    LOOP(i, (int)filelen)
                    {
                        ASSERT_RELEASE(filebuf[i] == decodebuf[i]);
                    }

        } // rans

        //---------------------------------------------

        lprintf("tANS:\n");

        std::pair<void*, int> tans_end;
        Timer::tsc_type encode_time_withtable = (uint64)-1;
        Timer::tsc_type encode_time = (uint64)-1;
        Timer::tsc_type decode_time_withtable = (uint64)-1;
        Timer::tsc_type decode_time = (uint64)-1;

        TIMING_LOOP_PRE(timing_repeats)
        {

            {
                Timer::tsc_type t1 = Timer::rdtsc();

                tans_encode_table entables;
                tans_make_encode_table(
                    &entables, L, alphabet, normalized_counts.data());

                Timer::tsc_type t2 = Timer::rdtsc();

                // tans_end =
                // tans_encode(filebuf,(ptrdiff_t)filelen,comp_end,&entables);
                tans_end = tans_encode_int2(
                    filebuf, (ptrdiff_t)filelen, comp_end, &entables);

                Timer::tsc_type t3 = Timer::rdtsc();
                encode_time_withtable = MIN(encode_time_withtable, (t3 - t1));
                encode_time = MIN(encode_time, (t3 - t2));
            }

            TrashTheCache();

            {
                Timer::tsc_type t1 = Timer::rdtsc();

                tans_decode_table detables;
                tans_make_decode_table(
                    &detables, L, alphabet, normalized_counts.data());

                Timer::tsc_type t2 = Timer::rdtsc();

                // tans_decode(tans_end.first,tans_end.second,decodebuf,(ptrdiff_t)filelen,&detables);
                tans_decode_int2(tans_end.first, tans_end.second, decodebuf,
                    (ptrdiff_t)filelen, &detables);

                Timer::tsc_type t3 = Timer::rdtsc();
                decode_time_withtable = MIN(decode_time_withtable, (t3 - t1));
                decode_time = MIN(decode_time, (t3 - t2));
            }
        }
        TIMING_LOOP_POST();

        ASSERT(*comp_end == 0xCC);
        // we wrote from bout_final_ptr[0] to comp_end[-1]
        //*comp_end = 0xEE;
        // bout_final_ptr[-1] = 0xEE;

        int64 wrote_bits
            = (int64)8 * (comp_end - (uint8*)tans_end.first) + tans_end.second;
        lprintf("wrote %a bpb\n", (float)wrote_bits / filelen);

        lprintf("ticks to encode: %.2f decode: %.2f\n",
            encode_time / (double)filelen, decode_time / (double)filelen);

        lprintf("mbps encode: %.2f decode: %.2f\n", (double)filelen
                / (1000 * 1000 * Timer::ConvertTicksToSeconds(encode_time)),
            (double)filelen
                / (1000 * 1000 * Timer::ConvertTicksToSeconds(decode_time)));

        lprintf("withtable ticks to encode: %.2f decode: %.2f\n",
            encode_time_withtable / (double)filelen,
            decode_time_withtable / (double)filelen);

        lprintf("withtable mbps encode: %.2f decode: %.2f\n",
            (double)filelen / (1000 * 1000 * Timer::ConvertTicksToSeconds(
                                                 encode_time_withtable)),
            (double)filelen / (1000 * 1000 * Timer::ConvertTicksToSeconds(
                                                 decode_time_withtable)));

        lprintf("build table ticks encode: %a decode: %a\n",
            (encode_time_withtable - encode_time),
            (decode_time_withtable - decode_time));

        //===============================================
        // check it :

        for
            LOOP(i, (int)filelen)
            {
                ASSERT_RELEASE(filebuf[i] == decodebuf[i]);
            }

        //===============================================

        CBFREE(filebuf);
        CBFREE(decodebuf);

        return wrote_bits;
}

int main(int argc, char* argv[])
{
    cb::vector<cb::String> files;

    for (int argi = 1; argi < argc; argi++) {
        if (NameIsDir(argv[argi])) {
            EnumFiles(argv[argi], false, &files);
        } else {
            files.push_back(cb::String(argv[argi]));
        }
    }

    int64 total = 0;

        for
            LOOPVEC(i, files)
            {
                const char* fname = files[i].CStr();

                lprintf("loading : %s\n", fname);

                total += TestFile(fname);
            }

        // lprintfvar(total);
        lprintf("total bytes out : %.2f\n", total / 8.0);

        return 0;
}

//=====================================================================
