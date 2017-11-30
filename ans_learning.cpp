/*************

ans_learning.cpp

by cbloom

v2 02-21-2014

---------------------------

incremental educational implementations of tANS

with logging so you can trace what's happening

---------------------------

please link to this :

http://cbloomrants.blogspot.com/2014/02/02-18-14-understanding-ans-conclusion.html

*******************/

#include "cblib/FileEnum.h"
#include "cblib/FileUtil.h"
#include "cblib/Log.h"
#include "cblib/LogUtil.h"
#include "cblib/Rand.h"
#include "cblib/inc.h"
#include "cblib/vector.h"

USE_CB

//===========================================================================================

struct sort_sym {
    // sort array of "sym" by comparing "rank"
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

/*
normalize_counts :
make {to} sum to "to_sum_desired"

normalized_counts = Fs
to_sum_desired = M
(here we will use M=L)
*/

static void normalize_counts(
    uint32* to, const uint32* from, int alphabet, int to_sum_desired)
{
    ASSERT_RELEASE(to_sum_desired > alphabet);

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

        // now apply correction to make output sum to "to_sum"
        int32 correction = to_sum_desired - to_sum;
        if (correction != 0) {
            lprintfvar(correction);
            int32 correction_sign = (correction > 0) ? 1 : -1;

            vector<sort_sym> heap;
            heap.reserve(alphabet);

                for
                    LOOP(i, alphabet)
                    {
                        if (from[i] == 0)
                            continue;
                        ASSERT(to[i] != 0);
                        if (to[i] > 1 || correction_sign == 1) {
                            float change
                                = (float)log(
                                      1.0 + correction_sign * (1.0 / to[i]))
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
                        float change
                            = (float)log(1.0 + correction_sign * (1.0 / to[i]))
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

// total code len from given frequencies
static double codelen_frequencies(const uint32* true_counts,
    uint32 true_count_sum, const uint32* normalized_counts,
    uint32 normalized_sum, int alphabet)
{
    ASSERT(cb::sum(normalized_counts, normalized_counts + alphabet)
        == normalized_sum);
    ASSERT(cb::sum(true_counts, true_counts + alphabet) == true_count_sum);

    double total_cl = 0;
        for
            LOOP(a, alphabet)
            {
                if (true_counts[a] == 0)
                    continue;
                ASSERT_RELEASE(normalized_counts[a] != 0);
                double cl = log2((double)normalized_sum / normalized_counts[a]);
                total_cl += cl * true_counts[a] / true_count_sum;
            }
        return total_cl;
}

//===========================================================================================
/*

make_sort

takes "normalized_counts"

produces an "output string"

this is the key step in tANS table construction

*/

void make_sort_truesort(int* sorted_syms, int sort_count,
    const uint32* normalized_counts, int alphabet)
{
    ASSERT((int)cb::sum(normalized_counts, normalized_counts + alphabet)
        == sort_count);

    vector<sort_sym> sort_syms;
    sort_syms.resize(sort_count);

    int s = 0;

        for
            LOOP(sym, alphabet)
            {
                uint32 count = normalized_counts[sym];
                if (count == 0)
                    continue;

                float invp = 1.f / count;

                // float base = 0.5f * invp; // Duda Precise
                // best if normalized_counts[] are true counts

                float base = invp; // 1.0 - better in practice

                for
                    LOOP(c, (int)count)
                    {
                        sort_syms[s].sym = sym;
                        sort_syms[s].rank = base + c * invp;
                        s++;
                    }
            }

        ASSERT_RELEASE(s == sort_count);

        // stable sort so rank ties are broken consistently
        std::stable_sort(sort_syms.begin(), sort_syms.end());

        // sort for lowest rank first :
        ASSERT_RELEASE(sort_syms[0].rank <= sort_syms[1].rank);

        for
            LOOP(s, sort_count) { sorted_syms[s] = sort_syms[s].sym; }
}

// sort by striding and looking for holes :
void make_sort_greedy(int* sorted_syms, int sort_count,
    const uint32* normalized_counts, int alphabet)
{
    ASSERT((int)cb::sum(normalized_counts, normalized_counts + alphabet)
        == sort_count);

    vector<sort_sym> sort_syms;
    sort_syms.resize(alphabet);

        for
            LOOP(s, sort_count) { sorted_syms[s] = -1; }

        // sort for least probable first :
        // ABCABAACBABACBAB
        for
            LOOP(a, alphabet)
            {
                sort_syms[a].sym = a;
                sort_syms[a].rank = (float)normalized_counts[a];
            }

        // sort for lps : 1803793.88 | 144683407.38
        // no sort : 1801351.13	 | 144685135.75
        // sort for mps : 1801187.75 | 144693738.13

        // sort for lowest rank first :
        // std::sort(sort_syms.begin(),sort_syms.end());
        // ASSERT_RELEASE( sort_syms[0].rank <= sort_syms[1].rank );

        for
            LOOP(a, alphabet)
            {
                int sym = sort_syms[a].sym;

                uint32 count = normalized_counts[sym];
                if (count == 0)
                    continue;

                uint32 step = (sort_count + (count / 2)) / count;

                uint32 first = step / 2;

                for
                    LOOP(c, (int)count)
                    {
                        uint32 base = first + step * c;
                        // int base = froundint(( c + 0.5 ) * sort_count /
                        // count);
                        // base = Clamp(base,0,sort_count-1);

                        if (sorted_syms[base] == -1) {
                            sorted_syms[base] = sym;
                        } else {
                            // search out from base :
                            for (int i = 0;; i++) {
                                // 0123 -> 1,-1,2,-2
                                int sign = (i & 1);
                                int mag = (i >> 1) + 1;
                                int offset = sign ? -mag : mag; // 1801351
                                // int offset = sign ? mag : -mag; // 1801423 ,
                                // slightly worse / doesn't matter
                                // offset = i; // 1801608 , not much worse

                                // modulo ? slightly better than reject
                                uint32 pos
                                    = (base + offset + sort_count) % sort_count;

                                // reject out of range :
                                // int pos = base + offset;
                                // if ( pos < 0 || pos >= sort_count ) continue;

                                if (sorted_syms[pos] == -1) {
                                    sorted_syms[pos] = sym;
                                    break;
                                }
                            }
                        }
                    }
            }
}

// simpler version of greedy insertion :
void make_sort_greedy_simple(int* sorted_syms, int sort_count,
    const uint32* normalized_counts, int alphabet)
{
    ASSERT((int)cb::sum(normalized_counts, normalized_counts + alphabet)
        == sort_count);

    // make all slots empty :
        for
            LOOP(s, sort_count) { sorted_syms[s] = -1; }

        for
            LOOP(a, alphabet)
            {
                uint32 count = normalized_counts[a];
                if (count == 0)
                    continue;

                uint32 step = (sort_count + (count / 2)) / count;
                uint32 first = step / 2;

                for
                    LOOP(c, (int)count)
                    {
                        uint32 slot = first + step * c;
                        // find an empty spot near slot :
                        for (;;) {
                            if (sorted_syms[slot] == -1) {
                                sorted_syms[slot] = a;
                                break;
                            }
                            slot = (slot + 1) % sort_count;
                        }
                    }
            }
}

int bitreverse(int val, int bits)
{
    int ret = 0;
    while (bits--) {
        ret <<= 1;
        ret |= (val & 1);
        val >>= 1;
    }
    return ret;
}

void make_sort_bitreverse(int* sorted_syms, int sort_count,
    const uint32* normalized_counts, int alphabet)
{
    ASSERT((int)cb::sum(normalized_counts, normalized_counts + alphabet)
        == sort_count);

    // only works for pow2 counts :
    int sort_count_bits = intlog2(sort_count);
    ASSERT_RELEASE(sort_count == 1 << sort_count_bits);

    int s = 0;
        for
            LOOP(a, alphabet)
            {
                int count = normalized_counts[a];

                for
                    LOOP(c, count)
                    {
                        // sorted_syms[s] = a; // 1824053.75

                        // 1805230.75
                        int sr = bitreverse(s, sort_count_bits);
                        sorted_syms[sr] = a;

                        s++;
                    }
            }
        ASSERT(s == sort_count);
}

// Yann's heuristic from FSE.c :
void make_sort_stepping_yann(int* sorted_syms, int sort_count,
    const uint32* normalized_counts, int alphabet)
{
    vector<int> symbolPos;
    symbolPos.resize(alphabet + 1);

    // symbol start positions
    symbolPos[0] = 0;
    for (int i = 1; i < alphabet; i++) {
        symbolPos[i] = symbolPos[i - 1] + normalized_counts[i - 1];
    }
    symbolPos[alphabet] = sort_count + 1;

    const int step = (sort_count >> 1) + (sort_count >> 3) + 1;

    memset(sorted_syms, 0xFF, sort_count * sizeof(int)); // debug check

    // Spread symbols
    int position = 0;
    int s = 0;
    while (!symbolPos[s + 1])
        s++;
    for (int i = 0; i < sort_count; i++) {
        sorted_syms[position] = s;
        while (i + 2 > symbolPos[s + 1])
            s++;
        position = (position + step) % sort_count;
    }

    // debug
    for
        LOOP(s, sort_count) { ASSERT_RELEASE(sorted_syms[s] != -1); }
}

// simpler version of Yann's way :
void make_sort_stepping(int* sorted_syms, int sort_count,
    const uint32* normalized_counts, int alphabet)
{
    ASSERT((int)cb::sum(normalized_counts, normalized_counts + alphabet)
        == sort_count);

    const int step = (sort_count >> 1) + (sort_count >> 3) + 1;

    int s = 0;
        for
            LOOP(a, alphabet)
            {
                int count = normalized_counts[a];

                for
                    LOOP(c, count)
                    {
                        sorted_syms[s] = a;

                        s = (s + step) % sort_count;
                    }
            }
}

#define make_sort make_sort_truesort // best

//===========================================================================================
/*

tans_tables1

most straightforward naive tans implementation

encode_table[] is a huge 2d array of state & symbol

*/

struct tans_tables1 {
    int L; // renorm is to [L,2L-1]
    int state_count, alphabet;

    // encode_table = C(s,x)
    vector<int> encode_table; // index by state & symbol, gives next state , use
    // [ x * alphabet + s ]
    vector<int> max_state; // index by symbol, gives max state that can encode
    // that symbol

    // decode_table = D(x)
    struct decode_entry {
        int next_state;
        int sym;
    };
    vector<decode_entry> decode_table; // index by state
};

void make_tables1(tans_tables1* tables, int L, const int* sorted_syms,
    int alphabet, const uint32* normalized_counts)
{
    ASSERT(
        cb::sum(normalized_counts, normalized_counts + alphabet) == (uint32)L);

    tables->L = L;
    int state_count = 2 * L;
    tables->state_count = state_count;
    tables->alphabet = alphabet;

    tables->encode_table.resize(state_count * alphabet, 0);
    tables->max_state.resize(alphabet, 0);

    // encoder precursor range Is is [Fs,2Fs-1]
    // set up max_state[] to be the first state we want to fill
        for
            LOOP(a, alphabet)
            {
                // -1 because we do a pre ++
                // first state for each symbol is Fs
                tables->max_state[a] = normalized_counts[a] - 1;
            }

        tables->decode_table.resize(state_count);

        // destination states are in [L,2L-1]
        for (int to_state = L; to_state < 2 * L; to_state++) {
            // to_state = x'

            // walking over the output string tells us our output symbol :
            // sym = s
            int sym = sorted_syms[to_state - L];

            // from_state is the next one for this symbol in Is
            // from_state = x
            tables->max_state[sym]++;
            int from_state = tables->max_state[sym];

            // fill C(s,x) = x'
            tables->encode_table[from_state * alphabet + sym] = to_state;

            // and D(x') = {s,x}
            ASSERT(tables->decode_table[to_state].next_state == 0);
            tables->decode_table[to_state].next_state = from_state;
            tables->decode_table[to_state].sym = sym;
        }

        // check streamability :
        for
            LOOP(a, alphabet)
            {
                if (normalized_counts[a] == 0)
                    continue;
                int m = tables->max_state[a];
                ASSERT_RELEASE(m == (int)normalized_counts[a] * 2 - 1);
                // m is the max state that has an encoder entry
                // step to the next state that doesn't (m+1)
                // then stream it down (/2) :
                int fm_state = (m + 1) >> 1;
                int to_state = tables->encode_table[fm_state * alphabet + a];
                ASSERT_RELEASE(to_state >= L && to_state < 2 * L);
            }
}

/*

tans_tables2

starting to optimize

fix the bloated table sizes

tans_tables1->decode_table is empty in [0,L-1],
so just leave that off

encode_table is packed into L entries
by finding the consecutive parts of the big 2d tans_tables1->encode_table

*/

struct tans_tables2 {
    int L; // renorm is to [L,2L-1]
    int state_count, alphabet;

    vector<int> encode_table_packed; // encode_table_packed[ x +
    // packed_table_offset[s] ]
    vector<int> packed_table_offset; // index by symbol
    vector<int> max_state; // index by symbol, gives max state that can encode
    // that symbol

    struct decode_entry {
        int next_state;
        int sym;
    };
    vector<decode_entry> decode_table_packed; // index by [state - L]
};

void make_tables2(tans_tables2* tables, int L, const int* sorted_syms,
    int alphabet, const uint32* normalized_counts)
{
    ASSERT(
        cb::sum(normalized_counts, normalized_counts + alphabet) == (uint32)L);

    tables->L = L;
    tables->state_count = 2 * L;
    tables->alphabet = alphabet;

    // only L entries per table here :
    tables->encode_table_packed.resize(L);
    tables->decode_table_packed.resize(L);

    tables->max_state.resize(alphabet, 0);
    tables->packed_table_offset.resize(alphabet, 0);

    // set up max_state[] to be the first state we want to fill
    uint32 cumulative_count = 0;
        for
            LOOP(a, alphabet)
            {
                // -1 because we do a pre ++
                uint32 c = normalized_counts[a];
                tables->max_state[a] = c - 1;

                // packed_table_offset for encode table
                // cumulative_count is the number of entries used so far
                // we will use c entries
                // subtract c from offset because from_state is in [c,2c-1]
                tables->packed_table_offset[a] = cumulative_count - c;
                cumulative_count += c;
            }
        ASSERT_RELEASE(cumulative_count == (uint32)L);

        // destination states are in [L,2L-1]
        for (int i = 0; i < L; i++) {
            int to_state = L + i;

            int sym = sorted_syms[i];
            tables->max_state[sym]++;
            int from_state = tables->max_state[sym];

            // packed_table_offset[sym] to find the linear portion of
            // encode_table_packed for sym
            tables->encode_table_packed[from_state
                + tables->packed_table_offset[sym]]
                = to_state;

            ASSERT(tables->decode_table_packed[i].next_state == 0);
            tables->decode_table_packed[i].next_state = from_state;
            tables->decode_table_packed[i].sym = sym;
        }
}

/**

tans_tables3

optimizing more

precompute the number of bits to renormalize at each state transition

**/

struct tans_tables3 {
    int L; // renorm is to [L,2L-1]
    int state_count, alphabet;

    struct encode_sym_data {
        int min_bits;
        int more_bits_threshold;
        int max_state;
        int table_offset;
    };

    vector<encode_sym_data> encode_sym_table; // index by symbol
    vector<int> encode_table_packed; // encode_table_packed[ x +
    // packed_table_offset[s] ]

    struct decode_entry {
        int next_state;
        int sym;
        int num_bits;
    };
    vector<decode_entry> decode_table_packed; // index by [state - L]
};

void make_tables3(tans_tables3* tables, int L, const int* sorted_syms,
    int alphabet, const uint32* normalized_counts)
{
    ASSERT(
        cb::sum(normalized_counts, normalized_counts + alphabet) == (uint32)L);

    tables->L = L;
    tables->state_count = 2 * L;
    tables->alphabet = alphabet;

    tables->encode_table_packed.resize(L);
    tables->decode_table_packed.resize(L);
    tables->encode_sym_table.resize(alphabet);

    // set up max_state[] to be the first state we want to fill
    uint32 cumulative_count = 0;
        for
            LOOP(a, alphabet)
            {
                int c = normalized_counts[a];
                if (c == 0)
                    continue;

                // -1 because we do a pre ++ below
                tables->encode_sym_table[a].max_state = c - 1;
                tables->encode_sym_table[a].table_offset = cumulative_count - c;

                //-----------------------------------------------------
                // encoder num_bits
                // we start at any state in [L,2L-1]
                // we need to reach [c,2c-1]
                // we can do it with "min_bits" bits output
                // + one more if state >= threshold

                int max_fm_state = 2 * c - 1;
                int state = L;
                // can use BSR/intlog2 obviously (see ans_fast)
                int num_bits = 0;
                while (state > max_fm_state) {
                    state >>= 1;
                    num_bits++;
                }
                tables->encode_sym_table[a].min_bits = num_bits;

                ASSERT((L >> num_bits) <= max_fm_state && (L >> num_bits) >= c);

                tables->encode_sym_table[a].more_bits_threshold
                    = (max_fm_state + 1) << num_bits;

                // obviously any state >= more_bits_threshold
                // when >>num_bits
                // will still be > max_fm_state
                // but one more bit is enough to get any state up to max in
                // range :
                ASSERT((2 * L - 1) >> (num_bits + 1) <= max_fm_state);

                cumulative_count += c;
            }
        ASSERT_RELEASE(cumulative_count == (uint32)L);

        // destination states are in [L,2L-1]
        for (int i = 0; i < L; i++) {
            int to_state = L + i;

            int sym = sorted_syms[i];

            tables->encode_sym_table[sym].max_state++;
            int from_state = tables->encode_sym_table[sym].max_state;

            tables->encode_table_packed[from_state
                + tables->encode_sym_table[sym].table_offset]
                = to_state;

            //-----------------------------------------
            // decoder numbits :
            // we take the state transition from "to_state" back to "from_state"
            // once we are in from_state we then need to read bits to get back
            // to >= L

            int decoder_state = from_state;
            int num_bits = 0;
            // can use BSR/intlog2 obviously (see ans_fast)
            while (decoder_state < L) {
                decoder_state <<= 1;
                num_bits++;
            }

            ASSERT(tables->decode_table_packed[i].next_state == 0);
            tables->decode_table_packed[i].next_state = from_state;
            tables->decode_table_packed[i].num_bits = num_bits;
            tables->decode_table_packed[i].sym = sym;
        }
}

//===========================================================================

void log_tables(const tans_tables1* tables)
{
    int max_encode_state = 0;
    int num_syms = 0;

    lprintf(" S ");
        for
            LOOP(a, tables->alphabet)
            {
                int m = tables->max_state[a];
                if (m == -1)
                    continue;
                max_encode_state = MAX(max_encode_state, m);
                lprintf("|%3d", a);
                num_syms++;
            }
        lprintf("\n");

        lprintf("---");
        for
            LOOP(a, num_syms) { lprintf("|---"); }
        lprintf("\n");

        for (int s = 1; s <= max_encode_state; s++) {
            bool any = false;
                for
                    LOOP(a, tables->alphabet)
                    {
                        int to_state
                            = tables->encode_table[s * tables->alphabet + a];
                        if (to_state != 0) {
                            any = true;
                            break;
                        }
                    }
                if (!any)
                    continue;

                lprintf("%3d", s);

                for
                    LOOP(a, tables->alphabet)
                    {
                        if (tables->max_state[a] == -1)
                            continue;

                        int to_state
                            = tables->encode_table[s * tables->alphabet + a];
                        if (to_state == 0)
                            lprintf("|   ");
                        else
                            lprintf("|%3d", to_state);
                    }
                lprintf("\n");
        }
}

//=================================================================================
/*

very simple bit input & output

note that bit in & out act in opposite orders
eg. like a stack pushing & popping

*/

struct bit_output {
    uint32 bits;
    int numbits;
    uint8* ptr;

    void start(uint8* to)
    {
        bits = 0;
        numbits = 0;
        ptr = to;
    }

    void put(uint32 val, int nb)
    {
        bits <<= nb;
        bits |= val;
        ASSERT((numbits + nb) <= 32);
        numbits += nb;
        while (numbits >= 8) {
            *ptr++ = (uint8)(bits >> (numbits - 8));
            numbits -= 8;
        }
        ASSERT(numbits >= 0 && numbits < 8);
    }

    int flush()
    {
        ASSERT(numbits >= 0 && numbits < 8);
        *ptr = (uint8)bits;
        return numbits;
    }
};

// bit_input reads bits in *reverse* from bit_output
struct bit_input {
    uint32 bits;
    int numbits;
    const uint8* ptr;

    void start(const uint8* last, int last_bits)
    {
        bits = *last;
        numbits = last_bits;
        ptr = last - 1;
    }

    uint32 get(int nb)
    {
        while (numbits < nb) {
            bits |= (*ptr--) << numbits;
            numbits += 8;
        }

        uint32 got = bits & ((1 << nb) - 1);

        bits >>= nb;
        numbits -= nb;

        return got;
    }
};

//=========================================================================
/*

tans_encode & tans_decode

from the various types of table

*/

void tans_encode(int& ref_state, bit_output& output, int sym,
    const tans_tables1* table, bool log)
{
    int state = ref_state;

    ASSERT(sym >= 0 && sym < table->alphabet);
    // x starts in I = [L,2L-1]
    ASSERT(state >= table->L && state < table->state_count);

    // encode renormalization
    // bring state down to the precursor range Is=[Fs,2Fs-1]
    int ms = table->max_state[sym];
    while (state > ms) {
        // output bits from state to reduce it :
        uint32 outbit = (state & 1);
        state >>= 1;
        output.put(outbit, 1);
        if (log)
            lprintf("%d+%d->", state, outbit);
    }

    // x' = C(s,x)
    int to_state = table->encode_table[state * table->alphabet + sym];

    // end state of encode must be in renormalization range :
    ASSERT(to_state >= table->L && to_state < table->state_count);

    ref_state = to_state;
}

int tans_decode(
    int& ref_state, bit_input& input, const tans_tables1* table, bool log)
{
    int state = ref_state;

    // first input bits to get state up to renormalization range :
    while (state < table->L) {
        uint32 inbit = input.get(1);
        state <<= 1;
        state |= inbit;
        if (log)
            lprintf("%d+%d->", inbit, state);
    }

    ASSERT(state >= table->L && state < table->state_count);

    // look up D(x)
    const tans_tables1::decode_entry& de = table->decode_table[state];

    int sym = de.sym;
    int to_state = de.next_state;

    ASSERT(to_state != 0 && to_state < table->state_count);

    ref_state = to_state;

    ASSERT(sym >= 0 && sym < table->alphabet);
    return sym;
}

void tans_encode(
    int& ref_state, bit_output& output, int sym, const tans_tables2* table)
{
    int state = ref_state;

    ASSERT(sym >= 0 && sym < table->alphabet);
    ASSERT(state >= 1 && state < table->state_count);

    int ms = table->max_state[sym];
    while (state > ms) {
        // output bits from state to reduce it :
        uint32 outbit = (state & 1);
        state >>= 1;
        output.put(outbit, 1);
    }

    int to_state
        = table->encode_table_packed[state + table->packed_table_offset[sym]];

    // end state of encode must be in renormalization range :
    ASSERT(to_state >= table->L && to_state < 2 * table->L);

    ref_state = to_state;
}

int tans_decode(int& ref_state, bit_input& input, const tans_tables2* table)
{
    int state = ref_state;

    ASSERT(state >= 1 && state < table->state_count);

    while (state < table->L) {
        uint32 inbit = input.get(1);
        state <<= 1;
        state |= inbit;
    }

    ASSERT(state >= (table->L) && state < 2 * table->L);

    const tans_tables2::decode_entry& de
        = table->decode_table_packed[state - table->L];

    int sym = de.sym;
    int to_state = de.next_state;

    ASSERT(to_state != 0 && to_state < table->state_count);

    ref_state = to_state;

    ASSERT(sym >= 0 && sym < table->alphabet);
    return sym;
}

void tans_encode(
    int& ref_state, bit_output& output, int sym, const tans_tables3* table)
{
    int state = ref_state;

    ASSERT(sym >= 0 && sym < table->alphabet);
    ASSERT(state >= 1 && state < table->state_count);

    // per-symbol encode info :
    const tans_tables3::encode_sym_data& esd = table->encode_sym_table[sym];

    // find number of bits to get into Is range :
    int nb = esd.min_bits;
    if (state >= esd.more_bits_threshold)
        nb++;

    output.put(state & ((1 << nb) - 1), nb);
    state >>= nb;

    ASSERT(state <= esd.max_state && state >= (esd.max_state / 2));

    int to_state = table->encode_table_packed[state + esd.table_offset];

    // end state of encode must be in renormalization range :
    ASSERT(to_state >= table->L && to_state < 2 * table->L);

    ref_state = to_state;
}

int tans_decode(int& ref_state, bit_input& input, const tans_tables3* table)
{
    int state = ref_state;

    // different convention here
    // state is renormalized already coming in :
    ASSERT(state >= (table->L) && state < 2 * table->L);

    // look up D(x) :
    const tans_tables3::decode_entry& de
        = table->decode_table_packed[state - table->L];

    int sym = de.sym;
    state = de.next_state;

    // renormalize from the transition on exit
    int nb = de.num_bits;
    state <<= nb;
    state |= input.get(nb);

    ASSERT(state >= table->L && state < table->state_count);

    ref_state = state;

    return sym;
}

//===================================================================
/*

test1

synthetic small alphabet test

lots of logging

*/

void test1()
{
    /*
    int L = 16;
    int alphabet = 3;
    //uint32 unnormalized_counts[3] = { 3,3,2 };
    uint32 unnormalized_counts[3] = { 7,6,3 };
    */
    int L = 16;
    int alphabet = 4;
    uint32 unnormalized_counts[4] = { 8, 5, 2, 1 };

    ASSERT_RELEASE(alphabet == ARRAY_SIZE(unnormalized_counts));

    vector<uint32> normalized_counts(alphabet, (uint32)0);
    normalize_counts(
        normalized_counts.data(), unnormalized_counts, alphabet, L);

    lprintfCIntArray(
        (int*)normalized_counts.data(), alphabet, "normalized_counts", 16, 4);

    vector<int> sorted_syms(L, 0);
    make_sort(sorted_syms.data(), L, normalized_counts.data(), alphabet);

    if (1) // log the sort
    {
        lprintf("sort:\n");
                for
                    LOOP(i, L)
                    {
                        int c = sorted_syms[i] + 'A';
                        lprintf("%c", c);
                    }
                lprintf("\n");
    }

    tans_tables1 tables;
    make_tables1(
        &tables, L, sorted_syms.data(), alphabet, normalized_counts.data());

    log_tables(&tables);

    int state_bits = intlog2ceil((float)tables.state_count);

    //-------------------------------------

    int test[10];
    int test_count = ARRAY_SIZE(test);

    vector<uint8> comp;
    comp.resize(test_count * 2 + 256);
    uint8* comp_start = comp.data();

    bit_output bout;
    bout.start(comp_start);

    int bout_final_bits;
    uint8* bout_final_ptr;

    {
        int state = L;

        for
            LOOPBACK(i, test_count)
            {
                // test[i] = irandmod(alphabet);
                // test[i] =
                // random_symbol_weighted(normalized_counts.data(),L,alphabet);
                test[i] = sorted_syms[irandmod(L)];
                lprintf("encode %d, state %d ->", test[i], state);
                tans_encode(state, bout, test[i], &tables, true);
                lprintf("%d\n", state);
            }

        bout.put(state, state_bits);

        bout_final_bits = bout.flush();
        bout_final_ptr = bout.ptr;
    }

    {
        bit_input bin;
        bin.start(bout_final_ptr, bout_final_bits);

        int state = bin.get(state_bits);

        for
            LOOP(i, test_count)
            {
                lprintf("decode state %d ->", state);
                int got = tans_decode(state, bin, &tables, true);
                lprintf("%d = %d\n", state, got);
                ASSERT_RELEASE(got == test[i]);
            }

        while (state < L) {
            state <<= 1;
            state |= bin.get(1);
        }

        ASSERT(state == L);
    }

    ptrdiff_t wrote_bits
        = 8 * (bout_final_ptr - comp_start + 1) + bout_final_bits;
    lprintf("wrote %a bits, = %a bps\n", wrote_bits,
        (float)wrote_bits / test_count);
}

//===========================================
/*

test on a file :

*/

int64 TestFile(const char* filename)
{
    int64 filelen;
    uint8* filebuf = (uint8*)cb::ReadWholeFile(filename, &filelen);
    if (filebuf == NULL) {
        lprintf("failed to load file\n");
        return 0;
    }

    lprintfvar(filelen);

    int L = 1024;
    // int L = 4096;
    // int L = 16;
    // int L = 64;

    int alphabet = 256;

    uint32 unnormalized_counts[256] = { 0 };

        for
            LOOP(i, filelen) { unnormalized_counts[filebuf[i]]++; }

        if (1) {
            while (alphabet > 1 && unnormalized_counts[alphabet - 1] == 0)
                alphabet--;

            lprintfvar(alphabet);
        }

        double H = entropy(unnormalized_counts, 256);
        lprintfvar(H);

        vector<uint32> normalized_counts;
        normalized_counts.resize(alphabet, (uint32)0);

        normalize_counts(
            normalized_counts.data(), unnormalized_counts, alphabet, L);

        ASSERT_RELEASE(
            cb::sum(normalized_counts.begin(), normalized_counts.end())
            == (uint32)L);

        double CL = codelen_frequencies(unnormalized_counts, (uint32)filelen,
            normalized_counts.data(), L, alphabet);
        lprintfvar(CL);

        vector<int> sorted_syms(L, 0);
        make_sort(sorted_syms.data(), L, normalized_counts.data(), alphabet);

        //---------------------------------------------------

        // tans_tables2 tables;
        // make_tables2(&tables,L,sorted_syms.data(),alphabet,normalized_counts.data());

        tans_tables3 tables;
        make_tables3(
            &tables, L, sorted_syms.data(), alphabet, normalized_counts.data());

        int state_bits = intlog2ceil((float)tables.state_count);

        //---------------------------------------------------

        vector<uint8> comp;
        comp.resize((int)filelen + 4096);
        uint8* comp_start = comp.data();

        bit_output bout;
        bout.start(comp_start);

        int bout_final_bits;
        uint8* bout_final_ptr;

        {
            int state = L;

        for
            LOOPBACK(i, (int)filelen)
            {
                ASSERT_RELEASE(normalized_counts[filebuf[i]] != 0);
                tans_encode(state, bout, filebuf[i], &tables);
            }

        bout.put(state, state_bits);

        bout_final_bits = bout.flush();
        bout_final_ptr = bout.ptr;
        }

        int64 wrote_bits
            = (int64)8 * (bout_final_ptr - comp_start + 1) + bout_final_bits;
        lprintf("wrote %a bpb\n", (float)wrote_bits / filelen);

        {
            bit_input bin;
            bin.start(bout_final_ptr, bout_final_bits);

            int state = bin.get(state_bits);

        for
            LOOP(i, (int)filelen)
            {
                int got = tans_decode(state, bin, &tables);
                ASSERT_RELEASE(got == filebuf[i]);
            }

        while (state < L) {
            state <<= 1;
            state |= bin.get(1);
        }

        ASSERT(state == L);
        }

        //===============================================

        CBFREE(filebuf);

        return wrote_bits;
}

int main(int argc, char* argv[])
{
    if (argc >= 2) // file or dir arg
    {
        cb::vector<cb::String> files;

        if (cb::NameIsDir(argv[1])) {
            cb::EnumFiles(argv[1], false, &files);
        } else {
            files.push_back(cb::String(argv[1]));
        }

        int64 total = 0;

                for
                    LOOPVEC(i, files)
                    {
                        const char* fname = files[i].CStr();

                        lprintf("loading : %s\n", fname);

                        total += TestFile(fname);
                    }

                lprintf("total bytes out : %.2f\n", total / 8.0);
    } else // no arg
    {
        test1();
        // test_stats();
    }

    return 0;
}