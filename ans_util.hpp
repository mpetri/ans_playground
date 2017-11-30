#pragma once

#include <algorithm>
#include <assert.h>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <numeric>

namespace consts {
const uint64_t SIGMA = 256;
const uint64_t TABLE_SHIFT = 12;
const uint64_t M = 1ULL << TABLE_SHIFT;
const uint8_t PADDING = 8;
}

/* do this so the loop can be nicely vectorized */
void histogram8(const uint8_t* in8, size_t n, uint64_t* freqs)
{
    uint64_t F1[consts::SIGMA] = { 0 };
    uint64_t F2[consts::SIGMA] = { 0 };
    uint64_t F3[consts::SIGMA] = { 0 };
    uint64_t F4[consts::SIGMA] = { 0 };
    uint64_t F5[consts::SIGMA] = { 0 };
    uint64_t F6[consts::SIGMA] = { 0 };
    uint64_t F7[consts::SIGMA] = { 0 };
    size_t n8 = n & ~7ULL;
    size_t i;
#pragma simd
    for (i = 0; i < n8; i += 8) {
        freqs[in8[i + 0]]++;
        F1[in8[i + 1]]++;
        F2[in8[i + 2]]++;
        F3[in8[i + 3]]++;
        F4[in8[i + 4]]++;
        F5[in8[i + 5]]++;
        F6[in8[i + 6]]++;
        F7[in8[i + 7]]++;
    }
    /* the last few */
    while (i < n)
        freqs[in8[i++]]++;

    for (i = 0; i < consts::SIGMA; i++)
        freqs[i] += F1[i] + F2[i] + F3[i] + F4[i] + F5[i] + F6[i] + F7[i];

    assert(std::accumulate(freqs, freqs + consts::SIGMA, 0) == n);
}

void histogram(const uint8_t* in8, size_t n, uint64_t* freqs)
{
    size_t i = 0;
    /* the last few */
    while (i < n)
        freqs[in8[i++]]++;

    assert(std::accumulate(freqs, freqs + consts::SIGMA, 0ULL) == n);
}

void normalize_freqs(uint64_t* freqs, size_t n)
{
    double ratio = double(consts::M) / double(n);
    int32_t cur_table_size = 0;
    uint64_t most_freq_sym = 0;
    int64_t highest_freq = 0;
    for (size_t i = 0; i < consts::SIGMA; i++) {
        if (freqs[i] == 0)
            continue;
        if (highest_freq < freqs[i]) {
            most_freq_sym = i;
            highest_freq = freqs[i];
        }
        // scale down the frequency
        freqs[i] = freqs[i] * ratio + 0.499;
        if (freqs[i] == 0)
            freqs[i] = 1;
        cur_table_size += freqs[i];
    }
    int32_t adjust = int32_t(consts::M) - cur_table_size;
    if (adjust > 0) { // table too small currently
        // we just modify the prob of the most frequent symbol
        // as it will have the least impact on the overall error
        freqs[most_freq_sym] += adjust;
    } else if (adjust < 0) { // table too big currently
        if (highest_freq > -adjust) {
            // can we just fix it up with by decreasing the freq of
            // most frequent symbol
            freqs[most_freq_sym] += adjust;
        } else {
            // we have to distribute the adjustment to multiple
            // symbols.
            adjust += freqs[most_freq_sym] - 1;
            freqs[most_freq_sym] = 1;
            for (size_t j = 0; adjust && j < consts::SIGMA; j++) {
                if (freqs[j] < 2)
                    continue;
                int d = freqs[j] > -adjust;
                int64_t m = d ? adjust : 1 - freqs[j];
                freqs[j] += m;
                adjust -= m;
            }
        }
    }
    assert(std::accumulate(freqs, freqs + consts::SIGMA, 0ULL) == consts::M);
}