#pragma once

#include <cstdint>

namespace consts {
uint64_t SIGMA = 256;
uint64_t TABLE_SHIFT = 12;
uint64_t TABLE_SIZE = 1ULL << TABLE_SHIFT;
uint8_t PADDING = 8;
}

/* do this so the loop can be nicely vectorized */
void histogram8(const uint8_t* in8, std::size_t n, uint64_t* freqs)
{
    uint64_t F1[consts::SIGMA] = { 0 };
    uint64_t F2[consts::SIGMA] = { 0 };
    uint64_t F3[consts::SIGMA] = { 0 };
    uint64_t F4[consts::SIGMA] = { 0 };
    uint64_t F5[consts::SIGMA] = { 0 };
    uint64_t F6[consts::SIGMA] = { 0 };
    uint64_t F7[consts::SIGMA] = { 0 };
    std::size_t n8 = n & ~7ULL;
    std::size_t i;
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
}

/* do this so the loop can be nicely vectorized */
void histogram(const uint8_t* in8, std::size_t n, uint64_t* freqs)
{
    std::size_t i = 0;
    /* the last few */
    while (i < n)
        freqs[in8[i++]]++;
}
