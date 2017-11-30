#include <cstdint>

#include "ans_util.hpp"
#include "util.hpp"

void ans_compress(const uint8_t* data, size_t n)
{
    benchmark bench_histogram8("ans_compress", n);

    uint64_t freqs[consts::SIGMA] = { 0 };
    histogram8(data, n, freqs);

    normalize_freqs(freqs, n);
}

int main(int argc, char const* argv[])
{
    std::string file_name = argv[1];

    auto dat = read_file_contents(file_name);

    ans_compress(dat.data(), dat.size());

    return 0;
}