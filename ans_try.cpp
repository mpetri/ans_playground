#include <cstdint>

#include "ans_util.hpp"
#include "util.hpp"

size_t ans_compress(const uint8_t* data, size_t n, uint8_t* out)
{
    size_t written_bytes = 0;
    benchmark bench_histogram8("ans_compress", n);

    uint64_t freqs[consts::SIGMA] = { 0 };
    histogram8(data, n, freqs);

    normalize_freqs(freqs, n);

    written_bytes += encode_prelude_rle_vbyte(freqs, out);

    return written_bytes;
}

int main(int argc, char const* argv[])
{
    std::string file_name = argv[1];

    auto dat = read_file_contents(file_name);

    std::vector<uint8_t> out(dat.size());
    auto wb = ans_compress(dat.data(), dat.size(), out.data());
    out.resize(wb);

    double cr = double(wb) / double(dat.size());
    std::cout << "compression ratio = " << cr << std::endl;

    return 0;
}