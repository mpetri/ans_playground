#include <cstdint>

#include "ans_util.hpp"
#include "util.hpp"

int main(int argc, char const* argv[])
{
    std::string file_name = argv[1];

    auto dat = read_file_contents(file_name);

    {
        uint64_t freqs[consts::SIGMA] = { 0 };
        benchmark bench_histogram8("histogram8", dat.size());
        histogram8(dat.data(), dat.size(), freqs);
    }

    {
        uint64_t freqs[consts::SIGMA] = { 0 };
        benchmark bench_histogram8("histogram", dat.size());
        histogram(dat.data(), dat.size(), freqs);
    }

    return 0;
}