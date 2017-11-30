#pragma once

#include <chrono>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>

using namespace std::chrono;

template<typename T>
std::string print_time_diff(T prior, T latter)
{
    namespace sc = std::chrono;
    auto time_ns = sc::duration_cast<sc::nanoseconds>(latter - prior).count();
   	double dtime_ns = time_ns;
    if(dtime_ns < 1000) {
    	return std::to_string(dtime_ns) + " nsec";
    }
    if(dtime_ns < 1000,000) {
    	return std::to_string(dtime_ns/1000) + " usec";
    }
    if(dtime_ns < (1000*1000*1000)) {
    	return std::to_string(dtime_ns/1000000) + " msec";
    }
    return std::to_string(dtime_ns/(1000*1000*1000)) + " sec";
}

struct timer {
    high_resolution_clock::time_point start;
    std::string name;
    timer(const std::string& _n)
        : name(_n)
    {
        std::cerr << "START(" << name << ")" << std::endl;
        start = high_resolution_clock::now();
    }
    ~timer()
    {
        auto stop = high_resolution_clock::now();
        std::cerr << "STOP(" << name << ") - "
                  << duration_cast<milliseconds>(stop - start).count() / 1000.0f
                  << " sec" << std::endl;
    }
};

struct benchmark {
    high_resolution_clock::time_point start;
    std::string name;
    size_t total_bytes;
    benchmark(const std::string& _n,size_t tb)
        : name(_n), total_bytes(tb)
    {
        start = high_resolution_clock::now();
    }
    ~benchmark()
    {
        auto stop = high_resolution_clock::now();
        auto diff = stop - start;
        double bytes_per_ns = double(total_bytes) / double(diff.count());
        double mib_per_sec = bytes_per_ns * 953.6743;
        std::cerr << "BENCH(" << name << ") - TIME: "
                  << print_time_diff(start,stop)
                  << " - SPEED: " 
                  << mib_per_sec << " MiB/s" << std::endl;
    }
};

std::vector<uint8_t>
read_file_contents(std::string& file_name)
{
	FILE* f = fopen(file_name.c_str(), "rb");
    if (!f) {
    	std::cerr << "Can't open file '" << file_name << "'" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::vector<uint8_t> content;
    auto cur = ftell(f);
    fseek(f, 0, SEEK_END);
    auto end = ftell(f);
    size_t file_size = (end - cur);
    benchmark bench_read("read_file",file_size);
    content.resize(file_size);
    fseek(f, cur, SEEK_SET);
    size_t ret = fread(content.data(), 1, file_size, f);
    if (ret != file_size) {
        std::cerr << "Error reading file '" << file_name << "'" << std::endl;
        exit(EXIT_FAILURE);
    }
    fclose(f);
    return content;
}