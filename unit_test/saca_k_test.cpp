#include <gtest/gtest.h>
#include <random>
#include <chrono>
#include <algorithm>
#include "saca_k.hpp"
// The test case is designed to reach no recursion
TEST(SACA_KTest, Random)
{
    std::string seq {"aatcgaaggtcgtaaggacacggttgagcgttcagcgtta"}; // include $
    std::vector<uint32_t> ans_sa {39, 13, 5, 0, 17, 19, 33, 26, 14, 6, 1, 18, 32, 3, 20, 10, 35, 28, 4, 16, 25, 34, 27, 15, 7, 21, 11, 8, 36, 29, 22, 38, 12, 31, 2, 9, 24, 37, 30, 23}; 
    std::vector<uint32_t> sa(ans_sa.size());

    std::transform(seq.begin(), seq.end(), seq.begin(),
        [](char base) 
        {
            switch (base)
            {
                case 'a': return 0;
                case 'c': return 1;
                case 'g': return 2;
                case 't': return 3;
                default: 
                    throw std::runtime_error(
                        "unknown character, ascii code: " 
                      + std::to_string(base));
            }
        });

    SACA_K<decltype(seq), decltype(sa)> sa_builder;
    sa_builder.build(seq, sa, 4);
    // std::for_each(sa.begin(), sa.end(), [](auto& elem){ std::cerr << elem << " "; });
    EXPECT_EQ(ans_sa, sa);
}

// The test case is designed to reach no recursion
TEST(SACA_KTest, RecursionLevel0)
{
    std::string seq {"bananaa"}; // include $
    std::vector<uint32_t> ans_sa {6, 5, 3, 1, 0, 4, 2};
    std::vector<uint32_t> sa(ans_sa.size());

    std::transform(seq.begin(), seq.end(), seq.begin(),
        [](char base) 
        {
            switch (base)
            {
                case 'a': return 0;
                case 'b': return 1;
                case 'n': return 2;
                default: 
                    throw std::runtime_error(
                        "unknown character, ascii code: " 
                      + std::to_string(base));
            }
        });

    SACA_K<decltype(seq), decltype(sa)> sa_builder;
    sa_builder.build(seq, sa, 3);
    EXPECT_EQ(ans_sa, sa);
}

// The test case is designed to reach 1 recursion
TEST(SACA_KTest, RecursionLevel1)
{
    std::string seq {"banaananana"}; // include $
    std::vector<uint32_t> ans_sa {10, 3, 8, 1, 6, 4, 0, 9, 2, 7, 5}; 
    std::vector<uint32_t> sa(ans_sa.size());

    std::transform(seq.begin(), seq.end(), seq.begin(),
        [](char base) 
        {
            switch (base)
            {
                case 'a': return 0;
                case 'b': return 1;
                case 'n': return 2;
                default: 
                    throw std::runtime_error(
                        "unknown character, ascii code: " 
                      + std::to_string(base));
            }
        });

    SACA_K<decltype(seq), decltype(sa)> sa_builder;
    sa_builder.build(seq, sa, 3);
    EXPECT_EQ(ans_sa, sa);
}

// A more complecate and complete test case
TEST(SACA_KTest, BasicTest)
{
    std::string seq{"TAAAGGGGCCCCCCAATATAATTTTGGGGCAAAGGGGCCCCCCAATAATTTTGGGGCAATAAAAAAATTTTTA"}; // the extra A denote $
    std::vector<uint32_t> ans_sa{72, 60, 61, 62, 63, 30, 1, 64, 31, 2, 57, 43, 14, 19, 46, 65, 32, 3, 58, 17, 44, 15, 20, 47, 66, 29, 56, 42, 13, 41, 12, 40, 11, 39, 10, 38, 9, 37, 8, 28, 55, 36, 7, 27, 54, 35, 6, 26, 53, 34, 5, 25, 52, 33, 4, 71, 59, 0, 18, 45, 16, 24, 51, 70, 23, 50, 69, 22, 49, 68, 21, 48, 67};

    std::transform(seq.begin(), seq.end(), seq.begin(),
        [](char base) 
        {
            switch (base)
            {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: 
                    throw std::runtime_error(
                        "unknown character, ascii code: " 
                      + std::to_string(base));
            }
        });

    int alphabet_size = 4; // A T G C
    std::vector<uint32_t> sa(ans_sa.size());
    SACA_K<decltype(seq), decltype(sa)> sa_builder;
    sa_builder.build(seq, sa, alphabet_size);
    EXPECT_EQ(ans_sa, sa);
}

TEST(SACA_KTest, PerformanceTest1MB)
{
    uint32_t length = 1024*1024;
    std::vector<char> seq(length);
    std::vector<uint32_t> sa(length);
    SACA_K<decltype(seq), decltype(sa)> sa_builder;

    // generate radom DNA sequence 
    // {'A', 'C', 'G', 'T'} => {0, 1, 2, 3};
    std::default_random_engine eng;
    std::uniform_int_distribution<char> dist(0, 3); 
    auto get_random_base = [&eng, &dist]() { return dist(eng); };
    // seq[length] denote $
    std::generate(seq.begin(), seq.end()-1, get_random_base);

    auto start = std::chrono::high_resolution_clock::now();
    sa_builder.build(seq, sa, 4);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "Random DNA seq size: 1MB\n"
              << "Suffix array construction time: " 
              << elapsed.count() << "s\n";
}

