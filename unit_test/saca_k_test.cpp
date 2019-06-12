#include <gtest/gtest.h>
#include <random>
#include <chrono>
#include <algorithm>
#include "saca_k.hpp"

template<class SEQ, class SA>
bool sa_is_correct(const SEQ& seq, const SA& sa, int n, int k)
{
    // check elem in sa is in range [0, n)
    for (const auto& elem : sa)
        if (elem < 0 || elem >= n)
        {
            std::cerr << "sa incorrect: out of range\n";
            return false;
        }

    // check first character
    for (auto i = 1; i < n; i++)
        if (seq[sa[i-1]] > seq[sa[i]])
        {
            std::cerr << "sa incorrect: first character\n";
            return false;
        }

    // check suffix:
    // For all character in k, if sa[a, b] contains the suffixes
    // starting with the character c, 
    // then sa[a]+1, sa[a+1]+1,..., sa[b]+1 occur in sa in this order

    SA bkt(k);
    for (auto i = n-1; i >= 0; i--)
        bkt[seq[sa[i]]] = i; // get head bucket
    bkt[seq[n-1]]++; // skip $

    for (const auto& elem : sa)
        if (elem >= 1)
        {
            auto c = seq[elem - 1];
            if ((sa[bkt[c]] + 1) != elem)
            {
                std::cerr << "sa incorrect: suffix\n";
                return false;
            }
            bkt[c]++;
        }

    return true;
}

TEST(SACA_KTest, Random)
{
    std::string seq {"aatcgaaggtcgtaaggacacggttgagcgttcagcgtta"}; // include $
    std::vector<uint32_t> sa(seq.size());

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
    EXPECT_TRUE(sa_is_correct(seq, sa, sa.size(), 4));
}

// The test case is designed to reach no recursion
TEST(SACA_KTest, RecursionLevel0)
{
    std::string seq {"bananaa"}; // include $
    std::vector<uint32_t> sa(seq.size());

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
    EXPECT_TRUE(sa_is_correct(seq, sa, sa.size(), 3));
}

// The test case is designed to reach 1 recursion
TEST(SACA_KTest, RecursionLevel1)
{
    std::string seq {"banaananana"}; // include $
    std::vector<uint32_t> sa(seq.size());

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
    EXPECT_TRUE(sa_is_correct(seq, sa, sa.size(), 3));
}

// A more complecate and complete test case
TEST(SACA_KTest, BasicTest)
{
    std::string seq{"TAAAGGGGCCCCCCAATATAATTTTGGGGCAAAGGGGCCCCCCAATAATTTTGGGGCAATAAAAAAATTTTTA"}; // the extra A denote $
    std::vector<uint32_t> sa(seq.size());
    int alphabet_size = 4; // A T G C

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

    SACA_K<decltype(seq), decltype(sa)> sa_builder;
    sa_builder.build(seq, sa, alphabet_size);
    EXPECT_TRUE(sa_is_correct(seq, sa, sa.size(), alphabet_size));
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
    std::generate(
        seq.begin()
      , seq.end()-1 // seq[length] denote $
      , [&eng, &dist](){ return dist(eng); });

    auto start = std::chrono::high_resolution_clock::now();
    sa_builder.build(seq, sa, 4);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "Random DNA seq size: 1MB\n"
              << "Suffix array construction time: " 
              << elapsed.count() << "s\n";
    EXPECT_TRUE(sa_is_correct(seq, sa, sa.size(), 4));
}

