#include <gtest/gtest.h>
#include <random>
#include <chrono>
#include "sais.hpp"

TEST(SAISTest, BasicTest)
{
    std::string seq{"TAAAGGGGCCCCCCAATATAATTTTGGGGCAAAGGGGCCCCCCAATAATTTTGGGGCAATAAAAAAATTTTTA"}; // the extra A denote $
    std::vector<std::size_t> ans_sa{72, 60, 61, 62, 63, 30, 1, 64, 31, 2, 57, 43, 14, 19, 46, 65, 32, 3, 58, 17, 44, 15, 20, 47, 66, 29, 56, 42, 13, 41, 12, 40, 11, 39, 10, 38, 9, 37, 8, 28, 55, 36, 7, 27, 54, 35, 6, 26, 53, 34, 5, 25, 52, 33, 4, 71, 59, 0, 18, 45, 16, 24, 51, 70, 23, 50, 69, 22, 49, 68, 21, 48, 67};

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

    SAIS sais;
    std::vector<std::size_t> sa(ans_sa.size());
    sais(seq, sa, 4);
    EXPECT_EQ(ans_sa, sa);
}

TEST(SAISTest, PerformanceTest)
{
    std::size_t length = 1024*1024;
    std::string seq;
    seq.reserve(length);
    std::vector<std::size_t> sa(length);
    SAIS sais;

    // generate radom DNA sequence 
    // {'A', 'C', 'G', 'T'} => {0, 1, 2, 3};
    std::default_random_engine eng;
    std::uniform_int_distribution<int> dist(0, 3); 
    for(std::size_t i = 0; i < length; i ++) {
        seq.push_back(dist(eng));
    }
    seq.push_back(0); // $

    auto start = std::chrono::high_resolution_clock::now();
    sais(seq, sa, 4);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "Random DNA seq size: 1MB\n"
              << "Suffix array construction time: " 
              << elapsed.count() << "s\n";
}
