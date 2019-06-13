#include <gtest/gtest.h>
#include <random>
#include <fstream>
#include <chrono>
#include <algorithm>
#include "saca_k.hpp"

template<class SEQ, class SA>
bool sa_is_correct(const SEQ& seq, const SA& sa, int k)
{
    int n = sa.size();
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

TEST(SACA_K, GetBuckets)
{
    SACA_K<std::string, std::vector<uint32_t>> sa_builder;
    std::vector<uint32_t> count {1, 2, 3, 4};
    decltype(count) bkt(count.size());
    decltype(count) ans_bkt_head {0, 1, 3, 6};
    decltype(count) ans_bkt_tail {0, 2, 5, 9};

    bool get_tail = false;
    sa_builder.get_buckets(count, bkt, get_tail);
    EXPECT_EQ(bkt, ans_bkt_head);

    std::fill(bkt.begin(), bkt.end(), 0);
    get_tail = true;
    sa_builder.get_buckets(count, bkt, get_tail);
    EXPECT_EQ(bkt, ans_bkt_tail);
}

TEST(SACA_K, GetLmsLen)
{
    // LMS can only be either $(sentinal) or can be expressed
    // as regex S+L+S
    std::vector<int> seq {1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0};
    SACA_K<decltype(seq), std::vector<uint32_t>> sa_builder;
    // $
    EXPECT_EQ(1, sa_builder.get_lms_len(
        seq.begin(), seq.size(), seq.size()-1));
    // SLS
    EXPECT_EQ(3, sa_builder.get_lms_len(
        seq.begin(), seq.size(), 1));
    // SSLS
    EXPECT_EQ(4, sa_builder.get_lms_len(
        seq.begin(), seq.size(), 3));
    // SLLS
    EXPECT_EQ(4, sa_builder.get_lms_len(
        seq.begin(), seq.size(), 6));
    // SSLLS
    EXPECT_EQ(5, sa_builder.get_lms_len(
        seq.begin(), seq.size(), 9));
}

TEST(SACA_K, GetSaOfLms)
{
    std::vector<int> seq {1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0};
    std::vector<uint32_t> sa(seq.size());
    SACA_K<decltype(seq), decltype(sa)> sa_builder;
    decltype(sa) sa_of_lms {13, 3, 9, 1, 6};
    uint32_t EMPTY { ((uint32_t)1) << (sizeof(uint32_t)*8-1) };
    int num_lms = sa_of_lms.size();

    // level0
    for (auto i = num_lms; i < seq.size(); i++)
        sa_of_lms.push_back(0);
    // assign sa1
    sa[0] = 4;
    sa[1] = 1;
    sa[2] = 3;
    sa[3] = 0;
    sa[4] = 2;
    sa_builder.get_sa_of_lms(
        seq.begin()
      , sa.begin()
      , sa.begin() + seq.size() - num_lms
      , seq.size()
      , num_lms
      , 0);
    EXPECT_EQ(sa, sa_of_lms);

    // level1
    std::fill(sa_of_lms.begin() + num_lms, sa_of_lms.end(), EMPTY);
    sa[0] = 4;
    sa[1] = 1;
    sa[2] = 3;
    sa[3] = 0;
    sa[4] = 2;
    sa_builder.get_sa_of_lms(
        seq.begin()
      , sa.begin()
      , sa.begin() + seq.size() - num_lms 
      , seq.size()
      , num_lms
      , 1);
    EXPECT_EQ(sa, sa_of_lms);
}

// TEST(SACA_K, InduceSAL0)
// {
//     // +induce_sal0(const SEQ_ITR seq,SaItr sa,std::vector<Index> & bkt,std::vector<Index> & count,Index n,bool suffix)
//
// }
//
// TEST(SACA_K, InduceSAS0)
// {
//     // +induce_sas0(const SEQ_ITR seq,SaItr sa,std::vector<Index> & bkt,std::vector<Index> & count,Index n,bool suffix)
//
// }
//
// TEST(SACA_K, InduceSAL1)
// {
//     // +induce_sal1(const SEQ_ITR seq,SaItr sa,Index n,bool suffix)
//
// }
//
// TEST(SACA_K, InduceSAS1)
// {
//     // +induce_sas1(const SEQ_ITR seq,SaItr sa,Index n,bool suffix)
//
// }

// TEST(SACA_K, NameSubstr)
// {
//     // +name_substr(const SEQ_ITR seq,SaItr sa,SaItr s1,Index n,Index m,Index n1,Index level)
//
// }

// TEST(SACA_K, PutLmsSubstr0)
// {
//     // +put_lms_substr0(const SEQ_ITR seq,SaItr sa,std::vector<Index> & bkt,std::vector<Index> & count,Index n)
//
// }
//
// TEST(SACA_K, PutLmsSubstr1)
// {
//     // +put_lms_substr1(const SEQ_ITR seq,SaItr sa,Index n)
//
// }
//
// TEST(SACA_K, PutSuffix0)
// {
//     // +put_suffix0(const SEQ_ITR seq,SaItr sa,std::vector<Index> & bkt,std::vector<Index> & count,Index n,Index n1)
//
// }
//
// TEST(SACA_K, PutSuffix1)
// {
//     // +put_suffix1(const SEQ_ITR seq,SaItr sa,Index n1)
//
// }

TEST(SACA_K, IntegrationOneLevel)
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
    EXPECT_TRUE(sa_is_correct(seq, sa, 3));
}

TEST(SACA_K, IntegrationThreeLevel)
{
    // level 0, size 5000
    // level 1, size 1425
    // level 2, size 467

    std::ifstream ifs("../unit_test/data/three-level-seq.txt");
    std::string seq;
    std::getline(ifs, seq);
    for (auto& elem : seq)
        elem -= 48;

    std::vector<uint32_t> sa(seq.size());
    SACA_K<decltype(seq), decltype(sa)> sa_builder;
    sa_builder.build(seq, sa, 4);
    EXPECT_TRUE(sa_is_correct(seq, sa, 4));
}

TEST(SACA_K, EasyTwoLevel0)
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
    EXPECT_TRUE(sa_is_correct(seq, sa, 4));
}

TEST(SACA_K, EasyTwoLevel1)
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
    EXPECT_TRUE(sa_is_correct(seq, sa, 3));
}

TEST(SACA_K, EasyTwoLevel2)
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
    EXPECT_TRUE(sa_is_correct(seq, sa, alphabet_size));
}

TEST(SACA_K, PerformanceTest1MB)
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
    EXPECT_TRUE(sa_is_correct(seq, sa, 4));
}
