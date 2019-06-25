#include <gtest/gtest.h>
#include "fm_index.hpp"
#include "saca_k.hpp"

using ::testing::TestWithParam;
using ::testing::Values;
using SeqType = std::string;
using IndexType = uint8_t;

class IntegrationTest : public TestWithParam<int>
{
  protected:
    void SetUp() override
    { sample_step = GetParam(); }

    // Sequence: repeated short and long lms. (repeated short lms
    //           will make sure it has hash collistion and enter 
    //           sorter, which can test hash and sorter)
    //
    // short(12 mer)
    // a:AATA * 3
    // b:AATTTTGGGGCA * 2
    // c:ATA *1
    //
    // long
    // d:AAAGGGGCCCCCCA * 2
    // e:AAAAAAATTTTT$ * 1
    //
    // final seq: Tdacbdabae
    // TAAAGGGGCCCCCCAATATAATTTTGGGGCAAAGGGGCCCCCCAATAATTTTGGGGCAATAAAAAAATTTTT$
    //  *            *  * *          *            *  *          *  *           * (lms)
    //
    // distinct lms location: 1, 17, 30, 46, 57, 60, 72
    // lms location: 1, 14, 17, 19, 30, 43, 46, 57, 60, 72
    // seq size: 73
    // number of lms: 9 + 1($) = 10
    // number of distinct lms: 5 + 1(repeated long lms: d) + 1($) = 7
    SeqType seq{"TAAAGGGGCCCCCCAATATAATTTTGGGGCAAAGGGGCCCCCCAATAATTTTGGGGCAATAAAAAAATTTTTA"}; // the extra A denote $
    std::vector<IndexType> sa{72, 60, 61, 62, 63, 30, 1, 64, 31, 2, 57, 43, 14, 19, 46, 65, 32, 3, 58, 17, 44, 15, 20, 47, 66, 29, 56, 42, 13, 41, 12, 40, 11, 39, 10, 38, 9, 37, 8, 28, 55, 36, 7, 27, 54, 35, 6, 26, 53, 34, 5, 25, 52, 33, 4, 71, 59, 0, 18, 45, 16, 24, 51, 70, 23, 50, 69, 22, 49, 68, 21, 48, 67};
    std::vector<std::vector<IndexType>> lf_map{
        {1, 25, 39, 55}
      , {1, 25, 39, 56}
      , {1, 25, 39, 57}
      , {2, 25, 39, 57}
      , {3, 25, 39, 57}
      , {4, 25, 39, 57}
      , {4, 26, 39, 57}
      , {4, 26, 39, 58}
      , {5, 26, 39, 58}
      , {6, 26, 39, 58}
      , {7, 26, 39, 58}
      , {7, 27, 39, 58}
      , {7, 28, 39, 58}
      , {7, 29, 39, 58}
      , {7, 29, 39, 59}
      , {7, 29, 39, 60}
      , {8, 29, 39, 60}
      , {9, 29, 39, 60}
      , {10, 29, 39, 60}
      , {11, 29, 39, 60}
      , {11, 29, 39, 61}
      , {12, 29, 39, 61}
      , {13, 29, 39, 61}
      , {14, 29, 39, 61}
      , {15, 29, 39, 61}
      , {16, 29, 39, 61}
      , {16, 29, 40, 61}
      , {16, 29, 41, 61}
      , {16, 30, 41, 61}
      , {16, 31, 41, 61}
      , {16, 32, 41, 61}
      , {16, 33, 41, 61}
      , {16, 34, 41, 61}
      , {16, 35, 41, 61}
      , {16, 36, 41, 61}
      , {16, 37, 41, 61}
      , {16, 38, 41, 61}
      , {16, 39, 41, 61}
      , {16, 39, 42, 61}
      , {16, 39, 43, 61}
      , {16, 39, 44, 61}
      , {16, 39, 45, 61}
      , {16, 39, 46, 61}
      , {16, 39, 47, 61}
      , {16, 39, 48, 61}
      , {16, 39, 49, 61}
      , {16, 39, 50, 61}
      , {16, 39, 51, 61}
      , {16, 39, 52, 61}
      , {16, 39, 53, 61}
      , {16, 39, 54, 61}
      , {16, 39, 55, 61}
      , {16, 39, 55, 62}
      , {16, 39, 55, 63}
      , {17, 39, 55, 63}
      , {18, 39, 55, 63}
      , {18, 39, 55, 64}
      , {19, 39, 55, 64}
      , {19, 39, 55, 64}
      , {20, 39, 55, 64}
      , {21, 39, 55, 64}
      , {22, 39, 55, 64}
      , {22, 39, 55, 65}
      , {22, 39, 55, 66}
      , {22, 39, 55, 67}
      , {22, 39, 55, 68}
      , {22, 39, 55, 69}
      , {22, 39, 55, 70}
      , {22, 39, 55, 71}
      , {22, 39, 55, 72}
      , {22, 39, 55, 73}
      , {23, 39, 55, 73}
      , {24, 39, 55, 73}
    }; // lf_mapping(index, char)
    int sample_step; 
};

TEST_P(IntegrationTest, Constructor)
{
    auto map = 
    [](char base) 
    {
        switch (base)
        {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: 
                std::cerr << "unknown character, ascii code: " 
                          << (int)base << std::endl;
                throw std::runtime_error("unknown character");
        }
    };

    FmIndex<SeqType, IndexType, 2, SACA_K> fm_index(
        seq, map, sample_step);

    for (auto i = 0; i < seq.size(); i++)
    {
        EXPECT_EQ(fm_index.get_location(i), sa[i]);
        EXPECT_EQ(fm_index.lf_mapping(i, 'A'), lf_map[i][0]);
        EXPECT_EQ(fm_index.lf_mapping(i, 'C'), lf_map[i][1]);
        EXPECT_EQ(fm_index.lf_mapping(i, 'G'), lf_map[i][2]);
        EXPECT_EQ(fm_index.lf_mapping(i, 'T'), lf_map[i][3]);
    }
}

// Parameterized test: pass in sample_step
INSTANTIATE_TEST_CASE_P(DifferentSampleRate, IntegrationTest
    , Values(1, 2, 4, 8, 16, 32));
