#include <gtest/gtest.h>
#include "sais.hpp"

using SeqType   = std::string;
using IndexType = uint32_t;
using SAType    = std::vector<IndexType>;
using Sais      = SAIS<SeqType, SAType::iterator>;

SeqType seq{"mississippii"}; // extra i denode $
SAType ans_sa{11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2};
auto mapper = 
[](char base) 
{
    switch (base)
    {
        case 'i': return 0;
        case 'm': return 1;
        case 'p': return 2;
        case 's': return 3;
        default: 
            std::cerr << "unknown character, ascii code: " 
                      << (int)base << std::endl;
            throw std::runtime_error("unknown character");
    }
};

TEST(SAISTest, Constructor)
{
    SAType sa(seq.size());
    Sais sais(seq, sa.begin(), 4, mapper);
    EXPECT_EQ(sa, ans_sa);
}
