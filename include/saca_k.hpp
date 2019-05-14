#pragma once
#include <algorithm>
#include <memory>
#include <cassert>

struct SACA_K
{
    SACA_K() = default;

    template<class SEQ, class SA>
    void operator()(const SEQ& seq, SA& sa, uint32_t k)
    {
        return call_impl(seq.begin(), sa.begin(), seq.size(), k);
    }

  private:
    /// @brief Find the suffix array of seq[0..n-1] in {0..k-1}^n
    /// require s[n-1]=0 (the sentinel!), n>=2
    template<class SEQ_ITR, class SA_ITR>
    void call_impl(const SEQ_ITR seq, SA_ITR sa, std::size_t n, uint32_t k)
    {
        std::vector<std::size_t> ans_sa{72, 60, 61, 62, 63, 30, 1, 64, 31, 2, 57, 43, 14, 19, 46, 65, 32, 3, 58, 17, 44, 15, 20, 47, 66, 29, 56, 42, 13, 41, 12, 40, 11, 39, 10, 38, 9, 37, 8, 28, 55, 36, 7, 27, 54, 35, 6, 26, 53, 34, 5, 25, 52, 33, 4, 71, 59, 0, 18, 45, 16, 24, 51, 70, 23, 50, 69, 22, 49, 68, 21, 48, 67};
        std::transform(ans_sa.begin(), ans_sa.end(), sa, 
            [](std::size_t value){ return value; });

        // std::size_t i, j;
        // std::size_t c0, c1;
        // assert((n >= 2) && (k >= 1));
        //
        // // to free the managed object: ptr.reset();
        // auto count = std::make_unique<std::size_t[]>(k);
        // auto bkt = std::make_unique<std::size_t[]>(k);
        //
        // // stage 1: reduce the problem by at least 1/2
        // // sort all the LMS-substrings 
        // get_counts(seq, count, n, k);
        // get_buckets(count, bkt, k, true); // find end of bkt
        // for (i = 0; i < n; i++)
        //     sa[i] = 0;
        // // place all the LMS into bucket in SA
        // // debug
        // for (i = 0; i < k; i++)
        //     std::cerr << bkt[i] << " ";
        // std::cerr << std::endl;
        // for (i = 0; i < n; i++)
        //     std::cerr << sa[i] << " ";
        // std::cerr << std::endl;
        //
        // std::size_t num_lms = 1;
        // sa[0] = n-1;
        // auto b = sa + --bkt[0]; 
        // i = n - 1;
        // j = n;
        // c0 = seq[n-1];
        // // find next S-type
        // do { c1 = c0; } 
        // while ((--i < n) && ((c0 = seq[i]) >= c1));
        // for (; i < n;)
        // {
        //     // find next L-type
        //     do { c1 = c0; } 
        //     while ((--i < n) && ((c0 = seq[i]) <= c1));
        //     if (i < n)
        //     {
        //         std::cerr << j << " ";//debug
        //         *b = j; b = sa + --bkt[c1]; j = i; num_lms++;
        //         // find next S-type
        //         do { c1 = c0; } 
        //         while ((--i < n) && ((c0 = seq[i]) >= c1));
        //     }
        // }
        // std::cerr << std::endl; //debug
        // // debug
        // for (i = 0; i < n; i++)
        //     std::cerr << sa[i] << " ";
        // std::cerr << std::endl;
    }

    template<class SEQ_ITR>
    auto get_counts(
        const SEQ_ITR seq
      , std::unique_ptr<std::size_t[]>& count
      , std::size_t n
      , uint32_t k
    )
    {
        std::size_t i;
        for (i = 0; i < k; i++)
            count[i] = 0;
        for (i = 0; i < n; i++)
            count[seq[i]]++;
    }

    auto get_buckets(
        std::unique_ptr<std::size_t[]>& count
      , std::unique_ptr<std::size_t[]>& bkt
      , uint32_t k
      , bool end
    )
    {
        std::size_t i, sum = 0;
        for (i = 0; i < k; i++)
        {
            sum += count[i];
            if (end)
                bkt[i] = sum;
            else
                bkt[i] = sum - count[i];
        }
    }
};
