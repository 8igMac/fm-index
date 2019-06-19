#pragma once
#include <type_traits>
#include <algorithm>

template<class SEQ, class SA>
class SACA_K
{
    using SaItr = typename SA::iterator;
    using Index = typename SA::value_type;
    using SignedIndex = std::make_signed_t<Index>;
    const Index EMPTY { ((Index)1) << (sizeof(Index)*8-1) };

  public:
    SACA_K() = default;

    void build(const SEQ& seq, SA& sa, Index k)
    {
        return call_impl(
            seq.begin()
          , sa.begin()
          , seq.size()
          , k
          , seq.size());
    }

  // private:
    /// @brief Find the suffix array of seq[0..n-1] in {0..k-1}^n
    /// require seq[n-1]=0 (the sentinel!), n>=2
    template<class SEQ_ITR>
    void call_impl(
        const SEQ_ITR seq
      , SaItr sa
      , Index n // seq and sa size
      , Index k // the number of character used
      , Index m // maxima available space
      , Index level = 0
    )
    {
        // stage 1: reduce the problem by at least 1/2
        std::vector<Index> bkt(k), count(k);
        if (level == 0)
        {
            // count the number of each character
            std::for_each(seq, seq+n,
                [&count](auto& chr){ count[chr]++; });

            put_lms_substr0(seq, sa, bkt, count, n);
            induce_sal0(seq, sa, bkt, count, n, false);
            induce_sas0(seq, sa, bkt, count, n, false);
        }
        else
        {
            put_lms_substr1(seq, sa, n);
            induce_sal1(seq, sa, n, false);
            induce_sas1(seq, sa, n, false);
        }

        // now, all the LMS-substrings are sorted and stored 
        // sparsely in SA.

        // compact all the sorted substrings into the first n1 items
        // of SA. 2*n1 must be not larger than n (proveable)
        Index n1 = 0;
        std::for_each(sa, sa+n, 
            [&n1, &sa](auto& sa_value)
            {
                if (static_cast<SignedIndex>(sa_value) > 0)
                    sa[n1++] = sa_value;
            });

        auto sa1 = sa; // sa1: the first n1 elements in sa
        auto s1 = sa + m - n1; // s1: the last n1 element in sa
        Index name_count = name_substr(seq, sa, s1, n, m, n1);

        // stage 2: solve the reduced problem
        // recurse if names are not yet unique
        if (name_count < n1)
            call_impl(s1, sa1, n1, 0, m-n1, level+1);
        else
            std::for_each(s1, s1+n1
                , [&sa1, i = 0](auto& elem) mutable { sa1[elem] = i++; });

        // stage 3: induce SA(S) from SA(S1)
        get_sa_of_lms(seq, sa, s1, n, n1, level);
        if (level == 0)
        {
            put_suffix0(seq, sa, bkt, count, n, n1);
            induce_sal0(seq, sa, bkt, count, n, true);
            induce_sas0(seq, sa, bkt, count, n, true);
        }
        else 
        {
            put_suffix1(seq, sa, n1);
            induce_sal1(seq, sa, n, true);
            induce_sas1(seq, sa, n, true);
        }
    }

    void get_buckets(
        const std::vector<Index>& count
      , std::vector<Index>& bkt
      , bool end
    )
    {
        std::transform(count.begin(), count.end(), bkt.begin()
          , [&end, sum = 0](auto& chr_count) mutable
            {
                sum += chr_count;
                return (end) ? sum - 1 : sum - chr_count;
            });
    }

    template<class SEQ_ITR>
    Index get_lms_len(
        const SEQ_ITR seq
      , Index n
      , Index pos
    )
    {
        if (pos == n-1)
            return 1;

        Index dist, i = 1;
        // 2 consecutive LMS substr must be in th form: S+L+S+L+S
        while (1)
        {
            // break when seq[pos+i-1] is the first L-type we met
            if (seq[pos+i] < seq[pos+i-1])
                break;
            i++;
        }

        while (1)
        {
            // break when seq[pos+i-1] is the first S-type we met
            if (pos+i > n-1 || seq[pos+i] > seq[pos+i-1])
                break;
            // record LMS
            if (pos+i == n-1 || seq[pos+i] < seq[pos+i-1])
                dist = i;
            i++;
        }

        return dist + 1;
    }

    template<class SEQ_ITR>
    Index name_substr(
        const SEQ_ITR seq
      , SaItr sa
      , SaItr s1
      , Index n
      , Index m
      , Index n1
    )
    {
        // init name array buffer
        std::fill(sa+n1, sa+n, EMPTY);

        // scan to compute the interim s1
        Index name, name_counter = 0;
        Index pre_pos, pre_len = 0;
        for (auto i = 0; i < n1; i++)
        {
            auto diff = false;
            Index pos = sa[i];
            Index len = get_lms_len(seq, n, pos);
            if (len != pre_len)
                diff = true;
            else
                for (auto d = 0; d < len; d++)
                    if (pos+d == n-1 || pre_pos+d == n-1 || 
                        (seq[pos+d] != seq[pre_pos+d]))
                    {
                        diff = true;
                        break;
                    }

            if (diff)
            {
                name = i; 
                name_counter++;
                sa[name] = 1; // a new name
                pre_pos = pos;
                pre_len = len;
            }
            else
                sa[name]++; // count this name
                
            sa[n1 + pos/2] = name;
        }

        // compact the interim s1 sparsely stored in sa[n1, n-1] into
        // sa[m-n1, m-1]
        for (auto i = n-1, j = m-1; i >= n1; i--)
            if (sa[i] != EMPTY)
                sa[j--] = sa[i];

        // rename each S-type character of the interim s1 as the end
        // of its bucket to produce the final s1
        bool pre_type, cur_type = true;
        for (auto i = n1-1; i > 0; i--) 
        {
            auto ch = s1[i];
            auto pre_ch = s1[i-1];
            pre_type = (pre_ch < ch || (pre_ch == ch && cur_type));
            if (pre_type)
                s1[i-1] += sa[s1[i-1]] - 1;
            cur_type = pre_type;
        }

        return name_counter;
    }

    template<class SEQ_ITR>
    void get_sa_of_lms(
        const SEQ_ITR seq
      , SaItr sa
      , SaItr s1
      , Index n
      , Index n1
      , Index level
    )
    {
        // put LMS into s1
        Index j = n1-1;
        s1[j--] = n-1;
        bool pre_type, cur_type = false; // seq[n-2] must be L-type
        for (auto i = n-2; i > 0; i--)
        {
            pre_type = (seq[i-1] < seq[i] || 
                (seq[i-1] == seq[i] && cur_type));
            if (cur_type && !pre_type)
                s1[j--] = i;
            cur_type = pre_type;
        }

        // get suffix array of LMS
        std::for_each(sa, sa+n1
          , [&s1](auto& sa_value){ 
          sa_value = s1[sa_value]; });

        // init sa[n1..n-1]
        if (level)
            std::fill(sa+n1, sa+n, EMPTY);
        else
            std::fill(sa+n1, sa+n, 0);
    }

    template<class SEQ_ITR>
    void put_lms_substr0(
        const SEQ_ITR seq
      , SaItr sa
      , std::vector<Index>& bkt
      , std::vector<Index>& count
      , Index n
    )
    {
        // find end of each bucket
        get_buckets(count, bkt, true);

        // clear sa
        std::fill(sa, sa+n, 0);

        bool pre_type, cur_type = false; // seq[n-2] must be L-type
        for (auto i = n-2; i > 0; i--)
        {
            pre_type = (seq[i-1] < seq[i] || 
                (seq[i-1] == seq[i] && cur_type));
            if (cur_type && !pre_type)
                sa[bkt[seq[i]]--] = i;
            cur_type = pre_type;
        }

        sa[0] = n-1; // set the single sentinel LMS substr
    }

    template<class SEQ_ITR>
    void induce_sal0(
        const SEQ_ITR seq
      , SaItr sa
      , std::vector<Index>& bkt
      , std::vector<Index>& count
      , Index n
      , bool suffix
    )
    {
        get_buckets(count, bkt, false); // find the head of bucket
        bkt[0]++; // skip $
        for (auto i = 0; i < n; i++)
            if (sa[i] > 0)
            {
                auto j = sa[i] - 1;
                if (seq[j] >= seq[j+1])
                {
                    sa[ bkt[seq[j]]++ ] = j;
                    if (!suffix && i>0)
                        sa[i] = 0;
                }
            }
    }

    template<class SEQ_ITR>
    void induce_sas0(
        const SEQ_ITR seq
      , SaItr sa
      , std::vector<Index>& bkt
      , std::vector<Index>& count
      , Index n
      , bool suffix
    )
    {
        get_buckets(count, bkt, true); // find the end of bucket
        for (auto i = n-1; i > 0; i--)
            if (sa[i] > 0)
            {
                auto j = sa[i] - 1;
                if (seq[j] < seq[j+1] || 
                    (seq[j] == seq[j+1] && bkt[seq[j]] < i))
                {
                    sa[ bkt[seq[j]]-- ] = j;
                    if (!suffix)
                        sa[i] = 0;
                }
            }
    }

    template<class SEQ_ITR>
    void put_suffix0(
        const SEQ_ITR seq
      , SaItr sa
      , std::vector<Index>& bkt
      , std::vector<Index>& count
      , Index n
      , Index n1
    )
    {
        get_buckets(count, bkt, true); // find the end of bucket

        // put the suffix into their bucket
        for (auto i = n1-1; i > 0; i--)
        {
            sa[ bkt[seq[sa[i]]]-- ] = sa[i];
            sa[i] = 0;
        }
        sa[0] = n-1; // set the single sentinel suffix
    }

    template<class SEQ_ITR>
    void put_lms_substr1(
        const SEQ_ITR seq
      , SaItr sa
      , Index n
    )
    {
        std::fill(sa, sa+n, EMPTY);

        Index c, pre_c = seq[n-2];
        bool type, pre_type = false;
        for (auto i = n-2; i > 0; i--)
        {
            c = pre_c;
            type = pre_type;
            pre_c = seq[i-1];
            pre_type = (pre_c < c) 
                || (pre_c == c && type);

            if (type && !pre_type)
            {
                if (static_cast<SignedIndex>(sa[c]) >= 0)
                {
                    // SA[c] is borrowed by the right
                    //   neighbor bucket.
                    // shift-right the items in the
                    //   right neighbor bucket.
                    Index cnt_pos;
                    // find counter
                    for (cnt_pos = c+1
                      ; static_cast<SignedIndex>(sa[cnt_pos]) >= 0
                      ; cnt_pos++); 
                    std::move_backward(
                        sa+c
                      , sa+cnt_pos
                      , sa+cnt_pos+1);
                    sa[c] = EMPTY;
                }

                auto d = static_cast<SignedIndex>(sa[c]);
                if (d == EMPTY) // sa[c] is empty
                {
                    if (sa[c-1] == EMPTY)
                    {
                        sa[c] = -1; // init the counter
                        sa[c-1] = i;
                    }
                    else
                        sa[c] = i; // a size-1 bucket
                }
                else // sa[c] is reused as a counter
                {
                    auto pos = c+d-1;
                    if (sa[pos] != EMPTY)
                    {
                        // we are running into the left
                        //   neighbor bucket.
                        // shift-right one step the items 
                        //   of bucket(SA, S, i).
                        pos++;
                        std::move_backward(sa+pos, sa+c, sa+c+1);
                    }
                    else
                        sa[c]--;
                    sa[pos] = i;
                }
            }
        }

        // scan to shift-right the items in each bucket
        //   with its head being reused as a counter.
        for (auto i = n-1; i > 0; i--)
        {
            auto j = static_cast<SignedIndex>(sa[i]);
            // is sa[i] a counter?
            if (j < 0 && j != EMPTY) 
            {
                std::move_backward(sa+i+j, sa+i, sa+i+1);
                sa[i+j] = EMPTY;
            }
        }

        // put the single sentinel LMS-substring.
        sa[0] = n-1;
    }

    template<class SEQ_ITR>
    void induce_sal1(
        const SEQ_ITR seq
      , SaItr sa
      , Index n
      , bool suffix
    )
    {
        for (auto i = 0, step = 1; i < n; i += step, step = 1)
        {
            auto j = sa[i] - 1;
            if (static_cast<SignedIndex>(sa[i]) <= 0)
                continue;
            auto c = seq[j], c1 = seq[j+1];
            bool is_L_type = (c >= c1);
            if (!is_L_type)
                continue;

            auto d = static_cast<SignedIndex>(sa[c]);
            if (d >= 0)
            {
                // SA[c] is borrowed by the left
                //   neighbor bucket.
                // shift-left the items in the
                //   left neighbor bucket.
                Index cnt_pos;
                // find counter
                for (cnt_pos = c-1
                  ; static_cast<SignedIndex>(sa[cnt_pos]) >= 0 
                 || sa[cnt_pos] == EMPTY
                  ; cnt_pos--);
                std::move(sa+cnt_pos+1, sa+c+1, sa+cnt_pos);
                if (cnt_pos < i)
                    step = 0; //TODO: why do we need step?
                d = EMPTY;
            }
            
            if (d == EMPTY) // sa[c] is empty
                if (c < n-1 && sa[c + 1] == EMPTY)
                {
                    sa[c] = -1; // init the counter
                    sa[c+1] = j;
                }
                else
                    sa[c] = j; // a size 1 bucket
            else // sa[c] is reused as a counter
            {
                auto pos = c-d+1;
                if (pos > n-1 || sa[pos] != EMPTY)
                {
                    // we are running into the right
                    //   neighbor bucket.
                    // shift-left one step the items
                    //   of bucket(SA, S, j).
                    std::move(sa+c+1, sa+pos, sa+c);
                    pos--;
                    if (c < i)
                        step = 0;
                }
                else
                    sa[c]--;
                sa[pos] = j;
            }

            decltype(c1) c2;
            // is seq[sa[i]] L-type?
            bool is_L1 = (j+1 < n-1) && 
                (c1 > (c2 = seq[j+2]) || (c1 == c2 && c1 < i));
            if ((!suffix || !is_L1) && i > 0)
            {
                auto i1 = (step == 0) ? i-1 : i;
                sa[i1] = EMPTY;
            }
        }

        // scan to shift-left the items in each bucket 
        //   with its head being reused as a counter.
        for (auto i = 1; i < n; i++)
        {
            auto j = static_cast<SignedIndex>(sa[i]);
            if (j < 0 && j != EMPTY) // is sa[i] a counter?
            {
                std::move(sa+i+1, sa+i+1-j, sa+i);
                sa[i-j] = EMPTY;
            }
        }
    }

    template<class SEQ_ITR>
    void induce_sas1(
        const SEQ_ITR seq
      , SaItr sa
      , Index n
      , bool suffix
    )
    {
        for (auto i = n-1, step = 1; i > 0; i -= step, step = 1)
        {
            auto j = sa[i]-1;
            if (static_cast<SignedIndex>(sa[i]) <= 0)
                continue;
            auto c = seq[j], c1 = seq[j+1];
            bool is_s_type = (c < c1) || (c == c1 && c > i);
            if (!is_s_type)
                continue;

            auto d = static_cast<SignedIndex>(sa[c]);
            if (d >= 0)
            {
                // SA[c] is borrowed by the right
                //   neighbor bucket.
                // shift-right the items in the
                //   right neighbor bucket.
                Index cnt_pos;
                // find counter
                for (cnt_pos = c+1
                  ; static_cast<SignedIndex>(sa[cnt_pos]) >= 0
                  ; cnt_pos++); 
                std::move_backward(
                    sa+c
                  , sa+cnt_pos
                  , sa+cnt_pos+1);
                if (cnt_pos > i)
                    step = 0;
                d = EMPTY;
            }

            if (d == EMPTY) // sa[c] is empty
            {
                if (sa[c-1] == EMPTY)
                {
                    sa[c] = -1; // init the counter
                    sa[c-1] = j;
                }
                else
                    sa[c] = j;
            }
            else // sa[c] is reused as a counter
            {
                auto pos = c+d-1;
                if (sa[pos] != EMPTY)
                {
                    // we are running into the left
                    //   neighbor bucket.
                    // shift-right one step the items 
                    //   of bucket(SA, S, j).
                    std::move_backward(
                        sa+pos+1
                      , sa+c
                      , sa+c+1);
                    pos++;
                    if (c > i)
                        step = 0;
                }
                else 
                    sa[c]--;
                sa[pos] = j;
            }

            if (!suffix)
            {
                auto i1 = (step == 0) ? i+1 : i;
                sa[i1] = EMPTY;
            }
        }

        // scan to shift-right the items in each bucket
        //   with its head being reused as a counter.
        if (!suffix)
            for (auto i = n-1; i > 0; i--)
            {
                auto j = static_cast<SignedIndex>(sa[i]);
                if (j < 0 && j != EMPTY) // is sa[i] a counter?
                {
                    std::move_backward(sa+i+j, sa+i, sa+i+1);
                    sa[i+j] = EMPTY;
                }
            }
    }

    template<class SEQ_ITR>
    void put_suffix1(
        const SEQ_ITR seq
      , SaItr sa
      , Index n1
    )
    {
        Index pos, cur, pre = -1;
        for (auto i = n1-1; i > 0; i--)
        {
            auto j = sa[i];
            sa[i] = EMPTY;
            cur = seq[j];
            if (cur != pre)
            {
                pre = cur;
                pos = cur;
            }
            sa[pos--] = j;
        }
    }
};
