#pragma once
#include <cmath>
#include <list>
#include <fstream>
#include <algorithm>
#include <functional>

template<
    typename SEQ
  , typename INDEX 
  , int BITS
  , template<typename, typename> typename SORTER
>
class FmIndex
{
    using CharType     = typename SEQ::value_type;
    using CTableType   = std::array<INDEX
                          , static_cast<int>(std::pow(2, BITS))>;
    using OccTableType = std::vector<CTableType>;
    using LocTableType = std::vector<
                            std::pair<INDEX, INDEX>>;
    
    /// @brief bwt of the orignal seq
    SEQ               bwt_;

    /// @brief Sampled mapping of index from bwt to original seq
    LocTableType      loc_table_;

    /// @brief Sampled occurence of each alphabet in bwt
    OccTableType      occ_table_ {{}}; 

    /// @brief Start position of each alphabet in first column
    CTableType        c_table_ {};

    /// @brief The index of $(sentinal) in bwt
    INDEX             primary_index_;

    /// @brief Sample rate, valid value are 2^n, n>=0 
    /// (representing 1, 1/2, 1/4, ...)
    INDEX             sample_rate_;

    /// @brief Bit vector, set to 1 if i-th bwt's suffix array 
    ///        location is stored.
    std::vector<bool> bwt_marked_;

    std::function<INDEX(CharType)> map_;

  public:
    /// @brief Build fm-index using bwt-isfm alogrithm
    /// @param seq Sequence, required $(smalest alphabet) be 
    ///        inserted at the end
    /// @param map Map alphabet to their rank
    /// @param step Sample rate, valid value are 2^n, n>=0
    template<class MAPPER>
    FmIndex (const SEQ& seq, MAPPER map, INDEX step = 1)
             : map_(map)
             , sample_rate_(step)
    {
        assert(!(sample_rate_ & (sample_rate_-1)));

        // Init member var and other param
        constexpr int alph_size = std::pow(2, BITS);
        constexpr int bit_mask  = alph_size - 1;
        constexpr int short_lms_len = 12;

        // Scan through seq (right-to-left) to identify L/S type 
        // (S-type set to true)
        std::vector<bool> type(seq.size());
        type[seq.size() - 1] = true; // $ is S-type
        type[seq.size() - 2] = false; 
        for (auto i = seq.size() - 3; ~i; i--)
            if (seq[i] < seq[i+1] || 
                (seq[i] == seq[i+1] && type[i+1]))
                type[i] = true;

        // Count the total number of each alphabet
        // $ is counted as the smallest alphabet
        for (auto i = 0; i < seq.size(); i++)
            c_table_[map_(seq[i])]++;
        // Calculate accumulative sum
        INDEX sum = 0;
        for (auto& i : c_table_)
        {
            std::swap(sum, i);
            sum += i;
        }
        
        // Calculate number of LMS
        INDEX lms_size = 0;
        for (auto i = 0; i < type.size(); i++)
            if (is_lms(i, type))
                lms_size++;

        ///////////////////////
        // Produce shorten seq
        ///////////////////////
        constexpr int hash_size = std::pow(alph_size, short_lms_len);
        std::vector<INDEX> hash_table(hash_size); // for short LMS 

        // Scan through seq (right-to-left) to extract LMS,
        // only store distinct LMS (all long LMS and distinct short 
        // LMS)
        std::vector<INDEX> lms(lms_size);
        uint32_t key = 0;
        INDEX lms_len = 0;
        INDEX distinct_lms_size = 0;
        for (auto i = seq.size() - 1; ~i; i--)
        {
            // record LMS substr's backward complement
            // ex: ATGC -> ~(CGTA) = GCAT
            // TODO: this can be speed up if we can extract the
            // compressed LMS substr at once
            uint32_t complement = ~map_(seq[i]) && bit_mask;
            key = (key << BITS) + complement;
            lms_len++;

            if (is_lms(i, type))
            {
                if (lms_len <= short_lms_len &&
                    hash_table[key] == 0)
                {
                    hash_table[key] = i;
                    lms[distinct_lms_size++] = i;
                }
                else if (lms_len > short_lms_len)
                    lms[distinct_lms_size++] = i;

                key = complement;
                lms_len = 1;
            }
        }
        std::reverse(lms.begin(), lms.begin() + distinct_lms_size);

        // INFO
        // for (const auto& i : seq)
        //     std::cerr << i;
        // std::cerr << std::endl;
        //
        // for (auto i = 0; i < type.size(); i++)
        //     std::cerr << type[i];
        // std::cerr << std::endl;
        //
        // std::cerr << "distinct lms location: ";
        // for (auto i = 0; i < distinct_lms_size; i++)
        //     std::cerr << (int)lms[i] << " ";
        // std::cerr << std::endl;
        //
        // std::cerr << "all lms location: ";
        // for (auto i = 0; i < type.size(); i++)
        //     if (is_lms(i, type))
        //         std::cerr << i << " ";
        // std::cerr << std::endl;
        //
        // std::cerr << "seq size: " << type.size()
        //           << ", LMS size: " << lms_size
        //           << ", Distinct LMS size: " << (int)distinct_lms_size
        //           << ", Remains: " 
        //           << (float)distinct_lms_size / (float)lms_size
        //           << std::endl;

        // Qsort distinct_LMS
        std::vector<INDEX> lms_sa(lms_size);
        for (auto i = 0; i < distinct_lms_size; i++)
            lms_sa[i] = i;

        // // debug: 1, 17, 30, 46, 57, 60, 72
        // std::cerr << "distinct lms index: ";
        // for (auto i = 0; i < distinct_lms_size; i++)
        //     std::cerr << (int)lms[i] << " ";
        // std::cerr << std::endl;

        std::sort(lms_sa.begin(), lms_sa.begin() + distinct_lms_size, 
            // TODO: this can be speedup, for example extract the 
            // compressed LMS substr at once as a integer and compare
            // them
            [&seq, &type, &lms, this](auto a, auto b)
            {
                auto a_pos = lms[a], b_pos = lms[b];

                // handle $
                if (a_pos == seq.size()-1)
                    return true;
                else if (b_pos == seq.size()-1)
                    return false;

                // compare first LMS char
                if (map_(seq[a_pos]) < map_(seq[b_pos]))
                    return true;
                else if (map_(seq[a_pos]) > map_(seq[b_pos]))
                    return false;

                a_pos++; b_pos++;
                while (a_pos < seq.size()-1 && b_pos < seq.size()-1)
                {
                    // compare char
                    if (map_(seq[a_pos]) < map_(seq[b_pos]))
                        return true;
                    else if (map_(seq[a_pos]) > map_(seq[b_pos]))
                        return false;

                    // compare type: if char is same, L < S
                    if (type[a_pos] != type[b_pos])
                       return !type[a_pos]; 

                    // if a reach LMS, b must reach LMS as well
                    // then they are equal
                    if (is_lms(a_pos, type))
                        return true;

                    a_pos++; b_pos++;
                }

                // either one of a or b reach $
                if (a_pos == seq.size() - 1)
                    return true;
                else 
                    return false;
            });

        // // debug: 6, 5, 2, 0, 4, 3, 1
        // std::cerr << "distinct lms sa: ";
        // for (auto i = 0; i < distinct_lms_size; i++)
        //     std::cerr << (int)lms_sa[i] << " ";
        // std::cerr << std::endl;

        // Assign name to sorted distinct LMS
        // (TODO: use more space efficient index like char)
        std::vector<INDEX> lms_name(distinct_lms_size);
        INDEX name = 0;
        auto sa_equal = 
            [&seq, &type, &lms, this](auto a, auto b)
            {
                auto a_pos = lms[a], b_pos = lms[b];

                // handle $
                if (a_pos == seq.size()-1 || b_pos == seq.size()-1)
                    return false;

                // compare first LMS char
                if (map_(seq[a_pos]) != map_(seq[b_pos]))
                    return false;

                a_pos++; b_pos++;
                while (a_pos < seq.size()-1 && b_pos < seq.size()-1)
                {
                    // compare char
                    if (map_(seq[a_pos]) != map_(seq[b_pos]) &&
                        type[a_pos] != type[b_pos])
                        return false;

                    // if a reach LMS, b must reach LMS as well
                    // then they are equal
                    if (is_lms(a_pos, type))
                        return true;

                    a_pos++; b_pos++;
                }
                   
                // either one of a or b reach $
                return false;
            };
        for (auto i = 0; i < distinct_lms_size; i++)
        {
            if (i != 0 && !sa_equal(lms_sa[i], lms_sa[i-1]))
                name++;
            lms_name[i] = name;
        }
        // // debug: 0, 1, 2, 2, 3, 4, 5
        // std::cerr << "assigned name: ";
        // for (auto i = 0; i < distinct_lms_size; i++)
        //     std::cerr << (int)lms_name[i] << " ";
        // std::cerr << std::endl;


        // Place named lms back to order 
        // TODO: inplace reorder 
        std::vector<INDEX> correct_order(distinct_lms_size);
        for (auto i = 0; i < distinct_lms_size; i++)
            correct_order[lms_sa[i]] = lms_name[i];
        // // debug: 2, 5, 2, 4, 3, 1, 0
        // std::cerr << "correct order: ";
        // for (auto i = 0; i < distinct_lms_size; i++)
        //     std::cerr << (int)correct_order[i] << " ";
        // std::cerr << std::endl;

        // Produce T1
        key = 0;
        lms_len = 0;
        for (auto i = seq.size() - 1
                , j = distinct_lms_size-1
                , k = lms_size-1; ~i; i--)
        {
            uint32_t complement = ~map_(seq[i]) && bit_mask;
            key = (key << BITS) + complement;
            lms_len++;

            if (is_lms(i, type))
            {
                if (i == lms[j]) // distinct LMS
                {
                    if (lms_len <= short_lms_len)
                        hash_table[key] = correct_order[j];
                        
                    lms[k] = correct_order[j];
                    j--;
                }
                else // short LMS that is not recorded
                    lms[k] = hash_table[key];

                k--;
                // reset hash key
                key = complement;
                lms_len = 1;
            }
        }
        // // debug: 2, 3, 5, 4, 2, 3, 4, 3, 1, 0
        // std::cerr << "T1: ";
        // for (auto i = 0; i < lms_size; i++)
        //     std::cerr << (int)lms[i] << " ";
        // std::cerr << std::endl;

        /////////////////////////////////////////
        // Produce LMS SA if name not yet unique
        /////////////////////////////////////////
        if (name < lms_size)
        {
            SORTER<decltype(lms), decltype(lms_sa)> sa_builder;
            sa_builder.build(lms, lms_sa, name+1);
        }
        // // debug: 9, 8, 4, 0, 7, 5, 1, 3, 6, 2
        // std::cerr << "lms sa(before): ";
        // for (auto i = 0; i < lms_size; i++)
        //     std::cerr << (int)lms_sa[i] << " ";
        // std::cerr << std::endl;

        //////////////////////////////////
        // Transform SA1 to T's position
        //////////////////////////////////
        // Get all LMS
        for (auto i = 0, j = 0; i < seq.size(); i++)
            if (is_lms(i, type))
                lms[j++] = i; 
        // Transform SA1 to T's position
        for (auto i = 0; i < lms_size; i++)
            lms_sa[i] = lms[lms_sa[i]];
        // // debug: 72, 60, 30, 1, 57, 43, 14, 19, 46, 17
        // std::cerr << "lms sa(after): ";
        // for (auto i = 0; i < lms_size; i++)
        //     std::cerr << (int)lms_sa[i] << " ";
        // std::cerr << std::endl;


        ///////////////
        // Induce sort
        ///////////////
        std::vector<bool> pos_bit(seq.size());
        {
            std::vector<std::list<INDEX>> LMS(alph_size);
            std::vector<std::list<INDEX>>   L(alph_size);
            std::vector<std::list<INDEX>>  LS(alph_size);
            std::vector<std::list<INDEX>>   S(alph_size);
            CTableType head, tail;
            bwt_.resize(seq.size());
            bwt_marked_.resize(seq.size());

            // init head, tail
            for (auto i = 0; i < c_table_.size(); i++)
            {
                if (i == 0) // skip 0, it is for $
                    head[i] = 1;
                else
                    head[i] = c_table_[i];

                if (i == c_table_.size() - 1)
                    tail[i] = seq.size() - 1;
                else 
                    tail[i] = c_table_[i+1] - 1;
            }

            // Put sorted LMS to correspond character bucket
            for (auto i = 0; i < lms_sa.size(); i++)
                LMS[map_(seq[lms_sa[i]])].push_back(lms_sa[i]);

            std::ofstream ofs("temp_file");
            // handle $ first
            {
                auto idx = LMS[0].front();
                LMS[0].pop_front();
                induce_l(idx, seq, L, LS, head, pos_bit, ofs, false);
            }

            // Left-to-right scan
            for (auto i = 0; i < alph_size; i++)
            {
                while (!L[i].empty())
                {
                    auto idx = L[i].front();
                    L[i].pop_front();
                    induce_l(idx, seq, L, LS, head, pos_bit, ofs, true);
                }
                while (!LMS[i].empty())
                {
                    auto idx = LMS[i].front();
                    LMS[i].pop_front();
                    induce_l(idx, seq, L, LS, head, pos_bit, ofs, false);
                }
            }

            // Right-to-left scan
            for (auto i = alph_size-1; ~i; i--)
            {
                while (!S[i].empty())
                {
                    auto idx = S[i].back();
                    S[i].pop_back();
                    induce_s(idx, seq, S, tail, pos_bit, ofs);
                }
                while (!LS[i].empty()) 
                {
                    auto idx = LS[i].back();
                    LS[i].pop_back();
                    induce_s(idx, seq, S, tail, pos_bit, ofs);
                }
            }
            ofs.close();
        }

        // // Get location table from file
        // std::ifstream ifs("temp_file");
        // while (!ifs.eof())
        // {
        //     INDEX bwt_pos, sa_pos;
        //     ifs >> bwt_pos >> sa_pos;
        //     std::cerr << (int)bwt_pos << " " << (int)sa_pos << std::endl;
        //     loc_table[bwt_pos / sample_rate_] = std::make_pair(bwt_pos, sa_pos);
        // }
        // ifs.close();

        // sort location_table
        std::sort(loc_table_.begin(), loc_table_.end(), 
            [](const auto& lhs, const auto& rhs)
            { return lhs.first < rhs.first; });

        // Calculate occ
        for (auto& i : c_table_)
            i = 0;
        for (auto i = 0; i < bwt_.size(); i++)
        {
            if (i != primary_index_)
                c_table_[map_(bwt_[i])]++;
            if ((i+1 & (sample_rate_ - 1)) == 0)
                occ_table_.emplace_back(c_table_);
        }

        // Caculate c_table
        sum = 1;
        for (auto& i : c_table_)
        {
            std::swap(sum, i);
            sum += i;
        } 
    }

    /// @brief Map the i-th elemnet in bwt to the original seq
    /// @param i i-th element in bwt
    /// @return Location in the original seq
    INDEX get_location(INDEX i) const
    { 
        INDEX step_count;
        for (step_count = 0; !bwt_marked_[i]; step_count++)
            i = lf_mapping(i, bwt_[i]);

        auto itr = std::lower_bound(
            loc_table_.begin(), loc_table_.end(), i,
            [](const std::pair<INDEX, INDEX>& lhs, INDEX rhs)
            { return lhs.first < rhs; });

        return itr->second + step_count;
    }

    /// @brief Map element in last column to first column 
    /// @param i Last column index
    /// @param c Character
    /// @return First column index
    INDEX lf_mapping(INDEX i, CharType c) const
    { 
        return c_table_[map_(c)] + get_occ(i, c);
    }
    
  private:
    bool is_lms(INDEX i, const std::vector<bool>& type) const
    {
        // type[0]'s previous is $, so it can not be LMS
        if (i != 0 && type[i] && !type[i-1])
            return true;
        else
            return false;
    }

    void induce_l(
        INDEX idx
      , const SEQ& seq
      , std::vector<std::list<INDEX>>& L
      , std::vector<std::list<INDEX>>& LS
      , CTableType& head
      , std::vector<bool>& pos_bit
      , std::ofstream& ofs
      , bool is_l_queue
    )
    {
        INDEX maxui = -1;
        auto idx_prev = (idx == 0) ? seq.size()-1 : idx-1;
        auto idx_pprev = (idx_prev == 0) ? seq.size()-1 : idx_prev-1;
        auto c = map_(seq[idx]);
        auto c_prev = map_(seq[idx_prev]);

        // record $ position in bwt as primary index
        if (idx_pprev == seq.size()-1)
            primary_index_ = head[c_prev];

        if (c <= c_prev)
        {
            L[c_prev].push_back(idx_prev);
            bwt_[head[c_prev]] = seq[idx_pprev];
            // sample if mod step is 0
            if ( (idx_prev & (sample_rate_ - 1)) == 0 ) 
            {
                bwt_marked_[head[c_prev]] = true;

                // ofs << head[c_prev] << " " << idx_prev
                //     << "\n";
                loc_table_.emplace_back(
                    std::make_pair(head[c_prev], idx_prev));
            }
            if (idx_prev > maxui)
                pos_bit[head[c_prev]] = true;
            head[c_prev]++;
        }
        else if (is_l_queue)
            LS[c].push_back(idx);
    }

    void induce_s(
        INDEX idx
      , const SEQ& seq
      , std::vector<std::list<INDEX>>& S
      , CTableType& tail
      , std::vector<bool>& pos_bit
      , std::ofstream& ofs
    )
    {
        // skip $ because $ is LMS
        if (idx == seq.size()-1)
            return;
    
        INDEX maxui = -1;
        auto idx_prev = (idx == 0) ? seq.size()-1 : idx-1;
        auto idx_pprev = (idx_prev == 0) ? seq.size()-1 : idx_prev-1;
        auto c = map_(seq[idx]);
        auto c_prev = map_(seq[idx_prev]);

        // $'s bwt_pos is 0
        auto bwt_pos = (idx_prev == seq.size()-1) ? 0 : tail[c_prev];
        if (c_prev <= c)
        {
            S[c_prev].push_front(idx_prev);
            bwt_[bwt_pos] = seq[idx_pprev];
            // sample if mod step is 0
            if ( (idx_prev & (sample_rate_ - 1)) == 0 ) 
            {
                bwt_marked_[bwt_pos] = true;

                // ofs << bwt_pos << " " << idx_prev 
                //     << "\n";
                loc_table_.emplace_back(
                    std::make_pair(bwt_pos, idx_prev));
            }
            if (idx_prev > maxui)
                pos_bit[bwt_pos] = true;

            if (bwt_pos != 0) // for $
                tail[c_prev]--;
        }
    }

    /// @brief Get occurence of character c uptile bwt index i 
    /// @param i Bwt index
    /// @param c Character
    /// @return Number of occurence
    INDEX get_occ(INDEX i, CharType c) const
    { 
        auto occ_lower_index = i / sample_rate_;
        auto occ_upper_index = occ_lower_index + 1;
        auto lower_offset = i & sample_rate_-1;
        INDEX c_count = 0;

        if (lower_offset <= sample_rate_ / 2 || 
            occ_upper_index == occ_table_.size())
        {
            auto lower_index = occ_lower_index * sample_rate_;
            for (auto j = lower_index; j < i; j++)
                if (bwt_[j] == c && j != primary_index_)
                    c_count++;

            return occ_table_[occ_lower_index][map_(c)] + c_count;
        }
        else
        {
            auto upper_index = occ_upper_index * sample_rate_;
            for (auto j = i; j < upper_index; j++)
                if (bwt_[j] == c && j != primary_index_)
                    c_count++;

            return occ_table_[occ_upper_index][map_(c)] - c_count;
        }
    }
};
