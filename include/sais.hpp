#pragma once
#include <cstdlib>
#include <iterator>

template <class StringType, class SuffixArrayIteratorType>
class SAIS
{
    using IndexType  = typename std::iterator_traits<SuffixArrayIteratorType>::value_type;
    using BucketType = IndexType*;

    StringType&     s;
    SuffixArrayIteratorType SA;
    IndexType       n, K;
    int             level;
    BucketType      C, B;
    unsigned char*  t;

    const IndexType EMPTY = -1;
    unsigned char mask[8] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
#define tget(i) ( (t[(i)>>3] & mask[(i)&7]) ? 1 : 0 )
#define tset(i, b) t[(i)>>3] = (b) ? (mask[(i)&7] | t[(i)>>3]) : ((~mask[(i)&7]) & t[(i)>>3])
#define isLMS(i) ((i)>0 && tget(i) && !tget(i-1))

    void getCounts()
    {
        IndexType i;
        for (i = 0; i < K; ++i) { C[i] = 0; } // clear all buckets
        for (i = 0; i < n; ++i) { ++C[s[i]]; } // compute the size of each bucket
    }

    // compute the head or end of each bucket
    void getBuckets(bool end)
    {
        IndexType i, sum = 0;
        if(end) { for (i = 0; i < K; ++i) { sum += C[i]; B[i] = sum - 1; } }
        else { for (i = 0; i < K; ++i) { sum += C[i]; B[i] = sum - C[i]; } }
    }

    // compute SAl
    void induceSAl(const unsigned char* t, int level)
    {
        IndexType i, j;
        getBuckets(false); // find heads of buckets
        if (level == 0) 
            B[0]++;
        for (i = 0; i < n; i++)
            if (SA[i] != (IndexType)-1)
            {
                j = SA[i] - 1;
                if (j != (IndexType)-1 && !tget(j)) 
                    SA[B[s[j]]++] = j;
            }
    }

    // compute SAs
    void induceSAs(const unsigned char* t)
    {
        IndexType i, j;
        getBuckets(true); // find ends of buckets
        for (i = n-1; i != (IndexType)-1; i--)
            if (SA[i] != (IndexType)-1)
            {
                j = SA[i] - 1;
                if (j != (IndexType)-1 && tget(j)) 
                    SA[B[s[j]]--] = j;
            }
    }

    // find the suffix array SA of s[0..n-1] in {0..K}^n;
    // require s[n-1]=0 (the virtual sentinel!), n>=2;
    // use a space of at most 6.25n+(1) for a constant alphabet;
    // level starts from 0.
public:
    explicit SAIS(StringType& s, SuffixArrayIteratorType SA, IndexType n, IndexType K, int level = 0): s(s), SA(SA), n(n), K(K), level(level)
    {
        IndexType i, j;
        t = new unsigned char[n/8 + 1]; // LS-type array in bits

        // stage 1: reduce the problem by at least 1/2

        // Classify the type of each character
        tset(n-2, 0);
        tset(n-1, 1); // the sentinel must be in s1, important!!!
        for (i = n-3; i != (IndexType)-1; i--) 
            tset(i, (s[i] < s[i+1] || (s[i] == s[i+1] && tget(i+1) == 1)) ? 1 : 0);

        C = new IndexType[K];
        B = new IndexType[K];

        // sort all the S-substrings
        getCounts();
        getBuckets(true); // find ends of buckets
        for (i = 0; i < n; i++) SA[i] = EMPTY;
        for (i = n-2; i != (IndexType)-1; i--)
            if (isLMS(i)) 
                SA[B[s[i]]--] = i;
        SA[0] = n-1; // set the single sentinel LMS-substring

        induceSAl(t, level);
        induceSAs(t);

        delete[](C);
        delete[](B);

        // compact all the sorted substrings into the first n1 items of s
        // 2*n1 must be not larger than n (proveable)
        IndexType n1 = 0;
        for (i = 0; i < n; i++)
            if (isLMS(SA[i]))
                SA[n1++] = SA[i];

        // Init the name array buffer
        for (i = n1; i < n; i++) SA[i] = EMPTY;
        // find the lexicographic names of all substrings
        IndexType name = 0, prev = -1;
        for (i = 0; i < n1; i++)
        {
            IndexType pos = SA[i];
            bool diff = false;
            for (IndexType d = 0; d < n; d++)
                if (prev == -1 || pos+d == n-1 || prev+d == n-1 || s[pos+d] != s[prev+d] || tget(pos+d) != tget(prev+d))
                { diff = true; break; }
                else
                if (d > 0 && (isLMS(pos+d) || isLMS(prev+d)))
                    break;

            if (diff)
            { name++; prev = pos; }
            pos = pos / 2;
            SA[n1 + pos] = name - 1;
        }
        for (i = n-1, j = n-1; i >= n1; i--)
            if (SA[i] != (IndexType)-1) 
                SA[j--] = SA[i];

        // s1 is done now

        SuffixArrayIteratorType SA1 = SA, s1 = SA + n - n1;

        // stage 2: solve the reduced problem

        // recurse if names are not yet unique
        if (name < n1)
        {
            SAIS<SuffixArrayIteratorType, SuffixArrayIteratorType>(s1, SA1, n1, name, level+1);
        }
        else
        {   // generate the suffix array of s1 directly
            for (i = 0; i < n1; i++) SA1[s1[i]] = i;
        }

        // stage 3: induce the result for the original problem

        C = new IndexType[K];
        B = new IndexType[K];

        // put all left-most S characters into their buckets
        getCounts();
        getBuckets(true); // find ends of buckets

        j = 0;
        for (i = 1; i < n; i++)
            if (isLMS(i)) 
                s1[j++] = i; // get p1
        for (i = 0; i < n1; i++) SA1[i] = s1[SA1[i]]; // get index in s1
        for (i = n1; i < n; i++) SA[i] = EMPTY; // init SA[n1..n-1]
        for (i = n1-1; i != (IndexType)-1; i--)
        {
            j = SA[i]; SA[i] = EMPTY;
            if (level == 0 && i == 0)
                SA[0] = n-1;
            else
                SA[B[s[j]]--] = j;
        }

        induceSAl(t, level);
        induceSAs(t);

        delete[](t);
        delete[](C);
        delete[](B);
    }
};

