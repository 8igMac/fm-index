# Fm-index
Index for fast query substring in large text.

## Dependencies
- gcc (c++11)
- cmake
- gtest
- Boost

## How to build

## How to run

## Todo list
- Basic interface
    - constructor 
        - sequence constructor
        - copy constructor
        - file constructor (input archive)
    - output archive
    - exact-match
- Architecture
- Test data (use file to scale)
- Implemnetation
    - fm-index from SAIS
    - blockwise-bwt
    - BWT-IS

## Reference
- blockwise-bwt: [Kärkkäinen J., Fast BWT in small space by blockwise suffix sorting, 2007](https://www.sciencedirect.com/science/article/pii/S0304397507005245)
- SAIS: [Nong G., et al. Two efficient algorithms for linear time suffix array construction, 2011](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5582081&tag=1)
- BWT-IS:[D. Okanohara and K. Sadakane, A linear-time Burrows-Wheeler transform using induced sorting, 2009](https://link.springer.com/chapter/10.1007/978-3-642-03784-9_9)
