# Fm-index
[![Build Status](https://travis-ci.org/8igMac/fm-index.svg?branch=master)](https://travis-ci.org/8igMac/fm-index)
[![codecov](https://codecov.io/gh/8igMac/fm-index/branch/master/graph/badge.svg)](https://codecov.io/gh/8igMac/fm-index)

Index for fast query substring in large text.

## Dependencies
- C++ 14
- cmake 3.12

## Build and run
- `mkdir -p build`
- `cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON`
- `cmake --build .`
- `ctest`

## Reference
- SACA-K: [Nong G. Practical linear-time O(1)-workspace suffix sorting for constant alphabets](https://dl.acm.org/citation.cfm?id=2493180)
- BWT-ISFM: [Elena Y. Practical Space-efficient Linear Time Construction of FM-index for Large Genomes](https://www.searchdl.org/Resources/Public/Conf/2018/BICOB/1034.pdf)
