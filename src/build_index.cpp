#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include "fm_index.hpp"
#include "saca_k.hpp"

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " FILE\n";
        return 1;
    }
    std::ifstream ifs(argv[1]);

    // Check file size
    int file_size;
    try {
        ifs.seekg(0, std::ios_base::end);
        file_size = ifs.tellg();
        std::cerr << "file size: " << file_size << ", ";
        ifs.seekg(0); // rewind
    } catch (const std::ios_base::failure& e)
    {
        std::cerr << "Can't seekg: " << e.what() 
                  << ", error code: " << e.code() << "\n";
    }

    // Read genome
    std::default_random_engine eng;
    std::uniform_int_distribution<int> dist(0, 3); 
    std::vector<char> seq;
    std::vector<char> char_set {'A', 'C', 'G', 'T'};
    seq.reserve(file_size);
    std::string buf;
    auto start = std::chrono::high_resolution_clock::now();
    while (std::getline(ifs, buf))
    {
        for (auto& chr : buf)
        {
            switch (chr)
            {
                case 'A': case 'a': seq.push_back('A'); break;
                case 'C': case 'c': seq.push_back('C'); break;
                case 'G': case 'g': seq.push_back('G'); break;
                case 'T': case 't': seq.push_back('T'); break;
                default: seq.push_back(char_set[dist(eng)]);
            }
        }
    }
    seq.push_back('A'); // $
    ifs.close();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "seq size: " << seq.size() << ", "
              << "file read time: " << elapsed.count() << "s\n";

    // debug
    for (auto& elem : seq)
        std::cerr << elem;
    std::cerr << std::endl;

    // construct fm-index
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
    start = std::chrono::high_resolution_clock::now();
    FmIndex<decltype(seq), uint32_t, 2, SACA_K> index(seq, map, 16);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cerr << "FmIndex construction time: " 
              << elapsed.count() << "s\n";

    return 0;
}
