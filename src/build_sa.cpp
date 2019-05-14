#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include "sais.hpp"
using Sorter = SAIS;

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
    seq.reserve(file_size);
    std::string buf;
    auto start = std::chrono::high_resolution_clock::now();
    while (std::getline(ifs, buf))
    {
        for (auto& chr : buf)
        {
            switch (chr)
            {
                case 'A': case 'a': seq.push_back(0); break;
                case 'C': case 'c': seq.push_back(1); break;
                case 'G': case 'g': seq.push_back(2); break;
                case 'T': case 't': seq.push_back(3); break;
                default: seq.push_back(dist(eng));
            }
        }
    }
    seq.push_back(0); // $
    ifs.close();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "seq size: " << seq.size() << ", "
              << "file read time: " << elapsed.count() << "s\n";

    // Construct suffix array 
    std::vector<std::size_t> sa(seq.size());
    Sorter sorter;
    start = std::chrono::high_resolution_clock::now();
    sorter(seq, sa, 4); 
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cerr << "Suffix array construction time: " 
              << elapsed.count() << "s\n";

    return 0;
}
