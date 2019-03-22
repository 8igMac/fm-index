#pragma once

template<class DERIVED>
struct FmIndex
{
    /// @brief Build index, call derived class implementation
    template<class... ARGS>
    void build(ARGS&&... args)
    {
        static_cast<DERIVED*>(this)->build(
            std::forward<ARGS>(args)...
        );
    }
    
    /// @brief Load index from file
    /// @param filename Input file
    void load(const std::string& filename)
    { 
        static_cast<DERIVED*>(this)->load(filename); 
    }
    
    /// @brief Save index to file
    /// @param filename Output file
    void save(const std::string& filename)
    { 
        static_cast<DERIVED*>(this)->save(filename); 
    }
    
    /// @brief Map element in last column to first column 
    /// @param i Last column index
    /// @param c Character
    /// @return First column index
    template<class INDEX, class CHAR>
    decltype(auto) lf_mapping(INDEX i, CHAR c)
    { 
        return static_cast<DERIVED*>(this)->get_c(c) +  
               static_cast<DERIVED*>(this)->get_occ(i, c); 
    }
    
    /// @brief Map the i-th elemnet in bwt to the original seq
    /// @param i i-th element in bwt
    /// @return Location in the original seq
    template<class INDEX>
    decltype(auto) get_location(INDEX i)
    { 
        return static_cast<DERIVED*>(this)->get_location(i); 
    }
};
