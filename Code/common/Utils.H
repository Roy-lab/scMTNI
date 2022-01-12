#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/* Assorted utility routines. */

class Utils
{
public:
    // Split on whitespace
    static 
    std::vector<std::string> split(std::string const);
    
    // Split on separator
    static 
    std::vector<std::string> split(std::string const, char const);

    // Split on line length, e.g., for FASTA output
    static
    std::vector<std::string> splitLines(std::string const, int const);

    // Trim whitespace

    static 
    inline std::string& ltrim(std::string& s) 
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    static 
    inline std::string& rtrim(std::string& s) 
    {
        s.erase(std::find_if(s.rbegin(), s.rend(),
                std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    static 
    std::string toupper(std::string s);

    static 
    std::string tolower(std::string s);

    static 
    inline std::string& trim(std::string& s) 
    {
        return ltrim(rtrim(s));
    }

    static
    bool hasMatchingPrefix(std::string s, std::vector<std::string> prefixes, bool caseSensitive=false);

    static
    std::string extension(const std::string filename);

    // Analogous to Python maketrans
    static 
    std::function<std::string(std::string)>
    maketrans(const std::string& from, const std::string& to);

    // Return the reverse complement of a sequence
    static 
    std::string revcomp(const std::string sequence);

    // Convert a string representation of a double to double. Handles out of range values.
    static
    double convertDouble(std::string val);

    // Convert A, C, G, T, N to the corresponding index
    static
    int convertBaseToIndex(const char base);

    static
    char convertIndexToBase(const int index);

    // Returns a vector of integers 0 - N-1 in random order
    static
    std::vector<int> randomSequence(int n);

    static
    bool contains(std::shared_ptr<std::vector<std::string>> container, std::string value);

    static
    bool contains(std::shared_ptr<std::unordered_map<std::string, double>> container, std::string value);

    static 
    void dumpVector(const std::vector<std::string> vec);
};

#endif /* UTILS_H */

