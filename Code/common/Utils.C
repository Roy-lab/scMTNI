#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "Utils.H"

using namespace std;

// Split on whitespace
vector<string>
Utils::split(string const s)
{
    stringstream ss(s);
    string token;
    vector<string> tokens;

    while(ss.good())
    {
        ss >> token;
        tokens.push_back(move(token));
    }
    return tokens;
}

// Split on separator
vector<string>
Utils::split(string const s, char const sep)
{
    stringstream ss(s);
    string token;
    vector<string> tokens;
    
    while (getline(ss, token, sep))
    {
        tokens.push_back(move(token));
    }
    return tokens;
}

// Split on line length, e.g., for FASTA output
vector<string> 
Utils::splitLines(string const s, int const len)
{
    vector<string> v;
    int i = 0;
    while(i < s.length())
    {
        int l = min(len, static_cast<int>(s.length()) - i);
        v.push_back(s.substr(i, l));
        i += l;
    }
    return v;
}

string
Utils::toupper(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { 
        return std::toupper(c); 
    });
    return s;
}

string
Utils::tolower(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { 
        return std::tolower(c); 
    });
    return s;
}

bool 
Utils::hasMatchingPrefix(string s, vector<string> prefixes, bool caseSensitive)
{
    if (!caseSensitive)
    {
        s = Utils::tolower(s);
    }
    for (auto p : prefixes)
    {
        if (!caseSensitive)
        {
            p = Utils::tolower(p);
        }
        if (s.rfind(p, 0) == 0) 
        {
            return true;
        }
    }
    return false;
}

string 
Utils::extension(string const filename) 
{
    auto pos = filename.find_last_of('.');
    if (pos == string::npos) return "";
    return filename.substr(pos + 1);
}

// Analogous to Python maketrans
function<string(string)>
Utils::maketrans(const string& from, const string& to)
{
    // Set up the translation map
    unordered_map<char, char> map;
    for (string::size_type i = 0; i != min(from.size(), to.size()); ++i)
    {
        map[from[i]] = to[i];
    }
    // Return the function that will do the translation
    return [=](string s) 
    {
        for (auto& c : s)
        {
            const auto mapped_c = map.find(c);
            if (mapped_c != map.end()) 
            {
                c = mapped_c -> second;
            }
        }
        return s;
    };
}

// Return the reverse complement of a sequence
string 
Utils::revcomp(string sequence)
{
    const auto translate = maketrans("AGCTagctWSKMYRBDHVnN", "TCGAtcgaWSMKRYVHDBnN");
    auto translated = translate(sequence);
    reverse(translated.begin(), translated.end());
    return translated;
}

// Convert a string representation of a double to double. Handles out of range values.
double 
Utils::convertDouble(string val)
{
    // The following are not as small or large as the Standard Library versions, but they are displayable in Excel
    const double myDBL_MIN = 1e-300;
    const double myDBL_MAX = 1e+300;

    try
    {
        return stod(val);
    }
    catch (out_of_range& oor)
    {
        auto n = val.find("e-");
        if (n != string::npos) return myDBL_MIN;
        n = val.find("e+");
        if (n != string::npos) return myDBL_MAX;
        throw out_of_range(val);
    }
}

// Convert A, C, G, T, N to the corresponding index.
// Index is 0-based.
int 
Utils::convertBaseToIndex(const char base)
{
    if (base == 'A' || base == 'a') return 0;
    if (base == 'C' || base == 'c') return 1;
    if (base == 'G' || base == 'g') return 2;
    if (base == 'T' || base == 't') return 3;
    return 4;
}

// Convert index to A, C, G, T, N. Index is 0-based.
char 
Utils::convertIndexToBase(const int index)
{
    switch(index)
    {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
    }
    return 'N';
}

// Returns a vector of integers 1 - N in random order
vector<int> 
Utils::randomSequence(int n)
{
    vector<int> seq;
    for (int i = 0; i < n; ++i)
    {
        seq.push_back(i);
    }

    random_device rd;
    mt19937 generator(rd());
 
    shuffle(seq.begin(), seq.end(), generator);

    return seq;
}

bool 
Utils::contains(shared_ptr<vector<string>> container, string value)
{
    auto iter = find(container->begin(), container->end(), value);
    return (iter != container->end());
}

bool 
Utils::contains(shared_ptr<unordered_map<string, double>> container, string value)
{
    auto iter = container->find(value);
    return (iter != container->end());
}

void 
Utils::dumpVector(vector<string> vec)
{
    for (string s : vec)
    {
        cout << s << endl;
    }
}
