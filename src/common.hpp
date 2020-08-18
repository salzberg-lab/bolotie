#pragma once

#include <cassert>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

static constexpr char    table_nt[4]   = { 'A', 'C', 'G', 'T' };
static constexpr uint8_t nt_table[256] =
    { 99, 99, 99, 99, 99, 99, 99, 99,   // 0
      99, 99, 99, 99, 99, 99, 99, 99,   // 8
      99, 99, 99, 99, 99, 99, 99, 99,   // 16
      99, 99, 99, 99, 99, 99, 99, 99,   // 24
      99, 99, 99, 99, 99, 99, 99, 99,   // 32
      99, 99, 99, 99, 99, 99, 99, 99,   // 40
      99, 99, 99, 99, 99, 99, 99, 99,   // 48
      99, 99, 99, 99, 99, 99, 99, 99,   // 56
      99,  0, 99,  1, 99, 99, 99, 2,    // 64
      99, 99, 99, 99, 99, 99, 99, 99,   // 72
      99, 99, 99, 99,  3, 99, 99, 99,   // 80
      99, 99, 99, 99, 99, 99, 99, 99,   // 88
      99,  0, 99,  1, 99, 99, 99, 2,    // 96
      99, 99, 99, 99, 99, 99, 99, 99,   // 104
      99, 99, 99, 99,  3, 99, 99, 99,   // 112
      99, 99, 99, 99, 99, 99, 99, 99,   // 120
      99, 99, 99, 99, 99, 99, 99, 99,   // 128
      99, 99, 99, 99, 99, 99, 99, 99,   // 136
      99, 99, 99, 99, 99, 99, 99, 99,   // 144
      99, 99, 99, 99, 99, 99, 99, 99,   // 152
      99, 99, 99, 99, 99, 99, 99, 99,   // 160
      99, 99, 99, 99, 99, 99, 99, 99,   // 168
      99, 99, 99, 99, 99, 99, 99, 99,   // 176
      99, 99, 99, 99, 99, 99, 99, 99,   // 184
      99, 99, 99, 99, 99, 99, 99, 99,   // 192
      99, 99, 99, 99, 99, 99, 99, 99,   // 200
      99, 99, 99, 99, 99, 99, 99, 99,   // 208
      99, 99, 99, 99, 99, 99, 99, 99,   // 216
      99, 99, 99, 99, 99, 99, 99, 99,   // 224
      99, 99, 99, 99, 99, 99, 99, 99,   // 232
      99, 99, 99, 99, 99, 99, 99, 99,   // 240
      99, 99, 99, 99, 99, 99, 99, 99 }; // 248

static constexpr bool IUPAC[256] =
        { 0, 0, 0, 0, 0, 0, 0, 0,   // 0
          0, 0, 0, 0, 0, 0, 0, 0,   // 8
          0, 0, 0, 0, 0, 0, 0, 0,   // 16
          0, 0, 0, 0, 0, 0, 0, 0,   // 24
          0, 0, 0, 0, 0, 0, 0, 0,   // 32
          0, 0, 0, 0, 0, 0, 0, 0,   // 40
          0, 0, 0, 0, 0, 0, 0, 0,   // 48
          0, 0, 0, 0, 0, 0, 0, 0,   // 56
          0, 1, 1, 1, 1, 0, 0, 1,   // 64 ABCDG
          1, 0, 0, 1, 0, 1, 1, 0,   // 72 HKMN
          0, 0, 1, 1, 1, 1, 1, 1,   // 80 RSTUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 88 Y
          0, 1, 1, 1, 1, 0, 0, 1,   // 96 ABCDG
          1, 0, 0, 1, 0, 1, 1, 0,   // 104 HKMN
          0, 0, 1, 1, 1, 1, 1, 1,   // 112 RSTUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 120 Y
          0, 0, 0, 0, 0, 0, 0, 0,   // 128
          0, 0, 0, 0, 0, 0, 0, 0,   // 136
          0, 0, 0, 0, 0, 0, 0, 0,   // 144
          0, 0, 0, 0, 0, 0, 0, 0,   // 152
          0, 0, 0, 0, 0, 0, 0, 0,   // 160
          0, 0, 0, 0, 0, 0, 0, 0,   // 168
          0, 0, 0, 0, 0, 0, 0, 0,   // 176
          0, 0, 0, 0, 0, 0, 0, 0,   // 184
          0, 0, 0, 0, 0, 0, 0, 0,   // 192
          0, 0, 0, 0, 0, 0, 0, 0,   // 200
          0, 0, 0, 0, 0, 0, 0, 0,   // 208
          0, 0, 0, 0, 0, 0, 0, 0,   // 216
          0, 0, 0, 0, 0, 0, 0, 0,   // 224
          0, 0, 0, 0, 0, 0, 0, 0,   // 232
          0, 0, 0, 0, 0, 0, 0, 0,   // 240
          0, 0, 0, 0, 0, 0, 0, 0 }; // 248

static constexpr bool IUPAC_noACGT[256] =
        { 0, 0, 0, 0, 0, 0, 0, 0,   // 0
          0, 0, 0, 0, 0, 0, 0, 0,   // 8
          0, 0, 0, 0, 0, 0, 0, 0,   // 16
          0, 0, 0, 0, 0, 0, 0, 0,   // 24
          0, 0, 0, 0, 0, 0, 0, 0,   // 32
          0, 0, 0, 0, 0, 0, 0, 0,   // 40
          0, 0, 0, 0, 0, 0, 0, 0,   // 48
          0, 0, 0, 0, 0, 0, 0, 0,   // 56
          0, 0, 1, 0, 1, 0, 0, 0,   // 64 BD
          1, 0, 0, 1, 0, 1, 1, 0,   // 72 HKMN
          0, 0, 1, 1, 1, 1, 1, 1,   // 80 RSUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 88 Y
          0, 1, 1, 1, 1, 0, 0, 1,   // 96 BD
          1, 0, 0, 1, 0, 1, 1, 0,   // 104 HKMN
          0, 0, 1, 1, 0, 1, 1, 1,   // 112 RSUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 120 Y
          0, 0, 0, 0, 0, 0, 0, 0,   // 128
          0, 0, 0, 0, 0, 0, 0, 0,   // 136
          0, 0, 0, 0, 0, 0, 0, 0,   // 144
          0, 0, 0, 0, 0, 0, 0, 0,   // 152
          0, 0, 0, 0, 0, 0, 0, 0,   // 160
          0, 0, 0, 0, 0, 0, 0, 0,   // 168
          0, 0, 0, 0, 0, 0, 0, 0,   // 176
          0, 0, 0, 0, 0, 0, 0, 0,   // 184
          0, 0, 0, 0, 0, 0, 0, 0,   // 192
          0, 0, 0, 0, 0, 0, 0, 0,   // 200
          0, 0, 0, 0, 0, 0, 0, 0,   // 208
          0, 0, 0, 0, 0, 0, 0, 0,   // 216
          0, 0, 0, 0, 0, 0, 0, 0,   // 224
          0, 0, 0, 0, 0, 0, 0, 0,   // 232
          0, 0, 0, 0, 0, 0, 0, 0,   // 240
          0, 0, 0, 0, 0, 0, 0, 0 }; // 248

static constexpr uint8_t ACGT[256] =
        { 0, 0, 0, 0, 0, 0, 0, 0,   // 0
          0, 0, 0, 0, 0, 0, 0, 0,   // 8
          0, 0, 0, 0, 0, 0, 0, 0,   // 16
          0, 0, 0, 0, 0, 0, 0, 0,   // 24
          0, 0, 0, 0, 0, 0, 0, 0,   // 32
          0, 0, 0, 0, 0, 0, 0, 0,   // 40
          0, 0, 0, 0, 0, 0, 0, 0,   // 48
          0, 0, 0, 0, 0, 0, 0, 0,   // 56
          0, 1, 0, 1, 0, 0, 0, 1,    // 64
          0, 0, 0, 0, 0, 0, 0, 0,   // 72
          0, 0, 0, 0, 1, 0, 0, 0,   // 80
          0, 0, 0, 0, 0, 0, 0, 0,   // 88
          0, 1, 0, 1, 0, 0, 0, 1,    // 96
          0, 0, 0, 0, 0, 0, 0, 0,   // 104
          0, 0, 0, 0, 1, 0, 0, 0,   // 112
          0, 0, 0, 0, 0, 0, 0, 0,   // 120
          0, 0, 0, 0, 0, 0, 0, 0,   // 128
          0, 0, 0, 0, 0, 0, 0, 0,   // 136
          0, 0, 0, 0, 0, 0, 0, 0,   // 144
          0, 0, 0, 0, 0, 0, 0, 0,   // 152
          0, 0, 0, 0, 0, 0, 0, 0,   // 160
          0, 0, 0, 0, 0, 0, 0, 0,   // 168
          0, 0, 0, 0, 0, 0, 0, 0,   // 176
          0, 0, 0, 0, 0, 0, 0, 0,   // 184
          0, 0, 0, 0, 0, 0, 0, 0,   // 192
          0, 0, 0, 0, 0, 0, 0, 0,   // 200
          0, 0, 0, 0, 0, 0, 0, 0,   // 208
          0, 0, 0, 0, 0, 0, 0, 0,   // 216
          0, 0, 0, 0, 0, 0, 0, 0,   // 224
          0, 0, 0, 0, 0, 0, 0, 0,   // 232
          0, 0, 0, 0, 0, 0, 0, 0,   // 240
          0, 0, 0, 0, 0, 0, 0, 0 }; // 248

// implicitly convert non ACGT characters to Ns when reading a fasta file
namespace seqan {
    template <>
    struct FastaIgnoreOrAssertFunctor_<Dna5>
    {
        typedef typename FastaIgnoreFunctor_<Dna5>::Type TIgnore;
        typedef TIgnore Type;
    };
}

void add_seq(std::vector<uint64_t>& sv, const std::string& s)
{
    int nc = 0; // number of characters
    int i  = 0;

    for (const auto c : s)
    {
        if (nc == 32) // full int - ca write
        {
            nc = 0;
            ++i;
        }
        sv[i] = sv[i] | (uint64_t)(nt_table[c] & 3) << (nc * 2);
        ++nc;
    }
}

void add_seq(std::vector<uint64_t>& sv, const seqan::Dna5String& s)
{
    int nc = 0; // number of characters
    int i  = 0;

    for (const auto c : s)
    {
        if (nc == 32) // full int - ca write
        {
            nc = 0;
            ++i;
        }
        sv[i] = sv[i] | (uint64_t)(seqan::ordValue(c) & 3) << (nc * 2); // TODO: assumes that A is 0
        ++nc;
    }
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& p)
{
    os << '(' << p.first << ',' << p.second << ')';
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << '[';
    for (size_t i = 0; i < v.size(); ++i)
    {
        os << v[i];
        if (i != v.size() - 1)
            os << ',';
    }
    os << ']';
    return os;
}

std::vector<uint64_t>& operator&=(std::vector<uint64_t>& lhs, const std::vector<uint64_t>& rhs)
{
    assert(lhs.size() == rhs.size());

    for (int i = 0; i < lhs.size(); ++i)
    {
        lhs[i] &= rhs[i];
    }
    return lhs;
}

std::vector<uint64_t> operator|(const std::vector<uint64_t>& lhs, const std::vector<uint64_t>& rhs)
{
    std::vector<uint64_t> res;
    res.reserve(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i] | rhs[i]);
    }
    return res;
}

std::vector<uint64_t> operator&(const std::vector<uint64_t>& lhs, const std::vector<uint64_t>& rhs)
{
    std::vector<uint64_t> res;
    res.reserve(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i] & rhs[i]);
    }
    return res;
}

std::vector<uint64_t> operator~(const std::vector<uint64_t>& v)
{
    std::vector<uint64_t> res;
    res.reserve(v.size());
    for (int i = 0; i < v.size(); ++i)
    {
        res.push_back(~v[i]);
    }
    return res;
}

uint64_t popcount(const std::vector<uint64_t>& v)
{
    uint64_t result = 0;
    for (int i = 0; i < v.size(); ++i)
    {
        result += __builtin_popcountll(v[i]);
    }
    return result;
}
