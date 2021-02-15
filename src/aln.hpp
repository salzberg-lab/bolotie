#include <iostream>
#include "ksw2.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <fstream>
#include <map>
#include <math.h>

#include "common.hpp"
#include "arg_parse.hpp"

#define CMATCH  0
#define CINS    1
#define CDEL    2
#define CREF_SKIP   3
#define CSOFT_CLIP  4
#define CHARD_CLIP  5
#define CPAD    6
#define CEQUAL  7
#define CDIFF   8

struct Parameters{
    int match = 5;
    int mismatch = -4;
    int gapopen = 12;
    int gapextend = 2;
    int trim_len = 100;
    int num_threads = 1;
    int min_len = 0;
    int max_len = INT_MAX;
    double fracn = FLT_MAX;
    int cigarops = INT_MAX;
};

struct Alignment{
    std::string ref_name = "";
    int start_clip = 0;
    int end_clip = 0;
    std::vector<std::tuple<int,uint8_t,uint8_t> > snvs;
    std::vector<std::pair<int,std::string> > inss;
    std::vector<std::pair<int,std::string> > dels;
    std::vector<int> ambigs; // ambiguous bases (non-ACGT) outside of soft clipping
    std::vector<std::pair<int,uint8_t> > cigar;
    bool pass = true;
};

class Finder{
public:
    explicit Finder() = default;
    ~Finder() = default;
    void set_scores(int sc_mch,int sc_mis,int gapo,int gape){
        this->sc_mch = sc_mch; // ?
        this->sc_mis = sc_mis; // mismatch?
        this->gapo = gapo; // gap open penalty
        this->gape = gape; // gap extend penalty
    }
    void set_trim_len(int len){
        this->trim_len = len;
    }
    bool find(std::string &t,std::string& q,ksw_extz_t& ez){
        align(t.c_str(),q.c_str(),ez);
        return true;
    }
    // extract:
    // 1. clipping with respect to the reference
    // 2. SNVs
    // 3. insertions
    // 4. deletions
    int parse_aln(ksw_extz_t& ez, Alignment& aln, std::string& t, std::string& q, const Parameters& params)
    {
        // print_cigar(ez);
        // first find the position of the last matched base on the query
        int last_m_q = ez.max_q;
        for (uint16_t c = ez.n_cigar - 1; c >= 0; --c)
        {
            int opcode =  ez.cigar[c] & 0xf;
            int oplen = ez.cigar[c] >> 4;

            if (opcode == CDEL)
                continue;
            else if (opcode == CINS)
                last_m_q-=oplen;
            else if (opcode == CMATCH)
                break;
            else
            {
                std::cerr << "unidentified cigar operation\n";
                exit(3);
            }
        }

        // need the template coordinate of the first match
        int cur_pos_q = 0;
        int cur_pos_t = 0;

        for (uint16_t c = 0;c < ez.n_cigar; ++c)
        {
            int opcode = ez.cigar[c] & 0xf;
            int oplen = ez.cigar[c] >> 4;

            if (opcode == CMATCH)
                break;
            else if (opcode == CINS)
                continue;
            else if (opcode == CDEL)
                cur_pos_t += oplen;
            else
            {
                std::cerr << "unidentified cigar operation\n";
                exit(3);
            }
        }

        int first_match_pos_t = cur_pos_t;
        int last_match_pos_t = ez.max_t;
        if (last_match_pos_t - first_match_pos_t < params.min_len)
            return -1;

        // next find the where to begin trimming the sequence
        int end_trim_start_pos = last_m_q-this->trim_len;
        int start_trim_end_pos = this->trim_len;

        // now can begin parsing
        cur_pos_q = 0;
        cur_pos_t = 0;

        for (uint16_t c = 0; c < ez.n_cigar; ++c)
        {
            int opcode = ez.cigar[c] & 0xf;
            int oplen = ez.cigar[c] >> 4;

            for(int i = 0; i < oplen; ++i)
            {
                if (cur_pos_q < start_trim_end_pos || cur_pos_q >= end_trim_start_pos) // within trimming range
                {
                    if (opcode == CMATCH)
                    {
                        if (!aln.cigar.empty() && aln.cigar.back().second == CMATCH)
                            ++aln.cigar.back().first;
                        else
                            aln.cigar.emplace_back(1, CMATCH);
                        ++cur_pos_q;
                        ++cur_pos_t;
                    }
                    else if (opcode == CDEL)
                    {
                        if (!aln.cigar.empty() && aln.cigar.back().second == CMATCH)
                            ++aln.cigar.back().first;
                        else
                            aln.cigar.emplace_back(1, CMATCH);
                        ++cur_pos_t;
                    }
                    else if (opcode == CINS)
                    {
                        ++cur_pos_q;
                        ++start_trim_end_pos;
                    }
                    else
                    {
                        std::cerr << "unidentified cigar operation\n";
                        exit(3);
                    }
                }
                else
                {
                    if (opcode == CMATCH)
                    {
                        if (q[cur_pos_q] != t[cur_pos_t])
                        {
                            aln.snvs.emplace_back(cur_pos_t, q[cur_pos_q], t[cur_pos_t]);
                            aln.cigar.emplace_back(1, 'Z');
                        }
                        else
                        {
                            if (!aln.cigar.empty() && aln.cigar.back().second == CMATCH)
                                ++aln.cigar.back().first;
                            else
                                aln.cigar.emplace_back(1, CMATCH);
                        }
                        ++cur_pos_q;
                        ++cur_pos_t;
                    }
                    else if (opcode == CINS)
                    {
                        if (!aln.cigar.empty() && aln.cigar.back().second == CINS)
                        {
                            ++aln.cigar.back().first;
                            aln.inss.back().second += q[cur_pos_q];
                        }
                        else
                        {
                            aln.cigar.emplace_back(1, CINS);
                            aln.inss.emplace_back(cur_pos_t, std::string(1, q[cur_pos_q]));
                        }
                        ++cur_pos_q;
                    }
                    else if (opcode == CDEL)
                    {
                        if (!aln.cigar.empty() && aln.cigar.back().second == CDEL)
                        {
                            ++aln.cigar.back().first;
                            aln.dels.back().second += t[cur_pos_t];
                        }
                        else
                        {
                            aln.cigar.emplace_back(1, CDEL);
                            aln.dels.emplace_back(cur_pos_t, std::string(1, t[cur_pos_t]));
                        }
                        ++cur_pos_t;
                    }
                    else
                    {
                        std::cerr << "unidentified cigar operation\n";
                        exit(3);
                    }
                }
            }
        }
        return 1;
    }

protected:
    void print_cigar(ksw_extz_t& ez)
    {
        for (int i = 0; i < ez.n_cigar; ++i)
        {
            printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
        }
        putchar('\n');
    }

    void align(const char *tseq, const char *qseq, ksw_extz_t& ez) // from the ksw2 example
    {
        int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
        int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
        int tl = strlen(tseq), ql = strlen(qseq);
        uint8_t *ts, *qs, c[256];

        memset(&ez, 0, sizeof(ksw_extz_t));
        memset(c, 4, 256); // TODO: make this global static const
        c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
        c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
        ts = (uint8_t*)malloc(tl);
        qs = (uint8_t*)malloc(ql);
        for (int i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
        for (int i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
        const int flag = KSW_EZ_EXTZ_ONLY;
        ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, flag, &ez);
        // ksw_gg2_sse(0,ql,qs,tl,ts,5,mat,gapo,gape,-1,&ez.m_cigar,&ez.n_cigar,&ez.cigar);
        // print_cigar(ez);
        free(ts);
        free(qs);
    }

    int sc_mch = 5, sc_mis = -4, gapo = 12, gape = 2;
    int trim_len = 0;

};

// clean a sequence - convert all to UPPER case and remove non-IUPAC codes
int clean_seq(std::string& seq, const Parameters& params)
{
    // remove ambiguity from the start
    for (auto itr = seq.begin(); itr != seq.end(); )
    {
        if (!ACGT[*itr])
            itr = seq.erase(itr);
        else // found the first ACGT base
            break;
    }

    // remove ambiguity from the end
    for (auto itr = seq.end() - 1; itr != seq.begin(); )
    {
        if (!ACGT[*itr])
        {
            itr = seq.erase(itr);
            --itr;
        }
        else // found the first ACGT base
        {
            break;
        }
    }

    int num_amb = 0;
    for (auto itr = seq.begin(); itr != seq.end(); )
    {
        if (!IUPAC[*itr]) // not a valid IUPAC code - skip
        {
            // ++num_amb;  // TODO: for ales. do we want to remove it?
            itr = seq.erase(itr);
            continue;
        }
        *itr = std::toupper(*itr);
        ++itr;
    }

    // now check length of the resulting sequence making sure it is within a threshold margin of the reference
    if (seq.size() < params.min_len || seq.size() > params.max_len)
        return -1;

    const float frac_amb = (float)num_amb / seq.size();
    if (frac_amb > params.fracn)
        return -1;

    return 1;
}

void write(const std::string& out_fname, const std::vector<Alignment>& alns, const Parameters& params)
{
    std::ostream *out_fp = &std::cout;
    std::ofstream file_out_fp;
    if (std::strcmp(out_fname.c_str(), "-") != 0)
    {
        file_out_fp.open(out_fname);
        out_fp = &file_out_fp;
    }

    for (const auto& aln : alns)
    {
        if (!aln.pass)
            continue;

        // TODO: remove 1s
        std::string res = aln.ref_name + "," + std::to_string(aln.ambigs.size()) + ",";
        for (const auto& v : aln.snvs)
            res += std::to_string(std::get<0>(v)) + std::string(1, std::get<1>(v)) + std::string(1, std::get<2>(v));

        res += ",";
        for (const auto& v : aln.dels)
            res += std::to_string(std::get<0>(v)) + std::get<1>(v);

        res += ",";
        for (const auto& v : aln.inss)
            res+=std::to_string(std::get<0>(v))+std::get<1>(v);

        res += ",";
        for (const auto& c : aln.cigar)
        {
            res += std::to_string(std::get<0>(c));
            switch (std::get<1>(c))
            {
                case CMATCH:
                    res += "M";
                    continue;
                case CDEL:
                    res += "D";
                    continue;
                case CINS:
                    res += "I";
                    continue;
                case 'Z':
                    res += "Z";
                    continue;
                default:
                    std::cerr << "unknown opcode: " << (int)std::get<1>(c) << '\n';
                    exit(3);
            }
        }
        res += "\n";

        out_fp->write(res.c_str(), res.size());
    }

    file_out_fp.close();
}

// write variants into a separate file
void write_vars(const std::string& out_fname, const std::vector<Alignment>& alns, const Parameters& params){
    std::ostream *out_fp = &std::cout;
    std::ofstream file_out_fp;

    if (std::strcmp(out_fname.c_str(), "-") != 0)
    {
        file_out_fp.open(out_fname);
        out_fp = &file_out_fp;
    }

    for(int i = 0; i < alns.size(); ++i)
    {
        if (!alns[i].pass)
            continue;

        std::string res = "";
        for (const auto& v : alns[i].snvs)
            res += alns[i].ref_name + ",S," + std::to_string(std::get<0>(v)) + "," + std::string(1,std::get<1>(v)) + "," + std::string(1,std::get<2>(v)) + "\n";

        for (const auto& v : alns[i].dels)
            res += alns[i].ref_name + ",D," + std::to_string(std::get<0>(v)) + "," + std::get<1>(v) + ",-\n";

        for (const auto& v : alns[i].inss)
            res += alns[i].ref_name + ",I," + std::to_string(std::get<0>(v)) + ",-," + std::get<1>(v) + "\n";

        out_fp->write(res.c_str(), res.size());
    }
    file_out_fp.close();
}

void run(Alignment& aln, std::string& ref_seq, std::string& query_seq,Parameters& params){
    ksw_extz_t ez;

    Finder fnd;
    fnd.set_scores(params.match, params.mismatch, params.gapopen, params.gapextend);
    fnd.set_trim_len(params.trim_len);

    fnd.find(ref_seq, query_seq, ez);

    if (ez.n_cigar == 0 || ez.n_cigar>params.cigarops)
        aln.pass = false;

    const int res = fnd.parse_aln(ez,aln,ref_seq,query_seq,params);
    if (res <= 0)
        aln.pass = false;

    free(ez.cigar);
}

int bolotie_aln(int argc, char **argv)
{
    enum Opt_ALN { REFERENCE   = 'x',
                   QUERIES     = 'i',
                   OUTPUT      = 'o',
                   NUM_THREADS = 't',
                   MATCH       = 'm',
                   MISMATCH    = 's',
                   GAPOPEN     = 'g',
                   GAPEXTEND   = 'e',
                   TRIMLEN     = 'l',
                   FRACN       = 'f',
                   FRACLEN     = 'n',
                   CIGAROPS    = 'c'};

    ArgParse args_aln("aln");

    args_aln.add_string(Opt_ALN::REFERENCE, "reference", "", "path to the reference sequence in FASTA format", true);
    args_aln.add_string(Opt_ALN::QUERIES, "query", "", "path to the FASTA file with query sequences", true);
    args_aln.add_int(Opt_ALN::NUM_THREADS, "threads", 1, "number of threads to use", false);
    args_aln.add_string(Opt_ALN::OUTPUT, "output", "-", "output path and name", false);
    args_aln.add_int(Opt_ALN::MATCH, "match", 5, "score to assign to matches", false);
    args_aln.add_int(Opt_ALN::MISMATCH, "mismatch", -4, "score to assign to mismatches", false);
    args_aln.add_int(Opt_ALN::GAPOPEN, "gapopen", 12, "gap open score", false);
    args_aln.add_int(Opt_ALN::GAPEXTEND, "gapextend", 2, "gap extension score", false);
    args_aln.add_int(Opt_ALN::TRIMLEN, "trimlen", 100, "number of query bases at 3' and 5' ends to be replaced with reference", false);
    args_aln.add_double(Opt_ALN::FRACN,"fracn",0.01,"fraction of the sequence that can be covered by ambiguous nucleotides",false);
    args_aln.add_double(Opt_ALN::FRACLEN,"fraclen",0.02,"permitted deviation from the reference sequence size (after trimming ambiguity from 5' and 3')",false);
    args_aln.add_int(Opt_ALN::CIGAROPS,"cigarops",25,"maximum number of cigar operations permitted for an alignment",false);

    args_aln.parse_args(argc, argv);

    Parameters params;
    params.match = args_aln.get_int(Opt_ALN::MATCH);
    params.mismatch = args_aln.get_int(Opt_ALN::MISMATCH);
    params.gapopen = args_aln.get_int(Opt_ALN::GAPOPEN);
    params.gapextend = args_aln.get_int(Opt_ALN::GAPEXTEND);
    params.trim_len = args_aln.get_int(Opt_ALN::TRIMLEN);
    params.num_threads = args_aln.get_int(Opt_ALN::NUM_THREADS);
    params.cigarops = args_aln.get_int(Opt_ALN::CIGAROPS);

    std::ifstream ref_source;
    ref_source.open(args_aln.get_string(Opt_ALN::REFERENCE));
    if (!ref_source)
    {
        std::cerr << "Can't open the reference fasta file!\n";
        exit(2);
    }

    std::ifstream query_source;
    query_source.open(args_aln.get_string(Opt_ALN::QUERIES));
    if (!ref_source)
    {
        std::cerr << "Can't open the query fasta file!\n";
        exit(2);
    }

    // load reference
    std::string ref_seqid, ref_seq, ref_line;
    while (std::getline(ref_source, ref_line))
    {
        if (ref_line[ref_line.size() - 1] == '\r')
        {
            ref_line.pop_back();
        }
        if (ref_line[0] == '>')
        {
            ref_seqid = ref_line.substr(1, ref_line.size());
            ref_seq = "";
        }
        else
        {
            ref_seq += ref_line;
        }
    }

    int clean_res;
    clean_res = clean_seq(ref_seq, params);
    if (clean_res <= 0)
    {
        std::cerr << "reference sequence does not pass cleaning\n";
        exit(3);
    }

    int seq_diff = (int)(ref_seq.size() * args_aln.get_double(Opt_ALN::FRACLEN)); // permitted number of bases in trimmed sequence by which it differs from the reference
    params.min_len = ref_seq.size() - seq_diff;
    params.max_len = ref_seq.size() + seq_diff;
    params.fracn = args_aln.get_double(Opt_ALN::FRACN);

    std::string query_seqid, query_seq, query_line;
    std::vector<Alignment> alns;
    std::vector<std::pair<int,std::string> > seqs; // holds enough data for the number of requested threads. 1st value - index in alns vector; 2nd - sequence
    alns.reserve(80000);
    seqs.reserve(80000);
    while (std::getline(query_source, query_line)) // parse fasta and add to sequence
    {
        if (query_line[query_line.size() - 1] == '\r')
            query_line.pop_back();

        if (query_line[0] == '>')
        {
            if (!query_seq.empty()) // not the first sequence - can align and write
            {
                clean_res = clean_seq(query_seq,params);
                if (clean_res > 0)
                {
                    seqs.emplace_back(alns.size()-1, query_seq);
                }

                if (seqs.size() == 1000)
                {
                    #pragma omp parallel for num_threads(params.num_threads)
                    for (int i = 0; i < seqs.size(); i++)
                    {
                        std::cout<<alns[seqs[i].first].ref_name<<std::endl;
                        run(alns[seqs[i].first], ref_seq, seqs[i].second,params);
                    }
                    seqs.clear();
                }
            }
            query_seq = "";
            query_seqid = query_line.substr(1, query_line.size());
            alns.emplace_back();
            alns.back().ref_name = query_seqid;
        }
        else{
            query_seq += query_line;
        }
    }

    ref_source.close();
    query_source.close();

    // run the last sequence
    clean_res = clean_seq(query_seq,params);
    if (clean_res > 0)
        seqs.emplace_back(alns.size() - 1, query_seq);

    #pragma omp parallel for num_threads(params.num_threads)
    for(int i = 0; i < seqs.size(); i++)
    {
        run(alns[seqs[i].first], ref_seq, seqs[i].second, params);
    }

    std::cerr << "writing output\n";
    write(args_aln.get_string(Opt_ALN::OUTPUT) + ".csv", alns, params);
    write_vars(args_aln.get_string(Opt_ALN::OUTPUT) + ".vars.csv", alns, params);

    return 0;
}
