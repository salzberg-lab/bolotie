#include <vector>
#include <fstream>
#include <unordered_set>
#include <regex>

#include "common.hpp"
#include "arg_parse.hpp"

struct Sequence
{
    uint16_t cluster_id;
    char date[11]; // YYYY-MM-DD plus null terminator
    std::vector<uint64_t> seq;
    std::vector<uint64_t> n_mask;

    Sequence(const uint16_t cluster_id, const char* date, const std::vector<uint64_t> & seq, const std::vector<uint64_t> & n_mask)
    {
        this->cluster_id = cluster_id;
        this->seq = seq;
        this->n_mask = n_mask;
        strcpy(this->date, date);
    }
};

// possible date formats observed (in 68973 sequences):
// *_YYYY-MM-DD (67218 seq)
// *_YYYY-MM (940 seq)
// *_YYYY (815 seq)

static std::regex regExDateYMD("_\\d{4}-\\d{2}-\\d{2}$");
static std::regex regExDateYM("_\\d{4}-\\d{2}$");
static std::regex regExDateY("_\\d{4}$");

inline void formatDate(std::string& date)
{
    std::smatch match;
    // remove leading underscore "_" and append earliest possible month/date (to be as conservative as possible)
    if (std::regex_search(date, match, regExDateYMD))
        date = std::string(match.str()).substr(1);
    else if (std::regex_search(date, match, regExDateYM))
        date = std::string(match.str()).substr(1) + "-31";
    else if (std::regex_search(date, match, regExDateY))
        date = std::string(match.str()).substr(1) + "-12-31";
}

inline double kimura(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b, const uint32_t nbr_n_bases)
{
    assert(a.size() == b.size());

    constexpr uint64_t bitmask = 0b0101010101010101010101010101010101010101010101010101010101010101; // every 2nd bit is a one

    const int seq_len         = a.size() * 32;
    int       num_mismatches  = 0;
    int       num_transitions = 0;

    for (int i = 0; i < a.size(); ++i)
    {
        const uint64_t c = a[i] ^ b[i];          // xor to get the differences for each bit
        uint64_t d = c | (c >> 1);               // will make sure that the 2nd / right bit is 1 if any of the two bits are 1 (i.e. chars don't match)
        d &= bitmask;                            // remove any bits on the left, so we can count the ones on the right
        num_mismatches += __builtin_popcountll(d); // # mismatches (transition + transversion)

        // a^b: 2nd bit has to be 0 (didn't change), 1st bit has to be 1
        //      hence: (~c) & (c >> 1). if it is a transition, the 2nd bit will be a 1
        const uint64_t e = ((~c) & (c >> 1)) & bitmask; // ignore the right bit (i.e., negate the bitmask)
        num_transitions += __builtin_popcountll(e); // # mismatches (transition + transversion)
    }

    num_mismatches += nbr_n_bases;
    // this will count Ns as transversions
    const int num_transversions = num_mismatches - num_transitions;

    const double frac_transitions   = (double)num_transitions / (double)seq_len;
    const double frac_transversions = (double)num_transversions / (double)seq_len;

    double res = (-0.5) * log(((1.0 - (2.0 * frac_transitions)) - frac_transversions) * sqrt(1.0 - (2.0 * frac_transversions)));
    if (res < 0)
        res = 0;

    return(res);
}

inline int readClustersFile(const std::string & path, std::unordered_map<std::string, Sequence> & seqs)
{
    std::ifstream source;
    source.open(path, std::ios_base::in);
    if (!source)
    {
        std::cerr << "Can't open cluster file!\n";
        exit(1);
    }

    uint16_t max_cluster_id = 0;

    std::string line;
    std::string seq_id;
    uint16_t c_id;
    while (getline(source, line))
    {
        std::istringstream in(line);
        in >> seq_id >> c_id;
        std::string seq_date = seq_id;
        formatDate(seq_date);

        seqs.emplace(seq_id, Sequence(c_id, seq_date.c_str(), std::vector<uint64_t>(), std::vector<uint64_t>()));
        if (c_id > max_cluster_id)
            max_cluster_id = c_id;
    }
    source.close();

    return max_cluster_id;
}

inline int readFastaFile(const std::string & path, std::unordered_map<std::string, Sequence> & seqs)
{
    seqan::SeqFileIn seqFileIn(path.c_str());
    bool first_seq = true;
    seqan::CharString fasta_id;
    seqan::Dna5String fasta_seq;
    int seq_length;
    while (!seqan::atEnd(seqFileIn))
    {
        seqan::readRecord(fasta_id, fasta_seq, seqFileIn);

        if (first_seq)
        {
            seq_length = seqan::length(fasta_seq);
            first_seq = false;
        }
        else if (seq_length != seqan::length(fasta_seq))
        {
            std::cerr << "Different sequence lengths.!\n"
                      << fasta_id << " is " << seqan::length(fasta_seq) << "bp long, "
                      << "but expected " << seq_length << "bp.\n";
            exit(2);
        }

        std::vector<uint64_t> int_seq;
        int_seq.resize(((seq_length - 1) / 32) + 1, 0);
        add_seq(int_seq, fasta_seq);

        std::vector<uint64_t> n_mask;
        n_mask.resize(((seq_length - 1) / 32) + 1, 0);
        for (uint32_t i = 0; i < length(fasta_seq); ++i)
        {
            // want to set all N to A (00), everything else to T (11)
            // if we then mask it with logical AND, N will be removed, everything else remains unchanged
            // we can overwrite fasta_seq since it is not needed beyong this point
            fasta_seq[i] = (fasta_seq[i] == seqan::Dna5('N') ? seqan::Dna5('A') : seqan::Dna5('T'));
        }
        add_seq(n_mask, fasta_seq);

        auto & s = seqs.at(std::string(seqan::toCString(fasta_id)));
        s.seq = std::move(int_seq);
        s.n_mask = std::move(n_mask);
    }

    return seq_length;
}

void readPathFile(const std::string & path, std::unordered_map<std::string, std::vector<std::pair<uint16_t, uint16_t> > > & paths)
{
    std::ifstream source;
    source.open(path, std::ios_base::in);
    if (!source)
    {
        std::cerr << "Can't open paths file!\n";
        exit(1);
    }

    double score;
    std::string line;
    std::string token;
    std::string seq_id;
    std::string path_string;
    while (getline(source, line))
    {
        std::istringstream in(line);
        in >> seq_id >> score >> path_string;

        std::istringstream path_ss(path_string);

        std::vector<std::pair<uint16_t, uint16_t> > paths_vector;
        while(std::getline(path_ss, token, ';'))
        {
            const uint32_t delimiter_pos = token.find(':');
            const uint32_t cluster_id = std::stoi(token.substr(0, delimiter_pos));
            const uint32_t length = std::stoi(token.substr(delimiter_pos + 1));
            paths_vector.emplace_back(cluster_id, length);
        }
        paths.emplace(seq_id, paths_vector);
    }
    source.close();
}

int bolotie_parents(int argc, char **argv)
{
    enum Opt_PARENTS
    {
        FASTA_IN    = 'f',
        CLUSTERS_IN = 'c',
        PATHS_IN    = 'p',
        OUTFNAME    = 'o'
    };

    ArgParse args_parents("parents");

    args_parents.add_string(Opt_PARENTS::FASTA_IN, "fasta", "", "path to the sequences in FASTA format", true);
    args_parents.add_string(Opt_PARENTS::CLUSTERS_IN, "clusters", "", "path to the cluster file", true);
    args_parents.add_string(Opt_PARENTS::PATHS_IN, "paths", "", "path to the paths file", true);
    args_parents.add_string(Opt_PARENTS::OUTFNAME,"out","-","path to the output file",false);
    args_parents.parse_args(argc, argv);

    std::ofstream     out_ss;
    if(args_parents.is_set(Opt_PARENTS::OUTFNAME))
    {
        out_ss.open(args_parents.get_string(Opt_PARENTS::OUTFNAME), std::ios_base::out);
        if (!out_ss)
        {
            std::cerr << "Can't open output file for writing parents!\n";
            exit(2);
        }
    }

    std::unordered_map<std::string, Sequence> seqs;
    const int max_cluster_id = readClustersFile(args_parents.get_string(Opt_PARENTS::CLUSTERS_IN), seqs);
    std::cerr << "\033[1;35mCluster file loaded ...\033[0m\n";

    const int seq_len = readFastaFile(args_parents.get_string(Opt_PARENTS::FASTA_IN), seqs);
    std::cerr << "\033[1;35mFasta file loaded ...\033[0m\n";

    std::unordered_map<std::string, std::vector<std::pair<uint16_t, uint16_t> > > paths;
    readPathFile(args_parents.get_string(Opt_PARENTS::PATHS_IN), paths);
    std::cerr << "\033[1;35mPaths file loaded ...\033[0m\n";

    std::vector<std::pair<std::string, double> > kimura_results;

    // unmasked[cluster_id] = vector of intervals that are unmasked
    // ranges are 0-based and half-open, i.e., [a,b)
    std::vector<std::pair<uint32_t, uint32_t> > unmasked;

    std::string out_str;
    bool write_str;
    for (const auto & p : paths)
    {
        write_str = true;
        out_str = p.first + ",";

        for (int c_id = 0; c_id < max_cluster_id + 1; ++c_id)
        {
            // convert to unmasked intervals
            uint32_t cum_bases = 0;
            for (const auto& entry : p.second)
            {
                if (entry.first == c_id) // unmask
                    unmasked.emplace_back(cum_bases, cum_bases + entry.second);
                cum_bases += entry.second;
            }

            if (unmasked.empty())
                continue; // this cluster does not occur in the recombinant -> skip cluster

            out_str += std::to_string(c_id) + "," + std::to_string(unmasked[0].first) + "," + std::to_string(unmasked[0].second) + ",";

            // create the mask
            std::string masked_str(seq_len, 'A'); // A == 00 (masked positions)
            for (const auto& u : unmasked)
            {
                for (uint32_t i = u.first; i < u.second; ++i)
                {
                    masked_str[i] = 'T'; // T == 11 (unmasked position)
                }
            }
            std::vector<uint64_t> masked_seq;
            masked_seq.resize(((seq_len - 1) / 32) + 1, 0);
            add_seq(masked_seq, masked_str); // transform std::string into bitmask

            const Sequence & recombinant = seqs.at(p.first);

            for (const auto& seq : seqs)
            {
                // wrong cluster or sequence newer than recombinant
                if (seq.second.cluster_id != c_id || strcmp(recombinant.date, seq.second.date) <= 0)
                    continue;

                // no fasta available
                if (seq.second.n_mask.size() == 0)
                    continue;

                std::vector<uint64_t> masked_seq_copy = masked_seq;
                masked_seq_copy &= seq.second.n_mask; // merge masking of cluster region and recombinants Ns and now: potential parent's Ns
                masked_seq_copy &= recombinant.n_mask; // merge masking of cluster region and recombinants Ns

                std::vector<uint64_t> recombinant_seq_copy = recombinant.seq; // copy on purpose (for masking)
                std::vector<uint64_t> seq_copy = seq.second.seq; // copy on purpose (for masking)

                // count Ns in unmasked region to correct kimura distance computed on DNA4
                // each character is either A (00) or T(11). no other patterns/chars in those masks:
                const std::vector<uint64_t> n_vector = ((~recombinant.n_mask | ~seq.second.n_mask) & masked_seq);
                const uint32_t nbr_n_bases = popcount(n_vector)/2;

                seq_copy &= masked_seq_copy; // mask the sequence
                recombinant_seq_copy &= masked_seq_copy; // mask the sequence

                const double k = kimura(seq_copy, recombinant_seq_copy, nbr_n_bases);
                kimura_results.emplace_back(seq.first, k);
            }

            std::sort(kimura_results.begin(), kimura_results.end(), [](const auto &l, const auto &r) {
                return l.second < r.second;
            });

            // output X most similar sequences
            // TODO: Ales wants me to check whether ABA recombinants are outputted correctly in csv, and not only AB recombinants
            if(kimura_results.empty())
            {
                write_str = false;
            }
            else{
                out_str += kimura_results[0].first + "," + std::to_string(kimura_results[0].second) + ",";
            }

            kimura_results.clear();
            unmasked.clear();
        }

        if (write_str)
        {
            out_str.pop_back();
            out_ss << out_str << '\n';
        }
    }

    out_ss.close();

    return 0;
}
