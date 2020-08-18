#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "arg_parse.hpp"

int trivialLog(const int clusters, const float p_stay, const float p_switch,
               std::vector<float>& probs, const std::string& id, std::stringstream& ss)
{
    const int length = probs.size() / clusters;

    for (int i = 0; i < probs.size(); ++i)
    {
        probs[i] = std::log(probs[i]);
    }

    // the first half stores the best probabilities/paths for ending in each cluster,
    // the second half is used as temporary values.
    std::vector<float> best_probs(2 * clusters);
    std::vector<std::vector<uint8_t> > best_paths(2 * clusters);

    for (int c = 0; c < clusters; ++c)
    {
        best_probs[c] = probs[c]; // first position
        best_paths[c].resize(length);
        best_paths[c][0] = c;     // we start in this cluster
    }

    for (int l = 1; l < length; ++l)
    {
        for (int c_to = 0; c_to < clusters; ++c_to)
        {
            for (int c_start = 0; c_start < clusters; ++c_start)
            {
                const float new_tmp = best_probs[c_start]
                                      + (c_start == c_to ? p_stay : p_switch)
                                      + probs[l * clusters + c_to];

                if (c_start == 0 || new_tmp > best_probs[clusters + c_to])
                {
                    best_paths[clusters + c_to]    = best_paths[c_start];
                    best_paths[clusters + c_to][l] = c_to;
                    best_probs[clusters + c_to]    = new_tmp;
                }
            }
        }
        // update probs and paths from tmp results
        for (int c = 0; c < clusters; ++c)
        {
            best_probs[c] = best_probs[clusters + c];
            best_paths[c] = best_paths[clusters + c];
        }
    }

    const uint8_t best_final_cluster =
        std::distance(best_probs.begin(),
                      std::max_element(best_probs.begin(), best_probs.begin() + clusters));

    // use this if you want to access the results outside of this function
    // result_log_prob = best_probs[best_final_cluster];
    // result_path = std::move(best_paths[best_final_cluster]);

    // output
    ss << id << '\t'
       << best_probs[best_final_cluster] << '\t';
    uint8_t lastCluster  = best_paths[best_final_cluster][0];
    int     occ          = 0;
    int     l            = 0;
    int     num_clusters = 0;
    for (; l < length; ++l)
    {
        while (lastCluster == best_paths[best_final_cluster][l] && l < length)
        {
            ++occ;
            ++l;
        }
        ss << (int)lastCluster << ':' << occ << ((l < length) ? ";" : "");
        lastCluster = best_paths[best_final_cluster][l];
        num_clusters++;
        occ = 1;
    }
    if (l == length)
    {
        ss << (int)lastCluster << ':' << 1;
    }

    ss << '\n';
    return(num_clusters);
}

int bolotie_find(int argc, char **argv)
{
    enum Opt_FIND { INDEX  = 'x',
                    INPUT  = 'i',
                    PROB   = 'p',
                    THREADS= 't',
                    OUTPUT = 'o',
                    RECOMB = 'r',
                    MIN_LEN= 'm' };

    ArgParse args_find("find");

    args_find.add_string(Opt_FIND::INDEX, "index", "", "path to the basename of the bolotie index.", true);
    args_find.add_string(Opt_FIND::INPUT, "input", "", "path to the FASTA-formatted file with query sequences", true);
    args_find.add_double(Opt_FIND::PROB, "prob", 0.9999, "probability of staying in the same cluster.", false);
    args_find.add_int(Opt_FIND::THREADS, "threads", 1, "number of threads to use", false);
    args_find.add_string(Opt_FIND::OUTPUT, "output", "", "output file where paths will be saved. If not specified - output will be directed to stdout", false);
    args_find.add_flag(Opt_FIND::RECOMB, "recomb", "Only write output for paths which cross multiple clusters", false);
    args_find.add_int(Opt_FIND::MIN_LEN, "minlen", 1, "minimum cumulative legth of each parent", false);

    args_find.parse_args(argc, argv);

    std::string cl = "bolotie find";
    for (int i = 0; i < argc; i++)
    {
        if (i == 0)
        {
            cl += argv[i];
        }
        else
        {
            cl += " ";
            cl += argv[i];
        }
    }

    // load the index
    std::vector<std::vector<std::vector<float> > > prob_table; // 1D - length of the genome; 2D - number of clusters; 3D - possible base values
    const std::string index_fname = args_find.get_string(Opt_FIND::INDEX);
    std::ifstream     index_source;
    index_source.open(index_fname, std::ios_base::in);
    if (!index_source)
    {
        std::cerr << "Can't open the index file!\n";
        exit(2);
    }

    std::string index_line, nt_vals, clu_val;
    while (std::getline(index_source, index_line))
    {
        std::istringstream line_ss(index_line);
        prob_table.push_back(std::vector<std::vector<float> >{});
        while (std::getline(line_ss, nt_vals, '\t'))
        {
            std::istringstream nt_ss(nt_vals);
            prob_table.back().push_back(std::vector<float>{});
            while (std::getline(nt_ss, clu_val, ','))
            {
                prob_table.back().back().push_back(std::stof(clu_val));
            }
        }
    }

    int seq_len  = prob_table.size();
    int clusters = prob_table[0][0].size();

    const std::string fasta_fname    = args_find.get_string(Opt_FIND::INPUT);
    double            p_stay         = args_find.get_double(Opt_FIND::PROB);
    double            p_switch       = 1 - p_stay;
    const int         threads        = args_find.get_int(Opt_FIND::THREADS);
    const bool        recomb_only    = args_find.get_flag(Opt_FIND::RECOMB);

    // transform probabilities to log space
    p_stay   = std::log(p_stay);
    p_switch = std::log(p_switch);

    // quickly open file to count number of clusters
    std::ifstream fasta_source;
    std::string   fasta_line;
    fasta_source.open(fasta_fname, std::ios_base::in);
    if (!fasta_source)
    {
        std::cerr << "Can't open file!\n";
        exit(2);
    }

    // probabilities of all clusters are stored in one vector.
    // the first `clusters` values are the first probabilities of each cluster.
    std::vector<std::string>         ids;
    std::vector<std::vector<float> > probs;

    probs.resize(1000); // compute this many sequences in parallel
    ids.reserve(probs.size());

    for (int i = 0; i < probs.size(); ++i)
    {
        probs[i].reserve(seq_len * clusters);
    }

    int  seq_id      = 0;
    bool first_line  = true;
    int  cur_seq_len = 0;

    // setup output stream - http://www.gotw.ca/gotw/048.htm
    std::ostream *out_fp = &std::cout;
    std::ofstream file_out_fp;
    if (args_find.is_set(Opt_FIND::OUTPUT))
    {
        file_out_fp.open(args_find.get_string(Opt_FIND::OUTPUT));
        out_fp = &file_out_fp;
    }

    while (getline(fasta_source, fasta_line))
    {
        if (fasta_line[0] == '>')
        {
            if (!first_line)
            {
                // finished reading one entire sequence
                ++seq_id;
                if (seq_id == probs.size())
                {
                    #pragma omp parallel for num_threads(threads)
                    for (int i = 0; i < seq_id; ++i)
                    {
                        // compute  and output
                        std::stringstream ss;
                        int num_clusters_in_path = trivialLog(clusters, p_stay, p_switch, probs[i], ids[i], ss);
                        probs[i].clear();

                        #pragma omp critical
                        {
                            if (!(recomb_only && num_clusters_in_path == 1))
                            {
                                out_fp->write(ss.str().c_str(), ss.str().size());
                            }
                        }
                    }
                    ids.clear();
                    seq_id = 0;
                }
                if (cur_seq_len != prob_table.size())
                {
                    std::cerr << "incompatible sequence lengths\n";
                    exit(3);
                }
                cur_seq_len = 0;
            }
            // continue with next sequence
            ids.push_back(fasta_line.substr(1, fasta_line.size()));
        }
        else
        {
            for (int i = 0; i < fasta_line.size(); i++)
            {
                if (ACGT[fasta_line[i]])
                {
                    int nt = nt_table[fasta_line[i]];
                    for (int c = 0; c < clusters; c++)
                    {
                        probs[seq_id].push_back(prob_table[cur_seq_len + i][nt][c]); // cur_seq_len+i is position of nt for multi-line fasta (i - pos within line)
                    }
                }
                else if (IUPAC[fasta_line[i]])
                {
                    for (int c = 0; c < clusters; c++)
                    {
                        probs[seq_id].push_back(1.0f / clusters); // cur_seq_len+i is position of nt for multi-line fasta (i - pos within line)
                    }
                }
                else
                {
                    std::cerr << "unknown nucleotide\n";
                    exit(3);
                }
            }
            cur_seq_len += fasta_line.size();
        }
        first_line = false;
    }

    // compute and output for last sequence in file
    // trivialLog(clusters, p_stay, p_switch, probs);
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < seq_id + 1; ++i) // NOTE: +1 because we didn't increment it since there was no new > to see
    {                                    // compute  and output
        std::stringstream ss;
        int num_clusters_in_path = trivialLog(clusters, p_stay, p_switch, probs[i], ids[i], ss);

        #pragma omp critical
        {
            if (!(recomb_only && num_clusters_in_path == 1))
            {
                out_fp->write(ss.str().c_str(), ss.str().size());
            }
        }
    }

    file_out_fp.close();
    fasta_source.close();

    return(0);
}
