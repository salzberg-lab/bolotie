#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>

#include "common.hpp"
#include "arg_parse.hpp"

int bolotie_build(int argc, char **argv)
{
    enum Opt_BUILD { INDEX = 'x',
        MSA   = 'i',
        CLUS  = 'c' };

    ArgParse args_build("build");

    args_build.add_string(Opt_BUILD::INDEX, "input", "", "path to the basename of the bolotie index.", true);
    args_build.add_string(Opt_BUILD::MSA, "msa", "", "path to the MSA file", true);
    args_build.add_string(Opt_BUILD::CLUS, "clu", "", "path the the clusters file", true);

    args_build.parse_args(argc, argv);

    std::string cl = "bolotie build";
    for (int i = 0; i < argc; ++i)
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

    // load clusters
    std::unordered_map<std::string, uint32_t> refname_uid;                            // holds reference name to unique integer id conversion table for fast downstream accesss
    std::vector<uint32_t> refid2cluid;                                      // position in the vector is the id of the reference and the value is the cluster assignment
    std::vector<int>      clusters;                                         // cluster sizes (the size of the vector is the number of clusters)
    std::ifstream         cluster_source;
    cluster_source.open(args_build.get_string(Opt_BUILD::CLUS), std::ios_base::in);
    if (!cluster_source)
    {
        std::cerr << "Can't open the cluster file!\n";
        exit(2);
    }
    std::string line;
    int         tab_pos = 0;
    std::string refname;
    int         clu_id   = -1;
    int         num_refs = 0;

    while (std::getline(cluster_source, line))
    {
        tab_pos = line.find('\t');
        refname = line.substr(0, tab_pos);
        clu_id  = std::stoi(line.substr(tab_pos, line.size()));
        if (clu_id < 0)
        {
            std::cerr << "incorrect cluster assignment\n";
            exit(3);
        }

        while (clusters.size() < (clu_id + 1))
        {
            clusters.push_back(0);
        }
        ++clusters[clu_id];

        const bool ru_result = refname_uid.emplace(refname, num_refs).second;
        if (!ru_result)
        {
            std::cerr << "duplicate reference names\n";
            exit(3);
        }

        // assign cluster to sequence
        refid2cluid.push_back(clu_id); // no need to reserve since no duplicate names and uid is incremental
        num_refs++;
    }
    cluster_source.close();

    // initialize the structure to hold counts of bases in each cluster
    std::vector<std::vector<std::vector<float> > > probs_tmp;    // 1D - length of the genome; 2D - number of clusters; 3D - possible base values
    probs_tmp.reserve(29003);
    std::vector<std::vector<std::vector<int> > > base_counts;    // same but only holds raw base counts per cluster
    base_counts.reserve(29003);


//    std::vector<std::vector<float> > clu_weights;
//    clu_weights.reserve(29003);
//// compute normalized weight of sequences for each cluster (num_seqs_cluster/total_num_seqs)
//    float norm_factor = (float)num_refs/(float)clusters.size();
//    std::vector<float> clu_seq_weight;
//    for(auto& clu : clusters){
//        clu_seq_weight.push_back(norm_factor/clu);
//    }

    std::ifstream msa_source;
    msa_source.open(args_build.get_string(Opt_BUILD::MSA), std::ios_base::in);
    if (!msa_source)
    {
        std::cerr << "Can't open the MSA file!\n";
        exit(2);
    }
    std::string fasta_line = "";
    int         seq_len = 0, cur_seq_len = 0;

    refname = "";
    int cur_cluid = -1;
    int refid = -1;

    int refno = 0;
    while (std::getline(msa_source, fasta_line))  // parse fasta and add to sequence
    {
        if (fasta_line[0] == '>')
        {
            if (!refno) // not first sequence
            {
                if (cur_seq_len != seq_len)
                {
                    std::cerr << "sequences of different lengths\n";
                    exit(3);
                }
            }
            else
            {
                seq_len = cur_seq_len;
            }
            cur_seq_len = 0;
            refname     = fasta_line.substr(1, fasta_line.size());
            auto ru_it = refname_uid.find(refname);
            if (ru_it == refname_uid.end())
            {
                std::cerr << "skipping unknown sequence: " << refname << '\n';
                continue;
            }
            refid     = ru_it->second;
            cur_cluid = refid2cluid[refid];
            refno++;
        }
        else
        {
            for (int i = 0; i < fasta_line.size(); i++)
            {
                if (refno == 1) // first reference - need to setup each position
                {
                    probs_tmp.emplace_back();
                    base_counts.emplace_back();
                    for (int c = 0; c < clusters.size(); c++)
                    {
                        probs_tmp.back().emplace_back(4, 0);
                        base_counts.back().emplace_back(4, 0);
                    }
                }
                if (ACGT[fasta_line[i]])
                {
                    int nt = nt_table[fasta_line[i]];
                    base_counts[cur_seq_len + i][cur_cluid][nt] += 1;
                }
                else if (!IUPAC[fasta_line[i]])
                {
                    std::cerr << "unknown nucleotide\n";
                    exit(3);
                }
            }

            cur_seq_len += fasta_line.size();
        }
    }
    msa_source.close();

    // TODO:
    //   add priors to the probabilities:
    //     1. fraction of sequences in cluster which have this variant
    //     2. transversion vs transition (one more likely to occur by chance since it is less likely to change amino acids)

    const std::string probs_fname = args_build.get_string(Opt_BUILD::INDEX) + ".probs";
    std::ofstream     probs_ss(probs_fname, std::ios_base::out);
    if (!probs_ss)
    {
        std::cerr << "Can't open index file for writing probabilities!\n";
        exit(2);
    }

    const std::string counts_fname = args_build.get_string(Opt_BUILD::INDEX) + ".counts";
    std::ofstream     counts_ss(counts_fname, std::ios_base::out);
    if (!counts_ss)
    {
        std::cerr << "Can't open index file for writing counts!\n";
        exit(2);
    }

    const std::string totals_fname = args_build.get_string(Opt_BUILD::INDEX) + ".totals";
    std::ofstream     totals_ss(totals_fname, std::ios_base::out);
    if (!totals_ss)
    {
        std::cerr << "Can't open index file for writing totals!\n";
        exit(2);
    }

    // now to compute conditional probabilities from the counts
    float min_prob       = 1.0 / (float)num_refs; // minimum possible probability value to use instead of 0
    int total_nt_raw_count = 0;
    float cond_prob      = 0;
    float prob_nt = 0;
    for (int p = 0; p < probs_tmp.size(); ++p) // for each position on the genome
    {
        // compute total number of sequences with non-ambiguous code at this position
        int total_acgt = 0;
        std::vector<int> clu_acgt(clusters.size(), 0);
        for (int c = 0; c < clusters.size(); c++){
            for (int nt = 0; nt < 4; nt++){
                total_acgt += base_counts[p][c][nt];
                clu_acgt[c] += base_counts[p][c][nt];
            }
        }

        // compute normalized weight of sequences for each cluster (num_seqs_cluster/total_num_seqs)
        float norm_factor = (float)total_acgt/(float)clusters.size();
        std::vector<float> clu_seq_weight;
        for(auto& clu : clu_acgt)
        {
            clu_seq_weight.push_back(norm_factor/clu);
        }

        for (int nt = 0; nt < 4; nt++)           // for each base
        {
            prob_nt = 0;
            total_nt_raw_count = 0;

            // get total number of times a base occurs across all clusters
            for (int c = 0; c < clusters.size(); c++)
            {
                prob_nt += base_counts[p][c][nt] * clu_seq_weight[c];
                total_nt_raw_count += base_counts[p][c][nt];
            }
            prob_nt = prob_nt == 0 ? 1.0 / (float)num_refs : prob_nt;
            for (int c = 0; c < clusters.size(); c++)
            {
                cond_prob = (base_counts[p][c][nt] * clu_seq_weight[c]) / prob_nt;
                probs_ss << (cond_prob == 0 ? min_prob : cond_prob) << ",";
                counts_ss << base_counts[p][c][nt] << ",";
            }
            probs_ss.seekp(-1, std::ios_base::end);
            counts_ss.seekp(-1, std::ios_base::end);
            probs_ss << "\t";
            counts_ss << "\t";
            totals_ss << total_nt_raw_count << "\t";
        }
        probs_ss.seekp(-1, std::ios_base::end);
        counts_ss.seekp(-1, std::ios_base::end);
        totals_ss.seekp(-1, std::ios_base::end);
        probs_ss << '\n';
        counts_ss << '\n';
        totals_ss << '\n';
    }

    probs_ss.close();
    counts_ss.close();
    totals_ss.close();

    return(0);
}
