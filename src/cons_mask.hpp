#include <iostream>
#include "ksw2.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <fstream>
#include <map>
#include <htslib/htslib/vcf.h>

#include "common.hpp"
#include "arg_parse.hpp"

// TODO: does not currently support indels
struct Variants{
public:
    explicit Variants() = default;
    ~Variants() = default;

    void load_vcf_mask(const std::string& vcf_fname,const std::string& ref_seqid){ // given a vcf file, load all SNPs to be used later as a mask when building a consensus sequence
        // http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html - thank you for many useful pointers on processing vcf with htslib
        htsFile * vcf_fp = bcf_open(vcf_fname.c_str(), "r");
        if (vcf_fp == NULL)
        {
            std::cerr<<"Can't open provided VCF file"<<std::endl;
            exit(2);
        }
        bcf_hdr_t *hdr = bcf_hdr_read(vcf_fp);

        const char **seqnames = NULL;
        int nseq = 0;
        seqnames = bcf_hdr_seqnames(hdr, &nseq);
        if (seqnames == NULL)
        {
            std::cerr<<"VCF appears empty"<<std::endl;
            exit(2);
        }

        bcf1_t *rec = bcf_init();
        if (rec == NULL)
        {
            std::cerr<<"couldn't initiate a record device"<<std::endl;
            exit(2);
        }

        this->mask.resize(29903,false);
        while (bcf_read(vcf_fp, hdr, rec) == 0)
        {
            if (bcf_is_snp(rec)) {
                if(this->mask.size()<=rec->pos)
                {
                    this->mask.resize(this->mask.size()+10000,false);
                }
                this->mask[rec->pos]=true;
            }
        }
        free(seqnames);
        bcf_hdr_destroy(hdr);
        bcf_close(vcf_fp);
        bcf_destroy(rec);

        this->mask_used = true;
    }

    void add_var(const std::string& ref_seqid, const int pos, const uint8_t q, const uint8_t t)
    {
        this->ri_it = this->ref2id.insert(std::make_pair(ref_seqid,max_ref_id));
        if (this->ri_it.second) // new reference
            max_ref_id++;

        this->ir_it = this->id2ref.emplace(this->ri_it.first->second, ref_seqid);

        this->vit = this->vars.emplace(Var(pos,q,t), 1);
        if (!this->vit.second)
            this->vit.first->second++;

        // add variant to the reference
        this->rv_it = this->ref2vars.emplace(this->ri_it.first->second, std::vector<Var>{});
        this->rv_it.first->second.emplace_back(pos, q, t);

        // now also add to the pos2var
        this->pv_it = this->pos2var.emplace(pos, std::array<int,5>{0,0,0,0,0});
        if (!ACGT[q]) // if any of the ambiguous codes
            this->pv_it.first->second[4]++;
        else
            this->pv_it.first->second[nt_table[q]]++;
    }

    void write_cons(const std::string& ref_seq, const std::string& out_fname, const int min_var_count, const int amb2ref_count, const bool dedup)
    {
        // first find all position which contain high-frequency variants
        std::set<int> pos_set; // positions in which high-frequency variants occur
        std::set<int>::iterator ps_it;

        std::set<int> amb2ref_pos_set; // positions in which any ambiguous variant needs to be replaced with the reference
        std::set<int>::iterator a2r_it;

        for (auto& uid : this->ref2vars) // over each reference
        {
            for (auto& var : uid.second) // for each variant in the sequence
            {
                if (!ACGT[std::get<1>(var)]) // if a variant has an ambiguous character at this position
                {
                    this->pv_it.first = this->pos2var.find(std::get<0>(var));
                    int max_acgt_vars = std::max({this->pv_it.first->second[0],
                                                     this->pv_it.first->second[1],
                                                     this->pv_it.first->second[2],
                                                     this->pv_it.first->second[3]});
                    if (max_acgt_vars < amb2ref_count) // if there are no high-frequency non-ambiguous variants at this position
                    {
                        amb2ref_pos_set.insert(std::get<0>(var));
                    }
                }

                this->vit.first = this->vars.find(var);
                if (this->vit.first->second >= min_var_count) // passes frequency threshold
                {
                    pos_set.insert(std::get<0>(var));
                }
            }
        }

        // next add all variants at those positions to the variant sets
        std::map<std::vector<Var>,std::vector<int> > var_sets; // key is the variant set; value holds ids of all other sequences with the same sets of references
        std::pair<std::map<std::vector<Var>,std::vector<int> >::iterator,bool> vs_it;

        for(auto& uid : this->ref2vars) // over each reference
        {
            std::vector<Var> cur_passing_vars;
            for(auto& var : uid.second) // for each variant in the sequence
            {
                // first check if the variant itself passes
                this->vit.first = this->vars.find(var);
                if(this->mask_used && this->mask[std::get<0>(var)]) // if position is masked - replace with N in consensus
                {
                    ps_it = pos_set.find(std::get<0>(var));
                    if (ps_it != pos_set.end()) // position passes frequency threshold
                    {
                        cur_passing_vars.emplace_back(std::get<0>(var), 'N', std::get<2>(var));
                    }
                }
                else if(this->vit.first->second >= min_var_count) // variant passes frequency threshold
                {
                    a2r_it = amb2ref_pos_set.find(std::get<0>(var));
                    if (a2r_it != amb2ref_pos_set.end()) // should be replaced with the reference allele
                    {
                        continue;
                    }
                    else
                    {
                        cur_passing_vars.push_back(var);
                    }
                }
                else // if the variant did not pass we consider the possibility that another variant at this position may exist. To avoid bias to the reference we replace base with N
                {
                    ps_it = pos_set.find(std::get<0>(var));
                    if (ps_it != pos_set.end()) // position passes frequency threshold
                    {
                        cur_passing_vars.emplace_back(std::get<0>(var), 'N', std::get<2>(var));
                    }
                }
            }
            vs_it = var_sets.emplace(cur_passing_vars, std::vector<int>{});
            vs_it.first->second.push_back(uid.first);
        }

        // now we can write the consensus sequences
        const std::string fasta_fname = out_fname + ".fa";
        std::ofstream fasta_ss(fasta_fname, std::ios_base::out);
        if (!fasta_ss)
        {
            std::cerr << "Can't open fasta file for writing!\n";
            exit(2);
        }

        for (auto& rvar : var_sets) // for a set of variants
        {
            for (auto& r : rvar.second) // for reference
            {
                this->ir_it.first = this->id2ref.find(r);
                if (this->ir_it.first == this->id2ref.end())
                {
                    std::cerr << "reference id not found\n";
                    exit(3);
                }
                fasta_ss << '>' << this->ir_it.first->second << '\n';
                int cur_ref_pos = 0;
                std::string cons_seq = "";
                int last_var_pos = -1;
                for (auto& var : rvar.first) // for each variant
                {
                    if (last_var_pos!=-1 && last_var_pos >= std::get<0>(var))
                    {
                        std::cerr << "variants in wrong order\n";
                        exit(-1);
                    }
                    last_var_pos = std::get<0>(var);
                    cons_seq += ref_seq.substr(cur_ref_pos, std::get<0>(var) - cur_ref_pos);
                    cons_seq += std::get<1>(var);
                    cur_ref_pos = std::get<0>(var) + 1;
                    if (ref_seq[std::get<0>(var)] != std::get<2>(var))
                    {
                        std::cerr << "specified variant does not match reference\n";
                        fasta_ss.close();
                        exit(3);
                    }
                    if (cons_seq[std::get<0>(var)] != std::get<1>(var))
                    {
                        std::cerr << "incorrectly merged variants\n";
                        fasta_ss.close();
                        exit(3);
                    }
                }
                // add the last remaining bit
                cons_seq += ref_seq.substr(cur_ref_pos, (ref_seq.size() - 1) - cur_ref_pos);
                fasta_ss << cons_seq << '\n';
                if (dedup) // terminate after writing the first of the duplicates
                    break;
            }
        }

        fasta_ss.close();

        // now we can write duplicate sequences
        const std::string dups_fname = out_fname + ".dups";
        std::ofstream dups_ss(dups_fname, std::ios_base::out);
        if (!dups_ss)
        {
            std::cerr << "Can't open dups file for writing!\n";
            exit(2);
        }

        dups_ss << "vars,rep,names,num_seqs\n";
        for (auto& rvar : var_sets)
        {
            for(auto& var : rvar.first)
            {
                dups_ss << std::get<0>(var) << std::get<1>(var) << std::get<2>(var);
            }
            dups_ss << ',' << rvar.second[0] << ',';
            for (auto& r : rvar.second)
            {
                this->ir_it.first = this->id2ref.find(r);
                dups_ss << this->ir_it.first->second << ';';
            }
            dups_ss.seekp(-1, std::ios_base::end);
            dups_ss << ',' << rvar.second.size() << '\n';
        }

        dups_ss.close();
    }

private:

    typedef std::tuple<int,uint8_t,uint8_t> Var; // position, query bases, reference bases
    std::map<Var,int> vars;
    std::pair<std::map<Var,int>::iterator,bool> vit;

    std::map<std::string,int> ref2id;
    std::pair<std::map<std::string,int>::iterator,bool> ri_it;
    std::map<int,std::string> id2ref;
    std::pair<std::map<int,std::string>::iterator,bool> ir_it;

    std::map<int,std::vector<Var>> ref2vars; // for each reference holds all variants
    std::pair<std::map<int,std::vector<Var> >::iterator,bool> rv_it;

    std::map<int,std::array<int,5> > pos2var; // position to counters of A:0,C:1,G:2,T:3.AMB:4
    std::pair<std::map<int,std::array<int,5> >::iterator,bool> pv_it;

    int max_ref_id = 0;

	bool mask_used = false;
    std::vector<bool> mask; // mask loaded from VCF/BCF

};

int bolotie_cons(int argc, char **argv){
    enum Opt_CONS { REFERENCE   = 'x',
                    VARIANTS    = 'i',
                    OUTPUT      = 'o',
                    NUM_THREADS = 't',
                    DEDUP       = 'd',
                    MINVARNUM   = 'm',
                    KEEPINDEL   = 'k',
                    NOAMB       = 'n',
                    AMB2REF     = 'r',
					MASK        = 'v'};

    ArgParse args_cons("cons");

    args_cons.add_string(Opt_CONS::REFERENCE, "reference", "", "path to the reference sequence in FASTA format", true);
    args_cons.add_string(Opt_CONS::VARIANTS, "query", "", "path to the variants file generated by the aln method", true);
    args_cons.add_int(Opt_CONS::NUM_THREADS, "threads", 1, "number of threads to use", false);
    args_cons.add_string(Opt_CONS::OUTPUT, "output", "-", "output path and name", true);
    args_cons.add_flag(Opt_CONS::DEDUP, "dedup", "only write unique sequences", false);
    args_cons.add_int(Opt_CONS::MINVARNUM, "min_var_num", 1, "minimum number of sequences in which a variant is to be present", false);
    args_cons.add_flag(Opt_CONS::KEEPINDEL, "keepindel", "keep insertions and deletions", false);
    args_cons.add_flag(Opt_CONS::NOAMB, "noamb", "replaces all ambiguous characters with reference alleles", false);
    args_cons.add_int(Opt_CONS::AMB2REF, "amb2ref", 0, "if at most in this many sequences variants are present at a position where an ambiguous base is observed, then the reference allele is written", false);
	args_cons.add_string(Opt_CONS::MASK,"mask","","VCF/BCF file listing variants to be masked (hidden) in the consensus sequence.",false);

    args_cons.parse_args(argc, argv);

    // load reference
    std::ifstream ref_source;
    ref_source.open(args_cons.get_string(Opt_CONS::REFERENCE));
    if (!ref_source)
    {
        std::cerr << "Can't open the reference fasta file!\n";
        exit(2);
    }

    std::string ref_seqid, ref_seq, ref_line;
    ref_seq.reserve(29003);
    while (std::getline(ref_source, ref_line)) // parse fasta and add to sequence
    {
        if (ref_line[ref_line.size()-1] == '\r')
        {
            ref_line.pop_back();
        }
        if (ref_line[0] == '>')
        {
            ref_seqid = ref_line.substr(1,ref_line.size());
            ref_seq = "";
        }
        else
        {
            ref_seq += ref_line;
        }
    }

    // load variants
    Variants vars;

	// if available - load masking data
    if(args_cons.is_set(Opt_CONS::MASK))
	{
        std::cerr<<"loading masking variant set"<<std::endl;
        vars.load_vcf_mask(args_cons.get_string(Opt_CONS::MASK),ref_seqid);
    }

    std::ifstream var_source;
    var_source.open(args_cons.get_string(Opt_CONS::VARIANTS));
    if (!var_source)
    {
        std::cerr << "Can't open the variants file\n";
        exit(2);
    }

    std::string var_line, var_seqid, var_type, var_pos, var_query, var_ref;
    while (std::getline(var_source, var_line))
    {
        // now read columns
        std::istringstream ss(var_line);

        std::getline(ss, var_seqid, ',');
        std::getline(ss, var_type, ',');

        if (var_type != "S") // if not SNV
            continue;

        std::getline(ss, var_pos, ',');
        std::getline(ss, var_query, ',');
        std::getline(ss, var_ref, ',');

        if (!(args_cons.is_set(Opt_CONS::NOAMB) && !ACGT[var_query[0]]))
            vars.add_var(var_seqid, std::stoi(var_pos), var_query[0], var_ref[0]);
    }

    std::transform(ref_seq.begin(),ref_seq.end(),ref_seq.begin(),
                   [](unsigned char c) -> unsigned char { return std::toupper(c); });
    ref_source.close();

    // write consensus
    const std::string out_base_fname = args_cons.get_string(Opt_CONS::OUTPUT);
    vars.write_cons(ref_seq, out_base_fname, args_cons.get_int(Opt_CONS::MINVARNUM), args_cons.get_int(Opt_CONS::AMB2REF), args_cons.get_flag(Opt_CONS::DEDUP));

    return 0;
}
