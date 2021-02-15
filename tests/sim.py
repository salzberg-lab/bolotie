#!/usr/bin/env python3

import pandas as pd
import subprocess
import random
import sys
import os
import argparse


# assessment of FPR on simulated data

# 1. get consensus for each clade as a founder population
# 2. for each clade iteratively create generations of sequences

def add_vars(var_list, seq,
             refseq):  # accepts dictionary to store data and the sequence to process and the reference for comparison
    assert len(seq) == len(refseq), "sequence and reference of different lengths"
    for i in range(len(seq)):
        if seq[i] in "ACGT":
            var_list[i][seq[i]] += 1
        else:
            var_list[i][refseq[i]] += 1


# only use bases which pass a minimum threshold (eg. 10%)

def get_random_nt(nts, min_freq):
    # find the base with minimum count above the random
    keys = []
    vals = []
    total = sum(nts.values())
    for k, v in nts.items():
        if v / total >= min_freq:
            keys.append(k)
            vals.append(v / total)
    return random.choices(population=keys, weights=vals, k=1)[0]


# now generate consensus and for each consensus sequence we can generate a collection of genomes for analysis

def mutate_write(outfp, cid, sid, seq, mut_rate):
    outfp.write(">" + str(cid) + "_" + str(sid) + "\n")
    sim_seq = ""
    # need the number of generations
    for nt in seq:
        p = random.uniform(0, 1)
        if p <= mut_rate:
            nt = random.choice("ACGT".replace(nt, ""))
        sim_seq += nt
    outfp.write(sim_seq + "\n")


def run(args):
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    run_dir = args.run_dir
    ref_fname = args.reference

    num_genomes_per_clade = args.num_genome_per_clade
    mut_rate = args.mut_rate
    min_allele_freq = args.min_allele_freq  # 10 percent minimum frequency of an allele

    # get consensus for each clade
    seqid2clu_fname = run_dir + "query.clus"
    seqid2clu = dict()
    clades = dict()

    with open(seqid2clu_fname, "r") as inFP:
        for line in inFP.readlines():
            seqid, clu = line.strip().split("\t")
            assert seqid not in seqid2clu, "duplicate seqids"
            seqid2clu[seqid] = clu
            if clu not in clades:
                clades[clu] = list()
            clades[clu].append(seqid)

    print("loaded " + str(len(seqid2clu)) + " sequence IDs with clades")
    print("loaded " + str(len(clades)) + " clades")

    # load the reference sequence
    num_refs = 0
    refseq = ""
    with open(ref_fname, "r") as inFP:
        for line in inFP.readlines():
            if line[0] == ">":
                num_refs += 1
                continue
            else:
                refseq += line.strip().upper()
    assert num_refs == 1, "incorrect number of sequences in the reference"

    # trim the last base of the reference because that's what bolotie does
    refseq = refseq[:-1]

    # compute the consensus for each clade
    gisaid_fasta = run_dir + "cons4parent.fa"  # using cons4parent as sequences are aligned, normalized and the same
    # length while preserving all meaningful variants

    clu_vars = dict()

    with open(gisaid_fasta, "r") as inFP:
        seqid = None
        cluid = None
        seq = ""
        for line in inFP.readlines():
            if line[0] == ">":
                if seqid is not None and cluid is not None:
                    add_vars(clu_vars[cluid], seq, refseq)
                seqid = line[1:].strip()
                cluid = seqid2clu[seqid]
                seq = ""
                if cluid not in clu_vars:
                    clu_vars[cluid] = [dict({"A": 0, "C": 0, "G": 0, "T": 0}) for _ in range(
                        len(refseq))]  # value is a dictionary of nucleotide counts for each base in the sequence
            else:
                seq += line.strip().upper()
        # process last sequence
        add_vars(clu_vars[cluid], seq, refseq)

    # compute not only the consensus, but also at least several sublineages (50% consensus,40% consensus,
    # 30% consensus,20% consensus, 10% consensus)

    # compute frequency of alleles

    var_df = pd.DataFrame(clu_vars["1"])
    var_df["total"] = var_df["A"] + var_df["C"] + var_df["G"] + var_df["T"]
    var_df["A_frac"] = var_df["A"] / var_df["total"]
    var_df["C_frac"] = var_df["C"] / var_df["total"]
    var_df["G_frac"] = var_df["G"] / var_df["total"]
    var_df["T_frac"] = var_df["T"] / var_df["total"]

    var_df.head()

    sim_seqs_fname = out_dir + "sim_seqs.fa"
    with open(sim_seqs_fname, "w+") as outFP:
        clade_seq = None
        for cladeid, nts in clu_vars.items():
            for sid in range(num_genomes_per_clade):
                if clade_seq is not None:
                    mutate_write(outFP, cladeid, sid, clade_seq, mut_rate)
                clade_seq = ""
                for nt in nts:
                    cons_nt = get_random_nt(nt, min_allele_freq)  # get base (weighted random)
                    clade_seq += cons_nt
    #                 if not refseq[len(clade_seq)-1]==clade_seq[len(clade_seq)-1]:
    #                     print(cladeid,sid,len(clade_seq)-1,refseq[len(clade_seq)-1],clade_seq[len(clade_seq)-1])

    # write out cluster information
    sim_seqs_fname = out_dir + "sim_seqs.fa"
    sim_clus_fname = out_dir + "sim_seqs.clus"

    with open(sim_clus_fname, "w+") as outFP:
        with open(sim_seqs_fname, "r") as inFP:
            for line in inFP.readlines():
                if line[0] == ">":
                    seqid = line[1:].strip()
                    cluid = seqid.split("_")[0]
                    outFP.write(seqid + "\t" + cluid + "\n")

    # run bolotie
    blt_cmd = [args.bolotie_run_py,
               "--input", out_dir + "sim_seqs.fa",
               "--reference", args.reference,
               "--clusters", out_dir + "sim_seqs.clus",
               "--outdir", out_dir,
               "--threads", args.threads,
               "--freq-threshold", args.freq_threshold,
               "--trim-len", args.trim_len,
               "--probability", args.probability,
               "--perc_ambiguous", args.perc_ambiguous]
    print(" ".join(blt_cmd))
    subprocess.call(blt_cmd)


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('--run_dir',
                        required=True,
                        type=str,
                        help="path to the directory in which bolotie was previously run and all outputs are located "
                             "and named according to the run.py wrapper script")
    parser.add_argument('--reference',
                        required=True,
                        type=str,
                        help="File containing the reference genome in FASTA format")
    parser.add_argument('--outdir',
                        required=True,
                        type=str,
                        help="File containing the reference genome in FASTA format")
    parser.add_argument('--bolotie_run_py',
                        required=True,
                        type=str,
                        help="full path to the bolotie run.py script")
    parser.add_argument("--threads",
                        required=False,
                        type=int,
                        default=1,
                        help="Number of threads to use for parallelized sections of the analysis")
    parser.add_argument("--freq-threshold",
                        required=False,
                        type=int,
                        default=100,
                        help="Variant frequency threshold")
    parser.add_argument("--trim-len",
                        required=False,
                        type=int,
                        default=200,
                        help="Number of bases at 3' and 5' ends of each sequence kept same as the reference for all "
                             "sequences.")
    parser.add_argument("--probability",
                        type=float,
                        required=False,
                        default=0.9999,
                        help="Probability of staying in the same state (i.e., clade) of the probability table")
    parser.add_argument("--perc_ambiguous",
                        type=float,
                        required=False,
                        default=1.0,
                        help="percent of non-ACGT characters allowed for sequences to be analyzed")
    parser.add_argument("--min-allele-freq",
                        type=float,
                        required=False,
                        default=0.1,
                        help="minimum allele frequency in the input data to consider for the simulation")
    parser.add_argument("--mut-rate",
                        type=float,
                        required=False,
                        default=0.0006,
                        help="mutation rate (number of substitutions per site per year")
    parser.add_argument("--num-genomes-per-clade",
                        type=int,
                        required=False,
                        default=25000,
                        help="number of genomes to simulate for each clade")

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
