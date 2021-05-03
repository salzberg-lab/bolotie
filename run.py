#!/usr/bin/env python3

import os
import sys
import time
import glob
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

# Common Declarations and Functions

all_stages = ["clean", "align", "consensus", "index", "path", "parents", "plot"]

acgt = {"A": 0,
        "C": 1,
        "G": 2,
        "T": 3}
tgca = ["A", "C", "G", "T"]


# https://stackoverflow.com/questions/5560248/programmatically-lighten-or-darken-a-hex-color-or-rgb-and-blend-colors
def shade_color(color, percent):
    color = color[1:]
    rgb = tuple(int(color[i:i + 2], 16) for i in (0, 2, 4))
    rgb = [c * (100 + percent) / 100 for c in rgb]
    rgb = [int(c) if c < 255 else 255 for c in rgb]
    return '#%02x%02x%02x' % tuple(rgb)


# load the probability table
def load_mat(prob_fname):
    # scan the file to get dimensions
    pos_count = 0
    clu_count = 0
    with open(prob_fname, "r") as inFP:
        for line in inFP.readlines():  # each line is a position
            pos_count += 1
            for bid, col in enumerate(line.split("\t")):  # each column is an ACGT
                clu_count = len(col.split(","))

    mat = np.zeros(shape=(pos_count, len(acgt), clu_count))

    pos = 0
    with open(prob_fname, "r") as inFP:
        for line in inFP.readlines():  # each line is a position
            for bid, col in enumerate(line.split("\t")):  # each column is an ACGT
                for cluid, prob in enumerate(col.split(",")):  # each value is a probability of a cluster
                    mat[pos, bid, cluid] = float(prob)
            pos += 1

    return mat


# load the probability table
def load_totals(counts_fname):
    # scan the file to get dimensions
    pos_count = 0
    with open(counts_fname, "r") as inFP:
        for _ in inFP.readlines():  # each line is a position
            pos_count += 1

    mat = np.zeros(shape=(pos_count, len(acgt)))

    pos = 0
    with open(counts_fname, "r") as inFP:
        for line in inFP.readlines():  # each line is a position
            for bid, total in enumerate(line.split("\t")):  # each column is an ACGT
                mat[pos, bid] = int(total)
            pos += 1

    return mat


# now load all sequences with best paths
def load_paths(path_fname, seqs_fname):
    paths_dict = dict()
    with open(path_fname, "r") as inFP:
        for line in inFP.readlines():
            name, prob, signature = line.strip().split("\t")
            assert name not in paths_dict, "duplicate names"
            paths_dict[name] = [signature]

    # now load sequences
    with open(seqs_fname, "r") as inFP:
        write_seq = False
        cur_name = ""
        for line in inFP.readlines():
            if line[0] == ">":
                if write_seq:
                    paths_dict[cur_name].append(cur_seq)
                cur_name = line[1:].strip()
                cur_seq = ""
                if cur_name in paths_dict:
                    write_seq = True
                else:
                    write_seq = False
            else:
                if write_seq:
                    cur_seq += line.strip()

    return paths_dict


def get_probs(seq, mat):
    probs = []
    for i, nt in enumerate(seq):
        probs.append(list())
        for cluid in range(mat.shape[2]):
            if nt.upper() not in "ACGT":
                probs[-1].append(1.0 / float(mat.shape[2]))
            else:
                probs[-1].append(mat[i, acgt[nt.upper()], cluid])
    return probs


def timer(start, end):
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)


def get_tmp_id(base_tmp_dir):
    id_val = 0
    tmp = base_tmp_dir
    while True:
        tmp = tmp + str(id_val)
        if not os.path.exists(tmp):
            return tmp
        else:
            id_val += 1


def get_stages(args):
    global all_stages

    run_stages = []

    found_first = False
    for stage in all_stages:
        if args.clusters is None and stage == "index":
            sys.stderr.write(
                "can not perform indexing, path finding, parents finding and plotting without provided clusters. You can use consensus calls to generate a tree and clusters\n")
            return run_stages

        if args.first_stage is None or args.first_stage == stage:
            found_first = True
        if found_first:
            run_stages.append(stage)
        if args.last_stage == stage:
            return run_stages

    return run_stages


def load_vars(vars_fname, clu_fname):
    vars_df = pd.read_csv(vars_fname, names=["seqid", "type", "pos", "query", "ref"])
    vars_df["pos"] = vars_df["pos"].astype(int)
    vars_df = vars_df[vars_df["type"] == "S"].reset_index(drop=True)
    vars_df["query"] = vars_df["query"].str.upper()
    vars_df = vars_df[vars_df["query"].isin(["A", "C", "G", "T"])].reset_index(drop=True)
    vars_df.drop("type", inplace=True, axis=1)
    vars_df["seqid"] = vars_df["seqid"].str.replace(" ", "_")

    clu_df = pd.read_csv(clu_fname, sep="\t", names=["seqid", "clu"])
    # merge
    vars_df = vars_df.merge(clu_df, how="left", on="seqid", indicator=True)
    vars_df.drop("_merge", inplace=True, axis=1)

    clu_grp_vars = dict()

    for clu, grp in vars_df[["clu", "pos", "query", "ref"]].groupby(by="clu"):
        tmp = pd.DataFrame(grp.groupby(by=["pos", "query"])["ref"].count()).reset_index()
        tmp.columns = ["pos", "var", "count"]
        tmp = tmp[tmp["count"] > 10].reset_index(drop=True)
        tmp.sort_values(by="pos", inplace=True, ascending=True)

        xs = tmp[["pos"]].drop_duplicates()["pos"].to_list()
        nt_var_dict = dict()
        for nt in "ACGT":  # for each base
            # build a list of counts - 0 must be placed if the nucleotide does not appear at this position
            tmp2 = tmp[["pos"]].drop_duplicates()
            tmp2 = tmp2.merge(tmp[tmp["var"] == nt].reset_index(drop=True), on="pos", how="left")
            tmp2.drop("var", axis=1, inplace=True)
            tmp2.fillna(0, inplace=True)
            tmp2.sort_values(by="pos", ascending=True, inplace=True)
            nt_var_dict[nt] = tmp2["count"].to_list()
        clu_grp_vars[clu] = [xs, nt_var_dict]

    return vars_df, clu_grp_vars


def plot_recomb(axes, seq, mat, counts, totals, clus_plot, clus_scatter, palette, thresh=0.1):
    probs = get_probs(seq, mat)
    probs = pd.DataFrame(probs)

    # now need to pull values for the bases for these clusters
    cur_xs = {key: list() for key in clus_plot}  # position
    cur_counts = {key: list() for key in clus_plot}  # fraction of sequences with this variant
    cur_probs = {key: list() for key in clus_plot}  # probabilities for each cluster
    base_prob = 1.0 / float(mat.shape[2])
    for i, nt in enumerate(seq):
        for clu in clus_plot:
            if nt.upper() not in "ACGT":
                prob = base_prob
            else:
                prob = mat[i][acgt[nt]][clu]
            if abs(prob - base_prob) >= 0.09:
                cur_probs[clu].append(prob)
                cur_xs[clu].append(i)
                cur_clu_count = int(counts[i][acgt[nt]][clu])
                cur_counts[clu].append(str(cur_clu_count))

    for rec_clu in clus_plot:
        axes.plot(probs[rec_clu], linewidth=10, color=palette.get(rec_clu, '#333333'))

    base_frac = 1.0 / float(len(clus_scatter))
    markers = ["o", "s", "v", "X", "^", "<", ">", "8", "p", "P", "*", "h", "H", "D", "d"]
    for i, rec_clu in enumerate(clus_scatter):
        small_size = (5 * len(clus_scatter)) - (i * 5)
        large_size = (150 * len(clus_scatter)) - (i * 150)
        mark_size = [large_size if x > thresh else small_size for x in list(abs(probs[rec_clu] - base_frac))]
        axes.scatter(list(probs.index), probs[rec_clu], s=mark_size, marker=markers[i],
                     color=palette.get(rec_clu, '#333333'), label=str(rec_clu))

        if rec_clu in clus_plot:
            res = pd.concat([pd.DataFrame(cur_xs[rec_clu]),
                             pd.DataFrame(cur_probs[rec_clu]),
                             pd.DataFrame(cur_counts[rec_clu])], axis=1)
            res.columns = ["xs", "prob", "counts"]

            idx = 0
            counts_tmp = []
            is_tmp = []
            js_tmp = []

            j = list(res["prob"])[0]
            i = list(res["xs"])[0]
            counts_tmp.append(int(list(res["counts"])[0]))
            is_tmp.append(i)
            js_tmp.append(j)
            prev_j = j
            prev_i = i
            idx += 1

            for i, j in zip(list(res["xs"])[1:], list(res["prob"])[1:]):  # TODO: "i" is already declared
                if max([j, prev_j]) - min([j, prev_j]) <= 0.05 and max([i, prev_i]) - min([i, prev_i]) <= 500:
                    counts_tmp.append(int(list(res["counts"])[idx]))
                    is_tmp.append(i)
                    js_tmp.append(j)
                    prev_j = j
                    prev_i = i
                else:
                    txt = str(int(np.mean(counts_tmp)))
                    if len(counts_tmp) > 1:
                        txt = txt + "(" + str(len(counts_tmp)) + ")"
                    axes.annotate(txt, xy=(np.mean(is_tmp), np.mean(js_tmp)), fontsize=30, rotation=0)

                    counts_tmp = [int(list(res["counts"])[idx])]
                    is_tmp = [i]
                    js_tmp = [j]
                    prev_j = j
                    prev_i = i
                idx += 1

            txt = str(int(np.mean(counts_tmp)))
            if len(counts_tmp) > 1:
                txt = txt + "(" + str(len(counts_tmp)) + ")"
            axes.annotate(txt, xy=(np.mean(is_tmp), np.mean(js_tmp)), fontsize=30, rotation=0)

    axes.legend(loc=7)


def load_parents(parents_fname, mat):
    parents_dict = dict()
    with open(parents_fname, "r") as inFP:
        for line in inFP.readlines():
            line_cols = line.strip().split(",")
            rec = line_cols[0]
            par1c = int(line_cols[1])
            par1 = line_cols[4]
            par2c = int(line_cols[6])
            par2 = line_cols[9]
            parents_dict[rec] = [0 for _ in range(mat.shape[2])]
            parents_dict[rec][par1c] = par1
            parents_dict[rec][par2c] = par2

    return parents_dict


# Main Methods

def clean(args):
    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nCleaning sequence names in all inputs\n")

    excluded_seqids = set() # list of sequence identifiers to be excluded based on the parameters
    included_seqids = set() # list of passable identifiers

    # read cluster data and find seqids to remove based on missing data
    with open(args.clusters,"r") as inFP:
        for line in inFP.readlines():
            line = line.replace("/", "_").replace("|", "_").replace(" ", "_").replace(",","_")
            cols = line.strip().split("\t")
            cluid = cols[-1].split(".")[0]
            if not cluid.isdigit():
                excluded_seqids.add(cols[0])

    out_fasta = od + "query.fasta"
    with open(out_fasta, "w+") as outFP:
        with open(args.input, "r") as inFP:
            seqid = None
            seq = ""
            for line in inFP.readlines():
                if line[0]==">":
                    if seqid is not None: # write out previous entry
                        num_acgt = seq.count("A")+seq.count("C")+seq.count("G")+seq.count("T")
                        perc_non_acgt = (float(len(seq)-num_acgt)/float(len(seq)))*100
                        if perc_non_acgt <= args.perc_ambiguous and seqid not in excluded_seqids:
                            included_seqids.add(seqid)
                            outFP.write(">"+seqid+"\n")
                            outFP.write(seq+"\n")
                        else:
                            excluded_seqids.add(seqid)

                    line = line.replace("/", "_").replace("|", "_").replace(" ", "_").replace(",","_")
                    seqid = line.strip()[1:]
                    seq = ""
                else:
                    seq += line.strip().upper()

            if seqid is not None: # write out last entry
                num_acgt = seq.count("A")+seq.count("C")+seq.count("G")+seq.count("T")
                perc_non_acgt = (float(len(seq)-num_acgt)/float(len(seq)))*100
                if perc_non_acgt <= args.perc_ambiguous:
                    included_seqids.add(seqid)
                    outFP.write(">"+seqid+"\n")
                    outFP.write(seq+"\n")
                else:
                    excluded_seqids.add(seqid)

    if args.clusters is not None:
        out_clusters = od + "query.clus"
        with open(out_clusters, "w+") as outFP:
            with open(args.clusters, "r") as inFP:
                for line in inFP.readlines():
                    line = line.replace("/", "_").replace("|", "_").replace(" ", "_").replace(",","_")
                    seqid = line.strip().split("\t")[0]
                    cluid = line.strip().split("\t")[-1].split(".")[0]
                    if seqid in included_seqids:
                        assert seqid not in excluded_seqids,"included and excluded"
                        outFP.write(seqid+"\t"+cluid+"\n")
                    else:
                        assert seqid in excluded_seqids,"not included and not excluded: "+str(seqid)

    sys.stderr.write("Done\n")
    return


def align(args, bolotie_path):
    od = args.outdir.rstrip("/") + "/"
    clean_input = od + "query.fasta"

    sys.stderr.write("\n==================\nAligning sequences\n")
    start_time = time.time()

    aln_cmd = [bolotie_path, "aln",
               "-x", args.reference,
               "-i", clean_input,
               "-o", od + "aln",
               "-t", str(args.threads),
               "-l", str(args.trim_len)]
    sys.stderr.write("command: " + " ".join(aln_cmd) + "\n")
    ret = subprocess.call(aln_cmd)
    if ret != 0:
        sys.stderr.write("Non-zero exit code in alignment")
        exit(1)

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done aligning sequences in " + exec_time + "\n")

    return


def consensus(args, bolotie_path):
    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nBuilding deduplicated consensus sequences\n")

    start_time = time.time()
    cons4tree_cmd = [bolotie_path, "cons",
                     "-x", args.reference,
                     "-i", od + "aln.vars.csv",
                     "-o", od + "cons4tree",
                     "-m", str(args.freq_threshold),
                     "-d"]
    sys.stderr.write("command: " + " ".join(cons4tree_cmd) + "\n")
    ret = subprocess.call(cons4tree_cmd)
    if ret != 0:
        sys.stderr.write("Non-zero exit code in consensus 4 tree")
        exit(1)

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done in " + exec_time + "\n")

    #########################################################
    #########################################################
    #########################################################

    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nBuilding consensus sequences with frequency cut\n")

    start_time = time.time()
    cons4prob_cmd = [bolotie_path, "cons",
                     "-x", args.reference,
                     "-i", od + "aln.vars.csv",
                     "-o", od + "cons4prob",
                     "-m", str(args.freq_threshold)]
    sys.stderr.write("command: " + " ".join(cons4prob_cmd) + "\n")
    ret = subprocess.call(cons4prob_cmd)
    if ret != 0:
        sys.stderr.write("Non-zero exit code in consensus for probability")
        exit(1)

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done in " + exec_time + "\n")

    #########################################################
    #########################################################
    #########################################################

    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nBuilding raw consensus sequences\n")

    start_time = time.time()
    cons4parent_cmd = [bolotie_path, "cons",
                       "-x", args.reference,
                       "-i", od + "aln.vars.csv",
                       "-o", od + "cons4parent"]
    sys.stderr.write("command: " + " ".join(cons4parent_cmd) + "\n")
    ret = subprocess.call(cons4parent_cmd)
    if ret != 0:
        sys.stderr.write("Non-zero exit code in consensus for parents")
        exit(1)

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done in " + exec_time + "\n")

    return


def index(args, bolotie_path):
    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nBuilding the conditional probability table\n")

    start_time = time.time()
    build_cmd = [bolotie_path, "build",
                 "-x", od + "probmat",
                 "-i", od + "cons4prob.fa",
                 "-c", od + "query.clus"]
    sys.stderr.write("command: " + " ".join(build_cmd) + "\n")
    ret = subprocess.call(build_cmd)
    if ret != 0:
        sys.stderr.write("Non-zero exit code in indexing")
        exit(1)

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done in " + exec_time + "\n")

    return


def paths(args, bolotie_path):
    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nSearching for recombinant paths\n")

    start_time = time.time()
    path_cmd = [bolotie_path, "find",
                "-x", od + "probmat.probs",
                "-i", od + "cons4prob.fa",
                "-p", str(args.probability),
                "-t", str(args.threads),
                "-o", od + "paths",
                "-r"]
    sys.stderr.write("command: " + " ".join(path_cmd) + "\n")
    ret = subprocess.call(path_cmd)
    if ret != 0:
        sys.stderr.write("Non-zero exit code in path finding")
        exit(1)

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done in " + exec_time + "\n")

    return


def parents(args, bolotie_path):
    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nSearching for parents\n")

    start_time = time.time()
    parents_cmd = [bolotie_path, "parents",
                   "-f", od + "cons4parent.fa",
                   "-c", od + "query.clus",
                   "-p", od + "paths",
                   "-o", od + "parents"]
    sys.stderr.write("command: " + " ".join(parents_cmd) + "\n")
    ret = subprocess.call(parents_cmd)
    if ret != 0:
        sys.stderr.write("Non-zero exit code in parents finding")
        exit(1)

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done in " + exec_time + "\n")

    return


def plot(args, figure_dir):
    od = args.outdir.rstrip("/") + "/"

    sys.stderr.write("\n==================\nPlotting figures\n")
    start_time = time.time()

    curpal = sns.color_palette("pastel")
    cold = {i: curpal.as_hex()[i] for i in range(len(curpal))}
    curpal_dark = sns.color_palette("colorblind")
    cold_dark = {i: curpal_dark.as_hex()[i] for i in range(len(curpal_dark))}
    light_shades = []
    for _, c in cold.items():
        light_shades.append(shade_color(c, 40))

    # save the colorblind color_palette
    cold_dark_df = pd.DataFrame.from_dict(cold_dark, orient="index").reset_index()
    cold_dark_df.columns = ["clu_id", "color"]
    cold_dark_df.to_csv(figure_dir + "palette.csv", index=False, header=False, sep="\t")

    mat = load_mat(od + "probmat.probs")
    counts = load_mat(od + "probmat.counts")
    totals = load_totals(od + "probmat.totals")
    paths = load_paths(od + "paths", od + "cons4prob.fa")

    vars_df, clu_grp_vars = load_vars(od + "aln.vars.csv", od + "query.clus")
    load_parents(od + "parents", mat)

    params = {'font.size': 30}
    plt.rcParams.update(params)

    for refname in paths:
        signature, seq = paths[refname]
        lens = [int(x.split(":")[1]) for x in signature.split(";")]

        rec_clus = [int(x.split(":")[0]) for x in signature.split(";")]

        fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(40, 8), sharex=True,
                                 gridspec_kw={'hspace': 0, 'height_ratios': [6, 1, 1, 1]})
        axes[0].set_ylabel("P(Clade | Base)")
        plt.setp(axes[1].get_yticklabels(), visible=False)
        axes[1].tick_params(axis='both', which='both', length=0)
        plt.setp(axes[2].get_yticklabels(), visible=False)
        axes[2].tick_params(axis='both', which='both', length=0)
        plt.setp(axes[3].get_yticklabels(), visible=False)
        axes[3].tick_params(axis='both', which='both', length=0)

        axes[1].set_ylabel("Clade " + str(rec_clus[0]), rotation=0, va="center")
        axes[2].set_ylabel("Recombinant", rotation=0, va="center")
        axes[3].set_ylabel("Clade " + str(rec_clus[1]), rotation=0, va="center")

        plot_recomb(axes[0], seq, mat, counts, totals, rec_clus, list(range(mat.shape[2])), cold_dark, thresh=0.1)

        # plot recomb variants
        rec_vars = vars_df[vars_df["seqid"] == refname].sort_values(by="pos")
        for nt in "ACGT":
            xs = rec_vars[rec_vars["query"] == nt]["pos"].to_list()
            heights = [1 for _ in range(len(xs))]
            axes[2].bar(x=xs, height=heights, width=200, color="black")

        # 1st should be consensus SNPS for the 1st cluster
        # 2nd should be bases of the 1st parent
        axes[1].set_facecolor(cold_dark.get(rec_clus[0], '#333333') + str("20"))
        xs = clu_grp_vars[rec_clus[0]][0]
        prevs = [0 for _ in range(len(xs))]
        for nt in "ACGT":
            cur_vals = clu_grp_vars[rec_clus[0]][1][nt]
            axes[1].bar(x=xs, height=cur_vals, width=200, bottom=prevs, color="black")
            for i in range(len(prevs)):
                prevs[i] += cur_vals[i]

        axes[3].set_facecolor(cold_dark.get(rec_clus[1], '#333333') + str("20"))
        xs = clu_grp_vars[rec_clus[1]][0]
        prevs = [0 for _ in range(len(xs))]
        for nt in "ACGT":
            cur_vals = clu_grp_vars[rec_clus[1]][1][nt]
            axes[3].bar(x=xs, height=cur_vals, width=200, bottom=prevs, color="black")
            for i in range(len(prevs)):
                prevs[i] += cur_vals[i]

        plt.suptitle(refname, fontsize=40)

        if args.vector:
            plt.savefig(figure_dir + refname + ".svg")
        else:
            plt.savefig(figure_dir + refname + ".png")

        plt.clf()
        plt.close('all')

    end_time = time.time()
    exec_time = timer(start_time, end_time)
    sys.stderr.write("Done plotting in " + exec_time + "\n")

    return


def wrapper(args):
    total_start_time = time.time()

    run_stages = get_stages(args)
    sys.stderr.write("Running stages: " + ", ".join(run_stages) + "\n")

    cur_path = os.path.dirname(os.path.abspath(__file__))
    bolotie_path = cur_path + "/bolotie"
    if not os.path.exists(bolotie_path):
        found = False
        for tmp_path in glob.glob(cur_path + "/*/*"):
            tmp_fname = tmp_path.split("/")[-1]
            if tmp_fname == "bolotie":
                bolotie_path = tmp_path
                found = True
                break
        if not found:
            sys.stderr.write(
                "Did not find bolotie executable path. Please make sure bolotie executable is located within the same directory as the run.py script or whithin a single directory level down.")
            sys.exit(2)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    tmp_dir = get_tmp_id(args.outdir.rstrip("/") + "/tmp")
    assert not os.path.exists(tmp_dir), "new tmp directory already exists"
    os.makedirs(tmp_dir)

    figure_dir = args.outdir.rstrip("/") + "/figures/"
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    if "clean" in run_stages:
        clean(args)
    if "align" in run_stages:
        align(args, bolotie_path)
    if "consensus" in run_stages:
        consensus(args, bolotie_path)
    if "index" in run_stages:
        index(args, bolotie_path)
    if "path" in run_stages:
        paths(args, bolotie_path)
    if "parents" in run_stages:
        parents(args, bolotie_path)
    if "plot" in run_stages:
        plot(args, figure_dir)

    if not args.keep_tmp:
        shutil.rmtree(tmp_dir)

    total_end_time = time.time()
    exec_time = timer(total_start_time, total_end_time)
    sys.stderr.write("Total runtime " + exec_time + "\n")

    return


def main(args):
    global all_stages
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('--input',
                        required=True,
                        type=str,
                        help="File containing all assembled genomes in FASTA format")
    parser.add_argument('--reference',
                        required=True,
                        type=str,
                        help="File containing the reference genome in FASTA format")
    parser.add_argument('--outdir',
                        required=True,
                        type=str,
                        help="Output directory in which all output and temporary data will be stored")
    parser.add_argument("--keep_tmp",
                        action='store_true',
                        help="If enabled - temporary data will not be removed")
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
                        help="Number of bases at 3' and 5' ends of each sequence kept same as the reference for all sequences.")
    parser.add_argument("--clusters",
                        required=False,
                        type=str,
                        help="File with sequence to cluster assignments in a tsv format: <seqname>\t<cluster>")
    parser.add_argument("--first-stage",
                        required=False,
                        type=str,
                        default="clean",
                        choices=all_stages,
                        help="Stage at which to start the run")
    parser.add_argument("--last-stage",
                        required=False,
                        type=str,
                        default="plot",
                        choices=all_stages,
                        help="Stage at which to end the run")
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
    parser.add_argument("--vector",
                        required=False,
                        action="store_true",
                        help="If enabled - all plots will be saved in SVG format. PNG images are stored by default")

    parser.set_defaults(func=wrapper)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
