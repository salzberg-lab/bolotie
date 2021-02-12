Bolotie: detecting recombinations in viruses using large data sets with high sequence similarity
================================================================================================

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://opensource.org/licenses/GPL-3.0
    :alt: GPLv3 License

.. contents::
   :local:
   :depth: 2

Introduction
^^^^^^^^^^^^

Bolotie is a new algorithm designed to conduct mutational analysis and to detect recombinant forms and other anomalies among a very large set of viral sequences.
methods implemented in Bolotie are designed such that novel sequences can also be analyzed efficiently without the need to rerun the entire protocol.

Ales Varabyou, Christopher Pockrandt, Steven L. Salzberg, Mihaela Pertea. **Rapid detection of inter-clade recombination in SARS-CoV-2 with Bolotie**. `BioRxiv`_, 2020.

Installation
^^^^^^^^^^^^

Building from source
""""""""""""""""""""

If you want to build it from source, we recommend cloning the git repository as shown below.

::

    $ git clone --recursive https://github.com/salzberg-lab/bolotie.git --recursive
    $ cmake -DCMAKE_BUILD_TYPE=Release .
    $ make -j4

If you are using a very old version of Git (< 1.6.5) the flag ``--recursive`` does not exist.
In this case you need to clone the submodule separately  (``git submodule update --init --recursive``).

**Requirements**

Operating System
  GNU/Linux

Architecture
  Intel/AMD platforms that support ``POPCNT``

Compiler
  GCC ≥ 4.9, Clang ≥ 3.8

Build system
  CMake ≥ 3.0

Language support
  C++14

Python3 requirements
::
    $ pip install -r requirements.txt

R dependencies:
  ape
  plyr

Getting started
^^^^^^^^^^^^^^^
A wrapper script `run.py` is provided to run all steps of the protocol.
This script, while convenient, is intended for replicability and testing and lacks some of the available features of Bolotie.
Due to additional external dependencies, tree-building and plotting is not included as part of the script but the tree plotting script
implemented in R is also provided and described below. Lastly, the clade assignment .tsv file requires clades IDs to be listed as unique integer IDs in incremental order.
The plotting function only supports upto 8 clades but Bolotie.

Arguments:
     -h, --help            show this help message and exit
     --input INPUT         File containing all assembled genomes in FASTA format
     --reference REFERENCE
                           File containing the reference genome in FASTA format
     --outdir OUTDIR       Output directory in which all output and temporary
                           data will be stored
     --keep_tmp            If enabled - temporary data will not be removed (default=False)
     --threads THREADS     Number of threads to use for parallelized sections of
                           the analysis. (default = 1)
     --freq-threshold FREQ_THRESHOLD
                           Variant frequency threshold. (default = 100)
     --trim-len TRIM_LEN   Number of bases at 3' and 5' ends of each sequence kept same as the reference for all sequences. (default = 200)
     --clusters CLUSTERS   File with sequence to cluster assignments in a tsv
                           format: <seqname> <cluster>
     --first-stage
                           Stage at which to start the run. Options: clean,align,consensus,index,path,parents,plot. (default = clean)
     --last-stage
                           Stage at which to end the run. Options: clean,align,consensus,index,path,parents,plot. (default = plot)
     --probability PROBABILITY
                           Probability of staying in the same state of the probability table. (default = 0.9999)

The script completes by creating several important files in the output directory:
 1. query.fasta - this file contains a copy of the input set of sequences with sequence names standardized in a way compatible with Bolotie, plotting scripts and any additional downstream methods such as IQ-TREE
 2. query.clus  - contains a copy of the clusters file standardized in the same way as the query.fasta
 3. aln.csv     - provides an extended cigar string for each sequences
 4. alb.vars.csv - individual variants detected from alignments
 5. cons4prob.fasta - file contains consensus sequences used to build the probability matrix. These sequences only contain the high-frequency variants.
 6. cons4parents.fasta - file contains consensus sequences used to search for parent sequences. These sequences contain all high and low-frequency variants detected by the aligner.
 7. probmat.probs - file contains the 3-dimensional probability table, where 1D is the positions along the reference genome; 2D corresponds to ACGT nucleotides and 3D contains individual probabilities for each cluster
 8. probmat.totals - contains a 2D table where rows are positions along the reference genome and columns are bases. Values correspond to the total number of times each base is seen at that position
 9. probmat.counts - file contains the 3-dimensional matrix where 1D is the positions along the reference genome; 2D corresponds to ACGT nucleotides and 3D contains individual counts for each cluster
 10. paths - recombinant and otherwise anomalous sequences detected by Bolotie
 11. parents - for each sequence in the "paths" file this file lists most likely parental sequences for each segment of a recombinant.

A plotting function `plot_utree.R` is used to produce a plot of an unrooted tree colored by clades. The script also labels any anomalous sequences detected by Bolotie.
To run the script use the following command:
::
    $ plot_utree.R treefile.nwk clusterfile.tsv palette.tsv output.png paths.
Treefile.nwk needs to be computed separately, but all other inputs can be found in the output directory of the `run.py`.


Data
^^^^

If you want to evaluate specific sequences for recombination events, you can download the probability table which we computed on 87,695 complete high-coverage genomes obtained from GISAID (September 2). This significantly reduces the running time from several hours to seconds.

`ftp://ftp.ccb.jhu.edu/pub/data/bolotie_sars_cov_2/ <ftp://ftp.ccb.jhu.edu/pub/data/bolotie_sars_cov_2/>`_
