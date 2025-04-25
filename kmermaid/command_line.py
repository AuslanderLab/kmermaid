####################################################
#                        env                       #  
####################################################

import pandas as pd
import numpy as np
import pickle
import sys
import os
from os.path import isdir,isfile
import timeit
from .kmermaid import proc_classify_fasta,proc_classify_fastq,validate_fasta,validate_fastq
import argparse
from pathlib import Path
import kmermaid



####################################################
#                     constants                    #
####################################################

DESCRIPTION = "kmermaid: Ultrafast functional annotations of shotgut metagenomic sequencing reads" \
              ""
mpath = kmermaid.__file__
PWD = mpath[:-11]


SEGMENT_LENGTH=1000000000000
MIN_LEN=99 ##Minimun sequence lentgh
K = 5


gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

basepairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}



####################################################
#               Terminal functions                 #
####################################################

def kmermaid_predict():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type = str, help="Input fasta or fastq", required=True,
    )
    parser.add_argument(
        "--output", type = str, default = os.path.join('outputs', 'kmermaid_out'), help="Output file name"
    )
    parser.add_argument(
        "--cluster_reps", type = str, default = os.path.join('db', 'cluster_names.pkl'), help="Path to cluster names (pkl file) **modifications of this argument is not generally supported and at user responsibility**"
    )
    parser.add_argument(
        "--trained_model", type = str, default = os.path.join('db', 'kmer_model.pkl'), help="Path to cluster (pkl file) **modifications of this argument is not generally supported and at user responsibility**"
    )
    parser.add_argument(
        "--append_path", action='store_true', help="Use this flag with slurm"
    )
    parser.add_argument(
        "--check_format", action = 'store_true', help="Validate format of input file. Turning off this option may lead to incorrect results if the fasta or fastq is not properly formatted"
    )
    args = parser.parse_args()


    input_path = args.input
    output_path = args.output
    reps_path = PWD+args.cluster_reps
    model_path = PWD+args.trained_model
    app_path = args.append_path
    check_fmt = args.check_format

    if input_path is None:
        print('Error: Syntax: mikclust COMMAND [OPTIONS]. To print help message: mikclust -h')
        sys.exit(1)

    if not isfile(input_path):
        print('kmermaid: error: the following arguments are not found: '+input_path)
        sys.exit(1)

    if not isdir(os.path.dirname(os.path.abspath(output_path))):
        print('kmermaid: warning: the following arguments are not found:' + output_path + ' output will be generated in outputs/kmermaid_out')
        if not isdir(os.path.join('outputs', 'miiklust_out')):
            os.makedirs(os.path.join('outputs', 'kmermaid_out'), exist_ok = True)
        output_path=os.path.join('outputs', 'kmermaid_out', 'kmermaid_out')

    if app_path:
        sys.path.append(os.getcwd())


    ##load things, should be moved somewhere else:
    print("Loading databases", file=sys.stderr)

    with open(reps_path, 'rb') as f:
        names_dict = pickle.load(f)

    ###get the trained kmer dictionaries

    with open(model_path, 'rb') as fp:
        mod = pickle.load(fp)

    print("Predicting ", file=sys.stderr)

    with open(input_path) as fitmp:
        line0=fitmp.readline()

    start_time = timeit.default_timer()
    if line0.startswith('>'):
        if check_fmt:
            validate_fasta(open(input_path))
        else: 
            print("Skipping input validation")
        proc_classify_fasta(open(input_path), output_path + ".tsv", nmd = names_dict, segment_lengths = SEGMENT_LENGTH, minlen = MIN_LEN, gcode = gencode, bpairs = basepairs, K = K, dc = mod, check_fmt = check_fmt)
    elif line0.startswith('@'):
        if check_fmt:
            validate_fastq(open(input_path))
        else:
            print("Skipping input validatation")
        proc_classify_fastq(open(input_path), output_path + ".tsv", nmd=names_dict, segment_lengths=SEGMENT_LENGTH,
                            minlen=MIN_LEN, gcode=gencode, bpairs=basepairs, K=K, dc=mod, check_fmt = check_fmt)
    if not line0.startswith('@') and not line0.startswith('>'):
        print('kmermaid: error: only fasta, fastq format supported, exiting...')
        sys.exit(1)

    end_time = timeit.default_timer()
    tot_time = end_time - start_time



