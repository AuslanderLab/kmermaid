import pandas as pd
import numpy as np
import pickle
import sys
from kmermaid import *

SEGMENT_LENGTH=1000000000000
MIN_LEN=99
K = 5
AAs = ['H', 'T', 'W', 'P', 'X', 'I', 'Y', 'R', 'Q', 'K', 'U', 'G', 'A', 'C', 'N', 'M', 'V', 'E', 'L', 'D', 'F', 'S']
RETRAIN = False

##load things, should be moved somewhere else:
print("Loading databases", file=sys.stderr)

with open('../db/clsd0.pkl', 'rb') as fp:
    clsd = pickle.load(fp)

bac_d = segment_fasta(open('../db/expanded_bac.fa'), SEGMENT_LENGTH)

k1 = list(clsd.keys())
v1 = list(clsd.values())
prot2clust = {}
for i in range(len(k1)):
    vi = v1[i]
    for j in range(len(vi)):
        prot2clust[vi[j]] = k1[i]

with open('../db/cluster_representative_sequences.fa') as f:
    l = f.readlines()

nm = [i for i in l if '>' in i]
nmd = {i.split(' ')[0].replace('>',''):i.split('.')[1].split('>')[0] for i in nm}

dc,dfreq = get_kmer_dict_cluster(bac_d, clsd, K)


ddc={}
for k in dc.keys():
    v=dc[k]
    v1=dfreq[k]
    ddc[k]={}
    for i in range(0,len(v)):
        ddc[k][v[i]]=v1[i]

pickle.dump( ddc, open( "../db/combine_retrained.pkl", "wb" ))
