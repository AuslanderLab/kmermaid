from collections import OrderedDict
from collections import Counter
from collections import defaultdict
import random
import sys
import timeit
import re
import threading
import time

######################### retrain_kmer.py ##########################

flatten = lambda l: [item for sublist in l for item in sublist]


def get_kmer_dict_cluster(bac_d, clsd, K):
    '''
    train k-mer model mapping reads to clusters based on average k-mer frequencies in each cluster
    :param bac_d:
    :param clsd:
    :param K:
    :return:
    '''
    dc = {}  # dictionary mapping kmers to cluster reps
    dfreq = {}  # dictionary mapping (kmer,cluster_rep) to score, may be two large and need changing strategy
    cluster_reps = list(clsd.keys())
    cnt=0
    for rep in cluster_reps: ##for every cluster
        print(cnt)
        prot_list = clsd[rep] ##protein list in cluster
        # all kmers of proteins in the cluster
        vkmer = [[bac_d[p][0][i:i + K] for i in range(len(bac_d[p][0]) - (K-1))] for p in prot_list]
        flat_kmr=flatten(vkmer) #flat list of all kmers in cluster
        count_km = Counter(flat_kmr)
        uniqe_kmers = set(flat_kmr)
        lp = len(prot_list)
        ##get cluster frequency of each kmer
        for uk in uniqe_kmers:
            # print(uk)
            val = count_km[uk]/lp
            if uk not in dc:
                dc[uk] = [rep]
                dfreq[uk] = [val]
            else:
                dc[uk].append(rep)
                dfreq[uk].append(val)

        cnt+=1
    return dc,dfreq

######################### get_kmer_script_fast.py ##########################

#@profile
def translate_frameshifted(sequence, gcode):
    """
    Translation of nucleotide to amino acid
    :param sequence: a section of nucleotide sequence
    :param gcode: gencode dictionary
    :return: amino acide sequence
    """
    translate = ''.join([gcode.get(sequence[3 * i:3 * i + 3]) for i in range(len(sequence) // 3)])
    return translate

#@profile
def reverse_complement(sequence, bpairs):
    """
    Genertate the reversed sequence
    :param sequence: a section of nucleotide sequence
    :param bpairs: basepairs dictionary
    :return: the reversed string of the nucleotide sequence
    """
    reversed_sequence = (sequence[::-1])
    rc = ''.join([bpairs.get(reversed_sequence[i]) for i in range(len(sequence))])
    return rc


def six_frame_trans(seq, gcode, bpairs):
    x1 = translate_frameshifted(seq, gcode)
    x2 = translate_frameshifted(seq[1:], gcode)
    x3 = translate_frameshifted(seq[2:], gcode)
    rc = reverse_complement(seq, bpairs)
    x4 = translate_frameshifted(rc, gcode)
    x5 = translate_frameshifted(rc[1:], gcode)
    x6 = translate_frameshifted(rc[2:], gcode)
    x = [x1, x2, x3, x4, x5, x6]
    return x
    
def handle_non_ATGC(sequence):
    """
    Handle non ATGCs.
    :param sequence: String input.
    :return: String output (only ATCGs), with randomly assigned bp to non-ATGCs.
    """
    ret = re.sub('[^ATCG]', random.choice(['A', 'C', 'G', 'T']), sequence)
    assert len(ret) == len(sequence)
    return ret

def proc_classify_fasta(fasta, outfilepath, nmd, segment_lengths, minlen, gcode, bpairs, K, dc, print_prog=True):
    """
    Parse Fasta into segments.
    :param fasta: File handle in Fasta format.
    :param segment_lengths: Length of segments for model.
    :return: Dictionary of Fasta name -> list of segments.
    """
    dck = set(dc.keys())
    seq_name = None
    seq_value = None
    fo = open(outfilepath,'w')
    fo.write('seq_name' + "\t" + 'cluste_rep' + "\t" + 'prot_name' + "\t"+'score'+"\n")
    count = 0
    for line in fasta:
        # print(line)
        line = line.strip()
        if line.startswith(">"):
            # Reached a new sequence in Fasta, if this is not the first sequence let's save the previous one
            if seq_name is not None:  # This is not the first sequence
                #assert seq_value is not None, seq_name

                if len(seq_value) > minlen:
                    p, s = map_sequence_to_prot(handle_non_ATGC(seq_value), gcode, bpairs, K, dc, dck)

                    if s>=3:
                        fo.write(seq_name + "\t" + p + "\t" + nmd[p] + "\t" + str(round(s, 2)) + "\n")


                count += 1

            seq_name = line.lstrip(">").split()[0]
            seq_value = ""
        else:
            seq_value += line

    count += 1
    ret_str = "Finished"
    print("Read Fasta entry {}, total {} entries read".format(seq_name, count), file=sys.stderr)
    fo.close()
    return ret_str


def proc_classify_fastq(fastq, outfilepath, nmd, segment_lengths, minlen, gcode, bpairs, K, dc, print_prog=True):
    """
    Parse Fasta into segments.
    :param fasta: File handle in Fasta format.
    :param segment_lengths: Length of segments for model.
    :return: Dictionary of Fasta name -> list of segments.
    """
    # ret = OrderedDict()
    dck = set(dc.keys())
    fo = open(outfilepath,'w')
    fo.write('seq_name' + "\t" + 'cluste_rep' + "\t" + 'prot_name' + "\t"+'score'+"\n")
    count = 0
    for line in fastq:
        if (count % 4)==0:
            seq_name = line.lstrip("@").split()[0]

        if (count % 4)==1:
            seq_value = line
            if len(seq_value) > minlen:
                p, s = map_sequence_to_prot(handle_non_ATGC(seq_value), gcode, bpairs, K, dc, dck)
                if s >= 3:
                    fo.write(seq_name + "\t" + p + "\t" + nmd[p] + "\t" + str(round(s, 2)) + "\n")

        count += 1

    ret_str = "Finished"
    print("Read Fastq entry {}, total {} entries read".format(seq_name, int(count/4)), file=sys.stderr)
    fo.close()
    return ret_str

def get_time():
   t = time.localtime()
   current_time = time.strftime("%m/%d/%Y, %I:%M:%S", t)
   print(str(current_time)+"\n")


def map_sequence_to_prot(seq, gcode, bpairs, K, dc, dck):
    fr0 = six_frame_trans(seq, gcode, bpairs)  ##maybe?Take only complete frames and remove stops
    fr = [i for i in fr0 if '_' not in i] ##fix - if all frames are truncated do not classify
    glob_max_s = 0
    glob_max_p = ""
    for v in fr:
        psd = {}
        #Appears that dc is of the form {kmer -> {protein id -> frequency}}
        items = [dc[j] for j in set([v[i:i + K] for i in range(len(v) - K)]) & dck] ##from kmer to prot,freqs
        for itm in items: ##items is the {k -> v} for kmers in the overlap
            for key, val in itm.items():
                if key in psd: ##key is prot id, val is freq
                    # basically this just provides more confidence when multiple kmers are present 
                    # in cluster and then assigns to the one with the highest composite score based
                    # on all clusters
                    s_new = psd[key] + val
                else:
                    s_new = val
                if s_new > glob_max_s:
                    glob_max_s = s_new
                    glob_max_p = key
                psd[key] = s_new
    return glob_max_p,glob_max_s

