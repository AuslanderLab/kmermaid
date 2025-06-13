### <i>k</i>Mermaid: Ultrafast metagenomic read assignment to protein clusters by hashing of amino-acid k-mer frequencies

## Overview 
This file describes the software package <i>k</i>mermaid [1], a k-mer based method for functional classification of metagenomic reads into protein clusters. 

Operating systems tested: Linux (CentOS version 7 and Ubuntu 20.04 LTS) and MacOS.

Python versions tested: Python 3.6, 3.7, 3.8 and 3.9

## Citation
[1] Anastasia Lucas, Daniel Schaffer, Jayamanna Wickramasinghe, Noam Auslander. kmermaid: Ultrafast functional classification of microbial reads

## Installation

Typical install time on a "normal" desktop computer: less than 30 minutes (depending on the number of packages already installed)


### Install with pip

To install the current version of this Github repo, run the following commands
```
git clone https://github.com/AuslanderLab/kmermaid.git 
cd kmermaid
pip install .
```

#### File dependncy: To download the larger files uploaded to this repo, install git-lfs https://git-lfs.com (otherwise larger files including the kmer model would not be downloded through git clone). Alternatively, download the files directly. For example, to download the kmer model to ```kmermaid/db/``` run:
```
cd kmermaid/db/
rm kmer_model.pkl
wget https://zenodo.org/records/15658544/files/kmer_model.pkl 
cd ../../
```

### Install with conda 
1- Get anaconda (64 bit)installer python3.x for linux : https://www.anaconda.com/download/#linux <br />
2- Run the installer : bash Anaconda3-2021.11-Linux-x86_64.sh, and follow the instructions to install anaconda at your preferred directory

#### Run the following commands: <br />
```
git clone https://github.com/AuslanderLab/kmermaid.git
cd kmermaid
conda create --name kmermaid python=3.8 pip
conda activate kmermaid
pip install .
```

#### To deactivate kmermaid environment: <br />
```
conda deactivate
```

## Usage
To use kmermaid, a user must provide an input fasta/fastq and is recommended to provide an output path:

### Running example (Demo):

Expected run time for demo on a "normal" desktop computer: less than 1 minute

a. To run with an example input fasta file (```inputs/reads_file.fa```) run

```
kmermaid --input inputs/reads_file.fa --output outputs/out_file
```

And evaluate the output file generated in ```outputs/``` using the expected output in  ```expected_output/expected_out_exmp.tsv```

To test if the above command worked as expected, run the additional command

```
 diff outputs/out_file.tsv expected_output/expected_out_exmp.tsv
```
The installation is correct if the above diff command retruns either no differences or small differences in the less significant digits.
 

b. To run the example remote homology sequences not classified with blastx as described in the manuscript (```inputs/remote_homology_sequences.txt```) run

```
kmermaid --input inputs/remote_homology_sequences.txt --output outputs/remote_ho
```


## Output
kmermaid output is the K-mer based cluster classification of each read which is a tab delimited file with the following columns:
1) seq_name-read id from input fasta 
2) cluster_rep-protein ID of cluster representative
3) prot_name- name of the protein 
4) score - confidence scores when above 3

for example, an output line from <i>k</i>mermaid:
`
WP_002358485.1_0	WP_002358485.1	lantipeptide cytolysin subunit CylL-L	23.00
`

### Parameters description:

| Parameter |     type      |           description            |       default       |
| :---: |:-------------:|:--------------------------------:|:-------------------:|
| input |  path (txt)   |    fasta or fastq input file     |          -          |
| output |  path (txt)   |       path to output file        | False (no argument) |
| cluster_reps |  path (txt)   | path to cluster names (pkl file) | False (no argument) |
| trained_model |  path (txt)   |    path to cluster (pkl file)    | False (no argument) |
| append_path | Boolean (0/1) |       use flag with slurm        |         `1`          |

## Retraining
To retrain <i>k</i>mermaid with a new clustered database of protein sequences, the retraining code is provided in ```kmermaid/_retrain_kmermaid_.py```

<b>WARNING:</b> retraining <i>k</i>mermaid  with a new clustered database is not recommended without sufficient benchmarking and evaluation, and is therefore at user responsibility.

## Contact

If you have any questions or encounter any difficulties, please create an issue on Github or email us at either Anastasia.Lucas@pennmedicine.upenn.edu or nauslander@wistar.org.


