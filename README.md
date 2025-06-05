# Extract-all-genomes

designed to compare sequence strings within large genome databases by percent identity and structure

genomes.py will download your chosen organisms genomes from ncbi GenBank

sequence_logo.py will compare the percent identity of user supplied sequence string within each genome, and the nt positions of these sequences

## Comparing secondary structures between genomes

create a covariance model using the program infernal.

- can supply multiple sequences
- create a .sto file with format:

  \# STOCKHOLM 1.0

  seq1 ATCG

  =GC SS_cons )..( \#-->secondary strucure plot denoted with ).


  //
- Build the cm model: using bash, an example command might be 'cmbuild xx.cm seed_alignment.sto'
- calibrate said model: 'cmcalibrate xx.cm'
- search ncbi built fasta file from genomes.py and output table 'cmsearch --tblout results.tbl xx.cm xx.fasta > results.txt'

compare_and_extract.py will extract the output results.tbl and take the sequences from every matched genome and provide a sequence logo plot to give a sense of the sequences which match a certain secondary structure within a genome - helpful to show if these structures share any regions of sequence similarity. 

### Requirements
local installations of **Biopython**, **clustalo**, **blast** and **infernal** 

available in most package managers for linux

python libraries for:

logomaker

matplotlib

os

tkinter

argparse

pandas

subprocess

os

tempfile

###### Many of these libraries are included in biopython. Install the others using 'pip install xxx' unless python is managed system wide, in which case you may need a virtual environment for these packages

#### Updates
2025/06/04: genomes.py rewritten to enable user prompts for organisms and email.
