# Extract-all-genomes

designed to compare sequence strings within large genome databases by percent identity and structure

genomes.py will download your chosen organisms genomes from ncbi GenBank

sequence_logo.py will compare the percent identity of user supplied sequence string

## Comparing secondary structures between genomes

create a covariance model using the program infernal.

- can supply multiple sequences
- create a .sto file with format:

  \# STOCKHOLM 1.0

  seq1 ATCG

  =GC SS_cons (. -->secondary strucure plot denoted with ).


  //
- using bash example command might be 'cmbuild xx.cm seed_alignment.sto
- calibrare model 'cmcalibrate xx.cm'
- search above ncbi built fasta file from genome.py and output table 'cmsearch --tblout results.tbl xx.cm xx.fasta > results.txt'

compare_and_extract.py will extract the output results.tbl and take the sequences from every matched genome and provide a sequence logo plot to give a sense of the sequences which match a certain secondary structure within a genome. 
