# Brown_2023

# Sensitized piRNA reporter identifies multiple RNA processing factors involved in piRNA-mediated gene silencing

This repository holds all scripts necessary to perform small RNA alignment and plot construction described in Brown et. al. 2023.

## Summary of scripts

### Preparation and requirements

To run the pipeline script, bowtie v1.2.1.1, bedtools v2.30.0, and R v4.0.3 must be installed.
The reference directory should contain subdirectories for the WS230 reference genome. All subdirectories should contain the WS230 genome in FASTA format, available from Wormbase.
A bowtie index is necessary for alignment to genome, splice junctions, structural RNAs and miRNA hairpins. To build all necessary bowtie index files, run
```

HOME=`pwd`
cd $HOME/reference/WS230/
path_to_genome_fasta=$HOME/reference/WS230/c_elegans.WS230.genomic.fa
bowtie-build ${path_to_genome_fasta} c_elegans.WS230.genomic
bowtie-build ce_WS230.coding_transcript.juncs.fa ce_WS230.coding_transcript.juncs
bowtie-build ce_WS230.rna.knownRNA.fa ce_WS230.rna.knownRNA
bowtie-build ce_hairpin.dna.fa ce_hairpin.dna

```

Compressed fastq.gz files must be located in the fastq directory.

### Alignment script

To map reads contained in the fastq directory and calculate normalized reads per million against genes, run bowtie_alignment.sh
```
bash ./bowtie_alignment.sh
```
The directory config should contain a TXT file with 3 comma separated fields fields: path to the fastq.gz file, name of the library, and the genome to use for alignment.

### piRNA precursor script

To find reads that correspond to piRNA precursor molecules, the precursor.pipeline.sh script must be run with relevant libraries
```
bash ./precursor.pipeline.sh JB.20181018_RRS.L4440
bash ./precursor.pipeline.sh JB.20181018_RRS.snpc4
bash ./precursor.pipeline.sh JB.20181018_RRS.L4440_ints
bash ./precursor.pipeline.sh JB.20181018_RRS.ints1
bash ./precursor.pipeline.sh JB.20181018_RRS.dic1
bash ./precursor.pipeline.sh JB.20201013_ev
bash ./precursor.pipeline.sh JB.20201013_ints1
```
The alignment script must be run prior to this step.

### Browser image and metagene construction

Run plots.R to use data generated in the bowtie_alignment.sh and precursor.pipeline.sh scripts to construct plots shown in Figure 2B, Figure S2C, Figure 3A, and Figure 3B.
This script depends on R libraries:
ggplot2
scales
reshape2
ggpubr
gridExtra
dplyr
rstatix

```
Rscript ./plots.R
```
