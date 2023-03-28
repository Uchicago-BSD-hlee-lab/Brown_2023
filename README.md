# Brown_2023: Sensitized piRNA reporter identifies multiple RNA processing factors involved in piRNA-mediated gene silencing

This repository holds all scripts necessary to perform small RNA alignment and plot construction described in Brown et. al. 2023.

## Summary of scripts

### Preparation and requirements

First, clone repository into working directory:

```
git clone https://github.com/jordan-scot-brown/Brown_2023.git
```

To run the pipeline script, bowtie v1.2.1.1, bedtools v2.30.0, R v4.0.3, and various R packages must be installed. If these requirements are met, then proceed to genome index building with bowtie. Otherwise, build or pull docker image for a complete environment which reproduces the environment used by the authors:

#### Docker image

To build a docker image from the Dockerfile in this repository, run the following in your Brown_2023/ directory:

```
docker build -t brown2023 -f Dockerfile ./
```

Then, to run the container:

```
docker run -it --name brown2023_c1 -p 8889:8889 -v $GIT_PATH/:/project brown2023
```

Alternatively, pull the docker image from DockerHub:

#### Install bowtie and sratoolkit

To fully reproduce the results shown in this paper, install bowtie 1.2.1.1. First, download binaries from https://sourceforge.net/projects/bowtie-bio/ and extract the folder in the same directory as you cloned this repository.

Additionally, sratoolkit is necessary to pull raw reads from NCBI. Download and extract sratoolkit binaries for the docker environment in the same directory:

```
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
```

To associate these binaries with the bash profile of the docker envrionment, add the following lines to /root/.bashrc and source:
```
export PATH=$PATH:/project/bowtie-1.2.1.1-linux-x86_64/bowtie-1.2.1.1
export PATH=$PATH:/project/sratoolkit.3.0.1-ubuntu64/bin
source /root/.bashrc
```

#### get fastq files generated in this paper from NCBI
```
mkdir fastq
cd fastq
fastq-dump --gzip SRR23256087
mv SRR23256087.fastq.gz ints1_capseq.fastq.gz
fastq-dump --gzip SRR23256088
mv SRR23256088.fastq.gz control_capseq.fastq.gz
fastq-dump --gzip SRR23256089
mv SRR23256089.fastq.gz dic1_sRNA.fastq.gz
fastq-dump --gzip SRR23256090
mv SRR23256090.fastq.gz ints1_sRNA.fastq.gz
fastq-dump --gzip SRR23256091
mv SRR23256091.fastq.gz control_ints1_dic1_sRNA.fastq.gz
fastq-dump --gzip SRR23256092
mv SRR23256092.fastq.gz snpc4_sRNA.fastq.gz
fastq-dump --gzip SRR23256093
mv SRR23256093.fastq.gz control_snpc4_sRNA.fastq.gz
fastq-dump --gzip SRR23256094
mv SRR23256094.fastq.gz prp17_sRNA.fastq.gz
fastq-dump --gzip SRR23256095
mv SRR23256095.fastq.gz npp7_sRNA.fastq.gz
fastq-dump --gzip SRR23256096
mv SRR23256096.fastq.gz control_npp7_prp17_sRNA.fastq.gz
```

SRA numbers almost completely mixed up!
As of March 2023, reassign names:

```
mv npp7_sRNA.fastq.gz control_capseq.fastq.gz
mv prp17_sRNA.fastq.gz npp7_sRNA.fastq.gz
mv control_snpc4_sRNA.fastq.gz prp17_sRNA.fastq.gz
mv snpc4_sRNA.fastq.gz control_snpc4_sRNA.fastq.gz
mv control_ints1_dic1_sRNA.fastq.gz snpc4_sRNA.fastq.gz
mv ints1_sRNA.fastq.gz control_ints1_dic1_sRNA.fastq.gz
mv dic1_sRNA.fastq.gz ints1_sRNA.fastq.gz
mv control_capseq.fastq.gz dic1_sRNA.fastq.gz
```

#### get fastq files for snpc-4 RNAi from Kasper et al 2014 from NCBI
```
fastq-dump --gzip SRR1054267
mv SRR1054267.fastq.gz Kasper2014_control.fastq.gz
fastq-dump --gzip SRR1054268
mv SRR1054268.fastq.gz Kasper2014_snpc4.fastq.gz
```
#### Build bowtie index for genome and known RNAs
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

Note: Compressed fastq.gz files must be located in the fastq directory.

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
