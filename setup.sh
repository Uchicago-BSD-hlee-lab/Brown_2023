#!/bin/bash

# get docker image from DockerHub and initialize a container brown2023_c1

# setup necessary packages not included in docker image

# bowtie
# add path to bash profile
# download bowtie/1.2.1.1 binaries from sourceforge project page:
# https://sourceforge.net/projects/bowtie-bio/
# extract in working directory
# export PATH=$PATH:$PWD/bowtie-1.2.1.1-linux-x86_64/bowtie-1.2.1.1

# sra-toolkit
# add path to bash profile
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
# export PATH=$PATH:$PWD/sratoolkit.3.0.1-ubuntu64/bin

# get fastq files generated in this paper from NCBI
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

# SRA numbers almost completely mixed up!
# As of March 2023, reassign names:
mv npp7_sRNA.fastq.gz control_capseq.fastq.gz
mv prp17_sRNA.fastq.gz npp7_sRNA.fastq.gz
mv control_snpc4_sRNA.fastq.gz prp17_sRNA.fastq.gz
mv snpc4_sRNA.fastq.gz control_snpc4_sRNA.fastq.gz
mv control_ints1_dic1_sRNA.fastq.gz snpc4_sRNA.fastq.gz
mv ints1_sRNA.fastq.gz control_ints1_dic1_sRNA.fastq.gz
mv dic1_sRNA.fastq.gz ints1_sRNA.fastq.gz
mv control_capseq.fastq.gz dic1_sRNA.fastq.gz

# get fastq files for snpc-4 RNAi from Kasper et al 2014 from NCBI
fastq-dump --gzip SRR1054267
mv SRR1054267.fastq.gz Kasper2014_control.fastq.gz
fastq-dump --gzip SRR1054268
mv SRR1054268.fastq.gz Kasper2014_snpc4.fastq.gz