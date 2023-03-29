#!/bin/bash

bowtie_cores=16
HOME=`pwd`

lib=$1
sgefile="./config/qsub.$lib.sh"

mkdir $HOME/termination/

echo """#!/bin/bash

mkdir $HOME/termination/$lib
cd $HOME/termination/$lib
cp $HOME/reference/WS230/ce_WS230.pirna.improve.fa ./

cp $HOME/results/$lib/$lib.inserts.gz ./
gunzip $lib.inserts.gz

cat $lib.inserts | awk '{print \">\" \$1 \"\n\" \$1; }' > $lib.fa
bowtie-build $lib.fa $lib.index --threads $bowtie_cores

bowtie $HOME/termination/$lib/$lib.index ce_WS230.pirna.improve.fa -p $bowtie_cores -v 0 --best --strata -a -f > $lib.pirna.bowtie
awk '{FS = \"\t\"} ; {print \$1 \"\t\" \$5 \"\t\" \$3; }' $lib.pirna.bowtie > $lib.bed
awk '{print \">\" \$1 \",\" \$5 \",\" \$3 \"\n\" \$3; }' $lib.pirna.bowtie > $lib.bed.fa

bowtie $HOME/reference/WS230/c_elegans.WS230.genomic $lib.bed.fa -p $bowtie_cores -v 0 --best --strata -a -f > $lib.precursor.v0.bowtie

# bowtie $HOME/reference/WS230/c_elegans.WS230.genomic $lib.bed.fa -p $bowtie_cores -v 1 --best --strata -a -f > $lib.precursor.v1.bowtie

cd $HOME
Rscript $HOME/scripts/construct_precursor_df.R $lib $bowtie_cores
""">$sgefile

echo $sgefile
bash $sgefile