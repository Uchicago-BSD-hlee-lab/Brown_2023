#!/bin/bash

bowtie_cores=16
HOME=`pwd`

lib_array=(JB.20181018_RRS.L4440 JB.20181018_RRS.snpc4 JB.20181018_RRS.L4440_ints JB.20181018_RRS.ints1 JB.20181018_RRS.dic1 JB.20201013_ev JB.20201013_ints1 Kasper.2014_control Kasper.2014_snpc4)

mkdir $HOME/termination/

precursor_pipeline () {
    local lib=$1

    mkdir $HOME/termination/$lib
    cd $HOME/termination/$lib
    cp $HOME/reference/WS230/ce_WS230.pirna.improve.fa ./

    cp $HOME/results/$lib/$lib.inserts.gz ./
    gunzip $lib.inserts.gz

    cat $lib.inserts | awk '{print ">" $1 "\n" $1; }' > $lib.fa
    bowtie-build $lib.fa $lib.index --threads $bowtie_cores

    bowtie $HOME/termination/$lib/$lib.index ce_WS230.pirna.improve.fa -p $bowtie_cores -v 0 --best --strata -a -f > $lib.pirna.bowtie
    awk '{FS = "\t"} ; {print $1 "\t" $5 "\t" $3; }' $lib.pirna.bowtie > $lib.bed
    awk '{print ">" $1 "," $5 "," $3 "\n" $3; }' $lib.pirna.bowtie > $lib.bed.fa

    bowtie $HOME/reference/WS230/c_elegans.WS230.genomic $lib.bed.fa -p $bowtie_cores -v 0 --best --strata -a -f > $lib.precursor.v0.bowtie

    # bowtie $HOME/reference/WS230/c_elegans.WS230.genomic $lib.bed.fa -p $bowtie_cores -v 1 --best --strata -a -f > $lib.precursor.v1.bowtie

    cd $HOME
    Rscript $HOME/scripts/construct_precursor_df.R $lib $bowtie_cores

}

for ((i=0;i<=$((${#lib_array[@]}-1));i++)); do
	precursor_pipeline ${lib_array[i]}
done
wait
