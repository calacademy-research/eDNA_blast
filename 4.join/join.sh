#!/usr/bin/env bash
loci=( 16S 18S COI 28S )
samples=( 'SouthShore_April' 'SouthShore_May' 'stream_May')

for sample in "${samples[@]}"
do

    for locus in "${loci[@]}"
    do
        echo "Processing ${sample}_${locus}"
        python ~/dbcAmplicons/bin/dbcAmplicons join -t 4 -O ./${sample}_${locus} -1 ../3.preprocess/${sample}.intermediate/${locus}_R1.fastq.gz -2 ../3.preprocess/${sample}.intermediate/${locus}_R2.fastq.gz > $sample\_$locus.log
    done
done



