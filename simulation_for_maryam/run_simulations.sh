#!/bin/bash

echo "Simuating data for Maryam"

mkdir data 2>/dev/null
echo "Distributing SVs..."
Rscript simulate_SVs.R 

echo "Expanding different VAFs..."
for vaf in 1 0.5 0.25 0.1
do
    for size in small large 
    do
        awk -v vaf="$vaf" '{print $0 "\t" vaf}' data/sv_file.${size}.txt > data/sv_config.vaf${vaf}.${size}.txt
    done
done

echo "Simulating cells with MC"
for vaf in 1 0.5 0.25 0.1
do
    for size in small large
    do
        for cov in 5 15
        do
            ls data/sv_config.vaf${vaf}.${size}.txt
            ../../strseq/build/mosaic simulate -w 50000000 -n 100 -c $cov -C $cov -o data/counts.cov${cov}.vaf${vaf}.${size}.txt.gz data/sv_config.vaf${vaf}.${size}.txt
        done
    done
done
