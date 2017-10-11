#!/bin/bash
module load Boost HTSlib

echo "Simuating data for Maryam"

mkdir data 2>/dev/null

if [ ! -f data/sv_file.small.txt ]
then
    echo "Distributing SVs..."
    Rscript simulate_SVs.R 
fi

for vaf in 1 0.5 0.25 0.1 0.08 0.05 0.02
do
    for size in small large 
    do
        if [ ! -f data/sv_config.vaf${vaf}.${size}.txt ]
        then
            echo "Expanding different VAFs..."
            awk -v vaf="$vaf" '{print $0 "\t" vaf}' data/sv_file.${size}.txt > data/sv_config.vaf${vaf}.${size}.txt
        fi
    done
done

echo "Simulating cells with MC"
for vaf in 1 0.5 0.25 0.1 0.08 0.05 0.02
do
    for size in small large
    do
        for cov in 5 15 45
        do
            for p in 0.8 0.3 0.1 0.03
            do
                if [ ! -f data/counts.cov${cov}.vaf${vaf}.${size}.p${p}.txt.gz ]
                then
                    echo -e "Size ${size},\tVAF ${vaf},\tcov ${cov}\tp $p"
                    ../../mosaicatcher_new/build/mosaic simulate -w 50000 -n 200 \
                        -p $p \
                        -c $cov -C $cov \
                        -o data/counts.cov${cov}.vaf${vaf}.${size}.p${p}.txt.gz \
                        -S data/sces.cov${cov}.vaf${vaf}.${size}.p${p}.txt \
                        -V data/variants.cov${cov}.vaf${vaf}.${size}.p${p}.txt \
                        -U data/breakpoints.cov${cov}.vaf${vaf}.${size}.p${p}.txt \
                        -P data/phases.cov${cov}.vaf${vaf}.${size}.p${p}.txt \
                        data/sv_config.vaf${vaf}.${size}.txt
					Rscript qc.R data/counts.cov${cov}.vaf${vaf}.${size}.p${p}.txt.gz data/plot.cov${cov}.vaf${vaf}.${size}.p${p}.pdf
                fi
            done
        done
    done
done

