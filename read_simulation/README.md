#

## Part 1: simulating from reference genome

### Mason

failed. Didn't even write a fastq file

### `wgsim`

```
R=/g/korbel/shared/datasets/refgenomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
wgsim -N 500000000 -1 64 -2 64 -h -S 1 $R  wgsim.1.fq wgsim.2.fq > wgsim.variants.txt
```


## Part 2: simulating from real Strand-seq data (trio data)

See `sample_from_real_data.

```
for x in data/*.bam; 
    do 
	y=${x#*/}; 
	echo "Processing $x"; 
	samtools view -f 65 -F 2304 $x \
	| awk '{print "@" $1 "\n" substr($10,1,64) "\n+\n" substr($11,1,64)}' \
	> ${y%.bam}.1st.trim64.fq
done
```

