library(data.table)
#args = commandArgs(trailingOnly = T) # args = c(80, "het_del", 400e3, 1, "test.out.txt")

SV_num    = 200
SV_types  = c("het_del", "hom_del", "het_dup", "hom_dup", "het_inv", "hom_inv", "inv_dup", "false_del")
SV_sizes  = c(50e3, 800e3)
SV_vafs   = c(0.1,0.2,0.5,1)
f_out     = "simulation.input.txt"

# GRCh38
chrom_sizes = data.table(
    chrom = c(
        "chr1",  "chr2",  "chr3",   "chr4",  "chr5",  "chr6",
        "chr7",  "chr8",  "chr9",  "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),
    size = c(
        248956422, 242193529, 198295559, 190214555,
        181538259, 170805979, 159345973, 145138636,
        138394717, 133797422, 135086622, 133275309,
        114364328, 107043718, 101991189,  90338345,
        83257441,  80373285,  58617616,  64444167,
        46709983, 50818468, 156040895, 57227415))


SVs = NULL
buffer = 1e6
iterations = 0
repeat {
    # sample chromosome
    chrom_ = sample(chrom_sizes$chrom, 1, prob = chrom_sizes$size / sum(chrom_sizes$size), replace = T)
    
    # sample SV size
    SV_size = as.integer(runif(1, SV_sizes[1], SV_sizes[2]))
        
    # sample position
    new_s = as.integer(runif(1, 1, chrom_sizes[chrom == chrom_]$size - SV_size))
    new_e = new_s + SV_size
    
    # sample SV_type
    SV_type = sample(SV_types,1)
    
    # sample SV VAF
    SV_vaf  = sample(SV_vafs, 1)
    
    # check overlap
    if (is.null(SVs) || all(new_e <= SVs$start - buffer | new_s > SVs$end + buffer) ) {
        SVs = rbind(SVs, data.table(chrom = chrom_, start = new_s, end = new_e, type = SV_type, vaf = SV_vaf))
    }
    iterations = iterations + 1
    if (nrow(SVs) == SV_num || iterations > 15*SV_num) break
}
iterations
SVs <- SVs[order(chrom, start),]
write.table(SVs, file = f_out, quote=F, row.names = F, col.names = F, sep = "\t")

# Then run:
# ../segmentation/scripts/simul -w 50000 -n 100 -o simulation.n200.txt.gz -c 4 -C 20 -s 3 -S simulation.sces.txt simulation.input.txt

# Then segment with MC:
# ../segmentation/scripts/segmentation -m 0.2 -o simulation.n100.segments_mc.txt2 simulation.n100.txt.gz

# Then reformat that table a bit and finally plot it:
# Rscript plot_segmentation.R simulation.n100.txt.gz simulation.input.txt simulation.n100.segments_mc.txt 50000 simulation.n100.segments_tA.pdf

