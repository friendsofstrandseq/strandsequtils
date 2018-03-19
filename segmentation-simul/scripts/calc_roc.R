library(data.table)
library(assertthat)

MAX_DIST_FACTOR = 1 # * window size
genome_sizes = data.table(chrom = paste0("chr", c(1:22,"X","Y")),
                          genome_size  = c( 248956422, 242193529, 198295559, 190214555,
                                     181538259, 170805979, 159345973, 145138636,
                                     138394717, 133797422, 135086622, 133275309,
                                     114364328, 107043718, 101991189,  90338345,
                                      83257441,  80373285,  58617616,  64444167,
                                      46709983,  50818468, 156040895,  57227415))


### 1. Read tables
f_bp = snakemake@input[["bp"]] # f_bp = "segmentation/cells-30_window-50000/cov-10_sces-4/class-het_del_size-400000_vaf-0.2.jointSeg.txt"
f_sv = snakemake@input[["sv"]] # f_sv = "svs/class-het_del_size-400000_vaf-0.2.txt"
windowsize = as.numeric(snakemake@params[["w"]])


bp = fread(f_bp)
colnames(bp) = c("k","bin","chrom")

sv = fread(f_sv)
colnames(sv) = c("chrom","start","end","type","vaf")


### 2. Keep only chromosomes with at least 25 breakpoints
good_chroms = bp[, max(k), by = chrom][V1>25,chrom]
bp <- bp[chrom %in% good_chroms]
sv <- sv[chrom %in% good_chroms]
message("[calc_roc.R] Window size ", windowsize, " (", nrow(sv), " SVs)")
message("[calc_roc.R] Only consider the chromosomes ", paste(good_chroms, collapse = " "))


### 3. Augment bp table with "bp_pos" and "bp_per_Mb"
bp[, bp_pos := bin * windowsize] # calculate window position
bp <- merge(bp, genome_sizes, by = "chrom") # annotate genome size
bp[, bp_per_Mb := round(k/genome_size * 1e6, 2)] # Go from "number segments" to "segments per Mb"
bp <- bp[, .SD[k == max(k),], by = .(chrom,bp_per_Mb)] #


### 4. Keep only the "bp_per_Mb" values for which all chromosomes exist.
good_perMb <- bp[, length(unique(chrom)), by = bp_per_Mb][V1 == length(unique(bp$chrom)), bp_per_Mb]
bp <- bp[bp_per_Mb %in% good_perMb]


### NEW calculation of Sensitivity (Recall) and Specificity ()
summary = NULL
for (bp_per_Mb_ in unique(bp$bp_per_Mb)) {
    
    # Sub-table for a given number of bp_per_Mb
    bpx = bp[bp_per_Mb == bp_per_Mb_,. (bpx.chrom = chrom,
                                        bpx.bp_per_Mb = bp_per_Mb,
                                        bpx.k = k,
                                        bpx.bin = bin,
                                        bpx.bp_pos = bp_pos)]

    # These are the SV breakpoints
    SVs = rbind(sv[, .(chrom, sv_pos = start)], sv[, .(chrom, sv_pos = end)])

    # calculate for each predicted breakpoint whether it matches a real one
    bpx[, bpx.true := nrow(SVs[chrom == bpx.chrom & abs(sv_pos - bpx.bp_pos) <= MAX_DIST_FACTOR*windowsize]), by = .(bpx.chrom,bpx.bp_pos)]
    assert_that(max(bpx$bpx.true)<=1)
    #message("[calc_roc.R] ", nrow(bpx[bpx.true>0]), " of ", nrow(bpx), " prediction are true")

    # Calculate for each real breakpoint whether it was found
    SVs[, found := nrow(bpx[bpx.chrom == chrom & abs(sv_pos - bpx.bp_pos) <= MAX_DIST_FACTOR*windowsize]), by = .(chrom,sv_pos)]
    #assert_that(max(SVs$found)<=1) % this is not always the case: it could be found by 2 breakpoints
    #message("[calc_roc.R] ", nrow(SVs[found>0]), " of ", nrow(SVs), " SV breakpoints were found")


    X = data.table(bp_per_Mb = bp_per_Mb_,
                   SV_found  = nrow(SVs[found>0]),
                   SV_total  = nrow(SVs),
                   bp_true   = nrow(bpx[bpx.true>0]),
                   bp_total  = nrow(bpx))
    summary = rbind(summary, X)
}

write.table(summary,
            file = paste0(snakemake@output[[1]]),
            row.names = F,
            quote = F,
            col.names = T,
            sep = "\t")
