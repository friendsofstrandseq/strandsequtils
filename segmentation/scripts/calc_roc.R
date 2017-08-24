library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly = T)
f_bp = args[1]
f_sv = args[2]
windowsize = as.integer(args[3])
f_out = args[4]

bp = fread(f_bp)
colnames(bp) = c("k","bin","chrom")

sv = fread(f_sv)
colnames(sv) = c("chrom","start","end","type","vaf")


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
chrom_sizes$num_bins = as.integer(floor(chrom_sizes$size / windowsize))




# Positions of SV breakpoints
sv_coords = rbind(sv[, .(chrom, sv_pos = start)], sv[, .(chrom, sv_pos = end)])
sv_coords = rbind(cbind(sv_coords, bin = floor(sv_coords$sv_pos/windowsize)),
                  cbind(sv_coords, bin = ceiling(sv_coords$sv_pos/windowsize)))

# Positions of breakpoints
bp$bp_pos = bp$bin * windowsize

# All bins
all_bins = chrom_sizes[, data.table(bin = 0:num_bins), by = chrom]
all_bins$pos = (all_bins$bin) * windowsize
all_bins <- merge(all_bins, sv_coords, by = c("chrom", "bin"), all.x = T)

# Annotate all bins by whether they match breakpoints and calculate ROC
ROC <- NULL
for (k_ in unique(bp$k)) {
    bkps <- merge(all_bins, bp[k==k_,], by = c("chrom", "bin"), all = T, allow.cartesian = T)
    test <- bkps[, 
                 .(total = length(pos),
                   sv = sum(!is.na(sv_pos)),
                   bp = sum(!is.na(bp_pos)),
                   tp = sum(!is.na(sv_pos) & !is.na(bp_pos)),
                   fp = sum(is.na(sv_pos) & !is.na(bp_pos))), 
                 by = chrom][order(chrom),]
    test <- test[, c("fp_rate","tp_rate") := .(fp / (total-sv), 2*tp/sv)][]
    ROC <- rbind(ROC, cbind(test, k = k_))
}

write.table(ROC, file = paste0(f_out, ".txt"), row.names = F, quote = F, col.names = T, sep = "\t")
p <- ggplot(ROC[chrom == "chr1",]) + 
    aes(fp_rate, tp_rate) + 
    geom_line() + 
    geom_label(aes(label = k)) +
    ggtitle(paste("ROC curve", "chr1"))
ggsave(p, filename = f_out, width=8, height = 6)


