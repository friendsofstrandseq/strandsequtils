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

# Positions of breakpoints
bp$bp_pos = bp$bin * windowsize


### NEW calculation of Sensitivity and Specificity
summary = NULL
for (k_ in unique(bp$k)) {
    
    bp_k = bp[k == k_,]
    bp_k$true = 0
    
    SVs = rbind(sv[, .(chrom, sv_pos = start)], sv[, .(chrom, sv_pos = end)])
    SVs$found = 0
    
    for (i in 1:nrow(bp_k)) {
        # true positives
        # How many SV breakpoints were identified by (at least 1) breakpoint?
        SVs[chrom == bp_k$chrom[i] & abs(sv_pos - bp_k$bp_pos[i]) <= 5*windowsize,]$found = 1
        
        # false positives
        # How many of breakpoints are not close to SVs
        if (nrow(SVs[chrom == bp_k$chrom[i] & abs(sv_pos - bp_k$bp_pos[i]) <= 5*windowsize,])>0) bp_k$true[i] = 1
    }
    summary_k = merge(SVs[, .(tp = sum(found)), by = chrom],
          bp_k[, .(bps = .N, fp = .N - sum(true)), by = chrom],
          by = "chrom")
    summary = rbind(summary, cbind(summary_k, k=k_))
}

# Add total number of 
summary <- merge(summary, sv[, .(true = 2*.N), by = chrom], all.x=T)
S = summary[chrom == "chr1", ]

write.table(S, file = paste0(f_out, ".txt"), row.names = F, quote = F, col.names = T, sep = "\t")
p <- ggplot(S) + 
    aes(fp, tp/true) + 
    geom_line() + 
    geom_label(aes(label = k)) +
    ggtitle(paste("ROC curve", "chr1"))
p
ggsave(p, filename = f_out, width=8, height = 6)



