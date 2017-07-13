library(data.table)
library(ggExtra)
library(ggplot2)
library(scales)

f_in = paste("/usr/local/bin/zcat",
             "/Volumes/korbel/meiers/projects/20170105_strandseq_code/mosaicatcher/data/RPE/BM160815_WT.500k.txt.gz")
d = fread(f_in)
d = merge(d, d[, .(total = sum(w+c)), by = cell], by = "cell")

# Sex ratio - but all cells are female 
#sex = d[,.(Y   = sum(w[chrom=="Y"]+c[chrom=="Y"]), 
#           X   = sum(w[chrom=="X"]+c[chrom=="X"]),
#           ten = sum(w[chrom=="10"]+c[chrom=="10"]),
#           total = sum(w+c)), by = cell]
#ggplot(sex[total > 1e5,]) + 
#    aes(X/total, Y/total) +
#    geom_point(alpha=0.3)
#ggplot(sex[total > 1e5,]) + 
#    aes(total, ten/total) + 
#    geom_point(alpha=0.3)



# Overview of crick ratio along chromosomes
p0 <- ggplot(d[total > 1e5 & chrom == "10",]) +
    aes(x = (start + end)/2, y = c/(c+w)) +
    geom_point(size = 0.1) +
    facet_grid(cell ~ .) +
    theme_bw() + 
    scale_x_continuous(breaks=pretty_breaks(15)) +
    scale_y_continuous(breaks=pretty_breaks(3)) +
    ggtitle("Shift in crick ratio on chr10")
ggsave(p0, filename = "rpe_translocation.wcfrac_overiew.pdf", width=10, height = 30)


# Watson crick median counts left. vs. right of translocation
lr = merge(d[total > 1e5 &start>=40e6 & end<=60e6 & chrom=="10", .(w_left = median(w), c_left = median(c), frac_left = median(c/(c+w))), by = cell],
           d[total > 1e5 &start>=62e6 & end<=82e6 & chrom=="10", .(w_right = median(w), c_right = median(c), frac_right = median(c/(c+w))), by = cell],
           by = "cell")
lr <- lr[, `X state q-arm` := d[total > 1e5 & start > 130e6 & chrom == "X", .(X_state = names(table(class))[which.max(table(class))]), by = cell]$X_state ]
lr <- lr[, `X state p-arm` := d[total > 1e5 & start < 25e6 & chrom == "X", .(X_state = names(table(class))[which.max(table(class))]), by = cell]$X_state ]
lr <- lr[, `chromosome state` := ifelse(c_left/(c_left+w_left) < 0.2, "CC", ifelse(w_left/(c_left+w_left)<0.2, "WW", "WC")) ][]

p1 = ggplot(lr) +
    aes(c_left/(c_left+w_left), c_right/(c_right+w_right), col = `X state p-arm`) + 
    geom_point(alpha = 0.5) + 
    xlab("Crick fraction chr10:40-60Mb") +
    ylab("Crick fraction chr10:62-82Mb") +
    theme_classic() +
    theme(legend.position = "bottom") +
    ggtitle("Strand state around translocation")
p1_ = ggMarginal(p1, type = "histogram", binwidth = 0.02, fill = "dodgerblue3", col = "dodgerblue3", alpha=0.3)
ggsave(p1_, filename = "rpe_translocation.local_state_change.left_arm.pdf", width=4,height=4)

p1 = ggplot(lr) +
    aes(c_left/(c_left+w_left), c_right/(c_right+w_right), col = `X state q-arm`) + 
    geom_point(alpha = 0.5) + 
    xlab("Crick fraction chr10:40-60Mb") +
    ylab("Crick fraction chr10:62-82Mb") +
    theme_classic() +
    theme(legend.position = "bottom") +
    ggtitle("Strand state around translocation")
p1_ = ggMarginal(p1, type = "histogram", binwidth = 0.02, fill = "dodgerblue3", col = "dodgerblue3", alpha=0.3)
ggsave(p1_, filename = "rpe_translocation.local_state_change.right_arm.pdf", width=4,height=4)
