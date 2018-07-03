
# README for SV_ConsistencyCheck_PlottinFunction.R

# plotting function that shows level of agreement for each SV call made by MosaiCatcher
# generates four plots showing "high" - "rare" SV calls located in a sample
# for each SV, plot the number and class of SV calls made at that locus
# for each SV, % of cell with same SV class (cell fraction) is calculated and the best-fit SV is shown (highest cell faction)

# Inputs
# fileLoc == PATH to MosaiCatcher output
# e.g. 
fileLoc <- ('./rpe_mosaiClassifier_20180615/sv_calls/C7_data/100000_fixed_norm.many/')
# File == specific file to plot 
# e.g. 
File <- ('simpleCalls_llr4.txt')

# to run:
SVplotting(fileLoc, File)
# output plots are saved in a new folder "SV_ConsistencyCheck"