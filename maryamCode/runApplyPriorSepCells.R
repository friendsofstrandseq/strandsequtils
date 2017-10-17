directory = "/local/data/maryam/codes/strand-seq/SVcalling/simulatedData/dataset2/"
priorFile = "/local/data/maryam/codes/strand-seq/SVcalling/priorProbsForGenotypes.txt"
#setwd(directory)
subdirs = system(paste0("ls ", directory, "newFormat"), intern = T)

for (subdir in subdirs)
{
  print(subdir)
  system(paste("Rscript applyPriorAndSeperateCellPobs.R", priorFile, directory, subdir))
}

print("done")
