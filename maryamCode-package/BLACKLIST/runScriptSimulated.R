directory = "/local/data/maryam/codes/strand-seq/SVcalling/simulatedData/dataset2/newFormat/"

for (dir in system(paste("ls", directory), intern = T))
{
  system(paste0("mkdir ",directory, dir, "/NBfitPlots"))
  system(paste0("Rscript runCalssificationOnSimulatedData.R ", directory, dir, "/ &disown"))
}
