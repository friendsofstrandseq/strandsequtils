library(data.table)
library(stringr)
require(scales)

format_Mb <- function(x) {
  y = as.character(x)
  y[x >=  1e6] = paste(round(x[x >= 1e6]/1e6,1), "Mb") 
  y[x <   1e6] = paste(round(x[x <  1e6]/1e3,0), "kb") 
  return(y)
}


read_files <- function(input_list) {
  
  D = NULL
  regex = "roc/cells-(\\d+)_window-(\\d+)/cov-(\\d+)_sces-(\\d+)/class-([a-zA-Z_]+)_size-(\\d+)_vaf-([0-9.]+)\\.([a-zA-Z]+)\\.txt"
  
  # Read all input files and extract variables from their names
  for (f in input_list) {
    
    # Get variables
    vars = str_match(f,regex)
    
    # Read table
    d = fread(f)
    d$n_cells  = as.integer(vars[2])
    d$window   = as.integer(vars[3])
    d$cov      = as.integer(vars[4])
    d$sces     = as.integer(vars[5])
    d$SV_class = vars[6]
    d$SV_size  = as.integer(vars[7])
    d$SV_vaf   = as.numeric(vars[8])
    d$Method   = vars[9]
    
    # Store in table
    D = rbind(D, d)
  }
  
  D[, SV_class := vapply(SV_class, function(x) {switch(x, 
                                                       het_dup = "Duplication (het)", 
                                                       hom_dup = "Duplication (hom)",
                                                       het_inv = "Inversion (het)",
                                                       hom_inv = "Inversion (hom)",
                                                       het_del = "Deletion (het)",
                                                       hom_del = "Deletion (hom)",
                                                       inv_dup = "Inv. duplication")}, "char")]
  D[, SV_size := factor(vapply(SV_size, format_Mb, "char"), 
                        levels = vapply(sort(unique(D$SV_size)), format_Mb, "char"), 
                        ordered = T)]
  D[, SV_vaf  := factor(vapply(SV_vaf, function(x){paste(x*100,"%")}, "char"),
                        levels = vapply(sort(unique(D$SV_vaf)), function(x){paste(x*100,"%")}, "char"),
                        ordered = T)]
  D
}


