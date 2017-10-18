# SCE-Utils
Utilities for finding and plotting SCEs in Strand-Seq data.

## Shiny
A web app for semi-manual SCE annotation.

## The Scripts `SCE.R` and `plot_all_into_pdf.R`
A script for finding potential sister chromatid exchanges, and another to plot the results into a pdf.

### Usage from command line:
```
Rscript SCE.R count.table.gz SCEs.txt
Rscript plot_all_into_pdf.R count.table.gz SCEs.txt output.pdf
```
