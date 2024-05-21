x <- new_xModel(".")

otype(x)
otype(x$GSE138309)

class(x)
class(x$GSE138309)

s3_class(x)
s3_class(x$GSE138309)


methods(print)





target_dir <- "E:/UniTo Drive/WORKS/0010 - Ongoing/Endothelion/data/in/Lines/hCMEC_D3" # @ HOME
target_dir <- "D:/UniTo Drive/WORKS/0010 - Ongoing/Endothelion/data/in/Lines/hCMEC_D3" # @ DBIOS



series_ID <- "GSE195781"

series_GSE195781 <- new_xSeries("GSE195781", target_dir)




hCMEC_D3 <- new_xModel(target_dir)


size.xSeries(hCMEC_D3$GSE138309)
size.xSeries(hCMEC_D3$GSE139133)
size.xSeries(hCMEC_D3$GSE195781)
length(hCMEC_D3$GSE195781)

a <- countMatrix.xSeries(hCMEC_D3$GSE138309)
a <- countMatrix.xSeries(hCMEC_D3$GSE138309, T)



series <- hCMEC_D3$GSE138309


b <- geneStats.xSeries(hCMEC_D3$GSE138309)
b <- geneStats.xSeries(hCMEC_D3$GSE138309, T)
b <- geneStats.xSeries(hCMEC_D3$GSE138309, annot = FALSE, robust = TRUE)
b <- geneStats.xSeries(hCMEC_D3$GSE138309, T, T)


d1 <- geneStats.xModel(hCMEC_D3)
d1 <- geneStats.xModel(hCMEC_D3, annot = T)
d2 <- geneStats.xModel(hCMEC_D3, descriptive = MEDIAN)
d3 <- geneStats.xModel(hCMEC_D3, descriptive = NWMEAN)




f1 <- keepRuns.xSeries(hCMEC_D3$GSE138309, extra == 1)
f2 <- keepRuns2.xSeries(series_GSE195781, "extra == 1")


  
  # to revert the effect of split(seq(nrow(meta_df))) and get the metadata dataframe back!! 
  unsplit(meta_list, seq(length(meta_list)))
  
  