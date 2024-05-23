


library(sloop) # Helpers for 'OOP' in R



ftype(geneStats)
ftype(geneStats.xModel)

geneStats.xModel
s3_get_method(geneStats.xModel) # this works even when method is not exported from a package

methods(geneStats)
s3_methods_class("xSeries")
s3_methods_class("xModel")

s3_methods_generic("countMatrix")
s3_methods_generic("geneStats")








target_dir <- "E:/UniTo Drive/WORKS/0010 - Ongoing/Endothelion/data/in/Lines/hCMEC_D3" # @ HOME
target_dir <- "D:/UniTo Drive/WORKS/0010 - Ongoing/Endothelion/data/in/Lines/hCMEC_D3" # @ DBIOS

series_ID <- "GSE195781"

# Build a xSeries object
series_GSE195781 <- new_xSeries(series_ID, target_dir)

# Build a xModel object
hCMEC_D3 <- new_xModel(target_dir)



s3_dispatch(geneStats(hCMEC_D3))
s3_dispatch(geneStats(hCMEC_D3$GSE138309))

s3_dispatch(`[`(hCMEC_D3))
s3_dispatch(`[`(hCMEC_D3$GSE138309))
s3_dispatch(`[`(hCMEC_D3$GSE138309$SRR10217414))


otype(hCMEC_D3)
otype(hCMEC_D3$GSE138309)
otype(hCMEC_D3$GSE138309$SRR10217414)

class(hCMEC_D3)
class(hCMEC_D3$GSE138309)
class(hCMEC_D3$GSE138309$SRR10217414)

s3_class(hCMEC_D3)
s3_class(hCMEC_D3$GSE138309)
s3_class(hCMEC_D3$GSE138309$SRR10217414)


# Dispatch to N_*.xSeries methods
N_series(hCMEC_D3$GSE138309)
N_selection(hCMEC_D3$GSE138309)
N_genome(hCMEC_D3$GSE138309)

# Dispatch to countMatrix.xSeries
countMatrix(hCMEC_D3$GSE138309) |> head()
countMatrix(hCMEC_D3$GSE138309, annot = T) |> head()

# Dispatch to geneStats.xSeries
geneStats(hCMEC_D3$GSE138309) |> head()
geneStats(hCMEC_D3$GSE138309, T) |> head()
geneStats(hCMEC_D3$GSE138309, annot = FALSE, robust = TRUE) |> head()
geneStats(hCMEC_D3$GSE138309, T, T) |> head()

# Dispatch to geneStats.xModel
geneStats(hCMEC_D3) |> head()
geneStats(hCMEC_D3, annot = T) |> head()
geneStats(hCMEC_D3, descriptive = MEDIAN) |> head()
geneStats(hCMEC_D3, descriptive = NWMEAN) |> head()

# No default defined
geneStats(hCMEC_D3$GSE138309$SRR10217414)


# Dispatch to keepRuns.xSeries
attributes(f1)
f1 <- keepRuns(hCMEC_D3$GSE138309, extra == 1)
attributes(f1)

# Cannot dispatch
f2 <- keepRuns2.xSeries(series_GSE195781, "extra == 1")






attributes(series_GSE195781)
attr(series_GSE195781, "pippo") <- "bestia"
attributes(series_GSE195781)
keepRuns(series_GSE195781, extra == 1) |> View()
keepRuns(series_GSE195781, extra == 1) |> attributes()





cerca tutti [[]] sostituibili con $
  
  %<>%



# to revert the effect of split(seq(nrow(meta_df))) and get the metadata dataframe back!! 
unsplit(meta_list, seq(length(meta_list)))
  
  