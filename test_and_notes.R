# This test is designed to be performed in interactive mode with manual
# inspection of the output. However, it shouldn't be too hard to automate it.

# --- Load packages ------------------------------------------------------------

library(httr)    # Tools for working with URLs and HTTP
library(stringr) # Common string operations
library(sloop)   # Helpers for 'OOP' in R

source("./seq_loader.R")

# --- Prepare the test file set ------------------------------------------------

# Check the operating system to set the right target_dir
os_type <- Sys.info()["sysname"]

if (os_type == "Windows") {
  target_dir <- Sys.getenv("USERPROFILE") %//% "Desktop/hCMEC_D3" 
} else if (os_type == "Linux") {
  target_dir <- Sys.getenv("HOME") %//% "hCMEC_D3" 
} else {
  stop("Running on an unknown OS\n")
}

if (!dir.exists(target_dir)) {dir.create(target_dir, recursive = TRUE)}

# Vector of URLs
URLs <- c(
  "https://zenodo.org/records/10960146/files/GSE138309_CountMatrix_genes_TPM.tsv",
  "https://zenodo.org/records/10960146/files/GSE139133_CountMatrix_genes_TPM.tsv",
  "https://zenodo.org/records/10960146/files/GSE195781_CountMatrix_genes_TPM.tsv",
  "https://zenodo.org/records/10960146/files/GSE205739_CountMatrix_genes_TPM.tsv",
  "https://zenodo.org/records/10960146/files/GSE76528_CountMatrix_genes_TPM.tsv",
  "https://zenodo.org/records/10960146/files/GSE138309_meta.csv",
  "https://zenodo.org/records/10960146/files/GSE139133_meta.csv",
  "https://zenodo.org/records/10960146/files/GSE195781_meta.csv",
  "https://zenodo.org/records/10960146/files/GSE205739_meta.csv",
  "https://zenodo.org/records/10960146/files/GSE76528_meta.csv")
# Download files
URLs |> lapply(\(URL) {
  target_dir %//% basename(URL) -> dest_file
  "Downloading: " %+% dest_file |> print(quote = FALSE)
  URL |> httr::GET(httr::write_disk(dest_file, overwrite = TRUE))
}) -> download_status

# Selective file remodeling to introduce some exceptions
# Introduce harmless modifications to a metadata filename
file.rename(target_dir %//% "GSE76528_meta.csv",
            target_dir %//% "GSE76528_MeTaDATA_xxx.csv")
# Introduce breaking modifications to a data filename (make pattern invalid)
file.rename(target_dir %//% "GSE195781_CountMatrix_genes_TPM.tsv",
            target_dir %//% "GSE195781_CounMatrix_genes_TPM.tsv")
# Simulate two metadata files (and no data) for a single Series
file.rename(target_dir %//% "GSE139133_CountMatrix_genes_TPM.tsv",
            target_dir %//% "GSE139133_metadata2.csv")
# Push gene ID column forward (to third position) in a Series count matrix
read.delim(target_dir %//% "GSE205739_CountMatrix_genes_TPM.tsv") |> 
  select(c(SYMBOL, GENENAME), everything()) |>
  write.table(target_dir %//% "GSE205739_CountMatrix_genes_TPM.tsv",
              quote = FALSE, sep = "\t", eol = "\n",
              row.names = FALSE, col.names = TRUE)
# Assign a bad quality score (extra == 0) to a Run of a Series
read.csv(target_dir %//% "GSE138309_meta.csv") |>
  mutate(extra = replace(extra, 6, 0)) |>
  write.table(target_dir %//% "GSE138309_meta.csv",
              quote = TRUE, sep = ",", eol = "\n",
              row.names = FALSE, col.names = TRUE)

# --- Let's test ---------------------------------------------------------------

# Determine function type
ftype(new_xModel)
ftype(geneStats)
ftype(geneStats.xModel)

# Show method definition (even when it is not exported from a package)
geneStats.xModel
s3_get_method(geneStats.xModel)

# List methods of a class
s3_methods_class("xSeries")
s3_methods_class("xModel")

# List methods for a generic
s3_methods_generic("countMatrix")
s3_methods_generic("geneStats")
methods(geneStats)

# Check `get_files()` behavior
# only 'GSE195781_CounMatrix_genes_TPM.tsv' should be excluded
target_dir |> get_files()

# Build an xSeries object
GSE138309 <- new_xSeries("GSE138309", target_dir)
# Explore Series (e.g., find megareads for a given Run)
GSE138309$SRR10217414$read_count/1e6

# Build a new xModel object and test `check_pairing()` behavior
# Two series should be excluded with the following warnings:
# Warning messages:
#  1: In FUN(X[[i]], ...) : Bad filename pair in series GSE139133: skip!
#  2: In FUN(X[[i]], ...) : Wrong number of files in series GSE195781: skip!
hCMEC_D3 <- new_xModel(target_dir)

# Check reordering/renaming of the gene ID column in `annotation` data.frame
read.xsv(target_dir %//% "GSE205739_CountMatrix_genes_TPM.tsv") |> colnames()
hCMEC_D3$GSE205739$annotation |> colnames()

# Test method dispatch
s3_dispatch(geneStats(hCMEC_D3))
s3_dispatch(geneStats(hCMEC_D3$GSE138309))

s3_dispatch(`[`(hCMEC_D3))
s3_dispatch(`[`(hCMEC_D3$GSE138309))
s3_dispatch(`[`(hCMEC_D3$GSE138309$SRR10217414))

# Check object and class types
otype(hCMEC_D3)
otype(hCMEC_D3$GSE138309)
otype(hCMEC_D3$GSE138309$SRR10217414)

class(hCMEC_D3)
class(hCMEC_D3$GSE138309)
class(hCMEC_D3$GSE138309$SRR10217414)

s3_class(hCMEC_D3)
s3_class(hCMEC_D3$GSE138309)
s3_class(hCMEC_D3$GSE138309$SRR10217414)

# ~~~~~~~~~~~~~~ #
#  Test methods  #
# ~~~~~~~~~~~~~~ #

# Dispatch to N_*.xSeries methods
N_series(hCMEC_D3$GSE138309)
N_selection(hCMEC_D3$GSE138309)
N_genome(hCMEC_D3$GSE138309)

# Dispatch to totalCounts.xSeries
totalCounts(hCMEC_D3$GSE138309)
hCMEC_D3 |> lapply(totalCounts)

# Dispatch to metaTable.xSeries
metaTable(hCMEC_D3$GSE138309)
hCMEC_D3 |> lapply(metaTable)

# Dispatch to countMatrix.xSeries
countMatrix(hCMEC_D3$GSE138309) |> head()
countMatrix(hCMEC_D3$GSE138309, annot = TRUE) |> head()

# Disassemble and reassemble an xSeries and check that it matches the original
to_xSeries(counts_df = countMatrix(hCMEC_D3$GSE138309, annot = TRUE),
           meta_df = metaTable(hCMEC_D3$GSE138309)) |> View()

# Dispatch to geneStats.xSeries
geneStats(hCMEC_D3$GSE138309) |> head()
geneStats(hCMEC_D3$GSE138309, T) |> head()
geneStats(hCMEC_D3$GSE138309, annot = FALSE, robust = TRUE) |> head()
geneStats(hCMEC_D3$GSE138309, T, T) |> head()

# Dispatch to geneStats.xModel
geneStats(hCMEC_D3) |> head()
geneStats(hCMEC_D3, annot = TRUE) |> head()
geneStats(hCMEC_D3, descriptive = MEDIAN) |> head()
geneStats(hCMEC_D3, descriptive = NWMEAN) |> head()

# No default defined
geneStats(hCMEC_D3$GSE138309$SRR10217414)

# Dispatch to pruneRuns.xSeries and check attribute preservation (by dispatching
# to extraction operator [...])
hCMEC_D3$GSE138309 |> attributes()
GSE138309_pruned <- pruneRuns(hCMEC_D3$GSE138309)
attributes(GSE138309_pruned)
# Notice that
GSE138309 |> metaTable()
GSE138309_pruned |> metaTable()

# Dispatch to pruneRuns.xModel and check attribute preservation
attributes(hCMEC_D3)
hCMEC_D3_pruned <- pruneRuns(hCMEC_D3)
attributes(hCMEC_D3_pruned)
attributes(hCMEC_D3_pruned$GSE138309)

# Dispatch to keepRuns.xSeries and check attribute preservation (by dispatching
# to extraction operator [...])
GSE138309_pruned |> attributes()
GSE138309_reduced <- keepRuns(GSE138309_pruned, "extra == 1")
attributes(GSE138309_reduced)

# Alternative NSE version (directly applied to unpruned series)
GSE138309_reduced2 <- keepRuns2.xSeries(hCMEC_D3$GSE138309, extra == 1)

# Dispatch to keepRuns.xModel and check attribute preservation
attributes(hCMEC_D3_pruned)
hCMEC_D3_reduced <- keepRuns(hCMEC_D3_pruned, "extra == 1")
attributes(hCMEC_D3_reduced)
attributes(hCMEC_D3_reduced$GSE138309)

# Dispatch to subsetGenes.xSeries and check attribute preservation
GSE138309 |> attributes()
gois <- c("CFTR", "AQP1", "TRPV1", "RYR1" ,"ADORA1")
GSE138309_filtered <- subsetGenes(hCMEC_D3$GSE138309, "SYMBOL", gois)
GSE138309_filtered |> attributes()

GSE138309 |> N_genome()
GSE138309_filtered |> N_genome()

# Dispatch to keepRuns.xModel and check attribute preservation
attributes(hCMEC_D3)
hCMEC_D3_filtered <- subsetGenes(hCMEC_D3, "SYMBOL", gois)
attributes(hCMEC_D3_filtered)
attributes(hCMEC_D3$GSE138309)
attributes(hCMEC_D3_filtered$GSE138309)


