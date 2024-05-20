# Constructors and methods for `xSeries` and `xModel` S3 classes
#
# Motivation
# ~~~~~~~~~~
# The `xSeries` S3 class attempts to give an effective representation of the
# hierarchical data structure whereby each transcriptomics study typically
# consists of a series of sequencing runs (usually each referring to a different
# RNA sample, coming from different biological sources).
# Built on this, the `xModel` S3 class addresses the need for grouping several
# independent scientific studies together for the purpose of meta-analysis or
# reanalysis.
# In this regard, we introduce the following
#
# Dictionary
# ~~~~~~~~~~
#   Run:     the set of all raw reads resulting from a single sequencing run,
#            referred to by the "Run" ENA accession `(E|D|S)RR[0-9]{6,}`.
#            Basically, each individual FASTQ file, or file pair in the case of
#            non-interleaved PE reads.
#   Sample:  the set of all runs from the same biological RNA sample, referred
#            to by the "BioSample" INSDC accession `SAM(E|D|N)[A-Z]?[0-9]+`, or
#            the "Samples" ENA accession `(E|D|S)RS[0-9]{6,}`, or the
#            corresponding GEO sample (alias) accession `GSM[0-9]+`. Most of the
#            times, runs and samples are the same things, however it could
#            happen that a single RNA sample is sequenced through multiple runs
#            (also referred to as 'technical replicates').
#   Series:  (also referred to as Project (INSDC), or Study (ENA)) the set of
#            all samples pertaining to the same experimental design, across all
#            conditions of interest to that particular scientific study or
#            research project. It can be referred to by the "BioProject" INSDC
#            accession `PRJ(E|D|N)[A-Z][0-9]+`, or the "Study" ENA accession
#            `(E|D|S)RP[0-9]{6,}`, or the corresponding GEO series (alias)
#            accession `GSE[0-9]+`, or the corresponding ArrayExpress experiment
#            (alias) accession `E-[A-Z]{4}-[0-9]+`.
#   Model:   a collection of studies (i.e., a set of Series) from the same
#            biological model.
#
# NOTE: since the ENA run accession is the only ID that is ensured to be unique
#       on a per-file basis (even in the case of technical replicates), it is
#       used here as the base reference for the construction of `xSeries`
#       class objects.
#       Importantly, ENA provides Archive-Generated FASTQs Files for download as
#       an effective  means of bringing consistency to data provided. By
#       imposing a level of uniformity on these files, users know what to expect
#       of them and may incorporate them into pipelines with minimal friction.
#
# Assumptions
# ~~~~~~~~~~~
# A number of (hopefully reasonable) assumptions upon file and data organization
# are made for object construction to be successful:
#
# 1 - Each series (or study) is represented by two CSV or TSV files, namely the
#     actual gene expression data (read counts) and related metadata for all the
#     samples/runs making up the series.
#
# 2 - All series related to the same model are stored in the same directory,
#     here referred to as the 'target directory' for xModel construction, and
#     each model lives in a separate target directory.
#
# 3 - Both data and metadata have file names starting with the same series
#     accession ID (both GEO and ENA are fine), followed by an underscore (`_`)
#     and either the distinctive string `CountMatrix` or `meta` (ignoring case)
#     if they are counts or metadata, respectively. Any additional characters
#     are allowed after this pattern, provided that the terminal file name
#     extension is either CSV or TSV (again ignoring case). Examples of good
#     file name patterns are:
#
#         GSE205739_CountMatrix_genes_TPM.tsv
#         PRJNA847413_countMATRIX.csv
#         GSE76528_meta.csv
#         PRJNA463482_metadata.csv
#
# . - Count tables have proper column names as header, among which a mandatory
#     ID column matching the regex `gene.*id` or `transcript.*id` (usually being
#     `gene_id` or `transcript_id`, respectively).
#
# . - Metadata files have proper column names as table header; among them only
#     `ena_run` is mandatory, as the only ID that is certainly unique to each
#     file.
#
# . - In count table, the names of the columns containing the counts of each
#     different run feature the corresponding `ena_run` ID (and possibly other
#     text strings).
#
# IMPLEMENTATION NOTE
# ~~~~~~~~~~~~~~~~~~~
# To make the code lighter, anonymous functions defined in *apply or pipes are
# NOT self-contained, but instead happily access variables from the outer scope.



# --- Package Dependencies -----------------------------------------------------

library(dplyr)    # `arrange()`, `select()`, `mutate()`
library(rlang)    # Injection operator `!!`
library(magrittr) # For pipe assignment operator %<>% and Aliases (equals())

# --- Internal Functions -------------------------------------------------------

# Automatically adapt file reading to CSV or TSV format
read.xsv <- function(file, header = TRUE) {
  if (grepl(".csv$", file, ignore.case = TRUE)) {
    read.csv(file, header = header)
  } else if (grepl(".tsv$", file, ignore.case = TRUE)) {
    read.delim(file, header = header)
  } else {
    stop("Non-compliant file extension.")
  }
}

# A Java-style binary operator (with 2 aliases) to concatenate strings in pipe
`%|+>%` <- `%+>%` <- `%+%` <- \(x,y)paste0(x,y)

# Takes in a series ID (e.g., PRJNA141411, or GSE29580), a file list (actually
# a character vector), and a vector of two patterns to match (one for the count
# matrix files and the other for metadata).
# Returns FALSE (i.e., do not skip that series) if both files of counts and
# metadata are within the file list. Returns TRUE (i.e., skip that series)
# otherwise.
check_filenames <- function(series_ID, files, pattern) {
  
  # Set "keep it" as default
  skip_this <- FALSE
  
  # Find file pair
  series_ID %+% "_" |> grep(files, value=T) -> run_pair
  # Check them
  if (length(run_pair) != 2) {
    "Wrong number of files in series " %+% series_ID %+% "... skip it!" |> warning()
    skip_this <- TRUE
  } else {
    # Check if the first `run_pair` element matches 'count' pattern AND the
    # second one matches 'meta' pattern (remember files are sorted). NOTE: the
    # logic below may not seem immediately self-evident, but it's fast and,
    # trust me, it works: this 'sapply' will return an identity matrix iff
    # condition is met. Try it.
    sapply(pattern, grepl, run_pair, ignore.case=T) |> equals(diag(2)) |> all() -> matching
    if (not(matching)) {
      "Bad filename pair in series " %+% series_ID %+% "... skip it!" |> warning()
      skip_this <- TRUE
    }
  }
  return(skip_this)
}
  
# --- Constructors -------------------------------------------------------------

# Create a new xSeries
new_xSeries <- function(series_ID, target_dir = ".") {

  # Get all file names from target directory
  list.files(path = target_dir, pattern = "\\.[ct]sv$") |> sort() -> files
  # Patterns to match
  pattern <- c(count = "_countmatrix.*", meta = "_meta.*")
  
  # Load data-metadata pair (also sort metadata by `ena_run`)
  file.path(target_dir,
            series_ID %+% pattern["count"] |> grep(files, ignore.case=T, value=T)
            ) |> read.xsv() -> counts_df
  file.path(target_dir,
            series_ID %+% pattern["meta"]  |> grep(files, ignore.case=T, value=T)
            ) |> read.xsv() |> arrange(ena_run) -> meta_df
  
  # Convert rows to (named) list
  meta_df |> split(seq(nrow(meta_df))) |> setNames(meta_df$ena_run) -> meta_list
  
  # Find gene ID column in `counts_df`
  "gene.*id|transcript.*id|ENSEMBL|ENSEMBLTRANS" |>
    grep(colnames(counts_df), ignore.case=T) -> ids_index
  
  # Build up a `series` object (list)
  lapply(meta_list, function(run) {
    # Look for Run's count data...
    run$ena_run |> grep(colnames(counts_df)) -> count_index
    # ...and add both counts (if present) and IDs to each Run-list as data frame
    counts_df |> select(IDs = !!ids_index, counts = !!count_index) |>
      list(genes=_) |> append(run, values=_)
    }) -> series

  # Add annotation to `series`
  meta_df$ena_run |> paste(collapse = "|") |> grep(colnames(counts_df), invert=T) -> annot_index
  counts_df |> select(!!annot_index) |> list(annotation=_) |> append(series, values=_) -> series
  
  # Any list knows all the names of the elements it contains (under its
  # attribute 'names'), but it doesn't know its own!
  # So, set as attribute for each run its own name (to access them later)
  sapply(names(series), function(name) {
    attr(series[[name]], "own_name") <- name
    return(series[[name]])}) -> series
  
  # Make `series` an S3 object (inheriting from class 'list') and return it
  structure(series, class = c("xSeries", "list"))
}

# Create a new xModel
new_xModel <- function(target_dir = ".") {

  # Get all file names from target directory
  list.files(path = target_dir, pattern = "\\.[ct]sv$") |> sort() -> files
  # Clean series IDs
  files |> sub("_.*$", "", x=_) |> sort() |> unique() -> series_IDs
  if (length(series_IDs) == 0) {
    "Cannot find expression series in " %+% target_dir %+% ". Stop constructor." |> stop()
  }
  
  # Patterns to match
  pattern <- c(count = "_countmatrix.*", meta = "_meta.*")
  
  # Check filenames for series not to include
  to_skip <- vector(mode = "logical", length = 0)
  sapply(series_IDs, check_filenames, files, pattern) -> to_skip
  series_IDs <- series_IDs[not(to_skip)]
  
  # Build up the xModel object
  lapply(series_IDs, new_xSeries, target_dir) |> setNames(series_IDs) -> base_list
  
  # Set as attribute for series its own name (to access them later)
  sapply(names(base_list), function(name) {
    attr(base_list[[name]], "own_name") <- name
    return(base_list[[name]])}) -> base_list
  
  # Make `base_list` an S3 object (inheriting from class 'list') and return it
  structure(base_list, class = c("xModel", "list"))
}



# --- Generics -----------------------------------------------------------------



# --- Methods ------------------------------------------------------------------

# Get the actual sample-size of a given xSeries
size.xSeries <- function(series) {
  series |> sapply(\(run) ncol(run$genes) == 2) |> unlist() |> sum()
}

# Get the read count matrix out of an xSeries object
countMatrix.xSeries <- function(series, annot = FALSE) {
  # Find run elements
  grep("(E|D|S)RR[0-9]{6,}", names(series), ignore.case=T) -> run_index
  
  # Extract counts, restore run ID names, then merge into one data frame
  series[run_index] |> lapply(function(run) {
    counts_df <- run$genes
    colnames(counts_df)[colnames(counts_df) == "counts"] <- run$ena_run
    counts_df
  }) |> Reduce(\(x, y) merge(x, y, by = "IDs", all = TRUE), x=_) -> count_matrix
  
  if (annot) {
    # Get annotation
    annot <- series$annotation
    "gene.*id|transcript.*id|ENSEMBL|ENSEMBLTRANS" |>
      grep(colnames(annot), ignore.case=T) -> ids_index
    count_matrix <- merge(annot, count_matrix,
                          by.x = ids_index, by.y = "IDs", all.y = TRUE)
  }
  return(count_matrix)
}

# Get gene-wise summary stats out of an xSeries object
geneStats.xSeries <- function(series, annot = FALSE, robust = FALSE) {
  
  # Get a log-transformed count matrix
  count_matrix <- countMatrix.xSeries(series, annot = annot)
  count_index <- sapply(count_matrix, is.numeric)
  count_matrix[,count_index] <- log2(count_matrix[,count_index] + 1)
  
  # Compute descriptive statistics and assemble results
  count_matrix[,!count_index, drop = FALSE] -> xSeries_stats
  xSeries_stats$Mean <- rowMeans(count_matrix[,count_index], na.rm=T)
  xSeries_stats$Std_Dev <- apply(count_matrix[,count_index], 1, sd, na.rm=T)
  xSeries_stats |> mutate(SEM = Std_Dev/sqrt(sum(count_index))) -> xSeries_stats
  # Add median and quartiles if robust = TRUE
  if (robust) {
    xSeries_stats$Median <- apply(count_matrix[,count_index], 1, median, na.rm=T)
    xSeries_stats$Q1 <- apply(count_matrix[,count_index], 1, quantile,
                              probs = 0.25, na.rm=T, names = FALSE)
    xSeries_stats$Q3 <- apply(count_matrix[,count_index], 1, quantile,
                              probs = 0.75, na.rm=T, names = FALSE)
    xSeries_stats |> mutate(IQR = Q3-Q1) -> xSeries_stats
  }
  return(xSeries_stats)
}

# Get gene-wise summary stats out of an xModel object
geneStats.xModel <- function(model, descriptive = MEAN, annot = FALSE) {
  
  # Store descriptive stats for each series into one data frame
  model |> lapply(function(series) {
    xSeries_stats <- geneStats.xSeries(series, annot = FALSE, robust = FALSE)
    colnames(xSeries_stats)[-1] <- colnames(xSeries_stats)[-1] %+% "_" %+% attr(series, "own_name")
    xSeries_stats
  }) |> Reduce(\(x, y) merge(x, y, by = 1, all = TRUE), x=_) -> large_stats
  
  # Set the list of effective sample sizes as attribute
  # ...to make them available to descriptive() function
  attr(large_stats, "sample_size") <- sapply(model, size.xSeries)
  
  # Compute model-level descriptive stats
  large_stats |> descriptive() -> xModel_stats
  
  if (annot) {
    # Implement annotation appending
    # Merge annotations from all the series by gene_ids (keeping all), then collapse by ','
    # possible different gene symbols or gene names associated with the same ENSGENE.
    # Finally merge this global annoation with xModel_stats
    #model |> lapply(\(series) series$annotaton) |>
    #  Reduce(\(x, y) merge(x, y, by = 1, all = TRUE), x=_)
  }
  return(xModel_stats)
}

MEAN <- function(large_stats) {
  large_stats |> colnames() |> grep("^Mean_", x=_) -> mean_index
  large_stats[,1, drop = FALSE] -> xModel_stats
  
  xModel_stats$Mean <- rowMeans(large_stats[,mean_index], na.rm=T)
  xModel_stats$Std_Dev <- apply(large_stats[,mean_index], 1, sd, na.rm=T)
  xModel_stats |> mutate(SEM = Std_Dev/sqrt(sum(mean_index)))
}

MEDIAN <- function(large_stats) {
  large_stats |> colnames() |> grep("^Mean_", x=_) -> mean_index
  large_stats[,1, drop = FALSE] -> xModel_stats
  
  xModel_stats$Median <- apply(large_stats[,mean_index], 1, median, na.rm=T)
  xModel_stats$Q1 <- apply(large_stats[,mean_index], 1, quantile,
                           probs = 0.25, na.rm=T, names = FALSE)
  xModel_stats$Q3 <- apply(large_stats[,mean_index], 1, quantile,
                           probs = 0.75, na.rm=T, names = FALSE)
  xModel_stats |> mutate(IQR = Q3-Q1)
}

# Sample size (n) weighted mean
NWMEAN <- function(large_stats) {
  large_stats |> colnames() |> grep("^Mean_", x=_) -> mean_index
  large_stats[,1, drop = FALSE] -> xModel_stats
  
  large_stats |> attr("sample_size") -> N
  large_stats[,mean_index] -> series_mean
  
  # Use mapply to apply grepl element-wise and check correspondence between data
  # and attributes for each series (it should never happen... but just in case)
  mapply(grepl, names(N), colnames(series_mean)) |> all() -> good
  if (!good) {
    stop("Series ID in 'large_stats' colnames do not match 'sample_size' attribute names.")
  }
  
  # Operator %*% is used for matrix multiplication or dot product of two vectors
  xModel_stats$Weighted_Mean <- apply(series_mean, 1, '%*%', N)/sum(N)
  # squared deviations (d2) and Bessel-corrected Weighted SD
  nwmeans <- xModel_stats$Weighted_Mean
  (series_mean - matrix(nwmeans,
                        nrow = length(nwmeans),
                        ncol = length(N), byrow = FALSE))^2 -> d2
  xModel_stats$Weighted_Std_Dev <- (apply(d2, 1, '%*%', N)/(sum(N)-1))^(0.5)
  # NOTE: There is no widely accepted definition of standard error of the
  #       weighted mean. For a discussion about this topic see, e.g.,
  # https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
  return(xModel_stats)
}






# Here `logic` is an unquoted logical expression to be used as filter criterium.
keepRuns.xSeries <- function(series, logic) {
  
  # Capture `logic` expression for Non-Standard Evaluation
  logic_call <- substitute(logic)
  
  # Find runs in `series` that match the `logic` condition
  series |> sapply(function(element){
    if(grepl("(E|D|S)RR[0-9]{6,}", element |> attr("own_name"))) {
      # Evaluate captured expression in the proper environment
      logic_call |> eval(envir = element)
    } else {TRUE}
  }) -> keep_these
  # Restore attributes (subsetting drops all, except names, dim and dimnames)
  # and return
  series[keep_these] |> structure(class = c("xModel", "list"))
}


# Here `logic` is a string, namely the double-quoted logical expression to be
# used as filter criterium.
# NOTE: this is an alternative backup version of the previous method, in case
#       'substitute()' is not successfully processed in some particular settings.
keepRuns2.xSeries <- function(series, logic) {
  
  # Find runs in `series` that match the `logic` condition
  series |> sapply(function(element){
    if(grepl("(E|D|S)RR[0-9]{6,}", element |> attr("own_name"))) {
      # Evaluate the string expression in the proper data environment
      logic |> parse(text=_) |> eval() |> with(element, expr=_)
    } else {TRUE}
  }) -> keep_these
  # Restore attributes (subsetting drops all, except names, dim and dimnames)
  # and return
  series[keep_these] |> structure(class = c("xModel", "list"))
}



