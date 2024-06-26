# Constructors and methods for `xSeries` and `xModel` S3 classes
#
# IMPLEMENTATION NOTES
# ~~~~~~~~~~~~~~~~~~~~
# To make the code lighter, anonymous functions defined in *apply or pipes are
# NOT self-contained, but instead happily access variables from the outer scope.
#
# Also elements from GLOBAL list are accessed without being passed as arguments
# (otherwise what global vars would they be?). However, such variables are never
# modified and the (uppercase) `GLOBAL` name should be enough to make their
# global nature explicit.
#
# Default methods are (still?) not implemented.

# --- Package Dependencies -----------------------------------------------------

library(dplyr)    # `arrange()`, `select()`, `filter()`, `mutate()`
library(rlang)    # `sym()` and Injection operator `!!`
library(magrittr) # For pipe assignment operator %<>% and Aliases (equals())
library(tidyr)    # `separate_rows()`

# --- Globals ------------------------------------------------------------------

GLOBAL <-
  list(filename = list(counts = "_countmatrix.*",
                       metadata = "_meta.*"),
       geneID_regex = "^IDs$|gene.*id|transcript.*id|ENSEMBL|ENSEMBLTRANS",
       run_regex = "(E|D|S)RR[0-9]{6,}")

# --- Internal Functions -------------------------------------------------------

# A Java-style binary operator to concatenate strings in pipe.
`%+%` <- \(x,y){paste0(x,y)}

# A Java-style binary operator to construct paths in a platform-independent way.
`%//%` <- \(x,y){file.path(x, y, fsep = .Platform$file.sep)}

# Tab Stop (borrowed from r4tcpl package)
tab <- function (word = "", sp = 7) {
  word <- toString(word)
  if (sp < nchar(word)) {
    warning("The word to include exceeds tab stop... can't align properly!")
  } else {
    word <- paste0(word, paste(rep(" ", sp - nchar(word)), collapse = ""))
  }
  return(word)
}

# Shorthand for case insensitive 'grep' and 'grepl'
grepi <- function(pattern, x, ...) {grep(pattern, x, ignore.case = TRUE, ...)}
grepli <- function(pattern, x, ...) {grepl(pattern, x, ignore.case = TRUE, ...)}

# Only called by constructors. Finds all the CSV and TSV files within the target
# directory. If any, performs a filename validity check and returns a character
# vector of suitable files. Stops execution (i.e., object construction) if none
# is found matching the required filename pattern. 
get_files <- function(target_dir) {
  # Scan target directory looking for valid filenames
  target_dir |> list.files(pattern = "\\.[ct]sv$", ignore.case=T) -> files
  GLOBAL$filename |> paste(collapse = "|") |>
    grepi(files, value=T) |> sort() -> files
  if (length(files) == 0) {
    "Can't find suitable CSV/TSV files in " %+% target_dir |> stop()
  }
  return(files)
}

# Reads a table from file, automatically adapting to CSV or TSV formats.
read.xsv <- function(file, header = TRUE, ...) {
  if (grepli(".csv$", file)) {read.csv(file, header = header, ...)}
  else if (grepli(".tsv$", file)) {read.delim(file, header = header, ...)}
  else {stop("Non-compliant file extension.")}
}

# Any named list knows the names of all the elements it contains (under its
# 'names' attribute), but it doesn't know its own (even when it has one, e.g.,
# because it is itself an element of a named list)! So, this function sets as
# attribute for each element of a named list its own name (to access it later).
set_own_names <- function(parent_list) {
  names(parent_list) |> lapply(function(element_name) {
    attr(parent_list[[element_name]], "own_name") <- element_name
    return(parent_list[[element_name]])
  }) |> setNames(names(parent_list)) # Also keep original names in parent list
}

# Takes in a series ID (e.g., PRJNA141411, or GSE29580), a file list (actually
# a character vector), and a vector or list of two patterns to match (one for
# the count matrix files and the other for metadata). Returns FALSE (i.e., do
# not skip that series) if both files of counts and metadata are within the file
# list. Returns TRUE (i.e., skip that series) otherwise.
check_pairing <- function(series_ID, files, pattern) {
  # Set "keep it" as default
  skip_this <- FALSE
  # Find file pair
  series_ID %+% "_" |> grepi(files, value=T) -> file_pair
  # Check them
  if (length(file_pair) != 2) {
    "Wrong number of files in series " %+% series_ID %+% ": skip!" |> warning()
    skip_this <- TRUE
  } else {
    # Check if the first `file_pair` element matches 'count' pattern AND the
    # second one matches 'meta' pattern (remember files are sorted). NOTE: the
    # logic below may not seem immediately self-evident, but it's fast and,
    # trust me, it works: this 'sapply' will return a 2-by-2 identity matrix iff
    # condition is met. Try it.
    pattern |> sapply(grepli, file_pair) |> equals(diag(2)) |> all() -> matching
    if (not(matching)) {
      "Bad filename pair in series " %+% series_ID %+% ": skip!" |> warning()
      skip_this <- TRUE
    }
  }
  return(skip_this)
}

# Utility function to spot Run objects
is_run <- function(element) {
  element |> attr("own_name") -> element_name
  if (element_name |> is.null()) {
    FALSE
  } else {
    GLOBAL$run_regex |> grepli(element_name)
  }
}

# Utility function that removes one element by name from a list, suite to work
# in pipe and preserve class attribute if applied to xSeries/xModel objects.
remove_element <- function(target_list, element) {
  target_list[!(names(target_list) == element)]
}

# --- Constructors -------------------------------------------------------------

# Create a new xSeries
new_xSeries <- function(series_ID, target_dir = ".") {

  # Get all CSV and TSV files from target directory and check pairing
  target_dir |> get_files() -> files
  series_ID |> check_pairing(files, GLOBAL$filename) |> ifelse(return(), NA) |>
    invisible()
  
  # Load data-metadata pair (also sort metadata by `ena_run`)
  series_ID %+% GLOBAL$filename$counts |> grepi(files, value=T) -> count_file
  target_dir %//% count_file |> read.xsv() -> counts_df
  
  series_ID %+% GLOBAL$filename$metadata |> grepi(files, value=T) -> meta_file
  target_dir %//% meta_file |> read.xsv() |> arrange(ena_run) -> meta_df
  
  # From two data frames to one xSeries object
  to_xSeries(counts_df, meta_df)
}

# Create an xSeries object from R variables
to_xSeries <- function(counts_df, meta_df) {
  # Convert rows to a (named) list
  meta_df |> split(seq(nrow(meta_df))) |> setNames(meta_df$ena_run) -> series
  
  # Find gene ID column in `counts_df`
  GLOBAL$geneID_regex |> grepi(colnames(counts_df), value=T) -> ids_header
  if (length(ids_header) != 1) {
    "Cannot find a suitable gene ID column in `counts_df`" |> stop()
  }
  
  # Add `genes` information to each Run in `series`
  series %<>% lapply(function(run) {
    # Look for Run's count data...
    # (search for "isolated" Run IDs, not to confuse, e.g., SRR123 with SRR1234)
    "(^|[^a-zA-Z0-9])" %+% run$ena_run %+% "($|[^a-zA-Z0-9])" |>
      grepi(colnames(counts_df), value=T) -> count_header
    # ...and add both counts (if present) and IDs to each Run, as data frame
    counts_df |> select(IDs = !!ids_header, counts = !!count_header) |>
      list(genes=_) |> append(run, values=_)
  })
  
  # Add annotation to `series` (and force IDs to first position)
  GLOBAL$run_regex |> grepi(colnames(counts_df), invert=T) -> annot_index
  counts_df |> select(!!annot_index) |>
    select(IDs = !!ids_header, everything()) |>
    list(annotation=_) |> append(series, values=_) -> series
  
  # Assign to each Run its own name
  series %<>% set_own_names()
  
  # Make `series` an S3 object (inheriting from class 'list') and return it
  structure(series, class = c("xSeries", "list"))
}

# Create a new xModel
new_xModel <- function(target_dir = ".") {

  # Get all CSV and TSV files from target directory
  target_dir |> get_files() -> files
  
  # Clean series IDs
  files |> sub("_.*$", "", x=_) |> sort() |> unique() -> series_IDs
  
  # Check filename pairing for series not to include
  to_skip <- vector(mode = "logical", length = 0)
  series_IDs |> sapply(check_pairing, files, GLOBAL$filename) -> to_skip
  series_IDs <- series_IDs[not(to_skip)]
  
  # Build up the xModel object
  series_IDs |> lapply(new_xSeries, target_dir) |> setNames(series_IDs) -> model
  
  # Assign to each Series its own name
  model %<>% set_own_names()
  # Also give a name to the whole model
  attr(model, "own_name") <- basename(target_dir)
  
  # Make `model` an S3 object (inheriting from class 'list') and return it
  structure(model, class = c("xModel", "list"))
}

# --- [ ... ] ------------------------------------------------------------------

# Special subsetting behavior for objects that preserves their attributes, given
# that standard subsetting drops all attributes except names, dim and dimnames.
# Obviously, in this case, the generic already exists, so we only need to
# provide the method.
`[.xSeries` <- function(series, i, ...) {
  # Store original attributes (for later restoring) and update names
  series |> attributes() -> bkp_attribs
  if(!is.null(bkp_attribs$names)) {
    bkp_attribs$names <- names(series)[i]
  }
  # Standard subset of list
  out <- unclass(series)[i]
  # Restore attributes (including 'class')
  attributes(out) <- bkp_attribs
  return(out)
}

`[.xModel` <- function(model, i, ...) {
  # Store original attributes (for later restoring) and update names
  model |> attributes() -> bkp_attribs
  if(!is.null(bkp_attribs$names)) {
    bkp_attribs$names <- names(model)[i]
  }
  # Standard subset of list
  out <- unclass(model)[i]
  # Restore attributes (including 'class')
  attributes(out) <- bkp_attribs
  return(out)
}

# --- N_series -----------------------------------------------------------------

N_series <- function(xSeries) {
  UseMethod("N_series")
}

# Returns the size of the complete Series
# (number of Runs in the metadata table)
N_series.xSeries <- function(series) {
  names(series) |> grepi(GLOBAL$run_regex, x=_) |> length()
}

# --- N_selection --------------------------------------------------------------

N_selection <- function(xSeries) {
  UseMethod("N_selection")
}

# Returns the actual size of the sub-Series of interest
# (number of Runs in the count matrix)
N_selection.xSeries <- function(series) {
  series |> sapply(\(run){ncol(run$genes) == 2}) |> unlist() |> sum()
}

# --- N_genome -----------------------------------------------------------------

N_genome <- function(xSeries) {
  UseMethod("N_genome")
}

# Returns the size of the genome screened within a given Series
N_genome.xSeries <- function(series){series$annotation |> nrow()}

# --- totalCounts --------------------------------------------------------------

totalCounts <- function(xSeries) {
  UseMethod("totalCounts")
}

# Return the sum of counts for each Run in a Series (useful to guess the metric
# and distinguish among raw counts, normalized counts, FPKM, TPM ...) 
totalCounts.xSeries <- function(series) {
  # Find Run elements and get their count sum
  GLOBAL$run_regex |> grepi(names(series)) -> run_index
  series[run_index] |> sapply(\(run) run$genes$counts |> sum())
}

# --- factTable ----------------------------------------------------------------

factTable <- function(xObject) {
  UseMethod("factTable")
}

# Print main facts about an xSeries object structure and contents
factTable.xSeries <- function(series) {
  cat(tab("xSeries", 8), ":", attr(series, "own_name"), "\n")
  cat(tab("Contents", 8), ":", N_series(series), "Runs + annotation\n")
  series |> sapply(\(element) {
    if(element |> is_run()) {
      cat(" - ", tab(attr(element, "own_name"), 12),
          "[ ", nrow(element$genes), " genes | ",
          tab(sprintf("%.2e", sum(element$genes$counts)), 9),
          "TOT counts ]\n", sep = "")
    } else {
      cat(" + ", tab(attr(element, "own_name"), 12),
          "[ ", nrow(element), " genes ]\n", sep = "")
      names(element) |> sapply(\(field) cat("    -", field, "\n"))
    }
  }) |> invisible()
}

# Print main facts about an xModel object structure and contents
factTable.xModel <- function(model) {
  cat(tab("xModel", 8), ":", attr(model, "own_name"), "\n")
  cat(tab("Contents", 8), ":", length(model), "xSeries\n")
  model |> sapply(\(series) {
    cat(" - ", tab(attr(series, "own_name"), 12),
        "[ ", tab(N_series(series), 3), "Runs | ",
        N_genome(series), " genes ]\n", sep = "")
  }) |> invisible()
}

# --- metaTable ----------------------------------------------------------------

metaTable <- function(xSeries) {
  UseMethod("metaTable")
}

# Rebuild the matadata table out of an xSeries object
metaTable.xSeries <- function(series) {
  series |> remove_element("annotation") |> lapply(remove_element, "genes") |>
    lapply(as.data.frame) |> lapply(\(run) {
      attr(run, "row.names") <- run$ena_run
      return(run)
    }) |> unsplit(seq(N_series(series)))
}

# --- countMatrix --------------------------------------------------------------

countMatrix <- function(xSeries, annot) {
  UseMethod("countMatrix")
}

# Rebuilds the read count matrix out of an xSeries object
countMatrix.xSeries <- function(series, annot = FALSE) {
  # Find Run elements
  GLOBAL$run_regex |> grepi(names(series)) -> run_index
  
  # Extract counts, set clean Run IDs as column names, merge into one data frame
  series[run_index] |> lapply(function(run) {
    counts_df <- run$genes
    colnames(counts_df)[colnames(counts_df) == "counts"] <- run$ena_run
    return(counts_df)
  }) |> Reduce(\(x,y){merge(x, y, by = "IDs", all = TRUE)}, x=_) -> count_matrix
  
  if (annot) {
    # Merge annotation with counts
    count_matrix <- merge(series$annotation, count_matrix, by = "IDs", all.y=T)
  }
  return(count_matrix)
}

# --- geneStats ----------------------------------------------------------------

geneStats <- function(xObject, ...) {
  UseMethod("geneStats")
}

# Computes gene-wise summary stats out of an xSeries object
# Value: returns a data frame with the extra-attribute 'metadata' about the Runs
#        actually present in the count matrix used for stats (i.e, selection).
#        Gene expression stats are always in log scale.
geneStats.xSeries <- function(series, annot = FALSE, robust = FALSE) {
  
  # Get a log-transformed count matrix
  series |> countMatrix(annot = annot) -> count_matrix
  count_matrix |> sapply(is.numeric) -> count_index
  count_matrix[,count_index] <- log2(count_matrix[,count_index] + 1)
  
  # Compute Series-level descriptive stats
  count_matrix[,!count_index, drop = FALSE] -> xSeries_stats
  count_matrix[,count_index] -> x
  x |> rowMeans(na.rm=T) -> xSeries_stats$Mean
  x |> apply(1, sd, na.rm=T) -> xSeries_stats$Std_Dev
  xSeries_stats %<>% mutate(SEM = Std_Dev/sqrt(sum(count_index)))
  
  # Add median and quartiles if robust = TRUE
  if (robust) {
    x |> apply(1, median, na.rm=T) -> xSeries_stats$Median
    x |> apply(1, quantile, probs = 0.25, na.rm=T, names=F) -> xSeries_stats$Q1
    x |> apply(1, quantile, probs = 0.75, na.rm=T, names=F) -> xSeries_stats$Q3
    xSeries_stats %<>% mutate(IQR = Q3-Q1)
  }
  
  # Add selection metadata as attribute
  metaTable(series)$ena_run %in% names(count_matrix) -> meta_index
  attr(xSeries_stats, "metadata") <- metaTable(series)[meta_index,]
  
  return(xSeries_stats)
}

# Computes gene-wise summary stats out of an xModel object
# Meta-analysis Inclusion Criteria (maic)
#   inclusive -> keep ALL genes, even if not present in every Series (union)
#   exclusive -> keep only genes that present in every Series (intersection)
geneStats.xModel <- function(model, descriptive = MEAN,
                             maic = "inclusive", annot = FALSE) {
  
  # Check if screened genomes are all equal across different Series
  model |> sapply(function(series) {
    setequal(series[[1]]$genes$IDs, model[[1]][[1]]$genes$IDs)
  }) |> all() -> equal_genomes
  if (not(equal_genomes)) {
    warning("xModel with different genome sizes." %+%
            "\n`maic` parameter is critical in determining final genome size.")
  }
  # Store descriptive stats for each Series into one data frame
  model |> lapply(function(series) {
    series |> geneStats.xSeries(annot = FALSE, robust = FALSE) -> xSeries_stats
    colnames(xSeries_stats)[-1] %+%
      "_" %+% attr(series, "own_name") -> colnames(xSeries_stats)[-1]
    return(xSeries_stats)
  }) |> Reduce(\(x,y) merge(x, y, by = "IDs",
                            all = ifelse(maic == "inclusive", T, F)),
               x=_) -> large_stats
  
  # Set the list of the actual number of Runs per Series as a Model attribute
  # ...to make them available to descriptive() function
  model |> sapply(N_selection) -> attr(large_stats, "selection_size")
  
  # Compute Model-level descriptive stats
  large_stats |> descriptive() -> xModel_stats
  
  # Model-level annotation synthesis
  if (annot) {
    warning("Model-level annotation synthesis coming soon...")
  }
  return(xModel_stats)
}

MEAN <- function(large_stats) {
  large_stats |> colnames() |> grepl("^Mean_", x=_) -> mean_index
  large_stats[,"IDs", drop = FALSE] -> xModel_stats
  large_stats[,mean_index] -> x
  x |> rowMeans(na.rm=T) -> xModel_stats$Mean
  x |> apply(1, sd, na.rm=T) -> xModel_stats$Std_Dev
  xModel_stats |> mutate(SEM = Std_Dev/sqrt(sum(mean_index)))
}

MEDIAN <- function(large_stats) {
  large_stats |> colnames() |> grepl("^Mean_", x=_) -> mean_index
  large_stats[,"IDs", drop = FALSE] -> xModel_stats
  large_stats[,mean_index] -> x
  x |> apply(1, median, na.rm=T) -> xModel_stats$Median
  x |> apply(1, quantile, probs = 0.25, na.rm=T, names=F) -> xModel_stats$Q1
  x |> apply(1, quantile, probs = 0.75, na.rm=T, names=F) -> xModel_stats$Q3
  xModel_stats |> mutate(IQR = Q3-Q1)
}

# Sample size (N) weighted mean
NWMEAN <- function(large_stats) {
  large_stats |> colnames() |> grepl("^Mean_", x=_) -> mean_index
  large_stats[,"IDs", drop = FALSE] -> xModel_stats
  large_stats[,mean_index] -> x
  large_stats |> attr("selection_size") -> N
  
  # Use mapply to apply grepl element-wise and check correspondence between data
  # and attributes for each series (it should always match... but just in case)
  mapply(grepl, names(N), colnames(x)) |> all() -> good
  if (!good) {
    stop("Internal error:\nSeries IDs in 'large_stats' colnames" %+%
           " do not match 'selection_size' attribute names.")
  }
  # Operator %*% is used for matrix multiplication or dot product of two vectors
  x |> apply(1, '%*%', N)/sum(N) -> xModel_stats$Weighted_Mean
  # squared deviations (d2) and Bessel-corrected Weighted SD
  xModel_stats$Weighted_Mean -> nwmeans
  (x - matrix(nwmeans, nrow = length(nwmeans), ncol = ncol(x)))^2 -> d2
  (d2 |> apply(1, '%*%', N)/(sum(N)-1))^(0.5) -> xModel_stats$Weighted_Std_Dev
  # NOTE: There is no widely accepted definition of standard error of the
  #       weighted mean. For a discussion about this topic see, e.g.,
  # https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
  return(xModel_stats)
}

# --- pruneRuns ----------------------------------------------------------------

pruneRuns <- function(xObject) {
  UseMethod("pruneRuns")
}

# Remove all Runs with no count data from an xSeries object
pruneRuns.xSeries <- function(series) {
  # Find Runs in `series` that have the `count` column in gene data frame
  series |> sapply(function(element) {
    if(element |> is_run()) {
      ncol(element$genes) == 2
    } else {TRUE}
  }) -> keep_these
  # Dispatch substetting to `[.xSeries`
  return(series[keep_these])
}

# Remove all Runs with no count data from each xSeries of an xModel object
pruneRuns.xModel <- function(model) {
  # Store original attributes (for later restoring)
  model |> attributes() -> bkp_attribs
  model |> lapply(pruneRuns.xSeries) -> out
  # Restore attributes (including 'class')
  attributes(out) <- bkp_attribs
  return(out)
}

# --- keepRuns -----------------------------------------------------------------

keepRuns <- function(xObject, logic) {
  UseMethod("keepRuns")
}

# Filter an xSeries object by selecting Runs that match a logical condition.
# Here `logic` is a string, namely the quoted logical expression to be used as
# filter criterium.
keepRuns.xSeries <- function(series, logic) {
  # Find Runs in `series` that match the `logic` condition
  series |> sapply(function(element) {
    if(element |> is_run()) {
      # Evaluate the string expression in the proper data environment
      logic |> parse(text=_) |> eval() |> with(element, expr=_)
    } else {TRUE}
  }) -> keep_these
  # Dispatch substetting to `[.xSeries`
  return(series[keep_these])
}

# NOTE: Currently with no generic!
# This is an alternative version of the previous method, the only difference
# being that this method uses Non-Standard Evaluation (NSE) to take an unquoted
# logical expression as filter criterium. Although working and more elegant, NSE
# complicated the implementation of this method a bit and I couldn't get it
# working properly when called from `keepRuns.xModel`.
keepRuns2.xSeries <- function(series, logic) {
  
  # Capture `logic` expression for Non-Standard Evaluation
  logic_call <- substitute(logic)
  
  # Find Runs in `series` that match the `logic` condition
  series |> sapply(function(element) {
    if(element |> is_run()) {
      # Evaluate captured expression in the proper environment
      logic_call |> eval(envir = element)
    } else {TRUE}
  }) -> keep_these
  # Dispatch substetting to `[.xSeries`
  return(series[keep_these])
}

# Filter an xModel object by selecting Runs that match a logical condition in
# each one of its xSeries elements.
# Here `logic` is a string, namely the quoted logical expression to be used as
# filter criterium.
keepRuns.xModel <- function(model, logic) {
  # Store original attributes (for later restoring)
  model |> attributes() -> bkp_attribs
  model |> lapply(keepRuns.xSeries, logic) -> out
  # Restore attributes (including 'class')
  attributes(out) <- bkp_attribs
  return(out)
}

# --- subsetGenes --------------------------------------------------------------

subsetGenes <- function(xObject, key, geneset) {
  UseMethod("subsetGenes")
}

# Filter all the elements of an xSeries (both Runs and annotation) based on a
# key and a list of genes of interest (GOIs). Actually, `geneset` can be a
# character vector, a n-by-1 matrix, or a n-by-1 dataframe.
subsetGenes.xSeries <- function(series, key, geneset) {
  # Get count matrix (with annotation), explode possibly collapsed entries
  # (e.g., many SYMBOLS per ENSG ID), then filter row-wise.
  geneset %<>% unlist()
  series |> countMatrix(annot = TRUE) |> separate_rows(!!key, sep = ",") |>
    as.data.frame() |> filter(!!sym(key) %in% geneset) -> sub_counts_df
  
  # Check for completeness
  if (setdiff(geneset, sub_counts_df[,key]) |> length() > 0) {
    cat("\nWARNING by ", attr(series, "own_name"), ":",
        "\n Can't find these Genes of Interest in count matrix:\n  ", sep = "")
    cat(setdiff(geneset, sub_counts_df[,key]), sep = "\n  ")
  }
  # Check for duplicate entries in `key` column (e.g., many ENSG IDs per SYMBOL)
  if (sub_counts_df[key] |> duplicated() |> sum() > 0) {
    "Duplicated '" %+% key %+% "' entries detected but not handled" |> warning()
  }

  # Reassemble the xSeries from filtered data
  to_xSeries(counts_df = sub_counts_df,
             meta_df = metaTable(series))
}

# Filter all the elements of all the xSeries of an xModel object.
subsetGenes.xModel <- function(model, key, geneset) {
  # Store original attributes (for later restoring)
  model |> attributes() -> bkp_attribs
  model |> lapply(subsetGenes.xSeries, key, geneset) -> out
  # Re-assign to each Run its own name
  out %<>% set_own_names()
  # Restore attributes (including 'class')
  attributes(out) <- bkp_attribs
  return(out)
}


