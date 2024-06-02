# SeqLoader Package
###### Constructors and methods for `xSeries` and `xModel` S3 classes

## Motivation
The S3 `xSeries` class attempts to provide a convenient R representation of a
typical RNA-Seq output consisting of a series of sequencing runs, usually each
of them referring to a different RNA sample, from a different biological source.
Building on this, the S3 `xModel` class addresses the need to group several
independent scientific studies for the purpose of meta-analysis or reanalysis.
To describe this hierarchical data structure, we will use a terminology
primarily modeled after the one adopted by the partners of the International
Nucleotide Sequence Database Collaboration (INSDC) in the context of the
Sequence Read Archive (SRA).
See [here](https://feat-fear.github.io/logSpace/docs/genomics/INSDC/) for a
friendly explanation of the INSDC and the related SRA data model.

## Hierarchy and Terminology
1. ___Run___: each sequencing run directly generating a single FASTQ file (or
   file pair in the case of non-interleaved PE reads). Within the INSDC
   framework, _Run metadata objects_ are referred to by an accession following
   the pattern `(E|D|S)RR[0-9]{6,}`.
1. ___Sample___: the set of all _Runs_ from the same biological RNA sample.
   The corresponding INSDC metadata objects are _Samples_ (accession
   `(E|D|S)RS[0-9]{6,}`) or, equivalently, _BioSamples_ (accession
   `SAM(E|D|N)[A-Z]?[0-9]+`). Even _GEO Sample_ IDs (accession `GSM[0-9]+`) can
   be found if data were originally brokered to INSDC’s SRA by GEO.
1. ___Series___: the set of all _Samples_/_Runs_ pertaining to a given
   biological research project or experimental design, across all conditions of
   interest to that particular study. This is roughly what INSDC refers to as
   _Study_ (accession `(E|D|S)RP[0-9]{6,}`) or, equivalently, _BioProject_
   (accession `PRJ(E|D|N)[A-Z][0-9]+`). If data were originally brokered by GEO
   or ArrayExpress, _GEO Series_ IDs (accession `GSE[0-9]+`) or _AE Experiment_
   IDs (accession `E-[A-Z]{4}-[0-9]+`) will be also present in SRA DBs and can
   be used as _Study_ alias accessions, respectively.
1. ___Model___: a collection of _Series_ dealing with the same biological model.
   Since this object is characteristic of the meta-analytic level, which by
   definition involves multiple independent studies, it has no counterpart in
   INSDC.

> [!NOTE]
> Most of the times, _Runs_ and _Samples_ are the same thing, however it could
> happen that a single RNA _Sample_ is sequenced through multiple _Runs_ (also
> referred to as 'technical replicates'). Thus, the only ID that is ensured to
> be unique on a per-file basis (even in the case of technical replicates), is
> the _Run_ accession. For this reason, it will be used here as the base
> reference for the construction of _xSeries_ and _xModel_ data structures.

## SeqLoader Package
### General
___SeqLoader___ implements S3 `xSeries` and `xModel` classes to provide an
integrated representation of RNA-Seq data (i.e., read counts), gene annotation,
and _Sample_/_Run_ metadata for _Series_ and _Models_, respectively.

Both `new_xSeries` and `new_xModel` ___SeqLoader___ constructors assume that
RNA-Seq low-level analysis (i.e., read alignment and transcript abundance
quantification) has already been performed and that the starting point is a
matrix of raw or normalized counts and their related metadata. This means that
`xSeries` and `xModel` objects can only be constructed _from files_ containing
such information and not by direct data input.

With the exclusion of a few methods designed to filter and subset structures
(`[...]`, `pruneRuns`, `keepRuns`), all the other methods of the class access
the data structures in read-only mode, leaving the object completely unchanged.
___SeqLoader___ is indeed designed to facilitate high-level analysis operations
(i.e., downstream of the count matrix) such as data representation or
descriptive or inferential statistical analysis, in both absolute and
differential expression contexts.  

### `xSeries` objects
These S3 objects are constructed from a gene expression matrix and a metadata
table related to a single study. Within the expression matrix, each row
represents a gene (or transcript), while each column refers to a different
sequencing _Run_, with the possible addition of one or more columns for gene
annotations (e.g., gene symbol, gene name, gene type, ...). In contrast, the
metadata table is assumed to dedicate one row for each _Run_ and one column for
each metadata of interest, among which a mandatory field for _Run_ accessions.
So, `xSeries` objects combine information from a typical gene expression matrix
(possibly including gene annotations) with metadata about individual _Runs_.

More technically, an `xSeries` object is a named list containing one element for
each _Run_ of the study in question, and an additional element for gene
annotations (which are common to all _Runs_). Each _Run_ element is itself a
list containing relevant metadata for the specific _Run_, and a `gene` data
frame with raw or normalized counts for each gene. The `annotation` element, on
the contrary, is a simple data frame (see tree graph below). Notice that,
_Samples_ are not explicitly represented in `xSeries` objects. However, a method
will be implemented to collapse _Runs_ to the samples they belong in the case of
technical replicates (see issue #4).

> [!IMPORTANT]
> For the purpose of a data re-analysis or meta-analysis, one may not
> necessarily be interested in _all_ the samples or experimental conditions
> included in the the original study. Rather, just a selection of them is
> usually considered. For this reason, recomputed count matrices could feature
> less columns (i.e., _Runs_) compared to the number of rows within the complete
> metadata table. By construction, `xSeries` objects are made up of all _Runs_
> for which there exist a metadata entry (i.e., usually all the _Runs_ from the
> original study), but only the ones included in the loaded count matrices will
> also feature a `counts` column in the `genes` data frame. In any case, one can
> always use the `pruneRuns` method to easily remove all the _Runs_ with no
> count data from _xSeries_ or _xModel_ objects, reducing them to the sole
> columns actually present in the count matrix used for object construction.

###  `xModel` objects
Structurally speaking, these S3 objects are just collections (lists) of many
`xSeries` related to the same biological model. Rather, the focus here is on
functions. Methods implementing generics for this class are primarily meant to
aid in meta-analytic data synthesis (see `geneStats` in particular).

Here is a graph of a generic _xModel_ made up of _n_ _xSeries_. The _n-th_
_Series_ consists of _m_ _Runs_ and the gene `annotation` element. Each _Run_
contains a `gene` dataframe for the `counts` and an arbitrary number of
metadata, including the mandatory `ena_run` ID.
```
Model    : xModel
 │
 ├──$ Series_1    : xSeries
 ├──$ Series_2    : xSeries
 ├──$ ...
 └──$ Series_n    : xSeries
    │
    ├──$ Run_1       : list
    ├──$ Run_2       : list
    ├──$ Run_3       : list
    ├──$ ...
    ├──$ Run_m       : list
    │  │
    │  ├──$ ena_sample_title  : chr
    │  ├──$ geo_series        : chr
    │  ├──$ geo_sample        : chr
    │  ├──$ ena_project       : chr
    │  ├──$ ena_sample        : chr
    │  ├──$ ena_run           : chr
    │  ├──$ read_count        : int
    │  ├──$ library_layout    : chr
    │  ├──$ extra             : int
    │  ├──$ ...
    │  └──$ genes             : data.frame
    │     ├──$ IDs               : chr
    │     └──$ counts            : num
    │
    └──$ annotation  : data.frame
       ├──$ IDs         : chr
       ├──$ SYMBOL      : chr
       ├──$ GENENAME    : chr
       ├──$ GENETYPE    : chr 
       └──$ ...
  ```

## Formal Requirements
A number of (hopefully reasonable) assumptions upon file and data organization
are made for object construction to be successful.
1. Each _Series_ (or _Study_) is represented by __two__ `CSV` or `TSV` files,
   containing respectively the read counts (along with possible gene
   annotation) and metadata for all the _Runs_ making up the series.
1. All _Series_ related to the same _Model_ are stored in the same directory,
   here referred to as the _target directory_ for __xModel__ construction.
   Consistently, each _Model_ is expected to live in a separate target
   directory.
1. Both data and metadata files have names starting with the __same__ _Series_
   accession ID (any of the INSDC BioProject, ENA Study, or GEO Series format
   is fine), followed by an underscore (`_`) and either the distinctive string
   `CountMatrix` or `meta` (ignoring case) if they are counts or metadata,
   respectively. Any additional characters after this pattern are allowed, as
   long as the terminal filename extension is either `CSV` or `TSV` (again
   ignoring case). Examples of well-formed filenames are:
   ```
   GSE205739_CountMatrix_genes_TPM.tsv
   PRJNA847413_countMATRIX.CSV
   PRJNA847413_meta.csv
   SRP379268_MeTaDaTa.TSV
   ```
1. Metadata tables have proper column names as header, among which only
   `ena_run` is mandatory, as the only ID that is certainly unique to each FASTQ
   file in the original study (i.e., _Series_).
1. Count tables have proper column names as header, including a mandatory gene
   ID column matching the regex `gene.*id|transcript.*id|ENSEMBL|ENSEMBLTRANS`
   (usually being `gene_id` or `transcript_id`).
1. Within the count table, each count-containing column is named using the
   corresponding `ena_run` ID in the header (along with other possible text
   strings).
1. Input counts from `CountMatrix` file can be either raw or normalized counts,
   but in any case they are supposed to be in linear scale (not log-transformed)

> [!TIP]
> Most of these requirements are automatically met if the _low-level_ RNA-Seq
> data analysis is performed using
> [___x.FASTQ___](https://github.com/TCP-Lab/x.FASTQ).
