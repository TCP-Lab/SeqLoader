# SeqLoader Package
###### Constructors and methods for `xSeries` and `xModel` S3 classes

## Motivation
The `xSeries` S3 class attempts to give a practical representation in R of a
typical RNA-Seq output consisting of a series of sequencing runs (usually each
of them referring to a different RNA sample, coming from different biological
sources). Built on this, the `xModel` S3 class addresses the need to collate
several independent scientific studies for the purpose of meta-analysis or
reanalysis. In order to describe such a hierarchical data structure, we will
use a terminology primarily modeled after the one adopted by all the partners of
The International Nucleotide Sequence Database Collaboration (INSDC) in the
context of the Sequence Read Archive (SRA).
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
	be used as _Study_ alias accession, respectively.
1. ___Model___: a collection of _Series_ dealing with the same biological model.
	Since this is the central object of a meta-analysis, it has no corresponding
	object in INSDC.

> [!NOTE]
> Most of the times, _Runs_ and _Samples_ are the same things, however it could
> happen that a single RNA _Sample_ is sequenced through multiple _Runs_ (also
> referred to as 'technical replicates'). Thus, the only ID that is ensured to
> be unique on a per-file basis (even in the case of technical replicates), is
> the _Run_ accession. For this reason, it will be used here as the base
> reference for the construction of _Series_ and _Models_ data structures.

## SeqLoader Package
___SeqLoader___ provides `xSeries` and `xModel` S3 classes to represent _Series_
and _Models_, respectively, as R objects that include both data and metadata.

___SeqLoader___ constructor methods assume that RNA-Seq low-level analysis
(i.e., read alignment and transcript abundance quantification) has already been
performed and that the starting point is a matrix of raw or normalized counts
and the related metadata.

### `xSeries` objects
It is an object of class S3 that is constructed from a gene expression matrix
and a metadata table. In the expression matrix, each row represents a gene (or
transcript), while each column refers to a different sequencing _Run_, with the
possible addition of one or more columns for gene annotations (e.g., gene
symbol, gene name, gene type, ...). In contrast, the metadata table is assumed
to dedicate one row for each _Run_ and one column for each metadata of interest,
among which a mandatory field for _Run_ accessions. The `xSeries` object then
combines information from a typical gene expression matrix (possibly complete
with annotations for genes) with metadata about individual _Runs_. Technically,
an `xSeries` object is an R list containing one element for each _Run_ of the
study in question, and an additional element containing gene annotations (which
are common to each _Run_). Each _Run_ element is itself a list containing
relevant metadata for the specific _Run_, and a `gene` data frame with raw or
normalized counts for each gene. The `annotation` element, on the contrary, is a
simple data frame.

> [!IMPORTANT]
> For the purpose of a data re-analysis or meta-analysis, one may not
> necessarily be interested in _all_ the samples or experimental conditions
> included in the the original study. Rather, just a selection of them is
> usually considered. For this reason, recomputed count matrices could feature
> less columns (i.e., _Runs_) compared to the number of rows of the metadata
> table. The `xSeries` object is made up of all _Runs_ for which there exist an
> entry within the metadata table (i.e., usually all the _Runs_ from the
> original study), but only those _Runs_ included in the loaded count matrix
> will also feature a `counts` column in the `genes` data frame within the
> `xSeries` representation.

###  `xModel` objects
...


## Formal Assumptions
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
	PRJNA847413_countMATRIX.csv
	PRJNA847413_meta.csv
	SRP379268_metadata.tsv
	```
1. Count tables have proper column names as header, including a mandatory gene
	ID column matching the regex `gene.*id|transcript.*id|ENSEMBL|ENSEMBLTRANS`
	(usually being `gene_id` or `transcript_id`).
1. Metadata tables have proper column names as header, among which only
	`ena_run` is mandatory, as the only ID that was certainly unique to each
	FASTQ file in the original study (i.e., _Series_).
1. Within the count table, each count-containing column is named using the
	corresponding `ena_run` ID in the header (along with other possible text
	strings).
1. Input counts from `CountMatrix` file can be either raw or normalized counts,
	but in any case they are supposed to be in linear scale (not log-transformed)

> [!IMPORTANT]
> ENA provides Archive-Generated FASTQs Files for download as an effective means of bringing consistency to data provided. By
   imposing a level of uniformity on these files, users know what to expect
   of them and may incorporate them into pipelines with minimal friction.


