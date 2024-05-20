# SeqLoader Package
###### Constructors and methods for `bioSeries` and `bioModel` S3 classes

## Motivation

The `bioSeries` S3 class attempts to give an effective representation of the
hierarchical data structure whereby each transcriptomics study typically
consists of a series of sequencing runs (usually each referring to a different
RNA sample, coming from different biological sources).
Built on this, the `bioModel` S3 class addresses the need for grouping several
independent scientific studies together for the purpose of meta-analysis or
reanalysis.
In this regard, we introduce the following dictionary.

## Dictionary

Run:     the set of all raw reads resulting from a single sequencing run,
        referred to by the "Run" ENA accession `(E|D|S)RR[0-9]{6,}`.
        Basically, each individual FASTQ file, or file pair in the case of
        non-interleaved PE reads.
Sample:  the set of all runs from the same biological RNA sample, referred
        to by the "BioSample" INSDC accession `SAM(E|D|N)[A-Z]?[0-9]+`, or
        the "Samples" ENA accession `(E|D|S)RS[0-9]{6,}`, or the
        corresponding GEO sample (alias) accession `GSM[0-9]+`. Most of the
        times, runs and samples are the same things, however it could
        happen that a single RNA sample is sequenced through multiple runs
        (also referred to as 'technical replicates').
Series:  (also referred to as Project (INSDC), or Study (ENA)) the set of
        all samples pertaining to the same experimental design, across all
        conditions of interest to that particular scientific study or
        research project. It can be referred to by the "BioProject" INSDC
        accession `PRJ(E|D|N)[A-Z][0-9]+`, or the "Study" ENA accession
        `(E|D|S)RP[0-9]{6,}`, or the corresponding GEO series (alias)
        accession `GSE[0-9]+`, or the corresponding ArrayExpress experiment
        (alias) accession `E-[A-Z]{4}-[0-9]+`.
Model:   a collection of studies (i.e., a set of Series) from the same
        biological model.

In the following, these terms are italicized when referring to bioSeries
or bioModel objects from SeqLoader package.

NOTE: since the ENA run accession is the only ID that is ensured to be unique
   on a per-file basis (even in the case of technical replicates), it is
   used here as the base reference for the construction of `bioSeries`
   class objects.
   Importantly, ENA provides Archive-Generated FASTQs Files for download as
   an effective  means of bringing consistency to data provided. By
   imposing a level of uniformity on these files, users know what to expect
   of them and may incorporate them into pipelines with minimal friction.






## Assumptions
A number of (hopefully reasonable) assumptions upon file and data organization
are made for object construction to be successful.
1. Each _Series_ (or _Study_) is represented by __two__ `CSV` or `TSV` files,
	containing respectively the read counts (along with possible gene
	annotation) and metadata for all the _Runs_ making up the series.
1. All _Series_ related to the same _Model_ are stored in the same directory,
	here referred to as the _target directory_ for __bioModel__ construction.
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

## `xSeries` objects
E' un oggetto di classe S3 che unisce le informazioni tipiche di una matrice di espressione annotata con i metadati relativi ai singoli campioni (o corse). In altre parole una oggetto xSeries contiene i dati di espressione (o le conte grezze) di tutti i campioni di interesse dello studio), metadati relativi ad ogni campione dlelo studio, e le annotazioni dei geni trasversali ad ogni campione. 
NOTA: l'oggetto 'serie' contine tutte le run per cui esiste un entri nel file di metadati
(tipicamente TUTTE le run dello studio), ma solo quelle con una corrispondente colonna nella matrice delle conte avranno un elemento serie_xxx$run_yyy$genes$counts non nullo!!

A list with  xxx elements, each of them being in turn a list of chars and integer 
and a dataframe containing gene IDs and possibly raw counts 
with attributes yyy

##  `xModel` objects
...


## IMPLEMENTATION NOTE

To make the code lighter, anonymous functions defined in \*apply or pipes are
NOT self-contained, but instead happily access variables from the outer scope.

