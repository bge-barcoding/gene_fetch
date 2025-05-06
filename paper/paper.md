---
title: 'Gene Fetch: A Python tool for high-throughput sequence retrieval from NCBI databases'
tags:
  - Python
  - Bioinformatics
  - Genomics
  - Biodiversity Genomics
  - Barcoding
authors:
  - name: Daniel A.J. Parsons
    orcid: 0000-0002-5246-0725
    affiliation: 1 
  - name: Benjamin W. Price
    orcid: 0000-0001-5497-4087
    affiliation: 1
affiliations:
 - name: Natural History Museum, Cromwell Road, London, SW7 5BD, United Kingdom
   index: 1
   ror: 000017519k585
date: 07 May 2025
bibliography: paper.bib
---



# Summary

`Gene Fetch` (https://github.com/bge-barcoding/gene_fetch) is an open-source Python tool that automates the retrieval of sequence data from the National Center for Biotechnology Information (NCBI) sequence databases, namely GenBank [@Benson:2012]. GenBank is one of three mirrored partner databases, along with the European Nucleotide Archive (ENA) [@OCathail:2024] and DNA DataBank of Japan (DDBJ) [@Fukuda:2020], that form the International Nucleotide Sequence Database Collaboration (INSDC) [@Karsch-Mizrachi:2024]. As of writing, GenBank comprises >3.6 billion nucleotide sequences from >500,000 formally described species, and exhibits biennial doubling [@Sayers:2023]. 

By streamlining the often tedious and error-prone task of sequence acquisition at scale, `Gene Fetch` addresses a critical need in biological research. It is capable of retrieving both protein and nucleotide sequences for a user-specified genetic marker from across the tree of life, including protein-coding genes (such as cytochrome oxidase subunits, NADH dehydrogenases, RuBisCO and matK), ribosomal RNA genes (like 16S, 18S, and 28S) and the internal transcribed spacer (ITS) regions. The tool accepts both NCBI taxonomic identifiers (taxids) and hierarchical taxonomic information as input. 

It can systematically traverse taxonomic hierarchies when target sequences are unavailable at the initially specified taxonomic rank (e.g., species → genus → family, etc), ensuring retrieval of the taxonomically closest available reference sequence. This functionality is particularly valuable for researchers working with non-model organisms or taxonomic groups with sparse genetic data, and serves to “fetch” the closest available sequence for each of the focal taxa. The integrated ‘batch’ mode processes multiple input taxa and retrieves the ‘best’ sequence per taxa, while the ‘single’ mode enables exhaustive searching and retrieval of sequence targets for a specified taxon. Collectively, these two modes facilitate more efficient retrieval of sequence data for genomic and phylogenetic studies across diverse taxa.




# Statement of need

Comparative genetic analyses within or between taxonomic groups often requires researchers to gather large numbers of gene or protein sequences from public repositories. Streamlined access to these invaluable databases is essential for comprehensive sequence data analysis, however, this task presents several significant challenges: 
(1) Limited database curation leading to inconsistent sequence annotation.
(2) Variable sequence representation across taxonomic groups on databases.
(3) Time-intensive manual retrieval processes that do not scale efficiently.

While several existing tools enable sequence retrieval from NCBI sequence repositories, they often require considerable bioinformatics expertise or are limited in scope. E-utilities like Entrez Direct [@Kans:2024] provide broad access to all NCBI databases via a collection of APIs but rely on manually constructed search terms and significant scripting skills, placing a high burden on the user to navigate database-specific syntax and structure. Similarly, other tools like CRABS: Creating Reference databases for Amplicon-Based Sequencing [@Jeunen:2023] and RESCRIPt: REference Sequence annotation and CuRatIon Pipeline [@Robeson:2021] do offer bulk, programmatic retrieval of sequences from several databases (e.g. NCBI, BOLD, etc.), yet they also require manual construction of complex NCBI search terms, operate on a single taxon at a time, and often necessitate substantial post-processing and filtering by the user.

In contrast, `Gene Fetch` offers an accessible, high-throughput solution that automates and simplifies sequence retrieval from NCBI databases, and unlike existing tools, requires no prior knowledge of NCBI query syntax. `Gene Fetch` integrates robust progress logging, error handling, checkpointing, and provides a structured output (including GenBank files), making it suited for researchers seeking reproducible, efficient, and biologically-informed sequence retrieval at scale.

`Gene Fetch` supports a "batch" and "single" query modes across both protein and nucleotide sequences, with automated CDS extraction, customisable length filtering, and fallback mechanisms for atypical GenBank annotations. Additionally, `Gene Fetch` can process variable GenBank features, such as complementary strands, joined sequences, and whole genome shotgun sequence entries. The tool isolates and extracts regions of interest from annotated features (e.g., COI from a mitogenome record). The tool also cross-validates retrieved NCBI taxonomy against input taxonomy, effectively preventing common issues such as taxonomic homonyms (identical names referring to different organisms across distinct branches of the tree of life) and eliminating unreliable database entries (including unverified records, sequences with missing data, or entries with substandard annotations).

`Gene Fetch` offers several key capabilities:
- Traverses taxonomic hierarchies from lower taxonomic ranks (e.g. species) to higher ranks (e.g. phylum) when sequences are unavailable at the requested level, documenting the matched taxonomic rank.
- Optimised for 15 common genetic markers at point of release, with curated synonyms for improved search specificity, including popular “barcoding” genes. Additional markers and synonyms can be specified by the user.
- Flexible inputs, accepting either taxids or taxonomic hierarchies. 
- Supports both "batch" mode (retrieving one "best" sequence per input taxon) and "single" mode (retrieving all available sequences for a single specified taxon).
- Extracts coding sequences (CDS) and rRNA features from larger records (e.g. mitogenomes and WGS records).
- Implements length-based sequence filtering before retrieval to fit the user’s requirements.
- Includes fallback mechanisms for handling inconsistently annotated entries (sequence feature extraction).
- Produces detailed logs of applied parameters and search progress.
- Checkpoint recovery to resume interrupted runs.
- Utilises batch processing and caching of retrieved taxonomic lineages to maximise efficiency.
- Optional retrieval of GenBank records corresponding to the fetched sequence data.



# Implementation

`Gene Fetch` is implemented in Python `(>=3.9)` and leverages two main libraries: Biopython [@Cock:2009], which, through subpackages of the Bio package (Bio.Entrez, Bio.Seq, Bio.SeqIO, and Bio.SeqRecord), provides the foundation for NCBI database access, sequence parsing and manipulation; and RateLimit, which manages NCBI API rate constraints (beyond those provided in Bio.Entrez) to prevent request throttling. 

The tool follows a modular, community standard design with four primary components: configuration manager (handles search parameters and target-specific search term generation), Entrez handler (manages API interactions with NCBI's Entrez system with comprehensive error handling), sequence processor (implements the core logic for extracting and validating target sequences), and output manager (controls file generation and reporting). `Gene Fetch` can process hundereds-thousands of samples (in 'batch' mode) with modest computational resources, as outlined in the benchmarking section of the GitHub repository.



# Availability

`Gene Fetch` is distributed as a Python package on [PyPi](https://pypi.org/project/gene-fetch/), with the source code, testing modules, and standalone script available under an MIT license through the [GitHub](https://github.com/bge-barcoding/gene_fetch) repository. The repository also holds comprehensive documentation, usage examples, installation instructions, and benchmarking statistics. A supplementary shell wrapper script is also provided for submitting a `Gene Fetch` job to a High-Performance Computing (HPC) cluster running SLURM. 



# Acknowledgments and contributions

Contributions following the project's guidelines are welcomed, and we encourage users to report issues through the GitHub issue tracker. We extend our thanks to Maria Kamouyiaros and Oliver White at the Natural History Museum in London for their contributions to the conceptualisation and testing of `Gene Fetch`. This work is supported by Biodiversity Genomics Europe. Biodiversity Genomics Europe (Grant no.101059492) is funded by Horizon Europe under the Biodiversity, Circular Economy and Environment call (REA.B.3); co-funded by the Swiss State Secretariat for Education, Research and Innovation (SERI) under contract numbers 22.00173 and 24.00054; and by the UK Research and Innovation (UKRI) under the Department for Business, Energy and Industrial Strategy’s Horizon Europe Guarantee Scheme.



# References

