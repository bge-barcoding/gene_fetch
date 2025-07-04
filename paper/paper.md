---	
title: 'Gene Fetch: A Python tool for sequence retrieval from GenBank across the tree of life'
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

`Gene Fetch` is an open-source Python tool that automates the retrieval of sequence data from the National Center for Biotechnology Information (NCBI) sequence databases, namely GenBank [@Benson:2012; @Sayers:2024]. GenBank is one of three mirrored partner databases which, along with the European Nucleotide Archive (ENA) [@OCathail2024] and DNA DataBank of Japan (DDBJ) [@Fukuda:2020], form the International Nucleotide Sequence Database Collaboration (INSDC) [@Karsch-Mizrachi:2024]. GenBank currently holds over 3.6 billion nucleotide sequences from over 500,000 formally described species, and exhibits biennial doubling [@Sayers:2023]. 

By streamlining the often tedious and error-prone task of sequence acquisition at scale, `Gene Fetch` addresses a critical need in biological research. It is capable of retrieving both protein and nucleotide sequences for a user-specified genetic marker from across the tree of life, including protein-coding genes (such as cytochrome oxidase subunits, NADH dehydrogenases, RuBisCO, and matK), ribosomal RNA genes (like 16S, 18S, and 28S) and the internal transcribed spacer (ITS) regions. The tool accepts both NCBI taxonomic identifiers (taxids) and hierarchical taxonomic information as input. 

It can systematically traverse taxonomic hierarchies when target sequences are unavailable at the initially specified taxonomic rank (e.g., species → genus → family, etc), documenting the matched taxonomic rank. This is especially valuable for researchers working with non-model organisms or taxonomic groups with sparse genetic data, facilitating retrieval of the taxonomically closest available sequence. The integrated 'batch' mode processes multiple input taxa and retrieve the single ‘best’ sequence per taxa, whilst 'single' mode exhaustively searches for all target sequences for a specified taxon. Collectively, these modes enable efficient retrieval of sequence data for genomic and phylogenetic studies across diverse taxa.




# Statement of need

Comparative genetic analyses within or between taxonomic groups often requires researchers to gather large numbers of gene or protein sequences from public repositories. Streamlined access to these invaluable databases is essential for comprehensive sequence data analysis, however, this task presents several significant challenges:
(1) Limited database curation leading to inconsistent sequence annotation.
(2) Variable sequence representation across taxonomic groups on databases.
(3) Time-intensive manual retrieval processes that do not scale efficiently.

While several existing tools enable sequence retrieval from NCBI sequence repositories, they often require considerable bioinformatics expertise, are limited in functionality or data scope, or are suited for slightly different applications. E-utilities like Entrez Direct [@Kans:2024] offer broad NCBI database access via several APIs but rely on manual NCBI search term construction and significant scripting expertise, burdening the user with navigating database-specific syntax and structure. Similarly, other tools like CRABS: Creating Reference databases for Amplicon-Based Sequencing [@Jeunen:2023] and RESCRIPt: REference Sequence annotation and CuRatIon Pipeline [@Robeson:2021] offer bulk, programmatic retrieval of sequences from several databases (e.g. NCBI, BOLD, etc.), yet they also require manual search terms construction, operate on a single taxon, and often require substantial post-processing by the user. NCBI Datasets [@Oleary2024], while facilitating web- and programmatic-access to NCBI sequence data, is restricted to the curated RefSeq database (a small subset of sequences available on GenBank), limited to species-level queries, and lacks sequence filtering and batch processing capabilities. 

In contrast, `Gene Fetch` offers an accessible, high-throughput solution that automates and simplifies sequence retrieval from GenBank, and requires no prior NCBI syntax knowledge. `Gene Fetch` integrates robust logging, error handling, checkpointing, and a standardised output format, making it suited for reproducible, efficient, and biologically-informed sequence retrieval at scale.

`Gene Fetch` supports 'batch' and 'single' query modes across both protein and nucleotide sequences, with automated CDS extraction, customisable length filtering, and fallback mechanisms for atypical GenBank annotations. It can also process variable GenBank features, including complementary strands, joined sequences, and whole genome shotgun entries, enabling extraction regions of interest from variable feature annotations (e.g., COI from mitogenome records). `Gene Fetch` cross-validates retrieved NCBI taxonomy against input taxonomy, preventing taxonomic homonyms matches (identical names referring to different organisms across the tree of life) and eliminating unreliable database entries (e.g., unverified records, missing data, or substandard annotations). At release, the tool is optimised for 16 common genetic markers, including “barcoding” genes, with curated synonyms for improved search specificity. Users can also specify additional markers, and optionally retrieve corresponding GenBank records for each fetched sequence. 



# Implementation

`Gene Fetch` is implemented in Python `(>=3.9)` and leverages two main libraries: Biopython [@Cock:2009], which, through subpackages of the Bio package (Bio.Entrez, Bio.Seq, Bio.SeqIO, and Bio.SeqRecord), provides the foundation for NCBI database access, sequence parsing and manipulation; and RateLimit, which manages NCBI API rate constraints (beyond those provided in Bio.Entrez) to prevent request throttling. 

The tool follows a modular design with four primary components: configuration manager (handles search parameters and target-specific search term generation), Entrez handler (manages NCBI API interactions with comprehensive error handling), sequence processor (implements core logic for sequence extraction and validation), and output manager (controls file generation and reporting). Detailed logs of parameter and search progress are produced, and checkpoint recovery enables interrupted runs to be resumed. `Gene Fetch` utilises batch processing and taxonomic lineage caching to maximise efficiency, and can process hundreds to thousands of samples (in 'batch' mode) with modest computational resources, as outlined in the GitHub repository.



# Availability

`Gene Fetch` is distributed as a Python package on [PyPi](https://pypi.org/project/gene-fetch/) and a [bioconda package](https://bioconda.github.io/recipes/gene-fetch/README.html), with the source code, testing modules, and a standalone script available under an MIT license through the [GitHub](https://github.com/bge-barcoding/gene_fetch) repository. A supplementary shell wrapper script is also provided in the GitHub repository for submitting a `Gene Fetch` job to a High-Performance Computing (HPC) cluster running the SLURM job scheduler. 



# Acknowledgments and contributions

Contributions following the project's guidelines are welcomed. Users are encouraged to report bugs and suggest enhancements through the GitHub issue tracker. We thank Maria Kamouyiaros and Oliver White (Natural History Museum, London) for their contributions to the conceptualisation and testing of `Gene Fetch`. This work is supported by Biodiversity Genomics Europe (Grant no.101059492), funded by Horizon Europe under the Biodiversity, Circular Economy and Environment call (REA.B.3); co-funded by the Swiss State Secretariat for Education, Research and Innovation (SERI) under contract numbers 22.00173 and 24.00054; and by the UK Research and Innovation (UKRI) under the Department for Business, Energy and Industrial Strategy’s Horizon Europe Guarantee Scheme.



# References

