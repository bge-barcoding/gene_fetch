def main():
    print("Starting gene_fetch.py")
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # Validate NCBI credentials first
    print(f"DEBUG: Validating NCBI credentials: email='{args.email}', api_key='{args.api_key}'")
    try:
        Config.validate_credentials(args.email, args.api_key)
        print("DEBUG: Credential validation passed")
    except ValueError as e:
        print(f"ERROR: Credential validation failed: {e}")
        print("DEBUG: About to call sys.exit(1)")
        sys.exit(1)
    print("DEBUG: Continuing after credential validation")

    gene_name = args.gene.lower()
    output_dir = Path(args.out)
    sequence_type = args.type.lower()
    save_genbank = args.genbank  # Get genbank flag

    # Setup output directory and logging
    make_out_dir(output_dir)
    setup_logging(output_dir)

    # Log if GenBank download is enabled
    if save_genbank:
        logger.info(
            "GenBank download mode enabled - will save .gb files in genbank/ subdirectory"
        )

    # Initialize components with required email/api_key
    # No try/catch needed here since we already validated credentials
    config = Config(email=args.email, api_key=args.api_key)

    # Always update thresholds based on user input, regardless of mode
    config.update_thresholds(args.protein_size, args.nucleotide_size)

    # In single-taxid mode, log use of user-specified thresholds
    if args.single:
        logger.info(
            f"Single-taxid mode activated: using protein size threshold {args.protein_size} and nucleotide size threshold {args.nucleotide_size}"
        )

    search_type = config.set_gene_search_term(gene_name)

    if sequence_type not in config.valid_sequence_types:
        print(
            f"Invalid sequence type. Choose from: {', '.join(config.valid_sequence_types)}"
        )
        sys.exit(1)

    logger.info(f"Using {search_type} search terms for {gene_name}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Sequence type: {sequence_type}")

    # Initialize remaining components
    entrez = EntrezHandler(config)
    processor = SequenceProcessor(config, entrez)

    # Check if in single-taxid mode
    if args.single:
        logger.info(f"Single-taxid mode activated for taxid: {args.single}")

        if args.max_sequences:
            logger.info(
                f"Maximum number of sequences to fetch: {args.max_sequences}"
            )
            if sequence_type == "both":
                logger.info(
                    "Note: The max_sequences limit will be applied separately to protein and nucleotide sequences"
                )

        process_single_taxid(
            taxid=args.single,
            gene_name=gene_name,
            sequence_type=sequence_type,
            processor=processor,
            output_dir=output_dir,
            max_sequences=args.max_sequences,
            save_genbank=save_genbank,  # Pass genbank flag
        )
        logger.info("Single taxid processing completed")
        sys.exit(0)
    elif args.max_sequences is not None:
        logger.warning(
            "--max-sequences parameter is ignored when not in single taxid mode"
        )

    # Create output manager
    output_manager = OutputManager(output_dir)

    # Check input file requirements
    if args.input_csv is None and args.input_taxonomy_csv is None:
        logger.error(
            "Error: Either input CSV file (-i/--in) or input taxonomy CSV file (-i2/--in2) must be provided"
        )
        sys.exit(1)

    # Process input samples.csv
    if args.input_csv:
        logger.info(
            f"Starting gene fetch for {gene_name} using taxids from {args.input_csv}"
        )
        process_taxid_csv(
            args.input_csv,
            gene_name,
            sequence_type,
            processor,
            output_manager,
            save_genbank,  # Pass genbank flag
        )

    # Process input samples_taxonomy.csv
    elif args.input_taxonomy_csv:
        logger.info(
            f"Starting gene fetch for {gene_name} using taxonomy from {args.input_taxonomy_csv}"
        )
        process_taxonomy_csv(
            args.input_taxonomy_csv,
            gene_name,
            sequence_type,
            processor,
            output_manager,
            entrez,
            save_genbank,
        )

    logger.info("***********************************************************")
    logger.info("              ? ? ? Gene fetch complete ? ? ?              ")
    logger.info("***********************************************************")
