DUPLICATION_BACKEND = config["duplication"].get("backend", "dupgen_finder_unique")
DUPLICATION_TARGETS = config["duplication"].get("targets") or [
    record["id"] for record in SPECIES_RECORDS if record.get("outgroup")
]
DUPLICATION_OPTIONAL = [config["duplication"].get("precomputed_table")] if config["duplication"].get("precomputed_table") else []


rule duplication_classification:
    input:
        members=MODULE_TARGETS["family"],
        maps=NORMALIZED_MAPS,
        proteins=NORMALIZED_PROTEINS,
        optional=DUPLICATION_OPTIONAL,
    output:
        modes=MODULE_TARGETS["duplication"],
        pairs=join_path(RESULTS, "08_duplication", "duplication_pairs.tsv"),
        xlsx=join_path(RESULTS, "08_duplication", "duplication.xlsx"),
        plot_pdf=join_path(RESULTS, "08_duplication", "duplication_mode_counts.pdf"),
        plot_png=join_path(RESULTS, "08_duplication", "duplication_mode_counts.png"),
    params:
        backend=DUPLICATION_BACKEND,
        targets=DUPLICATION_TARGETS,
        species_ids=SPECIES,
        species_records=SPECIES_RECORDS,
        separator=SEPARATOR,
        precomputed_table=config["duplication"].get("precomputed_table"),
        dupgen_executable=config["duplication"].get("dupgen_executable", "DupGen_finder-unique.pl"),
        diamond_evalue=config["duplication"].get("diamond_evalue", 1e-10),
        max_target_seqs=config["duplication"].get("max_target_seqs", 5),
        proximal_max_gene_distance=config["duplication"].get("proximal_max_gene_distance", 10),
        extra_args=config["duplication"].get("extra_args", []),
        work_dir=join_path(WORK, "08_duplication"),
        png_dpi=PNG_DPI,
    threads:
        min(32, int(RUN.get("cores", 16)))
    conda:
        "../envs/duplication.yaml"
    script:
        "../scripts/run_duplication.py"
