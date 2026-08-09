KAKS_PAIR_SOURCE = config.get("kaks", {}).get("pair_source", "both")
KAKS_REFERENCE_SPECIES = config.get("kaks", {}).get("reference_species") or (
    REPRESENTATIVES[0] if REPRESENTATIVES else SPECIES[0]
)
KAKS_MEMBERSHIP = [join_path(RESULTS, "06_pan_family", "family_hog_membership.tsv")] if KAKS_PAIR_SOURCE in {"orthology", "both"} else []
KAKS_DUPLICATION_PAIRS = [join_path(RESULTS, "08_duplication", "duplication_pairs.tsv")] if KAKS_PAIR_SOURCE in {"duplication", "both"} else []


rule calculate_kaks:
    input:
        proteins=join_path(RESULTS, "02_family", "family_proteins.fa"),
        cds=join_path(RESULTS, "02_family", "family_cds.fa"),
        membership=KAKS_MEMBERSHIP,
        duplication_pairs=KAKS_DUPLICATION_PAIRS,
    output:
        tsv=ensure(MODULE_TARGETS["kaks"], non_empty=True),
        xlsx=ensure(join_path(RESULTS, "09_kaks", "kaks_pairs.xlsx"), non_empty=True),
        plot_pdf=join_path(RESULTS, "09_kaks", "kaks_distribution.pdf"),
        plot_png=join_path(RESULTS, "09_kaks", "kaks_distribution.png"),
    params:
        pair_source=KAKS_PAIR_SOURCE,
        reference_species=KAKS_REFERENCE_SPECIES,
        separator=SEPARATOR,
        method=config.get("kaks", {}).get("method", "MA"),
        max_pairs_per_group=config.get("kaks", {}).get("max_pairs_per_group"),
        saturation_ks=float(config.get("kaks", {}).get("saturation_ks", 2.0)),
        workers=int(config.get("kaks", {}).get("workers", 4)),
        work_dir=join_path(WORK, "09_kaks", "pairs"),
        png_dpi=PNG_DPI,
    threads:
        min(int(config.get("kaks", {}).get("workers", 4)), int(RUN.get("cores", 16)))
    resources:
        mem_mb=16000,
        runtime=1440,
    conda:
        "../envs/kaks.yaml"
    retries:
        int(RUN.get("retries", 1))
    script:
        "../scripts/run_kaks.py"
