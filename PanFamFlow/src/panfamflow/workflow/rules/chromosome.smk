rule chromosome_distribution:
    input:
        members=MODULE_TARGETS["family"],
        genomes=[record["genome"] for record in SPECIES_RECORDS],
    output:
        distribution=MODULE_TARGETS["chromosome"],
        summary=join_path(RESULTS, "07_chromosome", "chromosome_summary.tsv"),
        xlsx=join_path(RESULTS, "07_chromosome", "chromosome_distribution.xlsx"),
        plot_pdf=join_path(RESULTS, "07_chromosome", "chromosome_distribution.pdf"),
        plot_png=join_path(RESULTS, "07_chromosome", "chromosome_distribution.png"),
    params:
        species_ids=SPECIES,
        representatives=REPRESENTATIVES,
        representative_only=config["chromosome"].get("representative_only", False),
        png_dpi=PNG_DPI,
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/chromosome_distribution.py"
