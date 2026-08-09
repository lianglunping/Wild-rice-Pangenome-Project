rule gene_structure_metrics:
    input:
        members=MODULE_TARGETS["family"],
        gff3s=NORMALIZED_GFFS,
        maps=NORMALIZED_MAPS,
    output:
        metrics=MODULE_TARGETS["gene_structure"],
        summary=join_path(RESULTS, "04_gene_structure", "gene_structure_summary.tsv"),
        xlsx=join_path(RESULTS, "04_gene_structure", "gene_structure.xlsx"),
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/extract_gene_structure.py"
