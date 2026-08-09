REPORT_DEPENDENCIES = [
    MODULE_TARGETS[module] for module in SELECTED_MODULES if module != "report"
]


rule integrated_report:
    input:
        REPORT_DEPENDENCIES
    output:
        index=MODULE_TARGETS["report"],
        master_tsv=join_path(RESULTS, "12_integrated", "master_gene_table.tsv"),
        master_xlsx=join_path(RESULTS, "12_integrated", "master_gene_table.xlsx"),
        manifest=join_path(RESULTS, "report", "result_manifest.tsv"),
        software_versions=join_path(RESULTS, "report", "software_versions.tsv"),
        run_info=join_path(RESULTS, "report", "run_info.json"),
    params:
        results_dir=RESULTS,
        selected_modules=",".join(SELECTED_MODULES),
        title=config.get("report", {}).get("title"),
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/build_report.py"
