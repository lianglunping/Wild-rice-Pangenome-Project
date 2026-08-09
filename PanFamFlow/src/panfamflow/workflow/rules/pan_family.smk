rule pan_family_classification:
    input:
        result_dir=join_path(RESULTS, "05_orthology", "orthofinder.result_dir.txt"),
        members=join_path(RESULTS, "02_family", "family_members.tsv"),
    output:
        classification=ensure(MODULE_TARGETS["pan_family"], non_empty=True),
        membership=ensure(
            join_path(RESULTS, "06_pan_family", "family_hog_membership.tsv"),
            non_empty=True,
        ),
        presence=ensure(
            join_path(RESULTS, "06_pan_family", "family_presence_absence.tsv"),
            non_empty=True,
        ),
        unassigned_members=ensure(
            join_path(RESULTS, "06_pan_family", "unassigned_family_members.tsv"),
            non_empty=True,
        ),
        rarefaction=ensure(
            join_path(RESULTS, "06_pan_family", "pan_family_rarefaction_iterations.tsv"),
            non_empty=True,
        ),
        rarefaction_summary=ensure(
            join_path(RESULTS, "06_pan_family", "pan_family_rarefaction_summary.tsv"),
            non_empty=True,
        ),
        xlsx=ensure(
            join_path(RESULTS, "06_pan_family", "pan_family_results.xlsx"),
            non_empty=True,
        ),
        class_plot_pdf=join_path(
            RESULTS, "06_pan_family", "pan_family_classification.pdf"
        ),
        class_plot_png=join_path(
            RESULTS, "06_pan_family", "pan_family_classification.png"
        ),
        rarefaction_plot_pdf=join_path(
            RESULTS, "06_pan_family", "pan_family_rarefaction.pdf"
        ),
        rarefaction_plot_png=join_path(
            RESULTS, "06_pan_family", "pan_family_rarefaction.png"
        ),
    params:
        hog_node=config.get("orthofinder", {}).get("hog_node", "auto"),
        species_ids=SPECIES,
        separator=SEPARATOR,
        core_min=float(config.get("pan_family", {}).get("core_min", 0.99)),
        soft_core_min=float(config.get("pan_family", {}).get("soft_core_min", 0.90)),
        shell_min=float(config.get("pan_family", {}).get("shell_min", 0.10)),
        rarefaction_iterations=int(
            config.get("pan_family", {}).get("rarefaction_iterations", 1000)
        ),
        max_exact_combinations=int(
            config.get("pan_family", {}).get("max_exact_combinations", 5000)
        ),
        seed=SEED,
        png_dpi=PNG_DPI,
    conda:
        "../envs/analysis.yaml"
    retries:
        int(RUN.get("retries", 1))
    script:
        "../scripts/parse_pan_family.py"
