rule orthofinder:
    input:
        proteomes=NORMALIZED_PROTEINS
    output:
        done=MODULE_TARGETS["orthology"],
        result_dir=join_path(RESULTS, "05_orthology", "orthofinder_result_dir.txt"),
    params:
        species_ids=SPECIES,
        search_threads=config["orthofinder"].get("search_threads", 32),
        analysis_threads=config["orthofinder"].get("analysis_threads", 8),
        extra_args=config["orthofinder"].get("extra_args", []),
        work_dir=join_path(WORK, "05_orthology"),
    threads:
        max(
            int(config["orthofinder"].get("search_threads", 32)),
            int(config["orthofinder"].get("analysis_threads", 8)),
        )
    log:
        stdout=join_path(LOGS, "05_orthology", "orthofinder.stdout.log"),
        stderr=join_path(LOGS, "05_orthology", "orthofinder.stderr.log"),
    conda:
        "../envs/orthology.yaml"
    script:
        "../scripts/run_orthofinder.py"
