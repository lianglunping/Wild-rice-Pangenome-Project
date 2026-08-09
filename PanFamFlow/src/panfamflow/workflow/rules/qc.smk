AUDIT_RECORDS = []
for record in SPECIES_RECORDS:
    AUDIT_RECORDS.extend(
        [
            {"role": "genome", "species_id": record["id"], "path": record["genome"]},
            {"role": "gff3", "species_id": record["id"], "path": record["gff3"]},
        ]
    )
    if record.get("protein"):
        AUDIT_RECORDS.append(
            {"role": "protein", "species_id": record["id"], "path": record["protein"]}
        )
    if record.get("cds"):
        AUDIT_RECORDS.append(
            {"role": "cds", "species_id": record["id"], "path": record["cds"]}
        )

if config["family"].get("precomputed_members"):
    AUDIT_RECORDS.append(
        {"role": "family_members", "species_id": "", "path": config["family"]["precomputed_members"]}
    )
if config["family"]["hmm"].get("enabled") and config["family"]["hmm"].get("hmm"):
    AUDIT_RECORDS.append(
        {"role": "hmm", "species_id": "", "path": config["family"]["hmm"]["hmm"]}
    )
if config["family"]["blast"].get("enabled") and config["family"]["blast"].get("reference_proteins"):
    AUDIT_RECORDS.append(
        {
            "role": "reference_proteins",
            "species_id": "",
            "path": config["family"]["blast"]["reference_proteins"],
        }
    )
for sample in config["inputs"].get("rnaseq_samples", []):
    AUDIT_RECORDS.append({"role": "fastq_r1", "species_id": sample["species_id"], "path": sample["r1"]})
    if sample.get("r2"):
        AUDIT_RECORDS.append({"role": "fastq_r2", "species_id": sample["species_id"], "path": sample["r2"]})
if config["inputs"].get("expression_matrix"):
    AUDIT_RECORDS.append(
        {"role": "expression_matrix", "species_id": "", "path": config["inputs"]["expression_matrix"]}
    )
if config["inputs"].get("sample_metadata"):
    AUDIT_RECORDS.append(
        {"role": "sample_metadata", "species_id": "", "path": config["inputs"]["sample_metadata"]}
    )

AUDIT_INPUTS = sorted({str(record["path"]) for record in AUDIT_RECORDS})
BUSCO_ENABLED = bool(config.get("qc", {}).get("busco", {}).get("enabled", False))
BUSCO_SUMMARIES = [join_path(RESULTS, "00_qc", "busco", f"{species}.busco.tsv") for species in SPECIES]


rule input_audit:
    input:
        AUDIT_INPUTS
    output:
        tsv=ensure(join_path(RESULTS, "00_qc", "input_audit.tsv"), non_empty=True),
        xlsx=ensure(join_path(RESULTS, "00_qc", "input_audit.xlsx"), non_empty=True),
        manifest=ensure(join_path(RESULTS, "00_qc", "input_manifest.json"), non_empty=True),
    params:
        records=AUDIT_RECORDS,
        calculate_sha256=bool(config.get("qc", {}).get("calculate_sha256", True)),
    log:
        join_path(LOGS, "qc", "input_audit.log")
    conda:
        "../envs/qc.yaml"
    retries:
        int(RUN.get("retries", 1))
    script:
        "../scripts/input_audit.py"


if BUSCO_ENABLED:
    rule busco_species:
        input:
            genome=lambda wildcards: species_field(wildcards, "genome"),
        output:
            summary=ensure(join_path(RESULTS, "00_qc", "busco", "{species}.busco.tsv"), non_empty=True),
        params:
            species=lambda wildcards: wildcards.species,
            lineage=lambda wildcards: species_field(wildcards, "busco_lineage"),
            mode=config["qc"]["busco"].get("mode", "genome"),
            offline=bool(config["qc"]["busco"].get("offline", False)),
            download_path=config["qc"]["busco"].get("download_path"),
            extra_args=config["qc"]["busco"].get("extra_args", []),
            work_dir=lambda wildcards: join_path(WORK, "00_qc", "busco", wildcards.species),
        threads:
            int(config["qc"]["busco"].get("threads", 8))
        log:
            stdout=join_path(LOGS, "qc", "busco", "{species}.stdout.log"),
            stderr=join_path(LOGS, "qc", "busco", "{species}.stderr.log"),
        conda:
            "../envs/busco.yaml"
        retries:
            int(RUN.get("retries", 1))
        script:
            "../scripts/run_busco.py"


rule collect_qc:
    input:
        audit=rules.input_audit.output.tsv,
        busco=BUSCO_SUMMARIES if BUSCO_ENABLED else [],
    output:
        busco_tsv=ensure(join_path(RESULTS, "00_qc", "busco_summary.tsv"), non_empty=True),
        busco_xlsx=ensure(join_path(RESULTS, "00_qc", "busco_summary.xlsx"), non_empty=True),
        done=ensure(MODULE_TARGETS["qc"], non_empty=True),
    params:
        busco_enabled=BUSCO_ENABLED,
    conda:
        "../envs/qc.yaml"
    retries:
        int(RUN.get("retries", 1))
    script:
        "../scripts/collect_qc.py"
