EXPRESSION_MODE = config["expression"].get("mode", "imported_matrix")
EXPRESSION_MATRIX = config["inputs"].get("expression_matrix") or join_path(
    WORK, "11_expression", "UNCONFIGURED_EXPRESSION_MATRIX.tsv"
)
EXPRESSION_OUTPUTS = {
    "matrix": MODULE_TARGETS["expression"],
    "long": join_path(RESULTS, "11_expression", "expression_long.tsv"),
    "summary": join_path(RESULTS, "11_expression", "expression_summary.tsv"),
    "xlsx": join_path(RESULTS, "11_expression", "expression.xlsx"),
    "plot_pdf": join_path(RESULTS, "11_expression", "expression_heatmap.pdf"),
    "plot_png": join_path(RESULTS, "11_expression", "expression_heatmap.png"),
}

if EXPRESSION_MODE == "imported_matrix":
    rule expression_import:
        input:
            members=MODULE_TARGETS["family"],
            matrix=EXPRESSION_MATRIX,
        output:
            matrix=EXPRESSION_OUTPUTS["matrix"],
            long=EXPRESSION_OUTPUTS["long"],
            summary=EXPRESSION_OUTPUTS["summary"],
            xlsx=EXPRESSION_OUTPUTS["xlsx"],
            plot_pdf=EXPRESSION_OUTPUTS["plot_pdf"],
            plot_png=EXPRESSION_OUTPUTS["plot_png"],
        params:
            separator=SEPARATOR,
            min_tpm_detected=config["expression"].get("min_tpm_detected", 1.0),
            heatmap_transform=config["expression"].get("heatmap_transform", "log2_tpm1_zscore"),
            png_dpi=PNG_DPI,
        conda:
            "../envs/expression.yaml"
        script:
            "../scripts/import_expression.py"
else:
    rule fastp_reads:
        input:
            r1=lambda wildcards: sample_field(wildcards, "r1"),
            r2=lambda wildcards: SAMPLE_BY_ID[wildcards.sample].get("r2") or sample_field(wildcards, "r1"),
        output:
            r1=join_path(WORK, "11_expression", "fastp", "{sample}.R1.clean.fastq.gz"),
            r2=join_path(WORK, "11_expression", "fastp", "{sample}.R2.clean.fastq.gz"),
            json=join_path(RESULTS, "11_expression", "fastp", "{sample}.json"),
            html=join_path(RESULTS, "11_expression", "fastp", "{sample}.html"),
        params:
            paired=lambda wildcards: bool(SAMPLE_BY_ID[wildcards.sample].get("r2")),
            extra_args=config["expression"].get("fastp_extra_args", []),
        threads:
            min(8, int(RUN.get("cores", 16)))
        log:
            stdout=join_path(LOGS, "11_expression", "fastp", "{sample}.stdout.log"),
            stderr=join_path(LOGS, "11_expression", "fastp", "{sample}.stderr.log"),
        conda:
            "../envs/expression.yaml"
        script:
            "../scripts/run_fastp.py"

    rule hisat2_index:
        input:
            genome=lambda wildcards: species_field(wildcards, "genome")
        output:
            done=join_path(WORK, "11_expression", "hisat2_index", "{species}", "index.done")
        params:
            prefix=lambda wildcards: join_path(WORK, "11_expression", "hisat2_index", wildcards.species, "genome")
        threads:
            min(16, int(RUN.get("cores", 16)))
        log:
            join_path(LOGS, "11_expression", "hisat2_index", "{species}.log")
        conda:
            "../envs/expression.yaml"
        shell:
            "mkdir -p $(dirname {params.prefix}) $(dirname {log}); "
            "rm -f {params.prefix}*.ht2 {params.prefix}*.ht2l; "
            "hisat2-build -p {threads} {input.genome} {params.prefix} > {log} 2>&1; touch {output.done}"

    rule align_rnaseq:
        input:
            index=lambda wildcards: join_path(
                WORK,
                "11_expression",
                "hisat2_index",
                SAMPLE_BY_ID[wildcards.sample]["species_id"],
                "index.done",
            ),
            r1=rules.fastp_reads.output.r1,
            r2=rules.fastp_reads.output.r2,
        output:
            bam=join_path(WORK, "11_expression", "alignment", "{sample}.sorted.bam"),
            bai=join_path(WORK, "11_expression", "alignment", "{sample}.sorted.bam.bai"),
        params:
            index_prefix=lambda wildcards: join_path(
                WORK,
                "11_expression",
                "hisat2_index",
                SAMPLE_BY_ID[wildcards.sample]["species_id"],
                "genome",
            ),
            paired=lambda wildcards: bool(SAMPLE_BY_ID[wildcards.sample].get("r2")),
            extra_args=config["expression"].get("hisat2_extra_args", []),
        threads:
            min(16, int(RUN.get("cores", 16)))
        log:
            hisat2=join_path(LOGS, "11_expression", "hisat2", "{sample}.log"),
            samtools=join_path(LOGS, "11_expression", "samtools", "{sample}.log"),
        conda:
            "../envs/expression.yaml"
        script:
            "../scripts/run_hisat2.py"

    rule stringtie_quantify:
        input:
            bam=rules.align_rnaseq.output.bam,
            gff3=lambda wildcards: join_path(
                RESULTS,
                "01_normalized",
                f"{SAMPLE_BY_ID[wildcards.sample]['species_id']}.canonical.gff3",
            ),
        output:
            gtf=join_path(WORK, "11_expression", "stringtie", "{sample}.gtf"),
            abundance=join_path(WORK, "11_expression", "stringtie", "{sample}.abundance.tsv"),
        params:
            strandedness=lambda wildcards: SAMPLE_BY_ID[wildcards.sample].get("strandedness", "unstranded"),
            extra_args=config["expression"].get("stringtie_extra_args", []),
        threads:
            min(8, int(RUN.get("cores", 16)))
        log:
            stdout=join_path(LOGS, "11_expression", "stringtie", "{sample}.stdout.log"),
            stderr=join_path(LOGS, "11_expression", "stringtie", "{sample}.stderr.log"),
        conda:
            "../envs/expression.yaml"
        script:
            "../scripts/run_stringtie.py"

    rule expression_combine_stringtie:
        input:
            members=MODULE_TARGETS["family"],
            maps=NORMALIZED_MAPS,
            abundance=expand(join_path(WORK, "11_expression", "stringtie", "{sample}.abundance.tsv"), sample=SAMPLES),
        output:
            matrix=EXPRESSION_OUTPUTS["matrix"],
            long=EXPRESSION_OUTPUTS["long"],
            summary=EXPRESSION_OUTPUTS["summary"],
            xlsx=EXPRESSION_OUTPUTS["xlsx"],
            plot_pdf=EXPRESSION_OUTPUTS["plot_pdf"],
            plot_png=EXPRESSION_OUTPUTS["plot_png"],
        params:
            sample_ids=SAMPLES,
            sample_species_ids=[SAMPLE_BY_ID[sample]["species_id"] for sample in SAMPLES],
            min_tpm_detected=config["expression"].get("min_tpm_detected", 1.0),
            heatmap_transform=config["expression"].get("heatmap_transform", "log2_tpm1_zscore"),
            png_dpi=PNG_DPI,
        conda:
            "../envs/expression.yaml"
        script:
            "../scripts/combine_stringtie.py"
