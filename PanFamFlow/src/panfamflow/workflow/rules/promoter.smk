PROMOTER_BACKEND = config.get("promoter", {}).get("backend", "fimo")
PROMOTER_MOTIF_DATABASE = config.get("promoter", {}).get("motif_database")
PROMOTER_CATEGORY_MAP = config.get("promoter", {}).get("category_map")
PROMOTER_PRECOMPUTED = config.get("promoter", {}).get("precomputed_table")


rule extract_family_promoters:
    input:
        members=MODULE_TARGETS["family"],
        maps=NORMALIZED_MAPS,
        gff3=NORMALIZED_GFFS,
        genomes=[record["genome"] for record in SPECIES_RECORDS],
    output:
        fasta=ensure(join_path(RESULTS, "10_promoter", "family_promoters.fa"), non_empty=True),
        coordinates=ensure(join_path(RESULTS, "10_promoter", "promoter_coordinates.tsv"), non_empty=True),
        coordinates_xlsx=ensure(join_path(RESULTS, "10_promoter", "promoter_coordinates.xlsx"), non_empty=True),
    params:
        species_ids=SPECIES,
        genome_paths=[record["genome"] for record in SPECIES_RECORDS],
        upstream_bp=int(config.get("promoter", {}).get("upstream_bp", 2000)),
        downstream_bp=int(config.get("promoter", {}).get("downstream_bp", 0)),
    conda:
        "../envs/promoter.yaml"
    retries:
        int(RUN.get("retries", 1))
    script:
        "../scripts/extract_promoters.py"


if PROMOTER_BACKEND == "fimo":
    rule scan_promoters_fimo:
        input:
            promoters=rules.extract_family_promoters.output.fasta,
            motifs=PROMOTER_MOTIF_DATABASE,
        output:
            tsv=ensure(join_path(WORK, "10_promoter", "fimo.tsv"), non_empty=True),
        params:
            threshold=float(config.get("promoter", {}).get("fimo_threshold", 1.0e-4)),
            outdir=join_path(WORK, "10_promoter", "fimo_run"),
        log:
            join_path(LOGS, "promoter", "fimo.log")
        conda:
            "../envs/promoter.yaml"
        retries:
            int(RUN.get("retries", 1))
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
            tmpdir="{params.outdir}.partial.$$"
            rm -rf "$tmpdir"
            fimo --thresh {params.threshold} --oc "$tmpdir" {input.motifs} {input.promoters} > {log} 2>&1
            test -s "$tmpdir/fimo.tsv"
            mv "$tmpdir/fimo.tsv" "{output.tsv}.partial.$$"
            rm -rf "$tmpdir"
            mv "{output.tsv}.partial.$$" {output.tsv}
            """
else:
    rule scan_promoters_fimo:
        output:
            tsv=join_path(WORK, "10_promoter", "fimo.tsv"),
        shell:
            "mkdir -p $(dirname {output.tsv}); : > {output.tsv}"


PROMOTER_OPTIONAL_INPUTS = [path for path in (PROMOTER_CATEGORY_MAP, PROMOTER_PRECOMPUTED) if path]


rule parse_promoter_elements:
    input:
        coordinates=rules.extract_family_promoters.output.coordinates,
        fimo=rules.scan_promoters_fimo.output.tsv,
        optional=PROMOTER_OPTIONAL_INPUTS,
    output:
        elements=ensure(MODULE_TARGETS["promoter"], non_empty=True),
        counts=ensure(join_path(RESULTS, "10_promoter", "promoter_element_counts.tsv"), non_empty=True),
        xlsx=ensure(join_path(RESULTS, "10_promoter", "promoter_elements.xlsx"), non_empty=True),
        plot_pdf=join_path(RESULTS, "10_promoter", "promoter_element_counts.pdf"),
        plot_png=join_path(RESULTS, "10_promoter", "promoter_element_counts.png"),
    params:
        backend=PROMOTER_BACKEND,
        category_map=PROMOTER_CATEGORY_MAP,
        precomputed_table=PROMOTER_PRECOMPUTED,
        top_n_elements=int(config.get("promoter", {}).get("top_n_elements", 20)),
        png_dpi=PNG_DPI,
    conda:
        "../envs/promoter.yaml"
    retries:
        int(RUN.get("retries", 1))
    script:
        "../scripts/parse_promoter_elements.py"
