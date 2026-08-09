rule normalize_species:
    input:
        genome=lambda wildcards: species_field(wildcards, "genome"),
        gff3=lambda wildcards: species_field(wildcards, "gff3"),
        qc=MODULE_TARGETS["qc"],
    output:
        gff3=join_path(RESULTS, "01_normalized", "{species}.canonical.gff3"),
        proteins=join_path(RESULTS, "01_normalized", "{species}.proteins.fa"),
        cds=join_path(RESULTS, "01_normalized", "{species}.cds.fa"),
        transcripts=join_path(RESULTS, "01_normalized", "{species}.transcripts.fa"),
        mapping=join_path(RESULTS, "01_normalized", "{species}.gene_transcript_map.tsv"),
        mapping_xlsx=join_path(RESULTS, "01_normalized", "{species}.gene_transcript_map.xlsx"),
    params:
        species=lambda wildcards: wildcards.species,
        species_name=lambda wildcards: SPECIES_BY_ID[wildcards.species]["name"],
        group=lambda wildcards: SPECIES_BY_ID[wildcards.species].get("group"),
        subfamily=lambda wildcards: SPECIES_BY_ID[wildcards.species].get("subfamily"),
        separator=SEPARATOR,
        work_dir=lambda wildcards: join_path(WORK, "01_normalized", wildcards.species),
    log:
        agat_stdout=join_path(LOGS, "01_normalized", "{species}.agat.stdout.log"),
        agat_stderr=join_path(LOGS, "01_normalized", "{species}.agat.stderr.log"),
        gffread_stdout=join_path(LOGS, "01_normalized", "{species}.gffread.stdout.log"),
        gffread_stderr=join_path(LOGS, "01_normalized", "{species}.gffread.stderr.log"),
    conda:
        "../envs/normalize.yaml"
    script:
        "../scripts/normalize_canonical.py"


rule normalized_complete:
    input:
        maps=NORMALIZED_MAPS,
        proteins=NORMALIZED_PROTEINS,
        cds=NORMALIZED_CDS,
        gff3=NORMALIZED_GFFS,
    output:
        done=MODULE_TARGETS["normalize"]
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/collect_normalized.py"
