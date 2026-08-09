rule concatenate_canonical_proteomes:
    input:
        NORMALIZED_PROTEINS
    output:
        join_path(WORK, "02_family", "all_canonical_proteins.fa")
    conda:
        "../envs/family.yaml"
    script:
        "../scripts/concat_fastas.py"


HMM_ENABLED = bool(config["family"].get("hmm", {}).get("enabled", False))
if HMM_ENABLED:
    HMM_OPTIONS = (
        "--cut_ga"
        if config["family"]["hmm"].get("cut_ga", False)
        else f"-E {config['family']['hmm']['evalue']} --domE {config['family']['hmm']['domain_evalue']}"
    )

    rule family_hmmsearch:
        input:
            proteins=rules.concatenate_canonical_proteomes.output[0],
            hmm=config["family"]["hmm"]["hmm"],
        output:
            domtbl=join_path(WORK, "02_family", "family.domtblout")
        params:
            options=HMM_OPTIONS
        threads:
            min(16, int(RUN.get("cores", 16)))
        log:
            join_path(LOGS, "02_family", "hmmsearch.log")
        conda:
            "../envs/family.yaml"
        shell:
            "mkdir -p $(dirname {output.domtbl}) $(dirname {log}); "
            "hmmsearch --cpu {threads} {params.options} --domtblout {output.domtbl} "
            "{input.hmm} {input.proteins} > {log} 2>&1"
else:
    rule family_hmmsearch:
        input:
            proteins=rules.concatenate_canonical_proteomes.output[0]
        output:
            domtbl=join_path(WORK, "02_family", "family.domtblout")
        log:
            join_path(LOGS, "02_family", "hmmsearch.log")
        shell:
            "mkdir -p $(dirname {output.domtbl}) $(dirname {log}); : > {output.domtbl}; "
            "echo 'HMM search disabled' > {log}"


BLAST_ENABLED = bool(config["family"].get("blast", {}).get("enabled", False))
if BLAST_ENABLED:
    rule family_blastp:
        input:
            proteins=rules.concatenate_canonical_proteomes.output[0],
            queries=config["family"]["blast"]["reference_proteins"],
        output:
            tsv=join_path(WORK, "02_family", "reference_vs_proteome.tsv")
        params:
            db=join_path(WORK, "02_family", "all_canonical_proteins"),
            evalue=config["family"]["blast"]["evalue"],
            max_target_seqs=config["family"]["blast"]["max_target_seqs"],
        threads:
            min(16, int(RUN.get("cores", 16)))
        log:
            join_path(LOGS, "02_family", "blastp.log")
        conda:
            "../envs/family.yaml"
        shell:
            "mkdir -p $(dirname {output.tsv}) $(dirname {log}); "
            "makeblastdb -in {input.proteins} -dbtype prot -out {params.db} > {log} 2>&1; "
            "blastp -query {input.queries} -db {params.db} -evalue {params.evalue} "
            "-max_target_seqs {params.max_target_seqs} -num_threads {threads} "
            "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' "
            "-out {output.tsv} >> {log} 2>&1"
else:
    rule family_blastp:
        input:
            proteins=rules.concatenate_canonical_proteomes.output[0]
        output:
            tsv=join_path(WORK, "02_family", "reference_vs_proteome.tsv")
        log:
            join_path(LOGS, "02_family", "blastp.log")
        shell:
            "mkdir -p $(dirname {output.tsv}) $(dirname {log}); : > {output.tsv}; "
            "echo 'BLASTP search disabled' > {log}"


FAMILY_OPTIONAL_INPUTS = [
    path
    for path in (
        config["family"].get("precomputed_members"),
        config["family"].get("subfamily_assignments"),
        config["family"].get("domain_validation_table"),
        config["family"].get("subcellular_localization_table"),
    )
    if path
]


rule combine_family_evidence:
    input:
        maps=NORMALIZED_MAPS,
        proteins=NORMALIZED_PROTEINS,
        cds=NORMALIZED_CDS,
        hmm=rules.family_hmmsearch.output.domtbl,
        blast=rules.family_blastp.output.tsv,
        optional=FAMILY_OPTIONAL_INPUTS,
    output:
        members=MODULE_TARGETS["family"],
        rejected=join_path(RESULTS, "02_family", "family_candidates_rejected.tsv"),
        xlsx=join_path(RESULTS, "02_family", "family_members.xlsx"),
        proteins=join_path(RESULTS, "02_family", "family_proteins.fa"),
        cds=join_path(RESULTS, "02_family", "family_cds.fa"),
        domains=join_path(RESULTS, "02_family", "family_domains.fa"),
    params:
        family_name=config["family"]["name"],
        separator=SEPARATOR,
        combine_evidence=config["family"]["combine_evidence"],
        calculate_properties=config["family"].get("calculate_protein_properties", True),
        hmm_cut_ga=config["family"]["hmm"].get("cut_ga", False),
        hmm_evalue=config["family"]["hmm"].get("evalue", 1e-5),
        hmm_domain_evalue=config["family"]["hmm"].get("domain_evalue", 1e-3),
        blast_evalue=config["family"]["blast"].get("evalue", 1e-5),
        blast_min_identity=config["family"]["blast"].get("min_identity", 30),
        blast_min_query_coverage=config["family"]["blast"].get("min_query_coverage", 50),
        precomputed_members=config["family"].get("precomputed_members"),
        subfamily_assignments=config["family"].get("subfamily_assignments"),
        domain_validation_table=config["family"].get("domain_validation_table"),
        subcellular_localization_table=config["family"].get("subcellular_localization_table"),
    conda:
        "../envs/family.yaml"
    script:
        "../scripts/combine_family_evidence.py"
