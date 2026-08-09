# PanFamFlow scope contract

## Primary objective

PanFamFlow answers questions about **one configured target gene family** across multiple assembled and annotated genomes or accessions:

- Which genes are credible family members?
- How are family members assigned to subfamilies and HOGs?
- Which target-family HOGs are Core, Soft-core, Shell or Cloud in the configured sample set?
- How do gene structure, duplication mode, selection pressure, promoter motifs and expression vary within that family?

## Required starting point

PanFamFlow starts from:

```text
assembled genome FASTA
+
matching GFF3 annotation
+
optional source protein/CDS FASTA
+
optional RNA-seq or expression matrix
```

It does not generate the genome assemblies or annotations.

## Explicit non-goals

The workflow does not:

1. assemble genomes from raw reads;
2. build a graph pangenome;
3. construct a whole-genome pan-reference;
4. call whole-genome structural variants or gene PAVs;
5. classify all whole-genome HOGs as a general pangenome product;
6. prove gene loss from annotation absence alone.

Tools such as minigraph, PGGB, Cactus, vg, PanGenome Research Toolkit or whole-genome SV callers are outside this repository's primary scope.

## Why OrthoFinder still uses full proteomes

Orthology inference benefits from complete canonical proteomes. OrthoFinder therefore runs on configured canonical proteomes, but the `pan_family` parser intersects the selected HOG table with `family_members.tsv`. HOGs without a target-family gene are discarded and never enter pan-family classification.

The output records:

```text
analysis_scope = TARGET_GENE_FAMILY_ONLY
analysis_unit = ORTHOFINDER_HOG
presence_basis = ANNOTATION_AND_HOG_MEMBERSHIP
absence_validation_status = NOT_GENOME_RESCUED
```

This prevents whole-genome HOG classification from being confused with target pan-gene-family analysis.

## Terminology

- **Family subfamily/clade**: a group defined from target-family phylogeny, reference genes and domain architecture.
- **HOG**: an OrthoFinder hierarchical orthogroup at a specified species-tree node.
- **Pan-family class**: occupancy class of a target-family HOG in the configured genomes/accessions.
- **Pan-locus**: a syntenically homologous locus across accessions. PanFamFlow v0.1.1 does not yet infer pan-loci automatically.
- **Validated absence**: absence supported by genome-level rescue checks and assessable sequence context. A zero in the current HOG matrix is not automatically validated absence.

## Acceptance guard

The configuration model rejects the former pre-release setting:

```yaml
pangenome:
  scope: whole_genome
```

Legacy `pangenome` naming is migrated only when `scope: target_family`; all new documentation and output use `pan_family`.

## Configuration guard

The canonical configuration fixes:

```yaml
project:
  analysis_scope: target_pan_gene_family
```

Any other value is rejected by the Pydantic schema. This field is included in both analysis and execution provenance fingerprints, so a scope change cannot be hidden as a resource-only change.
