# PanFamFlow decision log

## ADR-001 — Analysis scope

**Decision:** PanFamFlow analyzes one configured target gene family across multiple already assembled and annotated genomes or materials.

**Excluded:** genome assembly, graph-pangenome construction, genome-wide SV/PAV calling, and a general whole-genome pan-gene catalog.

**Reason:** the source analysis template is organized around identification and comparative analysis of one XX gene family, followed by family-level Core/Soft-core/Shell/Cloud, duplication, Ka/Ks, promoter and expression analyses. The expanded workflow document explicitly identified “target family only versus whole-genome analysis” as a variable requiring confirmation; the project owner selected target-family-only analysis.

## ADR-002 — Use of complete proteomes

**Decision:** OrthoFinder may use complete canonical proteomes because orthology inference requires genomic context. Its downstream HOG occupancy and pan-family classification are restricted to the configured target-family members.

**Consequence:** complete-proteome computation does not change the scientific unit of analysis into a whole-genome pangenome.

## ADR-003 — Terminology

The following identifiers are distinct and must not be used interchangeably:

- `subfamily_id`: a target-family phylogenetic/functional subgroup;
- `HOG_ID`: a hierarchical orthogroup at a selected species-tree node;
- `pan_family_class`: Core, Soft-core, Shell or Cloud occupancy class of a target-family HOG;
- `pan_gene_locus_id`: a material-level homologous locus requiring synteny or genome projection; this is not inferred from a tree clade alone.

The canonical module and output directory are `pan_family` and `06_pan_family`.

## ADR-004 — Absence interpretation

**Decision:** absence from an annotation or selected HOG table is not reported as validated gene loss.

**Required evidence for strong loss claims:** genome-level rescue/search, assessable syntenic interval, assembly-gap review and appropriate sequence/coverage evidence. Until then, the workflow records an annotation/HOG occupancy state and an explicit uncertainty flag.

## ADR-005 — Release evidence

Software CI, unit tests, toy DAG execution and smart-resume tests establish software behavior only. A manually reviewed real rice-family benchmark remains required before publication-grade biological claims or a production release.
