# Wild rice Pangenome Project

> **Feature-branch addition:** [`PanFamFlow/`](./PanFamFlow/) is a target **pan-gene-family analysis** workflow. It consumes already assembled and annotated genomes; it is not a genome-assembly or graph-pangenome construction pipeline.

## PanFamFlow: target pan-gene-family analysis

PanFamFlow identifies one configured gene family across multiple assembled genomes or materials, then integrates family-member evidence, family phylogeny, target-family HOG occupancy, Core/Soft-core/Shell/Cloud classification, gene structure, chromosome distribution, duplication, Ka/Ks, promoter motifs, expression and a traceable report.

It does **not** assemble genomes, construct graph pangenomes, call genome-wide SV/PAV, or classify every genome-wide HOG as a general pangenome product. Read the [workflow README](./PanFamFlow/README.md), [scope contract](./PanFamFlow/docs/SCOPE.md), [audit notes](./PanFamFlow/docs/AUDIT.md), [resume semantics](./PanFamFlow/docs/RESUME.md), and [validation status](./PanFamFlow/docs/VALIDATION.md).

## Original repository content

The original Wild Rice Pangenome Project code remains available in directories 1–8:

- [1. Genome assembly](./1.%20Genome_assembly/README.md)
- [2. Gene annotation](./2.%20Gene_annotation/README.md)
- [3. TE annotation](./3.%20TE_annotation/README.md)
- [4. Variation calling](./4.%20Variation_calling/README.md)
- [5. Pangenome analysis](./5.%20Pangenome_analysis/README.md)
- [6. Evolutionary analysis](./6.%20Evolutionary_analysis/README.md)
- [7. Domestication analysis](./7.%20Domestication_analysis/README.md)
- [8. Indica–japonica differentiation](./8.%20Indica-japonica_differentiated/README.md)

Original database: [RicePandb](http://ricepandb.ncgr.ac.cn). Original publication: Guo, D., Li, Y., Lu, H. et al. *A pangenome reference of wild and cultivated rice*. Nature (2025). https://doi.org/10.1038/s41586-025-08883-6.
