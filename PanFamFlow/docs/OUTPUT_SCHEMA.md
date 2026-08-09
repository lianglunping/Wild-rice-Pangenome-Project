# 关键输出表结构

## `family_members.tsv`

核心列：

- `species_id`, `gene_id`, `transcript_id`, `stable_id`
- `family_name`, `subfamily`
- `chromosome`, `gene_start`, `gene_end`, `strand`
- `evidence_hmm`, `evidence_blast`, `evidence_precomputed`
- HMM/BLAST best-hit fields
- `molecular_weight_da`, `theoretical_pi`, `protein_property_qc`
- optional domain/localization prefixed columns
- `decision`, `rejection_reason`

## `06_pan_family/pan_family_classification.tsv`

- `HOG_ID`, `hog_node`, `hog_node_status`
- `analysis_scope=TARGET_GENE_FAMILY_ONLY`
- `analysis_unit=ORTHOFINDER_HOG`
- `presence_basis=ANNOTATION_AND_HOG_MEMBERSHIP`
- `absence_validation_status=NOT_GENOME_RESCUED`
- `interpretation_flag`
- `species_occupancy`, `species_fraction`, `family_gene_count`
- copy-number metrics
- `pan_family_class`, `is_private`

## `06_pan_family/family_hog_membership.tsv`

- `HOG_ID`, `species_id`, `gene_id`, `stable_id`
- `copy_number_in_species`

## `06_pan_family/unassigned_family_members.tsv`

目标家族成员未出现在所选 HOG node 时保留原 family columns，并增加：

- `reason=NOT_FOUND_IN_SELECTED_HOG_NODE`
- `selected_hog_node`
- `hog_node_status`

## `duplication_mode.tsv`

- `species_id`, `gene_id`, `stable_id`
- `duplication_mode`, `partner_stable_ids`
- `outgroup`, `backend`

## `kaks_pairs.tsv`

- `stable_id_1`, `stable_id_2`
- `pair_type`, `group_id`, `pair_id`
- `Ka`, `Ks`, `Ka_Ks`, `method`
- `qc_status`, `qc_message`, `resumed_from_cache`

## `promoter_elements.tsv`

- `stable_id`, `element`, `major_class`, `subclass`
- FIMO coordinate/score/p/q-value fields when available
- `species_id`, `gene_id`, `promoter_length`, `promoter_qc`

## `expression_matrix.tsv`

前置 metadata 列为 `stable_id`, `species_id`, `gene_id`, `subfamily`，后续为 sample columns。导入矩阵路线另含 `expression_data_status`。

## `master_gene_table.tsv`

以 `stable_id` 为一行，整合 family 证据、gene structure、target-family HOG/pan-family class、chromosome、duplication、per-gene Ka/Ks summary、promoter summary 和 expression summary。

表结构可能在 minor release 增列，但 `stable_id` 语义保持稳定；破坏性改名需要 schema 版本更新。
