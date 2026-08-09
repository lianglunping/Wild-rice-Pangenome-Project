# 模块说明

## `qc`

输入：配置中与所选模块相关的 genome、GFF3、reference protein、HMM、motif、expression/FASTQ 等。

输出：

- `input_audit.tsv/.xlsx`
- `run_manifest.json`
- 可选 `busco_summary.tsv/.xlsx`

失败条件：文件不存在、为空、FASTA 无记录、GFF3 无法解析、BUSCO 运行失败。

## `normalize`

方法：AGAT `agat_sp_keep_longest_isoform.pl` + gffread。

输出：每物种 canonical GFF3、protein、CDS、transcript 和 `gene_transcript_map.tsv/.xlsx`。

关键 QC：stable ID 唯一、每个 gene 只保留一个 canonical transcript、CDS 长度模 3。

## `family`

方法：HMMER `hmmsearch`、BLASTP 或预计算成员表。HMM/BLAST 候选先分别筛选，再按配置合并。

输出：

- `family_members.tsv/.xlsx`
- `family_candidates_rejected.tsv`
- `family_proteins.fa`
- `family_cds.fa`
- `family_domains.fa`

可导入 subfamily、CDD/domain validation 和 subcellular localization 表。大规模项目只加载候选 family 的序列和 mapping 子集，避免把全体 proteome/CDS 同时载入内存。

## `phylogeny`

方法：MAFFT → ClipKIT → IQ-TREE。最少序列数由 `phylogeny.min_sequences` 控制。

输出：alignment、trimmed alignment、`.treefile` 和 IQ-TREE report。

## `gene_structure`

从 canonical GFF3 提取 gene length、exon count、CDS count/length、intron count/length 和 UTR 指标。

## `orthology`

方法：OrthoFinder 3。输入 canonical proteomes；调用 `-X` 保留 PanFamFlow stable ID。

输出：完整 OrthoFinder 工作目录位置和 completion JSON。正式发布结果应保存 `SpeciesTree_rooted_node_labels.txt` 并固定目标 HOG node。

## `pan_family`

读取所选 OrthoFinder `Phylogenetic_Hierarchical_Orthogroups/N*.tsv`，并与 `family_members.tsv` 取交集。只处理目标家族 HOG；不提供 whole-genome 模式，也不构建泛基因组。

输出：

- `pan_family_classification.tsv`
- `family_hog_membership.tsv`
- `family_presence_absence.tsv`
- `unassigned_family_members.tsv`
- rarefaction 原始迭代、汇总和图件

分类表使用 `HOG_ID` 与 `pan_family_class`，同时记录 `presence_basis` 和 `absence_validation_status`。矩阵中的 0 表示注释/HOG 层面未检出，不自动代表已验证的基因丢失。旧预发布模块名 `pangenome` 仅作为兼容别名。

## `chromosome`

读取 family 成员坐标和各物种 genome chromosome 长度。输出每个 gene 的位置、相对位置、coordinate QC，以及每物种/染色体计数与 density。

不同物种保持独立 panel，不直接拼接物理坐标。

## `duplication`

### `dupgen_finder_unique`

对配置 target/outgroup 组合生成 DupGen_finder 输入，使用 DIAMOND 生成相似性文件，再解析 WGD/Tandem/Proximal/Transposed/Dispersed/Singleton。

DupGen_finder 本体通过 `scripts/install_dupgen.sh` 非覆盖安装，并应记录源 commit。

### `precomputed`

导入包含 stable ID、duplication mode 和可选 partner 的表。若 partner 可解析，会同时生成 `duplication_pairs.tsv`，供 Ka/Ks 使用。

## `kaks`

pair source：

- `orthology`：目标 HOG 内 reference species 与其他物种的严格单拷贝 pair。
- `duplication`：duplication partner pair。
- `both`：合并并去重。

方法：protein MAFFT → PAL2NAL codon alignment → KaKs_Calculator。

输出包含失败信息和 QC flags；任何单 pair 失败不会静默删除。

## `promoter`

先根据 strand 提取 promoter：默认 upstream 2000 bp、downstream 0 bp；记录 chromosome 边界截断和邻近 gene overlap。

backend：

- `fimo`：扫描 versioned MEME motif database。
- `precomputed_plantcare`：导入 PlantCARE/其他来源表。

所有 ID 必须可映射到提取的 family promoter。

## `expression`

### `imported_matrix`

导入 wide matrix，筛选 family genes，输出 wide/long/summary 和模式热图。

### `fastq_stringtie`

fastp → HISAT2 → sorted BAM → StringTie `-e -A` → family TPM matrix。gene ID 映射按 sample 的 species scope 处理，允许不同物种存在相同原始 gene ID。

## `report`

整合 family、gene structure、pan_family、chromosome、duplication、Ka/Ks、promoter 和 expression 指标，输出 master table、结果 SHA256 manifest、Python package versions、run info 和静态 HTML gallery。

## Ka/Ks 环境说明

Bioconda 包名为 `kakscalculator2=2.0.1`，安装后提供 `KaKs_Calculator` 命令。流程不会使用不存在的 `kaks_calculator` 包名。
