# 设计依据与需求追踪

## 1. 来源

本版本根据以下两份项目来源文档实现：

1. `03-泛基因家族分析-模板.pdf`
2. `comparative_genomics_gene_family_workflow_20260807.md`

原始文档不随代码仓库公开；本文件只记录需求摘要、实现决策和解释边界。

## 2. 来源要求到实现的映射

| 来源要求 | PanFamFlow 实现 | 状态/边界 |
|---|---|---|
| BLASTP + HMMER/Pfam 识别家族成员 | `family` 模块；支持 union/intersection/hmm_only/blast_only | 已实现；所有未通过候选保留在拒绝表 |
| CDD 保守结构域验证 | `family.domain_validation_table` 导入 | v0.1 不自动提交 CDD |
| 理化性质 MW/pI | Biopython `ProteinAnalysis` | 已实现；非标准残基有 QC 标记 |
| WoLF PSORT 亚细胞定位 | `family.subcellular_localization_table` 导入 | v0.1 不自动提交网页服务 |
| 多序列比对和系统发育树 | MAFFT + ClipKIT + IQ-TREE | 已实现；presence/absence 聚类不称为系统发育树 |
| gene/exon/intron/UTR 结构 | GFF3 解析 | 已实现；依赖规范的 Parent/ID 关系 |
| Core/Soft-core/Shell/Cloud | 目标家族 HOG occupancy classification | 已实现；不对全基因组 HOG 分类 |
| 泛基因家族 rarefaction | target-family HOG exact/random subset sampling | 已实现；随机种子固定 |
| chromosome distribution | 每物种独立坐标与 density | 已实现；禁止跨物种坐标直接串联 |
| MCScanX / DupGen_finder 复制分类 | DupGen_finder-unique 主路线；预计算导入 | 已实现 alpha；MCScanX 独立 backend 待补 |
| ParaAT/KaKs_Calculator | MAFFT + PAL2NAL + KaKs_Calculator | 已实现 pairwise 路线；高 Ks 标记潜在饱和 |
| 2 kb promoter + PlantCARE | strand-aware promoter；FIMO；PlantCARE 表格导入 | 自动 FIMO 已实现；PlantCARE 网页自动化未实现 |
| fastp/HISAT2/StringTie | FASTQ 路线或矩阵导入 | TPM 已实现；DESeq2 contrasts 待补 |
| 结果图表与主表 | TSV/XLSX、PDF/PNG、master table、HTML report | 数据层和代表性图已实现；模板全部组合图待补 |

## 3. 对综合方案文档的工程化落实

### 3.1 Canonical transcript

每个 gene 默认保留最长 CDS transcript。流程使用 AGAT 选择 isoform，再用 gffread 从同一 genome/GFF3 组合生成 transcript、CDS 和 protein。下游全部使用 `SpeciesID__GeneID` 稳定 ID，避免原始注释中同名基因跨物种冲突。

### 3.2 OrthoFinder HOG 与目标家族限定

OrthoFinder 在 canonical proteomes 上获得 HOG 背景；`pan_family` 读取 `Phylogenetic_Hierarchical_Orthogroups/N*.tsv` 后与 `family_members.tsv` 取交集。流程不把旧式 `Orthogroups.tsv` 当作最终 HOG 结果，也不对全部 HOG 做 whole-genome 分类。OrthoFinder 调用加入 `-X`，因为 canonical protein 已经有全局唯一的物种前缀。

### 3.3 模块化 QC

- 配置 schema 严格拒绝未知字段。
- 输入文件检查存在性、非空、FASTA/GFF3 统计和 SHA256。
- 家族候选同时输出 PASS 与 REJECT 审计表。
- HOG node 自动选择发出 warning。
- promoter 边界和邻近基因重叠输出 QC。
- Ka/Ks 输出非有限值、Ks=0、潜在饱和和失败信息。

### 3.4 跨物种解释边界

- 染色体图按物种分面，不制造跨物种连续坐标轴。
- 表达模块保留原始 TPM，但报告不自动给出跨物种绝对高低结论。
- pairwise Ka/Ks 不自动解释为正选择证据。
- gene absence 不自动解释为真实丢失。

### 3.5 可复现输出

- `project.seed` 控制 random rarefaction。
- 输入审计写入绝对路径、大小和 SHA256。
- 结构化输出写 TSV + XLSX。
- 图件写 PDF + 高分辨率 PNG。
- `result_manifest.tsv` 对结果文件再次计算 SHA256。
- `run_info.json` 记录配置路径和模块闭包。

## 4. 验收口径

工程验收与生物学验收分开：

- 工程验收：schema、CLI、模块闭包、单元测试、静态检查、构建、Snakemake dry-run。
- 生物学验收：在真实数据上核对 canonical 数量、HMM/BLAST 候选、HOG node、关键 absence、复制分类、Ka/Ks QC、promoter 坐标和表达样本设计。

v0.1.1-alpha 只声明完成前者的一部分，不声明真实项目的生物学结论已经验证。
