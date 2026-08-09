## 项目范围硬约束

```yaml
project:
  analysis_scope: target_pan_gene_family
```

该字段只接受 `target_pan_gene_family`。它用于阻止配置在不显式报错的情况下漂移到基因组组装、图泛基因组构建或全基因组 HOG 分类。

# `config.yaml` 说明

## 1. 设计原则

`config.yaml` 是唯一需要用户日常编辑的文件。配置通过 Pydantic 严格校验：字段拼写错误、未知模块、重复物种 ID、错误 outgroup、非法阈值顺序和模块所需输入缺失都会在运行前报告。

初始化项目后，编辑器可读取 `.panfamflow/config.schema.json` 完成补全和类型检查。

## 2. 顶层字段

| 字段 | 含义 |
|---|---|
| `schema_version` | 当前固定为 `1.0` |
| `project` | 项目名、根目录、随机种子和输出路径 |
| `run` | 模块、资源、engine、profile 和 Snakemake 参数 |
| `inputs` | 物种、RNA-seq 样本、表达矩阵与样本元数据 |
| `qc` | SHA256 与 BUSCO |
| `canonical_transcript` | canonical transcript 规则与稳定 ID 分隔符 |
| `family` | HMM/BLAST/预计算证据与外部注释导入 |
| `phylogeny` | MAFFT/ClipKIT/IQ-TREE 参数 |
| `orthofinder` | HOG node 和线程 |
| `pan_family` | 目标家族 HOG 的四分类阈值与 rarefaction |
| `chromosome` | representative 筛选与 density window |
| `duplication` | DupGen_finder 或预计算 backend |
| `kaks` | pair source、reference species、方法和饱和阈值 |
| `promoter` | promoter 长度、FIMO/PlantCARE backend |
| `expression` | imported matrix 或 FASTQ/StringTie 路线 |
| `plot` | PDF/PNG 与 DPI |
| `report` | 报告标题和已有结果整合 |

## 3. 路径解析

除绝对路径外，所有路径均相对于 `project.root`。`project.root` 本身相对于配置文件所在目录解析。

```yaml
project:
  root: .
```

意味着 `data/Os/genome.fa` 对应 `config.yaml` 所在项目目录下的 `data/Os/genome.fa`。

## 4. 模块选择

```yaml
run:
  modules: [family, phylogeny, gene_structure]
```

CLI 会自动解析为 `qc → normalize → family → phylogeny/gene_structure`。

临时覆盖：

```bash
uv run panfamflow run -c config.yaml -m promoter
```

## 5. 物种记录

```yaml
inputs:
  species:
    - id: Os
      name: Oryza_sativa
      genome: data/Os/genome.fa
      gff3: data/Os/annotation.gff3
      protein: null
      cds: null
      group: Cultivated
      subfamily: Oryza
      representative: true
      outgroup: Og
      busco_lineage: poales_odb12
```

约束：

- `id` 必须全局唯一，不能包含保留分隔符 `__`。
- `outgroup` 必须引用已配置物种，且不能指向自身。
- DupGen_finder target 必须有 outgroup。
- BUSCO 启用时每个物种必须配置 lineage。

## 6. 家族证据

```yaml
family:
  name: GPAT
  combine_evidence: intersection
  hmm:
    enabled: true
    hmm: references/PF01553.hmm
    evalue: 1.0e-5
    domain_evalue: 1.0e-3
  blast:
    enabled: true
    reference_proteins: references/known_GPAT.pep.fa
    evalue: 1.0e-5
    min_identity: 30
    min_query_coverage: 50
```

- `union`：HMM 或 BLAST 任一通过。
- `intersection`：两类证据均通过。
- `hmm_only`：只使用 HMM。
- `blast_only`：只使用 BLAST。
- `precomputed_members` 非空时直接以外部成员表为准，仍校验 stable ID。

外部表应优先包含 `stable_id`；也可包含 `species_id + gene_id`。

## 7. 泛基因家族占有率阈值

```yaml
pan_family:
  core_min: 0.99
  soft_core_min: 0.90
  shell_min: 0.10
```

只对包含 `family_members.tsv` 目标家族成员的 HOG 分类：

- `fraction >= core_min`：Core
- `soft_core_min <= fraction < core_min`：Soft-core
- `shell_min <= fraction < soft_core_min`：Shell
- `< shell_min`：Cloud

阈值必须满足 `core_min >= soft_core_min >= shell_min`。PanFamFlow 不接受 `whole_genome` scope；旧 `pangenome` 字段只在 `scope: target_family` 时自动迁移。

## 8. expression 两种输入

### 导入矩阵

```yaml
inputs:
  expression_matrix: references/expression.tsv
expression:
  mode: imported_matrix
```

行 ID 优先使用 `stable_id`。若只使用 gene ID，必须在全体 family genes 中唯一，否则拒绝自动映射。

### FASTQ/StringTie

```yaml
inputs:
  rnaseq_samples:
    - id: Os_control_R1
      species_id: Os
      condition: Control
      tissue: Leaf
      replicate: 1
      strandedness: unstranded
      r1: data/rnaseq/Os_control_R1.fastq.gz
      r2: null
expression:
  mode: fastq_stringtie
```

v0.1 生成 TPM，不自动执行统计推断型差异表达。

## 9. 配置校验

```bash
uv run panfamflow validate -c config.yaml
uv run panfamflow validate -c config.yaml -m promoter
uv run panfamflow schema -o config.schema.json
```

错误导致退出码 2；warning 不阻止执行，但必须在正式分析记录中处理。

## 断点续跑相关字段（v0.1.1）

```yaml
run:
  resume_mode: smart
  keep_going: true
  rerun_incomplete: true
  latency_wait: 120
  retries: 1
  rerun_triggers: [mtime, input, params, code, software-env]
  printshellcmds: true
  show_failed_logs: true
```

- `resume_mode`：`smart`、`mtime_only` 或 `off`。推荐 `smart`。
- `keep_going`：某个 job 失败后，继续执行与它无依赖关系的 job；最终流程仍返回失败状态。
- `rerun_incomplete`：重跑被 Snakemake 标记为 incomplete 的 job。
- `latency_wait`：网络文件系统上等待输出可见的秒数。
- `retries`：单个 job 的额外自动重试次数。格式错误和配置错误不会因重试而自行修复，因此不宜设置过大。
- `rerun_triggers`：决定输出何时因输入、参数、代码、环境或 mtime 变化而失效。
- `printshellcmds`：记录实际命令。
- `show_failed_logs`：失败时显示规则日志。

查看状态：

```bash
uv run panfamflow status -c config.yaml
```

修复问题后恢复：

```bash
uv run panfamflow resume -c config.yaml
```

完整语义见 [RESUME.md](RESUME.md)。
