# 真实生物学 benchmark 启动门

## 1. 目的

该启动门用于判断一个真实水稻目标泛基因家族项目是否具备启动条件。它不会运行 HMMER、OrthoFinder、DupGen_finder 或 RNA-seq；它只审计输入、版本、预注册和人工基线，避免在口径未冻结时直接进入长时间计算。

启动门严格区分：

- **软件工程验证**：代码、类型、测试、构建、toy DAG 和恢复机制通过；
- **生物学验收**：真实 genome panel、目标家族、HOG 节点、缺失救援、复制外群、Ka/Ks pair 和表达设计经人工核验。

前者不能替代后者。

## 2. 为什么必须使用 assembled genomes

`pan_family` 的占有率单位是独立组装并注释的 genome/material。多个样本的 BAM、VCF 或 gVCF 即使来自不同材料，只要它们都比对到同一个参考坐标，就仍是 reference-aligned samples，不能直接替代多套 genome FASTA + GFF3 + protein + CDS。

否则会把以下现象混在一起：

- 参考偏倚；
- annotation absence；
- assembly gap；
- 真正的 gene absence；
- 重测序覆盖不足。

## 3. CLI

### 3.1 初始化

```bash
panfamflow benchmark init benchmarks/rice_pilot
```

初始化操作拒绝覆盖非空目录，并生成：

```text
benchmark.yaml
species.tsv
decision_log.tsv
README.zh-CN.md
manual_review/manual_truth_set.tsv
.panfamflow/benchmark.schema.json
```

### 3.2 审计

```bash
panfamflow benchmark audit \
  --manifest benchmarks/rice_pilot/benchmark.yaml \
  --output benchmarks/rice_pilot/audits/intake_001
```

存在阻断项时，审计仍会写出完整报告，然后以退出码 `2` 结束。规划阶段可使用：

```bash
panfamflow benchmark audit \
  --manifest benchmarks/rice_pilot/benchmark.yaml \
  --output benchmarks/rice_pilot/audits/intake_001 \
  --allow-blocked
```

`--allow-blocked` 只改变进程退出码，不改变报告中的 `BLOCKED` 状态。

## 4. 状态语义

| 状态 | 含义 | 是否可启动正式计算 |
|---|---|---:|
| `READY` | 无阻断项、无警告项 | 是 |
| `REVIEW` | 无阻断项，但存在需要人工接受的降级项 | 需签字后决定 |
| `BLOCKED` | 至少一个 P0 阻断项失败 | 否 |

## 5. P0 阻断项

### 5.1 目标家族

- `family.name` 不能是 `XX`、`TARGET_FAMILY`、`UNRESOLVED` 等占位符；
- `family.approval_state` 必须冻结为 `approved`；
- 至少配置一个 Pfam 或 InterPro ID；
- 至少有一个非空、可读的 HMM 或人工审阅参考蛋白 FASTA；
- 人工 truth set 的 species ID 必须闭合。

### 5.2 预注册阈值

`acceptance.approval_state` 必须在查看最终组间差异前冻结为 `approved`。阈值变更会改变 pass/fail，应记录在 `decision_log.tsv`。

### 5.3 Genome panel

- 纳入 5–10 个独立 assembled genomes；
- `species_id` 全局唯一，且不得包含稳定 ID 分隔符 `__`；
- 每个 genome 必须记录 assembly accession、annotation version 和 coordinate system；
- 恰好指定一个代表基因组；
- `reference_aligned_sample` 不能计入 genome panel。

### 5.4 四类输入与完整性

每个纳入 genome 必须同时具备：

```text
genome FASTA
GFF3
protein FASTA
CDS FASTA
```

每个文件必须：

1. 是常规文件且非空；
2. 通过轻量格式嗅探；
3. 在 `species.tsv` 登记预期 SHA256；
4. 现场计算的 SHA256 与预期值一致。

审计只读取文件，不修改输入。

### 5.5 人工正负例

`manual_truth_set.tsv` 至少包含预注册数量的：

- `POSITIVE`：确定的家族成员；
- `NEGATIVE`：相似但不属于目标家族的蛋白；
- 可选 `UNCERTAIN`、`NOT_ASSESSABLE`。

必需字段：

```text
species_id
gene_id
expected_status
evidence
reviewer
```

## 6. 警告项

默认情况下，非 chromosome-level assembly 产生 `WARN`，因为染色体定位、共线性和复制解释会降级。可在 `benchmark.yaml` 中将策略设为：

```yaml
acceptance:
  chromosome_level_policy: ignore  # 不提醒
  # 或 warn
  # 或 block
```

群体数少于 2 也产生警告；它不会阻止只做家族鉴定，但会阻止有意义的群体比较。

## 7. 输出

每次审计输出到新的目录，不覆盖历史结果：

| 文件 | 用途 |
|---|---|
| `benchmark_readiness.html` | 中文人工审阅入口 |
| `benchmark_readiness.md` | 跨会话机器可读摘要 |
| `benchmark_readiness.json` | 自动化状态和完整证据 |
| `benchmark_readiness.tsv` | 检查项明细 |
| `benchmark_readiness.xlsx` | summary/checks/species/input_files 工作表 |
| `species_snapshot.tsv` | 本次审计使用的 panel 快照 |
| `input_files.tsv` | expected/observed SHA256、大小和格式状态 |
| `SHA256SUMS.tsv` | 本次报告文件的完整性清单 |

## 8. READY 之后仍需完成的验收

启动门通过后，仍不能直接宣称生物学结论成立。正式 benchmark 还需完成：

1. Family member precision/recall 与 rejected/rescued/fragment 审阅；
2. MSA、模型、UFBoot、SH-aLRT 和亚家族边界审阅；
3. OrthoFinder HOG node 与 unassigned member 审阅；
4. `PRESENT`、`PRESENT_RESCUED`、`ABSENT_VALIDATED`、`UNCERTAIN` 状态机；
5. DupGen 外群与复制类型稳定性；
6. Ka/Ks pair provenance、无效 pair、Ks=0 和饱和审计；
7. motif 的长度/GC 背景与 enrichment；
8. RNA-seq reference bias、batch、生物学重复和 contrast；
9. 真实数据的首次执行、自动跳过、参数失效传播和局部恢复。
