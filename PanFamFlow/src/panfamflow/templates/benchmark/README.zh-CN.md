# PanFamFlow 真实水稻生物学 benchmark 工作区

此目录用于冻结 5–10 个**相互独立、已经组装并注释**的水稻基因组，以及一个经人工审阅的目标基因家族。这里不接受把同一参考坐标上的多个 BAM/VCF 样本当成多个 genome。

## 目录口径

- `benchmark.yaml`：目标家族、panel 规模和验收阈值；查看最终组间差异前冻结。
- `species.tsv`：每个 assembled genome 的版本、坐标体系、四类必需输入和 SHA256。
- `manual_review/manual_truth_set.tsv`：人工正例、负例、不确定例及证据。
- `decision_log.tsv`：影响结果口径的决策记录。
- `inputs/`：本地输入，仅放副本或只读挂载；不得覆盖原始数据。
- `references/`：版本化 HMM、参考蛋白及来源记录。
- `audits/`：每次审计输出到新的时间戳目录，不覆盖旧结果。

## 初始化后首先执行

```bash
panfamflow benchmark audit \
  --manifest benchmark.yaml \
  --output audits/intake_001 \
  --allow-blocked
```

`--allow-blocked` 只允许生成规划期报告，不会把 BLOCKED 解释为 READY。去掉该选项后，存在阻断项时命令以退出码 2 结束。

## 达到 READY 的必要条件

1. `family.approval_state: approved`，并冻结家族名称、Pfam/InterPro ID、HMM/参考蛋白。
2. `acceptance.approval_state: approved`，在查看最终差异前冻结阈值。
3. 纳入 5–10 个独立 assembled genomes，不能以同一参考上的重测序样本替代。
4. 每个 genome 均有同版本的 genome FASTA、GFF3、protein FASTA、CDS FASTA 和匹配 SHA256。
5. 物种主键唯一、版本与坐标体系明确、恰好一个代表基因组。
6. 人工正负例基线满足预注册最低数量。

软件 CI 通过仅说明代码门禁通过，不等于真实水稻生物学验收通过。
