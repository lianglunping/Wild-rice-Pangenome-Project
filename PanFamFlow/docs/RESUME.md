# 断点续跑、自动跳过与故障恢复

PanFamFlow v0.1.1 将 Snakemake 的增量执行作为默认运行语义。正常情况下，用户不需要手工指定“从第几步开始”。修复失败原因后重新执行同一条命令，工作流会重新计算 DAG，并只执行失败、不完整、缺失或过期的任务。

## 默认配置

```yaml
run:
  resume_mode: smart
  keep_going: true
  rerun_incomplete: true
  latency_wait: 120
  retries: 1
  rerun_triggers:
    - mtime
    - input
    - params
    - code
    - software-env
  printshellcmds: true
  show_failed_logs: true
```

`smart` 是推荐模式。`mtime_only` 只按文件时间判断，适合迁移旧结果时临时使用，但不会感知参数、代码和软件环境变化。`off` 关闭 PanFamFlow 自动追加的恢复参数；仍可使用 `panfamflow resume` 强制恢复。

## 标准恢复操作

```bash
uv run panfamflow status -c config.yaml
uv run panfamflow resume -c config.yaml
```

`retry` 是 `resume` 的语义别名：

```bash
uv run panfamflow retry -c config.yaml -m expression
```

不得通过 `touch`、手工制造 `.done` 或删除 `.snakemake/` 来伪造完成状态。`--unlock` 只应在确认没有其他 Snakemake 进程运行且锁确实来自异常退出时使用。

## 自动跳过判定

已完成任务仅在下列条件仍成立时跳过：

1. 声明的关键输出存在；
2. Snakemake 未将其标记为 incomplete；
3. 输入、参数、代码和规则级软件环境未发生需要重跑的变化；
4. 关键表或标记文件通过非空检查；
5. 大型软件自己的完成清单仍与输入签名一致。

仅改变 `cores`、`jobs`、队列、内存、运行时间或重试次数，不改变 biological analysis hash。改变 HMM/BLAST 阈值、canonical transcript 规则、HOG node、泛基因家族占有率阈值、Ka/Ks 配对口径、promoter 长度或表达设计，则改变 analysis hash，并触发相关模块及其下游更新。

## 原子输出

Python 表格、JSON、FASTA、HTML 和命令标准输出先写入同目录下的：

```text
.<filename>.partial.<pid>.<uuid>
```

只有命令成功、文件写完后才通过 `os.replace` 原子发布为正式文件。失败时正式结果不会被半成品覆盖，partial 文件保留用于审计。

## 软件级恢复

### IQ-TREE

`work/03_phylogeny/family.ckp.gz` 保留。输入比对和关键参数签名一致时，重新执行相同命令，IQ-TREE 自动读取 checkpoint。恢复路径不使用 `--redo`。签名变化时，旧 prefix 文件移动到 `work/03_phylogeny/stale/<UTC timestamp>/`。

### OrthoFinder

每套 proteome 内容、物种集合、参数和 OrthoFinder 版本形成一个输入签名：

```text
work/05_orthology/runs/<signature>/
```

完整结果直接复用；存在 `WorkingDirectory` 的未完成结果使用 `orthofinder -b` 尝试继续。输入签名变化时创建新 run，不删除旧 run。

### Ka/Ks

每个基因对使用稳定 pair signature，结果缓存于：

```text
work/09_kaks/pairs/pair_<signature>/result.json
```

完整且签名一致的基因对自动复用；失败基因对单独重算，不会使已成功基因对重新运行。

### BUSCO

BUSCO 输出目录包含 genome、lineage、参数和版本签名。完整 summary 直接复用；不完整目录移动到 `failed/` 后重新运行，不覆盖历史结果。

## Provenance

每次 run/resume/retry 在项目下写入：

```text
.panfamflow/provenance/
├── resolved_config.yaml
├── fingerprints.json
└── runs/<UTC run id>.json
```

其中记录：

- analysis config SHA256；
- execution config SHA256；
- 模块闭包；
- 实际 Snakemake 命令；
- resume mode、triggers、retries；
- 配置文件绝对路径和 UTC 时间。

## 故障注入验收

发布前至少验证：

| 场景 | 验收标准 |
|---|---|
| 中断单物种 HMMER | 其他物种成功任务保留；恢复时只补失败物种及下游 |
| 中断 IQ-TREE | checkpoint 保留并被相同签名重新使用 |
| 命令写出部分 TSV 后失败 | 正式 TSV 不出现，`.partial.*` 被保留 |
| 删除一个合法输出 | 只重建该输出及其下游 |
| 修改 HMM E-value | `qc/normalize` 跳过，`family` 及下游更新 |
| 只修改 cores/memory | biological outputs 不失效 |
| 一个 RNA-seq 样本失败 | 其他样本继续；恢复只补该样本和相关矩阵/contrast |
| 全部完成后再次运行 | Snakemake 报告无需执行或所有目标状态为 complete |
