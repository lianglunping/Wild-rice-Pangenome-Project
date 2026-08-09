# 可复现性与审计

## 1. 版本层

- Python CLI：`pyproject.toml` + `uv.lock`。
- Snakemake engine：`environment.yaml` 或 `environment-slurm.yaml`。
- 规则级工具：`src/panfamflow/workflow/envs/*.yaml`。
- DupGen_finder：安装后记录 `git rev-parse HEAD`。
- 代码：记录仓库 commit SHA 和 release/tag。

## 2. 数据层

`qc` 模块记录：

- 输入绝对路径
- 文件大小
- SHA256
- FASTA/GFF3 基本统计
- 审计状态和错误说明

不覆盖原始输入。需要替换数据时，建议创建新项目版本目录或更新路径并保留旧 config/manifest。

## 3. 参数层

所有研究阈值写入 `config.yaml`。随机步骤读取 `project.seed`，默认 `20260807`。

正式结果应保存：

```text
config.yaml
uv.lock
environment*.yaml
src/panfamflow/workflow/envs/*.yaml
results/00_qc/run_manifest.json
results/report/result_manifest.tsv
results/report/run_info.json
Git commit SHA
```

## 4. 输出层

- TSV 为机器可读主格式。
- XLSX 用于人工核查，不作为唯一数据源。
- PDF 为矢量图主格式。
- PNG 默认 600 dpi，用于预览或投稿系统。
- HTML report 只做索引，不替代表格。

## 5. 失败与不确定性

不应通过经验补齐缺失结果。以下情况必须明确标记：

- HOG node 仍为 auto。
- family candidate 只通过单一证据。
- promoter 发生 chromosome boundary truncation 或邻近 gene overlap。
- Ka/Ks 失败、Ks=0、非有限 ratio 或高 Ks。
- expression 样本缺失、无重复或跨物种绝对比较。
- absence 未经基因组层复核。

## Python 锁文件

Python 依赖由 `uv 0.12.3` 解析并写入 `uv.lock`。CI 使用 `uv lock --check` 和 `uv sync --locked --dev`，禁止静默更新。
