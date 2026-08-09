# 验证记录

## 声明规则

只有分支中真实存在 `PanFamFlow/`、源码可直接浏览，并且同一源码状态通过自动化检查后，才可标记为软件工程层面的已发布。临时 staging、压缩分片或外部归档不算源码发布。

## v0.1.1-alpha 软件工程验证

验证日期：2026-08-09

已在 GitHub Actions 的 Ubuntu 24.04 / Python 3.12 环境中完成：

```text
uv 0.12.3 lock and locked sync
Ruff lint and formatting
mypy strict checking
27 pytest tests
wheel and source-distribution build
PanFamFlow CLI version/list/validate/plan
Snakemake 9.25.1 toy DAG dry-run
actual toy QC first execution
second identical execution with unchanged hashes and mtimes
deletion of qc.done followed by downstream-only reconstruction
three terminal provenance records with COMPLETED / exit code 0
```

关键兼容性修复包括：

- 在 Snakemake 多值 `--rerun-triggers` 与工作流目标之间加入 `--` 终止符，避免目标被解析为触发器；
- 移除 `workflow/scripts/*.py` 中不兼容 Snakemake `script:` 包装器的 `from __future__ import annotations`；
- 将工作流脚本测试所需的 Biopython 与 Matplotlib 纳入锁定的开发依赖；
- 保留原子输出、运行指纹和最终 `COMPLETED`/`FAILED` provenance 状态。

首次执行、自动跳过和局部恢复已经由 GitHub Actions run `31304925922` 实际验证；该运行还生成了非空的 verified-source artifact。长期 CI 会在每次相关提交和 PR 更新时重复执行同类门禁。

清理后的精确提交 `60edf5504f0f2e7a80508ae5f84662deef6f5e37` 已由长期 CI push run `31305053638` 验证成功。

## v0.1.2 benchmark-gate 开发验证

在不修改 v0.1.1 基线的独立工作目录中，已使用 Python 3.13.5 执行：

```text
python -m compileall -q src tests
PYTHONPATH=src pytest -q
```

结果：

```text
33 passed
```

新增回归覆盖：

- benchmark 初始化拒绝覆盖非空目录；
- 完整双基因组 fixture 达到 `READY`；
- 缺失 genome/GFF3/protein/CDS 时保持 `BLOCKED`；
- 输入内容发生 SHA256 漂移时保持 `BLOCKED`；
- reference-aligned sample 不能替代 assembled genome；
- 中文 HTML、JSON、TSV、XLSX 和 SHA256 receipt 均非空；
- 规划期 `--allow-blocked` 不改变报告中的阻断状态。

Ruff、mypy、wheel/sdist 和 benchmark CLI smoke test 由新分支长期 CI 在精确提交上复核。

## 研究范围验证

标准配置固定为：

```yaml
project:
  analysis_scope: target_pan_gene_family

pan_family:
  # target-family HOG occupancy and classification
```

标准模块为 `pan_family`，不提供 `whole_genome` 分析范围。OrthoFinder 可在完整 canonical proteome 背景中推断 HOG，但最终占有率、Core/Soft-core/Shell/Cloud 分类和下游整合仅投影到已鉴定的目标家族成员。

## 尚未完成的生物学验收

软件工程验证不等于生物学结论已经得到验证。以下仍未完成：

- 5–10 个高质量水稻基因组及一个人工核验目标家族的端到端 benchmark；
- apparent absence 的 TBLASTN/miniprot、共线性区域和 assembly-gap 救援；
- OrthoFinder HOG node 的真实材料人工核验；
- DupGen 外群与复制类型稳定性评估；
- 大规模 Ka/Ks 分块、饱和过滤及统计独立性评估；
- 多材料 RNA-seq reference bias、批次和生物学重复验收；
- 来源模板全部图件的逐图复现。

因此当前版本保持 `v0.1.1-alpha`，PR 保持 Draft，不应被表述为生产级生物学结果流程。
