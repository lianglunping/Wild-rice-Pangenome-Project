# 验证记录

## 声明规则

只有分支中真实存在 `PanFamFlow/`、源码可直接浏览，并且同一 head SHA 的长期 CI 通过后，才可标记为已发布。临时 staging 或外部压缩包不算发布。

## v0.1.1-alpha 本地验证

验证日期：2026-08-09

已执行：

```bash
python -m compileall -q src tests
python -m pytest -q
```

当前结果：

```text
待本次修复后由实际测试日志更新
```

覆盖内容：

- 配置加载与重复物种 ID 拒绝
- 模块依赖展开及 Ka/Ks 动态依赖
- module-aware input validation
- CLI version/list/validate/init 非覆盖行为
- list-based Snakemake command 构造
- streaming FASTA parser/length audit
- input audit 脚本
- target pan-family HOG parser、unassigned-member audit 和 whole-genome scope rejection
- StringTie 跨物种同名 gene ID 映射
- OrthoFinder stable ID `-X` 保护

## 尚待远程 CI 确认

当前本地运行环境无法解析外部包索引，因此以下项目由 GitHub Actions 完成：

- `uv 0.12.3 lock` / `uv sync --locked`
- Ruff lint/format
- mypy
- wheel/sdist build
- Snakemake 9 parser + toy DAG dry-run

## 尚未完成的验收

- 未在真实植物多基因组输入上完成全模块端到端运行。
- 未验证所有第三方工具在每个目标 HPC 平台上的资源需求。
- 未把来源模板的全部图件逐图复现。
- 未把 PlantCARE/CDD/WoLF PSORT 网页服务自动化。

因此当前标签必须保持 alpha。
