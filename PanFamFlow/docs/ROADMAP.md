# Roadmap

## v0.1.x：工程稳定

- 已实现真实生物学 benchmark 的 fail-closed 输入与预注册启动门。
- 完成 5–10 个真实水稻 assembled genomes 和一个人工核验目标家族的端到端 pilot。
- 记录每个 rule 的真实软件版本和命令清单。
- 增加 Snakemake end-to-end toy execution（使用小型 mock/binary fixtures）。
- 增加 GFF3 方言、gzipped input 和异常 ID 测试。
- 固化 OrthoFinder HOG node 审阅流程。
- 对 DupGen_finder 输出格式做更多版本 fixture 测试。

## v0.2：分析扩展

- MCScanX 独立 duplication/synteny backend。
- family domain architecture 与 MEME sequence logo。
- codeml branch/site/branch-site 模型。
- DESeq2 contrasts、样本设计矩阵与多重检验。
- expression 元数据分组、species-within normalization 和响应方向整合。
- 重点 absence 的 TBLASTN/miniprot/共线性复核模块。
- 同种材料的 syntenic pan-locus 构建，将 HOG occupancy 与材料级 locus PAV 分离。

## v0.3：图形和报告

- 完整实现来源模板中的组合图。
- 可配置、跨图固定且色盲友好的颜色映射。
- family tree + gene structure + motif + expression 的整合图。
- 可选交互式 HTML 数据浏览器。

## 发布门槛

从 alpha 升级到 beta 前至少需要：

1. 两套独立真实数据完成端到端运行。
2. 所有主要模块有 fixture 或 integration test。
3. CI 在固定 lockfile 下通过 lint/type/test/build/dry-run。
4. 输出 schema、目标家族边界和 HOG/clade/pan-locus 术语经人工审阅。
5. 已知高风险模块有失败注入测试。
