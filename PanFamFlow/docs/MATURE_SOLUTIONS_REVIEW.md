# GitHub 成熟方案评估

评估目标不是寻找“代码看起来相似”的仓库，而是判断是否存在可直接满足以下组合需求的维护方案：植物多基因组、canonical transcript、HMM+BLAST family discovery、OrthoFinder 3 HOG、泛基因四分类、复制方式、Ka/Ks、promoter、expression、one-config、模块选择、rule-specific Conda 和现代 Python packaging。

## 1. Snaketool

仓库：`beardymcjohnface/Snaketool`

优点：提供 Snakemake workflow + Python CLI 的项目模板思想，适合借鉴 launcher、配置和包结构。

不足：它是通用工程模板，不提供本项目所需的植物基因家族和泛基因分析逻辑；旧模板结构也不直接覆盖 uv/Pydantic v2/Snakemake 9/HOG 需求。

决策：借鉴架构，不直接 fork。

## 2. orthosnake

仓库：`paraslonic/orthosnake`

优点：用 Snakemake 串联 annotation 与 OrthoFinder，结构简单。

不足：以 Prokka 和 `.fna` 为主，明显偏原核；覆盖范围主要是 annotation + orthology，不含植物 GFF3 canonical transcript、HOG node、复制、Ka/Ks、promoter 和 expression；公开更新时间较旧。

决策：不作为基线。

## 3. smsk_selection

仓库：`jlanga/smsk_selection`

优点：包含 OrthoFinder、正选择和 Conda/Snakemake 思路，功能比 orthosnake 更广。

不足：以 transcriptome/TransDecoder/FastCodeML 路线为主，目录和依赖结构较旧；不满足 one-config、OrthoFinder 3 HOG、植物复制分类、promoter 和表达整合要求。

决策：方法模块可参考，不直接 fork。

## 4. 官方算法软件

以下软件应直接复用而不是重写算法：

- OrthoFinder
- HMMER
- BLAST+
- MAFFT
- ClipKIT
- IQ-TREE
- AGAT/gffread
- DupGen_finder
- PAL2NAL/KaKs_Calculator
- MEME Suite/FIMO
- fastp/HISAT2/StringTie
- BUSCO

## 5. 最终方案

采用“自建工程层 + 复用发表算法软件”的组合：

- 工程层由 PanFamFlow 管理配置、DAG、ID、QC、输出 schema 和报告。
- 算法层由成熟软件执行。
- 与第三方工具的接口通过 TSV/FASTA/GFF3 和规则级 Conda 环境隔离。

这样比直接延续旧仓库更容易控制技术债务、方法更新、测试覆盖和跨项目复用。
