# HPC / SLURM 执行

## 1. 创建 engine

```bash
mamba env create -f environment-slurm.yaml
```

配置：

```yaml
run:
  cores: 64
  jobs: 100
  engine_runner: mamba
  engine_env: panfamflow-engine-slurm
  use_conda: true
  profile: /absolute/path/to/PanFamFlow/profiles/slurm
```

先验证和 dry-run：

```bash
uv run panfamflow validate -c config.yaml
uv run panfamflow run -c config.yaml --dry-run
```

## 2. 后台提交

建议创建非覆盖式日志目录：

```bash
mkdir -p submit_logs
```

示例 `submit_panfamflow.sbatch`：

```bash
#!/usr/bin/env bash
#SBATCH --job-name=panfamflow
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=submit_logs/%x.%j.out
#SBATCH --error=submit_logs/%x.%j.err

set -euo pipefail
cd /path/to/project
uv run panfamflow run -c config.yaml
```

提交：

```bash
sbatch submit_panfamflow.sbatch
```

不要在登录节点直接执行长时间全流程。

## 3. 分级退避监控

推荐监控节奏：

1. 前 3 分钟：每分钟检查一次。
2. 进入稳定运行后：每 10 分钟检查一次。
3. 连续 3 次状态正常后：每 30 分钟检查一次。
4. 发现失败、长时间无输出或资源异常时恢复高频检查。

最小检查命令：

```bash
squeue -j <JOB_ID> -o '%.18i %.9P %.20j %.2t %.10M %.10l %.6D %R'
sacct -j <JOB_ID> --format=JobID,State,Elapsed,MaxRSS,AllocCPUS,ExitCode
find logs -type f -mmin -30 -print | sort
tail -n 100 submit_logs/panfamflow.<JOB_ID>.err
```

异常报告至少记录：时间点、job/rule、日志证据、可能原因排序、最小修复动作和复现命令。

## 4. 失败恢复

Snakemake 默认启用 `--rerun-incomplete`。修正配置或环境后重新提交相同命令，已完成且输入未变化的结果不会重算。

锁残留：

```bash
uv run panfamflow run -c config.yaml --unlock
```

仅在确认没有其他 Snakemake 实例操作同一工作目录后使用。

## 5. 资源覆盖

`profiles/slurm/config.yaml` 提供默认资源和重点规则覆盖。可根据真实测试修改：

```yaml
set-resources:
  orthofinder:
    mem_mb: 128000
    runtime: 5760
  family_phylogeny:
    mem_mb: 64000
    runtime: 2880
```

正式批量运行前，应使用代表性物种/样本做 pilot，记录 MaxRSS 和 elapsed，再调整资源。
