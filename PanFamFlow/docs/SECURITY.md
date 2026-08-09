# 安全与数据边界

- 不要把未公开 genome、FASTQ、样本身份、访问令牌或服务器凭据提交到 GitHub。
- `data/`, `work/`, `results/`, `logs/`, `.snakemake/` 默认被 `.gitignore` 排除。
- `scripts/install_dupgen.sh` 在目标目录存在时拒绝覆盖。
- CLI 使用参数列表调用 subprocess，不拼接 shell 命令。
- 用户提供的 `extra_args` 会传递给相应工具；只应使用可信配置文件。
- 公开 issue/PR 不应附带受限数据或路径中包含的敏感标识。
