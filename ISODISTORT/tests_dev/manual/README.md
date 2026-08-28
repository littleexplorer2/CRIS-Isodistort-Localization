# tests_dev/manual — 手工/长时验证脚本

本目录不会被 pytest 收集；用于生成 CIF 与批量回归。

```text
python tests_dev/manual/make_cifs_30.py
python tests_dev/manual/fetch_cod_cifs.py
python tests_dev/manual/run_web.py spotcheck
python tests_dev/manual/run_web.py m134
python tests_dev/manual/run_web.py method2_ld
python tests_dev/manual/run_batch.py cif30
python tests_dev/manual/run_batch.py external
python tests_dev/manual/run_batch.py terminal
```

请使用项目 `.venv` 中的 Python。CIF 数据仍在 `tests_dev/cifs_30/` 与 `cifs_external/`。

## 用户需要配置 / 放置的路径

本目录脚本**沿用**上层 ISODISTORT 配置，一般不必单独改路径：

| 项目 | 说明 |
| --- | --- |
| **母相 / 批量 CIF** | 默认读 `tests_dev/cifs_30/`、`cifs_external/`，以及仓库 `experiment_data/`（只读）。换样本时把 CIF 放进对应测试数据目录，或改脚本参数里的路径 |
| **iso / WSL / output** | 与正式程序相同：`config/settings.yaml`、`isobyu/`，见 [ISODISTORT/README.md](../../README.md) §3.4 |
| **网页 spotcheck** | 依赖本机已能启动 `main_web.py`（端口见 `runtime.web_port`） |

跨项目总表见仓库根 [README.md](../../../README.md)。
