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
