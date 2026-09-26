# ISOVIZ_INPUT agent rules

适用于本目录；先读根 `AGENTS.md`。使用与环境限制见 `README.md`，路径与启动器
查找见 `config/settings.yaml`。本文件不记录项目进度或运行报告。

## 职责与边界

- 只负责读取 GD 振幅 CSV、写入子群 `.isoviz` 的 `amp` 并可选启动 IsoVIZ；
  不搜索子群、不比较 CIF、不拟合振幅。
- 唯一职责模块：CSV/amp 为 `amplitudes.py`，路径为 `paths.py`，启动器为
  `launcher.py`，配置为 `config_loader.py`，唯一入口 `main.py`。
- 禁止修改或提交用户 CSV、实验衍射数据、原始 `.isoviz`、桌面 GD 参考笔记本；
  写入结果使用系统临时文件，`input_content/` 不作为产品输出。
- IsoVIZ `amp` 使用 **Best Model Parameter**，不是 Normalized Amplitude。

## 实现门禁

- 不按材料、IR、模式名或某次振幅表硬编码；通用处理列名别名、标签匹配和
  `a1,a2,…` 顺序别名。
- CSV 与 `.isoviz` 路径显式传递；启动器顺序为 yaml 环境变量，再按 yaml
  目录/名称搜索。没有 Java/启动器是环境限制，不得假装窗口已打开。
- 接口、列名、查找规则变化须同步 README、yaml、提示文案和 fixtures 测试。

## 完成标准

- 运行 `.\.venv\Scripts\python.exe -m pytest ISOVIZ_INPUT\tests_dev -q --tb=line`。
- 默认门槛是 CSV 可解析、匹配报告正确、临时 `.isoviz` 可写、启动器路径可解析；
  只有用户要求时才启动 GUI。
- 执行 `git diff --check`，清理缓存/日志；行为改 README，默认值改 yaml。
