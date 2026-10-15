# ISOVIZ_INPUT agent rules

适用于本目录；先读根 `AGENTS.md`。使用与环境限制见 `README.md`，路径与启动器
查找见 `config/settings.yaml`。

## 职责与边界

- 只负责接收文件或粘贴文本形式的 GD 振幅 CSV 与子群 `.isoviz`、写入临时副本的
  `amp` 并可选启动 IsoVIZ；不搜索子群、不比较 CIF、不拟合振幅。
- `isoviz_input/` 固定三层；每个职责模块独立成目录，对外仅使用该目录
  `__init__.py` 的公开接口，不跨模块引用实现文件或私有符号：
  - `backend/`：`amplitudes/` 负责 CSV/amp；`inputs/` 负责文本输入；
    `validation/` 负责静态不变量；`paths/` 负责路径；`launcher/` 负责启动适配；
    `config_loader/` 读取唯一配置 `config/settings.yaml`。
  - `application/`：`pipeline/` 是唯一准备、验证、临时写入与启动管线；
    `input_selection/` 统一文件/文本选择接口；`environment/` 提供启动器查询。
  - `frontend/`：`gui/` 和 `cli/` 处理展示与交互；`entrypoint/` 选择运行模式。
- 调用方向为 `frontend -> application -> backend`。前端不得直接依赖底层或
  包根兼容文件；应用层不得导入前端；底层不得导入应用层或前端。不得在 GUI/CLI
  重复解析、写入或启动逻辑；包根同名文件仅作旧导入兼容，不新增业务实现。
- 根 `main.py` 仅引导至 `frontend.entrypoint`；默认 GUI，`--cli` 进入命令行。
  安装包命令 `isoviz-input` 使用同一入口。根 `main_requirement.py`、`.gitignore`
  和配置目录保留；wheel 中配置只能由 `setup.py` 从唯一配置生成，不另维护副本。
- 禁止修改或提交用户 CSV、实验衍射数据、原始 `.isoviz`、同级 GD 的参考笔记本；
  写入结果使用系统临时文件，`input_content/` 不作为产品输出。
- IsoVIZ `amp` 使用 **Best Model Parameter**，不是 Normalized Amplitude。

## 实现门禁

- 列名别名、标签匹配和 `a1,a2,…` 顺序别名必须通用处理。
- 文件上传与粘贴内容必须经过同一解析器和静态门禁；源 CSV 与源 `.isoviz` 均只读。
- CLI 的 CSV 与 `.isoviz` 路径显式传递；启动器顺序为 yaml 环境变量，再按 yaml
  目录/名称搜索。没有 Java/启动器时不得报告窗口已打开；进程分发也不得表述为
  窗口已打开或数据已识别。

## 完成标准

- 运行 `.\.venv\Scripts\python.exe -m pytest ISOVIZ_INPUT\tests_dev -q --tb=line`；
  入口与 wheel 验证依赖 `requirements-dev.txt`，构建测试在测试临时目录中进行。
- 默认门槛是 CSV 可解析、匹配报告正确、临时 `.isoviz` 可写、启动器路径可解析；
  只有用户要求时才启动 GUI。
- 真实 GUI 验收使用 `tests_dev/manual/validate_isoviz_gui.py --launch`；进程分发成功
  不等于窗口已打开或数据已识别，必须分别记录人工确认结果。
