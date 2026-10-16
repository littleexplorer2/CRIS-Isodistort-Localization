# ISOVIZ_INPUT agent rules

适用于本目录；先读根 `AGENTS.md`。使用与环境限制见 `README.md`，路径与启动器
查找见 `config/settings.yaml`。

## 职责与边界

- 只负责接收文件或粘贴文本形式的 GD 振幅 CSV 与子群 `.isoviz`、写入临时副本的
  `amp` 并可选启动 IsoVIZ；不搜索子群、不比较 CIF、不拟合振幅。
- 分层、公共入口与兼容约束执行根规则；`isoviz_input/` 内职责划分如下：
  - `backend/`：`amplitudes/` 负责 CSV/amp；`inputs/` 负责文本输入；
    `validation/` 负责静态不变量；`paths/` 负责路径；`launcher/` 负责启动适配；
    `config_loader/` 读取唯一配置 `config/settings.yaml`。
  - `application/`：`pipeline/` 是唯一准备、验证、临时写入与启动管线；
    `input_selection/` 统一文件/文本选择接口；`environment/` 提供启动器查询。
  - `frontend/`：`gui/` 和 `cli/` 处理展示与交互；`entrypoint/` 选择运行模式。
- GUI/CLI 必须调用同一应用管线，禁止另行实现解析、写入或启动逻辑。
- 根 `main.py` 仅引导至 `frontend.entrypoint`；默认 GUI，`--cli` 进入命令行。
  安装包命令 `isoviz-input` 使用同一入口。根 `main_requirement.py`、`.gitignore`
  和配置目录保留；wheel 中配置只能由 `setup.py` 从唯一配置生成，不另维护副本。
- 禁止修改或提交用户 CSV、实验衍射数据、原始 `.isoviz`、同级 GD 的参考笔记本；
  写入结果使用系统临时文件，`input_content/` 不作为产品输出。
- IsoVIZ `amp` 使用 **Best Model Parameter**，不是 Normalized Amplitude。
- 振幅约定、模式标签/顺序、CSV 到 `amp` 的映射变更属于科学语义修改，必须按根
  规则核对输入模型、模式基、单位及归一化合同；禁止用官网振幅或经验比例补齐、
  替换用户参数。历史振幅换算假设必须标注证据状态，不得作为已证明合同沿用。

## 实现门禁

- 列名别名、标签匹配和 `a1,a2,…` 顺序别名必须通用处理。
- 文件上传与粘贴内容必须经过同一解析器和静态门禁；源 CSV 与源 `.isoviz` 均只读。
- CLI 的 CSV 与 `.isoviz` 路径显式传递；启动器顺序为 yaml 环境变量，再按 yaml
  目录/名称搜索。没有 Java/启动器时不得报告窗口已打开；进程分发也不得表述为
  窗口已打开或数据已识别。

## 完成标准

- 修改代码或配置后，从 CRIS 根目录运行
  `.\.venv\Scripts\python.exe -m pytest ISOVIZ_INPUT\tests_dev -q --tb=line`；
  模块、入口、资源或打包变更必须验证入口与 wheel，依赖 `requirements-dev.txt`，
  构建测试必须在测试临时目录中进行。
- 默认门槛是 CSV 可解析、匹配报告正确、临时 `.isoviz` 可写、启动器路径可解析；
  只有用户要求时才启动 GUI。
- 真实 GUI 验收使用 `tests_dev/manual/validate_isoviz_gui.py --launch`，必须分别记录
  窗口打开与数据识别的人工确认结果。
