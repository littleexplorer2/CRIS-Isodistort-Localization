# ISODISTORT_VALIDATE agent rules

适用于本目录；先读根 `AGENTS.md`。使用与限制见 `README.md`，路径和容差见
`config/settings.yaml`。

## 职责与边界

- 本工具只比较本地/官网 CIF 的晶体学语义，不生成结构、不改黄金 CIF、不调用
  WSL/iso/VESTA/IsoVIZ。
- 可执行入口 `main.py` 只转发 `isodistort_validate.frontend.cli.main`；安装后的
  `isodistort-validate` 命令使用同一入口。
- `compare/`、`output_compare/`、`experiment_data/`、`webpage_info/` 只读对照且
  不提交；不要访问 ISODISTORT 私有对象或复制其计算逻辑。
- 配对身份是 `item/` 与 `true/` 的相同相对路径。默认做语义比较；
  `--strict` 仅用于排版调试，`byte_exact=false` 不单独导致失败。

## 三层结构与调用规则

- 源码唯一实现位于 `isodistort_validate/` 的三个层级，各职责模块独立成目录，
  模块的 `__init__.py` 是公共调用入口；后续修改沿用以下职责划分。
- 底层 `backend/`：`config/` 加载 yaml 与路径基准；`paths/` 固定目录、配对和
  路径诊断；`cif/` 只读解析与哈希；`comparison/` 原子匹配、晶格和元数据比较；
  `results/` 比较结果契约。底层不得导入应用层或前端。
- 应用层 `application/`：`single/` 编排单对文件检查、解析、比较和哈希校验；
  `batch/` 编排目录配对、manifest 与批量汇总；`workspace/` 向前端提供配置、
  配对状态和默认文件选择。仅调用底层或同层模块公共入口，不解析 CLI 参数。
- 前端 `frontend/cli/`：`menu.py` 交互菜单与分发；`compare.py`、`batch.py`
  各自处理参数和文本/JSON 输出。前端通过应用层调用，不直接导入底层。
- 跨模块只能导入目录公共入口，不直接导入其它模块的 `service.py` 或私有对象。
  同一职责目录内部允许使用实现文件；不要复制比较算法或另建配置来源。
- 旧 `compare_cif.py`、`batch_compare.py`、`compare_paths.py`、`config_loader.py`
  只保留导入/模块别名兼容转发，新生产代码不得依赖它们；旧路径与配置模块别名
  保持同一对象，避免临时覆盖路径后存在两份状态。
- `config/settings.yaml` 是唯一配置源；`setup.py` 仅在构建时生成 wheel 内资源，
  不在源码中维护第二份 yaml。安装包的资源路径不得改变默认比较目录的相对语义。
- `tests_dev/test_architecture.py` 检查层级方向、模块公共入口、兼容转发和导入隔离。

## 实现门禁

- 原子配对应基于物种、周期晶格度量、占据率/磁矩和空间群语义，并明确处理行序。
- 路径和默认容差只来自 yaml。

## 完成标准

- 运行 `.\.venv\Scripts\python.exe -m pytest ISODISTORT_VALIDATE\tests_dev -q --tb=line`。
- wheel 验证使用 `requirements-dev.txt` 中的构建依赖，必须离开源码目录验证
  安装包入口、配置资源与单对/批量调用；缺少构建依赖时明确记录跳过。
- 改比较语义时再运行 `main.py compare` / `main.py batch`；缺少成对用户文件时
  说明该项未执行。
