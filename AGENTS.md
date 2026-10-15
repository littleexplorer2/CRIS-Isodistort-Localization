# CRIS agent rules

进入子项目后还要读取该目录的 `agent.md`；使用方法读各自 `README.md`。

## 边界

- 永远只读：`experiment_data/`、`webpage_info/`、`output_compare/`、
  `ISODISTORT/resources/isobyu/`、仓库内 `GD/`。
- 与 CRIS 同级的 `GD/` 服从其自身规则；禁止修改 tianren 参考笔记本及
  `D:\OneDrive\...` 原始数据。
- 可改主体：`ISODISTORT/`（除 `resources/isobyu/`）、`ISODISTORT_VALIDATE/`、
  `ISOVIZ_INPUT/`。
- 三个子项目共用根目录内的实体 `.venv`，除非用户明确要求，不重建环境。普通
  Python 命令使用 `.\.venv\Scripts\python.exe`；本机从 OneDrive 内的解释器调用
  WSL 会触发 `Wsl/E_ACCESSDENIED`，需要 WSL 的 ISODISTORT 命令改用根目录
  `.\run_cris.ps1 <Python 参数>`，包括网页入口与相关测试；其依赖仍只从 `.venv`
  加载。项目仅使用该实体 `.venv`，不要另行创建环境。
- 根部署器与仓库编排测试统一放在 `tests_repo/`；三个子项目各自的 `tests/` 或
  `tests_dev/` 保持不变。运行 `.\run_cris.ps1 -m pytest tests_repo -q --tb=line`。
- 本机 GUI 快捷方式与可执行启动器统一放在 `local_tools/launchers/`，该目录内容不提交；
  `ISOVIZ` / `ISOVIZ_JAR` 环境变量仍可显式覆盖默认查找位置。

## DSH 诊断会话

- DSH 源码 checkout 为 `C:\Users\devou\deepseek-harness`。从 CRIS 根目录用
  `C:\Users\devou\deepseek-harness\node_modules\.bin\tsx.cmd --tsconfig
  C:\Users\devou\deepseek-harness\tsconfig.json
  C:\Users\devou\deepseek-harness\apps\cli\src\bin.ts --profile headless`
  启动，保证 Session 工作目录与只读诊断范围是本仓库；不使用 npm/npx 中的另一份
  DSH。多行提示词从 `prompt.md` 经 stdin 传入；不要把多行文本作为 `.cmd` 的位置
  参数，否则 Windows 包装器可能只传递第一行。
- 每次调用先建立
  `ISODISTORT/output/validation/dsh_sessions/<timestamp>-<slug>/`，至少保存
  `prompt.md`、`metadata.json`、原始 stdout 和 stderr；该目录不提交。提示词必须重申
  DSH 只读以及本文件中的受保护目录，禁止 DSH 修改或清理任何工作区文件。

## 文档职责

- 根 `README.md`：跨项目总览。
- 子项目 `README.md`：用户安装、使用、能力与限制。
- `ISODISTORT/resources/config/settings.yaml`：运行时默认值的事实来源。
- 三层模块目录、依赖方向与兼容边界见下节；可执行清单为
  `tests_repo/architecture_manifest.json`。
- `ISODISTORT/docs/DEVELOPMENT_PLAN.md`：唯一开发计划和待办，固定包含两个常驻分区
  “待修的程序漏洞”和“待执行的优化计划”；正文只维护当前与后续待办。
- `ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md`：当前验证缺口、在用证据入口与科学依据。
- 已完成的 debug、下载和修改历史无损压缩到 `ISODISTORT/docs/history/`；索引记录
  原路径、SHA-256 和恢复方法，归档须验证可恢复后再精简正文。有效操作说明和在用
  机器 manifest 保持可直接读取，不重复维护历史叙述。
- `ISODISTORT/docs/DOWNLOAD_CHECKLIST.md`：人工下载动作与目录。
- `ISODISTORT/docs/MANUAL_VALIDATION.md`：长时/手工验证命令。

## 三层架构基线（后续修改必须遵守）

- 调用方向固定为 `frontend → application → backend`。底层不依赖应用层或前端；
  应用层不依赖前端；前端只经应用层调用计算、文件工作流与运行环境查询。
- 每个模块必须有独立文件夹和 `__init__.py` 入口；模块内部可按职责拆文件。
  跨模块复用现有公开入口或其明确公开的子模块，不复制实现；状态通过调用合同传递。
  ISODISTORT 原有科学子模块的文件级接口登记在清单 `public_submodules` 中；
  少量既有下划线命名合同显式登记在 `legacy_scientific_interfaces`，仅保留原调用。
  未登记的实现文件和包内私有符号不得成为其他模块的新依赖；新功能从包公开入口提供。
- `backend` 保存科学算法、数据模型、格式读写和运行基础设施；`application` 保存
  会话、用例协调、事务/锁与结果合同；`frontend` 保存 HTTP/CLI/GUI 输入输出适配。
  `resources`、配置、文档、测试、输出不是额外的业务层。
- 新增、迁移或重命名模块时，同步可执行清单、调用方、打包配置、兼容入口、
  来源签名和对应文档。所有后续改动以以下已实现目录为基线。

```text
ISODISTORT/
  backend/
    wrappers/ models/ tables/ utils/ messages/
    input_cif/ method1/ method2/ method3/ method4/ export/
  application/
    session/                 # IsoDistort Python 会话与计算/导出协调
    workflows/               # 交互用例、revision、上传、导出任务与进度
  frontend/
    web/                     # HTTP、页面/静态资源、浏览器与标签页生命周期
ISODISTORT_VALIDATE/isodistort_validate/
  backend/
    config/ paths/ cif/ comparison/ results/
  application/
    single/ batch/ workspace/
  frontend/
    cli/
ISOVIZ_INPUT/isoviz_input/
  backend/
    amplitudes/ inputs/ validation/ paths/ launcher/ config_loader/
  application/
    pipeline/ input_selection/ environment/
  frontend/
    cli/ gui/ entrypoint/
```

- `ISODISTORT/backend/api/`、`ISODISTORT/features/`、`ISODISTORT/frontend/i18n/`
  及另两个包根的旧平铺模块只作兼容转发；不再加入实现，生产层不得经这些入口绕过
  分层。特别保留外部 GD 使用的 `backend.api` 与 `features.export.distortion_formats`；
  对应新入口为 `application.session` 与 `backend.export.distortion_formats`。
- 根 `setup_cris.py`、`run_cris.ps1` 与子项目 `main*.py` 保留安装/启动入口职责；
  不在入口复制领域算法或应用工作流。
- 架构门禁：`.\run_cris.ps1 -m pytest tests_repo/test_layer_boundaries.py -q --tb=line`。
  它检查模块登记、独立目录、单向导入与兼容层无业务实现；运行子项目规则指定的
  回归测试与打包验证后才能交付。静态架构检查不能替代科学结果或 GUI 验收。
