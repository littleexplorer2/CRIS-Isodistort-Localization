# CRIS agent rules

进入子项目后还要读取该目录的 `agent.md`；使用方法读各自 `README.md`。

## 边界

- 永远只读：`experiment_data/`、`webpage_info/`、`output_compare/`、
  `ISODISTORT/isobyu/`、仓库内 `GD/`。
- 桌面 `GD（未同步git）` 服从其自身规则；禁止修改 tianren 参考笔记本及
  `D:\OneDrive\...` 原始数据。
- 可改主体：`ISODISTORT/`（除 `isobyu/`）、`ISODISTORT_VALIDATE/`、
  `ISOVIZ_INPUT/`。
- 三个子项目共用根目录内的实体 `.venv`，除非用户明确要求，不重建环境。普通
  Python 命令使用 `.\.venv\Scripts\python.exe`；本机从 OneDrive 内的解释器调用
  WSL 会触发 `Wsl/E_ACCESSDENIED`，需要 WSL 的 ISODISTORT 命令改用根目录
  `.\run_cris.ps1 <Python 参数>`，包括 `main_web.py`、`main_terminal.py` 和相关测试；
  其依赖仍只从 `.venv` 加载。两个入口会拒绝已知必失败的直接 `.venv` 启动。

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
- `config/settings.yaml`：运行时默认值的事实来源。
- `ISODISTORT/docs/DEVELOPMENT_PLAN.md`：唯一开发计划和待办。
- `ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md`：已修问题、验证结论、机器报告索引。
- `ISODISTORT/docs/DOWNLOAD_CHECKLIST.md`：人工下载动作与目录。
- `ISODISTORT/docs/MANUAL_VALIDATION.md`：长时/手工验证命令。
