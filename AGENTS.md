# CRIS agent rules

进入子项目后还要读取该目录的 `agent.md`；使用方法读各自 `README.md`。

## 边界

- 永远只读：`experiment_data/`、`webpage_info/`、`output_compare/`、
  `ISODISTORT/isobyu/`、仓库内 `GD/`。
- 桌面 `GD（未同步git）` 服从其自身规则；禁止修改 tianren 参考笔记本及
  `D:\OneDrive\...` 原始数据。
- 可改主体：`ISODISTORT/`（除 `isobyu/`）、`ISODISTORT_VALIDATE/`、
  `ISOVIZ_INPUT/`。
- 三个子项目共用根目录 `.venv`；运行 Python 一律用
  `.\.venv\Scripts\python.exe`，除非用户明确要求，不重建环境。

## 文档职责

- 根 `README.md`：跨项目总览。
- 子项目 `README.md`：用户安装、使用、能力与限制。
- `config/settings.yaml`：运行时默认值的事实来源。
- `ISODISTORT/docs/DEVELOPMENT_PLAN.md`：唯一开发计划和待办。
- `ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md`：已修问题、验证结论、机器报告索引。
- `ISODISTORT/docs/DOWNLOAD_CHECKLIST.md`：人工下载动作与目录。
- `ISODISTORT/docs/MANUAL_VALIDATION.md`：长时/手工验证命令。
