# CRIS agent rules

本文件只保存长期稳定的仓库级约束。进入子项目后还要读取该目录的
`agent.md`；使用方法读各自 `README.md`，不要把进度或报告追加到本文件。

## 边界

- 永远只读：`experiment_data/`、`webpage_info/`、`output_compare/`、
  `ISODISTORT/isobyu/`、仓库内 `GD/`。
- 桌面 `GD（未同步git）` 服从其自身规则；禁止修改 tianren 参考笔记本及
  `D:\OneDrive\...` 原始数据。
- 可改主体：`ISODISTORT/`（除 `isobyu/`）、`ISODISTORT_VALIDATE/`、
  `ISOVIZ_INPUT/`。不要用修改黄金数据来让测试通过。
- 三个子项目共用根目录 `.venv`；运行 Python 一律用
  `.\.venv\Scripts\python.exe`，除非用户明确要求，不重建环境。

## 工作方式

1. 先读本文件、目标子项目 `agent.md`、`git status`，再按该文件的文档路由
   读取最少必要资料。
2. 修通用算法与接口，不按材料名、IR、文件名或官网个例硬编码答案。
3. 保持唯一责任模块；变更接口时同步全部调用者、测试、配置和用户文档。
4. 先用最小反例复现，再实现并跑按风险匹配的测试；事实不够时明确写
   `inconclusive`，不能把“程序运行结束”当作正确性证明。
5. 只在用户明确要求时提交/推送。提交前执行相关测试、`git diff --check`，
   清除缓存、日志、临时下载与生成物；不得提交受保护数据或用户输入。

## 文档职责

- 根 `README.md`：跨项目总览。
- 子项目 `README.md`：用户安装、使用、能力与限制。
- `config/settings.yaml`：运行时默认值的事实来源。
- `ISODISTORT/docs/DEVELOPMENT_PLAN.md`：唯一开发计划和待办。
- `ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md`：已修问题、验证结论、机器报告索引。
- `ISODISTORT/docs/DOWNLOAD_CHECKLIST.md`：人工下载动作与目录。
- `ISODISTORT/docs/MANUAL_VALIDATION.md`：长时/手工验证命令。

进度、批次数字、bug 历史、下载状态和运行报告只能写入上述对应文档；
`AGENTS.md` / `agent.md` 仅保留规则与引用。
