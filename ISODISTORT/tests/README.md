# ISODISTORT 功能检查

本目录按**被检验代码的责任模块**组织自动回归。每个生产模块只保留一个 `test_<模块>.py`；
用例及辅助函数直接在新文件中维护，不再加载旧测试脚本。

## 一条命令检查

在 PowerShell 中进入 CRIS 根目录，例如 `cd E:\CRIS`。Windows 下用现有 `run_cris.ps1`
调用实体 `.venv` 的依赖并连接 WSL；不需要另建测试环境。

完整检查：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests -q --tb=line
```

只检查某模块，例如 Method 2、会话或网页：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests/backend/test_method2.py -q --tb=line
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests/application/test_session.py -q --tb=line
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests/frontend/test_web.py -q --tb=line
```

将上面的文件名换为下列模块名即可。也可以把路径写成 `ISODISTORT/tests/backend` 来检查整个层。
不直接双击 `test_*.py`；它们由 pytest 收集并管理 fixture、参数化及运行隔离。

## 目录和归属

- `backend/test_wrappers.py`：本地计算套件、WSL 目录与进程、缓存。
- `backend/test_models.py`：共享模式、来源身份、相路径等数据合同。
- `backend/test_tables.py`：空间群、表示、k 点和元数据表。
- `backend/test_utils.py`：晶格与仿射工具、解析、配置、自检。
- `backend/test_messages.py`：共享英文文案合同。
- `backend/test_input_cif.py`：结构读取、对称性、晶胞变换与站点映射。
- `backend/test_method1.py`：特殊 k 搜索及支持边界。
- `backend/test_method2.py`：指定候选、调制与模式映射。
- `backend/test_method3.py`：任意 k/格、仿射搜索、耦合路线。
- `backend/test_method4.py`：分解、应变与生成回代。
- `backend/test_export.py`：CIF、IsoVIZ、TOPAS 和模式说明的格式合同。
- `application/test_session.py`：计算会话、状态、批量导出及发布事务。
- `application/test_workflows.py`：用例、上传、revision 和导出任务。
- `frontend/test_web.py`：HTTP、页面、浏览器及生命周期。

跨层架构检查、独立审计工具的回归及准备/诊断/比较脚本直接放在 `tests/` 根目录。
`conftest.py` 提供收集环境；`support_*.py` 在根目录保存跨文件共享的夹具、辅助构造和稳定路径，
避免复制或导入另一个测试脚本。测试数据放在 `fixtures/`（已有生成样本位于 `cifs_30/`），不是待执行脚本。
`validation/` 是已有的只读比较**实现包**，保留其三层模块和公共入口；它不是另一套测试脚本目录，
也不会启动本地生产计算。禁止把官网答案送回生产链。

## 手工工具与科学验收

原 `tests/manual/` 工具已移到这里，入口例如 `tests/verify_parametric_k_grid.py`、
`tests/validate_method_outputs.py`、`tests/prepare_method4_inputs.py`。
操作前阅读 [手工验证说明](../docs/MANUAL_VALIDATION.md)，使用新批次输出目录；不要直接运行
依赖已撤销冻结记录的历史长时批次。只读比较入口仍是 `-m tests.validation compare` / `batch`。

成功退出或自动回归通过不等于全材料科学验收。跳过项、缺失输入及未执行的 GUI/长时检查应单独记录；
具体科学合同、容差和本次迁移证据见 [验证报告](../docs/BUGFIX_VALIDATION_REPORT.md)。
