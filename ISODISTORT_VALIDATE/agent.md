# ISODISTORT_VALIDATE agent rules

适用于本目录；先读根 `AGENTS.md`。使用与限制见 `README.md`，路径和容差见
`config/settings.yaml`。本文件不记录项目进度或验证报告。

## 职责与边界

- 本工具只比较本地/官网 CIF 的晶体学语义，不生成结构、不改黄金 CIF、不调用
  WSL/iso/VESTA/IsoVIZ。
- 核心职责：单对比较 `compare_cif.py`，批量比较 `batch_compare.py`，固定目录
  与配对 `compare_paths.py`，配置 `config_loader.py`，唯一入口 `main.py`。
- `compare/`、`output_compare/`、`experiment_data/`、`webpage_info/` 只读对照且
  不提交；不要访问 ISODISTORT 私有对象或复制其计算逻辑。
- 配对身份是 `item/` 与 `true/` 的相同相对路径。默认做语义比较；
  `--strict` 仅用于排版调试，`byte_exact=false` 不单独导致失败。

## 实现门禁

- 不得按文件名、材料或 IR 硬编码 PASS。原子配对应基于物种、周期晶格度量、
  占据率/磁矩和空间群语义，并明确处理行序。
- 路径、默认容差只来自 yaml；接口/退出码变化须同步 CLI、菜单、README 和测试。
- 修复后用最小反例和成对 CIF 复验，不能修改 `true/` 来迁就算法。

## 完成标准

- 运行 `.\.venv\Scripts\python.exe -m pytest ISODISTORT_VALIDATE\tests_dev -q --tb=line`。
- 改比较语义时再运行 `main.py compare` / `main.py batch`；缺少成对用户文件时
  明确说明，不伪造样本。
- 执行 `git diff --check`，清理缓存/日志；行为改 README，默认值改 yaml。
