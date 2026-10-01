# ISODISTORT_VALIDATE agent rules

适用于本目录；先读根 `AGENTS.md`。使用与限制见 `README.md`，路径和容差见
`config/settings.yaml`。

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

- 原子配对应基于物种、周期晶格度量、占据率/磁矩和空间群语义，并明确处理行序。
- 路径和默认容差只来自 yaml。

## 完成标准

- 运行 `.\.venv\Scripts\python.exe -m pytest ISODISTORT_VALIDATE\tests_dev -q --tb=line`。
- 改比较语义时再运行 `main.py compare` / `main.py batch`；缺少成对用户文件时
  说明该项未执行。
