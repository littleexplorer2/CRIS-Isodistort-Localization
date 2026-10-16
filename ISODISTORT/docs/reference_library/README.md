# ISODISTORT 理论参考资料库

本库为 CRIS 的科学计算提供可追踪的文献、数学合同和证明入口。检索截止日为
2026-10-15，使用 Exa、Consensus、Undermind，并以作者官网、原始论文、出版社及
IUCr 标准复核。检索平台的摘要或回答本身不作为数学依据。

## 两份论文

- **补充证明**：[LaTeX 源文件](latex/unsupported_algorithm_proofs.tex) ·
  [阅读 PDF](latex/unsupported_algorithm_proofs.pdf)。讨论缺少直接实现依据的环节：
  精确稳定子与有理见证、平移格与谐波、周期匹配、模式身份、斜晶格应变、
  非公度微观模式构造和搜索完备性的条件。
- **完整能力的理论框架**：[LaTeX 源文件](latex/isodistort_theoretical_foundations.tex) ·
  [阅读 PDF](latex/isodistort_theoretical_foundations.pdf)。按官方帮助的能力范围解释
  空间群、物理表示、Method 1–4、调制、幅度、畴、子群树、衍射、导出与精修；
  包含定义、假设、定理、证明和参考文献。

两份文档是项目技术论证，尚未经独立同行评审。它们证明明确假设下的数学合同，
不声称形式化认证全部官网程序、外部二进制、表数据、GUI 或浮点容差分支。
已发表分类算法作为带来源的前提；无法确定的官网约定和实现条件单独列出，
不能用“主题已覆盖”代替“全部实现已验证”。

## 查找与审核

- [注释书目](notes/ANNOTATED_SOURCES.md)：逐条说明支持内容、实际核验材料及限制。
- [能力覆盖清单](coverage.json)：26 个主题与来源、证明标签、本地模块的对应关系。
- [信源和数学审核](audit/QUALITY_REVIEW.md)：来源取舍、冲突处理、未解决的证明义务。
- [理论与代码审查](audit/THEORY_CODE_AUDIT.md)：实际反例、已修复路径、能力限制及
  [逐项实现审查清单](audit/implementation_audit.json)；不将文献覆盖当作软件认证。
- [书目主数据](sources.json)：50 条来源的唯一可编辑元数据。
- [统一 BibTeX](latex/references.bib)：由主数据生成，两篇论文共用。
- [检索记录](audit/search_log.json)、[下载清单与 SHA-256](audit/download_manifest.json)、
  [当前文档检查结果](audit/validation.json)：分别记录发现、获取和验收证据。

建议先读官方帮助与 Campbell 2006，再读 Stokes 1991、逆 Landau 2002、
ISO-IR 2013 和超空间 2017；随后按需要读模式晶体学综述、磁超空间综述以及独立
AMPLIMODES/spglib 方法。扫描 PDF、仅摘要、出版社入口与取得的正文节选在书目中分开标记。

## 文件布局与获取边界

`pdf/` 保存合法公开获取的作者/机构 PDF，`snapshots/` 保存网页及抓取文本；
`notes/` 保存注释，`latex/` 保存两篇论文、共享排版和书目，`audit/` 保存来源与审核证据。
本次本地取得 26 份参考 PDF。其余来源保留权威链接及获取状态；403、证书错误或
未取得全书正文都没有被标为下载成功，也没有绕过访问控制。

第三方全文、网页快照、提取文本、排版缓存和渲染图在 `.gitignore` 中按路径排除；
它们仍在本地资料目录中。主数据、链接、注释、证明源文件、两份自写 PDF 和审核结果
可纳入版本管理。从远程重新获取仓库时，用下列脚本恢复可公开取得的资料；
同一 URL 内容可能更新，重新获取后的 SHA-256 必须复核。

## 更新与复现

从 CRIS 根目录运行；使用已有 `.venv`，不安装或重建项目环境：

```powershell
.\.venv\Scripts\python.exe ISODISTORT/docs/reference_library/audit/acquire_sources.py
# 也可只获取某些条目，例如 --ids Campbell2006 IsoIr
.\.venv\Scripts\python.exe ISODISTORT/docs/reference_library/audit/build_catalogue.py
& ISODISTORT/docs/reference_library/latex/build.ps1
.\.venv\Scripts\python.exe ISODISTORT/docs/reference_library/audit/check_library.py
```

编译需要现有 TeX Live 的 `latexmk`、XeLaTeX、BibTeX 和 Fandol 中文字体。
`build.ps1 -Latexmk 'D:\texlive\2026\bin\windows\latexmk.exe'` 可指定本机路径。
PDF 提取/检查需要 `pdftotext` 和 `pdfinfo` 在 PATH 上；获取脚本对扫描件不伪造 OCR。
数学例证检查另运行 `audit/check_mathematical_examples.py`，结果仅校验文中例证，
不替代一般证明或 CRIS 科学回归。排版中间文件和 TeX 缓存均保存在本资料库内。

新增计算或修改约定时，先更新相应来源/证明与 `coverage.json`，明确坐标、单位、
输入域、列身份及不变量，再在原模块做实现和科学验收。未来验收必须使用当时的源码签名；
本资料库中的源码定位摘要不是软件测试通过证据。
