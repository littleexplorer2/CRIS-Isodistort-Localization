# 说明性文档可恢复归档

## 2026-10-20：Method 2 母相上下文及新增列表登记前快照

[method2-context-before-20261020.zip](method2-context-before-20261020.zip)
保存第 9 版冻结记录及修订前清单、计划、报告、索引、机器记录和定向审计脚本，共 14 个文件。
ZIP SHA-256：`2ab2fa7b1ae73f4c85a639492bad54ded9d9ef829e58a7918dc1a9bdf46c9281`。
条目名为原仓库相对路径；`_RECOVERY_SHA256.json` 保存原路径、长度和 SHA-256。
已校验 ZIP CRC 及逐项解压字节。恢复到新的空目录后按该索引验证，
不要覆盖当前记录；旧“6 个待下载查询”状态由后继记录取代，不能沿用为当前状态。
官网页面及资源原件没有被改写或删除。

## 2026-10-20：Method 2 定向下载范围修订前快照

[method2-targeted-before-20261020.zip](method2-targeted-before-20261020.zip)
保存原 102 查询清单、14 个改动前的说明/manifest/索引/审计脚本文件及第 8 版冻结记录。
ZIP SHA-256：`c692330a6351d6aab7be3a2c8c80298855b935ab92fee702e6457286c43edbaf`。
条目名是原仓库相对路径；`_RECOVERY_SHA256.json` 保存原路径、长度、SHA-256。
已逐项读取解压字节并比对哈希，CRC 检查通过。恢复时解压到新的空目录，
按 `_RECOVERY_SHA256.json` 校验，不能覆盖当前记录或复活旧下载范围/运行签名。
当前范围只由 [下载清单](../DOWNLOAD_CHECKLIST.md) 及后继冻结记录决定。

当前待办见 [开发计划](../DEVELOPMENT_PLAN.md)。这里保存精简前工作树的完整原始字节，
包括当时未提交的文档修改；压缩文件是历史快照，不作为当前执行计划。
机器 manifest、原始官网资料与机器验证产物未被本次压缩改写。

## 2026-10-20：连续修复与架构复审前的记录

[architecture-loop-before-20261020.zip](architecture-loop-before-20261020.zip)
保存当时工作树中的 ISODISTORT README、计划、报告、本索引、架构合同、能力映射、
实现审计及第 7 版冻结记录，共 8 个文件，另附 `_SHA256.json`。
SHA-256：`18ab5b243671ab7f088ceb0d9c4fc6004f261b65434e6ed94e5ea7a621653cb5`。
条目名为仓库相对路径，已核对 ZIP CRC 与每项恢复字节/SHA-256。
恢复时使用下文 Python 方法解压到独立目录；不要覆盖当前记录或把旧检查结果当作
新源码的全量验证结果。

## 2026-10-20：Method 2 查询归属修复前的记录

[method2-query-before-20261020.zip](method2-query-before-20261020.zip)
保存当时工作树中的 ISODISTORT README、计划、报告、本索引、架构合同、能力映射、
实现审计及第 6 版冻结记录，共 8 个文件，另附 `_SHA256.json`。
SHA-256：`a311a113b71fc8c1f0dedd76021ce76e4f99f96fd296653f6ef586803647f01c`。
条目名为仓库相对路径；ZIP CRC 及每项恢复字节/SHA-256 均已核对。
按下文 Python 方法恢复到独立目录，不要覆盖当前文档。该快照保存修前计划与证据，
不代表修后源码已经完成全量或科学验证。

## 2026-10-20：三个已复现问题修复前的记录

[architecture-bugs-before-20261020.zip](architecture-bugs-before-20261020.zip)
保存当时工作树中的 ISODISTORT README、计划、报告、本索引、架构合同、能力映射、
实现审计及第 5 版冻结记录，共 8 个文件，另附 `_SHA256.json`。
SHA-256：`a80879d7c05a205e6c3abe14472b26adb3d776dc89d09acf092ca260c629e592`。
条目名为仓库相对路径，已逐条解压到内存核对字节/SHA-256，并通过 ZIP CRC 检查。
按下文 Python 方法恢复到独立目录，对照 `_SHA256.json`，不要覆盖当前文件；
归档保留旧待办状态与证据，不作为新源码的验证结论。

## 2026-10-20：剩余架构修改前的支撑记录

[architecture-followup-before-20261020.zip](architecture-followup-before-20261020.zip)
保存修订前的 ISODISTORT README、开发计划、验证报告、本索引、架构合同、
能力映射、实现审计及第 4 版清单冻结记录，共 8 个原始文件，另附 `_SHA256.json`。
SHA-256：`c65c83b359070f0b8bce50d9dd322dc587bc8347d64163af08099153125fd116`。
ZIP 条目名为仓库相对路径；创建后已逐条解压到内存核对字节与 SHA-256，CRC 自检通过。
按下文 Python 恢复方法解压到独立目录，并对照 `_SHA256.json`；不要覆盖当前文档。
这些记录保留当时的延期事项和旧验证边界，不代表修改后源码已完成全量验证。

## 2026-10-09 快照

- [四份说明文档](docs-before-compaction-20261009.zip)：75,079 字节。
  SHA-256：`30aedc230f611a0a724f53321b9d52b3b9af2f387fafc7794d4b26a4eaef9311`。
  ZIP 条目位于根层；恢复时对应仓库的 `ISODISTORT/docs/<条目名>`。
- [路由调整前原文](routing-before-compaction-20261009.zip)：包含仓库根 `AGENTS.md`
  和 `ISODISTORT/README.md`，ZIP 条目名分别为 `AGENTS.md`、`README.md`。
  SHA-256：`4f4198622e4e3f8e56346b0d60d313f21501518cfc19b82c4595c2e165291fbb`。

四份说明文档的原始 SHA-256（原文合计 160,980 字节）：

```text
8c8045277b0badc6464034789cac6839bd849939677dd61c7eac0318bd0ff9f3  DEVELOPMENT_PLAN.md
3a40b7c6484826495bb9317378fe4dae7ea25f91aedbbb9a602cace73e4339b6  BUGFIX_VALIDATION_REPORT.md
1cf3682e16286d66eadc8fe2b680f49fa2d00409d1ff3e80cbee5fa11e0715ca  DOWNLOAD_CHECKLIST.md
00115517119c27878a8e19525e914f17cb36ab4543e64448eda2cca7640ccc12  MANUAL_VALIDATION.md
```

已实际解压到独立临时目录，四份文件逐一核对 SHA-256；归档内容同时与编辑前读取的
原文一致。编码、换行和完整历史叙述均可恢复。此验证不涉及科学算法或网页运行。

## 2026-10-16 重置归档

用户要求先完成官网三母相有限全量输出的下载、再重新开始验证与 debug，因此旧的冻结范围与其
批次产物整体作废并归档。重置记录见
[2026-10-16 重置记录（已归档）](reset-records-20261019.zip)。

- [local-web-output-refresh-20261009.zip](local-web-output-refresh-20261009.zip)：
  11,335,017 字节，2,150 个条目。
  SHA-256：`935d62e3ce549a489b3cc773077bfafd2648c0c1616e2c49f365bd3a4ad8dead`。
  条目名即仓库相对路径，覆盖：
  `ISODISTORT/output/validation/web_output_refresh/20261009-local-web-output-refresh-v1/**`、
  冻结清单 `ISODISTORT/docs/manifests/local_web_output_refresh_manifest.json` 的冻结字节
  （SHA-256 `322a75e7d28bd9bd9516309602ff01b298c7faa77e1974b7b562273b2206728d`）、
  以及 `ISODISTORT/output/validation/local_web_output_refresh_manifest_audit_20261009.json`。
- 归档后已实际完整解压到独立临时目录，逐文件核对 SHA-256：2,150/2,150 一致，0 处不符；
  ZIP CRC 自检通过。原批次目录已在核对通过后从工作区移除；冻结清单原地保留为
  `superseded` 记录，其冻结字节只在归档中。
- 恢复方法（条目名含中文，请用 Python 而非 PowerShell `Expand-Archive`，后者会改写编码）：

```python
import zipfile
zipfile.ZipFile(r"ISODISTORT/docs/history/local-web-output-refresh-20261009.zip").extractall(r"<目标目录>")
```

## 2026-10-17：下载清单 v1 原文

下载目标迁移到既有 `output_compare` 树并补正 nmod 覆盖前，v1 官网下载清单的冻结字节归档为
[official-download-manifest-20261016.zip](official-download-manifest-20261016.zip)：9,230 字节，
条目 `ISODISTORT/docs/manifests/official_output_download_manifest.json`。
清单原文 SHA-256：`879c72e98ccdffe982857723b799e69a14b5841cd4734ac47d8778d2cc20a211`；
ZIP SHA-256：`0d5f55a3b2856db3fe2f0d713bfd35c59ffa3b77d69c915e269d4ad9ea846e14`。
创建后已解压到独立临时目录核对条目 SHA-256，CRC 自检通过。该 v1 清单从未执行下载，
只用于追溯范围与计数差异；现行清单见
[official_output_download_manifest.json](../manifests/official_output_download_manifest.json)。

## 查看或恢复

在仓库根目录运行 PowerShell，将文件解压到新建的临时目录查看，不覆盖当前文档：

```powershell
$restoreDir = Join-Path ([IO.Path]::GetTempPath()) ('cris-doc-history-' + [guid]::NewGuid().ToString('N'))
Expand-Archive -LiteralPath 'ISODISTORT/docs/history/docs-before-compaction-20261009.zip' -DestinationPath $restoreDir
Get-ChildItem -LiteralPath $restoreDir -File | Get-FileHash -Algorithm SHA256
$restoreDir
```

与上列摘要逐项比较后，即得到可逐字节还原的原文件。若需恢复某个历史段落，从解压文件
取回；若需回退整个文件，先保存当前版本，再按原路径映射复制。路由 ZIP 用相同方式
解压到另一个新目录，注意其中 README.md 的原位置是 `ISODISTORT/README.md`。
后续压缩创建新的日期/版本 ZIP 并追加索引，不能覆盖本快照。

## 2026-10-17：撤销 DS 范围判据

[ds-scope-before-correction-20261017.zip](ds-scope-before-correction-20261017.zip)
保存修正前的两个 DS 脚本、相关计划/清单/机器记录以及当时存在的网格和闭包报告，共 11 个
原始文件，另含逐文件 `SHA256.json`。ZIP SHA-256：
`1f549199c92a4efb7ae726626f61dc3b3f9d789031e60d94680f72f2215a0b7b`。
创建后逐条解码并重新计算 SHA-256，全部匹配；CRC 自检通过。
按上文 Python 方法解压到独立目录，再对照 SHA256.json；不能覆盖当前文件或恢复旧 checkpoint。
这些内容是**已失效的历史证据**，不是当前计算或下载范围的依据。

## 2026-10-19：全部旧工作记录与本地产物重置

- [旧记录原文](reset-records-20261019.zip)：20 个原始文件，SHA-256 `99a5b7d990cd44359e8b96d62c7e6752e80cfc0424ce12cbd0b6a774ae0c0ce0`。
- [旧本地产物](reset-local-outputs-20261019.zip)：14,754 个原始文件，SHA-256 `6fe92e6d070a8f200bf62633a4c0b7a0fce3eac9723a89d11b75d67ac989f939`。
- 两份 ZIP 均含 `_RESET_SHA256.json`，按仓库相对路径记录原文件名、大小和 SHA-256。
  创建后逐条解压读取、复算全部字节哈希，并通过 ZIP CRC 检查，之后才清理工作文件。
- 归档只作恢复，不是当前检验通过证据。后续从输出 ZIP 恢复了历史官网网页修复副本与上传输入，
  记录见 `output/validation/rebuild_20261019/recovered_references.json`；没有删除主官网原件。
- 旧冻结、checkpoint 和 debug 结论不得用于新批次；输入和理论资料的科学限制仍有效。

恢复到新的空目录（不能覆盖当前文件）：

```python
import zipfile
with zipfile.ZipFile(r"E:\CRIS\ISODISTORT\docs\history\reset-records-20261019.zip") as archive:
    archive.extractall(r"<新的空目录>")
```

按 `_RESET_SHA256.json` 校验后再选择需要的历史文件。产物 ZIP 同理。当前事实入口为
[新下载清单](../DOWNLOAD_CHECKLIST.md)、[新检验任务](../manifests/validation_checklist_20261019.json)。

## 2026-10-19：架构调整前的说明与冻结记录

[architecture-before-20261019.zip](architecture-before-20261019.zip) 保存当时的验证报告、
开发计划、手工操作、历史索引、ISODISTORT README 和 checklist_freeze_20261019.json，
共六份原始文件，另含 `_SHA256.json`。
ZIP SHA-256：`60a08de0f8cabea5aec31754ea7fb113110dd904c35156a28c6364c22e94f3ca`。
创建后逐条读取解压字节核对 SHA-256，全部一致；ZIP CRC 检查通过。
条目名即仓库相对路径，按上文 Python 方法解压到新的空目录，再按 `_SHA256.json` 复核。
旧冻结记录可恢复；架构调整只修订支撑文档/工具指纹，下载任务及原件身份继续使用现行冻结。
旧报告的全量测试中断/失败记录保留，不改写成当前通过证据。

## 2026-10-19：测试整合前的源码、说明和冻结记录

[test-layout-before-20261019.zip](test-layout-before-20261019.zip) 保存测试/审计源码、当前说明、
能力映射、规则和冻结记录等 145 个文件，另含 `_SHA256.json`。
ZIP SHA-256：`495c76d18c76294a3452554b2ee5ab3e042382213802335c074ff412cfe5e1d6`。
创建时逐文件读取解压字节、比对 SHA-256，并通过 ZIP CRC 检查；条目名就是原仓库相对路径。
恢复时按本页方法解压到新的空目录并检查 `_SHA256.json`，不得覆盖当前源码或恢复旧 checkpoint。
历史报告、code_context.json 和既有复用登记中的 generator 路径保留当时含义；当前脚本入口
以 `tests/README.md`、能力映射和新版冻结记录为准。位置迁移不会使历史科学结果获得新验证。

## 2026-10-19：安装入口和启动器目录修复前快照

[launcher-before-20261019.zip](launcher-before-20261019.zip) 保存 14 个改动前的源码、配置、
说明和第 3 版冻结记录，另含 `_SHA256.json`。ZIP SHA-256：
`5342bfbe223e32e976e708ba8e807d3b08e89cdfa5fa48b364b4fec21c7c194b`。
已通过 CRC 和逐项解压字节哈希核查；按原相对路径解压到新的空目录恢复，不覆盖当前文件。

`reset-local-outputs-20261019.zip` 是本机的大型派生产物恢复包，由 Git 忽略并保留在原位置；
远程源码仓库不包含该包、官网原件、实验输入或本机启动器。下载源码不能代替备份这些本地材料。
