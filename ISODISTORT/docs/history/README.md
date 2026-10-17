# 说明性文档可恢复归档

当前待办见 [开发计划](../DEVELOPMENT_PLAN.md)。这里保存精简前工作树的完整原始字节，
包括当时未提交的文档修改；压缩文件是历史快照，不作为当前执行计划。
机器 manifest、原始官网资料与机器验证产物未被本次压缩改写。

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
[validation_reset_20261016.json](../manifests/validation_reset_20261016.json)。

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
