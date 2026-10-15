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
