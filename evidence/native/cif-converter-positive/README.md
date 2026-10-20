# 合成 CIF 转换正例

本例实际通过冻结的 `cris_rea` 应用运行了 `findsym_cifinput`，再把转换器
**未编辑的原始 stdout 字节**交给 `findsym`。两次程序均退出 0，无包装层已识别的错误诊断。
人工审读确认转换器输出为 FINDSYM 关键词格式；FINDSYM 接受输入，输出 P1、
两个原子及完整结果 CIF。

输入是本次独立创建的合成格式测试，不是实验结构或官网输出。三斜晶格参数为
`4.123, 5.287, 6.431 Å; 78.2, 83.7, 74.6°`，明确声明 P1 和恒等操作 `x,y,z`，
并指定 Si、O 两个不同分数坐标位点。
[input_manifest.json](input_manifest.json) 在执行前固定了输入、12 个应用源码文件和官方程序指纹。
第二阶段的输入来源与哈希在 [findsym_input_manifest.json](findsym_input_manifest.json) 中单独固定。

[verification.json](verification.json) 记录：

- CIF 输入哈希与转换器实际输入一致，下游输入哈希与转换器 stdout 完全一致。
- 转换器保留了晶格、原子种类/数量、分数坐标、占据率、恒等操作和所声明的群。
  十进制参数用 `Decimal` 精确相等比较，没有放宽容差。
- 运行后应用源码与预先冻结的 12 个文件一致。
- FINDSYM stdout 明确报告 P1 并产生完整 CIF；这是此输入的格式接受证据，
  不把输入中已有的 P1/操作再次输出称为独立的群理论证明。

程序自己的容差由转换器输出：晶格 `1e-5 Å`、原子位置 `1e-3 Å`、
位置最大容差 `0.33`（这里只原样记录，未独立解码该额外字段的完整作用）、
占据率 `1e-3`。包装层没有修改这些字段。
stdout 主标题报告 FINDSYM 7.1.6，内嵌 CIF 注释报告 writer 7.1.3；
两处原文均保留，工具身份以二进制 SHA-256 为准。

两次完整运行分别保存在 ignored 的
`runs/479ce6fd-3d63-410b-bedc-68c9d658860e/` 和
`runs/0538df13-283f-4680-89b2-5ad58feb1c96/`。
本目录只保存输入及较小的 stdout/stderr、执行 JSON 和 evidence JSON，没有复制数据库。
没有读取官网案例结果，没有修改生产代码或原 CRIS 目录。

本例验收范围为单个明确合成输入的真实转换与格式接受，
不覆盖任意 CIF 语法、磁性结构、不同 setting 或完整 ISODISTORT 科学验收。
