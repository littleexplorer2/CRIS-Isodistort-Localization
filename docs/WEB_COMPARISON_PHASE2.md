# 同一独立母相的官网对照

本记录是一个有限样例的独立比较。结论为：同一自构母相、GM4− / P1 选择下，
官网与本地得到相同子群和位移空间。没有验证全部 k、全部空间群、逐列 IR 基
或当前官网应变约定等价。

## 输入冻结与隔离

母相是新项目自行构造的立方 BaTiO₃ 几何，a=4 Å，不是用户实验数据。
原始文件 802 字节，SHA-256：

    8bd7d4afaefcde08c6666a1467f24da4aa25659709bc3e38a8fbc3eff33e602a

本地已完成结果位于
`evidence/phase2/workflow/20261020T145121147627Z/cubic_polar_axis/`。
2026-10-20 14:56:02 UTC，先对该目录 15 个文件逐一记录字节数与 SHA-256，
再开始官网请求。冻结记录为
`evidence/phase2/web-comparison/20261020_cubic_polar_axis/local_freeze_manifest.json`。
每次官网阶段结束及最终比较前均重验这些指纹，没有重新计算或改写本地结果。

本次没有用旧 SrTiO₃ 页面代替新上传结果。新官网响应只供比较，保存于上述
`web-comparison` 子目录；产品模块未导入捕获或比较脚本。
冻结结果早于随后的 CIF 特殊位置导出修复，因此本记录也不代替修复后导出验证。

## 实际官网交互

通过官网返回的真实表单字段，以 HTTPS 请求执行：

1. GET `isodistort.php` 取得上传表单。
2. POST `isodistortuploadfile.php`，`multipart/form-data`，
   `input=uploadparentcif`，`toProcess` 携带冻结的原始 CIF 字节。
3. 官网返回带服务器临时文件引用的自动提交中间表单；原样 POST
   `isodistortform.php` 取得母相搜索页。
4. Method 2：选择显示为 `GM, k12 (0,0,0)` 的实际选项，随后选择 GM4−，
   选择 `P1 (a,0,0)`，进入畸变页。
5. `input=displaydistort`、`origintype=completemodesdetails` 取得完整模式页。

七次响应均为 HTTP 200；完整模式响应开始于 14:58:06 UTC，4303 字节。
目录保留请求体、字段字典、原始 HTML、表单解析结果及请求/响应元数据。
上传中的临时路径是官网返回值，不是本地访问的路径，不能证明服务器调用了
哪个进程。此新流程采用 Python HTTPS 表单重放，**没有冒称为新的 REA 浏览器
捕获**；先前 REA 观察的协议依据见 `WEB_PROTOCOL.md`。

阶段脚本最初没有处理自动提交中间页，并把 GM 显示文本误当成带星号的选项值；
两次本地选择检查因此拒绝继续，未发送错误 k 请求。随后根据实际表单修正选择。
`capture_script_history.json` 及三份脚本快照保留可重建并匹配原先 SHA 的版本，
没有把失败隐藏为官网成功。

## 精确子群比较与位移张成空间

官网报告：

    99 P4mm
    basis = {(0,1,0), (0,0,1), (1,0,0)}
    origin = (0,0,0)
    primitive size = 1
    group index = 6

将打印基矢解释为列后，这些字段与冻结本地结果逐项精确一致。
官网有四个位移列，分别按 A、B、O 站点表示分类；本地同为四列，按三个内部
位移方向及一个允许的均匀平移方向排列。二者列名和列向量不相同，比较对象是
它们张成的空间。

独立比较器 `tools/compare_web_phase2.py` 只读取冻结输入、冻结 JSON 和新官网
HTML，不导入产品计算模块。它从原始输入晶胞参数重建列晶格矩阵 A，从官网
读出列基 C，令 B=A C；按精确的超胞位置匹配原子，并将官网每个分数位移向量
转换成 `B v`。对转换后的 15×4 矩阵作 QR，比较两个正交投影：

    max_abs(Q_web Q_web^T − U_local U_local^T)
      = 3.3306690738754696e-16

交叉基矩阵的四个奇异值在双精度舍入范围内均为 1。本例检查界为
`2048 eps max(15,4)=6.821210263296962e-12`。这一数值结论针对官网**展示出来**
的模式行；它不测量服务器内部未打印系数的误差。所展示位移分量在本例为 0 或 1，
没有用官网列调整本地列或修改本地归一化。

独立负对照也能被区分：删除一个官网方向后投影差为约 0.5，交换父笛卡尔 x/y
方向后差为约 1，均远超检查界。机器记录为 `comparison.json` 和
`comparison_controls.json`。

## 官网 normfactor 与未覆盖范围

官网前三列的分数向量经 B 转换后范数为 4 Å，打印 normfactor 为 0.25000 Å⁻¹；
第四列范数为 `4 sqrt(2)` Å，按定义的倒数为 0.17677669529663687 Å⁻¹，
官网打印 0.17678。后者差 `3.3047033631e-6`，位于五位小数的半单位
`5e-6` 内。直接使用打印因子会得到范数约 1.0000186942，不能把显示舍入
误差当成算法不正交，也不能用该打印因子覆盖本地结果。

本例子群为原胞，常规胞与原胞 RSS 的区别恰好不显现；它不能验证心化子群
下的官方 As/Ap 因子。官网还显示两个应变模式，但此次只归档，未用它们来
证明现代官网有限应变或归一化合同。没有完成网站特有模式 CIF、畴、磁性、
占位、调制或不公度能力的等价验收。

## CIF 导出修复的独立复核

另一个不依赖官网的 P3 特殊位置反例，曾使旧 H-CIF 的精确读回从重数 1
变成 3。原失败完整保留于
`evidence/phase2/theory/review_export_roundtrip.json`。
主代理增加读回校验及必要时的 P1 全胞表示后，独立调用修后 `checked_cif`
得到 `P1_full_cell`，读回一个原子，分数坐标残差
`3.3306690738754696e-16`，晶格度量残差 0。新证据为
`review_export_roundtrip_fixed.json`；没有覆盖旧失败。
P1 在这里仅是避免有限小数特殊位置约束歧义的序列化表示，目标 H 仍在 JSON
和独立元数据中；不能称物理对称性已被判定为 P1。

实际命令分别为：

    python tools/capture_web_phase2_comparison.py freeze
    python tools/capture_web_phase2_comparison.py upload
    python tools/capture_web_phase2_comparison.py process
    python tools/capture_web_phase2_comparison.py k
    python tools/capture_web_phase2_comparison.py irrep
    python tools/capture_web_phase2_comparison.py opd
    python tools/capture_web_phase2_comparison.py complete
    python tools/compare_web_phase2.py
    python evidence/phase2/theory/review_export_roundtrip.py evidence/phase2/theory/review_export_roundtrip_fixed.json --checked

捕获脚本和比较脚本拒绝覆盖既有阶段产物；上述命令记录实际执行过程，不应
直接重复运行来覆盖本轮证据。没有运行原 CRIS 全量回归。
