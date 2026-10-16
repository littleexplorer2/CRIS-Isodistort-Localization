# 当前验证缺口与证据入口

待办与调度只维护在 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)；
下载动作见 [DOWNLOAD_CHECKLIST.md](DOWNLOAD_CHECKLIST.md)。
已完成 debug、下载、迁移和修改历史的完整原文见 [可恢复归档](history/README.md)。

## 集成验证与管线合同检查（2026-10-16）

验证公共入口为 `tests.validation`；生产模块不导入测试支持包。入口、默认配置、
源码签名、部署器和架构清单已按当前目录统一。原 CIF 比较器 AST 与迁移前相同；
Method 4 的两个数据合同及五个计算方法与迁移前快照 AST 相同，这只证明移动未改写
相应算法，不代表对全部科学输入的认证。公共兼容入口继续调用同一实现。

本轮涉及科学语义的共享合同如下：

- 心化矩阵采用精确有理数；惯用胞中原胞个数为 `1 / abs(det(P))`，原胞/惯用胞
  体积比与该指数的乘积为一。验证覆盖 P/A/B/C/I/F/R 的独立晶格重数；不引入
  材料相关比例或振幅拟合。
- 超胞采用行约定 `x_parent = x_child @ B`；物种和占据保持一致，位点及整数平移
  必须在给定笛卡尔误差球中唯一。容差单位为 Å，默认值来自配置
  `symmetry_cartesian_tolerance_angstrom`。重排和扩胞坐标守恒为正例；缺失站点、
  斜晶格和同站点多个周期像为反例。该合同仅适用于未畸变对应，不代替畸变结构匹配。
- 数学依据为资料库补充证明 `sec:shared-mapping`，对应 C02/C03；Method 4 归属
  与 C15 的来源链同步。这里未改变动力学或能量模型，能量守恒不是这次重构的验收项。
- 比较容差沿用测试 YAML 的晶格、周期分数坐标和标量三个域；非有限或负容差拒绝。
  语义、声明的参考哈希与严格字节条件逐项成立才通过；空批次不表示成功。
  比较器只读取冻结文件，验证参数不回传生产计算。

原始证据位于 `output/validation/pipeline_cleanup/`（不提交）：
`final-source-evidence.json` 保存源码/配置/依赖指纹和 AST 核对；`repo-final.xml`
记录仓库门禁及实际 wheel 检查 34 passed；目录清理后的 `directory-removal.xml`
记录门禁、集成比较和受影响配对测试 75 passed、1 skipped（缺少外部 CIF 配对）。
完整科学回归于 2026-10-17 完成：**1132 passed、3 skipped、0 failed**，耗时
20 分 18 秒，见同目录 `full-final.xml` 与 `full-final.log`。三个跳过项分别缺少
X4- P3 外部配对文件、需要原生 POSIX shell、缺少 EuAl4 LD1 `g=1/6` 本地子群数据库。
运行期间的目录清理只改动测试支持和路径引用，生产三层指纹未变；对应最终版本已由
上述 75 passed、1 skipped 的补充复测覆盖。源码指纹、Ruff F821/F401 和
`git diff --check` 终检通过。

资料完整性及引用检查通过，两篇 LaTeX 编译成功；新增证明页面的渲染检查范围见
资料库 `audit/DELIVERY_REVIEW.md`。未运行官网全量差分或人工 GUI 验收，既有科学
缺口仍保留；静态去重与调用扫描不构成对所有动态路径不存在死代码的形式化证明。

## 理论与代码一致性审查（2026-10-15）

本轮按文献和两份证明检查实际调用链，详细反例、推导、修复和限制统一见
[理论审查报告](reference_library/audit/THEORY_CODE_AUDIT.md)。
[实现清单](reference_library/audit/implementation_audit.json) 保存逐能力状态与最终检查数字；
原始日志、JUnit、k 参数回代及源码指纹在
`output/validation/theory_design_audit_20261015/`（相对于 ISODISTORT）。

已移除官网存档回填、占据率排序二分/夹取、以位移代替转动和按重数猜 Wyckoff 字母；
修复行基矢点群作用、精确晶格判定、分数显示和含常数/负号/耦合的波矢换算。
Wyckoff 标签增加实际 Hall/原点门禁，标准表种子经独立交叉检查；旧母胞映射的均匀
复制与最近邻补位已替换为既有完整子胞映射的覆盖性检查。
最小二乘仍用于有来源且满秩的固定模式基分解，它不是材料答案拟合。
未实现物理模型明确拒绝；未对全部二进制、表数据、晶体、官网版本或 GUI 作认证。
此前三层重构的测试数字只适用于其当时源码，不代替本轮结果。

## 三层重构验证（2026-10-15）

三个子项目统一为 `frontend → application → backend`，模块目录与维护约束以根
`AGENTS.md`、各子项目 `agent.md` 及 `tests_repo/architecture_manifest.json` 为准。
本轮仅修改 CRIS；选课工具只读参考。原有工作区改动保留。

- 科学算法与会话迁移：对照本轮修改前快照，31 个文件在规范导入替换后的 AST 相同。
  该检查证明迁移没有改写这些计算实现，不代表新增科学输入已验证。
- 首次 ISODISTORT 完整运行：1034 passed、4 failed、3 skipped。两项失败是本轮前已存在
  的导出状态报告断言遗漏（在修改前快照复现）；另两项是旧会话包与新包的 monkeypatch
  状态未共享。分别补齐报告内容/文件保持断言、改为同一模块对象兼容别名后，相关文件
  复测 53 passed/1 skipped 与 51 passed。最终完整复核为 **1040 passed、3 skipped、
  0 failed**，耗时 19 分 24 秒；证据为同目录 `full-final-results.xml` 与
  `full-final-pytest.log`。
  三个跳过项分别缺少 X4- P3 官网/本地配对文件、需要原生 POSIX shell，以及缺少
  EuAl4 LD1 `g=1/6` 本地子群数据库；本轮不生成可能耗时数小时的子群数据库。
- 仓库级门禁与实际 ISODISTORT wheel：32 passed；当时的比较工具：29 passed；
  ISOVIZ_INPUT：40 passed。当时的两个小工具均包含隔离源码后的实际 wheel/CLI 检查；GUI 未启动。
- 三个 wheel 均检查规范入口和资源；两个小工具的 YAML 由唯一源配置在构建时打包，
  自定义配置与包内配置分别保持正确的相对路径基准。旧入口的可变状态兼容有回归覆盖。
- 部署诊断：24 pass、2 warn、0 fail；两个 warn 分别为缺少成对比较 CIF 与振幅 CSV，
  没有将无用户输入视作真实使用验收成功。静态 Ruff F 与 `git diff --check` 通过。
- 手工审计源码签名已覆盖实际 `backend/`、`application/` 及独立比较器的实现，避免
  迁移后只签名旧导入壳。原历史报告继续只适用于当时的源码签名。
- 提交前补查手工网页下载器，将应用层纳入运行前后源码摘要，并新增 2026-10-15
  运行修订；旧冻结输入和旧修订未改写。回归以独立临时文件验证三层源码变化会改变
  摘要、字节码不会改变摘要、合同和父 manifest 漂移会被拒绝；没有启动浏览器下载。

本轮原始日志、JUnit XML、部署诊断及 AST 对照：
`../output/validation/layer_refactor_20261015/`（不提交）。打包测试所需的 setuptools/wheel
临时安装在 CRIS 内的 `ISOVIZ_INPUT/.test-tmp/build-deps`，共享 `.venv` 没有重建或改写；
长期开发依赖已在各项目开发依赖清单中登记。

未执行官网全量重新下载、全部科学输入扩展矩阵或 IsoVIZ/VESTA 人工 GUI 验收；
下面已有科学验证缺口继续保留。

## 当前需要核验

- 本地输出刷新：执行冻结清单的 74 项真实网页运行；冻结自检通过不等于网页或科学结果通过。
- NdNiO2 Y `a=1/3`：完成 nmod=0/1 全候选长时验收，逐项审查导出、跳过原因和网页状态。
- 官网补证：EuAl4 G05 auto-origin basis HTML；无法取得仍保留证据警告。
- 4310 Method 2 扩展及 Method 3/4：先定义输入，再执行；未定义范围不能称为漏下载。
- 泛用性、非 Γ 来源合同、混合占位、跨晶系、工具 GUI、性能与并发缺口：
  按开发计划逐项验收。

## 正在使用的证据入口

- 冻结输入：`manifests/local_web_output_refresh_manifest.json` 及其引用的 Method 3/4 manifest。
- 冻结自检：`../output/validation/local_web_output_refresh_manifest_audit_20261009.json`。
- 导出诊断：`../output/validation/ndnio2_export_audit_20261006/`。
- 新批次产物、报告和 checkpoint：`../output/validation/web_output_refresh/`。
- 既有机器报告仍在 `../output/validation/`；历史索引和适用源码签名可从归档恢复。

解释任何结果前核对输入、实际源码、配置、依赖、后端及结束签名；
签名漂移、工具缺失、GUI 未运行或证据不足须显式记录。
归档中的通过结论仅适用于其原有签名和范围。

## 主要科学依据

统一入口为[理论参考资料库](reference_library/README.md)，其中维护作者论文、独立方法、
IUCr 标准与教程、来源审核、能力映射及两份带引用的 LaTeX 论证。
书目元数据只在 `reference_library/sources.json` 维护；当前编译与资料完整性检查见
`reference_library/audit/validation.json`，数学/实现的适用边界见
`reference_library/audit/QUALITY_REVIEW.md`。资料整理未重跑本报告中的科学或 GUI 批次。

这些文献定义算法与不变量；官网结果只用于检验实现是否复现同一科学语义。
