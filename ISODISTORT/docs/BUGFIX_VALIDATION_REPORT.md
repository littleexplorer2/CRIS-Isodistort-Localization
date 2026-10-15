# 当前验证缺口与证据入口

待办与调度只维护在 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)；
下载动作见 [DOWNLOAD_CHECKLIST.md](DOWNLOAD_CHECKLIST.md)。
已完成 debug、下载、迁移和修改历史的完整原文见 [可恢复归档](history/README.md)。

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
- 仓库级门禁与实际 ISODISTORT wheel：32 passed；ISODISTORT_VALIDATE：29 passed；
  ISOVIZ_INPUT：40 passed。两个小工具均包含隔离源码后的实际 wheel/CLI 检查；GUI 未启动。
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

- [ISODISTORT Method 3 help](https://iso.byu.edu/isodistorthelp.php)
- [ISOTROPY user documentation](https://iso.byu.edu/isotropy_doc.php)
- Campbell, Stokes, Tanner & Hatch, *J. Appl. Cryst.* 39 (2006), [DOI 10.1107/S0021889806014075](https://doi.org/10.1107/S0021889806014075)
- Stokes & Campbell, *Acta Cryst.* A73 (2017), Appendix B, DOI `10.1107/S2053273316017629`
- Hatch & Stokes, *Phys. Rev. B* 65 (2002), DOI `10.1103/PhysRevB.65.014113`
- Stokes, van Orden & Campbell, *J. Appl. Cryst.* 49 (2016), DOI `10.1107/S160057671601311X`
- Wagner & Schönleber, *Acta Cryst.* B65 (2009), [superspace introduction](https://journals.iucr.org/b/issues/2009/03/00/bk5084/index.html)：有理调制、公度周期、cosine/sine 分量及相位约定。
- [IUCr CIF specifications](https://www.iucr.org/resources/cif/spec)：data block、loop 与结构数据的语法边界。
- [International Tables, Vol. C, §1.1](https://it.iucr.org/Cb/ch1o1v0001/)：原/惯用胞、中心化平移与倒格矢的定义及体积关系。
- [IUCr Superstructure definition](https://dictionary.iucr.org/Superstructure)：子胞平移群必须为母相平移群的子群。

这些文献定义算法与不变量；官网结果只用于检验实现是否复现同一科学语义。
