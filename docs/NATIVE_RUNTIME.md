# 本项目的 Linux 原生逆向工具环境

本项目采用可移走的工具目录 `tools/runtime/`；该目录已加入 `.gitignore`。
没有全局安装 Ghidra/JDK/Node，没有修改 Codex 的 MCP 注册，也没有修改 `HOME`。
这一环境仅用于静态分析，不是 CRIS 的生产计算依赖。

## 固定版本与下载身份

元数据于 2026-10-20 从官方来源查询；证据在 `evidence/runtime-plan/`。
共下载 808,082,961 字节的三个压缩包。

- **Ghidra 12.1.3**：
  [NSA 官方 release](https://github.com/NationalSecurityAgency/ghidra/releases/tag/Ghidra_12.1.3_build)，
  文件 `ghidra_12.1.3_PUBLIC_20260817.zip`，569,445,154 字节。
  SHA-256：`93a5d11a9ad510622acaaf908c556a7b9b764d338e78a7567f3689bf5081fd54`。
  对应 tag 的 `Ghidra/application.properties` 明确要求 Java 最低 21，编译目标也是 21；
  没有用当前 `master` 的不同 Java 要求替代这个固定 release 的要求。
- **Eclipse Temurin JDK 21.0.12.1+1，Linux x64 HotSpot**：
  [官方发布包](https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.12.1%2B1/OpenJDK21U-jdk_x64_linux_hotspot_21.0.12.1_1.tar.gz)，
  207,473,347 字节。
  SHA-256：`ce79869e1307ed8ee1e2baa86a412b1eb5b75d10a01006d788a6f968bcfaee94`。
- **Node.js 24.15.0，Linux x64**：
  [官方发布包](https://nodejs.org/dist/v24.15.0/node-v24.15.0-linux-x64.tar.xz)，31,164,460 字节。
  SHA-256：`472655581fb851559730c48763e0c9d3bc25975c59d518003fc0849d3e4ba0f6`，
  来源为同目录官方 `SHASUMS256.txt`；字节数通过该 URL 的 HTTP HEAD 确认。
- **REA npm 包 6.7.0**：
  [npm registry 元数据](https://registry.npmjs.org/rea-agents/6.7.0)声明公共 bin 为
  `scripts/rea.mjs`，Node 支持范围包含 24.15.0。
  包完整性值为
  `sha512-yhIx2ioJYYT1C/tONfYeitC6LhnKmCh+ZwRUOgYuhSc733lGhkkt73zhug9WmLgxoMdKuqIYGCoBfD4fXkOhFw==`。
  项目内使用固定版本安装；`package-lock.json` 位于 ignored 的 `tools/runtime/rea/`。

下载在 Windows 侧完成，以避免 WSL NAT 无法使用 Windows localhost proxy 的问题。
下载器逐包核对字节数及 SHA-256 后才把 `.part` 改为正式包名。
解压前确认所有归档成员均是相对路径、没有 `..` 逃逸；只向新的 `tools/runtime/unpacked`
目录写入。归档工具的版本与安装结果以 `evidence/runtime-plan/installation.json` 为准。
运行期版本核查和 REA 的实际工具握手是两个不同的验证层级。

## 目录与公共入口

```text
tools/runtime/
  archives/                         # 三个校验后的官方包
  unpacked/
    ghidra_12.1.3_PUBLIC/
    jdk-21.0.12.1+1/
    node-v24.15.0-linux-x64/
  rea/node_modules/rea-agents/
  npm-cache/
  temp/
  cache/
  config/
  run_rea.sh
```

公共 REA MCP 入口为：

```text
<runtime>/unpacked/node-v24.15.0-linux-x64/bin/node
<runtime>/rea/node_modules/rea-agents/scripts/rea.mjs
mcp
```

`package.json` 的 bin 声明和 `scripts/rea.mjs` 都已审读。该脚本在 `mcp` 模式直接加载
发行包的 `dist/mcp/main.js`，没有猜测内部执行入口，也没有使用 REA 安装/注册命令。

项目包装器 `tools/runtime/run_rea.sh` 接受原样的 REA 参数，进程级设置：

- `GHIDRA_INSTALL_DIR=<runtime>/unpacked/ghidra_12.1.3_PUBLIC`
- `JAVA_HOME=<runtime>/unpacked/jdk-21.0.12.1+1`
- `REA_ANALYSIS_PROVIDER=ghidra`
- `TMPDIR=<runtime>/temp/ghidra-projects-and-analysis-artifacts`
- `XDG_CACHE_HOME=<runtime>/cache`、`XDG_CONFIG_HOME=<runtime>/config`
- `JDK_JAVA_OPTIONS` 设置项目内 Java 用户目录和同一 TMPDIR

REA 包中 `dist/config/parseConfig.js` 确认前三个 provider/JDK 配置名。
`dist/ghidra/GhidraSessionRoot.js` 说明临时会话路径不能含以点开头的路径段；
本项目选取的 `tools/runtime/temp` 满足这一约束。
`dist/ghidra/GhidraJavaLaunch.js` 还会把具体会话的 Java home/temp 属性放到 JVM 参数里，
将会话产物进一步隔离。

Windows 挂载盘不支持本次 Ghidra 所需的 Unix-domain socket。首次分析已保留完整失败
记录；启动器随后采用 REA 自带的长路径回退：项目内保留分析文件，私有 Linux `/tmp`
目录仅放短命通信 socket，结束时清理，不存储研究结果。
精确根因、源码定位和实际绑定探针见 [GHIDRA_TRANSPORT_RECOVERY](GHIDRA_TRANSPORT_RECOVERY.md)。

WSL 调用示例（目录内执行不需要注册 MCP）：

```powershell
wsl.exe --exec sh /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/runtime/run_rea.sh mcp
```

`mcp` 使用标准输入/输出交换 JSON-RPC，不能用“进程成功启动”代替 `initialize`、
`tools/list`、`tools/call` 的真实握手。项目的 `tools/rea_mcp_batch.py` 保存这类证据。

## 依赖选择与限制

准备 REA 时在 Windows 侧运行了：

```powershell
$env:npm_config_cache = '<CRIS_REA>\tools\runtime\npm-cache'
$env:TEMP = '<CRIS_REA>\tools\runtime\temp'
$env:TMP = $env:TEMP
npm install --prefix '<CRIS_REA>\tools\runtime\rea' `
  --ignore-scripts --omit=optional --no-audit --no-fund --save-exact rea-agents@6.7.0
```

这会跳过安装脚本和可选依赖；因此它**不是所有 REA 功能的完整安装**。
本任务使用 Ghidra 静态原生分析，Frida、交互 PTY、浏览器等其它能力不在此安装的验收范围。
是否能在 Linux 上导入公共入口，以及 Ghidra provider 能否处理所选 ELF，
必须由实际 Linux 运行与 MCP 证据回答。

GitHub 重定向下载可能包含短时签名 query。本地 ignored 的下载日志可保留用于诊断，
提交的证据只记录固定官方 URL、版本、大小、哈希和结果，不提交临时签名 URL。

本文件是工具安装与可复现说明，不构成科学计算正确性证明。

## 本次 REA 实际结果

成功与失败的 MCP 逐条响应分别位于
`evidence/native/rea-linux-iso-recovered/` 和早期 `rea-linux-iso/`。
恢复后取得 overview 和十个目标函数的伪代码，`MAIN__` 请求仍失败。
随后导出 12 条 Evidence 并成功关闭目标、保存分析快照。
具体算法结论和限制只在 [NATIVE_ANALYSIS](NATIVE_ANALYSIS.md) 维护。

安装目录不进入 Git，但本次实际执行的启动器和独立 socket 探针副本保存在
`evidence/runtime-plan/run_rea.sh` 与 `probe_ghidra_socket.mjs`；
重放时把它们放回 `tools/runtime/`，不要直接在证据目录运行。
先恢复上述固定工具版本及目录，再运行：

```powershell
wsl.exe --exec python3 /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/rea_mcp_batch.py `
  --plan /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/rea_iso_plan.json `
  --output /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/evidence/native/rea-replay-new `
  -- sh /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/runtime/run_rea.sh mcp
```

重放前须将计划中的目标和导出路径改成该机器的实际项目路径，并给 bundle/snapshot
选择新的文件名；不能覆盖既有证据。计划保留了本次失败的 `MAIN__` 请求，用于准确
复述本次调用范围；该请求成功与否必须单独检查 `isError`，不能只看脚本退出码。
