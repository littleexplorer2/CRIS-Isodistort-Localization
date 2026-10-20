# REA / Ghidra 通信故障与隔离恢复

## 已观察到的故障

`evidence/native/rea-linux-iso/04-overview.json` 保存了原失败证据：Ghidra 对 ISO 二进制
完成静态分析后，`ReaGhidraBridge.serveUnixSocket` 在
`UnixDomainSockets.bind` 抛出 `Operation not supported`。这次没有成功建立分析桥，
不能把静态分析日志中的 “Analysis succeeded” 等同于 REA 查询成功。

原临时目录位于 Windows 文件系统的 WSL `/mnt/c/` 挂载中。现象与该文件系统不支持
所需 Unix socket 的限制一致；本次只修复通信地址的放置，不改变目标二进制或分析器。

## 实际审读范围

本项目隔离安装的 `rea-agents` 为 **6.7.0**。审读了以下发行包实现：

- `dist/ghidra/GhidraProviderClient.js:216`：Windows 使用
  `authenticated-loopback-tcp`，其他平台固定使用 `unix-socket`。
- `dist/config/environment.js:14` 起：未发现可选择 Ghidra transport 的公开环境项。
- `bridge/ghidra/ReaGhidraBridge.java:251` 起：底层桥同时实现 Unix socket 和
  经认证的回环 TCP，但不意味着 Linux provider 暴露了 TCP 配置入口。
- `dist/ghidra/GhidraTransport.js:15` 起：若 Linux Unix socket 路径长于 107 字节，
  **只为 socket** 创建 `/tmp/rea-ghidra-socket-*/s`，分析 runtimeRoot 不迁移。
- `dist/ghidra/GhidraClient.js:514` 起及 `:779` 起：分析目录与 socket 目录分开管理；
  关闭时分别清理，清理失败会保留失败状态。
- `dist/process/PrivateRuntimeRoot.js`：POSIX 私有目录通过 `mkdtemp` 创建。

因此不修改发行包、不伪装 Windows 平台，也不声称存在未发现的 TCP 开关。

## 本次隔离配置修改

仅修改未跟踪的 `tools/runtime/run_rea.sh`：

```sh
task_session_temp="$task_runtime/temp/ghidra-projects-and-analysis-artifacts"
mkdir -p "$task_session_temp" "$task_runtime/cache" "$task_runtime/config"
export TMPDIR="$task_session_temp"
export JDK_JAVA_OPTIONS="-Duser.home=$task_runtime/cache -Djava.io.tmpdir=$task_session_temp"
```

较长但有明确含义的项目内目录使 REA 既有长路径处理逻辑生效。本机候选 socket 路径
为 132 字节，超过 107 字节限制；REA 将其短暂 IPC socket 放入 Linux 原生 `/tmp`，
分析快照、项目、session 描述和日志仍使用桌面 `CRIS_REA` 内的 runtimeRoot。
没有更改 `HOME` 环境变量；原来的 Java `user.home` 配置保持原值。

`/tmp` 中的短命 IPC socket 不包含逆向分析结果，完成通信后由 REA 清理。
没有移动整个项目或官方数据到 Linux `/tmp`。

## 已执行的定向验证

未启动新的 Ghidra 分析或 MCP 批次，仅运行一个独立通信探针：

```powershell
wsl.exe --exec /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/runtime/unpacked/node-v24.15.0-linux-x64/bin/node /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/runtime/probe_ghidra_socket.mjs
wsl.exe --exec sh -n /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/runtime/run_rea.sh
```

两条命令退出码均为 0。真实记录：
`evidence/environment/ghidra-socket-transport-probe.json`。

探针调用发行包原有 `createGhidraEndpoint`，确认选出的地址为
`/tmp/rea-ghidra-socket-EZEiwY/s`、父目录权限 `0700`、socket 绑定成功且清理后目录
不存在。未复制/执行待逆向程序。脚本语法检查通过。

**证据范围：** 此处证明 socket 地址选择、绑定和清理可行，尚不等同于随后整个
Ghidra 分析批次成功。完整批次必须由主代理重新启动并保存其真实结果；失败批次证据
继续保留，不覆盖为成功。
