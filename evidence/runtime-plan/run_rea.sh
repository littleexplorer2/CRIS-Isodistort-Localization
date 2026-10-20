#!/bin/sh
set -eu
task_runtime=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
# Keep analysis artifacts here. REA's supported long-path socket allocator
# places only the Unix IPC socket in a private /tmp directory (>107 bytes),
# avoiding the lack of Unix socket support on the Windows DrvFS mount.
task_session_temp="$task_runtime/temp/ghidra-projects-and-analysis-artifacts"
mkdir -p "$task_session_temp" "$task_runtime/cache" "$task_runtime/config"
export GHIDRA_INSTALL_DIR="$task_runtime/unpacked/ghidra_12.1.3_PUBLIC"
export JAVA_HOME="$task_runtime/unpacked/jdk-21.0.12.1+1"
export TMPDIR="$task_session_temp"
export XDG_CACHE_HOME="$task_runtime/cache"
export XDG_CONFIG_HOME="$task_runtime/config"
export JDK_JAVA_OPTIONS="-Duser.home=$task_runtime/cache -Djava.io.tmpdir=$task_session_temp"
export REA_ANALYSIS_PROVIDER=ghidra
exec "$task_runtime/unpacked/node-v24.15.0-linux-x64/bin/node" \
  "$task_runtime/rea/node_modules/rea-agents/scripts/rea.mjs" "$@"
