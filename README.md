# `LD1_C1_gradient_descent_tianren.ipynb` 说明

本文件说明**原始参考笔记本** `LD1_C1_gradient_descent_tianren.ipynb` 做什么、怎么跑、内部怎样把衍射强度拟合成模式振幅。

- 本笔记本是只读对照：日常改代码请用同目录的 `GD_modified.ipynb`。
- 标题虽写 “EuAl4 version Fit modes”，文件名里的 **LD1 / C1** 才是这次拟合的不可约表示与子群。

环境与修改边界见同目录 [`agent.md`](agent.md)。CSV 自动导出只存在于 `GD_modified.ipynb`，**本笔记本不会写出桌面 CSV**。

---

## 1. 它在解决什么问题

实验上有一份 EuAl4 相关的衍射强度表（HKL + `intensity_exp`）。晶体学上，母相是四方 **I4/mmm**，晶胞大约 \(a=4.402\,\text{Å}\)、\(c=11.163\,\text{Å}\)。畸变用 ISODISTORT 的一组 **位移模式** 描述：每个模式有一个振幅，振幅改变原子坐标，从而改变结构因子 \(F(hkl)\) 和强度 \(I \propto |F|^2\)。

本笔记本把这件事做成可微优化：

| 机器学习用语 | 在本笔记本里 |
| --- | --- |
| 特征 | 实验测到的 \((h,k,l)\) |
| 标签 | 归一化后的实验强度 |
| 可学习参数 | `mode_num` 个模式振幅（本次为 48） |
| 前向网络 | 不是 MLP，而是物理结构因子 `fun_tf` |
| 损失 | 晶体学 R 因子（对 \(\sqrt{I}\) 的相对残差） |
| 优化器 | Adam，全数据一个 batch |

拟合完成后打印每个模式的归一化振幅（相对 `max_mode_amps`），并画出模拟强度 vs 实验强度。

---

## 2. 运行前准备

### 2.1 Python 环境

需要 **Python 3.12.x** 与 TensorFlow。在本 GD 目录：

```powershell
cd <桌面>\GD（未同步git）
python main_requirement.py
```

Jupyter 内核选：

```text
.\.venv\Scripts\python.exe
```

不要用 CRIS 仓库的 `.venv`（一般没有可用的 TensorFlow 3.12 环境）。

### 2.2 本机文件（笔记本写死的路径）

第一个代码单元格里的路径指向作者机器上的 OneDrive。你的电脑上必须能读到同样布局，或先改**工作副本** `GD_modified.ipynb` 里的路径（不要改本参考笔记本）。

配置常量：

```text
structure_type = 'C1'
irrep = 'LD1'
path = D:\OneDrive\PhD\Projects\2511HXRD\Refinement\LD1\C1
sub_name = path\LD1_C1
data_path = D:\OneDrive\PhD\Projects\2511HXRD\Data Process\table\All Combined.csv
```

本机若没有 `D:\OneDrive\...` 那份精修目录，用 CRIS `.venv` 从本地 ISODISTORT 生成一套**同样契约**的文件（不改 tianren 笔记本，也不覆盖 `D:\` 原件）：

```powershell
cd <桌面>\GD（未同步git）
C:\Users\devou\OneDrive\Desktop\CRIS\.venv\Scripts\python.exe isodistort_to_gd.py --irrep LD1 --opd C1 --k LD --k-params 1/6 --nmod 0
```

默认写出：

```text
generated\LD1_C1\LD1_C1_mario_mode_names.txt
generated\LD1_C1\LD1_C1_norm_factors.txt
generated\LD1_C1\LD1_C1_max_bound_vectors.txt
generated\LD1_C1\LD1_C1_displacive_modes.txt
generated\LD1_C1\LD1_C1_alris_functions.py
LD1_C1_alris_functions.py          # 与笔记本同目录，供 import
```

在 **工作副本** `GD_modified.ipynb` 里把 `path` 指到 `generated\LD1_C1`（`sub_name` 仍是 `{path}\{irrep}_{structure_type}`）。`All Combined.csv` 仍是实验衍射表，转换器不会生成。

`sub_name` 会拼出这些文件（缺一则后面单元格会失败）：

| 文件 | 用途 |
| --- | --- |
| `LD1_C1_mario_mode_names.txt` | 统计 `mode_num`（非空行数） |
| `LD1_C1_norm_factors.txt` | 振幅进入结构因子前的缩放 `norm_factors` |
| `LD1_C1_max_bound_vectors.txt` | 每个振幅的绝对上界，Keras `clip_by_value` 用 |
| `LD1_C1_displacive_modes.txt` | 空白分隔、无表头；第 2 列是模式全名 |
| `All Combined.csv` | 实验点，至少含列 `h, k, l, intensity_exp` |

另外必须能 `import LD1_C1_alris_functions`（与笔记本同目录，或已在 `PYTHONPATH`）。它提供：

- `transform_list_hkl_p63_p65`：HKL 网格变换
- `atom_position_list(*amplitudes)`：由模式振幅得到原子坐标
- `get_structure_factors`：结构因子
- `get_structure_factors_CO`：电荷有序相关的另一套结构因子（本笔记本里注释掉了）

### 2.3 建议的运行顺序

在 Jupyter 中 **Restart Kernel**，然后 **Run All**，或按单元格从上到下执行。后面的训练依赖前面的 `features` / `labels` / `fun_tf` / `hkl_indices`。

训练默认 `n_epochs = 1000`、全 batch，第一次会较慢；结构因子在 TensorFlow 图里，后面的 epoch 通常远快于第 1 步。

---

## 3. 单元格在干什么（按执行顺序）

### 3.1 配置与导入

- 设置 `KMP_DUPLICATE_LIB_OK`（Windows 上 MKL/OpenMP 重复加载）。
- `tf.keras.utils.set_random_seed(1)`，便于复现。
- `reload(LD1_C1_alris_functions)`，改完结构因子模块后可重新导入而不必重启内核。
- Debye–Waller 因子现在固定 `w = 0.001`（旧值留在注释里）。
- `intensity_normalization ≈ 26915.86`：后面用它把实验强度除到大约 \([0,1]\)。

### 3.2 HKL 网格与实验点对齐

`fn_all_hkl_list()` 生成完整计算网格：

- \(h \in [-10,10]\)，21 个点（步长 1）
- \(k \in [-4,4]\)，9 个点
- \(l \in [-10-1/6,\,10+1/6]\)，123 个点（步长 \(1/6\)，对应 LD 调制）

实验 CSV 里的 HKL 必须能落到这张网格上。`make_hkl_indices` 把 HKL 乘以 6 再取整来匹配（因为 \(l\) 以 \(1/6\) 为步），返回：

- `hkl_indices`：实验点在全网格中的下标，供 `tf.gather`
- 重排后的 HKL 与强度，保证与网格顺序一致

一次成功运行大约匹配 **4818** 个实验点。

### 3.3 特征与标签

- `features`：HKL 张量，形状 `(n_features, 3)`。
- `labels`：`intensity_exp / intensity_normalization`，再 `expand_dims` 成 `(n_features, 1)`。
- `labels_err`：按 \(1/I\) 做的权重（强度为 0 时用 `1e-9`）。训练函数接收了它，但 `model.fit` 里 `sample_weight` 被注释掉了，**当前损失不加权**。

### 3.4 前向：`fun_tf`

这是物理模型，不是全连接层。

1. 缓存全网格 HKL（经 `transform_list_hkl_p63_p65`），以及 \(|Q|\) 与 Debye–Waller \(\mathrm{e}^{-w|Q|^2}\)。\(l\) 在缓存里曾乘 6，算 \(|Q|\) 时要除回去。
2. 输入振幅 `pars` 先乘 `norm_factors`，再 `atom_position_list` 得到畸变结构。
3. `get_structure_factors` → \(I = |F|^2\) → 乘 DW → 除以网格上的最大强度，把模拟强度也收到约 \([0,1]\)。
4. `tf.gather(..., hkl_indices)` 只取出实验点上的强度，才能和标签比。

电荷有序版本 `get_structure_factors_CO` 留在注释中，默认不用。

### 3.5 Keras 层 `FunAsLayer`

- 权重 `param` 形状 `(mode_num,)`，`RandomNormal(mean=0, stddev=0.05, seed=2)`。
- `constraint`：每个分量裁剪到 \([-max_mode_amps, +max_mode_amps]\)，防止振幅超出物理上界。
- `call`：把当前权重交给 `fun_tf`。

`create_model` 用 `Input(shape=(3,))` 接这一层。HKL 作为输入只是为了走 Keras `fit` 的 API；真正用到的是层内部的 `hkl_indices` 与缓存网格。

### 3.6 损失：R 因子

\[
R = \frac{\sum |\sqrt{I_\mathrm{obs}}-\sqrt{I_\mathrm{calc}}|}{\sum \sqrt{I_\mathrm{obs}}}
\]

`tf.sqrt` 前用 `maximum(..., 1e-8)`，避免强度为 0 时梯度炸掉。`RFactorLoss` 与度量 `r_factor_metric` 是同一公式。笔记本注释提过 MSE，**当前实现用的是 R 因子，不是 MSE**。

### 3.7 训练循环

默认超参：

| 符号 | 值 | 含义 |
| --- | --- | --- |
| `lr` | `[2e-3]` | Adam 学习率列表，可试多档 |
| `n_epochs` | `1000` | 每轮 epoch |
| `n_iter` | `1` | 随机重启次数（局部极小多时可加大） |
| `batch_size` | 全部 4818 点 | 全批量，`shuffle=False` |

对每个学习率、每次重启：新建模型 → `train_single_run` → 按最后一步 loss 保留该学习率下最好的 `best_model_pars`。外层用 `best_rf_overall` 跨学习率再选一次，结果在 `best_pars_overall`。

注意：变量名 `best_rf_overall` 实际比较的是 **loss 标量**（与 R 因子同公式），不是另一套指标。一次完整训练后，日志里大约是 `Best overall loss / Best R-factor ≈ 0.442`。

### 3.8 打印最佳振幅

从 `{sub_name}_displacive_modes.txt` 读模式名，对每个 \(i\) 打印：

```text
Mode a{i+1}  {Mode Name} :  best_pars[i] / max_mode_amps[i]
```

这里打印的是**相对上界的归一化振幅**。`save_best_model_parameters.py`（给 `GD_modified.ipynb` 用）会同时写出未归一化的 `Best Model Parameter` 和 `Normalized Amplitude`；IsoVIZ 应使用前者（与进度条 `amp` 同一单位），不要误用 Normalized Amplitude。

随后画初始权重 vs 最终权重的直方图。

### 3.9 平面图与校验图

`alris_r_factor` 用 DataFrame 再算一遍 R 因子。然后：

- `plot_plane_sim_vs_exp_DIM3_CO`：每个固定 \(h\) 的 \((k,l)\) 平面，**去掉**整数 Bragg 点，看调制/卫星点；半圆左红（模拟）右蓝（实验）。
- `plot_plane_sim_vs_exp_DIM3_bragg`：只保留整数 \(k,l\) 的 Bragg 点。
- 校验图（parity）：横轴实验、纵轴模拟，加 \(y=x\)。
- 相对误差直方图。

\(h\) 从 0 画到 8。

---

## 4. 数据流（简图）

```text
All Combined.csv          LD1_C1_* 模式文件
        │                         │
        ▼                         ▼
  实验 HKL + 强度          mode_num / 上界 / 名字 / 归一化
        │                         │
        └────────► fun_tf ◄───────┘
                     │
                     │  I_sim(hkl) = |F(振幅)|^2 × DW
                     ▼
              R(I_exp, I_sim)
                     │
                     ▼
              Adam 更新 48 个振幅
                     │
                     ▼
              打印 Best model parameters
              （本笔记本到此为止，不写 CSV）
```

`GD_modified.ipynb` 在「ANALYZE MODE AMPLITUDES」之后会调用 `save_best_model_parameters(...)`，写到：

```text
<桌面>/Best_Model_Parameters/LD1_C1/LD1_C1_best_model_parameters.csv
```

再交给 CRIS 的 `ISOVIZ_INPUT/main.py` 写入 `.isoviz`。

---

## 5. 和 `GD_modified.ipynb` 的关系

| | `LD1_C1_gradient_descent_tianren.ipynb` | `GD_modified.ipynb` |
| --- | --- | --- |
| 角色 | 原始参考 | 工作副本 |
| 可否修改 | 否 | 是 |
| 拟合流程 | 本节描述的流程 | 同源，可继续改 |
| 导出 CSV | 无 | 分析振幅后自动导出 |

对照原理或超参时读本文件对应的笔记本即可；不要在参考本上试验。

---

## 6. 常见失败

1. **找不到 `LD1_C1_alris_functions`**：在本目录用 CRIS `.venv` 跑 `isodistort_to_gd.py`，或把已生成的 `.py` 放到本目录。
2. **`FileNotFoundError` on `D:\OneDrive\...`**：数据不在这台机器；只改 `GD_modified.ipynb` 的 `path` / `data_path`。
3. **内核无 TensorFlow**：先 `python main_requirement.py`，再选 GD `.venv`。
4. **HKL 匹配数为 0**：实验表的 \(l\) 必须能对上 \(1/6\) 网格；`SCALE_FACTOR = 6` 的取整要能对上。
5. **振幅全是 0 或顶在边界**：检查 `max_mode_amps` 与 `norm_factors` 是否与 `mode_num` 同长；学习率过大过小都会表现为不降或震荡。
6. **想把结果送进 IsoVIZ**：请运行 `GD_modified.ipynb` 的分析单元格，或在训练完成后对 `best_pars_overall` 调用 `save_best_model_parameters.py`；本参考笔记本不会自动保存。
