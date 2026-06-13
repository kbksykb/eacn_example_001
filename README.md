# Harmony + REAL 流水线使用说明

本仓库根目录提供两个独立入口脚本，对应 **任务一（benchmark 整合）** 与 **任务二（REAL 评估流水线）**。二者从同一份原始 `.h5ad` 出发，参数与输出目录分离，**互不覆盖**。

| 脚本 | 用途 | 核心代码 |
|------|------|----------|
| `run_task1_harmony.py` | HVG 2000 + PCA 50 + Harmony | `workspace/code/harness/` |
| `run_task2_harmony_real.py` | HVG 1000 + Leiden + Harmony 32 维 + REAL | `workspace/code/pilots/synth_pilot.py` |

辅助模块：`pipeline_common.py`（路径解析、参数默认值，无需直接运行）。

---

## 1. 环境准备

### 1.1 Conda / Python

推荐在仓库克隆目录下激活环境（服务器示例：`conda activate scenv`）。

### 1.2 依赖（任务一）

```bash
pip install anndata scanpy harmonypy==0.0.10
```

> Windows + Python 3.13 上较新的 `harmonypy` 可能编译失败，请固定 `0.0.10`。

### 1.3 额外依赖（任务二）

```bash
pip install torch pyarrow igraph leidenalg
```

任务二 REAL 的 OT 通道需要 `torch`；Leiden 聚类需要 `igraph` 与 `leidenalg`。

### 1.4 运行位置

**必须在仓库根目录**（与 `workspace/` 同级）执行，以便 `workspace.code.*` 包路径正确：

```bash
cd /path/to/eacn_example_001
python run_task1_harmony.py --input workspace/Galaxy1.h5ad
```

---

## 2. 输入文件放在哪里

输入可以是**任意路径**的 `.h5ad`，通过 `--input` 指定：

| 方式 | 示例 |
|------|------|
| 相对仓库根目录 | `--input workspace/Galaxy1.h5ad` |
| 推荐目录 | `--input workspace/data/input/MyDataset.h5ad` |
| 绝对路径 | `--input C:/data/MyDataset.h5ad` |

### 2.1 输入 AnnData 要求

- **必需**：`obs` 中含批次列，默认名为 `batch`（可用 `--batch-key` 修改）。
- **规模**：脚本在内存中读入整表；Galaxy1 示例为 5220×33694。
- **格式**：AnnData `.h5ad`（Scanpy 可读）。

### 2.2 推荐目录结构

```
eacn_example_001/
├── run_task1_harmony.py
├── run_task2_harmony_real.py
├── pipeline_common.py
├── workspace/
│   ├── Galaxy1.h5ad              # 可放这里（已在 .gitignore）
│   └── data/
│       └── input/                # 推荐统一放输入
│           └── MyDataset.h5ad
└── output/                       # 本机默认输出根目录
    ├── task1/
    └── task2/
```

大文件 `.h5ad` 通常不入库；请将数据放在本地 `workspace/` 或 `workspace/data/input/`。

---

## 3. 输出放在哪里

默认输出根目录：**`./output/`**（仓库根下）。

优先级：

1. 命令行 `--out-root /your/path`
2. 环境变量 `EACN_OUT_ROOT`
3. 默认 `<repo>/output/`

任务一与任务二分别写入 **`task1/`** 与 **`task2/`** 子目录，避免互相覆盖。

### 3.1 任务一输出树

```
<out-root>/task1/
├── integrations/harmony/<name>.h5ad    # obsm['X_integrated']，50 维（与 PCA 一致）
├── runtime/harmony_<name>.json       # harness 运行时指标
└── runtime/task1_<name>_manifest.json # 汇总清单（含输入/输出绝对路径）
```

**整合规格（harness）**：

- 2000 HVG（`HarnessConfig.n_top_genes`）
- PCA 50 维
- Harmony 作用于 `X_pca`，输出维数与 PCA 相同

### 3.2 任务二输出树

```
<out-root>/task2/
├── data/<name>_pre.h5ad                         # 预处理后（含 X_uncorrected_pca、motif_id_pre）
├── integrations/harmony/<name>.h5ad             # 32 维 X_integrated
├── detections/real/harmony_<name>.parquet       # 每个 motif 的 REAL 通道分数
└── runtime/harmony_real_<name>.json             # 汇总（LossRate@1/3、耗时等）
```

**流水线规格（synth_pilot）**：

- 1000 HVG + log 归一化 + scale
- PCA 50 → 存为 `X_uncorrected_pca`
- Leiden 过聚类（默认 resolution=3.0）→ `motif_id_pre`
- Harmony → PCA 投影至 **32 维**
- REAL 三通道 + Fisher 融合 → parquet

> **注意**：任务二**不会**读取任务一的 h5ad；两条 Harmony 实现与超参不同，属设计内行为。

### 3.3 服务器示例

在 8×A100 服务器上可将输出指到 shared 盘：

```bash
export EACN_OUT_ROOT=/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared
python run_task1_harmony.py --input workspace/data/input/Galaxy1.h5ad
python run_task2_harmony_real.py --input workspace/data/input/Galaxy1.h5ad
```

实际文件会落在 `$EACN_OUT_ROOT/task1/` 与 `$EACN_OUT_ROOT/task2/`。

---

## 4. 命令行参数详解

### 4.1 公共参数（两个脚本均支持）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--input` | **必填** | 输入 `.h5ad` 路径 |
| `--name` | 输入文件 stem | 输出文件名中的数据集名，如 `Galaxy1` |
| `--batch-key` | `batch` | `obs` 中批次列名 |
| `--out-root` | `$EACN_OUT_ROOT` 或 `./output` | 输出根目录 |
| `--seed` | `0` | 随机种子 |

### 4.2 任务二专用参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--leiden-resolution` | `3.0` | Leiden 分辨率，控制 motif 候选数量 |
| `--cell-type-key` | 无 | 可选；用于 parquet 中 `ground_truth_label` 列 |
| `--rare-label` | 无 | 与 `--cell-type-key` 联用，标记“真稀有”细胞类型字符串 |

未提供 `--cell-type-key` 时，`ground_truth_label` 全为 `False`（不影响 REAL 打分，仅便于事后对照）。

---

## 5. 使用示例

### 5.1 任务一 — 仅 Harmony 整合

```bash
python run_task1_harmony.py \
  --input workspace/Galaxy1.h5ad \
  --name Galaxy1 \
  --batch-key batch \
  --seed 0
```

成功后在终端打印 JSON manifest，并写入 `output/task1/`。

### 5.2 任务二 — Harmony + REAL

```bash
python run_task2_harmony_real.py \
  --input workspace/Galaxy1.h5ad \
  --name Galaxy1 \
  --batch-key batch \
  --leiden-resolution 3.0 \
  --seed 0
```

若数据含细胞类型注释且需标记 ground truth（可选）：

```bash
python run_task2_harmony_real.py \
  --input workspace/Galaxy1.h5ad \
  --cell-type-key celltype \
  --rare-label "SomeRareType"
```

### 5.3 自定义输出目录

```bash
python run_task1_harmony.py --input workspace/Galaxy1.h5ad --out-root D:/runs/eacn
python run_task2_harmony_real.py --input workspace/Galaxy1.h5ad --out-root D:/runs/eacn
```

---

## 6. 输出文件字段说明

### 6.1 任务一 `integrations/harmony/<name>.h5ad`

- `obsm['X_integrated']`：Harmony 校正后的嵌入（**50 维**）
- 其余 `obs` / `var` 与输入一致

### 6.2 任务二 `data/<name>_pre.h5ad`

- `obsm['X_uncorrected_pca']`：整合前 PCA（50 维）
- `obsm['X_pca']`：与上相同（供 neighbors/Leiden）
- `obs['motif_id_pre']`：Leiden 簇 ID

### 6.3 任务二 `integrations/harmony/<name>.h5ad`

- `obsm['X_integrated']`：**32 维** Harmony 嵌入
- `uns['method']` / `uns['method_version']` / `uns['pipeline']`

### 6.4 任务二 `detections/real/harmony_<name>.parquet`

每行一个 motif 候选，主要列包括：

- `candidate_id`, `abundance`, `loss_probability`
- `channel_mknn_purity_pre/post`, `channel_proc_displacement`, `channel_boot_stable`
- `channel_ot_stat`, `fisher_p_value`
- `ground_truth_label`（可选，见 §4.2）

### 6.5 运行时 JSON

- 任务一：`runtime/harmony_<name>.json`（harness）+ `task1_<name>_manifest.json`
- 任务二：`runtime/harmony_real_<name>.json`（含 `LossRate@1`、`LossRate@3`）

---

## 7. 任务一 vs 任务二对比

| 项目 | 任务一 | 任务二 |
|------|--------|--------|
| HVG 数量 | 2000 | 1000 |
| Harmony 输出维 | 50 | 32 |
| Leiden / REAL | 无 | 有 |
| 典型用途 | 与其他整合方法 benchmark 对齐 | 稀有亚群丢失风险检测 |
| 输出子目录 | `<out-root>/task1/` | `<out-root>/task2/` |

---

## 8. 常见问题

**Q：两个任务可以并行跑吗？**  
可以。输出目录已分离，只要 `--input` 指向同一份只读 h5ad 即可。

**Q：任务二能否接着读任务一的 h5ad？**  
当前设计**不支持**也**不需要**。任务二从原始 count 矩阵重新预处理。

**Q：Windows 上报 `ImportError: resource`？**  
`synth_pilot.py` 已对 `resource` 模块做可选导入；请更新到当前分支版本。

**Q：Harmonypy 安装失败？**  
使用 `pip install harmonypy==0.0.10`。

**Q：REAL 很慢？**  
OT 置换默认 200 次、CPU 设备；可在源码 `synth_pilot.run_real_channels` 中调整 `OTChannelConfig`（脚本层未暴露该参数，后续可扩展）。

---

## 9. 相关文档与代码

- 实验布局：`workspace/notes/experimental_plan_v1.md`
- Harness 设计：`workspace/code/harness/__init__.py`
- REAL 流水线：`workspace/code/pilots/synth_pilot.py`
- 服务器状态：`workspace/STATE.md`
