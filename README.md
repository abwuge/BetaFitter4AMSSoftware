# AMS TOF 粒子速度非线性重建

从四层 TOF 时间重建粒子在 AMS 中心或 TOF 入口处的速度。原 Nonlinear 将 TOF 与 Tracker 的实测沉积转换为传播能损；Hybrid 用材料图预测平均能损，并用中心化的实测沉积修正事件涨落。

## 编译与单文件运行

在 sdiat 的项目目录中：

```bash
./build.sh
./run.sh input.root output.root 0 1.9 1.0 all center
```

`build.sh` 和 `run.sh` 加载 AMSVX/ROOT5 环境。合并结果与绘图使用干净 shell 中的 ROOT6。编译可选 `./build.sh debug` 或 `./build.sh clean`。

工程依赖外部 AMSVX 软件、AMS 数据目录和材料图 `$AMSDATADIR/LAPP/dEdxPDF/g4mscan.root`。脚本当前从 `/home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh` 加载环境。输入为 `amstreea` ntuple；`input_Z*.list` 和 `config/*.sources` 保存生产输入路径。代码归档保留清单，原始 ROOT 数据与 AMS 软件在各自存储中。

`run.sh` 参数按以下顺序传入：

```text
input output [option=0] [tofScale=1] [trackerScale=1]
[scaleMode=all] [reference=center] [scaleConfig=none]
[hybridConfig=none] [materialStepCm=0]
```

| 参数 | 取值与作用 |
|---|---|
| `option` | `0` 重建；`1` 固定 MC beta 拟合共同 scale；`2` 计时；`3` 导出 beta 差；`4` 导出 Hybrid 训练样本；`5` 训练 Hybrid 参数；`-1/-2` 导出磁场/能损 |
| `tofScale`, `trackerScale` | 分别缩放原 Nonlinear 的 TOF、Tracker 能损。直接调用可执行文件时 TOF scale 默认是 `2`，`run.sh` 默认是 `1` |
| `scaleMode` | `all` 缩放 S1/S2/S3；`s1s2` 缩放 S1/S2；`s2` 仅缩放 S2，其余沉积系数固定为 1 |
| `reference` | `center` 重建中心 beta；`before_tof` 重建 TOF 入口 beta |
| `scaleConfig` | 三列 `charge tofScale trackerScale`，仅处理配置中的电荷 |
| `hybridConfig` | Hybrid 标定文件，`none` 关闭 Hybrid |
| `materialStepCm` | `0` 使用 coarse 传播；正数指定精细材料积分的最大 z 步长，正式 fine 配置使用 `1` cm |

输出 `betaTree` 包含 `linearBeta`、含 Tracker 的 `nonlinearBeta`、旧 TOF-only 输出 `nonlinearBetaLegacy`、`hybridBeta`、`hybridFitStatus`、`hybridChi2`、`hybridSteps`、`mcBeta`、`innerRigidity`、`betaRigidity`、`Z` 和 `mcCharge`。未启用 Hybrid 时 `hybridBeta` 为 NaN，状态为 -1。

## 当前 Hybrid 配置

```bash
./run.sh input.root beta.root 0 1 1 all center \
  config/energy_loss_scales_center_common_z1268.conf \
  config/hybrid_material_fine.conf 1
```

原 Nonlinear 的共同 TOF/Tracker scale 按 Z1/Z2/Z6/Z8 配置。Hybrid 的平均能损由包含入射电荷 Z² 的 Bethe 材料预测给出；区域模型为：

```text
loss = materialScale * materialMean
     + fluctuationScale * (measuredDeposit - meanDeposit)
```

向下粒子的四个区域按实际命中 TOF bar 的物理表面划分：

| 区域 | 范围 | 实测沉积变量 |
|---:|---|---|
| 0 | S1 入口 → S2 出口 | S1 + S2 |
| 1 | S2 出口 → Tracker L5 中心 z=1.698 cm | L2+L3+L4+½L5 |
| 2 | L5 中心 → S3 入口 | ½L5+L6+L7+L8 |
| 3 | S3 入口 → post-S4 acceptance plane z=-69.975 cm | S3 + S4 |

材料积分在原生材料网格及探测器节点处分步，每事件拟合中心逆 beta 与共同时间零点。区域参数来自独立训练。正式 fine 配置含 128 行参数，覆盖 4 电荷 × 4 区域 × 8 局部 beta 箱。

训练与重建入口：

```bash
./run.sh input.root training_part.root 4 1 1 all center none none 1
find hybrid_training_parts -name '*.root' | sort > training.list
./run.sh training.list config/hybrid_material_fine.conf 5
```

`hybridTraining` 保存材料预测、truth loss、实测沉积及边界审计字段。训练用 primary L5 truth；历史比较输出 `mcBeta` 使用 acceptance checkpoint 10。当前三方法结果沿用这一评价口径。

[Hybrid 方法、训练和结果](docs/hybrid-material-fine.md)给出完整定义。`config/hybrid_energy_loss_b1308_z268_20260904.conf` 为按相同物理边界重训的 coarse 配置，文件名保留原日期。

## 结果与分析入口

当前物理边界结果在 `macro/results/beta_comparison/hybrid_material_physical*`；独立验证合并数据在 `results/hybrid_physical_20260913/merged/`。旧 `hybrid_material_fine*` 图对应修正前固定边界版本。

共同有效 MC 样本 Z1/Z2/Z6/Z8 为 193,655 / 118,554 / 45,894 / 24,747。以 `1/beta_reco - 1/beta_MC` 为残差，Hybrid 在低 beta 样本中改善 Z1/Z2 的核心均值，Z6/Z8 出现轻微负向偏差。核心宽度与原 Nonlinear 接近。

使用 ROOT6 复现比较图：

```bash
root -l -b -q 'macro/plotChargeBetaComparison.C("results/hybrid_physical_20260913/merged/beta_Z2.root","results/hybrid_physical_20260913/merged/beta_Z6.root","results/hybrid_physical_20260913/merged/beta_Z8.root","macro/results/beta_comparison/hybrid_material_physical.pdf",100,"results/hybrid_physical_20260913/merged/beta_Z1.root",true)'
```

| 资料 | 内容 |
|---|---|
| `docs/energy-loss-scale-scope.md` | all、s1s2、s2 的同事件比较与累计能损解释 |
| `docs/hybrid-material-fine.md` | 正式物理边界、charge-aware Bethe、训练、数值收敛与性能表 |
| `docs/hybrid-energy-loss-b1308-z268-20260904.md` | coarse Hybrid 的方法与更新说明 |
| `research/zeta_study_20260722/` | ζ 的全局拟合与时间核分解 |
| `research/nonlinear_bias_decomposition_20260802/` | 传播和 TOF 时间两项相反偏差的分解 |
| `research/real_data_validation_20260722/` | 四层盲验证、独立时间对闭合、RICH 和留一验证 |

`research/` 收录研究源码、说明与小型数值摘要；生成的 ROOT、图像、二进制和批量运行日志保留在工作目录中。`ref/LxBetalhd.*` 为外部实现参考，运行使用 AMSVX 提供的接口。

旧式批量脚本 `run_local.sh` 保留原有八参数接口，可运行 Linear/Nonlinear；Hybrid 批量重建请逐文件调用上述 `run.sh`。独立的后续估计器在相邻 `betaEstimator` 项目中。
