# S1–S4 精细材料传播与中心 beta 对比

> 2026-09-13 已完成物理边界修复并重新训练。旧版本使用固定
> `65.72, 61.05, 1.69, -61.50, -70.10 cm`，会切进错层 TOF bar，
> 把部分 S2/S3 闪烁体误归入 Tracker 关联区域；该固定边界版本及其旧图只保留
> 为历史对照，不再作为当前模型。当前实现按**逐事件命中 bar 的真实表面**定义区域，
> 并让材料积分、MC truth 标签和实测沉积变量覆盖一致的物理范围。

## 物理边界定义与标签对齐

方向约定为向下粒子，沿运动方向 z 减小。`CC/Tofdbc.C` 中 TOF 闪烁体厚度
`Sci_t=1 cm`，错层 bar 的中心 `Sci_pz` 为：

| 层 | bar 中心 z（cm） | 向下入口 → 出口（cm） |
|---|---|---|
| S1 | 64.425 或 65.975 | 64.925→63.925 或 66.475→65.475 |
| S2 | 61.325 或 62.875 | 61.825→60.825 或 63.375→62.375 |
| S3 | -62.875 或 -61.325 | -62.375→-63.375 或 -60.825→-61.825 |
| S4 | -64.425 或 -65.975 | -63.925→-64.925 或 -65.475→-66.475 |

当前每个事件从 `tof_barid` 取得命中 bar，并通过 `TOFGeom::GetBarEdge` 计算表面。
五个区域边界保存到训练树 `segmentZ[5]`：

| 边界 | 物理面 | 可能的 z（cm） |
|---:|---|---|
| 0 | 命中 S1 的**入口面** | 64.925 或 66.475 |
| 1 | 命中 S2 的**出口面** | 60.825 或 62.375 |
| 2 | Tracker L5 参考面 | 1.698 |
| 3 | 命中 S3 的**入口面** | -62.375 或 -60.825 |
| 4 | G4 post-S4 acceptance plane | -69.975 |

四个区域因此定义为：

| 区域 | 传播范围 | 包含 / 排除 | 实测涨落变量 |
|---:|---|---|---|
| 0 | S1 入口 → S2 出口 | 完整 S1、层间被动材料、完整 S2 | S1+S2 edep |
| 1 | S2 出口 → L5 中心 | **不含 S2 闪烁体**；上半 Tracker/被动材料 | L2+L3+L4+½L5 edep |
| 2 | L5 中心 → S3 入口 | **不含 S3 闪烁体**；下半 Tracker/被动材料 | ½L5+L6+L7+L8 edep |
| 3 | S3 入口 → -69.975 plane | 完整 S3/S4，并包含 S4 后到该 plane 的被动材料 | S3+S4 edep |

region 3 **不能称为“S3 入口→S4 出口”**。G4 acceptance table 没有与所有错层
S4 bar 中心一一对应的统一末端 checkpoint；为避免用长距离插值冒充精确 MC 真值，
当前统一采用实际存在的 -69.975 cm post-S4 acceptance plane。S4 后的额外被动材料
由材料图均值承担，涨落观测量仍只来自 S3/S4 主动沉积。

训练要求四层 `mtofbar == tof_barid` 且 `mtofpri=1`，避免几何和 MC 沉积来自
不同 bar。S1/S2/S3 的 G4 center checkpoint 是 crossing 后的 post-step 状态，
因此代码使用 checkpoint 的真实 z 与该 bar 的 primary `mtofdep` 在 1 cm 闪烁体内
做局部线性修正，把动能移到对应入口/出口表面。L5 中心改用 primary
`TrMCCluster` 的 `mtrmom[4]`，并要求 `mtrcoo1[4][2]` 位于 1.698±0.03 cm。
末端 checkpoint 17 仍是 post-step 状态，因此 -69.975 cm 端点保留小于单个
transport step 的近似，但它位于同一 post-S4 被动区域，不再跨越 TOF/Tracker 边界。

`amsd69n` 的 `tof_pos` 是 `TOFGeom::GetMeanZ` 层平均平面，并非逐 bar 中心；
当前精细传播的 TOF 时间节点已改用实际命中 bar 的 `Sci_pz`，各层间路径长度仍用
cluster 位置得到的 `tof_leng` 路径差归一化。

本次材料传播仍沿 `g4mscan.root` 原生 z 网格逐步积分；四个区域只是标定与
实测涨落的统计粒度，并不意味着传播只经过四块均匀材料。

## 汇报主线与逐步解释

建议题目：**基于材料均值与实测涨落的 AMS 中心 beta 重建**。
本工作是在既有 nonlinear 重建框架中改造能损预测与传播积分，最终每个事件
仍只拟合中心逆 beta 与共同时间零点；标定参数来自独立 MC 训练，不在单个事件中重拟合。

### 第一步：明确原 nonlinear 的起点

图中的原 nonlinear 对应 `nonlinearBeta`，已经包含 TOF 与 Tracker 的
能量沉积。不要把它描述成“未考虑 Tracker”或“恒速直线拟合”。
它把主动探测器的实测沉积乘以经验 scale，作为传播中的能量损失：

$$
\Delta E_{\rm TOF,i}=\zeta_{\rm TOF}d_i,
\qquad
\Delta E_{\rm Tk,l}=\zeta_{\rm Tk}d_l.
$$

Tracker 层之间采用相应区段能量计算飞行时间，到层位置后更新能量。
因此原方法已经是非线性的：试探中心 beta 决定能量，能损又改变沿途速度，
最后用传播时间与 TOF 时间比较。当前基准配置在每个电荷内给 TOF/Tracker
使用同一个 scale；Z1/Z2/Z6/Z8 分别为 2.4765/2.5352/2.4372/2.1956。

它的限制是：主动硅层或闪烁体的 edep 同时承担了整个路径的平均能损代理和
事件涨落代理。被动支撑等材料对实际能损的贡献，只能间接被经验 scale 吸收。
原方法还有一个独立输出 `nonlinearBetaLegacy`，不是本次蓝色比较曲线。

### 第二步：让材料图承担平均能损

Hybrid 使用重建轨迹和内层刚度穿过九元素材料图，缓存每个小步的元素丰度。
在当前试探 beta 下，用既有 Bethe 计算得到小步材料预测，再按区域求和：

$$
\mu_s^{\rm map}(\beta,m,\mathbf g)
=\sum_{j\in s}\mu_j^{\rm map}(\beta,m,\mathbf g),
\qquad B_s=\kappa_s\mu_s^{\rm map}.
$$

这里的几何变量 $\mathbf g$ 包含事件的轨迹与路径；两个同 Z、同 beta 的
事件若穿过不同材料，材料均值仍可不同。材料缓存每个事件只建立一次；
每次试探 beta 时重新计算能损和时间，不重新拟合轨迹。

`MaterialMean` 现在把事件的入射电荷传给 `BetheEquation`，因此材料基准本身
已经包含 Bethe 方程的 projectile-charge $Z^2$ 因子。$\kappa$ 只继续承担
材料图、连续能损近似和训练标签之间的剩余定标，不再用于吸收从 Z1 到重核的
主导 $Z^2$ 依赖。旧的单位电荷标定参数与这一实现不兼容，必须重新训练。

### 第三步：实测 edep 只提供中心化涨落修正

对区域 $s$，定义实测沉积汇总 $d_s$，模型为：

$$
\widehat{\Delta E}_s
=\kappa_s\mu_s^{\rm map}+a_s(d_s-\bar d_s).
$$

| 配置字段 | 符号 | 作用 |
|---|---|---|
| `materialScale` | $\kappa_s$ | 定标材料预测的平均能损 |
| `meanDeposit` | $\bar d_s$ | 定义相应训练箱中的平均实测沉积 |
| `fluctuationScale` | $a_s$ | 将实测沉积相对均值的偏离映射到能损残差 |

减去 $\bar d_s$ 是为了避免把平均 edep 再加一次。这里“中心化”只保证
所选训练箱内的沉积偏离均值为零；换样本、换几何分布以及输出裁剪以后，
不能宣称该项的均值仍严格为零。

仅为说明公式：若材料项为 30 MeV，实测沉积比箱均值高 2 MeV，且
$a=0.5$，则预测为 31 MeV。此例不是实际 MC 事例。

### 第四步：四个区域与精细传播是两种粒度

区域边界不是固定常数，而由每个事件的命中 TOF bar 决定；完整定义见上面的
“物理边界定义与标签对齐”。region 1/2 明确从 S2 出口开始、到 S3 入口结束，
因此不再包含 TOF 闪烁体。L5 位于两个 Tracker 区域的共同边界，实测 L5 沉积
按 1/2 + 1/2 分给两侧。

四段是统计标定区域，不是四个独立 beta 拟合量，也不是传播只经过四块均匀材料。
真实时间积分仍在原生 1 cm 材料网格上进行，并在 Tracker 层、实际命中 bar
中心、L5 参考面及区域边界处分步；0.5 cm 仅作为数值收敛检查，不会创造
低于材料图原生分辨率的新结构信息。

### 第五步：把区域模型分配到小步

设当前小步 $j$ 属于区域 $s$，按当前局部 beta 计算：

$$
w_j=\frac{\mu_j^{\rm map}}{\mu_s^{\rm map}},
\qquad
\widehat{\Delta E}_j
=w_j\,\mathrm{clip}\!\left[
\kappa_s\mu_s^{\rm map}+a_s(d_s-\bar d_s),\ 0,\ 3\kappa_s\mu_s^{\rm map}
\right].
$$

若不触发裁剪，它等于 $\kappa_s\mu_j^{\rm map}+w_j a_s(d_s-\bar d_s)$。
材料多的小步分配较多修正，真空步只累计时间。
这是**区域涨落按材料权重分配**的假设，不是逐小步测到了真实涨落；
当前也不是只在主动探测器位置注入涨落。因为每步 beta 会更新，权重与
系数按该步局部 beta 重新求值，不应把整次传播理解成预先固定的四个损失数。

### 第六步：从中心 beta 向上下传播，预测四个 TOF 时间

对一个试探 $\beta_c$，用项目的按电荷平均质量假设得到：

$$
E_c=\frac{m}{\sqrt{1-\beta_c^2}},
\qquad \beta(E)=\sqrt{1-\frac{m^2}{E^2}}.
$$

中心往 S3/S4 传播时扣能量；中心往 S2/S1 反推时恢复能量。
每步先算一次初值能损，用中点能量重估能损，再计算该步时间：

$$
E_{j+1}=E_j\mp\widehat{\Delta E}_j(E_{\rm mid}),
\qquad
\Delta t_j=\pm\frac{\Delta l_j}{2c}
\left[\frac{1}{\beta(E_j)}+\frac{1}{\beta(E_{j+1})}\right].
$$

时间向下为正、向上为负。所有小步累计后得到 $T_i(\beta_c)$，预测 TOF
时间为 $t_i^{\rm pred}=\tau+T_i(\beta_c)$。小步时间路径按实测 TOF 路径差
归一化，能损影响速度；没有每次 trial 都重新拟合整条轨迹。

最终仍最小化：

$$
\chi^2(\beta_c,\tau)=
\sum_i\frac{[t_i^{\rm meas}-\tau-T_i(\beta_c)]^2}{\sigma_{t,i}^2}.
$$

单事件拟合只有 $1/\beta_c$ 与 $\tau$ 两个参数。材料参数来自预训练，
生产重建不读取 MC 真值来决定拟合结果；MC 真值只用于训练标签和评价。
训练材料项使用 MC 质量/局部 beta，重建使用项目的平均质量假设。

上述能量公式用于 beta<1 的试探。beta≥1 时沿用既有恒速预测分支，
允许 TOF 测量分辨率造成的超 1 拟合值；这不是本次新增物理模型。

几何参考面需区分训练与评价：精细 hybrid 的传播中心为 Tracker L5 的
z=1.698 cm，训练标签在这里使用 primary `TrMCCluster` 动量；当前 `betaTree.mcBeta`
为了与历史 linear/nonlinear 输出保持可比，评价仍沿用 acceptance checkpoint 10。
两者非常接近但并非同一个 transport 状态，因此不能把这项差异宣称为严格为零；
如需最终物理精度结论，应统一三种方法的 MC reference 定义后再做系统学比较。

### 第七步：三个能损参数具体如何训练

Option 4 构造与区域边界对齐的 MC 动能 $K=\sqrt{p^2+m^2}-m$，区域标签为
$y=K_{\rm in}-K_{\rm out}$，并导出材料预测 $x$、实测沉积 $d$、局部入口 beta。
S1/S2/S3 根据 bar 奇偶选择 G4 center checkpoint（3/4、5/6、14/15），再用
实际 post-step z 与 primary `mtofdep` 修正到闪烁体表面；L5 使用 primary
`mtrmom[4]`；末端使用 checkpoint 17。训练树同时保存 `segmentZ`、`tofBar`、
`mcTofBar`、`truthBoundaryKinetic`、`truthCheckpointZ` 和 `truthCheckpoint`，
便于逐事件审计边界与 truth 来源。

Option 5 对每个 Z、每个区域独立按局部 beta 分箱，当前各 8 箱。
在每个箱里，用 $q=y/x$ 的 median 与 $1.4826\,MAD(q)$ 构造宽度为
10 倍稳健 sigma 的窗口，压低极端损失尾对标定的影响。它是统计去尾，
不等于逐事件识别了核相互作用。

对保留下来的事件，代码依次计算：

$$
\kappa=\frac{\sum_i y_i}{\sum_i x_i},\qquad
\bar d=\frac{1}{N}\sum_i d_i,\qquad
r_i=y_i-\kappa x_i,
$$
$$
a_{\rm fit}=\frac{\sum_i(d_i-\bar d)r_i}{\sum_i(d_i-\bar d)^2}.
$$

只有沉积偏离与真实残差的相关系数 $\rho>0.1$ 时才保留该斜率，否则令
$a=0$。0.1 是效应大小的经验门槛，不是 p-value 或严格显著性检验。
因此本模型是分箱定标和一元残差回归，不是神经网络；$\kappa$ 也不是
联合最小二乘拟合出来的参数，它使用的是上面的均值比定义。

共 4 电荷 × 4 区域 × 8 beta 箱 = 128 组，每组 3 个预测参数。
计数、RMS 与相关系数列是诊断信息，不增加单事件拟合自由度。

### 第八步：把数值保护与物理模型改进分开

当前还有三项工程控制：区域预测裁剪至 $[0,3\kappa\mu^{\rm map}]$；
hybrid 的逆 beta 搜索区间限制到 linear 初值的 ±0.25，并与旧的
$[0.6,10]$ 区间取交集；精细传播出现不物理能量时返回大 chi2 惩罚。

原 nonlinear 的比较曲线没有同步使用这个较窄的搜索窗口。因此普通 RMS
或失败尾的差异含有优化器约束的作用，不能作为材料模型单独改善的证据。
已有结果是组合实现的比较；要分离材料均值、涨落与步长的独立贡献，还需
统一拟合控制后的消融试验，本报告不将尚未做的消融当作已有结论。

## 简单汇报的讲解顺序（约 6 页）

| 页 | 内容 | 可以直接讲的一句话 |
|---:|---|---|
| 1 | 目标与旧方法 | 我们保持中心 beta 的 TOF 拟合目标，改进从中心到四个 TOF 的能损与时间预测。 |
| 2 | Hybrid 公式与训练 | 用材料图描述平均能损，实测沉积只修正相对训练均值的涨落。 |
| 3 | 四区域与精细传播 | 四段按逐事件物理边界标定，实际在原生材料网格逐步更新能量和时间。 |
| 4 | 核心均值随 beta | 低 beta 的部分偏差减小，最低 beta 和重核过修正仍存在。 |
| 5 | 核心 sigma 与残差叠加 | 核心宽度相对原 nonlinear 没有一致改善，不能声称分辨率全面提升。 |
| 6 | 验证与后续 | 数值步长主体稳定；下一步需要统一拟合控制并做模型消融。 |

汇报用结果图已在文末列出，三页 PDF 的顺序为均值、sigma、低 beta 分布。

## 实现与参考定义

- 使用既有 `g4mscan.root` 九种元素材料密度与重建轨迹、内层刚度。
  按每个 z voxel 的实际边界积分，在 TOF 测时面、Tracker 层和中心参考面
  插入节点；每小步取五个横向采样点。原 API 会把端点取到 voxel 中心，
  所以不能通过反复调用短区间 API 得到等价的精细积分。
- 每步以局部能量计算材料能损，用中点预测修正能量，用端点逆 beta 的
  梯形积分计算飞行时间。沿中心向上恢复能量、向下扣除能量；时间积分在
  四个实际 TOF 测时面取值，S4 之后的传播不进入 TOF 时间。
- 同一区域的实测涨落按小步材料权重分配；被动结构没有新增独立的实测
  涨落变量。使用四区域 × 四电荷（Z1/Z2/Z6/Z8）× 八局部 beta 箱的 128 组标定。
- 小步路径按各 TOF 间的实测路径差归一化。Hybrid 训练的中心 truth 使用
  primary L5 `TrMCCluster`；为保持历史比较口径，当前输出树的 `mcBeta` 仍沿用
  acceptance checkpoint 10。`linearBeta` 与原 `nonlinearBeta` 仍以此前算法计算。
- 输出 `hybridFitStatus`、`hybridChi2` 和 `hybridSteps` 便于检查拟合。
  保留旧的有限拟合输出惯例；本报告的共同事件筛选不额外要求 Minuit status=0。

## 样本与复现

`config/hybrid_material_fine.sources` 固定了 Z1/Z2/Z6/Z8 各 24 个训练文件和
24 个独立验证文件，共 192 个完整 MC 文件。2026-09-13 的物理边界版本重新从
原始 MC 生成 Option 4 训练树，没有复用旧固定边界 training ROOT。

Option 4 导出的可审计样本数为：Z1 172,056、Z2 100,468、Z6 34,958、
Z8 17,946，共 325,428。Option 5 还要求所有局部 beta、材料预测、truth loss
和 measured deposit 满足训练有效性条件，最终进入标定的是 Z1 171,276、
Z2 100,468、Z6 34,958、Z8 17,946，共 **324,648** 条记录，生成
4 电荷 × 4 区域 × 8 beta 箱 = 128 行参数。验证文件未参与参数拟合。

96 个独立验证文件全部重新重建。要求重建电荷等于 MC 电荷，且 linear、原
nonlinear、Hybrid 三种 beta 都有限且为正后，共同样本数为：
Z1 **193,655**、Z2 **118,554**、Z6 **45,894**、Z8 **24,747**，共 **382,850**。

训练分片和 Option 5 使用项目 ROOT5 环境；结果合并和绘图使用干净 shell 的
ROOT 6.36.10。ROOT6 写出的合并文件不再交回 ROOT5 训练器。

```bash
./build.sh
./run.sh input.root training_part.root 4 1 1 all center none none 1
./run.sh training.list config/hybrid_material_fine.conf 5
./run.sh input.root beta.root 0 1 1 all center \
  config/energy_loss_scales_center_common_z1268.conf \
  config/hybrid_material_fine.conf 1

root -l -b -q 'macro/plotChargeBetaComparison.C("results/hybrid_physical_20260913/merged/beta_Z2.root","results/hybrid_physical_20260913/merged/beta_Z6.root","results/hybrid_physical_20260913/merged/beta_Z8.root","macro/results/beta_comparison/hybrid_material_physical.pdf",100,"results/hybrid_physical_20260913/merged/beta_Z1.root",true)'
```

第十个参数为最大材料 z 步长（cm）。正式 fine 配置使用 `1`；`0` 使用 coarse
传播并配套同样已按物理边界重训的 `config/hybrid_energy_loss_b1308_z268_20260904.conf`。
Z1 原 nonlinear 的共同 TOF/Tracker scale 仍为 `2.47651231289`；Z2/Z6/Z8
沿用 `config/energy_loss_scales_center_common_z1268.conf` 中既有值，这些 scale
只控制原 nonlinear 比较曲线，不进入 Hybrid 材料标定。

## 数值收敛

在每个 Z 的一个独立验证文件上，用同一物理边界和同一新配置分别运行
`materialStepCm=1` 与 `0.5`，比较两个结果都有限时的
`|1/beta_0.5cm - 1/beta_1cm|`：

| Z | 共同事件 | 绝对差中位数 | P95 | 最大绝对差 |
|---:|---:|---:|---:|---:|
| 1 | 8,254 | 4.48e-9 | 7.01e-7 | 0.0368 |
| 2 | 5,067 | 4.64e-9 | 7.35e-7 | 0.00188 |
| 6 | 1,548 | 4.94e-9 | 1.35e-6 | 0.000158 |
| 8 | 1,013 | 3.70e-9 | 5.87e-7 | 0.000621 |

主体分布对减半步长非常稳定，但 Z1 存在一个明显的极端敏感事件，因此只能说
**分布层面收敛**，不能声称每个事件都收敛。该检查验证数值积分稳定性，不代表
材料图具有低于原生 1 cm 网格的真实结构信息。

## 对比结果

残差定义为 `1/beta_reco - 1/beta_MC`。当前 `mcBeta` 为历史比较口径的
acceptance checkpoint 10；因此下面是重建性能比较，不应误写成“相对于精确 L5
同一物理面的绝对偏差”。随 beta 图每箱至少 100 个共同事件，低 beta 图使用
`beta_MC < 0.9`，并在稳健窗口内做高斯核心拟合。

| Z | 方法 | 低 beta 核心均值 mu | 核心 sigma |
|---:|---|---:|---:|
| 1 | Linear | 0.006501 ± 0.000155 | 0.025965 ± 0.000127 |
| 1 | 原 Nonlinear | 0.005771 ± 0.000152 | 0.025349 ± 0.000126 |
| 1 | Hybrid physical | **0.004899 ± 0.000153** | 0.025371 ± 0.000128 |
| 2 | Linear | 0.002673 ± 0.000107 | 0.014880 ± 0.000086 |
| 2 | 原 Nonlinear | 0.002165 ± 0.000103 | 0.014444 ± 0.000082 |
| 2 | Hybrid physical | **0.001023 ± 0.000104** | 0.014310 ± 0.000085 |
| 6 | Linear | 0.001036 ± 0.000128 | 0.009700 ± 0.000098 |
| 6 | 原 Nonlinear | 0.000717 ± 0.000123 | 0.009533 ± 0.000096 |
| 6 | Hybrid physical | **-0.000590 ± 0.000125** | 0.009550 ± 0.000099 |
| 8 | Linear | 0.000901 ± 0.000179 | 0.009545 ± 0.000138 |
| 8 | 原 Nonlinear | **0.000605 ± 0.000172** | 0.009385 ± 0.000130 |
| 8 | Hybrid physical | -0.000668 ± 0.000176 | 0.009457 ± 0.000138 |

物理边界 Hybrid 相对 linear 在四个 Z 上都把低 beta 核心均值向 0 拉近；相对原
nonlinear，Z1/Z2 明显改善，Z6 的绝对偏差略有改善但转为负向，Z8 则出现轻微
负向过修正且绝对核心偏差略大于原 nonlinear。核心 sigma 与原 nonlinear 接近，
不能声称分辨率全面提升。

普通全样本 RMS 还会混入 Minuit 失败尾和 Hybrid 专用 ±0.25 逆 beta 搜索窗口的
影响；因此本节主要使用核心拟合判断材料模型，不把尾部改善单独归因于边界或材料
模型。若要给出最终系统学结论，还需要统一三种方法的拟合窗口和 MC reference 面。

## 产物

当前物理边界版本：

- `config/hybrid_material_fine.conf`：正式 fine 标定，charge-aware Bethe，逐 bar 物理边界。
- `config/hybrid_energy_loss_b1308_z268_20260904.conf`：同样按物理边界重训的 coarse Z2/Z6/Z8 标定；文件名保留历史日期以避免破坏已有调用路径。
- `results/hybrid_physical_20260913/merged/beta_Z1.root`、`beta_Z2.root`、`beta_Z6.root`、`beta_Z8.root`：96 个独立 validation 文件按电荷合并后的重建结果。
- `macro/results/beta_comparison/hybrid_material_physical.pdf` / `.png`：三方法随 beta 的核心均值与 sigma。
- `macro/results/beta_comparison/hybrid_material_physical_mean.pdf` / `.png`：2×2 电荷偏差页。
- `macro/results/beta_comparison/hybrid_material_physical_sigma.pdf` / `.png`：2×2 电荷分辨率页。
- `macro/results/beta_comparison/hybrid_material_physical_low_beta.pdf` / `.png`：低 beta 分布页。
- `macro/results/beta_comparison/hybrid_material_physical_residuals.pdf` / `.png`：低 beta 残差叠加与高斯拟合。
- `macro/results/beta_comparison/hybrid_material_physical.root`：绘图直方图和拟合点。
- `macro/results/beta_comparison/hybrid_material_physical_report.pdf`：均值、sigma、低 beta 三页合并报告。
- `results/hybrid_physical_20260913/convergence/`：1 cm / 0.5 cm 数值收敛检查。

旧的 `hybrid_material_fine*` 图和旧 merged ROOT 对应修复前固定边界原型，只能作为
历史对照，不应再引用为当前模型结果。

