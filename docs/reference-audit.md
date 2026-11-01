# 实测 RICH / Tracker 参考量审计

审计范围：只读 ssh:sdiat 当前 producer 与 AMS library 源码；只读取 flight 文件
/lustre02/user/yachen/public/hxwu/amsd69n_getrigidity_dt_p000-009/1305853512_10.root 的 measured branches。
未读取 MC 文件、真值 branch 或任何真值评价结果。当前 producer 与实测文件 branch shapes 一致，但未取得该历史 flight production 对应 binary/source hash，不能把当前 checkout 当作已核实的完整生产 provenance。

## 主要结论

RICH beta 可以提供来自另一探测器的速度约束，基础 RICH quality flags 没有直接用 TOF beta。现成 SelectRich() 另外施加 TOF–RICH agreement cut，训练和独立验证不能原样复用。RICH 是下方 radiator 处的速度，不能无修正地作为 AMS 中心或 Tracker L5 处的速度标签。低 beta 的能损差和不完整接受度必须保留为限制。

特别重要：producer 在无环时没有清空全部 RICH quality fields；必须先检查有效 ring 和 track association。样本 7883 events 中，5358 个没有有效匹配 RICH 的事件仍带有 rich_good && rich_clean=true，不能只靠这两个 bool 判断。

## Branch 定义与形状

来源 /nas01/home/hxwu/AMSSoft/amsd69n/analysis_amsd69n.C:663–700、1238–1265、4371–4409。

| Branch | 形状 / 单位 | 实际定义 |
|---|---|---|
| rich_beta | float[3]，无量纲 | [0]=ring->getBeta()：BetaRefit × betaCorrection；[1]=ring->Beta：只使用 hit 的原始重建；[2]=BetaRefit：按每 hit 光电子数重加权的重建 |
| rich_ebeta | scalar float，beta 的绝对误差 | ring->ErrorBeta，旧重建误差估计 |
| rich_ebeta1 | scalar float，beta 的绝对误差 | ring->getBetaError()；AgL sqrt((0.00302035)^2/Used+(0.000247571)^2)，NaF 0.00774510/sqrt(Used) |
| rich_good / rich_clean | scalar bool | good = clean && uncorrected expected PE > 1e-6；clean = (Status&1)==0 |
| rich_NaF | scalar bool | (Status&2)!=0 |
| rich_hit | scalar int | ring->getHits()，原始 Used |
| rich_used | scalar int | ring->getUsedHits()，默认应用 PMT correction 后的 hit 数，可能比 rich_hit 小 |
| rich_pmt | scalar int | ring->getPMTs()：环中不同 PMT 的数量，跳过被带电粒子穿过的 hit |
| rich_stat | scalar int | RichRingR::Status bit field，不可整体与0比较；bit0=清理重建标志，bit1=NaF，另含关联和几何信息 |
| rich_pb | scalar float，概率 | ring->getProb()，环电荷方位分布的 Kolmogorov quality probability |
| rich_BetaConsistency | scalar float | abs(RichRingR::Beta - matched RichRingBR::Beta)，比较两种 RICH 算法；无标准 tracker 匹配解返回 -1 |
| rich_itrtrack | scalar int | ring->iTrTrack()，必须与当前 itrtrack 对应 |
| rich_pos | float[3]，cm | getTrackEmissionPoint() 返回 AMSTrPars[0:3]，即重建所用 radiator 处 track 坐标 |
| rich_theta | scalar float，rad | AMSTrPars[3]，track theta |
| rich_n | float[2] | [0]=getIndexUsed()，含 betaCorrection 的折射率；[1]=getRawIndexUsed() |
| rich_npe | float[3]，光电子数 | [0]=环收集 PE，默认应用 correction；[1]=该轨迹/速度下 Z=1 环期望 PE，默认 correction；[2]=事件全部 RICH hit 的 PE 和，排除被带电粒子穿过的 PMT |
| rich_PMTChargeConsistency | scalar float | hit 与 PMT charge consistency test |
| rich_cstat | scalar int | PMT corrections flag：-1 未做，0 完成，1 失败；不能凭字段名把0判失败 |

rich_npe[2] 不是“包含一切 PMT 的原始总 PE”。root.h:1907–1908 和 CC/root.C:8470–8478 明确排除了 crossed PMT。现成 cut 使用 rich_npe[0]/rich_npe[2]；应先要求分母有限且>0，不额外强制比例≤1。

library 定义来源：
/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/include/root.h:2040–2200；
同目录 CC/root.C:8860–8983、8997–9015、9108–9174；
同目录 CC/richrec.C:1031–1048。

getBeta() 默认 correctBias=false；getBeta 的索引/位置均匀性校正并不等于近阈值残余偏差已消失。CC/root.C:9108–9174 的专门 CalibrationBias 分支默认没有开启。其每事件校正没有直接使用 TOF 时间：CC/RichBeta.C:52、195 使用位置/方向，:381 的扩展项可取 Tracker charge；GetInnerQYJ 默认 beta=1。共享 Tracker 几何仍会使 TOF 路长误差和 RICH beta 误差相关，不能据“不同探测器”断言条件无偏噪声。

## 建议固定的独立参考质量选择

来自 /nas01/home/hxwu/AMSSoft/TofCalib/include/parameter.h:149–167，与 amsd69n/SelEvent14.h:4667–4705 的历史 RICH cut 配对。保留其 RICH 自身条件，移除 get_betadiff()，增加有效性下界：

- irich>=0，itrtrack>=0，rich_itrtrack==itrtrack；所有参与运算的 float 有限；rich_beta[0]>0。
- rich_good && rich_clean；rich_hit>=3，rich_used>=3，rich_pmt>=3，rich_pb>=0.01。
- rich_npe[2]>0 且 rich_npe[0]/rich_npe[2]>=0.4；0<=rich_PMTChargeConsistency<=15。
- 0<=rich_BetaConsistency<=0.01 (NaF) 或 <=0.005 (AgL)。下界是必需的，-1 表示没有可比较解。
- rich_npe[1]>=1 (NaF) 或 >=2 (AgL)；有限且>0 的 rich_ebeta1。可额外排除 rich_cstat==1，并单独报告 -1。
- NaF：max(abs(x),abs(y))<=17 cm。
- AgL：max(abs(x),abs(y))>=19 cm，x*x+y*y<=3500 cm²，rich_tile 不属于 {3,7,12,20,87,100,108}。bad-tile veto **仅适用于 AgL**。
- 不加任何 TOF–RICH beta difference cut，不用模型新 beta 做 selection。不要强制 beta_RICH<=1，否则会截断测量噪声并偏移高能参考。

上述数值是已有 calibration selection 的可追踪起点，不是已证明在本任务所有 charge 都最优的门槛。若更改门槛，用 flight holdout 稳定性决定并冻结，不按 MC 表现选择。

SelectRich() 的循环 selection 源码：
amsd69n/SelEvent14.h:599 定义 get_betadiff()=(tof_betah-rich_beta[0])/rich_beta[0]；
:4681 要求其绝对值<=0.06。RichRingR::getBetaConsistency() 本身不做这个比较；它只接受 RichRingBR Status%10==2 的标准 tracker 解，排除 TOF-track fallback，见 root.h:2350–2368 与 root.C:8983–8994。

## 几何与速度参考面

producer analysis_amsd69n.C:97 的 TOFZ={65.2,62.1,-62.1,-65.2} cm 明确标注 central cm；
:106 的 Tracker nominal z 为
{159.04,53.05,29.22,25.24,1.706,-2.292,-25.26,-29.25,-136.03} cm。

这些是标称层中心，**不是完整材料边界或 bar 包络**，不能据此构造材料分区而省略前/中/后定义。

本次唯一 flight 文件中：
- L5 track z 中位数 +1.7101 cm，L6 -2.2866 cm；z=0 的 AMS 中心与 L5 处并非同一点。
- AgL rich_pos z 中位数 -74.7050 cm。
- NaF rich_pos z 中位数 -75.5516 cm。
- 向下事件从 L5 / z=0 至 RICH 仍经过 L6–L8、下方两层 TOF 及相关支撑/被动材料。即使 RICH beta 本身准确，也通常低于上游中心 beta。

若目标是 u=1/beta_center，直接训练贴合 1/beta_RICH 会混入正的能损位移。需要明确 transport model / 可控高能近似 / 参考平面，而不能把 flight RICH 标签称为 center truth。AgL n约1.05，Cherenkov 阈值 beta约0.95；NaF n约1.334，阈值约0.75，且覆盖中心较小几何区域。低于NaF阈值没有这种实测标签，AgL/NaF 有标签子样本不能证明整个 acceptance 的全面提升。

## Rigidity

producer:450、999：tk_rigidity1 float[3][3][7]；
:2802–2917 定义 axes：
算法 [0] Choutko、[1] GBL、[2] Kalman；
alignment [0] PG、[1] PG+CIEMAT、[2] V6；
span [0] maximal、[1] inner、[2] inner+L1、[3] inner+L9、[4] full span、[5] inner up、[6] inner down。

tk_rigidity 为 scalar float，由 Fill():4654–4658 的 trk->GetRigidity() 填充，使用默认 fit。rigidity 单位 GV，tk_erigidity1 单位 GV^-1 且代表 sigma(1/R)：
AMSVX_TOFBetaFit_el9/include/TrTrack.h:105–107、625–630。

若需要可审计高 rigidity beta≈1 anchor，可固定例如 tk_rigidity1[0,2,1]，配套 tk_chis1[0,2,1,0:2] 与 tk_erigidity1[0,1]；后者只存 V6 的误差，不能直接假定它与 raw tk_rigidity 默认 fit 相配。
GetRigidity 默认 iz=-1 并非显式 center plane；API iz=1 对应 z=0，iz=2 对应 z=-70（仍不同于真实 radiator z），且源码只有特定 Kalman/DAF fit 会使用不同 PZRigidity。

Tracker charge reconstruction 使用 obeta=TOF beta 做能损修正；随后 mass≈(MHe/2)*Z 用于重拟合，因此 rigidity 不是完全没有 TOF 间接输入，但拟合没有直接传入 TOF 速度标签。高 rigidity 饱和参考受同位素 A/Z 和误重建/尾部影响，低 rigidity 直接换算 beta 还需要质量假设；不能当普适独立无误差真值。

## 小样本可用性

唯一文件共 7883 events，有效匹配 RICH 2380。按上述主要严格 cuts（含 geometry、ratio、consistency、error，未额外 cstat/tile-border margin）为874：
Z1=591、Z2=253、Z6=9、Z8=6。
该文件 Z6 与 Z8 的严格样本均为 AgL，beta范围分别0.98209–0.99888与0.98180–0.99780。
这只是容量与branch sanity check，没有做模型比较，也不能据一个文件估计整个 flight reference coverage。后续应固定并报告全样本 Z × beta × radiator × run 的数量与 coverage。

