<div align="center">
  <h1>AMS飞行时间探测器的粒子速度非线性重建方法研究</h1>
  
  [![English](https://badgen.net/badge/Language/English/blue?icon=github)](README_EN.md) [![简体中文](https://badgen.net/badge/语言/简体中文/red?icon=github)](README.md)
</div>

## 项目介绍

本项目旨在研究阿尔法磁谱仪（AMS）飞行时间探测器（TOF）中粒子速度的非线性重建方法。AMS是一个在国际空间站上运行的粒子物理实验，TOF探测器用于测量带电粒子的飞行方向和速度。然而，由于带电粒子在TOF材料中发生电离能量损失，导致粒子减速。如果使用线性函数拟合粒子的时间-空间关系，会导致粒子速度重建存在一定偏差，特别是在低能时更加明显。

本研究通过在粒子速度拟合过程中引入粒子能量损失项，进行非线性拟合，以减少速度重建的偏差并提高精度。这将有助于AMS更准确地鉴别核同位素。

## 功能特性

- 实现粒子运动模拟与能量损失计算
- 提供线性和非线性粒子速度重建方法
- 磁场数据可视化与分析
- 能量损失校正模型
- Beta（速度/光速）比较分析工具

## 软件架构

本项目基于ROOT框架开发，主要包含以下组件：

- `ParticleData`：粒子数据结构定义
- `ParticlePropagator`：粒子传播算法实现
- `BetaFitter`：粒子速度拟合核心类
- `macro`：数据分析和可视化脚本

## 安装与使用

### 系统要求
- Linux操作系统
- ROOT 5.34+（支持ROOT 6.x）
- 编译器支持C++11标准（gcc或icpx）

### 编译步骤

1. 标准编译：
```bash
./build.csh
```

2. 以调试模式编译：
```bash
./build.csh debug
```

3. 清理编译文件：
```bash
./build.csh clean
```

### 运行示例

```bash
./run.csh [输入文件] [输出文件] [参数]
```

或使用本地处理脚本：
```bash
./run_local.sh
```

## 数据分析

项目提供了多种数据分析和可视化工具：

### Beta比较分析

```bash
root -l 'macro/plotBetaComparison.C("test.root", "beta_comparison.pdf")'
```

### 能量损失分析

```bash
root -l 'macro/plotEnergyLoss.C("test_el.root", "energy_loss.pdf")'
```

### 磁场分析

```bash
root -l 'macro/plotMagneticField.C(0.0, 0.0, "test.root", "magnetic_field.pdf")'
```

## 结果分析

项目通过比较线性和非线性速度重建方法的残差分析，证明了非线性方法在β较小的情况下具有明显的精度提升。研究表明带电粒子在探测器材料中的能量损失会对速度测量产生系统性偏差，而通过引入能量损失校正项可以显著改善测量精度。