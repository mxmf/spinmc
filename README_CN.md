# SpinMC

SpinMC 是一个开源的经典自旋模型蒙特卡洛模拟程序包，专为研究磁性系统和统计物理中的相变而设计。

## ✨ 功能特点

- **多种自旋模型**：
  - Ising 模型
  - XY 模型
  - Heisenberg 模型

- **支持的相互作用**：
  - 交换耦合（各向同性）
  - 单离子各向异性能
  - [计划支持] 各向异性交换耦合
  - [计划支持] 外磁场

- **支持的算法**：
  - Metropolis
  - Wolff 簇算法
  - [计划支持] 并行回火（副本交换）

- **模拟能力**：
  - 能量、热容
  - 磁化强度、磁化率
  - 绝对磁化强度、绝对磁化率
  - 分组（子晶格）磁化强度与磁化率
  - [计划支持] 磁滞回线
  - [计划支持] Binder 累积量 (u4)

## 📦 安装

### 环境要求

- Python 3.9+（pip 安装时需要）
- Rust 1.88+（从源码编译时需要）

### 从 PyPI 安装

```bash
pip install spinmc
```

### 从预构建版本下载

1. 访问 Releases 页面
2. 下载适合你操作系统的包
3. 解压缩并运行可执行文件

### 从源码安装

```bash
git clone https://github.com/mxmf/spinmc.git
cd spinmc
cargo build --release
```

## 🚀 快速开始

1. 创建一个配置文件（例如`ising.toml`)

```toml
[grid]
dimensions = [50, 50, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, false]

[simulation]
initial_state = "random"
boltzmann_constant = 1
model = "ising"
equilibration_steps = 10000
measurement_steps = 100000
algorithm = "wolff"
num_threads = 10
[simulation.temperature_range]
start = 1
end = 3
step = 0.1

[output]
outfile = "result.txt"
energy = true
heat_capacity = true
magnetization = true
susceptibility = true
magnetization_abs = true
susceptibility_abs = true

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[0, -1, 0], [0, 1, 0], [-1, 0, 0], [1, 0, 0]]
strength = 1.0
```

2. 运行模拟

```bash
spinmc -i ising.toml
```

3. 模拟结果将保存到 `result.txt`, 包含你在配置文件中选择的观测量

4. 如果是通过 `Python` 安装 的 spinmc，可以通过如下命令绘图查看：

```bash
spinmc plot -i result.txt
```

## 📚 更多示例请查看[示例文件夹](examples)
