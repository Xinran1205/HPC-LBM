# 格子玻尔兹曼方法 (LBM) 优化项目

## 项目介绍
格子玻尔兹曼方法 (LBM) 是一种功能强大且灵活的方法，可用于模拟复杂的流体流动。然而，要充分利用其功能，优化底层代码以提高性能至关重要。本项目详细介绍了我的全面优化过程，该过程涵盖几个阶段：

## 优化阶段

### 串行优化
我首先从串行优化技术入手，为性能提升打下坚实的基础。我实施了循环融合和指针交换等技术来简化执行流程并降低开销。

### 向量化优化
在优化过程中，我深入分析了内存访问模式对整体性能的影响。为了提升数据访问速度并促进编译器进行有效的向量化，我首先对代码中的数据结构进行了调整，确保它们更适合向量操作。接着，实施内存对齐技术，这优化了数据加载过程。为了进一步提升向量化效率，我尝试使用了多种编译器及其不同版本，针对每种情况应用特定的编译器标志（flags）。通过详细分析各编译器提供的优化报告（optimization reports），我能够选择出最适合项目的编译策略，从而最大限度地提高了代码的执行效率和CPU并行处理能力。

### 使用 OpenMP 实现并行加速
我的优化工作最终利用了 OpenMP（一种并行编程模型）来加速计算。通过将代码从单核扩展到 28 核，我实现了性能和效率的显著提升，使我的 LBM 模拟更适用于复杂场景。

### 使用 MPI 实现分布式内存并行
本部分通过使用消息传递接口 (MPI) 扩展了 LBM 代码的优化。优化后的代码在 BlueCrystal 超级计算机上运行，利用了四个节点，每个节点配备 28 个内核。此部分基于之前优化的串行代码。与 OpenMP 中使用的共享内存模型不同，MPI 采用适用于多节点计算环境的分布式内存模型。每个节点都有自己独立的物理内存，数据交换通过网络或其他通信方法明确进行。我为每个进程实施了数据的分配和初始化，并且采取了负载平衡策略。此外，通过 Halo Exchange 策略，我优化了节点间的数据通信，在确保数据正确性的同时优化了内存存储。

## Results
<img src="/pic/1.png" alt="scalability" width="450" height="350">
<img src="/pic/2.png" alt="speedup" width="500" height="350">
