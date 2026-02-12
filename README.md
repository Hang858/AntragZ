# ANTRAG_Z: 基于整数运算的格基签名算法参考实现

## 1. 项目简介

**ANTRAG_Z** 是一种基于 Falcon 框架改进的抗量子密码数字签名算法。
本项目是该算法的 **C 语言参考实现**。

### 核心特性

* **签名验签过程完全无浮点**：移除了原 Falcon 算法中的所有浮点数运算，仅使用整数和定点数逻辑，但在密钥生成阶段仍有浮点运算。
* **具体算法逻辑参考算法设计文档**
---

## 2. 快速开始

本项目使用 `CMake` 进行构建，使用如下命令可输出调试信息。
```bash
CMake -CMAKE_BUILD_TYPE=Debug ..
```

### 编译

```bash
mkdir build
cd build
cmake ..
make

```

### 运行测试

编译完成后，可以运行测试程序以验证算法的正确性：

```bash
# 运行签名与验签的完整流程测试
./test/test_sign_impl

```

---

## 3. 代码结构导航

### 3.1 基础数学库 (`src/math/`)


* **`ntt.c` / `fft.c**`: 数论变换（NTT）实现。这是硬件加速的重点，用于加速多项式乘法。ntt.c 是搬运Falcon原有的算法中的ntt，主要用于NTRUSolve，目前项目中其他模块还未使用ntt进行优化，先用朴素方法实现的正确性版本。
* **`poly.c`**: 多项式加减乘除的基本操作，此文件中的多项式乘法均为朴素乘法，尚未优化。
* **`zint.c`**: 大整数运算逻辑。用于NTRUSolve，原来Falcon项目中搬运来的，其他模块未使用。

### 3.2 密钥生成 (`src/keygen.c`, `src/ntru/`, `src/prematrix/`)

对应设计文档中的 。

* **逻辑流**：`keygen.c` -> `ntru/keygen_ntru.c` (求解F,G) -> `prematrix/prematrix.c`。
* **硬件注意**：
* `src/prematrix/` 包含了 **MIGD** 和 **EIGD** 算法（`src/sampler/eigd.c`），用于生成预计算矩阵 。
* 这部分计算量较大，但通常只在生成密钥时运行一次。



### 3.3 签名生成 (`src/sign.c`, `src/sampler/`)

对应设计文档中的 ，分为**离线（Offline）**和**在线（Online）**两个阶段。

* **`sign.c`**: 顶层逻辑。
* `crypto_sign_signature()`: 包含完整的签名流程。


* **`src/sampler/zalcon_samp.c`**: **核心采样器**。
* 这是无浮点高斯采样的核心实现，替代了 Falcon 的 Fast Fourier Sampling，具体方法基于Zalcon。


* **`src/sampler/sample.c`**: 基础采样函数。

### 3.4 验签 (`src/verify.c`)

对应设计文档中的 。

<<<<<<< HEAD
* 逻辑相对简单，主要包含 Hash 计算、解压缩（Decompress）和范数检查， Hash计算，解压缩，与Falcon算法保持一致。
=======
* 逻辑相对简单，主要包含 Hash 计算、解压缩（Decompress）和范数检查。 Hash计算，解压缩，与Falcon算法保持一致。
>>>>>>> ef9671a (增加速度测试程序，验签替换ntt算子)

### 3.5 辅助工具 (`include/`)

* **`CDT_gen.py` / `calc_lw_params.py**`: Python 脚本。
* 用于生成采样所需的 **CDT 表（Cumulative Distribution Table）** 和其他常量参数。

---

## 5. 联系方式

* **算法实现**: 张航
<<<<<<< HEAD
* **邮箱**: zhanghang2024@iie.ac.cn
=======
* **邮箱**: zhanghang2024@iie.ac.cn
>>>>>>> ef9671a (增加速度测试程序，验签替换ntt算子)
