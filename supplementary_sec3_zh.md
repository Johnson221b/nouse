# 补充材料 3：空间耦合指标 S、C 的推导与证明

本节给出正文 §1 中两个空间耦合指标的严格数学推导、上下界证明，以及完全随机情形下的期望值，并附蒙特卡洛数值验证。

## 3.1 定义与记号

设研究区由 $N$ 个空间单元（CBG 或出租车分区）组成，其地理邻接关系构成图 $G$。定义行标准化 Queen 邻接权重矩阵 $W \in \mathbb{R}^{N \times N}$：

$$
W_{ij} =
\begin{cases}
\frac{1}{k_i}, & j \in \mathcal{N}(i),\\[4pt]
0, & \text{否则},
\end{cases}
$$

其中 $\mathcal{N}(i)$ 为单元 $i$ 的地理邻接集合，$k_i = |\mathcal{N}(i)|$ 为单元 $i$ 的度（邻接单元数）。由此 $W$ 行随机，即 $W\mathbf{1} = \mathbf{1}$，但一般**不**双随机。

设 $F \in \mathbb{R}_{\ge 0}^{N \times N}$ 为逐日 OD 流矩阵（$F_{ij}$ 为从 $i$ 到 $j$ 的流量）。将其归一化为流分布

$$
F_n = \frac{F}{\sum_{i,j} F_{ij}}, \qquad \sum_{i,j} (F_n)_{ij} = 1,
$$

并中心化（去除全局均匀基线）

$$
F_c = F_n - \tfrac{1}{N^2}\mathbf{1}\mathbf{1}^\top, \qquad \sum_{i,j} (F_c)_{ij} = 0 .
$$

记流入与流出占比向量

$$
p_{\text{out}} = F_n \mathbf{1}, \qquad p_{\text{in}} = F_n^\top \mathbf{1},
$$

其中 $(p_{\text{out}})_i$ 为从单元 $i$ 出发的流量占比，$(p_{\text{in}})_i$ 为到达单元 $i$ 的流量占比。对任意向量 $v$，定义 z-score 标准化 $z(v) = \frac{v - \bar{v}}{s}$，其中 $\bar{v} = \frac{1}{N}\sum_i v_i$，$s^2 = \frac{1}{N}\sum_i (v_i - \bar{v})^2$；于是 $\|z(v)\|_2^2 = N$。

正文中的两个空间耦合指标（即逐日流矩阵的 Moran 型指标）定义为：

$$
S = \frac{\|W F_c\|_F^2}{\|F_c\|_F^2} \qquad \text{（流级结构稳定性）},
\qquad
C = \frac{z(p_{\text{in}})^\top W\, z(p_{\text{out}})}{N} \qquad \text{（交叉流向耦合强度）},
$$

其中 $\|\cdot\|_F$ 为 Frobenius 范数。

## 3.2 构造动机

**S 的构造。** 记 $v = \text{vec}(F_c) \in \mathbb{R}^{N^2}$ 为 $F_c$ 的按列堆叠向量。由 Kronecker 乘积恒等式 $\text{vec}(W F_c) = (I_N \otimes W)\, v$，S 可改写为 Rayleigh 商

$$
S = \frac{v^\top (I_N \otimes W^\top W)\, v}{v^\top v}.
$$

矩阵 $W^\top W$ 半正定，故 $I_N \otimes W^\top W$ 半正定。$(W F_c)$ 的第 $i$ 行是单元 $i$ 各邻居流出占比的算术平均，因此 $\|W F_c\|_F^2$ 度量的是「经过一步地理平滑后的流矩阵能量」：若相邻单元的流出剖面彼此相似，平滑几乎不损失能量，S 接近其上界；若流在空间上完全无序，平滑会大幅抵消能量，S 接近其下界。因此 S 本质上是流矩阵在其空间低频子空间（沿图平滑方向）中的能量占比，即流级空间自相关强度。

**C 的构造。** C 是流入场 $p_{\text{in}}$ 与流出场 $p_{\text{out}}$ 之间的双变量 Moran's I：分子 $z(p_{\text{in}})^\top W z(p_{\text{out}})$ 度量「高流入单元的邻居是否同时高流出」。当流入与流出在空间上同向耦合（高流入区邻近高流出区）时 C 为正，反向耦合时为负，彼此独立时趋近于 0。

## 3.3 上下界

**命题 1（S 的上下界）。** 设 $\sigma_{\max}(W)$ 为 $W$ 的最大奇异值，则

$$
0 \le S \le \sigma_{\max}^2(W).
$$

*证明。* 因 $W^\top W$ 半正定，$I_N \otimes W^\top W$ 半正定，Rayleigh 商非负，故 $S \ge 0$。由 Courant–Fischer 极大极小定理，

$$
S \le \lambda_{\max}\!\big(I_N \otimes W^\top W\big) = \lambda_{\max}(W^\top W) = \sigma_{\max}^2(W),
$$

其中用到 $\lambda_{\max}(I_N \otimes M) = \lambda_{\max}(M)$ 与 $\lambda_{\max}(W^\top W) = \sigma_{\max}^2(W)$。$\blacksquare$

**推论 1。** 若 $W$ 双随机（例如正则图且权重对称化），则 $\sigma_{\max}(W) = 1$，故 $S \in [0,1]$。实际采用 Queen 邻接（仅行标准化、非列标准化），一般 $\sigma_{\max}(W) \ge 1$，故按一般情形 $0 \le S \le \sigma_{\max}^2(W)$ 表述。

*说明。* 对任意行随机 $W$，因 $W\mathbf{1}=\mathbf{1}$ 有 $\|W\mathbf{1}\|_2 = \sqrt{N} = \|\mathbf{1}\|_2$，故 $\sigma_{\max}(W) \ge 1$；等号当且仅当 $\mathbf{1}$ 同时为 $W^\top$ 的（近似）主导特征向量，即双随机情形。

**命题 2（C 的上下界）。**

$$
|C| \le \sigma_{\max}(W).
$$

*证明。* 由 Cauchy–Schwarz 不等式与谱范数定义，

$$
|C| = \frac{1}{N}\big|z(p_{\text{in}})^\top W z(p_{\text{out}})\big|
\le \frac{1}{N}\,\|z(p_{\text{in}})\|_2\,\|W z(p_{\text{out}})\|_2
\le \frac{1}{N}\,\|z(p_{\text{in}})\|_2\,\sigma_{\max}(W)\,\|z(p_{\text{out}})\|_2 .
$$

由 $\|z(p_{\text{in}})\|_2 = \|z(p_{\text{out}})\|_2 = \sqrt{N}$，得

$$
|C| \le \frac{1}{N}\sqrt{N}\,\sigma_{\max}(W)\,\sqrt{N} = \sigma_{\max}(W). \qquad \blacksquare
$$

**推论 2。** 若 $W$ 双随机，则 $|C| \le 1$，即 $C \in [-1,1]$。

## 3.4 完全随机情形的期望值

**命题 3（S 的零模型期望值）。** 在零假设（$F_c$ 各元素独立同分布、无空间结构）下，

$$
\mathbb{E}[S] = \frac{N\sum_{i=1}^N \frac{1}{k_i} - 1}{N^2 - 1} \;=\; \frac{1}{N}\sum_{i=1}^N \frac{1}{k_i} + O(N^{-2}).
$$

特别地，当图为 $k$-正则图（所有 $k_i = k$）时，

$$
\mathbb{E}[S] \;\xrightarrow{N\to\infty}\; \frac{1}{k}.
$$

*证明。* 记 $v = \text{vec}(F_c)$、$A = I_N \otimes W^\top W$，则 $S = \frac{v^\top A v}{v^\top v}$。中心化约束 $\sum_{ij}(F_c)_{ij}=0$ 使 $v$ 落在与 $\mathbf{1}_{N^2}$ 正交的 $d = N^2 - 1$ 维子空间 $V$ 上。在零假设下 $v$ 在该子空间内各向同性（方向均匀），利用单位球面上均匀向量 $u$ 的恒等式 $\mathbb{E}[uu^\top] = \frac{1}{d}P_V$（$P_V$ 为到 $V$ 的正交投影），得

$$
\mathbb{E}[S] = \frac{\text{tr}(P_V A)}{N^2 - 1}.
$$

记到 $\mathbf{1}_{N^2}$ 方向的投影 $P_\perp = \frac{1}{N^2}\mathbf{1}_{N^2}\mathbf{1}_{N^2}^\top$，则 $\text{tr}(P_V A) = \text{tr}(A) - \text{tr}(P_\perp A)$。因 $W$ 行随机，$W\mathbf{1}=\mathbf{1}$，故

$$
\text{tr}(P_\perp A) = \frac{1}{N^2}\,\mathbf{1}_{N^2}^\top A\, \mathbf{1}_{N^2}
= \frac{1}{N^2}\sum_{i=1}^{N} \mathbf{1}^\top W^\top W \mathbf{1}
= \frac{1}{N^2}\cdot N\,\|W\mathbf{1}\|_2^2 = 1 .
$$

另一方面，对行标准化权重矩阵 $W$，

$$
\text{tr}(W^\top W) = \sum_{i,j} W_{ij}^2 = \sum_{i=1}^N \sum_{j\in\mathcal{N}(i)} \frac{1}{k_i^2} = \sum_{i=1}^N k_i \cdot \frac{1}{k_i^2} = \sum_{i=1}^N \frac{1}{k_i},
$$

从而 $\text{tr}(A) = N\,\text{tr}(W^\top W) = N\sum_i \frac{1}{k_i}$。代入得

$$
\mathbb{E}[S] = \frac{N\sum_i \frac{1}{k_i} - 1}{N^2 - 1},
$$

其主项为 $\frac{1}{N}\sum_i \frac{1}{k_i}$，误差 $O(N^{-2})$。$k$-正则图时主项退化为 $\frac{1}{k}$。$\blacksquare$

**命题 4（C 的零模型期望值）。** 在零假设下，流入场与流出场彼此独立，故 $\mathbb{E}[z(p_{\text{in}})] = 0$ 与 $\mathbb{E}[z(p_{\text{out}})] = 0$，进而

$$
\mathbb{E}[C] = \frac{\mathbb{E}\big[z(p_{\text{in}})^\top W z(p_{\text{out}})\big]}{N} = 0 .
$$

*说明。* 该结果意味着：无空间结构时，S 并不趋于 0，而是趋于其几何期望 $\frac{1}{N}\sum_i \frac{1}{k_i}$（正则图 $1/k$）；C 则趋于 0。因此正文以 S 的偏离度（相对零模型）衡量结构失序、以 C 的正负衡量流入–流出耦合，均有明确的零模型基准。

## 3.5 蒙特卡洛验证

图 S7 以蒙特卡洛模拟验证上述理论值。在环图（$k=2$，理论 $\mathbb{E}[S]\to 1/2$）上，随节点数 $N$ 增大，模拟 S 的均值（200 次独立随机流矩阵）逐步收敛到虚线所示的 $1/2$，与命题 3 一致；同时 C 的均值收敛到 0，与命题 4 一致，验证了零模型基准的正确性。
