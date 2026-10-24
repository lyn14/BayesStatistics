# BayesStatistics

贝叶斯统计课程的学习与作业代码，按作业编号整理，包含共轭推断、后验预测和离散网格抽样等示例。

## 内容

| 目录 | 主要内容 |
| --- | --- |
| `hw2/` | Poisson–Gamma 共轭推断、负二项后验预测、数值积分验证，以及课程 PDF |
| `hw3/` | 通过离散网格抽样估计标准正态与均匀分布的前六阶原点矩，比较区间与网格点数的影响；附 PDF 和运行日志 |
| `hw4/` | 大鼠肿瘤 Beta–Binomial 层次模型：vague prior 与均匀先验下的网格计算及等高线图 |

## 大鼠肿瘤数据

`hw4/Hierarchical_Models_Vague.R` 中的数据来自 Gelman 等人的《Bayesian Data Analysis》第三版（BDA3）第 5 章表 5.1 的大鼠肿瘤示例，历史实验数据源自 Tarone（1982）。教材示例包含 70 组历史实验及一组新增实验，用于说明二项分布的层次贝叶斯建模。

- `rat.y[j]`：第 `j` 组实验中出现肿瘤的大鼠数量。
- `rat.n[j]`：第 `j` 组实验的大鼠总数。

两者是观测数据；`y[j] / n[j]` 是该组观察到的肿瘤发生比例。常用模型为 `y[j] ~ Binomial(n[j], theta[j])`，其中 `theta[j]` 才是待推断的该组肿瘤发生概率。

参考：[数据说明与文献出处](https://paulnorthrop.github.io/bang/reference/rat.html)；Tarone, R. E. (1982). *The use of historical information in testing for a trend in proportions*. Biometrics, 38, 215–220. [DOI](https://doi.org/10.2307/2530304)。

## hw4：大鼠肿瘤层次模型

使用上述 70 组历史实验数据，建立 `y[j] ~ Binomial(n[j], theta[j])`、`theta[j] ~ Beta(a, b)` 的层次模型。积分消去各组的 `theta[j]` 后，在 `u = log(a/b)`、`v = log(a+b)` 坐标下进行网格计算并绘制等高线。

- [Hierarchical_Models_Vague.R](hw4/Hierarchical_Models_Vague.R)：采用 BDA3 的 vague hyperprior，`p(a,b) ∝ (a+b)^(-5/2)`，并加入坐标变换的 Jacobian 项 `log(a) + log(b)`。
- [Hierarchical_Models_Uniform.R](hw4/Hierarchical_Models_Uniform.R)：采用独立均匀先验 `a ~ Uniform(0,40)`、`b ~ Uniform(0,100)`。
- 结果图：[vague prior](hw4/RatTumor_Vague.png)、[均匀先验](hw4/RatTumor_Uniform.png)。

等高线数值为相对于网格最大值的密度高度，未归一化为概率，也不代表可信区域的概率。当前均匀先验脚本未加入 Jacobian 项，绘制的是 `(a,b)` 后验密度在 `(u,v)` 网格上的取值；若要与 vague prior 图作为 `(u,v)` 后验密度比较，需先补入 `log(a) + log(b)`。

## 运行

Python 脚本建议使用 Python 3.9 或更高版本。在仓库根目录执行：

```bash
python -m pip install scipy sympy
python "hw2/Poisson-Gamma posterior.py"
python hw3/Normal_grid_sampling.py
python hw3/grid_moments.py
```

`hw3` 的两个脚本仅依赖 Python 标准库。网格抽样示例使用固定随机种子 `2026`，输出矩估计及其与理论值的误差。

hw4 使用 R 基础函数，无需额外安装 R 包。在 R / RStudio 中打开脚本运行，或在仓库根目录执行：

```bash
Rscript hw4/Hierarchical_Models_Vague.R
Rscript hw4/Hierarchical_Models_Uniform.R
```

交互式运行时显示等高线图；通过 `Rscript` 运行时，默认图形设备通常将图写入当前目录的 `Rplots.pdf`，连续运行可能覆盖该文件。脚本不会自动更新仓库中已有的 PNG 图。

本仓库用于课程学习与数值实验，代码随学习进度更新。
