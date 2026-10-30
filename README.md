# BayesStatistics

贝叶斯统计课程的学习与作业代码，按作业编号整理，包含共轭推断、层次模型、后验预测和离散网格抽样等示例。

## 内容

| 目录 | 主要内容 |
| --- | --- |
| `hw2/` | Poisson–Gamma 共轭推断、负二项后验预测、数值积分验证，以及课程 PDF |
| `hw3/` | 通过离散网格抽样估计标准正态与均匀分布的前六阶原点矩，比较区间与网格点数的影响；附 PDF 和运行日志 |
| `hw4/` | 大鼠肿瘤与 BDA3 Exercise 5.13 自行车比例的 Beta–Binomial 层次模型：先验比较、网格抽样、可信区间及收缩估计 |
| `hw5/` | 新组自行车数量的后验预测，以及交通总量的 Poisson–Gamma 层次模型、超先验比较与联合后验抽样 |

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

等高线数值为相对于网格最大值的密度高度，未归一化为概率，也不代表可信区域的概率。两份大鼠脚本均已加入坐标变换的 Jacobian 项 `log(a) + log(b)`。

## hw4：BDA3 Exercise 5.13 — 自行车比例

使用 10 组车辆计数数据，`bikes.y[j]` 为第 `j` 组自行车数量，`bikes.n[j]` 为该组车辆总数，`theta[j]` 为自行车所占比例。沿用 Beta–Binomial 层次模型，比较 vague hyperprior 与独立均匀超先验对后验推断的影响。

| 程序 | 内容 | 已保存的结果图 |
| --- | --- | --- |
| [Bikes_Vague.R](hw4/Bikes_Vague.R) | vague prior 下的超参数后验等高线、网格抽样及第 5 组比例的 95% 等尾可信区间 | [等高线](hw4/Bikes_Vague.png)、[区间与后验直方图](hw4/Interval_Vague.png) |
| [Bikes_Uniform.R](hw4/Bikes_Uniform.R) | 均匀先验下的对应计算 | [等高线](hw4/Bikes_Uniform.png)、[区间与后验直方图](hw4/Interval_Uniform.png) |
| [Bikes_Expectation.R](hw4/Bikes_Expectation.R) | 通过 1000 次超参数网格抽样，平均条件期望 `(a+y[j])/(a+b+n[j])`，估计各组后验均值并与 MLE `y[j]/n[j]` 比较；随机种子为 `123` | [后验均值与 MLE](hw4/Means.png) |


## hw5：新组预测与层次 Poisson 模型

延续 hw4 的自行车数据分析，预测新路口在 100 辆车中出现的自行车数量；另用 [Traffic_total.txt](hw5/Traffic_total.txt) 中的 10 组交通总量建立 Poisson–Gamma 层次模型。详细推导和已有数值结果见 [hw5/README.md](hw5/README.md)。

| 程序 | 内容 | 已保存的结果图 |
| --- | --- | --- |
| [Bikes_new_sample.R](hw5/Bikes_new_sample.R) | Beta–Binomial 模型的新组后验预测：抽取超参数、新组比例，再生成计数；与条件 Beta–Binomial 概率的蒙特卡洛平均比较，并计算 95% 预测区间 | [新组预测分布](hw5/posterior_predictive_new_group.png) |
| [Bikes_poisson.R](hw5/Bikes_poisson.R) | 平坦超先验下的交通计数模型，计算超参数边际后验等高线并进行 10,000 次网格抽样 | [等高线](hw5/hyperparameter_contours.png)、[散点图](hw5/hyperparameter_scatter.png) |
| [Bikes_poisson_vague.R](hw5/Bikes_poisson_vague.R) | 将超先验改为 `p(alpha,beta) ∝ (alpha+beta)^(-5/2)`，比较超参数后验 | [等高线](hw5/hyperparameter_contours_vague.png)、[散点图](hw5/hyperparameter_scatter_vague.png) |
| [Bikes_poisson_joint.R](hw5/Bikes_poisson_joint.R) | 抽取 1000 组超参数，再由条件后验 `Gamma(alpha+y[j], rate=beta+1)` 抽取各组强度，展示后验均值及 95% 区间 | [联合后验分解图](hw5/joint_theta_alpha_beta_decomp.png) |

Poisson 模型中的 `y[j]` 是交通总量，`theta[j]` 是对应观测时段的期望计数；与二项模型中“自行车数量”和“自行车比例”的含义不同。Gamma 分布采用 shape–rate 参数化。

实现备注：当前三份 Poisson 脚本直接用原尺度后验密度作为等距 `(log(alpha), log(beta))` 网格的抽样权重；若目标密度定义在 `(alpha,beta)` 上，抽样权重还需乘 Jacobian `alpha*beta`（对数权重增加 `log(alpha)+log(beta)`）。因此现有抽样摘要与区间需在修正后重新核验；仅用对数坐标展示原尺度密度等高线不需要这一调整。

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
Rscript hw4/Bikes_Expectation.R
```

修正上述数据引用后，也可分别运行 `Rscript hw4/Bikes_Vague.R` 与 `Rscript hw4/Bikes_Uniform.R`。

交互式运行时显示相应图形；通过 `Rscript` 运行时，默认图形设备通常将图写入当前目录的 `Rplots.pdf`，连续运行可能覆盖该文件。脚本不会自动更新仓库中已有的 PNG 图。

hw5 同样使用 R 基础函数。为使交通数据可被找到，在仓库根目录的终端中执行：

```bash
cd hw5
Rscript Bikes_new_sample.R
Rscript Bikes_poisson.R
Rscript Bikes_poisson_vague.R
Rscript Bikes_poisson_joint.R
```

`Bikes_new_sample.R` 将图发送到当前图形设备，不会自动保存已有 PNG；其未显式 `print()` 的摘要表达式可在 RStudio 中逐段执行查看。两份 Poisson 超参数脚本将 PNG 保存到当前工作目录，联合抽样脚本将分解图保存到脚本所在目录。

本仓库用于课程学习与数值实验，代码随学习进度更新。
