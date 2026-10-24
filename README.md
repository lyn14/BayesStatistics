# BayesStatistics

贝叶斯统计课程的学习与作业代码，按作业编号整理，包含共轭推断、后验预测和离散网格抽样等示例。

## 内容

| 目录 | 主要内容 |
| --- | --- |
| `hw2/` | Poisson–Gamma 共轭推断、负二项后验预测、数值积分验证，以及课程 PDF |
| `hw3/` | 通过离散网格抽样估计标准正态与均匀分布的前六阶原点矩，比较区间与网格点数的影响；附 PDF 和运行日志 |
| `hw4/` | 层次模型的 R 代码草稿，包含大鼠实验数据；函数尚未完成，暂不能直接运行 |

## 大鼠肿瘤数据

`hw4/Hierarchical_Models_lec.R` 中的数据来自 Gelman 等人的《Bayesian Data Analysis》第三版（BDA3）第 5 章表 5.1 的大鼠肿瘤示例，历史实验数据源自 Tarone（1982）。教材示例包含 70 组历史实验及一组新增实验，用于说明二项分布的层次贝叶斯建模。

- `rat.y[j]`：第 `j` 组实验中出现肿瘤的大鼠数量。
- `rat.n[j]`：第 `j` 组实验的大鼠总数。

两者是观测数据；`y[j] / n[j]` 是该组观察到的肿瘤发生比例。常用模型为 `y[j] ~ Binomial(n[j], theta[j])`，其中 `theta[j]` 才是待推断的该组肿瘤发生概率。

参考：[数据说明与文献出处](https://paulnorthrop.github.io/bang/reference/rat.html)；Tarone, R. E. (1982). *The use of historical information in testing for a trend in proportions*. Biometrics, 38, 215–220. [DOI](https://doi.org/10.2307/2530304)。

## 运行

Python 脚本建议使用 Python 3.9 或更高版本。在仓库根目录执行：

```bash
python -m pip install scipy sympy
python "hw2/Poisson-Gamma posterior.py"
python hw3/Normal_grid_sampling.py
python hw3/grid_moments.py
```

`hw3` 的两个脚本仅依赖 Python 标准库。网格抽样示例使用固定随机种子 `2026`，输出矩估计及其与理论值的误差。

本仓库用于课程学习与数值实验，代码随学习进度更新。
