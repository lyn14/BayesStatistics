"""用离散网格抽样估计标准正态分布的前六阶矩。"""

import random
from math import exp


SAMPLE_COUNT = 1000
GRID_RANGES = ((-0.1, 0.1), (-1.0, 1.0), (-3.0, 3.0), (-10.0, 10.0))
GRID_SIZES = (10, 100, 1000)
EXACT_MOMENTS = (0.0, 1.0, 0.0, 3.0, 0.0, 15.0)


def make_grid(lower: float, upper: float, point_count: int) -> tuple[list[float], list[float]]:
    """生成等距网格，并用标准正态密度近似各网格点的概率权重。"""
    if point_count < 2:
        raise ValueError("网格点数必须至少为 2。")

    step = (upper - lower) / (point_count - 1)
    points = [lower + index * step for index in range(point_count)]
    weights = [exp(-(point**2) / 2) for point in points]

    # 梯形求积对区间两端的权重减半。
    weights[0] *= 0.5
    weights[-1] *= 0.5
    return points, weights


def estimate_moments(
    lower: float,
    upper: float,
    point_count: int,
    sample_count: int,
    rng: random.Random,
) -> tuple[float, ...]:
    """从归一化的离散网格分布中有放回抽样，并估计一至六阶原点矩。"""
    if sample_count < 1:
        raise ValueError("样本数必须至少为 1。")

    points, weights = make_grid(lower, upper, point_count)
    samples = rng.choices(points, weights=weights, k=sample_count)
    return tuple(
        sum(sample**power for sample in samples) / sample_count
        for power in range(1, 7)
    )


def main() -> None:
    """运行全部区间和网格点数组合，并输出矩估计及平均绝对误差。"""
    print(f"每组随机样本数: {SAMPLE_COUNT}（随机种子: 2026）")
    print("区间       网格点       E[X]       E[X^2]     E[X^3]     E[X^4]     E[X^5]     E[X^6]     平均绝对误差")
    for lower, upper in GRID_RANGES:
        for point_count in GRID_SIZES:
            estimates = estimate_moments(
                lower, upper, point_count, SAMPLE_COUNT, random.Random(2026)
            )
            mean_absolute_error = sum(
                abs(estimate - exact)
                for estimate, exact in zip(estimates, EXACT_MOMENTS)
            ) / len(EXACT_MOMENTS)
            interval = f"[{lower:g}, {upper:g}]"
            values = " ".join(f"{estimate:10.4f}" for estimate in estimates)
            print(
                f"{interval:<11} {point_count:>4}       {values} "
                f"{mean_absolute_error:10.4f}"
            )

    print("\n标准正态分布的理论矩: E[X]=0, E[X^2]=1, E[X^3]=0, "
          "E[X^4]=3, E[X^5]=0, E[X^6]=15。")
    print("区间过窄会截掉尾部概率，尤其使高阶矩偏小；区间过宽而网格点数少会使网格过稀。")
    print("增加网格点数通常改善离散分布对区间内密度的近似，但不会消除 1000 次抽样的随机误差。")


if __name__ == "__main__":
    main()
