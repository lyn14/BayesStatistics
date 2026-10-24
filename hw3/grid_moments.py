"""用离散网格抽样估计 N(0,1) 与 Unif(-1,1) 的前六阶矩，并比较网格设置的影响。"""

import random
from collections.abc import Callable
from math import exp


SAMPLE_COUNT = 1000
RANDOM_SEED = 2026
GRID_RANGES = ((-0.1, 0.1), (-1.0, 1.0), (-3.0, 3.0), (-10.0, 10.0))
GRID_SIZES = (10, 100, 1000)


def normal_density(x: float) -> float:
    """标准正态密度的未归一化核（常数因子不影响离散权重）。"""
    return exp(-(x**2) / 2)


def uniform_density(x: float) -> float:
    """Unif(-1,1) 的密度：支撑集内为 1，支撑集外为 0。"""
    return 1.0 if -1.0 <= x <= 1.0 else 0.0


# (分布名, 密度函数, 一至六阶理论矩)
DISTRIBUTIONS: tuple[tuple[str, Callable[[float], float], tuple[float, ...]], ...] = (
    ("N(0,1)", normal_density, (0.0, 1.0, 0.0, 3.0, 0.0, 15.0)),
    # Unif(-1,1): 奇数阶矩为 0，偶数阶矩 E[X^k] = 1/(k+1)。
    ("Unif(-1,1)", uniform_density, (0.0, 1.0 / 3.0, 0.0, 1.0 / 5.0, 0.0, 1.0 / 7.0)),
)


def make_grid(
    lower: float,
    upper: float,
    point_count: int,
    density: Callable[[float], float],
) -> tuple[list[float], list[float]]:
    """生成等距网格，并用密度在网格点处的值作为概率权重。"""
    if point_count < 2:
        raise ValueError("网格点数必须至少为 2。")

    step = (upper - lower) / (point_count - 1)
    points = [lower + index * step for index in range(point_count)]
    weights = [density(point) for point in points]

    # 梯形求积对区间两端的权重减半。
    weights[0] *= 0.5
    weights[-1] *= 0.5
    if sum(weights) <= 0.0:
        raise ValueError("网格上总权重为 0：该区间没有覆盖分布的支撑集。")
    return points, weights


def estimate_moments(
    lower: float,
    upper: float,
    point_count: int,
    sample_count: int,
    density: Callable[[float], float],
    rng: random.Random,
) -> tuple[float, ...]:
    """从归一化的离散网格分布中有放回抽样，并估计一至六阶原点矩。"""
    if sample_count < 1:
        raise ValueError("样本数必须至少为 1。")

    points, weights = make_grid(lower, upper, point_count, density)
    samples = rng.choices(points, weights=weights, k=sample_count)
    return tuple(
        sum(sample**power for sample in samples) / sample_count
        for power in range(1, 7)
    )


def run_distribution(
    name: str,
    density: Callable[[float], float],
    exact_moments: tuple[float, ...],
) -> None:
    """对单个分布输出全部区间和网格点数组合的矩估计。"""
    exact_text = ", ".join(
        f"E[X^{power}]={exact:g}" for power, exact in enumerate(exact_moments, start=1)
    )
    print(f"\n=== {name}（理论矩: {exact_text}）===")
    print(
        "区间         网格点       E[X]       E[X^2]     E[X^3]     E[X^4]     "
        "E[X^5]     E[X^6]     平均绝对误差"
    )

    for lower, upper in GRID_RANGES:
        for point_count in GRID_SIZES:
            try:
                estimates = estimate_moments(
                    lower,
                    upper,
                    point_count,
                    SAMPLE_COUNT,
                    density,
                    random.Random(RANDOM_SEED),
                )
            except ValueError:
                # 区间内没有任何网格点落在分布支撑集上，无法构造离散近似。
                interval = f"[{lower:g}, {upper:g}]"
                print(
                    f"{interval:<13} {point_count:>4}       "
                    "（支撑集内无网格点，无法估计）"
                )
                continue
            mean_absolute_error = sum(
                abs(estimate - exact)
                for estimate, exact in zip(estimates, exact_moments)
            ) / len(exact_moments)
            interval = f"[{lower:g}, {upper:g}]"
            values = " ".join(f"{estimate:10.4f}" for estimate in estimates)
            print(
                f"{interval:<13} {point_count:>4}       {values} "
                f"{mean_absolute_error:10.4f}"
            )


def main() -> None:
    """运行两种分布的全部组合，并输出综合结论。"""
    print(f"每组随机样本数: {SAMPLE_COUNT}（随机种子: {RANDOM_SEED}）")
    for name, density, exact_moments in DISTRIBUTIONS:
        run_distribution(name, density, exact_moments)

    print("\n综合结论:")
    print("1. 区间必须覆盖分布的有效支撑。N(0,1) 的概率散布在整个实轴: [-3,3] 之外")
    print("   虽然只有约 0.27% 的概率，却贡献了约 1/4 的 E[X^6]，因此 [-3,3] 的")
    print("   E[X^6] 只有约 11 而非 15; Unif(-1,1) 的支撑恰为 [-1,1]，区间取")
    print("   [-1,1] 就没有截断偏差，更宽的区间不增加偏差，但会稀释网格点。")
    print("2. 网格点数决定支撑内的分辨率，关键是落在密度非零区域的点数: [-10,10]")
    print("   只有 10 个点时，Unif(-1,1) 在支撑内一个网格点都没有，程序直接报错;")
    print("   点数增加到 100/1000 后，支撑内的点变密，误差迅速下降。")
    print("3. 区间过窄造成的截断偏差无法靠增加点数弥补（[-0.1,0.1] 对两种分布")
    print("   都系统性低估各阶矩）; 点数过稀的问题则可以靠加多点或收窄区间解决。")
    print("4. 1000 次抽样带来量级约 1/sqrt(1000) 的随机误差，是两种分布共同的")
    print("   精度下限，单次运行的误差排名会受随机波动影响。总体原则是区间匹配")
    print("   有效支撑且点数尽量多: N(0,1) 用 [-10,10]x1000，Unif(-1,1) 用")
    print("   [-1,1]x1000（或更宽区间配 1000 点）。")


if __name__ == "__main__":
    main()
