from statistics import NormalDist
from scipy.stats import gamma

# --1 导入数据--
# 两列表格：第一列年份(1998-2017)，第二列事故数量(待填写)
Data = {
    1998: 20,  
    1999: 21,  # 事故数量
    2000: 18,  # 事故数量
    2001: 14,  # 事故数量
    2002: 14,  # 事故数量
    2003: 7,  # 事故数量
    2004: 9,  # 事故数量
    2005: 17,  # 事故数量
    2006: 23,
    2007: 11
}

# --2 构造共轭先验 --
# 认为出事故的概率比较小，采取泊松分布，共轭先验为Gamma分布.根据要求，置参数为2，2
alpha_prior, beta_prior = 2, 2   # Gamma(shape=2, rate=2)，密度 ∝ theta^(2-1) * exp(-2*theta)

# --3 贝叶斯推断 --
# 似然函数(泊松)：L(theta) ∝ prod_i exp(-theta) * theta^(n_i) / n_i!
#                              = exp(-N*theta) * theta^(sum n_i) / prod n_i!
# 先验(Gamma)：p(theta) ∝ theta^(2-1) * exp(-2*theta)
# 后验 ∝ L(theta) * p(theta)
#      ∝ exp{-theta(N+2)} * theta^(sum n_i + 1)   <- 后验核心项
# 识别为 Gamma 分布：shape = sum(n_i) + 2, rate = N + 2

N = len(Data)               # 数据总数（年数）
S = sum(Data.values())      # sigma_ni：各年事故数之和

alpha_post = S + alpha_prior    # 后验 shape = sum(n_i) + 2
beta_post  = N + beta_prior     # 后验 rate  = N + 2

print(f"N = {N}, sum(n_i) = {S}")
print(f"后验核心项: exp(-theta*{beta_post}) * theta^({S}+1)")
print(f"后验分布: Gamma(shape={alpha_post}, rate={beta_post})")
print(f"后验均值 E[theta|data] = {alpha_post}/{beta_post} = {alpha_post/beta_post:.4f}")

# --3.1 用 sympy 符号验证核心项 --
import sympy as sp

theta = sp.symbols('theta', positive=True)

# 似然 ∝ prod_i exp(-theta) * theta^(n_i)（省略常数 1/prod n_i!）
likelihood = sp.prod([sp.exp(-theta) * theta**ni for ni in Data.values()])
# 先验 ∝ theta^(2-1) * exp(-2*theta)
prior = theta**(alpha_prior - 1) * sp.exp(-beta_prior * theta)

posterior_kernel = sp.expand(likelihood * prior)
core_term = sp.exp(-theta * (N + beta_prior)) * theta**(S + alpha_prior - 1)

# 两者只差一个常数因子（1/prod n_i!），即成比例
ratio = sp.simplify(posterior_kernel / core_term)
print(f"后验核/核心项 = {ratio}  (与 theta 无关的常数，验证核心项成立)")

# --4 计算后验预测分布 n_tilde | n --
# 定义：p(n_tilde|n) = ∫ p(n_tilde|theta) · p(theta|n) dtheta
# 代入泊松似然与 Gamma 后验：
#   p(n_tilde|n) = ∫ [e^{-theta}·theta^{n_t}/n_t!] · [beta^alpha/Gamma(alpha)·theta^{alpha-1}·e^{-beta·theta}] dtheta
#                = Gamma(n_t+alpha)/(n_t!·Gamma(alpha)) · (beta/(beta+1))^alpha · (1/(beta+1))^{n_t}
# 识别为负二项分布：n_tilde | n ~ NegBinom(r=alpha_post, p=beta_post/(beta_post+1))
# scipy 参数化 nbinom.pmf(k, n, p)：均值 = n(1-p)/p = alpha_post/beta_post（恰为后验均值）

from scipy.stats import nbinom, poisson
from scipy import integrate
import math

r_pred = alpha_pred = alpha_post         # 负二项 size 参数 r = 156
p_pred = beta_post / (beta_post + 1)     # 成功概率 p = 12/13

print(f"\n后验预测分布: n_tilde|n ~ NegBinom(r={r_pred}, p={p_pred:.6f})")
print(f"预测均值 E[n_tilde|n] = r(1-p)/p = {r_pred*(1-p_pred)/p_pred:.4f}  (= 后验均值 {alpha_post/beta_post:.4f})")
print(f"预测方差 Var[n_tilde|n] = r(1-p)/p^2 = {r_pred*(1-p_pred)/p_pred**2:.4f}")

# 95% 预测区间
lo = nbinom.ppf(0.025, r_pred, p_pred)
hi = nbinom.ppf(0.975, r_pred, p_pred)
print(f"95% 预测区间: [{lo:.0f}, {hi:.0f}]")

# --4.1 数值积分验证闭式解 --
# p(n_t|n) = ∫ Poisson(n_t|theta) · Gamma(theta|alpha_post, beta_post) dtheta
def predictive_integral(n_t):
    integrand = lambda th: poisson.pmf(n_t, th) * gamma.pdf(th, a=alpha_post, scale=1/beta_post)
    val, _ = integrate.quad(integrand, 0, math.inf)
    return val

print("\n数值积分 vs 负二项闭式解：")
for n_t in [5, 10, 13, 20, 25]:
    p_closed = nbinom.pmf(n_t, r_pred, p_pred)
    p_num = predictive_integral(n_t)
    print(f"  n_tilde={n_t:3d}: 闭式={p_closed:.6f}, 数值积分={p_num:.6f}, 差={abs(p_closed-p_num):.2e}")

# --4.2 预测分布 PMF 一览（均值附近） --
print("\n后验预测分布 PMF（n_tilde = 7..20）：")
for n_t in range(7, 21):
    pmf_val = nbinom.pmf(n_t, r_pred, p_pred)
    print(f"  {n_t:3d}: {pmf_val:.4f} {'█' * int(round(pmf_val * 100))}")

# --5 计算题目要求的：来年事故数 > 20 的概率 --
# P(n_tilde > 20 | n) = ∫ P(n_tilde > 20 | theta) · p(theta|n) dtheta
# 方法1(闭式)：后验预测为负二项分布，直接用生存函数
#   P(n_tilde > 20 | n) = nbinom.sf(20, r_pred, p_pred)
# 方法2(数值积分)：∫ [1 - Poisson.CDF(20|theta)] · Gamma(theta|156,12) dtheta
#   其中 1 - Poisson.CDF(20|theta) = poisson.sf(20, theta) 是给定 theta 时泊松的尾部概率

threshold = 20

# 方法1：负二项生存函数（闭式解）
p_exceed_closed = nbinom.sf(threshold, r_pred, p_pred)

# 方法2：数值积分验证
integrand = lambda th: poisson.sf(threshold, th) * gamma.pdf(th, a=alpha_post, scale=1/beta_post)
p_exceed_num, _ = integrate.quad(integrand, 0, math.inf)

print(f"\n来年事故数 > {threshold} 的概率 P(n_tilde > {threshold} | n):")
print(f"  负二项闭式解 nbinom.sf: {p_exceed_closed:.6f}")
print(f"  数值积分验证:          {p_exceed_num:.6f}")
print(f"  两者之差:              {abs(p_exceed_closed - p_exceed_num):.2e}")
print(f"  即约 {p_exceed_closed*100:.2f}%")

# 对比：若只用后验均值 theta=13 的泊松分布（忽略参数不确定性）
p_exceed_poisson_at_mean = poisson.sf(threshold, alpha_post/beta_post)
print(f"\n[对比] 若忽略 theta 不确定性(固定 theta=13 的泊松): {p_exceed_poisson_at_mean:.6f}")
print("  负二项尾部更厚，因为额外叠加了 theta 本身的不确定性")
 



