# ==========================================
# 1. 准备数据（使用课件中的数据）
# ==========================================
bikes.y <- c(16,9,10,13,19,20,18,17,35,55)
bikes.n <- c(74,99,58,70,121,77,104,129,308,119)
J <- length(bikes.y)


# ==========================================
# 2. 定义两种超先验的对数后验函数
# ==========================================
# 老师的 Vague 先验：p(u,v) ∝ αβ(α+β)^{-5/2}
logpost_vague <- function(u, v) {
  a <- (exp(u) * exp(v))/(exp(u) + 1)
  b <- exp(v)/(exp(u) + 1)
  x <- J * (lgamma(a + b) - lgamma(a) - lgamma(b)) + log(a) + log(b) - 2.5 * log(a + b)
  for(i in 1:J) {
    x <- x + lgamma(a + bikes.y[i]) + lgamma(b + bikes.n[i] - bikes.y[i]) - lgamma(a + b + bikes.n[i])
  }
  return(x)
}

# 受限 Uniform 先验：α ~ U(0,40), β ~ U(0,100)
logpost_uniform <- function(u, v) {
  a <- (exp(u) * exp(v))/(exp(u) + 1)
  b <- exp(v)/(exp(u) + 1)
  
  # 用 dunif 计算先验的对数。超出边界的会自动变成 -Inf
  prior_log <- dunif(a, min = 0, max = 40, log = TRUE) + 
    dunif(b, min = 0, max = 100, log = TRUE)
  
  x <- J * (lgamma(a + b) - lgamma(a) - lgamma(b)) + prior_log
  for(i in 1:J) {
    x <- x + lgamma(a + bikes.y[i]) + lgamma(b + bikes.n[i] - bikes.y[i]) - lgamma(a + b + bikes.n[i])
  }
  return(x)
}


# ==========================================
# 3. 网格抽样与计算 E(θ_j | y)
# ==========================================
# 设定网格（和老师代码一致）
u_grid <- seq(-2.3, -1.3, length.out = 100) 
v_grid <- seq(1.0, 5.0, length.out = 200)   

# 定义一个辅助函数：从二维网格中按概率抽取 K 组 (α, β)
sample_grid <- function(u_grid, v_grid, logpost_func, K = 1000) {
  z <- outer(u_grid, v_grid, logpost_func) # 计算高度
  z <- exp(z - max(z))                     # 转化为相对概率
  probs <- as.vector(z) / sum(z)           # 归一化为总和=1的权重
  
  # 抽样
  idx <- sample(length(probs), size = K, replace = TRUE, prob = probs)
  
  # 还原行号和列号，并映射回 u 和 v
  n_rows <- length(u_grid)
  row_idx <- (idx - 1) %% n_rows + 1
  col_idx <- (idx - 1) %/% n_rows + 1
  
  u_samples <- u_grid[row_idx]
  v_samples <- v_grid[col_idx]
  
  # 换算回 α 和 β
  a_samples <- (exp(u_samples) * exp(v_samples)) / (exp(u_samples) + 1)
  b_samples <- exp(v_samples) / (exp(u_samples) + 1)
  
  return(data.frame(a = a_samples, b = b_samples))
}

# 计算 E(θ_j | y) 的函数
calc_E_theta <- function(samples, bikes.y, bikes.n) {
  J <- length(bikes.y)
  E_theta <- numeric(J)
  for (j in 1:J) {
    # 根据公式：E(θ_j | α, β, y) = (α + y_j) / (α + β + n_j)
    # 我们对抽取的 1000 个条件期望取平均
    cond_means <- (samples$a + bikes.y[j]) / (samples$a + samples$b + bikes.n[j])
    E_theta[j] <- mean(cond_means)
  }
  return(E_theta)
}

# 设置随机种子，确保结果可重复
set.seed(123)

# 抽取样本并计算
samples_vague <- sample_grid(u_grid, v_grid, logpost_vague, K = 1000)
E_theta_vague <- calc_E_theta(samples_vague, bikes.y, bikes.n)

samples_uniform <- sample_grid(u_grid, v_grid, logpost_uniform, K = 1000)
E_theta_uniform <- calc_E_theta(samples_uniform, bikes.y, bikes.n)

# 计算 MLE (y / n)
mle <- bikes.y / bikes.n


# ==========================================
# 4. 绘图 (对应课件第 41 页的图)
# ==========================================
# 画空图打底
plot(mle, E_theta_vague, col = "red", pch = 16, cex = 0.8,
     xlab = "MLE: y_i / n_i", ylab = "E(θ_i | y)", 
     main = "Posterior Mean vs MLE",
     ylim = c(0, max(E_theta_vague, E_theta_uniform) * 1.1),
     xlim = c(0, max(mle) * 1.1))

# 加上 Uniform 的黑点
points(mle, E_theta_uniform, col = "black", pch = 16, cex = 0.8)

# 添加 y = x 的对角线 (代表没有收缩的情况，即完全相信 MLE)
abline(0, 1, col = "gray", lty = 2)

# 添加图例
legend("topleft", legend = c("Vague", "Uniform"), 
       col = c("red", "black"), pch = 16, bty = "n")


