# ============================================================================
# Bikes_poisson_vague.R
# BDA3 第 5 章 · 层次 Poisson 模型 · 第 (b) 问
#
# 与 Bikes_poisson.R 的唯一区别是超先验：
#   平坦先验 p(alpha, beta) ∝ 1 不合适（p(phi | y) 的积分不收敛），
#   改用 BDA3 §5.3 的经典非信息先验
#       p(alpha, beta) ∝ (alpha + beta)^(-5/2).
#
# 模型：
#   y_j | theta_j ~ Poisson(theta_j),              j = 1, ..., J
#   theta_j | alpha, beta ~ Gamma(alpha, beta)
#   p(alpha, beta) ∝ (alpha + beta)^(-5/2)
#
# 把 theta_1, ..., theta_J 解析积掉（Gamma-Poisson 共轭）：
#
#   p(alpha, beta | y) ∝ p(alpha, beta) *
#        prod_j beta^alpha * Gamma(alpha + y_j) / (Gamma(alpha) * y_j! * (beta+1)^(alpha+y_j))
#
# 这正是第 (b) 问要画的超参数边际后验。它在 (alpha, beta) 上是二维的，
# 本脚本用 (log alpha, log beta) 上的规则网格数值计算它：
#   1) 画等高线（相对最高点的 5%, 15%, ..., 95%）；
#   2) 从网格离散后验中抽 10,000 个 (alpha, beta)，画散点图。
#
# 数据：工作目录下的 Total_traffic.txt（本机实际文件名为
#       Traffic_total.txt，脚本两个名字都试），10 个路口的车辆计数。
#
# 输出：hyperparameter_contours_vague.png, hyperparameter_scatter_vague.png
# ============================================================================


# ---------------------------------------------------------------------------
# 0. 读数据
#    文件里用中文全角逗号分隔，所以不按逗号切分，直接抽出所有数字，
#    这样半角/全角逗号、空格、换行都不影响。
# ---------------------------------------------------------------------------
read_counts <- function(candidates) {
  f <- candidates[file.exists(candidates)]
  if (!length(f)) {
    f <- list.files(".", pattern = "traffic.*\\.txt$", ignore.case = TRUE)
  }
  if (!length(f)) {
    stop("找不到数据文件（试过：", paste(candidates, collapse = ", "), "）")
  }
  txt <- paste(readLines(f[1], warn = FALSE), collapse = " ")
  num <- regmatches(txt, gregexpr("[0-9]+(\\.[0-9]+)?", txt))[[1]]
  if (!length(num)) stop("在 ", f[1], " 里没有找到数字")
  list(y = as.numeric(num), file = f[1])
}

dat <- read_counts(c("Total_traffic.txt", "Traffic_total.txt"))
y <- dat$y
J <- length(y)
if (J < 2) stop("至少需要两组数据")
if (any(y < 0)) stop("Poisson 计数不能为负")

cat("数据文件：", dat$file, "\n", sep = "")
cat("y =", paste(y, collapse = ", "), "\n")
cat("J =", J, "  sum(y) =", sum(y), "  mean(y) =", round(mean(y), 3), "\n")


# ---------------------------------------------------------------------------
# 1. 边际后验的对数（差一个常数）
#    p(alpha, beta | y) ∝ p(alpha, beta) * prod_j [ beta^alpha Gamma(alpha+y_j) /
#                                                   (Gamma(alpha) y_j! (beta+1)^(alpha+y_j)) ]
#    超先验：p(alpha, beta) ∝ (alpha + beta)^(-5/2)
# ---------------------------------------------------------------------------
log.hyperprior <- function(alpha, beta) -2.5 * log(alpha + beta)

log.post.ab <- function(alpha, beta) {
  log.hyperprior(alpha, beta) +
    J * alpha * log(beta) - J * lgamma(alpha) -
    (J * alpha + sum(y)) * log(beta + 1) +
    sum(lgamma(alpha + y)) - sum(lgamma(y + 1))
}


# ---------------------------------------------------------------------------
# 2. 在 (log alpha, log beta) 上建网格
#    后验在 log 尺度上近似为一条窄的对角线带状分布，log 网格比原始尺度
#    均匀得多。这一版先验把后验整体推向更小的 (alpha, beta)，所以网格
#    下界取得比 Bikes_poisson.R 更低，保证几乎不截断尾部质量。
# ---------------------------------------------------------------------------
a.range <- c(0.1, 40)
b.range <- c(0.0005, 0.5)
u.grid <- seq(log(a.range[1]), log(a.range[2]), length.out = 520)   # u = log alpha
v.grid <- seq(log(b.range[1]), log(b.range[2]), length.out = 520)   # v = log beta
nu <- length(u.grid)
nv <- length(v.grid)

# 与 dim(z) = c(nu, nv) 的列优先顺序对齐：u 变化最快，v 最慢
A <- exp(rep(u.grid, times = nv))
B <- exp(rep(v.grid, each = nu))

logz <- log.hyperprior(A, B) +
  J * A * log(B) - J * lgamma(A) - (J * A + sum(y)) * log(B + 1)
for (yj in y) logz <- logz + lgamma(A + yj)        # prod_j Gamma(alpha + y_j)
logz <- logz - sum(lgamma(y + 1))                  # 常数 y_j!，对后验形状无影响

dim(logz) <- c(nu, nv)
z2 <- exp(logz - max(logz))                        # 相对高度，最大值为 1


# ---------------------------------------------------------------------------
# 3. 众数（用 optim 在网格最大值附近精修）
# ---------------------------------------------------------------------------
ij.max <- which(z2 == max(z2), arr.ind = TRUE)[1, ]
fit <- optim(c(u.grid[ij.max[1]], v.grid[ij.max[2]]),
             function(p) -log.post.ab(exp(p[1]), exp(p[2])),
             method = "BFGS", control = list(maxit = 1000))
mode.a <- exp(fit$par[1])
mode.b <- exp(fit$par[2])


# ---------------------------------------------------------------------------
# 4. 从离散后验抽样
#    在 log 尺度上每个网格单元面积相同，所以抽中格点的概率 ∝ z2；
#    再在单元内加均匀抖动，近似连续后验（否则散点图会出现网格花纹）。
# ---------------------------------------------------------------------------
set.seed(2024)
ndraw <- 10000

cell <- sample.int(length(z2), ndraw, replace = TRUE, prob = as.vector(z2))
ij <- arrayInd(cell, dim(z2))
du <- u.grid[2] - u.grid[1]
dv <- v.grid[2] - v.grid[1]
u.draw <- u.grid[ij[, 1]] + runif(ndraw, -du / 2, du / 2)
v.draw <- v.grid[ij[, 2]] + runif(ndraw, -dv / 2, dv / 2)
a.draw <- exp(u.draw)
b.draw <- exp(v.draw)


# ---------------------------------------------------------------------------
# 5. 画图
#    画图范围随后验位置下移（先验把后验推到更小的 alpha, beta）。
# ---------------------------------------------------------------------------
a.lim <- c(0.25, 12)
b.lim <- c(0.0012, 0.12)
a.ticks <- c(0.3, 0.5, 1, 2, 5, 10)
b.ticks <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
b.labels <- c("0.002", "0.005", "0.01", "0.02", "0.05", "0.1")
lev <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)

# 两个图共用的坐标轴
axes.ab <- function() {
  axis(1, at = log(a.ticks), labels = a.ticks, tcl = -0.35, cex.axis = 0.95)
  axis(2, at = log(b.ticks), labels = b.labels, tcl = -0.35, cex.axis = 0.95)
  box()
  mtext(expression(alpha * " (log scale)"), side = 1, line = 3.0, cex = 1.05)
  mtext(expression(beta * " (log scale)"), side = 2, line = 3.4, cex = 1.05, las = 0)
}

draw_contours <- function() {
  par(mar = c(4.7, 5.0, 4.3, 1.6), las = 1, mgp = c(2.8, 0.6, 0),
      xaxs = "i", yaxs = "i")

  image(u.grid, v.grid, z2,
        col = gray.colors(64, start = 0.98, end = 0.68),
        xlim = log(a.lim), ylim = log(b.lim),
        axes = FALSE, xlab = "", ylab = "")
  contour(u.grid, v.grid, z2, levels = lev, add = TRUE,
          drawlabels = FALSE, col = "grey10", lwd = 1.15)
  axes.ab()

  mtext("Marginal posterior density of the hyperparameters", side = 3,
        line = 2.4, font = 2, cex = 1.05)
  mtext("contours at 5%, 15%, ..., 95% of the joint maximum", side = 3,
        line = 1.1, cex = 0.82)

  legend("topleft", bty = "n", cex = 0.82,
         pch = c(NA, 19), lty = c(1, NA), lwd = c(1.15, NA),
         col = c("grey10", 2),
         legend = c("posterior contours",
                    sprintf("mode (%.2f, %.3f)", mode.a, mode.b)))

  points(log(mode.a), log(mode.b), pch = 21, bg = 2, col = "white",
         cex = 1.4, lwd = 1.1)
}

draw_scatter <- function() {
  par(mar = c(4.7, 5.0, 4.3, 1.6), las = 1, mgp = c(2.8, 0.6, 0),
      xaxs = "i", yaxs = "i")

  plot(u.draw, v.draw, type = "n",
       xlim = log(a.lim), ylim = log(b.lim),
       axes = FALSE, xlab = "", ylab = "")
  points(u.draw, v.draw, pch = 16, cex = 0.35,
         col = adjustcolor("steelblue4", alpha.f = 0.10))
  contour(u.grid, v.grid, z2, levels = c(0.05, 0.25, 0.5, 0.75, 0.95),
          add = TRUE, drawlabels = FALSE, col = "grey25", lwd = 1.0)
  axes.ab()

  mtext("Posterior draws of the hyperparameters", side = 3,
        line = 2.4, font = 2, cex = 1.05)
  mtext(sprintf("%s simulations from the grid posterior (with within-cell jitter)",
                format(ndraw, big.mark = ",")),
        side = 3, line = 1.1, cex = 0.82)

  legend("topleft", bty = "n", cex = 0.82,
         pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 1.0),
         col = c(adjustcolor("steelblue4", 0.5), "grey25"),
         legend = c("posterior draws", "contours: 5%-95% of max"))

  points(log(mode.a), log(mode.b), pch = 21, bg = 2, col = "white",
         cex = 1.4, lwd = 1.1)
}

# 保存 PNG；交互式运行（RStudio）时再在屏幕上画一遍
save_fig <- function(file, draw, width = 1600, height = 1400) {
  png(file, width = width, height = height, res = 220, bg = "white")
  draw()
  dev.off()
  if (interactive()) draw()
  cat("已保存：", file, "\n", sep = "")
}

save_fig("hyperparameter_contours_vague.png", draw_contours)
save_fig("hyperparameter_scatter_vague.png", draw_scatter)


# ---------------------------------------------------------------------------
# 6. 数值摘要
# ---------------------------------------------------------------------------
sm <- function(x) c(mean = mean(x), sd = sd(x),
                    "2.5%" = unname(quantile(x, 0.025)),
                    "50%"  = unname(quantile(x, 0.50)),
                    "97.5%" = unname(quantile(x, 0.975)))

tab <- rbind(alpha = sm(a.draw),
             beta = sm(b.draw),
             "alpha/beta" = sm(a.draw / b.draw))

cat("\nHyperprior: p(alpha, beta) ∝ (alpha + beta)^(-5/2)\n")
cat("Posterior mode: alpha =", round(mode.a, 3),
    " beta =", round(mode.b, 4),
    " alpha/beta =", round(mode.a / mode.b, 1), "\n\n")
print(round(tab, 4))
cat("\ncor(log alpha, log beta) =", round(cor(u.draw, v.draw), 3), "\n")
cat("Draws outside the plotting window:",
    sum(a.draw < a.lim[1] | a.draw > a.lim[2] |
          b.draw < b.lim[1] | b.draw > b.lim[2]), "of", ndraw, "\n")
