# ============================================================================
# Bikes_poisson_joint.R
# BDA3 第 5 章 · 层次 Poisson 模型 · 第 (b) 问的后续：抽样 (theta, alpha, beta)
#
# 模型与超先验（与 Bikes_poisson_vague.R 一致）：
#   y_j | theta_j ~ Poisson(theta_j)
#   theta_j | alpha, beta ~ Gamma(alpha, beta)
#   p(alpha, beta) ∝ (alpha + beta)^(-5/2)
#
# 边际后验（theta 已积掉）：
#   p(alpha, beta | y) ∝ p(alpha, beta) *
#        prod_j beta^alpha Gamma(alpha + y_j) / (Gamma(alpha) y_j! (beta+1)^(alpha+y_j))
#
# 由 Gamma-Poisson 共轭，给定 (alpha, beta) 后每个 theta_j 的条件后验为
#   theta_j | alpha, beta, y_j ~ Gamma(alpha + y_j, beta + 1)
# 且 theta_1, ..., theta_J 在给定 (alpha, beta) 后相互独立。
#
# 抽样步骤（1000 轮）：
#   1. 从网格近似的 p(alpha, beta | y) 中抽 1000 个 (alpha, beta)；
#   2. 每一轮再抽 theta_j ~ Gamma(alpha + y_j, beta + 1)。
# 最后用一张二联图展示联合后验的分解：
#   (a) 超参数 p(alpha, beta | y) 的等高线；
#   (b) 10 个 theta_j 的 95% 后验区间与观测 y_j。
#
# 输出：joint_theta_alpha_beta_decomp.png
# ============================================================================


# ---------------------------------------------------------------------------
# 0. 读数据
#    RStudio 的 working directory 不一定是脚本所在目录，所以除了 getwd()，
#    还会在脚本所在目录（及其上一级）里找。
#    如果数据放在别处，直接把下面的 data.file 改成绝对路径即可。
# ---------------------------------------------------------------------------
data.file <- NULL            # put working directory here

script_dir <- function() {
  # 1) Rscript file.R
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) {
    return(dirname(normalizePath(sub("^--file=", "", m[1]), mustWork = FALSE)))
  }
  # 2) RStudio 里 Source
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  # 3) source("file.R")
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile, mustWork = FALSE)))
  }
  getwd()
}

find_data <- function(candidates) {
  here <- tryCatch(script_dir(), error = function(e) getwd())
  dirs <- unique(c(getwd(), here, dirname(here)))
  dirs <- dirs[dir.exists(dirs)]

  # 先按候选文件名找
  for (d in dirs) {
    hit <- file.path(d, candidates)
    hit <- hit[file.exists(hit)]
    if (length(hit)) return(normalizePath(hit[1], mustWork = FALSE))
  }
  # 再按 *traffic*.txt 模糊匹配（只搜当前目录和脚本目录，避免扫全盘）
  for (d in dirs) {
    hit <- list.files(d, pattern = "traffic.*\\.txt$", ignore.case = TRUE,
                      full.names = TRUE)
    if (length(hit)) return(normalizePath(hit[1], mustWork = FALSE))
  }
  # 最后在脚本目录下递归找（例如 data/ 子目录）
  hit <- list.files(here, pattern = "traffic.*\\.txt$", ignore.case = TRUE,
                    full.names = TRUE, recursive = TRUE)
  if (length(hit)) return(normalizePath(hit[1], mustWork = FALSE))
  NULL
}

read_counts <- function(candidates, file = NULL) {
  if (!is.null(file) && nzchar(file)) {
    if (!file.exists(file)) stop("data.file 指定的文件不存在：", file)
    f <- file
  } else {
    f <- find_data(candidates)
    if (is.null(f)) {
      stop("找不到数据文件（试过：", paste(candidates, collapse = ", "), "）。\n",
           "  当前 working directory：", getwd(), "\n",
           "  脚本所在目录：", tryCatch(script_dir(), error = function(e) "未知"), "\n",
           "  请把数据文件放到上面任一目录，或在脚本顶部设置 data.file <- \"绝对路径\"。")
    }
  }
  txt <- paste(readLines(f[1], warn = FALSE), collapse = " ")
  num <- regmatches(txt, gregexpr("[0-9]+(\\.[0-9]+)?", txt))[[1]]
  if (!length(num)) stop("在 ", f[1], " 里没有找到数字")
  list(y = as.numeric(num), file = f[1])
}

dat <- read_counts(c("Total_traffic.txt", "Traffic_total.txt"), data.file)
y <- dat$y
J <- length(y)
if (J < 2) stop("至少需要两组数据")

cat("数据文件：", dat$file, "\n", sep = "")
cat("y =", paste(y, collapse = ", "), "\n")
cat("J =", J, "  sum(y) =", sum(y), "  mean(y) =", round(mean(y), 3), "\n")


# ---------------------------------------------------------------------------
# 1. 超参数边际后验（与 Bikes_poisson_vague.R 相同）
# ---------------------------------------------------------------------------
log.hyperprior <- function(alpha, beta) -2.5 * log(alpha + beta)

log.post.ab <- function(alpha, beta) {
  log.hyperprior(alpha, beta) +
    J * alpha * log(beta) - J * lgamma(alpha) -
    (J * alpha + sum(y)) * log(beta + 1) +
    sum(lgamma(alpha + y)) - sum(lgamma(y + 1))
}

a.range <- c(0.1, 40)
b.range <- c(0.0005, 0.5)
u.grid <- seq(log(a.range[1]), log(a.range[2]), length.out = 520)
v.grid <- seq(log(b.range[1]), log(b.range[2]), length.out = 520)
nu <- length(u.grid)
nv <- length(v.grid)

A <- exp(rep(u.grid, times = nv))
B <- exp(rep(v.grid, each = nu))
logz <- log.hyperprior(A, B) +
  J * A * log(B) - J * lgamma(A) - (J * A + sum(y)) * log(B + 1)
for (yj in y) logz <- logz + lgamma(A + yj)
logz <- logz - sum(lgamma(y + 1))
dim(logz) <- c(nu, nv)
z2 <- exp(logz - max(logz))


# ---------------------------------------------------------------------------
# 2. 抽 1000 个 (alpha, beta)
#    网格在 (log alpha, log beta) 上等面积，故按 z2 抽格点并在格内抖动。
# ---------------------------------------------------------------------------
set.seed(2024)
ndraw <- 1000

cell <- sample.int(length(z2), ndraw, replace = TRUE, prob = as.vector(z2))
ij <- arrayInd(cell, dim(z2))
du <- u.grid[2] - u.grid[1]
dv <- v.grid[2] - v.grid[1]
u.draw <- u.grid[ij[, 1]] + runif(ndraw, -du / 2, du / 2)
v.draw <- v.grid[ij[, 2]] + runif(ndraw, -dv / 2, dv / 2)
a.draw <- exp(u.draw)
b.draw <- exp(v.draw)


# ---------------------------------------------------------------------------
# 3. 给定 (alpha, beta) 抽 theta_j
#    theta_j | alpha, beta, y_j ~ Gamma(shape = alpha + y_j, rate = beta + 1)
#    R 的 rgamma 第二个参数是 rate，所以 rate = beta + 1。
# ---------------------------------------------------------------------------
theta.mat <- sapply(y, function(yj) rgamma(ndraw, shape = a.draw + yj, rate = b.draw + 1))
# theta.mat 是 ndraw x J 矩阵，第 j 列是 theta_j 的 1000 个抽样


# ---------------------------------------------------------------------------
# 4. 分解图（替代 3D 散点）
#    p(theta, alpha, beta | y) = p(alpha, beta | y) * prod_j p(theta_j | alpha, beta, y)
#    (a) 超参数联合后验的等高线；(b) 10 个 theta_j 的 95% 区间与观测 y_j。
# ---------------------------------------------------------------------------
# 后验众数（在网格最大值附近用 optim 精修）
ij.max <- which(z2 == max(z2), arr.ind = TRUE)[1, ]
fit <- optim(c(u.grid[ij.max[1]], v.grid[ij.max[2]]),
             function(p) -log.post.ab(exp(p[1]), exp(p[2])),
             method = "BFGS", control = list(maxit = 1000))
mode.a <- exp(fit$par[1])
mode.b <- exp(fit$par[2])

a.lim <- c(0.25, 12)
b.lim <- c(0.0012, 0.12)
a.ticks <- c(0.3, 0.5, 1, 2, 5, 10)
b.ticks <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
b.labels <- c("0.002", "0.005", "0.01", "0.02", "0.05", "0.1")
lev <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)

theta.mean <- colMeans(theta.mat)
theta.lo <- apply(theta.mat, 2, quantile, probs = 0.025)
theta.hi <- apply(theta.mat, 2, quantile, probs = 0.975)
pop.mean <- mean(a.draw / b.draw)

draw_decomp <- function() {
  layout(matrix(1:2, 1), widths = c(1, 1.05))

  # (a) 超参数后验
  par(mar = c(4.4, 4.6, 3.6, 1.0), las = 1, mgp = c(2.6, 0.6, 0),
      xaxs = "i", yaxs = "i")
  image(u.grid, v.grid, z2,
        col = gray.colors(64, start = 0.98, end = 0.68),
        xlim = log(a.lim), ylim = log(b.lim),
        axes = FALSE, xlab = "", ylab = "")
  contour(u.grid, v.grid, z2, levels = lev, add = TRUE,
          drawlabels = FALSE, col = "grey10", lwd = 1.1)
  axis(1, at = log(a.ticks), labels = a.ticks, tcl = -0.3, cex.axis = 0.85)
  axis(2, at = log(b.ticks), labels = b.labels, tcl = -0.3, cex.axis = 0.85)
  box()
  mtext(expression(alpha * " (log scale)"), side = 1, line = 2.9, cex = 0.95)
  mtext(expression(beta * " (log scale)"), side = 2, line = 3.2, cex = 0.95, las = 0)
  mtext("(a) hyperparameter posterior p(alpha, beta | y)", side = 3,
        line = 1.2, font = 2, cex = 0.95)
  legend("topleft", bty = "n", cex = 0.78,
         pch = c(NA, 19), lty = c(1, NA), col = c("grey10", 2),
         legend = c("contours: 5%-95% of max",
                    sprintf("mode (%.2f, %.3f)", mode.a, mode.b)))
  points(log(mode.a), log(mode.b), pch = 21, bg = 2, col = "white", cex = 1.3)

  # (b) theta_j 后验区间
  par(mar = c(4.4, 4.4, 3.6, 1.0), las = 1, mgp = c(2.6, 0.6, 0),
      xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(0.5, J + 0.5), ylim = c(40, 360), log = "y",
       xaxt = "n", yaxt = "n", xlab = "location j", ylab = expression(theta[j]),
       cex.axis = 0.85, cex.lab = 0.95)
  axis(1, at = 1:J, labels = 1:J, tcl = -0.3, cex.axis = 0.85)
  axis(2, at = c(50, 75, 100, 150, 200, 300),
       labels = c(50, 75, 100, 150, 200, 300), las = 1, cex.axis = 0.85)
  abline(h = pop.mean, lty = 2, col = 2, lwd = 1.2)
  segments(1:J, theta.lo, 1:J, theta.hi, lwd = 2.2, col = "steelblue4")
  points(1:J, theta.mean, pch = 19, col = "steelblue4", cex = 1.1)
  points(1:J, y, pch = 4, col = "grey20", cex = 1.05, lwd = 1.2)
  mtext("(b) theta_j 95% intervals vs observed y_j", side = 3,
        line = 1.2, font = 2, cex = 0.95)
  legend("topleft", bty = "n", cex = 0.75,
         legend = c("95% credible interval", "posterior mean",
                    "observed y_j", sprintf("E[alpha/beta] = %.0f", pop.mean)),
         lty = c(1, NA, NA, 2), lwd = c(2.2, NA, NA, 1.2),
         pch = c(NA, 19, 4, NA),
         col = c("steelblue4", "steelblue4", "grey20", 2))
}

save_fig <- function(file, draw, width = 2000, height = 950) {
  png(file, width = width, height = height, res = 220, bg = "white")
  draw()
  dev.off()
  if (interactive()) draw()
  cat("已保存：", file, "\n", sep = "")
}

save_fig(file.path(script_dir(), "joint_theta_alpha_beta_decomp.png"), draw_decomp)


# ---------------------------------------------------------------------------
# 5. 数值摘要
# ---------------------------------------------------------------------------
cat("\nHyperparameters: mean alpha =", round(mean(a.draw), 3),
    " mean beta =", round(mean(b.draw), 4),
    " mean alpha/beta =", round(pop.mean, 1), "\n")

tab <- data.frame(
  j = 1:J,
  y = y,
  mean = theta.mean,
  sd = apply(theta.mat, 2, sd),
  "2.5%" = theta.lo,
  "50%" = apply(theta.mat, 2, quantile, probs = 0.5),
  "97.5%" = theta.hi,
  check.names = FALSE
)
cat("\nPosterior of theta_j (1000 draws):\n")
print(round(tab, 2), row.names = FALSE)

cat("\ncor(theta_j, alpha):", paste(round(apply(theta.mat, 2, cor, a.draw), 3),
                                   collapse = "  "), "\n")
cat("cor(theta_j, beta) :", paste(round(apply(theta.mat, 2, cor, b.draw), 3),
                                   collapse = "  "), "\n")
