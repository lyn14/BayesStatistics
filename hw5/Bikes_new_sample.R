# y: The number of bicycles in group;
# n: The number of total vehicles in group
bikes.y <- c(16,9,10,13,19,20,18,17,35,55)
bikes.n <- c(74,99,58,70,121,77,104,129,308,119)
logpost4 <- function(u, v) {a <- exp(u + v)/(exp(u) + 1);
b <- exp(v)/(exp(u) + 1);
J <- length(bikes.y); # Number of data
x <- J * (lgamma(a + b)- lgamma(a)- lgamma(b)) + log(a) + log(b)-2.5 * log(a + b); # Constants, Jaccobi and hyperprior
for (i in (1:J)){
  x <- x + lgamma(a + bikes.y[i]) + lgamma(b + bikes.n[i] - bikes.y[i]) - 
    lgamma(a + b + bikes.n[i]) # Adding likelihood
}
x # Return log p(theta, phi | y)
}
# Genebikese grid sampling
x4 <- c(-230:-50)/100 # set the grid values of u
y4 <- c(100:500)/100 # set the grid values of v
# evaluate grid
z <- outer(x4,y4,logpost4) # evaluate the key part of the log posterior p(u,v|y) for each grid (u,v)
# exponentiate log posterior (set maximum to be one)
z2 <- exp(z- max(z)) # transfer log posterior into posterior (the joint posterior)

contour(x4,y4,z2,levels=c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95))

# ---------------------------------------------------------------------------
# Posterior predictive distribution for a NEW group:  p(y.tilde | y)
# A new location from the same population, n.tilde = 100 vehicles:
#     phi | a, b ~ Beta(a, b)                 (phi = theta of the new group)
#     y.tilde | phi ~ Binomial(n.tilde, phi)
# ---------------------------------------------------------------------------
set.seed(2024)
n.tilde <- 100
ndraw <- 10000

# 1. normalise the grid, then draw (u, v) -- and hence (a, b) -- from it -----
p.joint <- z2/sum(z2)                              # joint posterior on the grid
idx <- sample(length(p.joint), ndraw, replace = TRUE, prob = p.joint)
ij <- arrayInd(idx, dim(p.joint))
u.s <- x4[ij[, 1]]                                 # draws of u = logit(a/(a+b))
v.s <- y4[ij[, 2]]                                 # draws of v = log(a+b)
a.s <- exp(u.s + v.s)/(exp(u.s) + 1)               # transform back to (a, b)
b.s <- exp(v.s)/(exp(u.s) + 1)

# 2. two-stage sampling from the posterior predictive ------------------------
phi <- rbeta(ndraw, a.s, b.s)                      # theta of the new group
ynew <- rbinom(ndraw, n.tilde, phi)                # bicycles out of 100 vehicles

# 3. exact form of the same distribution ------------------------------------
# Integrating phi out of the beta-binomial gives
#   p(y.tilde = k | a, b) = C(n,k) * B(a+k, b+n-k)/B(a,b),  then average over draws of (a,b)
k <- 0:n.tilde
logbb <- function(k, a, b) {
  lgamma(n.tilde + 1) - lgamma(k + 1) - lgamma(n.tilde - k + 1) +
    lgamma(a + b) - lgamma(a) - lgamma(b) +
    lgamma(a + k) + lgamma(b + n.tilde - k) - lgamma(a + b + n.tilde)
}
p.pred <- colMeans(exp(sapply(k, logbb, a = a.s, b = b.s)))
p.pred <- p.pred/sum(p.pred)                       # length 101, k = 0,...,100

# 4. picture ---------------------------------------------------------------
par(mar = c(6.5, 4.2, 4.2, 2.2))                    # extra room for the caption below
emp <- as.numeric(table(factor(ynew, levels = 0:n.tilde)))/ndraw   # empirical pmf
y.max <- max(emp, p.pred)                                          # tallest bar

hist(ynew, breaks = seq(-0.5, n.tilde + 0.5, 1), probability = TRUE,
     col = "grey85", border = "grey60", xlim = c(-0.5, 60.5), ylim = c(0, y.max * 1.45),
     main = "Posterior predictive, new group (n = 100 vehicles)",
     xlab = "y.tilde = number of bicycles", ylab = "probability")
lines(k, p.pred, type = "h", col = 2, lwd = 2)
legend("topright", c("predictive draws (histogram)", "beta-binomial pmf"),
       fill = c("grey85", NA), border = c("grey60", NA), lwd = c(NA, 2),
       col = c(NA, 2), bty = "n", cex = 0.85)

# 5. summaries -------------------------------------------------------------
round(c(mean = mean(ynew), sd = sd(ynew),
        q025 = as.numeric(quantile(ynew, 0.025)),
        q50 = as.numeric(quantile(ynew, 0.5)),
        q975 = as.numeric(quantile(ynew, 0.975))), 2)
round(p.pred[c(1, 11, 21, 31, 41)], 4)             # P(y.tilde = 0, 10, 20, 30, 40)
round(sum(p.pred[1:11]), 4)                        # P(y.tilde <= 10)

# 6. 95% credible interval at both sides of the peak ------------------------
oo <- order(p.pred, decreasing = TRUE)
sel <- oo[seq_len(which(cumsum(p.pred[oo]) >= 0.95)[1])]
k.lo <- min(k[sel])
k.hi <- max(k[sel])
mass.ci <- sum(p.pred[(k.lo + 1):(k.hi + 1)])      # mass actually inside
c(lower = k.lo, upper = k.hi, mass = round(mass.ci, 4))

# frame it on the figure: tinted band + dashed box round the whole region
rect(k.lo - 0.5, 0, k.hi + 0.5, y.max * 1.10,
     col = rgb(1, 0, 0, 0.08), border = 2, lty = 2, lwd = 1.3)
mtext(side = 1, line = 5.1, cex = 0.8, col = 2,
      text = sprintf("dashed box: 95%% peak-centred (HPD) interval [%d, %d], %.1f%% of the mass",
                     k.lo, k.hi, 100 * mass.ci))



