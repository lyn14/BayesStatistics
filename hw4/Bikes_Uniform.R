# y: The number of tumor-emerging bikess in group;
# n: The number of total bikess in group
bikes.y <- c(16,9,10,13,19,20,18,17,35,55)
bikes.n <- c(74,99,58,70,121,77,104,129,308,119)
logpost4 <- function(u, v) {a <- exp(u+v)/(exp(u)+1);
b <- exp(v)/(exp(u) + 1);
J <- length(bikes.y); # Number of data
logp_a <- dunif(a , min = 0 , max = 40 , log = TRUE);
logp_b <- dunif(b , min = 0 , max = 100 , log = TRUE);
log_prior <- logp_a + logp_b;
x <- J * (lgamma(a + b) - lgamma(a) - lgamma(b)) + log(a) + log(b) + log_prior;
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

m=1000
z2 <- z2/sum(z2) #get the value of joint posterior, where sum(z2) serves as normalized constant
p.x <- apply(z2, 1, sum) #get marginal posterior p(u|y) from joint posterior p(u,v|y) above
newx <- sample(x4, m, replace=T, prob=p.x) 
xid <- (newx + 2.31)*100 # round to integer
# sample v
newy <- rep(0, m)
for(i in(1:m)) newy[i] <- sample(y4, 1, prob=z2[xid[i],]/sum(z2[xid[i],]))
newa <- (exp(newx)*exp(newy))/(exp(newx)+1)
newb <- exp(newy)/(exp(newx)+1)

# 1. 设定你要预测的组（比如第 5 组）
j <- 5 
yj <- bikes.y[j]  
nj <- bikes.n[j]  

# 2. 根据抽出来的 alpha 和 beta 样本，进行第二步抽样，生成 1000 个 theta 样本
# 根据条件后验：theta | alpha, beta, y ~ Beta(alpha + yj, beta + nj - yj)
theta_samples <- rbeta(m, shape1 = newa + yj, shape2 = newb + nj - yj)

# 3. 计算这 1000 个 theta 样本的 2.5% 和 97.5% 分位数，这就是 95% 可信区间
ci_95 <- quantile(theta_samples, probs = c(0.025, 0.975))

# 4. 打印结果
print(ci_95)

# 顺便画个直方图看看后验分布的形状
hist(theta_samples, breaks = 30, main = paste("Posterior of theta_", j), 
     xlab = "theta", col = "lightblue", border = "white")
abline(v = ci_95, col = "red", lwd = 2, lty = 2)