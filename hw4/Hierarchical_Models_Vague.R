# y: The number of tumor-emerging rats in group;
# n: The number of total rats in group
rat.y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,
           2,1,5,2,5,3,2,7,7,3,3,2,9,10,4,4,4,4,4,4,
           4,10,4,4,4,5,11,12,5,5,6,5,6,6,6,6,16,15,15,9)
rat.n <- c(20,20,20,20,20,20,20,19,19,19,19,18,18,17,20,20,20,20,19,19,
           18,18,25,24,23,20,20,20,20,20,20,10,49,19,46,27,17,49,47,20,
           20,13,48,50,20,20,20,20,20,20,20,48,19,19,19,22,46,49,20,20,
           23,19,22,20,20,20,52,47,46,24)
logpost4 <- function(u, v) {a <- exp(u + v)/(exp(u) + 1);
b <- exp(v)/(exp(u) + 1);
J <- length(rat.y); # Number of data
x <- J * (lgamma(a + b)- lgamma(a)- lgamma(b)) + log(a) + log(b)-2.5 * log(a + b); # Constants, Jaccobi and hyperprior
for (i in (1:J)){
  x <- x + lgamma(a + rat.y[i]) + lgamma(b + rat.n[i] - rat.y[i]) - 
    lgamma(a + b + rat.n[i]) # Adding likelihood
}
x # Return log p(theta, phi | y)
}
# Generate grid sampling
x4 <- c(-230:-130)/100 # set the grid values of u
y4 <- c(100:500)/100 # set the grid values of v
# evaluate grid
z <- outer(x4,y4,logpost4) # evaluate the key part of the log posterior p(u,v|y) for each grid (u,v)
# exponentiate log posterior (set maximum to be one)
z2 <- exp(z- max(z)) # transfer log posterior into posterior (the joint posterior)

contour(x4,y4,z2,levels=c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95))