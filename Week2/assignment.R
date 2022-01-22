data <- read.csv("data.txt")
yt <- ts(data$y)
ts.plot(yt)
T=400
cat('Data contains ',T,' rows \n')

p=8
y=rev(yt[(p+1):T]) # response
X=t(matrix(yt[rev(rep((1:p),T-p)+rep((0:(T-p-1)),rep(p,T-p)))],p,T-p));
XtX=t(X)%*%X
XtX_inv=solve(XtX)
phi_MLE=XtX_inv%*%t(X)%*%y # MLE for phi
s2=sum((y - X%*%phi_MLE)^2)/(length(y) - p) #unbiased estimate for v

cat("\n MLE of conditional likelihood for phi: ", phi_MLE, "\n",
    "Estimate for v: ", s2, "\n")

n_sample=500 # posterior sample size
library(MASS)

## step 1: sample v from inverse gamma distribution
v_sample=1/rgamma(n_sample, (T-2*p)/2, sum((y-X%*%phi_MLE)^2)/2)

## step 2: sample phi conditional on v from normal distribution
phi_sample=matrix(0, nrow = n_sample, ncol = p)
for(i in 1:n_sample){
  phi_sample[i, ]=mvrnorm(1,phi_MLE,Sigma=v_sample[i]*XtX_inv)
}

## plot histogram of posterior samples of phi and nu
par(mfrow = c(3,3), cex.lab = 1.3)
for(i in 1:8){
  hist(phi_sample[, i],xlab=bquote('phi'~.(i)),main=bquote("Histogram of phi"~.(i)))
}
hist(v_sample,xlab=bquote('v'),main=bquote("Histogram of "~v))

roots=1/polyroot(c(1, -phi_MLE)) # compute reciprocal characteristic roots
r=Mod(roots)
# compute moduli of reciprocal roots
lambda=2*pi/Arg(roots) # compute periods of reciprocal roots
# print results modulus and frequency by decreasing order
print(cbind(r, abs(lambda))[order(r, decreasing=TRUE), ])
