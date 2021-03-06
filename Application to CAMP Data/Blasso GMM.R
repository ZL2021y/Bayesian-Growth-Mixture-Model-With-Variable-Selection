#rm(list = ls())
library(nlme)
library(lme4)
library(ggplot2)
library(mvtnorm)
library(MCMCpack)
library(mclust)
library(cluster) 
library(statmod)		# for invgauss
library(abind)	# abind()
library(label.switching)
 
bayes.lasso <- function(	mydat,				# trajectory data
			    	covar,				# covariate matrix
			    	num.cluster,			# num of clusters
				K = 2, 				# dimention of random effects: 1 - intercept; 2-slope for linear term; 3-slope for quadratic terms
			    	burn.in=1000, 
			    	thin=5, 
			    	per = 100, 
				max.iter = 2000, 
				lambda.prior = list(a.lambda=NULL,b.lambda=NULL),  
				myseed=2020){
set.seed(myseed)
y <- mydat$y
xt <- mydat$time 
N <- length(unique(mydat$id))
n.obs <- as.numeric(table(mydat$id))
mydat$ID <- mydat$id
id <- mydat$id
mydat <- mydat[order(mydat$ID,mydat$time),]  # Sort by ID and time
mydat_uq <- mydat[!duplicated(mydat$id, fromLast=TRUE),] # Keep last observation per ID
p <- dim(covar)[2] ; p		# total number of covariates
#--------------------------------------------------------------#
# starting values;
#--------------------------------------------------------------#
fit.lme <-lme(y ~ time +  I(time^2), random = ~ time| ID, data=mydat, control = list(maxIter = 500),method="ML");fit.lme 
# fit.lme <-lme(y ~ AGE +  I(AGE^2), random = ~ AGE| ID, data=mydat, control = list(maxIter = 500),method="ML");
# fixed.effect <- apply(coef(fit.lme),2,mean); fixed.effect
# random.effect <- coef(fit.lme) - matrix(fixed.effect,nrow=N,ncol=3,byrow=TRUE) ; random.effect # Need to use "byrow" options here
# cov(random.effect[,1:2])
# str(fit.lme)

fixed.effect <- apply(coef(fit.lme),2,mean); fixed.effect
random.effect <- coef(fit.lme) - matrix(fixed.effect,nrow=N,ncol=3,byrow=TRUE) ; random.effect # Need to use "byrow" options here

# starting value for cluster membership
fit.clust <- pam(random.effect,num.cluster)
cluster <- fit.clust$clustering; table(cluster)
mydf.clust <- data.frame(ID=1:length(cluster),cluster)
dat <- merge(mydat,mydf.clust,by="ID")
#------------------------------------------------------#
# starting value for regression coefficients
fixed.effect <- NULL
for (j in 1:num.cluster){
fit.lme <- lme(y ~  time +  I(time^2), random = ~ 1|ID, data=dat[dat$cluster==j,] )
fixed.effect <- rbind(fixed.effect,apply(coef(fit.lme),2,mean))
}
theta <- 0.01*as.matrix(random.effect) # intial values for random effect (all subjects): random intercept and slope

# avoid using "pi", which conflicts with the function sampling from "PG" distribution;
# instead use "ppi" here to represent the proportion between the two groups;
library(robustHD)
library(nnet)

starts <- list(beta0 = rep(0,num.cluster-1),beta = matrix(0,ncol=p,nrow=num.cluster-1),ga=fixed.effect,sigma.sq.e=rep(2,num.cluster))
beta0 <- c(starts$beta0); beta <- matrix(starts$beta,nrow=num.cluster-1)
zz <- cluster; table(zz)
# regression coeffcients
ga <- starts$ga ;ga 
# for residual variance
sigma.sq.e <- starts$sigma.sq.e  ; sigma.sq.e

# RANDOM EFFECT
sigma.sq.u <- sigma.sq.u.inv <-array(0,c(K,K,num.cluster))
for (j in 1:num.cluster){sigma.sq.u[,,j] <-  var(random.effect[cluster==j,1:K]);
				sigma.sq.u.inv[,,j] <- solve(sigma.sq.u[,,j]) }	 
#--------------------------------------------------------------#
# prior values (hyper-parameters)
#--------------------------------------------------------------#
#---- residual variables
a0 <- rep(3,num.cluster)
b0 <- rep(1e-2,num.cluster)

#----- Fixed Effect Variables
w0 <- matrix(0,nrow=num.cluster,ncol=dim(ga)[2])
omega0 <- array(diag(1e3,dim(ga)[2]),dim = c(dim(ga)[2],dim(ga)[2],num.cluster))

# hyper-parameters for Lasso
# prior for lambda2 
library(glmnet)
xx <- standardize(as.matrix(covar))
out.zz <- cluster
out.zz <- relevel(factor(out.zz), ref = max(out.zz))  # the last categories as reference category
set.seed(2020)
cvfit <- cv.glmnet(x = xx, y = out.zz, family = "multinomial");cvfit$lambda.min 

# According to park and Casella (2008): the prior density for lambda2 should approach 0 sufficiently fast as lambda2 -> \infty (to avoid mixing problems),
# but should be relatively flat and place high probablity near the maximum likelihood estimate
# Choosing inappropriate prior here will have mixing problem
#-------------------------------------------------------------
if(is.null(lambda.prior$a.lambda)==TRUE & is.null(lambda.prior$b.lambda)==TRUE){
	a.lambda <- 10
	b.lambda <- 1/cvfit$lambda.min} else{
	a.lambda <- lambda.prior$a.lambda
	b.lambda <- lambda.prior$b.lambda}
#---- Random effect variables
r0 <-  rep(4,num.cluster)
R0 <-  array(0,c(K,K,num.cluster))
for (j in 1:num.cluster){R0[,,j] <- diag(3,K)}

# initialize
lambda2 <- rep(0.01,num.cluster-1)
tau2 <- matrix(1,ncol=p,nrow=num.cluster - 1)
#--------------------------------------------------------------#
# Storing the sample;
ZZ <- NULL
GA <- NULL
SIGMA.SQ.E <- NULL
SIGMA.SQ.U <- NULL
BETA <- BETA0 <- NULL
T <- NULL

cat(paste(rep('-',60),sep='',collapse=''), '\n'); cat(paste(rep('-',60),sep='',collapse=''), '\n'); 
cat('Running LASSO Joint Modeling Method', '\n')
cat(paste(rep('-',60),sep='',collapse=''), '\n'); cat(paste(rep('-',60),sep='',collapse=''), '\n'); 

iter <- 0
begin = proc.time()[1]
repeat
  {
#----------------------------------------------#
# Analyzing the data with variable selection		
#----------------------------------------------#
dm <- matrix(model.matrix(~factor(zz))[,-1],ncol=num.cluster-1); colnames(dm)<- NULL
kappa <- dm-0.5		# first-category as the reference category
#-----------------------------------------------#
covar <- cbind(as.matrix(standardize(covar)))		# NOT INCLUDED INTERCEPT	 
xi.tmp <- matrix(0,ncol=num.cluster-1,nrow=N)
for (j in 1:(num.cluster-1)){xi.tmp[,j] <- beta0[j] + covar %*% beta[j,]}
w <- NULL
for (j in 1:(num.cluster-1)){
	w <- cbind(w,rpg.devroye.R(num=N, n=rep(1,N), z=xi.tmp[,j]- log(apply(exp(as.matrix(xi.tmp[,-j])),1,sum) + 1 ))$x)}
C <- NULL
C.sum <- 0
for (j in 1:(num.cluster-1)){C.sum <- C.sum + exp(beta0[j] + covar %*% beta[j,])}
for (j in 1:(num.cluster-1)){C <- cbind(C,log(C.sum- exp(beta0[j] + covar %*% beta[j,]) + 1 )) }

for (j in 1:(num.cluster-1)){
#--------------------------------------------------------------#
# Sample beta (WITHOUT THE INTERCEPT)
#--------------------------------------------------------------#
	Omega <- diag(w[,j]) 
	eta <- kappa[,j]/w[,j]
	D <- diag(tau2[j,])
	V <- solve(t(covar) %*% Omega %*% covar + solve(D,tol=1e-20))
	m <- V %*% (t(covar) %*% Omega %*% matrix(eta - beta0[j],nrow=N)) +  V %*% t(covar) %*% Omega %*% C[,j]
	beta[j,] <- rmvnorm(n = 1, mean= m, sigma =V)

#--------------------------------------------------------------#
# Sample beta0 (THE INTERCEPT)
#--------------------------------------------------------------#
	s <- sum(w[,j])
	nu <- as.vector(eta - covar %*% beta[j,]  + C[,j])*w[,j]
	beta0[j] <- rnorm(1,mean=sum(nu)/s,sd=sqrt(1/s)) }
#--------------------------------------------------------------#
# sample lambda2
#--------------------------------------------------------------#
for (j in 1:num.cluster-1) {lambda2[j] <- rgamma(1,shape=p+a.lambda,rate=sum(tau2)/2 + b.lambda)}
#--------------------------------------------------------------#
# sample tau2
#--------------------------------------------------------------#
eval <- function(m) class(try(solve(m, tol = 1e-20),silent=T))=="matrix"	# function to evaluate singularity issue
for (j in 1:(num.cluster-1)){
	repeat{
	for (s in 1:p) {tau2[j,s] <- 1/rinvgauss(n=1, mean=sqrt((lambda2[j])/beta[j,s]^2), shape=lambda2[j])}
		if (eval(diag(tau2[j,]))==TRUE) break}	# if invertiable, then break
	}	# lambda2 denote lambda square
#--------------------------------------------------------------#
# Sample zz_i
#--------------------------------------------------------------#
t <- matrix(0,ncol=num.cluster,nrow=N)
for (i in 1:N) {
deno <-  1 + sum(exp(beta0 + apply(beta,1,function(x) covar[i,]%*%x)))
ppi <- c(exp(beta0 +apply(beta,1,function(x) covar[i,]%*%x)),1)/deno
m <- cbind(1,xt[which(id==i)],I(xt[which(id==i)]^2))	# construct the basis
g <- matrix(apply(ga,1,function(x) x %*% t(m) + theta[i,] %*% t(m)),ncol=num.cluster)	# this is n_i times num.cluster dimention
f <- NULL
for (j in 1:num.cluster){f <- c(f,dmvnorm(t(y[which(id==i)]), mean=g[,j],sigma=diag(as.numeric(sigma.sq.e[j]),n.obs[i] )))}
t[i,] <- (ppi*f)/sum(ppi*f)
zz[i] <- which.is.max(t[i,])
}
#--------------------------------------------------------------#
# Sample sigmas (residual variances)
#--------------------------------------------------------------#
a0.tid <-matrix((2*a0 +  aggregate(n.obs,by=list(zz),FUN=sum)[,2])/2,ncol=num.cluster)
tmp <- rep(0,num.cluster)
for (i in 1:N) {
m <- cbind(1,xt[which(id==i)],I(xt[which(id==i)]^2))	# construct the basis
g <- matrix(apply(ga,1,function(x) x %*% t(m) + theta[i,] %*% t(m) ),ncol=num.cluster)	# this is n_i times num.cluster dimention
	for (j in 1:num.cluster){tmp[j] <- tmp[j] + (zz[i]==j)*norm(as.matrix(y[which(id==i)] - g[,j]),"f")^2}
}
b0.tid <- matrix(b0 + tmp/2,ncol=num.cluster)

# https://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters
sigma.sq.e <-  mapply(function(x,y) rinvgamma(1,shape=x,scale=y),a0.tid,b0.tid); sigma.sq.e 

#--------------------------------------------------------------#
# Sample Fixed Effect Coefficients
#--------------------------------------------------------------#
tmp <-  array(0,c(dim(ga)[2],dim(ga)[2],num.cluster))
tmpp <-  matrix(0,ncol=dim(ga)[2],nrow=num.cluster)
for (i in 1:N) {
m <- cbind(1,xt[which(id==i)],I(xt[which(id==i)]^2))	# construct the basis
  for(j in 1: num.cluster){
	tmp[,,j] <- tmp[,,j] + (zz[i]==j)*(t(m) %*% m) * (1/sigma.sq.e[j])
	tmpp[j,] <- tmpp[j,] + (zz[i]==j)*(t(m) %*% matrix(y[which(id==i)] - theta[i,] %*% t(m),ncol=1))* (1/sigma.sq.e[j])}
}
omega0.tid <-  array(0,c(dim(ga)[2],dim(ga)[2],num.cluster))
w0.tid  <-  matrix(0,ncol=dim(ga)[2],nrow=num.cluster)
for(j in 1: num.cluster){
	omega0.tid[,,j] <-  solve(solve(omega0[,,j]) + tmp[,,j])
	w0.tid[j,] <- omega0.tid[,,j] %*% (tmpp[j,])
	ga[j,] <- mvrnorm(n = 1, mu = w0.tid[j,], Sigma = omega0.tid[,,j])}
#--------------------------------------------------------------#
# Reorder the fixed effect coefficients
#--------------------------------------------------------------#
pred <- NULL
max.t <- max(mydat$time)
m <- cbind(1,max.t,max.t^2)	# construct the basis
for (j in 1:num.cluster) {
		# predicted mean curve
		pred.mean <- ga[j,]%*% t(m); pred <- rbind(pred,pred.mean)}
ga <- ga[order(pred,decreasing = TRUE),]
#--------------------------------------------------------------#
# Reorder the fixed effect coefficients
#--------------------------------------------------------------#
pred <- NULL
max.t <- max(mydat$time)
m <- cbind(1,max.t,max.t^2)	# construct the basis
for (j in 1:num.cluster) {
		# predicted mean curve
		pred.mean <- ga[j,]%*% t(m); pred <- rbind(pred,pred.mean)}
ga <- ga[order(pred,decreasing = TRUE),]
#--------------------------------------------------------------#
# Sample Sigma's (Random effect variances)
#--------------------------------------------------------------#
tmp <- sigma.sq.u.inv <- sigma.sq.u <-array(0,c(K,K,num.cluster))
for (i in 1:N) {
	for (j in 1:num.cluster){
		tmp[,,j] <- tmp[,,j] + (zz[i]==j)* theta[i,1:K] %*% t(theta[i,1:K])}}
R0.tid <- r0*R0 + tmp
for (j in 1:num.cluster){sigma.sq.u.inv[,,j] <-  rWishart(1,sum(zz==j) + r0[j], R0.tid[,,j])[,,1];
				 sigma.sq.u[,,j] <- solve(sigma.sq.u.inv[,,j])}

#--------------------------------------------------------------#
# Sample random effects via Multivariate Normal  distribution
#--------------------------------------------------------------#
tmp <-  array(0,c(K,K,num.cluster) )
tmpp <-  matrix(0,ncol=K,nrow=num.cluster)
for (i in 1:N){
	m <- cbind(1,xt[which(id==i)],I(xt[which(id==i)]^2))	# construct the basis
	mz <- cbind(1,xt[which(id==i)])	# construct the basis
	for (j in 1:num.cluster){
		tmp[,,j] <- tmp[,,j]  + (zz[i]==j)* (t(mz) %*% mz)*(1/sigma.sq.e[j]) 
		tmpp[j,] <- tmpp[j,] +  (zz[i]==j)*(t(mz) %*% matrix(y[which(id==i)] - ga[j,]%*%t(m),ncol=1))* (1/sigma.sq.e[j])}}
for (j in 1:num.cluster){
	sigma.tid <- solve(solve(sigma.sq.u[,,j]) + tmp[,,j])
	mu.tid <- sigma.tid %*% tmpp[j,]
	theta[zz==j,1:K] <- mvrnorm(n = sum(zz==j), mu = mu.tid, Sigma = sigma.tid)
}


# storing the sample;
if (iter >= burn.in & iter %% thin == 0){
SIGMA.SQ.E <- rbind(SIGMA.SQ.E,sigma.sq.e)
SIGMA.SQ.U <- abind(SIGMA.SQ.U,apply(sigma.sq.u,3,vech),along = 3)
T  <- abind(T,t,along = 3)
ZZ <- rbind(ZZ,zz)
GA  <- abind(GA,ga,along = 3)
BETA <- abind(BETA,beta,along = 3)
BETA0 <- rbind(BETA0,beta0)}
iter <- iter + 1

if(iter %% per == 0) {cat('iter =',iter,'\n');print(lambda2);print(beta0);print(beta)}
if(iter == max.iter) break
}
end = proc.time()[1]
cat('It took', end - begin, 'seconds\n')
time <-  end - begin
#----------------------------------------------------------------------------#
num.sample <- length(seq(burn.in + 1,iter,thin))
#----------------------------------------------------------------------------#
# Address Label Switching Using Stephens' algorithm
#----------------------------------------------------------------------------#
T.trans <- array(0,c(num.sample,N,num.cluster))
for (j in 1:num.cluster){ T.trans[,,j] <- t(T[,j,])}
out.relabel <- label.switching(method="STEPHENS",z=ZZ,K=num.cluster,p=T.trans)
mylabel <- out.relabel$permutations$STEPHENS
for (i in 1:num.sample){
	SIGMA.SQ.E[i,] <- SIGMA.SQ.E[i,mylabel[i,]]
	SIGMA.SQ.U[,,i] <- SIGMA.SQ.U[,mylabel[i,],i]
	GA[,,i] <- GA[mylabel[i,],,i]
	BETA[,,i] <- BETA[mylabel[i,-which.is.max(mylabel[i,])],,i]
	T[,,i] <- T[,mylabel[i,],i]}
cluster.marginal <- apply(apply(T,c(1,2),mean),1,which.is.max)
postprob <- apply(apply(T,c(1,2),mean),1,max)

# returning the final parameters;
list(GA=GA, SIGMA.SQ.E=SIGMA.SQ.E, SIGMA.SQ.U=SIGMA.SQ.U,BETA0=BETA0,BETA=BETA, cluster.marginal=cluster.marginal,postprob=postprob)
}
 