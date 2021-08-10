library(TMB)
library(Matrix)
compile("templates/cc.cpp")
dyn.load(dynlib("templates/cc"))
set.seed(123)

n <- 20
case_day = 2:n
control_days = cbind(1:(n-1),0:(n-2))
count <- rpois(length(case_day),10)

d0 <- 3
d1 <- 15
d2 <- 12


X <- matrix(rnorm(n*d0),n,d0)
A1 <- matrix(0,n,d1)
A1[cbind(1:n,sample(d1,n,T))] <- 1
A1 <- A1[,-(5:7)]
A2 <- matrix(0,n,d2)
A2[cbind(1:n,sample(d2,n,T))] <- 1
A2 <- A2[,-(5:7)]
A <- cbind(A1,A2)
A <- as(A,"dgTMatrix")
createD <-function(d,p){
  if(p==0) return(Diagonal(d,1))
  D <- Matrix::bandSparse(d,k =c(0,1),diagonals =list(rep(-1,d),rep(1,d-1)))[-d, ]
  if(p==1) return(D)
  else return(createD(d,p-1)[-1,-1] %*% D)
}
# Q1 <- as(diag(d1),"dgTMatrix")
# Q2 <- as(diag(d2),"dgTMatrix")
Q1 <- crossprod(createD(d1,3)[,-(5:7)])
Q2 <- crossprod(createD(d2,3)[,-(5:7)])
Q <- bdiag(Q1,Q2)
gamma_dims <- c(ncol(Q1),ncol(Q2))

beta_prec <- c(.001,.005,.01)
theta_prior_id <- c(1,2,2)
theta_hypers <- c(.01,.5,
                  1.1,2.1,
                  1.1,2.1)

beta <- rep(0,ncol(X))
gamma <- rep(0,ncol(A))
z <- rep(0,n)
theta <- c(1,2,3)

data <- list(count = count, case_day = case_day, control_days = control_days,
             X = X, A = A, Q = Q, gamma_dims = gamma_dims,
             beta_prec = beta_prec, theta_prior_id = theta_prior_id, theta_hypers = theta_hypers)


length(count); length(case_day); dim(control_days)
dim(X); dim(A); dim(Q); gamma_dims
length(beta_prec); length(theta_prior_id); length(theta_hypers)

parameters <- list(beta = beta, gamma = gamma, z = z, theta = theta)

obj <- MakeADFun(data, parameters, random = c("beta","gamma","z"), DLL="cc", hessian=T)

obj$fn(c(1,2,3))

opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
sdreport(obj)
