############################# Packages
library(MASS)
library(Matrix)
library(dplyr)

############################# Data preparation / Statistics summary 
# Computational methods
source("./Computing_Penalized_Pseudo_FIML_LargeBT_CovMatrix.R")

# Loading the data 
source("./revision2_dataprep1.R")

# Data management
y1.var <- apply(y1, c(3,1,2), t)
p.dim<-dim(y1.var)[1]
N.dim<-dim(y1.var)[2]
T.dim<-dim(y1.var)[3]

# Summary statistics
Z.summary<-summary.stats(y1.var)



####################### Optimization
# Combine into parameter vector
r.dim <- p.dim
m.dim <- 7

# Parameters initialization & boundaries
B<-diag(0.5*runif(m.dim,-1,1))
min.B <- rep(-1,m.dim^2) 
ini.B <- B%>%as.vector  
max.B <- rep(1,m.dim^2) 

#
min.ly <-vector("list", m.dim)
ini.ly <-vector("list", m.dim)
max.ly <-vector("list", m.dim)
for (m in c(1,m.dim-1,m.dim)){
  min.ly[[m]]<-c(1,-2,-2)
  ini.ly[[m]]<-rep(1,3)
  max.ly[[m]]<-c(1,2,2)
}

for (m in 2:(m.dim-2)){
  min.ly[[m]]<-c(1,-2)
  ini.ly[[m]]<-rep(1,2)
  max.ly[[m]]<-c(1,2)
}
min.Lambda <- min.ly %>% bdiag %>% as.vector
ini.Lambda <- ini.ly %>% bdiag %>% as.vector
max.Lambda <- max.ly %>% bdiag %>% as.vector

#
min.Sigeta <- matrix(-2,m.dim,m.dim)
diag(min.Sigeta) <- 0
ini.Sigeta <- (diag(m.dim)+1)
max.Sigeta <- matrix(2,m.dim,m.dim)

#
min.Sigeps <- matrix(0,r.dim,r.dim)
diag(min.Sigeps) <- 0
ini.Sigeps <- diag(r.dim)
max.Sigeps <- diag(rep(2,r.dim)) #matrix(0,r.dim,r.dim)

#
lower.b <- c(min.B,min.Lambda,
             min.Sigeta[lower.tri(min.Sigeta,diag=T)],
             min.Sigeps[lower.tri(min.Sigeps,diag=T)]#,
             # min.Sigmu[lower.tri(min.Sigmu,diag=T)]
)%>%as.vector

upper.b <- c(max.B,max.Lambda,
             max.Sigeta[lower.tri(max.Sigeta,diag=T)],
             max.Sigeps[lower.tri(max.Sigeps,diag=T)]#,
             # max.Sigmu[lower.tri(max.Sigmu,diag=T)]
)%>%as.vector

# 
theta0 <- c(ini.B*0,
            ini.Lambda,
            0.1*ini.Sigeta[lower.tri(ini.Sigeta,diag=T)],
            0.1*ini.Sigeps[lower.tri(ini.Sigeps,diag=T)]#,
            # 0.1*ini.Sigmu[lower.tri(ini.Sigmu,diag=T)]
)

# Adaptive/standard Lasso
FIML.est<-read.csv("FIMLest-SAM-data.csv")
est.inv <- 1/FIML.est[1:length(theta0),'x']
inf.vals <- ifelse(est.inv == Inf,0,1)
weights <- inf.vals*est.inv
weights[is.na(weights)]<-max(weights,na.rm=T)
lambda.w <- abs(weights)/sqrt(N.dim*T.dim)
lambda.w <- lambda.w%>%as.vector
print(round(lambda.w,2))

###### Using optim()
# Optimize Full-Information ML using L-BFGS-B
lBFGSb.penFIML <- optim(
  theta0, 
  fn=  log.penFIMLit.Y, 
  gr = ScoreY.penFIMLit, 
  method = "L-BFGS-B",
  hessian = TRUE,
  data = Z.summary,  
  m.dim = m.dim,
  r.dim = r.dim,
  N.dim = N.dim,
  T.dim = T.dim,
  pen.l = lambda.w,
  lower = lower.b,
  upper = upper.b
)
print(lBFGSb.penFIML)
Est <- lBFGSb.penFIML$par
write.csv(Est, file = "alasso1-FIMLest-SAM-data.csv", row.names = FALSE)

# Numerical derivation
Hess.theta <- lBFGSb.penFIML$hessian
write.csv(Hess.theta, file = "alasso1-FIMLHessian-SAM-data.csv", row.names = FALSE)

# Other quantities
Neg.logLik<-lBFGSb.penFIML$value
iter <-lBFGSb.penFIML$counts
CVce<-lBFGSb.penFIML$convergence
Quantities.penFIML <- c(
      Neg.logLik,iter['gradient'],CVce)
names(Quantities.penFIML)<-c("F_reg","gradient","Convergence")
write.csv(Quantities.penFIML, file = "alasso1-Quantities_from_optim-SAM-data.csv", row.names = FALSE)

