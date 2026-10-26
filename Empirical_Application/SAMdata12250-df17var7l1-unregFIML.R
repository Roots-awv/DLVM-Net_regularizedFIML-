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

###### Using optim()
# Optimize Full-Information ML using L-BFGS-B
lBFGSb.FIML <- optim(
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
    pen.l = 0,
    lower = lower.b,
    upper = upper.b,
    control = list(maxit = 500, factr = 1e7)
)
print(lBFGSb.FIML)
Est <- lBFGSb.FIML$par
write.csv(Est, file = "FIMLest-SAM-data.csv", row.names = FALSE)

# Numerical derivation
Hess.theta <- lBFGSb.FIML$hessian
write.csv(Hess.theta, file = "FIMLHessian-SAM-data.csv", row.names = FALSE)

# Other quantities
Neg.logLik<-lBFGSb.FIML$value
iter <-lBFGSb.FIML$counts
CVce<-lBFGSb.FIML$convergence
Quantities.FIML <- c(
      Neg.logLik,iter['gradient'],CVce)
names(Quantities.FIML)<-c("F_reg","gradient","Convergence")
write.csv(Quantities.FIML, file = "FIMLunreg-Quantities_from_optim-SAM-data.csv", row.names = FALSE)

