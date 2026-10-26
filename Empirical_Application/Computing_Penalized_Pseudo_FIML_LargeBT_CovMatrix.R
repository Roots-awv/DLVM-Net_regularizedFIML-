###########################################################
######################### Computational Methods / Functions
###########################################################

################### Loading packages
library(MASS)
library(Matrix)
library(dplyr)


################### Transformation matrix/vector 

# Identity matrix's vector column
e.vec <- function(m,n){
  e <- Matrix(rep(0,m),sparse=TRUE)
  e[n,1] <- 1
  return(e)
}

# Commutation matrim
K.mn <- function(m,n){
  K <- c()
  I.m<-Matrix(diag(m),sparse=TRUE)
  for (i in 1:n){
    K <- cbind(K,kronecker(I.m, e.vec(n,i)))
  }
  K.sparse<-Matrix(K,sparse=T)
  return(K.sparse)
}

# Duplication matrix
D.n <- function(n,sparse=T){
  D <- matrix(0,n^2,n*(n+1)*0.5)
  for(i in 1:n){
    for(j in 1:i){
      index.1 <- (j-1)*n+i
      index.2 <- (i-1)*n+j
      p.ij <- (n-0.5*j)*(j-1)+i
      delta.ineqj <- ifelse(i!=j,1,0) 
      D <- D+(e.vec(n^2,index.1)+
                e.vec(n^2,index.2)*delta.ineqj)%*%t(e.vec(n*(n+1)*0.5,p.ij))
    } 
  }
  D.sparse<-Matrix(D,sparse=T)
  return(D.sparse)
}

################### Parameters

# Parameter matrices
Param.ls<-function(x,m.dim,r.dim){
  # x<-theta
  B <- matrix(x[1:(m.dim^2)], m.dim, m.dim)%>%Matrix(sparse = TRUE)
  # 
  L.start <- (m.dim^2)+1
  L.end <- m.dim^2+m.dim*r.dim
  Lambda <- matrix(x[L.start:L.end], r.dim, m.dim)%>%Matrix(sparse = TRUE)
  # 
  Seta.start <- L.end + 1
  Seta.end <- L.end + m.dim * (m.dim + 1) * 0.5
  Seta.vec <- D.n(m.dim) %*% x[Seta.start:Seta.end]
  Sig.zeta <- matrix(Seta.vec, m.dim, m.dim)%>%Matrix(sparse = TRUE)
  #
  Seps.start<-Seta.end+1
  Seps.end<-Seta.end+r.dim*(r.dim+1)*0.5
  Seps.vec <- D.n(r.dim)%*%x[Seps.start:Seps.end]
  Sig.eps <- matrix(Seps.vec, r.dim, r.dim)%>%Matrix(sparse = TRUE)
  #
  # Smu.start <- Seps.end+1
  # Smu.end <- Seps.end+m.dim*(m.dim+1)*0.5
  # Smu.vec <- D.n(m.dim) %*% x[Smu.start:Smu.end]
  # Sig.mu <- matrix(Smu.vec,m.dim,m.dim)%>%Matrix(sparse = TRUE)
  #  
  ls.out<-list(
    B=B,Lambda=Lambda,Sig.zeta=Sig.zeta,Sig.eps=Sig.eps)#,Sig.mu=Sig.mu)
  return(ls.out)
}

##################### Covariance matrix

# Covariance matrix structure
Vcov.matrix <- function(Sig) {
  G.11<-G.12<-G.21<-G.22<-matrix(0,2,2)
  G.11[1,1]<-G.12[1,2]<-
    G.21[2,1]<-G.22[2,2]<-1
  Sig.m11 <- kronecker(G.11,Sig[[1]]) 
  Sig.m12 <- kronecker(G.12,Sig[[2]]) 
  Sig.m21 <- kronecker(G.21,Sig[[2]]) 
  Sig.m22 <- kronecker(G.22,Sig[[1]]) 
  Sig.m <- Sig.m11+Sig.m12+Sig.m21+Sig.m22 
  return(Sig.m)
}

# Splitting variance matrix
HxI.item.ij<-function(a,b,m){
  #a<-2;b<-1
  G.ij<-matrix(0,2,2)
  G.ij[a,b]<-1
  vec.G<-G.ij%>%as.vector
  H.term1<-kronecker(diag(2),K.mn(m,2))
  H.term2<-kronecker(vec.G,diag(m))
  H.term3<-kronecker(H.term1%*%H.term2,diag(m))
  return(H.term3)
}

# Latent variable covaraince matrix
Sigma.eta<-function(C,S,m){
  #
  inv.IminBxB<-solve(diag(m^2)-kronecker(C,C))
  vecSig.eta.zero<-inv.IminBxB%*%as.vector(S)  
  out.var<-matrix(vecSig.eta.zero,m,m)%>%Matrix(sparse = TRUE)
  # 
  return(out.var)
}

# Observed items covaraince matrix
Sigma.y<-function(x,m.dim,r.dim){
  # x<-theta0
  par.ls<-Param.ls(x,m.dim,r.dim)
  B <- par.ls$B
  Lambda <- par.ls$Lambda
  Sig.zeta <- par.ls$Sig.zeta
  #Sig.mu<-par.ls$Sig.mu
  Sig.eps <- par.ls$Sig.eps
  #
  Sig.eta.z<-Sigma.eta(B,Sig.zeta,m.dim)
  Sig.eta.zero<-Sig.eta.z#+Sig.mu
  Sig.eta.one<-B%*%Sig.eta.zero
  #
  Sig.y.zero<-Lambda%*%(Sig.eta.zero)%*%t(Lambda)+Sig.eps 
  Sig.y.one<-Lambda%*%Sig.eta.one%*%t(Lambda)
  #
  #diag(Lambda%*%(Sig.eta.zero-Sig.mu)%*%t(Lambda))/diag(Sig.y.zero)
  #
  out.list<-list(lag0=Sig.y.zero,lag1=Sig.y.one)
  return(out.list)
}

# Estimation of the Between level variation
Mu_u.hat <- function(invS,M,bar.Z){
  # Compute Eigen Decomposition
  MSM<-t(M)%*%invS%*%M
  # Compute Sigma Inverse
  eig_full.MSM <- eigen(MSM, symmetric = TRUE)
  eig_vals.MSM <- pmax(eig_full.MSM$values, 1e-4)
  inv.MSM <- eig_full.MSM$vectors%*%diag(1/eig_vals.MSM)%*%t(eig_full.MSM$vectors)
  Mu.hat<-inv.MSM%*%(t(M)%*%invS%*%bar.Z)
  #mu.hat<-Mu.hat%>%as.matrix
  return(Mu.hat%>%as.matrix)
}

Sigma_u.hat<-function(x,mean.T.dim,m.dim,r.dim,pen){
  # pen<-lambda.w[1]
  # x<-penFIML.est.lBFGSb
  # Initial parameters
  Sigma.y.ls <- Sigma.y(x,m.dim,r.dim)
  Sigma.hat <- rbind(cbind(Sigma.y.ls$lag0,Sigma.y.ls$lag1),
                     cbind(Sigma.y.ls$lag1,Sigma.y.ls$lag0))
  Sigma.hat <-(Sigma.hat+Sigma.hat)/2
  # Parameter estimates
  theta.hat<-x
  Params.theta<-Param.ls(theta.hat,m.dim,r.dim)
  Lambda.hat<-Params.theta$Lambda
  Sigeps.hat<-Params.theta$Sig.eps
  B.hat<-Params.theta$B
  # Compute Eigen Decomposition
  eig_full.hat <- eigen(Sigma.hat, symmetric = TRUE)
  eig_vals.hat <- pmax(eig_full.hat$values, 1e-4)  
  # Compute Sigma Inverse
  inv.Sigma.hat <- eig_full.hat$vectors%*%diag(1/eig_vals.hat)%*%t(eig_full.hat$vectors)
  M <- Matrix(rep(1,2)%>%kronecker(diag(r.dim)),sparse = TRUE)  
  # Compute Eigen Decomposition
  MSM<-t(M)%*%inv.Sigma.hat%*%M
  eig_full.MSM <- eigen(MSM, symmetric = TRUE)
  eig_vals.MSM <- pmax(eig_full.MSM$values, 1e-4)  
  # Compute Sigma Inverse
  inv.MSM <- eig_full.MSM$vectors%*%diag(1/eig_vals.MSM)%*%t(eig_full.MSM$vectors)
  # Random intercept
  Mu.hat<-inv.MSM%*%(t(M)%*%inv.Sigma.hat%*%mean.T.dim-pen)
  mu.hat<-Mu.hat%>%as.matrix
  mu.est<-solve(t(Lambda.hat)%*%(Lambda.hat))%*%t(Lambda.hat)%*%(mu.hat)
  # Between-person variances
  cov(t(as.matrix(mu.est)))#*(N.dim-1)/(N.dim)
  dim(mu.est)
  dim(mu.hat)
  P.Lambda <- solve(t(Lambda.hat)%*%Lambda.hat)%*%t(Lambda.hat)
  cov_mu <- P.Lambda%*%(cov(t(mu.hat)))%*%t(P.Lambda)
  (diag(m.dim)-B.hat)%*%cov_mu%*%(diag(m.dim)-t(B.hat))
  cov_mu-B.hat%*%cov_mu%*%t(B.hat)
  cov_mu/2
  return(cov_mu)
}
    
##### Full information maximum likelihood
log.penFIMLit.Y <- function(x, data, m.dim, r.dim, N.dim, T.dim, pen.l) {
  # x<-theta
  # data<-Z.summary
  # 
  # pen.l <-0 lambda.w
  mean.Y <- data[['mean']]
  Y.l0 <- data[['Y.l0']]
  Y.l1 <- data[['Y.l1']]
  Svcov.Z <- data[['VcoV']]
  # Compute Sigma.x 
  Sigma.ls <- Sigma.y(x, m.dim, r.dim)
  Sigma.x <- Vcov.matrix(Sigma.ls)
  # Sigma.x <- (Sigma.x + t(Sigma.x)) / 2
  # Sigma.x <- Sigma.x + diag(1e-6, nrow(Sigma.x))
  # Compute Eigen Decomposition
  eig_full <- eigen(Sigma.x, symmetric = TRUE)
  eig_vals <- pmax(eig_full$values, 1e-4)
  # Compute Sigma Inverse
  Sigma.x <- eig_full$vectors%*%diag(eig_vals)%*%t(eig_full$vectors)
  L.x <- chol(Sigma.x)
  Sigma.invx <- chol2inv(L.x)
  Linv.x<-chol(Sigma.invx)
  # Projection Matrix
  M <- Matrix(rep(1,2)%>%kronecker(diag(r.dim)),sparse = TRUE)
  P <- diag(r.dim*2) - M%*%solve(t(M)%*%M)%*%t(M)
  tY <- P%*%mean.Y
  # Observed data
  outer_full <- tcrossprod(tY)
  trace_term_matrix <- Sigma.invx%*%(Svcov.Z+outer_full)
  # Log 
  log.det.SS <- log(diag(Linv.x))
  # Penalty
  # Sigma.invdiag <- diag(Sigma.invx)
  # Sigma.invpen <- Sigma.invx - diag(Sigma.invdiag)
  penalty.term <- sum(abs(pen.l*x)) #*Sigma.invpen
  #
  LL.y <- 0
  for (n in 1:N.dim) {
    # 
    Y.l0.n <- Y.l0[, n, ]
    Y.l1.n <- Y.l1[, n, ]
    Y.l0l1.n <- rbind(Y.l0.n, Y.l1.n)
    # tr.terms <-0
    # log.terms <-0
    for (t in 1:T.dim) {
      #
      Y.l0l1.nt <- Y.l0l1.n[, t]
      obs_idx <- which(!is.na(Y.l0l1.nt))
      if (length(obs_idx) == 0) next
      tr.m <- trace_term_matrix[obs_idx, obs_idx]
      tr.terms <- sum(diag(tr.m))
      # pen.term <- sum(abs(penalty.term[obs_idx, obs_idx]))
      # Negative log-Lik / fit function
      log.v<-log.det.SS[obs_idx]
      log.terms <- 2*sum(log.v)
      LL.y <- LL.y-N.dim*log.terms+tr.terms#+pen.term
    }
  }
  fiml <- LL.y/(N.dim*T.dim)
  reg.fiml <- fiml+penalty.term
  reg.fiml
  return(reg.fiml)
}

#################### Score function
dvecSy.dTheta <- function(dvec.S0,dvec.S1,o.dim){
  H11xI<-HxI.item.ij(1,1,o.dim)
  H12xI<-HxI.item.ij(1,2,o.dim)
  H21xI<-HxI.item.ij(2,1,o.dim)
  H22xI<-HxI.item.ij(2,2,o.dim)
  #
  dvec.S.out<-(H11xI+H22xI)%*%dvec.S0+
    (H12xI+H21xI)%*%dvec.S1
  return(dvec.S.out)
}

### Covariance derivatives 
# WRT B
dLy.dB <- function(x,m.dim,r.dim){
  # x<-theta
  # data<-Z.summary
  #
  par.ls<-Param.ls(x,m.dim,r.dim)
  B<-par.ls$B
  Lambda<-par.ls$Lambda
  Sig.zeta<-par.ls$Sig.zeta
  # Sig.mu<-par.ls$Sig.mu
  Sig.eps<-par.ls$Sig.eps
  #
  # Sigma.ls <- Sigma.y(x,m.dim,r.dim)
  Sigma.eta.z<-Sigma.eta(B,Sig.zeta,m.dim)
  Sigma.eta.zero<-Sigma.eta.z#+Sig.mu
  Sigma.eta.one<-B%*%Sigma.eta.zero
  #
  inv.IminBB <- solve(diag(m.dim)-B%*%t(B))
  # inv.IminBB.B<-inv.IminBB%*%B
  # S0.inv.IminBB<-(Sig.zeta)%*%inv.IminBB
  # inv.IminBB.B.S0<-inv.IminBB%*%B%*%(Sig.zeta)
  kron.Szeta.I<-kronecker(Sig.zeta,diag(m.dim))
  kron2.inv.IminBB<-kronecker(inv.IminBB,inv.IminBB)
  #kron.tB.I<-kronecker(2*t(B),diag(m.dim)) 
  kron.B.I<-kronecker((B),diag(m.dim)) 
  kron.I.B<-kronecker(diag(m.dim),(B)) 
  Kkk<-K.mn(m.dim,m.dim)
  #IxKkk<-kronecker(diag(m.dim),Kkk)
  #
  #term1<-kronecker(inv.IminBB.B.S0,inv.IminBB)
  #term2<-kronecker(S0.inv.IminBB,inv.IminBB.B)
  #dSS.eta0<-((term1)+(term2))
  term1<-kron.Szeta.I%*%(kron2.inv.IminBB)
  dSS.eta0<-((term1)%*%(kron.I.B+Kkk%*%kron.I.B%*%Kkk))
  term3<-kronecker(t(Lambda),t(Lambda))
  dSS0.dB <-dSS.eta0%*%term3
  # B%*%Sigma.eta.zero%*%t(B)+Sig.zeta
  term4 <-kronecker(Sigma.eta.zero,diag(m.dim))
  term5 <-kronecker(t(B),diag(m.dim))
  term6 <- (term4+term5%*%dSS.eta0)
  dSS1.dB <- (term6)%*%(term3)
  #
  dSS.dB <- dvecSy.dTheta(t(dSS0.dB),t(dSS1.dB),r.dim)
  # Score 
  #matrix(c(0.75,-1.5,1.5,-0.75),m.dim,m.dim)
  # dL1<-dpL1.terms(x,data,m.dim,r.dim,N.dim,T.dim)
  # t(dSS.dB)%*%dl.term%>%as.vector
  # grad_num[1:4]
  return(t(dSS.dB))
}

# WRT Lambda
dLy.dLambda <- function(x,m.dim,r.dim){
  # x<-theta
  # data<-Z.summary
  #
  par.ls<-Param.ls(x,m.dim,r.dim)
  B<-par.ls$B
  Lambda<-par.ls$Lambda
  Sig.zeta<-par.ls$Sig.zeta
  Sig.mu<-par.ls$Sig.mu
  #
  Sigma.ls <- Sigma.y(x,m.dim,r.dim)
  Sigma.eta.z<-Sigma.eta(B,Sig.zeta,m.dim)
  Sigma.eta.zero<-Sigma.eta.z#+Sig.mu
  Sigma.eta.one<-B%*%Sigma.eta.zero
  #
  term1<-kronecker(
    Lambda%*%Sigma.eta.one,diag(r.dim))
  term2<-kronecker(
    diag(r.dim),Lambda%*%t(Sigma.eta.one))
  term3<-term2%*%K.mn(r.dim,m.dim)
  dSS0.dL <- term1+term3
  #
  term4<-kronecker(
    Lambda%*%Sigma.eta.one,diag(r.dim))
  term5<-kronecker(
    diag(r.dim),Lambda%*%t(Sigma.eta.one))
  term6<-term5%*%K.mn(r.dim,m.dim)
  dSS1.dL <- term4+term6
  #
  dSS.dL<-dvecSy.dTheta((dSS0.dL),(dSS1.dL),r.dim)
  # Score
  #dL1<-dpL1.terms(x,data,m.dim,r.dim,N.dim,T.dim)
  return(t(dSS.dL))
}

# WRT Sigma_eps
dLy.dSeps <- function(x,r.dim){
  #x<-theta
  #data<-Z.summary
  #
  dSS0.dS <- D.n(r.dim)%>%t
  dSS1.dS <- matrix(0,dim(dSS0.dS)[1],dim(dSS0.dS)[2])
  #
  dSS.dS<-dvecSy.dTheta(t(dSS0.dS),t(dSS1.dS),r.dim)
  # Score
  #dL1<-dpL1.terms(x,data,m.dim,r.dim,N.dim,T.dim)
  return(t(dSS.dS))
}

# WRT Sigma_zeta
dLy.dSzeta <- function(x,m.dim,r.dim){
  # x<-theta
  # data<-Z.summary
  #
  par.ls<-Param.ls(x,m.dim,r.dim)
  B <- par.ls$B
  Lambda <- par.ls$Lambda
  #
  inv.IminBxB<-solve(diag(m.dim^2)-kronecker(B,B))
  #
  term1 <- D.n(m.dim)%>%t
  term2<-term1%*%inv.IminBxB
  term3<-kronecker(t(Lambda),t(Lambda))
  dSS0.dS <- term1%*%term3
  #
  term4 <- kronecker(diag(m.dim),B)
  # dSS1.dS <- term1%*%term4%*%inv.IminBxB%*%term3
  dSS1.dS <- matrix(0,dim(dSS0.dS)[1],dim(dSS0.dS)[2])
  #
  dSS.dS<-dvecSy.dTheta(t(dSS0.dS),t(dSS1.dS),r.dim)
  # Score
  #dL1<-dpL1.terms(x,data,m.dim,r.dim,N.dim,T.dim)
  return(t(dSS.dS))
}

# WRT Sigma_mu
dLy.dSmu <- function(x,m.dim,r.dim){
  # x<-theta
  # data<-Z.summary
  #
  par.ls<-Param.ls(x,m.dim,r.dim)
  B <- par.ls$B
  Lambda <- par.ls$Lambda
  #
  inv.IminBxB<-solve(diag(m.dim^2)-kronecker(B,B))
  #
  term1 <- D.n(m.dim)%>%t
  term2 <- kronecker(t(Lambda),t(Lambda))
  term3<- kronecker(diag(m.dim),B)
  dSS0.dS <- term1%*%term2
  #
  dSS1.dS <-  matrix(0,dim(dSS0.dS)[1],dim(dSS0.dS)[2])
  #dSS1.dS <- term1%*%term3%*%term2
  # dim(dSS.dS)
  dSS.dS<-dvecSy.dTheta(t(dSS0.dS),t(dSS1.dS),r.dim)
  # Score
  #dL1<-dpL1.terms(x,data,m.dim,r.dim,N.dim,T.dim)
  return(t(dSS.dS))
}

###### Score function
ScoreY.penFIMLit<-function(x,data,m.dim,r.dim,N.dim,T.dim,pen.l){
  # x<-theta0
  # data<-Z.summary
  # 
  # pen.l <- lambda.w
  mean.Y <- data[['mean']]
  Y.l0 <- data[['Y.l0']]
  Y.l1 <- data[['Y.l1']]
  Svcov.Z <- data[['VcoV']]
  # Compute Sigma.x 
  Sigma.ls <- Sigma.y(x, m.dim, r.dim)
  Sigma.x <- Vcov.matrix(Sigma.ls)
  # Sigma.x <- (Sigma.x + t(Sigma.x)) / 2
  # Sigma.x <- Sigma.x + diag(1e-6, nrow(Sigma.x))
  # Compute Eigen Decomposition
  eig_full <- eigen(Sigma.x, symmetric = TRUE)
  eig_vals <- pmax(eig_full$values, 1e-4)
  # Compute Sigma Inverse
  Sigma.x <- eig_full$vectors%*%diag(eig_vals)%*%t(eig_full$vectors)
  L.x <- chol(Sigma.x)
  Sigma.invx <- chol2inv(L.x)
  vec.Sigma.invx <- Sigma.invx %>% as.vector
  # Projection Matrix
  M <- Matrix(rep(1,2)%>%kronecker(diag(r.dim)),sparse = TRUE)
  P <- diag(r.dim*2) - M%*%solve(t(M)%*%M)%*%t(M)
  tY <- P%*%mean.Y
  # Observed data
  outer_full <- tcrossprod(tY)
  terms.1<-Sigma.invx%*%(Svcov.Z+outer_full)
  terms.2<-terms.1%*%Sigma.invx
  vec.terms<-terms.2%>%as.vector
  #
  Score.B <- dLy.dB(x,m.dim,r.dim)
  Score.L <- dLy.dLambda(x,m.dim,r.dim)
  Score.Szeta <- dLy.dSzeta(x,m.dim,r.dim)
  Score.Seps <- dLy.dSeps(x,r.dim)
  #Score.Smu <- dLy.dSmu(x,m.dim,r.dim)
  # Penalty
  # Sigma.invdiag <- diag(Sigma.invx)
  # Sigma.invpen <- Sigma.invx - diag(Sigma.invdiag)
  # Sigma.penvec <- Sigma.invpen%>%as.vector
  penalty.vec <- pen.l*sign(x) #Sigma.penvec 
  #
  dL1.out <- rep(0,length(x))
  for (n in 1:N.dim) {
    # n<-80
    Y.l0.n <- Y.l0[,n,]
    Y.l1.n <- Y.l1[,n,]
    Y.l0l1.n <- rbind(Y.l0.n, Y.l1.n)
    # dL1.out.t <- rep(0,length(x))
    for (t in 1:T.dim) {
      # t<-49 
      Y.l0l1.nt <- Y.l0l1.n[, t]
      obs_idx <- which(!is.na(Y.l0l1.nt))
      #obs_idx+r.dim
      if (length(obs_idx) == 0) next
      # Compute vec indices for submatrix
      row.s <- rep(obs_idx,times=length(obs_idx))
      col.s <- rep(obs_idx,each =length(obs_idx)) 
      vec_indices <- row.s+(col.s-1)*2*r.dim
      # Subsets
      vec.1 <- vec.terms[vec_indices]#as.vector(terms.2)
      vec.2 <- vec.Sigma.invx[vec_indices]
      dl.term <- -vec.1+N.dim*vec.2
      #pen.term <- sign(penalty.vec[vec_indices])
      #
      Score.vect<- rbind(
        Score.B[,vec_indices],
        Score.L[,vec_indices],
        Score.Szeta[,vec_indices],
        Score.Seps[,vec_indices])#,
        #Score.Smu[,vec_indices])
      #Score.vect<-Score.vect%>%as.numeric
      dL1.out <- dL1.out+(dl.term)%*%t(Score.vect)
    }
  }
  dfiml <- as.vector(dL1.out/(N.dim*T.dim))
  dreg.fiml <- dfiml+penalty.vec
  dreg.fiml
  # grad_num #*T.dim
  # pdml
  return(dreg.fiml)
}

##### Data,stats & summary
summary.stats<-function(y.var){
  # Lagged variables
  y.var.lag1 <- array(0, dim = dim(y.var))
  y.var.lag1[,,2:T.dim] <- y.var[,,1:(T.dim-1)]
  
  # Sample mean vector
  bar.yl1<-apply(y.var.lag1,c(1,2),'mean',na.rm = TRUE)
  bar.yl0<-apply(y.var,c(1,2),'mean',na.rm = TRUE)
  dim(bar.yl0)
  bar.Z<-rbind(bar.yl1,bar.yl0)
  
  # Sample covariance matrix
  S.l1<-S.l0<-list()
  S.l1.sum<-S.l0.sum<-matrix(0,p.dim,p.dim)
  for(i in 1:N.dim){
    S.l1[[i]]<-(T.dim/(T.dim-1))*var(t(y.var.lag1[,i,]),t(y.var[,i,]),na.rm=TRUE)
    S.l0[[i]]<-(T.dim/(T.dim-1))*var(t(y.var[,i,]),na.rm=TRUE)
  }
  S.l1.sum <- Reduce(`+`, lapply(S.l1, function(mat) {
    if (anyNA(mat)) matrix(0, nrow = nrow(mat), ncol = ncol(mat)) else mat
  }))
  S.l0.sum <- Reduce(`+`, lapply(S.l0, function(mat) {
    if (anyNA(mat)) matrix(0, nrow = nrow(mat), ncol = ncol(mat)) else mat
  }))
  S.list<-list(S.l0.sum,S.l1.sum)
  S.data <- Vcov.matrix(S.list)
  
  # Summary statistics
  obs.summary<-list(
    Y.l0=y.var,Y.l1=y.var.lag1,mean=bar.Z,S.l0=S.l0,S.l1=S.l1,VcoV=S.data)
  return(obs.summary)  
}

