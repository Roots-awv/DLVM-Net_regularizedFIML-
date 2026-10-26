##################################################################
library("qgraph")
library("dplyr")
library("Matrix")
library("reshape2")
library("ggplot2")

##### Working directory
setwd("~/Documents/Research_Asa/ml-DLVM-Net/ml-DLVM-Net_brute")

#### Loading methods
source('./Computing_Penalized_Pseudo_FIML_LargeBT_CovMatrix.R')

#### Loading the data 
source("./revision2_dataprep1.R")

# Data management
y1.var <- apply(y1, c(3,1,2), t)
p.dim<-dim(y1.var)[1]
N.dim<-dim(y1.var)[2]
T.dim<-dim(y1.var)[3]
k.dim <- 7
na.all<-is.na(y1.var)*1
na.dim<-mean(na.all)

# Summary statistics
Z.summary<-summary.stats(y1.var)

#### Estimates and results
# Adaptive/standard Lasso
FIML.est<-read.csv("FIMLest-SAM-data.csv")
est.length<-dim(FIML.est)[1] -0.5*(k.dim*(k.dim+1))
est.inv<-1/FIML.est
inf.vals <- ifelse(est.inv == Inf,0,1)
weights <- inf.vals*est.inv
weights[is.na(weights)]<-max(weights,na.rm=T)

# AL-1
lambda.w.AL1 <- abs(weights)/sqrt(N.dim*T.dim)
lambda.w.AL1 <- lambda.w.AL1%>%as.vector
print(round(lambda.w.AL1$x,2))
lambda.alasso1 <- lambda.w.AL1$x

# AL-2 
lambda.w.AL2 <- abs(weights)*1/sqrt(N.dim*T.dim)
lambda.w.AL2 <- lambda.w.AL2%>%as.vector
lambda.w.AL2<-lambda.w.AL2$x
lambda.w.AL2[(1+k.dim^2):dim(FIML.est)[1]]<-sqrt(N.dim)
print(round(lambda.w.AL2,2))
lambda.alasso4 <- lambda.w.AL2

# AL-3 
lambda.w.AL3 <- abs(weights)/sqrt((N.dim+T.dim))
lambda.w.AL3 <- lambda.w.AL3%>%as.vector
print(round(lambda.w.AL3$x,2))
lambda.alasso5 <- lambda.w.AL3$x
length(lambda.alasso1)

# AL-4 
lambda.w.AL4 <- abs(weights)*5
lambda.w.AL4 <- lambda.w.AL4%>%as.vector
print(round(lambda.w.AL4$x,2))
lambda.alasso6 <- lambda.w.AL4$x

# AL-5 
lambda.w.AL5 <- abs(weights)*7.5
lambda.w.AL5 <- lambda.w.AL5%>%as.vector
print(round(lambda.w.AL5$x,2))
lambda.alasso7 <- lambda.w.AL5$x

# Lasso small
lambda.lasso1 <- rep(1/sqrt(N.dim*T.dim*p.dim*na.dim),length(lambda.alasso1))
lambda.lasso4 <- rep(sqrt(p.dim*0.5),length(lambda.alasso1))
lambda.lasso5 <- rep(sqrt(N.dim/T.dim),length(lambda.alasso1))
lambda.lasso6 <- rep(1,length(lambda.alasso1))
lambda.lasso7 <- rep(sqrt(1/T.dim),length(lambda.alasso1))
lambda.lasso9 <- rep(sqrt(T.dim/N.dim),length(lambda.alasso1))

# Lasso large
lambda.lasso10 <- rep(sqrt(N.dim),length(lambda.alasso1))
lambda.lasso11 <- rep(sqrt(N.dim/2),length(lambda.alasso1))
lambda.lasso12 <- rep(7.5,length(lambda.alasso1))
lambda.lasso13 <- rep(sqrt(T.dim),length(lambda.alasso1))
lambda.lasso14 <- rep(sqrt(T.dim*2),length(lambda.alasso1))
lambda.lasso15 <- rep(sqrt(T.dim/2),length(lambda.alasso1))

# All penalties
lambda.all <- cbind(
  lambda.alasso1,lambda.alasso4,lambda.alasso5,lambda.alasso6,lambda.alasso7,
  lambda.lasso1,lambda.lasso4,lambda.lasso5,lambda.lasso6,lambda.lasso7,lambda.lasso9,
  lambda.lasso10,lambda.lasso11,lambda.lasso12,lambda.lasso13,lambda.lasso14,lambda.lasso15,0
             )
#l0<-rep(0,est.length)

# Penalty names 
M <- Matrix(rep(1,2)%>%kronecker(diag(p.dim)),sparse = TRUE)

penalty.methods<-c(
  "alasso1","alasso4","alasso5","alasso6","alasso7",
  "lasso1","lasso4","lasso5","lasso6","lasso7","lasso9","lasso10",
  "lasso11","lasso12","lasso13","lasso14","lasso15","Unreg"
  )

penalty.names<-c("AL-1","AL-2","AL-3","AL-4","AL-5",
                 "L-1","L-2","L-3","L-4","L-5","L-6","L-7",
                 "L-8","L-9","L-10","L-11","L-12","Unreg")

id.B<-sample(1:49,49)
id.Szeta<-sample(50:(49+17*7),49)
#id.Smu<-sample((length(lambda.alasso1)-49):length(lambda.alasso1),10)
id.theta<-c(id.B, id.Szeta)
id.theta<-sort(id.theta)

for (pen in 1:length(penalty.methods)) {
  # Fisher information matrix
  Est <- read.csv(paste0(penalty.methods[pen],"-FIMLest-SAM-data.csv"))
  estimates <- Est %>% as.vector
  Hess.theta <- read.csv(paste0(penalty.methods[pen],"-FIMLHessian-SAM-data.csv"))
  Hess.theta<-as.matrix(Hess.theta)
  Hess.theta<-(Hess.theta+t(Hess.theta))/2

  # # Model selection  
  # df.free <- sum(abs(estimates$x) > 1e-5)
  # Neg.logLik <-  N.dim*T.dim * log.penFIMLit.Y(estimates$x,Z.summary,k.dim,p.dim,N.dim,T.dim,pen.l=lambda.all[pen])
  # AIC <- 2*(Neg.logLik + df.free)/(N.dim*T.dim)
  # BIC <- (2*Neg.logLik + df.free*log(N.dim*T.dim))/(N.dim*T.dim)
  # EBIC <- BIC + (log(choose(length(estimates$x),df.free)))/(N.dim*T.dim)
  # model.select<-list(c(df.free,Neg.logLik/(N.dim*T.dim),AIC,BIC,EBIC))
  # names(model.select)<-penalty.methods[pen]
  # select <- append(select, model.select)
  # #write.csv(model.select,file = paste0(pen,"Table-SAM-data.csv"), row.names = FALSE)
  
  # Invert
  Grad.theta <-ScoreY.penFIMLit(estimates$x,Z.summary,k.dim,p.dim,N.dim,T.dim,pen.l=lambda.all[pen])
  Hcal.theta<-Grad.theta%*%t(Grad.theta) 
  inv.Hess <- 2*solve(Hess.theta)

  # eigen_Hess <- eigen(Hess.theta, symmetric = TRUE)
  # eigenHess_values <- pmax(eigen_Hess$values, 0.025)  # Push eigenvalues up if too small  
  # Hess.theta <- eigen_Hess$vectors%*%diag(eigenHess_values)%*%t(eigen_Hess$vectors)
  # inv.Hess <- eigen_Hess$vectors%*%diag(1/eigenHess_values)%*%t(eigen_Hess$vectors)
  # L.Hess<-chol(Hess.theta)
  # inv.Hess <- chol2inv(L.Hess)
  
  Sigma.items <- estimates$x %>% Sigma.y(k.dim,p.dim) %>% Vcov.matrix
  Mu_u.est<-Mu_u.hat(solve(Sigma.items),M,Z.summary$mean) 
  Mu_eta.est<-rep(1,2)%>%kronecker(Mu_u.est)
  Grad.mu <- t(M) %*% solve(Sigma.items) %*% (Z.summary$mean-Mu_eta.est)
  Hcal.mu<-Grad.mu%*%t(Grad.mu)
  Hess.mu <- t(M) %*% solve(Sigma.items) %*% M
  inv.Hess_mu <- 2*solve(Hess.mu)
  
  # Variance-covariance matrix
  As.Var.W <- inv.Hess%*%Hcal.theta%*%t(inv.Hess)/(N.dim*T.dim)#2*inv.Hess/(N.dim*T.dim)
  As.Var.H <- 4*Hcal.theta/(N.dim*T.dim)#2*inv.Hess/(N.dim*T.dim)
  # max(As.Var)
  # min(As.Var)
  # diag(As.Var)
  As.Var_mu.W <- inv.Hess_mu%*%Hcal.mu%*%t(inv.Hess_mu)/(N.dim*T.dim)#/(N.dim*T.dim)#2*inv.Hess/(N.dim*T.dim)
  As.Var_mu.H <- 4*Hcal.mu/(N.dim*T.dim)#/(N.dim*T.dim)#2*inv.Hess/(N.dim*T.dim)
  # dim(As.Var_mu)
  
  # # Heatmap B
  # theta_names <- parse(text = paste0("theta[", id.theta, "]"))
  # As.Var.theta <- As.Var[id.theta,id.theta]
  # As.Var.melt <- As.Var.theta %>% melt
  # head(As.Var.melt)
  # As.Var.gg <- ggplot(As.Var.melt, aes(x = Var1, y = Var2, fill = value)) +
  #           geom_tile() + theme_minimal() + 
  #           theme(
  #             axis.title = element_blank(),
  #             axis.text = element_blank(),
  #             axis.ticks = element_blank(),
  #             plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
  #             plot.subtitle = element_text(hjust = 0.5, size = 12)
  #           ) +
  #           #scale_x_discrete(labels = theta_names) +
  #           scale_y_discrete(labels = theta_names,limits =rev) +
  #           scale_fill_viridis_c(option = 'D',name = "Values") +
  #           ggtitle(expression(As.Var(theta)),
  #                   subtitle = paste0("(Case of ",penalty.names[pen],")"))
  # ggsave(paste0("As.Var",penalty.methods[pen],".jpeg"),As.Var.gg)          
  
  # Standard errors
  SE.W <- sqrt(diag(As.Var.W))
  l_CI95.W<-Est+qnorm(0.025)*SE.W
  u_CI95.W<-Est+qnorm(0.975)*SE.W

  SE.H <- sqrt(diag(As.Var.H))
  l_CI95.H<-Est+qnorm(0.025)*SE.H
  u_CI95.H<-Est+qnorm(0.975)*SE.H
  
  mean_mueta <- apply(Mu_eta.est[1:17,],1,mean)
  SE_mu.W <- sqrt(diag(As.Var_mu.W))
  l_CI95_mu.W<-mean_mueta+qnorm(0.025)*SE_mu.W
  u_CI95_mu.W<-mean_mueta+qnorm(0.975)*SE_mu.W
  
  SE_mu.H <- sqrt(diag(As.Var_mu.H))
  l_CI95_mu.H<-mean_mueta+qnorm(0.025)*SE_mu.H
  u_CI95_mu.H<-mean_mueta+qnorm(0.975)*SE_mu.H
  
  # Table
  lBFGSb.FIML.results<-cbind(Est,SE.H,l_CI95.H,u_CI95.H,SE.W,l_CI95.W,u_CI95.W)
  table<-round(lBFGSb.FIML.results,4)%>%as.data.frame
  head(table)
  colnames(table)<-c('Est','SE_H','Lci_H','Uci_H','SE_W','Lci_W','Uci_W')
  write.csv(round(table,4),file = paste0(penalty.methods[pen],"Table-SAM-data.csv"), row.names = FALSE)
  
  lBFGSb.FIML_mu<-cbind(mean_mueta,SE_mu.H,l_CI95_mu.H,u_CI95_mu.H,SE_mu.W,l_CI95_mu.W,u_CI95_mu.W)
  table_mu<-round(lBFGSb.FIML_mu,4)%>%as.data.frame
  head(table_mu)
  colnames(table_mu)<-c('Est','SE_H','Lci_H','Uci_H','SE_W','Lci_W','Uci_W')
  write.csv(round(table_mu,4),file = paste0(penalty.methods[pen],"Table_mu-SAM-data.csv"), row.names = FALSE)
  
}

###########################
# Model selection  

ini.ly <-vector("list", k.dim)
for (m in c(1,k.dim-1,k.dim)){ini.ly[[m]]<-rep(1,3)}
for (m in 2:(k.dim-2)){ini.ly[[m]]<-rep(1,2)}
ini.Lambda <- ini.ly %>% bdiag %>% as.vector

select<-list()
values.select<-c("df.H","df.W","Fit","EBIC.H","EBIC.W")

for (pen in 1:length(penalty.methods)) {
  table <- read.csv(paste0(penalty.methods[pen],"Table-SAM-data.csv"))
  df.free.W<-(table$Lci_W < 0 & table$Uci_W > 0)%>%ifelse(0,1)%>%sum
  df.free.H<-(table$Lci_H < 0 & table$Uci_H > 0)%>%ifelse(0,1)%>%sum
  
  #df.free <- sum(abs(estimates$x) > 1e-2)
  Neg.logLik <-  N.dim*T.dim * log.penFIMLit.Y(table$Est,Z.summary,k.dim,p.dim,N.dim,T.dim,pen.l=0)
  # AIC <- 2*(Neg.logLik + df.free)/(N.dim*T.dim)
  BIC.H <- (2*Neg.logLik + df.free.H*log(N.dim*T.dim))/(N.dim*T.dim)
  EBIC.H <- BIC.H + (log(choose(length(table$Est),df.free.H)))/(N.dim*T.dim)
  BIC.W <- (2*Neg.logLik + df.free.W*log(N.dim*T.dim))/(N.dim*T.dim)
  EBIC.W <- BIC.W + (log(choose(length(table$Est),df.free.W)))/(N.dim*T.dim)
  model.select<-list(c(df.free.H,df.free.W,Neg.logLik/(N.dim*T.dim),EBIC.H,EBIC.W))
  names(model.select)<-penalty.methods[pen]
  select <- append(select, model.select)
}

selection0 <- data.frame(select)#/(N.dim*T.dim)
#lambda.choices <- c(1/sqrt(N.dim*T.dim),lambda.all)
selection0<- rbind(lambda.all[1,],selection0)%>%data.frame
row.names(selection0) <- c("Penalty",values.select)
dfree.c<-selection0['df',-15]
selection<-rbind(selection0[,-15])
colnames(selection) <- penalty.names[-16]
min(selection['Fit',])
min(selection['EBIC.H',])
min(selection['EBIC.W',])
selection.mix<-selection%>%round(2)%>%data.frame
#selection_sorted<-selection.mix[order(selection.mix$Penalty,decreasing = F),]

# Second best
selection.2nd.best<-selection.mix[,-9]
min(selection.2nd.best['Fit',])
min(selection.2nd.best['EBIC.H',])
min(selection.2nd.best['EBIC.W',])

# # Example data
# df <- data.frame(A = c(3,1,4), B = c(1,5,9), C = c(2,6,5))
# rownames(df) <- c("row1", "row2", "row3")
# 
# # Order columns by the values in "row2"
# df_ordered <- df[, order(df["row2", ])]
# df_ordered

# Latex Table
selection.L<-selection.mix[,-c(1:5,17)]
selection.aLasso<-selection.mix[,1:5]
selection.Unreg<-selection.mix[,17]
selection.Lasso<-selection.L[,order(unlist(selection.L['Penalty',]))]
Lasso.names<-paste0('L.',as.character(1:11))
colnames(selection.Lasso)<-Lasso.names
selection.forall<-cbind(selection.aLasso,selection.Lasso,selection.Unreg)
colnames(selection.forall)<-penalty.names[-17]

model.selection<-selection.forall['EBIC.W',]
print(model.selection)

abs.diff.EBIC<-abs(model.selection-model.selection$`L-4`)
#ggplot(abs.diff.EBIC,a)

#stargazer::stargazer(t(model.selection), summary = F)


####### Extracting Network

# is_valid_covariance <- function(Sigma, tol = 1e-8) {
#   Sigma<-S
#   is_finite <- all(is.finite(Sigma))
#   is_symmetric <- isSymmetric(Sigma, tol = tol)
#   eigenvalues <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
#   is_pos_semidefinite <- all(eigenvalues >= -tol)
#   
#   return(is_finite && is_symmetric && is_pos_semidefinite)
# }

cov.to.pcorr <- function(S) {
  S <- (S + t(S)) / 2  # Force symmetry
  eigen_S <- eigen(S, symmetric = TRUE)
  eigen_values <- pmax(eigen_S$values, 1e-3)  # Push eigenvalues up if too small  
  S <- eigen_S$vectors%*%diag(eigen_values)%*%t(eigen_S$vectors)
  S.inv <- eigen_S$vectors%*%diag(1/eigen_values)%*%t(eigen_S$vectors)
  # Correlation matrix
  D <- sqrt(diag(S))
  R <- diag(1/D)%*%S%*%diag(1/D)
  R <- (R + t(R)) / 2  # Symmetrize again±
  diag(R) <- 1
  # Partial correlation matrix
  Delta <- sqrt(diag(S.inv))
  pR <- diag(1/Delta)%*%S.inv%*%diag(1/Delta)
  pR <- (pR + t(pR)) / 2  # Symmetrize again
  diag(pR) <- 0
  return(-pR)
}

############# CI based
#labels<-c("content1","content2","anxiety","understand","stress","PA","NA")
labels<-c("SIC","CC","Anx","SAU","Str","PA","NA")

###############
Net<- 9
#Net1
Net_1.ini<-read.csv(paste0(penalty.methods[Net],"Table-SAM-data.csv"))

head(Net_1.ini)
Net_1<-Net_1.ini
Net_1$Est[Net_1$Lci_H < 0 & Net_1$Uci_H > 0] <- 0
cbind(Net_1.ini,Net_1$Est)


# Temporal network 
par.matrices_net1<-Param.ls(Net_1$Est,k.dim,p.dim)
tempo.net_net1<-par.matrices_net1$B%>%as.matrix
colnames(tempo.net_net1)<-labels
rownames(tempo.net_net1)<- labels
max(abs(tempo.net_net1))
min(abs(tempo.net_net1))
sum(abs(tempo.net_net1)!=0)

# Contemporaneous network
par.matrices_net1.ini<-Param.ls(Net_1.ini$Est,k.dim,p.dim)
par.matrices_net1.to0<-Param.ls(Net_1$Est,k.dim,p.dim)
Szeta.est_net1.ini<-par.matrices_net1.ini$Sig.zeta
Szeta.est_net1_to0<-(par.matrices_net1.to0$Sig.zeta==Szeta.est_net1.ini)*1
contempo.cor_net1<-Szeta.est_net1.ini%>%cov.to.pcorr
contempo.cor_net1_to0<-contempo.cor_net1*Szeta.est_net1_to0
max(abs(contempo.cor_net1))
sum(abs(contempo.cor_net1[lower.tri(contempo.cor_net1,diag = F)])>1e-5)

# Between network
Sigma_u.est_net1<-Sigma_u.hat(Net_1.ini$Est,Z.summary$mean,k.dim,p.dim,0)
between.cor_net1<-Sigma_u.est_net1%>%cov.to.pcorr
max(between.cor_net1)
sum(abs(between.cor_net1[lower.tri(between.cor_net1,diag = F)])>1e-5)

# Net2
Net_2.ini<-read.csv(paste0(penalty.methods[Net],"Table-SAM-data.csv"))

head(Net_2.ini)
# Net_2<-Net_2.ini
# Net_2$Est[Net_2$CI_lower < 0 & Net_2$CI_upper > 0] <- 0
# cbind(Net_2.ini,Net_2$Est)
Net_2<-Net_2.ini
Net_2$Est[Net_2$Lci_W < 0 & Net_2$Uci_W > 0] <- 0
#cbind(Net.ini,Net_2$Est)

# Temporal network 
par.matrices_net2<-Param.ls(Net_2$Est,k.dim,p.dim)
tempo.net_net2<-par.matrices_net2$B%>%as.matrix
colnames(tempo.net_net2)<-labels
rownames(tempo.net_net2)<- labels

# # Contemporaneous network
# par.matrices_net2.ini<-Param.ls(Net_2$Est,k.dim,p.dim)
# Szeta.est_net2_to0<-(par.matrices_net2$Sig.zeta<1e-5)*1
# Szeta.est_net2.ini<-par.matrices_net2.ini$Sig.zeta
# contempo.cor_net2<-Szeta.est_net2.ini%>%cov.to.pcorr
# contempo.cor_net2_to0<-Szeta.est_net2_to0*contempo.cor_net2
# max(abs(contempo.cor_net2))
# sum(abs(contempo.cor_net2[lower.tri(contempo.cor_net2,diag = F)])>0.1)

# Contemporaneous network
par.matrices_net2.ini<-Param.ls(Net_2.ini$Est,k.dim,p.dim)
par.matrices_net2.to0<-Param.ls(Net_2$Est,k.dim,p.dim)
Szeta.est_net2.ini<-par.matrices_net2.ini$Sig.zeta
Szeta.est_net2_to0<-(par.matrices_net2.to0$Sig.zeta==Szeta.est_net2.ini)*1
contempo.cor_net2<-Szeta.est_net2.ini%>%cov.to.pcorr
contempo.cor_net2_to0<-contempo.cor_net2*Szeta.est_net2_to0
max(abs(contempo.cor_net2))
sum(abs(contempo.cor_net2[lower.tri(contempo.cor_net2,diag = F)])>0.1)

# Between network
Sigma_u.est_net2<-Sigma_u.hat(Net_2.ini$Est,Z.summary$mean,k.dim,p.dim,0)
between.cor_net2<-Sigma_u.est_net2%>%cov.to.pcorr
max(between.cor_net2)
sum(abs(between.cor_net2)>0.1)
sum(abs(between.cor_net2[lower.tri(between.cor_net2,diag = F)])>0.1)

#### threshold

selected.est2 <- read.csv(paste0(penalty.methods[Net],"-FIMLest-SAM-data.csv"))
selected.estimates2 <- selected.est2 %>% as.vector
selected.table2 <-read.csv(paste0(penalty.methods[Net],"Table-SAM-data.csv"))

par.matrices2<-Param.ls(selected.estimates2$x,k.dim,p.dim)
tempo.net3<-par.matrices2$B%>%as.matrix
colnames(tempo.net3)<-labels
rownames(tempo.net3)<- labels
max(abs(tempo.net3))
min(abs(tempo.net3))
sum(abs(tempo.net3)>0.01)

Szeta.est2<-par.matrices2$Sig.zeta
contempo.cor3<-Szeta.est2%>%cov.to.pcorr
max(abs(contempo.cor3))
sum(abs(contempo.cor3)>0.01)
sum(abs(contempo.cor3[lower.tri(contempo.cor3,diag = F)])>0.01)

Sigma_u.est2<-Sigma_u.hat(selected.estimates2$x,Z.summary$mean,k.dim,p.dim,0)
between.cor3<-Sigma_u.est2%>%cov.to.pcorr
max(between.cor3)
sum(abs(between.cor3)>0.1)
sum(abs(between.cor3[lower.tri(between.cor3,diag = F)])>0.01)

######### Network analysis
### nb og edges

# Temporal net
tempo.H<-(tempo.net_net1!=0)*1
tempo.W<-(tempo.net_net2!=0)*1
tempo.01<-(tempo.net3!=0)*1
nb.tempo.in<-rbind(colSums(tempo.01),colSums(tempo.H),colSums(tempo.W))
row.names(nb.tempo.in)<-c('Threshold','OPG','Sandwich')
nb.tempo.in%>%colMeans%>%round%>%print

# Result: Anxeity and PA # the most central latent construct 
# when we only take into account nodes that received the most incoming arrows

nb.tempo.out<-rbind(rowSums(tempo.01),rowSums(tempo.H),rowSums(tempo.W))
row.names(nb.tempo.out)<-c('Threshold','OPG','Sandwich')
nb.tempo.out%>%colMeans%>%round%>%print
nb.tempo.out
abs(nb.tempo.out-nb.tempo.in)
# Result: Anxeity and PA # the most central latent construct 
# when we only take into account nodes that received the most incoming arrows

# Contempo net
contempo.H<-(abs(contempo.cor_net1_to0)>1e-5)*1
contempo.W<-(abs(contempo.cor_net2_to0)>1e-5)*1
contempo.01<-(abs(contempo.cor3)>1e-2)*1
nb.contempo<-rbind(colSums(contempo.01),colSums(contempo.H),colSums(contempo.W))
row.names(nb.contempo)<-c('Threshold','OPG','Sandwich')
colnames(nb.contempo)<-colnames(nb.tempo.in)
nb.contempo%>%colMeans%>%round%>%print
contempo.intersect<-contempo.01*contempo.H*contempo.W
colnames(contempo.intersect)<-colnames(nb.tempo.in)
rownames(contempo.intersect)<-colnames(nb.tempo.in)

# Result: By aggregating the result, it seemed that Content1 and NA was 
# the most central latent construct in the network

# Between net
between.H<-(abs(between.cor_net1)>0)*1
between.W<-(abs(between.cor_net2)>0)*1
between.01<-(abs(between.cor3)>1e-2)*1
nb.between<-rbind(colSums(between.01),colSums(between.H),colSums(between.W))
row.names(nb.between)<-c('Threshold','OPG','Sandwich')
nb.between%>%colMeans%>%round%>%print

###### Longest path
longest_directed_path <- function(A) {
  # Convert adjacency matrix to directed graph
  g <- igraph::graph_from_adjacency_matrix(A,mode = "directed",weighted = NULL,diag = FALSE)
  # Find all simple directed paths
  paths <- list()
  for (v in V(g)) {
    paths <- c(paths,all_simple_paths(g, from = v, mode = "out"))
  }
  # Length = number of edges
  lengths <- sapply(paths, function(x) length(x) - 1)
  # Longest path
  k <- which.max(lengths)
  list(path = paths[[k]],length = lengths[k])

  }

adj.net1<-(tempo.net_net1!=0)*1
adj.net2<-(tempo.net_net2!=0)*1
adj.net3<-(abs(tempo.net3)>1e-2)*1
dir_path.net1 <- longest_directed_path(adj.net1)
dir_path.net2 <- longest_directed_path(adj.net2)
dir_path.net3 <- longest_directed_path(adj.net3)
dir_path.net1$path
dir_path.net2$path
dir_path.net3$path

adj.net.inter<-adj.net1*adj.net2*adj.net3
dir_path.inter <- longest_directed_path(adj.net.inter)

# Visual 3
# width = 1054, height= 822

layout((matrix(1:9,3,3)))

qgraph(tempo.net3, layout = "circle", labels = labels,
       title = paste0("Temporal (Threshold)"), shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(paste("      ",B)),cex.main = 2,line=1.5)

qgraph(tempo.net_net2, layout = "circle", labels = labels,
       title = paste0("Temporal (Sandwich)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(paste("      ",B)),cex.main = 2,line=1.5)

qgraph(tempo.net_net1, layout = "circle", labels = labels,
       title = paste0("Temporal (OPG)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(paste("      ",B)),cex.main = 2,line=1.5)

qgraph(contempo.cor3, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (Threshold)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       theme = "TeamFortress",
       curve = 0,
       height=5,width=10)
box("figure")
title(expression(Omega^(zeta)), cex.main = 2,line=1.5)

qgraph(contempo.cor_net2_to0, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (Sandwich)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(zeta)), cex.main = 2,line=1.5)

qgraph(contempo.cor_net1_to0, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (OPG)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(zeta)), cex.main = 2,line=1.5)

### Between network
#layout(matrix(1:4,2,2))

qgraph(between.cor3, layout = "circle", labels = labels,
       title = paste0("Between (Threshold)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(mu)),cex.main = 2,line=1.5)

qgraph(between.cor_net2, layout = "circle", labels = labels,
       title = paste0("Between (Sandwich)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(mu)),cex.main = 2,line=1.5)

qgraph(between.cor_net1, layout = "circle", labels = labels,
       title = paste0("Between (OPG)"),shape = "circle",
       vsize = 15, vsize2 = 15,
       mar = rep(8,4), asize = 5,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       label.cex = 1.5,
       fade = FALSE,
       curve = 0,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(mu)),cex.main = 2,line=1.5)

# Witdth = 1125 ; Height = 1012.5


#############

mu.param <- 1:p.dim
Mod1.mueta <- read.csv(paste0(penalty.methods[Net1],"Table_mu-SAM-data.csv"))
Mod1.mueta.df<-data.frame(cbind(mu.param,Mod1.mueta),Method=penalty.names[Net1])
c(min(Mod1.mueta.df$Lci_W),max(Mod1.mueta.df$Uci_W))

Mod2.mueta <- read.csv(paste0(penalty.methods[Net2],"Table_mu-SAM-data.csv"))
Mod2.mueta.df<-data.frame(cbind(mu.param,Mod2.mueta),Method=penalty.names[Net2])
c(min(Mod2.mueta.df$Lci_W),max(Mod2.mueta.df$Uci_W))

mueta.df <- rbind(Mod1.mueta.df,Mod2.mueta.df)

mueta.gg <- ggplot(mueta.df, aes(x = mu.param, y = Est,color = Method,shape = Method,group = Method)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Lci_W, ymax = Uci_W), 
                width = 1,position = position_dodge(width = 0.5)) +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(0.5, 3.5))+
  scale_color_manual(values = c("#B36B00","#00724F")) +
  labs(
    x = "",
    y = expression(mu^(eta) ~ "estimates"),
    title = "(a)"
  ) + theme(
    legend.position = "none",
    axis.text.y = element_blank(),     # remove tick labels
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank(),     # remove axis title
    plot.title = element_text(hjust = 0.5)
  )

mueta.gg

#############
edge.labels<-matrix('',length(labels),length(labels))
for (node1 in 1:length(labels)) {
  for (node2 in 1:length(labels)) {
    edge.labels[node1,node2]<-paste0(labels[node1],'--',labels[node2])
  }
}
edge.rank <- edge.labels[lower.tri(edge.labels,diag=T)]

# Models
Mod1<-read.csv(paste0(penalty.methods[Net],"Table-SAM-data.csv"))
Mod2<-read.csv(paste0(penalty.methods[Net],"Table-SAM-data.csv"))

####### B
H.name<-c('Est','SE_H','Lci_H','Uci_H')
W.name<-c('Est','SE_W','Lci_W','Uci_W')
real.name<-c('Est','SE','Lci','Uci')

B.param<-matrix(1:(k.dim^2),k.dim,k.dim)
B.param<-B.param[lower.tri(B.param,diag=T)]

Mod1.B.df<-data.frame(cbind(edge.rank,Mod1[B.param,H.name]),Asymptotics="OOG")
colnames(Mod1.B.df)<-c('Edges',real.name,'Std. Error')
min(Mod1.B.df$Lci)
max(Mod1.B.df$Uci)

Mod2.B.df<-data.frame(cbind(edge.rank,Mod2[B.param,W.name]),Asymptotics='Sandwich')
colnames(Mod2.B.df)<-c('Edges',real.name,'Std. Error')
min(Mod2.B.df$Lci)
max(Mod2.B.df$Uci)

B.df<-rbind(Mod2.B.df,Mod1.B.df)
col.B<-names(B.df)
#names(B.df)<-c('Edges',col.B[-1])
B.df$Edges<-factor(B.df$Edges,levels = unique(B.df$Edges))

B.gg <- ggplot(B.df,  aes(x = Edges, y = Est,color =`Std. Error`,shape =`Std. Error`,group =`Std. Error`)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = Lci, ymax = Uci), 
                width = 1,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(-0.4, 0.65))+
  scale_color_manual(values = c("darkgray","black")) +  
  labs(
    x = "",
    y = expression(B ~ "estimates"),
    title = "(a)"
  ) +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank(),     # remove axis title
    plot.title = element_text(hjust = 0.5)
  ) 
B.gg

### Szeta
L.end <- k.dim^2+k.dim*p.dim
Szeta.param <- (L.end+1):(L.end + k.dim * (k.dim + 1) * 0.5)
Mod1.Szeta.df<-data.frame(cbind(edge.rank,Mod1[Szeta.param,H.name]),Asymptotics="OOG")
colnames(Mod1.Szeta.df)<-c('Edges',real.name,'Std. Error')
min(Mod1.Szeta.df$Lci)
max(Mod1.Szeta.df$Uci)

Mod2.Szeta.df<-data.frame(cbind(edge.rank,Mod2[Szeta.param,W.name]),Asymptotics='Sandwich')
colnames(Mod2.Szeta.df)<-c('Edges',real.name,'Std. Error')
min(Mod2.Szeta.df$Lci)
max(Mod2.Szeta.df$Uci)

Szeta.df<-rbind(Mod2.Szeta.df,Mod1.Szeta.df)
col.Szeta<-names(Szeta.df)
#names(Szeta.df)<-c('Edges',col.Szeta[-1])
Szeta.df$Edges<-factor(Szeta.df$Edges,levels = unique(Szeta.df$Edges))

Szeta.gg <- ggplot(Szeta.df, aes(x = Edges,  y = Est,color =`Std. Error`,shape=`Std. Error`,group=`Std. Error`)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = Lci, ymax = Uci), 
                width = 1,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(-2, 2))+
  scale_color_manual(values = c("darkgray","black")) +  
  
  labs(
    x = "",
    y = expression(Sigma^{zeta} ~ "estimates"),
    title = "(b)"
  ) + theme(
    axis.text.y = element_blank(),     # remove tick labels
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank(),     # remove axis title
    plot.title = element_text(hjust = 0.5)
  )

Szeta.gg

# ggtitle("Parameter Estimates with 95% Confidence Intervals")
# Mod1.gg2 <-Mod1.Lambda.gg + Mod1.Seps.gg + plot_layout(ncol=2)
library(patchwork)

# Mod1.gg1 <- B.gg + Szeta.gg + plot_layout(ncol=2) + 
#   plot_annotation(title = "AL-1 estimates with 95% Confidence Intervals",
#                   theme = theme(plot.title = element_text(hjust = 0.5)))
# Mod1.gg1 

CI.gg1 <-  B.gg + Szeta.gg + plot_layout(ncol=2) + 
  plot_annotation(title = "L-4 estimates with 95% Confidence Intervals",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
CI.gg1 

####################################################################################
#############

mu.param <- 1:p.dim
Mod1.mueta <- read.csv(paste0(penalty.methods[Net1],"Table_mu-SAM-data.csv"))
Mod1.mueta.df<-data.frame(cbind(mu.param,Mod1.mueta),Method=penalty.names[Net1])
c(min(Mod1.mueta.df$Lci_W),max(Mod1.mueta.df$Uci_W))

Mod2.mueta <- read.csv(paste0(penalty.methods[Net2],"Table_mu-SAM-data.csv"))
Mod2.mueta.df<-data.frame(cbind(mu.param,Mod2.mueta),Method=penalty.names[Net2])
c(min(Mod2.mueta.df$Lci_W),max(Mod2.mueta.df$Uci_W))

mueta.df <- rbind(Mod1.mueta.df,Mod2.mueta.df)

mueta.gg <- ggplot(mueta.df, aes(x = mu.param, y = Est,color = Method,shape = Method,group = Method)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Lci_W, ymax = Uci_W), 
                width = 1,position = position_dodge(width = 0.5)) +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(0.5, 3.5))+
  scale_color_manual(values = c("#B36B00","#00724F")) +
  labs(
    x = "",
    y = expression(mu^(eta) ~ "estimates"),
    title = "(a)"
  ) + theme(
    legend.position = "none",
    axis.text.y = element_blank(),     # remove tick labels
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank(),     # remove axis title
    plot.title = element_text(hjust = 0.5)
  )

mueta.gg

#############
edge.labels<-matrix('',length(labels),length(labels))
for (node1 in 1:length(labels)) {
  for (node2 in 1:length(labels)) {
    edge.labels[node1,node2]<-paste0(labels[node1],'--',labels[node2])
  }
}
edge.rank <- edge.labels[lower.tri(edge.labels,diag=T)]

# Models
Mod1<-read.csv(paste0(penalty.methods[Net1],"Table-SAM-data.csv"))
Mod2<-read.csv(paste0(penalty.methods[Net2],"Table-SAM-data.csv"))

####### B
B.param<-matrix(1:(k.dim^2),k.dim,k.dim)
B.param<-B.param[lower.tri(B.param,diag=T)]

Mod1.B.df<-data.frame(cbind(edge.rank,Mod1[B.param,]),Method=penalty.names[Net1])
min(Mod1.B.df$Lci_W)
max(Mod1.B.df$Uci_W)

Mod2.B.df<-data.frame(cbind(edge.rank,Mod2[B.param,]),Method=penalty.names[Net2])
min(Mod2.B.df$Lci_W)
max(Mod2.B.df$Uci_W)

B.df<-rbind(Mod1.B.df,Mod2.B.df)
col.B<-names(B.df)
names(B.df)<-c('Edges',col.B[-1])
B.df$Edges<-factor(B.df$Edges,levels = unique(B.df$Edges))

B.gg <- ggplot(B.df,  aes(x = Edges, y = Est,color = Method,shape = Method,group = Method)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = Lci_W, ymax = Uci_W), 
                width = 1,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(-0.6, 0.7))+
  scale_color_manual(values = c("#B36B00","#00724F")) +  
  labs(
    x = "",
    y = expression(B ~ "estimates"),
    title = "(a)"
  ) +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank(),     # remove axis title
    plot.title = element_text(hjust = 0.5)
  ) 
B.gg

### Szeta
L.end <- k.dim^2+k.dim*p.dim
Szeta.param <- (L.end+1):(L.end + k.dim * (k.dim + 1) * 0.5)
Mod1.Szeta.df<-data.frame(cbind(edge.rank,Mod1[Szeta.param,]),Method=penalty.names[Net1])
min(Mod1.Szeta.df$Lci_W)
max(Mod1.Szeta.df$Uci_W)

Mod2.Szeta.df<-data.frame(cbind(edge.rank,Mod2[Szeta.param,]),Method=penalty.names[Net2])
min(Mod2.Szeta.df$Lci_W)
max(Mod2.Szeta.df$Uci_W)

Szeta.df<-rbind(Mod1.Szeta.df,Mod2.Szeta.df)
col.Szeta<-names(Szeta.df)
names(Szeta.df)<-c('Edges',col.Szeta[-1])
Szeta.df$Edges<-factor(Szeta.df$Edges,levels = unique(Szeta.df$Edges))

Szeta.gg <- ggplot(Szeta.df, aes(x = Edges, y = Est,color = Method,shape = Method,group = Method)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = Lci_W, ymax = Uci_W), 
                width = 1,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(-2, 2))+
  scale_color_manual(values = c("#B36B00","#00724F")) + 
  labs(
    x = "",
    y = expression(Sigma^{zeta} ~ "estimates"),
    title = "(b)"
  ) + theme(
    axis.text.y = element_blank(),     # remove tick labels
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank(),     # remove axis title
    plot.title = element_text(hjust = 0.5)
  )

Szeta.gg

# ggtitle("Parameter Estimates with 95% Confidence Intervals")
# Mod1.gg2 <-Mod1.Lambda.gg + Mod1.Seps.gg + plot_layout(ncol=2)
library(patchwork)

# Mod1.gg1 <- B.gg + Szeta.gg + plot_layout(ncol=2) + 
#   plot_annotation(title = "AL-1 estimates with 95% Confidence Intervals",
#                   theme = theme(plot.title = element_text(hjust = 0.5)))
# Mod1.gg1 

CI.gg1 <-  B.gg + Szeta.gg + plot_layout(ncol=2) + 
  plot_annotation(title = "Estimates with 95% Confidence Intervals",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
CI.gg1 


###### Threshold rule
# Selected networks
Net1<- 8
Net2<- 9
#labels<-c("content1","content2","anxiety","understand","stress","PA","NA")

selected.est1 <- read.csv(paste0(penalty.methods[Net1],"-FIMLest-SAM-data.csv"))
selected.estimates1 <- selected.est1 %>% as.vector
selected.table1<-read.csv(paste0(penalty.methods[Net1],"Table-SAM-data.csv"))
selected.table1$CI_upper-selected.table1$CI_lower

selected.est2 <- read.csv(paste0(penalty.methods[Net2],"-FIMLest-SAM-data.csv"))
selected.estimates2 <- selected.est2 %>% as.vector
selected.table2 <-read.csv(paste0(penalty.methods[Net2],"Table-SAM-data.csv"))

# Temporal network 
par.matrices1<-Param.ls(selected.estimates1$x,k.dim,p.dim)
tempo.net1<-par.matrices1$B%>%as.matrix
colnames(tempo.net1)<-labels
rownames(tempo.net1)<- labels
max(abs(tempo.net1))
min(abs(tempo.net1))
sum(abs(tempo.net1)>0.01)

par.matrices2<-Param.ls(selected.estimates2$x,k.dim,p.dim)
tempo.net2<-par.matrices2$B%>%as.matrix
colnames(tempo.net2)<-labels
rownames(tempo.net2)<- labels
max(abs(tempo.net2))
min(abs(tempo.net2))
sum(abs(tempo.net2)>0.01)

# Contemporaneous network
Szeta.est1<-par.matrices1$Sig.zeta
contempo.cor1<-Szeta.est1%>%cov.to.pcorr
max(abs(contempo.cor1))
min(abs(contempo.cor1))
sum(abs(contempo.cor1[lower.tri(contempo.cor1,diag = F)])>0.01)

Szeta.est2<-par.matrices2$Sig.zeta
contempo.cor2<-Szeta.est2%>%cov.to.pcorr
max(abs(contempo.cor2))
sum(abs(contempo.cor2)>0.01)
sum(abs(contempo.cor2[lower.tri(contempo.cor2,diag = F)])>0.01)

# Contemporaneous network
Sigma_u.est1<-Sigma_u.hat(selected.estimates1$x,Z.summary$mean,k.dim,p.dim,0)
between.cor1<-Sigma_u.est1%>%cov.to.pcorr
max(between.cor1)
sum(abs(between.cor1)>0.01)
sum(abs(between.cor1[lower.tri(between.cor1,diag = F)])>0.01)

Sigma_u.est2<-Sigma_u.hat(selected.estimates2$x,Z.summary$mean,k.dim,p.dim,0)
between.cor2<-Sigma_u.est2%>%cov.to.pcorr
max(between.cor2)
sum(abs(between.cor2)>0.1)
sum(abs(between.cor2[lower.tri(between.cor2,diag = F)])>0.01)

colfac<-max(max(between.cor1),max(contempo.cor1),max(tempo.net1),
            max(between.cor2),max(contempo.cor2),max(tempo.net2))

# # Network viszualization
# # Visual 1
# labels<-c("content1","content2","anxiety","understand","stress","PA","NA")
# 
# qgraph(tempo.net, layout = "circle", labels = labels,
#        title = "Temporal", shape = "rectangle", 
#        vsize = 15, vsize2 = 15,
#        mar = rep(8,4), asize = 7,
#        fade=T,
#        theme = "colorblind",
#        filetype='jpg',
#        height=5,width=10,
#        filename="TemporalNet-FIML"
# )
# 
# qgraph(contempo.cor, layout = "circle", labels = labels,
#        title = "Contemporaneous",shape = "rectangle", 
#        vsize = 15, vsize2 = 15,
#        mar = rep(8,4),
#        fade=T,
#        theme = "colorblind",
#        filetype='jpg',
#        height=5,width=10,
#        filename="ContemporanousNet-FIML"
# )
# 
# qgraph(between.cor, layout = "circle", labels = labels,
#        title = "Between",shape = "rectangle", 
#        vsize = 15, vsize2 = 15,
#        mar = rep(8,4),
#        fade=T,
#        theme = "colorblind",
#        filetype='jpg',
#        height=5,width=10,
#        filename="BetweenNet-FIML"
# )

# Visual 2
# width = 877, height= 412
layout((matrix(1:6,2,3)))
qgraph(tempo.net1, layout = "circle", labels = labels,
       title = paste0("Temporal (",penalty.names[Net1],")"), shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 7,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(paste("      ",B)),cex.main = 2,line=1.5)

qgraph(tempo.net2, layout = "circle", labels = labels,
       title = paste0("Temporal (",penalty.names[Net2],")"), shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 7,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(paste("      ",B)),cex.main = 2,line=1.5)

qgraph(contempo.cor1, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (",penalty.names[Net1],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(zeta)), cex.main = 2,line=1.5)

qgraph(contempo.cor2, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (",penalty.names[Net2],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(zeta)), cex.main = 2,line=1.5)

qgraph(between.cor1, layout = "circle", labels = labels,
       title = paste0("Between (",penalty.names[Net1],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(mu)),cex.main = 2,line=1.5)

qgraph(between.cor2, layout = "circle", labels = labels,
       title = paste0("Between (",penalty.names[Net2],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-2,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(mu)),cex.main = 2,line=1.5)
# Witdth = 1125 ; Height = 675


############# CI based

###############
#Net1
Net1<- 9

Net_1.ini<-read.csv(paste0(penalty.methods[Net1],"Table-SAM-data.csv"))

head(Net.ini)
Net_1<-Net.ini
Net_1$Est[Net_1$Lci_H < 0 & Net_1$Uci_H > 0] <- 0
cbind(Net.ini,Net_1$Est)


# Temporal network 
par.matrices_net1<-Param.ls(Net_1$Est,k.dim,p.dim)
tempo.net_net1<-par.matrices_net1$B%>%as.matrix
colnames(tempo.net_net1)<-labels
rownames(tempo.net_net1)<- labels
max(abs(tempo.net_net1))
min(abs(tempo.net_net1))
sum(abs(tempo.net_net1)!=0)

# Contemporaneous network
par.matrices_net1.ini<-Param.ls(Net_1.ini$Est,k.dim,p.dim)
par.matrices_net1.to0<-Param.ls(Net_1$Est,k.dim,p.dim)
Szeta.est_net1.ini<-par.matrices_net1.ini$Sig.zeta
Szeta.est_net1_to0<-(par.matrices_net1.to0$Sig.zeta==Szeta.est_net1.ini)*1
contempo.cor_net1<-Szeta.est_net1.ini%>%cov.to.pcorr
contempo.cor_net1_to0<-contempo.cor_net1*Szeta.est_net1_to0
max(abs(contempo.cor_net1))
sum(abs(contempo.cor_net1[lower.tri(contempo.cor_net1,diag = F)])>1e-5)

# Between network
Sigma_u.est_net1<-Sigma_u.hat(Net_1.ini$Est,Z.summary$mean,k.dim,p.dim,0)
between.cor_net1<-Sigma_u.est_net1%>%cov.to.pcorr
max(between.cor_net1)
sum(abs(between.cor_net1[lower.tri(between.cor_net1,diag = F)])>1e-5)

# Net2
Net2<-8
Net_2.ini<-read.csv(paste0(penalty.methods[Net2],"Table-SAM-data.csv"))

head(Net_2.ini)
Net_2<-Net_2.ini
Net_2$Est[Net_2$Uci_H < 0 & Net_2$Uci_H > 0] <- 0
cbind(Net_2.ini,Net_2$Est)

# Temporal network 
par.matrices_net2<-Param.ls(Net_2$Est,k.dim,p.dim)
tempo.net_net2<-par.matrices_net2$B%>%as.matrix
colnames(tempo.net_net2)<-labels
rownames(tempo.net_net2)<- labels

# # Contemporaneous network
# par.matrices_net2.ini<-Param.ls(Net_2$Est,k.dim,p.dim)
# Szeta.est_net2_to0<-(par.matrices_net2$Sig.zeta<1e-5)*1
# Szeta.est_net2.ini<-par.matrices_net2.ini$Sig.zeta
# contempo.cor_net2<-Szeta.est_net2.ini%>%cov.to.pcorr
# contempo.cor_net2_to0<-Szeta.est_net2_to0*contempo.cor_net2
# max(abs(contempo.cor_net2))
# sum(abs(contempo.cor_net2[lower.tri(contempo.cor_net2,diag = F)])>0.1)

# Contemporaneous network
par.matrices_net2.ini<-Param.ls(Net_2.ini$Est,k.dim,p.dim)
par.matrices_net2.to0<-Param.ls(Net_2$Est,k.dim,p.dim)
Szeta.est_net2.ini<-par.matrices_net2.ini$Sig.zeta
Szeta.est_net2_to0<-(par.matrices_net2.to0$Sig.zeta==Szeta.est_net2.ini)*1
contempo.cor_net2<-Szeta.est_net2.ini%>%cov.to.pcorr
contempo.cor_net2_to0<-contempo.cor_net2*Szeta.est_net2_to0
max(abs(contempo.cor_net2))
sum(abs(contempo.cor_net2[lower.tri(contempo.cor_net2,diag = F)])>0.1)

# Between network
Sigma_u.est_net2<-Sigma_u.hat(Net_2.ini$Est,Z.summary$mean,k.dim,p.dim,0)
between.cor_net2<-Sigma_u.est_net2%>%cov.to.pcorr
max(between.cor_net2)
sum(abs(between.cor_net2)>0.1)
sum(abs(between.cor_net2[lower.tri(between.cor_net2,diag = F)])>0.1)


# Visual 3
# width = 1054, height= 822
layout((matrix(1:6,2,3)))
qgraph(tempo.net_net1, layout = "circle", labels = labels,
       title = paste0("Temporal (",penalty.names[Net1],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 7,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(paste("      ",B)),cex.main = 2,line=1.5)

qgraph(tempo.net_net2, layout = "circle", labels = labels,
       title = paste0("Temporal (",penalty.names[Net2],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 7,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(paste("      ",B)),cex.main = 2,line=1.5)

qgraph(contempo.cor_net1_to0, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (",penalty.names[Net1],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(zeta)), cex.main = 2,line=1.5)

qgraph(contempo.cor_net2_to0, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (",penalty.names[Net2],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(zeta)), cex.main = 2,line=1.5)

qgraph(between.cor_net1, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (",penalty.names[Net1],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(mu)),cex.main = 2,line=1.5)

qgraph(between.cor_net2, layout = "circle", labels = labels,
       title = paste0("Contemporaneous (",penalty.names[Net2],")"),shape = "circle",
       vsize = 15, vsize2 = 20,
       mar = rep(8,4), asize = 2,
       threshold =1e-5,
       colFactor=1,
       title.cex = 1.5,
       theme = "TeamFortress",
       height=5,width=10)
box("figure")
title(expression(Omega^(mu)),cex.main = 2,line=1.5)
# Witdth = 1125 ; Height = 675

####### Others 
S.end <- L.end + k.dim * (k.dim + 1) * 0.5
Seps.param <- (S.end+1):(S.end +k.dim*p.dim)
Mod1.Seps.df<-data.frame(cbind(Seps.param,Mod1[Seps.param,]))

Mod1.Seps.gg <- ggplot(Mod1.Seps.df, aes(x = Seps.param, y = Est)) +
  geom_point(size = 1, color = "black") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.25, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(-1, 1))+
  labs(
    y = expression(Sigma^(epsilon)),
    x = "",
    title = "(b)"
  ) + theme(
    axis.text.y = element_blank(),     # remove tick labels
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank()     # remove axis title
  )


Lambda.param <- (k.dim^2+1):(k.dim^2+k.dim*p.dim)
Mod1.Lambda.df<-data.frame(cbind(Lambda.param,Mod1[Lambda.param,]))

Mod1.Lambda.gg <- ggplot(Mod1.Lambda.df[Mod1.Lambda.df$Est!=1.0000,], aes(x = Lambda.param, y = Est)) +
  geom_point(size = 1, color = "black") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.25, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +   # <<--- transpose the plot
  scale_y_continuous(limits = c(-0.05, 1.1))+
  labs(
    y = expression(Lambda),
    x = "",
    title = "                         AL-1 estimates with 95% Confidence Intervals"
  ) + theme(
    axis.text.y = element_blank(),     # remove tick labels
    axis.ticks.y = element_blank(),    # remove tick marks
    axis.title.y = element_blank(),     # remove axis title
    plot.title = element_text(hjust = 0.5)
  )



