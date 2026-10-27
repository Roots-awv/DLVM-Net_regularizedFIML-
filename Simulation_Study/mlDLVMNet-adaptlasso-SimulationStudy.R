# library(optimParallel)
library(parallel)
library(MASS)
library(Matrix)
library(dplyr)

############################# Data preparation / Statistics summary 
# Loading the data / methods
source("../Computing_Penalized_Pseudo_FIML_LargeBT_CovMatrix.R")

y.simu.ls <- readRDS("Simulated_Data_rep200.rds") 

##### Dimension of the data
M.rep <- 1:25
N.ls<- c(75,100,125,150)
T.dim<- 25 # c(25,50,75)
p.ls<- c(5,10,15)
k.ls<- c(2,4,6)
miss <- 'y' #c('y','n')
item.3<-3 
item.2<-2
penalty.term<- 0.25#c(0,sqrt(mean(N.ls)/2))
cov.v <- 3

penFIML.tab_condition.list<- Quant.tab_condition.list<-list()


### Nb of individual
for(N.dim in N.ls){

    ### Number of items   
    for (p.dim in p.ls) {

        # N.dim<-100
        # T.dim<-50
        # p.dim<-10
        # k.dim<-4
        # miss<-'y'
        # penalty.term<-penalty.ls[2]
        
        penFIML.tab_replic.list<-Quant.tab_replic.list<-list()
  
        ### Replication
        for(m.rep in M.rep){
          
          # m.rep<-10
          # N.dim<-75
          # T.dim<-50
          # p.dim<-10
          # k.dim<-4
          # miss<-'y'
          # penalty.term<-penalty.ls[2]
          
          # Dimension of the latent space
          if (p.dim == p.ls[1]){
              k.dim = k.ls[1]
              k1.c = c(1)
              k2.c = c(2) 
          } else if (p.dim == p.ls[2]){
              k.dim = k.ls[2] 
              k1.c = c(1,4)
              k2.c = 2:3
              df.free = p.dim - 2
          } else if (p.dim == p.ls[3]){
              k.dim = k.ls[3] 
              k1.c = c(1,2,6)
              k2.c = 3:5
          }
          
          ### Loading files 
          lsdata.name <- paste0(
                              #"Simulateddata_files/",
                              "y_var",
                              "_M",as.character(m.rep),
                              "_N",as.character(N.dim),
                              "_T",as.character(T.dim),
                              "_p",as.character(p.dim),
                              "_CovCode",as.character(cov.v),
                              "_NA_",miss)
          y.var <- y.simu.ls[[lsdata.name]]
          dim(y.var)
          ### Vector of parameters
          fltheta<-paste0("theta_p_",as.character(p.dim),
                        "_CovCode",as.character(cov.v),".csv")
          theta <- read.csv(fltheta)
          # theta$x <- theta%>%round(4)
          
          # Summary statistics
          Z.summary<-summary.stats(y.var)
          
          ####################### Optimization
          # Combine into parameter vector
          r.dim <- p.dim
          m.dim <- k.dim
          
          # Parameters initialization & boundaries
          B<-diag(0.5*runif(m.dim,-1,1))
          min.B <- rep(-1,m.dim^2) 
          ini.B <- B%>%as.vector  
          max.B <- rep(1,m.dim^2) 
          
          #
          min.ly <-vector("list", m.dim)
          ini.ly <-vector("list", m.dim)
          max.ly <-vector("list", m.dim)
          for (k1 in k1.c){
            min.ly[[k1]]<-c(1,-2,-2)
            ini.ly[[k1]]<-rep(1,3)
            max.ly[[k1]]<-c(1,2,2)
          }
          
          for (k2 in k2.c){
            min.ly[[k2]] <- c(1,-2)
            ini.ly[[k2]] <- rep(1,2)
            max.ly[[k2]] <- c(1,2)
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
          est.file <- paste0(
                        "unregularized_FIML/estimate",
                        "_Penalty0",
                        "_M",as.character(m.rep),
                        "_N",as.character(N.dim),
                        "_T",as.character(T.dim),
                        "_p",as.character(p.dim),
                        "_NA_",miss,
                        ".csv")
          FIML.est<-read.csv(est.file)
          est.inv <- 1/FIML.est[1:length(theta0),'x']
          inf.vals <- ifelse(est.inv == Inf,0,1)
          weights <- inf.vals*est.inv
          weights[is.na(weights)]<- max(weights,na.rm=T)
          lambda.w <- abs(weights)*penalty.term #penalty.tuning
          lambda.w <- lambda.w%>%as.vector
          # lambda.w[1:(m.dim^2)] <- max(weights,na.rm=T)
          print(round(lambda.w,2))
           
          # Optimize Full-Information ML using L-BFGS-B
          start_time <- proc.time()
          lBFGSb.penFIML <- optim(
              theta0, 
              fn=  log.penFIMLit.Y, 
              gr = ScoreY.penFIMLit, 
              method = "L-BFGS-B",
              #hessian = TRUE,
              data = Z.summary,  
              m.dim = m.dim,
              r.dim = r.dim,
              N.dim = N.dim,
              T.dim = T.dim,
              pen.l = lambda.w,
              lower = lower.b,
              upper = upper.b,
              control = list(maxit = 500, factr = 1e7)
          )
          end_time <- proc.time()
          elapsed <- end_time - start_time
          print(lBFGSb.penFIML$message)
          print(elapsed)
          # Parameter Estimates 
          penFIML.est.lBFGSb <- lBFGSb.penFIML$par
          Sigma_u.est <- Sigma_u.hat(lBFGSb.penFIML$par,Z.summary$mean,m.dim,r.dim,lambda.w[1])
          Sigma_u.vec_est <- Sigma_u.est[lower.tri(Sigma_u.est,diag=T)]
          # remove_u<- length(penFIML.est.lBFGSb) - length(Sigma_u.vec_est) [1:remove_u]
          penFIML.est.lBFGSb.all <- c(penFIML.est.lBFGSb,Sigma_u.vec_est)
          round(penFIML.est.lBFGSb.all,3)
          # Bias and Deviation
          true_values<-theta$x
          penFIML.lBFGSb.bias <- penFIML.est.lBFGSb.all-true_values
          penFIML.lBFGSb.abs_bias <- abs(penFIML.lBFGSb.bias)
          penFIML.lBFGSb.sqrd_bias <- penFIML.lBFGSb.bias^2

          # Asymptotic inference
          # Invert
          Hess.theta<- penFIML.est.lBFGSb%>%length%>%diag
          Grad.theta <- ScoreY.penFIMLit(penFIML.est.lBFGSb,Z.summary,k.dim,p.dim,N.dim,T.dim,pen.l=lambda.w)
          Hcal.theta<-Grad.theta%*%t(Grad.theta)  
          inv.Hess <- 2*solve(Hess.theta)
          dim(Hess.theta)
          all(eigen(Hess.theta, symmetric = TRUE)$values > 0)
          
          # Variance-covariance matrix
          As.Var <- inv.Hess%*%Hcal.theta%*%t(inv.Hess)/(N.dim*T.dim)#2*inv.Hess/(N.dim*T.dim)
          max(As.Var)
          min(As.Var)
          dim(As.Var)

          # Standard errors
          SE_ <- sqrt(diag(As.Var))
          l_CI95<-penFIML.est.lBFGSb+qnorm(0.025)*SE_
          u_CI95<-penFIML.est.lBFGSb+qnorm(0.975)*SE_
          #cbind(penFIML.est.lBFGSb,SE_,l_CI95,u_CI95)%>%round(4)

          # Sensitivity and Specificity
          all.est<-Param.ls(penFIML.est.lBFGSb,m.dim,r.dim)          
          all.l_CI95<-Param.ls(l_CI95,m.dim,r.dim)
          all.u_CI95<-Param.ls(u_CI95,m.dim,r.dim)
          all.true<-Param.ls(theta$x,m.dim,r.dim)
          
          # 
          B.pen <- all.est$B%>%as.vector
          B.l_CI95 <- all.l_CI95$B%>%as.vector
          B.u_CI95 <- all.u_CI95$B%>%as.vector
          B.true <- all.true$B%>%as.vector
          # 
          Lambda.pen <- all.est$Lambda[which(abs(all.est$Lambda)!=0&abs(all.est$Lambda)!=1)]
          Lambda.l_CI95 <- all.l_CI95$Lambda[which(abs(all.est$Lambda)!=0&abs(all.est$Lambda)!=1)]
          Lambda.u_CI95 <- all.u_CI95$Lambda[which(abs(all.est$Lambda)!=0&abs(all.est$Lambda)!=1)]
          Lambda.true <- all.true$Lambda[which(abs(all.est$Lambda)!=0&abs(all.est$Lambda)!=1)]
          # 
          Sig.zeta.pen <- all.est$Sig.zeta[lower.tri(all.est$Sig.zeta,diag=T)]
          Sig.zeta.l_CI95 <- all.l_CI95$Sig.zeta[lower.tri(all.l_CI95$Sig.zeta,diag=T)]
          Sig.zeta.u_CI95 <- all.u_CI95$Sig.zeta[lower.tri(all.u_CI95$Sig.zeta,diag=T)]
          Sig.zeta.true <- all.true$Sig.zeta[lower.tri(all.true$Sig.zeta,diag=T)]
          # 
          Sig.eps.pen <- diag(all.est$Sig.eps)
          Sig.eps.l_CI95 <- diag(all.l_CI95$Sig.eps)
          Sig.eps.u_CI95 <- diag(all.u_CI95$Sig.eps)
          Sig.eps.true <- diag(all.true$Sig.eps)
          # 
          penalized.est<-c(B.pen,Lambda.pen,Sig.zeta.pen,Sig.eps.pen)
          lower.CI95<-c(B.l_CI95,Lambda.l_CI95,Sig.zeta.l_CI95,Sig.eps.l_CI95)
          upper.CI95<-c(B.u_CI95,Lambda.u_CI95,Sig.zeta.u_CI95,Sig.eps.u_CI95)
          true_values_reg <- c(B.true,Lambda.true,Sig.zeta.true,Sig.eps.true)

          # True zeros and true non-zeros 
          true_nonzero_indices <- which(true_values_reg != 0)
          true_zero_indices <- which(true_values_reg == 0)

          # Cutoff values based        
          estimated_zero_indices.0.01 <- which(abs(penalized.est)<0.01)
          estimated_nonzero_indices.0.01 <- which(abs(penalized.est)>=0.01)
  
          TP.0.01 <- sum(estimated_nonzero_indices.0.01 %in% true_nonzero_indices)          
          FP.0.01 <- sum(estimated_nonzero_indices.0.01 %in% true_zero_indices)
          TN.0.01 <- sum(estimated_zero_indices.0.01 %in% true_zero_indices)
          FN.0.01 <- sum(estimated_zero_indices.0.01 %in% true_nonzero_indices)
              
          FPR.0.01 <- FP.0.01/(FP.0.01+TN.0.01)          
          TNR.0.01 <- TN.0.01/(FP.0.01+TN.0.01)
          FDR.0.01 <- FP.0.01/(FP.0.01+TP.0.01)

          # CI 95% based  
          estimated_zero_indices.CI95 <- which(lower.CI95 <= 0 & upper.CI95 >= 0)
          estimated_nonzero_indices.CI95 <- which(lower.CI95 >= 0 | upper.CI95 <= 0)

          TP.CI95 <- sum(estimated_nonzero_indices.CI95 %in% true_nonzero_indices)          
          FP.CI95 <- sum(estimated_nonzero_indices.CI95 %in% true_zero_indices)
          TN.CI95 <- sum(estimated_zero_indices.CI95 %in% true_zero_indices)
          FN.CI95 <- sum(estimated_zero_indices.CI95 %in% true_nonzero_indices)
          
          FPR.CI95 <- FP.CI95/(FP.CI95+TN.CI95)          
          TNR.CI95 <- TN.CI95/(FP.CI95+TN.CI95)
          FDR.CI95 <- FP.CI95/(FP.CI95+TP.CI95)

          # Table
          lBFGSb.penFIML.results <- cbind(
                        true_values,
                        penFIML.est.lBFGSb.all,
                        # As.std_err.theta,
                        penFIML.lBFGSb.bias,
                        penFIML.lBFGSb.abs_bias,
                        penFIML.lBFGSb.sqrd_bias,
                        FPR.0.01,TNR.0.01,FDR.0.01,
                        FPR.CI95,TNR.CI95,FDR.CI95)
          colnames(lBFGSb.penFIML.results) <- c(
            "Theta","Est.","Bias","aBias","SSE",
            "FPR.0.01","TNR.0.01","FDR.0.01",
            "FPR.CI95","TNR.CI95","FDR.CI95"
          )
          
          ## Other quantities
          df.free <- sum(abs(penalized.est) > 1e-2)
          unreg.neglL<-N.dim*T.dim*log.penFIMLit.Y(penFIML.est.lBFGSb,Z.summary,m.dim,r.dim, N.dim, T.dim, pen.l=0)
          Neg.logLik<-lBFGSb.penFIML$value
          AIC <- 2*(unreg.neglL + df.free)
          BIC <- 2*unreg.neglL + df.free*log(N.dim*T.dim)
          EBIC <- BIC + log(choose(length(penalized.est),df.free))
          iter <-lBFGSb.penFIML$counts
          CVce<-lBFGSb.penFIML$convergence
          Quantities.penFIML <- c(
                Neg.logLik,unreg.neglL,AIC,BIC,EBIC,iter['gradient'],CVce,elapsed['elapsed'])
          names(Quantities.penFIML)<-c(
            "F_reg","F_ML","AIC","BIC","EBIC","gradient","Convergence","Elapsed")
          t(Quantities.penFIML)
          ## List of replication's results
          lBFGSb.penFIML.ls<-list(lBFGSb.penFIML.results)
          penFIML.tab_replic.list <- append(penFIML.tab_replic.list,lBFGSb.penFIML.ls)
          
          Quantities.penFIML.ls<-list(Quantities.penFIML)
          Quant.tab_replic.list<-append(Quant.tab_replic.list,Quantities.penFIML.ls)
          
        }
        
          # List of Tables per conditions
          ls.table.name <- paste0(
                        "Properties_Table_plus",
                        "_Penalty_aLasso3",
                        "_Rep",as.character(min(M.rep)),
                        "_",as.character(max(M.rep)),
                        "_N",as.character(N.dim),
                        "_T",as.character(T.dim),
                        "_p",as.character(p.dim),
                        "_NA_",miss,".csv") 

          penFIML.table.sum <- Reduce("+",penFIML.tab_replic.list)
          penFIML.table.mean <- penFIML.table.sum/length(M.rep)
          
          write.csv(penFIML.table.mean, file = ls.table.name)
          
          # List of "Quantities" per conditions
          ls.Quant.name <- paste0(
                        "Quantities",
                        "_Penalty_aLasso3",
                        "_Rep",as.character(min(M.rep)),
                        "_",as.character(max(M.rep)),
                        "_N",as.character(N.dim),
                        "_T",as.character(T.dim),
                        "_p",as.character(p.dim),
                        "_NA_",miss,".csv") 
                        
          Quant.table.sum <- Reduce("+",Quant.tab_replic.list)
          Quant.table.mean <- Quant.table.sum/length(M.rep)
          
          write.csv(Quant.table.mean, file = ls.Quant.name)

    }
}
#


# #### Saving the results
# saveRDS(penFIML.tab_condition.list, file = "Properties_Table_lasso1_rep1_rep40.rds")
# saveRDS(Quant.tab_condition.list, file = "OtherQuantities_lasso1_rep1_rep40.rds")
# 
# #### Read results 
# penFIML.tab_condition.list<-readRDS("Properties_Table_lasso1_rep1_rep40.rds")
# Quant.tab_condition.list<-readRDS("OtherQuantities_lasso1_rep1_rep40.rds")


