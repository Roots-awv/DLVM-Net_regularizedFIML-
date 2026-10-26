######################### Data generation process
library(MASS)
library(Matrix)
library(dplyr)
library(mvtnorm)

### Random number generationg
set.seed(070225)

##### Dimension of the data
M.rep <- 200
N.ls<- c(75,100,125,150)
T.ls<- c(25,50)
p.ls<- c(5,10,15)
k.ls<- c(2,4,6)
miss <- 'y'
cov.ls <- c(1,3)
item.3<-3 
item.2<-2

dgp.list <- list()

### Replication
for(m.rep in 1:M.rep){

  ### Nb of individual
  for(N.dim in N.ls){
        
    ### Length of time points
    for (T.dim in T.ls) {

      ### Number of items   
      for (p.dim in p.ls) {

        ### Covariance 
        for (cov.v in cov.ls) {
          
          # N.dim<-100
          # T.dim<-50
          # p.dim<-10
          # k.dim<-4
          # miss<-'y'
          
          # Dimension of the latent space
          if (p.dim == p.ls[1]){
            k.dim = k.ls[1]
            k1.c = c(1)
            k2.c = c(2)
            cov_ = matrix(1,k.dim,k.dim)
          } else if (p.dim == p.ls[2]){
            k.dim = k.ls[2] 
            k1.c = c(1,4)
            k2.c = 2:3
            cov_ = matrix(c(1,1,0,0,
                            1,1,0,1,
                            0,0,1,1,
                            0,1,1,1),
                      k.dim,k.dim)
          }  else if (p.dim == p.ls[3]){
            k.dim = k.ls[3] 
            k1.c = c(1,2,6)
            k2.c = 3:5
            cov_ = matrix(c(1,1,0,0,0,1,
                            1,1,0,1,0,0,
                            0,0,1,1,0,1,
                            0,1,1,1,0,0,
                            0,0,0,0,1,1,
                            1,0,1,0,1,1),
                      k.dim,k.dim)
          }
                    
          #### Parameters / population values
          ## For Measurement / Factor model
          # Intercept
          tau<-rep(0,p.dim)

          # loading factor
          lambda <-vector("list", k.dim)
          for (k1 in k1.c){lambda[[k1]]<-rep(1,item.3)}
          for (k2 in k2.c){lambda[[k2]]<-rep(1,item.2)}
          Lambda <- lambda %>% bdiag %>% as.matrix
          
          # Variance
          Sigma_eps <- diag(0.25,p.dim)
          
          ## For Structural model
          # Intercept
          Sigma_u <- 0.125*(diag(k.dim)* cov.v + cov_)
          mu.eta <- rmvnorm(N.dim,rep(0,k.dim),Sigma_u)
          
          # Slopes
          b <- rep(0.8,k.dim)#k.dim^2 %>% runif(-.8,.8) %>% round(1) 
          B <- b %>% diag(k.dim)#matrix(k.dim,k.dim)
          # B <- b %>% diag %>% diag#as.matrix
          # B[2,1]<--0
          
          # Variance
          Sigma_eta <-0.18*(diag(k.dim)+cov_)
          # 0.25/sqrt(0.75*0.75) 
          
          #### Structural model as VAR(1)
          zeta <- array(0,c(k.dim,N.dim,T.dim)) # resilduals
          eta <- array(0,c(k.dim,N.dim,T.dim)) # latent variables
          for(i in 1:N.dim){
            zeta[,i,1] <- rmvnorm(1,rep(0,k.dim),Sigma_eta)
            eta[,i,1] <- zeta[,i,1]
            for(t in 2:T.dim){
              zeta[,i,t] <- rmvnorm(1,rep(0,k.dim),Sigma_eta)
              eta[,i,t] <- B %*% (eta[,i,t-1]) + zeta[,i,t]
            }
          }
          
          #### Dynamic Factor model
          eps <- array(0,c(p.dim,N.dim,T.dim)) # residuals
          y.var <- array(0,c(p.dim,N.dim,T.dim)) # observed items
          for(i in 1:N.dim){
            for(t in 1:T.dim){
              eps[,i,t] <- rmvnorm(1,rep(0,p.dim),Sigma_eps)
              y.var[,i,t] <- tau + Lambda %*%mu.eta[i,]+ Lambda %*% eta[,i,t] + eps[,i,t]
            }
          }
          
          # Missing data
          N.miss<-sample(1:N.dim,N.dim/5)
          T.miss<-sample(1:T.dim,T.dim/5)
          p.miss<-sample(1:p.dim,3)
          if (miss=='y'){
            y.var[p.miss,N.miss,T.miss]<-NA
          } else{
            y.var <- y.var 
          }

          ### Saving files
          fldata.name <- paste0(
                              #"Simulateddata_files/",
                              "y_var",
                              "_M",as.character(m.rep),
                              "_N",as.character(N.dim),
                              "_T",as.character(T.dim),
                              "_p",as.character(p.dim),
                              "_CovCode",as.character(cov.v),
                              "_NA_",miss)
          # write.csv(y.var, file = fldata.name)
          y.var.list<-list(y.var)
          names(y.var.list)<-fldata.name
          dgp.list<-append(dgp.list,y.var.list)
          # print(fldata.name)
                  
        ### Vector of parameters
        theta <- c(B,Lambda,
                   Sigma_eta[lower.tri(Sigma_eta,diag=T)],
                   Sigma_eps[lower.tri(Sigma_eps,diag=T)],
                   Sigma_u[lower.tri(Sigma_u,diag=T)])
        fltheta<-paste0("theta_p_",as.character(p.dim),
                        "_CovCode",as.character(cov.v),".csv")
        write.csv(theta, file = fltheta)
        
        print(paste(m.rep, N.dim, p.dim, miss,cov.v))
        }
      }
    }
  }
}


saveRDS(dgp.list, file = "Simulated_Data_rep200.rds")
# names(dgp.list)
# dim(dgp.list[["y_var_M5_N150_T50_p10_NA_y"]])
