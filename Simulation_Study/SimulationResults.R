library(MASS)
library(dplyr)

#################################
##### Dimension of the data
M.rep <- cbind(1:25,26:50,51:75,76:100,101:125,126:150,151:175,176:200)
N.ls<- c(75,100,125,150)
T.dim<- 50 #c(25,50,75)
p.ls<- c(5,10,15)
k.ls<- c(2,4,6)
miss <- 'y'#c('y','n')
item.3<-3 
item.2<-2

penalization <- c("_aLasso","_aLasso2","7.5","0")

penFIML.tab_condition.list<- Quant.tab_condition.list<-list()

for(pen in penalization){

    ### Nb of individual
    for(N.dim in N.ls){

        ### Number of items   
        for (p.dim in p.ls) {

            # N.dim<-150
            # T.dim<-50
            # p.dim<-15
            # k.dim<-6
            # miss<-'y'
            # pen <- penalization[1]
            # penalty.term<-penalty.ls[2]

            penFIML.tab_condition.list<- Quant.tab_condition.list<-list()

            for(g in 1:8) {
                
                # g=1
                # List of Tables per conditions
                ls.table.name <- paste0(
                              "Properties_Table",
                              "_Penalty",as.character(pen),
                              "_Rep",as.character(min(M.rep[,g])),
                              "_",as.character(max(M.rep[,g])),
                              "_N",as.character(N.dim),
                              "_T",as.character(T.dim),
                              "_p",as.character(p.dim),
                              "_NA_",miss,".csv") 

                penFIML.table.mean <- read.csv(ls.table.name)
                penFIML.tab_mean.list <- list(penFIML.table.mean)

                penFIML.tab_condition.list <- append(
                    penFIML.tab_condition.list,penFIML.tab_mean.list)

                # List of "Quantities" per conditions
                ls.Quant.name <- paste0(
                              "Quantities",
                              "_Penalty",as.character(pen),
                              "_Rep",as.character(min(M.rep)),
                              "_",as.character(max(M.rep)),
                              "_N",as.character(N.dim),
                              "_T",as.character(T.dim),
                              "_p",as.character(p.dim),
                              "_NA_",miss,".csv") 

                Quant.table.mean <- read.csv(ls.table.name)
                Quant.tab_mean.list <- list(Quant.table.mean)

                Quant.tab_condition.list <- append(
                    Quant.tab_condition.list,Quant.tab_mean.list)
            }  

            penFIML.table.sum.ls <- Reduce("+",penFIML.tab_condition.list)
            penFIML.table.mean.ls <- penFIML.table.sum.ls/8
            round(penFIML.table.mean.ls[,-1],3)
            penFIML.file_name <- paste0(
                "Brute_Results-Penalty",as.character(pen),
                "_N",as.character(N.dim),
                "_T",as.character(T.dim),
                "_p",as.character(p.dim),".csv"
            )
            write.csv(penFIML.table.mean.ls[,-1],file=penFIML.file_name)

            Quant.table.sum.ls <- Reduce("+",Quant.tab_condition.list)
            Quant.table.mean.ls <- Quant.table.sum.ls/8
            round(Quant.table.mean.ls,3)
            Quant.file_name <- paste0(
                "Brute_Quantities-Penalty",as.character(pen),
                "_N",as.character(N.dim),
                "_T",as.character(T.dim),
                "_p",as.character(p.dim),".csv"
            )
            write.csv(Quant.table.mean.ls[,-1],file=Quant.file_name)
        }
    }
}
