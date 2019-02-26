
  install.packages("pROC")
  library(pROC)  

  install.packages("RCurl")                  
  require(RCurl)
   
  dat <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/logit.regression.rice/master/ricediversity.44k.germplasm3.csv"),header=T)   
  phe <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/logit.regression.rice/master/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"),header=T)


  dat <- read.csv("ricediversity.44k.germplasm3.csv",header=T)   
  phe <- read.table("RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",header=T,sep=",")
  inters <- intersect(dat$NSFTV.ID,phe$NSFTVID)
  dat <- dat[match(inters,dat$NSFTV.ID),]
  phe <- phe[match(inters,phe$NSFTVID),]  
  
  dat1 <- cbind(dat,phe)    
  dat2 <- dat1[which( dat1$Sub.population == "JAP" | dat1$Sub.population == "IND"),]  

buff <- NULL
res <- NULL
cnt <- 0
    for (kk in 1:100){

     pop <- sample(1:dim(dat2)[1],replace=F)   
     for (z in 0:9){   
        if (z < 9){
           st <- z * ceiling(dim(dat2)[1]/10) + 1
           en <- (z+1) * ceiling(dim(dat2)[1]/10)
        }   
        if (z == 9){
           st <- z * ceiling(dim(dat2)[1]/10) + 1
           en <- dim(dat2)[1]
         }
         tra_sam <- pop[-1* c(st:en)]
     
         dat22 <- dat2[tra_sam,]    
         fit <- glm(Sub.population ~ Panicle.number.per.plant + Seed.number.per.panicle + Florets.per.panicle + Panicle.fertility + Straighthead.suseptability + Blast.resistance + Protein.content,data=dat2,family=binomial(link = "logit"),control = list(maxit = 100))  
         summary(fit)


         logit <- function(x1, x2, x3, x4, x5, x6, x7){  
         res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] + x4*coefficients(fit)[5] + x5*coefficients(fit)[6] + x6*coefficients(fit)[7] + x7*coefficients(fit)[8] )))  
         return(res)
         }     

         val_sam <- c(1:dim(dat2)[1])[ -1 * tra_sam]
         logit.res <- logit(dat2$Panicle.number.per.plant[val_sam], dat2$Seed.number.per.panicle[val_sam], dat2$Florets.per.panicle[val_sam], dat2$Panicle.fertility[val_sam], dat2$Straighthead.suseptability[val_sam], dat2$Blast.resistance[val_sam], dat2$Protein.content[val_sam])

        rr <- cbind(logit.res,dat2$Sub.population[val_sam])
        sorted.logit.res <- rr[match(sort(logit.res),rr),]
        col.topo <- topo.colors(3)
        palette(col.topo)
        plot(sorted.logit.res[,1],col=sorted.logit.res[,2],pch=20)
        abline(h=0.5)


        total_tej <- which(dat2$Sub.population == "JAP")
        pred_tej  <- val_sam[which(logit.res >= .5)]
    
        total_ind <- which(dat2$Sub.population == "IND")
        pred_ind  <- val_sam[which(logit.res <  .5)]
        
        cnt <- cnt + 1
        res[cnt] <- ( length(intersect(total_tej,pred_tej)) + length(intersect(total_ind,pred_ind)) )/ length(which(!is.na(logit.res)))
      }
}

mean(res)




