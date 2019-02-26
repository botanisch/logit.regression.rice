
  install.packages("pROC")
  library(pROC)  

  install.packages("RCurl")                  
  require(RCurl)
   
  dat <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/logit.regression.rice/master/ricediversity.44k.germplasm3.csv"),header=T)   
  phe <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/logit.regression.rice/master/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"),header=T)
  inters <- intersect(dat$NSFTV.ID,phe$NSFTVID)
  dat <- dat[match(inters,dat$NSFTV.ID),]
  phe <- phe[match(inters,phe$NSFTVID),]  
  
  dat1 <- cbind(dat,phe)    
  dat2 <- dat1[which( dat1$Sub.population == "JAP" | dat1$Sub.population == "IND"),]  



#### one predictors    

 temp  <- NULL
 cnt   <- 0
 lines <- NULL

 for (j in 15:38){               
                                 
    phen1 <- dat2[,j] 
    fit <- glm(Sub.population ~ phen1 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1){  

      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] )))  
      return(res)
    }     
    logit.res <- logit(phen1)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    res <- roc(response, logit.res, plot=TRUE, direction="<")  

    cnt       <- cnt + 1
    temp[cnt] <- res$auc[1]
    lines[cnt] <- paste(res$auc[1], colnames(dat2)[j],sep="\t")
    cat(lines[cnt])
    cat("\n")     
 }

#  Classification accuracy with the parameters yielding the highest AUC
    phen1 <- dat2$Panicle.number.per.plant
    fit <- glm(Sub.population ~ phen1 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1){  
      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] )))  
      return(res)
    }     
    logit.res <- logit(phen1)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],ann=FALSE,col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    AUC <- roc(response, logit.res, plot=TRUE, direction="<")  
    AUC



#### two predictors   
 temp  <- NULL
 cnt   <- 0
 lines <- NULL

 for (j in 15:37){
 for (k in (j+1):38){                   
                                 
    phen1 <- dat2[,j]
    phen2 <- dat2[,k]
   
    fit <- glm(Sub.population ~ phen1 + phen2 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1,x2){  

      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] )))  
      return(res)
    }     
    logit.res <- logit(phen1,phen2)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)
     
    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    res <- roc(response, logit.res, plot=TRUE, direction="<") 

    cnt       <- cnt + 1
    temp[cnt] <- res$auc[1]
    lines[cnt] <- paste(res$auc[1], colnames(dat2)[j], colnames(dat2)[k],sep="\t")
    cat(lines[cnt])
    cat("\n") 
 }}

#  Classification accuracy with the parameters yielding the highest AUC
    phen1 <- dat2$Panicle.number.per.plant
    phen2 <- dat2$Brown.rice.seed.width
    
    fit <- glm(Sub.population ~ phen1 + phen2 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1, x2){  
      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] )))  
      return(res)
    }     
    logit.res <- logit(phen1,phen2)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],ann=FALSE,col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    AUC <- roc(response, logit.res, plot=TRUE, direction="<")  
    AUC


      
#### three predictors   
 temp  <- NULL
 cnt   <- 0
 lines <- NULL

 for (j in 15:36){
 for (k in (j+1):37){ 
 for (m in (k+1):38){                     
                                 
    phen1 <- dat2[,j]
    phen2 <- dat2[,k]
    phen3 <- dat2[,m]

   
    fit <- glm(Sub.population ~ phen1 + phen2 + phen3,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1,x2,x3){   

    res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] )))   
    return(res)
    }     
    logit.res <- logit(phen1,phen2,phen3)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    res <- roc(response, logit.res, plot=TRUE, direction="<")  

    cnt       <- cnt + 1
    temp[cnt] <- res$auc[1]
    lines[cnt] <- paste(res$auc[1], colnames(dat2)[j], colnames(dat2)[k], colnames(dat2)[m], sep="\t")
    cat(lines[cnt])
    cat("\n") 
 }}}

#  Classification accuracy with the parameters yielding the highest AUC  
    phen1 <- dat2$Panicle.number.per.plant
    phen2 <- dat2$Straighthead.suseptability
    phen3 <- dat2$Blast.resistance
        
    fit <- glm(Sub.population ~ phen1 + phen2 + phen3 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1, x2, x3){  
      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3]  + x3*coefficients(fit)[4] )))  
      return(res)
    }     
    logit.res <- logit(phen1,phen2,phen3)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],ann=FALSE,col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    AUC <- roc(response, logit.res, plot=TRUE, direction="<")  
    AUC



#### four predictors   
 temp  <- NULL
 cnt   <- 0
 lines <- NULL

 for (i in 15:35){
 for (j in (i+1):36){
 for (k in (j+1):37){ 
 for (m in (k+1):38){                     
                                 
    phen1 <- dat2[,i]
    phen2 <- dat2[,j]
    phen3 <- dat2[,k]
    phen4 <- dat2[,m]
    
    tryCatch( { 
      fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4,data=dat2,family=binomial(link = "logit"))
      summary(fit) 
    }, error = function(e) { cat("Error happen\n") })
    
    tryCatch( { 
       logit <- function(x1,x2,x3,x4){       
       res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] + x4*coefficients(fit)[5] ))); return(res) 
       }     
       logit.res <- logit(phen1,phen2,phen3,phen4)      
    }, error = function(e) { cat("Error happen\n") })

   tryCatch( { 
   response <- rep(0,dim(dat2)[1])
   response[ which(dat2$Sub.population == "JAP") ] <- 1
   response  <- response[-1 * which(is.na(logit.res) == T)]
   logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
   res <- roc(response, logit.res, plot=TRUE, direction="<")     
   }, error = function(e) { cat("Error happen\n") })

   cnt       <- cnt + 1
   temp[cnt] <- res$auc[1]
   lines[cnt] <- paste(res$auc[1], colnames(dat2)[i], colnames(dat2)[j], colnames(dat2)[k], colnames(dat2)[m], sep="\t")
   cat(lines[cnt])
   cat("\n") 
 }}}}

#  Classification accuracy with the parameters yielding the highest AUC   
    phen1 <- dat2$Panicle.number.per.plant
    phen2 <- dat2$Brown.rice.volume
    phen3 <- dat2$Straighthead.suseptability
    phen4 <- dat2$Blast.resistance
        
    fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1, x2, x3, x4){  
      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4]  + x4*coefficients(fit)[5] )))  
      return(res)
    }     
    logit.res <- logit(phen1,phen2,phen3,phen4)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],ann=FALSE,col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    AUC <- roc(response, logit.res, plot=TRUE, direction="<")  
    AUC



#### five predictors 
 temp  <- NULL
 cnt   <- 0
 lines <- NULL
 res   <- NULL
 for (i in 15:34){
 for (j in (i+1):35){
 for (k in (j+1):36){ 
 for (m in (k+1):37){
 for (n in (m+1):38){
                      
                                 
    phen1 <- dat2[,i]
    phen2 <- dat2[,j]
    phen3 <- dat2[,k]
    phen4 <- dat2[,m]
    phen5 <- dat2[,n]
        
    fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4 + phen5, data=dat2,family=binomial(link = "logit"))
    summary(fit) 

    logit <- function(x1,x2,x3,x4,x5){       
    res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] + x4*coefficients(fit)[5] + x5*coefficients(fit)[6] ))); return(res) 
    }     
    logit.res <- logit(phen1,phen2,phen3,phen4,phen5)      

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    res <- roc(response, logit.res, plot=TRUE, direction="<")     

    cnt       <- cnt + 1
    temp[cnt] <- res$auc[1]
    lines[cnt] <- paste(res$auc[1], colnames(dat2)[i], colnames(dat2)[j], colnames(dat2)[k], colnames(dat2)[m], colnames(dat2)[n], sep="\t")
    cat(lines[cnt])
    cat("\n") 
 }}}}}

#  Classification accuracy with the parameters yielding the highest AUC
    phen1 <- dat2$Panicle.number.per.plant
    phen2 <- dat2$Brown.rice.volume
    phen3 <- dat2$Straighthead.suseptability
    phen4 <- dat2$Blast.resistance
    
        
    fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1, x2, x3, x4){  
      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4]  + x4*coefficients(fit)[5] )))  
      return(res)
    }     
    logit.res <- logit(phen1,phen2,phen3,phen4)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],ann=FALSE,col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    AUC <- roc(response, logit.res, plot=TRUE, direction="<")  
    AUC



#### six predictors   

 temp  <- NULL
 cnt   <- 0
 lines <- NULL
 res   <- NULL
 for (i in 15:33){
 for (j in (i+1):34){
 for (k in (j+1):35){ 
 for (m in (k+1):36){
 for (n in (m+1):37){
 for (p in (n+1):38){ 
                                 
   phen1 <- dat2[,i]
   phen2 <- dat2[,j]
   phen3 <- dat2[,k]
   phen4 <- dat2[,m]
   phen5 <- dat2[,n]
   phen6 <- dat2[,p]       

   fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4 + phen5 + phen6, data=dat2,family=binomial(link = "logit"))
   summary(fit) 


   logit <- function(x1,x2,x3,x4,x5,x6){       
   res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] + x4*coefficients(fit)[5] + x5*coefficients(fit)[6] + x6*coefficients(fit)[7] ))); return(res) 
   }     
   logit.res <- logit(phen1,phen2,phen3,phen4,phen5,phen6)      

   response <- rep(0,dim(dat2)[1])
   response[ which(dat2$Sub.population == "JAP") ] <- 1
   response  <- response[-1 * which(is.na(logit.res) == T)]
   logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
   res <- roc(response, logit.res, plot=TRUE, direction="<")     


   cnt       <- cnt + 1
   temp[cnt] <- res$auc[1]
   lines[cnt] <- paste(res$auc[1], colnames(dat2)[i], colnames(dat2)[j], colnames(dat2)[k], colnames(dat2)[m], colnames(dat2)[n], colnames(dat2)[p], sep="\t")
   cat(lines[cnt])
   cat("\n") 
 }}}}}}

#  Classification accuracy with the parameters yielding the highest AUC    
    phen1 <- dat2$Panicle.number.per.plant
    phen2 <- dat2$Seed.number.per.panicle
    phen3 <- dat2$Florets.per.panicle
    phen4 <- dat2$Panicle.fertility 
    phen5 <- dat2$Straighthead.suseptability
    phen6 <- dat2$Blast.resistance


    fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4 + phen5 + phen6 ,data=dat2,family=binomial(link = "logit"))  
    summary(fit)
    logit <- function(x1, x2, x3, x4, x5, x6){  
      res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4]  + x4*coefficients(fit)[5]  + x5*coefficients(fit)[6]  + x6*coefficients(fit)[7] )))  
      return(res)
    }     
    logit.res <- logit(phen1,phen2,phen3,phen4,phen5,phen6)

    plot(sort(logit.res))
    rr <- cbind(logit.res,dat2$Sub.population)
    sorted.logit.res <- rr[match(sort(logit.res),rr),]
    col.topo <- topo.colors(3)
    palette(col.topo)
    plot(sorted.logit.res[,1],ann=FALSE,col=sorted.logit.res[,2],pch=20,main=colnames(dat2[j]))
    abline(h=0.5)

    response <- rep(0,dim(dat2)[1])
    response[ which(dat2$Sub.population == "JAP") ] <- 1
    response  <- response[-1 * which(is.na(logit.res) == T)]
    logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
    AUC <- roc(response, logit.res, plot=TRUE, direction="<")  
    AUC



#### seven predictors   

 temp  <- NULL
 cnt   <- 0
 lines <- NULL
 res   <- NULL
 for (i in 15:32){
 for (j in (i+1):33){
 for (k in (j+1):34){ 
 for (m in (k+1):35){
 for (n in (m+1):36){
 for (p in (n+1):37){ 
 for (q in (p+1):38){                     
            
   phen1 <- dat2[,i]
   phen2 <- dat2[,j]
   phen3 <- dat2[,k]
   phen4 <- dat2[,m]
   phen5 <- dat2[,n]
   phen6 <- dat2[,p]       
   phen7 <- dat2[,q]       


   fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4 + phen5 + phen6 + phen7, data=dat2,family=binomial(link = "logit"))
   summary(fit) 


   logit <- function(x1,x2,x3,x4,x5,x6,x7){       
   res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] + x4*coefficients(fit)[5] + x5*coefficients(fit)[6] + x6*coefficients(fit)[7]  + x7*coefficients(fit)[8] ))); return(res) 
   }     
   logit.res <- logit(phen1,phen2,phen3,phen4,phen5,phen6,phen7)      

   response <- rep(0,dim(dat2)[1])
   response[ which(dat2$Sub.population == "JAP") ] <- 1
   response  <- response[-1 * which(is.na(logit.res) == T)]
   logit.res <- logit.res[-1 * which(is.na(logit.res) == T)]
   res <- roc(response, logit.res, plot=TRUE, direction="<")     


   cnt       <- cnt + 1
   temp[cnt] <- res$auc[1]
   lines[cnt] <- paste(res$auc[1], colnames(dat2)[i], colnames(dat2)[j], colnames(dat2)[k], colnames(dat2)[m], colnames(dat2)[n], colnames(dat2)[p], colnames(dat2)[q], sep="\t")
   cat(lines[cnt])
   cat("\n") 
 }}}}}}}



