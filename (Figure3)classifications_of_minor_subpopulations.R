

install.packages("RCurl")                  
require(RCurl)
   
dat <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/logit.regression.rice/master/ricediversity.44k.germplasm3.csv"),header=T)   
phe <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/logit.regression.rice/master/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"),header=T)


inters <- intersect(dat$NSFTV.ID,phe$NSFTVID)
dat <- dat[match(inters,dat$NSFTV.ID),]
phe <- phe[match(inters,phe$NSFTVID),]  
  
dat1 <- cbind(dat,phe)    
dat2 <- dat1[which( dat1$Sub.population == "JAP" | dat1$Sub.population == "IND"),]  

phen1 <- dat2$Panicle.number.per.plant
phen2 <- dat2$Seed.number.per.panicle
phen3 <- dat2$Florets.per.panicle
phen4 <- dat2$Panicle.fertility 
phen5 <- dat2$Straighthead.suseptability
phen6 <- dat2$Blast.resistance  
phen7 <- dat2$Protein.content
fit <- glm(dat2$Sub.population ~ phen1 + phen2 + phen3 + phen4 + phen5 + phen6 + phen7,data=dat2,family=binomial(link = "logit"),control = list(maxit = 100))
summary(fit) 

logit <- function(x1,x2,x3,x4,x5,x6){       
res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] + x4*coefficients(fit)[5] + x5*coefficients(fit)[6] + x6*coefficients(fit)[7] + x7*coefficients(fit)[8]))); return(res) 
       } 


## Aus validation

dat3 <- dat1[which( dat1$Sub.population == "AUS"),]  

phen1 <- dat3$Panicle.number.per.plant
phen2 <- dat3$Seed.number.per.panicle
phen3 <- dat3$Florets.per.panicle
phen4 <- dat3$Panicle.fertility 
phen5 <- dat3$Straighthead.suseptability
phen6 <- dat3$Blast.resistance    
fit <- glm(Sub.population ~ phen1 + phen2 + phen3 + phen4 + phen5 + phen6, data=dat2,family=binomial(link = "logit"))
summary(fit) 


logit <- function(x1,x2,x3,x4,x5,x6){       
res <- 1/(1 + exp( -1* (coefficients(fit)[1] + x1*coefficients(fit)[2] + x2*coefficients(fit)[3] + x3*coefficients(fit)[4] + x4*coefficients(fit)[5] + x5*coefficients(fit)[6] + x6*coefficients(fit)[7] ))); return(res) 
       } 

logit.res <- logit(phen1,phen2,phen3,phen4,phen5,phen6)

res.tab1 <- cbind(as.character(dat3[,4]),logit.res)[which(is.na(logit.res) == FALSE),]
res.tab2 <- res.tab1[which(as.numeric(res.tab1[,2]) < .5),]
res.tab3 <- res.tab1[which(as.numeric(res.tab1[,2]) > .5),]

plot(sort(as.numeric(res.tab1[,2])),col="red",pch=20,ylab="indica / japonica")
abline(h=0.5)


## Admix validation

dat4 <- dat1[which( dat1$Sub.population == "ADMIX"),]  
phen1 <- dat4$Panicle.number.per.plant
phen2 <- dat4$Seed.number.per.panicle
phen3 <- dat4$Florets.per.panicle
phen4 <- dat4$Panicle.fertility 
phen5 <- dat4$Straighthead.suseptability
phen6 <- dat4$Blast.resistance    

logit.res <- logit(phen1,phen2,phen3,phen4,phen5,phen6)

res.tab1 <- cbind(as.character(dat4[,4]),logit.res)[which(is.na(logit.res) == FALSE),]
res.tab2 <- res.tab1[which(as.numeric(res.tab1[,2]) < .5),]
res.tab3 <- res.tab1[which(as.numeric(res.tab1[,2]) > .5),]

plot(sort(as.numeric(res.tab1[,2])),col="red",pch=20,ylab="indica / japonica")
abline(h=0.5)


## Aroma validation

dat5 <- dat1[which( dat1$Sub.population == "AROMATIC"),]  
phen1 <- dat5$Panicle.number.per.plant
phen2 <- dat5$Seed.number.per.panicle
phen3 <- dat5$Florets.per.panicle
phen4 <- dat5$Panicle.fertility 
phen5 <- dat5$Straighthead.suseptability
phen6 <- dat5$Blast.resistance    

logit.res <- logit(phen1,phen2,phen3,phen4,phen5,phen6)

res.tab1 <- cbind(as.character(dat5[,4]),logit.res)[which(is.na(logit.res) == FALSE),]
res.tab2 <- res.tab1[which(as.numeric(res.tab1[,2]) < .5),]
res.tab3 <- res.tab1[which(as.numeric(res.tab1[,2]) > .5),]

plot(sort(res.tab1[,2]),col="red",pch=20,ylab="indica / japonica")
abline(h=0.5)


