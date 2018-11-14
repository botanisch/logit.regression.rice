   
   install.packages("RCurl")                  
   require(RCurl)
   
   dat1 <- read.csv(text=getURL("https://raw.githubusercontent.com/bongsongkim/logit.regression.rice/master/IBS_matrix.csv"),row.names=1)   
   dat2 <- matrix(rep(2,dim(dat1)[1]^2),nrow=dim(dat1)[1]) 
   dat  <- dat2 - dat1

   germ <- read.csv("RiceDiversity.44K.germplasm.csv",header=T)
   subs <- germ[which(germ[,8] == "TEJ" | germ[,8] == "TRJ" | germ[,8] == "IND"),3]    
   subs <- paste("NSFTV_",subs,sep="")

   dat <- dat[which(colnames(dat) %in% subs),which(colnames(dat) %in% subs)]

   germ <- read.csv("RiceDiversity.44K.germplasm.csv",header=T)
   germ <- germ[which(germ[,8] == "TEJ" | germ[,8] == "TRJ" | germ[,8] == "IND"),]
   germ[,3] <- paste("NSFTV_",germ[,3],sep="")

   germ_lbl <- germ[,3]                              
   sorted   <- germ[match(colnames(dat), germ_lbl),] 
   names    <- paste("(",germ[match(colnames(dat), germ_lbl),8],") ",germ[match(colnames(dat), germ_lbl),4],sep="")
   colnames(dat) <- names
   dat <- t(dat)
   colnames(dat) <- names
   dat[upper.tri(dat)] <- NA  
   dat <- as.dist(dat, diag = TRUE)
   hc  <- hclust(dat) 
   plot(as.dendrogram(hc),horiz=T)

   jpeg("ind.jap.dendrogram.jpg",width = 1300, height = 5000) 
   par(lwd=5,cex=0.9, mar=c(3,0,0,9))
   plot(as.dendrogram(hc),horiz=T)
   dev.off()
 

