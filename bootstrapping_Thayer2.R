#Bootstrapping code
####### PURPOSE: Script to bootstrap on the data set to find the confidence interval for the threshold
#estimate (for Fig 2A)
####### INSTRUCTIONS: set the bootstrap sample size you want to do in the Fig2A_B.R script
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset
####### REQUIREMENTS: you need the script bootstrapping.R and to have the R packages 'gam',
#'changepoint', and 'time' installed
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Figs. 2A
#############################################################
### The global data set contains some holes in the times series. We first select the years we have both a
#breeding success and an abundance estimate for the prey

Run.Bootstrap <- function(t_glob){
  data <- t_glob[complete.cases(t_glob),]
  
  ids <- c(1:nrow(data))
  
  bird <- as.numeric(data[,2])  
  prey <- as.numeric(data[,3])
  req1 <- which(is.na(bird)==T)
  if (length(req1) >0) {bird<-bird[-req1];prey<-prey[-req1] }
  req2<-which(is.na(prey)==T)
  if (length(req2) >0) {bird<-bird[-req2];prey<-prey[-req2] }
  ####### RESAMPLING #########
  nb <- bootstrap_size
  thresh <- vector(length=nb) ## create a vector to be filled with the estimated threshold
  prev <- txtProgressBar(min=0,max=nb,style=3) ## just to display how long it remains for the bootstrap
  for (i in 1:nb) {
    setTxtProgressBar(prev,i)
    set.seed(i)
    bootsample <- sample(1:length(bird), replace=T)
    birdi <- bird[bootsample]
    preyi <- prey[bootsample]
    #idsi <- ids[bootsample]
    
    #dat2 <- data.frame(idsi,birdi,preyi)
    #dat2 <- dat2[order(dat2$idsi),]
    
    #######Computing threshold on new sample #######
    gam.object <- gam(birdi~s(preyi, k=smoother_size))#,data=dat2)
    prediction <- predict(gam.object, se.fit=T)
    fds <- cbind(prediction$fit,preyi,prediction$se.fit,0,0)
    fds2 <- order(fds[,2])
    titu <- spline(prediction$fit[fds2]~preyi[fds2],n=200)
    
    data <- cbind(titu$x[1:150],1,titu$y[1:150])
    
    x <- sum(gam.object$edf)
  
    if(round(x,4) < 2.1){
      thresh[i] <- NA
    }else{
      try({
        #change_point <- cpt.mean(data[,1],penalty="MBIC")
        change_point <- EnvCpt:::cpt.reg(data,penalty="SIC") ## Calculate the position of the threshold (see changepoint package)  ##cpt.reg
        cpts(change_point)
        #data[cpts(change_point),1]
        
      })
      #print(titu$x[cpts(change_point)[1]])
      thresh[i] <- titu$x[cpts(change_point)[1]]
    }
    
    
    remove(change_point)
  }
  close(prev)
  return(thresh)
}


### Reset the values to initial values


# data <- t_glob
# bird<-as.numeric(data[,2])
# prey<-as.numeric(data[,3])
# req1<-which(is.na(bird)==T)
# if (length(req1) >0){bird<-bird[-req1];prey<-prey[-req1] }
# req2<-which(is.na(prey)==T)
# if (length(req2) >0){bird<-bird[-req2];prey<-prey[-req2] }
# gam.object<-gam(bird~s(prey, smoother_size))
# prediction<-predict(gam.object, se.fit=T)
# fds<-cbind(prediction$fit,prey,prediction$se.fit,0,0)
# fds2<-order(fds[,2])
# titu <- spline(prediction$fit[fds2]~prey[fds2],n=200)
# 
# MaxPrey<-max(data[,3],na.rm=T)
# 
# MinPrey<-min(data[,3],na.rm=T)
# Quant<-quantile(thresh, c(0.1, 0.9),na.rm=T)
# PercMax<-round((max(Quant) + abs(MinPrey) )/(MaxPrey+abs(MinPrey)),2)*100
# PercMin<-round((min(Quant) + abs(MinPrey) )/(MaxPrey+abs(MinPrey)),2)*100
# 
# 
# print(paste("max value of threshold CV is ",as.character(max(Quant)),sep=""))
# print(paste("max value is ",PercMax,"% of max prey biomass (",as.character(MaxPrey),")",sep=""))
# print(paste("min value of threshold CV is ",as.character(min(Quant)),sep=""))
# print(paste("min value is ",PercMin,"% of max prey biomass (",as.character(MaxPrey),")",sep=""))


plotdata <- function(data,FigName,pngfigure,smoother_size,labs=FALSE,plotlines=TRUE){
  
  data <- data[complete.cases(data),]
  data$legend <- as.character(data$legend)
  predator <- data[,2]
  prey <- data[,3]
  
  gam.obj <- gam(predator~s(prey, k=smoother_size))
  
  
  if (pngfigure==1){
    pnam <- paste('../Figures/',FigName,sep='')
    png(pnam, width=1200, height=1200)
  }
  par(mar=c(7,6,3,1))
  
  MIN.x<-min(data[which(!is.na(data[,3])),3])
  MAX.x<-max(data[which(!is.na(data[,3])),3])
  MIN.y<-min(data[which(!is.na(data[,2])),2])
  MAX.y <- 3
  #MAX.y<-max(data[which(!is.na(data[,2])),2])
  
  prediction<-predict(gam.obj, se.fit=T)
  fds<-cbind(prediction$fit,prey,prediction$se.fit,0,0)
  fds2<-order(fds[,2])
  titu <- spline(prediction$fit[fds2]~prey[fds2],n=200)
  
  plot(min(prey),min(predator), col="transparent" ,main="",ylab="Breeding success (-)", 
       xlab="Prey abundance (-)", cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, 
       xlim=c(MIN.x,MAX.x),ylim=c(-2,MAX.y))
  lines(titu, lwd=2)
  xd<-c(fds[fds2,2],rev(fds[fds2,2]))
  up.sd<-prediction$fit[fds2]+1.96*prediction$se.fit[fds2]
  down.sd<-prediction$fit[fds2]-1.96*prediction$se.fit[fds2]
  yd<-c(up.sd,rev(down.sd))
  polygon(xd,yd,col="grey81",lty=1,border="transparent")
  lines(titu$y~titu$x, lwd=4)
  if(labs == TRUE){
    points(data[,2],data[,3],pch=16,col=data$colour,cex=1.5)  
  }else if(labs == FALSE){
    points(data[i,2]~data[i,3],pch=16, col='black', cex=1)
  }
  

  data2 <- cbind(titu$x[1:150],1,titu$y[1:150])
  
  tryCatch({
    #change_point <- cpt.mean(data2[,1],penalty="SIC0") 
    change_point <- EnvCpt:::cpt.reg(data2,penalty="SIC",pen.value=0) ## Calculate the position of the threshold (see changepoint package)  ##cpt.reg
  },error = function(e){
    #    ##cpt.reg
    change_point <- -10
  })
  
  #change_point <- cpt.mean(data2[,1],penalty="MBIC")     ##cpt.reg
  
  threshold <- titu$x[cpts(change_point)[1]]    ##
  if(plotlines == TRUE){
    abline(v=threshold, lwd=4, col="orange")  
  }
  

  percentage_threshold <- (threshold+abs(min(prey)))/(max(prey)+abs(min(prey)))
  
  if(labs == TRUE){
    if(plotlines == TRUE){
      text((threshold+0.3),(MAX.y-0.3), paste("Threshold: ",signif(threshold,2), "\ni.e.",
                                           paste(signif(percentage_threshold*100,2)),
                                           "% of max prey abundance"), pos=4, cex=2.5)
    }
    #text(2.5,2.5, paste("AIC=",signif(AIC(gam.obj),5),sep=" "), pos=4, cex=2.5)
    legend("bottomright", pch=16, col=unique(paste(data$colour)), legend=unique(data$legend),
          cex=2)
  }
 
  ### plot the confidence interval
  IC<-quantile(thresh, c(0.025, 0.975),na.rm=T)
  
  if(plotlines == TRUE){
    abline(v=IC[1], lwd=3, lty=2, col="black")
    abline(v=IC[2], lwd=3, lty=2, col="black")  
  }
  
  if (pngfigure==1){
    dev.off()
  }
  
  return(change_point)
}





######################################################################################
#### Calculate parametric AIC values for global models ####
###########################################################



# 
# nls_cont <- nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/(1024*20),
#                         printEval = FALSE, warnOnly = FALSE)




















