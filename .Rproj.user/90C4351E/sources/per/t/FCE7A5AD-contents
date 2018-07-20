###################################################
### Code from Cury et al ###
############################
normalized prey abundance (Fig. 2A) and Change in variance across the range of normalized food
abundance ranging from -1.5 to 2 standard deviations in 8 classes.
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset
####### REQUIREMENTS: you need the script bootstrapping.R and to have the R packages 'gam',
#######  'changepoint', and 'time' installed
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Figs. 2A and B
####### MISCELLANEOUS:
####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
#############################################################
#############################################################
### LOAD the required packages
#############################################################
library(gam)
library(changepoint) # package version 0.3
library(time)
#############################################################
### !!!! SET the path to the data (download the global_dataset.csv file and paste it in the repository you
### want to work in) !!!
  #############################################################
setwd('D:/Dropbox/BlackBawks/PROJECTS/Farallon/OneThird/Data')
#############################################################
### LOAD the data
#############################################################
t_glob <- read.table('./global_dataset.csv', header=F, sep=",")
##########################
### Script to boostrap
bootstrap_size <- 1000
smoother_size <- 3
source('../bootstrapping.R')
##########################
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
###########################################
####### GLOBAL DATASET FIGURE (Fig 2A) ###
###########################################
if (png_figure==1){
  png('./Fig2A.png', width=1200, height=1200)
}
par(mar=c(5,5,3,1))
plot(min(prey),min(bird), col="transparent" ,main="",ylab="Breeding success (-)", xlab="Prey
     abundance (-)", cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, xlim=c(-2,3.5), ylim=c(-3.5, 3))
lines(titu, lwd=2)
xd<-c(fds[fds2,2],rev(fds[fds2,2]))
up.sd<-prediction$fit[fds2]+1.96*prediction$se.fit[fds2]
down.sd<-prediction$fit[fds2]-1.96*prediction$se.fit[fds2]
yd<-c(up.sd,rev(down.sd))
polygon(xd,yd,col="grey81",lty=1,border="transparent")
lines(titu$y~titu$x, lwd=4)
for (i in 1:length(data[,1])){
  points(data[i,2]~data[i,3],pch=16, col=paste(data[i,5]), cex=2)
}



### Identify the threshold
data2 <- cbind(titu$y[1:150],1,titu$x[1:150])
change_point <- cpt.reg(data2,penalty="SIC")
threshold <- titu$x[cpts(change_point)[1]]
abline(v=threshold, lwd=4, col="orange")
percentage_threshold <- (threshold+abs(min(prey)))/(max(prey)+abs(min(prey)))
text(threshold,-3.2, paste(signif(threshold,2), "\ni.e.",paste(signif(percentage_threshold*100,2)),"% of
                           max prey abundance" ), pos=4, cex=2.5)
text(2.5,2.5, paste("AIC=",signif(AIC(gam.object),5),sep=" "), pos=4, cex=2.5)
legend("bottomright", pch=16, col=unique(paste(t_glob[,5])), legend=unique(t_glob[,4]), bty="n",
       cex=2)
### plot the confidence interval
IC<-quantile(thresh, c(0.025, 0.975))
abline(v=IC[1], lwd=3, lty=2, col="black")
abline(v=IC[2], lwd=3, lty=2, col="black")
if (png_figure==1){
  dev.off()
}
###########################################
########## VARIANCE FIGURE (Fig 2B) ######
###########################################
min_prey <- min(prey)
max_prey <- max(prey)
diff_prey <- max_prey-min_prey
variance <- 0
count <- 0
for (i in seq(-2,3, by=0.5)){
  j <- i+0.5
  variance <- c(variance,var(bird[which(prey>i & prey <= j )], na.rm=T))
  count <- c(count,length(bird[which(prey>i & prey <= j )]))
}
variance <- variance[-1]
count <- count[-1]
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  png('./Fig2B.png', width=1200, height=1200)
}
par(mar=c(6,6,1,1))
plot(smooth(variance[2:9])~seq(-1.5,2,by=0.5), type="l", ylab="Variance of seabird breeding success",
     xlab="Prey abundance (-)", cex.axis=3, cex.lab=3, lwd=3)
abline(v=threshold, lwd=3, col="orange")
legend("topright", legend=c("Threshold value as identified\n in the global analysis"), col="orange",
       lwd=3, bty="n", cex=3)
if (png_figure==1){
  dev.off()
}
####################################################
############### Test for difference in variance ##
####################################################
### Bartlett test for variance
first_var <- bird[which(prey < -0.092)]
second_var <- bird[which(prey > -0.092)]
var.test(first_var,second_var)
####### PURPOSE: Script to fit the parametric model (see Table S2) to the global dataset for each
ecosystem
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Fig. 2C
####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
#############################################################
#############################################################
### !!!! SET the path to the data !!!!
#############################################################
setwd('XXX')
#############################################################
### LOAD the data
#############################################################
t_glob <- read.table('./global_dataset.csv', header=F, sep=";")
#########################
##### CONTROL PARAMETERS
#########################
### for nls to enable larger number of iterations and decrease the minFactor
nls_cont <- nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/(1024*20),
                        printEval = FALSE, warnOnly = FALSE)
###############################
###### PARAMETER ESTIMATION ###
###############################
### Initialization for the first ecosystem
toto <- unique (t_glob[,4])
j <- 1; i <- 1
data <- subset(t_glob, V4==toto[j])
### The global data set contains some holes in the times series. We first select the years we have both a
breeding success and an abundance estimate for the prey
bird<-as.numeric(data[,2])
prey<-as.numeric(data[,3])
req1<-which(is.na(bird)==T)
if (length(req1) >0){bird<-bird[-req1];prey<-prey[-req1] }
req2<-which(is.na(prey)==T)
if (length(req2) >0){bird<-bird[-req2];prey<-prey[-req2] }
model3 <- nls(bird~c+a*(1-exp(-b*prey)), start=list(a=1, b=1, c=0),control=nls_cont)
x_prey <- seq(min(prey),max(prey),by=(max(prey)-min(prey))/200)
y_bird <- coef(model3)[3]+coef(model3)[1]*(1-exp(-coef(model3)[2]*x_prey))
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  png('./Fig2C.png', width=1200, height=1200)
}
par(mar=c(5,5,3,1))
plot(bird~prey, xlim=c(-2,3.5), ylim=c(-2,2) ,ylab="Breeding success (-)", xlab="Prey abundance (-)",
     cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, col="transparent", main="", cex.main=2.5)
lines(smooth(y_bird)~x_prey,col=paste(data[i,5]), lwd=4)
for (j in 2: length(toto) ) {
  data <- subset(t_glob, V4==toto[j])
  bird<-as.numeric(data[,2])
  prey<-as.numeric(data[,3])
  req1<-which(is.na(bird)==T)
  if (length(req1) >0){bird<-bird[-req1];prey<-prey[-req1] }
  req2<-which(is.na(prey)==T)
  if (length(req2) >0){bird<-bird[-req2];prey<-prey[-req2] }
  model3 <- nls(bird~c+a*(1-exp(-b*prey)), start=list(a=1, b=1, c=0),control=nls_cont)
  x_prey <- seq(min(prey),max(prey),by=(max(prey)-min(prey))/200)
  y_bird <- coef(model3)[3]+coef(model3)[1]*(1-exp(-coef(model3)[2]*x_prey))
  lines(smooth(y_bird)~x_prey,col=paste(data[i,5]), lwd=4)
}
legend("bottomright", pch=16, col=unique(paste(t_glob[,5])), legend=unique(t_glob[,4]), bty="n",
       cex=2)
if (png_figure==1){
  dev.off()


  
  
  
  ####### PURPOSE: Script to fit the parametric model to the global dataset for each species
  ####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset
  ####### INPUTS: the global_dataset.csv file
  ####### OUTPUTS: Fig. 2D
  ####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
  #############################################################
  ############################################################
  ### !!!! SET the path to the data (download the global_dataset.csv file and paste it in the repository you
  want to work in) !!!
  #############################################################
setwd('XXX')
#############################################################
### LOAD the data
#############################################################
t_glob <- read.table('./global_dataset.csv', header=F, sep=";")
#########################
##### CONTROL PARAMETERS
#########################
### for nls to enable larger number of iterations and decrease the minFactor
nls_cont <- nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/(1024*20),
                        printEval = FALSE, warnOnly = FALSE)
#########################
##### Model FITS
#########################
colors <- rainbow(18)
data <- t_glob
un_data <- unique(paste(t_glob[,6])) # Get the tag for each data sea
### Initialization for the first species
i <- 1
subset_data <- subset(t_glob,subset=c(t_glob[,6]==paste(un_data[i])))
### The global data set contains some holes in the times series. We first select the years we have both a
breeding success and an abundance estimate for the prey
bird<-as.numeric(subset_data[,2])
prey<-as.numeric(subset_data[,3])
req1<-which(is.na(bird)==T)
if (length(req1) >0){bird<-bird[-req1];prey<-prey[-req1] }
req2<-which(is.na(prey)==T)
if (length(req2) >0){bird<-bird[-req2];prey<-prey[-req2] }
model3 <- nls(bird~c+a*(1-exp(-b*prey)), start=list(a=1, b=1, c=0),control=nls_cont)
x_prey <- seq(min(prey),max(prey),by=(max(prey)-min(prey))/200)
y_bird <- coef(model3)[3]+coef(model3)[1]*(1-exp(-coef(model3)[2]*x_prey))
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  png('./Fig2D.png', width=1200, height=1200)
}
par(mar=c(5,5,3,1))
plot(bird~prey, xlim=c(-2,3.5), ylim=c(-2,2) ,main="",ylab="Breeding success (-)", xlab="Prey
abundance (-)", cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, col="transparent")
lines(smooth(y_bird)~x_prey,col=colors[i], lwd=4)
legend_text <- paste(subset_data[1,6])
counter <- 1
for (i in 2:length(un_data)){
  if (i<6) {i2 <- i} else {i2 <- i-1}
  if (i!=6){
    subset_data <- subset(t_glob,subset=c(t_glob[,6]==paste(un_data[i])))  
    bird<-as.numeric(subset_data[,2])
    prey<-as.numeric(subset_data[,3])
    req1<-which(is.na(bird)==T)
    if (length(req1) >0){bird<-bird[-req1];prey<-prey[-req1] }
    req2<-which(is.na(prey)==T)
    if (length(req2) >0){bird<-bird[-req2];prey<-prey[-req2] }
    model3 <- nls(bird~c+a*(1-exp(-b*prey)), start=list(a=1, b=1, c=0),control=nls_cont)
    x_prey <- seq(min(prey),max(prey),by=(max(prey)-min(prey))/200)
    y_bird <- coef(model3)[3]+coef(model3)[1]*(1-exp(-coef(model3)[2]*x_prey))
    lines(smooth(y_bird)~x_prey,col=colors[i2], lwd=4)
    legend_text <- c(legend_text,paste(subset_data[1,6]))
    counter<- c(counter,i2)
  }
}
legend_text <- paste(seq(1,13), paste(un_data[-6]))
legend("bottomright",cex=2, col=colors[counter], legend=legend_text, bty="n", lwd=2)
if (png_figure==1){
  dev.off()
}




####### PURPOSE: Script to fit the parametric model to the global dataset for each ecosystem with
the data overplotted
####### INSTRUCTIONS: run this script to
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset, install the
package splines, quantreg and VGAM
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Fig. 3
####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
#############################################################
#############################################################
### LOAD the required packages
#############################################################
library(gam)
library(changepoint)
library(time)
#############################################################
### SET the path to the data
#############################################################
setwd('XXX')
#############################################################
### LOAD the data
#############################################################
t_glob <- read.table('./global_dataset.csv', header=F, sep=";")
#########################
##### CONTROL PARAMETERS
#########################
### for nls to enable larger number of iterations and decrease the minFactor
nls_cont <- nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/(1024*20),
                        printEval = FALSE, warnOnly = FALSE)
### Initialization for the first ecosystem
toto <- unique (t_glob[,4])
j=1; i=1
data <- subset(t_glob, V4==toto[j])
### The global data set contains some holes in the times series. We first select the years we have both a
breeding success and an abundance estimate for the prey
bird<-as.numeric(data[,2])
prey<-as.numeric(data[,3])
req1<-which(is.na(bird)==T)
if (length(req1) >0){bird<-bird[-req1];prey<-prey[-req1] }
req2<-which(is.na(prey)==T)
if (length(req2) >0){bird<-bird[-req2];prey<-prey[-req2] }
model3 <- nls(bird~c+a*(1-exp(-b*prey)), start=list(a=1, b=1, c=0),control=nls_cont)
x_prey <- seq(min(prey),max(prey),by=(max(prey)-min(prey))/200)
y_bird <- coef(model3)[3]+coef(model3)[1]*(1-exp(-coef(model3)[2]*x_prey))
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  png('./Fig3.png', width=1200, height=1200)
}
par(mar=c(5,5,3,1), mfrow=c(3,3))
plot(bird~prey, xlim=c(-2,3.5), ylim=c(-2,2) ,ylab="Breeding success (-)", xlab="Prey abundance (-)",
     cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, col="transparent", main=toto[j], cex.main=2.5)
lines(smooth(y_bird)~x_prey,col=paste(data[i,5]), lwd=4)
points(bird~prey, pch=17,col=paste(data[i,5]), lwd=4)
for (j in 2: length(toto) ) {
  data <- subset(t_glob, V4==toto[j])
  bird<-as.numeric(data[,2])
  prey<-as.numeric(data[,3])
  req1<-which(is.na(bird)==T)
  if (length(req1) >0){bird<-bird[-req1];prey<-prey[-req1] }
  req2<-which(is.na(prey)==T)
  if (length(req2) >0){bird<-bird[-req2];prey<-prey[-req2] }
  model3 <- nls(bird~c+a*(1-exp(-b*prey)), start=list(a=1, b=1, c=0),control=nls_cont)
  x_prey <- seq(min(prey),max(prey),by=(max(prey)-min(prey))/200)
  y_bird <- coef(model3)[3]+coef(model3)[1]*(1-exp(-coef(model3)[2]*x_prey))
  plot(bird~prey, xlim=c(-2,3.5), ylim=c(-2,2) ,ylab="Breeding success (-)", xlab="Prey
abundance (-)", cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, col="transparent", main=toto[j],
       cex.main=2.5)
  lines(smooth(y_bird)~x_prey,col=paste(data[i,5]), lwd=4)
  points(bird~prey, pch=17,col=paste(data[i,5]), lwd=4)
}
legend("bottomright", pch=16, col=unique(paste(t_glob[,5])), legend=unique(t_glob[,4]), bty="n",
       cex=2)
if (png_figure==1){
  dev.off()
}  
  
  





####### PURPOSE: Script to estimate the minimum of number of years required to estimate the
threshold
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Fig. S1
####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
####### WARNINGS: this is a bootstrap procedure. It will repeat 100 times this procedure (100
bootstraps...) so get some coffee.
#############################################################
#############################################################
### LOAD the required packages
#############################################################
library(gam)
library(changepoint)
library(time)
#############################################################
### !!!! SET the path to the data (download the global_dataset.csv file and paste it in the repository you
want to work in) !!!
  #############################################################
setwd('XXX')
#############################################################
### LOAD the data
#############################################################
t_glob <- read.table('./global_dataset.csv', header=F, sep=";")
### GAM
nb_boot <- 100
year_est <-rep(0,nb_boot)
threshold <- array(0,dim=c(nb_boot,38,nb_boot))
years <- 1964:2010
for (j in 1:nb_boot){
  prev <- progressBar()
  for (year_size in 5:42){
    prev2 <- progressBar()
    for (i in 1:nb_boot){
      prev2 <- progressBar(i/nb_boot, prev2)
      year_start <- sample(years[1:(length(years)-year_size)], size=1)
      data <- subset(t_glob,subset=c(t_glob[,1]>=year_start &
                                       t_glob[,1]<(year_start+year_size)))
      bird <- as.numeric(data[,2])
      prey <- as.numeric(data[,3])
      req1 <- which(is.na(bird)==T)
      if (length(req1) >0) {bird<-bird[-req1];prey<-prey[-req1] }
      req2<-which(is.na(prey)==T)
      if (length(req2) >0) {bird<-bird[-req2];prey<-prey[-req2] }
      gam.object<-gam(bird~s(prey,3), na.action='na.omit')
      prediction<-predict.gam(gam.object, se.fit=T)
      fds<-cbind(prediction$fit,prey,prediction$se.fit,0,0)
      fds2<-order(fds[,2])
      titu <- spline(prediction$fit[fds2]~prey[fds2],n=200)
      data2 <- cbind(titu$y[1:150],1,titu$x[1:150])
      change_point <- cpt.reg(data2,penalty="SIC")
      threshold[i,(year_size-4),j] <- titu$x[cpts(change_point)[1]]
    }
  }
}
year_est <-rep(0,nb_boot)
year_mean <-rep(0,38)
year_2_5 <-rep(0,38)
year_97_5 <-rep(0,38)
for (j in 1:nb_boot){
  xd<-c(1:38,rev(1:38))
  up.sd<-smooth(apply(threshold[,,j],2, quantile, probs=c(0.9)))
  down.sd<-smooth(apply(threshold[,,j],2, quantile, probs=c(0.1)))
  yd<-c(up.sd,rev(down.sd))
  year_est[j] <- max(which(up.sd>0.13),which(down.sd< -0.3))
}
for (j in 1:38){
  year_mean[j] <- mean(threshold[,j,])
  year_2_5[j] <- quantile(threshold[,j,],probs=0.0225)
  year_97_5[j] <- quantile(threshold[,j,],probs=0.975)
}
save(threshold,file=paste("./threshold_sim.Rdata", sep=""))
save(threshold,file=paste("./threshold_sim",j,".Rdata", sep=""))
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  pdf('./FigS1.pdf')
}
plot(year_mean, xaxt="n", xlab="Length of the time series (years)", ylab="Estimated threshold",
     ylim=c(-1.5,1), type="l", lwd=3)
axis(side=1, at=1:42, labels=5:46)
xd<-c(1:38,rev(1:38))
yd<-c(year_97_5,rev(year_2_5))
polygon(xd,yd,col="grey81",lty=1,border="transparent")
abline(h=-0.092, lwd=3, col="dark red", lty=1)
lines(year_mean,lwd=2)
abline(v=(mean(year_est)-4), lwd=2, col="orange")
abline(v=(quantile(year_est, probs=c(0.025,0.975))-4), lwd=2, lty=2, col="orange")
legend("bottomright", legend=c("Estimated threshold (global dataset)", "Estimated threshold (subset)",
                               "Minimum number of years required \n to estimate the threshold","Confidence interval for \n the
minimum number of years required"), cex=1, bty="n", lwd=2, col=c("dark red","black", "orange",
                                                                 "orange" ), lty=c(1,1,1,2))
if (png_figure==1){
  dev.off()
}


####### PURPOSE: Script to estimate when the maximum prey biomass is observed in the global
dataset
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Fig. S2
####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
#############################################################
#############################################################
### SET the path to the data (download the global_dataset.csv file and paste it in the repository you
want to work in)
#############################################################
setwd('XXX')
#############################################################
### LOAD the data
#############################################################
data <- read.table('./global_dataset.csv', header=F, sep=";")
year <- rep(0,19)
for (i in 1:19){
  data2 <- subset(data, subset=c(data[,8]==i))
  iwh <- which(data2[,3]==max(data2[,3], na.rm=T))
  year[i] <- data2[iwh,1]-data2[1,1]+1
}
names <- c("Krill_biomass", rep("Sardine_anchovy_w_agulhas",3),
           rep("Sandeel",4),rep("Shetland_sandeel_TSB",3),"Euphausiid_availability_index",rep("biom",2),rep("
all_rockfish",3),rep("herring_vpa_0",2))
max_bio <- as.data.frame(cbind(year,names))
index <- max_bio[1:(length(max_bio[,1])-1),2]!=max_bio[2:length(max_bio[,1]),2]
year2 <- year[index][year[index]>0]
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  png('./FigS2.png')
}
par(mar=c(5,5,1,4))
titi<-lowess(density(year2)$y~density(year2)$x, f=0.56)
plot(titi, type="l", xlab="Years after the beginning of data collection", ylab="Probability",
     cex.axis=2.5, cex.lab=1.7, xlim=c(0,50), lwd=2)
abline(v=round(titi$x[which(titi$y==max(titi$y))]), lwd=3)
text(5.,0.006, paste(round(titi$x[which(titi$y==max(titi$y))])," years", sep=""), pos=1, cex=2)
if (png_figure==1){
  dev.off()
}





####### PURPOSE: Script to run the generalized additive mixed model (GAMM) analysis on the data
set
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset, install the
package gamm4
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Fig. S3
####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
#############################################################
#############################################################
### LOAD the required packages
#############################################################
library(gam)
library(changepoint)
library(time)
library(gamm4)
#############################################################
### !!!! SET the path to the data (download the global_dataset.csv file and paste it in the repository you
want to work in) !!!
  #############################################################
setwd('XXX')
#############################################################
### LOAD the data
#############################################################
t_glob <- read.table('./global_dataset.csv', header=F, sep=";")
data <- t_glob
### The global data set contains some holes in the times series. We first select the years we have both a
breeding success and an abundance estimate for the prey
bird <- as.numeric(data[,2])
prey <- as.numeric(data[,3])
Eco <- as.numeric(data[,4])
req1 <- which(is.na(bird)==T)
if (length(req1) >0) {bird<-bird[-req1];prey<-prey[-req1]; Eco<-Eco[-req1] }
req2<-which(is.na(prey)==T)
if (length(req2) >0) {bird<-bird[-req2];prey<-prey[-req2] ; Eco<-Eco[-req2]}
### GAMM
gam.object <- gamm4(bird~s(prey), random = ~ (1|Eco), na.action="na.omit")
prediction <- predict.gam(gam.object$gam, se.fit=T)
fds <- cbind(prediction$fit,prey,prediction$se.fit,0,0)
fds2 <- order(fds[,2])
titu <- spline(prediction$fit[fds2]~prey[fds2],n=200)
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  png('./FigS3.png', width=1200, height=1200)
}
plot(min(prey),min(bird), col="transparent" ,main="",ylab="Breeding success (-)", xlab="Prey
     abundance (-)", cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, xlim=c(-2,3.5), ylim=c(-3.5, 3))
lines(titu, lwd=2)
xd<-c(fds[fds2,2],rev(fds[fds2,2]))
up.sd<-prediction$fit[fds2]+1.96*prediction$se.fit[fds2]
down.sd<-prediction$fit[fds2]-1.96*prediction$se.fit[fds2]
yd<-c(up.sd,rev(down.sd))
polygon(xd,yd,col="grey81",lty=1,border="transparent")
lines(titu$y~titu$x, lwd=4)
for (i in 1:length(data[,1])){
  points(data[i,2]~data[i,3],pch=16, col=paste(data[i,5]), cex=2)
}
## Identify the threshold
data2 <- cbind(titu$y[1:150],1,titu$x[1:150])
change_point <- cpt.reg(data2,penalty="SIC")
threshold <- titu$x[cpts(change_point)[1]]
abline(v=threshold, lwd=4, col="orange")
percentage_threshold <- (threshold+abs(min(prey)))/(max(prey)+abs(min(prey)))
text(threshold,-3.2, paste(signif(threshold,2), "\ni.e.",paste(signif(percentage_threshold*100,2)),"% of
max prey abundance" ), pos=4, cex=2.5)
legend("bottomright", pch=16, col=unique(paste(t_glob[,5])), legend=unique(t_glob[,4]), bty="n",
       cex=2)
if (png_figure==1){
  dev.off()
}




####### PURPOSE: Script to run the quantile regression analyses on the global dataset
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset, install the
package splines, quantreg and VGAM
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Fig. S4
####### REQUIREMENTS: you need the script bootstrapping.R
####### AUTHORS: Sylvain Bonhommeau & Philippe Cury (Last update 23rd of August 2011)
#############################################################
#############################################################
### LOAD the required packages
#############################################################
library(gam)
library(changepoint)
library(time)
library(splines)
library(quantreg)
library(VGAM)
#############################################################
### !!!! SET the path to the data (download the global_dataset.csv file and paste it in the repository you
want to work in) !!!
  #############################################################
setwd('XXXX')
t_glob <- read.table('./global_dataset.csv', header=F, sep=";")
data <- t_glob
### The global data set contains some holes in the times series. We first select the years we have both a
breeding success and an abundance estimate for the prey
bird <- as.numeric(data[,2])
prey <- as.numeric(data[,3])
req1 <- which(is.na(bird)==T)
if (length(req1) >0) {bird<-bird[-req1];prey<-prey[-req1] }
req2<-which(is.na(prey)==T)
if (length(req2) >0) {bird<-bird[-req2];prey<-prey[-req2] }
##########################
### Script to boostrap to plot as in Fig 2A
bootstrap_size <- 1000
smoother_size <- 3
source('./bootstrapping.R')
### Classic GAM as for the Figure 1A of the manuscript
gam.object<-gam(bird~s(prey,3), na.action='na.omit')
prediction<-predict.gam(gam.object, se.fit=T)
fds<-cbind(prediction$fit,prey,prediction$se.fit,0,0)
fds2<-order(fds[,2])
titu <- spline(prediction$fit[fds2]~prey[fds2],n=200)
### Fit a "quantile" gam on the data (with an offset of 10 since vgam requires that the response variable
is positive (then offset substracted for the plot)
fit <- vgam(bird+10 ~ s(prey), lms.bcn(zero=c(1,3)),trac=TRUE)
###quantile regression with GAM
qtreg_gam <- qtplot.lmscreg(fit, percentiles=c(10,90), main="Quantiles",pch=16, las=1,
                            ylab="Bredding success (-)", xlab="Prey abundance (-)", lwd=2, lcol=4, yaxt="n")
quant_matrix <- cbind(prey,qtreg_gam$fitted.values-10)
order_quant_matrix <- order(quant_matrix[,1])
titu95 <- spline(quant_matrix[order_quant_matrix,3]~quant_matrix[order_quant_matrix,1],n=200)
titu5 <- spline(quant_matrix[order_quant_matrix,2]~quant_matrix[order_quant_matrix,1],n=200)
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 0
if (png_figure==1){
  png('./FigS4.png', width=1200, height=1200)
  
}
par(mar=c(5,5,3,1))
plot(titu95, lwd=2, col="dark red", type="l",main="",ylab="Breeding success (-)", xlab="Prey
     abundance (-)", cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, xlim=c(-2,3.5), ylim=c(-3.5, 3) )
lines(titu, lwd=2)
lines(titu5, lwd=2, col="dark blue")
xd<-c(fds[fds2,2],rev(fds[fds2,2]))
up.sd<-prediction$fit[fds2]+1.96*prediction$se.fit[fds2]
down.sd<-prediction$fit[fds2]-1.96*prediction$se.fit[fds2]
yd<-c(up.sd,rev(down.sd))
polygon(xd,yd,col="grey81",lty=1,border="transparent")
for (i in 1:length(data[,1])){
  points(data[i,2]~data[i,3],pch=16, col=paste(data[i,5]), cex=2)
}
lines(titu, lwd=2)
lines(titu5, lwd=2, col="dark blue")
lines(titu95, lwd=2, col="dark red")
### Identify the threshold for each GAM
data2 <- cbind(titu$y[1:150],1,titu$x[1:150])
change_point <- cpt.reg(data2,penalty="SIC")
threshold <- titu$x[cpts(change_point)[1]]
abline(v=threshold, lwd=4, col="orange")
data5 <- cbind(titu5$y[1:150],1,titu5$x[1:150])
change_point5 <- cpt.reg(data5,penalty="SIC")
threshold5 <- titu5$x[cpts(change_point5)[1]]
data95 <- cbind(titu95$y[1:150],1,titu95$x[1:150])
change_point95 <- cpt.reg(data95,penalty="SIC")
threshold95 <- titu95$x[cpts(change_point95)[1]]
### plot the confidence interval
IC<-quantile(thresh, c(0.025, 0.975))
abline(v=IC[1], lwd=3, lty=2, col="black")
abline(v=IC[2], lwd=3, lty=2, col="black")
legend("bottomright", pch=16, col=unique(paste(t_glob[,5])), legend=unique(t_glob[,4]), bty="n",
       cex=2)
legend("topleft", legend=c("Quantile GAM 90%", "GAM as in Fig 1A","Quantile GAM 10%",
                           "Confidence interval") , lwd=4, lty=c(1,1,1,2,1) ,col=c("dark red", "black", "dark blue", "black",
                                                                                   "orange"), bty="n", cex=2)
if (png_figure==1){
  dev.off()
}