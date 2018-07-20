###Normalized prey abundance (Fig. 2A) and Change in variance across the range of normalized food
###abundance ranging from -1.5 to 2 standard deviations in 8 classes.
####### INSTRUCTIONS: set the path where you downloaded and saved the global dataset
####### REQUIREMENTS: you need the script bootstrapping.R and to have the R packages 'gam',
#'changepoint'
####### INPUTS: the global_dataset.csv file
####### OUTPUTS: Fig. 2
#############################################################
### LOAD the required packages
#############################################################
library(gam)
library(changepoint) 
#library(time)
#############################################################
### !!!! SET the path to the data (download the global_dataset.csv file and paste it in the repository you
#want to work in) !!!
  #############################################################
setwd("D:/Dropbox/BlackBawks/PROJECTS/Farallon/OneThird/Data/") 
#############################################################
### LOAD the data
#############################################################
t_glob <- read.table('./Socal.thresholds.sealion.n.d.anch-biomass.n.d.1975.csv', header=T, sep=",")
#### data <- read.csv("data.seabird.prey.csv") # THIS IS THE COMMAND I USE TO OPEN MY  NEW DATASET.  DOES THE COMMAND ABOVE FROM THE SAMPLE SOMEHOW TURN THE LIST OF DATA IN THE CSV FILE INTO A TABLE?  HOW MANY COLUMNS, ROWS, OR DOES IT MATTER?
##########################
### Script to boostrap
bootstrap_size <- 1000
smoother_size <- 2
source('../bootstrapping_Thayer.R')
##########################
### Figure (set png_figure to 1) or not (any other value)
png_figure <- 1
###########################################
####### GLOBAL DATASET FIGURE (Fig 2A) ###
###########################################

################# Change './Fig2A.png' to whatever you want.. you could make it
### './Socal_treshold.png' or whatever else... 

if (png_figure==1){
  png('./Fig2A.png', width=1200, height=1200)
}
par(mar=c(7,6,3,1))

##################################
### GET LIMITS FOR THE FIGURE ####
##################################

MIN.x<-min(data[which(!is.na(data[,3])),3])
MAX.x<-max(data[which(!is.na(data[,3])),3])
MIN.y<-min(data[which(!is.na(data[,2])),2])
MAX.y<-max(data[which(!is.na(data[,2])),2])


plot(min(prey),min(bird), col="transparent" ,main="",ylab="Breeding success (-)", 
     xlab="Prey abundance (-)", cex.lab=2.5, cex.axis=2.5, pch=16, cex.main=2, 
     xlim=c(MIN.x,MAX.x),ylim=c(MIN.y,MAX.y))
lines(titu, lwd=2)
xd<-c(fds[fds2,2],rev(fds[fds2,2]))
up.sd<-prediction$fit[fds2]+1.96*prediction$se.fit[fds2]
down.sd<-prediction$fit[fds2]-1.96*prediction$se.fit[fds2]
yd<-c(up.sd,rev(down.sd))
polygon(xd,yd,col="grey81",lty=1,border="transparent")
lines(titu$y~titu$x, lwd=4)
for (i in 1:length(data[,1])){
  points(data[i,2]~data[i,3],pch=16, col='black', cex=2) #paste(data[i,5])
}







### Identify the threshold
data2 <- cbind(titu$y[1:150],1,titu$x[1:150])
change_point <- EnvCpt:::cpt.reg(data2,penalty="SIC") ## Calculate the position of the threshold (see changepoint package)  ##cpt.reg
change_point <- cpt.mean(data2[,1],penalty="SIC")     ##cpt.reg
threshold <- titu$x[cpts(change_point)[1]]    ##



abline(v=threshold, lwd=4, col="orange")
percentage_threshold <- (threshold+abs(min(prey)))/(max(prey)+abs(min(prey)))
text((threshold+0.3),(MAX.y-0.3), paste("Threshold: ",signif(threshold,2), "\ni.e.",
                                        paste(signif(percentage_threshold*100,2)),
                                              "% of max prey abundance"), pos=4, cex=2.5)
text(2.5,2.5, paste("AIC=",signif(AIC(gam.object),5),sep=" "), pos=4, cex=2.5)
legend("bottomright", pch=16, col=unique(paste(t_glob[,5])), legend=unique(t_glob[,4]),
       cex=2)
### plot the confidence interval
IC<-quantile(thresh, c(0.1, 0.9),na.rm=T)
abline(v=IC[1], lwd=3, lty=2, col="black")
abline(v=IC[2], lwd=3, lty=2, col="black")
if (png_figure==1){
  dev.off()
}


