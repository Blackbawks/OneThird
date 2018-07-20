######################################################################

#### Pull data together for modelling ####

library(dplyr)
library(foreach)
library(DescTools)
library(changepoint)
source('../bootstrapping_Thayer.R')

wd <- "D:/Dropbox/BlackBawks/PROJECTS/Farallon/OneThird/Data/"

setwd(wd)

Globaldat <- read.table('./Raw_data.csv', header=T, sep=",")


preylist <- c("anchovy.biomass","rockfish.jrssurvey",#"sardine.ichthyo",
              "sardine.biomass","squid.jrssurvey")#,"squid.ichthyo"


seabirdlist <- c("BRPE.SBI","BRPE.ANA","LETE.VB","BRAC.SFI","COMU.SFI",
                "PIGU.SFI","PECO.SFI","RHAU.SFI","RHAU.ANI",
                "BRAC.ANI","BRAC.ALZ","PECO.ANI")


sealionlist <- c("CASL_SBI","CASL.SCI","CASL.SMI","CASL.SNI")


salmon <- "SALM.SAC"


getData <- function(predator,prey){
  A1 <- data.frame(year = Globaldat$Year, 
                   prd = Globaldat[predator], 
                   pry = Globaldat[prey],
                   prdname=rep(predator,42),
                   pryname=rep(prey,42))
  return(A1)
}


seabirds <- foreach(i = preylist, .combine='rbind') %do% {

  X <- foreach(j=seabirdlist,.combine='rbind') %do% {
    A <- getData(j,i)
    
    names(A) <- c('year','prd','pry','prdname','pryname')
    return(A)
  }
  
  return(X)
  
}



sealions <- foreach(i = preylist, .combine='rbind') %do% {
  
  X <- foreach(j=sealionlist,.combine='rbind') %do% {
    A <- getData(j,i)
    names(A) <- c('year','prd','pry','prdname','pryname')
    return(A)
  }
  
  return(X)
  
}


salmon <- foreach(i = preylist, .combine='rbind') %do% {
  
  X <- foreach(j=salmon,.combine='rbind') %do% {
    A <- getData(j,i)
    names(A) <- c('year','prd','pry','prdname','pryname')
    return(A)
  }
  
  return(X)
  
}

############################################################################################################

#### First we run the GAMS with mgcv in order to get the EDF values out

# Abundance of I) anchovy prey relative to productivity of 
# a) Brandt’s cormorants at SE Farallon Island (SFI), 
# b) Año Nuevo Island (ANI) and 
# c) Alcatraz Island (ALZ), 
# d) common murre at SFI, 
# e) rhinoceros auklet at SFI and 
# f) ANI, 
# g) brown pelican at Santa Barbara Island (SBI) and 
# h) Anacapa Island (ANA), 
# i) least tern at Venice Beach (VB), and 
# j) California sea lions at San Miguel Island (SMI) 
# k) San Clemente Island (SCI), 
# l) SBI, and 
# m) San Nicolas Island (SNI); 
# 
# II) juvenile rockfish prey relative to productivity of 
# a) common murre at SFI, 
# b) rhinoceros auklet at SFI and 
# c) ANI, 
# d) pelagic cormorant at SFI and 
# e) ANI, 
# f) pigeon guillemot at SFI, and 
# g) Chinook salmon from the Central Valley fall run (CV). 
# 
# III) sardine prey relative to 
# a) productivity of sea lions at (SMI) and 
# b) Chinook salmon from CV; and 
# 
# IV) market squid prey relative to productivity of 
# a) common murre at SFI, 
# b) rhinoceros auklet at SFI and 
# c) ANI, and 
# d) sea lions at SMI.



detach("package:mgcv", unload=TRUE)
library(gam)
smoother_size=2


Ia <- getData('BRAC.SFI','anchovy.biomass')
Ib <- getData('BRAC.ANI','anchovy.biomass')
Ic <- getData('BRAC.ALZ','anchovy.biomass')
Id <- getData('COMU.SFI','anchovy.biomass')
Ie <- getData('RHAU.SFI','anchovy.biomass')
If <- getData('RHAU.ANI','anchovy.biomass')
Ig <- getData('BRPE.SBI','anchovy.biomass')
Ih <- getData('BRPE.ANA','anchovy.biomass')
Ii <- getData('LETE.VB','anchovy.biomass')
Ij <- getData('CASL.SMI','anchovy.biomass')
Ik <- getData('CASL.SCI','anchovy.biomass')
Il <- getData('CASL_SBI','anchovy.biomass')
Im <- getData('CASL.SNI','anchovy.biomass')

IIa <- getData('COMU.SFI','rockfish.jrssurvey')
IIb <- getData('RHAU.SFI','rockfish.jrssurvey')
IIc <- getData('RHAU.ANI','rockfish.jrssurvey')
IId <- getData('PECO.SFI','rockfish.jrssurvey')
IIe <- getData('PECO.ANI','rockfish.jrssurvey')
IIf <- getData('PIGU.SFI','rockfish.jrssurvey')
IIg <- getData('SALM.SAC','rockfish.jrssurvey')


IIIa <- getData('CASL.SMI','sardine.ichthyo')
IIIb <- getData('SALM.SAC','sardine.ichthyo')
IIIc <- getData('BRPE.SBI','sardine.ichthyo')
IIId <- getData('BRPE.ANA','sardine.ichthyo')


IVa <- getData('COMU.SFI','squid.ichthyo')
IVb <- getData('RHAU.SFI','squid.ichthyo')
IVc <- getData('RHAU.ANI','squid.ichthyo')
IVd <- getData('CASL.SMI','squid.ichthyo')



DATAS <- list(Ia,Ib,Ic,Id,Ie,If,Ig,Ih,Ii,Ij,Ik,Il,Im,IIa,IIb,IIc,IId,IIe,IIf,IIg,IIIa,IIIb,IVa,IVb,IVc,IVd)

#DATAS <- list(Ia,Ib,Ic,Id,Ie,If,Ig,Ih,Ii,Ij,Ik,Il,Im,IIa,IIb,IIc,IId,IIe,IIf,IIg,IIIa,IIIb,IVa,IVb,IVc,IVd,IIIa2,IIIb2,IVa2,IVb2,IVc2,IVd2)


DATAS <- list(IIIc,IIId,IIIc2,IIId2)

EDFval <- function(data){
  data <- data[complete.cases(data),]
  predator <- data[,2]
  prey <- data[,3]
  
  gam.obj <- gam(predator~s(prey, k=smoother_size))
  EDF <- sum(gam.obj$edf)
  
  pd <- names(data)[2]
  py <- names(data)[3]
  
  M <- data.frame(pd,py,as.character(EDF))
  names(M) <- c('pred','prey','edf')
  return(M)
}


OUT <- foreach(x = DATAS, .combine='rbind') %do% {
  return(EDFval(x))
}

#write.csv(OUT,'EDFs.csv',row.names=F)



bootstrap_size <- 1000
smoother_size <- 3



DATAS <- list(Ia,Ib,Ic,Id,Ie,If,Ig,Ih,Ii,Ij,Ik,Il,Im,IIa,IIb,IIc,IId,IIe,IIf,IIg,IIIa,IIIb,IVa,IVb,IVc,IVd)

FIGNAMES <- list("Ia",'Ib','Ic','Id','Ie','If','Ig','Ih','Ii','Ij','Ik','Il','Im',
                 'IIa','IIb','IIc','IId','IIe','IIf','IIg',
                 'IIIa','IIIb',
                 'IVa','IVb','IVc','IVd')


for(i in 1:length(DATAS)){
  t_glob <- DATAS[[i]]
  
  thresh <- Run.Bootstrap(t_glob)
  
  pNAME <- paste(FIGNAMES[i],'.png',sep='')
  
  
  data <- t_glob[complete.cases(t_glob),]
  predator <- data[,2]
  prey <- data[,3]
  
  gam.obj <- gam(predator~s(prey, k=smoother_size))
  EDF <- sum(gam.obj$edf)
  if(round(EDF,4) > 2){
    print(pNAME)
    plotdata(t_glob,pNAME,pngfigure=1,3)  
  }
  
  
    
}


t_glob <- seabirds
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'anchovy.biomass')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'red'
t_glob$colour[which(t_glob$pryname == 'sardine.biomass')] <- 'dark blue'
t_glob$colour[which(t_glob$pryname == 'squid.jrssurvey')] <- 'dark green'

t_glob$legend[which(t_glob$pryname == 'anchovy.biomass')] <- 'seabird-anchovy'
t_glob$legend[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'seabird-rockfish'
t_glob$legend[which(t_glob$pryname == 'sardine.biomass')] <- 'seabird-sardine'
t_glob$legend[which(t_glob$pryname == 'squid.jrssurvey')] <- 'seabird-squid'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'seabirds.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)





t_glob <- sealions
## remove rockfish...
t_glob <- tbl_df(sealions) %>% filter(pryname != 'rockfish.jrssurvey')
t_glob <- data.frame(t_glob)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'anchovy.biomass')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'sardine.biomass')] <- 'red'
t_glob$colour[which(t_glob$pryname == 'squid.jrssurvey')] <- 'dark green'

t_glob$legend[which(t_glob$pryname == 'anchovy.biomass')] <- 'sealion-anchovy'
t_glob$legend[which(t_glob$pryname == 'sardine.biomass')] <- 'sealion-sardine'
t_glob$legend[which(t_glob$pryname == 'squid.jrssurvey')] <- 'sealion-squid'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'sealions.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=2,labs=TRUE)





t_glob <- tbl_df(salmon) %>% filter(!pryname %in% c('squid.jrssurvey','anchovy.biomass'))
t_glob <- data.frame(t_glob)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'sardine.biomass')] <- 'red'

t_glob$legend[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'salmon-rockfish'
t_glob$legend[which(t_glob$pryname == 'sardine.biomass')] <- 'salmon-sardine'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'salmon.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)





















