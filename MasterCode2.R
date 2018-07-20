######################################################################

#### Pull data together for modelling ####

library(dplyr)
library(foreach)
library(DescTools)
library(changepoint)
library(mgcv)
wd <- "D:/Dropbox/BlackBawks/PROJECTS/Farallon/OneThird/Data/"

setwd(wd)

source('../bootstrapping_Thayer.R')
Globaldat <- read.table('./Raw_data.csv', header=T, sep=",")

'%!in%' <- function(x,y)!('%in%'(x,y))

getData <- function(predator,prey){
  A1 <- data.frame(year = Globaldat$Year, 
                   prd = Globaldat[predator], 
                   pry = Globaldat[prey],
                   prdname=rep(predator,42),
                   pryname=rep(prey,42))
  names(A1) <- c('year','prd','pry','prdname','pryname')
  return(A1)
}

bootstrap_size <- 1000
smoother_size <- 3
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
IIIe <- getData('CASL.SCI','sardine.ichthyo')
IIIf <- getData('CASL_SBI','sardine.ichthyo')
IIIg <- getData('CASL.SNI','sardine.ichthyo')


IVa <- getData('COMU.SFI','squid.ichthyo')
IVb <- getData('RHAU.SFI','squid.ichthyo')
IVc <- getData('RHAU.ANI','squid.ichthyo')
IVd <- getData('CASL.SMI','squid.ichthyo')
IVe <- getData('CASL.SCI','squid.ichthyo')
IVf <- getData('CASL_SBI','squid.ichthyo')
IVg <- getData('CASL.SNI','squid.ichthyo')


##### Seabirds, Sealions, Salmon (by prey)
## seabirds

sbrds <- list(Ia,Ib,Ic,Id,Ie,If,Ig,Ih,Ii,IIa,IIb,IIc,IId,IIe,IIf,IIIc,IIId,IVa,IVb,IVc)
seabirds <- do.call('rbind',sbrds)

t_glob <- seabirds
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'anchovy.biomass')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'red'
t_glob$colour[which(t_glob$pryname == 'sardine.ichthyo')] <- 'blue'
t_glob$colour[which(t_glob$pryname == 'squid.ichthyo')] <- 'green'

t_glob$legend[which(t_glob$pryname == 'anchovy.biomass')] <- 'anchovy'
t_glob$legend[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'rockfish'
t_glob$legend[which(t_glob$pryname == 'sardine.ichthyo')] <- 'sardine'
t_glob$legend[which(t_glob$pryname == 'squid.ichthyo')] <- 'squid'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'seabirds_test.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)



## sealions

selins <- list(Ij,Ik,Il,Im,IIIa,IIIe,IIIf,IIIg,IVd,IVe,IVf,IVg)
sealions <- do.call('rbind',selins)
t_glob <- sealions
## remove rockfish...
t_glob <- tbl_df(sealions) %>% filter(pryname != 'rockfish.jrssurvey')
t_glob <- data.frame(t_glob)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'anchovy.biomass')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'sardine.ichthyo')] <- 'blue'
t_glob$colour[which(t_glob$pryname == 'squid.ichthyo')] <- 'green'

t_glob$legend[which(t_glob$pryname == 'anchovy.biomass')] <- 'anchovy'
t_glob$legend[which(t_glob$pryname == 'sardine.ichthyo')] <- 'sardine'
t_glob$legend[which(t_glob$pryname == 'squid.ichthyo')] <- 'squid'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'sealions.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)

### Salmon

salms <- list(IIg,IIIb)
salmon <- do.call('rbind',salms)
t_glob <- tbl_df(salmon) %>% filter(!pryname %in% c('squid.jrssurvey','anchovy.biomass'))
t_glob <- data.frame(t_glob)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'red'
t_glob$colour[which(t_glob$pryname == 'sardine.ichthyo')] <- 'blue'

t_glob$legend[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'rockfish'
t_glob$legend[which(t_glob$pryname == 'sardine.ichthyo')] <- 'sardine'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'salmon.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE,plotlines=FALSE)


################################################################################################################

### Anchovy
acns <- list(Ia,Ib,Ic,Id,Ie,If,Ig,Ih,Ii,Ij,Ik,Il,Im)
anchovy <- do.call('rbind',acns)

t_glob <- data.frame(anchovy)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$prdname %in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'red'
t_glob$colour[which(t_glob$prdname %!in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'black'
  
t_glob$legend[which(t_glob$prdname %in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'sealions'
t_glob$legend[which(t_glob$prdname %!in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'seabirds'
  

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'anchovy.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)






### Rockfish
rcfs <- list(IIa,IIb,IIc,IId,IIe,IIf,IIg)
rockfish <- do.call('rbind',rcfs)

t_glob <- data.frame(rockfish)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$prdname %in% c('SALM.SAC'))] <- 'green'
t_glob$colour[which(t_glob$prdname %!in% c('SALM.SAC'))] <- 'black'

t_glob$legend[which(t_glob$prdname %in% c('SALM.SAC'))] <- 'salmon'
t_glob$legend[which(t_glob$prdname %!in% c('SALM.SAC'))] <- 'seabirds'


thresh <- Run.Bootstrap(t_glob)
pNAME <- 'rockfish.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)





### Sardine
srds <- list(IIIa,IIIb,IIIc,IIId,IIIe,IIIf,IIIg)
sardine <- do.call('rbind',srds)

t_glob <- data.frame(sardine)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$prdname %in% c('SALM.SAC'))] <- 'green'
t_glob$colour[which(t_glob$prdname %in% c('BRPE.SBI','BRPE.ANA'))] <- 'black'
t_glob$colour[which(t_glob$prdname %in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'red'

t_glob$legend[which(t_glob$prdname %in% c('SALM.SAC'))] <- 'salmon'
t_glob$legend[which(t_glob$prdname %in% c('BRPE.SBI','BRPE.ANA'))] <- 'seabirds'
t_glob$legend[which(t_glob$prdname %in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'sealions'


thresh <- Run.Bootstrap(t_glob)
pNAME <- 'sardine.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE,plotlines=FALSE)




### Squid
sqds <- list(IVa,IVb,IVc,IVd,IVe,IVf,IVg)
squid <- do.call('rbind',sqds)
t_glob <- data.frame(squid)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$prdname %in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'red'
t_glob$colour[which(t_glob$prdname %!in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'black'

t_glob$legend[which(t_glob$prdname %in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'sealions'
t_glob$legend[which(t_glob$prdname %!in% c('CASL.SMI','CASL.SCI','CASL_SBI','CASL.SNI'))] <- 'seabirds'


thresh <- Run.Bootstrap(t_glob)
pNAME <- 'squid.png'

out <- plotdata(t_glob,pNAME,pngfigure=0,smoother_size=3,labs=TRUE)






#######################################################################################################################
##### CENCAL



cencal <- list(Ia,Ib,Ic,Id,Ie,If,Ig,IIa,IIb,IIc,IId,IIe,IIf,IIg,IIIb,IVa,IVb,IVc)
cen.cal <- do.call('rbind',cencal)

t_glob <- data.frame(cen.cal)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'anchovy.biomass')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'red'
t_glob$colour[which(t_glob$pryname == 'sardine.ichthyo')] <- 'blue'
t_glob$colour[which(t_glob$pryname == 'squid.ichthyo')] <- 'green'

t_glob$legend[which(t_glob$pryname == 'anchovy.biomass')] <- 'anchovy'
t_glob$legend[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'rockfish'
t_glob$legend[which(t_glob$pryname == 'sardine.ichthyo')] <- 'sardine'
t_glob$legend[which(t_glob$pryname == 'squid.ichthyo')] <- 'squid'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'centralcal.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)


#################




socal <- list(Ih,Ii,Ij,Ik,Il,Im,IIIa,IIIc,IIId,IIIe,IIIf,IIIg,IVd,IVe,IVf,IVg)
so.cal <- do.call('rbind',socal)
t_glob <- data.frame(so.cal)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'anchovy.biomass')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'red'
t_glob$colour[which(t_glob$pryname == 'sardine.ichthyo')] <- 'blue'
t_glob$colour[which(t_glob$pryname == 'squid.ichthyo')] <- 'green'

t_glob$legend[which(t_glob$pryname == 'anchovy.biomass')] <- 'anchovy'
t_glob$legend[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'rockfish'
t_glob$legend[which(t_glob$pryname == 'sardine.ichthyo')] <- 'sardine'
t_glob$legend[which(t_glob$pryname == 'squid.ichthyo')] <- 'squid'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'southerncal.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)


###########################################################################################################

#### Global ####


DATAS <- list(Ia,Ib,Ic,Id,Ie,If,Ig,Ih,Ii,Ij,Ik,Il,Im,IIa,IIb,IIc,IId,IIe,IIf,IIg,IIIa,IIIb,
              IIIc,IIId,IIIe,IIIf,IIIg,IVa,IVb,IVc,IVd,IVe,IVf,IVg)

global <- do.call('rbind',DATAS)


t_glob <- data.frame(global)
t_glob$colour <- NA
t_glob$legend <- NA

t_glob$colour[which(t_glob$pryname == 'anchovy.biomass')] <- 'black'
t_glob$colour[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'red'
t_glob$colour[which(t_glob$pryname == 'sardine.ichthyo')] <- 'blue'
t_glob$colour[which(t_glob$pryname == 'squid.ichthyo')] <- 'green'

t_glob$legend[which(t_glob$pryname == 'anchovy.biomass')] <- 'anchovy'
t_glob$legend[which(t_glob$pryname == 'rockfish.jrssurvey')] <- 'rockfish'
t_glob$legend[which(t_glob$pryname == 'sardine.ichthyo')] <- 'sardine'
t_glob$legend[which(t_glob$pryname == 'squid.ichthyo')] <- 'squid'

thresh <- Run.Bootstrap(t_glob)
pNAME <- 'global.png'

out <- plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=TRUE)



#################################################################################################################################
#################################################################################################################################
#################################################################################################################################





############ This function will pull out the EDF value from the gam

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





############# This plots all the individual combinations 
### If the EDF value is a straight line (i.e. 2), then we don't plot that line

DATAS <- list(Ia,Ib,Ic,Id,Ie,If,Ig,Ih,Ii,Ij,Ik,Il,Im,IIa,IIb,IIc,IId,IIe,IIf,IIg,IIIa,IIIb,
              IIIc,IIId,IIIe,IIIf,IIIg,IVa,IVb,IVc,IVd)

FIGNAMES <- list("Ia",'Ib','Ic','Id','Ie','If','Ig','Ih','Ii','Ij','Ik','Il','Im',
                 'IIa','IIb','IIc','IId','IIe','IIf','IIg',
                 'IIIa','IIIb','IIIc','IIId','IIIe','IIIf','IIIg',
                 'IVa','IVb','IVc','IVd')


DATAS <- list(IVf,IVg)
FIGNAMES <-list('IVf','IVg')

for(i in 1:length(DATAS)){
  t_glob <- DATAS[[i]]
  t_glob$colour <- 'black'
  thresh <- Run.Bootstrap(t_glob)
  
  pNAME <- paste(FIGNAMES[i],'.png',sep='')
  
  
  data <- t_glob[complete.cases(t_glob),]
  predator <- data[,2]
  prey <- data[,3]
  
  gam.obj <- gam(predator~s(prey, k=smoother_size))
  EDF <- sum(gam.obj$edf)
  if(round(EDF,4) > 2){
    print(pNAME)
    plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=FALSE,pntsize=2)  
  }else{
    print(pNAME)
    plotdata(t_glob,pNAME,pngfigure=1,smoother_size=3,labs=FALSE,pntsize=2,plotlines=FALSE)  
  }
  
  
    
}


#### Grab EDF values for the 'global' data

EDFval(anchovy)
EDFval(sealions)
EDFval(so.cal)
EDFval(squid)
EDFval(global)

EDFval(sardine)
EDFval(seabirds)
EDFval(cen.cal)
EDFval(rockfish)
EDFval(salmon)


#################################################################################################################








