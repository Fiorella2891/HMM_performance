rm(list=ls())
#Load package
library (ctmm)
library(dplyr)
#Load example buffalo data

# KERNEL ESTIMATION
################################################################################
especie <- "SULA"
#POPULATION LEVEL
datos1 = read.csv(paste0("outputs/",especie,"_4.csv"))
d1 <- datos1 %>% 
  dplyr::select("ave", "ID","x","y","datetime","dives")
colnames(d1) <- c("ave","ID","longitude","latitude","timestamp","dives")
d1$dives[which(d1$dives > 0)] <- "dives" 

d2 <- datos1 %>%
  dplyr::select("ave","ID","x","y","datetime","model")
colnames(d2) <- c("ave","ID","longitude","latitude","timestamp","dives")
d2$dives[which(d2$dives == 2)] <- "foraging_HMM"

dTDR <- d1[which(d1$dives == "dives"),]
dHMM <- d2[which(d2$dives == "foraging_HMM"),]

# DIVES
dTDR$timestamp <- as.POSIXct(strptime(dTDR$timestamp,format="%Y-%m-%d %H:%M"),
                                tz="GMT")
viajes = unique(dTDR$ID)

ALL1 <- dTDR$timestamp[dTDR$ID == viajes[1]]
for (i in 2:length(viajes)){
  DATA <- dTDR$timestamp[dTDR$ID == viajes[i]]
  DATA <- DATA - DATA[1] # start times at zero
  DATA <- DATA + tail(ALL1,1) + (1 %#% 'week')
  print(range(DATA))# start next trip one week after the previous
  ALL1 <- c(ALL1,DATA)
} 

dTDR$newdate <- ALL1
solo_TDR <- dTDR
TDR <- dTDR[,c(3,4,6,7)]
colnames(TDR)[4] <- "timestamp"
t_TDR <- as.telemetry(TDR, projection = "+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

variograma <- variogram(t_TDR,res = 2)

#KDE and AKDE
M.IID_T <- ctmm.fit(t_TDR) # no autocorrelation timescales
m.ouf_T <- ctmm.guess(t_TDR,interactive=FALSE) # automated model guess
M.OUF_T <- ctmm.fit(t_TDR,m.ouf_T)

UD0_TDR <- akde(t_TDR,M.IID_T)
save(UD0_TDR, file = paste0("outputs/",especie,"_DIVE_KDE.Rdata"))
UD2w_TDR <- akde(t_TDR,M.OUF_T, weights=TRUE, trace = TRUE)
save(UD2w_TDR, file = paste0("outputs/",especie,"_DIVE_AKDE.Rdata"))

# HMM
dHMM$timestamp <- as.POSIXct(strptime(dHMM$timestamp,format="%Y-%m-%d %H:%M"),
                             tz="GMT")
viajes = unique(dHMM$ID)
ALL <- dHMM$timestamp[dHMM$ID == viajes[1]]
for (i in 2:length(viajes)){
  DATA <- dHMM$timestamp[dHMM$ID == viajes[i]]
  DATA <- DATA - DATA[1] # start times at zero
  DATA <- DATA + tail(ALL,1) + (1 %#% 'week') # start next trip one week after the previous
  ALL <- c(ALL,DATA)
} 

dHMM$newdate <- ALL
solo_HMM <- dHMM
HMM <- dHMM[,c(3,4,6,7)]
colnames(HMM)[4] <- "timestamp"
t_HMM <- as.telemetry(HMM, projection = "+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

variograma <- variogram(t_HMM)
plot(variograma,level=c(0.5,0.95))

# KDE and AKDE
M.IID_H <- ctmm.fit(t_HMM) # no autocorrelation timescales
m.ouf_H <- ctmm.guess(t_HMM,interactive=FALSE) # automated model guess
M.OUF_H <- ctmm.fit(t_HMM,m.ouf_H)

UD0_HMM <- akde(t_HMM,M.IID_H)
save(UD0_HMM, file = paste0("outputs/",especie,"_HMM_KDE.Rdata"))
UD2w_HMM <- akde(t_HMM, M.OUF_H, weights=TRUE, trace = TRUE)
save(UD2w_HMM, file = paste0("outputs/",especie,"_HMM_AKDE.Rdata"))

# Overlap
#KDE
x <- paste0("outputs/",especie,"_DIVE_KDE.Rdata")
y <- paste0("outputs/",especie,"_HMM_KDE.Rdata")
overlap_kde <- overlap_geral(x,y)

#AKDE
x <- paste0("outputs/",especie,"_DIVE_AKDE.Rdata")
y <- paste0("outputs/",especie,"_HMM_AKDE.Rdata")
overlap_akde <- overlap_geral(x,y)


# INDIVIDUAL LEVEL
# dives
viajeT <- unique(dTDR$ave)
viajeT <- viajeT[-(12:14)]

#HMM
viajeH <- unique(dHMM$ave)
viajeH <- viajeH[-(12:14)]
ov_akde <- NULL

for(i in viajeT){
#TDR
  data1 <- solo_TDR[which(solo_TDR$ID == i),]
  tdr <- data1[,c(1,3,4,7)]
  colnames(tdr)[4] <- "timestamp"
  t_tdr <- as.telemetry(tdr, projection = "+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  #KERNEL
  M.IID_T <- ctmm.fit(t_tdr) # no autocorrelation timescales

  UD0_tdr <- akde(t_tdr,M.IID_T)
  save(UD0_tdr, file = paste("outputs/individual/",especie,"_",i,"DIVE_KDE.Rdata"))
  rm(UD0_tdr)
#HMM
  data1 <- solo_HMM[which(solo_HMM$ID == i),]
  hmm <- data1[,c(1,3,4,7)]
  colnames(hmm)[4] <- "timestamp"
  t_hmm <- as.telemetry(hmm, projection = "+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  #KERNEL
  M.IID_T <- ctmm.fit(t_hmm) # no autocorrelation timescales

  UD0_hmm <- akde(t_hmm,M.IID_T)
  save(UD0_hmm, file = paste("outputs/individual/",especie,"_",i,"HMM_KDE.Rdata"))
  rm(UD0_hmm)
  print(i)
}

# TRIP LEVEL
# dives
viajeT <- unique(dTDR$ID)
viajeT <- viajeT[-(12:14)]

#HMM
viajeH <- unique(dHMM$ID)
viajeH <- viajeH[-(12:14)]
ov_akde <- NULL

for(i in viajeT){
  #TDR
  data1 <- solo_TDR[which(solo_TDR$ID == i),]
  tdr <- data1[,c(1,3,4,7)]
  colnames(tdr)[4] <- "timestamp"
  t_tdr <- as.telemetry(tdr, projection = "+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  #KERNEL
  M.IID_T <- ctmm.fit(t_tdr) # no autocorrelation timescales
  
  UD0_tdr <- akde(t_tdr,M.IID_T)
  save(UD0_tdr, file = paste("outputs/trip/",especie,"_",i,"DIVE_KDE.Rdata"))
  rm(UD0_tdr)
  #HMM
  data1 <- solo_HMM[which(solo_HMM$ID == i),]
  hmm <- data1[,c(1,3,4,7)]
  colnames(hmm)[4] <- "timestamp"
  t_hmm <- as.telemetry(hmm, projection = "+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  #KERNEL
  M.IID_T <- ctmm.fit(t_hmm) # no autocorrelation timescales
  
  UD0_hmm <- akde(t_hmm,M.IID_T)
  save(UD0_hmm, file = paste("outputs/trip/",especie_,i,"HMM_KDE.Rdata"))
  rm(UD0_hmm)
  print(i)
}

#Overlap
source("codes/tesis/funciones/overlap_geral_vi.R")

tipo <- "individual" # or trip
especie <- "PHAET"

files <- list.files(paste0("outputs/",tipo,"/"), full.names = T)
files <- files[which(grepl(especie,files))]

files_DIVE <- files[which(grepl("DIVE",files))]
files_HMM <- files[which(grepl("HMM",files))]

overlap_HMMenTDR <- NULL
overlap_TDRenHMM <- NULL

for(i in 1:length(files_HMM)){
  #  if(i == 5){next}
  print(i)
  x <- files_TDR[i]
  y <- files_HMM[i]
  name <- substr(files_DIVE[i],1,nchar(files_HMM[i])-20)
  result <-  overlap_geral_vi(x,y,name,especie,tipo)
  r1 <- data.frame(o50=result[[1]][2,1],o75=result[[2]][2,1],
                   o95=result[[3]][2,1], ID = names_TDR[i] )
  r2 <- data.frame(o50=result[[1]][1,2],o75=result[[2]][1,2],
                   o95=result[[3]][1,2], ID = names_HMM[i])
  overlap_HMMenTDR <- rbind(overlap_HMMenTDR,r1)
  overlap_TDRenHMM <- rbind(overlap_TDRenHMM,r2)
  remove(result)
} 

median(overlap_TDRenHMM$o50)

write.csv(x = overlap_HMMenTDR,paste0("outputs/",especie,"/",tipo,"/overlap_HMMenTDR.csv"))
write.csv(x = overlap_TDRenHMM,paste0("outputs/",especie,"/",tipo,"/overlap_TDRenHMM.csv"))
