rm(list=ls())
library('lubridate')
library("maps")
library("dplyr")
library("moveHMM")
library("factoextra")
library(caret)
library(mapproj)
library(mclust)
library("ggplot2")
source("codes/tesis/funciones/fun_parini.R")

# HIDDEN MARKOV MODELS
################################################################################
# 1. Phaethon aethereus
files <- list.files("data/PHAET/",full.names = T)

baseR <- NULL
for(i in 1:length(files)){

  data <- read.csv(files[i], header = T)
  data <- data[,-1]
  datab <- data[,c(1:5,9)]
    baseR = rbind(baseR,datab)
}

baseR$datetime <- as.POSIXct(strptime(baseR$datetime,
                                      format="%Y-%m-%d %H:%M:%S"),tz="GMT")
colnames(baseR) <- c("ave","ID","datetime","lon","lat","dives")

data1 <- baseR
data1$ID <- as.factor(data1$ID)
data.HMM<-prepData(data1, type="LL", coordNames=c("lon","lat"))

# Initial parameters
#Parametros Mclust
set.seed(123)
datakm = data.HMM[2:(nrow(data.HMM)-1),c(2:3)]
datakm <- datakm[-c(which(is.na(datakm$step)),which(is.na(datakm$angle))),]

n_estados <- 5 #From 2 to 5
in_par <- get_k_ip(datakm,n_estados)
stepPar0 <- c(in_par$step$mean, in_par$step$sd)
anglePar0 <- c(in_par$angle$mean, in_par$angle$sd)

# Values
ini_par = list(
  e2 = list(stepPar0 <-  c(11,   2,    10,   2),
            anglePar0 <- c( 0,   0,   30,   6)),
  e3 = list(stepPar0 <-  c(11, 1.4,  0.2, 2.3, 1.6, 0.1),
            anglePar0 <- c( 0,   0,    0,  38,   4,  23)),
  e4 = list(stepPar0 <-  c(11,   5,  0.1,   2,   4,   2,  0.1,    1),
            anglePar0 <- c( 0,   0,    0,   0,  30,   2,   10,    6)),
  e5 = list(stepPar0 <-  c(11,   5,  0.1,   2, 0.5,   4,    2,  0.1,  1, 0.5),
            anglePar0 <- c( 0,   0,    0,   0,   0,  30,    2,   10,  6,   5)))

s <- seq(2:5)
for(s in seq(2,5,1)){
stepPar0 <- unlist(ini_par[[s-1]][1])
anglePar0 <- unlist(ini_par[[s-1]][2])

model1 <- fitHMM(data = data.HMM, nbStates = s, stepPar0 = stepPar0,
                  anglePar0 = anglePar0, formula = ~1)

estado <- paste0(s,"_est")
print(paste(AIC(model1),s))

estados_vit <- viterbi(model1)
estados_vit <- as.factor(estados_vit)

model <- estados_vit
base_r <- cbind(data.HMM,model)

write.csv(base_r,paste0("outputs/PAET_",s,".csv"))

}

# 2. Sula leucogaster
files <- list.files("data/SULA/",full.names = T)

baseR <- NULL
for(i in 1:length(files)){
  
  data <- read.csv(files[i], header = T)
  data <- data[,-1]
  datab <- data[,c(1:5,9)]
  baseR = rbind(baseR,datab)
}

baseR$datetime <- as.POSIXct(strptime(baseR$datetime,
                                      format="%Y-%m-%d %H:%M:%S"),tz="GMT")
colnames(baseR) <- c("ave","ID","datetime","lon","lat","dives")

data1 <- baseR
data1$ID <- as.factor(data1$ID)
data.HMM<-prepData(data1, type="LL", coordNames=c("lon","lat"))

# Initial parameters
#Parametros Mclust
datakm = data.HMM[2:(nrow(data.HMM)-1),c(2:3)]
datakm <- datakm[-c(which(is.na(datakm$step)),which(is.na(datakm$angle))),]

set.seed(123)
obj <- Mclust(datakm)
in_par <- get_k_ip(datakm,5)
stepPar0 <- c(in_par$step$mean, in_par$step$sd)
anglePar0 <- c(in_par$angle$mean, in_par$angle$sd)

ini_par = list(
  e2 = list(stepPar0 <-  c(0.6, 0.15, 0.4,  0.15),
            anglePar0 <- c(  0,    0,  64,     3)),
  e3 = list(stepPar0 <-  c(0.7, 0.01, 0.3,   0.4, 0.01, 0.2),
            anglePar0 <- c(  0,    0,   0,    40,   10,  80)),
  e4 = list(stepPar0 <-  c(0.7, 0.03, 0.3, 0.001,  0.7,  0.02,  0.2,   0.001),
            anglePar0 <- c(  0,    0,   0,     0,   40,     7,   40,    5)),
  e5 = list(stepPar0 <-  c(0.7, 0.07, 0.4, 0.01, 0.5,  
                           0.3, 0.05, 0.23, 0.05, 0.2),
            anglePar0 <- c( 0,    0,   0,     0,   0,     
                            60, 1,  10,  40,  40)))

s <- seq(2:5)
for(s in seq(2,5,1)){
  stepPar0 <- unlist(ini_par[[s-1]][1])
  anglePar0 <- unlist(ini_par[[s-1]][2])
  
  model1 <- fitHMM(data = data.HMM, nbStates = s, stepPar0 = stepPar0,
                   anglePar0 = anglePar0, formula = ~1)
  
  estado <- paste0(s,"_est")
  
  print(paste(AIC(model1),s))
  estados_vit <- viterbi(model1)
  estados_vit <- as.factor(estados_vit)
  
  model <- estados_vit
  base_r <- cbind(data.HMM,model)
  write.csv(base_r, paste0("outputs/SULA_",s,".csv"))
}

# ASSESING PERFORMANCE BEST MODEL (AIC)
################################################################################
especie <- "SULA" #or PHAET
indi <- read.csv(paste0("outputs/",especie,"_4.csv"))
indi <- indi %>% 
  dplyr::select("ave","x","y","dives","model","ID","datetime")

indi$datetime <- as.POSIXct(strptime(indi$datetime,
                                     format="%Y-%m-%d %H:%M:%S"),tz="GMT")

# POPULATION LEVEL
indi$dives[indi$dives > 0] = 1
indi$dives[indi$dives == 0] = 2
indi$dives <- as.factor(indi$dives)
levels(indi$dives) <- c("Dives", "No dives")

indi$model[indi$model != 2] = 5
indi$model <- as.factor(indi$model)
levels(indi$model) <- c("Dives", "No dives")

a <- confusionMatrix(indi$model, indi$dives)
A <- a$table[1]
B <- a$table[3]
C <- a$table[2]
D <- a$table[4]
sensitivity = A/(A+C)
specificity = D/(B+D)
prevalence = (A+C)/(A+B+C+D)
PPV = (sensitivity * prevalence)/((sensitivity*prevalence) + ((1-specificity)*(1-prevalence)))
NPV = (specificity * (1-prevalence))/(((1-sensitivity)*prevalence) + ((specificity)*(1-prevalence)))
accuracy = (sensitivity+specificity)/2
precision = A/(A+B)
values <- data.frame(ave = indi$ID[1], PPV=PPV,NPV=NPV,accuracy=accuracy,
                     precision=precision,sensitivity=sensitivity,
                     specificity=specificity)
write.csv(values,paste0("outputs/",especie,"population_matrix_conf.csv"))

# INDIVIDUAL AND TRIP LEVELS
indi$ID <- as.factor(indi$ID)
indi$ave <- as.factor(indi$ave)

datos <- NULL
for(a in levels(indi$ave)){
  indi2 <- indi[which(indi$ave == a),] 
  a <- confusionMatrix(indi2$model, indi2$dives)
  A <- a$table[1]
  B <- a$table[3]
  C <- a$table[2]
  D <- a$table[4]
  sensitivity = A/(A+C)
  specificity = D/(B+D)
  prevalence = (A+C)/(A+B+C+D)
  PPV = (sensitivity * prevalence)/((sensitivity*prevalence) + ((1-specificity)*(1-prevalence)))
  NPV = (specificity * (1-prevalence))/(((1-sensitivity)*prevalence) + ((specificity)*(1-prevalence)))
  accuracy = (sensitivity+specificity)/2
  precision = A/(A+B)
  values <- data.frame(ave = indi2$ID[1], PPV=PPV,NPV=NPV,accuracy=accuracy,
                       precision=precision,sensitivity=sensitivity,
                       specificity=specificity)
  datos <- rbind(datos,values)
  
}
write.csv(datos,paste0("outputs/",especie,"_individual_matrix_conf.csv"))
write.csv(datos,paste0("outputs/",especie,"_trip_matrix_conf.csv"))

# Plots
# Individual
phaet <- read.csv("outputs/PAET_individual_matrix_conf.csv")
phaet <- phaet[,-1]
sula <- read.csv("outputs/SULA_individual_matrix_conf.csv")
sula <- sula[,-1]

phaet2 <- phaet %>% 
  dplyr::select(ave,accuracy,precision,sensitivity,specificity) %>% 
  gather(key = "accuracy:specificity", value = "indicador",-"ave")
colnames(phaet2) <- c("ave", "ind","value")
phaet2$especie <- "Red-billed tropicbird"

sula2 <- sula %>% 
  dplyr::select(ave,accuracy,precision,sensitivity,specificity) %>% 
  gather(key = "accuracy:specificity", value = "indicador",-"ave")
colnames(sula2) <- c("ave", "ind","value")
sula2$especie <- "Brown booby"

base <- rbind(phaet2,sula2)
base$ind <- as.factor(base$ind)
levels(base$ind) <- c('Acc','Pre','Sen','Spe')

gg_ind <- ggplot(data = base, aes(x = ind, y = value))+
  geom_boxplot(aes(fill = especie))+
  ylim(0,1)+
  xlab("")+
  ylab("Values")+
  theme_bw()+
  labs(fill = "")+
  ggplot2::annotate(geom="text", x=0.9, y=.99, label="Individual",
                    color=1, size = 3.5)+
  theme(axis.text.x = element_text(face = "bold", size = 9), legend.position = "none",
        axis.title.x = element_blank(), axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8))

#Trip
phaet <- read.csv("outputs/PAET_trip_matrix_conf.csv")
phaet <- phaet[,-1]
sula <- read.csv("outputs/SULA_trip_matrix_conf.csv")
sula <- sula[,-1]

phaet2 <- phaet %>% 
  dplyr::select(ave,accuracy,precision,sensitivity,specificity) %>% 
  gather(key = "accuracy:specificity", value = "indicador",-"ave")
colnames(phaet2) <- c("ave", "ind","value")
phaet2$especie <- "Red-billed tropicbird"

sula2 <- sula %>% 
  dplyr::select(ave,accuracy,precision,sensitivity,specificity) %>% 
  gather(key = "accuracy:specificity", value = "indicador",-"ave")
colnames(sula2) <- c("ave", "ind","value")
sula2$especie <- "Brown booby"

base <- rbind(phaet2,sula2)
base$ind <- as.factor(base$ind)
levels(base$ind) <- c('Acc','Pre','Sen','Spe')
gg_trips <- ggplot(data = base, aes(x = ind, y = value))+
  geom_boxplot(aes(fill = especie))+
  ylab("")+
  ylim(0,1)+
  theme_bw()+
  labs(fill = "")+
  ggplot2::annotate(geom="text", x=0.7, y=.99, label="Trip",
                    color=1, size = 3.5)+
  theme(axis.text.x = element_text(face = "bold", size = 9),
        legend.position = 'none',
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),axis.text.y = element_text(size = 7),
        legend.text = element_text(size = 7))
figure <- ggarrange(gg_ind, gg_trips,
                    ncol = 2, nrow = 1)

tiff(filename="figures/box_final.tiff",
     height=2600,width=5900,units='px',res=800,compression='lzw')
figure
dev.off()

