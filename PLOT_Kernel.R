rm(list = ls())
library(RColorBrewer)
library(sftrack)
library(sf)
library(rnaturalearth) #mapas mundo
library(plotly)
library(ggmap)
library(mapdata)
library(maps)
library("rnaturalearthdata")
library(lubridate)
library(rgdal)
library(sf)
library(ggplot2)
library(dplyr)
library(MASS)
library(contoureR)
library(tibble)
library(tidyverse)
library(ggsn)
library(raster)
library(ggpubr)
library(ggpolypath)

#PHAET
datos1 = read.csv("outputs/PAET_4.csv")
vijen <- unique(datos1$ave)

lev <- 0.95
################################################################################
#DIVE
# kde
load("outputs/PHAET_DIVE_KDE.Rdata")
idd_tdr <- UD0_TDR

rm(UD0_TDR)

x = idd_tdr@.Data[[11]]$x
y = idd_tdr@.Data[[11]]$y
z = idd_tdr@.Data[[9]]
xy <- expand.grid(x,y)

araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
a95 <- rasterToContour(araster, levels=c(lev))

coords.x1 <- NULL
coords.x2 <- NULL
for(i in 1:lengths(coordinates(a95))){
  coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
  coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
  coords.x1 <- c(coords.x1,coords.x1_1)
  coords.x2 <- c(coords.x2,coords.x2_1)
}

sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
iid95 <- coordinates(sp95)
iid95 <- as.data.frame(iid95)
iid95$z <- lev

#AKDE
load("outputs/PHAET_DIVE_AKDE.Rdata")
akde_tdr <- UD2w_TDR

x = akde_tdr@.Data[[11]]$x
y = akde_tdr@.Data[[11]]$y
z = akde_tdr@.Data[[9]]
xy <- expand.grid(x,y)

araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
a95 <- rasterToContour(araster, levels=c(lev))

coords.x1 <- NULL
coords.x2 <- NULL
for(i in 1:lengths(coordinates(a95))){
  coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
  coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
  coords.x1 <- c(coords.x1,coords.x1_1)
  coords.x2 <- c(coords.x2,coords.x2_1)
  }
plot(coords.x1,coords.x2)
sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
akdew95 <- coordinates(sp95)
akdew95 <- as.data.frame(akdew95)
akdew95$z <- lev

#HMM
load("outputs/PHAET_HMM_KDE.Rdata")
idd_HMM <- UD0_HMM

rm(UD0_HMM)

x = idd_HMM@.Data[[11]]$x
y = idd_HMM@.Data[[11]]$y
z = idd_HMM@.Data[[9]]
xy <- expand.grid(x,y)

araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
a95 <- rasterToContour(araster, levels=c(lev))

coords.x1 <- NULL
coords.x2 <- NULL
for(i in 1:lengths(coordinates(a95))){
  coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
  coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
  coords.x1 <- c(coords.x1,coords.x1_1)
  coords.x2 <- c(coords.x2,coords.x2_1)
}
plot(coords.x1,coords.x2)

sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
iidHMM95 <- coordinates(sp95)
iidHMM95 <- as.data.frame(iidHMM95)
iidHMM95$z <- lev

#AKDE
load("outputs/PHAET_HMM_AKDE.Rdata")
akde_HMM <- UD2w_HMM

x = akde_HMM@.Data[[11]]$x
y = akde_HMM@.Data[[11]]$y
z = akde_HMM@.Data[[9]]
xy <- expand.grid(x,y)

araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
a95 <- rasterToContour(araster, levels=c(lev))

coords.x1 <- NULL
coords.x2 <- NULL
for(i in 1:lengths(coordinates(a95))){
  coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
  coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
  coords.x1 <- c(coords.x1,coords.x1_1)
  coords.x2 <- c(coords.x2,coords.x2_1)
}
plot(coords.x1,coords.x2)
sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
akdewHMM95 <- coordinates(sp95)
akdewHMM95 <- as.data.frame(akdewHMM95)
akdewHMM95$z <- lev

# Puntos de HMM y TDR
d1 <- datos1 %>% 
  dplyr::select("ave","x","y","datetime","dives")
colnames(d1) <- c("ID","longitude","latitude","timestamp","dives")
d1$dives[which(d1$dives > 0)] <- "buceo_TDR" 
#d1 <- d1[which(d1$ID == vijen[i]),]

d2 <- datos1 %>%
  dplyr::select("ave","x","y","datetime","model")
colnames(d2) <- c("ID","longitude","latitude","timestamp","dives")
d2$dives[which(d2$dives == 2)] <- "forrajeo_HMM"

dTDR <- d1[which(d1$dives == "buceo_TDR"),]
dHMM <- d2[which(d2$dives == "forrajeo_HMM"),]

data1 <- iid95[which(iid95$z == as.character(lev)),]
data2 <- iidHMM95[which(iidHMM95$z == as.character(lev)),]
data3 <- akdew95[which(akdew95$z == as.character(lev)),]
data4 <- akdewHMM95[which(akdewHMM95$z == as.character(lev)),]

#Plots
#KDE
p_s <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-40,-35), ylim = c(-20,-15))+
  north(x.min = -40, x.max = -35,
        y.min = -20, y.max = -15,location = "topleft", symbol = 3)+
geom_polygon(data = data1, aes(x = coords.x1, y = coords.x2),
             alpha = 0.3, show.legend = F, fill = "blue")+
geom_polygon(data=data2,aes(x = coords.x1, y = coords.x2),
               alpha = 0.3, show.legend = F, fill = "red")+
geom_point(data = dHMM, aes(x = longitude, y = latitude),
             col = "red", size = 0.4)+
  geom_point(data = dTDR, aes(x = longitude, y = latitude),
             col = "blue", size = 0.2)+
  ylab("Latitude")+
  xlab("Longitude")+
  geom_polygon(data = islaT,aes(x = V1,y=V2, group = 1), colour = "black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 6),axis.title = element_text(size = 8),
        plot.title = element_text(size = 10))

#AKDE  
p_s2 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-40,-35), ylim = c(-20,-15))+
  ggsn::scalebar(x.min = -40, x.max = -35, border.size = 0.1,
                 y.min = -20, y.max = -15, dist = 50, dist_unit = "km",
                 transform = TRUE, model = "WGS84", st.size = 1.5)+
  geom_polygon(data = data3, aes(x = coords.x1, y = coords.x2),
               alpha = 0.3, show.legend = F, fill = "blue")+
  geom_polygon(data=data4,aes(x = coords.x1, y = coords.x2),
               alpha = 0.3, show.legend = F, fill = "red")+
  geom_point(data = dHMM, aes(x = longitude, y = latitude),
             col = "red", size = 0.4, show.legend = T)+
  geom_point(data = dTDR, aes(x = longitude, y = latitude),
             col = "blue", size = 0.2, show.legend = T)+
  ylab("Latitude")+
  xlab("Longitude")+
 # ggtitle(paste0("AKDE (95%)"))+
  geom_polygon(data = islaT,aes(x = V1,y=V2, group = 1), colour = "black")+
  scale_color_manual(limits=c("Dives", "HMM"),
                     values = c("blue","red")) +
  guides(colour = guide_legend(override.aes = list(pch = c(20, 20),
                                                        col = c("blue", "red"), 
                                                        cex = c(4),
                                                        alpha = 0.3)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = c(0.85,0.85), axis.text = element_text(size = 6),
        axis.title = element_text(size = 8), 
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.1, 'cm'))

figure <- ggarrange(p_s, p_s2 ,
                    ncol = 2, nrow = 1)
ggsave(plot = figure, width = 20,height = 9,units = "cm",
       filename = paste("figures/PAET/kernel_pob_2021.pdf")) 

print(i)

# SULA
datos1 = read.csv("outputs/SULA_4.csv")
vijen <- unique(datos1$ID)
lev <- 0.95

#DIVES
#KDE
  load("outputs/SULA_DIVE_KDE.Rdata")
  idd_tdr <- UD0_TDR

  rm(UD0_TDR)
  
  x = idd_tdr@.Data[[11]]$x
  y = idd_tdr@.Data[[11]]$y
  z = idd_tdr@.Data[[9]]
  xy <- expand.grid(x,y)
  
  araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
  a95 <- rasterToContour(araster, levels=c(lev))
  
  coords.x1 <- NULL
  coords.x2 <- NULL
  group <- NULL
  for(i in 1:lengths(coordinates(a95))){
    coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
    coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
    group <- c(group,rep(i,length(coords.x2_1)))
    coords.x1 <- c(coords.x1,coords.x1_1)
    coords.x2 <- c(coords.x2,coords.x2_1)
  }
  plot(coords.x1,coords.x2)
  
  sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
  sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
  s_iid95 <- coordinates(sp95)
  s_iid95 <- as.data.frame(s_iid95)
  s_iid95$z <- lev
  s_iid95$group <- group
  
#AKDE
  load("outputs/SULA_DIVE_AKDE.Rdata")
  akde_tdr <- UD2w_TDR
  
  x = akde_tdr@.Data[[11]]$x
  y = akde_tdr@.Data[[11]]$y
  z = akde_tdr@.Data[[9]]
  xy <- expand.grid(x,y)
  
  araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
  a95 <- rasterToContour(araster, levels=c(lev))
  
  coords.x1 <- NULL
  coords.x2 <- NULL
  group <- NULL
  for(i in 1:lengths(coordinates(a95))){
    coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
    coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
    group <- c(group, rep(i,length(coords.x2_1)))
    coords.x1 <- c(coords.x1,coords.x1_1)
    coords.x2 <- c(coords.x2,coords.x2_1)
  }
  plot(coords.x1,coords.x2)
  
  sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
  sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
  s_akdew95 <- coordinates(sp95)
  s_akdew95 <- as.data.frame(s_akdew95)
  s_akdew95$z <- lev
  s_akdew95$group <- group
  
#HMM
#KDE
  load("outputs/SULA_HMM_KDE.Rdata")
  idd_HMM <- UD0_HMM
  
  x = idd_HMM@.Data[[11]]$x
  y = idd_HMM@.Data[[11]]$y
  z = idd_HMM@.Data[[9]]
  xy <- expand.grid(x,y)
  
  araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
  a95 <- rasterToContour(araster, levels=c(lev))
  
  coords.x1 <- NULL
  coords.x2 <- NULL
  group <- NULL
  for(i in 1:lengths(coordinates(a95))){
    coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
    coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
    group <- c(group, rep(i,length(coords.x2_1)))
    coords.x1 <- c(coords.x1,coords.x1_1)
    coords.x2 <- c(coords.x2,coords.x2_1)
  }
  plot(coords.x1,coords.x2)
  
  sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
  sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
  s_iidHMM95 <- coordinates(sp95)
  s_iidHMM95 <- as.data.frame(s_iidHMM95)
  s_iidHMM95$z <- lev
  s_iidHMM95$group <- group
  
#AKDE
  load("outputs/SULA_HMM_AKDE.Rdata")
  akde_HMM <- UD2w_HMM
  
  x = akde_HMM@.Data[[11]]$x
  y = akde_HMM@.Data[[11]]$y
  z = akde_HMM@.Data[[9]]
  xy <- expand.grid(x,y)
  
  araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
  a95 <- rasterToContour(araster, levels=c(lev))
  
  coords.x1 <- NULL
  coords.x2 <- NULL
  group <- NULL
  for(i in 1:lengths(coordinates(a95))){
    coords.x1_1 = coordinates(a95)[[1]][[i]][,1]
    coords.x2_1 = coordinates(a95)[[1]][[i]][,2]
    coords.x1 <- c(coords.x1,coords.x1_1)
    coords.x2 <- c(coords.x2,coords.x2_1)
    group <- c(group, rep(i,length(coords.x2_1)))
  }
  plot(coords.x1,coords.x2)
  
  sp95 <- SpatialPoints(cbind(coords.x1,coords.x2), proj4string=CRS("+proj=utm +zone=25 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) 
  sp95 <- spTransform(sp95, CRS("+proj=longlat +datum=WGS84"))
  s_akdewHMM95 <- coordinates(sp95)
  s_akdewHMM95 <- as.data.frame(s_akdewHMM95)
  s_akdewHMM95$z <- lev
  s_akdewHMM95$group <- group
  
  # Puntos de HMM y TDR
  d3 <- datos1 %>% 
    dplyr::select("ID","x","y","datetime","dives")
  colnames(d3) <- c("ID","longitude","latitude","timestamp","dives")
  d3$dives[which(d3$dives > 0)] <- "buceo_TDR" 

  d4 <- datos1 %>%
    dplyr::select("ID","x","y","datetime","model")
  colnames(d4) <- c("ID","longitude","latitude","timestamp","dives")
  d4$dives[which(d4$dives == 2)] <- "forrajeo_HMM"

  dTDR_s <- d3[which(d3$dives == "buceo_TDR"),]
  dHMM_s <- d4[which(d4$dives == "forrajeo_HMM"),]
  
  #Dibujando Kernels
  data5 <- s_iid95[which(s_iid95$z == as.character(lev)),]
  data6 <- s_iidHMM95[which(s_iidHMM95$z == as.character(lev)),]
  data7 <- s_akdew95[which(s_akdew95$z == as.character(lev)),]
  data8 <- s_akdewHMM95[which(s_akdewHMM95$z == as.character(lev)),]
  
  #Plots
  p_su <- ggplot(data = world) +
    geom_sf() +
    coord_sf(c(-40,-37.5), ylim = c(-19,-16.5))+
    north(x.min = -40, x.max = -37.5,
          y.min = -19, y.max = -16.5,location = "topleft", symbol = 3)+
    geom_polypath(data = data5, aes(x = coords.x1, y = coords.x2, group = group),
                 alpha = 0.3, show.legend = F, fill = "blue")+
    geom_polypath(data=data6,aes(x = coords.x1, y = coords.x2, group = group),
                 alpha = 0.3, show.legend = F, fill = "red")+
    geom_point(data = dHMM_s, aes(x = longitude, y = latitude),
               col = "red", size = 0.4)+
    geom_point(data = dTDR_s, aes(x = longitude, y = latitude),
               col = "blue", size = 0.2)+
    ylab("Latitude")+
    xlab("Longitude")+
#    ggtitle(paste0("KDE (95%)"))+
    geom_polygon(data = islaT,aes(x = V1,y=V2, group = 1), colour = "black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 6),axis.title = element_text(size = 8),
          plot.title = element_text(size = 10))
  
  #AKDE
  p_su2 <- ggplot(data = world) +
    geom_sf() +
    coord_sf(c(-40,-37.5), ylim = c(-19,-16.5))+
    geom_polygon(data = data7, aes(x = coords.x1, y = coords.x2, group = group),
                 alpha = 0.3, show.legend = F, fill = "blue")+
    geom_polygon(data=data8,aes(x = coords.x1, y = coords.x2, group = group),
                 alpha = 0.3, show.legend = F, fill = "red")+
    geom_point(data = dHMM_s, aes(x = longitude, y = latitude),
               col = "red", size = 0.4, show.legend = T)+
    geom_point(data = dTDR_s, aes(x = longitude, y = latitude),
               col = "blue", size = 0.2, show.legend = T)+
    ggsn::scalebar(x.min = -40, x.max = -37.5, border.size = 0.1,
                   y.min = -19, y.max = -16.5, dist = 40, dist_unit = "km",
                   transform = TRUE, model = "WGS84", st.size = 1.5)+
    ylab("Latitude")+
    xlab("Longitude")+
    geom_polygon(data = islaT,aes(x = V1,y=V2, group = 1), colour = "black")+
    scale_color_manual(limits=c("Dives", "HMM"),
                       values = c("blue","red")) +
    guides(colour = guide_legend(override.aes = list(pch = c(20, 20),
                                                     col = c("blue", "red"), 
                                                     cex = c(4),
                                                     alpha = 0.3)))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          legend.position = 'none', axis.text = element_text(size = 6),
          axis.title = element_text(size = 8), 
          plot.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.1, 'cm'))
  figure <- ggarrange(p_s, p_s2 , p_su, p_su2,
                      ncol = 2, nrow = 2)

tiff(filename="figures/SULA/kernel_pob_2021.tiff",
       height=5600,width=5200,units='px',res=800,compression='lzw')
figure
dev.off()




 