overlap_geral <- function(x,y){

  idd_tdr <- get(load(x))
  idd_HMM <- get(load(y))
  
  # TDR
  x = idd_tdr@.Data[[11]]$x
  y = idd_tdr@.Data[[11]]$y
  z = idd_tdr@.Data[[9]]
  xy <- expand.grid(x,y)
  
  #HMM
  x2 = idd_HMM@.Data[[11]]$x
  y2 = idd_HMM@.Data[[11]]$y
  z2 = idd_HMM@.Data[[9]]
  xy2 <- expand.grid(x2,y2)
  
  ## pasando a raster
  araster = rasterFromXYZ(data.frame(x=xy$Var1,y=xy$Var2,z=matrix(z,ncol=1)))
  araster2 = rasterFromXYZ(data.frame(x=xy2$Var1,y=xy2$Var2,z=matrix(z2,ncol=1)))
  
  #Inspeccion
  # plot(araster2, xlim = c(-4e+05,2e+05))
  # contour(araster, add = T, level = c(0.95))
  # contour(araster2, add = T, level = c(0.95), col = "red")
  # #obteniendo contornos
  con <- c(0.5, 0.75, 0.95)
  
  #Overlapp a cada nivel
  ker_con <- list()
  for(f in 1:length(con)){
    # Contour 95
    aa <- rasterToContour(araster, levels = con[f])
    aa2 <- rasterToContour(araster2, levels = con[f])
    
    at1 <- NULL
    at2 <- NULL
    aint <- NULL
    
    for(d in 1:lengths(coordinates(aa))){
      for(a in 1:lengths(coordinates(aa2))){
    coords.x1 = coordinates(aa)[[1]][[d]][,1]
    coords.x2 = coordinates(aa)[[1]][[d]][,2]
    coords_t <- data.frame(x = coords.x1, y = coords.x2, ID = "TDR")
    
    coords.x12 = coordinates(aa2)[[1]][[a]][,1]
    coords.x22 = coordinates(aa2)[[1]][[a]][,2]
    coords_t2 <- data.frame(x = coords.x12, y = coords.x22, ID = "HMM")
    
    coords <- rbind(coords_t,coords_t2)
    polys <- lapply(unique(coords$ID), function(i) {
      Polygons(list(Polygon(coords[coords$ID==i, 1:2])), ID=i)
    })
    spa_polys <- SpatialPolygons(polys)
    
    areas_poly <- vector(length = length(spa_polys)) 
    for (x in seq_along(spa_polys)) areas_poly[x]<-area(spa_polys[x])
    ## areas for the overlapping regions:
    idx <- combn(seq_along(spa_polys),2)
    areas_intersect <- sapply(1:ncol(idx), function(x) {
      interseccion <- intersect(spa_polys[idx[1,x]], spa_polys[idx[2,x]]) 
      
      if(is.null(interseccion)){
        return(0)
      }else{
        return(area(interseccion))
      }
    })
    
    at1 <- c(at1,areas_poly[idx[1,1]]) #TDR
    at2 <- c(at2,areas_poly[idx[2,1]]) #HMM
    aint <- c(aint, areas_intersect)
    }
      }
    #####
    ## get overlaps in percentages:
    areas_intersect = sum(aint)
    areas_poly <- c(sum(unique(at1)),sum(unique(at2)))
    print(areas_poly)
    overlap_perc <- 
      round(do.call(cbind, lapply(seq_len(ncol(idx)), function(x)  
        rbind(
          areas_intersect[x] / areas_poly[idx[1,x]] * 100, 
          areas_intersect[x] / areas_poly[idx[2,x]] * 100
        )
      )), 2)
    
    ## into matrix form:
    m <- matrix(100, ncol=length(spa_polys), nrow=length(spa_polys))
    m[rbind(t(idx),t(idx)[,2:1])] <- as.vector(t(overlap_perc))
    
    ker_con[[f]] <- m
  }
  return(ker_con)
}
