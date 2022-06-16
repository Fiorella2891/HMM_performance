get_k_ip <- function(x,n){
  library(CircStats)
  dat_km = Mclust(x,n)
  x$cluster <- dat_km$classification
  kmean_step <- tapply(x$step , x$cluster, mean)
  ksd_step <- tapply(x$step , x$cluster, sd)
  kmean_angle <- tapply(x$angle , x$cluster, circ.mean)
  ksd_angle <- tapply(x$angle , x$cluster, circ.disp)
  ksd_angle<- 1/(unlist(ksd_angle)[seq(4,4*n,by=4)])
  par_k <- list(step = data.frame(mean = kmean_step, sd = ksd_step),
                angle = data.frame(mean = kmean_angle, sd = ksd_angle), 
                cluster = x$cluster)
  return(par_k)
}
