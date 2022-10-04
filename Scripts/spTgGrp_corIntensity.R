#' Species group correlation intensity
#' 
#' @param A data.frame  including id, X, Y and poc.
#' 
#' @return A data.frame 
#' 
#' @export
spTgGrp.corIntensity <- function(dat) {
  
  names(dat)[ncol(dat)] <- "biasVar"
  # Calculate distance decay from points
  wgt <- rast(dat[, c("X", "Y", "biasVar")], type = "xyz",
              crs = "EPSG:4326")
  wgt <- aggregate(wgt, 50, mean)
  wgt <- disagg(wgt, 50, method = "bilinear")
  wgt <- (wgt / minmax(wgt)[2,])

  # Create logistic decay
  dat$dist <- values(wgt)[,1]
  for (B in c(0.6, 0.8, 1)) {
    dat[paste0("wgt", sub("\\.", "", B))] <-
      scales::rescale(1 / (1 + exp((dat$dist - (B^3)) / -0.1)))
  }
  
  # B = 1
  # x <- seq(0,1, 0.1)
  # y <- 1 / (1 + exp((x - (B^3)) / -0.1))
  # ggplot() + geom_point(aes(x = x, y = y))
  # 
  # Add random
  dat$wgt0 <- 1
  # ggplot(dat, aes(x = X, y = Y)) +
  #   geom_raster(aes(fill = wgt1)) +
  #   geom_point(data = dat[dat$biasSeed == 1,]) +
  #   scale_fill_viridis_c() + theme_bw()
  dat <- dat[,!(names(dat) %in% c("biasVar", "biasSeed", "dist"))]
  
}
