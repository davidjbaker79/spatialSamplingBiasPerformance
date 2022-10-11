#' Species group correlation intensity
#' 
#' @param A data.frame  including id, X, Y and poc.
#' 
#' @return A data.frame 
#' 
#' @export
dat = sp_in[, c("id", "X", "Y", "poc")]
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

  # Add random
  dat$wgt0 <- 1
  dat <- dat[,!(names(dat) %in% c("biasVar", "biasSeed", "dist"))]
  
}
