#' Create spatial sampling bias surface for sampling based on 
#' distance decay to road features.
#' 
#' @param dat A data.frame containing distance decay values for grid cells to 
#' nearest road.
#' 
#' @return 
#' 
#' @export
samp.bias.feature.dd <- function(dat) {
  
  names(dat)[ncol(dat)] <- "biasVar"
  # Calculate distance decay from points
  # wgt <- rast(dat[, c("X", "Y", "biasVar")], type = "xyz",
  #             crs = "EPSG:4326")
  # dist <- terra::distance(wgt)
  # dist <- (1 - (dist / minmax(dist)[2,]))
  
  # Create logistic decay
  dat$dist <- dat$biasVar
  for (B in c(0.6, 0.8, 1)) {
    dat[paste0("wgt", sub("\\.", "", B))] <-
      scales::rescale(1 / (1 + exp((dat$dist - B) / -0.05)))
  }
  # Add random
  dat$wgt0 <- 1
  # ggplot(dat, aes(x = X, y = Y)) +
  #   geom_raster(aes(fill = wgt06)) +
  #   geom_point(data = dat[dat$biasSeed == 1,]) +
  #   scale_fill_viridis_b() + theme_bw()
  dat <- dat[,!(names(dat) %in% c("biasVar", "biasSeed", "dist"))]
  
}

#' Create spatial sampling bias surface for sampling background data based
#' buffers from occurrence records
#' 
#' @param spDat A data.frame with 
#' @param grd A data.frame with coordinates of study grid
#' @param bf_w A numeric giving the distance for the buffer
#' 
#' @return 
#' 
#' @export
samp.bias.occ.buff <- function(dat, grd, buf_w) {
  # Calculate distance decay from points
  grd_r <- rast(grd, type = "xyz", crs = "EPSG:4326")
  # sf to buffer
  occ_sf <- st_as_sf(dat, coords = c("X", "Y"), crs = "EPSG:4326")
  # Buffer
  occ_b <- st_buffer(occ_sf, buf_w)
  occ_b <- st_union(occ_b)
  # rasterize
  occ_r <- rasterize(vect(occ_b), grd_r)
  occ_r[is.na(occ_r)] <- 0
  occ_r <-  occ_r / ncell(occ_r)
  # Extract
  as.data.frame(occ_r)[, 1]#, xy = TRUE, cells = TRUE)
  
}






