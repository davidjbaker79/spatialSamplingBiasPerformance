#' Create linear features
#'
#' @param df A dataframe with X, Y, Var from \code{generate.enviro}.
#'
#' @import gdistance
#' @import terra
#'
#' @return Dataframe holding the generated environmental data.
create.linear.features <- function(dat) {
  gd <-
    terra::rast(dat, crs = "+proj=longlat +datum=WGS84", type = "xyz")
  N <- nrow(gd) / 5
  
  # Add linear features (e.g. major roads)
  landscape <- 1 - raster::raster(gd) # doesn't yet work using SpatRaster
  #plot(landscape)
  tl <- gdistance::transition(landscape, function(x)
    min(x), 16)
  path1 <-
    shortestPath(tl, c(0, N), c(N, 0), output = "SpatialLines")
  path2 <-
    shortestPath(tl, c(0, N * 0.5), c(N * 0.75, N), output = "SpatialLines")
  path3 <-
    shortestPath(tl, c(N * 0.25, 0), c(N, N * 0.75), output = "SpatialLines")
  path <- do.call(rbind, list(path1, path2, path3))
  path <- vect(path) #plot(gd[[1]]); lines(path, w)
  
  # Aggregate to speed up distance calculation
  gl <- aggregate(gd[[1]], 5)
  pathRast <- terra::rasterize(path, gl)
  pathdd <- terra::distance(pathRast)
  pathdd <- disagg(pathdd, 5, method = "bilinear")
  #
  plot(path)
  
  # Data frame
  pathdd_df <- data.frame(terra::xyFromCell(pathdd, 1:nrow(dat)),
                          terra::values(pathdd))
  names(pathdd_df) <- c("X", "Y", "pathdd")
  pathdd_df$id <- 1:nrow(pathdd_df)
  pathdd_df <- pathdd_df[, c("id", "X", "Y", "pathdd")]
  pathdd_df$pathdd <- round(scales::rescale(pathdd_df$pathdd), 3)
  pathdd_df
  
}

#' Create area features
#'
#' @param df A dataframe with X, Y, Var from \code{generate.enviro}.
#'
#' @import gdistance
#' @import terra
#'
#' @return Dataframe holding the generated environmental data.
#' dat = enviro[c("X","Y","V2")]
create.area.features <- function(dat, n = 5) {
  names(dat)[ncol(dat)] <- "biasVar"
  dat$biasSeed <- NA
  seeds <-
    sample(
      1:nrow(dat),
      size = n,
      replace = FALSE,
      prob = dat$biasVar
    )
  dat$biasSeed[seeds] <- 1
  gd <-
    terra::rast(dat[,c("X","Y","biasSeed")], crs = "+proj=longlat +datum=WGS84", type = "xyz")
  #gd <- terra::raster(gd)
  gd <- terra::buffer(gd, 100000)
  gd <- aggregate(gd, 5)
  gd[gd == 0] <- NA
  padd <- terra::distance(gd)
  padd <- disagg(padd, 5, method = "bilinear")
  plot(sqrt(padd))
  # Data frame
  padd_df <- data.frame(terra::xyFromCell(padd, 1:nrow(dat)),
                        terra::values(padd))
  names(padd_df) <- c("X", "Y", "padd")
  padd_df$id <- 1:nrow(padd_df)
  padd_df <- padd_df[, c("id", "X", "Y", "padd")]
  padd_df$padd <- round(scales::rescale(padd_df$padd), 3)
  padd_df
  
}
