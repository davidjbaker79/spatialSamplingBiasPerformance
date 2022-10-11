#' Create blocking (random, systematic, environmental)
#' 
#' @param dat A data.frame containing simulated enviro and species data.
#' @param block.type A character string giving the type of blocking 
#' required (random, spatial_dist, spatial_systematic, environmental).
#' @param n_folds A numeric giving the number of folds required.
#' 
#' @return A data.frame with blocks added.
#' 
#' @export
create.blocks <- function(dat, env, block.type, n_folds) {
  # Create blocking
  if (block.type == "random") {
    dat$block <-
      sample(1:n_folds, nrow(dat), replace =  TRUE)
  } else {
    # make a sf object from data.frame
    pa_dat <- sf::st_as_sf(dat[, c('X', 'Y', 'det')],
                           coords = c("X", "Y"),
                           crs = "EPSG:4326")
    dat_i <- dplyr::select(env, matches("X|Y|V"))
    env_dat <- rasterFromXYZ(dat_i)
    crs(env_dat) <- "EPSG:4326"
    
    if (block.type == "spatial_dist") {
      bk <- spatialBlock(
        speciesData = pa_dat,
        species = "det",
        theRange = 5000,
        k = n_folds,
        selection = "random",
        iteration = 100,
        numLimit = NULL,
        biomod2Format = FALSE,
      )
    }
    
    if (block.type == "spatial_sytematic") {
      bk <- spatialBlock(
        speciesData = pa_dat,
        species = "det",
        rasterLayer = env_dat,
        rows = 5,
        cols = 8,
        k = n_folds,
        selection = "systematic",
        biomod2Format = FALSE,
        xOffset = 0.1
      )
    }
    
    if (block.type == "enviro") {
      # environmental clustering
      bk <- envBlock(
        rasterLayer = env_dat,
        speciesData = pa_dat,
        species = "det",
        k = n_folds,
        standardization = "normal",
        rasterBlock = FALSE
      )
    }
    dat$block <- bk$foldID
    dat$block[dat$block %in%
                seq(1, n_folds)[bk$records$test_1 == 0]] <- 0
  }
  dat
}
