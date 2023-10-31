#' Create blocking (random, systematic, environmental)
#' 
#' Create blocking either randomly or using blockCV.
#' 
#' @param dat A data.frame containing simulated environmental and species data.
#' @param env A data.frame containing simulated environmental data
#' @param block.type A character string giving the type of blocking 
#' required (random, spatial_dist, spatial_systematic, environmental).
#' @param n_folds A numeric giving the number of folds required.
#' 
#' @import sf
#' @import dplyr
#' @import terra
#' @import blockCV
#' 
#' @return A data.frame with blocks added.
#' 
#' @export
#' dat = spBkDat
create.blocks <- function(dat, env, block.type, n_folds) {
  if (block.type == "random") {
    dat$block <-
      sample(1:n_folds, nrow(dat), replace =  TRUE)
  } else {
    # make a sf object from data.frame
    pa_dat <- sf::st_as_sf(dat[, c('X', 'Y', 'det')],
                           coords = c("X", "Y"),
                           crs = "EPSG:4326")
    dat_i <- dplyr::select(env, matches("X|Y|^V"))
    env_dat <- rast(dat_i, type = "xyz")
    crs(env_dat) <- "EPSG:4326"
  
    if (block.type == "spatial_sytematic") {
      bk <- cv_spatial(
        x = pa_dat,
        column = "det",
        hexagon = FALSE,
        r = env_dat,
        k = n_folds,
        size = 900000,
        selection = "systematic",
        biomod2 = TRUE)
      )
    }
    
    if (block.type == "enviro") {
      # environmental clustering
      bk <- cv_cluster(
        x = pa_dat,
        column = "det",
        r = env_dat,
        k = n_folds,
        scale = TRUE
      )
    }
    dat$block <- bk$folds_ids
    dat$block[dat$block %in%
                seq(1, n_folds)[bk$records$test_1 < 10]] <- 0
  }
  dat
}

dat_r <- rast(dat[,c(2:3,10)], type = "xyz")
plot(dat_r)

