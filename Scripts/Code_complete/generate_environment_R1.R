#' Generate environmental variables with spatial autocorrelation
#'
#' Simulate environmental variables for generating virtual species.
#'
#' @param n Gives the square dimensions of the model arena (n x n)
#' as the number of cells * 10
#' @param phi A vector defining the strength of autocorrelation. The length of
#' the vector determines the number of surfaces to return.
#'
#' @import mgcv
#' @import scales
#' @import Rfast
#'
#' @return Dataframe holding the generated environmental data.
#'
#' @export
generate.enviro <- function(n = 10, phi = c(0.1, 0.1)) {
  # Set up a square lattice region
  simgrid <- expand.grid(1:n, 1:n)
  n <- nrow(simgrid)
  
  # Set up distance matrix
  distance <- as.matrix(Rfast::Dist(simgrid))
  
  # Generate random variable
  env_df <- simgrid[, 1:2] - 0.5
  names(env_df) <- c("X", "Y")
  for (i in 1:length(phi)) {
    env_df[[paste0("V", i)]] <-
      scales::rescale(mgcv::rmvn(1, rep(0, n), exp(-phi[i] * distance)))
  }
  
  # Visualize results
  gd <- rast(env_df, type = "xyz",crs = "EPSG:4326")
  gd <- disagg(gd, 5, method = "bilinear")
  
  # Data frame
  env_df <- data.frame(crds(gd), terra::values(gd))
  names(env_df) <- c("X", "Y", names(gd))
  env_df$id <- 1:nrow(env_df)
  env_df <- env_df[, c("id", "X", "Y", paste0("V", 1:length(phi)))]
  
  # Return
  env_df
  
}
