#' Create virtual species using PCA species approach
#'
#' @param x A data.frame holding environmental data. 
#' @param pca_axes A numeric giving the number of PCA axes.
#' 
#' @return A data.frame with environmental data and the virtual species
#' 
#' @export
create.virtual.species <- function(x, pca_axes) {
  pca_res <-
    stats::prcomp(x,
                  center = TRUE,
                  scale = TRUE)
  pca_dat <- pca_res$x[, 1:pca_axes]
  # Gaussian responses to PCA
  gaus_m <- rnorm(pca_axes, 0, 0.5)
  gaus_sd <- abs(rnorm(pca_axes, 0, 0.5))
  for (i in 1:pca_axes) {
    pca_dat[, i] <-
      dnorm(pca_dat[, i], mean = gaus_m[i], sd = gaus_sd[i])
  }
  # Calculate suitability
  poc <- apply(pca_dat, 1, function(x)
    prod(x))
  poc <- round((poc - min(poc)) / (max(poc) - min(poc)), 4)
  # Convert to PA
  pa <- rbinom(n = length(poc),
               size = 1,
               prob = poc)
  # Join with XY coordinates
  out <- cbind(poc, pa)
}
