#' Create virtual species using using response curve approach
#'
#' Simulate virtual species based on responses to environmental variables.
#'
#' @param x A data.frame holding environmental data.
#' @param pca_axes A numeric giving the number of PCA axes.
#'
#' @return A data.frame with environmental data and the virtual species
#'
#' @export
create.virtual.species.pca <- function(x, pca_axes = 2) {
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

#' Virtual species generating function
#' 
#' Generate virtual species using response curve functions with parameters 
#' generated from random distributions.
#'
#' @param x A data.frame holding environmental data.
#'
#' @import scales
#'
#' @return data.frame holding the species probability of occurrence data.
#' 
#' @export
species.occ <- function(x) {
    # --- Logistic response to e.g. elevation
    beta_a <- 0.5 - runif(1, -0.5, 0.5)
    poc_V1 <-
      1 / (1 + exp((x$V1 - beta_a) / 0.25))# <- 1 * exp(-((x$V1 - 0.5) ^ 2 / (v1_nb ^ 2)))
    # poc_g <-
    #   sapply(seq(0, 1, 0.01), function(x) 1 / (1 + exp((x - beta_a) / 0.25)))
    # ggplot() +
    #   geom_line(aes(y = poc_g, x = seq(0, 1, 0.01))) +
    #   theme_classic() + ylim(0, 1) + xlab("V1") + ylab("Pocc")
    
    # --- Logistic response to e.g. habitat cover
    beta_b <- 0.5 - runif(1, -0.5, 0.5)
    poc_V2 <- 1 / (1 + exp((x$V2 - beta_b) / -0.15))
    # poc_g <-
    #   sapply(seq(0, 1, 0.01), function(x) 1 / (1 + exp((x - beta_b) / -0.15)))
    # ggplot() +
    #   geom_line(aes(y = poc_g, x = seq(0, 1, 0.01))) +
    #   theme_classic() + ylim(0, 1) + xlab("V1") + ylab("Pocc")
    
    # --- Gaussian response to e.g. latitude "temperature" gradient
    sigma_c <- runif(1, 0, 1)
    mu_c <- rnorm(1, 0.5, 0.1)
    poc_V3 <-
      (1 / (sigma_c * sqrt(2 * pi))) * exp(-((x$V3 - mu_c) ^ 2 / (2 * sigma_c ^ 2)))
    # poc_g <-
    #   sapply(seq(0, 1, 0.01), function(x)  ( 1 * exp(-((x - mu_c) ^ 2 / (2 * sigma_c ^ 2)))))
    # ggplot() +
    #   geom_line(aes(y = poc_g, x = seq(0, 1, 0.01))) +
    #   theme_classic() + ylim(0, 1) + xlab("V1") + ylab("Pocc")
    
    # --- Combine data
    nd <- cbind.data.frame(poc_V1, poc_V2, poc_V3)
    
    # --- Combine responses
    nd$poc <- round(nd$poc_V1 * nd$poc_V2 * nd$poc_V3, 2)
    
    # --- Rescale (x - from[1])/diff(from) * diff(to) + to[1]
    nd$poc <- scales::rescale(nd$poc, to = c(0, 1))
    
    # --- Convert to PA
    nd$pa <- rbinom(n = length(nd$poc),
                    size = 1,
                    prob = nd$poc)
    x_i <- nd
    
    # --- Out list
    param <-
      data.frame(
        variable = c("V1", "V2", "V3", "V3"),
        curve = c("Logistic", "Logistic", "Gaussian","Gaussian"),
        param = c("beta", "beta", "sigma", "mu"),
        value = c(beta_a, beta_b, sigma_c, mu_c)
      )
    outlist <-
      list(x_i, param, sum(nd$pa == 1), sum(nd$pa == 1) / nrow(nd))
   
    # ggplot(x_i[x_i$pa == 1, ]) +
    #   geom_point(aes(x = V1, y = V2, colour = poc, alpha = 0.1)) +
    #   scale_colour_viridis_c()
    #
    # ggplot() +
    #   geom_raster(aes(x = X, y = Y, fill = poc), data = x_i) +
    #   geom_point(aes(x = X, y = Y), size = 0.1,colour = "red", data= x_i[x_i$pa == 1, ]) +
    #   scale_fill_viridis_c()
    # 
  
}
