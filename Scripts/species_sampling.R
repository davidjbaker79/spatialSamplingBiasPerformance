#' Species virtual sampling
#'
#' @param sp_dat Data frame return by \code{species.occ}.
#' @param bias Vector of weights returned by \code{samp.bias.Landcov}.
#' @param nPres Numeric indicating the number of occurrence records to
#' be sampled.
#' @param detProb Numeric indicating the probability of detecting a
#' species on a survey visit.
#' @param fpr Numeric indicating the rate at which species are miss identified.
#' @param listProb Numeric indicating the probability that a complete list
#' is returned.
#'
#' @return Data frame with species det column added.
#'
#' @export
species.sampling <- function(sp_dat,
                             bias = 1,
                             nPres = 200,
                             detProb = 1,
                             fpr = 0,
                             listProb = 1) {
  # Join bias
  sp_dat <- cbind(sp_dat, bias)
  
  # Prevalence
  pr <- sp_dat[sp_dat$pa == 1,]
  prev <- nrow(pr) / nrow(sp_dat)
  n_absen <- nPres / prev
  
  # Assign nPres detections based on biased surface
  pr <-
    pr[sample(1:nrow(pr), nPres, prob = pr$bias, replace = TRUE), ]
  ab <- sp_dat[sp_dat$pa == 0,]
  ab <-
    ab[sample(1:nrow(ab),
              n_absen,
              prob = ab$bias,
              replace = TRUE), ]
  
  # Combine data
  sp_dat_s <- rbind(pr, ab)
  
  # Assign detections based on detection probability
  sp_dat_s$det <- sp_dat_s$pa * rbinom(nrow(sp_dat_s), 1, detProb)
  
  # Adjust detections based on probability of completing list
  sp_dat_s$det <- sp_dat_s$det * rbinom(nrow(sp_dat_s), 1, listProb)
  
  # Adjust detections based on probability of incorrect id
  det_id <- (1 - sp_dat_s$pa) * rbinom(nrow(sp_dat_s), 1, fpr)
  sp_dat_s$det <-  sp_dat_s$det + det_id
  
  # Return
  sp_dat_s
  
}
