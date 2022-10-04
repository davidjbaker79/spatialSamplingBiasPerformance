#' Process data for comparisons between corrected and uncorrected SDMs
#'
#' @param dat A data.frame containing model results.
#' @param mod.qual A character vector indicating whether to extract the best or
#' worse performing model
#' @param n_pres A numeric giving the results to be analysed in terms of the
#' number of presences.
#' @param metric A character string giving the test statistic.
#' @param ssbTrue A character string giving the true ssb to be analysed.
#' @param calStdEff Logical indicating whether to use standardised metrics
#' i.e. Hedge g.
#'
#' @return data.frame
#' 
#' @export
process.comparision.data <-
  function(dat,
           mod.qual = "best",
           n_pres,
           metric,
           ssbTrue,
           calStdEff = TRUE) {
    #--- First, filter and aggregate cross validation to m and sd
    dat_f <- dat[dat$nPres %in% n_pres,]
    dat_f <- dat_f[dat_f$nicheSsbCor == "none",]
    dat_f <- dat_f[dat_f$ssb_true == ssbTrue,]
    dat_f <- dat_f[dat_f$metric %in% metric,]
    dat_f <- dat_f[!is.na(dat_f$value),]
    dat_f <- dat_f %>%
      group_by(ssb_correct,
               spId,
               metric,
               ssb_true,
               nPres,
               testDat,
               ssb_meth,
               blkType ,
               nicheSsbCor) %>%
      summarise(
        m = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        N = n()
      ) %>%
      ungroup() %>%
      dplyr::select(!c(nicheSsbCor, ssb_true))
    
    #--- Identify the best performing model via internal cv
    best_int <- dat_f %>%
      filter(testDat == "int",
             ssb_correct != 0 )
    if (mod.qual == "best") {
      best_int <- best_int %>%
        group_by(spId, metric, nPres, testDat, ssb_meth, blkType) %>%
        slice(which.max(m))
    } else if (mod.qual == "worse") {
      best_int <- best_int %>%
        group_by(spId, metric, nPres, testDat, ssb_meth, blkType) %>%
        slice(which.min(m))
    } else {
      best_int <- best_int %>%
        group_by(spId, metric, nPres, testDat, ssb_meth, blkType) %>%
        slice_sample()
    }
    best_int <- best_int %>%
      rename(m_int = "m",
             sd_int = "sd",
             N_int = "N") %>%
      ungroup() %>%
      dplyr::select(!c(testDat))
    
    #--- Find corresponding independent evaluation
    best_ind <- dat_f %>%
      filter(testDat == "ind",
             ssb_correct != 0) %>%
      rename(m_ind = "m",
             sd_ind = "sd",
             N_ind = "N") %>%
      dplyr::select(!c(testDat)) %>%
      left_join(best_int,
                .,
                by =
                  c(
                    "spId",
                    "metric",
                    "nPres",
                    "ssb_meth",
                    "ssb_correct",
                    "blkType"
                  )) %>%
      dplyr::select(!c(m_int, sd_int, N_int)) %>%
      rename(m = "m_ind", sd = "sd_ind", N = "N_ind")
    
    #--- Best adjusted models
    best_int <- best_int %>%
      rename(m = "m_int", sd = "sd_int", N = "N_int") %>%
      mutate(testDat = "int")
    best_ind <- best_ind %>%
      mutate(testDat = "ind")
    best_adj <- bind_rows(best_int, best_ind) %>%
      mutate(ssb_correct = "Yes") %>%
      distinct()
    
    #--- No adjustment data
    no_Adj <- dat_f %>%
      filter(ssb_correct %in% c("wgt0", 0)) %>%
      mutate(ssb_correct = "No") %>%
      distinct()
    
    #--- Combine adjustment and no adjustment
    dat_comp <- bind_rows(best_adj, no_Adj)
    if (calStdEff) {
      dat_comp <- dat_comp %>%
        pivot_wider(names_from = ssb_correct,
                    values_from = c(m, sd, N)) %>%
        mutate(
          d_m = m_Yes - m_No,
          A = sd_Yes ^ 2 / N_Yes,
          B = sd_No ^ 2 / N_No,
          d_sd = sqrt(A + B)
        ) %>%
        dplyr::select(spId, metric, nPres, ssb_meth,
                      blkType, testDat, d_m, d_sd)
    }
    dat_comp
  }
