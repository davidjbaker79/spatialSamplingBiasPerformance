#' Collate and process simulation results
#' 
#' From simulation output folder, collate and process simulation results into 
#' data.frame.
#'
#' @param nicheCor Logical Whether or not niche is correlated with SSB
#'
#' @import snowfall
#' @import tidyverse
#' @import data.table
#' 
#' @return data.frame
#' 
#' @export
prepare_sim_data <- function(nicheCor = FALSE) {
  # -- Uncorrelated niches
  if (nicheCor) {
    pat <- "corNiche"
  } else{
    pat <-  "_.RData"
  }
  f_sim_res <-
    list.files(
      "Outputs/simulation_res_optimization_R1/",
      pattern = pat,
      full.names = TRUE
    )
  
  sfInit(parallel = TRUE, cpus = 7)
  sfExport("f_sim_res")
  f_sim_res_l <- sfLapply(f_sim_res, function(i) {
    print(i)
    possibleError <- tryCatch(
      sp_i  <- get(load(i)),
      error = function(e)
        e
    )
    if (inherits(possibleError, "error")) {
      file.remove(i)
      return(NULL)
    } else {
      return(sp_i)
    }
    
  })
  sim_res <- do.call(bind_rows, f_sim_res_l)
  
  #--- Process factor levels and names
  sim_res_i <- sim_res %>%
    # Split out the test data type and metric
    rowwise %>%
    mutate(
      testDat = strsplit(metric, split = "\\.")[[1]][2],
      metric = strsplit(metric, split = "\\.")[[1]][1],
    ) %>%
    # Set the factor levels that we need and recode names
    mutate(
      metric = factor(metric,
                      levels =  c("AUC", "boyce", "corP", "corS")),
      metric = fct_recode(
        metric,
        AUC = "AUC",
        Boyce = "boyce",
        "S (Pearson's)" = "corP",
        "S (Spearman's)" = "corS"
      ),
      ssb_true = factor(ssb_true,
                        levels = c("wgt0",
                                   "wgt06",
                                   "wgt08",
                                   "wgt1")),
      ssb_meth = factor(
        ssb_meth,
        levels = c("spTgGrp",
                   "featDD",
                   "occBuff",
                   "occFilter",
                   "covCond")
      ),
      testDat = factor(testDat,
                       levels = c("int",
                                  "ind")),
      blkType  = factor(blkType ,
                        levels = c(
                          "random",
                          "spatial_sytematic",
                          "enviro"
                        )),
      nPres = factor(nPres,
                     levels = c("50",
                                "100",
                                "200",
                                "400",
                                "800",
                                "1600")),
      
    ) %>%
    pivot_longer(cols = starts_with("V")) # Pivot to a longer table format
  
  #--- First, filter and aggregate cross validation to m and sd
  sim_res_i <- setDT(sim_res_i)[, .(
    m = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    N = .N
  ), by = .(ssb_correct,
            spId,
            prevelance,
            metric,
            ssb_true,
            nPres,
            testDat,
            ssb_meth,
            blkType)]
  
  
}

#' Process data for comparisons between corrected and uncorrected SDMs
#' 
#' Find 'best' (or other) model based on internal cv for each SSB method and
#' treatment from amongst the multiple implementations and then find the 
#' corresponding independent test statistic.
#'
#' @param dat A data.frame containing model results.
#' @param mod.qual A character vector indicating whether to extract the best or
#' worse performing model
#' @param calStdEff Logical indicating whether to use standardised metrics
#' i.e. Hedge g.
#'
#' @return data.frame
#' 
#' @export
process.comparision.data <-
  function(dat_f,
           mod.qual = "best",
           calStdEff = FALSE) {
    
    dat_f$prevelance <- dat_f$prevelance / 62500
    dat_f$prevCat <-
      cut(
        dat_f$prevelance,
        breaks = seq(0, 0.6, 0.1),
        labels = seq(0.1, 0.6, 0.1)
      )

    #--- Identify the best performing model via internal cv
    best_int <- dat_f %>%
      filter(testDat == "int",
             ssb_correct != 0 )
    if (mod.qual == "best") {
      best_int <- best_int %>%
        group_by(spId, metric, ssb_true, nPres, prevCat, testDat, ssb_meth, blkType) %>%
        slice(which.max(m))
    } else if (mod.qual == "worst") {
      best_int <- best_int %>%
        group_by(spId, metric, ssb_true, nPres,  prevCat, testDat, ssb_meth, blkType) %>%
        slice(which.min(m))
    } else {
      best_int <- best_int %>%
        group_by(spId, metric, ssb_true, nPres,  prevCat, testDat, ssb_meth, blkType) %>%
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
                    "prevCat",
                    "prevelance",
                    "metric",
                    "ssb_true", 
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
        dplyr::select(spId, metric, nPres, ssb_true, ssb_meth,
                      blkType, testDat, d_m, d_sd)
    }
    dat_comp
  }
