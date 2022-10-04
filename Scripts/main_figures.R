#' Prepared data for evaluation error vs predictive performance plots
#' 
#' @param dat_in Model predictive performance.
#' @param true_ssb Value of true SSB ("wgt0", "wgt06", "wgt08", "wgt1").
#' @param mod_pp Selecting the best, worse or random model.
#' 
#' @return ggplot object
#' 
#' @export
err_vs_pp_dat <-
  function(dat_in,
           true_ssb,
           mod_pp = c("best", "worse", "random")) {
    #--- Prepare best data
    d <-
      process.comparision.data(
        dat_in,
        mod_pp,
        n_pres = c("50", "200"),
        metric = c("AUC",  "Boyce", "S (Pearson's)", "S (Spearman's)"),
        ssbTrue = true_ssb,
        FALSE
      )
    
    #-- Mean dif. and statistical sig. between internal and independent data
    d_pp_cv <- d %>%
      filter(ssb_correct == "Yes") %>%
      pivot_wider(names_from = testDat,
                  values_from = c(m, sd, N)) %>%
      mutate(
        mdif = (m_int - m_ind),
        A = sd_int ^ 2 / N_int,
        B = sd_ind ^ 2 / N_ind,
        ds = sqrt(A + B),
        tt =  mdif / ds,
        df = (A + B) ^ 2 / (A ^ 2 / (N_int - 1) + B ^ 2 / (N_ind - 1)),
        p = 1 - pt(tt, df),
        pN =  N_int + N_ind
      ) %>%
      dplyr::select(!c(A, B, sd_int,  sd_ind)) %>%
      mutate(p_sig = ifelse(p < 0.05, 1, 0)) %>%
      group_by(metric, nPres, ssb_meth, blkType) %>%
      summarise(
        gm1 = weighted.mean(mdif, pN, na.rm = TRUE),
        # Grand mean weighted
        gsd1 = sqrt(
          weighted.mean(ds ^ 2 + mdif ^ 2, pN, na.rm = TRUE) -
            weighted.mean(mdif, pN, na.rm = TRUE) ^ 2
        ),
        # Grand sd weighted
        RMSE = rmse(m_ind, m_int),
        n_sp = n(),
        n_sig = sum(p_sig, na.rm = TRUE)
      ) %>%
      mutate(pSig_intvInd = n_sig / n_sp)
    
    #-- Model improvement (delta predictive performance) corrected vs uncorrected
    d_pp_m <-
      d %>%
      pivot_wider(names_from = ssb_correct,
                  values_from = c(m, sd, N)) %>%
      mutate(
        mdif = (m_Yes - m_No),
        A = sd_Yes ^ 2 / N_Yes,
        B = sd_No ^ 2 / N_No,
        ds = sqrt(A + B),
        tt =  mdif / ds,
        df = (A + B) ^ 2 / (A ^ 2 / (N_Yes - 1) + B ^ 2 / (N_No - 1)),
        p = 1 - pt(tt, df),
        pN =  N_Yes + N_No
      ) %>%
      dplyr::select(!c(A, B, df)) %>%
      mutate(p_sig = ifelse(p < 0.05, 1, 0),
             impr = ifelse(mdif > 0, 1, 0)) %>%
      group_by(metric, nPres, ssb_meth, blkType, testDat) %>%
      summarise(
        n_sp = n(),
        n_sig = sum(p_sig, na.rm = TRUE),
        propImp = sum(impr, na.rm = TRUE) / n_sp,
        gm2 = weighted.mean(mdif, pN, na.rm = TRUE),
        # Grand mean weighted
        gsd2 = sqrt(
          weighted.mean(ds ^ 2 + mdif ^ 2, pN, na.rm = TRUE) -
            weighted.mean(mdif, pN, na.rm = TRUE) ^ 2
        ),
        c95L_fx = quantile(mdif, prob = 0.025, na.rm = TRUE),
        c95U_fx = quantile(mdif, prob = 0.975, na.rm = TRUE),
      ) %>%
      mutate(pSig_Adj = n_sig / n_sp) %>%
      filter(testDat == "ind") %>%
      dplyr::select(!c(n_sp, n_sig))
    
    #-- Figure RMSE vs effect size
    pp_comp <- left_join(d_pp_cv, d_pp_m) %>%
      mutate(
        ssb_meth =
          fct_recode(
            ssb_meth,
            "Spp. target grp" = "spTgGrp",
            'Feat. distance decay' = "featDD",
            'Occurrence buffer' =    "occBuff",
            'Occurrence filter' =    "occFilter",
            'Covariate (cond. pred)' =   "covCond"
          ),
        blkType =
          fct_recode(
            blkType,
            'Random' = "random",
            'Spatial systematic' = "spatial_sytematic",
            'Environmental' = "enviro"
          ),
        nPres =
          fct_recode(nPres,
                     'N = 50' = "50",
                     'N = 200' = "200"),
      )
    pp_comp
  }

#' Create evaluation error vs predictive performance plots
#' 
#' @param pp_comp Output from err_vs_pp_dat. 
#' @param metricA Metric for top panel
#' @param metricB Metric for lower panel
#' 
#' @return ggplot object
#' 
#' @export
err_vs_pp_plot <- function(pp_comp, 
                           metricA = c("AUC", "Boyce"), 
                           metricB = c("S (Spearman's)","S (Pearson's)")) {
  
  a_lab = metricA
  b_lab = metricB
  if(b_lab == "S (Spearman's)") {
    b_lab <- "s"
  } else {
    b_lab <- "p"
  }
  
  figa <-  pp_comp %>%
    filter(metric == {{metricA}}) %>%
    ggplot() +
    geom_hline(
      yintercept = 0,
      size = 0.5,
      colour = "grey50",
      linetype = 2
    ) +
    geom_vline(
      xintercept = 0,
      size = 0.5,
      colour = "grey50",
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        x = gm1,
        ymin = gm2 - gsd2,
        ymax = gm2 + gsd2,
        colour = ssb_meth
      ),
      size = 0.5,
      width = 0.045
    ) +
    geom_errorbarh(
      aes(
        y = gm2,
        xmin = gm1 - gsd1,
        xmax = gm1 + gsd1,
        colour = ssb_meth
      ),
      #alpha = 0.5,
      size = 0.5,
      height = 0.01
    ) +
    geom_point(aes(x = gm1,
                   y = gm2    ,
                   #size = pSig_intvInd,
                   colour = ssb_meth),
               alpha = 1,
               size = 1.75) +
    facet_grid(nPres ~ blkType, scales = "fixed") +
    coord_flip() +
    theme_bw(base_size = 12) %+replace%
    theme(strip.background = element_rect(fill = "white")) +
    scale_color_tableau(name = "SSB correction\nmethod") +
    scale_size(range = c(1, 4),
               name = "Proportion \ntests sig. dif.\nbetween test\ndatasets") +
    xlab(bquote(atop(Delta*.(a_lab) ~ " (\U03BC \u00b11SD)", "Indep. vs. Intern. cv"))) +
    ylab(bquote(Delta*.(a_lab) ~ " (\U03BC \u00b11SD) with SSB adjustment"))
  
  figb <- pp_comp %>%
    filter(metric == {{metricB}}) %>% 
    ggplot() +
    geom_hline(
      yintercept = 0,
      size = 0.5,
      colour = "grey50",
      linetype = 2
    ) +
    geom_vline(
      xintercept = 0,
      size = 0.5,
      colour = "grey50",
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        x = gm1,
        ymin = gm2 - gsd2,
        ymax = gm2 + gsd2,
        colour = ssb_meth
      ),
      size = 0.5,
      width = 0.045
    ) +
    geom_errorbarh(
      aes(
        y = gm2,
        xmin = gm1 - gsd1,
        xmax = gm1 + gsd1,
        colour = ssb_meth
      ),
      #alpha = 0.5,
      size = 0.5,
      height = 0.01
    ) +
    geom_point(aes(x = gm1,
                   y = gm2    ,
                   #size = pSig_intvInd,
                   colour = ssb_meth),
               alpha = 1,
               size = 1.75) +
    facet_grid(nPres ~ blkType, scales = "fixed") +
    coord_flip() +
    theme_bw(base_size = 12) %+replace%
    theme(strip.background = element_rect(fill = "white")) +
    scale_color_tableau(name = "SSB correction\nmethod") +
    scale_size(range = c(1, 4),
               name = "Proportion \ntests sig. dif.\nbetween test\ndatasets") +
    xlab(bquote(atop(Delta*"r"[.(b_lab)] ~ " (\U03BC \u00b11SD)", 
                     "Indep. vs. Intern. cv"))) +
    ylab(bquote(Delta*"r"[.(b_lab)] ~ " (\U03BC \u00b11SD) with SSB adjustment"))
  
  #- Panel figure with legend
  fig <- figa / figb + 
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
    theme(legend.position = "bottom", legend.title = element_blank())
  
}


#' Cross-validation model improvement
#' 
#' 
cv_mod_improve <- function(dat_in, 
                           true_ssb, 
                           metricA = c("AUC", "Boyce"), 
                           metricB = c("S (Spearman's)","S (Pearson's)")) {
  
  a_lab = metricA
  b_lab = metricB
  if(b_lab == "S (Spearman's)") {
    b_lab <- "s"
  } else {
    b_lab <- "p"
  }
  
  #- Find best/worse model
  d_best <-
    process.comparision.data(
      dat_in,
      "best",
      n_pres = c("50", "200"),
      metric = c("AUC",  "Boyce", "S (Pearson's)", "S (Spearman's)"),
      ssbTrue = true_ssb,
      TRUE
    ) %>%
    mutate(cv_qual = "Best model identified using internal cv")
  d_worse <-
    process.comparision.data(
      dat_in,
      "worse",
      n_pres = c("50", "200"),
      metric = c("AUC",  "Boyce", "S (Pearson's)", "S (Spearman's)"),
      ssbTrue = true_ssb,
      TRUE
    ) %>%
    mutate(cv_qual = "Worse model identified using internal cv")
  d <- rbind(d_best, d_worse)
  
  #- Calculate delta between int and ind evaluation
  d <- d %>%
    filter(testDat == "ind") %>%
    group_by(metric, nPres, ssb_meth, blkType, cv_qual) %>%
    #- Calculate mean and 95% quantiles across treatments
    summarise(
      m = mean(d_m, na.rm = TRUE),
      l95 = quantile(d_m, prob = 0.025, na.rm = TRUE),
      u95 = quantile(d_m, prob = 0.975, na.rm = TRUE)
    ) %>%
    mutate(
      ssb_meth =
        fct_recode(
          ssb_meth,
          "Spp. target grp" = "spTgGrp",
          'Feat. distance decay' = "featDD",
          'Occurrence buffer' =    "occBuff",
          'Occurrence filter' =    "occFilter",
          'Covariate (cond. pred)' =   "covCond"
        ),
      blkType =
        fct_recode(
          blkType,
          'Random' = "random",
          'Spatial systematic' = "spatial_sytematic",
          'Environmental' = "enviro"
        ),
      nPres =
        fct_recode(nPres,
                   'N = 50' = "50",
                   'N = 200' = "200"),
    )
  
  #- Plot
  figa <- d %>% 
    filter(metric == 'AUC') %>%
    ggplot() +
    geom_point(aes(x = ssb_meth, 
                   y = m,
                   colour = cv_qual,
                   shape = cv_qual),
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = ssb_meth, 
                      ymax = u95, 
                      ymin = l95,
                      colour = cv_qual,
                      linetype = cv_qual),
                  position = position_dodge(width = 0.5),
                  width = 0.2) +
    facet_grid(nPres ~ blkType) +
    theme_bw(base_size = 12) %+replace%
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) +
    scale_color_tableau(name = "SSB correction\nmethod") +
    scale_shape(name = "SSB correction\nmethod") +
    scale_linetype(name = "SSB correction\nmethod") +
    ylab(bquote(Delta*.(a_lab) ~ " (\U03BC \u00b1 95% quantiles)"))
  
  figb <- d %>% 
    filter(metric == "S (Spearman's)") %>%
    ggplot() +
    geom_point(aes(x = ssb_meth, 
                   y = m,
                   colour = cv_qual,
                   shape = cv_qual),
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = ssb_meth, 
                      ymax = u95, 
                      ymin = l95,
                      colour = cv_qual,
                      linetype = cv_qual),
                  position = position_dodge(width = 0.5),
                  width = 0.2) +
    facet_grid(nPres ~ blkType) +
    theme_bw(base_size = 12) %+replace%
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 45, 
                                     hjust = 1, 
                                     vjust = 1),
          axis.title.x = element_blank()) +
    scale_color_tableau(name = "SSB correction\nmethod") +
    scale_shape(name = "SSB correction\nmethod") +
    scale_linetype(name = "SSB correction\nmethod") +
    ylab(bquote(Delta*"r"[.(b_lab)] ~ " (\U03BC \u00b1 95% quantiles)"))
    
  #- Panel figure with legend
  fig <- figa / figb + 
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
    theme(legend.position = "bottom", legend.title = element_blank())
  
}




