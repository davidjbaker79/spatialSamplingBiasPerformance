#' Prepare statistics for agreement vs predictive performance
#' 
#' Calculate performance statistics for specified subset of simulation 
#' data, including differences between models corrected and uncorrected for SSB
#' and differences in performance evaluation between internal and independent 
#' test data.
#' 
#' @param dat_in Model predictive performance.
#' @param true_ssb Value of true SSB ("wgt0", "wgt06", "wgt08", "wgt1").
#' 
#' @import ggplot2
#' @import ddplyr
#' 
#' @return data.frame object
#' 
#' @export
ag_vs_pp_dat <- function(dat_in, true_ssb) {

    dat_in <- dat_in %>% 
      filter(ssb_true == true_ssb)
    
    #-- Compare internal vs independent testing 
    d_pp_cv <-
      dat_in %>%
      pivot_wider(names_from = c(ssb_correct, testDat),
                  values_from = c(m, sd, N)) %>%
      dplyr::select(
        "spId",
        "nPres",
        "prevCat",
        "metric",
        "ssb_meth",
        "blkType",
        "m_Yes_int",
        "m_Yes_ind",
        "m_No_int",
        "m_No_ind"
      ) %>%
      mutate(mdif_int = (m_Yes_int - m_No_int),
             mdif_ind = (m_Yes_ind - m_No_ind),
             intIndAgg = ifelse(sign(mdif_int) != sign(mdif_ind), 0, 1),
             propDecr = ifelse(sign(mdif_ind) < 0, 1, 0),
             mdif_ind_Decr = mdif_ind * propDecr,
             propIncr = ifelse(sign(mdif_ind) > 0, 1, 0),
             mdif_ind_Incr = mdif_ind * propIncr) %>%
      na.omit() %>% 
      group_by(metric, nPres, prevCat, blkType, ssb_meth) %>%
      summarise(cor = cor(mdif_int, mdif_ind),
                intIndAgg = sum(intIndAgg)/n(),
                propDecr = sum(propDecr)/n(),
                propIncr = sum(propIncr)/n(),
                mdif_ind_m = mean(mdif_ind, na.rm = T),
                mdif_ind_Decr = mean(mdif_ind_Decr, na.rm = T),
                mdif_ind_Incr = mean(mdif_ind_Incr, na.rm = T))
    
    #-- Model improvement (delta predictive performance) corrected vs uncorrected
    # for both internal and independent.
    d_pp_m <-
      dat_in %>%
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
      group_by(metric, nPres, prevCat, ssb_meth, blkType, testDat) %>%
      summarise(
        n_sp = n(),
        n_sig = sum(p_sig, na.rm = TRUE),
        #propImp = sum(impr, na.rm = TRUE) / n_sp,
        # Grand mean weighted
        gm = weighted.mean(mdif, pN, na.rm = TRUE),
        # Grand sd weighted
        gsd = sqrt(
          weighted.mean(ds ^ 2 + mdif ^ 2, pN, na.rm = TRUE) -
            weighted.mean(mdif, pN, na.rm = TRUE) ^ 2
        ),
        c95L_fx = quantile(mdif, prob = 0.025, na.rm = TRUE),
        c95U_fx = quantile(mdif, prob = 0.975, na.rm = TRUE),
        c50L_fx = quantile(mdif, prob = 0.25, na.rm = TRUE),
        c50U_fx = quantile(mdif, prob = 0.75, na.rm = TRUE),
      ) %>%
      mutate(pSig_Adj = n_sig / n_sp) %>%
      filter(testDat == "ind") %>%
      dplyr::select(!c( n_sig))
    
    #-- Figure cv vs effect size
    pp_comp <- left_join(d_pp_cv, d_pp_m) %>%
      mutate(
        ssb_meth =
          fct_recode(
            ssb_meth,
            "AdjBkGrd(SppTgGrp)" = "spTgGrp",
            'AdjBkGrd(Feature)' = "featDD",
            'AdjBkGrd(Occurrence)' = "occBuff",
            'OccFilter(Geo)' = "occFilter",
            'Covariate(CondPred)' =  "covCond"
          ),
        blkType =
          fct_recode(
            blkType,
            'Random' = "random",
            'Spatial systematic' = "spatial_sytematic",
            'Environmental' = "enviro"
          ),
        nPres =
          fct_recode(
            nPres,
            'N = 50' = "50",
            'N = 100' = "100",
            'N = 200' = "200",
            'N = 400' = "400",
            'N = 800' = "800"
          ),
        prevCat = fct_recode(
          prevCat,
          'Prev. = 10%' = "0.1",
          'Prev. = 20%' = "0.2",
          'Prev. = 30%' = "0.3",
          'Prev. = 40%' = "0.4",
          'Prev. = 50%' = "0.5",
          'Prev. = 60%' = "0.6"
        )
      )
    pp_comp
  }

#' Create figure for comparison of agreement vs predictive performance
#' 
#' Takes outputs from \code{ag_vs_pp_dat} and creates facet plot of showing 
#' agreement vs predictive performance for simulation results.
#' 
#' @param pp_comp Output from \code{ag_vs_pp_dat} 
#' @param metricA A character vector giving the metric for plotting
#' @param N A character vector with required sample sizes.
#' @param Prev A character vector with required prevalence category.
#' 
#' @import ggplot2
#' @import ddplyr
#' 
#' @return list including ggplot object and summary data.frame.
#' 
#' @export
ag_vs_pp_plot <- function(pp_comp,
                           metricA = c("AUC", "Boyce", "S (Spearman's)", "S (Pearson's)"),
                           N = c("N = 100", "N = 800"),
                           Prev = c("Prev. = 10%", "Prev. = 40%")) {

  pp_comp <- pp_comp %>%  
    filter(metric == {
      {
        metricA
      }
    })  
  pp_comp_p <- pp_comp
  pp_comp_p <- pp_comp_p[pp_comp_p$nPres %in% N,]
  pp_comp_p <- pp_comp_p[pp_comp_p$prevCat %in% Prev,]
  pp_comp_p <- pp_comp_p %>% droplevels()
  
  figa <- pp_comp_p %>%
    filter(metric == {
      {
        metricA
      }
    }) %>%
    ggplot() +
    geom_vline(
      xintercept = 0,
      linewidth = 0.5,
      colour = "grey50",
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        y = intIndAgg,
        xmin = gm - gsd,
        xmax = gm + gsd,
        colour = ssb_meth
      ),
      linewidth = 0.5,
      width = 0.045
    ) +
    geom_point(aes(y = intIndAgg,
                   x = gm,
                   colour = ssb_meth),
               alpha = 1,
               size = 1.75) +
    facet_grid(nPres + prevCat ~ blkType, scales = "fixed") +
    ylim(0, 1.05) +
    coord_flip() +
    theme_bw(base_size = 12) %+replace%
    theme(strip.background = element_rect(fill = "white")) +
    scale_color_tableau(name = "SSB correction\nmethod")
  
  if (metricA == "AUC") {
    figa <- figa +
      ylab(bquote(
        atop(
          "Proportion agreement in direction of effect between",
          "internal vs independent evaluation of SSB (" * Delta * .("AUC") * ")"
        )
      )) +
      xlab(bquote(Delta * .("AUC") ~ " (\U03BC \u00b1 SD) with SSB adjustment"))
  }
  
  if (metricA == "Boyce") {
    figa <- figa +
      ylab(bquote(
        atop(
          "Proportion agreement in direction of effect between",
          "internal vs independent evaluation of SSB (" * Delta * .("Boyce") * ")"
        )
      )) +
      xlab(bquote(Delta * .("Boyce") ~ " (\U03BC \u00b1 SD) with SSB adjustment"))
  }
  
  if (metricA == "S (Spearman's)") {
    figa <- figa +
      ylab(bquote(
        atop(
          "Proportion agreement in direction of effect between",
          "internal vs independent evaluation of SSB (" * Delta * "r"[.("s")] * ")"
        )
      )) +
      xlab(bquote(Delta * "r"[.("s")] ~ " (\U03BC \u00b1 SD) with SSB adjustment"))
  }
  
  if (metricA == "S (Pearson's)") {
    figa <- figa +
      ylab(bquote(
        atop(
          "Proportion agreement in direction of effect between",
          "internal vs independent evaluation of SSB (" * Delta * "r"[.("p")] * ")"
        )
      )) +
      xlab(bquote(Delta * "r"[.("p")] ~ " (\U03BC \u00b1 SD) with SSB adjustment"))
  }
  
  list(figa, pp_comp)
  
}

#' Wrapper for ag_vs_pp_plot (fig. 2)
#' 
#' For running combinations of wgt and metric efficiently.
#' 
#' @param x Output from ag_vs_pp_dat. 
#' @param wgt Metric for top panel
#' @param metricA Metric for lower panel
#' @param corNiche
#' 
#' @return ggplot object
#' 
#' @export
create.fig2.wrapper <- function(x, wgt, metricA, corNiche) {
  
  pp_comp <- ag_vs_pp_dat(x, wgt)
  ssb_fig <- ag_vs_pp_plot(pp_comp, metricA)
  ssb_fig[[2]]$ssb_true <- wgt
  ssb_fig
}

#' Panel wrapper for fig. 2
#' 
#' Creates panels for fig. 2 and supplementary figures.
#' 
#' @param p1 ggplot for plot1
#' @param p2 ggplot for plot2
#' @param wgt A character string for labelling file
#' @param metric A character string for labelling file
#' 
#' @import patchwork
#' 
#' @return ggplot object
#' 
#' @export
panel_fig_wrapper <- function(p1, p2, corType, metric) {
  pp <- p1 / p2  +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'a', tag_suffix = ')') &
    theme(legend.position = "bottom", legend.title = element_blank())
  ggsave(
    filename = paste0("Outputs/Figures_R2/Figure_2_R2_", corType, "_", metric, ".png"),
    plot = pp,
    width = 10,
    height = 12
  )
  p1
}
  
#' Proportion increasing / decreasing by block and ssb strength
#' 
#' Extract statistics for proportion of test statistics increasing / decreasing
#' with SSB correction, as indicated by evaluation on independent test data.
#' 
#' @param dat A data.frame from \Code{ag_vs_pp_plot}.
#' 
#' @import dplyr (edit function for piping)
#' 
#' @return data.frame
#' 
#' @export
prop_inc_dec <- function(dat) {
  ssb_inc_dec_tab <- dat %>%
    mutate(
      propIncr = propIncr * 100,
      propDecr = propDecr * 100,
      blkType = sub(" ", "_", blkType)
    ) %>%
    group_by(ssb_meth, ssb_true, blkType) %>%
    summarise(
      propInc_m = round(mean(propIncr), 1),
      mdif_inc = round(mean(mdif_ind_Incr, na.rm = T), 3),
      propDec_m = round(mean(propDecr), 1),
      mdif_dec = round(mean(mdif_ind_Decr, na.rm = T), 3)
    ) %>%
    ungroup() %>%
    pivot_wider(
      names_from = c(blkType),
      values_from = c(propInc_m, mdif_inc, propDec_m, mdif_dec)
    )  %>%
    mutate(
      propIncRand = paste0(propInc_m_Random, "% (", mdif_inc_Random, ")"),
      propIncSyst = paste0(
        propInc_m_Spatial_systematic,
        "% (",
        mdif_inc_Spatial_systematic,
        ")"
      ),
      propIncEnv = paste0(propInc_m_Environmental, "% (", mdif_inc_Environmental, ")"),
      propDecRand = paste0(propDec_m_Random, "% (", mdif_dec_Random, ")"),
      propDecSyst = paste0(
        propDec_m_Spatial_systematic,
        "% (",
        mdif_dec_Spatial_systematic,
        ")"
      ),
      propDecEnv = paste0(propDec_m_Environmental, "% (", mdif_dec_Environmental, ")")
    ) %>%
    dplyr::select(
      ssb_meth,
      ssb_true,
      propIncRand,
      propIncSyst,
      propIncEnv,
      propDecRand,
      propDecSyst,
      propDecEnv
    )
  
}

#' Proportion agreement between internal and independent evaluation
#' 
#' Extract statistics for proportion agreement between internal and independent
#' evaluation on the direction of effect of SSB correction.
#' 
#' @param dat A data.frame from \Code{ag_vs_pp_plot}.
#' 
#' @import dplyr (edit function for piping)
#' 
#' @return data.frame
#' 
#' @export
prop_agree_indInt <- function(dat) {
  dat %>%
    mutate(intIndAgg = intIndAgg * 100,
           blkType = sub(" ", "_", blkType)) %>%
    group_by(ssb_meth, ssb_true, blkType) %>%
    summarise(intIndAgg = round(mean(intIndAgg), 1)) %>%
    ungroup() %>%
    pivot_wider(names_from = c(blkType),
                values_from = c(intIndAgg)) %>%
    dplyr::select(ssb_meth, 
                  ssb_true, 
                  Random, 
                  Spatial_systematic, 
                  Environmental)
    
}

#' Calculate and plot partial effects
#' 
#' Calculate and plot partial effects from models relating Δ[metric] to SSB 
#' strength (none, weak, moderate, strong), SSB correction method, number of 
#' presences using to build SDMs, and the species prevalence, and each two-way 
#' interactions between the correction method, prevalence and number of 
#' presences. Species ID was included as a random intercept. Models were fit 
#' using robust linear mixed-effects models, using the robustlmm package 
#' (Koller 2016), to account for the presence of extreme values in estimates 
#' of the mean effects.
#' 
#' @param mod_i A statistical model output.
#' 
#' @import ggeffects 
#' @import ggplot2
#' @import patchwork
#' 
#' @return data.frame
#' 
#' @export
plot_partial_effects <- function(mod_i) {

  # Predicted effect
  g1 <-
    ggpredict(mod_i,
              terms = c("nPres", "prevelance", "ssb_meth", "ssb_true"))
  
  # Create plot
  p <- plot(g1)
  
  # New labels
  hp_labeller <- as_labeller(
    c(
      "featDD" = 'AdjBkGrd(Feature)' ,
      "occBuff" = 'AdjBkGrd(Occurrence)',
      "occFilter" = 'OccFilter(Geo)',
      "covCond" = 'Covariate(CondPred)'
    )
  )
  
  # Create plot for each panel
  pZero <- p[[1]] +
    facet_grid(cols = vars(facet),
               labeller = labeller(facet = hp_labeller)) +
    scale_fill_colorblind(name = "Species\nprevalence") +
    scale_colour_colorblind(name = "Species\nprevalence") +
    theme_bw(base_size = 12) %+replace% theme(
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom"
    ) +
    ggtitle("A) SSB strength: None") +
    xlab("Number of presences") +
    ylab(bquote("Predicted " * Delta [.("[metric]")]))
  
  pWeak <- p[[2]] +
    facet_grid(cols = vars(facet),
               labeller = labeller(facet = hp_labeller)) +
    scale_fill_colorblind(name = "Species\nprevalence") +
    scale_colour_colorblind(name = "Species\nprevalence") +
    theme_bw(base_size = 12) %+replace% theme(
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom"
    ) +
    ggtitle("B) SSB strength: Weak") +
    xlab("Number of presences") +
    ylab(bquote("Predicted " * Delta [.("[metric]")]))
  
  pMode <- p[[3]] +
    facet_grid(cols = vars(facet),
               labeller = labeller(facet = hp_labeller)) +
    scale_fill_colorblind(name = "Species\nprevalence") +
    scale_colour_colorblind(name = "Species\nprevalence") +
    theme_bw(base_size = 12) %+replace% theme(
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom"
    ) +
    ggtitle("A) SSB strength: Moderate") +
    xlab("Number of presences") +
    ylab(bquote("Predicted " * Delta [.("[metric]")]))
  
  pStrg <- p[[4]] +
    facet_grid(cols = vars(facet),
               labeller = labeller(facet = hp_labeller)) +
    scale_fill_colorblind(name = "Species\nprevalence") +
    scale_colour_colorblind(name = "Species\nprevalence") +
    theme_bw(base_size = 12) %+replace% theme(
      plot.subtitle = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom"
    ) +
    ggtitle("A) SSB strength: Strong") +
    xlab("Number of presences") +
    ylab(bquote("Predicted " * Delta [.("[metric]")]))
  
  # Create panel
  pout <- (pZero / pWeak / pMode / pStrg)  +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  pout
  
}
