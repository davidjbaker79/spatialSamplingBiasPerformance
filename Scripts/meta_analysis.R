#------------------------------------------------------------------------------#
#
#     Spatial sampling bias systematic review and meta analysis
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Load literature extraction data
#------------------------------------------------------------------------------#

#' The systematic review was conducted using Web of Science and Scopus databases
#' on the 22/05/2022 using the following search string:
#'
#' ALL=(("species distribution*" OR SDM OR "environmental niche" OR ENM OR
#' "resource selection" OR "habitat selection" OR suitability OR occurrence)
#' AND ("presence-only" OR “presence data” OR "presence-background" OR
#' “pseudo absence” OR opportunistic OR “citizen science” OR “preferential”) AND
#' (spatial* OR geographic* OR bias* OR correct* OR adjust* OR account
#' OR control))

#--- Load data for full systematic review search results and filtering
sr_search_dat <-
  read_excel("Data/baker_sdm_ssb_literature_review_22_05_22.xlsx",
             sheet = "literatureFilter")

#--- Load data for data extraction from included studies
sr_dat <-
  read_excel(
    "Data/baker_sdm_ssb_literature_review_22_05_22.xlsx",
    sheet = "dataExtraction",
    na = "NA"
  )
names(sr_dat)
#---
#-- Classify bias correction methods into broad classes
#---

#-#-# Covert to factors and level
sr_dat <- sr_dat %>%
  mutate(
    bias_correction_method = factor(
      bias_correction_method,
      levels =
        c(
          "Species_grp_records",
          "Feature_bufferOrDD",
          "Effort_bufferOrDD",
          "Occurrence_bufferOrDD",
          "Occurrence_filter",
          "Covariate",
          "Covariate_condPred",
          "Model_other",
          "Occurrence_filter + Species_grp_records",
          "Occurrence_filter + Feature_bufferOrDD",
          "Occurrence_filter + Occurrence_bufferOrDD",
          "Occurrence_filter + Covariate",
          "Effort_bufferOrDD + Feature_bufferOrDD",
          "Occurrence_filter + Feature_bufferOrDD + Species_grp_records",
          "Unclear"
        )
    ),
    bias_correction_method = fct_recode(
      bias_correction_method,
      "Multiple" = "Occurrence_filter + Feature_bufferOrDD",
      "Multiple" = "Occurrence_filter + Occurrence_bufferOrDD",
      "Multiple" = "Occurrence_filter + Covariate",
      "Multiple" = "Effort_bufferOrDD + Feature_bufferOrDD",
      "Multiple" = "Occurrence_filter + Feature_bufferOrDD + Species_grp_records",
      "Multiple" = "Occurrence_filter + Species_grp_records"
    )
  )

#-------------------------------------------------------------------------------
# Summary statistics
#-------------------------------------------------------------------------------

#---
#-- Basic summary statistics for systematic review filtering process
#---

#--- Number of studies in initial search (non duplicated)
nrow(sr_search_dat) # 4150

#--- Number of studies filtered by title
sr_search_dat %>%
  filter(select_title == "Y") %>%
  nrow() # 1548

#--- Number of studies filtered by abstract
sr_search_dat %>%
  filter(select_abstract == "Y") %>%
  nrow() # 1053

#--- Number of studies filtered by content
sr_search_dat %>%
  filter(select_content == "Y") %>%
  nrow() # 220

#-------------------------------------------------------------------------------
# Summarise the frequency SSB correction method usage
#-------------------------------------------------------------------------------

#--- Summarise number of studies by method type
ssb_mthd_div <- sr_dat %>%
  dplyr::select(study_id, bias_correction_method) %>% # Class2 is the broad classes
  distinct() %>%
  group_by(bias_correction_method) %>%
  tally() %>%
  arrange(desc(n)) %>%
  mutate(bias_correction_method = fct_inorder(bias_correction_method)) %>%
  mutate(pc = n / sum(n)) %>%
  filter(bias_correction_method != "NA")  %>%
  mutate(
    bias_correction_method =
      fct_recode(
        bias_correction_method,
        "Spp. target grp" = "Species_grp_records",
        "Covariate (cond. pred.)" =  "Covariate_condPred",
        "Feature (buffer or dist. decay)" = "Feature_bufferOrDD",
        "Occurrence (buffer or dist. decay)" = "Occurrence_bufferOrDD",
        "Occurrence filter" = "Occurrence_filter",
        "Effort (buffer or dist. decay)" = "Effort_bufferOrDD",
        "Model-based (other)" = "Model_other"
      )
  )

#--- Create frequency histogram of SSB frequency
ssb_mthd_div_p <- ssb_mthd_div %>%
  ggplot() +
  geom_col(aes(x = bias_correction_method, y = n), alpha = 0.75) +
  theme_bw(base_size = 14) %+replace%
  theme(axis.text.x = element_text(
    angle = 55,
    vjust = 1,
    hjust = 1
  ),
  axis.title.x = element_blank()) +
  ylab("Frequency of use")

#--- Save figure 1
ggsave(
  filename = "Outputs/Figures/Figure_1_methodFreq.png",
  plot = ssb_mthd_div_p,
  width = 10,
  height = 6
)

#-------------------------------------------------------------------------------
# Meta-analysis: the effect of SSB correction on SDM performance
#-------------------------------------------------------------------------------

#---
#-- Prepare model data
#---

#-- Subset data to just those with comparisons between naive and bias corrected
biasCorrEval <- sr_dat %>%
  filter(comparison_to_naive == "Yes",
         useInMa == 1) %>%
  filter(test_dataset != "Simulated") %>%
  mutate(
    focal_m = as.numeric(focal_m),
    focal_sd = as.numeric(focal_sd),
    focal_n = as.numeric(focal_n),
    focal_lci = as.numeric(focal_lci),
    focal_uci = as.numeric(focal_uci),
    naïve_m = as.numeric(naïve_m),
    naïve_sd = as.numeric(naïve_sd),
    naïve_n = as.numeric(naïve_n),
    naïve_lci = as.numeric(naïve_lci),
    naïve_uci = as.numeric(naïve_uci)
  ) %>%
  mutate(
    test_dataset =
      recode(test_dataset,
             Independent = "Test data = Independent",
             Internal = "Test data = Internal")
  )

# Number of studies with data
length(unique(biasCorrEval$study_id))

# Number of independent
length(unique(biasCorrEval$study_id[biasCorrEval$test_dataset == "Test data = Independent"]))

# Number of internal
length(unique(biasCorrEval$study_id[biasCorrEval$test_dataset == "Test data = Internal"]))

# Number using AUC
length(unique(biasCorrEval$study_id[biasCorrEval$performance_metric == "AUC"]))

#---
#-- Calculate standardised mean differences
#---

#-- Set the columns by which to group test statistics
grpByCols <- c(
  "study_id",
  "bias_correction_method",
  "bias_correction_treatment",
  "sdm_method",
  "performance_metric",
  "test_dataset",
  "useInMa"
)

#--- Individual reported results from group assessments
#' Some results were given for individual models but want to aggregate these
#' to study / treatment level statistics for ease of comparison. Some were only
#' given as a single value (no replication or no variation reported). These will
#' be removed for effect size calculations, but leave in for now as we'll plot
#' then unstandardised first.

#-- Summarise at group_by levels
biasCorrEval_values <- biasCorrEval %>%
  filter(quant_metric == "value")
naïve_m <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(naïve_m = mean(naïve_m, na.rm = TRUE))
naïve_sd <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(naïve_sd = sd(naïve_m, na.rm = TRUE))
naïve_n <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(naïve_n = length(naïve_m))
focal_m <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(focal_m = mean(focal_m, na.rm = TRUE))
focal_sd <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(focal_sd = sd(focal_m, na.rm = TRUE))
focal_n <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(focal_n = length(focal_m))

#--- Meta data
metaDat_values <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(
    grain = mean(as.numeric(grain), na.rm = TRUE),
    min_number_rec = mean(as.numeric(min_number_rec), na.rm = TRUE),
    max_number_rec = mean(as.numeric(max_number_rec), na.rm = TRUE),
    focal_kingdom = paste(unique(focal_kingdom), collapse = "|")
  )

#--- Combine
biasCorrEval_values <- Reduce(
  function(...)
    left_join(...),
  list(focal_m, focal_sd, focal_n,
       naïve_m, naïve_sd, naïve_n)
)

#--- For results given as Quantile + Q1 and Q3 (i.e. from boxplots)
biasCorrEval_iqr2msd <- biasCorrEval %>%
  filter(quant_metric == "median_iqr") %>%
  rename(
    naïve_q1 = "naïve_lci",
    naïve_q3 = "naïve_uci",
    focal_q1 = "focal_lci",
    focal_q3 = "focal_uci"
  ) %>%
  mutate(
    naïve_m =  (naïve_q1 + naïve_m + naïve_q3) / 3,
    naïve_sd = (naïve_q3 - naïve_q1) / (2 * qnorm((0.75 * naïve_n - 0.125) / (naïve_n + 0.25), 0, 1
    )),
    focal_m =  (focal_q1 + focal_m + focal_q3) / 3,
    focal_sd = (focal_q3 - focal_q1) / (2 * qnorm((0.75 * focal_n - 0.125) / (focal_n + 0.25), 0, 1
    ))
  ) %>%
  dplyr::select(names(biasCorrEval_values))

metaDat_iqr2msd <- biasCorrEval %>%
  filter(quant_metric == "median_iqr") %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(
    grain = mean(as.numeric(grain), na.rm = TRUE),
    min_number_rec = mean(as.numeric(min_number_rec), na.rm = TRUE),
    max_number_rec = mean(as.numeric(max_number_rec), na.rm = TRUE),
    focal_kingdom = paste(unique(focal_kingdom), collapse = "|")
  )

#--- For results given as mean + sd
biasCorrEval_msd <- biasCorrEval %>%
  filter(quant_metric == "mean_sd") %>%
  dplyr::select(names(biasCorrEval_values))

metaDat_msd <- biasCorrEval %>%
  filter(quant_metric == "mean_sd") %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(
    grain = mean(as.numeric(grain), na.rm = TRUE),
    min_number_rec = mean(as.numeric(min_number_rec), na.rm = TRUE),
    max_number_rec = mean(as.numeric(max_number_rec), na.rm = TRUE),
    focal_kingdom = paste(unique(focal_kingdom), collapse = "|")
  )

#--- Combine
biasCorrEval_msd <-
  do.call(rbind,
          list(biasCorrEval_values, biasCorrEval_iqr2msd, biasCorrEval_msd)) %>%
  filter(focal_n > 2)

metaData <-
  do.call(rbind,
          list(metaDat_values, metaDat_iqr2msd, metaDat_msd))

#---
#-- Create a 1:1 plot of results with and without correction
#---

#' We want to include this figure before dropping data that can't be used in the
#' meta-analysis (i.e. point estimates or metrics with no variance)

#-- Individual points (i.e. with no +- SD)
indiv_points <- biasCorrEval_values %>%
  filter(is.na(focal_sd))
ns <-
  c(unique(biasCorrEval_msd$study_id),
    unique(indiv_points$study_id))
indiv_points2 <-
  biasCorrEval[biasCorrEval$study_id %in% setdiff(unique(biasCorrEval$study_id), ns), ]

#-- Check number of studies
length(unique(biasCorrEval$study_id))

#--- Plot all raw metric comparisons
all_data_comp_p <- biasCorrEval_msd %>%
  mutate(
    focal_1sd_u = focal_m + focal_sd,
    focal_1sd_l = focal_m - focal_sd,
    naive_1sd_u = naïve_m + naïve_sd,
    naive_1sd_l = naïve_m - naïve_sd,
    focal_1sd_u = ifelse(focal_1sd_u > 1, 1, focal_1sd_u),
    focal_1sd_l = ifelse(focal_1sd_l < 0, 0, focal_1sd_l),
    naive_1sd_u = ifelse(naive_1sd_u > 1, 1, naive_1sd_u),
    naive_1sd_l = ifelse(naive_1sd_l < 0, 0, naive_1sd_l)
  ) %>%
  mutate(
    test_dataset =
      recode(test_dataset,
             Independent = "Test data = Independent",
             Internal = "Test data = Internal")
  ) %>%
  ggplot(
    aes(
      x = naïve_m,
      y = focal_m,
      colour = bias_correction_method,
      shape = performance_metric
    ),
    alpha = 0.5
  ) +
  geom_errorbarh(aes(xmin = naive_1sd_l,
                     xmax = naive_1sd_u),
                 alpha = 0.2) +
  geom_errorbar(aes(ymin = focal_1sd_l,
                    ymax = focal_1sd_u),
                alpha = 0.2) +
  geom_point() +
  geom_point(
    aes(
      x = naïve_m,
      y = focal_m,
      colour = bias_correction_method,
      shape = performance_metric
    ),
    alpha = 0.5,
    data = indiv_points
  ) +
  geom_point(
    aes(
      x = naïve_m,
      y = focal_m,
      colour = bias_correction_method,
      shape = performance_metric
    ),
    alpha = 0.5,
    data = indiv_points2
  ) +
  geom_abline() +
  facet_wrap(~ test_dataset, nrow = 2) +
  theme_bw(base_size = 14) %+replace%
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "right") +
  scale_colour_colorblind(name = "Bias correction method") +
  scale_shape(name = "Test statistic") +
  coord_equal() +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab(expression(
    atop(
      "Test statistic (mean \u00B11sd)",
      "No spatial sampling bias correction"
    )
  )) +
  ylab(expression(
    atop(
      "Test statistic (mean \u00B11sd)",
      "Spatial sampling bias correction"
    )
  ))

#--- Save for supplementary material
ggsave(
  filename = "Outputs/Figures/Figure_S1.png",
  plot = all_data_comp_p,
  width = 8,
  height = 10
)

#---
#-- Calculate Hedge's g (standardised effect size)
#---

#' Hedge's g is a small sample bias corrected standardised effect size

#--- Quick sense check
biasCorrEval_msd %>%
  #filter(performance_metric == "AUC") %>%
  rowwise() %>%
  mutate(sc = ifelse(focal_m >= naïve_m, 1, 0)) %>%
  ungroup() %>%
  dplyr::select(sc) %>%
  table()
# Looks like more cases of bias correction making things worse.

#--- Calculate effect sizes: Cohen's d and Hedge's g
biasCorrEval_d <-
  biasCorrEval_msd  %>%
  filter(!is.na(naïve_m)) %>%
  mutate(m_diff = focal_m - naïve_m) %>%
  rowwise() %>%
  mutate(
    S_within = sqrt((((naïve_n - 1) * naïve_sd ^ 2
    ) +
      ((focal_n - 1) * focal_sd ^ 2
      )) / ((naïve_n - 1) + (focal_n - 1))),
    J = 1 - (3 / (4 * (naïve_n + focal_n - 2) - 1)),
    d = (m_diff / S_within),
    d = ifelse(is.na(d), m_diff, d),
    g = (m_diff / S_within) * J,
    g = ifelse(is.na(g), m_diff, g),
    Vd = (naïve_n + focal_n) / ((naïve_n * focal_n) + (d ^ 2 / (
      2 * (naïve_n + focal_n)
    ))),
    SEd = sqrt(Vd),
    invSEd = 1 / SEd,
    Vg = Vd * J ^ 2,
    SEg = sqrt(Vg),
    invSEg = 1 / SEg
  )

#---
#-- Estimated effect size (mean +/- 95CI) of mean difference
#---

#--- Subset to single metric per study
#' Because there are sometimes multiple metrics per study we just use one metric
#' in the meta-analysis. useInMa == 1 selects the primary metric per study,
#' which is AUC (i.e. the most common).
ssb_cor_d <- filter(biasCorrEval_d, useInMa == 1)

#--- Quick sense check
ssb_cor_d %>%
  rowwise() %>%
  mutate(sc = sign(g)) %>%
  ungroup() %>%
  group_by(test_dataset, bias_correction_method) %>%
  dplyr::select(sc) %>%
  table()

#--- Independent test data
#' Calculate for each study the grand mean and sd of Hedge's g to show study
#' specific mean and variation in effect sizes.
fx_d <- ssb_cor_d %>%
  rowwise %>%
  mutate(pN = focal_n + naïve_n) %>%
  group_by(test_dataset, bias_correction_method, study_id) %>%
  summarise(gm = weighted.mean(g, pN, na.rm = TRUE),
            gsd = sqrt(
              weighted.mean(Vg + g ^ 2, pN, na.rm = TRUE) -
                weighted.mean(g, pN, na.rm = TRUE) ^ 2
            )) %>%
  group_by(test_dataset) %>%
  mutate(
    id = seq_along(gm),
    # Create plotting order
    wd = ifelse(test_dataset == "Test data = Internal", 0.5, 0.2),
    # for CI width
    bias_correction_method =
      fct_recode(
        bias_correction_method,
        "Spp. target grp" = "Species_grp_records",
        "Covariate (cond. pred.)" =  "Covariate_condPred",
        "Feature (buffer or dist. decay)" = "Feature_bufferOrDD",
        "Occurrence (buffer or dist. decay)" = "Occurrence_bufferOrDD",
        "Occurrence filter" = "Occurrence_filter",
        "Effort (buffer or dist. decay)" = "Effort_bufferOrDD"
      )
  )

#---
#-- Plot the results for the independent evaluations
#---

#-- Filter independent data and create id label
fx_dat_ind <-
  fx_d %>% filter(test_dataset == "Test data = Independent") %>%
  mutate(id = rev(id))

#-- Produce text labels for points that have extreme range
ind_lab <- fx_dat_ind %>%
  filter(abs(gm) > 5) %>%
  mutate(
    l95 = gm - 1.96 * gsd,
    u95 = gm + 1.96 * gsd,
    gm = round(gm, 2),
    l95 =  round(l95, 2),
    u95 =  round(u95, 2),
    lab_i = paste0(gm, " (", l95, ", ", u95, ")"),
    y_pos = 14 * sign(gm),
    hj = ifelse(sign(gm) > 0, 1, 0)
  )

#---
#-- Plot the results for the independent evaluations
#---

#-- Filter independent data and create id label
fx_dat_ind <- fx_d %>% filter(test_dataset == "Test data = Independent") %>%
  mutate(id = rev(id)) 

#-- Produce text labels for points that have extreme range
ind_lab <- fx_dat_ind %>%
  filter(abs(gm) > 5) %>%
  mutate(
    l95 = gm - 1.96 * gsd,
    u95 = gm + 1.96 * gsd,
    gm = round(gm, 2),
    l95 =  round(l95, 2),
    u95 =  round(u95, 2),
    lab_i = paste0(gm, " (", l95, ", ", u95, ")"),
    y_pos = 14 * sign(gm),
    hj = ifelse(sign(gm) > 0, 1, 0)
  )

#-- Remove extreme points (which will be given text labels)
fx_dat_ind <- mutate(fx_dat_ind, gm = ifelse(abs(gm) > 5, NA, gm))

#-- Plot, adding in text for extreme points
fx_plot_ind <-  ggplot(fx_dat_ind) +
  geom_point(aes(x = id,
                 y = gm,
                 colour = bias_correction_method),) +
  geom_errorbar(aes(
    x = id,
    ymin = gm - 1.96 * gsd,
    ymax = gm + 1.96 * gsd,
    colour = bias_correction_method
  ),
  width = 0,
  size = 0.5) +
  # Add in text labels
  geom_text(
    aes(x = id,
        y = y_pos,
        label = lab_i,
        hjust = hj),
    data = ind_lab
  ) +
  ylim(-14, 14) + # set to be the same across independent and interval panels
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw(base_size = 14) %+replace%
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()#,legend.position = "none"
  ) +
  coord_flip() +
  scale_color_tableau(name = "SSB correction\nmethod", drop=FALSE) +
  scale_x_continuous(breaks = fx_dat_ind$id, 
                     labels = fx_dat_ind$bias_correction_method)

#---
#-- Plot the results for the internal evaluations
#---

#-- Filter internal data and create id label
fx_dat_int <- fx_d %>% filter(test_dataset == "Test data = Internal") %>%
  mutate(id = rev(id)) 

#-- Produce text labels for points that have extreme range
int_lab <- fx_dat_int %>%
  filter(abs(gm) > 6) %>%
  mutate(
    l95 = gm - 1.96 * gsd,
    u95 = gm + 1.96 * gsd,
    gm = round(gm, 2),
    l95 =  round(l95, 2),
    u95 =  round(u95, 2),
    lab_i = paste0(gm, " (", l95, ", ", u95, ")"),
    y_pos = 14 * sign(gm),
    hj = ifelse(sign(gm) > 0, 1, 0)
  )

#-- Remove extreme points (which will be given text labels)
fx_dat_int <- mutate(fx_dat_int, gm = ifelse(abs(gm) > 5, NA, gm))

#-- Plot, adding in text for extreme points
fx_plot_int <-
  ggplot(fx_dat_int) +
  geom_point(aes(x = id,
                 y = gm,
                 colour = bias_correction_method),) +
  geom_errorbar(aes(
    x = id,
    ymin = gm - 1.96 * gsd,
    ymax = gm + 1.96 * gsd,
    colour = bias_correction_method,
  ),
  width = 0,
  size = 0.5) +
  # Add in text labels
  geom_text(
    aes(x = id,
        y = y_pos,
        label = lab_i,
        hjust = hj),
    data = int_lab
  ) +
  ylim(-14, 14) +  # set to be the same across independent and interval panels
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw(base_size = 14) %+replace%
  theme(
    axis.title.y = element_blank()#,legend.position = "none"
  ) +
  coord_flip() +
  scale_color_tableau(name = "SSB correction\nmethod", drop=FALSE)  +
  scale_x_continuous(breaks = fx_dat_int$id, 
                     labels = fx_dat_int$bias_correction_method) +
  ylab(expression("Hedge's g \u00b195CI"))

#--- Create and save Figure 2 panel plot
fig_2 <- fx_plot_ind / fx_plot_int +
  plot_layout(guides = "collect", heights = c(1, 2.5)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
  theme(legend.position = "none", legend.title = element_blank())
#--- Save for supplementary material
ggsave(
  filename = "Outputs/Figures/Figure_2.png",
  plot = fig_2,
  width = 8,
  height = 10
)
                
#---
#-- Test mean difference between performance measured using dif. test dataset
#---

#--- Model testing effect for test dataset on Hedge's g
#' Weights are the inverse SE of Hedge's g and study_id is a random intercept
m1 <-
  brm(
    bf(g | weights(invSEg) ~ test_dataset + (1 | study_id)),
    data = ssb_cor_d,
    family = student(),
    chains = 4,
    iter = 60000,
    thin = 10,
    cores = 4,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 14)
  )
pp_check(m1)
pp_check(m1, type = "ecdf_overlay")
summary(m1)
fixef(m1)
#save(m1, file = "Outputs/brms_g_testdat.Rdata", compress = "xz")

#-- Hypothesis test of internal < 0
#load("Outputs/brms_g_testdat.Rdata")
hypothesis(m1, "test_datasetTestdataEQInternal < 0")

#-- Treatment effects
cond_eff_int <- conditional_effects(m1)[[1]]

#---
#-- Test mean difference between performance measured using dif. SDM methods
#---

#--- Prepare data (add in metadata)
ssb_d_met <- left_join(ssb_cor_d, metaData)

#--- Separate into model classes
ssb_d_met$modClass <- NA

#-- Maxent
mx <- c("Maxent", "MaxEnt", "Maxlike")
ssb_d_met$modClass[ssb_d_met$sdm_method %in% mx] <-
  "Maxent"

#-- GLM/GAM
glm_gam <-
  c(
    "GLM",
    "GLMM",
    "Lasso GLM",
    "ENFA",
    "GAM",
    "GAMM",
    "FDA",
    "Hierarchical model",
    "Occupancy-detection",
    "Multi-species hierarchical model"
  )
ssb_d_met$modClass[ssb_d_met$sdm_method %in% glm_gam] <-
  "GLM/GAM"

#-- Machine Learn./Classif.
ml <-
  c(
    "Elastic net",
    "ANN",
    "BRT",
    "CART",
    "CTA",
    "CTREE",
    "GARP",
    "GBM",
    "MARS",
    "KNN",
    "NNET",
    "SVM",
    "SRE",
    "Naïve bayes",
    "RF"
  )
ssb_d_met$modClass[ssb_d_met$sdm_method %in% ml] <-
  "Machine Learn./Classif."

#-- PPM
ppm <- c(
  "PPM",
  "DWPR",
  "Bayesian Image restoration",
  "GP",
  "CAR",
  "Spatial log-Gaussian cox process"
)
ssb_d_met$modClass[ssb_d_met$sdm_method %in% ppm] <-
  "Point Process Models"

#-- Multimodel
mm <- "ANN|BRT|FDA|GAM|GLM|MARS|Maxent|RF"
ssb_d_met$modClass[ssb_d_met$sdm_method == mm] <-
  "Ensemble"

#' Test the effect of SDM method on Hedge's g
m2 <-
  brm(
    bf(g | weights(invSEg) ~  modClass + (1 | study_id)),
    data = ssb_d_met,
    family = student(),
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.95,
                   max_treedepth = 14),
    seed = 123
  )
pp_check(m2)
pp_check(m2, type = "ecdf_overlay")
summary(m2)
fixef(m2)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------- END ------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#