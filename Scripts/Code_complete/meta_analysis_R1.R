#------------------------------------------------------------------------------#
#
#     Spatial sampling bias systematic review and meta analysis
#
#------------------------------------------------------------------------------#

# --- Load packages
library(tidyverse)
library(readr)
library(readxl)
library(openxlsx)
library(data.table)
library(patchwork)

#------------------------------------------------------------------------------#
# Load literature extraction data
#------------------------------------------------------------------------------#

#' The systematic review was conducted using Web of Science and Scopus databases
#' on the 13/02/2023 using the following search string:
#'
#' ALL=(("species distribution*" OR SDM OR "environmental niche" OR ENM OR 
#' "resource selection" OR "habitat selection" OR suitability OR occurrence) 
#' AND ("presence-only" OR “presence data” OR "presence-background" OR 
#' “pseudo absence” OR opportunistic OR “citizen science” OR preferential 
#' OR maxent OR biomod)) 

#-------------------------------------------------------------------------------
# Summary statistics
#-------------------------------------------------------------------------------

# --- Load data for full systematic review search results and filtering
sr_search_dat <-
  read_excel("Data/sysRev_databases_R1/dataset_bias_correction_systematic_review_R1_4_8_23.xlsx",
             sheet = "search_results_filter")
sr_search_dat$title <- tolower(sr_search_dat$title)
sr_search_dat <- sr_search_dat %>% 
  dplyr::select(study_id, select_title, select_abstract, select_content, comparision, comments)

#---
#-- Basic summary statistics for systematic review filtering process
#---

#--- Number of studies in initial search (non duplicated)
nrow(sr_search_dat) # 8570

#--- Number of studies filtered by title
sr_search_dat %>%
  filter(select_title == "Y") %>%
  nrow() # 5555

#--- Number of studies filtered by abstract
sr_search_dat %>%
  filter(select_abstract == "Y") %>%
  nrow() # 4964

#--- Number of studies filtered by content
sr_search_dat %>%
  filter(select_content == "Y") %>%
  nrow() # 935

#--- Number of studies filtered by content
sr_search_dat %>%
  filter(comparision == "Y") %>%
  nrow() # 73


#-------------------------------------------------------------------------------
# Meta-analysis: the effect of SSB correction on SDM performance
#-------------------------------------------------------------------------------

# -- Load new search
sr_dat <-
  read_excel("Data/sysRev_databases_R1/dataset_bias_correction_systematic_review_R1_4_8_23.xlsx",
             sheet = "data_extraction")
sr_dat$bibtexkey <- gsub(" ", "", sr_dat$bibtexkey )

unique(sr_dat$bias_correction_method)

#-#-# Covert to factors and level
sr_dat <- sr_dat %>%
  mutate(
    bias_correction_method = factor(
      bias_correction_method,
      levels =
        c(
          "Species_grp_records",
          "Feature_bufferOrDD",   
          "Covariate",
          "Covariate_condPred",
          "Occurrence_filter_spa + Occurrence_bufferOrDD",
          "Occurrence_filter_spa + Species_grp_records",
          "Occurrence_bufferOrDD",
          "Occurrence_filter_spa",    
          "Unclear",
          "Occurrence_filter_spa + Covariate",
          "Effort_bufferOrDD",
          "Occurrence_filter_env",
          "Effort_bufferOrDD + Feature_bufferOrDD",
          "Occurrence_filter_spa + Feature_bufferOrDD + Species_grp_records",
          "Occurrence_bufferOrDD + Covariate", 
          "Occurrence_filter_spa + Occurrence_bufferOrDD + Covariate",
          "Feature_bufferOrDD + Occurrence_filter",
          "ZI"
        )
    ),
    bias_correction_method = fct_recode(
      bias_correction_method,
      "Multiple" = "Occurrence_filter_spa + Occurrence_bufferOrDD",
      "Multiple" = "Occurrence_filter_spa + Species_grp_records",
      "Multiple" = "Occurrence_filter_spa + Covariate",
      "Multiple" = "Effort_bufferOrDD + Feature_bufferOrDD",
      "Multiple" = "Occurrence_filter_spa + Feature_bufferOrDD + Species_grp_records",
      "Multiple" = "Occurrence_filter_spa + Species_grp_records",
      "Multiple" = "Occurrence_bufferOrDD + Covariate",
      "Multiple" = "Occurrence_filter_spa + Occurrence_bufferOrDD + Covariate",
      "Multiple" = "Feature_bufferOrDD + Occurrence_filter"
    ),
    bias_correction_method =
      fct_recode(
        bias_correction_method,
        "OccFilter(Geo)" = "Occurrence_filter_spa",
        "OccFilter(Env)" = "Occurrence_filter_env",
        "AdjBkGrd(SppTgGrp)" = "Species_grp_records",
        "AdjBkGrd(Occurrence)" = "Occurrence_bufferOrDD",
        "AdjBkGrd(Feature)" = "Feature_bufferOrDD",
        "AdjBkGrd(Effort)" = "Effort_bufferOrDD",
        "Covariate(CondPred)" =  "Covariate_condPred",
        "Zero-inflated" = "ZI"
      )
  )

#---
#-- Prepare model data
#---

#-- Subset data to just those with comparisons between uncorrected and corrected
biasCorrEval <- sr_dat %>%
  filter(comparison_to_naive == "Yes",
         useInMa == 1) %>%
  filter(test_dataset != "Simulated") %>%
  mutate(
    corrected_m = as.numeric(corrected_m),
    corrected_sd = as.numeric(corrected_sd),
    corrected_n = as.numeric(corrected_n),
    corrected_lci = as.numeric(corrected_lci),
    corrected_uci = as.numeric(corrected_uci),
    uncorrected_m = as.numeric(uncorrected_m),
    uncorrected_sd = as.numeric(uncorrected_sd),
    uncorrected_n = as.numeric(uncorrected_n),
    uncorrected_lci = as.numeric(uncorrected_lci),
    uncorrected_uci = as.numeric(uncorrected_uci)
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
uncorrected_m <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(uncorrected_m = mean(uncorrected_m, na.rm = TRUE))
uncorrected_sd <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(uncorrected_sd = sd(uncorrected_m, na.rm = TRUE))
uncorrected_n <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(uncorrected_n = length(uncorrected_m))
corrected_m <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(corrected_m = mean(corrected_m, na.rm = TRUE))
corrected_sd <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(corrected_sd = sd(corrected_m, na.rm = TRUE))
corrected_n <- biasCorrEval_values %>%
  group_by(across(all_of(grpByCols))) %>%
  summarise(corrected_n = length(corrected_m))

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
  list(corrected_m, corrected_sd, corrected_n,
       uncorrected_m, uncorrected_sd, uncorrected_n)
)

#--- For results given as Quantile + Q1 and Q3 (i.e. from boxplots)
biasCorrEval_iqr2msd <- biasCorrEval %>%
  filter(quant_metric == "median_iqr") %>%
  rename(
    uncorrected_q1 = "uncorrected_lci",
    uncorrected_q3 = "uncorrected_uci",
    corrected_q1 = "corrected_lci",
    corrected_q3 = "corrected_uci"
  ) %>%
  mutate(
    uncorrected_m =  (uncorrected_q1 + uncorrected_m + uncorrected_q3) / 3,
    uncorrected_sd = (uncorrected_q3 - uncorrected_q1) / (2 * qnorm((0.75 * uncorrected_n - 0.125) / (uncorrected_n + 0.25), 0, 1
    )),
    corrected_m =  (corrected_q1 + corrected_m + corrected_q3) / 3,
    corrected_sd = (corrected_q3 - corrected_q1) / (2 * qnorm((0.75 * corrected_n - 0.125) / (corrected_n + 0.25), 0, 1
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
  filter(corrected_n > 2)

metaData <-
  do.call(rbind,
          list(metaDat_values, metaDat_iqr2msd, metaDat_msd))

summary(unique(metaData$max_number_rec))

#---
#-- Create a 1:1 plot of results with and without correction
#---

#' We want to include this figure before dropping data that can't be used in the
#' meta-analysis (i.e. point estimates or metrics with no variance)

#-- Individual points (i.e. with no +- SD)
indiv_points <- biasCorrEval_values %>%
  filter(is.na(corrected_sd))
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
    corrected_1sd_u = corrected_m + corrected_sd,
    corrected_1sd_l = corrected_m - corrected_sd,
    uncorrected_1sd_u = uncorrected_m + uncorrected_sd,
    uncorrected_1sd_l = uncorrected_m - uncorrected_sd,
    corrected_1sd_u = ifelse(corrected_1sd_u > 1, 1, corrected_1sd_u),
    corrected_1sd_l = ifelse(corrected_1sd_l < 0, 0, corrected_1sd_l),
    uncorrected_1sd_u = ifelse(uncorrected_1sd_u > 1, 1, uncorrected_1sd_u),
    uncorrected_1sd_l = ifelse(uncorrected_1sd_l < 0, 0, uncorrected_1sd_l)
  ) %>%
  mutate(
    test_dataset =
      recode(test_dataset,
             Independent = "Test data = Independent",
             Internal = "Test data = Internal")
  ) %>%
  ggplot(
    aes(
      x = uncorrected_m,
      y = corrected_m,
      colour = bias_correction_method,
      shape = performance_metric
    ),
    alpha = 0.5
  ) +
  geom_errorbarh(aes(xmin = uncorrected_1sd_l,
                     xmax = uncorrected_1sd_u),
                 alpha = 0.2) +
  geom_errorbar(aes(ymin = corrected_1sd_l,
                    ymax = corrected_1sd_u),
                alpha = 0.2) +
  geom_point() +
  geom_point(
    aes(
      x = uncorrected_m,
      y = corrected_m,
      colour = bias_correction_method,
      shape = performance_metric
    ),
    alpha = 0.5,
    data = indiv_points
  ) +
  geom_point(
    aes(
      x = uncorrected_m,
      y = corrected_m,
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
  scale_colour_brewer(palette = "Set3", name = "Bias correction method") +
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
  filename = "Outputs/Figures_R2/Figure_S1.png",
  plot = all_data_comp_p,
  width = 10,
  height = 11
)

#---
#-- Calculate Hedge's g (standardised effect size)
#---

#' Hedge's g is a small sample bias corrected standardised effect size

#--- Quick sense check
biasCorrEval_msd %>%
  #filter(performance_metric == "AUC") %>%
  rowwise() %>%
  mutate(sc = ifelse(corrected_m >= uncorrected_m, 1, 0)) %>%
  ungroup() %>%
  dplyr::select(sc) %>%
  table()
# Looks like more cases of bias correction making things worse.

#--- Calculate effect sizes: Cohen's d and Hedge's g
biasCorrEval_d <-
  biasCorrEval_msd  %>%
  filter(!is.na(uncorrected_m)) %>%
  mutate(m_diff = corrected_m - uncorrected_m) %>%
  rowwise() %>%
  mutate(
    S_within = sqrt((((uncorrected_n - 1) * uncorrected_sd ^ 2
    ) +
      ((corrected_n - 1) * corrected_sd ^ 2
      )) / ((uncorrected_n - 1) + (corrected_n - 1))),
    J = 1 - (3 / (4 * (uncorrected_n + corrected_n - 2) - 1)),
    d = (m_diff / S_within),
    d = ifelse(is.na(d), m_diff, d),
    g = (m_diff / S_within) * J,
    g = ifelse(is.na(g), m_diff, g),
    Vd = (uncorrected_n + corrected_n) / ((uncorrected_n * corrected_n) + (d ^ 2 / (
      2 * (uncorrected_n + corrected_n)
    ))),
    SEd = sqrt(Vd),
    invSEd = 1 / SEd,
    Vg = Vd * J ^ 2,
    SEg = sqrt(Vg),
    invSEg = 1 / SEg
  )

# --- AUC mean difference - independent
ind_auc <- biasCorrEval_d %>% 
  filter(test_dataset == "Test data = Independent",
         performance_metric == "AUC")
summary(ind_auc$m_diff)
ind_cor <- biasCorrEval_d %>% 
  filter(test_dataset == "Test data = Independent",
         performance_metric == "COR")
summary(ind_cor$m_diff)

# --- AUC mean difference - internal
int_auc <- biasCorrEval_d %>% 
  filter(test_dataset == "Test data = Internal",
         performance_metric == "AUC")
summary(int_auc$m_diff)

ind_cor <- biasCorrEval_d %>% 
  filter(test_dataset == "Test data = Internal",
         performance_metric == "COR")
summary(ind_cor$m_diff)

#---
#-- Estimated effect size (mean +/- 95CI) of mean difference
#---

#--- Subset to single metric per study
#' Because there are sometimes multiple metrics per study we just use one metric
#' in the meta-analysis. useInMa == 1 selects the primary metric per study,
#' which is AUC (i.e. the most common).
ssb_cor_d <- filter(biasCorrEval_d, useInMa == 1)

# -- 
s_id <- unique(ssb_cor_d$study_id)
study_index <- sr_dat %>% 
  filter(study_id %in% {{s_id}}) %>% 
  dplyr::select(study_id, bibtexkey, doi) %>% 
  dplyr::distinct()
write.table(study_index, "Outputs/study_index_figure_S2.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

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
  mutate(pN = corrected_n + uncorrected_n) %>%
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
   bias_correction_method =
    fct_relevel(
       bias_correction_method,
       "OccFilter(Geo)",
       "OccFilter(Env)" ,
       "AdjBkGrd(SppTargetGrp)",
       "AdjBkGrd(occurrence_weight_or_buffer)",
       "AdjBkGrd(feature_weight_or_buffer)",
       "AdjBkGrd(effort_weight_or_buffer)",
       "Covariate(CondPred)",
       "Covariate",
       "Multiple",
       "Zero-inflated",
       "Unclear"
     )
 )


#---
#-- Plot the results for the independent evaluations
#---

#-- Filter independent data and create id label
fx_dat_ind <-
  fx_d %>% filter(test_dataset == "Test data = Independent") %>%
  mutate(id = rev(id))

# Number of independent
length(unique(fx_dat_ind$study_id[fx_dat_ind$test_dataset == "Test data = Independent"]))

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
  filter(abs(gm) > 7) %>%
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
fx_dat_ind <- mutate(fx_dat_ind, gm = ifelse(abs(gm) > 7, NA, gm))

#-- Plot, adding in text for extreme points
fx_plot_ind <-  ggplot(fx_dat_ind) +
  geom_errorbar(aes(
    x = id,
    ymin = gm - 1.96 * gsd,
    ymax = gm + 1.96 * gsd,
    colour = bias_correction_method
  ),
  width = 0,
  linewidth = 0.5) +
  geom_point(aes(x = id,
                 y = gm),
             colour = "black") +
  # Add in text labels
  geom_text(
    aes(x = id,
        y = y_pos,
        label = lab_i,
        hjust = hj),
    data = ind_lab
  ) +
  #For supp mat add in study id labels
  geom_text(
    aes(x = id,
        y = -16,
        label = study_id,
        hjust = 0.5,
        vjust = 0),
    size = 2
  ) +
  ylim(-16, 16) + # set to be the same across independent and interval panels
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw(base_size = 14) %+replace%
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()#,legend.position = "none"
  ) +
  coord_flip() +
  scale_color_brewer(palette ="Paired" ,name = "SSB correction\nmethod", drop=FALSE) +
  scale_x_continuous(breaks = fx_dat_ind$id, 
                     labels = fx_dat_ind$bias_correction_method)

#---
#-- Plot the results for the internal evaluations
#---

#-- Filter internal data and create id label
fx_dat_int <- fx_d %>% filter(test_dataset == "Test data = Internal") %>%
  mutate(id = rev(id))  

# Number of internal
length(unique(fx_dat_int$study_id[fx_dat_int$test_dataset == "Test data = Internal"]))


#-- Produce text labels for points that have extreme range
int_lab <- fx_dat_int %>%
  filter(abs(gm) > 7) %>%
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
fx_dat_int <- mutate(fx_dat_int, gm = ifelse(abs(gm) > 7, NA, gm))

#-- Plot, adding in text for extreme points
fx_plot_int <-
  ggplot(fx_dat_int) +
  geom_errorbar(aes(
    x = id,
    ymin = gm - 1.96 * gsd,
    ymax = gm + 1.96 * gsd,
    colour = bias_correction_method,
  ),
  width = 0,
  size = 0.5) +
  geom_point(aes(x = id,
                 y = gm),
             colour = "black") +
  # Add in text labels
  geom_text(
    aes(x = id,
        y = y_pos,
        label = lab_i,
        hjust = hj),
    data = int_lab
  ) + #For supp mat add in study id labels
  geom_text(
    aes(x = id,
        y = -16,
        label = study_id,
        hjust = 0.5,
        vjust = 0),
    size = 2
  ) +
  ylim(-16, 16) +  # set to be the same across independent and interval panels
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw(base_size = 14) %+replace%
  theme(
    axis.title.y = element_blank()#,legend.position = "none"
  ) +
  coord_flip() +
  scale_color_brewer(palette ="Paired" ,name = "SSB correction\nmethod", drop=FALSE) +
  scale_x_continuous(breaks = fx_dat_int$id, 
                     labels = fx_dat_int$bias_correction_method) +
  ylab(expression("Grand Mean Hedges' g \u00b195CI"))

#--- Create and save Figure 2 panel plot
fig_1 <- fx_plot_ind / fx_plot_int +
  plot_layout(guides = "collect", heights = c(1, 4)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
  theme(legend.position = "none", legend.title = element_blank())

#--- Save for supplementary material
ggsave(
  filename = "Outputs/Figures_R2/Figure_S2_metaAnalysis.png",
  plot = fig_1,
  width = 8,
  height = 12
)
   
#---
#-- Test mean overall between hedge's g
#---

#--- Independent - Model testing effect for test dataset on Hedge's g
#' Weights are the inverse SE of Hedge's g and study_id is a random intercept
ind_dat <- ssb_cor_d %>% filter(test_dataset == "Test data = Independent")

m1_ind <-
  brm(
    bf(g | weights(invSEg) ~ 1 + (1 | study_id)),
    data = ind_dat,
    family = student(),
    chains = 4,
    iter = 12000,
    thin = 10,
    cores = 4,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 14)
  )
# pp_check(m1)
# pp_check(m1, type = "ecdf_overlay")
# summary(m1_ind)
# fixef(m1_ind)
hypothesis(m1_ind, "Intercept > 0")
save(m1_ind, file = "Outputs/brms_g_testdat_Independent.Rdata", compress = "xz")


#--- Internal - Model testing effect for test dataset on Hedge's g
#' Weights are the inverse SE of Hedge's g and study_id is a random intercept
int_dat <- ssb_cor_d %>% filter(test_dataset == "Test data = Internal")
m1_int <-
  brm(
    bf(g | weights(invSEg) ~ 1 + (1 | study_id)),
    data = int_dat,
    family = student(),
    chains = 4,
    iter = 12000,
    thin = 10,
    cores = 4,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 14)
  )
# pp_check(m1)
# pp_check(m1, type = "ecdf_overlay")
# summary(m1)
# fixef(m1_int)
hypothesis(m1_int, "Intercept > 0")
save(m1_int, file = "Outputs/brms_g_testdat_Internal.Rdata", compress = "xz")
             
#---
#-- Test mean difference between performance measured using dif. test dataset
#---

#--- Model testing effect for test dataset on Hedge's g
#' Weights are the inverse SE of Hedge's g and study_id is a random intercept
m1_test <-
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
# pp_check(m1)
# pp_check(m1, type = "ecdf_overlay")
# summary(m1)
# fixef(m1)
save(m1_test, file = "Outputs/brms_g_testdat.Rdata", compress = "xz")

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
    "uncorrected bayes",
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
    iter = 12000,
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