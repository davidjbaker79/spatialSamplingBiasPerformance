#------------------------------------------------------------------------------#
#
#   Simulation of the spatial sampling bias analysis
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# (1) Simulate 100 environments (we do this first as this is slow(ish))
#------------------------------------------------------------------------------#

#--- Functions
source("Scripts/generate_environment_R1.R")
source("Scripts/create_enviro_features_R1.R")

#--- Run N = 100
#' Environment is created and then the variables altered to create more 
#' realistic variables
#' V1 is given a stronger gradient
#' V2 is given an upper 'elevation' cut off
#' V3 is given a 'latitude' gradient
#' 'Roads' are created through the landscape
for(i in 1:100) {
  
  # Create environmental data
  enviro <- generate.enviro(n = 50, phi = c(0.1, 0.75, 0.5))
  
  # Create a stronger "elevation" gradient
  enviro$V1 <- round(1 / (1 + exp((enviro$V1 - 0.5) / -0.1)), 3)
  #ggplot(enviro) + geom_raster(aes(x = X, y = Y, fill = V1))
  
  # Create "habitat" layer with upper elevation cut off
  enviro$V2 <-
    round(enviro$V2 * (1 / (1 + exp((enviro$V1 - 0.5) / 0.25
    ))), 3)
  #ggplot(enviro) + geom_raster(aes(x=X,y=Y, fill=V2))
  
  # Create gradient (e.g. latitude temperature gradient)
  enviro$V3 <- enviro$V3 * scales::rescale(enviro$X * enviro$Y, c(0.1,0.9))
  #ggplot(enviro) + geom_raster(aes(x=X,y=Y, fill=V3))
  
  # Create linear feature through "elevation"
  l.feat <- create.linear.features(enviro[c("X", "Y", "V1")])
  
  # Create area features (e.g. PAs) in reference to "forest"
  a.feat <- create.area.features(enviro[c("X", "Y", "V2")])
  
  # Combine
  enviro <-
    Reduce(function(...)
      left_join(...), list(enviro, l.feat, a.feat))
  
  save(
    enviro,
    file = paste0("Outputs/simulated_environments/env_", i, "_R1.RData"),
    compress = "xz"
  )
  
}

#enviro <- get(load(paste0("Outputs/simulated_environments/env_", 1, "_R1.RData")))



# --- Measure SAC
i = 4
for (i  in 1:100) {
  
  # Load environmental data
  env_i <- 
    get(load(paste0("Outputs/simulated_environments/env_", i, "_R1.RData")))
  env_i <- rast(env_i[,c(2:3,5)], type = "xyz")
  crs(env_i) <- "EPSG:4326"
  cv_sac <- cv_spatial_autocor(
    r = env_i,
    num_sample = 5000L,
    deg_to_metre = 111325,
    plot = TRUE,
    progress = TRUE
  )
  
  sb1 <- cv_spatial(
    x = pa_dat,
    column = "det",
    size = cv_sac$range,
    hexagon = FALSE,
    k = 10,
    selection = "systematic",
    offset = 1,
    iteration = 50
  )
  
}

#------------------------------------------------------------------------------#
# (2) Simulate 1000 virtual species 
#------------------------------------------------------------------------------#

# --- Simulate 100 species (different niche breadths) per 100 environments
for(i in 1:100) { 
  # -- Load environmental data
  enviro_i <-
    get(load(
      paste0("Outputs/simulated_environments/env_", i, "_R1.RData")
    ))
  for(j in 1:100) { print(i * j)
    # -- Simulate species with different niche breadths
    vs_i <- species.occ(enviro_i)
    save(
      vs_i,
      file = paste0("Outputs/virtual_species/vspp_", j, "_env_", i, ".RData"),
      compress = "xz"
    )
  }
}

x_i <- get(load("Outputs/virtual_species/vspp_2_env_1.RData"))[[1]]
x_i <- left_join(enviro_i[,c("id", "X", "Y")], x_i)
# ggplot() +
#   geom_raster(aes(x = X, y = Y, fill = poc), data = x_i) +
#   geom_point(aes(x = X, y = Y), size = 0.1,colour = "red", data= x_i[x_i$pa == 1, ]) +
#   scale_fill_viridis_c()

# --- Subsample to get species with different types of niche breadths

# - Get prevelance for each virtual species
f <- list.files("Outputs/virtual_species/", full.names = TRUE)
spp_prev <- lapply(f, function(i){ print(i)
  f_prev <- get(load(i))[[4]]
  f_str <- strsplit(i, split = "_")[[1]]
  env <- sub(".RData","",f_str[length(f_str)])
  spp <- f_str[3]
  data.frame(env = env, spp = spp, prev = f_prev)
})
spp_prev_df <- do.call(rbind, spp_prev)

# -- Calculate sample weightings based on number of species with similar prev
spp_prev_df <- spp_prev_df %>% 
  mutate(prev_cat = plyr::round_any(prev, .05))
spp_prev_df <- spp_prev_df[spp_prev_df$prev_cat <= 0.5,]

# -- Random sample 50 virtual species from each prevalence category
spp_prev_select <- spp_prev_df %>% 
  group_by(prev_cat) %>% 
  sample_n(50)
save(spp_prev_select, file = "Outputs/species_prev_selection.RData")

#------------------------------------------------------------------------------#
#  (2) Second, run simulation models of SSB correction
#------------------------------------------------------------------------------#

#' Generate bias to linear and site features.
#' Test correction based on:
#' 1) degree to which bias correction approaches the true bias 
#' 2) how correlated the species niches is to the true bias
#' 3) sample size of presences
#' 4) how strong the sampling bias is compare to the correction


#--- Source functions
source("Scripts/Code_complete/simulation_function_R1.R")
source("Scripts/Code_complete/spatial_sampling_bias_R1.R")
source("Scripts/Code_complete/species_sampling_R1.R")
source("Scripts/Code_complete/generate_virtual_species_R1.R")
source("Scripts/Code_complete/create_blocking_R1.R")
source("Scripts/Code_complete/spatial_filtering_R1.R")
#source("Scripts/Code_complete/spTgGrp_corIntensity.R")

#--- Combinations of models to run
comb <- expand.grid(
  #- True sampling bias (0 = no bias -> 1 = strong bias)
  trueBias = c("wgt1"),
  # #- Is the species niche correlated to the SSB
  # nicheSsbCor = c("none"),
  #- Number of presences
  nPres = c(100),
  #- Type of spatial blocking
  block.type = c("spatial_sytematic"),
  #- Type of SSB correction method
  ssb_meth =  "occFilter"#c( "occBuff", "occFilter", "covCond", "featDD")
)

# --- Load sampled species
#load( "Outputs/species_prev_selection.RData")

#--- Run in parallel using snowfall

#-- Set up parallel
sfInit( parallel = TRUE, cpus = 4)

#- Load libraries
sfLibrary(terra)
sfLibrary(sf)
sfLibrary(tidyverse)
sfLibrary(Rfast)
sfLibrary(gtools)
sfLibrary(mgcv)
sfLibrary(spatstat)
sfLibrary(Metrics)
sfLibrary(blockCV)
sfLibrary(raster)
sfLibrary(stats)
sfLibrary(enmSdm)
sfLibrary(gdistance)
sfLibrary(speedglm)

#- Export data
sfExport(
  list = c(
    "spp_prev_select",
    "comb",
    "run_ssb_sim",
    "create.blocks",
    "filter.occurrences",
    "species.occ",
    "samp.bias.feature.dd",
    "species.sampling",
    "samp.bias.occ.buff"
  )
)

#--- Run simulation in parallel j = 1; i = 1
sfLapply(1:nrow(comb), function(i) {
  
  #- Run for N = 100 virtual species / virtual environment j =1
  for (j in 1:nrow(splkup)) {
    
    print(i)
    print(j)
    
    # - Species and env
    spp_j <- spp_prev_select$spp[j]
    env_j <- spp_prev_select$env[j]
    fname <- paste0("vspp_", spp_j, "_env_", env_j, ".RData")
    
    # - Load species data
    spDat_j <- get(load(paste0("Outputs/virtual_species/", fname)))[[1]]
    spDat_j <- spDat_j[,c("poc", "pa")] 
    
    # - Load environmental data scenario
    enviro_j <-
      get(load(
        paste0("Outputs/simulated_environments/env_", env_j, "_R1.RData")
      ))
    
    # - Join species occurrence data to environmental
    spEnv_j <- cbind(enviro_j, spDat_j)
    
    # ggplot() +
    #   geom_raster(aes(x = X, y = Y, fill =  poc), data = spEnv_j) +
    #   geom_point(aes(x = X, y = Y), size = 0.1,colour = "red", data= spEnv_j[spEnv_j$pa == 1, ]) +
    #   scale_fill_viridis_c()
    
    #- File output name
    f_out <- paste0(
      "Outputs/simulation_res_optimization_R1/",
      "ppm_sim_",
      comb[i, 4],
      "_",
      comb[i, 1],
      "_",
      comb[i, 2],
      "_",
      comb[i, 3],
      "_",
      spp_j,
      "_",
      env_j,
      "_",
      ".RData"
    )
    
    #- Run simulation (if scenario is not already run)
    if (!file.exists(f_out)) {
      out_i <- run_ssb_sim(
        sppEnv = spEnv_j,
        nPres = comb[i, "nPres"],
        n_folds = 10,
        trueBias = comb[i, "trueBias"],
        ssb_meth = comb[i, "ssb_meth"],
        modVar = c("V1", "V2", "V3"),
        corNiche = FALSE,
        block.type = comb[i, "block.type"],
        spp = spp_j,
        env = env_j
      )
      #- Save results if they are not NULL (i.e. if the simulation ran)
      if (!is.null(out_i)) {
        save(out_i, file = f_out, compress = "xz")
      }
    }
  }
  return(NULL)
})

#- Stop sf
sfStop()

#------------------------------------------------------------------------------#
# (3) Analyse simulation results
#------------------------------------------------------------------------------#

#--- Source functions
source("Scripts/process_comparision_data_R1.R")
source("Scripts/main_figures_R1.R")

#' (Q1) Which SSB correction method and evaluation approach provides the most
#' accurate measure of independent cross-validation test score?
#' (Q2) Which SSB correction method most consistently provides gains in model 
#' predictive performance?

#--- Load in all results files

# -- Uncorrelated niches
sim_res_uncor <- prepare_sim_data(nicheCor = FALSE)
# - Save (it takes time to process)
save(sim_res_uncor, 
     file = "Outputs/Simulation_results_compiled_unCorNiche_R1.Rdata", 
     compress = "xz")
#--- Prepare best data
sim_comp_uncor <- process.comparision.data(sim_res_uncor, mod.qual = "best")
save(sim_comp_uncor, 
     file = "Outputs/Simulation_results_comp_unCorNiche_R1.Rdata", 
     compress = "xz")

# -- Correlated niches
sim_res_cor <- prepare_sim_data(nicheCor = TRUE)
# - Save (it takes time to process)
save(sim_res_cor, 
     file = "Outputs/Simulation_results_compiled_corNiche_R1.Rdata", 
     compress = "xz")
#--- Prepare best data
#load("Outputs/Simulation_results_compiled_corNiche_R1.Rdata")
sim_comp_cor <- process.comparision.data(sim_res_cor, mod.qual = "best")
save(sim_comp_cor, 
     file = "Outputs/Simulation_results_comp_corNiche_R1.Rdata", 
     compress = "xz")

#---
#-- [Q1 & Q2] Which SSB correction method and evaluation approach provides the most
#-- accurate measure of independent cross-validation test score and the 
#-- which method produced the best model performance over uncorrected models?
#---

#' Which approach provides the best solution to SSB correction as assessed using
#' internal cross-validation given the need to minimise errors relative to 
#' independent cross-validation and provide the best SSB correction possible for 
#' the data.

# ---
# --- Uncorrelated niche
# ---

sim_comp_uncor <- get(load("Outputs/Simulation_results_comp_unCorNiche_R1.Rdata")) %>% 
  filter(ssb_meth  != "spTgGrp") %>% 
  droplevels()

# - AUC (nicheCor = F)
auc_wgt1_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt1", "AUC", corNiche = FALSE)
auc_wgt08_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt08", "AUC", corNiche = FALSE)
auc_wgt06_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt06", "AUC", corNiche = FALSE)
auc_wgt0_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt0", "AUC", corNiche = FALSE)
# - boyce (nicheCor = F)
boyce_wgt1_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt1", "Boyce", corNiche = FALSE)
boyce_wgt08_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt08", "Boyce", corNiche = FALSE)
boyce_wgt06_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt06", "Boyce", corNiche = FALSE)
boyce_wgt0_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt0", "Boyce", corNiche = FALSE)
# - Spearman's correlation (nicheCor = F)
rs_wgt1_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt1", "S (Spearman's)", corNiche = FALSE)
rs_wgt08_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt08", "S (Spearman's)", corNiche = FALSE)
rs_wgt06_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt06", "S (Spearman's)", corNiche = FALSE)
rs_wgt0_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt0", "S (Spearman's)", corNiche = FALSE)
# - Pearson's correlation (nicheCor = F)
rp_wgt1_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt1", "S (Pearson's)", corNiche = FALSE)
rp_wgt08_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt08", "S (Pearson's)", corNiche = FALSE)
rp_wgt06_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt06", "S (Pearson's)", corNiche = FALSE)
rp_wgt0_uncor <- create.fig2.wrapper(sim_comp_uncor, "wgt0", "S (Pearson's)", corNiche = FALSE)


# -- Panel figure
# - AUC - Uncorrelated
auc_uc_strong <- panel_fig_wrapper(auc_wgt08_uncor[[1]], auc_wgt1_uncor[[1]], "uncor", "auc_strong")
auc_uc_weak <- panel_fig_wrapper(auc_wgt06_uncor[[1]], auc_wgt0_uncor[[1]], "uncor", "auc_weak") 
# - Boyce - Uncorrelated
boyce_uc <- panel_fig_wrapper(boyce_wgt08_uncor[[1]], boyce_wgt1_uncor[[1]], "uncor", "boyce_strong") 
boyce_uc_weak <- panel_fig_wrapper(boyce_wgt06_uncor[[1]], boyce_wgt0_uncor[[1]], "uncor", "boyce_weak") 
# - Spearman's correlation  - Uncorrelated
rs_uc <- panel_fig_wrapper(rs_wgt08_uncor[[1]], rs_wgt1_uncor[[1]], "uncor", "rs_strong") 
rs_uc_weak <- panel_fig_wrapper(rs_wgt06_uncor[[1]], rs_wgt0_uncor[[1]], "uncor", "rs_weak") 
# - Pearson's correlation  - Uncorrelated
rp_uc <- panel_fig_wrapper(rp_wgt08_uncor[[1]], rp_wgt1_uncor[[1]], "uncor", "rp_strong") 
rp_uc_weak <- panel_fig_wrapper(rp_wgt06_uncor[[1]], rp_wgt0_uncor[[1]], "uncor", "rp_weak") 

# ---
# ---  Summary statistic
# ---

# -- Proportion decrease, difference for decrease and agree int/ind

# - AUC summary effect sizes
auc_all <- 
  rbind(auc_wgt1_uncor[[2]],auc_wgt08_uncor[[2]], 
        auc_wgt06_uncor[[2]],auc_wgt0_uncor[[2]])
# - Boyce summary effect sizes
boyce_all <- 
  rbind(boyce_wgt1_uncor[[2]],boyce_wgt08_uncor[[2]],
        boyce_wgt06_uncor[[2]],boyce_wgt0_uncor[[2]])
# - Spearman's summary effect sizes
rs_all <- 
  rbind(rs_wgt1_uncor[[2]],rs_wgt08_uncor[[2]],
        rs_wgt06_uncor[[2]],rs_wgt0_uncor[[2]])
# - Pearson's summary effect sizes
rp_all <- 
  rbind(rp_wgt1_uncor[[2]],rp_wgt08_uncor[[2]],
        rp_wgt06_uncor[[2]],rp_wgt0_uncor[[2]])

# -- Summarise overall effects for AUC
auc_all %>%
  group_by(ssb_meth, ssb_true, blkType) %>%
  summarise(
    mdif_ind_m = mean(abs(mdif_ind_m), na.rm = T),
    propInc_m = mean(propIncr),
    mdif_inc = mean(mdif_ind_Incr, na.rm = T),
    propDec_m = mean(propDecr),
    mdif_dec = mean(mdif_ind_Decr, na.rm = T),
    intIndAgg = mean(intIndAgg)
  )

# --- Tables of proportion increasing / decreasing by block and ssb strength
# - AUC
prop_inc_dec_auc <- prop_inc_dec(auc_all) 
write.csv(prop_inc_dec_auc, "Outputs/Figures_R2/Table_incDec_prop_auc.csv")
# - Boyce
prop_inc_dec_boyce <- prop_inc_dec(boyce_all) 
write.csv(prop_inc_dec_boyce, "Outputs/Figures_R2/Table_incDec_prop_boyce.csv")
# - Spearman's
prop_inc_dec_rs <- prop_inc_dec(rs_all) 
write.csv(prop_inc_dec_rs, "Outputs/Figures_R2/Table_incDec_prop_rs.csv")
# - Pearson's
prop_inc_dec_rp <- prop_inc_dec(rp_all) 
write.csv(prop_inc_dec_rp, "Outputs/Figures_R2/Table_incDec_prop_rp.csv")

# --- Tables proportion agreement internal vs independent on direction effect
# - AUC
prAg_indInt_auc <- prop_agree_indInt(auc_all)
write.csv(prAg_indInt_auc, "Outputs/Figures_R2/Table_intInd_prAg_auc.csv")
# - Boyce
prAg_indInt_boyce <- prop_agree_indInt(boyce_all)
write.csv(prAg_indInt_boyce, "Outputs/Figures_R2/Table_intInd_prAg_boyce.csv")
# - Spearman's
prAg_indInt_rs <- prop_agree_indInt(rs_all) 
write.csv(prAg_indInt_rs, "Outputs/Figures_R2/Table_intInd_prAg_rs.csv")
# - Pearson's
prAg_indInt_rp <- prop_agree_indInt(rp_all) 
write.csv(prAg_indInt_rp, "Outputs/Figures_R2/Table_intInd_prAg_rp.csv")

#---
#-- Does internal cv idenfify best / avoid worst model?
#---

# --- Consistency of finding best model

sim_res_uncor <- sim_res_uncor %>%
  filter(ssb_meth != "spTgGrp") %>% 
  mutate(
    prevelance = prevelance / 62500,
    prevCat =
      cut(
        prevelance,
        breaks = seq(0, 0.6, 0.1),
        labels = seq(0.1, 0.6, 0.1)
      ))

# - Best model independent
best_int <- sim_res_uncor %>%
  filter(testDat == "int") %>%
  group_by(spId,
           metric,
           ssb_true,
           nPres,
           prevCat,
           testDat,
           ssb_meth,
           blkType) %>%
  slice(which.max(m)) %>%
  rename(m_int_b = "m") %>%
  ungroup() %>% 
  dplyr::select(
    ssb_correct,
    spId,
    prevelance,
    prevCat,
    metric,
    ssb_true,
    nPres,
    ssb_meth,
    blkType,
    m_int_b
  )
best_int <- sim_res_uncor %>%
  filter(testDat == "ind") %>%
  left_join(
    best_int,
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
      )
  ) %>%
  rename(m_ind_bInt = "m") %>%
  dplyr::select(spId,
                prevelance,
                prevCat,
                metric,
                ssb_true,
                nPres,
                ssb_meth,
                blkType,
                m_int_b,
                m_ind_bInt) 

# - Best model independent
best_ind <- sim_res_uncor %>%
  filter(testDat == "ind") %>%
  group_by(spId,
           metric,
           ssb_true,
           nPres,
           prevCat,
           testDat,
           ssb_meth,
           blkType) %>%
  slice(which.max(m)) %>%
  rename(m_ind_b = "m") %>%
  ungroup() %>% 
  dplyr::select(
    spId,
    prevelance,
    prevCat,
    metric,
    ssb_true,
    nPres,
    ssb_meth,
    blkType,
    m_ind_b
  )

# - Worst model independent
worst_ind <- sim_res_uncor %>%
  filter(testDat == "ind",
         !is.na(m)) %>%
  group_by(spId,
           metric,
           ssb_true,
           nPres,
           prevCat,
           testDat,
           ssb_meth,
           blkType) %>%
  slice(which.min(m)) %>%
  rename(m_ind_w = "m") %>%
  ungroup() %>% 
  dplyr::select(
    spId,
    prevelance,
    prevCat,
    metric,
    ssb_true,
    nPres,
    ssb_meth,
    blkType,
    m_ind_w
  )

# - Combine and then join with int
ind_uncor_mods <- full_join(best_ind, worst_ind) 
tune_uncor_all <- 
  best_int %>% 
  left_join(ind_uncor_mods) %>% 
  mutate(b = ifelse(m_ind_bInt == m_ind_b, 1, 0),
         w = ifelse(m_ind_bInt > m_ind_w, 1, 0))

# --- Summarise proportion ids best and better than worse
bw <- tune_uncor_all %>% 
  #filter(metric == "S (Pearson's)") %>% 
  group_by(metric, ssb_meth, ssb_true, blkType) %>% 
  summarise(w = round((sum(w, na.rm = T) / n()) * 100, 1),
            b = round((sum(b, na.rm = T) / n()) * 100, 1)) %>% 
  ungroup() %>%
  pivot_wider(
    names_from = c(blkType),
    values_from = c(w,b)
  ) 
bw
write.csv(bw, "Outputs/Figures_R2/Table_CV_tuning_proportions.csv")
mean(bw$b_random)
mean(bw$b_spatial_sytematic)
mean(bw$b_enviro)
mean(bw$w_random)
mean(bw$w_spatial_sytematic)
mean(bw$w_enviro)

#---
#-- Effects of nPres and Prevalence on delta[metric]
#---

# -- Robust linear mixed model

# - Combinations to run
combs <-
  expand.grid(c("AUC", "Boyce", "S (Pearson's)", "S (Spearman's)"),
              c("uncor"),
              c("random", "enviro", "spatial_sytematic"))

# -- Process uncorrelated
load("Outputs/Simulation_results_comp_unCorNiche_R1.Rdata")
sim_uncor <- sim_comp_uncor %>%
  mutate(corType = "uncor") %>%
  pivot_wider(names_from = ssb_correct,
              values_from = c(m, sd, N)) %>%
  mutate(mdif = (m_Yes - m_No)) %>%
  dplyr::select(spId,
                metric,
                ssb_true,
                nPres,
                ssb_meth,
                blkType,
                prevelance,
                testDat,
                corType,
                mdif,
                m_No) %>% 
  mutate(nPres = as.numeric(as.character(nPres))) %>% 
  filter(ssb_true %in% c("wgt0","wgt06","wgt08","wgt1"))

# - Set up parallel
sfInit( parallel = TRUE, cpus = 6)

# - Load libraries
sfLibrary(robustlmm)

# - Export data
sfExport(
  list = c(
    "sim_res",
    "combs"
  )
)

# - Run model
sfLapply(1:nrow(combs), function(i) {
  sim_res_b_i <- 
    sim_res[sim_res$metric == combs[i,1] & sim_res$corType == combs[i,2] &
              sim_res$blkType == combs[i,3],]
  mod_i <-
    rlmer(
      mdif ~
        #m_No +
        ssb_true +
        ssb_meth +
        nPres +
        prevelance +
        ssb_meth:nPres +
        ssb_meth:prevelance +
        nPres:prevelance +
        (1 | spId),
      data = sim_res_b_i
    )
  if(combs[i,1] == "AUC") {lab_i <- "auc"}
  if(combs[i,1] == "S (Pearson's)") {lab_i <- "rp"}
  if(combs[i,1] == "S (Spearman's)") {lab_i <- "rs"}
  if(combs[i,1] == "Boyce") {lab_i <- "boyce"}
  corType <- combs[i,2]
  blkType <- combs[i,3]
  save(
    mod_i,
    file = paste0(
      "Outputs/models_R1/mod_typePerformance_",
      lab_i,
      "_",
      corType,
      "_",
      blkType,
      ".RData"
    ),
    compress = "xz"
  )
})

# --- 
# -- Marginal effects figures
# ---

# - AUC Random
mod_auc_r <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_auc_uncor_random.RData")))
p_auc_r <- plot_partial_effects(mod_auc_r)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_auc_r.png",
  plot = p_auc_r,
  width = 8,
  height = 12
)

# - AUC Spatial systematic
mod_auc_s <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_auc_uncor_spatial_sytematic.RData")))
p_auc_s <- plot_partial_effects(mod_auc_s)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_auc_s.png",
  plot = p_auc_s,
  width = 8,
  height = 12
)

# - AUC environmental
mod_auc_e <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_auc_uncor_enviro.RData")))
p_auc_e <- plot_partial_effects(mod_auc_e)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_auc_e.png",
  plot = p_auc_e,
  width = 8,
  height = 12
)

# - boyce Random
mod_boyce_r <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_boyce_uncor_random.RData")))
p_boyce_r <- plot_partial_effects(mod_boyce_r)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_boyce_r.png",
  plot = p_boyce_r,
  width = 8,
  height = 12
)

# - boyce Spatial systematic
mod_boyce_s <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_boyce_uncor_spatial_sytematic.RData")))
p_boyce_s <- plot_partial_effects(mod_boyce_s)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_boyce_s.png",
  plot = p_boyce_s,
  width = 8,
  height = 12
)

# - Boyce environmental
mod_boyce_e <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_boyce_uncor_enviro.RData")))
p_boyce_e <- plot_partial_effects(mod_boyce_e)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_boyce_e.png",
  plot = p_boyce_e,
  width = 8,
  height = 12
)


# - rs Random
mod_rs_r <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_rs_uncor_random.RData")))
p_rs_r <- plot_partial_effects(mod_rs_r)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_rs_r.png",
  plot = p_rs_r,
  width = 8,
  height = 12
)

# - rs Spatial systematic
mod_rs_s <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_rs_uncor_spatial_sytematic.RData")))
p_rs_s <- plot_partial_effects(mod_rs_s)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_rs_s.png",
  plot = p_rs_s,
  width = 8,
  height = 12
)

# - rs environmental
mod_rs_e <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_rs_uncor_enviro.RData")))
p_rs_e <- plot_partial_effects(mod_rs_e)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_rs_e.png",
  plot = p_rs_e,
  width = 8,
  height = 12
)

# - rp Random
mod_rp_r <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_rp_uncor_random.RData")))
p_rp_r <- plot_partial_effects(mod_rp_r)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_rp_r.png",
  plot = p_rp_r,
  width = 8,
  height = 12
)

# - rp Spatial systematic
mod_rp_s <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_rp_uncor_spatial_sytematic.RData")))
p_rp_s <- plot_partial_effects(mod_rp_s)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_rp_s.png",
  plot = p_rp_s,
  width = 8,
  height = 12
)

# - rp environmental
mod_rp_e <-
  get(load(paste0("Outputs/models_R1/mod_typePerformance_rp_uncor_enviro.RData")))
p_rp_e <- plot_partial_effects(mod_rp_e)
ggsave(
  file = "Outputs/Figures_R2/marEff_p_rp_e.png",
  plot = p_rp_e,
  width = 8,
  height = 12
)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------END---------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#