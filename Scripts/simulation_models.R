#------------------------------------------------------------------------------#
#
#   Simulation of the spatial sampling bias analysis
#
#------------------------------------------------------------------------------#

#' 1) Simulate 100 virtual environments.
#' 2) Simulate virtual species on these environments, simulation sampling and
#'    run Poisson Point Process models using different SSB correction methods. 
#'    Evaluate on withheld and independent data (i.e. sampled strategically).
#' 3) 

#------------------------------------------------------------------------------#
# (1) Simulate 100 environments (we do this first as this is slow(ish))
#------------------------------------------------------------------------------#

#--- Functions
source("Scripts/generate_environment.R")
source("Scripts/create_enviro_features.R")

#--- Run N = 100
#' Environment is created and then the variables altered to create more 
#' realistic variables
#' V1 is given a stronger gradient
#' V2 is given an upper 'elevation' cut off
#' V3 is given a 'latitude' gradient
#' V4 is negatively correlated with distance from 'road'
#' 'Roads' are created through the landscape
for(i in 1:100) {
  
  # Create environmental data
  enviro <- generate.enviro(n = 50, phi = c(0.1, 1, 0.25, 0.25))
  
  # Create a stronger gradient
  enviro$V1 <- round(1 / (1 + exp((enviro$V1 - 0.5) / -0.1)), 3)
  
  # Create upper elevation cut off for species
  enviro$V2 <-
    round(enviro$V2 * (1 / (1 + exp((enviro$V1 - 0.5) / 0.1
    ))), 3)
  
  # Create gradient (e.g. latitude)
  enviro$V3 <- enviro$V3 * scales::rescale(enviro$X * enviro$Y)
  
  # Create linear feature through "elevation"
  l.feat <- create.linear.features(enviro[c("X", "Y", "V1")])
  
  # Create area features (e.g. PAs) in reference to "forest"
  a.feat <- create.area.features(enviro[c("X", "Y", "V2")])
  
  # Combine
  enviro <-
    Reduce(function(...)
      left_join(...), list(enviro, l.feat, a.feat))
  
  # Make V3 negatively correlated with road features
  enviro$V4 <-
    round(enviro$V4 * (1 / (1 + exp((enviro$pathdd - 0.1) / -0.1
    ))), 3)
  
  save(
    enviro,
    file = paste0("Outputs/simulated_environments/env_", i, ".RData"),
    compress = "xz"
  )
  
}

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
source("Scripts/simulation_function.R")
source("Scripts/spatial_sampling_bias.R")
source("Scripts/species_sampling.R")
source("Scripts/generate_virtual_species.R")
source("Scripts/create_blocking.R")
source("Scripts/spatial_filtering.R")
source("Scripts/spTgGrp_corIntensity.R")

#--- Combinations of models to run
comb <- expand.grid(
  #- True sampling bias (0 = no bias -> 1 = strong bias)
  trueBias = c("wgt0", "wgt06", "wgt08", "wgt1"),
  #- Is the species niche correlated to the SSB
  nicheSsbCor = c("none"),
  #- Number of presences
  nPres = c(50, 200),
  #- Type of spatial blocking
  block.type = c("random", "spatial_sytematic", "enviro"),
  #- Type of SSB correction method
  ssb_meth = c("spTgGrp", "featDD", "occBuff", "occFilter", "covCond")
)


#--- Run in parallel using snowfall

#-- Set up parallel
sfInit( parallel = TRUE, cpus = 7 )

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
    "comb",
    "run_ssb_sim",
    "create.blocks",
    "filter.occurrences",
    "species.occ",
    "create.virtual.species",
    "samp.bias.feature.dd",
    "species.sampling",
    "samp.bias.occ.buff",
    "spTgGrp.corIntensity"
  )
)
j = 1
#--- Run simulation in parallel
sfLapply(1:nrow(comb), function(i) {
  
  #- Run for N = 100 virtual species / virtual environment
  for (j in 1:100) {
    
    #- Load environmental data scenario
    enviro_j <-
      get(load(
        paste0("Outputs/simulated_environments/env_", j, ".RData")
      ))
    
    #- File output name
    f_name <- paste0(
      "Outputs/simulation_res_optimization/",
      "ppm_sim_",
      comb[i, 4],
      "_",
      comb[i, 1],
      "_",
      comb[i, 2],
      "_",
      comb[i, 3],
      "_",
      j,
      "_",
      comb[i, 5],
      ".RData"
    )
    
    #- Run simulation (if scenario is not already run)
    if (!file.exists(f_name)) {
      out_i <- run_ssb_sim(
        enviro = enviro_j,
        nPres = comb[i, 3],
        n_folds = 10,
        trueBias = comb[i, 1],
        nicheSsbCor = comb[i, 2],
        ssb_meth = comb[i, 5],
        modVar = c("V1", "V2", "V4"),
        block.type = comb[i, 4],
        env_i = j
      )
      #- Save results if they are not NULL (i.e. if the simulation ran)
      if (!is.null(out_i)) {
        save(out_i, file = f_name, compress = "xz")
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
source("Scripts/process_comparision_data.R")
source("Scripts/main_figures.R")

#' (Q1) Which SSB correction method and evaluation approach provides the most
#' accurate measure of independent cross-validation test score?
#' (Q2) Which SSB correction method most consistently provides gains in model 
#' predictive performance?
#' (Q3) How does optimisation affect predictive performance by SSB correction?

#--- Load in all results files
sim_res <-
  list.files("Outputs/simulation_res_optimization/", full.names = TRUE) %>%
  lapply(., function(x) {
    get(load(x))
  }) %>%
  do.call(rbind, .) %>%
  filter(blkType %in% c("random", "spatial_sytematic", "enviro"),
         nPres %in% c(50, 200))

#--- Process factor levels and names
sim_res <- sim_res %>%
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
    nicheSsbCor = factor(nicheSsbCor,
                         levels = c("none",
                                    "partial",
                                    "full")),
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
                      levels = c("random",
                                 "spatial_sytematic",
                                 "enviro")),
    nPres = factor(nPres,
                   levels = c("50",
                              "200")),
    
  ) %>%
  pivot_longer(cols = starts_with("V")) # Pivot to a longer table format

#--- Save (it takes time to process)
save(sim_res, 
     file = "Outputs/Simulation_results_compiled_all.Rdata", 
     compress = "xz")
#load("Outputs/Simulation_results_compiled_all.Rdata")

#---
#-- [Q1 & Q2] Which SSB correction method and evaluation approach provides the most
#-- accurate measure of independent cross-validation test score and the 
#-- which method produced the best model performance over uncorrected models?
#---

#' Which approach provides the best solution to SSB correction as assessed using
#' internal cross-validation given the need to minimise errors relative to 
#' independent cross-validation and provide the best SSB correction possible for 
#' the data.

#--- Real SSB = zero (wgt0)

pp_comp_0 <- err_vs_pp_dat(sim_res, "wgt0", "best")

#- Fig 3 - SSB 0
ssb0_fig3 <- err_vs_pp_plot(pp_comp_0, "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_3_ssb0.png",
  plot = ssb0_fig3,
  width = 9,
  height = 8
)

#- Fig S3 - SSB 0
ssb0_figS3 <- err_vs_pp_plot(pp_comp_0, "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S3_ssb0.png",
  plot = ssb0_figS3,
  width = 9,
  height = 8
)

#--- Real SSB = weak (wgt06)
pp_comp_06 <- err_vs_pp_dat(sim_res, "wgt06", "best")

#- Fig 3 - SSB 06
ssb06_fig3 <- err_vs_pp_plot(pp_comp_06, "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_3_ssb06.png", 
  plot = ssb06_fig3, width = 9, height = 8
)

#- Fig S3 - SSB 06
ssb06_figS3 <- err_vs_pp_plot(pp_comp_06, "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S3_ssb06.png", 
  plot = ssb06_figS3, width = 9, height = 8
)

#--- Real SSB = moderate (wgt08)

pp_comp_08 <- err_vs_pp_dat(sim_res, "wgt08", "best")

#-Summarise AUC
pp_comp_08 %>% 
  filter(metric == "S (Spearman's)",
         blkType == "Random") %>%
  print(n =20) %>%
  dplyr::select(nPres, ssb_meth, blkType, gm1, gsd1, pSig_intvInd )

#- Fig 3 - SSB 08
ssb08_fig3 <- err_vs_pp_plot(pp_comp_08, "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_3_ssb08.png", 
  plot = ssb08_fig3, width = 9, height = 8
)

#- Fig S3 - SSB 08
ssb08_figS3 <- err_vs_pp_plot(pp_comp_08, "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S3_ssb08.png", 
  plot = ssb08_figS3, width = 9, height = 8
)

#--- Real SSB = strong (wgt1)
pp_comp_1 <- err_vs_pp_dat(sim_res, "wgt1", "best")

#- Fig 3 - SSB 1
ssb1_fig3 <- err_vs_pp_plot(pp_comp_1, "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_3_ssb1.png", 
  plot = ssb1_fig3, width = 9, height = 8
)

#- Fig S3 - SSB 1
ssb1_figS3 <- err_vs_pp_plot(pp_comp_1, "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S3_ssb1.png", 
  plot = ssb1_figS3, width = 9, height = 8
)

#---
#-- Effect of model optimisation via cross-validation, comparing corrected and 
#-- uncorrected models
#---

#--- Real SSB = zero (wgt0)
#- Fig 4 - SSB 0
ssb0_fig4 <- cv_mod_improve(sim_res, "wgt0", "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_4_ssb0.png", 
  plot = ssb0_fig4, width = 9, height = 8
)
#- Fig S4 - SSB 0
ssb0_figS4 <- cv_mod_improve(sim_res, "wgt0", "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S4_ssb0.png", 
  plot = ssb0_figS4, width = 9, height = 8
)

#--- Real SSB = weak (wgt06)
#- Fig 4 - SSB 06
ssb06_fig4 <- cv_mod_improve(sim_res, "wgt06", "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_4_ssb06.png", 
  plot = ssb06_fig4, width = 9, height = 8
)
#- Fig S4 - SSB 06
ssb06_figS4 <- cv_mod_improve(sim_res, "wgt06", "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S4_ssb06.png", 
  plot = ssb06_figS4, width = 9, height = 8
)

#--- Real SSB = moderate (wgt08)
#- Fig 4 - SSB 08
ssb08_fig4 <- cv_mod_improve(sim_res, "wgt08", "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_4_ssb08.png", 
  plot = ssb08_fig4, width = 9, height = 8
)
#- Fig S4 - SSB 08
ssb08_figS4 <- cv_mod_improve(sim_res, "wgt08", "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S4_ssb08.png", 
  plot = ssb08_figS4, width = 9, height = 8
)

#--- Real SSB = strong (wgt1)
#- Fig 4 - SSB 1
ssb1_fig4 <- cv_mod_improve(sim_res, "wgt1", "AUC", "S (Spearman's)")
ggsave(
  filename = "Outputs/Figures/Figure_4_ssb1.png", 
  plot = ssb1_fig4, width = 9, height = 8
)
#- Fig S4 - SSB 108
ssb1_figS4 <- cv_mod_improve(sim_res, "wgt1", "Boyce", "S (Pearson's)")
ggsave(
  filename = "Outputs/Figures/Figure_S4_ssb1.png", 
  plot = ssb1_figS4, width = 9, height = 8
)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------END---------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
