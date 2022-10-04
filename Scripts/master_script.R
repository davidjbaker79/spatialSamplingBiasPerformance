#------------------------------------------------------------------------------#
#
#   Identifying effective strategies for correcting spatial sampling bias in 
#   species distribution models without independent test data 
#
#   By David Baker (Mar 2022)
#
#------------------------------------------------------------------------------#

#' This project aims to determine reliable methods for correcting for Spatial
#' sampling bias (SSB) in presence-only species data using presence-background 
#' species distribution models (SDMs).
#' 
#' The first part of the analysis focuses on previously published comparisons
#' of model performance extracted from a systematic review of the literature.
#' 
#' The meta-analysis of these results leads onto a simulation model to determine
#' which methods of SSB correction and cross-validation procedure provide the 
#' most reliable measure of true model performance and which methods produce the
#' best predictive performance.

#------------------------------------------------------------------------------#
# Setup
#------------------------------------------------------------------------------#

#--- Increase memory size
memory.size(1000000)

#--- Set seed
set.seed(17)

#--- Load packages

#- Workflow
library(tidyverse)

#-# Data handling 
library(readxl)
library(bib2df)
library(gtools)
library(dplyr)
library(scales)

#- Plotting
library(patchwork)
library(ggthemes)
library(colorBlindness)

#- Modelling
library(brms)
library(mgcv)
library(speedglm)
library(Metrics)
library(ecospat)
library(blockCV)
#remotes::install_github('adamlilith/enmSdm', dependencies=TRUE)
library(enmSdm)

#- Computing 
library(snowfall)

#- Spatial
library(terra)
library(sf)
library(Rfast)
library(gdistance)
library(raster)


#------------------------------------------------------------------------------#
# Analysis 
#------------------------------------------------------------------------------#

#' This analysis is carried out in two main scripts (with associated functions):

#- 1) Meta-analysis of systematic review data
#' Run meta-analysis script "meta_analysis.R"
#' This produces statistics for systematic review, fig 1, fig 2, fig S1, and 
#' inferential models (Table 2).

#- 2) Simulation models
#' Run simulation models "simulation_models.R"
#' This produces simulation results and analyses the outputs, producing fig 3 
#' and fig 4, and supplementary figures

#------------------------------------------------------------------------------#
#-------------------------------------END--------------------------------------#
#------------------------------------------------------------------------------#