#' Simulation of virtual species distribution, sampling and modelling
#'
#' @param enviro data.frame of virtual environment data.
#' @param nPres Numeric specifying the number of presences.
#' @param n_folds Numeric specifying the number of cross-validation folds.
#' @param trueBias Character string specifying the amount of spatial sampling
#' bias ("wgt06", "wgt08", "wgt1", "wgt0").
#' @param ssb_meth Character string indicating the type of spatial sampling bias. 
#' One of "spTgGrp", "featDD", "occBuff", "occFilter", or "covCond".
#' @param block.type A character string indicating the type of blocking to be 
#' used in model evaluation. One of "random", "spatial_sytematic", 
#' "spatial_dist", or "enviro".
#' @param detProb A numeric indicating the detection probability for the virtual 
#' species. Not for use in this analysis and should be left at 1. 
#' @param falPosRat A numeric indicating the false positive rate for the virtual 
#' species. Not for use in this analysis and should be left at 0. 
#' @param listInclusion A numeric indicating the probability of being included in 
#' a species list. Not for use in this analysis and should be left at 1. 
#' @param effort A numeric indicating the sampling effort. Not for use in this 
#' analysis and should be left at 1. 
#' @param pca_axes A numeric indicating the number of pca axes to use when 
#' creating virtual species.
#' @param modVar A vector of character strings indicating the names of the 
#' variables to be used as predictors in the SDMs.
#' @param env_i A numeric giving the environmental data ID number.
#' 
#' @return list 
#' 
#' @export
# 
run_ssb_sim <-
  function(enviro,
           nPres = 100,
           n_folds = 10,
           trueBias = c("wgt06", "wgt08", "wgt1", "wgt0"),
           ssb_meth = c("spTgGrp", "featDD", "occBuff", "occFilter", "covCond"),
           block.type = c("random", "spatial_sytematic", "spatial_dist", "enviro"),
           detProb = 1,
           falPosRat = 0,
           listInclusion = 1,
           effort = 1,
           pca_axes = 3,
           modVar = c("V1", "V2", "V4"),
           env_i) {
    #--- Create species using PCA species approach
    sp_in <-
      create.virtual.species(enviro[, c("V1", "V2", "V3", "V4")], pca_axes)
    sp_in <- cbind(enviro, sp_in)
    
    #--- Adjust the road dd bias (wrong sign)
    sp_in$pathdd <- 1 - (sp_in$pathdd)
    
    #--- Modelling section
    if (sum(sp_in$pa) > 200) {    # Don't run for very small range species
      
      #- Create SSB based on "feature" distance decay
      sp_in$biasVar <- sp_in$pathdd
      samp_bias <-
        samp.bias.feature.dd(sp_in[, c("id", "X", "Y", "biasVar")])
      #' This creates a column for each of wgt0, wgt06, wgt08, wgt1 giving the
      #' probability of sampling for four different strengths of SSB, with wgt0
      #' meaning no bias and wgt1 indicating a strong bias.
      
      #- Correlation weights for spTgGrp
      if (agrepl(ssb_meth, pattern = "spTgGrp")) {
        cor_wgh_df <-
          spTgGrp.corIntensity(sp_in[, c("id", "X", "Y", "poc")])
      }
      
      #- Names of trials (i.e. coding different strengths of SSB correction)
      trials <-
        data.frame(
          ssb_meth = c("spTgGrp", "featDD", "occBuff", "occFilter", "covCond"),
          rep1 = c("wgt06", "wgt06", 1000000, 2, 1),
          rep2 = c("wgt08", "wgt08", 500000, 5, 2),
          rep3 = c("wgt1", "wgt1",   250000, 10, 3),
          rep4 = c("wgt0", "wgt0", 0, 0, 0)
        )
      trials <- trials[trials$ssb_meth == ssb_meth, 2:5]
      
      #- Run multiple correction strengths for cv x = trials[1,3]
      rep_all <- lapply(trials, function(x) {
        print(x)
        
        #- Sample species records
        spSamp <- species.sampling(
          sp_dat = sp_in,
          bias = samp_bias[[trueBias]], # This get the true bias
          nPres = nPres,
          detProb = detProb,
          fpr = falPosRat,
          listProb = listInclusion
        )
  
        #- Summary data
        n_sites <- length(unique(spSamp$id))
        n_visits <- nrow(spSamp)
        n_det <- sum(spSamp$det)
        n_occ <- sum(spSamp$pa)
        
        #- Group vars
        groupVars <-
          lapply(c('id', 'X', 'Y', 'poc', 'pathdd', modVar), as.symbol)
        spDat <- spSamp %>%
          group_by(!!!groupVars) %>%
          summarise(det = max(det)) %>%
          filter(det == 1) %>%
          as.data.frame()
        
        #- Select covariate polynomial from 1 to 3
        if (ssb_meth == "covCond" & x != 0) {
          modVar <- c(modVar, paste0("poly(pathdd,", x, ")"))
        }
        
        #- Occurrence filtering with multiple filter distances
        if (ssb_meth == "occFilter" & x != "wgt0") {
          spDat <- filter.occurrences(enviro, spDat, as.numeric(x))
        }
        
        #- Create occurrence buffers
        samp_bias_x <- NULL
        if (ssb_meth == "occBuff" & x != 0) {
            samp_bias_x <-
              samp.bias.occ.buff(spDat[, c("X", "Y", "det")],
                                 enviro[, c("X", "Y", "id")],
                                 buf_w = as.numeric(x))
        }
        
        #- Create feature distance decay 
        if (ssb_meth == "featDD" & x != "wgt0") {
          samp_bias_x <- samp_bias[[x]]
        }
        
        #- Target group
        if (agrepl(ssb_meth, pattern = "spTgGrp")) {
          cor_wgh <- cor_wgh_df[[x]] 
          richness <- rpois(nrow(sp_in), 20 * cor_wgh)
          samp_bias_x <-  richness * samp_bias[[trueBias]]
          samp_bias_x <- scales::rescale(samp_bias_x)
        }
       
        #- Select background sample
        bkgd <- sp_in[sample(
          1:nrow(sp_in),
          size = 10000,
          replace = TRUE,
          prob = samp_bias_x
        ),]
        bkgd <- dplyr::select(bkgd, !!!groupVars)
        bkgd$det <- 0
        
        #- Combine
        spBkDat <- rbind(spDat, bkgd)

        # Weight
        detN <- nrow(spBkDat[spBkDat$det == 1,])
        bkgN <- nrow(spBkDat[spBkDat$det == 0,])
        wgt <- c(rep(1, detN), rep(detN / bkgN, bkgN))
        
        # Number of data points used to fit model
        n_mod_fit <- nrow(spBkDat)
        
        # Create blocking n_folds = 10
        spBkDat <-
          create.blocks(spBkDat, sp_in, block.type, n_folds)
        
        # Output data.frame
        AUC.int <- rep(NA, n_folds)
        corP.int <- rep(NA, n_folds)
        corS.int <- rep(NA, n_folds)
        boyce.int <- rep(NA, n_folds)
        AUC.ind <- rep(NA, n_folds)
        corP.ind <- rep(NA, n_folds)
        corS.ind <- rep(NA, n_folds)
        boyce.ind <- rep(NA, n_folds)
        
        # SDM fit using down-weighted Poisson regression
        blocks <- unique(spBkDat$block)
        blocks <- blocks[blocks != 0]
        for (blk in  blocks) {
          train.blk <- spBkDat[spBkDat$block != blk, ]
          test.blk <- spBkDat[spBkDat$block == blk, ]
          test.blk.ind <-
            sp_in[sample(1:nrow(sp_in), nrow(test.blk), replace = FALSE), ]
          
          #- Assign weightings
          p.wt <- rep(1.e-6, nrow(train.blk))
          p.wt[train.blk$det == 0] <-
            nrow(enviro) / sum(train.blk$det)
          train.blk$det <- train.blk$det / p.wt
          
          #- Select covariate data
          if (ssb_meth == "covCond" & x != 0) {
            modVar <- c(modVar, paste0("poly(pathdd,", x, ")"))
          }
          
          #- Model structure
          modForm <-
            as.formula(paste("det~", paste(modVar, collapse = "+")))
          
          modfit <- speedglm(
            modForm,
            family = poisson(),
            data = train.blk,
            weights = p.wt
          )
          
          # Condition on pathdd mean for covariate SSB correction method
          if (ssb_meth == "covCond") {
            test.blk$pathdd <- mean(spBkDat$pathdd)
            test.blk.ind$pathdd <- mean(spBkDat$pathdd)
          }
          
          # Internal CV
          fit.int <-
            predict(modfit, newdata = test.blk, type = "response")
          
          #- AUC
          AUC.int[blk] <- Metrics::auc(test.blk$det, fit.int)
          corP.int[blk] <- cor(x = test.blk$det,
                               y = fit.int,
                               method = "pearson")[[1]]
          corS.int[blk] <- cor(x = test.blk$det,
                               y = fit.int,
                               method = "spearman")[[1]]
          boyce.int[blk] <-
            enmSdm::contBoyce(fit.int[test.blk$det == 1],
                              fit.int[test.blk$det == 0])
          
          # Independent CV
          fit.ind <-
            predict(modfit, newdata = test.blk.ind, type = "response")
          #- AUC
          AUC.ind[blk] <- Metrics::auc(test.blk.ind$pa, fit.ind)
          corP.ind[blk] <- cor(x = test.blk.ind$pa,
                               y = fit.ind,
                               method = "pearson")[[1]]
          corS.ind[blk] <- cor(x = test.blk.ind$pa,
                               y = fit.ind,
                               method = "spearman")[[1]]
          boyce.ind[blk] <-
            enmSdm::contBoyce(fit.ind[test.blk.ind$pa == 1],
                              fit.ind[test.blk.ind$pa == 0])
          
        }
        
        resOut <- do.call(
          rbind,
          list(
            AUC.int,
            corP.int,
            corS.int,
            boyce.int,
            AUC.ind,
            corP.ind,
            corS.ind,
            boyce.ind
          )
        )
        resOut <- as.data.frame(resOut)
        resOut$metric <- c(
          "AUC.int",
          "corP.int",
          "corS.int",
          "boyce.int",
          "AUC.ind",
          "corP.ind",
          "corS.ind",
          "boyce.ind"
        )
        resOut$ssb_correct <- x
        resOut$ssb_true <- trueBias
        resOut$nicheSsbCor <- nicheSsbCor
        resOut$nPres <- nPres
        resOut$nPresFit <- sum(spBkDat$det)
        resOut$ssb_meth <- ssb_meth
        resOut
        
      })
      rep_all <- do.call(rbind, rep_all)
      
      rep_all$spId <- env_i
      rep_all$blkType = block.type
      rep_all$prevelance <- sum(sp_in$pa)
      rep_all
      
    }
  }
