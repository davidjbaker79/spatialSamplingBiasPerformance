#' Spatial thinning of occurrence data
#'
#' @param x A data.frame produced by \code{species.sampling} from within
#' \code{run_exp_i}
#' @param y A data.frame with species' detection.
#' @param fdis A numeric giving the degree of filtering.
#'
#' @import data.table
#' @import dplyr
#'
#' @return Dataframe with spatial adjusted data
#' 
#' @export
filter.occurrences <- function(x, y, fdis) {
  x_i <- rast(x[, c("X", "Y", "id")], type = "xyz", crs = "EPSG:4326")
  x_i <- setValues(x_i, NA)
  x_i[y$id] <- 1
  x_i <- terra::aggregate(x_i, fdis, sum, na.rm = TRUE)
  x_i <- disagg(x_i, fdis)
  wdf <- terra::as.data.frame(x_i, cells = TRUE)
  wdf$wgt <-  1 / wdf$id
  names(wdf) <- c("id", "N", "wgt")
  wdf <- left_join(wdf, y)
  wdf <- na.omit(wdf)
  #- Thin data
  wdf <- wdf[rbinom(1:nrow(wdf), size = 1, prob = wdf$wgt) == 1,]
  wdf <- wdf[, names(y)]
  # x_ii <- rast(wdf[, c("X", "Y", "det")], type = "xyz",
  # crs = "EPSG:4326")
  # plot(x_ii)

}
