#' Estimate Changes in Precipitation Using CMIP Model Scenarios
#'
#' This function calculates changes in precipitation using CMIP model scenarios, comparing future projections with current conditions.
#'
#' @param ROI A `SpatVector` representing the region of interest.
#' @param res Numeric; the resolution of the data. Possible values are `2.5`, `5`, or `10`.
#' @param ssp Character; the Shared Socio-economic Pathway code. Options are `"126"`, `"245"`, `"370"`, or `"585"`.
#' @param model Character; the climate model to use. See details at \url{https://www.worldclim.org/data/cmip6/cmip6climate.html}.
#' @param year Character; the future time period to evaluate. Must be one of `"2021-2040"`, `"2041-2060"`, or `"2061-2080"`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{abs_chang}{A raster representing the absolute change in precipitation (mm).}
#'   \item{rel_change}{A raster representing the relative change in precipitation, expressed as a fraction of the current precipitation.}
#' }
#'
#' @details
#' This function downloads current and future precipitation data from the WorldClim database and CMIP6 projections.
#' The data is clipped to the provided region of interest (`ROI`), and annual changes in precipitation are calculated.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' ROI <- vect("path_to_shapefile.shp")
#' result <- myfunction_Prec_Change(ROI, res = 5, ssp = "245", model = "BCC-CSM2-MR", year = "2041-2060")
#' plot(result$abs_chang)
#' plot(result$rel_change)
#' }
myfunction_Prec_Change <- function(ROI, res, ssp, model, year) {
  require(geodata)
  require(terra)
  #### download data
  current <- worldclim_global(var = "prec", res = res)
  future <- cmip6_world(model, ssp, year, var = "prec", res = res)
  
  #### clip to ROI
  current <- current %>%
    terra::crop(ROI) %>%
    terra::mask(ROI)
  future <- future %>%
    terra::crop(ROI) %>%
    terra::mask(ROI)
  
  #### Calculate annual change
  delta_annual <- sum(future - current)
  delta_prec <- sum(future - current) / sum(current)
  
  return(list(abs_chang = delta_annual, rel_change = delta_prec))
}