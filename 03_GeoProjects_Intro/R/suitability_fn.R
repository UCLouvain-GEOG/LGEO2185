## ---- Suitability --------
#' Suitability Map Modeling
#'
#' This set of functions models species habitat suitability using environmental and occurrence data.
#' It predicts suitability for present and future climatic conditions.
#'
#' @param region A spatial object defining the region of interest.
#' @param genus The genus of the species.
#' @param species The species name.
#' @param res Resolution of environmental data. Must be one of 2.5, 5, 10
#' @param ssp Shared Socioeconomic Pathway scenario. Must be one of '126', '245', '370', '585'
#' @param model Climatic model for future prediction. See \code{geodata::cmip6_world()}
#' @param time A character vector of time periods. Must be one or several of '1970-2000', '2021-2040', '2041-2060', or '2061-2080'
#' @return A list containing suitability maps and input_data.
#' @import geodata terra leaflet leafsync
#' @examples
#' library("terra")
#' library("geodata")
#' library("here")
#' options(geodata_default_path = here())
#' region <- gadm(country = "MEX", level = 0, path = here())
#' output <- suitability_map(
#'     region = region, genus = "ambystoma",
#'     species = "mexicanum", res = 10,
#'     model = "ACCESS-CM2", ssp = "126"
#' )
#' r <- output$suitability_map
#' r
#' @importFrom assertthat assert_that
#' @export
suitability_map <- function(region, genus, species, res, ssp, model, time = c("1970-2000", "2021-2040", "2041-2060", "2061-2080")) {
    ### Validate inputs ###
    assertthat::assert_that(!is.null(getOption("geodata_default_path")), msg = "Setup your geodata folder using: options(geodata_default_path = 'mypath')")

    # STEP 1
    predictors <- get_predictors(region, res, ssp, model, time)
    # STEP 2
    occ <- get_occurrence(region, genus, species)
    # STEP 3
    occurrence_df <- create_presence_absence_df(occ, predictors[1])
    # STEP 4
    model_fit <- fit_suitability_model(occurrence_df)
    # STEP 5
    maps <- predict_suitability(model_fit, predictors)
    return(list(
        "suitability_map" = maps,
        "input_data" = list(
            "occurrence" = occ, "predictors" = predictors,
            "model_fit" = model_fit, "sampled_df" = occurrence_df
        )
    ))
}

#' Download and prepare occurrence data
#'
#' @param region The region of interest, a \code{SpatVector}
#' @param genus The genus of the species, a \code{character}
#' @param species The species name, \code{character}
#' @param path output path
#' @return A spatial object with occurrence data.
#' @import geodata terra
#' @importFrom assertthat assert_that
#' @export
get_occurrence <- function(region, genus, species, path = getOption("geodata_default_path")) {
    assertthat::assert_that(!is.null(path), msg = "Setup your geodata folder using: options(geodata_default_path = 'mypath')")
    assertthat::assert_that(is.character(genus), msg = "Genus must be a character string")
    assertthat::assert_that(is.character(species), msg = "Species must be a character string")
    assertthat::assert_that(inherits(region, "SpatVector"), msg = "Region must be a 'SpatVector' object.")

    occ <- geodata::sp_occurrence(genus = genus, species = species, ext = region, sp = TRUE, download = TRUE, geo = TRUE, path = path)
    occ <- terra::vect(occ)
    crs(occ) <- terra::crs(region)

    if (length(occ) < 40) {
        stop("Not enough occurrences (ie < 40) in the selected region to model the suitability map.")
    }

    occ$pres <- 1
    return(occ)
}

#' Create environmental raster stack
#'
#' @param region The region of interest, a \code{SpatVector}
#' @param res Resolution of environmental data.
#' @param ssp Shared Socioeconomic Pathway scenario. Must be one of '126', '245', '370', '585'
#' @param model Climatic model for future prediction.
#' @param time A vector of time period. One or several of "1970-2000", "2021-2040", "2041-2060", or "2061-2080"
#' @param path output path
#' @return A list of environmental raster layers for different time periods.
#' @import terra
#' @importFrom geodata worldclim_global cmip6_world
#' @importFrom assertthat assert_that
#' @export
get_predictors <- function(region, res, ssp, model, time = c("1970-2000", "2021-2040", "2041-2060", "2061-2080"), path = getOption("geodata_default_path")) {
    assertthat::assert_that(!is.null(path), msg = "Setup your geodata folder using: options(geodata_default_path = 'mypath')")
    assertthat::assert_that(res %in% c(2.5, 5, 10), msg = "Resolution must be 2.5, 5, or 10")
    assertthat::assert_that(ssp %in% c("126", "245", "370", "585"), msg = "ssp must be '126', '245', '370', '585'")
    assertthat::assert_that(all(time %in% c("1970-2000", "2021-2040", "2041-2060", "2061-2080")), msg = "time should be one of '1970-2000', '2021-2040', '2041-2060', or '2061-2080'")
    assertthat::assert_that(inherits(region, "SpatVector"), msg = "Region must be a 'SpatVector' object.")

    # Retrieve all climate data
    climate_data <- purrr::map(time, function(year) {
        if (year == "1970-2000") {
            x <- geodata::worldclim_global(var = "bio", res = res, ssp = ssp, path = path)
            names(x) <- paste0("bio_", 1:nlyr(x))
            x
        } else {
            x <- geodata::cmip6_world(model = model, ssp = ssp, time = year, var = "bioc", res = res, path = path)
            names(x) <- paste0("bio_", 1:nlyr(x))
            x
        }
    })

    # mask & crop
    climate_data <- purrr::map(climate_data, ~ .x %>%
        terra::crop(region) %>%
        terra::mask(region))

    names(climate_data) <- time

    # Create a SpatDataSet (ie holding a fourth dimension)
    if (length(climate_data) > 1) {
        climate_data <- terra::sds(climate_data)
    } else {
        climate_data <- terra::rast(climate_data)
    }

    return(climate_data)
}

#' Generate presence-absence data and extract environmental variables
#'
#' @param df \code{SpatVect} with occurrence data.
#' @param r \code{SpatRaster} with environmental covariates
#' @param seed seed for sampling
#' @return A \code{data.frame} with presence-absence data and extracted environmental variables.
#' @import terra
#' @export
#' @examples
#' library("terra")
#' library("geodata")
#' library("here")
#' options(geodata_default_path = here())
#' data(occurrences)
#' region <- gadm(country = "MEX", level = 0, path = here())
#' predictors <- get_predictors(
#'     region = region, res = 10,
#'     ssp = "126", model = "ACCESS-CM2",
#'     time = "1970-2000"
#' )
#' occurrence_df <- create_presence_absence_df(df = vect(occurrences), r = predictors)
create_presence_absence_df <- function(df, r, seed = 1975) {
    set.seed(seed)

    absence <- terra::spatSample(r, size = length(df), method = "random", na.rm = T)
    absence$pres <- 0
    presence <- terra::extract(r, df, df = TRUE)
    presence$pres <- 1
    occurrence_df <- dplyr::bind_rows(absence, presence) %>% dplyr::select(-ID)
    return(occurrence_df)
}

#' Fit a suitability model
#'
#' @param df \code{data.frame} containing presence-absence and environmental variables.
#' @return A fitted logistic regression model.
#' @importFrom stats glm
#' @export
fit_suitability_model <- function(df) {
    model_fit <- glm(pres ~ ., data = df, family = "binomial")
    return(model_fit)
}

#' Predict habitat suitability
#'
#' @param model_fit A fitted model.
#' @param r List of environmental data for different time periods.
#' @return Suitability rasters.
#' @import terra
#' @export
predict_suitability <- function(model_fit, r) {
    if (inherits(r, "SpatRasterDataset")) {
        r <- purrr::map(r, ~ predict(.x, model_fit, type = "response"))
        r <- terra::rast(r) # collapse
    } else {
        r <- predict(r, model_fit, type = "response")
    }
    names(r) <- paste("Suitability", names(r))
    return(r)
}

#' Plot suitability rasters using leaflet
#'
#' @param r Suitability rasters.
#' @return leaflet plots
#' @import terra leaflet leafsync
#' @export
plot_leaflet <- function(r) {
    p <- purrr::map2(as.list(r), names(r), ~ .plot_leaflet(.x, .y))
    return(leafsync::sync(p))
}

# helper
.plot_leaflet <- function(r, title, cols = c("#c45741", "#FFFFCC", "#0C2C84")) {
    pal <- leaflet::colorNumeric(cols,
        domain = c(-1, 1),
        na.color = "transparent"
    )
    leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(r, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(
            pal = pal, values = terra::values(r),
            title = title
        )
}
