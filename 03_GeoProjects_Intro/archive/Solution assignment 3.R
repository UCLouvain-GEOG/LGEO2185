setwd("/Users/kvanoost/Documents/UCL/Cours/GIS/work")
options(geodata_default_path = "/Users/kvanoost/Documents/UCL/Cours/GIS/geodata")
geodata_path()

library(ggplot2) ;library(lattice)
library(latticeExtra);library(terra);library(geodata)
library(rasterVis); library(RColorBrewer); library(dplyr)

DRC <- gadm(country = "COD", level = 0, resolution = 1)

myfunction_Future_ensemble_Precipitation <- function(res,ssp,models,year,ROI) {
  #### Description: Function to estimate changes in precipitation using CMIP model scenarios
  #### Input: ROI (SpatVector), res (resolution 2.5 5 or 10), 
  #### ssp (Shared Socio-economic Pathway code: "126", "245", "370" or "585".) model (see https://www.worldclim.org/data/cmip6/cmip6climate.html)
  #### year (future year/period One of "2021-2040", "2041-2060", or "2061-2080")
  #### model he possible CMIP6 models are One of “ACCESS-CM2”, “ACCESS-ESM1-5”, “AWI-CM-1-1-MR”, “BCC-CSM2-MR”, “CanESM5”, “CanESM5-CanOE”, “CMCC-ESM2”, “CNRM-CM6-1”, “CNRM-CM6-1-HR”, “CNRM-ESM2-1”, “EC-Earth3-Veg”, “EC-Earth3-Veg-LR”, “FIO-ESM-2-0”, “GFDL-ESM4”, “GISS-E2-1-G”, “GISS-E2-1-H”, “HadGEM3-GC31-LL”, “INM-CM4-8”, “INM-CM5-0”, “IPSL-CM6A-LR”, “MIROC-ES2L”, “MIROC6”, “MPI-ESM1-2-HR”, “MPI-ESM1-2-LR”, “MRI-ESM2-0”, “UKESM1-0-LL”. 
  #### Output: a list with two rasters in unit mm with absolute & relative change relative to current precipitation
  
  #### download data and clip to ROI
  
  current <- (worldclim_global(var="prec",res=res)) %>% crop(ROI) %>% mask(ROI) %>% sum
  future <- rast()
  for (i in models) {
    print(i)
    add(future) <- (cmip6_world(i, ssp, year, var="prec", res=res,)) %>% crop(ROI) %>% mask(ROI) %>% sum
  }
  
  future <- mean(future)
  
  #### Calculate annual change
  delta_annual <- sum(future -current)
  delta_prec <- sum(future - current)/sum(current)
  return(list(abs_chang = delta_annual, rel_change = delta_prec))  
  
}
myfunction_Future_ensemble_Precipitation_fast <- function(res,ssp,models,year,ROI) {
  #### Description: Function to estimate changes in precipitation using CMIP model scenarios
  #### Input: ROI (SpatVector), res (resolution 2.5 5 or 10), 
  #### ssp (Shared Socio-economic Pathway code: "126", "245", "370" or "585".) model (see https://www.worldclim.org/data/cmip6/cmip6climate.html)
  #### year (future year/period One of "2021-2040", "2041-2060", or "2061-2080")
  #### model he possible CMIP6 models are One of “ACCESS-CM2”, “ACCESS-ESM1-5”, “AWI-CM-1-1-MR”, “BCC-CSM2-MR”, “CanESM5”, “CanESM5-CanOE”, “CMCC-ESM2”, “CNRM-CM6-1”, “CNRM-CM6-1-HR”, “CNRM-ESM2-1”, “EC-Earth3-Veg”, “EC-Earth3-Veg-LR”, “FIO-ESM-2-0”, “GFDL-ESM4”, “GISS-E2-1-G”, “GISS-E2-1-H”, “HadGEM3-GC31-LL”, “INM-CM4-8”, “INM-CM5-0”, “IPSL-CM6A-LR”, “MIROC-ES2L”, “MIROC6”, “MPI-ESM1-2-HR”, “MPI-ESM1-2-LR”, “MRI-ESM2-0”, “UKESM1-0-LL”. 
  #### Output: a list with two rasters in unit mm with absolute & relative change relative to current precipitation
  
  #### download data and clip to ROI
  
  current <- (worldclim_global(var="prec",res=res)) %>% crop(ROI) %>% mask(ROI) %>% sum
  future <- rast(nlyrs=length(models))
  for (i in models) {
    print(i)
    add(future) <- (cmip6_world(i, ssp, year, var="prec", res=res,)) %>% crop(ROI) %>% mask(ROI) %>% sum
  }
  
  future <- mean(future)
  
  #### Calculate annual change
  delta_annual <- sum(future -current)
  delta_prec <- sum(future - current)/sum(current)
  return(list(abs_chang = delta_annual, rel_change = delta_prec))  
  
}
DRC <- gadm(country="COD", level=0)
models <- c("ACCESS-CM2", "HadGEM3-GC31-LL")

### 
start_time <- Sys.time()
precip_change<-(myfunction_Future_ensemble_Precipitation(res=5,ssp="585",models=models,year="2041-2060",DRC))
end_time <- Sys.time()
end_time - start_time

####
start_time <- Sys.time()
precip_change<-(myfunction_Future_ensemble_Precipitation_fast(res=5,ssp="585",models=models,year="2041-2060",DRC))
end_time <- Sys.time()
end_time - start_time

cols <- colorRampPalette(brewer.pal(9,"Blues"))
m_abs <- mapview(precip_change$abs_chang, col.regions=cols, na.color = "transparent")
m_rel <- mapview(precip_change$rel_chang, col.regions=cols, na.color = "transparent")
leafsync::sync(m_abs, m_rel)

start_time <- Sys.time()

end_time <- Sys.time()
end_time - start_time