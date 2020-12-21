library(GeoLight)
library(tidyverse)
library(lubridate)

get_zenith <- function(inputs, params){
  s <- solar(inputs$time)
  zenith(s, params$lon, params$lat)
}

day_LAI <- function (inputs, params){
  yday <- yday(inputs$time)
  if (yday < params$leaf_out) { # before leaves are out LAI is 0
    return(0)
  }
  if (yday >= params$leaf_out & yday < params$leaf_full ) {
    ndays = params$leaf_full - params$leaf_out # n days of the transition
    return(params$max_LAI * (yday - params$leaf_out) / ndays)
  }
  if (yday >= params$leaf_full & yday < params$leaf_fall ) {
    return(params$max_LAI)
  }
  if (yday >= params$leaf_fall & yday < params$leaf_fall_complete ) {
    ndays = params$leaf_fall_complete - params$leaf_fall # n days of the transition
    return(params$max_LAI * (params$leaf_fall_complete - yday) / ndays)
  }
  if (yday >= params$leaf_fall_complete){
    return(0)
  }
}







