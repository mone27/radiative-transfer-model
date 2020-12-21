library(GeoLight)
s <- solar(as.POSIXct("2020-12-21 12:00:00"))

params <-  list(lat= 51.099,
              lon = 10.426)

zenith(s, params$lon, params$lat)

get_zenith <- function(inputs, params){
  s <- solar(inputs$time)
  zenith(s, params$lon, params$lat)
}
