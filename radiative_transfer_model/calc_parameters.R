library(pracma) # for rad2deg and deg2rad
library(GeoLight)
library(lubridate)


#' Simple model to get current LAI
#'
#' Use a simple model to get the LAI of the different day of the year.
#' It has 4 phases:
#'  - Winter: from leaf_fall_complete to leaf_out -> LAI = 0
#'  - Spring: from leaf_out to leaf_full -> linear growth from 0 to max_LAI
#'  - Summer: from leaf_full to leaf_fall -> LAI = max_LAI
#'  - Fall: from leaf_fall to leaf_fall_complete -> linear decrease from max_LAI to 0
#' The 4 paramenters (leaf_out...) are the day of the year
#'
#' @param time a datetime object
#' @param params list with the paramenters
#' @return LAI Leaf Area Index value for the day of the year
get_day_LAI <- function (datetime, max_LAI, min_LAI=0, leaf_out, leaf_full, leaf_fall, leaf_fall_complete){
  yday <- yday(datetime)
  if (yday < leaf_out) { # before leaves are out LAI is min
    return(min_LAI)
  }
  if (yday >= leaf_out & yday < leaf_full ) {
    ndays <-  leaf_full - leaf_out # n days of the transition
    return((max_LAI - min_LAI) * (yday - leaf_out) / ndays + min_LAI)
  }
  if (yday >= leaf_full & yday < leaf_fall ) {
    return(max_LAI)
  }
  if (yday >= leaf_fall & yday < leaf_fall_complete ) {
    ndays <-  leaf_fall_complete - leaf_fall # n days of the transition
    return((max_LAI - min_LAI) * (leaf_fall_complete - yday) / ndays + min_LAI)
  }
  if (yday >= leaf_fall_complete){
    return(min_LAI)
  }
}

#' Solar zenith from datetime and geographical coordinates
get_zenith <- function(time, params){
  s <- solar(time)
  Z <- zenith(s, params$lon, params$lat)
  min(90, Z) # the zenith function returns negative zenith during the night
}

#'  All the following function assumes a SPHERICAL leaves distribution
#' Chapter 2.2


#' Direct beam extiction coefficient
#' @param zenith in degrees
#' @return Kb
get_Kb <-function(zenith){
        # Eq. 14.29
        Kb <- 0.5/cos(deg2rad(zenith)) # extinction coefficient
        return(Kb)
}

#' Diffuse beam extiction coefficient
#' @param LAI
#' @return Kd
get_Kd <- function (LAI){
        G_z <-  0.5

        # Eq. 14.33
        td <-  0
        for (z in seq(0, pi / 4, pi / 18)){ # make 9 steps from 0 till π/2
             td <- td + exp( - G_z / cos(z) * LAI)*sin(z)*cos(z)*(pi / 18)
        }

        # Eq 14.34
        Kd <- -log(2 * td)/LAI
        return(Kd)
}

#' Fraction of diffuse light scattered backward
#' @param params list of rho_lead and and tau_leaf
#' @return beta
get_Beta <- function (params) {
        # Derived from equations 14.81 following the book approximation for sperical distribution
        Beta <- ( 0.625 * params$rho_leaf +  0.375 * params$tau_leaf ) / (params$rho_leaf + params$tau_leaf)
        return(Beta)

}

#' Fraction of direct light scattered backward
#' @param zenith in degrees
#' @param params list of Kb, Kd and omega_leaf
#' @return beta0
get_Beta_0 <- function (zenith, params){


        # Eq. 14.31
        ross <- 0
        phi_1 <- 0.5 - 0.633 * ross - 0.333 * (ross)^2
        phi_2 <- 0.877 * (1 - 2 * phi_1 )

        G_mu <- 0.5 #mu is cos(Z) but G(Z) for sperical leaves distribution is costant
        mu <- cos(deg2rad(zenith))

        # Equation 14.84

        #defining commonly used terms
        mphi_1 <- mu * phi_1
        mphi_2 <- mu * phi_2

        a_s <- ((params$omega_leaf / 2) * (G_mu) / (G_mu + mphi_2) *
                (1 - (mphi_1/(G_mu + mphi_2) * log((G_mu + mphi_1 + mphi_2) / mphi_1))))

        Beta_0 <-  (((params$Kb + params$Kd) / params$Kb) * a_s ) / params$omega_leaf
        return(Beta_0)
}

#' Updated the given paramenters list by updating with the new values depending on provided time and LAI.
update_parameters <- function (datetime, LAI, params){
  zenith <- get_zenith(datetime, params)
  params$Kb <- get_Kb(zenith)
  params$Kd <- get_Kd(LAI)
  params$beta <- get_Beta(params)
  params$beta0 <- get_Beta_0(zenith, params)

  return(params)
}