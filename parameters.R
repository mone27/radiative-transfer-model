library(pracma) # for rad2deg and deg2rad
library(GeoLight)
library(lubridate)

#' Solar zenith from datetime and geographical coordinates
get_zenith <- function(time, params){
  s <- solar(time)
  zenith(s, params$lon, params$lat)
}

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
#' @return LAI LAI value for the day of the year
get_day_LAI <- function (time, params){
  yday <- yday(time)
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

#'  All the following function assumes a SPHERICAL leaves distribution
#' Chapter 2.2

#' Direct beam extiction coefficient
get_Kb <-function(zenith){
        # Eq. 14.29
        Kb <- 0.5/cos(deg2rad(zenith)) # extinction coefficient
        return(list(Kb = Kb))
}

#' Diffuse beam extiction coefficient
get_Kd <- function (LAI){
        G_z <-  0.5

        # Eq. 14.33
        td <-  0
        for (z in seq(0, pi / 4, pi / 18)){ # make 9 steps from 0 till π/2
             td <- td + exp( - G_z / cos(z) * LAI)*sin(z)*cos(z)*(pi / 18)
        }

        # Eq 14.34
        Kd <- -log(2 * td)/LAI
        return(list(Kd=Kd))
}

#' Fraction of diffuse light scattered backward
get_Beta <- function (params) {
        # Derived from equations 14.81 following the book approximation for sperical distribution
        Beta <- ( 0.625 * params$rho_leaf +  0.375 * params$tau_leaf ) / (params$rho_leaf + params$tau_leaf)
        return(list(Beta = Beta))

}

#' Fraction of direct light scattered backward
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

        Beta_0 <-  (((Kb + Kb) / Kb) * a_s ) / params$omega_leaf
        return(list(Beta_0 = Beta_0))
}