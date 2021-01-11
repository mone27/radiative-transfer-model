#' Entry point for the radiative transfer model

library(tidyverse)
source("radiative_transfer_model/shortwave.R")
source("radiative_transfer_model/longwave.R")
source("radiative_transfer_model/calc_parameters.R")


#' Radiative transfer model step
#'
#' This the core routine of the radiative trasfer model. It calls all the models function
#'
#' @param input A data frame row (or list) containing at least the following elements
#' - datetime
#' - sw_sky_b direct beam shortwave radiation incoming
#' - sw_sky_d diffuse shortwave radiation incoming
#' - lw_sky longwave radiation incoming
#'
#' @param p a list of the model parameters containing at least the following elements
#' max LAI value in the summer
#' min_LAI min value of LAI during winter, it is an aproximation that consider the total Plant Area Index as LAI
#' leaf_out day leaves start in spring
#' leaf_full day leaves reach max LAI
#' leaf_fall day leaves start to fall
#' leaf_fall_complete day all leaves are fallen
#'
#' lat latidude
#' lon longitude
#'
#' rho_leaf Reflencance of leaf
#' tau_leaf trasmissivity of leaf
#' omega_leaf scattering coefficient of leaf
#' clump_OMEGA canopy clumping coefficient
#' alb_soil_b soil albedo direct beam
#' alb_soil_d soil albedo diffuse
#'
#' em_leaf emittivity of leaves
#' em_soil emittivity of soil
#'
#' @return One row data Dataframe with
#' TODO document here

radiative_transfer_model_step <- function(input, p){
    LAI <- get_day_LAI(input$datetime, p$max_LAI, p$min_LAI, p$leaf_out, p$leaf_full, p$leaf_fall, p$leaf_fall_complete)

    zenith <- get_zenith(input$datetime, p$lat, p$lon) # should be 15 mins earlier because is a better average value of the half an hour interval
    Kb <- get_Kb(zenith)
    Kd <- get_Kd(LAI)
    beta <- get_beta(p$rho_leaf, p$tau_leaf)
    beta0 <- get_beta0(zenith, Kb, Kd, p$omega_leaf)

    shortwave <- shortwave_radiation(input$sw_sky_b, input$sw_sky_d, LAI, Kd, Kd, beta, beta0 , p$omega_leaf,
                                     p$clump_OMEGA, p$alb_soil_b, p$alb_soil_d)
    longwave <- longwave_radiation(input$lw_sky, LAI, input$t_leaf, input$t_soil, Kb, Kd, p$em_leaf, p$em_soil)

    # values calculated during model run outputed to give more info about the model
    interm_params <- list(LAI=LAI, Kb=Kb, Kd=Kd, beta=beta, beta0 = beta0)

    return(data.frame(c(shortwave, longwave, interm_params)))

}




#' Runs radiative transfer model over input data
#'
#' @seealso radiative_transfer_model_step
radiative_transfer_over_input <- function (input, params) {

  # Initiate an empty dataframe to hold output
  len <- nrow(input)
  out <- tibble(
        LAI = double(len),
        ic = double(len),       # total absorbed shortwave
        ic_sun = double(len),   # absorbed shortwave from sunlit leaves
        ic_sha = double(len),   # absorbed shortwave from shaded leaves
        ig = double(len),       # shortwave absorbed by the soil
        i_up = double(len),     # shortwave emitted above the canopy
        i_down = double(len),   # shortwave reaching the soil
        lc = double(len),       # total absorbed longwave from canopy
        lg = double(len),       # total absorbed longwave from soil
        lc_sun = double(len),   # absorbed longwave from sunlit leaves
        lc_sha = double(len),   # absorbed longwave from shaded leaves
        l_up = double(len),     # emitted longwave above the canopy
        l_down = double(len),    # longwave reaching the soil
        Kb = double(len),
        kd = double(len),
        beta = double(len),
        beta0 = double(len)
  )

  print("Running radiative transfer model")
  pb <- txtProgressBar(1, len)
  for (i in seq_len(nrow(input))){
    out_step <- radiative_transfer_model_step(input[i,], params)
    out[i,] <-  out_step

    setTxtProgressBar(pb, i)
  }

  return(out)
}

