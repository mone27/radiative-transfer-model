## simplified longwave model


#' Calculate the longwave radiation absorbed by the canopy with sunlit and shaded components
#' @param lw_sky_d longwave radiation above the canopy in W m-2
#' @param t_leaf temperature leaves (in Kelvin)
#' @param t_soil temperature of soil (in Kelvin)
#' @param LAI Leaf Area Index
#' @param params list containing teh values of Kd, em_leaf,
#'
#' @return list of lc, lg, lc_sun, lc_sha, l_up, l_down
longwave_radiation <- function(lw_sky, LAI, t_leaf, t_soil, params){
  ## --
  lw_soil_emit <- params$em_soil * params$sigma * t_soil^4
  lw_leaf_emit <- params$em_leaf * params$sigma  * t_leaf^4

  ## equations 14.134
  lw_down_trans <- function(x)
    lw_sky * (1 - params$em_leaf * (1-exp(-params$Kd * x)))
  lw_down_emit <- function(x)
    lw_leaf_emit *(1-exp(-params$Kd * x))

  lw_down <- function(LAI) lw_down_emit(LAI) + lw_down_trans(LAI)


  ## equations 14.135

  lw_up_trans <- function(x) {
    lw_soil_emit * (1 - params$em_leaf * (1-exp(-params$Kd * (LAI - x))))
  }

  lw_up_emit <- function(x)
    lw_leaf_emit * (1-exp(-params$Kd * (LAI - x)))

  lw_up <- function (x) lw_up_trans(x) + lw_up_emit(x)
  ## equation 14.137
  perc_abs <-  1-exp(-params$Kd * LAI) # amount absorbed
  lc <- perc_abs * (params$em_leaf * (lw_sky + lw_soil_emit)
         - 2 * lw_leaf_emit)

  ## equation 14.138
  lg<- lw_down(LAI) - lw_soil_emit # Lw adboserbed by the soil


  ### --- Sunlit and shaded leaves ---

  # Sunlit Eq. 14.140

  lc_sun <-
  (
      (
        (params$em_leaf * (lw_sky - params$sigma * t_leaf ^ 4 ) * params$Kd )
        / (params$Kd + params$Kb)
        * (1 - exp(-(params$Kd + params$Kb) * LAI))
      )
      +
      (
        (params$em_leaf *(lw_soil_emit - params$sigma * t_leaf ^ 4 ) * params$Kd)
        / (params$Kd - params$Kb)
        * (exp(-params$Kb * LAI)  - exp(-params$Kd * LAI))
      )
  )

  # Shaded Eq. 14.141

  lc_sha <- lc - lc_sun

  return(list(lc = lc,                        # Lw radiation absorbed by the canopy
              lg = lg,                        # Lw radiation absorbed by the soil
              lc_sun = lc_sun,                # Lw radiation absorbed by the sunlit canopy
              lc_sha = lc_sha,                # Lw radiation absorbed by the shaded canopy
              lw_up = lw_up(0),             # Lw emitted into the sky
              lw_down = lw_down(LAI)           # Lw reaching the soil
  ))
}
