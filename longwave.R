## simplified longwave model

## NEED TO TEST THIS with the Matlab code

longwave_radiation <- function(inputs, params){
  
  ## equations 14.134
  lw_down_trans <- function(x) 
    inputs$lw_sky * (1 - params$em_leaf * (1-exp(-params$Kd * x)))
  lw_down_emit <- function(x)
    params$em_leaf * params$sigma * inputs$t_leaf^4 *(1-exp(-params$Kd * x))
  
  lw_down <- function(LAI) lw_down_emit(LAI) + lw_down_trans(LAI)
  
  ## --
  lw_up_soil <- params$em_soil * params$sigma * inputs$t_soil^4
  
  ## equations 14.135
  
  lw_up_trans <- function(x) {
    l_upw_soil * (1 - params$em_leaf * (1-exp(-params$Kd * (inputs$LAI - x))))
  }
  
  lw_leaf_emit <- params$em_leaf * params$sigma * inputs$t_leaf^4
  lw_up_emit <- function(x)
    l_leaf_emit * (1-exp(-params$Kd * (inputs$LAI - x)))
  
  ## equation 14.137
  perc_abs <-  1-exp(-params$Kd * inputs$LAI) # amount absorbed
  lc <- params$em_leaf * (inputs$lw_sky + lw_up_soil) * perc_abs
         - 2 * lw_leaf_emit * perc_abs
  
  ## equation 14.138
  lc_soil <- lw_down(inputs$LAI) - lw_up_soil
  
  return(list(lc = lc, lc_soil = lc_soil))
}


params <- list(
  rho_leaf = 0.10,                    # Leaf reflectance
  tau_leaf = 0.05,                    # Leaf transmittance
  Kb = 0.58,                          # Direct beam extinction coefficient
  Kd = 0.70,                          # Diffuse extinction coefficient
  beta = 0.54,                        # Upscatter parameter for diffuse radiation
  beta0 = 0.46,                       # Upscatter parameter for direct beam radiation
  clump_OMEGA = 0.75,                 # Clumping coefficient 
  alb_soil_b = 0.1,                   # Soil albedo (direct)  
  alb_soil_d = 0.1,                   # Soil albedo (diffuse)
  em_leaf = .98,                      # Emissivity leaves
  em_soil = 1.00,                     # Emissivity soil
  sigma = 5.67e-08                    # Stefan-Boltzan costant (TODO: Unit)
  
)

params$omega_leaf = params$rho_leaf + params$tau_leaf   # Leaf scattering coefficient


inputs <- list(
  sw_sky_b = 0.8,             # Short wave (sw) direct beam radiation (W m^-2)
  sw_sky_d = 0.2,             # Short wave (sw) diffuse radiation (W m^-2)
  lw_sky = 400,               # Long wave (sw) diffuse radiation (W m^-2)
  t_leaf = 273.15 + 25,
  t_soil =273.15 + 20,
  LAI = 6                     # Leaf Area Index
)
longwave_radiation(inputs, params)
