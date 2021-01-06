## solving for direct/diffuse shortwave using the two stream approximation

direct_beam_radiation <- function(inputs, params){
  
  # defining common terms between direct/diffuse
  # --- Common terms: Eqs. (14.87) - (14.91)
  
  b <- (1 - (1 - params$beta) * params$omega_leaf) * params$Kd
  c <- params$beta * params$omega_leaf * params$Kd
  h <- sqrt(b*b - c*c)
  u <- (h - b - c) / (2 * h)
  v <- (h + b + c) / (2 * h)
  g1 <- ((params$beta0 * params$Kb - b * params$beta0 - c * (1 - params$beta0))
         * params$omega_leaf * params$Kb * inputs$sw_sky_b / (h^2 - params$Kb^2))
  g2 <- ((1 - params$beta0) * params$Kb + c * params$beta0 + b * (1 - params$beta0)) * params$omega_leaf * params$Kb * inputs$sw_sky_b / (h*h - params$Kb^2)
  
  # --- Exponential functions of leaf area
  
  s1 <- function(x) exp(-h * params$clump_OMEGA * x);
  s2 <- function(x) exp(-params$Kb * params$clump_OMEGA * x)
  
  
  
  # --- Direct beam solution
  # n1 (Eq. 14.92) and n2 (14.93)
  
  
  num1 <- v * (g1 + g2 * params$alb_soil_b + params$alb_soil_b * inputs$sw_sky_b) * s2(inputs$LAI)
  num2 <- g2 * (u + v * params$alb_soil_b) * s1(inputs$LAI)
  den1 <- v * (v + u * params$alb_soil_b) / s1(inputs$LAI)
  den2 <- u * (u + v * params$alb_soil_b) * s1(inputs$LAI)
  n2b <- (num1 - num2) / (den1 - den2)
  n1b <- (g2 - n2b * u) / v
  
  # Scattered direct beam fluxes:
  # iupwb - direct beam flux scattered upward above cumulative LAI (W/m2); Eq. (14.94)
  # idwnb - direct beam flux scattered downward below cumulative LAI (W/m2); Eq. (14.95)
  
  i_upw_b <-  function(x) -g1 * s2(x) + n1b * u * s1(x) + n2b * v / s1(x)
  i_dwn_b <-  function(x)  g2 * s2(x) - n1b * v * s1(x) - n2b * u / s1(x)
  
  # icb - direct beam flux absorbed by canopy (W/m2); Eq. (14.97)
  
  ic_b <- inputs$sw_sky_b * (1 - s2(inputs$LAI)) - i_upw_b(0) + i_upw_b(inputs$LAI) - i_dwn_b(inputs$LAI)
  
  # icsunb - direct beam flux absorbed by sunlit canopy (W/m2); Eq. (14.114)
  # icshab - direct beam flux absorbed by shaded canopy (W/m2); Eq. (14.115)
  
  a1b <- -g1 *      (1 - s2(inputs$LAI)*s2(inputs$LAI)) / (2 * params$Kb) + 
    n1b * u * (1 - s2(inputs$LAI)*s1(inputs$LAI)) / (params$Kb + h) + n2b * v * (1 - s2(inputs$LAI)/s1(inputs$LAI)) / (params$Kb - h)
  a2b <-  g2 *      (1 - s2(inputs$LAI)*s2(inputs$LAI)) / (2 * params$Kb) -
    n1b * v * (1 - s2(inputs$LAI)*s1(inputs$LAI)) / (params$Kb + h) - n2b * u * (1 - s2(inputs$LAI)/s1(inputs$LAI)) / (params$Kb - h)
  
  ic_sun_b <- (1 - params$omega_leaf) * ((1 - s2(inputs$LAI)) * inputs$sw_sky_b + params$Kd * (a1b + a2b) * params$clump_OMEGA)
  ic_sha_b <- ic_b - ic_sun_b
  
  return(list(ic_b = ic_b, ic_sun_b=ic_sun_b, ic_sha_b=ic_sha_b))
}

diffuse_radiation <- function(inputs, params){
  
  # defining common terms between direct/diffuse
  # --- Common terms: Eqs. (14.87) - (14.91)
  
  b <- (1 - (1 - params$beta) * params$omega_leaf) * params$Kd
  c <- params$beta * params$omega_leaf * params$Kd
  h <- sqrt(b*b - c*c)
  u <- (h - b - c) / (2 * h)
  v <- (h + b + c) / (2 * h)
  g1 <- ((params$beta0 * params$Kb - b * params$beta0 - c * (1 - params$beta0))
         * params$omega_leaf * params$Kb * inputs$sw_sky_b / (h^2 - params$Kb^2))
  g2 <- ((1 - params$beta0) * params$Kb + c * params$beta0 + b * (1 - params$beta0)) * params$omega_leaf * params$Kb * inputs$sw_sky_b / (h*h - params$Kb^2)
  
  # --- Exponential functions of leaf area
  
  s1 <- function(x) exp(-h * params$clump_OMEGA * x);
  s2 <- function(x) exp(-params$Kb * params$clump_OMEGA * x)
  # --- Diffuse solution
  
  # n1 (Eq. 14.99) and n2 (14.100)
  
  num <- inputs$sw_sky_d * (u + v * params$alb_soil_d) * s1(inputs$LAI)
  den1 <- v * (v + u * params$alb_soil_d) / s1(inputs$LAI)
  den2 <- u * (u + v * params$alb_soil_d) * s1(inputs$LAI)
  n2d <- num / (den1 - den2)
  n1d <- -(inputs$sw_sky_d + n2d * u) / v
  
  # Scattered diffuse fluxes:
  # iupwd - diffuse flux scattered upward above cumulative LAI (W/m2); Eq. (14.101)
  # idwnd - diffuse flux scattered downward below cumulative LAI (W/m2); Eq. (14.102)
  # and their derivatives with respect to LAI
  
  i_upw_d <- function(x)  n1d * u * s1(x) + n2d * v / s1(x)
  i_dwn_d <- function(x) -n1d * v * s1(x) - n2d * u / s1(x)
  
  
  # icd - diffuse flux absorbed by canopy (W/m2); Eq. (14.104)
  
  ic_d <- inputs$sw_sky_d - i_upw_d(0) + i_upw_d(inputs$LAI) - i_dwn_d(inputs$LAI)
  
  # icsund - diffuse flux absorbed by sunlit canopy (W/m2); Eq. (14.118)
  # icshad - diffuse flux absorbed by shaded canopy (W/m2); Eq. (14.119)
  
  a1d <-  n1d * u * (1 - s2(inputs$LAI)*s1(inputs$LAI)) / (params$Kb + h) + n2d * v * (1 - s2(inputs$LAI)/s1(inputs$LAI)) / (params$Kb - h)
  a2d <- -n1d * v * (1 - s2(inputs$LAI)*s1(inputs$LAI)) / (params$Kb + h) - n2d * u * (1 - s2(inputs$LAI)/s1(inputs$LAI)) / (params$Kb - h)
  
  ic_sun_d <- (1 - params$omega_leaf) * params$Kd * (a1d + a2d) * params$clump_OMEGA
  ic_sha_d <- ic_d - ic_sun_d
  
  return(list(ic_d = ic_d, ic_sun_d = ic_sun_d, ic_sha_d = ic_sha_d))
}


shortwave_radiation <- function(inputs, params){
 
  ib <- direct_beam_radiation(inputs, params)
  id <- diffuse_radiation(inputs, params)
  
  ic <- ib$ic_b + id$ic_d
  ic_sun <- ib$ic_sun_b + id$ic_sun_d
  ic_sha <- ib$ic_sha_b + id$ic_sha_d
  
  return(list(ic = ic, ic_sun = ic_sun, ic_sha = ic_sha))
  
}



