#' Define the parameters for the Hainich site TODO specify origin of each value
params <- list(
  rho_leaf = 0.057,                   # Leaf reflectance
  tau_leaf = 0.048,                   # Leaf transmittance
  # Now calculated by the model
  #Kb = 0.58,                          # Direct beam extinction coefficient
  #Kd = 0.70,                          # Diffuse extinction coefficient
  #beta = 0.54,                        # Upscatter parameter for diffuse radiation
  #beta0 = 0.46,                       # Upscatter parameter for direct beam radiation
  clump_OMEGA = 0.75,                 # Clumping coefficient
  alb_soil_b = 0.1,                   # Soil albedo (direct)
  alb_soil_d = 0.1,                   # Soil albedo (diffuse)
  em_leaf = .97,                      # Emissivity leaves
  em_soil = 1.00,                     # Emissivity soil
  sigma = 5.67e-08,                   # Stefan-Boltzmann constant (TODO: Unit)

  lat = 51.099,                       # Latitude for Hainich
  lon = 10.426,                       # Longitude for Hainich
  leaf_out = 110,                     # start of leaf out
  leaf_full = 170,                    # day of full leaves
  leaf_fall = 280,                    # day of leaf fall
  leaf_fall_complete = 300,           # day of leaf fall end
  max_LAI = 5                         # Max LAI in summer
)

params$omega_leaf = params$rho_leaf + params$tau_leaf   # Leaf scattering coefficient
