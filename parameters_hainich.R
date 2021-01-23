#' Define the parameters for the Hainich site TODO specify origin of each value
params <- tibble(
  # NB those are test parameters that make the model work. Still need to validate that they make sense
  rho_leaf = 0.4,                     # Leaf reflectance
  tau_leaf = 0.1,                     # Leaf transmittance
  # omega_leaf =  rho_leaf + tau_leaf,  # Leaf scattering coefficient

  clump_OMEGA = 1,                    # Clumping coefficient
  alb_soil_b = 0.1,                   # Soil albedo (direct)
  alb_soil_d = 0.1,                   # Soil albedo (diffuse)
  em_leaf = .97,                      # Emissivity leaves
  em_soil = 1.00,                     # Emissivity soil
  # sigma = 5.67e-08,                   # Stefan-Boltzmann constant (TODO: Unit)

  lat = 51.099,                       # Latitude for Hainich
  lon = 10.426,                       # Longitude for Hainich
  leaf_out = 110,                     # start of leaf out
  leaf_full = 170,                    # day of full leaves
  leaf_fall = 280,                    # day of leaf fall
  leaf_fall_complete = 300,           # day of leaf fall end
  max_LAI = 5,                         # Max LAI in summer
  min_LAI = 1
)

params$omega_leaf = params$rho_leaf + params$tau_leaf   # Leaf scattering coefficient
