library(readxl)
read_fluxnet_data <- function(path) {
  read_excel(path, na = "-9999") %>%
    select("Date/Time", SW_IN_F, SW_DIF, SW_OUT, LW_IN_F, LW_OUT, TS_F_MDS_1, TA_F, NIGHT) %>% # Columns interesting for radiative model
    rename(datetime = "Date/Time") %>%
    transmute(datetime = force_tz(datetime, "Etc/GMT+1" ),
              sw_sky_d = SW_DIF,
              sw_sky_b = SW_IN_F - sw_sky_d,  # SW_IN_F is the total radiation direct + diffuse
              lw_sky = LW_IN_F,
              t_leaf = 273.15 + TA_F,       # Temperature of air aproximation for leaf T for now
              t_soil = 273.15 + TS_F_MDS_1, # Temp soil at 2cm depth
              sw_out = SW_OUT,
              lw_out = LW_OUT,
              night = NIGHT

    )
}

#' Use this function to read the preprocessed data with the zenith added
read_fluxnet_zenith <- function(path) {
  read_csv(path, ) %>%
    select("Date/Time", SW_IN_F, SW_DIF, SW_OUT, LW_IN_F, LW_OUT, TS_F_MDS_1, TA_F, NIGHT, zenith) %>% # Columns interesting for radiative model
    rename(datetime = "Date/Time") %>%
    transmute(datetime = force_tz(datetime, "Etc/GMT+1" ),
              sw_sky = SW_IN_F,
              sw_sky_d = SW_DIF,
              sw_sky_b = SW_IN_F - sw_sky_d,  # SW_IN_F is the total radiation direct + diffuse
              lw_sky = LW_IN_F,
              t_leaf = 273.15 + TA_F,       # Temperature of air aproximation for leaf T for now
              t_soil = 273.15 + TS_F_MDS_1, # Temp soil at 2cm depth
              sw_out = SW_OUT,
              lw_out = LW_OUT,
              night = NIGHT,
              zenith = zenith

    )
}
