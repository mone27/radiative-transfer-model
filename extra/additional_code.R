### This code has been written as part of the model development but then it ended up being not necessary


# Sperical leaf angle distribution
leaf_angle_pdf <- function(theta) {
        # theta in degrees
        # Eq. from table 2.1
        sin(deg2rad(theta))
}

cumulative_leaf_angle <- function (theta){
        # cumulative sphetical leaf angle distribution - theta in degrees
        cos(deg2rad(theta))
}

# generic leaf angle distribution
G_z_theta <- function (Z, theta){
        # Eq. 14.25
        theta <- deg2rad(theta)
        Z <- deg2rad(Z)
        if (theta < (deg2rad(90) - Z)) {
                cos(Z) * cos(theta)
        }
        else {
                a <- cos(Z) * cos(theta)
                b <- sin(Z) * sin(theta)
                c <- sqrt(sin(theta) ^ 2 - cos(Z) ^ 2)
                2 / pi * (c + a * asin(a/b) )
        }
}

# generic leaf angle distribution
G <- function(Z) {
        ## Eq. 14.26 and Eq. 2.15
        G <- 0
        for (theta in seq(5, 85, 10)) { # 10 degrees increments from 5 to 85
           G <- G + G_z_theta(Z, theta) * ( cumulative_leaf_angle(theta - 5) - cumulative_leaf_angle(theta + 5))
        }
        return(G)
}



library(testthat)
library(tidyverse)


test_that("spherical leaves distribution discrete cumulative function", {
  ## Using the numbers from page 31 tests the correct implementation of Eq. 2.15
  exp_vals <- c(0.9848078, 0.9254166, 0.8137977, 0.6634139, 0.4924039, 0.3213938, 0.1710101, 0.05939117, 1.063288e-17)
  cal_vals <- c()
  for (theta in seq(5, 85, 10)) { # 10 degrees increments from 5 to 85
    c( cal_vals, cumulative_leaf_angle(theta - 5) - cumulative_leaf_angle(theta + 5))
  }
  expect_true(all(near(exp_vals, cal_vals)))
})

