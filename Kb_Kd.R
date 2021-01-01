library(pracma) # for rad2deg and deg2rad
# we are assuming a spherical leaves distribution

get_Kb <-function(zenith){
        Kb <- 0.5/cos(zenith) # extinction coefficient
        return(list(Kb = Kb))
}


get_Kd <- function (LAI){
        G_z <-  0.5

        # Eq. ??
        td <-  0
        for (z in seq(0, pi / 4, pi / 18)){ # make 9 steps from 0 till π/2
             td <- td + exp( - G_z / cos(z) * LAI)*sin(z)*cos(z)*(pi / 18)
        }
        Kd <- -log(2 * td)/LAI
        return(list(Kd=Kd))
}

get_B <- function (params) {
        # Derived from equations 14.81 following the book approximation for sperical distribution
        Beta <- ( 0.625 * params$rho_leaf +  0.375 * params$tau_leaf ) / (params$rho_leaf + params$tau_leaf)

}

get_B0 <- function (zenith, params){
        attach(params)

        ross <- 0
        phi_1 <- 0.5 - 0.633 * ross - 0.333 * (ross)^2
        phi_2 <- 0.877(1 - 2 * phi_1 )

        G_mu <- 0.5
        mu <- cos(deg2rad(zenith))

        # Equation 14.84

        #defining commonly used terms
        mphi_1 <- mu * phi_1
        mphi_2 <- mu * phi_2

        # part 1 of the equation a(s) TO BE continuted
        a_s <- ((omega_leaf / 2) * (G_mu) / (G_mu + mphi_2) *
                (1 - (mphi_1/(G_mu + mphi_2) * log((G_mu + mphi_1 + mphi_2) / mphi_1))))

        B0 <-  (((Kb + Kb) / Kb) * a_s ) / omega_leaf

}