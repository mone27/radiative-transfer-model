diffuse_transmittance <-function(inputs, params){
  
        z <- inputs$zenith_angle
        Kb <- 0.5/cos(z) # extinction coefficient

        zi <- z*0.01745
        L <-inputs$LAI


        b_ext <- G(zi)/cos(zi)
     ### We don't need this code for now as we are assuming xl=0 or spherical leaves distribution
        # in R you cannot create a function in this way G(zi)
        G(zi) <- fi1+fi2*cos(zi)                  # projection of leaf area

        fi1 <- 0.5-0.633*xl-0.33*xl^2
        xl <-params$xl                            # ross index
        fi2 <- 0.877(1-2*fi1)
        E <- exp(-b_ext*L)

        ##
     ####  end

       # you need to make a loop an sum all you can't sum the numbers from 0 to 9
        Td <- 2*[sum(1:9)*E*sin(zi)*cos(zi)*zi] # diffuse radiation transmission

        ## square brackets are not valid syntax you need to always use normal brackets
        Kd <- [-log(Td)/L]                      #effective extinction coefficient


}
