**Start of the journal**

________________
*2021/01/05*

so far the longwave and shortwave model are working.
Work in progress for the calculation of the parameters.

________________
*2021/01/06*

If LAI is equal to `0` Kd is `Inf`, so need to find a way to properly deal with it.  
Same thing for Kb is the zenith is `90` you get `8.16562e+15` so not too okay.

`Beta` and `Beta_0` values seems realistic, but they'll need more verification.

Finished to code `"parameters.R"` 

----------------
*2020/01/07*

Progress with model run. Still need to debug some issues.
Added docstring and minor file structure refactor

-----------
*2020/01/08*

Fixed error in model running and added output vars to shortwave (ig, i_up, i_down)
First evaluation of model results.

Still need to double check that the equation are right (eg. no Kd instead of Kb)

----------------
*2020/01/09*
New model run gives way worse results that previous model run, with in theory the same parameters.
Need to investigate

--------------
*2020/01/11*
Fixes minor errors in shortwave.
Found a negative value in sw_sky_b from the input data. Investigate

-------------
*2020/01/11*
Refactor to have explicit arguments (instead of params) in the submodules functions (found a typo in this way :D)
LAI cannot go to 0 in winter, so fixing the issues of Inf Kb and NA outputs

Found issue in model as `out$ic, out$ic_sun + out$ic_sha)` is not true contrary to the theory.
Need to check this.
Also would be a good idea to make some test to check the energy balance of the model
Fixed issue above the order of output column was wrong, now ic is correct.

