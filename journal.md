---
*2021/01/05*

**Start of the journal**

so far the longwave and shortwave model are working.
Work in progress for the calculation of the parameters.

---
*2021/01/06*

If LAI is equal to `0` Kd is `Inf`, so need to find a way to properly deal with it.  
Same thing for Kb is the zenith is `90` you get `8.16562e+15` so not too okay.

`Beta` and `Beta_0` values seems realistic, but they'll need more verification.

Finished to code `"parameters.R"` 