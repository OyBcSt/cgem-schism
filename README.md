# cgem-schism

This repo is CGEM code that is being reformatted into a SCHISM module.  (Sorry, no actual SCHISM here!)

It can be used to run and visualize CGEM equations in a 1D water column with no boundary conditions.

So far:
- fixed temperature and salinity
- RAD calculated by solar zenith angle (lat/lon/time)
- Cloern is there but I didn't test it (check DailyRad_init)
- no fluxes tested (yet)
- no sinking yet...(this used to be handled in the Vmixing routine)

And this will never have:
- transport.  That's SCHISM's job.

## Compile and run CGEM and make basic plots
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lisalenorelowe/cgem-schism.git/HEAD?labpath=cgem.ipynb)


## Run CGEM with different temperatures
Shows how to change an input with f90nml command, then rerun the code to check results.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lisalenorelowe/cgem-schism.git/HEAD?labpath=cgem_testing.ipynb)

## Check all variables after code changes
This notebook prints each variable at each k layer

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lisalenorelowe/cgem-schism.git/HEAD?labpath=cgem_check.ipynb)

If you happen to be browsing and see something odd, please open an Issue and let us know.  Thanks!
