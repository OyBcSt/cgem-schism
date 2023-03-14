# cgem-schism

This repo is CGEM code that is being reformatted into a SCHISM module.  (Sorry, no actual SCHISM here!)

It can be used to run and visualize CGEM equations in a 1D water column with no boundary conditions.

Right now:
- Temperature and salinity are fixed, set in grid.nml
- RAD calculated by solar zenith angle (lat/lon/time)
- Cloern is there but I didn't test it (need to check DailyRad_init)
- no fluxes implemented yet
- no sinking yet...(this used to be handled in the Adv3D routine)

And this will never have:
- transport.  That's SCHISM's job.

This is not a production repo!  The purpose is for testing code changes.  The general workflow is to `git commit` locally and only do a `push` to GitHub when everything is working, then after the push, make necessary changes to the Python/Notebook.  That narrows the **OMG it's broken!** window a bit, but not completely.  If we're broken, please come back later. 

## Compile and run CGEM and make basic plots
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lisalenorelowe/cgem-schism.git/HEAD?labpath=cgem.ipynb)


## Run CGEM with different temperatures
Shows how to change an input with f90nml command, then rerun the code to check results.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lisalenorelowe/cgem-schism.git/HEAD?labpath=cgem_testing.ipynb)

## Check all variables after code changes
This notebook prints each variable at each k layer

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lisalenorelowe/cgem-schism.git/HEAD?labpath=cgem_check.ipynb)

## Check Phytoplankton dynamics with different parameters
This notebook shows how to change parameters and makes plots of each phytoplankton group at layer k=1.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lisalenorelowe/cgem-schism.git/HEAD?labpath=cgem_A6.ipynb)

## That's it for now!
If you happen to be browsing and see something odd, please open an Issue and let us know.  

If you browse, and make changes, figure out some good parameters, or develop new Python code to look at stuff, SAVE your work locally! The Notebooks, nml files, cgem.py, everything new thing you created will disappear forever when Binder stops.

Thanks!
