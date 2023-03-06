import f90nml
from opentea.noob.asciigraph import nob_asciigraph
nml = f90nml.read('cgem.nml')
nml_dict = nml.todict()
print(nob_asciigraph(nml_dict)) # show output
