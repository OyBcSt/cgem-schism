import subprocess
import numpy as np
import matplotlib.pyplot as plt

def cgem_plot1D(which_var):
    print("Plotting CGEM variable",which_var)
    results = subprocess.run(['./CGEM.exe',which_var],stdout=subprocess.PIPE, text=True)
    ar = results.stdout.splitlines()
    result = np.array(list(map(str.strip,ar))).astype(float)
    plt.plot(result,linestyle='dotted')
    plt.show()
