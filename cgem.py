import subprocess
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def cgem_plot1D(which_var):
    print("Plotting CGEM variable",which_var)
    results = subprocess.run(['./CGEM.exe',which_var],stdout=subprocess.PIPE, text=True)
    ar = results.stdout.splitlines()
    result = np.array(list(map(str.strip,ar))).astype(float)
    plt.plot(result,linestyle='dotted')
    plt.show()
    
def cgem_getvar(which_var):
    print("Calculating CGEM variable",which_var)
    results = subprocess.run(['./CGEM.exe',which_var],stdout=subprocess.PIPE, text=True)
    ar = results.stdout.splitlines()
    result = np.array(list(map(str.strip,ar))).astype(float)
    return result
        
def cgem_tstart(grid):
    iYrS = grid.get('time').get('iyrs')
    iMonS = grid.get('time').get('imons')
    iDayS = grid.get('time').get('idays')
    iHrS = grid.get('time').get('ihrs')
    iMinS = grid.get('time').get('imins')
    iSecS = grid.get('time').get('isecs')
    T = datetime(iYrS, iMonS, iDayS, iHrS, iMinS, iSecS)
    return T

def cgem_tend(grid):
    iYrE = grid.get('time').get('iyre')
    iMonE = grid.get('time').get('imone')
    iDayE = grid.get('time').get('idaye')
    iHrE = grid.get('time').get('ihre')
    iMinE = grid.get('time').get('imine')
    iSecE = grid.get('time').get('isece')
    T = datetime(iYrE, iMonE, iDayE, iHrE, iMinE, iSecE)
    return T

def cgem_timearray(var,grid):
    Tstart = cgem_tstart(grid)
    dtout = grid.get('time').get('dT')
    dt = timedelta(seconds=dtout)
    res = []
    for x in range (0, len(var)):
        res.append(Tstart+x*dt)
    return res

