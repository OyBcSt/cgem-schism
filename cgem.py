import subprocess
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def cgem_plot1D(grid,which_var):
    print("Plotting CGEM variable",which_var)
    results = subprocess.run(['./CGEM.exe',which_var],stdout=subprocess.PIPE, text=True)
    ar = results.stdout.splitlines()
    result = np.array(list(map(str.strip,ar))).astype(float)
    time = cgem_timearray(result,grid)
    fig, ax = plt.subplots(figsize=(15, 3))
    ax.plot(time,result,label=which_var)
    ax.legend(loc='upper left')
    
def cgem_getvar(which_var):
    print("Calculating CGEM variable",which_var)
    results = subprocess.run(['./CGEM.exe',which_var],stdout=subprocess.PIPE, text=True)
    ar = results.stdout.splitlines()
    result = np.array(list(map(str.strip,ar))).astype(float)
    return result
        
def cgem_tstart(grid):
    iYrS = grid.get('hydro').get('iyrs')
    iMonS = grid.get('hydro').get('imons')
    iDayS = grid.get('hydro').get('idays')
    iHrS = grid.get('hydro').get('ihrs')
    iMinS = grid.get('hydro').get('imins')
    iSecS = grid.get('hydro').get('isecs')
    T = datetime(iYrS, iMonS, iDayS, iHrS, iMinS, iSecS)
    return T

def cgem_tend(grid):
    iYrE = grid.get('hydro').get('iyre')
    iMonE = grid.get('hydro').get('imone')
    iDayE = grid.get('hydro').get('idaye')
    iHrE = grid.get('hydro').get('ihre')
    iMinE = grid.get('hydro').get('imine')
    iSecE = grid.get('hydro').get('isece')
    T = datetime(iYrE, iMonE, iDayE, iHrE, iMinE, iSecE)
    return T

def cgem_timearray(var,grid):
    Tstart = cgem_tstart(grid)
    dtout = grid.get('hydro').get('dt')
    dt = timedelta(seconds=dtout)
    res = []
    for x in range (0, len(var)):
        res.append(Tstart+x*dt)
    return res

