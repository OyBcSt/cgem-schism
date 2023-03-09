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
    
def cgem_getvar(which_var):
    print("Calculating CGEM variable",which_var)
    results = subprocess.run(['./CGEM.exe',which_var],stdout=subprocess.PIPE, text=True)
    ar = results.stdout.splitlines()
    result = np.array(list(map(str.strip,ar))).astype(float)
    return result
        
def cgem_tstart(grid):
    iYrS = grid.get('grid').get('iyrs')
    iMonS = grid.get('grid').get('imons')
    iDayS = grid.get('grid').get('idays')
    iHrS = grid.get('grid').get('ihrs')
    iMinS = grid.get('grid').get('imins')
    iSecS = grid.get('grid').get('isecs')
    T = datetime.datetime(iYrS, iMonS, iDayS, iHrS, iMinS, iSecS)
    return T

def cgem_tend(grid):
    iYrE = grid.get('grid').get('iyre')
    iMonE = grid.get('grid').get('imone')
    iDayE = grid.get('grid').get('idaye')
    iHrE = grid.get('grid').get('ihre')
    iMinE = grid.get('grid').get('imine')
    iSecE = grid.get('grid').get('isece')
    dtout = grid.get('grid').get('dt_out')
    return T

def cgem_timearray(var,grid):
    Tstart = cgem_tstart(grid)
    dtout = grid.get('grid').get('dt_out')
    dt = datetime.timedelta(seconds=dtout)
    res = []
    for x in range (0, len(A)):
        res.append(Tstart+x*dt)
    return res

