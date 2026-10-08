# Description
#------------
# This script visualizes the most important SOLPS 
# simulation quantities based on state files provided by the 
# simulation. No attempt at recomputation of derived quantities is made
# as this will depend on the version of SOLPS and the available 
# physics therein. 

from python_source import solps_types as st 
from python_source import Datahandler as dh
from python_source import Plotter as pl 
from matplotlib import pyplot as plt
import numpy as np 
import sys 
import copy
import os 

# Initialize
#-----------
# Set initial paths
statefile = 'b2fstate' # state
gridfile = 'traduit.out.b2us' # grid
plasmafile = 'b2fplasmf'
gmtryfile = 'b2fgmtry'

# Set options
plotGuardCells = False 
plotstate = True  
plotresiduals = False 
plottransport = True 
plotmetrics = True 
nsmax = 2 # maximal number of species to plot 

# Read command line arguments
narg = len(sys.argv)
for i in range(1, narg):
    # First argument is state file
    if (i == 1):
        statefile = sys.argv[i]
    elif (i == 2):
        gridfile = sys.argv[i]
    elif (i == 3):
        plasmafile = sys.argv[i]
    elif (i == 4):
        gmtryfile = sys.argv[i]

simdir = os.getcwd()
statedir = simdir + '/' + statefile
plasmadir = simdir + '/' + plasmafile
griddir = simdir + '/' + gridfile
gmtrydir = simdir + '/' + gmtryfile

# Print out paths 
print('SolpsPostProcess: reading state file from: ' + statedir)
print('SolpsPostProcess: reading grid file from: ' + griddir)
print('SolpsPostProcess: reading b2fplasmf file from: ' + plasmadir)
print('SolpsPostProcess: reading b2fgmtry file from: ' + gmtrydir)


# Initialize objects
state = st.PlasmaState()
const = st.PhysicalConstants()

# Load state 
state.ReadB2fstatefile(statedir)

# Try loading b2fplasmf with residuals, transport coefficients, ...
try:
    state.ReadB2fplasmfFile(plasmadir)
except:
    plotresiduals = False
    plottransport = False 
    print("SolpsPostProcess: could not read b2fplasmf, not plotting residuals")

# Try loading b2fgmtry to be able to plot guard cell quantities
try:
    grid = dh.ReadGridFromB2fgmtryus(gmtrydir)
except:
    plotmetrics = False 
    try: 
        # Try loading through traduit file
        grid = dh.ReadTraduitOutB2us(griddir)
        print("SolpsPostProcess: could not read b2fgmtry, reading from traduit file")
    except:
        raise ValueError(("SolpsPostProcess: could not read b2fgmtry nor traduit file"))     

# Figure counter
fignum = 0

# Visualize grid
#---------------
# Face labels
pl.PlotGridFaceLabels(grid, fignum)
pl.PlotGridFaces(grid, fignum)
# Set title and other descriptors
thisaxes = plt.gca()
thisaxes.set_title('Grid face labels')
thisaxes.set_xlabel('x [m]')
thisaxes.set_ylabel('y [m]')
thisaxes.legend(loc='upper right')
fignum = fignum + 1

# Visualize plasma state
#-----------------------
# Determine plotting range
plotns = state.ns 
if plotns > nsmax:
    plotns = nsmax 
if plotstate:
    # Density
    for i in range(0, plotns):
        pl.PlotCellBasedQuantity2D(grid, np.log10(state.na[0:grid.cell.ntot, i]), fignum, doguardcells=plotGuardCells)
        thisaxes = plt.gca()
        thisaxes.set_title('log10(density) [# m^-3], species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 

    # Velocity
    for i in range(0, plotns):
        pl.PlotCellBasedQuantity2D(grid, state.ua[0:grid.cell.ntot, i], fignum, doguardcells=plotGuardCells)
        thisaxes = plt.gca()
        thisaxes.set_title('Parallel velocity [m s^-1], species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 

    # Temperature 
    pl.PlotCellBasedQuantity2D(grid, state.ti[0:grid.cell.ntot]*const.eVperJoule, fignum, doguardcells=plotGuardCells)
    thisaxes = plt.gca()
    thisaxes.set_title('Ion temperature [eV], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1
    pl.PlotCellBasedQuantity2D(grid, state.te[0:grid.cell.ntot]*const.eVperJoule, fignum, doguardcells=plotGuardCells)
    thisaxes = plt.gca()
    thisaxes.set_title('Electron temperature [eV], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1
    pl.PlotCellBasedQuantity2D(grid, state.tn[0:grid.cell.ntot]*const.eVperJoule, fignum, doguardcells=plotGuardCells)
    thisaxes = plt.gca()
    thisaxes.set_title('Neutral temperature [eV], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1

    # Potential 
    pl.PlotCellBasedQuantity2D(grid, state.po[0:grid.cell.ntot], fignum, doguardcells=plotGuardCells)
    thisaxes = plt.gca()
    thisaxes.set_title('Potential [V]')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1

# Residuals
#----------
if plotresiduals:
    # Density
    for i in range(0, plotns):
        pl.PlotCellBasedQuantity2D(grid, (state.resco[0:grid.cell.ntot, i]), fignum, doguardcells=plotGuardCells)
        thisaxes = plt.gca()
        thisaxes.set_title('Density residual [s^-1], species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 
    
    # Momentum
    for i in range(0, plotns):
        pl.PlotCellBasedQuantity2D(grid, (state.resmo[0:grid.cell.ntot, i]), fignum, doguardcells=plotGuardCells)
        thisaxes = plt.gca()
        thisaxes.set_title('Momentum residual [s^-1], species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 

    # Temperatures
    pl.PlotCellBasedQuantity2D(grid, state.reshi[0:grid.cell.ntot]*const.eVperJoule, fignum, doguardcells=plotGuardCells)
    thisaxes = plt.gca()
    thisaxes.set_title('Ion temperature residual [s^-1], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1
    pl.PlotCellBasedQuantity2D(grid, state.reshe[0:grid.cell.ntot]*const.eVperJoule, fignum, doguardcells=plotGuardCells)
    thisaxes = plt.gca()
    thisaxes.set_title('Electron temperature residual [s^-1], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1
    pl.PlotCellBasedQuantity2D(grid, state.reshn[0:grid.cell.ntot]*const.eVperJoule, fignum, doguardcells=plotGuardCells)
    thisaxes = plt.gca()
    thisaxes.set_title('Neutral temperature residual [s^-1], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1


# Transport
#----------
if plottransport:
    # cvsa
    for i in range(0, plotns):
        pl.PlotFaceBasedQuantity2D(grid, state.cvsa[0:grid.face.ntot, 0, i], fignum)
        pl.PlotGridFaces(grid, fignum)
        thisaxes = plt.gca()
        thisaxes.set_title('cvsa_pol on faces, species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 

        pl.PlotFaceBasedQuantity2D(grid, state.cvsa[0:grid.face.ntot, 1, i], fignum)
        pl.PlotGridFaces(grid, fignum)
        thisaxes = plt.gca()
        thisaxes.set_title('cvsa_rad on faces, species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 

# Grid metrics
#-------------
if plotmetrics:
    # Cell centers
    pl.PlotGridFaces(grid, fignum)
    pl.PlotPoints2D(grid.cell.x, grid.cell.y, fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('Cell centers')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1 

    # fcQalf
    pl.PlotFaceBasedQuantity2D(grid, grid.face.qalf[:, 0], fignum)
    pl.PlotPoints2D(grid.face.x[grid.face.qalf[:, 0] == 0], grid.face.y[grid.face.qalf[:, 0] == 0], fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('fcQalf[:, 0] on faces')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1 

    pl.PlotFaceBasedQuantity2D(grid, grid.face.qalf[:, 1], fignum)
    pl.PlotPoints2D(grid.face.x[grid.face.qalf[:, 1] == 0], grid.face.y[grid.face.qalf[:, 1] == 0], fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('fcQalf[:, 1] on faces')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1 

    # fcQbet
    pl.PlotFaceBasedQuantity2D(grid, grid.face.qbet[:, 0], fignum)
    pl.PlotPoints2D(grid.face.x[grid.face.qbet[:, 0] == 0], grid.face.y[grid.face.qbet[:, 0] == 0], fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('fcQbet[:, 0] on faces')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1 

    pl.PlotFaceBasedQuantity2D(grid, grid.face.qbet[:, 1], fignum)
    pl.PlotPoints2D(grid.face.x[grid.face.qbet[:, 1] == 0], grid.face.y[grid.face.qbet[:, 1] == 0], fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('fcQbet[:, 1] on faces')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1 

    
pl.ShowFigures()