# Main plotter class
import matplotlib as mpl
mpl.use('TkAgg') # set the gui backend for the cluster...
from matplotlib import pyplot as plt
import numpy as np
from . import Datahandler as dh
import time
from . import goat_types as gt
import copy

#==========================================================================#
#                                                                          #
#                              GLOBAL VARIABLES                            #
#                                                                          #
#==========================================================================#

# Special characters
#-------------------
filesep = '/' # file separator

# Data
#-----
GUARD_CELL_WIDTH = 0.1 # width of guard cell for plotting

# Enable gui events
#------------------

# Paths & files
#--------------
# file where grid vertices are stored in [ID, x, y] format
gridverticesfile = 'vertices_init.dat' # all grid vertices (initial coordinates)

# Grid deformation cost function
cfvLRvertexpairsfile = 'costfunction_vertexpairs_LR.dat'
cfvLRdesiredbiasfile = 'costfunctionLR_desiredbias.dat'
cfvLRbiasfile = 'costfunctionLR_desiredbias.dat'
cfvLRradvertexpairsfile = 'costfunction_vertexpairs_LRrad.dat'
cfvFADvertexpairsfile = 'costfunction_vertexpairs_FAD.dat'
cfvFADweightsfile = 'costfunctionFAD_weights.dat'
cfvFAvertexpairsfile = 'costfunction_vertexpairs_FA.dat'
cfvFAweightsfile = 'costfunctionFA_weights.dat'
cfvPRPBbiasfile = 'costfunctionPRPB_desiredbias.dat'
cfvPRPBweightfile = 'costfunctionPRPB_weight.dat'

cfvPRPBfspairsfile = 'costfunctionPRPB_fspairs.dat'
cfvPRPBvertfile = 'costfunction_vertexpairs_PRPB.dat'

cfvLRvalverticesfileit = 'cfv_lr_val_vertices.dat'
cfvFADvalverticesfileit = 'cfv_fad_val_vertices.dat'
cfvLRradvalverticesfileit = 'cfv_lrrad_val_vertices.dat'
cfvFAvalverticesfileit = 'cfv_fa_val_vertices.dat'
cfvPRPBvalverticesfileit = 'cfv_prpb_val_vertices.dat'


# Grid deformation equality constraints
bndconverticesfile = 'con_bnd_vertices.dat' # vertices constrained by boundary constraints
bndconverticesvalfile = 'con_bnd_val_vertices.dat'
ffconverticesfile = 'con_ff_vertices.dat' # vertices constrained by flux function constraints
ffconverticesvalfile = 'con_ff_val_vertices.dat'
ffconverticesfilestem = 'con_ff_vertices_'
xpconverticesfile = 'con_xp_vertices.dat' # vertices constrained as x-points
elconvertexpairsfile = 'con_el_vertices.dat' # vertex pairs constrained for edge lengths
orthconvertexpairsfile = 'con_orth_vertices.dat' # vertex pairs constrained for orthogonality
orthconverticesvalfile = 'con_orth_val_vertices.dat'
ffvconverticesfile = 'con_ffv_vertices.dat' # vertices constrained by constant flux value for flux surface 

# Grid deformation inequality constraints
convertpairsLFpol = 'con_lf_vpairspol.dat'
convertpairsLFrad = 'con_lf_vpairsrad.dat'
convertpairsLFves = 'con_lf_vpairsves.dat'
convertpairsLFpolit = 'con_lf_vpairspol_iterate.dat'
convertpairsLFradit = 'con_lf_vpairsrad_iterate.dat'
convertpairsLFvesit = 'con_lf_vpairsves_iterate.dat'


# file where grid cells are stored in polygon format ([x, y] with blank
# lines between polygons)
gridcellsfile = 'cells_init.dat'
gridcellsiteratefile = 'cells_iterate.dat'

# Files with optimization history
goathistoryfile = 'goat_optimization_history.dat'
shapeopthistoryfile = 'so_optimization_history.dat'

# Shape optimization paths
vesselpolygonfile = 'vesselpolygon.dat'
fvpfile = 'so_con_fvp_vertices.dat' # fixed vessel points file
fvffile = 'so_con_fvf_vertices.dat' # fixed vessel flux file
origvesselpolygonfile = 'vesselpolygon_orig_so.dat' # original/initial vessel file
currentvesselpolygonfile = 'vesselpolygon_iterate_so.dat' # current/new vessel file

# Grid generation paths
tm_beforecellsfile = 'topomesh_beforecells.dat'
tmfile = 'topomesh.dat'

#==========================================================================#
#                                                                          #
#                          SPECIFIC PLOTTING ROUTINES                      #
#                                                                          #
#==========================================================================#

#--------------------------------------------------------------------------#
#                                   Grid                                   #
#--------------------------------------------------------------------------#


def PlotGridCellsFromFile(dirpath, fignum):
    # Description
    #------------
    # Plot the grid cells once by reading in the grid cell polygon data as
    # written out in the cells.dat file

    # Get filepath
    filepath = dirpath + filesep + gridcellsfile
    vesselfilepath = dirpath + filesep + vesselpolygonfile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)
    vesselvals = dh.GetPolygonCoordinates(vesselfilepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='',
        linewidth=0.1, label='Grid cells')
    PlotPolygons2D(vesselvals[:, 0], vesselvals[:, 1], fignum, color='b', marker='',
        linewidth=0.1, label='Vessel')

    # Set axes
    SetAxesLimits2D(plt.gca(), vesselvals[:, 0], vesselvals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Grid cells')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotGridCellsIterate(dirpath, fignum):
    # Description
    #------------
    # Plot the grid cells once by reading in the grid cell polygon data as
    # written out in the cells.dat file

    # Get filepath
    filepath = dirpath + filesep + gridcellsiteratefile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='',
        linewidth=0.1, label='Grid cells')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Grid cells')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def MonitorGrid(datadir, num, pausetime, maxruntime):
    # Description
    #------------
    # This routine makes a plot that is continuously updated.

    # Time interval to replot [s]
    starttime = time.time()

    # Prepare for gui event loop
    plt.ion()

    # Plot the data
    PlotGridCellsIterate(datadir, 1)
    thisfig = plt.gcf()

    # Loop until time has passed
    while (time.time() - starttime <= maxruntime):
        try: 
            # Plot the data
            PlotGridCellsIterate(datadir, 1)

            # Plot constraints
            PlotActiveLinefoldingConstraintsIterate(datadir, 1)

            # Draw
            thisfig.canvas.draw()
            thisfig.canvas.flush_events()

            # Pause
            time.sleep(pausetime)

            # Clear figure
            ClearCurrentAxes()
        except: 
            time.sleep(pausetime)

def PlotGridFaceLabels(grid, fignum):
    # Description
    #------------
    # Plot facelabels of an unstructured simulation grid

    # Compute face locations
    xf = 0.5*(grid.vert.x[grid.face.v1-1] + grid.vert.x[grid.face.v2-1])
    yf = 0.5*(grid.vert.y[grid.face.v1-1] + grid.vert.y[grid.face.v2-1])

    # Only plot non-zero labels
    ind = np.nonzero(grid.face.label)

    # Plot
    PlotPoints2DWithID(xf[ind], yf[ind], grid.face.label[ind], fignum) 

def PlotGridFaceRegions(grid, fignum):
    # Description
    #------------
    # Plot facelabels of an unstructured simulation grid

    # Compute face locations
    xf = 0.5*(grid.vert.x[grid.face.v1-1] + grid.vert.x[grid.face.v2-1])
    yf = 0.5*(grid.vert.y[grid.face.v1-1] + grid.vert.y[grid.face.v2-1])

    # Only plot non-zero labels
    ind = np.nonzero(grid.face.region)

    # Plot
    PlotPoints2DWithID(xf[ind], yf[ind], grid.face.region[ind], fignum) 

def PlotGridVertFieldlineID(grid, fignum, vID=None):
    # Description
    #------------
    # plot vertex fieldline IDs (may be heavy)

    if vID == None:
        tvID = grid.vert.fieldlineID
        tx = grid.vert.x
        ty = grid.vert.y
    else:
        hasID = np.zeros(grid.vert.ntot)
        for i in vID:
            hasID = np.where(grid.vert.fieldlineID == i, 1, hasID)
        this = np.array(hasID > 0)
        tx = grid.vert.x[this]
        ty = grid.vert.y[this]
        tvID = grid.vert.fieldlineID[this]
        

    # Plot
    PlotPoints2DWithID(tx, ty, tvID, fignum) 

def PlotGridCellArea(grid, fignum, doinverse=False):
    # Description
    #------------
    # Make a plot of the cell area of the grid, just to detect small
    # cells. 

    # Calculate surface area
    surfA = grid.ComputeCellSurfaceArea()

    # Plot
    
    
    if doinverse:
        PlotCellBasedQuantity2D(grid, np.log10(1/np.abs(surfA)), fignum)
        thisaxes = plt.gca()
        thisaxes.set_title('Log10 of grid cell surface area (inverse) [m^-2]')
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
    else:
        PlotCellBasedQuantity2D(grid, np.log10(np.abs(surfA)), fignum)
        thisaxes = plt.gca()
        thisaxes.set_title('Log10 of grid cell surface area [m^2]')
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')

    # Plot minimal location
    loc = np.argmin(abs(surfA))
    PlotPoints2D(grid.cell.x[loc], grid.cell.y[loc], fignum, color='r', marker='o')
    print ("PlotGridCellArea: smallest cell located at ", grid.cell.x[loc], grid.cell.y[loc], " with surface area: ", np.min(abs(surfA)))

#--------------------------------------------------------------------------#
#                              Grid Optimization                           #
#--------------------------------------------------------------------------#
def PlotLRCostfunctionVertexPairs(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs in a quiver-like plot that are considered 
    # in the length ratio costfunction (or any derivate of that costfunction)

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + cfvLRvertexpairsfile

    # Get the vertex pairs
    vals = dh.GetVertexPairsPairCoordinates(vertfilepath)

    # Plot in a quiver style
    PlotPoints2D(0.5*(vals[:, 0] + vals[:, 2]), 0.5*(vals[:, 1] + vals[:, 3]), fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotPoints2D(0.5*(vals[:, 4] + vals[:, 6]), 0.5*(vals[:, 5] + vals[:, 7]), fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotVertexPairs2DQuiver(0.5*(vals[:, 0] + vals[:, 2]), 0.5*(vals[:, 1] + vals[:, 3]), 
        0.5*(vals[:, 4] + vals[:, 6]), 0.5*(vals[:, 5] + vals[:, 7]), 
        fignum, color='b')
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Length ratio costfunction pairs of vertex pairs')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotLRCostfunctionDesiredBias(dirpath, fignum):
    # Description
    #------------
    # Plot the desired bias for the length ratio cost function
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + cfvLRdesiredbiasfile

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Length ratio costfunction desired bias')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotLRCostfunctionValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + cfvLRvalverticesfileit

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Length ratio costfunction value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    
def PlotLRradCostfunctionVertexPairs(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs in a quiver-like plot that are considered 
    # in the length ratio costfunction (or any derivate of that costfunction)

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + cfvLRradvertexpairsfile

    # Get the vertex pairs
    vals = dh.GetVertexPairsPairCoordinates(vertfilepath)

    # Plot in a quiver style
    PlotPoints2D(0.5*(vals[:, 0] + vals[:, 2]), 0.5*(vals[:, 1] + vals[:, 3]), fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotPoints2D(0.5*(vals[:, 4] + vals[:, 6]), 0.5*(vals[:, 5] + vals[:, 7]), fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotVertexPairs2DQuiver(0.5*(vals[:, 0] + vals[:, 2]), 0.5*(vals[:, 1] + vals[:, 3]), 
        0.5*(vals[:, 4] + vals[:, 6]), 0.5*(vals[:, 5] + vals[:, 7]), 
        fignum, color='b')
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Radial length ratio costfunction pairs of vertex pairs')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotLRradCostfunctionValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + cfvLRradvalverticesfileit

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Radial length ratio costfunction value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFACostfunctionVertexPairs(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs in a quiver-like plot that are considered 
    # in the face angle costfunction (or any derivate of that costfunction)

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + cfvFAvertexpairsfile
    fffilepath = dirpath + filesep + ffconverticesfilestem

    # Get the vertex pairs
    vals = dh.GetVertexPairCoordinates(vertfilepath)

    # Plot in a quiver style
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotPoints2D(vals[:, 2], vals[:, 3], fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotVertexPairs2DQuiver(vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3], 
        fignum, color='b')
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Face angle costfunction vertex pairs')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFACostfunctionValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + cfvFAvalverticesfileit

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Face angle costfunction value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFADCostfunctionVertexPairs(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs in a quiver-like plot that are considered 
    # in the face angle difference costfunction (or any derivate of that costfunction)

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + cfvFADvertexpairsfile

    # Get the vertex pairs
    vals = dh.GetVertexPairsPairCoordinates(vertfilepath)

    # Plot in a quiver style
    PlotPoints2D(0.5*(vals[:, 0] + vals[:, 2]), 0.5*(vals[:, 1] + vals[:, 3]), fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotPoints2D(0.5*(vals[:, 4] + vals[:, 6]), 0.5*(vals[:, 5] + vals[:, 7]), fignum, color='b', marker='o',
        facecolors='none', label='Vertices')
    PlotVertexPairs2DQuiver(0.5*(vals[:, 0] + vals[:, 2]), 0.5*(vals[:, 1] + vals[:, 3]), 
        0.5*(vals[:, 4] + vals[:, 6]), 0.5*(vals[:, 5] + vals[:, 7]), 
        fignum, color='b')
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Face angle difference costfunction pairs of vertex pairs')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFADCostfunctionValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + cfvFADvalverticesfileit

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Face angle difference costfunction value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotPRPBCostfunctionVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the PRPB vertices

    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    prpbvertfilepath = dirpath + filesep + cfvPRPBvertfile

    # Get the data
    vals = dh.GetVertexCoordinates(vertfilepath)
    valscon = dh.GetVertexCoordinates(prpbvertfilepath)

    # Plot the data
    PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
        marker='+', label='Flux surface vertices')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vertices of flux surfaces considered in PRPB cost function')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotPRPBCostfunctionValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + cfvPRPBvalverticesfileit

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Psi ratio costfunction value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFluxfunctionConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vertices that have their flux values constrained. Vertices
    # should be stored in ID, x, y format.

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    fffilepath = dirpath + filesep + ffconverticesfilestem

    # Get the grid data
    vals = dh.GetVertexCoordinates(vertfilepath)

    # Special points
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'sp.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='*', label='Special vertices')
        valscon = dh.GetVertexCoordinates(fffilepath+'spfs.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='+', label='Special flux surface vertices')
    except:
        print('could not print special points of flux function constraints')

    # Tangency points
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'tp.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='s', label='Tangency points')
    except:
        print('could not print tangency points of flux function constraints')

    # Fixed points
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'fp.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='g',
                 marker='*', label='Fixed points')
    except:
        print('could not print fixed points of flux function constraints')

    # Flux surfaces
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'fs.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
                 marker='+', label='Flux surfaces')
    except:
        print('could not print flux surface points of flux function constraints')
    

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Flux function constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFluxFunctionConstraintValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + ffconverticesvalfile

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Flux function constraint value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotBoundaryConstraintValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + bndconverticesvalfile

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Boundary constraint value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotBoundaryConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # This routine plots the vertices that are constrained by boundary
    # constraints. The data of all vertex coordinates should be stored in
    # a vertices.dat file, together with the boundary vertex data in a
    # con_bnd_vertices.dat file. Both should be present in the folder
    # specified in 'dirpath'.

    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    bndfilepath = dirpath + filesep + bndconverticesfile

    # Get the data
    vals = dh.GetVertexCoordinates(vertfilepath)
    valscon = dh.GetVertexCoordinates(bndfilepath)

    # Plot the data
    PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
        marker='+', label='Constrained vertices')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vertices constrained with boundary constraints')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotBoundaryConstraintPLF(dirpath, fignum):
    # Description
    #------------
    # Plot the polygon levelset function of the boundary constraints
    # Set the filepaths
    datafilepath = dirpath + filesep + bndconPLFfile

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum, levels=25)

def PlotXPointConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vertices that are x-points and are constrained as such

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    confilepath = dirpath + filesep + xpconverticesfile

    # Get the data
    vals = dh.GetVertexCoordinates(vertfilepath)
    valscon = dh.GetVertexCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
                 marker='+', label='Constrained vertices')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('X-point constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotEdgelengthsConstraintEdges(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs that belong to the edgelengths.

    # Set filepaths
    # Set the filepaths
    cellfilepath = dirpath + filesep + gridcellsfile
    confilepath = dirpath + filesep + elconvertexpairsfile

    # Get the data
    vals = dh.GetPolygonCoordinates(cellfilepath)
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='b',
        marker='+', label='Constrained edges')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Edge length constrained edges')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFixedFluxFunctionConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vertices that belong to flux surfaces of which the 
    # flux value is being constrained. 

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    confilepath = dirpath + filesep + ffvconverticesfile

    # Get the data
    vals = dh.GetVertexCoordinates(vertfilepath)
    valscon = dh.GetVertexCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
                 marker='+', label='Constrained vertices')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Fixed flux value for flux surfaces constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotOrthogonalityConstraintEdges(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs that belong to the edgelengths.

    # Set filepaths
    # Set the filepaths
    cellfilepath = dirpath + filesep + gridcellsfile
    confilepath = dirpath + filesep + orthconvertexpairsfile

    # Get the data
    vals = dh.GetPolygonCoordinates(cellfilepath)
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='b',
        marker='o', label='Constrained edges')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Orthogonality constrained edges')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotOrthogonalityConstraintValueAtVertices(dirpath, fignum):
    # Description
    #------------
    # Visualize the cost function value per vertex
    # Set filepaths
    # Set the filepaths
    datafilepath = dirpath + filesep + orthconverticesvalfile

    # Plot 
    Plot2DSurfaceDataContourf(datafilepath, fignum)
    
    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Orthogonality constraint value per vertex')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotLinefoldingConstraintEdges(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs that belong to the linefolding constraints.
    # This includes poloidal, radial, and vessel line folding 
    # constraints. 

    # General
    #--------
    # Underlying grid
    cellfilepath = dirpath + filesep + gridcellsfile
    vals = dh.GetPolygonCoordinates(cellfilepath)

    # Poloidal
    #---------
    # Set filepaths
    confilepath = dirpath + filesep + convertpairsLFpol

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='b',
        marker='o', label='Constrained edges')
    
    # Radial
    #-------
    confilepath = dirpath + filesep + convertpairsLFrad

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='g',
        marker='x', label='Constrained edges')
    
    # Vessel
    #-------
    confilepath = dirpath + filesep + convertpairsLFves

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='k',
        marker='+', label='Constrained edges')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Linefolding constrained edges')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotActiveLinefoldingConstraintsIterate(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs that belong to the linefolding constraints.
    # This includes poloidal, radial, and vessel line folding 
    # constraints. 

    # Poloidal
    #---------
    # Set filepaths
    confilepath = dirpath + filesep + convertpairsLFpolit

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    if len(valscon) > 0: 
        PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
            0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='b',
            marker='o', label='Constrained edges')
    
    # Radial
    #-------
    confilepath = dirpath + filesep + convertpairsLFradit

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    if len(valscon) > 0: 
        PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
            0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='g',
            marker='x', label='Constrained edges')
    
    # Vessel
    #-------
    confilepath = dirpath + filesep + convertpairsLFvesit

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data -  if available
    if len(valscon) > 0: 
        PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
            0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='k',
            marker='+', label='Constrained edges')

    # Set axes
    #SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.legend(loc='upper right')

def PlotGoatOptimizationHistory(dirpath, fignum):
    # Description
    #------------
    # Plot the goat optimization history on a log(value)-iterate plot.
    # This includes the cost function, (max) value of 
    # constraints, the convergence history (infinity norm of 
    # lagrangian gradient), and the linesearch step length 

    # Read data
    #----------
    historypath = dirpath + filesep + goathistoryfile
    [valnames, vals] = dh.ReadGeneralColumnwiseFloatData(historypath)

    # Hedge for negative numbers
    for j in range(len(vals[1,:])):
        for i in range(len(vals[:, j])):
            if (vals[i, j] <= 0):
                vals[i, j] = np.nan

    # Set figure
    #-----------
    plt.figure(fignum)

    # Plot data
    #----------
    # First entry should be iteration counter, then convnorm, dphi, L, 
    # J, max(G), max(H), rxf, step, tol, ...

    # Convnorm
    plt.plot(vals[:, 0], vals[:, 1], 'rx-', label='max(abs(grad L))')
    
    # Cost function
    plt.plot(vals[:, 0], vals[:, 4], 'bx-', label='J')

    # Equality constraints
    plt.plot(vals[:, 0], vals[:, 5], 'gx-', label='max(G)')

    # Inequality constraints
    plt.plot(vals[:, 0], vals[:, 6], 'mx-', label='max(H)')

    # Line search step length
    plt.plot(vals[:, 0], vals[:, 8], 'kx-', label='alpha_ls')

    # Set figure data
    #----------------
    # Set axes
    SetAxesLimitsLogplot(plt.gca(), vals[:, 0], vals[:, [1, 4, 5, 6, 8]])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Goat grid deformation convergence history')
    thisaxes.set_xlabel('iteration number')
    thisaxes.set_ylabel('value')
    thisaxes.set_yscale('log')
    thisaxes.legend(loc='upper right')

#--------------------------------------------------------------------------#
#                             Grid Generation                              #
#--------------------------------------------------------------------------#

# Topological mesh plotting from file
def PlotTopologicalMeshFromFile(dirpath, fignum):
    filepath = dirpath + filesep + tmfile
    topomesh = dh.ReadTopomeshFile(dirpath)
    PlotTopologicalMesh(topomesh, fignum)

# Topological mesh plotting
def PlotTopologicalMesh(topomesh, fignum, plotvertID=True):
    # Description
    #------------
    # Plot the topological mesh 

    # Plot the vertices
    #------------------
    # Plot per vertex type
    regular = np.where(topomesh.vert.type == gt.TMvertexregularID) 
    bnd = np.where(topomesh.vert.type == gt.TMvertexbndID)
    split = np.where(topomesh.vert.type == gt.TMvertexsplitID) 
    minima = np.where(topomesh.vert.type == gt.TMvertexminID)
    saddle = np.where(topomesh.vert.type == gt.TMvertexsaddleID)
    maxima = np.where(topomesh.vert.type == gt.TMvertexmaxID)
    tp1 = np.where(topomesh.vert.type == gt.TMvertextp1ID)
    tp2 = np.where(topomesh.vert.type == gt.TMvertextp2ID)
    if plotvertID:
        PlotPoints2DWithID(topomesh.vert.x[maxima], topomesh.vert.y[maxima], topomesh.vert.ID[maxima], fignum, color='b', 
            marker='o', label='field maxima')
        PlotPoints2DWithID(topomesh.vert.x[saddle], topomesh.vert.y[saddle], topomesh.vert.ID[saddle], fignum, color='b', 
            marker='x', label='field saddle')
        PlotPoints2DWithID(topomesh.vert.x[minima], topomesh.vert.y[minima], topomesh.vert.ID[minima], fignum, color='b', 
            marker='s', label='field minima')
        PlotPoints2DWithID(topomesh.vert.x[tp1], topomesh.vert.y[tp1], topomesh.vert.ID[tp1], fignum, color='r', 
            marker='o', label='field tangency point type 1')
        PlotPoints2DWithID(topomesh.vert.x[tp2], topomesh.vert.y[tp2], topomesh.vert.ID[tp2], fignum, color='r', 
            marker='x', label='field tangency point type 2')
        PlotPoints2DWithID(topomesh.vert.x[regular], topomesh.vert.y[regular], topomesh.vert.ID[regular], fignum, color='k', 
            marker='.', label='regular vertex')
        PlotPoints2DWithID(topomesh.vert.x[bnd], topomesh.vert.y[bnd], topomesh.vert.ID[bnd], fignum, color='k', 
            marker='.', label='bnd vertex')
        PlotPoints2DWithID(topomesh.vert.x[split], topomesh.vert.y[split], topomesh.vert.ID[split], fignum, color='g', 
            marker='d', label='splitted face vertex')
    
    PlotPoints2D(topomesh.vert.x[maxima], topomesh.vert.y[maxima],  fignum, color='b', 
        marker='o', label='field maxima')
    PlotPoints2D(topomesh.vert.x[saddle], topomesh.vert.y[saddle],  fignum, color='b', 
        marker='x', label='field saddle')
    PlotPoints2D(topomesh.vert.x[minima], topomesh.vert.y[minima],  fignum, color='b', 
        marker='s', label='field minima')
    PlotPoints2D(topomesh.vert.x[tp1], topomesh.vert.y[tp1],  fignum, color='r', 
        marker='o', label='field tangency point type 1')
    PlotPoints2D(topomesh.vert.x[tp2], topomesh.vert.y[tp2],  fignum, color='r', 
        marker='x', label='field tangency point type 2')
    PlotPoints2D(topomesh.vert.x[regular], topomesh.vert.y[regular], fignum, color='k', 
        marker='.', label='regular vertex')
    PlotPoints2D(topomesh.vert.x[bnd], topomesh.vert.y[bnd],  fignum, color='k', 
        marker='.', label='bnd vertex')
    PlotPoints2D(topomesh.vert.x[split], topomesh.vert.y[split],  fignum, color='g', 
        marker='d', label='splitted face vertex')
    
    # Store bounds
    xb = [np.min(topomesh.vert.x), np.max(topomesh.vert.x)]
    yb = [np.min(topomesh.vert.y), np.max(topomesh.vert.y)]
    
    # Plot the faces
    #---------------
    # Plot per face type
    pf = np.where(topomesh.face.type == gt.TMfacepolID)
    rf = np.where(topomesh.face.type == gt.TMfaceradID)
    bndf = np.where(topomesh.face.type == gt.TMfacebndID)
    sepf = np.where(topomesh.face.type == gt.TMfacesepID)
    coref = np.where(topomesh.face.type == gt.TMfacecoreID)
    PFf = np.where(topomesh.face.type ==  gt.TMfacePFID)
    SOLf = np.where(topomesh.face.type ==  gt.TMfaceSOLID)
    albndf = np.where(topomesh.face.type ==  gt.TMfacealbndID)

    for i in pf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='r')
    for i in rf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='g')
    for i in bndf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='k')
    for i in sepf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='m')
    for i in coref[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='r')
    for i in PFf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='r')
    for i in SOLf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='r')
    for i in albndf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='g')
        
    # Check bounds
    for i in np.arange(0, topomesh.face.ntot, 1):
        xb[0] = np.min([xb[0], np.min(topomesh.face.data[i].x)])
        xb[1] = np.max([xb[1], np.max(topomesh.face.data[i].x)])
        yb[0] = np.min([yb[0], np.min(topomesh.face.data[i].y)])
        yb[1] = np.max([yb[1], np.max(topomesh.face.data[i].y)])

    # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

    # thisaxes.legend(loc='upper right')

# Topological cell plotting
def PlotTopologicalMeshCells(topomesh, fignum, **plotargs):
    # Description
    #------------
    # Separate routine to plot topological mesh cells in order not to
    # overburden the standard plot
    maxcellvert = 1000
    thisrange = np.arange(0, topomesh.cell.ntot, 1)
    for i in thisrange:
        index = np.arange(0, len(topomesh.cell.data[i].x), len(topomesh.cell.data[i].x)/maxcellvert, dtype=int)
        #PlotPolygons2D(topomesh.cell.data[i].x[index], topomesh.cell.data[i].y[index], fignum)
        PlotPolygons2D(topomesh.cell.data[i].x, topomesh.cell.data[i].y, fignum, **plotargs)

# Single topological cell, filled 
def PlotSingleTopologicalMeshCellFilled(topomesh, fignum, cellID, **plotargs):
    # Description
    #------------
    # Plot a single cell of which the surface will be filled according
    # to the arguments given in **plotargs

    ind = cellID-1 # zero-based indexing
    PlotFilledPolygons2D(topomesh.cell.data[ind].x, topomesh.cell.data[ind].y, fignum, **plotargs)


# Topological mesh tube plotting
def PlotTopologicalMeshTubes(topomesh, fignum):
    # Description
    #------------
    # This routine plots the topological mesh tubes by plotting their 
    # cells and plotting the tube number on each radial face (as these
    # should be unique for each tube)
    
    # Plot cells
    PlotTopologicalMeshCells(topomesh, fignum, color='k')

    # Loop over all tubes and plot the numbers at approx. the middle of 
    # the face
    for i in np.arange(0, topomesh.ft.ntot, 1):
        # Get the tube faces
        tf = topomesh.ft.GetFace(i)

        # Plot tube ID for each face
        xf = np.zeros(len(tf), dtype=float)
        yf = np.zeros(len(tf), dtype=float)
        ftID = np.zeros(len(tf), dtype=int)
        cc = 0 
        for j in tf:
            txf, tyf = topomesh.face.data[j-1].InterpolateCoordinates(0.5)
            xf[cc] = txf 
            yf[cc] = tyf
            ftID[cc] = topomesh.ft.ID[i]
            cc = cc + 1
        PlotPoints2DWithID(xf, yf, ftID, fignum, color='r')


# Grid generation data plotting: face vertex distributions
def PlotGGTMDataFaceVertexDistribution(ggtmdata, topomesh, fignum):
    # Plot all vertices that were distributed on faces only
   
    # Plot per face type
    pf = np.where(topomesh.face.type == gt.TMfacepolID)
    rf = np.where(topomesh.face.type == gt.TMfaceradID)
    bndf = np.where(topomesh.face.type == gt.TMfacebndID)
    sepf = np.where(topomesh.face.type == gt.TMfacesepID)

    for i in pf[0]:
        if not (ggtmdata.face[i].ID == 0):
            PlotPolygons2D(ggtmdata.face[i].x, ggtmdata.face[i].y, 
                fignum, color='r', marker='.')
    for i in rf[0]:
        if not (ggtmdata.face[i].ID == 0):
            PlotPolygons2D(ggtmdata.face[i].x, ggtmdata.face[i].y, 
                fignum, color='g', marker='.')
    for i in bndf[0]:
        if not (ggtmdata.face[i].ID == 0):
            PlotPolygons2D(ggtmdata.face[i].x, ggtmdata.face[i].y, 
                fignum, color='k', marker='.')
    for i in sepf[0]:
        if not (ggtmdata.face[i].ID == 0):
            PlotPolygons2D(ggtmdata.face[i].x, ggtmdata.face[i].y, 
                fignum, color='m', marker='.')
        
    # Check bounds
    xb = [np.inf, -np.inf]
    yb = [np.inf, -np.inf]
    for i in np.arange(0, topomesh.face.ntot, 1):
        xb[0] = np.min([xb[0], np.min(topomesh.face.data[i].x)])
        xb[1] = np.max([xb[1], np.max(topomesh.face.data[i].x)])
        yb[0] = np.min([yb[0], np.min(topomesh.face.data[i].y)])
        yb[1] = np.max([yb[1], np.max(topomesh.face.data[i].y)])

    # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('topomesh face vertex distribution')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

# Grid generation data plotting: cell vertex distributions
def PlotGGTMDataCellVertexDistribution(ggtmdata, topomesh, fignum):
    # Plot all vertices that were distributed on cells only

    # Initialize plotting bounds
    xb = [np.inf, -np.inf]
    yb = [np.inf, -np.inf]

    # Loop over all cells
    for thiscell in ggtmdata.cell:
        # Plot high field line and low field line
        nt = len(thiscell.tubes)-1
        PlotPolygons2DQuiver(thiscell.tubes[0].hfline.x, thiscell.tubes[0].hfline.y, fignum, color='r')
        PlotPolygons2DQuiver(thiscell.tubes[nt].lfline.x, thiscell.tubes[nt].lfline.y, fignum, color='g')

        # Plot all other lines
        for thistube in thiscell.tubes: 
            PlotPolygons2DQuiver(thistube.hfline.x, thistube.hfline.y, fignum, color='b')
            PlotPolygons2DQuiver(thistube.lfline.x, thistube.lfline.y, fignum, color='b')
        
    # Check bounds
    
    for i in np.arange(0, topomesh.face.ntot, 1):
        xb[0] = np.min([xb[0], np.min(topomesh.face.data[i].x)])
        xb[1] = np.max([xb[1], np.max(topomesh.face.data[i].x)])
        yb[0] = np.min([yb[0], np.min(topomesh.face.data[i].y)])
        yb[1] = np.max([yb[1], np.max(topomesh.face.data[i].y)])

    # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('topomesh face vertex distribution')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

# Grid generation data plotting: vertices
def PlotGridVertices(grid, fignum):
    # Plot vertices only
    
    # Initialize plotting bounds
    xb = [min(grid.vert.x), max(grid.vert.x)]
    yb = [min(grid.vert.y), max(grid.vert.y)]

    # Make plot
    PlotPoints2D(grid.vert.x, grid.vert.y, fignum, color='k', marker='.')

    # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

# Grid generation data plotting: faces
def PlotGridFaces(grid, fignum):
    # Plot faces only

    # Initialize plotting bounds
    xb = [min(grid.vert.x), max(grid.vert.x)]
    yb = [min(grid.vert.y), max(grid.vert.y)]

    # Construct aligned and non-aligned face coordinates
    try:
        aligned = np.nonzero(grid.face.aligned)
        aligned = aligned[0]
        nonaligned = np.nonzero(grid.face.aligned-1) # a bit of a hack but eh
        nonaligned = nonaligned[0]
    except:
        aligned = np.arange(0, grid.face.ntot, 1, dtype=int)
        nonaligned = aligned
    xfal = np.zeros(len(aligned)*3, dtype=float)
    yfal = np.zeros(len(aligned)*3, dtype=float)
    xfnonal = np.zeros(len(nonaligned)*3, dtype=float)
    yfnonal = np.zeros(len(nonaligned)*3, dtype=float)

    cc = 0
    for i in aligned: 
        xfal[3*cc] = grid.vert.x[grid.face.v1[i]-1]
        xfal[3*cc+1] = grid.vert.x[grid.face.v2[i]-1]
        xfal[3*cc+2] = np.nan
        yfal[3*cc] = grid.vert.y[grid.face.v1[i]-1]
        yfal[3*cc+1] = grid.vert.y[grid.face.v2[i]-1]
        yfal[3*cc+2] = np.nan
        cc = cc + 1 

    cc = 0
    for i in nonaligned: 
        xfnonal[3*cc] = grid.vert.x[grid.face.v1[i]-1]
        xfnonal[3*cc+1] = grid.vert.x[grid.face.v2[i]-1]
        xfnonal[3*cc+2] = np.nan
        yfnonal[3*cc] = grid.vert.y[grid.face.v1[i]-1]
        yfnonal[3*cc+1] = grid.vert.y[grid.face.v2[i]-1]
        yfnonal[3*cc+2] = np.nan
        cc = cc + 1 
    
    # Plot faces
    PlotPolygons2D(xfal, yfal, fignum, color='r', marker='', linewidth=0.25)
    PlotPolygons2D(xfnonal, yfnonal, fignum, color='k', marker='', linewidth=0.25)

    # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

# Grid generation data plotting: cells
def PlotGridCells(grid, fignum):
    # Plot cells only

    # Initialize plotting bounds
    xb = [min(grid.vert.x), max(grid.vert.x)]
    yb = [min(grid.vert.y), max(grid.vert.y)]

    # Construct cell coordinates
    xc = np.zeros(grid.cell.nvert + 2*grid.cell.ntot, dtype=float)
    yc = np.zeros(grid.cell.nvert + 2*grid.cell.ntot, dtype=float)

    counter = 0
    for i in np.arange(0, grid.cell.ntot): 
        nvc = grid.cell.vp2[i]
        tv = grid.cell.GetVert(i)-1
        xc[counter:counter+nvc] = grid.vert.x[tv]
        xc[counter+nvc] = grid.vert.x[tv[0]]
        xc[counter+nvc+1] = np.nan
        yc[counter:counter+nvc] = grid.vert.y[tv]
        yc[counter+nvc] = grid.vert.y[tv[0]]
        yc[counter+nvc+1] = np.nan

        counter = counter + nvc + 2
    
    # Plot
    PlotPolygons2D(xc, yc, fignum, color='k', marker='', linewidth=0.25)

     # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

# Grid generation data plotting: cells and indicated aligned faces
def PlotGridCellsAlignedFaces(grid, fignum):
    # Plot faces only
    PlotGridFaces(grid, fignum)

    # Plot vertexID
    #for i in np.arange(0, grid.vert.ntot, 1):
    #    PlotPoints2D(grid.vert.x[i], grid.vert.y[i], fignum, marker='.', color='r', label=str(grid.vert.ID[i]))


    # Plot aligned faces
    xf = 0.5*(grid.vert.x[grid.face.v1-1] + grid.vert.x[grid.face.v2-1])
    yf = 0.5*(grid.vert.y[grid.face.v1-1] + grid.vert.y[grid.face.v2-1]) 
    for i in np.arange(0, grid.face.ntot, 1): 
        if (grid.face.aligned[i] == 1) :
            PlotPoints2D(xf[i], yf[i], fignum, marker='.', color='b')

    # Plot cells without aligned faces
    #for i in np.arange(0, grid.cell.ntot, 1):
    #    s = grid.cell.faceP1[i]-1
    #    n = grid.cell.faceP2[i]
    #    if (n == 1) :
    #        cells = grid.cell.face[s]
    #    else: 
    #        cells = grid.cell.face[s:s+n]
    #    if all(j == 0 for j in grid.face.aligned[fcs]):
    #       PlotPoints2D(grid.cell.x[i], grid.cell.y[i], fignum, marker='.', color='g') 


# Grid generation data plotting: cells and indicated boundary faces
def PlotGridCellsBoundaryFaces(grid, fignum):
    # Plot faces only
    PlotGridFaces(grid, fignum)

    # Plot vertexID
    #for i in np.arange(0, grid.vert.ntot, 1):
    #    PlotPoints2D(grid.vert.x[i], grid.vert.y[i], fignum, marker='.', color='r', label=str(grid.vert.ID[i]))


    # Plot aligned faces
    xf = 0.5*(grid.vert.x[grid.face.v1-1] + grid.vert.x[grid.face.v2-1])
    yf = 0.5*(grid.vert.y[grid.face.v1-1] + grid.vert.y[grid.face.v2-1]) 
    for i in np.arange(0, grid.face.ntot, 1):
        if (not(grid.face.label[i] ==  0)) :
            PlotPoints2D(xf[i], yf[i], fignum, marker='.', color='b')

def PlotGridCellCutcells(grid,cctria,cctrapsP1,cctrapsP2,cctraps, fignum ):
    # Plot cells only
    PlotGridCells(grid, fignum)

    # Plot cutcells
    for i in np.arange(0,len(cctria), 1):
        PlotPoints2D(grid.cell.x[cctria[i]-1], grid.cell.y[cctria[i]-1], fignum, marker='.', color='b')

    for i in np.arange(0,len(cctrapsP1), 1):
        s = cctrapsP1[i]-1
        n = cctrapsP2[i]
        if (n == 1) :
            cells = cctraps[s]
        else: 
            cells = cctraps[s:s+n]
        PlotPoints2D(grid.cell.x[cells-1], grid.cell.y[cells-1], fignum, marker='.', color='g')

def PlotGridCellsFsFc(grid, fsFc, fsFcP1, fsFcP2, fignum):
    # Plot cells only
    PlotGridFaces(grid, fignum)

    # Plot cutcells
    xf = 0.5*(grid.vert.x[grid.face.v1-1] + grid.vert.x[grid.face.v2-1])
    yf = 0.5*(grid.vert.y[grid.face.v1-1] + grid.vert.y[grid.face.v2-1]) 
    for i in np.arange(0,len(fsFcP1), 1):
        s = fsFcP1[i]-1
        n = fsFcP2[i]
        if (n == 1) :
            faces = fsFc[s]
        else: 
            faces = fsFc[s:s+n]
        PlotPoints2D(xf[faces-1], yf[faces-1], fignum, marker='.', color='g')    

def PlotGridCellsFtCv(grid, ftcv, ftcvP1, ftcvP2, fignum):
    # Plot cells only
    PlotGridFaces(grid, fignum)

    # Plot cutcells 
    color_array = ['g', 'b', 'r', 'b']
    for i in np.arange(0,len(ftcvP1), 1):
        t = i % 4
        s = ftcvP1[i]-1
        n = ftcvP2[i]
        if (n == 1) :
            cells = ftcv[s]
        else: 
            cells = ftcv[s:s+n]
        PlotPoints2D(grid.cell.x[cells-1], grid.cell.y[cells-1], fignum, marker='.', color=color_array[t])

def PlotGridCellValue(grid, array, threshold, fignum):
    # Plot cells only
    PlotGridCells(grid, fignum)

    for i in np.arange(0, grid.cell.ntot, 1):
        if (array[i] >  threshold) :
            PlotPoints2D(grid.cell.x[i], grid.cell.y[i], fignum, marker='.', color='b')

def PlotGridCellArray(grid, array, fignum):
    # Plot cells only
    PlotGridCells(grid, fignum)  

    # Plot cutcells
    for i in np.arange(0,len(array), 1):
        PlotPoints2D(grid.cell.x[array[i]-1], grid.cell.y[array[i]-1], fignum, marker='.', color='b') 

def PlotGridCellVertArray(grid, array, fignum):
    # Plot cells only
    PlotGridCells(grid, fignum)  

    # Plot verts
    for i in np.arange(0,len(array), 1):
        PlotPoints2D(grid.vert.x[array[i]-1], grid.vert.y[array[i]-1], fignum, marker='.', color='b') 

# Grid topological data
def PlotGridTopologicalData(grid, fignum):
    # Description
    #------------
    # Plot a grid and highlight the topological data 
    
    # Plot the grid
    PlotGridFaces(grid, fignum)

    # Plot x-points
    if (grid.topodata.nxp > 0):
        xpID = grid.topodata.XpointID-1
        PlotPoints2D(grid.vert.x[xpID], grid.vert.y[xpID], fignum, marker='x', color='r', label='x-points')

    # Plot o-points
    if (grid.topodata.nop > 0):
        xpID = grid.topodata.OpointID-1
        PlotPoints2D(grid.vert.x[xpID], grid.vert.y[xpID], fignum, marker='o', color='r', label='o-points')

    # Plot strike points
    if (grid.topodata.nsp > 0):
        xpID = grid.topodata.SpointID-1
        PlotPoints2D(grid.vert.x[xpID], grid.vert.y[xpID], fignum, marker='^', color='r', label='strike points')

    # Plot tangency points
    if (grid.topodata.ntp > 0):
        xpID = grid.topodata.TpointID-1
        PlotPoints2D(grid.vert.x[xpID], grid.vert.y[xpID], fignum, marker='s', color='r', label='tangency points')

    # Plot divertor targets
    if (grid.topodata.ndiv > 0):
        xf = 0.5*(grid.vert.x[grid.face.v1-1] + grid.vert.x[grid.face.v2-1])
        yf = 0.5*(grid.vert.y[grid.face.v1-1] + grid.vert.y[grid.face.v2-1])
        for i in np.arange(0, grid.topodata.ndiv, 1):
            tdivfaces = grid.topodata.GetDivFace(i)
            PlotPolygons2D(xf[tdivfaces-1], yf[tdivfaces-1], fignum, label='divertor target ' + str(i))

# Grid generation data plotting: cells
def PlotTriaCells(tria, fignum):
    # Plot cells only

    # Initialize plotting bounds
    xb = [min(tria.vert.x), max(tria.vert.x)]
    yb = [min(tria.vert.y), max(tria.vert.y)]

    # Construct cell coordinates
    xc = np.zeros(tria.cell.nvert + 2*tria.cell.ntot, dtype=float)
    yc = np.zeros(tria.cell.nvert + 2*tria.cell.ntot, dtype=float)

    counter = 0
    for i in np.arange(0, tria.cell.ntot): 
        tv = tria.cell.vert[i,:]-1
        xc[counter:counter+3] = tria.vert.x[tv]
        xc[counter+3] = tria.vert.x[tv[0]]
        xc[counter+4] = np.nan
        yc[counter:counter+3] = tria.vert.y[tv]
        yc[counter+3] = tria.vert.y[tv[0]]
        yc[counter+4] = np.nan

        counter = counter + 5
    
    # Plot
    PlotPolygons2D(xc, yc, fignum, color='k', marker='', linewidth=0.25)

     # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

def PlotTriaCellVertArray(tria, array, fignum):
    # Plot cells only
    PlotTriaCells(tria, fignum)

    # Plot verts
    for i in np.arange(0,len(array), 1):
        PlotPoints2D(tria.vert.x[array[i]-1], tria.vert.y[array[i]-1], fignum, marker='.', color='b')     

#--------------------------------------------------------------------------#
#                             Shape Optimization                           #
#--------------------------------------------------------------------------#

def PlotFixedVesselPointsConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel vertices that have their position constrained. 

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + origvesselpolygonfile
    filepath = dirpath + filesep + fvpfile 

    # Get the grid data
    vals = dh.GetPolygonCoordinates(vertfilepath)
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
        facecolors='none', label='Vertices')

    # Special points
    try: 
        valscon = dh.GetVertexCoordinates(filepath)
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='*', label='Fixed vessel points')
    except:
        print('could not print fixed vessel point vertex constraints')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Fixed vessel points constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFixedVesselFluxConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel vertices that have their flux value constrained. 

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + origvesselpolygonfile
    filepath = dirpath + filesep + fvffile

    # Get the grid data
    vals = dh.GetPolygonCoordinates(vertfilepath)
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
        facecolors='none', label='Vertices')

    # Special points
    try: 
        valscon = dh.GetVertexCoordinates(filepath)
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='*', label='Fixed vessel points')
    except:
        print('could not print fixed vessel point vertex constraints')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Fixed vessel flux constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotShapeOptimizationHistory(dirpath, fignum):
        # Description
    #------------
    # Plot the goat optimization history on a log(value)-iterate plot.
    # This includes the cost function, (max) value of 
    # constraints, the convergence history (infinity norm of 
    # lagrangian gradient), and the linesearch step length 

    # Read data
    #----------
    historypath = dirpath + filesep + shapeopthistoryfile
    [valnames, vals] = dh.ReadGeneralColumnwiseFloatData(historypath)

    # Hedge for negative numbers
    for j in range(len(vals[1,:])):
        for i in range(len(vals[:, j])):
            if (vals[i, j] <= 0):
                vals[i, j] = np.nan

    # Set figure
    #-----------
    plt.figure(fignum)

    # Plot data
    #----------
    # First entry should be iteration counter, then convnorm, dphi, L, 
    # J, max(G), max(H), rxf, step, tol, ...

    # Convnorm
    plt.plot(vals[:, 0], vals[:, 1], 'rx-', label='max(abs(grad L))')
    
    # Cost function
    plt.plot(vals[:, 0], vals[:, 4], 'bx-', label='J')

    # Equality constraints
    plt.plot(vals[:, 0], vals[:, 5], 'gx-', label='max(G)')

    # Inequality constraints
    plt.plot(vals[:, 0], vals[:, 6], 'mx-', label='max(H)')

    # Line search step length
    plt.plot(vals[:, 0], vals[:, 8], 'kx-', label='alpha_ls')

    # Set figure data
    #----------------
    # Set axes
    SetAxesLimitsLogplot(plt.gca(), vals[:, 0], vals[:, [1, 4, 5, 6, 8]])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Shape optimization convergence history')
    thisaxes.set_xlabel('iteration number')
    thisaxes.set_ylabel('value')
    thisaxes.set_yscale('log')
    thisaxes.legend(loc='upper right')

def PlotVesselDisplacement(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel displacement between original and new vessel

    # Get filepath
    filepathorig = dirpath + filesep + origvesselpolygonfile 
    filepathnew = dirpath + filesep + vesselpolygonfile 

    # Get data
    valsorig = dh.GetPolygonCoordinates(filepathorig)
    valsnew = dh.GetPolygonCoordinates(filepathnew)

    # Plot both polygons
    PlotPolygons2D(valsorig[:, 0], valsorig[:, 1], fignum, color='b', marker='x',
        label='Original vessel polygon')
    PlotPolygons2D(valsnew[:, 0], valsnew[:, 1], fignum, color='r', marker='o',
        label='New vessel polygon')

    # Plot displacement vector
    d = valsnew - valsorig 
    plt.quiver(valsorig[:, 0], valsorig[:, 1], d[:, 0], d[:, 1], color='g', angles='xy', scale_units='xy', scale=1)

    # Set axes
    SetAxesLimits2D(plt.gca(), valsnew[:, 0], valsnew[:, 1])
    
#--------------------------------------------------------------------------#
#                                 Vessel                                   #
#--------------------------------------------------------------------------#

def PlotInitialVesselPolygon(dirpath, fignum):
    # Description
    # ------------
    # Plot the initial vessel polygon

    # Get filepath
    filepath = dirpath + filesep + origvesselpolygonfile 

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='b', marker='x',
                   label='Vessel polygon')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vessel polygon')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotVesselPolygon(dirpath, fignum):
    # Description
    # ------------
    # Plot the grid cells once by reading in the grid cell polygon data as
    # written out in the cells.dat file

    # Get filepath
    filepath = dirpath + filesep + vesselpolygonfile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
                   label='Vessel polygon')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vessel polygon')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotVesselIterate(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel polygon obtained during a shape iteration
    # Get filepath
    filepath = dirpath + filesep + currentvesselpolygonfile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
                   label='Vessel polygon')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vessel polygon')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def MonitorGridAndVessel(datadir, num, pausetime, maxruntime):
    # Description
    #------------
    # This routine makes a plot that is continuously updated of the
    # current grid and vessel

    # Time interval to replot [s]
    starttime = time.time()

    # Prepare for gui event loop
    plt.ion()

    # Plot the data
    PlotGridCellsIterate(datadir, 1)
    PlotVesselIterate(datadir, 1)
    thisfig = plt.gcf()

    # Loop until time has passed
    while (time.time() - starttime <= maxruntime):
        try: 
            # Plot the data
            PlotGridCellsIterate(datadir, 1)
            PlotVesselIterate(datadir, 1)

            # Draw
            thisfig.canvas.draw()
            thisfig.canvas.flush_events()

            # Pause
            time.sleep(pausetime)

            # Clear figure
            ClearCurrentAxes()
        except: 
            time.sleep(pausetime)

def PlotStructure(structure, fignum, **plotargs):
    # Description
    #------------
    # Plot the structure polygons
    minx = np.inf
    maxx = -np.inf
    miny = np.inf
    maxy = -np.inf
    cc = 1
    for i in structure:
        PlotPolygons2D(i.x, i.y, fignum, **plotargs, label='structure ' + str(i.ID))
        ID = np.arange(1, len(i.x)+1, 1) 
        PlotPoints2DWithID(i.x, i.y, ID, fignum, **plotargs)
        minx = np.min([minx, np.min(i.x)])
        maxx = np.max([maxx, np.max(i.x)])
        miny = np.min([miny, np.min(i.y)])
        maxy = np.max([maxy, np.max(i.y)])
        cc = cc + 1

    # Set limits
    # Set axes
    SetAxesLimits2D(plt.gca(), [minx, maxx], [miny, maxy])
    thisaxes = plt.gca()
    thisaxes.legend()

def PlotStructureFromFile(dirpath, fignum):
    # Read the structure
    structure = dh.ReadStructureFile(dirpath)

    # Plot
    PlotStructure(structure, fignum, color='k', linewidth='2')

#--------------------------------------------------------------------------#
#                             Magnetic field                               #
#--------------------------------------------------------------------------#
def PlotMagneticFieldFromRZPsi(dirpath, fignum, plottype, **plotargs):
    # Description
    #------------
    # Plot magnetic field that is stored in an rzpsi.dat file format. 
    # Use 'plottype' to determine what kind of plot to make (with
    # arguments defined in plotargs):
    # 'contour'     classic contour plot
    # 'contourf'    filled contours
    # 'quiver'      quiver plot of Bx, By (approximated using FD). Note
    #               that these arrows are normalized
    
    # Load values
    [R, Z, Psi] = dh.ReadRZPsiFile(dirpath)

    if plottype == 'contour':
        PlotStructured2DContour(R, Z, Psi, fignum, **plotargs)
    elif plottype == 'contourf':
        PlotStructured2DContourf(R, Z, Psi, fignum, **plotargs)
    elif plottype == 'quiver':
        bx = np.gradient(Psi.transpose(), Z, axis=0)
        by = np.gradient(Psi.transpose(), R, axis=1)
        bx = -bx
        bn = np.sqrt(bx**2 + by**2)
        bx = np.divide(bx, bn) 
        by = np.divide(by, bn)
        PlotStructured2DQuiver(R, Z, bx, by, fignum, **plotargs)
    else:
        raise('PlotMagneticfieldFromRZPsi: unknown method')

#--------------------------------------------------------------------------#
#                               2D surface plots                           #
#--------------------------------------------------------------------------#
    
# Wrappers 
def Plot2DSurfaceData(filepath, fignum):
    # Description
    #------------
    # Read and plot 2D data as patchplot

    # Read
    data = dh.GetGeneral2DSurfacedata(filepath)

    # Plot
    PlotGeneral2DSurface(data[:, 0], data[:, 1], data[:, 2], fignum)

def Plot2DSurfaceDataContourf(filepath, fignum, **plotargs):
    # Description
    #------------
    # Read and plot 2D data with filled contours

    # Read
    data = dh.GetGeneral2DSurfaceData(filepath)

    # Plot
    PlotGeneral2DContourf(data[:, 0], data[:, 1], data[:, 2], fignum, **plotargs)

    # Add colorbar

def Plot2DSurfaceDataContour(filepath, fignum, **plotargs):
    # Description
    #------------
    # Read and plot 2D data with filled contours

    # Read
    data = dh.GetGeneral2DSurfaceData(filepath)

    # Plot
    PlotGeneral2DContour(data[:, 0], data[:, 1], data[:, 2], fignum, **plotargs)

    # Add colorbar

#--------------------------------------------------------------------------#
#                                 Interpolant                              #
#--------------------------------------------------------------------------#

def VisualizeGridInterpolant2D(interp, fignum):
    # Description
    #------------
    # Simple routine to visualize the grid interpolant evaluated on 
    # some predefined raster

    # Construct grid vectors
    xb = np.array([np.min(interp.grid.vert.x), np.max(interp.grid.vert.x)])
    yb = np.array([np.min(interp.grid.vert.y), np.max(interp.grid.vert.y)]) 
    resx = 100
    resy = 100
    x = np.linspace(xb[0], xb[1], resx)
    y = np.linspace(yb[0], yb[1], resy)
    val = np.zeros((len(x), len(y)), dtype=float)
    for i in np.arange(0, len(x), 1): 
        for j in np.arange(0, len(y), 1):
            val[i, j] = interp.Evaluate(x[i], y[j])

    # Visualize
    PlotStructured2DSurface(x, y, val, fignum)


#==========================================================================#
#                                                                          #
#                           GENERAL PLOTTING ROUTINES                      #
#                                                                          #
#==========================================================================#
def PlotPolygons2D(x, y, fignum, **plotargs):
    # General polygon plotter. x and y should be np.arrays containing the
    # coordinates to plot. Fignum should contain the figure number on which
    # to plot the data.

    # Set the current figure
    plt.figure(fignum)

    # Plot the data as a polygon plot
    #fig, ax = plt.subplots()
    plt.plot(x, y, **plotargs)
    
    #myrange = range(len(x))
    #for i, txt in enumerate(myrange):
    #    ax.annotate(txt, (x[i], y[i]))
    #plt.draw()

def PlotFilledPolygons2D(x, y, fignum, **plotargs):
    # General polygon plotter. x and y should be np.arrays containing the
    # coordinates to plot. Fignum should contain the figure number on which
    # to plot the data.

    # Set the current figure
    plt.figure(fignum)

    # Plot the data as a polygon plot
    #fig, ax = plt.subplots()
    plt.fill(x, y, **plotargs)
    
    #myrange = range(len(x))
    #for i, txt in enumerate(myrange):
    #    ax.annotate(txt, (x[i], y[i]))
    #plt.draw()

def PlotPolygons2DQuiver(x, y, fignum, **plotargs):
    # Same as PlotPolygons2D, but now we plot arrows between the 
    # different nodes according to the polygon orientation

    # Set the current figure
    plt.figure(fignum)

    # Plot the data as a polygon plot
    #fig, ax = plt.subplots()
    plt.plot(x, y, **plotargs)

    # Plot displacement vector
    dx = -(x[0:len(x)-1] - x[1:len(x)])
    dy = -(y[0:len(y)-1] - y[1:len(y)])
    plt.quiver(x[0:len(x)-1], y[0:len(y)-1], dx, dy, 
        **plotargs, angles='xy', scale_units='xy', scale=1)

def PlotVertexPairs2DQuiver(x1, y1, x2, y2, fignum, **plotargs):
    # Plot an arrow between vertex pairs defined by (x1, y1), (x2, y2)

    # Set the current figure
    plt.figure(fignum)

    # Plot displacement vector
    dx = (x2 - x1)
    dy = (y2 - y1)
    plt.quiver(x1, y1, dx, dy, 
        **plotargs, angles='xy', scale_units='xy', scale=1)

def PlotPolygonData(filepath, fignum, **plotargs):
    # Read data
    data = dh.GetPolygonCoordinates(filepath)

    # Plot
    plt.figure(fignum)
    PlotPolygons2D(data[:, 0], data[:, 1], fignum, **plotargs)

def PlotVertexPairsData2D(filepath, fignum, **plotargs):
    # Read data
    data = dh.GetVertexPairCoordinates(filepath)

    # Compute faces
    xf = 0.5*(data[:, 0] + data[:, 2])
    yf = 0.5*(data[:, 1] + data[:, 3])

    # Plot
    PlotPoints2D(xf, yf, fignum, **plotargs)

def PlotPoints2D(x, y, fignum, **plotargs):
    # General point plotter.

    # Set the current figure
    plt.figure(fignum)
    plt.scatter(x, y, **plotargs)
    plt.draw()
    
def PlotPoints2DWithID(x, y, ID, fignum, **plotargs):
    # Set the current figure
    fig = plt.figure(fignum)
    plt.scatter(x, y, **plotargs)
    ax = fig.axes
    k = 0
    for i, txt in enumerate(ID):
        ax[0].text(x[k], y[k], str(txt), size=10)
        k = k + 1
    plt.draw()

def PlotGeneral2DSurface(x, y, z, fignum, **plotargs):
    # General z = f(x, y) surface plotter - may be expensive since
    # a triangulation is created under the hood. 

    # Set the current figure
    plt.figure(fignum)
    plt.tripcolor(x, y, z, **plotargs)
    plt.draw()

def PlotGeneral2DContourf(x, y, z, fignum, **plotargs):
    # Same as PlotGeneral2DSurface, but now we plot filled contours
    
    # Set the current figure
    fig = plt.figure(fignum)
    CS = plt.tricontourf(x, y, z, **plotargs)
    cbar = fig.colorbar(CS)
    plt.draw()

def PlotGeneral2DContour(x, y, z, fignum, **plotargs):
    # Same as PlotGeneral2DSurface, but now we plot filled contours
    
    # Set the current figure
    fig = plt.figure(fignum)
    CS = plt.tricontour(x, y, z, **plotargs)
    cbar = fig.colorbar(CS)
    plt.draw()

def PlotGeneral2DPatch(verts, val, fignum):
    # Create a patchplot by creating a polygon collection

    # Set the current figure
    fig = plt.figure(fignum)

    # Create the polygon collection
    poly = mpl.collections.PolyCollection(verts, closed=True, array=val, edgecolor='k', linewidth=0.0, cmap='viridis')
    
    # Get the current axes
    ax = fig.axes 
    if len(ax) == 0:
        ax = fig.add_subplot()
    else:
        ax = ax[0] 

    # Add the collection
    ax.add_collection(poly, autolim=True)
    ax.autoscale_view()

    # Set colorbar
    fig.colorbar(poly, ax=ax)

def PlotStructured2DSurface(x, y, val, fignum, **plotargs):
    # Description
    #------------
    # Plot structured 2D data given by the values 'val', where val[i, j]
    # represents the value at x[i], y[j]. Note that under the hood, we 
    # simply call the more general plotting version. This is likely 
    # inefficient but I don't care

    # Construct plotting points
    xp = np.zeros(val.size)
    yp = np.zeros(val.size)
    vp = np.zeros(val.size)
    k = 0
    for j in np.arange(val.shape[1]):
        for i in np.arange(val.shape[0]):
            xp[k] = x[i]
            yp[k] = y[j]
            vp[k] = val[i, j]
            k = k + 1

    # Plot
    PlotGeneral2DSurface(xp, yp, vp, fignum, **plotargs)

def PlotStructured2DContour(x, y, val, fignum, **plotargs):
    # Description
    #------------
    # Plot structured 2D data given by the values 'val', where val[i, j]
    # represents the value at x[i], y[j]. Note that under the hood, we 
    # simply call the more general plotting version. This is likely 
    # inefficient but I don't care

    # Construct plotting points
    xp = np.zeros(val.size)
    yp = np.zeros(val.size)
    vp = np.zeros(val.size)
    k = 0
    for j in np.arange(val.shape[1]):
        for i in np.arange(val.shape[0]):
            xp[k] = x[i]
            yp[k] = y[j]
            vp[k] = val[i, j]
            k = k + 1

    # Plot
    PlotGeneral2DContour(xp, yp, vp, fignum, **plotargs)

def PlotStructured2DContourf(x, y, val, fignum, **plotargs):
    # Description
    #------------
    # Plot structured 2D data given by the values 'val', where val[i, j]
    # represents the value at x[i], y[j]. Note that under the hood, we 
    # simply call the more general plotting version. This is likely 
    # inefficient but I don't care

    # Construct plotting points
    xp = np.zeros(val.size)
    yp = np.zeros(val.size)
    vp = np.zeros(val.size)
    k = 0
    for j in np.arange(val.shape[1]):
        for i in np.arange(val.shape[0]):
            xp[k] = x[i]
            yp[k] = y[j]
            vp[k] = val[i, j]
            k = k + 1

    # Plot
    PlotGeneral2DContourf(xp, yp, vp, fignum, **plotargs)

def PlotStructured2DQuiver(x, y, vx, vy, fignum, **plotargs):
    # Description
    #------------
    # Make a quiver plot from structured data. X, y should give the 
    # axis spacing, vx, vy should be nx-by-ny arrays containing the 
    # x and y vector components. 

    # Set the current figure
    plt.figure(fignum)

    # Plot
    plt.quiver(x, y, vx, vy, 
        **plotargs, angles='xy', scale_units='xy') # , angles='xy', scale_units='xy', scale=1

#==========================================================================#
#                                                                          #
#                                 AUXILIARY                                #
#                                                                          #
#==========================================================================#
def ShowFigures(*args, **kwargs):
    # Just a wrapper for plt.show()
    plt.show(*args, **kwargs)

def ClearCurrentAxes():
    # Just a wrapper for plt.clf()
    plt.clf()

def SetAxesLimits2D(thisaxes, xdata, ydata):
    # Automatically set the axes limits based on the x and y figure data.
    # Also applies true scaling

    # Compute data limits - ignore NaNs
    maxx = np.nanmax(xdata)
    minx = np.nanmin(xdata)
    maxy = np.nanmax(ydata)
    miny = np.nanmin(ydata)

    # Compute ranges and add some buffer to the sides
    xrange = maxx - minx
    yrange = maxy - miny
    xbuffer = 0.05*xrange
    ybuffer = 0.05*yrange

    # Compute limits
    xlim = [minx - xbuffer, maxx + xbuffer]
    ylim = [miny - ybuffer, maxy + ybuffer]

    # Set limits
    thisaxes.set_xlim(xlim)
    thisaxes.set_ylim(ylim)

    # Set proper scaling
    thisaxes.set_aspect(1)

def SetAxesLimitsLogplot(thisaxes, xdata, ydata):
    # Automatically set the axes limits based on the x and y figure data.
    # Also applies true scaling

    # Compute data limits - ignore NaNs
    maxx = np.nanmax(xdata)
    minx = np.nanmin(xdata)
    maxy = np.nanmax(ydata)
    miny = np.nanmin(ydata)


    # Compute limits
    xlim = [minx, maxx]
    ylim = [miny, maxy]

    # Set limits
    thisaxes.set_xlim(xlim)
    thisaxes.set_ylim(ylim)

    # Set proper scaling
    thisaxes.set_aspect(1)

def GetColorsFromValue(val, minval, maxval):
    # Set colormap
    colormaptype = 'viridis'
    cm = mpl.colormaps[colormaptype]

    # Compute the color values
    norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    valc = norm(val)
    this = mpl.cm.ScalarMappable(norm=norm, cmap=cm)
    col = this.to_rgba(valc, alpha=None, bytes=False, norm=True)

    return col

def PlotCellBasedQuantity2D(grid, val, fignum, bounds=[-np.inf, np.inf], doguardcells=True):
    # Description
    #------------
    # Make a patchplot of a cell based quantity. If guard cells are 
    # present, we need to create additional guard vertices in order to
    # visualize the values there. These guard cells will not be be drawn
    # to scale but will have the width defined in GUARD_CELL_WIDTH.
    # It is assumed in this case that face and cell coordinates are
    # available. 

    # Check
    if (len(val) != grid.cell.ntot):
        raise ValueError('PlotCellBasedQuantity2D: ' \
            'value length is not equal to number of grid cells')
        
    # Construct cell polygon collection
    verts = []
    if (hasattr(grid.cell, 'nci')):
        # Interior cells
        for i in np.arange(0, grid.cell.nci): 
            nvc = grid.cell.vp2[i]
            tv = grid.cell.GetVert(i)-1
            
            verts.append(list(zip(grid.vert.x[tv], grid.vert.y[tv])))

        # Guard cells
        if doguardcells:
            for i in np.arange(grid.cell.nci, grid.cell.ntot): 
                # Get the guard cell face
                tf = grid.cell.GetFace(i)
                if len(tf) == 1:
                    pass 
                else:
                    raise ValueError('PlotCellBasedQuantity2D: guard cell does not ' \
                    'have exactly one face, unexpected.') 
                tfind = tf[0]-1

                # Construct the cell connector
                nb1 = grid.face.nb1[tfind]
                nb2 = grid.face.nb2[tfind]
                if nb1-1 == i: 
                    nbind = nb2-1
                elif nb2-1 == i:
                    nbind = nb1-1
                else:
                    raise ValueError('PlotCellBasedQuantity2D: guard face ' \
                    'does not have the guard cell as neighbour, unexpected.')
                ccx = grid.face.x[tfind] - grid.cell.x[nbind]
                ccy = grid.face.y[tfind] - grid.cell.y[nbind]
                ccn = np.sqrt(ccx**2 + ccy**2)
                ccx = ccx/ccn 
                ccy = ccy/ccn 

                # Construct the outward facing face normal
                tfv1ind = grid.face.v1[tfind]-1
                tfv2ind = grid.face.v2[tfind]-1
                nx = -(grid.vert.y[tfv2ind] - grid.vert.y[tfv1ind])
                ny = grid.vert.x[tfv2ind] - grid.vert.x[tfv1ind]
                nn = np.sqrt(nx**2 + ny**2)
                nx = nx/nn 
                ny = ny/nn
                if (nx*ccx + ny*ccy) < 0:
                    nx = -nx 
                    ny = -ny

                # Get the boundary vertices from the guard cell
                nvc = grid.cell.vp2[i]
                tv = grid.cell.GetVert(i)-1

                # Project the vertices along the connector vector
                tvx = grid.vert.x[tv]
                tvy = grid.vert.y[tv]
                tgvx = tvx + nx*GUARD_CELL_WIDTH 
                tgvy = tvy + ny*GUARD_CELL_WIDTH
                tx = np.append(tvx, np.flipud(tgvx))
                ty = np.append(tvy, np.flipud(tgvy))
                
                # Add vertex coordinates
                verts.append(list(zip(tx, ty)))

        

    else:
        for i in np.arange(0, grid.cell.ntot): 
            nvc = grid.cell.vp2[i]
            tv = grid.cell.GetVert(i)-1
            
            verts.append(list(zip(grid.vert.x[tv], grid.vert.y[tv])))

    # Make patchplot
    thisval = copy.deepcopy(val)
    thisval = np.where(thisval > bounds[1], np.repeat(bounds[1], len(val), 0), thisval)
    thisval = np.where(thisval < bounds[0], np.repeat(bounds[0], len(val), 0), thisval)
    PlotGeneral2DPatch(verts, thisval, fignum)

    # Set axes
    SetAxesLimits2D(plt.gca(), grid.cell.x, grid.cell.y)

def PlotVertBasedQuantity2D(grid, val, fignum):

# Description
    #------------
    # Make a patchplot of a cell based quantity

    # Check
    if (len(val) != grid.vert.ntot):
        raise ValueError('PlotVertBasedQuantity2D: ' \
            'value length is not equal to number of grid vertices')
        
    
    # Construct cell polygon collection
    verts = []
    cval = np.zeros(grid.cell.ntot, dtype=float)

    counter = 0
    for i in np.arange(0, grid.cell.ntot): 
        nvc = grid.cell.vp2[i]
        tv = grid.cell.GetVert(i)-1
        
        verts.append(list(zip(grid.vert.x[tv], grid.vert.y[tv])))
        cval = verts.append
        counter = counter + nvc + 2

    # Make patchplot
    PlotGeneral2DPatch(verts, val, fignum)

    # Set axes
    SetAxesLimits2D(plt.gca(), grid.cell.x, grid.cell.y)

def PlotTMCellBasedQuantity(topomesh, val, fignum):
    # Description
    #------------
    # Make a patchplot of a cell based quantity

    # Check
    if (len(val) != topomesh.cell.ntot):
        raise ValueError('PlotCellBasedQuantity2D: ' \
            'value length is not equal to number of grid cells')
        
    
    # Construct cell polygon collection
    verts = []

    for i in np.arange(0, topomesh.cell.ntot): 
        
        verts.append(list(zip(topomesh.cell.data[i].x, topomesh.cell.data[i].y)))

    # Make patchplot
    PlotGeneral2DPatch(verts, val, fignum)

    # Set axes
    SetAxesLimits2D(plt.gca(), topomesh.vert.x, topomesh.vert.y)





