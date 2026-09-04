# This module handles plotting data. It reads in data and returns the data
# in proper structures.
import numpy as np
import os 
from . import goat_types as gt

#--------------------------------------------------------------------------#
#                                   goat                                   #
#--------------------------------------------------------------------------#
def GetFolder():
	# Description
	#------------

	# Check
	# Set base_dir, should not be changed
	bd = './'

	# Get all subdirectories in base_dir
	subfolders = [f for f in os.listdir(bd) if os.path.isdir(os.path.join(bd, f))]

	# Check if there are any subfolders
	if not subfolders:
		print("No subfolders found in Runs/. Exiting.")
		exit()

	# Print the subfolders with numbers
	print("Select a folder:")
	for idx, folder in enumerate(subfolders, start=1):
		print(f"{idx}) {folder}")

	# User input selection
	while True:
		try:
			choice = int(input("Enter the number of the folder: "))
			if 1 <= choice <= len(subfolders):
				sd = subfolders[choice - 1]
				break
			else:
				print("Invalid number. Please enter a valid option.")
		except ValueError:
			print("Invalid input. Please enter a number.")

	# Print and store the selected folder name
	print(f"You selected: {sd}")

	# Set base_dir, should not be changed
	base_dir = sd

	# Get all subdirectories in base_dir
	subfolders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]

	# Check if there are any subfolders
	if not subfolders:
		print("No subfolders found in Runs/. Exiting.")
		exit()

	# Print the subfolders with numbers
	print("Select a folder:")
	for idx, folder in enumerate(subfolders, start=1):
		print(f"{idx}) {folder}")

	# User input selection
	while True:
		try:
			choice = int(input("Enter the number of the folder: "))
			if 1 <= choice <= len(subfolders):
				s_dir = subfolders[choice - 1]
				break
			else:
				print("Invalid number. Please enter a valid option.")
		except ValueError:
			print("Invalid input. Please enter a number.")

	# Print and store the selected folder name
	print(f"You selected: {s_dir}")

	datadir = './' + sd + '/' + s_dir
	print(f'Final dir: {datadir}')
	return datadir

def GetDataDirectory():
    # Description
    #------------
    # This routine determines based on the environment variables which
    # data directory has to be considered. For goat not coupled to 
    # solps, this is simply the current directory '.', when coupled to
    # solps, all output should be written to './output'.

    # Check
    if 'SOLPSTOP' in os.environ:
        datadir = './output'
    else:
        datadir = './output'

    return datadir

def GetVertexCoordinates(filepath):
    # Description
    #------------
    # This routine reads a formatted file where the vertex IDs and
    # coordinates (2D) are stored as [ID, x, y] columns. White/empty lines
    # are ignored.

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 3])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            pass
            # Don't read in
        else:
            # Read
            vals[cc, 0:3] = np.fromstring(i, dtype=float, count=3, sep=' ')
            cc = cc + 1

    # Return values
    return vals[0:cc, 1:3]

def GetVertexCoordinatesWithID(filepath):
    # Description
    #------------
    # This routine reads a formatted file where the vertex IDs and
    # coordinates (2D) are stored as [ID, x, y] columns. White/empty lines
    # are ignored.

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 2], dtype=float)
    IDs = np.zeros([len(alllines), 1], dtype=int)

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            pass
            # Don't read in
        else:
            # Read
            values = i.split()
            IDs[cc] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
            vals[cc, 0] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
            vals[cc, 1] = np.fromstring(values[2], dtype=float, count=1, sep=' ')
            cc = cc + 1

    # Return values
    return IDs, vals

def GetPolygonCoordinates(filepath):
    # Description
    #------------
    # This routine reads a formatted file where the polygon
    # coordinates (2D) are stored as [x, y] columns. White/empty lines
    # indicate the end of a polygon and are padded with NaNs in the
    # data.

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 2])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            vals[cc, :] = np.nan
        else:
            # Read
            vals[cc, 0:2] = np.fromstring(i, dtype=float, count=2, sep=' ')

        # Update counter
        cc = cc + 1

    # Return values
    return vals[0:cc, 0:2]

def GetVertexPairCoordinates(filepath):
    # Description
    #------------
    # Read in a vertex pair format file where the vertices are stored in
    # ID1, ID2, x1, x2, y1, y2 format. The values that are returned are
    # in [x1, y1, x2, y2] format. Empty lines are skipped

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    if len(alllines) < 2:
        vals = np.zeros([0, 0])
        return vals 

    # Retrieve sizes
    dim = np.fromstring(alllines[1], dtype=float, sep=' ')
    del alllines[0]
    npoints = np.floor(len(dim)/3).astype(int)
    vals = np.zeros([len(alllines), len(dim)])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            # don't read
            pass
        else:
            # Read
            vals[cc, 0:6] = np.fromstring(i, dtype=float, count=6, sep=' ')

        # Update counter
        cc = cc + 1

    # Return values
    returnvec = np.zeros([2*npoints], dtype=int)
    for i in np.arange(0, npoints, 1):
        returnvec[2*i] = npoints + i 
        returnvec[2*i+1] = 2*npoints + i 
    return vals[0:cc, returnvec]

def GetVertexPairsPairCoordinates(filepath):
    # Description
    #------------
    # Read in a pair of vertex pairs format file where the vertices are stored in
    # ID1, ID2, ID3, ID4,, x1, x2, x3, x4, y1, y2, y3, y4 format. The values that are returned are
    # in [x1, y1, x2, y2, x3, y3, x4, y4] format. Empty lines are skipped

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Retrieve sizes
    dim = np.fromstring(alllines[1], dtype=float, sep=' ')
    del alllines[0]
    npoints = np.floor(len(dim)/3).astype(int)
    vals = np.zeros([len(alllines), len(dim)])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            # don't read
            pass
        else:
            # Read
            vals[cc, 0:12] = np.fromstring(i, dtype=float, count=12, sep=' ')

        # Update counter
        cc = cc + 1

    # Return values
    returnvec = np.zeros([2*npoints], dtype=int)
    for i in np.arange(0, npoints, 1):
        returnvec[2*i] = npoints + i 
        returnvec[2*i+1] = 2*npoints + i 
    return vals[0:cc, returnvec]

def GetGeneral2DSurfaceData(filepath):
    # Description
    #------------
    # Read in x, y, z data where x and y are assumed to be spatial 
    # coordinates and z is some field. Blank lines are allowed but 
    # skipped. It is assumed that the first line is a header. 

    # Read in vertex coordinates from the file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 3])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            pass
            # Don't read in
        else:
            # Read
            vals[cc, 0:3] = np.fromstring(i, dtype=float, count=3, sep=' ')
            cc = cc + 1

    # Return values
    return vals[0:cc, 0:3]

def ReadGeneralColumnwiseFloatData(filepath):
    # Description
    #------------
    # Read in any columnwise organized float data. It is assumed that 
    # the first line is a header with the same amount of columns as the
    # data later on, and that the string of each header describes 
    # the data. 

    # Read data
    #----------
    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Split the header
    valnames = alllines[0].split()

    # Delete the header
    del alllines[0]

    # Count the number of columns
    ncol = len(valnames)

    # Initialize
    vals = np.zeros([len(alllines), ncol])

    # Read data
    cc = 0
    for i in alllines:
        if i == '\n':  # empty string
            pass
            # Don't read in
        else:
            # Read
            vals[cc, 0:ncol] = np.fromstring(i, dtype=float, count=ncol, sep=' ')
            cc = cc + 1
        
    # Return 
    return valnames, vals

def ReadTopomeshFile(filepath):
    # Description
    #------------
    # Read in the data of a topological mesh by reading the topological
    # mesh file (format specific here)
    topomesh = gt.Topomesh()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read in vertex data
    #--------------------
    # Get the 'vertices' header
    i = 0
    while i < len(alllines): 
        if "vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read total amount of vertices
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=1, sep=' ')
    topomesh.vert.ntot = ntot[0]
    
    # Initialize fields
    topomesh.vert.Initialize(topomesh.vert.ntot)

    # Skip header
    i = i + 2

    # Read vertex data
    for j in np.arange(0, topomesh.vert.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.vert.ID[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        topomesh.vert.x[j] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
        topomesh.vert.y[j] = np.fromstring(values[2], dtype=float, count=1, sep=' ')
        topomesh.vert.type[j] = np.fromstring(values[3], dtype=int, count=1, sep=' ')
        topomesh.vert.fsID[j] = np.fromstring(values[4], dtype=int, count=1, sep=' ')
        topomesh.vert.fval[j] = np.fromstring(values[5], dtype=float, count=1, sep=' ')
        topomesh.vert.BV[j] = np.fromstring(values[6], dtype=int, count=1, sep=' ')

    # Read face data
    #---------------
    # Get the 'faces' header
    while i < len(alllines): 
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read total amount of faces
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=1, sep=' ')
    topomesh.face.ntot = ntot[0]
    
    # Initialize fields
    topomesh.face.Initialize(topomesh.face.ntot)

    # Skip header
    i = i + 2

    # Read vertex data
    for j in np.arange(0, topomesh.face.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.face.ID[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        topomesh.face.fsID[j] = np.fromstring(values[1], dtype=int, count=1, sep=' ')
        topomesh.face.type[j] = np.fromstring(values[2], dtype=int, count=1, sep=' ')
        topomesh.face.vert[j, 0] = np.fromstring(values[3], dtype=int, count=1, sep=' ')
        topomesh.face.vert[j, 1] = np.fromstring(values[4], dtype=int, count=1, sep=' ')
        topomesh.face.BF[j] = np.fromstring(values[5], dtype=int, count=1, sep=' ')
        topomesh.face.nc[j] = np.fromstring(values[6], dtype=int, count=1, sep=' ')

    # Update file position
    i = i + topomesh.face.ntot-1

    # Read face coordinates
    #----------------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "face coordinates" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read coordinates
    allfound = False 
    counter = 0
    while (i < len(alllines)) and (not allfound):
        # Read until we find 'face'
        if "face" in alllines[i]:
            # Get face ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1
            if (counter == topomesh.face.ntot):
                allfound = True 

            # Update position
            i = i + 1

            # Read nc coordinates
            topomesh.face.data[fID].Initialize(topomesh.face.nc[fID])
            for k in np.arange(0, topomesh.face.nc[fID], 1):
                values = alllines[i+k].split()
                this = np.fromstring(values[0], dtype=float, count=1)
                topomesh.face.data[fID].x[k] = np.fromstring(values[0], dtype=float, count=1, sep=' ')
                topomesh.face.data[fID].y[k] = np.fromstring(values[1], dtype=float, count=1, sep=' ')

            # Compute metrics
            topomesh.face.data[fID].ComputeMetrics()

            # Update position
            i = i + topomesh.face.nc[fID]
        else: 
            i = i + 1 

    # Read cell data
    #---------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1
    
    # Read total amount of cells
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=3, sep=' ')
    topomesh.cell.ntot = ntot[0]
    topomesh.cell.nvert = ntot[1]
    topomesh.cell.nface = ntot[2]
    
    # Initialize fields
    topomesh.cell.Initialize(topomesh.cell.ntot, topomesh.cell.nvert, topomesh.cell.nface)

    # Skip header
    i = i + 2

    # Read vertex data
    for j in np.arange(0, topomesh.cell.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.cell.ID[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        topomesh.cell.nc[j] = np.fromstring(values[1], dtype=int, count=1, sep=' ')

    # Update file position
    i = i + topomesh.cell.ntot-1

    # Get header position
    i = 0
    while i < len(alllines): 
        if "cell vertexlist" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read vertex data
    for j in np.arange(0, topomesh.cell.nvert, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.cell.vert[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')

    # Get header position
    i = 0
    while i < len(alllines): 
        if "cell vertex pointer" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read vertex data
    for j in np.arange(0, topomesh.cell.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        tv = np.fromstring(values[0], dtype=int, count=2, sep=' ')
        topomesh.cell.vertP[j, 0] = tv[0] - 1 # zero-based indexing!
        topomesh.cell.vertP[j, 1] = tv[1]

    # Read cell coordinates
    #----------------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "cell polygons" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read coordinates
    allfound = False 
    if topomesh.cell.ntot == 0:
        allfound = True
    counter = 0
    while (i < len(alllines)) and (not allfound):
        # Read until we find 'cell'
        if "cell" in alllines[i]:
            # Get cell ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1
            if (counter == topomesh.cell.ntot):
                allfound = True 

            # Update position
            i = i + 1

            # Read nc coordinates
            topomesh.cell.data[fID].Initialize(topomesh.cell.nc[fID])
            for k in np.arange(0, topomesh.cell.nc[fID], 1):
                values = alllines[i+k].split()
                this = np.fromstring(values[0], dtype=float, count=1)
                topomesh.cell.data[fID].x[k] = np.fromstring(values[0], dtype=float, count=1, sep=' ')
                topomesh.cell.data[fID].y[k] = np.fromstring(values[1], dtype=float, count=1, sep=' ')

            # Update position
            i = i + topomesh.cell.nc[fID]
        else: 
            i = i + 1 

    # Read flux surfaces
    #-------------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "flux surfaces" in alllines[i]:
            break 
        else: 
            i = i + 1
    
    # Read total amount of flux surfaces
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=1, sep=' ')
    topomesh.fs.ntot = ntot[0]
    
    # Initialize fields
    topomesh.fs.Initialize(topomesh.fs.ntot)

    # Read data if available
    if topomesh.fs.ntot > 0:
        # Skip header
        i = i + 2

        # Read data
        for j in np.arange(0, topomesh.fs.ntot, 1):
            # Split the string
            values = alllines[i+j].split()

            # Add the values
            topomesh.fs.ID[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
            topomesh.fs.psi[j] = np.fromstring(values[1], dtype=float, count=1, sep=' ')

    # Read flux tubes
    #----------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "tubes" in alllines[i]:
            break 
        else: 
            i = i + 1
    
    # Read total amount of flux tubes and the size of face and cell
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=3, sep=' ')
    topomesh.ft.ntot = ntot[0]
    topomesh.ft.nface = ntot[1]
    topomesh.ft.ncell = ntot[2]
    
    # Initialize fields
    topomesh.ft.Initialize(topomesh.ft.ntot, topomesh.ft.nface, topomesh.ft.ncell)

    # Read data if available
    if topomesh.ft.ntot > 0:

        # Set ID 
        topomesh.ft.ID = np.arange(1, topomesh.ft.ntot+1, 1)

        # Skip header
        i = i + 2

        # Read facelist
        for j in np.arange(0, topomesh.ft.nface, 1):
            # Split the string
            values = alllines[i].split()

            # Add the values
            topomesh.ft.face[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')

            # Update counter
            i = i + 1

        # Read face pointer 
        i = i + 1
        for j in np.arange(0, topomesh.ft.ntot, 1):
            # Split the string
            values = alllines[i].split()

            # Add the values
            topomesh.ft.faceP[j, 0] = np.fromstring(values[0], dtype=int, count=1, sep=' ')-1
            topomesh.ft.faceP[j, 1] = np.fromstring(values[1], dtype=int, count=1, sep=' ')

            # Update counter
            i = i + 1

        # Read cell list
        i = i + 1
        for j in np.arange(0, topomesh.ft.ncell, 1):
            # Split the string
            values = alllines[i].split()

            # Add the values
            topomesh.ft.cell[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')

            # Update counter
            i = i + 1

        # Read cell pointer 
        i = i + 1
        for j in np.arange(0, topomesh.ft.ntot, 1):
            # Split the string
            values = alllines[i].split()

            # Add the values
            topomesh.ft.cellP[j, 0] = np.fromstring(values[0], dtype=int, count=1, sep=' ')-1
            topomesh.ft.cellP[j, 1] = np.fromstring(values[1], dtype=int, count=1, sep=' ')

            # Update counter
            i = i + 1

    # Return 
    return topomesh 

def ReadGGTMDataFile(filepath):
    # Description
    #------------
    # Read in the grid generation data of a topological mesh by reading 
    # the ggtmdata file (format specific here)
    ggtmdata = gt.GGTMData()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read basic face and cell data
    i = 0
    while i < len(alllines):
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1
    i = i + 1
    values = alllines[i].split()
    nftot = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nf = np.fromstring(values[1], dtype=int, count=1, sep =' ')
    nf = nf[0]

    i = 0
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1
    i = i + 1
    values = alllines[i].split()
    nctot = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nc = np.fromstring(values[1], dtype=int, count=1, sep =' ')
    nc = nc[0]

    # Initialize
    ggtmdata.Initialize(nftot[0], nctot[0])

    # Read in face data
    #------------------
    # Get the face header
    i = 0 # reset counter to start from beginning
    while i < len(alllines):
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read until ID header is found
    while i < len(alllines):
        if "ID, nc" in alllines[i]:
            break 
        else:
            i = i + 1

    # Skip header
    i = i + 1

    # Read face data 
    for j in np.arange(0, nf):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        ID = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        ID = ID[0]
        nv = np.fromstring(values[1], dtype=int, count=1, sep=' ')
        nv = nv[0]
        ggtmdata.face[ID-1].nv = nv
        ggtmdata.face[ID-1].ID = ID

    # Update position
    i = i + nf - 1

    # Get the face coordinates
    while i < len(alllines): 
        if "face coordinates" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read 
    counter = 0
    while (i < len(alllines)) and (counter < nf):
        # Read until we find 'face'
        if "face" in alllines[i]:
            # Get face ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1

            # Update position
            i = i + 1

            # Read nc coordinates
            ggtmdata.face[fID].Initialize(ggtmdata.face[fID].nv, ggtmdata.face[fID].ID)
            for k in np.arange(0, ggtmdata.face[fID].nv, 1):
                values = alllines[i+k].split()
                ggtmdata.face[fID].vert[k] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                ggtmdata.face[fID].x[k] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
                ggtmdata.face[fID].y[k] = np.fromstring(values[2], dtype=float, count=1, sep=' ')

            # Update position
            i = i + ggtmdata.face[fID].nv
        else: 
            i = i + 1 

    # Read cell data
    #---------------
    # Get the cell header
    i = 0 # reset counter to start from beginning
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read until ID header is found
    while i < len(alllines):
        if "ID, srf, erf, cell line size" in alllines[i]:
            break 
        else:
            i = i + 1

    # Skip header
    i = i + 1

    # Read cell data 
    for j in np.arange(0, nc):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        ID = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        ID = ID[0]
        ggtmdata.cell[ID-1].ID = ID
        srf = np.fromstring(values[1], dtype=int, count=1, sep=' ')
        erf = np.fromstring(values[2], dtype=int, count=1, sep=' ')
        nl = np.fromstring(values[3], dtype=int, count=1, sep=' ')
        ggtmdata.cell[ID-1].srf = srf[0]
        ggtmdata.cell[ID-1].erf = erf[0]
        ggtmdata.cell[ID-1].nl = nl[0]
        
    # Update file position
    i = i + nc-1

    # Get header position
    while i < len(alllines): 
        if "cell lines" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read line data
    counter = 0
    while (i < len(alllines)) and (counter < nc):
        # Read until we find 'cell'
        if "cell" in alllines[i]:
            # Get cell ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1

            # Read nl lines
            ggtmdata.cell[fID].Initialize(ggtmdata.cell[fID].nl-2, ggtmdata.cell[fID].ID)
            for k in np.arange(0, ggtmdata.cell[fID].nl+2, 1):
                # Update position
                i = i + 1

                # Read size of line
                values = alllines[i].split()

                # Initialize line
                nv = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                nv = nv[0]
                line = gt.GGTMLine()
                line.Initialize(nv)

                # Read 
                for cc in np.arange(0, nv, 1): 
                    values = alllines[i+cc+1].split()
                    line.vert[cc] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                    line.x[cc] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
                    line.y[cc] = np.fromstring(values[2], dtype=float, count=1, sep=' ')

                # Add line
                if k == 0:
                    # First line -> hfline
                    ggtmdata.cell[fID].hfline = line
                elif k == ggtmdata.cell[fID].nl+1:
                    # Last line -> lfline
                    ggtmdata.cell[fID].lfline = line
                else:
                    # intermediate line
                    ggtmdata.cell[fID].lines[k-1] = line

                # Update position
                i = i + nv

        else: 
            i = i + 1 

    # Get header position
    while i < len(alllines): 
        if "cell tubes" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read tube data
    counter = 0
    while (i < len(alllines)) and (counter < nc):
        # Read until we find 'cell'
        if "cell" in alllines[i]:
            # Get cell ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1
            i = i + 1

            # Read nl+1 tube hf and lf lines
            for k in np.arange(0, ggtmdata.cell[fID].nl+1, 1):
                # Update position
                i = i + 1

                # Read size of line
                values = alllines[i].split()

                # Initialize hfline
                nvhf = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                nvhf = nvhf[0]
                xhf = np.zeros(nvhf, dtype=float)
                yhf = np.zeros(nvhf, dtype=float)
                vhf = np.zeros(nvhf, dtype=int)

                # Read 
                for cc in np.arange(0, nvhf, 1): 
                    values = alllines[i+cc+1].split()
                    vhf[cc] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                    xhf[cc] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
                    yhf[cc] = np.fromstring(values[2], dtype=float, count=1, sep=' ')

                # Update position
                i = i + nvhf + 2
                
                # Initialize lfline
                values = alllines[i].split()
                nvlf = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                nvlf = nvlf[0]
                xlf = np.zeros(nvlf, dtype=float)
                ylf = np.zeros(nvlf, dtype=float)
                vlf = np.zeros(nvlf, dtype=int)

                # Read 
                for cc in np.arange(0, nvlf, 1): 
                    values = alllines[i+cc+1].split()
                    vlf[cc] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                    xlf[cc] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
                    ylf[cc] = np.fromstring(values[2], dtype=float, count=1, sep=' ')

                # Add tube
                ggtmdata.cell[fID].AddTube(k, xhf, yhf, vhf, xlf, ylf, vlf)

                # Update position
                i = i + nvlf + 1

        else: 
            i = i + 1 

    # Return
    return ggtmdata

def ReadGGGridDataFile(filepath):
    # Description
    #------------
    # This routine reads in the grid data from an intermediate grid 
    # generator grid object. This is not the same as the full grid
    # that is used in the grid deformation module.

    # Initialize
    grid = gt.GGGrid()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read vertex data
    #-----------------
    i = 0
    while i < len(alllines):
        if "vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of vertices
    values = alllines[i].split()
    nv = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nv = nv[0]

    # Initialize vertex structure
    grid.vert.Initialize(nv)

    # Read vertex data
    while i < len(alllines):
        if "ID, x, y, fieldlineID" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nv):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        x = np.fromstring(values[1], dtype=float, count=1, sep =' ')
        y = np.fromstring(values[2], dtype=float, count=1, sep =' ')
        fID = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        grid.vert.ID[ID-1] = ID
        grid.vert.x[ID-1] = x[0]
        grid.vert.y[ID-1] = y[0]
        grid.vert.fieldlineID[ID-1] = fID[0] 

    # Read in face data
    #------------------
    i = 0
    while i < len(alllines):
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of faces
    values = alllines[i].split()
    nf = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nf = nf[0]

    # Initialize vertex structure
    grid.face.Initialize(nf)

    # Read vertex data
    while i < len(alllines):
        if "ID, v1, v2, label, region" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nf):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        v1 = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        v2 = np.fromstring(values[2], dtype=int, count=1, sep =' ')
        label = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        region = np.fromstring(values[4], dtype=int, count=1, sep =' ')
        grid.face.ID[ID-1] = ID
        grid.face.v1[ID-1] = v1[0]
        grid.face.v2[ID-1] = v2[0]
        grid.face.label[ID-1] = label[0] 
        grid.face.region[ID-1] = region[0]

    # Read in cell data
    #------------------
    i = 0
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of cells and cell vertices
    values = alllines[i].split()
    nc = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    ncv = np.fromstring(values[1], dtype=int, count=1, sep =' ')
    nc = nc[0]
    ncv = ncv[0]

    # Initialize vertex structure
    grid.cell.Initialize(nc, ncv)

    # Read cell data
    while i < len(alllines):
        if "ID, vp1, vp2, region" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nc):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        v1 = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        v2 = np.fromstring(values[2], dtype=int, count=1, sep =' ')
        region = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        grid.cell.ID[ID-1] = ID
        grid.cell.vp1[ID-1] = v1[0]-1 # Need to account for 0-based indexing
        grid.cell.vp2[ID-1] = v2[0]
        grid.cell.region[ID-1] = region[0]

    # Read cell vertices
    while i < len(alllines):
        if "cell vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read vertices
    for j in np.arange(0, ncv):
        values = alllines[i+j].split()
        tv = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        tv = tv[0]
        grid.cell.vert[j] = tv

    # Return
    return grid

def ReadGAIntegerArrayFile(filepath):
    # Description
    # -----------
    # This routine read a file with an array as writting by WriteArray
    # in gamod_types.F90.

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read array data
    #----------------
    i = 0
    while i < len(alllines):
        if "Elements" in alllines[i]:
            break
        else:
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of elements
    values = alllines[i].split()
    n_el = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    n_el = n_el[0]

    # Read vertex data
    while i < len(alllines):
        if "ID val(ID)" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    ar = np.zeros(n_el, dtype=int)
    for j in np.arange(0, n_el):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        el = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        ar[ID-1] = el[0]

    return ar

def ReadGARealArrayFile(filepath):
    # Description
    # -----------
    # This routine read a file with an array as writting by WriteArray
    # in gamod_types.F90.

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read array data
    #----------------
    i = 0
    while i < len(alllines):
        if "Elements" in alllines[i]:
            break
        else:
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of elements
    values = alllines[i].split()
    n_el = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    n_el = n_el[0]

    # Read vertex data
    while i < len(alllines):
        if "ID val(ID)" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    ar = np.zeros(n_el, dtype=float)
    for j in np.arange(0, n_el):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        el = np.fromstring(values[1], dtype=float, count=1, sep =' ')
        ar[ID-1] = el[0]

    return ar

def ReadGAGridDataFile(filepath):
    # Description
    #------------
    # This routine reads in the grid data from an 
    # adaptation grid object. This is not the same as the full grid
    # that is used in the grid deformation module.

    # Initialize
    grid = gt.GAGrid()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read vertex data
    #-----------------
    i = 0
    while i < len(alllines):
        if "vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of vertices
    values = alllines[i].split()
    nv = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nv = nv[0]

    # Initialize vertex structure
    grid.vert.Initialize(nv)

    # Read vertex data
    while i < len(alllines):
        if "ID, x, y" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nv):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        x = np.fromstring(values[1], dtype=float, count=1, sep =' ')
        y = np.fromstring(values[2], dtype=float, count=1, sep =' ')
        grid.vert.ID[ID-1] = ID
        grid.vert.x[ID-1] = x[0]
        grid.vert.y[ID-1] = y[0] 

    # Read in face data
    #------------------
    i = 0
    while i < len(alllines):
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of faces
    values = alllines[i].split()
    nf = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nf = nf[0]

    # Initialize vertex structure
    grid.face.Initialize(nf)

    # Read vertex data
    while i < len(alllines):
        if "ID, v1, v2, label, aligned" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nf):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        v1 = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        v2 = np.fromstring(values[2], dtype=int, count=1, sep =' ')
        label = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        aligned = np.fromstring(values[4], dtype=int, count=1, sep =' ')
        grid.face.ID[ID-1] = ID
        grid.face.v1[ID-1] = v1[0]
        grid.face.v2[ID-1] = v2[0]
        grid.face.label[ID-1] = label[0] 
        grid.face.aligned[ID-1] = aligned[0]

    # Read in cell data
    #------------------
    i = 0
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of cells and cell vertices
    values = alllines[i].split()
    nc = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    ncv = np.fromstring(values[1], dtype=int, count=1, sep =' ')
    nc = nc[0]
    ncv = ncv[0]

    # Initialize vertex structure
    grid.cell.Initialize(nc, ncv)

    # Read cell data
    while i < len(alllines):
        if "ID, vp1, vp2, region, x, y" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nc):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        v1 = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        v2 = np.fromstring(values[2], dtype=int, count=1, sep =' ')
        region = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        x = np.fromstring(values[4], dtype=float, count=1, sep =' ')
        y = np.fromstring(values[5], dtype=float, count=1, sep =' ')
        grid.cell.ID[ID-1] = ID
        grid.cell.vp1[ID-1] = v1[0]-1 # Need to account for 0-based indexing
        grid.cell.vp2[ID-1] = v2[0]
        grid.cell.region[ID-1] = region[0]
        grid.cell.x[ID-1] = x[0]
        grid.cell.y[ID-1] = y[0]

    # Read cell vertices
    while i < len(alllines):
        if "cell vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read vertices
    for j in np.arange(0, ncv):
        values = alllines[i+j].split()
        tv = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        tv = tv[0]
        grid.cell.vert[j] = tv

    # Return
    return grid

def ReadTraduitOutB2us(filepath):
    # Description
    #------------
    # This routine reads in the grid data in unstructured traduit file
    # into a 'Grid' object. This does not contain any guard cells yet

    # Initialize
    grid = gt.Grid()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read version to determine what data to read in 
    temp = alllines[0].split(); temp = temp[0]
    version = temp[7:17]
    topodataversion = '03.002.001'

    if version >= topodataversion:
        hasTopologicalData = True
    else: 
        hasTopologicalData = False 
    

    # Read dimensions
    #----------------
    # Start at top of file
    i = 0
    while i < len(alllines):
        if "nCi,nFc,nVx,nCg,nFs,nFt" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read dimensions
    values = alllines[i].split()
    nc = np.fromstring(values[0], dtype=int, count=1, sep =' '); nc = nc[0]
    nf = np.fromstring(values[1], dtype=int, count=1, sep =' '); nf = nf[0]
    nv = np.fromstring(values[2], dtype=int, count=1, sep =' '); nv = nv[0]
    nfs = np.fromstring(values[4], dtype=int, count=1, sep =' '); nfs = nfs[0]
    nft = np.fromstring(values[5], dtype=int, count=1, sep =' '); nft = nft[0]

    # Read secondary dimensions (continue from previous line)
    while i < len(alllines):
        if "nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read dimensions
    values = alllines[i].split()
    ncv = np.fromstring(values[0], dtype=int, count=1, sep =' '); ncv = ncv[0]
    ncf = np.fromstring(values[1], dtype=int, count=1, sep =' '); ncf = ncf[0]

    # Initialize the grid (except flux surfaces etc - later on)
    grid.vert.Initialize(nv)
    grid.face.Initialize(nf)
    grid.cell.Initialize(nc, 0, ncv, ncf)
    grid.ft.Initialize(nft)
    grid.fs.Initialize(nfs)

    # Initialize vertex structure
    grid.vert.Initialize(nv)

    # Read topological data
    #----------------------
    if hasTopologicalData:
        while i < len(alllines):
            if "topoflag" in alllines[i]:
                break 
            else: 
                i = i + 1

        # Read in topoflag
        i = i + 1
        values = alllines[i].split()
        topoflag = np.fromstring(values[0], dtype=int, count=1, sep =' '); topoflag = topoflag[0]
        i = i + 1

        # Read in numbers
        i = i + 1
        values = alllines[i].split()
        nX = np.fromstring(values[0], dtype=int, count=1, sep =' '); nX = nX[0]
        nO = np.fromstring(values[1], dtype=int, count=1, sep =' '); nO = nO[0]
        nS = np.fromstring(values[2], dtype=int, count=1, sep =' '); nS = nS[0]
        nT = np.fromstring(values[3], dtype=int, count=1, sep =' '); nT = nT[0]
        nDiv = np.fromstring(values[4], dtype=int, count=1, sep =' '); nDiv = nDiv[0]
        nDivFc = np.fromstring(values[5], dtype=int, count=1, sep =' '); nDivFc = nDivFc[0]
        i = i + 1

        # Initialize
        grid.topodata.Initialize(nX, nO, nS, nT, nDiv, nDivFc, topoflag)

        # Read x-point data
        i = i + 1
        for j in np.arange(0, nX, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]
            isprimaryxp = np.fromstring(values[1], dtype=int, count=1, sep =' '); isprimaryxp = isprimaryxp[0]

            grid.topodata.XpointID[j] = xind 
            grid.topodata.isprimaryxp[j] = isprimaryxp 
        i = i + nX 

        # Read O-point data
        i = i + 1
        for j in np.arange(0, nO, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]

            grid.topodata.OpointID[j] = xind 
        i = i + nO

        # Read S-point data
        i = i + 1
        for j in np.arange(0, nS, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]
            spointxpID = np.fromstring(values[1], dtype=int, count=1, sep =' '); spointxpID = spointxpID[0]
            spointdivID = np.fromstring(values[2], dtype=int, count=1, sep =' '); spointdivID = spointdivID[0]

            grid.topodata.SpointID[j] = xind 
            grid.topodata.spointdivID[j] = spointdivID 
            grid.topodata.spointxpID[j] = spointxpID 
        i = i + nS

        # Read T-point data
        i = i + 1
        for j in np.arange(0, nT, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]
            tpointdivID = np.fromstring(values[1], dtype=int, count=1, sep =' '); tpointdivID = tpointdivID[0]

            grid.topodata.TpointID[j] = xind 
            grid.topodata.tpointdivID[j] = tpointdivID 
        i = i + nT

        # Read divertor face pointer
        i = i + 1
        for j in np.arange(0, nDiv, 1):
            values = alllines[i+j].split()
            divFcP1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); divFcP1 = divFcP1[0]
            divFcP2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); divFcP2 = divFcP2[0]

            grid.topodata.divFcP1[j] = divFcP1-1 # account for zero-based indexing 
            grid.topodata.divFcP2[j] = divFcP2 
        i = i + nDiv 

        # Read divertor face list
        i = i + 1
        k = 0
        while k < nDivFc:
            values = alllines[i].split()
            for j in np.arange(0, len(values)):
                ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
                grid.topodata.divFc[k] = ID
                k = k + 1
            i = i + 1
    else:
        # Simply initialize to zero
        grid.topodata.Initialize(0, 0, 0, 0, 0, 0, 0)
        

    # Read vertex data
    #-----------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: Vx vxX vxY vxPsi vxBx vxBy vxFfbz" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nv):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        x = np.fromstring(values[1], dtype=float, count=1, sep =' '); x = x[0]
        y = np.fromstring(values[2], dtype=float, count=1, sep =' '); y = y[0]
        psi = np.fromstring(values[3], dtype=float, count=1, sep =' '); psi = psi[0]
        bx = np.fromstring(values[4], dtype=float, count=1, sep =' '); bx = bx[0]
        by = np.fromstring(values[5], dtype=float, count=1, sep =' '); by = by[0]
        ffbz = np.fromstring(values[6], dtype=int, count=1, sep =' '); ffbz = ffbz[0]
        grid.vert.ID[ID-1] = ID
        grid.vert.x[ID-1] = x
        grid.vert.y[ID-1] = y
        grid.vert.bx[ID-1] = bx
        grid.vert.by[ID-1] = by
        grid.vert.psi[ID-1] = psi
        grid.vert.ffbz[ID-1] = ffbz

    # Read cell data
    #---------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: cv cvVxP(:,1) cvVxP(:,2) cvX cvY psi bp bt cflags(:) cvReg cvFt" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nc):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        vp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); vp1 = vp1[0]
        vp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); vp2 = vp2[0]
        x = np.fromstring(values[3], dtype=float, count=1, sep =' '); x = x[0]
        y = np.fromstring(values[4], dtype=float, count=1, sep =' '); y = y[0]
        psi = np.fromstring(values[5], dtype=float, count=1, sep =' '); psi = psi[0]
        bp = np.fromstring(values[6], dtype=float, count=1, sep =' '); bp = bp[0]
        bt = np.fromstring(values[7], dtype=float, count=1, sep =' '); bt = bt[0]
        cflags = np.fromstring(values[8], dtype=int, count=1, sep =' '); cflags = cflags[0]
        cvreg = np.fromstring(values[9], dtype=int, count=1, sep =' '); cvreg = cvreg[0]
        cvft = np.fromstring(values[10], dtype=int, count=1, sep =' '); cvft = cvft[0]
        
        grid.cell.ID[ID-1] = ID
        ID = ID-1
        grid.cell.vp1[ID] = vp1-1 # Need to account for 0-based indexing
        grid.cell.vp2[ID] = vp2
        grid.cell.x[ID] = x
        grid.cell.y[ID] = y
        grid.cell.psi[ID] = psi
        grid.cell.bp[ID] = bp
        grid.cell.bt[ID] = bt
        grid.cell.cflags[ID] = cflags
        grid.cell.region[ID] = cvreg
        grid.cell.ft[ID] = cvft

    # Update i 
    i = i + nc

    # Vertex and face pointers are the same
    grid.cell.fp1 = grid.cell.vp1
    grid.cell.fp2 = grid.cell.vp2

    # Skip next line for reading cell vertices
    i = i + 1

    # Read cell vertices
    k = 0
    while k < ncv:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.cell.vert[k] = ID
            k = k + 1
        i = i + 1

    # Skip header
    i = i + 1

    # Read cell faces
    k = 0
    while k < ncf:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.cell.face[k] = ID
            k = k + 1
        i = i + 1

    # Read face data
    #---------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: fc fcVx(:,1) fcVx(:,2) fcLbl fcReg fcAligned" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Start reading
    i = i + 1
    for j in np.arange(0, nf):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        vp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); vp1 = vp1[0]
        vp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); vp2 = vp2[0]
        label = np.fromstring(values[3], dtype=int, count=1, sep =' '); label = label[0]
        region = np.fromstring(values[4], dtype=int, count=1, sep =' '); region = region[0]
        aligned = np.fromstring(values[5], dtype=int, count=1, sep =' '); aligned = aligned[0]

        grid.face.ID[ID-1] = ID
        ID = ID-1
        grid.face.v1[ID] = vp1
        grid.face.v2[ID] = vp2
        grid.face.label[ID] = label
        grid.face.region[ID] = region
        grid.face.aligned[ID] = aligned

    # Update i 
    i = i + nf

    # Flux tube
    #----------
    # Continue from last line
    while i < len(alllines):
        if "*cf: ft ftCvP(:,1) ftCvP(:,2) ftFcP(:,1) ftFcP(:,2) ftReg" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nft):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        cp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); cp1 = cp1[0]
        cp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); cp2 = cp2[0]
        fp1 = np.fromstring(values[3], dtype=int, count=1, sep =' '); fp1 = fp1[0]
        fp2 = np.fromstring(values[4], dtype=int, count=1, sep =' '); fp2 = fp2[0]
        region = np.fromstring(values[5], dtype=int, count=1, sep =' '); region = region[0]

        grid.ft.ID[ID-1] = ID
        ID = ID-1
        grid.ft.cp1[ID] = cp1-1 # Need to account for 0-based indexing
        grid.ft.cp2[ID] = cp2
        grid.ft.fp1[ID] = fp1-1 # Need to account for 0-based indexing
        grid.ft.fp2[ID] = fp2
        grid.ft.region[ID] = region

    # Update
    i = i + nft 

    # Compute and initialize face and cell data
    nftf = grid.ft.fp1[nft-1]+grid.ft.fp2[nft-1]
    nftc = grid.ft.cp1[nft-1]+grid.ft.cp2[nft-1]
    grid.ft.InitializeCellData(nftc)
    grid.ft.InitializeFaceData(nftf)

    # Skip header and read
    i = i + 1
    k = 0
    while k < nftc:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.ft.cell[k] = ID
            k = k + 1
        i = i + 1
    

    # Update and read
    i = i + 1
    k = 0
    while k < nftf:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.ft.face[k] = ID
            k = k + 1
        i = i + 1
    

    # Flux surfaces
    #--------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: fs fsFcP(:,1) fsFcP(:,2) fsPsi" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nfs):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        fp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); fp1 = fp1[0]
        fp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); fp2 = fp2[0]
        psi = np.fromstring(values[3], dtype=float, count=1, sep =' '); psi = psi[0]

        grid.fs.ID[ID-1] = ID
        ID = ID-1
        grid.fs.fp1[ID] = fp1-1 # Need to account for 0-based indexing
        grid.fs.fp2[ID] = fp2
        grid.fs.psi[ID] = psi

    # Update
    i = i + nfs 

    # Compute and initialize face data
    nfsf = grid.fs.fp1[nfs-1]+grid.fs.fp2[nfs-1]
    grid.fs.InitializeFaceData(nfsf)

    # Skip header and read
    i = i + 1
    k = 0
    while k < nfsf:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.fs.face[k] = ID
            k = k + 1
        i = i + 1

    # Construct vertex field line ID
    for j in np.arange(0, nfs):
        tf = grid.fs.face[grid.fs.fp1[j]:grid.fs.fp1[j]+grid.fs.fp2[j]]
        tv1 = grid.face.v1[tf-1]
        tv2 = grid.face.v2[tf-1]
        for k in tv1:
            grid.vert.fieldlineID[k-1] = j+1
        for k in tv2:
            grid.vert.fieldlineID[k-1] = j+1

    # Compute grid interconnections
    grid.ComputeInterconnections()

    # Compute face coordinates
    grid.ComputeFaceCoordinates()

    # Compute cell coordinates
    grid.ComputeCellCoordinates()
    
    # Return
    return grid

def ReadGridFromB2fgmtryus(filepath):
    # Description
    #------------
    # This routine reads in the grid data in unstructured b2fgmtry file
    # format. Besides the usual grid data (without guard cells), also
    # grid data of guard cells, and coordinate data of faces and cells
    # are loaded. 
    # Note: not all data available in b2fgmtry is loaded in, only the 
    # data required for post-processing purposes

    # Initialize
    grid = gt.Grid()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read version to determine what data to read in 
    temp = alllines[0].split(); temp = temp[0]
    version = temp[7:17]
    topodataversion = '03.002.001'

    if version >= topodataversion:
        hasTopologicalData = True
    else: 
        hasTopologicalData = False 
    

    # Read dimensions
    #----------------
    # Start at top of file
    i = 0
    while i < len(alllines):
        if "nCi,nCg,nCv,nFc,nVx,nFs,nFt" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read dimensions
    values = alllines[i].split()
    nci = np.fromstring(values[0], dtype=int, count=1, sep =' '); nci = nci[0]
    ncg = np.fromstring(values[1], dtype=int, count=1, sep =' '); ncg = ncg[0]
    nc = np.fromstring(values[2], dtype=int, count=1, sep =' '); nc = nc[0]
    nf = np.fromstring(values[3], dtype=int, count=1, sep =' '); nf = nf[0]
    nv = np.fromstring(values[4], dtype=int, count=1, sep =' '); nv = nv[0]
    nfs = np.fromstring(values[5], dtype=int, count=1, sep =' '); nfs = nfs[0]
    nft = np.fromstring(values[6], dtype=int, count=1, sep =' '); nft = nft[0]

    # Read secondary dimensions (continue from previous line)
    while i < len(alllines):
        if "nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmx" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read dimensions
    values = alllines[i].split()
    ncv = np.fromstring(values[0], dtype=int, count=1, sep =' '); ncv = ncv[0]
    ncf = np.fromstring(values[1], dtype=int, count=1, sep =' '); ncf = ncf[0]

    # Initialize the grid (except flux surfaces etc - later on)
    grid.vert.Initialize(nv)
    grid.face.Initialize(nf)
    grid.cell.Initialize(nci, ncg, ncv, ncf)
    grid.ft.Initialize(nft)
    grid.fs.Initialize(nfs)

    # Initialize vertex structure
    grid.vert.Initialize(nv)

    # Read topological data
    #----------------------
    # First only header, other information is further downstream
    if hasTopologicalData:
        while i < len(alllines):
            if "topoflag" in alllines[i]:
                break 
            else: 
                i = i + 1

        # Read in topoflag
        i = i + 1
        values = alllines[i].split()
        topoflag = np.fromstring(values[0], dtype=int, count=1, sep =' '); topoflag = topoflag[0]
        i = i + 1

        # Read in numbers
        i = i + 1
        values = alllines[i].split()
        nX = np.fromstring(values[0], dtype=int, count=1, sep =' '); nX = nX[0]
        nO = np.fromstring(values[1], dtype=int, count=1, sep =' '); nO = nO[0]
        nS = np.fromstring(values[2], dtype=int, count=1, sep =' '); nS = nS[0]
        nT = np.fromstring(values[3], dtype=int, count=1, sep =' '); nT = nT[0]
        nDiv = np.fromstring(values[4], dtype=int, count=1, sep =' '); nDiv = nDiv[0]
        nDivFc = np.fromstring(values[5], dtype=int, count=1, sep =' '); nDivFc = nDivFc[0]
        i = i + 1   
        
    else:
        # Simply initialize to zero
        nX = 0 
        nO = 0
        nS = 0
        nT = 0
        nDiv = 0
        nDivFc = 0
        topoflag = 0

    # Initialize
    grid.topodata.Initialize(nX, nO, nS, nT, nDiv, nDivFc, topoflag)        

    # Read interconnection data
    #--------------------------
    # Read
    cvFcP, i = ReadSOLPSField('cvFcP', alllines, i, nc*2, 'int', np.array([nc, 2]))
    cvFc,  i = ReadSOLPSField('cvFc', alllines, i, ncf, 'int', np.array([ncf]))
    fcCv,  i = ReadSOLPSField('fcCv', alllines, i, nf*2, 'int', np.array([nf, 2]))
    fcVx,  i = ReadSOLPSField('fcVx', alllines, i, nf*2, 'int', np.array([nf, 2]))
    cvVxP, i = ReadSOLPSField('cvVxP', alllines, i, nc*2, 'int', np.array([nc, 2]))
    cvVx,  i = ReadSOLPSField('cvVx', alllines, i, ncv, 'int', np.array([ncv]))
    # vxFcP: not read
    # vxFc: not read
    # vxCvP: not read
    # vxCv: not read
    ftCvP, i = ReadSOLPSField('ftCvP', alllines, i, nft*2, 'int', np.array([nft, 2]))
    ftCv,  i = ReadSOLPSField('ftCv', alllines, i, nc, 'int', np.array([nc]))
    ftFcP, i = ReadSOLPSField('ftFcP', alllines, i, nft*2, 'int', np.array([nft, 2]))
    ftFc,  i = ReadSOLPSField('ftFc', alllines, i, nf, 'int', np.array([nf])) # in this case equal to amount of faces...
    cvFt,  i = ReadSOLPSField('cvFt', alllines, i, nc, 'int', np.array([nc]))
    fsFcP, i = ReadSOLPSField('fsFcP', alllines, i, nfs*2, 'int', np.array([nfs, 2]))
    fsFc,  i = ReadSOLPSField('fsFc', alllines, i, nf, 'int', np.array([nf]))
    fcReg, i = ReadSOLPSField('fcReg', alllines, i, nf, 'int', np.array([nf]))
    cvReg, i = ReadSOLPSField('cvReg', alllines, i, nc, 'int', np.array([nc]))
    ftReg, i = ReadSOLPSField('ftReg', alllines, i, nft, 'int', np.array([nft]))
    # intcellP, intcellR, imap*, icorn: not read
    fcLbl, i = ReadSOLPSField('fcLbl', alllines, i, nf, 'int', np.array([nf]))
    cvLbl, i = ReadSOLPSField('cvLbl', alllines, i, nc, 'int', np.array([nc]))
    # ftLbl: not read

    # Add data
    grid.cell.fp1 = cvFcP[:, 0]-1
    grid.cell.fp2 = cvFcP[:, 1]
    grid.cell.vp1 = cvVxP[:, 0]-1
    grid.cell.vp2 = cvVxP[:, 1]
    grid.cell.face  = cvFc 
    grid.cell.vert  = cvVx
    grid.cell.ft = cvFt 
    grid.cell.region = cvReg 
    grid.cell.cflags = cvLbl

    grid.face.nb1 = fcCv[:, 0]
    grid.face.nb2 = fcCv[:, 1]
    grid.face.v1 = fcVx[:, 0]
    grid.face.v2 = fcVx[:, 1]
    grid.face.region = fcReg
    grid.face.label = fcLbl

    grid.ft.cp1 = ftCvP[:, 0]-1
    grid.ft.cp2 = ftCvP[:, 1]
    grid.ft.cell = ftCv 
    grid.ft.fp1 = ftFcP[:, 0]-1
    grid.ft.fp2 = ftFcP[:, 1]
    grid.ft.face = ftFc 
    grid.ft.region = ftReg 

    grid.fs.fp1 = fsFcP[:, 0]-1
    grid.fs.fp2 = fsFcP[:, 1]
    grid.fs.face = fsFc

    
    # Read additional topological data
    #---------------------------------
    if hasTopologicalData:
        # Read x-point data
        while i < len(alllines):
            if "*cf:    Vx isPrimary fsID" in alllines[i]:
                break 
            else: 
                i = i + 1
        i = i + 1
        for j in np.arange(0, nX, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]
            isprimaryxp = np.fromstring(values[1], dtype=int, count=1, sep =' '); isprimaryxp = isprimaryxp[0]

            grid.topodata.XpointID[j] = xind 
            grid.topodata.isprimaryxp[j] = isprimaryxp 
        i = i + nX 

        # Read S-point data
        i = i + 1
        for j in np.arange(0, nS, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]
            spointxpID = np.fromstring(values[1], dtype=int, count=1, sep =' '); spointxpID = spointxpID[0]
            spointdivID = np.fromstring(values[2], dtype=int, count=1, sep =' '); spointdivID = spointdivID[0]

            grid.topodata.SpointID[j] = xind 
            grid.topodata.spointdivID[j] = spointdivID 
            grid.topodata.spointxpID[j] = spointxpID 
        i = i + nS

        # Read O-point data
        i = i + 1
        for j in np.arange(0, nO, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]

            grid.topodata.OpointID[j] = xind 
        i = i + nO

        

        # Read T-point data - still some things going wrong at b2fgmtry side...
        #i = i + 1
        for j in np.arange(0, nT, 1):
            values = alllines[i+j].split()
            xind = np.fromstring(values[0], dtype=int, count=1, sep =' '); xind = xind[0]
            tpointdivID = np.fromstring(values[1], dtype=int, count=1, sep =' '); tpointdivID = tpointdivID[0]

            grid.topodata.TpointID[j] = xind 
            grid.topodata.tpointdivID[j] = tpointdivID 
        i = i + nT

        # Read divertor face pointer
        i = i + 1
        for j in np.arange(0, nDiv, 1):
            values = alllines[i+j].split()
            divFcP1 = np.fromstring(values[0], dtype=int, count=1, sep =' '); divFcP1 = divFcP1[0]
            divFcP2 = np.fromstring(values[1], dtype=int, count=1, sep =' '); divFcP2 = divFcP2[0]

            grid.topodata.divFcP1[j] = divFcP1-1 # account for zero-based indexing 
            grid.topodata.divFcP2[j] = divFcP2 
        i = i + nDiv 

        # Read divertor face list
        i = i + 1
        k = 0
        while k < nDivFc:
            values = alllines[i].split()
            for j in np.arange(0, len(values)):
                ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
                grid.topodata.divFc[k] = ID
                k = k + 1
            i = i + 1

    # Read metrics
    #-------------
    cvBb, i = ReadSOLPSField('cvBb', alllines, i, 4*nc, 'float', np.array([nc, 4]))
    cvEb, i = ReadSOLPSField('cvEb', alllines, i, 3*nc, 'float', np.array([nc, 3]))
    cvX,  i = ReadSOLPSField('cvX', alllines, i, nc, 'float', np.array([nc]))
    cvY,  i = ReadSOLPSField('cvY', alllines, i, nc, 'float', np.array([nc]))
    cvSz,  i = ReadSOLPSField('cvSz', alllines, i, nc, 'float', np.array([nc]))
    cvHz,  i = ReadSOLPSField('cvHz', alllines, i, nc, 'float', np.array([nc]))
    cvHx,  i = ReadSOLPSField('cvHx', alllines, i, nc, 'float', np.array([nc]))
    cvQgam,  i = ReadSOLPSField('cvQgam', alllines, i, 2*nc, 'float', np.array([nc, 2]))
    cvVol,  i = ReadSOLPSField('cvVol', alllines, i, nc, 'float', np.array([nc]))

    fcBb, i = ReadSOLPSField('fcBb', alllines, i, 4*nf, 'float', np.array([nf, 4]))
    fcS, i = ReadSOLPSField('fcS', alllines, i, nf, 'float', np.array([nf]))
    fcHc, i = ReadSOLPSField('fcHc', alllines, i, 2*nf, 'float', np.array([nf, 2]))
    fcHt, i = ReadSOLPSField('fcHt', alllines, i, nf, 'float', np.array([nf]))
    fcQgam,  i = ReadSOLPSField('fcQgam', alllines, i, 2*nf, 'float', np.array([nf, 2]))
    fcQalf,  i = ReadSOLPSField('fcQalf', alllines, i, 2*nf, 'float', np.array([nf, 2]))
    fcQbet,  i = ReadSOLPSField('fcQbet', alllines, i, 2*nf, 'float', np.array([nf, 2]))
    fcPbs, i = ReadSOLPSField('fcPbs', alllines, i, nf, 'float', np.array([nf]))

    vxBb, i = ReadSOLPSField('vxBb', alllines, i, 4*nv, 'float', np.array([nv, 4]))
    vxX, i = ReadSOLPSField('vxX', alllines, i, nv, 'float', np.array([nv]))
    vxY, i = ReadSOLPSField('vxY', alllines, i, nv, 'float', np.array([nv]))
    vxFfbz, i = ReadSOLPSField('vxFfbz', alllines, i, nv, 'float', np.array([nv]))
    vxfsPsi, i = ReadSOLPSField('vxfsPsi', alllines, i, nv, 'float', np.array([nv]))

    # cvConn, fsPsi: skipped

    grid.cell.x = cvX 
    grid.cell.y = cvY
    grid.cell.bb = cvBb 
    grid.cell.eb = cvEb
    grid.cell.sz = cvSz
    grid.cell.hz = cvHz
    grid.cell.hx = cvHx    
    grid.cell.qgam = cvQgam
    grid.cell.vol = cvVol 

    grid.face.bb = fcBb
    grid.face.s = fcS 
    grid.face.hc = fcHc
    grid.face.ht = fcHt
    grid.face.qgam = fcQgam
    grid.face.qbet = fcQbet
    grid.face.qalf = fcQalf
    grid.face.pbs = fcPbs

    grid.vert.bb = vxBb 
    grid.vert.x = vxX
    grid.vert.y = vxY 
    grid.vert.ffbz = vxFfbz
    grid.vert.fspsi = vxfsPsi 


    # Compute grid interconnections
    grid.ComputeInterconnections()

    # Compute face coordinates
    grid.ComputeFaceCoordinates()
    
    # Return
    return grid

def ReadStructureFile(filepath):
    # Description
    #------------
    # This function reads a structure.dat file and returns the structure

    # Open file
    thisfile = open(filepath)

    # Read all lines
    alllines = thisfile.readlines()

    # Read total number of structures
    values = alllines[0].split()
    ns = np.fromstring(values[0], dtype=int, count=1, sep =' '); ns = ns[0]

    # Initialize
    structure = [gt.Structure() for i in np.arange(0, ns, 1)]

    # Read
    lind = 2 # skip header
    for i in np.arange(0, ns, 1):
        # Read structure header to get the structure ID
        values = alllines[lind].split()
        sID = np.fromstring(values[1], dtype=int, count=1, sep =' '); sID = sID[0]

        # Update line index
        lind = lind + 1

        # Read number of structure coordinates
        values = alllines[lind].split()
        nc = np.fromstring(values[0], dtype=int, count=1, sep =' '); nc = nc[0]

        # Read coordinates
        lind = lind + 1
        tx = np.zeros(abs(nc))
        ty = np.zeros(abs(nc))
        for j in np.arange(0, abs(nc), 1):
            # Read
            values = alllines[lind].split()
            xp = np.fromstring(values[0], dtype=float, count=1, sep =' '); xp = xp[0]
            yp = np.fromstring(values[1], dtype=float, count=1, sep =' '); yp = yp[0]
            tx[j] = xp
            ty[j] = yp 

            # Update counter
            lind = lind + 1
        
        # Make structure
        structure[i].Initialize(abs(nc), tx, ty, sID)

    # Return
    return structure

def WriteStructureFile(dirpath, structure):
    # Description
    #------------
    # write a set of structures with specified, x, y coordinates and
    # number of points n (and indication if open or closed by sign of
    # n, negative: open, positive: closed) to a file specified by 
    # dirpath/structure.dat (only dirpath should be given)

    thisfile = open(dirpath + r"/structure.dat", "w")

    # Write total number of structures
    thisfile.write("             " + str(len(structure)) + "\n")

    # Write header
    thisfile.write("$structures\n")

    # Loop over all structures
    k = 1
    for i in structure:
        # Write structure header
        thisfile.write("Structure    " + str(k) + "\n")

        # Write number
        thisfile.write("             " + str(i.n) + "\n")
        
        # Loop over coordinates
        for j in np.arange(0, abs(i.n), 1):
            thisfile.write(str(i.x[j]) + "    " + str(i.y[j]) + " \n")

        # Update counter
        k = k + 1

    # Write final end
    thisfile.write("$end\n")

    # Close the file
    thisfile.close()

def ReadTriaGrid(filepath):
    # Description
    #------------
    # This routine reads in the triangulation grid data writen by 
    # the subroutine VisualizeTriangulation

    tria = gt.TriaGrid()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines() 

    # Read vertex data
    #-----------------
    i = 0
    while i < len(alllines):
        if "vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of vertices
    values = alllines[i].split()
    nv = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nv = nv[0]

    # Initialize vertex structure
    tria.vert.Initialize(nv)

    # Read vertex data
    while i < len(alllines):
        if "ID, x, y" in alllines[i]:
            break 
        else: 
            i = i + 1  

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nv):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        x = np.fromstring(values[1], dtype=float, count=1, sep =' ')
        y = np.fromstring(values[2], dtype=float, count=1, sep =' ')
        tria.vert.ID[ID-1] = ID
        tria.vert.x[ID-1] = x[0]
        tria.vert.y[ID-1] = y[0] 
        if x[0] == -1:
            print('value of -1 found at')
            print(ID)

    # Read in cell data
    #------------------
    i = 0
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1    

    # Read number of cells and cell vertices
    values = alllines[i].split()
    nc = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nc = nc[0]

    # Initialize vertex structure
    tria.cell.Initialize(nc)

    # Read cell data
    while i < len(alllines):
        if "ID, v1, v2, v3" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1    

    # Start reading
    for j in np.arange(0, nc):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        v1 = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        v2 = np.fromstring(values[2], dtype=int, count=1, sep =' ')
        v3 = np.fromstring(values[3], dtype=int, count=1, sep =' ')
      
        tria.cell.ID[ID-1] = ID
        tria.cell.vert[ID-1, 0] = v1[0] # Need to account for 0-based indexing
        tria.cell.vert[ID-1, 1] = v2[0]
        tria.cell.vert[ID-1, 2] = v3[0]

    # Return
    return tria


#--------------------------------------------------------------------------#
#                           MAGNETIC FIELD                                 #
#--------------------------------------------------------------------------#

def ReadRZPsiFile(filepath):
    # Description
    #------------
    # Read in an rzpsi file and return the coordinate vectors R, Z and
    # the values Psi, where Psi[i, j] yields the Psi value at R[i], Z[j]

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # First line is empty
    lind = 1

    # Read R, Z sizes
    values = alllines[lind].split()
    nR = np.fromstring(values[0], dtype=int, count=1, sep =' '); nR = nR[0]
    nZ = np.fromstring(values[1], dtype=int, count=1, sep =' '); nZ = nZ[0]

    # Skip the next two lines
    lind = lind + 2

    # Read R coordinate
    lind = lind + 1
    R = np.zeros(nR)
    k = 0
    while k < nR:
        # Read line
        values = alllines[lind].split()
        for j in values:
            R[k] = j 
            k = k + 1

        # Go to next line
        lind = lind + 1

    # Skip the next two lines
    lind = lind + 2

    # Read Z coordinate
    Z = np.zeros(nZ)
    k = 0
    while k < nZ:
        # Read line
        values = alllines[lind].split()
        for j in values:
            Z[k] = j 
            k = k + 1

        # Go to next line
        lind = lind + 1

    # Skip next two lines
    lind = lind + 2

    # Read Psi values
    Psi = np.zeros([nR, nZ])
    i = 0
    j = 0
    k = 0
    while k < nR*nZ:
        # Read
        values = alllines[lind].split() 
        for p in values:
            # Add value
            Psi[i, j] = p
            
            # Determine new i, j indices
            if (i == nR-1):
                i = 0
                j = j + 1
            else:
                i = i + 1

            k = k + 1

        # Update couner
        lind = lind + 1

    # Return
    return [R, Z, Psi]

def WriteRZPsiFile(dirpath, R, Z, Psi):
    # Description
    #------------
    # Write out an rzpsi.dat file to the filepath specified as dirpath
    # (should contain the full name)

    # Define number of values per line
    nvpl = 6

    # Open file
    thisfile = open(dirpath, "w")

    # Write initial white line
    thisfile.write("\n")

    # Write sizes
    thisfile.write(str(R.size) + "    " + str(Z.size) + "\n")

    # R coordinate data
    thisfile.write("$r \n")
    thisfile.write("nr=" + str(R.size) + "\n")
    k = 0
    while k < R.size:
        thisstr = ""
        if k + nvpl < R.size:
            tv = np.zeros(nvpl)
        else:
            tv = np.zeros(R.size-k)
        for i in np.arange(0, tv.size, 1):
            thisstr = thisstr + np.format_float_scientific(R[k], precision=8) + "  "
            k = k + 1
        thisstr = thisstr + "\n"
        thisfile.write(thisstr)

    # Z coordinate data
    thisfile.write("$z \n")
    thisfile.write("nz=" + str(Z.size) + "\n")
    k = 0
    while k < Z.size:
        thisstr = ""
        if k + nvpl < Z.size:
            tv = np.zeros(nvpl)
        else:
            tv = np.zeros(Z.size-k)
        for i in np.arange(0, tv.size, 1):
            thisstr = thisstr + np.format_float_scientific(Z[k], precision=8) + "  "
            k = k + 1
        thisstr = thisstr + "\n"
        thisfile.write(thisstr)

    # Psi data
    thisfile.write("\n")
    thisfile.write("$psi \n")
    thispsi = np.reshape(Psi, Psi.size, order='F')
    k = 0
    while k < thispsi.size:
        thisstr = ""
        if k + nvpl < thispsi.size:
            tv = np.zeros(nvpl)
        else:
            tv = np.zeros(thispsi.size-k)
        for i in np.arange(0, tv.size, 1):
           thisstr = thisstr + np.format_float_scientific(thispsi[k], precision=8) + "  "
           k = k + 1
        thisstr = thisstr + "\n"
        thisfile.write(thisstr)

    # Close the file
    thisfile.close()      

def ReadRZPsiFromEqdskFile(filepath):
    # Description
    #------------
    # Read the R, Z, Psi values (in Weber) from an eqdsk file and 
    # convert them to classic Weber/rad by dividing through 2*pi. We 
    # have to exploit the convention that each float has 16 characters,
    # with no whitespace between (except if the first character, which
    # indicates the sign, is positive, then there is a space). 

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Initialize line position
    i = 0

    # Read first line: comment, dummy, nR, nZ 
    values = alllines[i].split()
    nR = np.fromstring(values[-2], dtype=int, count=1, sep =' '); nR = nR[0]
    nZ = np.fromstring(values[-1], dtype=int, count=1, sep =' '); nZ = nZ[0]
    i = i + 1

    # Read next numbers
    values = SplitEqdskFileLine(alllines[i])
    rboxlength = np.fromstring(values[0], dtype=float, count=1, sep =' ')[0] # length of box in R direction
    zboxlength = np.fromstring(values[1], dtype=float, count=1, sep =' ')[0] # length of box in Z direction
    R0 = np.fromstring(values[2], dtype=float, count=1, sep =' ')[0]  # Major radius
    rboxleft = np.fromstring(values[3], dtype=float, count=1, sep =' ')[0] # Start of R coordinate
    zboxmid = np.fromstring(values[4], dtype=float, count=1, sep =' ')[0] # Mid of Z coordinate
    i = i + 1

    # Construct coordinates
    R = rboxleft + rboxlength/(nR-1)*np.arange(0, nR, 1, dtype=float)
    Z = zboxmid - zboxlength/2.0 + zboxlength/(nZ-1)*np.arange(0, nZ, 1, dtype=float)

    # Read the magnetic axit
    values = SplitEqdskFileLine(alllines[i])
    Rpsi0 = np.fromstring(values[0], dtype=float, count=1, sep =' ')[0] #
    Zpsi0 = np.fromstring(values[1], dtype=float, count=1, sep =' ')[0] # 
    PsiaxisVs = np.fromstring(values[2], dtype=float, count=1, sep =' ')[0]
    PsiedgeVs = np.fromstring(values[3], dtype=float, count=1, sep =' ')[0]
    BtoratR0 = np.fromstring(values[4], dtype=float, count=1, sep =' ')[0]
    i = i + 1

    # Skip two lines - currents and stuff
    i = i + 2

    # Skip a lot of profiles
    i = i + int(np.ceil(4*nR/5.0)) # because we know there are five numbers per line

    # Read psi in V*s = Wb
    # Read Psi values
    Psi = np.zeros([nR, nZ])
    ii = 0
    j = 0
    k = 0
    while k < nR*nZ:
        # Read
        values = SplitEqdskFileLine(alllines[i])
        for p in values:
            # Add value
            Psi[ii, j] = p
            
            # Determine new i, j indices
            if (ii == nR-1):
                ii = 0
                j = j + 1
            else:
                ii = ii + 1

            k = k + 1

        # Update counter
        i = i + 1

    # Divide by 2*Pi
    Psi = Psi/(2.0*np.pi)

    # Return
    return R, Z, Psi

def SplitEqdskFileLine(line):
    # Description
    #------------
    # This routine splits a line from an eqdsk file (except for the header)
    # by simply taking for each value the first 16 characters. It also 
    # assumes that there are 5*16 characters per line. Very ugly. 
    # Probably not general at all.
    values = []
    cc = 16
    while (cc < len(line)):
        values.append(line[cc-16:cc].split()[0])
        cc = cc + 16

    return values

def ReadRZPsiFromCSV(Rfilepath, Zfilepath, Psifilepath, separator):
    # Description
    #------------
    # Read the rzpsi values from three different csv files where R, Z 
    # and psi are given all as nZ-by-nR matrices (so nZ rows, nR columns) 
    # These are returned, however, as vectors and as a nR-by-nZ array
    import csv 

    # Open psi file
    thisfile = open(Psifilepath)

    # Read lines
    alllines = thisfile.readlines()

    # Determine dimensions
    nZ = len(alllines)
    nR = len(alllines[0].split(separator))

    # Initialize
    psi = np.zeros((nR, nZ), dtype=float)
    R = np.zeros(nR, dtype=float)
    Z = np.zeros(nZ, dtype=float)

    # Read psi
    for j in np.arange(0, nZ):
        temp = alllines[j].split(separator)
        for i in np.arange(0, nR):
            psi[i, j] = np.fromstring(temp[i], dtype=float, count=1, sep =' ')

    # Close file
    thisfile.close()

    # Open R file
    thisfile = open(Rfilepath)

    # Read lines
    alllines = thisfile.readlines()

    # Determine dimensions
    tempnZ = len(alllines)
    tempnR = len(alllines[0].split(separator))
    if (tempnZ == nR and tempnR == 1):
        for i in np.arange(0, nR):
            R[i] = np.fromstring(alllines[i], dtype=float, count=1, sep=' ')
    elif ((tempnZ == nZ and tempnR == nR) or (tempnZ == nZ and tempnR == 1)):
        temp = alllines[1].split(separator)
        for i in np.arange(0, nR):
            R[i] = np.fromstring(temp[i], dtype=float, count=1, sep=' ')
    else: 
        print('Trying transposing R coordinates')
        alllines = [' '.join(row) for row in zip(*[line.split(separator) for line in alllines])]
        tempnZ = len(alllines)
        tempnR = len(alllines[0].split(separator))
        if (tempnZ == nR and tempnR == 1):
            for i in np.arange(0, nR):
                R[i] = np.fromstring(alllines[i], dtype=float, count=1, sep=' ')
        elif ((tempnZ == nZ and tempnR == nR) or (tempnZ == nZ and tempnR == 1)):
            temp = alllines[1].split(separator)
            for i in np.arange(0, nR):
                R[i] = np.fromstring(temp[i], dtype=float, count=1, sep=' ')
        else:
            # Shouldn't be happening
            raise ValueError('Inconsistent dimensions of R coordinates') 

    
    # close file
    thisfile.close()

    # Open Z file
    thisfile = open(Zfilepath)

    # Read lines
    alllines = thisfile.readlines()

    # Determine dimensions
    tempnZ = len(alllines)
    tempnR = len(alllines[0].split(separator))
    if ((tempnZ == nZ and tempnR == nR) or (tempnZ == nZ and tempnR == 1)):
        for i in np.arange(0, nZ):
            Z[i] = np.fromstring(alllines[i], dtype=float, count=1, sep=' ')
    elif (tempnZ == 1 and tempnR == nZ):
        temp = alllines[1].split(separator)
        for i in np.arange(0, nZ):
            Z[i] = np.fromstring(temp[i], dtype=float, count=1, sep=' ')
    else: 
        print('Trying transposing Z coordinates')
        alllines = [' '.join(row) for row in zip(*[line.split(separator) for line in alllines])]
        tempnZ = len(alllines)
        tempnR = len(alllines[0].split(separator))
        if ((tempnZ == nZ and tempnR == nR) or (tempnZ == nZ and tempnR == 1)):
            for i in np.arange(0, nZ):
                Z[i] = np.fromstring(alllines[i], dtype=float, count=1, sep=' ')
        elif (tempnZ == 1 and tempnR == nZ):
            temp = alllines[1].split(separator)
            for i in np.arange(0, nZ):
                Z[i] = np.fromstring(temp[i], dtype=float, count=1, sep=' ')
        else: 
            # Shouldn't be happening
            raise ValueError('Inconsistent dimensions of Z coordinates')

    
    # close file
    thisfile.close()
    


    return [R, Z, psi, nZ, nR]

#--------------------------------------------------------------------------#
#                                SOLPS                                     #
#--------------------------------------------------------------------------#   

# Read field 
def ReadSOLPSField(fieldname, lines, startindex, nval, fieldtype, dims):

    # Description
    #------------
    # Read in a field with the name 'fieldname' that has an amount of 
    # values 'nval'. 'startindex' is the starting line index in the 
    # lines. At exit, 'startindex' will be the last line read in. 

    # Initialize (always skip header)
    if startindex >= len(lines):
        startindex = 1
    elif startindex < 0:
        startindex = 1
    i = startindex 
    if fieldtype == 'float':
        vals = np.zeros(nval, dtype=float)
    elif fieldtype == 'int':
        vals = np.zeros(nval, dtype=int)
    else: 
        raise Exception("ReadField: unknown fieldtype option")

    # Look for the field name
    isfound = False 
    while i < len(lines):
        if '*cf' in lines[i]:
            values = lines[i].split()
            for j in values: 
                if fieldname == j:
                    isfound = True 
            if isfound:
                break 
            else:
                i = i + 1
        else: 
            i = i + 1
    if not isfound:
        i = 1 
        while i < startindex: 
            if '*cf' in lines[i]:
                values = lines[i].split()
                for j in values: 
                    if fieldname == j:
                        isfound = True 
                if isfound:
                    break 
                else: 
                    i = i + 1
            else: 
                i = i + 1
    
    # Check if the field was found
    if not isfound:
        print("ReadField: could not find field: '" + fieldname + "' in "
            "in file. Returning zeros...")
    else: 
        # Read in the field
        i = i + 1 # Skip header
        k = 0
        while k < nval:
            # Read line
            values = lines[i].split()
            for j in values:
                vals[k] = j 
                k = k + 1

            # Go to next line
            i = i + 1

    # Reshape the field to the desired dimensions
    out = np.reshape(vals, dims, order='F')

    # Return 
    startindex = i
    return out, startindex # need to return startindex, since it is immutable and therefore isn't returned...

# Write field 
def WriteSOLPSField(filehandle, vals, fieldname, fieldtype):
    # Description
    #------------
    # Write out a scalar, unstructured field with the name 'fieldname' 
    # into the file opened by the filehandle. This includes the classical
    # field name header, followed by the solps format of fields. 
    # For fields of multiple species, the values are appended columnwise
    # as expected.

    # Write header
    nval = vals.size
    headerstring = ( "*cf:".ljust(8, ' ') + fieldtype.ljust(8, ' ') + 
        str(nval).ljust(12, ' ') + " ".ljust(4, ' ') + fieldname )
    filehandle.write(headerstring + '\n')
    val = np.reshape(vals, nval, order='F')
    nvalsperline = 6 # max. number of values on each line

    # Write data
    if fieldtype == 'real':
        k = 0
        while k < nval:
            valind = 0 
            datastring = ""
            while (valind < nvalsperline) and (k+valind < nval):
                datastring = datastring + np.format_float_scientific(val[k+valind], 
                    unique=False, precision=13, pad_left=3, exp_digits=3, trim='k')
                valind = valind + 1
            filehandle.write(datastring + ' \n')
            k = k + nvalsperline # no problem if k > nval 
    elif fieldtype == 'int':
        k = 0
        while k < nval:
            valind = 0 
            datastring = ""
            while (valind < nvalsperline) and (k+valind < nval):
                datastring = datastring + str(val[k+valind]).rjust(11, " ")
                valind = valind + 1
            filehandle.write(datastring + ' \n')
            k = k + nvalsperline # no problem if k > nval 
    else:
        raise Exception("WriteSOLPSField: unknown data format: " + fieldtype)

# Read fort.78 file
def ReadFort78File(filepath):
    # Description
    #------------
    # This routine reads a formatted fort.78 file where the polygon
    # coordinates (2D) are stored as [x, y] columns and are prepended
    # with the amount of points. Polygons in output are separated by
    # NaNs

    # Note: original values are in cm but are converted to m here!

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0:1]
    vals = np.zeros([len(alllines), 2])

    # Read in vertex data
    cc = 0
    for i in alllines:
        # Read in values
        temp = np.fromstring(i, dtype=float, count=-1, sep=' ')

        # Check the size
        if len(temp) != 3:
            # Insert NaN
            vals[cc, :] = np.nan
        else:
            # Read
            vals[cc, 0:2] = np.fromstring(i, dtype=float, count=2, sep=' ')

        # Update counter
        cc = cc + 1

    # Change units from cm to m
    vals = vals/100.0

    # Return values
    return vals[0:cc, 0:2]
