# This module defines some convenient types for post-processing that
# mimick the types of goat 
import numpy as np 
import shapely
from shapely import geometry

#======================================================================#
#                                                                      #
#                       GRID GENERATOR TYPES                           #
#                                                                      #
#======================================================================#

#----------------------------------------------------------------------#
#                           DEFINITIONS                                #
#----------------------------------------------------------------------#
# Vessel parts definitions
targetID = 1
coreID = 2
outerboundaryID = 3
vesselID = 4
interiorID = 5

# Definitions for grid generation
# Topological mesh vertex IDs (field minimum, saddle point, maximum: 1, 2, 3,
# field & boundary tangency point types 1 and 2: 4, 5, boundary vertex
# that is no tangency point or whatever: 6, vertex constructed by 
# splitting face: -1)
TMvertexminID = 1
TMvertexsaddleID = 2
TMvertexmaxID = 3
TMvertextp1ID = 4
TMvertextp2ID = 5
TMvertexbndID = 6
TMvertexsplitID = -1
TMvertexregularID = 0

# Topological mesh boundary IDs (1: radial face, 2: poloidal face, &
# 3: boundary face)
TMfaceradID = 1
TMfacepolID = 2
TMfacebndID = 3
TMfacesepID = 4
TMfacecoreID = 5
TMfacePFID = 6
TMfacealbndID = 7
TMfaceSOLID = 8

#----------------------------------------------------------------------#
#                               I/O                                    #
#----------------------------------------------------------------------#

# Structure format
class Structure:
    def __init__(self):
        # Total number of coordinates (negative if open contour)
        self.n = 0

        # Coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # ID
        self.ID = 0

    def Initialize(self, ntot, x, y, ID):
        assert ntot == len(x)
        assert ntot == len(y)
        self.n = -ntot 
        if x[0] == x[len(x)-1] and y[0] == y[len(y)-1]:
            self.n =  ntot 
        self.x = x 
        self.y = y
        self.ID = ID

#----------------------------------------------------------------------#
#                        TOPOLOGICAL MESH                              #
#----------------------------------------------------------------------#

# Topological mesh vertices
class TopomeshVert:
    # Definition
    def __init__(self):
        # total number of vertices
        self.ntot = 0

        # Coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # ID
        self.ID = np.zeros(0, dtype=int)
        self.type = np.zeros(0, dtype=int)
        self.fsID = np.zeros(0, dtype=int)

        # Field value
        self.fval = np.zeros(0, dtype=float)

        # Boundary vertex
        self.BV = np.zeros(0, dtype=int)

    # Initialization
    def Initialize(self, ntot):
        # total number of vertices
        self.ntot = ntot

        # Coordinates
        self.x = np.zeros(ntot, dtype=float)
        self.y = np.zeros(ntot, dtype=float)

        # ID
        self.ID = np.zeros(ntot, dtype=int)
        self.type = np.zeros(ntot, dtype=int)
        self.fsID = np.zeros(ntot, dtype=int)

        # Field value
        self.fval = np.zeros(ntot, dtype=float)

        # Boundary vertex
        self.BV = np.zeros(ntot, dtype=int)


# Topological mesh face data
class TopomeshFacedata: 
    # Definition
    def __init__(self):
        # coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)
        self.dx = np.zeros(0, dtype=float)
        self.dy = np.zeros(0, dtype=float)
        self.dl = np.zeros(0, dtype=float)
        self.dlsum = np.zeros(0, dtype=float)
        self.L = 0

    # Initializer
    def Initialize(self, nc):
        # Initialize coordinates
        self.x = np.zeros(nc, dtype=float)
        self.y = np.zeros(nc, dtype=float)
        self.ComputeMetrics()

    def ComputeMetrics(self):
        self.dx = np.diff(self.x)
        self.dy = np.diff(self.y)
        self.dl = np.sqrt(self.dx**2 + self.dy**2)
        self.dlsum = np.cumsum(np.append(np.array([0]), (self.dl)))
        self.L = np.sum(self.dl)
        
    # Coordinate interpolator
    def InterpolateCoordinates(self, frac):
        # Checks
        if (frac < 0.0) or (frac > 1.0):
            raise ValueError("InterpolateCoordinates: frac should be between 0.0" \
                "(start of line) and 1.0 (end of line)")
        
        # Interpolate
        lcoord = frac*self.L 
        xc = np.interp(lcoord, self.dlsum, self.x)
        yc = np.interp(lcoord, self.dlsum, self.y)
        return xc, yc

# Topological mesh faces
class TopomeshFace:
    # Definition
    def __init__(self):
        # total number of faces
        self.ntot = 0

        # Coordinates
        self.data = TopomeshFacedata()

        # ID
        self.ID = np.zeros(0, dtype=int)

        # Field value
        self.fval = np.zeros(0, dtype=float)

        # Boundary face
        self.BF = np.zeros(0, dtype=int)

        # Flux surface ID
        self.fsID = np.zeros(0, dtype=int)
        
        # Type
        self.type = np.zeros(0, dtype=int)
        self.vert = np.zeros((0, 2), dtype=int)

        # Number of coordinates
        self.nc = np.zeros(0, dtype=int)

    # Initialization
    def Initialize(self, ntot):
        # total number of faces
        self.ntot = ntot

        # Coordinates
        self.data = [TopomeshFacedata() for i in range(ntot)]

        # ID
        self.ID = np.zeros(ntot, dtype=int)

        # Field value
        self.fval = np.zeros(ntot, dtype=float)

        # Boundary face
        self.BF = np.zeros(ntot, dtype=int)

        # Flux surface ID
        self.fsID = np.zeros(ntot, dtype=int)
        
        # Type
        self.type = np.zeros(ntot, dtype=int)
        self.vert = np.zeros((ntot, 2), dtype=int)

        # Number of coordinates
        self.nc = np.zeros(ntot, dtype=int)


    def AddFaceCoordinates(self, faceindex, xf, yf):
        # Add face coordinates
        self.data[faceindex].x = xf 
        self.data[faceindex].y = yf

# Topological mesh cells
class TopomeshCell:
    # Definition
    def __init__(self):
        # Total number of cells
        self.ntot = 0
        self.nvert = 0
        self.nface = 0

        # ID
        self.ID = np.zeros(0, dtype=int)

        # Number of coordinates
        self.nc = np.zeros(0, dtype=int)

        # Vertices
        self.vert = np.zeros(0, dtype=int)
        self.vertP = np.zeros((0, 2), dtype=int)

        # Faces
        self.face = np.zeros(0, dtype=int)
        self.faceP = np.zeros((0, 2), dtype=int)

        # Coordinates
        self.data = TopomeshFacedata() 

    # Initialization
    def Initialize(self, ntot, nvert, nface):
        # Total number of cells
        self.ntot = ntot
        self.nvert = nvert 
        self.nface = nface

        # ID
        self.ID = np.zeros(ntot, dtype=int)

        # Number of polygon coordinates
        self.nc = np.zeros(ntot, dtype=int)

        # Vertices
        self.vert = np.zeros(nvert, dtype=int)
        self.vertP = np.zeros((ntot, 2), dtype=int)

        # Faces
        self.face = np.zeros(nface, dtype=int)
        self.faceP = np.zeros((ntot, 2), dtype=int)

        # Coordinates
        self.data = [TopomeshFacedata() for i in range(ntot)]

    def AddCellCoordinates(self, cellindex, xc, yc):
        # Add face coordinates
        self.data[cellindex].x = xc 
        self.data[cellindex].y = yc

# Topological mesh flux surface
class TopomeshFluxSurface:
    def __init__(self):
        # Dimensions
        self.ntot = 0
        
        # IDs etc
        self.ID = np.zeros(self.ntot, dtype=int)
        self.psi = np.zeros(self.ntot, dtype=float)

    def Initialize(self, ntot):
        # Dimensions
        self.ntot = ntot
        
        # IDs etc
        self.ID = np.zeros(self.ntot, dtype=int)
        self.psi = np.zeros(self.ntot, dtype=float)

# Topological mesh tubes
class TopomeshFluxTube:
    def __init__(self):
        # Dimensions
        self.ntot   = 0
        self.nface  = 0
        self.ncell  = 0

        # ID
        self.ID     = np.zeros(self.ntot, dtype=int)

        # Faces
        self.face   = np.zeros(self.nface, dtype=int)
        self.faceP = np.zeros((self.ntot, 2), dtype=int)

        # Cells
        self.cell   = np.zeros(self.ncell, dtype=int)
        self.cellP  = np.zeros((self.ntot, 2), dtype=int)

    def Initialize(self, ntot, nface, ncell):
        # Dimensions
        self.ntot   = ntot 
        self.nface  = nface 
        self.ncell  = ncell 

        # ID
        self.ID     = np.zeros(self.ntot, dtype=int)

        # Faces
        self.face   = np.zeros(self.nface, dtype=int)
        self.faceP = np.zeros((self.ntot, 2), dtype=int)

        # Cells
        self.cell   = np.zeros(self.ncell, dtype=int)
        self.cellP  = np.zeros((self.ntot, 2), dtype=int)

    # Face getter
    def GetFace(self, i):
        return self.face[self.faceP[i, 0]:self.faceP[i, 0]+self.faceP[i, 1]]
    
    # Cell getter
    def GetCell(self, i):
        return self.cell[self.cellP[i, 0]:self.cellP[i, 0]+self.cellP[i, 1]]
        

# Topological mesh
class Topomesh:
    def __init__(self):
        # Fields
        self.vert   = TopomeshVert()
        self.face   = TopomeshFace() 
        self.cell   = TopomeshCell()
        self.fs     = TopomeshFluxSurface()
        self.ft     = TopomeshFluxTube()
        
#----------------------------------------------------------------------#
#                         GRID GENERATOR                               #
#----------------------------------------------------------------------#

# GGTM lines
class GGTMLine:
    # Definition
    def __init__(self):
        # coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Vertices
        self.vert = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nv):
        # Initialize coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # Vertices
        self.vert = np.zeros(nv, dtype=int)

    # Coordinate addition
    def AddCoordinates(self, xc, yc, vc):
        self.x = xc 
        self.y = yc 
        self.vert = vc

# GGTM tubes
class GGTMTube: 
    # Definition
    def __init__(self):
        # Lines
        self.hfline = GGTMLine()
        self.lfline = GGTMLine()

    # Initializer
    def Initialize(self, hfvert, hfx, hfy, lfvert, lfx, lfy):
        # Initialize lines
        self.hfline.AddCoordinates(hfx, hfy, hfvert)
        self.lfline.AddCoordinates(lfx, lfy, lfvert)

# GGTM faces
class GGTMFace:
    # Definition
    def __init__(self):
        # coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Vertices
        self.vert = np.zeros(0, dtype=int)

        # Number
        self.nv = 0

        # ID
        self.ID = 0

    # Initializer
    def Initialize(self, nv, ID):
        # Initialize coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # Vertices
        self.vert = np.zeros(nv, dtype=int)

        # Number
        self.nv = nv

        # ID
        self.ID = ID

# GGTM cells
class GGTMCell:
    # Definition
    def __init__(self):

        # hfline, lfline
        self.hfline = GGTMLine()
        self.lfline = GGTMLine()

        # Start and end radial face IDs
        self.srf = 0
        self.erf = 0

        # Lines
        self.nl = 0
        self.lines = GGTMLine()
        self.tubes = GGTMTube()

        # ID
        self.ID = 0

    # Initialization
    def Initialize(self, nl, ID):

        # Coordinates
        self.nl = nl
        self.lines = [GGTMLine() for i in range(nl)]
        self.tubes = [GGTMTube() for i in range(nl+1)] # need to account for lfline and hfline being stored differently

        # ID
        self.ID = ID

    # Adding line coordinates
    def AddLineCoordinates(self, lineindex, xl, yl, vl):
        # Add face coordinates
        self.lines[lineindex].AddCoordinates(xl, yl, vl)

    # Adding tube lines
    def AddTube(self, tubeindex, hfx, hfy, hfvert, lfx, lfy, lfvert):
        # Add tube line coordinates
        self.tubes[tubeindex].Initialize(hfvert, hfx, hfy, lfvert, lfx, lfy)

# GGTM data class
class GGTMData:
    # Init
    def __init__(self):
        # Fields
        self.face = GGTMFace() 
        self.cell = GGTMCell()

    # Initializer
    def Initialize(self, nf, nc):
        # Fields
        self.face = [GGTMFace() for i in range(nf)]
        self.cell = [GGTMCell() for i in range(nc)]
        
# Grid vertices
class GGVert:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        nv = self.ntot

        # coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(nv, dtype=int)

    # Initializer
    def Initialize(self, nv):
        # Number
        self.ntot = nv 

        # Coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(nv, dtype=int)

# Grid faces
class GGFace:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0 

        # Vertices
        self.v1 = np.zeros(0, dtype=int)
        self.v2 = np.zeros(0, dtype=int)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Label
        self.label = np.zeros(0, dtype=int)
        self.region = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nf):
        # Number
        self.ntot = nf

        # Vertices
        self.v1 = np.zeros(nf, dtype=int)
        self.v2 = np.zeros(nf, dtype=int)

        # ID
        self.ID =  np.zeros(nf, dtype=int)

        # Label
        self.label = np.zeros(nf, dtype=int)
        self.region = np.zeros(nf, dtype=int)

# Grid cells
class GGCell:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nvert = 0
        nc = self.ntot 

        # Vertex pointer
        self.vp1 = np.zeros(nc, dtype=int)
        self.vp2 = np.zeros(nc, dtype=int)
        
        # Vertices
        self.vert = np.zeros(nc, dtype=int)

        # ID
        self.ID =  np.zeros(nc, dtype=int)

        # Region
        self.region = np.zeros(nc, dtype=int)

    # Initializer
    def Initialize(self, nc, ncv):
        # Number
        self.ntot = nc
        self.nvert = ncv

        # Vertex pointer
        self.vp1 = np.zeros(nc, dtype=int)
        self.vp2 = np.zeros(nc, dtype=int)
        
        # Vertices
        self.vert = np.zeros(ncv, dtype=int)

        # ID
        self.ID =  np.zeros(nc, dtype=int)

        # Region
        self.region = np.zeros(nc, dtype=int)

    # Vertex getter
    def GetVert(self, i):
        return self.vert[self.vp1[i]:self.vp1[i]+self.vp2[i]]

# Grid
class GGGrid:
    # Init
    def __init__(self):
        # Fields
        self.vert = GGVert()
        self.face = GGFace() 
        self.cell = GGCell()

#----------------------------------------------------------------------#
#                         GRID ADAPTATION                              #
#----------------------------------------------------------------------#
# GAGrid vertices
class GAVert:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        nv = self.ntot

        # coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(nv, dtype=int)

    # Initializer
    def Initialize(self, nv):
        # Number
        self.ntot = nv 

        # Coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(nv, dtype=int)

# GAGrid faces
class GAFace:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0 

        # Vertices
        self.v1 = np.zeros(0, dtype=int)
        self.v2 = np.zeros(0, dtype=int)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Label
        self.label = np.zeros(0, dtype=int)
        self.region = np.zeros(0, dtype=int)
        self.aligned = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nf):
        # Number
        self.ntot = nf

        # Vertices
        self.v1 = np.zeros(nf, dtype=int)
        self.v2 = np.zeros(nf, dtype=int)

        # ID
        self.ID =  np.zeros(nf, dtype=int)

        # Label
        self.label = np.zeros(nf, dtype=int)
        self.region = np.zeros(nf, dtype=int)
        self.aligned = np.zeros(nf, dtype=int)

# GAGrid cells
class GACell:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nvert = 0
        nc = self.ntot 

        # Vertex pointer
        self.vp1 = np.zeros(nc, dtype=int)
        self.vp2 = np.zeros(nc, dtype=int)
        
        # Vertices
        self.vert = np.zeros(nc, dtype=int)

        # ID
        self.ID =  np.zeros(nc, dtype=int)

        # Region
        self.region = np.zeros(nc, dtype=int)

        # Coordinates
        self.x  = np.zeros(nc, dtype=float)
        self.y  = np.zeros(nc, dtype=float)

    # Initializer
    def Initialize(self, nc, ncv):
        # Number
        self.ntot = nc
        self.nvert = ncv

        # Vertex pointer
        self.vp1 = np.zeros(nc, dtype=int)
        self.vp2 = np.zeros(nc, dtype=int)
        
        # Vertices
        self.vert = np.zeros(ncv, dtype=int)

        # ID
        self.ID =  np.zeros(nc, dtype=int)

        # Region
        self.region = np.zeros(nc, dtype=int)

        # Coordinates
        self.x  = np.zeros(nc, dtype=float)
        self.y  = np.zeros(nc, dtype=float)        

    # Vertex getter
    def GetVert(self, i):
        return self.vert[self.vp1[i]:self.vp1[i]+self.vp2[i]]

class GAGrid:
    # Init
    def __init__(self):
        # Fields
        self.vert = GAVert()
        self.face = GAFace()
        self.cell = GACell()

#----------------------------------------------------------------------#
#                         TRIANGULATION                                #
#----------------------------------------------------------------------#
# triangulation vert
class TriaVert:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        nv = self.ntot

        # coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)


    # Initializer
    def Initialize(self, nv):
        # Number
        self.ntot = nv 

        # Coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

# triangulation cell        
class TriaCell:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nvert = 0
        nc = self.ntot

        # cvert
        shape = (nc, 3)
        self.vert = np.zeros(shape, dtype=int)

        # ID
        self.ID = np.zeros(nc, dtype=int)

    # Initializer
    def Initialize(self, nc):
         # Number
        self.ntot = nc
        self.nvert = nc*3

        # vert
        shape = (nc, 3)
        self.vert = np.zeros(shape, dtype=int)

        # ID
        self.ID = np.zeros(nc, dtype=int)   

class TriaGrid:
    # Init
    def __init__(self):
        # Fields
        self.vert = TriaVert()
        self.cell = TriaCell()        


#----------------------------------------------------------------------#
#                        SIMULATION GRID                               #
#----------------------------------------------------------------------#
        
# Grid vertices
class Vert:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        nv = self.ntot

        # Coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Magnetic field
        self.psi = np.zeros(0, dtype=float)
        self.bxv = np.zeros(0, dtype=float)
        self.byv = np.zeros(0, dtype=float)
        self.ffbz = np.zeros(0, dtype=float)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(0, dtype=int)

        # Metrics
        self.bb = np.zeros((nv, 4), dtype=float)
        self.ffbz = np.zeros(nv, dtype=float)
        self.fspsi = np.zeros(nv, dtype=float)

    # Initializer
    def Initialize(self, nv):
        # Number
        self.ntot = nv 

        # Coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # Magnetic field
        self.psi = np.zeros(nv, dtype=float)
        self.bx = np.zeros(nv, dtype=float)
        self.by = np.zeros(nv, dtype=float)
        self.ffbz = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(nv, dtype=int)

        # Cells (unsorted)
        self.cp1 = np.zeros(nv, dtype=int)
        self.cp2 = np.zeros(nv, dtype=int)
        self.cell = np.zeros(0, dtype=int) # To be determined in grid interconnections

        # Metrics
        self.bb = np.zeros((nv, 4), dtype=float)
        self.ffbz = np.zeros(nv, dtype=float)
        self.fspsi = np.zeros(nv, dtype=float)

    def GetCell(self, i):
        return self.cell[self.cp1[i]:self.cp1[i]+self.cp2[i]]

# Grid faces
class Face:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0 
        nf = self.ntot

        # Vertices
        self.v1 = np.zeros(0, dtype=int)
        self.v2 = np.zeros(0, dtype=int)

        # Cell neighbours
        self.nb1 = np.zeros(0, dtype=int)
        self.nb2 = np.zeros(0, dtype=int)

        # Labels
        self.label = np.zeros(0, dtype=int)
        self.region = np.zeros(0, dtype=int)
        
        # Magnetic field
        self.aligned = np.zeros(0, dtype=int)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Metrics
        self.bb = np.zeros((nf, 4), dtype=float)
        self.s = np.zeros(nf, dtype=float)
        self.hc = np.zeros((nf, 4), dtype=float)
        self.ht = np.zeros(nf, dtype=float)
        self.qgam = np.zeros((nf, 2), dtype=float)
        self.qalf = np.zeros((nf, 2), dtype=float)
        self.qbet = np.zeros((nf, 2), dtype=float)
        self.pbs = np.zeros(nf, dtype=float)

    # Initializer
    def Initialize(self, nf):
        # Number
        self.ntot = nf

        # Vertices
        self.v1 = np.zeros(nf, dtype=int)
        self.v2 = np.zeros(nf, dtype=int)

        # Cells
        self.nb1 = np.zeros(nf, dtype=int)
        self.nb2 = np.zeros(nf, dtype=int)

        # Labels
        self.label = np.zeros(nf, dtype=int)
        self.region = np.zeros(nf, dtype=int)
        
        # Magnetic field
        self.aligned = np.zeros(nf, dtype=int)

        # ID
        self.ID =  np.zeros(nf, dtype=int)

        # Coordinates
        self.x = np.zeros(nf, dtype=float)
        self.y = np.zeros(nf, dtype=float)

        # Metrics
        self.bb = np.zeros((nf, 4), dtype=float)
        self.s = np.zeros(nf, dtype=float)
        self.hc = np.zeros((nf, 4), dtype=float)
        self.ht = np.zeros(nf, dtype=float)
        self.qgam = np.zeros((nf, 2), dtype=float)
        self.qalf = np.zeros((nf, 2), dtype=float)
        self.qbet = np.zeros((nf, 2), dtype=float)
        self.pbs = np.zeros(nf, dtype=float)

# Grid cells
class Cell:
    # Definition
    def __init__(self):
        # Number
        self.nvert = 0
        self.nface = 0
        self.ncg = 0 # number of guard cells 
        self.nci = 0 # number of internal (non-guard) cells
        self.ntot = self.ncg + self.nci 
        nc = self.ntot

        # Vertex pointer
        self.vp1 = np.zeros(0, dtype=int)
        self.vp2 = np.zeros(0, dtype=int)

        # Face pointer
        self.fp1 = np.zeros(0, dtype=int)
        self.fp2 = np.zeros(0, dtype=int)

        # Neighbour pointer
        self.nbp1 = np.zeros(0, dtype=int)
        self.nbp2 = np.zeros(0, dtype=int)
        
        # Vertices
        self.vert = np.zeros(0, dtype=int)

        # Faces
        self.face = np.zeros(0, dtype=int)

        # Neighbours
        self.nb = np.zeros(0, dtype=int)

        # Coordinates
        self.x  = np.zeros(0, dtype=float)
        self.y  = np.zeros(0, dtype=float)

        # Magnetic Field
        self.bt = np.zeros(0, dtype=float)
        self.bp = np.zeros(0, dtype=float)
        self.psi   = np.zeros(0, dtype=float)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Region etc
        self.cflags = np.zeros(0, dtype=int)
        self.region = np.zeros(0, dtype=int)
        self.ft     = np.zeros(0, dtype=int)

        # Metrics
        self.bb = np.zeros((nc, 4), dtype=float)
        self.zb = np.zeros((nc, 3), dtype=float)
        self.sz = np.zeros(nc, dtype=int)
        self.hz = np.zeros(nc, dtype=int)
        self.hx = np.zeros(nc, dtype=int)
        self.qgam = np.zeros((nc, 2), dtype=int)
        self.vol = np.zeros(nc, dtype=int)

    # Initializer
    def Initialize(self, nci, ncg, ncv, ncf):
        # Number
        self.nci = nci 
        self.ncg = ncg 
        self.nvert = ncv
        self.nface = ncf
        nc = self.nci + self.ncg 
        self.ntot = nc

        # Vertex pointer
        self.vp1 = np.zeros(nc, dtype=int)
        self.vp2 = np.zeros(nc, dtype=int)

        # Face pointer
        self.fp1 = np.zeros(nc, dtype=int)
        self.fp2 = np.zeros(nc, dtype=int)

        # Neighbour pointer 
        self.nbp1 = np.zeros(nc, dtype=int)
        self.nbp2 = np.zeros(nc, dtype=int)
        
        # Vertices
        self.vert = np.zeros(ncv, dtype=int)

        # Faces
        self.face = np.zeros(ncf, dtype=int)

        # Neighbours
        self.nb = np.zeros(ncf, dtype=int)

        # Coordinates
        self.x  = np.zeros(nc, dtype=float)
        self.y  = np.zeros(nc, dtype=float)

        # Magnetic Field
        self.bt = np.zeros(nc, dtype=float)
        self.bp = np.zeros(nc, dtype=float)
        self.psi   = np.zeros(nc, dtype=float)

        # ID
        self.ID =  np.zeros(nc, dtype=int)

        # Region etc
        self.cflags = np.zeros(nc, dtype=int)
        self.region = np.zeros(nc, dtype=int)
        self.ft     = np.zeros(nc, dtype=int)

        # Metrics
        self.bb = np.zeros((nc, 4), dtype=float)
        self.zb = np.zeros((nc, 3), dtype=float)
        self.sz = np.zeros(nc, dtype=int)
        self.hz = np.zeros(nc, dtype=int)
        self.hx = np.zeros(nc, dtype=int)
        self.qgam = np.zeros((nc, 2), dtype=int)
        self.vol = np.zeros(nc, dtype=int)

    # Vertex getter
    def GetVert(self, i):
        return self.vert[self.vp1[i]:self.vp1[i]+self.vp2[i]]
    
    # Face getter
    def GetFace(self, i):
        return self.face[self.fp1[i]:self.fp1[i]+self.fp2[i]]

    # Cell neighbour getter
    def GetNeig(self, i):
        return self.nb[self.nbp1[i]:self.nbp1[i]+self.nbp2[i]]
        

# Grid flux surfaces
class FluxSurf:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nface = 0

        # Faces
        self.fp1 = np.zeros(0, dtype=int)
        self.fp2 = np.zeros(0, dtype=int)
        self.psi = np.zeros(0, dtype=float)
        self.face = np.zeros(0, dtype=int)

        # ID
        self.ID = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nfs):
        # Number
        self.ntot = nfs 

        # Faces
        self.fp1 = np.zeros(nfs, dtype=int)
        self.fp2 = np.zeros(nfs, dtype=int)

        # ID
        self.ID = np.zeros(nfs, dtype=int)

        # Magnetic field
        self.psi = np.zeros(nfs, dtype=float)

    def InitializeFaceData(self, nfsf):
        self.nface = nfsf 
        self.face = np.zeros(nfsf, dtype=int)

    # Face getter
    def GetFace(self, i):
        return self.face[self.fp1[i]:self.fp1[i]+self.fp2[i]]
    
# Grid flux tubes
class FluxTube:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nface = 0
        self.ncell = 0

        # Faces
        self.fp1 = np.zeros(0, dtype=int)
        self.fp2 = np.zeros(0, dtype=int)
        self.face = np.zeros(0, dtype=int)

        # Cells
        self.cp1 = np.zeros(0, dtype=int)
        self.cp2 = np.zeros(0, dtype=int)
        self.cell = np.zeros(0, dtype=int)

        # Region
        self.region = np.zeros(0, dtype=int)

        # ID
        self.ID = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nft):
        # Number
        self.ntot = nft

        # Faces
        self.fp1 = np.zeros(nft, dtype=int)
        self.fp2 = np.zeros(nft, dtype=int)

        # Cells
        self.cp1 = np.zeros(nft, dtype=int)
        self.cp2 = np.zeros(nft, dtype=int)

        # Region
        self.region = np.zeros(nft, dtype=int)

        # ID
        self.ID = np.zeros(nft, dtype=int)

    def InitializeFaceData(self, nftf):
        self.nface = nftf
        self.face = np.zeros(nftf, dtype=int)

    def InitializeCellData(self, nftc):
        self.ncell = nftc
        self.cell = np.zeros(nftc, dtype=int)

    # Face getter
    def GetFace(self, i):
        return self.face[self.fp1[i]:self.fp1[i]+self.fp2[i]]
    
    # Cell getter
    def GetCell(self, i):
        return self.cell[self.cp1[i]:self.cp1[i]+self.cp2[i]]

# Grid topological data
class TopoData: 
    # Definition
    def __init__(self):
        # Number
        self.nxp = 0
        self.nop = 0
        self.nsp = 0
        self.ntp = 0
        self.ndiv = 0
        self.ndivFc = 0
        self.topoID = 0

        # Points
        self.XpointID = np.zeros(0, dtype=int)
        self.OpointID = np.zeros(0, dtype=int)
        self.SpointID = np.zeros(0, dtype=int)
        self.TpointID = np.zeros(0, dtype=int)

        # Point data
        self.isprimaryxp = np.zeros(0, dtype=int)
        self.spointxpID = np.zeros(0, dtype=int)
        self.spointdivID = np.zeros(0, dtype=int)
        self.tpointdivID = np.zeros(0, dtype=int)

        # Divertor face data
        self.divFcP1 = np.zeros(0, dtype=int)
        self.divFcP2 = np.zeros(0, dtype=int)
        self.divFc = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nX, nO, nS, nT, nDiv, nDivFc, topoID):
        # Number
        self.nxp = nX
        self.nop = nO
        self.nsp = nS
        self.ntp = nT
        self.ndiv = nDiv
        self.ndivFc = nDivFc
        self.topoID = topoID

        # Points
        self.XpointID = np.zeros(nX, dtype=int)
        self.OpointID = np.zeros(nO, dtype=int)
        self.SpointID = np.zeros(nS, dtype=int)
        self.TpointID = np.zeros(nT, dtype=int)

        # Point data
        self.isprimaryxp = np.zeros(nX, dtype=int)
        self.spointxpID = np.zeros(nS, dtype=int)
        self.spointdivID = np.zeros(nS, dtype=int)
        self.tpointdivID = np.zeros(nT, dtype=int)

        # Divertor face data
        self.divFcP1 = np.zeros(nDiv, dtype=int)
        self.divFcP2 = np.zeros(nDiv, dtype=int)
        self.divFc = np.zeros(nDivFc, dtype=int)

    # Divertor face getter
    def GetDivFace(self, i):
        return self.divFc[self.divFcP1[i]:self.divFcP1[i]+self.divFcP2[i]]

# Grid
class Grid:
    # Init
    def __init__(self):
        # Fields
        self.vert = Vert()
        self.face = Face() 
        self.cell = Cell()
        self.fs = FluxSurf()
        self.ft = FluxTube()
        self.topodata = TopoData()

    # Interconnection computation
    def ComputeInterconnections(self):
        # Computes additional grid interconnection data, assuming
        # all other necessary fields were read in

        # Face cells
        self.face.nb1[:] = 0 
        self.face.nb2[:] = 0
        for i in np.arange(0, self.cell.ntot, 1):
            # Get cell faces
            tf = self.cell.GetFace(i)

            # Loop over each face
            for j in np.arange(0, len(tf), 1):
                if self.face.nb1[tf[j]-1] == 0:
                    self.face.nb1[tf[j]-1] = i+1 
                elif self.face.nb2[tf[j]-1] == 0:
                    self.face.nb2[tf[j]-1] = i+1 
                else:
                    raise('ComputeInterconnections: face has more than two cell neighbours')

        # Cell neighbours
        cc = 0 
        for i in np.arange(0, self.cell.ntot, 1):
            # Get cell faces
            tf = self.cell.GetFace(i)

            # Loop over each face
            for j in np.arange(0, len(tf), 1):
                if (self.face.nb1[tf[j]-1] == i+1 and self.face.nb2[tf[j]-1] != 0):
                    self.cell.nb[cc] = self.face.nb2[tf[j]-1]
                    self.cell.nbp2[i] = self.cell.nbp2[i] + 1
                    cc = cc + 1
                elif (self.face.nb2[tf[j]-1] == i+1 and self.face.nb1[tf[j]-1] != 0):
                    self.cell.nb[cc] = self.face.nb1[tf[j]-1]
                    cc = cc + 1
                    self.cell.nbp2[i] = self.cell.nbp2[i] + 1

        # Construct pointer
        self.cell.nbp1[0] = 0 # Account for zero-based indexing
        for i in np.arange(1, self.cell.ntot, 1):
            self.cell.nbp1[i] = self.cell.nbp1[i-1] + self.cell.nbp2[i-1]

        # Vertex pointer setup (assumed initialized)
        for i in np.arange(0, self.cell.ntot, 1):
            # Get vertices of this cell
            tv = self.cell.GetVert(i)

            # Update counter
            self.vert.cp2[tv-1] = self.vert.cp2[tv-1] + 1

        for i in np.arange(1, self.vert.ntot, 1):
            self.vert.cp1[i] = self.vert.cp1[i-1] + self.vert.cp2[i-1]

        # Vertex cell list construction
        vcounter = np.zeros(self.vert.ntot, dtype=int)
        self.vert.cell = np.zeros(np.sum(self.vert.cp2), dtype=int)
        for i in np.arange(0, self.cell.ntot, 1):
            #  Get the cell vertices
            tv = self.cell.GetVert(i)

            # Set the vertex cells
            for j in np.arange(0, len(tv), 1):
                ind = self.vert.cp1[tv[j]-1] + vcounter[tv[j]-1]
                self.vert.cell[ind] = i+1 
                vcounter[tv[j]-1] = vcounter[tv[j]-1] + 1

    # Check if point lies in cell domain
    def InCell(self, i, xq, yq):
        # Get cell vertices
        tv = self.cell.GetVert(i)

        # Get coordinates
        vx = self.vert.x[tv-1]
        vy = self.vert.y[tv-1]

        # Build points
        pts = np.zeros((len(vx), 2), dtype=float)
        for j in np.arange(0, len(tv), 1):
            pts[j, 0] = vx[j]
            pts[j, 1] = vy[j]

        # Create a polygon
        pol = geometry.Polygon(pts)
        pt = geometry.Point(xq, yq)
        
        # Check if point in polygon
        inpolyg = pol.contains(pt)

        # Check if point on boundary
        bnd = shapely.boundary(pol)
        onboundary = bnd.contains(pt)
        isincell = inpolyg or onboundary

        return isincell

    # Compute grid metrics
    def ComputeFaceCoordinates(self):
        # Compute face center coordinates
        for i in np.arange(0, self.face.ntot, 1):
            tv = np.array([self.face.v1[i], self.face.v2[i]])
            self.face.x[i] = np.mean(self.vert.x[tv-1])
            self.face.y[i] = np.mean(self.vert.y[tv-1])

    def ComputeCellCoordinates(self):
        # Compute cell center coordinates - only for internal cells, 
        # guard cells should be given or will be zero
        for i in np.arange(0, self.cell.nci, 1): 
            tv = self.cell.GetVert(i)
            self.cell.x[i] = np.mean(self.vert.x[tv-1])
            self.cell.y[i] = np.mean(self.vert.y[tv-1])

    def ComputeCellSurfaceArea(self):
        surfA = np.zeros(self.cell.ntot, dtype=float)
        for i in np.arange(0, self.cell.ntot):
             tv = self.cell.GetVert(i)
             xv = self.vert.x[tv-1]
             yv = self.vert.y[tv-1]
             surfA[i] = Compute2DSurfaceArea(xv, yv)
        return surfA


#----------------------------------------------------------------------#
#                                METRICS                               #
#----------------------------------------------------------------------#

def Compute2DSurfaceArea(x, y):
    # Description
    #------------
    # Compute the (signed) 2D surface area underneath a polygon defined by the 
    # x, y coordinate sequence. We do trapezoidal 

    if (x[0] == x[len(x)-1]) and  (y[0] == y[len(y)-1]):
        pass
    else:
        x = np.append(x, x[0])
        y = np.append(y, y[0])

    dx = x[1:len(x)] - x[0:len(x)-1]
    yf = 0.5*(y[0:len(y)-1] + y[1:len(y)])
    surfA = np.sum(dx*yf)
    return surfA

#----------------------------------------------------------------------#
#                        GENERAL 2D INTERPOLANT                        #
#----------------------------------------------------------------------#

# General 2D unstructured interpolant based on grid
class GridInterpolant2D:
    def __init__(self):
        # Grid description
        self.grid = Grid()

        # Points
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # (normalized) field for curvilinear method evaluated at cell centers
        self.bx = np.zeros(0, dtype=float)
        self.by = np.zeros(0, dtype=float)

        # Values at points
        self.f = np.zeros(0, dtype=float)

        # Gradient at points (doesn't have to be in x, y direction)
        self.dfdx = np.zeros(0, dtype=float)
        self.dfdy = np.zeros(0, dtype=float)

        # Length method
        self.method = ''

    def Initialize(self, grid):
        self.grid = grid 
        self.x = np.zeros(grid.cell.ntot, dtype=float)
        self.y = np.zeros(grid.cell.ntot, dtype=float)
        self.bx = np.zeros(grid.cell.ntot, dtype=float)
        self.by = np.zeros(grid.cell.ntot, dtype=float)
        self.f = np.zeros(grid.cell.ntot, dtype=float) 
        self.dfdx = np.zeros(grid.cell.ntot, dtype=float)
        self.dfdy = np.zeros(grid.cell.ntot, dtype=float)
        self.method = ''

    def Construct(self, grid, cellvalues, method):
        # Construct the interpolant based on the values defined 
        # in the cell centers. 

        # Initialize
        self.Initialize(grid)

        # checks
        if (not (len(cellvalues) == grid.cell.ntot)):
            raise('Construction of GridInterpolant2D: cell values do '
                  + 'not have the same dimension as grid cells')
        
        # Add basic values
        self.f = cellvalues 
        self.method = method 
        self.grid = grid 

        # Compute cell center coordinates and add
        for i in np.arange(0, grid.cell.ntot, 1):
            tv = grid.cell.GetVert(i)
            self.x[i] = np.mean(grid.vert.x[tv-1])
            self.y[i] = np.mean(grid.vert.y[tv-1])

            tbx = -grid.vert.by[tv-1] # actually bx is dpsidx etc 
            tby = grid.vert.bx[tv-1] 

            self.bx[i] = np.mean(tbx) # May need a better way to determine this...
            self.by[i] = np.mean(tby)

        # Normalize
        bn = np.sqrt(self.bx**2 + self.by**2)
        self.bx = self.bx/bn 
        self.by = self.by/bn

        # Compute derived quantities
        match self.method: 
            case 'cartesian':
                # Reconstruct gradient based on cartesian coordinates
                for i in np.arange(0, grid.cell.ntot, 1):
                    # Get cell neighbours (9point)
                    tv = grid.cell.GetVert(i)
                    tc = np.zeros(0, dtype=int)
                    for j in tv:
                        tc = np.append(tc, grid.vert.GetCell(j-1))
                    tc = np.unique(tc)

                    # Compute distances
                    tdx = self.x[tc-1] - self.x[i]
                    tdy = self.y[tc-1] - self.y[i]

                    # Compute rhs
                    b = self.f[i] - self.f[tc-1]

                    # Compute lhs 
                    a = np.zeros((len(tc), 2)) # not square per se
                    for j in np.arange(0, len(tc)):
                        a[j, :] = [tdx[j], tdy[j]]

                    A = np.matmul(np.transpose(a), a)

                    # Solve
                    x = np.linalg.solve(A, np.matmul(np.transpose(a), b))

                    # Add
                    self.dfdx[i] = x[0]
                    self.dfdy[i] = x[1]


                    # Compute coefficients by unweighted least squares
            case 'curvilinear':
                # Reconstruct gradient based on curvilinear coordinates
                for i in np.arange(0, grid.cell.ntot, 1):
                    # Get cell neighbours (9point)
                    tv = grid.cell.GetVert(i)
                    tc = np.zeros(0, dtype=int)
                    for j in tv:
                        tc = np.append(tc, grid.vert.GetCell(j-1))
                    tc = np.unique(tc)

                    # Compute distances
                    tdx = self.x[tc-1] - self.x[i]
                    tdy = self.y[tc-1] - self.y[i]
                    tdtheta = tdx*self.bx[tc-1] + tdy*self.by[tc-1]
                    tdr = -tdx*self.by[tc-1] + tdy*self.bx[tc-1]

                    # Compute rhs
                    b = self.f[i] - self.f[tc-1]

                    # Compute lhs 
                    a = np.zeros((len(tc), 2)) # not square per se
                    for j in np.arange(0, len(tc)):
                        a[j, :] = [tdtheta[j], tdr[j]]

                    A = np.matmul(np.transpose(a), a)

                    # Solve
                    x = np.linalg.solve(A, np.matmul(np.transpose(a), b))

                    # Add
                    self.dfdx[i] = x[0]
                    self.dfdy[i] = x[1]
            case _: 
                raise('ConstructGridInterpolant2D: unknown method')

    def Evaluate(self, xq, yq):
        # Description
        #------------
        # Evaluate the interpolant at xq, yq coordinates and return the
        # values in vq

        # Checks
        if (not (len(xq) == len(yq) )):
            raise('EvaluateGridInterpolant2D: incompatible sizes of query point coordinates')
        
        # Initialize
        vq = np.zeros(len(xq), dtype=float)

        # Check which distance metric to take
        match self.method:

            case 'cartesian': # based on regular x, y coordinates

                # Check for each query point which cell lies closest
                for i in np.arange(0, len(xq), 1):
                    # Unpack
                    txq = xq[i]
                    tyq = yq[i]

                    # Compute distance
                    dist = np.sqrt((self.x - txq)**2 + (self.y - tyq)**2) 
                    mindistind = np.argmin(dist)

                    # Check if point in polygon, if not: set dist to nan and
                    # keep looking
                    isincell = self.grid.InCell(mindistind, txq, tyq)
                    dist[mindistind] = np.inf 
                    while ((not isincell) and (not all(dist == np.inf))):
                        mindistind = np.argmin(dist)
                        isincell = self.grid.InCell(mindistind, txq, tyq)
                        dist[mindistind] = np.inf 
                    
                    if (all(dist == np.inf) and (not isincell)):
                        vq[i] = np.nan # set to nan and return
                    else:

                        dx = txq - self.x[mindistind]
                        dy = tyq - self.y[mindistind] 

                        # Compute value
                        vq[i] = self.f[mindistind] + self.dfdx[mindistind]*dx + self.dfdy[mindistind]*dy

            case 'curvilinear': # based on radial/poloidal distance

                # Check for each query point which cell lies closest
                for i in np.arange(0, len(xq), 1):
                    # Unpack
                    txq = xq[i]
                    tyq = yq[i]

                    # Compute distance
                    dist = np.sqrt((self.x - txq)**2 + (self.y - tyq)**2) 
                    mindistind = np.argmin(dist)
                    dx = txq - self.x[mindistind]
                    dy = tyq - self.y[mindistind] 

                    # Check if point in polygon, if not: set dist to nan and
                    # keep looking
                    isincell = self.grid.InCell(mindistind, txq, tyq)
                    dist[mindistind] = np.inf 
                    while ((not isincell) and (not all(dist == np.inf))):
                        mindistind = np.argmin(dist)
                        isincell = self.grid.InCell(mindistind, txq, tyq)
                        dist[mindistind] = np.inf 

                    if (all(dist == np.inf) and (not isincell)):
                        vq[i] = np.nan # set to nan and return
                    else:
                        # Compute poloidal (parallel) and radial (perpendicular) distance
                        dtheta = dx*self.bx[mindistind] + dy*self.by[mindistind]
                        dr = -dx*self.by[mindistind] + dy*self.bx[mindistind]

                        # Compute value
                        vq[i] = self.f[mindistind] + self.dfdx[mindistind]*dtheta + self.dfdy[mindistind]*dr
                    
            case _: 
                raise('EvaluateGridInterpolant2D: unknown method')

        # Return 
        return vq
