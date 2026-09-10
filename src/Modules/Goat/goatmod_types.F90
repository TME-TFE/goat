!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all the type definitions used in the goat. 
! See the option setting routines in goatmod_userinput
! for an explanation of the different fields and options. Additionally,
! precision types are set here that can (should) be used in all 
! subroutines. 

! Note: default accessibility is set to public, since all these types 
! and corresponding functions are free to be used. 

module goatmod_types

    ! Initialize
    !============
    use goatmod_userinput
    use mod_precision
    use mod_polygon
    use mod_inputfileparser
    use mod_plotter
    use mod_dynamicarrays
    use mod_sort, only: Unique
    use Interpolant
    use PolygonShapeFunction
    use PolygonLevelsetFunction2D

    ! The usual
    implicit none
    save
    public 

    ! All private types
    
    ! All private functions

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                               Grid                               !
    !------------------------------------------------------------------!

    ! Main grid substructures
    !========================
    ! Vertex structure
    type VertexUDT

        ! Description
        !============
        ! Fields:
        ! - x, y            : coordinates [m]
        ! - BV              : logical index to indicate if vertex lies
        !                   on boundary
        ! - fieldlineID     : ID of the field line the vertex lies on
        ! - ntot            : total number of vertices  
        ! - faceP           : ntot-by-2 array containing in the first 
        !                   column the starting index and in the second
        !                   column the number of faces of the cell 
        !                   (for querying the faces of the cell
        !                   stored in vert%facelist)            
        ! - facelist        : list of vertex faces, to be queried as 
        !                   vert%facelist(vert%faceP(i,1):
        !                   vert%faceP(i,1)+vert%faceP(i,2)-1) 
        ! - nfacelist       : length of facelist   
        ! - cellP           : similar to faceP, but for vertex cells
        ! - celllist        : similar to facelist, but for vertex cells
        ! - ncellist        : length of celllist    
        ! - neiglist        : list of neighbouring vertices (nneiglist-
        !                   by-1)
        ! - nneiglist       : dimension of neiglist (scalar)
        ! - neigP           : ntot-by-2 array, analogous to faceP and 
        !                   cellP, but for neiglist.
        ! - psi             : psi values at vertex locations    
        ! - bx, by          : magnetic field vector at vertex locations  
        ! - ffbz            : ???

        ! Coordinates
        real(R8), allocatable               :: x(:), y(:)

        ! Logicals and indices
        logical, allocatable                :: BV(:)
        integer(I8), allocatable            :: fieldlineID(:) 
        integer(I8)                         :: ntot = 0

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: face(:)
        integer(I8)                         :: nface = 0

        integer(I8), allocatable            :: cellP(:,:)
        integer(I8), allocatable            :: cell(:)
        integer(I8)                         :: ncell = 0

        integer(I8), allocatable            :: neigP(:,:)
        integer(I8), allocatable            :: neig(:)
        integer(I8)                         :: nneig = 0

        ! Other data
        real(R8), allocatable, dimension(:) :: bx, by, psi, ffbz 

    end type

    ! Face structure
    type FaceUDT

        ! Description
        !============
        ! Fields:
        ! - ntot            : total number of faces
        ! - vert            : set of (two) vertices belonging to that 
        !                   face
        ! - cell            : list of cell neighbours (max 2)
        ! - cellP           : pointer for list of cell neighbours
        ! - neig            : cell neighbours of face
        ! - BF              : logical index that is true if the face
        !                   is a boundary face
        ! - label           : face labels
        ! - TMfacelabel     : face label corresponding to topomesh ID (only when GG is used)
        ! - reg             : face regions
        ! - aligned         : integer that is 1 if the face is aligned


        ! Logicals and indices
        integer(I8), allocatable            :: vert(:,:)

        integer(I8), allocatable            :: cell(:), label(:), &
            reg(:), aligned(:), TMfacelabel(:)
        integer(I8), allocatable            :: cellP(:, :)
        integer(I8)                         :: ncell = 0

        integer(I8)                         :: ntot = 0
        logical, allocatable                :: BF(:)

    end type


    ! Cell structure
    type CellUDT

        ! Description
        !============
        ! Fields:   
        ! - ntot            : total number of cells
        ! - vertP           : ntot-by-2 array containing in the first 
        !                   column the starting index and in the second
        !                   column the number of vertices of the cell 
        !                   (for querying the vertices of the cell
        !                   stored in cells%vert)            
        ! - vert            : list of cell vertices, to be queried as 
        !                   cells%vert(cells%vertP(i,1):
        !                   cells%vertP(i,1)+cells%vertP(i,2)-1) 
        ! - nvert          : length of vert
        ! - faceP           : similar to vertP, but for faces
        ! - face            : similar to vert, but for faces
        ! - nface          : similar to nvert, but for faces
        ! - GC              : ntot-by-1 logical vector indicating if the
        !                   cell is a guard cell


        ! Logicals and indices
        integer(I8), allocatable            :: vertP(:,:)
        integer(I8), allocatable            :: vert(:)
        integer(I8)                         :: nvert = 0

        integer(I8), allocatable            :: faceP(:,:)
        integer(I8), allocatable            :: face(:)
        integer(I8)                         :: nface = 0

        logical, allocatable                :: GC(:)

        integer(I8)                         :: ntot = 0, ngc
        
        real(R8), allocatable, dimension(:) :: psi, bp, bt, x, y
        integer(I8), allocatable, dimension(:)  :: cflags, reg, &
            ft
    end type


    ! Boundary structure
    type BndUDT

        ! Description
        !============
        ! Fields:   
        ! - nvert:                      total number of vertices
        ! - vert:                       all vertices belonging to that 
        !                               boundary
        ! - nface:                      total number of faces 
        ! - face:                       all face indices 
        ! - ID:                         boundary ID (see below) that 
        !                               identifies the boundary type.
        ! 
        ! See also the interface routine InterfaceBoundaryMapping in 
        ! gdmod_interfaces on how the boundary IDs are mapped given 
        ! input from the grid generator. 

        ! Note: the amount of vertices is per definition the amount of
        ! faces minus one, also for closed boundaries. In the latter 
        ! case, the begin and end vertex index are the same. This 
        ! assumes that the boundaries should be polygons that either 
        ! close perfectly on themselves, or are simple polygons (i.e. 
        ! no polygon vertex should have more than three faces it belongs
        ! to).

        ! Boundary ID meaning
        !====================
        ! GRID DEFORMATION              PHYSICAL MEANING
        !
        ! 1                             target plate (inner)
        ! 2                             target plate (outer)
        ! 3                             private flux
        ! 4                             core boundary
        ! 5                             outermost flux surf.


        ! Boundary vertices
        integer(I8)                         :: nvert = 0 ! simply nfaces-1 actually
        integer(I8), allocatable            :: vert(:)   

        ! Boundary faces
        integer(I8)                         :: nface = 0
        integer(I8), allocatable            :: face(:)             
        
        ! Boundary ID
        integer(I8)                         :: ID

    end type

    ! Additional grid data structures
    !================================
    ! Used to collect derived grid data, e.g. which link the grid cells 
    ! flux tubes. Other output data from the grid generator that does 
    ! not fit within the classical cell/faces/vert structures can be 
    ! added as well through user defined types. 

    ! Flux data
    type FluxDataUDT

        ! Description
        !============
        ! This data type contains all information regarding flux tubes 
        ! and flux surfaces, and the cells associated with these. 
        ! Fields:
        !
        ! - nFt:            : total number of flux tubes (scalar, to be 
        !                   read in)
        ! - nFs:            : total number of flux surfaces (scalar, to be
        !                   read in)
        ! - fluxtubecellsP  : nFt-by-2 array where the first index 
        !                   is the start index in the fluxtubecells 
        !                   array, and the second the amount of cells of
        !                   the flux tube.
        ! - fluxtubecells   : nCv(number of cells)-by-1 array containing 
        !                   the cell numbers that correspond to flux 
        !                   tubes. 
        ! - fluxtubefacesP  : nFt-by-2 array where the first index 
        !                   is the start index in the fluxtubefaces 
        !                   array, and the second the amount of faces of
        !                   the flux tube.
        ! - fluxtubefsIDs   : flux surface IDs that bound the flux tube
        !                   (nFt-by-2)
        ! - fluxtubefaces   : nFv(number of faces)-by-1 array containing 
        !                   the face numbers that correspond to flux 
        !                   tubes. 
        ! - fluxsurfacefacesP  : nFs-by-2 array where the first index 
        !                   is the start index in the fluxsurfacefaces
        !                   array, and the second the amount of faces of
        !                   the flux surface.
        ! - fluxsurfacefaces   : nFv(number of faces)-by-1 array containing 
        !                   the face numbers that correspond to flux 
        !                   surfaces. 
        ! - fluxsurfaceID   : nv(number of vertices)-by-1 array 
        !                   containing the flux tube ID of each vertex

        ! Logicals and indices
        integer(I8)                         :: nFt = 0
        integer(I8)                         :: nFs = 0

        ! Arrays, flux tube data
        integer(I8), allocatable            :: fluxtubecellsP(:,:)
        integer(I8), allocatable            :: fluxtubecells(:)
        integer(I8), allocatable            :: fluxtubefacesP(:,:)
        integer(I8), allocatable            :: fluxtubefaces(:)
        integer(I8), allocatable            :: fluxtubefsIDs(:, :), &
            fluxtuberegID(:)
        logical, allocatable                :: isclosedft(:)

        ! Arrays, flux surface data
        integer(I8), allocatable            :: fluxsurfacefacesP(:,:)
        integer(I8), allocatable            :: fluxsurfacefaces(:)
        integer(I8), allocatable            :: fluxsurfaceID(:)
        integer(I8), allocatable            :: fluxsurfaceneig(:), fluxsurfaceneigP(:, :)
        real(R8), allocatable               :: fluxsurfacepsi(:)

        integer(I8), allocatable            :: fluxsurfacevertsP(:,:)
        integer(I8), allocatable            :: fluxsurfaceverts(:)


    end type

    ! Structured grid data (to be removed in the future)
    type StructuredGridDataUDT

        ! Description
        !============
        ! Data type to collect all data related to a possibly initial 
        ! structured grid that served as basis for the current 
        ! unstructured grid data. Only saved for backward compatibility 
        ! reasons. Hopefully deleted in the future. 
        ! Fields:
        !
        ! - isClassicalGrid         : integer indicating whether the
        !                           whether the grid was a classic 
        !                           structured grid
        ! - nx, ny                  : dimensions of the original 
        !                           structured grid
        ! - nncut                   : number of cuts

        ! Logicals & scalars
        integer(I8)             :: nx, ny, nncut
        integer(I4)             :: isClassicalGrid

    end type

    ! Goat grid generator data
    type GoatGGDataUDT

        ! Description
        !============
        ! Type that contains any useful grid generator data from goat
        ! that does not fall in the scope of other GridDataUDT 
        ! fields. For example this can contain topomesh data of grid 
        ! faces/vertices (TMfacetype/TMverttype), which faces are 
        ! boundary layer faces, etc. 

        ! Current fields: 
        ! - TMfacetype:     topological mesh face type (0 if not on topomesh face)
        ! - TMfacevert:     same, but for vertices
        ! - BLind:          boundary layer number (0 if not a boundary 
        !                   layer face, otherwise 1, 2, ... where 1 is closest to vessel)    
        ! - facelabelsGG, facelabelsGD:     mapping from GG to GD face labels

        ! Fields
        integer(I8), allocatable, dimension(:)  :: TMfacetype, TMverttype, &   
            BLind, facelabelsGG, facelabelsGD

    end type

    ! Grid data 
    type GridDataUDT

        ! Description
        !===========
        ! Data type to collect all other grid related data which is not 
        ! linked to the core data, i.e. the cells, vertices, and faces.
        ! This can be additional info from the grid generator on flux 
        ! tubes, for example, and interconnection data thereof (e.g. 
        ! the cells that belong to a certain flux tube). 
        ! Fields:
        !
        ! - fluxdata            : UDT with all flux data such as flux
        !                       flux tube data, flux surfaces, ... 
        ! - goatggdata          : additional goat grid generator data
        !                       (if available)
        ! - sglegacy            : data from legacy structured grids
        ! - OMPcell, OMPface    : cells and faces belonging to outer mid
        !                       plane
        ! - IMPcell, IMPface    : same, but inner mid plane
        ! - OMPr, OMPz          : points defining line segment of OMP
        ! - IMPr, IMPz          : same but for inner mid plane
        ! - topoflag            : flag indicating the topological mesh 
        !                       type (see also mod_definitions)
        ! - xpointID            : array containing all X-point vertex IDs
        ! - isprimaryxp         : integer array which is 1 if the 
        !                       x-point is a primary x-point (i.e. it 
        !                       lies on a separatrix that is connected
        !                       directly to a core region), otherwise 0
        ! - nxp                 : number of x-points
        ! - spointID            : array containing all strike point vertex IDs
        ! - nsp                 : number of strike points
        ! - opointID            : O point IDs (in the grid)
        ! - nop                 : number of o points
        ! - tpointID            : Tangency point IDs (but here tangency 
        !                       points are only those tangency points 
        !                       that have a closed contour in the grid!)
        ! - ntp                 : number of t points
        ! - spointxpID          : x-point vertex number of the x-point 
        !                       the strike point belongs to
        ! - ndiv                : number of divertor targets/plates
        ! - spointdivID         : list of targets that belong to strike 
        !                       points (one target each)
        ! - tpointdivID         : same as above but for tangency points
        ! - ndivFc              : total number of divertor faces
        ! - divFc               : list of divertor faces
        ! - divFcP              : pointer for the list of divertor faces
        ! - sepID               : provides the fieldlineID of the separatrices
        ! - nsep                : number of separatrices    

        ! Flux data
        type(FluxDataUDT)           :: fluxdata

        ! Legacy data of structured grid
        type(StructuredGridDataUDT) :: sglegacy

        ! Goat grid generator data
        logical                     :: hasGoatGGData
        type(GoatGGDataUDT)         :: goatggdata 

        ! Topological mesh type
        integer(I8)                             :: topoflag
        
        ! OMP & IMP
        integer(I8), allocatable, dimension(:)  :: OMPcell, OMPface, &
            IMPcell, IMPface
        integer(I8)                             :: nOMPcell, nOMPface, &
            nIMPcell, nIMPface
        real(R8), dimension(1:2)                :: OMPr, OMPz, IMPr, &
            IMPz

        ! X-point(s), strike points, o points, ...
        integer(I8), allocatable, dimension(:)  :: xpointID, &
            spointID, opointID, tpointID, isprimaryxp, spointxpID, &
            divFc, spointdivID, tpointdivID
        integer(I8), allocatable, dimension(:, :)   :: divFcP
        integer(I8)                             :: nxp, nsp, nop, ntp, &
            ndiv, ndivFc

        ! Separatrices
        integer(I8), allocatable                :: sepID(:)
        integer(I8)                             :: nsep

    end type

    ! Main grid structure
    !====================
    type GridUDT

        ! Description
        !============
        ! Data structure containing all the grid data substructures. 
        ! Fields:
        !
        ! - vert            : see type definition for description
        ! - face           
        ! - cell
        ! - data
        ! - bnd
        !
        ! Note: the bnd substructure has to be allocated separately

        ! Vertices
        type(VertexUDT)                     :: vert

        ! Faces
        type(FaceUDT)                       :: face

        ! Cells
        type(CellUDT)                       :: cell

        ! Additional data
        type(GridDataUDT)                   :: data

        ! Boundaries
        type(BndUDT), allocatable           :: bnd(:)

    contains

        procedure  :: WriteData         => WriteGridData

    end type

    !------------------------------------------------------------------!
    !                               Numerics                           !
    !------------------------------------------------------------------!

    ! Interpolation
    !==============

    !------------------------------------------------------------------!
    !                           Magnetic field                         !
    !------------------------------------------------------------------!

    ! Magnetic field
    !===============
    type MagneticFieldUDT

        ! Description
        !============
        ! Data structure containing all the magnetic field data. The 
        ! following fields should be present:
        !
        ! - interpolant:    see type definition for explanation
        ! - nR:             number of sample points in the R-direction
        ! - nZ:             number of sample points in the Z-direction
        ! - R:              nR-by-1 array with R-coordinates
        ! - Z:              nZ-by-1 array with Z-coordinates
        ! - Psi:            nR-by-nZ array with magnetic flux values
        ! - RBtor:          product of major radius and toroidal field (constant)

        ! Coordinates
        integer(I8)                 :: nR = 0 
        integer(I8)                 :: nZ = 0
        real(R8), allocatable       :: R(:)
        real(R8), allocatable       :: Z(:)
        real(R8), allocatable       :: Psi(:,:)
        real(R8)                    :: RBtor 
        
        ! Interpolant
        type(StructuredInterpolant2DUDT)    :: interp

    end type

    !------------------------------------------------------------------!
    !                          DivGeo coupling                         !
    !------------------------------------------------------------------!
    type DivGeoDataUDT

        ! Description
        !============
        ! DivGeo data structure to manipulate output obtained from .dgo
        ! file. This typically includes elements that form structure or
        ! vessel polygons etc. From this file, inputs may be retrieved
        ! if it is available. Note that this object does not read in 
        ! magnetic field data, as this is defined in the .equ and should
        ! be read in through the magnetic field routines.

        ! Currently, the following data is read in :
        ! - Element point coordinates pairs P1, P2 (x, y, z, assumed in 
        !   csv-like format - but we don't store z)
        ! - Element face labels (fcLbl)
        ! - Vessel elements (elvessel)

        ! We indicate whether data is available using logicals
        
        integer(I8)                             :: nel, nv
        integer(I8), allocatable, dimension(:)  :: elID, elfcLbl, elv1, &
            elv2, elvessel, trimark
        real(R8), allocatable, dimension(:)     :: elvx, elvy

        logical     :: hasEl, hasElFcLbl, hasElVessel, hasElVoid

    contains 

        ! Allocation
        procedure :: Allocate   => AllocateDGData
        procedure :: Deallocate => DeallocateDGData

        ! Initialization
        procedure :: Initialize => InitializeDGData

        ! Reader
        procedure :: Read       => ReadDGFile 

        ! Structure extraction 
        procedure :: ExtractStructures  => ExtractDGStructures

        ! Vessel structure extraction
        procedure :: ExtractVesselStructures => ExtractDGVesselStructures

        ! Triangulation structure extraction
        procedure :: ExtractTriangulationStructures  => ExtractDGTriangulationStructures

    end type

    !------------------------------------------------------------------!
    !                            Environment                           !
    !------------------------------------------------------------------!
    ! Vessel structures
    !==================
    type VesselStructureUDT

        ! Description
        !============
        ! User defined type for the different substructures of the 
        ! vessel. These are actually polygons (so with x, y coordinates)
        ! that are used to build up the vessel structure afterwards. 
        ! The following fields are present:
        !
        ! - np:         number of points in the polygon
        ! - x, y:       polygon coordinates
        ! - isclosed:   logical to indicate if the polygon should close
        !               upon itself
        ! - id:         identifier (integer number)
        ! - label:      (optional) additional label. Zero by default, 
        !               unless otherwise specified. 

        ! Coordinates
        integer(I8)                         :: np
        real(R8), allocatable               :: x(:), y(:)

        ! Logicals
        logical                             :: isclosed

        ! ID
        integer(I4)                         :: ID, label

    end type

    ! State 
    type StateUDT

        ! Description
        !============
        ! User defined type containing simulation data (SOLPS-ITER) on the given grid
        ! The following field are present:
        !
        ! - na:     ion density
        ! - ne:     electron density
        ! - ua:     parallel velocity
        ! - uadia:  diamagnetic velocity
        ! - te:     electron density
        ! - ti:     ion temperature
        ! - tn:     neutral temperature
        ! - po:     potential
        ! - kt:     ??
        ! - zt:     ??

        ! - resco:  residual of continuity equation
        ! - reshe:  residual of electron temperature equation
        ! - reshi:  residual of ion temperature equation
        ! - reshn:  residual of neutral temperature equation
        ! - resmo:  residual of ion momentum equation
        ! - resmt:  residual of total momentum equation
        ! - respo:  residual of potential equation
        ! - reskt:  ??
        ! - reszt:  ??

        ! State 
        integer(I8)                         :: nc, nf, ns
        real(R8), allocatable               :: na(:,:), ne(:), ua(:,:), uadia(:,:,:), & 
                                                te(:), ti(:), tn(:), po(:), kt(:), zt(:)

        ! Residual
        real(R8), allocatable               :: resco (:,:), reshe(:), reshi(:), reshn(:), &
                                                resmo(:,:), resmt(:), respo(:), reskt(:), reszt(:)

    end type

    ! Vessel
    !=======
    type VesselUDT
        
        ! Description
        !============
        ! User defined structure to keep track of all vessel related 
        ! data that is required for the grid optimization. Currently, 
        ! this includes the following structures and fields:
        !
        ! - np:         number of points in the boundary polygon
        ! - x, y:       x, y-coordinates (np-by-1 array)
        ! - TPind:      integer array (np-by-1) indicating on which TP
        !               a node is located
        ! - ntp:        number of all target plate indices
        ! - allTPind:   a list (ntp-by-1) of all target plate indices
        ! - geom:       the vessel geometry being considered (character)
        ! - nstructures:number of sub-structures of the vessel
        ! - structures : nstructures-by-1 array of structures
        ! - targetpolygons:     closed polygon representation of each 
        !                       vessel target plate
        ! - plfvessel:  polygon levelset representation of the closed 
        !               vessel polygon (arbitrary method of approximation)
        ! - exactplfvessel: exact representation of the closed vessel 
        !                   polygon, possibly with refinement of this polygon
        !                   (<...>_noref is the same, but with no refinement, so 
        !                   vertices are original structure polygon vertices except 
        !                   at intersections of structures)

        ! Target plates
        integer(I4)                         :: ntp = 0
        integer(I4), allocatable            :: allTPind(:)
        integer(I4), allocatable            :: TPind(:)

        ! Structures
        integer(I4)                         :: nstructures = 0
        type(VesselStructureUDT), allocatable       :: structures(:)
        type(PolygonSetUDT)                 :: polygonset, targetps 
        type(PolygonUDT), allocatable       :: targetpolygons(:)
        
        ! Shape function
        class(PolygonShapeFunctionUDT), allocatable     :: psf
        class(PolygonLevelsetFunction2DUDT), allocatable   :: plfvessel, &
            plftarget
        type(PolygonLevelsetFunction2DClosedExactUDT)   :: exactplfvessel, &
            exactplfvessel_noref

    contains 

        ! Update vessel description using coordinates
        procedure :: UpdateVesselCoordinates

        ! Vessel coordinate getter
        procedure :: GetVesselCoordinates

        ! Vessel vertex pairs getter
        procedure :: GetVesselVertexPairs

        ! Vessel structure label getters
        procedure :: GetVesselStructureLabelsOnPoints, &
            GetVesselStructureLabelsOnEdges

        ! Vessel polygonset writer (structure.dat format)
        procedure :: WriteVesselPolygon

    end type

    ! Environment
    !============
    type EnvironmentUDT

        ! Description
        !============
        ! Overarching type that stores all other structures etc which 
        ! may be needed for grid optimization, and which are not 
        ! related to the grid or the magnetic field. Currently, 
        ! the vessel structure is stored here that is used to construct
        ! the plasma grid (named 'vessel') and the vessel structure
        ! that should be used to determine the triangulation 
        ! (triangulationvessel). It is assumed that the triangulation 
        ! vessel structure is compatible with the plasma edge grid
        ! vessel structure (they may be the same, or the triangulation
        ! vessel structure may hold additional void regions and exclude
        ! vessel structures with negative label, which are assumed to be
        ! 'fake' vessel structures used to bound the plasma grid)

        ! Note: the routine to set up the vessel is currently a 
        ! standalone routine. Should we include it here as a 
        ! method of the vessel structure?

        type(VesselUDT)                 :: vessel
        type(VesselUDT)                 :: triangulationvessel
        type(StateUDT)                  :: SOLPSstate    

    contains

    end type
    
    !==================================================================!
    !                                                                  !
    !                             INTERFACES                           !
    !                                                                  !
    !==================================================================!

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    contains 


    !------------------------------------------------------------------!
    !                               Grid                               !
    !------------------------------------------------------------------!

    ! Allocation 
    !===========
    ! Main grid structure
    subroutine AllocateGrid(grid)

        ! Description
        !============
        ! Allocate the fields in the vertex, faces, and cell structures
        ! of the grid. 
        !
        ! Note: the bnd structure has to be allocated separately, since 
        ! the size of this structure may vary, depending on the rest of
        ! the grid. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)       :: grid

        ! Return if already allocated
        if (allocated(grid%vert%x)) return

        ! Allocate vertex data
        call AllocateVertices(grid%vert)
        
        ! Allocate face data
        call AllocateFaces(grid%face)

        ! Allocate cell data
        call AllocateCells(grid%cell)

        ! Allocate other data structures
        call AllocateGridData(grid%data,grid)


    end subroutine

    ! Constructor
    subroutine ConstructGrid(grid, gridoptions)
    
        ! Description
        !============
        ! This function further processes the grid and constructs additional 
        ! data needed for most GOAT drivers. 
    
        ! Initialize
        !===========
    
        ! The usual
        implicit none
    
        ! Arguments
        type(GridUDT)                           :: grid
        type(GridOptionsUDT)                    :: gridoptions
    
        ! Construct grid
        !===============
        select case (gridoptions%type)
    
        case ('plasma')
    
            ! Default plasma edge grid, check the provided format
            select case (gridoptions%readmeth)
    
            case ('b2fgmtry_us')
    
                ! Extract the grid data structures
                call ExtractGridData(grid, 'b2fgmtry_us', gridoptions)
    
            case ('traduitb2us')
    
                ! Extract the grid data structures 
                call ExtractGridData(grid, 'traduitb2us', gridoptions)
    
            case default 
    
                call gdErrorHandler('unknown plasma grid input type')
    
            end select
    
    
        case default
    
            call gdErrorHandler('unknown grid type')
    
        end select
    
        ! Compute grid interconnections
        !==============================
        call ComputeGridInterconnections(grid)
    
    end subroutine

    ! Readers
    subroutine ReadTraduitUS(grid, filepath)

        ! Description
        !============
        ! Read in the necessary grid data from a traduit.out.b2us file.
    
    
        ! Initialize
        !===========
        ! Declare modules
        use mod_readwrite
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Main variables
        type(GridUDT)               :: grid
    
        integer                     :: filespec
        integer(I8)                 :: idum(0:9), idum2(1)    
        integer(I8)                 :: nc, nf, nv, nfsFc, nftCv, nftFc
        character(*)               :: filepath
    
        logical                     :: reachedeof, readTopologicalData 
    
        character(:), allocatable   :: chardummy   ! dummy array
        integer(I8), allocatable    :: cdummy(:,:) ! dummy array
        integer(I8), allocatable    :: fdummy(:,:) ! dummy array
        integer(I8), allocatable    :: vdummy(:,:) ! dummy array
        integer(I8), allocatable    :: ftdummy(:) ! dummy array
    
        real(R8), allocatable       :: cdummyr(:,:) ! dummy array
        real(R8), allocatable       :: fdummyr(:,:) ! dummy array
        real(R8), allocatable       :: vdummyr(:,:) ! dummy array
        real(R8), allocatable       :: fsdummyr(:) ! dummy array
        real(R8), allocatable       :: facelistdummy(:) ! dummy array
        integer(I8), allocatable    :: ftCvdum(:), ftFcdum(:), fsFcdum(:)
    
        real(R8), allocatable       :: vdata(:, :), cdata(:, :), &
            fsdata(:, :)
            
        integer(I8), allocatable    :: vlist(:), clist(:), &
            flist(:), ftlist(:), fslist(:), cdatai1(:, :), &
            fdatai(:, :), ftdatai(:, :), fsdatai(:, :), cdatai2(:, :)
    
        ! Loop
        integer(I8)                 :: i
    
        ! Data
        data filespec /60/
        
        ! Read grid dimensions & allocate
        !================================
        ! Open the file
        print *, 'reading grid in traduit format from file: ' // filepath
        open(unit = filespec, file = filepath)
        rewind(filespec)
    
        ! First, read the header with the version
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then
            call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
        end if 
    
        ! Check the version to determine what to read in 
        readTopologicalData = .false. 
        if (chardummy(8:17) >= '03.002.002') then
            ! Topological data should be present
            readTopologicalData = .true.
        end if 
        ! call cfverr(filespec,b2fgmtryversion)
    
        ! Primary array dimensions
        call ReadArray(filespec, 6, idum, 'nCi,nFc,nVx,nCg,nFs,nFt')
        nc = idum(0) ! note: only reading in actual cells, no guard cells
        nf = idum(1)
        nv = idum(2)
    
        ! Add to grid
        grid%cell%ntot         = nc
        grid%face%ntot         = nf
        grid%vert%ntot          = nv
        grid%data%fluxdata%nFs  = idum(4)
        grid%data%fluxdata%nFt  = idum(5)
    
        ! Secondary array dimensions
        call ReadArray (filespec,5,idum,'nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc')
    
        ! Add to grid
        grid%cell%nvert = idum(0)
        grid%cell%nface = idum(1)
        grid%vert%ncell  = idum(3)
        grid%vert%nface  = idum(4)
    
        ! Initialize topological data
        if (allocated(grid%data%isprimaryxp)) deallocate(grid%data%isprimaryxp)
        if (allocated(grid%data%xpointID)) deallocate(grid%data%xpointID)
        if (allocated(grid%data%opointID)) deallocate(grid%data%opointID)
        if (allocated(grid%data%spointID)) deallocate(grid%data%spointID)
        if (allocated(grid%data%divFcP)) deallocate(grid%data%divFcP)
        if (allocated(grid%data%divFc)) deallocate(grid%data%divFc)
        if (allocated(grid%data%spointdivID)) deallocate(grid%data%spointdivID)
        if (allocated(grid%data%tpointdivID)) deallocate(grid%data%tpointdivID)

        ! Read data for structured grid remapping (to be deleted in future)
        call ReadArray (filespec,1,idum2,'isClassicalGrid') 
        grid%data%sglegacy%isClassicalGrid = int(idum2(1), I4) ! cast to I4 type
        if (grid%data%sglegacy%isClassicalGrid == 1) then 
            call ReadArray (filespec,3,idum,'nx,ny,nncut') ! this seems to be wrongly formatted for now - to be checked in the future
            grid%data%sglegacy%nx = idum(0)
            grid%data%sglegacy%ny = idum(1)
        end if
        
        ! Read topological data
        if (readTopologicalData) then 
            ! Read topological mesh flag
            call ReadArray(filespec, 1, idum, 'topoflag')
            grid%data%topoflag = idum(0)
    
            ! Read number of topological points
            call ReadArray(filespec, 6, idum, 'nX,nO,nS,nT,nDiv,nDivFc')
            grid%data%nxp = idum(0)
            grid%data%nop = idum(1)
            grid%data%nsp = idum(2)
            grid%data%ntp = idum(3)
            grid%data%ndiv = idum(4)
            grid%data%ndivFc = idum(5)
            grid%data%nsep = 0
    
            ! Allocate
            allocate(grid%data%xpointID(idum(0)), grid%data%opointID(idum(1)), &
                grid%data%spointID(idum(2)), grid%data%tpointID(idum(3)), &
                grid%data%isprimaryxp(idum(0)), &
                grid%data%divFcP(grid%data%ndiv, 2), grid%data%divFc(grid%data%ndivFc), &
                grid%data%spointdivID(grid%data%nsp), grid%data%tpointdivID(grid%data%ntp), &
                grid%data%spointxpID(grid%data%nsp), grid%data%sepID(grid%data%nsep))
    
            ! Read X-point data
            call ReadSingleLine(filespec, chardummy, reachedeof) ! header
            do i = 1, grid%data%nxp  
                ! Read 
                read(filespec, *) grid%data%xpointID(i), grid%data%isprimaryxp(i), &
                    idum(0)
            end do
    
            ! Read O-point data
            call ReadSingleLine(filespec, chardummy, reachedeof) ! header
            do i = 1, grid%data%nop  
                ! Read 
                read(filespec, *) grid%data%opointID(i)
            end do
    
            ! Read strike point data
            call ReadSingleLine(filespec, chardummy, reachedeof) ! header
            do i = 1, grid%data%nsp  
                ! Read 
                read(filespec, *) grid%data%spointID(i), grid%data%spointxpID(i), &
                    grid%data%spointdivID(i), idum(0)
            end do
    
            ! Read tangency point data
            call ReadSingleLine(filespec, chardummy, reachedeof) ! header
            do i = 1, grid%data%ntp  
                ! Read 
                read(filespec, *) grid%data%tpointID(i), grid%data%tpointdivID(i), &
                    idum(0)
            end do
    
            ! Read divertor data
            call ReadSingleLine(filespec, chardummy, reachedeof) ! header
            do i = 1, grid%data%ndiv  
                ! Read 
                read(filespec, *) idum(0), grid%data%divFcP(i, 1), &
                    grid%data%divFcP(i, 2)
            end do
    
            ! Read divertor face data
            call ReadArray(filespec, grid%data%ndivFc, grid%data%divFc, 'divFc')
    
        else 
            ! Initialize to zero
            grid%data%topoflag = 0
            grid%data%nxp = 0
            grid%data%nop = 0
            grid%data%nsp = 0
            grid%data%ntp = 0
            grid%data%ndiv = 0
            grid%data%ndivFc = 0
            grid%data%nsep = 0
    
            ! Allocate
            allocate(grid%data%xpointID(0), grid%data%opointID(0), &
                grid%data%spointID(0), grid%data%isprimaryxp(0), &
                grid%data%divFcP(0, 2), grid%data%divFc(0), &
                grid%data%spointdivID(0), grid%data%tpointdivID(0), &
                grid%data%sepID(0))
    
        end if 
    
        ! Allocate grid
        call AllocateGrid(grid)
    
        ! Read grid interconnection data
        !===============================
        ! Initialize
        !-----------
        ! Allocate dumies
        allocate(cdummy(nc,4))
        allocate(fdummy(nf,4))
        allocate(vdummy(nv,4))
        allocate(cdummyr(nc,4))
        allocate(fdummyr(nf,4))
        allocate(vdummyr(nv,4))
        allocate(fsdummyr(grid%data%fluxdata%nFs))
        allocate(facelistdummy(grid%cell%nface))
        allocate(ftdummy(grid%data%fluxdata%nFt))
        
        ! Attempt to read in OMP/IMP data 
        !--------------------------------
        ! OMP
        call ReadUntilFound(filespec, 'OMPr', reachedeof)
        if (reachedeof) then 
            ! Not found, issue warning and rewind
            print *, 'ReadTraduitUS: could not find field "OMPr", setting to zero...'
            rewind(filespec)
        else 
            ! Found, read in coordinates (assumed two)
            backspace(filespec)
            call ReadArray(filespec, 2, grid%data%OMPr, 'OMPr')
        end if 
        call ReadUntilFound(filespec, 'OMPz', reachedeof)
        if (reachedeof) then 
            ! Not found, issue warning and rewind
            print *, 'ReadTraduitUS: could not find field "OMPz", setting to zero...'
            rewind(filespec)
        else 
            ! Found, read in coordinates (assumed two)
            backspace(filespec)
            call ReadArray(filespec, 2, grid%data%OMPz, 'OMPz')
        end if 
    
        ! IMP
        call ReadUntilFound(filespec, 'IMPr', reachedeof)
        if (reachedeof) then 
            ! Not found, issue warning and rewind
            print *, 'ReadTraduitUS: could not find field "IMPr", setting to zero...'
            rewind(filespec)
        else 
            ! Found, read in coordinates (assumed two)
            backspace(filespec)
            call ReadArray(filespec, 2, grid%data%IMPr, 'IMPr')
        end if 
        call ReadUntilFound(filespec, 'IMPz', reachedeof)
        if (reachedeof) then 
            ! Not found, issue warning and rewind
            print *, 'ReadTraduitUS: could not find field "IMPz", setting to zero...'
            rewind(filespec)
        else 
            ! Found, read in coordinates (assumed two)
            backspace(filespec)
            call ReadArray(filespec, 2, grid%data%IMPz, 'IMPz')
        end if 
    
        ! Rewind and reset to expected location
        rewind(filespec)
        call ReadUntilFound(filespec, 'vxFfbz', reachedeof)
        if (reachedeof) then 
            ! Something wrong, exit
            call gdErrorHandler('ReadTraduitUS: could not reset file, check file content')
        else
            backspace(filespec)
        end if
    
        ! Vertex data
        !------------
        ! Here we simply loop ourselves, as ReadArray etc
        ! do not seem to be capable to load in this data (?)
        allocate(vdata(nv, 6), vlist(nv))
        
        ! Skip the header
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
        end if 
    
        ! Read in data
        do i = 1, nv 
            ! Read 
            read(filespec, *) vlist(i), vdata(i,:) ! apparently this works fine...
        end do
    
        ! Add to grid
        grid%vert%x(vlist)      = vdata(:, 1)
        grid%vert%y(vlist)      = vdata(:, 2)
        grid%vert%psi(vlist)    = vdata(:, 3)
        grid%vert%bx(vlist)     = vdata(:, 4)
        grid%vert%by(vlist)     = vdata(:, 5)
        grid%vert%ffbz(vlist)   = vdata(:, 6)
    
        ! Cell data
        !----------
        allocate(cdatai1(nc, 2), cdata(nc, 5), cdatai2(nc, 3), clist(nc))
    
        ! Skip the header
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
        end if 
    
        ! Read in data
        do i = 1, nc 
            ! Read 
            read(filespec, *) clist(i), cdatai1(i,:), cdata(i, :), cdatai2(i, :) ! apparently this works fine...
        end do 
    
        ! Add to grid
        grid%cell%vertP(clist, 1)          = cdatai1(:, 1)
        grid%cell%vertP(clist, 2)          = cdatai1(:, 2)
        grid%cell%reg(clist)               = cdatai2(:, 2)
        grid%cell%ft(clist)                = cdatai2(:, 3)
        grid%cell%cflags(clist)            = cdatai2(:, 1)
        grid%cell%x(clist)                 = cdata(:, 1)
        grid%cell%y(clist)                 = cdata(:, 2)
        grid%cell%psi(clist)               = cdata(:, 3)
        grid%cell%bp(clist)                = cdata(:, 4)
        grid%cell%bt(clist)                = cdata(:, 5)
    
        ! Cell vertices and faces
        call ReadArray (filespec, grid%cell%nvert, grid%cell%vert,  'cvVx')
        call ReadArray (filespec, grid%cell%nface, grid%cell%face,  'cvFc')
    
        ! Faces
        !------
        allocate(fdatai(nf, 5), flist(nf))
    
        ! Skip the header
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
        end if 
    
        ! Read in data
        do i = 1, nf 
            ! Read 
            read(filespec, *) flist(i), fdatai(i,:) ! apparently this works fine...
        end do 
    
        ! Add to grid
        grid%face%vert(flist, 1)            = fdatai(:, 1)
        grid%face%vert(flist, 2)            = fdatai(:, 2)
        grid%face%label(flist)              = fdatai(:, 3)
        grid%face%reg(flist)                = fdatai(:, 4)
        grid%face%aligned                   = fdatai(:, 5) 
        grid%face%TMfacelabel               = grid%face%label
    
        ! Flux tubes
        !-----------
        allocate(ftdatai(grid%data%fluxdata%nFt, 5), &
            ftlist(grid%data%fluxdata%nFt))
    
        ! Skip the header
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
        end if 
    
        ! Read in data
        do i = 1, grid%data%fluxdata%nFt 
            ! Read 
            read(filespec, *) ftlist(i), ftdatai(i,:) ! apparently this works fine...
        end do
    
        ! Compute the number of ftCv and ftFc that are actually there
        nftCv = ftdatai(grid%data%fluxdata%nFt, 2) + &
            ftdatai(grid%data%fluxdata%nFt, 1)-1
        nftFc = ftdatai(grid%data%fluxdata%nFt, 4) + &
            ftdatai(grid%data%fluxdata%nFt, 3)-1
        
        ! Read ftCv, ftFc
        allocate(ftCvdum(nftCv), ftFcdum(nftFc))
        call ReadArray (filespec, nftCv, ftCvdum, 'ftCv')
        call ReadArray (filespec, nftFc, ftFcdum, 'ftFc')
        
    
        ! Add to grid
        grid%data%fluxdata%fluxtubecells(1:nftCv) = ftCvdum
        grid%data%fluxdata%fluxtubefaces(1:nftFc) = ftFcdum
        grid%data%fluxdata%fluxtubecellsP(ftlist, 1) = ftdatai(:, 1)
        grid%data%fluxdata%fluxtubecellsP(ftlist, 2) = ftdatai(:, 2)
        grid%data%fluxdata%fluxtubefacesP(ftlist, 1) = ftdatai(:, 3)
        grid%data%fluxdata%fluxtubefacesP(ftlist, 2) = ftdatai(:, 4)
        grid%data%fluxdata%fluxtuberegID(ftlist) = ftdatai(:, 5)
    
        ! Flux surfaces
        !--------------
        allocate(fsdatai(grid%data%fluxdata%nFs, 2), &
            fsdata(grid%data%fluxdata%nFs, 1), &
            fslist(grid%data%fluxdata%nFs)) 
    
        ! Skip the header
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTraduitUS: reached EOF prematurely')
        end if 
    
        ! Read in data
        do i = 1, grid%data%fluxdata%nFs 
            ! Read 
            read(filespec, *) fslist(i), fsdatai(i,:), fsdata(i, :) ! apparently this works fine...
        end do
    
        ! Compute the number of fsFc
        nfsFc = fsdatai(grid%data%fluxdata%nFs, 2) + &
            fsdatai(grid%data%fluxdata%nFs, 1)-1
    
        ! Read fsFc
        allocate(fsFcdum(nfsFc))
        call ReadArray (filespec, nfsFc, fsFcdum ,  'fsFc')
    
        ! Add to grid
        grid%data%fluxdata%fluxsurfacefaces(1:nfsFc) = fsFcdum
        grid%data%fluxdata%fluxsurfacefacesP(fslist, 1) = fsdatai(:, 1)
        grid%data%fluxdata%fluxsurfacefacesP(fslist, 2) = fsdatai(:, 2)
        grid%data%fluxdata%fluxsurfacepsi(fslist)       = fsdata(:, 1)
        ! fsdata(:, 1) contains psi value of flux surfaces, to be read in in the future?
        
        ! Housekeeping
        !=============
        deallocate(vdata, cdata, fsdata, vlist, flist, clist, ftlist, &
            fslist, cdatai1, cdatai2, fdatai, ftdatai, fsdatai, ftCvdum, &
            ftFcdum, fsFcdum)
        close(filespec)
    
    end subroutine

    subroutine ReadB2fgmtryUS(grid, filepath)

        ! Description
        !============
        ! Read in the necessary grid data from a b2fgmtry_us file for the 
        ! grid deformation module. This routine can serve as wrapper in the 
        ! future to call the dedicated reading and writing routines 
        ! present in b2mod_geo for example. 
        
        ! Notes
        !======
        ! Note 1: typically, there are still spurious x-point vertices and 
        ! guard cells in the b2fgmtry format. These are removed here. To 
        ! this end, we first read in the grid in a classical way. Then, 
        ! we check which vertices and cells have to be removed and 
        ! construct a grid without guard cells etc. It is assumed that all
        ! guard cells are added at the end of the cell list
    
        ! Note 2: contrary to the traduit.out.b2us file, we cannot assume
        ! that the cell vertices are ordened either clockwise or 
        ! counterclockwise. This necessitates an additional sorting loop
        ! to comply with the assumptions made in GOAT. (we could at some
        ! point include this sorting in the grid interconnection computation
        ! routine, but this would lead to unnecessary overhead for traduit
        ! files)
    
        ! Initialize
        !===========
        ! Declare modules
        use mod_polygon
        use mod_readwrite
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Main variables
        type(GridUDT)               :: grid
        character(*), intent(in)    :: filepath
    
        ! Auxiliary
        logical                     :: reachedeof, readTopologicalData
        integer(I8)                 :: idum(0:9), idum2(1), filespec   
        integer(I8)                 :: nc,nf,nv ! total number of cells, faces, vertices
    
        character(120)              :: chardummy   ! dummy array
        character(:), allocatable   :: chardummy2
        integer(I8), allocatable    :: cdummy(:,:), cdummy2(:) ! dummy array
        integer(I8), allocatable    :: fdummy(:,:), fdummy2(:, :) ! dummy array
        integer(I8), allocatable    :: vdummy(:,:) ! dummy array
        integer(I8), allocatable    :: ftdummy(:) ! dummy array
    
        real(R8), allocatable       :: fcQalf(:, :) ! to reconstruct aligned faces
        real(R8), allocatable       :: cdummyr(:,:) ! dummy array
        real(R8), allocatable       :: fdummyr(:,:) ! dummy array
        real(R8), allocatable       :: vdummyr(:,:) ! dummy array
        real(R8), allocatable       :: fsdummyr(:) ! dummy array
        real(R8), allocatable       :: facedummy(:) ! dummy array
        integer(I8), allocatable       :: n2dummy(:) ! dummy array
        integer(I8), allocatable       :: nxdummy(:) ! dummy array        
        real(R8), allocatable       :: vx(:), vy(:), vxbx(:), vxby(:), &
            vffbz(:), vpsi(:)
        integer(I8), allocatable    :: vfieldlineID(:), vface(:), &
            vcell(:), vcell1(:), cellft(:), cellreg(:), &
            ftcells(:), cellv(:), cellf(:)
        logical, allocatable        :: vBV(:)
        integer(I8)                 :: n2,nx  ! legacy structured data
    
        integer(I8)                                 :: ngv, ngc 
        integer(I8), allocatable, dimension(:)      :: vertmap, cellmap, &
            sortindex, tcf, tv, tftc
        integer(I8), allocatable, dimension(:, :)   :: tempvertfaceP, &
            tempvertcellP, tempcellfaceP, tempcellvertP, tcfv
    
        logical, allocatable, dimension(:)      :: isnoghostvert, keepvertface, &
            keepvertcell, isnoguardcell, keepcellface, keepcellvert, &
            ispolygonstart, isbranchingpolygon
    
        ! Loop
        integer(I8)                 :: i, j, k
        
        ! Data
        data filespec /60/
        
        ! Read (inlucding guard cells)
        !=============================
        ! Open the file
        print *, 'reading grid in b2gmtry_us format from file: ' // filepath
        open(unit = filespec, file = filepath)
    
        ! First, read the header with the version - just skip it...
        call ReadSingleLine(filespec, chardummy2, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadB2fgmtryUS: reached EOF prematurely')
        end if 

        ! Check the version to determine what to read in 
        readTopologicalData = .false. 
        if (chardummy2(8:17) >= '03.002.002') then
            ! Topological data should be present
            readTopologicalData = .true.
        end if 
    
        ! Primary array dimensions
        call ReadArray (filespec,7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
        ngc = int(idum(1), I8)
        nc = int(idum(2), I8) ! make sure to cast to correct type
        nf = int(idum(3), I8)
        nv = int(idum(4), I8)
    
        ! Add to grid
        grid%cell%ntot         = nc
        grid%face%ntot         = nf
        grid%vert%ntot          = nv
        grid%data%fluxdata%nFs  = idum(5)
        grid%data%fluxdata%nFt  = idum(6)
    
        ! Secondary array dimensions
        call ReadArray (filespec,5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
    
        ! Add to grid
        grid%cell%nvert = idum(0)
        grid%cell%nface = idum(1)
        grid%vert%ncell  = idum(2)
        grid%vert%nface  = idum(3)
    
        ! Allocate grid
        call AllocateGrid(grid)

        ! Initialize topological data
        if (allocated(grid%data%isprimaryxp)) deallocate(grid%data%isprimaryxp)
        if (allocated(grid%data%xpointID)) deallocate(grid%data%xpointID)
        if (allocated(grid%data%opointID)) deallocate(grid%data%opointID)
        if (allocated(grid%data%spointID)) deallocate(grid%data%spointID)
        if (allocated(grid%data%divFcP)) deallocate(grid%data%divFcP)
        if (allocated(grid%data%divFc)) deallocate(grid%data%divFc)
        if (allocated(grid%data%spointdivID)) deallocate(grid%data%spointdivID)
        if (allocated(grid%data%tpointdivID)) deallocate(grid%data%tpointdivID)

        ! Check
        if (.not. allocated(grid%vert%face)) then 
            allocate(grid%vert%face(grid%vert%nface))
        end if 
        if (.not. allocated(grid%vert%cell)) then 
            allocate(grid%vert%cell(grid%vert%ncell))
        end if 
    
        ! Allocate dumies
        allocate(cdummy(nc,4), cdummy2(nc))
        allocate(fdummy(nf,4))
        allocate(vdummy(nv,4))
        allocate(cdummyr(nc,4))
        allocate(fdummyr(nf,4))
        allocate(vdummyr(nv,4))
        allocate(fsdummyr(grid%data%fluxdata%nFs))
        allocate(facedummy(grid%cell%nface))
        allocate(ftdummy(grid%data%fluxdata%nFt))
        allocate(fcQalf(nf, 2))
    
        ! Read data for structured grid remapping (to be deleted in future)
        call ReadArray (filespec,1,idum2,'isClassicalGrid')
        grid%data%sglegacy%isClassicalGrid = int(idum2(1), I4)
        call ReadArray (filespec,3,idum,'nx,ny,nncut')
        grid%data%sglegacy%nx       = idum(0)
        grid%data%sglegacy%ny       = idum(1)
        grid%data%sglegacy%nncut = idum(2)

        ! Read topological data
        if (readTopologicalData) then 
            ! Read topological mesh flag
            call ReadArray(filespec, 1, idum, 'topoflag')
            grid%data%topoflag = idum(0)
    
            ! Read number of topological points
            call ReadArray(filespec, 6, idum, 'nX,nO,nS,nT,nDiv,nDivFc')
            grid%data%nxp = idum(0)
            grid%data%nop = idum(1)
            grid%data%nsp = idum(2)
            grid%data%ntp = idum(3)
            grid%data%ndiv = idum(4)
            grid%data%ndivFc = idum(5)
            grid%data%nsep = 0
    
        else 
            ! Initialize to zero
            grid%data%topoflag = 0
            grid%data%nxp = 0
            grid%data%nop = 0
            grid%data%nsp = 0
            grid%data%ntp = 0
            grid%data%ndiv = 0
            grid%data%ndivFc = 0
            grid%data%nsep = 0
    
        end if 

        ! Allocate
        allocate(grid%data%xpointID(idum(0)), grid%data%opointID(idum(1)), &
            grid%data%spointID(idum(2)), grid%data%tpointID(idum(3)), &
            grid%data%isprimaryxp(idum(0)), &
            grid%data%divFcP(grid%data%ndiv, 2), grid%data%divFc(grid%data%ndivFc), &
            grid%data%spointdivID(grid%data%nsp), grid%data%tpointdivID(grid%data%ntp), &
            grid%data%spointxpID(grid%data%nsp), grid%data%sepID(grid%data%nsep))
        
        ! Read data that is not used
        call ReadSingleLine(filespec, chardummy2, reachedeof)
        call ReadArray (filespec,1,idum, 'isymm')
    
        ! Add grid mapping data
        allocate(fdummy2(nf, 2))
        call ReadArray (filespec, nc*2,    grid%cell%faceP, 'cvFcP')
        call ReadArray (filespec, grid%cell%nface, grid%cell%face,  'cvFc')
        call ReadArray (filespec, nf*2,    fdummy2,  'fcCv')
        call ReadArray (filespec, nf*2,    grid%face%vert,  'fcVx')
        call ReadArray (filespec, nc*2,    grid%cell%vertP, 'cvVxP')
        call ReadArray (filespec, grid%cell%nvert, grid%cell%vert,  'cvVx')
        call ReadArray (filespec, nv*2,     grid%vert%faceP, 'vxFcP')
        call ReadArray (filespec, grid%vert%nface,  grid%vert%face,  'vxFc')
        call ReadArray (filespec, nv*2,     grid%vert%cellP, 'vxCvP')
        call ReadArray (filespec, grid%vert%ncell,  grid%vert%cell,  'vxCv')
        call ReadArray (filespec, grid%data%fluxdata%nFt*2,   grid%data%fluxdata%fluxtubecellsP, 'ftCvP')
        call ReadArray (filespec, nc,     grid%data%fluxdata%fluxtubecells,  'ftCv')
        call ReadArray (filespec, grid%data%fluxdata%nFt*2,   grid%data%fluxdata%fluxtubefacesP, 'ftFcP')
        call ReadArray (filespec, nf,     grid%data%fluxdata%fluxtubefaces,  'ftFc')
        call ReadArray (filespec, nc,     grid%cell%ft,  'cvFt')   
        call ReadArray (filespec, grid%data%fluxdata%nFs*2,   grid%data%fluxdata%fluxsurfacefacesP, 'fsFcP')
        call ReadArray (filespec, nf,     grid%data%fluxdata%fluxsurfacefaces,  'fsFc')
        call ReadArray (filespec, nf,     grid%face%reg, 'fcReg')
        call ReadArray (filespec, nc,     grid%cell%reg, 'cvReg')
        call ReadArray (filespec, grid%data%fluxdata%nFt,     grid%data%fluxdata%fluxtuberegID, 'ftReg')
        call ReadArray (filespec, grid%cell%nface,  facedummy,'intcellP') ! not used
        call ReadArray (filespec, grid%cell%nface,  facedummy,'intcellR') ! not used
        deallocate(fdummy2)
    
        ! Add additional data from underlying structured grid (to be removed)
        n2 = 0
        nx = 0
        if (grid%data%sglegacy%isClassicalGrid == 1) n2 = (grid%data%sglegacy%nx+2)*(grid%data%sglegacy%ny+2)
        if (grid%data%sglegacy%isClassicalGrid == 1) nx =  grid%data%sglegacy%nx
    
        ! For now, put it into dummies
        allocate(n2dummy(n2))
        allocate(nxdummy(nx))
        if (grid%data%sglegacy%isClassicalGrid == 1) then 
            call ReadArray (filespec, n2, n2dummy,  'imapCv')
            call ReadArray (filespec, n2, n2dummy, 'imapFcx')
            call ReadArray (filespec, n2, n2dummy, 'imapFcy')
            call ReadArray (filespec, n2, n2dummy,  'imapVx')
            call ReadArray (filespec, nx, nxdummy, 'icornVx')
        end if 
    
        ! Read some labels, ignore
        call ReadArray (filespec, nf, grid%face%label,'fcLbl')
        call ReadArray (filespec, nc, cdummy(:, 1),'cvLbl')
        call ReadArray (filespec, grid%data%fluxdata%nFt, ftdummy,'ftLbl')

        ! Read topological data
        if (readTopologicalData) then 
            ! Read X-point data
            call ReadSingleLine(filespec, chardummy2, reachedeof) ! header
            do i = 1, grid%data%nxp  
                ! Read 
                read(filespec, *) grid%data%xpointID(i), grid%data%isprimaryxp(i), &
                    idum(0)
            end do
    
            ! Read O-point data - not yet written in b2fgmtry
            !call ReadSingleLine(filespec, chardummy2, reachedeof) ! header
            !do i = 1, grid%data%nop  
            !    ! Read 
            !    read(filespec, *) grid%data%opointID(i)
            !end do
    
            ! Read strike point data
            call ReadSingleLine(filespec, chardummy2, reachedeof) ! header
            do i = 1, grid%data%nsp  
                ! Read 
                read(filespec, *) grid%data%spointID(i), grid%data%spointxpID(i), &
                    grid%data%spointdivID(i), idum(0)
            end do
    
            ! Read tangency point data
            call ReadSingleLine(filespec, chardummy2, reachedeof) ! header
            do i = 1, grid%data%ntp  
                ! Read 
                read(filespec, *) grid%data%tpointID(i), grid%data%tpointdivID(i), &
                    idum(0)
            end do
    
            ! Read divertor data
            call ReadSingleLine(filespec, chardummy2, reachedeof) ! header
            do i = 1, grid%data%ndiv  
                ! Read 
                read(filespec, *) idum(0), grid%data%divFcP(i, 1), &
                    grid%data%divFcP(i, 2)
            end do
    
            ! Read divertor face data
            call ReadArray(filespec, grid%data%ndivFc, grid%data%divFc, 'divFc')
        end if 
    
        ! Add geometry data
        ! cell data - ignore all
        call ReadArray (filespec, nc*4, cdummyr(:,1:4),   'cvBb')
        call ReadArray (filespec, nc*3, cdummyr(:,1:3),   'cvEb')
        call ReadArray (filespec, nc,   cdummyr(:,1),    'cvX')
        call ReadArray (filespec, nc,   cdummyr(:,1),    'cvY')
        call ReadArray (filespec, nc,   cdummyr(:,1),   'cvSz')
        call ReadArray (filespec, nc,   cdummyr(:,1),   'cvHz')
        call ReadArray (filespec, nc,   cdummyr(:,1),   'cvHx') !WG temp!
        call ReadArray (filespec, nc*2, cdummyr(:,1:2), 'cvQgam')
        call ReadArray (filespec, nc,   cdummyr(:,1),  'cvVol')
    
        ! face quantities - ignore all except fcQalf
        call ReadArray (filespec, nf*4, fdummyr(:,1:4),   'fcBb')
        call ReadArray (filespec, nf,   fdummyr(:,1),    'fcS')
        call ReadArray (filespec, nf*2, fdummyr(:,1:2),   'fcHc')
        call ReadArray (filespec, nf,   fdummyr(:,1),   'fcHt')
        call ReadArray (filespec, nf*2, fdummyr(:,1:2), 'fcQgam')
        call ReadArray (filespec, nf*2, fcQalf, 'fcQalf')
        call ReadArray (filespec, nf*2, fdummyr(:,1:2), 'fcQbet')
        call ReadArray (filespec, nf,   fdummyr(:,1),  'fcPbs')
    
        ! Determine aligned faces as those faces with cos(alpha) = 0
        ! But do it with a very small tolerance to avoid numerical garbage
        where (fcQalf(:, 1) < 1e-10) 
            grid%face%aligned = 1
        elsewhere
            grid%face%aligned = 0
        end where

        ! vertex quantities - only keep coordinates (and ffbz)
        call ReadArray (filespec, nv*4, vdummyr(:,1:4),   'vxBb')
        call ReadArray (filespec, nv,   grid%vert%x,    'vxX')
        call ReadArray (filespec, nv,   grid%vert%y,    'vxY')
        call ReadArray (filespec, nv,   grid%vert%ffbz, 'vxFfbz')
        call ReadArray (filespec, nv,   vdummyr(:,1), 'vxFpsi')
    
        ! flux surface quantities
        call ReadArray (filespec, nc,   cdummyr(:,1), 'cvConn')
        call ReadArray (filespec, grid%data%fluxdata%nFs,   &
            grid%data%fluxdata%fluxsurfacepsi, 'fsPsi')
    
        ! Eliminate 'ghost' vertices
        !===========================
        ! Allocate & initialize
        allocate(isnoghostvert(grid%vert%ntot))
        isnoghostvert(:) = .false.
    
        ! Loop
        do j = 1, 2
            do i = 1, grid%face%ntot
                isnoghostvert(grid%face%vert(i,j)) = .true.
            enddo
        enddo
    
        ! Count number of ghost vertices
        ngv = count(.not. isnoghostvert) ! total number of ghost vertices
    
        ! Construct vertex mapping for vertex arrays
        allocate(vertmap(grid%vert%ntot))
        vertmap = 0
        k = 0
        do i = 1, grid%vert%ntot
            ! Skip if ghostvert
            if (isnoghostvert(i)) then 
    
                ! Update k
                k = k + 1
    
                ! Set mapping
                vertmap(i) = k 
            end if
        end do
    
        ! Rebuild vertex fields
        grid%vert%ntot = grid%vert%ntot - ngv
        vx = grid%vert%x
        vy = grid%vert%y
        vfieldlineID = grid%vert%fieldlineID
        vxbx = grid%vert%bx
        vxby = grid%vert%by
        vBV = grid%vert%BV
        vffbz = grid%vert%ffbz
        vpsi = grid%vert%psi
        deallocate(grid%vert%x, grid%vert%y, grid%vert%fieldlineID, grid%vert%bx, &
        grid%vert%by, grid%vert%BV, grid%vert%ffbz, grid%vert%psi)
        allocate(grid%vert%x(grid%vert%ntot), grid%vert%y(grid%vert%ntot), &
        grid%vert%fieldlineID(grid%vert%ntot), grid%vert%bx(grid%vert%ntot), &
        grid%vert%by(grid%vert%ntot), grid%vert%BV(grid%vert%ntot), &
        grid%vert%ffbz(grid%vert%ntot), grid%vert%psi(grid%vert%ntot))

        grid%vert%x = pack(vx, isnoghostvert)
        grid%vert%y = pack(vy, isnoghostvert)
        grid%vert%fieldlineID = pack(vfieldlineID, isnoghostvert)
        grid%vert%bx = pack(vxbx, isnoghostvert)
        grid%vert%by = pack(vxby, isnoghostvert)
        grid%vert%BV = pack(vBV, isnoghostvert)
        grid%vert%ffbz = pack(vffbz, isnoghostvert)
        grid%vert%psi = pack(vpsi, isnoghostvert)
        if (allocated(grid%vert%neigP)) then 
            deallocate(grid%vert%neigP)
            allocate(grid%vert%neigP(grid%vert%ntot, 2))
        end if 
    
        ! Rebuild face fields
        grid%face%vert(:, 1) = vertmap(grid%face%vert(:, 1))
        grid%face%vert(:, 2) = vertmap(grid%face%vert(:, 2))
    
        ! Rebuild cell fields
        grid%cell%vert = vertmap(grid%cell%vert)
    
        ! Sanity checks
        if (any(grid%cell%vert == 0)) then 
            ! Zero vertices in cells detected, shouldn't happen
            call gdErrorHandler('ReadB2fgmtryUS: cells with zero vertex ID ' // & 
                'detected after removing ghost vertices. Check input grid')
        end if 
        if (any(any(grid%face%vert == 0, 1))) then 
            ! Zero vertices in faces detected, shouldn't happen
            call gdErrorHandler('ReadB2fgmtryUS: faces with zero vertex ID' // & 
                'detected after removing ghost vertices. Check input grid')
        end if 
    
        ! Rebuild vertex interconnection data
        allocate(keepvertface(grid%vert%nface), keepvertcell(grid%vert%ncell), &
            tempvertfaceP(grid%vert%ntot, 2), tempvertcellP(grid%vert%ntot, 2))
        keepvertface = .true. 
        keepvertcell = .true. 
        where (grid%vert%face == 0) keepvertface = .false. ! eliminate zero values 
        where (grid%vert%cell == 0) keepvertcell = .false.
        do i = 1, grid%vert%ntot + ngv
            if (.not. isnoghostvert(i)) then 
                ! Set to false
                keepvertface(grid%vert%faceP(i, 1):grid%vert%faceP(i, 1)+grid%vert%faceP(i, 2)-1) = .false.
                keepvertcell(grid%vert%cellP(i, 1):grid%vert%cellP(i, 1)+grid%vert%cellP(i, 2)-1) = .false.
            end if
        end do 
        tempvertfaceP(:, 1) = pack(grid%vert%faceP(:, 1), isnoghostvert)
        tempvertfaceP(:, 2) = pack(grid%vert%faceP(:, 2), isnoghostvert)
        tempvertcellP(:, 1) = pack(grid%vert%cellP(:, 1), isnoghostvert)
        tempvertcellP(:, 2) = pack(grid%vert%cellP(:, 2), isnoghostvert)
        grid%vert%faceP = tempvertfaceP 
        grid%vert%faceP(1, 1)  = 1
        grid%vert%cellP = tempvertcellP 
        grid%vert%cellP(1, 1) = 1
        do i = 2, grid%vert%ntot
            grid%vert%faceP(i, 1) = grid%vert%faceP(i-1, 1) + grid%vert%faceP(i-1, 2)
            grid%vert%cellP(i, 1) = grid%vert%cellP(i-1, 1) + grid%vert%cellP(i-1, 2)
        end do 
        vface = grid%vert%face
        vcell = grid%vert%cell
        deallocate(grid%vert%face, grid%vert%cell)
        allocate(grid%vert%face(count(keepvertface)), grid%vert%cell(count(keepvertcell)))
        grid%vert%face = pack(vface, keepvertface)
        grid%vert%cell = pack(vcell, keepvertcell)
        grid%vert%nface = size(grid%vert%face)
        grid%vert%ncell = size(grid%vert%cell)
    
        ! Eliminate guard cells
        !======================
        ! Allocate & initialize
        allocate(isnoguardcell(grid%cell%ntot))
        isnoguardcell(:) = .false.
        isnoguardcell(1:grid%cell%ntot-ngc) = .true.
    
        ! Construct cell mapping for cell arrays
        allocate(cellmap(grid%cell%ntot))
        cellmap = 0
        k = 0
        do i = 1, grid%cell%ntot
            ! Skip if ghostvert
            if (isnoguardcell(i)) then 
    
                ! Update k
                k = k + 1
    
                ! Set mapping
                cellmap(i) = k 
            end if
        end do
    
        ! Rebuild vertex fields
        vcell1 = grid%vert%cell
        deallocate(grid%vert%cell)
        allocate(grid%vert%cell(count(isnoguardcell(vcell1))))
        grid%vert%cell = pack(vcell1, isnoguardcell(vcell1))
        grid%vert%cell = cellmap(grid%vert%cell)
    
        ! Rebuild face fields (nothing to be done)
    
        ! Rebuild cell fields
        grid%cell%ntot = grid%cell%ntot - ngc 
        cellft = grid%cell%ft
        cellreg = grid%cell%reg
        deallocate(grid%cell%ft, grid%cell%reg, grid%cell%psi, grid%cell%bp, &
            grid%cell%bt, grid%cell%x, grid%cell%y, grid%cell%cflags)
        allocate(grid%cell%ft(count(isnoguardcell)), grid%cell%reg(count(isnoguardcell)), grid%cell%psi(grid%cell%ntot), &
            grid%cell%bp(grid%cell%ntot), grid%cell%bt(grid%cell%ntot), grid%cell%x(grid%cell%ntot), grid%cell%y(grid%cell%ntot), &
            grid%cell%cflags(grid%cell%ntot))
        grid%cell%ft = pack(cellft, isnoguardcell)
        grid%cell%reg = pack(cellreg, isnoguardcell)
        grid%cell%psi = 0
        grid%cell%bp = 0
        grid%cell%bt = 0
        grid%cell%x = 0
        grid%cell%y = 0 
        grid%cell%cflags = 0
        
        
        ! Rebuild flux tube data
        ftcells = grid%data%fluxdata%fluxtubecells
        deallocate(grid%data%fluxdata%fluxtubecells)
        allocate(grid%data%fluxdata%fluxtubecells(count(ftcells /= 0)))
        grid%data%fluxdata%fluxtubecells = pack(ftcells, ftcells /= 0)        
        grid%data%fluxdata%fluxtubecells = cellmap(grid%data%fluxdata%fluxtubecells)
        do i = 1, grid%data%fluxdata%nft 
            tftc = GetFTCell(grid%data%fluxdata, i)
            grid%data%fluxdata%fluxtubecellsP(i, 2) = grid%data%fluxdata%fluxtubecellsP(i, 2) - count(tftc == 0)
        end do 
        grid%data%fluxdata%fluxtubecellsP(1, 1) = 1
        grid%data%fluxdata%fluxtubefacesP(1, 1) = 1 ! hedge for junk here from input
        do i = 2, grid%data%fluxdata%nft 
            grid%data%fluxdata%fluxtubecellsP(i, 1) = &
                grid%data%fluxdata%fluxtubecellsP(i-1, 1) + grid%data%fluxdata%fluxtubecellsP(i-1, 2)
            grid%data%fluxdata%fluxtubefacesP(i, 1) = &
                grid%data%fluxdata%fluxtubefacesP(i-1, 1) + grid%data%fluxdata%fluxtubefacesP(i-1, 2)
        end do 
    
    
        ! Rebuild cell interconnection data
        allocate(keepcellface(grid%cell%nface), keepcellvert(grid%cell%nvert), &
            tempcellfaceP(grid%cell%ntot, 2), tempcellvertP(grid%cell%ntot, 2))
        keepcellface = .true. 
        keepcellvert = .true.
        where (grid%cell%face == 0) keepcellface = .false.
        where (grid%cell%vert == 0) keepcellvert = .false. 
        do i = 1, grid%cell%ntot + ngc
            if (.not. isnoguardcell(i)) then 
                ! Set to false
                keepcellface(grid%cell%faceP(i, 1):grid%cell%faceP(i, 1)+grid%cell%faceP(i, 2)-1) = .false.
                keepcellvert(grid%cell%vertP(i, 1):grid%cell%vertP(i, 1)+grid%cell%vertP(i, 2)-1) = .false.
            end if
        end do 
        tempcellfaceP(:, 1) = pack(grid%cell%faceP(:, 1), isnoguardcell)
        tempcellfaceP(:, 2) = pack(grid%cell%faceP(:, 2), isnoguardcell)
        tempcellvertP(:, 1) = pack(grid%cell%vertP(:, 1), isnoguardcell)
        tempcellvertP(:, 2) = pack(grid%cell%vertP(:, 2), isnoguardcell)
        grid%cell%faceP = tempcellfaceP 
        grid%cell%faceP(1, 1)  = 1
        grid%cell%vertP = tempcellvertP 
        grid%cell%vertP(1, 1) = 1
        do i = 2, grid%cell%ntot
            grid%cell%faceP(i, 1) = grid%cell%faceP(i-1, 1) + grid%cell%faceP(i-1, 2)
            grid%cell%vertP(i, 1) = grid%cell%vertP(i-1, 1) + grid%cell%vertP(i-1, 2)
        end do 
        cellf = grid%cell%face
        cellv = grid%cell%vert
        deallocate(grid%cell%face, grid%cell%vert)
        allocate(grid%cell%face(count(keepcellface)), grid%cell%vert(count(keepcellvert)))
        grid%cell%face = pack(cellf, keepcellface)
        grid%cell%vert = pack(cellv, keepcellvert)
        grid%cell%nface = size(grid%cell%face)
        grid%cell%nvert = size(grid%cell%vert)
    
        ! Reorden cell vertices
        !======================
        ! Since we know the cell faces from the format, we can do this 
        ! easily by using the SortPolygonEdges and ExtractPolygonVertices 
        ! routines
        do i = 1, grid%cell%ntot 
            ! Get cell faces
            tcf = GetCellFace(grid%cell, i)
    
            ! Get vertices
            tcfv = grid%face%vert(tcf, :)
    
            ! Sort
            allocate(sortindex(size(tcf)), ispolygonstart(size(tcf)), isbranchingpolygon(size(tcf)))
            call SortPolygonEdges(tcfv, size(tcf), sortindex, ispolygonstart, &
                isbranchingpolygon)
            tcfv = tcfv(sortindex, :)
    
            ! Check
            if ( (count(ispolygonstart) /= 1)) then 
                call gdErrorHandler('ReadB2fgmtry: multiple polygons detected ' // &
                    'for single cell, not supported. Check grid input')
            end if 
    
            ! Extract vertices
            allocate(tv(size(tcf)+1))
            call ExtractPolygonVertices(tcfv, size(tcf), tv)
    
            ! Check
            if (tv(1) /= tv(size(tcf)+1)) then 
                call gdErrorHandler('ReadB2fgmtry: cell vertices do not form ' // &
                    'closed polygon, not supported. Check grid input')
            end if 
    
            ! Add to grid vertices
            grid%cell%vert(grid%cell%vertP(i, 1):grid%cell%vertP(i, 1)+grid%cell%vertP(i, 2)-1) = &
                tv(1:size(tcf))
    
            ! Housekeeping
            deallocate(tv, sortindex, ispolygonstart, isbranchingpolygon)
    
        end do
    
        ! Housekeeping
        !=============
        ! Close file
        close(filespec)
    
    end subroutine

    subroutine ReadGoatGGData(grid, filepath)

        ! Description
        !============
        ! This routine reads in additional data provided by the goat
        ! grid generator in goatggdata.dat, if it is available. 

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: grid 
        character(*), intent(in)            :: filepath 

        ! Auxiliary
        integer                             :: fu 
        logical                             :: fileExists, reachedeof 
        character(:), allocatable           :: thisline
        integer(I8)                         :: tmptype, tmpBLind, &
            tfID, nfacelabelmappingGGtoGD

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Check existence in current directory
        print *, 'ReadGoatGGData: looking for ' // filepath // ' in current folder'
        inquire(file=filepath, exist=fileExists)
        if (.not. fileExists) then 
            ! Check existence in output directory
            print *, 'ReadGoatGGData: looking for ' // filepath // ' in output folder'
            inquire(file=plotdir // filesepchar // filepath, exist=fileExists)
            if (fileExists) then 
                ! Open file
                open (action='read', file=plotdir // filesepchar // filepath, newunit=fu, &
                    status='unknown')
            end if 
        else
            ! Open file
            open (action='read', file=filepath, newunit=fu, &
                status='unknown')
        end if 

        ! Check if we can continue
        if (.not. fileExists) then 
            ! Set to false and exit
            print *, 'ReadGoatGGData: could not find file, not reading data...'
            grid%data%hasGoatGGData = .false. 
            return 
        else 
            ! Set to true and continue
            grid%data%hasGoatGGData = .true. 
        end if 

        ! Initialize arrays
        if (allocated(grid%data%goatggdata%TMfacetype)) deallocate(grid%data%goatggdata%TMfacetype)
        if (allocated(grid%data%goatggdata%TMverttype)) deallocate(grid%data%goatggdata%TMverttype)
        if (allocated(grid%data%goatggdata%BLind)) deallocate(grid%data%goatggdata%BLind)
        allocate(grid%data%goatggdata%TMfacetype(grid%face%ntot), &
            grid%data%goatggdata%BLind(grid%face%ntot), grid%data%goatggdata%TMverttype(grid%vert%ntot))

        ! Read
        !=====
        ! Skip header
        call ReadSingleLine(fu, thisline, reachedeof)

        ! Read face data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, grid%face%ntot
            read(fu, *) tfID, tmptype, tmpBLind
            grid%data%goatggdata%TMfacetype(tfID)    = tmptype
            grid%data%goatggdata%BLind(tfID)         = tmpBLind
        end do 

        ! Read vertex data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, grid%vert%ntot
            read(fu, *) tfID, tmptype
            grid%data%goatggdata%TMverttype(tfID)    = tmptype
        end do 

        ! Read face label mappings
        call ReadSingleLine(fu, thisline, reachedeof) ! header
        read (fu, *) nfacelabelmappingGGtoGD 
        allocate(grid%data%goatggdata%facelabelsGG(nfacelabelmappingGGtoGD))
        grid%data%goatggdata%facelabelsGD = grid%data%goatggdata%facelabelsGG
        call ReadSingleLine(fu, thisline, reachedeof) ! header
        do i = 1, nfacelabelmappingGGtoGD
            read(fu, *) grid%data%goatggdata%facelabelsGG(i)
        end do 
        call ReadSingleLine(fu, thisline, reachedeof) ! header
        do i = 1, nfacelabelmappingGGtoGD
            read(fu, *) grid%data%goatggdata%facelabelsGD(i)
        end do 

        ! Close file
        close(fu)


    end subroutine

    ! Writers
    subroutine WriteGOAT(goatoptions, grid, magneticField, environment)

        ! Description
        !============
        ! This routine writes out the grid data from GOAT to the output file
        ! specified in the goatoptions structure. The format is an 
        ! unstructured traduit.out.b2us file 
    
        ! Note: the format madness that you will find in this file is due to
        ! the overly pedantic read/write routines of SOLPS, which fail if 
        ! you dare to even add/remove a single whitespace. If, hopefully
        ! somewhere in the near future, this would get revision, one can
        ! get rid of all the 'fmt' business and just put 'write(fu, *)' 
        ! probably everywhere. 
    
        ! Modules
        !========
        use goatmod_userinput
        use mod_std_formatspecs
        use mod_constants
        use mod_definitions
    
        ! The usual
        !==========
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        type(GoatoptionsUDT), intent(in)        :: goatoptions 
        type(GridUDT), intent(inout)            :: grid 
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(EnvironmentUDT), intent(in)        :: environment
    
        ! Auxiliary
        character(:), allocatable               :: version, tempstring, fmt, &
            gridversion
        integer                                 :: fu 
    
        ! Loop
        integer(I8)                             :: i 
    
        ! Initialize
        !===========
        ! Check if we need to write
        if (.not. goatoptions%write_final) then 
            ! Return
            return 
        end if 
    
        ! Open file, overwrite if existing
        open (action='write', file=trim(goatoptions%writefilepath), newunit=fu, &
            status='unknown')
    
        ! Recompute data
        call UpdateGridData(grid, magneticField, environment)
    
        ! Write grid (unstructured traduit format)
        !=========================================
        ! Associate
        associate(&
            mf              => magneticField%interp,    &
            xv              => grid%vert%x,         &
            yv              => grid%vert%y,         &
            nv              => grid%vert%ntot,      &
            psiv            => grid%vert%psi,       &
            Bxv             => grid%vert%bx,        &
            Byv             => grid%vert%by,        &
            ffbzv           => grid%vert%ffbz,      &
            nvertcell       => grid%vert%ncell,     &
            nvertface       => grid%vert%nface,     &
            fieldlineID     => grid%vert%fieldlineID,   &
            isprimaryxp     => grid%data%isprimaryxp,   &
            SpointxpID      => grid%data%spointxpID,    &
            SpointdivID     => grid%data%spointdivID,   &
            Xpoint          => grid%data%xpointID,  &
            Opoint          => grid%data%opointID,  &
            Spoint          => grid%data%spointID,  &
            Tpoint          => grid%data%tpointID,  &
            TpointdivID     => grid%data%tpointdivID,   &
            topoflag        => grid%data%topoflag,  &
            nX              => grid%data%nxp,       &
            nO              => grid%data%nop,       &
            nS              => grid%data%nsp,       &
            nT              => grid%data%ntp,       &
            nDiv            => grid%data%ndiv,      &
            ndivFc          => grid%data%ndivFc,    &
            divFc           => grid%data%divFc,     &
            divFcP          => grid%data%divFcP,    &
            nf              => grid%face%ntot,      &
            facevert        => grid%face%vert,      &
            labelf          => grid%face%label,     &
            regf            => grid%face%reg,       &
            alignedf        => grid%face%aligned,   &
            fcOMP           => grid%data%OMPface,   &
            nfcOMP          => grid%data%nOMPface,  &    
            fcIMP           => grid%data%IMPface,   &
            nfcIMP          => grid%data%nIMPface,  &
            nFmxCv          => grid%face%ncell,     &
            nc              => grid%cell%ntot,      &
            ngc             => grid%cell%ngc,       & 
            cellvert        => grid%cell%vert,      &
            cellvertP       => grid%cell%vertP,     &
            xc              => grid%cell%x,         &
            yc              => grid%cell%y,         &
            psic            => grid%cell%psi,       &
            bpc             => grid%cell%bp,        &
            btc             => grid%cell%bt,        &
            flagc           => grid%cell%cflags,    &
            regc            => grid%cell%reg,       &
            cellface        => grid%cell%face,      &
            ncellvert       => grid%cell%nvert,     &
            ncellface       => grid%cell%nface,     &
            cellft          => grid%cell%ft,        &    
            cvOMP           => grid%data%OMPcell,       &
            ncvOMP          => grid%data%nOMPcell,      &
            cvIMP           => grid%data%IMPcell,       &
            ncvIMP          => grid%data%nIMPcell,      &
            nft             => grid%data%fluxdata%nft,                  &
            ftcell          => grid%data%fluxdata%fluxtubecells,        &
            ftcellP         => grid%data%fluxdata%fluxtubecellsP,       &
            ftface          => grid%data%fluxdata%fluxtubefaces,        &
            ftfaceP         => grid%data%fluxdata%fluxtubefacesP,       &
            ftreg           => grid%data%fluxdata%fluxtuberegID,        &
            nfs             => grid%data%fluxdata%nfs,                  &
            fsface          => grid%data%fluxdata%fluxsurfacefaces,     &
            fsfaceP         => grid%data%fluxdata%fluxsurfacefacesP,    &
            fspsi           => grid%data%fluxdata%fluxsurfacepsi,       &
            sgnx            => grid%data%sglegacy%nx,       &
            sgny            => grid%data%sglegacy%ny,       &
            sgncut          => grid%data%sglegacy%nncut,    &
            isClassicalGrid => grid%data%sglegacy%isclassicalgrid       &
            )
    
        ! Write
        !======
        ! Version
        if (goatoptions%write_topologicaldata) then 
            gridversion = '03.002.002'
        else
            gridversion = '03.002.000'
        end if 
        version = 'VERSION' // gridversion // ' traduit.out.b2us'
        write(fu, '(a)') version 
    
        ! General grid information
        tempstring = '*cf:    int                6    nCi,nFc,nVx,nCg,nFs,nFt'
        write(fu, '(a)' ) tempstring
        fmt = '(6'//Ifm//')'
        write(fu, fmt) nc, nf, nv, ngc, nfs, nft
    
        tempstring = '*cf:    int                5    nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc'
        write(fu, '(a)' ) tempstring
    
        fmt = '(5'//Ifm//')'
        write(fu, fmt) ncellvert, ncellface, nFmxCv, nvertcell, nvertface 
        tempstring = '*cf:    int                1    isClassicalGrid'
        write(fu, '(a)' ) tempstring
    
        fmt = '('//Ifm//')'
        write(fu, fmt) isClassicalGrid
    
        if (isClassicalGrid == 1) then 
            tempstring = '*cf:    int                3    nx,ny,nncut'
            write(fu, '(a)' ) tempstring
            fmt = '(3'//Ifm//')'
            write(fu, fmt) sgnx, sgny, sgncut
        end if 
        
        if (goatoptions%write_topologicaldata) then 
            ! Topological flag
            tempstring = '*cf:    int                1    topoflag'
            write(fu, '(a)' ) tempstring
            fmt = '(1'//Ifm//')'
            write(fu, fmt) topoflag
    
            ! Topological points information
            tempstring = '*cf:    int                6    nX,nO,nS,nT,nDiv,nDivFc'
            write(fu, '(a)' ) tempstring
            fmt = '(6'//Ifm//')'
            write(fu, fmt) nX, nO, nS, nT, nDiv, nDivFc
    
            ! X-point data
            tempstring = '*cf: Vx isprimary fsID'
            write(fu, '(a)' ) tempstring 
            fmt = '(3'//Ifm//')'
            do i = 1, nX 
                write(fu, fmt) Xpoint(i), isprimaryxp(i), fieldlineID(Xpoint(i))
            end do
    
            ! O-point data
            tempstring = '*cf: Vx'
            write(fu, '(a)' ) tempstring 
            do i = 1, nO 
                fmt = '('//Ifm//')' 
                write (fu, fmt) Opoint(i)
            end do 
    
            ! S-point data
            tempstring = '*cf: Vx xpID divID fsID'
            write(fu, '(a)' ) tempstring 
            do i = 1, nS 
                fmt = '(4'//Ifm//')' 
                write (fu, fmt) Spoint(i), SpointxpID(i), SpointdivID(i), fieldlineID(Spoint(i))
            end do
    
            ! T-point data
            tempstring = '*cf: Vx divID fsID'
            write(fu, '(a)' ) tempstring 
            do i = 1, nT 
                fmt = '(3'//Ifm//')' 
                write (fu, fmt) Tpoint(i), TpointdivID(i), fieldlineID(Tpoint(i))
            end do
    
            ! Divertor target data
            tempstring  = '*cf: div divFcP(:,1) divFcP(:,2)'
            write(fu, '(a)' ) tempstring 
            do i = 1, nDiv
                fmt = '(3'//Ifm//')'
                write (fu, fmt) i, divFcP(i, 1), divFcP(i, 2)
            end do
    
            ! Divertor face list
            tempstring = repeat(' ', 1000)
            write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', ndivFc, 'divFc'
            fmt = repeat(' ', 1000)
            write(fmt, *) '(', 12, '(', Ifm, '))'
            write(fu, fmt) divFc
        end if 
    
        ! Vertex information
        tempstring = '*cf: Vx vxX vxY vxPsi vxBx vxBy vxFfbz'
        write(fu, '(a)' ) tempstring 
        fmt = '(' // Ifm // ',' // &
     &          repeat(spacefm // ',' // Rfm // ',', 5) // &
     &          spacefm // ',' // Rfm // ')'
        do i = 1, nv 
            write(fu, fmt) i, xv(i), yv(i), psiv(i), Bxv(i), Byv(i), ffbzv(i)
        end do
    
        ! Cell information
        tempstring = '*cf: cv cvVxP(:,1) cvVxP(:,2) cvX cvY psi bp bt cflags(:) cvReg cvFt'
        write(fu, '(a)' ) tempstring 
        fmt = '(' // &
     &          repeat(Ifm // ',' // spacefm // ',', 3) // &
     &          repeat(Rfm // ',' // spacefm // ',', 5) // &
     &          repeat(Ifm // ',' // spacefm // ',', 2) // &
     &          Ifm // ')'
        do i = 1, nc 
            write(fu, fmt) i, cellvertP(i, 1), cellvertP(i, 2), xc(i), yc(i), &
                psic(i), bpc(i), btc(i), flagc(i), regc(i), cellft(i)
        end do
        fmt = repeat(' ', 1000)
        write(fmt, *) '(', 12, '(', Ifm, '))'
        write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ncellvert, 'cvVx'
        write(fu, fmt) cellvert
    
        fmt = repeat(' ', 1000)
        write(fmt, *) '(', 12, '(', Ifm, '))'
        write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ncellface, 'cvFc'
        write(fu, fmt) cellface
    
        if (goatoptions%write_OMPdata) then 
            fmt = repeat(' ', 1000)
            write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', ncvOMP, 'cvOMP'
            write(fmt, *) '(', 12, '(', Ifm, '))'
            write(fu, fmt) cvOMP
            write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', ncvIMP, 'cvIMP'
            write(fmt, *) '(', 12, '(', Ifm, '))'
            write(fu, fmt) cvIMP
        end if 
    
        ! Face information
        tempstring = '*cf: fc fcVx(:,1) fcVx(:,2) fcLbl fcReg fcAligned'
        write(fu, '(a)' ) tempstring
        fmt = '(' // &
     &          repeat(Ifm // ',' // spacefm // ',', 5) // &
     &          Ifm // ')'
        do i = 1, nf 
              write(fu, fmt) i, facevert(i, 1), facevert(i, 2), labelf(i), &
                 regf(i), alignedf(i)
        end do
    
        if (goatoptions%write_OMPdata) then 
            tempstring = repeat(' ', 1000)
            write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', nfcOMP, 'fcOMP'
            fmt = repeat(' ', 1000)
            write(fmt, *) '(', 12, '(', Ifm, '))'
            write(fu, fmt) fcOMP
            write(fu, '(2a8,i12,4x,a5)') '*cf:    ','int     ', nfcIMP, 'fcIMP'
            fmt = repeat(' ', 1000)
            write(fmt, *) '(', 12, '(', Ifm, '))'
            write(fu, fmt) fcIMP
        end if
    
        ! Flux tube information
        tempstring = '*cf: ft ftCvP(:,1) ftCvP(:,2) ftFcP(:,1) ftFcP(:,2) ftReg'
        write(fu, '(a)' ) tempstring
        fmt = '(' // &
     &      repeat(Ifm // ',' // spacefm // ',', 5) // &
     &      Ifm // ')'
        do i = 1, nft
              write(fu, fmt) i, ftcellP(i, 1), ftcellP(i, 2), ftfaceP(i, 1), &
                 ftfaceP(i, 2), ftreg(i)
        end do
    
        tempstring = repeat(' ', 1000)
        write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ftcellP(nft, 1)+ftcellP(nft, 2)-1, 'ftCv'
        fmt = repeat(' ', 1000)
        write(fmt, *) '(', 12, '(', Ifm, '))'
        write(fu, fmt) ftcell(1:ftcellP(nft, 1)+ftcellP(nft, 2)-1) 
        write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', ftfaceP(nft, 1)+ftfaceP(nft, 2)-1, 'ftFc'
        fmt = repeat(' ', 1000)
        write(fmt, *) '(', 12, '(', Ifm, '))'  
        write(fu, fmt) ftface(1:ftfaceP(nft, 1)+ftfaceP(nft, 2)-1)
    
        ! Flux surface information
        tempstring = '*cf: fs fsFcP(:,1) fsFcP(:,2) fsPsi'
        write(fu, '(a)' ) tempstring
        fmt = '('//repeat(Ifm, 4) //')'
        do i = 1, nfs
              write(fu, *) i, fsfaceP(i, 1), fsfaceP(i, 2), fspsi(i)
        end do
    
        tempstring = repeat(' ', 1000)
        write(fu, '(2a8,i12,4x,a4)') '*cf:    ','int     ', fsfaceP(nfs, 1)+fsfaceP(nfs, 2)-1, 'fsFc'
        fmt = repeat(' ', 1000)
        write(fmt, *) '(', 12, '(', Ifm, '))'  
        write(fu, fmt) fsface(1:fsfaceP(nfs, 1)+fsfaceP(nfs, 2)-1)
    
        ! Housekeeping
        !=============
        ! Close file
        close (fu)
    
        ! Additional data output
        !=======================
        ! Write vessel
        call environment%vessel%polygonset%WriteData(goatoptions%writefilepath // '_vesselpolygonset')

        ! Write additional goat data
        if (grid%data%hasGoatGGData) then 
            call WriteGoatGGData(grid, 'goatggdata.dat')
        end if 
    
        ! Write grid in .ogr format
        open (action='write', file=goatoptions%writefilepath // '.ogr', newunit=fu, &
            status='unknown')
    
        do i = 1, nf
            ! Coordinates in mm!
            write (fu, *) xv(facevert(i, 1))*1000, yv(facevert(i, 1))*1000
            write (fu, *) xv(facevert(i, 2))*1000, yv(facevert(i, 2))*1000
            write (fu, *) ' '
        end do 
    
        close (fu)
        ! End associate
        end associate
    
    
    end subroutine

    subroutine WriteGoatGGData(grid, filepath)

        ! Description
        !============
        ! This routine writes additional data provided by the goat
        ! grid generator in goatggdata.dat, if it is available. 

        ! Modules
        !========
        use mod_definitions, only: goatversion

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in)           :: grid 
        character(*), intent(in)            :: filepath 

        ! Auxiliary
        integer                             :: fu 

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Check if it is available
        if (.not. grid%data%hasGoatGGData) then 
            print *, 'WriteGoatGGData: data not available, not writing file'
            return 
        end if 

        ! Open file
        open (action='write', file=plotdir // filesepchar // filepath, newunit=fu, &
            status='unknown')

        ! Unpack
        associate(&
            face            => grid%face,                   &
            vert            => grid%vert,                   &
            facelabelsGG    => grid%data%goatggdata%facelabelsGG,   &
            facelabelsGD    => grid%data%goatggdata%facelabelsGD,   &
            TMfacetype      => grid%data%goatggdata%TMfacetype,  &
            TMverttype      => grid%data%goatggdata%TMverttype,  &
            BLind           => grid%data%goatggdata%BLind        &
            )

        ! Write
        !======
        ! Write version
        write (fu, *) goatversion 
        
        ! Write face data 
        write (fu, *) 'faces: ID, TMfacetype, BLind'
        do i = 1, face%ntot
            write (fu, *) i, TMfacetype(i), BLind(i)
        end do 

        ! Write vertex data
        write (fu, *) 'vertices: ID, TMverttype'
        do i = 1, vert%ntot
            write (fu, *) i, TMverttype(i)
        end do 

        ! Write face label mappings
        write (fu, *) 'nfacelabelmappingGGtoGD'
        write (fu, *) size(facelabelsGG)
        write (fu, *) 'facelabelmappingGG'
        do i = 1, size(facelabelsGG)
            write(fu, *) facelabelsGG(i)
        end do 
        write (fu, *) 'facelabelmappingGD'
        do i = 1, size(facelabelsGD)
            write(fu, *) facelabelsGD(i)
        end do 

        ! Housekeeping
        !=============
        end associate

        ! Close file
        close(fu)

    end subroutine

    subroutine WriteGridData(grid, filename)

        ! Description
        !============
        ! Writing out the GridUDT to plot it.

        ! 'vertices'
        ! <vert%ntot>
        ! 'ID, x, y'
        ! <ID, x, y>
        ! 'faces'
        ! <face%ntot> 
        ! 'ID, v1, v2, label'
        ! <ID, v1, v2, label>
        ! 'cells'
        ! <cell%ntot, cell%nvert> 
        ! 'ID, vp1, vp2, region>'
        ! <ID, vp1, vp2, region>
        ! 'cell vertices'
        ! <cell%vert> 

        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only : filesepchar


        ! Arguments
        class(GridUDT)                        :: grid
        character(*), intent(in)                :: filename 

        ! Auxiliary
        integer                                 :: fu
        real(R8), allocatable, dimension(:)     :: x, y, cx, cy
        integer(I8), allocatable, dimension(:)  :: v1, v2, region, &   
            label, vc, aligned
        character(:), allocatable               :: dir

        ! Loop
        integer(I8)                             :: i


        ! Initialize
        !===========
        ! Unpack
        associate(&
            f       => grid%face,   &
            c       => grid%cell,   &
            v       => grid%vert    &
        )

        ! Construct writing directory
        dir = plotdir // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Write header
        write(fu, *) 'VERSION3.00.00'

        ! Write vertex data
        !==================
        ! Unpack
        x = v%x
        y = v%y

        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) v%ntot 

        ! Vertex data
        write (fu, *) 'ID, x, y'
        do i = 1, v%ntot 
            write (fu, *) i, x(i), y(i)
        end do 

        ! Write face data
        !================
        ! Unpack
        v1 = f%vert(:,1)
        v2 = f%vert(:,2)
        label = f%label
        aligned = f%aligned

        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) f%ntot

        ! Face data
        write (fu, *) 'ID, v1, v2, label, aligned'
        do i = 1, f%ntot
            write (fu, *) i, v1(i), v2(i), label(i), aligned(i)
        end do 

        ! Write cell data
        !================
        ! Unpack
        vc = c%vert
        v1 = c%vertP(:,1)
        v2 = c%vertP(:,2)
        region = c%reg
        cx = c%x
        cy = c%y

        ! Number of cells
        write (fu, *) 'cells'
        write (fu, *) c%ntot, size(vc)

        ! Cell data
        write (fu, *) 'ID, vp1, vp2, region, x, y'
        do i = 1, c%ntot
            write (fu, *) i, v1(i), v2(i), region(i), cx(i), cy(i)
        end do 

        ! Cell vertices
        write (fu, *) 'cell vertices'
        do i = 1, size(vc)
            write (fu, *) vc(i)
        end do 

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)

    end subroutine

    ! Data extraction
    subroutine ExtractGridData(grid, meth, gridoptions)

        ! Description
        !============
        ! Process the grid data given in 'grid' by adding the necessary 
        ! additional datastructures, checking for ghost vertices, ... 
    
        ! Notes
        !======
        ! Note 1: it is assumed that the grid that is read does not contain
        ! any guard cells yet. Otherwise this routine has to be extended. 
        ! Additionally, it is assumed that no 'ghost' vertices (i.e. vertices 
        ! without any other connection to the grid) are present. 
    
        ! Initialize
        !===========
        ! Declare modules
        use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: grid
        type(GridOptionsUDT), intent(inout) :: gridoptions
        character(len=*), intent(in)        :: meth
    
        ! Loop variables
        integer(I8)                 :: i, j, k, iFT, ib, il 
    
        ! Auxiliary variables (
        type(FluxDataUDT)           :: fluxdata
    
        integer(I8)                 :: itf, ntf, nbnd, nfpb, nseg, &
            nlabels 
        integer(I8), allocatable    :: tf(:), tfv(:,:)
        integer(I8), allocatable    :: &
            sortindex(:), temparray(:,:), tempfaces(:), segstart(:), &
            gglabels(:), gdlabels(:)
    
        logical, allocatable        :: mask(:), ispolygonstart(:), &
            isbranchingpolygon(:)
    
        integer(I8), allocatable    :: facevec(:) ! simply 1:grid%face%ntot
        integer(I8)                 :: segend
    
        ! Plotting
        real(R8)                    :: NaN
        integer(I8)                 :: ntemp
        real(R8), allocatable       :: tempx(:), tempy(:) ! coordinates
        logical                     :: makedebugplots = .false. 
    
        ! Initialize
        !===========
        ! Initialize
        grid%vert%fieldlineID = 0
        grid%data%fluxdata%fluxsurfaceID = 0
    
        ! Check method
        !=============
        select case(meth)
    
        case ('traduitb2us', 'b2fgmtry_us')
    
            ! Process the grid
            !=================
            ! Extract the fieldline ID's of the vertices based on the flux 
            ! surface data
    
            ! Allocate
            fluxdata = grid%data%fluxdata
            allocate(tf(maxval(fluxdata%fluxsurfacefacesP(:,2),1)))
            allocate(tfv(maxval(fluxdata%fluxsurfacefacesP(:,2),1),2))
    
            ! Loop
            do iFT = 1, fluxdata%nFs
               ! Unpack
                itf = fluxdata%fluxsurfacefacesP(iFT,1); ! start index
                ntf = fluxdata%fluxsurfacefacesP(iFT,2); ! number of faces 
            
                ! Extract faces
                tf(1:ntf) = fluxdata%fluxsurfacefaces(itf:itf+ntf-1)
    
                ! Extract vertices of these faces
                tfv(1:ntf,:) = grid%face%vert(tf(1:ntf),:)
    
                ! Set the flux tube index
                do j = 1, 2
                    do i = 1, ntf 
                        grid%data%fluxdata%fluxsurfaceID(tfv(i,j)) = iFT
                        grid%vert%fieldlineID(tfv(i,j)) = iFT
                    enddo
                enddo
            enddo
    
            ! Deallocate
            deallocate(tf)
            deallocate(tfv)
    
            ! Extract boundaries
            !===================
            ! Check if any goat grid generator data is available. If so, 
            ! determine mapping automatically 
            if (grid%data%hasGoatGGData) then 
                ! Print a message that we're overriding the user defined mapping
                print *, 'ExtractGridData: goat grid generator data found, '  // & 
                    'overwriting face label mapping'
                print *, 'facelabelmappingGG: ', grid%data%goatggdata%facelabelsGG
                print *, 'facelabelmappingGD: ', grid%data%goatggdata%facelabelsGD

                ! Overwrite
                gridoptions%facelabelmappingGG = grid%data%goatggdata%facelabelsGG
                gridoptions%facelabelmappingGD = grid%data%goatggdata%facelabelsGD
            end if 

            ! Get the mapping between boundary labels 
            gglabels = gridoptions%facelabelmappingGG
            gdlabels = gridoptions%facelabelmappingGD 
    
            ! Substitute labels
            do i = 1, size(gridoptions%facelabelsubfrom)
                where (grid%face%label == gridoptions%facelabelsubfrom(i)) &
                    grid%face%label = gridoptions%facelabelsubto(i)
            end do
    
            ! Loop over all face labels (not regions here!) to precompute
            ! number of grid boundaries (can be more/less)
            allocate(mask(grid%face%ntot))
            allocate(facevec(grid%face%ntot))
            facevec(:) = (/(i, i=1,grid%face%ntot,1)/)
            nlabels = size(gglabels) 
            nbnd = 0 ! number of boundaries
            do il = 1, nlabels 
                ! Get the faces of this boundary
                mask(:) = grid%face%label == gglabels(il);
                nfpb = count(mask)
    
                ! Check
                if (nfpb == 0) then 
                    ! this will not become a boundary, skip rest of loop
                    cycle
                end if
    
                ! Extract faces
                allocate(tempfaces(nfpb))
                tempfaces = pack(facevec, mask) 
    
                ! Determine number of boundaries by sorting
                allocate(sortindex(nfpb))
                allocate(ispolygonstart(nfpb), isbranchingpolygon(nfpb))
                allocate(temparray(nfpb, 2))
                temparray(:, :) = grid%face%vert(tempfaces, :)
                call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart, &
                    isbranchingpolygon)
                nbnd = nbnd + count(ispolygonstart)
    
                ! Housekeeping
                deallocate(tempfaces, sortindex, ispolygonstart, temparray, isbranchingpolygon)
    
            end do 
    
            ! Extract boundaries
            allocate(grid%bnd(nbnd))
            ib = 0 ! boundary counter
            do il = 1, nlabels
                
    
                ! Get the faces of this boundary
                mask(:) = grid%face%label == gglabels(il);
                nfpb = count(mask)
    
                ! Check
                if (nfpb == 0) then 
                    ! Don't add as a boundary, skip rest of loop
                    cycle 
                end if
    
                ! Sort the boundary faces
                allocate(sortindex(nfpb), ispolygonstart(nfpb), &
                    tempfaces(nfpb), temparray(nfpb,2), isbranchingpolygon(nfpb))
    
                tempfaces = pack(facevec, mask)
                temparray(:,:) = grid%face%vert(tempfaces, :)
                call SortPolygonEdges(temparray, nfpb, sortindex, ispolygonstart, &
                    isbranchingpolygon)
                tempfaces = tempfaces(sortindex)
                
                ! Loop over all found boundary segments
                nseg = count(ispolygonstart)
                allocate(segstart(nseg))
                segstart = pack([(k, k = 1, nfpb)], ispolygonstart)
                do j = 1, nseg
                    ! Update the boundary counter
                    ib = ib + 1
                    
                    ! Compute end segment index
                    if (j < nseg) then 
                        segend = segstart(j+1)-1
                    else
                        segend = nfpb
                    end if 
    
                    ! Add the boundary ID
                    grid%bnd(ib)%ID = gdlabels(il)
                    grid%bnd(ib)%nface = segend-segstart(j)+1
    
                    ! Allocate this boundary
                    call AllocateBnd(grid%bnd(ib))
    
                    ! Add faces (already sorted before)
                    grid%bnd(ib)%face(:) = tempfaces(segstart(j):segend)
    
                    ! Extract vertices
                    call ExtractPolygonVertices( & 
                        grid%face%vert(grid%bnd(ib)%face,:), &
                        grid%bnd(ib)%nface, grid%bnd(ib)%vert)
                end do 
    
                ! Deallocate
                deallocate(sortindex, ispolygonstart, temparray, tempfaces, &
                    segstart, isbranchingpolygon)
            end do
    
            ! Make a plot to check
            if (makedebugplots) then
                ntemp = 0
                do ib = 1, nbnd
                    ntemp = ntemp + grid%bnd(ib)%nvert + 1
                end do 
                allocate(tempx(ntemp))
                allocate(tempy(ntemp))
                NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)
                ntemp = 1
                do ib = 1, nbnd
                    !print *, grid%bnd(ib)%vert
                    tempx(ntemp:ntemp+grid%bnd(ib)%nvert-1) = &
                       grid%vert%x(grid%bnd(ib)%vert)
                    tempy(ntemp:ntemp+grid%bnd(ib)%nvert-1) = &
                        grid%vert%y(grid%bnd(ib)%vert)
                    tempx(ntemp+grid%bnd(ib)%nvert) = NaN
                    tempy(ntemp+grid%bnd(ib)%nvert) = NaN
                    ntemp = ntemp + grid%bnd(ib)%nvert + 1
                end do
    
                ntemp = size(tempx)
                call Plot2DPolygon(tempx, tempy, ntemp, '-p')
                deallocate(tempx)
                deallocate(tempy)
            end if
    
    
        case default
    
            ! Unknown extraction case, throw error
            call gdErrorHandler('ExtractGridData: unknown method')
    
        end select 
        
    end subroutine

    ! Vertex substructure
    subroutine AllocateVertices(vert)

        ! Description
        !============
        ! Allocate the fields in the vertex structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of vertices

        ! Note: derived quantities, such as the vertex neighbours, 
        ! cannot be allocated beforehand. These should be done 
        ! separately, e.g. in the ComputeGridInterconnections routine


        ! The usual
        implicit none

        ! Declare variables
        type(VertexUDT)       :: vert

        ! Allocate
        !=========
        ! Vertex data
        allocate(vert%x(vert%ntot))
        allocate(vert%y(vert%ntot))
        allocate(vert%BV(vert%ntot))
        allocate(vert%fieldlineID(vert%ntot))
        allocate(vert%bx(vert%ntot), vert%by(vert%ntot), &
            vert%ffbz(vert%ntot), vert%psi(vert%ntot))

        ! Face data
        allocate(vert%faceP(vert%ntot,2))
        !allocate(vert%face(vert%nface)) 

        ! Cell data
        allocate(vert%cellP(vert%ntot,2))
        !allocate(vert%cell(vert%ncell))

        ! Neighbour data
        allocate(vert%neigP(vert%ntot, 2))
        !allocate(vert%neig(vert%nneig))

    end subroutine

    ! Face substructure
    subroutine AllocateFaces(face)

        ! Description
        !============
        ! Allocate the fields in the faces structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of faces

        ! Note: the cell structure cannot be allocated a priori - see 
        ! ComputeGridInterconnections


        ! The usual
        implicit none

        ! Declare variables
        type(FaceUDT)       :: face

        ! Allocate
        !=========
        ! Vertices
        allocate(face%vert(face%ntot,2))

        ! Cells
        allocate(face%cellP(face%ntot,2)) 
        ! allocate(face%cell(face%ncell))
        allocate(face%BF(face%ntot))

        ! Other
        allocate(face%aligned(face%ntot), face%label(face%ntot), &
            face%reg(face%ntot), face%TMfacelabel(face%ntot))

    end subroutine

    ! Cell substruture
    subroutine AllocateCells(cell)

        ! Description
        !============
        ! Allocate the fields in the cells structure. At least the 
        ! following scalar fields have to be present:
        !
        ! - ntot        : total number of cells
        ! - nvert       : length of cells%vert
        ! - nfacelist   : length of cells%face
        ! 
        ! The following fields are allocated: 

        ! Note: roughly the same conventions on data structures. fields
        ! 


        ! The usual
        implicit none

        ! Declare variables
        type(CellUDT)       :: cell

        ! Allocate
        !=========
        ! Guard cells
        allocate(cell%GC(cell%ntot))

        ! Vertex values
        allocate(cell%vertP(cell%ntot,2))
        allocate(cell%vert(cell%nvert)) 

        ! Face values
        allocate(cell%faceP(cell%ntot,2))
        allocate(cell%face(cell%nface))

        ! Others
        allocate(cell%psi(cell%ntot), cell%bp(cell%ntot), &
            cell%bt(cell%ntot), cell%x(cell%ntot), cell%y(cell%ntot), &
            cell%cflags(cell%ntot), cell%reg(cell%ntot), cell%ft(cell%ntot))

        ! Initialize
        cell%psi    = 0
        cell%bp     = 0
        cell%bt     = 0
        cell%x      = 0
        cell%y      = 0
        cell%cflags = 0
        cell%reg    = 0
        cell%ft     = 0

    end subroutine

    ! Boundary substructure
    subroutine AllocateBnd(bnd)

        ! Description
        !============
        ! Allocate the boundary fields. It is assumed that the following
        ! fields are present:
        !
        ! - nface:          number of boundary faces
        !
        ! Note: the amount of vertices is per definition the amount of
        ! faces plus one, also for closed boundaries. In the latter 
        ! case, the begin and end vertex index are the same. This 
        ! assumes that the boundaries should be polygons that either 
        ! close perfectly on themselves, or are simple polygons (i.e. 
        ! no polygon vertex should have more than three faces it belongs
        ! to)
        
        ! The usual
        implicit none

        ! Declare variables
        type(BndUDT)       :: bnd

        ! Allocate
        !=========
        bnd%nvert = bnd%nface+1
        allocate(bnd%face(bnd%nface))
        allocate(bnd%vert(bnd%nvert))
        

    end subroutine
 

    ! Additional grid data structures
    !================================
    ! Data substructure
    subroutine AllocateGridData(data,grid)

        ! Description
        !============
        ! Allocate the structures inside the grid data structures. It 
        ! is assumed that the main grid structures (cells/faces/vert) 
        ! are already properly initialized. 

        ! The usual
        implicit none

        ! Declare variables
        type(GridDataUDT)       :: data
        type(GridUDT)           :: grid

        ! Flux data
        call AllocateFluxData(data%fluxdata,grid)

    end subroutine

    ! Flux data substructure
    subroutine AllocateFluxData(fluxdata,grid)

        ! Description
        !============
        ! Allocate the flux data structure. At least the following 
        ! fields have to be initialized (see definition of FluxDataUDT
        ! for explanation): 
        !
        ! - nFs
        ! - nFt
        ! - grid%cell%ntot
        ! - grid%face%ntot

        ! The usual
        implicit none

        ! Declare variables
        type(FluxDataUDT)       :: fluxdata
        type(GridUDT)           :: grid

        ! Allocate
        !=========
        ! Flux tube data
        allocate(fluxdata%fluxtubecellsP(fluxdata%nFt,2))
        allocate(fluxdata%fluxtubecells(grid%cell%ntot)) 
        allocate(fluxdata%fluxtubefacesP(fluxdata%nFt,2))
        allocate(fluxdata%fluxtubefaces(grid%face%ntot)) 
        allocate(fluxdata%fluxtuberegID(fluxdata%nFt))

        ! Flux surface data
        allocate(fluxdata%fluxsurfacefacesP(fluxdata%nFs,2))
        allocate(fluxdata%fluxsurfacefaces(grid%face%ntot))
        allocate(fluxdata%fluxsurfacevertsP(fluxdata%nFs,2))
        allocate(fluxdata%fluxsurfaceverts(grid%vert%ntot))
        allocate(fluxdata%fluxsurfaceID(grid%vert%ntot))
        allocate(fluxdata%fluxsurfacepsi(fluxdata%nFs))

    end subroutine

    ! Deallocation
    !=============
    ! Vertex substructure
    subroutine DeallocateVertices(vert)

        ! Description
        !============
        ! Deallocate the fields in the vert structure

        ! The usual
        implicit none

        ! Declare variables
        type(VertexUDT)       :: vert

        ! Deallocate
        !=========
        ! Vertex data
        deallocate(vert%x)
        deallocate(vert%y)
        deallocate(vert%BV)
        deallocate(vert%fieldlineID)

        ! Face data
        deallocate(vert%faceP)
        deallocate(vert%face) 

        ! Cell data
        deallocate(vert%cellP)
        deallocate(vert%cell)

        ! Neighbour data
        deallocate(vert%neigP)
        deallocate(vert%neig)

    end subroutine

    ! Face structure
    subroutine DeallocateFaces(face)

        ! Description
        !============
        ! Deallocate the fields in the face structure

        ! The usual
        implicit none

        ! Declare variables
        type(FaceUDT)       :: face

        ! Deallocate
        !===========
        deallocate(face%vert)
        deallocate(face%cell)
        deallocate(face%cellP)
        deallocate(face%label)
        deallocate(face%reg)
        deallocate(face%aligned)
        deallocate(face%TMfacelabel)
        if (allocated(face%cellP)) deallocate(face%cellP)
        if (allocated(face%BF)) deallocate(face%BF)


    end subroutine

    ! Cell structure
    subroutine DeallocateCells(cell)

        ! Description
        !============
        ! Deallocate the fields in the cell structure

        ! The usual
        implicit none

        ! Declare variables
        type(CellUDT)       :: cell

        ! Deallocate
        !===========
        deallocate(cell%vertP)
        deallocate(cell%vert)
        deallocate(cell%faceP)
        deallocate(cell%face)
        deallocate(cell%GC)
        deallocate(cell%psi)
        deallocate(cell%bp)
        deallocate(cell%bt)
        deallocate(cell%x)
        deallocate(cell%y)
        deallocate(cell%cflags)
        deallocate(cell%reg)
        deallocate(cell%ft)

    end subroutine

    ! Grid data structure
    subroutine DeallocateGridData(data)

        ! Description
        !============
        ! Deallocate the fields in the grid data structure

        ! The usual
        implicit none

        ! Declare variables
        type(GridDataUDT)       :: data

        ! Deallocate
        !===========
        if (allocated(data%OMPcell))        deallocate(data%OMPcell)
        if (allocated(data%OMPface))        deallocate(data%OMPface)
        if (allocated(data%IMPcell))        deallocate(data%IMPcell)
        if (allocated(data%IMPface))        deallocate(data%IMPface)
        if (allocated(data%xpointID))       deallocate(data%xpointID)
        if (allocated(data%spointID))       deallocate(data%spointID)
        if (allocated(data%opointID))       deallocate(data%opointID)
        if (allocated(data%tpointID))       deallocate(data%tpointID)
        if (allocated(data%isprimaryxp))    deallocate(data%isprimaryxp)
        if (allocated(data%spointxpID))     deallocate(data%spointxpID)
        if (allocated(data%divFc))          deallocate(data%divFc)
        if (allocated(data%spointdivID))    deallocate(data%spointdivID)
        if (allocated(data%tpointdivID))    deallocate(data%tpointdivID)
        if (allocated(data%divFcP))         deallocate(data%divFcP)
        if (allocated(data%sepID))          deallocate(data%sepID)
        
        call DeallocateFluxData(data%fluxdata)

    end subroutine

    subroutine DeallocateFluxData(fluxdata)

        ! Description
        !============
        ! Deallocate the fields in the flux data structure

        ! The usual
        implicit none

        ! Declare variables
        type(FluxDataUDT)       :: fluxdata

        ! Deallocate
        !===========        
        if (allocated(fluxdata%fluxtubecellsP))     deallocate(fluxdata%fluxtubecellsP)
        if (allocated(fluxdata%fluxtubecells))      deallocate(fluxdata%fluxtubecells)
        if (allocated(fluxdata%fluxtubefacesP))     deallocate(fluxdata%fluxtubefacesP)
        if (allocated(fluxdata%fluxtubefaces))      deallocate(fluxdata%fluxtubefaces)
        if (allocated(fluxdata%fluxtubefsIDs))      deallocate(fluxdata%fluxtubefsIDs)
        if (allocated(fluxdata%fluxtuberegID))      deallocate(fluxdata%fluxtuberegID)
        if (allocated(fluxdata%isclosedft))         deallocate(fluxdata%isclosedft)
        if (allocated(fluxdata%fluxsurfacefacesP))  deallocate(fluxdata%fluxsurfacefacesP)
        if (allocated(fluxdata%fluxsurfacefaces))   deallocate(fluxdata%fluxsurfacefaces)
        if (allocated(fluxdata%fluxsurfaceID))      deallocate(fluxdata%fluxsurfaceID)
        if (allocated(fluxdata%fluxsurfaceneig))    deallocate(fluxdata%fluxsurfaceneig)
        if (allocated(fluxdata%fluxsurfaceneigP))   deallocate(fluxdata%fluxsurfaceneigP)
        if (allocated(fluxdata%fluxsurfacepsi))     deallocate(fluxdata%fluxsurfacepsi)
        if (allocated(fluxdata%fluxsurfacevertsP))  deallocate(fluxdata%fluxsurfacevertsP)
        if (allocated(fluxdata%fluxsurfaceverts))   deallocate(fluxdata%fluxsurfaceverts)

    end subroutine

    ! Grid structure
    subroutine DeallocateGrid(grid)

        ! Description
        !============
        ! Deallocate the fields in the grid structure

        ! The usual
        implicit none

        ! Declare variables
        type(GridUDT)       :: grid 
        
        ! Deallocate substructures
        call DeallocateVertices(grid%vert)
        call DeallocateFaces(grid%face)
        call DeallocateCells(grid%cell)
        call DeallocateGridData(grid%data)
        deallocate(grid%bnd)

    end subroutine

    ! Basic operations
    !=================
    ! Basic pointer getter
    function GetArrayFromPointer(p1, p2, list, i) result(array)
        !  Description
        !=============
        ! General array getter from pointer input (only for integers)

        ! Declare variables
        !==================
        integer(I8), intent(in)                 :: p1(:), p2(:), list(:)
        integer(I8), intent(in)                 :: i 
        integer(I8), allocatable                :: array(:) 

        ! Get array
        array = list(p1(i):(p1(i) + p2(i) - 1))

    end function

    ! Get face of vertex set
    subroutine MapVertexPairToFace(v1,v2,fvert,nfvert,faceID)

        ! Description
        !============
        ! This routine maps a vertex pair, given by the vertex indices v1 
        ! and v2, to the face ID that has these vertices (i.e. it retrieves
        ! the face ID with the same vertices). The sequence of the vertices 
        ! doesn't matter. Caution: the routine does not check for multiple 
        ! occurences and simply returns the first face id that matches!
    
        ! Notes
        !======
        ! Note 1: altough the routine works without sorting the face vertex
        ! array fvert, performance can be improved by doing this in a pre-
        ! processing step, e.g. setting fvert(:,1) = minval(fvert,2) and
        ! fvert(:,2) = maxval(fvert,2) - this is pseudocode!
    
        ! Note 2: it is assumed that v1 and v2 are scalar integers. 
    
        ! Note 3: if no face vertex pair matches, faceID is set to -1.
    
        ! Algorithm
        !==========
        ! First, we loop and check if the pair [v1 v2] is found to be equal 
        ! to any fvert(i,:). If none are found, then we loop to see if we 
        ! find [v2 v1]. If no pairs are found, faceID is set to -1. 
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        integer(I8)                         :: v1, v2, nfvert
        integer(I8), dimension(nfvert,2)    :: fvert
        integer(I8)                         :: faceID
    
        ! Loop variables
        integer(I8)                         :: k
        logical                             :: notfound
    
        ! Auxiliary variables 
    
        ! Initialize
        !===========
        ! Default faceID output
        faceID = -1
    
        ! Iteration counter
        k = 1
    
        ! Loop variable
        notfound = .true.
    
        ! Loop
        !=====
        ! Is [v1 v2] a pair? 
        do while (notfound .and. k <= nfvert)
    
            ! Check the first vertex
            if (fvert(k,1) == v1) then 
                ! Check the second vertex
                if (fvert(k,2) == v2) then 
                    ! Found edge, exit
                    faceID = k
                    notfound = .false.
                end if
            end if
    
            ! Update counter
            k = k+1;
    
        end do
    
        ! Reset k
        k = 1
    
        ! Is [v2 v1] a pair?
        do while (notfound .and. k <= nfvert)
    
            ! Check the second vertex
            if (fvert(k,2) == v1) then 
                ! Check the first vertex
                if (fvert(k,1) == v2) then 
                    ! Found edge, exit
                    faceID = k
                    notfound = .false.
                end if
            end if
    
            ! Update counter
            k = k+1;
    
        end do
    
    
    
    
    
    end subroutine

    ! Get cells of a vertex
    function GetVertCell(vert, i) result(res)
        integer(I8)                 :: i 
        type(VertexUDT)             :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%cell(vert%cellP(i, 1):(vert%cellP(i, 1) + vert%cellP(i, 2) - 1))
    end function

    ! Get faces of a vertex
    function GetVertFace(vert, i) result(res)
        integer(I8)                 :: i 
        type(VertexUDT)             :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%face(vert%faceP(i, 1):(vert%faceP(i, 1) + vert%faceP(i, 2) - 1))
    end function

    ! Get neig of a vertex
    function GetVertNeig(vert, i) result(res)
        integer(I8)                 :: i 
        type(VertexUDT)             :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%neig(vert%neigP(i, 1):(vert%neigP(i, 1) + vert%neigP(i, 2) - 1))
    end function

    ! Get cells of a face
    function GetFaceCell(face, i) result(res)
        integer(I8)                 :: i 
        type(FaceUDT)               :: face 
        integer(I8), allocatable    :: res(:)
        res = face%cell(face%cellP(i, 1):(face%cellP(i, 1) + face%cellP(i, 2) - 1))
    end function

    ! Get faces of a cell
    function GetCellFace(cell, i) result(res)
        integer(I8)                 :: i 
        type(CellUDT)               :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%face(cell%faceP(i, 1):(cell%faceP(i, 1) + cell%faceP(i, 2) - 1))
    end function

    ! Get vertices of a cell
    function GetCellVert(cell, i) result(res)
        integer(I8)                 :: i 
        type(CellUDT)               :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%vert(cell%vertP(i, 1):(cell%vertP(i, 1) + cell%vertP(i, 2) - 1))
    end function

    ! Get faces of a flux surface
    function GetFSFace(fd, i) result(res)
        integer(I8)                 :: i 
        type(FluxDataUDT)           :: fd 
        integer(I8), allocatable    :: res(:)
        res = fd%fluxsurfacefaces(fd%fluxsurfacefacesP(i, 1):&
            (fd%fluxsurfacefacesP(i, 1) + fd%fluxsurfacefacesP(i, 2) - 1))
    end function

    ! Get cells of a flux tube
    function GetFTCell(fd, i) result(res)
        integer(I8)                 :: i 
        type(FluxDataUDT)           :: fd 
        integer(I8), allocatable    :: res(:)
        res = fd%fluxtubecells(fd%fluxtubecellsP(i, 1):&
            (fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1))
    end function

    ! Get faces of a flux tube
    function GetFTFace(fd, i) result(res)
        integer(I8)                 :: i 
        type(FluxDataUDT)           :: fd 
        integer(I8), allocatable    :: res(:)
        res = fd%fluxtubefaces(fd%fluxtubefacesP(i, 1):&
            (fd%fluxtubefacesP(i, 1) + fd%fluxtubefacesP(i, 2) - 1))
    end function

    ! Interconnections
    subroutine ComputeGridInterconnections(grid)

        ! Description
        !============
        ! Compute additional grid topology data starting from basic grid
        ! information. The following fields are required (see gdmod_types
        ! for an explanation of the different fields)
        !
        ! Vertices: ntot, BV
        ! Faces: ntot, vert
        ! Cells: ntot, vertP, vert, nvert 
        !   
        ! All other fields are recomputed based on this basic 
        ! interconnection data.             
    
        ! Notes
        !======
        ! Note 1: guard cells are identified by checking which cells only 
        ! have two vertices. It
    
        ! Note 2: it is assumed that all fields of which the size is 
        ! predetermined by the number of vertices/faces/cells are already
        ! allocated. Other fields are allocated or overwritten. 
    
        ! Note 3: the original cell vertices are replaced by a sorted list.
        ! This list is generated based on the faces of the cell and their
        ! corresponding vertices. It is assumed that each cell has the same
        ! amount of faces as it has vertices.
    
        ! Note 4: the vertex neighbours and vertex cells are ordened in 
        ! clockwise or counter-clockwise manner. This sorting assumes that
        ! each vertex has more than two cells - otherwise, it doesn't really
        ! matter how it is sorted. If it has more than two cells, the cells
        ! have to be properly connected, in the sense that each cell should
        ! only have one neighbour. 
    
        ! Initialize
        !===========
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT)               :: grid
    
        ! Loop variables
        integer(I8)                 :: i, j, k, m, q
    
        ! Auxiliary variables 
        type(VertexUDT)             :: v 
        type(FaceUDT)               :: f
        type(CellUDT)               :: c
    
        integer(I8)                 :: nc, nv, nf, sp, ep, tcs, thiscell, &
            nvc, ntv, tf, ind, tfv(1:2), tc, tcn, ncf, nextcell, tcf, vc, &
            ntcf2, nb1, nb2, nfc, fcc
        integer(I8), allocatable    :: fcount(:), vcount(:), tv(:), &
            tempfcell(:, :), localID(:), fc(:), fn(:), allvertcells(:), &
            allfv(:), tcf2(:), nb(:), vcc(:), indv(:)
    
        logical, allocatable        :: cellfound(:), allfvind(:)
        logical                     :: startcellnotfound, allnotfound, &
            accountforGC
    
        ! Unpack & initialize
        !==================
        ! Data structures 
        v = grid%vert
        f = grid%face
        c = grid%cell  
    
        ! Checks
        if (size(f%vert,2) /= 2) then
            ! This should actually never happen, but just to be sure
            call gdErrorHandler('ComputeGridInterconnections: too many vertices per face')
        end if
        if (any(any(f%vert == 0,2),1)) then
            ! This could be the case when the face vertices are not correctly determined
            call gdErrorHandler('ComputeGridInterconnections: some vertex indices are zero in faces, not supported')
        end if
    
        ! Quantities
        nv = v%ntot
        nf = f%ntot
        nc = c%ntot
    
        ! Counters
        allocate(fcount(nf))
        allocate(vcount(nv))
    
        ! Face neighbours, vertex cells, cell faces
        !==========================================
        ! Check allocation
        if (allocated(c%face)) then
            ! Print out warning
            print *, 'ComputeGridInterconnections: recomputing cell%face'
    
            ! Deallocate
            deallocate(c%face)
        end if
        if (allocated(v%cell)) then 
            ! print out warning
            print *, 'ComputeGridInterconnections: recomputing vert%cell'
    
            ! Deallocate
            deallocate(v%cell)
        end if
        if (allocated(f%cell)) then 
            ! Print out warning
            print *, 'ComputeGridInterconnections: recomputing face%cell'
    
            ! Deallocate
            deallocate(f%cell)
        end if
        if (allocated(v%cell)) then 
            ! Print out warning
            print *, 'ComputeGridInterconnections: recomputing vert%cell'
    
            ! Deallocate
            deallocate(v%cell)
        end if
    
        ! Check which cells are guard cells
        c%GC(:) = .false. ! initialize to false
        do i = 1, nc 
            if (c%vertP(i,2) <= 2) then
                ! Guard cell found, set true
                c%GC(i) = .true.
            end if
        end do
    
        ! Check if we have to account for guard cells
        accountforGC = .false.
        if (any(c%GC)) then
            accountforGC = .true.
        end if 
    
        ! Make the pointer list for cell faces
        c%faceP = 0
        c%faceP(1,1) = 1
        if (c%GC(1)) then 
            ! Only one face
            c%faceP(1,2) = 1
        else
            ! Same amount of faces as vertices 
            c%faceP(1,2) = c%vertP(1,2) 
        end if
    
        do i = 2, nc
            ! Pointer
            c%faceP(i,1) = c%faceP(i-1,1) + c%faceP(i-1,2)
    
            ! Number of faces
            if (c%GC(i)) then 
                c%faceP(i,2) = 1
            else
                c%faceP(i,2) = c%vertP(i,2)
            end if
    
        end do
    
        ! Compute the total faces in cells%face
        c%nface = sum(c%faceP(:,2)) 
    
        ! Allocate
        allocate(c%face(c%nface))
        c%face = 0
    
        ! Make the pointer list for vertex cell neighbours
        v%cellP = 0
        do i = 1, nc
            ! Get the number of vertices of the current cell
            ntv = c%vertP(i,2)
    
            ! Allocate
            allocate(tv(ntv))
    
            ! Get the vertices of the current cell
            tv = c%vert(c%vertP(i,1):(c%vertP(i,1)+ntv-1))
    
            ! Update the counter
            v%cellP(tv,2) = v%cellP(tv,2)+1
    
            ! Deallocate
            deallocate(tv)
    
        end do
    
        ! Compute the total number of cells in vert%cell
        v%ncell = sum(v%cellP(:,2))
    
        ! Allocate
        allocate(v%cell(v%ncell))
        v%cell = 0
    
        ! Construct the pointer
        v%cellP(1,1) = 1
        do i = 2, nv 
            v%cellP(i,1) = v%cellP(i-1,1) + v%cellP(i-1,2)
        end do
    
        ! Set vertex cell and cell face and face neighbours
        fcount = 1
        vcount = 1
        allocate(tempfcell(f%ntot, 2))
        tempfcell = 0
        !!$omp parallel do default(none) schedule(guided) if(.not. omp_in_parallel()) &
        !!$omp private(i, j, k, tv, ntv, tf, ind) &
        !!$omp shared(v, f, c, nv, nf, nc, fcount, tempfcell, vcount) 
        do i = 1, nc ! loop over all cells
            ! Get vertices of this cell
            tv = GetCellVert(c, i)
            ntv = c%vertP(i,2)
    
            ! Loop over all vertices
            do j = 1, ntv
    
                ! Get the faces of the current cell
                ! Get the next face 
                if (j == ntv) then
                    k = 1
                else
                    k = j+1
                end if 
    
                call MapVertexPairToFace(tv(j), tv(k), f%vert, nf, tf)
    
                ! Check
                if (tf < 0) then
                    ! This shouldn't happen
                    print *, 'vertex pair: ', tv(j), tv(k)
                    call gdErrorHandler(&
                        'ComputeGridInterconnections: could not find '&
                        // 'matching face for this vertex pair')
                else
    
                    ! Add the current cell index as neighbour to this face
                    if (.not. ((j > 1) .and. c%GC(i))) then ! hedge for guard cells
                        ! Sanity check
                        if (fcount(tf) > 2) then
                            print *, 'face ID: ', tf
                            print *, 'neighbours: ', tempfcell(tf,:), i
                            print *, 'vert: ', f%vert(tf, 1), f%vert(tf, 2)
                            call gdErrorHandler(& 
                            'ComputeGridInterconnections: too many ' &
                                // 'neighbours for this face')
                        end if
                    
                        !$omp critical
                        tempfcell(tf,fcount(tf)) = i
    
                        ! Update fcount
                        fcount(tf) = fcount(tf)+1
                        !$omp end critical
                    end if
    
                    ! Add the current cell to the j'th vertex
                    !$omp critical
                    ind = v%cellP(tv(j),1) + vcount(tv(j)) - 1
                    if (ind > v%ncell) then
                        call gdErrorHandler('unknown error')
                    end if
                    v%cell(ind) = i
    
                    ! Update vcount
                    vcount(tv(j)) = vcount(tv(j))+1
                    
                    ! Add cell face
                    if (.not. ((j > 1) .and. c%GC(i))) then ! hedge for guard cells
                        ind = c%faceP(i,1) + j - 1
                        if (ind > c%nface) then
                            call gdErrorHandler('unknown error')
                        end if
                        c%face(ind) = tf
                    end if
                    !$omp end critical
                end if
    
            end do
    
            ! Deallocate
            deallocate(tv)
    
        end do
        !!$omp end parallel do
    
        ! Construct cell arrays for faces and vertices
        fcount = fcount-1
        f%cellP(:, 2) = fcount 
        f%cellP(1, 1) = 1
        f%ncell = sum(f%cellP(:, 2))
        allocate(f%cell(f%ncell))
        fcc = 0
        nfc = fcount(1)
        f%cell(fcc+1:fcc+nfc) = tempfcell(1, 1:nfc)
        fcc = fcc + nfc
        do i = 2, f%ntot
            f%cellP(i, 1) = f%cellP(i-1, 2) + f%cellP(i-1, 1)
            nfc = fcount(i)
            f%cell(fcc+1:fcc+nfc) = tempfcell(i, 1:nfc)
            fcc = fcc + nfc
        end do
    
        ! Housekeeping
        deallocate(vcount, fcount)
        
        ! Logicals
        !=========
        ! Based on the interconnections computed above, the logical indices
        ! (boundary vertices => BV, boundary faces => BF) etc can be 
        ! computed. Boundary faces are either: 1) faces with a single 
        ! neighbour (if no guard cells are present), 2) faces with a guard
        ! cell neighbour (if guard cells are present)
    
        ! Compute BF & BV
        f%BF(:) = .false. ! initialize
        v%BV(:) = .false. ! initialize
        do i = 1, f%ntot
    
            ! Get face neighbours
            nb  = GetFaceCell(f, i)
            nb1 = 0 ! Initialize to be sure
            nb2 = 0
            if (size(nb, 1) == 1) then 
                nb1 = nb(1)
                nb2 = 0
            elseif (size(nb, 1) == 2) then   
                nb1 = nb(1)
                nb2 = nb(2)
            else 
                ! This shouldn't happen
                call gdErrorHandler('More than two cells detected as face neighbours, check grid interconnection')
            end if 
            !deallocate(nb)
    
            ! Sanity checks
            if ( (nb1 == 0) .or. (accountforGC .and. (nb2 == 0)) ) then 
                ! Call error, this should not be the case
                print *, 'face ID: ', i
                call gdErrorHandler('ComputeGridInterconnections:' &
                // 'could not find neighbours for this face')
            end if
    
            ! Set BF & BV logical
            if (accountforGC) then
                if (c%GC(nb2) .or. c%GC(nb1) ) then
                    ! This is a boundary face
                    f%BF(i) = .true.
                    v%BV(f%vert(i,:)) = .true.
                end if
            else
                if (nb2 == 0) then
                    ! This is a boundary face
                    f%BF(i) = .true.
                    v%BV(f%vert(i,:)) = .true.
                end if
            end if
        end do
    
    
        ! Additional vertex interconnections
        !===================================
        ! Vertex neighbours
        ! Check allocation
        if (.not. allocated(v%neigP)) then
            allocate(v%neigP(v%ntot,2))
        end if
    
        ! Compute how much neighbours each vertex has by looping over faces
        v%neigP = 0
        do i = 1, f%ntot
            ! Get face vertices
            tfv = f%vert(i,:)
    
            ! Update counters
            v%neigP(tfv(1),2) = v%neigP(tfv(1),2) + 1
            v%neigP(tfv(2),2) = v%neigP(tfv(2),2) + 1
        end do
    
        ! Compute the first column of neigP
        v%neigP(1,1) = 1
        do i = 2, v%ntot
            v%neigP(i,1) = v%neigP(i-1,1) + v%neigP(i-1,2)
        end do
    
        ! Check the neiglist allocation
        f%ncell = sum(v%neigP(:,2))
        if (allocated(v%neig)) then 
            ! This shouldn't be the case, but we hedge for it. Reallocate
            deallocate(v%neig)
            allocate(v%neig(f%ncell))
        else
            allocate(v%neig(f%ncell))
        end if
        v%neig = 0
    
        ! Compute the neiglist
        allocate(vcount(v%ntot))
        vcount = 0
        do i = 1, f%ntot
            ! Get the face vertices
            tfv = f%vert(i,:)
    
            ! Add the neighbours
            sp = v%neigP(tfv(1),1) + vcount(tfv(1))
            v%neig(sp) = tfv(2)
            sp = v%neigP(tfv(2),1) + vcount(tfv(2))
            v%neig(sp) = tfv(1)
    
            ! Update counter
            vcount(tfv) = vcount(tfv) + 1
    
        end do
        deallocate(vcount)
    
        ! Compute vertex faces
        v%faceP = 0
        do i = 1, f%ntot
            v%faceP(f%vert(i, :), 2) = v%faceP(f%vert(i, :), 2) + 1;
        end do
        v%nface = sum(v%faceP(:, 2))
        v%faceP(1, 1) = 1
        do i = 2, v%ntot 
            v%faceP(i, 1) = v%faceP(i-1, 1) + v%faceP(i-1, 2)
        end do 
        if (allocated(v%face)) then 
            deallocate(v%face)
        end if 
        allocate(v%face(v%nface), vcc(v%ntot)) 
        vcc(:) = 0
        do i = 1, f%ntot
            tfv = f%vert(i, :)
            indv = vcc(tfv) + v%faceP(tfv, 1)
            v%face(indv) = i
            vcc(tfv) = vcc(tfv) + 1
        end do 
    
    
        ! Sort vertex neighbours and vertex cells
        !========================================
        
        ! Loop over all vertices
        !!$omp parallel do default(none) schedule(guided) if(.not. omp_in_parallel()) &
        !!$omp private (i, j, tcs, nvc, cellfound, localID, allvertcells, &
        !!$omp sp, ep, tcf, startcellnotfound, thiscell, m, tc, fc, ncf, &
        !!$omp k, tcn, fn, q, vc, allnotfound, nextcell, allfv, allfvind, &
        !!$omp tcf2, ntcf2, tfv) &
        !!$omp shared (v, f, c, accountforGC)
        do i = 1, v%ntot ! v%ntot
    
            ! Check how many distinct cell sequences there are by checking 
            ! the amount of faces of this vertex and the amount of cells. 
            ! If the difference is zero -> internal vertex (check if this is
            ! true, otherwise throw error), starting cell doesn't matter. 
            ! If one or higher -> only possible if boundary vertex 
            ! (check this). Difference indicates the amount of cell 
            ! sequences.
    
            ! Note: we have to account for guard cells here. Normally, if 
            ! guard cells are present, each boundary vertex should have a
            ! multiple of two amount of guard cells (more than 2 only holds
            ! for boundary vertices that connect two domains in a pointwise
            ! way (or more domains)).
    
            ! Number of cell sequences
            tcs = 0
            tcs = v%neigP(i, 2) - v%cellP(i, 2)
    
            ! Account for guard cells, if any
            if (accountforGC) then 
                tcs = tcs + count(c%GC(v%cell(v%cellP(i, 1):v%cellP(i, 2)+v%cellP(i, 1)-1)))
            end if
    
            ! Check
            if ( (tcs == 0) .and. (v%BV(i)) ) then 
                ! Throw error - this vertex should be an internal vertex,
                ! not a boundary one. Check the grid interconnection. 
                print *, 'vertex ID: ', i
                call gdErrorHandler(&
                    'ComputeGridInterconnections: this vertex is ' &
                    // 'falsely classified as boundary vertex')
            end if
            if ( (tcs > 0) .and. (.not. v%BV(i)) ) then 
                ! Throw error - this vertex should be a boundary vertex.
                ! Check the grid interconnection. 
                print *, 'vertex ID: ', i
                call gdErrorHandler(&
                    'ComputeGridInterconnections: this vertex is ' &
                    // 'falsely classified as internal vertex')
            end if
    
            ! Issue warning when tcs is higher than 1 - only occurs if two
            ! boundaries touch each other in one point, which should hardly
            ! ever be the case. 
            if (tcs > 1) then 
                print *, 'WARNING: boundary vertex ', i, &
                    ' found with multiple cell sequences - check interconnection'
            end if 
    
            ! Set to one if equal to zero (for internal vertices)
            if (tcs == 0) then
                tcs = 1
            end if
    
            ! Loop over all cell sequences
            nvc = v%cellP(i,2)
            allocate(cellfound(nvc))
            allocate(localID(nvc))
            allocate(allvertcells(nvc))
            cellfound(:) = .false.
            sp = v%cellP(i,1)
            ep = v%cellP(i,1) + nvc-1
            allvertcells = v%cell(sp:ep)
    
            tcf = 0
            do j = 1, tcs
    
                ! Get the first cell and its faces
                startcellnotfound = .true.
                thiscell = 0
                if (.not. v%BV(i)) then 
                    ! Internal cell, simply take the first cell, since only 
                    ! one sequence
                    sp = v%cellP(i,1)
                    thiscell = v%cell(sp)
    
    
                    ! Update 
                    cellfound(1) = .true.
                    startcellnotfound = .false.
    
                else
                    ! Boundary cell, search for cell with only one (or none)
                    ! faces in common with the other cells
    
                    ! Find cell with one or zero common faces with the other
                    ! cells
                    
                    m = 1
                    do while (startcellnotfound .and. (m <= nvc))
    
                        ! Get the current cell
                        if (cellfound(m)) then
                            ! Update m and skip the rest
                            m = m + 1
                            cycle
                        end if
    
                        ! Get the faces of the current cell
                        ! tc = allremcells(m) ! current cell
                        tc = allvertcells(m)
                        fc = GetCellFace(c, tc)
                        !sp = c%faceP(tc,1)
                        !ep = c%faceP(tc,1) + c%faceP(tc,2) - 1 
                        !allocate(fc(c%faceP(tc,2)))
                        !fc = c%face(sp:ep)
    
                        ! Check how many common faces there are with faces
                        ! of the other remaining cells
                        ncf = 0 ! number of common faces
                        do k = 1, nvc
                            ! Skip if k == m or if cell is already found
                            if ((k == m) .or. cellfound(k)) then 
                                cycle 
                            end if
    
                            ! Get the faces of the next cell
                            tcn = allvertcells(k)
                            sp = c%faceP(tcn,1)
                            ep = c%faceP(tcn,1) + c%faceP(tcn,2) - 1 
                            allocate(fn(c%faceP(tcn,2)))
                            fn = c%face(sp:ep)
    
                            ! Count how many faces they have in common
                            do q = 1, c%faceP(tc,2)
                                ncf = ncf + count(fc(q) == fn)
                            end do
    
                            ! Deallocate
                            deallocate(fn)
    
                        end do
    
                        ! Deallocate
                        deallocate(fc)
    
                        ! Check how many common faces there are
                        !print *, ncf
                        if (ncf < 2) then 
                            ! Cell found
                            startcellnotfound = .false.
    
                            ! Set current cell
                            sp = v%cellP(i,1) + m - 1
                            thiscell = v%cell(sp)
    
                            ! Update logicals
                            ! cellfound(thislocalID(m)) = .true.
                            cellfound(m) = .true.
    
                            ! Set the remaining cells
                            !allocate(tcrem(nremcells-1))
                            !tcrem = allremcells((/(q, q = i, m-1),(q, q = m+1,nremcells)/))
    
    
                        else
                            ! Update m
                            m = m + 1
    
                        end if
    
                    end do
    
                    ! Check
                    if (startcellnotfound) then
                        ! Should have found starting cell, throw error. 
                        ! Something is wrong then with the interconnection,
                        ! and it could be that this boundary vertex should 
                        ! in fact not be a boundary vertex. 
                        print *, 'vertex ID: ', i 
                        call gdErrorHandler('ComputeGridInterconnections: could not find starting cell for this boundary vertex')
    
                    end if 
    
    
                end if
    
                ! Find the sequence
                !==================
                ! Now that the first cell of the sequence has been found, we
                ! need to extract the sequence itself. 
    
                ! Loop over the remaining cells
                k = 1
                vc = 0 ! vertex counter
                allnotfound = .true. 
                
                do while (k <= nvc)
    
                    ! Set the faces of the current cell
                    allocate(fc(c%faceP(thiscell,2)))
                    sp = c%faceP(thiscell, 1)
                    ep = sp + c%faceP(thiscell,2)-1
                    fc = c%face(sp:ep)
    
                    ! Take the next cell
                    if (cellfound(k)) then
                        ! Update k, cycle
                        k = k + 1
                        deallocate(fc)
                        cycle
                    else
                        nextcell = allvertcells(k)
                    end if
                    
                    ! Get cell faces
                    allocate(fn(c%faceP(nextcell,2)))
                    sp = c%faceP(nextcell, 1)
                    ep = sp + c%faceP(nextcell,2)-1
                    fn = c%face(sp:ep)
    
                    ! Check for common faces - should only be one
                    ncf = 0
                    do q = 1, c%faceP(thiscell, 2)
                        ncf = ncf + count(fc(q) == fn)
                        if ((count(fc(q) == fn)) == 1) then
                            ! Set the current face
                            tcf = fc(q)
                        end if
                    end do
    
                    ! Sanity check
                    if (ncf > 1) then
                        ! Throw error: the next cell is connected with too
                        ! many faces to the previous one, which shouldn't be
                        ! possible in any decent grid. Though in some 
                        ! exceptional cases, this may be allowed, we don't 
                        ! support that here. 
                        call gdErrorHandler('ComputeGridInterconnections: '&
                            // 'too many common faces detected when ' &
                            // 'sorting vertex neighbours, check grid ' &
                            // 'interconnection')
                    end if
    
                    ! Check if a common face was found
                    if (ncf == 1) then 
                        ! Add the first vertex if this is the first found 
                        ! cell
                        if (vc == 0) then
                            ! Get other face with this vertex
    
                            ! Get the faces that have the current vertex
                            if (c%GC(thiscell)) then
                                ! Hedge for guard cells: only one face
                                allocate(allfv(1))
                                allocate(allfvind(1))
                                allocate(tcf2(1))
                                tcf2(1) = fc(1)
                            else
                                allocate(allfvind(c%faceP(thiscell,2)))
                                allfvind = (f%vert(fc, 1) == i) .or. &
                                    (f%vert(fc, 2) == i)
                                allocate(allfv(count(allfvind)))
                                allfv = pack(fc, allfvind)
                                allocate(tcf2(count(allfvind)))
                                tcf2 = pack(allfv, allfv .ne. tcf)
                            end if
                            
                            ! Check 
                            ntcf2 = size(tcf2)
                            if (ntcf2 .ne. 1) then
                                call gdErrorHandler( &
                                'ComputeGridInterconnections: could not ' &
                                // 'find a single common face, check grid' &
                                //' interconnections')
                            end if
    
                            ! Add cell and vertex neighbour - neighbour only
                            ! if the current cell is not a guard cell
                            
                            sp = v%cellP(i,1)
                            !$omp critical
                            v%cell(sp) = thiscell 
                            if (.not. c%GC(thiscell)) then
                                sp = v%neigP(i,1)
                                tfv = f%vert(tcf2(1),:)
                                if (tfv(1) == i) then
                                    v%neig(sp) = tfv(2)
                                else
                                    v%neig(sp) = tfv(1)
                                end if
    
                                ! Update vc
                                vc = vc+1
                            end if
                            !$omp end critical
    
                            ! Housekeeping
                            deallocate(allfvind)
                            deallocate(allfv)
                            deallocate(tcf2)
    
                        end if 
    
                        ! Add the cell and vertex neighbour
                        sp = v%cellP(i,1) + vc
                        !$omp critical
                        v%cell(sp) = nextcell
                        sp = v%neigP(i,1) + vc
                        tfv = f%vert(tcf,:)
                        if (tfv(1) == i) then 
                            v%neig(sp) = tfv(2)
                        else
                            v%neig(sp) = tfv(1)
                        end if
                        !$omp end critical
    
                        ! Update counter
                        vc = vc + 1
    
                        ! Update logicals
                        cellfound(k) = .true.
    
                        ! Reset iterator
                        k = 1
    
                        ! Set current cell
                        thiscell = nextcell
    
                    else
                        ! Update the counter
                        k = k + 1
                    end if
    
                    ! Housekeeping
                    deallocate(fc, fn)
    
                end do
    
                ! Add last vertex
                !================
                ! Only in case of boundary vertices
                if ( (v%BV(i)) .and. (v%cellP(i,2) > 1) ) then
    
                    ! First, update vc - is now one-off
                    ! vc = vc - 1
    
                    ! Sanity check
                    if (vc > v%neigP(i,2)) then
                        ! Throw error. This can happen if the supposed 
                        ! boundary vertex is no actual boundary vertex.
                        call gdErrorHandler('ComputeGridInterconnections: ' &
                            // 'supposed boundary vertex may be internal ' &
                            // 'vertex. Check grid interconnection')
                    end if
    
                    ! Set the faces of the current cell
                    allocate(fc(c%faceP(thiscell,2)))
                    sp = c%faceP(thiscell, 1)
                    ep = sp + c%faceP(thiscell,2)-1
                    fc = c%face(sp:ep)
    
                    ! Add the vertex neighbour
                    if (c%GC(thiscell)) then
                        ! Hedge for guard cells: only one face
                        allocate(tcf2(1))
                        allocate(allfv(1))
                        allocate(allfvind(1))
                        tcf2(1) = fc(1)
                    else
                        allocate(allfvind(c%faceP(thiscell,2)))
                        allfvind = (f%vert(fc, 1) == i) .or. &
                            (f%vert(fc, 2) == i)
                        allocate(allfv(count(allfvind)))
                        allfv = pack(fc, allfvind)
                        tcf2 = pack(allfv, allfv .ne. tcf)
    
                    end if
                    
                    ! Check 
                    ntcf2 = size(tcf2)
                    if (ntcf2 .ne. 1) then
                        call gdErrorHandler( &
                        'ComputeGridInterconnections: could not ' &
                        // 'find a single common face, check grid' &
                        //' interconnections')
                    end if
    
                    ! Add cell and vertex neighbour
                    sp = v%neigP(i,1) + vc
                    tfv = f%vert(tcf2(1),:)
                    !$omp critical
                    if (tfv(1) == i) then
                        v%neig(sp) = tfv(2)
                    else
                        v%neig(sp) = tfv(1)
                    end if
                    !$omp end critical
    
                    ! Housekeeping
                    deallocate(allfvind)
                    deallocate(allfv)
                    deallocate(fc, tcf2)
                    
                end if
            end do
    
            ! Sanity check
            if (any(.not. cellfound)) then
                ! Not all cells were ordened, throw error
                print *, 'vertex ID: ', i
                call gdErrorHandler('ComputeGridInterconnections: could not sort all cells, check grid interconnection')
            end if 
    
            ! Housekeeping
            deallocate(cellfound)
            deallocate(localID)
            deallocate(allvertcells)
    
        end do
        !!$omp end parallel do 
    
        ! Add to grid
        !============
        grid%cell   = c 
        grid%face   = f 
        grid%vert    = v
    
    
    end subroutine


    !------------------------------------------------------------------!
    !                             Numerics                             !
    !------------------------------------------------------------------!

    ! Grid data updating
    subroutine UpdateGridData(grid, magneticField, environment)

        ! Description
        !============
        ! Update grid data such as metrics, magnetic field at vertices,
        ! ... that depends on the magnetic field and possibly on 
        ! the environment. Assumes all other fields are properly 
        ! initialized and allocated. 

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: grid 
        type(MagneticFieldUDT), intent(in)  :: magneticField 
        type(EnvironmentUDT), intent(in)    :: environment

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tv
        real(R8), allocatable, dimension(:)     :: tpsi
        
        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Associate
        associate(mf        => magneticField%interp,    &
                nc          => grid%cell%ntot,          &
                nfs         => grid%data%fluxdata%nfs)

        ! Compute
        !========
        ! Magnetic field at vertices
        call mf%Evaluate(grid%vert%x, grid%vert%y, 0, 0, grid%vert%psi)
        call mf%Evaluate(grid%vert%x, grid%vert%y, 1, 0, grid%vert%bx) 
        call mf%Evaluate(grid%vert%x, grid%vert%y, 0, 1, grid%vert%by) 

        ! Cell centers
        do i = 1, nc 
            tv = GetCellVert(grid%cell, i)
            grid%cell%x(i) = sum(grid%vert%x(tv))/size(tv)
            grid%cell%y(i) = sum(grid%vert%y(tv))/size(tv)
        end do 

        ! Amount of guard cells (not present here, simply amount of boundary faces)
        grid%cell%ngc = count(grid%face%BF)

        ! Magnetic field at cell centers
        call mf%Evaluate(grid%cell%x, grid%cell%y, 0, 0, grid%cell%psi)
        do i = 1, nc
            tv = GetCellVert(grid%cell, i)
            grid%cell%bp(i) = -sum(sqrt(grid%vert%bx(tv)**2 + grid%vert%by(tv)**2))/(size(tv)*2*pi_R8*grid%cell%x(i))
            grid%cell%bt(i) = sum(grid%vert%ffbz(tv))/(size(tv)*2*pi_R8*grid%cell%x(i))
        end do

        ! Flux surface psi values
        do i = 1, nfs 
            ! Get psi values
            allocate(tpsi(count(grid%vert%fieldlineID == i)))
            tpsi = pack(grid%vert%psi, grid%vert%fieldlineID == i)

            ! Average and add
            grid%data%fluxdata%fluxsurfacepsi(i) = sum(tpsi)/size(tpsi)
            deallocate(tpsi)

        end do

        ! Housekeeping
        !=============
        end associate


    end subroutine

    !------------------------------------------------------------------!
    !                          Magnetic field                          !
    !------------------------------------------------------------------!

    ! Magnetic field
    !===============
    ! Allocation
    subroutine AllocateMagneticField(magneticField)

        ! Description
        !============
        ! Allocate the magnetic field properties. At least the fields nR
        ! and nZ should be present. 

        ! The usual
        implicit none

        ! Declare variables
        type(MagneticFieldUDT)          :: magneticField

        ! Allocate
        !=========
        allocate(magneticField%R(magneticField%nR))
        allocate(magneticField%Z(magneticField%nZ))
        allocate(magneticField%Psi(magneticField%nR, magneticField%nZ))

    end subroutine

    ! Constructor
    subroutine ConstructMagneticField(magneticField, mfoptions)

        ! Description
        !============
        ! Construct a representation of the magnetic field based on given
        ! magnetic field flux data (or other data). It is assumed that the 
        ! basic magnetic field data is already available in the 
        ! magneticField structure
    
        ! Notes
        !======
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        type(MagneticFieldOptionsUDT)       :: mfoptions 
        type(MagneticFieldUDT)              :: magneticField
    
        ! Auxiliary 
        real(R8)                                :: dR, dZ 
        real(R8), allocatable, dimension(:)     :: x, y, temp
        real(R8), allocatable, dimension(: ,:)  :: val, tempval  
        type(StructuredInterpolant2DUDT)    :: interp
    
        ! Loop
        integer(I8)                         :: i, k 
    
    
        ! Main program
        !=============
        ! Check if we want to reinterpolate on a different resolution
        if (mfoptions%reinterpolate) then 
            ! Compute field sizes
            dR = magneticField%R(magneticField%nR) - magneticField%R(1)
            dZ = magneticField%Z(magneticField%nZ) - magneticField%Z(1)
    
            ! Compute reinterpolation vectors 
            x = [(k, k = 0, mfoptions%resx)]*dR/real(mfoptions%resx, kind=R8) + magneticField%R(1)
            x(1) = magneticField%R(1) 
            x(size(x)) = magneticField%R(magneticField%nR)
    
            y = [(k, k = 0, mfoptions%resy)]*dZ/real(mfoptions%resy, kind=R8) + magneticField%Z(1)
            y(1) = magneticField%Z(1) 
            y(size(y)) = magneticField%Z(magneticField%nZ)
    
            ! Interpolate
            val = magneticField%Psi 
            allocate(tempval(size(x), magneticField%nZ))
            do i = 1, magneticField%nZ 
                call Interpolate1D(x, temp, magneticField%R, val(:, i))
                tempval(:, i) = temp
            end do 
            val = tempval 
            deallocate(tempval)
            allocate(tempval(size(x), size(y)))
            do i = 1, size(x) 
                call Interpolate1D(y, temp, magneticField%Z, val(i, :))
                tempval(i, :) = temp
            end do 
            val = tempval 
        else
            val = magneticField%Psi 
            x = magneticField%R 
            y = magneticField%Z 
        end if 
    
        ! Construct interpolant representation
        call interp%SetParameters(mfoptions%interpmeth, mfoptions%interpC, &
            mfoptions%interpM)
        call interp%Construct(x, y, val)
    
        ! Add interpolant to the magnetic field
        magneticField%interp = interp
    
        ! Visualize
        call magneticField%interp%Visualize('magneticfield_visualization')
        call magneticField%interp%Visualize('magneticfield_dpsidx_visualization', &
            xderivin=1, yderivin=0)
        call magneticField%interp%Visualize('magneticfield_dpsidy_visualization', &
            xderivin=0, yderivin=1)
    
    end subroutine
    
    ! Main reader
    subroutine ReadMagneticField(magneticField, mfoptions, filepath)

        ! Description
        !============
        ! Read in the magnetic field data. It is assumed that the data is 
        ! given as flux values on a structured (but not necessarily
        ! equidistant) grid. 
    
        ! Notes
        !======
    
        ! Initialize
        !===========
        ! Declare modules
        use goatmod_userinput
        
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        integer                         :: filespecifier
        type(MagneticFieldOptionsUDT)   :: mfoptions 
        type(MagneticFieldUDT)          :: magneticField
        character(*)                    :: filepath
    
        ! Loop variables
    
        ! Auxiliary variables 
    
        ! Data
        data filespecifier /60/
    
        ! Main program
        !=============
        ! Open the file
        print *, 'reading magnetic field data from file: ' // filepath
        open(unit = filespecifier, file = filepath)
    
        ! Check how to read the file
        select case(mfoptions%readmeth)
    
        case ('rzpsi','default')
    
            ! Read in a classic rzpsi file 
            call Readrzpsi(filespecifier, magneticField)

            ! Add additional data - needs to be specified separately
            magneticField%RBtor = mfoptions%RBtor
    
        case ('equ')
    
            ! Read the equ file using CARRE subroutines - wrapper here
            call Readequ(filespecifier, magneticField)

            ! RBtor computed from equilibrium file
    
        case default
    
            ! Unknown reading method, throw error
            call gdErrorHandler('ReadMagneticField: unknown reading method')
    
        end select

        ! Print RBtor value
        print *, 'ReadMagneticField: RBtor [m T] = ', magneticField%RBtor
    
        
    
        ! Housekeeping
        !=============
        close(filespecifier)
    
    
    end subroutine

    ! Read from rzpsi file (assumed units in [Wb])
    subroutine Readrzpsi(filespecifier, magneticField)

        ! Description
        !============
        ! Read in the rzpsi data and add those to the magnetic field. It is 
        ! assumed that the file first gives the r coordinates, then z 
        ! coordinates, and finally the psi values in a 2D structured way. 
        ! The following format is assumed: 
        ! - first, we should encounter 'nr=xxx', where 'xxx' is an integer
        ! that specified the number of R-coordinates
        ! - then, the R-coordinates are specified as reals
    
        ! Notes
        !======
    
        ! Initialize
        !===========
        ! Declare modules
        use goatmod_userinput
        use mod_inputfileparser
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        integer                         :: filespecifier
        type(MagneticFieldUDT)          :: magneticField
    
        ! Loop variables
    
        ! Auxiliary variables 
        real(R8), allocatable       :: R(:), Z(:), psi(:, :)
    
        integer(I8)                 :: nr, nz ! real number of points in r, z directions
        integer(I8)                 :: nval 
        integer(I8), allocatable, dimension(:)  :: vali 
        
        character(:), allocatable   :: thisline ! temporary variable for line
        logical                     :: iseof
    
        ! Read R-coordinates
        !===================
        ! Initialize
        iseof = .false. 

        ! Read until we encounter 'nr'
        call ReadUntilFound(filespecifier, 'nr', iseof)
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find r-coordinates in rzpsi file')
        end if 

        ! Read nr
        backspace(filespecifier)
        call ReadSingleLine(filespecifier, thisline, iseof)
        call ReadIntegersFromString(thisline, vali, nval)
        if (nval /= 1) then 
            call gdErrorHandler('Readrzpsi: could not read number of r coordinates')
        end if 
        nr = vali(1)

        ! Read r coordinates
        allocate(R(nr))
        R = 0.0_R8
        call ReadArray(filespecifier, nr, R)
    
        ! Read Z-coordinates
        !===================
        ! Read until we encounter 'nr'
        call ReadUntilFound(filespecifier, 'nz', iseof)
        if (iseof) then 
            ! Throw error
            call gdErrorHandler('Readrzpsi: could not find z-coordinates in rzpsi file')
        end if 

        ! Read nz
        backspace(filespecifier)
        call ReadSingleLine(filespecifier, thisline, iseof)
        call ReadIntegersFromString(thisline, vali, nval)
        if (nval /= 1) then 
            call gdErrorHandler('Readrzpsi: could not read number of z coordinates')
        end if 
        nz = vali(1)

        ! Read z coordinates
        allocate(Z(nz))
        Z = 0.0_R8
        call ReadArray(filespecifier, nz, Z)

        ! Read Psi values
        !================
        ! Allocate
        allocate(psi(nr, nz))
        psi = 0.0_R8
        call ReadArray(filespecifier, nr*nz, psi, 'psi')
    
        ! Add to magnetic field
        !======================
        ! Allocate the magnetic field 
        magneticField%nR = nr
        magneticField%nZ = nz
        call AllocateMagneticField(magneticField)
    
        ! Add data
        magneticField%R = R
        magneticField%Z = Z
        magneticField%Psi = psi
    
    end subroutine

    ! Read from equ file (assumed units in [wb/rad], rescaled to [Wb])
    subroutine Readequ(filespecifier, magneticField)

        ! Description
        !============
        ! Read in the equ data and add those to the magnetic field. 
    
        ! Notes
        !======
    
        ! Initialize
        !===========
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        integer                         :: filespecifier
        type(MagneticFieldUDT)          :: magneticField
    
        ! Loop variables
    
        ! Auxiliary variables 
        real(R8)  , allocatable     :: R(:), Z(:) ! coordinates of magnetic field grid
        real(R8)                    :: btf, rtf ! toroidal field strength and radius
        real(R8)                    :: psib ! psi at boundary (?)
        integer(I8)                 :: nr, nz ! real number of points in r, z directions
        integer(I8)                 :: nval
        integer(I8), allocatable, dimension(:)  :: vali 
        real(R8), allocatable, dimension(:)     :: valr
        real(R8), allocatable       :: pfm(:, :) ! temporary psi value array
        logical                     :: reachedeof
        integer(I8)                 :: ind
        character(:), allocatable   :: thisline, valname
    
        ! Main program
        !=============
        ! Find start position for reading (i.e. jm = ...)
        reachedeof = .false. 
        do while (.true.)
            ! Read line
            call ReadSingleLine(filespecifier, thisline, reachedeof)

            ! Check for eof
            if (reachedeof) then 
                exit 
            end if 

            ! Check if it's a comment line
            ind = index(thisline, ':=')
            if (ind /= 0) then 
                ! Comment line, cycle
                cycle 
            end if 

            ! Check if it's not a comment line
            ind = index(thisline, '=')
            if (ind == 0) then 
                ! Comment or empty line, cycle
                cycle 
            end if 

            ! If we got here, we're at the start. Exit the loop
            backspace(filespecifier)
            exit 
        end do 

        ! Check
        if (reachedeof) then 
            call gdErrorHandler('Readequ: could not read in header, reached EOF')
        end if 

        ! Read from header
        valname = 'jm'
        call ReadSingleLine(filespecifier, thisline, reachedeof)
        if (index(thisline, valname) == 0) call gdErrorHandler(&
            'Readequ: ' // valname // ' not at expected location')
        call ReadIntegersFromString(thisline, vali, nval)
        if (nval /= 1) call gdErrorHandler(&
            'Readequ: could not read in value of ' // valname)
        nr = vali(1)

        valname = 'km'
        call ReadSingleLine(filespecifier, thisline, reachedeof)
        if (index(thisline, valname) == 0) call gdErrorHandler(&
            'Readequ: ' // valname // ' not at expected location')
        call ReadIntegersFromString(thisline, vali, nval)
        if (nval /= 1) call gdErrorHandler(&
            'Readequ: could not read in value of ' // valname)
        nz = vali(1)

        valname = 'psib'
        call ReadSingleLine(filespecifier, thisline, reachedeof)
        if (index(thisline, valname) == 0) call gdErrorHandler(&
            'Readequ: ' // valname // ' not at expected location')
        call ReadRealsFromString(thisline, valr, nval)
        if (nval /= 1) call gdErrorHandler(&
            'Readequ: could not read in value of ' // valname)
        psib = valr(1)

        valname = 'btf'
        call ReadSingleLine(filespecifier, thisline, reachedeof)
        if (index(thisline, valname) == 0) call gdErrorHandler(&
            'Readequ: ' // valname // ' not at expected location')
        call ReadRealsFromString(thisline, valr, nval)
        if (nval /= 1) call gdErrorHandler(&
            'Readequ: could not read in value of ' // valname)
        btf = valr(1)

        valname = 'rtf'
        call ReadSingleLine(filespecifier, thisline, reachedeof)
        if (index(thisline, valname) == 0) call gdErrorHandler(&
            'Readequ: ' // valname // ' not at expected location')
        call ReadRealsFromString(thisline, valr, nval)
        if (nval /= 1) call gdErrorHandler(&
            'Readequ: could not read in value of ' // valname)
        rtf = valr(1)

        ! Initialize
        allocate(R(nr), Z(nz), pfm(nr, nz))
        R = 0.0_R8 
        Z = 0.0_R8 
        pfm = 0.0_R8

        ! Read r values
        call ReadArray(filespecifier, nr, R, 'r(1:jm)')

        ! Read z values
        call ReadArray(filespecifier, nz, Z, 'z(1:km)')

        ! Read psi values
        call ReadArray(filespecifier, nr*nz, pfm, 'psi(j,k)')
    
        ! Allocate the magnetic field 
        magneticField%nR = nr
        magneticField%nZ = nz
        call AllocateMagneticField(magneticField)
    
        ! Add data
        magneticField%R = R(1:nr)
        magneticField%Z = Z(1:nz)
        magneticField%Psi = pfm(1:nr, 1:nz)*2*pi_R8 ! convert from Wb/rad to Wb
        magneticField%RBtor = btf*rtf
    
    end subroutine

    !------------------------------------------------------------------!
    !                               State                              !
    !------------------------------------------------------------------!

    ! State
    !======
    subroutine AllocateState(state, nc, nf, ns)

        ! Description
        !============
        ! Allocate state field

        ! Declare variables
        !==================
        ! Arguments
        type(StateUDT), intent(inout)   :: state
        integer(I8), intent(in)         :: nc, nf, ns

        if (.not.allocated(state%na)) then

            ! State variables
            state%nc = nc
            state%nf = nf
            state%ns = ns
            allocate(state%na(nc, ns))
            allocate(state%ne(nc))
            allocate(state%ua(nc, ns))
            allocate(state%uadia(nf, 2, ns))
            allocate(state%te(nc))
            allocate(state%ti(nc))
            allocate(state%tn(nc))
            allocate(state%po(nc))
            allocate(state%kt(nc))
            allocate(state%zt(nc))
            state%na = 0.0_R8
            state%ne = 0.0_R8
            state%ua = 0.0_R8
            state%uadia = 0.0_R8
            state%te = 0.0_R8
            state%ti = 0.0_R8
            state%tn = 0.0_R8
            state%po = 0.0_R8
            state%kt = 0.0_R8
            state%zt = 0.0_R8

            ! Residuals
            allocate(state%resco(nc, ns))
            allocate(state%resmo(nc, ns))
            allocate(state%resmt(nc))
            allocate(state%reshe(nc))
            allocate(state%reshi(nc))
            allocate(state%reshn(nc))
            allocate(state%respo(nc))
            allocate(state%reskt(nc))
            allocate(state%reszt(nc))
            state%resco = 0.0_R8
            state%resmo = 0.0_R8
            state%resmt = 0.0_R8
            state%reshe = 0.0_R8
            state%reshi = 0.0_R8
            state%reshn = 0.0_R8
            state%respo = 0.0_R8
            state%reskt = 0.0_R8
            state%reszt = 0.0_R8

        end if
    end subroutine
    ! Main reader
    subroutine ReadState(state, options)

        ! Description
        !============
        ! Read in the state data. 

        ! Declare variables
        !==================
        type(StateUDT)                  :: state
        type(EnvironmentOptionsUDT)     :: options 
        
        ! Loop variables
    
        ! Auxiliary variables
    
        ! Main program
        !=============
        ! Open the file
        if (options%readstate) then 
            
            ! Check how to read the file
            select case (options%readstatemeth)

            case ('b2fstate')

                ! Read b2fstate file
                call ReadB2fstate(state, options%statefilepath)

            case ('b2fplasmf')

                ! Read b2fplasmf file
                call ReadB2fplasmf(state, options%statefilepath)

            case default

                ! Unknown reading method, throw error
                call gdErrorHandler('ReadState: unknown reading method') 

            end select

        end if

    end subroutine

    ! Specific readers
    subroutine ReadB2fstate(state, filepath)

        ! Description
        !============
        ! Reading b2fstate file to extract state information

        ! Declare variables
        !==================
        ! Arguments
        type(StateUDT), intent(inout)   :: state
        character(:), allocatable       :: filepath    

        ! Auxiliary
        integer(I8) :: nc, nf, ns, filespec, idum(0:9)
        character(:), allocatable   :: chardummy   ! dummy array
        logical :: reachedeof

        ! Data
        data filespec /60/

        ! Read grid dimensions & allocate
        !================================
        ! Open the file
        print *, 'reading state from file: ' // filepath
        open(unit = filespec, file = filepath)

        ! First, read the header with the version
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadB2fstate: reached EOF prematurely')
        end if  
        
        ! Check the version to determine what to read in
        if (chardummy(8:17) >= '03.002.000') then

            ! Primary array dimensions
            call ReadArray (filespec,3,idum,'nCv,nFc,ns')
            nc = idum(0) ! note: only reading in actual cells, no guard cells
            nf = idum(1)
            ns = idum(2)

        else

            call gdErrorHandler('ReadState: structured format not supported')

        end if

        ! Allocate state based on that information
        call AllocateState(state, nc, nf, ns)

        ! Read charges - not needed
        !call ReadArray(filespec, state%zamin)

        ! Read state variables
        !---------------------
        call ReadUntilFound(filespec, 'na', reachedeof)
        if (reachedeof) then
            ! Not found, issue warning and rewind
            print *, 'ReadState: could not find field "na", setting to zero'
            rewind(filespec)
        else 
            ! Found, read state
            backspace(filespec)
            call ReadArray(filespec, nc*ns,    state%na, 'na')
            call ReadArray(filespec, nc,       state%ne, 'ne')
            call ReadArray(filespec, nc*ns,    state%ua, 'ua')
            call ReadArray(filespec, nf*2*ns,  state%uadia, 'uadia')
            call ReadArray(filespec, nc,       state%te, 'te')
            call ReadArray(filespec, nc,       state%ti, 'ti')
            call ReadArray(filespec, nc,       state%tn, 'tn')
            call ReadArray(filespec, nc,       state%po, 'po')
            call ReadArray(filespec, nc,       state%kt, 'kt')
            call ReadArray(filespec, nc,       state%zt, 'zt')
            
        end if

        ! Read residuals
        !---------------
        print *, 'ReadState: reading from b2fstate file, so no residuals'

        ! Housekeeping
        close(filespec)

    end subroutine

    subroutine ReadB2fplasmf(state, filepath)

        ! Description
        !============
        ! Read b2fplasmf file to extract state information

        ! Declare variables
        !==================
        ! Arguments
        type(StateUDT), intent(inout)   :: state
        character(:), allocatable       :: filepath 

        ! Auxiliary
        integer(I8) :: nc, nf, ns, filespec, idum(0:9)
        character(:), allocatable   :: chardummy   ! dummy array
        logical :: reachedeof

        ! Data
        data filespec /60/

        ! Read grid dimensions & allocate
        !================================
        ! Open the file
        print *, 'reading state from file: ' // filepath
        open(unit = filespec, file = filepath)

        ! First, read the header with the version
        call ReadSingleLine(filespec, chardummy, reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadB2fstate: reached EOF prematurely')
        end if  
        
        ! Check the version to determine what to read in
        if (chardummy(8:17) >= '03.002.000') then

            ! Primary array dimensions
            call ReadArray (filespec,3,idum,'nCv,nFc,ns')
            nc = idum(0) ! note: only reading in actual cells, no guard cells
            nf = idum(1)
            ns = idum(2)

        else

            call gdErrorHandler('ReadState: structured format not supported')

        end if

        ! Allocate state based on that information
        call AllocateState(state, nc, nf, ns)

        ! Read state variables
        !---------------------
        call ReadUntilFound(filespec, 'na', reachedeof)
        if (reachedeof) then
            ! Not found, issue warning and rewind
            print *, 'ReadState: could not find field "na", setting to zero'
            rewind(filespec)
        else 
            ! Found, read state
            backspace(filespec)
            call ReadArray(filespec, nc*ns,    state%na, 'na')
            call ReadArray(filespec, nc,       state%ne, 'ne')
            call ReadArray(filespec, nc*ns,    state%ua, 'ua')
            call ReadArray(filespec, nf*2*ns,  state%uadia, 'uadia')
            call ReadArray(filespec, nc,       state%te, 'te')
            call ReadArray(filespec, nc,       state%ti, 'ti')
            call ReadArray(filespec, nc,       state%tn, 'tn')
            call ReadArray(filespec, nc,       state%po, 'po')
            call ReadArray(filespec, nc,       state%kt, 'kt')
            call ReadArray(filespec, nc,       state%zt, 'zt')

        end if   
        
        ! Read residuals
        !---------------
        call ReadUntilFound(filespec, 'resco', reachedeof)
        if (reachedeof) then
            ! Not found, issue warning and rewind
            print *, 'ReadState: could not find field "resco", setting to zero'
            rewind(filespec)
        else

            ! Found, read residuals
            backspace(filespec)            
            call ReadArray(filespec, nc*ns,    state%resco, 'resco')
            call ReadArray(filespec, nc,       state%reshe, 'reshe')
            call ReadArray(filespec, nc,       state%reshi, 'reshi')
            call ReadArray(filespec, nc,       state%reshn, 'reshn')
            call ReadArray(filespec, nc*ns,    state%resmo, 'resmo')
            call ReadArray(filespec, nc,       state%resmt, 'resmt')
            call ReadArray(filespec, nc,       state%respo, 'respo')
            call ReadArray(filespec, nc,       state%reskt, 'reskt')
            call ReadArray(filespec, nc,       state%reszt, 'reszt')

        endif

        ! Housekeeping
        close(filespec)

    end subroutine

    !------------------------------------------------------------------!
    !                            Environment                           !
    !------------------------------------------------------------------!
    ! Environment
    !============
    ! General constructor
    subroutine ConstructEnvironment(environment, environmentoptions)

        ! Description
        !============
        ! Construct the environment for the grid optimization. Here, all 
        ! other peripheral structures, vessel data, ... should be stored
        ! that may be needed in goat (apart from the 
        ! magnetic field, which, due to it's main role in the grid 
        ! generation, has it's own structure). It goes without saying that 
        ! this routine is heavily case dependent. 
    
        ! Notes
        !======
        ! The vessel is added and also state information for a simulation.
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        type(EnvironmentOptionsUDT)       :: environmentoptions 
        type(EnvironmentUDT)              :: environment
    
        ! Loop variables
    
        ! Auxiliary variables 
        integer                             :: filespecifier
    
        ! Additional environment structures
        type(VesselOptionsUDT)              :: vesseloptions
    
        ! Data
        data filespecifier /60/
    
        ! Initialize
        !===========
        ! Check which type of environment to construct
        select case (trim(environmentoptions%type))
    
        case ('vessel')
    
            ! Read in vessel options
            vesseloptions%inputfilepath = environmentoptions%inputfilepath
            call vesseloptions%Set()
    
            ! Read in the vessel structure
            vesseloptions%filepath = environmentoptions%vesselfilepath
            call ReadVessel(filespecifier, environment%vessel, &
                environment%triangulationvessel, vesseloptions)
    
            ! Extract the vessel data
            call ExtractVesselData(environment%vessel, vesseloptions)
            
            ! Extract triangulation vessel data (options assumed to be 
            ! checked in ReadVessel!)
            call ExtractVesselData(environment%triangulationvessel, vesseloptions) 
    
        case default
    
            ! Unknown case, throw error
            call gdErrorHandler('Unknown environment type')
    
        end select

        ! State information
        call ReadState(environment%SOLPSstate, environmentoptions)
    
    end subroutine

    ! Vessel structures
    !==================
    ! Allocation
    subroutine AllocateVesselStructure(vesselstructure)

        ! Description
        !============
        ! Routine to allocate the vessel structure quantities assuming
        ! that np is given.

        ! The usual
        implicit none

        ! Declare variables
        type(VesselStructureUDT)          :: vesselstructure

        ! Allocate
        !=========
        allocate(vesselstructure%x(vesselstructure%np))
        allocate(vesselstructure%y(vesselstructure%np))

    end subroutine

    ! Reader (structure.dat or vessel.dat)
    ! Reading vessel data from structure.dat file
    subroutine read_structure(filespecifier, vessel, vesseloptions)

        ! Description
        !============
        ! Read the vessel data from a structure.dat file. This file 
        ! should be strictly formatted (though not checked) as follows:
        ! - The first line should contain the number of structures 
        ! - The next two lines are skipped
        ! - Starting from the fourth line, the structures follow. Each 
        ! structure should first start with a header (which is ignored)
        ! and subsequently with the number of points in the 
        ! structure, where the sign of this number indicates whether the 
        ! polygon should be closed (negative) or open (positive).

        ! Notes:
        ! The vessel should have an allocatable array of substructures
        ! of the type 'VesselStructure'. This array is allocated while
        ! reading in the number of vessel structures.

        ! Declare variables
        !==================
        ! Arguments
        type(VesselUDT)                         :: vessel
        type(VesselOptionsUDT), intent(in)      :: vesseloptions
        integer(I8), intent(in)                 :: filespecifier

        ! Loop variables
        integer(I8)                             :: i, j

        ! Auxiliary variables
        integer(I4)                             :: nstruct, npoints
        character(C32)                          :: dummy

        ! Initialize
        !===========
        ! Print from where we're reading
        print *, 'reading vessel from file: ', vesseloptions%filepath

        ! Open the file
        open(unit = filespecifier, file = vesseloptions%filepath)

        ! Read
        !=====
        ! Read the amount of structures
        read(filespecifier, *) nstruct
        print *, 'there are ', nstruct, ' structures present'
        print *, 'vessel structure ID |vessel structure label | ' // &
            'number of points | is closed'

        ! Allocate
        vessel%nstructures = nstruct
        allocate(vessel%structures(nstruct))

        ! Skip the next line
        read(filespecifier, *)

        ! Read in structures
        do i = 1, nstruct
            ! Read the header (structure <structurelabel>)
            vessel%structures(i)%ID = int(i, kind=I4) 
            read(filespecifier, *) dummy, vessel%structures(i)%label

            ! Read the number of points
            read(filespecifier, *) npoints

            ! Allocate
            vessel%structures(i)%np = abs(npoints)
            vessel%structures(i)%isclosed = (npoints .gt. 0)
            call AllocateVesselStructure(vessel%structures(i))

            ! Print
            print *, vessel%structures(i)%ID,  vessel%structures(i)%label, &
                vessel%structures(i)%np, vessel%structures(i)%isclosed

            ! Read coordinates
            do j = 1, vessel%structures(i)%np
                ! First x, then y coordinate
                read(filespecifier, *) vessel%structures(i)%x(j), vessel%structures(i)%y(j)
            end do
        end do

        ! Close the file
        !===============
        close(filespecifier)
        

    end subroutine

    ! Vessel
    !=======
    ! Allocation
    subroutine AllocateVessel(vessel)

        ! Description
        !============
        ! Allocate the vessel quantities assuming that np and ntp are 
        ! given. 

        ! The usual
        implicit none

        ! Declare variables
        type(VesselUDT)          :: vessel

        ! Allocate
        !=========
        !allocate(vessel%x(vessel%nv))
        !allocate(vessel%y(vessel%nv))
        !allocate(vessel%TPind(vessel%nv))
        !allocate(vessel%allTPind(vessel%ntv))

    end subroutine

    ! Deallocation
    subroutine DeallocateVessel(vessel)

        ! Description
        !============
        ! Deallocate the vessel quantities.

        ! The usual
        implicit none

        ! Declare variables
        type(VesselUDT)          :: vessel

        ! Allocate
        !=========
        deallocate(vessel%TPind)
        deallocate(vessel%allTPind)

    end subroutine

    ! Readers
    subroutine ReadVessel(filespecifier, vessel, triangulationvessel, vesseloptions)

        ! Description
        !============
        ! Read in the vessel data. 
    
        ! Notes
        !======
    
        ! Initialize
        !===========
        ! Declare modules
        use goatmod_userinput
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        integer                         :: filespecifier
        type(VesselOptionsUDT)          :: vesseloptions 
        type(VesselUDT)                 :: vessel, triangulationvessel
        type(DivGeoDataUDT)             :: dgdata
    
        ! Loop variables
    
        ! Auxiliary variables 
        type(VesselStructureUDT), allocatable   :: structures(:), &
            triangulationstructures(:)
        integer(I8)                             :: flag
    
        ! Main program
        !=============
        ! Check how to read the vessel structure
        select case(vesseloptions%readmeth)
            
        case ('read_structure')
    
            ! Read in the separate vessel structures
            call read_structure(filespecifier, vessel, vesseloptions)
    
            ! No additional information on triangulation vessel, so will
            ! be the same as the original vessel 
            triangulationvessel = vessel

        case ('read_dg')

            ! Check options: no structures are allowed to be deleted to
            ! avoid issues with void polygon! These structures should then 
            ! not be included in the divgeo setup...
        
            ! Check options
            if (size(vesseloptions%exclude) > 0) then 
                ! Throw hard error for now
                call gdErrorHandler('ReadVessel: when reading dgo output, ' // & 
                    'vessel parts may not be excluded using the exclude option ' // & 
                    'of goat. If structures should be excluded, do this in ' // & 
                    'the divgeo setup.')
            end if 

            ! Read from DivGeo file
            call dgdata%Initialize()
            call dgdata%Read(vesseloptions%filepath)

            ! Read plasma vessel structures
            call dgdata%ExtractVesselStructures(structures, flag)

            ! Check
            if (flag > 0) then 
                call gdErrorHandler('ReadVessel: could not read in ' // & 
                    'vessel structure data from DivGeo file')
            end if 

            ! Initialize vessel
            vessel%nstructures = int(size(structures), kind=I4)
            vessel%structures = structures

            ! Read triangulation structures
            call dgdata%ExtractTriangulationStructures(triangulationstructures, flag)
            if (flag > 0) then 
                call gdErrorHandler('ReadVessel: could not read in ' // & 
                    'triangulation vessel structure data from DivGeo file')
            end if

            ! Initialize triangulation vessel
            triangulationvessel%nstructures = int(size(triangulationstructures), kind=I4)
            triangulationvessel%structures = triangulationstructures
    
        case default
    
            ! Throw error
            call gdErrorHandler('Unknown vessel reading method')
    
        end select
    
    end subroutine

    ! Writers
    subroutine WriteVesselPolygon(vessel, filename)

        ! Description
        !============
        ! This routine writes out the vessel polygon (set) in 
        ! structure.dat format.

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)            :: vessel 
        character(*), intent(in)    :: filename

        ! Auxiliary
        type(VesselStructureUDT), allocatable, dimension(:)     :: structures

        ! Convert to structures
        !======================
        call GetStructuresFromPolygonset(vessel%polygonset, structures)
        call WriteStructureFile(filename, structures)

    end subroutine

    ! Data addition
    subroutine ExtractVesselData(vessel, vesseloptions)

        ! Description
        !============
        ! Process the vessel structure to provide the necessary data for the
        ! grid deformation. Here, we construct a polygon set representation 
        ! of the vessel to easier manipulate the vessel structure. It is 
        ! assumed that the vessel structure contains segments that 
        ! properly intersect each other such that a single closed polygon
        ! can be extracted (nested polygons are supported as long as they 
        ! are simply nested, i.e. the interior of the domain can be 
        ! distinguished). 
    
        ! Notes
        !======
        ! Note 1: (20/03/2024) Added a separate polygon representation 
        ! for each target structure. Each structure that has to be a target
        ! should be identified by the field 'TP' in the vesseloptions 
        ! structure. Each integer i corresponds to the i-th vessel structure
        ! present in vessel%structures. Given that the aim of this 
        ! additional structure is to easily determine distributions etc
        ! around these structures, we artifically close the polygon (if not
        ! already closed), since this allows more types of distribution
        ! functions to be used. Note that self-intersecting polygons are not 
        ! allowed for targets and will cause an error. 
    
        ! Note 2: (27/05/2024) Added labels in polygon set representation with 
        ! the following format: label(:, 1), label(:, 2) are vessel structure
        ! ID(s) of the vertex (max. 2 allowed per vertex), label(:, 3)  is 
        ! the unique vertex ID (may, after intersections, not go from 1 to 
        ! number of vertices due to exclusion of vertices)

        ! Note 3: (1/10/2025) Added additional polygon 'level' label in 
        ! labels(:, 4), which indicates the nestedness level of the 
        ! polygon. 
    
        ! Initialize
        !===========
        ! Declare modules
        use PolygonLevelsetFunction2D
        use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    
        ! The usual
        implicit none 
    
        ! Declare variables
        !==================
        ! Arguments
        type(VesselUDT), intent(inout)      :: vessel
        type(VesselOptionsUDT), intent(in)  :: vesseloptions
    
        ! Loop variables
        integer(I8)                         :: i, j, k, status, ks, ke, pc
    
        ! Auxiliary variables 
        type(PolygonSetUDT)                 :: tempps, polygonset_noref
        type(PolygonUDT)                    :: temppol
        integer(I8)                         :: nvp, cc, nexcl, nTP, tne, &
            vID, tl, flag
        integer(I8), allocatable            :: tv(:), templabels(:, :), &
            vesselIDmap(:), uniquelabels(:)
        real(R8), allocatable               :: tempx(:), &
            tempy(:), tx(:), ty(:), xp(:), yp(:), tvx(:), tvy(:), tnxp(:), &
            tnyp(:), tnnp(:), tnx(:), tny(:), tnn(:)
    
        class(PLF2DOptionsUDT), allocatable :: plfoptions
    
        ! Plotting
    
        ! Set NaN
        real(R8)                            :: NaN
    
        ! Initialize
        !===========
        ! Set NaN
        NaN = IEEE_VALUE(nan, IEEE_QUIET_NAN)
        
        ! Set vertex ID
        vID = 0
    
        ! Associate
        associate(&
            nvs         => vessel%nstructures,  &
            vs          => vessel%structures)
    
        ! Vessel polygon representation
        !==============================
        ! Count total number of points currently in the structures
        nvp = 0
        nexcl = 0
        do i = 1, nvs
            if (.not. any(vesseloptions%exclude == i)) then 
                nvp = vessel%structures(i)%np + nvp 
            else 
                nexcl = nexcl + 1
            end if 
        end do
    
        ! Allocate (account for NaNs)
        allocate(vesselIDmap(nvs - nexcl))
        allocate(tempx(nvp+nvs-1-nexcl), tempy(nvp+nvs-1-nexcl))
        allocate(templabels(size(tempx), 4))
        templabels = 0
        
        ! Set vertices
        cc = 0
        pc = 0
        do i = 1, nvs 
            if (.not. any(vesseloptions%exclude == i)) then 
                ! Mapping 
                pc = pc + 1
                vesselIDmap(pc) = vs(i)%label ! now propagating label instead of ID
    
                ! Coordinates
                tempx(cc+1:cc+vs(i)%np) = vs(i)%x
                tempy(cc+1:cc+vs(i)%np) = vs(i)%y
                
                ! Labels
                templabels(cc+1:cc+vs(i)%np, 1) = pc !vs(i)%ID ! Initial label is polygon-based, will be mapped back later
                templabels(cc+1:cc+vs(i)%np, 3) = [(k, k = vID+1, vID+vs(i)%np)]
    
                ! Update counters
                cc = cc + vs(i)%np
                vID = vID + vs(i)%np 
    
                ! Add NaN
                if (cc .ne. (size(tempx))) then 
                    tempx(cc+1) = NaN 
                    tempy(cc+1) = NaN 
                    templabels(cc+1, :) = 0
                    cc = cc + 1
                end if 
            end if
        end do
    
        ! Construct temporary vessel polygon set
        call tempps%Construct(tempx, tempy, templabels)
    
        ! Construct vessel polygon set 
        call ConstructVesselPolygonSet(vessel, tempps)

        ! Store the original, non-refined set
        polygonset_noref = vessel%polygonset
    
        ! Refine vessel polygonset
        !=========================
        if (vesseloptions%refine) then 
            ! Call refiner 
            call vessel%polygonset%Refine(vesseloptions%maxdist, vesseloptions%minreffac)
    
            ! Update labels - normally, labels are kept the same on original
            ! nodes but zero elsewhere. 
            do i = 1, vessel%polygonset%np
    
                ! Associate for ease
                associate(tp    => vessel%polygonset%polygons(i))
    
                ! Get polygon vertices
                tv = tp%vert
    
                ! Adjust labels by checking common non-zero labels between
                ! original vertices. Normally, first vertex should always be
                ! an original vertex
                ke = 1
                ks = 0
                do while (.true.)
                    ! Starting index
                    ks = ke 
    
                    ! Sanity check
                    if (all(tp%labels(tv(ks), :) == 0_i8)) then 
                        call gdErrorHandler('ExtractVesselData: original vertex ' // & 
                            'has no labels assigned after refinement, unexpected')
                    end if 
    
                    ! Check exit conditions
                    if (ks >= size(tv)) then 
                        exit 
                    end if 
    
                    ! Ending index
                    ke = findloc(tp%labels(tv(ks+1:tp%ne+1), 3) /= 0, .true., 1, back=.false.) + ks
    
                    ! Sanity check
                    if (ke == ks) then 
                        call gdErrorHandler('ExtractVesselData: could not ' // & 
                            'find ending vertex, unexpected')
                    end if 
    
                    if (all(tp%labels(tv(ke), :) == 0_i8)) then 
                        call gdErrorHandler('ExtractVesselData: original vertex ' // & 
                            'has no labels assigned after refinement, unexpected')
                    end if 
    
                    ! Check for common labels
                    associate(&
                        sl      => tp%labels(tv(ks), 1:2),  &
                        el      => tp%labels(tv(ke), 1:2)   &
                        )
    
                    ! Sanity checks
                    if (all(sl == 0) .or. all(el == 0)) then 
                        call gdErrorHandler('ExtractVesselData: original vertex ' // & 
                            'has no target labels after refinement, unexpected')
                    end if 
    
                    tl = 0
                    if (sl(2) == 0 .and. el(2) == 0) then
                        ! Check
                        if (sl(1) /= el(1) .and. (ke - ks > 1)) then ! fine if there are no refined vertices in between
                            ! Weird, but continue
                            print *, 'ExtractVesselData: no common labels ' // & 
                                'between subsequent vertices found, unexpected. ' // & 
                                'Taking first vertex label...' 
                        end if 
    
                        ! Assign
                        tl = sl(1)
                    elseif ((sl(2) == 0) .and. (el(2) /= 0)) then 
                        ! Check
                        if (all(sl(1) /= el) .and. (ke - ks > 1)) then ! fine if there are no refined vertices in between
                            ! Weird, but continue
                            print *, 'ExtractVesselData: no common labels ' // & 
                                'between subsequent vertices found, unexpected. ' // & 
                                'Taking first vertex label...' 
                        end if 
    
                        ! Assign
                        tl = sl(1)
                    elseif ((sl(2) /= 0) .and. (el(2) == 0)) then 
                            ! Check
                            if (all(el(1) /= sl) .and. (ke - ks > 1)) then ! fine if there are no refined vertices in between
                                ! Weird, but continue
                                print *, 'ExtractVesselData: no common labels ' // & 
                                    'between subsequent vertices found, unexpected. ' // & 
                                    'Taking first vertex label...' 
                            end if 
        
                            ! Assign
                            tl = el(1)
                    else
                        ! Both are non-zero, need to check common one
                        if (any(sl(1) == el) .and. any(sl(2) == el)) then 
                            if (sl(1) /= sl(2)) then 
                                ! Cannot distinguish, issue warning
                                print *, 'ExtractVesselData: both labels ' // & 
                                    'between subsequent vertices are present. ' // & 
                                    'Cannot distinguish, taking first vertex label...' 
                            end if 
                            tl = sl(1)
                        elseif (any(sl(1) == el) .and. .not. any(sl(2) == el)) then 
                            tl = sl(1)
                        elseif (any(sl(2) == el) .and. .not. any(sl(1) == el)) then 
                            tl = sl(2)
                        else
                            ! Could not find a common label - unexpected
                            print *, 'ExtractVesselData: no common labels ' // & 
                                'between subsequent vertices found, unexpected. ' // & 
                                'Taking first vertex label...' 
                            tl = sl(1)
                        end if 
    
                    end if 
    
                    ! Assign labels
                    tp%labels(tv(ks+1:ke-1), 1) = tl 

                    ! Housekeeping
                    end associate
                end do 

                ! Check level - should be the same everywhere (if any value
                ! is a non-zero value, )
                if (any(tp%labels(:, 4) /= 0)) then 
                    ! Get all unique labels
                    call Unique(tp%labels(:, 4), uniquelabels)

                    ! Sanity check: only one non-zero value should be found
                    if (count(uniquelabels /= 0) > 1) then
                        call gdErrorHandler('ExtractVesselData: polygon ' // & 
                            'has multiple levels assigned after refinement, ' // & 
                            'this should not happen.')
                    end if 

                    ! Get that value and assign label to all vertices
                    uniquelabels = pack(uniquelabels, uniquelabels /= 0)
                    tp%labels(:, 4) = uniquelabels(1)
                end if 
            
                ! Housekeeping
                end associate
            end do
    
            ! Test orientation
            call vessel%polygonset%OrientNestedClosedPolygons(flag)
    
            ! Check
            if (flag .ne. 0) then  
                ! Throw error
                print *, 'flag: ', flag
                call gdErrorHandler('ConstructVesselPolygon: could not ' // & 
                    'orient polygons, OrientNestedClosedPolygons exited with flag above')
            end if 
    
            ! Write data
            call vessel%polygonset%WriteData('vesselpolygon')
        end if 
    
        ! Remap the labels
        !=================
        do i = 1, vessel%polygonset%np 
            ! Initialize
            associate(tp    => vessel%polygonset%polygons(i)) 
    
            ! Remap if nonzero
            where (tp%labels(:, 1) /= 0) tp%labels(:, 1) = vesselIDmap(tp%labels(:, 1))
            where (tp%labels(:, 2) /= 0) tp%labels(:, 2) = vesselIDmap(tp%labels(:, 2))
    
            ! Add vertex IDs 
            do j = 1, tp%nv
                if (tp%labels(j, 3) == 0) then 
                    vID = vID + 1
                    tp%labels(j, 3) = vID 
                end if 
            end do 
    
            ! Housekeeping
            end associate 
        end do
        do i = 1, polygonset_noref%np 
            ! Initialize
            associate(tp    => polygonset_noref%polygons(i)) 
    
            ! Remap if nonzero
            where (tp%labels(:, 1) /= 0) tp%labels(:, 1) = vesselIDmap(tp%labels(:, 1))
            where (tp%labels(:, 2) /= 0) tp%labels(:, 2) = vesselIDmap(tp%labels(:, 2))
    
            ! Add vertex IDs 
            do j = 1, tp%nv
                if (tp%labels(j, 3) == 0) then 
                    vID = vID + 1
                    tp%labels(j, 3) = vID 
                end if 
            end do 
    
            ! Housekeeping
            end associate 
        end do

        ! Target polygon representation
        !==============================
        ! Check how many targets there are
        nTP = size(vesseloptions%TP)
    
        ! If there are targets, check that TP does not exceed the number of
        ! structures
        if (nTP > 0) then 
            if (any(vesseloptions%TP > nvs)) then 
                ! Throw error
                call gdErrorHandler('ExtractVesselData: some target plate ' // &
                    'indices are larger than number of vessel structures, ' // &
                    'check input')
            end if
        end if
    
    
        ! Allocate
        allocate(vessel%targetpolygons(nTP))
    
        ! Construct polygons
        do i = 1, nTP
    
            ! Construct test polygon
            call temppol%Construct(vs(vesseloptions%TP(i))%x, vs(vesseloptions%TP(i))%y)
    
            ! Check
            if (temppol%IsSelfIntersectingPolygon()) then 
                ! Throw error
                call gdErrorHandler('ExtractVesselData: target polygon is self intersecting, not supported')
            end if 
    
            ! Check if we need to close the polygon
            if (.not. temppol%isclosed) then 
                ! Display
                print *, 'ExtractVesselData: closing target polygon with vessel structure ID: ', vesseloptions%TP(i) 
    
                ! Get points in vertex order
                tne = temppol%ne 
                tv = temppol%vert 
                tvx = temppol%x(tv)
                tvy = temppol%y(tv)
    
                ! Compute normals in points
                tnx = -(tvy(2:tne+1) - tvy(1:tne))
                tny = (tvx(2:tne+1) - tvx(1:tne))
                tnn = sqrt(tnx**2 + tny**2)
                tnx = tnx/tnn 
                tny = tny/tnn
                tnxp = [tnx(1), 0.5*(tnx(1:tne-1)+tnx(2:tne)), tnx(tne)]
                tnyp = [tny(1), 0.5*(tny(1:tne-1)+tny(2:tne)), tny(tne)]
                tnnp = sqrt(tnxp**2 + tnyp**2)
                tnxp = tnxp/tnnp 
                tnyp = tnyp/tnnp 
    
                ! Shift the points slightly 
                tx = [tvx, tvx(size(tv)-1:1:-1)+1e-5*tnxp(size(tv)-1:1:-1), tvx(1)]
                ty = [tvy, tvy(size(tv)-1:1:-1)+1e-5*tnyp(size(tv)-1:1:-1), tvy(1)]
    
                call Write2DPolygonData(tx, ty, 'testpolyg')
    
                ! Close the polygon
                !tv2 = [tv, tv(size(tv)-1:1:-1)]
    
                ! Get new coordinates
                !tx = temppol%x(tv2)
                !ty = temppol%y(tv2)
    
                ! Construct the polygon
                call temppol%Deallocate()
                call temppol%Construct(tx, ty)
    
            end if 
    
            ! Assign
            vessel%targetpolygons(i) = temppol
    
            ! Deallocate test polygon
            call temppol%Deallocate()
    
    
        end do 
    
        ! Construct target polygonset
        call vessel%targetps%Construct(vessel%targetpolygons)
    
        ! Construct polygon representations
        !==================================
        ! Check how to construct
        select case (trim(vesseloptions%shapemeth))
    
        case ('polygon')
    
            ! Exact polygon representation
            allocate(PLF2DGeneralOptionsUDT::plfoptions)
    
            ! Set options (nothing to do here)
    
        case ('closedpolygon_exact')
    
            ! Exact representation of closed polygon
            allocate(PLF2DClosedExactOptionsUDT::plfoptions)
    
            ! Set options (nothing to do here)
    
        case ('closedpolygon_smoothapproximation')
    
            ! Approximate representation of closed polygon
            allocate(PLF2DClosedApproximationOptionsUDT::plfoptions, stat=status)
    
    
            ! Set options
            select type (plfoptions)
    
            type is (PLF2DClosedApproximationOptionsUDT) 
    
                ! Interpolation settings
                plfoptions%meth = 'uniformgrid'
                plfoptions%resx = vesseloptions%resx
                plfoptions%resy = vesseloptions%resy
                plfoptions%offsetx = vesseloptions%offsetfracx
                plfoptions%offsety = vesseloptions%offsetfracy
                plfoptions%C = vesseloptions%interpC
                plfoptions%M = vesseloptions%interpM
                plfoptions%xrange = vesseloptions%xrange
                plfoptions%yrange = vesseloptions%yrange
    
                if (size(plfoptions%xrange, 1) < 2) then 
                    ! Reset
                    call vessel%polygonset%GetVertices(xp, yp)
                    plfoptions%xrange = xp 
                    plfoptions%yrange = yp
                end if
    
            end select
    
        case default 
    
            ! Throw error
            call gdErrorHandler('ExtractVesselData: vessel shapemeth not implemented')
    
        end select
    
        ! Construct
        call InitializePolygonLevelsetFunction2D(vessel%plfvessel, vessel%polygonset, plfoptions)
        call InitializePolygonLevelsetFunction2D(vessel%plftarget, vessel%targetps, plfoptions)
        call vessel%exactplfvessel%Initialize(vessel%polygonset)
        call vessel%exactplfvessel_noref%Initialize(polygonset_noref)

        ! Write
        call vessel%exactplfvessel%VisualizeLabel('levelset_label_1', labelindin=1)
        call vessel%exactplfvessel%VisualizeLabel('levelset_label_3', labelindin=3)
        call vessel%exactplfvessel%VisualizeLabel('levelset_label_4', labelindin=4)
    
        end associate 
    
    end subroutine

    ! Vessel polygon construction
    subroutine ConstructVesselPolygonSet(vessel, ps, compress)

        ! Description
        !============
        ! ConstructVesselPolygonSet constructs a closed vessel polygon set (possibly
        ! consisting of multiple closed polygons) starting from an initial set of
        ! polygons, given in the oldvessel input argument (see later for input
        ! specifications). The routine has support for multiple open and closed polygon
        ! boundaries, though the closed polygons (or the closed polygons resulting
        ! from combining multiple open polygons) should not intersect themselves
        ! with other open or closed polygons. Otherwise, an error will be thrown.
        ! See the notes for more explanation on how closed polygons are identified
        ! and how open polygons are merged to create a closed one. 
    
        ! Algorithm
        !==========
        ! 1) Check which structures are closed. For these closed structures, check
        ! if they have any intersections with other structures -> throw error if
        ! this is the case
        ! 2) Compute intersections of the remaining open polygons with all other
        ! polygons. If not exactly two intersections are found, an error is thrown,
        ! as this is indicates that the input assumptions are violated.
        ! 3) Construct the different closed polygon parts by looking at the
        ! intersections of each polygon.
        ! 4) Construct the new vessel polygon(s) by concatenating the polygon
        ! segments that were found. This is also the stage where target plate
        ! indices etc are constructed. 
    
        ! Notes
        !------
        ! Note 1: closed polygons are identified by comparing their end points. 
        ! These should be the same up to 100*macheps (to get rid of possible 
        ! numerical artefacts). This is measured in terms of Euclidean distance
        ! between the points (i.e. dist <= 100*macheps
    
        ! Note 2: open polygons can either truly intersect with each other at
        ! certain edges or coincide in one of their end vertices. The latter is
        ! checked by taking the distance between the end points and comparing
        ! to 100*macheps (same criterion as above for determining closed polygons).
        ! For polygon intersections, the end sections are cut off. Intersection
        ! computation is done using the 'intersections' routine
    
        ! Note 3: the resulting vessel polygon(s) is(are) arbitrarily sorted
        ! (either clockwise or counterclockwise). 
    
        ! Note 4: we propagate the structure number(s) of each vertex and 
        ! the vertex IDs as they were present in the original vessel 
        ! structure. Note that vertices are only allowed up to two structure
        ! IDs and only one vertex ID. Vertex IDs will be removed since 
        ! either the structures intersect exactly in the end points and only
        ! one node is retained, or the structures intersect somewhere along
        ! the polygon and the 'dangling' nodes are removed. It is assumed 
        ! that this information is already available through the polygon
        ! labels field, which should be setup in ExtractVesselData.F90. 

        ! Note 5: sometimes it may be desireable to 'compress' the 
        ! initial polygon set, i.e. remove any redundant straight edges
        ! that do not influence the actual form of the polygon. In that
        ! case, the optional argument 'compress' should be parsed and 
        ! set to true. Machine precision effects may arise though. 
    
        ! Modules
        !========
        implicit none
    
        ! Declare variables
        !==================
        ! Arguments
        type(VesselUDT)             :: vessel
        type(PolygonSetUDT)         :: ps 
        logical, intent(in), optional   :: compress
    
        ! Auxiliary
        integer(I8)                 :: ni, nfinpol, thisp, firstpolygon, &
            c1, c2, nvest, nvv, tempnv, nextp, sv, ev, flag, vID, indis, indie
        logical                     :: polygonnotfound, doflip, docompression
        real(R8)                    :: nan, xs, ys, xe, ye 
    
        real(R8), allocatable       :: xi(:), yi(:), tempx(:), tempy(:), &
            xv(:), yv(:)
        integer(I8)                 :: pis(1:2), pie(1:2)
        integer(I8), allocatable    :: p1(:), p2(:), s1(:), s2(:), &
            polcat(:, :), npol(:), remp(:), tempp(:), indi(:, :), &
            si(:, :), pi(:, :), ci(:), templabels(:, :), &
            labelsv(:, :), polygonlevels(:)
        logical, allocatable        :: notfound(:)
    
        character(:), allocatable   :: vesselpath
    
        ! Loop
        integer(I8)                 :: i, j, k
    
        ! Initialize
        !===========
        ! Checks 
        if (present(compress)) then
            docompression = compress
        else
            docompression = .false.
        end if 
        
        ! If desired, compress the given polygon set 
        if (docompression) then 
            call ps%Compress()
        end if 
    
        ! Set NaN
        nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    
        ! Get maximal vertex ID label - if available
        vID = 0
        if (size(ps%polygons(1)%labels, 2) >= 3) then 
            do i = 1, ps%np 
                vID = max(maxval(ps%polygons(i)%labels(:, 3)), vID)
            end do 
        else
            ! Throw warning
            print *, 'ConstructVesselPolygonSet: polygons do not have ' // & 
                'all labels, constructing vessel structure and vertex labels ' // & 
                'based on polygons'
            do i = 1, ps%np
                deallocate(ps%polygons(i)%labels)
                allocate(ps%polygons(i)%labels(ps%polygons(i)%nv, 4))
                ps%polygons(i)%labels(:, 1) = i 
                ps%polygons(i)%labels(:, 2) = 0
                ps%polygons(i)%labels(:, 3) = [(k, k = vID+1, vID+ps%polygons(i)%nv)]
                vID = vID + ps%polygons(i)%nv 
            end do 
        end if 
    
        ! Determine intersections
        !========================
        ! Compute all intersections
        call ps%SelfIntersections(xi, yi, p1, p2, s1, s2)
        
        ! Check if each polygon has either exactly zero or two intersections
        do i = 1, ps%np 
            ! Count
            ni = count(i == p1) + count(i == p2)
    
            ! Check
            if ( ni == 1 ) then 
                ! This shouldn't be happening with the right input format
                ! -> probably something wrong on the input side 
                print *, 'polygon number (excluded polygons not accounted):'
                print *, i
                call gdErrorHandler('ConstructVesselPolygon: only one ' &
                    // 'intersection found for this polygon, check input')
            elseif ( ni > 2 ) then 
                ! This shouldn't be happening with the right input format
                ! -> probably something wrong on the input side 
                print *, 'polygon number (excluded polygons not accounted):'
                print *, i
                call gdErrorHandler('ConstructVesselPolygon: too many ' &
                    // 'intersections found for this polygon, check input')
            elseif ( ni == 0 ) then 
                if (.not. ps%polygons(i)%isclosed) then 
                    ! Open polygons should intersect twice, this shouldn't 
                    ! be happening
                    print *, 'polygon number (excluded polygons not accounted):'
                    print *, i
                    call gdErrorHandler('ConstructVesselPolygon: not enough ' &
                        // 'intersections found for this polygon, check input')
                end if 
            elseif ( (ni .ne. 0) .and. ps%polygons(i)%isclosed) then 
                ! Closed polygons should not intersect with itself (this is
                ! topologically possible but is not allowed here)
                print *, 'polygon number (excluded polygons not accounted):'
                print *, i
                print *, ps%polygons(i)%edges
                call gdErrorHandler('ConstructVesselPolygon: closed ' &
                    // 'polygon found that intersects with itself, check input')
    
            end if
    
        end do
    
        ! Compute closed polygons
        !========================
        ! First, we check which polygons form one single closed polygon and
        ! how many closed polygons there are. 
    
        ! Allocate
        allocate(notfound(ps%np), polcat(ps%np, ps%np), &
            npol(ps%np))
    
        ! Initialize
        notfound(:) = .true. 
        polcat(:, :) = 0 ! concatenation structure: each column holds the polygon indices of a closed polygon
        npol(:) = 0 ! counter for how many polygons belong to each final polygon
        nfinpol = 0 ! number of final polygons
    
        ! Loop
        nfinpol = 0
        do while (any(notfound))
            ! Get the next polygon piece
            allocate(remp(count(notfound))) 
            remp = pack([(k, k=1, ps%np)], notfound)
            thisp = remp(1) 
    
            ! If polygon is closed -> add
            if (ps%polygons(thisp)%isclosed) then 
    
                ! Update total number of closed polygons
                nfinpol = nfinpol + 1
               
                ! Set logical to false 
                notfound(thisp) = .false.
    
                ! Update polygon structure counter
                npol(nfinpol) = npol(nfinpol) + 1 
    
                ! Add to polygon structure
                polcat(npol(nfinpol), nfinpol) = thisp 
    
            else 
                ! Polygon is not closed, need to loop untill we find all
                ! polygons for this closed polygon
    
                ! Update total number of closed polygons
                nfinpol = nfinpol + 1
    
                ! Set this polygon as the first
                firstpolygon = thisp 
    
                ! Loop
                polygonnotfound = .true. 
                do while (polygonnotfound) 
                    
                    ! Save the current polygon
                    npol(nfinpol) = npol(nfinpol) + 1 
                    polcat(npol(nfinpol), nfinpol) = thisp
                    notfound(thisp) = .false. 
    
                    ! Get next polygon candidates
                    c2 = count( (thisp == p1) .and. (notfound(p2)) ) ! don't include polygons that were found already
                    c1 = count( (thisp == p2) .and. (notfound(p1)) ) 
    
                    ! Check the amount of candidates
                    if ( (c2 + c1) == 1) then 
                        ! Just one candidate - ok 
                        if (c1 == 0) then 
                            ! Update thisp 
                            allocate(tempp(c2))
                            tempp = pack(p2, (thisp == p1) .and. (notfound(p2)))
                            thisp = tempp(1)
                            deallocate(tempp)
                        else 
                            ! Update thisp 
                            allocate(tempp(c1))
                            tempp = pack(p1, (thisp == p2) .and. (notfound(p1)))
                            thisp = tempp(1)
                            deallocate(tempp)
                        end if  
                    elseif ( (c2 + c1) == 2) then 
                        ! Two candidates - only possible if this is the 
                        ! first polygon.
                        if (thisp .ne. firstpolygon) then  
                            ! Throw error
                            print *, 'polygon number: ', thisp
                            call gdErrorHandler('ConstructVesselPolygon: ' &
                                // 'multiple polygons found for polygon ' &
                                // 'that is not the starting polygon. ' &
                                // 'Check vessel input') 
                        end if 
    
                        ! Simply take one of the two
                        if (c1 > 0) then 
                            ! Update thisp 
                            allocate(tempp(c1))
                            tempp = pack(p1, (thisp == p2) .and. (notfound(p1)))
                            thisp = tempp(1)
                            deallocate(tempp)
                        else 
                            ! Update thisp 
                            allocate(tempp(c2))
                            tempp = pack(p2, (thisp == p1) .and. (notfound(p2)))
                            thisp = tempp(1)
                            deallocate(tempp)
                        end if 
                    elseif ( (c2 + c1) == 0) then 
                        ! We should have reached the end of the polygon here
                        ! - need to check if the next polygon here would be 
                        ! the first one again
                        c2 = count( (thisp == p1) .and. (p2 == firstpolygon) )
                        c1 = count( (thisp == p2) .and. (p1 == firstpolygon) ) 
                        if ( (c1 + c2) >= 1) then ! might be 2 as well in the case of exactly two polygons
                            ! Ok, exit the loop
                            polygonnotfound = .false. 
                        else 
                            ! Not OK, throw error - this means an open 
                            ! vessel structure was found
                            call gdErrorHandler('ConstructVesselPolygon: ' &
                                // 'Open vessel structure detected, ' &
                                // 'check vessel input') 
                        end if 
                    elseif ( (c2 + c1) > 2) then 
                        ! This should actually not happen anymore since this
                        ! is checked upfront, but still we handle this error
                        ! should it occur due to unknown causes. 
                        call gdErrorHandler('ConstructVesselPolygon: ' &
                            // 'more than two polygons found that ' &
                            // 'intersect with this polygon, ' &
                            // 'check vessel input')
                    end if
                end do
            end if 
    
            ! Housekeeping
            deallocate(remp)
        end do 
    
        ! Construct new vessel
        !=====================
        ! Here, we 'merge' all the polygon pieces together for each closed
        ! polygon and reconstruct the vessel polygon set. To this end, we
        ! must make sure that the polygon pieces are properly trimmed at the
        ! intersection points and that the order of the points is correct as
        ! well. Since it is tedious (but not impossible) to precompute the 
        ! final total number of polygon vertices in order to allocate the 
        ! coordinate arrays, we make an overestimation and allocate too 
        ! much memory to trim later. 
        !
        ! The segments are sorted as follows: for each segment, we check
        ! which intersection corresponds to the intersection with the next
        ! polygon in the polcat structure. If that intersection has the 
        ! highest segment index of both intersections in that polygon, the 
        ! polygon is correctly oriented. Otherwise, we need to reverse the
        ! vertex order of this polygon before adding it to the vessel 
        ! structure. 
    
        ! First, cast the intersection data in easier to navigate 
        ! structures. We can safely assume here that there are either two 
        ! or zero intersections per polygon. We store:
        ! - pi: polygon indices of the intersections
        ! - si: section indices of the intersection of the polygon
        ! - indi: the index of the intersection (to retrieve the x- and y- 
        ! coordinates of this intersection as xi(indi) ) 
        ! - ci is a counter
        allocate(indi(2, ps%np), si(2, ps%np), pi(2, ps%np), ci(ps%np))
        indi(:, :)  = 0
        si(:, :)    = 0
        pi(:, :)    = 0
        ci(:)       = 0
    
        do i = 1, size(xi)
            ! First polygon contribution
            thisp = p1(i)
            ci(thisp) = ci(thisp) + 1
            indi(ci(thisp), thisp) = i
            si(ci(thisp), thisp) = s1(i)
            pi(ci(thisp), thisp) = p2(i)
            
            ! Second polygon contribution
            thisp = p2(i)
            ci(thisp) = ci(thisp) + 1
            indi(ci(thisp), thisp) = i
            si(ci(thisp), thisp) = s2(i)
            pi(ci(thisp), thisp) = p1(i)
        end do
    
        ! Make estimate of number of vessel vertices
        nvest = 0
        do i = 1, ps%np 
            nvest = nvest + ps%polygons(i)%nv 
        end do 
        nvest = nvest*2 ! factor 2 just to be sure 
    
        ! Allocate
        allocate(tempx(nvest), tempy(nvest), templabels(nvest, 4))
    
        ! Loop 
        nvv = 0 ! vessel vertex counter
        do i = 1, nfinpol 
            ! Check number of polygon pieces
            if (npol(i) == 1) then ! Closed polygon piece
                ! Get the polygon number
                thisp = polcat(1, i)
    
                ! Get number of vertices of this polygon
                ! Note: we need to work through the %vert structure, which 
                ! has size ne+1! Otherwise, results may be wrong 
                tempnv = ps%polygons(thisp)%ne+1
    
                ! Add coordinates
                tempx(nvv+1:nvv+tempnv) = ps%polygons(thisp)%x(ps%polygons(thisp)%vert) 
                tempy(nvv+1:nvv+tempnv) = ps%polygons(thisp)%y(ps%polygons(thisp)%vert)  
    
                ! Add labels
                templabels(nvv+1:nvv+tempnv, :) = ps%polygons(thisp)%labels(ps%polygons(thisp)%vert, :) 
    
                ! Update counter
                nvv = nvv + tempnv 
    
            elseif (npol(i) == 2) then  ! special case, needs special treatment
    
                ! Print warning, this code part is not yet thoroughly tested
                print *, 'ConstructVesselPolygonSet: code part for two segments is not ' // & 
                    'yet properly tested'
    
                ! Get polygons
                thisp = polcat(1, i)
                nextp = polcat(2, i)
    
                ! First polygon
                !--------------
                ! Get polygons
                thisp = polcat(1, i)
                nextp = polcat(2, i)
                doflip = .false. 
    
                ! First intersection, check si
                indis = indi(2, thisp)
                indie = indi(1, thisp)
                if (si(1, thisp) >= si(2, thisp) ) then 
                    ! Don't flip
                    xs = xi(indis) ! start 
                    ys = yi(indis) 
                    pis = [p1(indis), p2(indis)]
                    xe = xi(indie) ! end 
                    ye = yi(indie) 
                    pie = [p1(indie), p2(indie)]
                    sv = si(2, thisp) + 1
                    ev = si(1, thisp) 
                else 
                    ! Flip
                    doflip = .true. 
                    xs = xi(indis) ! start 
                    ys = yi(indis) 
                    pis = [p1(indis), p2(indis)]
                    xe = xi(indie) ! end 
                    ye = yi(indie) 
                    pie = [p1(indie), p2(indie)]
                    sv = si(1, thisp) + 1
                    ev = si(2, thisp) 
                end if 
    
                ! Add
                    ! Starting point
                tempx(nvv+1) = xs 
                tempy(nvv+1) = ys 
                templabels(nvv+1, 1:2) = pis
                templabels(nvv+1, 3) = vID+1
    
                ! Update counters
                nvv     = nvv + 1
                vID     = vID + 1
    
                ! Polygon points
                tempnv = ev - sv + 1 
                
                if (doflip) then 
                    tempx(nvv+1:nvv+tempnv) = &
                        ps%polygons(thisp)%x(ps%polygons(thisp)%vert(ev:sv:-1))
                    tempy(nvv+1:nvv+tempnv) = &
                        ps%polygons(thisp)%y(ps%polygons(thisp)%vert(ev:sv:-1))
                    templabels(nvv+1:nvv+tempnv, :) = &
                        ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(ev:sv:-1), :)
                else 
                    tempx(nvv+1:nvv+tempnv) = ps%polygons(thisp)%x(ps%polygons(thisp)%vert(sv:ev))
                    tempy(nvv+1:nvv+tempnv) = ps%polygons(thisp)%y(ps%polygons(thisp)%vert(sv:ev))
                    templabels(nvv+1:nvv+tempnv, :) = ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(sv:ev), :)
                end if
    
                ! Hedge for zero/negative length
                if (tempnv < 0) then 
                    tempnv = 0
                end if
                nvv = nvv+tempnv
    
                ! End point
                tempx(nvv+1) = xe
                tempy(nvv+1) = ye 
                templabels(nvv+1, 1:2) = pie
                templabels(nvv+1, 3) = vID+1
    
                ! Update counters
                nvv = nvv + 1 
                vID = vID + 1
    
                ! Second polygon
                !---------------
                ! Get polygons
                thisp = polcat(2, i)
                nextp = polcat(1, i)
                doflip = .false. 
    
                ! Check intersection index
                if (indie == indi(1, thisp)) then 
                    indis = indi(1, thisp)
                    indie = indi(2, thisp)
    
                    if (si(1, thisp) <= si(2, thisp) ) then 
                        ! Don't flip
                        xs = xi(indis) ! start 
                        ys = yi(indis) 
                        pis = [p1(indis), p2(indis)]
                        xe = xi(indie) ! end 
                        ye = yi(indie) 
                        pie = [p1(indie), p2(indie)]
                        sv = si(1, thisp) + 1
                        ev = si(2, thisp) 
                    else 
                        ! Flip
                        doflip = .true. 
                        xs = xi(indis) ! start 
                        ys = yi(indis) 
                        pis = [p1(indis), p2(indis)]
                        xe = xi(indie) ! end 
                        ye = yi(indie) 
                        pie = [p1(indie), p2(indie)]
                        sv = si(2, thisp) + 1
                        ev = si(1, thisp) 
                    end if 
    
                elseif (indis == indi(1, thisp)) then 
                    indis = indi(2, thisp)
                    indie = indi(1, thisp)
    
                    if (si(1, thisp) >= si(2, thisp) ) then 
                        ! Don't flip
                        xs = xi(indis) ! start 
                        ys = yi(indis) 
                        pis = [p1(indis), p2(indis)]
                        xe = xi(indie) ! end 
                        ye = yi(indie) 
                        pie = [p1(indie), p2(indie)]
                        sv = si(2, thisp) + 1
                        ev = si(1, thisp) 
                    else 
                        ! Flip
                        doflip = .true. 
                        xs = xi(indis) ! start 
                        ys = yi(indis) 
                        pis = [p1(indis), p2(indis)]
                        xe = xi(indie) ! end 
                        ye = yi(indie) 
                        pie = [p1(indie), p2(indie)]
                        sv = si(1, thisp) + 1
                        ev = si(2, thisp) 
                    end if 
                else
                    ! This shouldn't happen
                    call gdErrorHandler('Unknown error')
                end if 
    
                ! Add
                ! Starting point
                tempx(nvv+1) = xs 
                tempy(nvv+1) = ys 
                templabels(nvv+1, 1:2) = pis
                templabels(nvv+1, 3) = vID+1
    
                ! Update counters
                nvv     = nvv + 1
                vID     = vID + 1
    
                ! Polygon points
                tempnv = ev - sv + 1 
                
                if (doflip) then 
                    tempx(nvv+1:nvv+tempnv) = &
                        ps%polygons(thisp)%x(ps%polygons(thisp)%vert(ev:sv:-1))
                    tempy(nvv+1:nvv+tempnv) = &
                        ps%polygons(thisp)%y(ps%polygons(thisp)%vert(ev:sv:-1))
                    templabels(nvv+1:nvv+tempnv, :) = &
                        ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(ev:sv:-1), :)
                else 
                    tempx(nvv+1:nvv+tempnv) = ps%polygons(thisp)%x(ps%polygons(thisp)%vert(sv:ev))
                    tempy(nvv+1:nvv+tempnv) = ps%polygons(thisp)%y(ps%polygons(thisp)%vert(sv:ev))
                    templabels(nvv+1:nvv+tempnv, :) = ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(sv:ev), :)
                end if
    
                ! Hedge for zero/negative length
                if (tempnv < 0) then 
                    tempnv = 0
                end if
                nvv = nvv+tempnv
    
                ! End point
                tempx(nvv+1) = xe
                tempy(nvv+1) = ye 
                templabels(nvv+1, 1:2) = pie
                templabels(nvv+1, 3) = vID+1
    
                ! Update counters
                nvv = nvv + 1 
                vID = vID + 1
    
    
    
    
            else ! Loop over all segments
                do j = 1, npol(i)
                    ! Get the current polygon
                    thisp = polcat(j, i) 
    
                    ! Check if the intersection with the largest si 
                    ! corresponds to the intersection with the next polygon.
                    ! If so, add the polygon as is. Otherwise, flip the 
                    ! coordinates.
    
                    ! Get the next polygon
                    if (j < npol(i)) then 
                        nextp = polcat(j+1, i)
                    else 
                        nextp = polcat(1, i)
                    end if 
    
                    ! Get start and end indices and check if we need to 
                    ! flip. Set start and end coordinates 
                    doflip = .false. 
                    if (pi(1, thisp) == nextp) then 
                        ! First intersection, check si
                        if (si(1, thisp) >= si(2, thisp) ) then 
                            ! Don't flip
                            xs = xi(indi(2, thisp)) ! start 
                            ys = yi(indi(2, thisp)) 
                            pis = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                            xe = xi(indi(1, thisp)) ! end 
                            ye = yi(indi(1, thisp)) 
                            pie = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                            sv = si(2, thisp) + 1
                            ev = si(1, thisp) 
                        else 
                            ! Flip
                            doflip = .true. 
                            xs = xi(indi(2, thisp)) ! start 
                            ys = yi(indi(2, thisp)) 
                            pis = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                            xe = xi(indi(1, thisp)) ! end 
                            ye = yi(indi(1, thisp)) 
                            pie = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                            sv = si(1, thisp) + 1
                            ev = si(2, thisp) 
                        end if 
                    elseif ( pi(2, thisp) == nextp) then 
                        ! Second intersection, check si
                        if (si(2, thisp) >= si(1, thisp) ) then 
                            ! Don't flip
                            xs = xi(indi(1, thisp)) ! start 
                            ys = yi(indi(1, thisp)) 
                            pis = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                            xe = xi(indi(2, thisp)) ! end 
                            ye = yi(indi(2, thisp)) 
                            pie = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                            sv = si(1, thisp) + 1 
                            ev = si(2, thisp) 
                        else 
                            ! Flip
                            doflip = .true. 
                            xs = xi(indi(1, thisp)) ! start 
                            ys = yi(indi(1, thisp)) 
                            pis = [p1(indi(1, thisp)), p2(indi(1, thisp))]
                            xe = xi(indi(2, thisp)) ! end 
                            ye = yi(indi(2, thisp))
                            pie = [p1(indi(2, thisp)), p2(indi(2, thisp))]
                            sv = si(2, thisp) + 1
                            ev = si(1, thisp) 
                        end if 
                    else 
                        ! Unexpected error
                        call gdErrorHandler('ConstructVesselPolygon: ' &
                            // 'encountered unexpected error when ' &
                            // 'constructing full vessel polygon')
                    end if 
    
                    ! Add
                    ! Starting point
                    tempx(nvv+1) = xs 
                    tempy(nvv+1) = ys 
                    templabels(nvv+1, 1:2) = pis
                    templabels(nvv+1, 3) = vID+1
    
                    ! Update counters
                    nvv     = nvv + 1
                    vID     = vID + 1
    
                    ! Polygon points
                    tempnv = ev - sv + 1 
                    
                    if (doflip) then 
                        tempx(nvv+1:nvv+tempnv) = &
                            ps%polygons(thisp)%x(ps%polygons(thisp)%vert(ev:sv:-1))
                        tempy(nvv+1:nvv+tempnv) = &
                            ps%polygons(thisp)%y(ps%polygons(thisp)%vert(ev:sv:-1))
                        templabels(nvv+1:nvv+tempnv, :) = &
                            ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(ev:sv:-1), :)
                    else 
                        tempx(nvv+1:nvv+tempnv) = ps%polygons(thisp)%x(ps%polygons(thisp)%vert(sv:ev))
                        tempy(nvv+1:nvv+tempnv) = ps%polygons(thisp)%y(ps%polygons(thisp)%vert(sv:ev))
                        templabels(nvv+1:nvv+tempnv, :) = ps%polygons(thisp)%labels(ps%polygons(thisp)%vert(sv:ev), :)
                    end if
    
                    ! Hedge for zero/negative length
                    if (tempnv < 0) then 
                        tempnv = 0
                    end if
                    nvv = nvv+tempnv
    
                    ! End point
                    tempx(nvv+1) = xe
                    tempy(nvv+1) = ye 
                    templabels(nvv+1, 1:2) = pie
                    templabels(nvv+1, 3) = vID+1
    
                    ! Update counters
                    nvv = nvv + 1 
                    vID = vID + 1
    
                end do 
            end if
    
            ! Add nan (unless if last polygon)
            if (i .ne. nfinpol) then 
                tempx(nvv+1) = nan 
                tempy(nvv+1) = nan
                templabels(nvv+1, :) = 0 
                nvv = nvv+1
            end if
        end do
    
        ! Construct vessel polygon set
        xv = tempx(1:nvv)
        yv = tempy(1:nvv)
        labelsv = templabels(1:nvv, :)
        call vessel%polygonset%Construct(xv, yv, labelsv)
    
        ! Test orientation
        call vessel%polygonset%OrientNestedClosedPolygons(flag, polygonlevels)
    
        ! Check
        if (flag .ne. 0) then  
            ! Throw error
            print *, 'flag: ', flag
            call gdErrorHandler('ConstructVesselPolygon: could not ' // &
                'orient polygons, OrientNestedClosedPolygons exited with ' // &
                'flag above')
        end if 

        ! Add polygonlevels to polygons
        do i = 1, vessel%polygonset%np
            vessel%polygonset%polygons(i)%labels(:, 4) = polygonlevels(i)
        end do 
    
        ! Write data
        vesselpath = 'vesselpolygon'
        call vessel%polygonset%WriteData(vesselpath)
    
        ! Housekeeping
        !=============
        deallocate(indi, si, pi, ci, tempx, tempy) 
    
    end subroutine

    ! Getting coordinates (assumed already initialized)
    subroutine GetVesselCoordinates(vessel, xv, yv)

        ! Description
        !============
        ! This routine returns the vessel coordinates in the same order
        ! as they are updated. Useful for e.g. shape optimization 
        ! purposes. 

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)                    :: vessel 
        real(R8), allocatable, intent(out)  :: xv(:), yv(:)

        ! Auxiliary
        integer(I8)                     :: nv

        ! Loop 
        integer(I8)                     :: i, cc

        ! Compute total number of vertices
        nv = 0
        do i = 1, vessel%polygonset%np
            nv = nv + vessel%polygonset%polygons(i)%nv
        end do

        ! Extract
        cc = 0
        if (allocated(xv)) then 
            deallocate(xv)
        end if 
        if (allocated(yv)) then 
            deallocate(yv)
        end if 
        allocate(xv(nv), yv(nv))
        do i = 1, vessel%polygonset%np 
            xv(cc+1:cc+vessel%polygonset%polygons(i)%nv) = vessel%polygonset%polygons(i)%x
            yv(cc+1:cc+vessel%polygonset%polygons(i)%nv) = vessel%polygonset%polygons(i)%y

            ! Update counter
            cc = cc + vessel%polygonset%polygons(i)%nv
        end do


    end subroutine

    ! Get vessel labels on query points
    subroutine GetVesselStructureLabelsOnPoints(vessel, xq, yq, labels)

        ! Description
        !============
        ! Dedicated routine to get the vessel structure labels on 
        ! arbitrary query points. This routine is more involved than
        ! just calling the EvaluateLabel routine of the exact polygon
        ! representation, because that one may be inaccurate near vessel
        ! structure transitions. Here, we evaluate first the edges/vertices
        ! on which the query points lie, and then check the labels of 
        ! the vertices (of those edges). The two first labels should 
        ! indicate the structure labels (the second one is zero if it 
        ! lies only on that structure). 

        ! Note: two label values are possible, since query points may
        ! lie on an edge and that edge may hold multiple labels, for
        ! which we cannot distinguish in this routine (this should be 
        ! dealt with upstream...)

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)                        :: vessel 
        real(R8), dimension(:), intent(in)      :: xq, yq 
        integer(I8), allocatable, dimension(:, :)   :: labels 

        ! Auxiliary
        integer(I8)                             :: nq, v1, v2
        integer(I8), allocatable, dimension(:)  :: edgeID, vertID
        integer(I8), allocatable, dimension(:, :)   :: templabels

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Unpack number of query points
        nq = size(xq)

        ! Sanity check
        if (size(xq) /= size(yq)) then 
            call gdErrorHandler('GetVesselStructureLabelsOnPoints: size mismatch ' // & 
                'between xq, yq')
        end if 

        ! Check labels status
        if (allocated(labels)) then 
            if (size(labels, 1) /= size(xq)) then 
                deallocate(labels)
                allocate(labels(nq, 2))
            end if 
        else
            allocate(labels(nq, 2))
        end if 
        
        ! Initialize
        labels = 0

        ! Compute labels
        !===============
        ! Labels are assumed to be given in first and second column
        call vessel%exactplfvessel_noref%EvaluateLabel(xq, yq, templabels, &
            edgeIDopt=edgeID, vertIDopt=vertID)

        ! Sanity check
        if (size(templabels, 2) < 2) then 
            call gdErrorHandler('GetVesselStructureLabelsOnPoints: expected at ' // &
                'least two labels per vertex but found less. Unexpected')
        end if 

        ! Unpack
        associate(&
            vlabels     => vessel%exactplfvessel_noref%vertlabel,   &
            vp1         => vessel%exactplfvessel_noref%vp1,         &
            vp2         => vessel%exactplfvessel_noref%vp2          &
            )

        ! Loop over all query points
        do i = 1, nq
            ! First, check if it lies on an edge (typically most likely 
            ! scenario)
            if (edgeID(i) /= 0) then
                ! Take the vertices of this edge
                v1 = vp1(edgeID(i))
                v2 = vp2(edgeID(i))
                
                ! Sanity checks
                if (vlabels(v1, 1) == 0 .or. vlabels(v2, 1) == 0) then 
                    call gdErrorHandler('GetVesselStructureLabelsOnPoints: ' // &
                        'vertices detected that do not have structure ' // & 
                        'labels, unexpected')
                end if 

                ! Check label occurrence
                if (vlabels(v1, 2) == 0) then 
                    ! Normally, label should simply be vlabels(v1, 1), but 
                    ! we do some sanity checks
                    if (vlabels(v2, 2) == 0) then 
                        if (vlabels(v1, 1) /= vlabels(v2, 1)) then 
                            ! Different labels, unexpected - throw warning
                            print *, 'GetVesselStructureLabelsOnPoints: each vertex ' // & 
                                'belongs to a single yet different structure, ' // & 
                                'unexpected. Taking first structure value...'
                        end if
                    else
                        if (.not. any(vlabels(v1, 1) == vlabels(v2, 1:2))) then 
                            ! Different labels, unexpected - throw warning
                            print *, 'GetVesselStructureLabelsOnPoints: each vertex ' // & 
                                'belongs to a single yet different structure, ' // & 
                                'unexpected. Taking first structure value...'
                        end if 
                    end if 

                    ! Add
                    labels(i, 1) = vlabels(v1, 1)
                else 
                    ! Need to check which one to take 
                    if (vlabels(v2, 2) == 0) then 
                        ! Normally, label should simply be vlabels(v2, 1), but 
                        ! we do some sanity checks
                        if (.not. any(vlabels(v2, 1) == vlabels(v1, 1:2))) then 
                            ! Different labels, unexpected - throw warning
                            print *, 'GetVesselStructureLabelsOnPoints: each vertex ' // & 
                                'belongs to a single yet different structure, ' // & 
                                'unexpected. Taking first structure value...'
                        end if 

                        ! Add
                        labels(i, 1) = vlabels(v2, 1)
                    else
                        ! Both are non-negative: need to check all combinations
                        if (any(vlabels(v1, 1) == vlabels(v2, 1:2)) .and. &
                            any(vlabels(v1, 2) == vlabels(v2, 1:2))) then 
                            ! Cannot distinguish between both possibilities, unexpected
                            print *, 'GetVesselStructureLabelsOnPoints: each vertex ' // & 
                                'belongs to two structures, cannot discriminate ' // & 
                                'which one to take - unexpected. Taking both structure values...'   
                            labels(i, 1:2) = vlabels(v1, 1:2)
                        elseif (any(vlabels(v1, 1) == vlabels(v2, 1:2)) .and. &
                            .not. any(vlabels(v1, 2) == vlabels(v2, 1:2))) then 
                            labels(i, 1) = vlabels(v1, 1)
                        elseif (any(vlabels(v1, 2) == vlabels(v2, 1:2)) .and. &
                            .not. any(vlabels(v1, 1) == vlabels(v2, 1:2))) then 
                            labels(i, 1) = vlabels(v1, 2)
                        else
                            ! No agreement found - unexpected    
                            print *, 'GetVesselStructureLabelsOnPoints: each vertex ' // & 
                                'belongs to multiple yet different structures, ' // & 
                                'unexpected. Taking first structure value...'
                            labels(i, 1) = vlabels(v1, 1)
                        end if 
                    end if 
                end if 
            end if  

            ! Otherwise, check the vertex
            if (vertID(i) /= 0) then 
                ! Simply return the vertex labels
                labels(i, 1:2) = vlabels(vertID(i), 1:2)
            end if
        end do 

        ! Housekeeping
        end associate 

    end subroutine

    ! Get vessel labels on edges
    subroutine GetVesselStructureLabelsOnEdges(vessel, xq1, yq1, &
        xq2, yq2, labels)

        ! Description
        !============
        ! Same routine as the '...OnPoints' one, but now the query 
        ! coordinates are assumed to form edges (e.g. face vertex 
        ! coordinates of grid faces). Here, the labels are first defined
        ! for each edge vertex. Then, these labels are compared such 
        ! that a single label is returned per edge (if possible). 

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)                        :: vessel 
        real(R8), dimension(:), intent(in)      :: xq1, xq2, yq1, yq2 
        integer(I8), allocatable, dimension(:)  :: labels 

        ! Auxiliary
        integer(I8)                             :: nq
        integer(I8), allocatable, dimension(:, :)   :: labelsv1, &
            labelsv2 

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Number of query points
        nq = size(xq1, 1)

        ! Sanity checks
        if ((nq /= size(yq1)) .or. (nq /= size(yq1)) .or. (nq /= size(yq1))) then 
            call gdErrorHandler('GetVesselSTructureLabelsOnEdges: ' // &
                'incompatible dimensions of input coordinates')
        end if 

        ! Edge labels
        if (allocated(labels)) then 
            if (size(labels) /= nq) then 
                deallocate(labels)
                allocate(labels(nq)) 
            end if
        else
            allocate(labels(nq))
        end if 
        labels = 0

        ! Vertex labels
        call vessel%GetVesselStructureLabelsOnPoints(xq1, yq1, &
            labelsv1)
        call vessel%GetVesselStructureLabelsOnPoints(xq2, yq2, &
            labelsv2)

        ! Compute edge labels
        !====================
        do i = 1, nq
            ! Check labels
            ! Check label occurrence
            if (labelsv1(i, 2) == 0) then 
                ! Normally, label should simply be labelsv1(i, 1), but 
                ! we do some sanity checks
                if (labelsv2(i, 2) == 0) then 
                    if (labelsv1(i, 1) /= labelsv2(i, 1)) then 
                        ! Different labels, unexpected - throw warning
                        print *, 'GetVesselStructureLabelsOnEdges: each edge vertex ' // & 
                            'belongs to a single yet different structure, ' // & 
                            'unexpected. Taking first structure value...'
                    end if
                else
                    if (.not. any(labelsv1(i, 1) == labelsv2(i, 1:2))) then 
                        ! Different labels, unexpected - throw warning
                        print *, 'GetVesselStructureLabelsOnEdges: each edge vertex ' // & 
                            'belongs to a single yet different structure, ' // & 
                            'unexpected. Taking first structure value...'
                    end if 
                end if 

                ! Add
                labels(i) = labelsv1(i, 1)
            else 
                ! Need to check which one to take 
                if (labelsv2(i, 2) == 0) then 
                    ! Normally, label should simply be labelsv2(i, 1), but 
                    ! we do some sanity checks
                    if (.not. any(labelsv2(i, 1) == labelsv1(i, 1:2))) then 
                        ! Different labels, unexpected - throw warning
                        print *, 'GetVesselStructureLabelsOnEdges: each edge vertex ' // & 
                            'belongs to a single yet different structure, ' // & 
                            'unexpected. Taking first structure value...'
                    end if 

                    ! Add
                    labels(i) = labelsv2(i, 1)
                else
                    ! Both are non-negative: need to check all combinations
                    if (any(labelsv1(i, 1) == labelsv2(i, 1:2)) .and. &
                        any(labelsv1(i, 2) == labelsv2(i, 1:2))) then 
                        ! Cannot distinguish between both possibilities, unexpected
                        print *, 'GetVesselStructureLabelsOnEdges: each edge vertex ' // & 
                            'belongs to two structures, cannot discriminate ' // & 
                            'which one to take - unexpected. Taking first structure value...'   
                        labels(i) = labelsv1(i, 1)
                    elseif (any(labelsv1(i, 1) == labelsv2(i, 1:2)) .and. &
                        .not. any(labelsv1(i, 2) == labelsv2(i, 1:2))) then 
                        labels(i) = labelsv1(i, 1)
                    elseif (any(labelsv1(i, 2) == labelsv2(i, 1:2)) .and. &
                        .not. any(labelsv1(i, 1) == labelsv2(i, 1:2))) then 
                        labels(i) = labelsv1(i, 2)
                    else
                        ! No agreement found - unexpected    
                        print *, 'GetVesselStructureLabelsOnEdges: each edge vertex '  // & 
                            'belongs to multiple yet different structures, ' // & 
                            'unexpected. Taking first structure value...'
                        labels(i) = labelsv1(i, 1)
                    end if 
                end if 
            end if   
        end do 

    end subroutine

    ! Updating (assumed already initialized)
    subroutine UpdateVesselCoordinates(vessel, xv, yv)

        ! Description
        !============
        ! This routine updates the vessel coordinates and all other 
        ! dependent structures according to the new vessel coordinates
        ! given by the xv, yv pairs. It is assumed that xv, yv has an 
        ! equal number of entries as there are vertices in the vessel
        ! polygon structure. 

        ! Notes
        !======
        ! Note 1: the targets are NOT updated yet!
        ! Note 2: we must construct the polygon sets using the polygon
        ! vertices, not simply passing the new set of coordinates! This
        ! will lead to an open polygon by default...

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)                :: vessel 
        real(R8), intent(in)            :: xv(:), yv(:)

        ! Auxiliary
        integer(I8)                     :: flag, npv, npvtot
        real(R8), allocatable           :: xvp(:), yvp(:)
        character(:), allocatable       :: vesselpath 

        ! Loop 
        integer(I8)                     :: i, k 

        ! Check
        !======
        npvtot = 0
        do i = 1, vessel%polygonset%np 
            npvtot = npvtot + vessel%polygonset%polygons(i)%nv 
        end do
        if ( (npvtot /= size(xv, 1)) .or. (size(xv, 1) /= size(yv, 1)) ) then 
            ! Incompatible dimensions
            call gdErrorHandler('UpdateVesselCoordinates: incompatible ' // &
                'dimensions of new coordinates and original vessel polygon')
        end if  

        ! Adjust coordinates
        !===================
        ! Vessel polygon vertices
        k = 0
        do i = 1, vessel%polygonset%np
            ! Get number of vertices of this polygon
            npv = vessel%polygonset%polygons(i)%nv 

            ! Get the coordinates
            xvp = xv(k+1:k+npv)
            yvp = yv(k+1:k+npv)

            ! Assign - don't reconstruct, may alter vertex oder!
            vessel%polygonset%polygons(i)%x = xvp 
            vessel%polygonset%polygons(i)%y = yvp

            ! H
            !call vessel%polygonset%polygons(i)%Construct(&
            !   xv(vessel%polygonset%polygons(i)%vert), &
            !    yv(vessel%polygonset%polygons(i)%vert), &
            !    vessel%polygonset%polygons(i)%labels(vessel%polygonset%polygons(i)%vert,:))
            
            ! Update counter
            k = k + npv 
        end do

        ! Reconstruct polygonset (just in case, shouldn't be necessary)
        call vessel%polygonset%Construct(vessel%polygonset%polygons)

        ! Test orientation
        call vessel%polygonset%OrientNestedClosedPolygons(flag)

        ! Check if vessel polygon is still a closed and simple polygon
        if (flag .ne. 0) then  
            ! Throw error, but do not exit program - may be dealt with upstream
            print *, 'flag: ', flag
            call gdErrorHandler('UpdateVesselCoordinates: could not orient ' // &
                'polygons, OrientNestedClosedPolygons exited with flag above.' // &
                'Not updating polygon levelset function any further', severityin=0)
            return 
        end if 

        ! Write data
        vesselpath = 'vesselpolygon'
        call vessel%polygonset%WriteData(vesselpath)

        ! Adjust vessel description
        !==========================
        call vessel%plfvessel%Initialize(vessel%polygonset)

    end subroutine

    !------------------------------------------------------------------!
    !                         DivGeo coupling                          !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeDGData(dgdata)

        ! Description
        !============
        ! Initialize and allocate all data to default values

        ! Declare variables
        !==================
        ! Arguments
        class(DivGeoDataUDT)        :: dgdata 

        ! Initialize
        !===========
        ! Logicals
        dgdata%hasEl            = .false. 
        dgdata%hasElFcLbl       = .false.
        dgdata%hasElVessel      = .false. 
        dgdata%hasElVoid        = .false. 

        ! Elements & points
        dgdata%nel      = 0
        dgdata%nv       = 0 

        ! Allocate
        !=========
        call dgdata%Allocate()
        
        
    end subroutine

    ! Allocation
    subroutine AllocateDGData(dgdata)

        ! Description
        !============
        ! Allocate the divgeo data object, assuming that all relevant 
        ! sizes have been initialized.

        ! Declare variables
        !==================
        ! Arguments
        class(DivGeoDataUDT)        :: dgdata 

        ! First deallocate
        !=================
        ! Soft deallocation
        call dgdata%Deallocate()
        
        ! Then allocate
        !==============
        ! Element & points
        allocate(dgdata%elID(dgdata%nel), dgdata%elfcLbl(dgdata%nel), &
            dgdata%elvx(dgdata%nv), dgdata%elvy(dgdata%nv), &
            dgdata%elv1(dgdata%nel), dgdata%elv2(dgdata%nel))
        
        ! Vessel elements
        allocate(dgdata%elvessel(dgdata%nel)) 

    end subroutine

    ! Deallocation
    subroutine DeallocateDGData(dgdata)

        ! Description
        !============
        ! Deallocate the divgeo data object (soft)
        ! Declare variables
        !==================
        ! Arguments
        class(DivGeoDataUDT)        :: dgdata 

        ! Soft deallocation
        !==================
        ! Element & points
        if (allocated(dgdata%elID)) deallocate(dgdata%elID)
        if (allocated(dgdata%elfcLbl)) deallocate(dgdata%elfcLbl)
        if (allocated(dgdata%elv1)) deallocate(dgdata%elv1)
        if (allocated(dgdata%elv2)) deallocate(dgdata%elv2)
        if (allocated(dgdata%elvx)) deallocate(dgdata%elvx)
        if (allocated(dgdata%elvy)) deallocate(dgdata%elvy)

        ! Vessel elements
        if (allocated(dgdata%elvessel)) deallocate(dgdata%elvessel)

    end subroutine

    ! DivGeo .dgo reader
    subroutine ReadDGFile(dgdata, filepath)

        ! Description
        !============
        ! This routine allocates and reads in data from a .dgo file. The
        ! following fields are read in:
        ! - p1: list of first vertices of elements (x, y, z coordinates)
        ! - p2: list of second vertices of elements (x, y, z coordinates)
        ! - fcLbl: list of element labels
        ! - vess_elm: list of vessel element indices

        ! From this data, different other base quantities are derived, 
        ! such as the element vertices (as a ne-by-2 array) and the 
        ! vertices themselves. If data is not present, then the 
        ! corresponding logicals are set to false. 

        ! Declare variables
        !==================
        ! Arguments
        class(DivGeoDataUDT)            :: dgdata
        character(*), intent(in)        :: filepath 

        ! Auxiliary
        integer(I8)                             :: np1, np2, nfcLbl
        integer(I8), allocatable, dimension(:)  :: tempi, fcLbl, v1, v2, &
            worki 
        logical                                 :: reachedeof, hasp1, &
            hasp2, islegal 
        real(R8), allocatable, dimension(:)     :: tempr, x1, x2, &
            y1, y2, z1, z2, x, y
        character(:), allocatable       :: thisline

        ! Loop

        ! Data
        integer         :: fid
        data fid /70/
        
        ! Initialize
        !===========
        ! Initialize structure to defaults
        call dgdata%Initialize()

        ! Open file
        !==========
        ! Print from where we're reading
        print *, 'reading divgeo output from file: ', filepath

        ! Open the file
        open(unit = fid, file = filepath)

        ! Read file
        !==========
        ! Read elements & points
        !-----------------------
        ! Initialize
        hasp1 = .true. 
        hasp2 = .true. 
        np1 = 0
        np2 = 0

        ! Read until we find 'p1'
        call ReadUntilMatchFound(fid, 'p1', ' ', reachedeof)
        if (reachedeof) then 
            ! Set to false
            hasp1 = .false. 
        else
            ! Initialize
            allocate(x1(0), y1(0), z1(0))

            ! Read in data
            do while (.true.)
                ! Read in the next triplet
                call ReadSingleLine(fid, thisline, reachedeof)

                ! Check for EOF
                if (reachedeof) then 
                    exit 
                end if 

                ! Check if we can extract coordinates
                call ExtractRealFromString1D(thisline, tempr, islegal)

                ! Check if the read was legal
                if (.not. islegal) then 
                    exit
                end if 
                if (size(tempr) /= 3) then 
                    exit
                end if 

                ! Add the coordinates
                x1 = [x1, tempr(1)]
                y1 = [y1, tempr(2)]
                z1 = [z1, tempr(3)]

            end do 

            ! Determine size
            np1 = size(x1)
        end if 

        ! Rewind the file
        rewind(fid)

        ! Read until we find 'p2'
        call ReadUntilMatchFound(fid, 'p2', ' ', reachedeof)
        if (reachedeof) then 
            ! Set to false
            hasp2 = .false. 
        else
            ! Initialize
            allocate(x2(0), y2(0), z2(0))

            ! Read in data
            do while (.true.)
                ! Read in the next triplet
                call ReadSingleLine(fid, thisline, reachedeof)

                ! Check for EOF
                if (reachedeof) then 
                    exit 
                end if 

                ! Check if we can extract coordinates
                call ExtractRealFromString1D(thisline, tempr, islegal)

                ! Check if the read was legal
                if (.not. islegal) then 
                    exit
                end if 
                if (size(tempr) /= 3) then 
                    exit
                end if 

                ! Add the coordinates
                x2 = [x2, tempr(1)]
                y2 = [y2, tempr(2)]
                z2 = [z2, tempr(3)]

            end do 

            ! Determine size
            np2 = size(x2)
        end if 

        ! Sanity checks
        if (.not. hasp1 .and. .not. hasp2) then 
            ! Print
            print *, 'ReadDGData: element data not available'

            ! Set logicals
            dgdata%hasEl = .false.

        elseif (hasp1 .and. hasp2) then 
            ! Check if dimensions are consistent
            if (np1 /= np2) then 
                call gdErrorHandler('ReadDGData: inconsistent number of ' // & 
                    'p1 and p2, check input file')
            else
                ! Set logicals
                dgdata%hasEl = .true. 

                ! Set number of elements
                dgdata%nel = np1 
                print *, 'ReadDGData: number of elements = ', dgdata%nel
            end if 
        else
            ! Inconsistent 
            call gdErrorHandler('ReadDGData: only one set of point ' // & 
                'coordinates found, not supported. Check input file')
        end if 

        ! Read fcLbl
        !-----------
        ! Only if elements are present
        dgdata%haselfcLbl = .false.
        if (dgdata%hasEl) then 
            ! Rewind the file
            rewind(fid)

            ! Read until we find 'fclbl'
            dgdata%haselfcLbl = .true.
            call ReadUntilMatchFound(fid, 'fclbl', ' ', reachedeof)
            if (reachedeof) then 
                ! Set to false
                dgdata%haselfcLbl = .false. 
            else
                ! Initialize
                allocate(fcLbl(0))

                ! Read in data
                do while (.true.)
                    ! Read in the next triplet
                    call ReadSingleLine(fid, thisline, reachedeof)

                    ! Check for EOF
                    if (reachedeof) then 
                        exit 
                    end if 

                    ! Check if we can extract 
                    call ExtractIntegerFromString1D(thisline, tempi, islegal)

                    ! Check if the read was legal
                    if (.not. islegal) then 
                        exit
                    end if 
                    if (size(tempi) /= 1) then ! Expected size 1
                        exit
                    end if 

                    ! Add the coordinates
                    fcLbl = [fcLbl, tempi]

                end do 

                ! Determine size
                nfcLbl = size(fcLbl)

                ! Sanity checks
                if (nfcLbl /= dgdata%nel) then 
                    call gdErrorHandler('ReadDGData: number of fcLbl does ' // & 
                        'not equal number of elements, unexpected')
                end if 

                ! Add
                dgdata%elfcLbl = fcLbl

            end if 
        end if 

        ! Print 
        if (dgdata%haselfcLbl) then 
            print *, 'ReadDGData: fcLbl read in'
        else
            print *, 'ReadDGData: could not read fcLbl'
        end if 

        ! Read vessel elements
        !---------------------
        ! Only if elements are present
        if (dgdata%hasEl) then 
            ! Rewind the file
            rewind(fid)

            ! Read until we find 'vess_elm'
            dgdata%haselvessel = .true.
            call ReadUntilMatchFound(fid, 'vess_elm', ' ', reachedeof)
            if (reachedeof) then 
                ! Set to false
                dgdata%haselvessel = .false. 
            else
                ! Initialize
                allocate(worki(0))

                ! Read in data
                do while (.true.)
                    ! Read in the next triplet
                    call ReadSingleLine(fid, thisline, reachedeof)

                    ! Check for EOF
                    if (reachedeof) then 
                        exit 
                    end if 

                    ! Check if we can extract 
                    call ExtractIntegerFromString1D(thisline, tempi, islegal)

                    ! Check if the read was legal
                    if (.not. islegal) then 
                        exit
                    end if 
                    if (size(tempi) /= 1) then ! Expected size 1
                        exit
                    end if 

                    ! Add the coordinates
                    worki = [worki, tempi]

                end do 

                ! Add
                dgdata%elvessel = worki

                ! Housekeeping
                deallocate(worki)

            end if 
        end if 

        ! Print 
        if (dgdata%haselvessel) then 
            print *, 'ReadDGData: vessel elements read in'
        else
            print *, 'ReadDGData: could not read vessel elements'
        end if 

        ! Read void elements
        !-------------------
        ! Only if elements are present
        if (dgdata%hasEl) then 
            ! Rewind the file
            rewind(fid)

            ! Read until we find 'trimark'
            dgdata%haselvoid = .true.
            call ReadUntilMatchFound(fid, 'trimark', ' ', reachedeof)
            if (reachedeof) then 
                ! Set to false
                dgdata%haselvoid = .false. 
            else
                ! Initialize
                allocate(worki(0))

                ! Read in data
                do while (.true.)
                    ! Read in the next triplet
                    call ReadSingleLine(fid, thisline, reachedeof)

                    ! Check for EOF
                    if (reachedeof) then 
                        exit 
                    end if 

                    ! Check if we can extract 
                    call ExtractIntegerFromString1D(thisline, tempi, islegal)

                    ! Check if the read was legal
                    if (.not. islegal) then 
                        exit
                    end if 
                    if (size(tempi) /= 1) then ! Expected size 1
                        exit
                    end if 

                    ! Add the coordinates
                    worki = [worki, tempi]

                end do 

                ! Add
                dgdata%trimark = worki

                ! Housekeeping
                deallocate(worki)

            end if 
        end if 

        ! Print 
        if (dgdata%haselvessel) then 
            print *, 'ReadDGData: void elements read in'
        else
            print *, 'ReadDGData: could not read void elements'
        end if 


        ! Housekeeping
        !=============
        ! Close the file
        close(fid)

        ! Add additional data
        !====================
        ! If elements are present, we extract the element vertices
        ! and the local vertex numbering. These are stored in x, y and 
        ! v1, v2 arrays
        if (dgdata%hasEl) then 
            ! Extract
            call ExtractEdgesFromCoordinates(x1, y1, x2, y2, &
                v1, v2, x, y)

            ! Add
            dgdata%elvx = x 
            dgdata%elvy = y 
            dgdata%elv1 = v1 
            dgdata%elv2 = v2

            ! rescale x, y (originally in mm)
            dgdata%elvx = dgdata%elvx/1e3_R8
            dgdata%elvy = dgdata%elvy/1e3_R8
        end if 



    end subroutine

    ! Structure extraction
    subroutine ExtractDGStructures(dgdata, structures, elementIDs, flag)

        ! Description
        !============
        ! This routine extracts structure data from the dgdata
        ! structure, assuming that all required fields have been 
        ! properly initialized. This is only possible if elements are
        ! present. If in addition labels are present, then the
        ! structures are further subdivided based on these labels. The 
        ! flag is zero if the structures were succesfully extracted, 
        ! otherwise it has a positive value. Structure IDs are set 
        ! equal to the face label found in DivGeo if present, otherwise
        ! it simply goes from 1 to the number of structures.

        ! Note: since in general an arbitrary number of elements may 
        ! coincide, we leave it up to the user to specify which subset 
        ! of elements (defined as the array 'elementIDs') that should 
        ! be taken to form structures. It is then assumed that all these
        ! elements form simple (open or closed) non-branching polygons. 

        ! Note: under the hood this routine uses the polygon edge 
        ! sorter to determine all structures etc. 

        ! Declare variables
        !==================
        ! Arguments
        class(DivGeoDataUDT)            :: dgdata 
        type(VesselStructureUDT), allocatable, intent(out)  :: &
            structures(:)
        integer(I8), dimension(:), intent(in)               :: elementIDs
        integer(I8), intent(out)        :: flag 

        ! Auxiliary
        integer(I8)                             :: si, ei, nel
        integer(I8), allocatable, dimension(:)  :: fcLbl, &
            sortindex, pv 
        integer(I8), allocatable                :: pein(:, :), &
            sortededges(:, :)
        logical, allocatable, dimension(:)      :: ispolygonstart, &
            isbranchingpolygon
        
        ! Loop
        integer(I8)                             :: k 

        ! Initialize
        !===========
        ! Set flag to success
        flag = 0

        ! Unpack
        nel = size(elementIDs)

        ! Hedge for trivial case
        if (size(elementIDs) == 0) then 
            allocate(structures(0))
            return 
        end if

        ! Check if elements are present
        if (.not. dgdata%hasEl) then 
            print *, 'ExtractDGStructures: no elements present, ' // & 
                'cannot extract structures. Returning...'
            flag = 1
            return 
        end if

        ! Check if labels are present
        if (.not. dgdata%hasElFcLbl) then 
            print *, 'ExtractDGStructures: no element fcLbl present, ' // & 
                'extracting based on element data only and setting label ' // & 
                'equal to structure ID' 
            allocate(fcLbl(dgdata%nel))
            fcLbl = 0
        else
            ! Get labels
            fcLbl = dgdata%elfcLbl
        end if
        
        ! Extract structures
        !===================
        ! Sort polygon edges
        allocate(pein(nel, 2))
        pein(:, 1) = dgdata%elv1(elementIDs)
        pein(:, 2) = dgdata%elv2(elementIDs)
        allocate(sortindex(nel), ispolygonstart(nel), &
            isbranchingpolygon(nel))
        call SortPolygonEdges(pein, nel, sortindex, &
            ispolygonstart, isbranchingpolygon)
        sortededges = pein(sortindex, :)

        ! Sanity checks
        if (any(isbranchingpolygon)) then 
            ! This is not supported
            flag = 2
            print *, 'ExtractDGStructures: branching polygons present, ' // & 
                'cannot proceed. Returning...'
            return 
        end if 

        ! Extract polygons
        allocate(structures(count(ispolygonstart)))
        k = 0
        si = 0
        ei = 0
        do while (k < count(ispolygonstart)) 
            ! Update counter
            k = k + 1

            ! Find the next polygon indices
            si = findloc(ispolygonstart(ei+1:), .true., 1, back=.false.) + ei 
            ei = findloc(ispolygonstart(si+1:), .true., 1, back=.false.) + si
            if (ei == si) then ! no further start found, so we reached the end
                ei = size(ispolygonstart)
            else
                ei = ei - 1
            end if   

            ! Get the vertices
            call ExtractPolygonVertices(sortededges(si:ei, :), ei-si+1, pv)

            ! Add to the structure
            structures(k)%x = dgdata%elvx(pv)
            structures(k)%y = dgdata%elvy(pv)
            structures(k)%np = size(structures(k)%x)
            if (pv(1) == pv(size(pv))) then 
                structures(k)%isclosed = .true.
            else
                structures(k)%isclosed = .false. 
            end if 
            structures(k)%ID = int(k, kind=I4) 
            
            ! Check label
            if (dgdata%hasElFcLbl) then 
                structures(k)%label = int(fcLbl(elementIDs(sortindex(si))), kind=I4)
                if (.not. all(fcLbl(elementIDs(sortindex(si:ei))) == &
                    fcLbl(elementIDs(sortindex(si))))) then 
                    print *, 'ExtractDGStructures: found multiple labels for ' // & 
                        'structure ', k, ' taking first label: ', structures(k)%label 
                end if 
            else
                structures(k)%label = structures(k)%ID 
            end if 
        end do 

    end subroutine

    ! Vessel structure extraction
    subroutine ExtractDGVesselStructures(dgdata, structures, flag)

        ! Description
        !============
        ! This routine extracts structures that are part of the vessel
        ! structure. These structures are defined as structures with a 
        ! non-zero label. The routine extracts the vessel structures by
        ! checking for each unique vessel label (propagated to 
        ! structure.label) which elements belong to it, and how many
        ! polygons (vessel structures) that can be constructed from it.
        ! Normally, only one, non-brancching, non-selfintersecting 
        ! polygon should emerge, otherwise we throw an error.
         
        ! Note: vessel elements are assumed to hold a non-zero label
        ! (it can be negative). Furthermore, vessel polygons are 
        ! assumed to be non-branching, simple polygons, and should have
        ! the same label on all elements

        ! Declare variables
        !==================
        ! Arguments
        class(DivGeoDataUDT)            :: dgdata
        type(VesselStructureUDT), allocatable, intent(out)  :: structures(:) 
        integer(I8), intent(out)        :: flag

        ! Auxiliary
        integer(I8)                 :: nvs, tfcLbl 
        integer(I8), allocatable, dimension(:)  :: fcLblu, vessfcLbl, &
            tel
        type(VesselStructureUDT), allocatable, dimension(:)     :: &
            tempstructures

        ! Loop
        integer(I8)                 :: i

        ! Extract structures
        !===================
        ! Set to success
        flag = 0

        ! Check if vessel elements are present
        if (.not. dgdata%hasElVessel) then 
            print *, 'ExtractDGVesselStructures: no vessel element data ' // &
                'present, cannot proceed. Returning...'
            flag = 1
            return 
        end if 
        if (.not. dgdata%hasElFcLbl) then 
            print *, 'ExtractDGVesselStructures: fcLbl not present, ' // &
                'which is required for vessel data. Returning...'
            flag  = 2
            return 
        end if 
        
        ! Check if all vessel labels are non-zero
        if (any(dgdata%elfcLbl(dgdata%elvessel) == 0)) then 
            print *, 'ExtractDGVesselStructures: vessel segments with ' // &
                'zero label detected, not supported. Check dg setup. ' // & 
                'Returning...'
            flag = 3
            return 
        end if

        ! Compute unique number of face labels
        vessfcLbl = dgdata%elfcLbl(dgdata%elvessel)
        call Unique(vessfcLbl, fcLblu)
        nvs = size(fcLblu)

        ! Allocate
        allocate(structures(nvs))

        ! Determine structures
        do i = 1, nvs 
            ! Get labels for this unique label
            tfcLbl = fcLblu(i)

            ! Check which elements have this label
            allocate(tel(count(vessfcLbl == tfcLbl)))
            tel = pack(dgdata%elvessel, vessfcLbl == tfcLbl)

            ! Extract the structures
            call dgdata%ExtractStructures(tempstructures, tel, flag)

            ! Checks
            if (flag /= 0) then 
                print *, 'ExtractDGVesselStructures: could not extract ' // & 
                    'structures for label', fcLblu(i), ', returning...'
                return
            end if 
            if (size(tempstructures) /= 1) then 
                print *, 'ExtractDGVesselStructures: did not find exactly ' // &
                    'one structure for label', fcLblu(i), ', returning...'
                flag = 4
                return 
            end if
            
            ! Add
            structures(i) = tempstructures(1)
            structures(i)%ID = int(i, kind=I4) 
            structures(i)%label = int(fcLblu(i), kind=I4)

            ! Housekeeping
            deallocate(tel)
        end do 

        ! Write
        call WriteStructureFile('vessel_dgo', structures)

    end subroutine
 
    ! Void structure extraction
    subroutine ExtractDGTriangulationStructures(dgdata, structures, flag)

        ! Description
        !============
        ! This routine extracts structures that are part of the structure
        ! that should be considered for triangulation (and EIRENE) by
        ! an external triangulation routine. Additional elements may
        ! be introduced (these elements may overlap with vessel elements 
        ! as long as they form a set of closed non-intersecting
        ! polygons). Additionally, elements with negative fcLbl are 
        ! ignored as they are assumed to be 'fake' vessel boundaries 
        ! used to limit the plasma. 

        ! Declare variables
        !==================
        ! Arguments
        class(DivGeoDataUDT)            :: dgdata
        type(VesselStructureUDT), allocatable, intent(out)  :: structures(:) 
        integer(I8), intent(out)        :: flag

        ! Auxiliary
        integer(I8)                 :: nvs, tfcLbl, flagv 
        integer(I8), allocatable, dimension(:)  :: fcLblu, vessfcLbl, &
            tel, temptrimark 
        type(VesselStructureUDT), allocatable, dimension(:)     :: &
            tempstructures, vesselstructures

        ! Loop
        integer(I8)                 :: i, j, k, sID 

        ! Initialize
        !===========
        ! Set to success
        flag = 0

        ! Check if vessel elements are present
        if (.not. dgdata%hasElVessel) then 
            print *, 'ExtractDGTriangulationStructures: no vessel element data ' // &
                'present, cannot proceed. Returning...'
            flag = 1
            return 
        end if 
        if (.not. dgdata%hasElFcLbl) then 
            print *, 'ExtractDGTriangulationStructures: fcLbl not present, ' // &
                'which is required for vessel data. Returning...'
            flag  = 2
            return 
        end if 
        if (.not. dgdata%hasElVoid) then 
            print *, 'ExtractDGTriangulationStructures: void element data not present, ' // &
                'which is required for triangulation data. Returning...'
            flag = 3
        end if 
        
        ! Check if all vessel labels are non-zero
        if (any(dgdata%elfcLbl(dgdata%elvessel) == 0)) then 
            print *, 'ExtractDGTriangulationStructures: vessel segments with ' // &
                'zero label detected, not supported. Check dg setup. ' // & 
                'Returning...'
            flag = 3
            return 
        end if

        ! Extract standard vessel structures
        !===================================
        ! Extract vessel structures using dedicated routine
        call ExtractDGVesselStructures(dgdata, vesselstructures, flagv)

        ! Sanity check
        if (flagv /= 0) then 
            print *, 'flag of vessel structure extraction = ', flagv 
            print *, 'ExtractDGTriangulationStructures: could not extract ' // & 
                'vessel structures. Returning...'
            flag = 5
            return 
        end if 

        ! Extract additional void structures
        !===================================
        ! Set any trimark of vessel edges to zero - should not be added
        ! twice 
        temptrimark = dgdata%trimark 
        temptrimark(dgdata%elvessel) = 0

        ! Compute unique number of polygon labels
        allocate(vessfcLbl(count(temptrimark /= 0)))
        vessfcLbl = pack(temptrimark, temptrimark /= 0)
        call Unique(vessfcLbl, fcLblu)
        nvs = size(fcLblu)

        ! Allocate
        allocate(structures(0)) ! not a priori known how many polygons we'll have

        ! Determine structures
        sID = 0
        do i = 1, nvs 
            ! Get labels for this unique label
            tfcLbl = fcLblu(i)

            ! Check which elements have this label
            allocate(tel(count(temptrimark == tfcLbl)))
            tel = pack([(k, k = 1, dgdata%nel)], temptrimark == tfcLbl)

            ! Extract the structures
            call dgdata%ExtractStructures(tempstructures, tel, flag)

            ! Checks
            if (flag /= 0) then 
                print *, 'ExtractDGVesselStructures: could not extract ' // & 
                    'structures for label', fcLblu(i), ', returning...'
                return
            end if 

            ! Adjust ID
            do j = 1, size(tempstructures)
                sID = sID + 1
                tempstructures(j)%ID = int(sID, kind=I4)
                tempstructures(j)%label = int(fcLblu(i), kind=I4)
            end do 
            
            ! Add
            structures = [structures, tempstructures]

            ! Housekeeping
            deallocate(tel)
        end do 

        ! Append vessel structures with non-negative label (and overwrite
        ! their ID)
        do i  = 1, size(vesselstructures)
            if (vesselstructures(i)%label > 0) then 
                ! Update counter
                sID = sID + 1 

                ! Append structure
                structures = [structures, vesselstructures(i)]

                ! Update ID
                structures(size(structures))%ID = int(sID, kind=I4)
            end if 
        end do

        ! Write
        call WriteStructureFile('vessel_tria_dgo', structures)

    end subroutine

    ! Structure file writing
    subroutine WriteStructureFile(filename, structures)

        ! Description
        !============
        ! This routine writes out a structure file (filename specified
        ! through filename variable) containing the structures given 
        ! in 'structures' in a structure.dat file format.

        ! Declare modules
        !================
        use mod_std_formatspecs

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)        :: filename 
        type(VesselStructureUDT), intent(in)    :: structures(:)
        
        ! Auxiliary
        integer(I8)                     :: ns, fu 
        character(:), allocatable       :: fmt 

        ! Loop
        integer(I8)                     :: i, j

        ! Initialize
        !===========
        ! Open file 
        open (action='write', file=trim(filename // '.dat'), newunit=fu, &
            status='unknown')

        ! Compute size
        ns = size(structures)

        ! Write number of structures
        write(fu, *) ns 
        write(fu, *) '$structures'

        ! Set writing format
        fmt = '(' // Rfm // ',' // spacefm // ','// Rfm // ')'

        ! Write structure data
        !=====================
        do i = 1, ns
            ! Write structure header
            write (fu, *) 'Structure ', structures(i)%label

            ! Write structure polygon
            write (fu, *) structures(i)%np 
            do j = 1, structures(i)%np 
                write (fu, fmt) structures(i)%x(j), structures(i)%y(j)
            end do  
        end do

        ! Housekeeping
        close(fu)

    end subroutine

    !------------------------------------------------------------------!
    !                    Auxiliary data extraction                     !
    !------------------------------------------------------------------!
    ! Data extraction for grid generation
    subroutine ExtractGGData(magneticField, environment, options)

        ! Description
        !============
        ! Extract all necessary data to run the grid generator driver. This includes
        ! the magnetic field and environment. The options structure
        ! should contain the necessary paths to read in these quantities. 

        ! Declare variables
        !==================
        ! Arguments
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(GoatoptionsUDT)                :: options 

        ! Auxiliary 
        type(MagneticFieldOptionsUDT)       :: mfoptions 
        type(EnvironmentOptionsUDT)         :: environmentoptions 

        integer                 :: filespecifier(0:2)
        
        ! Declare data
        data filespecifier /50, 51, 52/ ! numbers for file specification

        ! Initialize
        !===========
        ! Set the path - the same as for general goat options 
        environmentoptions%inputfilepath    = options%inputfilepath
        mfoptions%inputfilepath             = options%inputfilepath 
        
        ! Set the options
        call environmentoptions%Set()
        call mfoptions%Set()

        ! Reset vessel reading 
        environmentoptions%vesselfilepath = options%structurefilepath

        ! Read data
        !==========
        ! Read magnetic field
        call ReadMagneticField(magneticField, mfoptions, options%magneticfieldfilepath)

        ! Read additional data
        !=====================
        call ConstructMagneticField(magneticField, mfoptions) 
        call ConstructEnvironment(environment, environmentoptions) 

        ! Housekeeping
        !=============


    end subroutine

    ! Data extraction for other goat drivers (GD, GA, ...)
    subroutine ExtractGoatData(grid, magneticField, environment, options)

        ! Description
        !============
        ! Extract all necessary data to run the goat drivers. This includes
        ! the grid, magnetic field, and environment. The options structure
        ! should contain the necessary paths to read in these quantities. 
    
        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT)                       :: grid 
        type(MagneticFieldUDT)              :: magneticField 
        type(EnvironmentUDT)                :: environment 
        type(GoatoptionsUDT)                :: options 
    
        ! Auxiliary 
        type(GridOptionsUDT)                :: gridoptions
        type(MagneticFieldOptionsUDT)       :: mfoptions 
        type(EnvironmentOptionsUDT)         :: environmentoptions 
    
        integer                 :: filespecifier(0:2)
        
        ! Declare data
        data filespecifier /50, 51, 52/ ! numbers for file specification
    
        ! Initialize
        !===========
        ! Set the path - the same as for general goat options 
        gridoptions%inputfilepath           = options%inputfilepath 
        environmentoptions%inputfilepath    = options%inputfilepath
        mfoptions%inputfilepath             = options%inputfilepath 
        
        ! Set the options
        call gridoptions%Set()
        call environmentoptions%Set()
        call mfoptions%Set()
    
        ! Reset label mappings for grid options
        gridoptions%facelabelmappingGD      = options%facelabelmappingGD
        gridoptions%facelabelmappingGG      = options%facelabelmappingGG
        gridoptions%facelabelsubfrom        = options%facelabelsubfrom
        gridoptions%facelabelsubto          = options%facelabelsubto
    
        ! Reset vessel reading 
        environmentoptions%vesselfilepath   = options%structurefilepath

        ! Carry over some environment options
        options%readstate                   = environmentoptions%readstate
        options%readstatemeth               = environmentoptions%readstatemeth
    
        ! Read data
        !==========
        ! Read grid
        select case (gridoptions%readmeth)
    
        case ('traduitb2us')
    
            ! Unstructured traduit file reading
            call ReadTraduitUS(grid, options%gridfilepath)
    
        case ('b2fgmtry_us')
    
            ! Structured traduit file reading
            call ReadB2fgmtryUS(grid, options%gridfilepath)
    
        case default
    
            call gdErrorHandler('ExtractGoatData: unknown readtype option')
    
        end select 
    
        ! Read magnetic field
        call ReadMagneticField(magneticField, mfoptions, options%magneticfieldfilepath)
    
        ! Read additional data
        !=====================
        call ReadGoatGGData(grid, 'goatggdata.dat')
        call ConstructGrid(grid, gridoptions)
        call ConstructMagneticField(magneticField, mfoptions) 
        call ConstructEnvironment(environment, environmentoptions) 
    
        ! Housekeeping
        !=============
    
    
    end subroutine
    
    ! Vessel vertex pairs
    subroutine GetVesselVertexPairs(vessel, vpairs, structureIDs, vertIDs)

        ! Description
        !============
        ! This routine returns the vessel vertex pairs (i.e. the 
        ! subsequent edge pairs) in a nvpairs-by-3 structure. The 
        ! edges are formed by v1v2 and v2v3 (in that specific order). 
        ! This routine is useful for many shape optimization related 
        ! cost functions and constraints. Only vertex pairs of the 
        ! specified structure IDs or vertices are taken (for the latter,
        ! the vertex pair of v2, i.e. v1v2 v2v3 is stored). Note that 
        ! the first and last pair may not be part of the specified 
        ! structure. 

        ! Declare variables
        !==================
        ! Arguments
        class(VesselUDT)                        :: vessel 
        integer(I8), allocatable, intent(out)   :: vpairs(:, :)
        integer(I8), intent(in)                 :: structureIDs(:), &
            vertIDs(:)

        ! Auxiliary
        logical, allocatable, dimension(:)      :: includevert, &
            ispolygonstart, ispolygonend
        integer(I8)                             :: nv, nvpairs, psind, &
            prevv, ind
        integer(I8), allocatable                :: labels(:, :)
        integer(I8), allocatable, dimension(:)  :: pID, vID
        real(R8), allocatable, dimension(:)     :: xv, yv

        ! Loop
        integer(I8)                             :: i, cc
 
        ! Initialize
        !===========
        ! Associate
        associate(ps        => vessel%polygonset)

        ! Determine vertex pairs
        !=======================
        ! Get all vertices and vertex labels
        call ps%GetLabels(labels)
        call ps%GetVertices(xv, yv, pID)
        call ps%GetVertices(vID)
        nv = size(labels, 1)

        ! Check if contributions should be included
        allocate(includevert(nv))
        includevert = .false. 

        ! Constrain per vessel structure (label 1 and 2)
        do i = 1, size(structureIDs)
            ! Unpack ID
            associate(tID       => structureIDs(i))

            ! Check vertices
            where ( (labels(:, 1) == tID) .or. (labels(:, 2) == tID) ) &
                includevert = .true. 

            ! Housekeeping
            end associate
        end do

        ! Constrain per vertex ID
        do i = 1, size(vertIDs)
            ! Unpack ID
            associate(tID       => vertIDs(i))

            ! Check vertices
            where( (labels(:, 3) == tID)) includevert = .true. 

            ! Housekeeping
            end associate
        end do

        ! Check starting points of polygon
        ispolygonstart = [.true., pID(2:) - pID(:nv-1) /= 0]
        ispolygonend = [ispolygonstart(2:), .true.]

        ! Initialize vertex pairs 
        nvpairs = count(includevert .and. (.not. ispolygonend))
        allocate(vpairs(nvpairs, 3))

        ! Loop
        cc = 0
        psind = 1
        do i = 1, size(vID)
            if (includevert(vID(i))) then 

                ! Check if we should include the pair (need to account
                ! for starting vertex appearing twice in vID)
                if (ispolygonend(i)) then 
                    ! Update psind
                    psind = psind + ps%polygons(pID(i))%ne + 1
                    cycle 
                end if 

                ! Update counter
                cc = cc + 1

                ! Set current vertex
                vpairs(cc, 2) = vID(i) 

                ! Check
                if (ispolygonstart(i)) then
                    ! Previous vertex ID should be current polygon start 
                    ! plus number of edges minus 1
                    ind = psind + ps%polygons(pID(i))%ne - 1
                    prevv = vID(ind) 
                else 
                    prevv = vID(i-1)
                end if 

                ! Add
                vpairs(cc, 1) = prevv 
                vpairs(cc, 3) = vID(i+1)
            end if 
        end do

        ! Sanity check
        if (cc /= nvpairs) then 
            ! Should be a bug
            call gdErrorHandler('GetVesselVertexPairs: this is a bug')
        end if 

        ! Housekeeping
        end associate

    end subroutine

    ! Structures from polygonset
    subroutine GetStructuresFromPolygonset(polygonset, structures)

        ! Description
        !============
        ! Simple routine that converts each polygon of a polygon set
        ! to a structure datatype. Useful to write out polygonset data
        ! in structure.dat format for plotting. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT), intent(in)        :: polygonset 
        type(VesselStructureUDT), allocatable, intent(out)  :: structures(:)

        ! Auxiliary

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Allocate
        allocate(structures(polygonset%np))

        ! Construct structures
        !=====================
        do i = 1, polygonset%np
            ! Associate for ease
            associate(pol => polygonset%polygons(i))

            ! Construct structure
            structures(i)%np = pol%ne+1
            structures(i)%x = pol%x(pol%vert)
            structures(i)%y = pol%y(pol%vert)
            structures(i)%isclosed = pol%isclosed
            structures(i)%ID = int(i, kind=I4)
            structures(i)%label = 0
            
            ! Housekeeping
            end associate
        end do 
    end subroutine

end module
